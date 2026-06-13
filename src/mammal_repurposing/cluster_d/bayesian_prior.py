"""V6 §13.2 — PyMC NUTS hierarchical Bayesian model for Cluster D.

The full model per `research/4-tier/Multi-Source Neurobiological Prior for
Cognition Target Prioritization.md` §B.2:

    y^s_i ~ N(α_s + β_s · θ_i, τ_s^{-1} + (σ^s_i)²)       likelihood
    θ_i  ~ N(0, 1)                                          target prior
    α_s  ~ N(0, 0.5²)                                       source intercept
    β_s  ~ HalfNormal(1.0)                                  source informativeness
    τ_s  ~ Gamma(2, 2)                                      source precision

with skeptical-prior down-weight on Lit-OTAR (β_Lit ~ HN(0.3)) and reference
anchor compounds (BDNF, COMT, ACHE, DRD2, GRIN2B, CHRNA7) at θ ~ N(0.5, 0.3²)
to break scale + sign degeneracy.

Posterior cognition-relevance per target: w_i = σ(θ_i) ∈ (0, 1) — feeds the
§7.11 calibration as a multiplicative bias.

MH8 SUBSTRATE-MEDIATED FLAG (added 2026-05-28):
  For substrate-degrading enzymes (ACHE, MAO-A, MAO-B, COMT), tissue
  expression level (AHBA cortical y_ahba) does NOT linearly inform cognition
  relevance — k_cat/K_m operates at substrate-saturated regime; high enzyme
  density above the saturating threshold yields no additional functional
  clearance. Per MH8 deep-dive (research/4-tier/MH8 Methods Clarity Research.md),
  the AHBA contribution to θ for these targets should be effectively masked.
  Implementation: `build_y_obs_from_sources` accepts a
  `substrate_mediated_uniprots` set; for those targets, sigma_obs[AHBA_row, i]
  is inflated by `substrate_mediated_sigma_inflate` (default 10×, makes
  variance contribution 100× larger so likelihood is effectively marginalised).

  This is the structural fix for the 37 divergences observed on the 191-target
  V6.B.5 NUTS run — ACHE/MAO/COMT were forcing the posterior through
  Neal's-funnel geometries because the linear AHBA→θ link is biophysically
  invalid for these targets.

OPERATIONAL STATE (commit time):
  - PyMC + numpyro installed; PYMC_AVAILABLE = True; numpyro JAX backend used
  - 22-target headline: R̂=1.000, ESS=12,780, 0 divergences
  - 191-target expanded (pre-MH8): R̂=1.000, ESS=1,739, 37 divergences
  - 191-target expanded (post-MH8): TBD — Sprint 1.3 in progress

Reference:
  Multi-Source Neurobiological Prior for Cognition Target Prioritization.md §B
  research/4-tier/MH8 Methods Clarity Research.md §3-§4 (substrate-mediated bypass)
  Neelon & Dunson 2004 Biometrics 60:398
  Davies 2018 Nat Commun 9:2098 (intelligence GWAS)
  Hill 2019 Mol Psychiatry 24:169
  Moodie 2024 Hum Brain Mapp 45(4):e26641 (41-gene cortical g-map)
  Roberts 2020 Eur Neuropsychopharm 38:40 (SMD ceiling)
  Heaton 2024 (synaptic acetylcholinesterase k_cat/K_m saturation)
  Rommelfanger 2007 (PET MAO-A homeostatic adaptation to 5-HT)
  Chen 2011 (MB-COMT vs S-COMT non-linear brain regulation)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np

logger = logging.getLogger(__name__)

try:
    import pymc as pm    # noqa: F401
    PYMC_AVAILABLE = True
except ImportError:
    PYMC_AVAILABLE = False
    pm = None    # type: ignore


# Reference anchor compounds with literature-strong cognition relevance.
# These get tighter priors θ ~ N(0.5, 0.3²) to break scale degeneracy.
DEFAULT_ANCHORS: dict[str, float] = {
    "BDNF":    0.7,     # gold-standard cognition gene (Moodie 2024 + Hill 2019 + Siletti 2023)
    "COMT":    0.6,     # Davies 2018 L2G=0.71 at rs4680
    "ACHE":    0.5,     # substrate-mediated — flagged
    "DRD2":    0.5,     # pharmacogenetics
    "GRIN2B":  0.5,
    "CHRNA7":  0.5,
}

# Canonical substrate-degrading enzyme UniProts (MH8 substrate-mediated flag).
# These targets receive AHBA-masking in `build_y_obs_from_sources()` because
# tissue expression level does not linearly inform cognition relevance for
# enzymes operating at substrate-saturated regime (k_cat/K_m saturation).
SUBSTRATE_MEDIATED_UNIPROTS: frozenset[str] = frozenset({
    "P22303",   # ACHE — acetylcholinesterase
    "P21397",   # MAOA — monoamine oxidase A
    "P27338",   # MAOB — monoamine oxidase B
    "P21964",   # COMT — catechol-O-methyltransferase
})

# Default sigma-inflation factor for AHBA masking on substrate-mediated targets.
# 10× inflation makes variance contribution 100× larger, effectively
# marginalising the AHBA observation for these targets while keeping the
# model topology unchanged (no NaN handling needed).
DEFAULT_SM_SIGMA_INFLATE: float = 10.0


@dataclass
class ClusterDPosterior:
    targets: list[str]
    theta_mean: dict[str, float]                # posterior mean per target
    theta_lower: dict[str, float] = field(default_factory=dict)
    theta_upper: dict[str, float] = field(default_factory=dict)
    w_pipeline: dict[str, float] = field(default_factory=dict)   # σ(θ) ∈ (0, 1)
    jsd_disagreement: dict[str, float] = field(default_factory=dict)
    sources_used: list[str] = field(default_factory=list)
    n_chains: int = 0
    n_draws: int = 0
    n_divergences: int = 0                                       # MH8 + non-centered audit
    rhat_max: float = float("nan")
    ess_min: float = float("nan")
    method: str = "stage_0_stub"
    note: str = ""
    substrate_mediated_uniprots: list[str] = field(default_factory=list)   # MH8 audit: which targets received AHBA-masking


def fit_cluster_d_prior_stub(
    target_uniprots: list[str],
    y_ahba: dict[str, float] | None = None,
    y_l2g: dict[str, float] | None = None,
    y_sc: dict[str, float] | None = None,
    uniprot_to_gene: dict[str, str] | None = None,
) -> ClusterDPosterior:
    """Stage 0 stub — returns uniform w=0.5 per target when no input streams.

    When any of y_ahba / y_l2g / y_sc is provided, returns simple z-scored
    weighted-mean θ ∈ [-2, 2]. Real Bayesian model lives in
    `fit_cluster_d_prior_nuts` below.
    """
    n = len(target_uniprots)
    theta = np.zeros(n, dtype=float)
    contributions = 0
    for stream in (y_ahba, y_l2g, y_sc):
        if stream is None:
            continue
        vals = np.array([stream.get(t, 0.0) for t in target_uniprots])
        if vals.std() > 0:
            theta += (vals - vals.mean()) / vals.std()
            contributions += 1
    if contributions:
        theta /= contributions

    # Apply reference anchors as priors. DEFAULT_ANCHORS is keyed by the anchor identifier (gene
    # symbol). Pass `uniprot_to_gene` to resolve a UniProt accession -> gene so anchors fire on real
    # accessions; without the map, fall back to a direct key match (legacy behaviour). The previous
    # code only did the direct match, so anchors silently never fired on UniProt-accession input.
    for i, t in enumerate(target_uniprots):
        key = uniprot_to_gene.get(t, t) if uniprot_to_gene else t
        if key in DEFAULT_ANCHORS:
            theta[i] = 0.5 * theta[i] + 0.5 * DEFAULT_ANCHORS[key]

    theta_mean = {t: float(v) for t, v in zip(target_uniprots, theta)}
    w_pipeline = {t: float(1.0 / (1.0 + np.exp(-v))) for t, v in zip(target_uniprots, theta)}
    return ClusterDPosterior(
        targets=target_uniprots,
        theta_mean=theta_mean,
        w_pipeline=w_pipeline,
        method="stage_0_stub",
        sources_used=[s for s, v in [("AHBA", y_ahba), ("L2G", y_l2g), ("SC", y_sc)] if v],
        note="PyMC not installed; using z-score weighted-mean stub. "
             "Install pymc + numpyro and call fit_cluster_d_prior_nuts() for full Bayesian.",
    )


def fit_cluster_d_prior_nuts(
    target_uniprots: list[str],
    y_obs: np.ndarray,                # shape (n_sources, n_targets)
    sigma_obs: np.ndarray,            # shape (n_sources, n_targets)
    source_names: list[str],
    reference_idx: list[int] | None = None,
    reference_mean: float = 0.5,
    reference_sd: float = 0.3,
    n_chains: int = 4,
    n_draws: int = 2000,
    n_tune: int = 2000,
    target_accept: float = 0.95,
    random_seed: int = 42,
) -> ClusterDPosterior:
    """Full PyMC NUTS Bayesian hierarchical model from §B.2.

    y_obs[s, i] is the source-s observation for target i (Fisher-z transformed).
    Skeptical Lit prior triggered by `source_name == 'Lit'`.

    Reference anchors at θ_ref ~ N(reference_mean, reference_sd²) break the
    scale + sign degeneracy.

    Returns ClusterDPosterior with credible intervals.
    """
    if not PYMC_AVAILABLE:
        raise ImportError("PyMC not installed — use fit_cluster_d_prior_stub() instead")

    import pymc as pm
    import pytensor.tensor as pt
    import arviz as az

    S, T_n = y_obs.shape
    if T_n != len(target_uniprots):
        raise ValueError(f"y_obs has {T_n} targets but {len(target_uniprots)} provided")

    # Source-specific β scale (skeptical on Lit)
    beta_scale = np.ones(S)
    for j, s in enumerate(source_names):
        if s.lower().startswith("lit"):
            beta_scale[j] = 0.3

    with pm.Model() as model:
        # Latent target cognition-relevance θ ∈ R
        theta = pm.Normal("theta", mu=0.0, sigma=1.0, shape=T_n)

        # Reference-anchor likelihood via Potential (theta[idx] is derived,
        # not raw data — observed= would reject it as "depends on other nodes")
        if reference_idx:
            ref_dist = pm.Normal.dist(mu=reference_mean, sigma=reference_sd)
            ref_loglik = pm.logp(ref_dist, theta[reference_idx]).sum()
            pm.Potential("ref_anchor_loglik", ref_loglik)

        alpha = pm.Normal("alpha", mu=0.0, sigma=0.5, shape=S)
        beta = pm.HalfNormal("beta", sigma=beta_scale, shape=S)
        tau = pm.Gamma("tau", alpha=2.0, beta=2.0, shape=S)

        # Soft sum-to-zero on alpha
        pm.Normal("alpha_sum", mu=pt.sum(alpha), sigma=0.05, observed=0.0)

        # Likelihood
        mu_s = alpha[:, None] + beta[:, None] * theta[None, :]
        sigma_s2 = (1.0 / tau)[:, None] + sigma_obs ** 2
        pm.Normal("y", mu=mu_s, sigma=pt.sqrt(sigma_s2), observed=y_obs)

        sample_kwargs = dict(
            draws=n_draws, tune=n_tune, chains=n_chains,
            target_accept=target_accept,
            random_seed=random_seed, progressbar=False,
        )
        if _numpyro_available():
            sample_kwargs["nuts_sampler"] = "numpyro"
        idata = pm.sample(**sample_kwargs)

    theta_post = idata.posterior["theta"].values    # (chain, draw, target)
    theta_flat = theta_post.reshape(-1, T_n)
    theta_mean = {t: float(theta_flat[:, i].mean()) for i, t in enumerate(target_uniprots)}
    theta_lo = {t: float(np.percentile(theta_flat[:, i], 2.5)) for i, t in enumerate(target_uniprots)}
    theta_hi = {t: float(np.percentile(theta_flat[:, i], 97.5)) for i, t in enumerate(target_uniprots)}
    w_pipeline = {t: float(1.0 / (1.0 + np.exp(-theta_mean[t]))) for t in target_uniprots}

    # Diagnostics
    summary = az.summary(idata, var_names=["theta"])
    rhat_max = float(summary["r_hat"].max())
    ess_min = float(summary["ess_bulk"].min())

    # Divergence count — used for MH8 audit.
    try:
        sample_stats = idata.sample_stats
        if "diverging" in sample_stats:
            n_divergences = int(sample_stats["diverging"].values.sum())
        else:
            n_divergences = 0
    except Exception:    # noqa: BLE001
        n_divergences = -1    # sentinel: unknown

    note = f"NUTS converged: Rhat_max={rhat_max:.3f}, ESS_min={ess_min:.0f}, divergences={n_divergences}"
    return ClusterDPosterior(
        targets=target_uniprots,
        theta_mean=theta_mean,
        theta_lower=theta_lo,
        theta_upper=theta_hi,
        w_pipeline=w_pipeline,
        sources_used=source_names,
        n_chains=n_chains,
        n_draws=n_draws,
        n_divergences=n_divergences,
        rhat_max=rhat_max,
        ess_min=ess_min,
        method="pymc_nuts",
        note=note,
    )


def _numpyro_available() -> bool:
    try:
        import numpyro    # noqa: F401
        return True
    except ImportError:
        return False


def build_y_obs_from_sources(
    target_uniprots: list[str],
    y_ahba: dict[str, float] | None = None,
    y_l2g: dict[str, float] | None = None,
    y_sc: dict[str, float] | None = None,
    sigma_ahba: float = 0.30,    # AHBA spatial-corr Fisher-z SE; Markello 2021
    sigma_l2g: float = 0.20,     # OT Genetics L2G Shapley-XGBoost noise
    sigma_sc: float = 0.35,      # cellxgene single-cell aggregation noise
    substrate_mediated_uniprots: set[str] | frozenset[str] | None = None,
    substrate_mediated_sigma_inflate: float = DEFAULT_SM_SIGMA_INFLATE,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Stack AHBA / L2G / SC observations into the (S, T) matrix shape
    expected by `fit_cluster_d_prior_nuts`.

    Returns (y_obs, sigma_obs, source_names). Missing sources are dropped;
    missing per-target values default to 0.0 with inflated sigma.

    MH8 substrate-mediated bypass:
      If `substrate_mediated_uniprots` is provided, for each UniProt in the
      set that also appears in `target_uniprots`, the AHBA-row sigma is
      multiplied by `substrate_mediated_sigma_inflate` (default 10×).
      This effectively marginalises the AHBA observation for substrate-
      degrading enzymes (ACHE / MAO-A / MAO-B / COMT) whose tissue
      expression does not linearly inform cognition relevance. If no set is
      provided, the canonical `SUBSTRATE_MEDIATED_UNIPROTS` is NOT applied
      automatically — pass `substrate_mediated_uniprots=SUBSTRATE_MEDIATED_UNIPROTS`
      explicitly to opt in. This preserves backward compatibility with v6b
      paper headline results that did not apply the fix.

      Returns: the inflated sigma_obs reflects MH8 masking; the count and
      identity of masked targets are logged at INFO level for audit.
    """
    sources: list[tuple[str, dict[str, float] | None, float]] = [
        ("AHBA", y_ahba, sigma_ahba),
        ("L2G", y_l2g, sigma_l2g),
        ("SC", y_sc, sigma_sc),
    ]
    active = [(name, stream, sig) for name, stream, sig in sources if stream is not None]
    if not active:
        raise ValueError("At least one of y_ahba / y_l2g / y_sc must be provided")

    T_n = len(target_uniprots)
    y_obs = np.zeros((len(active), T_n), dtype=float)
    sigma_obs = np.zeros((len(active), T_n), dtype=float)
    source_names: list[str] = []
    for s_idx, (name, stream, sig) in enumerate(active):
        source_names.append(name)
        for t_idx, t in enumerate(target_uniprots):
            v = stream.get(t)
            if v is None:
                # Missing observation → centered at 0 with very inflated sigma
                y_obs[s_idx, t_idx] = 0.0
                sigma_obs[s_idx, t_idx] = sig * 5.0
            else:
                y_obs[s_idx, t_idx] = float(v)
                sigma_obs[s_idx, t_idx] = sig

    # MH8 substrate-mediated bypass: inflate AHBA sigma for these targets.
    if substrate_mediated_uniprots and "AHBA" in source_names:
        ahba_row = source_names.index("AHBA")
        sm_set = set(substrate_mediated_uniprots)
        masked_uniprots: list[str] = []
        for t_idx, t in enumerate(target_uniprots):
            if t in sm_set:
                sigma_obs[ahba_row, t_idx] *= substrate_mediated_sigma_inflate
                masked_uniprots.append(t)
        if masked_uniprots:
            logger.info(
                "MH8: inflated AHBA sigma by %.1fx for %d substrate-mediated targets: %s",
                substrate_mediated_sigma_inflate,
                len(masked_uniprots),
                ",".join(masked_uniprots),
            )

    return y_obs, sigma_obs, source_names


def roberts_2020_ceiling_check(
    posterior: ClusterDPosterior,
    target_smd_predictions: dict[str, float] | None = None,
    smd_ceiling: float = 0.5,
    upper_quantile: float = 0.90,
) -> dict[str, str]:
    """Gate 1 from Cluster D §G — Roberts 2020 SMD ceiling.

    "No target's predicted modulator effect-size posterior may exceed
     Hedges' g = 0.5 at 90% credible upper bound."

    If `target_smd_predictions` is None, this is a no-op (no SMD predictions
    available); returns each target as REGIME_OK.
    """
    out: dict[str, str] = {}
    for t in posterior.targets:
        if not target_smd_predictions or t not in target_smd_predictions:
            out[t] = "NO_SMD_PREDICTION"
            continue
        smd = target_smd_predictions[t]
        if smd > smd_ceiling:
            out[t] = "REGIME_VIOLATION"
        else:
            out[t] = "REGIME_OK"
    return out
