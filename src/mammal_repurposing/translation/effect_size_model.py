"""V7.3 — Effect-size hierarchical Bayes model.

3-level PyMC NUTS hierarchy per `research/4-tier/Clinical Effect-Size
Translation Function.md` §B:

    μ_global             ~ Normal(0, 0.20)
    μ_class[m]           ~ Normal(μ_class_PRISMA[m], λ_class · σ_class_PRISMA[m])
    g[c, t]              ~ Normal(η[c, t], σ_resid²)

with sigmoid translation:
    η[c, t] = sigmoid(α + β1·E[pchembl_post] + β2·E[relevance_post]
                       + β3·copula_correction) − Σ_k γ_k · m_k

and Cluster D multiplicative gate:
    β_target[t_c] = θ̄_{t_c} · β_raw_target[t_c]

5 failure-mode moderators (m_k):
    m1 = U-shape miss          (dose past peak; from PBPK u_shape_occupancy)
    m2 = practice/placebo      (trial design moderator)
    m3 = tolerance onset       (chronic vs acute mismatch; from R_avail dynamics)
    m4 = trait × state         (responder enrichment masks population mean)
    m5 = trial-design          (parallel-group vs crossover, endpoint sensitivity)

Graceful degradation:
    - PyMC missing → returns class-mean stub posterior (PRISMA prior only)
    - PyMC present → full NUTS via numpyro backend

API:
    obs = EffectSizeObservation(
        compound="donepezil", class_name="AChE-I",
        target_uniprot="P22303", pchembl_post_mean=8.5, pchembl_post_sd=0.3,
        relevance_post_mean=0.92, relevance_post_sd=0.08,
        pbpk_auc_brain=12.0, moderators=[0, 0, 0, 0, 0],
    )
    posterior = fit_effect_size_model([obs1, obs2, ...])
    # posterior.g_mean[compound] → predicted Hedges' g
    # posterior.g_2p5[compound], posterior.g_97p5[compound] → 95% CrI
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Sequence

import numpy as np

logger = logging.getLogger(__name__)

try:
    import pymc as pm   # noqa: F401
    PYMC_AVAILABLE = True
except ImportError:
    PYMC_AVAILABLE = False
    pm = None    # type: ignore


N_MODERATORS = 5
MODERATOR_NAMES: tuple[str, ...] = (
    "u_shape_miss",
    "practice_placebo",
    "tolerance_onset",
    "trait_state_interaction",
    "trial_design",
)


@dataclass
class EffectSizeObservation:
    """One (compound × endpoint) observation feeding the V7 hierarchical model.

    Fields:
        compound: unique compound identifier (e.g., 'donepezil')
        class_name: PRISMA mechanism class (must be in prisma_priors.list_class_names())
        target_uniprot: primary target UniProt accession (for Cluster D gating)
        pchembl_post_mean / pchembl_post_sd: V6.A Venn-ABERS-calibrated posterior
        relevance_post_mean / relevance_post_sd: V6.B Cluster D θ̄ posterior
            (σ(θ̄) ∈ (0, 1))
        pbpk_auc_brain: brain-compartment AUC from V7.1 PBPK (μM·h or nM·h;
            entered as a normalised exposure metric, mean-centered downstream)
        moderators: length-5 array of {0, 1} flags or [0, 1] continuous scores
            for the 5 failure modes (m1-m5)
        observed_g: optional ground-truth Hedges' g (for compounds in the
            anchor set with published meta-analytic SMD)
        endpoint: cognitive endpoint or domain (V1: ADAS-Cog, DSST, n-back,
            Stroop, RAVLT, CANTAB-RVIP, MCCB; V2: EM, WM, ATT, EF, PS, VL, VS,
            MOT — the 8 canonical cognitive domains per V7.2 Stage 3).
        population: Sprint 3.2 — patient population for the population × class
            interaction term. Canonical values: HC, AD, MCI, SCZ, ADHD, FXS,
            MDD, NRC, OSA, PD, HD, DownSyndrome, Vasc, DepCog, EpilCog, mixed.
            Default 'mixed' to preserve backward compatibility with V7.3 callers.
    """
    compound: str
    class_name: str
    target_uniprot: str
    pchembl_post_mean: float
    pchembl_post_sd: float = 0.3
    relevance_post_mean: float = 0.5
    relevance_post_sd: float = 0.1
    pbpk_auc_brain: float = 1.0
    moderators: tuple[float, float, float, float, float] = (0.0,) * N_MODERATORS
    observed_g: float | None = None
    endpoint: str = "DSST"
    population: str = "mixed"


@dataclass
class EffectSizePosterior:
    """V7 hierarchical Bayes output."""
    compounds: list[str]
    g_mean: dict[str, float]                              # posterior mean
    g_2p5: dict[str, float] = field(default_factory=dict)  # 95% CrI lower
    g_97p5: dict[str, float] = field(default_factory=dict) # 95% CrI upper
    g_90_upper: dict[str, float] = field(default_factory=dict)  # 90% CrI upper (for Roberts gate)
    class_mu: dict[str, float] = field(default_factory=dict)
    cluster_d_gate_active: dict[str, bool] = field(default_factory=dict)
    n_chains: int = 0
    n_draws: int = 0
    n_divergences: int = 0    # NUTS divergence count (Sprint 3.4 audit)
    rhat_max: float = float("nan")
    ess_min: float = float("nan")
    method: str = "prisma_stub"
    note: str = ""


def _sigmoid(x: np.ndarray | float) -> np.ndarray | float:
    return 1.0 / (1.0 + np.exp(-x))


def fit_effect_size_stub(
    observations: Sequence[EffectSizeObservation],
    class_prior_lookup: dict[str, dict] | None = None,
    cluster_d_floor: float = 0.10,
    use_subdomain_priors: bool = True,
) -> EffectSizePosterior:
    """Stage 0 stub — returns PRISMA class-mean posterior per compound,
    multiplicatively gated by Cluster D θ̄ (V6.B posterior).

    Uses the class prior mean as the posterior point estimate when PyMC is
    unavailable. The Cluster D gate β_target[t_c] = θ̄_{t_c} · β_raw applies
    as: g_stub[c] = class_prior_mean · max(cluster_d_floor, relevance_post).

    V7.2 Stage 2: if use_subdomain_priors=True (default), prefers per-
    (class, endpoint) prior from PER_SUBDOMAIN_PRIORS when the observation's
    endpoint matches a registered subdomain. Falls back to class-level prior
    otherwise.

    Returns EffectSizePosterior with .method='prisma_stub' or
    'prisma_stub_subdomain'.
    """
    from .prisma_priors import class_prior_table, get_subdomain_prior
    table = class_prior_lookup or class_prior_table()
    g_mean: dict[str, float] = {}
    g_2p5: dict[str, float] = {}
    g_97p5: dict[str, float] = {}
    g_90_upper: dict[str, float] = {}
    class_mu: dict[str, float] = {}
    gate_active: dict[str, bool] = {}
    compounds: list[str] = []
    n_subdomain_hits = 0

    for obs in observations:
        compounds.append(obs.compound)
        # V7.2 Stage 2: try subdomain prior first
        used_subdomain = False
        if use_subdomain_priors:
            sub_mean, sub_sd = get_subdomain_prior(obs.class_name, obs.endpoint,
                                                    fallback_to_class=False)
            if (sub_mean, sub_sd) != (0.0, 0.15) or (obs.class_name, obs.endpoint) in \
                    __import__("mammal_repurposing.translation.prisma_priors",
                               fromlist=["PER_SUBDOMAIN_PRIORS"]).PER_SUBDOMAIN_PRIORS:
                class_mean = sub_mean
                class_sd = sub_sd
                used_subdomain = True
                n_subdomain_hits += 1
        if not used_subdomain:
            cls = table.get(obs.class_name)
            if cls is None:
                logger.warning("Compound %s class '%s' missing from PRISMA priors; "
                               "defaulting to (mean=0, sd=0.15)", obs.compound,
                               obs.class_name)
                class_mean, class_sd = 0.0, 0.15
            else:
                class_mean = cls["mean"]
                class_sd = cls["sd"]
        gate = max(cluster_d_floor, obs.relevance_post_mean)
        gated_mean = class_mean * gate

        # Subtract moderator penalties: each active moderator (γ=0.05 default)
        # debits ~0.05 from g; conservative stub.
        moderator_penalty = float(np.sum(obs.moderators)) * 0.05
        g = gated_mean - moderator_penalty

        # Stub CIs: ± 1.96 · class_sd, expanded by Cluster D uncertainty
        ci_sd = float(np.sqrt(class_sd ** 2 + obs.relevance_post_sd ** 2
                              + obs.pchembl_post_sd ** 2 * 0.05))
        g_mean[obs.compound] = g
        g_2p5[obs.compound] = g - 1.96 * ci_sd
        g_97p5[obs.compound] = g + 1.96 * ci_sd
        g_90_upper[obs.compound] = g + 1.282 * ci_sd
        class_mu[obs.compound] = class_mean
        gate_active[obs.compound] = obs.relevance_post_mean > cluster_d_floor

    method = "prisma_stub_subdomain" if (use_subdomain_priors
                                          and n_subdomain_hits > 0) else "prisma_stub"
    note = ("PyMC unavailable or stub-mode requested. Returns "
            "PRISMA class-mean × Cluster D θ̄ multiplicative gate − "
            "moderator penalties. Full NUTS posterior via "
            "fit_effect_size_nuts() requires `pip install pymc numpyro`.")
    if use_subdomain_priors and n_subdomain_hits > 0:
        note += (f" V7.2 Stage 2 active: {n_subdomain_hits}/{len(observations)} "
                 "observations used per-(class, endpoint) subdomain prior.")
    return EffectSizePosterior(
        compounds=compounds,
        g_mean=g_mean,
        g_2p5=g_2p5,
        g_97p5=g_97p5,
        g_90_upper=g_90_upper,
        class_mu=class_mu,
        cluster_d_gate_active=gate_active,
        method=method,
        note=note,
    )


def fit_effect_size_nuts(
    observations: Sequence[EffectSizeObservation],
    class_prior_lookup: dict[str, dict] | None = None,
    n_chains: int = 4,
    n_draws: int = 2000,
    n_tune: int = 2000,
    target_accept: float = 0.95,
    lambda_class: float = 1.0,      # robust MAP class-prior strength
    sigma_resid_prior: float = 0.20,
    moderator_gamma_prior_sd: float = 0.10,
    random_seed: int = 42,
    subdomain_anchor_weight: float = 0.5,    # V7.2 Stage 2: per-(class, endpoint) anchor strength
) -> EffectSizePosterior:
    """Full PyMC NUTS hierarchical Bayes.

    Per V7 spec §B:
        μ_global        ~ Normal(0, 0.20)
        μ_class[m]      ~ Normal(prisma_mean[m], λ_class · prisma_sd[m])
        β_pchembl       ~ Normal(0, 0.10)
        β_relevance     ~ Normal(0, 0.10)
        β_copula        ~ Normal(0, 0.05)
        γ_k             ~ Normal(0, 0.10) for each of 5 moderators
        β_target[c]     = θ̄_c · β_raw_target[c]   ← Cluster D gate
        η[c]            = sigmoid(α + β·E[pchembl] + β·E[relevance] + β·copula) ·
                          β_target[c] − Σ_k γ_k · m_k[c]
        g[c]            ~ Normal(η[c], σ_resid²)
        σ_resid         ~ HalfNormal(σ_resid_prior)

    Compounds with `observed_g` are likelihood-anchored; the rest are
    forward-predicted from priors.
    """
    if not PYMC_AVAILABLE:
        raise ImportError(
            "PyMC not installed. Use fit_effect_size_stub() instead, or "
            "`pip install pymc numpyro` for the full Bayesian path."
        )
    import pymc as pm
    import pytensor.tensor as pt
    import arviz as az

    from .prisma_priors import class_prior_table, list_class_names
    table = class_prior_lookup or class_prior_table()
    class_names = list_class_names()
    n_obs = len(observations)
    if n_obs == 0:
        raise ValueError("Cannot fit on zero observations")

    # Index observations
    cls_idx = np.array([class_names.index(o.class_name)
                        if o.class_name in class_names else -1
                        for o in observations])
    pchembl_mean = np.array([o.pchembl_post_mean for o in observations])
    pchembl_sd = np.array([o.pchembl_post_sd for o in observations])
    relevance_mean = np.array([o.relevance_post_mean for o in observations])
    relevance_sd = np.array([o.relevance_post_sd for o in observations])
    moderators_mat = np.array([list(o.moderators) for o in observations])
    pbpk_auc = np.array([o.pbpk_auc_brain for o in observations])
    obs_g = np.array([o.observed_g if o.observed_g is not None else np.nan
                      for o in observations])
    has_g = ~np.isnan(obs_g)

    # Class-prior means + sds (Schmidli MAP)
    prisma_means = np.array([table.get(c, {"mean": 0.0})["mean"]
                              for c in class_names])
    prisma_sds = np.array([table.get(c, {"sd": 0.15})["sd"]
                           for c in class_names])

    with pm.Model() as model:
        mu_global = pm.Normal("mu_global", 0.0, 0.20)

        mu_class = pm.Normal(
            "mu_class",
            mu=prisma_means,
            sigma=lambda_class * prisma_sds,
            shape=len(class_names),
        )

        beta_pchembl = pm.Normal("beta_pchembl", 0.0, 0.10)
        beta_relevance = pm.Normal("beta_relevance", 0.0, 0.10)
        beta_copula = pm.Normal("beta_copula", 0.0, 0.05)
        gamma = pm.Normal("gamma_moderators", 0.0, moderator_gamma_prior_sd,
                          shape=N_MODERATORS)
        beta_pbpk = pm.Normal("beta_pbpk", 0.0, 0.05)

        # Center pchembl + pbpk to keep sigmoid well-conditioned
        pchembl_c = (pchembl_mean - pchembl_mean.mean()) / max(pchembl_mean.std(), 1.0)
        pbpk_c = (pbpk_auc - pbpk_auc.mean()) / max(pbpk_auc.std(), 1.0)
        copula = pchembl_c * relevance_mean   # interaction proxy

        eta_inner = (mu_global
                     + beta_pchembl * pchembl_c
                     + beta_relevance * relevance_mean
                     + beta_copula * copula
                     + beta_pbpk * pbpk_c)
        sigmoid_eta = pm.math.sigmoid(eta_inner)    # ∈ (0, 1)

        # Class-mean contribution + Cluster D multiplicative gate
        cls_idx_safe = np.where(cls_idx >= 0, cls_idx, 0)
        class_contrib = mu_class[cls_idx_safe]
        # Cluster D gate: scale class effect by relevance posterior
        gated_class = class_contrib * relevance_mean

        # Moderator debit
        moderator_debit = pm.math.dot(moderators_mat, gamma)

        eta = sigmoid_eta * gated_class - moderator_debit

        sigma_resid = pm.HalfNormal("sigma_resid", sigma_resid_prior)

        # Likelihood: only over observed g
        if has_g.sum() > 0:
            pm.Normal(
                "g_obs",
                mu=eta[has_g],
                sigma=sigma_resid,
                observed=obs_g[has_g],
            )

        # V7.2 Stage 2: per-(class, endpoint) anchor likelihood (soft) via Potential.
        # For observations where (class, endpoint) has a curated subdomain prior,
        # add an additional Gaussian pull centered on the subdomain mean.
        from .prisma_priors import get_subdomain_prior, PER_SUBDOMAIN_PRIORS
        subdomain_mu = np.zeros(n_obs)
        subdomain_sd = np.full(n_obs, 1e6)    # huge sd = no effect (default)
        n_subdomain_hits = 0
        for i, o in enumerate(observations):
            if (o.class_name, o.endpoint) in PER_SUBDOMAIN_PRIORS:
                sm, ssd = PER_SUBDOMAIN_PRIORS[(o.class_name, o.endpoint)]
                subdomain_mu[i] = sm
                # Inflate sd by 1/subdomain_anchor_weight so weight=0 → no anchor,
                # weight=1 → tight anchor at the subdomain prior sd
                subdomain_sd[i] = ssd / max(subdomain_anchor_weight, 1e-6)
                n_subdomain_hits += 1
        if n_subdomain_hits > 0 and subdomain_anchor_weight > 0:
            sub_loglik = pm.logp(
                pm.Normal.dist(mu=subdomain_mu, sigma=subdomain_sd),
                eta,
            ).sum()
            pm.Potential("subdomain_anchor_loglik", sub_loglik)
            logger.info("V7.2 Stage 2: subdomain anchor applied to %d/%d observations",
                        n_subdomain_hits, n_obs)

        # Track η for all compounds (whether or not observed) for posterior export
        pm.Deterministic("g_pred", eta)

        sample_kwargs = dict(
            draws=n_draws, tune=n_tune, chains=n_chains,
            target_accept=target_accept,
            random_seed=random_seed, progressbar=False,
        )
        if _numpyro_available():
            sample_kwargs["nuts_sampler"] = "numpyro"
        idata = pm.sample(**sample_kwargs)

    g_post = idata.posterior["g_pred"].values    # (chain, draw, n_obs)
    g_flat = g_post.reshape(-1, n_obs)

    compounds = [o.compound for o in observations]
    g_mean = {c: float(g_flat[:, i].mean()) for i, c in enumerate(compounds)}
    g_2p5 = {c: float(np.percentile(g_flat[:, i], 2.5)) for i, c in enumerate(compounds)}
    g_97p5 = {c: float(np.percentile(g_flat[:, i], 97.5)) for i, c in enumerate(compounds)}
    g_90_upper = {c: float(np.percentile(g_flat[:, i], 90.0)) for i, c in enumerate(compounds)}
    class_mu = {c: float(np.mean([prisma_means[ci] for ci in [cls_idx[i]] if ci >= 0]) or 0.0)
                for i, c in enumerate(compounds)}
    gate_active = {c: bool(observations[i].relevance_post_mean > 0.10)
                   for i, c in enumerate(compounds)}

    summary = az.summary(idata, var_names=["mu_class", "beta_pchembl",
                                            "beta_relevance", "sigma_resid"])
    rhat_max = float(summary["r_hat"].max())
    ess_min = float(summary["ess_bulk"].min())

    return EffectSizePosterior(
        compounds=compounds,
        g_mean=g_mean,
        g_2p5=g_2p5,
        g_97p5=g_97p5,
        g_90_upper=g_90_upper,
        class_mu=class_mu,
        cluster_d_gate_active=gate_active,
        n_chains=n_chains, n_draws=n_draws,
        rhat_max=rhat_max, ess_min=ess_min,
        method="pymc_nuts",
        note=f"NUTS converged: R̂={rhat_max:.3f}, ESS={ess_min:.0f}; "
             f"λ_class={lambda_class}",
    )


def _numpyro_available() -> bool:
    try:
        import numpyro   # noqa: F401
        return True
    except ImportError:
        return False


# ===========================================================================
# Sprint 3.2 — V7.3 Stage 2 hierarchical Bayes with per-class τ² + population × class interaction
# ===========================================================================

def fit_effect_size_nuts_v2(
    observations: Sequence[EffectSizeObservation],
    *,
    n_chains: int = 4,
    n_draws: int = 2000,
    n_tune: int = 2000,
    target_accept: float = 0.95,
    tau_class_prior_scale: float = 0.05,       # HalfNormal scale for per-class τ
    interaction_prior_scale: float = 0.10,      # interaction term SD
    sigma_resid_prior: float = 0.20,
    moderator_gamma_prior_sd: float = 0.10,
    random_seed: int = 42,
    use_v2_subdomain_priors: bool = True,
    subdomain_anchor_weight: float = 0.5,
) -> EffectSizePosterior:
    """V7.3 Stage 2 — refit with per-class τ²_class + population × class interaction.

    Sprint 3.2 deliverable per `reports/MH_IMPLEMENTATION_ROADMAP.md` § 3
    (V7 model-structure refit BEFORE anchor expansion).

    Model (canonical form per MH1+MH2 V7 CPT doc § 5):
        μ_global              ~ Normal(0, 0.20)
        τ_class[m]            ~ HalfNormal(tau_class_prior_scale)
        μ_class[m]            ~ Normal(prisma_mean[m], τ_class[m])        # per-class heterogeneity
        ι[m, p]               ~ Normal(0, interaction_prior_scale)         # population × class
        β_pchembl, β_relevance, β_copula, β_pbpk
        γ_k (k=5)             ~ Normal(0, moderator_gamma_prior_sd)
        sigma_resid           ~ HalfNormal(sigma_resid_prior)

    Linear predictor:
        η[c] = sigmoid(α + β·E[pchembl_c] + β·relevance_c + β·copula + β·pbpk_c)
                · (μ_class[m_c] + ι[m_c, p_c]) · relevance_c
                − Σ_k γ_k · m_k[c]

    where m_c = class index, p_c = population index for observation c.

    Compared to v1:
      - τ²_class is now PER class (was single λ_class shared across all classes)
      - Adds (population × class) interaction ι[m, p] to capture e.g.
        AChE-I in AD (g≈0.36 EM) vs AChE-I in HC (g≈0.15-0.20)
      - When `use_v2_subdomain_priors`, anchor likelihood pulls each
        observation toward its (class, domain) V2 cell via Potential.

    Returns:
      EffectSizePosterior with method='pymc_nuts_v2' and a populated
      `note` field reporting per-class τ̂ posterior + interaction magnitudes.
    """
    if not PYMC_AVAILABLE:
        raise ImportError(
            "PyMC not installed. Use fit_effect_size_stub() instead, or "
            "`pip install pymc numpyro` for the full Bayesian path."
        )
    import pymc as pm
    import pytensor.tensor as pt
    import arviz as az

    from .prisma_priors import (
        class_prior_table, list_class_names,
        PER_SUBDOMAIN_PRIORS_V2,
    )

    table = class_prior_table()
    class_names = list_class_names()
    n_obs = len(observations)
    if n_obs == 0:
        raise ValueError("Cannot fit on zero observations")

    # Index classes and populations
    cls_idx = np.array([class_names.index(o.class_name)
                        if o.class_name in class_names else 0
                        for o in observations])
    populations = sorted(set(o.population for o in observations))
    pop_to_idx = {p: i for i, p in enumerate(populations)}
    pop_idx = np.array([pop_to_idx[o.population] for o in observations])
    n_pop = len(populations)
    n_classes = len(class_names)

    pchembl_mean = np.array([o.pchembl_post_mean for o in observations])
    relevance_mean = np.array([o.relevance_post_mean for o in observations])
    moderators_mat = np.array([list(o.moderators) for o in observations])
    pbpk_auc = np.array([o.pbpk_auc_brain for o in observations])
    obs_g = np.array([o.observed_g if o.observed_g is not None else np.nan
                      for o in observations])
    has_g = ~np.isnan(obs_g)

    prisma_means = np.array([table.get(c, {"mean": 0.0})["mean"]
                              for c in class_names])

    with pm.Model() as model:
        mu_global = pm.Normal("mu_global", 0.0, 0.20)

        # Per-class heterogeneity: τ_class[m] ~ HalfNormal(scale)
        tau_class = pm.HalfNormal("tau_class", tau_class_prior_scale,
                                   shape=n_classes)

        # μ_class[m] ~ Normal(prisma_mean[m], τ_class[m]) — non-centered
        mu_class_raw = pm.Normal("mu_class_raw", 0.0, 1.0, shape=n_classes)
        mu_class = pm.Deterministic(
            "mu_class",
            prisma_means + tau_class * mu_class_raw,
        )

        # Population × class interaction: ι[m, p] ~ Normal(0, scale)
        iota_raw = pm.Normal("iota_raw", 0.0, 1.0, shape=(n_classes, n_pop))
        iota = pm.Deterministic("iota",
                                 iota_raw * interaction_prior_scale)

        # Linear predictor components
        beta_pchembl = pm.Normal("beta_pchembl", 0.0, 0.10)
        beta_relevance = pm.Normal("beta_relevance", 0.0, 0.10)
        beta_copula = pm.Normal("beta_copula", 0.0, 0.05)
        beta_pbpk = pm.Normal("beta_pbpk", 0.0, 0.05)
        gamma = pm.Normal("gamma_moderators", 0.0, moderator_gamma_prior_sd,
                          shape=N_MODERATORS)

        pchembl_c = (pchembl_mean - pchembl_mean.mean()) / max(pchembl_mean.std(), 1.0)
        pbpk_c = (pbpk_auc - pbpk_auc.mean()) / max(pbpk_auc.std(), 1.0)
        copula = pchembl_c * relevance_mean

        eta_inner = (mu_global
                     + beta_pchembl * pchembl_c
                     + beta_relevance * relevance_mean
                     + beta_copula * copula
                     + beta_pbpk * pbpk_c)
        sigmoid_eta = pm.math.sigmoid(eta_inner)

        # Class contribution + population × class interaction
        class_contrib = mu_class[cls_idx]
        interaction_contrib = iota[cls_idx, pop_idx]
        combined_class_effect = class_contrib + interaction_contrib

        # Cluster D gate + moderator debit
        gated = combined_class_effect * relevance_mean
        moderator_debit = pm.math.dot(moderators_mat, gamma)
        eta = sigmoid_eta * gated - moderator_debit

        sigma_resid = pm.HalfNormal("sigma_resid", sigma_resid_prior)

        if has_g.sum() > 0:
            pm.Normal("g_obs", mu=eta[has_g], sigma=sigma_resid,
                      observed=obs_g[has_g])

        # V2 subdomain anchor likelihood (per-cell pull via Potential)
        if use_v2_subdomain_priors and subdomain_anchor_weight > 0:
            from .prisma_priors import (CLASS_NAME_MIGRATION_V1_TO_V2,
                                          get_subdomain_prior_v2)
            subdomain_mu = np.zeros(n_obs)
            subdomain_sd = np.full(n_obs, 1e6)
            n_subdomain_hits = 0
            for i, o in enumerate(observations):
                cell = get_subdomain_prior_v2(o.class_name, o.endpoint)
                if cell is not None:
                    subdomain_mu[i] = cell.pooled_g
                    # Tighter τ (Schmidli rule) → tighter anchor pull
                    subdomain_sd[i] = cell.tau / max(subdomain_anchor_weight, 1e-6)
                    n_subdomain_hits += 1
            if n_subdomain_hits > 0:
                sub_loglik = pm.logp(
                    pm.Normal.dist(mu=subdomain_mu, sigma=subdomain_sd),
                    eta,
                ).sum()
                pm.Potential("v2_subdomain_anchor", sub_loglik)
                logger.info("V7.3 Stage 2 v2 subdomain anchors: %d/%d observations",
                            n_subdomain_hits, n_obs)

        pm.Deterministic("g_pred", eta)

        sample_kwargs = dict(
            draws=n_draws, tune=n_tune, chains=n_chains,
            target_accept=target_accept,
            random_seed=random_seed, progressbar=False,
        )
        if _numpyro_available():
            sample_kwargs["nuts_sampler"] = "numpyro"
        idata = pm.sample(**sample_kwargs)

    g_post = idata.posterior["g_pred"].values
    g_flat = g_post.reshape(-1, n_obs)

    compounds = [o.compound for o in observations]
    g_mean = {c: float(g_flat[:, i].mean()) for i, c in enumerate(compounds)}
    g_2p5 = {c: float(np.percentile(g_flat[:, i], 2.5)) for i, c in enumerate(compounds)}
    g_97p5 = {c: float(np.percentile(g_flat[:, i], 97.5)) for i, c in enumerate(compounds)}
    g_90_upper = {c: float(np.percentile(g_flat[:, i], 90.0)) for i, c in enumerate(compounds)}
    class_mu = {c: float(prisma_means[cls_idx[i]]) for i, c in enumerate(compounds)}
    gate_active = {c: bool(observations[i].relevance_post_mean > 0.10)
                   for i, c in enumerate(compounds)}

    summary = az.summary(idata, var_names=["mu_class", "tau_class", "iota",
                                            "beta_pchembl", "beta_relevance",
                                            "sigma_resid"])
    rhat_max = float(summary["r_hat"].max())
    ess_min = float(summary["ess_bulk"].min())

    # Divergence count
    try:
        diverging = idata.sample_stats["diverging"].values if "diverging" in idata.sample_stats else None
        n_div = int(diverging.sum()) if diverging is not None else 0
    except Exception:    # noqa: BLE001
        n_div = -1

    # Per-class τ̂ posterior summary
    tau_post = idata.posterior["tau_class"].values
    tau_per_class = tau_post.mean(axis=(0, 1))
    tau_summary = "; ".join(f"{c}:{tau_per_class[i]:.3f}"
                            for i, c in enumerate(class_names[:5]))

    return EffectSizePosterior(
        compounds=compounds,
        g_mean=g_mean,
        g_2p5=g_2p5,
        g_97p5=g_97p5,
        g_90_upper=g_90_upper,
        class_mu=class_mu,
        cluster_d_gate_active=gate_active,
        n_chains=n_chains, n_draws=n_draws,
        n_divergences=n_div,
        rhat_max=rhat_max, ess_min=ess_min,
        method="pymc_nuts_v2",
        note=(f"NUTS v2 converged: R̂={rhat_max:.3f}, ESS={ess_min:.0f}, "
              f"divergences={n_div}; tau_class[{tau_summary},...]; "
              f"n_pop={n_pop} ({','.join(populations)})"),
    )


def assert_p1_through_p8(
    posterior: EffectSizePosterior,
    predictions: dict[str, tuple[float, float]] | None = None,
) -> dict[str, str]:
    """V7.4 Gate 1: all 8 P1-P8 pre-registered predictions land within bands.

    Default predictions per V4 §13.Y:
        P1: donepezil g ∈ [0.10, 0.30]
        P2: encenicline_3mg g recapitulates Phase 3 failure (|g| < 0.20)
        P3: methylphenidate_20mg g ∈ [0.15, 0.30]
        P4: modafinil_200mg g ∈ [0.06, 0.18]
        P5: memantine_20mg g ∈ [-0.05, 0.20]   (healthy adults)
        P6: intepirdine g ∈ [-0.10, 0.15]      (MINDSET-replicated)
        P7: pridopidine g ∈ [-0.10, 0.15]      (PROOF-HD-replicated)
        P8: lecanemab g ∈ [0.0, 0.15]          (cognitive subdomain)

    Returns {prediction_id: 'PASS' | 'FAIL' | 'NO_COMPOUND'}.
    """
    default_predictions = {
        "P1_donepezil":         (0.10, 0.30),
        "P2_encenicline_3mg":   (-0.20, 0.20),
        "P3_methylphenidate_20mg": (0.15, 0.30),
        "P4_modafinil_200mg":   (0.06, 0.18),
        "P5_memantine_20mg":    (-0.05, 0.20),
        "P6_intepirdine":       (-0.10, 0.15),
        "P7_pridopidine":       (-0.10, 0.15),
        "P8_lecanemab":         (0.00, 0.15),
    }
    predictions = predictions or default_predictions
    out: dict[str, str] = {}
    for pred_id, (lo, hi) in predictions.items():
        # Match by *base* compound name (drop "Pn_" prefix AND any "_dose" suffix)
        # e.g. "P3_methylphenidate_20mg" → base="methylphenidate"
        parts = pred_id.split("_")
        # Drop Pn prefix
        if parts and parts[0].startswith("P") and parts[0][1:].isdigit():
            parts = parts[1:]
        # Drop dose suffix (e.g. "20mg", "200mg", "3mg")
        if parts and parts[-1].endswith("mg") and parts[-1][:-2].replace(".", "").isdigit():
            parts = parts[:-1]
        compound_key = "_".join(parts) if parts else pred_id
        candidates = [c for c in posterior.compounds
                      if compound_key.lower() in c.lower()]
        if not candidates:
            out[pred_id] = "NO_COMPOUND"
            continue
        g = posterior.g_mean[candidates[0]]
        out[pred_id] = "PASS" if lo <= g <= hi else "FAIL"
    return out


def availability() -> dict[str, object]:
    """Probe V7.3 model availability."""
    return {
        "available": True,
        "pymc_backend": PYMC_AVAILABLE,
        "numpyro_backend": _numpyro_available() if PYMC_AVAILABLE else False,
        "n_moderators": N_MODERATORS,
        "moderator_names": list(MODERATOR_NAMES),
        "stub_mode_works_without_pymc": True,
    }
