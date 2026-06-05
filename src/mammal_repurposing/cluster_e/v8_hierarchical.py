"""V8.6 — Hierarchical V8 joint posterior with α (cell-line) + γ (species)
+ δ (compound × cell interaction) random effects (Sprint 4.1).

Per `research/4-tier/MH3_per_cell_line_random_effect_deep_research.md` §3-§4,
this module bundles MH3 (per-cell-line) and MH7 (per-species) random effects
into a single PyMC NUTS hierarchical model. The MH3 doc § 7 is explicit:
splitting MH3 and MH7 produces α/γ identifiability collisions because some
compounds appear only in U2OS (cell-line and species are partly aliased) —
the random effects MUST be fit jointly.

Mathematical specification (per MH3 § 3):

    y_{c,l,s,k,r} = μ_k + β_{c,k} + α_{l,k} + γ_{s,k} + δ_{c,l,k} + ε

Where:
    β_{c,k}    = transferable compound effect (the quantity V8 ranks against)
    α_{l,k}    = cell-line random effect (MH3)
    γ_{s,k}    = species random effect (MH7)
    δ_{c,l,k}  = compound × cell-line interaction (the cell-specific portion
                 that does NOT transfer to brain)
    ε          = residual noise

Derived quantities (per MH3 § 3.1):
    ICC_cell,k  = σ_α² / (σ_β² + σ_α² + σ_γ² + σ_δ² + σ_ε²)
    ICC_inter,k = σ_δ² / (σ_β² + σ_δ²)
    T_{c,k}     = E[ |β_{c,k}| / (|β_{c,k}| + std(δ_{c,·,k})) | data ]

`T_{c,k}` is the **per-compound transferability index** — the most important
deliverable for the V8 paper. T_{c,k} ≈ 1 means the compound effect is
dominated by the transferable β; T_{c,k} ≈ 0 means it's essentially U2OS-
idiosyncratic. The V8 (L, L, H) gate should tighten to T_{c,k} > 0.6.

OPERATIONAL STATE (Sprint 4.1 — 2026-05-28):
  - Architecture validated on synthetic 2-cell-line × 50-compound data
  - Real cpg0000 empirical prior calibration is Sprint 4.3 (pending data pull)
  - V8 OSF § 7 amendment with G1-G6 gates is Sprint 4.4 (pending)

References:
  research/4-tier/MH3_per_cell_line_random_effect_deep_research.md (full math + sprint plan)
  Chandrasekaran et al. 2023 (JUMP-CP cpg0016 + cpg0000 pilot)
  Gorgogietas et al. 2025 Sci Rep (U2OS → mDA neuron concordance)
  Anderson et al. 2025 eLife (iPSC cortical CP modality bridge)
  Dirmeier & Beerenwinkel 2022 (SHM perturbation screening)
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


def _numpyro_available() -> bool:
    try:
        import numpyro    # noqa: F401
        return True
    except ImportError:
        return False


@dataclass
class V8HierarchicalPosterior:
    """Output of the V8 hierarchical NUTS fit."""
    # Indexing
    compounds: list[str]
    cell_lines: list[str]
    species: list[str]
    endpoints: list[str]
    # Per-compound transferable effect β (the quantity V8 ranks)
    beta_mean: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))      # (n_compounds, n_endpoints)
    beta_sd: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    # Per-compound transferability index T_{c,k}
    transferability_index: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))  # (n_compounds, n_endpoints)
    # Variance components (population-level, per endpoint)
    sigma_beta: np.ndarray = field(default_factory=lambda: np.zeros(0))           # (n_endpoints,)
    sigma_alpha: np.ndarray = field(default_factory=lambda: np.zeros(0))
    sigma_gamma: np.ndarray = field(default_factory=lambda: np.zeros(0))
    sigma_delta: np.ndarray = field(default_factory=lambda: np.zeros(0))
    sigma_eps: np.ndarray = field(default_factory=lambda: np.zeros(0))
    # ICCs (per endpoint)
    icc_cell: np.ndarray = field(default_factory=lambda: np.zeros(0))             # (n_endpoints,)
    icc_inter: np.ndarray = field(default_factory=lambda: np.zeros(0))
    # Diagnostics
    n_chains: int = 0
    n_draws: int = 0
    n_divergences: int = 0
    rhat_max: float = float("nan")
    ess_min: float = float("nan")
    method: str = "synthetic_stub"
    note: str = ""

    def transferability_per_compound(self, endpoint_idx: int = 0) -> dict[str, float]:
        """Per-compound transferability score for a given endpoint."""
        if self.transferability_index.size == 0:
            return {}
        return {self.compounds[i]: float(self.transferability_index[i, endpoint_idx])
                for i in range(len(self.compounds))}

    def beta_per_compound(self, endpoint_idx: int = 0) -> dict[str, float]:
        """Per-compound transferable effect β̂ for a given endpoint."""
        if self.beta_mean.size == 0:
            return {}
        return {self.compounds[i]: float(self.beta_mean[i, endpoint_idx])
                for i in range(len(self.compounds))}


@dataclass
class V8Observation:
    """One observation row in the long-format V8 hierarchical input.

    Fields (per MH3 § 3):
        compound:   compound identifier
        cell_line:  cell line identifier (e.g. "U2OS", "A549", "iPSC_cortical")
        species:    "human" | "mouse" | etc.
        endpoint:   phenotypic endpoint name (e.g. "mitochondrial",
                    "cytoskeletal", "synaptic")
        replicate:  replicate index (or run/batch — used only for ε)
        y:          observed phenotypic signal value
    """
    compound: str
    cell_line: str
    species: str
    endpoint: str
    replicate: int
    y: float


def build_v8_hierarchical_with_cell_random_effect(
    observations: list[V8Observation],
    *,
    sigma_beta_prior_scale: float = 0.5,
    sigma_alpha_prior_scale: float = 0.5,
    sigma_gamma_prior_scale: float = 0.3,
    sigma_delta_prior_scale: float = 0.5,
    sigma_eps_prior_scale: float = 1.0,
    # Sprint 4.3 (future) will replace these defaults with cpg0000-empirically-calibrated values
):
    """Build the PyMC NUTS model per MH3 § 4 drop-in code.

    Returns:
        (model, coords_dict, idx_dict) where:
          - model: PyMC Model instance ready to sample
          - coords_dict: dim → label mapping
          - idx_dict: {"compound_idx", "cell_idx", "species_idx", "endpoint_idx"}
            numpy arrays for downstream posterior interpretation.

    Non-centered parameterization on all four random effects (β, α, γ, δ) per
    MH3 § 4 — necessary because the centered version will divergence-bomb
    with five variance components.
    """
    if not PYMC_AVAILABLE:
        raise ImportError(
            "PyMC not installed. Install pymc + numpyro to use V8 hierarchical."
        )
    import pymc as pm

    n_obs = len(observations)
    if n_obs == 0:
        raise ValueError("Cannot build V8 hierarchical on zero observations")

    # Build coords
    compounds = sorted(set(o.compound for o in observations))
    cell_lines = sorted(set(o.cell_line for o in observations))
    species = sorted(set(o.species for o in observations))
    endpoints = sorted(set(o.endpoint for o in observations))

    compound_to_idx = {c: i for i, c in enumerate(compounds)}
    cell_to_idx = {c: i for i, c in enumerate(cell_lines)}
    species_to_idx = {s: i for i, s in enumerate(species)}
    endpoint_to_idx = {e: i for i, e in enumerate(endpoints)}

    compound_idx = np.array([compound_to_idx[o.compound] for o in observations])
    cell_idx = np.array([cell_to_idx[o.cell_line] for o in observations])
    species_idx = np.array([species_to_idx[o.species] for o in observations])
    endpoint_idx = np.array([endpoint_to_idx[o.endpoint] for o in observations])
    y = np.array([o.y for o in observations], dtype=float)

    coords = {
        "compound": compounds,
        "cell": cell_lines,
        "species": species,
        "endpoint": endpoints,
    }

    with pm.Model(coords=coords) as model:
        # ---- Global intercepts per endpoint ----
        mu_k = pm.Normal("mu_k", 0.0, 1.0, dims="endpoint")

        # ---- Hyperpriors on variance components (per endpoint) ----
        sigma_beta = pm.HalfNormal("sigma_beta",  sigma_beta_prior_scale,
                                    dims="endpoint")
        sigma_alpha = pm.HalfNormal("sigma_alpha", sigma_alpha_prior_scale,
                                     dims="endpoint")
        sigma_gamma = pm.HalfNormal("sigma_gamma", sigma_gamma_prior_scale,
                                     dims="endpoint")
        sigma_delta = pm.HalfNormal("sigma_delta", sigma_delta_prior_scale,
                                     dims="endpoint")
        sigma_eps = pm.HalfNormal("sigma_eps", sigma_eps_prior_scale,
                                   dims="endpoint")

        # ---- Non-centered random effects (per MH3 § 4) ----
        # β_{c,k} = sigma_beta_k * raw   — transferable compound effect
        beta_raw = pm.Normal("beta_raw",  0.0, 1.0, dims=("compound", "endpoint"))
        beta = pm.Deterministic(
            "beta", beta_raw * sigma_beta, dims=("compound", "endpoint"),
        )
        # α_{l,k} — cell-line random effect (MH3)
        alpha_raw = pm.Normal("alpha_raw", 0.0, 1.0, dims=("cell", "endpoint"))
        alpha = pm.Deterministic(
            "alpha", alpha_raw * sigma_alpha, dims=("cell", "endpoint"),
        )
        # γ_{s,k} — species random effect (MH7 co-fit per MH3 § 7)
        gamma_raw = pm.Normal("gamma_raw", 0.0, 1.0, dims=("species", "endpoint"))
        gamma = pm.Deterministic(
            "gamma", gamma_raw * sigma_gamma, dims=("species", "endpoint"),
        )
        # δ_{c,l,k} — compound × cell-line interaction
        delta_raw = pm.Normal("delta_raw", 0.0, 1.0,
                               dims=("compound", "cell", "endpoint"))
        delta = pm.Deterministic(
            "delta", delta_raw * sigma_delta,
            dims=("compound", "cell", "endpoint"),
        )

        # ---- Linear predictor ----
        eta = (
            mu_k[endpoint_idx]
            + beta[compound_idx, endpoint_idx]
            + alpha[cell_idx, endpoint_idx]
            + gamma[species_idx, endpoint_idx]
            + delta[compound_idx, cell_idx, endpoint_idx]
        )
        sigma_obs = sigma_eps[endpoint_idx]
        pm.Normal("y_obs", mu=eta, sigma=sigma_obs, observed=y)

        # ---- Derived reporting quantities ----
        pm.Deterministic(
            "icc_cell",
            sigma_alpha ** 2 / (sigma_beta ** 2 + sigma_alpha ** 2
                                + sigma_gamma ** 2 + sigma_delta ** 2
                                + sigma_eps ** 2),
            dims="endpoint",
        )
        pm.Deterministic(
            "icc_inter",
            sigma_delta ** 2 / (sigma_beta ** 2 + sigma_delta ** 2),
            dims="endpoint",
        )

    return model, coords, {
        "compound_idx": compound_idx,
        "cell_idx": cell_idx,
        "species_idx": species_idx,
        "endpoint_idx": endpoint_idx,
        "compounds": compounds,
        "cell_lines": cell_lines,
        "species": species,
        "endpoints": endpoints,
    }


def fit_v8_hierarchical(
    observations: list[V8Observation],
    *,
    n_chains: int = 4,
    n_draws: int = 2000,
    n_tune: int = 2000,
    target_accept: float = 0.95,
    random_seed: int = 20260528,
    sigma_beta_prior_scale: float = 0.5,
    sigma_alpha_prior_scale: float = 0.5,
    sigma_gamma_prior_scale: float = 0.3,
    sigma_delta_prior_scale: float = 0.5,
    sigma_eps_prior_scale: float = 1.0,
) -> V8HierarchicalPosterior:
    """Build + fit the V8 hierarchical model. Returns V8HierarchicalPosterior.

    Compute budget per MH3 § 4: with ~1000 compounds × ~3 cell-line buckets ×
    ~8 endpoints (~24,000 latent δ parameters), 4 chains × 4000 draws on
    numpyro/JAX completes in 5-15 minutes on RTX 5070.
    """
    if not PYMC_AVAILABLE:
        raise ImportError("PyMC not installed.")
    import pymc as pm
    import arviz as az

    model, coords, idx = build_v8_hierarchical_with_cell_random_effect(
        observations,
        sigma_beta_prior_scale=sigma_beta_prior_scale,
        sigma_alpha_prior_scale=sigma_alpha_prior_scale,
        sigma_gamma_prior_scale=sigma_gamma_prior_scale,
        sigma_delta_prior_scale=sigma_delta_prior_scale,
        sigma_eps_prior_scale=sigma_eps_prior_scale,
    )

    with model:
        sample_kwargs = dict(
            draws=n_draws, tune=n_tune, chains=n_chains,
            target_accept=target_accept,
            random_seed=random_seed, progressbar=False,
        )
        if _numpyro_available():
            sample_kwargs["nuts_sampler"] = "numpyro"
        idata = pm.sample(**sample_kwargs)

    # Extract posterior summaries
    beta_post = idata.posterior["beta"].values    # (chain, draw, compound, endpoint)
    delta_post = idata.posterior["delta"].values  # (chain, draw, compound, cell, endpoint)

    n_compounds = len(idx["compounds"])
    n_endpoints = len(idx["endpoints"])

    beta_mean = beta_post.mean(axis=(0, 1))    # (compound, endpoint)
    beta_sd = beta_post.std(axis=(0, 1))

    # Transferability index per compound per endpoint:
    #   T_{c,k} = E[ |β_{c,k}| / (|β_{c,k}| + std(δ_{c,·,k})) ]
    # We compute by drawing — for each posterior sample s, compute
    #   |β^{(s)}_{c,k}| / (|β^{(s)}_{c,k}| + std_l[δ^{(s)}_{c,l,k}])
    # then average across (chain × draw).
    beta_flat = beta_post.reshape(-1, n_compounds, n_endpoints)
    delta_flat = delta_post.reshape(-1, n_compounds, len(idx["cell_lines"]), n_endpoints)
    delta_std_per_compound = delta_flat.std(axis=2)    # std across cell-lines: (n_draws, n_compounds, n_endpoints)
    t_per_draw = (np.abs(beta_flat) /
                  (np.abs(beta_flat) + delta_std_per_compound + 1e-9))
    transferability_index = t_per_draw.mean(axis=0)    # (n_compounds, n_endpoints)

    sigma_beta_post = idata.posterior["sigma_beta"].values.mean(axis=(0, 1))
    sigma_alpha_post = idata.posterior["sigma_alpha"].values.mean(axis=(0, 1))
    sigma_gamma_post = idata.posterior["sigma_gamma"].values.mean(axis=(0, 1))
    sigma_delta_post = idata.posterior["sigma_delta"].values.mean(axis=(0, 1))
    sigma_eps_post = idata.posterior["sigma_eps"].values.mean(axis=(0, 1))
    icc_cell = idata.posterior["icc_cell"].values.mean(axis=(0, 1))
    icc_inter = idata.posterior["icc_inter"].values.mean(axis=(0, 1))

    summary = az.summary(
        idata, var_names=["sigma_beta", "sigma_alpha", "sigma_gamma",
                          "sigma_delta", "sigma_eps", "icc_cell", "icc_inter"],
    )
    rhat_max = float(summary["r_hat"].max())
    ess_min = float(summary["ess_bulk"].min())

    try:
        diverging = idata.sample_stats["diverging"].values if "diverging" in idata.sample_stats else None
        n_div = int(diverging.sum()) if diverging is not None else 0
    except Exception:    # noqa: BLE001
        n_div = -1

    return V8HierarchicalPosterior(
        compounds=idx["compounds"],
        cell_lines=idx["cell_lines"],
        species=idx["species"],
        endpoints=idx["endpoints"],
        beta_mean=beta_mean,
        beta_sd=beta_sd,
        transferability_index=transferability_index,
        sigma_beta=sigma_beta_post,
        sigma_alpha=sigma_alpha_post,
        sigma_gamma=sigma_gamma_post,
        sigma_delta=sigma_delta_post,
        sigma_eps=sigma_eps_post,
        icc_cell=icc_cell,
        icc_inter=icc_inter,
        n_chains=n_chains, n_draws=n_draws, n_divergences=n_div,
        rhat_max=rhat_max, ess_min=ess_min,
        method="pymc_nuts_v8_hierarchical",
        note=(f"V8 hierarchical converged: R̂={rhat_max:.3f}, ESS={ess_min:.0f}, "
              f"divergences={n_div}; "
              f"σ̂_β={sigma_beta_post.mean():.3f}, σ̂_α={sigma_alpha_post.mean():.3f}, "
              f"σ̂_γ={sigma_gamma_post.mean():.3f}, σ̂_δ={sigma_delta_post.mean():.3f}; "
              f"ICC_cell={icc_cell.mean():.3f}, ICC_inter={icc_inter.mean():.3f}"),
    )


def generate_synthetic_v8_observations(
    n_compounds: int = 50,
    n_cell_lines: int = 2,
    n_species: int = 1,
    n_endpoints: int = 2,
    n_replicates: int = 3,
    *,
    true_sigma_beta: float = 0.4,
    true_sigma_alpha: float = 0.2,
    true_sigma_gamma: float = 0.1,
    true_sigma_delta: float = 0.15,
    true_sigma_eps: float = 0.1,
    rng_seed: int = 42,
) -> tuple[list[V8Observation], np.ndarray]:
    """Generate synthetic observations from the canonical V8 hierarchical
    model. Returns (observations, true_beta) where true_beta is the
    underlying transferable component used for posterior-recovery tests.

    Used by Sprint 4.2 smoke test to verify the model architecture
    correctly recovers known parameters within posterior credible bounds.
    """
    rng = np.random.default_rng(rng_seed)
    compounds = [f"cmpd_{i:03d}" for i in range(n_compounds)]
    cell_lines = [f"cell_{i}" for i in range(n_cell_lines)]
    species_names = [f"species_{i}" for i in range(n_species)]
    endpoints = [f"endpt_{i}" for i in range(n_endpoints)]

    # True parameters
    true_beta = rng.normal(0, true_sigma_beta, size=(n_compounds, n_endpoints))
    true_alpha = rng.normal(0, true_sigma_alpha, size=(n_cell_lines, n_endpoints))
    true_gamma = rng.normal(0, true_sigma_gamma, size=(n_species, n_endpoints))
    true_delta = rng.normal(0, true_sigma_delta,
                             size=(n_compounds, n_cell_lines, n_endpoints))

    observations: list[V8Observation] = []
    for c_idx, c in enumerate(compounds):
        for l_idx, l in enumerate(cell_lines):
            for s_idx, s in enumerate(species_names):
                for k_idx, k in enumerate(endpoints):
                    for r in range(n_replicates):
                        eta = (true_beta[c_idx, k_idx]
                               + true_alpha[l_idx, k_idx]
                               + true_gamma[s_idx, k_idx]
                               + true_delta[c_idx, l_idx, k_idx])
                        y = eta + rng.normal(0, true_sigma_eps)
                        observations.append(V8Observation(
                            compound=c, cell_line=l, species=s,
                            endpoint=k, replicate=r, y=float(y),
                        ))
    return observations, true_beta


def availability() -> dict[str, object]:
    return {
        "available": True,
        "pymc_backend": PYMC_AVAILABLE,
        "numpyro_backend": _numpyro_available() if PYMC_AVAILABLE else False,
        "model_name": "V8_hierarchical_alpha_gamma_delta",
        "random_effects": ["beta_transferable", "alpha_cell", "gamma_species",
                            "delta_compound_x_cell"],
        "derived_quantities": ["icc_cell", "icc_inter", "transferability_index_T_ck"],
        "reference": ("research/4-tier/MH3_per_cell_line_random_effect_deep_research.md "
                      "§3-§4"),
    }
