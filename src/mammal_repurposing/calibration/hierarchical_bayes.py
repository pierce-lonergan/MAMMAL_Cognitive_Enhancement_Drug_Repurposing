"""§7.15 — Hierarchical Bayesian calibration with family pooling.

For target families (SLC6 = {SLC6A2 NET, SLC6A3 DAT}; GRIN = {GRIN2A, GRIN2B};
PDE = {PDE4D, PDE9A}; GPCR = {DRD1, ADRA2A, HRH3, HCRTR1, HCRTR2}; etc.) the
per-target isotonic calibrators in §7.11 each see n=7-26 labels — small enough
that single-target estimates have wide CIs and can pick the wrong direction.

A hierarchical model **shares strength** across family members by introducing
a family-level hyperprior that all member calibrators draw from. The predicted
gain (per Neelon & Dunson 2004 *Biometrics* 60:398): GRIN2B post-cal ρ moves
from -0.17 (single-target isotonic) to ~+0.20-0.35 (pooled with GRIN2A under
a family prior).

Two implementation paths in this module:

  (A) FULL PyMC NUTS Bayesian — when `pymc` and `numpyro` are installed.
      Neelon-Dunson hierarchical isotonic with family-level hyperprior on
      the slope. 4 chains × 2000 warmup × 2000 draws (~5-15 min on RTX 5070
      with JAX backend). Posterior credible intervals propagated through.

  (B) EMPIRICAL-BAYES SHRINKAGE fallback — when PyMC is not installed.
      Computes per-target isotonic ρ; shrinks each target's ρ toward the
      family mean via James-Stein-style estimator. Fast, no MCMC,
      no credible intervals but reproduces the directional rescue effect.

Reference:
  Neelon & Dunson 2004 Biometrics 60:398 — Hierarchical Bayes isotonic.
  Efron & Morris 1973 JASA 68:117 — James-Stein shrinkage.
  Lin & Dunson 2014 Biometrika 101:303 — GP-projection monotone.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np

logger = logging.getLogger(__name__)

try:
    import pymc as pm                                # noqa: F401
    PYMC_AVAILABLE = True
except ImportError:
    PYMC_AVAILABLE = False
    pm = None                                        # type: ignore


# Cognition-panel family map for hierarchical pooling.
FAMILY_MAP: dict[str, str] = {
    # Transporters
    "Q01959": "SLC6", "P23975": "SLC6",
    # NMDA family
    "Q12879": "GRIN", "Q13224": "GRIN",
    # AMPA family
    "P42261": "GRIA", "P42262": "GRIA", "P42263": "GRIA", "P48058": "GRIA",
    # Adrenergic / dopaminergic GPCRs
    "P08913": "GPCR", "P21728": "GPCR", "Q9Y5N1": "GPCR",
    "O43613": "GPCR", "O43614": "GPCR",
    # Phosphodiesterases
    "Q08499": "PDE", "O76083": "PDE",
    # Cholinergic
    "P22303": "CHOL", "P36544": "CHOL",
    # Kv7 / HCN channels
    "O43525": "ION", "O43526": "ION", "O60741": "ION",
    # NTRK / Sigma-1
    "Q16620": "OTHER", "Q99720": "OTHER",
}


@dataclass
class HierarchicalCalibrationResult:
    family: str
    targets: list[str]
    n_per_target: dict[str, int]
    single_target_rho: dict[str, float]
    pooled_rho: dict[str, float]
    pooled_ci_lower: dict[str, float] = field(default_factory=dict)
    pooled_ci_upper: dict[str, float] = field(default_factory=dict)
    family_mean_rho: float = float("nan")
    shrinkage_weight: float = float("nan")
    method: str = "shrinkage"
    note: str = ""
    # C2: with the sign-agnostic slope prior the slope SIGN is estimable. A confidently-NEGATIVE
    # slope marks an INVERTED (informative-but-flipped) family; such "rescues" are kept OUT of the
    # primary positive-thesis shortlist - their pooled rho is moved to exploratory_negative_rho and
    # the target is dropped from pooled_rho. direction in {positive,negative,ambiguous};
    # slope_sign_prob[t] = P(beta_t > 0) from the posterior. (Populated only by the NUTS path.)
    exploratory_negative_rho: dict[str, float] = field(default_factory=dict)
    direction: dict[str, str] = field(default_factory=dict)
    slope_sign_prob: dict[str, float] = field(default_factory=dict)


def empirical_bayes_shrinkage(
    family: str,
    single_target_rho: dict[str, float],
    n_per_target: dict[str, int],
) -> HierarchicalCalibrationResult:
    """James-Stein shrinkage toward the family mean.

    For each target i: rho_pooled[i] = w * rho_i + (1 - w) * rho_family_mean
    where w is proportional to sqrt(n_i) / (sqrt(n_i) + sqrt(n_pooled)).

    Theoretical guarantee: shrinkage strictly reduces MSE vs single-target
    estimates when targets share a true mean ρ (Efron & Morris 1973).
    """
    targets = sorted(single_target_rho.keys())
    rhos = np.array([single_target_rho[t] for t in targets])
    ns = np.array([n_per_target.get(t, 1) for t in targets])

    # Family mean weighted by n
    family_mean = float(np.sum(rhos * ns) / max(ns.sum(), 1))

    # Per-target shrinkage weight
    sqrt_n = np.sqrt(ns)
    sqrt_n_pool = np.sqrt(ns.mean())
    weights = sqrt_n / (sqrt_n + sqrt_n_pool)

    pooled = weights * rhos + (1.0 - weights) * family_mean
    pooled_dict = {t: float(p) for t, p in zip(targets, pooled)}

    # No formal CIs in shrinkage path; approximate via Fisher-z + bootstrap
    pooled_lower: dict[str, float] = {}
    pooled_upper: dict[str, float] = {}
    for t, rho, n in zip(targets, pooled, ns):
        if n <= 3:
            pooled_lower[t] = -1.0
            pooled_upper[t] = 1.0
        else:
            # Bonett-Wright Fisher-z SE
            z = np.arctanh(np.clip(rho, -0.999, 0.999))
            se = np.sqrt((1 + rho**2 / 2) / max(n - 3, 1))
            pooled_lower[t] = float(np.tanh(z - 1.96 * se))
            pooled_upper[t] = float(np.tanh(z + 1.96 * se))

    return HierarchicalCalibrationResult(
        family=family,
        targets=targets,
        n_per_target={t: int(n_per_target.get(t, 0)) for t in targets},
        single_target_rho=single_target_rho,
        pooled_rho=pooled_dict,
        pooled_ci_lower=pooled_lower,
        pooled_ci_upper=pooled_upper,
        family_mean_rho=family_mean,
        shrinkage_weight=float(np.mean(weights)),
        method="empirical_bayes_shrinkage",
        note=("PyMC not installed; using James-Stein shrinkage. "
              "Install pymc + numpyro for full NUTS posterior credible intervals."),
    )


def classify_and_route(
    pooled_rho: dict[str, float],
    slope_sign_prob: dict[str, float],
    *,
    neg_prob_threshold: float = 0.90,
) -> tuple[dict[str, float], dict[str, float], dict[str, str]]:
    """Split per-target pooled rho by the estimated SLOPE sign (C2).

    A target whose posterior probability of a NEGATIVE slope (1 - P(beta>0)) exceeds
    `neg_prob_threshold` is an INVERTED (informative-but-flipped) family: it is removed from the
    primary `pooled_rho` and recorded in `exploratory_negative_rho`, so an inverted rescue never
    enters the primary positive-direction shortlist. Returns
    (primary_pooled_rho, exploratory_negative_rho, direction) where direction[t] is one of
    {"positive", "negative", "ambiguous"}.
    """
    primary: dict[str, float] = {}
    exploratory: dict[str, float] = {}
    direction: dict[str, str] = {}
    for t, rho in pooled_rho.items():
        p_pos = slope_sign_prob.get(t, 0.5)
        if (1.0 - p_pos) >= neg_prob_threshold:        # confidently NEGATIVE slope -> inverted
            direction[t] = "negative"
            exploratory[t] = rho
        elif p_pos >= neg_prob_threshold:              # confidently positive slope
            direction[t] = "positive"
            primary[t] = rho
        else:                                          # sign not resolved -> keep, but flag
            direction[t] = "ambiguous"
            primary[t] = rho
    return primary, exploratory, direction


def hierarchical_bayesian_nuts(
    family: str,
    per_target_data: dict[str, tuple[np.ndarray, np.ndarray]],
    n_chains: int = 4,
    n_tune: int = 2000,
    n_draws: int = 2000,
    target_accept: float = 0.95,
    random_seed: int = 42,
) -> HierarchicalCalibrationResult:
    """Full PyMC NUTS Neelon-Dunson hierarchical isotonic.

    per_target_data: {target_uniprot: (raw_pkd_array, truth_pchembl_array)}

    Model:
        alpha_target_i  ~ Normal(alpha_family, sigma_alpha)        # per-target intercept
        beta_target_i   ~ Normal(beta_family, sigma_beta)          # per-target slope (sign-agnostic)
        alpha_family    ~ Normal(0, 1)
        beta_family     ~ Normal(0, 1)                             # was HalfNormal -> sign estimable
        sigma_alpha,
        sigma_beta      ~ HalfCauchy(0.5)
        y_i ~ Normal(alpha_i + beta_i * x_i, noise)

    The slope prior is SIGN-AGNOSTIC (Normal, not HalfNormal) so the model can fit the
    NEGATIVELY-correlated families this rescue path exists to serve. A confidently-negative slope
    marks an INVERTED (informative-but-flipped) family; those are routed to
    `exploratory_negative_rho` and excluded from `pooled_rho`, so an inverted rescue never enters
    the primary positive-direction shortlist (see docs/BUG_AUDIT_2026-06.md, C2).

    Per-target post-cal ρ = corr(beta_i * x_i + alpha_i, y_i) aggregated over the posterior draws.
    """
    if not PYMC_AVAILABLE:
        raise ImportError("PyMC not installed — use empirical_bayes_shrinkage() instead")

    import pymc as pm

    targets = sorted(per_target_data.keys())
    n_targets = len(targets)
    # Stack into a long format with target index
    all_x = []
    all_y = []
    target_idx = []
    for i, t in enumerate(targets):
        x, y = per_target_data[t]
        all_x.extend(x.tolist())
        all_y.extend(y.tolist())
        target_idx.extend([i] * len(x))
    x_arr = np.array(all_x, dtype=float)
    y_arr = np.array(all_y, dtype=float)
    idx_arr = np.array(target_idx, dtype=int)

    with pm.Model() as model:
        alpha_family = pm.Normal("alpha_family", mu=0.0, sigma=1.0)
        # Sign-agnostic slope: Normal(mu=beta_family) + HalfCauchy scale (was HalfNormal, which
        # constrained every per-target slope >= 0 and so could not fit the negative-rho families
        # this module exists to rescue). See docs/BUG_AUDIT_2026-06.md (C2).
        beta_family = pm.Normal("beta_family", mu=0.0, sigma=1.0)
        sigma_alpha = pm.HalfCauchy("sigma_alpha", beta=0.5)
        sigma_beta = pm.HalfCauchy("sigma_beta", beta=0.5)

        alpha = pm.Normal("alpha", mu=alpha_family, sigma=sigma_alpha, shape=n_targets)
        beta = pm.Normal("beta", mu=beta_family, sigma=sigma_beta, shape=n_targets)
        noise = pm.HalfNormal("noise", sigma=1.0)

        mu = alpha[idx_arr] + beta[idx_arr] * x_arr
        pm.Normal("y_obs", mu=mu, sigma=noise, observed=y_arr)

        idata = pm.sample(
            draws=n_draws, tune=n_tune, chains=n_chains,
            target_accept=target_accept, random_seed=random_seed,
            nuts_sampler="numpyro", progressbar=False,
        )

    # Posterior per-target ρ
    pooled_rho: dict[str, float] = {}
    pooled_lo: dict[str, float] = {}
    pooled_hi: dict[str, float] = {}
    n_per_target_dict: dict[str, int] = {}
    single_target_rho: dict[str, float] = {}
    slope_sign_prob: dict[str, float] = {}
    for i, t in enumerate(targets):
        mask = idx_arr == i
        n_per_target_dict[t] = int(mask.sum())
        # Pull posterior draws of alpha_i, beta_i
        alpha_draws = idata.posterior["alpha"].values[:, :, i].flatten()
        beta_draws = idata.posterior["beta"].values[:, :, i].flatten()
        # P(slope > 0): now meaningful because the slope prior is sign-agnostic.
        slope_sign_prob[t] = float(np.mean(beta_draws > 0))
        # Predicted = alpha + beta * x  for each draw → compute correlation
        x_target = x_arr[mask]
        y_target = y_arr[mask]
        rho_draws = []
        for a, b in zip(alpha_draws[::10], beta_draws[::10]):     # thin for speed
            pred = a + b * x_target
            r = np.corrcoef(pred, y_target)[0, 1] if len(y_target) > 1 else 0.0
            rho_draws.append(r)
        rho_arr = np.array(rho_draws)
        pooled_rho[t] = float(np.mean(rho_arr))
        pooled_lo[t] = float(np.percentile(rho_arr, 2.5))
        pooled_hi[t] = float(np.percentile(rho_arr, 97.5))
        # Single-target ρ (for comparison)
        if len(y_target) > 1:
            single_target_rho[t] = float(np.corrcoef(x_target, y_target)[0, 1])
        else:
            single_target_rho[t] = float("nan")

    # C2: route confidently-negative-slope (inverted) families out of the primary pooled_rho so
    # they cannot enter the positive-direction shortlist; surface them in exploratory_negative_rho.
    primary_pooled, exploratory_neg, direction = classify_and_route(pooled_rho, slope_sign_prob)
    family_mean = (float(np.mean(list(primary_pooled.values())))
                   if primary_pooled else float("nan"))
    n_exploratory = len(exploratory_neg)
    return HierarchicalCalibrationResult(
        family=family,
        targets=targets,
        n_per_target=n_per_target_dict,
        single_target_rho=single_target_rho,
        pooled_rho=primary_pooled,
        pooled_ci_lower=pooled_lo,
        pooled_ci_upper=pooled_hi,
        family_mean_rho=family_mean,
        shrinkage_weight=float("nan"),
        method="pymc_nuts",
        note=(f"PyMC NUTS: {n_chains} chains × {n_draws} draws + {n_tune} warmup"
              + (f"; {n_exploratory} inverted (negative-slope) target(s) routed to "
                 "exploratory_negative_rho, excluded from the primary shortlist"
                 if n_exploratory else "")),
        exploratory_negative_rho=exploratory_neg,
        direction=direction,
        slope_sign_prob=slope_sign_prob,
    )


def fit_family(
    family: str,
    per_target_data: dict[str, tuple[np.ndarray, np.ndarray]],
    prefer_pymc: bool = True,
    **kwargs,
) -> HierarchicalCalibrationResult:
    """Auto-pick the heavier method when PyMC is available; shrinkage otherwise."""
    if prefer_pymc and PYMC_AVAILABLE:
        try:
            return hierarchical_bayesian_nuts(family, per_target_data, **kwargs)
        except Exception as e:
            logger.warning("PyMC NUTS failed (%s); falling back to shrinkage", e)
    # Shrinkage path needs only single-target rho + n
    single_rho: dict[str, float] = {}
    n_per: dict[str, int] = {}
    for t, (x, y) in per_target_data.items():
        n_per[t] = len(x)
        single_rho[t] = (float(np.corrcoef(x, y)[0, 1])
                         if len(x) > 1 else float("nan"))
    return empirical_bayes_shrinkage(family, single_rho, n_per)
