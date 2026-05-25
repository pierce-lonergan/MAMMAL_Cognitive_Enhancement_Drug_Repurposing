"""Bonett-Wright Fisher-z CI + permutation tests for Spearman ρ.

The research doc's blocking concern: ρ = -0.71 at n=26 is real (95% CI
[-0.86, -0.45]); but ρ = -0.35 at n=8 (GRIN2A) has CI [-0.83, +0.43] —
indistinguishable from zero. Don't invest engineering time on noise.

Reference: Bonett & Wright 2000, *Psychometrika* 65(1):23-28.
Spearman ρ permutation: standard non-parametric two-sided test.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


@dataclass
class PowerResult:
    target_uniprot: str
    target_gene: str
    n: int
    rho: float
    fisher_z_ci_low: float
    fisher_z_ci_high: float
    perm_p_value: float        # two-sided
    perm_ci_low: float         # 2.5th percentile of null
    perm_ci_high: float        # 97.5th percentile of null
    distinguishable_from_zero: bool
    verdict: str               # REAL | MARGINAL | NOISE


def fisher_z_ci_spearman(rho: float, n: int, alpha: float = 0.05) -> tuple[float, float]:
    """Bonett-Wright Fisher-z CI for Spearman ρ.

    Uses the variance correction from Bonett & Wright 2000 (Eq 6):
    Var(z) = (1 + ρ²/2) / (n - 3) — accounts for the slight underestimation
    in the classic Fisher-z formula at finite n.
    """
    if abs(rho) >= 1.0 or n < 4:
        return float("nan"), float("nan")
    # Fisher z transform
    z = 0.5 * math.log((1 + rho) / (1 - rho))
    # Bonett-Wright variance
    se = math.sqrt((1 + rho * rho / 2) / (n - 3))
    z_crit = 1.959963984540054  # 2-sided 95%
    z_lo = z - z_crit * se
    z_hi = z + z_crit * se
    # Back-transform
    rho_lo = (math.exp(2 * z_lo) - 1) / (math.exp(2 * z_lo) + 1)
    rho_hi = (math.exp(2 * z_hi) - 1) / (math.exp(2 * z_hi) + 1)
    return rho_lo, rho_hi


def permutation_test_rho(
    x: np.ndarray,
    y: np.ndarray,
    n_iter: int = 10000,
    seed: int = 42,
) -> tuple[float, float, float]:
    """Two-sided permutation test for Spearman ρ.

    Returns: (observed_rho, p_value, ci_low, ci_high) — strict null built by
    shuffling x; CIs are the 2.5/97.5 percentiles of the null distribution.
    """
    rng = np.random.default_rng(seed)
    if len(x) < 4:
        return float("nan"), float("nan"), float("nan")
    obs = spearmanr(x, y)[0]
    if math.isnan(obs):
        return float("nan"), float("nan"), float("nan")
    null = np.empty(n_iter, dtype=np.float64)
    for i in range(n_iter):
        null[i] = spearmanr(rng.permutation(x), y)[0]
    # Two-sided p-value
    p = (np.sum(np.abs(null) >= abs(obs)) + 1) / (n_iter + 1)
    ci_lo = float(np.percentile(null, 2.5))
    ci_hi = float(np.percentile(null, 97.5))
    return float(p), ci_lo, ci_hi


def evaluate(
    target_uniprot: str,
    target_gene: str,
    pred: np.ndarray,
    truth: np.ndarray,
    *,
    n_perm: int = 10000,
    seed: int = 42,
) -> PowerResult:
    """Full power analysis for one target."""
    if len(pred) < 4 or len(truth) < 4 or len(pred) != len(truth):
        return PowerResult(
            target_uniprot=target_uniprot, target_gene=target_gene,
            n=len(pred), rho=float("nan"),
            fisher_z_ci_low=float("nan"), fisher_z_ci_high=float("nan"),
            perm_p_value=float("nan"), perm_ci_low=float("nan"), perm_ci_high=float("nan"),
            distinguishable_from_zero=False, verdict="INSUFFICIENT_N",
        )

    rho_obs, _ = spearmanr(pred, truth)
    fz_lo, fz_hi = fisher_z_ci_spearman(rho_obs, len(pred))
    p, perm_lo, perm_hi = permutation_test_rho(pred, truth, n_iter=n_perm, seed=seed)

    # Distinguishable from zero — CI excludes 0
    distinguishable = not (fz_lo <= 0 <= fz_hi)
    if distinguishable and p < 0.05:
        verdict = "REAL"
    elif (fz_lo <= 0 <= fz_hi) and abs(rho_obs) > 0.2:
        verdict = "MARGINAL"     # large rho but wide CI — n too small
    else:
        verdict = "NOISE"

    return PowerResult(
        target_uniprot=target_uniprot, target_gene=target_gene,
        n=len(pred), rho=float(rho_obs),
        fisher_z_ci_low=fz_lo, fisher_z_ci_high=fz_hi,
        perm_p_value=float(p),
        perm_ci_low=perm_lo, perm_ci_high=perm_hi,
        distinguishable_from_zero=distinguishable, verdict=verdict,
    )


def power_panel(joined_records: dict[str, dict]) -> pd.DataFrame:
    """Compute power per target.

    joined_records: {target_uniprot: {"gene": str, "pred": np.array, "truth": np.array}}
    """
    rows = []
    for uniprot, payload in joined_records.items():
        r = evaluate(
            uniprot, payload.get("gene", "?"),
            payload["pred"], payload["truth"],
        )
        rows.append({
            "target_uniprot": r.target_uniprot,
            "gene": r.target_gene,
            "n": r.n,
            "rho": r.rho,
            "fisher_ci": f"[{r.fisher_z_ci_low:+.2f}, {r.fisher_z_ci_high:+.2f}]",
            "perm_p": r.perm_p_value,
            "perm_null_ci": f"[{r.perm_ci_low:+.2f}, {r.perm_ci_high:+.2f}]",
            "distinguishable_from_zero": r.distinguishable_from_zero,
            "verdict": r.verdict,
        })
    return pd.DataFrame(rows)
