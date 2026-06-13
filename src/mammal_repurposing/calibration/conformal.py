"""§7.12 — Inductive conformal prediction per target.

Wraps the §7.11 isotonic calibrators with an inductive conformal layer that
emits per-(compound, target) prediction intervals with explicit marginal
coverage guarantees.

The recipe (Vovk-Gammerman split-conformal):
    1. Hold out a calibration fold from the ChEMBL ground-truth at target t.
    2. Fit isotonic on the training fold.
    3. Compute residuals on calibration fold: |y_cal − isotonic(x_cal)|.
    4. q_alpha = ⌈(n_cal + 1)(1 − alpha)⌉-th order statistic of residuals.
    5. At inference: predicted pKd ± q_alpha is a marginally-valid 1−alpha
       prediction interval.

Default alpha = 0.20 → 80% marginal coverage (looser than usual to reflect
the n=7-26 per-target regime — at higher confidence levels the interval
becomes uninformative).

When the per-target calibration set is too small (n < 10) we fall back to a
LOCO-residual interval: use the leave-one-compound-out residuals already
computed by the §7.11 diagnostics module.

Reference:
  Vovk, Gammerman & Shafer 2005. _Algorithmic Learning in a Random World_.
  Papadopoulos et al. 2007 IEEE TPAMI 29:9 — inductive conformal predictors.
  Mervin et al. 2020 J Chem Inf Model 60:4546 — Venn-ABERS (alt UQ for DTI).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
from sklearn.isotonic import IsotonicRegression

logger = logging.getLogger(__name__)


@dataclass
class ConformalCalibratorResult:
    target_uniprot: str
    n_train: int
    n_cal: int
    alpha: float
    q_alpha: float                      # half-width of the prediction interval
    empirical_coverage: float | None    # measured on the cal fold (sanity)
    isotonic: IsotonicRegression
    raw_min: float
    raw_max: float


def fit_inductive_conformal(
    raw_pkd: np.ndarray,
    truth_pchembl: np.ndarray,
    target_uniprot: str,
    alpha: float = 0.20,
    cal_frac: float = 0.30,
    seed: int = 42,
    direction: str | bool = "auto",
) -> ConformalCalibratorResult:
    """Split-conformal isotonic.

    Args:
        raw_pkd:        MAMMAL predicted_pkd (training + calibration).
        truth_pchembl:  ChEMBL pchembl ground truth (aligned with raw_pkd).
        target_uniprot: target identifier (book-keeping only).
        alpha:          target miscoverage rate; 1−alpha is the nominal coverage.
        cal_frac:       fraction of the data to hold out as the calibration
                        fold; default 0.30.
        direction:      isotonic direction; 'auto' is recommended for n ≥ 15.
    """
    n = len(raw_pkd)
    if n < 10:
        raise ValueError(f"Need ≥10 points for split-conformal; got n={n} at {target_uniprot}")

    rng = np.random.default_rng(seed)
    idx = np.arange(n)
    rng.shuffle(idx)
    n_cal = max(5, int(round(cal_frac * n)))
    cal_idx, train_idx = idx[:n_cal], idx[n_cal:]

    iso = IsotonicRegression(
        increasing="auto" if direction == "auto" else bool(direction),
        out_of_bounds="clip",
        y_min=2.0,
        y_max=11.0,
    )
    iso.fit(raw_pkd[train_idx], truth_pchembl[train_idx])

    # Residuals on cal fold
    cal_pred = iso.predict(raw_pkd[cal_idx])
    residuals = np.abs(truth_pchembl[cal_idx] - cal_pred)
    # q_alpha = ceil((n_cal + 1) * (1 - alpha)) order statistic. When that rank exceeds n_cal the
    # finite-sample (1-alpha) quantile does not exist, so the valid conformal interval is +inf
    # (abstain). Clamping to the max residual instead (the old behaviour) silently UNDER-covers at
    # tight alpha / small n_cal. The finite-rank path is unchanged.
    n_cal_actual = len(residuals)
    rank = int(np.ceil((n_cal_actual + 1) * (1 - alpha)))
    sorted_res = np.sort(residuals)
    q_alpha = float("inf") if rank > n_cal_actual else float(sorted_res[max(1, rank) - 1])

    # Empirical coverage on the SAME cal fold (sanity, not validity)
    in_interval = residuals <= q_alpha
    emp_cov = float(in_interval.mean()) if len(in_interval) else None

    return ConformalCalibratorResult(
        target_uniprot=target_uniprot,
        n_train=int(len(train_idx)),
        n_cal=int(n_cal_actual),
        alpha=alpha,
        q_alpha=q_alpha,
        empirical_coverage=emp_cov,
        isotonic=iso,
        raw_min=float(raw_pkd.min()),
        raw_max=float(raw_pkd.max()),
    )


def predict_with_interval(
    result: ConformalCalibratorResult,
    raw_pkd: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Returns (point_estimate, lower_bound, upper_bound)."""
    point = result.isotonic.predict(raw_pkd)
    lower = np.clip(point - result.q_alpha, 2.0, 11.0)
    upper = np.clip(point + result.q_alpha, 2.0, 11.0)
    return point, lower, upper


def loco_residuals_fallback(
    raw_pkd: np.ndarray,
    truth_pchembl: np.ndarray,
    direction: str | bool = "auto",
) -> np.ndarray:
    """Leave-one-out residuals — fallback when n < 10 is too small to split.

    Returns the array of |y_i − iso_{-i}(x_i)| values.
    """
    n = len(raw_pkd)
    residuals = np.empty(n, dtype=float)
    for i in range(n):
        mask = np.arange(n) != i
        iso = IsotonicRegression(
            increasing="auto" if direction == "auto" else bool(direction),
            out_of_bounds="clip",
            y_min=2.0,
            y_max=11.0,
        )
        iso.fit(raw_pkd[mask], truth_pchembl[mask])
        residuals[i] = abs(truth_pchembl[i] - iso.predict([raw_pkd[i]])[0])
    return residuals


def q_alpha_from_loco(residuals: np.ndarray, alpha: float = 0.20) -> float:
    """Compute q_alpha from a LOCO residual array (fallback path)."""
    if residuals.size == 0:
        return float("nan")
    n = residuals.size
    rank = int(np.ceil((n + 1) * (1 - alpha)))
    if rank > n:
        return float("inf")     # (1-alpha) order statistic does not exist -> abstain
    return float(np.sort(residuals)[max(1, rank) - 1])
