"""LOCO + bootstrap diagnostics shared across calibration methods.

LOCO (leave-one-compound-out): for each calibration point i, refit the
calibrator on the other n-1 points and predict on i. Gives an honest
out-of-fold prediction without losing any data.

Bootstrap CI on Spearman ρ: 1,000 resamples with replacement; report
2.5/97.5 percentiles of the post-cal ρ.
"""

from __future__ import annotations

import numpy as np
from scipy.stats import spearmanr


def loco_predictions(
    raw: np.ndarray,
    truth: np.ndarray,
    fit_predict_fn,
) -> np.ndarray:
    """Leave-one-out predictions.

    Args:
        raw: 1D array of raw predictor (e.g. MAMMAL pKd)
        truth: 1D array of ground truth
        fit_predict_fn: callable (X_train, y_train, X_test) -> y_pred[len(X_test)]

    Returns:
        1D array of LOCO out-of-fold predictions (same shape as raw).
    """
    raw = np.asarray(raw, dtype=float)
    truth = np.asarray(truth, dtype=float)
    n = len(raw)
    out = np.zeros(n, dtype=float)
    for i in range(n):
        mask = np.arange(n) != i
        try:
            pred = fit_predict_fn(raw[mask], truth[mask], np.array([raw[i]]))
            out[i] = float(pred[0])
        except Exception:
            out[i] = float("nan")
    return out


def bootstrap_loco_rho(
    raw: np.ndarray,
    truth: np.ndarray,
    fit_predict_fn,
    n_iter: int = 1000,
    seed: int = 42,
    ci_low: float = 2.5,
    ci_high: float = 97.5,
) -> tuple[float, float, float, float]:
    """Out-of-bag bootstrap of the post-cal Spearman rho.

    Delegates to the shared leakage-safe primitive (validation.folding.oob_bootstrap_rho): each
    iteration fits the calibrator on the in-bag resample and scores rho on the OUT-OF-BAG points
    ONLY. The previous version refit on the resample but scored on the FULL original raw array, so
    ~63% of every resample was in-sample -> an optimistic ci_low that fed the ship/escalate Tier
    gate (pure-null data could pass `ci_low > 0`). See docs/BUG_AUDIT_2026-06.md (C1).

    Returns: (point_insample, mean_oob_rho, ci_low, ci_high). The point estimate is the in-sample
    fit-on-full anchor; the bootstrap mean + CI are now out-of-bag.
    """
    from ..validation.folding import oob_bootstrap_rho
    return oob_bootstrap_rho(
        raw, truth, fit_predict_fn,
        n_iter=n_iter, seed=seed, ci_low=ci_low, ci_high=ci_high,
    )


def summarise_metrics(
    raw: np.ndarray,
    truth: np.ndarray,
    loco_pred: np.ndarray,
) -> dict:
    """Compute LOCO RMSE + Spearman ρ + Pearson r. Always returns all keys."""
    mask = ~(np.isnan(raw) | np.isnan(truth) | np.isnan(loco_pred))
    raw_v, truth_v, pred_v = raw[mask], truth[mask], loco_pred[mask]
    out = {
        "loco_rmse": float("nan"), "loco_rho": float("nan"),
        "loco_pearson": float("nan"),
        "raw_rmse": float("nan"), "raw_rho": float("nan"),
        "n_used": int(len(truth_v)),
    }
    if len(truth_v) >= 4:
        out["loco_rmse"] = float(np.sqrt(np.mean((pred_v - truth_v) ** 2)))
        out["loco_rho"] = float(spearmanr(pred_v, truth_v)[0])
        out["loco_pearson"] = float(np.corrcoef(pred_v, truth_v)[0, 1])
        out["raw_rmse"] = float(np.sqrt(np.mean((raw_v - truth_v) ** 2)))
        out["raw_rho"] = float(spearmanr(raw_v, truth_v)[0])
    return out
