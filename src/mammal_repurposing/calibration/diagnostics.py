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
    """Bootstrap the post-cal Spearman ρ.

    For each bootstrap iteration:
      - resample with replacement, refit the calibrator on the resample,
        predict on the original raw values, compute Spearman ρ vs truth.

    Returns: (point_estimate, mean_boot_rho, ci_low, ci_high). Point estimate
    is the on-fit-data ρ; bootstrap statistics characterise its variability.
    """
    raw = np.asarray(raw, dtype=float)
    truth = np.asarray(truth, dtype=float)
    n = len(raw)
    if n < 4:
        return float("nan"), float("nan"), float("nan"), float("nan")

    # Point estimate via fit-on-full + predict-on-full (in-sample ρ)
    try:
        pred_full = fit_predict_fn(raw, truth, raw)
        point = float(spearmanr(pred_full, truth)[0])
    except Exception:
        point = float("nan")

    rng = np.random.default_rng(seed)
    rhos = []
    for _ in range(n_iter):
        idx = rng.integers(0, n, n)
        if len(np.unique(idx)) < 4:
            continue
        try:
            pred = fit_predict_fn(raw[idx], truth[idx], raw)
            r, _ = spearmanr(pred, truth)
            if not np.isnan(r):
                rhos.append(r)
        except Exception:
            continue

    if not rhos:
        return point, float("nan"), float("nan"), float("nan")

    return (
        point,
        float(np.mean(rhos)),
        float(np.percentile(rhos, ci_low)),
        float(np.percentile(rhos, ci_high)),
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
