"""Leakage-safe folding primitives (see docs/BUG_AUDIT_2026-06.md, systemic finding).

Several evaluation/calibration paths computed their self-assessment statistics on IN-SAMPLE
(pre-split) data -- discretizing labels, imputing, or scoring on the full array BEFORE the
train/test or bootstrap split. That inflates the metric optimistically and, worse, feeds the
automated ship/escalate gates and published validation numbers. This module centralizes the
"fit on TRAIN, apply/score on HELD-OUT" contract so every evaluation folds the same correct way:

  - fit_quantile_edges / apply_quantile_edges : discretize using edges fit on train only.
  - oob_bootstrap_rho                          : bootstrap a Spearman rho scored OUT-OF-BAG only.

The companion regression test (tests/test_folding.py) feeds RANDOM labels and asserts each of these
collapses to ~0 / a zero-spanning CI -- the property the contaminated versions lacked (a flexible
calibrator scored in-sample could report a confident positive rho on pure noise).
"""
from __future__ import annotations

import numpy as np
from scipy.stats import spearmanr


def fit_quantile_edges(values: np.ndarray, n_bins: int) -> np.ndarray:
    """Interior quantile cut-points fit on TRAINING values only.

    Returns the (n_bins - 1) interior quantile edges (deduplicated). Pair with
    apply_quantile_edges so a train/test discretization never lets the held-out labels define the
    training buckets. Non-finite values are dropped before fitting.
    """
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0 or n_bins < 2:
        return np.array([], dtype=float)
    qs = np.linspace(0.0, 1.0, n_bins + 1)[1:-1]   # interior quantiles only
    return np.unique(np.quantile(values, qs))


def apply_quantile_edges(values: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """Map values to integer buckets 0..len(edges) using edges from fit_quantile_edges.

    NaN maps to bucket 0 (lowest) so integer-label consumers (e.g. LightGBM gain labels) never see
    a NaN. With k interior edges this yields k+1 monotone buckets.
    """
    values = np.asarray(values, dtype=float)
    edges = np.asarray(edges, dtype=float)
    buckets = np.digitize(np.nan_to_num(values, nan=-np.inf), edges)
    return buckets.astype(int)


def oob_bootstrap_rho(
    raw: np.ndarray,
    truth: np.ndarray,
    fit_predict_fn,
    *,
    n_iter: int = 1000,
    seed: int = 42,
    ci_low: float = 2.5,
    ci_high: float = 97.5,
    min_oob: int = 4,
) -> tuple[float, float, float, float]:
    """Out-of-bag bootstrap of post-calibration Spearman rho.

    Each iteration draws an in-bag resample (with replacement), fits
    ``fit_predict_fn(raw[inbag], truth[inbag], raw[oob])`` and scores Spearman rho on the
    OUT-OF-BAG points ONLY -- never on the points used to fit. The contaminated predecessor scored
    on the FULL original array, so ~63% of every "bootstrap" sample was in-sample; that inflated
    ``ci_low`` and let pure-null data pass a ``ci_low > 0`` ship gate.

    Args:
        raw, truth: 1D arrays (predictor, ground truth).
        fit_predict_fn: callable (X_train, y_train, X_eval) -> y_pred[len(X_eval)].
    Returns:
        (point_insample, mean_oob_rho, ci_low_pct, ci_high_pct). The point estimate is the
        in-sample fit-on-full anchor; the mean and CI are out-of-bag.
    """
    raw = np.asarray(raw, dtype=float)
    truth = np.asarray(truth, dtype=float)
    n = len(raw)
    if n < 4:
        return float("nan"), float("nan"), float("nan"), float("nan")

    try:
        point = float(spearmanr(fit_predict_fn(raw, truth, raw), truth)[0])
    except Exception:
        point = float("nan")

    rng = np.random.default_rng(seed)
    all_idx = np.arange(n)
    rhos: list[float] = []
    for _ in range(n_iter):
        inbag = rng.integers(0, n, n)
        uniq = np.unique(inbag)
        oob = np.setdiff1d(all_idx, uniq, assume_unique=True)
        if oob.size < min_oob or uniq.size < min_oob:
            continue
        try:
            pred = fit_predict_fn(raw[inbag], truth[inbag], raw[oob])
            r, _ = spearmanr(pred, truth[oob])
            if not np.isnan(r):
                rhos.append(float(r))
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
