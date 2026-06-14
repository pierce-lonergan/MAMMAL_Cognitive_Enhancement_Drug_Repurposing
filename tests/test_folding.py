"""Regression tests for the leakage-safe folding primitives (docs/BUG_AUDIT_2026-06.md).

The decisive property the contaminated evaluation paths lacked: on RANDOM labels, an honest
out-of-bag / fit-on-train metric must NOT declare signal. These tests lock that in.
"""
from __future__ import annotations

import numpy as np

from mammal_repurposing.validation.folding import (
    apply_quantile_edges,
    fit_quantile_edges,
    oob_bootstrap_rho,
)


def _identity(x_tr, y_tr, x_eval):
    """A trivial monotone 'calibrator' that returns the predictor unchanged."""
    return np.asarray(x_eval, dtype=float)


def _overfit(x_tr, y_tr, x_eval):
    """A calibrator that memorizes the train mapping (nearest-neighbour in x). In-sample this
    reproduces y exactly (rho=1); out-of-bag it cannot, so it exposes optimistic in-sample CIs."""
    x_tr = np.asarray(x_tr, dtype=float)
    y_tr = np.asarray(y_tr, dtype=float)
    x_eval = np.asarray(x_eval, dtype=float)
    idx = np.array([int(np.argmin(np.abs(x_tr - v))) for v in x_eval])
    return y_tr[idx]


def test_oob_bootstrap_null_labels_ci_spans_zero():
    rng = np.random.default_rng(0)
    raw = rng.normal(size=80)
    truth = rng.normal(size=80)   # independent of raw -> no real signal
    _, _, ci_lo, ci_hi = oob_bootstrap_rho(raw, truth, _identity, n_iter=600, seed=1)
    assert ci_lo <= 0.0 <= ci_hi, f"null OOB CI must span 0, got [{ci_lo:.3f}, {ci_hi:.3f}]"


def test_oob_bootstrap_overfit_calibrator_does_not_pass_on_noise():
    """The key contamination guard: a memorizing calibrator scores rho=1 IN-SAMPLE on random data,
    but its OUT-OF-BAG CI must still span 0 (it has no real signal to generalize)."""
    rng = np.random.default_rng(3)
    raw = rng.normal(size=80)
    truth = rng.normal(size=80)
    point, _, ci_lo, _ = oob_bootstrap_rho(raw, truth, _overfit, n_iter=600, seed=2)
    assert point > 0.9          # in-sample anchor is inflated by the memorization
    assert ci_lo <= 0.0         # but the honest OOB CI does not declare signal


def test_oob_bootstrap_detects_real_signal():
    rng = np.random.default_rng(0)
    raw = rng.normal(size=100)
    truth = raw + 0.3 * rng.normal(size=100)   # strong monotone signal
    _, _, ci_lo, _ = oob_bootstrap_rho(raw, truth, _identity, n_iter=600, seed=1)
    assert ci_lo > 0.0, "real monotone signal should give an OOB CI excluding 0"


def test_quantile_edges_fit_on_train_apply_to_test():
    train = np.arange(100, dtype=float)
    edges = fit_quantile_edges(train, n_bins=5)
    assert len(edges) == 4
    buckets = apply_quantile_edges(np.array([-5.0, 10.0, 50.0, 95.0, 200.0]), edges)
    assert buckets.min() == 0 and buckets.max() == 4
    assert list(buckets) == sorted(buckets)   # monotone in the input


def test_quantile_edges_handle_nan_and_degenerate():
    assert fit_quantile_edges(np.array([]), 5).size == 0
    assert fit_quantile_edges(np.array([1.0, 2.0]), 1).size == 0
    # NaN inputs map to the lowest bucket, never raise
    out = apply_quantile_edges(np.array([np.nan, 1.0]), np.array([0.5]))
    assert out[0] == 0
