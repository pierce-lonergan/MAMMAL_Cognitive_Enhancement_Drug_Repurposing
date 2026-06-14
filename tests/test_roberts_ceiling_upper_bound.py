"""C5 regression: roberts_2020_ceiling_check must use the `upper_quantile` credible UPPER BOUND
when a per-target SD is supplied, not a point comparison (BUG_AUDIT_2026-06.md, C5).
"""
from __future__ import annotations

from mammal_repurposing.cluster_d.bayesian_prior import (
    fit_cluster_d_prior_stub,
    roberts_2020_ceiling_check,
)


def test_roberts_ceiling_uses_upper_quantile_when_sd_given():
    post = fit_cluster_d_prior_stub(["T1"])
    preds = {"T1": 0.40}          # point is BELOW the 0.5 ceiling...
    sds = {"T1": 0.20}            # ...but 90% upper bound = 0.40 + 1.2816*0.20 = 0.656 > 0.5
    out = roberts_2020_ceiling_check(post, preds, target_smd_sd=sds)
    # FAILS on the pre-fix code (it ignored sd/upper_quantile -> REGIME_OK); PASSES post-fix.
    assert out["T1"] == "REGIME_VIOLATION"


def test_roberts_ceiling_point_only_is_backward_compatible():
    post = fit_cluster_d_prior_stub(["T1"])
    # Without an SD the check degrades to the point comparison (0.40 < 0.5 -> OK); existing
    # point-only callers/tests are unchanged.
    assert roberts_2020_ceiling_check(post, {"T1": 0.40})["T1"] == "REGIME_OK"
    assert roberts_2020_ceiling_check(post, {"T1": 0.70})["T1"] == "REGIME_VIOLATION"
    assert roberts_2020_ceiling_check(post, None)["T1"] == "NO_SMD_PREDICTION"
