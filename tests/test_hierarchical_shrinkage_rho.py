"""B7 regression: fit_family's shipped SHRINKAGE path must compute single_target_rho with Spearman
(rank) correlation -- the framework convention -- not Pearson (BUG_AUDIT_2026-06.md, B7).
"""
from __future__ import annotations

import numpy as np
from scipy.stats import spearmanr

from mammal_repurposing.calibration.hierarchical_bayes import fit_family


def test_shrinkage_single_rho_is_spearman_not_pearson():
    x = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
    y = x ** 3                                       # strictly monotone, nonlinear
    res = fit_family("F", {"T": (x, y)}, prefer_pymc=False)   # force the shrinkage path
    sp = float(spearmanr(x, y)[0])                   # == 1.0 (perfect rank agreement)
    pe = float(np.corrcoef(x, y)[0, 1])              # < 1.0 (nonlinear)
    assert abs(sp - pe) > 0.02, "test data must distinguish Spearman from Pearson"
    # FAILS on the pre-fix code (Pearson pe ~0.93); PASSES post-fix (Spearman 1.0).
    assert abs(res.single_target_rho["T"] - sp) < 1e-9
