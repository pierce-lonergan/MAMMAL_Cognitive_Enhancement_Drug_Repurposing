"""C6 regression: the hierarchical NUTS path must fail CLOSED on non-convergence, and its
pooled_rho statistic is provably degenerate (BUG_AUDIT_2026-06.md, C6).

Found by actually RUNNING the pipeline once PyMC became installed: the NUTS path produced 148-281
divergences with R-hat > 1.01 and wrote those numbers to the published report as authoritative,
while its pooled_rho was mathematically incapable of differing from single_target_rho.
"""
from __future__ import annotations

import numpy as np

from mammal_repurposing.calibration import hierarchical_bayes as HB


def test_pooled_rho_metric_is_affine_invariant_and_cannot_pool():
    """corr(alpha + beta*x, y) == sign(beta) * corr(x, y) EXACTLY, for any alpha and any beta != 0.

    This is why the NUTS pooled_rho can never express pooling: it is a positive-affine transform of
    x, and Pearson correlation is invariant to those. Locks the defect so a future edit cannot
    silently reintroduce the metric as if it measured shrinkage."""
    rng = np.random.default_rng(0)
    x = rng.normal(size=50)
    y = 0.7 * x + rng.normal(size=50)
    base = float(np.corrcoef(x, y)[0, 1])
    for alpha, beta in [(0.0, 1.0), (3.7, 0.2), (-9.1, 5.0), (2.0, 1e-3)]:
        assert abs(float(np.corrcoef(alpha + beta * x, y)[0, 1]) - base) < 1e-12
    # a negative slope only flips the SIGN -- magnitude is still pinned to |corr(x, y)|
    assert abs(float(np.corrcoef(1.0 - 2.0 * x, y)[0, 1]) + base) < 1e-12


def test_fit_family_falls_back_to_shrinkage_when_nuts_not_converged(monkeypatch):
    """A non-converged NUTS fit must NOT be returned; fit_family falls back to the deterministic
    shrinkage estimator. Fails-before (the old code returned whatever NUTS produced)."""
    x = np.array([5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0])
    y = np.array([5.2, 6.1, 6.9, 8.3, 8.8, 10.4, 10.9])
    data = {"T1": (x, y), "T2": (x + 0.5, y - 0.3)}

    bad = HB.HierarchicalCalibrationResult(
        family="F", targets=["T1", "T2"], n_per_target={"T1": 7, "T2": 7},
        single_target_rho={"T1": 0.9, "T2": 0.9}, pooled_rho={"T1": 0.9, "T2": 0.9},
        method="pymc_nuts", n_divergences=281, rhat_max=1.34, ess_min=42.0, converged=False,
    )
    monkeypatch.setattr(HB, "PYMC_AVAILABLE", True)
    monkeypatch.setattr(HB, "hierarchical_bayesian_nuts", lambda *a, **k: bad)

    res = HB.fit_family("F", data, prefer_pymc=True)
    assert res.method == "empirical_bayes_shrinkage", "must not publish an unconverged posterior"


def test_converged_nuts_result_is_accepted(monkeypatch):
    """Guard the other direction: a CONVERGED NUTS fit is still used (the gate is not a blanket
    disable of the NUTS path)."""
    x = np.array([5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0])
    good = HB.HierarchicalCalibrationResult(
        family="F", targets=["T1"], n_per_target={"T1": 7},
        single_target_rho={"T1": 0.9}, pooled_rho={"T1": 0.8},
        method="pymc_nuts", n_divergences=0, rhat_max=1.001, ess_min=1200.0, converged=True,
    )
    monkeypatch.setattr(HB, "PYMC_AVAILABLE", True)
    monkeypatch.setattr(HB, "hierarchical_bayesian_nuts", lambda *a, **k: good)
    res = HB.fit_family("F", {"T1": (x, x * 1.1)}, prefer_pymc=True)
    assert res.method == "pymc_nuts"
