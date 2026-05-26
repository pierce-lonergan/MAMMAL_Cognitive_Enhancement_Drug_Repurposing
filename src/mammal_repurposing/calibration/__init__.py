"""§7.11 — Per-target post-hoc calibration for MAMMAL DTI predictions.

Per research/4-tier/Isotonic-PerTarget-Calibration.md.

Modules:
  isotonic.py     — sklearn IsotonicRegression with LOCO + bootstrap CI
  beta_cal.py     — betacal.BetaCalibration wrapper (small-n parametric)
  router.py       — decision tree per §1D of the research doc
  diagnostics.py  — LOCO ρ, Brier, reliability metrics
  (hierarchical.py — PyMC SLC6/GRIN pool, deferred)
  (venn_abers.py  — binary-discretised mode, deferred)

The elegant insight: IsotonicRegression(increasing='auto') naturally
absorbs sign inversion at MAMMAL_ONLY_INVERTED targets, replacing the
awkward weight=0.30 hack in weights_calibrated.yaml.
"""

from .isotonic import (
    IsotonicCalibrationResult,
    fit_isotonic_per_target,
    fit_isotonic_with_diagnostics,
)
from .beta_cal import BetaCalibrationResult, fit_beta_per_target
from .router import (
    CalibratorChoice,
    decide_calibrator,
    DECISION_MATRIX,
)
from .diagnostics import bootstrap_loco_rho, loco_predictions

__all__ = [
    "IsotonicCalibrationResult",
    "fit_isotonic_per_target",
    "fit_isotonic_with_diagnostics",
    "BetaCalibrationResult",
    "fit_beta_per_target",
    "CalibratorChoice",
    "decide_calibrator",
    "DECISION_MATRIX",
    "bootstrap_loco_rho",
    "loco_predictions",
]
