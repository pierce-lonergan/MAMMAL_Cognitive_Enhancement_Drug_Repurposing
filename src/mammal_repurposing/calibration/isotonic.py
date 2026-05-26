"""Classical isotonic regression (sklearn) per-target calibrator.

Per research/4-tier/Isotonic-PerTarget-Calibration.md §1A.

Critical design: IsotonicRegression(increasing='auto') runs PAVA in both
directions and selects the lower-SSE fit. This naturally absorbs sign
inversion at MAMMAL_ONLY_INVERTED targets (SLC6A3, SLC6A2) without
requiring an external sign-flip hack.

Below n=15 the 'auto' selection becomes unstable; the router (§7.11
router.py) forces direction from empirical Spearman sign when n<15.
"""

from __future__ import annotations

import pickle
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from sklearn.isotonic import IsotonicRegression

from .diagnostics import bootstrap_loco_rho, loco_predictions, summarise_metrics


# pKd guardrails — anything outside [2, 11] is biologically implausible
# (Y_MIN/Y_MAX clip the isotonic output range).
PKD_MIN = 2.0
PKD_MAX = 11.0


@dataclass
class IsotonicCalibrationResult:
    target_uniprot: str
    target_gene: str
    n: int
    direction: str                       # 'increasing' | 'decreasing' | 'auto'
    raw_rho: float
    raw_rmse: float
    loco_rho: float                      # held-out per-compound
    loco_rmse: float
    loco_pearson: float
    in_sample_rho: float                 # fit-on-full, predict-on-full
    boot_mean_rho: float
    boot_ci_low: float
    boot_ci_high: float
    spans_zero: bool                     # CI 95% lower bound <= 0
    fitted_calibrator: IsotonicRegression | None = None
    pickle_path: str | None = None


def _make_iso(direction: str | bool) -> IsotonicRegression:
    inc = "auto" if direction == "auto" or direction is None else bool(direction)
    return IsotonicRegression(
        increasing=inc,
        out_of_bounds="clip",
        y_min=PKD_MIN,
        y_max=PKD_MAX,
    )


def fit_isotonic_per_target(
    raw_pkd: np.ndarray,
    truth_pchembl: np.ndarray,
    direction: str | bool = "auto",
) -> IsotonicRegression:
    """Fit a single isotonic regressor. Returns the fitted object only."""
    iso = _make_iso(direction)
    iso.fit(raw_pkd, truth_pchembl)
    return iso


def fit_isotonic_with_diagnostics(
    target_uniprot: str,
    target_gene: str,
    raw_pkd: np.ndarray,
    truth_pchembl: np.ndarray,
    direction: str | bool = "auto",
    n_bootstrap: int = 1000,
    pickle_dir: Path | None = None,
) -> IsotonicCalibrationResult:
    """Full diagnostic pipeline: fit + LOCO + bootstrap CI + optional pickle.

    Returns IsotonicCalibrationResult ready to populate the router's
    `weights_calibrated.yaml` schema.
    """
    raw_pkd = np.asarray(raw_pkd, dtype=float)
    truth_pchembl = np.asarray(truth_pchembl, dtype=float)
    mask = ~(np.isnan(raw_pkd) | np.isnan(truth_pchembl))
    raw_v = raw_pkd[mask]
    truth_v = truth_pchembl[mask]
    n = len(raw_v)

    if n < 4:
        return IsotonicCalibrationResult(
            target_uniprot=target_uniprot, target_gene=target_gene,
            n=n, direction=str(direction),
            raw_rho=float("nan"), raw_rmse=float("nan"),
            loco_rho=float("nan"), loco_rmse=float("nan"), loco_pearson=float("nan"),
            in_sample_rho=float("nan"),
            boot_mean_rho=float("nan"),
            boot_ci_low=float("nan"), boot_ci_high=float("nan"),
            spans_zero=True, fitted_calibrator=None,
        )

    iso = _make_iso(direction)
    iso.fit(raw_v, truth_v)

    # Inferred direction after auto — read sklearn's attribute
    inferred = "increasing" if getattr(iso, "increasing_", True) else "decreasing"
    if direction != "auto" and direction is not None:
        inferred = "increasing" if direction is True else "decreasing"

    # LOCO predictions
    def _fit_predict(x_tr, y_tr, x_te):
        iso_i = _make_iso(direction)
        iso_i.fit(x_tr, y_tr)
        return iso_i.predict(x_te)

    loco_pred = loco_predictions(raw_v, truth_v, _fit_predict)
    metrics = summarise_metrics(raw_v, truth_v, loco_pred)

    # Bootstrap CI
    point, boot_mean, ci_lo, ci_hi = bootstrap_loco_rho(
        raw_v, truth_v, _fit_predict, n_iter=n_bootstrap,
    )

    # Optional pickle
    pickle_path = None
    if pickle_dir is not None:
        pickle_dir.mkdir(parents=True, exist_ok=True)
        pickle_path = pickle_dir / f"{target_uniprot}.pkl"
        with open(pickle_path, "wb") as fh:
            pickle.dump({
                "iso": iso, "direction": inferred,
                "n": n, "raw_min": float(raw_v.min()),
                "raw_max": float(raw_v.max()),
            }, fh)

    return IsotonicCalibrationResult(
        target_uniprot=target_uniprot,
        target_gene=target_gene,
        n=n,
        direction=inferred,
        raw_rho=metrics["raw_rho"],
        raw_rmse=metrics["raw_rmse"],
        loco_rho=metrics["loco_rho"],
        loco_rmse=metrics["loco_rmse"],
        loco_pearson=metrics["loco_pearson"],
        in_sample_rho=point,
        boot_mean_rho=boot_mean,
        boot_ci_low=ci_lo,
        boot_ci_high=ci_hi,
        spans_zero=bool(np.isnan(ci_lo) or ci_lo <= 0),
        fitted_calibrator=iso,
        pickle_path=str(pickle_path) if pickle_path else None,
    )


def predict_calibrated(
    iso: IsotonicRegression,
    raw_pkd: np.ndarray | float,
) -> np.ndarray:
    """Apply a fitted isotonic calibrator. Wraps sklearn predict for the
    common single-value case."""
    arr = np.atleast_1d(np.asarray(raw_pkd, dtype=float))
    return iso.predict(arr)
