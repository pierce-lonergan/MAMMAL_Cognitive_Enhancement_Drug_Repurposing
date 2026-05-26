"""Beta-calibration (Kull, Silva Filho & Flach 2017 AISTATS) wrapper.

⚠️ NOT USABLE for the v3 use case as currently implemented.

Per research/4-tier/Isotonic-PerTarget-Calibration.md §1C, beta-calibration
was proposed as the parametric small-n fallback. However the `betacal`
package (v1.1.0) is internally a binary-classification probability
calibrator (wraps sklearn's LogisticRegression) — it cannot regress on
continuous pchembl targets. The research doc's example code is mistaken.

This module remains as a reference + binary-mode entry point. To use:
either (a) discretise pchembl at 6.0 (1 µM) and treat as binary, or (b)
substitute a different parametric monotone smoother (e.g., a 3-parameter
piecewise-linear monotone spline).

For v1 we ship isotonic-only; beta is deferred to v2.
"""

from __future__ import annotations

import pickle
from dataclasses import dataclass
from pathlib import Path

import numpy as np

try:
    from betacal import BetaCalibration as _BetaCal
    BETACAL_AVAILABLE = True
except ImportError:
    _BetaCal = None
    BETACAL_AVAILABLE = False

from .diagnostics import bootstrap_loco_rho, loco_predictions, summarise_metrics

# pKd range — used to rescale to [0, 1] required by betacal.
PKD_MIN = 2.0
PKD_MAX = 11.0


@dataclass
class BetaCalibrationResult:
    target_uniprot: str
    target_gene: str
    n: int
    parameters: str                     # 'abm' / 'am' / 'ab'
    raw_rho: float
    raw_rmse: float
    loco_rho: float
    loco_rmse: float
    loco_pearson: float
    in_sample_rho: float
    boot_mean_rho: float
    boot_ci_low: float
    boot_ci_high: float
    spans_zero: bool
    rescale_lo: float                   # raw_pkd min (for fwd transform)
    rescale_hi: float                   # raw_pkd max
    pickle_path: str | None = None


def _rescale(x: np.ndarray, lo: float, hi: float) -> np.ndarray:
    """Map raw_pkd to [0, 1] via min-max with calibration-set bounds."""
    return np.clip((x - lo) / (hi - lo + 1e-9), 0.0, 1.0)


def _inv_rescale(s: np.ndarray) -> np.ndarray:
    """Map calibrated [0,1] back to pKd range [PKD_MIN, PKD_MAX]."""
    return s * (PKD_MAX - PKD_MIN) + PKD_MIN


def fit_beta_per_target(
    target_uniprot: str,
    target_gene: str,
    raw_pkd: np.ndarray,
    truth_pchembl: np.ndarray,
    parameters: str = "abm",            # 'abm' (full 3-param) | 'am' | 'ab'
    n_bootstrap: int = 1000,
    pickle_dir: Path | None = None,
) -> BetaCalibrationResult:
    """Fit beta-calibration with diagnostics. Returns BetaCalibrationResult."""
    if not BETACAL_AVAILABLE:
        raise ImportError(
            "betacal not installed. `pip install betacal` to enable."
        )

    raw_pkd = np.asarray(raw_pkd, dtype=float)
    truth_pchembl = np.asarray(truth_pchembl, dtype=float)
    mask = ~(np.isnan(raw_pkd) | np.isnan(truth_pchembl))
    raw_v = raw_pkd[mask]
    truth_v = truth_pchembl[mask]
    n = len(raw_v)

    if n < 4:
        return BetaCalibrationResult(
            target_uniprot=target_uniprot, target_gene=target_gene,
            n=n, parameters=parameters,
            raw_rho=float("nan"), raw_rmse=float("nan"),
            loco_rho=float("nan"), loco_rmse=float("nan"), loco_pearson=float("nan"),
            in_sample_rho=float("nan"),
            boot_mean_rho=float("nan"),
            boot_ci_low=float("nan"), boot_ci_high=float("nan"),
            spans_zero=True,
            rescale_lo=float("nan"), rescale_hi=float("nan"),
        )

    lo, hi = float(raw_v.min() - 0.1), float(raw_v.max() + 0.1)

    def _fit_predict(x_tr, y_tr, x_te):
        s_tr = _rescale(x_tr, lo, hi).reshape(-1, 1)
        y_tr_n = (y_tr - PKD_MIN) / (PKD_MAX - PKD_MIN)
        y_tr_n = np.clip(y_tr_n, 1e-6, 1 - 1e-6)
        bc = _BetaCal(parameters=parameters)
        bc.fit(s_tr, y_tr_n)
        s_te = _rescale(x_te, lo, hi).reshape(-1, 1)
        out_norm = bc.predict(s_te)
        return _inv_rescale(out_norm)

    loco_pred = loco_predictions(raw_v, truth_v, _fit_predict)
    metrics = summarise_metrics(raw_v, truth_v, loco_pred)
    point, boot_mean, ci_lo, ci_hi = bootstrap_loco_rho(
        raw_v, truth_v, _fit_predict, n_iter=n_bootstrap,
    )

    # Pickle the fitted calibrator (full data, no LOCO)
    pickle_path = None
    if pickle_dir is not None:
        bc_full = _BetaCal(parameters=parameters)
        s_full = _rescale(raw_v, lo, hi).reshape(-1, 1)
        y_full_n = np.clip(
            (truth_v - PKD_MIN) / (PKD_MAX - PKD_MIN), 1e-6, 1 - 1e-6,
        )
        bc_full.fit(s_full, y_full_n)
        pickle_dir.mkdir(parents=True, exist_ok=True)
        pickle_path = pickle_dir / f"{target_uniprot}.pkl"
        with open(pickle_path, "wb") as fh:
            pickle.dump({
                "bc": bc_full, "parameters": parameters,
                "rescale_lo": lo, "rescale_hi": hi, "n": n,
            }, fh)

    return BetaCalibrationResult(
        target_uniprot=target_uniprot,
        target_gene=target_gene,
        n=n, parameters=parameters,
        raw_rho=metrics["raw_rho"], raw_rmse=metrics["raw_rmse"],
        loco_rho=metrics["loco_rho"], loco_rmse=metrics["loco_rmse"],
        loco_pearson=metrics["loco_pearson"],
        in_sample_rho=point, boot_mean_rho=boot_mean,
        boot_ci_low=ci_lo, boot_ci_high=ci_hi,
        spans_zero=bool(np.isnan(ci_lo) or ci_lo <= 0),
        rescale_lo=lo, rescale_hi=hi,
        pickle_path=str(pickle_path) if pickle_path else None,
    )


def predict_calibrated_beta(
    payload: dict,
    raw_pkd: np.ndarray | float,
) -> np.ndarray:
    """Apply a pickled beta calibrator. `payload` matches the pickle dict
    layout from fit_beta_per_target."""
    arr = np.atleast_1d(np.asarray(raw_pkd, dtype=float))
    s = _rescale(arr, payload["rescale_lo"], payload["rescale_hi"]).reshape(-1, 1)
    out_norm = payload["bc"].predict(s)
    return _inv_rescale(out_norm)
