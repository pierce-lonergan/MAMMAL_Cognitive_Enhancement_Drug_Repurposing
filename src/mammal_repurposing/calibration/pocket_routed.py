"""§8.14 — Pocket-class-routed isotonic calibration.

For targets with multiple pocket classes (SLC6A3 has S1 orthosteric +
S2 vestibule allosteric; ACHE has CAS + PAS; CHRNA7 has Type-I + Type-II
PAMs), the single per-target isotonic calibrator from §7.11 assumes one
monotone MAMMAL-vs-truth relationship. If the relationship differs by
pocket class, this assumption is violated and the calibrator becomes a
weighted average that fits neither sub-class well.

This module fits per-(target, pocket_class) isotonic calibrators when pose
data is available, falls back to the per-target calibrator otherwise.

Reference:
  Pocket-Conditioned-Boltz2.md §6.4 — pocket-routed calibration spec.
  Cheng et al. 2020 J Chem Inf Model — DAT S2 vestibule discovery.
  Nielsen et al. 2024 Nature — hDAT cocaine cryo-EM at S1.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.isotonic import IsotonicRegression

from .isotonic import _make_iso, PKD_MIN, PKD_MAX

logger = logging.getLogger(__name__)


@dataclass
class PocketRoutedCalibrator:
    target_uniprot: str
    by_pocket_class: dict[str, IsotonicRegression]
    fallback: IsotonicRegression | None
    n_by_pocket: dict[str, int]
    raw_min_by_pocket: dict[str, float]
    raw_max_by_pocket: dict[str, float]


def fit_pocket_routed(
    target_uniprot: str,
    raw_pkd: np.ndarray,
    truth_pchembl: np.ndarray,
    pocket_classes: np.ndarray,
    min_n_per_pocket: int = 5,
    direction: str | bool = "auto",
) -> PocketRoutedCalibrator:
    """Fit a separate isotonic per pocket_class; fall back to a global
    isotonic when a pocket class has < `min_n_per_pocket` samples.

    Args:
        raw_pkd, truth_pchembl: aligned arrays of length n.
        pocket_classes: per-sample pocket class label (e.g. "orthosteric",
            "allosteric_known"). Pass empty strings when unknown.
        min_n_per_pocket: below this, fall back to the global calibrator.
    """
    raw = np.asarray(raw_pkd, dtype=float)
    truth = np.asarray(truth_pchembl, dtype=float)
    pcls = np.asarray(pocket_classes)
    by_pocket: dict[str, IsotonicRegression] = {}
    n_by_pocket: dict[str, int] = {}
    raw_min: dict[str, float] = {}
    raw_max: dict[str, float] = {}

    # Per-pocket fits
    for cls in sorted(set(pcls.tolist())):
        if not cls:
            continue
        mask = pcls == cls
        n_cls = int(mask.sum())
        n_by_pocket[cls] = n_cls
        if n_cls < min_n_per_pocket:
            logger.info("  %s/%s: only %d samples; falling back", target_uniprot, cls, n_cls)
            continue
        iso = _make_iso(direction)
        iso.fit(raw[mask], truth[mask])
        by_pocket[cls] = iso
        raw_min[cls] = float(raw[mask].min())
        raw_max[cls] = float(raw[mask].max())
        logger.info("  %s/%s: fit isotonic on n=%d", target_uniprot, cls, n_cls)

    # Global fallback always available
    fallback = _make_iso(direction)
    fallback.fit(raw, truth)

    return PocketRoutedCalibrator(
        target_uniprot=target_uniprot,
        by_pocket_class=by_pocket,
        fallback=fallback,
        n_by_pocket=n_by_pocket,
        raw_min_by_pocket=raw_min,
        raw_max_by_pocket=raw_max,
    )


def predict_with_routing(
    calibrator: PocketRoutedCalibrator,
    raw_pkd: np.ndarray,
    pocket_classes: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Returns (predictions, used_calibrator_tag) where used_calibrator_tag
    is "pocket:<cls>" if a pocket-specific isotonic was used or "fallback"."""
    raw = np.asarray(raw_pkd, dtype=float)
    pcls = np.asarray(pocket_classes)
    out = np.empty_like(raw)
    tags = np.empty(len(raw), dtype=object)
    for i, (x, c) in enumerate(zip(raw, pcls)):
        if c and c in calibrator.by_pocket_class:
            out[i] = float(calibrator.by_pocket_class[c].predict([x])[0])
            tags[i] = f"pocket:{c}"
        elif calibrator.fallback is not None:
            out[i] = float(calibrator.fallback.predict([x])[0])
            tags[i] = "fallback"
        else:
            out[i] = float("nan")
            tags[i] = "no_calibrator"
    return out, tags


def evaluate_routing_lift(
    raw_pkd: np.ndarray,
    truth_pchembl: np.ndarray,
    pocket_classes: np.ndarray,
    direction: str | bool = "auto",
) -> dict:
    """Measure whether per-pocket isotonic outperforms global isotonic on
    a sum-of-squared-residuals basis. Hypothesis: routing gives lower SSR
    when the per-pocket relationships truly differ.

    Returns dict with: global_ssr, routed_ssr, lift_pct.
    """
    raw = np.asarray(raw_pkd, dtype=float)
    truth = np.asarray(truth_pchembl, dtype=float)
    pcls = np.asarray(pocket_classes)
    if len(raw) < 10:
        return {"global_ssr": float("nan"), "routed_ssr": float("nan"),
                "lift_pct": float("nan"),
                "note": "n < 10; insufficient"}
    # Global
    global_iso = _make_iso(direction)
    global_iso.fit(raw, truth)
    global_pred = global_iso.predict(raw)
    global_ssr = float(np.sum((truth - global_pred) ** 2))
    # Routed (with same min_n_per_pocket policy)
    cal = fit_pocket_routed("__eval__", raw, truth, pcls, direction=direction)
    routed_pred, _ = predict_with_routing(cal, raw, pcls)
    routed_ssr = float(np.sum((truth - routed_pred) ** 2))
    lift_pct = 100.0 * (global_ssr - routed_ssr) / max(global_ssr, 1e-9)
    return {
        "global_ssr": global_ssr,
        "routed_ssr": routed_ssr,
        "lift_pct": lift_pct,
        "n_total": int(len(raw)),
        "n_by_pocket": cal.n_by_pocket,
    }
