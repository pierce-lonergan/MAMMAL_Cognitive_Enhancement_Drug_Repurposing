"""Tests for the persistence-target DTI module (PERSEUS roadmap #2).

The calibration math (AUROC / permutation-p / Youden) and the substrate-hypothesis
aggregation are pure-Python, so these run with no GPU. The honesty-critical assertions are:
calibration GATES the substrate read, and only an ABLATIVE engagement promotes to durable.
"""
from __future__ import annotations

import math
from pathlib import Path

import pytest

from mammal_repurposing.engine.persistence_dti import (
    PanelTarget, auroc, calibrate_target, load_panel, permutation_p,
    substrate_hypothesis, youden_threshold,
)

ROOT = Path(__file__).resolve().parents[1]
INTERIM = ROOT / "data" / "interim"


# --------------------------------------------------------------------------
# calibration math
# --------------------------------------------------------------------------

def test_auroc_separable_reversed_tied_and_nan():
    assert auroc([8, 9, 10], [1, 2, 3]) == 1.0
    assert auroc([1, 2, 3], [8, 9, 10]) == 0.0
    assert auroc([5, 5, 5], [5, 5, 5]) == 0.5          # all ties -> 0.5
    assert auroc([9, float("nan"), 8], [1, 2]) == 1.0  # NaN dropped, still separable
    assert math.isnan(auroc([], [1, 2]))


def test_youden_threshold_between_classes():
    t = youden_threshold([8, 9, 10], [1, 2, 3])
    assert 3 < t <= 8     # a clean cut sits above the negatives


def test_calibrate_pass_fail_and_min_pos():
    # clean separation, n_pos>=3 -> PASS, and engagers clear the specificity threshold so
    # the channel is usable for single-compound inference
    good = calibrate_target([8, 8.5, 9, 9.5], [1, 2, 3, 4, 5])
    assert good["passed"] is True and good["auroc"] == 1.0 and good["perm_p"] < 0.05
    assert good["inference_usable"] is True and good["sensitivity_at_threshold"] == 1.0
    # fully overlapping -> fail both gates
    bad = calibrate_target([4, 5, 6], [4, 5, 6])
    assert bad["passed"] is False and bad["inference_usable"] is False
    # high AUROC but too few positives -> fail on the n_pos guard
    thin = calibrate_target([9, 10], [1, 2, 3, 4])
    assert thin["passed"] is False and thin["n_pos_scored"] == 2


def test_ranks_well_but_not_inference_usable():
    # engagers out-RANK negatives (good AUROC) but sit INSIDE the negative band, so a
    # specificity-first single-compound threshold cannot separate them: passed yet unusable.
    # negatives 1..20 (so the 0.95 quantile ~ 19); engagers 15..18 rank above ~80% of the
    # negatives but none clear the top-5% negative tail.
    neg = list(range(1, 21))
    pos = [15, 16, 17, 18]
    cal = calibrate_target(pos, neg)
    assert cal["auroc"] > 0.7 and cal["passed"] is True   # ranking gate PASSES
    assert cal["sensitivity_at_threshold"] == 0.0
    assert cal["inference_usable"] is False               # but unusable per-compound


def test_permutation_p_is_high_for_no_separation():
    p = permutation_p([5, 5, 5], [5, 5, 5], n_perm=500)
    assert p > 0.5   # no real separation -> not significant


# --------------------------------------------------------------------------
# substrate-hypothesis aggregation (the head)
# --------------------------------------------------------------------------

_PANEL = {
    "BCL2": PanelTarget("BCL2", "P10415", "", "ablative", True),
    "HDAC1": PanelTarget("HDAC1", "Q13547", "", "capability", False),
    "NTRK2": PanelTarget("NTRK2", "Q16620", "", "plasticity_window", False),
}
# BCL2 + HDAC1 calibrated + usable; NTRK2 NOT calibrated (passed False)
_CALIB = {
    "BCL2": {"threshold": 6.0, "passed": True, "inference_usable": True},
    "HDAC1": {"threshold": 6.0, "passed": True, "inference_usable": True},
    "NTRK2": {"threshold": 6.0, "passed": False, "inference_usable": False},
}


def test_ablative_engagement_promotes_durable():
    h = substrate_hypothesis({"BCL2": 7.5, "HDAC1": 3.0, "NTRK2": 2.0}, _PANEL, _CALIB)
    assert h.promotes_durable is True
    assert h.substrate_hypothesis == "ablative"
    assert [e["gene"] for e in h.engaged] == ["BCL2"]


def test_capability_engagement_is_flag_not_durable():
    h = substrate_hypothesis({"BCL2": 2.0, "HDAC1": 8.0, "NTRK2": 2.0}, _PANEL, _CALIB)
    assert h.promotes_durable is False
    assert h.substrate_hypothesis == "capability"
    assert "capability:HDAC1" in h.capability_flags


def test_uncalibrated_target_is_abstained_not_trusted():
    # NTRK2 "looks engaged" (pkd above threshold) but its channel failed calibration ->
    # it must be IGNORED, contributing nothing, and recorded as abstained.
    h = substrate_hypothesis({"BCL2": 2.0, "HDAC1": 2.0, "NTRK2": 9.0}, _PANEL, _CALIB)
    assert h.substrate_hypothesis is None and h.promotes_durable is False
    assert h.abstained_targets == ["NTRK2"]
    assert "UN-calibrated" in h.note


def test_passed_but_not_inference_usable_is_abstained():
    # a channel can pass the AUROC/ranking gate yet be unusable for single-compound calls;
    # a crossing on such a channel must be IGNORED (abstained), never trusted.
    panel = {"HDAC1": PanelTarget("HDAC1", "Q13547", "", "capability", False)}
    calib = {"HDAC1": {"threshold": 6.0, "passed": True, "inference_usable": False}}
    h = substrate_hypothesis({"HDAC1": 9.0}, panel, calib)
    assert h.substrate_hypothesis is None and h.abstained_targets == ["HDAC1"]


def test_below_threshold_is_not_engaged():
    h = substrate_hypothesis({"BCL2": 5.0, "HDAC1": 5.5}, _PANEL, _CALIB)  # both < 6.0
    assert h.substrate_hypothesis is None and not h.engaged


def test_ablative_wins_when_both_engaged():
    # both an ablative and a capability target engage -> top tier is ablative, durable True
    h = substrate_hypothesis({"BCL2": 7.0, "HDAC1": 7.0}, _PANEL, _CALIB)
    assert h.substrate_hypothesis == "ablative" and h.promotes_durable is True
    assert "capability:HDAC1" in h.capability_flags   # still flagged as a capability


# --------------------------------------------------------------------------
# real panel loads
# --------------------------------------------------------------------------

@pytest.mark.skipif(not (INTERIM / "persistence_targets.csv").exists(),
                    reason="persistence target panel not fetched")
def test_real_panel_tiers():
    panel = load_panel(INTERIM / "persistence_targets.csv")
    assert {"BCL2", "BCL2L1", "HDAC1", "DNMT1", "KEAP1", "NTRK2"} <= set(panel)
    assert panel["BCL2"].tier == "ablative" and panel["BCL2"].promotes_durable is True
    assert panel["HDAC1"].tier == "capability" and panel["HDAC1"].promotes_durable is False
    assert panel["NTRK2"].tier == "plasticity_window"
