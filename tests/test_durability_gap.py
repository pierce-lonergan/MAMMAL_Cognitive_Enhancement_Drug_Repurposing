"""Locks the durability gap (scripts/123): the target cell for a PERMANENT cognitive enhancer is
empty from both directions, and that is what gates the programme.

If any of these assertions ever fails, it means the project has acquired the ground truth it
currently lacks -- which is a MAJOR result and must be reported, not silently absorbed.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "123_durability_gap.py"
ACUTE = ROOT / "data" / "raw" / "healthy_adult_cognition_ledger.csv"
PERSIST = ROOT / "data" / "raw" / "persistence_positive_ledger.csv"
pytestmark = pytest.mark.skipif(not (ACUTE.exists() and PERSIST.exists()), reason="ledgers absent")


def _mod():
    spec = importlib.util.spec_from_file_location("dg", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_acute_ledger_carries_no_durability_information():
    """The ledger that establishes WHO enhances cognition in healthy adults has no post-washout
    column, so it cannot support a durability claim even in principle."""
    import re
    cols = pd.read_csv(ACUTE).columns
    dur = [c for c in cols if re.search(r"washout|durab|persist|retention|follow", c, re.I)]
    assert dur == [], (
        f"acute ledger gained durability columns {dur} -- if real post-washout data has been "
        "curated, scripts/123 and the durability-gap report must be updated")


def test_durable_cognitive_enhancement_in_healthy_adults_has_zero_examples():
    """THE decisive fact. Every durable cognitive entry is a deficit population, where the mechanism
    is restoration of a degraded system rather than enhancement of an intact one."""
    m = _mod()
    per = pd.read_csv(PERSIST)
    per["_is_cog"] = per.apply(m.is_cognition, axis=1)
    per["_healthy"] = ~per.apply(m.classify_population, axis=1).str.startswith("patient")
    target = per[per["_is_cog"] & per["_healthy"]]
    assert len(target) == 0, (
        "the target cell is no longer empty -- a durable cognitive gain in healthy adults would be "
        f"a major finding: {target['compound'].tolist()}")
    # the cell is empty for a substantive reason, not because either axis is empty
    assert per["_is_cog"].sum() >= 1, "there ARE cognition entries (they are just all clinical)"
    assert per["_healthy"].sum() >= 1, "there ARE healthy entries (they are just not cognition)"


def test_population_heuristic_is_conservative_toward_healthy():
    """The classifier only calls an entry healthy when it finds NO clinical marker, so it can only
    OVER-count the healthy cell. The emptiness conclusion is therefore robust to its errors."""
    m = _mod()
    assert m.classify_population({"persistence_design": "randomised in vascular dementia"}
                                 ).startswith("patient")
    assert m.classify_population({"persistence_design": "healthy volunteers, double-blind"}
                                 ) == "no_clinical_marker_found"


def test_every_healthy_brain_mechanism_class_is_symptomatic():
    """The only non-symptomatic class is disease-modifying and patient-scoped: every class acting on
    an intact brain reverses on washout by construction."""
    axis = pd.read_csv(ROOT / "data" / "raw" / "persistence_axis_classes.csv")
    nonsympt = axis[axis["persistence_status"] != "symptomatic"]
    assert (nonsympt["persistence_status"] == "disease_modifying_patients").all(), (
        "a non-symptomatic, non-patient persistence class has appeared -- that would be the first "
        "mechanistic route to durable enhancement in a healthy brain and must be reported")
