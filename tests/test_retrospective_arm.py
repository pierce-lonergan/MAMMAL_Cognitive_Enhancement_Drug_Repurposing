"""Invariants for the retrospective arm.

The arm exists because the prospective registry cannot be scored. Its whole value depends on three
things staying true: that it is never pooled with the prospective arm, that nothing enters it
without the attribution check the retraction taught, and that a single row is never allowed to look
like evidence.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
ARM = ROOT / "data" / "raw" / "retrospective_arm.csv"
PROSPECTIVE = ROOT / "data" / "raw" / "prospective_healthy_adult.csv"


@pytest.fixture(scope="module")
def arm():
    return pd.read_csv(ARM)


def test_the_two_arms_share_no_registration(arm):
    """Pooling a retrodiction with a prediction is the error that made the patient arm's
    '2/2 correct, genuine out-of-sample' claim wrong. They must not even overlap."""
    pro = pd.read_csv(PROSPECTIVE)
    overlap = set(arm["identifier"].astype(str)) & set(pro["identifier"].astype(str))
    assert not overlap, f"same registration in both arms: {overlap}"


def test_the_retracted_registration_is_not_here(arm):
    """NCT07469852's paper belonged to a different study. It must not reappear via this door."""
    assert "NCT07469852" not in set(arm["identifier"].astype(str))


def test_no_systematic_review_is_used_as_a_test_of_a_meta_analytic_estimate(arm):
    """Two reviews of the same literature re-pool overlapping primary studies, so a new
    meta-analysis is not an independent observation of a meta-analytic estimate."""
    ids = arm["identifier"].astype(str)
    assert not ids.str.startswith(("CRD", "INPLASY")).any()


def test_every_scoreable_row_states_how_it_was_attributed(arm):
    """Confirming a source exists is not confirming it belongs to the registration."""
    sco = arm[arm["scoreable"]]
    assert len(sco) >= 1
    for _, r in sco.iterrows():
        a = str(r["attribution"])
        assert len(a) > 40, "attribution must be argued, not asserted"
        assert str(r["source"]).strip()
        assert str(r["quote"]).strip(), "a scoreable row must carry the verbatim result"


def test_unscoreable_rows_carry_no_outcome_and_a_reason(arm):
    uns = arm[~arm["scoreable"]]
    assert uns["actual_outcome"].fillna("").astype(str).str.strip().replace("nan", "").eq("").all()
    assert uns["caveat"].astype(str).str.len().gt(40).all()


def test_a_single_row_cannot_claim_to_beat_the_constant(arm):
    """The arm's one scoreable row favours the rule. At n = 1 that must not reach significance."""
    from mammal_repurposing.reporting.prospective import score_healthy_adult
    sc = score_healthy_adult(arm)
    assert sc["n_scored"] >= 1
    assert sc["beats_constant"] is False, "n=1 must never clear the bar"
    assert sc["p_value"] >= 0.05


def test_posted_results_rows_are_unresolved_for_the_stated_reason(arm):
    """Posting results is not posting an answer: all three posted group means with no analysis."""
    posted = arm[arm["channel"] == "registry_posted_results"]
    assert len(posted) >= 1
    assert (posted["outcome"] == "UNRESOLVED").all()
    assert posted["caveat"].str.contains("statistical analyses", case=False).all()
