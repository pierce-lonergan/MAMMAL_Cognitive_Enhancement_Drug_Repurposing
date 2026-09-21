"""Tests for the pre-registered prospective-prediction instrument."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from mammal_repurposing.reporting import prospective as P

ROOT = Path(__file__).resolve().parents[1]
CSV = ROOT / "data" / "raw" / "prospective_predictions.csv"


def _toy():
    return pd.DataFrame({
        "drug": ["a", "b", "c"],
        "mechanism_class": ["GlyT1", "PDE4", "M1_M4"],
        "indication": ["CIAS", "FXS", "AD"],
        "predicted_outcome": ["FAILURE", "SUCCESS", "SUCCESS"],
        "actual_outcome": ["FAILURE", None, None],
        "status": ["RESOLVED", "PENDING", "PENDING"],
    })


def test_score_resolved_counts_only_resolved():
    sc = P.score_resolved(_toy())
    assert sc["n_resolved"] == 1
    assert sc["n_correct"] == 1
    assert sc["accuracy"] == pytest.approx(1.0)


def test_score_resolved_detects_miss():
    df = _toy()
    df.loc[0, "actual_outcome"] = "SUCCESS"   # predicted FAILURE -> miss
    sc = P.score_resolved(df)
    assert sc["n_correct"] == 0
    assert sc["accuracy"] == pytest.approx(0.0)


def test_summary_partitions():
    s = P.summary(_toy())
    assert s["n_total"] == 3 and s["n_pending"] == 2 and s["n_resolved"] == 1


def test_cognition_primary_classifier():
    # cognitive batteries -> cognition-primary
    assert P.is_cognition_primary("Change From Baseline in the MCCB composite")
    assert P.is_cognition_primary("ADAS-Cog13 at week 26")
    assert P.is_cognition_primary("Severe Impairment Battery (SIB)")
    # safety / psychosis / behaviour primaries -> NOT cognition-primary
    assert not P.is_cognition_primary("Number of Treatment-Emergent Adverse Events")
    assert not P.is_cognition_primary("Change in PANSS Total Score")
    assert not P.is_cognition_primary("Aberrant Behavior Checklist")
    assert not P.is_cognition_primary("")
    assert not P.is_cognition_primary(None)


@pytest.mark.skipif(not CSV.exists(), reason="prospective CSV absent")
def test_real_prospective_predictions_resolved_correct():
    """The two resolved NMDA-coagonist-enhancer predictions (iclepertin GlyT1,
    luvadaxistat DAAO) were both FAILURE and both predicted FAILURE."""
    df = P.load_prospective(CSV)
    sc = P.score_resolved(df)
    assert sc["n_resolved"] >= 2
    assert sc["accuracy"] == pytest.approx(1.0)        # all resolved predicted right
    # there must be genuinely pending falsifiable predictions too
    assert P.summary(df)["n_pending"] >= 3


# --- the constant-baseline guard (added 2026-09-20) -------------------------------------------
# The registry reported "2/2 correct, accuracy 100%" as out-of-sample evidence when both resolved
# predictions were FAILURE, the majority outcome, which a constant predictor also gets right. These
# pin the correction so a headline accuracy can never again be quoted without its baseline.

def _frame(rows):
    import pandas as pd
    return pd.DataFrame(rows)


def test_agreeing_with_the_majority_is_not_evidence():
    from mammal_repurposing.reporting.prospective import score_resolved
    df = _frame([
        dict(drug="a", predicted_outcome="FAILURE", actual_outcome="FAILURE", status="RESOLVED"),
        dict(drug="b", predicted_outcome="FAILURE", actual_outcome="FAILURE", status="RESOLVED"),
    ])
    sc = score_resolved(df)
    assert sc["accuracy"] == 1.0
    assert sc["baseline_accuracy"] == 1.0
    assert sc["beats_baseline"] is False, "a constant predictor scores identically here"
    assert sc["n_informative"] == 0, "neither row departs from the majority"
    assert sc["p_value"] == 1.0


def test_departing_from_the_majority_and_being_right_IS_evidence():
    from mammal_repurposing.reporting.prospective import score_resolved
    rows = [dict(drug=f"f{i}", predicted_outcome="FAILURE", actual_outcome="FAILURE",
                 status="RESOLVED") for i in range(6)]
    rows += [dict(drug=f"s{i}", predicted_outcome="SUCCESS", actual_outcome="SUCCESS",
                  status="RESOLVED") for i in range(4)]
    sc = score_resolved(rows and _frame(rows))
    assert sc["majority_outcome"] == "FAILURE"
    assert sc["n_informative"] == 4
    assert sc["n_informative_right"] == 4
    assert sc["beats_baseline"] is True
    assert sc["p_value"] < 0.05


def test_being_wrong_on_the_departures_is_punished():
    from mammal_repurposing.reporting.prospective import score_resolved
    rows = [dict(drug=f"f{i}", predicted_outcome="FAILURE", actual_outcome="FAILURE",
                 status="RESOLVED") for i in range(6)]
    rows += [dict(drug=f"s{i}", predicted_outcome="SUCCESS", actual_outcome="FAILURE",
                  status="RESOLVED") for i in range(4)]
    sc = score_resolved(_frame(rows))
    assert sc["n_informative"] == 4
    assert sc["n_informative_right"] == 0
    assert sc["beats_baseline"] is False


def test_pending_informativeness_counts_real_bets():
    from mammal_repurposing.reporting.prospective import pending_informativeness
    df = _frame([
        dict(drug="a", predicted_outcome="FAILURE", actual_outcome="FAILURE", status="RESOLVED"),
        dict(drug="b", predicted_outcome="SUCCESS", actual_outcome=None, status="PENDING"),
        dict(drug="c", predicted_outcome="FAILURE", actual_outcome=None, status="PENDING"),
    ])
    pi = pending_informativeness(df)
    assert pi["assumed_majority"] == "FAILURE"
    assert pi["n_pending"] == 2
    assert pi["n_pending_informative"] == 1, "only the SUCCESS bet can be lost"
    assert pi["drugs"] == ["b"]


def test_live_registry_has_no_discriminative_evidence_yet():
    """Pins the actual state as of 2026-09-20, so a future change has to be deliberate."""
    from pathlib import Path
    from mammal_repurposing.reporting.prospective import (
        load_prospective, pending_informativeness, score_resolved,
    )
    root = Path(__file__).resolve().parents[1]
    df = load_prospective(root / "data" / "raw" / "prospective_predictions.csv")
    sc, pi = score_resolved(df), pending_informativeness(df)
    assert sc["n_informative"] == 0
    assert sc["beats_baseline"] is False
    assert pi["n_pending_informative"] == pi["n_pending"] > 0, (
        "every outstanding prediction should still be a real bet")


def test_no_resolved_rows_does_not_crash_or_claim_a_baseline():
    from mammal_repurposing.reporting.prospective import score_resolved
    sc = score_resolved(_frame([
        dict(drug="a", predicted_outcome="SUCCESS", actual_outcome=None, status="PENDING")]))
    assert sc["n_resolved"] == 0
    assert sc["n_informative"] == 0
    assert sc["beats_baseline"] is False
