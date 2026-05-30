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
