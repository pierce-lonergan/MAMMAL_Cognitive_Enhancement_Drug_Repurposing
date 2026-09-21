"""Tests for the prospective trial-watch engine (src/.../reporting/trial_watch.py).

Uses the real committed ledgers + prospective seed (small, deterministic).
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.reporting import trial_watch as tw  # noqa: E402

LEDGERS = [
    ROOT / "data" / "raw" / "clinical_outcomes_ledger.csv",
    ROOT / "data" / "raw" / "clinical_outcomes_ledger_EXTENSION.csv",
    ROOT / "data" / "raw" / "clinical_outcomes_ledger_CTGOV.csv",
]
PROSPECTIVE = ROOT / "data" / "raw" / "prospective_predictions.csv"


@pytest.fixture(scope="module")
def ledger():
    return tw.load_combined_ledger(LEDGERS)


@pytest.fixture(scope="module")
def prospective():
    return pd.read_csv(PROSPECTIVE, comment="#")


def test_combined_ledger_loads_and_dedups(ledger):
    assert len(ledger) == 47
    assert ledger["mechanism_class"].nunique() == 20
    assert set(ledger["label"].unique()) == {0, 1}
    # base success rate (15 successes / 47)
    assert abs(float(ledger["label"].mean()) - 15 / 47) < 1e-9


def test_norm_drug_strips_parenthetical():
    assert tw._norm_drug("xanomeline-trospium (KarXT)") == "xanomeline-trospium"
    assert tw._norm_drug("  Donepezil  ") == "donepezil"


def test_class_success_table_separates_success_from_failure(ledger):
    tbl = tw.class_success_table(ledger).set_index("mechanism_class")
    # outcome-pure success classes shrink high, failure classes shrink low
    assert tbl.loc["AChE_inhibitor", "p_shrunk"] > 0.7
    assert tbl.loc["catecholaminergic_ADHD", "p_shrunk"] > 0.8
    assert tbl.loc["alpha7_nAChR", "p_shrunk"] < 0.1
    assert tbl.loc["GlyT1_NMDA_coagonist", "p_shrunk"] < 0.15
    # raw rates are pure 0/1 (the famous outcome homogeneity)
    assert tbl.loc["AChE_inhibitor", "p_raw"] == 1.0
    assert tbl.loc["alpha7_nAChR", "p_raw"] == 0.0


def test_predict_high_confidence_class_failure(ledger):
    p = tw.predict_for("iclepertin", "GlyT1_NMDA_coagonist", "CIAS", ledger)
    assert p.predicted_outcome == "FAILURE"
    assert p.confidence == tw.CONF_HIGH
    assert p.evidence_level == "class"
    assert p.n_evidence == 2          # bitopertin + MK-5757 (iclepertin held out)
    assert p.p_success < 0.5


def test_predict_leaves_the_drug_itself_out(ledger):
    # bitopertin must be predicted from its siblings, never from its own row
    p = tw.predict_for("bitopertin", "GlyT1_NMDA_coagonist", "schizophrenia", ledger)
    assert p.predicted_outcome == "FAILURE"
    assert p.n_evidence == 2          # iclepertin + MK-5757
    assert "bitopertin" not in p.basis.lower()


def test_predict_superclass_carries_axis_evidence(ledger):
    # luvadaxistat (DAAO) has no direct class members; the NMDA-coagonist
    # super-class carries the GlyT1 failure evidence
    p = tw.predict_for("luvadaxistat", "DAAO_NMDA_coagonist", "CIAS", ledger)
    assert p.predicted_outcome == "FAILURE"
    assert p.evidence_level == "superclass"
    assert p.n_evidence >= 2


def test_predict_same_drug_new_indication_via_normalization(ledger):
    # KarXT alias must resolve to the ledger's xanomeline-trospium (same drug,
    # new indication), not be treated as a fresh super-class borrow
    p = tw.predict_for("xanomeline-trospium (KarXT)", "M1_M4_agonist",
                       "AD-cognition", ledger)
    assert p.evidence_level == "base_rate"
    assert "same-drug" in p.basis.lower()
    assert p.confidence == tw.CONF_LOW


def test_predict_abstains_on_unknown_class(ledger):
    p = tw.predict_for("madeupdrug", "TOTALLY_NOVEL_CLASS", "AD", ledger)
    assert p.confidence == tw.CONF_ABSTAIN
    assert p.n_evidence == 0
    assert abs(p.p_success - float(ledger["label"].mean())) < 1e-3


def test_roundtrip_class_purity_every_drug_predicts_itself(ledger):
    """Leave-one-drug-out over the whole ledger recovers every drug's own
    outcome, because the mechanism classes are outcome-pure. This is the
    class-prognostic result reproduced by the forward engine."""
    correct = 0
    for _, r in ledger.iterrows():
        p = tw.predict_for(r["compound"], r["mechanism_class"],
                           r.get("indication", ""), ledger)
        if p.predicted_outcome == r["clinical_outcome"]:
            correct += 1
    assert correct == len(ledger)     # 47/47


def test_build_registry_engine_matches_frozen(ledger, prospective):
    reg = tw.build_registry(prospective, ledger)
    assert len(reg) == len(prospective)
    # the automated engine reproduces every frozen hand prediction
    assert bool(reg["engine_matches_frozen"].all())


def test_score_registry_resolved_predictions(ledger, prospective):
    reg = tw.build_registry(prospective, ledger)
    sc = tw.score_registry(reg)
    assert sc["n_resolved"] == 2          # iclepertin + luvadaxistat
    assert sc["accuracy"] == 1.0          # both correctly called FAILURE
    assert sc["brier"] < 0.05             # well calibrated
    # AUROC undefined until both a success and a failure resolve
    assert sc["auroc"] != sc["auroc"]     # nan
    assert sc["n_failure"] == 2 and sc["n_success"] == 0


def test_score_registry_empty_is_graceful():
    empty = pd.DataFrame({"status": [], "actual_outcome": [], "p_success": [],
                          "predicted_outcome": [], "confidence": []})
    sc = tw.score_registry(empty)
    assert sc["n_resolved"] == 0


# --- the Brier baseline guard (added 2026-09-20) ----------------------------------------------
# The report benchmarked Brier against "0.25 = no-skill at base rate 0.5" while the ledger's own
# base rate is 0.319. A no-skill reference taken from the wrong base rate flatters or punishes the
# engine arbitrarily, so the comparator is now the ledger's real rate and these pin that.

def _reg(rows):
    import pandas as pd
    return pd.DataFrame(rows)


def test_omitting_the_base_rate_gives_nan_not_a_silent_default():
    from mammal_repurposing.reporting.trial_watch import score_registry
    sc = score_registry(_reg([
        dict(status="RESOLVED", actual_outcome="FAILURE", predicted_outcome="FAILURE",
             p_success=0.1, confidence="HIGH"),
    ]))
    assert sc["brier_base_rate"] != sc["brier_base_rate"], "expected NaN, not a 0.5 default"
    assert sc["beats_base_rate_brier"] is False, "cannot claim to beat an uncomputed baseline"


def test_engine_beating_the_constant_predictor_is_detected():
    from mammal_repurposing.reporting.trial_watch import score_registry
    rows = [dict(status="RESOLVED", actual_outcome="FAILURE", predicted_outcome="FAILURE",
                 p_success=0.05, confidence="HIGH") for _ in range(4)]
    sc = score_registry(_reg(rows), base_rate=0.319)
    assert sc["base_rate"] == 0.319
    assert sc["brier"] < sc["brier_base_rate"]
    assert sc["beats_base_rate_brier"] is True


def test_a_confidently_wrong_engine_loses_to_the_constant_predictor():
    from mammal_repurposing.reporting.trial_watch import score_registry
    rows = [dict(status="RESOLVED", actual_outcome="FAILURE", predicted_outcome="SUCCESS",
                 p_success=0.95, confidence="HIGH") for _ in range(4)]
    sc = score_registry(_reg(rows), base_rate=0.319)
    assert sc["brier"] > sc["brier_base_rate"]
    assert sc["beats_base_rate_brier"] is False


def test_base_rate_comparator_is_not_hardcoded_to_a_half():
    """A base rate far from 0.5 must move the comparator, or the fix did nothing."""
    from mammal_repurposing.reporting.trial_watch import score_registry
    rows = [dict(status="RESOLVED", actual_outcome="FAILURE", predicted_outcome="FAILURE",
                 p_success=0.2, confidence="HIGH") for _ in range(4)]
    lo = score_registry(_reg(rows), base_rate=0.10)["brier_base_rate"]
    hi = score_registry(_reg(rows), base_rate=0.90)["brier_base_rate"]
    assert lo < hi, "a lower base rate must score better against all-failure outcomes"
    assert abs(lo - 0.25) > 1e-9 and abs(hi - 0.25) > 1e-9
