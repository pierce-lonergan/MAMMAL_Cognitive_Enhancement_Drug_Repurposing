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


# --- the structural audit (added 2026-09-20) ---------------------------------------------------
# Found by re-verifying every NCT by hand: one row pointed at a study of abdominal massage for
# dysmenorrhea, two "pending" trials had completed in 2025, one row's primary was a
# negative-symptom scale, and both resolved predictions were dated after their own readouts.
# Accuracy showed none of it.

def test_prediction_dated_after_readout_is_flagged():
    from mammal_repurposing.reporting.prospective import audit_registry
    iss = audit_registry(_frame([dict(
        drug="x", nct="NCT1", status="RESOLVED", primary_endpoint="MCCB composite",
        readout_year=2024.0, prediction_date="2026-05-30", nct_verified="Y")]))
    assert any(i["code"] == "PREDICTION_AFTER_READOUT" for i in iss)


def test_prediction_dated_before_readout_is_clean():
    from mammal_repurposing.reporting.prospective import audit_registry
    iss = audit_registry(_frame([dict(
        drug="x", nct="NCT1", status="RESOLVED", primary_endpoint="MCCB composite",
        readout_year=2027.0, prediction_date="2026-05-30", nct_verified="Y")]))
    assert not any(i["code"] == "PREDICTION_AFTER_READOUT" for i in iss)


def test_non_cognition_primary_is_flagged_on_pending_rows_too():
    """A category error that has not been cashed out yet is still a category error."""
    from mammal_repurposing.reporting.prospective import audit_registry
    iss = audit_registry(_frame([dict(
        drug="x", nct="NCT1", status="PENDING", primary_endpoint="PANSS-negative",
        readout_year=float("nan"), prediction_date="2026-05-30", nct_verified="Y")]))
    assert any(i["code"] == "PRIMARY_NOT_COGNITION" for i in iss)


def test_row_that_looks_clean_but_registry_disagrees_is_flagged():
    """The worse case: the row's own label passes, the actual registered primary does not."""
    from mammal_repurposing.reporting.prospective import audit_registry
    iss = audit_registry(_frame([dict(
        drug="x", nct="NCT1", status="PENDING", primary_endpoint="cognition (exploratory)",
        ctgov_primary="adverse events; PK; PANSS total", readout_year=float("nan"),
        prediction_date="2026-05-30", nct_verified="Y")]))
    assert any(i["code"] == "PRIMARY_NOT_COGNITION_PER_REGISTRY" for i in iss)


def test_pending_row_whose_trial_already_completed_is_flagged():
    from mammal_repurposing.reporting.prospective import audit_registry
    iss = audit_registry(_frame([dict(
        drug="x", nct="NCT1", status="PENDING", primary_endpoint="NIH Toolbox cognition",
        ctgov_status="COMPLETED", ctgov_completion="2025-07-18", ctgov_primary="NIH Toolbox",
        readout_year=float("nan"), prediction_date="2026-05-30", nct_verified="Y")]))
    assert any(i["code"] == "STATUS_STALE" for i in iss)


def test_still_recruiting_pending_row_is_not_flagged_stale():
    from mammal_repurposing.reporting.prospective import audit_registry
    iss = audit_registry(_frame([dict(
        drug="x", nct="NCT1", status="PENDING", primary_endpoint="ADAS-Cog11",
        ctgov_status="RECRUITING", ctgov_completion="2028-09-11", ctgov_primary="ADAS-Cog11",
        readout_year=float("nan"), prediction_date="2026-05-30", nct_verified="Y")]))
    assert not any(i["code"] == "STATUS_STALE" for i in iss)
    assert iss == [], "a fully clean row should raise nothing at all"


def test_live_registry_has_exactly_one_structurally_clean_row():
    """Pins the state found on 2026-09-20: only MINDSET 2 survives the audit."""
    from pathlib import Path
    from mammal_repurposing.reporting.prospective import audit_registry, load_prospective
    root = Path(__file__).resolve().parents[1]
    df = load_prospective(root / "data" / "raw" / "prospective_predictions.csv")
    iss = audit_registry(df)
    clean = set(df["nct"].astype(str)) - {str(i["nct"]) for i in iss}
    assert clean == {"NCT06976203"}, f"expected only MINDSET 2 clean, got {clean}"


def test_no_row_still_carries_the_bad_dysmenorrhea_nct():
    from pathlib import Path
    from mammal_repurposing.reporting.prospective import load_prospective
    root = Path(__file__).resolve().parents[1]
    df = load_prospective(root / "data" / "raw" / "prospective_predictions.csv")
    assert "NCT03821207" not in set(df["nct"].astype(str)), (
        "NCT03821207 resolves to an abdominal-massage dysmenorrhea study, not luvadaxistat")
    assert "NCT03382639" in set(df["nct"].astype(str))


# --- the pre-registered healthy-adult forward rule (2026-09-20) --------------------------------

def _led(rows):
    import pandas as pd
    return pd.DataFrame(rows)


def test_rule_calls_positive_only_when_interval_clears_zero_and_the_floor():
    from mammal_repurposing.reporting.prospective import predict_healthy_adult
    led = _led([dict(compound="x", representative_g=0.28, ci_lo=0.21, ci_hi=0.36,
                     citation_short="X 2020")])
    assert predict_healthy_adult("x", led)["prediction"] == "POSITIVE"


def test_real_but_trivial_effect_is_called_null_not_positive():
    """An interval above zero with g under the floor must NOT become a POSITIVE call."""
    from mammal_repurposing.reporting.prospective import predict_healthy_adult
    led = _led([dict(compound="x", representative_g=0.12, ci_lo=0.02, ci_hi=0.21,
                     citation_short="X 2020")])
    assert predict_healthy_adult("x", led)["prediction"] == "NULL"


def test_interval_spanning_zero_is_null_and_below_zero_is_negative():
    from mammal_repurposing.reporting.prospective import predict_healthy_adult
    led = _led([
        dict(compound="span", representative_g=0.05, ci_lo=-0.11, ci_hi=0.19, citation_short=""),
        dict(compound="neg", representative_g=-0.40, ci_lo=-0.70, ci_hi=-0.10, citation_short=""),
    ])
    assert predict_healthy_adult("span", led)["prediction"] == "NULL"
    assert predict_healthy_adult("neg", led)["prediction"] == "NEGATIVE"


def test_missing_compound_and_missing_interval_both_abstain():
    from mammal_repurposing.reporting.prospective import predict_healthy_adult
    led = _led([dict(compound="noci", representative_g=0.30, ci_lo=float("nan"),
                     ci_hi=float("nan"), citation_short="")])
    assert predict_healthy_adult("noci", led)["prediction"] == "ABSTAIN"
    assert predict_healthy_adult("absent", led)["prediction"] == "ABSTAIN"


def test_lookup_is_case_and_whitespace_insensitive():
    from mammal_repurposing.reporting.prospective import predict_healthy_adult
    led = _led([dict(compound="Caffeine", representative_g=0.28, ci_lo=0.21, ci_hi=0.36,
                     citation_short="")])
    assert predict_healthy_adult("  caffeine ", led)["prediction"] == "POSITIVE"


def test_every_prediction_carries_its_own_reason():
    from mammal_repurposing.reporting.prospective import predict_healthy_adult
    led = _led([dict(compound="x", representative_g=0.28, ci_lo=0.21, ci_hi=0.36,
                     citation_short="X 2020")])
    for c in ("x", "absent"):
        assert predict_healthy_adult(c, led)["reason"].strip(), "a bare call is not disputable"


def test_baseline_is_the_constant_null_predictor_not_a_coin_flip():
    from mammal_repurposing.reporting.prospective import healthy_adult_baseline
    led = _led([dict(compound=f"p{i}", ci_lo=0.1, ci_hi=0.3) for i in range(3)]
               + [dict(compound=f"n{i}", ci_lo=-0.2, ci_hi=0.2) for i in range(7)])
    b = healthy_adult_baseline(led)
    assert b["n_with_interval"] == 10 and b["n_positive"] == 3 and b["n_null"] == 7
    assert b["constant_call"] == "NULL"
    assert b["constant_accuracy"] == 0.7, "must reflect the real mix, not 0.5"


def test_live_ledger_baseline_is_a_strong_opponent():
    """Pins that the constant-NULL predictor is hard to beat, which is the whole point."""
    from pathlib import Path
    import pandas as pd
    from mammal_repurposing.reporting.prospective import healthy_adult_baseline
    root = Path(__file__).resolve().parents[1]
    led = pd.read_csv(root / "data" / "raw" / "healthy_adult_cognition_ledger.csv")
    b = healthy_adult_baseline(led)
    assert b["n_with_interval"] >= 30
    assert 0.2 < b["constant_accuracy"] < 0.8


# --- the healthy-adult forward registry (2026-09-20) -------------------------------------------

def _ha():
    from pathlib import Path
    import pandas as pd
    root = Path(__file__).resolve().parents[1]
    return pd.read_csv(root / "data" / "raw" / "prospective_healthy_adult.csv")


def _s135():
    import importlib.util
    from pathlib import Path
    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location(
        "s135", root / "scripts" / "135_prospective_healthy_adult.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_every_row_is_pending_and_unscored():
    """The registry is only worth anything if nothing in it has been filled in after the fact."""
    df = _ha()
    assert (df["status"] == "PENDING").all()
    assert df["actual_outcome"].fillna("").eq("").all()


def test_every_row_carries_a_reason_and_a_frozen_date():
    df = _ha()
    assert df["prediction_reason"].fillna("").str.strip().ne("").all()
    assert (df["prediction_date"] == "2026-09-20").all()
    assert (df["rule_committed_at"] == "2baedf8").all()


def test_microdosing_maps_to_the_microdosing_ledger_row():
    """The full-dose psilocybin estimate is -1.13; microdosing is -0.34. Conflating them was a bug."""
    m = _s135()
    assert m.COMPOUND_MAP["NCT07449351"] == "psilocybin_lsd_microdosing"


def test_non_cognitive_primaries_are_excluded_not_predicted():
    """Imaging, self-report and BOLD primaries cannot grade a cognition claim."""
    m, df = _s135(), _ha()
    for ident in ("NCT06367738", "NCT06768944", "NCT06041048"):
        assert ident in m.ENDPOINT_EXCLUDE
        row = df[df["identifier"] == ident].iloc[0]
        assert row["prediction"] == "ABSTAIN"
        assert row["prediction_basis"] == "endpoint_not_cognitive"


def test_post_acute_endpoints_are_excluded_from_an_acute_rule():
    """The ledger's g is on-drug; predicting a post-washout readout from it is a category error."""
    m, df = _s135(), _ha()
    assert "NCT06692192" in m.TIMING_EXCLUDE
    row = df[df["identifier"] == "NCT06692192"].iloc[0]
    assert row["prediction"] == "ABSTAIN"
    assert row["prediction_basis"] == "post_acute_endpoint"


def test_no_confirmed_registration_is_silently_dropped():
    """Every confirmed row must appear with either a call or a recorded abstention reason."""
    from pathlib import Path
    import pandas as pd
    root = Path(__file__).resolve().parents[1]
    pend = pd.read_csv(root / "data" / "raw" / "pending_readouts_2026-09.csv")
    conf = set(pend[pend["verdict"] == "CONFIRMED"]["identifier"].astype(str))
    assert conf == set(_ha()["identifier"].astype(str))


def test_secondary_endpoint_rows_are_tagged_for_stratification():
    df = _ha()
    tagged = df[df["endpoint_tier"] == "secondary"]
    assert len(tagged) >= 3
    assert (tagged["prediction"] != "ABSTAIN").any(), "tagging should not silently abstain them"


def test_abstentions_outnumber_calls_and_that_is_expected():
    """A registry that forces a call on every row is manufacturing a track record."""
    df = _ha()
    assert (df["prediction"] == "ABSTAIN").sum() > (df["prediction"] != "ABSTAIN").sum() / 2


def test_the_overturned_psychosis_trial_never_returns():
    """NCT07231497 is COSTEP 1, single-arm, unmasked, medicated psychosis. It passed the sweep's
    verifier because the registry's own healthyVolunteers flag contradicts its eligibility text."""
    m, df = _s135(), _ha()
    assert "NCT07231497" not in m.COMPOUND_MAP
    assert "NCT07231497" not in set(df["identifier"].astype(str))


def test_the_overturn_is_recorded_not_silently_dropped():
    from pathlib import Path
    import pandas as pd
    root = Path(__file__).resolve().parents[1]
    pend = pd.read_csv(root / "data" / "raw" / "pending_readouts_2026-09.csv")
    row = pend[pend["identifier"] == "NCT07231497"].iloc[0]
    assert row["verdict"] == "REJECT_WRONG_POPULATION"
    assert "OVERTURNED" in str(row["post_sweep_correction"])


def test_modafinil_has_no_live_prediction():
    """Both modafinil readouts fell: one overturned on population, one on a BOLD primary."""
    df = _ha()
    assert df[df["ledger_compound"] == "modafinil"].empty
