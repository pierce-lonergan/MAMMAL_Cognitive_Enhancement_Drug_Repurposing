"""The healthy-adult ledger must satisfy its own stated contract (validation/ledger_guard.py)."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from mammal_repurposing.validation.ledger_guard import validate_ledger

LEDGER = Path(__file__).resolve().parents[1] / "data" / "raw" / "healthy_adult_cognition_ledger.csv"


def _base_row(**kw):
    row = {"compound": "x", "mechanism_class": "m", "supergroup": "nonstimulant",
           "primary_domain": "memory", "representative_g": 0.3, "ci_lo": 0.1, "ci_hi": 0.5,
           "n_studies": 5, "population": "healthy_young", "enhances_healthy_young": 1,
           "evidence_tier": "clean_MA", "citation_short": "A 2024 J", "pmid_doi": "PMID:1",
           "robustness": "", "verified": "yes"}
    row.update(kw)
    return pd.DataFrame([row])


def test_clean_row_passes():
    assert validate_ledger(_base_row()) == []


def test_label_rule_conflict_is_flagged():
    """The l-theanine class of defect: CI excludes 0 but the row is labelled null."""
    v = validate_ledger(_base_row(enhances_healthy_young=0))
    assert any(x.rule == "label_rule_conflict" for x in v)


def test_unjustified_conflict_is_an_error_but_justified_is_a_warning():
    unjust = validate_ledger(_base_row(enhances_healthy_young=0, robustness=""))
    assert any(x.rule == "label_rule_conflict" and x.severity == "error" for x in unjust)
    just = validate_ledger(_base_row(
        enhances_healthy_young=0,
        robustness="only choice-RT significant (faster); simple-RT null - borderline call"))
    conf = [x for x in just if x.rule == "label_rule_conflict"]
    assert conf and conf[0].severity == "warn"


def test_evidence_row_requires_provenance_and_effect():
    assert any(x.rule == "missing_provenance" for x in validate_ledger(_base_row(pmid_doi="")))
    assert any(x.rule == "missing_effect_size"
               for x in validate_ledger(_base_row(representative_g=None)))


def test_missing_ci_is_warned_not_fatal():
    v = validate_ledger(_base_row(ci_lo=None, ci_hi=None))
    assert any(x.rule == "missing_ci" and x.severity == "warn" for x in v)
    assert not any(x.severity == "error" for x in v)


def test_absent_tier_may_not_carry_an_effect_or_label():
    v = validate_ledger(_base_row(evidence_tier="absent"))
    rules = {x.rule for x in v}
    assert "absent_row_has_effect" in rules and "absent_row_has_label" in rules


def test_patient_population_is_an_error_but_healthy_variants_are_fine():
    """The guard must catch a CLINICAL sample without crying wolf on the ledger's own descriptive
    vocabulary for healthy samples (a whitelist previously false-alarmed on 7 legitimate rows)."""
    for ok in ["healthy_nonSD", "rested_healthy", "nonsmokers", "healthy_unimpaired",
               "mixed_stressed", "healthy_young"]:
        assert not any(x.rule == "patient_population"
                       for x in validate_ledger(_base_row(population=ok))), ok
    for bad in ["alzheimer_patients", "MCI", "adhd_adults", "schizophrenia"]:
        v = validate_ledger(_base_row(population=bad))
        assert any(x.rule == "patient_population" and x.severity == "error" for x in v), bad


def test_duplicate_compounds_are_an_error():
    df = pd.concat([_base_row(), _base_row()], ignore_index=True)
    assert any(x.rule == "duplicate_compound" for x in validate_ledger(df))


@pytest.mark.skipif(not LEDGER.exists(), reason="ledger absent")
def test_live_ledger_has_no_error_severity_violations():
    """The shipped ledger must be free of contract BREACHES. Warnings (missing CIs, the one
    justified label departure) are expected and are surfaced by scripts/121, not suppressed here."""
    v = validate_ledger(pd.read_csv(LEDGER))
    errors = [x for x in v if x.severity == "error"]
    assert not errors, "live ledger has error-severity violations:\n" + "\n".join(str(e) for e in errors)
