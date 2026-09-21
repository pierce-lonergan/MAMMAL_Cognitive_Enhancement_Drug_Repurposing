"""Pin the three filters that decide what enters the retrospective arm.

All three were wrong on the first pass and each let through a different kind of junk:

  1. The cognitive-primary filter. ClinicalTrials.gov expands search terms, so a query for
     "attention" returns "Plasma Cotinine CONCENTRATION". Three of the first four candidates were
     pharmacokinetic studies that matched only through that expansion.
  2. The population filter. A registry healthyVolunteers flag is set true when a study enrols a
     healthy CONTROL arm beside patients, so two schizophrenia trials reached the usable list.
  3. The circularity guard's meta-analysis years, which were hardcoded from memory and wrong for
     five of fifteen compounds, including oxytocin (2013, written as 2025).
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def m():
    spec = importlib.util.spec_from_file_location(
        "s137", ROOT / "scripts" / "137_posted_results_arm.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# --- 1. the term-expansion trap ---------------------------------------------------------------

@pytest.mark.parametrize("text", [
    "Plasma Cotinine Concentration",
    "Area Under the Curve From Time 0 to the Last Quantifiable Sample, AUC(0-t)",
    "Measurement of Maximum Serum Concentration (Cmax)",
    "Differences in the Geometric Mean Ratio (GMR) of the Peak Concentration (Cmax)",
    "Urinary NNAL Concentration",
    "Number of Participants With Treatment-Emergent Adverse Events",
])
def test_pharmacokinetic_endpoints_are_rejected(m, text):
    """Every one of these was returned by a query for cognitive terms."""
    assert not m.is_cognitive_primary(text)


@pytest.mark.parametrize("text", [
    "Modified Suicide Stroop Task",
    "Reaction Time | Signal Detection Performance",
    "Change from baseline in CANTAB Paired Associates Learning",
    "NIH Toolbox Cognition Battery composite",
    "Psychomotor Vigilance Task mean reaction time",
    "Trail Making Test Part B time to completion",
])
def test_real_cognitive_endpoints_are_kept(m, text):
    assert m.is_cognitive_primary(text)


def test_a_cognitive_word_inside_a_pk_endpoint_still_loses(m):
    """Mentioning a battery does not rescue an exposure endpoint."""
    assert not m.is_cognitive_primary(
        "Plasma concentration of drug, with Stroop administered as a safety check")


# --- 2. inclusion versus exclusion ------------------------------------------------------------

def test_diagnosis_in_inclusion_criteria_disqualifies(m):
    assert m.patient_population(
        "The Effects of Nicotine on Cognition in Schizophrenia",
        "Inclusion Criteria: Patients: DSM IV diagnosis of schizophrenia. "
        "Exclusion Criteria: none") == "Schizophrenia"


def test_same_diagnosis_in_exclusion_criteria_does_not(m):
    """'No history of psychosis' is what a HEALTHY study looks like."""
    assert m.patient_population(
        "Nicotinic Modulation of the Default Network",
        "Inclusion Criteria: Age 21 through 50. "
        "Exclusion Criteria: history of psychotic disorder, dementia, bipolar disorder") == ""


def test_diagnosis_in_the_title_disqualifies_whatever_the_criteria_say(m):
    assert m.patient_population("Cognition in Mild Cognitive Impairment", "Inclusion: adults")


def test_inclusion_text_stops_at_the_exclusion_heading(m):
    t = m.inclusion_text("Inclusion Criteria: healthy adults EXCLUSION CRITERIA: schizophrenia")
    assert "schizophrenia" not in t.lower()
    assert "healthy adults" in t


def test_eligibility_with_no_exclusion_heading_is_searched_whole(m):
    assert m.patient_population("A study", "Participants with bipolar disorder aged 30-50")


# --- 3. the circularity guard -----------------------------------------------------------------

def test_meta_analysis_year_is_read_from_the_ledger_not_a_literal(m):
    """Hardcoding these got five of fifteen wrong. They must come from the ledger row."""
    led = pd.read_csv(ROOT / "data" / "raw" / "healthy_adult_cognition_ledger.csv")
    assert m.ledger_ma_year(led, "oxytocin") == 2013, "was hardcoded as 2025"
    assert m.ledger_ma_year(led, "dehydration") == 2018, "was hardcoded as 2025"
    assert m.ledger_ma_year(led, "nicotine") == 2010
    assert m.ledger_ma_year(led, "caffeine") == 2025
    assert m.ledger_ma_year(led, "not_a_compound") is None


def test_search_terms_cover_only_discordant_generating_compounds(m):
    """A compound whose rule call is NULL_EFFECT cannot produce an informative pair, so searching
    for it is wasted effort. None should appear in the target list."""
    from mammal_repurposing.reporting.prospective import predict_healthy_adult
    led = pd.read_csv(ROOT / "data" / "raw" / "healthy_adult_cognition_ledger.csv")
    for compound in m.SEARCH_TERMS:
        call = predict_healthy_adult(compound, led)["prediction"]
        assert call in ("POSITIVE", "NEGATIVE"), f"{compound} calls {call}"


# --- the live output --------------------------------------------------------------------------

def test_no_usable_row_is_a_patient_study_or_a_pk_study():
    d = pd.read_csv(ROOT / "data" / "raw" / "posted_results_candidates.csv")
    usable = d[d["usable"]]
    assert len(usable) > 0
    assert usable["patient_population"].fillna("").eq("").all()
    assert usable["cognitive_primary"].all()
    assert usable["independent_of_ledger"].all()


def test_the_two_schizophrenia_trials_never_return():
    """Both reached the usable list on the first pass with healthyVolunteers = true."""
    d = pd.read_csv(ROOT / "data" / "raw" / "posted_results_candidates.csv")
    for nct in ("NCT03838484", "NCT00383747"):
        row = d[d["nct"] == nct]
        assert len(row) == 1
        assert not bool(row.iloc[0]["usable"])
        assert str(row.iloc[0]["patient_population"]).lower().startswith("schizophren")
