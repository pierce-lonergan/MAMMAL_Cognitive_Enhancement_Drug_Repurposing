"""Robustness state of the healthy-adult headline, AFTER the 2026-07 verified ledger expansion.

HISTORY -- THESE ASSERTIONS WERE REVERSED BY BETTER DATA, DELIBERATELY.
At n=11 (pre-expansion) this file locked three findings that made the headline look fragile:
  R1  n_studies (a pure power proxy) BEAT the stimulant gate, AUROC 0.88 vs 0.86  -> confound
  R2  the gate lost significance under the stated CI rule, p 0.046 -> 0.176       -> not robust
  R3  ZERO of 7 labelled nulls had actually been refuted                          -> no closure
Expanding the primary set to n=19 with independently citation-verified meta-analyses REFUTED R1 and
R2 and largely resolved R3. The power confound was an artifact of the tiny sample. The tests below
now lock the CURRENT state; the old expectations are recorded above so the reversal is auditable and
nobody re-derives the superseded claim.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "data" / "raw" / "healthy_adult_cognition_ledger.csv"
SCRIPT = ROOT / "scripts" / "121_healthy_adult_robustness.py"
pytest.importorskip("pandas")


def _mod():
    spec = importlib.util.spec_from_file_location("har", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _primary():
    import pandas as pd
    return _mod().clean_ma(pd.read_csv(LEDGER))


def test_primary_set_excludes_impairment_exposures():
    """The primary set must contain only genuine ENHANCER CANDIDATES. Pooling in impairment
    exposures (acute alcohol, dehydration, daytime melatonin, acute psilocybin) would let a
    classifier score well trivially on easy negatives and inflate every AUROC."""
    p = _primary()
    assert "candidate_enhancer" in p.columns
    assert (p["candidate_enhancer"] != 0).all()
    for impairer in ["alcohol_acute", "dehydration", "melatonin", "psilocybin"]:
        assert impairer not in set(p["compound"]), f"{impairer} must not enter the primary set"
    # The guard is on SILENT shrinkage. The set went 19 -> 18 on 2026-09-20 when interval recovery
    # found that omega_3's source pools cognitively normal older adults TOGETHER WITH MCI
    # participants, which fails this ledger's own clean-healthy rule. That row was never eligible,
    # so removing it is a correction rather than a loss, and it is recorded in
    # data/raw/provenance/interval_recovery_2026-09.json and in PREREG_DEVIATIONS_2026-06.md.
    # The guard is on SILENT change, not on change. Every move of this number is recorded in
    # docs/PREREG_DEVIATIONS_2026-06.md: 19 -> 18 when omega_3 was found to pool MCI participants,
    # then 18 -> 21 when three independently-verified healthy-adult rows were curated in.
    assert len(p) == 21, (
        f"primary set is {len(p)}, expected 21. Any further change must be justified in the "
        "deviations ledger before this number is edited.")
    assert "omega_3" not in set(p["compound"]), (
        "omega_3 pools MCI participants and cannot be tiered clean_MA")
    # non-independence: caffeine combinations must never be pooled alongside caffeine itself
    for combo in ["theanine_plus_caffeine", "caffeine_plus_taurine"]:
        assert combo not in set(p["compound"]), (
            f"{combo} contains caffeine, which is already a labelled enhancer here; pooling both "
            "counts one effect twice and inflates every AUROC")
    assert "caffeine" in set(p["compound"]), "the independent caffeine row itself must remain"


def test_power_confound_does_not_survive_the_larger_sample():
    """R1 REVERSED. At n=11 the power proxy beat the biology gate (0.88 vs 0.86). At n=19 it
    collapses to chance, so study volume does NOT explain the stimulant gate."""
    r1 = _mod().r1_power_confound(_primary())
    assert r1["au_k"] < r1["au_stim"], "power proxy must no longer beat the gate"
    assert r1["p_k"] > 0.05, "power proxy must no longer be a significant predictor"
    assert r1["p_stim"] < 0.05, "the stimulant gate remains nominally significant"


def test_the_headline_does_not_survive_its_own_sensitivity_analysis():
    """R2, REVERSED TWICE, and the second reversal is the one that matters.

    At n = 11 the stimulant gate did NOT survive re-labelling l-theanine per the ledger's stated CI
    rule (p = 0.176). At n = 19 it DID (p = 0.0456), and this test asserted that. On 2026-09-20 the
    set became n = 18, because interval recovery found that omega_3's source pools cognitively
    normal older adults together with MCI participants and it was never eligible for the clean-MA
    tier. Removing that single mis-tiered row moves the sensitivity p from 0.0456 to 0.0562, and the
    headline stops surviving.

    The conclusion is not that the gate is dead. It is that a headline whose survival turns on the
    tiering of one row out of nineteen is not robust, which is the same thing
    ledger_resolving_power_v1.md says from the other direction: at this n, a chance ranker scores
    AUROC 0.18 to 0.80.

    This test now pins the REVERSAL, so that a future change which restores significance has to be
    deliberate and has to say why.
    """
    _, sens = _mod().r2_label_rule_consistency(_primary())
    assert sens["p_shipped"] < 0.05, "as shipped, the gate is still nominally significant"
    assert sens["p_rule"] > 0.05, (
        f"the headline must be reported as NOT surviving its own stated-rule sensitivity "
        f"(p_rule = {sens['p_rule']:.4f}). If this starts passing again, the deviations ledger "
        "needs a new entry explaining what changed.")
    assert sens["au_rule"] < sens["au_shipped"], "and the AUROC still drops under the re-labelling"


def test_l_theanine_remains_the_unique_label_rule_conflict():
    rows, sens = _mod().r2_label_rule_consistency(_primary())
    assert sens["conflicts"] == ["l_theanine"]
    assert sum(1 for r in rows if r["conflict"]) == 1


def test_expansion_actually_refuted_compounds():
    """R3 RESOLVED. At n=11 not one null was refutable (0 of 7). The expansion brought real CIs, so
    several compounds are now genuinely ruled out as meaningful (g>=0.2) enhancers."""
    r3 = _mod().r3_evidence_of_absence(_primary())
    verdicts = [r["verdict"] for r in r3]
    assert verdicts.count("REFUTED") >= 5, "the expansion must close real doors"
    refuted = {r["compound"] for r in r3 if r["verdict"] == "REFUTED"}
    # microdosing is the headline closure: heavily marketed as an enhancer, pooled effect negative
    assert "psilocybin_lsd_microdosing" in refuted
    for r in r3:
        if r["verdict"] == "REFUTED":
            assert r["ci_hi"] < 0.20


def test_enhancer_set_is_no_longer_exclusively_stimulant():
    """The pre-expansion claim "the clean enhancers are exactly the four acute CNS stimulants" is
    now FALSIFIED BY THE DATA (not merely by a sensitivity analysis): a non-stimulant clears the
    stated CI-excludes-zero bar."""
    p = _primary()
    enh = p[p["enhances_healthy_young"] == 1]
    assert len(enh) >= 5
    assert (enh["supergroup"] != "stimulant").any(), "expected at least one non-stimulant enhancer"
