"""Locks the robustness audit of the healthy-adult headline (scripts/121).

The axis report's headline — "the only separator is a coarse acute-CNS-stimulant gate, AUROC 0.86,
permutation p = 0.046" — is FRAGILE in three specific, reproducible ways. These tests pin each one
so the caveat cannot quietly disappear from the record.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
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


def test_power_proxy_matches_or_beats_the_biology_gate():
    """R1: n_studies carries ZERO biology, yet predicts the label at least as well as the stimulant
    gate — so the gate is not clean evidence of a pharmacological signal."""
    r1 = _mod().r1_power_confound(_primary())
    assert r1["au_k"] >= r1["au_stim"], (
        f"power proxy AUROC {r1['au_k']:.2f} should match/beat the stimulant gate "
        f"{r1['au_stim']:.2f}; if this flips, the confound argument must be revisited")
    # labelled enhancers are the more heavily studied compounds
    assert r1["k_med_enh"] > r1["k_med_null"]


def test_l_theanine_is_the_unique_label_rule_conflict():
    """R2: exactly one compound's assigned label disagrees with the ledger's own stated
    'CI excludes 0' rule, and it is the sole non-stimulant that would qualify."""
    rows, sens = _mod().r2_label_rule_consistency(_primary())
    assert sens["conflicts"] == ["l_theanine"]
    conflicted = [r for r in rows if r["conflict"]]
    assert len(conflicted) == 1 and conflicted[0]["supergroup"] != "stimulant"
    # it is a real positive point estimate, not a rounding artifact
    assert conflicted[0]["ci_lo"] > 0 and conflicted[0]["g"] > 0.3


def test_headline_significance_does_not_survive_the_stated_rule():
    """R2: the ONE significant healthy-adult result loses significance under a single defensible
    re-reading, and a non-stimulant then appears in the enhancer set."""
    _, sens = _mod().r2_label_rule_consistency(_primary())
    assert sens["p_shipped"] < 0.05, "as shipped the gate is nominally significant"
    assert sens["p_rule"] > 0.05, "under the stated CI rule it must NOT remain significant"
    assert sens["au_rule"] < sens["au_shipped"]
    assert "l_theanine" in sens["nonstim_enhancers_under_rule"]


def test_no_null_has_actually_been_refuted():
    """R3: a 'null' whose CI still admits g>=0.2 is under-powered, not refuted.

    The measured state is stronger than 'most': at the g>=0.2 threshold NOT ONE labelled null has
    been genuinely refuted — every one is either inconclusive or has no CI recorded at all. So the
    healthy-adult evidence base does not currently support ANY claim of the form 'compound X has
    been ruled out as a meaningful enhancer'."""
    r3 = _mod().r3_evidence_of_absence(_primary())
    verdicts = [r["verdict"] for r in r3]
    assert verdicts.count("REFUTED") == 0, (
        "if a null ever becomes genuinely refuted, this claim must be updated in the report too")
    assert verdicts.count("INCONCLUSIVE") >= 3
    assert verdicts.count("NO CI RECORDED") >= 3
    for r in r3:
        if r["verdict"] == "INCONCLUSIVE":
            assert np.isnan(r["ci_hi"]) or r["ci_hi"] >= 0.20
