"""Tests for the F2 novel-compound onboarding engine.

RDKit + numpy/pandas. Unit tests for the assignment / prior / abstention logic on
synthetic data, plus a real-ledger regression that locks the F2 headline (~0.97
leave-one-compound-out class recovery on the expanded exemplar base).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("rdkit")

from mammal_repurposing.validation.novel_compound import (  # noqa: E402
    MIN_CLASS_N, TAU_OOD,
    _profile_stats, build_class_priors, build_exemplars, build_profile_centroids,
    is_allosteric_class, load_profiles, loco_class_recovery, profile_class_scores,
    score_catalogue, score_compound,
)

ROOT = Path(__file__).resolve().parents[1]


# --------------------------------------------------------------------------
# synthetic fixtures - four structurally separable classes
# --------------------------------------------------------------------------

@pytest.fixture
def synth():
    ledger = pd.DataFrame({
        "compound": ["a", "b", "c", "d", "e", "f", "g"],
        "mechanism_class": ["alc", "alc", "acid", "acid",
                            "alpha7_nAChR", "alpha7_nAChR", "single"],
        "clinical_g": [0.5, 0.4, 0.0, 0.0, 0.3, 0.2, 0.6],
        "label": [1, 1, 0, 0, 1, 1, 1],
    })
    smiles = pd.DataFrame({
        "compound": ["a", "b", "c", "d", "e", "f", "g"],
        "smiles": ["CCO", "CCCO", "CC(=O)O", "CCC(=O)O",
                   "C1CCNCC1", "C1CCN(C)CC1", "CCCCCCCCCC"],
    })
    return ledger, smiles


# --------------------------------------------------------------------------
# allosteric flag
# --------------------------------------------------------------------------

def test_is_allosteric_class():
    assert is_allosteric_class("alpha7_nAChR")
    assert is_allosteric_class("AMPA_PAM")
    assert is_allosteric_class("mGluR")
    assert not is_allosteric_class("AChE_inhibitor")
    assert not is_allosteric_class("catecholaminergic_ADHD")


# --------------------------------------------------------------------------
# class priors
# --------------------------------------------------------------------------

def test_build_class_priors(synth):
    ledger, _ = synth
    pr = build_class_priors(ledger, n_boot=200, seed=0)
    assert set(pr) == {"alc", "acid", "alpha7_nAChR", "single"}
    # CrI brackets the point estimate for a multi-member class
    assert pr["alc"].g_ci_lo <= pr["alc"].prior_g <= pr["alc"].g_ci_hi
    # success class vs failure class ordering
    assert pr["alc"].p_success > pr["acid"].p_success
    assert pr["alc"].prior_g > pr["acid"].prior_g
    # singleton -> degenerate (zero-width) CrI, flagged thin via n
    assert pr["single"].n == 1
    assert pr["single"].g_ci_lo == pr["single"].g_ci_hi
    # allosteric flag carried
    assert pr["alpha7_nAChR"].allosteric and not pr["alc"].allosteric


# --------------------------------------------------------------------------
# exemplar library
# --------------------------------------------------------------------------

def test_build_exemplars_groups_and_skips(synth):
    ledger, smiles = synth
    ex = build_exemplars(ledger, smiles)
    assert ex.keys() == {"alc", "acid", "alpha7_nAChR", "single"}
    assert len(ex["alc"]) == 2 and len(ex["single"]) == 1
    # a compound with no SMILES is skipped, not crashed
    led2 = pd.concat([ledger, pd.DataFrame([{
        "compound": "z", "mechanism_class": "alc", "clinical_g": 0.1, "label": 0}])],
        ignore_index=True)
    ex2 = build_exemplars(led2, smiles)  # 'z' has no SMILES row
    assert len(ex2["alc"]) == 2


# --------------------------------------------------------------------------
# scoring + abstention guardrails
# --------------------------------------------------------------------------

def test_clean_route_high(synth):
    ledger, smiles = synth
    ex, pr = build_exemplars(ledger, smiles), build_class_priors(ledger, n_boot=200)
    # query identical to an alcohol exemplar -> routes to 'alc', non-allosteric
    s = score_compound("q", "CCO", ex, pr)
    assert s.assigned_class == "alc"
    assert s.similarity == pytest.approx(1.0, abs=1e-6)
    assert s.tier == "HIGH"
    assert s.predicted_outcome == "SUCCESS"
    assert np.isfinite(s.prior_g)


def test_allosteric_downgrade(synth):
    ledger, smiles = synth
    ex, pr = build_exemplars(ledger, smiles), build_class_priors(ledger, n_boot=200)
    # query identical to an alpha7 (allosteric) exemplar -> capped at MED
    s = score_compound("q", "C1CCNCC1", ex, pr)
    assert s.assigned_class == "alpha7_nAChR"
    assert s.allosteric_flag is True
    assert s.tier == "MED"
    assert "allosteric" in s.reason.lower()


def test_thin_prior_low(synth):
    ledger, smiles = synth
    ex, pr = build_exemplars(ledger, smiles), build_class_priors(ledger, n_boot=200)
    # query identical to the singleton-class member -> LOW (n < MIN_CLASS_N)
    s = score_compound("q", "CCCCCCCCCC", ex, pr)
    assert s.assigned_class == "single"
    assert s.tier == "LOW"
    assert pr["single"].n < MIN_CLASS_N


def test_out_of_manifold_abstains(synth):
    ledger, smiles = synth
    ex, pr = build_exemplars(ledger, smiles), build_class_priors(ledger, n_boot=200)
    # force OOD by a near-1.0 floor: a non-exact alcohol (sim < 1.0) now abstains
    s = score_compound("q", "CCCCCO", ex, pr, tau_ood=0.999)
    assert s.tier == "ABSTAIN"
    assert "out-of-manifold" in s.reason

    # this query is genuinely far from every exemplar at the default floor
    s2 = score_compound("q", "O=S(=O)(O)O", ex, pr)  # sulfuric acid
    assert s2.tier == "ABSTAIN" and s2.assigned_class is None


def test_unparseable_and_empty(synth):
    ledger, smiles = synth
    ex, pr = build_exemplars(ledger, smiles), build_class_priors(ledger, n_boot=200)
    assert score_compound("q", "not_a_smiles!!", ex, pr).reason == "unparseable SMILES"
    assert score_compound("q", "CCO", {}, pr).reason == "no exemplar library"


def test_score_catalogue_sorts_routed_first(synth):
    ledger, smiles = synth
    ex, pr = build_exemplars(ledger, smiles), build_class_priors(ledger, n_boot=200)
    cat = pd.DataFrame({"id": ["x", "y"], "smiles": ["O=S(=O)(O)O", "CCO"]})
    out = score_catalogue(cat, ex, pr)
    assert list(out["query_id"])[0] == "y"          # routed before abstained
    assert out.iloc[-1]["tier"] == "ABSTAIN"


def test_loco_returns_sane_structure(synth):
    ledger, smiles = synth
    rec = loco_class_recovery(ledger, smiles)
    assert set(rec) >= {"n_evaluable", "n_routed", "top1_acc", "abstain_rate", "detail"}
    assert rec["n_evaluable"] >= 1
    assert 0.0 <= rec["abstain_rate"] <= 1.0


# --------------------------------------------------------------------------
# DTI-profile class signal (F2 spec signal a)
# --------------------------------------------------------------------------

def test_load_profiles_pivots_and_imputes():
    long = pd.DataFrame({
        "compound": ["a", "a", "a", "b", "b", "b"],
        "target_uniprot": ["T1", "T2", "T3", "T1", "T2", "T3"],
        "predicted_pkd": [9.0, 5.0, np.nan, 8.0, 5.0, 7.0],  # one NaN -> imputed
    })
    profiles, order = load_profiles(long)
    assert set(order) == {"T1", "T2", "T3"}
    assert set(profiles) == {"a", "b"}
    assert all(np.isfinite(v).all() for v in profiles.values())  # NaN imputed


def test_profile_centroids_discriminate_classes():
    # class A prefers T1, class B prefers T3; a query like A must score A > B
    profiles = {"a": np.array([9., 5., 5.]), "b": np.array([8., 5., 5.]),
                "c": np.array([5., 5., 9.]), "d": np.array([5., 5., 8.])}
    class_of = {"a": "A", "b": "A", "c": "B", "d": "B"}
    mu, sd = _profile_stats(profiles)
    cent = build_profile_centroids(profiles, class_of, mu, sd)
    assert set(cent) == {"A", "B"}
    scores = profile_class_scores(np.array([9., 5., 5.]), cent, mu, sd)
    assert scores["A"] > scores["B"]
    # a T3-preferring query routes to B
    scores_b = profile_class_scores(np.array([5., 5., 9.]), cent, mu, sd)
    assert scores_b["B"] > scores_b["A"]


def test_build_profile_centroids_excludes_holdout():
    profiles = {"a": np.array([9., 5.]), "b": np.array([8., 5.]), "c": np.array([5., 9.])}
    class_of = {"a": "A", "b": "A", "c": "B"}
    mu, sd = _profile_stats(profiles)
    cent = build_profile_centroids(profiles, class_of, mu, sd, exclude="a")
    # A's centroid now rests on 'b' only; 'a' did not leak in
    assert np.allclose(cent["A"], (profiles["b"] - mu) / sd)


# --------------------------------------------------------------------------
# real-ledger regression: the F2 headline
# --------------------------------------------------------------------------

_RAW = ROOT / "data" / "raw"
_LEDGERS = [_RAW / "clinical_outcomes_ledger.csv",
            _RAW / "clinical_outcomes_ledger_EXTENSION.csv",
            _RAW / "clinical_outcomes_ledger_CTGOV.csv",
            _RAW / "clinical_outcomes_ledger_RESEARCH.csv"]
_SMILES = _RAW / "ledger_compound_smiles.csv"
_HAVE = _SMILES.exists() and all(p.exists() for p in _LEDGERS)


@pytest.mark.skipif(not _HAVE, reason="ledger / SMILES CSVs not present")
def test_f2_headline_class_recovery():
    from mammal_repurposing.reporting.trial_watch import load_combined_ledger
    led = load_combined_ledger(_LEDGERS)
    smi = pd.read_csv(_SMILES)[["compound", "smiles"]]
    ex = build_exemplars(led, smi)
    # expanded exemplar base covers most classes
    assert len(ex) >= 20
    assert sum(1 for v in ex.values() if len(v) >= 2) >= 15
    rec = loco_class_recovery(led, smi)
    # structure recovers the true mechanism class for routed compounds
    assert rec["top1_acc"] >= 0.90
    # abstention is substantial (the guardrail) but not total
    assert 0.3 < rec["abstain_rate"] < 0.85
    assert rec["n_routed"] >= 20
