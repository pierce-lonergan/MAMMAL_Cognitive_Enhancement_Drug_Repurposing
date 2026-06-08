"""Tests for the persistence evaluation metrics (over-claim, coverage-accuracy, budget).

Pure numpy - no engine - so these are fast and deterministic.
"""
from __future__ import annotations

from mammal_repurposing.validation.persistence_eval import (
    coverage_accuracy_curve, evaluate, label_budget, over_claims,
)


def test_over_claims_directional():
    # asserting MORE durability than the label supports is an over-claim
    assert over_claims("CANDIDATE_MECHANISTIC", "not_persistent") is True
    assert over_claims("DEMONSTRATED_HEALTHY", "not_persistent") is True
    assert over_claims("CONTESTED", "not_persistent") is True       # equivocal > negative
    # asserting equal-or-less durability is fine
    assert over_claims("NULL_SYMPTOMATIC", "not_persistent") is False
    assert over_claims("CONTESTED", "contested") is False
    assert over_claims("DISEASE_MODIFYING_PATIENTS", "disease_modifying_patients") is False
    assert over_claims("NULL_SYMPTOMATIC", "contested") is False    # under-claim allowed


def test_evaluate_counts_overclaims():
    recs = [
        {"compound": "a", "mechanism_class": "X", "persistence_verdict": "NULL_SYMPTOMATIC",
         "persistence_label": "not_persistent"},
        {"compound": "b", "mechanism_class": "X", "persistence_verdict": "CANDIDATE_MECHANISTIC",
         "persistence_label": "not_persistent"},   # over-claim
        {"compound": "c", "mechanism_class": "Y", "persistence_verdict": "CONTESTED",
         "persistence_label": "contested"},
    ]
    ev = evaluate(recs)
    assert ev["n"] == 3 and ev["n_over_claims"] == 1
    assert ev["over_claimers"] == ["b"]
    assert ev["per_mechanism"]["X"] == {"over": 1, "n": 2}


def test_coverage_accuracy_curve_monotone_coverage():
    recs = [
        {"persistence_verdict": "CONTESTED", "persistence_label": "contested", "_r": 7},
        {"persistence_verdict": "CANDIDATE_MECHANISTIC", "persistence_label": "not_persistent", "_r": 2},
        {"persistence_verdict": "NULL_SYMPTOMATIC", "persistence_label": "not_persistent", "_r": 0},
    ]
    curve = coverage_accuracy_curve(recs, lambda r: r["_r"])
    covs = [c["coverage"] for c in curve]
    assert covs == sorted(covs, reverse=True)   # coverage non-increasing as threshold rises
    assert all(0.0 <= c["coverage"] <= 1.0 for c in curve)


def test_label_budget_positive_and_scales_with_rarity():
    assert label_budget(prior=0.01) > label_budget(prior=0.1) > 0
