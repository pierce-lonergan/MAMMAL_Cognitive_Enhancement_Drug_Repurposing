"""Tests for the persistence evaluation metrics (over-claim, coverage-accuracy, budget).

Pure numpy - no engine - so these are fast and deterministic.
"""
from __future__ import annotations

from mammal_repurposing.validation.persistence_eval import (
    coverage_accuracy_curve, evaluate, label_budget, over_claims, sensitivity,
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


def test_sensitivity_recall_over_verified_positives():
    recs = [
        # flagged: PERSEUS asserts durability (level >= 1)
        {"compound": "a", "persistence_verdict": "CANDIDATE_MECHANISTIC", "domain": "neuroplasticity"},
        {"compound": "b", "persistence_verdict": "WINDOW_CONDITIONAL", "domain": "neuroplasticity"},
        # missed: PERSEUS dismisses as null / abstain
        {"compound": "c", "persistence_verdict": "NULL_SYMPTOMATIC", "domain": "cognition"},
        {"compound": "d", "persistence_verdict": "ABSTAIN", "domain": "mood"},
    ]
    s = sensitivity(recs)
    assert s["n"] == 4 and s["n_flagged"] == 2 and s["sensitivity"] == 0.5
    assert set(s["flagged"]) == {"a", "b"} and set(s["missed"]) == {"c", "d"}
    assert s["by_domain"]["neuroplasticity"] == {"flagged": 2, "n": 2}
    assert s["by_domain"]["cognition"] == {"flagged": 0, "n": 1}


def test_sensitivity_empty_is_nan():
    s = sensitivity([])
    assert s["n"] == 0 and s["sensitivity"] != s["sensitivity"]   # nan


# --- rigorous empty-positive PU evaluation (Gap 4) ---

def test_recall_ci_jeffreys_brackets_point_and_is_wide_at_small_n():
    from mammal_repurposing.validation.persistence_pu_eval import recall_ci
    r = recall_ci(7, 13)
    assert abs(r["recall"] - 7 / 13) < 1e-9
    assert r["lo"] < r["recall"] < r["hi"]
    assert (r["hi"] - r["lo"]) > 0.4          # n=13 -> honestly wide interval


def test_fpr_zero_has_nonzero_upper_bound():
    from mammal_repurposing.validation.persistence_pu_eval import fpr_ci
    f = fpr_ci(0, 15)
    assert f["fpr"] == 0.0 and f["lo"] == 0.0 and 0.1 < f["hi"] < 0.3   # 0/15 -> upper ~0.18


def test_ppv_curve_rises_with_prior_and_handles_zero_fpr():
    from mammal_repurposing.validation.persistence_pu_eval import ppv_curve
    c = ppv_curve(0.54, 0.05, priors=(0.01, 0.03))
    assert c[1]["ppv"] > c[0]["ppv"]          # higher prior -> higher PPV
    z = ppv_curve(0.54, 0.0, priors=(0.01,))
    assert z[0]["ppv"] == 1.0                 # zero FPR -> PPV 1 (degenerate, by construction)


def test_wilson_and_jeffreys_agree_roughly():
    from mammal_repurposing.validation.persistence_pu_eval import jeffreys_ci, wilson_ci
    wlo, whi = wilson_ci(7, 13)
    jlo, jhi = jeffreys_ci(7, 13)
    assert abs(wlo - jlo) < 0.12 and abs(whi - jhi) < 0.12


def test_grouped_lomo_per_mechanism_recall():
    from mammal_repurposing.validation.persistence_pu_eval import grouped_lomo
    recs = [
        {"compound": "psilocin", "mechanism_class": "serotonergic", "flagged": True},
        {"compound": "lsd", "mechanism_class": "serotonergic", "flagged": True},
        {"compound": "ketamine", "mechanism_class": "nmda", "flagged": False},
        {"compound": "scopolamine", "mechanism_class": "muscarinic", "flagged": False},
    ]
    g = grouped_lomo(recs)
    assert g["per_mechanism"]["serotonergic"]["recall"] == 1.0
    assert g["per_mechanism"]["nmda"]["recall"] == 0.0
    assert g["covered_mechanisms"] == ["serotonergic"] and g["fitted_model"] is False


def test_label_shift_transport_prior_correction():
    from mammal_repurposing.validation.persistence_pu_eval import label_shift_transport
    # zero FPR -> perfect precision; the expected confusion scales with the deployment prior
    t = label_shift_transport(sens=0.5, fpr=0.0, deploy_prior=0.01, n_screen=10000)
    assert t["tp"] == 50 and t["fp"] == 0 and t["ppv"] == 1.0
    # nonzero FPR at a 1% base rate sinks PPV even at decent sensitivity (the rare-disease trap)
    t2 = label_shift_transport(sens=0.5, fpr=0.15, deploy_prior=0.01, n_screen=10000)
    assert t2["fp"] > t2["tp"] and t2["ppv"] < 0.05
