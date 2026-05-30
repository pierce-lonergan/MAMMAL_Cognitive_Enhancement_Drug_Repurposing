"""Gap 3 — Tests for the leakage-audited retrospective clinical validation.

Locks: AUROC correctness, leave-one-compound-out self-exclusion (no leakage),
leave-one-class-out class removal, metric sanity, and the real-data headline
contrast (class track-record discriminates; target affinity does not).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mammal_repurposing.validation import retrospective as R

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "data" / "raw" / "clinical_outcomes_ledger.csv"
PRED = ROOT / "data" / "results" / "v2" / "retrospective_validation_predictions.parquet"


# ---------------------------------------------------------------------------
# AUROC correctness (numpy-only implementation)
# ---------------------------------------------------------------------------

def test_auroc_perfect_separation():
    scores = np.array([0.1, 0.2, 0.3, 0.8, 0.9, 1.0])
    labels = np.array([0, 0, 0, 1, 1, 1])
    assert R.auroc(scores, labels) == 1.0


def test_auroc_perfect_inversion():
    scores = np.array([0.9, 0.8, 0.7, 0.2, 0.1, 0.0])
    labels = np.array([0, 0, 0, 1, 1, 1])
    assert R.auroc(scores, labels) == 0.0


def test_auroc_ties_handled():
    scores = np.array([0.5, 0.5, 0.5, 0.5])
    labels = np.array([0, 1, 0, 1])
    # all tied → AUROC 0.5
    assert abs(R.auroc(scores, labels) - 0.5) < 1e-9


def test_auroc_matches_known_value():
    # classic small example
    scores = np.array([0.1, 0.4, 0.35, 0.8])
    labels = np.array([0, 0, 1, 1])
    # pairs: (s=0.35 vs 0.1)>, (0.35 vs 0.4)<, (0.8 vs 0.1)>, (0.8 vs 0.4)> → 3/4
    assert abs(R.auroc(scores, labels) - 0.75) < 1e-9


# ---------------------------------------------------------------------------
# Leakage discipline
# ---------------------------------------------------------------------------

@pytest.fixture
def mini_ledger():
    return pd.DataFrame({
        "compound": ["a1", "a2", "a3", "b1", "b2"],
        "mechanism_class": ["A", "A", "A", "B", "B"],
        "target_uniprot": ["P1", "P1", "P1", "P2", "P2"],
        "indication": ["x"] * 5,
        "clinical_outcome": ["SUCCESS", "SUCCESS", "SUCCESS", "FAILURE", "FAILURE"],
        "clinical_g": [0.4, 0.4, 0.4, 0.0, 0.0],
        "label": [1, 1, 1, 0, 0],
        "compound_lower": ["a1", "a2", "a3", "b1", "b2"],
    })


def test_class_loco_excludes_self(mini_ledger):
    pred = R.class_loco_g(mini_ledger, shrinkage_k0=0.0)  # no shrink → pure sibling mean
    # a1 predicted from a2,a3 only (both 0.4) → 0.4, NOT influenced by its own value
    assert abs(pred["a1"] - 0.4) < 1e-9
    # b1 predicted from b2 only (0.0) → 0.0
    assert abs(pred["b1"] - 0.0) < 1e-9


def test_class_loco_shrinks_toward_global(mini_ledger):
    pred = R.class_loco_g(mini_ledger, shrinkage_k0=1.0)
    gmean = mini_ledger["clinical_g"].mean()  # 0.24
    # a1: (2*0.4 + 1*0.24)/3 = 0.3467
    assert abs(pred["a1"] - (2 * 0.4 + gmean) / 3) < 1e-6
    # all predictions strictly between sibling mean and global mean
    assert 0.24 < pred["a1"] < 0.4


def test_leave_one_class_out_removes_whole_class(mini_ledger):
    pred = R.leave_one_class_out_g(mini_ledger)
    # a1 predicted from class B only (all 0.0) → 0.0
    assert abs(pred["a1"] - 0.0) < 1e-9
    # b1 predicted from class A only (all 0.4) → 0.4
    assert abs(pred["b1"] - 0.4) < 1e-9


def test_class_loco_discriminates_mini(mini_ledger):
    pred = R.class_loco_g(mini_ledger)
    s = np.array([pred[c] for c in mini_ledger["compound"]])
    y = mini_ledger["label"].to_numpy()
    assert R.auroc(s, y) == 1.0  # class-homogeneous → perfect


def test_leave_one_class_out_inverts_mini(mini_ledger):
    pred = R.leave_one_class_out_g(mini_ledger)
    s = np.array([pred[c] for c in mini_ledger["compound"]])
    y = mini_ledger["label"].to_numpy()
    # successes get class-B (low) prediction, failures get class-A (high) → inverted
    assert R.auroc(s, y) == 0.0


# ---------------------------------------------------------------------------
# Metric helpers
# ---------------------------------------------------------------------------

def test_permutation_p_significant_for_perfect():
    scores = np.array([0.1, 0.2, 0.3, 0.8, 0.9, 1.0])
    labels = np.array([0, 0, 0, 1, 1, 1])
    p = R.permutation_p(scores, labels, n_perm=2000, seed=1)
    assert p < 0.10


def test_failure_recall(mini_ledger):
    pred = R.class_loco_g(mini_ledger)
    recall, flagged = R.failure_recall(pred, mini_ledger, threshold=0.2)
    assert recall == 1.0  # both failures predicted < 0.2
    assert set(flagged) == {"b1", "b2"}


# ---------------------------------------------------------------------------
# Real-data integration (if artifacts present)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not LEDGER.exists(), reason="clinical ledger not built")
def test_real_ledger_loads_balanced():
    led = R.load_clinical_ledger(LEDGER)
    assert len(led) >= 25
    assert led["label"].sum() >= 8           # successes
    assert (1 - led["label"]).sum() >= 8     # failures
    assert led["mechanism_class"].nunique() >= 8


@pytest.mark.skipif(not LEDGER.exists(), reason="clinical ledger not built")
def test_real_headline_contrast():
    """The core scientific claim: class track-record (P2) discriminates clinical
    outcome far better than target affinity / relevance (P1)."""
    led = R.load_clinical_ledger(LEDGER)
    p2 = R.class_loco_g(led)
    s2 = np.array([p2[c] for c in led["compound"]])
    y = led["label"].to_numpy()
    auroc_p2 = R.auroc(s2, y)
    assert auroc_p2 >= 0.85, f"class-LOCO AUROC {auroc_p2} should be high"
    # all 31 correctly classified at the global-mean threshold
    thr = float(led["clinical_g"].mean())
    correct = sum(((p2[c] >= thr) == bool(row.label))
                  for c, (_, row) in zip(led["compound"], led.iterrows()))
    assert correct >= int(0.9 * len(led))


# ---------------------------------------------------------------------------
# Gap 6 — external benchmark helpers
# ---------------------------------------------------------------------------

def test_paired_auroc_bootstrap_detects_winner():
    # a ranks perfectly, b is inverted -> delta strongly positive
    a = np.array([0.1, 0.2, 0.3, 0.8, 0.9, 1.0])
    b = np.array([0.9, 0.8, 0.7, 0.2, 0.1, 0.0])
    y = np.array([0, 0, 0, 1, 1, 1])
    d = R.paired_auroc_bootstrap(a, b, y, n_boot=1000, seed=1)
    assert d["delta"] == pytest.approx(1.0, abs=1e-9)   # 1.0 - 0.0
    assert d["p_a_gt_b"] > 0.95


def test_target_popularity_score_log_records():
    led = pd.DataFrame({
        "compound": ["d1", "d2"], "mechanism_class": ["A", "B"],
        "target_uniprot": ["P1", "P2"], "indication": ["x", "y"],
        "clinical_outcome": ["SUCCESS", "FAILURE"], "clinical_g": [0.4, 0.0],
        "label": [1, 0], "compound_lower": ["d1", "d2"],
    })
    chembl = pd.DataFrame({
        "target_uniprot": ["P1", "P1", "P2"],
        "n_records": [99.0, 0.0, 9.0],   # P1 -> 99, P2 -> 9
    })
    s = R.target_popularity_score(led, chembl)
    assert s["d1"] == pytest.approx(np.log10(100.0))
    assert s["d2"] == pytest.approx(np.log10(10.0))


@pytest.mark.skipif(not LEDGER.exists(), reason="ledger absent")
def test_class_beats_leakagefree_target_predictors():
    """Gap-6 headline: the class track record out-ranks the genuinely
    leakage-free target-centric paradigms (affinity + genetic relevance)."""
    led = R.load_clinical_ledger(LEDGER)
    p2 = R.class_loco_g(led)
    s2 = np.array([p2[c] for c in led["compound"]])
    auroc_class = R.auroc(s2, led["label"].to_numpy())
    v6b_path = ROOT / "data" / "results" / "v2" / "cluster_d_posterior_expanded_v2_mh8_ta99.parquet"
    if v6b_path.exists():
        v6b = pd.read_parquet(v6b_path)
        rel = R.target_relevance_score(led, v6b)
        rows = led[led["compound"].isin(rel)]
        au_rel = R.auroc(np.array([rel[c] for c in rows["compound"]]),
                         rows["label"].to_numpy())
        assert auroc_class > au_rel        # class >> genetics
    assert auroc_class >= 0.95


@pytest.mark.skipif(not PRED.exists(), reason="predictions parquet not generated")
def test_famous_failures_flagged():
    pred = pd.read_parquet(PRED)
    famous = {"encenicline", "idalopirdine", "intepirdine", "pomaglumetad",
              "PF-04447943", "SUVN-502", "ABT-126", "TC-5619", "MK-0249"}
    present = pred[pred["compound"].isin(famous)]
    # every famous failure present must be predicted FAILURE
    assert (present["p2_predicted_outcome"] == "FAILURE").all()
    assert (present["clinical_outcome"] == "FAILURE").all()


# ---------------------------------------------------------------------------
# Review-driven additions: class-level CI, network + structure comparators
# ---------------------------------------------------------------------------

def test_class_cluster_bootstrap_pure_classes_stays_perfect(mini_ledger):
    """With outcome-pure classes the class-level (cluster) bootstrap CI is
    degenerate at [1, 1] — it does not widen, because resampling classes cannot
    break a perfect separation when every class is uniformly one outcome."""
    res = R.class_cluster_bootstrap_auroc(mini_ledger, n_boot=300, seed=0)
    assert res["auroc"] == pytest.approx(1.0)
    assert res["ci_lo"] == pytest.approx(1.0)
    assert res["ci_hi"] == pytest.approx(1.0)
    assert res["n_classes"] == 2


def test_class_cluster_bootstrap_widens_with_mixed_class():
    """A genuinely mixed class injects between-class variance, so the class-level
    CI must widen below 1.0 — proving the estimator is not trivially pinned."""
    led = pd.DataFrame({
        "compound": ["a1", "a2", "b1", "b2", "c1", "c2"],
        "mechanism_class": ["A", "A", "B", "B", "C", "C"],
        "target_uniprot": ["P1"] * 6, "indication": ["x"] * 6,
        "clinical_outcome": ["SUCCESS", "SUCCESS", "FAILURE", "FAILURE",
                             "SUCCESS", "FAILURE"],   # class C is MIXED
        "clinical_g": [0.5, 0.5, 0.0, 0.0, 0.5, 0.0],
        "label": [1, 1, 0, 0, 1, 0],
        "compound_lower": ["a1", "a2", "b1", "b2", "c1", "c2"],
    })
    res = R.class_cluster_bootstrap_auroc(led, n_boot=800, seed=1)
    assert res["ci_lo"] < 1.0          # mixed class => genuine uncertainty


def test_kg_network_score_maps_ppr(mini_ledger):
    kg = pd.DataFrame({
        "compound_name": ["A1", "B1", "zzz"],   # case-insensitive join
        "kg_ppr_sum": [0.9, 0.1, 0.5],
    })
    led = mini_ledger.copy()
    led["compound"] = ["A1", "x2", "x3", "B1", "x5"]
    led["compound_lower"] = led["compound"].str.lower()
    s = R.kg_network_score(led, kg)
    assert s["A1"] == pytest.approx(0.9)
    assert s["B1"] == pytest.approx(0.1)
    assert "x2" not in s               # absent from kg table -> dropped


def test_structure_nn_success_score_prefers_success_neighbour():
    """The LOO structure score = max Tanimoto to another SUCCESS drug. A drug
    chemically identical to a known success should score ~1.0."""
    pytest.importorskip("rdkit")
    led = pd.DataFrame({
        "compound": ["aspirin", "aspirin_twin", "octane"],
        "mechanism_class": ["A", "A", "B"],
        "target_uniprot": ["P1", "P1", "P2"], "indication": ["x"] * 3,
        "clinical_outcome": ["SUCCESS", "FAILURE", "FAILURE"],
        "clinical_g": [0.5, 0.0, 0.0], "label": [1, 0, 0],
        "compound_lower": ["aspirin", "aspirin_twin", "octane"],
    })
    comp = pd.DataFrame({
        "name": ["aspirin", "aspirin_twin", "octane"],
        "smiles": ["CC(=O)Oc1ccccc1C(=O)O",      # aspirin
                   "CC(=O)Oc1ccccc1C(=O)O",      # identical twin (a FAILURE)
                   "CCCCCCCC"],                  # dissimilar
    })
    s = R.structure_nn_success_score(led, comp)
    # the failure twin is identical to the lone success -> ~1.0
    assert s["aspirin_twin"] == pytest.approx(1.0, abs=1e-6)
    # octane has no resemblance to the success -> low
    assert s["octane"] < 0.3


# ---------------------------------------------------------------------------
# Review-3 additions: temporal validation + taxonomy sensitivity
# ---------------------------------------------------------------------------

def _temporal_ledger():
    # class A succeeds (2000, 2005); class B fails (2001, 2006) — later members
    # are predictable from earlier same-class members
    return pd.DataFrame({
        "compound": ["a_old", "b_old", "a_new", "b_new"],
        "mechanism_class": ["A", "B", "A", "B"],
        "target_uniprot": ["P1", "P2", "P1", "P2"], "indication": ["x"] * 4,
        "clinical_outcome": ["SUCCESS", "FAILURE", "SUCCESS", "FAILURE"],
        "clinical_g": [0.5, 0.0, 0.5, 0.0], "label": [1, 0, 1, 0],
        "readout_year": [2000, 2001, 2005, 2006],
        "compound_lower": ["a_old", "b_old", "a_new", "b_new"],
    })


def test_temporal_holdout_predicts_later_from_earlier():
    led = _temporal_ledger()
    r = R.temporal_holdout_auroc(led, 2002)
    assert r["n_train"] == 2 and r["n_test"] == 2
    assert r["auroc"] == pytest.approx(1.0)   # earlier classes predict later members
    assert r["coverage"] == pytest.approx(1.0)


def test_prequential_excludes_first_of_class():
    led = _temporal_ledger()
    pq = R.prequential_class_loco(led)
    # only the 2nd member of each class has a strictly-earlier same-class sibling
    assert pq["n_informed"] == 2
    assert pq["auroc_informed"] == pytest.approx(1.0)
    tab = pq["table"]
    assert not bool(tab[tab["compound"] == "a_old"]["informed"].iloc[0])
    assert bool(tab[tab["compound"] == "a_new"]["informed"].iloc[0])


def test_auroc_under_coarse_taxonomy_can_drop(mini_ledger):
    # identity taxonomy reproduces the perfect separation
    ident = dict(zip(mini_ledger["mechanism_class"], mini_ledger["mechanism_class"]))
    assert R.auroc_under_taxonomy(mini_ledger, ident) == pytest.approx(1.0)
    # collapsing the two (pure, opposite-outcome) classes into one destroys it
    lumped = {"A": "all", "B": "all"}
    assert R.auroc_under_taxonomy(mini_ledger, lumped) < 1.0


def test_taxonomy_perturbation_observed_above_null(mini_ledger):
    res = R.taxonomy_perturbation_test(mini_ledger, n_perm=500, seed=0)
    assert res["observed"] == pytest.approx(1.0)
    assert res["null_mean"] < res["observed"]      # random grouping is worse
    assert 0.0 <= res["frac_reaching_observed"] <= 1.0


@pytest.mark.skipif(not LEDGER.exists(), reason="ledger absent")
def test_real_taxonomy_bracket():
    """Headline round-3: real mechanism-class taxonomy >> coarse >> random."""
    led = R.load_clinical_ledger(LEDGER)
    medium = R.auroc_under_taxonomy(led, dict(zip(led["mechanism_class"],
                                                  led["mechanism_class"])))
    coarse = R.auroc_under_taxonomy(led, R.COARSE_SYSTEM_MAP)
    pert = R.taxonomy_perturbation_test(led, n_perm=500, seed=0)
    assert medium == pytest.approx(1.0)
    assert coarse < 0.8                          # lumping mechanisms hurts
    assert pert["null_mean"] < 0.7               # random is near chance
    assert pert["frac_reaching_observed"] < 0.01  # ~0/N reach 1.0
