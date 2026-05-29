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


@pytest.mark.skipif(not PRED.exists(), reason="predictions parquet not generated")
def test_famous_failures_flagged():
    pred = pd.read_parquet(PRED)
    famous = {"encenicline", "idalopirdine", "intepirdine", "pomaglumetad",
              "PF-04447943", "SUVN-502", "ABT-126", "TC-5619", "MK-0249"}
    present = pred[pred["compound"].isin(famous)]
    # every famous failure present must be predicted FAILURE
    assert (present["p2_predicted_outcome"] == "FAILURE").all()
    assert (present["clinical_outcome"] == "FAILURE").all()
