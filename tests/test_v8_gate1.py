"""Tests for the reusable V8 Gate 1 (src/.../cluster_e/gate1.py)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.cluster_e import gate1  # noqa: E402

pytestmark = pytest.mark.skipif(not gate1.SKLEARN_AVAILABLE,
                                reason="scikit-learn not installed")


def test_gate1_verdict_bands():
    assert gate1.gate1_verdict(0.60, 0.50) == "PASS"
    assert gate1.gate1_verdict(0.50, 0.40) == "PASS"
    assert gate1.gate1_verdict(0.40, 0.30) == "DEGRADE"
    assert gate1.gate1_verdict(0.30, 0.25) == "DEGRADE"
    assert gate1.gate1_verdict(0.29, 0.24) == "FAIL"
    assert gate1.gate1_verdict(0.10, 0.10) == "FAIL"
    # high AMI but low ARI does not pass
    assert gate1.gate1_verdict(0.55, 0.10) != "PASS"


def test_compound_consensus_averages_replicates():
    # two compounds, 3 + 2 replicate signatures
    X = np.array([[1.0, 1.0], [3.0, 3.0], [2.0, 2.0],   # compound A -> mean (2,2)
                  [10.0, 0.0], [0.0, 10.0]])            # compound B -> mean (5,5)
    keys = ["A", "A", "A", "B", "B"]
    Xc, uniq = gate1.compound_consensus(X, keys)
    assert uniq == ["A", "B"]
    assert np.allclose(Xc[0], [2.0, 2.0])
    assert np.allclose(Xc[1], [5.0, 5.0])


def test_cluster_and_score_recovers_separable_structure():
    rng = np.random.default_rng(0)
    # 3 well-separated classes in 5-d
    centroids = np.eye(3, 5) * 6.0
    X, y = [], []
    for c in range(3):
        for _ in range(12):
            X.append(centroids[c] + rng.normal(0, 0.3, 5))
            y.append(c)
    res = gate1.cluster_and_score(np.array(X), np.array(y),
                                  method="agglomerative", n_clusters=3)
    assert res["ami"] > 0.8
    assert gate1.gate1_verdict(res["ami"], res["ari"]) == "PASS"


def test_cluster_and_score_does_not_inflate_on_noise():
    rng = np.random.default_rng(1)
    # labels are random vs the data -> no recoverable structure
    X = rng.normal(0, 1, (36, 5))
    y = rng.integers(0, 3, 36)
    res = gate1.cluster_and_score(X, y, method="agglomerative", n_clusters=3)
    assert res["ami"] < 0.3
    assert gate1.gate1_verdict(res["ami"], res["ari"]) == "FAIL"
