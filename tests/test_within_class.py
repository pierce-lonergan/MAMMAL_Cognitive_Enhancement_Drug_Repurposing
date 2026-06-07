"""Tests for the F1 within-class resolution harness.

The statistical core (variance decomposition, pooled within-class Spearman,
within-class permutation, leave-one-compound-out MAE) is numpy/scipy only, so
these run in CI. A planted within-class signal must be detected; pure noise must
not be; and the real ledger must reproduce the headline "class is the resolution
limit" decomposition. RDKit-dependent feature builders are tested separately and
skip gracefully when RDKit is absent.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mammal_repurposing.validation.within_class import (
    variance_decomposition,
    within_class_spearman,
    loco_within_class_mae,
    _avg_ranks,
    rdkit_descriptors,
    class_centroid_tanimoto,
)

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "data" / "raw" / "clinical_outcomes_ledger.csv"


def _frame(classes, values, feature=None):
    rows = []
    for ci, (g_list) in enumerate(values):
        for j, g in enumerate(g_list):
            row = {"compound": f"c{ci}_{j}", "mechanism_class": classes[ci],
                   "clinical_g": g}
            if feature is not None:
                row["feat"] = feature[ci][j]
            rows.append(row)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# variance decomposition
# --------------------------------------------------------------------------

def test_variance_all_between():
    # each class constant within, classes differ -> all variance between
    df = _frame(["A", "B", "C"], [[0.4, 0.4], [0.2, 0.2], [0.0, 0.0]])
    vd = variance_decomposition(df)
    assert vd.frac_between == pytest.approx(1.0, abs=1e-9)
    assert vd.frac_within == pytest.approx(0.0, abs=1e-9)
    assert vd.icc1 == pytest.approx(1.0, abs=1e-6)


def test_variance_all_within():
    # one class, all spread -> nothing between
    df = _frame(["A"], [[0.0, 0.2, 0.4, 0.6]])
    vd = variance_decomposition(df)
    assert vd.frac_between == pytest.approx(0.0, abs=1e-9)
    assert vd.frac_within == pytest.approx(1.0, abs=1e-9)


def test_variance_mixed_monotone_in_separation():
    tight = _frame(["A", "B"], [[0.39, 0.41], [0.01, -0.01]])
    loose = _frame(["A", "B"], [[0.1, 0.7], [-0.3, 0.3]])
    assert variance_decomposition(tight).frac_between > variance_decomposition(loose).frac_between


# --------------------------------------------------------------------------
# pooled within-class Spearman + permutation
# --------------------------------------------------------------------------

def test_within_spearman_planted_signal_detected():
    # within every class, feature ranks value perfectly -> pooled rho = +1
    classes = ["A", "B", "C", "D"]
    vals = [[1.0, 2.0, 3.0, 4.0]] * 4
    feats = [[10.0, 20.0, 30.0, 40.0]] * 4   # monotone with value, different scale
    df = _frame(classes, vals, feats)
    r = within_class_spearman(df, "feat", n_perm=500, n_boot=200, seed=0)
    assert r.rho == pytest.approx(1.0, abs=1e-9)
    assert r.n_classes == 4
    assert r.perm_p < 0.05           # planted signal is significant
    assert r.ci_lo > 0.0             # CI excludes 0


def test_within_spearman_noise_not_detected():
    # half the classes align, half anti-align -> pooled rho cancels to ~0
    classes = ["A", "B", "C", "D"]
    vals = [[1.0, 2.0, 3.0, 4.0]] * 4
    feats = [[1.0, 2.0, 3.0, 4.0], [1.0, 2.0, 3.0, 4.0],
             [4.0, 3.0, 2.0, 1.0], [4.0, 3.0, 2.0, 1.0]]
    df = _frame(classes, vals, feats)
    r = within_class_spearman(df, "feat", n_perm=500, n_boot=200, seed=0)
    assert abs(r.rho) < 0.2
    assert r.perm_p > 0.2            # not significant


def test_within_spearman_skips_constant_and_singleton():
    # one singleton class + one flat class -> neither contributes
    df = _frame(["S", "F", "G"], [[0.3], [0.2, 0.2, 0.2], [1.0, 2.0]],
                [[9.0], [1.0, 2.0, 3.0], [5.0, 6.0]])
    r = within_class_spearman(df, "feat", n_perm=200, n_boot=100, seed=0)
    assert r.n_classes == 1          # only class G (non-constant, n>=2) is used


def test_permutation_p_is_valid_probability():
    df = _frame(["A", "B"], [[1.0, 2.0, 3.0], [1.0, 2.0, 3.0]],
                [[1.0, 2.0, 3.0], [3.0, 2.0, 1.0]])
    r = within_class_spearman(df, "feat", n_perm=300, n_boot=100, seed=1)
    assert 0.0 < r.perm_p <= 1.0


# --------------------------------------------------------------------------
# leave-one-compound-out MAE
# --------------------------------------------------------------------------

def test_loco_planted_feature_beats_class_mean():
    # within-class g is a clean linear function of feature -> augmenting helps
    classes = ["A", "B", "C"]
    base = {"A": 0.4, "B": 0.2, "C": 0.0}
    vals, feats = [], []
    for c in classes:
        f = [0.0, 1.0, 2.0, 3.0]
        vals.append([base[c] + 0.05 * x for x in f])   # g = base + 0.05*feat
        feats.append(f)
    df = _frame(classes, vals, feats)
    lo = loco_within_class_mae(df, "feat")
    assert lo.delta_mae > 0.0          # augmented beats the class mean
    assert lo.n_adjusted > 0


def test_loco_noise_does_not_beat_class_mean():
    rng = np.random.default_rng(0)
    classes = ["A", "B", "C"]
    base = {"A": 0.4, "B": 0.2, "C": 0.0}
    vals, feats = [], []
    for c in classes:
        vals.append([base[c]] * 4)                      # flat within class
        feats.append(list(rng.normal(size=4)))          # random feature
    df = _frame(classes, vals, feats)
    lo = loco_within_class_mae(df, "feat")
    assert lo.delta_mae <= 1e-9        # no improvement over the class mean


# --------------------------------------------------------------------------
# rank helper
# --------------------------------------------------------------------------

def test_avg_ranks_matches_scipy():
    sp = pytest.importorskip("scipy.stats")
    for a in ([3, 1, 2], [1, 1, 2, 2], [5.0], [2, 2, 2], [-1, 0, 1, 1, 1]):
        got = _avg_ranks(np.array(a, dtype=float))
        exp = sp.rankdata(a, method="average")
        assert np.allclose(got, exp)


# --------------------------------------------------------------------------
# real-ledger regression: the headline F1 result
# --------------------------------------------------------------------------

def test_real_ledger_class_is_resolution_limit():
    if not LEDGER.exists():
        pytest.skip("ledger CSV not present")
    led = pd.read_csv(LEDGER, comment="#")
    led = led[led["clinical_outcome"].isin(["SUCCESS", "FAILURE"])].copy()
    vd = variance_decomposition(led)
    # The headline F1 finding: the overwhelming majority of clinical-g variance
    # is BETWEEN mechanism classes; class membership all but determines g.
    assert vd.frac_between > 0.90
    assert vd.icc1 > 0.85
    assert vd.frac_within < 0.10


# --------------------------------------------------------------------------
# RDKit feature builders (skip if rdkit absent)
# --------------------------------------------------------------------------

def test_rdkit_descriptors_optional():
    pytest.importorskip("rdkit")
    d = rdkit_descriptors("CC(=O)Oc1ccccc1C(=O)O")   # aspirin
    assert d and d["mw"] == pytest.approx(180.16, abs=0.5)
    assert "cns_mpo" in d and 0.0 <= d["cns_mpo"] <= 5.0
    assert rdkit_descriptors("not_a_smiles") == {}


def test_class_centroid_tanimoto_optional():
    pytest.importorskip("rdkit")
    smi = {"a1": "c1ccccc1", "a2": "c1ccccc1C", "b1": "CCO"}
    cls = {"a1": "A", "a2": "A", "b1": "B"}
    out = class_centroid_tanimoto(smi, cls)
    assert 0.0 <= out["a1"] <= 1.0          # has a same-class peer (a2)
    assert np.isnan(out["b1"])              # singleton class -> NaN
