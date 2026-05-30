"""Gap 4 — Tests for the allosteric learn-to-rank head.

Locks: physicochemical featurisation, feature-table assembly + imputation,
within-target Spearman, the quantified MAMMAL-flatness finding, and the
real-data headline (the fused head beats MAMMAL-alone on the held-out
allosteric benchmark).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("sklearn")
pytest.importorskip("rdkit")

from mammal_repurposing.cluster_a import allosteric_ltr as A

ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / "data" / "raw" / "allosteric_benchmark.csv"
CHEMBL = ROOT / "data" / "results" / "chembl_evidence.parquet"
DTI = ROOT / "data" / "results" / "dti_scores.parquet"
TANI = ROOT / "data" / "results" / "v2" / "disagreement_signal.parquet"


# ---------------------------------------------------------------------------
# Featurisation
# ---------------------------------------------------------------------------

def test_physchem_parses_known_drug():
    d = A.physchem_descriptors("O=C(OC1CC2CCC(C1)[NH+]2C)c1ccccc1")  # ~atropine-like
    assert np.isfinite(d["mw"]) and d["mw"] > 100
    assert set(A.PHYSCHEM_COLS) <= set(d)


def test_physchem_bad_smiles_nan():
    d = A.physchem_descriptors("not_a_smiles((")
    assert all(np.isnan(d[c]) for c in A.PHYSCHEM_COLS)


def test_build_feature_table_imputes():
    pairs = pd.DataFrame({
        "compound_name": ["a", "b"],
        "target_uniprot": ["P1", "P1"],
        "smiles": ["CCO", "c1ccccc1"],
    })
    mam = pd.DataFrame({"compound_name": ["a"], "target_uniprot": ["P1"],
                        "predicted_pkd": [6.0]})  # b missing -> imputed
    ft = A.build_feature_table(pairs, mammal=mam)
    assert set(A.FUSION_FEATURES) <= set(ft.columns)
    assert ft["mammal_pkd"].notna().all()        # imputed
    assert ft.loc[ft.compound_name == "a", "mammal_pkd"].iloc[0] == 6.0


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

def test_within_target_spearman_perfect():
    df = pd.DataFrame({
        "target_uniprot": ["P1"] * 4,
        "score": [1.0, 2.0, 3.0, 4.0],
        "label": [10.0, 20.0, 30.0, 40.0],
    })
    pooled, per = A.within_target_spearman(df, "score", "label")
    assert abs(pooled - 1.0) < 1e-9
    assert abs(per["P1"] - 1.0) < 1e-9


def test_within_target_spearman_skips_small_groups():
    df = pd.DataFrame({
        "target_uniprot": ["P1", "P1", "P2", "P2"],   # both groups n=2 < min_n
        "score": [1.0, 2.0, 1.0, 2.0],
        "label": [1.0, 2.0, 2.0, 1.0],
    })
    pooled, per = A.within_target_spearman(df, "score", "label", min_n=3)
    assert per == {} and np.isnan(pooled)


def test_mammal_flatness_detects_flat():
    df = pd.DataFrame({
        "target_uniprot": ["P1"] * 5,
        "mammal_pkd": [6.60, 6.61, 6.60, 6.62, 6.61],   # ~flat
    })
    fl = A.mammal_flatness(df)
    assert fl["P1"] < 0.05


# ---------------------------------------------------------------------------
# Real-data headline
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not (BENCH.exists() and CHEMBL.exists() and DTI.exists()),
                    reason="benchmark / chembl / dti artifacts absent")
def test_fused_beats_mammal_on_benchmark():
    dti = pd.read_parquet(DTI)
    tani = pd.read_parquet(TANI) if TANI.exists() else None
    bench = pd.read_csv(BENCH)
    bench["pact"] = 9.0 - np.log10(bench["measured_activity_nm"].astype(float))
    bf = A.build_feature_table(bench[["compound_name", "target_uniprot", "smiles"]],
                               mammal=dti, tanimoto=tani)
    bf = bf.merge(bench[["compound_name", "target_uniprot", "pact"]],
                  on=["compound_name", "target_uniprot"], how="left")

    ch = pd.read_parquet(CHEMBL)
    ch = ch[ch["best_pchembl"].notna()].copy()
    ch["pact"] = ch["best_pchembl"].astype(float)
    ch = ch[~ch["compound_name"].str.lower().isin(set(bench["compound_name"].str.lower()))]
    tf = A.build_feature_table(ch[["compound_name", "target_uniprot", "smiles"]],
                               mammal=dti, tanimoto=tani)
    tf = tf.merge(ch[["compound_name", "target_uniprot", "pact"]],
                  on=["compound_name", "target_uniprot"], how="left")
    tf = tf[tf["pact"].notna()].reset_index(drop=True)

    res = A.evaluate(tf, bf, label_col="pact", seed=0)
    # MAMMAL is flat within target (structural blindness)
    assert np.nanmean(list(res.flatness.values())) < 0.1
    # MAMMAL-alone ranking is ~chance; the fusion recovers a real ranking
    assert res.pooled_rho["mammal_only"] < 0.15
    assert res.pooled_rho["fused_ltr"] > 0.30
    assert res.pooled_rho["fused_ltr"] > res.pooled_rho["mammal_only"]


@pytest.mark.skipif(not (CHEMBL.exists() and DTI.exists()),
                    reason="chembl / dti artifacts absent")
def test_loto_scale_fused_beats_mammal():
    """Publication-strength scale-up: leave-one-target-out CV on the 297 real
    ChEMBL pChEMBL pairs (21 targets). The fused head must beat MAMMAL-alone
    within-target ranking across independent held-out targets."""
    dti = pd.read_parquet(DTI)
    tani = pd.read_parquet(TANI) if TANI.exists() else None
    ch = pd.read_parquet(CHEMBL)
    ch = ch[ch["best_pchembl"].notna()].copy()
    ch["pact"] = ch["best_pchembl"].astype(float)
    feat = A.build_feature_table(ch[["compound_name", "target_uniprot", "smiles"]],
                                 mammal=dti, tanimoto=tani)
    feat = feat.merge(ch[["compound_name", "target_uniprot", "pact"]],
                      on=["compound_name", "target_uniprot"], how="left")
    feat = feat[feat["pact"].notna()].reset_index(drop=True)
    res = A.loto_evaluate(feat, label_col="pact", seed=0)
    assert res.n_eval >= 250
    assert int(res.feature_importance["n_folds"]) >= 15
    assert res.pooled_rho["mammal_only"] < 0.2          # MAMMAL ~chance at scale
    assert res.pooled_rho["fused_ltr"] > 0.45           # fusion strong
    assert res.pooled_rho["fused_ltr"] > res.pooled_rho["mammal_only"]
