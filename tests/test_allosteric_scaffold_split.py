"""Gap 4 robustness — tests for the scaffold-split evaluation.

A split that silently leaks is the exact failure this whole analysis exists to
measure, so the leak-proofness of the split is pinned here rather than left to
inspection. Locks:

  * scaffold keying, including BOTH sentinel cases (unparseable, acyclic) and the
    fact that sentinels are UNIQUE per row, which is what stops them acting as a
    shared pseudo-scaffold;
  * train/test disjointness on scaffold key AND on InChIKey, the property whose
    absence is the leak;
  * that the scaffold filter actually bites on the real data (a split that
    removes nothing would pass a naive disjointness test trivially);
  * deterministic, RNG-free blocked-CV assignment with no scaffold crossing folds;
  * bootstrap determinism under the pre-registered seed, and the NaN-fold
    handling that keeps one degenerate block from voiding a resample;
  * the pre-registered band boundaries, evaluated at their exact edges.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("sklearn")
pytest.importorskip("rdkit")

from mammal_repurposing.cluster_a import allosteric_scaffold_split as S

ROOT = Path(__file__).resolve().parents[1]
CHEMBL = ROOT / "data" / "results" / "chembl_evidence.parquet"

# Two donepezil-like benzylpiperidines that differ only in substitution: the
# generic Murcko scaffold is deliberately coarse enough to merge them.
BENZYLPIPERIDINE_A = "O=C1c2cc(OC)c(OC)cc2CC1CC1CCN(Cc2ccccc2)CC1"
BENZYLPIPERIDINE_B = "O=C1c2cc(O)c(O)cc2CC1CC1CCN(Cc2ccccc2)CC1"


# ---------------------------------------------------------------------------
# Scaffold keying and its two sentinel cases
# ---------------------------------------------------------------------------

def test_scaffold_key_is_canonical_and_groups_analogues():
    ka = S.scaffold_key(BENZYLPIPERIDINE_A, "a")
    kb = S.scaffold_key(BENZYLPIPERIDINE_B, "b")
    assert ka and not S.is_sentinel(ka)
    assert ka == kb, "generic Murcko must merge the same skeleton across substituents"


def test_scaffold_key_distinguishes_different_skeletons():
    assert S.scaffold_key("c1ccccc1", "x") != S.scaffold_key("c1ccc2ccccc2c1", "y")


def test_unparseable_smiles_gets_unique_sentinel():
    k1 = S.scaffold_key("not_a_smiles((", "ROW1")
    k2 = S.scaffold_key("also))bad", "ROW2")
    assert k1.startswith("UNPARSEABLE::") and k2.startswith("UNPARSEABLE::")
    assert S.is_sentinel(k1) and S.is_sentinel(k2)
    # UNIQUE, not a shared bucket — otherwise every unparseable row would be
    # treated as one analogue series and strip training data for no reason.
    assert k1 != k2


def test_acyclic_molecule_gets_unique_sentinel():
    ach = S.scaffold_key("CC(=O)OCC[N+](C)(C)C", "acetylcholine")   # acetylcholine
    glu = S.scaffold_key("N[C@@H](CCC(=O)O)C(=O)O", "glutamate")    # glutamate
    assert ach.startswith("ACYCLIC::") and glu.startswith("ACYCLIC::")
    assert ach != glu, "acetylcholine and glutamate are not the same series"


def test_sentinels_never_match_a_real_scaffold():
    real = S.scaffold_key("c1ccccc1", "benzene")
    assert not S.is_sentinel(real)
    for s in (S.scaffold_key("bad((", "r"), S.scaffold_key("CCO", "r")):
        assert s != real


# ---------------------------------------------------------------------------
# The split itself — the property whose absence is the leak
# ---------------------------------------------------------------------------

def _toy() -> pd.DataFrame:
    """Two targets sharing a scaffold and sharing a compound — both leak paths."""
    df = pd.DataFrame({
        "compound_name": ["a", "b", "c", "d", "e", "f"],
        "target_uniprot": ["T1", "T1", "T1", "T2", "T2", "T2"],
        "inchikey": ["K_A", "K_B", "K_C", "K_A", "K_E", "K_F"],
        "smiles": [BENZYLPIPERIDINE_A, "c1ccccc1", "CCO",
                   BENZYLPIPERIDINE_A, BENZYLPIPERIDINE_B, "c1ccc2ccccc2c1"],
        "pact": [8.0, 7.0, 6.0, 9.0, 5.0, 7.5],
    })
    return S.assign_scaffold_keys(df)


def test_scaffold_disjoint_train_shares_no_scaffold_or_compound_with_test():
    df = _toy()
    test_idx = df.index[df["target_uniprot"] == "T2"]
    train = df[S.scaffold_disjoint_train_mask(df, test_idx)]
    test = df.loc[test_idx]

    real_test_scaffolds = {k for k in test["scaffold_key"] if not S.is_sentinel(k)}
    assert not (set(train["scaffold_key"]) & real_test_scaffolds)
    assert not (set(train["inchikey"]) & set(test["inchikey"]))
    assert "T2" not in set(train["target_uniprot"])


def test_scaffold_filter_actually_removes_rows_the_loto_split_keeps():
    """The leak, made concrete: leave-one-target-out keeps compound `a` (same
    InChIKey as a test row) and `b`-vs-`e` scaffold mates; Arm B must drop them."""
    df = _toy()
    test_idx = df.index[df["target_uniprot"] == "T2"]
    loto = df[df["target_uniprot"] != "T2"]
    scaf = df[S.scaffold_disjoint_train_mask(df, test_idx)]
    assert len(scaf) < len(loto), "scaffold filter must bite, or it is testing nothing"
    # `a` shares BOTH an InChIKey and a scaffold with test rows.
    assert "a" in set(loto["compound_name"])
    assert "a" not in set(scaf["compound_name"])


def test_acyclic_sentinel_does_not_strip_unrelated_acyclic_training_rows():
    """Regression guard: bucketing all acyclic molecules under one key would make
    ethanol and acetylcholine the same series and delete training rows for a leak
    that does not exist."""
    df = pd.DataFrame({
        "compound_name": ["ethanol", "acetylcholine"],
        "target_uniprot": ["T1", "T2"],
        "inchikey": ["K_ETOH", "K_ACH"],
        "smiles": ["CCO", "CC(=O)OCC[N+](C)(C)C"],
        "pact": [5.0, 6.0],
    })
    df = S.assign_scaffold_keys(df)
    test_idx = df.index[df["target_uniprot"] == "T2"]
    train = df[S.scaffold_disjoint_train_mask(df, test_idx)]
    assert set(train["compound_name"]) == {"ethanol"}


@pytest.mark.skipif(not CHEMBL.exists(), reason="chembl_evidence.parquet absent")
def test_real_data_split_is_disjoint_on_every_fold():
    ch = pd.read_parquet(CHEMBL)
    ch = ch[ch["best_pchembl"].notna()].copy()
    ch["target_uniprot"] = ch["target_uniprot"].astype(str)
    ch["pact"] = ch["best_pchembl"].astype(float)
    df = S.assign_scaffold_keys(ch.reset_index(drop=True))

    folds = S.target_folds(df, min_n=4)
    assert len(folds) >= 12, "the pre-registered headline floor needs >= 12 folds"

    bit = 0
    for t in folds:
        test_idx = df.index[df["target_uniprot"] == t]
        train = df[S.scaffold_disjoint_train_mask(df, test_idx)]
        test = df.loc[test_idx]
        real = {k for k in test["scaffold_key"] if not S.is_sentinel(k)}
        assert not (set(train["scaffold_key"]) & real), f"scaffold leak at {t}"
        assert not (set(train["inchikey"]) & set(test["inchikey"])), f"compound leak at {t}"
        assert len(train) >= S.MIN_TRAIN_ROWS_PER_FOLD, f"fold {t} starved"
        if len(train) < len(df) - len(test):
            bit += 1
    assert bit > 0, "on real data the scaffold filter must remove something somewhere"


# ---------------------------------------------------------------------------
# Blocked CV (Arm D)
# ---------------------------------------------------------------------------

def test_blocked_cv_keeps_scaffolds_within_one_fold_and_is_deterministic():
    df = _toy()
    a1 = S.blocked_cv_assignment(df, k=3)
    a2 = S.blocked_cv_assignment(df, k=3)
    assert a1.equals(a2), "assignment must be RNG-free and reproducible"
    for key, g in df.groupby("scaffold_key"):
        assert a1.loc[g.index].nunique() == 1, f"scaffold {key} crosses folds"


# ---------------------------------------------------------------------------
# Pooling, bootstrap and bands
# ---------------------------------------------------------------------------

def test_pooled_rho_is_sample_size_weighted_and_skips_undefined_folds():
    tab = pd.DataFrame({"fold": ["a", "b", "c"], "rho": [1.0, 0.0, np.nan], "n": [30, 10, 100]})
    # The NaN fold must not drag the mean toward zero, nor count in the weights.
    assert S.pooled_rho(tab) == pytest.approx((1.0 * 30 + 0.0 * 10) / 40)
    assert S.n_contributing_folds(tab) == 2


def test_bootstrap_is_deterministic_under_the_preregistered_seed():
    i1 = S.bootstrap_indices(19)
    i2 = S.bootstrap_indices(19)
    assert np.array_equal(i1, i2)
    assert i1.shape == (S.BOOTSTRAP_B, 19)
    assert S.BOOTSTRAP_SEED == 20260824 and S.BOOTSTRAP_B == 2000


def test_bootstrap_tolerates_folds_where_a_block_is_undefined():
    """A block that is degenerate on some folds (Boltz, wherever a target has no
    coverage) must pool over the rest, not void every resample containing one."""
    tab = pd.DataFrame({"fold": list("abcde"), "rho": [0.5, np.nan, 0.5, np.nan, 0.5],
                        "n": [10, 10, 10, 10, 10]})
    d = S.bootstrap_pooled(tab, S.bootstrap_indices(5, b=200))
    assert np.isfinite(d).mean() > 0.95
    assert np.nanmax(np.abs(d[np.isfinite(d)] - 0.5)) < 1e-9


def test_percentile_ci_is_the_preregistered_95_percent_interval():
    draws = np.linspace(0.0, 1.0, 10001)
    lo, hi = S.percentile_ci(draws)
    assert lo == pytest.approx(0.025, abs=1e-3)
    assert hi == pytest.approx(0.975, abs=1e-3)


@pytest.mark.parametrize("rho,lo,ref,expected", [
    (0.50, 0.30, 0.60, "SURVIVES"),    # clears all three clauses
    (0.45, 0.25, 0.60, "SURVIVES"),    # exactly on both edges, drop exactly 0.15
    (0.44, 0.30, 0.50, "DEGRADES"),    # below the 0.45 bar
    (0.50, 0.30, 0.70, "DEGRADES"),    # high absolute, but paired drop > 0.15
    (0.50, 0.20, 0.55, "DEGRADES"),    # lower bound under 0.25
    (0.25, 0.10, 0.30, "DEGRADES"),    # exactly on the FAIL floor -> not FAILS
    (0.24, 0.10, 0.30, "FAILS"),       # under the floor
    (0.60, 0.00, 0.65, "FAILS"),       # lower bound touches 0
])
def test_band_verdict_matches_the_preregistered_boundaries(rho, lo, ref, expected):
    verdict, _why = S.band_verdict(rho, lo, ref)
    assert verdict == expected


def test_block_interpretation_fires_the_preregistered_margin_rule():
    # Tanimoto within 0.05 of the fusion -> the "matches" headline is required.
    notes = S.block_interpretation({"F_fusion": 0.42, "T_tanimoto": 0.39, "F_minus_T": 0.40})
    assert any("matches the fusion" in n for n in notes)
    # Tanimoto above the fusion -> the stronger "net negative" reading.
    notes = S.block_interpretation({"F_fusion": 0.42, "T_tanimoto": 0.53, "F_minus_T": 0.40})
    assert any("EXCEEDS" in n for n in notes)
    # Fusion collapsing without Tanimoto -> the dependence must be quantified.
    notes = S.block_interpretation({"F_fusion": 0.42, "T_tanimoto": 0.53, "F_minus_T": 0.017})
    assert any("F-minus-T" in n for n in notes)
