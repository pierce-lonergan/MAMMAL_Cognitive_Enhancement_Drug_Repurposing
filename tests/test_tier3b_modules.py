"""Pytest coverage for the second Tier-3 sprint (commit 458881c onward):

  - §8.2 DrugComb fetcher (mocked HTTP; degradation path)
  - LambdaMART meta-ranker (synthetic features + NDCG sanity)
  - §7.15 hierarchical Bayes shrinkage fallback
  - §8.14 pocket-routed isotonic
  - §8.9 ANI-2x stub
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


# ---------------------------------------------------------------------------
# §8.2 DrugComb
# ---------------------------------------------------------------------------
class TestDrugComb:
    def test_graceful_when_unreachable(self):
        """fetch_drug_metadata returns empty metadata when both hosts fail."""
        from mammal_repurposing.fetchers.drugcomb import fetch_drug_metadata
        # Both DRUGCOMB_HOSTS are unreachable in this env — should NOT raise
        meta = fetch_drug_metadata("nonexistent_compound_xyz_12345")
        assert meta.name == "nonexistent_compound_xyz_12345"
        assert meta.drugcomb_id is None

    def test_summarise_synergy_empty(self):
        from mammal_repurposing.fetchers.drugcomb import summarise_synergy
        s = summarise_synergy(pd.DataFrame(), drug_a="x", partner_name="y")
        assert s.n_pairs == 0
        assert s.best_bliss is None

    def test_summarise_synergy_with_data(self):
        from mammal_repurposing.fetchers.drugcomb import summarise_synergy
        combos = pd.DataFrame({
            "drug_a": ["donepezil", "donepezil", "donepezil"],
            "drug_b": ["memantine", "memantine", "rivastigmine"],
            "cell_line": ["A", "B", "A"],
            "study": ["S1", "S1", "S2"],
            "synergy_bliss": [5.0, 8.0, 2.0],
            "synergy_loewe": [3.0, 6.0, 1.0],
        })
        s = summarise_synergy(combos, drug_a="donepezil", partner_name="memantine")
        assert s.n_pairs == 2
        assert s.best_bliss == 8.0
        assert s.best_loewe == 6.0
        assert s.median_bliss == 6.5


# ---------------------------------------------------------------------------
# LambdaMART meta-ranker
# ---------------------------------------------------------------------------
class TestLambdaMart:
    def test_discretize_pchembl(self):
        from mammal_repurposing.fusion.lambdamart_meta import discretize_pchembl
        # 25 evenly-spaced values → 5 buckets of 5 each
        pv = np.linspace(5.0, 10.0, 25)
        labels = discretize_pchembl(pv, n_buckets=5)
        # All buckets should be present
        assert set(labels) == {0, 1, 2, 3, 4}
        # Lowest values → 0; highest → 4
        assert labels[0] == 0
        assert labels[-1] == 4

    def test_ndcg_at_k_per_query(self):
        """Perfect ranking should give NDCG=1.0; reversed should give < 1.0."""
        from mammal_repurposing.fusion.lambdamart_meta import _ndcg_at_k_per_query
        df = pd.DataFrame({
            "target_uniprot": ["X"] * 5,
            "__gain": [4, 3, 2, 1, 0],
            "__pred": [5.0, 4.0, 3.0, 2.0, 1.0],   # perfect descending order
        })
        ndcg = _ndcg_at_k_per_query(df, k=5)
        assert abs(ndcg - 1.0) < 1e-6

        df["__pred"] = [1.0, 2.0, 3.0, 4.0, 5.0]    # reversed
        ndcg_bad = _ndcg_at_k_per_query(df, k=5)
        assert ndcg_bad < 1.0

    def test_fit_lambdamart_synthetic(self):
        """Train on a tiny synthetic frame; should produce a usable booster."""
        from mammal_repurposing.fusion.lambdamart_meta import (
            FEATURE_COLUMNS, fit_lambdamart, LambdaMartConfig,
        )
        # 4 targets × 10 compounds each = 40 rows; synthetic features
        rng = np.random.default_rng(0)
        rows = []
        for u in ["U1", "U2", "U3", "U4"]:
            for c in range(10):
                feats = rng.uniform(0, 1, len(FEATURE_COLUMNS))
                # Label correlated with first feature
                pchembl = 5.0 + 4.0 * feats[0] + rng.normal(0, 0.3)
                row = {
                    "compound_name": f"{u}_c{c}",
                    "target_uniprot": u,
                    "best_pchembl": pchembl,
                }
                for f, v in zip(FEATURE_COLUMNS, feats):
                    row[f] = v
                rows.append(row)
        df = pd.DataFrame(rows)
        res = fit_lambdamart(df, config=LambdaMartConfig(n_estimators=20))
        assert res.booster is not None
        assert res.n_train > 0
        assert res.n_test > 0
        # NDCG should be > 0 (synthetic signal exists)
        if res.test_ndcg_at_25 is not None:
            assert res.test_ndcg_at_25 > 0


# ---------------------------------------------------------------------------
# §7.15 Hierarchical Bayes
# ---------------------------------------------------------------------------
class TestHierarchicalBayes:
    def test_shrinkage_pulls_toward_family_mean(self):
        from mammal_repurposing.calibration.hierarchical_bayes import (
            empirical_bayes_shrinkage,
        )
        single = {"A": 0.8, "B": -0.4}
        n_per = {"A": 5, "B": 5}
        res = empirical_bayes_shrinkage("TEST", single, n_per)
        # Family mean = (0.8*5 + -0.4*5) / 10 = 0.2
        assert abs(res.family_mean_rho - 0.2) < 1e-6
        # Pooled A should be between single (0.8) and family mean (0.2)
        assert 0.2 <= res.pooled_rho["A"] <= 0.8
        # Pooled B should be between -0.4 and 0.2
        assert -0.4 <= res.pooled_rho["B"] <= 0.2

    def test_shrinkage_with_unbalanced_n(self):
        """Larger-n targets shrink less (the J-S idea)."""
        from mammal_repurposing.calibration.hierarchical_bayes import (
            empirical_bayes_shrinkage,
        )
        single = {"A": 0.8, "B": 0.0}
        n_per = {"A": 100, "B": 5}
        res = empirical_bayes_shrinkage("TEST", single, n_per)
        # A has more data → shrinks less → pooled rho closer to 0.8
        # B has less data → shrinks more → pooled rho closer to family mean
        shift_A = abs(res.pooled_rho["A"] - single["A"])
        shift_B = abs(res.pooled_rho["B"] - single["B"])
        assert shift_A < shift_B

    def test_fit_family_with_synthetic_data(self):
        from mammal_repurposing.calibration.hierarchical_bayes import fit_family
        rng = np.random.default_rng(7)
        data = {
            "T1": (rng.uniform(5, 9, 30), rng.uniform(5, 9, 30)),
            "T2": (rng.uniform(5, 9, 20), rng.uniform(5, 9, 20)),
        }
        res = fit_family("TEST", data, prefer_pymc=False)
        # Both targets present in output
        assert set(res.targets) == {"T1", "T2"}
        assert res.method == "empirical_bayes_shrinkage"


# ---------------------------------------------------------------------------
# §8.14 Pocket-routed isotonic
# ---------------------------------------------------------------------------
class TestPocketRouted:
    def test_routing_lift_on_distinct_pockets(self):
        """When two pockets have genuinely different slopes, routed beats global."""
        from mammal_repurposing.calibration.pocket_routed import (
            evaluate_routing_lift,
        )
        rng = np.random.default_rng(42)
        # S1: y = x; S2: y = -x + 20
        x1 = rng.uniform(5, 9, 20)
        y1 = x1 + rng.normal(0, 0.2, 20)
        x2 = rng.uniform(5, 9, 20)
        y2 = -x2 + 20 + rng.normal(0, 0.2, 20)
        raw = np.concatenate([x1, x2])
        truth = np.concatenate([y1, y2])
        pcls = np.array(["S1"] * 20 + ["S2"] * 20)
        lift = evaluate_routing_lift(raw, truth, pcls)
        # Per-pocket isotonic should significantly improve SSR over global
        assert lift["lift_pct"] > 10.0
        assert lift["routed_ssr"] < lift["global_ssr"]

    def test_predict_with_routing_falls_back(self):
        """Compounds with unknown pocket class get fallback calibrator output."""
        from mammal_repurposing.calibration.pocket_routed import (
            fit_pocket_routed, predict_with_routing,
        )
        rng = np.random.default_rng(0)
        x = rng.uniform(5, 9, 20)
        y = x + rng.normal(0, 0.1, 20)
        pcls = np.array(["S1"] * 10 + ["S2"] * 10)
        cal = fit_pocket_routed("TEST", x, y, pcls)
        # Predict on mixed: known + unknown
        preds, tags = predict_with_routing(cal,
                                            np.array([7.0, 7.0, 7.0]),
                                            np.array(["S1", "S2", "UNKNOWN"]))
        assert tags[0] == "pocket:S1"
        assert tags[1] == "pocket:S2"
        assert tags[2] == "fallback"

    def test_below_min_n_per_pocket_falls_back(self):
        from mammal_repurposing.calibration.pocket_routed import fit_pocket_routed
        x = np.array([5.0, 6.0, 7.0, 8.0, 9.0, 10.0])
        y = np.array([5.0, 6.0, 7.0, 8.0, 9.0, 10.0])
        pcls = np.array(["S1"] * 6)    # only one pocket, all in it
        cal = fit_pocket_routed("TEST", x, y, pcls, min_n_per_pocket=10)
        # S1 has 6 < 10 → no per-pocket fit; only fallback exists
        assert "S1" not in cal.by_pocket_class
        assert cal.fallback is not None


# ---------------------------------------------------------------------------
# §8.9 ANI-2x stub
# ---------------------------------------------------------------------------
class TestAniValidator:
    def test_stub_returns_deterministic_energy(self):
        from mammal_repurposing.scoring.ani2x_validator import (
            _compute_energy_stub, ANI_AVAILABLE,
        )
        # Stub formula: -50 kcal/mol per heavy atom
        e = _compute_energy_stub(["C", "C", "N", "O"], np.zeros((4, 3)))
        assert e == -200.0

    def test_validate_pose_set_runs(self, tmp_path):
        from mammal_repurposing.scoring.ani2x_validator import (
            summarise_validation, validate_pose_set,
        )
        # Mini mmCIF with 3 heavy ligand atoms
        cif = """data_t
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
HETATM 1 C C1 LIG L . 0.0 0.0 0.0 1.0 0.0
HETATM 2 N N1 LIG L . 1.5 0.0 0.0 1.0 0.0
HETATM 3 O O1 LIG L . 0.0 1.5 0.0 1.0 0.0
"""
        p = tmp_path / "pose.cif"
        p.write_text(cif, encoding="utf-8")
        results = validate_pose_set([p, p])
        assert len(results) == 2
        for r in results:
            assert r.n_atoms == 3
            assert r.status in ("STUB", "OK")    # depends on ANI install
        summary = summarise_validation(results)
        assert summary["n_total"] == 2
