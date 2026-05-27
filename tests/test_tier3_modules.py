"""Pytest coverage for the Tier-3 sprint modules:

  - §8.0a Pareto NSGA-III non-dominated sort (fusion/pareto.py)
  - §7.12 Conformal prediction split (calibration/conformal.py)
  - §7.13 Scaffold-aware AL re-ranker (diagnostics/scaffold_aware_al.py)
  - §8.6 Brain-region annotator (analysis/brain_region.py)
  - Hypothesis-validation harness components (scripts/41)
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
# §8.0a Pareto NSGA-III
# ---------------------------------------------------------------------------
class TestParetoSort:
    def test_dominates_correctness(self):
        from mammal_repurposing.fusion.pareto import _dominates
        a = np.array([1.0, 1.0, 1.0])
        b = np.array([0.5, 0.5, 0.5])
        c = np.array([1.0, 1.0, 1.0])
        d = np.array([2.0, 0.0, 0.0])
        assert _dominates(a, b)         # strict domination
        assert not _dominates(b, a)
        assert not _dominates(a, c)     # equal — no strict dominance
        assert not _dominates(a, d)     # incomparable
        assert not _dominates(d, a)

    def test_pareto_front_has_no_dominator(self):
        """Brute-force: every rank-0 point is non-dominated by every other point."""
        from mammal_repurposing.fusion.pareto import non_dominated_sort
        # 10 random points in 3D
        rng = np.random.default_rng(7)
        A = rng.uniform(0, 1, size=(10, 3))
        ranks = non_dominated_sort(A)
        rank0 = [i for i, r in enumerate(ranks) if r == 0]
        # For every rank-0 point i, no other point j dominates it
        for i in rank0:
            for j in range(len(A)):
                if i == j:
                    continue
                # j dominates i means j >= i on all, > on some
                if np.all(A[j] >= A[i]) and np.any(A[j] > A[i]):
                    pytest.fail(f"rank-0 point {i} is dominated by {j}")

    def test_pareto_rank_monotone(self):
        """Rank should be consistent: each frontier strictly dominated only by lower-rank."""
        from mammal_repurposing.fusion.pareto import non_dominated_sort
        # Construct: (3,3) is the only undominated point. (3,2) and (2,3) are
        # dominated by (3,3) only — rank 1. (2,2) is dominated by (3,2) + (2,3)
        # (and (3,3)) — rank 2. And so on.
        A = np.array([
            [3, 3],                      # rank 0
            [3, 2], [2, 3],              # rank 1
            [2, 2], [2, 1], [1, 2],      # rank 2
            [1, 1], [1, 0], [0, 1],      # rank ≥ 3
        ], dtype=float)
        ranks = non_dominated_sort(A)
        assert ranks[0] == 0
        assert ranks[1] == 1
        assert ranks[2] == 1
        assert all(r >= 2 for r in ranks[3:6])
        assert all(r >= 3 for r in ranks[6:])

    def test_crowding_distance_extremes_are_inf(self):
        from mammal_repurposing.fusion.pareto import crowding_distance
        A = np.array([[0, 0], [1, 1], [2, 2], [3, 3]], dtype=float)
        cd = crowding_distance(A)
        # Extremes (sorted endpoints) get inf in each axis
        assert np.isinf(cd[0])
        assert np.isinf(cd[3])
        # Middle points get finite values
        assert np.isfinite(cd[1])
        assert np.isfinite(cd[2])

    def test_hypervolume_mc_monotone(self):
        """Adding a dominant point should not decrease hypervolume."""
        from mammal_repurposing.fusion.pareto import hypervolume_mc
        front1 = np.array([[1.0, 1.0]])
        front2 = np.array([[1.0, 1.0], [0.5, 1.5]])    # adds another non-dominated point
        ref = np.array([0.0, 0.0])
        hv1 = hypervolume_mc(front1, ref, n_samples=20_000, seed=0)
        hv2 = hypervolume_mc(front2, ref, n_samples=20_000, seed=0)
        assert hv2 >= hv1 - 0.02     # Monte Carlo noise tolerance

    def test_rank_pareto_populates_columns(self):
        """rank_pareto adds the expected columns even with minimal input."""
        from mammal_repurposing.fusion.pareto import ParetoConfig, rank_pareto
        df = pd.DataFrame({
            "compound_name": ["a", "b", "c"],
            "rrf_score": [0.9, 0.5, 0.7],
        })
        out = rank_pareto(df, ParetoConfig())
        assert "pareto_rank" in out.columns
        assert "crowding_distance" in out.columns
        assert "_axis_efficacy_rrf" in out.columns


# ---------------------------------------------------------------------------
# §7.12 Conformal prediction
# ---------------------------------------------------------------------------
class TestConformalCalibration:
    def test_fit_and_predict_roundtrip(self):
        from mammal_repurposing.calibration.conformal import (
            fit_inductive_conformal, predict_with_interval,
        )
        rng = np.random.default_rng(0)
        n = 100
        x = rng.uniform(4, 9, size=n)
        # Monotone-increasing relationship with noise
        y = 1.2 * x + rng.normal(0, 0.5, size=n)
        res = fit_inductive_conformal(x, y, "TEST", alpha=0.20, seed=0)
        assert res.q_alpha > 0
        assert res.n_cal > 0
        # Empirical coverage on cal fold should be ≥ nominal (80%)
        assert res.empirical_coverage >= 0.70

        # Predicted intervals contain the point estimate
        point, lo, hi = predict_with_interval(res, x[:10])
        assert np.all(lo <= point)
        assert np.all(point <= hi)
        assert np.all(hi - lo > 0)

    def test_insufficient_n_raises(self):
        from mammal_repurposing.calibration.conformal import fit_inductive_conformal
        with pytest.raises(ValueError, match="Need ≥10"):
            fit_inductive_conformal(np.array([1.0, 2.0]), np.array([1.0, 2.0]),
                                    "TEST", alpha=0.2)

    def test_q_alpha_from_loco(self):
        from mammal_repurposing.calibration.conformal import q_alpha_from_loco
        residuals = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
        # alpha=0.20 → rank = ceil(11 * 0.80) = 9 → 9th order stat = 0.9
        q = q_alpha_from_loco(residuals, alpha=0.20)
        assert abs(q - 0.9) < 1e-6


# ---------------------------------------------------------------------------
# §7.13 Scaffold-aware AL
# ---------------------------------------------------------------------------
class TestScaffoldAwareAL:
    def test_undersampled_scaffold_gets_bonus(self):
        from mammal_repurposing.diagnostics.scaffold_aware_al import (
            ScaffoldAwareConfig, rank_with_scaffold_bonus,
        )
        # Construct: 5 compounds with same scaffold (benzene rings) + 1 singleton
        df = pd.DataFrame({
            "compound_name": ["c1", "c2", "c3", "c4", "c5", "rare"],
            "rrf_score":     [0.9,   0.85, 0.80, 0.75, 0.70, 0.60],
            "compound_smiles": [
                "c1ccccc1C",       # benzene + methyl
                "c1ccccc1CC",      # benzene + ethyl
                "c1ccccc1CCC",     # benzene + propyl
                "c1ccccc1CCCC",    # benzene + butyl
                "c1ccccc1CCCCC",   # benzene + pentyl  — share benzene Murcko
                "C1CCCCCCCCCCC1",  # cyclododecane — unique Murcko
            ],
        })
        out = rank_with_scaffold_bonus(df, ScaffoldAwareConfig(alpha=0.5))
        rare = out[out["compound_name"] == "rare"].iloc[0]
        c1 = out[out["compound_name"] == "c1"].iloc[0]
        # rare has lower RRF but its scaffold density = 1 (singleton); the
        # exploration bonus should give it a higher AL score than expected
        # from RRF alone.
        assert rare["scaffold_density"] == 1
        assert c1["scaffold_density"] >= 4     # benzene scaffold shared by ≥4
        assert rare["scaffold_exploration_bonus"] > c1["scaffold_exploration_bonus"]

    def test_diversity_evaluation(self):
        from mammal_repurposing.diagnostics.scaffold_aware_al import evaluate_diversity
        baseline = pd.DataFrame({
            "compound_smiles": ["c1ccccc1", "c1ccccc1C", "c1ccccc1CC"],   # all benzene
        })
        al = pd.DataFrame({
            "compound_smiles": ["c1ccccc1", "C1CCCCC1", "C1CCCNC1"],      # 3 distinct
        })
        d = evaluate_diversity(baseline, al)
        assert d["baseline_n_distinct_scaffolds"] == 1
        assert d["al_n_distinct_scaffolds"] == 3
        assert d["delta"] == 2


# ---------------------------------------------------------------------------
# §8.6 Brain-region annotator
# ---------------------------------------------------------------------------
class TestBrainRegion:
    def test_panel_complete(self):
        from mammal_repurposing.analysis.brain_region import BRAIN_REGION_BIAS
        assert len(BRAIN_REGION_BIAS) == 22

    def test_annotate_known_target(self):
        from mammal_repurposing.analysis.brain_region import annotate
        df = pd.DataFrame([
            {"compound_name": "methylphenidate", "mammal_best_target": "Q01959"},
            {"compound_name": "modafinil",       "mammal_best_target": "Q01959"},
            {"compound_name": "ex1",             "mammal_best_target": "UNKNOWN"},
        ])
        out = annotate(df)
        assert out.iloc[0]["brain_bias"] == "brainstem"     # DAT → midbrain DA neurons
        assert "VTA" in out.iloc[0]["brain_primary_region"]
        # Unknown target gives empty annotation
        assert out.iloc[2]["brain_bias"] == ""

    def test_summary_counts_covers_all_categories(self):
        from mammal_repurposing.analysis.brain_region import summary_counts
        c = summary_counts()
        assert sum(c.values()) == 22
        # Expect ≥3 distinct categories present
        assert len(c) >= 3
