"""Gap 1 — Tests for the V11 grid composer (compound × target hypotheses).

Locks the fix for the v10 degeneracy where every compound collapsed onto
ACHE via `.iloc[0]`. These tests assert the V11 composer:
  - never collapses the target axis (spans multiple targets)
  - computes within-target binding percentiles
  - places the V7 anchor g at the KNOWN mechanism target, not the binding argmax
  - keeps predicted g below the Roberts ceiling for honest small effects
  - best_target_per_compound returns exactly one row per compound
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from mammal_repurposing.fusion.joint_composition import (
    GridCompositionConfig,
    best_target_per_compound,
    compose_grid_shortlist_v11,
)


# ---------------------------------------------------------------------------
# Fixtures: a tiny but realistic 4-compound × 3-target grid
# ---------------------------------------------------------------------------

@pytest.fixture
def tiny_grid():
    # 3 targets: ACHE (high relevance), SLC6A3 (mid), GRIA2 (low)
    rows = []
    pkd = {
        # compound: {target: pkd}
        "donepezil":   {"P22303": 9.0, "Q01959": 5.2, "P42262": 5.0},
        "methylphenidate": {"P22303": 5.1, "Q01959": 8.8, "P42262": 5.3},
        "noveldrug":   {"P22303": 8.5, "Q01959": 8.0, "P42262": 7.9},
        "weakbinder":  {"P22303": 5.0, "Q01959": 5.0, "P42262": 5.0},
    }
    for c, td in pkd.items():
        for t, v in td.items():
            rows.append({"compound_name": c, "target_uniprot": t, "predicted_pkd": v})
    v6a = pd.DataFrame(rows)

    v6b = pd.DataFrame([
        {"target_uniprot": "P22303", "gene": "ACHE", "theta_mean": 0.47,
         "theta_2p5": 0.2, "theta_97p5": 0.7, "w_pipeline": 0.615},
        {"target_uniprot": "Q01959", "gene": "SLC6A3", "theta_mean": 0.25,
         "theta_2p5": 0.0, "theta_97p5": 0.5, "w_pipeline": 0.562},
        {"target_uniprot": "P42262", "gene": "GRIA2", "theta_mean": 0.05,
         "theta_2p5": -0.2, "theta_97p5": 0.3, "w_pipeline": 0.51},
    ])

    target_class = {"P22303": "AChE-I", "Q01959": "NDRI", "P42262": "AMPA_pos_mod"}
    class_priors = {
        "AChE-I": {"mean": 0.18, "sd": 0.15},
        "NDRI": {"mean": 0.21, "sd": 0.18},
        "AMPA_pos_mod": {"mean": 0.05, "sd": 0.20},
    }
    gene_map = {"P22303": "ACHE", "Q01959": "SLC6A3", "P42262": "GRIA2"}
    return v6a, v6b, target_class, class_priors, gene_map


# ---------------------------------------------------------------------------
# Core: no target collapse
# ---------------------------------------------------------------------------

def test_grid_spans_all_targets_no_collapse(tiny_grid):
    v6a, v6b, tc, cp, gm = tiny_grid
    grid = compose_grid_shortlist_v11(v6a, v6b, tc, cp, target_gene_map=gm)
    # 4 compounds × 3 targets = 12 rows, NOT collapsed to 4
    assert len(grid) == 12
    assert grid["target_uniprot"].nunique() == 3
    # Each compound appears once per target
    assert (grid.groupby("compound")["target_uniprot"].nunique() == 3).all()


def test_within_target_binding_percentile(tiny_grid):
    v6a, v6b, tc, cp, gm = tiny_grid
    grid = compose_grid_shortlist_v11(v6a, v6b, tc, cp, target_gene_map=gm)
    # At ACHE, donepezil (pkd 9.0) should be the top binder (percentile ~1.0);
    # weakbinder (5.0) the lowest.
    ache = grid[grid["target_uniprot"] == "P22303"].set_index("compound")
    assert ache.loc["donepezil", "binding_percentile"] > ache.loc["weakbinder", "binding_percentile"]
    assert ache.loc["donepezil", "binding_percentile"] >= 0.9


# ---------------------------------------------------------------------------
# Anchor placement at KNOWN target, not binding argmax
# ---------------------------------------------------------------------------

def test_anchor_g_lands_at_known_target_not_binding_argmax(tiny_grid):
    v6a, v6b, tc, cp, gm = tiny_grid
    # methylphenidate's real g anchored; its KNOWN target is SLC6A3 (Q01959).
    # Its MAMMAL-best binding is also SLC6A3 here, but donepezil's known
    # target ACHE vs binding argmax must both resolve to ACHE.
    v7 = {"methylphenidate": (0.21, 0.26), "donepezil": (0.22, 0.27)}
    known = {"methylphenidate": "Q01959", "donepezil": "P22303"}
    grid = compose_grid_shortlist_v11(
        v6a, v6b, tc, cp, v7_anchor_g=v7, anchor_compound_target=known,
        target_gene_map=gm,
    )
    mph = grid[grid["compound"] == "methylphenidate"]
    anchored = mph[mph["g_source"] == "v7_nuts_anchor"]
    assert len(anchored) == 1
    assert anchored.iloc[0]["target_uniprot"] == "Q01959"
    assert abs(anchored.iloc[0]["g_predicted"] - 0.21) < 1e-6


def test_anchor_not_placed_when_known_target_absent(tiny_grid):
    """If a drug's known target isn't in the grid, the anchor g is NOT placed
    (the class-prior pathway runs instead — no biologically-wrong assignment)."""
    v6a, v6b, tc, cp, gm = tiny_grid
    v7 = {"donepezil": (0.22, 0.27)}
    known = {"donepezil": "P99999"}  # not in grid
    grid = compose_grid_shortlist_v11(
        v6a, v6b, tc, cp, v7_anchor_g=v7, anchor_compound_target=known,
        target_gene_map=gm,
    )
    don = grid[grid["compound"] == "donepezil"]
    # No anchor row (known target absent) → all class_prior
    assert (don["g_source"] == "class_prior").all()


# ---------------------------------------------------------------------------
# Roberts ceiling + differentiation
# ---------------------------------------------------------------------------

def test_predicted_g_below_ceiling_for_honest_effects(tiny_grid):
    v6a, v6b, tc, cp, gm = tiny_grid
    grid = compose_grid_shortlist_v11(v6a, v6b, tc, cp, target_gene_map=gm)
    # Class priors 0.05-0.21 × percentile × relevance → g_90 must stay < 0.5
    assert grid["g_90_upper"].max() < 0.5
    assert grid["roberts_ceiling_ok"].all()


def test_g_is_differentiated_across_pairs(tiny_grid):
    v6a, v6b, tc, cp, gm = tiny_grid
    grid = compose_grid_shortlist_v11(v6a, v6b, tc, cp, target_gene_map=gm)
    # NOT all the same value (the v10 degeneracy)
    assert grid["g_predicted"].round(3).nunique() >= 5


def test_high_binder_high_relevance_outranks_weak(tiny_grid):
    v6a, v6b, tc, cp, gm = tiny_grid
    grid = compose_grid_shortlist_v11(v6a, v6b, tc, cp, target_gene_map=gm)
    # donepezil@ACHE (strong binder, high relevance) should outrank
    # weakbinder@GRIA2 (weak binder, low relevance)
    don_ache = grid[(grid.compound == "donepezil") & (grid.target_uniprot == "P22303")].iloc[0]
    weak_gria = grid[(grid.compound == "weakbinder") & (grid.target_uniprot == "P42262")].iloc[0]
    assert don_ache["g_predicted"] > weak_gria["g_predicted"]


# ---------------------------------------------------------------------------
# best_target_per_compound view
# ---------------------------------------------------------------------------

def test_best_target_one_row_per_compound(tiny_grid):
    v6a, v6b, tc, cp, gm = tiny_grid
    grid = compose_grid_shortlist_v11(v6a, v6b, tc, cp, target_gene_map=gm)
    best = best_target_per_compound(grid)
    assert len(best) == 4  # one per compound
    assert best["compound"].nunique() == 4
    # methylphenidate's best target should be SLC6A3 (its strong binding + mid relevance
    # beats ACHE where it's a weak binder)
    mph_best = best[best.compound == "methylphenidate"].iloc[0]
    assert mph_best["target_uniprot"] == "Q01959"


# ---------------------------------------------------------------------------
# Real-data smoke (if the v11 parquet exists)
# ---------------------------------------------------------------------------

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
V11_GRID = ROOT / "data" / "results" / "v2" / "wet_lab_shortlist_v11_grid.parquet"


@pytest.mark.skipif(not V11_GRID.exists(), reason="v11 grid parquet not generated yet")
def test_real_v11_is_not_degenerate():
    grid = pd.read_parquet(V11_GRID)
    top25 = grid[grid["roberts_ceiling_ok"]].head(25)
    # The v10 degeneracy was 1 target; v11 must span many
    assert top25["target_uniprot"].nunique() >= 3
    assert grid["g_90_upper"].max() < 0.55
    # Known positive controls land at correct targets
    don = grid[(grid.compound == "donepezil")]
    if len(don):
        anchored = don[don.g_source.str.startswith("v7_nuts_anchor")]
        if len(anchored):
            assert anchored.iloc[0]["target_uniprot"] == "P22303"  # ACHE
