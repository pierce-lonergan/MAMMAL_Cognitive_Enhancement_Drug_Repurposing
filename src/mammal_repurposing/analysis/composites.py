"""Cognitive composite scoring panels.

Each composite is a weighted aggregate over a subset of the target panel.
Per the research doc Section 2 Step 6, we define three proxy panels — working
memory, processing speed, learning rate — plus a global composite (mean of three).

Composite math:
    1. For each (compound, target) score, compute z-score WITHIN the target's
       distribution (so panels with naturally-higher pKd targets don't dominate).
    2. Weighted mean of z-scores per panel.
    3. Global composite = unweighted mean of the three panel scores.

We also emit polypharm metrics:
    - polypharm_breadth   : # panel targets with raw pKd > 6.0
    - polypharm_weighted_score : sum of (pKd - 6.0) over panel targets with pKd > 6
"""

from __future__ import annotations

import pandas as pd

from mammal_repurposing.config import POLYPHARM_PKD_THRESHOLD


# Panels are keyed by gene symbol (we resolve gene -> uniprot from the targets df).
COMPOSITE_PANELS: dict[str, dict[str, float]] = {
    "working_memory": {
        "HRH3": 1.0, "ADRA2A": 1.0, "DRD1": 1.0,
        "CHRNA7": 1.0, "GRIN2B": 1.0,
    },
    "processing_speed": {
        "SLC6A3": 1.0, "SLC6A2": 1.0, "HRH3": 0.8,
        "HCRTR1": 0.6, "HCRTR2": 0.6,
    },
    "learning_rate": {
        "GRIA1": 1.0, "GRIA2": 1.0, "GRIA3": 1.0, "GRIA4": 1.0,
        "GRIN2A": 1.0, "NTRK2": 1.0, "PDE4D": 0.8,
    },
}


def _zscore_within_target(scores: pd.DataFrame) -> pd.DataFrame:
    """Add a ``pkd_z`` column with z-score of predicted_pkd within each target."""
    scores = scores.copy()
    g = scores.groupby("target_uniprot")["predicted_pkd"]
    scores["pkd_z"] = (scores["predicted_pkd"] - g.transform("mean")) / g.transform("std")
    # If std=0 (single value per target — degenerate), set z=0
    scores["pkd_z"] = scores["pkd_z"].fillna(0.0)
    return scores


def _gene_to_uniprot(targets: pd.DataFrame) -> dict[str, str]:
    return dict(zip(targets["gene"], targets["uniprot"]))


def compute_panel_score(
    scores_z: pd.DataFrame,
    panel: dict[str, float],
    gene_to_uniprot: dict[str, str],
) -> pd.Series:
    """Weighted-mean z-score per compound for one panel.

    Returns a Series indexed by compound_name.
    """
    weights_by_uniprot: dict[str, float] = {}
    for gene, w in panel.items():
        u = gene_to_uniprot.get(gene)
        if u is None:
            continue
        weights_by_uniprot[u] = w

    if not weights_by_uniprot:
        return pd.Series(dtype=float)

    sub = scores_z[scores_z["target_uniprot"].isin(weights_by_uniprot)].copy()
    sub["weight"] = sub["target_uniprot"].map(weights_by_uniprot)
    sub["weighted_z"] = sub["pkd_z"] * sub["weight"]

    per_compound = sub.groupby("compound_name").agg(
        sum_wz=("weighted_z", "sum"),
        sum_w=("weight", "sum"),
    )
    return (per_compound["sum_wz"] / per_compound["sum_w"]).rename("score")


def compute_polypharm(
    scores: pd.DataFrame,
    *,
    threshold: float = POLYPHARM_PKD_THRESHOLD,
) -> pd.DataFrame:
    """Per-compound polypharm breadth + weighted score (raw pKd, not z-scored)."""
    hits = scores[scores["predicted_pkd"] > threshold].copy()
    if hits.empty:
        return pd.DataFrame(columns=["compound_name", "polypharm_breadth",
                                     "polypharm_weighted_score"])

    hits["excess"] = hits["predicted_pkd"] - threshold
    return (
        hits.groupby("compound_name")
        .agg(
            polypharm_breadth=("target_uniprot", "nunique"),
            polypharm_weighted_score=("excess", "sum"),
        )
        .reset_index()
    )


def compute_composites(
    scores: pd.DataFrame,
    targets: pd.DataFrame,
    *,
    threshold: float = POLYPHARM_PKD_THRESHOLD,
) -> pd.DataFrame:
    """Per-compound DataFrame with all three panel scores + global composite
    + polypharm metrics."""
    scores_z = _zscore_within_target(scores)
    gene_map = _gene_to_uniprot(targets)

    panel_scores = {
        name: compute_panel_score(scores_z, panel, gene_map)
        for name, panel in COMPOSITE_PANELS.items()
    }

    df = pd.DataFrame(panel_scores).reset_index().rename(columns={"index": "compound_name"})
    # Global composite = unweighted mean across the three panels
    df["global_composite"] = df[list(COMPOSITE_PANELS.keys())].mean(axis=1)

    pp = compute_polypharm(scores, threshold=threshold)
    df = df.merge(pp, on="compound_name", how="left")
    df["polypharm_breadth"] = df["polypharm_breadth"].fillna(0).astype(int)
    df["polypharm_weighted_score"] = df["polypharm_weighted_score"].fillna(0.0)

    return df.sort_values("global_composite", ascending=False).reset_index(drop=True)
