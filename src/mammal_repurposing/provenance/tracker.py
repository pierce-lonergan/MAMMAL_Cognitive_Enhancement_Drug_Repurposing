"""Per-compound provenance: one row capturing cluster ranks + gate status.

Schema (parquet):
    compound_name | smiles | evidence_tier |
    gate_status   | gates_failed | gates_flagged | regulatory_bypass |
    admet_score   |
    mammal_best_target | mammal_best_pkd | mammal_best_rank | mammal_polypharm_n |
    boltzina_best_target | boltzina_best_logic50 | boltzina_best_binder_prob | boltzina_best_rank |
    txgnn_max_indication | txgnn_max_contraindication | txgnn_max_rank |
    kg_path_count | kg_side_effect_count |
    rrf_score | rrf_rank | n_clusters_contributing
"""

from __future__ import annotations

import pandas as pd


def build_provenance(
    *,
    compounds: pd.DataFrame,
    gates: pd.DataFrame,
    mammal_scores: pd.DataFrame | None = None,
    boltzina_scores: pd.DataFrame | None = None,
    txgnn_scores: pd.DataFrame | None = None,
    kg_scores: pd.DataFrame | None = None,
    rrf_ranking: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Join every cluster's per-compound view into a single provenance table."""
    base = compounds[["name", "smiles", "evidence_tier"]].rename(
        columns={"name": "compound_name"}
    )

    out = base.merge(
        gates[["compound_name", "gate_status", "gates_failed", "gates_flagged",
               "regulatory_bypass", "admet_score"]],
        on="compound_name", how="left",
    )

    if mammal_scores is not None:
        m_per_compound = (
            mammal_scores.sort_values(["compound_name", "predicted_pkd"],
                                       ascending=[True, False])
            .groupby("compound_name")
            .agg(
                mammal_best_target=("target_uniprot", "first"),
                mammal_best_pkd=("predicted_pkd", "max"),
                mammal_polypharm_n=("target_uniprot", "nunique"),
            ).reset_index()
        )
        out = out.merge(m_per_compound, on="compound_name", how="left")

    if boltzina_scores is not None and not boltzina_scores.empty:
        b_per_compound = (
            boltzina_scores.sort_values(["compound_name", "binder_prob"],
                                         ascending=[True, False])
            .groupby("compound_name")
            .agg(
                boltzina_best_target=("target_uniprot", "first"),
                boltzina_best_logic50=("log_ic50", "min"),
                boltzina_best_binder_prob=("binder_prob", "max"),
            ).reset_index()
        )
        out = out.merge(b_per_compound, on="compound_name", how="left")

    if txgnn_scores is not None and not txgnn_scores.empty:
        out = out.merge(
            txgnn_scores[["compound_name", "indication_score", "contraindication_score"]]
            .rename(columns={"indication_score": "txgnn_max_indication",
                              "contraindication_score": "txgnn_max_contraindication"}),
            on="compound_name", how="left",
        )

    if kg_scores is not None and not kg_scores.empty:
        out = out.merge(
            kg_scores[["compound_name", "path_count", "side_effect_count"]]
            .rename(columns={"path_count": "kg_path_count",
                              "side_effect_count": "kg_side_effect_count"}),
            on="compound_name", how="left",
        )

    if rrf_ranking is not None and not rrf_ranking.empty:
        rrf_min = rrf_ranking[["compound_name", "per_compound_rrf"]].rename(
            columns={"per_compound_rrf": "rrf_score"}
        ).copy()
        rrf_min["rrf_rank"] = rrf_min["rrf_score"].rank(method="min", ascending=False).astype(int)
        out = out.merge(rrf_min, on="compound_name", how="left")

    contributing_cols = [
        ("mammal", "mammal_best_pkd"),
        ("boltzina", "boltzina_best_binder_prob"),
        ("txgnn", "txgnn_max_indication"),
        ("admet", "admet_score"),
        ("kg", "kg_path_count"),
    ]
    out["n_clusters_contributing"] = sum(
        out[c].notna().astype(int) if c in out.columns else 0
        for _, c in contributing_cols
    )

    return out
