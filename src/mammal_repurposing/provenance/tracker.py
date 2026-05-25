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
        # Real schema: affinity_probability_binary (0-1, higher = better) and
        # affinity_pred_value (logIC50 in µM; lower = better). Older smoke files
        # may still have the legacy binder_prob/log_ic50 names.
        bz = boltzina_scores.copy()
        prob_col = ("affinity_probability_binary"
                    if "affinity_probability_binary" in bz.columns
                    else "binder_prob")
        logic50_col = ("affinity_pred_value"
                       if "affinity_pred_value" in bz.columns
                       else "log_ic50")
        b_per_compound = (
            bz.sort_values(["compound_name", prob_col], ascending=[True, False])
            .groupby("compound_name")
            .agg(
                boltzina_best_target=("target_uniprot", "first"),
                boltzina_best_logic50=(logic50_col, "min"),
                boltzina_best_binder_prob=(prob_col, "max"),
            ).reset_index()
        )
        out = out.merge(b_per_compound, on="compound_name", how="left")

    if txgnn_scores is not None and not txgnn_scores.empty:
        # Real Cluster C parquet uses txgnn_mean_p_indication /
        # txgnn_mean_p_contraindication. Older stubs use indication_score /
        # contraindication_score. Pick whichever pair exists.
        ind_col = ("txgnn_mean_p_indication"
                   if "txgnn_mean_p_indication" in txgnn_scores.columns
                   else "indication_score")
        contra_col = ("txgnn_mean_p_contraindication"
                      if "txgnn_mean_p_contraindication" in txgnn_scores.columns
                      else "contraindication_score")
        keep = ["compound_name"]
        rename: dict[str, str] = {}
        if ind_col in txgnn_scores.columns:
            keep.append(ind_col)
            rename[ind_col] = "txgnn_max_indication"
        if contra_col in txgnn_scores.columns:
            keep.append(contra_col)
            rename[contra_col] = "txgnn_max_contraindication"
        if len(keep) > 1:
            out = out.merge(txgnn_scores[keep].rename(columns=rename),
                            on="compound_name", how="left")

    if kg_scores is not None and not kg_scores.empty:
        # Real Cluster C parquet uses kg_ppr_sum / kg_n_targets_reachable.
        # Older stubs use path_count / side_effect_count.
        pc_col = ("kg_ppr_sum" if "kg_ppr_sum" in kg_scores.columns
                  else ("path_count" if "path_count" in kg_scores.columns else None))
        se_col = ("kg_n_targets_reachable"
                  if "kg_n_targets_reachable" in kg_scores.columns
                  else ("side_effect_count"
                        if "side_effect_count" in kg_scores.columns else None))
        keep = ["compound_name"]
        rename = {}
        if pc_col:
            keep.append(pc_col)
            rename[pc_col] = "kg_path_count"
        if se_col:
            keep.append(se_col)
            rename[se_col] = "kg_side_effect_count"
        if len(keep) > 1:
            out = out.merge(kg_scores[keep].rename(columns=rename),
                            on="compound_name", how="left")

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
