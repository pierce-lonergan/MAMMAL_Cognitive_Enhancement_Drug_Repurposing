"""Faceted shortlist generator — by-mechanism-class + targeted-pair.

Per research/4-tier/Graczyk-Style ... .md §2. The point is to dissolve the
v3 HRH3-23/25 single-target lock-in by producing 8 mechanism-class top-5
tables + 9 targeted-pair top-5 tables, with cross-facet provenance so a
reviewer sees that donepezil appears in (cholinergic, CHRNA7+ACHE,
SIGMAR1+NTRK2) without triple-counting it as three independent hits.

Within each facet, ranking is:
    composite_facet = rrf_efficacy_normalized * (1 + alpha * gini_in_facet)
with alpha = 0.3 (small enough that polypharm survives, big enough that
mono-selective hits get a deserved boost).

This module is pure dataframe assembly — no ML.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from ..selectivity.gini_scorecard import MECHANISM_CLASS_TARGETS


# Targeted-pair facets per research doc §2.2. We hard-code the 9 pre-registered
# pairs to avoid combinatorial enumeration over C(22, 2) = 231 pairs.
TARGETED_PAIRS: dict[str, tuple[list[str], str]] = {
    "CHRNA7+ACHE":     (["CHRNA7", "ACHE"],
                        "Galantamine-class dual cholinergic"),
    "PDE4D+CHRNA7":    (["PDE4D", "CHRNA7"],
                        "cAMP + cholinergic LTP/spine convergence"),
    "HRH3+DRD1":       (["HRH3", "DRD1"],
                        "Dual aminergic for processing speed"),
    "GRIA+PDE4D":      (["GRIA1", "GRIA2", "GRIA3", "GRIA4", "PDE4D"],
                        "LTP via AMPA + cAMP convergence"),
    "SIGMAR1+NTRK2":   (["SIGMAR1", "NTRK2"],
                        "Neuroprotection axis (ANAVEX + 7,8-DHF combo logic)"),
    "DAT+NET":         (["SLC6A3", "SLC6A2"],
                        "Dual reuptake avoiding SERT (solriamfetol phenotype)"),
    "HCN1+KCNQ":       (["HCN1", "KCNQ2", "KCNQ3"],
                        "Intrinsic excitability tuning"),
    "GRIN2A_pref":     (["GRIN2A"],
                        "GluN2A-preferring procognitive PAMs"),
    "HCRTR1+DRD1":     (["HCRTR1", "DRD1"],
                        "Motivation/arousal axis"),
}


@dataclass
class FacetEntry:
    facet_type: str       # mechanism_class | targeted_pair
    facet_name: str
    facet_rank: int
    compound_name: str
    composite_score: float
    gini: float
    s_10x: int
    selectivity_category: str
    top_target: str
    notes: str


def _compute_gini_in_facet(
    pkd_in_facet: list[float],
    log_window: float = 1.0,
) -> float:
    """Local Gini computed on the subset of panel members IN this facet.

    For a mechanism-class facet with N targets, Gini-in-facet measures
    selectivity WITHIN the class — high value means one class member dominates.
    """
    arr = np.array([v for v in pkd_in_facet if not pd.isna(v)], dtype=float)
    if arr.size < 2:
        return 0.0
    m = arr.min()
    if m < 0:
        arr = arr - m
    arr = arr + 1e-9
    arr = np.sort(arr)
    n = arr.size
    idx = np.arange(1, n + 1)
    return float(((2 * idx - n - 1) * arr).sum() / (n * arr.sum()))


def build_by_class_facets(
    selectivity_df: pd.DataFrame,
    dti_grid: pd.DataFrame,
    target_gene_col: str = "target_gene",
    compound_col: str = "compound_name",
    pkd_col: str = "predicted_pkd",
    composite_col: str = "rrf_score",
    top_n: int = 5,
    alpha: float = 0.3,
    direction_notes: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Build the 8 mechanism-class facets, top-N each.

    Args:
        selectivity_df: output of score_panel() — one row per compound.
        dti_grid: long-format MAMMAL DTI grid (compound × target × pkd).
        composite_col: column on selectivity_df OR a sibling DataFrame
            containing the efficacy/RRF score per compound. If absent on
            selectivity_df, it'll fall back to top_target_pkd.
        direction_notes: optional dict {facet_name: note} for direction
            warnings (e.g. orexin antagonist vs cognition).
    """
    direction_notes = direction_notes or {}
    pivot = dti_grid.pivot_table(
        index=compound_col, columns=target_gene_col, values=pkd_col, aggfunc="first",
    )

    rows: list[FacetEntry] = []
    for facet_name, facet_genes in MECHANISM_CLASS_TARGETS.items():
        per_compound = []
        for compound in selectivity_df["compound_name"]:
            if compound not in pivot.index:
                continue
            pkd_in_facet = [pivot.loc[compound].get(g, np.nan) for g in facet_genes]
            pkd_in_facet = [v for v in pkd_in_facet if not pd.isna(v)]
            if not pkd_in_facet:
                continue
            max_in_facet = float(np.max(pkd_in_facet))
            mean_in_facet = float(np.mean(pkd_in_facet))
            gini_in_facet = _compute_gini_in_facet(pkd_in_facet)

            sel_row = selectivity_df.loc[selectivity_df["compound_name"] == compound].iloc[0]
            # Facet-internal composite: how strongly does this compound bind THIS
            # facet's targets, weighted by within-facet concentration. We use
            # max_in_facet (the top facet target's score) — NOT the global RRF —
            # because we want compounds that DOMINATE in this facet to surface.
            global_efficacy = sel_row.get(composite_col, np.nan)
            top_target = sel_row.get("top_target", "")
            efficacy = max_in_facet * (1 + alpha * gini_in_facet)

            # CLASS-CONCENTRATION BONUS (research doc §2.1): "bonus rewards compounds
            # whose Gini concentrates *on this class* (top-2 panel members are both
            # in this class)". We implement as a multiplicative boost when the
            # compound's panel-wide top target is IN this facet — i.e., this
            # compound's PRIMARY mechanism IS this class.
            class_match_bonus = 1.0
            if top_target in facet_genes:
                class_match_bonus = 1.10                 # +10% if top target is in-class
                # +5% more if second_target is also in-class
                second = sel_row.get("second_target", "")
                if second in facet_genes:
                    class_match_bonus = 1.15
            efficacy = efficacy * class_match_bonus

            # Optional tiebreaker via global score (small contribution)
            if not pd.isna(global_efficacy):
                efficacy = efficacy + 0.01 * float(global_efficacy)

            per_compound.append({
                "compound": compound,
                "composite": float(efficacy),
                "gini_in_facet": gini_in_facet,
                "max_pkd_in_facet": max_in_facet,
                "mean_pkd_in_facet": mean_in_facet,
                "gini": sel_row.get("gini", float("nan")),
                "s_10x": sel_row.get("s_10x", 0),
                "selectivity_category": sel_row.get("selectivity_category", "unknown"),
                "top_target": sel_row.get("top_target", ""),
            })

        per_compound.sort(key=lambda r: r["composite"], reverse=True)
        for rank, r in enumerate(per_compound[:top_n], 1):
            rows.append(FacetEntry(
                facet_type="mechanism_class",
                facet_name=facet_name,
                facet_rank=rank,
                compound_name=r["compound"],
                composite_score=r["composite"],
                gini=r["gini"],
                s_10x=int(r["s_10x"]),
                selectivity_category=r["selectivity_category"],
                top_target=r["top_target"],
                notes=direction_notes.get(facet_name, ""),
            ))
    return pd.DataFrame([r.__dict__ for r in rows])


def build_targeted_pair_facets(
    selectivity_df: pd.DataFrame,
    dti_grid: pd.DataFrame,
    target_gene_col: str = "target_gene",
    compound_col: str = "compound_name",
    pkd_col: str = "predicted_pkd",
    top_n: int = 5,
    lambda_off_pair: float = 0.5,
) -> pd.DataFrame:
    """Build the 9 pre-registered targeted-pair facets.

    Within-pair ranking: `(pkd_A + pkd_B) − λ * max(pkd_off_pair)`.
    """
    pivot = dti_grid.pivot_table(
        index=compound_col, columns=target_gene_col, values=pkd_col, aggfunc="first",
    )

    panel_set = set()
    for genes in MECHANISM_CLASS_TARGETS.values():
        panel_set.update(genes)

    rows: list[FacetEntry] = []
    for facet_name, (facet_genes, hypothesis) in TARGETED_PAIRS.items():
        per_compound = []
        off_genes = list(panel_set - set(facet_genes))
        for compound in selectivity_df["compound_name"]:
            if compound not in pivot.index:
                continue
            in_pair = np.array(
                [pivot.loc[compound].get(g, np.nan) for g in facet_genes], dtype=float,
            )
            in_pair = in_pair[~np.isnan(in_pair)]
            if in_pair.size == 0:
                continue
            off_pair = np.array(
                [pivot.loc[compound].get(g, np.nan) for g in off_genes], dtype=float,
            )
            off_pair = off_pair[~np.isnan(off_pair)]
            score = float(in_pair.sum()) - lambda_off_pair * (
                float(off_pair.max()) if off_pair.size > 0 else 0.0
            )

            sel_row = selectivity_df.loc[selectivity_df["compound_name"] == compound].iloc[0]
            per_compound.append({
                "compound": compound,
                "composite": score,
                "gini": sel_row.get("gini", float("nan")),
                "s_10x": sel_row.get("s_10x", 0),
                "selectivity_category": sel_row.get("selectivity_category", "unknown"),
                "top_target": sel_row.get("top_target", ""),
                "hypothesis": hypothesis,
            })

        per_compound.sort(key=lambda r: r["composite"], reverse=True)
        for rank, r in enumerate(per_compound[:top_n], 1):
            rows.append(FacetEntry(
                facet_type="targeted_pair",
                facet_name=facet_name,
                facet_rank=rank,
                compound_name=r["compound"],
                composite_score=r["composite"],
                gini=r["gini"],
                s_10x=int(r["s_10x"]),
                selectivity_category=r["selectivity_category"],
                top_target=r["top_target"],
                notes=r["hypothesis"],
            ))
    return pd.DataFrame([r.__dict__ for r in rows])


def compute_cross_facet_provenance(faceted_df: pd.DataFrame) -> pd.DataFrame:
    """For each compound, list every facet it appears in.

    Returns a DataFrame: compound_name, cross_facet_count, cross_facet_list.
    Used to stop a reviewer from triple-counting one compound as if it were
    three independent hits.
    """
    if faceted_df.empty:
        return pd.DataFrame(columns=["compound_name", "cross_facet_count", "cross_facet_list"])

    rows = []
    for compound, sub in faceted_df.groupby("compound_name"):
        labels = [f"{r['facet_name']} #{int(r['facet_rank'])}"
                  for _, r in sub.iterrows()]
        rows.append({
            "compound_name": str(compound),
            "cross_facet_count": int(len(labels)),
            "cross_facet_list": "; ".join(labels),
        })
    return pd.DataFrame(rows).sort_values("cross_facet_count", ascending=False)
