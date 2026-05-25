"""Polypharmacology aggregator: count how many panel targets each compound
binds with predicted pKd above a threshold."""

from __future__ import annotations

import pandas as pd

from mammal_repurposing.config import POLYPHARM_PKD_THRESHOLD


def compute_polypharm(
    scores: pd.DataFrame,
    *,
    threshold: float = POLYPHARM_PKD_THRESHOLD,
) -> pd.DataFrame:
    """Return a per-compound table sorted by number of targets hit.

    Args:
        scores: DataFrame from ``data/results/dti_scores.parquet`` with columns
            ``compound_name``, ``target_uniprot``, ``predicted_pkd``.
        threshold: pKd above which a (compound, target) pair counts as a "hit".

    Returns:
        DataFrame with columns: compound_name, n_hits, mean_pkd_hits,
        max_pkd, targets_hit (semicolon-separated UniProt list).
    """
    hits = scores[scores["predicted_pkd"] > threshold].copy()

    if hits.empty:
        return pd.DataFrame(
            columns=["compound_name", "n_hits", "mean_pkd_hits", "max_pkd", "targets_hit"]
        )

    agg = (
        hits.groupby("compound_name")
        .agg(
            n_hits=("target_uniprot", "count"),
            mean_pkd_hits=("predicted_pkd", "mean"),
            max_pkd=("predicted_pkd", "max"),
            targets_hit=("target_uniprot", lambda s: ";".join(sorted(set(s)))),
        )
        .reset_index()
        .sort_values(["n_hits", "mean_pkd_hits"], ascending=[False, False])
        .reset_index(drop=True)
    )
    return agg
