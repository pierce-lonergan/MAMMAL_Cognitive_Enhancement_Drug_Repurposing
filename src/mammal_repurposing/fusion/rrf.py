"""Reciprocal Rank Fusion (RRF) — default fusion method for v2.

Math:  RRF_score(item) = Σ_k  weight_k / (k_const + rank_k(item))

where rank_k is the 1-indexed rank from ranker k (1 = best) and ``k_const``
is the Cormack et al. 2009 default of 60 (insensitive in [30, 80]).

Properties:
    - Rank-only: immune to score-distribution heterogeneity across rankers
      (Boltzina logIC50 is signed real; MAMMAL pKd is real; TxGNN is [0,1];
       ADMET-AI is [0,1]). Mixing these with CombSUM requires Platt scaling
       first; RRF doesn't.
    - Tolerant of missing rankers: if a compound has no rank from one ranker,
      its contribution from that ranker is zero (configurable).
    - Per-cluster weights supported via ``weights`` argument.

Source: Cormack, Clarke & Buettcher. "Reciprocal Rank Fusion Outperforms
Condorcet and Individual Rank Learning Methods." SIGIR 2009.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class RankerInput:
    """One ranker's contribution to RRF.

    Attributes:
        name: ranker identifier (e.g. "cluster_a_mammal", "cluster_b_admet").
        scores: Series indexed by item key, values are scores where HIGHER = BETTER.
                Items not in the index are treated as "not ranked by this model".
        ascending: if True, smaller scores are better (e.g. logIC50). Default False.
        weight: per-ranker weight applied to its 1/(k+rank) contribution.
    """

    name: str
    scores: pd.Series
    ascending: bool = False
    weight: float = 1.0


def rrf(
    rankers: list[RankerInput],
    *,
    k_const: int = 60,
    missing_rank_strategy: str = "skip",  # or "worst"
) -> pd.DataFrame:
    """Fuse multiple ranker scores into a single RRF-scored DataFrame.

    Args:
        rankers: list of RankerInput objects.
        k_const: RRF constant (Cormack default 60).
        missing_rank_strategy: how to treat items not ranked by a given ranker.
            "skip" = contribute 0 (recommended default for sparse-coverage rankers).
            "worst" = assign rank = N+1 where N is that ranker's coverage.

    Returns:
        DataFrame indexed by item, columns:
            rrf_score                    — sum of weighted contributions
            n_rankers_contributing       — count of rankers that placed this item
            rank_<name>                  — rank from each ranker (NaN if missing)
            contribution_<name>          — weighted 1/(k+rank) from each ranker
        sorted descending by rrf_score.
    """
    if not rankers:
        raise ValueError("RRF requires at least one ranker.")

    # Union of all item indices across rankers
    all_items = pd.Index([], dtype=object)
    for r in rankers:
        all_items = all_items.union(pd.Index(r.scores.dropna().index))

    rrf_score = pd.Series(0.0, index=all_items, name="rrf_score")
    n_contrib = pd.Series(0, index=all_items, name="n_rankers_contributing")
    rank_cols: dict[str, pd.Series] = {}
    contrib_cols: dict[str, pd.Series] = {}

    for r in rankers:
        s = r.scores.dropna()
        if s.empty:
            logger.warning("Ranker %s has no items; skipping.", r.name)
            continue
        # Higher score = better rank => sort descending then rank
        ranks = s.rank(method="min", ascending=r.ascending)
        # If we want missing items to get worst rank, fill from ranks.max()+1
        if missing_rank_strategy == "worst":
            worst = float(ranks.max()) + 1.0
            ranks = ranks.reindex(all_items, fill_value=worst)
        else:
            ranks = ranks.reindex(all_items)

        rank_cols[f"rank_{r.name}"] = ranks
        contrib = (r.weight / (k_const + ranks)).fillna(0.0)
        contrib_cols[f"contribution_{r.name}"] = contrib
        rrf_score = rrf_score.add(contrib, fill_value=0.0)
        n_contrib = n_contrib.add((~ranks.isna()).astype(int), fill_value=0)

    out = pd.DataFrame({
        **{name: col for name, col in rank_cols.items()},
        **{name: col for name, col in contrib_cols.items()},
        "rrf_score": rrf_score,
        "n_rankers_contributing": n_contrib.astype(int),
    })
    out = out.sort_values("rrf_score", ascending=False)
    out.index.name = "item"
    return out.reset_index()


def rrf_per_target_then_compound(
    long_scores: pd.DataFrame,
    *,
    item_col: str = "compound_name",
    target_col: str = "target_uniprot",
    score_col: str = "predicted_pkd",
    ranker_col: str = "ranker_name",
    ascending: bool = False,
    weight_map: dict[str, float] | None = None,
    k_const: int = 60,
) -> pd.DataFrame:
    """Per-target RRF over (target, compound) pairs, then aggregate per compound.

    Long-format input: rows are (target, compound, ranker_name, score).
    For each (target, compound) pair across rankers, compute the per-pair RRF.
    Then aggregate to per-compound by summing per-target RRF scores.

    Returns DataFrame with compound_name + per_compound_rrf + n_targets_supporting.
    """
    weight_map = weight_map or {}

    per_pair_rows: list[dict] = []
    for target, target_group in long_scores.groupby(target_col):
        rankers: list[RankerInput] = []
        for ranker_name, ranker_group in target_group.groupby(ranker_col):
            s = ranker_group.set_index(item_col)[score_col]
            rankers.append(RankerInput(
                name=str(ranker_name),
                scores=s,
                ascending=ascending,
                weight=weight_map.get(str(ranker_name), 1.0),
            ))
        if not rankers:
            continue
        pair_rrf = rrf(rankers, k_const=k_const)
        pair_rrf[target_col] = target
        per_pair_rows.append(pair_rrf)

    if not per_pair_rows:
        return pd.DataFrame(columns=["compound_name", "per_compound_rrf"])

    per_pair = pd.concat(per_pair_rows, ignore_index=True).rename(
        columns={"item": item_col}
    )
    agg = per_pair.groupby(item_col).agg(
        per_compound_rrf=("rrf_score", "sum"),
        n_targets_supporting=("rrf_score", "count"),
        mean_per_target_rrf=("rrf_score", "mean"),
    ).reset_index()
    return agg.sort_values("per_compound_rrf", ascending=False).reset_index(drop=True)
