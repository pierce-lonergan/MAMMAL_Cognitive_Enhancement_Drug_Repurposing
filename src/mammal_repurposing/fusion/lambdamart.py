"""LambdaMART promotion path (LightGBM) — used once ≥20 labeled positives exist.

Per the research doc §5: don't build this for v2 unless empirical validation
demands more. RRF is near-ceiling at low label counts. Keep this module as a
skeleton with the canonical feature schema; activate after labels accumulate.

Reference: Burges CJC. "From RankNet to LambdaRank to LambdaMART: An Overview."
Microsoft Research Tech Report MSR-TR-2010-82.

Promotion gate: ``len(positive_labels) >= weights.fusion.lambdamart_min_labels``
(default 20). Until then, the orchestrator routes to ``fusion.rrf`` only.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import pandas as pd

logger = logging.getLogger(__name__)


# Canonical feature schema for the meta-ranker.
# These columns are produced by the v2 phase scripts; this is the contract.
LAMBDAMART_FEATURES: list[str] = [
    "mammal_pkd",
    "mammal_polypharm_n",
    "boltzina_log_ic50",
    "boltzina_binder_prob",
    "admet_bbb",
    "admet_pgp_substrate",
    "admet_herg",
    "admet_dili",
    "admet_cyp3a4",
    "admet_cyp2d6",
    "admet_caco2_logpapp",
    "admet_score_composite",
    "txgnn_indication",
    "txgnn_contraindication",
    "kg_path_count_log",
    "kg_side_effect_penalty",
    "esm2_cos_to_canonical_binder",
    "morgan_cos_to_canonical_binder",
]


@dataclass
class LambdaMARTConfig:
    n_estimators: int = 300
    num_leaves: int = 31
    learning_rate: float = 0.05
    label_gain: list[int] = field(default_factory=lambda: [0, 1, 3])
    early_stopping_rounds: int | None = None
    feature_cols: list[str] = field(default_factory=lambda: LAMBDAMART_FEATURES.copy())


def fit_lambdamart(
    feature_df: pd.DataFrame,
    *,
    label_col: str = "label",
    group_col: str = "target_uniprot",
    config: LambdaMARTConfig | None = None,
):
    """Fit a LightGBM LambdaMART ranker. Requires lightgbm at runtime."""
    try:
        import lightgbm as lgb  # noqa: PLC0415
    except ImportError as e:
        raise ImportError(
            "lightgbm not installed. Run `pip install lightgbm` first."
        ) from e

    cfg = config or LambdaMARTConfig()
    df = feature_df.dropna(subset=cfg.feature_cols + [label_col, group_col]).copy()
    df = df.sort_values(group_col).reset_index(drop=True)

    groups = df.groupby(group_col).size().tolist()
    X = df[cfg.feature_cols].values
    y = df[label_col].astype(int).values

    ranker = lgb.LGBMRanker(
        objective="lambdarank",
        metric="ndcg",
        label_gain=cfg.label_gain,
        n_estimators=cfg.n_estimators,
        num_leaves=cfg.num_leaves,
        learning_rate=cfg.learning_rate,
    )
    ranker.fit(X, y, group=groups)
    logger.info("Fitted LambdaMART on %d rows across %d groups.",
                len(df), len(groups))
    return ranker


def predict_lambdamart(ranker, feature_df: pd.DataFrame,
                       *, feature_cols: list[str] | None = None) -> pd.Series:
    """Return predicted ranker scores; higher = more likely to be a hit."""
    cols = feature_cols or LAMBDAMART_FEATURES
    X = feature_df[cols].fillna(0.0).values
    preds = ranker.predict(X)
    return pd.Series(preds, index=feature_df.index, name="lambdamart_score")


def should_promote_to_lambdamart(n_positive_labels: int, min_required: int = 20) -> bool:
    """Gate condition for the RRF -> LambdaMART promotion (research doc §5)."""
    return n_positive_labels >= min_required
