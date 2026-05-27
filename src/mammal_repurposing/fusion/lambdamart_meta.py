"""LambdaMART meta-ranker — supervised on the 275 CORROBORATED ChEMBL labels.

Per V4 plan: with ≥20 labels per query (target), a learning-to-rank model
becomes appropriate. We use LightGBM's `lambdarank` objective (the canonical
LambdaMART implementation; Burges 2010 *Learning to Rank using Gradient
Descent and Beyond*; Wang et al. 2018 ICTIR *The LambdaLoss Framework*).

The training set comes from `data/results/chembl_evidence.parquet`:
    - 275 (compound, target) pairs with status='CORROBORATED'
    - `best_pchembl` is the per-pair relevance score
    - We discretize to 5 NDCG buckets (gain 0..4) by pchembl quintile

Features per (compound, target):
    - mammal_predicted_pkd (raw) and calibrated_pkd
    - cluster_a_tanimoto_score
    - cluster_a_boltzina_aff (when available; NaN otherwise)
    - cluster_b_admet_score
    - cluster_b_moa_score   (when MoA ranker run; else 0.5 default)
    - n_tier_1_liability   (from §8.0b-zn)
    - n_tier_2_liability
    - per-target Spearman ρ from Phase A.7 (target-level feature)

Hypothesis: LambdaMART NDCG@25 on held-out CORROBORATED labels ≥ RRF NDCG@25.

Reference:
  Burges 2010 — RankNet/LambdaRank/LambdaMART overview.
  Wang et al. 2018 ICTIR — LambdaLoss framework.
  Ke et al. 2017 NIPS — LightGBM.
"""

from __future__ import annotations

import logging
import pickle
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

try:
    import lightgbm as lgb
    HAS_LIGHTGBM = True
except ImportError:
    HAS_LIGHTGBM = False
    lgb = None    # type: ignore


@dataclass
class LambdaMartConfig:
    n_estimators: int = 200
    learning_rate: float = 0.05
    num_leaves: int = 31
    min_data_in_leaf: int = 2
    label_gain: list[int] = field(default_factory=lambda: [0, 1, 3, 7, 15])
    eval_at: list[int] = field(default_factory=lambda: [5, 10, 25])
    random_state: int = 42


def discretize_pchembl(pchembl_values: np.ndarray, n_buckets: int = 5) -> np.ndarray:
    """Discretize continuous pchembl into 0..n_buckets-1 NDCG gains.

    Uses quintile boundaries on the corroborated label distribution so each
    bucket has approximately equal mass. Returns int32 array.
    """
    pv = np.asarray(pchembl_values, dtype=float)
    valid = ~np.isnan(pv)
    if not valid.any():
        return np.zeros_like(pv, dtype=np.int32)
    # Quantile boundaries from the valid subset
    qs = np.linspace(0, 1, n_buckets + 1)
    edges = np.nanquantile(pv[valid], qs)
    # `digitize` returns 1..n; convert to 0..n-1
    labels = np.zeros_like(pv, dtype=np.int32)
    labels[valid] = np.clip(np.digitize(pv[valid], edges[1:-1], right=True),
                            0, n_buckets - 1)
    return labels


def build_feature_frame(
    evidence: pd.DataFrame,
    dti_calibrated: pd.DataFrame | None = None,
    boltzina: pd.DataFrame | None = None,
    admet: pd.DataFrame | None = None,
    moa: pd.DataFrame | None = None,
    liability_gates: pd.DataFrame | None = None,
    per_target_rho: dict[str, float] | None = None,
) -> pd.DataFrame:
    """Build the per-(compound, target) feature matrix from evidence + cluster outputs."""
    ev = evidence.copy()
    # Join calibrated MAMMAL (raw + calibrated)
    if dti_calibrated is not None and len(dti_calibrated):
        j = dti_calibrated.rename(
            columns={"predicted_pkd": "raw_pkd_dti",
                     "calibrated_pkd": "calibrated_pkd_dti"}
        )[["target_uniprot", "compound_name", "raw_pkd_dti",
           "calibrated_pkd_dti"]]
        ev = ev.merge(j, on=["target_uniprot", "compound_name"], how="left")
    # Boltzina affinity
    if boltzina is not None and len(boltzina):
        b_col = ("affinity_probability_binary"
                 if "affinity_probability_binary" in boltzina.columns
                 else "binder_prob")
        j = boltzina[["target_uniprot", "compound_name", b_col]].rename(
            columns={b_col: "boltzina_aff"}
        )
        ev = ev.merge(j, on=["target_uniprot", "compound_name"], how="left")
    # ADMET — per-compound (broadcast)
    if admet is not None and len(admet):
        a_col = ("admet_score" if "admet_score" in admet.columns else None)
        if a_col:
            j = admet[["compound_name", a_col]].rename(columns={a_col: "admet_score"})
            j["_jk"] = j["compound_name"].str.lower().str.strip()
            ev["_jk"] = ev["compound_name"].str.lower().str.strip()
            ev = ev.merge(j.drop(columns=["compound_name"]), on="_jk", how="left").drop(columns=["_jk"])
    # MoA per (compound, target)
    if moa is not None and len(moa):
        j = moa[["target_uniprot", "compound_name", "predicted_pkd"]].rename(
            columns={"predicted_pkd": "moa_score"}
        )
        ev = ev.merge(j, on=["target_uniprot", "compound_name"], how="left")
    # Liability counts per-compound
    if liability_gates is not None and len(liability_gates):
        j = liability_gates[["compound_name", "n_tier_1", "n_tier_2", "n_tier_3"]].rename(
            columns={"n_tier_1": "n_tier_1_liab", "n_tier_2": "n_tier_2_liab",
                     "n_tier_3": "n_tier_3_liab"}
        )
        j["_jk"] = j["compound_name"].str.lower().str.strip()
        ev["_jk"] = ev["compound_name"].str.lower().str.strip()
        ev = ev.merge(j.drop(columns=["compound_name"]), on="_jk", how="left").drop(columns=["_jk"])
    # Per-target Spearman ρ from Phase A.7
    if per_target_rho:
        ev["target_phase_a7_rho"] = ev["target_uniprot"].map(per_target_rho).fillna(0.0)
    return ev


FEATURE_COLUMNS = [
    "raw_pkd_dti",
    "calibrated_pkd_dti",
    "boltzina_aff",
    "admet_score",
    "moa_score",
    "n_tier_1_liab",
    "n_tier_2_liab",
    "n_tier_3_liab",
    "target_phase_a7_rho",
]


@dataclass
class LambdaMartTrainResult:
    booster: object
    n_train: int
    n_test: int
    train_targets: list[str]
    test_targets: list[str]
    test_ndcg_at_25: float | None
    feature_importance: dict[str, float]


def train_test_split_by_target(
    feature_frame: pd.DataFrame,
    test_frac: float = 0.25,
    seed: int = 42,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Split by TARGET (not by row) so test queries are held out entirely."""
    rng = np.random.default_rng(seed)
    targets = feature_frame["target_uniprot"].unique()
    rng.shuffle(targets)
    n_test = max(1, int(round(test_frac * len(targets))))
    test_targets = set(targets[:n_test])
    train = feature_frame[~feature_frame["target_uniprot"].isin(test_targets)].copy()
    test = feature_frame[feature_frame["target_uniprot"].isin(test_targets)].copy()
    return train, test


def fit_lambdamart(
    feature_frame: pd.DataFrame,
    label_col: str = "best_pchembl",
    feature_cols: list[str] | None = None,
    config: LambdaMartConfig | None = None,
    test_frac: float = 0.25,
    seed: int = 42,
) -> LambdaMartTrainResult:
    """Fit LightGBM LambdaMART (`lambdarank` objective) on the evidence frame.

    The frame must already contain the feature columns AND the label column.
    Groups are per `target_uniprot`. Test queries are held out as whole
    targets so we measure target-novel generalisation.

    Returns the trained Booster + held-out NDCG@25.
    """
    if not HAS_LIGHTGBM:
        raise ImportError("lightgbm is required for LambdaMART fit; install via `pip install lightgbm`")
    cfg = config or LambdaMartConfig()
    feats = feature_cols or FEATURE_COLUMNS

    # Drop rows with missing label
    df = feature_frame.dropna(subset=[label_col]).copy()
    # Replace remaining feature NaNs with median per column
    for c in feats:
        if c not in df.columns:
            df[c] = 0.0
        else:
            df[c] = df[c].fillna(df[c].median() if df[c].notna().any() else 0.0)

    # Discretize label to NDCG gain
    df["__gain"] = discretize_pchembl(df[label_col].to_numpy(), n_buckets=5)

    train, test = train_test_split_by_target(df, test_frac=test_frac, seed=seed)
    if len(train) == 0:
        raise ValueError("train split empty after target hold-out")

    # Sort train by group so LightGBM's group array works
    train = train.sort_values("target_uniprot").reset_index(drop=True)
    train_groups = train.groupby("target_uniprot").size().tolist()
    X_train = train[feats].to_numpy(dtype=float)
    y_train = train["__gain"].to_numpy(dtype=np.int32)

    lgb_train = lgb.Dataset(X_train, label=y_train, group=train_groups, feature_name=feats)
    params = dict(
        objective="lambdarank",
        learning_rate=cfg.learning_rate,
        num_leaves=cfg.num_leaves,
        min_data_in_leaf=cfg.min_data_in_leaf,
        label_gain=cfg.label_gain,
        eval_at=cfg.eval_at,
        verbose=-1,
        random_state=cfg.random_state,
    )
    booster = lgb.train(params, lgb_train, num_boost_round=cfg.n_estimators)

    # Held-out NDCG@25 on test targets
    ndcg25 = None
    if len(test):
        test = test.sort_values("target_uniprot").reset_index(drop=True)
        X_test = test[feats].to_numpy(dtype=float)
        test["__pred"] = booster.predict(X_test)
        ndcg25 = float(_ndcg_at_k_per_query(test, k=25))

    importance = dict(zip(
        booster.feature_name(),
        booster.feature_importance(importance_type="gain").tolist(),
    ))

    return LambdaMartTrainResult(
        booster=booster,
        n_train=int(len(train)),
        n_test=int(len(test)),
        train_targets=sorted(train["target_uniprot"].unique().tolist()),
        test_targets=sorted(test["target_uniprot"].unique().tolist()),
        test_ndcg_at_25=ndcg25,
        feature_importance=importance,
    )


def _dcg_at_k(gains: np.ndarray, k: int) -> float:
    g = gains[:k]
    if len(g) == 0:
        return 0.0
    discounts = 1.0 / np.log2(np.arange(2, len(g) + 2))
    return float(np.sum((np.power(2.0, g) - 1.0) * discounts))


def _ndcg_at_k_per_query(
    df: pd.DataFrame,
    k: int,
    score_col: str = "__pred",
    label_col: str = "__gain",
    group_col: str = "target_uniprot",
) -> float:
    """Mean NDCG@k across queries (one query per target_uniprot)."""
    out = []
    for _, sub in df.groupby(group_col):
        # Predicted ranking
        order = np.argsort(-sub[score_col].to_numpy())
        gains_pred = sub[label_col].to_numpy(dtype=int)[order]
        gains_ideal = np.sort(sub[label_col].to_numpy(dtype=int))[::-1]
        dcg = _dcg_at_k(gains_pred, k)
        idcg = _dcg_at_k(gains_ideal, k)
        if idcg == 0:
            continue
        out.append(dcg / idcg)
    return float(np.mean(out)) if out else 0.0


def ndcg_baseline_vs_lambdamart(
    feature_frame: pd.DataFrame,
    rrf_score_col: str = "raw_pkd_dti",   # baseline = single feature, e.g. MAMMAL pkd
    k: int = 25,
    label_col: str = "best_pchembl",
) -> dict:
    """Compute NDCG@k for the baseline single-feature ranking and report."""
    df = feature_frame.dropna(subset=[label_col]).copy()
    df["__gain"] = discretize_pchembl(df[label_col].to_numpy(), n_buckets=5)
    df["__pred"] = df[rrf_score_col].fillna(0.0)
    base_ndcg = _ndcg_at_k_per_query(df, k=k)
    return {"baseline_ndcg_at_k": base_ndcg, "k": k, "baseline_col": rrf_score_col}


def save_booster(result: LambdaMartTrainResult, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as f:
        pickle.dump({
            "booster": result.booster.model_to_string(),
            "n_train": result.n_train,
            "n_test": result.n_test,
            "test_ndcg_at_25": result.test_ndcg_at_25,
            "feature_importance": result.feature_importance,
        }, f)
