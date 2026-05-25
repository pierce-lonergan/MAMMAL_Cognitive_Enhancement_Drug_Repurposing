"""Prior-collapse sanity check — the FIRST diagnostic.

If MAMMAL predictions cluster tightly around norm_y_mean = 5.79
(BindingDB training prior) and have std << norm_y_std = 1.34, the
model is returning its prior, not learning a target-specific signal.

Any downstream rank analysis (Spearman ρ, etc.) computed on a collapsed
distribution is statistics on noise.

Source: Shoshan et al. 2026 npj Drug Discovery / arXiv:2410.22367 —
MAMMAL DTI head normalization constants `norm_y_mean=5.79384684128215,
norm_y_std=1.33808027428196`. The HuggingFace model card publishes these.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

NORM_Y_MEAN = 5.79384684128215
NORM_Y_STD = 1.33808027428196

# How much spread we want vs the training SD for predictions to be considered
# "non-collapsed." 0.5 = predictions span at least half the training SD.
NON_COLLAPSE_RATIO = 0.5
# IQR / SD ratio — if the bulk of predictions lives in a narrow window, the
# tails carry all the discrimination and ranking is fragile.
TIGHT_IQR_THRESHOLD = 0.25 * NORM_Y_STD


@dataclass
class PriorCollapseResult:
    target_uniprot: str
    target_gene: str
    n_predictions: int
    pred_mean: float
    pred_std: float
    pred_min: float
    pred_max: float
    pred_iqr: float
    distance_from_prior_mean: float          # |pred_mean - 5.79|
    collapse_ratio: float                    # NORM_Y_STD / pred_std (higher = worse)
    fraction_within_0_2_of_mean: float       # 99% collapse marker
    verdict: str                             # SEVERE | MODERATE | OK


def evaluate_target(target_uniprot: str, scores_df: pd.DataFrame) -> PriorCollapseResult:
    sub = scores_df[scores_df["target_uniprot"] == target_uniprot]
    if sub.empty:
        raise ValueError(f"No predictions for {target_uniprot}")

    pkd = sub["predicted_pkd"].to_numpy()
    gene = sub["target_gene"].iloc[0] if "target_gene" in sub.columns else "?"

    pred_mean = float(np.mean(pkd))
    pred_std = float(np.std(pkd, ddof=1)) if len(pkd) > 1 else 0.0
    pred_iqr = float(np.percentile(pkd, 75) - np.percentile(pkd, 25))
    pred_min = float(pkd.min())
    pred_max = float(pkd.max())
    distance_from_prior = abs(pred_mean - NORM_Y_MEAN)
    collapse_ratio = NORM_Y_STD / pred_std if pred_std > 0 else float("inf")
    fraction_within = float(np.mean(np.abs(pkd - pred_mean) <= 0.2))

    # Verdict — three tiers
    if collapse_ratio >= 10:
        verdict = "SEVERE"          # pred std < 1/10 of training SD; ranks are noise
    elif collapse_ratio >= 4:
        verdict = "MODERATE"        # constrained but discriminating
    else:
        verdict = "OK"

    return PriorCollapseResult(
        target_uniprot=target_uniprot,
        target_gene=gene,
        n_predictions=len(sub),
        pred_mean=pred_mean,
        pred_std=pred_std,
        pred_min=pred_min,
        pred_max=pred_max,
        pred_iqr=pred_iqr,
        distance_from_prior_mean=distance_from_prior,
        collapse_ratio=collapse_ratio,
        fraction_within_0_2_of_mean=fraction_within,
        verdict=verdict,
    )


def evaluate_panel(scores_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for u in scores_df["target_uniprot"].unique():
        r = evaluate_target(u, scores_df)
        rows.append({
            "target_uniprot": r.target_uniprot,
            "target_gene": r.target_gene,
            "n": r.n_predictions,
            "pred_mean": r.pred_mean,
            "pred_std": r.pred_std,
            "pred_range": r.pred_max - r.pred_min,
            "iqr": r.pred_iqr,
            "distance_from_prior_mean": r.distance_from_prior_mean,
            "collapse_ratio": r.collapse_ratio,
            "frac_within_0.2": r.fraction_within_0_2_of_mean,
            "verdict": r.verdict,
        })
    return pd.DataFrame(rows).sort_values("collapse_ratio", ascending=False)
