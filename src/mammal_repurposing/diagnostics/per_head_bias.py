"""V6 §13.1 — Per-head bias decomposition for the Multi-Head DTI ensemble.

For each DTI head (MAMMAL, Tanimoto, MMAtt-DTA, PSICHIC, BALM) at each target,
compute the 4-dimensional bias signature vector b_k(t) = (PC_k, SN_k, OOD_k, CT_k):

  PC_k  — prior-collapse ratio = σ(predictions) / σ(training labels)
          SEVERE: < 0.3   MODERATE: 0.3–0.5   ACCEPTABLE: > 0.5
  SN_k  — scaffold-novelty bias = Spearman ρ between predictions and
          Morgan-FP T_max(query, training_actives)
          ≥ 0.6 = essentially a similarity searcher  ≤ 0.3 = scaffold-independent
  OOD_k — fraction of query compounds beyond the head's training distribution
          (Mahalanobis threshold at 99th percentile of training embeddings)
  CT_k  — calibration tier (A/B/C/D per §7.11) — already computed for MAMMAL

The 5-head × 22-target × 4-feature tensor feeds the §13.1 trust matrix T(t, k).

Reference:
  Multi Head DTI.md §2 — bias decomposition framework
  Shoshan et al. 2026 arXiv:2410.22367 — MAMMAL norm_y_mean/std reference
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

logger = logging.getLogger(__name__)


# Tier mapping for CT_k → numeric trust score
TIER_TO_SCORE: dict[str, float] = {"A": 1.0, "B": 0.7, "C": 0.4, "D": 0.1, "": 0.5}


@dataclass
class HeadBiasSignature:
    head: str
    target_uniprot: str
    n_predictions: int
    pc_ratio: float            # prior-collapse ratio
    sn_rho: float              # scaffold-novelty bias
    ood_fraction: float        # fraction beyond training distribution
    calibration_tier: str      # A/B/C/D
    trust_score: float = field(default=0.0)

    @property
    def pc_severity(self) -> str:
        if self.pc_ratio < 0.3:
            return "SEVERE"
        if self.pc_ratio < 0.5:
            return "MODERATE"
        return "ACCEPTABLE"


def compute_pc_ratio(
    predictions: np.ndarray,
    training_label_std: float = 1.34,    # MAMMAL training pchembl SD reference
) -> float:
    """σ(predictions) / σ(training_labels) — Multi Head DTI.md §2.2."""
    arr = np.asarray(predictions, dtype=float)
    arr = arr[~np.isnan(arr)]
    if len(arr) < 2:
        return float("nan")
    return float(arr.std() / training_label_std)


def compute_sn_bias(
    predictions: np.ndarray,
    tanimoto_to_actives: np.ndarray,
) -> float:
    """Spearman ρ between head predictions and Tanimoto-to-training-actives.

    SN ≥ 0.6 → head behaves like a similarity searcher
    SN ≤ 0.3 → head is scaffold-independent
    """
    p = np.asarray(predictions, dtype=float)
    t = np.asarray(tanimoto_to_actives, dtype=float)
    mask = ~(np.isnan(p) | np.isnan(t))
    if mask.sum() < 5:
        return float("nan")
    rho, _ = spearmanr(p[mask], t[mask])
    return float(rho) if rho is not None else float("nan")


def compute_ood_fraction(
    query_embeddings: np.ndarray,
    training_embeddings: np.ndarray,
    percentile: float = 99.0,
) -> float:
    """Fraction of query points beyond the training distribution's
    `percentile`-th Mahalanobis-distance threshold.

    Per Multi Head DTI.md §4 (eMOSAIC), the per-head OOD score is the
    min-Mahalanobis to the training cluster centroids. Here we use a
    simplified single-cluster variant: distance to training mean under
    pooled covariance.
    """
    Q = np.atleast_2d(query_embeddings)
    T = np.atleast_2d(training_embeddings)
    if T.shape[0] < 5 or Q.shape[0] < 1:
        return float("nan")
    mu = T.mean(axis=0)
    # atleast_2d: a single-feature embedding makes np.cov return a 0-d scalar, so cov.shape[0]
    # below would IndexError.
    cov = np.atleast_2d(np.cov(T, rowvar=False))
    # Regularised inverse (covariance often near-singular for small n)
    cov_reg = cov + 1e-4 * np.eye(cov.shape[0])
    try:
        cov_inv = np.linalg.inv(cov_reg)
    except np.linalg.LinAlgError:
        return float("nan")
    # Per-training-point distances → threshold = `percentile`-th
    train_d = np.array([
        float((row - mu) @ cov_inv @ (row - mu)) for row in T
    ])
    threshold = float(np.percentile(train_d, percentile))
    # Per-query distance
    query_d = np.array([
        float((row - mu) @ cov_inv @ (row - mu)) for row in Q
    ])
    return float((query_d > threshold).mean())


def compute_signature(
    head: str,
    target_uniprot: str,
    predictions: np.ndarray,
    tanimoto_to_actives: np.ndarray | None = None,
    query_embeddings: np.ndarray | None = None,
    training_embeddings: np.ndarray | None = None,
    calibration_tier: str = "",
    training_label_std: float = 1.34,
) -> HeadBiasSignature:
    """Compute the full (PC, SN, OOD, CT) signature for one (head, target)."""
    pc = compute_pc_ratio(predictions, training_label_std)
    sn = (compute_sn_bias(predictions, tanimoto_to_actives)
          if tanimoto_to_actives is not None else float("nan"))
    ood = (compute_ood_fraction(query_embeddings, training_embeddings)
           if (query_embeddings is not None and training_embeddings is not None)
           else float("nan"))
    return HeadBiasSignature(
        head=head,
        target_uniprot=target_uniprot,
        n_predictions=int(np.sum(~np.isnan(predictions))),
        pc_ratio=pc,
        sn_rho=sn,
        ood_fraction=ood,
        calibration_tier=calibration_tier,
    )


def build_trust_matrix(
    signatures: Iterable[HeadBiasSignature],
    alpha: float = 0.5,      # prior-collapse weight
    beta: float = 0.15,      # scaffold-novelty weight
    gamma: float = 0.15,     # OOD weight
    delta: float = 0.2,      # calibration-tier weight
    temperature: float = 0.5,
    floor: float = 0.02,
    ceiling: float = 0.7,
) -> pd.DataFrame:
    """Build the target-by-head trust matrix T(t, k) ∈ [floor, ceiling].

    Per Multi Head DTI.md §2.4: per-head trust score
        s_k(t) = α·(1 − |PC_k − 1|) + β·(1 − SN_k) + γ·(1 − OOD_k) + δ·tier(CT_k)
    then softmax over heads per target with temperature τ, clipped to
    [floor, ceiling].

    Returns DataFrame indexed by target_uniprot, columns = heads, values = T(t,k).
    Rows sum to ≤ 1 (after clipping; renormalised).
    """
    sig_df = pd.DataFrame([{
        "head": s.head, "target": s.target_uniprot,
        "pc": s.pc_ratio, "sn": s.sn_rho, "ood": s.ood_fraction,
        "tier": s.calibration_tier,
    } for s in signatures])
    if sig_df.empty:
        return pd.DataFrame()

    # Per-feature score component
    def _safe(v: float, default: float) -> float:
        return float(default if (v is None or np.isnan(v)) else v)

    sig_df["s_raw"] = sig_df.apply(lambda r: (
        alpha * (1.0 - abs(_safe(r["pc"], 0.5) - 1.0))
        + beta * (1.0 - _safe(r["sn"], 0.5))
        + gamma * (1.0 - _safe(r["ood"], 0.5))
        + delta * TIER_TO_SCORE.get(r["tier"], 0.5)
    ), axis=1)

    # Per-target softmax
    out_rows: list[dict] = []
    for target, sub in sig_df.groupby("target"):
        z = sub["s_raw"].to_numpy() / temperature
        z = z - z.max()
        w = np.exp(z)
        w = w / w.sum()
        # Clip + renormalise
        w = np.clip(w, floor, ceiling)
        w = w / w.sum()
        row = {"target": target}
        for h, wi in zip(sub["head"], w):
            row[h] = float(wi)
        out_rows.append(row)
    out_df = pd.DataFrame(out_rows).set_index("target")
    return out_df
