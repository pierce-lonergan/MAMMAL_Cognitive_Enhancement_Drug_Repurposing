"""V6 §13.1 — Bayesian per-target router for the Multi-Head DTI ensemble.

Per Multi Head DTI.md §3: ensemble prediction
    ŷ(q, t) = Σ_k w_k(q, t) · y_k(q, t)
with
    w_k(q, t) ∝ T(t, k) · g(MD_k(q)) · h(σ_k(q, t))

where:
  T(t, k)     — trust matrix from §13.1 bias decomposition (diagnostics/per_head_bias.py)
  g(MD)       — OOD gate: exp(−max(0, MD − MD*)² / 2) with MD* = 99th-pctl
                of training Mahalanobis distances (diagnostics/ood_emosaic.py)
  h(σ)        — confidence gate: 1 / (σ² + ε), ε = 0.01 pchembl units

The router emits per-(query, target) (point estimate, lower, upper, w_vector)
with credible intervals propagated via Monte Carlo over per-head Venn-ABERS
predictive distributions (or Gaussian approximation when VA isn't available).

This is the V6 SKELETON — when 5 heads land (MAMMAL + Tanimoto already
shipped, MMAtt-DTA / PSICHIC / BALM pending), the router operationalises
immediately.

Reference:
  Multi Head DTI.md §3 — router specification + identifiability theorem
  Park et al. 2024 bioRxiv 2024.08.06.606753 — EnsDTI gating prior art
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class RouterPrediction:
    compound: str
    target: str
    y_hat: float                       # ensemble point estimate
    lo: float                          # 95% credible interval lower
    hi: float                          # 95% credible interval upper
    weights: dict[str, float]          # per-head weight in this prediction
    head_predictions: dict[str, float] # per-head point estimate
    head_sigmas: dict[str, float]      # per-head std
    n_active_heads: int                # heads with non-zero weight
    note: str = ""


def gate_ood(
    md: float,
    md_threshold: float,
    decay_scale: float = 1.0,
) -> float:
    """OOD gate g(MD) = exp(−max(0, MD − MD*)² / 2σ²).

    md_threshold is the 99th percentile of training Mahalanobis distances.
    Beyond threshold, weight decays exponentially with `decay_scale`.
    """
    if np.isnan(md):
        return 1.0    # missing OOD info → don't penalise
    excess = max(0.0, md - md_threshold)
    return float(np.exp(-(excess ** 2) / (2 * decay_scale ** 2)))


def gate_confidence(sigma: float, epsilon: float = 0.01) -> float:
    """Confidence gate h(σ) = 1 / (σ² + ε). Tighter σ → higher weight."""
    if np.isnan(sigma):
        return 1.0
    return float(1.0 / (sigma ** 2 + epsilon))


def route_one(
    compound: str,
    target: str,
    head_predictions: dict[str, float],
    head_sigmas: dict[str, float] | None,
    trust_row: dict[str, float],
    head_mds: dict[str, float] | None = None,
    md_threshold: float = 10.0,
    epsilon: float = 0.01,
) -> RouterPrediction:
    """Compute the ensemble prediction for one (compound, target).

    Args:
        head_predictions: {head: y_k(q, t)} per-head point estimates
        head_sigmas: {head: σ_k(q, t)} per-head standard deviations; None
            collapses to equal weighting
        trust_row: {head: T(t, k)} the trust-matrix row for this target
        head_mds: {head: Mahalanobis distance} for OOD gating; None disables
    """
    sigmas = head_sigmas or {}
    mds = head_mds or {}
    raw_weights: dict[str, float] = {}
    for head, y in head_predictions.items():
        if np.isnan(y):
            raw_weights[head] = 0.0
            continue
        T = float(trust_row.get(head, 0.0))
        g = gate_ood(mds.get(head, float("nan")), md_threshold)
        h = gate_confidence(sigmas.get(head, float("nan")), epsilon)
        raw_weights[head] = T * g * h

    total = sum(raw_weights.values())
    if total <= 0:
        return RouterPrediction(
            compound=compound, target=target,
            y_hat=float("nan"), lo=float("nan"), hi=float("nan"),
            weights={h: 0.0 for h in head_predictions},
            head_predictions=head_predictions,
            head_sigmas=sigmas,
            n_active_heads=0,
            note="all heads gated to zero",
        )
    norm_weights = {h: w / total for h, w in raw_weights.items()}
    y_hat = sum(norm_weights[h] * y for h, y in head_predictions.items()
                if not np.isnan(y))
    # Gaussian-approximated 95% CI under independence assumption
    var = sum(
        (norm_weights[h] ** 2) * (sigmas.get(h, 0.5) ** 2)
        for h in head_predictions
    )
    se = float(np.sqrt(var))
    return RouterPrediction(
        compound=compound, target=target,
        y_hat=float(y_hat),
        lo=float(y_hat - 1.96 * se),
        hi=float(y_hat + 1.96 * se),
        weights=norm_weights,
        head_predictions=head_predictions,
        head_sigmas=sigmas,
        n_active_heads=int(sum(1 for w in norm_weights.values() if w > 0.01)),
    )


def route_batch(
    long_predictions: pd.DataFrame,
    trust_matrix: pd.DataFrame,
    compound_col: str = "compound_name",
    target_col: str = "target_uniprot",
    head_col: str = "head",
    pred_col: str = "prediction",
    sigma_col: str | None = None,
    md_col: str | None = None,
    md_threshold: float = 10.0,
) -> pd.DataFrame:
    """Vectorised router over a long-format (compound, target, head, prediction).

    Returns a wide-format DataFrame: (compound, target, y_hat, lo, hi, w_<head>).
    """
    rows: list[dict] = []
    for (compound, target), sub in long_predictions.groupby(
            [compound_col, target_col]):
        head_preds = dict(zip(sub[head_col], sub[pred_col]))
        head_sigmas = (dict(zip(sub[head_col], sub[sigma_col]))
                       if sigma_col and sigma_col in sub.columns else None)
        head_mds = (dict(zip(sub[head_col], sub[md_col]))
                    if md_col and md_col in sub.columns else None)
        trust_row = (trust_matrix.loc[target].to_dict()
                     if target in trust_matrix.index else {})
        r = route_one(compound, target,
                      head_preds, head_sigmas, trust_row,
                      head_mds=head_mds, md_threshold=md_threshold)
        row = {
            compound_col: compound, target_col: target,
            "y_hat": r.y_hat, "lo": r.lo, "hi": r.hi,
            "n_active_heads": r.n_active_heads,
            "note": r.note,
        }
        for h, w in r.weights.items():
            row[f"w_{h}"] = w
        rows.append(row)
    return pd.DataFrame(rows)


def identifiability_diagnostic(
    trust_matrix: pd.DataFrame,
    per_target_n: dict[str, int],
    sigma_pchembl: float = 1.34,
    ci_width_target: float = 0.1,
) -> pd.DataFrame:
    """Per Multi Head DTI.md §3.4 identifiability theorem.

    Minimum n* for CI on w_k(t) to be < ci_width wide:
        n* ≈ 4·σ_pchembl² / ci_width²
    With v4 cognition panel n=7-26, this is FAR above the actual n. Honest
    interpretation: per-target router weights are not data-identified —
    they ARE priors, not posteriors.
    """
    n_star = 4 * (sigma_pchembl ** 2) / (ci_width_target ** 2)
    rows = []
    for target in trust_matrix.index:
        n_actual = per_target_n.get(target, 0)
        rows.append({
            "target": target,
            "n_actual": n_actual,
            "n_star_for_id": int(n_star),
            "identifiable": n_actual >= n_star,
            "trust_row_max": float(trust_matrix.loc[target].max()),
            "trust_row_min": float(trust_matrix.loc[target].min()),
            "trust_row_entropy": float(-np.sum(
                trust_matrix.loc[target].clip(1e-9, 1) *
                np.log(trust_matrix.loc[target].clip(1e-9, 1))
            )),
        })
    return pd.DataFrame(rows)
