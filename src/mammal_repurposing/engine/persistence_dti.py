"""L3-structural: read the persistence SUBSTRATE of any chemical from its predicted
engagement of a curated panel of substrate-defining targets.

PERSEUS's persistence head was, until now, curation-only: a compound with no entry in the
mechanism axis could only ABSTAIN. This module makes the substrate read STRUCTURE-COMPUTABLE
- score a compound against a small panel of targets whose engagement *defines* a substrate
class (senolytic BCL2/BCL-xL -> ablative; HDAC/DNMT/G9a/KEAP1 -> capability; TrkB ->
plasticity_window) and let the persistence head emit a substrate HYPOTHESIS from that.

The honesty discipline that makes this defensible (MAMMAL is a weak class router, so a raw
DTI score is NOT trustworthy on its own):

  1. CALIBRATION-GATED. A target may only contribute a substrate read if MAMMAL has been
     shown to RANK that target's known engagers above matched non-engagers - measured by
     per-target AUROC with a permutation p-value (scripts/104). A target that fails
     calibration contributes NOTHING; the compound abstains on that channel.
  2. ABLATIVE-ONLY PROMOTION. Only an ablative (senolytic) engagement is durable by
     construction. capability (epigenetic / NRF2) and plasticity_window engagements are
     reported as capability flags / window hypotheses, never auto-promoted to durable -
     consistent with the v2 reversibility design and the L4 permissive/instructive firewall.

This module is pure-Python for the calibration math + aggregation (testable with no GPU);
the single GPU entry point is :func:`score_compound_against_panel`.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

# substrate ordering: only ablative is durable-by-construction
TIER_RANK: dict[str, int] = {"capability": 1, "plasticity_window": 2, "ablative": 3}
DURABLE_TIER = "ablative"

DEFAULT_MIN_AUROC = 0.70
DEFAULT_MIN_POS = 3
DEFAULT_PERM_P = 0.05


# --------------------------------------------------------------------------
# calibration math (pure; NaN-dropping; small-n aware)
# --------------------------------------------------------------------------

def _clean(xs) -> np.ndarray:
    a = np.asarray(xs, dtype=float)
    return a[~np.isnan(a)]


def auroc(pos_scores, neg_scores) -> float:
    """Rank-based AUROC = P(score(pos) > score(neg)), ties counted as 0.5.

    NaNs are dropped first. Returns nan if either class is empty after cleaning.
    """
    pos, neg = _clean(pos_scores), _clean(neg_scores)
    if pos.size == 0 or neg.size == 0:
        return float("nan")
    pooled = np.concatenate([pos, neg])
    ranks = pd.Series(pooled).rank(method="average").to_numpy()
    r_pos = ranks[: pos.size].sum()
    u = r_pos - pos.size * (pos.size + 1) / 2.0
    return float(u / (pos.size * neg.size))


def permutation_p(pos_scores, neg_scores, *, n_perm: int = 2000, seed: int = 0) -> float:
    """One-sided permutation p-value for AUROC > 0.5 (label-shuffle null).

    Honest guard against small-n flukes: with 3-6 positives an AUROC of 1.0 is not rare by
    chance, so we report how often a random relabelling reaches the observed AUROC.
    """
    pos, neg = _clean(pos_scores), _clean(neg_scores)
    if pos.size == 0 or neg.size == 0:
        return float("nan")
    obs = auroc(pos, neg)
    pooled = np.concatenate([pos, neg])
    n_pos = pos.size
    rng = np.random.default_rng(seed)
    ge = 1  # +1 (observed) for an unbiased, never-zero estimate
    for _ in range(n_perm):
        perm = rng.permutation(pooled)
        if auroc(perm[:n_pos], perm[n_pos:]) >= obs - 1e-12:
            ge += 1
    return ge / (n_perm + 1)


def youden_threshold(pos_scores, neg_scores) -> float:
    """Score threshold t that maximises TPR(>=t) - FPR(>=t) (Youden's J).

    Used at inference: a compound's pKd at a calibrated target must clear this to count as
    an engagement. Returns nan if a class is empty.
    """
    pos, neg = _clean(pos_scores), _clean(neg_scores)
    if pos.size == 0 or neg.size == 0:
        return float("nan")
    cuts = np.unique(np.concatenate([pos, neg]))
    best_t, best_j = float(cuts[0]), -1.0
    for t in cuts:
        tpr = float((pos >= t).mean())
        fpr = float((neg >= t).mean())
        if tpr - fpr > best_j:
            best_j, best_t = tpr - fpr, float(t)
    return best_t


def specificity_threshold(pos_scores, neg_scores, *, target_fpr: float = 0.05) -> float:
    """Specificity-first threshold for SINGLE-COMPOUND inference: the score above which at
    most ``target_fpr`` of known non-engagers fall (the (1-target_fpr) quantile of negatives).

    Youden balances TPR/FPR and is fine for population ranking, but MAMMAL's pKd outputs are
    compressed, so a Youden cut lands inside the non-engager cloud and the channel fires for
    almost everyone. An abstain-by-default system needs the conservative cut instead: a novel
    compound must out-score essentially all known non-engagers to count as engaged.
    """
    pos, neg = _clean(pos_scores), _clean(neg_scores)
    if neg.size == 0:
        return float("nan")
    return float(np.quantile(neg, 1.0 - target_fpr))


def calibrate_target(
    pos_scores, neg_scores, *,
    min_auroc: float = DEFAULT_MIN_AUROC,
    min_pos: int = DEFAULT_MIN_POS,
    max_perm_p: float = DEFAULT_PERM_P,
    target_fpr: float = 0.05,
    min_sensitivity: float = 0.5,
    n_perm: int = 2000,
    seed: int = 0,
) -> dict:
    """Per-target calibration verdict with TWO gates:

      * ``passed`` - RANKING is valid: AUROC >= min_auroc AND permutation-p < max_perm_p AND
        enough scored positives. (Necessary, not sufficient for single-compound calls.)
      * ``inference_usable`` - the specificity-first threshold actually SEPARATES: at least
        ``min_sensitivity`` of engagers clear the (1-target_fpr) non-engager quantile. A
        channel can rank well yet be unusable per-compound when the score bands overlap.

    ``threshold`` (consumed at inference) is the specificity-first cut, NOT Youden.
    """
    pos, neg = _clean(pos_scores), _clean(neg_scores)
    a = auroc(pos, neg)
    p = permutation_p(pos, neg, n_perm=n_perm, seed=seed)
    spec_thr = specificity_threshold(pos, neg, target_fpr=target_fpr)
    sens = float((pos >= spec_thr).mean()) if pos.size and not math.isnan(spec_thr) else float("nan")
    passed = bool(
        pos.size >= min_pos
        and not math.isnan(a) and a >= min_auroc
        and not math.isnan(p) and p < max_perm_p
    )
    inference_usable = bool(passed and not math.isnan(sens) and sens >= min_sensitivity)
    return {
        "auroc": a, "perm_p": p,
        "n_pos_scored": int(pos.size), "n_neg_scored": int(neg.size),
        "threshold": spec_thr, "youden_threshold": youden_threshold(pos, neg),
        "sensitivity_at_threshold": sens,
        "passed": passed, "inference_usable": inference_usable,
    }


# --------------------------------------------------------------------------
# panel + calibration loaders
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class PanelTarget:
    gene: str
    uniprot: str
    sequence: str
    tier: str
    promotes_durable: bool


def load_panel(targets_csv: str | Path) -> dict[str, PanelTarget]:
    df = pd.read_csv(targets_csv)
    out: dict[str, PanelTarget] = {}
    for _, r in df.iterrows():
        out[r["gene"]] = PanelTarget(
            gene=r["gene"], uniprot=r["uniprot"],
            sequence=r.get("sequence", ""), tier=r["substrate_tier_implied"],
            promotes_durable=bool(r["promotes_durable"]),
        )
    return out


def load_calibration(calib_json: str | Path) -> dict[str, dict]:
    with open(calib_json, encoding="utf-8") as fh:
        return json.load(fh)["per_target"]


# --------------------------------------------------------------------------
# substrate-hypothesis aggregation (the head)
# --------------------------------------------------------------------------

@dataclass
class SubstrateHypothesis:
    substrate_hypothesis: str | None       # highest engaged tier, calibration-gated; None = abstain
    promotes_durable: bool                  # True only if an ablative target engaged + passed
    engaged: list[dict] = field(default_factory=list)   # [{gene, tier, pkd, threshold, promotes_durable}]
    capability_flags: list[str] = field(default_factory=list)
    abstained_targets: list[str] = field(default_factory=list)  # engaged-looking but calibration-failed
    note: str = ""


def substrate_hypothesis(
    compound_scores: dict[str, float],
    panel: dict[str, PanelTarget],
    calibration: dict[str, dict],
) -> SubstrateHypothesis:
    """Aggregate per-target DTI scores into a substrate HYPOTHESIS.

    A target contributes only if (a) it passed calibration and (b) the compound's pKd clears
    that target's calibrated Youden threshold. Durable promotion happens only for an ablative
    (senolytic) engagement. Everything else is a capability flag / window hypothesis.
    """
    engaged, flags, abstained = [], [], []
    for gene, pkd in compound_scores.items():
        if gene not in panel or pkd is None or (isinstance(pkd, float) and math.isnan(pkd)):
            continue
        cal = calibration.get(gene, {})
        thr = cal.get("threshold")
        looks_engaged = thr is not None and not math.isnan(thr) and pkd >= thr
        if not looks_engaged:
            continue
        # a channel must rank validly (passed) AND separate per-compound (inference_usable)
        if not (cal.get("passed", False) and cal.get("inference_usable", False)):
            abstained.append(gene)        # MAMMAL "engages" but the channel is untrustworthy
            continue
        t = panel[gene]
        engaged.append({"gene": gene, "tier": t.tier, "pkd": float(pkd),
                        "threshold": float(thr), "promotes_durable": t.promotes_durable})
        if t.tier != DURABLE_TIER:
            flags.append(f"{t.tier}:{gene}")

    if not engaged:
        note = ("no calibration-trusted substrate engagement"
                + (f"; {len(abstained)} engagement(s) on UN-calibrated targets ignored"
                   if abstained else ""))
        return SubstrateHypothesis(None, False, [], flags, abstained, note)

    top = max(engaged, key=lambda e: TIER_RANK.get(e["tier"], 0))
    durable = any(e["promotes_durable"] for e in engaged)
    note = (f"durable (ablative) engagement: {[e['gene'] for e in engaged if e['promotes_durable']]}"
            if durable else
            f"capability/window hypothesis only (no ablative engagement); flags={flags}")
    return SubstrateHypothesis(top["tier"], durable, engaged, flags, abstained, note)


# --------------------------------------------------------------------------
# GPU entry point (single compound vs the whole panel)
# --------------------------------------------------------------------------

def score_compound_against_panel(model, tokenizer, smiles: str,
                                 panel: dict[str, PanelTarget]) -> dict[str, float]:
    """MAMMAL-score one compound against every panel target. Returns {gene: pKd} (NaN on
    per-pair tokenizer/forward failure). This is the only function that needs the GPU."""
    from mammal_repurposing.scoring.dti import score_batch_safe  # noqa: PLC0415

    genes = list(panel)
    pairs = [(panel[g].sequence, smiles) for g in genes]
    ids = [f"{g}|persist" for g in genes]
    pkds = score_batch_safe(model, tokenizer, pairs, sample_ids=ids)
    return dict(zip(genes, pkds))
