"""V6.A.5 — Multi-head disagreement axis (generalises §8.15 from 2 → N heads).

§8.15 (`scripts/35_v3_disagreement_signal.py`) computes per-compound
rank disagreement between MAMMAL and Tanimoto. V6.A.5 generalises:

  - For each (compound, target) pair, compute the rank under EVERY head
  - Pairwise Kendall-τ between rank vectors per target
  - Per-compound disagreement entropy across heads
  - Facet-tag per compound: {novel_scaffold | activity_cliff | ood | noise}

Per Multi Head DTI.md §6, disagreement IS the discovery signal — compounds
where heads disagree are exactly the wet-lab-priority candidates that
single-head ranking would miss.

Heads (currently 3): MAMMAL_cal, Tanimoto, PrimeKG_PPR. When V6.A.1
activates MMAtt-DTA, this extends to 4 with no script changes.

Output:
  data/results/v2/disagreement_axis_v1.parquet
  reports/pipeline/disagreement_axis_v1.md
"""

from __future__ import annotations

import argparse
import logging
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kendalltau

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("v6_disagreement")

# Per Multi Head DTI.md §6.2 facet thresholds
NOVEL_SCAFFOLD_DELTA = 50      # rank delta threshold
ACTIVITY_CLIFF_DELTA = 50
HIGH_DISAGREEMENT_ENTROPY = 0.7


def _normalise_to_rank(values: np.ndarray) -> np.ndarray:
    """Returns 1-based ranks (lower = higher score). NaN-safe."""
    arr = np.asarray(values, dtype=float)
    # Argsort descending (higher score → lower rank)
    valid = ~np.isnan(arr)
    out = np.full(arr.shape, np.nan)
    if valid.sum() == 0:
        return out
    order = np.argsort(-arr[valid], kind="stable")
    ranks = np.empty(valid.sum())
    ranks[order] = np.arange(1, valid.sum() + 1)
    out[valid] = ranks
    return out


def _rank_entropy(rank_row: np.ndarray, n_compounds: int) -> float:
    """Normalised entropy of a compound's ranks across heads.

    rank_row: 1D array of ranks (one per head). NaN-skip.
    Higher entropy ≈ heads disagree more on this compound's position.
    """
    valid = rank_row[~np.isnan(rank_row)]
    if len(valid) < 2:
        return float("nan")
    # Convert ranks to position [0, 1] (rank 1 → 0, last → 1)
    pos = (valid - 1) / max(n_compounds - 1, 1)
    # Simple proxy: std of positions, normalised by max possible (0.5)
    return float(min(pos.std() / 0.5, 1.0))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path,
                        default=ROOT / "data/results/v2/disagreement_axis_v1.parquet")
    parser.add_argument("--report", type=Path,
                        default=ROOT / "reports/pipeline/disagreement_axis_v1.md")
    args = parser.parse_args()

    # Load each head's predictions per (compound, target)
    dti_cal = pd.read_parquet(ROOT / "data/results/dti_scores_calibrated.parquet")
    mammal_long = dti_cal[["target_uniprot", "compound_name", "calibrated_pkd"]].rename(
        columns={"calibrated_pkd": "score"}
    )
    mammal_long["head"] = "MAMMAL_cal"

    disag_path = ROOT / "data/results/v2/disagreement_signal.parquet"
    tanimoto_long = pd.DataFrame()
    if disag_path.exists():
        d = pd.read_parquet(disag_path)
        tanimoto_long = d[["target_uniprot", "compound_name", "tanimoto_score"]].rename(
            columns={"tanimoto_score": "score"}
        )
        tanimoto_long["head"] = "Tanimoto"

    kg = pd.read_parquet(ROOT / "data/results/v2/kg_scores.parquet")
    # PrimeKG PPR is per-compound; broadcast per target
    targets_uni = mammal_long["target_uniprot"].unique()
    kg_long_rows = []
    for u in targets_uni:
        for _, r in kg.iterrows():
            kg_long_rows.append({
                "target_uniprot": u,
                "compound_name": r["compound_name"],
                "score": float(r["kg_ppr_sum"]),
                "head": "PrimeKG_PPR",
            })
    kg_long = pd.DataFrame(kg_long_rows)

    # V6.A.1 — MMAtt-DTA per-(compound, target) predictions
    mmatt_path = ROOT / "data/results/v2/mmatt_dta_predictions.parquet"
    mmatt_long = pd.DataFrame()
    if mmatt_path.exists():
        mm = pd.read_parquet(mmatt_path)
        # Join SMILES → compound_name from compounds parquet
        compounds = pd.read_parquet(ROOT / "data/interim/compounds.parquet")
        smi_to_name = dict(zip(compounds["smiles"], compounds["name"]))
        mm["compound_name"] = mm["smiles"].map(smi_to_name)
        mm = mm.dropna(subset=["compound_name"])
        mmatt_long = mm[["uniprot_id", "compound_name", "prediction"]].rename(
            columns={"uniprot_id": "target_uniprot", "prediction": "score"}
        )
        mmatt_long["head"] = "MMAtt_DTA"

    all_long = pd.concat(
        [mammal_long, tanimoto_long, kg_long, mmatt_long], ignore_index=True
    )
    heads = sorted(all_long["head"].unique())
    logger.info("Heads available: %s", heads)

    # Per-target: compute per-head ranks + pairwise Kendall-τ
    pairwise_rows: list[dict] = []
    per_compound_rows: list[dict] = []

    for target_uni, sub in all_long.groupby("target_uniprot"):
        head_to_scores: dict[str, dict[str, float]] = {}
        for h, g in sub.groupby("head"):
            head_to_scores[h] = dict(zip(g["compound_name"], g["score"]))

        all_compounds = sorted(set().union(*[d.keys() for d in head_to_scores.values()]))
        n_c = len(all_compounds)

        # Per-head rank vector aligned by compound_name
        rank_matrix = np.full((n_c, len(heads)), np.nan)
        for hi, h in enumerate(heads):
            scores_dict = head_to_scores.get(h, {})
            scores = np.array([scores_dict.get(c, np.nan) for c in all_compounds])
            rank_matrix[:, hi] = _normalise_to_rank(scores)

        # Pairwise Kendall-τ per target
        for i, hi in combinations(range(len(heads)), 2):
            ri = rank_matrix[:, i]
            rj = rank_matrix[:, hi]
            mask = ~(np.isnan(ri) | np.isnan(rj))
            if mask.sum() < 3:
                continue
            tau, _ = kendalltau(ri[mask], rj[mask])
            pairwise_rows.append({
                "target_uniprot": target_uni,
                "head_a": heads[i], "head_b": heads[hi],
                "kendall_tau": float(tau) if tau is not None else float("nan"),
                "n_compounds": int(mask.sum()),
            })

        # Per-compound disagreement entropy + max pairwise rank delta
        for ci, cname in enumerate(all_compounds):
            row = rank_matrix[ci]
            ent = _rank_entropy(row, n_c)
            valid = row[~np.isnan(row)]
            max_delta = int(valid.max() - valid.min()) if len(valid) >= 2 else 0
            # Facet tag
            if ent >= HIGH_DISAGREEMENT_ENTROPY and max_delta >= NOVEL_SCAFFOLD_DELTA:
                tag = "high_information_value"
            elif max_delta >= NOVEL_SCAFFOLD_DELTA:
                tag = "moderate_disagreement"
            else:
                tag = "agree"
            per_compound_rows.append({
                "target_uniprot": target_uni,
                "compound_name": cname,
                "disagreement_entropy": ent,
                "max_rank_delta": max_delta,
                "n_heads": int((~np.isnan(row)).sum()),
                "facet_tag": tag,
            })

    pairwise_df = pd.DataFrame(pairwise_rows)
    per_compound_df = pd.DataFrame(per_compound_rows)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    per_compound_df.to_parquet(args.out, index=False)
    logger.info("Wrote %s (%d rows)", args.out, len(per_compound_df))

    # Markdown report
    L: list[str] = []
    L.append("# Multi-Head Disagreement Axis v1 (V6.A.5)")
    L.append("")
    L.append(f"Pairwise Kendall-τ rank correlations + per-compound disagreement "
             f"entropy across {len(heads)} heads: **{', '.join(heads)}**. "
             "Generalises §8.15 (which was MAMMAL-vs-Tanimoto only) to N heads.")
    L.append("")
    L.append("## Pairwise Kendall-τ (target-by-target)")
    L.append("")
    L.append("τ = +1 means heads rank-agree perfectly; τ = -1 means they "
             "rank-anti-agree; τ ≈ 0 means independent rankings (the "
             "interesting case for ensemble information).")
    L.append("")
    L.append("| Target | MAMMAL_cal ↔ Tanimoto | MAMMAL_cal ↔ PrimeKG_PPR | Tanimoto ↔ PrimeKG_PPR |")
    L.append("|---|---|---|---|")
    targets_seen = pairwise_df["target_uniprot"].unique()
    for t in sorted(targets_seen):
        pivot = pairwise_df[pairwise_df.target_uniprot == t]
        get_pair = lambda a, b: pivot[
            (((pivot.head_a == a) & (pivot.head_b == b))
             | ((pivot.head_a == b) & (pivot.head_b == a)))
        ]["kendall_tau"].iloc[0] if len(pivot[
            (((pivot.head_a == a) & (pivot.head_b == b))
             | ((pivot.head_a == b) & (pivot.head_b == a)))
        ]) else float("nan")
        m_t = get_pair("MAMMAL_cal", "Tanimoto")
        m_k = get_pair("MAMMAL_cal", "PrimeKG_PPR")
        t_k = get_pair("Tanimoto", "PrimeKG_PPR")
        def fmt(x):
            return f"{x:+.2f}" if not pd.isna(x) else "—"
        L.append(f"| {t} | {fmt(m_t)} | {fmt(m_k)} | {fmt(t_k)} |")
    L.append("")

    # Aggregate τ stats
    L.append("## Pairwise τ summary across all targets")
    L.append("")
    by_pair = pairwise_df.groupby(["head_a", "head_b"])["kendall_tau"].describe().round(3)
    L.append("```")
    L.append(by_pair.to_string())
    L.append("```")
    L.append("")

    L.append("## Per-compound facet-tag distribution")
    L.append("")
    facet_counts = per_compound_df["facet_tag"].value_counts().to_dict()
    L.append("| Tag | Count | % |")
    L.append("|---|---|---|")
    n_total = max(len(per_compound_df), 1)
    for tag, n in sorted(facet_counts.items(), key=lambda kv: -kv[1]):
        L.append(f"| {tag} | {n} | {100 * n / n_total:.1f}% |")
    L.append("")

    # Top high-information-value compounds
    L.append("## Top 30 high-information-value (compound, target) pairs")
    L.append("")
    L.append("These are compounds where heads disagree maximally — exactly the "
             "wet-lab-priority candidates that single-head ranking would miss.")
    L.append("")
    hi = per_compound_df[per_compound_df.facet_tag == "high_information_value"]\
        .sort_values("max_rank_delta", ascending=False).head(30)
    L.append("| Compound | Target | Entropy | Rank Δ | N heads |")
    L.append("|---|---|---|---|---|")
    for _, r in hi.iterrows():
        L.append(f"| {r['compound_name']} | {r['target_uniprot']} | "
                 f"{r['disagreement_entropy']:.2f} | "
                 f"{r['max_rank_delta']} | {r['n_heads']} |")
    L.append("")

    L.append("## When V6.A.1 activates (MMAtt-DTA + PSICHIC + BALM)")
    L.append("")
    L.append("This script re-runs unchanged with 5+ heads. The interesting "
             "axis becomes the cross-head consensus on novel-scaffold "
             "compounds: where 4 of 5 heads agree but Tanimoto dissents, "
             "that's a novel-scaffold discovery candidate.")
    L.append("")
    L.append("---")
    L.append("")
    L.append("Generated by `scripts/51_v6_multihead_disagreement.py`.")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(L), encoding="utf-8")
    logger.info("Wrote %s", args.report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
