"""§8.0a — Render a Pareto-restructured wet-lab shortlist.

Composes:
  - v7 RRF ranking                          (efficacy axis)
  - combined_gates (liability counts)        (safety axis)
  - selectivity_v6_tanimoto_4metrics PI      (selectivity axis)
  - clinical_trials_v1 ip_status             (IP-freedom axis)
  - nootropic_similarity_v1 novelty_tag      (scaffold-novelty axis)

Outputs:
  data/results/v2/pareto_ranking_v1.parquet
  reports/pareto_ranking_v1.md

The Pareto front (rank 0) is the actionable wet-lab set: every compound on
the front is non-dominated across all 5 axes. Rank 1+ are progressively less
optimal frontiers — useful for triage when you need to expand the candidate
pool.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.fusion.pareto import (  # noqa: E402
    ParetoConfig, hypervolume_mc, rank_pareto,
)

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("v5_pareto")

DEFAULT_RANKING = ROOT / "data" / "results" / "v2" / "final_ranking_v7_moa.parquet"
DEFAULT_GATES = ROOT / "data" / "results" / "v2" / "combined_gates.parquet"
DEFAULT_SELECTIVITY = ROOT / "data" / "results" / "v2" / "selectivity_scores_tanimoto_v6_metrics.parquet"
DEFAULT_TRIALS = ROOT / "data" / "results" / "v2" / "clinical_trials_v1.parquet"
DEFAULT_NOOTROPIC = ROOT / "data" / "results" / "v2" / "nootropic_similarity_v1.parquet"
DEFAULT_OUT = ROOT / "data" / "results" / "v2" / "pareto_ranking_v1.parquet"
DEFAULT_REPORT = ROOT / "reports" / "pareto_ranking_v1.md"


def _safe_merge(left: pd.DataFrame, right: pd.DataFrame, on_col: str,
                 cols_to_take: list[str]) -> pd.DataFrame:
    if right is None or right.empty:
        for c in cols_to_take:
            if c not in left.columns:
                left[c] = pd.NA
        return left
    keep = [on_col] + [c for c in cols_to_take if c in right.columns]
    sub = right[keep].copy()
    sub["_jk"] = sub[on_col].str.lower().str.strip()
    left["_jk"] = left["compound_name"].str.lower().str.strip()
    out = left.merge(sub.drop(columns=[on_col]), on="_jk", how="left")
    return out.drop(columns=["_jk"])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ranking", type=Path, default=DEFAULT_RANKING)
    parser.add_argument("--gates", type=Path, default=DEFAULT_GATES)
    parser.add_argument("--selectivity", type=Path, default=DEFAULT_SELECTIVITY)
    parser.add_argument("--trials", type=Path, default=DEFAULT_TRIALS)
    parser.add_argument("--nootropic", type=Path, default=DEFAULT_NOOTROPIC)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--pass-only", action="store_true",
                        help="Restrict to final_status==PASS before Pareto sort.")
    parser.add_argument("--top-n", type=int, default=30,
                        help="Show top N from front in the report.")
    args = parser.parse_args()

    rk = pd.read_parquet(args.ranking)
    gates = pd.read_parquet(args.gates) if args.gates.exists() else None
    sel = pd.read_parquet(args.selectivity) if args.selectivity.exists() else None
    trials = pd.read_parquet(args.trials) if args.trials.exists() else None
    noot = pd.read_parquet(args.nootropic) if args.nootropic.exists() else None
    logger.info("Loaded: ranking=%d, gates=%s, sel=%s, trials=%s, noot=%s",
                len(rk),
                len(gates) if gates is not None else None,
                len(sel) if sel is not None else None,
                len(trials) if trials is not None else None,
                len(noot) if noot is not None else None)

    df = rk.copy()
    df = _safe_merge(df, gates, "compound_name",
                     ["final_status", "n_tier_1", "n_tier_2", "n_tier_3"])
    df = _safe_merge(df, sel, "compound_name",
                     ["partition_index_top", "gini", "entropy", "s_10x",
                      "selectivity_category"])
    df = _safe_merge(df, trials, "compound_name",
                     ["ip_status", "n_trials", "latest_phase"])
    df = _safe_merge(df, noot, "compound_name",
                     ["nearest_nootropic", "nearest_nootropic_tanimoto",
                      "nootropic_novelty_tag"])

    if args.pass_only and "final_status" in df.columns:
        before = len(df)
        df = df[df["final_status"] == "PASS"].reset_index(drop=True)
        logger.info("PASS-only filter: %d → %d", before, len(df))

    # Fill missing axis inputs
    for c, fill in [
        ("n_tier_1", 0), ("n_tier_2", 0), ("partition_index_top", 0.0),
        ("ip_status", ""), ("nootropic_novelty_tag", "unknown"),
    ]:
        if c not in df.columns:
            df[c] = fill
        else:
            df[c] = df[c].fillna(fill)

    df = rank_pareto(df, ParetoConfig())
    logger.info("Pareto sorted: %d frontiers (rank 0..%d)",
                df["pareto_rank"].nunique(), int(df["pareto_rank"].max()))

    # Hypervolume of the front (rank 0)
    A = df[[f"_axis_{a}" for a in
            ("efficacy_rrf", "safety_neg_liability", "selectivity_pi",
             "ip_novelty", "scaffold_novelty")]].to_numpy()
    front_idx = df.index[df["pareto_rank"] == 0].to_numpy()
    front = A[front_idx]
    if len(front):
        ref_point = A.min(axis=0)
        hv = hypervolume_mc(front, ref_point, n_samples=200_000)
    else:
        hv = 0.0

    # Save parquet
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(args.out, index=False)
    logger.info("Wrote %s", args.out)

    # Report
    L: list[str] = []
    L.append("# Pareto-Restructured Shortlist v1 (§8.0a)")
    L.append("")
    L.append("Five-axis non-dominated sort over the V7 wet-lab handoff. Every "
             "compound on **frontier rank 0** is non-dominated — there is no "
             "compound that beats it on every axis.")
    L.append("")
    L.append("**Axes** (all directions = HIGHER is better):")
    L.append("- `efficacy_rrf` — v7 RRF score from 5-cluster fusion (MAMMAL+Tanimoto+Boltzina+ADMET+MoA)")
    L.append("- `safety_neg_liability` — −(n_tier_1 + 0.5×n_tier_2) from §8.0b-zn liability panel")
    L.append("- `selectivity_pi` — top-target Cheng 2010 Partition Index from §7.4 v2")
    L.append("- `ip_novelty` — CTgov §8.3 IP-status mapped to [0,1]: approved=0.0, investigational=0.4, early=0.7, none=1.0")
    L.append("- `scaffold_novelty` — §8.10 nootropic-similarity tag: analog=0.0, intermediate=0.5, novel_scaffold=1.0")
    L.append("")
    L.append(f"**Pareto front (rank 0)**: {int((df['pareto_rank']==0).sum())} compounds")
    L.append(f"**Hypervolume (MC, 200k samples)**: {hv:.4f}")
    L.append(f"**Frontier count**: {int(df['pareto_rank'].nunique())} (max rank = {int(df['pareto_rank'].max())})")
    L.append("")

    L.append(f"## Pareto front (rank 0) — top {args.top_n} by RRF")
    L.append("")
    front_df = df[df["pareto_rank"] == 0].sort_values(
        "rrf_score", ascending=False).head(args.top_n)
    L.append("| # | Compound | Tier | RRF | Safety | PI | IP | Novelty |")
    L.append("|---|---|---|---|---|---|---|---|")
    for i, r in front_df.iterrows():
        L.append(f"| {i+1} | {r['compound_name']} | "
                 f"{r.get('evidence_tier','?')} | "
                 f"{r['rrf_score']:.3f} | "
                 f"{r['_axis_safety_neg_liability']:+.1f} | "
                 f"{r.get('partition_index_top', 0):.2f} | "
                 f"{r.get('ip_status','')} | "
                 f"`{r.get('nootropic_novelty_tag','')}` |")
    L.append("")

    L.append("## Rank-1 frontier (one peel below Pareto) — top 15")
    L.append("")
    r1 = df[df["pareto_rank"] == 1].sort_values("rrf_score", ascending=False).head(15)
    if len(r1):
        L.append("| # | Compound | Tier | RRF | Safety | PI | IP | Novelty |")
        L.append("|---|---|---|---|---|---|---|---|")
        for i, r in r1.iterrows():
            L.append(f"| {i+1} | {r['compound_name']} | "
                     f"{r.get('evidence_tier','?')} | "
                     f"{r['rrf_score']:.3f} | "
                     f"{r['_axis_safety_neg_liability']:+.1f} | "
                     f"{r.get('partition_index_top', 0):.2f} | "
                     f"{r.get('ip_status','')} | "
                     f"`{r.get('nootropic_novelty_tag','')}` |")
    L.append("")

    L.append("## Frontier-size distribution")
    L.append("")
    counts = df["pareto_rank"].value_counts().sort_index()
    L.append("| Rank | N compounds |")
    L.append("|---|---|")
    for r, n in counts.items():
        L.append(f"| {int(r)} | {int(n)} |")
    L.append("")

    L.append("## Interpretation")
    L.append("")
    L.append("- **Front (rank 0)** is the *production* wet-lab handoff for "
             "constrained-budget runs. Every compound here strictly dominates "
             "nothing on the front.")
    L.append("- **Crowding distance** (saved as `crowding_distance` column in "
             "the parquet) within rank 0 ranks compounds by frontier diversity "
             "— pick the top-k by crowding when the front is too large.")
    L.append("- **Rank 1-2** are useful when expanding the candidate pool: "
             "these are compounds beaten on exactly one or two axes.")
    L.append("- **Hypervolume** measures the volume of objective space "
             "dominated by the front; higher = better coverage. Use as a "
             "regression metric when re-running with new evidence.")
    L.append("")
    L.append("---")
    L.append("")
    L.append(f"Generated by `scripts/42_v5_pareto_shortlist.py`. "
             f"Input ranking: `{args.ranking.name}`. PASS-only={args.pass_only}.")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(L), encoding="utf-8")
    logger.info("Wrote %s", args.report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
