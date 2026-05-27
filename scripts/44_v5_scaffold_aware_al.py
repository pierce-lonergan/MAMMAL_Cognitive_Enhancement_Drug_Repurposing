"""§7.13 — Run scaffold-aware AL re-ranking + diversity hypothesis test.

Composes v7 RRF ranking with the scaffold-aware AL re-ranker. Hypothesis:
the AL top-25 contains more distinct Murcko scaffolds than the baseline RRF
top-25. PASS if delta ≥ 2 distinct scaffolds.

Outputs:
    data/results/v2/scaffold_aware_v1.parquet
    reports/scaffold_aware_v1.md
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.diagnostics.scaffold_aware_al import (  # noqa: E402
    ScaffoldAwareConfig, evaluate_diversity, rank_with_scaffold_bonus,
)

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("v5_scaffold_al")

DEFAULT_RANKING = ROOT / "data" / "results" / "v2" / "final_ranking_v7_moa.parquet"
DEFAULT_OUT = ROOT / "data" / "results" / "v2" / "scaffold_aware_v1.parquet"
DEFAULT_REPORT = ROOT / "reports" / "scaffold_aware_v1.md"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ranking", type=Path, default=DEFAULT_RANKING)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--alpha", type=float, default=0.7,
                        help="Weight on efficacy vs exploration bonus.")
    parser.add_argument("--top-k", type=int, default=25)
    parser.add_argument("--pass-only", action="store_true",
                        help="Restrict to final_status==PASS first.")
    args = parser.parse_args()

    rk = pd.read_parquet(args.ranking)
    if args.pass_only:
        gates = pd.read_parquet(ROOT / "data" / "results" / "v2" / "combined_gates.parquet")
        keep = set(gates[gates["final_status"] == "PASS"]["compound_name"].str.lower().str.strip())
        before = len(rk)
        rk = rk[rk["compound_name"].str.lower().str.strip().isin(keep)].reset_index(drop=True)
        logger.info("PASS-only filter: %d → %d", before, len(rk))

    if "compound_smiles" not in rk.columns or rk["compound_smiles"].isna().all():
        compounds = pd.read_parquet(ROOT / "data" / "interim" / "compounds.parquet")
        smi_map = dict(zip(compounds["name"].str.lower().str.strip(),
                           compounds["smiles"]))
        rk = rk.copy()
        rk["compound_smiles"] = rk["compound_name"].str.lower().str.strip().map(smi_map)

    cfg = ScaffoldAwareConfig(alpha=args.alpha)
    reranked = rank_with_scaffold_bonus(rk, cfg)
    reranked.to_parquet(args.out, index=False)
    logger.info("Wrote %s", args.out)

    # Diversity hypothesis test
    baseline_top = rk.sort_values("rrf_score", ascending=False).head(args.top_k)
    al_top = reranked.head(args.top_k)
    div = evaluate_diversity(baseline_top, al_top)
    hypothesis_pass = div["delta"] >= 2

    # Report
    L: list[str] = []
    L.append(f"# Scaffold-Aware Active Learning v1 (§7.13) — α = {args.alpha}")
    L.append("")
    L.append("For every compound, computes the Bemis-Murcko scaffold and "
             "applies an exploration bonus inversely proportional to "
             "scaffold density in the library. The final AL score is "
             f"α · normalized_RRF + (1−α) · scaffold_bonus with α={args.alpha}.")
    L.append("")
    L.append("**Hypothesis**: AL top-k contains more distinct Murcko scaffolds "
             "than baseline RRF top-k (≥2 more).")
    L.append("")
    L.append(f"**Verdict**: {'PASS' if hypothesis_pass else 'DEGRADE'} "
             f"— baseline = {div['baseline_n_distinct_scaffolds']} distinct scaffolds; "
             f"AL = {div['al_n_distinct_scaffolds']} distinct scaffolds; "
             f"Δ = +{div['delta']}; AL-only scaffolds = {div['al_only_scaffolds']}.")
    L.append("")

    L.append(f"## Baseline top-{args.top_k} (by RRF)")
    L.append("")
    L.append("| # | Compound | Tier | RRF | Scaffold density |")
    L.append("|---|---|---|---|---|")
    baseline_with_scaf = rank_with_scaffold_bonus(baseline_top.copy(), cfg)
    for _, r in baseline_with_scaf.iterrows():
        L.append(f"| {int(r['al_rank'])} | {r['compound_name']} | "
                 f"{r.get('evidence_tier','?')} | {r['rrf_score']:.3f} | "
                 f"{int(r['scaffold_density'])} |")
    L.append("")

    L.append(f"## AL-reranked top-{args.top_k}")
    L.append("")
    L.append("| AL rank | Compound | Tier | RRF | Norm RRF | Scaf density | AL score |")
    L.append("|---|---|---|---|---|---|---|")
    for _, r in al_top.iterrows():
        L.append(f"| {int(r['al_rank'])} | {r['compound_name']} | "
                 f"{r.get('evidence_tier','?')} | "
                 f"{r['rrf_score']:.3f} | "
                 f"{r['normalized_rrf']:.2f} | "
                 f"{int(r['scaffold_density'])} | "
                 f"{r['al_score']:.3f} |")
    L.append("")

    # AL-only compounds: in AL top-k but NOT in baseline top-k
    baseline_names = set(baseline_top["compound_name"].str.lower())
    al_only = al_top[~al_top["compound_name"].str.lower().isin(baseline_names)]
    L.append(f"## AL surfaced new compounds (n={len(al_only)})")
    L.append("")
    L.append("These compounds entered the top-k via the scaffold bonus — "
             "their RRF ranks were below baseline but they sit in "
             "undersampled scaffold buckets.")
    L.append("")
    if len(al_only):
        L.append("| Compound | RRF | Scaffold density | AL score |")
        L.append("|---|---|---|---|")
        for _, r in al_only.iterrows():
            L.append(f"| {r['compound_name']} | {r['rrf_score']:.3f} | "
                     f"{int(r['scaffold_density'])} | {r['al_score']:.3f} |")
    L.append("")

    L.append("---")
    L.append("")
    L.append(f"Generated by `scripts/44_v5_scaffold_aware_al.py` "
             f"(α={args.alpha}, top_k={args.top_k}).")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(L), encoding="utf-8")
    logger.info("Wrote %s. Hypothesis %s: baseline=%d distinct vs AL=%d distinct, Δ=+%d",
                args.report,
                "PASS" if hypothesis_pass else "DEGRADE",
                div['baseline_n_distinct_scaffolds'],
                div['al_n_distinct_scaffolds'],
                div['delta'])
    return 0 if hypothesis_pass else 2


if __name__ == "__main__":
    raise SystemExit(main())
