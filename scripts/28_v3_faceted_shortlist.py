"""V3 §8.1 — Multi-class faceted wet-lab shortlist.

8 mechanism-class facets + 9 targeted-pair facets, top-5 each, with cross-
facet provenance to prevent triple-counting. Reads:
  - data/results/v2/selectivity_scores.parquet (from scripts/27)
  - data/results/v2/final_ranking_calibrated_v4_tanimoto.parquet (from 15)
  - data/results/dti_scores.parquet (the MAMMAL grid)

Validates gates G3-G6:
  G3 (faceted CHRNA7): top-5 must contain TC-5619 AND encenicline
  G4 (faceted ACHE):   top-5 must contain donepezil AND galantamine
  G5 (faceted HRH3):   top-5 must contain pitolisant
  G6 (cross-facet):    pitolisant in HRH3 facet but NOT in HRH3+DRD1 pair

Output: reports/wet_lab_shortlist_v4_faceted.md +
        data/results/v2/faceted_shortlist.parquet
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

from mammal_repurposing.config import DTI_SCORES_PARQUET, RESULTS_DIR, TARGETS_PARQUET  # noqa: E402
from mammal_repurposing.fusion.faceted_shortlist import (  # noqa: E402
    TARGETED_PAIRS,
    build_by_class_facets,
    build_targeted_pair_facets,
    compute_cross_facet_provenance,
)
from mammal_repurposing.selectivity.gini_scorecard import MECHANISM_CLASS_TARGETS  # noqa: E402

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("v3_faceted")

V2_DIR = RESULTS_DIR / "v2"
DIRECTION_NOTES = {
    "orexinergic": "⚠ FDA-approved orexin drugs (suvorexant, lemborexant) are "
                   "antagonists for sleep — opposite of procognitive direction. "
                   "Surfaced for review but tagged WRONG_DIRECTION_FOR_COGNITION.",
}

GATE_G3 = {"cholinergic": ["tc-5619", "encenicline"]}      # CHRNA7 + ACHE class
GATE_G4 = {"cholinergic": ["donepezil", "galantamine"]}    # same class, AChE inhibitors
GATE_G5 = {"histaminergic": ["pitolisant"]}                # HRH3 class
GATE_G6_HRH3DRD1_EXCLUDE = "pitolisant"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selectivity", type=Path,
                        default=V2_DIR / "selectivity_scores.parquet")
    parser.add_argument("--final-ranking", type=Path,
                        default=V2_DIR / "final_ranking_calibrated_v4_tanimoto.parquet")
    parser.add_argument("--dti", type=Path, default=DTI_SCORES_PARQUET)
    parser.add_argument("--targets", type=Path, default=TARGETS_PARQUET)
    parser.add_argument("--out", type=Path,
                        default=V2_DIR / "faceted_shortlist.parquet")
    parser.add_argument("--report", type=Path,
                        default=ROOT / "reports" / "wet_lab_shortlist_v4_faceted.md")
    parser.add_argument("--top-n", type=int, default=5)
    args = parser.parse_args()

    sel_df = pd.read_parquet(args.selectivity)
    dti = pd.read_parquet(args.dti)
    targets = pd.read_parquet(args.targets)
    uniprot_to_gene = dict(zip(targets["uniprot"], targets["gene"]))
    if "target_gene" not in dti.columns:
        dti = dti.assign(target_gene=dti["target_uniprot"].map(uniprot_to_gene))

    # Pull RRF score (efficacy proxy) from the v4 calibrated ranking if available
    if args.final_ranking.exists():
        final_df = pd.read_parquet(args.final_ranking)
        rrf = final_df[["compound_name", "rrf_score"]].drop_duplicates("compound_name")
        sel_df = sel_df.merge(rrf, on="compound_name", how="left")
        logger.info("Joined RRF scores from %s (%d compounds matched).",
                    args.final_ranking, sel_df["rrf_score"].notna().sum())
    else:
        logger.warning("No final ranking parquet at %s — using top_target_pkd as efficacy.",
                       args.final_ranking)
        sel_df["rrf_score"] = sel_df["top_target_pkd"]

    composite_col = "rrf_score"

    # --- Build facets ----------------------------------------------------------
    logger.info("Building 8 mechanism-class facets (top-%d each) ...", args.top_n)
    by_class = build_by_class_facets(
        sel_df, dti, top_n=args.top_n, direction_notes=DIRECTION_NOTES,
        composite_col=composite_col,
    )
    logger.info("  → %d entries.", len(by_class))

    logger.info("Building 9 targeted-pair facets (top-%d each) ...", args.top_n)
    by_pair = build_targeted_pair_facets(sel_df, dti, top_n=args.top_n)
    logger.info("  → %d entries.", len(by_pair))

    faceted = pd.concat([by_class, by_pair], ignore_index=True)
    cross_prov = compute_cross_facet_provenance(faceted)
    faceted_with_prov = faceted.merge(
        cross_prov.rename(columns={"cross_facet_list": "_cross_facet_list",
                                    "cross_facet_count": "_cross_facet_count"}),
        on="compound_name", how="left",
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    faceted_with_prov.to_parquet(args.out, index=False)
    logger.info("Wrote %s.", args.out)

    # --- Validation gates G3–G6 ------------------------------------------------
    gate_results = []
    for gate_id, gate_dict in [("G3", GATE_G3), ("G4", GATE_G4), ("G5", GATE_G5)]:
        for facet_name, must_contain in gate_dict.items():
            sub = by_class[by_class["facet_name"] == facet_name]
            top_names = [n.lower() for n in sub["compound_name"]]
            present = [n for n in must_contain if n in top_names]
            passed = len(present) == len(must_contain)
            gate_results.append({
                "gate": gate_id,
                "facet": facet_name,
                "must_contain": ", ".join(must_contain),
                "present": ", ".join(present),
                "passed": passed,
            })
    # G6: pitolisant in HRH3 facet but NOT in HRH3+DRD1 pair
    hrh3_top = [n.lower() for n in
                by_class[by_class["facet_name"] == "histaminergic"]["compound_name"]]
    pair_top = [n.lower() for n in
                by_pair[by_pair["facet_name"] == "HRH3+DRD1"]["compound_name"]]
    pit_in_hrh3 = GATE_G6_HRH3DRD1_EXCLUDE in hrh3_top
    pit_in_pair = GATE_G6_HRH3DRD1_EXCLUDE in pair_top
    g6_passed = pit_in_hrh3 and not pit_in_pair
    gate_results.append({
        "gate": "G6", "facet": "HRH3 / HRH3+DRD1 hygiene",
        "must_contain": "pitolisant in HRH3 facet AND not in HRH3+DRD1 pair",
        "present": (f"HRH3={pit_in_hrh3}, HRH3+DRD1={pit_in_pair}"),
        "passed": g6_passed,
    })
    gate_df = pd.DataFrame(gate_results)

    # --- Render markdown -------------------------------------------------------
    L: list[str] = []
    L.append("# Wet-Lab Shortlist v4 — Multi-Class Faceted (§8.1)")
    L.append("")
    L.append("Per research/4-tier/Graczyk-Style ... .md §2. Top-5 per facet across "
             "8 mechanism classes + 9 targeted pairs, with cross-facet provenance.")
    L.append("")
    L.append("This dissolves the v3 HRH3-23/25 lock-in into a structured, "
             "mechanism-orthogonal shortlist medicinal chemists can triage.")
    L.append("")

    # Gates
    L.append("## Validation gates (G3–G6)")
    L.append("")
    L.append("| Gate | Facet | Must contain | Present | Passed |")
    L.append("|---|---|---|---|---|")
    for _, r in gate_df.iterrows():
        status = "✅" if r["passed"] else "❌"
        L.append(f"| {r['gate']} | {r['facet']} | {r['must_contain']} | "
                 f"{r['present']} | {status} |")
    L.append("")

    # Mechanism-class facets
    L.append("## Mechanism-class facets (top-5 each)")
    L.append("")
    for class_name in MECHANISM_CLASS_TARGETS.keys():
        sub = by_class[by_class["facet_name"] == class_name]
        if sub.empty:
            continue
        L.append(f"### {class_name}")
        if class_name in DIRECTION_NOTES:
            L.append(f"_{DIRECTION_NOTES[class_name]}_")
        L.append("")
        L.append("| # | Compound | Composite | Gini | S(10x) | Category | Top target | "
                 "Cross-facet count |")
        L.append("|---|---|---|---|---|---|---|---|")
        for _, r in sub.iterrows():
            cf_count = int(faceted_with_prov.loc[
                faceted_with_prov["compound_name"] == r["compound_name"],
                "_cross_facet_count"].iloc[0])
            L.append(f"| {int(r['facet_rank'])} | {r['compound_name']} | "
                     f"{r['composite_score']:.3f} | {r['gini']:.2f} | "
                     f"{int(r['s_10x'])} | `{r['selectivity_category']}` | "
                     f"{r['top_target']} | {cf_count} |")
        L.append("")

    # Targeted-pair facets
    L.append("## Targeted-pair facets (top-5 each)")
    L.append("")
    for pair_name, (genes, hypothesis) in TARGETED_PAIRS.items():
        sub = by_pair[by_pair["facet_name"] == pair_name]
        if sub.empty:
            continue
        L.append(f"### {pair_name} ({', '.join(genes)})")
        L.append(f"_{hypothesis}_")
        L.append("")
        if len(sub) == 0:
            L.append("_(no compounds — facet empty; emptiness is the finding)_")
            L.append("")
            continue
        L.append("| # | Compound | Composite | Gini | Category | Top target | Cross-facet |")
        L.append("|---|---|---|---|---|---|---|")
        for _, r in sub.iterrows():
            cf_count = int(faceted_with_prov.loc[
                faceted_with_prov["compound_name"] == r["compound_name"],
                "_cross_facet_count"].iloc[0])
            L.append(f"| {int(r['facet_rank'])} | {r['compound_name']} | "
                     f"{r['composite_score']:.3f} | {r['gini']:.2f} | "
                     f"`{r['selectivity_category']}` | {r['top_target']} | {cf_count} |")
        L.append("")

    # Cross-facet champions — compounds appearing in 3+ facets
    L.append("## Cross-facet champions (≥3 facets)")
    L.append("")
    champs = cross_prov[cross_prov["cross_facet_count"] >= 3]
    L.append("| Compound | # facets | Facets |")
    L.append("|---|---|---|")
    for _, r in champs.iterrows():
        L.append(f"| {r['compound_name']} | {int(r['cross_facet_count'])} | "
                 f"{r['cross_facet_list']} |")
    L.append("")

    L.append("---")
    L.append("")
    L.append("Generated by `scripts/28_v3_faceted_shortlist.py`.")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(L), encoding="utf-8")
    logger.info("Wrote %s.", args.report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
