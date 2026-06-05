"""V6.B.5 — Generate the ~210-target GWAS-anchored cognition panel.

Stub-mode driver consuming src/mammal_repurposing/cluster_d/panel_expansion.py
build_expanded_panel(). Validates that:
  - 22-target V6.B panel is a strict subset
  - All UniProts unique
  - Each target has ≥1 inclusion rule (l2g / magma_p / ahba_cortical /
    sc_zscore / lit_otar / v6b_panel_22_anchor / v5_liability_panel_44)

Outputs:
  data/results/v2/panel_expanded_v1.parquet
  reports/pipeline/panel_expansion_v1.md

Real-mode (V6.B.5 Stage 2-3) replaces the hand-curated expansion list with
live OT Genetics L2G + cellxgene-census + Moodie 2024 g-cortical alignment
+ Lit-OTAR (Kafkas 2024) query results.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("v6b5_panel_expand")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path,
                        default=ROOT / "data" / "results" / "v2"
                        / "panel_expanded_v1.parquet")
    parser.add_argument("--report", type=Path,
                        default=ROOT / "reports" / "pipeline" / "panel_expansion_v1.md")
    args = parser.parse_args()

    from mammal_repurposing.cluster_d.panel_expansion import (
        build_expanded_panel, validate_panel, availability,
    )

    df = build_expanded_panel()
    val = validate_panel(df)
    avail = availability()

    logger.info("Expanded panel: %d total targets (22-panel subset: %d/%d ok=%s; "
                "liability subset: %d; unique uniprots: %d ok=%s)",
                val["n_total"], val["n_in_v6b_panel_22"],
                val["v6b_panel_22_expected"], val["v6b_panel_22_subset_ok"],
                val["n_in_liability_panel_44"], val["n_unique_uniprot"],
                val["dedup_ok"])

    # Persist
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(args.out, index=False)
    logger.info("Wrote %s (%d rows)", args.out, len(df))

    # Report
    L: list[str] = []
    L.append("# V6.B.5 Panel Expansion v1 (Stage 1, stub-mode)")
    L.append("")
    L.append("Hand-curated approximation of the GWAS+AHBA+SC+Lit-OTAR live "
             "query result per Cluster D §F. Real-mode V6.B.5 Stage 2-3 "
             "replaces this with live OT Genetics + cellxgene + Moodie 2024 "
             "+ Lit-OTAR queries.")
    L.append("")
    L.append("## Summary")
    L.append("")
    L.append(f"- **Total panel targets**: {val['n_total']}")
    L.append(f"- **22-target V6.B subset**: {val['n_in_v6b_panel_22']} / 22 "
             f"({'✅' if val['v6b_panel_22_subset_ok'] else '❌'})")
    L.append(f"- **Liability panel members**: {val['n_in_liability_panel_44']}")
    L.append(f"- **Unique UniProts**: {val['n_unique_uniprot']} "
             f"({'✅' if val['dedup_ok'] else '❌'})")
    L.append("")
    L.append("## Inclusion rules fired")
    L.append("")
    L.append("| Rule | Count |")
    L.append("|---|---|")
    for rule, n in sorted(val["rules_distribution"].items(),
                          key=lambda x: -x[1]):
        L.append(f"| {rule} | {n} |")
    L.append("")
    L.append("## Sample top-30 targets")
    L.append("")
    L.append("| UniProt | Gene | Rules | In V6.B-22 | In Liability-44 |")
    L.append("|---|---|---|---|---|")
    for _, r in df.head(30).iterrows():
        v22 = "✅" if r["in_v6b_panel_22"] else "—"
        v44 = "✅" if r["in_liability_panel_44"] else "—"
        rules_short = r["rules_fired"][:50] + ("…" if len(r["rules_fired"]) > 50 else "")
        L.append(f"| {r['uniprot']} | {r['gene_symbol']} | {rules_short} | "
                 f"{v22} | {v44} |")
    L.append("")
    L.append("## Honest caveats")
    L.append("")
    L.append("- This is the V6.B.5 Stage 1 STUB. The expansion list is "
             "hand-curated to approximate what a live GWAS L2G + MAGMA + "
             "AHBA (Moodie 2024 g-cortical) + cellxgene single-cell + "
             "Lit-OTAR (Kafkas 2024) query would return.")
    L.append("- Real-mode V6.B.5 Stage 2 requires the OT Genetics L2G live "
             "fetch to land (currently network-blocked in sandbox).")
    L.append("- Real-mode V6.B.5 Stage 3 requires cellxgene-census brain "
             "slice + Moodie 2024 g-cortical alignment.")
    L.append("- The ~210-target count is per Cluster D §F target estimate; "
             "actual real-mode count depends on threshold choices "
             "(L2G ≥ 0.2 vs 0.3, MAGMA α=0.05 vs 0.001, etc.).")
    L.append("- V6.B.3 PyMC NUTS run on this expanded panel is a V6.B.5 "
             "Stage 4 deliverable — gene-level T≈15,000 requires sparse "
             "approximation per V6.B plan.")
    L.append("")
    L.append("---")
    L.append("")
    L.append("Generated by `scripts/61_v6b5_panel_expand.py`. V6.B.5 Stage 1 "
             "stub-mode validation of the panel expansion architecture.")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(L), encoding="utf-8")
    logger.info("Wrote %s", args.report)

    # Exit: 0 if all sanity checks pass, 1 otherwise
    if not val["v6b_panel_22_subset_ok"] or not val["dedup_ok"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
