"""§8.13 — Pocket-class-conditioned liability gating demo.

**Status**: gate logic shipped; live operation awaits §7.17 pose-saving Boltz
wrapper. This script demonstrates the composition end-to-end using SYNTHETIC
pose classifications derived from the 7-target curated pocket DB (§7.5).

What it shows:
  1. Load the §8.0b liability gates (z-norm mode) — currently CUT=14 / FLAG=21 / PASS=80.
  2. Build a synthetic pose-classification table at the targets we have pocket
     DB coverage for (CHRNA7, ACHE, HRH3, DRD1, PDE4D, SIGMAR1, GRIN2B). For each
     CUT compound we test three scenarios: pose at orthosteric, pose at
     allosteric_known, pose at surface_artifact.
  3. Run pocket_aware_liability_gate() and report how many CUTs get demoted
     to FLAG when their pose is allosteric_known vs orthosteric.

When §7.17 ships, replace the synthetic pose-class column with the real
classifier output from `pockets/pocket_classifier.py:classify_pose(pose_xyz, ...)`.
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

from mammal_repurposing.gates.liability_panel import (  # noqa: E402
    POCKET_AWARE_DEMOTABLE, pocket_aware_liability_gate,
)

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("v5_pocket_liability")

DEFAULT_GATES = ROOT / "data" / "results" / "v2" / "liability_gates.parquet"
DEFAULT_GATES_OUT = ROOT / "data" / "results" / "v2" / "liability_gates_pocket_aware.parquet"
DEFAULT_REPORT = ROOT / "reports" / "liability_pocket_aware_v1.md"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gates", type=Path, default=DEFAULT_GATES)
    parser.add_argument("--out", type=Path, default=DEFAULT_GATES_OUT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument(
        "--scenario", choices=["orthosteric", "allosteric_known", "surface_artifact"],
        default="allosteric_known",
        help="Synthetic pose scenario for the demo (default: allosteric_known, "
             "the case where demotion fires).",
    )
    args = parser.parse_args()

    gates = pd.read_parquet(args.gates)
    logger.info("Loaded liability gates: %d rows (CUT=%d, FLAG=%d, PASS=%d)",
                len(gates),
                int((gates["liability_status"] == "CUT").sum()),
                int((gates["liability_status"] == "FLAG").sum()),
                int((gates["liability_status"] == "PASS").sum()))

    # Build synthetic pose-classifications: every CUT compound is assigned the
    # `scenario` pocket_class at every gene in its Tier 1 hits that is in
    # POCKET_AWARE_DEMOTABLE.
    rows: list[dict] = []
    for _, r in gates.iterrows():
        if r["liability_status"] != "CUT":
            continue
        for gene in str(r.get("tier_1_hits", "") or "").split(";"):
            if not gene:
                continue
            if gene not in POCKET_AWARE_DEMOTABLE:
                continue
            rows.append({
                "compound_name": r["compound_name"],
                "target_gene": gene,
                "pocket_class": args.scenario,
            })
    pose_df = pd.DataFrame(rows)
    logger.info("Synthetic pose-classification table: %d rows (scenario=%s)",
                len(pose_df), args.scenario)

    augmented = pocket_aware_liability_gate(gates, pose_df)
    augmented.to_parquet(args.out, index=False)
    logger.info("Wrote %s", args.out)

    # Summarise the demotion effect
    before = augmented["liability_status"].value_counts().to_dict()
    after = augmented["liability_status_pocket_aware"].value_counts().to_dict()

    L: list[str] = []
    L.append("# Pocket-Conditioned Liability Gate v1 (§8.13)")
    L.append("")
    L.append(f"Demo run with **synthetic pose scenario = `{args.scenario}`**. "
             "Real-grid operation awaits §7.17 pose-saving Boltz wrapper.")
    L.append("")
    L.append("Per the V4 §8.13 design + research/4-tier/archived/analysis-notes/"
             "Pocket-Conditioned-Boltz2.md §3.3, the absolute-mode §8.0b CUT "
             "is too aggressive when the predicted pose binds OUTSIDE the "
             "orthosteric pocket. This gate applies literature-grounded "
             "demotion rules:")
    L.append("")
    L.append("- **5-HT2B**: Roth 2007 valvulopathy class warning applies to "
             "orthosteric agonists (fen-phen / pergolide / cabergoline pattern); "
             "allosteric NAMs don't trigger.")
    L.append("- **hERG**: Dumotier & Urban 2024 — central-pore Y652/F656/T623 "
             "binding is the classical block; allosteric or vestibule binding "
             "is materially lower risk.")
    L.append("- **HRH1**: Gray 2015 anticholinergic dementia risk is "
             "orthosteric-antagonist mediated; allosteric H1 ligands have no "
             "documented cognition-impairing precedent.")
    L.append("- **CB1**: Topol 2010 CRESCENDO rimonabant neuropsych AEs are "
             "orthosteric-only; allosteric NAMs have different clinical "
             "risk profile.")
    L.append("- **CHRM1 / OPRM1 / MAOA**: demoted by the same logic.")
    L.append("")
    L.append("## Demotable pocket classes per gene")
    L.append("")
    L.append("| Gene | Demotable pocket classes |")
    L.append("|---|---|")
    for g, classes in POCKET_AWARE_DEMOTABLE.items():
        L.append(f"| {g} | {', '.join(sorted(classes))} |")
    L.append("")
    L.append("## Verdict shift (synthetic demo)")
    L.append("")
    L.append("| Status | Before | After (`pocket_aware`) | Δ |")
    L.append("|---|---|---|---|")
    for s in ("CUT", "FLAG", "PASS"):
        b = before.get(s, 0)
        a = after.get(s, 0)
        L.append(f"| {s} | {b} | {a} | {a - b:+d} |")
    L.append("")
    L.append(f"_Compounds touched_: {int((augmented['n_pocket_demoted'] > 0).sum())} "
             f"with non-zero demotions.")
    L.append("")

    demoted = augmented[augmented["n_pocket_demoted"] > 0]
    if len(demoted):
        L.append("## Per-compound demotion detail")
        L.append("")
        L.append("| Compound | Original | Pocket-aware | Demotions |")
        L.append("|---|---|---|---|")
        for _, r in demoted.iterrows():
            L.append(f"| {r['compound_name']} | {r['liability_status']} | "
                     f"{r['liability_status_pocket_aware']} | "
                     f"{r['pocket_demotions']} |")
        L.append("")

    L.append("## How to operationalise")
    L.append("")
    L.append("This script's `pose_df` is **synthetic** — every CUT compound is "
             "assigned the scenario pocket-class. To run the gate against real "
             "Boltz-2 poses:")
    L.append("")
    L.append("```python")
    L.append("# 1. Save mmCIF poses during the Boltz sweep (§7.17 unblocks this).")
    L.append("# 2. Extract heavy-atom centroid xyz per pose.")
    L.append("# 3. For each (compound, target_gene) classify via §7.5:")
    L.append("from mammal_repurposing.pockets.pocket_classifier import classify_pose")
    L.append("from mammal_repurposing.pockets.pocket_database import load_pocket_database")
    L.append("db = load_pocket_database('data/pockets/centroids/')")
    L.append("pose_df = (pose_centroids_df.apply(lambda r: ")
    L.append("    classify_pose(r[['x','y','z']].values, r['target_gene'], db),")
    L.append("    axis=1))")
    L.append("# 4. Then call pocket_aware_liability_gate(liability_gates_df, pose_df).")
    L.append("```")
    L.append("")
    L.append("Currently only 7 of 22 cognition targets have curated centroids "
             "(CHRNA7, ACHE, HRH3, DRD1, PDE4D, SIGMAR1, GRIN2B). To extend to "
             "the full 44-target liability panel, either: (a) reuse a single "
             "canonical orthosteric pocket per liability target via UniProt → "
             "reference PDB lookup, OR (b) extend the curated DB to 44 targets "
             "(~2 weeks of centroid curation).")
    L.append("")
    L.append("---")
    L.append("")
    L.append("Generated by `scripts/39_v5_pocket_conditional_liability.py` "
             f"(scenario={args.scenario}).")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(L), encoding="utf-8")
    logger.info("Wrote %s", args.report)

    logger.info("Before: %s; After: %s", before, after)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
