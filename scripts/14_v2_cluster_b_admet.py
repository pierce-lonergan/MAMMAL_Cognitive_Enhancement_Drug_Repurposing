"""V2 Cluster B — Score every compound with ADMET-AI and apply hard gates.

Replaces the v1 BBBP/ClinTox heads (which we found to be poorly calibrated).
ADMET-AI runs on CPU; safe to schedule in parallel with GPU stages.

Output:
    data/results/v2/admet_predictions.parquet  # 41 endpoints per SMILES (raw)
    data/results/v2/admet_gates.parquet        # gate status + composite ADMET_score
    data/results/v2/admet_gate_summary.md      # human-readable report
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

from mammal_repurposing.cluster_b.admet_ai_runner import predict_admet_for_compounds  # noqa: E402
from mammal_repurposing.config import COMPOUNDS_PARQUET, RESULTS_DIR, ensure_dirs  # noqa: E402
from mammal_repurposing.gates.admet_gates import (  # noqa: E402
    apply_gates,
    validate_positive_controls,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("v2_cluster_b")

V2_RESULTS_DIR = RESULTS_DIR / "v2"
DEFAULT_PRED_OUT = V2_RESULTS_DIR / "admet_predictions.parquet"
DEFAULT_GATE_OUT = V2_RESULTS_DIR / "admet_gates.parquet"
DEFAULT_REPORT_OUT = V2_RESULTS_DIR / "admet_gate_summary.md"


def _render_summary(gates: pd.DataFrame, controls_ok: bool, failures: list[str]) -> str:
    n_total = len(gates)
    n_pass = (gates["gate_status"] == "PASS").sum()
    n_flag = (gates["gate_status"] == "FLAG").sum()
    n_cut = (gates["gate_status"] == "CUT").sum()

    lines = []
    lines.append("# V2 Cluster B — ADMET Gate Summary")
    lines.append("")
    lines.append(f"- Total compounds scored: **{n_total}**")
    lines.append(f"- PASS: **{n_pass}** ({100*n_pass/n_total:.1f}%)")
    lines.append(f"- FLAG: **{n_flag}** ({100*n_flag/n_total:.1f}%)")
    lines.append(f"- CUT:  **{n_cut}** ({100*n_cut/n_total:.1f}%)")
    lines.append("")
    lines.append(f"- Positive-control validation: {'✅ all pass BBB' if controls_ok else '❌ FAILURES'}")
    if not controls_ok:
        lines.append("")
        lines.append("**Failures**:")
        for f in failures:
            lines.append(f"  - {f}")
    lines.append("")

    # Top cuts by frequency
    cut_compounds = gates[gates["gate_status"] == "CUT"]
    if not cut_compounds.empty:
        lines.append("## Top reasons for CUT")
        lines.append("")
        all_failures = cut_compounds["gates_failed"].str.split(";").explode()
        # Extract gate name (before "=")
        all_failures = all_failures[all_failures != ""].str.split("=").str[0]
        reason_counts = all_failures.value_counts()
        lines.append("| Gate | # cuts |")
        lines.append("|---|---|")
        for gate, n in reason_counts.items():
            lines.append(f"| {gate} | {n} |")
        lines.append("")

    lines.append("## ADMET-score distribution (passes only)")
    lines.append("")
    passing = gates[gates["gate_status"] != "CUT"]
    if not passing.empty:
        stats = passing["admet_score"].describe()
        lines.append("| stat | value |")
        lines.append("|---|---|")
        for k in ["min", "25%", "50%", "75%", "max", "mean", "std"]:
            v = stats.get(k)
            if v is not None:
                lines.append(f"| {k} | {v:.3f} |")
        lines.append("")

    lines.append("## Top 15 by ADMET score (passes only)")
    lines.append("")
    if not passing.empty:
        top = passing.sort_values("admet_score", ascending=False).head(15)
        lines.append("| compound | admet_score | gate_status | flags |")
        lines.append("|---|---|---|---|")
        for _, r in top.iterrows():
            flags = r.get("gates_flagged") or ""
            lines.append(f"| {r['compound_name']} | {r['admet_score']:.3f} | {r['gate_status']} | {flags} |")

    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compounds", type=Path, default=COMPOUNDS_PARQUET)
    parser.add_argument("--pred-out", type=Path, default=DEFAULT_PRED_OUT)
    parser.add_argument("--gates-out", type=Path, default=DEFAULT_GATE_OUT)
    parser.add_argument("--report-out", type=Path, default=DEFAULT_REPORT_OUT)
    parser.add_argument("--no-cache", action="store_true",
                        help="Disable per-SMILES ADMET cache.")
    parser.add_argument("--limit", type=int, default=None,
                        help="Only score first N compounds (smoke test).")
    args = parser.parse_args()

    ensure_dirs()
    V2_RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    if not args.compounds.exists():
        logger.error("Compounds parquet not found: %s", args.compounds)
        return 1

    compounds = pd.read_parquet(args.compounds)
    if args.limit:
        compounds = compounds.head(args.limit)
    logger.info("Scoring ADMET for %d compounds ...", len(compounds))

    preds = predict_admet_for_compounds(compounds, use_cache=not args.no_cache)
    preds.to_parquet(args.pred_out, index=False)
    logger.info("Wrote %d ADMET predictions (%d columns) -> %s",
                len(preds), len(preds.columns), args.pred_out)

    # Pass compounds_df so apply_gates can grant regulatory bypass to approved drugs.
    gates = apply_gates(preds, compounds_df=compounds)
    gates.to_parquet(args.gates_out, index=False)
    logger.info("Applied gates: %d PASS / %d FLAG / %d CUT",
                (gates["gate_status"] == "PASS").sum(),
                (gates["gate_status"] == "FLAG").sum(),
                (gates["gate_status"] == "CUT").sum())

    controls_ok, failures = validate_positive_controls(gates)
    if not controls_ok:
        logger.warning("Positive-control validation FAILED: %s", failures)
    else:
        logger.info("Positive-control BBB validation: ALL PASS")

    args.report_out.write_text(_render_summary(gates, controls_ok, failures),
                                encoding="utf-8")
    logger.info("Wrote gate summary -> %s", args.report_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
