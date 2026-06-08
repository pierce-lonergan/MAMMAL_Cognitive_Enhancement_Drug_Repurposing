"""Stage 4 of the persistence-target DTI module: MW-RESIDUALIZED re-calibration.

scripts/104 showed MAMMAL's pKd is heavily molecular-weight-driven (corr(MW,pKd)~0.6-0.9),
and that with size-matched negatives only the two ablative BH3-mimetic channels survive. This
asks the deeper question: is that a "MAMMAL has NO size-independent persistence signal"
result, or was a real signal merely MASKED by size? We de-confound by residualizing each
score against a size->score line fit on the non-engagers, then re-calibrate on the residuals.

A genuine, size-independent binder scores ABOVE what its weight predicts (positive residual).
If a channel that failed raw is RESCUED here, MAMMAL had hidden signal; if it still fails,
the negative is real (no size-independent signal). CPU only - reuses the scores from 104.

CAVEAT (reported per target): residualizing engagers whose MW lies OUTSIDE the non-engager
range is extrapolation - trustworthy mainly where engager/non-engager weights overlap.

Writes data/results/persistence_dti_calibration_mwresid.json +
reports/pipeline/persistence_dti_mwresidual.md.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from mammal_repurposing.engine.persistence_dti import (
    calibrate_target, load_panel, mw_baseline, mw_residualize,
)

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("persistence_dti_mwresid")

ROOT = Path(__file__).resolve().parents[1]
INTERIM = ROOT / "data" / "interim"
RESULTS = ROOT / "data" / "results"
PANEL = INTERIM / "persistence_targets.csv"
SCORES = RESULTS / "persistence_dti_scores.csv"
RAW_CALIB = RESULTS / "persistence_dti_calibration.json"
OUT = RESULTS / "persistence_dti_calibration_mwresid.json"
REPORT = ROOT / "reports" / "pipeline" / "persistence_dti_mwresidual.md"


def main() -> int:
    panel = load_panel(PANEL)
    scored = pd.read_csv(SCORES)
    raw = json.load(open(RAW_CALIB, encoding="utf-8"))["per_target"]

    rows, rescued, survived = [], [], []
    for gene, tgt in panel.items():
        sub = scored[scored["target_gene"] == gene]
        pos = sub[sub["role"] == "engager"].dropna(subset=["predicted_pkd", "mw"])
        neg = sub[sub["role"] == "non_engager"].dropna(subset=["predicted_pkd", "mw"])
        base = mw_baseline(neg["predicted_pkd"], neg["mw"])
        raw_auroc = raw.get(gene, {}).get("auroc")
        raw_pass = raw.get(gene, {}).get("passed", False)
        if base is None or pos.empty:
            rows.append({"gene": gene, "tier": tgt.tier, "raw_auroc": raw_auroc,
                         "raw_pass": raw_pass, "resid_auroc": None, "resid_pass": False,
                         "engager_in_range_frac": None, "slope": None})
            continue
        rp = mw_residualize(pos["predicted_pkd"], pos["mw"], base)
        rn = mw_residualize(neg["predicted_pkd"], neg["mw"], base)
        cal = calibrate_target(rp.tolist(), rn.tolist())
        lo, hi = float(neg["mw"].min()), float(neg["mw"].max())
        in_range = float(((pos["mw"] >= lo) & (pos["mw"] <= hi)).mean())
        cal.update({"tier": tgt.tier, "promotes_durable": tgt.promotes_durable,
                    "raw_auroc": raw_auroc, "raw_passed": raw_pass,
                    "engager_in_range_frac": round(in_range, 2),
                    "mw_slope": round(base[0], 4)})
        rows.append({"gene": gene, "tier": tgt.tier, "raw_auroc": raw_auroc,
                     "raw_pass": raw_pass, "resid_auroc": cal["auroc"],
                     "resid_pass": cal["passed"], "engager_in_range_frac": round(in_range, 2),
                     "slope": round(base[0], 4)})
        if cal["passed"] and not raw_pass:
            rescued.append(gene)
        if cal["passed"] and raw_pass:
            survived.append(gene)

    per_target = {r["gene"]: r for r in rows}
    OUT.write_text(json.dumps({"per_target": per_target,
                               "rescued_by_residualization": rescued,
                               "survived_residualization": survived}, indent=2),
                   encoding="utf-8")
    write_report(rows, rescued, survived)
    L.info("MW-residualized: rescued=%s survived=%s", rescued or "none", survived or "none")
    return 0


def write_report(rows, rescued, survived) -> None:
    def f(x):
        return "n/a" if x is None or (isinstance(x, float) and x != x) else f"{x:.2f}"
    Ls = ["# Persistence-target DTI module - MW-residualized re-calibration", "",
          "Does MAMMAL have ANY molecular-size-INDEPENDENT persistence-substrate signal? Each "
          "score is residualized against a size->score line fit on the non-engagers, then the "
          "channel is re-calibrated on residuals (a genuine binder scores above its size-"
          "expected pKd). Reproduced by `scripts/106_persistence_dti_mwresidual.py`.", "",
          f"**Rescued by de-confounding (failed raw, pass residualized): "
          f"{', '.join(rescued) or 'NONE'}.** "
          f"Survived (passed both): {', '.join(survived) or 'NONE'}.", "",
          "| target | tier | raw AUROC | raw PASS | residualized AUROC | resid PASS | "
          "engager-in-MW-range | size slope |",
          "|---|---|---|---|---|---|---|---|"]
    for r in sorted(rows, key=lambda x: -(x["resid_auroc"] or -1)):
        slope_s = "n/a" if r["slope"] is None else f"{r['slope']:.4f}"
        Ls.append(f"| {r['gene']} | {r['tier']} | {f(r['raw_auroc'])} | "
                  f"{'PASS' if r['raw_pass'] else 'fail'} | {f(r['resid_auroc'])} | "
                  f"{'**PASS**' if r['resid_pass'] else 'fail'} | "
                  f"{f(r['engager_in_range_frac'])} | {slope_s} |")
    Ls += ["", "## Reading", ""]
    if not rescued and not survived:
        Ls.append("After removing the molecular-weight confound, NO channel separates engagers "
                  "from non-engagers - the size-matched negative result was REAL: MAMMAL has no "
                  "size-independent persistence-substrate signal on this panel. The persistence "
                  "head correctly remains abstain-by-default.")
    else:
        Ls.append("Channels passing AFTER residualization carry size-INDEPENDENT signal and are "
                  "the trustworthy substrate channels; any 'rescued' channel was real but masked "
                  "by size in the raw calibration. The `engager-in-MW-range` column flags "
                  "validity: a low value means the residualized AUROC relies on extrapolating "
                  "the size line beyond the non-engager weight range (treat with caution; the "
                  "BH3-mimetics are far larger than the negatives).")
    Ls.append("")
    REPORT.write_text("\n".join(Ls), encoding="utf-8")
    L.info("Wrote %s", REPORT)


if __name__ == "__main__":
    raise SystemExit(main())
