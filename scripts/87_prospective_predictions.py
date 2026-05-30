"""Emit the pre-registered prospective-prediction report (review-3 items 2 + 9).

Turns the class prior into falsifiable, time-stamped predictions for REAL ongoing
cognition trials, and scores the ones that have already read out. Output feeds the
OSF pre-registration.

Output: reports/prospective_predictions_v1.md

Usage:
  python scripts/87_prospective_predictions.py
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("prospective")


def main() -> int:
    from mammal_repurposing.reporting import prospective as P

    df = P.load_prospective(ROOT / "data" / "raw" / "prospective_predictions.csv")
    sc = P.score_resolved(df)
    sm = P.summary(df)

    L = []
    L.append("# Pre-registered prospective class predictions")
    L.append("")
    L.append("Falsifiable, time-stamped predictions for **real** ongoing cognition trials, "
             "each following only from the drug's mechanism-class historical track record "
             "(the validated class prior). This is the forward test a retrospective AUROC "
             "cannot be. Frozen 2026-05-30; reconcile NCTs before OSF lock. Reproduced by "
             "`scripts/87_prospective_predictions.py`.")
    L.append("")
    L.append(f"**{sm['n_total']} predictions** across {len(sm['classes'])} mechanism "
             f"classes: {sm['n_resolved']} already resolved, {sm['n_pending']} pending.")
    L.append("")
    L.append("## Resolved since the ledger was curated (out-of-sample confirmations)")
    L.append("")
    if sc["n_resolved"]:
        L.append(f"**{sc['n_correct']} / {sc['n_resolved']} correct** "
                 f"(accuracy {sc['accuracy']:.0%}). These trials read out *after* the "
                 f"31-drug ledger was frozen, so they are genuine out-of-sample tests of "
                 f"class predictions.")
        L.append("")
        L.append("| Drug | class | indication | predicted | actual | ✓ | basis |")
        L.append("|---|---|---|---|---|---|---|")
        for _, r in sc["rows"].iterrows():
            ok = "✓" if r["correct"] else "✗"
            L.append(f"| {r['drug']} | {r['mechanism_class']} | {r['indication']} | "
                     f"{r['predicted_outcome']} | {r['actual_outcome']} | {ok} | "
                     f"{r['prediction_basis']} |")
        L.append("")
        L.append("The headline confirmation is **iclepertin** (GlyT1; CONNEX Phase 3, 2025): "
                 "the NMDA-coagonist-enhancer class had already failed in cognition "
                 "(bitopertin), so the class prior predicted FAILURE — and CONNEX failed its "
                 "MCCB cognition primary, with the programme discontinued. The same axis "
                 "(DAAO inhibitor luvadaxistat) also missed its primary and was halted in "
                 "2024. A class-history prediction made from pre-2020 precedent was thus "
                 "confirmed by 2024–2025 readouts.")
    else:
        L.append("_No resolved predictions yet._")
    L.append("")
    L.append("## Pending — genuinely prospective, falsifiable")
    L.append("")
    L.append("| Drug | trial (NCT) | class | indication | predicted | g range | basis |")
    L.append("|---|---|---|---|---|---|---|")
    for _, r in df[df["status"] == "PENDING"].iterrows():
        L.append(f"| {r['drug']} | {r['trial_program']} ({r['nct']}) | "
                 f"{r['mechanism_class']} | {r['indication']} | "
                 f"**{r['predicted_outcome']}** | {r['pred_g_low']:+.1f}–{r['pred_g_high']:+.1f} | "
                 f"{r['prediction_basis']} |")
    L.append("")
    L.append("Falsification rule (pre-specified): a PENDING prediction is **correct** if the "
             "trial's primary cognitive endpoint outcome (met / not-met on the pre-registered "
             "primary) matches `predicted_outcome`, and **wrong** otherwise. The honest "
             "counter-signal already on record — emraclidine's M4 EMPOWER Phase 2 miss (2024, "
             "a psychosis endpoint) — is retained as a caution that the muscarinic class is "
             "becoming drug-level heterogeneous; it tempers the M1/M4 SUCCESS predictions and "
             "is itself a falsifiable bet that KarXT-class agents fare better on cognition "
             "than emraclidine did on psychosis.")
    L.append("")
    L.append("## Why this matters")
    L.append("")
    L.append("A retrospective AUROC of 1.00 on a curated ledger is, by construction, a "
             "look-up. These predictions are not: they are named drugs, named trials, named "
             "endpoints, and a frozen date. If in 12–18 months the PDE4 (zatolmilast) and "
             "M1/M4 (KarXT-AD) bets resolve as predicted while the NMDA-coagonist axis keeps "
             "failing, the class-prognostic prior will have earned the word *predicts* in a "
             "way no retrospective analysis can. If they do not, that is recorded here against "
             "the method.")
    L.append("")
    L.append("Generated by `scripts/87_prospective_predictions.py`.")
    out = ROOT / "reports" / "prospective_predictions_v1.md"
    out.write_text("\n".join(L), encoding="utf-8")
    logger.info("Resolved: %d/%d correct (%.0f%%); %d pending",
                sc["n_correct"], sc["n_resolved"],
                100 * sc["accuracy"] if sc["n_resolved"] else 0, sm["n_pending"])
    logger.info("Wrote %s", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
