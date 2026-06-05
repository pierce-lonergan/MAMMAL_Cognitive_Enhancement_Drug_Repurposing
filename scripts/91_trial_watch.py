"""Prospective trial-watch: score the standing forward-prediction registry.

Derives the calibrated per-mechanism-class SUCCESS prior from the combined
clinical-outcome ledger (base + extension + CT.gov, n=47), applies it to the
ongoing cognition trials in data/raw/prospective_predictions.csv with each trial
drug held out, and scores the predictions that have read out so far.

Outputs:
  reports/pipeline/trial_watch_v1.md   - the prospective scorecard
  data/raw/trial_watch_registry.csv    - the locked, auto-scored registry

Usage:
  python scripts/91_trial_watch.py
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.reporting.trial_watch import (  # noqa: E402
    load_combined_ledger, class_success_table, build_registry, score_registry,
)

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
logger = logging.getLogger("trial_watch")

ROOT = Path(__file__).resolve().parents[1]
LEDGERS = [
    ROOT / "data" / "raw" / "clinical_outcomes_ledger.csv",
    ROOT / "data" / "raw" / "clinical_outcomes_ledger_EXTENSION.csv",
    ROOT / "data" / "raw" / "clinical_outcomes_ledger_CTGOV.csv",
]
PROSPECTIVE = ROOT / "data" / "raw" / "prospective_predictions.csv"


def _fmt_ci(ci) -> str:
    lo, hi = ci
    if lo != lo:  # nan
        return "n/a"
    return f"[{lo:.2f}, {hi:.2f}]"


def main() -> int:
    ledger = load_combined_ledger(LEDGERS)
    base = float(ledger["label"].mean())
    logger.info("combined ledger: %d drugs, %d classes, base success rate %.3f",
                len(ledger), ledger["mechanism_class"].nunique(), base)

    prior = class_success_table(ledger)
    prospective = pd.read_csv(PROSPECTIVE, comment="#")
    registry = build_registry(prospective, ledger)
    sc = score_registry(registry)

    reg_out = ROOT / "data" / "raw" / "trial_watch_registry.csv"
    registry.to_csv(reg_out, index=False)
    logger.info("wrote registry: %s (%d trials)", reg_out, len(registry))

    # ---- report ----------------------------------------------------------
    L: list[str] = []
    L.append("# Prospective trial-watch (v1)\n")
    L.append("A standing forward-prediction system. The calibrated per-mechanism-"
             "class SUCCESS prior is derived from the combined clinical-outcome "
             f"ledger (n={len(ledger)} drugs, {ledger['mechanism_class'].nunique()} "
             "classes, base success rate "
             f"{base:.2f}); each ongoing cognition trial inherits its class history "
             "with the trial drug held out, so no drug predicts its own outcome. "
             "Unlike the retrospective AUROC, these are falsifiable forward calls on "
             "named trials, checkable as they read out.\n")

    # scorecard
    L.append("## Prospective scorecard (RESOLVED trials)\n")
    if sc["n_resolved"] == 0:
        L.append("No predictions have resolved yet.\n")
    else:
        au = sc["auroc"]
        au_txt = "n/a (needs >=1 success and >=1 failure resolved)" if au != au \
            else f"{au:.2f} {_fmt_ci(sc.get('auroc_ci', (float('nan'), float('nan'))))}"
        L.append(f"- Resolved: **{sc['n_resolved']}** "
                 f"({sc['n_success']} success, {sc['n_failure']} failure)")
        L.append(f"- Accuracy (predicted vs actual): "
                 f"**{sc['accuracy']*100:.0f}%** "
                 f"({sum(v['correct'] for v in sc['by_confidence'].values())}"
                 f"/{sc['n_resolved']})")
        L.append(f"- Prospective AUROC: **{au_txt}**")
        L.append(f"- Brier score: **{sc['brier']:.3f}** "
                 "(lower is better; 0.25 = no-skill at base rate 0.5)")
        if sc["by_confidence"]:
            tiers = ", ".join(f"{c} {v['correct']}/{v['n']}"
                              for c, v in sorted(sc["by_confidence"].items()))
            L.append(f"- By confidence tier: {tiers}")
        L.append("")

    # registry table
    L.append("## Registry: locked forward predictions\n")
    L.append("| Drug | Trial | Class | P(success) | Call | Conf | Evidence | "
             "Status | Actual | Match |")
    L.append("|---|---|---|---|---|---|---|---|---|---|")
    for _, r in registry.iterrows():
        actual = r.get("actual_outcome")
        actual = "" if pd.isna(actual) else str(actual)
        match = ""
        if str(r.get("status")) == "RESOLVED" and actual:
            match = "OK" if str(r["predicted_outcome"]).strip() == actual.strip() else "MISS"
        L.append(
            f"| {r['drug']} | {r.get('trial_program','')} | "
            f"`{r['mechanism_class']}` | {r['p_success']:.2f} | "
            f"{r['predicted_outcome']} | {r['confidence']} | "
            f"{r['n_evidence']} ({r['evidence_level']}) | {r.get('status','')} | "
            f"{actual} | {match} |")
    L.append("")

    # engine-vs-frozen agreement
    agree = int(registry["engine_matches_frozen"].sum())
    L.append("## Engine vs frozen (hand) predictions\n")
    L.append(f"The automated class-prior engine reproduces the originally frozen "
             f"hand predictions on **{agree}/{len(registry)}** trials. "
             "Disagreements are cases where the hand prediction reasoned by "
             "mechanistic analogy beyond direct class evidence; the engine is "
             "deliberately more conservative and tags those LOW confidence or "
             "ABSTAIN.\n")

    # the prior
    L.append("## The class-prior table (calibrated success rate per class)\n")
    L.append("| Mechanism class | n | successes | raw rate | shrunk P(success) |")
    L.append("|---|---|---|---|---|")
    for _, r in prior.iterrows():
        L.append(f"| `{r['mechanism_class']}` | {int(r['n'])} | "
                 f"{int(r['n_success'])} | {r['p_raw']:.2f} | {r['p_shrunk']:.2f} |")
    L.append("")

    L.append("## Honest scope\n")
    L.append("- The prior makes a clean, falsifiable call only where a class has "
             ">=2 prior members (HIGH). Singletons, same-drug continuations, and "
             "shared-axis (super-class) borrows are flagged MED/LOW so a reader can "
             "see exactly how much evidence each call rests on.")
    L.append("- Mechanistic super-classes are pre-specified by shared pharmacology "
             "(for example GlyT1 and DAAO inhibitors both enhance NMDA "
             "co-agonism), never by outcome, and never lump axes with opposite "
             "track records.")
    L.append("- Prospective AUROC is reported only once both a success and a "
             "failure have resolved; until then accuracy and Brier carry the "
             "score. The registry is the accruing record: re-run as trials read "
             "out.")
    L.append("")

    out = ROOT / "reports" / "pipeline" / "trial_watch_v1.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(L), encoding="utf-8")
    logger.info("wrote report: %s", out)

    # console summary
    print("\n=== trial-watch ===")
    print(f"ledger n={len(ledger)}  classes={ledger['mechanism_class'].nunique()}  "
          f"base={base:.3f}")
    print(f"registry trials={len(registry)}  resolved={sc['n_resolved']}  "
          f"engine==frozen={agree}/{len(registry)}")
    if sc["n_resolved"]:
        print(f"accuracy={sc['accuracy']*100:.0f}%  brier={sc['brier']:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
