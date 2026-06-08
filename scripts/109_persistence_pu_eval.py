"""Rigorous bidirectional empty-positive evaluation of the PERSEUS persistence head (Gap 4).

Scores BOTH the verified positive ledger (recall) and the negative-control panel (FPR) through
PERSEUS with a single "flagged" rule (verdict asserts any durability, level >= 1), then reports:
  - sensitivity with a JEFFREYS 95% CI (small-sample-correct, not Wald);
  - FPR with a Jeffreys CI;
  - PPV as a CURVE over prior pi in [0.005, 0.03], at both the point and Jeffreys-upper FPR.
This is the publishable bidirectional metric no current repurposing persistence predictor reports.

CPU. Writes reports/pipeline/perseus_pu_eval_v1.md.
"""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from mammal_repurposing.engine.perseus import PerseusEngine, score_frame
from mammal_repurposing.validation.persistence_eval import VERDICT_DURABILITY
from mammal_repurposing.validation.persistence_pu_eval import evaluate

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("perseus_pu_eval")

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
LEDGERS = [RAW / "clinical_outcomes_ledger.csv", RAW / "clinical_outcomes_ledger_EXTENSION.csv",
           RAW / "clinical_outcomes_ledger_CTGOV.csv", RAW / "clinical_outcomes_ledger_RESEARCH.csv"]
SMILES = RAW / "ledger_compound_smiles.csv"
POS = RAW / "persistence_positive_ledger.csv"
NEG = RAW / "perseus_negative_controls.csv"
REPORT = ROOT / "reports" / "pipeline" / "perseus_pu_eval_v1.md"


def _flagged(scored: pd.DataFrame) -> int:
    return int(sum(VERDICT_DURABILITY.get(v, 0) >= 1 for v in scored["persistence_verdict"]))


def main() -> int:
    if not (POS.exists() and NEG.exists()):
        L.warning("Need both %s and %s", POS, NEG)
        return 0
    eng = PerseusEngine(LEDGERS, SMILES, RAW / "persistence_axis_classes.csv",
                        RAW / "persistence_axis_overrides.csv")

    pos = pd.read_csv(POS)
    pos = pos[(pos["is_small_molecule"].astype(str).str.lower().isin(["true", "1", "yes"]))
              & pos["smiles"].notna()].copy()
    pos_scored = score_frame(eng, pos.rename(columns={"compound": "query_id"}), dedup_salts=False)
    n_flagged, n_pos = _flagged(pos_scored), len(pos_scored)

    neg = pd.read_csv(NEG)
    neg_scored = score_frame(eng, neg, dedup_salts=False)
    n_fp, n_neg = _flagged(neg_scored), len(neg_scored)

    ev = evaluate(n_flagged, n_pos, n_fp, n_neg)
    rec, fpr = ev["recall"], ev["fpr"]

    Ls = ["# PERSEUS persistence - rigorous empty-positive evaluation (Gap 4)", "",
          "Bidirectional small-sample metrics on the verified positive ledger + negative-control "
          "panel, using Jeffreys intervals (not Wald) and a PPV-vs-prior curve. Reproduced by "
          "`scripts/109_persistence_pu_eval.py`.", "",
          f"## Sensitivity (recall): **{rec['recall']:.2f}** "
          f"(Jeffreys 95% CI {rec['lo']:.2f}-{rec['hi']:.2f}), {rec['n_flagged']}/{rec['n']} "
          "verified positives flagged (durability verdict >= 1).", "",
          f"## FPR (negative controls): **{fpr['fpr']:.2f}** "
          f"(Jeffreys 95% CI {fpr['lo']:.2f}-{fpr['hi']:.2f}), {fpr['n_fp']}/{fpr['n_neg']}.", "",
          "## PPV across an externally supplied prior (PPV = pi*S / (pi*S + (1-pi)*FPR))", "",
          "| prior pi | PPV @ point FPR | PPV @ Jeffreys-upper FPR |", "|---|---|---|"]
    up = {d["prior"]: d["ppv"] for d in ev["ppv_at_upper_fpr"]}
    for d in ev["ppv_at_point_fpr"]:
        pv_p = "n/a" if d["ppv"] != d["ppv"] else f"{d['ppv']:.2f}"
        pv_u = up.get(d["prior"])
        pv_u = "n/a" if pv_u is None or pv_u != pv_u else f"{pv_u:.2f}"
        Ls.append(f"| {d['prior']:.3f} | {pv_p} | {pv_u} |")
    Ls += ["", "## Caveats (load-bearing)", "",
           "- " + ev["scar_sar_caveat"],
           "- At n=" + str(rec["n"]) + " the recall CI is wide by construction (~+/-0.25); this "
           "is the honest precision limit, reported per Brown-Cai-DasGupta (Jeffreys), and is why "
           "the old Wald label_budget framing was replaced.",
           "- 'Flagged' counts WINDOW_CONDITIONAL (a permissive plasticity window) as durability "
           ">= 1; recall is therefore recall on the SEROTONERGIC-psychoplastogen sub-class the L4 "
           "window covers, not on every durable mechanism (NMDA/TrkB-TMD/GABA-A are off-channel).", ""]
    REPORT.write_text("\n".join(Ls), encoding="utf-8")
    L.info("Wrote %s", REPORT)
    L.info("PU-eval: recall %.2f (%.2f-%.2f) %d/%d | FPR %.2f (%.2f-%.2f) %d/%d",
           rec["recall"], rec["lo"], rec["hi"], n_flagged, n_pos,
           fpr["fpr"], fpr["lo"], fpr["hi"], n_fp, n_neg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
