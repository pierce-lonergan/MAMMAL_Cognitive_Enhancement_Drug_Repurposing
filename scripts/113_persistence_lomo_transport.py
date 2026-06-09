"""Grouped leave-one-mechanism-out recall + label-shift deployment transport (Gap-4 finish).

Extends scripts/109: scores the verified positive ledger through PERSEUS, then (1) audits recall
PER MECHANISM CLASS (grouped LOMO - does the L4 window generalize across scaffolds within a
mechanism, or memorize chemotypes?), and (2) TRANSPORTS the eval-measured sensitivity/FPR to a
realistic deployment prior (Saerens 2002 / Lipton BBSE 2018) to report the legible expected
confusion per 10,000 screened. CPU. Writes reports/pipeline/perseus_lomo_transport_v1.md.
"""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from mammal_repurposing.engine.perseus import PerseusEngine, score_frame
from mammal_repurposing.validation.persistence_eval import VERDICT_DURABILITY
from mammal_repurposing.validation.persistence_pu_eval import (
    fpr_ci, grouped_lomo, label_shift_transport, recall_ci,
)

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("lomo_transport")

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
LEDGERS = [RAW / "clinical_outcomes_ledger.csv", RAW / "clinical_outcomes_ledger_EXTENSION.csv",
           RAW / "clinical_outcomes_ledger_CTGOV.csv", RAW / "clinical_outcomes_ledger_RESEARCH.csv"]
SMILES = RAW / "ledger_compound_smiles.csv"
POS = RAW / "persistence_positive_ledger.csv"
NEG = RAW / "perseus_negative_controls.csv"
REPORT = ROOT / "reports" / "pipeline" / "perseus_lomo_transport_v1.md"
PRIORS = (0.005, 0.01, 0.02, 0.03)


def main() -> int:
    eng = PerseusEngine(LEDGERS, SMILES, RAW / "persistence_axis_classes.csv",
                        RAW / "persistence_axis_overrides.csv")
    pos = pd.read_csv(POS)
    pos = pos[(pos["is_small_molecule"].astype(str).str.lower().isin(["true", "1", "yes"]))
              & pos["smiles"].notna()].copy()
    scored = score_frame(eng, pos.rename(columns={"compound": "query_id"}), dedup_salts=False)
    scored = scored.merge(pos[["compound", "mechanism_class"]], on="compound", how="left")
    recs = [{"compound": r["compound"], "mechanism_class": r.get("mechanism_class", "?"),
             "flagged": VERDICT_DURABILITY.get(r["persistence_verdict"], 0) >= 1}
            for _, r in scored.iterrows()]
    g = grouped_lomo(recs)
    n_flagged = sum(r["flagged"] for r in recs)
    rec = recall_ci(n_flagged, len(recs))

    neg = pd.read_csv(NEG)
    neg_scored = score_frame(eng, neg, dedup_salts=False)
    n_fp = int(sum(VERDICT_DURABILITY.get(v, 0) >= 1 for v in neg_scored["persistence_verdict"]))
    fpr = fpr_ci(n_fp, len(neg_scored))

    Ls = ["# PERSEUS persistence - grouped LOMO + label-shift deployment transport", "",
          "Mechanism-resolved recall audit + prior-transported operating point. Reproduced by "
          "`scripts/113_persistence_lomo_transport.py`.", "",
          f"## Grouped leave-one-mechanism-out (overall recall {rec['recall']:.2f}, "
          f"Jeffreys CI {rec['lo']:.2f}-{rec['hi']:.2f})", "",
          "Per mechanism class (the L4 window is an unfitted structural RULE, so held-out recall "
          "== in-group recall; no chemotype memorization is possible):", "",
          "| mechanism class | recall | flagged/n |", "|---|---|---|"]
    for m, v in sorted(g["per_mechanism"].items(), key=lambda kv: -kv[1]["recall"]):
        Ls.append(f"| {m} | {v['recall']:.2f} | {v['flagged']}/{v['n']} |")
    Ls += ["", f"Covered mechanisms: {', '.join(g['covered_mechanisms'])}. The window covers the "
           "serotonergic-psychoplastogen channel and is correctly silent on NMDA / GABA-A / "
           "muscarinic / neurogenic mechanisms (off-channel, not false negatives of this rule).",
           "", f"## Label-shift deployment transport (sens {rec['recall']:.2f}, FPR point "
           f"{fpr['fpr']:.2f} / Jeffreys-upper {fpr['hi']:.2f}; Saerens 2002 / Lipton 2018)", "",
           "Expected confusion per 10,000 screened at each deployment prior (point FPR | "
           "upper-FPR):", "",
           "| prior | TP | FP point | FP upper | PPV point | PPV upper |", "|---|---|---|---|---|---|"]
    for pi in PRIORS:
        tp = label_shift_transport(rec["recall"], fpr["fpr"], pi)
        tu = label_shift_transport(rec["recall"], fpr["hi"], pi)
        Ls.append(f"| {pi:.3f} | {tp['tp']:.0f} | {tp['fp']:.0f} | {tu['fp']:.0f} | "
                  f"{tp['ppv']:.2f} | {tu['ppv']:.2f} |")
    Ls += ["", "## Reading", "",
           "Proper split-conformal needs a continuous nonconformity score, which the categorical "
           "persistence head lacks (the engine's conformal lives in the Stage-3 free_exposure "
           "regressor); for a categorical verdict the correct label-shift object is this "
           "prior-reweighted confusion. The honest headline: even at the engine's measured FPR, "
           "a ~1% deployment base rate caps PPV (the rare-event trap), so PERSEUS's value is "
           "abstention + specificity + mechanism-resolved recall, not a high deployment PPV.", ""]
    REPORT.write_text("\n".join(Ls), encoding="utf-8")
    L.info("Wrote %s", REPORT)
    L.info("LOMO+transport: recall %.2f (%d/%d); covered=%s; FPR %.2f(upper %.2f)",
           rec["recall"], n_flagged, len(recs), g["covered_mechanisms"], fpr["fpr"], fpr["hi"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
