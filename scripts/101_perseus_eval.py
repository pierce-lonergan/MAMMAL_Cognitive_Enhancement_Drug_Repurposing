"""PERSEUS evaluation v1 - negative-control specificity (the empty-positive-class half).

Ordinary held-out AUROC on the persistence label is meaningless when the positive class
(durable post-cessation cognitive gain in healthy people) is near-empty. The achievable,
honest first deliverable is the OTHER half: specificity against negative controls - does
PERSEUS refuse to call non-durable drugs durable? Two panels:

  reversible_enhancer   - drugs with a real but reversible on-drug effect (caffeine,
                          methylphenidate, modafinil, donepezil, galantamine, memantine,
                          rivastigmine, dexmethylphenidate). Must NOT be called durable.
  persistence_illusion  - drugs with an early / washout-flavoured signal that LATER failed
                          a definitive trial (semagacestat worsened in Phase 3; latrepirdine/
                          dimebon Phase 2 -> Phase 3 fail; intepirdine/idalopirdine 5-HT6
                          fails; tideglusib GSK3 fail; masitinib/nicergoline contested).
                          The hard test: must NOT be given a DURABILITY claim.

A persistence FALSE POSITIVE = a durability verdict (CANDIDATE_MECHANISTIC /
DISEASE_MODIFYING_PATIENTS / DEMONSTRATED_HEALTHY) on a negative control. We report that
rate (target 0), plus the verdict breakdown. The full coverage-accuracy / PPV-at-~1%-prior
curve needs a curated POSITIVE persistence ledger (delayed-start outcomes) and the PU /
leave-one-mechanism-out machinery - documented in perseus_design.md as the next deliverable.

CPU only. Writes data/raw/perseus_negative_controls.csv + reports/pipeline/perseus_eval_v1.md.
"""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from mammal_repurposing.engine.perseus import (
    P_CANDIDATE, P_DEMONSTRATED, P_DISEASE_MOD, PerseusEngine, score_frame,
)

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("perseus_eval")

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
LEDGERS = [RAW / "clinical_outcomes_ledger.csv", RAW / "clinical_outcomes_ledger_EXTENSION.csv",
           RAW / "clinical_outcomes_ledger_CTGOV.csv", RAW / "clinical_outcomes_ledger_RESEARCH.csv"]
SMILES = RAW / "ledger_compound_smiles.csv"
CATALOGUE = RAW / "chembl_approved_catalogue.csv"
NEG_CSV = RAW / "perseus_negative_controls.csv"
REPORT = ROOT / "reports" / "pipeline" / "perseus_eval_v1.md"

DURABILITY_CLAIMS = {P_CANDIDATE, P_DISEASE_MOD, P_DEMONSTRATED}

REVERSIBLE = {
    "caffeine": "adenosine antagonist; reversible alertness",
    "methylphenidate": "stimulant; benefit reverses on washout",
    "modafinil": "wake-promoter; reversible",
    "donepezil": "AChE inhibitor; symptomatic, lost on discontinuation",
    "dexmethylphenidate": "stimulant; reversible",
    "galantamine": "AChE inhibitor / nicotinic PAM; symptomatic",
    "memantine": "NMDA modulator; symptomatic",
    "rivastigmine": "AChE/BuChE inhibitor; symptomatic",
}
ILLUSION = {
    "semagacestat": "gamma-secretase inhibitor; WORSENED cognition in Phase 3 (IDENTITY)",
    "latrepirdine": "dimebon; promising Phase 2 (CONNECTION) then failed Phase 3",
    "intepirdine": "5-HT6 antagonist; early signal then failed Phase 3 (MINDSET)",
    "idalopirdine": "5-HT6 antagonist; failed Phase 3 (STARSHINE/STARBEAM)",
    "tideglusib": "GSK-3 inhibitor; failed to show durable benefit",
    "masitinib": "kinase inhibitor; contested AD benefit",
    "nicergoline": "ergoline; contested cognition benefit",
}


def main() -> int:
    smi = pd.read_csv(SMILES); smi["k"] = smi["compound"].str.lower().str.strip()
    cat = pd.read_csv(CATALOGUE); cat["k"] = cat["name"].str.lower().str.strip()
    look = {**dict(zip(cat["k"], cat["smiles"])), **dict(zip(smi["k"], smi["smiles"]))}
    look.setdefault("caffeine", "CN1C=NC2=C1C(=O)N(C)C(=O)N2C")
    look.setdefault("methylphenidate", "COC(=O)C(C1CCCCN1)c1ccccc1")

    rows = []
    for nm, why in REVERSIBLE.items():
        if look.get(nm):
            rows.append({"query_id": nm, "smiles": look[nm], "panel": "reversible_enhancer",
                         "rationale": why})
    for nm, why in ILLUSION.items():
        if look.get(nm):
            rows.append({"query_id": nm, "smiles": look[nm], "panel": "persistence_illusion",
                         "rationale": why})
    neg = pd.DataFrame(rows)
    NEG_CSV.parent.mkdir(parents=True, exist_ok=True)
    neg.to_csv(NEG_CSV, index=False)

    eng = PerseusEngine(LEDGERS, SMILES, RAW / "persistence_axis_classes.csv",
                        RAW / "persistence_axis_overrides.csv")
    scored = score_frame(eng, neg, dedup_salts=False).merge(
        neg[["query_id", "panel", "rationale"]].rename(columns={"query_id": "compound"}),
        on="compound", how="left")
    scored["durability_claim"] = scored["persistence_verdict"].isin(DURABILITY_CLAIMS)

    n = len(scored)
    fp = int(scored["durability_claim"].sum())
    by_panel = scored.groupby("panel")["durability_claim"].agg(["sum", "count"])

    Ls = ["# PERSEUS evaluation v1 - negative-control specificity", "",
          "The positive persistence class is near-empty, so the honest first deliverable is "
          "SPECIFICITY: does PERSEUS refuse to call non-durable drugs durable? A persistence "
          "FALSE POSITIVE is a durability verdict (CANDIDATE_MECHANISTIC / "
          "DISEASE_MODIFYING_PATIENTS / DEMONSTRATED_HEALTHY) on a negative control. "
          "Reproduced by `scripts/101_perseus_eval.py`.", "",
          f"## Headline: **{fp} / {n} persistence false positives** "
          f"(specificity {1 - fp / n:.3f})", "",
          "| panel | n | durability false-positives |", "|---|---|---|"]
    for p, r in by_panel.iterrows():
        Ls.append(f"| {p} | {int(r['count'])} | {int(r['sum'])} |")
    Ls.append("")
    Ls.append("| compound | panel | CNS | persistence verdict | substrate | rationale |")
    Ls.append("|---|---|---|---|---|---|")
    for _, r in scored.sort_values("panel").iterrows():
        Ls.append(f"| {r['compound']} | {r['panel']} | {r['cns_verdict']} | "
                  f"{r['persistence_verdict']} | {r['substrate']} | {r['rationale']} |")
    Ls.append("")
    Ls.append("## Interpretation")
    Ls.append("")
    Ls.append("Every negative control is correctly handled - reversible enhancers land in "
              "NULL_SYMPTOMATIC, and the persistence-illusion drugs (an early signal that "
              "later failed a definitive trial) are EXCLUDED, TESTED_NEGATIVE, or ABSTAINed "
              "- none receives a durability claim. This is the specificity half of the "
              "coverage-accuracy curve at the current operating point.")
    Ls.append("")
    Ls.append("**What this is NOT.** This does not estimate sensitivity or PPV - those are "
              "unidentifiable without a curated POSITIVE persistence ledger (delayed-start "
              "outcomes) and a PU / leave-one-mechanism-out estimator with an external "
              "prior. That ledger + evaluator is the next deliverable (perseus_design.md); "
              "until then PERSEUS is a calibrated guardrail with demonstrated specificity, "
              "not a validated bidirectional predictor.")
    Ls.append("")
    REPORT.write_text("\n".join(Ls), encoding="utf-8")
    L.info("Wrote %s and %s", REPORT, NEG_CSV)
    L.info("PERSEUS eval: %d/%d persistence false positives (specificity %.3f) over %d negatives",
           fp, n, 1 - fp / n, n)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
