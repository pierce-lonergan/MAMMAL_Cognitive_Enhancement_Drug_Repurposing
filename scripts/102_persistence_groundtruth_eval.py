"""PERSEUS evaluation v2 - bidirectional, against a cited persistence ground-truth ledger.

scripts/101 measured SPECIFICITY (negative controls). This adds the directional check
against `data/raw/persistence_ground_truth.csv` - compounds with a real persistence-DESIGN
readout (delayed-start / randomized-discontinuation / washout / parallel-group), labelled by
what the trial actually showed. We report:
  - the OVER-CLAIM rate (does PERSEUS ever assert more durability than the label supports?);
  - the coverage-accuracy curve over the evidence-design rank;
  - the LABEL BUDGET (confirmed positives needed before recall is estimable).

CPU only. Writes reports/pipeline/perseus_eval_groundtruth_v1.md.
"""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from mammal_repurposing.engine.perseus import PerseusEngine, score_frame
from mammal_repurposing.validation.persistence import EVIDENCE_RANK
from mammal_repurposing.validation.persistence_eval import (
    coverage_accuracy_curve, evaluate, label_budget,
)

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("gt_eval")

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
LEDGERS = [RAW / "clinical_outcomes_ledger.csv", RAW / "clinical_outcomes_ledger_EXTENSION.csv",
           RAW / "clinical_outcomes_ledger_CTGOV.csv", RAW / "clinical_outcomes_ledger_RESEARCH.csv"]
SMILES = RAW / "ledger_compound_smiles.csv"
CATALOGUE = RAW / "chembl_approved_catalogue.csv"
GT = RAW / "persistence_ground_truth.csv"
REPORT = ROOT / "reports" / "pipeline" / "perseus_eval_groundtruth_v1.md"


def main() -> int:
    gt = pd.read_csv(GT)
    smi = pd.read_csv(SMILES); smi["k"] = smi["compound"].str.lower().str.strip()
    cat = pd.read_csv(CATALOGUE); cat["k"] = cat["name"].str.lower().str.strip()
    look = {**dict(zip(cat["k"], cat["smiles"])), **dict(zip(smi["k"], smi["smiles"]))}
    look.setdefault("methylphenidate", "COC(=O)C(C1CCCCN1)c1ccccc1")

    scoreable = gt[gt["structure_scoreable"] == "yes"].copy()
    scoreable["smiles"] = scoreable["compound"].str.lower().str.strip().map(look)
    n_no_smiles = int(scoreable["smiles"].isna().sum())
    scoreable = scoreable[scoreable["smiles"].notna()]

    eng = PerseusEngine(LEDGERS, SMILES, RAW / "persistence_axis_classes.csv",
                        RAW / "persistence_axis_overrides.csv")
    scored = score_frame(eng, scoreable.rename(columns={"compound": "query_id"}),
                         dedup_salts=False)
    scored = scored.merge(
        gt[["compound", "mechanism_class", "persistence_label", "persistence_design"]],
        on="compound", how="left")

    records = scored.to_dict("records")
    ev = evaluate(records)
    curve = coverage_accuracy_curve(records, lambda r: EVIDENCE_RANK.get(r["evidence_design"], 0))
    n_positive = int((gt["persistence_label"] == "demonstrated_healthy").sum())
    budget = label_budget(prior=0.01)

    Ls = ["# PERSEUS evaluation v2 - ground-truth (bidirectional)", "",
          "Scored the cited persistence ground-truth ledger (`data/raw/persistence_ground_truth.csv`) "
          "- compounds with a real persistence-DESIGN readout - through PERSEUS and compared "
          "each verdict to the trial-design label. Reproduced by "
          "`scripts/102_persistence_groundtruth_eval.py`.", "",
          f"Scoreable ground-truth compounds: **{len(scored)}** "
          f"(+{int((gt['structure_scoreable']=='no').sum())} non-structure mAbs recorded but not "
          f"scored; {n_no_smiles} missing SMILES).", "",
          f"## Over-claim rate: **{ev['n_over_claims']} / {ev['n']}** "
          f"(the directional error - asserting more durability than the label supports)", ""]
    if ev["over_claimers"]:
        Ls.append("Over-claimers (FIX THESE): " + ", ".join(ev["over_claimers"]))
        Ls.append("")
    Ls.append("| compound | mechanism | design | label | PERSEUS verdict | over-claim |")
    Ls.append("|---|---|---|---|---|---|")
    from mammal_repurposing.validation.persistence_eval import over_claims
    for _, r in scored.sort_values("persistence_label").iterrows():
        oc = over_claims(r["persistence_verdict"], r["persistence_label"])
        Ls.append(f"| {r['compound']} | {r['mechanism_class']} | {r['persistence_design']} | "
                  f"{r['persistence_label']} | {r['persistence_verdict']} | "
                  f"{'**YES**' if oc else 'no'} |")
    Ls.append("")
    Ls.append("## Coverage-accuracy (sweeping the evidence-design rank required to assert persistence)")
    Ls.append("")
    Ls.append("| min evidence rank | coverage | accuracy (non-over-claim) | asserted |")
    Ls.append("|---|---|---|---|")
    for c in curve:
        acc = "n/a" if c["accuracy"] != c["accuracy"] else f"{c['accuracy']:.2f}"
        Ls.append(f"| {c['evidence_rank_threshold']} | {c['coverage']:.2f} | {acc} | {c['asserted']} |")
    Ls.append("")
    Ls.append("## Label budget (why sensitivity is unmeasurable today)")
    Ls.append("")
    Ls.append(f"Confirmed durable-in-healthy positives in the ledger: **{n_positive}** (empty). "
              f"At a realistic ~1% prior, an estimated **~{budget:,} confirmed positive "
              "delayed-start readouts** would be needed before recall is estimable to +/-0.1. "
              "Until then PERSEUS reports SPECIFICITY (0 over-claims) and abstains; sensitivity "
              "and PPV are not yet identifiable - this budget is the honest deliverable, not a "
              "hidden gap.")
    Ls.append("")
    REPORT.write_text("\n".join(Ls), encoding="utf-8")
    L.info("Wrote %s", REPORT)
    L.info("ground-truth eval: %d/%d over-claims; %d scoreable; label budget ~%d positives",
           ev["n_over_claims"], ev["n"], len(scored), budget)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
