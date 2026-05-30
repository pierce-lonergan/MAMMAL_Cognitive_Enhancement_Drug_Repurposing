"""Gap 5 — Generate clinician-legible evidence dossiers.

Produces a one-page GRADE-style evidence card per (compound, indication) for a
curated, story-telling set: SUCCESS exemplars (donepezil, memantine,
methylphenidate, pitolisant), FAILURE-class drugs the dossier should correctly
DOWN-grade and flag (idalopirdine, encenicline), and a repurposing pick that
exercises the allosteric-reliability caveat.

Output: reports/clinician_dossiers_v1.md

Usage:
  python scripts/80_clinician_dossier.py
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

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("clinician_dossier")

# (compound, indication, target_uniprot) — curated to span the evidence spectrum
DOSSIER_SET = [
    ("donepezil", "AD", "P22303"),
    ("memantine", "AD", "Q13224"),
    ("methylphenidate", "ADHD", "Q01959"),
    ("pitolisant", "narcolepsy", "Q9Y5N1"),
    ("idalopirdine", "AD", "P50406"),     # 5-HT6 FAILURE class — should down-grade + warn
    ("encenicline", "CIAS", "P36544"),    # α7 FAILURE class — should down-grade + warn
    ("galantamine", "AD", "P22303"),      # AChE-I + α7 PAM
]

# Targets where MAMMAL's sequence-only binding is structurally unreliable
# (allosteric / transporter — the Gap-4 finding).
ALLOSTERIC_BLIND = {"P36544", "Q08499", "O76083", "P42261", "P42262", "P42263",
                    "P48058", "Q01959", "P23975", "P50406", "P08908", "Q13639"}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ledger", type=Path,
                    default=ROOT / "data" / "raw" / "clinical_outcomes_ledger.csv")
    ap.add_argument("--anchors", type=Path,
                    default=ROOT / "data" / "raw" / "modulator_anchors_seed.csv")
    ap.add_argument("--liability", type=Path,
                    default=ROOT / "data" / "results" / "liability_dti.parquet")
    ap.add_argument("--panel", type=Path,
                    default=ROOT / "data" / "interim" / "targets.parquet")
    ap.add_argument("--report", type=Path,
                    default=ROOT / "reports" / "clinician_dossiers_v1.md")
    args = ap.parse_args()

    from mammal_repurposing.validation import retrospective as R
    from mammal_repurposing.validation import disease_reframe as D
    from mammal_repurposing.reporting import clinician_dossier as Dx

    ledger = R.load_clinical_ledger(args.ledger)
    anchors = pd.read_csv(args.anchors, comment="#")
    liability = pd.read_parquet(args.liability) if args.liability.exists() else None
    panel = pd.read_parquet(args.panel)
    panel_u = set(panel["uniprot"].astype(str))
    gene = dict(zip(panel["uniprot"].astype(str), panel["gene"]))
    # extend gene map with off-panel targets (HTR6 etc.) from the ledger
    for _, r in ledger.iterrows():
        gene.setdefault(str(r["target_uniprot"]), str(r["target_uniprot"]))

    evidence = D.load_disease_evidence(ledger, anchors)
    priors_by_disease = {d: D.build_disease_class_priors(d, evidence)
                         for d in ("AD", "CIAS", "ADHD", "narcolepsy", "FXS")}

    # anchor lookup by base compound name (lower)
    anchors["base"] = anchors["compound"].str.split("_").str[0].str.lower()

    cards = []
    for compound, indication, uniprot in DOSSIER_SET:
        mech = D.TARGET_TO_MECHCLASS.get(uniprot, "unknown")
        priors = priors_by_disease.get(indication, {})
        arow = None
        match = anchors[anchors["base"] == compound.lower()]
        if len(match):
            arow = match.iloc[0]
        # indication-matched ledger row (real pivotal outcome) takes priority
        lrow = None
        lmatch = ledger[(ledger["compound"].str.lower() == compound.lower())
                        & ledger["indication"].apply(
                            lambda s: D.disease_match(str(s), indication))]
        if len(lmatch):
            lrow = lmatch.iloc[0]
        card = Dx.build_dossier(
            compound, indication,
            ledger=ledger, disease_priors=priors, anchor_row=arow, ledger_row=lrow,
            liability_df=liability, panel_uniprots=panel_u,
            binding_reliable=(uniprot not in ALLOSTERIC_BLIND),
            target_gene=gene.get(uniprot, uniprot), mechanism_class=mech)
        cards.append(card)

    L = ["# Clinician Evidence Dossiers (Gap 5)", "",
         "One-page, GRADE-style evidence cards — the single artifact a clinician "
         "reads instead of the full report suite. Each distils the pipeline's real "
         "outputs into: predicted cognition effect size + credible interval, "
         "Cochrane-GRADE evidence quality with explicit up/down-grade reasons, the "
         "mechanism-class pivotal-trial track record (the Gap-3 prognostic signal), "
         "predicted off-target liability flags, the provenance trail, and explicit "
         "failure-mode caveats.", "",
         "> Effect sizes are predicted cognition Hedges' g bounded by the Roberts 2020 "
         "ceiling; off-target flags are model-predicted (MAMMAL DTI) and unvalidated. "
         "This is a triage aid, not a prescribing guide.", "",
         "---", ""]
    for card in cards:
        L.append(Dx.render_card_md(card))
        L.append("---")
        L.append("")
    L.append("Generated by `scripts/80_clinician_dossier.py` via "
             "`reporting/clinician_dossier.py`.")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(L), encoding="utf-8")

    logger.info("=" * 64)
    for card in cards:
        logger.info("  %-16s %-11s GRADE=%-9s class=%-8s g=%.2f",
                    card.compound, card.indication, card.grade,
                    card.class_verdict, card.g)
    logger.info("Wrote %s (%d dossiers)", args.report, len(cards))
    logger.info("=" * 64)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
