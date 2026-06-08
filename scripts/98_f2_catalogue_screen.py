"""F2 capstone - catalogue-scale prospective screen.

Run the F2 novel-compound engine over the full ChEMBL approved-drug catalogue
(~3.4k drugs across every therapeutic area) and surface the actionable output:
approved drugs that are structurally members of a STRONG-PRECEDENT cognition
mechanism class (catecholaminergic_ADHD / wake_promoting / AChE_inhibitor - the
classes whose clinical-g prior is positive and whose success rate >= 0.5) yet are
NOT in our cognition-outcome ledger. Each is a structure-grounded, prior-quantified
repurposing hypothesis: "this drug looks like a [class] member, and that class has
clinical cognition precedent - has it been tried for cognition?".

Honest scope. (1) Routing is by 2D structure (F2's validated signal; the GPU
DTI-profile was a tested negative). (2) "Not in our ledger" is not proof a drug was
never trialled for cognition - this is a HYPOTHESIS-GENERATION screen; prior-trial
verification (the trial-watch system) is the follow-up. (3) The predicted g is the
class prior, a model output, never a measured outcome. Nothing here is fabricated.

CPU only (RDKit + numpy/pandas). Reads data/raw/chembl_approved_catalogue.csv
(from scripts/_fetch_chembl_approved.py). Writes the shortlist CSV + report.
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from mammal_repurposing.reporting.trial_watch import _norm_drug, load_combined_ledger
from mammal_repurposing.validation.novel_compound import (
    MIN_CLASS_N, build_class_priors, build_exemplars, score_catalogue,
)

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("f2_cat")

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
LEDGERS = [RAW / "clinical_outcomes_ledger.csv",
           RAW / "clinical_outcomes_ledger_EXTENSION.csv",
           RAW / "clinical_outcomes_ledger_CTGOV.csv",
           RAW / "clinical_outcomes_ledger_RESEARCH.csv"]
SMILES = RAW / "ledger_compound_smiles.csv"
CATALOGUE = RAW / "chembl_approved_catalogue.csv"
REPORT = ROOT / "reports" / "pipeline" / "f2_catalogue_shortlist_v1.md"
SHORTLIST = ROOT / "reports" / "pipeline" / "f2_catalogue_shortlist.csv"


def main() -> int:
    led = load_combined_ledger(LEDGERS)
    smi = pd.read_csv(SMILES)[["compound", "smiles"]]
    ex = build_exemplars(led, smi)
    priors = build_class_priors(led, n_boot=2000, seed=0)

    # strong-precedent classes: positive prior g, success rate >= 0.5, real prior
    strong = {c for c, p in priors.items()
              if p.p_success >= 0.5 and p.prior_g > 0 and p.n >= MIN_CLASS_N}
    L.info("strong-precedent classes: %s", sorted(strong))

    cat = pd.read_csv(CATALOGUE)
    cat["name"] = cat["name"].astype(str)
    # collapse salt/combination forms: one row per unique PARENT structure
    # (the fetcher already reduced each to its largest organic fragment), keeping
    # the shortest name (usually the base drug, e.g. AMPHETAMINE over ... SULFATE).
    n_raw = len(cat)
    cat = (cat.assign(_l=cat["name"].str.len()).sort_values("_l")
           .drop_duplicates("smiles", keep="first").drop(columns="_l"))
    L.info("catalogue %d rows -> %d unique parent structures", n_raw, len(cat))
    # light, uniform druglikeness gate: drop tiny/excipient-like molecules (heavy
    # atoms < 12, ~MW < 150) such as benzyl alcohol that match a class only by
    # chemical coincidence. Applied to all compounds equally (not cherry-picking).
    from rdkit import Chem  # noqa: PLC0415

    def _heavy(s):
        m = Chem.MolFromSmiles(str(s))
        return m.GetNumHeavyAtoms() if m is not None else 0

    n_pre = len(cat)
    cat = cat[cat["smiles"].map(_heavy) >= 12].copy()
    L.info("druglikeness gate (>=12 heavy atoms): %d -> %d", n_pre, len(cat))
    # screen only drugs NOT already in our cognition-outcome ledger
    in_ledger = set(led["compound_norm"])
    cat["_norm"] = cat["name"].map(_norm_drug)
    cat_novel = cat[~cat["_norm"].isin(in_ledger)].drop_duplicates("_norm")
    L.info("%d unique structures -> %d not in cognition ledger", len(cat), len(cat_novel))

    scored = score_catalogue(cat_novel.rename(columns={"name": "id"})[["id", "smiles"]],
                             ex, priors)
    scored = scored.merge(cat_novel[["name", "chembl_id"]].rename(columns={"name": "query_id"}),
                          on="query_id", how="left")

    n_routed = int((scored["tier"] != "ABSTAIN").sum())
    # the discovery shortlist
    short = scored[(scored["tier"].isin(["HIGH", "MED"]))
                   & (scored["assigned_class"].isin(strong))
                   & (scored["predicted_outcome"] == "SUCCESS")].copy()
    short = short.sort_values(["prior_g", "similarity"], ascending=False).reset_index(drop=True)
    SHORTLIST.parent.mkdir(parents=True, exist_ok=True)
    short[["query_id", "chembl_id", "assigned_class", "tier", "similarity",
           "prior_g", "g_ci_lo", "g_ci_hi", "p_success", "scaffold_hit",
           "reason"]].to_csv(SHORTLIST, index=False)

    # per-class breakdown of the shortlist
    by_class = short.groupby("assigned_class").size().sort_values(ascending=False)

    # ---- report ----
    Ls = ["# F2 capstone - catalogue-scale repurposing screen", "",
          "**The actionable output of F2.** Every approved drug in ChEMBL (max_phase=4) "
          "routed through the novel-compound engine; the shortlist is the approved drugs "
          "that are structurally members of a STRONG-PRECEDENT cognition class but are not "
          "in our cognition-outcome ledger. Reproduced by `scripts/_fetch_chembl_approved.py` "
          "+ `scripts/98_f2_catalogue_screen.py`.", "",
          f"- Catalogue: **{n_raw}** approved-drug rows (ChEMBL max_phase=4, RDKit-parsed) "
          f"-> **{len(cat)}** unique drug-like parent structures (salts/combinations "
          f"collapsed, >=12 heavy atoms). Removing {len(cat) - len(cat_novel)} already in "
          f"the cognition ledger leaves **{len(cat_novel)}** screened.",
          f"- Routed (not abstained): **{n_routed}** ({100*n_routed/len(cat_novel):.0f}%); "
          f"the rest are out-of-manifold for the cognition exemplars and correctly abstain.",
          f"- Strong-precedent classes (prior g>0, success>=0.5, n>={MIN_CLASS_N}): "
          f"**{', '.join(sorted(strong))}**.",
          f"- **Shortlist: {len(short)} repurposing hypotheses** "
          f"(HIGH/MED, predicted SUCCESS). Full CSV: `reports/pipeline/f2_catalogue_shortlist.csv`.",
          ""]
    Ls.append("Per class: " + ", ".join(f"{c} ({n})" for c, n in by_class.items()) + ".")
    Ls.append("")
    Ls.append(f"## Top {min(40, len(short))} hypotheses (by class prior g)")
    Ls.append("")
    Ls.append("| drug | class | tier | sim | predicted g [90% CrI] | P(success) | scaffold |")
    Ls.append("|---|---|---|---|---|---|---|")
    for _, r in short.head(40).iterrows():
        g = (f"{r['prior_g']:+.2f} [{r['g_ci_lo']:+.2f}, {r['g_ci_hi']:+.2f}]"
             if np.isfinite(r["prior_g"]) else "-")
        sc = "yes" if r["scaffold_hit"] else "-"
        Ls.append(f"| {r['query_id']} | {r['assigned_class']} | {r['tier']} | "
                  f"{r['similarity']:.2f} | {g} | {r['p_success']:.2f} | {sc} |")
    Ls.append("")
    Ls.append("## Honest scope")
    Ls.append("")
    Ls.append("- Structure-based routing (F2's validated signal; the MAMMAL DTI-profile "
              "was a tested negative, `f2_profile_vs_structure_v1.md`).")
    Ls.append("- \"Not in our ledger\" is NOT proof a drug was never trialled for cognition. "
              "This is hypothesis generation; each hit needs prior-trial verification "
              "(the trial-watch system) before it is a genuine novel-repurposing claim.")
    Ls.append("- The predicted g is the assigned class's prior - a model output, not a "
              "measured outcome. Many hits will be near-analogs of the class exemplars "
              "(e.g. other sympathomimetics); the value is the ranked, prior-quantified, "
              "structure-grounded surface for triage.")
    Ls.append("")

    REPORT.write_text("\n".join(Ls), encoding="utf-8")
    L.info("Wrote %s and %s", REPORT, SHORTLIST)
    L.info("F2 capstone: %d screened -> %d routed -> %d shortlist hypotheses across %d classes",
           len(cat_novel), n_routed, len(short), len(by_class))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
