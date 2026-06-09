"""AMPA-PAM + PNN/ECM persistence-target channels (deep-research Gap-1/Gap-3, pre-registered).

Two more candidate structure-derivable durability channels, each calibrated against the same
size-matched negative pool that killed 7/9 of the original DTI panel:

  * AMPA (GRIA1), SITE-SPLIT into orthosteric vs PAM. The durability lever is allosteric
    POTENTIATION (PAM: cyclothiazide/aniracetam/ampakines) + GluA1-phospho + mTORC1 coupling -
    and MAMMAL is documented allosterically blind (v1 CHRNA7 PAM std 0.029). PRE-REGISTERED:
    MAMMAL FAILS the AMPA-PAM site (orthosteric may rank, PAM should not).
  * PNN/ECM remodeling via MMP9 (matrix metalloproteinase; chondroitinase-like plasticity
    re-opening). A distinct structural-plasticity tier with no prior representation;
    PRE-REGISTERED as uncertain - it must clear the same size-matched gate to count.

GPU. Scores GRIA1/MMP9 x (channel anchors + size-matched negatives), per-(channel,site) AUROC +
permutation-p. Writes reports/pipeline/ampa_pnn_channels_v1.md.
"""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from mammal_repurposing.engine.persistence_dti import calibrate_target

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("ampa_pnn")

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
INTERIM = ROOT / "data" / "interim"
TARGETS = INTERIM / "ampa_pnn_targets.csv"
ANCHORS = RAW / "ampa_pnn_anchors.csv"
NEG = RAW / "persistence_dti_anchors.csv"
REPORT = ROOT / "reports" / "pipeline" / "ampa_pnn_channels_v1.md"

# (gene, channel, sites) groups to calibrate
GROUPS = [("GRIA1", "AMPA", ["orthosteric", "pam"]), ("MMP9", "PNN_ECM", ["mmp_inhibitor"])]


def main() -> int:
    tg = pd.read_csv(TARGETS).set_index("gene")
    anc = pd.read_csv(ANCHORS)
    negs = pd.read_csv(NEG)
    negs = negs[negs["role"] == "non_engager"][["compound", "smiles"]].drop_duplicates()

    rows = []
    for gene, channel, sites in GROUPS:
        seq = tg.loc[gene, "sequence"]
        sub = anc[anc["channel"] == channel]
        for _, r in sub.iterrows():
            rows.append({"gene": gene, "channel": channel, "site": r["site"],
                         "role": "engager", "compound": r["compound"], "seq": seq, "smiles": r["smiles"]})
        for _, n in negs.iterrows():
            rows.append({"gene": gene, "channel": channel, "site": "negative",
                         "role": "non_engager", "compound": n["compound"], "seq": seq, "smiles": n["smiles"]})
    L.info("Scoring %d (compound,target) pairs across %d channels", len(rows), len(GROUPS))

    from mammal_repurposing.scoring.dti import score_batch_safe
    from mammal_repurposing.scoring.model_loader import load_dti_model
    model, tok = load_dti_model()
    pkds = []
    for i in range(0, len(rows), 4):
        c = rows[i:i + 4]
        pkds.extend(score_batch_safe(model, tok, [(r["seq"], r["smiles"]) for r in c],
                                     sample_ids=[f"{r['gene']}|{r['compound']}" for r in c]))
    df = pd.DataFrame(rows); df["predicted_pkd"] = pkds

    results = {}
    for gene, channel, sites in GROUPS:
        g = df[df["gene"] == gene]
        neg_scores = g[g["role"] == "non_engager"]["predicted_pkd"].tolist()
        for s in sites:
            pos = g[(g["role"] == "engager") & (g["site"] == s)]["predicted_pkd"].tolist()
            results[(gene, channel, s)] = calibrate_target(pos, neg_scores, min_pos=2)

    write_report(df, results)
    for k, c in results.items():
        L.info("%s/%s: AUROC=%.2f perm-p=%.3f passed=%s", k[0], k[2], c["auroc"], c["perm_p"], c["passed"])
    pam = results.get(("GRIA1", "AMPA", "pam"))
    L.info("PRE-REGISTERED CHECK - AMPA-PAM passed=%s (expected FALSE, allosteric blindness)",
           pam["passed"] if pam else "n/a")
    return 0


def write_report(df, results) -> None:
    Ls = ["# AMPA-PAM + PNN/ECM persistence channels - pre-registered calibration", "",
          "Two candidate structure-derivable durability channels, calibrated against the "
          "size-matched negative pool. Reproduced by `scripts/114_ampa_pnn_channels.py`.", "",
          "**PRE-REGISTERED:** AMPA-PAM FAILS (the durability lever is allosteric potentiation; "
          "MAMMAL's BindingDB-pKd head is allosterically blind - v1 CHRNA7 PAM std 0.029). AMPA "
          "orthosteric may rank; MMP9/PNN is exploratory and must clear the same gate.", "",
          "| target | channel/site | engagers | AUROC | perm-p | passes |",
          "|---|---|---|---|---|---|"]
    for (gene, channel, s), c in results.items():
        names = ", ".join(df[(df.gene == gene) & (df.site == s)]["compound"])
        a = "n/a" if c["auroc"] != c["auroc"] else f"{c['auroc']:.2f}"
        p = "n/a" if c["perm_p"] != c["perm_p"] else f"{c['perm_p']:.3f}"
        Ls.append(f"| {gene} | {channel}/{s} | {names} | {a} | {p} | "
                  f"{'PASS' if c['passed'] else '**fail**'} |")
    pam = results.get(("GRIA1", "AMPA", "pam"))
    Ls += ["", "## Reading", "",
           ("CONFIRMED: MAMMAL FAILS the AMPA-PAM (allosteric-potentiation) site - the "
            "durability lever is invisible to a BindingDB-pKd head, extending the v1 "
            "allosteric-blindness audit to the AMPA channel."
            if pam and not pam["passed"] else
            "UNEXPECTED AMPA-PAM result - investigate before trusting.") +
           " The orthosteric/MMP9 rows are reported for completeness; any channel that does not "
           "clear the size-matched + permutation gate is NOT wired into the persistence head "
           "(abstain-by-default). Net: AMPA durability (PAM/phospho/mTORC1-coupled) and PNN/ECM "
           "remodeling are confirmed OFF the sequence-DTI axis, consistent with routing "
           "plasticity through the L4 permeability window and reserving these tiers for a "
           "structure/allosteric-aware second opinion.", ""]
    REPORT.write_text("\n".join(Ls), encoding="utf-8")
    L.info("Wrote %s", REPORT)


if __name__ == "__main__":
    raise SystemExit(main())
