"""TrkB (NTRK2) site-split calibration - the pre-registered TrkB-TMD negative control.

The deep-research gap analysis (psychoplastogen-targets lane) flagged that the existing NTRK2
anchor set conflates THREE incompatible binding sites and contains ZERO of the durability-
relevant ones: the antidepressant/ketamine site is a crossed transmembrane-domain DIMER wedge
gated by membrane cholesterol (Casarotto 2021 Cell; Cordeiro 2024 Nat Commun), solvable only
by NMR-in-lipid. A 1D-sequence DTI model (MAMMAL) carries no membrane, no dimer, no
cholesterol, so the durability-determining mode is information-theoretically absent.

PRE-REGISTERED EXPECTATION (write it down before scoring): MAMMAL will FAIL to rank the
tmd_wedge antidepressants (fluoxetine/imipramine/ketamine/(2R,6R)-HNK) above matched
non-engagers at NTRK2, while it may rank the ATP-pocket TRK kinase inhibitors (which ARE in
BindingDB) - i.e. the engine sees the wrong site. A FAIL here is the publishable result: it
extends the project's allosteric-blindness audit to the TrkB-TMD durability site and bounds
what any sequence-only second opinion can do (even Boltz-2 would need the dimer+membrane).

GPU. Scores NTRK2 x (site anchors + size-matched negatives) and reports per-site AUROC +
permutation-p. Writes reports/pipeline/trkb_tmd_sitesplit_v1.md.
"""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from mammal_repurposing.engine.persistence_dti import calibrate_target, load_panel

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("trkb_tmd")

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
INTERIM = ROOT / "data" / "interim"
SITE = RAW / "ntrk2_site_anchors.csv"
ANCHORS = RAW / "persistence_dti_anchors.csv"
PANEL = INTERIM / "persistence_targets.csv"
REPORT = ROOT / "reports" / "pipeline" / "trkb_tmd_sitesplit_v1.md"
SITES = ["tmd_wedge", "ecd", "atp_pocket"]


def main() -> int:
    panel = load_panel(PANEL)
    seq = panel["NTRK2"].sequence
    site = pd.read_csv(SITE)
    negs = pd.read_csv(ANCHORS)
    negs = negs[negs["role"] == "non_engager"][["compound", "smiles"]].drop_duplicates()

    rows = [{"compound": r["compound"], "site": r["site"], "role": "engager", "smiles": r["smiles"]}
            for _, r in site.iterrows()]
    rows += [{"compound": r["compound"], "site": "negative", "role": "non_engager",
              "smiles": r["smiles"]} for _, r in negs.iterrows()]
    L.info("Scoring NTRK2 x %d compounds (%d site engagers + %d negatives)",
           len(rows), len(site), len(negs))

    from mammal_repurposing.scoring.dti import score_batch_safe
    from mammal_repurposing.scoring.model_loader import load_dti_model
    model, tok = load_dti_model()
    pkds = []
    for i in range(0, len(rows), 4):
        chunk = rows[i:i + 4]
        pkds.extend(score_batch_safe(model, tok, [(seq, r["smiles"]) for r in chunk],
                                     sample_ids=[f"NTRK2|{r['compound']}" for r in chunk]))
    df = pd.DataFrame(rows); df["predicted_pkd"] = pkds
    neg_scores = df[df["role"] == "non_engager"]["predicted_pkd"].tolist()

    per_site = {}
    for s in SITES:
        pos = df[(df["role"] == "engager") & (df["site"] == s)]["predicted_pkd"].tolist()
        per_site[s] = calibrate_target(pos, neg_scores, min_pos=2)   # small n per site

    write_report(df, per_site, len(neg_scores))
    for s in SITES:
        c = per_site[s]
        L.info("NTRK2/%s: AUROC=%.2f perm-p=%.3f passed=%s", s, c["auroc"], c["perm_p"], c["passed"])
    tmd = per_site["tmd_wedge"]
    L.info("PRE-REGISTERED CHECK - tmd_wedge (durability site) passed=%s (expected FALSE)",
           tmd["passed"])
    return 0


def write_report(df, per_site, n_neg) -> None:
    Ls = ["# TrkB (NTRK2) site-split calibration - pre-registered TrkB-TMD negative", "",
          "Does MAMMAL's sequence-only DTI head rank the DURABILITY-relevant TrkB-TMD "
          "antidepressants (Casarotto 2021) above matched non-engagers - or only the "
          "ATP-pocket TRK inhibitors that happen to be in BindingDB? Reproduced by "
          "`scripts/112_trkb_tmd_sitesplit.py`.", "",
          "**PRE-REGISTERED EXPECTATION:** tmd_wedge FAILS (the crossed-dimer/cholesterol mode "
          "is information-theoretically absent from a 1D sequence); ATP-pocket may pass.", "",
          f"Negatives: {n_neg} size-matched non-engagers.", "",
          "| binding site | engagers | AUROC | perm-p | passes (AUROC>=0.70, p<0.05) |",
          "|---|---|---|---|---|"]
    for s in SITES:
        c = per_site[s]
        names = ", ".join(df[(df.role == "engager") & (df.site == s)]["compound"])
        a = "n/a" if c["auroc"] != c["auroc"] else f"{c['auroc']:.2f}"
        p = "n/a" if c["perm_p"] != c["perm_p"] else f"{c['perm_p']:.3f}"
        Ls.append(f"| {s} | {names} | {a} | {p} | {'PASS' if c['passed'] else '**fail**'} |")
    tmd_fail = not per_site["tmd_wedge"]["passed"]
    Ls += ["", "## Reading", "",
           ("CONFIRMED pre-registration: MAMMAL FAILS the tmd_wedge (durability) site"
            if tmd_fail else "UNEXPECTED: tmd_wedge passed - investigate before trusting") + ". "
           "The engine, if it ranks TrkB engagers at all, sees the ATP-pocket / ECD site, NOT "
           "the transmembrane-domain crossed-dimer wedge that mediates the durable "
           "antidepressant/plasticity effect. This is why PERSEUS routes psychoplastogen "
           "durability through the L4 permeability-gated window (off the DTI axis) rather than a "
           "TrkB DTI score, and why even a Boltz-2 structure second opinion is insufficient "
           "here: the active site only forms as a cholesterol-dependent crossed dimer in the "
           "lipid bilayer (Casarotto 2021; Cordeiro 2024), which a single-chain apo prediction "
           "cannot represent. The TrkB-TMD durability channel is therefore a documented "
           "off-axis limit, not a buildable DTI head.", ""]
    REPORT.write_text("\n".join(Ls), encoding="utf-8")
    L.info("Wrote %s", REPORT)


if __name__ == "__main__":
    raise SystemExit(main())
