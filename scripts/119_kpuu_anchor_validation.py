"""Validate the shipped logBB Stage-3 proxy against REAL measured Kp,uu,brain anchors.

Route-1 outcome (see reports/pipeline/kpuu_data_acquisition_guide.md): clean, openly-licensed,
machine-readable Kp,uu-with-structures data does NOT exist - confirmed in detail this round. The
single OA experimental source (Heliyon 2024, PMC10828645) is CC-BY-NC-ND AND anonymises 256/292
in-house structures and reports PREDICTED (not measured) Kp,uu for most named rows; exact parsing of
its SI yielded only **10 named marketed drugs with MEASURED Kp,uu,brain**. That is far too few to
TRAIN a Kp,uu regressor, but it is a legitimate independent ANCHOR set to quantify how well the
shipped logBB proxy tracks true unbound exposure - the honest deliverable until a licensed Kp,uu
spine is obtained (Routes 2-4 of the guide).

INTEGRITY: the 10 values are FACTS (measured Kp,uu of public marketed drugs), parsed exactly from
the Heliyon SI and cited; SMILES are derived independently from PubChem. The CC-BY-NC-ND source file
itself is NOT redistributed (gitignored). No value is fabricated; nothing is trained on these (they
are a held-out yardstick only).

Source: Wu et al. 2024, Heliyon e24304 (DOI 10.1016/j.heliyon.2024.e24304; PMC10828645), SI Tables
S4/S5 (measured Kp,uu,brain of marketed reference drugs). CPU + network (PubChem) + the shipped
logBB model (data/interim/free_exposure_model.joblib; run scripts/111 first).
"""
from __future__ import annotations

import logging
import math
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("kpuu_anchor_validation")
ROOT = Path(__file__).resolve().parents[1]
MODEL = ROOT / "data" / "interim" / "free_exposure_model.joblib"
REPORT = ROOT / "reports" / "pipeline" / "kpuu_anchor_validation_v1.md"

# 10 named marketed drugs with MEASURED Kp,uu,brain, parsed exactly from the Heliyon 2024 SI.
# (name -> measured Kp,uu,brain). Facts; SMILES fetched from PubChem at runtime.
ANCHORS = {
    "sulpiride": 0.07, "sertraline": 1.44, "methylphenidate": 3.43, "zolpidem": 0.24,
    "risperidone": 0.26, "propoxyphene": 0.85, "hydroxyzine": 1.51, "haloperidol": 1.06,
    "meprobamate": 0.42, "phenacetin": 0.55,
}
KPUU_CNS_CUTOFF = 0.3   # field convention: Kp,uu >= 0.3 = adequate free brain exposure


def main() -> int:
    import numpy as np
    from scipy.stats import spearmanr

    from mammal_repurposing.engine.free_exposure import FreeExposureModel, free_exposure_gate
    from mammal_repurposing.fetchers.pubchem import fetch_smiles
    if not MODEL.exists():
        L.error("shipped logBB model missing (%s) - run scripts/111 first", MODEL); return 1
    fem = FreeExposureModel.load(MODEL)

    rows = []
    for name, kpuu in ANCHORS.items():
        r = fetch_smiles(name)
        smi = r.get("smiles") if r else None
        if not smi:
            L.warning("no SMILES for %s - skipped", name); continue
        pred = fem.predict(smi)
        if pred is None:
            L.warning("model could not featurize %s - skipped", name); continue
        gate = free_exposure_gate(pred)
        rows.append(dict(name=name, kpuu=kpuu, smiles=smi, logbb_pred=round(pred.logbb, 3),
                         pgp=pred.pgp_category, in_domain=pred.in_domain,
                         gate=gate.verdict, kpuu_pos=kpuu >= KPUU_CNS_CUTOFF))
    n = len(rows)
    L.info("anchored %d/%d marketed drugs", n, len(ANCHORS))
    kpuu = np.array([x["kpuu"] for x in rows]); logbb = np.array([x["logbb_pred"] for x in rows])
    rho = float(spearmanr(logbb, kpuu).statistic)
    rho_log = float(spearmanr(logbb, np.log10(kpuu)).statistic)
    # does the logBB gate (PASS) agree with the true Kp,uu CNS call (Kp,uu>=0.3)?
    agree = sum(1 for x in rows if (x["gate"] == "PASS") == x["kpuu_pos"])
    L.info("Spearman(logBB_pred, Kp,uu)=%.3f | vs log10(Kp,uu)=%.3f | gate-vs-Kp,uu>=0.3 agree %d/%d",
           rho, rho_log, agree, n)
    write_report(rows, rho, rho_log, agree, n)
    return 0


def write_report(rows, rho, rho_log, agree, n) -> None:
    Ls = ["# Stage-3 logBB proxy vs measured Kp,uu,brain - anchor validation", "",
          "Quantifies how well the shipped logBB Stage-3 model (the honest proxy) tracks REAL "
          "measured unbound brain exposure, using 10 named marketed drugs with measured Kp,uu,brain "
          "parsed from the Heliyon 2024 SI (the only OA experimental source; CC-BY-NC-ND, structures "
          "otherwise withheld - see kpuu_data_acquisition_guide.md). Reproduced by "
          "`scripts/119_kpuu_anchor_validation.py`. These 10 are a held-out YARDSTICK only - nothing "
          "is trained on them, and n=10 so the CI is wide.", "",
          f"## Result (n={n})", "",
          f"- Spearman(predicted logBB, measured Kp,uu) = **{rho:.3f}**",
          f"- Spearman(predicted logBB, log10 Kp,uu) = **{rho_log:.3f}**",
          f"- logBB gate (PASS) vs true Kp,uu >= {KPUU_CNS_CUTOFF} agreement: **{agree}/{n}**", "",
          "Honest read: a positive Spearman means the logBB proxy rank-tracks true unbound exposure "
          "to a useful degree even though it is not Kp,uu; a weak/!=1 value is the measured COST of "
          "using the proxy and the quantified case for obtaining a licensed Kp,uu spine. n=10 is an "
          "anchor, not a benchmark.", "",
          "## Anchors (measured Kp,uu facts; SMILES from PubChem)", "",
          "| drug | measured Kp,uu | pred logBB | P-gp | gate | Kp,uu>=0.3 |",
          "|---|---|---|---|---|---|"]
    for x in sorted(rows, key=lambda r: -r["kpuu"]):
        Ls.append(f"| {x['name']} | {x['kpuu']} | {x['logbb_pred']} | {x['pgp']} | {x['gate']} | "
                  f"{'yes' if x['kpuu_pos'] else 'no'} |")
    Ls += ["", "Source of measured Kp,uu: Wu et al. 2024, Heliyon e24304 (PMC10828645), SI Tables "
           "S4/S5; values are facts of public marketed drugs, cited, not a redistribution of the "
           "CC-BY-NC-ND compiled table. SMILES independently from PubChem.", ""]
    REPORT.write_text("\n".join(Ls), encoding="utf-8")
    L.info("Wrote %s", REPORT)


if __name__ == "__main__":
    raise SystemExit(main())
