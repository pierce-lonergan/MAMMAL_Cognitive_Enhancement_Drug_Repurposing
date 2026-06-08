"""Stage 2 of the persistence-target DTI module (PERSEUS roadmap #2): the CALIBRATION run.

THE question this answers: can MAMMAL's DTI head actually rank known engagers of each
persistence-substrate target above matched non-engagers? Only targets that pass may later
contribute a substrate read for a novel compound. A target that fails is dropped - and a
panel that mostly fails is itself an honest, publishable result (it would justify PERSEUS's
abstain-by-default persistence head and the project's documented finding that MAMMAL is a
weak class router).

GPU. Scores anchors x panel via the shared DTI wrapper, then per-target AUROC +
permutation-p + Youden threshold (engine.persistence_dti). Writes:
  - data/results/persistence_dti_scores.csv
  - data/results/persistence_dti_calibration.json
  - reports/pipeline/persistence_dti_calibration.md

Run with the MAMMAL venv:
  .venv-mammal/Scripts/python.exe scripts/104_persistence_dti_calibrate.py
"""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import pandas as pd

from mammal_repurposing.engine.persistence_dti import (
    DEFAULT_MIN_AUROC, DEFAULT_MIN_POS, DEFAULT_PERM_P, calibrate_target, load_panel,
)

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("persistence_dti_calib")

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
INTERIM = ROOT / "data" / "interim"
RESULTS = ROOT / "data" / "results"
TARGETS = INTERIM / "persistence_targets.csv"
ANCHORS = RAW / "persistence_dti_anchors.csv"
SCORES_OUT = RESULTS / "persistence_dti_scores.csv"
CALIB_OUT = RESULTS / "persistence_dti_calibration.json"
REPORT = ROOT / "reports" / "pipeline" / "persistence_dti_calibration.md"


def build_grid(panel, anchors: pd.DataFrame) -> list[dict]:
    """Per target T: positives = engagers of T; negatives = the shared non-engager pool."""
    negs = anchors[anchors["role"] == "non_engager"][["compound", "smiles"]].drop_duplicates()
    rows = []
    for gene, tgt in panel.items():
        for _, e in anchors[(anchors["role"] == "engager") & (anchors["target_gene"] == gene)].iterrows():
            rows.append({"target_gene": gene, "compound": e["compound"],
                         "role": "engager", "smiles": e["smiles"], "seq": tgt.sequence})
        for _, n in negs.iterrows():
            rows.append({"target_gene": gene, "compound": n["compound"],
                         "role": "non_engager", "smiles": n["smiles"], "seq": tgt.sequence})
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch-size", type=int, default=4,
                    help="small default: panel has 1200-1600 AA sequences (DNMT1/HDAC6/EHMT2)")
    args = ap.parse_args()

    panel = load_panel(TARGETS)
    anchors = pd.read_csv(ANCHORS)
    grid = build_grid(panel, anchors)
    L.info("Scoring grid: %d (compound,target) pairs across %d targets", len(grid), len(panel))

    from mammal_repurposing.scoring.dti import score_batch_safe
    from mammal_repurposing.scoring.model_loader import load_dti_model
    model, tok = load_dti_model()

    pkds: list[float] = []
    for i in range(0, len(grid), args.batch_size):
        chunk = grid[i:i + args.batch_size]
        pairs = [(r["seq"], r["smiles"]) for r in chunk]
        ids = [f"{r['target_gene']}|{r['compound']}" for r in chunk]
        pkds.extend(score_batch_safe(model, tok, pairs, sample_ids=ids))
        L.info("  scored %d/%d", min(i + args.batch_size, len(grid)), len(grid))

    scored = pd.DataFrame(grid)
    scored["predicted_pkd"] = pkds
    scored = scored.drop(columns="seq")
    mw = dict(zip(anchors["compound"], anchors["mw"])) if "mw" in anchors else {}
    scored["mw"] = scored["compound"].map(mw)
    RESULTS.mkdir(parents=True, exist_ok=True)
    scored.to_csv(SCORES_OUT, index=False)
    L.info("Wrote %s", SCORES_OUT)

    # molecular-size confound: how much of the predicted pKd is just molecular weight?
    # measured over the NON-engager pool (which now spans 129-671 Da by design).
    import numpy as np
    def _confound(df):
        d = df.dropna(subset=["predicted_pkd", "mw"])
        if len(d) < 3 or d["mw"].nunique() < 2:
            return None
        return round(float(np.corrcoef(d["mw"], d["predicted_pkd"])[0, 1]), 3)
    neg_all = scored[scored["role"] == "non_engager"]
    size_confound_overall = _confound(neg_all)
    neg_mw = neg_all["mw"].dropna()
    mw_lo, mw_hi = (int(neg_mw.min()), int(neg_mw.max())) if len(neg_mw) else (0, 0)

    per_target, n_nan = {}, int(scored["predicted_pkd"].isna().sum())
    for gene, tgt in panel.items():
        sub = scored[scored["target_gene"] == gene]
        pos = sub[sub["role"] == "engager"]["predicted_pkd"].tolist()
        neg = sub[sub["role"] == "non_engager"]["predicted_pkd"].tolist()
        cal = calibrate_target(pos, neg)
        cal.update({"tier": tgt.tier, "promotes_durable": tgt.promotes_durable,
                    "n_pos_total": len(pos),
                    "size_confound_r": _confound(sub[sub["role"] == "non_engager"])})
        per_target[gene] = cal

    passed = [g for g, c in per_target.items() if c["passed"]]
    calib = {
        "meta": {"min_auroc": DEFAULT_MIN_AUROC, "min_pos": DEFAULT_MIN_POS,
                 "max_perm_p": DEFAULT_PERM_P, "n_perm": 2000},
        "per_target": per_target,
        "summary": {"n_targets": len(panel), "n_passed": len(passed), "passed_genes": passed,
                    "n_pairs": len(grid), "n_nan_pairs": n_nan,
                    "size_confound_r": size_confound_overall,
                    "n_non_engagers": int(neg_all["compound"].nunique()),
                    "non_engager_mw_lo": mw_lo, "non_engager_mw_hi": mw_hi},
    }
    CALIB_OUT.write_text(json.dumps(calib, indent=2), encoding="utf-8")
    L.info("Wrote %s", CALIB_OUT)

    write_report(panel, per_target, calib["summary"])
    L.info("CALIBRATION: %d/%d targets pass (AUROC>=%.2f & perm-p<%.2f & n_pos>=%d): %s",
           len(passed), len(panel), DEFAULT_MIN_AUROC, DEFAULT_PERM_P, DEFAULT_MIN_POS, passed)
    return 0


def write_report(panel, per_target, summary) -> None:
    n_pass = summary["n_passed"]
    Ls = ["# Persistence-target DTI module - calibration", "",
          "Can MAMMAL's DTI head rank known ENGAGERS of each persistence-substrate target "
          "above matched non-engagers? A target may contribute a substrate read for a novel "
          "compound ONLY if it passes here. Reproduced by "
          "`scripts/104_persistence_dti_calibrate.py`.", "",
          f"**Headline: {n_pass}/{summary['n_targets']} targets pass** "
          f"(AUROC>=0.70 AND permutation-p<0.05 AND >=3 scored engagers). "
          f"Passing targets: {', '.join(summary['passed_genes']) or 'NONE'}.", "",
          f"Scored {summary['n_pairs']} (compound,target) pairs"
          + (f"; {summary['n_nan_pairs']} returned NaN (tokenizer/forward overflow on long "
             "sequences)." if summary["n_nan_pairs"] else "; no NaN pairs."), "",
          f"**Molecular-size confound:** corr(MW, predicted pKd) over the "
          f"{summary.get('n_non_engagers', '?')} size-matched non-engagers "
          f"({summary.get('non_engager_mw_lo', '?')}-{summary.get('non_engager_mw_hi', '?')} Da) "
          f"= **{summary.get('size_confound_r', 'n/a')}**. MAMMAL's pKd is "
          "substantially molecular-weight-driven; the negative pool is SIZE-MATCHED on "
          "purpose so a channel cannot pass just by scoring big molecules high.", "",
          "| target | tier | durable? | n engagers (scored/total) | AUROC | perm-p | "
          "sens@thr | size-r | PASS | usable |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for gene, c in sorted(per_target.items(), key=lambda kv: -(kv[1]["auroc"] if kv[1]["auroc"] == kv[1]["auroc"] else -1)):
        a = "n/a" if c["auroc"] != c["auroc"] else f"{c['auroc']:.2f}"
        p = "n/a" if c["perm_p"] != c["perm_p"] else f"{c['perm_p']:.3f}"
        sens = c.get("sensitivity_at_threshold")
        sens_s = "n/a" if sens is None or sens != sens else f"{sens:.2f}"
        sr = c.get("size_confound_r")
        sr_s = "n/a" if sr is None else f"{sr:.2f}"
        Ls.append(f"| {gene} | {c['tier']} | {'yes' if c['promotes_durable'] else 'no'} | "
                  f"{c['n_pos_scored']}/{c['n_pos_total']} | {a} | {p} | {sens_s} | {sr_s} | "
                  f"{'**PASS**' if c['passed'] else 'fail'} | "
                  f"{'yes' if c.get('inference_usable') else 'no'} |")
    Ls += ["", "## Interpretation", ""]
    if n_pass == 0:
        Ls.append("MAMMAL does NOT reliably route engagement of ANY persistence-substrate "
                  "target. This is an honest negative: it confirms the project's finding that "
                  "MAMMAL is a weak class router and JUSTIFIES PERSEUS's abstain-by-default "
                  "persistence head - the structure-computable substrate channel is not yet "
                  "trustworthy, so the engine correctly withholds rather than guesses.")
    else:
        Ls.append(f"{n_pass} target(s) pass calibration and may contribute a substrate read. "
                  "Only ablative (senolytic) passes promote toward DURABLE; capability/window "
                  "passes are reported as hypotheses (capability flags), never auto-promoted. "
                  "Every contribution is gated on the per-target threshold + PASS flag above, "
                  "so an un-calibrated channel is ignored at inference, not trusted.")
    Ls += ["", "Calibration is consumed by `engine/persistence_dti.py:substrate_hypothesis` "
           "(`data/results/persistence_dti_calibration.json`).", ""]
    REPORT.write_text("\n".join(Ls), encoding="utf-8")
    L.info("Wrote %s", REPORT)


if __name__ == "__main__":
    raise SystemExit(main())
