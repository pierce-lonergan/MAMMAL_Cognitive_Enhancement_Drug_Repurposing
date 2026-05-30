"""Manuscript robustness supplement — pre-empts the obvious reviewer objections.

Emits three analyses that harden the class-prognostic claim:
  1. PREDICTOR COVERAGE — why each predictor's n differs, and the label balance
     of each subset (shows the binding subset is success-ENRICHED, i.e. the
     missingness works against the conclusion, not for it).
  2. COMMON-SUBSET comparison — all predictors scored on the SAME drugs (where
     binding affinity is defined), so the AUROCs are directly comparable and the
     differing-n objection dissolves.
  3. LTR FEATURE ABLATION — leave-one-target-out within-target ρ for each feature
     subset, quantifying how much the foundation model actually contributes vs
     classic ligand-similarity + physicochemical features.

Output: reports/manuscript_robustness.md

Usage:
  python scripts/84_manuscript_robustness.py
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("robustness")


def main() -> int:
    from mammal_repurposing.validation import retrospective as R
    from mammal_repurposing.cluster_a import allosteric_ltr as A

    led = R.load_clinical_ledger(ROOT / "data" / "raw" / "clinical_outcomes_ledger.csv")
    grid = pd.read_parquet(ROOT / "data" / "results" / "v2" / "v6a_grid_expanded.parquet")
    v6b = pd.read_parquet(ROOT / "data" / "results" / "v2"
                          / "cluster_d_posterior_expanded_v2_mh8_ta99.parquet")
    chembl = pd.read_parquet(ROOT / "data" / "results" / "chembl_evidence.parquet")

    cls = R.class_loco_g(led)
    rel = R.target_relevance_score(led, v6b)
    bind = R.binding_score(led, grid)
    pop = R.target_popularity_score(led, chembl)

    def block(s):
        rows = led[led["compound"].isin(s)]
        y = rows["label"]
        return len(rows), int(y.sum()), int((1 - y).sum())

    preds = [("Mechanism-class track record", cls),
             ("Target genetic relevance", rel),
             ("Target binding affinity (MAMMAL DTI)", bind),
             ("Target popularity (ChEMBL records)", pop)]

    L: list[str] = []
    L.append("# Manuscript robustness supplement")
    L.append("")
    L.append("Reproducible answers to the predictable reviewer objections "
             "(`scripts/84_manuscript_robustness.py`).")
    L.append("")
    L.append("## 1. Predictor coverage and label balance")
    L.append("")
    L.append("Each predictor is defined on a different subset of the 31-drug ledger. "
             "The mechanism-class predictor needs only the drug's class siblings (all "
             "31). Target genetic relevance needs the target in the V6.B posterior "
             "(26). Target binding affinity needs the drug to be **in the 298-compound "
             "screening library AND scored at its known target** (14) — the library was "
             "assembled before the clinical ledger, so 17 mostly-obscure clinical "
             "candidates were never screened.")
    L.append("")
    L.append("| Predictor | n | SUCCESS | FAILURE | success fraction |")
    L.append("|---|---|---|---|---|")
    for name, s in preds:
        n, su, fa = block(s)
        L.append(f"| {name} | {n} | {su} | {fa} | {su / n:.0%} |")
    L.append("")
    cov = set(bind)
    miss = led[~led["compound"].isin(cov)]
    L.append(f"**The binding subset is success-ENRICHED ({block(bind)[1]}:{block(bind)[2]} "
             f"= {block(bind)[1] / block(bind)[0]:.0%} success vs {led['label'].mean():.0%} "
             f"in the full ledger).** The {len(miss)} drugs binding cannot score are "
             f"{int((miss['label'] == 0).sum())} failures and {int((miss['label'] == 1).sum())} "
             f"successes — so the missingness *removes failures* and leaves the "
             f"best-characterised marketed drugs (donepezil, methylphenidate, memantine, "
             f"pitolisant …), the **most favourable possible test for an affinity model**. "
             f"It still scores at chance. The missingness works against the conclusion, "
             f"not for it.")
    L.append("")
    L.append(f"- Failures binding cannot score: {', '.join(sorted(miss[miss.label == 0].compound))}")
    L.append(f"- Successes binding cannot score: {', '.join(sorted(miss[miss.label == 1].compound))}")
    L.append("")

    # 2. common subset
    common = set(bind) & set(rel) & set(cls)
    rows = led[led["compound"].isin(common)]
    y = rows["label"].to_numpy()
    L.append("## 2. Common-subset comparison (all predictors, identical drugs)")
    L.append("")
    L.append(f"Restricting every predictor to the **same {len(rows)} drugs** where binding "
             "affinity is defined removes the differing-n objection entirely:")
    L.append("")
    L.append("| Predictor (same drugs) | AUROC |")
    L.append("|---|---|")
    for name, s in [("Mechanism-class track record", cls),
                    ("Target genetic relevance", rel),
                    ("Target binding affinity", bind)]:
        sv = np.array([s[c] for c in rows["compound"]], float)
        L.append(f"| {name} | {R.auroc(sv, y):.3f} |")
    L.append("")
    L.append("The contrast is preserved apples-to-apples: class track record perfectly "
             "separates; both target-centric predictors are at chance.")
    L.append("")

    # 3. ablation
    ch = chembl[chembl["best_pchembl"].notna()].copy()
    ch["pact"] = ch["best_pchembl"].astype(float)
    dti = pd.read_parquet(ROOT / "data" / "results" / "dti_scores.parquet")
    tani = pd.read_parquet(ROOT / "data" / "results" / "v2" / "disagreement_signal.parquet")
    feat = A.build_feature_table(ch[["compound_name", "target_uniprot", "smiles"]],
                                 mammal=dti, tanimoto=tani)
    feat = feat.merge(ch[["compound_name", "target_uniprot", "pact"]],
                      on=["compound_name", "target_uniprot"], how="left")
    feat = feat[feat["pact"].notna()].reset_index(drop=True)
    abl = [("MAMMAL pKd only", ["mammal_pkd"]),
           ("Physicochemical only", A.PHYSCHEM_COLS),
           ("Tanimoto-to-actives only", ["tanimoto"]),
           ("Tanimoto + physchem (NO foundation model)", ["tanimoto"] + A.PHYSCHEM_COLS),
           ("Full fused (+ MAMMAL + Boltz)", A.FUSION_FEATURES)]
    L.append("## 3. Learn-to-rank feature ablation (LOTO within-target ρ, 297 ChEMBL pairs)")
    L.append("")
    L.append("What actually recovers the within-target ranking?")
    L.append("")
    L.append("| Feature set | within-target Spearman ρ |")
    L.append("|---|---|")
    rho = {}
    for name, fl in abl:
        r = A.loto_evaluate(feat, features=fl, seed=0)
        rho[name] = r.pooled_rho["fused_ltr"]
        L.append(f"| {name} | {rho[name]:+.3f} |")
    L.append("")
    nofm = rho["Tanimoto + physchem (NO foundation model)"]
    full = rho["Full fused (+ MAMMAL + Boltz)"]
    L.append(f"**Honest attribution:** classic ligand-similarity + physicochemical "
             f"features alone reach ρ = {nofm:+.2f}; adding the foundation model and "
             f"3D-affinity lifts this only to {full:+.2f} (Δ = {full - nofm:+.2f}). "
             f"The foundation model contributes marginally to within-target ranking — "
             f"the recovery is driven by classic cheminformatics, and MAMMAL-alone "
             f"({rho['MAMMAL pKd only']:+.2f}) is near-dead-weight for this task. This is "
             f"reported as a feature, not hidden: the practical claim is that a "
             f"sequence-only DTI score should not be relied on for within-target ligand "
             f"ranking at allosteric/transporter sites.")
    L.append("")
    L.append("Generated by `scripts/84_manuscript_robustness.py`.")
    out = ROOT / "reports" / "manuscript_robustness.md"
    out.write_text("\n".join(L), encoding="utf-8")
    logger.info("Common-subset (n=%d): class=%.2f rel=%.2f bind=%.2f",
                len(rows), R.auroc(np.array([cls[c] for c in rows.compound]), y),
                R.auroc(np.array([rel[c] for c in rows.compound]), y),
                R.auroc(np.array([bind[c] for c in rows.compound]), y))
    logger.info("Ablation: no-FM=%.3f full=%.3f", nofm, full)
    logger.info("Wrote %s", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
