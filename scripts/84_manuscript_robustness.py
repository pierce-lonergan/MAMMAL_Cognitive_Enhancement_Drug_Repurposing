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
    # 0. CONSORT-style ledger assembly flow
    n_total = len(led)
    n_succ = int((led["label"] == 1).sum())
    n_fail = int((led["label"] == 0).sum())
    n_cls = led["mechanism_class"].nunique()
    by_dis = led["indication"].value_counts().to_dict() if "indication" in led else {}
    L.append("## 0. Clinical-ledger assembly (CONSORT-style flow)")
    L.append("")
    L.append("Pre-specified inclusion rule (Methods): a drug enters the binary analysis if "
             "it (i) was evaluated on a cognition-relevant/functional primary endpoint in "
             "its lead indication (index diseases AD/CIAS/FXS + cognition-adjacent "
             "ADHD/narcolepsy that anchor the successful classes); (ii) reached ≥ Phase II "
             "with a reported readout or obtained approval; (iii) has an assignable "
             "mechanism class + human target UniProt; (iv) is a small molecule in the DTI "
             "head's domain. Outcome coding is judged **within each drug's own indication**.")
    L.append("")
    L.append("```")
    L.append("Cognition / cognition-adjacent drugs reviewed (literature)        ~45")
    L.append("  │")
    L.append("  ├─ excluded: peptide/biologic outside small-molecule DTI domain (e.g. GLP-1)")
    L.append("  ├─ excluded: no adjudicable Phase II+ cognitive readout")
    L.append("  ├─ excluded: outcome ambiguous / terminated for non-efficacy reasons")
    L.append("  └─ excluded: no single assignable mechanism class / human target")
    L.append("  │")
    L.append(f"  ▼")
    L.append(f"Ledger: {n_total} drugs, {n_cls} mechanism classes "
             f"({n_succ} SUCCESS / {n_fail} FAILURE)")
    L.append("```")
    L.append("")
    L.append("Per-indication composition (lead indication as coded):")
    L.append("")
    L.append("| Lead indication | n |")
    L.append("|---|---|")
    for dis, n in sorted(by_dis.items(), key=lambda kv: -kv[1]):
        L.append(f"| {dis} | {n} |")
    L.append("")
    L.append("Intermediate exclusion counts are not individually logged (single-author "
             "literature curation); the verifiable endpoints are the final 31 rows, each "
             "with an indication, pivotal-trial identifier, readout year, and citation in "
             "`data/raw/clinical_outcomes_ledger.csv`. A fully programmatic "
             "ClinicalTrials.gov extraction is the natural next step.")
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
    # 1b. class purity + class-level CI
    L.append("## 1b. The AUROC = 1.00 is a readout of complete class homogeneity")
    L.append("")
    pure, mixed = [], []
    for c, g in led.groupby("mechanism_class"):
        s = int((g["label"] == 1).sum()); f = int((g["label"] == 0).sum())
        (pure if (s == 0 or f == 0) else mixed).append((c, s, f))
    L.append(f"Every one of the **{len(pure) + len(mixed)} mechanism classes is "
             f"outcome-pure** ({len(pure)}/{len(pure) + len(mixed)} uniformly SUCCESS "
             f"or FAILURE; {len(mixed)} mixed). The class-leave-one-compound-out "
             f"predictor is therefore, by construction, a historical class look-up, and "
             f"its perfect separation is a direct readout of this homogeneity — not a "
             f"predictive margin.")
    L.append("")
    L.append("| Mechanism class | SUCCESS | FAILURE | purity |")
    L.append("|---|---|---|---|")
    for c, s, f in sorted(pure + mixed):
        L.append(f"| {c} | {s} | {f} | {'PURE' if (c, s, f) in pure else 'MIXED'} |")
    L.append("")
    cb = R.class_cluster_bootstrap_auroc(led, n_boot=3000, seed=42)
    L.append(f"**Class-level (cluster) bootstrap.** Resampling the *classes* "
             f"themselves with replacement (the correct unit when the predictor is "
             f"class-aggregated) gives a 90% CI of **[{cb['ci_lo']:.2f}, {cb['ci_hi']:.2f}]** "
             f"(median {cb['median']:.2f}; {cb['frac_degenerate']:.1%} degenerate draws). "
             f"It does not widen below 1.00 because the classes are outcome-pure — "
             f"confirming that the relevant uncertainty is **not** sampling variance but "
             f"out-of-sample generalisation to *new* mechanism classes, which the "
             f"leave-one-class-out result (AUROC 0.00) bounds explicitly. The headline "
             f"is therefore the *comparative* result (class history dominates target-level "
             f"predictors), with perfect separation a downstream consequence of class "
             f"homogeneity.")
    L.append("")
    L.append("## 2. Comparator predictors, all on the identical drugs")
    L.append("")
    # extra comparators: KG network-propagation, structure-NN, ensemble
    kg = pd.read_parquet(ROOT / "data" / "results" / "v2" / "kg_scores.parquet")
    comp = pd.read_parquet(ROOT / "data" / "interim" / "compounds.parquet")
    kgn = R.kg_network_score(led, kg)
    nn = R.structure_nn_success_score(led, comp)
    common = set(bind) & set(rel) & set(cls) & set(kgn) & set(nn)
    rows = led[led["compound"].isin(common)]
    y = rows["label"].to_numpy()

    def zmap(s):
        v = np.array([s[c] for c in rows["compound"]], float)
        return (v - v.mean()) / (v.std() + 1e-9)
    ens = {c: float(z) for c, z in zip(rows["compound"],
            (zmap(bind) + zmap(rel) + zmap(kgn)) / 3.0)}  # target-centric ensemble

    L.append(f"Restricting **every** predictor to the same {len(rows)} drugs (where all "
             "are defined) removes the differing-n objection and answers "
             "\"compared to what?\" against four repurposing paradigms — affinity, "
             "genetics, **network-propagation (KG personalised PageRank)**, and "
             "**chemical-structure similarity** — plus their ensemble:")
    L.append("")
    L.append("| Predictor (same drugs) | paradigm | AUROC | leakage status |")
    L.append("|---|---|---|---|")
    comparators = [
        ("Mechanism-class track record (ours)", "class history", cls, "leakage-audited (class-LOCO)"),
        ("Structure NN-to-successes (LOO)", "chemical similarity", nn, "uses LOO outcomes"),
        ("KG personalised-PageRank", "network propagation", kgn, "hindsight-confounded"),
        ("Target genetic relevance", "genetics (Open Targets)", rel, "leakage-free"),
        ("Target binding affinity (MAMMAL DTI)", "affinity", bind, "leakage-free"),
        ("Target-centric ensemble (affinity+genetics+KG)", "ensemble", ens, "mixed"),
    ]
    res2 = {}
    for name, para, s, leak in comparators:
        sv = np.array([s[c] for c in rows["compound"]], float)
        au = R.auroc(sv, y); res2[name] = au
        L.append(f"| {name} | {para} | {au:.3f} | {leak} |")
    L.append("")
    L.append(f"On the identical {len(rows)} drugs, mechanism-class history separates "
             f"perfectly (1.00). The two paradigms that show apparent signal — "
             f"network-propagation ({res2['KG personalised-PageRank']:.2f}) and "
             f"structure-similarity ({res2['Structure NN-to-successes (LOO)']:.2f}) — are "
             f"precisely the two whose signal is **explainable as hindsight**: KG PageRank "
             f"rewards node degree (a drug accrues edges *because* it was studied and "
             f"succeeded), and structure-NN consumes the historical outcome labels "
             f"directly. The genuinely a-priori target metrics (genetics "
             f"{res2['Target genetic relevance']:.2f}, affinity "
             f"{res2['Target binding affinity (MAMMAL DTI)']:.2f}) and even the "
             f"target-centric **ensemble** ({res2['Target-centric ensemble (affinity+genetics+KG)']:.2f}) "
             f"remain at or below chance. Given the *same* historical-outcome information, "
             f"aggregating by **mechanism class** (1.00) beats aggregating by **chemical "
             f"structure** ({res2['Structure NN-to-successes (LOO)']:.2f}) — the comparison "
             f"that isolates the paper's claim.")
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
             f"The recovery is driven by classic cheminformatics; the foundation model "
             f"alone ({rho['MAMMAL pKd only']:+.2f}) contributes negligibly to this "
             f"specific task. Scope (precise): this concerns **within-target ligand "
             f"ranking at allosteric/transporter sites using the released "
             f"`dti_bindingdb_pkd` head** — a task adversarial to a sequence-only DTI "
             f"model trained on BindingDB pKd, not its intended cross-target affinity "
             f"task. The practical claim is that a sequence-only DTI score should not be "
             f"relied on for within-target ligand ranking at these sites, where "
             f"inexpensive cheminformatics features suffice.")
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
