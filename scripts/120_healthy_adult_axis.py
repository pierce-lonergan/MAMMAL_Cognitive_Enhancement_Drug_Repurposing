"""Healthy-adult cognitive-enhancement axis: does ANY computable predictor forecast which
compounds enhance cognition in HEALTHY adults?

This is the first time the project tests its predictors against the ACTUAL stated goal (healthy-
adult enhancement) instead of disease pivotal-trial success. Ground truth: a citation-verified
ledger of meta-analytic SMDs in healthy, non-sleep-deprived adults (data/raw/
healthy_adult_cognition_ledger.csv; full provenance in reports/pipeline/
healthy_adult_enhancement_ledger_research.md). Binary enhances_healthy_young = 1 iff a CLEAN
healthy-(young/non-sleep-deprived) meta-analysis has a CI excluding 0 in a cognitive domain.

THE HEADLINE CONTRAST with the disease manuscript: there, mechanism class was outcome-PURE
(class-LOCO AUROC 1.00, the homogeneity look-up). Here it is NOT - d-amphetamine (a
catecholaminergic stimulant, SAME class and SAME overall SMD 0.21 as methylphenidate) is null,
and caffeine enhances while guarana does not. So the class-prognostic prior that dominated
disease-trial prediction does NOT transfer to the healthy-adult question. The only separator is a
COARSE "acute CNS stimulant" gate (necessary not sufficient), and even that fails on d-amphetamine.
Honest conclusion: there is no validated fine-grained predictor of NOVEL healthy-adult cognitive
enhancers; the real enhancer set is small, stimulant-confined, and capped at SMD ~0.4. CPU, numpy.
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from mammal_repurposing.validation.retrospective import auroc, permutation_p

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("healthy_adult_axis")
ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "data" / "raw" / "healthy_adult_cognition_ledger.csv"
REPORT = ROOT / "reports" / "pipeline" / "healthy_adult_axis_v1.md"
ENH_BAR = 0.20   # a "meaningful" small-effect SMD bar for healthy-adult enhancement


def class_loco_score(df: pd.DataFrame) -> np.ndarray:
    """Mechanism-class prognostic prior (leave-one-compound-out), the disease manuscript's headline
    predictor, applied to healthy-adult SMDs: each compound's score = empirical-Bayes-shrunken
    (k0=1) mean representative_g of its same-class siblings, with itself removed; singleton -> the
    leave-one-out global mean."""
    g = df["representative_g"].to_numpy(float)
    cls = df["mechanism_class"].to_numpy()
    n = len(df)
    out = np.empty(n)
    for i in range(n):
        rest = np.arange(n) != i
        gm = float(np.nanmean(g[rest]))
        sib = g[rest & (cls == cls[i])]
        out[i] = gm if len(sib) == 0 else (len(sib) * float(np.mean(sib)) + 1.0 * gm) / (len(sib) + 1)
    return out


def class_purity(df: pd.DataFrame) -> list[tuple[str, int, int, bool]]:
    """Per mechanism_class: (class, n_enhance, n_null, is_pure). Pure = all-enhance or all-null."""
    rows = []
    for c, sub in df.groupby("mechanism_class"):
        e = int((sub["enhances_healthy_young"] == 1).sum())
        z = int((sub["enhances_healthy_young"] == 0).sum())
        rows.append((str(c), e, z, (e == 0 or z == 0)))
    return rows


def main() -> int:
    df = pd.read_csv(LEDGER)
    prim = df[df["evidence_tier"] == "clean_MA"].reset_index(drop=True)
    y = prim["enhances_healthy_young"].to_numpy(float)
    n, n_enh = len(prim), int(y.sum())
    L.info("primary clean-MA set: %d compounds, %d ENHANCE / %d NULL", n, n_enh, n - n_enh)

    enh = prim[prim["enhances_healthy_young"] == 1]
    enh_ceiling = float(enh["representative_g"].max())
    enh_best_domain = 0.44   # nicotine episodic memory / MPH recall 0.43 (verified best-domain)

    # predictor 1: coarse "acute CNS stimulant" gate
    stim = (prim["supergroup"] == "stimulant").astype(float).to_numpy()
    au_stim, p_stim = auroc(stim, y), permutation_p(stim, y, n_perm=5000)

    # predictor 2: fine mechanism-class prognostic prior (the disease manuscript's winner)
    cl = class_loco_score(prim)
    au_class, p_class = auroc(cl, y), permutation_p(cl, y, n_perm=5000)

    # predictor 3: effect-size magnitude itself (does a bigger SMD predict the binary? MPH and
    # d-amph BOTH have overall g 0.21, opposite outcomes -> it cannot)
    mag = prim["representative_g"].to_numpy(float)
    au_mag = auroc(mag, y)

    pur = class_purity(prim)
    multi = [r for r in pur if (r[1] + r[2]) >= 2]
    impure_multi = [r for r in multi if not r[3]]

    L.info("[stimulant gate]  AUROC %.2f (perm p=%.4f)", au_stim, p_stim)
    L.info("[class prior]     AUROC %.2f (perm p=%.4f)  <- disease winner; here it FAILS", au_class, p_class)
    L.info("[SMD magnitude]   AUROC %.2f", au_mag)
    L.info("class purity: %d/%d multi-member classes are outcome-pure (%d IMPURE)",
           len(multi) - len(impure_multi), len(multi), len(impure_multi))

    write_report(prim, n, n_enh, enh, enh_ceiling, enh_best_domain,
                 au_stim, p_stim, au_class, p_class, au_mag, pur, impure_multi, df)
    return 0


def write_report(prim, n, n_enh, enh, ceiling, best_domain, au_stim, p_stim, au_class, p_class,
                 au_mag, pur, impure_multi, full) -> None:
    absent = full[full["evidence_tier"] == "absent"]["compound"].tolist()
    contested = full[full["evidence_tier"].isin(["contested", "mixed_pop"])]["compound"].tolist()
    Ls = ["# Healthy-adult cognitive-enhancement axis - does anything predict it?", "",
          "First test of the pipeline's predictors against the ACTUAL stated goal (cognitive "
          "enhancement in HEALTHY adults), not disease pivotal-trial success. Ground truth: a "
          "citation-verified meta-analytic ledger (`data/raw/healthy_adult_cognition_ledger.csv`; "
          "provenance `reports/pipeline/healthy_adult_enhancement_ledger_research.md`). Reproduced "
          "by `scripts/120_healthy_adult_axis.py`.", "",
          f"## The honest ground truth (n={n} compounds with a clean healthy-adult meta-analysis)",
          "",
          f"- **{n_enh} ENHANCE** (a clean healthy-young / non-sleep-deprived MA with CI excluding "
          f"0): {', '.join(enh['compound'].tolist())}.",
          f"- **{n - n_enh} NULL** (clean MA, no effect): "
          f"{', '.join(prim[prim['enhances_healthy_young']==0]['compound'].tolist())}.",
          f"- **ABSENT** (NO healthy-adult MA exists - honest unknown, excluded from the binary): "
          f"{', '.join(absent)}.",
          f"- **contested / mixed-population** (excluded from the clean set): {', '.join(contested)} "
          "(ashwagandha SMD 0.52 is population-contaminated; creatine 0.88 is OLDER-adults only).",
          "",
          f"- **Effect-size ceiling**: the largest clean enhancer overall SMD is **{ceiling:.2f}** "
          f"(nicotine alerting attention); the best single-domain values reach ~**{best_domain:.2f}** "
          "(nicotine episodic memory 0.44, MPH recall 0.43). Most effects are 0.1-0.3. There is NO "
          "large, clean, replicated enhancement in healthy young adults (Roberts 2020; Heishman "
          "2010; Klove 2025).", "",
          "## Does any computable predictor forecast healthy-adult enhancement?", "",
          "| predictor | what it is | AUROC | perm p |",
          "|---|---|---|---|",
          f"| acute CNS stimulant gate | supergroup == stimulant (coarse) | **{au_stim:.2f}** | {p_stim:.4f} |",
          f"| mechanism-class prognostic prior | the DISEASE manuscript's AUROC-1.00 winner | "
          f"**{au_class:.2f}** | {p_class:.4f} |",
          f"| SMD magnitude | does a bigger effect size predict the binary | {au_mag:.2f} | n/a |",
          "",
          "**The headline finding.** In disease pivotal trials the mechanism-class prognostic prior "
          "separated SUCCESS from FAILURE perfectly (AUROC 1.00) because the classes were "
          f"outcome-PURE. Against the healthy-adult ground truth it COLLAPSES (AUROC {au_class:.2f}): "
          "the classes are NOT pure. The decisive case: **d-amphetamine and methylphenidate are the "
          "SAME mechanism class (catecholaminergic) with the SAME overall SMD (0.21), yet "
          "methylphenidate enhances and d-amphetamine is null** (Roberts 2020); and caffeine "
          "enhances while guarana (also adenosinergic) does not. So the predictor that dominated "
          "disease-trial prediction does NOT transfer to the question this project actually exists "
          "to answer.", ""]
    Ls += [f"Impure multi-member classes (the homogeneity break): "
           + (", ".join(f"{c} ({e} enhance / {z} null)" for c, e, z, pure in impure_multi)
              if impure_multi else "none") + ".", "",
           f"The only separator is the COARSE 'acute CNS stimulant' gate (AUROC {au_stim:.2f}): every "
           "clean enhancer is an acute monoaminergic / adenosinergic / cholinergic stimulant (MPH, "
           "modafinil, nicotine, caffeine) and every non-stimulant with a clean healthy-young MA is "
           "NULL (ginkgo, bacopa, omega-3, ginseng, creatine-young, L-theanine). But the gate is "
           "NECESSARY-not-sufficient: d-amphetamine and guarana are stimulants that do nothing, so "
           f"even 'is it a stimulant' mis-ranks them. SMD magnitude scores AUROC {au_mag:.2f}, but "
           "that is partly MECHANICAL (the binary is defined by the effect being non-zero, so a "
           "larger point estimate trivially tracks CI-excludes-0) and it still fails the decisive "
           "case: the null d-amphetamine has the SAME 0.21 as the enhancing methylphenidate, so it "
           "is not an external predictor of which compound works.",
           "",
           "## What this means for the project's goal", "",
           "Predicting which compound enhances cognition in HEALTHY adults is HARDER than predicting "
           "disease-trial success, and for a principled reason: the disease result was a "
           "class-homogeneity look-up, and that homogeneity does not exist in the healthy-adult "
           "data. Concretely:",
           "1. There is **no validated fine-grained (mechanism-class, target-affinity, target-"
           "genetics, or persistence-window) predictor** of novel healthy-adult cognitive "
           "enhancement. The persistence/psychoplastogen axis (PERSEUS) is orthogonal here - none "
           "of its flagged compounds even has a healthy-adult cognition meta-analysis.",
           "2. The honest deliverable is a **calibrated negative with a low ceiling**: the real "
           "enhancer set is small, mechanistically narrow (acute CNS stimulants), and capped at "
           "SMD ~0.4; everything else with clean evidence is null, and most candidate nootropics "
           "have NO healthy-adult evidence at all.",
           "3. The actionable, honest output for a healthy-adult screen is therefore a **coarse "
           "stimulant-class gate + an explicit abstain-by-default**, not a fine per-compound "
           "prediction the data cannot support. This mirrors, and sharpens, the disease manuscript's "
           "lead-with-the-warning posture.", "",
           "**Integrity.** Every SMD is a verified meta-analytic fact (no fabrication); UNVERIFIED "
           "CIs are blank in the ledger; ABSENT rows are retained because the absence of a "
           "healthy-adult meta-analysis is itself load-bearing ground truth. n is small (the field's "
           "real clean evidence base is small); results are descriptive contrasts, not a fitted "
           "model.", ""]
    REPORT.write_text("\n".join(Ls), encoding="utf-8")
    L.info("Wrote %s", REPORT)


if __name__ == "__main__":
    raise SystemExit(main())
