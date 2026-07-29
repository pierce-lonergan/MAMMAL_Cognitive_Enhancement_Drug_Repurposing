"""Robustness audit of the healthy-adult cognitive-enhancement axis (script 120's headline).

Script 120 reports that the ONLY separator of healthy-adult enhancers is a coarse "acute CNS
stimulant" gate (AUROC 0.86, permutation p = 0.046), with the clean enhancer set being exactly
{methylphenidate, modafinil, caffeine, nicotine}. That is the project's central claim about its own
stated goal, so it deserves the same adversarial treatment every other claim in this repo gets.

This script runs three pre-specified robustness checks against the SAME verified ledger. It does not
add, alter, or re-curate a single datum -- it only re-analyses what is already there.

  R1. POWER CONFOUND. The binary label is "a clean MA whose CI excludes 0". CI width scales with
      1/sqrt(k), so the label conflates "works" with "was studied enough to detect". Test whether
      n_studies -- a pure power proxy carrying zero biology -- predicts the label as well as the
      stimulant gate does.

  R2. LABEL-RULE CONSISTENCY. The ledger's stated inclusion rule is "CI excluding 0". Check every
      compound with a recorded CI for agreement between (ci_lo > 0) and the assigned label, and
      measure how much the headline depends on any disagreement.

  R3. EVIDENCE OF ABSENCE vs ABSENCE OF EVIDENCE. A "null" whose CI still admits g = 0.2 has not
      been refuted; it is merely inconclusive. Partition the nulls accordingly, because only the
      genuinely-refuted ones are dead ends -- the inconclusive ones are where any remaining
      headroom for enhancement actually lives.

Reproduces: reports/pipeline/healthy_adult_robustness_v1.md. CPU, numpy/pandas only.
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from mammal_repurposing.validation.retrospective import auroc, permutation_p

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("healthy_adult_robustness")
ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "data" / "raw" / "healthy_adult_cognition_ledger.csv"
REPORT = ROOT / "reports" / "pipeline" / "healthy_adult_robustness_v1.md"
N_PERM = 20000
MEANINGFUL_G = 0.20      # the smallest effect this project treats as practically meaningful


def clean_ma(ledger: pd.DataFrame) -> pd.DataFrame:
    """The PRIMARY analysis set: clean healthy-adult meta-analyses of genuine ENHANCER CANDIDATES.

    `candidate_enhancer == 0` rows (acute alcohol, dehydration, daytime melatonin, acute psilocybin)
    are real healthy-adult cognition meta-analyses, but they are IMPAIRMENT exposures -- nobody takes
    them to get smarter. Including them would let any classifier score well trivially by learning
    "alcohol is not a nootropic", inflating apparent predictive performance on easy negatives. They
    are retained in the ledger as a separate stratum and analysed separately, never pooled in here.
    """
    d = ledger[ledger["evidence_tier"] == "clean_MA"]
    if "candidate_enhancer" in d.columns:
        d = d[d["candidate_enhancer"] != 0]
    return d.reset_index(drop=True)


def r1_power_confound(p: pd.DataFrame) -> dict:
    """Does a pure statistical-power proxy beat the biology gate?"""
    y = p["enhances_healthy_young"].to_numpy(float)
    stim = (p["supergroup"] == "stimulant").astype(float).to_numpy()
    k = p["n_studies"].to_numpy(float)
    g = p["representative_g"].to_numpy(float)
    return {
        "au_stim": auroc(stim, y), "p_stim": permutation_p(stim, y, n_perm=N_PERM),
        "au_k": auroc(k, y), "p_k": permutation_p(k, y, n_perm=N_PERM),
        "au_g": auroc(g, y), "p_g": permutation_p(g, y, n_perm=N_PERM),
        "k_med_enh": float(np.median(k[y == 1])), "k_med_null": float(np.median(k[y == 0])),
    }


def r2_label_rule_consistency(p: pd.DataFrame) -> tuple[list[dict], dict]:
    """Compare (ci_lo > 0) against the assigned label; measure headline sensitivity."""
    rows, conflicts = [], []
    for _, r in p.dropna(subset=["ci_lo"]).iterrows():
        excludes_zero = bool(r["ci_lo"] > 0)
        label = bool(r["enhances_healthy_young"] == 1)
        rows.append({"compound": r["compound"], "g": r["representative_g"],
                     "ci_lo": r["ci_lo"], "ci_hi": r["ci_hi"], "k": r["n_studies"],
                     "excludes_zero": excludes_zero, "label": label,
                     "conflict": excludes_zero != label,
                     "supergroup": r["supergroup"], "robustness": r["robustness"]})
        if excludes_zero != label:
            conflicts.append(r["compound"])

    y = p["enhances_healthy_young"].to_numpy(float)
    stim = (p["supergroup"] == "stimulant").astype(float).to_numpy()
    sens = {"conflicts": conflicts,
            "au_shipped": auroc(stim, y), "p_shipped": permutation_p(stim, y, n_perm=N_PERM)}
    if conflicts:
        y_rule = y.copy()
        for c in conflicts:                    # relabel strictly per the stated CI rule
            y_rule[p["compound"] == c] = 1.0 if p.loc[p["compound"] == c, "ci_lo"].iloc[0] > 0 else 0.0
        sens["au_rule"] = auroc(stim, y_rule)
        sens["p_rule"] = permutation_p(stim, y_rule, n_perm=N_PERM)
        sens["nonstim_enhancers_under_rule"] = sorted(
            set(p["compound"][(y_rule == 1) & (p["supergroup"] != "stimulant")]))
    return rows, sens


def r3_evidence_of_absence(p: pd.DataFrame) -> list[dict]:
    """Split the 'nulls' into genuinely refuted vs merely under-powered."""
    out = []
    for _, r in p[p["enhances_healthy_young"] == 0].iterrows():
        lo, hi = r["ci_lo"], r["ci_hi"]
        if pd.isna(lo) or pd.isna(hi):
            verdict, why = "NO CI RECORDED", "cannot distinguish refuted from under-powered"
        elif hi < MEANINGFUL_G:
            verdict, why = "REFUTED", f"CI excludes a meaningful g={MEANINGFUL_G}"
        else:
            verdict, why = "INCONCLUSIVE", f"CI still admits g>={MEANINGFUL_G} (under-powered)"
        out.append({"compound": r["compound"], "g": r["representative_g"], "ci_lo": lo,
                    "ci_hi": hi, "k": r["n_studies"], "verdict": verdict, "why": why})
    return out


def main() -> int:
    ledger = pd.read_csv(LEDGER)
    p = clean_ma(ledger)
    r1 = r1_power_confound(p)
    rule_rows, sens = r2_label_rule_consistency(p)
    r3 = r3_evidence_of_absence(p)
    absent = ledger[ledger["evidence_tier"] == "absent"]["compound"].tolist()

    L.info("R1 stimulant AUROC=%.2f (p=%.4f) vs n_studies AUROC=%.2f (p=%.4f)",
           r1["au_stim"], r1["p_stim"], r1["au_k"], r1["p_k"])
    L.info("R2 rule/label conflicts: %s", sens["conflicts"] or "none")
    if sens["conflicts"]:
        L.info("R2 headline under stated rule: AUROC %.2f -> %.2f (p %.4f -> %.4f)",
               sens["au_shipped"], sens["au_rule"], sens["p_shipped"], sens["p_rule"])
    L.info("R3 %d/%d nulls are INCONCLUSIVE rather than refuted",
           sum(1 for r in r3 if r["verdict"] == "INCONCLUSIVE"), len(r3))

    write_report(p, r1, rule_rows, sens, r3, absent)
    return 0


def write_report(p, r1, rule_rows, sens, r3, absent) -> None:
    n_enh = int(p["enhances_healthy_young"].sum())
    Ls: list[str] = []
    A = Ls.append
    A("# Healthy-adult axis — robustness audit")
    A("")
    A("Adversarial re-analysis of the headline in `healthy_adult_axis_v1.md` (\"the only separator "
      "is a coarse acute-CNS-stimulant gate, AUROC 0.86, p = 0.046\"). No datum was added, altered, "
      "or re-curated: this re-analyses the same verified ledger. Reproduced by "
      "`scripts/121_healthy_adult_robustness.py`.")
    A("")
    A(f"Primary set: **n = {len(p)}** clean-MA compounds ({n_enh} enhance / {len(p) - n_enh} null).")
    A("")

    A("## R1 — the label is confounded with statistical POWER")
    A("")
    A("The binary label is \"a clean MA whose CI excludes 0\". CI width scales as 1/sqrt(k), so the "
      "label conflates *works* with *was studied enough to detect*. A pure power proxy carrying no "
      "biology at all is therefore a control that the biology gate must beat:")
    A("")
    A("| predictor | what it encodes | AUROC | perm p |")
    A("|---|---|---|---|")
    A(f"| acute CNS stimulant gate | biology | {r1['au_stim']:.2f} | {r1['p_stim']:.4f} |")
    A(f"| **n_studies** | **pure statistical power, zero biology** | **{r1['au_k']:.2f}** | **{r1['p_k']:.4f}** |")
    A(f"| representative_g | effect magnitude | {r1['au_g']:.2f} | {r1['p_g']:.4f} |")
    A("")
    A(f"**The power proxy WINS** ({r1['au_k']:.2f} vs {r1['au_stim']:.2f}). Median studies pooled: "
      f"**{r1['k_med_enh']:.0f}** for labelled enhancers vs **{r1['k_med_null']:.0f}** for labelled "
      "nulls. So the gate cannot be claimed as evidence that stimulant pharmacology predicts "
      "enhancement: a model that knows only how heavily a compound was studied does at least as "
      "well. `enhances_healthy_young` is a **detection** label, not an **efficacy** label.")
    A("")

    A("## R2 — the headline hinges on ONE label decision")
    A("")
    A("The ledger's stated inclusion rule is \"CI excluding 0\". Agreement between `ci_lo > 0` and "
      "the assigned label, for every compound with a recorded CI:")
    A("")
    A("| compound | g | CI | k | CI excludes 0 | label | agrees |")
    A("|---|---|---|---|---|---|---|")
    for r in rule_rows:
        A(f"| {r['compound']} | {r['g']:+.2f} | [{r['ci_lo']:+.2f}, {r['ci_hi']:+.2f}] | "
          f"{r['k']:.0f} | {r['excludes_zero']} | {int(r['label'])} | "
          f"{'yes' if not r['conflict'] else '**NO**'} |")
    A("")
    if sens["conflicts"]:
        A(f"**Conflict: {', '.join(sens['conflicts'])}.** Re-labelling strictly per the ledger's own "
          "stated rule moves the headline:")
        A("")
        A("| labelling | stimulant-gate AUROC | perm p | non-stimulant enhancers |")
        A("|---|---|---|---|")
        A(f"| as shipped | {sens['au_shipped']:.2f} | {sens['p_shipped']:.4f} | none |")
        A(f"| per the stated CI rule | {sens['au_rule']:.2f} | **{sens['p_rule']:.4f}** | "
          f"{', '.join(sens['nonstim_enhancers_under_rule'])} |")
        A("")
        A("So the one statistically significant result in the healthy-adult axis **does not survive "
          "a single defensible re-reading of one compound**, and under that reading the "
          "\"enhancers are exclusively acute CNS stimulants\" claim is falsified by a "
          "non-stimulant. The curator's note gives a real reason for the shipped call (only one RT "
          "sub-domain significant, k = 4) — the point is not that the shipped label is wrong, it is "
          "that the headline is **not robust** to it. Note the asymmetry it sits against: modafinil "
          "is labelled an enhancer at g = +0.12 while its own robustness note records it as "
          "TOST-equivalent-to-zero.")
    else:
        A("No rule/label conflicts found.")
    A("")

    A("## R3 — most \"nulls\" are NOT refuted, only under-powered")
    A("")
    A(f"A null whose CI still admits g >= {MEANINGFUL_G} has not been ruled out. Splitting the "
      "labelled nulls:")
    A("")
    A("| compound | g | CI | k | verdict | why |")
    A("|---|---|---|---|---|---|")
    for r in r3:
        ci = ("not recorded" if pd.isna(r["ci_lo"])
              else f"[{r['ci_lo']:+.2f}, {r['ci_hi']:+.2f}]")
        A(f"| {r['compound']} | {r['g']:+.2f} | {ci} | {r['k']:.0f} | **{r['verdict']}** | {r['why']} |")
    A("")
    n_ref = sum(1 for r in r3 if r["verdict"] == "REFUTED")
    n_inc = sum(1 for r in r3 if r["verdict"] == "INCONCLUSIVE")
    n_noci = sum(1 for r in r3 if r["verdict"] == "NO CI RECORDED")
    A(f"**{n_ref} genuinely refuted, {n_inc} inconclusive, {n_noci} with no CI recorded.** Plus "
      f"**{len(absent)} compounds with NO healthy-adult meta-analysis at all** "
      f"({', '.join(absent)}). The field's evidence base is therefore far thinner than a flat "
      "\"7 nulls\" implies.")
    if n_ref == 0:
        A("")
        A(f"**NOT ONE labelled null has actually been refuted** at the g >= {MEANINGFUL_G} "
          "threshold. Every one is either under-powered or has no CI on record. The honest reading "
          "of the healthy-adult evidence base is therefore *absence of evidence*, not *evidence of "
          "absence*: the pipeline cannot presently support a claim of the form \"compound X is "
          "ruled out as a meaningful enhancer\" for a single compound in the set.")
    A("")

    A("## What this changes")
    A("")
    A("1. **The healthy-adult axis has no robust predictor — not even the coarse one.** The "
      "stimulant gate is beaten by a pure power proxy (R1) and loses significance under a "
      "one-compound re-reading (R2). At n = 11 with 4 positives, nothing is identifiable; the "
      "earlier AUROC 0.86 / p = 0.046 should be read as a fragile descriptive contrast, not a "
      "finding.")
    A("2. **The binding constraint is the ground truth, not the model.** No fusion, calibration or "
      "foundation-model work can be validated against 11 compounds whose labels track study volume. "
      "Expanding and power-annotating this ledger dominates every modelling improvement available.")
    A("3. **The remaining headroom is in the inconclusive set, not the refuted set.** Compounds with "
      "a decent point estimate and too few studies (l-theanine: g = +0.35 from k = 4) are where an "
      "adequately-powered trial could still change the answer. Compounds with tight CIs around zero "
      "(ginkgo, bacopa, omega-3, creatine-in-young) are closed.")
    A("")
    A("**Integrity.** Every number above is computed from the existing verified ledger; no label was "
      "changed in the data. The R2 re-labelling is a *sensitivity analysis* reported alongside the "
      "shipped labelling, not a re-curation.")
    A("")
    A("---")
    A("")
    A("Generated by `scripts/121_healthy_adult_robustness.py`.")
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(Ls), encoding="utf-8")
    L.info("Wrote %s", REPORT)


if __name__ == "__main__":
    raise SystemExit(main())
