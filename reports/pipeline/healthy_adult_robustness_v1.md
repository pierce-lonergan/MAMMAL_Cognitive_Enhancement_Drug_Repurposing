# Healthy-adult axis — robustness audit

Adversarial re-analysis of the headline in `healthy_adult_axis_v1.md` ("the only separator is a coarse acute-CNS-stimulant gate, AUROC 0.86, p = 0.046"). No datum was added, altered, or re-curated: this re-analyses the same verified ledger. Reproduced by `scripts/121_healthy_adult_robustness.py`.

Primary set: **n = 11** clean-MA compounds (4 enhance / 7 null).

## R1 — the label is confounded with statistical POWER

The binary label is "a clean MA whose CI excludes 0". CI width scales as 1/sqrt(k), so the label conflates *works* with *was studied enough to detect*. A pure power proxy carrying no biology at all is therefore a control that the biology gate must beat:

| predictor | what it encodes | AUROC | perm p |
|---|---|---|---|
| acute CNS stimulant gate | biology | 0.86 | 0.0456 |
| **n_studies** | **pure statistical power, zero biology** | **0.88** | **0.0223** |
| representative_g | effect magnitude | 0.80 | 0.0676 |

**The power proxy WINS** (0.88 vs 0.86). Median studies pooled: **19** for labelled enhancers vs **9** for labelled nulls. So the gate cannot be claimed as evidence that stimulant pharmacology predicts enhancement: a model that knows only how heavily a compound was studied does at least as well. `enhances_healthy_young` is a **detection** label, not an **efficacy** label.

## R2 — the headline hinges on ONE label decision

The ledger's stated inclusion rule is "CI excluding 0". Agreement between `ci_lo > 0` and the assigned label, for every compound with a recorded CI:

| compound | g | CI | k | CI excludes 0 | label | agrees |
|---|---|---|---|---|---|---|
| methylphenidate | +0.21 | [+0.09, +0.32] | 24 | True | 1 | yes |
| modafinil | +0.12 | [+0.02, +0.21] | 14 | True | 1 | yes |
| dextroamphetamine | +0.21 | [-0.06, +0.47] | 10 | False | 0 | yes |
| nicotine | +0.34 | [+0.18, +0.50] | 9 | True | 1 | yes |
| l_theanine | +0.35 | [+0.10, +0.61] | 4 | True | 0 | **NO** |
| creatine | +0.03 | [-0.14, +0.20] | 8 | False | 0 | yes |

**Conflict: l_theanine.** Re-labelling strictly per the ledger's own stated rule moves the headline:

| labelling | stimulant-gate AUROC | perm p | non-stimulant enhancers |
|---|---|---|---|
| as shipped | 0.86 | 0.0456 | none |
| per the stated CI rule | 0.73 | **0.1758** | l_theanine |

So the one statistically significant result in the healthy-adult axis **does not survive a single defensible re-reading of one compound**, and under that reading the "enhancers are exclusively acute CNS stimulants" claim is falsified by a non-stimulant. The curator's note gives a real reason for the shipped call (only one RT sub-domain significant, k = 4) — the point is not that the shipped label is wrong, it is that the headline is **not robust** to it. Note the asymmetry it sits against: modafinil is labelled an enhancer at g = +0.12 while its own robustness note records it as TOST-equivalent-to-zero.

## R3 — most "nulls" are NOT refuted, only under-powered

A null whose CI still admits g >= 0.2 has not been ruled out. Splitting the labelled nulls:

| compound | g | CI | k | verdict | why |
|---|---|---|---|---|---|
| dextroamphetamine | +0.21 | [-0.06, +0.47] | 10 | **INCONCLUSIVE** | CI still admits g>=0.2 (under-powered) |
| guarana | +0.08 | not recorded | 8 | **NO CI RECORDED** | cannot distinguish refuted from under-powered |
| l_theanine | +0.35 | [+0.10, +0.61] | 4 | **INCONCLUSIVE** | CI still admits g>=0.2 (under-powered) |
| ginkgo_biloba | -0.04 | not recorded | 13 | **NO CI RECORDED** | cannot distinguish refuted from under-powered |
| bacopa_monnieri | +0.00 | not recorded | 9 | **NO CI RECORDED** | cannot distinguish refuted from under-powered |
| omega_3 | +0.00 | not recorded | 11 | **NO CI RECORDED** | cannot distinguish refuted from under-powered |
| creatine | +0.03 | [-0.14, +0.20] | 8 | **INCONCLUSIVE** | CI still admits g>=0.2 (under-powered) |

**0 genuinely refuted, 3 inconclusive, 4 with no CI recorded.** Plus **6 compounds with NO healthy-adult meta-analysis at all** (tyrosine, rhodiola_rosea, citicoline, piracetam, phosphatidylserine, vinpocetine). The field's evidence base is therefore far thinner than a flat "7 nulls" implies.

**NOT ONE labelled null has actually been refuted** at the g >= 0.2 threshold. Every one is either under-powered or has no CI on record. The honest reading of the healthy-adult evidence base is therefore *absence of evidence*, not *evidence of absence*: the pipeline cannot presently support a claim of the form "compound X is ruled out as a meaningful enhancer" for a single compound in the set.

## What this changes

1. **The healthy-adult axis has no robust predictor — not even the coarse one.** The stimulant gate is beaten by a pure power proxy (R1) and loses significance under a one-compound re-reading (R2). At n = 11 with 4 positives, nothing is identifiable; the earlier AUROC 0.86 / p = 0.046 should be read as a fragile descriptive contrast, not a finding.
2. **The binding constraint is the ground truth, not the model.** No fusion, calibration or foundation-model work can be validated against 11 compounds whose labels track study volume. Expanding and power-annotating this ledger dominates every modelling improvement available.
3. **The remaining headroom is in the inconclusive set, not the refuted set.** Compounds with a decent point estimate and too few studies (l-theanine: g = +0.35 from k = 4) are where an adequately-powered trial could still change the answer. Compounds with tight CIs around zero (ginkgo, bacopa, omega-3, creatine-in-young) are closed.

**Integrity.** Every number above is computed from the existing verified ledger; no label was changed in the data. The R2 re-labelling is a *sensitivity analysis* reported alongside the shipped labelling, not a re-curation.

---

Generated by `scripts/121_healthy_adult_robustness.py`.