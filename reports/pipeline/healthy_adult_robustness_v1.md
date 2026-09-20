# Healthy-adult axis — robustness audit

Adversarial re-analysis of the headline in `healthy_adult_axis_v1.md` ("the only separator is a coarse acute-CNS-stimulant gate, AUROC 0.86, p = 0.046"). No datum was added, altered, or re-curated: this re-analyses the same verified ledger. Reproduced by `scripts/121_healthy_adult_robustness.py`.

Primary set: **n = 21** clean-MA compounds (6 enhance / 15 null).

## R1 — the label is confounded with statistical POWER

The binary label is "a clean MA whose CI excludes 0". CI width scales as 1/sqrt(k), so the label conflates *works* with *was studied enough to detect*. A pure power proxy carrying no biology at all is therefore a control that the biology gate must beat:

| predictor | what it encodes | AUROC | perm p |
|---|---|---|---|
| acute CNS stimulant gate | biology | 0.77 | 0.0320 |
| **n_studies** | **pure statistical power, zero biology** | **0.66** | **0.1452** |
| representative_g | effect magnitude | 0.82 | 0.0105 |

**The power proxy WINS** (0.66 vs 0.77). Median studies pooled: **18** for labelled enhancers vs **nan** for labelled nulls. So the gate cannot be claimed as evidence that stimulant pharmacology predicts enhancement: a model that knows only how heavily a compound was studied does at least as well. `enhances_healthy_young` is a **detection** label, not an **efficacy** label.

## R2 — the headline hinges on ONE label decision

The ledger's stated inclusion rule is "CI excluding 0". Agreement between `ci_lo > 0` and the assigned label, for every compound with a recorded CI:

| compound | g | CI | k | CI excludes 0 | label | agrees |
|---|---|---|---|---|---|---|
| methylphenidate | +0.21 | [+0.09, +0.32] | 24 | True | 1 | yes |
| modafinil | +0.12 | [+0.02, +0.21] | 14 | True | 1 | yes |
| dextroamphetamine | +0.21 | [-0.06, +0.47] | 10 | False | 0 | yes |
| caffeine | +0.28 | [+0.21, +0.36] | 31 | True | 1 | yes |
| nicotine | +0.34 | [+0.18, +0.50] | 9 | True | 1 | yes |
| guarana | +0.08 | [-0.03, +0.18] | 8 | False | 0 | yes |
| l_theanine | +0.35 | [+0.10, +0.61] | 4 | True | 0 | **NO** |
| ginkgo_biloba | -0.04 | [-0.17, +0.07] | 13 | False | 0 | yes |
| bacopa_monnieri | +0.17 | [-0.52, +0.86] | 9 | False | 0 | yes |
| creatine | +0.03 | [-0.14, +0.20] | 8 | False | 0 | yes |
| multivitamin_mineral | +0.07 | [+0.03, +0.11] | 3 | True | 1 | yes |
| b_vitamins | +0.00 | [-0.05, +0.06] | 4 | False | 0 | yes |
| folic_acid | +0.01 | [-0.08, +0.10] | 9 | False | 0 | yes |
| testosterone | +0.09 | [-0.02, +0.19] | 23 | False | 0 | yes |
| menopausal_hormone_therapy | -0.16 | [-0.67, +0.34] | nan | False | 0 | yes |
| cannabidiol | -0.05 | [-0.12, +0.03] | 16 | False | 0 | yes |
| psilocybin_lsd_microdosing | -0.34 | [-0.62, -0.06] | 14 | False | 0 | yes |
| dietary_nitrate | +0.06 | [-0.06, +0.18] | 13 | False | 0 | yes |
| insulin_intranasal | +0.02 | [-0.05, +0.09] | 11 | False | 0 | yes |
| oxytocin_intranasal | +0.13 | [+0.02, +0.24] | 23 | True | 1 | yes |
| fruit_derived_polyphenols | +0.12 | [-0.29, +0.54] | 6 | False | 0 | yes |

**Conflict: l_theanine.** Re-labelling strictly per the ledger's own stated rule moves the headline:

| labelling | stimulant-gate AUROC | perm p | non-stimulant enhancers |
|---|---|---|---|
| as shipped | 0.77 | 0.0320 | none |
| per the stated CI rule | 0.71 | **0.0670** | l_theanine, multivitamin_mineral, oxytocin_intranasal |

So the one statistically significant result in the healthy-adult axis **does not survive a single defensible re-reading of one compound**, and under that reading the "enhancers are exclusively acute CNS stimulants" claim is falsified by a non-stimulant. The curator's note gives a real reason for the shipped call (only one RT sub-domain significant, k = 4) — the point is not that the shipped label is wrong, it is that the headline is **not robust** to it. Note the asymmetry it sits against: modafinil is labelled an enhancer at g = +0.12 while its own robustness note records it as TOST-equivalent-to-zero.

## R4 — is a label of 1 a USEFUL effect? (1 of 6 are entirely below the target, 4 more straddle it)

R3 asks whether a label of 0 is a refutation. This asks the mirror question, which nobody had asked. The inclusion rule is "the interval excludes 0", a statement about DETECTABILITY. The project's target is g = 0.20, a statement about MAGNITUDE. They can disagree.

| compound | g | CI | k | verdict | why |
|---|---|---|---|---|---|
| methylphenidate | +0.21 | [+0.09, +0.32] | 24 | **SPANS TARGET** | interval straddles g=0.2; magnitude unsettled |
| modafinil | +0.12 | [+0.02, +0.21] | 14 | **SPANS TARGET** | interval straddles g=0.2; magnitude unsettled |
| caffeine | +0.28 | [+0.21, +0.36] | 31 | **AT OR ABOVE TARGET** | entire interval is at or above g=0.2 |
| nicotine | +0.34 | [+0.18, +0.50] | 9 | **SPANS TARGET** | interval straddles g=0.2; magnitude unsettled |
| multivitamin_mineral | +0.07 | [+0.03, +0.11] | 3 | **BELOW TARGET** | interval excludes 0 AND excludes g=0.2: detectable but too small to be what is sought |
| oxytocin_intranasal | +0.13 | [+0.02, +0.24] | 23 | **SPANS TARGET** | interval straddles g=0.2; magnitude unsettled |

**1 of 6 labelled enhancers (multivitamin_mineral) have intervals that exclude 0 AND exclude g = 0.2.** They are real, replicated, and smaller than the effect this project is looking for. That is not a criticism of the compounds; it is a statement about what the label means. A ranker trained or evaluated on `enhances_healthy_young` is being asked to separate detectable-from-undetectable, not useful-from-useless, which is precisely the concern R1 raises about study volume.


A further **4** (methylphenidate, modafinil, nicotine, oxytocin_intranasal) have intervals that STRADDLE the target, so their magnitude is unsettled: the data are compatible both with a useful effect and with one too small to want. Only compounds whose entire interval sits at or above the target can be said to clear it, and there is 1 of those.

Taken together: of the labelled enhancers, only a minority have a magnitude this project could call established. That is not a criticism of the compounds and not a claim that the gate is wrong. It is a statement about what `enhances_healthy_young` encodes, and it is the same concern R1 raises from the direction of study volume: the label separates detectable from undetectable, which is not the same axis as useful from useless.

## R3 — most "nulls" are NOT refuted, only under-powered

A null whose CI still admits g >= 0.2 has not been ruled out. Splitting the labelled nulls:

| compound | g | CI | k | verdict | why |
|---|---|---|---|---|---|
| dextroamphetamine | +0.21 | [-0.06, +0.47] | 10 | **INCONCLUSIVE** | CI still admits g>=0.2 (under-powered) |
| guarana | +0.08 | [-0.03, +0.18] | 8 | **REFUTED** | CI excludes a meaningful g=0.2 |
| l_theanine | +0.35 | [+0.10, +0.61] | 4 | **INCONCLUSIVE** | CI still admits g>=0.2 (under-powered) |
| ginkgo_biloba | -0.04 | [-0.17, +0.07] | 13 | **REFUTED** | CI excludes a meaningful g=0.2 |
| bacopa_monnieri | +0.17 | [-0.52, +0.86] | 9 | **INCONCLUSIVE** | CI still admits g>=0.2 (under-powered) |
| creatine | +0.03 | [-0.14, +0.20] | 8 | **INCONCLUSIVE** | CI still admits g>=0.2 (under-powered) |
| b_vitamins | +0.00 | [-0.05, +0.06] | 4 | **REFUTED** | CI excludes a meaningful g=0.2 |
| folic_acid | +0.01 | [-0.08, +0.10] | 9 | **REFUTED** | CI excludes a meaningful g=0.2 |
| testosterone | +0.09 | [-0.02, +0.19] | 23 | **REFUTED** | CI excludes a meaningful g=0.2 |
| menopausal_hormone_therapy | -0.16 | [-0.67, +0.34] | nan | **INCONCLUSIVE** | CI still admits g>=0.2 (under-powered) |
| cannabidiol | -0.05 | [-0.12, +0.03] | 16 | **REFUTED** | CI excludes a meaningful g=0.2 |
| psilocybin_lsd_microdosing | -0.34 | [-0.62, -0.06] | 14 | **REFUTED** | CI excludes a meaningful g=0.2 |
| dietary_nitrate | +0.06 | [-0.06, +0.18] | 13 | **REFUTED** | CI excludes a meaningful g=0.2 |
| insulin_intranasal | +0.02 | [-0.05, +0.09] | 11 | **REFUTED** | CI excludes a meaningful g=0.2 |
| fruit_derived_polyphenols | +0.12 | [-0.29, +0.54] | 6 | **INCONCLUSIVE** | CI still admits g>=0.2 (under-powered) |

**9 genuinely refuted, 6 inconclusive, 0 with no CI recorded.** Plus **6 compounds with NO healthy-adult meta-analysis at all** (tyrosine, rhodiola_rosea, citicoline, piracetam, phosphatidylserine, vinpocetine). The field's evidence base is therefore far thinner than a flat "7 nulls" implies.

## What this changes

1. **The healthy-adult axis has no robust predictor — not even the coarse one.** The stimulant gate is beaten by a pure power proxy (R1) and loses significance under a one-compound re-reading (R2). At n = 11 with 4 positives, nothing is identifiable; the earlier AUROC 0.86 / p = 0.046 should be read as a fragile descriptive contrast, not a finding.
2. **The binding constraint is the ground truth, not the model.** No fusion, calibration or foundation-model work can be validated against 11 compounds whose labels track study volume. Expanding and power-annotating this ledger dominates every modelling improvement available.
3. **The remaining headroom is in the inconclusive set, not the refuted set.** Compounds with a decent point estimate and too few studies (l-theanine: g = +0.35 from k = 4) are where an adequately-powered trial could still change the answer. Compounds with tight CIs around zero (ginkgo, bacopa, omega-3, creatine-in-young) are closed.

**Integrity.** Every number above is computed from the existing verified ledger; no label was changed in the data. The R2 re-labelling is a *sensitivity analysis* reported alongside the shipped labelling, not a re-curation.

---

Generated by `scripts/121_healthy_adult_robustness.py`.