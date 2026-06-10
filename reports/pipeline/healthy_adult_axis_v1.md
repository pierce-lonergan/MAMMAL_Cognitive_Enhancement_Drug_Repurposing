# Healthy-adult cognitive-enhancement axis - does anything predict it?

First test of the pipeline's predictors against the ACTUAL stated goal (cognitive enhancement in HEALTHY adults), not disease pivotal-trial success. Ground truth: a citation-verified meta-analytic ledger (`data/raw/healthy_adult_cognition_ledger.csv`; provenance `reports/pipeline/healthy_adult_enhancement_ledger_research.md`). Reproduced by `scripts/120_healthy_adult_axis.py`.

## The honest ground truth (n=11 compounds with a clean healthy-adult meta-analysis)

- **4 ENHANCE** (a clean healthy-young / non-sleep-deprived MA with CI excluding 0): methylphenidate, modafinil, caffeine, nicotine.
- **7 NULL** (clean MA, no effect): dextroamphetamine, guarana, l_theanine, ginkgo_biloba, bacopa_monnieri, omega_3, creatine.
- **ABSENT** (NO healthy-adult MA exists - honest unknown, excluded from the binary): tyrosine, rhodiola_rosea, citicoline, piracetam, phosphatidylserine, vinpocetine.
- **contested / mixed-population** (excluded from the clean set): panax_ginseng, ashwagandha (ashwagandha SMD 0.52 is population-contaminated; creatine 0.88 is OLDER-adults only).

- **Effect-size ceiling**: the largest clean enhancer overall SMD is **0.34** (nicotine alerting attention); the best single-domain values reach ~**0.44** (nicotine episodic memory 0.44, MPH recall 0.43). Most effects are 0.1-0.3. There is NO large, clean, replicated enhancement in healthy young adults (Roberts 2020; Heishman 2010; Klove 2025).

## Does any computable predictor forecast healthy-adult enhancement?

| predictor | what it is | AUROC | perm p |
|---|---|---|---|
| acute CNS stimulant gate | supergroup == stimulant (coarse) | **0.86** | 0.0454 |
| mechanism-class prognostic prior | the DISEASE manuscript's AUROC-1.00 winner | **0.55** | 0.4183 |
| SMD magnitude | does a bigger effect size predict the binary | 0.80 | n/a |

**The headline finding.** In disease pivotal trials the mechanism-class prognostic prior separated SUCCESS from FAILURE perfectly (AUROC 1.00) because the classes were outcome-PURE. Against the healthy-adult ground truth it COLLAPSES (AUROC 0.55): the classes are NOT pure. The decisive case: **d-amphetamine and methylphenidate are the SAME mechanism class (catecholaminergic) with the SAME overall SMD (0.21), yet methylphenidate enhances and d-amphetamine is null** (Roberts 2020); and caffeine enhances while guarana (also adenosinergic) does not. So the predictor that dominated disease-trial prediction does NOT transfer to the question this project actually exists to answer.

Impure multi-member classes (the homogeneity break): adenosinergic (1 enhance / 1 null), catecholaminergic (2 enhance / 1 null).

The only separator is the COARSE 'acute CNS stimulant' gate (AUROC 0.86): every clean enhancer is an acute monoaminergic / adenosinergic / cholinergic stimulant (MPH, modafinil, nicotine, caffeine) and every non-stimulant with a clean healthy-young MA is NULL (ginkgo, bacopa, omega-3, ginseng, creatine-young, L-theanine). But the gate is NECESSARY-not-sufficient: d-amphetamine and guarana are stimulants that do nothing, so even 'is it a stimulant' mis-ranks them. SMD magnitude scores AUROC 0.80, but that is partly MECHANICAL (the binary is defined by the effect being non-zero, so a larger point estimate trivially tracks CI-excludes-0) and it still fails the decisive case: the null d-amphetamine has the SAME 0.21 as the enhancing methylphenidate, so it is not an external predictor of which compound works.

## What this means for the project's goal

Predicting which compound enhances cognition in HEALTHY adults is HARDER than predicting disease-trial success, and for a principled reason: the disease result was a class-homogeneity look-up, and that homogeneity does not exist in the healthy-adult data. Concretely:
1. There is **no validated fine-grained (mechanism-class, target-affinity, target-genetics, or persistence-window) predictor** of novel healthy-adult cognitive enhancement. The persistence/psychoplastogen axis (PERSEUS) is orthogonal here - none of its flagged compounds even has a healthy-adult cognition meta-analysis.
2. The honest deliverable is a **calibrated negative with a low ceiling**: the real enhancer set is small, mechanistically narrow (acute CNS stimulants), and capped at SMD ~0.4; everything else with clean evidence is null, and most candidate nootropics have NO healthy-adult evidence at all.
3. The actionable, honest output for a healthy-adult screen is therefore a **coarse stimulant-class gate + an explicit abstain-by-default**, not a fine per-compound prediction the data cannot support. This mirrors, and sharpens, the disease manuscript's lead-with-the-warning posture.

**Integrity.** Every SMD is a verified meta-analytic fact (no fabrication); UNVERIFIED CIs are blank in the ledger; ABSENT rows are retained because the absence of a healthy-adult meta-analysis is itself load-bearing ground truth. n is small (the field's real clean evidence base is small); results are descriptive contrasts, not a fitted model.
