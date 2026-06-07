# F3 - Ledger scaling, per-domain structure, and the power roadmap

**Questions.** (1) Does the class-separation result survive scaling the leakage-audited ledger from n=31 to the cited n=47? (2) Is the class-success pattern consistent across cognitive domains? (3) How large must the ledger get for the F1 within-class test to become conclusive? Real cited ledgers only (base + EXTENSION + CT.gov); no fabricated outcomes. Reproduced by `scripts/94_ledger_scaling.py`.

## 1. Scaling trajectory

| Step | n | classes | outcome-pure | class-LOCO AUROC | perm p | % var between-class | ICC(1) |
|---|---|---|---|---|---|---|---|
| base (frozen 31) | 31 | 11 | 11/11 (100%) | 1.000 | 0.0002 | 96.5% | 0.951 |
| + EXTENSION | 42 | 17 | 17/17 (100%) | 0.990 | 0.0002 | 96.9% | 0.952 |
| + CT.gov (unbiased) | 47 | 20 | 20/20 (100%) | 0.967 | 0.0002 | 97.0% | 0.952 |
| + web-researched (verified) | 125 | 49 | 40/49 (82%) | 0.766 | 0.0002 | 66.3% | 0.456 |

Through the cited ledgers (n=47, 20 classes) the pattern is PRESERVED: class-LOCO AUROC 1.000 -> 0.967, classes stay 100% outcome-pure, and 97% of clinical-*g* variance remains between-class (ICC 0.95) - not a small-n artifact of the original 31. The web-researched step (n=125) is RESEARCH-GRADE and shows the pattern is scale-sensitive once classes are fully populated (raw AUROC 0.766); see section 3b for the sensitivity decomposition.

## 2. Per-domain structure

Each drug assigned its pivotal endpoint's primary cognitive domain. AUROC is the class-LOCO separation within that domain (computed only where the domain holds both outcomes and >=2 classes).

| Cognitive domain | n | success | failure | classes | within-domain AUROC |
|---|---|---|---|---|---|
| global_amnestic | 59 | 14 | 45 | 27 | 0.679 |
| scz_composite_battery | 18 | 1 | 17 | 10 | 0.059 |
| other | 14 | 4 | 10 | 13 | 0.000 |
| functional_composite | 10 | 3 | 7 | 8 | 0.762 |
| wakefulness | 6 | 3 | 3 | 3 | 1.000 |
| adhd_symptom | 5 | 5 | 0 | 1 | n/a (single-outcome) |
| episodic_memory | 4 | 1 | 3 | 3 | 0.000 |
| processing_speed | 3 | 1 | 2 | 3 | 0.000 |
| psychosis_secondary | 3 | 1 | 2 | 3 | 0.000 |
| working_memory | 2 | 0 | 2 | 2 | n/a (single-outcome) |
| executive_attention | 1 | 0 | 1 | 1 | n/a (single-outcome) |

Most drugs sit in a global-amnestic (AD: ADAS-Cog) or schizophrenia composite (MCCB) endpoint, so the current ledger supports domain *stratification* but not fine per-(drug, domain) *g* decomposition - the pivotal trials report one global/composite endpoint, not domain sub-scores. Splitting a drug's effect across working-memory / processing-speed / episodic-memory requires curating trial secondary analyses (the remaining F3 curation; schema below).

## 3. Power roadmap (the actionable output)

- Current pooled within-class points (members of multi-member, g-varying classes): **74** across 18 classes (avg 4.1/class).

| Target within-class rho | effective points needed | x current | implied total ledger n |
|---|---|---|---|
| 0.30 | 85 | 1.1x | ~144 |
| 0.40 | 47 | 0.6x | ~80 |
| 0.50 | 30 | 0.4x | ~51 |

To make the F1 within-class test conclusive at a moderate effect (rho=0.4), the ledger needs roughly **80 drugs** that land in multi-member, g-varying (i.e. SUCCESS) classes - several times the current 74 pooled points. Failure classes (all g~0) add purity evidence but zero within-class power, so the binding curation target is **more multi-member SUCCESS classes with genuine within-class g spread**, ideally with per-domain sub-scores.

## 3b. Research-batch scaling sensitivity (n=125, RESEARCH-GRADE)

The web-researched batch (independently existence-verified) takes the ledger to n=125 across 49 classes. It is RESEARCH-GRADE: the SUCCESS/FAILURE boundary for 14 old/controversial drugs is genuinely disputed (flagged in the provenance) and the agents' class vocabulary was harmonized. The frozen base-31 + EXTENSION + CT.gov analysis is unchanged; this is a sensitivity probe, not a headline.

| Scenario | n | class-LOCO AUROC |
|---|---|---|
| full (raw research-grade) | 125 | 0.766 |
| multi-member classes only | 103 | 0.808 |
| borderline successes -> FAILURE (conservative) | 125 | 0.906 |
| borderline successes dropped | 111 | 0.916 |

**Interpretation.** The raw AUROC falls to 0.77, but the drop is mostly the 14 controversial SUCCESS codings: under conservative handling the class signal holds at ~0.91 - still far above the leakage-free target-level predictors (affinity 0.47, genetics 0.59). 22 singleton classes (no siblings for leave-one-compound-out) add a further structural ~0.04. The genuinely robust mixed-outcome class is anti-amyloid mAbs (lecanemab/donanemab succeed where earlier anti-Abeta mAbs failed) - a real boundary on broad-mechanism purity. Net: class-history prognosis substantially SURVIVES scaling (~0.91 at n=125 under conservative coding), but the perfect 1.00 at n=31 was partly a sparse-sampling / selection effect.

Mixed-outcome classes at this n (S/F): AChE_inhibitor 7/3; AMPA_PAM 1/8; DAAO_inhibitor 1/1; GSK3_inhibitor 1/1; H3_cognition 1/2; amyloid_aggregation_inhibitor 1/2; anti_amyloid_beta_mab 3/5; serotonin_5HT1A_partial_agonist 1/2; tyrosine_kinase_inhibitor 1/2.

These rows require human adjudication (esp. the borderline successes and the AChE safety-vs-efficacy failures) before informing any published claim; see `data/raw/clinical_outcomes_ledger_RESEARCH_provenance.csv`.

## 4. Verdict and remaining curation

- **Scaling**: robust through the cited ledgers (n=47: AUROC 0.967, 100% pure, ICC 0.95). The research-grade n=125 step shows scale-sensitivity (raw 0.77; ~0.91 under conservative coding); see 3b.
- **Per-domain**: supported as stratification on the real pivotal endpoints; fine per-domain *g* needs sub-score curation.
- **F1 power**: needs ~80 drugs (rho=0.4 target) concentrated in SUCCESS classes to become conclusive.

The remaining step is genuine literature curation (real drugs, real trials, real adjudicated outcomes, real cited effect sizes) - it is deliberately NOT auto-generated here, to protect the ledger's integrity. The curation protocol and the per-domain schema live in `docs/LEDGER_CURATION.md`; `load_all_ledgers()` validates and ingests any additional cited ledger CSV that follows the schema.
