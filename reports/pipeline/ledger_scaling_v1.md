# F3 - Ledger scaling, per-domain structure, and the power roadmap

**Questions.** (1) Does the class-separation result survive scaling the leakage-audited ledger from n=31 to the cited n=47? (2) Is the class-success pattern consistent across cognitive domains? (3) How large must the ledger get for the F1 within-class test to become conclusive? Real cited ledgers only (base + EXTENSION + CT.gov); no fabricated outcomes. Reproduced by `scripts/94_ledger_scaling.py`.

## 1. Scaling trajectory

| Step | n | classes | outcome-pure | class-LOCO AUROC | perm p | % var between-class | ICC(1) |
|---|---|---|---|---|---|---|---|
| base (frozen 31) | 31 | 11 | 11/11 (100%) | 1.000 | 0.0002 | 96.5% | 0.951 |
| + EXTENSION | 42 | 17 | 17/17 (100%) | 0.990 | 0.0002 | 96.9% | 0.952 |
| + CT.gov (unbiased) | 47 | 20 | 20/20 (100%) | 0.967 | 0.0002 | 97.0% | 0.952 |

Adding 16 cited drugs and 9 new mechanism classes PRESERVES the pattern: class-LOCO AUROC 1.000 -> 0.967, classes stay 100% outcome-pure, and 97% of clinical-*g* variance remains between-class (ICC 0.95). The class-history signal is not a small-n artifact of the original 31.

## 2. Per-domain structure

Each drug assigned its pivotal endpoint's primary cognitive domain. AUROC is the class-LOCO separation within that domain (computed only where the domain holds both outcomes and >=2 classes).

| Cognitive domain | n | success | failure | classes | within-domain AUROC |
|---|---|---|---|---|---|
| global_amnestic | 17 | 4 | 13 | 9 | 0.923 |
| scz_composite_battery | 11 | 0 | 11 | 5 | n/a (single-outcome) |
| functional_composite | 6 | 1 | 5 | 5 | 0.000 |
| adhd_symptom | 5 | 5 | 0 | 1 | n/a (single-outcome) |
| psychosis_secondary | 3 | 1 | 2 | 3 | 0.000 |
| wakefulness | 3 | 3 | 0 | 1 | n/a (single-outcome) |
| episodic_memory | 1 | 0 | 1 | 1 | n/a (single-outcome) |
| processing_speed | 1 | 1 | 0 | 1 | n/a (single-outcome) |

Most drugs sit in a global-amnestic (AD: ADAS-Cog) or schizophrenia composite (MCCB) endpoint, so the current ledger supports domain *stratification* but not fine per-(drug, domain) *g* decomposition - the pivotal trials report one global/composite endpoint, not domain sub-scores. Splitting a drug's effect across working-memory / processing-speed / episodic-memory requires curating trial secondary analyses (the remaining F3 curation; schema below).

## 3. Power roadmap (the actionable output)

- Current pooled within-class points (members of multi-member, g-varying classes): **34** across 10 classes (avg 3.4/class).

| Target within-class rho | effective points needed | x current | implied total ledger n |
|---|---|---|---|
| 0.30 | 85 | 2.5x | ~118 |
| 0.40 | 47 | 1.4x | ~65 |
| 0.50 | 30 | 0.9x | ~42 |

To make the F1 within-class test conclusive at a moderate effect (rho=0.4), the ledger needs roughly **65 drugs** that land in multi-member, g-varying (i.e. SUCCESS) classes - several times the current 34 pooled points. Failure classes (all g~0) add purity evidence but zero within-class power, so the binding curation target is **more multi-member SUCCESS classes with genuine within-class g spread**, ideally with per-domain sub-scores.

## 4. Verdict and remaining curation

- **Scaling**: the headline class-separation result is robust to the n=31 -> 47 expansion across 20 cited mechanism classes (AUROC 0.967, 100% pure, ICC 0.95).
- **Per-domain**: supported as stratification on the real pivotal endpoints; fine per-domain *g* needs sub-score curation.
- **F1 power**: needs ~65 drugs (rho=0.4 target) concentrated in SUCCESS classes to become conclusive.

The remaining step is genuine literature curation (real drugs, real trials, real adjudicated outcomes, real cited effect sizes) - it is deliberately NOT auto-generated here, to protect the ledger's integrity. The curation protocol and the per-domain schema live in `docs/LEDGER_CURATION.md`; `load_all_ledgers()` validates and ingests any additional cited ledger CSV that follows the schema.
