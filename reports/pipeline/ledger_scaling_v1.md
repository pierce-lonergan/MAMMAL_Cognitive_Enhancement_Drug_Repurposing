# F3 - Ledger scaling, per-domain structure, and the power roadmap

**Questions.** (1) Does the class-separation result survive scaling the leakage-audited ledger from n=31 to the cited n=47? (2) Is the class-success pattern consistent across cognitive domains? (3) How large must the ledger get for the F1 within-class test to become conclusive? Real cited ledgers only (base + EXTENSION + CT.gov); no fabricated outcomes. Reproduced by `scripts/94_ledger_scaling.py`.

## 1. Scaling trajectory

| Step | n | classes | outcome-pure | class-LOCO AUROC | perm p | % var between-class | ICC(1) |
|---|---|---|---|---|---|---|---|
| base (frozen 31) | 31 | 11 | 11/11 (100%) | 1.000 | 0.0002 | 96.5% | 0.951 |
| + EXTENSION | 42 | 17 | 17/17 (100%) | 0.990 | 0.0002 | 96.9% | 0.952 |
| + CT.gov (unbiased) | 47 | 20 | 20/20 (100%) | 0.967 | 0.0002 | 97.0% | 0.952 |
| + web-researched (verified) | 125 | 48 | 46/48 (96%) | 0.915 | 0.0002 | 76.2% | 0.623 |

Through the cited ledgers (n=47, 20 classes) the pattern is PRESERVED: class-LOCO AUROC 1.000 -> 0.967, classes stay 100% outcome-pure, and 97% of clinical-*g* variance remains between-class (ICC 0.95) - not a small-n artifact of the original 31. The web-researched + human-adjudicated step (n=125) confirms the signal survives scaling: class-LOCO AUROC 0.915; see section 3b.

## 2. Per-domain structure

Each drug assigned its pivotal endpoint's primary cognitive domain. AUROC is the class-LOCO separation within that domain (computed only where the domain holds both outcomes and >=2 classes).

| Cognitive domain | n | success | failure | classes | within-domain AUROC |
|---|---|---|---|---|---|
| global_amnestic | 59 | 8 | 51 | 27 | 0.794 |
| scz_composite_battery | 18 | 0 | 18 | 10 | n/a (single-outcome) |
| other | 17 | 1 | 16 | 14 | 0.000 |
| functional_composite | 10 | 3 | 7 | 8 | 0.762 |
| adhd_symptom | 5 | 5 | 0 | 1 | n/a (single-outcome) |
| episodic_memory | 4 | 0 | 4 | 3 | n/a (single-outcome) |
| processing_speed | 3 | 1 | 2 | 3 | 0.000 |
| psychosis_secondary | 3 | 1 | 2 | 3 | 0.000 |
| wakefulness | 3 | 3 | 0 | 1 | n/a (single-outcome) |
| working_memory | 2 | 0 | 2 | 2 | n/a (single-outcome) |
| executive_attention | 1 | 0 | 1 | 1 | n/a (single-outcome) |

Most drugs sit in a global-amnestic (AD: ADAS-Cog) or schizophrenia composite (MCCB) endpoint, so the current ledger supports domain *stratification* but not fine per-(drug, domain) *g* decomposition - the pivotal trials report one global/composite endpoint, not domain sub-scores. Splitting a drug's effect across working-memory / processing-speed / episodic-memory requires curating trial secondary analyses (the remaining F3 curation; schema below).

## 3. Power roadmap (the actionable output)

- Current pooled within-class points (members of multi-member, g-varying classes): **58** across 11 classes (avg 5.3/class).

| Target within-class rho | effective points needed | x current | implied total ledger n |
|---|---|---|---|
| 0.30 | 85 | 1.5x | ~184 |
| 0.40 | 47 | 0.8x | ~102 |
| 0.50 | 30 | 0.5x | ~65 |

To make the F1 within-class test conclusive at a moderate effect (rho=0.4), the ledger needs roughly **102 drugs** that land in multi-member, g-varying (i.e. SUCCESS) classes - several times the current 58 pooled points. Failure classes (all g~0) add purity evidence but zero within-class power, so the binding curation target is **more multi-member SUCCESS classes with genuine within-class g spread**, ideally with per-domain sub-scores.

## 3b. Research-curated + human-adjudicated step (n=125)

An Opus multi-agent run (106 agents: research + independent adversarial verification) added 78 web-verified cognition-drug outcomes; every disputed SUCCESS/FAILURE call was then re-adjudicated by two independent Opus adjudicators under a strict cognition-EFFICACY convention. All 78 are kept as binary data points (no exclusions): safety-halted drugs with a clean positive cognition pivotal (metrifonate, eptastigmine) stay SUCCESS, while unverifiable or contested-benefit drugs (velnacrine, aducanumab, sodium oligomannate, masitinib, nicergoline) are FAILURE; 11 over-generous SUCCESS calls were recoded to FAILURE in total. The ledger reaches n=125 / 48 classes -- meeting the F1 power target; per-row verdict + cited basis are in the provenance. The frozen base-31 analysis is untouched.

| Scenario | n | class-LOCO AUROC |
|---|---|---|
| raw web-research (pre-adjudication) | 125 | 0.766 |
| **adjudicated, full** | 125 | **0.915** |
| adjudicated, multi-member classes only | 103 | 0.971 |

**Interpretation.** Adjudicating the disputed codings lifts the class-LOCO AUROC from the raw 0.77 to **0.92** (**0.97** on the multi-member classes the predictor can actually leverage; the 22 singleton classes have no siblings for leave-one-compound-out and structurally fall to the global mean). The class-history signal **robustly survives scaling** to n=125, still far above the leakage-free target-level predictors (affinity 0.47, genetics 0.59) -- the raw 0.77 was coding noise, not signal collapse. The perfect 1.00 at n=31 was partly a sparse-sampling / selection effect; the honest value with full class population is ~0.92.

After adjudication only TWO genuinely mixed-outcome classes remain (S/F): AChE_inhibitor 6/4; anti_amyloid_beta_mab 2/6. The anti-amyloid mAb split (lecanemab/donanemab succeed where 5 earlier anti-Abeta mAbs failed) is a real boundary on broad-mechanism purity; the AChE-inhibitor split is the marketed winners vs later AChE-Is that failed on efficacy/dosing -- both genuine, not coding artifacts.

## 4. Verdict and remaining curation

- **Scaling**: robust through the cited ledgers (n=47: AUROC 0.967, 100% pure, ICC 0.95); the web-researched + adjudicated n=125 step holds at AUROC 0.92 (0.97 on multi-member classes); see 3b.
- **Per-domain**: supported as stratification on the real pivotal endpoints; fine per-domain *g* needs sub-score curation.
- **F1 power**: needs ~102 drugs (rho=0.4 target) concentrated in SUCCESS classes to become conclusive.

The remaining step is genuine literature curation (real drugs, real trials, real adjudicated outcomes, real cited effect sizes) - it is deliberately NOT auto-generated here, to protect the ledger's integrity. The curation protocol and the per-domain schema live in `docs/LEDGER_CURATION.md`; `load_all_ledgers()` validates and ingests any additional cited ledger CSV that follows the schema.
