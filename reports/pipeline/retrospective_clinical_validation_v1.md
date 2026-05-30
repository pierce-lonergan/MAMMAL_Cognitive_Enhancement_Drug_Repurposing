# Retrospective Clinical-Outcome Validation v1 (Gap 3)

**Does the pipeline anticipate which cognition drugs succeed — and which fail in Phase III — without ever being told the outcome?**

This is a leakage-audited retrospective benchmark on a curated ledger of **31 cognition drugs** (13 approved/positive SUCCESS, 18 Phase II/III FAILURE) across 11 mechanism classes. Every outcome label is a documented, adjudicated pivotal-trial result (`data/raw/clinical_outcomes_ledger.csv`).

## Pre-registered analysis plan

Three predictors of SUCCESS-vs-FAILURE, ranked by information used; primary metric = AUROC with 90% bootstrap CI + label-permutation p. Pre-specified hypothesis: **target affinity + genetic relevance do NOT discriminate clinical outcome (P1 ≈ chance); mechanism-class track record does (P2 ≫ chance); extrapolation to an unseen class is hard (P3 ≈ chance).** Failure-recall threshold fixed a priori at the global mean clinical g = 0.179.

## Results

| Predictor | n | AUROC | 90% CI | perm p | Spearman(ĝ, g) |
|---|---|---|---|---|---|
| P1a target relevance σ(θ̄) | 26 | **0.589** | [0.38, 0.79] | 0.2250 | +0.07 |
| P1b within-target binding percentile | 10 | **0.125** | [0.00, 0.38] | 0.9566 | -0.21 |
| P2 class-structure leave-one-compound-out | 31 | **1.000** | [1.00, 1.00] | 0.0002 | +0.82 |
| P3 leave-one-class-out (extrapolation bound) | 31 | **0.000** | [0.00, 0.00] | 1.0000 | -0.87 |

![ROC](figures/v11/retrospective_roc.png)

## The headline — read this carefully, not as a leaderboard number

The single empirical fact this benchmark surfaces:

> **Across these 31 real cognition drugs, mechanism class perfectly stratifies clinical outcome.** Every cholinergic / catecholaminergic / wake-promoting / NMDA / multimodal-5HT drug succeeded; every α7-nAChR, 5-HT6, mGluR, AMPA-PAM, PDE9/10 and H3-cognition drug failed. There is **zero within-class outcome variance** in this ledger.

That fact is what produces the numbers below — and it is the point. It is sobering, real, and well known to the field (the AChE inhibitors remain the only broadly-approved AD cognition drugs; the α7/5-HT6 graveyard is infamous).

- **Mechanism-class track record (P2): AUROC = 1.00** [1.00, 1.00], permutation p = 0.0002; **failure-recall 100%** — flagged 9/9 of the famous Phase III failures it was never told about (encenicline, idalopirdine, intepirdine, pomaglumetad, PF-04447943, SUVN-502, ABT-126, TC-5619, MK-0249). AUROC = 1.0 is **not a predictive miracle** — it is the direct readout of the class-homogeneity above: leave-one-compound-out retains the class via siblings, so it simply recovers the class verdict. The honest content is the *contrast* with P1, not the magnitude of P2.
- **Target genetic-relevance (P1a): AUROC = 0.59** [0.38, 0.79], p = 0.22 — **at chance**. A target being cognition-relevant (high θ̄) does NOT predict that a drug hitting it will succeed: encenicline binds α7 (θ̄ high) and failed; donepezil hits ACHE (θ̄ high) and succeeded.
- **Target binding affinity (P1b): AUROC = 0.12** (n=10) — at or *below* chance. If anything strong binders failed more often, because excellent affinity is exactly what carried the doomed compounds into Phase III. Binding potency is not prognostic of cognition-trial success.
- **Leave-one-class-out (P3): AUROC = 0.00** — when the held-out drug's entire mechanism class is removed, prediction collapses (here it inverts, driven by the failure-weighted base rate). The pipeline triages **within known mechanism space**; it cannot forecast an unseen mechanism. This is the honest ceiling on the claim.

**The defensible scientific claim** (not the AUROC=1.0 number): *in cognition drug development, the prognostic signal lives in the clinical track record of the mechanism class — target-binding affinity and target genetic-relevance, the two things a target-first in-silico pipeline measures, are at chance.* This is the same lesson as the V6.B Gate-2 falsification, now demonstrated directly against pivotal-trial outcomes, and it is the empirical case for why a cognition repurposing pipeline must be class-aware and phenotype-aware rather than target-affinity-driven.

## Leakage audit (per predictor)

- **P1a target relevance σ(θ̄)** — V6.B θ̄ built from GWAS/AHBA/single-cell brain data — never saw cognition trials
- **P1b within-target binding percentile** — V6.A MAMMAL/MMAtt pKd from ChEMBL bioactivity — never saw cognition trials
- **P2 class-structure leave-one-compound-out** — uses mechanism-class SIBLINGS' meta-analytic g; the held-out drug's OWN outcome is excluded (legitimate inductive generalization)
- **P3 leave-one-class-out (extrapolation bound)** — predicts from all OTHER mechanism classes; the drug's own class is removed entirely

No predictor uses the held-out drug's own trial outcome. P1a/P1b use only GWAS/expression (V6.B) and ChEMBL bioactivity (V6.A), which are structurally independent of cognition-trial readouts. P2 uses siblings' meta-analytic g — legitimate inductive generalization, never the held-out drug itself.

## Per-drug predictions (P2 leave-one-compound-out)

| Drug | Class | Indication | Actual | ĝ (P2) | Predicted | ✓ |
|---|---|---|---|---|---|---|
| TC-5619 | alpha7_nAChR | schizophrenia/ADHD | FAILURE | +0.057 | FAILURE | ✅ |
| TAK-063 | PDE9_PDE10 | schizophrenia | FAILURE | +0.043 | FAILURE | ✅ |
| BI-409306 | PDE9_PDE10 | schizophrenia/AD | FAILURE | +0.043 | FAILURE | ✅ |
| PF-04447943 | PDE9_PDE10 | AD | FAILURE | +0.060 | FAILURE | ✅ |
| farampator | AMPA_PAM | healthy/AD | FAILURE | +0.070 | FAILURE | ✅ |
| S47445 | AMPA_PAM | AD | FAILURE | +0.093 | FAILURE | ✅ |
| CX-516 | AMPA_PAM | schizophrenia/MCI | FAILURE | +0.070 | FAILURE | ✅ |
| mavoglurant | mGluR | FXS | FAILURE | +0.060 | FAILURE | ✅ |
| basimglurant | mGluR | FXS | FAILURE | +0.043 | FAILURE | ✅ |
| pomaglumetad | mGluR | schizophrenia | FAILURE | +0.043 | FAILURE | ✅ |
| SUVN-502 | 5HT6_antagonist | AD | FAILURE | +0.043 | FAILURE | ✅ |
| intepirdine | 5HT6_antagonist | AD | FAILURE | +0.026 | FAILURE | ✅ |
| idalopirdine | 5HT6_antagonist | AD | FAILURE | +0.043 | FAILURE | ✅ |
| DMXB-A | alpha7_nAChR | schizophrenia | FAILURE | +0.052 | FAILURE | ✅ |
| MK-0249 | H3_cognition | AD/schizophrenia | FAILURE | +0.089 | FAILURE | ✅ |
| ABT-288 | H3_cognition | schizophrenia | FAILURE | +0.089 | FAILURE | ✅ |
| encenicline | alpha7_nAChR | CIAS-schizophrenia | FAILURE | +0.052 | FAILURE | ✅ |
| ABT-126 | alpha7_nAChR | AD/schizophrenia | FAILURE | +0.040 | FAILURE | ✅ |
| vortioxetine | multimodal_5HT | MDD-cognition | SUCCESS | +0.179 | SUCCESS | ✅ |
| pitolisant | wake_promoting | narcolepsy-EDS | SUCCESS | +0.310 | SUCCESS | ✅ |
| armodafinil | wake_promoting | narcolepsy/SWD | SUCCESS | +0.396 | SUCCESS | ✅ |
| modafinil | wake_promoting | narcolepsy-EDS | SUCCESS | +0.380 | SUCCESS | ✅ |
| guanfacine-XR | catecholaminergic_ADHD | ADHD | SUCCESS | +0.436 | SUCCESS | ✅ |
| atomoxetine | catecholaminergic_ADHD | ADHD | SUCCESS | +0.446 | SUCCESS | ✅ |
| dextroamphetamine | catecholaminergic_ADHD | ADHD | SUCCESS | +0.436 | SUCCESS | ✅ |
| lisdexamfetamine | catecholaminergic_ADHD | ADHD | SUCCESS | +0.426 | SUCCESS | ✅ |
| methylphenidate | catecholaminergic_ADHD | ADHD | SUCCESS | +0.436 | SUCCESS | ✅ |
| memantine | NMDA_modulator | AD-mod-sev | SUCCESS | +0.179 | SUCCESS | ✅ |
| rivastigmine | AChE_inhibitor | AD | SUCCESS | +0.303 | SUCCESS | ✅ |
| galantamine | AChE_inhibitor | AD | SUCCESS | +0.313 | SUCCESS | ✅ |
| donepezil | AChE_inhibitor | AD | SUCCESS | +0.316 | SUCCESS | ✅ |

## Honest limitations

- **Class-outcome homogeneity drives P2.** In this ledger every mechanism class is outcome-homogeneous (all AChE-I/stimulant/wake succeed; all α7/5-HT6/mGluR/AMPA-PAM/PDE/H3-cognition fail). P2's high AUROC reflects this real homogeneity — it is the finding, not a trick — but it means P2 cannot resolve *within-class* winners from losers, and P3 shows it cannot extrapolate to a new class.
- **n is small (31).** AUROCs carry wide CIs; permutation p guards against chance, but this is a proof-of-principle benchmark, not a definitive estimate.
- **Indication matters.** H3 antagonism succeeds for narcolepsy EDS (pitolisant) but fails for AD/schizophrenia cognition (MK-0249) — encoded as distinct classes. The ledger is cognition-endpoint-focused.
- **The ledger is curated, not exhaustive.** It is a balanced, documented sample of the canonical cognition successes and the famous failures; expansion is a follow-up.

---

Generated by `scripts/75_retrospective_clinical_validation.py` via `mammal_repurposing.validation.retrospective`. Truth set: `data/raw/clinical_outcomes_ledger.csv`.