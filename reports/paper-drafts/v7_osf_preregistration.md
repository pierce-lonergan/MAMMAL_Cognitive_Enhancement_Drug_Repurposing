# V7 — OSF Pre-Registration: Clinical Effect-Size Translation Function

**Project**: MAMMAL Cognitive Enhancement Drug Repurposing (V7 layer)
**Lead**: Pierce Lonergan
**Target venue**: *Clinical Pharmacology & Therapeutics* (Wiley, IF 7.3); fallback *CPT: Pharmacometrics & Systems Pharmacology* (CPT:PSP, IF 4.2)
**OSF lock date**: TBD (lock BEFORE unblinding the held-out anchor set)
**Pre-registration template**: OSF.io + AsPredicted.org
**Companion design doc**: `design/architecture-and-plans/V4_STATUS_AND_FORWARD_PLAN.md` §13.Y
**Research source docs**:
- `research/4-tier/Clinical Effect-Size Translation Function.md` (V7 spec)
- `research/4-tier/Clinical Effect-Size Translation Function A Methodology Pre-Registration for Bayesian Cognition-Enhancement Drug Repurposing.md` (V7 companion)

---

## 1. Hypothesis

V7 translates the in-silico stack (V6.A pchembl posteriors + V6.B Cluster D θ̄ posteriors + V7.1 PBPK exposure) into a predicted healthy-adult cognition Hedges' *g* with credible intervals, calibrated against the Roberts CA, Jones A, Sumnall H, Gage SH, Montgomery C 2020 *Eur Neuropsychopharmacol* 38:40-62 ceiling.

**Primary hypothesis**: a 3-level hierarchical Bayesian model with PRISMA-anchored 12-class meta-analytic priors, a Cluster D multiplicative gate `β_target[t_c] = θ̄_{t_c} · β_raw_target[t_c]`, and a 5-failure-mode moderator (m1 U-shape, m2 practice/placebo, m3 tolerance, m4 trait×state, m5 trial-design) predicts published healthy-adult Hedges' *g* with MAE < 0.15 on a 15-compound held-out anchor set.

**Secondary hypothesis**: no compound's posterior 90% credible upper bound for *g* exceeds 0.50 (the Roberts 2020 ceiling) — Gate 2 (HARD).

**Tertiary hypothesis**: the model correctly recapitulates failure of α7 partial agonists (encenicline) at high desensitisation-prone doses via the tolerance-kinetics moderator m3 — P2 of the 8 pre-registered predictions.

---

## 2. Pre-registered model specification

### 2.1 Three-level hierarchy

```
μ_global             ~ Normal(0, 0.20)
μ_class[m]           ~ Normal(μ_class_PRISMA[m], λ_class · σ_class_PRISMA[m])
η[c, t]              = sigmoid(α + β1·E[pchembl_post]
                                + β2·E[relevance_post]
                                + β3·copula_correction)
                       · class_contribution[c]
                       − Σ_k γ_k · m_k[c]
g[c, t]              ~ Normal(η[c, t], σ_resid²)
```

with the **Cluster D multiplicative gate**:
```
β_target[t_c] = θ̄_{t_c} · β_raw_target[t_c]
```

### 2.2 Hyperprior choices

| Hyperprior | Distribution | Source |
|---|---|---|
| μ_global | Normal(0, 0.20) | Mild zero-prior on the population mean g |
| λ_class | sweep ∈ {0.10, 0.30, 1.00, 3.00} | Schmidli 2014 robust MAP weight |
| β_pchembl | Normal(0, 0.10) | Weak prior; data-driven β-pchembl |
| β_relevance | Normal(0, 0.10) | Weak prior; θ̄ contribution scaled |
| β_copula | Normal(0, 0.05) | Gaussian-copula correction for pchembl × relevance |
| γ_k (5 moderators) | Normal(0, 0.10) per k | Weak; data identifies which moderators matter |
| σ_resid | HalfNormal(0.20) | Residual SD ≤ Roberts ceiling/2 |

### 2.3 PRISMA 12-class meta-analytic priors

**LOCKED**: the 12 mechanism classes + their (mean, sd, n_trials, n_subjects, peak_subdomain_g, representative_drug, citation) are frozen in `src/mammal_repurposing/translation/prisma_priors.py::PRISMA_CLASS_PRIORS` at commit time. Classes are:

1. **AChE-I** (donepezil, n=8 trials, prior_g=0.18, peak=0.31 on delayed recall)
2. **wake_promoting** (modafinil, n=14, prior_g=0.12, peak=0.30 on vigilance)
3. **NDRI** (methylphenidate, n=12, prior_g=0.21, peak=0.43 on delayed recall)
4. **NRI** (atomoxetine, n=6, prior_g=0.10, peak=0.20 on response inhibition)
5. **NMDA_antagonist** (memantine, n=5, prior_g=0.05, peak=0.15 on learning)
6. **multimodal_5HT** (vortioxetine, n=4, prior_g=0.12, peak=0.25 on processing speed)
7. **alpha2A_agonist** (guanfacine, n=5, prior_g=0.15, peak=0.28 on working memory)
8. **A2A_antagonist** (caffeine, n=22, prior_g=0.20, peak=0.40 on vigilance)
9. **AMPA_pos_mod** (piracetam, n=6, prior_g=0.05, peak=0.15 on declarative memory)
10. **creatine** (n=10, prior_g=0.08, peak=0.20 in low-baseline subgroups)
11. **omega3** (EPA_DHA, n=18, prior_g=0.07, peak=0.15 on episodic memory)
12. **minocycline** (n=4, prior_g=0.05, peak=0.10 on working memory)

All 12 classes have `peak_subdomain_g ≤ 0.50` (Roberts ceiling sanity check at the prior level).

### 2.4 5 failure-mode moderators

```
m1 = U-shape miss          (dose past peak; from PBPK u_shape_occupancy)
m2 = practice/placebo      (trial design moderator)
m3 = tolerance onset       (chronic vs acute mismatch; from R_avail dynamics)
m4 = trait × state         (responder enrichment masks population mean)
m5 = trial-design          (parallel-group vs crossover, endpoint sensitivity)
```

Each moderator: binary {0, 1} flag OR continuous [0, 1] score per (compound, trial).

### 2.5 PBPK 9-compartment

**LOCKED**: 9 compartments (gut/plasma/peripheral/cortex/striatum/hippocampus/basal-forebrain/brainstem/CSF) per `src/mammal_repurposing/translation/pbpk.py::COMPARTMENTS`. PET anchors:

| Drug | Receptor | Compartment | Expected peak occupancy | Reference |
|---|---|---|---|---|
| Donepezil | AChE | cortex | 19.1% | Bohnen 2005 *Neurology* 64:1037 |
| MPH 20mg | DAT | striatum | 54% | Volkow 1998 *Am J Psych* 155:1325 |
| Haloperidol 2mg | D2 | striatum | 65% | Kapur 2000 *Am J Psych* 157:514 |

### 2.6 Sampling configuration

- Sampler: PyMC NUTS via numpyro JAX backend (preferred) or PyMC default
- Chains: 4
- Tune: 2000
- Draws: 2000
- target_accept: 0.95
- Random seed: 42

---

## 3. Eight pre-registered predictions (P1–P8)

| # | Compound | Endpoint | Pre-registered band | Source |
|---|---|---|---|---|
| P1 | Donepezil | ADAS-Cog | g ∈ [0.10, 0.30] | Birks 2018 Cochrane CD001190 |
| P2 | Encenicline 3mg | MCCB | |g| < 0.20 (Phase 3 failure recapitulated) | Keefe 2015 *Neuropsychopharmacology* 40:3053; Brannan 2019 Phase 3 (Alzforum/FORUM Mar 2016 topline) |
| P3 | Methylphenidate 20mg | DSST | g ∈ [0.15, 0.30] | Roberts 2020 SMD=0.21 (MPH overall) |
| P4 | Modafinil 200mg | n-back | g ∈ [0.06, 0.18] | Roberts 2020 SMD=0.12 (modafinil overall) |
| P5 | Memantine 20mg | RAVLT | g ∈ [-0.05, 0.20] | Repantis 2010 healthy-adult meta |
| P6 | Intepirdine | ADAS-Cog | g ∈ [-0.10, 0.15] | Lang 2021 MINDSET (null) |
| P7 | Pridopidine | cUHDRS proxy | g ∈ [-0.10, 0.15] | Reilmann 2025 PROOF-HD (null) |
| P8 | Lecanemab | CDR-SB (cognitive subdomain) | g ∈ [0.0, 0.15] | Aβ-mAb published cognitive subdomain effects |

**Falsifier**: if ≥3 of P1–P8 fail their pre-registered bands, the V7 hierarchical translation framework is **falsified for healthy-adult cognition**. Publishable contribution then becomes the *negative-result paper* in CPT:PSP — "PBPK + receptor-occupancy hierarchical Bayes cannot translate in-silico DTI rankings to healthy-adult Hedges' g, because the Roberts 2020 ceiling is at the floor of the model's discriminative resolution."

---

## 4. Four validation gates

### 4.1 Gate 1 (HARD): P1–P8 prediction-band recovery

**Threshold**: ≥6 of 8 P1–P8 predictions land within their pre-registered bands (at least 75% PASS rate).

**Falsifier**: ≤2 of 8 PASS → Gate 1 FAIL → CPT:PSP negative-result fallback triggered.

### 4.2 Gate 2 (HARD): Roberts 2020 SMD ceiling

**Threshold**: no compound's posterior 90% credible upper bound for *g* exceeds 0.50 (the Roberts 2020 ceiling).

**Falsifier**: any compound predicts g₉₀ > 0.50 → Gate 2 FAIL → model overconfident → priors need tightening.

### 4.3 Gate 3: MAE on held-out anchor set

**Threshold**: mean absolute error on a 15-compound held-out anchor set (curated independently of training) < 0.15.

**Falsifier**: MAE > 0.25 → Gate 3 FAIL → translation function not predictive at clinically-meaningful resolution.

### 4.4 Gate 4: per-endpoint calibration

**Threshold**: per-endpoint calibration plot 90% CrI coverage ≥ 85% across 6 endpoints (ADAS-Cog, DSST, n-back, Stroop, RAVLT, CANTAB-RVIP).

**Falsifier**: coverage < 70% on ≥2 endpoints → Gate 4 FAIL → CrIs miscalibrated.

---

## 5. Sensitivity analyses (also pre-registered)

| Parameter | Sweep values | What it tests |
|---|---|---|
| λ_class | {0.10, 0.30, 1.00, 3.00} | Schmidli 2014 robust MAP weight |
| σ_resid prior | {0.10, 0.20, 0.30} | Tightness of residual prior |
| γ_k prior SD | {0.05, 0.10, 0.20} | Moderator-effect strength |
| Sub-divided PRISMA priors | with/without subgroup splits | Whether age, baseline subgroups matter |

Each sweep reports top-50 shortlist overlap; PASS if ≥60% of top-50 invariant across the sweep.

---

## 6. Held-out anchor set (locked)

**15 compounds curated from Roberts 2020 + Cochrane + MetaPsy + literature** for the Gate 3 MAE test:

Defined in `src/mammal_repurposing/cluster_d/validation_gates.py::REFERENCE_COMPOUND_SMD`:

donepezil, galantamine, rivastigmine, memantine, methylphenidate, d-amphetamine, modafinil, atomoxetine, varenicline, caffeine, encenicline, intepirdine, pridopidine, vortioxetine, guanfacine.

This set is **frozen at OSF lock time**. The held-out is computed via leave-one-out: for each anchor compound, refit V7 without that compound's published *g*, then predict and compute residual.

---

## 7. Falsifiability + pre-registration timing

OSF.io project will be locked **BEFORE**:

1. Running the V7 NUTS on the full V6.A + V6.B + V7 + V8 pipeline (any compound used for training cannot also be in the held-out)
2. Looking at any Gate-1/2/3/4 result
3. Adjusting any hyperprior in response to V7 output

The pre-registration lock includes:
- The 12-class PRISMA priors (their exact mean + sd + n_trials)
- The 5-moderator structure
- The 8 P1–P8 predictions with their bands
- The 4 validation-gate thresholds
- The 15-compound held-out anchor set
- The PBPK 9-compartment + 3 PET anchors

---

## 8. Publication plan

**Primary**: *Clinical Pharmacology & Therapeutics* (Wiley, IF 7.3) full paper.
- Title (draft): "PBPK-anchored hierarchical Bayes translation of in-silico DTI rankings to predicted healthy-adult cognition Hedges' *g*, validated against pre-registered Phase-3 outcomes including encenicline."

**Fallback**: *CPT: Pharmacometrics & Systems Pharmacology* (Wiley, IF 4.2) negative-result paper.
- Title (draft): "Healthy-adult cognition Hedges' *g* prediction from in-silico DTI rankings is bounded by the Roberts 2020 ceiling: a pre-registered falsification."

Both submissions will include:
- The OSF lock URL (proof of pre-registration)
- The 15-compound anchor set parquet
- All `src/mammal_repurposing/translation/*` code under permissive open-source license
- Replication instructions for the full V6.A → V6.B → V7 → V8 chain

---

## 9. Caveats + limitations

1. **Stub vs NUTS**: this pre-registration covers the *full PyMC NUTS path* (`fit_effect_size_nuts`). The earlier stub mode (`fit_effect_size_stub`) is a sanity-test fallback that does NOT count as a NUTS posterior. Gate 3 + 4 require NUTS draws.
2. **Anchor-set sparsity**: 15 compounds is small. Confidence intervals on Gate 3 MAE will be wide. We pre-register MAE < 0.15 as the threshold despite this; the Gate 1 P1–P8 outcome is the more discriminating test.
3. **Roberts ceiling g=0.50**: the highest single-subdomain significant effect Roberts 2020 reports is MPH delayed recall g=0.43. Overall mean SMDs are 0.12 (modafinil) – 0.21 (MPH). V7 predicting g > 0.50 with high confidence would invalidate Roberts 2020 itself, which is implausible.
4. **No new wet-lab data**: V7 is entirely in-silico. The Hedges' *g* predictions are made against the published meta-analytic literature; no prospective trial validation.
5. **PRISMA prior curation may evolve**: the 12-class priors are V7.2 Stage 1. Stage 2 will refine per-subdomain (working memory, processing speed, attention, declarative memory) breakdowns. Any post-lock prior refinement requires a V7.2 Stage 2 version bump and a new OSF pre-registration sub-page.

---

## 10. OSF + AsPredicted template fields

| OSF field | Value |
|---|---|
| Title | V7 Clinical Effect-Size Translation Function (MAMMAL Cognitive Enhancement Drug Repurposing) |
| Authors | Pierce Lonergan + Claude Opus 4.7 (1M context) |
| Lock timestamp | TBD |
| Public release | Upon publication |
| Data sources | Roberts 2020, Cochrane systematic reviews, MetaPsy.org, published PET occupancy studies |
| Hypothesis 1 | (See §1) |
| Sample size | V6.A panel = 298 compounds × 22 cognition targets; V7 anchor = 15 compounds |
| Outcomes | Predicted Hedges' *g* per (compound, endpoint) with 95% CrI |
| Analysis plan | (See §2-§7) |
| Pre-registered predictions | P1–P8 (See §3) |
| Falsification thresholds | (See §4) |
| Sensitivity analyses | (See §5) |

---

**This pre-registration is V7.5 Stage 1.** It locks the model and validation gates BEFORE the V7.5 Stage 2 (executed-on-real-data) NUTS run. Any deviation from this pre-registration in the final V7 paper must be flagged as a *post-hoc* analysis.

---

*Generated by `reports/paper-drafts/v7_osf_preregistration.md`. Companion to `design/architecture-and-plans/V4_STATUS_AND_FORWARD_PLAN.md` §13.Y and `design/architecture-and-plans/V6_ARCHITECTURE_PLAN.md` V7 section.*
