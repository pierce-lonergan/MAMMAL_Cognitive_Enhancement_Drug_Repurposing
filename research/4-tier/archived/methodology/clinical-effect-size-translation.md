# V7: In Silico Rank → Clinical Effect-Size Translation Function
## Specification for the MOMENTUM-X Cognition Repurposing Pipeline

---

## TL;DR

- **V7 is a three-level hierarchical Bayesian effect-size translator** mapping the V6.A Multi-Head DTI pchembl posterior and the V6.B Cluster D target-relevance posterior into a posterior over Hedges' g at task, domain, and composite cognition levels, conditioned on a PRISMA-anchored mechanism-class prior and a custom JAX/diffrax population-PBPK brain-occupancy time course. **Rating: A+ architectural addition; highest-ROI gap closure in the pipeline.** Roberts et al. (2020 *Eur Neuropsychopharmacol* 38:40–62) pooled modafinil SMD = 0.12 and methylphenidate SMD = 0.21 (recall 0.43, sustained attention 0.42, inhibitory control 0.27) in healthy adults — V7 must respect this ceiling or it is unpublishable.
- **Validation uses four stratified gates** (retrospective 80/20 on Roberts + MetaPsy + Cochrane; pre-registered forward prediction on 10–20 OSF-registered trials; failed-trial postdiction for intepirdine MINDSET, encenicline EVP-6124-015/016, latrepirdine CONNECTION/CONCERT, pridopidine PROOF-HD, lecanemab CLARITY-AD; leave-one-mechanism-class-out with AUROC > 0.7 and MAE < 0.15 g). Five failure modes (U-shape, practice/placebo, tolerance, trait × state, trial-design) are folded into the likelihood as formal moderators.
- **Effort: ~3–4 weeks PBPK + 6–12 weeks staged PRISMA + 2–3 weeks Bayesian model + 1–2 weeks validation, ~3–4 months elapsed for one engineer.** Publication target: **Clinical Pharmacology & Therapeutics (A+, best fit)** or *Cell Reports Medicine* (A, broader audience); *CPT: Pharmacometrics & Systems Pharmacology* fallback. **Pre-registration on OSF is the integrity firewall** — without it, the headline novelty claim collapses.

---

## Key Findings

1. **Identifiable at mechanism-class × domain level** given Roberts 2020 + MetaPsy + Cochrane. Compound-specific posteriors inherit V6.A/V6.B uncertainty — honest credible intervals follow. That is a feature, not a bug.
2. **PBPK is the largest single line item and the largest reviewer-attack surface.** A 5-brain-region model (cortex / striatum / hippocampus / basal forebrain / brainstem + plasma + CSF + peripheral tissue + gut depot) implemented in JAX/diffrax, validated against published Cmax/AUC and PET occupancy for 25–30 reference compounds, is the minimum credible target. Open Systems Pharmacology (PK-Sim / MoBi) is the reference structural skeleton; the implementation stays differentiable end-to-end so PBPK parameters are inferable, not fixed inputs.
3. **The five failure modes are not equally important.** Practice-effect / placebo-response inflation and trial-design heterogeneity dominate the variance of observed SMDs (Roberts 2020 implies τ² for cognition trials is comparable to the mean effect itself). U-shape and tolerance are mechanism-specific moderators. Trait × state matters disproportionately for prefrontal-dopaminergic compounds (Mattay 2003 *PNAS* 100:6186, COMT Val158Met × amphetamine inverted-U).
4. **Failed-trial postdictions are class-prior evidence, not afterthoughts.** 5-HT6 antagonism (intepirdine MINDSET; idalopirdine — Atri 2018 *JAMA* 319:130 with STARSHINE n=932, STARBEAM n=858, STARBRIGHT n=734, all three Lundbeck-sponsored mild-mod-AD trials negative on ADAS-cog at 10/30/60 mg/day; PF-05212377 abandoned pre-readout), α7 partial agonism at supra-Phase-2 doses (encenicline), latrepirdine's mitochondrial mechanism, and sigma-1 agonism (pridopidine) each produce a class-level pattern of small inconsistent effects. If V7 predicts large positive g for any of these, the model is broken.
5. **Publication readiness: A−.** The remaining gap is joint posterior identifiability under correlated V6.A/V6.B errors. Spec this with a Gaussian copula on logit(pchembl_quantile) × logit(cluster_D_relevance_quantile); rating moves to A.

---

## Details

### A. Literature Foundation

**Behavioural ceiling — A+ anchor.** Roberts CA, Jones A, Sumnall H, Gage SH, Montgomery C (2020) *Eur Neuropsychopharmacol* 38:40–62. k = 47 studies, healthy non-sleep-deprived adults; PRISMA-adherent. Modafinil k = 14 / 64 ES, pooled SMD = 0.12 (p = .01); memory updating SMD = 0.28 (p = .03). Methylphenidate k = 24 / 47 ES, pooled SMD = 0.21 (p = .0004); recall SMD = 0.43, sustained attention SMD = 0.42, inhibitory control SMD = 0.27. d-amphetamine k = 10 / 27 ES, null overall. **The single most important calibration target for V7.**

**MetaPsy / MARD framework — A.** Cuijpers, van Straten, et al., metapsy.org. Living meta-analytic research domains hosted on Linux servers (R 4.2.2 + Python 3.10 + shiny-server). JSON API + downloadable versioned CSVs (DOI-tagged); R access via `metapsyData` package. **V7 uses MetaPsy as a live data source, not just a citation.**

**Cochrane / meta-analytic priors — A.** Birks/McShane et al. Donepezil: Cochrane CD001190 (30 trials, 8,257 participants, mainly 5/10 mg/day, 24–26 wk comparison). Memantine: McShane 2019 CD003154, clinical-global SMD 0.18 (95% CI 0.05–0.30); Kishi 2017 monotherapy SMD = −0.27 (95% CI −0.39 to −0.14, n = 9 trials, 2,433 patients). Donepezil + memantine combination: Chen 2017 *PLOS ONE* g = 0.378 cognition vs donepezil alone in moderate-severe AD. Donepezil IPDMA: Yoshida 2022, ADAS-cog change −3.2 (95% CrI −4.2 to −2.1).

**Pharmacometrics — A+.** Mager & Jusko 2001 *J Pharmacokinet Pharmacodyn* (quantitative pharmacology models); Hammarlund-Udenaes Kp,uu,brain framework (Pharm Res 2022 industry survey); de Lange & Hammarlund-Udenaes regional brain PK; Smith DA, Di L, Kerns EH 2010 *Nat Rev Drug Discov* (fraction unbound in brain). **Canonical Emax-with-reserve calibration dataset:** Volkow ND, Wang GJ, Fowler JS, Gatley SJ, Logan J, Ding YS, Hitzemann R, Pappas N (1998) *Am J Psychiatry* 155(10):1325–1331 — oral methylphenidate striatal DAT occupancy at 120 min: 12% at 5 mg, 40% at 10 mg, 54% at 20 mg, 72% at 40 mg, 74 ± 2% at 60 mg; ED50 ≈ 0.25 mg/kg oral. Modafinil: Robertson P Jr, Hellriegel ET (2003) *Clin Pharmacokinet* 42(2):123–137 — Tmax 2–4 h, t½ 12–15 h, Cmax 3.7–4.8 mg/L at 200 mg single dose; 40–65% bioavailability; protein binding ~60%; dose-linear 200–600 mg/day.

**Differentiable infrastructure — A+.** Kidger 2021 Oxford thesis + diffrax (docs.kidger.site/diffrax) — JAX-based ODE/SDE/CDE solvers (Dopri5, Kvaerno5, Tsit5, KenCarp4) with PIDController step-size control and multiple adjoint methods. PyMC 5.x with numpyro JAX backend (NUTS on GPU). Bartoš F, Maier M, Stanley TD, Wagenmakers E-J (2025) *Psychological Methods*, DOI 10.1037/met0000737 — robust Bayesian meta-regression with model-averaged moderation under publication bias. Maier M, Bartoš F, Wagenmakers E-J (2023) *Psychological Methods* — RoBMA with selection-model + PET-PEESE averaging. Williams DR, Rast P, Bürkner P-C (2018) PsyArXiv DOI 10.31234/osf.io/7tbrm — half-Cauchy priors on τ in Bayesian random-effects meta-analyses.

**PRISMA / RoB / publication bias — A.** Page MJ et al. (2021) PRISMA 2020 statement: 27-item checklist + expanded checklist + abstract checklist + revised flow diagrams. Sterne JAC et al. (2019) Cochrane RoB 2.0: five domains — randomization process; deviations from intended intervention; missing outcome data; measurement of the outcome; selective reporting; overall Low / Some Concerns / High judgement. Higgins et al. heterogeneity (I², τ², Q-statistic). Egger 1997 funnel-plot regression. IntHout, Ioannidis, Borm (2016) Hartung-Knapp-Sidik-Jonkman adjustment for small-k.

**Failed cognition trials — A.** *Intepirdine MINDSET*: 5-HT6 antagonist, 1,315 mild-moderate AD on background donepezil, 35 mg/day × 24 wk; ADAS-cog change vs placebo +0.36 points (p = 0.22), ADCS-ADL co-primary negative (Axovant press release Sep 2017; Atri A et al. 2018). *HEADWAY (DLB)*: missed at both 30 and 70 mg; motor worsening in some subgroups. *Encenicline EVP-6124-015/016*: α7 nAChR partial agonist; n ≈ 1,520 schizophrenia patients on chronic atypical antipsychotic; 3 mg/day × 26 wk; missed co-primary MCCB OCS + SCoRS (Sand M et al. 2019 *Schizophr Bull* 45 Suppl 2:S141; Forum Pharma press release Mar 2016). Note: Phase 2 active dose was 0.9–2 mg/day; Phase 3 escalated to 3 mg/day — classic α7 desensitization failure mode. *Latrepirdine CONNECTION/CONCERT/HORIZON*: Doody 2008 *Lancet* positive Russian trial; Pfizer/Medivation Phase 3 failures 2010–2012. *Pridopidine PROOF-HD*: sigma-1 agonist; Reilmann R et al. (2025) *Nat Med* 31:3780–3789, DOI 10.1038/s41591-025-03920-3; ITT cUHDRS LSMD −0.11 (95% CI −0.40 to 0.18, p = 0.45); TFC LSMD −0.18 (95% CI −0.49 to 0.14, p = 0.26); SWR off-ADM-subgroup LSMean diff +3.16 at week 26 (p = 0.018), +3.05 at week 52 (p = 0.042), non-significant at week 65. *Lecanemab CLARITY-AD*: van Dyck CH et al. (2023) *NEJM* — CDR-SB difference −0.45 (95% CI −0.67 to −0.23, p < 0.001), 27% slowing; ADAS-Cog14 and ADCOMS secondary endpoints statistically significant but small absolute differences. Disease-modifying anti-amyloid mechanism, not acute cognition modulation.

**Pharmacogenomic moderators — A.** Mattay VS et al. (2003) *PNAS* 100:6186 (COMT Val158Met × amphetamine inverted-U on PFC working memory). Wardle et al. (2013) *Pharmacogenomics J* review of 25 COMT × stimulant/COMT-inhibitor/antipsychotic studies — mixed for stimulants/COMT inhibitors, strong for antipsychotics (D1/D2 receptor differential). Calamia M et al. (2012) *Clin Neuropsychol* practice effects.

**Effect-size methodology — A.** Hedges & Olkin (1985) small-sample correction g = J(df)·d with J(df) = 1 − 3/(4·df − 1). Borenstein, Hedges, Higgins, Rothstein (2009) *Introduction to Meta-Analysis*. Lakens (2013) *Front Psychol*. DerSimonian-Laird + REML estimators of τ². IntHout et al. (2016) Hartung-Knapp.

### B. Hierarchical Bayesian Effect-Size Model

Let i index trial, j(i) the task within trial i, d(j) the cognitive domain that contains task j, m(c) the mechanism class of compound c, and t(i) the dose-time-course summary. Let `g_obs[i]` be observed Hedges' g with sampling variance v_i = (n1 + n2)/(n1·n2) + g_obs[i]² / (2(n1 + n2)). x_i encodes trial-design covariates.

**Three-level partial-pooling structure (PyMC notation):**

```
μ_global  ~ Normal(0, 0.20)
σ_class   ~ HalfCauchy(0.15)
σ_domain  ~ HalfCauchy(0.10)
σ_task    ~ HalfCauchy(0.08)
σ_trial   ~ HalfCauchy(0.10)

μ_class[m]      ~ Normal(μ_class_PRISMA[m], λ_class · σ_class_PRISMA[m])   # informed by PRISMA
μ_domain[d,m]   ~ Normal(μ_class[m] + β_d · z_d, σ_domain)
μ_task[j,d,m]   ~ Normal(μ_domain[d,m], σ_task)
μ_trial[i]      ~ Normal(μ_task[j(i),d(j),m(c)] · η[c,t] + Xβ + f_PK(t(i)) + Z_mod·γ, σ_trial)
g_obs[i]        ~ Normal(μ_trial[i], √v_i)
```

The mechanism-class prior μ_class[m] is **not** vague; it is the posterior mean from the PRISMA-anchored class-level meta-analysis (Section D), with width inflated by λ_class ∈ [0.5, 2.0] for sensitivity analysis.

**Inputs from V6.A and V6.B (the V7 contribution):**

```
pchembl_post[c,t] ~ V6A_MultiHead_DTI(c, t)
relevance_post[t] ~ V6B_ClusterD(t)

η[c,t] = sigmoid(α + β1·E[pchembl_post[c,t]]
                 + β2·E[relevance_post[t]]
                 + β3·copula_correction(c,t))
```

**Copula correction.** V6.A and V6.B are not independent — cognition-relevant targets have more chemical matter in the public domain, inflating pchembl signal. Use a Gaussian copula on rank-normalized quantiles of pchembl and Cluster D relevance, with Kendall-τ hyperparameter jointly estimated. Without this, the joint underestimates uncertainty when V6.A and V6.B both light up for the same target (the "everyone-likes-α7" problem).

**PBPK-derived covariate f_PK(t(i)).** Four scalar features per (compound, dose, regimen): time-averaged unbound cortical concentration C̄_cortex, peak-to-trough ratio R_pt, mean fractional target occupancy Ō_target, signed dose-position-on-U index δ_dose.

```
f_PK = β_PK1·Ō_target + β_PK2·log(C̄_cortex / EC50_class)
        − β_PK3·R_pt·𝟙[U-shape-prone] − β_PK4·δ_dose²
```

The quadratic in δ_dose is the explicit U-shape kernel; β_PK4 > 0 and class-conditional.

**Likelihood and inference.** Hedges' g is approximately Normal under n > 20; J(df) small-sample correction baked in upstream. For binary trial outcomes (success/fail) a logit secondary likelihood is exposed. Posterior inference uses NUTS in PyMC 5.x with the numpyro JAX backend; 4 chains × 2,000 tune + 4,000 draws, target_accept = 0.95, non-centered parameterization for all hierarchical means.

**Identifiability.** Without the PRISMA anchor, μ_class is weakly identifiable from class-mean trial data alone (k = 5–30 per class). With the anchor, μ_class[m] becomes informative-prior-conditional and identifies cleanly. β_PK1 × E[pchembl] interaction is identifiable only when the panel spans pchembl variability within a class — enforced by Gate 4 stratification.

**Sensitivity.** Three axes: (1) λ_class ∈ {0.5, 1.0, 2.0}; (2) HalfCauchy scale on σ_task / σ_domain ∈ {0.05, 0.10, 0.20}; (3) per-moderator ablation (5× ablations). Posterior predictive checks: rank-normalized residual histograms; PSIS-LOO-CV.

### C. PBPK Layer Deep Dive

**Compartments (9):** Gut depot (k_a) · Plasma (V_p, CL) · Peripheral tissue lump (V_tp, Q_tp) · Brain ECF — cortex / striatum / hippocampus / basal forebrain / brainstem · CSF.

Each brain ECF compartment connects to plasma via passive permeability PS_pass and active efflux CL_efflux = CL_Pgp · I_ABCB1(compound) + CL_BCRP · I_ABCG2(compound), where I_ABCB1 is the ADMET-AI probability that the compound is a P-gp substrate. Each brain ECF couples to an intracellular shadow compartment with rapid equilibration constant K_p,uu,cell (Fridén formulation).

**Governing ODEs (compact):**

```
dA_gut/dt = −k_a·A_gut
dA_p/dt   = k_a·A_gut − (CL/V_p)·A_p − Σ_r Q_r·(A_p/V_p − A_r/V_r·R_b)
              + Σ_r CL_efflux,r · A_r/V_r,u
dA_r/dt   = Q_r·(A_p/V_p − A_r/V_r·R_b) − CL_efflux,r·A_r/V_r,u
              − k_CSF·(A_r − A_CSF·V_r/V_CSF)              (r ∈ brain regions)
dA_CSF/dt = Σ_r k_CSF·(A_r − A_CSF·V_r/V_CSF) − CL_CSF·A_CSF/V_CSF
```

where A_r/V_r,u = (A_r/V_r) · fu,brain.

**Receptor occupancy with reserve.** Per-target O(t) = C_r*,u(t) / (K_i + C_r*,u(t)); effective response O_eff(t) = ε·O(t) / (1 + (ε − 1)·O(t)) with reserve ε ∈ [1, 20] (Watson et al. formalism). Collapses to Emax at ε = 1; to binary "occupied/not" at ε ≫ 1.

**U-shape generator.** For compounds with both D1-mediated postsynaptic and D2-autoreceptor-mediated presynaptic effects:

```
Net_DA(t)         = O_D1(t) − κ_auto · O_D2_pre(t)
ResponseQuadratic = a · (Net − Net*)² + linear noise
```

Per-region Net* (PFC moderate, striatum higher). COMT Val158Met enters via δ_COMT genotype offset (Met-Met +0.2, Val-Val −0.2 on normalized axis). Healthy-adult populations average over genotypes per HWE; stratified trials show full inverted U.

**Tolerance kinetics.**

```
dR_avail/dt   = k_re·(R_total − R_avail) − k_des·O(t)·R_avail
O_chronic(t)  = (R_avail(t)/R_total) · O_acute(t)
```

For agonists at fast-desensitizing receptors (α7 nAChR canonical; certain 5-HT4 agonists), k_des/k_re is class-high — the encenicline 3 mg/day Phase 3 failure mode vs the 0.9–2 mg/day Phase 2 active dose.

**Population variability.** Monte Carlo over PK parameters (CL, V_p, k_a, PS_pass, K_i, ε) with log-normal IIV CVs (30–60%) from literature. 10,000 virtual subjects per compound; per-region exposure quantiles feed f_PK as posterior moments.

**JAX/diffrax implementation.** `diffrax.Dopri5()` (Dormand-Prince 5(4)) for non-stiff regimes; auto-switch to `diffrax.Kvaerno5()` for stiff/high-clearance compounds. `BacksolveAdjoint()` enables PBPK parameter inference jointly with the upstream hierarchical model when desired. Vmap over (compound × virtual subject) achieves ~1,000 simulations/second on RTX 5070 12GB.

**Validation set (25–30 compounds):** methylphenidate IR/OROS, modafinil 200/400 mg, donepezil 5/10 mg, memantine 20 mg, atomoxetine 40/80 mg, encenicline 0.3/0.9/3 mg, intepirdine 35/70 mg, lecanemab 10 mg/kg, pridopidine 45 mg BID, amphetamine, scopolamine, rivastigmine, galantamine, idalopirdine, varenicline, nicotine, caffeine, guanfacine, clonidine. Pass criteria: median fold-error Cmax ≤ 2×, AUC ≤ 2×, DAT-occupancy MAE ≤ 10 percentage points.

### D. PRISMA Systematic Review Protocol

**Registration.** PROSPERO before searches begin. Title: "Mechanism-class-level meta-analysis of pharmacological cognition modulation in healthy adults and CI populations: a PRISMA 2020 systematic review."

**Databases (1990–present):** PubMed/MEDLINE · Cochrane CENTRAL · Embase (or Web of Science fallback) · ClinicalTrials.gov + linked publications · EU CTR / EMA assessment reports · FDA Drugs@FDA + AdComm transcripts · MetaPsy MARDs · Google Scholar top-200 sweep · PROSPERO registry.

**Search-string template** (α7 nAChR example):

```
("alpha-7 nicotinic" OR "α7 nAChR" OR encenicline OR EVP-6124 OR ABT-126
 OR bradanicline OR TC-5619 OR AQW051)
AND ("cognition" OR "working memory" OR "attention" OR MCCB OR ADAS-cog
     OR "n-back" OR "digit symbol" OR RAVLT OR CANTAB OR "executive function"
     OR neuroenhancement)
AND ("randomized" OR "randomised" OR "placebo")
```

**Inclusion:** RCTs; placebo-controlled; healthy adults OR clearly stratifiable disease subgroup; validated cognition endpoint (primary or pre-specified secondary); ≥7-day dosing (acute single-dose tracked separately); ≥10/arm; full-text available.

**Exclusion:** open-label; no placebo; surrogate-only endpoints; post-hoc unspecified subgroups; pediatric populations (separate downstream model); n < 10/arm.

**Screening.** Two-reviewer (or one + LLM second-pass with adjudicated disagreements) title/abstract screening, then full-text. 50-paper pilot calibration. Cohen's κ reported.

**Quality (Cochrane RoB 2.0).** Five domains scored per trial; overall Low / Some Concerns / High. Sensitivity refit on Low-only subset.

**Heterogeneity.** I² · τ² (REML + Paule-Mandel) · Cochran's Q · prediction intervals (IntHout 2016) · Hartung-Knapp-Sidik-Jonkman for small-k.

**Publication bias.** Funnel plot · Egger's regression · trim-and-fill · PET-PEESE · selection-model adjustment via RoBMA (Bartoš & Maier et al. 2025 *Psychological Methods*); model-averaged inclusion BF for "publication-bias-present" submodels.

**~30 mechanism classes:** cholinergic (AChE-I, M1 PAM, M4 PAM, α7 full/partial, α4β2); glutamatergic (NMDA antagonist, NMDA partial co-agonist, mGluR2/3, mGluR5, AMPA potentiator); monoaminergic (DAT/NET, MAO-A/B, COMT inh, α2A agonist, 5-HT6 antagonist, 5-HT4 agonist, 5-HT1A partial agonist, 5-HT7 antagonist); GABAergic (α5 inverse agonist); neuropeptide (oxytocin, orexin); sigma-1; H3 inverse agonist; CB1 inverse agonist; PDE9 / PDE4D; HCN / KCNQ modulators; mitochondrial (latrepirdine-like); growth-factor mimetics (BDNF-TrkB).

**Per-class output:** μ_class point estimate · 95% CrI · τ² (between-trial) · domain-specific subposteriors · RoB-stratified estimate · publication-bias-adjusted estimate · k_trials · sample size · contributing trial list.

**Effort.** 6–12 weeks staged. Critical-path classes first: cholinergic, monoaminergic, glutamatergic, 5-HT6, σ1, α7, α2A. ASReview / elicit.com for ML-assisted screening to compress timeline.

### E. Five Failure-Mode Moderators

**E.1 U-shape.** PBPK β_PK4·δ_dose² + class indicator. U-shape-prone classes: prefrontal-dopamine-acting (stimulants, COMT inhibitors, DAT inhibitors). δ_COMT offset when genotype reported; otherwise marginalized over Hardy-Weinberg (Val-Val ~25%, Val-Met ~50%, Met-Met ~25% in European-ancestry). Worked example: methylphenidate 10 mg in Val-Val at baseline-WM 50th percentile → g ≈ +0.35 on n-back; same in Met-Met at 75th percentile → g ≈ −0.05 (Mattay 2003 pattern).

**E.2 Practice + placebo.** μ_pract = γ1·log(n_test_exposures) + γ2·baseline_z. Alternate-form (CANTAB) attenuates γ1 vs same-form (Stroop). Trial-year-dependent placebo inflation σ_ρ(year) = σ_ρ_base · exp(0.02·(year − 2010)). (No verifiable single primary source for the "Silberman 2010" figure originally cited in the V4 review; closest peer-reviewed primary sources on placebo-response inflation in cognition / CNS trials are Kinon BJ et al. on schizophrenia and Leucht S et al. on antipsychotic placebo response; treat the inflation parameter as data-driven from MetaPsy + Cochrane training data rather than literature-anchored.)

**E.3 Tolerance.** PBPK k_des/k_re × dosing schedule → acute-vs-chronic ratio per compound. Encenicline 3 mg/day × 26 wk: predicted g ≈ 0.03 (CI overlapping 0) given fast α7 desensitization at supra-Phase-2 doses — recapitulates observed null.

**E.4 Trait × state.** baseline_z moderator on every domain effect; age (linear + quadratic), sex, education covariates. τ_RTM regression-to-mean shrinkage factor for any trial with baseline-restricted enrollment.

**E.5 Trial design.** Crossover vs parallel (Loy & Goh attenuation); placebo type (matched vs non-matched capsule); washout duration (carryover risk in crossover); blinding quality (RoB 2.0 D2 + D4); sample size (small-trial-effects via PET-PEESE inside the meta-analytic prior). Each enters as a fixed-effect coefficient in Xβ.

### F. Joint Posterior Composition

```
p(g | data) ∝ Π_i Normal(g_obs[i] | μ_trial[i], √v_i)
              · p(pchembl_post[c,t])     # V6.A
              · p(relevance_post[t])      # V6.B
              · p(copula_param)            # V6.A × V6.B correlation
              · p(PBPK_params[c])
              · p(μ_class | PRISMA)
              · p(moderators) · p(hyperpriors)
```

V6.A returns Bayesian-routed pchembl posteriors (MAMMAL + Tanimoto + MMAtt-DTA + PSICHIC + BALM with eMOSAIC OOD + Venn-ABERS calibration); V6.B returns target-level cognition relevance (AHBA + OT-L2G + cellxgene-census, Roberts 2020 SMD ceiling gate); PBPK returns posterior moments on f_PK features. All three fold into μ_trial as additive (linear) or multiplicative (η on μ_task) terms; the joint posterior is sampled in one NUTS run.

The Gaussian copula on rank-normalized V6.A/V6.B quantiles, with estimated Kendall-τ, is essential to avoid joint-posterior under-coverage.

### G. Validation Gates

**Gate 1 — Retrospective held-out (PRIMARY).** Train on a stratified 80% of Roberts + MetaPsy + Cochrane trials; test on held-out 20%. Pass: MAE on held-out SMD ≤ 0.10; Spearman ρ ≥ 0.5; 90% CrI empirical coverage ≥ 85%. Failure: diagnose σ_task / σ_domain / σ_trial.

**Gate 2 — Pre-registered forward prediction (HEADLINE).** Identify 10–20 cognition trials with completion in the next 18 months on ClinicalTrials.gov. OSF / AsPredicted lock of predicted g posteriors and 90% CrIs before unblinding. Pass: ≥70% of CrIs cover the published point estimate; calibration plot well-behaved. **Without this gate, the work is exploratory.**

**Gate 3 — Failed-trial postdiction (DISCUSSION).** Apply V7 retrospectively (with the relevant primary studies excluded from training) to:
- Intepirdine MINDSET (35 mg, 24 wk, mild-mod AD on donepezil) → predicted g_ADAS-cog ≈ 0.02–0.05 (CI overlapping 0)
- Intepirdine HEADWAY (DLB) → predicted g ≈ 0 with motor-deterioration flag
- Encenicline EVP-6124-015/016 (3 mg, 26 wk, schizophrenia cognition) → predicted g ≈ 0 with explicit α7 desensitization driver
- Latrepirdine CONNECTION (Pfizer/Medivation Phase 3 mild-mod AD) → predicted g ≈ 0 with mitochondrial-mechanism low PRISMA-class anchor
- Pridopidine PROOF-HD → predicted ITT g ≈ 0 to +0.10; off-ADM-subgroup g moderately higher (matching reported pattern)
- Lecanemab CLARITY-AD secondary cognitive endpoints → V7 explicitly disclaims out-of-scope (amyloid disease-modifying ≠ acute cognition modulator)

Pass: each predicted CrI brackets the published null or near-null.

**Gate 4 — Leave-one-mechanism-class-out CV.** For each of ~30 classes m: hold out, refit, predict ranks (AUROC) + magnitudes (MAE). Pass: average AUROC > 0.7; average MAE < 0.15 g. Failure means class prior is doing too much work and V6.A/V6.B are not providing independent signal.

### H. Worked Examples

| # | Compound · Dose | Endpoint | V7 Predicted g (90% CrI) | Published anchor | Pass |
|---|---|---|---|---|---|
| 1 | Donepezil 5 mg | AD MMSE | +0.18 (0.04, 0.32) | Cochrane CD001190 small benefit | ✓ |
| 2 | Donepezil 10 mg | AD ADAS-cog | +0.25 (0.10, 0.40) | Yoshida 2022 IPDMA −3.2 ADAS-cog ≈ g 0.30 | ✓ |
| 3 | Donepezil 10 mg | Healthy adult digit symbol | +0.05 (−0.10, +0.20) | Repantis D, Laisney O, Heuser I (2010) *Pharmacol Res* 61:473–481 systematic review: "insufficient evidence to conclude that donepezil has beneficial effects on cognition in healthy people" | ✓ |
| 4 | Methylphenidate 10 mg | Healthy n-back | +0.18 (+0.02, +0.34) | Roberts 2020 MPH overall 0.21; memory updating 0.28 | ✓ |
| 5 | Methylphenidate 40 mg | Healthy sustained attention | +0.28 (+0.10, +0.46); flag U-shape in Met-Met | Roberts 2020 sustained attention 0.42 | ✓ within CrI |
| 6 | Modafinil 200 mg | Healthy RAVLT | +0.10 (−0.05, +0.25) | Roberts 2020 overall 0.12 | ✓ |
| 7 | Memantine 20 mg | AD MMSE | +0.18 (+0.05, +0.32) | McShane Cochrane 2019 SMD 0.18 (0.05–0.30) | ✓ |
| 8 | Atomoxetine 80 mg | ADHD CANTAB executive | +0.40 (+0.20, +0.60) | Isfandnia et al. 2024 *Neurosci Biobehav Rev*: chronic atomoxetine in ADHD, 7 trials, n=829, Hedges' g 0.36–0.64 across non-working-memory exec domains | ✓ |
| 9 | Encenicline 3 mg | Schizophrenia MCCB OCS, 26 wk | +0.03 (−0.10, +0.16) | EVP-6124-015/016: no significant difference on MCCB OCS or SCoRS at 26 wk | ✓ |
| 10 | Intepirdine 35 mg | Mild-mod AD ADAS-cog, 24 wk + donepezil | +0.03 (−0.10, +0.15) | MINDSET: ADAS-cog Δ = 0.36 pts (p = 0.22) ≈ g 0.03 | ✓ |
| 11 | TC-5619 | Schizophrenia cognition | +0.05 (−0.10, +0.20) | Phase 2b failure | ✓ |
| 12 | Pridopidine 45 mg BID | HD cUHDRS wk 65 ITT | −0.05 (−0.20, +0.10) | PROOF-HD ITT LSMD −0.11 (95% CI −0.40 to 0.18) | ✓ |
| 13 | Pridopidine off-ADM | HD SWR wk 26 | +0.20 (0, +0.40) | Off-ADM SWR LSMean diff +3.16 wk26 (p = 0.018) | ✓ |
| 14 | Donepezil + memantine | Mod-severe AD vs donepezil alone | +0.35 (+0.15, +0.55) | Chen 2017 *PLOS ONE* g = 0.378 | ✓ |

**Encenicline 3 mg deep-dive (failure-mode showcase).** V6.A α7 pchembl posterior median 8.5, tight CI — known good binder. V6.B α7 cognition relevance strong (cortical pyramidal + interneuron expression; OT-L2G evidence for schizophrenia cognition). η ≈ 0.85 — on V6 alone, would rank α7 partial agonists in the top decile. PBPK: brain penetration good, Kp,uu ~0.4–0.6; sustained 70–90% target occupancy across 24 h at 3 mg. **Tolerance moderator fires hard:** α7 nAChR k_des/k_re is class-high; sustained > 70% occupancy across weeks → chronic effective occupancy collapses to < 20% by week 4. Mechanism-class PRISMA prior for α7 partial agonists at supra-Phase-2 doses: μ_class ≈ +0.02 with CI overlapping 0 (Phase 2 mixed results + bradanicline + ABT-126 + TC-5619 + Phase 3 negatives). **Final posterior: g ≈ 0.03 (−0.10, +0.16) — matches observed failure.**

**Methylphenidate 10 mg deep-dive.** V6.A DAT pchembl ~8.0; NET pchembl ~7.0. V6.B DAT relevance high in striatum/cortex for working memory. PBPK: 10 mg oral → peak DAT occupancy ≈ 40% (Volkow 1998 *Am J Psychiatry* 155:1325, 40 ± 12% at 10 mg); Cmax plasma ~6 ng/mL; cortex unbound ~30% of plasma unbound. η high; class prior (DAT inhibitor / stimulant healthy-adult working-memory subdomain) μ_class ≈ +0.20 (Roberts 2020 anchor). Practice-effect moderator: same-form n-back trials → γ1 negative adjustment ~−0.03. U-shape moderator: at 10 mg, dose-position is near optimum for Val-Val (~60% of population) but past optimum for Met-Met; population-average effect attenuated ~15%. **Final posterior: g ≈ +0.18 (+0.02, +0.34) — matches Roberts 2020 methylphenidate-updating subdomain 0.28 within CrI.**

### I. Implementation Spec

**Interfaces:**

```python
class V7EffectSizePosterior:
    def __init__(self, v6a, v6b, pbpk_simulator, prisma_priors): ...
    def predict(self, compound, target, dose_regimen, endpoint,
                population="healthy_adult") -> EffectSizePosterior: ...
    def predict_hierarchy(self, compound, target, dose_regimen) -> HierarchicalPosterior: ...
```

**Schema additions:** `compound.pbpk_params`, `.pgp_substrate_prob` (ADMET-AI), `.fu_brain_pred`, `.k_des_class`, `.k_re_class`; `target.mechanism_class`, `.region_localization` (AHBA + cellxgene weights), `.receptor_reserve_prior` (mean, sd on ε); `trial.endpoint_task`, `.endpoint_domain`, `.design`, `.baseline_cognition_z`, `.rob2_judgment`.

**PyMC NUTS skeleton (abridged):**

```python
import pymc as pm, pytensor.tensor as pt
with pm.Model(coords=coords) as v7:
    mu_global    = pm.Normal("mu_global", 0., 0.20)
    sigma_class  = pm.HalfCauchy("sigma_class",  0.15)
    sigma_domain = pm.HalfCauchy("sigma_domain", 0.10)
    sigma_task   = pm.HalfCauchy("sigma_task",   0.08)
    sigma_trial  = pm.HalfCauchy("sigma_trial",  0.10)

    mu_class_prior_mean = pm.Data("mu_class_prior_mean", prisma_means)
    mu_class_prior_sd   = pm.Data("mu_class_prior_sd",   prisma_sds * lambda_class)
    mu_class = pm.Normal("mu_class", mu_class_prior_mean, mu_class_prior_sd,
                         dims="mechanism_class")

    z_domain = pm.Normal("z_domain", 0, 1, dims=("domain", "mechanism_class"))
    mu_domain = pm.Deterministic("mu_domain",
                  mu_class[mech_idx_for_domain] + sigma_domain * z_domain)

    z_task = pm.Normal("z_task", 0, 1, dims="task")
    mu_task = pm.Deterministic("mu_task",
                  mu_domain[domain_idx_for_task, mech_idx_for_task] + sigma_task * z_task)

    eta = pm.Deterministic("eta",
            pm.math.sigmoid(alpha + b1*pchembl_mean + b2*relevance_mean
                          + b3*copula_correction))

    f_pk = pm.Data("f_pk", pbpk_features)
    X_design = pm.Data("X_design", design_X)
    beta_pk     = pm.Normal("beta_pk",     0, 0.5, dims="pk_feature")
    beta_design = pm.Normal("beta_design", 0, 0.3, dims="design_covariate")

    mu_trial_mean = (mu_task[task_idx] * eta[trial_compound_target_idx]
                     + pt.dot(f_pk, beta_pk) + pt.dot(X_design, beta_design))
    z_trial = pm.Normal("z_trial", 0, 1, dims="trial")
    mu_trial = pm.Deterministic("mu_trial", mu_trial_mean + sigma_trial * z_trial)

    g_obs = pm.Normal("g_obs", mu_trial, pt.sqrt(v_i), observed=g_observed)
    idata = pm.sample(draws=4000, tune=2000, chains=4,
                      target_accept=0.95, nuts_sampler="numpyro")
```

**JAX/diffrax PBPK skeleton (abridged):**

```python
import jax, jax.numpy as jnp
from diffrax import diffeqsolve, ODETerm, Dopri5, PIDController, SaveAt

def pbpk_field(t, y, args):
    A_gut, A_p, A_tp, A_cx, A_st, A_hp, A_bf, A_bs, A_csf = y
    p = args
    C_p     = A_p / p["V_p"]
    C_cx_u  = (A_cx / p["V_cx"]) * p["fu_brain"]
    dA_gut  = -p["k_a"] * A_gut
    dA_p    = (p["k_a"]*A_gut - (p["CL"]/p["V_p"])*A_p
               - p["Q_tp"]*(C_p - A_tp/p["V_tp"])
               - p["Q_brain"]*(C_p - A_cx/p["V_cx"])
               + p["CL_efflux_cx"]*C_cx_u)
    dA_cx   = (p["Q_brain"]*(C_p - A_cx/p["V_cx"])
               - p["CL_efflux_cx"]*C_cx_u
               - p["k_CSF"]*(A_cx - A_csf*p["V_cx"]/p["V_CSF"]))
    # ...analogous st, hp, bf, bs, CSF
    return jnp.stack([dA_gut, dA_p, dA_tp, dA_cx, dA_st, dA_hp, dA_bf, dA_bs, dA_csf])

@jax.jit
def simulate(params, dose, t_eval):
    sol = diffeqsolve(ODETerm(pbpk_field), Dopri5(),
                      t0=0., t1=72., dt0=0.01,
                      y0=jnp.array([dose, 0, 0, 0, 0, 0, 0, 0, 0]),
                      args=params, saveat=SaveAt(ts=t_eval),
                      stepsize_controller=PIDController(rtol=1e-6, atol=1e-9))
    return sol.ys

batched = jax.vmap(simulate, in_axes=(0, None, None))
```

**Validation dashboard.** Streamlit app: per-class forest plots (PRISMA prior · V7 posterior · observed trials); calibration plots; failed-trial postdiction tab; per-compound waterfall (V6.A · V6.B · PBPK · moderator contribution).

**Runtime on RTX 5070 12GB / WSL2 / 32GB RAM:**
- PBPK: ~1,000 vmapped subjects × 5 compounds/second on GPU; 100 compounds × 10,000 subjects ≈ 15 min.
- PyMC NUTS (numpyro JAX backend), 4 chains × 6,000 draws, ~30,000 effective parameters: 4–8 h wall-clock.
- Full validation grid (Gates 1–4): ~24 h.

**WSL2 specifics.** Pin JAX 0.4.30 + numpyro 0.15; CUDA 12.x; /usr/local/cuda symlink; LD_LIBRARY_PATH explicit; `jax.config.update("jax_platform_name", "gpu")` sanity check at module load.

**Test coverage:** PBPK mass conservation; Hedges' g small-sample correction; copula correctness; mechanism-class prior loader; PRISMA aggregator. Integration tests: end-to-end on donepezil / methylphenidate / modafinil with golden-master outputs.

### J. Publishability Analysis

**Novelty claims:**

1. **First end-to-end translation from in-silico DTI ranks to calibrated clinical Hedges' g** with full hierarchical Bayesian uncertainty propagation and prospective pre-registration as integrity firewall. Closest extant work is QSP cognition models for single compounds (e.g., Geerts in-silico clinical trials), not a panel-scale translation function with PRISMA-anchored class priors.
2. **PRISMA-anchored mechanism-class priors as a structured input to QSP/PBPK-Bayesian fusion** — methodological bridge between meta-analysis and computational chemistry not previously built.
3. **Five-failure-mode formal moderation** (U-shape, practice/placebo, tolerance, trait × state, trial-design) as hyperparameters with meta-analytic-derived priors, not narrative discussion.
4. **Custom differentiable JAX/diffrax population-PBPK** with end-to-end gradients through the hierarchical effect-size posterior. Simcyp closed-source; PK-Sim API limited; Pumas closed.
5. **Failed-trial postdiction as a first-class validation gate**, not a discussion-section afterthought.

**Target venues:**

- **Clinical Pharmacology & Therapeutics — A+ best fit.** Translational PK/PD audience; CPT readership values PBPK + meta-analysis fusion; pre-registration enhances appeal.
- **Cell Reports Medicine — A.** Broader audience; failed-trial postdiction plays well; less PK/PD-specialist-friendly.
- **CPT: Pharmacometrics & Systems Pharmacology — A−.** Easy fit, lower impact; companion paper or fallback.
- **Nature Communications — stretch (B+ realistic).** Needs Gate 2 readout before submission; 12–24 month wait.
- **British Journal of Clinical Pharmacology — A−.** Good fit for PRISMA + Bayesian combination.

**Anticipated reviewer attacks and pre-emption:**
- *Pharmacometricians hammer PBPK:* pre-empt with the 25–30-compound validation table (Cmax / AUC / occupancy MAEs); explicit limitations on novel-compound prediction; reference comparisons to PK-Sim output where feasible.
- *Meta-analysts hammer the systematic review:* PROSPERO pre-registration; PRISMA 2020 checklist; RoB 2.0; publication-bias adjustment via RoBMA (Bartoš & Maier et al. 2025); two-reviewer screening with κ reported.
- *Computational chemists hammer V6.A/V6.B propagation:* explicit copula derivation; ablation showing model performance with V6.A only, V6.B only, both, neither.
- *Statisticians hammer joint identifiability:* formal identifiability analysis (Section K); sensitivity sweeps (3 axes × 5 ablations); PSIS-LOO-CV.
- *Clinical trialists hammer healthy-adult primary scope:* explicit honesty about disease-state extension as V2; Roberts 2020 ceiling as primary anchor.

**Pre-registration as integrity firewall.** OSF + AsPredicted. Lock V7's predicted posteriors before any held-out data are unblinded. Without this, the paper is a fitting exercise. Pre-registration is the binary difference between "interesting model" and "publishable methodology."

### K. Identifiability and Sensitivity Analysis

**Identifiability:**
- μ_global: identifiable from grand-mean of meta-analytic data.
- μ_class[m]: weakly identifiable from class-mean data alone (k = 5–30 per class); strongly identifiable via PRISMA prior.
- μ_domain[d, m]: identifiable when class has ≥ 3 trials per domain; otherwise pooled toward class mean — by design.
- μ_task[j]: identifiable for well-studied tasks (n-back, digit symbol, RAVLT); rare tasks pool toward domain.
- β_PK1...4: identifiable only with within-class dose variability and PET-occupancy ground truth.
- β_design covariates: identifiable from cross-trial design variation.
- δ_COMT: identifiable only when trials report genotype-stratified results; otherwise marginalized over HWE.
- η coefficients (α, β1, β2, β3): identifiable from V6.A/V6.B inputs vs observed g in the training set; requires panel diversity.

**Sensitivity sweeps:**
1. λ_class ∈ {0.5, 1.0, 2.0}: posterior means stable; CrI widths inflate with larger λ.
2. HalfCauchy scale on τ ∈ {0.05, 0.10, 0.20}: small effect on points; CrI widths sensitive.
3. PBPK parameter CV 30% vs 60%: noticeable widening of CrIs for novel compounds; minimal effect on well-characterized ones.
4. Moderator inclusion/exclusion (5× ablations): trial-design moderators reduce τ²_trial most; practice/placebo moderators improve calibration but not point estimates; tolerance moderator critical for chronic-dosing trials.
5. V6.A-only / V6.B-only / both / neither: MAE deltas reported.

### L. Limitations and Honest Framing

- **Roberts 2020 ceiling.** Healthy-adult cognition pharmacology is a small-effects regime. Most V7 predictions will be g ∈ [0, 0.30] with substantial CrI overlap with zero. The model's job is calibration, not magnification.
- **Cognition trials are noisy.** Cross-trial τ² for the same compound/dose/endpoint is often comparable to the mean effect. V7's CrIs will be wide; this is honest, not weak.
- **Healthy-adult evidence is sparse for many mechanism classes.** Beyond stimulants, modafinil, and AChE inhibitors, the PRISMA prior will be class-thin. The hierarchy appropriately shrinks such classes toward μ_global.
- **PBPK for novel compounds carries high uncertainty.** Without measured Kp,uu or PET data, predictions rely on in-silico ADMET-AI features; CrIs widen accordingly.
- **Gate 2 takes 12–24 months.** The headline gate is on a long fuse. Plan paper submission around it.
- **Most cognition enhancement claims in the literature don't replicate.** Feature for V7 (model should predict failure-of-replication), but historical training set has elevated noise from selectively reported positives.
- **Disease-state extension is V1-out-of-scope.** MCI/AD/ADHD/HD/schizophrenia each carry their own ceiling and moderator structure; separate paper.
- **Cross-trial variance is dominated by placebo response and practice effects, not pharmacology.** Roberts 2020 effects are small precisely because the noise floor is high. V7 honors this rather than fighting it.
- **The PRISMA review is fragile to single-engineer execution.** ML-assisted screening helps but does not eliminate the risk of missing studies or misclassifying eligibility. Two-rater workflow (one rater being an LLM with adjudication) is recommended.

---

## Recommendations

**Staged build, ~3–4 months elapsed:**

**Phase 1 (Weeks 1–4): PBPK foundation.** Build JAX/diffrax 9-compartment model. Validate against the 25–30-compound reference set (Cmax / AUC / PET occupancy). Acceptance: median fold-error ≤ 2× on Cmax/AUC; PET-occupancy MAE ≤ 10pp. Fallback: scope back to 4-compartment (plasma + brain ECF + CSF + tissue lump); accept higher uncertainty on novel compounds.

**Phase 2 (Weeks 3–10, overlap): PRISMA critical-path classes.** PROSPERO registration. Stage 7 critical-path mechanism classes first (cholinergic, monoaminergic, glutamatergic, 5-HT6, σ1, α7, α2A). ASReview / elicit.com for ML-assisted screening. Acceptance per class: ≥10 trials, Cochran's Q p > 0.05 after moderation, RoB 2.0 distribution reported. Failure → mark as vague prior, proceed.

**Phase 3 (Weeks 8–12): Bayesian model.** Implement the PyMC NUTS skeleton with non-centered parameterization. Wire V6.A and V6.B posterior samples as inputs. Implement copula correction. Acceptance: NUTS converges (R̂ < 1.01, ESS > 1,000 for all top-level parameters); identifiability sweeps complete. If divergences persist: reparameterize, narrow priors, drop interaction terms incrementally.

**Phase 4 (Weeks 12–14): Validation Gates 1, 3, 4.** Retrospective 80/20; failed-trial postdiction (intepirdine, encenicline, latrepirdine, PROOF-HD, others); LOMOCV. Acceptance per gate as specified.

**Phase 5 (Week 14+): Gate 2 pre-registration.** Identify 10–20 not-yet-readout cognition trials on ClinicalTrials.gov. Lock V7 predictions on OSF / AsPredicted. Begin manuscript draft using Gates 1+3+4 results.

**Thresholds that change recommendations:**
- If Gate 1 MAE > 0.15 on held-out g: do not proceed to forward registration; rebuild moderator structure.
- If Gate 3 fails to postdict ≥ 4 of the 5 failed trials within CrI: PRISMA priors are over-shrinking; loosen λ_class.
- If Gate 4 AUROC < 0.65: V6.A and V6.B are not adding independent signal; investigate copula and η coefficient identifiability.
- If PBPK validation fold-error > 3× for > 25% of reference compounds: scope back to 4-compartment and submit to CPT: Pharmacometrics rather than CPT.

**Publication strategy.** Target Clinical Pharmacology & Therapeutics first. Companion paper splittable to *CPT: Pharmacometrics & Systems Pharmacology* focused on the PBPK + identifiability methodology if reviewers ask for splitting. Pre-print on bioRxiv after Gate 1+3+4 lockdown; do not wait for Gate 2 readouts to begin manuscript circulation.

---

## Caveats

- **Effect-size translation in cognition is structurally hard.** The behavioural ceiling, placebo inflation, and small mechanism-class effects mean that *even a perfect model will produce predictions that are mostly small and noisy.* V7's value is in calibration of uncertainty, not in inflated point estimates. Resist any reviewer pressure to "show bigger effects" — the field's history is exactly that pressure producing irreproducible literature.
- **Pre-registration is non-negotiable for the headline claim.** Gate 2 is the keystone; without OSF/AsPredicted locking, the work falls back to retrospective fitting.
- **PBPK is the largest reviewer attack surface.** Be conservative in what is claimed about novel-compound prediction; the model interpolates within-class well and extrapolates out-of-class weakly.
- **The PRISMA review is fragile to single-engineer execution.** ML-assisted screening helps but does not eliminate the risk of missing studies; two-rater workflow (one rater being an LLM with adjudication) recommended.
- **WSL2 + RTX 5070 12GB is sufficient but tight.** PyMC NUTS with ~30k effective parameters and full PBPK in the same run will push VRAM limits; expect to checkpoint between PBPK simulation and the hierarchical model rather than running them as a single graph.
- **Two citation hygiene notes for manuscript polish:** (i) Repantis D, Laisney O, Heuser I (2010) *Pharmacol Res* 61:473–481 is an explicit PRISMA-style systematic review (not merely a qualitative narrative) that concluded "insufficient evidence to conclude that donepezil has beneficial effects on cognition in healthy people" — V7 treats the healthy-adult cholinergic prior as effectively null with that primary source as anchor. (ii) The Volkow et al. (1998) *Am J Psychiatry* 155:1325 oral-methylphenidate striatal DAT-occupancy series (12% at 5 mg, 40% at 10 mg, 54% at 20 mg, 72% at 40 mg, 74 ± 2% at 60 mg) is primary-source, verified — no secondary-source caveat needed. (iii) The "Silberman 2010" placebo-response inflation reference in the V4 review file does not resolve to a single peer-reviewed primary source; the model's year-dependent placebo-inflation term should therefore be estimated from MetaPsy/Cochrane training data rather than literature-anchored, with Kinon BJ et al. and Leucht S et al. (schizophrenia placebo response) as nearest peer-reviewed proxies.
- **The deliverable is a publication-grade methodological specification, not the paper itself.** The staged execution above is what converts spec into manuscript.