# In Silico Rank → Clinical Effect-Size Translation Function: A Methodology Pre-Registration for Bayesian Cognition-Enhancement Drug Repurposing

## TL;DR
- The V4–V6 pipeline predicts target-level ranks but cannot yet emit per-(compound, dose, endpoint) predicted Hedges' g with credible intervals; this gap separates a wet-lab shortlist from a clinical-trial prior, and no published methodology currently performs this translation with calibration error (the closest analogues are HINT binary trial-outcome classifiers, Wong/Siah/Lo ML approval-probability models, and QSP mechanistic PK/PD — none produce SMD outputs).
- We propose a hierarchical Bayesian model `g_hat ~ N(α_class + f(occupancy) · β_target + γ_endpoint, σ²)` anchored on Roberts 2020 (k=47 trials; modafinil pooled SMD=0.12, MPH SMD=0.21, d-amph null) and Birks 2018 Cochrane (donepezil 10 mg ADAS-Cog MD=−2.67 [−3.31, −2.02], back-calculable to SMD≈−0.30), with occupancy mapped from V6.A pchembl posteriors via Hill curves calibrated against PET-validated anchors (donepezil 19.1% cortical AChE inhibition → cognitive change in Bohnen 2005; haloperidol 60–80% D2 occupancy → response in Kapur 2000), gated by V6.B θ̄ posteriors, and validated by leave-one-trial-out coverage on the MetaPsy + Cochrane corpus.
- The 24-week build plan ships in three Phases (D.1 prior elicitation → D.2 PK/PD assembly → D.3 hierarchical model fit → D.4 LOO validation), with pre-committed falsification thresholds (calibration slope ∈ [0.8, 1.2]; 90% CI coverage ≥ 85%; must beat the mechanism-class-mean baseline by ΔELPD > 4 SE), and a negative-result publication path to CPT:PSP if the translation function fails. Cell Reports Medicine is the target venue if validation succeeds; CPT:PSP is the fallback regardless.

## Key Findings

**(F1) The methodological gap is real and quantifiable.** Across a focused search of in silico-to-clinical translation pipelines, no published methodology emits a continuous predicted Hedges' g with credible intervals from upstream binding/affinity predictions. Fu et al. 2022 (HINT, *Patterns* 3(4):100445; DOI 10.1016/j.patter.2022.100445) predict binary phase-success with F1 = 0.665/0.620/0.847 for phases I/II/III, but never effect magnitude. Wong, Siah & Lo (*Biostatistics* 20(2):273–286, 2019; DOI 10.1093/biostatistics/kxx069) provide approval-probability AUCs of 0.78 (Ph2→App) and 0.81 (Ph3→App), again binary. Quantitative systems pharmacology (Geerts et al. *Schizophr Bull* 44:S221, 2018; In Silico Biosciences cortico-striatal-thalamic loop) translates 32-target binding profiles into predicted EPS phenotypes blinded against 1,124 South London/Maudsley CRIS patients but does not output Hedges' g. The proposed translation function therefore occupies an open methodological niche.

**(F2) CNS attrition makes the translation function high-value.** BIO/QLS/Informa 2021 (12,728 phase transitions, 2011–2020) report Phase I→II 47.7%, II→III 26.8%, III→NDA 53.1%, NDA→Approval 86.7% for neurology (compounded LOA from Phase I = ~5.9%, well below all-indications 7.9%). Psychiatry LOA = 7.3%. Cummings et al. 2014 *Alz Res Ther* 6(4):37 (244 AD compounds, 413 trials, 2002–2012) reports a 99.6% failure rate for AD drug development — one approval (memantine, 2003) in a decade. The Phase II→III bottleneck (26.8%) is the dominant choke point; a calibrated SMD-predictor is precisely the tool that could redirect Phase II go/no-go decisions.

**(F3) Roberts 2020 anchors the healthy-adult ceiling but does not span mechanism space.** Roberts CA et al. *Eur Neuropsychopharmacol* 38:40–62 (2020); k=47 studies. Modafinil overall SMD = 0.12 (p=0.01); memory updating SMD = 0.28 (p=0.03). MPH overall SMD = 0.21 (p=0.0004); recall 0.43, sustained attention 0.42, inhibitory control 0.27. D-amphetamine null. This frames Cluster D's hard ceiling at SMD ≈ 0.5 in best sub-domains. For other mechanism classes the field offers a heterogeneous library of pooled estimates: caffeine accuracy g = 0.27 and reaction-time g = 0.28 in 31 RCTs / n=1,455 rested healthy adults (Kløve & Petersen 2025 *Psychopharmacology* 242(9):1909–1930; DOI 10.1007/s00213-025-06775-1); nicotine attention/short-term memory small-to-moderate (Heishman 2010 *Psychopharmacology* 210:453–69, 41 RCTs); vortioxetine DSST SMD = 0.325 (95% CI 0.120–0.529) in MDD (Baune 2018 *Int J Neuropsychopharmacol* 21:97–107).

**(F4) Birks 2018 donepezil and memantine MAs report mean differences, not SMDs, requiring back-calculation.** Donepezil 10 mg (Birks & Harvey 2018, Cochrane CD001190.pub3): ADAS-Cog MD = −2.67 [−3.31, −2.02], 1,130 participants; MMSE MD = +1.05 [0.73, 1.37], 1,757 participants; SIB MD = +5.92 [4.53, 7.31]. Memantine (Winblad/McShane Cochrane): ADAS-Cog SMD = −0.21 [−0.34, −0.08], 6 studies; moderate-severe AD SMD = −0.29 [−0.54, −0.03]. Back-calculation against typical ADAS-Cog SD ≈ 8–10 in mild-moderate AD yields donepezil SMD ≈ −0.27 to −0.33. These provide MCI/AD anchors complementary to the Roberts 2020 healthy-adult anchors. Note that in any-stage dementia (18 RCTs, n=5,948), Sheikh & Ammar 2024 (*Front Neurosci* 18:1398952) report donepezil 10 mg vs placebo MMSE g = 2.27 [1.25, 3.29] with no significant ADAS-Cog effect — the apparent inconsistency reflects scale-direction conventions and dementia severity heterogeneity, and should be treated cautiously as a class outlier.

**(F5) PET-validated occupancy→effect anchors exist for at least three target families.** Donepezil: Bohnen et al. 2005 *J Neurol Neurosurg Psychiatry* 76:315–319 reports mean cortical AChE inhibition of 19.1% (SD 9.4%) after 12-week donepezil 5 mg, correlating with Stroop interference improvement (R² = 0.59, p<0.01), but not with primary memory tests; donepezil also achieves ~60% (5 mg) and ~75% (10 mg) sigma-1 occupancy (Ishikawa et al. 2009 PMID 19573265). Haloperidol D2: Kapur et al. 2000 *Am J Psychiatry* 157:514–520 shows clinical response plateau begins at ~65% striatal D2 occupancy, with EPS rising above ~78%. These furnish the empirical priors for the Hill function f(occupancy → effect) that the translation function requires.

**(F6) Dopaminergic effects are inverted-U, not sigmoidal — and the model must encode this.** Cools & D'Esposito 2011 *Biol Psychiatry* 69:e113–e125 establishes inverted-U dependence of cognition on PFC dopamine, with effect direction conditional on baseline DA tone (responders vs non-responders by COMT genotype). Narayanan and colleagues (2022, *Behavioral Neuroscience* 136(3):207; DOI 10.1037/bne0000512; PMC9364670) quantified this with a 75-study meta-analysis across rodents, non-human primates, and humans: "10% of the variance in working memory behavior was explained by manipulations of prefrontal dopamine, and 26% of the variance was explained by prefrontal D1DR manipulations." The translation function therefore needs a class-conditional response surface: sigmoidal/Hill for cholinergic, glutamatergic, monoaminergic; inverted-U with baseline covariate for dopaminergic (and arguably noradrenergic — Arnsten guanfacine literature).

**(F7) ML-predicted BBB scores have known calibration limits; Kp,uu is the right currency.** Kim et al. 2022 *Front Pharmacol* 13:1040838 (FDA QSAR BBB models) and the BBB-uncertainty-quantification work in Chemical Information & Modeling demonstrate that logBB-trained classifiers transfer poorly to Kp,uu (the unbound brain:plasma ratio that the free-drug hypothesis requires). The recommended practice is to (i) train/use ADMET-AI for logBB classification as a coarse filter, (ii) layer a regression model — for example Hu, Jiang, Li, Zhang, Wu, Zhang & Zhuang 2025 *Drug Delivery* (DOI 10.1080/10717544.2025.2585612) which reports "a robust correlation between MDR1-derived Papp(A-B) and Kp,uu,brain (R = 0.8886), with the remaining 21 compounds validating predictive accuracy (≤2-fold error)" across 20 training + 21 validation drugs, or the LeiCNS-PK3.0 integration QSPR — for Kp,uu point estimates, (iii) propagate ML uncertainty by Monte Carlo sampling brain:plasma ratios from posterior predictive distributions.

**(F8) MBMA is the right inferential frame; it has CNS precedent.** Mandema, Salinger, Ahn and successors established model-based meta-analysis (MBMA) as the canonical pharmacometric tool for synthesizing summary-level dose-response/time-course data across trials (Mawdsley/Bujkiewicz 2016 *CPT:PSP* 5:393–401 MBNMA framework; Boucher/Bennetts 2018 *Stat Med* time-course MBNMA; Bachhav et al. 2025 *CPT:PSP* anti-Aβ mAb MBMA characterizing aducanumab/lecanemab/donanemab amyloid removal and ARIA-E together). The proposed model is essentially an MBMA whose dose-response surface is informed by upstream V6.A occupancy posteriors rather than dose alone — a hybrid the literature has not previously executed for cognition endpoints.

## Details

### 1. Mechanism-class meta-analytic SMD priors

The proposed model uses a four-level partial-pooling hierarchy: **mechanism class → molecular target → compound → endpoint**, with weakly-informative priors on variance components and class means anchored to the meta-analyses below. All SMDs/Hedges' g are reported for cognition endpoints; where the source reports raw MD, we provide back-calculation guidance.

| Class | Population | Pooled estimate | 95% CI | k / N | Source |
|---|---|---|---|---|---|
| **Cholinergic — AChE-i (donepezil 10 mg)** | Mild-mod AD | ADAS-Cog MD = −2.67 | [−3.31, −2.02] | 5 RCTs / 1,130 | Birks & Harvey 2018 Cochrane CD001190.pub3 |
| Cholinergic — donepezil 10 mg | Severe AD | SIB MD = +5.92 | [4.53, 7.31] | 5 RCTs / 1,348 | Birks & Harvey 2018 |
| Cholinergic — donepezil 10 mg | Any-stage dementia (MMSE) | Hedges' g = 2.27; ADAS-Cog non-significant | [1.25, 3.29] | 18 RCTs / 5,948 | Sheikh & Ammar 2024 *Front Neurosci* 18:1398952 (treat as outlier — see F4 caveat) |
| Cholinergic — nicotine acute | Healthy adults, 9 domains | Small-moderate for attention, motor, episodic/working memory | reported per-domain | 41 RCTs | Heishman et al. 2010 *Psychopharmacology* 210:453–69 |
| **Dopaminergic — modafinil acute** | Healthy non-sleep-deprived adults | Overall SMD = 0.12; memory updating 0.28 | p=0.01 | 14 studies / 64 ES | Roberts 2020 *Eur Neuropsychopharmacol* 38:40–62 |
| Dopaminergic — methylphenidate | Healthy adults | Overall 0.21; recall 0.43; sustained attention 0.42; inhibitory control 0.27 | p=0.0004 | 24 studies / 47 ES | Roberts 2020 |
| Dopaminergic — d-amphetamine | Healthy adults | Null | — | 10 studies / 27 ES | Roberts 2020 |
| Dopaminergic — atomoxetine | Adult ADHD core symptoms | SMD = −0.40 (clinician); −0.33 (patient) | — | 12 RCTs / 3,375 | Cunill et al. 2013 *Pharmacoepidemiol Drug Saf* 22:961–9 |
| Dopaminergic — atomoxetine | Adult ADHD network MA | self-rated SMD = −0.38 [−0.56, −0.21], 95% PI [−0.82, 0.06] | — | — | Cortese et al. 2025 *Lancet Psychiatry* (DOI: 10.1016/S2215-0366(24)00360-2) |
| Dopaminergic — tolcapone (COMT-i) | Healthy adults | **No formal pooled SMD published** | — | Apud 2007, Bhakta 2017, narrative reviews only | gap flagged |
| **Glutamatergic — memantine** | Mild-mod AD | ADAS-Cog SMD = −0.21 | [−0.34, −0.08] | 6 RCTs | Winblad 2007 Cochrane meta-analysis |
| Glutamatergic — memantine | Moderate-severe AD | SMD = −0.29 | [−0.54, −0.03] | trials per Cochrane | Cochrane Memantine review |
| **Monoaminergic / 5-HT — vortioxetine** | MDD (DSST) | SMD = 0.325 | [0.120, 0.529] | 12 RCTs (network MA) | Baune et al. 2018 *Int J Neuropsychopharmacol* 21:97–107 |
| Adrenergic α2A — guanfacine | ADHD pediatric / monkey/human PFC | Improves spatial WM; no formal pooled human-adult SMD | — | Arnsten serial primate studies; pediatric MAs | Avery 2000, Arnsten 2002 *J Neurosci* 22:8771 |
| Adenosinergic — caffeine acute | Healthy rested adults | g = 0.27 (accuracy); g = 0.28 (RT) | dose-response confirmed | 31 RCTs / 1,455 | Kløve & Petersen 2025 *Psychopharmacology* 242(9):1909–1930 |
| Peptide/cognitive enhancer — piracetam | Older adults with cognitive impairment | OR (CGI-C) = 3.55 [2.45, 5.16]; no continuous SMD | I² high | 19 DB-PC RCTs / 1,489 | Waegemans et al. 2002 *Dement Geriatr Cogn Disord* 13:217–224 |
| Mitochondrial — creatine | Healthy adults | Short-term memory, intelligence/reasoning improved; quantitative SMD only in updated meta-analysis | — | 6 RCTs / 281 (Avgerinos); updated Prokopidis 2023 *Nutr Rev* 81:416–427 | Avgerinos 2018 *Exp Gerontol* 108:166–73 |
| Mitochondrial — omega-3 PUFA | Healthy older adults global cognition | Pooled SMD = −0.02 [−0.07, 0.04] (essentially null); MMSE umbrella SMD = 0.16 [0.01, 0.32] | — | 11 RCTs (global); 14 RCTs / 26,881 (MMSE umbrella) | Cooper 2015 *J Psychopharmacol* 29:753–63; Nutrients 2025 |
| Anti-inflammatory — minocycline (schiz) | Schizophrenia | PANSS-total SMD = −0.59 to −0.64; executive function SMD = 0.22 [0.01, 0.44] | — | 6–8 RCTs / 215–548 | Solmi 2017 *CNS Spectr* 22:415–26; Xiang 2017 *Eur Neuropsychopharmacol* 27:8–18 |

**Heterogeneity and publication-bias adjustments.** Most cholinergic and dopaminergic MAs report I² in the 40–75% range, motivating random-effects pooling with τ² priors (e.g., half-normal(0, 0.25²) on τ as recommended by Williams, Rast & Bürkner 2018, *PsyArXiv*). PET-PEESE and Egger-style adjustments are recommended where k ≥ 10. The Roberts 2020 PRISMA pipeline already implements trim-and-fill; the memantine and donepezil Cochrane reviews use the standard Cochrane RoB 2.0 tool.

**Within-class heterogeneity** is non-trivial. The Roberts 2020 result shows methylphenidate effects differ by ~2× across sustained attention (0.42) vs working memory (smaller, not significant). The Bhakta 2017 tolcapone replication failed despite Apud 2007's success — a textbook reminder that compound-level pooling that ignores COMT genotype, baseline DA tone, or task structure is misleading. The model accordingly partial-pools by compound and endpoint, with the class mean acting only as the highest-level prior anchor.

### 2. PK/PD models for brain exposure

The translation function requires four nested PK/PD layers:

1. **Plasma exposure (Cmax, AUC)** from a 1- or 2-compartment population PK model parameterized from published popPK studies. For donepezil: Cmax ≈ 15–30 ng/mL at 5–10 mg/day steady state, t½ ~70 h, V/F ~12 L/kg, F = 100%; Reyes 2002 *Br J Clin Pharmacol*.
2. **Brain:plasma ratio (Kp,uu)** sampled from ML/QSAR posteriors (see §3) or experimental microdialysis estimates where available. For CNS-active drugs the median Kp,uu is ≈ 0.3–1.0 (Friden 2009 *Drug Metab Dispos* 37:1226).
3. **Target occupancy (RO)** via the Hill equation: RO(t) = E_max · C_brain,unbound(t)^n / (EC50^n + C_brain,unbound(t)^n). EC50 derived from PET dose-occupancy data where possible. PET-validated anchors:
   - **Donepezil → AChE:** Bohnen 2005 reports cortical AChE inhibition mean 19.1% (SD 9.4%) at 5 mg/day, with stronger inhibition (39 ± 5%) at PET k₃ measurement (Shinotoh 2001 *Neurology* 56:408–410); Ota 2010 *Clin Neuropharmacol* 33:74 estimates plasma IC50 for cortical AChE inhibition.
   - **Donepezil → sigma-1:** Ishikawa 2009 (PMID 19573265) — 5 mg yields ~60% sigma-1 occupancy, 10 mg ~75%, with [11C]SA4503.
   - **Antipsychotics → D2:** Kapur 2000 *Am J Psychiatry* 157:514–520 — therapeutic D2 occupancy 60–80%; EPS threshold ~78%; prolactin threshold ~72%; the Fitzgerald 2000 *Neuropsychopharmacology* 22:19–26 model relates plasma haloperidol to D2 RO with an E_max-Hill structure.
4. **Effect-size mapping** RO → SMD via class-conditional response surface. Cholinergic/glutamatergic/serotonergic: monotonic sigmoid. Dopaminergic/noradrenergic: inverted-U (Cools & D'Esposito 2011; Williams & Goldman-Rakic 1995 D1 inverted-U primate WM). The inverted-U is encoded as g(RO; θ_optim, σ_width) = β · exp(−(RO − θ_optim)² / 2σ_width²) with baseline-DA covariate centering.

**Cmax vs AUC** debate: for cognition endpoints, single-dose crossover designs (most of Roberts 2020) are driven by Cmax-coupled brain occupancy at test time; chronic dosing trials (most of Birks 2018) are AUC- and steady-state-driven. The translation function must condition on trial design.

**Inter-individual variability sources** to model as random effects: CYP2D6/3A4 status (donepezil, atomoxetine), COMT Val158Met (tolcapone, dopaminergics — Apud 2007, Bhakta 2017), age (BBB integrity, hepatic clearance), ApoE4 (AD relevance, AChE-i response heterogeneity), education/baseline cognitive performance (ceiling effects).

**Tooling:** Population PK via mrgsolve (Baron et al.) or PKPDsim within an R pipeline; Bayesian fit via brms or Stan + Torsten (Margossian 2022); for Julia users Pumas.jl provides integrated NLME-Bayesian workflow (Tarek et al. 2023 arXiv:2304.04752) with the only software outside Stan that does fully Bayesian NLME pharmacometric inference with the breadth Pumas does. NONMEM remains the regulatory standard. Given the V6 stack (PyMC + JAX), the pragmatic choice is PyMC with custom ODE integration or PKPDsim hand-off, accepting a ~3× compute cost vs Stan/Torsten.

### 3. ADMET-AI BBB calibration

ADMET-AI (Swanson et al. 2024 *Bioinformatics* 40:btae416) outputs 41 endpoints including BBB permeability classification. Best-practice calibration:

- Use ADMET-AI BBB score as a Boolean filter (BBB+ vs BBB−) with threshold tuned to balance recall/precision against a held-out Kp,uu ground-truth set (Loryan 2015 *AAPS J*; Friden 2009 *DMD* 37:1226).
- For Kp,uu point estimation, layer a regression model trained on Kp,uu (not logBB). Hu et al. 2025 *Drug Delivery* (DOI 10.1080/10717544.2025.2585612) demonstrated R = 0.8886 between MDR1-derived P_app(A-B) and Kp,uu across 20 training + 21 validation drugs (≤ 2-fold prediction error in the held-out set). The Chen 2024 *J Cheminform* 16:42 graph-neural-network Kp,uu predictor with uncertainty estimation is the current SOTA. The Saleh 2022 LeiCNS-PK3.0 integration (Saleh et al. *Pharm Res*) is the canonical PBPK template.
- Propagate uncertainty by Monte Carlo: at each MCMC iteration draw Kp,uu ~ posterior_from_ML(SMILES), feed into the PK/PD ODE, propagate to occupancy and to predicted Hedges' g.

### 4. Hierarchical Bayesian effect-size model — formal specification

For trial *i* of compound *c*, dose *d_i*, on endpoint *e_i*, in population *p_i*, with mechanism class *m_c*:

```
g_observed_i ~ N(g_predicted_i, SE_i² + σ_residual²)
g_predicted_i = α_class[m_c] + β_target[t_c] · f_class(RO_i; θ_class) + γ_endpoint[e_i] + δ_population[p_i]
RO_i ~ HillEquation(C_brain_unbound_i, EC50_target, n_Hill)
C_brain_unbound_i ~ popPK(dose_i, covariates_i) × Kp,uu_c
Kp,uu_c ~ MLposterior(SMILES_c)         # from ADMET-AI / Kp,uu regressor
β_target[t_c] ~ Normal(0, σ_β)
α_class[m] ~ Normal(α_class_prior_m, σ_class)   # Robust MAP prior from meta-analyses
σ_β, σ_class, σ_residual ~ HalfNormal(0, 0.25)
```

with **robust MAP priors** on `α_class_prior_m` derived per Schmidli et al. 2014 *Biometrics* 70:1023–1032 — a mixture of an informative component (the published meta-analytic posterior for class *m*) and a weakly-informative robust component (Normal(0, 1)) to guard against prior-data conflict. The mixture weight *w* is set per Schmidli's recommendations (typically 0.5, sensitivity-tested in [0.2, 0.8]).

**Cluster D θ̄ gating:** The Cluster D cognition-relevance posterior θ̄_t enters as a multiplicative gate on β_target: `β_target[t_c] = θ̄_{t_c} · β_raw_target[t_c]`. Targets with θ̄ near zero contribute negligibly to predicted g regardless of binding affinity — the V6.B neurobiological prior is encoded directly into translation.

**For dopaminergic / noradrenergic compounds**, swap f_class for an inverted-U:
`f_dopamine(RO) = β · exp(−(RO − θ_optim[baseline_DA])² / 2σ_width²)`
with baseline_DA modeled as a population covariate (default population mean for de novo predictions; specific draws when COMT genotype or age is known).

**Sampling:** PyMC NUTS with JAX backend (already in V6.B stack), target_accept = 0.95, 4 chains × 4,000 draws after 2,000 warmup. Expected runtime on the requester's RTX 5070: ~6–12 h for a model with ~210 targets × ~50 compounds × ~10 endpoints in trial-summary mode.

### 5. Endpoint-specific considerations

**Digit-symbol-substitution test (DSST)** is the canonical cognition-enhancement endpoint because (i) it integrates processing speed, attention, working memory, executive function; (ii) high test-retest reliability (r = 0.69–0.82 across visits; Hilliard 2022 *JMIR Ment Health* 10:e33871); (iii) sensitive to drug effects (vortioxetine SMD = 0.325; methylphenidate SMD ≈ 0.65 in TBI per Huang 2019 *Brain Sci* 9:291); (iv) widely used in MetaPsy and Cochrane databases. Population-level normative SD ≈ 11–17 in healthy older adults (CHS data, Rosano 2016 *Age Ageing* 45:688); ~20 in younger adults (WAIS-IV normative).

**Other endpoints in the panel:** n-back (working memory; high variance, ceiling-prone in healthy), Stroop (interference; sensitive to AChE-i per Bohnen 2005), RAVLT (episodic memory; AD-sensitive), CANTAB RVIP (sustained attention; modafinil-sensitive per Battleday 2015 *Eur Neuropsychopharmacol* 25:1865), and the MCCB battery for schizophrenia-related cognition (Bhakta 2017 used MCCB global composite for tolcapone).

**Endpoint-to-endpoint mapping** uses the network meta-analytic equivalence approach: when MAs report only one endpoint per compound (e.g., DSST for vortioxetine), assume cross-endpoint Pearson correlations of 0.5–0.7 based on Knowles 2005 *Neuropsychologia* DSST/Stroop convergence, and treat γ_endpoint as a fixed-effect offset per endpoint with informative prior centered on the literature mean offset.

**Population SDs** required for back-calculating raw MD → SMD: ADAS-Cog in mild-mod AD ≈ 8–10; MMSE ≈ 3–5 in same population; DSST ≈ 11–20 depending on age band.

**Crossover designs** require correlation-corrected SE per Borenstein 2009; Roberts 2020 handles this explicitly (501 within-subject vs 144 between-subject participants for MPH).

### 6. Uncertainty propagation

The four uncertainty layers, in order:
1. **Measurement noise** (trial-level SE_i, propagated as the data likelihood);
2. **Target uncertainty** (V6.A.4 Venn-ABERS calibrated pchembl posteriors per (compound, target));
3. **PK/PD uncertainty** (Kp,uu draws from ML posterior; popPK random effects; EC50/Hill_n posterior from PET-anchored Bayesian fits);
4. **Mechanism-class prior uncertainty** (the robust MAP mixture variance, including the heavy-tailed robust component).

All four are integrated by Hamiltonian Monte Carlo within a single joint posterior. Output is the posterior predictive distribution of g for novel (compound, dose, endpoint) tuples, summarized as median + 90% HDI.

**Diagnostics:** posterior predictive checks (Gabry et al. 2019 *JRSS-A*), R-hat < 1.01, ESS > 400, PSIS-LOO with shape parameter k_hat < 0.7 for ≥ 95% of observations (Vehtari 2017 *Stat Comput* 27:1413–1432), prequential PIT calibration (Czado 2009 *Biometrics*). Cross-validation via Bürkner's `projpred` for variable selection on which mechanism subclasses to retain as separate hierarchical levels.

### 7. Validation strategy

**Primary validation:** leave-one-trial-out (LOTO) cross-validation across the assembled training corpus. Calibration plot of predicted vs observed Hedges' g. Pre-committed thresholds:

- **Calibration slope ∈ [0.8, 1.2]** (intercept ≈ 0)
- **90% credible interval coverage ≥ 85%** (allowing modest under-coverage from random-effects shrinkage)
- **RMSE on held-out g < 0.15** (below the median within-class τ ≈ 0.2)
- **ELPD-LOO advantage over the mechanism-class-mean naive baseline: ΔELPD > 4 × SE(ΔELPD)** (Sivula et al. 2022 *Bayesian Anal*)

**Secondary validation:** held-out validation against compounds not in training set — e.g., train on all-but-donepezil, predict donepezil at MCI 10 mg on ADAS-Cog, compare against Birks 2018 Cochrane MD = −2.67 (back-calculated SMD ≈ −0.30). Pre-register the exact predicted point estimate and 90% HDI before unblinding.

**Tertiary validation:** prospective predictions for compounds with known Phase 2 readouts not yet integrated (e.g., recent vortioxetine extension trials, atomoxetine in MCI patients), and falsifiability check at one-year horizon.

**Small-N challenge:** typically 1–5 trials per compound at the (compound, dose, endpoint) cell. Partial pooling via the four-level hierarchy is the only feasible response; the model cannot generate calibrated novel-compound g without the class-level borrowing strength.

**MetaPsy** (https://www.metapsy.org/) provides rectangular Hedges'-g data for psychotherapy and some pharmacological trials in MDD; the metapsyTools R package (Cuijpers/Harrer) computes effect sizes from per-trial raw data following Borenstein 2009. It is **usable for the depression-cognition arm (vortioxetine, SSRIs, SNRIs)** but does **not** cover most cholinergic, dopaminergic, or glutamatergic cognition trials — those require manual extraction from Roberts 2020 supplementary tables, Birks 2018 Cochrane forest plots, the AlzForum Therapeutics database, and ClinicalTrials.gov MetaTrial results.

### 8. Related work / state of the art

- **HINT** (Fu et al. 2022 *Patterns* 3(4):100445): hierarchical interaction network predicts binary trial outcomes (F1 0.665/0.620/0.847 across phases). Predicts pass/fail; does not output Hedges' g.
- **Wong/Siah/Lo ML extension** (Lo et al. 2018 SSRN; *Biostatistics* 2019): AUC 0.78/0.81 for Phase 2/3 → approval. Approval probability, not effect size.
- **Quantitative systems pharmacology** (Geerts et al. 2018 *Schizophr Bull* 44:S221; Spiros 2017): cortico-striato-thalamic loop QSP model with 32 CNS targets translates binding profiles to EPS predictions in 1,124 patients on 772 unique antipsychotic combinations from CRIS database. Closest published methodology, but outputs phenotype-incidence-rate, not SMD.
- **MBMA in CNS:** Bachhav et al. 2025 *CPT:PSP* anti-Aβ mAb model (aducanumab + lecanemab + donanemab) for amyloid + ARIA-E. Mawdsley 2016 *CPT:PSP* dose-response MBNMA. Boucher & Bennetts 2018 *Stat Med* time-course MBNMA.
- **Translational PK/PD bridging:** Danhof et al. (multiple) on mechanism-based PK/PD; Hutmacher 2008 *Br J Clin Pharmacol* for placebo-response modeling.
- **No published methodology** translates in silico DTI affinity into predicted clinical SMD with calibration error. This is the contribution.

### 9. Identifiability and prior sensitivity

**Identifiability concerns** are real. With k=47 Roberts 2020 trials anchoring a four-level hierarchy (class → target → compound → endpoint), per-cell N is sparse and variance components are weakly identified. Mitigations:

- **Half-normal(0, 0.25)** priors on all τ (between-trial heterogeneity), justified by typical cognition-MA τ ≈ 0.1–0.3 (Williams/Rast/Bürkner 2018).
- **Partial pooling** is essential — fully unpooled per-compound estimation is impossible at N=1–3 trials/compound; the class prior is what gives the model leverage.
- **Robust MAP mixtures** (Schmidli 2014) on class means: 80% weight on the meta-analytic informative prior, 20% weight on Normal(0, 1) — protects against the case where a novel target inherits an inappropriate class-mean from the literature.
- **Sensitivity analyses to pre-register:** (i) re-fit with τ ~ Half-Cauchy(0, 0.5); (ii) re-fit with non-informative α_class ~ Normal(0, 5); (iii) drop the Cluster D θ̄ gate; (iv) drop the inverted-U for dopaminergics (linear fallback); (v) re-fit with leave-one-class-out.
- **Effective sample size monitoring**: if any α_class has ESS < 400 after 16,000 post-warmup draws, restrict that class hierarchy to fixed-effect rather than random-effect.

Heavy-tailed outliers (e.g., the Apud 2007 tolcapone result, which Bhakta 2017 failed to replicate) motivate Student-t(df=4) residuals at the trial-level likelihood rather than Normal.

### 10. Publication strategy

- **Top target: Cell Reports Medicine.** Recent fits include the Geerts QSP cohorts and translational drug-repurposing methodology pieces. Sells if (i) LOTO validation meets pre-committed thresholds, (ii) at least one prospective Phase-2 prediction comes within 0.1 SMD of observed, (iii) the open-source software release is professional-grade.
- **Strong alternative: Clinical Pharmacology & Therapeutics (CPT).** The 2018 Wagner et al. "dynamic map" framework and FDA MIDD guidance map directly to this scope. CPT is more methods-friendly than Cell Reports Medicine and slightly lower bar for novelty in methodology-only papers.
- **Methods-specialist fallbacks:** *CPT: Pharmacometrics & Systems Pharmacology (PSP)* — natural home given the MBMA + Bayesian pharmacometrics core; Bachhav 2025, Mawdsley 2016, Yu 2026 all in this venue. *Journal of Pharmacokinetics and Pharmacodynamics* — quantitative methods focus. *AAPS Journal* — broader but accessible.
- **Negative-result path:** if LOTO fails the pre-committed thresholds, the paper still reports the framework, the calibration failure, the diagnosed reason (class-prior misspecification, identifiability collapse, occupancy-effect curve mis-shape, etc.) and is publishable in *CPT:PSP* or *PLOS Computational Biology* as a methodology + failed-validation report. The pre-registration on OSF + the negative result is a contribution.

**To clear the Cell Reports Medicine bar:** the validation must include at least one prospective head-to-head against a recent (post-training cutoff) trial readout, with the prediction registered to OSF before unblinding. Anti-amyloid mAbs (lecanemab, donanemab) post-2024 readouts provide a natural prospective test. Pooled CDR-SB SMD ≈ −0.49 [−0.67, −0.30] for combined lecanemab+donanemab (Avgerinos 2024 *Sci Rep* 14:25741); the donepezil 10 mg ADAS-Cog Cochrane MD of −2.67 (back-calculated SMD ≈ −0.30) provides a stationary anchor.

### 11. Practical implementation

**Software stack:**
- PyMC 5.x with `numpyro` JAX backend (matches V6.B) — primary inference engine.
- For PK ODEs: `diffrax` JAX-native ODE solver (Kidger et al.) integrated as a PyMC `Op`, or hand-off to PKPDsim/mrgsolve for popPK simulation with results imported as covariates.
- For BBB/Kp,uu: ADMET-AI v1.x + a custom Kp,uu regressor (training on Loryan 2015 + Friden 2009 + Hu 2025 data).
- For MetaPsy: metapsyTools R + metapsyData (Harrer/Cuijpers) for the depression/cognition arm.
- For reporting/diagnostics: `arviz`, `posterior`, `loo`, `projpred`.

**Compute budget:** RTX 5070 (12 GB VRAM, ~20 TFLOPS FP32). Expected timings:
- 210-target model, trial-summary mode (k≈200 trials across all classes): ~6–12 h for full Bayesian fit with diffrax-backed PK ODEs.
- LOTO validation: ~k×0.5 h with warm-started posteriors. Budget 5–7 days total for full validation.
- ML BBB/Kp,uu retraining: ~1 day on GPU.

**Data assembly pipeline (the V6.D.0 step):**
1. Extract trial-level g + dose + compound + endpoint + population from Roberts 2020 supplementary table, Birks 2018 Cochrane Review Manager XML, Winblad/McShane memantine Cochrane, Heishman 2010, Battleday 2015, Repantis 2010, Baune 2018 vortioxetine, Solmi 2017 minocycline, Cunill 2013/Cortese 2025 atomoxetine, Kløve & Petersen 2025 caffeine, Avgerinos 2018/Prokopidis 2023 creatine, Cooper 2015 omega-3, plus MetaPsy depression-cognition. Manual transcription into a normalized Parquet schema with columns: `trial_id, compound_inchikey, dose_mg, dose_freq, n_treatment, n_control, endpoint_normalized, population_normalized, hedges_g, hedges_g_se, design (parallel/crossover), duration_weeks, mean_age, percent_female, baseline_severity_normalized`.
2. Crosswalk endpoints to a normalized set (ADAS-Cog, MMSE, DSST, n-back, Stroop, RAVLT, CANTAB-RVIP, MCCB-composite, ...) via a manual mapping table.
3. Crosswalk populations to {healthy adult, MCI, mild-mod AD, mod-sev AD, schizophrenia, MDD, adult ADHD, child ADHD, older healthy adult}.
4. Assemble PET-occupancy training data from Bohnen 2005, Ishikawa 2009, Kapur 2000, Ota 2010, Shinotoh 2001 supplementary, plus the Lassen-plot literature for newer compounds.
5. Pull V6.A pchembl posteriors and V6.B θ̄ posteriors for the same compound × target cells.

### 12. Falsifiability and pre-committed predictions

The following predictions are committed to OSF pre-registration before any model fit:

- **(P1) Donepezil 10 mg, mild-mod AD, ADAS-Cog** (held out from training): predicted SMD must fall in (−0.45, −0.15) and 90% HDI must cover Birks 2018 back-calculated SMD ≈ −0.30. *Falsifier: missed entirely (CI outside the interval).*
- **(P2) Modafinil 200 mg, healthy adults, executive function**: predicted g must fall in (0.05, 0.30) and 90% HDI cover Roberts 2020 modafinil-overall SMD 0.12. *Falsifier: HDI excludes 0.12.*
- **(P3) Methylphenidate 20 mg, healthy adults, sustained attention**: predicted g must fall in (0.25, 0.60) and 90% HDI cover Roberts 2020 SMD 0.42. *Falsifier: HDI excludes 0.42.*
- **(P4) Vortioxetine 20 mg, MDD, DSST**: predicted g must fall in (0.15, 0.50) and 90% HDI cover Baune 2018 SMD 0.325. *Falsifier: HDI excludes 0.325.*
- **(P5) Caffeine 200 mg, healthy adults, attention/RT**: predicted g must fall in (0.15, 0.40) and 90% HDI cover Kløve & Petersen 2025 g 0.27/0.28. *Falsifier: HDI excludes both.*
- **(P6) Negative control (d-amphetamine, healthy adults, any endpoint)**: predicted g 90% HDI must include 0 (per Roberts 2020 null). *Falsifier: HDI excludes 0 systematically.*
- **(P7) Calibration slope, all predictions, vs all observed**: must be in [0.8, 1.2]. *Falsifier: outside.*
- **(P8) ELPD-LOO vs class-mean baseline**: ΔELPD > 4 × SE(ΔELPD). *Falsifier: model fails to beat the simple baseline → publish as negative result in CPT:PSP.*

**Negative-result publication path:** if (P7) or (P8) fails, the paper documents the framework, the data assembly, the diagnostic story, and submits to *CPT:PSP* as a methods + null-validation report. This is itself a contribution to the field because no prior published methodology has even attempted this translation with calibration error.

## Recommendations

**Stage 1 — Prior elicitation (Weeks 1–6, V6.D.1).** Build the mechanism-class meta-analytic prior table (extending the table in §1). Manually extract trial-level data from Roberts 2020, Birks 2018 Cochrane, memantine Cochrane, Heishman 2010, Baune 2018, Cunill/Cortese atomoxetine, Solmi/Xiang minocycline. Compute Hedges' g per trial. Fit class-level random-effects models in `brms` or `metafor`. Output: Parquet table of (class, target, compound, endpoint, population, g, SE, design) and a posterior for α_class[m] per class. **Decision gate:** if the class-mean posteriors are unstable (90% HDI width > 0.4 for any class), expand the source MA library before proceeding.

**Stage 2 — PK/PD layer assembly (Weeks 5–12, V6.D.2).** Pull popPK parameter posteriors from published popPK studies for the 10–20 highest-priority compounds (donepezil, modafinil, MPH, atomoxetine, memantine, vortioxetine, caffeine, nicotine, guanfacine, tolcapone, etc.). Assemble PET dose-occupancy datasets from Bohnen 2005, Ishikawa 2009, Kapur 2000, Ota 2010, varenicline α4β2 occupancy (Lotfipour 2013), and others. Fit Hill curves per (target, compound) in PyMC. Train the Kp,uu regressor on Loryan/Friden/Hu data with conformal uncertainty. **Decision gate:** if no occupancy-effect anchor exists for a given mechanism class (e.g., adenosinergic — no caffeine PET A2A occupancy-effect curve), document and either commission a literature-derived prior or restrict the class to a flat dose-response.

**Stage 3 — Hierarchical model fit (Weeks 12–18, V6.D.3).** Integrate Stage 1 priors + Stage 2 PK/PD + V6.A.4 pchembl posteriors + V6.B θ̄ posteriors into the full Bayesian model. Run NUTS in PyMC with diffrax PK ODEs. Diagnostics: R-hat, ESS, PSIS-LOO, posterior predictive checks. **Decision gate:** R-hat < 1.01 and ESS > 400 for all class-level parameters; otherwise reparameterize.

**Stage 4 — Validation (Weeks 18–24, V6.D.4).** LOTO cross-validation; pre-registered prospective predictions (P1–P8); calibration plots; coverage diagnostics. Pre-register all of this on OSF before running the validation. **Decision gate:** if (P7), (P8), or the LOTO RMSE thresholds fail, pivot to negative-result publication in *CPT:PSP*. If they pass, proceed to manuscript drafting for Cell Reports Medicine.

**Stage 5 (stretch, Weeks 24–32) — Prospective external validation.** Wait for a post-cutoff readout (lecanemab open-label-extension cognition, donanemab subset, a recent Phase 2 atomoxetine-in-MCI readout, etc.) and report against the OSF-locked predictions. If the prospective validation also passes, Cell Reports Medicine becomes a credible submission target.

**Thresholds that change recommendations:**
- If any class has fewer than ~3 published meta-analyses worth of trials (e.g., orexinergic, cannabinoid, peptide), collapse it into a sibling class or drop from the predictor space rather than fit a poorly-identified hierarchy.
- If V6.A pchembl posteriors and V6.B θ̄ posteriors are themselves not yet validated (per the V6.A.4 / V6.B.0 Gate criteria), defer Stage 3 until upstream gates pass.
- If GPU memory becomes a bottleneck at 210 targets × ~50 compounds × ~10 endpoints, profile and prune to the top 50 most-trial-supported targets first; expand once the small-N model is calibrated.

## Caveats

**(C1) Identifiability with k=47 + sparse-N cells is genuinely hard.** The model leans heavily on the robust MAP priors and class-level partial pooling; if those priors are misspecified (wrong mechanism class assignment, ignored within-class heterogeneity like the methylphenidate working-memory vs attention split), predictions inherit the bias. Sensitivity analyses (§9) are mandatory rather than optional.

**(C2) The mechanism-class taxonomy is an analyst choice, not a discovery.** A compound like donepezil is cholinergic primary (AChE-i) but also sigma-1 agonist (60–75% occupancy at 5–10 mg per Ishikawa 2009). Multi-class compounds violate the single-class prior. The model handles this by allowing β_target to differ per target with weak class-level shrinkage, but the analyst must decide class assignment.

**(C3) PET occupancy → cognition mapping is sparse outside cholinergic and dopaminergic.** For glutamatergic (memantine NMDA partial-block), serotonergic (5-HT modulators), adenosinergic (caffeine A2A), and α2A-noradrenergic (guanfacine) compounds, the PET dose-occupancy literature is thin or absent. The class-specific occupancy → effect curve will be priored almost entirely from theory/preclinical for these classes, with high posterior uncertainty.

**(C4) ADMET-AI and ML Kp,uu predictors have applicability domain limits.** Niraparib, alectinib, and other novel scaffolds show good Kp,uu predictions; truly novel chemistry outside training-set chemical space inherits high but possibly miscalibrated uncertainty. Conformal-prediction wrappers (Chen 2024) help but do not fully solve the OOD problem.

**(C5) Roberts 2020 covers acute healthy-adult administration only.** Most clinically interesting predictions (donepezil chronic in MCI/AD; atomoxetine chronic in adult ADHD) are chronic-dosing trials that Roberts 2020 does not cover. The model integrates Birks 2018 and other chronic-dosing MAs explicitly, but the resulting class-mean prior is a weighted mix of acute/chronic; the population × duration interaction is a known source of unmodeled heterogeneity.

**(C6) The inverted-U for dopamine is parameterized from primate + small-N human studies (Cools & D'Esposito 2011; Narayanan et al. 2022 *Behav Neurosci* 136(3):207).** The optimal-DA point θ_optim and width σ_width are weakly identified empirically and may differ across cognitive domains (e.g., WM-updating vs WM-maintenance) and individuals (COMT genotype). The model treats these as priors but cannot estimate them precisely without explicit baseline-DA measurement, which the V4–V6 stack does not provide. Class-level prior on θ_optim is anchored at ~70% in the absence of better data, with broad uncertainty. The Narayanan meta-analysis sets a ceiling on the expected explained variance: PFC-DA manipulations explain ~10% of WM-behavior variance, D1DR manipulations ~26%, implying that even with perfect occupancy estimation a large unmodeled residual will remain.

**(C7) The 99.6% AD drug-development failure rate** (Cummings 2014) means that any class-prior anchored on AD trials is selecting on the small successful tail. The robust MAP mixture partially handles this; users should not over-interpret the model's confidence in novel AD predictions.

**(C8) MetaPsy does not cover most of the cognition-enhancement literature.** Its strength is psychotherapy and antidepressant trials in MDD; only the vortioxetine/SSRI arm of this proposal benefits from automated MetaPsy ingestion. Manual extraction is required for the rest.

**(C9) The translation function does not replace wet-lab validation.** It produces a calibrated prior over expected clinical effect-size for a (compound, dose, endpoint) cell. Wet-lab assay confirmation of the upstream pchembl prediction is still required before a clinical trial. The contribution is to refine the rank → prior step, not to skip experimentation.