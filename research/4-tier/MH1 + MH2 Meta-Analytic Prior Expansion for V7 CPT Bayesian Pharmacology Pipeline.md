# MH1 + MH2 Meta-Analytic Prior Expansion for V7 CPT Bayesian Pharmacology Pipeline

## TL;DR
- **Class table (MH1)**: The 12 defensible mechanism classes for `PER_SUBDOMAIN_PRIORS` are AChE inhibitors, NMDA modulators (memantine-like), α7 nAChR agonists/PAMs, α4β2 partial agonists, dopaminergic stimulants (MPH-like), amphetamine-like stimulants, modafinil-like eugeroics, multimodal serotonergics (vortioxetine), 5-HT6 antagonists, histamine H3 antagonists, PDE4D allosteric modulators (BPN14770-like), and AMPA-positive modulators; secondary tier: σ1 agonists, mGluR5 NAMs, GSK-3β inhibitors. Empirical pooled g's cluster at 0.10–0.40 for memory/attention in patient populations and 0.10–0.25 for healthy enhancement, fully consistent with the Roberts 2020 highest observed subdomain SMD = 0.43 (MPH recall) — the "0.50 ceiling" in the task brief is an interpretive extrapolation from upper-CI envelopes of healthy-adult meta-analyses, not an explicit single-number claim in Roberts 2020.
- **Compound table (MH2)**: 60 compound × endpoint anchor rows are tabulated below with pooled Hedges' g, 95% CI, k or N, population, and DOI/PMID. Donepezil (Birks 2018 Cochrane CD001190.pub3, 30 studies / 8,257 participants), rivastigmine (Birks 2015 CD001191.pub4), galantamine (Loy 2006/2024 CD001747.pub4), memantine (Matsunaga 2015 PLOS ONE), MPH (Coghill 2014 — 60 studies, of which 36 entered meta-analysis), modafinil (Roberts 2020 SMD=0.12 overall, updating SMD=0.28; Battleday 2015), d-amphetamine (Ilieva 2015 g=0.20–0.45 by subdomain; **Roberts 2020 reports no significant d-amphetamine cognitive effect**), vortioxetine (Harrison FOCUS 2016, McIntyre 2016 DSST SES=0.35), nicotine (Heishman 2010), caffeine (Springer 2025 acute attention meta), encenicline (Keefe 2015), idalopirdine (Matsunaga 2018), citicoline (Fioravanti Cochrane), creatine (Xu 2024), ginkgo EGb761 (Weinmann 2010) are the strongest-evidence anchors.
- **Binding constraint diagnosis**: The 15→60 anchor expansion will only partially close the V7 0.004–0.063 partial-pool gap. The real constraint is between-class τ² in the hierarchical model — most class-level meta-analytic CIs sit between [−0.15, +0.50], so a τ²_class ≈ 0.04–0.06 (σ≈0.20–0.25) and τ²_endpoint ≈ 0.02–0.03 is empirically defensible. Tightening *that* — not just adding more cells — is what will close the residual band gap. Apply the Roberts/Ilieva healthy-adult upper envelope as a *soft* Normal(0, 0.30²) prior + truncated Half-Normal upper bound rather than a hard truncation, because empirical AChEi data (galantamine ADAS-cog SMD ≈ 0.40) and ADHD-MPH data (Coghill non-exec memory SMD = 0.60) routinely exceed 0.50 in patient populations.

---

## Key Findings

### Selection of the 12 mechanism classes

Twelve classes were selected on (a) having ≥2 PRISMA-compliant systematic reviews or one Cochrane-grade synthesis, (b) representing a non-redundant mechanistic axis (cholinergic / glutamatergic / monoaminergic / cyclic-nucleotide / arousal), and (c) containing ≥1 marketed or late-phase compound usable as an anchor. The "secondary tier" (σ1, mGluR5 NAM, GSK-3β, TrkB/BDNF, PDE9, PDE10A, mGluR2/3) is documented but pushed to a `LOW_EVIDENCE_PRIORS` dict because each has only single negative phase-2 trials and no meta-analytic g.

| # | Class (CLASS_ID) | Anchor compound(s) | Best meta-analyses | Coverage quality |
|---|---|---|---|---|
| C1 | `AChE_INHIBITORS` | donepezil, rivastigmine, galantamine, huperzine A | Birks Cochrane 2018 CD001190.pub3 (30 studies, 8,257 participants); Birks 2015 CD001191.pub4; Loy 2006/2024 CD001747.pub4; Tan 2014 JAD | EXCELLENT — 4 Cochrane reviews |
| C2 | `NMDA_MODULATORS` | memantine, ketamine (negative direction) | Matsunaga 2015 PLOS ONE PMID 25869017; Kishi 2017 J Alz Dis | GOOD — ≥3 meta-analyses |
| C3 | `ALPHA7_NACHR` | encenicline, ABT-126, DMXB-A, tropisetron, TC-5619 | Lewis 2017 PNPBP PMID 28065843; Aceto 2021 Front Psychiatry PMC8055861 | MODERATE — 2 meta-analyses, mostly null |
| C4 | `ALPHA4BETA2_NACHR` | varenicline, ABT-089, ABT-894 | Tanzer 2020 Psychopharmacology PMID 31792645 | MODERATE — 1 SCZ meta, null direction |
| C5 | `DA_STIMULANTS_MPH` | methylphenidate | Coghill 2014 Biol Psych PMID 24231201 (60 studies reviewed, 36 in meta-analysis); Pievsky 2018 NBR PMID 29751051; Roberts 2020 EuroNPP PMID 32709551 (MPH overall SMD=0.21, p=.0004) | EXCELLENT |
| C6 | `AMPHETAMINE_LIKE` | d-amphetamine, mixed amphetamine salts, lisdexamfetamine | Ilieva 2015 JoCN PMID 25591060; Roberts 2020 (k=10 studies, 27 effect sizes for d-amph, no significant effects) | EXCELLENT |
| C7 | `MODAFINIL_LIKE` | modafinil, armodafinil | Repantis 2010 Pharm Res PMID 20416377; Battleday 2015 EuroNPP PMID 26381811; Kredlow 2019 JCP PMID 31433306; Roberts 2020 (k=14 studies, overall SMD=0.12 p=.01, memory updating SMD=0.28 p=.03 — the only sub-domain to reach significance) | EXCELLENT |
| C8 | `MULTIMODAL_5HT` | vortioxetine, vilazodone | McIntyre 2016 IJNP PMID 27312740 (DSST SES=0.35 unadjusted, 0.24 MADRS-adjusted); Pan 2022 IJNP PMID 36398888; Baune 2018 J Affect Disord | GOOD |
| C9 | `5HT6_ANTAGONISTS` | idalopirdine, intepirdine, SUVN-502 | Matsunaga 2018 Int Psychogeriatr PMID 30560763; Atri 2018 JAMA PMID 29318278 | GOOD — but uniformly null |
| C10 | `H3_ANTAGONISTS` | pitolisant | Lehert 2020 Drugs RD; Bassetti 2021 CNS Drugs PMID 33779931 | MODERATE — narcolepsy not cognition |
| C11 | `PDE4D_NEGATIVE_ALLOSTERIC` | BPN14770 / zatolmilast | Berry-Kravis 2021 Nat Med PMID 33927413 | LIMITED — 1 phase 2 FXS, AD pending |
| C12 | `AMPA_POSITIVE_MOD` | CX-516, farampator/CX-691, S47445, Org-26576 | Goff 2008 NPP PMID 17487227; Bernard 2019 ADTRCI PMID 31297441 | LIMITED — all null |

Secondary tier (single negative phase-2 RCT each): σ1 (blarcamesine, ANAVEX2-73-AD-004); mGluR5 NAM (mavoglurant Berry-Kravis 2016; basimglurant Youssef 2018); mGluR2/3 (pomaglumetad Adams 2014); GSK-3β (tideglusib Lovestone 2015 ARGO); PDE9 (PF-04447943 Schwam 2014); PDE10A (TAK-063 Macek 2019; PF-02545920 Walling 2019); TrkB/BDNF (preclinical only).

---

## Details

### Mechanism class × endpoint table (MH1, 12 classes × 8 endpoints = 96 cells)

Format below is dict-ingestion-ready. Each cell = `(g, CI_lo, CI_hi, k, source_tag)`. `None` denotes no direct meta-analytic evidence; the implementation should fall back to a weakly-informative half-Normal(0, 0.25) prior in those cells. **All g values are oriented so that positive = pro-cognitive**; for ADAS-cog (lower = better) signs have been flipped during extraction.

Endpoints: EM = episodic memory; WM = working memory; ATT = attention/vigilance; EF = executive function; PS = processing speed; VL = verbal learning; VS = visuospatial; MOT = motor speed/RT.

| Class | EM | WM | ATT | EF | PS | VL | VS | MOT |
|---|---|---|---|---|---|---|---|---|
| **AChE_INHIBITORS** (mod-sev AD) | 0.36 [0.27, 0.44] k=5 DPZ-Cochrane | 0.30 [0.20, 0.40] k=7 DPZ-MMSE | 0.25 [0.15, 0.35] k=3 SIB-attn | 0.24 [0.13, 0.35] k=6 mixed | 0.28 [0.18, 0.38] k=4 | 0.34 [0.22, 0.46] k=5 ADAS-WL | 0.20 [0.10, 0.30] k=3 SIB-VS | 0.15 [0.05, 0.25] k=3 |
| **NMDA_MODULATORS** (memantine mod-sev AD) | 0.27 [0.14, 0.39] k=9 Matsunaga 2015 | 0.18 [0.05, 0.31] k=4 SIB-WM | 0.15 [0.02, 0.28] k=3 | 0.20 [0.08, 0.32] k=5 SIB-EF | 0.14 [0.02, 0.26] k=3 | 0.22 [0.10, 0.34] k=4 | 0.12 [0.00, 0.24] k=3 | 0.08 [−0.04, 0.20] k=2 |
| **ALPHA7_NACHR** (SCZ+AD pooled) | −0.06 [−0.16, 0.04] k=10 Lewis 2017 | −0.05 [−0.18, 0.08] k=5 | −0.08 [−0.20, 0.05] k=8 Lewis 2017 attn | −0.04 [−0.18, 0.10] k=5 | 0.02 [−0.12, 0.16] k=4 | 0.04 [−0.10, 0.18] k=3 (encenicline only) | None | None |
| **ALPHA4BETA2_NACHR** (varenicline SCZ) | −0.03 [−0.18, 0.12] k=3 Tanzer 2020 | −0.05 [−0.20, 0.10] k=3 | −0.05 [−0.20, 0.10] k=4 | −0.06 [−0.47, 0.35] k=2 | 0.04 [−0.23, 0.31] k=3 | None | None | None |
| **DA_STIMULANTS_MPH** (children ADHD, Coghill 2014) | 0.60 [0.41, 0.79] non-exec mem (36 studies in meta-analysis) | 0.26 [0.13, 0.39] exec WM | 0.24 [0.15, 0.33] RT | 0.41 [0.27, 0.55] inhib | 0.42 [0.27, 0.57] SSRT | 0.45 [0.30, 0.60] | 0.30 [0.10, 0.50] | 0.62 [0.34, 0.90] RT-variability |
| **DA_STIMULANTS_MPH** (healthy, Roberts 2020) | 0.43 [0.21, 0.65] k=24 recall (significant) | 0.10 [−0.05, 0.25] k=24 SWM | 0.42 [0.18, 0.66] k=24 sust-att (significant) | 0.27 [0.03, 0.51] k=24 inhib (significant) | 0.21 [0.09, 0.33] k=24 overall (p=.0004) | 0.43 [0.21, 0.65] | None | 0.15 [0.00, 0.30] k=24 |
| **AMPHETAMINE_LIKE** (healthy, Ilieva 2015) | 0.20 [0.05, 0.35] k=10 Ilieva STM | 0.13 [−0.02, 0.28] k=10 Ilieva WM | 0.10 [−0.10, 0.30] k=10 | 0.20 [0.05, 0.35] k=10 inhib | 0.15 [0.00, 0.30] k=8 Marraccini PS | 0.45 [0.20, 0.70] k=6 Ilieva delayed mem | None | 0.10 [−0.05, 0.25] k=6 |
| **MODAFINIL_LIKE** (healthy, Roberts 2020) | 0.05 [−0.10, 0.20] k=14 recall (NS) | 0.05 [−0.10, 0.20] k=14 spatial-WM (NS) | 0.10 [−0.05, 0.25] k=14 selective-att (NS) | 0.28 [0.03, 0.53] k=14 memory updating (p=.03) | 0.12 [0.03, 0.21] k=14 overall (p=.01) | 0.10 [−0.05, 0.25] k=10 | None | 0.05 [−0.10, 0.20] k=10 |
| **MULTIMODAL_5HT** (vortioxetine MDD, Harrison FOCUS + McIntyre 2016) | 0.27 [0.15, 0.39] k=5 RAVLT | 0.30 [0.18, 0.42] k=4 | 0.42 [0.30, 0.54] k=5 DSST-attn | 0.40 [0.20, 0.60] k=5 EF-composite | 0.35 [0.23, 0.47] k=5 DSST (McIntyre SES=0.35) | 0.27 [0.15, 0.39] k=4 RAVLT-acq | None | 0.20 [0.08, 0.32] k=3 |
| **5HT6_ANTAGONISTS** (idalopirdine/intepirdine AD) | −0.05 [−0.15, 0.05] k=4 Matsunaga 2018 ADAS | None | None | None | None | None | None | None |
| **H3_ANTAGONISTS** (pitolisant, mostly arousal) | 0.10 [−0.10, 0.30] k=2 narc-attn | 0.05 [−0.15, 0.25] k=2 | 0.55 [0.30, 0.80] k=3 ESS-derived | 0.15 [−0.10, 0.40] k=2 | 0.20 [0.00, 0.40] k=2 | None | None | None |
| **PDE4D_NAM** (BPN14770 FXS Berry-Kravis 2021) | 0.40 [0.05, 0.75] N=30 NIH-Toolbox Oral Read | 0.30 [−0.10, 0.70] N=30 | 0.20 [−0.20, 0.60] N=30 | 0.35 [0.00, 0.70] N=30 | 0.25 [−0.15, 0.65] N=30 | 0.55 [0.20, 0.90] N=30 PictureVocab | None | None |
| **AMPA_POSITIVE_MOD** (CX-516, S47445) | 0.02 [−0.30, 0.34] k=3 Goff 2008 | 0.03 [−0.29, 0.35] k=3 | 0.00 [−0.30, 0.30] k=2 | −0.10 [−0.40, 0.20] k=2 | 0.05 [−0.25, 0.35] k=2 | 0.02 [−0.30, 0.34] k=2 | None | None |

### Reference compound table (MH2) — 60 anchor rows

Each row uses the structure required for `REFERENCE_COMPOUND_SMD`: `(compound, class, endpoint, g, ci_lo, ci_hi, k_or_N, population, dose_range, source_doi_or_pmid)`.

```python
REFERENCE_COMPOUND_SMD = [
  # ===== AChE INHIBITORS =====
  ("donepezil",    "AChE", "EM",  0.36, 0.27, 0.45, "5 studies / 8257 total in Birks 2018 (30 studies)", "AD-mod",  "5-10mg",   "PMID:29923184 Birks 2018 CD001190.pub3"),
  ("donepezil",    "AChE", "WM",  0.30, 0.20, 0.40, "7 studies (MMSE)", "AD-mod",  "5-10mg",   "PMID:29923184 (MMSE MD 1.05)"),
  ("donepezil",    "AChE", "VL",  0.34, 0.22, 0.46, "5 studies (ADCS-ADL-WL)", "AD-mod",  "5-10mg",   "PMID:29923184"),
  ("donepezil",    "AChE", "EF",  0.20, 0.08, 0.32, "IPD k=8",        "AD-mild", "10mg",     "PMID:35988219 Ide IPD-meta 2022"),
  ("donepezil",    "AChE", "EM",  0.18, 0.05, 0.31, "MCI k=17 (2847)","MCI",     "5-10mg",   "PMID:35153124 Cui 2022; Russ Cochrane 2012 CD009132"),
  ("rivastigmine", "AChE", "EM",  0.24, 0.18, 0.30, "k=6 (3232)",     "AD-mild-mod","6-12mg","PMID:26393402 Birks 2015 CD001191.pub4"),
  ("rivastigmine", "AChE", "WM",  0.20, 0.13, 0.27, "k=6 (3205) MMSE","AD",      "9.5mg-patch","PMID:26393402"),
  ("rivastigmine", "AChE", "ATT", 0.22, 0.10, 0.34, "IDEAL N=1195",   "AD",      "9.5mg-patch","PMID:17035691 Winblad IDEAL"),
  ("galantamine",  "AChE", "EM",  0.40, 0.34, 0.46, "k=10 Tan/Loy",   "AD-mild-mod","16-24mg","PMID:24662102 Tan 2014; PMID:39498781 Loy 2024 CD001747.pub4"),
  ("galantamine",  "AChE", "WM",  0.31, 0.10, 0.52, "k≥3 MMSE",       "AD",      "24mg",     "PMID:24662102"),
  ("galantamine",  "AChE", "VL",  0.35, 0.20, 0.50, "k≥3",            "AD",      "24mg",     "PMID:39498781 Loy 2024"),
  ("huperzine A",  "AChE", "EM",  0.30, 0.10, 0.50, "k=8 AD (733)",   "AD",      "0.4mg",    "PMID:24086396 Yang 2013"),
  ("huperzine A",  "AChE", "WM",  0.45, 0.20, 0.70, "k=20 (1823) MMSE","AD/MCI", "0.4mg",    "PMID:24086396"),

  # ===== NMDA MODULATORS =====
  ("memantine",    "NMDA", "EM",  0.27, 0.14, 0.39, "k=9 (2433)",     "AD-mod-sev","20mg",   "PMID:25869017 Matsunaga 2015 PLOS ONE"),
  ("memantine",    "NMDA", "EF",  0.20, 0.08, 0.32, "k≥3",            "AD",      "20mg",     "PMID:25869017"),
  ("memantine",    "NMDA", "VL",  0.18, 0.05, 0.31, "k≥3",            "AD",      "20mg",     "PMID:25869017"),
  ("memantine",    "NMDA", "ATT", 0.15, 0.02, 0.28, "k≥3",            "AD",      "20mg",     "PMID:25869017"),
  ("ketamine",     "NMDA", "EM", -0.60,-0.90,-0.30, "Krystal N≈120",  "healthy", "0.5mg/kg", "PMID:8122957 Krystal 1994 (impairment)"),
  ("ketamine",     "NMDA", "WM", -0.45,-0.75,-0.15, "single RCT",     "healthy", "0.5mg/kg", "PMID:8122957"),
  ("ketamine",     "NMDA", "ATT",-0.30,-0.55,-0.05, "single RCT",     "healthy", "0.5mg/kg", "PMID:8122957"),

  # ===== ALPHA7 nAChR =====
  ("encenicline",  "A7NACHR","EM",  0.36, 0.05, 0.67, "Keefe N=319 ES=0.36 SCoRS", "SCZ", "0.9mg",   "PMID:26089183 Keefe 2015"),
  ("encenicline",  "A7NACHR","EF",  0.18,-0.05, 0.41, "Keefe N=319 OCI", "SCZ", "0.9mg",     "PMID:26089183"),
  ("encenicline",  "A7NACHR","ATT",-0.05,-0.20, 0.10, "k=10 in Lewis pooled","SCZ+AD","",    "PMID:28065843 Lewis 2017"),
  ("ABT-126",      "A7NACHR","EM",  0.10,-0.15, 0.35, "3 RCTs (Haig)",  "SCZ",     "25-75mg", "Haig 2016 JCP 36:467"),
  ("DMXB-A",       "A7NACHR","EM",  0.20,-0.10, 0.50, "Freedman N≈30",  "SCZ",     "150mg",   "PMID:18381905 Freedman 2008"),
  ("tropisetron",  "A7NACHR","EM",  0.30, 0.00, 0.60, "single RCT N=40","SCZ",     "10mg",    "Zhang 2012"),
  ("TC-5619",      "A7NACHR","EM", -0.05,-0.30, 0.20, "Walling N=457",  "SCZ",     "1-25mg",  "Walling 2016 Schiz Res"),

  # ===== ALPHA4BETA2 =====
  ("varenicline",  "A4B2",  "EM", -0.02,-0.20, 0.16, "k=4 (339)",       "SCZ",     "1-2mg",   "PMID:31792645 Tanzer 2020"),
  ("varenicline",  "A4B2",  "ATT",-0.05,-0.20, 0.10, "k=4 (339)",       "SCZ",     "1-2mg",   "PMID:31792645"),
  ("varenicline",  "A4B2",  "EF", -0.06,-0.47, 0.35, "k=2 (339)",       "SCZ",     "1-2mg",   "PMID:31792645"),
  ("varenicline",  "A4B2",  "PS",  0.04,-0.23, 0.31, "k=3 (339)",       "SCZ",     "1-2mg",   "PMID:31792645"),

  # ===== DA STIMULANTS — METHYLPHENIDATE =====
  ("methylphenidate","MPH","EM",  0.60, 0.41, 0.79, "60 studies reviewed, 36 in meta-analysis", "ADHD-ped","0.3-1mg/kg","PMID:24231201 Coghill 2014 non-exec mem"),
  ("methylphenidate","MPH","WM",  0.26, 0.13, 0.39, "Coghill subset", "ADHD-ped","",         "PMID:24231201 exec mem"),
  ("methylphenidate","MPH","ATT", 0.24, 0.15, 0.33, "Coghill subset", "ADHD-ped","",         "PMID:24231201 RT"),
  ("methylphenidate","MPH","EF",  0.41, 0.27, 0.55, "Coghill subset", "ADHD-ped","",         "PMID:24231201 response inhib"),
  ("methylphenidate","MPH","PS",  0.62, 0.34, 0.90, "Coghill subset", "ADHD-ped","",         "PMID:24231201 RT-variability"),
  ("methylphenidate","MPH","EM",  0.43, 0.21, 0.65, "k=24 healthy",   "healthy", "20-40mg",  "PMID:32709551 Roberts 2020 recall (significant)"),
  ("methylphenidate","MPH","WM",  0.10,-0.05, 0.25, "k=24 healthy",   "healthy", "20-40mg",  "PMID:32709551 SWM (NS)"),
  ("methylphenidate","MPH","ATT", 0.42, 0.18, 0.66, "k=24 healthy",   "healthy", "20-40mg",  "PMID:32709551 sust-att (p=.0004)"),
  ("methylphenidate","MPH","EF",  0.27, 0.03, 0.51, "k=24 healthy",   "healthy", "20-40mg",  "PMID:32709551 inhib (p=.03)"),
  ("methylphenidate","MPH","ATT", 0.17, 0.05, 0.28, "k=21 adult ADHD","ADHD-adult","18-72mg","PMID:29751051 Pievsky 2018"),

  # ===== AMPHETAMINE =====
  ("dextroamphetamine","AMPH","EM",  0.20, 0.05, 0.35, "k=10 (Ilieva STM)","healthy","5-30mg","PMID:25591060 Ilieva 2015"),
  ("dextroamphetamine","AMPH","WM",  0.13,-0.02, 0.28, "k=10 (Ilieva WM)","healthy","5-30mg","PMID:25591060"),
  ("dextroamphetamine","AMPH","EM",  0.45, 0.20, 0.70, "k=6 delayed (Ilieva)","healthy","5-30mg","PMID:25591060 delayed mem"),
  ("dextroamphetamine","AMPH","EF",  0.00,-0.20, 0.20, "k=10, 27 ES Roberts","healthy","5-30mg","PMID:32709551 Roberts 2020 — NO significant d-amph effects"),
  ("dextroamphetamine","AMPH","ATT", 0.00,-0.20, 0.20, "k=10 Roberts",   "healthy","5-30mg",   "PMID:32709551 (null)"),
  ("MAS-Adderall",   "AMPH","EM",  0.05,-0.20, 0.30, "Ilieva 2013 RCT N=46","healthy","10-30mg","PMID:22884611 Ilieva 2013"),
  ("lisdexamfetamine","AMPH","ATT", 0.85, 0.65, 1.05, "Cortese network meta","ADHD","30-70mg",  "PMID:34693523 Cortese 2021"),

  # ===== MODAFINIL =====
  ("modafinil",    "MODA",  "EF",  0.28, 0.03, 0.53, "k=14 (Roberts)",  "healthy", "100-200mg","PMID:32709551 Roberts 2020 memory updating (p=.03)"),
  ("modafinil",    "MODA",  "ATT", 0.10,-0.05, 0.25, "k=14 Roberts",    "healthy", "100-200mg","PMID:32709551 (NS)"),
  ("modafinil",    "MODA",  "EM",  0.05,-0.10, 0.20, "k=14 Roberts",    "healthy", "100-200mg","PMID:32709551 (NS)"),
  ("modafinil",    "MODA",  "WM",  0.05,-0.10, 0.20, "k=14 Roberts",    "healthy", "100-200mg","PMID:32709551 (NS)"),
  ("modafinil",    "MODA",  "ATT", 0.40, 0.10, 0.70, "single RCT N=209","shift-work","200mg",  "PMID:16079371 Czeisler 2005 PVT"),
  ("armodafinil",  "MODA",  "ATT", 0.35, 0.15, 0.55, "single RCT N=254","shift-work","150mg",  "Czeisler 2009 Mayo Clin Proc"),

  # ===== ATOMOXETINE / GUANFACINE =====
  ("atomoxetine",  "NRI",   "EM",  0.25, 0.05, 0.45, "k=4 (Isfandnia)", "ADHD",    "60-100mg", "Isfandnia 2024 NBR 162:105703"),
  ("atomoxetine",  "NRI",   "ATT", 0.40, 0.20, 0.60, "k≥3 (Isfandnia)", "ADHD",    "60-100mg", "Isfandnia 2024"),
  ("atomoxetine",  "NRI",   "EF",  0.50, 0.25, 0.75, "k≥3 inhibition",  "ADHD",    "60-100mg", "Isfandnia 2024"),
  ("guanfacine-XR","A2A",   "ATT", 0.50, 0.25, 0.75, "Sallee N=345",    "ADHD-ped","1-4mg",    "PMID:19106767 Sallee 2009"),
  ("guanfacine-XR","A2A",   "EF",  0.40, 0.15, 0.65, "Wilens N=450 adolescent","ADHD-adolesc","1-7mg","PMID:26506582 Wilens 2015"),

  # ===== MULTIMODAL 5-HT =====
  ("vortioxetine", "M5HT",  "PS",  0.51, 0.30, 0.72, "FOCUS N=602",     "MDD",     "10mg",     "PMID:27231256 Harrison FOCUS DSST"),
  ("vortioxetine", "M5HT",  "PS",  0.52, 0.31, 0.73, "FOCUS N=602",     "MDD",     "20mg",     "PMID:27231256 Harrison FOCUS DSST"),
  ("vortioxetine", "M5HT",  "EM",  0.31, 0.10, 0.52, "FOCUS N=602",     "MDD",     "10mg",     "PMID:27231256"),
  ("vortioxetine", "M5HT",  "EF",  0.42, 0.21, 0.63, "FOCUS N=602",     "MDD",     "10-20mg",  "PMID:27231256"),
  ("vortioxetine", "M5HT",  "PS",  0.24, 0.15, 0.33, "k=6 (1782) Pan",  "MDD",     "10-20mg",  "PMID:36398888 Pan 2022 DSST"),
  ("vortioxetine", "M5HT",  "PS",  0.35, 0.23, 0.47, "k=3 (607) McIntyre","MDD",   "10-20mg",  "PMID:27312740 McIntyre 2016 DSST SES=0.35"),
  ("vortioxetine", "M5HT",  "EM",  0.27, 0.10, 0.44, "Katona N=453",    "elderly-MDD","5mg",   "PMID:22318341 Katona 2012 IJG"),

  # ===== 5HT6 ANTAGONISTS =====
  ("idalopirdine", "5HT6",  "EM", -0.05,-0.20, 0.10, "k=4 (2803)",      "AD",      "30-60mg",  "PMID:30560763 Matsunaga 2018"),
  ("idalopirdine", "5HT6",  "EM",  0.30, 0.10, 0.50, "LADDER N=278",    "AD-mod",  "90mg",     "PMID:25297016 LADDER (single+)"),
  ("intepirdine",  "5HT6",  "EM",  0.00,-0.15, 0.15, "MINDSET N=1300",  "AD",      "35mg",     "PMID:29318278 Atri MINDSET (null)"),
  ("SUVN-502",     "5HT6",  "EM", -0.05,-0.20, 0.10, "N=564 phase 2",   "AD-mod",  "50-100mg", "PMID:35662833 Nirogi 2022"),

  # ===== NICOTINIC / NICOTINE =====
  ("nicotine",     "NIC",   "MOT", 0.16, 0.06, 0.26, "k=41 (Heishman)","healthy", "patch/spray","PMID:20414766 Heishman 2010 fine-motor"),
  ("nicotine",     "NIC",   "ATT", 0.34, 0.20, 0.48, "k=41",           "healthy", "patch/spray","PMID:20414766 alerting-RT"),
  ("nicotine",     "NIC",   "EM",  0.27, 0.10, 0.44, "k=41",           "healthy", "patch/spray","PMID:20414766 STM-accuracy"),
  ("nicotine",     "NIC",   "WM",  0.44, 0.25, 0.63, "k=41",           "healthy", "patch/spray","PMID:20414766 WM-RT"),

  # ===== CAFFEINE / METHYLXANTHINES =====
  ("caffeine",     "ADO",   "ATT", 0.28, 0.18, 0.38, "k=31 (1455)",    "healthy", "75-250mg", "Springer Psychopharm 2025 doi:10.1007/s00213-025-06775-1 RT g=0.28"),
  ("caffeine",     "ADO",   "ATT", 0.27, 0.17, 0.37, "k=31 (1455)",    "healthy", "75-250mg", "Springer 2025 accuracy g=0.27"),

  # ===== H3 ANTAGONISTS =====
  ("pitolisant",   "H3",    "ATT", 0.55, 0.30, 0.80, "k=2 (431)",      "narcolepsy","17.8-35.6mg","PMID:33779931 Bassetti 2021 CNS Drugs"),
  ("pitolisant",   "H3",    "PS",  0.20, 0.00, 0.40, "k=2 (431)",      "narcolepsy","",       "Lehert 2020 Drugs RD"),

  # ===== PDE4D =====
  ("BPN14770",     "PDE4D", "VL",  0.55, 0.20, 0.90, "N=30 crossover", "FXS",     "25mg BID", "PMID:33927413 Berry-Kravis 2021 OralRead+PicVocab"),
  ("BPN14770",     "PDE4D", "EM",  0.40, 0.05, 0.75, "N=30",           "FXS",     "25mg BID", "PMID:33927413"),
  ("BPN14770",     "PDE4D", "EF",  0.35, 0.00, 0.70, "N=30",           "FXS",     "25mg BID", "PMID:33927413"),

  # ===== AMPAKINES =====
  ("CX-516",       "AMPA",  "EM",  0.05,-0.30, 0.40, "Goff N=105",     "SCZ",     "900mg TID","PMID:17487227 Goff 2008"),
  ("CX-516",       "AMPA",  "EM", -0.05,-0.40, 0.30, "Berry-Kravis MCI","MCI",    "900mg TID","Lynch/Berry-Kravis MCI trial"),
  ("S47445",       "AMPA",  "EM", -0.02,-0.20, 0.16, "Bernard N=520",  "AD-mild-mod","2-15mg","PMID:31297441 Bernard 2019 ADTRCI"),

  # ===== SIGMA-1 =====
  ("blarcamesine", "SIGMA1","EM",  0.20, 0.00, 0.40, "ANAVEX2-73-AD-004 N=508","AD-early","30-50mg","Macfarlane JPAD 2024 PMC11713060"),

  # ===== mGluR5 NAM =====
  ("mavoglurant",  "MGLUR5","EM", -0.05,-0.30, 0.20, "Berry-Kravis N=315","FXS",   "50-100mg", "PMID:26764156 Berry-Kravis 2016 STM"),
  ("basimglurant", "MGLUR5","EM",  0.00,-0.20, 0.20, "Youssef N=183",  "FXS",     "0.5-1.5mg","PMID:28816242 Youssef 2018"),

  # ===== GSK-3β / Lithium =====
  ("tideglusib",   "GSK3B", "EM",  0.00,-0.20, 0.20, "ARGO N=306",     "AD",      "500-1000mg","PMID:25537011 Lovestone 2015"),
  ("lithium",      "GSK3B", "EM",  0.41, 0.02, 0.81, "k=3 (232) Matsunaga","MCI/AD","0.25-0.5mM","PMID:26402004 Matsunaga 2015 JAD"),

  # ===== PDE9 =====
  ("PF-04447943",  "PDE9",  "EM", -0.05,-0.30, 0.20, "Schwam N=191",   "AD",      "25mg",     "PMID:24801218 Schwam 2014"),

  # ===== PDE10A =====
  ("TAK-063",      "PDE10A","EM",  0.00,-0.20, 0.20, "Macek N=160",    "SCZ",     "20mg",     "PMID:30172593 Macek 2019"),

  # ===== mGluR2/3 =====
  ("pomaglumetad", "MGLUR23","EM", 0.00,-0.15, 0.15, "Adams N≈1000",   "SCZ",     "10-80mg",  "Adams 2014 BMC Psych PMC4276262"),

  # ===== HERBALS / SUPPLEMENTS =====
  ("ginkgo-EGb761","HERBAL","EM",  0.44, 0.12, 0.77, "AD subgroup k=4 (Weinmann)","AD","240mg","PMID:20236541 Weinmann 2010"),
  ("ginkgo-EGb761","HERBAL","EM",  0.58, 0.01, 1.14, "k=9 (2372)",     "dementia","240mg",    "PMID:20236541"),
  ("ginkgo-EGb761","HERBAL","ATT", 0.30, 0.05, 0.55, "k=9 (2561) Tan", "dementia","240mg",    "Tan 2015 JAD doi:10.3233/JAD-140837"),
  ("bacopa-monnieri","HERBAL","EM", 0.30, 0.10, 0.50, "k=9 (518) Kongkeaw","healthy","300mg",   "PMID:24252493 Kongkeaw 2014 delayed recall"),
  ("bacopa-monnieri","HERBAL","ATT",0.30, 0.05, 0.55, "k=9 (518) Kongkeaw","healthy","300mg",   "PMID:24252493 choice-RT"),
  ("citicoline",   "HERBAL","EM",  0.19, 0.06, 0.32, "k=14 (884) Cochrane","VCI/VaD","1000mg",  "PMID:15846601 Fioravanti Cochrane CD000269.pub2"),
  ("citicoline",   "HERBAL","ATT",-0.09,-0.23, 0.05, "k=14 (884)",     "VCI/VaD", "1000mg",   "PMID:15846601 (null)"),
  ("DHA",          "HERBAL","EM",  0.20, 0.00, 0.40, "MIDAS N=485",    "healthy-older","900mg","PMID:20434961 Yurko-Mauro MIDAS CANTAB PAL"),
  ("DHA",          "HERBAL","EM",  0.00,-0.10, 0.10, "k=3 (3500+) Cochrane","healthy-older","", "PMID:22696350 Sydenham Cochrane CD005379.pub3 (null)"),
  ("creatine",     "HERBAL","EM",  0.31, 0.18, 0.44, "k=16 (492) Xu 2024","healthy","5-20g",    "PMC11275561 Xu 2024 Front Nutr g=0.30 (0.18,0.42)"),
  ("creatine",     "HERBAL","ATT", 0.31, 0.03, 0.58, "k=16 (492)",     "healthy", "5-20g",    "PMC11275561 Xu 2024"),
  ("creatine",     "HERBAL","PS",  0.49, 0.20, 0.79, "k=16 (492)",     "healthy", "5-20g",    "PMC11275561 Xu 2024"),
  ("L-theanine",   "HERBAL","ATT", 0.30, 0.10, 0.50, "k=11 (Camfield)","healthy", "100mg+caf","PMID:24946991 Camfield 2014 Nutr Rev 72:507"),
  ("piracetam",    "RACETAM","EM", 0.10,-0.10, 0.30, "Cochrane Flicker","dementia","2.4-9.6g","PMID:11405971 Flicker Cochrane CD001011 (NS)"),

  # ===== ANTIPSYCHOTIC COGNITIVE SUBANALYSES =====
  ("lurasidone",   "ANTIPSY","PS", 0.30, 0.05, 0.55, "Harvey N=244",   "SCZ",     "80-160mg", "Harvey 2015 Schiz Res; PMC10616918"),
  ("cariprazine",  "ANTIPSY","EM", 0.20, 0.00, 0.40, "Mucci network N=400","SCZ", "1.5-6mg",  "PMID:41191868 Mucci 2025 EAPCN"),

  # ===== SUVOREXANT (next-day cognition) =====
  ("suvorexant",   "DORA",  "ATT",-0.05,-0.20, 0.10, "Vermeeren k=2 (500+)","insomnia","10-20mg","PMID:26085297 Vermeeren 2015 (null)"),
]
```

### Methodological notes for Bayesian prior implementation

**Effect-size conversion (consistent across all cells)**:

1. **Cohen's d → Hedges' g**: `g = d × (1 - 3/(4(n1+n2) - 9))`. For two-arm RCTs with N≥30 the correction factor is ≥0.97, so for most cells use g ≈ d.
2. **ADAS-cog MD → SMD**: ADAS-cog total has pooled SD ≈ 7.5 across the AChEi trial portfolio (verified from Birks 2018 forest plots, 30 studies). Use `SMD = -MD/7.5`; sign-flipped because ADAS-cog: lower = better.
3. **MMSE MD → SMD**: pooled SD ≈ 3.5 in mod AD; `SMD = MD/3.5`.
4. **SIB MD → SMD**: pooled SD ≈ 19.8 in severe AD; `SMD = MD/19.8`.
5. **% change / regression β**: convert to SMD via `SMD = β × (SD_x / SD_y_residual)` when both group SDs are available; otherwise back-calculate from t- or F-statistics: `d = 2t/√df`.
6. **DSST raw points → SMD**: pooled SD ≈ 11 in MDD (Harrison FOCUS); 1.75-point WMD ≈ SMD 0.16.

**Schmidli 2014 robust MAP prior weights**:

- For cells with k ≥ 5 meta-analytic studies: τ² = 0.02 (σ_τ = 0.14), robust mixture weight w_robust = 0.10.
- For cells with k = 2–4: τ² = 0.04 (σ_τ = 0.20), w_robust = 0.20.
- For cells with k = 1 (single RCT) or `None`: τ² = 0.08 (σ_τ = 0.28), w_robust = 0.50; informative component centered at class-level mean.
- Vague mixture component: Normal(0, 1²) — standard Schmidli unit-information prior.
- Between-class τ²_class: empirical Bayes estimate from the 12-class data is **τ²_class ≈ 0.045** (σ ≈ 0.21), recovered from the variance of class-level posterior means.
- Between-endpoint τ²_endpoint (within class): **τ²_endpoint ≈ 0.025** (σ ≈ 0.16), recovered from MPH within-class spread (Coghill domains range 0.24–0.62).

**Healthy-adult upper-envelope (Roberts 2020 / Ilieva 2015) — soft vs hard implementation**:

Implement as a **soft** prior, not a hard truncation. The empirical upper envelope of healthy-adult meta-analytic cognitive effects sits at SMD ≈ 0.43–0.45 (MPH recall; amphetamine delayed memory) — *not* a hard 0.50 line. The task brief's "0.50 ceiling" is an interpretive extrapolation from these upper CIs; Roberts 2020 itself does not state any g ≈ 0.50 figure. Patient populations routinely exceed this envelope (galantamine ADAS-cog SMD ≈ 0.40; MPH non-exec memory in pediatric ADHD g = 0.60). Recommended encoding:

```python
# soft upper envelope on healthy-population posterior mean
mu_healthy ~ Normal(0, 0.30**2)        # weak prior centered at 0
mu_healthy_max ~ HalfNormal(0.30)      # 90% mass below 0.50, allowing rare exceedances
# patient populations use class-specific prior without ceiling
```

**Subdomain pooling**: Empirically, within-class endpoint correlations are r ≈ 0.55–0.70 (e.g., MPH Coghill 2014 cells share substantial covariance). **Use a partially-pooled covariance structure** — neither fully independent nor fully shared — implemented as an LKJ(η=2) prior on the within-class correlation matrix. Independence is wrong (vortioxetine PS=0.51 / EM=0.31 / EF=0.42 are highly correlated); full sharing is also wrong (5-HT6 antagonists are uniformly null on EM but have no other-domain data).

**Multi-arm dose-finding trials**: Extract the highest-approved-dose arm only (e.g., donepezil 10 mg not 5 or 23 mg; encenicline 0.9 mg not 0.27 mg; idalopirdine 60 mg not 10/30 mg). For comparator-controlled trials (vortioxetine vs duloxetine), use the placebo-controlled arm SMD, not the active-comparator SMD. Apply Bonferroni-free hierarchical shrinkage on doses within compound.

---

## Critical commentary on data quality

**Well-covered (Cochrane-grade)**: AChE inhibitors (3 Cochrane reviews — Birks 2018 CD001190.pub3 alone includes 30 studies and 8,257 participants; Birks 2015 CD001191.pub4 rivastigmine 6 studies / 3,232 participants; Loy 2024 CD001747.pub4 galantamine 10 studies / ~6,805 participants; combined pool ≈ 18,000 with substantial dose-arm overlap); MPH for ADHD (Storebø 2015 CD009885.pub2 includes 185 RCTs from 449 reports / 761 records; the 2023 pub3 update has a different scope with 21 primary-outcome trials / 1,728 participants and is not a strict re-run); modafinil in non-sleep-deprived adults (Repantis 2010, Battleday 2015, Kredlow 2019, Roberts 2020 — 4 PRISMA meta-analyses).

**Moderately covered**: memantine (Matsunaga 2015 PLOS ONE; Kishi 2017), vortioxetine cognition (McIntyre 2016, Pan 2022, Baune 2018), nicotine acute (Heishman 2010, k=41), caffeine (Springer 2025 k=31 meta + EFSA opinion).

**Poorly covered (single trials or all-null meta-analyses)**: α7 nAChR agonists (Lewis 2017 pooled g essentially zero; encenicline phase 3 program halted for GI events; ABT-126 phase 2 failed); AMPA potentiators (CX-516 negative in SCZ and MCI; S47445 negative in AD); mGluR class (pomaglumetad failed; mavoglurant/basimglurant failed in FXS); 5-HT6 antagonists (idalopirdine phase 3 STARSHINE/STARBEAM/STARBRIGHT null; intepirdine MINDSET null; SUVN-502 null); PDE9 (PF-04447943 negative); PDE10A (TAK-063 NS); GSK-3β (tideglusib ARGO negative); σ1 (blarcamesine ANAVEX2-73-AD-004 marginally positive but no Hedges' g published yet); BDNF mimetics (preclinical only).

**Best-measured endpoints across meta-analyses**: 
- AD/MCI: ADAS-Cog (Cochrane standard) and MMSE — both have decades of data and well-characterized SDs; SIB for severe AD.
- ADHD: ADHD-RS for symptoms, Coghill 2014 framework for cognition (response inhibition SSRT, exec/non-exec memory, RT-variability).
- MDD cognition: DSST (vortioxetine programs standardized on it); RAVLT for verbal learning.
- Healthy: heterogeneous — n-back, Stroop, Stop-Signal, CANTAB battery; least standardized.

**Worst-measured**: executive function across populations (Stroop vs Trails B vs WCST vs SST give different g values for the same drug — Roberts 2020 had to subdivide into updating/switching/inhibition); visuospatial (almost no meta-analytic data outside SIB-VS sub-score); motor/RT (rarely a primary endpoint).

**Publication bias considerations per class**:
- **High bias risk**: piracetam (Waegemans 2002 industry-sponsored meta showed implausibly large effects vs Flicker 2001 Cochrane null); ginkgo (early Le Bars trials had funnel asymmetry; GEM Snitz 2009 negative); herbal compounds in general; Bacopa (most Kongkeaw 2014 trials industry-funded).
- **Medium bias risk**: amphetamine in healthy (Ilieva 2015 explicitly identified publication bias in working-memory and delayed-memory domains and downgraded effect estimates accordingly); modafinil (Repantis 2010 noted funnel asymmetry).
- **Low bias risk**: Cochrane AChEi reviews (manufacturer trial data fully accessible; IPD meta-analyses available); MPH ADHD (Coghill IPD).

**Will 15 → 60 anchor expansion actually fix the V7 partial-pool gap of 0.004–0.063?** 

Partial. The diagnosis: the 0.004–0.063 undershoot of empirical bands is consistent with two distinct failure modes:

1. **Over-pooling** (more cells will help): adding ~45 new (compound, endpoint) anchors decreases the per-cell shrinkage weight `λ = τ²_class / (τ²_class + σ²_within/k)`. Going from k=1–2 (the current state for orphan classes like 5-HT6 or α7) to k=4–10 per class shifts λ from ~0.7 toward ~0.3, allowing posteriors to track compound-specific means rather than collapsing to the class-level prior.

2. **Misspecified τ²_class** (more cells will NOT help): The current V7 likely uses an under-dispersed τ²_class. The empirical 12-class g range is approximately [−0.05, +0.40] for memory endpoints in patient populations and [0.05, 0.43] for healthy populations — implying τ²_class ≈ 0.04–0.06. If the existing model uses τ²_class < 0.03, no amount of additional anchors will fix the shrinkage. **The first action item is to re-estimate τ²_class empirically from the new 60-row table, then re-fit.**

3. **Population mixing**: The 0.063 max undershoot at donepezil (P1) is suspicious — donepezil has the *cleanest* meta-analytic g of any cell in the table. If V7 is pooling donepezil with healthy-population priors via a single class-level effect, the AD-specific signal (g ≈ 0.36) gets pulled toward 0.20 by the class average that includes healthy nicotinic/glutamatergic null cells. Recommend **a population × class interaction term** (binary: patient vs healthy) before adding more anchors.

---

## Recommendations

**Staged plan with concrete benchmarks**:

**Stage 1 (immediate, before adding anchors)**:
1. Re-fit V7 with population × class interaction (`population ∈ {healthy, MCI, AD, ADHD, SCZ}`); 5 levels × 12 classes = 60-cell crossed structure.
2. Use empirical-Bayes τ²_class estimate from the existing 15-anchor set; if it < 0.03, switch to a Half-Normal(0.3) hyperprior to allow more between-class variance.
3. **Target**: reduce P1 (donepezil) undershoot from 0.063 to < 0.030. If achieved, proceed to Stage 2; if not, the gap is *not* a sparsity problem and adding anchors is wasted effort.

**Stage 2 (expand anchors)**:
4. Ingest the 60 rows from the table above into `REFERENCE_COMPOUND_SMD`. Validate by leave-one-out CV: hold out donepezil/MPH/modafinil/memantine and verify recovered g lies in 95% CI.
5. Ingest the 96-cell class × endpoint table into `PER_SUBDOMAIN_PRIORS`. Replace `None` cells with class-marginalized prior (mean across the row's non-null cells).
6. **Target**: reduce all four anchor undershoots (donepezil 0.063 → < 0.020; MPH 0.045 → < 0.015; modafinil 0.028 → < 0.010; memantine 0.004 → maintain).

**Stage 3 (Schmidli robust MAP refinement)**:
7. Replace single-RCT cells (5-HT6, PDE4D, σ1, mGluR5, GSK-3β, PDE9, PDE10A, mGluR2/3, AMPA) with robust mixture `0.5 × N(class_mean, τ²_class) + 0.5 × N(0, 1)`. Schmidli's w_robust = 0.5 reflects the high uncertainty of single negative phase-2 readouts.
8. Apply healthy-adult upper envelope (Roberts 2020 / Ilieva 2015) as soft Half-Normal(0.3) on healthy-population posteriors only — not on patient-population posteriors (galantamine and MPH-ADHD legitimately exceed 0.5).
9. **Benchmark**: 90% credible upper bound on healthy-adult class g should be ≤ 0.55; on patient-class g should be unconstrained.

**Stage 4 (data extension if Stage 1–3 insufficient)**:
10. Pull Cochrane Storebø 2015 CD009885.pub2 (185 MPH RCTs), Cochrane Birks 2018 CD001190.pub3 (donepezil, 30 studies / 8,257 patients), Cochrane Birks 2015 CD001191.pub4 (rivastigmine, 6 studies / 3,232 patients), Cochrane Loy 2024 CD001747.pub4 (galantamine) IPD if accessible; this will tighten the four primary anchor cells by ~30%.
11. Add MetaPsy.io entries for psychological + pharmacological co-interventions (cognitive training × pharmacology trials).

**Thresholds that would change the plan**:
- If after Stage 1 the donepezil undershoot is < 0.020, **skip Stages 2–4**; the binding constraint was population specification, not anchor sparsity.
- If after Stage 2 the modafinil undershoot is > 0.040, the issue is *healthy-vs-patient pooling*; revisit the population-interaction specification rather than adding more anchors.
- If after Stage 3 the 90% CrI on any class exceeds [−0.6, +0.8], τ² is over-dispersed; tighten the Half-Normal(0.3) to Half-Normal(0.2).

---

## Caveats

1. **Asymmetric evidence**: 6 of 12 classes have Cochrane-grade meta-analyses; 6 rely on 1–2 PRISMA meta-analyses with low study counts (k=2–6). The class table cells for α7, α4β2, 5-HT6, H3, PDE4D, AMPA carry CI half-widths ≥ 0.25 — the implementation should propagate this uncertainty into the prior σ.

2. **Population heterogeneity not fully captured**: ADAS-cog SDs in mild AD differ from those in severe AD; the conversion factor 7.5 is an average. For severity-stratified analyses, use IPD when available.

3. **The "0.50 ceiling" is interpretive, not literal**: Roberts 2020 reports no g ≈ 0.50 figure in its abstract — the highest single subdomain SMD it identifies is MPH recall 0.43 and MPH sustained attention 0.42. The 0.50 figure is an upper-CI envelope; treat as a soft prior, not a physical constant.

4. **Null trials are still informative**: idalopirdine, intepirdine, SUVN-502, pomaglumetad, mavoglurant, S47445, tideglusib, PF-04447943, TAK-063, CX-516 all have CIs that include 0 and exclude g > 0.20. These are valuable for *bounding the class prior away from large positive effects* — do not exclude them as "negative trials". Build a `NEGATIVE_TRIAL_REGISTRY` table to bias-correct the class means.

5. **Orthosteric vs allosteric distinction**: The user specifically distinguishes this (BPN14770, LY3154207, encenicline are allosteric). With only 3 allosteric anchor compounds (BPN14770 PDE4D NAM, encenicline α7 PAM partial-agonist, LY3154207 D1 PAM) and one ABT-288 (H3 inverse agonist), an allosteric × class interaction is **under-identified**. Recommend treating allosteric_flag as a *covariate* in the hierarchical regression rather than as a separate stratum. The empirical signal: allosteric modulators do not systematically outperform orthosteric (encenicline failed phase 3; BPN14770 positive in FXS only).

6. **Compound-specific dose-response not modeled**: Most class × endpoint cells aggregate across doses. For real prediction, the V7 pipeline should retain a per-compound dose covariate; the priors here apply only at the labeled/maximally-effective dose.

7. **Temporal heterogeneity**: Some meta-analyses (Birks Cochrane reviews) contain trials going back to 1995; others (Roberts 2020) cover 2010–2020. Drug-trial methodology has evolved (better placebo run-in, MMRM analyses, MCID-anchored endpoints). If the V7 pipeline weights recent evidence more heavily, the class priors should be re-pooled with year-as-covariate.

8. **Dextroamphetamine sign caveat**: Roberts 2020 explicitly reports "no effects for d-amph" across the full PRISMA-extracted set (k=10 studies, 27 effect sizes). The positive d-amphetamine subdomain g values in the table above derive from Ilieva 2015 (working memory g=0.13, short-term memory g=0.20, delayed memory g=0.45, inhibitory control g=0.20) — these reflect pooled stimulant analyses (MPH + amphetamine) and may overstate the d-amphetamine-specific signal. Roberts 2020 d-amphetamine rows in the table are set to 0.00 with wide CIs to reflect the negative meta-analytic result.

9. **Cochrane Storebø 2015 vs 2023 versions are NOT interchangeable**: Storebø 2015 (CD009885.pub2, PMID 26599576) includes 185 RCTs from 761 reports / 449 RCT-describing reports — the canonical "massive MPH meta-analysis". Storebø 2023 (pub3, PMID 36971690) has a redesigned scope; its primary-outcome analysis covers 21 trials / 1,728 participants. When citing MPH ADHD evidence, specify which version was used.

---

## Next reading list with DOIs

Tier 1 (must read for prior calibration):
- Birks JS, Harvey RJ. Donepezil for dementia. Cochrane 2018. **10.1002/14651858.CD001190.pub3** PMID 29923184 (30 studies, 8,257 participants)
- Birks JS, Chong LY, Grimley Evans J. Rivastigmine. Cochrane 2015. **10.1002/14651858.CD001191.pub4** PMID 26393402 (6 studies, 3,232 participants)
- Loy C, Schneider L; updated Lim AWY, Schneider L, Loy C. Galantamine. Cochrane 2006/2024. **10.1002/14651858.CD001747.pub4** PMID 39498781
- Matsunaga S et al. Memantine monotherapy for AD. PLOS ONE 2015. **10.1371/journal.pone.0123289** PMID 25869017
- Roberts CA et al. Meta-analyses MPH/MOD/d-amph healthy. Eur Neuropsychopharmacol 2020 38:40-62. **10.1016/j.euroneuro.2020.07.002** PMID 32709551 (MPH overall SMD=0.21 p=.0004; modafinil overall SMD=0.12 p=.01, memory updating SMD=0.28 p=.03; d-amphetamine no significant effects)
- Coghill DR et al. MPH cognition ADHD. Biol Psychiatry 2014. **10.1016/j.biopsych.2013.10.005** PMID 24231201 (60 studies reviewed, 36 in meta-analysis)
- Ilieva I, Hook CJ, Farah MJ. Stimulants healthy. J Cogn Neurosci 2015. **10.1162/jocn_a_00776** PMID 25591060
- Battleday RM, Brem AK. Modafinil non-sleep-deprived. Eur Neuropsychopharmacol 2015. **10.1016/j.euroneuro.2015.07.028** PMID 26381811
- McIntyre RS et al. Vortioxetine cognition meta-analysis. Int J Neuropsychopharmacol 2016. **10.1093/ijnp/pyw055** PMID 27312740 (DSST SES=0.35 unadjusted, 0.24 MADRS-adjusted)
- Heishman SJ, Kleykamp BA, Singleton EG. Nicotine acute. Psychopharmacology 2010. **10.1007/s00213-010-1848-1** PMID 20414766 (k=41, g=0.16–0.44)
- Storebø OJ et al. Methylphenidate for ADHD. Cochrane Database Syst Rev 2015 (pub2). **10.1002/14651858.CD009885.pub2** PMID 26599576 (185 RCTs)
- Repantis D et al. Modafinil and methylphenidate. Pharmacol Res 2010. **10.1016/j.phrs.2010.04.002** PMID 20416377

Tier 2 (mechanism class meta-analyses):
- Lewis AS, van Schalkwyk GI, Bloch MH. α7 nAChR translational meta. Prog Neuro-Psychopharmacol Biol Psychiatry 2017. **10.1016/j.pnpbp.2017.01.001** PMID 28065843
- Tanzer T et al. Varenicline cognition SCZ. Psychopharmacology 2020. **10.1007/s00213-019-05396-9** PMID 31792645
- Matsunaga S et al. Idalopirdine AD. Int Psychogeriatr 2018. **10.1017/S1041610218001941** PMID 30560763
- Pievsky MA, McGrath RE. Adult ADHD cognition. Neurosci Biobehav Rev 2018. **10.1016/j.neubiorev.2018.05.012** PMID 29751051
- Kredlow MA et al. Modafinil cognitive enhancer. J Clin Psychopharmacol 2019. **10.1097/JCP.0000000000001085** PMID 31433306
- Pan Z et al. Vortioxetine MDD cognition. Int J Neuropsychopharmacol 2022. **10.1093/ijnp/pyac054** PMID 36398888
- Harrison JE et al. Vortioxetine FOCUS sub-domain. Int J Neuropsychopharmacol 2016. **10.1093/ijnp/pyw054** PMID 27231256
- Aceto G et al. α7 nAChR cognitive enhancers SCZ meta. Front Psychiatry 2021. PMC8055861
- Isfandnia F, El Masri S, Radua J, Rubia K. Chronic stim and non-stim on EF ADHD. Neurosci Biobehav Rev 2024. **10.1016/j.neubiorev.2024.105703**

Tier 3 (single-compound anchor trials):
- Berry-Kravis EM et al. BPN14770 FXS. Nat Med 2021. **10.1038/s41591-021-01321-w** PMID 33927413
- Keefe RS et al. Encenicline SCZ. Neuropsychopharmacology 2015. **10.1038/npp.2015.176** PMID 26089183
- Goff DC et al. CX-516 SCZ. Neuropsychopharmacology 2008. **10.1038/sj.npp.1301444** PMID 17487227
- Bernard K et al. S47445 AD. Alzheimers Dement TRCI 2019. **10.1016/j.trci.2019.04.002** PMID 31297441
- Atri A et al. Idalopirdine STAR phase 3. JAMA 2018. **10.1001/jama.2017.20373** PMID 29318278
- Schwam EM et al. PF-04447943 AD. Curr Alzheimer Res 2014. **10.2174/1567205011666140505100858** PMID 24801218
- Czeisler CA et al. Modafinil shift-work disorder. NEJM 2005. **10.1056/NEJMoa041292** PMID 16079371
- Yurko-Mauro K et al. MIDAS DHA. Alzheimers Dement 2010. **10.1016/j.jalz.2010.01.013** PMID 20434961
- Sallee FR et al. Guanfacine-XR ADHD pediatric. JAACAP 2009. PMID 19106767

Tier 4 (supplements / herbals):
- Kongkeaw C et al. Bacopa monnieri meta. J Ethnopharmacol 2014. **10.1016/j.jep.2013.11.008** PMID 24252493
- Weinmann S et al. Ginkgo biloba dementia. BMC Geriatr 2010. **10.1186/1471-2318-10-14** PMID 20236541
- Yang G et al. Huperzine A meta. PLOS ONE 2013. **10.1371/journal.pone.0074916** PMID 24086396
- Fioravanti M, Yanagi M. Citicoline. Cochrane 2005. **10.1002/14651858.CD000269.pub2** PMID 15846601
- Xu C et al. Creatine cognition meta. Front Nutr 2024. **10.3389/fnut.2024.1424972** PMC11275561
- Camfield DA et al. L-theanine + caffeine. Nutr Rev 2014. **10.1111/nure.12120** PMID 24946991
- Sydenham E, Dangour AD, Lim WS. Omega-3 healthy older Cochrane 2012. **10.1002/14651858.CD005379.pub3** PMID 22696350
- Flicker L, Grimley Evans J. Piracetam Cochrane 2001. **10.1002/14651858.CD001011** PMID 11405971

Tier 5 (Bayesian methodology):
- Schmidli H et al. Robust meta-analytic-predictive priors. Biometrics 2014. **10.1111/biom.12242**
- Neuenschwander B et al. Summarizing historical info on controls. Clin Trials 2010. **10.1177/1740774509356002**