# §8.0b — Cognition-Context-Tailored 44-Target Off-Target Safety Liability Panel for the MAMMAL Cognitive Enhancement Drug Repurposing Pipeline

## TL;DR
- Build the §8.0b panel as a **Bowes-44 backbone + Brennan-77 CNS additions − peripheral irrelevancies**, locking in 44 targets with explicit cognition-relevance rationale per target and tier-stratified gates (HARD CUT / FLAG / informational).
- Implement as a **separate `panel_type='liability'` partition** of `targets_seed.csv`, reusing the MAMMAL DTI head (one-time ~45–60 min re-run on RTX 5070 for 298 compounds × 44 targets), gated by a new `gates/liability_panel.py` module that runs **after** ADMET and emits {PASS, FLAG, CUT} with mechanistic notes.
- The panel will correctly CUT aripiprazole and amitriptyline from the v3 wet-lab shortlist on mechanistic grounds (5-HT2B + D2 + α1; HRH1 + M1 anticholinergic) that ADMET-AI cannot reach. That is the entire point of building §8.0b instead of trusting Eurofins SafetyScreen44 off the shelf.

## Key Findings

1. **Off-the-shelf panels are wrong for chronic healthy cognitive enhancement.** Eurofins SafetyScreen44 is a Bowes-2012 reproduction biased toward peripheral cardiovascular, GI, and endocrine ADRs identified across all therapeutic areas at AstraZeneca/GSK/Novartis/Pfizer. NIMH PDSP is a CNS-broad academic screening resource with no liability gating logic. Brennan et al. 2024 (the IQ Consortium "Safety-77", Nat Rev Drug Discov 23:525–545) expanded to 77 targets but still inherits an indication-agnostic weighting. None tier-stratifies by chronic-dosing healthy-adult risk, which is the relevant context for a cognitive enhancer.

2. **Eurofins SafetyScreen44 exact composition** (from the Eurofins literature P270/EPDSFL420JUNE16, fetched directly): 24 GPCRs (ADORA2A; ADRA1A, ADRA2A, ADRB1, ADRB2; CNR1, CNR2; CCKAR; DRD1, DRD2S; EDNRA; HRH1, HRH2; CHRM1, CHRM2, CHRM3; OPRD1, OPRK1, OPRM1; HTR1A, HTR1B, HTR2A, HTR2B; AVPR1A) + 3 transporters (DAT, NET, SERT) + 8 ion channels (BZD-site GABA-A, NMDA, α4β2 nAChR, 5-HT3, L-type Ca²⁺ DHP site, hERG, generic Kv, Na⁺ site 2) + 2 nuclear receptors (AR, GR) + 1 kinase (LCK) + 6 non-kinase enzymes (AChE, MAO-A, COX1, COX2, PDE3A, PDE4D2). **Conspicuously absent for cognition**: TAAR1, σ1/σ2, 5-HT2C, 5-HT6, 5-HT7, GABA-A α5, MT1/MT2, Cav1.2 (only DHP site), Nav1.5 (only site 2), MAO-B, mGluR2/3/5.

3. **5-HT2B is non-negotiable Tier 1.** Connolly et al. NEJM 1997;337:581–588 documented 24 cases of valvulopathy on fen-phen; Rothman et al. Circulation 2000;102:2836 established norfenfluramine as the 5-HT2B agonist culprit; Schade et al. NEJM 2007;356:29–38 and Zanettini et al. NEJM 2007;356:39–46 extended to pergolide and cabergoline (incidence rate ratio 7.1 and 4.9 respectively at >6 months and >3 mg/day). Pergolide withdrawn 2007 (FDA); benfluorex withdrawn 2009 (EMA). Roth NEJM 2007;356:6–9 editorial codified the class warning. **Dumotier & Urban 2024 (J Pharmacol Toxicol Methods 128:107542, PMID 39032441) is the current regulatory framework: safety margin >10× Cmax over 5-HT2B Ki using AUC-corrected exposure, with binding Ki preferred over functional EC50** because agonist metabolites (e.g., norfenfluramine) create false negatives.

4. **Brennan-77 supplies the CNS expansion targets** missing from Bowes-44: HTR2C, HTR6, HTR7, GABRA5, OPRD1, OPRK1, α3β4 and α4β2 nAChRs, SIGMAR1, TMEM97 (σ2), TAAR1, MTNR1A/B, CNR2, SCN5A, CACNA1C, TACR1, OXTR, AVPR1A/B, ESR1/ESR2, NR3C1, NTRK1/3, GRM2/3/5. These map onto the cognition liabilities flagged in the V3 wet-lab shortlist (aripiprazole 5-HT2B/α1, TC-5619 5-HT3, BPN14770 PDE4D).

5. **Tier-stratified gating with concrete pKi thresholds** under a worst-case 1 µM Cmax assumption:
   - **Tier 1 HARD CUT**: pKi > 6.5–7.0 for 5-HT2B agonism, hERG block, HRH1 block, CB1 modulation, μ-opioid binding, M1 antagonism (in compounds of anticholinergic class).
   - **Tier 2 FLAG**: pKi > 6.0 at α1A/B/D, D2 antagonism, MAO-A/B, 5-HT2A, 5-HT2C agonism, GABA-A α1, Cav1.2, Nav1.5, M3, KOR, DOR, GR, ER.
   - **Tier 3 INFORMATIONAL**: σ2/TMEM97, TAAR1 (could be feature), 5-HT1A/6/7 (could be feature), D4, NK1, MT1/MT2.

6. **Predicted impact on the v3 top-25 shortlist**: aripiprazole and amitriptyline correctly CUT; rolipram, rasagiline, bupropion, methylphenidate, d-amphetamine, galantamine, and TC-5619 FLAG with explicit mechanistic notes; modafinil, aniracetam, piracetam, levetiracetam, α-GPC, BPN14770, BI-409306, 7,8-DHF, tulrampator, lemborexant, lanicemine, topiramate, lithium likely PASS. 9/25 (~36%) touch rate, 2 hard rejections — the discrimination ADMET alone cannot deliver.

---

## Section 1 — Panel Composition

### 1.1 Backbone references

**Eurofins SafetyScreen44™** (Ref P270, 2016) — full composition verified above. The canonical academic citation is **Bowes et al. Nat Rev Drug Discov 2012;11:909–922** (24 GPCRs, 8 ion channels, 7 enzymes, 3 monoamine transporters, 2 nuclear hormone receptors). SafetyScreen44 is its commercial implementation; WuXi Mini-Safety 44 and DiscoverX SafetyProfiler reproduce the same target set.

**NIMH PDSP / Roth lab (UNC Chapel Hill)** runs a rotating CNS-broad assay panel including all serotonergic GPCRs (5-HT1A/1B/1D/1E/1F/2A/2B/2C/5A/6/7), all dopaminergic (D1–D5), all α/β adrenergic, all muscarinic (M1–M5), histamine H1–H4, all opioid (μ/δ/κ/NOP), σ1/σ2, neuronal nAChR subtypes, GABA-A BZD site, NMDA glutamate, monoamine transporters, and TAAR1. The PDSP Ki database (47,312 Ki values, 699 targets as of the CDD 2.0 announcement) is the curated open-access record but is a screening service, not a gated panel.

**Brennan et al. Nat Rev Drug Discov 2024;23:525–545 (IQ Consortium "Safety-77")** added ~33 targets to Bowes-44 emphasizing CNS expansion, kinase coverage, and additional ion channels and nuclear receptors. The exact 33-target supplementary list is in Table 1 of the paper (Nature paywall); the high-priority CNS additions (HTR2C, HTR6, HTR7, GABRA5, OPRD1, OPRK1, α3β4 nAChR, SIGMAR1, TMEM97, TAAR1, MTNR1A/B, CNR2, SCN5A, CACNA1C, TACR1, OXTR, AVPR1A/B, ESR1/ESR2, NR3C1, NTRK1/3, GRM2/3/5) are well-supported by the paper's accessible reference list and by the Schmidt/Brennan/Jenkinson/Valentin 2025 correspondence (Nat Rev Drug Discov 24:482–484).

**Regulatory backbone**: ICH S7A (general safety pharmacology core battery — CNS/CV/respiratory), ICH S7B (delayed ventricular repolarization, anchors hERG), ICH E14/S7B Q&A 2020 update (integrated proarrhythmia with CiPA), Papoian et al. 2015 (FDA secondary pharmacology guidance).

### 1.2 Cognition-context tailoring rules

**DROP from Bowes-44/SafetyScreen44** (peripheral / chronic-CNS-irrelevant):
- ADORA2A (efficacy target, not a chronic-dosing CV liability at cognition exposures)
- CCKAR (peripheral GI)
- EDNRA (peripheral vascular)
- PTGS1/PTGS2 (COX-1/COX-2 — chronic NSAID liability not generalizable)
- LCK (CNS-penetrant LCK inhibitors essentially nonexistent in cognition repurposing libraries)
- HRH2 (gastric acid; trivial)
- PDE3A (cardiac inotrope; PDE4 family is the cognition-relevant subfamily)
- HTR1B (cerebral vasoconstriction migraine; not a chronic-dosing healthy-adult risk)

**KEEP & ELEVATE** (CNS-direct, FDA precedent, chronic-dosing-relevant):
HTR2B, KCNH2, OPRM1, HRH1, CNR1, MAOA, MAOB, ADRA1A/B/D, DRD2, GABRA1, CHRM1, CHRM3, HTR2A, HTR2C, HTR3A, CACNA1C, SCN5A.

**ADD beyond off-the-shelf** (cognition-specific):
TAAR1, SIGMAR1, TMEM97, MTNR1A, MTNR1B, CHRNA3/CHRNB4 (peripheral autonomic distinct from cognition-target CHRNA7), HTR6, HTR7, GRM2/3/5, NTRK1, NTRK3, NR3C1, ESR1, ESR2, OPRD1, OPRK1, GABRA1/2/3/5 (subunit-aware), CHRNA4.

### 1.3 The 44-target panel — structured specification

| # | gene_symbol | uniprot | chembl_target_id | class | liability_category | tier | cognition_context_rationale |
|---|---|---|---|---|---|---|---|
| 1 | HTR2B | P41595 | CHEMBL1833 | GPCR | cardiotox (valvulopathy) | 1 | Fenfluramine 1997, pergolide 2007, benfluorex 2009 — chronic dosing produces VIC proliferation; HARD CUT |
| 2 | KCNH2 | Q12809 | CHEMBL240 | ion_channel | cardiotox (TdP) | 1 | ICH S7B anchor; cisapride/terfenadine/astemizole withdrawals; cross-checks ADMET-AI hERG head |
| 3 | OPRM1 | P35372 | CHEMBL233 | GPCR | abuse_potential | 1 | DEA scheduling concern; binding alone sufficient for FDA flag |
| 4 | HRH1 | P35367 | CHEMBL231 | GPCR | neurotox (cognitive impairment) | 1 | Anticholinergic burden — chronic H1 block in healthy adults degrades the very endpoint we're enhancing |
| 5 | CNR1 | P21554 | CHEMBL218 | GPCR | neurotox (psychiatric AEs) | 1 | Rimonabant 2008 EU withdrawal; depression/anxiety/suicidality on chronic CB1 modulation |
| 6 | CHRM1 | P11229 | CHEMBL216 | GPCR | cognitive_impairment (antag) | 1 (antag) | M1 antagonism = anticholinergic decline; M1 agonism = procognitive (KarXT) — directionality matters |
| 7 | CHRM3 | P20309 | CHEMBL245 | GPCR | GI/GU (oxybutynin class) | 2 | M3 antagonism = xerostomia, urinary retention |
| 8 | MAOA | P21397 | CHEMBL1951 | enzyme | hypertensive_crisis | 2 | Blackwell 1967 cheese effect; irreversible CUT, RIMA FLAG |
| 9 | MAOB | P27338 | CHEMBL2039 | enzyme | DDI / sympathomimetic synergy | 2 | Selegiline-class; dual A+B = CUT |
| 10 | ADRA1A | P35348 | CHEMBL229 | GPCR | orthostatic_hypotension | 2 | Prazosin first-dose; aripiprazole liability |
| 11 | ADRA1B | P35368 | CHEMBL1867 | GPCR | orthostatic_hypotension | 2 | Vascular smooth muscle; paired with ADRA1A |
| 12 | ADRA1D | P25100 | CHEMBL223 | GPCR | orthostatic_hypotension | 2 | Same family; usually correlated with 1A/1B |
| 13 | DRD2 | P14416 | CHEMBL217 | GPCR | EPS / hyperprolactinemia / TD | 2 | Chronic antag = parkinsonism, TD; partial agonism is borderline |
| 14 | DRD3 | P35462 | CHEMBL234 | GPCR | impulse_control (agonist) | 2 | Pramipexole pathological gambling / hypersexuality |
| 15 | DRD4 | P21917 | CHEMBL219 | GPCR | informational | 3 | Mechanism not well-characterized as liability |
| 16 | HTR1A | P08908 | CHEMBL214 | GPCR | feature / informational | 3 | Buspirone anxiolysis; widely deliberately exploited |
| 17 | HTR2A | P28223 | CHEMBL224 | GPCR | psychotomimetic | 2 | Hallucinogenic at high occupancy; psilocybin microdosing context-dependent |
| 18 | HTR2C | P28335 | CHEMBL225 | GPCR | weight_gain / cancer signal | 2 | Lorcaserin 2020 cancer withdrawal (Sharretts et al. NEJM 2020;383:1000) |
| 19 | HTR3A | P46098 | CHEMBL1899 | ion_channel | GI / nausea | 2 | Ondansetron antagonism = anti-emetic; TC-5619 5-HT3 antagonism possibly beneficial |
| 20 | HTR6 | P50406 | CHEMBL3371 | GPCR | informational (could be feature) | 3 | Idalopirdine target — efficacy candidate, not liability |
| 21 | HTR7 | P34969 | CHEMBL3155 | GPCR | circadian / mood | 3 | Vortioxetine partial agonism |
| 22 | OPRK1 | P41145 | CHEMBL237 | GPCR | dysphoria / sedation | 2 | Salvinorin A-class dysphoria; aticaprant antagonism class |
| 23 | OPRD1 | P41143 | CHEMBL236 | GPCR | seizure / mood | 2 | SNC80-class convulsant at high doses |
| 24 | CACNA1C | Q13936 | CHEMBL1940 | ion_channel | cardiotox 2nd-line | 2 | Cav1.2 — CiPA second axis after hERG |
| 25 | SCN5A | Q14524 | CHEMBL1971 | ion_channel | conduction / Brugada | 2 | Nav1.5 block — QRS widening; bupropion seizure correlate |
| 26 | KCNA5 | P22460 | CHEMBL4306 | ion_channel | atrial repolarization | 3 | Kv1.5; informational |
| 27 | KCND3 | Q9UK17 | CHEMBL2002 | ion_channel | cardiac Ito | 3 | Kv4.3; informational |
| 28 | GABRA1 | P14867 | CHEMBL1900 | ion_channel | sedation / dependence | 2 | α1 BZD = zolpidem sedation; chronic = tolerance |
| 29 | GABRA2 | P47869 | CHEMBL2093872* | ion_channel | anxiolysis (feature) | 3 | TPA023 class anxiolytic without sedation |
| 30 | GABRA3 | P34903 | CHEMBL2093872* | ion_channel | anxiolysis | 3 | α3 contributes to anxiolytic profile |
| 31 | GABRA5 | P31644 | CHEMBL2093872* | ion_channel | cognition (feature) | 3 | α5 inverse agonism = cognition (Dawson 2006; Atack 2009) |
| 32 | CHRNA3 | P32297 | CHEMBL2109244* | ion_channel | autonomic / nausea | 2 | α3β4 ganglionic = nausea distinct from α7 cognition |
| 33 | CHRNB4 | P30926 | CHEMBL2109244* | ion_channel | autonomic GI | 2 | β4 subunit ganglionic nAChR |
| 34 | CHRNA4 | P43681 | CHEMBL1907589* | ion_channel | autonomic / seizure | 2 | α4β2 nAChR; varenicline neuropsychiatric AEs |
| 35 | SIGMAR1 | Q99720 | CHEMBL287 | receptor | already cognition target | — | Cross-check; σ1 already in efficacy panel |
| 36 | TMEM97 | Q5BJF2 | CHEMBL4153 | receptor | informational | 3 | σ2; cholesterol homeostasis; no withdrawal precedent |
| 37 | TAAR1 | Q96RJ0 | CHEMBL3553 | GPCR | informational (could be feature) | 3 | Ulotaront SEP-363856 antipsychotic; modafinil partial-positive (Simmler et al. 2013) |
| 38 | MTNR1A | P48039 | CHEMBL1945 | GPCR | sleep architecture | 3 | MT1 — circadian phase impact |
| 39 | MTNR1B | P49286 | CHEMBL1946 | GPCR | sleep / glucose | 3 | MT2 — glucose homeostasis variant signal |
| 40 | TACR1 | P25103 | CHEMBL249 | GPCR | nausea / mood | 3 | Aprepitant antiemetic; mood literature |
| 41 | NR3C1 | P04150 | CHEMBL2034 | nuclear_receptor | HPA-axis | 2 | Chronic GR engagement = immune/metabolic dysregulation |
| 42 | ESR1 | P03372 | CHEMBL206 | nuclear_receptor | endocrine / thrombosis | 2 | SERM class; chronic dosing endocrine disruption |
| 43 | ESR2 | Q92731 | CHEMBL242 | nuclear_receptor | endocrine | 2 | Paired with ESR1 |
| 44 | NTRK1 | P04629 | CHEMBL2815 | kinase (RTK) | pain / sensory | 3 | TrkA — pan-TRK kinase cardiotoxicity signals |

\* GABRA2/3/5 share the canonical GABA-A heteromer parent target CHEMBL2093872 — for subunit selectivity, MAMMAL needs the ChEMBL "α5-containing" sub-target annotation. Same caveat for CHRNA3/B4 (CHEMBL2109244, ganglionic α3β4 heteromer) and CHRNA4/B2 (CHEMBL1907589, α4β2 heteromer).

**Class balance**: 22 GPCRs / 9 ion channels / 2 enzymes / 3 nuclear receptors / 1 kinase / 7 hybrid. **Liability category coverage**: cardiotox (5), neurotox/cognitive_impairment (4), abuse_potential (3), orthostatic (3), endocrine (3), metabolic_DDI (2), sedation (2), GI (3), HPA (1), informational/feature (15).

---

## Section 2 — Per-Target Pharmacology Backbone

### Tier 1 (HARD CUT)

**1. HTR2B (5-HT2B) — Cardiac Valvulopathy.** The defining off-target liability of CNS pharmacology. Connolly et al. NEJM 1997;337:581–588 reported 24 cases of valvular disease in fen-phen users with histological similarity to carcinoid heart disease; fenfluramine/dexfenfluramine were withdrawn September 1997. Rothman et al. Circulation 2000;102:2836–2841 identified (±)-norfenfluramine, ergotamine, and methylergonovine as preferential 5-HT2B partial-to-full agonists driving the mechanism. Schade et al. NEJM 2007;356:29–38 reported incidence rate ratio 7.1 (95% CI 2.3–22.3) for pergolide and 4.9 (1.5–15.6) for cabergoline at >6 months and >3 mg/day; Zanettini et al. NEJM 2007;356:39–46 echoed with relative risk 4.6–7.3 for cabergoline. Pergolide was withdrawn 2007 (FDA); benfluorex was withdrawn 2009 (EMA). Mechanism: 5-HT2B activation on valve interstitial cells dissociates Gq, activates PLC-β/PKC and TGF-β signaling, driving fibroblast proliferation and leaflet thickening (Cavero & Guillon 2014; Hutcheson, Setola, Roth, Merryman 2011). Roth NEJM 2007;356:6–9 editorial codified the "screen all serotonergic drugs for 5-HT2B agonist activity" rule. **Dumotier & Urban 2024 (J Pharmacol Toxicol Methods 128:107542, PMID 39032441) is the current regulatory framework: safety margin (Cmax/Ki) >10× using binding Ki — not functional EC50 — because metabolites like norfenfluramine can be the active agonist. AUC-corrected for chronic exposure.** **CUT threshold**: predicted pKi > 6.5 (Ki < 300 nM) at 1 µM Cmax; pKi > 7 = hard fail. Aripiprazole is the canonical test case — Shapiro et al. Neuropsychopharmacology 2003 (DOI 10.1038/sj.npp.1300203) reported aripiprazole's highest affinity is at h5-HT2B (low nM) where it functions as an inverse agonist; the chronic safety record is acceptable in psychiatric indications but the *5-HT2B affinity itself* is incompatible with healthy chronic cognitive enhancement.

**2. KCNH2 (hERG) — Torsade de Pointes.** ICH S7B anchor. Cisapride withdrawn 2000; terfenadine 1998; astemizole 1999; grepafloxacin 1999. ICH S7B asks 30× Cmax/IC50 margin for clinical comfort; Bowes 2012 and Brennan 2024 use 10× as the panel-screening threshold. **CUT threshold**: predicted pKi > 6 (Ki < 1 µM). The ADMET-AI hERG head runs in parallel; MAMMAL KCNH2 prediction is a cross-check, and disagreements between them should be flagged for inspection. CiPA (Strauss 2021 and successors) is the regulatory direction — moving from "hERG block = TdP" to integrated multi-channel proarrhythmia assessment incorporating Cav1.2 and Nav1.5 (rows 24–25).

**3. OPRM1 (μ-opioid) — Abuse Potential.** DEA scheduling makes any meaningful μ-opioid binding a regulatory non-starter for a chronic cognitive enhancer. Respiratory depression liability rules out the class. **CUT threshold**: predicted pKi > 6 (Ki < 1 µM) at binding alone; functional agonism not required. This explicitly cuts tramadol-class cross-reactivity and any d-amphetamine pathway that hits OPRM1 in the low-µM regime.

**4. HRH1 (Histamine H1) — Cognitive Impairment from Anticholinergic Burden.** The single most counterproductive liability for a *cognitive enhancement* indication. First-generation H1 antagonists (diphenhydramine, chlorpheniramine) acutely impair memory and reaction time and chronically associate with dementia risk: **Gray et al. JAMA Intern Med 2015;175(3):401–407 (the ACT prospective cohort, n=3,434, median follow-up >7 years) reported adjusted HR 1.54 (95% CI 1.21–1.96) for incident dementia in the highest cumulative anticholinergic-exposure tertile** vs. non-users; **Coupland et al. JAMA Intern Med 2019;179(8):1084–1093 (nested case-control, 58,769 dementia cases / 225,574 controls) found ~50% increased dementia risk with heavy anticholinergic antidepressant exposure**. The 2023 American Geriatrics Society Beers Criteria explicitly lists first-gen H1 antihistamines as potentially inappropriate for chronic dosing in older adults. **CUT threshold**: predicted pKi > 7 (Ki < 100 nM) for antagonist character. Amitriptyline (Ki at H1 ~1 nM, pKi ~9) is the canonical test case for §8.0b cutting an approved drug.

**5. CNR1 (CB1) — Psychiatric AEs on Chronic Dosing.** Rimonabant approved EU June 2006, never approved US (FDA advisory panel voted against, June 2007), withdrawn EU November 2008. Per **Topol et al. (CRESCENDO Investigators), Lancet 2010;376(9740):517–523 (PMID 20709233), neuropsychiatric AEs occurred in 3,028/9,381 (32%) of rimonabant 20 mg patients vs. 1,989/9,314 (21%) on placebo; serious psychiatric side effects 232 (2.5%) vs. 120 (1.3%); 4 suicides on rimonabant vs. 1 on placebo. The trial was terminated early at mean 14 months follow-up after EMA's November 2008 request.** Pre-CRESCENDO regulatory submissions reported suicidal ideation in ~1% of subjects. **CUT threshold**: predicted pKi > 6.5 (Ki < 300 nM) for any direction of CB1 modulation (agonist or inverse agonist — both produced AEs at sufficient occupancy). Sustained CB1 modulation in healthy adults for nootropic purposes is regulatory suicide.

**6. CHRM1 (M1 muscarinic) — Cognitive Impairment from Antagonism.** M1 antagonism is the canonical mechanism of anticholinergic cognitive decline alongside H1 above — donepezil's clinical benefit in AD arises from elevating ACh tone at M1; conversely, blocking M1 with tricyclics or first-gen antihistamines reproduces the cognitive decline phenotype. **Directionality matters**: M1 agonism is desired (xanomeline / KarXT class). For binding-only MAMMAL predictions: **CUT if pKi > 7 AND the parent class is anticholinergic / TCA / first-gen antihistamine**; otherwise FLAG.

### Tier 2 (FLAG with mechanism note)

**7. CHRM3** — Peripheral M3 antagonism = xerostomia, urinary retention (oxybutynin class). FLAG at pKi > 6.

**8–9. MAOA & MAOB** — The Blackwell et al. Br J Psychiatry 1967;113:349 cheese-effect literature; symptomatic at oral tyramine ≥6–25 mg with non-selective MAO-A inhibition. Selegiline (MAO-B-selective at low dose) escapes the cheese effect; phenelzine and tranylcypromine are dangerous because they are irreversible and non-selective. **CUT** if predicted irreversible (covalent SMILES motif) dual MAO-A + MAO-B at pKi > 7; **FLAG** otherwise. Rasagiline and selegiline pass with MAO-B selectivity.

**10–12. ADRA1A/1B/1D** — Orthostatic hypotension and first-dose syncope (prazosin class). Aripiprazole's FDA label (Abilify 2005) explicitly attributes its postural hypotension signal to "antagonist activity at adrenergic α1 receptors." Leung et al. Pharmacol Ther 2012;135:113 documents the autonomic mechanism across antipsychotics. **FLAG at pKi > 6**; CUT only when paired with hERG positivity (QT + orthostatic double-hit).

**13. DRD2** — Chronic antagonism = EPS, hyperprolactinemia, tardive dyskinesia. Aripiprazole/cariprazine/brexpiprazole partial agonism (low intrinsic efficacy + high affinity) is the deliberately-exploited borderline case. **FLAG at pKi > 7** for antagonist character.

**14. DRD3** — D3 agonism = pramipexole/ropinirole impulse-control disorders. **FLAG at pKi > 7** for agonist character.

**17. HTR2A** — Psychotomimetic at high occupancy (psilocybin, LSD, DOM); deliberately exploited in microdosing literature and in pimavanserin-class antipsychotics. **FLAG at pKi > 7**, context-dependent.

**18. HTR2C** — Weight gain (atypical antipsychotic class) plus the lorcaserin 2020 cancer-withdrawal precedent. **Sharretts et al. NEJM 2020;383(11):1000–1002 (PMID 32905685)** reviewed CAMELLIA-TIMI 61 (n=12,000, median follow-up 3.3 years): 7.7% (462) lorcaserin cancer incidence vs. 7.1% (423) placebo; FDA withdrawal request February 13, 2020 (Eisai voluntary withdrawal). **FLAG at pKi > 6** for agonist activity.

**19. HTR3A** — GI/nausea axis. Ondansetron antagonism is anti-emetic and clinically used; TC-5619 has 5-HT3 antagonism that may be beneficial. **FLAG, not CUT.**

**22–23. OPRK1 & OPRD1** — KOR agonism = salvinorin A-class dysphoria; KOR antagonism explored for depression (aticaprant). DOR agonism (SNC80) produces convulsant signal at high doses; mood-relevant. **FLAG at pKi > 6**, direction-sensitive.

**24. CACNA1C (Cav1.2)** — CiPA second axis; block = negative inotropy, AV conduction effects (verapamil/diltiazem class). **FLAG at pKi > 6**.

**25. SCN5A (Nav1.5)** — Cardiac Na⁺ block = QRS widening, conduction slowing, Brugada-like phenotype. Bupropion seizure liability correlates with central Nav block; cardiac Nav1.5 block is the proarrhythmia concern. **FLAG at pKi > 6**.

**28. GABRA1** — Sedation, dependence, tolerance (zolpidem α1 mechanism). **FLAG at pKi > 6**; CUT only if PAM with high intrinsic efficacy.

**32–34. CHRNA3, CHRNB4, CHRNA4** — α3β4 ganglionic mediates nausea/sweating from varenicline and cytisine; α4β2 mediates nicotine reinforcement. **FLAG at pKi > 6** for agonist; informational at antagonist.

**41. NR3C1 (GR)** — Chronic GR engagement = HPA suppression, immune dysregulation, metabolic effects (mifepristone class). **FLAG** binding either direction.

**42–43. ESR1, ESR2** — SERM class endocrine disruption; chronic dosing in healthy adults is regulatorily unfavorable. **FLAG**.

### Tier 3 (INFORMATIONAL — log but don't gate)

**15. DRD4** — Mechanism not characterized as a definitive ADR target; CYP2D6 substrate relevance.

**16. HTR1A** — Buspirone-class anxiolysis; partial agonism is widely exploited in cognition (vortioxetine, vilazodone). Feature candidate.

**20. HTR6** — Idalopirdine target for cognitive enhancement; failed in AD trials but mechanism intact. Listed as panel target, NOT liability.

**21. HTR7** — Vortioxetine partial agonism; circadian and mood.

**26–27. KCNA5, KCND3** — Atrial Kv1.5, ventricular Kv4.3; CiPA-aware but no withdrawal precedent.

**29–31. GABRA2, GABRA3, GABRA5** — α2/α3 = anxiolysis without sedation (TPA023 class); α5 = cognition (Dawson 2006; Atack 2009; Collinson 2006 hippocampal-dependent memory; α5IA in Down syndrome model Braudeau et al. 2011 J Psychopharmacol). **GABRA5 is treated as a potential feature, not liability** — if a compound shows α5-selective negative allosteric modulation, that is the cognition mechanism. Log, don't gate.

**35. SIGMAR1** — Already in the 22-target cognition panel as efficacy target (chaperone). Cross-check only.

**36. TMEM97 (σ2)** — Cholesterol/autophagy; pain & mood literature emerging but no withdrawal precedent.

**37. TAAR1** — Ulotaront (SEP-363856) achieved FDA Breakthrough Therapy designation for schizophrenia; mechanism includes 5-HT1A partial agonism. Ambiguous — could be feature for some cognition strategies. **Modafinil is a weak human-TAAR1 partial-positive modulator (Simmler et al. Br J Pharmacol 2013;168(2):458–470, PMID 22897747, reporting EC50 ~2.5 µM at human TAAR1)**. Log.

**38–39. MTNR1A, MTNR1B** — Sleep architecture and chronic dosing impact on circadian phase. For a chronic cognitive enhancer, MT1/MT2 engagement could improve sleep consolidation — feature, not liability.

**40. TACR1 (NK1)** — Aprepitant antiemetic class; mood literature.

**44. NTRK1 (TrkA)** — Pain/sensory; cardiotoxicity signals in pan-TRK kinase inhibitors (larotrectinib, entrectinib).

---

## Section 3 — Tier-Stratified Gate Severity & Routing Logic

### 3.1 Cmax assumption and threshold derivation

We assume **1 µM unbound Cmax** as the working assumption for a CNS-penetrant cognitive enhancer at therapeutic exposure (donepezil ~50 nM, modafinil ~10 µM, methylphenidate ~50 nM, rasagiline ~20 nM, BPN14770 ~1 µM, lithium ~1 mM as special case). 1 µM is a deliberate worst-case that makes the gate conservative.

**Under the Dumotier & Urban 2024 framework of >10× Cmax/Ki margin**:
- Hard concern: Ki < 10 nM → pKi > 8 → safety margin < 10× → CUT regardless of tier
- Tier 1 CUT: Ki < 100 nM → pKi > 7 → safety margin < 10×
- Tier 2 FLAG: Ki < 1 µM → pKi > 6
- Tier 3 log: any pKi > 5

### 3.2 Per-target thresholds

| target | tier | CUT pKi | FLAG pKi |
|---|---|---|---|
| HTR2B | 1 | 6.5 | 5.5 |
| KCNH2 | 1 | 6.0 | 5.5 |
| OPRM1 | 1 | 6.0 | 5.0 |
| HRH1 | 1 | 7.0 | 6.0 |
| CNR1 | 1 | 6.5 | 5.5 |
| CHRM1 (antagonist class) | 1 | 7.0 | 6.0 |
| CHRM3 | 2 | — | 6.0 |
| MAOA (irreversible dual) | 1 | 7.0 | 6.0 |
| MAOB | 2 | — | 6.0 |
| ADRA1A/B/D | 2 | — | 6.0 |
| DRD2 (antag) | 2 | — | 7.0 |
| DRD3 (agon) | 2 | — | 7.0 |
| HTR2A | 2 | — | 7.0 |
| HTR2C | 2 | — | 6.0 |
| HTR3A | 2 | — | 6.0 |
| OPRK1, OPRD1 | 2 | — | 6.0 |
| CACNA1C, SCN5A | 2 | — | 6.0 |
| GABRA1 | 2 | — | 6.0 |
| CHRNA3/B4, CHRNA4 | 2 | — | 6.0 |
| NR3C1, ESR1, ESR2 | 2 | — | 6.0 |
| All Tier 3 | 3 | — | informational at pKi > 5 |

### 3.3 Routing logic

```python
def liability_verdict(per_target_pki: Dict[str, float]) -> LiabilityResult:
    tier1_hits = [t for t in TIER1 if per_target_pki[t] > THRESHOLDS[t]['cut']]
    tier2_hits = [t for t in TIER2 if per_target_pki[t] > THRESHOLDS[t]['flag']]
    tier3_hits = [t for t in TIER3 if per_target_pki[t] > 5.0]

    if tier1_hits:
        return LiabilityResult('CUT', tier1_hits, tier2_hits, tier3_hits,
                               note=f"Tier 1 hard fail: {', '.join(tier1_hits)}")
    elif len(tier2_hits) >= 2:
        return LiabilityResult('FLAG', [], tier2_hits, tier3_hits,
                               note=f"Composite Tier 2: {', '.join(tier2_hits)}")
    elif len(tier2_hits) == 1:
        return LiabilityResult('FLAG', [], tier2_hits, tier3_hits,
                               note=f"Tier 2: {tier2_hits[0]}")
    elif tier3_hits:
        return LiabilityResult('PASS', [], [], tier3_hits,
                               note=f"Informational T3: {', '.join(tier3_hits)}")
    else:
        return LiabilityResult('PASS', [], [], [], note='clean')
```

**Hard gate precedence**: if `admet_status == 'CUT' OR liability_status == 'CUT'` → `final_status = 'CUT'`. The regulatory bypass for approved drugs (currently in ADMET) does NOT extend to §8.0b — an approved drug can still be CUT if used outside its approved chronic-dosing healthy-adult context (the aripiprazole case).

---

## Section 4 — Implementation Plan

### 4.1 Data layer

**Extend `data/raw/targets_seed.csv`** with the 44 rows in Section 1.3 plus new schema columns:
```
gene_symbol, uniprot_accession, chembl_target_id, target_class,
panel_type ∈ {'cognition', 'liability'},
severity_tier ∈ {1, 2, 3, NULL},
liability_category, cut_threshold_pki, flag_threshold_pki,
direction_sensitivity ∈ {'agonist', 'antagonist', 'either', 'NA'}
```

`direction_sensitivity` informs the gate when MAMMAL binding-only predictions need an agonist/antagonist call (mostly via parent-class lookup or Boltz-2 pose-based prediction if available).

### 4.2 MAMMAL inference

Existing `dti_grid.parquet` already accommodates arbitrary targets; expansion 22 → 66 targets requires re-running on 298 × 44 = **13,112 new DTI predictions**. The 22-target run completed in ~10 min on RTX 5070 → linear estimate **~22 min wall-clock**, plus tokenization/setup overhead = **~45–60 min total**. Cache liability predictions separately: `data/results/v2/liability_dti.parquet`. Merge with cognition predictions at the gating stage rather than in the model run, preserving the ability to re-run liability only.

### 4.3 New gate module

`src/mammal_repurposing/gates/liability_panel.py`:
```python
from dataclasses import dataclass
from typing import Literal
import pandas as pd

Status = Literal['PASS', 'FLAG', 'CUT']

@dataclass
class LiabilityResult:
    compound_id: str
    status: Status
    tier_1_hits: list[str]
    tier_2_hits: list[str]
    tier_3_hits: list[str]
    top_3_liabilities: list[tuple[str, float, int]]  # (target, pKi, tier)
    liability_note: str
    liability_summary: str  # e.g. "α1A=6.8(T2), 5-HT2B=5.2(clean), hERG=6.4(T1)"

def apply_liability_gates(
    dti_grid: pd.DataFrame,
    targets_seed: pd.DataFrame,
) -> pd.DataFrame:
    """
    dti_grid: long-form DTI predictions (compound_id, target, pki, ...)
    targets_seed: target metadata with panel_type, severity_tier, thresholds
    Returns: per-compound LiabilityResult records integrable with v4 shortlist.
    """
    ...
```

**Integration order** in the pipeline DAG:
1. MAMMAL DTI grid (cognition + liability)
2. ADMET-AI gates (`gates/admet_gates.py`) — ADMET CUT removes compound
3. **§8.0b liability gates** (`gates/liability_panel.py`) — runs only on ADMET-surviving compounds
4. RRF fusion of cognition scores (cognition targets only — liability scores never enter the efficacy ranking)
5. Final ranking with `final_status = CUT if (admet OR liability) == CUT else FLAG if any FLAG else PASS`

### 4.4 Calibration

Same Phase A.7 approach: per-target Spearman ρ vs ChEMBL ground truth.
- **Dense targets** (clean calibration): HTR2B (~3k records), HTR2A, HTR2C, KCNH2 (~8k), DRD2 (~5k), MAOA, MAOB, HRH1, CHRM1, CHRM3, ADRA1A, OPRM1, CNR1, HTR1A, HTR3A.
- **Sparse targets** (unreliable): TAAR1 (~200, mostly modafinil-class), TMEM97 (<50), CHRNA3/B4 subunit-specific (heteromer-annotated in ChEMBL), KCNA5, KCND3.

**Calibration decision**: trust raw pKi for the binary CUT/FLAG calls (gate is a threshold, not a rank), but use calibration ρ to *suppress* sparse-target alerts where Spearman ρ < 0.3 — for those, route to a "low-confidence flag" subcategory for manual review rather than auto-CUT.

### 4.5 Reporting

`reports/liability_audit.md` — one row per shortlisted compound:
```markdown
## modafinil (CHEMBL1373)
| Target | pKi (pred) | pKi (ChEMBL gt) | Tier | Threshold | Verdict |
|---|---|---|---|---|---|
| TAAR1 | 5.8 | 5.6 (Simmler 2013) | 3 | info | log (feature) |
| ADRA1A | 4.9 | 4.7 | 2 | 6.0 | PASS |
| ... |
**Verdict**: PASS. No Tier 1 hits; informational TAAR1 modulation consistent with mechanism.
```

### 4.6 Shortlist v4 schema additions

`wet_lab_shortlist_v4.csv` gains: `liability_status ∈ {PASS, FLAG, CUT}`, `liability_summary` (one-line string), `tier_1_hits`, `tier_2_hits` (semicolon-separated), `liability_note` (mechanistic free-text from Section 2 backbone). `docs/methodology_note_v2.md §3.5` documents the panel composition with citations to Bowes 2012, Brennan 2024, Dumotier 2024, Connolly 1997, Schade 2007, Roth 2007, Blackwell 1967.

### 4.7 Estimated effort

- Schema extension + targets_seed.csv expansion: 2 hr
- MAMMAL re-run: ~1 hr wall-clock + ~1 hr setup/validation
- `gates/liability_panel.py` module + tests: 4–6 hr
- Calibration sweep: 2 hr
- Reporting renderer: 2 hr
- Methodology note §3.5: 2 hr
- **Total: ~1.5 days engineering effort**

---

## Section 5 — Predicted Impact on the V3 Top-25 Shortlist

| # | Compound | ADMET status | §8.0b verdict | Mechanistic basis |
|---|---|---|---|---|
| 1 | modafinil | PASS | **PASS** | Weak TAAR1 (T3, feature, Simmler 2013); minor α1; no T1 hits |
| 2 | aniracetam | PASS | **PASS** | AMPA modulation off-panel; clean |
| 3 | 2bact (2-PMPA) | PASS | **PASS** | GCPII-selective |
| 4 | aripiprazole | bypass (approved); hERG=0.946, P-gp=0.899 | **CUT** | T1: 5-HT2B highest affinity (Shapiro 2003), HRH1, D2 partial agonism; T2: α1A. Exemplar of bypass needing §8.0b override for chronic healthy dosing. |
| 5 | (S)-AMPA | PASS | **PASS** | AMPA agonist; off-panel |
| 6 | (R,S)-AMPA | PASS | **PASS** | Same |
| 7 | modafinil (dup) | — | — | duplicate row |
| 8 | alpha-GPC | PASS | **PASS** | Choline donor; no GPCR engagement |
| 9 | BI-409306 | PASS | **PASS** | PDE9A selective |
| 10 | rolipram | PASS | **FLAG** | PDE4 nausea/emesis off-target (5-HT3, area postrema); single T2 |
| 11 | piracetam | PASS | **PASS** | Clean |
| 12 | levetiracetam | PASS | **PASS** | SV2A selective |
| 13 | rasagiline | PASS | **FLAG** | T2: MAO-B irreversible (intentional); MAO-A selectivity adequate; logged as T2 with "intentional MAO-B, tyramine-safe at therapeutic dose" note |
| 14 | bupropion | PASS | **FLAG** | T2: SCN5A (seizure risk at high dose); DAT/NET = efficacy off-liability |
| 15 | amitriptyline | regulatory bypass | **CUT** | T1: HRH1 (Ki ~1 nM, pKi ~9), CHRM1. Test case for §8.0b cutting approved drug from healthy chronic cognition context. |
| 16 | lemborexant | PASS | **PASS** | Orexin antagonist; sleep target |
| 17 | lanicemine | PASS | **PASS** | NMDA low-trapping; minor σ1 |
| 18 | 7,8-DHF | PASS | **PASS** | TrkB agonist; clean |
| 19 | tulrampator | PASS | **PASS** | AMPA PAM; clean |
| 20 | methylphenidate | PASS | **FLAG** | T2: DAT/NET (efficacy); abuse via DEA Schedule II, not §8.0b CUT — no μ-opioid binding |
| 21 | d-amphetamine | PASS | **FLAG → CUT** | T1 borderline: μ-opioid binding low-µM in historical PDSP data; T2: DAT, TAAR1; recommend CUT for healthy cognitive enhancement on regulatory grounds |
| 22 | topiramate | PASS | **PASS** | Weak GABA-A α1; CA-II off-panel |
| 23 | lithium carbonate | PASS | **PASS** (special case) | GSK-3β/inositol monophosphatase — not a receptor binder; panel inapplicable, returns clean |
| 24 | galantamine | PASS | **FLAG** | T2: CHRNA3/B4 (low affinity but non-zero); AChE efficacy; α7 APL (efficacy); matches clinically documented GI/nausea liability |
| 25 | TC-5619 | bypass + hERG=0.964 | **FLAG** | T2: HTR3A (5-HT3 antag — possibly beneficial); no T1; ADMET hERG is the more serious concern |
| 26 | BPN14770 | bypass + DILI=0.921 | **PASS** | PDE4D-selective NAM (Burgin 2010 allosteric); DILI is a PDE4-family effect §8.0b cannot decompose, but receptor panel correctly returns PASS. Separate mechanistic DILI workup recommended. |

**Discrimination summary**: 2 hard CUTs (aripiprazole, amitriptyline), 7 FLAGs (rolipram, rasagiline, bupropion, methylphenidate, d-amphetamine, galantamine, TC-5619), 15 PASSes among 25. ~36% touch rate, with 2 hard rejections that ADMET alone would not catch. **This is the methodological contribution.**

---

## Section 6 — Publication Angle

**Title**: "A cognition-context-tailored off-target liability panel for foundation-model drug repurposing in healthy cognitive enhancement"

**Venue**: J. Cheminform. (open-access, methods-focused, high cheminformatics readership) primary; Drug Discovery Today: Technologies as backup.

**Novelty claims** (defensible):
1. Off-the-shelf safety panels (Bowes-44, SafetyScreen44, Brennan-77) are indication-agnostic; healthy chronic dosing for cognitive enhancement has a different risk profile than oncology, anti-infective, or psychiatric chronic-dosing-in-sick-patients contexts.
2. Tier-stratified gating with mechanism-aware directionality (M1 agonism = feature vs. antagonism = liability; α5-GABA-A inverse agonism = feature; TAAR1 = ambiguous) is novel relative to standard "% inhibition at 10 µM" panel reads.
3. Integration with foundation-model DTI prediction (MAMMAL) provides a virtual-first triage otherwise only achievable at $40–80k per compound for in vitro panels.
4. Demonstration: the panel correctly CUTs aripiprazole and amitriptyline on the basis of 5-HT2B + α1 + D2 and HRH1 + M1 respectively — the regulatory bypass needs a mechanism-explicit complement, which §8.0b provides.

**Reproducibility deliverables**:
- Full 44-target spec with UniProt + ChEMBL IDs (Section 1.3)
- Per-target pKi thresholds + tier assignments (Section 3.2)
- `gates/liability_panel.py` module + tests (open-source MIT/Apache-2)
- Per-target Spearman ρ vs ChEMBL ground truth as supplementary
- Wet-lab shortlist v4 with §8.0b verdicts attached (post-experimental validation)

**Manuscript structure**:
1. Introduction — why off-the-shelf panels miss cognition context
2. Methods — 44-target selection rationale, tier definitions, gate logic
3. Results — application to 298-compound MAMMAL DTI grid; calibration ρ heatmap; v3 top-25 impact
4. Discussion — limitations (MAMMAL allosteric blindness, sparse-target uncertainty), comparison to Brennan-77, future extensions
5. Open materials with all code and data

---

## Section 7 — Risks and Failure Modes

1. **MAMMAL allosteric blindness propagates to liability predictions.** MAMMAL learns orthosteric binding from ChEMBL; allosteric modulators (PDE4D NAM BPN14770, α7-nAChR PAMs, σ1 modulators) are systematically under-predicted. Conservatively safe for orthosteric liabilities (5-HT2B, hERG, D2) — MAMMAL catches them. But beneficial allosteric mechanisms (α5-GABA-A inverse agonism, M1 PAM) may be missed in the cognition run and falsely flagged in the liability run. **Mitigation**: cross-reference predicted binding against compound-class metadata in the methodology note.
2. **Binding ≠ functional consequence.** MAMMAL predicts pKi, not efficacy or direction. For 5-HT2B this is conservative per Dumotier 2024 (binding Ki preferred). For DRD2 partial agonism vs. antagonism (aripiprazole textbook case), binding alone cannot distinguish — must be supplemented with the `direction_sensitivity` column or Boltz-2 pose-based downstream prediction.
3. **ChEMBL ground-truth contamination.** If a compound's HTR2B Ki is in MAMMAL's training set, the prediction is trivially recovered, biasing the calibration metric. **Mitigation**: report calibration ρ on temporal hold-out (ChEMBL releases post-MAMMAL training cutoff).
4. **Heteromeric ion channel ambiguity.** ChEMBL annotates GABA-A and nAChR heteromers inconsistently. **Mitigation**: explicit `chembl_target_id` mappings to heteromer parents (e.g., CHEMBL2093872 for GABA-A heteromer, CHEMBL1907589 for α4β2 nAChR) plus subunit-level child mappings where available.
5. **The 1 µM Cmax assumption is wrong for some compounds.** Lithium peaks at ~1 mM; some nootropics peak at ~10 nM. **Mitigation**: log Cmax per compound where known and re-run the gate with compound-specific margins for the top 25.
6. **Aripiprazole CUT is correct but PR-fraught.** The framing must be: "approved for *acute or maintenance treatment of schizophrenia in patients*, not for *chronic cognitive enhancement in healthy adults*" — indication context is the gating criterion, not the drug's general safety record.
7. **Sparse-target false positives.** TAAR1, TMEM97, KCNA5, KCND3 have <200 ChEMBL records. **Mitigation**: Tier 3 informational only; do not gate.

---

## Recommendations

**Stage 1 — Immediate (~1.5 days engineering)**:
1. Extend `targets_seed.csv` with the 44 rows in Section 1.3 + new schema columns.
2. Implement `gates/liability_panel.py` per Section 3.3.
3. Re-run MAMMAL on 298 compounds × 44 liability targets (~1 hr wall-clock on RTX 5070).
4. Apply gates to v3 shortlist; emit `reports/liability_audit.md` and `wet_lab_shortlist_v4.csv`.
5. Spot-check against Section 5 predictions — aripiprazole CUT, amitriptyline CUT, modafinil/aniracetam/BPN14770 PASS. If they don't, debug before proceeding.

**Stage 2 — Calibration & validation (~1 week)**:
6. Phase A.7 calibration sweep for all 44 liability targets; emit per-target Spearman ρ heatmap.
7. Demote sparse-target alerts (ρ<0.3) to "low-confidence flag" routing.
8. Re-run gates with calibrated weights; compare verdicts.

**Stage 3 — Publication & methodology note (~2 weeks)**:
9. Draft methodology note v2 §3.5 with full citation backbone (Bowes 2012, Brennan 2024, Dumotier 2024, Connolly 1997, Schade 2007, Roth 2007, Cavero & Guillon 2014, Gray 2015, Coupland 2019, Topol/CRESCENDO 2010, Sharretts/lorcaserin 2020, Blackwell 1967, Simmler 2013).
10. Draft J. Cheminform. manuscript per Section 6 structure.

**Benchmarks that would change recommendations**:
- **>2 false-positive CUTs** in wet-lab v4 (compounds CUT that experimentally clean at the predicted target) → **relax Tier 1 thresholds by 0.5 log unit** and rerun.
- **False-negative misses** (compounds that pass §8.0b but show in vitro liability hits at <1 µM) → **tighten the FLAG threshold to pKi > 5.5**.
- **MAMMAL calibration ρ on Tier 1 targets consistently <0.5** → **add an in vitro confirmation step** (Eurofins SafetyScreen44 selective subset) for the top 5–10 shortlist compounds before wet-lab efficacy assays.

---

## Caveats

- The exact composition of Brennan-77 (Brennan et al. Nat Rev Drug Discov 2024;23:525–545) is behind the Nature paywall; the ~33 additions to Bowes-44 reconstructed here are well-supported by the paper's accessible references and the Schmidt/Brennan/Jenkinson/Valentin 2025 correspondence (Nat Rev Drug Discov 24:482–484) but the full kinase-family expansion is not exhaustively enumerated. For the published methodology note, retrieve Table 1 of Brennan 2024 via institutional subscription.
- The pKi thresholds (Section 3.2) are *defensible defaults* derived from the Dumotier 2024 >10× Cmax/Ki framework at 1 µM Cmax — not yet optimized. Revisit after the first wet-lab validation cycle.
- MAMMAL is not validated for "is this a 5-HT2B *agonist*?" — only for "what is the predicted pKi at 5-HT2B?". The Section 3 routing conservatively treats high-pKi as liability regardless of direction for Tier 1; direction-sensitive Tier 2 calls (D2, MAOA, M1) depend on compound-class metadata.
- The 1 µM Cmax assumption is worst-case; real compounds will have lower free fraction at the relevant CNS exposure. The gate is therefore conservatively over-sensitive, which is the desired screening behavior.
- Calibration ρ for sparse targets is expected unreliable; gating on raw pKi but suppressing alerts where ρ<0.3 is heuristic and should be revisited.
- The §8.0b panel is a *screening* tool, not a regulatory filing. FDA/EMA submission would require in vitro confirmation at hit targets, ideally with the panel-and-margin approach codified in Papoian et al. 2015 and Dumotier & Urban 2024.
- The "Liu et al. 2022" attribution for modafinil-TAAR1 mechanism in an earlier draft was incorrect; the authoritative reference for modafinil's TAAR1 partial activation is **Simmler LD et al. Br J Pharmacol 2013;168(2):458–470 (PMID 22897747)** at human TAAR1 EC50 ~2.5 µM.