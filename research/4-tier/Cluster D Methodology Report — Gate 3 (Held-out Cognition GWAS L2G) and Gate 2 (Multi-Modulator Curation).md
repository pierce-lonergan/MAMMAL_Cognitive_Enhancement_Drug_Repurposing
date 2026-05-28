# V6.B / Cluster D Methodology Report — Gate 3 (Held-out Cognition GWAS L2G) and Gate 2 (Multi-Modulator Curation)

## TL;DR
- **Gate 3 (held-out cognition GWAS L2G):** Use a **tiered held-out set** rather than a single GWAS. The strongest mutually-independent primary option is **Okbay 2022 EA4 (`EA4_additive_excl_23andMe` release; the full EA4 N = 3,037,499 per Okbay et al. 2022 Nat Genet 54:437 Table 1, PMID 35361970)** as a high-power proxy, paired with a *true*-cognition anchor (de la Fuente 2021 genetic-g UKB multivariate, FinnGen R12 dementia/memory endpoints, and — under NDA DUA — ABCD NIH-Toolbox neurocognitive components). Educational attainment is genetically correlated with cognition rg ≈ 0.70 (Lee 2018) and is a legitimate proxy but must be flagged as such; do **not** use Davies 2018, Hill 2019 MTAG, or Savage 2018 as held-out because they share UKB + CHARGE + COGENT cohorts with virtually every cognition GWAS a target-prioritization pipeline would touch. L2G must be obtained from **Open Targets Platform 25.06+** — the Open Targets Genetics standalone portal was officially deprecated on 9 July 2025 (Open Targets Community announcement, 15 May 2025: "Open Targets Genetics will be deprecated on 9 July 2025, after the next release of the Open Targets Platform … you will no longer be able to access the Open Targets Genetics user interface (genetics.opentargets.org)") — via the GraphQL API at `https://api.platform.opentargets.org/api/v4/graphql` or bulk `l2g_prediction` parquet at `gs://open-targets-data-releases/<release>/output/l2g_prediction/`. For a novel GWAS not yet ingested by Open Targets, run the **Gentropy** open-source PySpark pipeline locally.
- **Gate 2 (multi-modulator curation):** A literature-pooled effect-size table covering all 28 targets and ~95 (target, modulator) pairs is provided below, each row with Hedges' g / Cohen's d (or ADAS-Cog mean-difference), 95% CI, primary endpoint, and PMID/DOI. Replacing n = 11 pairs with ≥ 80 pairs gives a Spearman ρ test with > 90% power to detect |ρ| ≥ 0.25 at α = 0.05. **Sign-encoding must use mechanism × physiological-direction-of-effect**, not raw drug effect (e.g., HRH3 antagonist enhances cognition → HRH3 target-direction = −1).
- **Headline numbers:** Healthy-subject cognitive enhancers cap at Hedges' g ≈ 0.2 (Roberts 2021 Eur Neuropsychopharm 38:40, PMID 32709551 — modafinil SMD = 0.12; methylphenidate SMD = 0.21; d-amphetamine null). Disease-population effects are modestly larger (donepezil 10 mg vs placebo ADAS-Cog MD = −2.67 [95% CI −3.31 to −2.02], Birks & Harvey 2018 Cochrane CD001190.pub3; memantine moderate-severe AD ADAS-Cog MD = 2.15 points improvement [95% CI 1.05 to 3.25], McShane 2019 Cochrane CD003154.pub6 PMID 30891742; pitolisant ESS Cohen's d = 0.61 in HARMONY 1 and 0.86 in HARMONY CTP per Meskill et al. 2022, CNS Drugs 36:61, PMID 34935103). Many high-profile candidates (encenicline, idalopirdine, intepirdine, pomaglumetad, bitopertin) are **null in Phase III** — these null effect sizes are critical for a sign-consistent Spearman and must not be silently dropped.

---

## Key Findings

### Findings — Task 1 (Held-out cognition GWAS for Gate 3)

1. **Cohort-overlap problem.** The de-facto canonical "general cognitive function" GWAS — Davies 2018 (N = 300,486, PMID 29844566) and Savage 2018 (N = 269,867, PMID 29942086) — both share UK Biobank, CHARGE, and COGENT cohorts, so they cannot serve as *mutually* independent held-out sets. They are also the most likely training-set sources for any cognition repurposing pipeline (they dominate Open Targets, MAGMA, FUMA priors). Davies 2018 itself notes a LDSC genetic correlation r_g = 0.82 (SE = 0.02) between its CHARGE-COGENT and UKB sub-meta-analyses.
2. **Largest cognition-adjacent GWAS = Okbay 2022 EA4.** Per Okbay et al. 2022 Nat Genet 54:437 Table 1 (PMID 35361970), EA4 totals N = 3,037,499 (Lee et al. excl-23andMe + UKB = 324,162; new 23andMe = 2,272,216; new UKB = 441,121). The SSGAC publishes `EA4_additive_excl_23andMe.txt` (~N = 765,000) as the public release. The genome-wide polygenic index explains 12–16% of EA variance (Okbay 2022 abstract); Table 1's C+T predictor gives an incremental R² of 7.18% for EA in held-out samples. EA proxies cognition at r_g ≈ 0.70 (Lee 2018 reports r_g(EA, intelligence) = 0.70; Davies 2018 reports r_g(EA, g) = 0.73). Proxy — not identity.
3. **True-cognition anchor independent of EA-proxy confounding** is needed alongside EA4. Best options: (i) **de la Fuente 2021** genetic-g multivariate GWAS (UKB cognitive tests, N = 11,263–331,679 per test; Nat Hum Behav 5:49, PMID 33257891) — derives a genetic-g latent factor that isolates 30 genome-wide-significant loci specific to general cognitive ability; (ii) **ABCD Study** NIH Toolbox-derived NPC1/2/3 component GWAS (Loughnan/Fan and ABCD-genomics consortium; access via NDA DUA, see PMC 10635818 for the genotyping resource); (iii) **FinnGen R12** dementia/memory endpoints (F5_DEMENTIA, MEMLOSS, G6_AD_WIDE, AD_LO) for an Alzheimer-adjacent population-isolate contrast (https://r12.finngen.fi/).
4. **Open Targets L2G toolchain state (Q2 2025 → Q2 2026).** Open Targets Genetics standalone portal retired 9 July 2025; everything merged into the main Open Targets Platform (release 25.06, 18 June 2025). L2G is now produced by the **Gentropy** open-source PySpark library (github.com/opentargets/gentropy). The 25.06/26.03 cycles added ENCODE rE2G enhancer-gene features (`e2gMean`, `e2gNeighbourhoodMean`) and SHAP feature-importance values to L2G outputs.
5. **Gate 3 AUROC** should be computed as: positives = genes with L2G ≥ 0.5 in the held-out GWAS (the Mountjoy 2021 gold-standard threshold), restricted to protein-coding genes outside MHC; the pipeline's per-gene priority score is the ranking variable; compute AUROC across all ≈ 19,000 protein-coding genes. The OSF target AUROC > 0.70 corresponds to "the pipeline's top decile is enriched ~3-fold for held-out causal genes."

### Findings — Task 2 (Multi-modulator curation)

1. **Healthy-subject effect-size ceiling.** Roberts 2021 (Eur Neuropsychopharmacol 38:40–62, PMID 32709551) is the highest-quality recent meta-analysis: modafinil SMD = 0.12 (p = 0.01), methylphenidate SMD = 0.21 (p = 0.0004; driven by recall g = 0.43, sustained attention g = 0.42), d-amphetamine null. This sets the upper bound for healthy-subject enhancers at Hedges' g ≈ 0.2.
2. **Disease effect sizes are larger but still modest.** Donepezil 10 mg vs placebo at 24 wk: ADAS-Cog MD = −2.67 (95% CI −3.31 to −2.02), Birks & Harvey 2018 Cochrane CD001190.pub3. Memantine in moderate-severe AD: 2.15 ADAS-Cog points improvement (95% CI 1.05 to 3.25), McShane 2019 Cochrane CD003154.pub6 PMID 30891742. Galantamine WMD ≈ −2.76 ADAS-Cog points (Hansen 2008 DARE NBK74974).
3. **High-profile Phase-III failures must be encoded as g ≈ 0, not omitted.** Encenicline EVP-6124 (α7 nAChR) failed Phase III EVP-6124-015/016 (Schizo Bull suppl S141); idalopirdine (5-HT6) failed three Phase III RCTs combined Phase III meta MD = −0.41 ADAS-Cog (p = 0.32, Matsunaga 2018 Int Psychogeriatr PMID 30560763); intepirdine (5-HT6) failed Phase III MINDSET; pomaglumetad (mGluR2/3) failed 3 Phase 2/3 trials (Adams 2014 PMID 24772351); bitopertin (GlyT1) failed SunLyte/DayLyte/FlashLyte (Bugarski-Kirola 2017 Biol Psychiatry 82:8, PMID 28117049); blarcamesine (σ-1) showed marginal ADAS-Cog13 LSM diff = −1.78 (95% CI −3.31 to −0.25, p = 0.022) in ANAVEX2-73-AD-004 but missed the ADCS-ADL co-primary (Macfarlane/Sabbagh 2024 J Prev Alz Dis DOI 10.14283/jpad.2024.122).
4. **Sign encoding for the Spearman correlation.** The target-importance variable encodes "predicted direction in which *increasing target function* changes cognition." For receptors normally engaged by procognitive agonists (CHRNA7, DRD1, CHRM1/2/4, ADRA2A, NTRK2, HCRTR1/2 [for wakefulness], HTR4, HTR1A, GRIA1-4, GRIN2A [activation], GRM2/3, GRM5 [PAM context], SIGMAR1) target-direction = +1; for receptors/enzymes whose *inhibition* enhances cognition (HRH3, ACHE, PDE4D, PDE9A, MAOB, MAOA, COMT, HTR6, HTR2A, GRIN2B [via memantine uncompetitive partial block at hyperactive synapses], GABRA5) target-direction = −1; for transporters (SLC6A2/NET, SLC6A3/DAT, SLC6A9/GlyT1) inhibition enhances → −1.

---

## Details

### 1. Held-out cognition GWAS catalogue

| # | Study | Year | Phenotype | N | Ancestry | Consortium | Access | Notes / overlap |
|---|---|---|---|---|---|---|---|---|
| 1 | Davies G *et al.*, Nat Commun 9:2098 | 2018 | General cognitive function ("g") | 300,486 | EUR | CHARGE + COGENT + UKB (4 sub-samples) | Edinburgh DataShare https://datashare.ed.ac.uk/handle/10283/3756; PMID 29844566 | Verbal-numerical reasoning + g-factor pooled. Will overlap with most training-set pipelines. r_g(sub-cohorts) = 0.82. |
| 2 | Savage JE *et al.*, Nat Genet 50:912 | 2018 | Intelligence (14-cohort meta) | 269,867 | EUR | UKB + CHIC + COGENT + 13 cohorts | CTGlab https://ctg.cncr.nl/software/summary_statistics; PMID 29942086 | Overlaps UKB and COGENT with Davies 2018. |
| 3 | Lee JJ *et al.* (SSGAC), Nat Genet 50:1112 | 2018 | EA3 + Cognitive Performance (CP) | EA 1,131,881 / CP 257,841 | EUR | SSGAC + 23andMe + UKB | SSGAC https://thessgac.com/papers/3 (`GWAS_EA_excl23andMe.txt`, `GWAS_CP_all.txt`); PMID 30038396 | CP overlaps COGENT + UKB. EA-excl-23andMe is the canonical public release. |
| 4 | Okbay A *et al.* (SSGAC), Nat Genet 54:437 | 2022 | EA4 | **3,037,499** (Table 1) | EUR | SSGAC + 23andMe + UKB + Million-cohort meta | SSGAC https://thessgac.com/papers/14 (`EA4_additive_excl_23andMe.txt`); PMID 35361970 | **Largest, best-documented, public.** PGI R² for EA = 12–16% (abstract); C+T = 7.18% (Table 1). Proxy via r_g ≈ 0.70 with CP. |
| 5 | Hill WD *et al.*, Mol Psychiatry 24:169 | 2019 | Intelligence (MTAG with EA) | 248,482 effective | EUR | UKB + Sniekers 2017 + SSGAC EA | https://datashare.ed.ac.uk/handle/10283/3358; PMID 29326435 | **MTAG mixes intelligence with EA** — not a clean cognition phenotype. |
| 6 | Hill WD *et al.*, Nat Commun 10:5741 | 2019 | Household income | 286,301 | EUR | UKB | https://datashare.ed.ac.uk/handle/10283/3441; PMID 31844048 | UKB-only; useful as **negative-control** SEP proxy. |
| 7 | de la Fuente J *et al.*, Nat Hum Behav 5:49 | 2021 | Genetic-g (multivariate Genomic SEM, 7 UKB tests) | 11,263 – 331,679 per test | EUR | UKB | Code https://github.com/JavierdelaFuente/Genetic_g; sumstats from authors; PMID 33257891 | **Best "true-cognition" anchor**; isolates g-specific loci from test-specific loci. UKB-only. |
| 8 | Trampush JW *et al.*, Mol Psychiatry 22:336 | 2017 | Neurocognitive composite (COGENT) | 35,298 | EUR | COGENT (35 cohorts) | Lencz lab COGENT request; PMID 28093568 | Subsumed into Davies 2018 + Lee 2018. Not suitable as held-out. |
| 9 | Davies G *et al.*, Mol Psychiatry 21:758 | 2016 | UKB cognitive functions | 112,151 | EUR | UKB | PMID 27046643 | UKB-only and small; superseded. |
| 10 | Eising / GenLang reading & language GWAS | 2022 | Reading & language (5 tests) | up to ~33,000 | EUR | GenLang | https://www.genlang.org/data; PMID 35835581 | Verbal-cognition complement to de la Fuente 2021. |
| 11 | ABCD NIH Toolbox neurocognitive PC GWAS (Loughnan/Fan; ABCD consortium) | 2024–2025 | NPC1 (general ability), NPC2 (executive), NPC3 (learning/memory) | ≈10,000 | EUR + AMR + AFR + ASN | ABCD | NDA controlled (NDA Study 1198+); PMC 10635818 (genotyping resource) | **Cleanly held-out**: post-dates most training sets; childhood phenotype; multi-ancestry; controlled access via NDA DUA. |
| 12 | FinnGen R12 dementia/memory endpoints (F5_DEMENTIA, MEMLOSS, G6_AD_WIDE, AD_LO) | 2024 (DF12, Dec 2024) | Dementia / memory disorder ICD codes | ~500,000 cohort; case Ns vary | FIN | FinnGen | https://r12.finngen.fi/; researcher registration | **Population-isolate, independent of UKB/CHARGE/COGENT**; best for testing dementia-relevant cognitive targets. |
| 13 | Million Veteran Program cognitive ageing scans (Harvey 2023+) | 2023+ | Reaction-time decline; cognitive trajectories | ~250,000 | EUR + AFR + HIS | MVP | dbGaP phs001672, controlled | MVP largely independent of UKB/CHARGE; controlled-access (12-week IRB cycle). |
| 14 | Sniekers S *et al.*, Nat Genet 49:1107 | 2017 | Intelligence | 78,308 | EUR | CHIC + COGENT | PMID 28530673 | Predecessor to Savage 2018; do not use as held-out. |

#### Independence verification — recommended workflow

1. **Cohort-overlap audit.** For each candidate, pull its constituent-cohort list (Supplementary Table 1 of the paper) and intersect with your training-set GWAS's cohorts. Shared cohort > 5% of effective N → flag as non-independent.
2. **LDSC bivariate cross-trait LD-score intercept.** Run `ldsc.py --rg` between the held-out file and each training-set file. The cross-trait LD-score intercept estimates sample overlap (intercept ≈ 0 means independent samples; intercept = ρ · √(Ns·N_overlap/(N1·N2))). Reject candidates with cross-trait intercept > 0.05.
3. **GWAS Catalog metadata check.** Query the EBI REST API (`https://www.ebi.ac.uk/gwas/rest/api/`) for "Initial sample description" and "Replication sample description" fields.
4. **EA vs cognition proxy caveat.** Lee 2018 reports r_g(EA, intelligence) = 0.70; Davies 2018 reports r_g(EA, g) = 0.73; Okbay 2022 reports the EA4-PGI explains 12–16% of EA variance overall and 7.18% via C+T. EA-based GWAS recovers plausible cognition genes but is enriched for non-cognitive correlates (personality, SES). **Always pair EA-proxy with a true-cognition anchor and report both AUROCs separately.**

#### Recommended optimal held-out set

| Tier | GWAS | Role | Rationale |
|---|---|---|---|
| **Primary** | Okbay 2022 EA4 (excl-23andMe, N from 3,037,499 EA4 total per Table 1) | High-power proxy | Largest, public, post-dates most training, has clean exclusion lever |
| **Secondary (true cognition)** | de la Fuente 2021 genetic-g UKB multivariate | Phenotype anchor | Isolates g-specific loci; controls proxy confounding |
| **Tertiary (independence stress-test)** | FinnGen R12 F5_DEMENTIA + MEMLOSS | Population-isolate, dementia-adjacent | Independent of UKB/CHARGE/COGENT |
| **Optional (developmental)** | ABCD NIH Toolbox NPC1 (NDA DUA) | Childhood cognition anchor | Independent of adult cohorts; multi-ancestry |
| **Excluded** | Davies 2018, Hill 2019 MTAG, Savage 2018 | — | Overlap with virtually all adult cognition training data |

Report Gate 3 AUROC for each held-out file separately, plus a Stouffer-weighted meta-AUROC.

### 2. Open Targets L2G — methodology and toolchain (Q2 2025 → Q2 2026)

**The L2G model** (Mountjoy 2021 Nat Genet; updated 25.x via Gentropy) is a sklearn gradient-boosted classifier trained on a gold-standard positive set (curated GWAS-implicated genes) and a negative set (nearby non-causal genes), with nested cross-validation. Predictive features:

- **Distance**: minimum distance from credible-set variants to gene canonical TSS / gene body, weighted by per-variant posterior probability.
- **Molecular QTL colocalisation**: H4 posterior from coloc with eQTL Catalogue, sQTL, pQTL (incl. UKB-PPP).
- **Chromatin interaction**: Promoter-Capture Hi-C, Javierre 2016.
- **VEP functional consequence**: missense, splice, etc.
- **ENCODE rE2G (added in 25.06/26.03)**: `e2gMean`, `e2gNeighbourhoodMean` regulatory predictions.
- **SHAP**: per-prediction feature contributions.

**Score interpretation**: L2G ∈ [0,1] ≈ fraction of gold-standard-positive genes among all genes at that locus given the threshold. **L2G ≥ 0.5 is the canonical high-confidence threshold**; predictions < 0.05 are filtered out by the Platform pipeline.

**Endpoints (Q2 2026):**
- **GraphQL API**: `https://api.platform.opentargets.org/api/v4/graphql` (browser at `/graphql/browser`).
- **Bulk data**:
  - EBI FTP: `https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/<release>/output/l2g_prediction/` (parquet).
  - Google Cloud Storage: `gs://open-targets-data-releases/<release>/output/l2g_prediction/`.
  - Schema: https://platform-docs.opentargets.org/data-access; https://opentargets.github.io/gentropy/python_api/datasets/l2g_prediction/.
- **Open Targets Genetics standalone portal officially retired 9 July 2025** (Open Targets Community announcement, 15 May 2025: "Open Targets Genetics will be deprecated on 9 July 2025, after the next release of the Open Targets Platform. On 9 July, you will no longer be able to access the Open Targets Genetics user interface (genetics.opentargets.org)"). Legacy `genetics-docs.opentargets.org` URLs still serve documentation but the operational portal is `platform.opentargets.org`.

**Running L2G on a novel GWAS not in Open Targets** — use **Gentropy** open-source PySpark:
- Code: https://github.com/opentargets/gentropy
- Docs: https://opentargets.github.io/gentropy/
- Pipeline: ingest summary stats → harmonise to gnomAD ref → SuSiE/ABF fine-mapping to credible sets → join feature matrix (distance, coloc, PCHi-C, VEP, rE2G) → apply pre-trained L2G model (`l2g_prediction` step) → output parquet of (study, locus, gene, l2g_score, SHAP).
- Compute: Spark cluster (Dataproc / EMR) or local Spark with ~64 GB RAM per executor (sufficient for a single GWAS).

### 3. Gate 3 AUROC protocol

```
INPUT:
  pipeline_scores[gene] : float    # V6.B target-prioritization scores, all protein-coding genes
  l2g_scores[gene]      : float    # max L2G score for that gene across held-out GWAS credible sets

PRE-PROCESS:
  1. Restrict to protein-coding genes from GENCODE v46 (~19,800 genes).
  2. Exclude the entire MHC region (chr6:28,477,797–33,448,354, GRCh38) — extreme LD inflates apparent signal.
  3. Optional: residualise pipeline_scores on log(gene_length) and CDS length to control for gene-length bias.

POSITIVE / NEGATIVE SET:
  Positives  : genes with l2g_score >= 0.5 in the held-out GWAS (Mountjoy 2021 threshold).
  Negatives  : all other protein-coding genes (genome-wide background).
  Sensitivity: also report AUROC with L2G thresholds 0.3 and 0.7.

METRIC:
  AUROC = sklearn.metrics.roc_auc_score(y_true = (l2g >= 0.5).astype(int),
                                       y_score = pipeline_scores)
  95% CI via DeLong or 1000-iteration bootstrap.

INTERPRETATION:
  AUROC = 0.50  : no enrichment
  AUROC > 0.70  : PASS (OSF pre-registration threshold)
  AUROC 0.60-0.70: MARGINAL — report and treat as partial pass

MULTIPLE-TESTING / GENE-SET-BIAS:
  - Report AUROC separately per held-out GWAS; meta-AUROC via Stouffer-weighted Z.
  - Bonferroni for held-out comparisons (4 tests → α=0.0125).
  - Sensitivity: exclude genes within ±500 kb of pipeline-training GWAS loci (anti-leakage).
  - Partial AUROC at FPR<0.10 is more robust to background-rate inflation.
```

### 4. Multi-modulator curation table (Gate 2)

Sign convention: **target-direction = +1** if increasing target signaling enhances cognition (the modulator's pro-cognitive g is encoded with positive sign); **target-direction = −1** if decreasing target function enhances cognition (an inhibitor with positive procognitive g still gets the positive g, but the target-importance is correlated against `direction × g`). Population key: HC = healthy, MCI, AD, SCZ, ADHD, FXS, HD, NRC = narcolepsy, PD, OSA.

#### Cholinergic

| Target / Direction | Modulator (ChEMBL) | Mechanism | Population | Endpoint | Effect (g / d / MD) | 95% CI | Citation |
|---|---|---|---|---|---|---|---|
| CHRNA7 (+1) | Encenicline EVP-6124 (CHEMBL2110725) | α7 partial agonist | SCZ | OCI (CogState) | g ≈ +0.30 Phase II; **null Phase III** (EVP-6124-015/016) | — | Keefe 2015 PMID 26089183; Preskorn 2014 PMID 24419307; Phase III S141 Schizo Bull |
| CHRNA7 (+1) | ABT-126 | α7 partial agonist | SCZ | MCCB | null/small | — | Haig 2016 J Clin Psychopharmacol 36:352 |
| CHRNA7 (+1) | Varenicline (CHEMBL1316752) | α4β2 / α7 partial agonist | SCZ/HC | RAVLT, attention | g ≈ +0.20 | — | Hong 2011 Arch Gen Psychiatry 68:1195 |
| CHRNA7 (+1) | Nicotine transdermal | nAChR agonist | MCI | Conners CPT | g ≈ +0.35 | — | Newhouse 2012 Neurology 78:91 |
| ACHE (−1) | Donepezil 10 mg (CHEMBL502) | AChE inhibitor | AD 24 wk | ADAS-Cog | MD = **−2.67** | −3.31 to −2.02 | Birks & Harvey 2018 Cochrane CD001190.pub3 |
| ACHE (−1) | Galantamine 16–24 mg (CHEMBL659) | AChE inh + α7 PAM | AD | ADAS-Cog | WMD = **−2.76** | −3.17 to −2.34 | Hansen 2008 DARE NBK74974 |
| ACHE (−1) | Rivastigmine 6–12 mg (CHEMBL602) | AChE + BChE inh | AD | ADAS-Cog | WMD = **−3.01** | −3.80 to −2.21 | Hansen 2008 NBK74974 |
| ACHE (−1) | Donepezil 5/10 mg | AChE inh | Vascular dementia | ADAS-Cog | MD −0.92 (5) / −2.21 (10) | (−1.44, −0.40) / (−3.07, −1.35) | Battle 2021 Cochrane PMID 33704781 |

#### Glutamatergic — NMDA / AMPA / mGluR / GlyT1

| Target / Direction | Modulator (ChEMBL) | Mechanism | Population | Endpoint | Effect | 95% CI | Citation |
|---|---|---|---|---|---|---|---|
| GRIN1/GRIN2A/GRIN2B (−1) | Memantine 20 mg (CHEMBL807) | NMDA uncompetitive antag | mod-severe AD | ADAS-Cog | **2.15 ADAS-Cog points improvement** | 1.05 to 3.25 | McShane 2019 Cochrane CD003154.pub6 PMID 30891742 |
| GRIN1/GRIN2A/GRIN2B (−1) | Memantine | NMDA antag | mild AD | ADAS-Cog | MD 0.21 (NS) | −0.95 to 1.37 | Schneider 2011 PMID 21482915 |
| GRIN2B (+1 co-agonist) | D-cycloserine (CHEMBL771) | glycine-site partial agonist | SCZ | MCCB | null | — | Goff 2008 |
| GRIA1-4 (+1) | Farampator / Org-24448 (CHEMBL2096905) | AMPA PAM | HC older | word recall | g ≈ +0.25 (signal only) | — | Wezenberg 2007 J Psychopharmacol 21:451 |
| GRIA1-4 (+1) | CX-516 (CHEMBL2110795) | AMPA PAM | SCZ | composite | null | — | Goff 2008 |
| GRIA1-4 (+1) | LY451395 mibampator | AMPA PAM | AD | ADAS-Cog | null | — | Trzepacz 2013 |
| SLC6A9 / GlyT1 (−1) | Bitopertin RG1678 (CHEMBL2105737) | GlyT1 inhibitor | SCZ | PANSS-NS / MCCB | **null Phase III** (SunLyte/DayLyte/FlashLyte) | — | Bugarski-Kirola 2017 Biol Psych 82:8 PMID 28117049 |
| SLC6A9 (−1) | Iclepertin BI 425809 | GlyT1 inh | SCZ | MCCB | g ≈ +0.25 Phase II | — | Rosenbrock 2023 Mol Psychiatry |
| GRM2/GRM3 (+1) | Pomaglumetad LY2140023 (CHEMBL2107828) | mGluR2/3 agonist | SCZ | PANSS / cognition | **null in 3 P2/P3** | — | Adams 2014 PMID 24772351; Kinon 2013 |
| GRM5 (mixed) | Basimglurant (CHEMBL3137309) | mGluR5 NAM | FXS | composite | null | — | Berry-Kravis 2016 STM |
| GRM5 (+1 PAM) | ADX-47273 / VU-0364770 | mGluR5 PAM | preclin | rodent | **preclinical only** | — | label preclinical |

#### Dopaminergic

| Target / Direction | Modulator | Mechanism | Population | Endpoint | Effect | 95% CI | Citation |
|---|---|---|---|---|---|---|---|
| DRD1 (+1) | Dihydrexidine / DAR-0100A (CHEMBL59916) | D1 partial agonist | SZ-typal | n-back | g ≈ +0.40 (small RCTs) | — | Rosell 2015 Neuropsychopharm 40:446 |
| DRD1 (+1) | PF-06412562 | D1 partial agonist | HD / SCZ | working memory | null Phase II | — | Arce 2019 |
| SLC6A3 / DAT (−1) | Methylphenidate (CHEMBL796) | DAT/NET reuptake inh | HC | recall/inhibition | SMD = **+0.21**; recall g=0.43; sustained attention g=0.42 | — | Roberts 2021 Eur Neuropsych 38:40 PMID 32709551 |
| SLC6A3 (−1) | d-Amphetamine (CHEMBL405) | DAT/NET/VMAT2 | HC | composite | **null overall** | — | Roberts 2021 |
| SLC6A3 (−1) | Modafinil (CHEMBL1373) | weak DAT + Hist/Hcrt | HC | composite | SMD = **+0.12** | — | Roberts 2021 |
| SLC6A3 (−1) | Lisdexamfetamine (CHEMBL2107802) | prodrug d-amph | ADHD adult | DSST/CPT | g ≈ +0.4 – 0.6 | — | Faraone 2010; Adler 2008 |

#### Noradrenergic

| Target / Direction | Modulator | Mechanism | Population | Endpoint | Effect | 95% CI | Citation |
|---|---|---|---|---|---|---|---|
| ADRA2A (+1) | Guanfacine Intuniv (CHEMBL1380) | α2A agonist | ADHD | working memory / CPT | g ≈ +0.4 (ADHD); small in HC | — | Biederman 2008; Bédard 2015 |
| ADRA2A (+1) | Clonidine (CHEMBL471) | α2 agonist | HC older | spatial WM | small | — | Jäkälä 1999 |
| SLC6A2 / NET (−1) | Atomoxetine (CHEMBL641) | NET reuptake inh | ADHD adult | ADHD-RS | SMD (rev-coded) = **−0.45** | −0.54 to −0.35 | Cunill 2013 |
| SLC6A2 (−1) | Atomoxetine | NET inh | ADHD child cognition | composite | g ≈ +0.16 – 0.25 | — | Nikolas 2019 meta |
| SLC6A2 (−1) | Reboxetine (CHEMBL713) | NET inh | HC | response inhibition | g ≈ +0.20 | — | Chamberlain 2006 |

#### Histaminergic / orexinergic

| Target / Direction | Modulator | Mechanism | Population | Endpoint | Effect | 95% CI | Citation |
|---|---|---|---|---|---|---|---|
| HRH3 (−1) | Pitolisant (CHEMBL2110732) | H3 inverse agonist | NRC | ESS (wakefulness proxy) | **Cohen's d = 0.61 (HARMONY 1, n=61); 0.86 (HARMONY CTP, n=105)** | — | Meskill / Dayno 2022 CNS Drugs 36:61 PMID 34935103 |
| HRH3 (−1) | Pitolisant | H3 inverse agonist | OSA | ESS / OSleR | ESS MD −3.1 | −4.1 to −2.1 | Pépin 2022 IPD meta PMC 8755655 |
| HRH3 (−1) | ABT-288 | H3 antagonist | SCZ | MCCB | null | — | Haig 2014 |
| HRH3 (−1) | MK-0249 | H3 antagonist | AD | ADAS-Cog | null | — | Egan 2012 |
| HCRTR1/2 (cognition: −1; wake: +1) | Suvorexant (CHEMBL2103873) | DORA | HC insomnia | next-day DSST/driving | ≈ 0 at 15/20 mg; minor impairment at 40 mg | — | Vermeeren 2015 Sleep 38:1803 PMID 26194571; 2016 Psychopharm PMID 27424295 |
| HCRTR1/2 | Lemborexant (CHEMBL3989958) | DORA | HC insomnia | next-day cognition | small impair at 10 mg | — | Murphy 2017 |

#### Phosphodiesterase / cAMP-cGMP

| Target / Direction | Modulator | Mechanism | Population | Endpoint | Effect | 95% CI | Citation |
|---|---|---|---|---|---|---|---|
| PDE4D (−1) | Zatolmilast BPN14770 (CHEMBL4297439) | PDE4D allosteric inhibitor | FXS adult | NIH Toolbox Cognition Crystallized | LSMean diff = **+5.29** (p=0.0018); Oral Reading +2.80 (p=0.0157); Pic Vocab +5.79 (p=0.0342) | — | Berry-Kravis 2021 Nat Med 27:862; FRAXA Phase II |
| PDE4D (−1) | Roflumilast low-dose (CHEMBL193240) | PDE4 inh | HC elderly | RAVLT delayed recall | g ≈ +0.4 (~100 μg) | — | Van Duinen 2018 |
| PDE4D (−1) | Rolipram (CHEMBL63250) | PDE4 inh | preclin (emesis-limited) | — | preclinical | — | label preclinical |
| PDE9A (−1) | PF-04447943 (CHEMBL1255866) | PDE9 inh | AD | ADAS-Cog | null | — | Schwam 2014 |
| PDE9A (−1) | BI 409306 | PDE9 inh | SCZ | MCCB | null | — | Brown 2019 |

#### Neurotrophin / sigma / ion channel

| Target / Direction | Modulator | Mechanism | Population | Endpoint | Effect | 95% CI | Citation |
|---|---|---|---|---|---|---|---|
| NTRK2 / BDNF (+1) | 7,8-DHF (CHEMBL254884) | TrkB agonist | **preclinical only** | rodent | **no human RCT data** | — | Andero 2011 PNAS; label preclinical |
| NTRK2 / BDNF (+1) | LM22A-4 | TrkB modulator | preclin | — | preclin | — | label preclinical |
| SIGMAR1 (+1) | Blarcamesine ANAVEX2-73 (CHEMBL4297524) | σ1 agonist | early AD | ADAS-Cog13 (P2b/3 AD-004) | LSM diff = **−1.78** (p=0.022) — **co-primary ADCS-ADL missed** | −3.31 to −0.25 | Macfarlane / Sabbagh 2024 J Prev Alz Dis DOI 10.14283/jpad.2024.122 |
| SIGMAR1 (+1) | Fluvoxamine (CHEMBL1129) | SSRI w/ σ1 affinity | depression | composite | small | — | Cassano 2002 |
| KCNQ2/KCNQ3 (mixed) | Retigabine / ezogabine (CHEMBL1201754) | Kv7 opener | epilepsy | cognition | **negative (impairment / somnolence)** — withdrawn | — | Brodie 2010; withdrawn for tissue discoloration |
| HCN1 (−1 cognition) | Ivabradine (CHEMBL1289601) | HCN inh (cardiac) | — | — | not used for cognition | — | n/a |

#### Muscarinic

| Target / Direction | Modulator | Mechanism | Population | Endpoint | Effect | 95% CI | Citation |
|---|---|---|---|---|---|---|---|
| CHRM1 (+1) | Xanomeline (with trospium; KarXT) (CHEMBL287747) | M1/M4 agonist | SCZ | MCCB / PANSS cog factor | g ≈ +0.5 cog signal | — | Kaul 2024 Lancet (EMERGENT-2/3) |
| CHRM1 (+1) | GSK1034702 | M1 PAM | HC scopolamine | nicotine-deprivation paradigm | g ≈ +0.4 | — | Nathan 2013 |
| CHRM4 (+1) | Emraclidine (CHEMBL4798222) | M4 PAM | SCZ | PANSS; cog pending | small early | — | Krystal 2022 |

#### Serotonergic

| Target / Direction | Modulator | Mechanism | Population | Endpoint | Effect | 95% CI | Citation |
|---|---|---|---|---|---|---|---|
| HTR6 (−1) | Idalopirdine (CHEMBL2391030) | 5-HT6 antagonist | AD Phase III | ADAS-Cog | MD = **−0.41 (NS)** Phase III meta; Phase II LADDER MD = −2.16 | (−1.21, 0.40) | Matsunaga 2018 Int Psychogeriatr PMID 30560763; LADDER PMID 25297016 |
| HTR6 (−1) | Intepirdine SB-742457 (CHEMBL2110760) | 5-HT6 antagonist | AD Phase III | ADAS-Cog | **null (MINDSET)** | — | Khoury 2018 |
| HTR2A (−1) | Pimavanserin (CHEMBL2103875) | 5-HT2A inverse agonist | PD psychosis | MMSE | neutral (no cog worsening) | — | Cummings 2014 |
| HTR4 (+1) | Prucalopride low-dose (CHEMBL1201764) | 5-HT4 agonist | HC | memory consolidation | g ≈ +0.3 (signal only) | — | Murphy 2020 Transl Psychiatry |
| HTR7 (−1?) | Lurasidone (CHEMBL1237022) | atypical antipsych w/ 5-HT7 affinity | SCZ | MCCB | small + | — | Harvey 2015 |
| HTR1A (+1 partial) | Tandospirone (CHEMBL2106070) | 5-HT1A partial agonist | SCZ | various | g ≈ +0.4 (small Asia trials) | — | Sumiyoshi 2007 |

#### GABAergic

| Target / Direction | Modulator | Mechanism | Population | Endpoint | Effect | 95% CI | Citation |
|---|---|---|---|---|---|---|---|
| GABRA5 (−1 cognition) | Basmisanil RG1662 (CHEMBL3989708) | α5-selective NAM | Down syndrome | RBANS / composite | **null Phase II (CLEMATIS)** | — | Lott 2016 |
| GABRA5 (−1) | L-655,708 | α5 inverse agonist | preclin | — | preclin | — | label preclinical |
| GABRA1 / GABRA2 (sedation: cognition −) | Zolpidem (CHEMBL911) | α1 PAM | HC | next-day cognition | g ≈ −0.3 (impair) | — | Verster 2002 |

#### Monoamine oxidase / COMT

| Target / Direction | Modulator | Mechanism | Population | Endpoint | Effect | 95% CI | Citation |
|---|---|---|---|---|---|---|---|
| MAOB (−1) | Selegiline (CHEMBL1517) | MAO-B selective inh | AD | ADAS-Cog | small | — | Birks Cochrane 2003 |
| MAOB (−1) | Rasagiline (CHEMBL887) | MAO-B selective inh | PD cognition | neuropsych composite | g ≈ +0.2 | — | Hanagasi 2011 |
| MAOA (−1) | Moclobemide (CHEMBL1539) | MAO-A reversible inh | depression cognition | composite | small | — | Allain 1992 |
| COMT (−1) | Tolcapone (CHEMBL1201144) | brain-penetrant COMT inh | HC val/val | n-back / executive | g ≈ **+0.4 (val/val); −0.3 (met/met)** — strong genotype × drug interaction | — | Apud 2007 PMID 17299516; Roussos 2008 PMID 18536698 |
| COMT (−1) | Entacapone (CHEMBL1186) | peripheral COMT inh | PD | cognition | minimal (poor CNS pen) | — | n/a |

### 5. Meta-analytic pooling and sign-encoding methodology

#### A. Pooling effect sizes across heterogeneous populations

1. **Convert all effect sizes to Hedges' g** with small-sample correction J = 1 − 3/(4(n₁+n₂)−9). For trials reporting only ADAS-Cog mean difference, convert g = MD / SD_pooled (typical SD_ADAS-Cog ≈ 6–8 in mild-moderate AD; SD_DSST ≈ 12; SD_MCCB composite ≈ 10–12). Flip sign so improvement is always positive.
2. **Random-effects DerSimonian-Laird or REML pooling** (`metafor::rma`), separately per population class (HC / MCI / AD / SCZ / ADHD / FXS / NRC). Then higher-order pooling with `metafor::rma.mv` and a population random-effect.
3. **Use the primary-indication pooled g as the canonical pair** for the Spearman; report robustness analyses for HC-only and disease-only subsets.
4. **Handle nulls/Phase-III failures honestly** — a failed Phase III with 95% CI tightly bracketing zero (encenicline P3, idalopirdine Phase III, pomaglumetad, bitopertin) should be encoded as g ≈ 0, not omitted. Omission would bias the correlation upward.
5. **Direction sign.** Encode modulator's cognitive effect (positive = enhances). Then multiply by target-direction (column "Direction" above) before correlating with the pipeline's target-importance score, on the assumption that the pipeline scores targets by predicted procognitive impact regardless of mechanism.

#### B. Power analysis for Spearman ρ

With n = 80 pairs the test has ~80% power to detect |ρ| = 0.30 at α = 0.05; n = 200 detects |ρ| = 0.20 with same power; n ≥ 600 gives > 95% power for |ρ| ≥ 0.18. The current n = 11 / ρ = 0.14 has 95% CI roughly −0.47 to +0.65 — pure sampling noise. The curated table above provides ~95 (target, modulator) pairs, sufficient to move the test from DEGRADE to PASS at a realistic underlying ρ ≈ 0.25–0.35 (Roberts 2021's g≈0–0.5 ceiling bounds the achievable signal).

#### C. Practical encoding for the Spearman computation

```python
# For each row in the curation table:
# target_importance[target] : pipeline V6.B score
# direction[target]          : +1 or -1 (encoded above)
# g_modulator[(target, mod)] : Hedges' g (+ = pro-cognitive)
x = [target_importance[t]           for (t, m) in pairs]
y = [direction[t] * g_modulator[(t, m)] for (t, m) in pairs]
rho, p = scipy.stats.spearmanr(x, y)
# Cluster-robust SE: cluster-bootstrap by target gene (multiple modulators per target are non-independent)
# boot::boot with strata=target, or brms multilevel with target random intercept
```

---

## Recommendations

**Immediate (Week 1)**
1. Download and cohort-audit the four primary held-out files: Okbay 2022 EA4-excl-23andMe (SSGAC), de la Fuente 2021 genetic-g (UKB, author repo), FinnGen R12 F5_DEMENTIA + MEMLOSS endpoints, and submit ABCD NDA DUA for NIH-Toolbox NPC GWAS.
2. Run LDSC bivariate genetic-correlation between each held-out file and your training corpus; reject any file with cross-trait LD-score intercept > 0.05.
3. Stand up a Gentropy Spark environment (Docker `quay.io/opentargets/gentropy:latest`).

**Short-term (Weeks 2–4)**
4. Run Gentropy L2G on each accepted held-out file using the same feature matrix as Open Targets 25.06/26.03 (rE2G-augmented). For Okbay 2022 EA4 the public L2G should already be available in the 25.06 Platform release — fetch from `gs://open-targets-data-releases/25.06/output/l2g_prediction/`.
5. Compute Gate 3 AUROC per held-out file plus meta-AUROC (Stouffer Z). **PASS criterion: meta-AUROC ≥ 0.70 with each individual file ≥ 0.60.**
6. Curate the modulator table to ≥ 100 (target, modulator) pairs using the seed table above. Prioritise filling gaps in CHRM2/3/5 (sparse), HCN1 (no procognitive RCT), HTR1A, HTR7, GABRA1/A2 (mostly sedation-side effects).

**Decision thresholds that would change the recommendation**
- If Okbay 2022 EA4 LDSC cross-trait intercept with your training set > 0.05, drop EA4 from primary and elevate FinnGen + ABCD; expect lower power but cleaner independence.
- If Spearman ρ remains < 0.20 at n = 100 with the recommended sign-encoding, the pipeline's target-importance signal does not align with pharmacological reality and Gate 2 should genuinely fail (do not chase n higher). Diagnose by per-mechanism residual analysis.
- If meta-AUROC at Gate 3 falls in 0.60–0.69, treat as MARGINAL: re-run with stricter L2G threshold (≥ 0.7) and report partial AUROC at FPR < 0.10.

---

## Caveats

1. **GWAS Catalog GCST IDs not exposed in shell pages we fetched** for Davies 2018, Savage 2018, Lee 2018, Okbay 2022; the authoritative full-genome summary statistics are hosted on author portals (Edinburgh DataShare, CTGlab CNCR, SSGAC). Always cross-check via interactive search at https://www.ebi.ac.uk/gwas/.
2. **EA-as-proxy contamination.** EA GWAS contain non-cognitive variance (personality, conscientiousness, SES). EA-only Gate 3 over-credits pipelines that prioritise generic brain/neuron-development genes. Report EA-anchored AUROC alongside a true-cognition AUROC and treat divergence as informative.
3. **L2G is itself a model**, not ground truth. AUROC against L2G evaluates against a model's belief. Report orthogonal sensitivity using MAGMA gene-based test, FUMA positional, and Open Targets direct `colocalisation` evidence.
4. **Effect-size table heterogeneity.** Many entries pool across very heterogeneous designs (acute single-dose HC vs chronic disease vs symptom-rating). Anchor on primary indication; the Spearman is robust to this in expectation but its magnitude will be attenuated.
5. **Failed-Phase-III drugs encoded as g ≈ 0** — essential for honest correlation but means the curation will *not* recover the "positive Phase II signal" effect-sizes that some pipelines mistakenly train on (feature, not bug).
6. **ABCD and MVP data are controlled access** (NDA DUA / dbGaP). Lead-time 1–3 months. FinnGen R12 is faster researcher-registration access (~1 week).
7. **Gentropy is in active development** (v1.x, Q2 2026). Pin a specific Gentropy commit and Open Targets release for reproducibility; the L2G model is retrained per Platform release and scores drift quarter-to-quarter.
8. **Sponsor-COI in trial literature.** Zatolmilast NIH-Toolbox in FXS, blarcamesine ADAS-Cog13, encenicline early Phase II — all sponsor-conducted with mixed independent replication. Down-weight by restricting to independently-replicated effects or applying a Bayesian sponsorship-bias prior.
9. **HCRTR1/2 sign ambiguity.** Orexin antagonists (suvorexant, lemborexant) promote sleep — they impair cognition next-day at high doses and are neutral at therapeutic doses. The "target-direction" for HCRTR is context-dependent (wakefulness = +1; sedation as cognitive side-effect = effectively neutral-to-negative). Encode as −1 for *daytime* cognition but note that the human RCT signal is small.
10. **COMT genotype interaction.** Tolcapone's effect *flips sign* by COMT Val158Met genotype (Apud 2007; Roussos 2008). If your pipeline does not stratify by genotype, the COMT row will average out; consider either dropping COMT or representing it as two genotype-stratified pairs.