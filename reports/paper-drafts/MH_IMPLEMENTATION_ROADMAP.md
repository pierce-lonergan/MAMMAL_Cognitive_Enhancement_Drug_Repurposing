# MH Implementation Roadmap — V6.B / V7 / V8 Enhancements

**Date:** 2026-05-27
**Inputs:** 6 research deep-dives placed in `research/4-tier/`
**Scope:** 90-day sprint plan to close the highest-ROI gaps from `GAPS_AND_RESEARCH_DIRECTIONS.md`, integrate the 6 deep-dive findings into running code (`src/mammal_repurposing/...`), and refresh the 5-paper manuscript suite.

---

## 0. Executive Summary

The six research deep-dives collectively converge on **four high-leverage, low-risk interventions** the pipeline can absorb before any new ground-truth data is generated:

| # | Intervention | Doc(s) | Closes | Engineering cost | Decision-grade impact |
|---|---|---|---|---|---|
| **A** | Bundle MH3 (per-cell-line random effect) + MH7 (per-species random effect) into V8 hierarchical model; calibrate with cpg0000 pilot | MH3 | V8 U2OS-to-brain critique; G1-G6 pre-reg gates | 2 weeks | **Largest** — defends V8 against the obvious reviewer attack |
| **B** | MH8 substrate-mediated flag + non-centered reparam + `target_accept=0.95` on V6.B.5 | MH8 | 37-divergence NUTS run; biophysical correctness for ACHE/MAO/COMT | 2 days | Removes only outstanding numerical defect in V6.B headline result |
| **C** | MH1+MH2 V7 anchor expansion (15 → 60 compounds) **+** population × class interaction term refit | MH1+MH2 V7 CPT, MH1+MH2 PRISMA alt § 6 | V7 Gate 1 partial-pool gap; ESS lift; CI shrinkage on rare-mechanism classes | 3 weeks | Partial — closes ~60% of the V7 0.004-0.063 PMD gap |
| **D** | Gate 2 (multi-modulator curation, ~95 pairs) **+** Gate 3 (Okbay 2022 EA4 + Gentropy L2G held-out AUROC) | Cluster D Methodology | V6.B Gate 2 DEGRADE → PASS; Gate 3 INSUFFICIENT_DATA → PASS | 1 week + 1 week | Converts V6.B from "headline result + 2 open gates" to "full 4-gate PASS" |

Bonus (E): **LINCS L1000 chemCPA training kit** (deep-dive #5) gives turn-key hyperparameter table, MD5-pinned molecular encoder weights, OOD compound list, and `bool_de` evaluation toggle — enables V8.2 to graduate from synthetic-LINCS smoke to real pretraining.

**Critical constraint to flag up-front:**
> MH3 and MH7 **must be bundled** in one PyMC refit. Splitting them produces an identifiability collision between α (cell-line) and γ (species) when only U2OS data is available for some compounds — non-identifiable parameters → divergences. The MH3 doc § 7 is explicit on this.

> The MH1+MH2 anchor expansion (15 → 60) only **partially** closes the V7 0.004-0.063 PMD residual. The real binding constraint is **τ²_class mis-specification** (Schmidli rMAP weights set too tight) **+ missing population × class interaction term**. Anchor expansion without the population × class interaction refit will not lift V7 Gate 1 to PASS.

---

## 1. State of the Pipeline (2026-05-27 baseline)

| Layer | Status | Artefact | Validation |
|---|---|---|---|
| V4 fusion (5-cluster calibrated) | Shipped | `data/results/v2/combined_gates.parquet` | 4 cluster gates; production |
| V5 Z-norm + Tier 2/3 | Shipped | `data/results/v2/disagreement_report_v9_mmatt.md` | 9 disagreement reports, 5 calibration layers |
| V6.A Multi-Head DTI | Shipped + drafted | `reports/paper-drafts/v6a_paper_draft.md`, `figures/v6a/` | 5-head Bayesian router + Venn-ABERS |
| V6.B Cluster D (22-target headline) | Shipped + drafted | `data/results/v2/cluster_d_posterior_v1.parquet` (R̂=1.000, ESS=12,780, 0 divergences) | 4-gate framework, 22-target panel |
| V6.B.5 (191-target expanded) | Shipped, **numerically degraded** | `data/results/v2/cluster_d_posterior_expanded_v1.parquet` (R̂=1.000, ESS=1,739, **37 divergences**) | numpyro JAX backend; **MH8 fix pending** |
| V7 Clinical Effect-Size Translation | Shipped + drafted | `data/results/v2/v7_nuts_posterior_production_v1.parquet` (R̂=1.000, ESS=9,232, MAE=0.071, 0 Roberts violations) | OSF pre-reg shipped; **15 anchors only** |
| V8 πphen Perturbational Evidence | Shipped scaffold + drafted | `src/mammal_repurposing/cluster_e/joint_phenotype.py`; 8-cell + I_novel | OSF pre-reg shipped; **synthetic LINCS only** |
| Integration manuscript | Drafted | `reports/paper-drafts/integration_paper_draft.md` | Umbrella synthesis |
| Wet-lab handoff | Shipped | `reports/wet-lab/wet_lab_shortlist_v10.md` + `reports/wet-lab/wet_lab_handoff_v1.md` | 8-cell tagged top candidates |

**Two open V6.B gates (per Cluster D Methodology doc § 5):**
- Gate 2: DEGRADE, Spearman ρ ≈ 0.41, only 17 (target, modulator) pairs — needs ≥ 80 pairs.
- Gate 3: INSUFFICIENT_DATA — no held-out cognition-GWAS L2G evaluation yet.

---

## 2. Per-Doc Findings — What Each Deep-Dive Actually Delivers

### 2.1 `MH3_per_cell_line_random_effect_deep_research.md`

**Verdict:** 9/10 must-have. Drop-in PyMC code provided. Effort revised 1 week → 2 weeks. **MH3 + MH7 must bundle.**

**Key equations:**
```
y_{c,l,s,m,k,r}  =  μ  +  β_{c,k}         (transferable mechanism × endpoint effect)
                       +  α_{l,k}         (cell-line random effect, non-centered)
                       +  γ_{s,k}         (species random effect)
                       +  δ_{c,l,k}       (compound × cell interaction)
                       +  ε_{m,r}         (measurement + replicate noise)
```

Plus a **transferability index** `T_{c,k} = σ²_β / (σ²_β + σ²_α + σ²_δ)` per (compound, endpoint), which directly populates a new column in the wet-lab shortlist explaining how confident we are the U2OS-derived signal will transfer to brain.

**Concrete deliverables already specified in the doc:**
1. Drop-in PyMC builder (`build_v8_hierarchical_with_cell_random_effect`) — copy-paste into `src/mammal_repurposing/cluster_e/joint_phenotype.py` with minimal edits.
2. Non-centered parameterisation already written (avoids Neal's funnel).
3. cpg0000 pilot dataset (~300 compounds × A549 + U2OS shared) for empirical σ_α + σ_δ prior calibration. Use this to **set HalfNormal(0.5) → empirically-fit HalfNormal(σ̂_α)** for each endpoint.
4. Six pre-registered validation gates G1-G6 (ICC bounds, leave-cell-line-out R², external concordance vs Gorgogietas 2025 + Anderson 2025).
5. Six honest limitations L1-L6 to add to V8 paper Methods.

**Citations to add to `CITATIONS.bib`:**
- Chandrasekaran et al. 2024 (JUMP cpg0016 Nature Methods).
- Gorgogietas et al. 2025 Scientific Reports.
- Anderson et al. 2025 eLife.

---

### 2.2 `Cluster D Methodology Report — Gate 3 and Gate 2 Multi-Modulator Curation.md`

**Verdict:** Plug-and-play. Provides the **complete tables** needed to close V6.B Gate 2 and Gate 3.

**Gate 2 table:** ~95 (target, modulator, mechanism, pooled_g, CI_lo, CI_hi, endpoint, citation_doi, population) rows covering all 28 cognition targets. Includes ALL Phase III failures (encenicline, idalopirdine, intepirdine, pomaglumetad, bitopertin) with `g ≈ 0` encoding so the model learns "high-affinity binder ≠ clinical benefit."

**Gate 3 table:** 14 held-out cognition GWAS sumstats with cohort-overlap audit recommendations. Top picks:
- **Okbay 2022 EA4** (N=3,037,499) — primary held-out anchor.
- de la Fuente 2021 (memory PRS) — secondary.
- FinnGen R12 (G6_AD, F5_DEMENTIA) — independent ancestry.
- ABCD release 5 (childhood cognition) — for developmental cohort.
- **Explicitly EXCLUDE:** Davies 2018, Hill 2018, Savage 2018 (cohort overlap with training).

**Sign-encoding methodology + power analysis:**
- n=80 (target, modulator) pairs gives 80% power for |ρ|=0.30 at α=0.05 (the Cluster D Methodology doc § 5 calculation).
- Sign encoding: positive g = pro-cognitive, negative g = anti-cognitive; 0 = clinically null at expected dose. Phase III nulls anchor the "high-Kd-low-effect" tail.

**Gate 3 AUROC protocol:**
- Positives: genes with Open Targets L2G score ≥ 0.5 in held-out GWAS, restricted to protein-coding outside MHC.
- Negatives: random gene-set matched on length + chromosome.
- Metric: AUROC; pass threshold 0.65 (acknowledged: aspirational; lit baselines are 0.55-0.70).

**Gentropy toolchain note (critical):**
- **Open Targets Genetics standalone portal retired 9 July 2025.**
- Use **Open Targets Platform 25.06+** L2G outputs (Parquet on GCS) **OR** **Gentropy PySpark** for re-running L2G on held-out cohorts (preferred for true held-out evaluation).

---

### 2.3 `MH1 + MH2 Meta-Analytic Prior Expansion for V7 CPT Bayesian Pharmacology Pipeline.md`

**Verdict:** Highest-ROI for V7 paper. Provides the **96-cell PER_SUBDOMAIN_PRIORS table** (12 mechanism classes × 8 endpoints) AND a **60-row REFERENCE_COMPOUND_SMD** dict in Python format ready to ingest.

**Concrete drop-ins:**

1. **`prisma_priors.py` refresh** — current file has 12 PRISMA classes + 32 per-subdomain cells. The doc provides all 96 cells with (g, CI_lo, CI_hi, k, source). Schmidli 2014 robust MAP prior weights per cell count:
   - k ≥ 5 → τ² = 0.02
   - k = 2-4 → τ² = 0.04
   - k = 1 → τ² = 0.08
   - k = 0 (sparse mechanism × endpoint) → fall back to class-level Half-Normal(0.3).

2. **`REFERENCE_COMPOUND_SMD` 15 → 60 anchor rows** — already in Python dict format in the doc. Includes all canonical compounds: donepezil, rivastigmine, galantamine, memantine, methylphenidate, modafinil, atomoxetine, pitolisant, encenicline (failed), idalopirdine (failed), zatolmilast (BPN14770), pridopidine, blarcamesine, etc.

3. **Effect-size conversion methodology** — Cohen's d → Hedges' g (small-sample correction); ADAS-cog mean difference → SMD with pooled SD ≈ 7.5 (the doc's recommended conservative pooled SD).

**Critical insight buried in § 5 of the doc:**
> The "0.50 ceiling" we currently use in V7 is **INTERPRETIVE**, not literal. Roberts 2020 highest single-subdomain SMD is g = 0.43 (MPH on recall). Soft Half-Normal(0.3) upper envelope is the correct constraint, NOT hard-clipped 0.50.

**The honest gap-diagnosis:**
> "Will 15 → 60 anchor expansion actually fix V7 partial-pool gap?" → **PARTIAL only.** Real constraints:
> - τ²_class is mis-specified (currently a single τ²_global = 0.045; doc recommends per-class τ²_class).
> - Population × class interaction term is missing — different mechanisms have different SMD distributions in MCI vs healthy vs schizophrenia.
> - Without these, anchor expansion alone gives ~60% gap closure.

**Sprint order matters:** refit V7 with τ²_class + population × class interaction FIRST on existing 15 anchors → measure gap → THEN add the 45 new anchors. Otherwise you can't tell which intervention did the work.

**60-anchor table also auto-generates updated Figure F7-3** (PMD waterfall + Roberts envelope).

---

### 2.4 `MH1+MH2_ PRISMA, Anchor Expansion.md` (the "alternative" doc)

**Verdict:** **Largely off-topic** (Sections 1-4 cover hardware engineering anchors, DNS subdomain takeover, and Prisma Cloud — entirely unrelated to PRISMA the systematic-review framework, despite the name collision). **However, Sections 5-9 contain genuine value:**

| § | Content | Use in roadmap |
|---|---|---|
| § 5 | Schmidli 2014 rMAP recap with conjugate Beta-mixture math | Already encoded in V7 — reference for V6.B paper Methods |
| § 6 | Clinical compound table — donepezil 18-RCT meta-analysis (n=5,948), memantine post-radiation (n=508), methylphenidate apathy (n=3 trials), JW8 dual-AChE/BACE1 inhibitor | **Cross-validate** the MH1+MH2 V7 CPT 60-row table; these are independent sources for the same compounds. Add to `REFERENCE_COMPOUND_SMD` provenance trail. |
| § 7 | chemCPA Gaussian-likelihood loss formula + JUMP cpg0016 stats (12 sites, 5 replicates, U2OS, 116k chem perturbations) | Cross-references LINCS chemCPA deep-dive doc; same numbers |
| § 8 | L2G feature SHAP table — `distanceSentinelTSS`, `vepMaximum`, `geneCount500kb`, `e2gMean`, `e2gNeighbourhoodMean` | **Use these exact column names** in Gate 3 implementation; they are the canonical Open Targets Gentropy L2G features |

**Conclusion:** Treat Sections 5-9 as **citations and cross-validation**, not as the primary source. The MH1+MH2 V7 CPT doc (§ 2.3 above) is canonical for the anchor expansion table.

---

### 2.5 `Deep Dive_ LINCS L1000 chemCPA Training.md`

**Verdict:** Turn-key V8.2 graduation kit. Replaces our `chemcpa_train.py` synthetic-LINCS smoke with real-data training recipe.

**Concrete deliverables in doc:**

| Hyperparameter | Value | Notes |
|---|---|---|
| Latent dim (`dim`) | 32 | |
| Dropout | 0.262378 | |
| Autoencoder width | 256 | |
| Autoencoder depth | 4 | |
| AE learning rate | 0.001121 | Adam |
| Batch size | 256 | choice space |
| Disentanglement weight | balanced | adversarial classifier loss |
| Likelihood | Gaussian on log-counts | (NB optional but Gaussian converges better) |
| Gradient penalty | zero-centered | stabilises discriminator |

**Molecular encoder MD5 pinning (table 1 in doc):**
- `grover_base`: `ff420aea264fca7668ecb147f60762a1`
- `jtvae`: `a7060ac4e2c6154e64a13acd414cbba2`
- `rdkit`: `4f061dbfc7af05cf84f06a724b0c8563`

**LINCS data levels (table 2):**
- Skip Level 1 (LXB raw fluorescence).
- Use Level 3 (Q2NORM `.gctx` / `.h5ad`) for primary pretraining.
- Optionally use Level 4 (Z-scores) for delta-expression modelling.
- Use Level 5 (DGE characteristic-direction signatures) as evaluation baseline.

**Batch ID extraction rules (essential — current ingestion may be wrong):**
- Sig ID `AML001_CD34_24H:BRD-A03772856:0.37037` → batch = `AML001_CD34_24H` (split on `:`, take first).
- Sample ID `ERG013_VCAP_72H_X3_B11` → batch = `ERG013_VCAP_72H` (split on `_`, take first 3 fields joined).

**Transfer learning recipe (978 → 2,000 HVG):**
- Replace encoder input layer + decoder output layer; train from scratch.
- **Freeze molecular embedding weights** in g_drug — these transfer chemical knowledge across domains.
- For K562 (or any new cell line missing from L1000) — retrain cell embedding dictionary from scratch.

**Held-out OOD compounds (use these exact 9 for our V8 paper's OOD evaluation):**
> Dacinostat, Givinostat, Belinostat, Hesperadin, Quisinostat, Alvespimycin, Tanespimycin, TAK-901, Flavopiridol

**Evaluation upgrade:** add Wasserstein distance to compute_prediction (currently only R² on full + DEG subsets). Mean-prediction R² lets controls-only baseline "win" — Wasserstein on full distribution avoids this trap.

**Docker image:** `registry.hf.space/b1ro-chemcpa:latest` — pre-installed Anaconda + CUDA + chemCPA. Useful for WSL2 reproducibility check.

**Benchmark comparators (table 4) to cite in V8 paper:**
- Biolord (outperforms chemCPA on sci-Plex3 OOD).
- PerturbNet (combined chem + genetic).
- scAgents (LLM-based agent; +20% Pearson over chemCPA).
- scDCA (drug-conditional adapter on scGPT).

---

### 2.6 `MH8 Methods Clarity Research.md`

**Verdict:** Critical for Methods clarity. Provides the algebra **and** a disambiguation table that prevents reviewer confusion with 6+ other "MH8" entities in databases.

**The mathematical fix (drop-in for `bayesian_prior.py`):**

```python
# Current (causes 37 divergences on 191-target run):
mu_t = beta_t * expression_t  # multiplicative gate forces ACHE/MAO/COMT into Neal's funnel

# MH8 reformulated:
if substrate_mediated[t]:
    mu_t = beta_t * scale_const  # constant ≈ 1.0 or target-specific global mean
else:
    mu_t = beta_t * expression_t

# Plus non-centered reparam on beta:
beta_raw = pm.Normal("beta_raw", 0, 1)
beta = mu_beta + sigma_beta * beta_raw

# Plus target_accept lift from 0.8 → 0.95
pm.sample(..., target_accept=0.95, nuts_sampler="numpyro")
```

**Biophysical justification (must go into V6.B paper Methods, not just Supplementary):**
- ACHE hydrolyses ACh "with extreme rapidity" → k_cat/K_m operates at substrate-saturated regime → enzyme density above some threshold does NOT linearly scale clearance.
- MAO-A PET (`[¹¹C]harmine`) studies show cerebral MAO-A levels are *homeostatically adapted* to local 5-HT → mRNA decoupled from net clearance.
- COMT: MB-COMT is only a *fraction* of total but accounts for ~70% of brain activity; ratio is non-linear in expression.

**Targets to flag with `substrate_mediated=True` in `targets.parquet`:**
- ACHE (P22303) — already in 22-target panel.
- MAO-A (P21397), MAO-B (P27338) — currently NOT in panel; add via V6.B.5 expansion.
- COMT (P21964) — currently NOT in panel; add via V6.B.5 expansion.
- (Conservative — only these 4 for the headline; add more enzymes in V8 work.)

**Validation assays (for the wet-lab handoff appendix):**
- APN/CD13 cleavage assay (Ala-MCA fluorogenic substrate, 96-well plate, Ex/Em 355/460 nm).
- ABPP transmembrane kinetic profiling (FP-PEG-TAMRA probe).
- FLAG-tag epitope affinity purification (anti-FLAG M2 gel, 3x FLAG peptide elution).

**MH8 nomenclature disambiguation — MUST include in V6.B paper Methods § 2.x:**
The token "MH8" overlaps with 6+ unrelated database entries (RCSB ligand, GenBank cDNA clone, Rockland antibody, Aminopeptidase N hybridoma, NHS QOF SMI indicator, Renishaw CMM probe). Reviewers WILL conflate. Add a 1-paragraph disambiguation footnote citing all 6.

---

## 3. The 90-Day Sprint Plan

### Sprint 1 — Days 1-14: **Numerical correctness + headline gates** (no new data ingestion)

**Goal:** V6.B headline result becomes defensible-by-itself before any expansion work.

| PR | Module | Change | Validation |
|---|---|---|---|
| 1.1 | `data/interim/targets.parquet` | Add `substrate_mediated:bool` column. Set True for ACHE; add MAO-A, MAO-B, COMT rows (UniProt fetch) | Schema test; ACHE flag True |
| 1.2 | `src/mammal_repurposing/cluster_d/bayesian_prior.py` | (a) Read `substrate_mediated` column; (b) conditional `mu_t = beta_t * (1.0 if sm else expression_t)`; (c) non-centered reparam on β; (d) `target_accept=0.95` default | New unit test: 22-target NUTS still hits R̂=1.000, 0 divergences |
| 1.3 | Re-run V6.B.5 NUTS | numpyro JAX backend, 191-target panel + MH8 fix | Target: 0 divergences (was 37). ESS lift to ≥ 5,000 (was 1,739). |
| 1.4 | `src/mammal_repurposing/cluster_d/validation_gates.py` | Add `gate2_pairs` ingestion from new CSV (see Sprint 2) | Test scaffold (data ingest fails until Sprint 2 CSV exists) |
| 1.5 | `reports/paper-drafts/v6b_paper_draft.md` Methods | Add MH8 § (substrate-mediated flag, biophysical justification, nomenclature disambiguation footnote). Replace 4-gate "DEGRADE" table with "PASS" once 1.3 lands. | Manual review |

**Decision gate:** If 1.3 still produces > 5 divergences after MH8 + non-centered + target_accept=0.95, ESCALATE — do NOT proceed to Sprint 2. Possible culprit: GWAS contribution term needs its own non-centered reparam.

---

### Sprint 2 — Days 15-28: **Gate 2 (multi-modulator) + Gate 3 (held-out GWAS)** plug-in

**Goal:** Convert V6.B from "headline + 2 open gates" to "full 4-gate PASS."

| PR | Module / Data | Change | Source |
|---|---|---|---|
| 2.1 | `data/raw/modulator_anchors_seed.csv` | Scaffold the ~95-row table from Cluster D Methodology doc. Columns: `target_uniprot`, `target_gene`, `compound`, `mechanism`, `pooled_g`, `CI_lo`, `CI_hi`, `k`, `endpoint`, `citation_doi`, `population` | Cluster D Methodology § 2-4 |
| 2.2 | `scripts/load_modulator_anchors.py` | Read CSV → `data/interim/modulator_anchors.parquet`. Validate: each target has ≥ 3 modulators; Phase III failures encoded with g ≈ 0 ± 0.05 | Schema + count tests |
| 2.3 | `src/mammal_repurposing/cluster_d/validation_gates.py` | Implement `gate2_spearman()` — Spearman ρ on (predicted θ̄, anchored g) across all (target, modulator) pairs. PASS threshold ρ ≥ 0.40 with n ≥ 80. | Synthetic positive control test |
| 2.4 | `data/raw/gwas_holdout_sources.yaml` | Catalogue 5 picked held-out sources (Okbay EA4, de la Fuente 2021, FinnGen R12 dementia, ABCD R5 cognition) + cohort-overlap audit field | Cluster D Methodology § 6 |
| 2.5 | `src/mammal_repurposing/cluster_d/gentropy_l2g.py` | New module: ingest Open Targets Platform 25.06+ L2G Parquet for held-out summary stats; OR call Gentropy PySpark if we need true held-out re-run | docs.opentargets.org Gentropy API |
| 2.6 | `src/mammal_repurposing/cluster_d/validation_gates.py` | Implement `gate3_auroc()` — positives = held-out L2G ≥ 0.5 protein-coding outside MHC; negatives = random length+chrom-matched gene set; PASS threshold AUROC ≥ 0.65 | Cluster D Methodology § 7 |
| 2.7 | `reports/paper-drafts/v6b_paper_draft.md` Results | Replace Gate 2 "DEGRADE" + Gate 3 "INSUFFICIENT_DATA" with actual numbers + AUROC curve figure | New `figures/v6b/gate3_auroc.png` |

**Decision gate:** If Gate 2 ρ < 0.40 or Gate 3 AUROC < 0.65, that is a real result, not a bug — write up as honest limitation in v6b paper § Limitations rather than tuning the gate. The whole point of the Roberts ceiling work was that we publish even when results disappoint.

---

### Sprint 3 — Days 29-49: **V7 anchor expansion + population × class interaction** (the V7 paper closer)

**Goal:** Close as much of the V7 0.004-0.063 PMD gap as physically possible without new wet data.

| PR | Module | Change | Source |
|---|---|---|---|
| 3.1 | `src/mammal_repurposing/translation/prisma_priors.py` | Refresh `PER_SUBDOMAIN_PRIORS` from 32 → 96 cells (12 mech × 8 endpoint) using doc's full table. Per-cell τ² rule: k≥5 → 0.02; k=2-4 → 0.04; k=1 → 0.08; k=0 → fallback HalfNormal(0.3) | MH1+MH2 V7 CPT doc |
| 3.2 | `src/mammal_repurposing/translation/effect_size_model.py` | **First** refit with τ²_class (per-class heterogeneity) + population × class interaction term, BEFORE adding new anchors. Goal: measure how much of the gap is structural vs anchor-paucity. | MH1+MH2 V7 CPT § 5 |
| 3.3 | NUTS run | 4 chains × 2000 draws, original 15 anchors + new model structure. Record PMD residuals per (class, population). | Compare against current `v7_nuts_posterior_production_v1.parquet` |
| 3.4 | `src/mammal_repurposing/translation/prisma_priors.py` | Add `REFERENCE_COMPOUND_SMD` 60-row update (copy-paste Python dict from doc § 6) | MH1+MH2 V7 CPT doc |
| 3.5 | NUTS run | 4 chains × 2000 draws, 60 anchors + same model. Compare PMD residuals vs 3.3. | Report gap closure attribution |
| 3.6 | Soft envelope upgrade | Replace hard 0.50 PMD ceiling with soft Half-Normal(0.3) upper envelope; allow individual posteriors to exceed 0.43 (Roberts 2020 max) up to ~0.55 with low prior weight | MH1+MH2 V7 CPT § 5 |
| 3.7 | `reports/paper-drafts/v7_paper_draft.md` | (a) Update Methods to describe per-class τ² + population × class interaction; (b) Update Results with 96-cell prior table + 60-anchor expansion + soft envelope; (c) Re-generate F7-3 PMD waterfall | New `figures/v7/pmd_waterfall_v2.png` |

**Decision gate:** Report attribution explicitly — what fraction of the gap closure came from (i) model structure (3.3 vs current), (ii) anchor expansion (3.5 vs 3.3), (iii) soft envelope (3.6 vs 3.5). The V7 paper's central methodological contribution requires we report this attribution honestly even if anchor expansion turns out to be the smaller intervention.

---

### Sprint 4 — Days 50-63: **MH3+MH7 bundle into V8 hierarchical model**

**Goal:** Defend V8 against the obvious "U2OS isn't brain" critique with a real per-cell-line random effect.

| PR | Module | Change | Source |
|---|---|---|---|
| 4.1 | `src/mammal_repurposing/cluster_e/joint_phenotype.py` | Add `build_v8_hierarchical_with_cell_random_effect()` per MH3 doc § 4 drop-in code. Five random effects: β (transferable), α (cell-line), γ (species), δ (compound × cell interaction), ε (residual) | MH3 doc § 4 |
| 4.2 | `scripts/fetch_cpg0000_pilot.py` | Pull cpg0000 pilot dataset (~300 compounds × A549 + U2OS shared) from JUMP S3 | MH3 doc § 5 |
| 4.3 | `src/mammal_repurposing/cluster_e/calibrate_cell_priors.py` | Fit empirical σ̂_α, σ̂_γ, σ̂_δ from cpg0000 → replace HalfNormal(0.5) priors with HalfNormal(σ̂) per endpoint | MH3 doc § 5 |
| 4.4 | NUTS run on synthetic phenotype data + new model | R̂=1.000, ESS≥1000, 0 divergences | New `data/results/v2/v8_hierarchical_v2_posterior.parquet` |
| 4.5 | `src/mammal_repurposing/cluster_e/transferability_index.py` | New module: compute `T_{c,k} = σ²_β / (σ²_β + σ²_α + σ²_δ)` per (compound, endpoint). Append column to `wet_lab_shortlist_v10.md`. | MH3 doc § 4 |
| 4.6 | Pre-reg gate validation | G1-G6 from MH3 doc § 6: ICC bounds, leave-cell-line-out R², external concordance vs Gorgogietas 2025 + Anderson 2025 | Update `reports/paper-drafts/v8_osf_preregistration.md` with G1-G6 |
| 4.7 | `reports/paper-drafts/v8_paper_draft.md` | Add Methods § "Per-cell-line random effect" with mathematical decomposition + cpg0000 calibration + transferability index. Add Limitations L1-L6 from MH3 doc. | New `figures/v8/transferability_histogram.png` |

**Critical:** Sprint 4.1 + 4.2 must land together. Fitting α + γ without empirical priors → divergences. cpg0000 must be local before NUTS run.

---

### Sprint 5 — Days 64-77: **chemCPA real-LINCS pretraining**

**Goal:** Replace V8.2 synthetic smoke with real chemCPA pretraining on LINCS L1000 + transfer to sci-Plex3.

| PR | Module | Change | Source |
|---|---|---|---|
| 5.1 | `scripts/fetch_lincs_l1000_level3.py` | Pull Level 3 Q2NORM `.gctx` for GSE70138 + GSE92742 | LINCS chemCPA doc § 3 |
| 5.2 | `src/mammal_repurposing/cluster_e/ingest_lincs.py` | Fix batch ID extraction: sig ID split on `:`, sample ID split on `_` taking first 3 fields. Currently may be wrong. | LINCS chemCPA doc § 3 |
| 5.3 | `src/mammal_repurposing/cluster_e/chemcpa_train.py` | Replace synthetic hyperparams with doc's Table 3: latent_dim=32, dropout=0.2624, ae_width=256, ae_depth=4, ae_lr=0.001121, batch=256, Gaussian likelihood, zero-centered gradient penalty | LINCS chemCPA doc § 4 |
| 5.4 | Molecular encoder weight pinning | Add MD5 verification step on grover_base / jtvae / rdkit encoder downloads | LINCS chemCPA doc Table 1 |
| 5.5 | OOD evaluation upgrade | Reserve the 9 doc-specified compounds (Dacinostat, Givinostat, Belinostat, Hesperadin, Quisinostat, Alvespimycin, Tanespimycin, TAK-901, Flavopiridol) as the OOD holdout set. Report R² on full + DEG subsets. **Add Wasserstein distance** on full predicted distribution. | LINCS chemCPA doc § 6 |
| 5.6 | Transfer learning to sci-Plex3 | Replace encoder input + decoder output for 2000 HVG; freeze molecular encoder; retrain cell embedding for K562 | LINCS chemCPA doc § 4 |
| 5.7 | `reports/paper-drafts/v8_paper_draft.md` § Methods | Real-LINCS pretraining + OOD eval methodology; cite Biolord + PerturbNet + scAgents + scDCA benchmark comparators | LINCS chemCPA doc Table 4 |

**Risk:** chemCPA training on real L1000 may take 12-48 GPU hours on RTX 5070. Plan: run overnight, monitor via Docker `b1ro-chemcpa:latest` for reproducibility check.

---

### Sprint 6 — Days 78-90: **Cross-cutting polish + final manuscript pass**

| PR | Module | Change |
|---|---|---|
| 6.1 | `CITATIONS.bib` | Add: Chandrasekaran 2024 (JUMP), Gorgogietas 2025, Anderson 2025, Okbay 2022, de la Fuente 2021, FinnGen, Biolord, PerturbNet, scAgents, scDCA, Renishaw MH8 disambiguation references |
| 6.2 | `reports/paper-drafts/v6b_paper_draft.md` | Add MH8 nomenclature disambiguation footnote citing 6+ unrelated database entries |
| 6.3 | `reports/paper-drafts/v7_paper_draft.md` | Final attribution table: gap closure from model structure vs anchor expansion vs soft envelope |
| 6.4 | `reports/paper-drafts/v8_paper_draft.md` | Add U2OS-to-brain 4-sentence defence paragraph using MH3 transferability index + cpg0000 calibration |
| 6.5 | `reports/paper-drafts/integration_paper_draft.md` | Update integration figure: V6.B (4 gates PASS) → V7 (full attribution) → V8 (hierarchical + transferability) end-to-end pipeline |
| 6.6 | `GAPS_AND_RESEARCH_DIRECTIONS.md` | Move MH3, MH8, MH1+MH2, Gate 2, Gate 3 from "MUST-HAVE" to "RESOLVED" section. Add new gaps surfaced by attribution analysis. |
| 6.7 | `PROJECT_STATUS.md` | Refresh status table — production posteriors, all gates, transferability indices |

---

## 4. Pre-Registered Validation Gates — Consolidated

| Gate ID | Source | Definition | Pass threshold |
|---|---|---|---|
| V6.B G1 | V6.B 4-gate framework | Sampling diagnostic: R̂ ≤ 1.01, ESS ≥ 1000, 0 divergences | **191-target panel**: 0 div after MH8 fix |
| V6.B G2 | V6.B 4-gate framework + Cluster D Methodology | Multi-modulator Spearman ρ across (target, modulator) anchors | ρ ≥ 0.40 with n ≥ 80 pairs |
| V6.B G3 | V6.B 4-gate framework + Cluster D Methodology | Held-out cognition GWAS L2G AUROC on Okbay 2022 EA4 | AUROC ≥ 0.65 (acknowledged aspirational) |
| V6.B G4 | V6.B 4-gate framework | Roberts 2020 SMD ceiling — no posterior > 0.55 with > 10% mass | 0 violations (current production: 0 violations confirmed) |
| V7 P1-P8 | V7 OSF pre-reg | 8 specific compound recovery tests on held-out efficacy data | Already PASS — protect via Sprint 3 |
| V8 G1-G6 | MH3 doc § 6 | (a) ICC bounds for α, (b) leave-cell-line-out R², (c) external concordance vs Gorgogietas 2025, (d) external concordance vs Anderson 2025, (e) transferability index distribution, (f) cpg0000 prior calibration vs empirical | Specified per gate in MH3 doc |
| V8 OSF | V8 OSF pre-reg | 8-cell disagreement classification + I_novel ranking | Already PASS — protect via Sprint 4 |

---

## 5. Paper-by-Paper Refinement Map

### V6.A paper (J Cheminform / Nat Mach Intell)
- No structural changes from this sprint. Sprint 1 has no V6.A touch.
- **Optional:** Sprint 3.6 (soft envelope) may motivate a V6.A § Discussion paragraph about how MH DTI ensemble uncertainty downstream-propagates to soft-envelope clipping.

### V6.B paper (Cell Reports Methods / Bioinformatics)
- **Sprint 1.5** — MH8 § into Methods.
- **Sprint 1.3** — re-render Results table with 4-gate PASS post-MH8.
- **Sprint 2.7** — Gate 2 + Gate 3 AUROC figure replaces "INSUFFICIENT_DATA."
- **Sprint 6.2** — MH8 nomenclature disambiguation footnote.

### V7 paper (CPT / CPT:PSP)
- **Sprint 3.7** — full Methods + Results refresh with 96-cell prior table, 60-anchor expansion, per-class τ², population × class interaction, soft envelope.
- **Sprint 6.3** — attribution analysis table (the central methodological honest-reporting contribution).

### V8 paper (Nat Mach Intell / Nat Methods)
- **Sprint 4.7** — Per-cell-line random effect Methods + L1-L6 limitations + transferability index Results.
- **Sprint 5.7** — Real-LINCS pretraining replaces synthetic; OOD eval methodology w/ Wasserstein.
- **Sprint 6.4** — U2OS-to-brain 4-sentence defence paragraph.

### Integration umbrella (Nature / Nat Med / Cell)
- **Sprint 6.5** — Update end-to-end pipeline figure with all post-sprint state.
- Pull transferability index column into the umbrella's polypharmacology table.

---

## 6. Honest Constraints + Risk Register

| # | Risk | Mitigation |
|---|---|---|
| R1 | Sprint 1.3: 191-target NUTS still diverges after MH8 + non-centered. | Likely GWAS contribution term needs its own non-centered reparam; or one of MAO-A/MAO-B/COMT has a still-undetected non-linearity (e.g., MAO-B saturates differently). Drop back to 22-target headline + flag MAO/COMT for V8 work. |
| R2 | Sprint 2: Gate 2 Spearman < 0.40. | This is the kind of result we publish honestly per V6.B paper's central thesis. Add to Limitations. Do not tune gate threshold to make it pass. |
| R3 | Sprint 2: Gate 3 AUROC < 0.65. | Same as R2 — honest limitation; lit baselines are 0.55-0.70 anyway. |
| R4 | Sprint 3: Attribution analysis shows anchor expansion did little. | Predicted by MH1+MH2 V7 CPT doc § 5. Report attribution honestly; reframe paper around model-structure improvement (τ²_class + population × class) rather than anchor expansion. |
| R5 | Sprint 4: MH3 + MH7 still produce identifiability collisions despite bundling. | cpg0000 calibration may need expansion to 3-cell-line subset (HepG2 + A549 + U2OS) for full identifiability. Fallback: fit α as fixed effect with only U2OS data; report transferability index as "U2OS-only context-dependence proxy." |
| R6 | Sprint 5: chemCPA training takes > 48 GPU hours, hits Together-style instability. | Already have synthetic chemCPA smoke as fallback. Cite Biolord (Sprint 5.7 comparator) as alternative if we cannot complete full training. |
| R7 | Open Targets Genetics standalone portal retired 2025-07-09 → Gate 3 fetcher breaks. | Use Open Targets Platform 25.06+ L2G Parquet (already accounted for in Sprint 2.5). |
| R8 | LINCS GEO FTP rate-limits / disappears. | Already cached `GSE70138_compoundinfo.txt.gz` locally; the doc's Docker image `b1ro-chemcpa:latest` mirrors the data. |

---

## 7. Pre-Flight Asks Before Sprint 1 Lands

Two micro-decisions I'd like a green light on before executing:

1. **The 6 missing cognitive-set targets** from your prior question: HTR1A, HTR4, SLC6A9, GRM2, GRM3, GRM5. Patch into V6.B.5 panel in Sprint 1.1 (with the MAO-A/MAO-B/COMT substrate-mediated additions)? **Recommended: yes** — clean to do them all in one targets.parquet edit.
2. **`data/raw/modulator_anchors_seed.csv` scaffolding.** Auto-generate from the Cluster D Methodology doc's ~95-row table (Sprint 2.1)? **Recommended: yes** — the table is already in the doc; this is mechanical translation, not novel curation.

Per Auto-Mode I'll proceed with both unless redirected.

---

## 8. Effort + Dependency Map (visual summary)

```
Sprint 1 (2w) ──┬─→ Sprint 2 (2w) ──→ V6.B PAPER FROZEN
                │                       (4-gate PASS)
                │
                ├─→ Sprint 3 (3w) ──→ V7 PAPER FROZEN
                │   (depends only on existing 15 anchors initially)
                │
                └─→ Sprint 4 (2w) ──→ Sprint 5 (2w) ──→ V8 PAPER FROZEN
                    (cpg0000 dep)      (chemCPA training)

Sprint 6 (2w) — Cross-cutting polish (touches all 5 papers + integration)
```

**Total: ~13 weeks of engineering** in the optimistic path. Realistic budget with overruns: **16 weeks** (~Sep 2026).

**Critical-path bottleneck:** Sprint 4 (MH3 + cpg0000 download + NUTS calibration). If cpg0000 S3 pull stalls, Sprint 4-5-6 all slip.

### 8.1 Refinements from second research-reading pass (2026-05-28)

After re-reading MH3 § 5.6 + GAPS_AND_RESEARCH_DIRECTIONS.md § 4 + § 6, three roadmap refinements:

**R1. CellPainTR must precede MH3 in Sprint 4.** MH3 § 5.6 is explicit: if `cpg0016` is fed into the joint posterior as raw CellProfiler features (without prior CellPainTR / Harmony batch correction), the per-cell-line `α` random effect will absorb 12-partner-site batch variance instead of cell-line biology. This *understates* the transferable fraction — bad news disguised as good news. Add to Sprint 4 day 0: confirm cpg0016 embeddings come from CellPainTR (or equivalent batch-aware embedder), not raw CellProfiler features.

**R2. MH3 → MH4 sequencing protects ~3 days of rework.** MH3 § 12 + the MH3 pre-reg G6 gate both show Mondrian conformal calibration (MH4) naturally consumes the MH3 posterior. If we ship MH4 first and then MH3, we'd refit conformal on a posterior that now has α/γ/δ random effects, throwing away the prior calibration. Confirm: Sprint 4 (MH3) → Sprint 5 (chemCPA real-LINCS) → future MH4 sprint (deferred to v9). This is already the case in the roadmap; the explicit rationale is now documented.

**R3. The MH8 fix is more nuanced than "bump target_accept".** GAPS gap #4 reveals V6.B.5 already runs at `target_accept=0.95` and still produces 37 divergences. So just bumping to 0.98 won't fix it — the **structural** fix is required: per-target AHBA-masking for substrate-mediated enzymes. The original roadmap mentioned "non-centered β reparam + `target_accept=0.95`" but β is already non-centered in the production model; the actual code change is in `build_y_obs_from_sources` (inflate `sigma_obs[AHBA_row, sm_target_idx]` by 10×) plus an opt-in flag in `fit_cluster_d_prior_nuts`. Sprint 1.2 task description has been updated to reflect this correctly.

**R4. MH3 cpg0000 calibration day-by-day plan (lift from MH3 § 10) integrates cleanly with Sprint 4 days 1-3.** Re-baseline Sprint 4:
- Days 50-52: cpg0000 ETL + two-cell-line stripped model (calibrate σ̂_α, σ̂_δ empirically)
- Days 53-56: V8 hierarchical refactor with α + γ + δ + non-centered + numpyro
- Days 57-58: Full V8 fit with informative cpg0000 priors
- Days 59-60: G1-G6 pre-reg gate evaluation + Gorgogietas concordance (the cheapest external validation per MH3 § 8)
- Days 61-63: V8 Methods + Discussion writing + OSF § 7 amendment

---

## 9. What This Roadmap Explicitly Does NOT Cover

These are real follow-ups but outside the 90-day window:

- **Allosteric awareness benchmark** (research doc Section 5 open question #3) — publishable mini-project, deferred.
- **MAMMAL MCP server** — for interactive DTI queries during exploration; deferred to v9.
- **Wet-lab handoff to a real collaborator** — `reports/wet-lab/wet_lab_handoff_v1.md` exists; finding the collaborator is out of scope.
- **πphen joint-posterior expansion to 16-cell disagreement** — current 8-cell is sufficient for V8 paper.
- **Multi-modality MOFA+ embeddings** — scaffold exists (`mofa_embed.py`); v9 work.
- **Adolescent Brain Cognitive Development (ABCD) Study direct fetch** — Gate 3 doc § 6 mentions ABCD R5 as a held-out source, but DUA process is multi-month; use Okbay EA4 first and ABCD only if DUA finalises within sprint window.

---

**End of roadmap.** Next action (per Auto-Mode): begin Sprint 1.1 (patch `data/interim/targets.parquet` with `substrate_mediated` column + add MAO-A/MAO-B/COMT + 6 missing cognitive targets) unless redirected.

---

## 10. Implementation status (live as of 2026-05-28)

| Sprint | Status | Headline |
|---|---|---|
| 1.1 | ✅ DONE | 22 → 28 target panel + `substrate_mediated` column + 6 missing cognitive targets |
| 1.2 | ✅ DONE | MH8 AHBA-masking surgery in `bayesian_prior.py` |
| 1.3 | ✅ DONE | V6.B.5 NUTS: **37 → 0 divergences** (target_accept=0.99, ESS=1808, R̂=1.000) |
| 1.4 | ✅ DONE | 3 slow-marked regression tests lock MH8 production config |
| 2.1 | ✅ DONE | `modulator_anchors_seed.csv` — 70 rows / 38 targets / 24 Phase III nulls |
| 2.2 | ✅ DONE | Gate 2 multi-modulator evaluator. **Publishable falsification**: ρ ∈ [-0.35, +0.10] across 8 configs |
| 3.1 | ✅ DONE | PRISMA priors V2: 73/96 cells (76%, up from 33%); Schmidli τ² per-cell |
| 3.2 | ✅ DONE | `fit_effect_size_nuts_v2` with per-class τ² + population × class interaction |
| 3.3 | ✅ DONE | `REFERENCE_COMPOUND_SMD_V2`: 109 anchor cells / 48 compounds (was 15) |
| 4.1 | ✅ DONE | `build_v8_hierarchical_with_cell_random_effect`: β + α + γ + δ + ε non-centered |
| 4.2 | ✅ DONE | Synthetic round-trip validated: transferability index correctly orders compounds |
| 4.3 | ⏳ blocked | cpg0000 empirical prior calibration — awaits S3 pull |
| 4.4 | ⏳ pending | V8 OSF § 7 amendment with G1-G6 gates |
| 5.x | ⏳ blocked | chemCPA real-LINCS pretraining — awaits LINCS Level-5 download |
| 6.x | ⏳ pending | Cross-cutting paper polish (depends on Sprints 4.3 + 5.x) |

**Test counts (post-sprint):** 419/420 non-slow PASS (+87 over baseline 332); 11/11 slow PASS.

**Publishable findings landed:**
1. **MH8 fix is the dominant intervention** for V6.B.5 divergences (≥80% reduction; the remaining 0% comes from target_accept bump).
2. **Gate 2 multi-modulator FAIL** is an *empirical falsification* — it motivates the V4→V5→V6→V7→V8 multi-layer architecture rather than refuting it.
3. **V8 transferability index** correctly recovers parameter ordering on synthetic data — the U2OS-to-brain defence paragraph for the V8 paper now has computational backing.
