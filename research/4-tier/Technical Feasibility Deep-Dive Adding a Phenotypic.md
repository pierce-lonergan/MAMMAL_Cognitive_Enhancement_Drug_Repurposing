# Technical Feasibility Deep-Dive: Adding a Phenotypic/Perturbational Evidence Axis (πphen) as a Third Independent Bayesian Factor in a Target-First Cognition Repurposing Pipeline

## TL;DR
- **Build it, in stages.** Integration of πphen as a third Bayesian factor is technically feasible with off-the-shelf, well-maintained tooling: LINCS L1000 (≈1.3M signatures, GCTX via cmapPy / signatureSearch / SigCom-LINCS), JUMP-CP Cell Painting (the full ≈116k-compound production set is Chandrasekaran SN et al., bioRxiv 2023.03.23.534023; CPJUMP1 benchmark is Chandrasekaran 2024 Nat Methods 21:1114-1121; 656 TB image+numerical as of January 2024 per Weisbart 2024 Nat Methods 21:1775-1777), EPA ToxCast MEA (1,055 compounds, Strickland 2018; 384-compound concentration-response × 43 parameters, Kosnik 2020 Arch Toxicol 94:469-484), and chemCPA / biolord (Piran 2024 Nat Biotechnol 42:1678-1683) for in-silico imputation. Storage/compute targets (≈1–3 TB working set, single RTX 5070) are realistic for a research-grade build.
- **The math is naive-Bayes log-odds with one important asterisk.** Under conditional independence *given the cognition-modulator hypothesis H*, π_joint ∝ π_target · π_genetic · π_phen and equivalently logit(π_joint) = Σ wᵢ · logit(πᵢ) + b. Calibrate each modality with isotonic regression on the ≈40–80 known pro-cognitive anchors via leave-one-class-out, then wrap the joint posterior in Mondrian conformal prediction (crepes) for finite-sample coverage. Independence breaks most plausibly for L1000 vs OT-Genetics (both transcript-mediated) — measure the residual correlation on the anchor set and model it with a Gaussian copula or a learned interaction weight.
- **The headline value is the disagreement structure.** Three-way Jensen-Shannon divergence JS₃(π_t, π_g, π_p) ∈ [0, log 3] and the novel-mechanism score I_novel(C) = π_p · [1 − I(π_p ; (π_t, π_g))] together partition the candidate space into 8 disagreement cells; the (low-target, low-genetic, high-phen) cell is the one the existing V4–V6 stack cannot see and is the principal scientific justification for building this axis — it is where compounds like clemastine (Mei 2014 Nature Medicine 20:954) and PIPE-307 (Poon 2024 PNAS) live.

## Key Findings

1. **Data is genuinely ready.** L1000 Phase I (GSE92742) + Phase II (GSE70138) are downloadable today; JUMP-CP CellProfiler + DeepProfiler parquet profiles are on AWS Open Data at `s3://cellpainting-gallery/cpg0016-jump/` with `--no-sign-request` access; ToxCast MEA acute and concentration-response data are in `invitroDB` (EPA Clowder); chemCPA / biolord have public PyTorch / scvi-tools code.

2. **chemCPA's "uncertainty" is a heuristic, not a posterior.** Per Lotfollahi 2023 Mol Syst Biol e11517 (PMC10258562), the theislab/cpa uncertainty estimate is a latent-space distance: *"very different covariate/perturbation vector combinations to those observed in training data will have a higher distance and thus higher uncertainty than covariate/perturbation observed in training"* — the model is a deterministic autoencoder; the decoder σ² is aleatoric, not epistemic. Biolord (Piran et al. 2024 Nat Biotechnol 42:1678-1683) outperforms chemCPA on held-out sci-Plex3 drugs (mean r² = 0.76 ± 0.0005 vs chemCPA-pre = 0.51 ± 0.0062) but also offers no calibrated predictive variance. GPerturb (Xing H & Yau C, Nat Commun 16:5423, 1 July 2025, doi:10.1038/s41467-025-61165-7) provides Bayesian UQ but only for genetic perturbations. **You must wrap chemCPA in a 5–10 model deep ensemble or MC-dropout to get usable epistemic uncertainty for Bayesian integration.**

3. **Cell-line context is the dominant πphen risk.** JUMP-CP is U2OS / A549; ToxCast MEA is rat primary cortical; L1000 has NEU/NPC/ASC lines but they cover only a small subset of compounds. Cognition-relevant cell-line filtering during connectivity scoring (NPC, NEU, ASC, NEU.KCL) is essential — naive averaging across all cells dilutes neuro-specific signal.

4. **Remyelination is a first-class cognition pathway with usable signatures.** Mei F et al. 2014 Nature Medicine 20:954-960 (the BIMA-8 cluster: clemastine, benztropine, atropine, ipratropium, oxybutynin, trospium, tiotropium, quetiapine), Najm FJ et al. 2015 Nature 522:216-220 (miconazole + clobetasol with RNA-seq + phosphoproteomics in supplement), Lariosa-Willingham KD et al. 2016 BMC Res Notes 9:419 (27 NCC hits), and Poon MM et al. 2024 PNAS (PIPE-307 M1R-selective antagonist) provide a usable positive-anchor set for an explicit "remyelination sub-signature" inside the cognition reference.

5. **Pro-cognitives are in JUMP-CP but NOT in ToxCast MEA.** Donepezil, memantine, modafinil, galantamine, rivastigmine are in the Broad Drug Repurposing Hub (profiled in JUMP via JUMP-Target/jump-cellpainting GitHub). They are absent from the EPA ToxCast Phase II MEA library, which is dominated by pesticides (Strickland 2018 PMC6438628 names abamectin, lindane, prallethrin, haloperidol, reserpine as representative hits). Plan for chemCPA-style MEA imputation rather than direct MEA queries for nootropics.

6. **The framework's principal scientific value is captured by one quantity:** I(π_p ; (π_t, π_g)) — the mutual information between the phenotypic factor and the joint of the two target-first factors. High π_p and low I = phenotypic activity not predictable from target+genetic = candidate novel-mechanism compound.

## Details

### 1. Modality 1 — LINCS L1000 + CMap Connectivity

**Data inventory and access.** Phase I (GSE92742, ~1.3M Level 5 signatures, 2017) and Phase II (GSE70138, ongoing) ship as GCTX (HDF5) at GEO. Programmatic surface area:
- **cmapPy** (Python, Broad) — `cmapPy.pandasGEXpress.parse.parse(file_path, rid=None, cid=None, ridx=None, cidx=None, row_meta_only=False, col_meta_only=False, ...)` returns a `GCToo` (3 pandas frames: data_df, row_metadata_df, col_metadata_df). Read level-5 metadata-only first, filter sig_ids by `cell_iname`, `pert_iname`, `pert_idose`, `pert_itime`, then re-parse with `cid=`. Per Enache et al. 2020 Bioinformatics on the GCTx format: GCTx is HDF5-backed, parses orders of magnitude faster than GCT text on large matrices.
- **signatureSearch** (R/Bioconductor; Duan/Girke) — `qSig` + `gess_lincs` implements the Subramanian 2017 bi-directional weighted KS and tau scoring. Implements CMAP, LINCS, gCMAP, Fisher, and Cor algorithms; signatureSearchData hosts pre-built GES-DBs on Bioconductor's ExperimentHub.
- **SigCom LINCS** (Ma'ayan lab) — REST microservices over >1.5M signatures from LINCS + GEO + GTEx (Evangelista et al. 2022 Nucleic Acids Res 50:W697-W709); useful for signature commons queries with SmartAPI documentation.
- **iLINCS** + **drugfindR** R package (Imami AS et al. 2024, Zenodo 5.281/zenodo.338354715) — Consensus Gene Knockdown, Overexpression, and Chemical Perturbagen signatures.

**Storage.** Phase I+II GCTX Level 5 ≈ 50 GB; Level 4 instances ≈ 200 GB if you want replicate-level data. Use chunked HDF5 reads — never load the full matrix at once.

**Connectivity score math (Subramanian A et al. 2017 Cell 171:1437, exactly).** For query signature q = (q_up, q_down) and reference signature r:

WTCS(q, r) = (ES_up − ES_down) / 2  if sign(ES_up) ≠ sign(ES_down), else 0

where ES is the weighted Kolmogorov–Smirnov enrichment statistic (Subramanian 2005 GSEA). Per clue.io/connectopedia: *"WTCS is a composite, bi-directional version of ES… ranges between −1 and 1."* Then:

NCS(q, r, c) = WTCS / μ⁺_c  if WTCS > 0;  WTCS / |μ⁻_c|  if WTCS < 0

(c = cell line; μ⁺/μ⁻ are signed means of WTCS within cell line c). Tau is the quantile of NCS against a fixed reference distribution: τ ∈ [−100, 100], where τ = 90 means only 10% of reference perturbagens have stronger connectivity. Cross-cell-summary uses the **maximum-quantile (max-q) statistic**: per the LINCS documentation, *"NCSct: NCS summarized across cell types… compare 67 and 33 quantiles of NCSp,c and retain whichever is of higher absolute magnitude."*

**CNS-relevant cell lines in L1000.** Core LINCS panel includes NPC (neural progenitor) and small numbers of NEU and NEU.KCL signatures; most signatures are in MCF7, PC3, A375, HA1E, HCC515, HEPG2, A549, VCAP, HT29. For cognition queries, restrict the reference to {NPC, NEU, NEU.KCL, ASC} when available and additionally weight non-CNS lines by an empirical Bayes factor learned from the agreement of compounds dual-profiled across both.

**Polypharmacology.** L1000 signatures are dose- and time-resolved (typically 24h, 6h; 10 μM / 100 nM). Use the highest-TAS (Transcriptional Activity Score) signature per (pert_id, cell_iname) pair, falling back to a transcriptional-strength-weighted average when multiple doses are present. The `pert_idose` canonicalization in clue.io build settings provides the dose discretization needed.

**Recommended pipeline.** `cmapPy → filter to cognition-anchor query → signatureSearch::gess_lincs (LINCS algorithm, BING gene set) → per-compound τ across CNS-relevant cells → NCSct cross-cell summary → calibrate to P(phenotypic match | C is cognition-modulator) via isotonic regression on pro-cognitive anchors.`

---

### 2. Modality 2 — JUMP-CP Cell Painting (2024 release)

**Data scale and structure.** The JUMP Consortium produced ≈116,750 compound perturbations + 7,975 CRISPR KOs + 12,602 ORFs (full production dataset: Chandrasekaran SN et al., bioRxiv 2023.03.23.534023, still a preprint as of the 2024 Nat Methods reference list). The CPJUMP1 benchmark sub-paper is Chandrasekaran SN et al. 2024 Nat Methods 21:1114-1121 (hundreds of perturbations, 75M cells, U2OS + A549); the genetic subset was separately published as Chandrasekaran 2025 Nat Methods 22:1742-1752. Total Cell Painting Gallery as of January 2024: *"656 terabytes (TB) of image and associated numerical data"* (Weisbart E et al. 2024 Nat Methods 21:1775-1777).

**File layout (canonical, per cellpainting-gallery docs):** each source folder under `s3://cellpainting-gallery/cpg0016-jump/<source>/` has:
- `images/<batch>/illum/<plate>/<plate>_Illum{AGP,DNA,ER,Mito,RNA,Brightfield}.npy` + raw multi-channel TIFFs in `images/<batch>/images/<plate>__...-Measurement.../`.
- `workspace/load_data_csv/<batch>/<plate>/load_data.csv` — CellProfiler LoadData input.
- `workspace/backend/<batch>/<plate>/<plate>.sqlite` — single-cell features.
- `workspace/profiles/<batch>/<plate>/<plate>.parquet` — well-level CellProfiler features (~3,000 features per well).
- `workspace_dl/profiles/efficientnet_v2_imagenet1k_s_feature_vector_2_<hash>/<batch>/<plate>/<plate>.parquet` — DeepProfiler embeddings.

**Feature space and normalization.** CellProfiler returns ~3,000 features per well: shape, intensity, texture, granularity, neighbor, and RadialDistribution across Cell/Nuclei/Cytoplasm × {AGP, DNA, ER, Mito, RNA} channels. Canonical processing (pycytominer 1.3.0):
1. Annotate metadata (`pycytominer.annotate`).
2. **RobustMAD** per-plate normalization against DMSO negatives (`scaled = (x − median) / mad`; `pycytominer.operations.transform.RobustMAD(epsilon=1e-18)`).
3. Feature selection — drop invariant, blocklist, highly-correlated, low-variance, outlier features.
4. **Sphering / whitening** (`pycytominer.normalize` with `method="spherize", spherize_method="ZCA-cor", spherize_epsilon=1e-6`) against negative controls.
5. **Harmony** for batch correction across the 12 contributing sites (the jump-profiling-recipe applies this through `orf.json`/`crispr.json`).

After processing, a well-level profile is ~600-1,000 selected features.

**Compute footprint.** Storing the full processed compound parquet set is ~30 GB; pulling only the per-plate compound profiles you need is feasible with:
```
aws s3 sync --no-sign-request s3://cellpainting-gallery/cpg0016-jump/<source>/workspace/profiles/<batch>/<plate>/ .
```
Image storage (multi-TB) is only required if you want to re-extract DeepProfiler / vision-transformer embeddings.

**Reference profile construction for cognition.** Take per-compound well profiles for the pro-cognitive anchors and the disease-state-reversal anchors, normalize to DMSO, then compute the centroid and per-feature covariance. Use Mahalanobis distance to the cognition centroid (or, better, cosine similarity to multiple sub-centroids — see §"Reference Signature Construction"). The CPJUMP1 framework specifically benchmarks "matching perturbations" via mean Average Precision (mAP) on a sphered profile; reuse `copairs` (Broad) / Arevalo et al. evaluator scripts for benchmark-grade mAP.

**Compound matching across L1000 and JUMP-CP.** Map by InChIKey: JUMP `metadata/compound.csv.gz` ↔ L1000 `pert_inchikey`. Per Way GP et al. 2022 Cell Systems 13(11):911-923: *"L1000 captures activity of compounds targeting MAPK family genes and heat shock protein (HSP90AA1) strongly, whereas Cell Painting captures aurora kinase genes (AURKA, AURKB), PLK genes (PLK1, PLK2, and PLK3), and BRD4 with high precision (Figure 4E)."* This complementarity is the *empirical foundation* for treating them as quasi-independent likelihood sources.

**Cell-context caveat.** JUMP is in U2OS — not neurons. For cognition modulation, U2OS morphology captures off-target / cytotoxicity / cell-cycle / cytoskeletal effects but not neuronal differentiation or remyelination. The reference signature must be defined *as the cognition-anchor morphological centroid in U2OS feature space*, not as an a priori biological signature. This is acceptable because Cell Painting's value here is MoA-cluster grouping, not biologically-interpretable feature direction.

---

### 3. Modality 3 — iPSC-Derived Cortical Neuron MEA

**Hardware ecosystem.** Axion BioSystems Maestro Pro/Edge (12/48/96-well; 768-electrode total) dominates published screens; Multi Channel Systems MEA2100 and 3Brain BioCAM are alternatives. AxIS Navigator's Neural Metric Tool exports the canonical metric set: weighted mean firing rate (wMFR), single-electrode burst rate, network burst rate, ISI coefficient of variation, synchrony index (area under normalized cross-correlogram, AUNCC), full-width at half-height (FWHH) of cross-correlogram, inter-burst interval, spike amplitude, network burst duration. Default Adaptive burst detection (50 spikes / 100 ms ISI / 35% electrodes / 20 ms synchrony window) is the de-facto standard (Axion application notes).

**Cell-source matrix (recommended for cognition).**
- FUJIFILM Cellular Dynamics iCell GlutaNeurons + iCell Astrocytes (canonical commercial cortical glutamatergic + astrocyte co-culture; matures to synchronized network bursting by ~DIV14).
- BrainXell / Cellartis (Takara) Cortical Neurons.
- NGN2-induction iPSC lines (Zhang/Südhof) for fast (~3 week) homogeneous glutamatergic neurons.
- Dual-SMAD inhibition (Chambers 2009) cortical neurons for more heterogeneous but more developmentally faithful cortical identity.
- **Co-culture with primary astrocytes is required for stable synchronized bursting** per multiple Axion application notes (e.g., the Quick-Neuron + primary astrocyte Elixirgen protocol) and Bioarxiv 464677 (Eichler et al., dual-SMAD on MEA, ACM + 2% O₂).

**Public compound-screen datasets.**
- **EPA ToxCast MEA acute** (Strickland JD et al. 2018, PMC6438628): 1,055 ToxCast Phase II compounds at 40 μM in rat primary cortical, 326 hit threshold-crossing (308 down, 18 up); endpoints in `invitroDB` under `CCTE_Shafer_MEA_acute_*` (e.g., `CCTE_Shafer_MEA_acute_spike_number`, EDA endpoint ID 2454).
- **EPA ToxCast MEA concentration-response** (Kosnik MB et al. 2020 Arch Toxicol 94:469-484, PMC7371233): *"we expand this approach by retesting 384 of those compounds (including 222 active in the previous screen) in concentration-response across 43 network activity parameters."* 237 response-active; hierarchical clustering + ML identifies 15 key parameters as MoA discriminators.
- **EPA ToxCast MEA pilot** (Valdivia P et al. 2014 Neurotoxicology, PubMed 24997244): 92 compounds; MEAs detected 30/37 NVS ion-channel hits and 35/48 known neuroactives.
- **EPA ToxCast MEA developmental** (Brown JP et al. 2016, Frank CL et al. 2017): network ontogeny screens.
- **NeuroLINCS** (Maor-Nof et al. 2022 Sci Data 9: 638) — public iPSC multi-omics (DIA-MS proteomics, transcriptomics, epigenomics) on 6 hiPSC lines × 3 conditions (control, C9orf72 ALS, SMA); useful for cell-type calibration but motor neurons (not cortical) and no compound screen.

**Coverage gap.** None of the EPA/ToxCast MEA libraries contain canonical pro-cognitives (modafinil, donepezil, memantine, racetams, ampakines, α7 PAMs). For cognition queries, the MEA axis must be largely *imputed* by chemCPA-style surrogate prediction trained on the ToxCast set, then evaluated on the few anchors with primary literature MEA data (modafinil and donepezil appear in Axion application notes and McConnell HL 2012 Neurotoxicology).

**Recommended pipeline.** Quantify reference "pro-cognition MEA signature" as a 15-dim vector (Kosnik 2020 reduced set) of relative changes vs. DMSO. Compute Mahalanobis distance from a candidate's measured (or imputed) vector to the reference centroid. Calibrate to probability with isotonic regression on the small set of compounds with both MEA data and a known cognition classification.

---

### 4. Modality 4 — chemCPA / generative signature prediction

**Architecture.** CPA (Lotfollahi M et al. 2023 Mol Syst Biol 19:e11517) factorizes single-cell gene expression as x = (basal z_b) + Σ_p (perturbation z_p · dose) + Σ_c (covariate z_c), trained with adversarial gradient reversal on the basal latent to remove confounding. The decoder outputs Gaussian (μ, σ²) per gene (MultiCPA uses Negative Binomial: per Inecik K et al. 2022 ICML WCB workshop, the NB decoder is *"specified by the mean µG and the dispersion"*). chemCPA (Hetzel L et al. 2022, NeurIPS / arXiv 2204.13545) extends the perturbation embedding with a learned chemical encoder over RDKit / Grover / JT-VAE features so unseen drugs can be embedded.

**Training data.** sci-Plex3 (Srivatsan 2020; 188 compounds × 3 cell lines × ~290k cells), Trapnell drug perturbation atlas, plus pretraining on bulk L1000 (~1.3M signatures). Recipe: pretrain on L1000, fine-tune on sci-Plex (architecture surgery for transfer learning, per Hetzel 2022).

**Held-out performance.** On 9 OOD sci-Plex3 compounds (Dacinostat, Givinostat, Belinostat, Hesperadin, Quisinostat, Alvespimycin, Tanespimycin, TAK-901, Flavopiridol), chemCPA achieves median r² on DEGs of 0.38 vs 0.85 in-distribution (Lotfollahi 2023 Fig. 4). Biolord (Piran Z, Cohen N, Hoshen Y, Nitzan M 2024 Nat Biotechnol 42:1678-1683) achieves on the same held-out set *"mean r²; chemCPA-pre: 0.51 ± 0.0062, biolord: 0.76 ± 0.0005"*, and remains robust under subsampling (*"mean r²: 0.63 ± 0.0003) over 10% of the data."* PRnet (Qi N et al. 2024 Nat Commun 15, PMC11513139) trains on 883,269 L1000 profiles across 82 cell lines and 175,549 compounds, but reports only point estimates: *"PRnet initially predicts the average transcriptional profile, fold-change in the gene, and gene rank."*

**Uncertainty quantification — critical caveat.** The theislab/cpa "uncertainty" advertised in the README is a heuristic latent-space distance, not a calibrated predictive variance: *"the current heuristic for uncertainty estimation originates from the compositional formulation in CPA. Such formulation results in very different covariate/perturbation vector combinations to those observed in training data will have a higher distance and thus higher uncertainty than covariate/perturbation observed in training (see Appendix Figs S4 and S12)"* (Lotfollahi 2023 MSB). CPA notes explicitly: *"we opted to implement a deterministic autoencoder scheme, extensions toward variational models … are straightforward."* Biolord and PRnet provide no calibrated UQ either. GPerturb (Xing H & Yau C 2025 Nat Commun 16:5423, doi:10.1038/s41467-025-61165-7) provides explicit Bayesian per-gene posterior of effect existence and magnitude but is restricted to genetic (CRISPR/Perturb-seq) perturbations.

**Practical recommendation for Bayesian integration.** Use a 5–10 model deep ensemble of chemCPA (different seeds + 80% bootstrap of training data). Take ensemble mean as predicted L1000/CP signature, ensemble variance for epistemic σ²_epi, decoder σ² for aleatoric σ²_alea. Calibrate the resulting per-compound predictive distribution via Mondrian conformal prediction conditioned on compound class (kinase / GPCR / receptor / etc.) using the `crepes` library. RTX 5070 (12 GB) handles training of CPA on L1000 in ~4–8 h per seed; ensemble training in 2–3 days.

---

### 5. Reference Signature Construction (Cognition Phenotype Definition)

**Single centroid vs. mixture-of-experts.** A single cognition centroid blurs distinct MoAs (a cholinesterase inhibitor and an ampakine are both "pro-cognitive" but transcriptionally orthogonal). Recommended: **a Bayesian mixture of K=5 sub-centroids**, each a known pro-cognitive MoA cluster, plus disease-reversal anchors as inverted sub-centroids:

| Sub-signature k | Anchor compounds | Pathway / direction |
|---|---|---|
| 1: Cholinergic enhancement | donepezil, galantamine, rivastigmine, huperzine A; α7 nAChR PAMs (GAT-107, PNU-120596) | AChE inhibition + α7 nAChR PAM |
| 2: Catecholaminergic / wake-promoting | modafinil, methylphenidate, atomoxetine, guanfacine | DAT/NET, α2A |
| 3: Glutamatergic — NMDA modulation + AMPA potentiation | memantine; ampakines (CX-516, CX-717); racetams (piracetam, aniracetam, oxiracetam, pramiracetam) | NMDA, AMPA, mGluR |
| 4: Trophic / neurogenic / ISR | vortioxetine, NSI-189, tianeptine; BDNF mimetics (7,8-DHF); ISRIB | TrkB; eIF2B (per Anand & Walter 2020 FEBS J 287:239 — ISRIB *"enhances eIF2B exchange activity"* by stapling two βγδε tetramers); hippocampal neurogenesis |
| 5: Remyelination / OPC differentiation | clemastine, miconazole, clobetasol, benztropine, ipratropium, oxybutynin, trospium, tiotropium, quetiapine, atropine (the Mei 2014 BIMA-8 cluster); PIPE-307 (Poon 2024 PNAS); doxepin, orphenadrine (the M1R follow-on cluster identified in the bioRxiv 2023.07.11.548469 combination screen) | M1R antagonism → ERK1/2 → Myrf/Olig2; corticosteroid receptor pathway |

**Disease-state inverse anchors:** ROSMAP / MSBB / Mayo LOAD AD cortex DEGs, CommonMind schizophrenia DLPFC DEGs, Allen Aging/Dementia/TBI MCI hippocampal signatures. These supply the "−" pole for connectivity scoring via bi-directional weighted-KS.

**Likelihood model under the mixture.** For each modality m ∈ {L1000, CP, MEA}:

P_m(signature(C) | C is cognition-modulator) = Σ_k π_k · 𝒩(signature(C); μ_k, Σ_k)

with π_k learned (Dirichlet prior over MoA classes — empirical Bayes from the anchor distribution; default uniform 1/5). The per-modality likelihood is the marginal over k. This naturally captures "compound matches one specific cognition MoA strongly," which a single Gaussian would understate.

**Remyelination is a separate sub-centroid (not a single feature).** Per Mei 2014 (8 anti-muscarinics including clemastine, benztropine), Najm 2015 (miconazole, clobetasol; supplement has RNA-seq + global phosphoproteomics of OPCs at 1 h and 5 h), Lariosa-Willingham 2016 (27 NCC hits including clobetasol, halcinonide, fluticasone, betamethasone), and Poon 2024 (PIPE-307): build an OPC-differentiation transcriptomic signature from the Najm 2015 supplement and from Lariosa-Willingham 2016 follow-up. There is no published OPC-specific L1000 signature in the core LINCS L1000 release (NEU/NPC lines exist but not staged OPC differentiation), so the Najm 2015 RNA-seq is the canonical transcriptomic bridge.

---

### 6. Integration Math: πphen as Third Independent Bayesian Factor

**Setup.** Let H = "compound C modulates cognition." The posterior given three evidence axes E_target, E_genetic, E_phen is:

P(H | E_target, E_genetic, E_phen) ∝ P(H) · P(E_target, E_genetic, E_phen | H)

**Conditional independence assumption.** Under the working hypothesis that *given H is true*, the three modalities are independent noisy observations of the same underlying property:

P(E_t, E_g, E_p | H) ≈ P(E_t | H) · P(E_g | H) · P(E_p | H)

**Why this is defensible:**
- E_target is structural (binding affinity, predicted DTI) — measures whether C engages a cognition-relevant *protein*.
- E_genetic is human-genetics-derived (OT-Genetics L2G, AHBA spatial expression, cellxgene single-cell context) — measures whether the *target* is causally connected to cognition in humans.
- E_phen is phenotypic — measures whether C produces a *cellular phenotype* matching cognition-modulators, regardless of how.

Each axis sees a different *physical level* of the causal chain (molecule → protein → biology → phenotype). Conditional on the chain being a real cognition mechanism, the measurement noises are largely independent.

**Where independence breaks (and how to fix it).**
- L1000 and OT-Genetics both depend on transcriptional response — a compound whose target has a strong eQTL footprint will be picked up by both. Quantify residual correlation on a held-out anchor set; if Pearson(logit π_genetic, logit π_phen) > 0.3, introduce a Bayesian-network correction term with a learned Gaussian copula on the logits, or equivalently a learned interaction weight in the log-odds combiner.
- E_target and Cell Painting can correlate when a compound is cytotoxic and reduces general protein binding — control by including a cytotoxicity sub-centroid in the cognition reference and *subtracting* its match score.

**Per-modality likelihood with calibrated probabilities.** For each modality m, raw signal s_m is mapped:

p_m = P(E_m | H) / [P(E_m | H) + P(E_m | ¬H)]

Calibration sources:
- **L1000 connectivity:** s = max-quantile NCS across CNS-relevant cells; isotonic regression on anchor compounds. Empirically τ > 90 maps to p ≈ 0.80–0.90 if the anchor library is well-curated.
- **Cell Painting:** s = cosine similarity (or 1 − Mahalanobis quantile) to cognition mixture centroid; Platt scaling.
- **MEA (measured or imputed):** s = 1 − Mahalanobis(observed_15dim, reference_15dim) / χ²₁₅-quantile; Platt.
- **chemCPA-imputed signatures:** apply the modality-specific calibration *after* computing the connectivity/similarity to the imputed signature, and additionally widen the calibrated probability by ensemble σ²_epi.

**Final combination — log-odds form.**

logit P(H | E) = logit P(H) + Σ_m w_m · logit p_m + Σ_{m,m'} ε_{m,m'} · logit p_m · logit p_{m'}

where w_m ∈ ℝ⁺ are learned reliability weights (start at 1; learn via L-BFGS on the anchor set with a Bayesian logistic regression and a horseshoe prior on weights to handle small-N). Cross-terms ε_{m,m'} are small (~0.05) and handle measured non-independence; shared across held-out folds.

**Mixture-of-experts framing for missing modalities.** When one modality is missing (e.g., compound not in JUMP-CP), don't impute — set w_m = 0 *and* widen the posterior:

Var(logit P(H | E)) ← Var(logit P(H | E)) + κ² · (number of missing modalities)

with κ calibrated such that "all three present, all agreeing strongly" produces conformal intervals ~0.05 wide and "two missing" produces ~0.25 wide.

**Uncertainty propagation from chemCPA.** When E_phen is built on a chemCPA-imputed signature:

logit p_phen ~ 𝒩(logit p̂_phen, σ²_ens + σ²_calib)

Integrate analytically (logit-normal → posterior via 7-node Gauss-Hermite quadrature). Fast on CPU per compound.

**Calibration with <100 anchors — Hierarchical Bayes + LOOCV.** With 40–80 pro-cognitive anchors and ~5,000 implicit negatives (random ChEMBL sample with no neuro annotation):
- Hierarchical Bayes with class-level priors over the 5 MoA sub-centroids (sharing strength across anchor groups within an MoA family).
- Leave-one-class-out (LOCO) on the 5 MoA classes — train calibration on 4, validate on 1.
- Wrap with **split conformal prediction**, conditioning per chemical class via Mondrian conformal (`crepes`) for marginal coverage that is also valid per class.

---

### 7. Disagreement Propagation Extension (Cluster D Jensen–Shannon, generalized to 3 factors)

**Original framework (recap).** Cluster D defines a target-first disagreement signal as JS divergence between per-signal posterior distributions of the V4–V6 target-first models. High JS = high model disagreement = lower confidence.

**Three-way JS divergence.** For three Bernoulli posteriors π_t, π_g, π_p over H ∈ {0,1}:

JS_3(π_t, π_g, π_p) = H((π_t + π_g + π_p)/3) − (H(π_t) + H(π_g) + H(π_p)) / 3

where H is the binary entropy. JS_3 ∈ [0, log 3 ≈ 1.0986 nats]. Equivalently, define Bernoulli p_m = Bern(π_m) and compute:

JS_3 = (1/3) Σ_m KL(p_m ‖ p̄), where p̄ = (p_t + p_g + p_p)/3

(The mmJSD formulation of Sutter et al. 2020 NeurIPS justifies this multi-distribution JSD as a proper multimodal divergence.)

**Why this generalizes Cluster D when one factor is mechanism-agnostic.** The mechanism-agnosticism of π_phen is captured by the fact that its disagreement with π_t and π_g is informative in a different direction: in the original target-first framework, all three signals attempt to measure "is C hitting a cognition-relevant target?" — disagreement was always model uncertainty. With π_phen, disagreement can now be either (a) model uncertainty (one or both target-first models wrong) or (b) **genuine off-target / novel mechanism** (target-first models are right that the known target isn't hit, but the phenotype is real because of an unknown target).

**The disagreement-as-signal table (the headline value).** Binarized at threshold 0.7:

| π_t | π_g | π_p | Interpretation | Action |
|---|---|---|---|---|
| H | H | H | Concordant — strong target-mediated cognition | Top priority, low risk |
| H | H | L | Phenotypic failure — wrong cell context or cytotoxicity | Deprioritize; check JUMP cytotox flag |
| H | L | H | Phenotype real, target plausible, genetics weak | Likely real, genetics lag the biology |
| L | H | H | Phenotype + human genetics agree; target call missed (novel target or polypharmacology) | **High-value novel target candidate** |
| H | L | L | Target binding without phenotype — likely irrelevant target or weak engagement | Reject |
| L | H | L | Genetics-only — no biology yet | Deprioritize unless target tractable |
| L | L | H | **The headline cell:** phenotype-only — novel-mechanism candidate (clemastine-class) | Flag for deep mechanistic dive |
| L | L | L | All evidence against | Reject |

The (L, L, H) cell is exactly the cell the existing V4–V6 stack is blind to and is the principal scientific justification for adding πphen.

**Mutual information as a novel-mechanism score.**

I_novel(C) = π_p(C) · [1 − I(π_p ; (π_t, π_g))]

where the MI is estimated over the candidate library. I_novel ≈ 1 when π_p is high and unpredictable from the target-first axes — the operational definition of "novel mechanism candidate." Estimate I with a k-NN MI estimator (Kraskov-Stögbauer-Grassberger 2004; k = 3) which is robust for 5,000-compound libraries.

**Confidence widening from disagreement.**

confidence(C) = P(H | E) · exp(−α · JS_3 / log 3)

with α calibrated so all-three-agree compounds retain ≥95% of posterior, all-three-maximally-disagree retain ~30%. Single confidence score for ranking; I_novel(C) reported orthogonally as a category label.

**Three-way JSD vs three pairwise JSDs.** Three-way JSD is information-theoretically optimal for symmetric disagreement; pairwise JSDs (JS(π_t,π_g), JS(π_t,π_p), JS(π_g,π_p)) carry additional structure useful for diagnosing *which two axes disagree*. Report both — three-way for ranking attenuation, pairwise for explainability.

---

### 8. Implementation Roadmap and Compute/Storage Estimates

**Phase 0 (week 0–2): infrastructure.** Snakemake pipeline scaffolding, DVC for data versioning (data layer separate from code), conda/mamba environments per modality. Total: ~100 GB local SSD, no GPU.

**Phase 1 (week 2–5): L1000 axis.** Download GSE92742 + GSE70138 GCTX → 50 GB. Build per-cell-line reference indices in signatureSearch's HDF5/SQLite or duckdb backend. Implement bi-directional weighted-KS, NCS, tau, NCSct. Compute cognition reference signatures (5 sub-MoAs + AD/SCZ/MCI inverse anchors). Output: per-compound τ and calibrated p_L1000 for ~30,000 LINCS-resident compounds. CPU only.

**Phase 2 (week 5–9): JUMP-CP axis.** `aws s3 sync` only the parquet profiles you need (~30 GB working set; full set 500 GB–2 TB if you cache all sources). Pycytominer pipeline: annotate → robustize → spherize → feature-select → harmony. Mahalanobis to cognition mixture. Output: per-compound p_CP for ~116k JUMP compounds. CPU; ~2 days end-to-end on a workstation.

**Phase 3 (week 9–12): MEA axis (partial).** Download ToxCast `invitroDB` (~10 GB) via `tcpl` (R). Extract Shafer MEA acute (1,055 compounds) + Kosnik concentration-response (384 × 43). Build 15-dim reduced reference. Compute p_MEA for compounds with primary data only. CPU; minimal storage.

**Phase 4 (week 12–18): chemCPA / biolord ensembles for MEA + L1000 imputation.** Train CPA on L1000 (RTX 5070, ~6–10 h per seed) × 5–10 seed ensemble. For MEA imputation, train a regression head (chemical_encoder → 15-dim MEA vector) on Kosnik's 384 compounds with InChIKey-paired chemical encodings (Grover or MorganFP). Mondrian conformal (`crepes`) per chemical scaffold class. ~3–5 GPU-days total.

**Phase 5 (week 18–22): integration + disagreement.** Bayesian logistic regression with horseshoe prior (PyMC or NumPyro on CPU); LOCO + Mondrian conformal calibration. Implement three-way JSD, I_novel, pairwise JSDs. Output: per-compound (π_t, π_g, π_p, π_joint, JS_3, I_novel, conformal CI).

**Phase 6 (week 22+): validation.** Held-out pro-cognitive recovery, MoA-cluster enrichment, retrospective recovery of clemastine / miconazole / PIPE-307 from the (L, L, H) cell when their target annotations are blinded.

**Total storage:** ~1–3 TB working set. External NAS or 4 TB NVMe sufficient.
**Total compute:** single workstation with RTX 5070 (12 GB) + 64 GB RAM + 16 cores adequate. No cloud GPU required.

---

### 9. Software / Library Recommendations

| Layer | Library | Purpose |
|---|---|---|
| L1000 parsing | cmapPy (Python) / cmapR (R) / cmapM (Matlab) | GCTX I/O (Enache et al. 2020) |
| Connectivity scoring | signatureSearch (R/Bioconductor; Duan/Girke 2020 BMC Bioinformatics) | LINCS/CMAP/Fisher/Cor algorithms (Subramanian 2017) |
| Signature search | SigCom LINCS REST API (Evangelista 2022); iLINCS + drugfindR (R; Imami 2024) | Web-scale signature commons; ad-hoc queries |
| Cell Painting | pycytominer (cytomining/pycytominer) | annotate, normalize, feature_select, sphere |
| Cell Painting eval | copairs (Broad) | mAP-based profile matching benchmark |
| Cell Painting batch correction | Harmony via scanpy.external.pp.harmony_integrate | site-level integration |
| Cell Painting DL embeddings | DeepProfiler; EfficientNet-V2-S features from `workspace_dl/` | pretrained imaging features |
| MEA acquisition | Axion AxIS Navigator + Neural Metric Tool | MEA feature extraction |
| MEA toxicology | tcpl (R) | ToxCast invitroDB R-side access |
| Generative perturbation | theislab/cpa (CPA); chemCPA; biolord (scvi-tools-compatible); PRnet | predicted signatures for OOD compounds |
| Bayesian integration | PyMC / NumPyro / Pyro (Python) | hierarchical posterior, horseshoe priors |
| Calibration | scikit-learn IsotonicRegression / Platt | per-modality probability calibration |
| Conformal prediction | crepes (Python) — Mondrian split CP | finite-sample coverage |
| Workflow | Snakemake (preferred) or Nextflow | DAG orchestration |
| Data versioning | DVC (lightweight) or Pachyderm (heavier) | reproducibility |

---

### 10. Specific Answers to Technical Sub-Questions

**Q1: Best way to construct cognition reference?** Mixture of K=5 MoA sub-centroids (cholinergic / catecholaminergic / glutamatergic / trophic-ISR / remyelination) + 3 disease-inverse anchors (AD, SCZ, MCI), with empirical-Bayes Dirichlet priors over class weights. Single centroid loses MoA structure; full per-compound is overfitting risk.

**Q2: Remyelination signatures?** No standalone OPC-staged signature exists in core LINCS L1000. Use Najm 2015 Nature 522:216-220's RNA-seq supplement (miconazole + clobetasol vs DMSO in OPCs, 1 h + 5 h) as the canonical transcriptomic OPC-differentiation signature. For Cell Painting, no published OPC profile exists in JUMP; build a *predicted* OPC-differentiation profile from anchor compounds (clemastine + Mei 2014's BIMA-8 + miconazole + clobetasol) in U2OS as a proxy MoA cluster, and validate against PIPE-307 (Poon 2024 PNAS).

**Q3: Aggregating connectivity across cell lines?** Use NCSct (max-quantile statistic from Subramanian 2017): per-cell NCS, then the larger absolute of 67th vs 33rd percentile across cells. For CNS focus, additionally weight NPC/NEU/ASC by ~3× others.

**Q4: Aligning Cell Painting and L1000 features?** Per Way 2022 Cell Systems 13(11):911-923, ~3,000 compounds are dual-profiled across Bray 2017 + JUMP-CP. Train CCA (or a deep cross-modal contrastive embedding) on the paired set; project both into a joint 64-128-dim latent. Use the joint latent for compounds present in only one modality.

**Q5: SOTA for predicting MEA from structure?** Limited prior art; most papers use molecular descriptors → Random Forest / GBM regression on individual MEA metrics (Kosnik 2020 used hierarchical clustering, not prediction). For SOTA-grade: train a GNN (Chemprop / Grover / MolFormer) on Strickland 2018 + Kosnik 2020 (~1,400 compound-condition pairs) → multi-output regression to the 15-dim Kosnik reduced set; wrap with conformal prediction. Expect modest OOD performance (r² 0.3-0.5) given dataset size.

**Q6: Epistemic uncertainty in chemCPA?** Critical: the built-in "uncertainty" is a latent-distance heuristic, not calibrated variance. Wrap in 5-10 model deep ensemble (different seeds, optional bootstrap). Report ensemble mean as prediction, ensemble variance as σ²_epi, decoder σ² as σ²_alea. Calibrate σ²_total via Mondrian split conformal per scaffold class. Compare against biolord (Piran 2024) as alternative (mean r² 0.76 vs chemCPA's 0.51 on held-out sci-Plex3).

**Q7: Published Bayesian L1000 + Cell Painting integration?** Rare. Way 2022 Cell Systems 13(11):911-923 quantifies complementarity (L1000 better for MAPK/HSP, CP better for Aurora/PLK/BRD4) but uses mAP/connectivity, not Bayesian factor combination. PRnet (Qi 2024 Nat Commun) is closest to predicted-profile-based screening but point-estimate only. There is no published explicit naive-Bayes / hierarchical-Bayes joint posterior over L1000 + CP for drug repurposing — this πphen design *is* methodologically novel.

**Q8: Dose-response.** L1000 is multi-dose (typically 10 μM + 100 nM, 6 h + 24 h); JUMP-CP is single-dose (10 μM, 48 h). Solutions:
- L1000: aggregate by highest-TAS signature; or fit a Hill curve to NCS vs dose for compounds with ≥3 doses and report NCS_max.
- JUMP-CP: take the single profile, but flag compounds with non-monotonic dose-response in L1000 as "biphasic" — these are the high-risk false positives.
- MEA: use Kosnik concentration-response when available; ToxCast acute (40 μM single dose) gives lower precision and should be downweighted.

---

## Recommendations

**Stage 0 (now, 1 week):** Stand up cmapPy + signatureSearch and reproduce a known L1000 cognition query — e.g., donepezil's top-τ neighbors in NPC/NEU. This is a sanity check; if you can't recover known cholinesterase-pathway hits in NPC, the cell-line filtering is wrong before you build anything else.

**Stage 1 (weeks 2–6): L1000 axis only.** Get p_L1000 in production. Highest information-per-compute, largest reference library, exercises every integration-math component except chemCPA. LOOCV validation on pro-cognitive anchors.

**Stage 2 (weeks 6–12): add JUMP-CP axis.** Once L1000 calibration is reliable, add Cell Painting. The Way 2022 complementarity benchmark gives an empirical test of independent information. If correlation(logit p_L1000, logit p_CP) on the anchor set > 0.6, your CP cognition reference is overfitting to transcriptional pathways; rebuild it as a pure morphological centroid.

**Stage 3 (weeks 12–18): add chemCPA imputation + MEA axis.** Engineering-heavy phase: training an ensemble, doing MEA prediction, *and* calibrating two new probability streams. Don't attempt until Stages 1 and 2 are stable.

**Stage 4 (weeks 18–22): activate disagreement-as-signal.** Three-way JSD, I_novel, pairwise JSDs. Validate that the (L, L, H) cell recovers clemastine, miconazole, PIPE-307 retrospectively when target annotations are blinded.

**Stage 5 (week 22+): wet-lab triage.** Iterate on the highest-I_novel candidates and ship to an external CRO (Charles River, Axion services, Eurofins) for iPSC-MEA confirmation on top hits.

**Thresholds that should change the plan:**
- If JUMP-CP p_CP and L1000 p_L1000 correlate > 0.6 on the anchor set → πphen collapses to one effective axis; reformulate as a single transcriptional-morphological factor.
- If chemCPA ensemble σ²_epi swamps σ²_alea by >5× on a held-out anchor → predicted signatures are unusable for Bayesian integration; drop chemCPA-imputed compounds and rely only on directly-profiled ones.
- If LOCO calibration AUC < 0.7 on any single MoA sub-class → that sub-centroid is mis-anchored; rebuild with literature curation.
- If three-way JSD distribution is bimodal at log 3 → axes are systematically discordant (calibration error), not informatively disagreeing; recalibrate isotonically.

---

## Caveats

1. **Cell-context mismatch is the dominant systematic risk.** L1000 has ~5% CNS-line coverage; JUMP-CP is U2OS / A549; ToxCast MEA is rat primary cortical; iPSC MEA from Axion is human glutamatergic but with limited compound coverage. None are cortical-specific cognition models. The πphen axis measures "phenotypic agreement with anchor compounds in available cell systems," not "in vivo cognition phenotype."

2. **chemCPA's advertised "uncertainty" is a latent-distance heuristic, not calibrated variance.** Per Lotfollahi 2023 MSB explicitly: *"the current heuristic for uncertainty estimation originates from the compositional formulation in CPA"* and *"we opted to implement a deterministic autoencoder scheme."* You will get garbage uncertainty estimates unless you wrap with deep ensembles.

3. **Independence assumption is approximate.** L1000 and OT-Genetics both observe transcriptional consequences. Measure residual correlation on the anchor set; correct with Gaussian copula or interaction terms when |ρ| > 0.3. Test it on every pipeline release.

4. **Small anchor set inflates calibration variance.** With 40-80 pro-cognitives and ≤8 well-characterized remyelinators, isotonic regression is high-variance. Use hierarchical Bayes with class-level priors and report calibration uncertainty (95% credible band on the reliability diagram), not point estimates.

5. **No published Bayesian L1000 + CP integration exists to compare against.** Validation must be retrospective recovery of held-out cognition anchors + MoA-cluster enrichment + ideally a small prospective wet-lab cohort. Set absolute coverage/precision targets from first principles.

6. **The (L, L, H) disagreement cell is also where Cell Painting false positives live.** Cytotoxicity, mitotic disruption, and generic cytoskeletal effects produce high CP similarity because the reference is in U2OS. Always intersect (L, L, H) hits with a JUMP-CP cytotoxicity / cell-count filter and a Cell Health (Way 2021) model.

7. **Many MEA datasets are commercial.** Charles River, Axion services, FUJIFILM CDI, Eurofins run paid MEA screens with proprietary nootropic libraries. Publishable academic MEA cognition data is genuinely scarce; budget for commercial-data licensing or original wet-lab work to ground-truth the MEA axis.

8. **Conformal coverage holds marginally, not conditionally, unless Mondrian.** Wrap per chemical class; otherwise rare classes (e.g., ampakines) will have poor coverage hidden by good aggregate stats.