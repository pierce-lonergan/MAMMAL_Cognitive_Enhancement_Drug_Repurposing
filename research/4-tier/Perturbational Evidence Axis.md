# V8 — Phenotypic / Perturbational Evidence Axis
**Architecture Specification for the MOMENTUM-X Cognition-Repurposing Pipeline**

## TL;DR
- **V8 is the integration test to V6.A/B's unit tests.** A target-agnostic, multi-modal phenotypic prior (LINCS L1000 + JUMP-CP + iPSC-neuron MEA/snRNA-seq + chemCPA-imputed signatures) is jointly fused with V6.A binding, V6.B target relevance, and V7 effect size in a single conditionally-dependent hierarchical Bayesian model. Rating: **A+ architectural novelty, A− methodological risk, B+ wet-validation risk (mitigated by pre-registration).** The U2OS-to-brain transfer is the single biggest theoretical leap and is handled by an explicit cell-line random effect plus a `λ_phen` weight that degrades gracefully when Gate 1 fails.
- **Validation is gated and pre-registered on OSF before unblinding.** Gate 1 (mechanism-class recovery, AMI > 0.5 / ARI > 0.4) is *primary* — it tests whether the phenotypic prior is real signal at all. Gate 2 (held-out Hedges' *g* prediction, MAE < 0.20) is *secondary*. Gate 3 (nootropic-anchor connectivity sanity check, ≥ 7/9 expected nearest-neighbors) is a *go/no-go sanity*. Failures degrade weighting but do not necessarily kill the axis — three-tier PASS/DEGRADE/FAIL bands.
- **Publishability target: Nature Machine Intelligence (A, realistic); Nature Methods (A+, stretch contingent on Gate 1 AMI ≥ 0.6).** The novelty claim is narrow and defensible: *first* multi-modal phenotypic prior (transcriptomic + morphological + electrophysiological + generative-imputed) coupled to a target-first Bayesian repurposing posterior with explicit conditional-dependence structure and four-axis disagreement-as-signal facet tagging.

---

## Key Findings

1. **The structural gap V8 closes is real and named.** V6.A (Multi-Head DTI) and V6.B (Cluster D neurobiological prior) are both *target-first*: they presuppose a named target. They cannot see (i) polypharmacology, (ii) cryptic off-target activity, (iii) novel-scaffold mechanism-of-action, or (iv) cellular-state effects that no single target captures. The Lamb 2006 *Science* / Subramanian 2017 *Cell* CMap/L1000 paradigm and the Bray 2016 *Nat. Protoc.* / Chandrasekaran 2024 *Nat. Methods* Cell Painting paradigm both demonstrate that signature-matching recovers mechanism-of-action *without target hypothesis*. V8 imports that capability and constrains it Bayesianly with the V6 priors.

2. **Six-modality scope is justified, but unequal in weight.** L1000 (Subramanian et al. *Cell* 2017) provides **1,319,138 L1000 profiles from 42,080 perturbagens (19,811 small-molecule compounds, 18,493 shRNAs, 3,462 cDNAs, and 314 biologics), consolidated into 473,647 replicate-collapsed signatures across 9 core cell lines for systematic profiling** — the transcriptomic workhorse. JUMP-CP cpg0016 (Chandrasekaran et al. *Nat. Methods* 2024) supplies **136,000 chemical and genetic perturbations spanning 116,750 unique compounds, over-expression of 12,602 genes, and knockout of 7,975 genes by CRISPR-Cas9, all in human osteosarcoma (U2OS) cells** — the morphological prior. iPSC-neuron MEA + snRNA-seq from Brennand, Pașca, Studer, and FUJIFILM CDI lines is the neuron-context prior — patchy but irreplaceable. cellxgene-census brain-organoid slice is a normalization scaffold. PsychENCODE / BICCN / BrainSpan supply reference brain transcriptomes. chemCPA (Hetzel et al. NeurIPS 2022) extends coverage to the full ~50 K-candidate × ~210-target library via generative imputation.

3. **Joint posterior must be conditionally dependent, not independent.** Phenotype is downstream of binding and downstream of cognition-relevance; pretending it is independent and multiplying likelihoods is statistically wrong and produces overconfidence. The correct decomposition is `p(θ_B, θ_T, θ_P, θ_E | data) ∝ p(θ_T) · p(θ_B | θ_T) · p(θ_P | θ_B, θ_T) · p(θ_E | θ_B, θ_T, θ_P, PBPK, moderators)` — a four-level hierarchical model where phenotype is the *integration test* fusing binding × target into a single observable cellular phenotype, and effect size is the final translational layer.

4. **chemCPA is the load-bearing methodological choice and its failure modes are known.** The Hetzel 2022 RDKit-pretrained chemCPA on the 9-OOD sci-Plex3 benchmark achieves mean R² ≈ 0.69 (all genes) / 0.47 (DEGs) at the 10 µM dose (Supp. Tables 6–7); the more conservative, externally-cited benchmark from the Biolord paper (Piran et al. *Nat. Biotechnol.* 2024) reports **chemCPA-pre mean R² = 0.51 ± 0.0062 vs. Biolord 0.76 ± 0.0005** averaged across cell lines and dosages. This is *good enough* to function as a downweighted prior but *not good enough* to treat as a measurement. Imputation uncertainty must propagate into the joint posterior via a chemCPA-specific noise inflation term.

5. **The mechanism-discovery clustering validation is the right primary gate.** AMI > 0.5 / ARI > 0.4 thresholds against the V7 PRISMA ~30-class mechanism taxonomy are a strong test because the clustering is *unsupervised* on phenotypic embeddings while the labels are *external* (chemistry + literature, not phenotype). This is the principled answer to the "phenotype is noise" reviewer objection.

6. **The Gate-2 effect-size prediction will be statistically thin.** Healthy-adult cognition Hedges' *g* is small: Roberts CA, Jones A, Sumnall H, Gage SH, Montgomery C, *Eur. Neuropsychopharmacol.* 38:40–62 (2020), reports overall modafinil **SMD = 0.12 (p = .01)**, overall methylphenidate **SMD = 0.21 (p = .0004)**, and the maximum significant subdomain effect of **SMD = 0.43 (p = 0.0002)** for MPH on delayed recall — that 0.43 is the realistic ceiling. With ~50–100 anchor compounds and ceiling effects this low, regression signal is constrained. MAE < 0.20 is reasonable; 90% CrI coverage ≥ 85% is the more defensible metric.

7. **Encenicline is the canonical "phenotypic prior saves us from a target-true failure" case.** α7 nAChR partial agonist with valid V6.A binding; Phase 2 (Keefe et al. *Neuropsychopharmacology* 40:3053–3060, 2015) gave Cohen's d = 0.257 in the 0.27 mg arm (p = 0.034) and 0.093 in the 0.9 mg arm (p = 0.255), and Phase 3 missed co-primary endpoints in two schizophrenia trials per FORUM Pharmaceuticals' March 2016 topline announcement (Brannan et al. *Schizophr. Bull.* 45(S2):S141 abstract; Alzforum therapeutics database). Phenotypic signature in LINCS/JUMP is expected to look inert relative to active enhancers. Intepirdine MINDSET (Lang FM, Mo Y, Sabbagh M, et al. *Alzheimer's Dement.: TRCI* 7(1):e12136, 2021) reported in *n* = 1,315: ADAS-Cog adjusted mean difference −0.36 (p = 0.225) and ADCS-ADL −0.09 (p = 0.826). Pridopidine PROOF-HD (Reilmann R et al. *Nat. Med.* 2025, DOI 10.1038/s41591-025-03920-3): TFC LS mean diff −0.18 (p = 0.26), cUHDRS −0.11 (p = 0.45). These are V8 test cases for `target-true.phenotype-failed` facet detection.

8. **U2OS-to-brain is the elephant.** Addressed by (a) per-cell-line random effect in the joint embedding; (b) explicit downweighting of cell-line-only phenotypic evidence in iPSC-deficient subspaces; (c) Gate 1 stratification reporting per-modality contribution; (d) honest framing in Limitations.

---

## Details

### A. Literature Foundation

#### A.1 Connectivity Map lineage
- **Lamb J, Crawford ED, Peck D, Modell JW, Blat IC, Wrobel MJ, Lerner J, Brunet J-P, Subramanian A, Ross KN, Reich M, Hieronymus H, Wei G, Armstrong SA, Haggarty SJ, Clemons PA, Wei R, Carr SA, Lander ES, Golub TR.** *Science* 313(5795):1929–1935 (2006). DOI: 10.1126/science.1132939. The original 164-perturbagen × 3-cell-line CMap establishing signature-matching for MOA discovery.
- **Subramanian A, Narayan R, Corsello SM, Peck DD, Natoli TE, Lu X, Gould J, Davis JF, Tubelli AA, Asiedu JK, Lahr DL, Hirschman JE, Liu Z, Donahue M, Julian B, Khan M, Wadden D, Smith IC, Lam D, Liberzon A, Toder C, Bagul M, Orzechowski M, Enache OM, Piccioni F, Johnson SA, Lyons NJ, Berger AH, Shamji AF, Brooks AN, Vrcic A, Flynn C, Rosains J, Takeda DY, Hu R, Davison D, Lamb J, Ardlie K, Hogstrom L, Greenside P, Gray NS, Clemons PA, Silver S, Wu X, Zhao W-N, Read-Button W, Wu X, Haggarty SJ, Ronco LV, Boehm JS, Schreiber SL, Doench JG, Bittker JA, Root DE, Wong B, Golub TR.** *Cell* 171(6):1437–1452.e17 (2017). DOI: 10.1016/j.cell.2017.10.049. Introduces τ, WTCS, FDR-based connectivity scoring across 1.3 M L1000 profiles.
- **Sirota M, Dudley JT, Kim J, Chiang AP, Morgan AA, Sweet-Cordero A, Sage J, Butte AJ.** *Sci. Transl. Med.* 3(96):96ra77 (2011). Signature-reversal repurposing — cimetidine for NSCLC.
- **Zhang S-D, Gant TW.** *BMC Bioinformatics* 9:258 (2008). sscMap (related to the Iorio 2010 connectivity-score methodology referenced in the design brief).
- **Cheng J, Yang L, Kumar V, Agarwal P.** *Genome Med.* 6(12):95 (2014). Systematic CMap evaluation.

#### A.2 Cell Painting lineage
- **Bray M-A, Singh S, Han H, Davis CT, Borgeson B, Hartland C, Kost-Alimova M, Gustafsdottir SM, Gibson CC, Carpenter AE.** *Nat. Protoc.* 11(9):1757–1774 (2016). DOI: 10.1038/nprot.2016.105. Six dyes, five channels, eight cellular components, ~1,500 features.
- **Cimini BA et al.** *Nat. Protoc.* 18:1981–2013 (2023). Optimized Cell Painting v2.5.
- **Caicedo JC et al.** *Nat. Methods* 14(9):849–863 (2017). Data-analysis strategies for image-based profiling.
- **Chandrasekaran SN, Ackerman J, et al.** *Nat. Methods* 21(6):1114–1121 (2024). DOI: 10.1038/s41592-024-02241-6. CPJUMP1 and cpg0016 specifications.
- **Weisbart E, Kumar A, Arevalo J, Carpenter AE, Cimini BA, Singh S.** *Nat. Methods* 21:1775–1777 (2024). Cell Painting Gallery — **"As of May 2024, the Cell Painting Gallery holds 688 terabytes (TB) of image and associated numerical data."**
- **Moshkov N, Bornholdt M, Benoit S, Smith M, McQuin C, Goodman A, Senft RA, Han Y, Babadi M, Horvath P, Cimini BA, Carpenter AE, Singh S, Caicedo JC.** *Nat. Commun.* 15:1594 (2024). DOI: 10.1038/s41467-024-45999-1. DeepProfiler CellPainting_CNN.
- **Sypetkowski M et al.** Phenom-1 / OpenPhenom / DINOv2-for-Cell-Painting baselines (Recursion 2024).
- **Hofmarcher M et al.** *J. Chem. Inf. Model.* 59:1163–1171 (2019). Deep features for Cell Painting.

#### A.3 iPSC neuron lineage
- **Brennand KJ, Simone A, Jou J, Gelboin-Burkhart C, Tran N, Sangar S, Li Y, Mu Y, Chen G, Yu D, McCarthy S, Sebat J, Gage FH.** *Nature* 473(7346):221–225 (2011). hiPSC schizophrenia model.
- **Wen Z, Nguyen HN, Guo Z, et al.** *Nature* 515:414–418 (2014). DISC1 iPSC synaptic dysfunction.
- **Hoffman GE, Hartley BJ, Flaherty E, et al. (Brennand lab).** *Nat. Commun.* 8:2225 (2017). hiPSC-NPC/neuron transcriptional signatures concordant with post-mortem brain.
- **Hartley BJ, Brennand KJ et al.** *PNAS* 119(11):e2109395119 (2022). Important correction: this paper used **patch-clamp electrophysiology (not MEA)** on hiPSC-derived cortical neurons; it found that sEPSC amplitude predicted Wisconsin Card Sorting Test (WCST) performance in SCZ patients with **R = −0.93, P = 2.94 × 10⁻⁵**. MEA results from other groups (Frank, Odawara, Hyysalo) constitute the complementary MEA-based literature.
- **Paşca SP.** *Nature* 553:437–445 (2018). Cortical organoid review.
- **Yoon S-J et al. (Paşca lab).** *Nat. Methods* 16:75–78 (2019). Reliable cortical organoid generation.
- **Revah O et al. (Paşca lab).** *Nature* 610:319–326 (2022). Transplanted cortical organoid maturation.
- **Frank CL, Brown JP, Wallace K, Mundy WR, Shafer TJ.** *Toxicol. Sci.* 160:121–135 (2017). MEA for developmental neurotoxicity.
- **Hyysalo A et al.** *Stem Cell Res.* 24:118–127 (2017). hiPSC-neuron MEA characterization.
- **Odawara A, Katoh H, Matsuda N, Suzuki I.** *Sci. Rep.* 6:26181 (2016). Physiological maturation on MEAs.
- **Studer lab** directed differentiation protocols (Kriks et al. 2011; Chambers et al. 2009 dual-SMAD).
- **FUJIFILM Cellular Dynamics** iCell GlutaNeurons / iCell DopaNeurons / iCell GABA Neurons; **AxoSim**, **NeuCyte**, **BrainXell** commercial datasets on Synapse.org / PsychENCODE.

#### A.4 Generative perturbation modeling
- **Hetzel L, Böhm S, Kilbertus N, Günnemann S, Lotfollahi M, Theis FJ.** NeurIPS 2022. chemCPA.
- **Lotfollahi M, Klimovskaia Susmelj A, De Donno C, Hetzel L, Ji Y, Ibarra IL, Srivatsan SR, Naghipourfar M, Daza RM, Martin B, Shendure J, McFaline-Figueroa JL, Boyeau P, Wolf FA, Yakubova N, Günnemann S, Trapnell C, Lopez-Paz D, Theis FJ.** *Mol. Syst. Biol.* 19:e11517 (2023). CPA.
- **Piran Z, Cohen N, Hoshen Y, Nitzan M.** *Nat. Biotechnol.* 42(11):1678–1683 (2024). DOI: 10.1038/s41587-023-02079-x. Biolord; reports the canonical externally-benchmarked figure **chemCPA-pre 0.51 ± 0.0062 vs. Biolord 0.76 ± 0.0005**.
- **Lotfollahi M, Wolf FA, Theis FJ.** *Nat. Methods* 16:715–721 (2019). scGen.
- **Roohani Y, Huang K, Leskovec J.** *Nat. Biotechnol.* 42:927–935 (2024). GEARS.
- **Lopez R, Regier J, Cole MB, Jordan MI, Yosef N.** *Nat. Methods* 15:1053–1058 (2018). scVI.
- **Xu C, Lopez R, Mehlman E, Regier J, Jordan MI, Yosef N.** *Mol. Syst. Biol.* 17:e9620 (2021). scANVI.

#### A.5 Multi-view learning
- **Andrew G, Arora R, Bilmes J, Livescu K.** Deep Canonical Correlation Analysis. *ICML* PMLR 28(3):1247–1255 (2013). (ICML 2023 Test-of-Time runner-up.)
- **Argelaguet R, Arnol D, Bredikhin D, Deloro Y, Velten B, Marioni JC, Stegle O.** *Genome Biol.* 21:111 (2020). MOFA+.
- **Argelaguet R, Velten B, Arnol D, Dietrich S, Zenz T, Marioni JC, Buettner F, Huber W, Stegle O.** *Mol. Syst. Biol.* 14:e8124 (2018). MOFA v1.
- **Lock EF, Hoadley KA, Marron JS, Nobel AB.** *Ann. Appl. Stat.* 7:523–542 (2013). JIVE.
- **Lee C, van der Schaar M.** Multi-view VAE (MVAE) for missing-modality learning.

#### A.6 Pre-registration & reference outcomes
- **OSF.io** (Center for Open Science) and **AsPredicted.org** — pre-registration sinks.
- **Roberts CA, Jones A, Sumnall H, Gage SH, Montgomery C.** *Eur. Neuropsychopharmacol.* 38:40–62 (2020). DOI: 10.1016/j.euroneuro.2020.07.002.
- **Keefe RSE et al.** *Neuropsychopharmacology* 40(13):3053–3060 (2015). DOI: 10.1038/npp.2015.176.
- **Lang FM, Mo Y, Sabbagh M, Solomon P, Boada M, Jones RW, Frisoni GB, Grimmer T, Dubois B, Harnett M, Friedhoff SR, Coslett S, Cummings JL.** *Alzheimer's Dement.: TRCI* 7(1):e12136 (2021). MINDSET intepirdine — primary endpoints did not separate from placebo.
- **Reilmann R, McGarry A, Grachev ID, et al.** *Lancet Neurol.* 18(2):165–176 (2019). PRIDE-HD.
- **Reilmann R et al.** *Nat. Med.* 2025 (DOI 10.1038/s41591-025-03920-3). PROOF-HD pridopidine.

---

### B. Methodology — Multi-Modal Phenotypic Prior Architecture

#### B.1 Data harmonization stack (engineering metaphor: build a 4-language ABI)

| Modality | Native unit | Native dim | Harmonization step |
|---|---|---|---|
| L1000 | Level-5 z-score (MODZ) | 12,328 inferred / 978 landmark | Rank-normalize → WTCS query; project to z-scored landmark or BING space |
| JUMP-CP CellProfiler | Float feature vector | ~3,200 (after FS ~700) | Sphering + per-plate normalization + ComBat by source/site |
| JUMP-CP DeepProfiler | EfficientNetB0 block6a_activation | 672-dim single-cell (per Moshkov 2024) / aggregate per-well | Sphering to negative-control covariance; treatment-level consensus |
| JUMP-CP DINOv2 | ViT-S/16 patch tokens | 384 (or 768 pooled) | Same sphering pipeline |
| iPSC-neuron MEA | Spike-train metrics | ~25 (MFR, WMFR, BR, NBF, NBD, ISI CV, IBI, synchrony) | Log-transform, per-plate baseline subtraction, donor random effect |
| iPSC-neuron snRNA-seq | UMI counts | ~30 K genes | scVI/scANVI integration with cellxgene-census brain organoid prior; pseudobulk per-compound |
| chemCPA latent | Latent perturbation vector | 128–256 | Trained jointly on L1000 + JUMP-CP via architecture surgery |

**LINCS L1000 ingestion:** GEO **GSE92742** (Phase 1, Level-5 file `GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz` — 473,647 signatures × 12,328 genes) and **GSE70138** (Phase 2, `GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz` — 118,050 sigs × 12,328 genes); the current clue.io beta release is `level3_beta_all_n3026460x12328.gctx`. Access via `cmapPy.pandasGEXpress.parse`; compute WTCS (weighted-connectivity τ) per the Subramanian 2017 STAR Methods. Filter to cell lines with adequate coverage (A375, MCF7, PC3, VCAP, HA1E, HCC515, HEPG2, HT29, plus neural-lineage lines NPC, NEU, SHSY5Y where present). For cognition relevance, **upweight** the neural-lineage cell lines using V6.B target relevance as a prior covariate.

**JUMP-CP ingestion:** S3 sync from `s3://cellpainting-gallery/cpg0016-jump` (no-sign-request, AWS Registry of Open Data). The CPJUMP1 benchmark subset (303 compounds, 160 paired genes per the `jump-cellpainting/2024_Chandrasekaran_NatureMethods` repo) is the first sanity drop. The full compound set (116,750 compounds, ≈ 75 M single cells; cpg0016 estimated at ~115 TB of images + numerical data per Chandrasekaran 2023 bioRxiv) is **not** downloaded raw — pull only pre-computed pycytominer consensus profiles (sphered, feature-selected) and DeepProfiler well-level 672-dim embeddings, totaling tens of GB rather than ~115 TB. Use `pycytominer.normalize` and `pycytominer.feature_select` as canonical loaders. Apply Harmony or ComBat batch correction across sources S1–S13 (source nuisance, plate nested, well-type biological covariate).

**iPSC-neuron MEA ingestion:** Aggregate published Axion Maestro datasets — Frank 2017 *Tox. Sci.*, Odawara 2016 *Sci. Rep.*, Hyysalo 2017 *Stem Cell Res.* — plus FCDI- and BrainXell-derived MEA panels on Synapse.org PsychENCODE. Extract canonical features: mean firing rate (MFR), weighted MFR (WMFR), burst rate, ISI CV, IBI, network burst frequency (NBF), network burst duration (NBD), synchrony index. Per-compound aggregation = median across wells + IQR as uncertainty.

**iPSC-neuron scRNA-seq ingestion:** Pull from cellxgene-census the brain-organoid and iPSC-neuron slices (already loaded for V6.B). Integrate donors and batches via scVI / scANVI (Lopez 2018; Xu 2021). Pseudobulk per (compound × cell-type × time) for the sparse subset of perturbed iPSC neurons. Most iPSC evidence will be **observational baseline** rather than perturbation; treat as a normalization scaffold and donor-variation prior.

**Reference brain transcriptomes:** PsychENCODE (Synapse.org syn21557948 et seq.), BrainSpan, BICCN. Use these for the *cell-line-to-brain* projection — anchor U2OS baseline and U2OS+compound transcriptomic projections on brain reference axes.

#### B.2 chemCPA generative imputation methodology

chemCPA (Hetzel et al. NeurIPS 2022) is selected over GEARS or scGen because it (i) explicitly accepts unseen compound chemical structures via a pretrained molecule encoder, (ii) has well-characterized OOD validation on sci-Plex3, and (iii) supports bulk-to-single-cell architecture surgery — i.e., pretrained on bulk LINCS perturbation effects and fine-tuned on single-cell data, which exactly matches our integration use case.

**Architecture:**

z_basal = E_cell(x_control)

z_pert  = M(G(SMILES)) · S(d)

x̂ = D(z_basal + z_pert + Σ_c E_c(covariate_c))

with G = frozen molecule encoder (**use RDKit Morgan-FP-MLP pretrained**, as the best validated OOD R² in Hetzel 2022 Supp. Tables 6–7); M maps chemical embedding → perturbation latent; S is the amortized dose scaler; E_c are covariate encoders (cell line, time); D is the decoder. Adversarial discriminators on z_basal enforce disentanglement.

**Training objective:**

L = L_recon(x̂, x) + λ_adv · L_adv(disc(z_basal) → pert) + λ_KL · L_KL

**Validation regime (this pipeline):**
- *Random held-out (80/20):* sanity baseline.
- *Scaffold-held-out:* Bemis-Murcko scaffold split, Tanimoto < 0.5 to training set.
- *Leave-one-mechanism-class-out (LOMCO):* train on stimulants + cholinergic; test on glutamatergic. Hardest generalization test.
- *Anchor benchmark:* sci-Plex3, 9-OOD compounds at 10 µM — Hetzel 2022 ceiling: chemCPA-RDKit-pretrained mean R²(all) ≈ 0.69 / R²(DEGs) ≈ 0.47. Externally-benchmarked under broader conditions (Piran et al. *Nat. Biotechnol.* 2024 Biolord paper): **chemCPA-pre mean R² = 0.51 ± 0.0062 vs Biolord 0.76 ± 0.0005**. Target for this pipeline's broader-mix retraining: mean R²(all) ≥ 0.50 and R²(DEGs) ≥ 0.30.

**Decision:** Use chemCPA as the **primary imputer**, with Biolord as a secondary cross-check on the highest-uncertainty compounds. Inflate the chemCPA-imputed signature noise term by a learned scaling factor τ_chemCPA estimated on held-out compounds. Compounds with max-Tanimoto-to-train < 0.3 → flagged `phenotype.imputed.low_confidence` and downweighted.

#### B.3 Multi-modal joint embedding — choice and justification

| Method | Pros | Cons | Rating |
|---|---|---|---|
| Deep CCA (Andrew 2013 ICML) | Maximizes cross-view correlation; clean nonlinear CCA extension | Two-view native; DGCCA extends to many views but is finicky | **B+** |
| MOFA+ (Argelaguet 2020 *Genome Biol.*) | Bayesian, sparse, interpretable factors with per-modality variance decomposition; native missing-data handling | Linear in observation space; may underfit highly nonlinear morphology | **A** |
| Multi-view VAE / Biolord | Nonlinear, generative, supports counterfactual generation, disentanglement | More hyperparameters; identifiability harder | **A−** |

**Choice: MOFA+ as the primary integration backbone, with Biolord as a secondary nonlinear refinement.** (i) MOFA+ provides explicit per-modality variance decomposition — *exactly* what we need to defend U2OS-to-brain transfer; if a factor explains > 70% of variance in JUMP-CP but ~0% in iPSC-MEA, it is flagged morphology-only and downweighted for cognition relevance. (ii) MOFA+ natively handles missing modalities (most compounds will have L1000 but not iPSC-MEA, and chemCPA-imputed signatures should be a clearly-labeled view). (iii) MOFA+ scales to ~200 K compound × ~50 latent factor on RTX 5070 + 32 GB RAM. (iv) Biolord is layered on top as a nonlinear corrector for residual variance in the highest-importance factors, giving counterfactual generation for "what would compound X's signature look like in cell line Y" queries.

**Implementation:** `mofapy2` (Python) with views = {L1000_zscore, CP_CellProfiler, CP_DeepProfiler, CP_DINO, MEA_features, snRNA_pseudobulk, chemCPA_latent}; groups = {neural_lineage, non_neural_lineage, imputed}; K = 30–50 with ARD sparsity priors. Output: per-compound K-dim factor vector + per-view variance attribution.

**Similarity metrics:** L1000 → WTCS τ. JUMP-CP → cosine on sphered DeepProfiler 672-dim embeddings. MEA → cosine on log-transformed feature vector. snRNA → cosine on scVI latent. Joint → cosine on MOFA+ factor vector (the primary signal entering V8 likelihood).

---

### C. Joint Hierarchical Bayesian Model with Conditional Dependencies

#### C.1 DAG

```
                   cognition_priors (V6.B: AHBA + OT-Genetics + cellxgene)
                            │
                            ▼
                    target_relevance θ_T  ←── GWAS-anchored ~210 targets
                       ╱          ╲
                      ▼            ▼
              binding θ_B          (covariates: PBPK, dose, exposure)
                │     ╲            ╱
                ▼      ▼          ▼
            phenotype θ_P  ───────────► effect_size θ_E (Hedges' g)
            (V8 NEW)                    (V7 PRISMA-anchored)
```

#### C.2 Decomposition

p(θ_T, θ_B, θ_P, θ_E | D) ∝ p(θ_T | D_V6B) · p(θ_B | θ_T, D_V6A) · p(θ_P | θ_B, θ_T, D_V8) · p(θ_E | θ_B, θ_T, θ_P, PBPK, m, D_V7)

**Phenotype likelihood.** For compound c with MOFA+ factor φ_c ∈ ℝ^K, V6.B target relevance r_c, V6.A binding b_c:

φ_c | θ_B, θ_T ∼ 𝒩( A·b_c + B·r_c + C·(b_c ⊗ r_c) , Σ_φ(τ_chemCPA, τ_cellline) )

- A captures the binding-driven phenotypic axis (AChE-I binding → AChE-I-class L1000/CP signature).
- B captures target-relevance-driven cellular state.
- C is the **interaction term** — cross-product capturing conditional dependence (AChE-I + cognition-relevance ⇒ a *specific* AChE-I-cognition signature).
- Σ_φ inflates with τ_chemCPA for imputed signatures and τ_cellline for non-neural cell lines.

**Effect-size likelihood (V7 extended).**

g_c | θ_B, θ_T, θ_P, m, PBPK ∼ 𝒩( μ_mech(θ_T) + β_B·b_c + β_P·φ_c + β_PK·AUC_c − Σ_k γ_k·m_k , σ_g² )

The β_P·φ_c term is what V8 contributes to V7 — a phenotypic adjustment to the mechanism-class g prior.

#### C.3 Direct conditional priors vs. copulas

**Direct conditional priors (chosen):** A, B, C estimated jointly; conditional dependence captured in the mean function. Identifiable because we have ≥ 100 K compounds with signatures vs. ~210 targets × ~30 mechanism classes — degrees of freedom favor identifiability.

**Copula alternative (rejected for V8.1, queued for V8.2):** Gaussian copula on (θ_B, θ_T, θ_P) with marginals from V6.A/V6.B/V8 individually, plus an empirically-estimated coupling correlation R. Cleaner separation of concerns but harder to debug; rejected for first publication.

#### C.4 PyMC implementation skeleton

```python
import pymc as pm
import pytensor.tensor as pt

with pm.Model(coords={
    "compound": compound_ids,
    "target":   target_ids,
    "factor":   range(K),
    "mech_cls": mechanism_class_ids,
}) as v8_model:

    # V6.B target relevance (frozen mean/sd from V6.B posterior)
    mu_T   = pm.Data("mu_T_post",   v6b_mean)
    sig_T  = pm.Data("sig_T_post",  v6b_sd)
    theta_T = pm.Normal("theta_T", mu=mu_T, sigma=sig_T, dims=("compound","target"))

    # V6.A binding
    mu_B  = pm.Data("mu_B_post",  v6a_mean)
    sig_B = pm.Data("sig_B_post", v6a_sd)
    theta_B = pm.Normal("theta_B", mu=mu_B, sigma=sig_B, dims=("compound","target"))

    # V8 phenotype likelihood
    A = pm.Normal("A", 0, 1.0, dims=("factor","target"))
    B = pm.Normal("B", 0, 1.0, dims=("factor","target"))
    C = pm.Normal("C", 0, 0.5, dims=("factor","target"))  # interaction (tighter prior)
    mean_phi = (pt.dot(theta_B, A.T)
                + pt.dot(theta_T, B.T)
                + pt.dot(theta_B * theta_T, C.T))

    tau_chemCPA = pm.Data("tau_chemCPA", chemcpa_uncertainty)
    tau_cellline = pm.Data("tau_cellline", cellline_uncertainty)
    sigma_phi = pm.HalfNormal("sigma_phi", 1.0) * tau_chemCPA * tau_cellline

    phi_obs = pm.Normal("phi", mu=mean_phi, sigma=sigma_phi,
                        observed=mofa_factors, dims=("compound","factor"))

    # V7 effect size extended
    mu_mech = pm.Data("mu_mech", v7_mechclass_prior_mean)
    beta_B  = pm.Normal("beta_B", 0.0, 0.3)
    beta_P  = pm.Normal("beta_P", 0.0, 0.3)
    beta_PK = pm.Normal("beta_PK", 0.0, 0.3)
    gamma   = pm.Normal("gamma_moderators", 0.0, 0.3, shape=(5,))

    g_mean = (mu_mech
              + beta_B * theta_B.sum(axis=-1)
              + beta_P * mean_phi.sum(axis=-1)
              - pt.dot(moderator_matrix, gamma)
              + beta_PK * pbpk_auc)

    sigma_g = pm.HalfNormal("sigma_g", 0.15)
    g_obs = pm.Normal("g", mu=g_mean, sigma=sigma_g,
                      observed=hedges_g, dims="compound_known_g")

    trace = pm.sample(2000, tune=2000, chains=4,
                      nuts_sampler="numpyro", target_accept=0.95)
```

Runtime estimate on RTX 5070 12 GB + 32 GB RAM + WSL2: ~8–16 h with numpyro JAX backend, ~50 K compounds × 30 factors, partial pooling on mechanism class. Memory-dominant; use `pm.Data` for large fixed matrices and a minibatched ADVI sanity-pass first.

#### C.5 Identifiability analysis
- A, B: identifiable; cross-compound variability in (θ_B, θ_T) pins down rows.
- C (interaction): **weakly identifiable** — tighter prior σ = 0.5; sensitivity sweep on prior.
- β_B, β_P (translation weights): identifiable only via ~50–100 compounds with known g — binding constraint.
- τ_chemCPA, τ_cellline: identifiable from out-of-fold validation.

#### C.6 Sensitivity sweep — λ_phen
Posterior at λ_phen ∈ {0, 0.25, 0.5, 1.0, 2.0}; report top-50 shortlist overlap. **Acceptance: ≥ 60% overlap with λ ∈ [0.5, 1.0] (stable inference), ≤ 90% (V8 is doing something).**

#### C.7 Worked example — "strong phenotype, weak V6.A"
Novel scaffold (Tanimoto < 0.3 to any ChEMBL ligand): V6.A binding posterior centered at pchembl 5.5 ± 0.6 (weak); V6.B target relevance uninformative; V8 chemCPA-imputed L1000 signature shows WTCS τ = 0.83 to donepezil and 0.71 to galantamine; JUMP-CP DeepProfiler cosine 0.68 to AChE-I cluster centroid; MOFA+ factor 7 (AChE-I) loads heavily.

Joint posterior: θ_B pulls upward toward AChE-I via the C interaction; θ_T pulls toward AChE-related relevance; β_P · φ_c adds positive shift to V7 g mean. **Facet tag: `weak_binder.strong_phenotype.AChE-I-like`.** Without V8, this compound is ranked low; with V8, it ranks in the top 5%.

---

### D. Disagreement-as-Signal Facet Taxonomy (Four-Axis)

#### D.1 Axes

| Axis | Source | Quantity | Uncertainty |
|---|---|---|---|
| 1. Binding (θ_B) | V6.A ensemble | pchembl posterior mean | Posterior SD + ensemble disagreement |
| 2. Target relevance (θ_T) | V6.B AHBA prior | Cognition relevance score | Posterior SD + Roberts SMD ceiling gate |
| 3. Effect size (θ_E) | V7 hierarchical | Hedges' g posterior | 95% CrI width |
| 4. Phenotype (θ_P) | V8 MOFA+ | Phenotypic match score (cosine to mechanism-class centroid) | MOFA+ factor variance + chemCPA noise |

JSD(P_i, P_j) = ½ KL(P_i‖M) + ½ KL(P_j‖M), with M = (P_i + P_j)/2; six pairwise JSDs per compound → 4-element "axis isolation" vector.

#### D.2 Facet-tag taxonomy

| Facet tag | Pattern | Interpretation | Priority |
|---|---|---|---|
| `agreement.all_high` | All 4 axes high | Canonical positive (AChE-I-class for AD) | A+ |
| `strong_binder.weak_phenotype` | High θ_B, low θ_P | Off-pathway target, wrong functional state, or compensatory | B |
| `weak_binder.strong_phenotype` | Low θ_B, high θ_P | Polypharmacology / off-target / novel MOA | **A+** |
| `strong_BTP.weak_E` | High B, T, P; low g | Mechanism validated but PK/PD-limited | A |
| `phenotype_only` | Only θ_P high | Pure phenotypic discovery — orphan signal | **A+** |
| `binding_only` | Only θ_B high | Computational binding artifact / non-functional binder | C |
| `relevance_only` | Only θ_T high | Genetic relevance without compound | not actionable |
| `target-true.phenotype-failed` | High B+T, inert P | **Encenicline pattern** — target hypothesis without functional consequence | C (V8 saves us) |

#### D.3 Wet-shortlist priority score
prio_c = w_1·E[θ_E^(c)] + w_2·log(1 + Σ_k JSD_k^(c)) − w_3·CrI-width(θ_E^(c)) + w_4·novelty(c)

Disagreement JSDs *boost* priority, consistent with the V6 facet-tag philosophy: disagreement is information.

---

### E. Gate 1 — Mechanism-Discovery Validation (PRIMARY)

#### E.1 Reference truth — ~30 PRISMA mechanism classes

**Cholinergic:** AChE-I, M1 PAM, M4 PAM, α7 nAChR agonist/PAM, α4β2 partial agonist.
**Glutamatergic:** NMDA antagonist (uncompetitive / partial-channel), mGluR2/3, mGluR5 NAM/PAM, AMPA potentiator.
**Monoaminergic:** DAT inhibitor, NET inhibitor, MAO-A/B, COMT, α2A, 5-HT6 antagonist, 5-HT4 agonist, 5-HT1A partial agonist, 5-HT7 antagonist.
**GABAergic:** α5 inverse agonist.
**Other:** Neuropeptide (oxytocin, vasopressin), σ1 agonist, H3 antagonist, CB1 modulator, PDE9 inhibitor, PDE4D inhibitor, HCN/KCNQ modulator, mitochondrial enhancer, growth-factor mimetic (BDNF), erythropoietin-class, hormonal (estradiol, pregnenolone).

#### E.2 Clustering protocol
1. MOFA+ joint factor vector for all compounds.
2. SNN graph (k = 15, Jaccard).
3. **Leiden** (primary; γ ∈ {0.4, 0.6, 0.8, 1.0, 1.2}, pick best AMI on held-out half).
4. **HDBSCAN** (secondary; min_cluster_size ∈ {15, 25, 50}).
5. Compute AMI, ARI, V-measure, Fowlkes-Mallows.

#### E.3 Pre-registered thresholds

| Band | AMI | ARI | Action |
|---|---|---|---|
| **PASS** | ≥ 0.5 | ≥ 0.4 | Enter joint posterior at λ_phen = 1.0 |
| **DEGRADE** | [0.3, 0.5) | [0.25, 0.4) | Enter at λ_phen = 0.5; flag in report |
| **FAIL** | < 0.3 | < 0.25 | Do not enter joint posterior; publish negative |

#### E.4 Stratified evaluation
- **Per mechanism class.** Expected: stimulants (DAT/NET) and AChE-I cluster cleanly; novel classes (σ1, mitochondrial, growth-factor mimetic) struggle.
- **By data modality.** L1000-only vs JUMP-CP-only vs iPSC-MEA-only vs MOFA+ joint. *If joint AMI does not exceed best single modality by ≥ 0.05, multi-modal architecture is not earning its complexity.*
- **chemCPA-imputed vs observed.** If imputed-only AMI ≪ observed-only, inflate τ_chemCPA.

#### E.5 Pre-registration
OSF.io pre-registration **before unblinding** mechanism-class labels. Locked: K, Leiden γ space, AMI/ARI bands, mechanism list, stratification scheme. This is the integrity firewall against the "clusters defined by mechanism class" circularity objection — labels come from chemistry + literature, not phenotype.

---

### F. Gate 2 — Held-Out *g* Prediction (SECONDARY)

#### F.1 Training set
~50–100 compounds with usable healthy-adult Hedges' *g* from V7 PRISMA priors (Roberts 2020 + MetaPsy + Cochrane).

#### F.2 Model
Features: MOFA+ K-dim factor + chemCPA-imputation flag + cell-line-mix indicator. Model: **Gaussian process regression** with Matérn-5/2 kernel (transparent uncertainty at low data volume). Bayesian NN (numpyro + flax, ~50 K params, HMC) as robustness check.

#### F.3 Split & metrics
- 80/20 stratified by mechanism class; **leave-one-mechanism-class-out** (LOMCO) as the harder generalization test.
- Metrics: MAE on held-out g, calibration plot (predicted 90% CrI coverage), Spearman ρ rank order.

#### F.4 Thresholds

| Band | MAE | 90% CrI coverage | Action |
|---|---|---|---|
| **PASS** | < 0.20 | ≥ 85% | V8 → V7 effect-size translation validated |
| **DEGRADE** | [0.20, 0.35] | [70%, 85%) | Use V8 → V7 but flag β_P uncertainty |
| **FAIL** | > 0.35 | < 70% | β_P = 0; V8 contributes only to ranking |

Honest framing: Roberts 2020 ceiling is g = 0.43 (MPH delayed recall); overall mean SMDs are 0.12 (modafinil) – 0.21 (MPH). With this small dynamic range and ~50–100 anchors, MAE < 0.20 is non-trivial — we may legitimately land in DEGRADE.

---

### G. Gate 3 — Connectivity vs. Nootropics Sanity Check

#### G.1 Reference set (9+1)
**Positive:** donepezil, memantine, modafinil, methylphenidate, atomoxetine, varenicline, galantamine, rivastigmine, donepezil+memantine combo.
**Negative control:** encenicline.

#### G.2 Expected nearest-neighbor structure

| Pair | Expectation | Mechanism |
|---|---|---|
| Donepezil ↔ Galantamine | Tight (top-10) | Both AChE-I; galantamine also α7 PAM |
| Donepezil ↔ Rivastigmine | Tight (top-10) | Both AChE-I |
| MPH ↔ Atomoxetine | Tight (top-20) | DAT+NET vs NET-selective |
| MPH ↔ Modafinil | Moderate (top-50) | DAT vs DAT-low-affinity / wake-promoting |
| Modafinil ↔ Armodafinil/Caffeine | Tight (top-20) | Wake-promoting class |
| Memantine ↔ Ketamine | Tight (top-20) | NMDA antagonists |
| Memantine ↔ Amantadine | Tight (top-20) | NMDA + dopaminergic |
| Varenicline ↔ Galantamine (α7) | Moderate (top-100) | α7 nAChR ligands |
| Donepezil+Memantine ↔ Donepezil | Tight | Combination → donepezil-dominant + memantine offset |
| **Encenicline ↔ active cluster** | **Distant (> 500 NN rank)** | **Negative control — Phase 3 inert** |

#### G.3 Pass threshold
≥ 7/9 expected positive-NN relationships hold *and* encenicline is ranked > 500 from active-cluster centroid. Single failure (e.g., donepezil ≁ galantamine) → L1000+JUMP-CP integration is broken.

---

### H. chemCPA Generative Imputation (detailed)

#### H.1 Training data
- **Bulk pretraining:** LINCS L1000 Level-5 MODZ (GSE92742 + GSE70138 + clue.io beta), filtered to ≥ 3 replicates per (compound, cell line, dose) and τ-quality > 0.5.
- **Single-cell fine-tuning:** sci-Plex3 (Srivatsan 2020) + any LINCS single-cell extension; JUMP-CP DeepProfiler embeddings via a second decoder head (architecture surgery).
- **Held-out:** (a) ≥ 100 compounds with Tanimoto < 0.5 to training; (b) entire glutamatergic class as LOMCO.

#### H.2 Architecture & training
- Molecule encoder: RDKit Morgan-FP-1024 → 256-MLP → 128-d perturbation latent (RDKit-pretrained per Hetzel 2022 best benchmark).
- Cell-line encoder: 1-hot → 64-d embedding; dose, time continuous.
- Decoder: 2-layer MLP → 977-gene reconstruction (LINCS landmark) or 12,328 (BING inferred).
- Adversarial discriminator on basal latent → perturbation identity (gradient reversal); λ_adv = 1.0.
- AdamW lr=1e-4, weight decay 0.01, 200 epochs, early stopping on held-out R²_DEGs.

#### H.3 Validation targets
- Random-split: R² ≥ 0.7 / DEGs ≥ 0.5.
- Scaffold-split: R² ≥ 0.5 / DEGs ≥ 0.3.
- LOMCO: R² ≥ 0.3 / DEGs ≥ 0.15.
- Anchor: sci-Plex3 9-OOD at 10 µM, target R²(all) ≥ 0.65 / DEGs ≥ 0.40 (≤ Hetzel 2022 ceiling of 0.69/0.47 to account for broader training mixture).

#### H.4 Failure modes
- `chemCPA.low_neighbor` (max-Tanimoto-to-train < 0.3): τ_chemCPA × 3.
- `chemCPA.polypharm_risk` (≥ 3 high-confidence ChEMBL targets pchembl ≥ 7): τ_chemCPA × 2 + polypharm-specific MOFA+ factor.
- `chemCPA.outside_chembl` (Lipinski violations ≥ 2): τ_chemCPA × 5; heavily downweighted.

#### H.5 Alternative architectures

| Method | Strengths | Weaknesses | Verdict |
|---|---|---|---|
| **chemCPA (Hetzel 2022)** | OOD compound support; bulk-to-sc surgery; benchmarked on sci-Plex3 | Mono-target bias; weak on novel mechanism class | **Primary** |
| CPA (Lotfollahi 2023) | Cleaner combinatorial dose-time disentanglement | No OOD compound capability | Sensitivity-only |
| GEARS (Roohani 2024) | Knowledge-graph prior; multigene | Designed for genetic perturbations | Reject |
| scGen (Lotfollahi 2019) | Simple latent arithmetic; transparent | No native unseen-compound support | Cross-check |
| Biolord (Piran 2024) | Best published OOD R² (~0.76 ± 0.0005) vs. chemCPA-pre (0.51 ± 0.0062) on sci-Plex3 | Newer, less broadly benchmarked | **Secondary** |

---

### I. Worked Examples

#### I.1 Donepezil (positive control)
- V6.A binding: AChE pchembl ~8.5 ± 0.3.
- V6.B target relevance: AChE cognition relevance ~0.92.
- V7 effect size: Hedges' g posterior ~0.18 ± 0.10 (consistent with the literature meta-analysis).
- V8 phenotype: WTCS τ to galantamine = 0.81, rivastigmine = 0.78; MOFA+ factor 7 (AChE-I) loads strongly.
- **Facet tag:** `agreement.all_high`. **A+**.

#### I.2 Encenicline (negative control — the case for V8)
- V6.A: α7 nAChR pchembl ~7.8.
- V6.B: α7 cognition relevance moderate (~0.55).
- V7: Hedges' g posterior ~0 ± 0.15 — Phase 2 (Keefe 2015) d=0.257 (0.27 mg, p=0.034) did not replicate at 0.9 mg (d=0.093, p=0.255), and FORUM Pharmaceuticals' March 2016 topline announcement (per Alzforum) reported both Phase 3 schizophrenia trials missed co-primary endpoints.
- V8 phenotype: Expected WTCS to active AChE-I cluster < 0.2; JUMP-CP DeepProfiler cosine to AChE-I centroid < 0.3; MOFA+ joint factor vector does NOT load on AChE-I/cognition factors.
- **Facet tag:** `target-true.phenotype-failed`. **A+ — V8 saves us from a false positive even before V7 is integrated.**

#### I.3 Intepirdine
- V6.A: 5-HT6 pchembl ~8.0.
- V6.B: 5-HT6 cognition relevance moderate.
- V7: g ≈ 0 — MINDSET (Lang et al. *Alz Dem TRCI* 2021, n=1,315) reported ADAS-Cog adj mean diff −0.36 (p=0.225) and ADCS-ADL −0.09 (p=0.826).
- V8 phenotype: Predicted weak clustering with 5-HT modulators; not in the active-cognition-enhancer set.
- **Facet tag:** `target-true.phenotype-failed` (5-HT6 class).

#### I.4 Lithium (test case for `weak_binder.strong_phenotype`)
- V6.A: No clear cognition target with high pchembl (GSK3β pIC50 ~5–6).
- V6.B: Weak (Li⁺ not GWAS-anchored to cognition).
- V8 phenotype: Broad GSK3β downstream + Wnt-pathway signature; distinctive LINCS signature. JSD(θ_B, θ_P) high.
- **Facet tag:** `weak_binder.strong_phenotype.GSK3β-Wnt`. **A — illustrates V8 polypharmacology capture.**

#### I.5 Selegiline
- V6.A: MAO-B pchembl ~7.5.
- V6.B: MAO-B cognition relevance moderate.
- V7: g ~0.10–0.20 (modest).
- V8 phenotype: Clusters with monoaminergic class (rasagiline, deprenyl).
- **Facet tag:** `agreement.all_medium`. **A class-consistent.**

#### I.6 Novel scaffold A (chemCPA imputation works)
- V6.A pchembl ~7.0 on a novel cognition target via MAMMAL/MMAtt-DTA (Tanimoto < 0.4 to known ligands).
- V8: chemCPA-imputed signature clusters with M1 PAM compounds; flag `chemCPA.imputed.confidence_medium`.
- Joint posterior: phenotype evidence amplifies binding evidence → effect-size posterior shifts up. Wet-shortlist priority high.

#### I.7 Novel scaffold B (DeepProfiler captures what Tanimoto misses)
- Tanimoto < 0.3 to all known active cognition compounds; strong JUMP-CP morphology signature matching AChE-I cluster.
- Canonical Cell Painting use case (Bray 2016 et seq.). **V8's central pitch.**

#### I.8 Aripiprazole (polypharmacology)
- V6 ranks low (D2 partial agonist, sedation, metabolic liability).
- V8 phenotype broad — multiple MOFA+ factors load (D2, 5-HT2A, 5-HT1A); many JSDs.
- **Facet tag:** `polypharmacology.broad_phenotype`. **B for direct cognition use; illustrative.**

#### I.9 Pridopidine (σ1, complex)
- V6.A: σ1 binding moderate.
- V6.B: σ1 cognition relevance uncertain.
- V7: PROOF-HD (Reilmann et al. *Nat. Med.* 2025) TFC LS mean diff −0.18 (p=0.26), cUHDRS −0.11 (p=0.45) — g ≈ 0.
- V8 phenotype: σ1 agonists have distinctive but small signature; MOFA+ factor for σ1 may or may not separate.
- **Facet tag:** `target-true.phenotype-weak.σ1` (if σ1 factor emerges) or `phenotype_inconclusive`.

#### I.10 Caffeine (sanity-check positive)
- Should cluster with stimulants in V8.
- **Facet tag:** `agreement.stimulant_class`.

---

### J. Implementation Spec (~30%)

#### J.1 Data ingestion pipeline

```
data/
├── lincs_l1000/
│   ├── GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz   (~7-8 GB)
│   ├── GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz
│   └── compoundinfo_beta.txt, geneinfo_beta.txt, siginfo_beta.txt
├── jump_cp/
│   ├── cpg0016/deepprofiler_well_consensus.parquet (672-dim, ~20-40 GB)
│   ├── cpg0016/cellprofiler_well_consensus.parquet (~10-15 GB)
│   └── cpg0016/dinov2_well_consensus.parquet       (~10 GB)
├── ipsc_neuron/
│   ├── mea_features.parquet         (per-compound aggregated)
│   └── snrnaseq_pseudobulk.h5ad     (cellxgene-census brain organoid slice)
├── reference/
│   ├── psychencode_brainspan.h5ad
│   └── ahba_v6b_anchors.parquet     (V6.B already loaded)
└── chemcpa/
    ├── trained_model.pt              (~200 MB)
    ├── imputed_signatures.parquet    (~50K cmpd × 977 landmark)
    └── imputation_uncertainty.parquet
```

Total disk: ~55 GB local cache. DuckDB queries over parquet for shortlist filtering; cmapPy for gctx; scanpy + anndata for h5ad.

#### J.2 PyMC V7→V8 extension

```python
class V8PhenotypicJointModel(V7HierarchicalBayesModel):
    def __init__(self, mofa_factors, chemcpa_uncertainty,
                 cellline_uncertainty, **v7_kwargs):
        super().__init__(**v7_kwargs)
        self.mofa = mofa_factors
        self.tau_chemcpa = chemcpa_uncertainty
        self.tau_cellline = cellline_uncertainty
        self.K = mofa_factors.shape[1]

    def build_phenotype_node(self, model):
        with model:
            A = pm.Normal("A", 0, 1, shape=(self.K, self.n_targets))
            B = pm.Normal("B", 0, 1, shape=(self.K, self.n_targets))
            C = pm.Normal("C", 0, 0.5, shape=(self.K, self.n_targets))
            mu_phi = (pt.dot(self.theta_B, A.T)
                      + pt.dot(self.theta_T, B.T)
                      + pt.dot(self.theta_B * self.theta_T, C.T))
            sigma_phi = pm.HalfNormal("sigma_phi", 1.0)
            phi = pm.Normal("phi",
                            mu=mu_phi,
                            sigma=sigma_phi * self.tau_chemcpa * self.tau_cellline,
                            observed=self.mofa)
            return phi

    def extend_effect_size_node(self, model):
        with model:
            beta_P = pm.Normal("beta_P", 0.0, 0.3)
            self.g_mean = self.g_mean + beta_P * pt.sum(self.mofa, axis=-1)
```

#### J.3 Runtime budget on RTX 5070 12 GB + WSL2

| Stage | Time | Disk | RAM | GPU VRAM |
|---|---|---|---|---|
| L1000 ingest + WTCS index | 1–2 h | ~10 GB | 8 GB | — |
| JUMP-CP DeepProfiler/CP download | 2–4 h (network) | ~30 GB | 4 GB | — |
| iPSC MEA + snRNA-seq ingest | ~1 h | ~10 GB | 6 GB | — |
| chemCPA training (sci-Plex3 + LINCS subset) | 4–8 h | <2 GB | 12 GB | 9–11 GB |
| chemCPA inference (50K compounds) | ~1 h | <1 GB | 4 GB | 4 GB |
| MOFA+ fit (50K cmpd × 7 views) | 2–4 h | <1 GB | 24 GB (CPU) | — |
| PyMC NUTS V8 joint (numpyro) | 8–16 h | <2 GB | 12 GB | 8 GB |
| Gate 1 clustering + AMI | <1 h | — | 6 GB | — |
| Gate 2 GP held-out | <1 h | — | 4 GB | <2 GB |
| **Total wall-clock** | **~24–36 h** | **~55 GB** | **peak 24 GB** | **peak 11 GB** |

**WSL2 file-system note:** large parquet writes are noticeably slower on the WSL2/Windows interop boundary. Keep working data inside the Linux filesystem (`/home`), not `/mnt/c`. MOFA+ at 24 GB peak is close to the 32 GB ceiling; consider `mofapy2 gpu_mode=False, sparse_data=True` and minibatching if it OOMs.

#### J.4 Test coverage plan
- Unit tests per ingester (LINCS, JUMP-CP, MEA, snRNA, chemCPA).
- Integration tests: held-out 100-compound corpus with all four axes; assert facet tags match expected.
- Snapshot tests on canonical references (donepezil, MPH, encenicline) — if facet tags flip between commits, abort.
- Property test: λ_phen = 0 posterior must reduce to V7-only posterior (sanity guard).

---

### K. Publishability Analysis

#### K.1 Novelty claim
*First* hierarchical Bayesian repurposing pipeline integrating a **multi-modal phenotypic prior** (transcriptomic L1000 + morphological JUMP-CP + electrophysiological iPSC-MEA + generative-imputed chemCPA) with a **target-first axis** (Multi-Head DTI + AHBA-anchored cognition prior + PRISMA mechanism-class effect-size translation) using **explicit conditional-dependence structure** and **four-axis disagreement-as-signal** facet tagging.

Prior art:
- L1000 + morphology integration (Way 2022; Haghighi et al. 2022 *Nat. Methods*) — feature-level, no Bayesian fusion.
- Bayesian repurposing (DeepPurpose, MOLI, MOLTPiper) — target-first, no phenotypic prior.
- Phenotype-driven repurposing (Sirota 2011, Iorio 2018, Subramanian 2017) — signature-only, no target prior fusion.
- Multi-modal embeddings (MOFA+, Biolord) — methodology, not coupled to clinical effect-size translation.

The integration with explicit conditional dependence + effect-size translation is, to literature search, novel.

#### K.2 Target venue ratings

| Venue | Realistic? | Required Gate 1 AMI | Verdict |
|---|---|---|---|
| **Nature Methods** | Stretch | ≥ 0.6 | A+ stretch (best fit for methodology novelty) |
| **Nature Machine Intelligence** | Yes | ≥ 0.5 | **A — best realistic fit** |
| Cell Reports Methods | Yes | ≥ 0.5 | A |
| Cell Systems | Yes | ≥ 0.4 | A− |
| Bioinformatics | Almost certain | ≥ 0.3 | A− |
| Nature Communications | Stretch | ≥ 0.5 + wet validation | B+ stretch (no wet-lab is hard) |

#### K.3 Anticipated reviewer objections & pre-emption

| Objection | Pre-emption |
|---|---|
| "Phenotypic data is noisy" | Gate 1 AMI/ARI; OSF pre-registration; stratified AMI by mechanism class |
| "L1000 has batch effects" | Level-5 MODZ already batch-corrected; per-plate τ > 0.5 filter; report sensitivity |
| "JUMP-CP is U2OS — irrelevant to cognition" | **The elephant.** Cell-line random effect; per-modality variance attribution; downweight if Gate 1 AMI low on neural subset; honest framing |
| "chemCPA may hallucinate" | τ_chemCPA propagation; LOMCO benchmark; low-confidence flag |
| "Mechanism-class clustering is circular" | Labels external (chemistry + literature); LOMCO test with held-out class; pre-registration |
| "Multi-modal embedding has identifiability issues" | MOFA+ ARD priors + variance decomposition table; Biolord cross-check |
| "β_P overconfident" | g ceiling 0.43; σ_g prior < ceiling/2; calibration plot |

#### K.4 Pre-registration as integrity firewall
OSF.io project locks before unblinding:
- K = 30 primary; sweep {20, 30, 40, 50}.
- Leiden γ ∈ {0.4–1.2}.
- AMI cuts 0.5 / 0.3.
- ARI cuts 0.4 / 0.25.
- Mechanism class list (~30 classes).
- Reference nootropic set (9 + 1 encenicline negative).
- Gate 2 MAE cut 0.20; CI coverage cut 85%.

---

### L. Identifiability and Sensitivity Analysis

- **L.1** A, B identifiable; C (interaction) weakly identifiable → tight prior σ=0.5 + sensitivity sweep.
- **L.2** chemCPA imputation uncertainty propagation via τ ∈ {0.5, 1, 2, 4}.
- **L.3** Embedding sensitivity: MOFA+ vs Biolord vs Deep CCA; report joint-vs-best-single AMI delta.
- **L.4** λ_phen sweep ∈ {0, 0.25, 0.5, 1, 2}; top-50 shortlist overlap.
- **L.5** Leiden γ ∈ {0.4–1.2}; HDBSCAN min_cluster_size ∈ {15, 25, 50}; robust band.
- **L.6** Cell-line transfer (U2OS → brain). Engineering response: (i) MOFA+ per-modality variance decomposition reports what % of compound variance is explained by JUMP-CP-only factors vs L1000 vs MEA/snRNA; (ii) Gate 1 AMI three ways — all-modality joint, neural-only (L1000 NEU/NPC + MEA + snRNA), non-neural (JUMP-CP U2OS); (iii) If non-neural-only AMI ≥ 0.5, U2OS Cell Painting transfers meaningfully — publish the finding; (iv) If non-neural < 0.3 alone but joint > 0.5, U2OS contributes as chemistry-anchored consistency check, not brain proxy — frame accordingly.

---

### M. Limitations and Honest Framing

1. **L1000 cell-line bias.** Mostly cancer lines; NEU/NPC L1000 thin (few cell lines × dozens of compounds). Mitigation: cell-line random effect; downweight non-neural for cognition lookups.

2. **JUMP-CP U2OS bias — the elephant.** U2OS, A549, HUVEC are not neurons. Defensible only if (a) factor-decomposition shows morphology captures target-class signal rather than cell context; (b) Gate 1 AMI on JUMP-CP-only embeddings exceeds 0.3 for AChE-I, NMDA antagonist, and DAT/NET classes (canonical Cell Painting wins). If those classes fail in U2OS, JUMP-CP enters at λ_phen = 0.25 instead of 1.0.

3. **iPSC neuron data is patchy and donor-heterogeneous.** Protocol-dependent (Brennand, Pașca, Studer, Mariani — different maturity, regional identity, donor pool). Donor random effect; protocol covariate. Note: Hartley & Brennand 2022 *PNAS* used **patch-clamp electrophysiology** (sEPSC amplitude predicted WCST with R = −0.93, p = 2.94×10⁻⁵) **not MEA** — the MEA literature comes from Frank 2017, Odawara 2016, Hyysalo 2017, and toxicology consortia; cite each precisely.

4. **chemCPA scaffold-extrapolation failure.** RDKit Morgan-FP has known weakness for novel scaffolds outside Lipinski-compliant chemistry; flag and downweight.

5. **Phenotypic signatures are MOA proxies, not mechanism proofs.** A compound that "looks AChE-I-class" in L1000 may not engage AChE; it may engage a downstream pathway producing similar transcription. Standard CMap caveat (Lamb 2006); acceptable for prioritization, not mechanism claims.

6. **Mechanism-class clustering circularity.** Handled by LOMCO + pre-registration.

7. **Healthy-adult cognition effects are small.** Roberts 2020 ceiling g ≈ 0.43; overall means 0.12–0.21. V8's resolution to distinguish active from inactive may be saturated by SNR; Gate 2 MAE < 0.20 may be hard.

8. **No wet-lab dollars; no prospective validation.** All claims in-silico against published anchors (Roberts 2020, MetaPsy, Cochrane). Pipeline produces a *prioritized shortlist*, not a validated discovery.

9. **Species mismatch.** Mouse-vs-human iPSC-MEA; species in LINCS minor but present. Translational caveat.

10. **Pre-registration mitigates but doesn't eliminate circular-reasoning risk.** A reviewer may claim mechanism-class labels are loosely derived from phenotype-adjacent literature. Honest framing in Discussion.

---

## Recommendations

**Stage 1 — Foundation build (weeks 1–3, ~80 h engineering).**
1. Pull LINCS L1000 GSE92742 + GSE70138 + clue.io beta; build cmapPy WTCS index.
2. Pull JUMP-CP cpg0016 DeepProfiler + CellProfiler consensus profiles only (skip raw images; cpg0016 raw images are ~115 TB per Chandrasekaran 2023 bioRxiv, well over local budget).
3. Aggregate iPSC-neuron MEA + snRNA from Frank 2017, Odawara 2016, Hyysalo 2017, and Synapse PsychENCODE.
4. **Decision gate:** verify LINCS+JUMP-CP coverage of V6 shortlist. If ≥ 60% has direct signatures, chemCPA is supplementary; if < 60%, chemCPA is load-bearing.

**Stage 2 — chemCPA training & validation (weeks 3–4, ~40 h GPU).**
1. Train chemCPA-RDKit-pretrained on LINCS + sci-Plex3.
2. Validate against 9-OOD sci-Plex3 — target R²(all) ≥ 0.65 / DEGs ≥ 0.40 (vs Hetzel 2022 ceiling 0.69/0.47 and Biolord-paper-cited cross-condition mean 0.51 ± 0.0062 for chemCPA-pre).
3. LOMCO on glutamatergic class — R² ≥ 0.30 / DEGs ≥ 0.15.
4. **Decision gate:** if LOMCO R² < 0.20, downweight imputed evidence 3× and flag publication scope.

**Stage 3 — MOFA+ joint embedding (weeks 4–5, ~20 h CPU).**
1. Fit MOFA+ K=30 across 7 views.
2. Variance attribution table — per-factor per-view contribution.
3. Sphering and batch correction sanity (donepezil cluster check on raw MOFA+ vectors).

**Stage 4 — Gate 1 pre-registration & evaluation (weeks 5–6).**
1. Lock predictions on OSF.
2. Compute AMI/ARI; stratified table.
3. **Decision gate:** AMI ≥ 0.5 → PASS (λ_phen=1.0). [0.3, 0.5) → DEGRADE (λ_phen=0.5). < 0.3 → FAIL, negative-result paper.

**Stage 5 — Joint posterior PyMC (weeks 6–8, ~16 h GPU).**
1. Extend V7 PyMC with V8 phenotype node + extended g-mean.
2. Sample with numpyro, target_accept=0.95.
3. Sensitivity sweep λ_phen ∈ {0, 0.25, 0.5, 1, 2}.

**Stage 6 — Gate 2 + Gate 3 + paper draft (weeks 8–10).**
1. Gate 2 held-out g GP regression.
2. Gate 3 nootropic anchors NN check.
3. Worked examples Section I — 10 compounds.
4. Draft for Nature Machine Intelligence.

**Threshold-changing benchmarks:**
- Gate 1 AMI ≥ 0.65 → upgrade target to Nature Methods.
- Gate 1 AMI < 0.3 → degrade to Bioinformatics or write as negative-result paper for Cell Reports Methods.
- chemCPA LOMCO R² < 0.15 → drop "full ~50K coverage" claim; restrict to compounds with direct signatures (the LINCS ∩ JUMP-CP overlap of which only ~5–10 K are V6-relevant).

---

## Caveats

- LINCS clue.io beta releases iterate (`level3_beta_all_n3026460x12328` is current as of writing); verify QC at ingestion.
- cpg0016 was still under final QC at Chandrasekaran 2024 publication: the JUMP-CP README states "we don't recommend performing any analysis with the principal dataset [until] the full QC of the dataset is complete." Verify QC status before final results.
- The Hetzel 2022 chemCPA OOD numbers (R² ≈ 0.69 all genes / 0.47 DEGs at 10 µM) are at the highest dose in sci-Plex3 only; the externally-cited cross-condition mean R² = 0.51 ± 0.0062 (Biolord paper, Piran 2024 *Nat. Biotechnol.*) is more conservative.
- The Roberts 2020 ceiling g ≈ 0.43 is the maximum *significant* subdomain effect (MPH delayed recall); overall mean SMDs are 0.12 (modafinil) and 0.21 (MPH) — Gate 2 thresholds set against realistic ceiling, not maximum.
- The Encenicline Phase 3 failure is documented via Alzforum's reporting of FORUM Pharmaceuticals' March 2016 topline press release; the primary peer-reviewed Phase 3 record is Brannan et al. *Schizophr. Bull.* 45(Suppl_2):S141 (2019 conference abstract). FORUM Pharmaceuticals subsequently wound down; no full Phase 3 manuscript was published to a PubMed-indexed journal.
- The Cell Painting Gallery total is **688 TB as of May 2024** (Weisbart et al.), of which cpg0016-jump is approximately 115 TB of images + numerical data (Chandrasekaran 2023 bioRxiv estimate).
- Hartley & Brennand 2022 *PNAS* (cognition prediction from iPSC neurons) was patch-clamp not MEA; cite that paper as evidence for the *electrophysiology → cognition* link but cite Frank 2017 / Odawara 2016 / Hyysalo 2017 for the MEA-specific protocols used in V8.
- The "first multi-modal phenotypic prior integrated with target-first Bayesian repurposing" novelty claim must be re-verified at submission against Recursion, Insitro, and Genentech proprietary disclosures — what matters for novelty is peer-reviewed published precedent, but the field is moving fast.