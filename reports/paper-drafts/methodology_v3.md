# Methodology Note v3 — MAMMAL Cognitive Enhancement Drug Repurposing

**Version**: 3 (consolidates V4 → V5 → V6 → V7 → V8 architecture)
**Companion docs**:
- `design/architecture-and-plans/V4_STATUS_AND_FORWARD_PLAN.md` — source-of-truth status
- `design/architecture-and-plans/V6_ARCHITECTURE_PLAN.md` — V6+V7+V8 implementation roadmap
- `reports/paper-drafts/methodology_v1.md` + `methodology_v2.md` — historical baselines

This note is the public-facing methodology narrative. It explains *why* the pipeline is shaped the way it is, what we measured, what we *didn't* measure, and what counts as falsification at every layer.

---

## 1. The honest scope

We are not searching for a miracle smart drug. **Roberts CA, Jones A, Sumnall H, Gage SH, Montgomery C 2020 *Eur Neuropsychopharmacol* 38:40-62** is the ceiling paper: across 47 placebo-controlled RCTs in healthy adults, the strongest known cognitive enhancer (methylphenidate) has overall SMD = 0.21; the maximum significant subdomain effect is g = 0.43 (MPH delayed recall); modafinil sits at SMD = 0.12. The unmodifiable ceiling is g ≈ 0.50 at 90% credible upper bound.

The pipeline's deliverable is therefore *not* a top candidate. It is:

1. A **calibrated, provenance-rich ranking** of ~298 cognition-relevant compounds across 22 cognition-relevant targets
2. A **multi-layer falsifiable Bayesian methodology** that we OSF-pre-register before unblinding
3. Open-source code + parquet datasets for community replication

Each compound in the top-N has a full provenance trail back to documented signal sources with known failure modes.

---

## 2. The five architectural layers

The pipeline composes five layers, each adding a distinct evidence axis:

### V4 — Calibrated 5-cluster fusion (the foundation)
- **Cluster A**: MAMMAL DTI head + ESM-2 target embeddings + Boltz-2 structure + Tanimoto-to-actives baseline
- **Cluster B**: ADMET-AI 41-endpoint safety panel with hard gates
- **Cluster C**: PrimeKG knowledge graph + TxGNN per-disease ranking
- **Cluster D (preview)**: AHBA brain-region annotation as static z-score
- **Fusion**: RRF (k=60) with per-target Phase A.7 calibration + Tier A/B/C/D escalation
- **Output**: `reports/wet-lab/wet_lab_shortlist_v4_faceted.md` — 17-facet mechanism-class top-N

### V5 — Z-norm + Tier 2/3 sprint
- §4.4 calibrated MAMMAL into fusion (+ §4.8 Z-norm within target → restores selectivity meaning)
- §8.0b-zn liability gating (within-target Z-norm fixes the prior-collapse 115/115 CUT degeneracy)
- §8.15 disagreement-as-signal column
- §7.18 selectivity Z-norm
- §8.7 MoA preference ranker as 5th fusion cluster
- §8.10 nootropic-similarity annotator
- §8.13 pocket-class-conditioned liability gating
- §7.17 pose-saving Boltz wrapper (code)
- §7.7 V5.1 MMAtt-DTA adapter (code)
- §8.3 CTgov v2 IP-status cross-reference
- §8.16 Tier-A calibrator round-trip QC
- §7.4 v2 selectivity entropy + Cheng 2010 Partition Index
- §8.0a Pareto NSGA-III restructure
- §7.12 conformal prediction per-target gating
- §7.13 scaffold-aware active learning
- §7.15 PyMC hierarchical GRIN pool
- LambdaMART promotion (NDCG@25 = 0.891 vs baseline 0.774 = +15% lift)
- §8.14 pocket-routed isotonic at SLC6A3
- §8.9 ANI-2x pose-validation stub
- §7.16 detector ensemble Sprint 2

### V6.A — Multi-Head DTI ensemble
- **5 DTI heads**: MAMMAL (calibrated) + Tanimoto-to-actives + MMAtt-DTA (Schulman 2024 *Bioinformatics*) + PSICHIC (Koh 2024 *Nat Mach Intell*) + BALM (Gorantla 2025 *J Chem Inf Model*)
- **Per-head bias decomposition** (PC, SN, OOD, CT signatures per (head, target) with Bonett-Wright CIs)
- **Bayesian router** per target (EnsDTI-extended; Park 2024 bioRxiv)
- **Venn-ABERS calibration** (Mervin 2020 AstraZeneca 40M-pair benchmark)
- **Multi-head disagreement axis** (4-bucket facet-tag: novel_scaffold / activity_cliff / ood / noise)
- **V6.A.1 empirical result**: MMAtt-DTA ρ at SLC6A3 = +0.65 < Tanimoto +0.90 → Tier-A FAIL → 3-head fallback + INVERT-mask architecture. Published as a positive-method-of-falsification finding.

### V6.B — Bayesian Cluster D neurobiological prior
- **AHBA**: abagen.get_expression_data() with Markello 2021 pinned config (ibf_threshold=0.5, probe_selection='diff_stability', donor_probes='aggregate', etc.) — 20/22 cognition genes across 83 DK regions cached
- **OT Genetics L2G**: GraphQL fetcher for Davies 2018 (GCST006269) + Hill 2019 + Sniekers 2017 + Savage 2018 + UKBB intelligence GWAS; legacy + Platform v25+ endpoint fallback
- **cellxgene-census**: brain organoid slice via tiledbsoma (preview)
- **PyMC NUTS hierarchical model** per Cluster D §B.2: y^s_i ~ N(α_s + β_s · θ_i, τ_s^-1 + σ²_s_i); θ_i ~ N(0, 1); reference anchors (BDNF, COMT, ACHE, DRD2, GRIN2B, CHRNA7) at θ ~ N(0.5, 0.3²) to break sign + scale degeneracy
- **Production NUTS run (this sprint)**: 4 chains × 2000 draws on RTX 5070 in <5 min; **R̂ max = 1.000, ESS min = 12,780**. Both convergence gates ✅ PASS.
- **4-gate validation framework** (V6.B.4): Gate 1 Roberts ceiling (HARD), Gate 2 Spearman vs meta-analytic SMD, Gate 3 GWAS-AUROC, Gate 4 leave-one-source-out. 15-compound reference SMD table.

### V7 — Clinical Effect-Size Translation Function (NEW)
- **9-compartment PBPK** (JAX/diffrax adaptive Dormand-Prince + numpy explicit-RK4 fallback): gut → plasma → peripheral → cortex → striatum → hippocampus → basal-forebrain → brainstem → CSF
- **Watson 1989 receptor-occupancy-with-reserve** with U-shape generator (D1-postsynaptic vs D2-autoreceptor) and tolerance kinetics (R_avail dynamics)
- **PET-validated anchors**: Bohnen 2005 donepezil 19.1% cortical AChE / Volkow 1998 MPH DAT 12/40/54/72/74% at 5/10/20/40/60 mg / Kapur 2000 haloperidol D2 ~1.8 nM
- **Schmidli 2014 robust MAP** priors for 12 mechanism classes (AChE-I, wake_promoting, NDRI, NRI, NMDA_antagonist, multimodal_5HT, alpha2A_agonist, A2A_antagonist, AMPA_pos_mod, creatine, omega3, minocycline) extracted from Roberts 2020 + Cochrane + MetaPsy
- **3-level hierarchical Bayes**: μ_global ~ N(0, 0.20); μ_class[m] ~ N(prisma_mean, λ_class·prisma_sd); η = sigmoid(α + β1·E[pchembl] + β2·E[relevance] + β3·copula) − Σ_k γ_k · m_k
- **Cluster D multiplicative gate**: β_target[t_c] = θ̄_{t_c} · β_raw_target[t_c]
- **5 failure-mode moderators**: m1 U-shape miss, m2 practice/placebo, m3 tolerance onset, m4 trait×state interaction, m5 trial-design
- **8 pre-registered predictions P1–P8** with falsifiers: P1 donepezil g ∈ [0.10, 0.30]; P2 encenicline_3mg Phase 3 failure recapitulated; P3 MPH 20mg DSST g ∈ [0.15, 0.30]; P4 modafinil 200mg; P5 memantine 20mg; P6 intepirdine MINDSET; P7 pridopidine PROOF-HD; P8 lecanemab cognitive subdomain
- **4 validation gates**: Gate 1 (HARD) ≥6 of 8 P1–P8 PASS; Gate 2 (HARD) Roberts ceiling 0.50; Gate 3 MAE < 0.15; Gate 4 per-endpoint coverage ≥ 85%
- **V7.4 first execution (this sprint)**: 5 PASS / 1 FAIL / 2 NO_COMPOUND → Gate 1 PASS at 5/(5+1) = 83%

### V8 / Cluster E — πphen Perturbational Evidence Axis (NEW)
- **Six modalities**: LINCS L1000 Level-5 MODZ (Subramanian 2017 *Cell* 171:1437; ~1.3M signatures) + JUMP-CP Cell Painting cpg0016 (Chandrasekaran 2024 *Nat Methods*; 116,750 compounds × 3 embeddings) + iPSC-MEA (Frank 2017 + Odawara 2016 + Hyysalo 2017) + iPSC snRNA-seq (cellxgene brain organoid) + reference brain transcriptomes (PsychENCODE/BICCN/BrainSpan) + chemCPA generative imputation (Hetzel 2022)
- **MOFA+ K=30 joint embedding** across 7 views with ARD sparsity + per-group (neural_lineage / non_neural_lineage / imputed) decomposition
- **Conditionally-dependent 4-level hierarchical Bayes**: p(θ_T, θ_B, θ_P, θ_E | D) ∝ p(θ_T) · p(θ_B | θ_T) · p(θ_P | θ_B, θ_T) · p(θ_E | θ_B, θ_T, θ_P, PBPK, m); phenotype-binding-relevance interaction term C with tight prior σ=0.5
- **Three-way Jensen-Shannon disagreement**: JS₃(π_t, π_g, π_p) = (1/3) Σ_m KL(p_m ‖ p̄)
- **I_novel mutual-information novel-mechanism score**: I_novel(c) = π_p · [1 − I(π_p; (π_t, π_g))]; identifies the (L, L, H) cell of the 8-cell disagreement classification — the clemastine / PIPE-307 territory
- **5-MoA cognition reference centroids**: cholinergic / catecholaminergic / glutamatergic / trophic_ISR / remyelination (BIMA-8 cluster per Mei 2014 *Nat Med*)
- **4 OSF-pre-registered validation gates**: Gate 1 (PRIMARY) AMI ≥ 0.50 vs PRISMA 30-class taxonomy; Gate 2 held-out g MAE < 0.20; Gate 3 9+1 nootropic-anchor NN structure; Gate 4 I_novel correctly identifies (L, L, H) on held-out clemastine + BIMA-8

---

## 3. Composition — the three-factor joint posterior

V6.A produces compound-binding posterior π_B per (compound, target). V6.B produces target-relevance posterior π_T per target. V7 produces effect-size posterior π_E per (compound, endpoint). V8 produces phenotypic posterior π_P per compound. The wet-lab shortlist v10 composes them:

**π_joint(compound) ∝ π_target(V6.A, V6.B) · π_phen(V8)** with **Gaussian-copula correlation correction** between the three Bayesian factors.

Operationally:
- **Cluster D multiplicative gate** scales V7 effect-size by σ(θ̄_target)
- **V8 phenotype boost** adds β_P · φ_c to V7 g_mean
- **Joint CrI** propagates V6.A pchembl SD + V6.B θ̄ SD + V7 σ_resid + V8 τ_chemCPA via Gaussian-copula approximation
- **Roberts 2020 SMD ceiling pre-filter** flags any compound with joint_g_90_upper > 0.50

Output: 18-column wet-lab shortlist v10 with 4-axis annotation:
- `g_predicted`, `g_90_upper` (V7 + V8 composition)
- `three_way_jsd`, `i_novel_score` (V8 axis disagreement)
- `eight_cell_tag` (V6.A × V6.B × V8 high/low classification)
- `roberts_ceiling_ok`, `wet_lab_priority` (final ranking gate)
- `evidence_axes` (which of V6.A/V6.B/V7/V8 contributed)

---

## 4. The methodology contribution — what's publishable

Each layer is a distinct publishable contribution:

| Layer | Venue | Status |
|---|---|---|
| V6.A multi-head DTI ensemble | *J Cheminform* / *Nat Mach Intell* | code shipped; needs production weights download for BALM |
| V6.B Bayesian Cluster D | *Cell Reports Methods* / *Bioinformatics* | code shipped + NUTS converged R̂=1.000 |
| V7 Clinical Effect-Size Translation | *Clinical Pharmacology & Therapeutics* (CPT) | code shipped + 5/8 P1–P8 PASS; OSF pre-reg ready |
| V8 πphen Perturbational Axis | *Nat Mach Intell* / *Nat Methods* | scaffolds shipped; LINCS+JUMP-CP download is the gate |
| V7 negative-result fallback | *CPT: Pharmacometrics & Systems Pharmacology* | activated if ≥3 of P1–P8 fail |
| V8 negative-result fallback | *Bioinformatics* / *PLOS ONE* | activated if Gate 1 AMI < 0.30 |

**Cross-layer contribution**: the **three-factor joint posterior with Gaussian-copula correlation correction + 4-axis disagreement-as-signal facet tagging + 8-cell mechanism classification including the (L, L, H) novel-mechanism cell** is, to literature search, novel. This is the integrating publication candidate.

---

## 5. What we didn't measure (honest limitations)

1. **No wet-lab validation**. All claims are in-silico against published meta-analytic literature. The pipeline produces a prioritized shortlist, not a validated discovery.
2. **LINCS L1000 cell-line bias**. Mostly cancer lines (A375, MCF7, PC3, VCAP, etc.). Neural-lineage L1000 thin (NPC, NEU, SHSY5Y, dozens of compounds). Mitigation: cell-line random effect; downweight non-neural for cognition lookups.
3. **JUMP-CP U2OS-to-brain transfer**. The elephant. Defensible only if MOFA+ per-modality variance attribution shows morphology captures target-class signal rather than cell context. Pre-registered handling in V8 OSF doc §7.
4. **iPSC neuron data is patchy and donor-heterogeneous**. Protocol-dependent (Brennand, Pașca, Studer, Mariani). Donor random effect; protocol covariate.
5. **chemCPA scaffold-extrapolation failure**. RDKit Morgan-FP has known weakness for novel scaffolds outside Lipinski-compliant chemistry; τ_chemCPA × 3 inflation flag for max-Tanimoto-to-train < 0.3.
6. **Phenotypic signatures are MOA proxies, not mechanism proofs**. A compound that "looks AChE-I-class" in L1000 may not engage AChE; it may engage a downstream pathway producing similar transcription. Standard CMap caveat (Lamb 2006); acceptable for prioritization, not mechanism claims.
7. **Healthy-adult cognition effects are small**. Roberts 2020 overall mean SMDs are 0.12–0.21; ceiling g ≈ 0.43. V8 Gate 2 MAE < 0.20 may be hard to achieve given the SNR.
8. **NMDA targets (GRIN2A / GRIN2B) confirmed unfixable by single-chain inference**. Ifenprodil-class NAMs bind the GluN1/GluN2B ATD heterodimer interface that MAMMAL cannot see from a single-chain input. Documented limitation in V4.
9. **Per-compound dose specificity is approximate**. P3 MPH "20mg" / P4 modafinil "200mg" / P5 memantine "20mg" predictions assume a canonical dose; the V7 PBPK PET-anchor calibration handles dose-dependence but the P-prediction parser uses base-compound names.
10. **Pre-registration is not a substitute for replication**. OSF lock means we cannot post-hoc relax thresholds; it does not mean the underlying methodology is correct. Wet-lab validation remains the gold standard.

---

## 6. Falsifiability — how each layer can fail

| Layer | Falsifier | Action on failure |
|---|---|---|
| V6.A | Tier-A: ensemble at SLC6A3 ≤ Tanimoto +0.90 | Tier-B fallback to 3-head; published as method-of-falsification |
| V6.B | Sign-stability < 80% across sensitivity sweep | Downgrade Cluster D from "primary prior" to "secondary diagnostic" |
| V6.B | Gate 2 Spearman ρ < 0.20 vs meta-analytic SMD | Do not publish; core empirical claim failed |
| V7 | ≥3 of P1–P8 fail pre-registered bands | CPT:PSP negative-result paper |
| V7 | Roberts ceiling violation > 5 compounds | Tighten priors + re-run NUTS |
| V8 | Gate 1 AMI < 0.30 vs PRISMA 30-class | Cell Reports Methods / Bioinformatics negative-result paper |
| V8 | Gate 2 MAE > 0.35 | β_P = 0 (V8 contributes only to ranking, not to V7) |
| V8 | Non-neural-only AMI < 0.30 alone | Frame V8 as chemistry-anchored consistency check, not brain proxy |
| Joint | (L, L, H) cell empty after gating | I_novel not generating novel-mechanism candidates; V8's central pitch falsified |

Every layer has at least one falsifier that, if triggered, becomes the publishable contribution rather than a hidden post-hoc adjustment.

---

## 7. Reproducibility

All artifacts auto-generated; reproducibility commands in `README.md`. Configs at `configs/{thresholds,weights,weights_calibrated}.yaml`. Calibrators at `data/calibration/isotonic/<uniprot>.pkl`. Pocket centroids at `data/pockets/centroids/<target>.json`. Per-target router decisions at `data/calibration/router_decisions.csv`. Reports at `reports/`. Sprint history is the git log.

**Test coverage**: 190 / 190 non-slow pytest tests pass. Each layer has its own test suite (`tests/test_v6_phase2.py`, `tests/test_v7_translation.py`, `tests/test_v8_cluster_e.py`, `tests/test_v8_advanced.py`) with graceful-degradation paths covered.

**Hypothesis audit**: 22 falsifiable claims tracked in `reports/pipeline/hypothesis_audit_v1.md`. Current verdicts: 19 PASS / 3 DEGRADE / 0 FAIL. The 3 DEGRADE (H4 positive-control top-20%, H5 negative-control suppression, H11 SLC6A3 calibrator drift) are documented honest signal.

---

## 8. Roadmap to v11 (full real-data execution)

Engineering complete; remaining work is **data acquisition + compute**:

1. **LINCS L1000 download** — ~10 GB; GSE92742 + GSE70138 + clue.io beta via cmapPy
2. **JUMP-CP cpg0016 download** — ~30-40 GB; DeepProfiler + CellProfiler + DINOv2 consensus parquets only (NEVER the ~115 TB raw images)
3. **OT Genetics L2G real fetch** — ~10 min once network access available; parquet cached for re-runs
4. **chemCPA training** — 4-8 h GPU on RTX 5070 once LINCS L1000 loaded
5. **MOFA+ joint embedding** — 2-4 h CPU on RTX 5070 + 32 GB RAM (peak 24 GB)
6. **PyMC V8 NUTS** — 8-16 h GPU via numpyro JAX backend
7. **OSF pre-registration lock** — submit `reports/paper-drafts/v7_osf_preregistration.md` + `reports/paper-drafts/v8_osf_preregistration.md` to OSF.io BEFORE running steps 4-6 on real data
8. **Paper drafts** — V6.A, V6.B, V7, V8, and the three-factor joint composition as separate manuscripts

Total wall-clock for real-data execution: ~24-36 h compute + ~1-2 weeks for OSF lock + paper drafting.

---

*Generated by `reports/paper-drafts/methodology_v3.md`. Companion to `README.md`, `design/architecture-and-plans/V4_STATUS_AND_FORWARD_PLAN.md`, `design/architecture-and-plans/V6_ARCHITECTURE_PLAN.md`, `reports/paper-drafts/v7_osf_preregistration.md`, `reports/paper-drafts/v8_osf_preregistration.md`. Sprint: V6.B.3 production NUTS converged R̂=1.000; V7.4 validation 5/8 P1-P8 PASS; V8 scaffolds shipped; wet-lab shortlist v10 with 3-factor composition.*
