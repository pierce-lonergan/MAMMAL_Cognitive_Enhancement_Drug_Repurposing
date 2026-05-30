# V8 / Cluster E — OSF Pre-Registration: πphen Perturbational Evidence Axis

**Project**: MAMMAL Cognitive Enhancement Drug Repurposing (V8 / Cluster E layer)
**Lead**: Pierce Lonergan
**Target venue**: *Nature Machine Intelligence* (A realistic at Gate 1 AMI ≥ 0.5); stretch *Nature Methods* (A+ at AMI ≥ 0.6); fallback *Cell Reports Methods* / *Cell Systems* / *Bioinformatics*
**OSF lock date**: TBD (lock BEFORE unblinding the mechanism-class labels)
**Pre-registration template**: OSF.io + AsPredicted.org
**Companion design doc**: `design/architecture-and-plans/V4_STATUS_AND_FORWARD_PLAN.md` §13.Z
**Research source docs**:
- `research/4-tier/Perturbational Evidence Axis.md` (V8 spec)
- `research/4-tier/Technical Feasibility Deep-Dive Adding a Phenotypic.md` (V8 companion)

---

## 1. Hypothesis

V8 introduces πphen, a **target-agnostic phenotypic evidence axis** parallel to V6.A (target-binding) and V6.B (target-relevance). The joint posterior over predicted cognition Hedges' *g* becomes:

**π_joint(compound) ∝ π_target(V6.A, V6.B) · π_phen(V8)** with Gaussian-copula correlation correction.

**Primary hypothesis** (Gate 1, PRIMARY): the MOFA+ K=30 joint factor embedding across 7 multi-modal views (LINCS L1000 + JUMP-CP Cell Painting CellProfiler + DeepProfiler + DINOv2 + iPSC-MEA + snRNA + chemCPA-imputed) recovers a PRISMA-anchored ~30-class mechanism taxonomy with **AMI ≥ 0.5 and ARI ≥ 0.4**.

**Secondary hypothesis** (Gate 2): held-out Hedges' *g* prediction via GP regression on MOFA+ factors achieves **MAE < 0.20 with 90% CrI coverage ≥ 85%** on a ~50–100-compound held-out anchor set.

**Tertiary hypothesis** (Gate 3): the 9+1 nootropic-anchor nearest-neighbor structure holds — ≥7 of 9 expected positive pairings (donepezil↔galantamine, MPH↔atomoxetine, etc.) AND encenicline ranks > 500 from the active cognition-enhancer cluster centroid.

**Quaternary hypothesis** (Gate 4): the I_novel mutual-information score correctly identifies the (L, L, H) "novel-mechanism" 8-cell on the held-out clemastine + BIMA-8 cluster — ≥6 of 8 anchor compounds in the top-5% I_novel rank.

---

## 2. Pre-registered model specification

### 2.1 Six-modality view stack (LOCKED)

| View | Source | Dimensionality | Citation |
|---|---|---|---|
| L1000_zscore | LINCS GSE92742 + GSE70138 + clue.io beta | 977 landmark (or 12,328 BING) | Subramanian 2017 *Cell* 171:1437 |
| CP_CellProfiler | JUMP-CP cpg0016 | ~700 (post feature-select) | Chandrasekaran 2024 *Nat Methods* 21:1114 |
| CP_DeepProfiler | JUMP-CP cpg0016 | 672 | Moshkov 2024 *Nat Commun* 15:1594 |
| CP_DINO | JUMP-CP cpg0016 | 384 | Sypetkowski 2024 |
| MEA_features | Frank 2017 + Odawara 2016 + Hyysalo 2017 + PsychENCODE | 25 | Frank 2017 *Toxicol Sci* 160:121 |
| snRNA_pseudobulk | cellxgene-census brain organoid slice | 1000 (scVI latent) | Lopez 2018 *Nat Methods* 15:1053 |
| chemCPA_latent | RDKit-Morgan-FP pretrained chemCPA | 128 | Hetzel 2022 NeurIPS |

### 2.2 MOFA+ K=30 joint embedding (LOCKED)

**Sampling sweep over K**: {20, 30, 40, 50}. **Primary**: K=30. Reported in the paper: per-factor per-view variance attribution table.

**Groups** for ARD sparsity: {neural_lineage, non_neural_lineage, imputed}.

**Implementation**: `mofapy2` (Python) with `ard_per_factor=True` + `spikeslab_weights=True`. Fallback: numpy SVD on z-scored concatenated views.

### 2.3 chemCPA generative imputation (LOCKED)

**Architecture** (Hetzel 2022): RDKit-Morgan-FP-MLP molecule encoder (1024→256→128-d) + cell-line embedding (64-d) + dose scaler + 977-gene decoder + adversarial discriminator.

**Training data**: LINCS L1000 Level-5 MODZ + sci-Plex3 + JUMP-CP DeepProfiler architecture-surgery head.

**Held-out validation** (per V8 plan §H.3):
- Random 80/20 — R² ≥ 0.70 / DEGs ≥ 0.50
- Scaffold-split — R² ≥ 0.50 / DEGs ≥ 0.30
- LOMCO (Leave-One-Mechanism-Class-Out) — R² ≥ 0.30 / DEGs ≥ 0.15
- sci-Plex3 9-OOD anchor — R²(all) ≥ 0.50, R²(DEGs) ≥ 0.30 (vs Hetzel 2022 ceiling 0.69/0.47; Piran 2024 cross-condition mean 0.51 ± 0.0062)

**Uncertainty inflation**: τ_chemCPA × 3 for max-Tanimoto-to-train < 0.3 (flag = `chemCPA.imputed.low_confidence`).

### 2.4 Joint posterior PyMC NUTS specification (LOCKED)

Per V8 plan §C:

```
p(θ_T, θ_B, θ_P, θ_E | D)
  ∝ p(θ_T | D_V6B)
    · p(θ_B | θ_T, D_V6A)
    · p(θ_P | θ_B, θ_T, D_V8)
    · p(θ_E | θ_B, θ_T, θ_P, PBPK, m, D_V7)
```

Phenotype likelihood:
```
φ_c | θ_B, θ_T ∼ N(A·b_c + B·r_c + C·(b_c ⊗ r_c),
                   Σ_φ(τ_chemCPA, τ_cellline))
```

**Hyperprior choices** (LOCKED):

| Hyperprior | Distribution |
|---|---|
| A (binding axis) | Normal(0, 1.0) shape=(K, n_targets) |
| B (relevance axis) | Normal(0, 1.0) shape=(K, n_targets) |
| C (interaction) | Normal(0, 0.5) shape=(K, n_targets) — TIGHTER prior; weakly identifiable |
| σ_φ baseline | HalfNormal(1.0) |
| β_P (effect-size translation weight) | Normal(0, 0.3) |

**Sampling**: 4 chains × 2000 tune × 2000 draws, target_accept=0.95, numpyro JAX backend.

### 2.5 8-cell disagreement classification (LOCKED)

3-bit (target, genetic, phenotype) → 8 cells:

| Bits | Tag | Interpretation |
|---|---|---|
| (1, 1, 1) | agreement.all_high | Canonical positive (donepezil, MPH) |
| (1, 1, 0) | target_true.phenotype_failed | Encenicline / intepirdine / pridopidine |
| (1, 0, 1) | target.phenotype | Binding + functional, no genetics |
| (1, 0, 0) | target_only | Binding artifact / off-pathway |
| (0, 1, 1) | genetic.phenotype | Genetic + functional, no good binder |
| (0, 1, 0) | genetic_only | GWAS but no actionable binder |
| **(0, 0, 1)** | **phenotype_only.novel_mechanism** | **Clemastine territory** ★ |
| (0, 0, 0) | no_evidence | Nothing |

### 2.6 I_novel novel-mechanism score (LOCKED)

```
I_novel(compound) = π_p · [1 − I(π_p ; (π_t, π_g))]
                    / τ_chemCPA
```

where I(·;·) is the mutual information between V8 phenotype and the joint (target, genetic) axes (approximated via Pearson correlation magnitude). Compounds with high I_novel are the V8 publishable novel-mechanism candidates.

---

## 3. Reference mechanism class taxonomy (LOCKED)

**~30 mechanism classes** per V8 plan §E.1:

**Cholinergic** (5): AChE-I, M1 PAM, M4 PAM, α7 nAChR agonist/PAM, α4β2 partial agonist.

**Glutamatergic** (4): NMDA antagonist (uncompetitive / partial-channel), mGluR2/3, mGluR5 NAM/PAM, AMPA potentiator.

**Monoaminergic** (9): DAT inhibitor, NET inhibitor, MAO-A/B, COMT, α2A, 5-HT6 antagonist, 5-HT4 agonist, 5-HT1A partial agonist, 5-HT7 antagonist.

**GABAergic** (1): α5 inverse agonist.

**Other** (~11): Neuropeptide (oxytocin, vasopressin), σ1 agonist, H3 antagonist, CB1 modulator, PDE9 inhibitor, PDE4D inhibitor, HCN/KCNQ modulator, mitochondrial enhancer, growth-factor mimetic (BDNF), erythropoietin-class, hormonal (estradiol, pregnenolone).

**Total**: ~30. Frozen at OSF lock time.

---

## 4. 9+1 nootropic-anchor reference set (LOCKED)

**Positive set** (9 expected NN pairings):

| Pair | Expected NN distance | Mechanism |
|---|---|---|
| Donepezil ↔ Galantamine | top-10 | Both AChE-I |
| Donepezil ↔ Rivastigmine | top-10 | Both AChE-I |
| MPH ↔ Atomoxetine | top-20 | DAT+NET vs NET-selective |
| MPH ↔ Modafinil | top-50 | DAT vs DAT-low-affinity |
| Modafinil ↔ Armodafinil/Caffeine | top-20 | Wake-promoting class |
| Memantine ↔ Ketamine | top-20 | NMDA antagonists |
| Memantine ↔ Amantadine | top-20 | NMDA + dopaminergic |
| Varenicline ↔ Galantamine | top-100 | α7 nAChR ligands |
| Donepezil+Memantine combo ↔ Donepezil | tight | Combination → donepezil-dominant |

**Negative control** (1): encenicline must rank > 500 from the active-cognition cluster centroid (the cardinal "target-true.phenotype-failed" anchor).

---

## 5. Four validation gates

### 5.1 Gate 1 (PRIMARY): mechanism-class recovery via clustering

**Clustering protocol** (LOCKED):
- SNN graph: k=15 nearest neighbors, Jaccard similarity
- **Leiden** (primary): γ ∈ {0.4, 0.6, 0.8, 1.0, 1.2}; pick best AMI on held-out half
- **HDBSCAN** (secondary): min_cluster_size ∈ {15, 25, 50}
- Metrics: AMI, ARI, V-measure, Fowlkes-Mallows

**Thresholds**:

| Band | AMI | ARI | Action |
|---|---|---|---|
| **PASS** | ≥ 0.50 | ≥ 0.40 | Enter joint posterior at λ_phen = 1.0 |
| **DEGRADE** | [0.30, 0.50) | [0.25, 0.40) | Enter at λ_phen = 0.5; flag in report |
| **FAIL** | < 0.30 | < 0.25 | Publish negative result |

**Stratification** (also reported):
- Per mechanism class (which classes cluster cleanly, which fail)
- By data modality (L1000-only vs JUMP-CP-only vs MEA-only vs MOFA+ joint) — joint AMI must exceed best single modality by ≥0.05
- chemCPA-imputed vs observed (inflate τ_chemCPA if imputed AMI << observed)

### 5.2 Gate 2 (SECONDARY): held-out g prediction

**Model**: Gaussian process regression (Matérn-5/2 kernel) on MOFA+ K-dim factors. Bayesian NN (numpyro + flax, ~50K params, HMC) as robustness check.

**Split**: 80/20 stratified by mechanism class; LOMCO (Leave-One-Mechanism-Class-Out) as the harder generalization test.

**Thresholds**:

| Band | MAE | 90% CrI coverage | Action |
|---|---|---|---|
| **PASS** | < 0.20 | ≥ 85% | V8 → V7 effect-size translation validated |
| **DEGRADE** | [0.20, 0.35] | [70%, 85%) | Use V8 → V7 but flag β_P uncertainty |
| **FAIL** | > 0.35 | < 70% | β_P = 0 (V8 ranking only) |

### 5.3 Gate 3 (SANITY): nootropic-anchor connectivity

**Threshold**: ≥ 7 of 9 expected positive NNs hold AND encenicline ranks > 500 from active-cluster centroid.

**Failure on single canonical pair** (e.g. donepezil ≁ galantamine) → integration broken → halt.

### 5.4 Gate 4 (NOVELTY): I_novel + (L, L, H) cell identification

**Threshold**: ≥ 6 of 8 held-out novel-mechanism anchors (clemastine + BIMA-8 cluster: benztropine + atropine + ipratropium + oxybutynin + trospium + tiotropium + quetiapine + PIPE-307) in top-5% I_novel rank.

---

## 6. Sensitivity analyses (pre-registered)

| Parameter | Sweep values | What it tests |
|---|---|---|
| K (MOFA+ factors) | {20, 30, 40, 50} | Latent-dim invariance |
| Leiden γ | {0.4, 0.6, 0.8, 1.0, 1.2} | Resolution sensitivity |
| HDBSCAN min_cluster_size | {15, 25, 50} | Cluster-size robustness |
| λ_phen | {0, 0.25, 0.5, 1.0, 2.0} | Phenotype-axis weight sweep |
| Embedding method | MOFA+ vs Biolord vs Deep CCA | Architecture invariance |
| τ_chemCPA | {0.5, 1, 2, 4} | Imputation-uncertainty robustness |

**Acceptance**: top-50 shortlist overlap ≥ 60% across λ_phen ∈ [0.5, 1.0]; ≤ 90% (V8 is doing something).

---

## 7. Cell-line transfer (U2OS-to-brain) — pre-registered handling

The single biggest theoretical leap is using U2OS Cell Painting morphology as a proxy for brain biology. Engineering response (pre-registered):

1. **MOFA+ per-modality variance decomposition** — report % of compound variance explained by JUMP-CP-only factors vs L1000 vs MEA/snRNA.
2. **Gate 1 AMI three ways** — all-modality joint, neural-only (L1000 NEU/NPC + MEA + snRNA), non-neural (JUMP-CP U2OS).
3. **If non-neural-only AMI ≥ 0.5** → U2OS Cell Painting transfers meaningfully → publish the finding.
4. **If non-neural < 0.3 alone but joint > 0.5** → U2OS contributes as chemistry-anchored consistency check, not brain proxy → frame accordingly.
5. **Per-cell-line random effect** in joint posterior to soak up cell-context confound.
6. **Honest framing in Discussion / Limitations**.

---

## 8. Falsifiability + pre-registration timing

OSF.io project will be locked **BEFORE**:

1. Running MOFA+ on the full 7-view stack
2. Computing AMI/ARI against any mechanism-class labels
3. Looking at any Gate-1/2/3/4 result
4. Examining the I_novel ranking
5. Identifying the (L, L, H) cell members

The pre-registration lock includes:
- K = 30 primary; sweep {20, 30, 40, 50}
- Leiden γ ∈ {0.4–1.2}; HDBSCAN min_cluster_size ∈ {15, 25, 50}
- AMI cuts 0.5 / 0.3
- ARI cuts 0.4 / 0.25
- ~30-class mechanism list
- 9+1 nootropic-anchor set
- Gate 2 MAE cut 0.20; CI coverage cut 85%
- Gate 4 I_novel threshold 6/8 in top-5%

---

## 9. Publication plan

**Primary** (A realistic, AMI ≥ 0.50): *Nature Machine Intelligence* methodology paper.
- Title (draft): "πphen: a target-agnostic multi-modal phenotypic prior for Bayesian cognition-enhancement drug repurposing."

**Stretch** (A+, AMI ≥ 0.60): *Nature Methods*.
- Title (draft): "Multi-modal phenotypic integration of LINCS L1000 + JUMP-CP Cell Painting + iPSC-MEA + chemCPA imputation for novel-mechanism drug discovery."

**Fallback** (Gate 1 DEGRADE): *Cell Reports Methods* / *Cell Systems* / *Bioinformatics*.

**Negative-result fallback** (Gate 1 FAIL, AMI < 0.30): *Bioinformatics* short report or PLOS ONE.
- Title (draft): "Multi-modal phenotypic integration fails to recover mechanism-class structure for cognition-enhancement compounds: a pre-registered falsification."

All submissions will include:
- The OSF lock URL (proof of pre-registration)
- The 9+1 anchor set + 30-class taxonomy parquet
- All `src/mammal_repurposing/cluster_e/*` code under permissive open-source license
- Replication instructions for the full V6.A → V6.B → V7 → V8 chain

---

## 10. Caveats + limitations

1. **L1000 cell-line bias** — mostly cancer lines; NEU/NPC L1000 thin. Mitigation: cell-line random effect; downweight non-neural for cognition lookups.

2. **JUMP-CP U2OS bias** — the elephant. Defensible only if (a) factor decomposition shows morphology captures target-class signal rather than cell context; (b) Gate 1 AMI on JUMP-CP-only embeddings exceeds 0.3 for AChE-I, NMDA antagonist, and DAT/NET classes (canonical Cell Painting wins). If those classes fail in U2OS, JUMP-CP enters at λ_phen = 0.25 instead of 1.0.

3. **iPSC neuron data is patchy and donor-heterogeneous** — protocol-dependent (Brennand, Pașca, Studer, Mariani). Donor random effect; protocol covariate.

4. **chemCPA scaffold-extrapolation failure** — RDKit Morgan-FP has known weakness for novel scaffolds outside Lipinski-compliant chemistry; flag and downweight.

5. **Phenotypic signatures are MOA proxies, not mechanism proofs** — a compound that "looks AChE-I-class" in L1000 may not engage AChE; it may engage a downstream pathway producing similar transcription. Standard CMap caveat (Lamb 2006); acceptable for prioritization, not mechanism claims.

6. **Mechanism-class clustering circularity** — handled by LOMCO + pre-registration. A reviewer may still claim mechanism-class labels are loosely derived from phenotype-adjacent literature. Honest framing in Discussion.

7. **Healthy-adult cognition effects are small** — Roberts 2020 ceiling g ≈ 0.43 (MPH delayed recall); overall means 0.12–0.21. V8's resolution to distinguish active from inactive may be saturated by SNR; Gate 2 MAE < 0.20 may be hard.

8. **No wet-lab dollars; no prospective validation** — all claims in-silico against published anchors. Pipeline produces a *prioritized shortlist*, not a validated discovery.

9. **Species mismatch** — mouse-vs-human iPSC-MEA; species in LINCS minor but present.

10. **Compute envelope** — ~24-36 h wall-clock on RTX 5070 12 GB + 32 GB RAM + WSL2 (per V8 plan §J.3). LINCS+JUMP-CP raw data ~55 GB cache; the ~115 TB cpg0016 raw images are NEVER downloaded.

---

## 11. OSF + AsPredicted template fields

| OSF field | Value |
|---|---|
| Title | V8 πphen Perturbational Evidence Axis (MAMMAL Cognitive Enhancement Drug Repurposing) |
| Authors | Pierce Lonergan + Claude Opus 4.7 (1M context) |
| Lock timestamp | TBD |
| Public release | Upon publication |
| Data sources | GEO GSE92742 + GSE70138 + clue.io beta; AWS s3://cellpainting-gallery/cpg0016-jump; Frank 2017 + Odawara 2016 + Hyysalo 2017 MEA datasets; cellxgene-census brain organoid slice |
| Hypothesis 1 (Gate 1 PRIMARY) | (See §1) |
| Sample size | ~50K compounds via chemCPA imputation; ~50–100 anchor compounds for Gate 2 |
| Outcomes | AMI/ARI vs PRISMA labels; held-out Hedges' g MAE + 90% CrI coverage; nootropic-anchor NN structure; I_novel score for novel-mechanism identification |
| Analysis plan | (See §2-§7) |
| Pre-registered predictions | Gate 1 thresholds + 9+1 anchor set + ~30-class taxonomy + I_novel novel-mechanism floor |
| Falsification thresholds | (See §5) |
| Sensitivity analyses | (See §6) |

---

**This pre-registration is V8.4 Stage 1.** It locks the MOFA+ K, Leiden γ space, AMI/ARI bands, mechanism class list, nootropic anchor set, MAE/CrI thresholds, and I_novel novel-mechanism gate BEFORE the V8.4 Stage 2 (executed-on-real-data) MOFA+ fit + clustering + AMI computation. Any deviation from this pre-registration in the final V8 paper must be flagged as a *post-hoc* analysis.

---

*Generated by `reports/paper-drafts/v8_osf_preregistration.md`. Companion to `design/architecture-and-plans/V4_STATUS_AND_FORWARD_PLAN.md` §13.Z and `design/architecture-and-plans/V6_ARCHITECTURE_PLAN.md` V8 section.*
