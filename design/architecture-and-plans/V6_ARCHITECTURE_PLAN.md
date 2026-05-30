# V6 / V7 / V8 Architecture & Phased Implementation Plan

**Status**: live source-of-truth for the V6 → V7 → V8 design arc. Companion to
`design/architecture-and-plans/V4_STATUS_AND_FORWARD_PLAN.md` §13 (V5/V6 Path Forward) + §13.Y (V7)
+ §13.Z (V8). Concrete implementation roadmap for **four** workstreams:

- **V6.A — Multi Head DTI ensemble** (~12 weeks) per
  `research/4-tier/Multi Head DTI.md`
- **V6.B — Bayesian Cluster D neurobiological prior** (~16 weeks) per
  `research/4-tier/Multi-Source Neurobiological Prior for Cognition Target Prioritization.md`
- **V7 — Clinical Effect-Size Translation Function** (~3-4 months) per
  `research/4-tier/Clinical Effect-Size Translation Function.md` + the
  companion Pre-Registration Methodology doc. *Downstream consumer of V6.A
  pchembl posteriors + V6.B θ̄ via multiplicative Cluster D gate.*
- **V8 / Cluster E — πphen Perturbational Evidence Axis** (~22 weeks) per
  `research/4-tier/Perturbational Evidence Axis.md` + the Technical
  Feasibility Deep-Dive companion. *Parallel third Bayesian factor; target-
  agnostic phenotypic prior to V6's target-first axes.*

The V6 scaffolds for V6.A + V6.B are already shipped in `src/mammal_repurposing/`
(diagnostics/per_head_bias.py, fusion/bayesian_router.py,
cluster_d/bayesian_prior.py, cluster_d/data_fetchers.py, cluster_a/
{mmatt_dta_adapter.py, psichic_adapter.py}, calibration/venn_abers.py). The
implementation work below operationalises them with real heads + real data,
then layers V7 (translation) and V8 (πphen) on top.

---

## V6.A — Multi Head DTI Ensemble (12 weeks)

### Goal
Replace the V5 2-head DTI signal (MAMMAL calibrated + Tanimoto) with a
5-head mixture (MAMMAL + Tanimoto + MMAtt-DTA + PSICHIC + BALM) with
explicit bias decomposition, per-target Bayesian routing, eMOSAIC OOD
gating, calibrated uncertainty propagation, and a disagreement-as-signal
discovery facet.

### Pre-committed Tier-A criterion
At SLC6A3 + SLC6A2 the ensemble must beat the +0.90/+0.91 Tanimoto floor by
≥0.01 each AND not regress at SLC6A3. Failure → **Tier-B fallback**:
production stays at 3-head (MAMMAL + Tanimoto + Cluster D); the negative
finding is the publishable contribution.

### ⚠️ V6.A.1 EMPIRICAL RESULT (this sprint): Tier-A FAILS at SLC6A3

- **Measured MMAtt-DTA ρ at SLC6A3 = +0.65** (pre-committed +0.78; Tanimoto +0.90)
- **Tier-A criterion: FAIL** — MMAtt-DTA does not beat Tanimoto at transporter
- **Tier-B fallback triggered**: 3-head ensemble (MAMMAL+Tanimoto+PrimeKG)
  is production; MMAtt added as 4th *conditional* ranker via INVERT mask
  (kept at 13 targets where ρ > +0.15; dropped at 6 INVERT targets)
- **GPCR/PDE wins are real**: HRH3 +0.82, HCRTR2 +0.70, PDE4D +0.39 —
  superfamily-conditional architecture works as Schulman 2024 claimed
- **Reframed publishable contribution**: per-target Bayesian router is now
  empirically necessary (uniform ensembling would degrade SLC6A2/ADRA2A/CHRNA7)
- See `reports/pipeline/mmatt_dta_activation_v1.md` for the full empirical table

### Phased plan

| Phase | Weeks | Deliverable | Validation |
|---|---|---|---|
| **V6.A.1 Heads installed** | 1-3 | MMAtt-DTA (pip + Zenodo weights ~2 GB; adapter shipped at `cluster_a/mmatt_dta_adapter.py`); PSICHIC; BALM (ESM-2 + ChemBERTa-2, ~3 GB) | Each head produces per-target ρ vs ChEMBL pchembl≥8 truth |
| **V6.A.2 Bias decomposition** | 4-5 | Wire `diagnostics/per_head_bias.py` to compute PC_k, SN_k, OOD_k, CT_k per (head, target). Bonett-Wright CIs | 5-head × 22-target trust matrix T(t, k) ∈ [0.02, 0.7] with row entropy logged |
| **V6.A.3 Bayesian router** | 6-7 | `fusion/bayesian_router.py` is shipped; activate by passing real T(t, k); add identifiability diagnostic report | Per-target router weights logged; identifiability theorem confirms n*=720 >> v4 n=7-26 (priors, not posteriors) |
| **V6.A.4 Calibrated uncertainty** | 8 | Per-head Venn-ABERS (Mervin 2020 J Chem Inf Model 60:4546). Cross-head correlation matrix Σ_kk' from 133-tuple calibration set. Replace router's Gaussian-CI with VA MC propagation | CI width inflation factor √(1+(K-1)·r̄) ≈ 1.41 confirmed |
| **V6.A.5 Disagreement facet** | 9-10 | Extend `35_v3_disagreement_signal.py` to multi-head: pairwise Kendall-τ + rank-distance + facet-tag {novel_scaffold / activity_cliff / ood / noise} | `reports/pipeline/disagreement_axis_v1.md` per-compound bucket |
| **V6.A.6 Validation + paper** | 11-12 | Run hypothesis audit; Tier-A criterion check. If PASS → J Cheminform / Nat Mach Intell draft. If FAIL → publish the negative result + 3-head fallback architecture | Pre-committed predictions in §13.1 |

### Pre-committed predictions (per Multi Head DTI.md §0)

| Target (n) | MAMMAL cal. | Tanimoto | MMAtt-DTA | PSICHIC | BALM | Ensemble (router) |
|---|---|---|---|---|---|---|
| SLC6A3 (n=26) | −0.70 | +0.90 | +0.78 | +0.74 | +0.62 | **+0.91 [+0.81,+0.96]** |
| SLC6A2 (n=23) | −0.60 | +0.91 | +0.80 | +0.75 | +0.65 | **+0.92 [+0.82,+0.96]** |
| ACHE (n=24) | +0.24 | +0.81 | +0.72 | +0.78 | +0.55 | +0.84 [+0.66,+0.93] |
| DRD1 (n=21) | +0.29 | +0.85 | +0.85 | +0.84 | +0.60 | +0.88 [+0.72,+0.95] |
| HCRTR1 (n=18) | +0.37 | +0.78 | +0.80 | +0.82 | +0.55 | +0.84 [+0.62,+0.94] |

### Dependencies + risks

- **MMAtt-DTA install**: manual `git clone` + Zenodo weights download
  (~2 GB). Documented in `cluster_a/mmatt_dta_adapter.py`. **Risk**: weights
  Zenodo DOI may rot — pin commit hash.
- **PSICHIC + BALM**: pip-installable but BALM needs ESM-2 cache.
  **Risk**: BALM may not reach Tier-A at SLC6A3 (predicted +0.62) — fallback
  is to use BALM as a tiebreaker only.
- **Compute**: 5 heads × 12k library × 22 targets at inference ≈ 1-2 hr on
  RTX 5070. No training required.

### Falsifiability fallback

If MMAtt-DTA / PSICHIC / BALM cannot beat Tanimoto +0.90 at the transporters,
the negative result is publishable as a methodology contribution. The
architecture stays at the 3-head (MAMMAL + Tanimoto + Cluster D)
configuration, and Cluster D's behavioural anchor (Roberts 2020 ceiling)
becomes the primary V6 deliverable.

---

## V6.B — Bayesian Cluster D Neurobiological Prior (16 weeks)

### Goal
First **behavioural anchor** in the pipeline. Full PyMC NUTS hierarchical
model over (AHBA, OT Genetics L2G, cellxgene-census single-cell, Lit-OTAR)
with explicit credible intervals, Jensen-Shannon disagreement axis, and a
hard Roberts 2020 SMD ceiling gate.

### Goal
- Replace implicit "cognition = binding proxy" with three independent
  neurobiological evidence streams + behavioural validation gate
- Expand panel from 22 to ~210 GWAS-anchored targets
- Provide posterior credible intervals on every target's cognition relevance

### Pre-committed verdict structure (per Cluster D §H)

| Target | y^AHBA | y^L2G | y^SC | θ̄ [90% HDI] | D_i | Verdict |
|---|---|---|---|---|---|---|
| BDNF | +0.65 | +0.55 | +0.70 | +0.78 [0.62, 0.93] | 0.08 | Three-source agreement |
| HTR2A | +0.55 | +0.10 | +0.30 | +0.35 [0.05, 0.65] | **0.62** | **High-disagreement positive — exactly the framework's target** |
| CHRNA7 | +0.18 | +0.05 | +0.45 | +0.22 [-0.05, +0.50] | 0.48 | SC drives; cortical AHBA undersamples α7's hippocampal niche |
| ACHE | +0.05 | +0.10 | +0.10 | +0.10 [-0.10, +0.30] | 0.10 | **Substrate-mediated flag — framework limitation** |

### Phased plan

| Stage | Weeks | Deliverable | Validation gate |
|---|---|---|---|
| **V6.B.1 Foundation** | 1-3 | `abagen.get_expression_data()` (pinned: ibf_threshold=0.5, probe_selection='diff_stability', donor_probes='aggregate', etc.) → AHBA cache. `BrainSMASH` 10k surrogates. OT Genetics L2G GraphQL pull (Davies 2018, Hill 2019, Sniekers 2017, Savage 2018, UKBB). cellxgene-census brain slice cached (Siletti 2023 + Mathys 2019 + Allen). `cluster_d/data_fetchers.py` scaffolds (shipped) wired to real APIs | **BDNF positive with three-source agreement** at θ̄ ≥ +0.5 |
| **V6.B.2 Panel expansion** | 4-6 | Generate ~210-target panel per §F (GWAS L2G≥0.2 OR MAGMA p<2.7e-6 OR AHBA \|r\|>0.3 BrainSMASH-corrected OR cell-type z>2 OR Lit-OTAR≥0.5). Validate existing 22-target panel + 44-target liability panel are strict subsets | **≥80% of new targets have a published modulator chemotype in ChEMBL** |
| **V6.B.3 Bayesian model** | 7-9 | Implement PyMC NUTS per §B.2 (already in `cluster_d/bayesian_prior.py::fit_cluster_d_prior_nuts`). 4 chains × 2000 warmup × 2000 draws on RTX 5070 via numpyro backend. Sensitivity sweep over θ / β / τ priors + Lit weight + reference anchors | **R̂ < 1.01, ESS > 400 per θ_i, zero divergences at target_accept=0.95; sign-stability >90% across sweep** |
| **V6.B.4 Validation gates** | 10-12 | **Gate 1 (HARD) Roberts SMD ceiling**: no target's predicted modulator effect-size posterior > Hedges' g = 0.5 at 90% credible upper bound. **Gate 2 Spearman**: per-target θ̄ correlates with meta-analytic SMD (Spearman ρ > 0.3) across ~15 reference compounds. **Gate 3 GWAS held-out**: AUROC > 0.7 on ABCD + CAC held-out. **Gate 4 leave-one-source-out**: Spearman ρ > 0.2 in all 3 folds | **All four gates pass**. If Gate 2 fails, audit substrate-mediated tagging |
| **V6.B.5 §7.11 integration + paper** | 13-16 | Plug into calibration via w^final_i = w^cal_i · σ(θ_i^post) · (1 + γ/(1 + HDI_width)). Re-run downstream Pareto. Paper draft | **Cell Reports Methods** (A+ fit) or **Bioinformatics** (A fit) |

### Threshold-driven rules
Per Cluster D §Recommendations:
- Sign-stability < 80% → downgrade from "primary prior" to "secondary diagnostic"
- Gate 2 Spearman ρ < 0.2 → do not publish; core empirical claim failed
- Median D_i > 0.7 → sources too inconsistent; fall back to single-source priors and publish disagreement as the main finding
- 90% HDI width > 0.6 for >50% of targets → framework inconclusive; expand reference-anchor set

### Dependencies + risks

- **abagen + BrainSMASH** (~1.4 GB AHBA download via abagen.fetch_microarray + brainsmash variogram). Risk: AHBA only has 6 donors; right hemisphere n=2 → leave-one-donor-out sensitivity reported.
- **OT Genetics L2G GraphQL** is rate-limited (~10 qps); cache locally. Endpoint at api.genetics.opentargets.org/graphql.
- **cellxgene-census** is network-bound (~10 GB local cache for human brain slice). tiledbsoma + cellxgene-census versions must match.
- **PyMC NUTS + numpyro JAX backend**: heavy install (~1 GB). RTX 5070 fine for T=210; gene-level T≈15,000 requires sparse approximation (out of V6 scope).

### Critical citation correction (already applied in V4 doc + this plan)
"Mansuri 2024 41-gene cognition map" was a misattribution. The correct
citation is **Moodie JE, Harris SE, Harris MA, Buchanan CR, Davies G, et al.
2024 *Hum Brain Mapp* 45(4):e26641** (doi:10.1002/hbm.26641, PMID 38488470,
PMC10941541). Internal references corrected throughout V4 doc + Appendix A.9.

---

## V6.A × V6.B Composition (the V6 joint posterior)

Per V4 §13.3, the V6 joint posterior over (compound, target) pairs is:

p(cognitive_relevance(q, t)) ∝ π(t | cognition) · Σ_k w_k(t) · F_k(q, t)

where:
- **π(t | cognition)** = Cluster D posterior at target t (from V6.B)
- **w_k(t)** = Multi Head DTI router weight from V6.A.3 trust matrix
- **F_k(q, t)** = Venn-ABERS-calibrated predictive distribution from head k

**Cluster D and the cross-DTI ensemble are independent factors** — additive
evidence assembly, not multiplicative double-counting.

The composition produces the V6 wet-lab shortlist:
1. Ranked by joint posterior mean with credible intervals
2. Pre-filtered by Roberts 2020 SMD ceiling
3. Annotated with both:
   - disagreement-axis facet-tag (V6.A.5)
   - Cluster D D_i Jensen-Shannon disagreement (V6.B.3)
4. Wet-lab priority = (high cross-DTI disagreement) × (high Cluster D posterior)

These are the high-information-value candidates that justify wet-lab spend.

---

## V7 — Clinical Effect-Size Translation Function (~3-4 months)

### Goal
The first **translational head** in the pipeline. Translate the V6.A pchembl
posteriors + V6.B Cluster D θ̄ posteriors + PBPK exposure into a *predicted*
healthy-adult cognition Hedges' *g* with credible intervals. Every prior
layer ranks compounds; V7 predicts the magnitude of the effect a wet-lab
experiment would actually measure.

### V7 in one sentence
**g_predicted(compound, endpoint) ~ Normal(η[compound, target] − Σ_k γ_k m_k, σ²)**
where η = sigmoid(α + β1·E[pchembl_post] + β2·E[relevance_post] + β3·copula_correction),
**β_target[t_c] = θ̄_{t_c} · β_raw_target[t_c]** (Cluster D multiplicative gate),
and the 5 failure-mode moderators m_k cover U-shape, practice/placebo,
tolerance, trait×state, trial-design.

### Pre-committed Tier-A criterion
8 pre-registered predictions P1–P8 (donepezil, encenicline-3mg-failure,
MPH-DSST, modafinil-200mg, memantine-20mg, intepirdine-MINDSET, pridopidine-
PROOF-HD, lecanemab-cognitive-subdomain) must all land within pre-registered
posterior bands. If ≥3 fail → **CPT:PSP negative-result paper** fallback.

### Phased plan

| Stage | Wk | Deliverable | Validation |
|---|---|---|---|
| **V7.1 PBPK foundation** | 1-4 | `src/mammal_repurposing/translation/pbpk.py` — JAX/diffrax 9-compartment ODE solver (gut → plasma → peripheral → cortex → striatum → hippocampus → basal-forebrain → brainstem → CSF). Watson 1989 receptor-occupancy-with-reserve formalism. U-shape generator (D1-postsynaptic vs D2-autoreceptor). Tolerance kinetics (R_avail dynamics) | All 3 PET anchors reproduced within 1σ: donepezil 19.1% cortical AChE (Bohnen 2005); MPH 12/40/54/72/74% DAT at 5/10/20/40/60mg (Volkow 1998); haloperidol D2 EC50 ~1.8 nM (Kapur 2000) |
| **V7.2 Class priors** | 5-8 | `src/mammal_repurposing/translation/prisma_priors.py` — Schmidli 2014 robust MAP priors for 12 mechanism classes (donepezil, modafinil, MPH, atomoxetine, memantine, vortioxetine, guanfacine, caffeine, piracetam, creatine, omega-3, minocycline) extracted from Roberts 2020 + MetaPsy + Cochrane | Each of 12 classes has prior mean ± σ from ≥3 trials; coverage table |
| **V7.3 Hierarchical Bayes** | 9-14 | `src/mammal_repurposing/translation/effect_size_model.py` — PyMC 3-level model: μ_global ~ N(0, 0.20); μ_class[m] ~ N(μ_class_PRISMA[m], λ_class·σ_class_PRISMA[m]); η = sigmoid(α + β1·E[pchembl_post] + β2·E[relevance_post] + β3·copula_correction); β_target gated by Cluster D θ̄; 5 moderators; numpyro backend | R̂ < 1.01, ESS > 400 per θ; sensitivity sweep over λ_class ∈ {0.1, 0.3, 1.0, 3.0}, robust MAP weight, moderator strengths |
| **V7.4 Validation gates** | 15-20 | **Gate 1 (HARD)**: All 8 P1–P8 posterior bands land. **Gate 2**: Roberts 2020 SMD ceiling — no compound's posterior 90% credible upper bound exceeds Hedges' g = 0.50. **Gate 3**: MAE on held-out anchor set (15 compounds) < 0.15. **Gate 4**: per-endpoint calibration plot 90% CrI coverage ≥ 85% across 6 endpoints (ADAS-Cog, DSST, n-back, Stroop, RAVLT, CANTAB-RVIP) | All 4 gates pass; otherwise CPT:PSP negative-result framing |
| **V7.5 OSF pre-registration + paper** | 21-24 | OSF.io lock on (priors, moderators, P1–P8, validation thresholds) BEFORE unblinding the held-out anchor set; pre-registered MS deposited on bioRxiv; submission to **Clinical Pharmacology & Therapeutics** (Wiley, IF 7.3) | OSF DOI minted; CPT submission ID |

### Dependencies + risks

- **JAX + diffrax + numpyro** stack; PyMC 5.x with numpyro backend. ~1 GB
  install. Mostly tested in V6.B.3, but V7 PBPK adds new ODE solver
  surface area.
- **PRISMA meta-analytic prior curation** is the load-bearing manual step:
  ~2 weeks for one engineer to extract 12-class SMD priors from
  Roberts 2020 + MetaPsy + Cochrane with proper variance components.
- **Risk**: V7 cannot fire until **both** V6.A.4 Venn-ABERS posteriors and
  V6.B.3 PyMC NUTS θ̄ posteriors exist as proper Bayesian objects (means + sds
  or full draws). V6.A.4 is shipped this sprint; V6.B.3 is queued.
- **Risk**: Roberts 2020 ceiling (g = 0.43 maximum significant subdomain
  effect for MPH delayed recall) is so low that V7's discriminative
  resolution may be at the noise floor. This is honest — the negative-result
  paper IS the publishable contribution if that's the empirical outcome.

### Falsifiability fallback
If ≥3 of P1–P8 fail, the V7 hierarchical translation is falsified for
healthy-adult cognition. Publishable contribution becomes the *negative
result paper* in **CPT: Pharmacometrics & Systems Pharmacology** (Wiley,
IF 4.2) — "PBPK + receptor-occupancy hierarchical Bayes cannot translate
in-silico DTI rankings to healthy-adult Hedges' g, because the Roberts
2020 ceiling is at the floor of the model's discriminative resolution."
Pre-registration protects against post-hoc retreat.

---

## V8 / Cluster E — πphen Perturbational Evidence Axis (~22 weeks)

### Goal
The **third Bayesian factor** parallel to V6.A (target-binding) and V6.B
(target-relevance). V8 introduces πphen, a *target-agnostic* phenotypic
evidence axis built from LINCS L1000 + JUMP-CP Cell Painting + iPSC-neuron
MEA / snRNA-seq + chemCPA generative imputation. The joint posterior
becomes:

**π_joint ∝ π_target(V6.A, V6.B) · π_phen(V8)**

with **Gaussian-copula correlation correction** between the three factors
(V6.A pchembl ↔ V6.B relevance ↔ V8 phenotype share training-data variance;
ignoring it produces overconfidence).

### V8 in one sentence
**Encenicline** is `target-true.phenotype-failed` — V6.A says binds α7
(pchembl ~7.8), V6.B Cluster D says α7 cognition relevance moderate (~0.55),
V8 says inert phenotypic signature (WTCS to active AChE-I cluster < 0.2,
JUMP-CP DeepProfiler cosine to AChE-I centroid < 0.3). V6 alone cannot
catch this; V8 does. This is the single motivating case.

### Pre-committed Tier-A criterion
Gate 1 (PRIMARY): mechanism-class recovery vs PRISMA ~30-class taxonomy with
**AMI ≥ 0.5, ARI ≥ 0.4**. OSF-pre-registered before unblinding mechanism
labels.

### Phased plan

| Stage | Wk | Deliverable | Validation gate |
|---|---|---|---|
| **V8.1 Data ingestion** | 1-3 | LINCS L1000 ingest (GSE92742 + GSE70138 + clue.io beta via cmapPy WTCS index, ~10 GB). JUMP-CP cpg0016 S3 sync — DeepProfiler + CellProfiler + DINOv2 consensus profiles only, ~30-40 GB (NOT the ~115 TB raw images). iPSC-MEA aggregation from Frank 2017, Odawara 2016, Hyysalo 2017. cellxgene-census brain slice (already V6.B asset, no incremental cost) | ≥60% V6 shortlist coverage in LINCS ∩ JUMP-CP overlap; if <60% chemCPA imputation becomes load-bearing |
| **V8.2 chemCPA training** | 3-4 | `src/mammal_repurposing/cluster_e/chemcpa_train.py` — RDKit-Morgan-FP-pretrained chemCPA on LINCS + sci-Plex3 (Hetzel 2022 architecture; 4-8 h GPU); LOMCO benchmark on glutamatergic class | R²(all) ≥ 0.50, R²(DEGs) ≥ 0.30 on broader-mix retraining (vs Hetzel 2022 ceiling 0.69/0.47, Piran 2024 cross-condition mean 0.51 ± 0.0062) |
| **V8.3 MOFA+ joint embedding** | 4-5 | `src/mammal_repurposing/cluster_e/mofa_embed.py` — mofapy2 K=30 across 7 views {L1000, CP_CellProfiler, CP_DeepProfiler, CP_DINO, MEA, snRNA, chemCPA_latent}; per-factor per-view variance attribution table (defends U2OS-to-brain transfer) | Joint AMI exceeds best single modality by ≥0.05 (multi-modal architecture earning its complexity) |
| **V8.4 Gate 1 pre-registration** | 5-6 | OSF.io lock on (K=30, Leiden γ space, AMI/ARI bands, 30-class mechanism list, 9+1 nootropic-anchor set) BEFORE unblinding. Leiden + HDBSCAN clustering; AMI/ARI computation; stratified per-class + per-modality | **Gate 1 PASS** at AMI ≥ 0.5, ARI ≥ 0.4; **DEGRADE** at [0.3, 0.5) → λ_phen=0.5; **FAIL** at <0.3 → negative-result paper |
| **V8.5 Joint posterior PyMC** | 6-8 | `src/mammal_repurposing/cluster_e/joint_phenotype.py` — extends V7 PyMC with V8 phenotype node: φ_c | θ_B, θ_T ~ Normal(A·b_c + B·r_c + C·(b_c⊗r_c), Σ_φ(τ_chemCPA, τ_cellline)). Extended g-mean with β_P·φ_c. Gaussian-copula correlation correction. Sensitivity sweep λ_phen ∈ {0, 0.25, 0.5, 1.0, 2.0} | Top-50 shortlist overlap ∈ [60%, 90%] across λ ∈ [0.5, 1.0]; reduces to V7-only at λ_phen=0 (sanity guard) |
| **V8.6 Gate 2 + 3 + 4 + paper** | 8-10 | GP regression (Matérn-5/2) on MOFA+ factors for held-out g prediction (Gate 2). Nootropic-anchor NN check on 9+1 set (Gate 3). I_novel novel-mechanism score correctly identifies (L, L, H) cell on held-out clemastine + BIMA-8 cluster (Gate 4). Mondrian conformal calibration via `crepes`. | All 4 gates pass; draft for **Nature Machine Intelligence** (A realistic) or **Nature Methods** (A+ stretch at AMI ≥ 0.6) |

### Five-MoA cognition reference centroids (per Technical Feasibility Deep-Dive)
K=5 sub-centroids in MOFA+ space — the cognition reference shape against
which compounds are scored for connectivity:

- **Cholinergic**: donepezil, galantamine, rivastigmine
- **Catecholaminergic**: MPH, atomoxetine, modafinil, d-amphetamine
- **Glutamatergic**: memantine, ketamine, riluzole
- **Trophic / ISR**: ISRIB, DNL343, 7,8-DHF, LM22A-4
- **Remyelination**: Mei 2014 BIMA-8 (clemastine + benztropine + atropine +
  ipratropium + oxybutynin + trospium + tiotropium + quetiapine); Najm 2015
  RNA-seq supplement; PIPE-307

Compounds clustering near the remyelination centroid with (L, L, H) profile
(low V6.A binding, low V6.B target relevance, high V8 phenotypic match)
are the **clemastine-class novel-mechanism candidates** — V8's central
pitch.

### 8-cell disagreement table
high/low × {target, genetic, phenotype}:
- (H, H, H) = canonical positive (donepezil, MPH)
- (H, H, L) = `target-true.phenotype-failed` (encenicline, intepirdine, pridopidine)
- **(L, L, H) = novel-mechanism territory** (clemastine, PIPE-307, BIMA-8)
- (L, H, L) = genetic-relevance-only (not actionable)
- ...

I_novel(compound) = π_p · [1 − I(π_p ; (π_t, π_g))] — mutual-information
novel-mechanism score; high when phenotype is informative AND target-genetic
axes are uninformative or independent.

### Dependencies + risks

- **Working set**: ~55 GB local cache (LINCS L1000 ~10 GB + JUMP-CP
  DeepProfiler+CellProfiler ~30-40 GB + iPSC-MEA ~10 GB + chemCPA models
  <2 GB). cpg0016 full image dataset (~115 TB) is **never downloaded**.
- **MOFA+** at K=30 across 7 views ~24 GB RAM peak; close to 32 GB ceiling.
  Mitigation: `mofapy2 sparse_data=True` + minibatching.
- **PyMC NUTS** on V8 joint at ~50K compounds × 30 factors ~8-16 h GPU via
  numpyro JAX backend.
- **U2OS-to-brain transfer** (the elephant): per-cell-line random effect +
  explicit downweighting of JUMP-CP-only factors in iPSC-deficient
  subspaces + Gate 1 stratified per modality.
- **chemCPA hallucination risk** for scaffolds with max-Tanimoto-to-train <
  0.3: τ_chemCPA × 3 inflation; flagged as `chemCPA.imputed.low_confidence`.

### Compute envelope
~24-36 h wall-clock (per Perturbational Evidence Axis.md §J.3) — RTX 5070
12 GB sufficient. Peak 24 GB RAM, 11 GB VRAM.

### Falsifiability fallback
- Gate 1 AMI < 0.3 → publish negative result for **Cell Reports Methods** or
  **Bioinformatics**.
- Gate 2 MAE > 0.35 → β_P = 0 (V8 contributes only to ranking, not to V7).
- Non-neural-only AMI < 0.3 → frame V8 as chemistry-anchored consistency
  check rather than brain proxy.

---

## V6 × V7 × V8 composition (three-factor joint posterior)

When V6.A + V6.B + V7 + V8 all ship, the joint posterior over Hedges' g per
(compound, target) is:

**p(g | compound, target, evidence)**
**∝ p_V8(g | compound, φ) · p_V7(g | compound, target, θ̄, b, PBPK) · p_V6.B(θ̄ | target) · p_V6.A(b | compound, target)**

with **Gaussian-copula correlation correction** between the three Bayesian
factors. The composition produces the **V8 wet-lab shortlist** (the
production deliverable once V8 lands):

1. Ranked by joint posterior mean Hedges' g with 95% CrI
2. Pre-filtered by Roberts 2020 SMD ceiling (no g > 0.50 at 90% upper CrI)
3. Annotated with:
   - V6.A multi-head disagreement-axis facet-tag (novel_scaffold /
     activity_cliff / ood / noise)
   - V6.B Cluster D D_i Jensen-Shannon disagreement
   - V8 three-way JSD I_novel novel-mechanism score
4. 8-cell disagreement classification: (H, H, H), (H, H, L), **(L, L, H)**
   [novel mechanism], etc.
5. Wet-lab priority = high I_novel × high V6.B posterior × high V8
   phenotypic match × passing Roberts ceiling

These are the **clemastine-class** high-information-value candidates that
justify wet-lab spend.

---

## Resource allocation decision tree

If a **single research-engineer-month** is available right now:
1. **BALM adapter** (V6.A.1 phase 3, 2-3 days) — completes the Multi-Head DTI
   activation. Already shipped: MMAtt-DTA + PSICHIC; BALM is the third
   non-MAMMAL head. Tier-A criterion at SLC6A3 was empirically falsified by
   V6.A.1; BALM may still earn an INVERT-mask slot at specific superfamilies.
2. **OT Genetics L2G fetcher + PyMC NUTS** (V6.B Stages 2-3, ~1-2 weeks).
   Activates the Cluster D Bayesian prior — V6.B foundation already shipped
   (abagen AHBA cache via `scripts/54_v6b_cluster_d_foundation.py`). The
   `cluster_d/bayesian_prior.py::fit_cluster_d_prior_nuts` scaffold fires the
   moment real (AHBA, L2G, SC) observations arrive.
3. **TxGNN API rewrite** (V4 last Tier-1 item; ~1 day). Switches from broken
   `predict_indication(drug, disease)` to public `predict_disease(idx)` API.
   Adds Cluster C TxGNN as a 5th fusion cluster.

If **2-3 months** are available — ship V6.A full + V6.B foundation. Push
production wet-lab shortlist to v10 with all 4 DTI heads + Cluster D prior.

If **4-6 months** are available — ship V6.A + V6.B in full, then **begin V7
PBPK foundation + class priors** (V7.1 + V7.2). V7 cannot fire without V6.A.4
Venn-ABERS posteriors (shipped) + V6.B.3 PyMC NUTS posteriors (V6.B.3).

If **8-10 months** are available — ship V6 + V7 in full. V7 paper draft for
Clinical Pharmacology & Therapeutics. Begin V8.1 data ingestion + V8.2
chemCPA training in parallel.

If **12-15 months** are available — ship V6 + V7 + V8 in full. Three papers,
three venues:
- V6.A → J Cheminform / Nat Mach Intell
- V6.B → Cell Reports Methods / Bioinformatics
- V7 → Clinical Pharmacology & Therapeutics (CPT, IF 7.3) or CPT:PSP (negative-result fallback)
- V8 → Nature Machine Intelligence (A realistic) / Nature Methods (A+ stretch)

If **no engineer time** is available — the current
`reports/wet-lab/wet_lab_shortlist_v6_full.md` is the production deliverable.
43 PASS compounds with all V4 + V5 gates flowing through. V6.A.1 MMAtt-DTA
empirical activation + V6.B.1 Cluster D foundation shipped this sprint and
ready to bolt on when engineering capacity resumes.

---

## V6 + V7 + V8 timeline summary

| Track | Wks | Effort | Dependencies | Output |
|---|---|---|---|---|
| V6.A: Multi Head DTI ensemble | 12 | 5 heads + bias decomposition + Bayesian router + eMOSAIC OOD + Venn-ABERS + disagreement facet + validation | MMAtt-DTA Zenodo download (shipped); PSICHIC adapter (shipped); BALM (pending phase 3) | `fusion/bayesian_router.py` (shipped) + `diagnostics/per_head_bias.py` (shipped) + `cluster_a/{mmatt_dta_adapter, psichic_adapter, balm_adapter}.py` + paper draft (J Cheminform / Nat Mach Intell) |
| V6.B: Bayesian Cluster D prior | 16 | abagen + BrainSMASH + OT L2G + cellxgene + PyMC NUTS + 4-gate validation + §7.11 integration + paper | abagen/BrainSMASH/PyMC installs (shipped); cellxgene brain slice ~10 GB; OT Genetics GraphQL fetcher (V6.B.2 pending) | `cluster_d/{bayesian_prior, data_fetchers}.py` (shipped scaffolds) + AHBA cache `data/results/v2/ahba_expression_v1.parquet` (shipped) + 5 validation reports + paper draft (Cell Reports Methods / Bioinformatics) |
| V6 Composition | 4 | Joint-posterior plumbing + wet-lab shortlist re-render | V6.A + V6.B both shipped | `reports/wet_lab_shortlist_v9_joint.md` |
| **V7: Clinical Effect-Size Translation** | **16** | **JAX/diffrax PBPK + PRISMA-anchored 3-level hierarchical Bayes + 5 failure-mode moderators + 4 validation gates + OSF pre-reg + paper** | **V6.A.4 Venn-ABERS posteriors (shipped); V6.B.3 PyMC θ̄ posteriors (V6.B.3 pending); PRISMA prior curation (~2 wk manual)** | **`translation/{pbpk, prisma_priors, effect_size_model}.py` + OSF pre-registration + CPT paper draft** |
| **V8 / Cluster E: πphen Perturbational Axis** | **22** | **LINCS L1000 + JUMP-CP + iPSC-MEA + chemCPA imputation + MOFA+ + conditionally-dependent 4-level hierarchical Bayes + 4 validation gates + OSF pre-reg + paper** | **V6.A + V6.B + V7 all shipped; ~55 GB local cache; mofapy2; cmapPy; pycytominer; numpyro JAX backend** | **`cluster_e/{ingest_lincs, ingest_jumpcp, chemcpa_train, mofa_embed, joint_phenotype}.py` + OSF pre-registration + Nat Mach Intell paper draft** |

**Total V6 + V7 + V8**: ~70 weeks (~16-18 months) of focused engineering across
all four workstreams. Four distinct papers, four distinct validation regimes,
four distinct publication venues. The shortlist that lands at the end is the
first cognition-enhancement candidate set in the literature with:

- formal credible intervals on every compound's rank (V6.A + V6.B + V8),
- behavioural validation gate (Roberts 2020 SMD ceiling, V6.B + V7),
- multi-head ensemble with disagreement-as-signal discovery axis (V6.A.5),
- three-axis JSD with I_novel novel-mechanism score (V8),
- 8-cell disagreement classification surfacing (L, L, H) clemastine-class candidates,
- mechanism + liability gating (V4 + V5),
- predicted Hedges' g per compound × endpoint with PBPK-grounded receptor occupancy (V7),
- AND a per-compound provenance trail back to documented signal sources with
  known failure modes.

That candidate set, not the next nootropic, is the contribution.

---

## V6 + V7 + V8 scaffold inventory

### Already shipped (V6.A + V6.B foundation)

| Module | Lines | Purpose | Status |
|---|---|---|---|
| `cluster_a/mmatt_dta_adapter.py` | ~200 | MMAtt-DTA adapter with 22-target superfamily map | LIVE (V6.A.1 activated) |
| `cluster_a/psichic_adapter.py` | ~190 | PSICHIC subprocess adapter (Koh 2024 Nat Mach Intell) with `_find_psichic_repo` env probe | LIVE (V6.A.1 phase 2) |
| `calibration/venn_abers.py` | ~200 | Venn-ABERS regressor + correlated_mc_intervals Gaussian-copula MC propagation | LIVE (V6.A.4) |
| `diagnostics/per_head_bias.py` | ~180 | PC/SN/OOD/CT signature computation + trust matrix builder | LIVE |
| `fusion/bayesian_router.py` | ~210 | Per-target router + OOD + confidence gates + identifiability diag | LIVE |
| `cluster_d/bayesian_prior.py` | ~250 | Stage-0 stub + PyMC NUTS full Bayesian model (Neelon-Dunson) | SCAFFOLD; activates with real (AHBA, L2G, SC) observations |
| `cluster_d/data_fetchers.py` | ~150 | AHBA / OT Genetics / cellxgene adapters with availability probes | LIVE for AHBA (`scripts/54_v6b_cluster_d_foundation.py` shipped); OT L2G pending |
| `analysis/brain_region.py` | ~140 | Static 22-target brain-bias map (V5 fallback / V6 preview) | LIVE |

**V6 LOC total**: ~1,520 lines of scaffold + tests across the V6 architectural footprint (was 1,130 in earlier snapshot; +PSICHIC adapter, +Venn-ABERS calibration).

### Pending V6 phase 3 / V7 / V8 modules

| Module (planned) | Purpose | When |
|---|---|---|
| `cluster_a/balm_adapter.py` | BALM fine-tuned ESM-2 + ChemBERTa-2 head (Gorantla 2025 JCIM 65(22):12279) | V6.A.1 phase 3 (this sprint) |
| `cluster_c/txgnn.py` (rewrite) | TxGNN with public `predict_disease(idx)` API for cognition EFO anchors | V6 Cluster C activation (this sprint) |
| `translation/pbpk.py` | JAX/diffrax 9-compartment PBPK ODE solver | V7.1 |
| `translation/prisma_priors.py` | Schmidli 2014 robust MAP priors for 12 PRISMA classes | V7.2 |
| `translation/effect_size_model.py` | PyMC 3-level hierarchical Bayes with Cluster D θ̄ multiplicative gate + 5 moderators | V7.3 |
| `cluster_e/ingest_lincs.py` | LINCS L1000 cmapPy WTCS index builder | V8.1 |
| `cluster_e/ingest_jumpcp.py` | JUMP-CP cpg0016 S3 sync + pycytominer normalization | V8.1 |
| `cluster_e/chemcpa_train.py` | RDKit-pretrained chemCPA on LINCS + sci-Plex3 | V8.2 |
| `cluster_e/mofa_embed.py` | MOFA+ K=30 joint embedding across 7 views | V8.3 |
| `cluster_e/joint_phenotype.py` | V7+V8 joint posterior with Gaussian-copula correlation correction | V8.5 |

### Activation triggers

- Set `MMATT_DTA_ROOT` env var → MMAtt-DTA adapter active *(SHIPPED)*
- Set `PSICHIC_ROOT` env var + matching conda env → PSICHIC adapter active *(SHIPPED)*
- `pip install pymc numpyro` → Cluster D Bayesian path active *(SHIPPED)*
- `pip install abagen brainsmash` → AHBA real-mode active *(SHIPPED via scripts/54)*
- Set `CRYPTOBENCH_HOME` / `POCKETMINER_HOME` → detector ensemble Sprint 2 active
- Set `BALM_WEIGHTS_DIR` → BALM adapter active *(pending V6.A.1 phase 3)*
- `pip install jax diffrax` → V7 PBPK active *(pending V7.1)*
- `pip install mofapy2 cmappy pycytominer` → V8 ingestion active *(pending V8.1)*

---

## How this doc was written

V6/V7/V8 plan synthesized from:
- `research/4-tier/Multi Head DTI.md` (~60 KB) — V6.A spec
- `research/4-tier/Multi-Source Neurobiological Prior for Cognition Target Prioritization.md` (~39 KB) — V6.B Cluster D Bayesian model
- `research/4-tier/Clinical Effect-Size Translation Function.md` (~32 KB) — V7 PBPK + hierarchical Bayes spec
- `research/4-tier/Clinical Effect-Size Translation Function A Methodology Pre-Registration for Bayesian Cognition-Enhancement Drug Repurposing.md` (~18 KB) — V7 companion: OSF pre-registration + P1–P8 falsifiers
- `research/4-tier/Perturbational Evidence Axis.md` (~50 KB) — V8 Cluster E LINCS + JUMP-CP + iPSC + chemCPA + MOFA+ + conditionally-dependent hierarchical Bayes
- `research/4-tier/Technical Feasibility Deep-Dive Adding a Phenotypic.md` (~24 KB) — V8 feasibility companion: I_novel score + 8-cell disagreement table + Mondrian conformal
- V4 plan §13 V5/V6 Path Forward + §13.Y V7 + §13.Z V8 + §13.W three-factor joint posterior (already integrated)
- This sprint's V6.A.1 + V6.A.4 + V6.B.1 + V7/V8 architecture-planning work

When V7/V8 work begins, the next assistant should:
1. Re-read all 6 research deep-dives in full
2. Activate the scaffolds via the env var / install instructions above
3. Finish V6.A (BALM adapter) + V6.B (OT Genetics + PyMC NUTS) FIRST — they
   are the prerequisites for V7 (which consumes both posteriors)
4. Run V7.1 PBPK foundation NEXT — it's the highest-leverage testable
   first slice of V7 (PET-anchored receptor occupancy is well-validated
   territory)
5. Run V8.1 LINCS+JUMP-CP ingestion in parallel with V7.1 — they share no
   dependencies and use disjoint compute resources
6. Maintain `reports/pipeline/hypothesis_audit_v1.md` as the standing falsifiability check
7. Honor OSF pre-registration timing for V7.5 + V8.4 — they must lock
   priors and validation thresholds BEFORE unblinding the held-out anchor
   sets

Generated 2026-05-26 (V6 scaffold commit), extended 2026-05-27 (V7 + V8
architecture planning sprint).
