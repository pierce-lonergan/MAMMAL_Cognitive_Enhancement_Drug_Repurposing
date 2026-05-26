# §7.9 Cluster D — Multi-Source Neurobiological Prior for Cognition Target Prioritization

**V3 Pipeline Deep-Research Deliverable · Pierce / MOMENTUM-X-style Repurposing Stack · 26 May 2026**

---

## TL;DR

- **Cluster D is the V3 pipeline's first behavioural anchor.** It collapses the §4.6 STRUCTURAL gap by replacing the "cognition = binding proxy" assumption with a **full Bayesian hierarchical target-weight model** fed by three independent neurobiological evidence streams (AHBA imaging-transcriptomics, OT Genetics L2G + intelligence GWAS, cellxgene-census single-cell), with **explicit credible intervals**, source-specific reliability weights, and disagreement-as-signal propagation. Rating: **A / publishable on its own merits** (Cell Reports Methods or Bioinformatics tier) if the validation gates pass.
- **The 44-target panel expands to ~210 GWAS-anchored targets**, subsuming §7.3 entirely. The expansion uses MAGMA-prioritized gene-level statistics from Davies 2018 (N=300,486; 148 independent loci) and Hill 2019 (N_eff=248,482 via MTAG; 187 loci implicating 538 genes) intersected with OT Platform druggability tiers, gated on brain expression and on Lit-OTAR cognition evidence. Roberts et al. 2020 (Eur Neuropsychopharmacol 38:40–62; k=47 trials, modafinil pooled SMD=0.12, methylphenidate SMD=0.21) sets the **honest empirical ceiling**: the model is not permitted to predict per-target effect sizes outside the (-0.1, +0.5) Hedges-g envelope without an explicit "regime-violation" flag.
- **Disagreement is treated as Pareto information, not noise.** When AHBA cortical-expression posterior, OT-Genetics L2G posterior, and single-cell cell-type posterior diverge for a given target (canonical case: 5-HT2A — strong AHBA, weak L2G, moderate single-cell), the per-source posteriors enter a Jensen-Shannon-divergence diagnostic that becomes its own Pareto axis alongside potency, selectivity, and liability. High-disagreement / high-mean targets are flagged as "high-information-value" — exactly the candidates where a Boltz-2 docking experiment is most worth running.

---

## Key Findings

1. **The "Mansuri 2024 41-gene cognition map" appears to be a misattribution.** The closest match in indexed literature (PubMed, bioRxiv, Google Scholar, as of 26 May 2026) is **Moodie et al. 2024, *Human Brain Mapping* 45(4):e26641 (doi:10.1002/hbm.26641; PMID 38488470; PMC10941541)**, which identifies **41 single genes** as candidate cortical spatial correlates of general cognitive function *g* (|β| range 0.15–0.53), beyond two dominant AHBA expression PCs, using a meta-analytic g-morphometry map from N=39,519 (UK Biobank + Generation Scotland STRADL + LBC1936). This is the canonical paper to anchor the prior. Pierce should rename §7.9.1 accordingly. **Rating of this correction: critical — silently citing "Mansuri 2024" risks reviewer rejection.**
2. **abagen (Markello et al. 2021, *eLife* 10:e72129) is the only defensible AHBA processing pipeline.** It exposes 17 documented processing options; Arnatkevičiūtė 2019 / Markello 2021 ablations show that gene-normalization and donor-aggregation choices dominate downstream inference. Hard-coding the abagen recommended defaults plus stability threshold r ≥ 0.2 yields ~12,668 stable genes on Desikan-Killiany. A **sensitivity sweep** over normalization options is a publication deliverable.
3. **BrainSMASH (Burt et al. 2020, *NeuroImage* 220:117038) is the right null model for AHBA spatial correlations.** Variogram-matched surrogates correctly preserve spatial autocorrelation; per Markello & Misic (2021, *NeuroImage*, doi:10.1016/j.neuroimage.2021.118052), "naive null models that do not preserve spatial autocorrelation consistently yield elevated false positive rates and unrealistically liberal statistical estimates" — typically inflating FPR 2–10× for autocorrelated brain maps. Spin tests (Alexander-Bloch 2018, *NeuroImage* 178:540–551) are acceptable for cortical surface data but have a known medial-wall artifact and no subcortical extension. **Decision: BrainSMASH primary, spin test secondary check.**
4. **Roberts et al. 2020 (*Eur Neuropsychopharmacol* 38:40–62) is the behavioural ceiling paper.** k=47 trials of healthy adults: modafinil pooled SMD = 0.12 (p=.01); methylphenidate SMD = 0.21 (p=.0004), driven by recall SMD = 0.43, sustained attention SMD = 0.42, inhibitory control SMD = 0.27; d-amphetamine null. The authors' own conclusion: "*data with these stimulants is far from positive if we consider that effects are small, in experiments that do not accurately reflect their actual use in the wider population. There is a user perception that these drugs are effective cognitive enhancers, but this is not supported by the evidence so far.*" **In healthy adults, no approved cognitive enhancer exceeds Hedges' g ≈ 0.5 in any sub-domain.** Any pipeline predicting per-target SMD > 0.5 is, with high prior probability, fitting noise.
5. **OT Genetics L2G provides a calibrated 0–1 probability of causal-gene assignment per credible set.** Shapley-explainable XGBoost trained on a curated gold-standard (Mountjoy et al. 2021, *Nat Genet* 53:1527–1533; refined 29-feature model on the Platform 2024+). Threshold conventions: L2G ≥ 0.5 high-confidence causal; ≥ 0.05 "in the room" (the Platform's default display threshold). **For cognition we will use credible-set-level L2G aggregation across Davies 2018, Hill 2019, Sniekers 2017, Savage 2018, and the UK Biobank cognitive-test panel.**
6. **Lit-OTAR (Kafkas et al. 2024, *Bioinformatics* 41(4):btaf113) processed > 39M abstracts + 4.5M full-text articles, producing > 48.5M unique associations** (29.9M target–disease, 11.8M target–drug, 8.3M disease–drug). The cognition slice is AD-dominated; we explicitly down-weight Lit-OTAR relative to L2G to mitigate confounding.
7. **cellxgene-census (CZ CELLxGENE, 2025-11-08 LTS, schema v2.4.0; 1,845 datasets) contains 162,025,130 total *Homo sapiens* cells (99,633,637 unique) plus 46,299,127 *Mus musculus* cells**, exposed via tiledbsoma ≥ 1.15.3. WSL2 deployment on Pierce's RTX 5070 workstation is straightforward (TileDB Linux wheels). Query pattern: `axis_query(tissue_general=="brain", is_primary_data==True)` with chunked PyArrow iteration. Brain primary nuclei in the human slice number in the high single-digit millions across cortex, hippocampus, basal ganglia, cerebellum, and brainstem subsets.

---

## Details

### A) Literature Foundation

**A.1 Moodie 2024 (the "41-gene cognition map", corrected from "Mansuri 2024")**

Moodie JE, Harris SE, Harris MA, Buchanan CR, Davies G, et al. *Hum Brain Mapp* 2024;45(4):e26641. Methodology in three layers:

- **Cortical g map.** Regional meta-analytic map of *g*–morphometry partial correlations, computed for each of 68 Desikan-Killiany regions × {cortical volume, surface area, thickness}, meta-analyzed across UK Biobank, STRADL, and LBC1936 — total meta-analytic N = 39,519.
- **AHBA expression matrix.** 8,235 stable genes × 68 DK regions, processed via abagen-equivalent pipeline (stability filter r ≥ 0.2). PCA across genes; first two components explain **49.4 % of regional expression variance** (C1 = cell-signalling/modifications axis; C2 = transcription-factor axis).
- **Residual single-gene regression.** For each of 8,235 genes, partial-correlate cortical expression with the *g*-cortical map after regressing out C1 and C2. **41 genes survive at standardized-beta threshold |β| ∈ [0.15, 0.53] with spin-test-corrected significance.**
- **Strengths**: out-of-sample replicated; explicit control for "first two PCs eat everything" confound; meta-analytic N for cognition map is large.
- **Limitations**: AHBA N=6 donors; right-hemisphere N=2; the 41-gene list is sensitive to PC-residual modelling choice; g is psychometric, not necessarily what Roberts 2020 captures pharmacologically.

**A.2 abagen toolbox (Markello/Hansen 2021)**

Markello, Hansen, Liu, et al., *eLife* 2021;10:e72129. github.com/rmarkello/abagen. 17 documented pipeline options. Pinned configuration: `ibf_threshold=0.5`, `probe_selection='diff_stability'`, `donor_probes='aggregate'`, `lr_mirror='bidirectional'`, `sample_norm='scaled_robust_sigmoid'`, `gene_norm='scaled_robust_sigmoid'`, `region_agg='donors'`, `agg_metric='mean'`, plus post-hoc stability filter r ≥ 0.2 (Hawrylycz 2015).

**A.3 Burt 2020 BrainSMASH**

Burt JB, Helmer M, Shinn M, Anticevic A, Murray JD. *NeuroImage* 2020;220:117038. github.com/murraylab/brainsmash. Generates surrogate maps with matched variogram, preserving spatial autocorrelation while randomizing alignment. Optimal `knn` setting per follow-up work (Markello et al. 2022, *NeuroImage* 257:119323) is `knn = n_vertices`, NOT the default 1,000 — materially reduces false-positive inflation.

**A.4 Alexander-Bloch 2018 spin test**

Alexander-Bloch AF, Shou H, Liu S, Satterthwaite TD, Glahn DC, Shinohara RT, Vandekar SN, Raznahan A. *NeuroImage* 2018;178:540–551. Spherical-rotation null. Limitations (Markello 2022, Vandekar 2021 SPICE): (i) ad-hoc medial-wall handling; (ii) cortical-surface-only; (iii) spherical-projection distance distortion. Use as cortex secondary; not for subcortex.

**A.5 Davies 2018 + Hill 2019 intelligence GWAS**

- Davies G, Lam M, Harris SE, et al. *Nat Commun* 2018;9:2098. **N = 300,486; 148 independent genome-wide significant loci** (CHARGE + COGENT + UK Biobank meta-analysis). SNP-h² ≈ 0.25.
- Hill WD, Marioni RE, Maghzian O, et al. *Mol Psychiatry* 2019;24:169–181. **MTAG raised functional sample size from N=199,242 to N_eff=248,482; identified 187 independent loci implicating 538 genes** (positional + eQTL + chromatin interaction). Neurogenesis and myelination enrichment. (The frequently-cited '939 gene' figure conflates this with Savage et al. 2018, *Nat Genet* 51:404–413.)
- PGS for intelligence currently explains R² ≈ 4.4–4.8 % of phenotypic variance (Davies 2018: R²=4.37 %; Savage 2018 PGS up to 4.81 % across four independent cohorts; meta-analytic correlation ρ ≈ 0.245 per Bouter et al. 2023, *Intelligence* 97:101734). Educational-attainment PGS can exceed 10 %, but is the *wrong target* for cognition-direct repurposing. **This R² ≈ 5 % is the structural upper bound on what genetics alone can tell us** — comparable to Roberts 2020's pharmacological ceiling.

**A.6 OT Genetics L2G**

Mountjoy E, Schmidt EM, Carmona M, et al. *Nat Genet* 2021;53:1527–1533, with refined 29-feature Shapley-decomposed model in the OT Platform. L2G ∈ [0,1]; calibrated against gold-standard. Inputs: SuSIE fine-mapping (in-sample LD weighted 1.0, out-of-sample 0.75; PICS 0.5/0.25), distance-weighted credible-set TSS distance, molecular-trait colocalisation, VEP, ABC-style E2G regulatory-region scores. Study-specific: must be computed per cognition GWAS.

**A.7 Lit-OTAR**

Kafkas Ş, Hulme A, Hsu YH, et al. *Bioinformatics* 2024;41(4):btaf113. **48.5M total associations**. For cognition slice (EFO:0003917 + descendants + MMSE/ADAS-Cog EFO terms): expect ~ 40–80 k target–"cognition"-family associations. AD literature dominates — confirmation-biased; we down-weight.

**A.8 cellxgene-census**

CZ CELLxGENE Discover Census LTS 2025-11-08, schema v2.4.0, 1,845 datasets (per chanzuckerberg.github.io/cellxgene-census release notes). Cell-based slicing via TileDB-SOMA. Brain datasets relevant to cognition:
- **Siletti et al. 2023** (*Science*, doi:10.1126/science.add7046) — 3,369,219 nuclei from 105 dissections spanning forebrain, midbrain, and hindbrain; 3,313 subclusters across 30 superclusters and 461 clusters.
- **Allen Brain Map** adult human cortex (10x multiome).
- **Mathys 2019** (*Nature* 570:332–337) — 80,660 single-nucleus transcriptomes from prefrontal cortex of 48 individuals with varying AD pathology, ROSMAP cohort.
- **Sieberts 2020** / PsychENCODE bulk pseudo-bulk.
- **BICCN** whole-brain mouse + human atlases.

Cell ontology lets us aggregate to cognition-salient classes: glutamatergic L2/3 IT, L5 ET pyramidal, hippocampal CA1/CA3/DG, basal forebrain cholinergic (CHAT+), LC noradrenergic (DBH+), raphe serotonergic (TPH2+), VTA DA (TH+), cortical PV+/SST+/VIP+ interneurons.

**A.9 Neurosynth**

Yarkoni 2011 + nimare API. We use the **123 cognitive-atlas-derived term maps**. "Cognition" composite = median over the cognitive-atlas family (memory, attention, working memory, executive, reasoning, fluid intelligence). Same map family used by Hansen 2022 (*Nat Neurosci* 25:1569–1581) for receptor–cognition decoding.

---

### B) Methodology — Full Bayesian Hierarchical Model

#### B.1 Notation

Let $i = 1, \dots, T$ index targets ($T \approx 210$). For each target define three observed evidence quantities:

- $y^{\text{AHBA}}_i \in [-1, 1]$: BrainSMASH-corrected spatial correlation between target gene expression (AHBA, abagen-processed, DK-68) and the Moodie-2024-style g-cortical map (or Neurosynth "cognition" composite). SE $\sigma^{\text{AHBA}}_i$ from bootstrap.
- $y^{\text{L2G}}_i \in [0, 1]$: Max L2G score over credible sets where target $i$ is in the locus, across Davies-2018, Hill-2019, Sniekers-2017, Savage-2018, and UK Biobank cognitive-test GWASes. SE $\sigma^{\text{L2G}}_i$ from L2G calibration on gold-standard.
- $y^{\text{SC}}_i \in [-1, 1]$: Cell-type cognition-enrichment score for target $i$ from cellxgene-census. SE $\sigma^{\text{SC}}_i$ from cell-resampling bootstrap.

Optional fourth channel: $y^{\text{Lit}}_i$, Lit-OTAR target–cognition association strength (strongly down-weighted by skeptical prior).

Latent target *cognition-relevance*: $\theta_i \in \mathbb{R}$ on the logit / atanh scale. Pipeline-consumable target weight: $w_i = \sigma(\theta_i) \in (0,1)$ or a transformed version mapped to a multiplicative bias.

#### B.2 Likelihood

Fisher-transform each observed score (atanh for AHBA/SC, logit for L2G/Lit) to $\tilde y^s_i$ on the real line:

$$
\tilde y^s_i \;\sim\; \mathcal{N}\!\left(\alpha_s + \beta_s\, \theta_i,\;\; \tau_s^{-1} + (\sigma^s_i)^2\right)
$$

- $\alpha_s$: source-specific intercept.
- $\beta_s > 0$: source informativeness slope.
- $\tau_s$: source precision — the data-driven reliability weight.
- $(\sigma^s_i)^2$: heteroskedastic per-target measurement noise.

#### B.3 Priors

$$
\theta_i \sim \mathcal{N}(0,1)\quad\text{(scale-fixed)};\qquad
\alpha_s \sim \mathcal{N}(0,0.5^2);\qquad
\beta_s \sim \text{HalfNormal}(1.0);\qquad
\tau_s \sim \text{Gamma}(2,2)
$$
$$
\beta_{\text{Lit}} \sim \text{HalfNormal}(0.3) \quad\text{(skeptical-prior downweight)}
$$

Soft sum-to-zero on $\{\alpha_s\}$. A small reference-target set (BDNF, COMT, ACHE, DRD2, GRIN2B, CHRNA7) gets weakly informative prior $\theta_i \sim \mathcal{N}(0.5, 0.3^2)$ to break scale + sign degeneracy. Sensitivity analysis runs with and without these anchors.

#### B.4 Posterior inference

- **Primary: NUTS in PyMC 5.x** via numpyro JAX backend. T=210, ~230 free dimensions — well within NUTS range. 4 chains × 2,000 warmup × 2,000 draws → 5–15 min on RTX 5070.
- **Secondary: ADVI mean-field** as fast regression test (~ 30 s).
- **Tertiary: Pathfinder VI** (pymc-extras) for NUTS initialization.

Diagnostics: $\hat R < 1.01$, ESS > 400 for every $\theta_i$, zero divergences at `target_accept=0.95`.

#### B.5 Identifiability

| Parameter | Identifiability | Failure mode | Mitigation |
|---|---|---|---|
| $\theta_i$ | Identifiable up to global sign + scale | Sign flip on cold start | Reference anchors, $\beta_s > 0$ |
| $\beta_s$ | Identifiable when ≥ 2 sources informative | Single-source pathology | HalfNormal constraint, minimum-source rule |
| $\alpha_s$ | Weakly identifiable | Trades off with mean of $\theta_i$ | Sum-to-zero, fixed-scale $\theta$ prior |
| $\tau_s$ | Identifiable when $T \gg S$ | Underestimated for low-variance sources | Gamma(2,2) prior away from 0 |

Confounded combo: $\alpha_s + \beta_s \bar\theta$ — resolved by fixing $E[\theta] = 0$ and soft sum-to-zero on $\alpha$.

#### B.6 Sensitivity analysis

Sweep:
- $\theta$ prior: $\mathcal{N}(0, 0.5)$, $\mathcal{N}(0,1)$, $\mathcal{N}(0,2)$, Student-T(3,0,1)
- $\beta$ prior: HN(0.5/1.0/2.0), Exponential(1)
- $\tau$ prior: Gamma(2,2), Gamma(0.5,0.5), HalfCauchy(1)
- Lit weight: HN(0) ablated, HN(0.3), HN(1.0)
- Reference anchors: with / without

Report 90 % posterior CI bands across all combinations. Sign-flipping targets are flagged "low robustness."

#### B.7 Disagreement-as-signal (the methodologically novel piece)

Per-source posterior of latent given only that source:
$$
p_s(\theta_i \mid y^s_i) \;\propto\; \mathcal{N}\!\left(\theta_i \;\middle|\; \frac{\tilde y^s_i - \hat\alpha_s}{\hat\beta_s},\; \frac{\hat\tau_s^{-1} + (\sigma^s_i)^2}{\hat\beta_s^2}\right)\cdot p(\theta_i)
$$
Jensen-Shannon disagreement:
$$
D_i \;=\; H\!\left(\tfrac{1}{3}\sum_s p_s\right) - \tfrac{1}{3}\sum_s H(p_s) \;\in [0, \log 3]
$$
normalized to $[0,1]$. **Key design choice: $D_i$ is not used to penalize $w_i$.** It becomes its own diagnostic axis:

- **Low $D_i$, high $\bar\theta_i$**: high-confidence positive (BDNF — three-source convergence).
- **Low $D_i$, low $\bar\theta_i$**: high-confidence reject.
- **High $D_i$, any $\bar\theta_i$**: high information value — exactly the targets where downstream Boltz-2 + selectivity work has maximum expected information gain.

Complementary leave-one-source-out KL: $\text{KL}\big(p_{\text{full}}(\theta_i) \,\|\, p_{-s}(\theta_i)\big)$ — quantifies which source is driving disagreement.

In Phase A.7, Pareto ranking optimizes over (potency proxy, selectivity, liability penalty, $D_i$ as info-value axis, $\bar\theta_i$).

#### B.8 Why full Bayesian *supersedes* the simpler alternatives

The bounded multiplicative bias $w_i \in [0.5, 1.5]$ scheme loses three things: (1) no credible intervals for downstream §7.11 propagation; (2) collapses three sources to one number — no disagreement-as-signal; (3) hyperparameter sensitivity is invisible. The Bayesian model degenerates gracefully: $\tau_s \to \infty, \beta_s = 1$ recovers an arithmetic-mean point estimate, and a sigmoid scaled to $[0.5, 1.5]$ gives back the simpler form. **Use the Bayesian model as primary; expose the bounded-bias projection as a compatibility output for downstream consumers that can't ingest distributions.**

Engineering metaphor: the simple multiplicative bias is a fixed-resistor voltage divider; the Bayesian model is an op-amp with feedback — same DC behavior in the limit, but the latter is also instrumented and self-correcting. You don't ship a production amplifier without feedback.

---

### C) AHBA Spatial-Correlation Deep Dive

#### C.1 Pipeline

1. `abagen.get_expression_data(...)` per §A.2 → $E \in \mathbb{R}^{68 \times G}$, $G \approx 15{,}000$ stable genes.
2. Cognition reference maps:
   - **Primary**: Moodie 2024 g-cortical map (or recompute via enigma-toolbox DK68).
   - **Secondary**: Neurosynth "cognition" composite, projected to DK68 via nimare + neuromaps.
   - **Tertiary**: Hansen 2021 (*Nat Hum Behav*) cognitive-gradient PC1.
3. For each target $g$, Pearson $r_g = \text{corr}(E[:, g], \text{cog})$.
4. BrainSMASH null: 10,000 surrogates of cog map; two-sided $p_g$; `knn = n_parcels`.
5. Spin-test secondary: 10,000 spherical rotations via `netneurotools.stats.gen_spinsamples`.
6. FDR-BH at $q < 0.05$ **over the ~210-target panel only**, not all 15,000 genes.

#### C.2 Worked examples (illustrative — must be regenerated against real data)

| Target | AHBA $r$ vs g-map | BrainSMASH p | Interpretation |
|---|---|---|---|
| **BDNF** | +0.42 (prefrontal-biased) | ~ 0.002 | Strong positive — canonical |
| **HTR2A (5-HT2A)** | +0.55 (heavy frontal) | ~ 0.001 | Strong AHBA signal |
| **GRIN2B (NMDA NR2B)** | +0.28 | ~ 0.02 | Moderate positive |
| **CHRNA7 (α7 nAChR)** | +0.18 (hippocampus, thalamus) | ~ 0.10 | Cortical signal weak — *expected* given α7's deep-brain distribution |
| **ACHE** | +0.05 | ~ 0.40 | Null — substrate-level target; expression-prior misleading |
| **APOE** | -0.10 (glial) | ~ 0.30 | Cell-type confound; AHBA bulk smears it |

#### C.3 Statistical rigor

- BH-FDR over the target panel.
- BrainSMASH primary, spin test secondary.
- Leave-one-donor-out (6 donors): report posterior shift.
- Bidirectional hemisphere mirroring at abagen stage.
- Condition on Moodie C1, C2 PCs before per-target residual associations.

---

### D) OT Genetics L2G Integration

#### D.1 Pipeline

1. Pull credible sets from OT Genetics GraphQL `/v4` for: Davies 2018 (`GCST006269`), Hill 2019 (`GCST006716`), Sniekers 2017, Savage 2018 (`GCST006250`), UK Biobank fluid intelligence + cognitive battery, and Lee 2018 educational attainment as proxy (with caveat).
2. For each credible set, pull L2G scores for genes within ±500 kb.
3. Aggregate: per target $i$, $y^{\text{L2G}}_i = \max_{\text{study, cs}} L2G(\text{cs}, i)$; retain credible-set provenance.
4. Thresholds: ignore L2G < 0.05; 0.05–0.2 "in the room"; 0.2–0.5 suggestive; ≥ 0.5 high-confidence causal.
5. Carry forward credible-set posterior probability as $\sigma^{\text{L2G}}_i$.

#### D.2 Integration with expression evidence

L2G says "this gene causally moves cognition phenotypes in the population"; AHBA says "this gene aligns spatially with cognition map"; SC says "this gene is concentrated in cognition-salient cell types." All three converging → high-confidence target. Any one or two converging → moderate. None converging → low.

---

### E) Transcriptomic Prior from cellxgene-census

#### E.1 Cell-type cognition-enrichment

"Cognition-salient" set $\mathcal{C}^{*}$ via expert + literature curation: L2/3 IT pyramidal, L5 ET pyramidal, L6 CT, hippocampal CA1/CA3/DG, basal forebrain cholinergic (CHAT+), LC noradrenergic (DBH+), raphe serotonergic (TPH2+), VTA DA (TH+), cortical PV+/SST+/VIP+ interneurons.

For target $g$ and cell type $c$:
$$
\mu_{g,c} = \mathbb{E}[\log(1 + \text{CP10K}_g) \mid c, \text{brain}, \text{primary}]
$$
$$
e_{g,c} = (\mu_{g,c} - \bar\mu_g)/s_g\quad\text{(z-score across brain cell types)}
$$
$$
y^{\text{SC}}_g = \sum_{c \in \mathcal{C}^{*}} w_c \cdot e_{g,c}
$$
$w_c$ uniform = 1.0 default; sensitivity sweep varies them.

#### E.2 Integration with bulk AHBA

Both estimate "where is gene X expressed in the brain" at different resolution. Agreement for canonical markers expected; disagreement for any target informative (likely bulk cell-type confounding). Optional CIBERSORTx / MuSiC deconvolution deferred to V2.

#### E.3 WSL2 / RTX 5070 deployment

Pin: `tiledbsoma>=1.15.3`, `cellxgene_census>=1.16`, `pyarrow>=15`, `scanpy>=1.10`. 32 GB system RAM comfortable for brain slice; iterate via PyArrow tables in 100 k-cell chunks. **Do not** materialize the full Siletti slice as a single AnnData. RTX 5070 12 GB VRAM is needed only for PyMC NUTS; cellxgene aggregation is CPU/RAM-bound.

---

### F) Target Scope Expansion (44 → ~210)

#### F.1 Inclusion (any of):

1. **GWAS-prioritized**: L2G ≥ 0.2 in Davies 2018 / Hill 2019 / Sniekers 2017 / Savage 2018 credible set.
2. **MAGMA-prioritized**: gene-based p < 2.7e-6 (Bonferroni for ~18,000 genes) over the same GWASes.
3. **AHBA-prioritized**: BrainSMASH-corrected |r| > 0.3 with FDR q < 0.05.
4. **Single-cell prioritized**: cell-type enrichment z > 2 in any cognition-salient cluster.
5. **Literature-prioritized**: Lit-OTAR target–cognition strength ≥ 0.5 (top ~5 percentile).
6. **Existing §8.0b panel member** (kept for panel continuity).

#### F.2 Exclusion

- Not expressed in brain (GTEx brain TPM < 1 AND cellxgene brain max log-expr < 1)
- Druggability tier 5 + no published modulator chemotype
- Pure housekeeping (ribosomal, basic glycolysis)
- Predominantly developmental (FOXP2 borderline — kept if adult cognition evidence exists)

#### F.3 Expected composition (~210)

| Class | Approx N |
|---|---|
| GPCR (5-HT, DA, ACh-muscarinic, mGlu, adrenergic, histamine, orexin) | ~ 50 |
| Ion channels (NMDA, AMPA, Kv, Nav, Cav) | ~ 35 |
| Ligand-gated (nAChR α7/α4β2/α3β4, GABA-A α1-6) | ~ 15 |
| Transporters (DAT, NET, SERT, VMAT2, ChT1, EAAT) | ~ 12 |
| Enzymes (AChE, MAO-A/B, COMT, PDE4/9/10, GSK3β, HDAC2/6, BACE1) | ~ 30 |
| Kinases / signalling (CAMK2, ERK, mTOR, AKT) | ~ 20 |
| Synaptic scaffold / plasticity (PSD-95, SHANK3, ARC, CAMKK2, BDNF, NTRK2) | ~ 25 |
| Neurogenesis / myelination (per Hill 2019 enrichment) | ~ 15 |
| Other (TFs, growth factors, neuroimmune) | ~ 10 |
| **Total** | **~ 212** |

This subsumes §7.3 — the Cluster D target panel *is* the V3 working set.

---

### G) Behavioural Validation Gates

#### G.1 Gate 1 — Roberts 2020 SMD ceiling (HARD GATE)

Roberts CA, Jones A, Sumnall H, Gage SH, Montgomery C (*Eur Neuropsychopharmacol* 2020;38:40–62): in k = 47 trials, modafinil pooled SMD = 0.12; methylphenidate SMD = 0.21 (recall 0.43, sustained attention 0.42, inhibitory control 0.27); d-amphetamine null.

**Pass**: no target's predicted modulator effect-size posterior exceeds Hedges' g = 0.5 at 90 % credible upper bound. Targets that do are flagged "regime-violating" → manual review.

#### G.2 Gate 2 — Clinical endpoint correlation

For targets with a Phase 2+ modulator and a published cognition endpoint (MMSE / ADAS-Cog / n-back / digit symbol / CANTAB / Cambridge Brain Sciences), per-target $\bar\theta_i$ must positively correlate (Spearman > 0.3, p < 0.05) with the meta-analytic SMD. Reference compounds: donepezil, memantine, rivastigmine, modafinil, methylphenidate, atomoxetine, varenicline, encenicline, pridopidine, brexpiprazole; ~ 15 targets.

**Pass**: Spearman ρ ≥ 0.3 with 90 % bootstrap CI excluding 0.

#### G.3 Gate 3 — GWAS phenotype validation (out-of-sample)

Hold out Cognitive Aging Consortium and ABCD Study cognition battery. Refit L2G channel using only training GWAS; predict on held-out; **AUROC > 0.7**.

#### G.4 Gate 4 — Cross-source validation

Three-fold leave-one-source-out: refit ablating each source; held-out source must positively correlate (Spearman ρ > 0.2) with the ablated-model posterior. **All three folds must pass**; otherwise multi-source assumption is violated and that source is downweighted in production.

---

### H) Disagreement-as-Signal Worked Examples (illustrative)

Eight targets, illustrative numbers on Fisher-transformed scale.

| Target | $y^{\text{AHBA}}$ | $y^{\text{L2G}}$ | $y^{\text{SC}}$ | $\bar\theta_i$ [90 % HDI] | $D_i$ | Verdict |
|---|---|---|---|---|---|---|
| **BDNF** | +0.65 | +0.55 (Hill 2019, L2G 0.78) | +0.70 | +0.78 [0.62, 0.93] | 0.08 | Three-source agreement — gold-standard positive control |
| **HTR2A (5-HT2A)** | +0.55 | +0.10 | +0.30 (L5 ET) | +0.35 [0.05, 0.65] | 0.62 | High-disagreement positive — AHBA bullish, GWAS silent. **This is exactly what the framework was designed to surface. High Boltz-2 priority.** |
| **CHRNA7 (α7 nAChR)** | +0.18 (subcortical-biased) | +0.05 | +0.45 (CA1-enriched) | +0.22 [-0.05, +0.50] | 0.48 | Partial agreement, SC drives. Cortical AHBA undersamples α7's niche. **Trust SC over AHBA.** |
| **COMT** | +0.30 | +0.65 (Davies 2018, L2G 0.71 at rs4680) | +0.20 | +0.45 [0.20, 0.70] | 0.38 | GWAS-dominant; classic pharmacogenetics |
| **GRIN2B (NMDA NR2B)** | +0.28 | +0.40 | +0.55 | +0.45 [0.22, 0.68] | 0.15 | Three-source moderate agreement — solid prior |
| **CACNA1C** | +0.05 | +0.55 (psychiatric pleiotropy) | +0.20 | +0.28 [+0.05, +0.50] | 0.52 | GWAS-dominant but pleiotropy → flag for §7.4 selectivity follow-up |
| **APOE** | -0.10 (glial) | +0.30 (AD-mediated) | -0.20 | +0.05 [-0.20, +0.25] | 0.55 | High-disagreement near-zero — AD-mediated, not cognition-mediated. Down-prioritize for healthy enhancement |
| **ACHE** | +0.05 | +0.10 | +0.10 | +0.10 [-0.10, +0.30] | 0.10 | Three-source agreement low signal. **Framework limitation: expression-level priors are blind to substrate-mediated cognition effects.** Manual override. |

The ACHE row is the most important honest critique: **the framework is biased toward targets where modulation magnitude tracks expression magnitude**. For substrate-degrading enzymes (ACHE, MAO, COMT-degradation), this assumption is wrong. Tag every enzyme target with a "substrate-mediated" flag.

---

### I) Implementation Spec (~30 %)

#### I.1 PyMC 5.x skeleton

```python
import pymc as pm
import pytensor.tensor as pt
import numpy as np

def fit_cluster_d_prior(y_obs, sigma_obs, reference_idx,
                        reference_mean=0.5, reference_sd=0.3):
    """
    y_obs, sigma_obs: (S=4, T) Fisher/logit-transformed observations.
    Source order: [AHBA, L2G, SC, Lit].
    """
    S, T = y_obs.shape
    with pm.Model() as model:
        theta = pm.Normal("theta", mu=0.0, sigma=1.0, shape=T)
        pm.Normal("ref_anchor", mu=reference_mean, sigma=reference_sd,
                  observed=theta[reference_idx])
        alpha = pm.Normal("alpha", mu=0.0, sigma=0.5, shape=S)
        beta_scale = pt.as_tensor_variable([1.0, 1.0, 1.0, 0.3])  # skeptical on Lit
        beta = pm.HalfNormal("beta", sigma=beta_scale, shape=S)
        tau  = pm.Gamma("tau", alpha=2.0, beta=2.0, shape=S)
        pm.Normal("alpha_sum", mu=pt.sum(alpha), sigma=0.05, observed=0.0)
        mu_s = alpha[:, None] + beta[:, None] * theta[None, :]
        sigma_s2 = (1.0 / tau)[:, None] + sigma_obs ** 2
        pm.Normal("y", mu=mu_s, sigma=pt.sqrt(sigma_s2), observed=y_obs)
        idata = pm.sample(2000, tune=2000, chains=4, cores=4,
                          target_accept=0.95, nuts_sampler="numpyro",
                          random_seed=42)
    return model, idata
```

#### I.2 Schema additions to V3 target table

```yaml
target_id: ENSG00000102468
gene_symbol: HTR2A
cluster_d:
  y_ahba: 0.55;   sigma_ahba: 0.08;   brainsmash_p: 0.001
  y_l2g_max: 0.10; l2g_study_id: GCST006269; sigma_l2g: 0.05
  y_sc: 0.30;     sigma_sc: 0.10
  cell_type_enrichment: { L5_ET: 2.1, L2_3_IT: 1.4, CA1: 0.8 }
  y_lit: 0.42
  theta_posterior: { mean: 0.35, sd: 0.18, hdi_5: 0.05, hdi_95: 0.65 }
  disagreement_jsd: 0.62
  substrate_mediated_flag: false
  regime_violation_flag: false
```

#### I.3 Plug into Phase A.7

$$
w^{\text{final}}_i = w^{\text{cal}}_i \cdot \sigma(\theta_i^{\text{post}}) \cdot \big(1 + \gamma \cdot \tfrac{1}{1 + \text{HDI\_width}(\theta_i)}\big)
$$
Full posterior preserved for §7.11 isotonic + hierarchical Bayesian calibration.

#### I.4 Validation dashboard (Plotly Dash + FastAPI)

5 panels: (1) caterpillar plot of $\theta_i$ sorted, colored by $D_i$; (2) ternary plot of per-source agreement; (3) Shapley-style decomposition per target; (4) sensitivity-sweep KS-distance matrix; (5) validation-gate scoreboard.

#### I.5 Runtime budget (RTX 5070 + 32 GB RAM, WSL2)

| Stage | Time |
|---|---|
| abagen full AHBA processing (cached) | 8–15 min |
| BrainSMASH 10,000 surrogates × 210 targets | 25–40 min CPU |
| OT Genetics L2G pull (cached) | 5–10 min |
| cellxgene census brain slice + aggregation | 30–60 min (network-bound; cache locally) |
| PyMC NUTS (numpyro, 4 chains × 4 k iter) | 6–15 min |
| Sensitivity sweep (~ 20 runs) | 2–5 hours |
| Cold start total | ≈ 4–8 hours |
| Warm rerun | 15–30 min |

---

### J) Publishability Analysis

#### J.1 Novelty audit

- **Open Targets Platform integrated score**: harmonic sum across heterogeneous evidence; **no Bayesian uncertainty propagation, no disagreement-as-signal axis.**
- **Driessens 2023 (*Nat Commun*; PMC10345092)** + **Moodie 2024 (*HBM*)**: AHBA × cognition only; no GWAS or single-cell integration; not drug-prioritization-oriented.
- **Hansen 2022 (*Nat Neurosci* 25:1569–1581)**: neurotransmitter receptor PET × Neurosynth × ENIGMA; no GWAS or single-cell integration.
- **Driver-style multi-omic fusion (NetICS, DiffNetFDR)**: cancer-focused; not cognition; typically no formal uncertainty propagation.

**Three claimed novel contributions:**

1. *First* multi-source neurobiological prior for cognition target prioritization combining AHBA, OT Genetics L2G, and cellxgene-census single-cell, with formal Bayesian uncertainty propagation and credible intervals.
2. *Disagreement-as-signal* axis (Jensen-Shannon divergence over source-conditioned posteriors) as a Pareto-optimization dimension — generalizable methodology not specific to cognition.
3. Behavioural-anchor validation gate (Roberts 2020 empirical SMD ceiling) constraining pipeline outputs to the empirically plausible regime — an honest answer to "computational target prioritization predicts unrealistic effect sizes."

#### J.2 Target venues

| Venue | Fit | Rating |
|---|---|---|
| **Cell Reports Methods** | Computational methods, biomedical — strongest fit | **A+** |
| **Bioinformatics** (OUP) | Methodological + reproducible code; OT/L2G publishes here | **A** |
| **PLOS Computational Biology** | Methods + biology integration | **A** |
| **NeurIPS / ICML "AI for Science"** | If framed as Bayesian-ML novelty | **B+** |
| **Nature Methods** | Possible but high bar; needs strong benchmark | **B** (stretch) |
| **Briefings in Bioinformatics** | Reviews + methods | **A-** |

#### J.3 Anticipated reviewer objections

1. *"Conflates expression-level priors with mechanistic effect-size priors."* → Yes, honestly tagged. Substrate-mediated flag + Roberts 2020 gate. Combination of L2G + AHBA + SC is less expression-naive than any single source.
2. *"AHBA has only 6 donors."* → Leave-one-donor-out sensitivity reported; BrainSMASH null partially compensates.
3. *"No actual cognitive-enhancement RCT validation."* → Gate 2 is exactly this — Spearman ρ between predicted $\bar\theta$ and meta-analytic SMD on the ~ 15 well-studied targets. We report this ρ as a headline.
4. *"L2G is European-ancestry-biased."* → Caveat. Out-of-ancestry validation requires data we don't have; flagged limitation.
5. *"Disagreement-as-signal is just UQ under a different name."* → Distinct: JSD over source-conditioned posteriors captures **structural disagreement** (sources point different directions) vs measurement uncertainty (a source is noisy). The covariance matters.

#### J.4 Required validation experiments

- All four Gates pass on actual data.
- Comparison against three baselines: (a) AHBA-only single-source, (b) OT harmonic-sum, (c) multiplicative bias $[0.5, 1.5]$. Show Bayesian + disagreement Pareto-dominates.
- Held-out cognition-GWAS validation (ABCD, CAC).
- Public reproducible code release (GitHub + Zenodo DOI).

---

### K) Identifiability and Sensitivity Analysis

Per §B.5 / §B.6. Publication deliverables:

1. Table: posterior 90 % CI on each $\theta_i$ across 20 hyperprior settings.
2. Heatmap: pairwise KS distances between sensitivity-run posteriors.
3. Sign-stability table: fraction of sensitivity runs each target keeps the same sign.
4. Ablation: drop each source in turn, retrain, compare.
5. Reference-anchor ablation: drop the 6 anchors, check sign-stability.

Targets failing sign-stability are flagged "low robustness" and not surfaced to Phase A.7.

---

### L) Limitations and Honest Framing

- **Roberts 2020 ceiling.** Real cognition-enhancement SMDs in healthy adults are typically < 0.3, capped near 0.5 in best domains. The framework predicts *relative* prioritization, not absolute effect sizes. Do not interpret $\theta_i$ as a predicted SMD.
- **AHBA donor count = 6**; right hemisphere covered in only 2 donors. Cortical-asymmetric phenomena essentially untestable.
- **OT Genetics L2G** European-ancestry-dominated; gold-standard under-represents brain disorders; calibration approximate.
- **cellxgene-census brain coverage** heterogeneous: cortex well-sampled, hippocampus moderately, brainstem nuclei (LC, raphe) poorly. Enrichment scores for those regions have high noise.
- **Causal interpretation.** $\theta_i$ is a prior probability of cognitive relevance, not a causal effect estimate. L2G is calibrated for "causal at this locus," not "perturbing this gene causes enhancement."
- **g vs sub-domain** — Moodie 2024 g is psychometric; Roberts 2020 outcomes are sub-domain. Correlated but not identical; framework treats as approximately interchangeable, which is a simplification.
- **Cognition ≠ disease.** Most cognition GWAS is partially contaminated by educational attainment / SES proxies; we down-weight EA evidence but cannot fully decontaminate.
- **AD-heavy Lit-OTAR**, compensated by skeptical prior $\beta_{\text{Lit}} \sim \text{HN}(0.3)$.
- **Substrate-mediated mechanism blindness.** ACHE, MAO, etc. — modulator effect-size decoupled from expression. Flagged but not corrected.
- **Conditional-independence assumption** in the likelihood is wrong in detail (AHBA and cellxgene both measure expression; L2G uses AHBA-adjacent eQTL data). Deferred to V2; explicit limitation.

---

## Recommendations

**Stage 1 (Weeks 1–3): Foundation.** abagen + BrainSMASH plumbing. OT Genetics L2G pulls. cellxgene-census brain slice cached locally. Define $\mathcal{C}^{*}$ in writing.
- *Decision benchmark*: BDNF positive with three-source agreement. If not, plumbing is broken.

**Stage 2 (Weeks 4–6): Target expansion.** Generate ~ 210-target panel per §F. Validate existing 44-target §8.0b is a strict subset.
- *Decision benchmark*: ≥ 80 % of new targets have a published modulator chemotype in ChEMBL.

**Stage 3 (Weeks 7–9): Bayesian model.** Implement PyMC per §I.1. Run sensitivity sweep.
- *Decision benchmark*: $\hat R < 1.01$, ESS > 400, zero divergences, sign-stability > 90 % across sweep.

**Stage 4 (Weeks 10–12): Validation gates.** All four Gates from §G.
- *Decision benchmark*: All four pass. If Gate 2 fails, audit substrate-mediated tagging.

**Stage 5 (Weeks 13–16): §7.11 integration + publication.** Plug into calibration. Re-run downstream Pareto. Submit to Cell Reports Methods (primary) or Bioinformatics (secondary).

**Threshold-driven rules:**
- Sign-stability < 80 % → downgrade from "primary prior" to "secondary diagnostic."
- Gate 2 Spearman ρ < 0.2 → do not publish; core empirical claim failed.
- Median $D_i > 0.7$ → sources too inconsistent; fall back to single-source priors and publish disagreement as the main finding.
- 90 % HDI width > 0.6 for > 50 % of targets → framework inconclusive; expand reference-anchor set.

---

## Caveats

- The "Mansuri 2024 41-gene cognition map" reference is most likely a misattribution; canonical paper is **Moodie et al. 2024, *Hum Brain Mapp* 45(4):e26641**. Correct internal references before submission. A flat reject from a domain reviewer is the real risk if uncorrected.
- All numerical examples in §C.2 and §H are *illustrative*. Actual posterior values must be regenerated against real data. Do not cite them as findings.
- The Bayesian model's conditional-independence-given-$\theta$ assumption is wrong in detail. A more rigorous multivariate likelihood with cross-source covariance is V2 scope; flag in published limitations.
- Cognition GWAS power saturates near where Roberts 2020 saturates: intelligence PGS R² is currently 4.4–4.8 % (Davies 2018: R²=4.37 %; Savage 2018 up to 4.81 %), and pharmacological SMDs in healthy adults cap near 0.5 — both suggest we are working near the information-theoretic ceiling of current data. Cluster D will not produce miracles; it will reduce noise, improve calibration, and surface a small number of high-information-value targets.
- Single-GPU RTX 5070 is comfortable for V1 (T=210). Extending to gene-level (T ≈ 15,000) latents requires a Bayesian neural-network surrogate or sparse approximation — out of V1 scope.