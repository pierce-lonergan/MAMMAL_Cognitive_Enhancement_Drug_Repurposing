# V6 Architecture & Phased Implementation Plan

**Status**: live source-of-truth for V6 design. Companion to
`design/V4_STATUS_AND_FORWARD_PLAN.md` §13 (V5/V6 Path Forward). Concrete
implementation roadmap for the two V6 priorities:

- **V6.A — Multi Head DTI ensemble** (~12 weeks) per
  `research/4-tier/Multi Head DTI.md`
- **V6.B — Bayesian Cluster D neurobiological prior** (~16 weeks) per
  `research/4-tier/Multi-Source Neurobiological Prior for Cognition Target Prioritization.md`

The V6 scaffolds for both are already shipped in `src/mammal_repurposing/`
(diagnostics/per_head_bias.py, fusion/bayesian_router.py,
cluster_d/bayesian_prior.py, cluster_d/data_fetchers.py). The implementation
work below operationalises them with real heads + real data.

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

### Phased plan

| Phase | Weeks | Deliverable | Validation |
|---|---|---|---|
| **V6.A.1 Heads installed** | 1-3 | MMAtt-DTA (pip + Zenodo weights ~2 GB; adapter shipped at `cluster_a/mmatt_dta_adapter.py`); PSICHIC; BALM (ESM-2 + ChemBERTa-2, ~3 GB) | Each head produces per-target ρ vs ChEMBL pchembl≥8 truth |
| **V6.A.2 Bias decomposition** | 4-5 | Wire `diagnostics/per_head_bias.py` to compute PC_k, SN_k, OOD_k, CT_k per (head, target). Bonett-Wright CIs | 5-head × 22-target trust matrix T(t, k) ∈ [0.02, 0.7] with row entropy logged |
| **V6.A.3 Bayesian router** | 6-7 | `fusion/bayesian_router.py` is shipped; activate by passing real T(t, k); add identifiability diagnostic report | Per-target router weights logged; identifiability theorem confirms n*=720 >> v4 n=7-26 (priors, not posteriors) |
| **V6.A.4 Calibrated uncertainty** | 8 | Per-head Venn-ABERS (Mervin 2020 J Chem Inf Model 60:4546). Cross-head correlation matrix Σ_kk' from 133-tuple calibration set. Replace router's Gaussian-CI with VA MC propagation | CI width inflation factor √(1+(K-1)·r̄) ≈ 1.41 confirmed |
| **V6.A.5 Disagreement facet** | 9-10 | Extend `35_v3_disagreement_signal.py` to multi-head: pairwise Kendall-τ + rank-distance + facet-tag {novel_scaffold / activity_cliff / ood / noise} | `reports/disagreement_axis_v1.md` per-compound bucket |
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

## V6.A × V6.B Composition (the joint posterior)

Per V4 §13.3, the joint posterior over (compound, target) pairs is:

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

## Resource allocation decision tree

If a **single research-engineer-month** is available between now and V5 launch:
1. **MMAtt-DTA head** (V6.A.1 first slice, 2-3 days) — single biggest disambiguation between "Tanimoto is the right baseline" and "modern DTI heads can beat it." Settles the empirical question.
2. **Pose-saving Boltz wrapper** (V4 §7.17 — code shipped; ~6-10h pose-only Boltz re-run on WSL2 GPU). Operationalises §7.5 on the live grid and unblocks §8.13 pocket-conditioned liability gating.
3. **Cluster C TxGNN run** (V4 last Tier-1 item; ~1 day setup including DGL wheel pin + PrimeKG download + 1hr run). Adds the 5th cluster to RRF.

If **2-3 months** are available — ship the full V6.A.1-A.4 Multi Head DTI core (heads + bias decomposition + Bayesian router + calibration). Disagreement facet (V6.A.5) and publication (V6.A.6) come naturally.

If **6+ months** are available — ship V6.A in full, then begin V6.B Cluster D. The V6.B Bayesian model requires V6.A's calibrated uncertainty propagation as input.

If **no engineer time** is available — the current `reports/wet_lab_shortlist_v6_full.md` is the production deliverable. 43 PASS compounds with all V4 + V5 gates flowing through (calibrated MAMMAL + Z-norm + Tanimoto + ADMET + MoA + §8.0b-zn liability + Pareto + scaffold-AL + nootropic-similarity + CTgov IP). The contribution is honest and defensible as-is.

---

## V6 timeline summary

| Track | Wks | Effort | Dependencies | Output |
|---|---|---|---|---|
| V6.A: Multi Head DTI ensemble | 12 | 5 heads + bias decomposition + Bayesian router + eMOSAIC OOD + Venn-ABERS + disagreement facet + validation | None (heads are pip-installable; MMAtt-DTA needs Zenodo download) | `fusion/bayesian_router.py` (shipped scaffold) + `diagnostics/per_head_bias.py` (shipped scaffold) + `diagnostics/ood_emosaic.py` (pending) + `diagnostics/disagreement_axis.py` (pending; §8.15 is partial) + paper draft (J Cheminform / Nat Mach Intell) |
| V6.B: Bayesian Cluster D prior | 16 | abagen + BrainSMASH + OT L2G + cellxgene + PyMC NUTS + 4-gate validation + §7.11 integration + paper | abagen / BrainSMASH / PyMC installs; cellxgene-census brain slice (~10 GB local cache) | `cluster_d/bayesian_prior.py` (shipped scaffold) + `cluster_d/data_fetchers.py` (shipped scaffold) + 5 validation reports + paper draft (Cell Reports Methods / Bioinformatics) |
| Composition | 4 | Joint-posterior plumbing + V7 wet-lab shortlist re-render | V6.A + V6.B both shipped | `reports/wet_lab_shortlist_v7_joint.md` |

**Total V6 (A + B + Composition)**: ~32 weeks (~8 months) of focused engineering.

Two distinct papers, two distinct validation regimes, two distinct
publication venues. The shortlist that lands at the end is the first
cognition-enhancement candidate set in the literature with:
- formal credible intervals on every compound's rank,
- behavioural validation gate (Roberts 2020 SMD ceiling),
- multi-head ensemble with disagreement-as-signal discovery axis,
- mechanism + liability gating,
- AND a per-compound provenance trail back to documented signal sources with known failure modes.

That candidate set, not the next nootropic, is the contribution.

---

## V6 scaffold inventory (already shipped)

Pre-V6 plumbing already in `src/mammal_repurposing/`:

| Module | Lines | Purpose |
|---|---|---|
| `cluster_a/mmatt_dta_adapter.py` | ~200 | MMAtt-DTA adapter with 22-target superfamily map |
| `diagnostics/per_head_bias.py` | ~180 | PC/SN/OOD/CT signature computation + trust matrix builder |
| `fusion/bayesian_router.py` | ~210 | Per-target router + OOD + confidence gates + identifiability diag |
| `cluster_d/bayesian_prior.py` | ~250 | Stage-0 stub + PyMC NUTS full Bayesian model (Neelon-Dunson) |
| `cluster_d/data_fetchers.py` | ~150 | AHBA / OT Genetics / cellxgene adapters with availability probes |
| `analysis/brain_region.py` | ~140 | Static 22-target brain-bias map (V5 fallback / V6 preview) |

**LOC total**: ~1,130 lines of scaffold + tests across the V6 architectural footprint.

When the heads + data sources arrive, the scaffold activates via:
- Set `MMATT_DTA_ROOT` env var → MMAtt-DTA adapter active
- `pip install pymc numpyro` → Cluster D Bayesian path active
- `pip install abagen brainsmash` → AHBA real-mode active
- Set `CRYPTOBENCH_HOME` / `POCKETMINER_HOME` → detector ensemble Sprint 2 active

---

## How this doc was written

V6 plan synthesized from:
- `research/4-tier/Multi Head DTI.md` (~60 KB) — pre-committed Multi Head DTI spec
- `research/4-tier/Multi-Source Neurobiological Prior for Cognition Target Prioritization.md` (~39 KB) — pre-committed Cluster D Bayesian model
- V4 plan §13 V5/V6 Path Forward (already integrated)
- This sprint's V6 scaffolding work (commits 458881c + 4efcbea + the Tier-3b push)

When V6 work begins, the next assistant should:
1. Re-read both research deep-dives in full
2. Activate the scaffolds via the env var / install instructions above
3. Run V6.A.1 (MMAtt-DTA) FIRST — it's the smallest-effort, highest-leverage step
4. Use this doc's phase tables as the project plan
5. Maintain `reports/hypothesis_audit_v1.md` as the standing falsifiability check

Generated 2026-05-26 alongside the V6 scaffold commit.
