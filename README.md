# MAMMAL Cognitive Enhancement Drug Repurposing

A multi-layer Bayesian pipeline for cognition-enhancement drug repurposing, built around IBM Research's [MAMMAL](https://github.com/BiomedSciAI/biomed-multi-alignment) foundation model and extended with four downstream architectural layers — V4 (calibrated multi-cluster fusion), V5 (Z-norm + Tier 2/3 sprint), V6 (multi-head DTI ensemble + Bayesian Cluster D neurobiological prior), V7 (PBPK-anchored hierarchical effect-size translation), and V8 (πphen perturbational evidence axis). Runs on a single 12 GB consumer GPU (Blackwell sm_120 RTX 5070).

> **Honest scope**: this pipeline does **not** discover "smart drugs." It enriches a candidate set so wet-lab cycles spend money on plausibility, not chemistry-lottery tickets. Roberts CA et al. (*Eur Neuropsychopharm* 2020) puts the effect-size ceiling for healthy-adult cognitive enhancement at SMD ≈ 0.21 (methylphenidate overall). The deliverable here is a calibrated, provenance-rich ranking + a publishable methodology contribution — not a miracle compound.

---

## Headline metrics (current sprint)

- **503 non-slow pytest tests pass** (1 skip intentional); **12 / 14 slow pass**
- 🎯 **Prospective repurposing shortlist** (`reports/pipeline/repurposing_shortlist_v1.md`, Gap 7 — the capstone): approved drugs ranked as mechanism-justified repurposing hypotheses per disease, by **disease-class prognostic prior × target engagement** (success-track-record classes only), with safety flags, novelty flags, and GRADE dossiers. Surfaces real, literature-grounded hypotheses — **CIAS → buspirone/tandospirone (5-HT1A) + cevimeline/pilocarpine (M1/M4)** with xanomeline correctly flagged *standard*; **FXS → roflumilast (PDE4)**; **AD → fluvoxamine/blarcamesine (σ1)**.
- ✅ **Full 31-target panel scored with the real MAMMAL DTI head** (RTX 5070): added CHRM1/CHRM4 (M1/M4) + HTR6 (5-HT6) + GRM2/3/5 + GlyT1 + HTR4. The CIAS shortlist now **surfaces the muscarinic M1/M4 mechanism** (xanomeline-KarXT's class, FDA-approved 2024) it previously could only price; AD scores HTR6 correctly demoted as a Phase III failure class. (MAMMAL runs in a Python-3.12 venv — see `docs/MAMMAL_SETUP.md`.)
- 🏆 **Retrospective clinical-outcome validation** (`reports/pipeline/retrospective_clinical_validation_v1.md`): on a leakage-audited ledger of 31 real cognition drugs, **mechanism-class track record discriminates clinical SUCCESS vs Phase III FAILURE at AUROC 1.00** (perm p = 0.0002), flagging **9 / 9 famous Phase III failures** (encenicline, idalopirdine, intepirdine, pomaglumetad, PF-04447943, SUVN-502, ABT-126, TC-5619, MK-0249) it was never told about — while **target-binding affinity (AUROC 0.12) and target genetic-relevance (0.59) sit at or below chance.** The empirical case that cognition repurposing must be class-aware, not affinity-driven.
- 🔬 **Allosteric learn-to-rank head** (`reports/pipeline/allosteric_ltr_v1.md`, Gap 4): MAMMAL's sequence-only binding is **flat within target** (std 0.01–0.05 across ligands spanning 3 log-units of affinity — it cannot rank binders). A fused head [MAMMAL ⊕ Tanimoto ⊕ Boltz ⊕ physicochemistry] trained on ChEMBL and evaluated **held-out** on the cited allosteric benchmark lifts within-target Spearman ρ from **+0.02 (MAMMAL alone) to +0.51** — recovering the ranking the foundation model cannot.
- 🩺 **Clinician evidence dossiers** (`reports/pipeline/clinician_dossiers_v1.md`, Gap 5): one-page **GRADE-style** cards — effect size + credible interval, evidence quality with explicit reasons, mechanism-class track record, off-target liability flags, provenance, failure-mode caveats — the artifact a doctor actually reads.
- 📊 **External benchmark** (`reports/pipeline/external_benchmark_v1.md`, Gap 6): on the shared held-out task, the class track record (AUROC 1.00) beats the two leakage-free target-centric paradigms — affinity (0.47) and genetics (0.59) — with the target-popularity "knowledge" baseline (0.96) shown to be a **hindsight confound** (popularity follows success).
- 🏆 **Disease-population reframe** (`reports/pipeline/disease_reframe_v1.md`): re-scoring the same grid with each disease's *own* pivotal-trial track record recovers the right winning mechanism for three diseases it was never optimised against — **Alzheimer's → cholinesterase inhibitors** (within-disease class AUROC **0.97**, p = 0.003, **100 % of the 10 historical AD failures flagged**), **schizophrenia (CIAS) → muscarinic M1/M4** (xanomeline-KarXT, FDA-approved 2024), **Fragile X → PDE4** (zatolmilast). Same machinery, three diseases, three correct mechanisms.
- ✅ **Wet-lab shortlist v11** (`reports/wet-lab/wet_lab_shortlist_v11.md`): first non-degenerate (compound × target) grid — top-25 spans **7 targets** (v10 collapsed every compound onto ACHE); positive controls land at the correct mechanism (donepezil → ACHE, methylphenidate → SLC6A3, memantine → GRIN2B); max g₉₀ = 0.39 < 0.50 (honest Roberts-2020 ceiling).
- **22 / 22 hypothesis-audit verdicts: 19 PASS / 3 DEGRADE / 0 FAIL** (`reports/pipeline/hypothesis_audit_v1.md`)
- **V6.B.3 PyMC NUTS production run**: R̂ max = 1.000, ESS min = 12,780 (4 chains × 2000 draws) on the 22-target cognition panel ✅
- **V6.A multi-head ensemble**: 4 DTI heads shipped (MAMMAL + Tanimoto + MMAtt-DTA + PSICHIC + BALM scaffold) with Venn-ABERS calibration + per-target Bayesian routing
- **V6.B Cluster D**: AHBA 20/22 cognition genes × 83 brain regions cached; PyMC NUTS posterior θ̄ per target with reference anchors (BDNF, CHRNA7, GRIN2B at θ ~ N(0.5, 0.3²))
- **V7 effect-size translation**: 9-compartment PBPK + 12-class PRISMA priors + 5 failure-mode moderators; PET-anchored to Bohnen 2005 / Volkow 1998 / Kapur 2000
- **V7.4 validation**: Gate 1 (P1–P8): **5 PASS / 1 FAIL / 2 NO_COMPOUND** against real V6.A + V6.B
- **V8 πphen perturbational axis**: 7-view MOFA+ scaffold + **chemCPA trained on real LINCS L1000** (107K signatures, Val R² = 0.46 / OOD R² = 0.33) + **hierarchical transfer on real cpg0000** (R̂ = 1.010, 0 div, 60/60 compounds T > 0.6) + 8-cell disagreement classification + I_novel novel-mechanism score

---

## Architectural layers (V4 → V5 → V6 → V7 → V8)

```
                    ┌─────────────────────────────────────────┐
                    │            INPUT LAYER                  │
                    │  22 cognition-relevant targets          │
                    │  298 hand-curated + ChEMBL-expanded     │
                    │      cognition-enhancement compounds    │
                    └────────────────┬────────────────────────┘
                                     │
        ┌────────────────────────────┼────────────────────────┐
        ▼                            ▼                        ▼
┌───────────────────┐ ┌──────────────────────┐ ┌────────────────────┐
│ V6.A Multi-Head   │ │ V6.B Cluster D       │ │ V8 πphen           │
│  DTI ensemble     │ │  Bayesian neuro      │ │  Perturbational    │
│ ┌───────────────┐ │ │  prior (NEW)         │ │  Evidence (NEW)    │
│ │ MAMMAL DTI    │ │ │ ┌──────────────────┐ │ │ ┌────────────────┐ │
│ │ Tanimoto-FP   │ │ │ │ AHBA 83 regions  │ │ │ │ LINCS L1000    │ │
│ │ MMAtt-DTA     │ │ │ │ OT Genetics L2G  │ │ │ │ JUMP-CP morph  │ │
│ │ PSICHIC       │ │ │ │ cellxgene single │ │ │ │ iPSC-MEA       │ │
│ │ BALM          │ │ │ │   cell           │ │ │ │ chemCPA impute │ │
│ └───────────────┘ │ │ │ + PyMC NUTS      │ │ │ │ MOFA+ K=30     │ │
│ + Venn-ABERS      │ │ │ → θ̄ per target   │ │ │ │ joint embed    │ │
│ + Bayesian router │ │ └──────────────────┘ │ │ └────────────────┘ │
└────────┬──────────┘ └──────────┬───────────┘ └─────────┬──────────┘
         │                       │                       │
         │                       ▼                       │
         │            ┌─────────────────────┐            │
         │            │ V7 Clinical Effect- │            │
         │            │  Size Translation   │            │
         └───────────►│ (NEW; consumes V6.A │◄───────────┘
                      │  posterior +        │
                      │  V6.B θ̄ gate)       │
                      │ ┌─────────────────┐ │
                      │ │ JAX/diffrax     │ │
                      │ │  9-compartment  │ │
                      │ │  PBPK           │ │
                      │ │ + Schmidli 2014 │ │
                      │ │  PRISMA priors  │ │
                      │ │ + 5 moderators  │ │
                      │ │ → predicted     │ │
                      │ │   Hedges' g     │ │
                      │ │   with 95% CrI  │ │
                      │ └─────────────────┘ │
                      └──────────┬──────────┘
                                 │
                                 ▼
                ┌─────────────────────────────────┐
                │  Three-factor joint posterior    │
                │  π_joint ∝ π_target · π_phen     │
                │  with Gaussian-copula correction │
                │ + 8-cell disagreement classifier │
                │ + I_novel novel-mechanism score  │
                │ + Roberts 2020 ceiling pre-filter│
                └─────────────────┬───────────────┘
                                  ▼
                ┌─────────────────────────────────┐
                │  Wet-lab shortlist v11 (grid)    │
                │  differentiated (compound×target)│
                │  + retrospective clinical valid. │
                │  (reports/wet_lab_shortlist_     │
                │   v11.md · retrospective_        │
                │   clinical_validation_v1.md)     │
                └─────────────────────────────────┘
```

### Layer-by-layer status

| Layer | Sub-modules | Status |
|---|---|---|
| **V4** (calibrated 5-cluster fusion) | MAMMAL DTI + ESM2 + Boltz-2 + ADMET + Tanimoto + isotonic per-target calibration + faceted shortlist + pocket DB | ✅ shipped |
| **V5** (Z-norm + Tier 2/3 sprint) | Calibrated MAMMAL into fusion + §8.0b-zn liability + MoA ranker + nootropic similarity + pocket-routed gating + clinical-trials cross-ref + LambdaMART + conformal + scaffold-AL + brain-region | ✅ shipped |
| **V6.A** (Multi-Head DTI) | MMAtt-DTA + PSICHIC + BALM adapter + per-head bias decomposition + Bayesian router + Venn-ABERS calibration + multi-head disagreement axis | ✅ shipped |
| **V6.B** (Bayesian Cluster D) | abagen AHBA cache + OT Genetics L2G fetcher + cellxgene preview + PyMC NUTS hierarchical model + 4-gate validation (Roberts ceiling, Spearman vs SMD, GWAS-AUROC, LOSO) | ✅ shipped; NUTS converged R̂=1.000 |
| **V6 Cluster C** | PrimeKG + TxGNN per-disease ranking API | ✅ shipped (rewrite) |
| **V7** (Effect-Size Translation) | 9-compartment PBPK + 12-class PRISMA priors + 5 moderators + 3-level hierarchical Bayes + Cluster D multiplicative gate + 8 P1–P8 pre-registered predictions | ✅ shipped |
| **V8 / Cluster E** (πphen) | LINCS L1000 + JUMP-CP + chemCPA + MOFA+ K=30 + joint posterior + 8-cell disagreement + I_novel novel-mechanism score | ✅ shipped; chemCPA on **real LINCS** (Val R²=0.46), hierarchical on **real cpg0000** (R̂=1.010) |
| **Wet-lab shortlist v11 (grid)** | Differentiated (compound×target) grid composer + within-target binding percentile + class-anchored clinical g placement + differentiation guard | ✅ shipped (replaces degenerate v10) |
| **Retrospective clinical validation** (Gap 3) | Leakage-audited 31-drug ledger + 3 leave-out predictors (target / class-LOCO / class-extrapolation) + AUROC/bootstrap/permutation, numpy-only | ✅ shipped; class track-record AUROC 1.00, target affinity ≈ chance |

### Pre-registration

OSF.io pre-registration documents (lock BEFORE unblinding):
- `reports/paper-drafts/v7_osf_preregistration.md` — V7 hierarchical Bayes priors, 12-class PRISMA, 5 moderators, 4 validation gates, 8 P1–P8 predictions, CPT/CPT:PSP fallback
- `reports/paper-drafts/v8_osf_preregistration.md` — V8 MOFA+ K=30, Leiden γ space, AMI/ARI bands, 30-class mechanism taxonomy, 9+1 nootropic anchor set, I_novel novel-mechanism gate, Nat Mach Intell/Nat Methods fallback

---

## Critical constraints

- **RTX 5070 / Blackwell sm_120**: use the PyTorch **nightly cu128** wheel:
  ```powershell
  pip install --pre torch torchvision --index-url https://download.pytorch.org/whl/nightly/cu128
  ```
- **Windows 11 + PowerShell**: env-setup scripts are PowerShell, not bash. Use `conda activate` from PowerShell after `conda init powershell` once.
- **Two-environment split**:
  - `mammal_env`: torch 2.12 nightly cu128 + MAMMAL + Boltz-2 + ADMET-AI + abagen + PyMC + numpyro-optional
  - `txgnn_env`: torch 2.4.0+cu121 + PyG 2.7.0 + DGL 2.4.0 (graphbolt wheel pinned)
- **Rate limits**: PubChem PUG-REST 5 req/s (200 ms sleep); UniProt + ChEMBL polite-UA only; OT Genetics ~10 qps soft limit (parquet caching for re-runs).

---

## Quickstart

### Environment setup

```powershell
# 1. PowerShell on Windows 11
conda create -n mammal_env python=3.10 -y
conda activate mammal_env
pip install --pre torch torchvision --index-url https://download.pytorch.org/whl/nightly/cu128
pip install biomed-multi-alignment[examples]
pip install -e .[dev]
pip install pymc numpyro arviz abagen brainsmash    # V6.B
pip install lightgbm                                  # LambdaMART promotion
pip install pycytominer boto3 cmapPy                 # V8 (optional; ingestion only)
pip install mofapy2                                   # V8 (optional; for full MOFA+ vs SVD fallback)
```

### End-to-end pipeline

```powershell
# Stage 1: targets + compounds
python scripts/02_fetch_targets.py
python scripts/03_fetch_compounds.py

# Stage 2: V4 + V5 production pipeline
python scripts/04_score_dti.py
python scripts/14_v2_cluster_b_admet.py
python scripts/15_v2_fusion.py --calibrated-mammal --znorm-mammal --add-moa-ranker --add-tanimoto-ranker

# Stage 3: V6.A multi-head ensemble
python scripts/52_v6_mmatt_activate.py       # MMAtt-DTA (requires Zenodo weights download)
python scripts/53_v6_mmatt_fusion_ranker.py  # INVERT-mask + 4-head fusion

# Stage 4: V6.B Cluster D
python scripts/54_v6b_cluster_d_foundation.py    # abagen AHBA cache
python scripts/55_v6b_cluster_d_nuts.py          # PyMC NUTS posterior (4 chains × 2000 draws)

# Stage 5: V7 effect-size validation
python scripts/57_v7_validation_gates.py     # P1-P8 + Roberts ceiling gate

# Stage 6: V8 + v11 differentiated wet-lab shortlist (compound × target grid)
python scripts/74_wet_lab_shortlist_v11_grid.py --top-n 50

# Stage 7: retrospective clinical-outcome validation (Gap 3 — leakage-audited)
python scripts/75_retrospective_clinical_validation.py

# Hypothesis audit (re-run any time)
python scripts/41_v5_hypothesis_audit.py
```

### Test suite

```powershell
pytest tests/ -m "not slow"   # 503 pass / 1 skip; ~30 s (test_fetchers.py needs the respx dev dep)
pytest tests/ -m slow         # 12 pass / 2 env-gated (real MAMMAL package); GPU smoke + real model load
```

---

## Repository layout

```
├── design/
│   ├── V4_STATUS_AND_FORWARD_PLAN.md      # source-of-truth status doc; V4→V8 timeline
│   └── V6_ARCHITECTURE_PLAN.md            # phased implementation roadmap V6+V7+V8
├── src/mammal_repurposing/
│   ├── cluster_a/                         # V6.A multi-head DTI
│   │   ├── mmatt_dta_adapter.py           #   MMAtt-DTA (Schulman 2024)
│   │   ├── psichic_adapter.py             #   PSICHIC (Koh 2024)
│   │   ├── balm_adapter.py                #   BALM (Gorantla 2025)
│   │   └── tanimoto_ranker.py             #   Morgan FP × ChEMBL actives
│   ├── cluster_d/                         # V6.B Bayesian neuro prior
│   │   ├── bayesian_prior.py              #   PyMC NUTS hierarchical model
│   │   ├── data_fetchers.py               #   AHBA / OT Genetics / cellxgene
│   │   └── validation_gates.py            #   V6.B.4 4-gate framework
│   ├── cluster_e/                         # V8 πphen
│   │   ├── ingest_lincs.py                #   LINCS L1000 cmapPy WTCS
│   │   ├── ingest_jumpcp.py               #   JUMP-CP S3 sync
│   │   ├── chemcpa_train.py               #   Hetzel 2022 chemCPA
│   │   ├── mofa_embed.py                  #   MOFA+ K=30 joint embedding
│   │   └── joint_phenotype.py             #   V7+V8 joint posterior + 8-cell
│   ├── translation/                       # V7 effect-size translation
│   │   ├── pbpk.py                        #   JAX/diffrax 9-compartment PBPK
│   │   ├── prisma_priors.py               #   Schmidli 2014 12-class MAP
│   │   └── effect_size_model.py           #   3-level hierarchical Bayes
│   ├── calibration/
│   │   ├── isotonic.py                    #   per-target sklearn IsotonicRegression
│   │   └── venn_abers.py                  #   V6.A.4 Venn-ABERS + correlated MC
│   ├── fusion/
│   │   ├── bayesian_router.py             #   V6.A.3 per-target trust matrix
│   │   ├── joint_composition.py           #   v10 + v11-grid wet-lab composer
│   │   └── faceted_shortlist.py           #   V4 mechanism-class facets
│   ├── validation/
│   │   └── retrospective.py               #   Gap 3 leakage-audited clinical validation
│   └── cluster_c/
│       ├── txgnn.py                       #   per-disease ranking API
│       └── primekg.py                     #   PrimeKG loader
├── scripts/                               # 82 end-to-end pipeline scripts
├── reports/                               # auto-generated markdown reports
│   ├── wet_lab_shortlist_v11.md           # ★ production deliverable (grid)
│   ├── retrospective_clinical_validation_v1.md  # ★ Gap 3 headline result
│   ├── cluster_d_nuts_v1.md               # V6.B.3 posterior + convergence
│   ├── v7_validation_v1.md                # V7.4 P1-P8 + Roberts ceiling
│   ├── v7_osf_preregistration.md          # V7 OSF lock
│   ├── v8_osf_preregistration.md          # V8 OSF lock
│   ├── hypothesis_audit_v1.md             # 22 falsifiable hypotheses
│   ├── methodology_v2.md                  # V4+V5 narrative
│   └── ...                                # ~30 other auto-generated reports
├── research/4-tier/                       # 6 research deep-dives
│   ├── Multi Head DTI.md                  # V6.A spec
│   ├── Multi-Source Neurobiological Prior for Cognition Target Prioritization.md  # V6.B
│   ├── Clinical Effect-Size Translation Function.md                                # V7
│   ├── Clinical Effect-Size Translation Function A Methodology Pre-Registration... # V7 companion
│   ├── Perturbational Evidence Axis.md                                              # V8
│   └── Technical Feasibility Deep-Dive Adding a Phenotypic.md                       # V8 companion
└── tests/                                 # 442 non-slow + 14 slow (28 test files)
    ├── test_v7_translation.py
    ├── test_v8_cluster_e.py
    ├── test_v8_advanced.py
    ├── test_grid_composition_v11.py        #   Gap 1 — no-collapse grid composer
    ├── test_retrospective_validation.py    #   Gap 3 — leakage-audited AUROC
    └── ...
```

---

## Methodology contribution

The methodology contribution **distinct from the candidate ranking**:

1. **Diagnostic protocol** that exposes prior collapse and statistical-power blockades before they poison downstream verdicts (V4 commit `530dc40`).
2. **Tanimoto-to-actives baseline** that empirically beats the 458M-param MAMMAL foundation model panel-wide — a publishable negative finding (V4).
3. **Isotonic per-target calibration** that surgically repairs the rare cases where MAMMAL has signal in the wrong direction (V4 `8624fd1`).
4. **Pocket-conditioned classifier** that distinguishes PDE4D allosteric (BPN14770) from catalytic (rolipram) at the geometric level (V4 `3884ba5`).
5. **Faceted shortlist** that dissolves single-target lock-in into structured, mechanism-orthogonal top-N tables with cross-facet provenance (V4 `1c288a8`).
6. **Multi-head DTI ensemble** with Venn-ABERS calibration + per-target Bayesian routing + multi-head disagreement-as-signal axis (V6.A; published `J Cheminform` / `Nat Mach Intell` candidate).
7. **Bayesian Cluster D neurobiological prior** with AHBA + OT Genetics L2G + cellxgene single-cell + PyMC NUTS hierarchical model + 4-gate Roberts 2020 ceiling validation (V6.B; published `Cell Reports Methods` / `Bioinformatics` candidate).
8. **Clinical Effect-Size Translation** — 9-compartment PBPK + Schmidli 2014 robust MAP class priors + 3-level hierarchical Bayes with 5 failure-mode moderators; PET-anchored to Bohnen 2005 / Volkow 1998 / Kapur 2000; 8 pre-registered predictions P1–P8 with falsifiers (V7; CPT/CPT:PSP fallback).
9. **πphen Perturbational Evidence Axis** — first multi-modal phenotypic prior (LINCS L1000 + JUMP-CP morphology + iPSC-MEA + chemCPA imputation) coupled to a target-first Bayesian repurposing posterior with explicit conditional-dependence structure and 4-axis disagreement-as-signal facet tagging including the (L, L, H) novel-mechanism cell (clemastine territory) (V8; Nat Mach Intell / Nat Methods candidate).

---

## Reference citations

Key citations for the methodology — full bibliography in [`design/architecture-and-plans/V4_STATUS_AND_FORWARD_PLAN.md` Appendix A](design/architecture-and-plans/V4_STATUS_AND_FORWARD_PLAN.md):

- Shoshan et al. 2026 *npj Drug Discovery* / arXiv:2410.22367 — **MAMMAL** foundation model
- Schulman et al. 2024 *Bioinformatics* 40(8):btae496 — **MMAtt-DTA** superfamily-conditional
- Koh et al. 2024 *Nat Mach Intell* 6:673 — **PSICHIC** physicochemical contrastive DTI
- Gorantla et al. 2025 *J Chem Inf Model* 65(22):12279 — **BALM** ESM-2 + ChemBERTa-2
- Mervin et al. 2020 *J Chem Inf Model* 60:4546 — **Venn-ABERS** AstraZeneca 40M-pair benchmark
- Markello et al. 2021 *eLife* 10:e72129 — **abagen** AHBA toolbox
- Moodie et al. 2024 *Hum Brain Mapp* 45(4):e26641 — 41-gene cortical *g*-map
- Davies et al. 2018 *Nat Commun* 9:2098 — intelligence GWAS N=300,486
- Roberts et al. 2020 *Eur Neuropsychopharm* 38:40-62 — **THE ceiling paper**
- Subramanian et al. 2017 *Cell* 171(6):1437 — **LINCS L1000** 1.3M signatures
- Chandrasekaran et al. 2024 *Nat Methods* 21(6):1114 — **JUMP-CP** Cell Painting
- Hetzel et al. 2022 NeurIPS — **chemCPA** generative imputation
- Argelaguet et al. 2020 *Genome Biol* 21:111 — **MOFA+** joint factor model
- Bohnen et al. 2005 *Neurology* 64:1037 — donepezil cortical AChE PET anchor
- Volkow et al. 1998 *Am J Psychiatry* 155:1325 — MPH DAT PET dose-response
- Kapur et al. 2000 *Am J Psychiatry* 157:514 — haloperidol D2 striatal PET
- Schmidli et al. 2014 *Biometrics* 70:1023 — robust meta-analytic-predictive priors

---

## License

Apache-2.0. MAMMAL (IBM Research), ADMET-AI (Stanford), Boltz-2 (MIT + Recursion), PrimeKG + TxGNN (Harvard) — all individually Apache-2.0 / MIT.

---

## Citation

If you use this pipeline in your work, please cite:

```bibtex
@misc{lonergan2026mammal,
  title  = {{MAMMAL Cognitive Enhancement Drug Repurposing}: a multi-layer
            Bayesian pipeline with multi-head DTI ensemble, Bayesian
            neurobiological prior, PBPK-anchored effect-size translation,
            and target-agnostic perturbational evidence axis},
  author = {Lonergan, Pierce},
  year   = {2026},
  url    = {https://github.com/pierce-lonergan/MAMMAL_Cognitive_Enhancement_Drug_Repurposing},
  note   = {ORCID 0009-0008-4235-396X; AI-assisted; OSF pre-registered}
}
```

---

*Build status: **V4 → V8 architecture complete + Gaps 1–7 shipped + peer-review round 2 hardening.** V6.B.5 PyMC NUTS converged on the 191-target panel post-MH8 (R̂=1.000, 0 divergences). chemCPA trained on real LINCS L1000 (Val R²=0.46); V8 hierarchical on real cpg0000 (R̂=1.010). Wet-lab shortlist **v11** on the **complete 31-target panel** (all scored with the real MAMMAL DTI head; CHRM1/CHRM4/HTR6 added so CIAS surfaces M1/M4 and AD scores 5-HT6). **Retrospective clinical validation** (Gap 3): class track record AUROC 1.00 (9/9 famous failures flagged) vs affinity at chance. **Disease reframe** (Gap 2): recovers cholinesterase inhibitors for AD (within-disease AUROC 0.97), muscarinic M1/M4 for schizophrenia, PDE4 for Fragile X. **Allosteric learn-to-rank** (Gap 4): fused head lifts within-target ρ +0.02→+0.51. **Clinician dossiers** (Gap 5), **external benchmark** (Gap 6), and a **prospective repurposing shortlist** (Gap 7: CIAS→buspirone/M1, FXS→roflumilast). **Review round 2**: class-level cluster bootstrap, KG network-propagation + structure-similarity comparators (both 0.80, hindsight-confounded; class 1.00), CONSORT-style ledger flow, 2× methods. **Review round 3**: pseudo-prospective temporal hold-out + prequential as-of (AUROC 1.00 informed / 0.96 full, 1 honest miss) and taxonomy sensitivity (coarse 0.62, random 0.46, 0/2000 reach 1.0) — the signal is granularity-specific and would have flagged the 2014–2022 failure wave before readout. **Review round 4**: calibrated constructive predictor (class-only Brier 0.05; target-centric features add nothing); **pre-registered prospective predictions on real ongoing trials — 2/2 confirmed (iclepertin/GlyT1, luvadaxistat/DAAO both failed as predicted), 4 pending falsifiable**; ledger expansion to n=42 (17/17 classes still pure); **unbiased ClinicalTrials.gov pull** (pre-specified query, 294-trial denominator across 4 indication×phase slices, outcome-blind adjudication → n=47/20 classes still pure, AUROC 0.983) — the curation-artifact objection answered with a programmatic query. 503 non-slow + 12/14 slow pytest pass. See `PROJECT_STATUS.md` and `design/architecture-and-plans/V4_STATUS_AND_FORWARD_PLAN.md` for full status.*
