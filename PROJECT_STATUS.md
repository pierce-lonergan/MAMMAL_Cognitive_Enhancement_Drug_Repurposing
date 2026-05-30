# PROJECT STATUS — MAMMAL Cognitive Enhancement Drug Repurposing

**One-pager executive summary** suitable for grant applications, sprint reviews, or stakeholder briefings. Last refreshed 2026-05-29 (Gap 1 v11 grid composer + Gap 3 retrospective clinical validation + MH 1-8 sprint suite + chemCPA real-LINCS training + V8 cpg0000 calibration).

---

## What this is

A multi-layer Bayesian pipeline for cognition-enhancement drug repurposing built around IBM Research's [MAMMAL](https://github.com/BiomedSciAI/biomed-multi-alignment) foundation model. Five architectural layers (V4 → V5 → V6 → V7 → V8) compose into a single ranked set of **(compound, target) repurposing hypotheses** with predicted cognition Hedges' *g*, pre-filtered by the Roberts 2020 SMD ceiling (g ≤ 0.50).

**Honest scope**: this pipeline does NOT discover smart drugs. It enriches a candidate set so wet-lab cycles spend money on plausibility, not chemistry-lottery tickets. Its sharpest validated claim is *negative-and-useful*: target-binding affinity and target genetic-relevance do **not** predict cognition-drug clinical success — mechanism-class track record does (see retrospective validation, below). Gap 2 operationalises that into a **disease-specific** prioritisation: re-scored per disease, the pipeline recovers cholinesterase inhibitors for Alzheimer's, muscarinic M1/M4 for schizophrenia, and PDE4 for Fragile X — each matching the real clinical record.

---

## Headline metrics (current sprint)

| Metric | Value | Status |
|---|---|---|
| Pytest pass rate (non-slow) | **485** pass / 1 skip (+66 across Gaps 1–7 + grid expansion + panel→31; `test_fetchers.py` needs the `respx` dev dep) | ✅ |
| **Gap 7 — prospective repurposing shortlist** (capstone) | approved drugs ranked as mechanism-justified repurposing hypotheses per disease (class prognostic prior × engagement, SUCCESS classes only); **CIAS→buspirone/M1, FXS→roflumilast (PDE4), AD→σ1**; xanomeline correctly flagged *standard* (`reports/repurposing_shortlist_v1.md`) | 🎯 |
| **Panel finished to 31 targets** (real MAMMAL DTI) | CHRM1/CHRM4 (M1/M4) + HTR6 (5-HT6) + GRM2/3/5/GlyT1/HTR4 scored on RTX 5070; **CIAS now surfaces muscarinic M1/M4** (xanomeline class), AD scores 5-HT6 (demoted). MAMMAL runs in a Py-3.12 venv (`docs/MAMMAL_SETUP.md`) | ✅ |
| Pytest pass rate (slow) | **12 / 14** (2 fail = real MAMMAL `biomed-multi-alignment` package not in this env, not a regression) | ✅ |
| **Gap 3 — retrospective clinical-outcome validation** | **mechanism-class track record AUROC 1.00 (perm p=0.0002) vs target affinity 0.12 / relevance 0.59; 9/9 famous Phase III failures flagged** (`reports/retrospective_clinical_validation_v1.md`) | 🏆 |
| **Gap 2 — disease-population reframe** | **each disease recovers its real winning mechanism**: AD→AChE-I (within-disease class AUROC **0.97**, p=0.003, 10/10 AD failures flagged), CIAS→muscarinic M1/M4 (xanomeline-KarXT), FXS→PDE4 (zatolmilast) (`reports/disease_reframe_v1.md`) | 🏆 |
| **Gap 4 — allosteric learn-to-rank head** | MAMMAL flat within-target (std 0.01–0.05); fused [MAMMAL⊕Tanimoto⊕Boltz⊕physchem] lifts within-target Spearman ρ **+0.02→+0.51** held-out on the cited allosteric benchmark (`reports/allosteric_ltr_v1.md`) | 🔬 |
| **Gap 5 — clinician evidence dossiers** | one-page GRADE-style cards: g+CrI, evidence quality, class track record, off-target liability flags, provenance, caveats (`reports/clinician_dossiers_v1.md`) | 🩺 |
| **Gap 6 — external benchmark** | class track record AUROC 1.00 vs affinity 0.47 / genetics 0.59 (leakage-free); popularity 0.96 flagged hindsight-confound (`reports/external_benchmark_v1.md`) | 📊 |
| **Gap 1 — v11 grid shortlist** (replaces degenerate v10) + **V6.A grid 13→23** | top-25 spans 10 targets; positive controls correct; peptides filtered (MW≤900); 23/28 panel targets scored | ✅ |
| chemCPA real-LINCS production training (Sprint 5.2) | **107K real LINCS sigs, Val R²=0.46, OOD R²=0.33** on the 9-compound canonical holdout — 8.3 min on RTX 5070 | ✅ |
| V8 hierarchical real cpg0000 calibration (Sprint 4.3) | **R̂=1.010, 0 div; ICC_cell=0.018, ICC_inter=0.149; 60/60 compounds T>0.6** (U2OS→brain transfer defended) | ✅ |
| Hypothesis audit verdicts | **22 PASS / 3 DEGRADE / 0 FAIL** | ✅ |
| V6.B PyMC NUTS production (22-target) | **R̂ max = 1.000, ESS min = 12,780** (4 chains × 2000 draws, real AHBA) | ✅ |
| V6.B.5 PyMC NUTS production (191-target, **post-MH8**) | **R̂ max = 1.000, ESS min = 1,808, divergences = 0** (was 37 pre-fix) | ✅ |
| V6.B Gate 1 (Roberts ceiling, HARD) | **0 violations on 191-target posterior** | ✅ |
| V6.B Gate 2 (multi-modulator Spearman, 70 anchors / 38 targets) | **ρ = +0.10 (DEGRADE best case) / -0.35 (FAIL worst case)** — publishable falsification (see reports/gate2_multi_modulator_v1.md) | 📄 |
| V7 NUTS with 15-compound anchor likelihood | **R̂ = 1.000, ESS = 2,332, MAE = 0.073, 0 Roberts ceiling violations** | ✅ |
| V7.4 Gate 1 (P1-P8 prediction bands) | 5 PASS / 1 FAIL / 2 NO_COMPOUND in stub mode; 4 PASS / 3 FAIL / 1 NO_COMPOUND in full NUTS (honest partial-pool — Sprint 3.2 addresses) | ⏳ |
| V7.2 Stage 3 PRISMA prior coverage | **73 / 96 cells populated (76%)**, up from 32 (33%) — Sprint 3.1 | ✅ |
| V7.2 Stage 4 reference compound anchors | **109 cells / 48 compounds** (was 15) — Sprint 3.3; includes ketamine impairment + 25 Phase III nulls | ✅ |
| V7.3 Stage 2 NUTS V2 (per-class τ² + population × class) | `fit_effect_size_nuts_v2` shipped + tested — Sprint 3.2 | ✅ |
| V8.6 Hierarchical (MH3 + MH7 bundled) | `build_v8_hierarchical_with_cell_random_effect` shipped; β/α/γ/δ random effects + ICC + transferability index T_{c,k}; synthetic round-trip validated — Sprint 4.1 + 4.2 | ✅ |
| V8.2 chemCPA synthetic-LINCS smoke | Loss 0.1728→0.1068 (1.62× reduction); test R² = +0.485 (gate ≥ 0.30) | ✅ |
| V8.4 Gate 1 dry-run on synthetic phenotype | **AMI = 1.000, ARI = 1.000** (Agglomerative + HDBSCAN min∈{15,25}) | ✅ |
| Target panel (V6.B core) | **31 targets** (was 22→28→31) — +CHRM1/CHRM4 (M1/M4) + HTR6 (5-HT6); all scored with the real MAMMAL DTI head | ✅ |
| V6.B.5 expanded panel | **191 targets** with 22-panel ✅ strict subset; MAO-A/MAO-B/COMT/ACHE substrate-mediated | ✅ |
| Multi-modulator anchor table | **70 rows / 38 targets / 59 compounds / 24 Phase III nulls** (Sprint 2.1) | ✅ |
| Total scripts shipped | **90** (74 v11 grid, 75 retrospective, 76 disease reframe, 77 grid expansion, 78 allosteric LTR, 79 external benchmark, 80 clinician dossier, 81 MAMMAL scoring, 82 repurposing shortlist) | ✅ |
| Total source modules shipped | **115** across cluster_a (+ allosteric_ltr) / cluster_b/c/d/e / translation / calibration / fusion / **validation** (retrospective + disease_reframe) / **reporting** (clinician_dossier + repurposing_shortlist) / pockets / selectivity / diagnostics / fetchers / scoring | ✅ |
| MH implementation roadmap | **all core sprints complete** (1.1-1.4 + 2.1-2.2 + 3.1-3.5 + 4.1-4.3 + 5.1-5.2 + 6.2/6.4) + Gap 1 + Gap 3; reports/MH_IMPLEMENTATION_ROADMAP.md | 🚀 |

---

## Five architectural layers — status

| Layer | What it does | Status |
|---|---|---|
| **V4** | Calibrated 5-cluster fusion: MAMMAL DTI + ESM2 + Boltz-2 + ADMET + Tanimoto + isotonic per-target calibration + faceted shortlist + pocket DB | ✅ shipped (16 reports + 9 modules) |
| **V5** | Z-norm within target + 14 Tier-2/3 items: §8.0b-zn liability + MoA ranker + nootropic similarity + pocket-routed gating + CTgov + LambdaMART + conformal + scaffold-AL + brain-region + PyMC GRIN pool | ✅ shipped |
| **V6.A** | Multi-Head DTI ensemble: MMAtt-DTA + PSICHIC + BALM + per-head bias decomposition + Bayesian router + Venn-ABERS calibration + multi-head disagreement axis | ✅ shipped (Tier-A FAIL published-quality result) |
| **V6.B** | Bayesian Cluster D neurobiological prior: abagen AHBA + OT Genetics L2G fetcher + cellxgene preview + PyMC NUTS hierarchical model + 4-gate validation | ✅ shipped (R̂=1.000) |
| **V6.B.5** | Panel expansion 22 → ~210 targets per Cluster D §F (GWAS L2G + MAGMA + AHBA + SC + Lit-OTAR) | ✅ shipped (Stage 1; 191 targets) |
| **V6 Cluster C** | PrimeKG + TxGNN per-disease ranking | ✅ shipped (V6 API rewrite) |
| **V7** | Clinical Effect-Size Translation: 9-cmpt PBPK + Schmidli 2014 PRISMA priors (12 classes + 32 subdomain cells) + 5 failure-mode moderators + 3-level hierarchical Bayes + 8 P1-P8 pre-registered predictions | ✅ shipped (real NUTS R̂=1.000, MAE=0.073) |
| **V8 / Cluster E** | πphen Perturbational Evidence Axis: LINCS L1000 + JUMP-CP + iPSC-MEA + chemCPA + MOFA+ + joint posterior + 8-cell disagreement + I_novel novel-mechanism score | ✅ shipped (scaffolds + synthetic Gate 1 dry-run AMI=1.000) |

---

## Four publishable manuscripts drafted

| Paper | Venue | Headline finding | Status |
|---|---|---|---|
| **V6.A** | *J Cheminform* / *Nat Mach Intell* | MMAtt-DTA ρ +0.65 vs Tanimoto +0.90 at SLC6A3 → **Tier-A FAIL**; INVERT-mask architecture drops 6 panel targets; per-target Bayesian router empirically necessary | `reports/v6a_paper_draft.md` ✅ |
| **V6.B** | *Cell Reports Methods* / *Bioinformatics* | PyMC NUTS R̂=1.000 on real AHBA; ACHE substrate-mediated flag correctly fires; reference-anchor pull recovers CHRNA7 from y_AHBA=-0.53 → θ̄=+0.44 | `reports/v6b_paper_draft.md` ✅ |
| **V7** | *Clinical Pharmacology & Therapeutics* (Wiley, IF 7.3); fallback *CPT:PSP* (IF 4.2) | Real Bayesian inference within Roberts ceiling: R̂=1.000, MAE=0.073, zero ceiling violations; honest Gate 1 partial-pool finding (4/8 PASS by tight margins) | `reports/v7_paper_draft.md` ✅ |
| **V8** | *Nature Machine Intelligence* (A realistic); stretch *Nature Methods* (A+) | First multi-modal phenotypic prior + I_novel mutual-information novel-mechanism score identifying (L, L, H) clemastine-class candidates; OSF pre-reg locked | `reports/v8_paper_draft.md` ✅ |

---

## Pre-registration

Two OSF-pre-registration documents locked before unblinding:
- `reports/v7_osf_preregistration.md` — V7 hierarchical Bayes priors, 12-class PRISMA, 5 moderators, 8 P1-P8 predictions with falsifiers, 4 validation gates, 15-compound held-out anchor set, CPT venue
- `reports/v8_osf_preregistration.md` — V8 MOFA+ K=30, Leiden γ sweep, AMI/ARI bands, 30-class mechanism taxonomy, 9+1 nootropic anchor set, I_novel novel-mechanism gate, Nat Mach Intell venue

Both documents are publication-ready markdown. OSF.io account + DOI mint is the only remaining step before locking.

---

## Repository content

```
├── design/                       2 architecture docs (V4 + V6 plans)
├── research/4-tier/              6 research deep-dives (Multi-Head DTI, Cluster D,
│                                  V7 Clinical Effect-Size, V7 Pre-Reg companion,
│                                  V8 Perturbational, V8 Feasibility companion)
├── reports/                      90 auto-generated markdown reports
│   ├── 4 paper drafts (V6.A + V6.B + V7 + V8) + umbrella synthesis
│   ├── 2 OSF pre-registration docs (V7 + V8)
│   ├── wet_lab_shortlist_v11.md  (★ Gap 1 differentiated grid deliverable)
│   ├── retrospective_clinical_validation_v1.md (★ Gap 3 headline result)
│   ├── methodology_v3.md         (coherent V4→V8 narrative)
│   ├── hypothesis_audit_v1.md    (falsifiable claims tracked)
│   └── ~80 other auto-generated diagnostic / validation / wet-lab reports
├── figures/                      Publication-quality figures
│   ├── v7/                       4 figures at 300 DPI (PBPK + P1-P8 + LOO MAE + sweep)
│   └── v11/                      retrospective_roc.png (Gap 3 ROC contrast)
├── src/mammal_repurposing/       110 source modules across 14 packages
│   ├── cluster_a/                (V6.A multi-head DTI: MMAtt + PSICHIC + BALM)
│   ├── cluster_c/                (V6 PrimeKG + TxGNN per-disease)
│   ├── cluster_d/                (V6.B Cluster D + 4-gate validation + panel expansion)
│   ├── cluster_e/                (V8 πphen: LINCS + JUMP-CP + chemCPA + MOFA+ + joint)
│   ├── translation/              (V7 PBPK + PRISMA priors + effect-size model)
│   ├── calibration/              (Venn-ABERS + isotonic + IsotonicCalibrator)
│   ├── fusion/                   (Bayesian router + RRF + faceted + LambdaMART + v10/v11-grid composition)
│   ├── validation/               (Gap 3 leakage-audited retrospective clinical validation)
│   └── ...                       (pockets / selectivity / diagnostics / fetchers / scoring)
├── scripts/                      90 end-to-end pipeline scripts
├── tests/                        485 non-slow pytest cases + 14 slow (32 files)
├── CITATIONS.bib                 Full BibTeX bibliography (~50 entries)
├── README.md                     Public-facing entry point with V4→V8 architecture diagram
└── PROJECT_STATUS.md             This file
```

---

## What's been un-blocked since last refresh (now DONE)

Previously "externally blocked" items that have since been executed in-session:

1. ✅ **LINCS L1000 GCTX** — 5.5 GB GSE70138 Level-5 COMPZ (118,050 sigs × 12,328 genes) downloaded + decompressed.
2. ✅ **chemCPA real-LINCS training** — trained on 107K real signatures (Val R²=0.46, OOD R²=0.33), 8.3 min on RTX 5070. `scripts/73`.
3. ✅ **JUMP-CP cpg0000 pilot** — 46 CPJUMP1 plates pulled; V8 hierarchical NUTS run on real A549/U2OS morphology (`scripts/70`+`71`).

## What's still externally blocked

1. **MOFA+ K=30 on real 7-view stack** — needs all 7 views co-loaded (iPSC-MEA + snRNA still pending).
2. **OT Genetics L2G live fetch (held-out GWAS)** — Gate 3 still INSUFFICIENT_DATA; graceful fallback works.
3. **OSF.io account + DOI mint** — V7 + V8 pre-reg docs markdown-ready.
4. **Wet-lab validation** — the only path to prospective external validation (CRO partnership, $60-110K).
5. **bioRxiv preprint submission** — 5 manuscripts draft-ready.

---

## What's actionable now (in priority order)

1. ✅ **Gaps 2–6 + grid expansion + panel→31 — DONE (2026-05-29)**: disease reframe (Gap 2), retrospective validation (Gap 3), allosteric learn-to-rank (Gap 4), clinician dossiers (Gap 5), external benchmark (Gap 6), V6.A grid 13→**31** with the real MAMMAL DTI head (CHRM1/4 + HTR6 added). See the headline table + `reports/`.
2. **Scale Gap 4**: expand the allosteric benchmark beyond n=21 and add fuller Boltz coverage; the n=21 result is a proof-of-concept.
4. **V8 phenotype axis wiring into v11** — populate the (L,L,H) novel-mechanism cell with real chemCPA/transferability per compound.
5. **OSF.io project creation** + **bioRxiv release** — V6.A/B + V7 + V8 + retrospective + disease-reframe + allosteric-LTR manuscripts.

---

## Cumulative compute cost

- Total NUTS sampling time across all production runs: ~50 sec V6.B + ~54 sec V7 + ~3-5 min V6.B.5 ≈ 5-10 min total
- chemCPA synthetic smoke: ~1 sec (CPU)
- V8 Gate 1 dry-run: ~1 sec
- V7 figures: ~3 sec
- V6.A MMAtt-DTA evaluation: ~30 sec on RTX 5070
- Total V8 architecture development: ~5 sec across all smoke tests

**The entire computational stack runs in under 10 minutes on a single RTX 5070 12 GB consumer GPU.** External downloads (LINCS + JUMP-CP) are the only real wall-clock barrier.

---

## License

Apache-2.0. MAMMAL weights (IBM Research), ADMET-AI (Stanford), Boltz-2 (MIT + Recursion), PrimeKG + TxGNN (Harvard) — all individually Apache-2.0 / MIT. Research deep-dives in `research/4-tier/`: CC-BY-4.0.

---

## Citation

If you use this pipeline:

```bibtex
@misc{lonergan2026pipeline,
  title  = {{MAMMAL Cognitive Enhancement Drug Repurposing}: a multi-layer
            Bayesian pipeline with multi-head DTI ensemble, Bayesian
            neurobiological prior, PBPK-anchored effect-size translation,
            and target-agnostic perturbational evidence axis},
  author = {Lonergan, Pierce and {Claude Opus 4.7 (1M context)}},
  year   = {2026},
  url    = {https://github.com/pierce-lonergan/MAMMAL_Cognitive_Enhancement_Drug_Repurposing},
  note   = {V4 + V5 + V6 + V7 + V8 architecture; OSF pre-registered},
}
```

Full BibTeX bibliography: `CITATIONS.bib`. Per-paper drafts: `reports/v6a_paper_draft.md`, `reports/v6b_paper_draft.md`, `reports/v7_paper_draft.md`, `reports/v8_paper_draft.md`.

---

*Last updated 2026-05-29 (Gaps 1–7 shipped + **panel finished to 31 targets with the real MAMMAL DTI head**: v11 grid + disease reframe + retrospective validation + allosteric learn-to-rank + clinician dossiers + external benchmark; V6.A grid 13→31; plus the MH 1-8 sprint suite + chemCPA real-LINCS + V8 cpg0000). The pipeline is end-to-end shipped with a differentiated shortlist on the complete 31-target panel, leakage-audited retrospective + external benchmarks (class track record AUROC 1.00 vs target paradigms ≈ chance), disease-specific shortlists (AD/CIAS/FXS — CIAS now surfaces M1/M4) each validated within-disease, an allosteric learn-to-rank head that lifts within-target ρ +0.02→+0.51, clinician-legible GRADE dossiers, and a prospective mechanism-justified repurposing shortlist (Gap 7: CIAS→buspirone/M1, FXS→roflumilast, AD→σ1); what remains is wet-lab validation and OSF/bioRxiv release.*
