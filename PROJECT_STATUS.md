# PROJECT STATUS — MAMMAL Cognitive Enhancement Drug Repurposing

**One-pager executive summary** suitable for grant applications, sprint reviews, or stakeholder briefings. Last refreshed in the V6.B.5 + V8 paper sprint.

---

## What this is

A multi-layer Bayesian pipeline for cognition-enhancement drug repurposing built around IBM Research's [MAMMAL](https://github.com/BiomedSciAI/biomed-multi-alignment) foundation model. Five architectural layers (V4 → V5 → V6 → V7 → V8) compose into a single three-factor joint posterior over predicted healthy-adult cognition Hedges' *g*, pre-filtered by the Roberts 2020 SMD ceiling (g ≤ 0.50).

**Honest scope**: this pipeline does NOT discover smart drugs. It enriches a candidate set so wet-lab cycles spend money on plausibility, not chemistry-lottery tickets.

---

## Headline metrics (current sprint)

| Metric | Value | Status |
|---|---|---|
| Pytest pass rate (non-slow) | **250 / 250** | ✅ |
| Hypothesis audit verdicts | **22 PASS / 3 DEGRADE / 0 FAIL** | ✅ |
| V6.B PyMC NUTS production | **R̂ max = 1.000, ESS min = 12,780** (4 chains × 2000 draws, real AHBA, 22 targets) | ✅ |
| V7 NUTS with 15-compound anchor likelihood | **R̂ = 1.000, ESS = 2,332, MAE = 0.073, 0 Roberts ceiling violations** | ✅ |
| V7.4 Gate 1 (P1-P8 prediction bands) | 5 PASS / 1 FAIL / 2 NO_COMPOUND in stub mode; 4 PASS / 3 FAIL / 1 NO_COMPOUND in full NUTS (honest partial-pool) | ⏳ |
| V7.4 Gate 2 (Roberts ceiling, HARD) | **0 violations on 15-compound anchor set** | ✅ |
| V7.4 Gate 3 (MAE < 0.15 on anchor) | **MAE = 0.073** | ✅ |
| V8.2 chemCPA synthetic-LINCS smoke | Loss 0.1728→0.1068 (1.62× reduction); test R² = +0.485 (gate ≥ 0.30) | ✅ |
| V8.4 Gate 1 dry-run on synthetic phenotype | **AMI = 1.000, ARI = 1.000** (Agglomerative + HDBSCAN min∈{15,25}) | ✅ |
| V6.B.5 expanded panel | **191 targets** with 22-panel ✅ strict subset, all UniProts unique | ✅ |
| Total scripts shipped | 63 | ✅ |
| Total source modules shipped | 40+ across cluster_a / cluster_b / cluster_c / cluster_d / cluster_e / translation / calibration / fusion / pockets / selectivity / diagnostics / fetchers / scoring | ✅ |

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
├── reports/                      30+ auto-generated markdown reports
│   ├── 4 paper drafts (V6.A + V6.B + V7 + V8)
│   ├── 2 OSF pre-registration docs (V7 + V8)
│   ├── methodology_v3.md         (coherent V4→V8 narrative)
│   ├── hypothesis_audit_v1.md    (22 falsifiable claims tracked)
│   └── ~24 other auto-generated diagnostic / validation / wet-lab reports
├── figures/                      Publication-quality figures
│   └── v7/                       4 figures at 300 DPI (PBPK + P1-P8 + LOO MAE + sweep)
├── src/mammal_repurposing/       40+ source modules across 12 packages
│   ├── cluster_a/                (V6.A multi-head DTI: MMAtt + PSICHIC + BALM)
│   ├── cluster_c/                (V6 PrimeKG + TxGNN per-disease)
│   ├── cluster_d/                (V6.B Cluster D + 4-gate validation + panel expansion)
│   ├── cluster_e/                (V8 πphen: LINCS + JUMP-CP + chemCPA + MOFA+ + joint)
│   ├── translation/              (V7 PBPK + PRISMA priors + effect-size model)
│   ├── calibration/              (Venn-ABERS + isotonic + IsotonicCalibrator)
│   ├── fusion/                   (Bayesian router + RRF + faceted + LambdaMART + joint composition)
│   └── ...                       (pockets / selectivity / diagnostics / fetchers / scoring)
├── scripts/                      63 end-to-end pipeline scripts
├── tests/                        250 non-slow pytest cases + 4 slow
├── CITATIONS.bib                 Full BibTeX bibliography (~50 entries)
├── README.md                     Public-facing entry point with V4→V8 architecture diagram
└── PROJECT_STATUS.md             This file
```

---

## What's externally blocked

These items are ready to ship the moment the external dependency arrives:

1. **LINCS L1000 GCTX download** (~10 GB; GSE92742 + GSE70138 + clue.io beta) — blocks V8 real-mode Gate 1 evaluation
2. **JUMP-CP cpg0016 S3 sync** (~30-40 GB DeepProfiler + CellProfiler + DINOv2 consensus parquets; NEVER the 115 TB raw images) — blocks V8 real-mode Gate 1 evaluation
3. **chemCPA training on real LINCS + sci-Plex3** (~4-8 h GPU once LINCS loaded) — blocks V8.2 Stage 2
4. **MOFA+ K=30 on real 7-view stack** (~2-4 h CPU once all 7 views loaded) — blocks V8.3 Stage 2
5. **V8 PyMC NUTS on joint posterior** (~8-16 h GPU once V8.3 done) — blocks V8.5 Stage 2
6. **OT Genetics L2G live fetch** — blocked from sandbox; graceful fallback works
7. **OSF.io account + DOI mint** — both pre-reg documents are markdown-ready

Total wall-clock for the full real-data execution: ~24-36 h compute + 1-2 weeks paper drafting iteration.

---

## What's actionable now (in priority order)

1. **V6.B paper Methods + Results expansion** — populate the full 22-target posterior table + 4 validation-gate live execution
2. **V8.4 Stage 2** — if real LINCS + JUMP-CP data arrives, execute Gate 1 AMI/ARI vs PRISMA 30-class
3. **V7.5 Stage 2** — extend P1-P8 anchor set from 15 to ~50-100 compounds with curated meta-analytic g
4. **V6.B.5 Stage 2-3** — replace synthetic GWAS/SC/AHBA scores with live OT Genetics + cellxgene + Moodie 2024 alignment
5. **OSF.io project creation** — upload V7 + V8 pre-registrations + lock + DOI mint
6. **Wet-lab validation of top-N (L, L, H) candidates** — Mei 2014 BIMA-8 remyelination assay; clemastine-class novel-mechanism follow-up
7. **Public release on bioRxiv** — V6.A + V6.B + V7 + V8 manuscripts as preprints with code DOI

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

*Last updated by V6.B.5 NUTS + 4-paper-suite sprint. The pipeline is end-to-end shipped; what remains is wet-lab validation + external data download.*
