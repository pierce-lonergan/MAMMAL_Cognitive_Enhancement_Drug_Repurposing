# MAMMAL Cognitive Enhancement Drug Repurposing

An in silico drug-repurposing pipeline built on IBM Research's [MAMMAL](https://github.com/BiomedSciAI/biomed-multi-alignment) foundation model. Scores candidate compounds against a curated panel of cognition-relevant human protein targets and ranks them by predicted binding affinity (pKd).

## Honest Scope

MAMMAL has no cognition-specific training data and no behavioral benchmark. It is a **target-affinity oracle**, not a "working memory enhancer" predictor. This pipeline:

- Uses the released `dti_bindingdb_pkd` head to score (compound, target) pairs.
- Ranks compounds within a hand-curated 22-target panel covering cholinergic, glutamatergic, dopaminergic, noradrenergic, histaminergic, orexinergic, PDE, BDNF/TrkB, sigma-1, and ion-channel mechanisms.
- Sanity-checks predictions against positive controls (donepezil/AChE, methylphenidate/DAT, pitolisant/HRH3, etc.) before trusting any output.

**Treat predictions as rank-ordering signals, not as Kd estimates.** This is a prefilter for downstream wet-lab or literature work, not a substitute for it.

See [`research/compass_artifact_wf-...md`](research/compass_artifact_wf-5221caf1-01b1-4d5c-8b44-e6657ccb9c8e_text_markdown.md) for the full technical deep-dive that motivates this design.

## Full Pipeline (Phases 0-5)

**MVP (Phase 0) — runs the core MAMMAL DTI scoring + decision gate.**

```
Phase 0.1: Environment setup        → scripts/01_setup_env.ps1
Phase 0.2: Fetch target sequences   → scripts/02_fetch_targets.py
Phase 0.3: Build compound library   → scripts/03_fetch_compounds.py
Phase 0.4: DTI scoring              → scripts/04_score_dti.py
Phase 0.5: Sanity check (gate)      → scripts/05_sanity_check.py    ★ go/no-go gate
```

**Phase 1 — Triangulation stack (parallel-safe, all required).**

```
Phase 1.1: BBBP + ClinTox heads     → scripts/06_score_aux_heads.py
Phase 1.2: ChEMBL ground truth      → scripts/07_chembl_evidence.py
Phase 1.3: OpenTargets context      → scripts/08_opentargets_context.py
```

**Phase 2-5 — Composites, calibration, audit, benchmark, shortlist.**

```
Phase 2:   Cognitive composites     → scripts/09_cognitive_composites.py
Phase 3.1: Platt calibration        → scripts/10_calibration.py    ★ rank vs absolute gate
Phase 3.2: Allosteric audit         → scripts/11_allosteric_audit.py
Phase 4.1: Allosteric benchmark     → scripts/12_allosteric_benchmark.py (publishable methods note)
Phase 5:   Wet-lab shortlist        → scripts/13_wet_lab_shortlist.py   ★ final actionable artifact
```

Deferred (v3): LINCS L1000 signature reversal (requires clue.io key), MCP server integration, ISR fine-tune (Phase 4.2), L1000 reference signature (Phase 4.3).

## Quickstart

### 1. Setup (one-time, ~10-30 min)

Requires conda/miniconda on Windows 11. From PowerShell in this repo:

```powershell
.\scripts\01_setup_env.ps1
```

This creates a `mammal_env` conda environment, installs **PyTorch nightly cu128** (required for Blackwell sm_120 on the RTX 5070), installs `biomed-multi-alignment[examples]`, installs this package in editable mode, and triggers the ~1.8 GB MAMMAL weight download.

### 2. Run the pipeline

```powershell
conda activate mammal_env

python scripts/02_fetch_targets.py
python scripts/03_fetch_compounds.py
python scripts/04_score_dti.py
python scripts/05_sanity_check.py
```

Outputs land in `data/results/`:
- `dti_scores.parquet` — full grid: ~6,600 rows of (target, compound, predicted_pKd)
- `sanity_report.md` — positive-control validation, polypharmacology leaderboard

### 3. Validate via smoke notebook (optional but recommended)

```powershell
jupyter notebook notebooks/00_smoke_test.ipynb
```

Runs IBM's reference DTI inference call locally to confirm env + GPU + weights work.

## CLI

After install, a `mammal-repurposing` console script is available:

```powershell
mammal-repurposing fetch-targets
mammal-repurposing fetch-compounds
mammal-repurposing score
mammal-repurposing sanity
```

## Project Layout

```
.
├── data/
│   ├── raw/           # hand-curated seed CSVs (committed)
│   ├── interim/       # enriched target/compound parquets (not committed)
│   └── results/       # scoring output + sanity report (not committed)
├── src/mammal_repurposing/
│   ├── config.py      # model IDs, normalization constants, target panel
│   ├── fetchers/      # uniprot, pubchem, chembl
│   ├── scoring/       # MAMMAL inference (model_loader, dti, runner)
│   ├── analysis/      # sanity, polypharmacology
│   └── cli.py
├── scripts/           # numbered one-shot pipeline scripts
├── notebooks/         # smoke test + exploratory notebooks
├── tests/
└── research/          # background research deep-dive
```

## Decision Gate

The sanity check is a **go/no-go gate**, not just a status report. If positive controls (donepezil at ACHE, methylphenidate at SLC6A3, pitolisant at HRH3) do not rank in the top 20% of their target's score distribution, the script exits nonzero and you should **stop and debug prompt formatting before trusting any other output**.

## License

Apache-2.0. MAMMAL itself is also Apache-2.0 (IBM Research).
