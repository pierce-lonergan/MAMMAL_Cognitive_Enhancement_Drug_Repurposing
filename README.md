# MAMMAL Cognitive Enhancement Drug Repurposing

A four-cluster hybrid pipeline for drug repurposing on cognition-relevant targets, built around IBM Research's [MAMMAL](https://github.com/BiomedSciAI/biomed-multi-alignment) foundation model and augmented with structure-aware affinity (Boltz-2), comprehensive ADMET (ADMET-AI), and mechanism reasoning (PrimeKG + TxGNN). Runs on a single 12 GB consumer GPU.

## Honest scope

This pipeline does **not** discover "smart drugs." It enriches a candidate set so wet-lab cycles spend money on plausibility, not chemistry-lottery tickets. Roberts CA et al. (*Eur Neuropsychopharm* 2020, PMID 32709551) puts the effect-size ceiling for healthy-adult cognitive enhancement at SMD ≈ 0.21 for the strongest known agents (methylphenidate). The deliverable here is a calibrated, provenance-rich ranking + a publishable methodology contribution — not a miracle compound.

## Current state (v3 in progress)

**Sprint commit log**: 3baf422 → ed4761c → 7467ac1 → 8f75a32 → … (`git log --oneline`)

| Component | Status | Evidence |
|---|---|---|
| **Cluster A.1** — MAMMAL DTI head | ✅ live | 6,556 pairs scored; 4/7 v1 positive controls top-20% |
| **Cluster A.2** — ESM2-650M target embeddings | ✅ live | 22 targets cached; paralog cos = 0.997 |
| **Cluster A.3** — Boltz-2 affinity (Windows) | ✅ live (slow) | ~149 s/pair, `--no_kernels` fallback |
| **Cluster A.3** — Boltz-2 affinity (WSL2 + cuEquiv kernels) | ✅ live (fast) | ~23 s/pair, 5.6× speedup confirmed |
| **Cluster B** — ADMET-AI 41 endpoints + hard gates | ✅ live | 53 PASS / 64 FLAG / 181 CUT of 298 compounds |
| **Cluster C.1** — PrimeKG loader | 🔄 sprint Phase B (next) | code stub + cognition virtual anchor (5 EFO union) ready |
| **Cluster C.2** — TxGNN scorer | 🔄 sprint Phase B (next) | stub ready; install path via WSL2 + PyG cu128 verified |
| **Fusion** — RRF over MAMMAL + ADMET | ✅ live (2-of-4 clusters) | `scripts/15_v2_fusion.py` |
| **Provenance** — disagreement archetypes + funnel narrative | ✅ live | `provenance/*.py`, renders on every fusion run |
| **Phase 0.5 CHRNA7 rescue gate** | ✅ **PASSED** | TC-5619 100%, encenicline 80% (vs v1 19%, 7%) |
| **Phase 0.4 full Boltz sweep** | 🔄 **running overnight** in WSL2 | 1,165 pairs ETA ~10 h; check via `scripts/_wsl2_sweep_status.sh` |
| **Phase 3.1 calibration** | 🔄 sprint Phase A (now) | gated on T2 (local ChEMBL SQLite mirror) |

**The single most important empirical result so far**: at CHRNA7, Boltz-2's pose-conditioned affinity rescues the α7 nAChR positive allosteric modulators that MAMMAL's sequence-only head buries in the bottom quartile. See [`design/PHASE_0_5_DECISION_RECORD.md`](design/PHASE_0_5_DECISION_RECORD.md).

## The hybrid architecture in one diagram

```
                    ┌────────────────────────┐
                    │   Input layer          │
                    │  - 22-target panel     │
                    │  - 298 compounds       │
                    └──────────┬─────────────┘
                               │
       ┌───────────────────────┼───────────────────────┐
       ▼                       ▼                       ▼
┌─────────────────┐ ┌─────────────────────┐ ┌───────────────────┐
│ CLUSTER A       │ │ CLUSTER B           │ │ CLUSTER C         │
│ Structure       │ │ ADMET / safety      │ │ Mechanism / KG    │
│ (Windows + WSL2)│ │ (Windows, CPU)      │ │ (WSL2 — sprint B) │
│                 │ │                     │ │                   │
│ ESM2-650M ──┐   │ │ ADMET-AI 41 EP      │ │ PrimeKG subgraph  │
│             ├── │ │  ├ BBB / hERG / DILI│ │ TxGNN zero-shot   │
│ MAMMAL DTI ─┘   │ │  ├ P-gp / CYP3A4    │ │  vs cognition     │
│      │          │ │  └ Ames / Caco-2    │ │  virtual anchor   │
│      ▼          │ │                     │ │  (MCI ∪ AD ∪ ADHD │
│ Boltz-2 affin.  │ │  HARD GATES         │ │   ∪ FXS ∪ narco)  │
│ (Boltzina mode) │ │  ↓ drops 60% before │ │                   │
│   ↳ WSL2 +      │ │     rank fusion     │ │                   │
│     cuEquivar.  │ │                     │ │                   │
└────────┬────────┘ └──────────┬──────────┘ └─────────┬─────────┘
         │                     │                      │
         └──────────────┬──────┴──────────────────────┘
                        ▼
            ┌────────────────────────┐
            │  RRF FUSION (k=60)     │
            │  4 rankers, per-target │
            │  calibrated weights    │
            └──────────┬─────────────┘
                       ▼
       ┌────────────────────────────────┐
       │  V3 OUTPUTS                    │
       │  - final_ranking_calibrated    │
       │  - provenance per candidate    │
       │  - disagreement archetypes     │
       │  - funnel narrative (top 20)   │
       │  - wet-lab shortlist           │
       └────────────────────────────────┘
```

## Quickstart paths (pick one)

### A. Run only the live pieces (~2 hrs cold cache, < 10 min warm)

```powershell
# One-time env (Windows + conda)
.\scripts\01_setup_env.ps1
conda activate mammal_env

# V1 + v2 live stack
python scripts/02_fetch_targets.py
python scripts/03_fetch_compounds.py
python scripts/04_score_dti.py
python scripts/05_sanity_check.py
python scripts/14_v2_cluster_b_admet.py
python scripts/15_v2_fusion.py
```

Outputs: `data/results/v2/admet_gates.parquet`, `rrf_ranking.parquet`, `final_ranking.parquet`, `reports/wet_lab_shortlist.md`.

### B. Add the structural cluster (overnight WSL2 sweep)

```powershell
# One-time WSL2 setup
wsl -d Ubuntu -u root -- bash scripts/_wsl2_setup_boltz.sh

# Launch overnight sweep (detached; runs ~10 h independent of your shell)
wsl -d Ubuntu -u root -- bash -c "source /root/mammal_env/bin/activate && nohup python /mnt/c/.../scripts/_wsl2_boltz_full_sweep.py > /tmp/wsl2_sweep.log 2>&1 &"

# Check progress any time
wsl -d Ubuntu -u root -- bash scripts/_wsl2_sweep_status.sh

# After ~10 hours, re-run fusion with the 3rd cluster
python scripts/15_v2_fusion.py
```

### C. Sprint: ChEMBL SQLite mirror + calibration + Cluster C + 4-cluster fusion

In flight; see [`design/V3_ATTACK_PLAN.md`](design/V3_ATTACK_PLAN.md) and the sprint phase docs (`design/PHASE_*.md`).

## Where to look for what

| Question | File |
|---|---|
| What's the architecture? | [`design/V2_HYBRID_ARCHITECTURE.md`](design/V2_HYBRID_ARCHITECTURE.md) |
| What's working / what's broken / what's the WSL2 plan? | [`design/V2_STATUS_AND_FORWARD_PLAN.md`](design/V2_STATUS_AND_FORWARD_PLAN.md) |
| What's the next 4 days of work? | [`design/V3_ATTACK_PLAN.md`](design/V3_ATTACK_PLAN.md) |
| Did the CHRNA7 rescue work? | [`design/PHASE_0_5_DECISION_RECORD.md`](design/PHASE_0_5_DECISION_RECORD.md) |
| Did WSL2 + cuEquivariance kernels work on RTX 5070? | [`design/WSL2_VALIDATION_RECORD.md`](design/WSL2_VALIDATION_RECORD.md) |
| Were the ChEMBL target IDs right? | [`design/T1_CHEMBL_AUDIT_VERDICT.md`](design/T1_CHEMBL_AUDIT_VERDICT.md) |
| Original research deep-dive | [`research/compass_artifact_wf-*.md`](research/) |
| Hybrid architecture research | [`research/Hybrid Architecture for MAMMAL-Based...md`](research/) |
| Deep-research on tier-4 issues | [`research/4-tier/`](research/4-tier/) (5 files) |

## Project layout

```
.
├── README.md                       (this file)
├── design/                         (architecture + decision records)
├── research/                       (background deep-dives, gitignored size)
├── configs/                        (thresholds.yaml, weights.yaml, anchor.yaml)
├── src/mammal_repurposing/
│   ├── cluster_a/                  (ESM2, Boltz-2, Boltzina)
│   ├── cluster_b/                  (ADMET-AI)
│   ├── cluster_c/                  (PrimeKG, TxGNN, cognition_anchor — sprint B)
│   ├── fetchers/                   (uniprot, pubchem, chembl, chembl_sqlite — sprint A)
│   ├── gates/                      (admet_gates with regulatory bypass)
│   ├── fusion/                     (rrf, lambdamart, calibration)
│   ├── provenance/                 (tracker, disagreement_report, narrative)
│   ├── analysis/                   (sanity, polypharm, composites, filters)
│   ├── scoring/                    (MAMMAL DTI, model_loader, runner)
│   └── _compat.py                  (Windows tiledbsoma stub)
├── scripts/
│   ├── 01-13_*.py / .ps1           (v1 + v2 pipeline, in order)
│   ├── 14_v2_cluster_b_admet.py    (ADMET-AI scoring)
│   ├── 15_v2_fusion.py             (RRF + provenance + narrative)
│   ├── 16_v2_esm2_embed.py         (ESM2 target cache)
│   ├── 17_v2_boltz_smoke.py        (single-pair Boltz validation)
│   ├── 18_v2_boltzina_sweep.py     (Windows full sweep, slow)
│   ├── 19_v2_allosteric_gate.py    (Phase 0.5 gate)
│   ├── 20_v3_audit_chembl_targets  (T1 audit)
│   └── _wsl2_*.{sh,py}             (WSL2-side scripts)
├── data/
│   ├── raw/                        (committed: targets_seed.csv, compounds_seed.csv, etc.)
│   ├── interim/                    (gitignored: fetched parquets)
│   ├── results/v2/                 (gitignored: scoring output)
│   ├── cache/                      (gitignored: per-pair memoization)
│   └── kg/primekg/                 (gitignored: ~2 GB download, sprint B)
├── reports/                        (committed: gate reports, audits, methodology)
├── notebooks/00_smoke_test.ipynb
└── tests/                          (respx-mocked fetchers, slow-marked scoring)
```

## Decision gates that matter

1. **Sanity gate** (v1, `scripts/05_sanity_check.py`) — peptide pollution + structural bias inflate negative controls; relaxed via the [`analysis/filters.py`](src/mammal_repurposing/analysis/filters.py) exclusion list
2. **ADMET hard gates** ([`gates/admet_gates.py`](src/mammal_repurposing/gates/admet_gates.py)) — physical kill criteria; regulatory bypass for already-approved drugs prevents over-aggressive filtering of clinically-validated CNS drugs
3. **Phase 0.5 CHRNA7 rescue** (`scripts/19_v2_allosteric_gate.py`) — ✅ passed; structural cluster justified
4. **Phase 3.1 calibration** (sprint Phase A, in flight) — per-target per-cluster Spearman ρ vs ChEMBL ground truth; output drives `configs/weights_calibrated.yaml`
5. **Cluster C positive-control rescue** (sprint Phase B gate) — donepezil, BPN14770, pitolisant must reach top decile by TxGNN indication score against the cognition virtual anchor

## Hardware

| Component | Spec |
|---|---|
| GPU | NVIDIA RTX 5070 (Blackwell sm_120, 12 GB VRAM) |
| Driver | 591.59 (or later) |
| OS | Windows 11 |
| WSL2 distro | Ubuntu 24.04 (for Boltz, PyG, future cellxgene-census) |
| Disk | ~50 GB (model weights ~10 GB, ChEMBL SQLite ~12 GB, PrimeKG ~2 GB, caches) |
| RAM | ≥16 GB (PrimeKG load needs 8+ GB during ingest) |

## License

Apache-2.0. MAMMAL (IBM Research), ADMET-AI (Stanford), Boltz-2 (MIT + Recursion), PrimeKG + TxGNN (Harvard) — all individually Apache-2.0 / MIT.

---

*Build status: V3 sprint in progress. Phase 0.5 ✅ passed; overnight WSL2 Boltz sweep running; ChEMBL SQLite mirror + per-target calibration starting now.*
