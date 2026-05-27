# MAMMAL Cognitive Enhancement Drug Repurposing

A four-cluster hybrid pipeline for drug repurposing on cognition-relevant targets, built around IBM Research's [MAMMAL](https://github.com/BiomedSciAI/biomed-multi-alignment) foundation model and augmented with structure-aware affinity (Boltz-2), comprehensive ADMET (ADMET-AI), and mechanism reasoning (PrimeKG + TxGNN). Runs on a single 12 GB consumer GPU.

## Honest scope

This pipeline does **not** discover "smart drugs." It enriches a candidate set so wet-lab cycles spend money on plausibility, not chemistry-lottery tickets. Roberts CA et al. (*Eur Neuropsychopharm* 2020, PMID 32709551) puts the effect-size ceiling for healthy-adult cognitive enhancement at SMD ≈ 0.21 for the strongest known agents (methylphenidate). The deliverable here is a calibrated, provenance-rich ranking + a publishable methodology contribution — not a miracle compound.

## Current state (v3 sprint Phases A, C, D, E shipped)

**Sprint commit log**: 3baf422 → ed4761c → 7467ac1 → 8f75a32 → 9d40e19 → 42b6597 → 1551a3e → 7c1d55e → b716ec7 → 9f800f8 → 5b415b5 → 5d7f67d (`git log --oneline`)

| Component | Status | Evidence |
|---|---|---|
| **Cluster A.1** — MAMMAL DTI head | ✅ live | 6,556 pairs scored; calibrated ρ vs ChEMBL per target (see below) |
| **Cluster A.2** — ESM2-650M target embeddings | ✅ live | 22 targets cached; paralog cos = 0.997 |
| **Cluster A.3** — Boltz-2 affinity (Windows) | ✅ live (slow) | ~149 s/pair, `--no_kernels` fallback |
| **Cluster A.3** — Boltz-2 affinity (WSL2 + cuEquiv kernels) | ✅ live (fast) | ~23 s/pair, 5.6× speedup confirmed |
| **Cluster B** — ADMET-AI 41 endpoints + hard gates | ✅ live | 53 PASS / 64 FLAG / 181 CUT of 298 compounds |
| **Cluster C.1** — PrimeKG loader (real, igraph) | ✅ code live | runs in WSL2 `txgnn_env` venv; download + setup scripted |
| **Cluster C.2** — TxGNN scorer (real) | ✅ code live | cognition virtual anchor (5 EFO union) wired |
| **ChEMBL 36 SQLite mirror** | ✅ live | local mirror via `chembl_downloader`; A.5 PASS (19/20 vs REST) |
| **Phase A.7 Calibration linchpin** | ✅ **shipped** | `reports/calibration_report.md`, `configs/weights_calibrated.yaml` |
| **Phase C 4-cluster RRF (calibrated + uncalibrated)** | ✅ live | `scripts/15_v2_fusion.py --out-suffix _{un,}calibrated` |
| **Phase D calibration diff** | ✅ shipped | `reports/fusion_calibration_diff.md` (Spearman ρ = +0.994) |
| **Phase E methodology note v1** | ✅ shipped | `reports/methodology_v1.md` |
| **Wet-lab shortlist v3 (4-cluster scorecards)** | ✅ shipped | `reports/wet_lab_shortlist_v3.md` |
| **Phase 0.5 CHRNA7 rescue gate** | ✅ **PASSED** | TC-5619 100%, encenicline 80% (vs v1 19%, 7%) |
| **Phase A.4 full backstop (6,556 pairs)** | 🔄 running | SQLite, ~75 min ETA; `data/results/chembl_evidence.parquet` |
| **Phase 0.4 full Boltz overnight sweep** | 🔄 **running** in WSL2 | 1,165 pairs ETA ~24 h; check via `scripts/_wsl2_sweep_status.sh` |
| **Phase B PrimeKG+TxGNN run** | ⏳ blocked | waiting for overnight Boltz sweep to free the GPU |

**The single most important empirical result so far**: at CHRNA7, Boltz-2's pose-conditioned affinity rescues the α7 nAChR positive allosteric modulators that MAMMAL's sequence-only head buries in the bottom quartile. See [`design/PHASE_0_5_DECISION_RECORD.md`](design/PHASE_0_5_DECISION_RECORD.md).

**The single most important Phase A.7 finding**: only **2 of 22 panel targets** (DRD1 ρ=+0.31, HCRTR1 ρ=+0.37) have MAMMAL ρ ≥ 0.30 against ChEMBL pchembl. **4 are inverted** (SLC6A3/DAT ρ=-0.71, SLC6A2/NET ρ=-0.53, GRIN2A ρ=-0.35, GRIN2B ρ=-0.30). The remaining 14 are weakly informative (0 ≤ ρ < 0.30). See [`reports/calibration_report.md`](reports/calibration_report.md) and the limitations section of [`reports/methodology_v1.md`](reports/methodology_v1.md).

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

### C. Full sprint reproduction (ChEMBL SQLite + calibration + 4-cluster fusion)

```powershell
# One-time: download ChEMBL 36 SQLite (~4 GB compressed, ~13 GB extracted, ~15 min)
python -c "import chembl_downloader; chembl_downloader.download_extract_sqlite()"

# Phase A.5 — verify SQLite vs REST agreement (gate: must PASS)
python scripts/_v3_sqlite_vs_rest_smoke.py

# Phase A.6 — clean T1 audit
python scripts/24_v3_audit_chembl_targets_sqlite.py

# Phase A.7 — calibration linchpin (THE unlock)
python scripts/22_v3_calibration.py

# Phase A.4 — full backstop (~75 min for 6,556 pairs)
python scripts/21_v3_chembl_evidence_sqlite.py --all-pairs

# Phase C — both fusion passes
python scripts/15_v2_fusion.py --calibrated-weights NONE --out-suffix _uncalibrated
python scripts/15_v2_fusion.py --out-suffix _calibrated

# Phase D — diff calibrated vs uncalibrated
python scripts/25_v3_fusion_diff.py

# Phase C completion — v3 wet-lab shortlist with 4-cluster scorecards
python scripts/26_v3_wet_lab_shortlist.py
```

See also: [`design/V3_ATTACK_PLAN.md`](design/V3_ATTACK_PLAN.md) and the sprint phase docs (`design/PHASE_*.md`).

## Where to look for what

| Question | File |
|---|---|
| What does this pipeline DO and NOT do? | [`reports/methodology_v2.md`](reports/methodology_v2.md) — V4+V5 coherent narrative (supersedes v1) |
| **What's the v4 state + forward plan?** | [`design/V4_STATUS_AND_FORWARD_PLAN.md`](design/V4_STATUS_AND_FORWARD_PLAN.md) |
| What was the v3 state (historical)? | [`design/V3_STATUS_AND_FORWARD_PLAN.md`](design/V3_STATUS_AND_FORWARD_PLAN.md) |
| What's the (v2) architecture skeleton? | [`design/V2_HYBRID_ARCHITECTURE.md`](design/V2_HYBRID_ARCHITECTURE.md) |
| What was working/broken at v2 (historical)? | [`design/V2_STATUS_AND_FORWARD_PLAN.md`](design/V2_STATUS_AND_FORWARD_PLAN.md) |
| What was the next 4 days of v3 sprint work? | [`design/V3_ATTACK_PLAN.md`](design/V3_ATTACK_PLAN.md) |
| Did the CHRNA7 rescue work? | [`design/PHASE_0_5_DECISION_RECORD.md`](design/PHASE_0_5_DECISION_RECORD.md) |
| Did WSL2 + cuEquivariance kernels work on RTX 5070? | [`design/WSL2_VALIDATION_RECORD.md`](design/WSL2_VALIDATION_RECORD.md) |
| Were the ChEMBL target IDs right? | [`design/T1_CHEMBL_AUDIT_VERDICT.md`](design/T1_CHEMBL_AUDIT_VERDICT.md), [`reports/chembl_target_id_audit_sqlite.md`](reports/chembl_target_id_audit_sqlite.md) |
| Per-target calibration linchpin (V3) | [`reports/calibration_report.md`](reports/calibration_report.md) |
| **THE Tanimoto-beats-MAMMAL breakthrough (V4)** | [`reports/tanimoto_baseline_v1.md`](reports/tanimoto_baseline_v1.md) |
| Diagnostic protocol (prior-collapse + power + 5 lateral) | [`reports/diagnostics_v1.md`](reports/diagnostics_v1.md) |
| §7.11 isotonic per-target calibration (V4) | [`reports/calibration_comparison_v1.md`](reports/calibration_comparison_v1.md) |
| §7.5 pocket-conditioned MVP (13/13 gates) | [`reports/pocket_database_v1.md`](reports/pocket_database_v1.md) |
| §8.1 multi-class faceted shortlist (V4 deliverable) | [`reports/wet_lab_shortlist_v4_faceted.md`](reports/wet_lab_shortlist_v4_faceted.md) |
| **§8.0b-zn liability audit (V5 transition)** | [`reports/liability_audit_v1.md`](reports/liability_audit_v1.md) (z-norm), [`reports/liability_audit_v1_absolute.md`](reports/liability_audit_v1_absolute.md) (calibration-failure evidence) |
| **§8.15 Tanimoto-vs-MAMMAL disagreement signal (V5 transition)** | [`reports/disagreement_signal_v1.md`](reports/disagreement_signal_v1.md) |
| **THE V5 production wet-lab handoff** | [`reports/wet_lab_shortlist_v6_full.md`](reports/wet_lab_shortlist_v6_full.md) — calibrated MAMMAL + Z-norm + §8.0b-zn + faceted; 43 PASS / 60 FLAG / 195 CUT |
| **§8.10 nootropic-similarity annotation** (V5+ sprint) | [`reports/nootropic_similarity_v1.md`](reports/nootropic_similarity_v1.md) |
| **§8.13 pocket-conditioned liability composition** (V5+ sprint) | [`reports/liability_pocket_aware_v1.md`](reports/liability_pocket_aware_v1.md) |
| **§8.3 ClinicalTrials.gov IP-status cross-reference** (V5+ sprint) | [`reports/clinical_trials_v1.md`](reports/clinical_trials_v1.md) |
| **§8.16 calibrator round-trip QC** (V5+ sprint) | [`reports/calibrator_qc_v1.md`](reports/calibrator_qc_v1.md) |
| **§7.4 v2 selectivity 4-metrics (Gini + S10× + Entropy + PI)** | [`reports/selectivity_v6_tanimoto_4metrics.md`](reports/selectivity_v6_tanimoto_4metrics.md) |
| **§14 Hypothesis Validation Ledger** (Tier-3 sprint) — falsifiable claim audit | [`reports/hypothesis_audit_v1.md`](reports/hypothesis_audit_v1.md) |
| **§8.0a Pareto NSGA-III 5-axis shortlist** (Tier-3 sprint) | [`reports/pareto_ranking_v1.md`](reports/pareto_ranking_v1.md) |
| **§7.12 Conformal prediction per-target** (Tier-3 sprint) | [`reports/conformal_calibration_v1.md`](reports/conformal_calibration_v1.md) |
| **§7.13 Scaffold-aware AL re-ranking** (Tier-3 sprint) | [`reports/scaffold_aware_v1.md`](reports/scaffold_aware_v1.md) |
| **§8.6 Brain-region annotation** (Tier-3 sprint, V6 Cluster D preview) | [`reports/brain_region_v1.md`](reports/brain_region_v1.md) |
| **LambdaMART meta-ranker** (Tier-3b sprint, NDCG@25 = 0.891 vs 0.774 baseline) | [`reports/lambdamart_meta_v1.md`](reports/lambdamart_meta_v1.md) |
| **§7.15 Hierarchical Bayes (GRIN pool)** (Tier-3b sprint) | [`reports/hierarchical_bayes_v1.md`](reports/hierarchical_bayes_v1.md) |
| **§8.14 Pocket-routed isotonic** (Tier-3b sprint, S1 vs S2 demo) | [`reports/pocket_routed_calibration_v1.md`](reports/pocket_routed_calibration_v1.md) |
| **§8.2 DrugComb combinations** (Tier-3b sprint, fallback table) | [`reports/drugcomb_combinations_v1.md`](reports/drugcomb_combinations_v1.md) |
| **V6 architecture plan** (Multi Head DTI + Cluster D) | [`design/V6_ARCHITECTURE_PLAN.md`](design/V6_ARCHITECTURE_PLAN.md) |
| What did calibration change in the rankings? | [`reports/fusion_calibration_diff.md`](reports/fusion_calibration_diff.md), [`reports/fusion_tanimoto_addition_diff.md`](reports/fusion_tanimoto_addition_diff.md) |
| Phase A.5 SQLite agreement gate | [`reports/sqlite_vs_rest_smoke.md`](reports/sqlite_vs_rest_smoke.md) |
| Original research deep-dive | [`research/compass_artifact_wf-*.md`](research/) |
| Hybrid architecture research | [`research/Hybrid Architecture for MAMMAL-Based...md`](research/) |
| Deep-research on tier-4 issues (7 V3+/V5+ docs) | [`research/4-tier/`](research/4-tier/) — includes `Multi Head DTI.md` (V5 priority) and `Multi-Source Neurobiological Prior...md` (V6 priority); 5 shipped docs are in `archived/` |

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
