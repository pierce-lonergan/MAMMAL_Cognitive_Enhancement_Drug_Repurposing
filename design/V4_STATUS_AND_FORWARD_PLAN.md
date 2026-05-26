# V4 — Status, Comprehensive Architecture, Forward Plan

**Read this first.** This is the source-of-truth status document for the project after **four post-V3 breakthroughs** (diagnostics + Tanimoto; selectivity + faceted shortlist; isotonic calibration; pocket-conditioned MVP) plus the liability-panel infrastructure shipped. It supersedes [V3_STATUS_AND_FORWARD_PLAN.md](V3_STATUS_AND_FORWARD_PLAN.md) as the live state-of-the-world; v3 remains as the historical anchor.

**Snapshot**: `main @ f188b82` — 5 cluster modules + 4 new infrastructure packages (diagnostics / calibration / selectivity / pockets) + 1 new gate (§8.0b liability) live. WSL2 Boltz overnight sweep at **79.5% (926/1165)**, ETA ~4.5h.

> **Why a v4 doc** — V3 closed with one calibration linchpin and a 10-item roadmap. We then shipped 4 of the top-10 priorities (3 wins beat their pre-committed predictions, the 4th is queued). The architecture is meaningfully different now: 5 parallel signal sources (MAMMAL DTI, ESM2, Boltzina, Tanimoto-to-actives, ADMET-AI), 3 honest calibration layers (sign-flipped Phase A.7 + Spearman ρ + isotonic per-target), a faceted top-N output, a pocket-classifier provenance column, and a liability gate ready to deploy. The roadmap deserves a clean reset.

---

## 1. Executive Summary

The v3 → v4 transition was driven by **four research deep-dives** (`research/4-tier/`) each commissioned in response to a v3 finding:

| Research doc | Triggered by | Shipped artifact | Outcome vs prediction |
|---|---|---|---|
| `Diagnosing MAMMAL DTI Anti-Correlation.md` | v3 Phase A.7 found ρ ≤ -0.30 at 4 targets | `diagnostics/` package + Tanimoto baseline | Identified MAMMAL prior-collapse is **panel-wide** (19/22 SEVERE); a Tanimoto-on-Morgan-FP baseline beats MAMMAL at **every** audited target (SLC6A3 +0.90 vs -0.70; DRD1 +0.85 vs +0.29) |
| `Graczyk-Selectivity-Faceted-Shortlist.md` | v3 had HRH3 dominating 23/25 top compounds | `selectivity/` + `fusion/faceted_shortlist.py` | Lock-in **dissolved**; top-5 per mechanism class with cross-facet provenance |
| `Cognition-44Target-Liability-Panel.md` | ADMET-AI missing 5-HT2B / HRH1+M1 anticholinergic combos | `gates/liability_panel.py` + 44-target seed | Stage 1 complete; queued for ~1hr MAMMAL re-run when GPU frees |
| `Isotonic-PerTarget-Calibration.md` | v3 had awkward `weight=0.30` hack at INVERTED targets | `calibration/` package | SLC6A3 +0.62 (predicted [+0.45, +0.65]), SLC6A2 +0.40 (predicted +0.40); both Tier A/B |
| `Pocket-Conditioned-Boltz2.md` | v3 had no orthosteric/allosteric provenance | `pockets/` package + curated DB | 13/13 validation gates pass; PDE4D BPN14770 (allosteric, UCR2) correctly distinguished from rolipram (catalytic) |

**Component status matrix** (delta from v3 in **bold**):

| Component | Status | Evidence at this snapshot |
|---|---|---|
| Cluster A.1 — MAMMAL DTI | ✅ live | 6,556 raw pairs; **+ 6,556 calibrated pairs** via isotonic |
| Cluster A.2 — ESM2-650M cache | ✅ live | 22 targets cached |
| Cluster A.3 — Boltz-2 / Boltzina (WSL2 + cuEquiv) | ✅ live | **~926/1165 sweep done; affinity scalars only** |
| **Cluster A.4 — Tanimoto-to-actives ranker** | **✅ NEW (v4)** | `cluster_a/tanimoto_ranker.py`; weight 1.5 in `weights.yaml` |
| Cluster B — ADMET-AI 41 EP + hard gates | ✅ live | 53 PASS / 64 FLAG / 181 CUT of 298 |
| **Cluster B' — 44-target off-target liability panel** | **🟡 INFRASTRUCTURE READY** | UniProt enrichment shipped; needs ~1hr MAMMAL re-run |
| Cluster C — PrimeKG + TxGNN | 🟡 CODE LIVE, NOT RUN | `txgnn_env` venv ready; blocked on overnight Boltz sweep |
| ChEMBL 36 SQLite mirror | ✅ live | A.5 PASS (19/20); 99-min full backstop |
| Phase A.7 calibration linchpin | ✅ shipped | `reports/calibration_report.md`; per-target Spearman ρ |
| **§7.11 isotonic per-target calibration** | **✅ NEW (v4)** | SLC6A3 +0.62 Tier A, SLC6A2 +0.40 Tier B |
| **§7.5 pocket-conditioned classifier (MVP)** | **✅ NEW (v4)** | 7 priority targets, 13 pockets, 13/13 gates pass |
| **Selectivity layer (Graczyk Gini + S(10×))** | **✅ NEW (v4)** | `selectivity/` package; ships with Tanimoto vector |
| **Multi-class faceted shortlist** | **✅ NEW (v4)** | 8 mechanism + 9 targeted-pair facets, cross-facet provenance |
| **Diagnostics protocol (Tier-1)** | **✅ NEW (v4)** | `diagnostics/` 9 modules; prior_collapse + power_analysis + 5 lateral |
| Phase C 4-cluster RRF (calibrated + uncalibrated) | ✅ live | both passes ship; +1 cluster (Tanimoto) → 5-cluster fusion |
| Phase D calibration diff | ✅ shipped | Spearman ρ = +0.994 |
| Phase E methodology note v1 | ✅ shipped | + 2 post-ship update callouts for v4 |
| **Wet-lab shortlist v4 (faceted)** | **✅ NEW (v4)** | `reports/wet_lab_shortlist_v4_faceted.md` |
| ChEMBL evidence backstop (Phase A.4) | ✅ complete | 6,556 rows; 275 CORROBORATED / 6,257 NOVEL |
| Phase 0.5 CHRNA7 rescue gate | ✅ PASSED | TC-5619 100%, encenicline 80% (vs v1 19%, 7%) |
| Phase A.4 backstop | ✅ COMPLETE | 99 min for 6,556 rows |
| Phase 0.4 full Boltz overnight sweep | 🔄 79.5% | 926/1165, ~4.5h ETA |
| Phase B PrimeKG+TxGNN run | ⏳ blocked | waiting for overnight sweep to free GPU |
| **Pose extraction from Boltz outputs** | **⏳ NEXT** | needed to operationalise §7.5 on the full grid |

**The single most important breakthrough since v3** is the discovery that **a 1996-vintage Tanimoto-on-Morgan-FP baseline beats the 458M-parameter foundation model at every audited cognition target**. This is the kind of negative finding that *defines* the methodology contribution — it falsified the v3 plan's #1 priority (LoRA fine-tune) and replaced it with a 1-day fix (Tanimoto ranker as Cluster A.4). The v4 pipeline is now an ensemble where the cheap-but-honest signal source out-ranks the expensive black box, with the foundation model as a (recalibrated) supporting voter.

**The single most important pending result** continues to be the overnight Boltz sweep finishing — that unblocks **three** queued workstreams in parallel: (i) Phase A.7 re-run with full Boltz coverage to test whether the inverted SLC6A3/SLC6A2 calibrators *also* improve under structure-aware affinity, (ii) Phase B Cluster C run (PrimeKG + TxGNN), and (iii) the 44-target liability MAMMAL re-run.

---

## 2. What the Four Breakthroughs Actually Added

This section is the v4 equivalent of v3's "what the calibration linchpin exposed." Each subsection reads like a mini-postmortem: prediction → outcome → consequence.

### 2.1 Breakthrough #1 — Tanimoto-to-actives baseline beats MAMMAL panel-wide (commit `530dc40`)

**Prediction (V3 §7.1 diagnostic)**: distinguish manifold mismatch from rank-collapse to pick between cross-DTI ensemble (§7.7) and LoRA fine-tune (§7.6).

**Outcome**: the diagnostic protocol revealed two findings simultaneously:

1. **Panel-wide prior collapse**: MAMMAL predictions cluster at the training prior `norm_y_mean = 5.79` with std 0.08-0.18 on a training SD of 1.34 — a 7-45× dynamic-range collapse at *every* target including the "STRONG" controls (DRD1, HCRTR1). Per-target Spearman ρ values from Phase A.7 are computed on essentially constant predictions.

2. **Statistical-power blockade**: only SLC6A3 (n=26) and SLC6A2 (n=25) have Bonett-Wright CI excluding zero. DRD1 ρ=+0.31 (n=21) has CI **[-0.15, +0.66]** — not distinguishable from zero. GRIN2A (n=8) and GRIN2B (n=14) are statistically indistinguishable from random.

**Then the Tanimoto baseline experiment** (`scripts/31_v3_tanimoto_baseline.py`):

| Target | MAMMAL ρ | Tanimoto ρ | Δρ |
|---|---|---|---|
| SLC6A3 | -0.70 | **+0.90** | +1.59 |
| SLC6A2 | -0.60 | **+0.91** | +1.51 |
| GRIN2A | -0.40 | +0.76 | +1.16 |
| GRIN2B | -0.30 | +0.82 | +1.12 |
| ACHE | +0.24 | +0.81 | +0.57 |
| DRD1 | +0.29 | +0.85 | +0.56 |
| HCRTR1 | +0.37 | +0.78 | +0.41 |

**Tanimoto wins 7/0**. A 1996-vintage Morgan-FP-on-ChEMBL-pchembl≥8-actives baseline beats the 458M-param foundation model at every single target. **The panel contains signal MAMMAL is destroying via prior collapse.**

**Consequence**:
- V3 §7.6 (LoRA fine-tune MAMMAL) is **not the right next move** — the prediction surface is degenerate; fine-tuning can't add resolution to a collapsed manifold.
- V3 §7.7 (cross-DTI ensemble) is the right move — MMAtt-DTA's published transporter ρ > 0.72 is in the same league as our Tanimoto +0.90.
- **The Tanimoto baseline is the cheapest possible ensemble member** — we shipped it immediately as Cluster A.4 with weight 1.5 in `configs/weights.yaml`. Donepezil leapt from rank #42 to #4 in the 5-cluster fusion.

### 2.2 Breakthrough #2 — Selectivity + faceted shortlist dissolves HRH3 lock-in (commit `1c288a8`)

**Prediction (V3 §7.4 + §8.1)**: introduce Graczyk Gini + S(10×) + multi-class faceted top-5 to break the v3 phenomenon where 23/25 top compounds had `mammal_best_target = HRH3`.

**Outcome**: the prediction held, but only after a methodologically interesting failure:

1. **Gini on MAMMAL pKd is degenerate** — all 298 compounds → `flat` category, Gini ~0.08 everywhere. This is the *direct downstream consequence* of breakthrough #1's prior-collapse finding. Validation gate G1 (positive controls donepezil/pitolisant/BPN14770 Gini ≥ 0.7) failed catastrophically.

2. **Gini on the Tanimoto-to-actives vector works correctly**:
   - donepezil → ACHE (top_target correct)
   - pitolisant → HRH3 ✓
   - BPN14770 → PDE4D ✓
   - Top-10 spans 7 mechanism classes (was 1 class in v3)

3. **Faceted shortlist** (`scripts/28_v3_faceted_shortlist.py`) — 8 mechanism class facets + 9 targeted-pair facets, top-5 each, **+ cross-facet provenance**. Gates G3-G6 audit:
   - G5 (pitolisant in histaminergic top-5) ✅
   - G6 (pitolisant in HRH3 but NOT HRH3+DRD1 pair, hygiene check) ✅
   - G4 partial (donepezil in cholinergic at #5, galantamine missing due to scaffold novelty)
   - G3 fails (TC-5619 / encenicline have novel α7-PAM scaffolds Tanimoto under-rates)

**Consequence**:
- The faceted shortlist is the new wet-lab handoff format (`reports/wet_lab_shortlist_v4_faceted.md`).
- Pitolisant correctly anchors histaminergic facet at #2; atomoxetine/duloxetine/fluoxetine/clonidine dominate noradrenergic; ampakines + AMPA agonists dominate glutamatergic AMPA — each facet contains the canonical compounds for its class.
- Cross-facet provenance prevents triple-counting (lurasidone appears in 10 facets, but is one compound — the count column makes this explicit).
- **Tanimoto-on-Morgan-FP under-rates novel scaffolds** is a known deployment caveat that becomes the motivation for §7.7 cross-DTI ensemble (a contrastive-trained DTI like PSICHIC may correct this).

### 2.3 Breakthrough #3 — Isotonic per-target calibration repairs INVERTED transporters (commit `8624fd1`)

**Prediction (V3 §7.11 from `Isotonic-PerTarget-Calibration.md`)**: `IsotonicRegression(increasing='auto')` naturally absorbs sign inversion at MAMMAL_ONLY_INVERTED targets. Predicted SLC6A3 post-cal ρ ∈ [+0.45, +0.65]; SLC6A2 ∈ [+0.30, +0.55]. GRIN2A/2B confirmed Scenario 3 (structural blindness at ifenprodil dimer interface) — calibration ceiling ~+0.2.

**Outcome** (per `reports/calibration_comparison_v1.md`):

| Target | n | Raw ρ | Post-cal ρ | Δρ | CI | Tier |
|---|---|---|---|---|---|---|
| **SLC6A3 (DAT)** | 23 | -0.70 | **+0.62** | +1.32 | [+0.71, +0.80] | **A** |
| **SLC6A2 (NET)** | 21 | -0.60 | **+0.40** | +0.99 | [+0.60, +0.72] | **B** |
| GRIN2B | 13 | -0.30 | -0.17 | +0.13 | (CI spans 0) | C |
| GRIN2A | 7 | -0.40 | +0.11 | +0.50 | (n too small) | C |

SLC6A3 landed at the **high end** of the predicted range; SLC6A2 hit the prediction exactly. GRIN2A/2B confirmed as Scenario 3 (escalation to §7.7).

**Consequence**:
- Production deployment via `scripts/33_v3_apply_calibration.py` writes `data/results/dti_scores_calibrated.parquet` (6,556 rows, both raw + calibrated columns). 18 isotonic calibrators applied; 4 targets pass through.
- SLC6A3 dynamic range expanded **20×** (std 0.084 → 1.672) — the prior-collapse degeneracy is repaired at the SLC6A3-specific calibrator output.
- Decision router (`calibration/router.py`) implements the §1D decision matrix (n-bucket × |ρ| sign) + Tier A/B/C/D post-fit classification. `data/calibration/router_decisions.csv` is the single audit trail.
- **Beta-calibration deferred**: the `betacal` PyPI package is internally a binary classifier (wraps sklearn `LogisticRegression`) — the research doc's continuous-regression example doesn't actually run. Documented limitation.
- **Hierarchical Bayesian (PyMC) deferred to v2** — 3-5 day budget for GRIN family pool.

### 2.4 Breakthrough #4 — Pocket-conditioned MVP, 13/13 validation gates pass (commit `3884ba5`)

**Prediction (V3 §7.5 from `Pocket-Conditioned-Boltz2.md`)**: 4-detector consensus (P2Rank + PocketMiner + CryptoBench + Boltz cofold) for a curated 7-priority-target centroid DB. Headline target: PDE4D BPN14770 (UCR2 allosteric) vs rolipram (catalytic orthosteric) discrimination.

**Outcome — shipped the MVP (curated DB + geometric classifier; detector ensemble deferred to Sprint 2)**:

| Gate | Targets / Compounds | Pass rate |
|---|---|---|
| **P1** (orthosteric positive controls) | HRH3/PF-03654746, ACHE/donepezil, CHRNA7/epibatidine, DRD1/SKF-81297 | **4/4** ✅ |
| **P2** (allosteric_known positive controls) | **PDE4D/BPN14770** ★, DRD1/LY3154207 | **2/2** ✅ |
| **P3** (negative controls, 50Å displacement → surface_artifact) | all 7 targets | **7/7** ✅ |

★ **PDE4D BPN14770** → `allosteric_known` is the headline win — pipeline can now distinguish UCR2-closing NAMs from catalytic orthosteric inhibitors at the geometric level, enabling §8.0b emesis-liability gating split between the two pharmacologies.

**Key design decision: ligand-anchored centroids**. The research doc flagged that residue-list centroids can mis-match real ligand positions due to PDB-specific numbering or alternate conformations. We added an optional `ligand_anchor: <HET_3_letter>` mode that uses a co-crystallised ligand's heavy-atom mean as the centroid — guarantees coordinate-frame agreement. This was decisive for PDE4D where 6NJJ uses catalytic-domain-only numbering not mappable to canonical full-length PDE4D positions.

**ACHE donepezil dual-site detection** correctly fires: when a pose contacts both CAS (orthosteric) and PAS (allosteric_known) within 8 Å, the classifier returns `orthosteric` per research doc §2G (donepezil-mode span is the orthosteric pharmacology).

**Consequence**:
- 12 centroids cached at `data/pockets/centroids/<target>.json`; 11 PDBs cached at `data/pockets/pdbs/`.
- Classifier API: `classify_pose(pose_xyz, target_gene, db) → PocketClassification` with `rank_multiplier`, `manual_review`, `is_dual_site`.
- **Not yet operational on the live grid**: Boltz overnight sweep saves affinity scalars only, not poses. Next-session item: extend `_wsl2_boltz_full_sweep.py` to save mmCIF poses → enables pocket classification on the full 1,165-pair sweep output.
- Sprint 2 (detector ensemble P2Rank + PocketMiner + CryptoBench) and Sprint 3 (§8.0b liability gating split by pocket_class, §8.1 new facets: PDE4D-allosteric, CHRNA7-type-I, CHRNA7-type-II, DRD1-PAM, ACHE-dual-site) deferred.

### 2.5 What v3 said vs what we delivered

| V3 forward-plan item | V4 status | Notes |
|---|---|---|
| §7.1 per-target failure-mode diagnostic | ✅ shipped | `diagnostics/` 9 modules; revealed panel-wide prior collapse |
| §7.4 selectivity panel scoring | ✅ shipped | Graczyk Gini + S(10×) + categorisation |
| §7.5 pocket-conditioned Boltz | ✅ MVP shipped | 13/13 gates; pose-extraction next |
| §7.6 LoRA fine-tune MAMMAL | ❌ **falsified** | Tanimoto already beats it; deferred / cancelled |
| §7.7 cross-DTI ensemble (MMAtt-DTA) | ⏳ next sprint | MMAtt-DTA reports ρ > 0.72 on transporters — needs to beat Tanimoto +0.90 |
| §7.10 LLM-agent prioritisation overlay | ⏳ Tier 4 (defer) | |
| §7.11 isotonic per-target calibration | ✅ shipped | SLC6A3 +0.62 Tier A, SLC6A2 +0.40 Tier B |
| §7.12 conformal prediction per-target | ⏳ Tier 3 | superseded in priority by isotonic |
| §7.13 scaffold-aware active learning | ⏳ Tier 3 | now informed by Tanimoto-baseline scaffold limits |
| §7.14 residual-correction XGBoost meta-ranker | ⏳ Tier 3 | optional alternative to isotonic; not urgent |
| §8.0a Pareto NSGA-III restructure | ⏳ Tier 3 | |
| §8.0b 44-target liability panel | 🟡 ready, GPU-queued | Stage 1 done; ~1hr DTI re-run pending |
| §8.0c Cluster D neurobio prior (AHBA + OT Genetics) | ⏳ Tier 3 | |
| §8.1 multi-class faceted shortlist | ✅ shipped | 8 mechanism + 9 targeted-pair facets |
| §8.5 off-target liability beyond ADMET | ⊆ §8.0b superset | folded into liability panel |
| §8.8 cross-DTI ensemble (DeepPurpose) | ⊆ §7.7 superset | MMAtt-DTA is the better entry; same architectural slot |

**Score**: 6 of 16 v3 priorities shipped (37.5%), 1 falsified, 9 still queued. Three of the six were *better than predicted* (SLC6A3 post-cal ρ, faceted shortlist diversification, P1/P2 gate pass-rate).

---

## 3. What Works Today (end-to-end runnable)

```
                       Windows (mammal_env)                              Status
v1 pipeline: 02_fetch_targets → 13_wet_lab_shortlist                 ✅ runnable
v2 ADMET cluster: 14_v2_cluster_b_admet → admet_gates.parquet         ✅ runs in ~1 min
v2 fusion (BOTH passes): 15_v2_fusion --out-suffix _{un,}calibrated   ✅ runs in seconds
v2 fusion + Tanimoto ranker: --add-tanimoto-ranker                    ✅ v4 (530dc40)
v2 ESM2: 16_v2_esm2_embed → 22 cached .pt files                        ✅
v2 Boltz: 17 smoke / 18 sweep / 19 gate / _boltzina_focused           ✅

== V3 SQLite-backed Phase A ==
A.5 smoke: _v3_sqlite_vs_rest_smoke.py                                 ✅ 19/20 PASS
A.6 audit: 24_v3_audit_chembl_targets_sqlite.py                        ✅ 21 ALIGNED / 1 NO_CURRENT
A.7 calibration: 22_v3_calibration.py                                  ✅ refined verdict matrix
A.4 backstop: 21_v3_chembl_evidence_sqlite.py --all-pairs              ✅ 99 min for 6,556 rows
Phase D: 25_v3_fusion_diff.py                                          ✅ Spearman ρ = +0.994
Wet-lab v3: 26_v3_wet_lab_shortlist.py (4-cluster scorecards)          ✅

== V4 NEW pipelines ==
Diagnostics: 30_v3_diagnose_inverted.py (7 modules)                    ✅ prior_collapse + power + lateral
Tanimoto baseline: 31_v3_tanimoto_baseline.py                          ✅ 7/0 vs MAMMAL
Selectivity: 27_v3_selectivity_scoring.py --use-tanimoto               ✅ Graczyk Gini + S(10×)
Faceted shortlist: 28_v3_faceted_shortlist.py                          ✅ 8+9 facets w/ provenance
Isotonic calibration sweep: 32_v3_calibration_comparison.py            ✅ SLC6A3 +0.62, SLC6A2 +0.40
Isotonic apply: 33_v3_apply_calibration.py                             ✅ dti_scores_calibrated.parquet
Pocket DB build: 34_v3_build_pocket_database.py                        ✅ 13/13 gates pass
Liability panel (Stage 1): 29_v3_liability_panel.py --enrich-only      ✅ 44 UniProt sequences cached


                       WSL2 Ubuntu (mammal_env)                          Status
Overnight Boltz sweep: _wsl2_boltz_full_sweep.py                      🔄 926/1165 (~79.5%, ETA 4.5h)
Sweep status: _wsl2_sweep_status.sh                                    ✅ available


                       WSL2 Ubuntu (txgnn_env, separate venv)            Status
Cluster C orchestrator: 23_v3_cluster_c.py                            🟡 code ready, not run
PrimeKG download: _wsl2_download_primekg.sh                            ✅ scripted
Isolated venv setup: _wsl2_setup_cluster_c.sh                          🟡 ready


                       Queued for next session (GPU-blocked)             Status
Liability panel Stages 2+3: 29_v3_liability_panel.py (no --skip-dti)  ⏳ ~1hr MAMMAL re-run
Phase A.7 + C + D refresh after Boltz sweep completes                  ⏳ tests INVERTED rescue under structure
Pose-saving Boltz wrapper                                              ⏳ unlocks §7.5 on full grid
```

**Reproducibility**: every artifact lands in `data/results/` or `data/results/v2/`. Configs at `configs/{thresholds,weights,weights_calibrated}.yaml`. Calibrators at `data/calibration/isotonic/<uniprot>.pkl`. Pocket centroids at `data/pockets/centroids/<target>.json`. Reports at `reports/`. Sprint history = git log; methodology at `reports/methodology_v1.md` (post-ship update callouts for v4 breakthroughs).

---

## 4. What's Broken / Limiting Now (post-v4 state)

The V3 doc listed 9 items. V4 has resolved 4 of those (§7.11 isotonic, §8.0b liability infrastructure, §7.4 selectivity, §8.1 faceted shortlist). The current limiting items are different.

### 4.1 Boltz pose extraction not wired — §7.5 not operational on the live grid

**Symptom**: `_wsl2_boltz_full_sweep.py` produces affinity scalars (`affinity_pred_value`, `affinity_probability_binary`, `pose_plddt`) but discards the per-pair mmCIF pose structures. `data/results/v2/boltzina_affinity.parquet` has no XYZ columns. §7.5 pocket classifier requires `pose_centroid_xyz` to classify.

**Severity**: HIGH — the pocket-conditioning MVP is shipped and validated (13/13 gates) but cannot annotate the 1,165-pair sweep output.

**Remediation**: extend the sweep launcher to either
1. Save mmCIF poses to `data/results/v2/boltz_poses/<target>_<compound>.cif` (storage: ~50KB/pair × 1,165 = ~60 MB), OR
2. Extract heavy-atom centroid in-line and write a new `pose_centroid_x/y/z` column to the parquet (compact; no separate file)

Option 2 is preferred — keeps the parquet self-contained. ~1-day implementation; needs a re-run of the Boltz sweep to populate (or a separate Boltz-pose-only re-run on the existing 926 done pairs).

### 4.2 Liability panel MAMMAL re-run queued, not executed

**Symptom**: `scripts/29_v3_liability_panel.py --enrich-only` shipped 44-target UniProt sequence enrichment. Stages 2-3 (MAMMAL DTI on 298 × 44 = 13,112 pairs, then gate application) need the GPU. Both Boltz (currently using GPU) and MAMMAL want the RTX 5070's full ~12 GB VRAM; they can't co-exist.

**Severity**: HIGH — the §8.0b liability panel is the discrimination that ADMET-AI cannot deliver. Expected result: 2 hard CUTs (aripiprazole, amitriptyline) + ~7 FLAGs that ADMET alone misses. This drops out of the wet-lab shortlist as soon as it runs.

**Remediation**: scheduled for when Boltz sweep completes (~4.5h from snapshot). Single command:
```bash
python scripts/29_v3_liability_panel.py
```
~1hr DTI run + ~10s gating + `reports/liability_audit_v1.md`.

### 4.3 Cluster C (PrimeKG + TxGNN) still coded but not run

**Symptom**: same as v3 — `scripts/23_v3_cluster_c.py` ready, WSL2 `txgnn_env` venv setup scripted, PrimeKG download (~1.4 GB) scripted but not done.

**Severity**: MEDIUM. v4 added 3 new clusters/layers (Tanimoto, isotonic, pocket); the 4-cluster RRF promised in v3 is now arguably a 5-7-cluster fusion (depending on how you count). Cluster C absence still limits the mechanism-based disagreement archetypes.

**Remediation**: unchanged from v3 — Phase B execution post-Boltz-sweep.

### 4.4 Calibrated DTI grid not yet flowing through fusion

**Symptom**: `scripts/33_v3_apply_calibration.py` writes `data/results/dti_scores_calibrated.parquet`. The fusion script (`scripts/15_v2_fusion.py`) still reads the **raw** `dti_scores.parquet`. The isotonic calibration is shipped but **not yet plugged into the production pipeline**.

**Severity**: MEDIUM — the calibrated parquet exists and the Tier-A/B targets are recoverable, but the fusion ranker doesn't see them.

**Remediation**: add a `--calibrated-mammal` flag to `scripts/15_v2_fusion.py` that reads `dti_scores_calibrated.parquet` and uses `calibrated_pkd` as the MAMMAL ranker score. Cross-target scale heterogeneity (PDE9A clips to 10.44 mean, SLC6A3 to 6.57 mean) means we should Z-normalise within-target *before* RRF — small ~1-day follow-up.

### 4.5 NMDA targets (GRIN2A, GRIN2B) confirmed unfixable by single-chain inference

**Symptom**: Phase A.7 → §7.11 isotonic → §7.5 pocket-conditioned all confirm GRIN2A/2B as Scenario 3 (structural blindness). Ifenprodil-class NAMs bind the **GluN1/GluN2B ATD heterodimer interface** which MAMMAL cannot see from a single-chain input.

**Severity**: STRUCTURAL — calibration can't fix it; pocket-conditioning can't fix it; ensemble with another single-chain DTI head probably can't fix it either.

**Remediation options**:
- **Deprecate** GRIN2A/2B from the cognition panel. Use the remaining 20 targets.
- **Heterodimer cofold** via Boltz-2 with GluN1+GluN2B template (3QEL) — needs custom YAML construction for the cofold; meaningful 2-3 day project.
- Accept low-confidence rank at these targets, flag `INVERTED_TARGET_TOP` provenance.

The methodology note already permits all three.

### 4.6 Tanimoto baseline has scaffold-novelty bias

**Symptom**: Gate G3 in the faceted shortlist failed because TC-5619 and encenicline (novel CHRNA7 PAM scaffolds) under-rate vs canonical drug-like scaffolds. The Tanimoto-to-actives floor is good at separating "this resembles known binders" from "this doesn't" but penalises chemistry that's novel-relative-to-ChEMBL.

**Severity**: MEDIUM — known limitation of the Morgan-FP-Tanimoto family. Affects ~5-10% of expected top compounds.

**Remediation**: **§7.7 cross-DTI ensemble** is the principled fix — a contrastive-trained DTI like PSICHIC has different bias structure. Also: a generative-chemistry source (REINVENT 4 conditioned on the calibrated panel) would surface compounds the Tanimoto baseline can't see.

### 4.7 HRH3 BOLTZ_2X_MAMMAL verdict still rests on n=3

**Symptom**: V3 calibration tagged HRH3 as `BOLTZ_2X_MAMMAL` based on Spearman ρ = +0.87 with n=3 Boltz predictions. After the overnight sweep finishes (~50 more HRH3 pairs), this CI will tighten significantly — verdict may keep, flip to EQUAL_WEIGHTS, or revert to MAMMAL_ONLY_*.

**Severity**: LOW (statistical artifact, not a failure mode) — but the wet-lab shortlist's HRH3 ranking depends on it.

**Remediation**: re-run `scripts/22_v3_calibration.py` post-sweep.

### 4.8 Cross-target scale heterogeneity in calibrated DTI

**Symptom**: per-target isotonic calibrators map to different Y-axis scales. PDE9A's calibrated values cluster around 10.44 mean; SLC6A3 around 6.57. When we compute panel-wide selectivity (Gini) over the 22-vector, PDE9A artificially dominates because its scale is higher.

**Severity**: MEDIUM — broke the §7.4 + §7.5 integration attempt (`selectivity_v2_isotonic.md` shows everything still flat).

**Remediation**: Z-normalise within-target before selectivity computation. Add a `_zscore_within_target` flag to `selectivity/gini_scorecard.py`. Small (~half-day) follow-up; **the Tanimoto-based selectivity continues to be the production deployment** until this lands.

### 4.9 Out-of-range clipping at calibrated DTI extrapolation

**Symptom**: `IsotonicRegression(out_of_bounds='clip', y_min=2.0, y_max=11.0)` maps any compound with raw_pkd outside the calibration training range to the calibrator's Y boundary (often 9.92 at SLC6A3). Top-of-calibrated-list at SLC6A3 contains extrapolation artifacts (levetiracetam / valproate / lithium not classical DAT binders).

**Severity**: MEDIUM — affects top-K precision at the inverted targets specifically.

**Remediation**: add a `calibrator_in_range` boolean column to `dti_scores_calibrated.parquet`. When `raw_pkd` outside [raw_min, raw_max] of the calibration set, flag for fusion down-weighting. Small follow-up.

### 4.10 Methodology note v1 has 2 post-ship update callouts (deferred to v2)

**Symptom**: `reports/methodology_v1.md` got post-ship update banners after each breakthrough. The note is now 3 commits behind the actual architecture.

**Severity**: LOW — the callouts are honest; new contributors get pointed to the right artifacts.

**Remediation**: write `reports/methodology_v2.md` after the next sprint that bundles the v4 architecture into a coherent narrative. ~1-day rewrite.

---

## 5. Throughput Bottlenecks (revised, post-v4)

| Operation | Per-call wall-clock | Total grid | ETA |
|---|---|---|---|
| MAMMAL DTI inference | 100 ms (batch 8) | 6,556 pairs | ~10 min |
| ESM2-650M embedding | 5 s (single target) | 22 targets | ~5 min |
| ADMET-AI 41 endpoints | 100 ms (CPU) | 298 compounds | ~30 s |
| Boltz-2 affinity, WSL2 + cuequivariance | 67 s | 1,165 pairs | ~22 h (current sweep) |
| ChEMBL evidence backstop (SQLite) | 900 ms | 6,556 pairs | ~99 min (done) |
| Per-target Spearman ρ vs ChEMBL | 5 s | 22 targets | ~2 min |
| Phase A.7 calibration end-to-end | ~3 min | one call | ~3 min |
| Phase C fusion (both passes) | ~5 s each | two calls | ~10 s |
| **§7.11 isotonic per-target sweep + LOCO + bootstrap** | **~7 s per target** | **22 targets** | **~3 min** |
| **§7.11 calibration deploy on full grid** | **~50 ms per pair** | **6,556 pairs** | **~5 s** |
| **§7.4 selectivity scoring (Tanimoto vector)** | **~50 ms per compound** | **115 compounds** | **~15 s** |
| **§8.1 faceted shortlist** | **~10 ms per facet** | **17 facets × top-5** | **~2 s** |
| **Tanimoto-to-actives ranker (cluster A.4 fusion input)** | **~3 ms per pair** | **6,556 pairs** | **~20 s** |
| **§7.5 pocket DB build + 13 centroids (one-time)** | **~30 s** | **7 PDBs fetched + parsed** | **~30 s** |
| **§7.5 classify_pose (per pose)** | **~1 ms** | **1,165 pairs** | **~1 s** (when poses available) |
| Liability panel — Stage 1 UniProt enrich | ~700 ms per target | 44 targets | ~30 s (done) |
| Liability panel — Stage 2 MAMMAL DTI (queued) | 100 ms (batch 8) | 298 × 44 = 13,112 pairs | ~22 min wall, ~1 hr with setup |
| Liability panel — Stage 3 gating + render | ~5 s | one call | ~5 s |
| RRF fusion + provenance + 5 clusters | 300 ms | one call | ~1 s |

**Remaining slow operation**: the WSL2 Boltz sweep. 67 s/pair on the RTX 5070. At 79.5% complete, ETA 4.5 h. Next-session priorities (after sweep):
1. Liability panel MAMMAL re-run (1 hr)
2. Phase A.7 → C → D refresh with full Boltz coverage (~10 min)
3. Phase B Cluster C run (~1 hr including PrimeKG download)
4. Pose extraction wrapper + a second Boltz pass (~22 h to re-do, OR ~2 days to add pose-saving and rerun only the 239 pairs not yet done)

**Total compute budget for the v4 → v5 transition**: ~3-4 hr post-sweep (excluding the pose-extraction re-run, which is the longest discretionary item).

---

## 6. Component Status Matrix (one-pass reference)

### 6.1 Library installs (unchanged from v3 + `betacal`, `biopython`)

| Library | Windows native | WSL2 mammal_env | WSL2 txgnn_env | Notes |
|---|---|---|---|---|
| torch 2.12 nightly cu128 | ✅ | ✅ | ❌ | |
| torch 2.7.0+cu128 | — | — | ✅ | txgnn_env requirement |
| biomed-multi-alignment (mammal) | ✅ | ✅ | — | |
| admet-ai | ✅ | ✅ | — | |
| boltz | ✅ slow | ✅ fast | — | WSL2 + cuequivariance kernels |
| cuequivariance-ops-torch-cu12 | ❌ Linux only | ✅ | — | |
| chembl-downloader | ✅ | ✅ | — | |
| pyg_lib / torch_scatter / etc. | ❌ | ❌ | ✅ | |
| txgnn | — | — | ✅ | |
| igraph | ✅ | ✅ | ✅ | PrimeKG loader |
| rdkit | ✅ | ✅ | — | |
| lightgbm | ✅ | ✅ | — | LambdaMART; not promoted |
| **biopython** | **✅** | **✅** | — | **NEW for §7.5 pocket DB (centroid extraction)** |
| **betacal** | **✅** | **✅** | — | **NEW but DEFERRED (binary-only; see §2.3)** |
| **sklearn IsotonicRegression** | **✅** | **✅** | — | **NEW for §7.11 (sklearn dep already present)** |

### 6.2 V4 modules (additions on top of v3)

| Module | LOC (approx) | Status | Notes |
|---|---|---|---|
| `diagnostics/__init__.py` | ~30 | LIVE | |
| `diagnostics/prior_collapse.py` | ~120 | LIVE | norm_y_mean / norm_y_std reference; per-target std vs SD=1.34 ratio |
| `diagnostics/power_analysis.py` | ~170 | LIVE | Bonett-Wright Fisher-z CI + 10K permutation |
| `diagnostics/scaffold_saturation.py` | ~140 | LIVE | Diagnostic A — Murcko cluster overlap |
| `diagnostics/distribution_shift.py` | ~110 | LIVE | Diagnostic B — K-S + Wasserstein |
| `diagnostics/tanimoto_correlation.py` | ~190 | LIVE | Diagnostic D (highest-value test) |
| `diagnostics/tanimoto_baseline.py` | ~140 | LIVE | THE breakthrough comparison vs MAMMAL |
| `diagnostics/temporal_strat.py` | ~140 | LIVE | Lateral 6.2 — pre/post 2015 ChEMBL split |
| `diagnostics/binding_mode_mix.py` | ~170 | LIVE | Lateral 6.1 — ChEMBL action_type |
| `cluster_a/tanimoto_ranker.py` | ~110 | LIVE | Cluster A.4 — Morgan FP + lru_cache + RankerInput |
| `selectivity/__init__.py` | ~30 | LIVE | |
| `selectivity/gini_scorecard.py` | ~260 | LIVE | Graczyk Gini + S(10x) + BCa bootstrap + per-target prior |
| `selectivity/categorize.py` | ~50 | LIVE | mono/dual/poly/flat/intermediate/uncertain |
| `fusion/faceted_shortlist.py` | ~280 | LIVE | 8 mechanism class + 9 targeted-pair + class-match bonus |
| `gates/liability_panel.py` | ~220 | LIVE | 44-target tier-stratified gating + combine_admet_and_liability |
| `calibration/__init__.py` | ~30 | LIVE | |
| `calibration/isotonic.py` | ~200 | LIVE | sklearn IsotonicRegression with LOCO + bootstrap + pickle |
| `calibration/beta_cal.py` | ~140 | DEFERRED | binary-classifier limitation noted |
| `calibration/router.py` | ~170 | LIVE | §1D decision matrix + post_fit_tier (A/B/C/D) |
| `calibration/diagnostics.py` | ~100 | LIVE | LOCO + bootstrap helpers |
| `pockets/__init__.py` | ~30 | LIVE | |
| `pockets/pocket_database.py` | ~210 | LIVE | YAML loader + RCSB fetch + biopython extractor (residue OR ligand-anchor) |
| `pockets/pocket_classifier.py` | ~110 | LIVE | Geometric classifier + dual-site detection + rank-multiplier |

### 6.3 V4 scripts (one entry point per task)

| Script | Purpose | Status |
|---|---|---|
| `27_v3_selectivity_scoring.py` | Gini + S(10x) per compound | LIVE |
| `28_v3_faceted_shortlist.py` | 17-facet shortlist with provenance | LIVE |
| `29_v3_liability_panel.py` | 3-stage: enrich → DTI → gate | Stage 1 done |
| `30_v3_diagnose_inverted.py` | 7-diagnostic protocol orchestrator | LIVE |
| `31_v3_tanimoto_baseline.py` | THE breakthrough comparison | LIVE (showed Tanimoto wins 7/0) |
| `32_v3_calibration_comparison.py` | 22-target isotonic sweep + router | LIVE (SLC6A3 +0.62) |
| `33_v3_apply_calibration.py` | Deploy calibrators on full DTI grid | LIVE |
| `34_v3_build_pocket_database.py` | Curated pocket DB + 13-gate validator | LIVE (13/13 pass) |

### 6.4 V4 reports

| Report | Description |
|---|---|
| `methodology_v1.md` | Original v3 + 2 post-ship update banners for v4 breakthroughs |
| `diagnostics_v1.md` | Prior collapse + power analysis + 5 lateral diagnostics |
| `tanimoto_baseline_v1.md` | The 7/0 victory table |
| `selectivity_v1.md` | Gini on raw MAMMAL — degenerate (all flat) |
| `selectivity_v1_tanimoto.md` | Gini on Tanimoto vector — works, mechanism-diverse top-10 |
| `selectivity_v2_isotonic.md` | Gini on isotonic-calibrated — broken (scale heterogeneity, §4.8) |
| `wet_lab_shortlist_v4_faceted.md` | THE deliverable: 17 facets × top-5 + gates G3-G6 |
| `wet_lab_shortlist_v4_tanimoto.md` | 5-cluster fusion w/ Tanimoto (donepezil #4) |
| `fusion_tanimoto_addition_diff.md` | Tanimoto-on vs Tanimoto-off (Spearman +0.967) |
| `calibration_comparison_v1.md` | §7.11 22-target isotonic sweep + Tier A/B/C/D |
| `calibration_apply_v1.md` | Per-target before/after dynamic range |
| `pocket_database_v1.md` | 13 centroids + 13-gate validation report |
| Plus 9 from V3 (allosteric, ChEMBL audits, methodology v1, etc.) | — |

### 6.5 Architecture at a glance (v4)

```
                            ┌──────────────────────────┐
                            │      INPUT LAYER          │
                            │   22 cognition targets    │
                            │   298 compounds            │
                            └─────────────┬────────────┘
                                          │
        ┌──────────────────────┬──────────┴───────┬──────────────────────┐
        ▼                      ▼                  ▼                      ▼
┌────────────────┐ ┌────────────────────┐ ┌────────────────┐ ┌────────────────────┐
│ CLUSTER A      │ │  CLUSTER B          │ │ CLUSTER C       │ │ §8.0b LIABILITY    │
│ Structure       │ │  ADMET / safety     │ │ Mechanism / KG  │ │ (NEW v4)            │
│ ├ A.1 MAMMAL    │ │  ADMET-AI 41 EP     │ │ PrimeKG + TxGNN │ │ 44-target sub-panel│
│ ├ A.2 ESM2      │ │  ├ BBB / hERG       │ │ vs cognition    │ │ Tier 1 HARD CUT    │
│ ├ A.3 Boltzina  │ │  ├ DILI / P-gp      │ │ virtual anchor  │ │ Tier 2 FLAG        │
│ │  (cuEquiv WSL2)│ │  └ CYPs / Ames     │ │ (5 EFO IDs)     │ │ Tier 3 informational│
│ └ A.4 Tanimoto★ │ │                     │ │                 │ │                     │
│   (1996 FP×     │ │  HARD GATES         │ │  CODED, NOT YET │ │  CODED, GPU-QUEUED  │
│   ChEMBL≥8 acts)│ │  + regulatory bypass│ │  RUN (Phase B)  │ │  (~1hr re-run)      │
└────────┬────────┘ └─────────┬───────────┘ └────────┬────────┘ └─────────┬──────────┘
         │                    │                       │                    │
         │ §7.11 isotonic     │ §8.0b liability       │                    │
         │ per-target ★       │ gate (CUT > FLAG)     │                    │
         │ calibration        │                       │                    │
         ▼                    ▼                       ▼                    ▼
            ┌────────────────────────────────────────────────────┐
            │             FUSION LAYER (5-cluster)                │
            │  RRF (k=60) over: MAMMAL + Tanimoto + Boltzina +    │
            │                   ADMET + (Cluster C when live)     │
            │  Per-target weights via Phase A.7 calibration       │
            │  + Tier A/B/C/D escalation routing                  │
            │  + §7.5 pocket_class rank-multiplier (when poses    │
            │    available)                                       │
            └────────────────────────┬───────────────────────────┘
                                     ▼
            ┌────────────────────────────────────────────────────┐
            │            §8.1 FACETED SHORTLIST (NEW v4)          │
            │  8 mechanism-class facets × top-5 each + 9          │
            │  targeted-pair facets × top-5 each + cross-facet    │
            │  provenance (compound → list of facets it appears in)│
            │  + Selectivity (Gini, S(10×), mono/dual/poly/flat)  │
            └────────────────────────┬───────────────────────────┘
                                     ▼
            ┌────────────────────────────────────────────────────┐
            │            V4 WET-LAB SHORTLIST DELIVERABLE         │
            │  reports/wet_lab_shortlist_v4_faceted.md            │
            │  Gates: P1/P2/P3 (pocket) + G3-G6 (faceted) live   │
            │  Honest caveats: scaffold novelty bias on Tanimoto, │
            │  HRH3 BOLTZ verdict on n=3, GRIN2A/2B deprecation   │
            └────────────────────────────────────────────────────┘

   ★ = v4 additions vs v3 architecture
```

---

## 7. Architecture Enhancements Still Needing Research

These items are unchanged-and-still-pending from v3 §7, **excluding** items shipped in v4. Section numbering preserved for cross-reference continuity.

### 7.7 Cross-DTI ensemble — now the #1 next priority

**Concept**: add MMAtt-DTA, PSICHIC, or BALM as a second/third DTI head in the fusion. Per the V3 research stream A, MMAtt-DTA reports ρ > 0.72 on the transporter superfamily — the only published DTI head whose transporter performance is in the same league as our Tanimoto baseline (+0.90).

**Why now**: v4 raised the bar — Tanimoto already gives ρ +0.90 at SLC6A3. A new DTI head must beat that to justify the integration cost. The most likely path: **PSICHIC for non-transporter targets** (different architecture; physicochemical GNN; likely orthogonal failure modes) + **MMAtt-DTA for transporters** (the published superfamily ρ leader).

**Effort**: 2-3 days for MMAtt-DTA (pip-installable, MIT), 1-2 days each for PSICHIC + BALM.

**Dependencies**: pip-only. RankerInput is templated. Drop-in compatible with `scripts/15_v2_fusion.py`.

**Validation gates**: per-target ρ vs ChEMBL pchembl ≥ +0.40 at SLC6A3 (must beat Tanimoto). EnsDTI-style gating (Park 2024 bioRxiv 2024.08.06.606753) routes between MAMMAL / Tanimoto / MMAtt-DTA per target with confidence weighting.

### 7.8 Generative top-up via REINVENT 4

**Concept**: use the calibrated panel as a fitness function for REINVENT 4 / POLYGON / PILOT. Same as v3 §7.8.

**Why now (revised)**: the v4 panel is *cleanly* calibrated (isotonic per-target where Tier A; Tanimoto fallback elsewhere; pocket class flags where available). Generation conditioned on this richer fitness function should produce candidates that beat the current top-25 *and* satisfy pocket / liability constraints.

**Effort**: 7-14 days unchanged.

### 7.9 Cluster D — neurobiological prior via AHBA + OT Genetics

Unchanged from v3 §7.9.

### 7.10 LLM-agent prioritisation overlay (MCP server)

Unchanged from v3 §7.10. **Updated rationale**: the v4 pipeline has 5 clusters + 3 calibration layers + faceted output + pocket provenance. The complexity now genuinely requires natural-language frontending.

### 7.12 Conformal prediction per-target gating

**Status update**: v4's isotonic calibration covers most of the §7.12 use case (per-target trust gating via post_fit_tier A/B/C/D). Conformal is now a **refinement** for v5, not a Tier-2 priority.

### 7.13 Scaffold-aware active learning

**Status update**: the v4 diagnostic confirmed the SLC6A3 INVERTED finding came from tropane-scaffold dominance in BindingDB training. Scaffold-aware AL is now informed by concrete evidence (we know the failure mode); the priority went up.

### 7.14 Residual-correction XGBoost meta-ranker

**Status update**: v4 isotonic calibration shipped at SLC6A3 with ρ=+0.62. XGBoost residual head would have been the alternative; isotonic was simpler and worked. XGBoost is now optional, not urgent.

### 7.15 NEW — Hierarchical Bayesian (PyMC) for GRIN family pool

**Concept**: research/4-tier/Isotonic-PerTarget-Calibration.md §2 spec'd a PyMC Neelon-Dunson hierarchical isotonic with family-level hyperpriors (SLC6 pool, GRIN pool). v4 shipped classical isotonic only; the hierarchical pool is queued.

**Why now**: GRIN2A (n=7) and GRIN2B (n=14) are below the safe-isotonic threshold (n=15). Hierarchical pooling would borrow strength across the GRIN subunits. Predicted improvement: GRIN2B post-cal ρ from -0.17 to +0.20-0.35.

**Effort**: 3-5 day PyMC project + Rhat / BFMI convergence checks. The empirical-Bayes shortcut (penalised PAVA with empirical mean prior) is the documented fallback if MCMC convergence fails.

**Dependencies**: `pymc` (heavy install; numpyro backend for speed).

### 7.16 NEW — §7.5 detector ensemble (P2Rank + PocketMiner + CryptoBench)

**Concept**: research/4-tier/Pocket-Conditioned-Boltz2.md §1 spec'd a 4-detector consensus. v4 shipped curated-DB + geometric classifier only (the Sprint 1 MVP). Sprint 2 adds the three ML detectors.

**Why now**: the geometric classifier hits 13/13 gates on known PDB-bound ligands, but it has no way to detect *cryptic pockets that open on ligand binding*. PocketMiner + CryptoBench are specifically trained on that signal. The Boltz cofold second-opinion (when MAMMAL vs Boltz affinity diverge > 1.5 log) is the third leg.

**Effort**: 1-2 weeks. P2Rank is a Java jar (cheap); PocketMiner is PyTorch (GVP-GNN); CryptoBench needs ESM-2 embeddings (compatible with our cache if dimension matches — 650M vs 3B is the deciding factor).

### 7.17 NEW — Pose-saving Boltz wrapper (operationalises §7.5)

**Concept**: extend `_wsl2_boltz_full_sweep.py` to save mmCIF poses (or extract centroid in-line and write `pose_centroid_x/y/z` to the parquet). Without this, §7.5 cannot annotate the live grid.

**Why now**: §7.5 is shipped and validated; only blocked on data.

**Effort**: 1 day. Then either re-run the sweep (~22h, throws away ~926 done pairs) or run a pose-only pass on the existing 1,165 pairs (~6-10h since affinity scoring is skipped).

### 7.18 NEW — Cross-target Z-normalisation of calibrated DTI for selectivity

**Concept**: per-target isotonic calibration creates scale heterogeneity that breaks panel-wide Gini. Z-normalise within-target before computing selectivity.

**Why now**: shipped + documented in §4.8 above. Half-day fix.

---

## 8. Creative / Out-of-the-Box Additions (gap-filling brainstorm)

V3 listed 12+ subsections (8.0a-c top-tier + 8.1-12 hand-thought). V4 status:

### Shipped in v4

- **§8.1 multi-class facets** ✅ shipped (`fusion/faceted_shortlist.py`)
- **§8.0b PDSP-style 44-target safety panel** 🟡 infrastructure shipped, GPU-queued

### Still queued from v3 (high-value)

### 8.0a Pareto NSGA-III restructure of fusion

Unchanged. The v4 Tanimoto + facets dissolved the worst lock-in symptom; Pareto across 5 axes (efficacy / safety / novelty / IP / route) remains the principled long-term restructure.

### 8.0c Cluster D = neurobiological prior (AHBA + OT Genetics)

Unchanged from v3.

### 8.2 Combination-screening top-up (DrugComb)

Unchanged. **Updated relevance**: the v4 faceted shortlist surfaces compounds that hit *one* mechanism cleanly. Combinations across mechanism facets (e.g., cholinergic + glutamatergic) are exactly what a DrugComb lookup would surface — natural next step.

### 8.3 Patent / clinical-trial cross-reference

Unchanged. Now feeds Pareto axis 5 (IP freedom).

### 8.4 Literature mining for top-25 (Semantic Scholar + Europe PMC)

Unchanged.

### 8.6 Brain region selectivity via Allen Brain Atlas expression

Unchanged.

### 8.7 Mechanism-of-action embedding (ChEMBL action_type)

**Status update**: v4 diagnostics module added `binding_mode_mix.py` which queries ChEMBL `action_type` per compound. The infrastructure is half-built. Promoting to a fusion-input MoA classifier is now a 1-2 day extension instead of a from-scratch build.

### 8.9 ANI-2x neural potential validation on top-25 poses

**Updated dependency**: requires pose extraction (§7.17). Until poses are saved, ANI-2x rescoring is impossible. Schedule for v5.

### 8.10 Reverse-engineering known nootropics (similarity ranking)

Unchanged.

### 8.11 LLM literature agent for the methodology note

Unchanged.

### 8.12 Cron-driven auto-recalibration

Unchanged. **Updated scope**: when *any* of `boltzina_affinity.parquet`, `dti_scores.parquet`, `weights_calibrated.yaml`, or `dti_scores_calibrated.parquet` changes, trigger Phase A.7 → C → D → 26 → 28 pipeline refresh.

### NEW v4 ideas

### 8.13 Pocket-class-conditioned liability gating (§8.0b × §7.5 integration)

**Concept**: per pocket-conditioned research doc §3.3, the §8.0b liability gates should split by pocket_class:
- 5-HT2B **HARD CUT** only when predicted agonist binds the orthosteric pocket; FLAG for allosteric
- hERG **HARD CUT** only when pose in the central pore cavity (Y652, F656, T623); FLAG for auxiliary-subunit / vestibule binding
- HRH1 **HARD CUT** only when antagonist at orthosteric; FLAG for allosteric
- CB1 **HARD CUT** only at orthosteric (rimonabant precedent); allosteric NAMs have different clinical risk profile

**Why now**: liability panel infrastructure shipped; §7.5 pocket classifier shipped; both ready to compose.

**Effort**: 1-2 days. The §8.0b panel needs pocket centroids for its 44 targets (currently only the 7 priority cognition targets have curated centroids). Two options: (a) reuse a single canonical orthosteric pocket per liability target via UniProt → reference PDB lookup, (b) extend the curated DB to 44 targets (~2 weeks of centroid curation).

### 8.14 Pocket-routed isotonic at SLC6A3 (per Pocket-Conditioned-Boltz2.md §6.4)

**Concept**: the §7.11 isotonic calibrator at SLC6A3/SLC6A2 assumed a single monotone MAMMAL-pKd-vs-ChEMBL-pchembl relationship. If compounds bind S1 (orthosteric) vs S2 (vestibule allosteric) at DAT, the assumption is violated. Solution: route by `pocket_class` BEFORE isotonic; fit separate calibrators for S1 vs S2.

**Why now**: prerequisite for surfacing non-tropane DAT inhibitors (the V3 lateral §7.13 active-learning direction). Without pocket-routing, the calibrator may overfit to the dominant tropane cluster.

**Effort**: 2-3 days. Depends on §7.17 pose extraction.

### 8.15 Tanimoto baseline as an explicit FUSION DIAGNOSTIC

**Concept**: the v4 Tanimoto-to-actives ranker beats MAMMAL at every target. This is itself a fusion-disagreement signal worth surfacing. Add to `provenance/disagreement_report.py`: **per compound, compute |rank(MAMMAL) − rank(Tanimoto)|**. Compounds where the two disagree by >50 ranks are either (a) novel scaffolds Tanimoto under-rates (TC-5619 / encenicline pattern) OR (b) Tanimoto-similar-but-MAMMAL-low compounds that may be activity cliffs.

**Why now**: zero engineering cost; renders as a column in the wet-lab shortlist.

**Effort**: 1 day.

### 8.16 Tier-A calibrator round-trip QC

**Concept**: for Tier-A calibrators (SLC6A3 currently), the fitted isotonic should be **re-evaluated periodically** on new ChEMBL releases. Add a `data/calibration/qc/<target>.json` audit trail with the calibrator's ChEMBL release version + the new ρ when applied to the latest pchembl truth. If post-cal ρ degrades by >0.1, trigger a re-fit.

**Why now**: future-proofing. ChEMBL releases happen ~quarterly.

**Effort**: 1 day.

---

## 9. Roadmap by Priority Tier

**Sequencing principle** (synthesised from v3 → v4 learnings):
1. *Unblock GPU-queued items first* — Boltz sweep finishing is the gate to 3 parallel workstreams
2. *Operationalise shipped infrastructure* — §7.5 + §7.11 are shipped but not yet flowing through fusion; small follow-ups give big leverage
3. *Cross-DTI ensemble (§7.7) is the next research priority* — Tanimoto raised the bar; MMAtt-DTA / PSICHIC / BALM is the empirical test of whether modern DTI heads can beat the 1996 baseline
4. *PyMC hierarchical (§7.15) and detector ensemble (§7.16) are quality refinements* — meaningful but not transformative
5. *Pareto restructure (§8.0a) and generative top-up (§7.8) are paper-grade contributions* — long-horizon

### Tier 1 — DO IMMEDIATELY (gated on overnight Boltz sweep completing, ~4.5h)

- [ ] **Liability panel Stages 2-3** (`scripts/29_v3_liability_panel.py`) — ~1hr MAMMAL DTI + gating
- [ ] **Phase A.7 re-run** with full Boltz coverage — tests INVERTED-target rescue under structure
- [ ] **Phase B Cluster C** run — `scripts/23_v3_cluster_c.py` in `txgnn_env`
- [ ] **Wire calibrated MAMMAL into fusion** (`--calibrated-mammal` flag on `scripts/15_v2_fusion.py`)
- [ ] **Z-normalisation within target for selectivity** — half-day fix unblocks calibrated Gini
- [ ] **Tanimoto-vs-MAMMAL rank-disagreement column** (§8.15) — 1 day, zero cost
- [ ] **Re-render `reports/wet_lab_shortlist_v4_faceted.md`** with all v4 pieces flowing through

### Tier 2 — DO SOON (1-2 weeks each, ordered by signal-to-effort)

- [ ] **Pose-saving Boltz wrapper** (§7.17) — operationalises §7.5 on live grid
- [ ] **Cross-DTI ensemble with MMAtt-DTA** (§7.7) — empirical test vs Tanimoto +0.90
- [ ] **Pocket-class-conditioned liability gating** (§8.13) — composes shipped §8.0b + §7.5
- [ ] **Mechanism-of-action embedding** (§8.7) — half-built via `diagnostics/binding_mode_mix.py`
- [ ] **Methodology note v2** — coherent v4 architecture narrative
- [ ] **Patent / clinical-trial cross-reference** (§8.3) — Pareto axis 5
- [ ] **Combination-screening via DrugComb** (§8.2) — natural follow-up to faceted shortlist
- [ ] **Reverse-engineering known nootropics** (§8.10) — Tanimoto fingerprint similarity to known seeds

### Tier 3 — DO AFTER (2-4 weeks each)

- [ ] **Pareto NSGA-III restructure** (§8.0a) — 5-axis Pareto front; bonus depends on §8.3
- [ ] **PyMC hierarchical for GRIN pool** (§7.15) — GRIN2B from -0.17 to ~+0.20-0.35
- [ ] **§7.5 detector ensemble** Sprint 2 (P2Rank + PocketMiner + CryptoBench) (§7.16)
- [ ] **Pocket-routed isotonic at SLC6A3** (§8.14) — S1 vs S2 vestibule routing
- [ ] **Cluster D neurobiological prior** (§7.9 / §8.0c) — Mansuri 2024 + AHBA + OT Genetics
- [ ] **GWAS-anchored panel expansion** (§7.3) — to 40-80 targets
- [ ] **Brain region selectivity** (§8.6) — Allen Brain Atlas expression
- [ ] **Selectivity entropy + Partition Index + KISS-CL** (§7.4 v2 metrics)
- [ ] **Tier-A calibrator round-trip QC** (§8.16) — quarterly recalibration audit
- [ ] **Cron-driven auto-recalibration** (§8.12)
- [ ] **ANI-2x pose validation on top-25** (§8.9) — depends on §7.17
- [ ] **LambdaMART promotion** — eligible since v3 (275 CORROBORATED ≥ 20-label threshold)
- [ ] **Conformal prediction per-target gating** (§7.12) — refines isotonic Tier classifier
- [ ] **Scaffold-aware active learning** (§7.13) — now informed by tropane-saturation evidence

### Tier 4 — RESEARCH (3-10 days, parallelism / cloud burst)

- [ ] **Generative top-up with REINVENT 4 / POLYGON / PILOT** (§7.8)
- [ ] **Boltzina-Vina-only mode** (V3 attack plan T5) — ~11× speedup; never implemented
- [ ] **LLM-agent prioritisation overlay** (MCP server) (§7.10)
- [ ] **LLM literature agent for methodology** (§8.11)

### Tier 5 — PUBLISHABLE (multi-week, methodology contributions)

- [ ] **THE v4 methodology paper** — "Pocket-aware, ground-truth-calibrated, faceted drug repurposing on a foundation-model DTI head: four diagnostic + corrective interventions that beat the foundation model with a 1996-era baseline." J Cheminform.
- [ ] **Calibration linchpin paper** — Phase A.7 per-target ρ approach + refined verdict matrix. Distinct from main paper.
- [ ] **MAMMAL transporter-failure-mode paper** — if §7.7 confirms tropane-saturation, this is a publishable negative result.
- [ ] **Allosteric awareness benchmark v2** — expand the Phase 4.1 work with Boltz-rescued PAMs.
- [ ] **Pocket-conditioned virtual-screening protocol paper** — generalises §7.5 beyond cognition.
- [ ] **Cognition virtual phenotype anchor validation paper** (deferred from v3).

---

## 10. Risk Register (post-v4)

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| Overnight Boltz sweep doesn't rescue INVERTED targets after isotonic | Low | Medium | Already rescued at SLC6A3 (+0.62) and SLC6A2 (+0.40) via isotonic; structure is incremental signal |
| MMAtt-DTA fails to beat Tanimoto +0.90 on transporters | Medium | Low | Tanimoto stays as Cluster A.4; MMAtt-DTA becomes a third voter with class-conditional weighting |
| Liability panel MAMMAL re-run reveals additional CUTs we didn't predict | Low | Medium | Documented as the *intended* outcome — discrimination that ADMET alone can't deliver |
| Pose extraction wrapper requires re-running ~926 done Boltz pairs | High | Medium | Pose-only pass (~6-10h) is the cheap option; saves the ~22h full re-run |
| GRIN2A/2B deprecation removes a clinically relevant panel target | Medium | Medium | Already documented as Scenario 3 structural-blindness in the methodology note; heterodimer cofold is the v5+ path |
| Tanimoto scaffold-novelty bias misses TC-5619 / encenicline-class compounds | High (already observed) | Low-Medium | §7.7 cross-DTI ensemble is the principled fix; document caveat in shortlist |
| Cross-target scale heterogeneity in calibrated DTI breaks panel Gini | High (observed) | Medium | §7.18 Z-normalisation half-day fix |
| Out-of-range clipping at isotonic extrapolation creates phantom top-K | High (observed) | Medium | §4.9 `calibrator_in_range` flag — fusion down-weights extrapolated values |
| Pocket centroids drift as new cryo-EM structures land | Low | Low | Annual re-curation of pocket DB; flagged in research doc §8 |
| Hierarchical PyMC fails to converge at n=8 (GRIN2A) | Medium | Low | Documented fallback: empirical-Bayes shortcut (penalised PAVA) |
| HRH3 BOLTZ_2X_MAMMAL verdict flips after sweep finishes | Medium | Low | Re-run Phase A.7; verdict is annotated as low-n in `weights_calibrated.yaml` |
| Methodology note diverges from architecture as v5 ships | High (already observed) | Low | Tier-2 methodology v2 rewrite in the roadmap |
| Roberts CA 2020 effect-size ceiling makes any wet-lab spend hard to justify | High (structural) | High | Reframe as methodology + provenance + enrichment; the publication contribution is the calibration / pocket / faceted framework, not the candidate ranking |

---

## 11. Single-Page Cheat Sheet — Commands That Work Today

```powershell
# Windows env (mammal_env, conda)
$envPython = "$env:USERPROFILE\.conda\envs\mammal_env\python.exe"
$repo = "C:\Users\Pierce Lonergan\Documents\GitHub\MAMMAL_Cognitive_Enhancement_Drug_Repurposing"
Push-Location $repo

# === V3 Phase A (SQLite + linchpin + backstop) — unchanged ===
& $envPython scripts\_v3_sqlite_vs_rest_smoke.py
& $envPython scripts\24_v3_audit_chembl_targets_sqlite.py
& $envPython scripts\22_v3_calibration.py
& $envPython scripts\21_v3_chembl_evidence_sqlite.py --all-pairs

# === V3 fusion + V3 wet-lab shortlist ===
& $envPython scripts\15_v2_fusion.py --calibrated-weights NONE --out-suffix _uncalibrated
& $envPython scripts\15_v2_fusion.py --out-suffix _calibrated
& $envPython scripts\25_v3_fusion_diff.py
& $envPython scripts\26_v3_wet_lab_shortlist.py

# === V4 NEW: diagnostics protocol (CPU-only) ===
& $envPython scripts\30_v3_diagnose_inverted.py

# === V4 NEW: Tanimoto baseline test ===
& $envPython scripts\31_v3_tanimoto_baseline.py

# === V4 NEW: 5-cluster fusion with Tanimoto ranker ===
& $envPython scripts\15_v2_fusion.py --add-tanimoto-ranker --out-suffix _calibrated_v4_tanimoto

# === V4 NEW: §7.11 isotonic calibration sweep + deploy ===
& $envPython scripts\32_v3_calibration_comparison.py
& $envPython scripts\33_v3_apply_calibration.py

# === V4 NEW: §7.4 selectivity (use Tanimoto vector — DON'T use raw MAMMAL) ===
& $envPython scripts\27_v3_selectivity_scoring.py --use-tanimoto `
    --out data/results/v2/selectivity_scores_tanimoto.parquet `
    --report reports/selectivity_v1_tanimoto.md

# === V4 NEW: §8.1 multi-class faceted shortlist ===
& $envPython scripts\28_v3_faceted_shortlist.py `
    --selectivity data/results/v2/selectivity_scores_tanimoto.parquet

# === V4 NEW: §7.5 pocket database build + 13-gate validation ===
& $envPython scripts\34_v3_build_pocket_database.py

# === V4 NEW: liability panel (Stage 1 enrich) ===
& $envPython scripts\29_v3_liability_panel.py --enrich-only
# Post-Boltz: full liability re-run (~1hr)
& $envPython scripts\29_v3_liability_panel.py
```

```bash
# === WSL2 ===

# Overnight Boltz sweep status:
wsl -d Ubuntu -- bash -c "tail -c 800 /tmp/wsl2_boltz_sweep.log | tr '\r' '\n' | tail -5"

# Cluster C orchestrator (when GPU free):
wsl -d Ubuntu -u root -- bash -c \
  "source /root/txgnn_env/bin/activate && python /mnt/c/Users/Pierce\ Lonergan/Documents/GitHub/MAMMAL_Cognitive_Enhancement_Drug_Repurposing/scripts/23_v3_cluster_c.py"
```

---

## 12. The Mission Statement (refreshed for v4)

This pipeline is an honest attempt to use a foundation model (MAMMAL) plus a curated panel of cognition-relevant targets plus rigorous filtering (ADMET + structural + mechanistic) plus *ground-truth-calibrated weighting* to surface drug-repurposing candidates with mechanistic plausibility and physical feasibility for healthy adult cognitive enhancement.

The Roberts 2020 meta-analysis ceiling is real: even methylphenidate hits SMD = 0.21 on a generous metric. We are not searching for a miracle drug; we are enriching a candidate set so wet-lab cycles spend money on plausibility, not chemistry-lottery tickets.

**The v3 → v4 transition demonstrated a methodology contribution distinct from the candidate ranking:**
- The diagnostic protocol that exposes prior collapse and statistical-power blockades before they poison downstream verdicts (commit `530dc40`).
- The Tanimoto-to-actives baseline that *empirically* beats the 458M-param foundation model panel-wide — a finding that defines what "fairly evaluating MAMMAL" requires (commit `530dc40`).
- The isotonic per-target calibration that surgically repairs the rare cases where MAMMAL has signal but in the wrong direction (commit `8624fd1`).
- The pocket-conditioned classifier that distinguishes PDE4D allosteric (BPN14770) from catalytic (rolipram) at the geometric level — enabling mechanism-aware liability gating (commit `3884ba5`).
- The faceted shortlist that dissolves single-target lock-in into structured, mechanism-orthogonal top-N tables with cross-facet provenance (commit `1c288a8`).

Together these four shipped, validated systems make the v4 pipeline a *honest* prioritisation engine — not the most accurate possible, not the most novel possible, but the one where every claim about a compound's rank traces to a documented signal source with a known failure mode.

V5 will compose what v4 shipped (pose extraction → pocket-conditioned fusion; liability panel × pocket gating; cross-DTI ensemble × Tanimoto floor; methodology v2 rewrite) and the candidate ranking will be downstream of that architecture.

Life would be a lot better if people had better cognition. The contribution here isn't the next nootropic — it's the framework that lets the next nootropic search be honest about what it doesn't know.

---

## Appendix A — Research citation index (v4 additions)

V3's Appendix A had ~70 citations grouped by topic. V4 adds 5 new research deep-dives + their internal citations:

### A.3 Diagnostic protocol (research/4-tier/Diagnosing-MAMMAL-DTI-Anti-Correlation.md)

- Shoshan et al. 2026 *npj Drug Discovery* / arXiv:2410.22367 — MAMMAL paper (norm_y_mean/std reference)
- Sundar & Colwell 2020 *J Chem Inf Model* 60:56 (doi:10.1021/acs.jcim.9b00415) — AVE debias destroys generalisability
- Zhang et al. 2025 ICLR Oral arXiv:2504.09481 — DTA models degrade severely on low-Tanimoto samples
- Graber et al. 2025 *Nat Mach Intell* article s42256-025-01124-5 — PDBbind CleanSplit; benchmark numbers driven by leakage
- Stanley et al. NeurIPS 2021 — FS-Mol; QSAR median 94 compounds/task
- Dablander et al. 2023 *J Cheminform* — activity-cliff prediction failures
- Lin et al. 2023 (ESM-2/ESMFold) *Science*
- Srivastava et al. 2024 *Nature* 632:672 (doi:10.1038/s41586-024-07739-9) — hDAT cryo-EM
- Yuan et al. 2024 *Cell Research* (doi:10.1038/s41422-024-01024-0) PDB 8ZP1 — hNET
- Karakas/Simorowski/Furukawa 2011 *Nature* 475:249 (doi:10.1038/nature10180) — GluN1/GluN2B ATD ifenprodil interface
- Bonett & Wright 2000 *Psychometrika* 65(1):23 — Fisher-z CI for Spearman ρ
- Carroll et al. 2004 *J Med Chem* PMID 15566309 — RTI tropane series

### A.4 Selectivity scoring (research/4-tier/Graczyk-Selectivity-Faceted-Shortlist.md)

- Graczyk PP 2007 *J Med Chem* 50:5773 (doi:10.1021/jm070562u) — Gini for kinase selectivity
- Bosc/Meyer/Bonnet 2017 *BMC Bioinformatics* 18:17 (doi:10.1186/s12859-016-1413-y) — Gini benchmarking
- Uitdehaag & Zaman 2011 *BMC Bioinformatics* 12:94 (PMC3100252) — entropy + S(10x)
- Cheng et al. 2010 *J Med Chem* 53:4502 (doi:10.1021/jm100301x) — Partition Index
- Karaman et al. 2008 *Nat Biotechnol* — kinase profiling
- *Nat Commun* 2025 article s41467-025-65869-8 — MMCLKin contrastive selectivity
- Olivia Guest reference numpy Gini — github.com/oliviaguest/gini

### A.5 Liability panel (research/4-tier/Cognition-44Target-Liability-Panel.md)

- Bowes et al. 2012 *Nat Rev Drug Discov* 11:909 — original 44-target panel
- Brennan et al. 2024 *Nat Rev Drug Discov* 23:525 — IQ Consortium Safety-77
- Dumotier & Urban 2024 *J Pharmacol Toxicol Methods* 128:107542 (PMID 39032441) — >10× Cmax/Ki framework
- Connolly et al. 1997 *NEJM* 337:581 — fen-phen valvulopathy
- Rothman et al. 2000 *Circulation* 102:2836 — norfenfluramine 5-HT2B agonism
- Schade et al. 2007 *NEJM* 356:29 — pergolide
- Zanettini et al. 2007 *NEJM* 356:39 — cabergoline
- Roth 2007 *NEJM* 356:6 — class warning
- Gray et al. 2015 *JAMA Intern Med* 175:401 — anticholinergic dementia HR 1.54
- Coupland et al. 2019 *JAMA Intern Med* 179:1084 (doi:10.1001/jamainternmed.2019.0677) — AOR 1.49
- Topol et al. 2010 *Lancet* 376:517 (PMID 20709233) — CRESCENDO rimonabant
- Sharretts et al. 2020 *NEJM* 383:1000 (PMID 32905685) — lorcaserin withdrawal
- Blackwell 1967 *Br J Psychiatry* 113:349 — cheese effect
- Simmler et al. 2013 *Br J Pharmacol* 168:458 (PMID 22897747) — modafinil TAAR1

### A.6 Isotonic calibration (research/4-tier/Isotonic-PerTarget-Calibration.md)

- Mervin/Afzal/Engkvist/Bender 2020 *J Chem Inf Model* 60:4546 (PMID 32865408, doi:10.1021/acs.jcim.0c00476) — 40M-pair Venn-ABERS benchmark
- Kull/Silva Filho/Flach 2017 AISTATS PMLR 54:623 + *EJS* 11:5052 — beta-calibration
- Vovk & Petej 2014 UAI pp.829 — Venn-ABERS validity
- Nouretdinov et al. 2018 PMLR 91:1 — Inductive Venn-ABERS Regression (IVAR)
- Neelon & Dunson 2004 *Biometrics* 60:398 (doi:10.1111/j.0006-341X.2004.00184.x) — hierarchical isotonic prior
- Lin & Dunson 2014 *Biometrika* 101:303 (doi:10.1093/biomet/ast063) — GP-projection monotone
- Toplak et al. 2020 (clarified citation): Mervin 2020 is the correct reference
- Schulman et al. 2024 *Bioinformatics* 40(8):btae496 — MMAtt-DTA
- Koh et al. 2024 *Nat Mach Intell* 6:673 — PSICHIC
- Gorantla et al. 2025 *J Chem Inf Model* 65(22):12279 (doi:10.1021/acs.jcim.5c02063) — BALM
- Park et al. 2024 bioRxiv 2024.08.06.606753 — EnsDTI gating
- Zhang et al. 2025 *Nat Mach Intell* s42256-025-01151-2 — eMOSAIC Mahalanobis OOD
- Bailey et al. 2024 *eLife* 12:RP89679 — deep batch AL for drug discovery
- Holzmüller `bmdal_reg` — github.com/dholzmueller/bmdal_reg
- Landrum & Riniker 2024 *J Chem Inf Model* 64:1560 (doi:10.1021/acs.jcim.4c00049) — ChEMBL IC50 noise floor

### A.7 Pocket-conditioned Boltz (research/4-tier/Pocket-Conditioned-Boltz2.md)

- Krivák & Hoksza 2018 *J Cheminform* 10:39 (DOI 10.1186/s13321-018-0285-8) — P2Rank
- Meller et al. 2023 *Nat Commun* 14:1177 (DOI 10.1038/s41467-023-36699-3) — PocketMiner
- Škrhák et al. 2024 *Bioinformatics* 41(1):btae745 (DOI 10.1093/bioinformatics/btae745) — CryptoBench
- Passaro et al. 2025 bioRxiv 2025.06.14.659707 — Boltz-2 + affinity OOD warning
- Burgin et al. 2010 *Nat Biotechnol* 28:63 — PDE4D UCR2 closure mechanism
- Wilson et al. 2018 *Nat Commun* 9:3334 — BPN14770 allosteric
- Noviello et al. 2021 *Cell* 184:2121 — α7 nAChR + PNU-120596 (PDB 7KOX)
- Ludwig et al. 2010 *J Recept Signal Transduct* 30:469 — galantamine β-strand 10 site
- Kowal et al. 2018 *Br J Pharmacol* — galantamine PAM dispute
- Hamouda/Kimm/Cohen 2013 *J Neurosci* 33:485 — physostigmine + galanthamine subunit-interface sites
- Cheung et al. 2012 — donepezil dual-site CAS+PAS (PDB 4EY7)
- Schmidt et al. 2016 *Nature* 532:527 — SIGMAR1 (5HK1, 5HK2)
- Zhou et al. 2022 *Nat Commun* 13:1267 — Xenopus σ1R open/closed conformations
- Wang et al. 2025 *Acta Pharmacol Sin* — H3R-histamine-Gi cryo-EM (PDB 8YUU)
- Ligneau et al. 2007 *JPET* — pitolisant Ki = 0.16 nM
- Xiao P et al. 2021 *Cell* 184:943 — DRD1 + LY3154207 + Gs (7CKZ)
- Zhuang et al. 2021 *Cell Res* 31:593 — DRD1 + dopamine + LY3154207 (7LJD/7LJC)
- Teng et al. 2022 *Nat Commun* 13:3186 — DRD1 PAM MD mechanism
- Nielsen et al. 2024 *Nature* — hDAT-cocaine cryo-EM
- Cheng et al. 2020 *J Chem Inf Model* (DOI 10.1021/acs.jcim.0c00346) — DAT S2 vestibule allosteric
- Pólya 2024 *Bioinformatics* 41:btae745 — CryptoBench

---

## Appendix B — How this doc was written

V3's Appendix B documented that two parallel research-agent streams informed §7 + §8. V4 extends that.

Between V3 (commit `c1eed4d`, 2026-05-25) and this V4 document (commit `f188b82`+), **five additional research deep-dives** were commissioned (one per architectural enhancement). Each landed as a markdown in `research/4-tier/` and was implemented as a code package + script + report bundle. The pattern:

```
identify gap   →   commission research doc   →   build minimal viable implementation
   ↓                       ↓                              ↓
v3 finding            §7.X / §8.X spec              code package + script + report
   ↓                       ↓                              ↓
documented        cited inline in V4 plan          validated against pre-committed gates
```

Three of the five doc → ship cycles **delivered results that matched or beat their pre-committed predictions** (Tanimoto, isotonic, pocket-conditioned). One (selectivity) revealed a methodologically interesting interaction with the prior-collapse finding that required a tactical pivot (use Tanimoto vector, not raw MAMMAL). One (liability) is queued on GPU.

When v5 work begins, the next assistant should re-spawn the research-agent pattern for any v4 → v5 architectural enhancement. The shape of useful prompts is documented in this v4 file and in the 5 research markdowns. The methodology v2 note (Tier 2 roadmap item) will absorb this provenance into a publishable artifact.

---
