# V3 — Status, Calibration Findings, Forward Plan

**Read this first.** This is the source-of-truth status document for the project after the V3 Post-Phase-0.5 Aggressive Sprint shipped Phases A, C, D, E. The architecture spec still lives in [V2_HYBRID_ARCHITECTURE.md](V2_HYBRID_ARCHITECTURE.md); the original v3 plan in [V3_ATTACK_PLAN.md](V3_ATTACK_PLAN.md); this doc supersedes [V2_STATUS_AND_FORWARD_PLAN.md](V2_STATUS_AND_FORWARD_PLAN.md) as the live state-of-the-world.

**Snapshot**: `main @ 530dc40` — diagnostics + Tanimoto ranker shipped.

> **🔥 POST-SPRINT BREAKTHROUGH (commit `530dc40`)** — *Diagnostic Protocol +
> Tanimoto Baseline*: a 1996-vintage Tanimoto-on-Morgan-FP baseline beats
> MAMMAL at **every audited cognition target** (SLC6A3 ρ +0.90 vs -0.70; DRD1
> ρ +0.85 vs +0.29; ACHE ρ +0.81 vs +0.24). MAMMAL prior collapse is
> panel-wide (19/22 SEVERE at >10× collapse vs training SD). The Phase A.7
> ρ values are statistics on noise. **Cluster A.4 (Tanimoto-to-actives) has
> been added as a real ranker in fusion** — donepezil leapt from rank #42 to
> #4 in the 5-cluster output. This RE-PRIORITISES §7.6 (LoRA fine-tune is
> NOT the right next move) and PROMOTES §7.7 (cross-DTI ensemble with
> MMAtt-DTA) to Tier 1. See `reports/diagnostics_v1.md`,
> `reports/tanimoto_baseline_v1.md`, `reports/fusion_tanimoto_addition_diff.md`.

> **🔥 POST-SPRINT BREAKTHROUGH #2 (commit `b3fddfa`)** — *Selectivity layer
> + faceted shortlist + 44-target liability panel*. The HRH3-23/25 lock-in
> from the v3 calibrated output is **DISSOLVED**. New `wet_lab_shortlist_v4_faceted.md`
> surfaces top-5 per mechanism class across 8 classes + 9 targeted pairs,
> with cross-facet provenance. Pitolisant correctly ranks #2 in histaminergic
> facet; donepezil in cholinergic top-5; atomoxetine/fluoxetine/duloxetine
> in noradrenergic; (s)-AMPA and ampakines dominate glutamatergic AMPA.
> The 44-target off-target liability panel (§8.0b) is built and ready — Stage 1
> (UniProt enrichment) shipped; Stages 2+3 (MAMMAL DTI + gating) queued for
> when WSL2 Boltz frees the GPU. Expected effect: 2 hard CUTs (aripiprazole,
> amitriptyline on 5-HT2B+HRH1) + ~7 FLAGs that ADMET-AI alone misses.
> See `reports/wet_lab_shortlist_v4_faceted.md`, `reports/selectivity_v1_tanimoto.md`.

> **🔥 POST-SPRINT BREAKTHROUGH #3 (commit `d0b4bd7`)** — *§7.11 Isotonic
> per-target post-hoc calibration*. `IsotonicRegression(increasing='auto')`
> naturally absorbs sign inversion at MAMMAL_ONLY_INVERTED targets,
> replacing the awkward `weight=0.30` hack. **Headline**: SLC6A3 raw ρ=-0.70 →
> post-cal ρ=**+0.62** (Δρ=+1.32, CI=[+0.71, +0.80], Tier A); SLC6A2 -0.60 →
> +0.40 (Tier B). The research doc predicted [+0.45, +0.65] for SLC6A3; we
> landed at the high end. GRIN2A/2B confirmed as Scenario 3 (structural
> blindness at ifenprodil dimer-interface) — calibration can't fix it; they
> escalate to §7.7 cross-DTI ensemble. Decision router (`calibration/router.py`)
> implements §1D matrix + Tier A/B/C/D post-fit classifier. Calibrated DTI
> grid at `data/results/dti_scores_calibrated.parquet`. Beta-calibration
> deferred (betacal package is binary-classifier-only); PyMC hierarchical
> deferred to v2. See `reports/calibration_comparison_v1.md`,
> `reports/calibration_apply_v1.md`, `data/calibration/router_decisions.csv`.

> **🔥 POST-SPRINT BREAKTHROUGH #4 (this commit)** — *§7.5 Pocket-conditioned
> Boltz-2 MVP: curated centroid database + geometric classifier, all 13
> validation gates pass*. 7 priority targets (CHRNA7 / GRIN2B / PDE4D /
> SIGMAR1 / HRH3 / DRD1 / ACHE), 13 pockets curated from RCSB PDBs via
> biopython, with both residue-derived AND ligand-anchored centroid modes.
> The headline win is **PDE4D BPN14770 (allosteric, UCR2) vs rolipram-class
> (orthosteric, catalytic) discrimination at the geometric level** — enables
> §8.0b emesis-liability gating split. ACHE donepezil correctly classified
> as orthosteric via dual-site (CAS+PAS span) detection. P1=4/4 orthosteric
> + P2=2/2 allosteric + P3=7/7 negative-control = 13/13 gates pass.
> Sprint 2 (P2Rank + PocketMiner + CryptoBench consensus) and Sprint 3
> (§8.0b liability split + §8.1 new facets) deferred. Pose-conditioning
> goes operational once a pose-saving Boltz wrapper is added. See
> `reports/pocket_database_v1.md`, `data/pockets/pocket_database.yml`,
> `src/mammal_repurposing/pockets/`.

---

## 1. Executive Summary

The 4-cluster hybrid is shipped end-to-end. Calibration vs ChEMBL ground truth is the linchpin we built; it exposed a finding the rest of v3 must respond to.

**The headline finding** (Phase A.7 calibration report):

> Of 22 cognition panel targets, only **2** have MAMMAL Spearman ρ ≥ +0.30 against ChEMBL pchembl_value (DRD1 ρ=+0.31, HCRTR1 ρ=+0.37). **14** are weakly informative (0 ≤ ρ < 0.30). **4 are anti-correlated** at ρ ≤ -0.30 — SLC6A3 (DAT, ρ=-0.71, n=26), SLC6A2 (NET, ρ=-0.53, n=25), GRIN2A (ρ=-0.35), GRIN2B (ρ=-0.30). KCNQ3 has no joinable cluster data.

This is not a pipeline failure — this is the *first time we measured* per-target signal quality against ground truth instead of trusting the foundation model's BindingDB-derived prior. The pipeline ships correctly; the calibrated weights now down-weight 18/22 targets and de-weight 4 to 0.30. The question for v3 is: **what do we do about the structural bias the calibration exposed?**

| Component | Status | Evidence at this snapshot |
|---|---|---|
| Cluster A.1 — MAMMAL DTI | ✅ LIVE | 6,556 pairs; per-target ρ in `reports/calibration_report.md` |
| Cluster A.2 — ESM2-650M cache | ✅ LIVE | 22 targets cached |
| Cluster A.3 — Boltz-2 / Boltzina (WSL2 + cuequivariance) | ✅ LIVE | 92 pairs done; **overnight sweep 22% of 1,165** |
| Cluster B — ADMET-AI 41 EP + hard gates | ✅ LIVE | 53 PASS / 64 FLAG / 181 CUT |
| Cluster C — PrimeKG + TxGNN | 🟡 CODE LIVE, NOT RUN | WSL2 venv ready; blocked on overnight sweep freeing GPU |
| ChEMBL 36 SQLite mirror | ✅ LIVE | `~/.data/chembl/36/chembl_36.db`; A.5 PASS 19/20 |
| Phase A.7 Calibration linchpin | ✅ SHIPPED | `reports/calibration_report.md`, `configs/weights_calibrated.yaml` |
| Phase C 4-cluster RRF (calibrated + uncalibrated) | ✅ LIVE | both passes produce versioned parquets |
| Phase D calibrated-vs-uncalibrated diff | ✅ SHIPPED | `reports/fusion_calibration_diff.md` (ρ = +0.994) |
| Phase E methodology note v1 | ✅ SHIPPED | `reports/methodology_v1.md` |
| ChEMBL evidence backstop (Phase A.4) | ✅ COMPLETE | 6,556 rows; 275 CORROBORATED / 18 AMBIGUOUS / 6 CONTRADICTED / 6,257 NOVEL |
| Wet-lab shortlist v3 (4-cluster scorecards) | ✅ SHIPPED | `reports/wet_lab_shortlist_v3.md` |

**The single most important pending result**: when the overnight Boltz sweep completes (~17h from this snapshot), re-running Phase A.7 → C → D will tell us whether Boltzina's structure-aware affinity *rescues* the 4 INVERTED-MAMMAL targets (flipping them to `BOLTZ_2X_MAMMAL`) or *confirms* they are structurally bad targets (flipping them to `DE_WEIGHT_TARGET`). That distinction reshapes the next sprint.

---

## 2. What the Calibration Linchpin Actually Exposed

The Phase A.7 report (`reports/calibration_report.md`) is the new ground truth about pipeline trustworthiness. This section unpacks it because the rest of v3 hangs off it.

### 2.1 The four MAMMAL_ONLY_INVERTED targets

| Target | UniProt | ρ vs ChEMBL | n joined preds | Truth records | Mechanism notes |
|---|---|---|---|---|---|
| SLC6A3 | Q01959 | **-0.71** | 26 | 3,391 | DAT, cocaine/methylphenidate target |
| SLC6A2 | P23975 | **-0.53** | 25 | 3,635 | NET, atomoxetine target |
| GRIN2B | Q13224 | -0.30 | 14 | 3,228 | NMDA NR2B subunit |
| GRIN2A | Q12879 | -0.35 | 8 | 239 | NMDA NR2A subunit |

These are not rare/edge-case targets — they are the **most clinically validated cognition targets in the panel** (modafinil, methylphenidate, d-amphetamine all hit SLC6A3 + SLC6A2; memantine works at GRIN2B). MAMMAL ranking their top candidates as worst binders is the single most important known failure mode of the v3 pipeline as it ships.

**Three hypotheses** (the methodology note tags this as not-yet-explained):

1. **BindingDB sampling bias toward narrow chemical space**: BindingDB has dense coverage of high-affinity transporter inhibitors clustered in the tropane / phenethylamine scaffolds. MAMMAL learns "structurally similar to known cocaine analog" → "high pKd," but pKd magnitude within that cluster is invariant to subtle changes that ChEMBL pchembl distinguishes. The signal exists; the regression target is collapsed.

2. **Tokenisation artefacts near the substrate-binding pocket**: the SLC6 family has a tightly conserved substrate-binding region. If MAMMAL's protein tokenisation truncates or mis-aligns that region, predictions become rank-uninformative across compounds.

3. **ChEMBL vs BindingDB chemical-space mismatch**: BindingDB and ChEMBL overlap in compound coverage but not perfectly. The compounds with high ChEMBL pchembl might be systematically *under-sampled* in BindingDB at SLC6A3/SLC6A2.

We don't know which is right. Diagnosing this is v3 priority.

### 2.2 The two MAMMAL_ONLY_STRONG targets

DRD1 (ρ=+0.31, n=21) and HCRTR1 (ρ=+0.37, n=6). Both are GPCRs with well-characterised orthosteric pockets and dense BindingDB coverage of selective antagonists. These are the targets where MAMMAL is doing what it was trained to do.

**Implication for the panel**: DRD1 + HCRTR1 are the only two targets the v3 wet-lab shortlist can lean on at default MAMMAL weight. Compounds whose `mammal_best_target` is one of these get the strongest signal-to-noise; compounds whose `mammal_best_target` is one of the inverted four are essentially being scored by a coin flip with a malicious thumb on it.

### 2.3 The one BOLTZ_2X_MAMMAL target (HRH3) — and why it's fragile

HRH3 (ρ_B = +0.87 with **n=3** vs ρ_M = -0.14 with n=12). On the strength of n=3, calibration flipped HRH3 to "trust Boltz 2×." This is the only target where the v3 calibration's verdict matrix produced a non-default override toward structure.

The Phase D diff shows the consequence: 23 of 25 top compounds in the calibrated shortlist have `mammal_best_target = HRH3`, because it's the only un-down-weighted target — everything else got pushed down by the WEAK / INVERTED override. So the v3 shortlist is in practice mostly an "H3 antagonist + ADMET-clean" ranking.

This is fragile. n=3 is statistically thin; the overnight sweep adding ~50 more HRH3 Boltz predictions might keep, flip, or destroy this verdict. Re-run priority is high.

### 2.4 What the ChEMBL backstop (Phase A.4) confirmed independently

Phase A.4 wrote 6,556 (target, compound) evidence rows. Status breakdown:
- **NOVEL**: 6,257 (95.4%) — most pairs have no ChEMBL bioactivity record
- **CORROBORATED** (pchembl ≥ 6.0, ~1 µM): **275** (4.2%)
- **AMBIGUOUS** (5 ≤ pchembl < 6): 18 (0.3%)
- **CONTRADICTED** (pchembl < 5): 6 (0.1%)

Notable: 4 PDE9A binders at pchembl 11.0 (~10 pM), fenpropimorph @ SIGMAR1 pchembl 11.0, piclamilast @ PDE4D pchembl 10.7. **These are real, externally-validated tight binders our pipeline surfaces** — even at targets where per-target Spearman ρ is weak. The two statements ("MAMMAL ρ at this target is low" and "MAMMAL's top picks at this target are real binders") are both true. ρ is a *correlation* statistic; the pipeline can be rank-informative at the top end while uncorrelated in aggregate. The Phase A.7 linchpin caught this nuance; the wet-lab shortlist v3 now surfaces the CORROBORATED count per compound so a human can see it.

### 2.5 The 14 MAMMAL_ONLY_WEAK targets — the silent majority

ACHE (the donepezil target) at ρ=+0.20 with n=10. CHRNA7 at ρ=-0.01 with n=10. GRIA1-4 family at ρ=+0.10 to -0.24. The cluster of targets the v1 sanity gate "passed" because positive controls ranked well within our 298-compound library are exactly the targets the calibration linchpin marks as poorly informative across all of ChEMBL.

**Both statements remain true**. The sanity gate measured rank percentile within a curated library. The calibration measures Spearman ρ across all ChEMBL records. They answer different questions. The v3 pipeline preserves both verdicts and labels them clearly — this is what the methodology note is for.

---

## 3. What Works Today (end-to-end runnable)

```
                       Windows (mammal_env)                              Status
v1 pipeline: 02_fetch_targets → 13_wet_lab_shortlist                 ✅ committed, runnable
v2 ADMET cluster: 14_v2_cluster_b_admet → admet_gates.parquet         ✅ runs in ~1 min
v2 fusion (BOTH passes): 15_v2_fusion --out-suffix _{un,}calibrated   ✅ runs in seconds
v2 ESM2: 16_v2_esm2_embed → 22 cached .pt files                        ✅ runs in ~5 min
v2 Boltz: 17 smoke / 18 sweep / 19 gate / _boltzina_focused           ✅ all runnable
v3 ChEMBL SQLite mirror: chembl_downloader + chembl_sqlite.py         ✅ live (A.5 PASS)
v3 Phase A.4 backstop: 21_v3_chembl_evidence_sqlite.py --all-pairs    ✅ 99 min for 6,556 pairs
v3 Phase A.6 audit: 24_v3_audit_chembl_targets_sqlite.py              ✅ 21 ALIGNED / 1 NO_CURRENT
v3 Phase A.7 calibration: 22_v3_calibration.py                         ✅ report + per-target weights
v3 Phase D diff: 25_v3_fusion_diff.py                                  ✅ Spearman ρ = +0.994
v3 wet-lab shortlist v3: 26_v3_wet_lab_shortlist.py                    ✅ 4-cluster scorecards
A.5 smoke gate: _v3_sqlite_vs_rest_smoke.py                            ✅ 19/20 agreement


                       WSL2 Ubuntu (mammal_env)                          Status
Overnight Boltz sweep: _wsl2_boltz_full_sweep.py                      🔄 RUNNING (PID 304, 22%)
Sweep status: _wsl2_sweep_status.sh                                    ✅ available


                       WSL2 Ubuntu (txgnn_env, separate venv)            Status
Cluster C orchestrator: 23_v3_cluster_c.py                            🟡 CODE READY, NOT RUN
PrimeKG download: _wsl2_download_primekg.sh                            ✅ scripted (download pending)
Isolated venv: _wsl2_setup_cluster_c.sh                                🟡 setup script ready
```

**Reproducibility**: every artifact lands in `data/results/v2/` or `data/results/`. Configs at `configs/{thresholds,weights,weights_calibrated}.yaml`. Module API at `src/mammal_repurposing/{cluster_a,cluster_b,cluster_c,gates,fusion,provenance,pipeline,fetchers}/`. Reports at `reports/`. Sprint history is the git log; methodology note at `reports/methodology_v1.md`.

---

## 4. What's Broken / Limiting Now (post-sprint state)

The V2 plan listed library-install pain (cuequivariance, PyG, tiledbsoma, ChEMBL REST). All four are resolved. The new limiting items are different.

### 4.1 MAMMAL_ONLY_INVERTED at monoamine transporters — the central problem

**Symptom**: ρ ≤ -0.30 at SLC6A3, SLC6A2, GRIN2A, GRIN2B. Calibrated fusion down-weights to 0.30 with a `_note: "INVERTED — manual review"`. This is honest but not solved.

**Severity**: HIGH — these are the most clinically validated cognition targets in the panel.

**Remediation (v3+)**: any combination of
- **Wait for Boltz coverage** (overnight sweep) — see if structure-aware ρ rescues those targets to BOLTZ_2X_MAMMAL.
- **Diagnose via residual analysis**: which compounds in our 298-library are at SLC6A3/SLC6A2? Are they clustered in tropane scaffolds (cocaine-like)? Build a small structural-similarity audit (see §7.1).
- **Ensemble with a non-MAMMAL DTI** (TANKBind, GraphDTA, MolTrans) — see §8 stream 3.
- **LoRA fine-tune MAMMAL on a transporter-specific corpus** — risky (3-5 days, may not improve), high-leverage if it works.
- **Acquire targeted ChEMBL records** — but the n is already ≥ 25 at those two targets; data acquisition won't fix bias.

### 4.2 Boltz coverage is partial → most calibration verdicts are MAMMAL_ONLY_*

**Symptom**: only HRH3 has Boltz n ≥ 3 in Phase A.7's calibration. 17 verdicts that should be 2-cluster (`BOLTZ_2X_MAMMAL` / `MAMMAL_2X_BOLTZ` / `EQUAL_WEIGHTS` / `DE_WEIGHT_TARGET`) are forced to `MAMMAL_ONLY_*`.

**Severity**: HIGH (blocking) — calibration's whole point is comparison.

**Remediation**: wait for the overnight WSL2 sweep (~17h ETA). Then re-run Phase A.7 → C → D. The expected delta is large; most current verdicts will move.

### 4.3 Cluster C (PrimeKG + TxGNN) coded but not run

**Symptom**: `scripts/23_v3_cluster_c.py` exists; the WSL2 `txgnn_env` venv setup is scripted but not executed; PrimeKG download (~1.4 GB) is scripted but not done.

**Severity**: MEDIUM. Cluster C absence means the 4-cluster RRF is currently a 3-cluster (MAMMAL + Boltzina + ADMET, with TxGNN/PrimeKG defaulting to NaN in the rrf module). The mechanism-based disagreement archetypes the v2 design promised are mostly dormant.

**Remediation**: scheduled as Phase B in the original V3_ATTACK_PLAN. Blocked on the overnight Boltz sweep finishing so the GPU is free. After that:
1. `bash scripts/_wsl2_setup_cluster_c.sh` (~30 min)
2. `bash scripts/_wsl2_download_primekg.sh` (~10 min)
3. `python scripts/23_v3_cluster_c.py` (~1 hour for 298 compounds)
4. Re-run fusion both passes; Cluster C parquet now feeds 4 rankers.

### 4.4 HRH3 dominates the calibrated shortlist (single-target lock-in)

**Symptom**: 23 of 25 top compounds in `wet_lab_shortlist_v3.md` show `mammal_best_target = HRH3` because HRH3 was the only target NOT down-weighted by the v1 calibration. The shortlist degenerated into an "H3 antagonist + ADMET-clean" filter.

**Severity**: MEDIUM. The shortlist is still defensible (pitolisant, modafinil, methylphenidate, donepezil, galantamine all in top-25; reasonable cognition priors) but it's not surfacing diversity across mechanism classes.

**Remediation**:
- After Boltz coverage expands, more targets get default weights, dominance loosens.
- **Selectivity-aware reranking** (see §8 stream 1) — surface compounds that have a *flat* profile across the panel (polypharm) vs `HRH3-dominant` ones.
- **Per-mechanism-class top-N**: instead of one big top-25, surface top-5 per mechanism (cholinergic, glutamatergic, dopaminergic, etc.) so monoaminergic targets aren't silenced.

### 4.5 The 22-target panel itself is curated, not validated against cognition GWAS

**Symptom**: the panel came from the research doc + clinical priors. We don't know how it overlaps with:
- UK Biobank cognition GWAS hits (e.g., Lam et al. 2017 *Cell Reports*)
- AD genetics (e.g., Bellenguez et al. 2022 *Nat Genet*, 75 loci)
- iPSC-neuron functional screens (e.g., the FENS / iScience datasets)
- Allen Brain Atlas expression patterns for "cognition-relevant" cortical regions

**Severity**: MEDIUM. The panel might be missing high-value targets (e.g., APOE, SORL1, TREM2 for AD-cognition; KCNN3 for synaptic plasticity; CTNNB1 / WNT signalling for hippocampal LTP) or over-weighting clinical-precedent targets that aren't actually cognition-causal.

**Remediation**: see §8 stream 5 + §7.3.

### 4.6 No behavioural anchor — cognition treated purely as binding proxy

**Symptom**: the cognition virtual phenotype anchor is 5 disease EFO IDs in TxGNN. There's no link from any predicted score to a measurable cognitive endpoint (working memory delta, processing speed, attention span). The pipeline says nothing about Roberts CA 2020's SMD = 0.21 ceiling beyond mentioning it in the methodology note.

**Severity**: STRUCTURAL — this is the fundamental honesty limit of the pipeline.

**Remediation**: any quantitative endpoint requires either (a) a target → endpoint mapping (heuristic, possible now), or (b) a clinical/behavioural meta-analysis backbone (Frumkin AM 2022 *J Psychopharm* on stimulants; the cognitive aging consortium; ABCD Study data). See §8 stream 5.

### 4.7 LambdaMART is now eligible but not promoted

**Symptom**: `fusion/lambdamart.py` is gated on ≥20 labels per target. We now have 275 CORROBORATED records globally; per-target it's often ≥20 (especially at the well-covered targets like ACHE 19,486 activities, SLC6A3 13,293, HCRTR1 24,283).

**Severity**: LOW (RRF works) — but LambdaMART would let calibrated weights become continuous per-(target, compound) rather than piecewise per-target.

**Remediation**: scheduled v3 work — binarise pchembl at ≥ 6.0 vs < 6.0, train LightGBM LambdaRanker on (compound morgan FP + ESM2 embedding + per-cluster scores) features. Per-target rankers; voted into the fusion.

### 4.8 Sweep runs use the GPU exclusively — Cluster C work serialised

**Symptom**: Boltz sweep is single-GPU on RTX 5070 12 GB. Cluster C (PrimeKG load = ~2 GB RAM; TxGNN inference = ~1-2 GB VRAM) can't run in parallel without contention.

**Severity**: LOW (scheduling, not correctness).

**Remediation**: scheduling — run Cluster C in the gaps between Boltz batches, or split the GPU via MIG (not supported on RTX 5070), or accept overnight-serialisation. For v3+: cloud GPU burst (1 × A100 for the Boltz sweep, freeing the 5070 for parallel Cluster C work) is a 1-day-$30 option if priority demands it.

### 4.9 Reports are stale-vs-snapshot — no provenance hash in the markdown

**Symptom**: each report says "Generated by `scripts/22_v3_calibration.py`" but doesn't pin the input parquet hashes or the git SHA at generation time. When we re-run after Boltz finishes, the old reports become silently stale.

**Severity**: LOW. Easy to fix in the report renderers — add `_meta: {git_sha, generated_at, input_hashes}` blocks.

**Remediation**: small follow-up; one PR across all report-writing scripts.

---

## 5. Throughput Bottlenecks (revised, post-A.4)

| Operation | Per-call wall-clock | Total grid | ETA |
|---|---|---|---|
| MAMMAL DTI inference | 100 ms (batch 8) | 6,556 pairs | ~10 min |
| ESM2-650M embedding | 5 s (single target) | 22 targets | ~5 min |
| ADMET-AI 41 endpoints | 100 ms (CPU) | 298 compounds | ~30 s |
| Boltz-2 affinity, WSL2 + cuequivariance kernels | ~67 s | 1,165 pairs (current sweep) | ~22 h (running) |
| Boltz-2 affinity, Windows (no kernels) | 150 s | n/a (deprecated path) | — |
| ChEMBL evidence backstop (SQLite mirror) | ~900 ms incl. RDKit InChIKey | 6,556 pairs | ~99 min (done) |
| ChEMBL per_target_pchembl_records query | ~5 s | 22 targets | ~2 min |
| Phase A.7 calibration end-to-end | ~3 min | one call | ~3 min |
| Phase C fusion (both passes) | ~5 s each | two calls | ~10 s |
| Phase D diff | ~1 s | one call | ~1 s |
| Wet-lab shortlist v3 render | ~1 s | one call (top 25) | ~1 s |
| RRF fusion + provenance | 200 ms | one call over full grid | ~1 s |

**Remaining slow operation**: the WSL2 Boltz sweep. At 67 s/pair on a 12 GB RTX 5070, expanding from 1,165 to the full 6,556-pair grid would take ~5 days. Realistic options:
1. **Stay scoped** — the 1,165 pairs is top-N per ADMET-surviving compound (10 best targets per compound). That's the right scope.
2. **Speedup** — Boltzina-Vina-only (cluster A.2 variant) reduces per-pair to ~10-15 s (research doc §1.2, never implemented). Would compress the sweep to ~3 hours.
3. **Cloud burst** — A100 80GB at ~$2/hour finishes the same 1,165 sweep in ~6 hours, $12. Trivial budget.

**Phase A.4 was a 470× speedup vs the original REST grinder** (78h projected → 99 min actual). Storage cost: ~13 GB for the extracted SQLite mirror. Pays for itself permanently.

---

## 6. Component Status Matrix (one-pass reference)

### 6.1 Library installs (post-v3)

| Library | Windows native | WSL2 mammal_env | WSL2 txgnn_env | Notes |
|---|---|---|---|---|
| torch 2.12 nightly cu128 | ✅ | ✅ | ❌ | mammal_env stays on 2.12; txgnn_env must use 2.7 for PyG |
| torch 2.7.0+cu128 | — | — | ✅ | txgnn_env requirement (PyG wheel target) |
| biomed-multi-alignment (mammal) | ✅ | ✅ | — | |
| admet-ai | ✅ | ✅ | — | |
| boltz | ✅ slow | ✅ fast | — | WSL2 + cuequivariance native kernels work on Blackwell |
| cuequivariance-ops-torch-cu12 | ❌ Linux only | ✅ | — | confirmed via overnight sweep |
| chembl-downloader | ✅ | ✅ | — | `latest()` returns version string; use `download_extract_sqlite()` for Path |
| pyg_lib / torch_scatter / etc. | ❌ no nightly wheel | ❌ ABI mismatch with torch 2.12 | ✅ | Separate venv required |
| txgnn | — | — | ✅ | install via `pip install git+https://github.com/mims-harvard/TxGNN.git` |
| igraph | ✅ | ✅ | ✅ | PrimeKG loader, 10× faster than networkx for 4M edges |
| rdkit | ✅ | ✅ | — | for InChIKey calc in calibration join + SMILES→InChIKey |
| networkx | ✅ | ✅ | ✅ | only used for PrimeKG fallback |
| lightgbm | ✅ | ✅ | — | LambdaMART; not yet promoted |

### 6.2 V3 modules (additions on top of v2)

| Module | LOC (approx) | Status | Test coverage | Notes |
|---|---|---|---|---|
| `fetchers/chembl_sqlite.py` | ~330 | LIVE | smoke (Phase A.5) | salt-form parent resolution via molecule_hierarchy JOIN |
| `cluster_c/primekg.py` (real) | ~190 | LIVE (code) | none | igraph-backed; PPR scoring; deferred until Phase B runs |
| `cluster_c/txgnn.py` (real) | ~120 | LIVE (code) | none | 5-anchor weighted mean over EFO union |
| `fusion/rrf.py` v3 | ~165 (+25 vs v2) | LIVE | implicit | new `per_target_weights` kwarg |
| `provenance/tracker.py` v3 | ~140 (+30 vs v2) | LIVE | implicit | defensive on Boltzina/TxGNN/KG column-name variants |
| `scripts/21_v3_chembl_evidence_sqlite.py` | ~120 | LIVE | A.4 complete | SQLite backstop, --all-pairs flag |
| `scripts/22_v3_calibration.py` | ~340 | LIVE | none | linchpin; refined verdict matrix in v3 |
| `scripts/23_v3_cluster_c.py` | ~170 | CODE, NOT RUN | none | WSL2 txgnn_env entry point |
| `scripts/24_v3_audit_chembl_targets_sqlite.py` | ~180 | LIVE | A.6 done (21 ALIGNED) | clean SQLite re-run of T1 audit |
| `scripts/25_v3_fusion_diff.py` | ~190 | LIVE | none | Phase D output |
| `scripts/26_v3_wet_lab_shortlist.py` | ~210 | LIVE | none | per-cluster scorecards for top-N |
| `configs/weights_calibrated.yaml` | 67 lines | GENERATED | regenerate when Boltz expands | per-target overrides from Phase A.7 |

---

## 7. Architecture Enhancements Needing Research

These are the items the user asked for. Each has: concept, why now, research gap that needs commissioning, effort tier, dependencies.

### 7.1 Per-target failure-mode diagnostic ("why is MAMMAL inverted at DAT?")

**Concept**: a notebook + small module that for each MAMMAL_ONLY_INVERTED target produces:
- A scatter plot of MAMMAL pKd vs ChEMBL pchembl for every compound with both
- A structural-similarity audit (Tanimoto on Morgan FP + **Murcko scaffold cluster assignment**) showing how our 298-library distributes across BindingDB's typical chemical clusters for that target
- An ESM2-residue-attention heatmap on the target's substrate-binding region — is MAMMAL "looking" at the right pocket?
- A pchembl-distribution comparison between ChEMBL's full record and BindingDB's filtered subset for that target

**Strong prior hypothesis** (from research stream A): the failure is **manifold mismatch from training-set saturation**, not epistemic uncertainty.

For DAT (SLC6A3) specifically, BindingDB is saturated with two scaffold families:
- 3β-(4-substituted-phenyl)tropane-2β-carboxylic acid methyl esters — Carroll 2004 *J Med Chem* / PubMed 15566309: β-CIT, β-CFT, RTI-series cocaine analogs
- Methylphenidate / phenidate analogs (Wikipedia "methylphenidate analogues" lineage)

When BindingDB compresses high-affinity DAT binders into one scaffold cluster, MAMMAL learns "tropane-like → bind" but **loses pKd dynamic range within that cluster**. ChEMBL pchembl distinguishes within-cluster potency the model can't see. Result: rank inversion.

**Concrete diagnostic gate**: compute Murcko scaffolds for the 26 SLC6A3 compounds in our DTI grid. If **< 30% match the tropane/phenidate family**, MAMMAL is extrapolating outside its training manifold — the anti-correlation is a domain-shift signature → ensemble with MMAtt-DTA (whose training distribution differs) is the fix. If **> 70% match**, the problem is rank-resolution loss within a saturated cluster → fine-tuning or scaffold-stratified data acquisition is the fix.

**Why now**: this single diagnostic decides between LoRA fine-tuning (§7.6) vs cross-DTI ensembling (§7.7) vs targeted data acquisition (§7.13).

**Closest published interpretability analogues**: Lin ESM3 2024; Sgarbossa 2024 on per-residue attention; CryptoBench (Pólya 2024 *Bioinformatics* 41) for pocket-region attribution.

**Effort**: 1-2 days (most components exist; RDKit Murcko + scatter + ESM2 attention extraction).

**Dependencies**: none — uses existing parquets.

### 7.2 Sign-correction layer (lossy but explicit)

**Concept**: an optional `--invert-flagged-targets` mode for fusion that multiplies MAMMAL's per-target ranks by -1 at the 4 INVERTED targets. Outputs into `rrf_ranking_calibrated_sign_corrected.parquet` clearly separate from the honest pass. The methodology note explicitly does NOT auto-invert because the sign could flip with more data, so we make it a one-flag deliberate choice.

**Why now**: gives us a "what if we *did* trust the negative correlation?" stress-test for the wet-lab shortlist. If sign-correction surfaces drastically different top-25, that's information about how much the INVERTED targets are driving current ranking decisions.

**Research gap**: minor — sign-correction is well-understood; the framing is calibration-policy, not novel.

**Effort**: half-day.

**Dependencies**: none.

### 7.3 Cognition panel expansion: GWAS-anchored target prioritisation

**Concept**: pull the top 100 cognition GWAS hits (UK Biobank, Davies et al. 2018 *Mol Psychiatry* PMID 29186139; Lee et al. 2018 *Nat Genet* PMID 30038396 on educational attainment as proxy; Hill et al. 2019 *Nat Commun* on general cognition). Score panel coverage. Surface high-evidence non-panel targets and add them in v4 with confidence stratification.

**Why now**: the current 22-target panel was curated from clinical precedent. GWAS gives an orthogonal target-discovery anchor that biases less toward "drugs that already exist."

**Research gap**: GWAS → druggable target translation. Pharmaprojects-style annotations are paid; OpenTargets has free target-tractability scores. The GWAS-Catalog has open EFO-mapped hits.

**Effort**: 2-3 days (data pull + scoring + integration into `targets_seed.csv`).

**Dependencies**: requires a "tractability" annotation layer — could pull from OpenTargets free GraphQL.

### 7.4 Selectivity panel scoring layer

**Concept**: for each compound, compute a "selectivity profile" across the 22 targets and surface selectivity *categories* (mono-selective, dual-selective, polypharmacological, panel-flat).

**Concrete metrics** (informed by research stream B — there is no consensus CNS-panel equivalent of the kinome's Karaman/Davis percent-control, but the toolkit transfers cleanly):

| Metric | Citation | What it gives us |
|---|---|---|
| **Selectivity entropy** | Uitdehaag & Zaman 2011 *PLOS ONE*, PMC3100252 | Shannon entropy over normalized 1/Kd vector. Low = selective, high = promiscuous. One numpy function. |
| **Partition Index** | Cheng et al. 2010 *J Med Chem*, doi 10.1021/jm100301x | Thermodynamically-grounded fraction of compound bound to target X at equilibrium across the panel. Better than entropy when one "intended" target exists. |
| **Gini + S(10x) tandem** | Graczyk-style scorecard, BMC Bioinformatics 2017 doi 10.1186/s12859-016-1413-y | Two numbers: Gini coefficient of panel affinities + fraction of panel within 10× of best target. Industry-friendly; reviewers expect both. |
| **KISS-CL** | *Nat Commun* 2025 article 65869-8 | Contrastive embedding trained to separate same-scaffold compounds with divergent panel profiles. SOTA on kinase panel selectivity. Pilot as re-scorer on MAMMAL embeddings. |

**Cleanest published "selectivity > target-X" definition**: S(10x) ≥ 0.30 with target-X ranked #1 by pKd, **and** entropy ≥ median of polypharmacology training set (Uitdehaag 2011 conventions). Reproduce with `scipy.stats.entropy` + a one-line fraction.

**Why now**: cognition drugs split cleanly along this axis. Donepezil is mono-selective (AChE). Modafinil is dual (DAT + HRH3). Galantamine is polypharm (AChE + α7 PAM + others). The current pipeline collapses this into `mammal_best_target` and `mammal_polypharm_n`; the structured selectivity profile lets the methodology note say something quantitative about diversification + lets §8.1 (multi-class top-N) surface mono- vs poly-selective top compounds separately.

**Effort**: 2 days (compute + integrate into wet-lab shortlist v3). KISS-CL pilot adds 3-5 days.

**Dependencies**: none — uses DTI grid as-is.

### 7.5 Pocket-conditioned Boltz inference (allosteric vs orthosteric)

**Concept**: before / after running Boltz on a (target, compound) pair, identify which pocket the predicted pose lands in. Surface "binds orthosteric vs allosteric vs cryptic vs surface artefact" as a provenance flag.

**Concrete approach** (research stream B):

1. **First-pass orthosteric + cryptic detection**: P2Rank (Krivak & Hoksza, github.com/rdk/p2rank, Java single-jar) for orthosteric pockets, PocketMiner (Meller 2023 *Nat Commun* article 36699-3) for cryptic. Tag each Boltz pose with the nearest predicted pocket centroid.

2. **CryptoBench classifier** (Pólya 2024 *Bioinformatics* 41 article btae745) — ESM2-based cryptic-site classifier that beats P2Rank and PocketMiner. We already have ESM2-650M embeddings cached (`data/cache/esm2/`); retraining the head against CryptoBench labels is ~2h on the 5070. Bonus: lets us add cryptic-site flags to the panel without re-embedding.

3. **AlphaFold3 / Boltz-2 cofolding for cryptic discovery** (Bryant 2023 *Nat Commun* PMC10373493): cofold ligand + target; cryptic pockets that only open upon ligand binding are revealed. Use as **second opinion** for any compound where MAMMAL pKd and Boltz affinity diverge by > 1.5 log units — i.e. the disagreement archetypes already tracked in `provenance/disagreement_report.py`.

4. **Provenance flag schema** (add to wet-lab shortlist v3):
```
pocket_class ∈ {orthosteric, allosteric_known, allosteric_putative,
                cryptic_predicted, surface_artifact, no_pocket_match}
```
`surface_artifact` and `no_pocket_match` are **auto-demoted** in the wet-lab shortlist.

**Why now**: the Phase 0.5 CHRNA7 rescue showed Boltz finds the PAM site implicitly — but we can't *prove* it's the PAM site without pocket annotation. For wet-lab handoff, "this compound is predicted to bind the allosteric site, not the orthosteric site" is a clinically actionable distinction. This is the single highest-impact addition for downstream trust — orthosteric vs allosteric is the question wet-lab will ask first.

**Effort**: 3-5 days (P2Rank install + pocket DB curation for 22 targets + CryptoBench classifier + Boltz integration).

**Dependencies**: P2Rank (open-source, Java; runs in WSL2). Per-target manual pocket DB is curation work (PDB lookup per target's known orthosteric + any annotated allosteric sites). CryptoBench needs the ESM2 cache (✅ exists).

### 7.6 LoRA fine-tune MAMMAL on cognition-DTI corpus

**Concept**: per the V3 attack plan T11. Use the now-validated 275 CORROBORATED ChEMBL pairs + the 22 panel targets to construct a fine-tuning corpus. LoRA on MAMMAL's DTI head + encoder. ~3-5 day training budget on RTX 5070.

**Why now**: the Phase A.7 calibration finally gives us the metric we'd train against. Without per-target Spearman ρ we couldn't have measured improvement.

**Concrete approach (informed by research stream A)**: copy the **BALM** pattern (Gorantla et al., 2024, bioRxiv 2024.11.01.621495; github.com/meyresearch/BALM, Apache-2.0). BALM freezes ESM-2 + ChemBERTa-2 and injects two PEFT adapters — **LoKr** (Kronecker LoRA) on the protein tower and **LoHa** (Hadamard low-rank) on the ligand tower — into a shared cosine-similarity space. It beats full fine-tuning on every BindingDB split.

For MAMMAL specifically:
- Keep T5 encoder-decoder frozen; attach LoRA adapters at **decoder cross-attention Q/V projections** (where target+compound tokens interact)
- **Rank r=8 default**, sweep r ∈ {4, 8, 16} on held-out per-target ρ
- Schmirler et al. 2024 *PNAS* on protein-LM PEFT, ESM-LoRA-Gly 2025, SeqProFT 2024 all converge on r=4-8 being the sweet spot
- One MAMMAL LoRA reproduction exists in the wild (`Kymi808/mammal-lora-bbbp` GitHub; matches Shoshan 2025 BBBP within 3σ) but is poorly documented

**Data-floor estimates** (no SOTA — open question; published anchors):
- FS-Mol (Stanley 2021): 32-64 active examples per task lifts ROC-AUC by 5-10 points
- BALM: +0.05 to +0.15 Pearson on ~3K BindingDB-Kd pairs per target class
- **Pragmatic floor: ~50 records/target to move ρ off noisy baseline; ~200/target to push from ρ ≤ -0.30 into ρ ≥ +0.30**. Our 275 CORROBORATED + ~1,400 truth records on the top 4 worst-calibrated targets is right in the meaningful band.

**Catastrophic forgetting mitigation** (transferable from LLM literature):
- LoRA inherently minimises forgetting (base weights frozen)
- **Replay**: mix 10-30% of original BindingDB samples into each fine-tune batch (Aljanabi & Khaleel 2025, arXiv:2501.13669)
- **L2-SP** / **EWC** (Kirkpatrick 2017) for full fine-tuning paths
- **Adapter routing per target class** (Qiao 2024 NeurIPS, "Learn more, but bother less") — different LoRA per protein superfamily, gated at inference. Maps directly to "different adapter for SLC6A* than for GRIN2*"

**ChEMBL corpus prep** (consolidated from Bento 2020 *J Cheminform* 12:51; Landrum 2024 update):
1. Standardise: RDKit `chembl_structure_pipeline.standardize_mol` → `get_parent_mol` (salt strip → largest fragment → neutralise charge → canonical tautomer)
2. Filter: `assay_type='B'`, **`confidence_score ≥ 8`** (bump from our current ≥ 7 to remove "homologous protein" mappings), `standard_relation='='` only (drop censored), pchembl not null
3. Aggregate: mean pchembl per (parent InChIKey × target ChEMBL_ID); drop pairs with intra-pair stdev > 1.0 log unit
4. Threshold-rebalance by potency deciles when splitting (BASE webservice, Song 2024 *BMC Bioinformatics* 25)
5. **Murcko-scaffold split** train/cal/test (not random) + a time-split holdout for temporal drift

**Effort**: 5-10 days end-to-end (corpus prep + LoRA training + per-target re-calibration + comparison).

**Dependencies**: T1 audit complete (✅ done). Tightening `confidence_score` filter from 7 to 8 is a one-line change in `chembl_sqlite.py`; check first whether it leaves enough rows.

### 7.7 Cross-DTI ensemble (orthogonal voters into RRF)

**Concept**: add a second + third DTI head as additional rankers in the fusion. The calibration linchpin's finding (MAMMAL ρ ≤ -0.30 at 4 targets) is exactly the case where ensembling helps most. If MAMMAL is biased a particular way at SLC6A3, a graph-attention or contrastive DTI may not share that bias.

**Concrete candidate panel** (all open-source, permissive license, fits 12 GB RTX 5070 at inference):

| Method | Citation | License | Inductive bias | Why orthogonal to MAMMAL |
|---|---|---|---|---|
| **MMAtt-DTA** ★ | Schulman 2024 *Bioinformatics* 40:btae496 | MIT | Attention-Transformer + LASSO-selected multimodal descriptors | **Reports ρ > 0.72 on transporter superfamily** — directly addresses our SLC6A3/SLC6A2 failure |
| **PSICHIC** | Koh 2024 *Nat Mach Intell* 6:673 | MIT | PhysicoChemical GNN; sequence-only, no 3D; learns interaction fingerprints | GNN vs MAMMAL's T5 — different architecture entirely |
| **BALM** | Gorantla 2024 bioRxiv 2024.11.01.621495 | Apache-2.0 | ESM-2 + ChemBERTa-2 with LoKr/LoHa PEFT, cosine-similarity readout | Separately-pretrained protein/ligand towers, no joint-pretraining bias |
| **ConPLex** | Singh 2023 *PNAS* 120:e2220778120 | MIT | ProtBERT + Morgan FP, contrastive with decoys; <2 GB VRAM | Contrastive-trained, decoy-aware |
| **Komet** | Guichaoua 2024 *J Chem Inf Model* 64:6938 | BSD-3 | Kronecker chemogenomic kernel with Nyström approximation; 1M+ LCIdb pairs | Classical kernel — robust to MAMMAL's failure modes |
| **GEMS** | Lemm 2024 bioRxiv 2024.12.09.627482 | Apache-2.0 | GNN + LM embeddings on PDBbind CleanSplit | Structure-aware, uses 3D pocket info |
| **DrugCLIP** | Gao 2024 arXiv:2310.06367 | MIT | Contrastive ligand-pocket co-embedding; retrieval-trained | Ranks via dense retrieval, not regression |
| **MGNDTI** | Peng 2024 *JCIM* PMID 39137398 | MIT | Multimodal representation with explicit gating | Gating vs MAMMAL's fused-embedding |
| **PLAPT** | Rose 2024 bioRxiv 2024.02.08.575577 | MIT | ProtBERT + ChemBERTa + FC head; ~1 GB VRAM | Simplest baseline, no domain-specific pretraining |

★ = Top recommendation. **MMAtt-DTA is the single strongest single addition if transporter ρ is the metric we want to fix.**

**MoE gating to route the ensemble**:
- **EnsDTI** (Park 2024, bioRxiv 2024.08.06.606753): freezes K DTI experts, trains a softmax gating network on prediction probabilities. Error < 1% when gating confidence > 60%. Exactly the architecture we want — wrap MAMMAL + MMAtt-DTA + PSICHIC + fallback, learn per-target routing
- **MoSE** (Yan 2025 arXiv:2503.15796 AAAI 2025): intrinsic/extrinsic experts with adaptive sample-level fusion; +53% on data-scarce benchmarks

**Per-compound OOD detection**:
- **eMOSAIC** (Zhang 2025 *Nat Mach Intell* article s42256-025-01151-2; bioRxiv 2024.01.05.574359): Mahalanobis distance in embedding space scores OOD-ness on multi-target binding models. Use as a confidence flag in the wet-lab shortlist — if compound is OOD vs all training distributions, flag for manual review

**Effort**: 2-3 days for MMAtt-DTA (first additional ranker, pip-installable), 1-2 days each for PSICHIC + BALM. RankerInput is already templated so each new method drops in without fusion changes.

**Dependencies**: pip-only for MMAtt-DTA and PSICHIC. BALM needs the PEFT install but no model surgery.

### 7.8 Compound-library expansion via REINVENT / POLYGON / PILOT

**Concept**: use the calibrated panel as a fitness function for generative chemistry. The current 298-compound library is a known-nootropic prior; generation breaks out of that prior, conditioned on the calibrated multi-target objective.

**Concrete platforms** (research stream B):

| Platform | Citation | Best fit for | Notes |
|---|---|---|---|
| **REINVENT 4** | AstraZeneca, Loeffler 2024 *J Cheminform* article PMC10882833 | General-purpose multi-property RL | De-facto open-source standard; scoring function is user-supplied via custom Python component; multi-target conditioning is standard recipe |
| **POLYGON** | *Nat Commun* 2024 article 47120-y | Dual-target ligand RL | 82.5% polypharmacology recall; wet-lab validated on MEK1/mTOR. Closest published precedent for our 22-target objective |
| **PILOT** | *Chem Sci* 2024 article d4sc03523b | Pocket-conditioned diffusion | Equivariant diffusion, multi-objective via importance sampling. Better than REINVENT when structures exist (we do, via Boltz) |
| **Pocket2Mol-RL** | Deargen 2024, github.com/deargen/Pocket2Mol_RL_public | Custom-reward RL on pocket-aware generator | Minimal code to attach our fused score |

**Realistic wall-clock on RTX 5070 (12 GB)**:
- REINVENT 4 RNN sampling: ~200 mol/s
- Panel scoring (MAMMAL inference 22 targets): ~50 mol/s after batching
- ADMET-AI: ~30 mol/s
- **Generating 1,000 candidates that beat top-25 RRF**: 4-8 hours per RL epoch, 3-5 epochs needed = **one weekend**
- PILOT diffusion: 5-10× slower per molecule but better pose-quality. Budget overnight for 1,000.

**Why now**: the v3 calibration finally gives us a confidence-weighted panel score. Random library expansion is wasteful; generative-conditional sampling targets the right region. Also: even the current top-25 is mostly known nootropics — generative would surface scaffolds we've never considered.

**Effort**: 7-14 days (REINVENT integration + multi-target conditioning + ADMET-in-loop + re-scoring + IP filter).

**Dependencies**: pip-installable. Tie IP filter into §8.3 patent cross-reference so generated candidates with strong existing IP are auto-deprioritised before display.

### 7.9 Behavioural / EEG / fMRI biomarker integration as a *target prior* (Cluster D)

**Concept**: instead of pretending the pipeline predicts behavioural outcomes, use behavioural/EEG/fMRI biomarker data to **prioritise targets** in the panel. The pipeline still predicts binding; the *weighting* of targets reflects biological plausibility for cognition.

**Concrete datasets + recipes** (research stream B):

1. **Mansuri et al. 2024** (PMC10941541) — identified 41 individual cortical genes whose AHBA expression maps spatially correlate with intelligence-related fMRI. Drop our 22 panel targets onto this map → cognition-anchored prior on which targets to upweight. Currently targets are weighted by ChEMBL data density only; this gives a biological-relevance prior.

2. **UK Biobank multimodal MRI marker of cognition** (eLife reviewed preprint 2025, elifesciences.org/reviewed-preprints/108109) — 48% of g-factor variance explained by structural + functional MRI; published target-set imputation directly usable.

3. **Open Targets Genetics + UK Biobank GWAS cognition portal** (genetics.opentargets.org) — L2G score for each panel target on cognitive-trait GWAS hits; penalise targets with strong loss-of-function intolerance + cognitive impairment phenotype. Single GraphQL endpoint.

4. **Lit-OTAR** (Sangkuhl 2024, PMC11978389): 48M target-disease-drug associations extracted from full text. Use for both target weighting AND per-compound prior-evidence flag.

5. **Hansen 2020 AHBA + Neurosynth pipeline** (biorxiv 2020.07.16.203026) — well-cited spatial-correlation recipe, ~200 lines of Python.

**Concrete addition**: new **Cluster D = "neurobiological prior"** = per-target weight ∈ [0.5, 1.5] derived from (AHBA spatial corr to cognition fMRI + OT Genetics L2G + lesion-impairment penalty). Feeds into RRF as target weights, multiplied into the existing per-target calibrated overrides from Phase A.7.

**Why now**: this is the principled way to add a behavioural signal without overclaiming. Mansuri's 41 genes are public and downloadable today.

**Effort**: 5-10 days (AHBA download + spatial correlation + Open Targets Genetics pull + integration).

**Dependencies**: HTTP + abagen library for AHBA + Neurosynth API.

### 7.10 LLM-agent prioritisation overlay

**Concept**: wrap the entire pipeline as MCP tools an LLM agent can call. The agent receives a freeform query ("find compounds that enhance working memory via prefrontal-cortex-dependent mechanisms without DAT activity") and orchestrates appropriate filters + rankers + literature retrieval. Output: a 1-2 paragraph candidate brief with citations.

**Why now**: the pipeline has matured to the point where it has too many knobs for non-experts to use. An LLM agent that knows the methodology note can serve as the natural-language frontend.

**Templates**: AgentD blueprint (Liu 2025 *JCIM* doi 10.1021/acs.jcim.5c02454); PharmaSwarm (Chen 2025 arXiv:2504.17967) — both benchmarked Claude-3.7-Sonnet > GPT-4o on tool-orchestration for drug discovery.

**Effort**: 7-14 days (MCP server + tool surface + agent prompt + integration testing).

**Dependencies**: MCP server scaffolding; methodology note as system prompt anchor; the calibration_report.md as a "what targets you can / cannot trust" runbook.

### 7.11 Isotonic per-target post-hoc calibration (the elegant sign-correction)

**Concept**: instead of multiplying inverted targets' MAMMAL scores by -1 (lossy + sign-fragile), fit a **monotone isotonic regression** of (MAMMAL pKd → ChEMBL pchembl) per target on the calibration set. The fit is monotone but flexible — if the true MAMMAL-vs-truth map is monotonically *decreasing* on DAT, the isotonic fit recovers that ordering and the sign-flip falls out as a property of the fit, not an arbitrary policy.

**Concrete approach** (research stream A):
- **Toplak 2020** (*J Chem Inf Model* 60:3829, PMID 32865408) compared Platt, isotonic, and Venn-ABERS on 40M ligand-target pairs across 2,112 targets. **Venn-ABERS won by Brier score across all base models and CV splits**
- Isotonic is right for our case: monotone but flexible, doesn't require parametric form
- **Key caveat**: with n=8 (GRIN2A) to n=26 (SLC6A3), classical isotonic overfits. Solutions:
  - Bayesian isotonic prior (Neelon & Dunson 2004)
  - Pooled-target hierarchical isotonic (no published DTI implementation — open research)
  - **Beta-calibration** (Kull 2017) — parametric middle ground that handles ≤ 30 calibration points more gracefully

**Why now**: the cleanest, lowest-risk path to "use MAMMAL at INVERTED targets in a defensible way." Outputs `data/calibration/per_target_isotonic_*.pkl`; fusion has a `--use-isotonic-calibration` flag that swaps raw MAMMAL pKd for the isotonic-transformed value.

**Effort**: 2-3 days (scikit-learn `IsotonicRegression` + Bayesian beta-calibration via `betacal` package + per-target fit + integration).

**Dependencies**: existing calibration parquet; sklearn + betacal (pip).

### 7.12 Conformal prediction per-target gating (route MAMMAL out at bad targets)

**Concept**: rather than fight MAMMAL's per-target failures, **gate** MAMMAL out at the targets where its conformal-predicted intervals are uninformative. The fusion's per-target weights become a *function* of the prediction uncertainty, not a static yaml.

**Concrete approach** (research stream A):
- **Rakhshaninejad 2025** (arXiv:2505.18890) — cluster-conditioned conformal prediction on DTI; nonconformity-clustered CP gives the tightest subgroup intervals on KIBA splits
- **Bosc 2023** (*J Cheminform* 15:79, PMC10457664) — dynamic applicability domain / per-target local CP, construct calibration set per query from k-nearest neighbors in chemical space
- Per-target Expected Calibration Error (ECE) tracked across new data → gate MAMMAL out automatically at targets where ECE exceeds threshold
- Routes to fallback (cross-DTI ensemble §7.7) when MAMMAL is uninformative

**Why now**: gives the pipeline *automatic* per-target trust decisions — when overnight Boltz data arrives, the conformal layer notices that MAMMAL's uncertainty interval at SLC6A3 is much wider than at DRD1 and routes accordingly. No human re-calibration step needed.

**Effort**: 5-7 days (per-target CP fit + ECE monitoring + gating logic + integration into fusion).

**Dependencies**: `nonconformist` or `mapie` pip libraries; existing calibration parquet.

### 7.13 Scaffold-aware active learning for calibration-set expansion

**Concept**: if §7.1 diagnostic confirms tropane-scaffold saturation as the cause of MAMMAL_ONLY_INVERTED at DAT, the fix is **not** more compounds from the same scaffold — it's compounds from *underrepresented* scaffolds in the training distribution. Use scaffold-stratified Tanimoto-MaxMin acquisition (not BatchBALD) to add 50-200 compounds to the calibration set per worst target.

**Concrete approach** (research stream A):
- **Murcko-scaffold cluster** the 26 SLC6A3 compounds + all ChEMBL records for the target
- **MaxMin pick** (RDKit `SimDivFilters.MaxMinPicker`) across underrepresented scaffolds
- Within each scaffold cluster, **GP-UCB on (MAMMAL pKd, fingerprint) features** to pick informative points
- Re-fit isotonic (§7.11) or residual head (§7.14) per target

**Bailey 2024** (eLife 12:RP89679; biorxiv 2023.07.26.550653) benchmarks: batched greedy + diversity-weighted uncertainty beat BatchBALD on wall-clock at equivalent regret in drug-discovery campaigns. **Holzmüller `bmdal_reg`** (github.com/dholzmueller/bmdal_reg) is the canonical regression AL library (BAIT, LCMD, BALD-regression).

**Bioactivity similarity index (Skinnider 2025, *Front Bioinform* 5 doi 10.3389/fbinf.2025.1695353)** outperforms Tanimoto for chemistry-relevance ordering.

**Why now**: when the diagnostic in §7.1 fires, this is the cheapest data-acquisition strategy — uses *existing* ChEMBL records rather than asking for new experiments.

**Effort**: 3-5 days (scaffold clustering + MaxMin picker + GP-UCB + re-calibration loop).

**Dependencies**: RDKit (✅ installed); ChEMBL SQLite (✅ live).

### 7.14 Residual-correction head (XGBoost meta-ranker)

**Concept**: instead of fine-tuning MAMMAL, train a small residual head on `(MAMMAL_pKd, Morgan FP, target one-hot, ESM2 embedding) → ChEMBL pchembl`. The model learns per-target offset + per-fingerprint-region correction in one model. Cheapest path; no MAMMAL surgery.

**Concrete approach**:
- XGBoost or 2-layer MLP
- Per-(target, compound) features: MAMMAL pKd + Boltz binder_prob + ADMET-AI score + ESM2 cosine to nearest panel-target paralog + Morgan FP (radius 2, 2048 bits)
- Target: ChEMBL pchembl (from Phase A.4 backstop)
- Train on the 275 CORROBORATED + 18 AMBIGUOUS rows; cross-validate per target
- Output a *residual* added to the fusion's RRF score

**Pattern reference**: no published DTI residual-correction head, but well-established in molecular property prediction (MLT-LE, Ivanov & Polykovskiy 2022 arXiv:2209.06274; ResDTA, Rahman 2023 arXiv:2303.11434). **eMOSAIC** (Zhang 2025 *Nat Mach Intell*) uses Mahalanobis distance for OOD-aware residual confidence.

**Why now**: easiest possible meta-ranker; LambdaMART is one such residual-correction approach with rank-loss training, but a regression XGBoost on continuous pchembl is simpler and more interpretable. SHAP values on the trained model surface "MAMMAL is wrong when [feature]" insights.

**Effort**: 2-4 days end-to-end.

**Dependencies**: XGBoost (pip); existing parquets.

---

## 8. Creative / Out-of-the-Box Additions (gap-filling brainstorm)

Each subsection connects to a specific failure mode of the current pipeline. The first three (8.0a-c) are the *highest-impact* picks from the research streams; the rest (8.1+) are the original hand-thought brainstorm.

### 8.0a Restructure RRF as a Pareto front (NSGA-III axes)

**Gap**: RRF k=60 is a single-scalarizer — it hides Pareto-equivalent alternatives. A compound that wins by 0.01 RRF over another might be Pareto-dominated on safety / IP / route. The wet-lab shortlist should reveal the frontier, not collapse it.

**Concept** (research stream B): explicit non-dominated sorting via NSGA-III (Deb 2014, pymoo.org/algorithms/moo/nsga3.html) over **5 axes**:
1. Calibrated RRF score (the existing fused efficacy proxy)
2. ADMET violation count (B cluster gate stats)
3. Novelty: Tanimoto distance to nearest known-nootropic seed compound
4. Route / oral bioavailability (Caco-2 logPapp + B/B ratio prior)
5. IP freedom (from §8.3 patent cross-reference)

Output: top-50 as a **Pareto frontier**, not a ranked list. Surface "Pareto-optimal" tier (no compound dominates it across all 5 axes) and "near-Pareto" (within ε on any axis) separately.

**Anchors**:
- **DrugEx v3 / Pareto-MOEA** (Liu 2021 PMC8588612)
- **PMMG** (Adv Sci 2025 article 2410640): Pareto + MCTS over 7 objectives simultaneously, ~1-in-2 generated molecules satisfies all
- **Pareto MCTS for target-aware generation** (PNAS Nexus 2024 PMC11368924)
- **Honest landscape note**: even Nat Med 2024 TxGNN collapses to ranked list at the end. No one is doing Pareto well in repurposing today; we'd be among the first.

**Effort**: 5-7 days (pymoo NSGA-III wrap + tier renderer + integration into wet-lab shortlist v4).

**Dependencies**: pymoo (pip); novelty + IP axes need §8.3 + Bemis-Murcko computation.

### 8.0b PDSP-style 44-target off-target safety panel (extending ADMET)

**Gap**: ADMET-AI covers hERG, DILI, P-gp, CYPs. It does **not** cover:
- **5-HT2B** (valvulopathy — withdrew fenfluramine, pergolide, cabergoline)
- 5-HT2A halluc liability
- MAO-A / MAO-B irreversibility proxies
- α1-adrenergic (orthostatic hypotension)
- μ-opioid / κ-opioid (abuse potential, dysphoria)

**Concept**: extend the MAMMAL DTI grid to a **secondary 44-target liability sub-panel** (PDSP-inspired). Outputs become hard gates parallel to ADMET-AI. A compound that passes hERG/BBB/DILI but is a 5-HT2B agonist is auto-flagged.

**Threshold gate**: Dumotier 2024 (PMID 39032441) reaffirms regulatory ask = 10× Cmax/Ki margin.

**Why now**: this is the **single biggest risk-reduction we can buy for ~1 day of work**. The wet-lab shortlist currently has a hERG-shaped hole the size of a heart valve.

**Effort**: 1 day (extend `targets_seed.csv` with the 44-target liability sub-panel; re-run MAMMAL DTI; add gates in `gates/admet_gates.py`).

**Dependencies**: GPU for one ~1-hour MAMMAL re-run on the liability panel.

### 8.0c Cluster D = "neurobiological prior" via AHBA / Open Targets Genetics

See §7.9 for the full description. Listed here as a top-tier creative pick because it covers a *gap not addressed by anything else in the pipeline*: cognition-relevance of targets independent of binding data. Effectively a 5th cluster.

---

### 8.1 Selectivity-aware reranking → multi-class top-N

**Gap**: HRH3 single-target lock-in in the v3 shortlist (§4.4).

**Concept**: instead of one big top-25, surface multiple top-N tables per panel-class:
- Top-5 *mono-selective* compounds per panel target (selectivity score > threshold)
- Top-5 *polypharmacological* compounds (≥ 3 panel targets, balanced)
- Top-5 *novel-scaffold* compounds (Bemis-Murcko cluster not in known-nootropic seed)
- Top-5 by *Pareto frontier* across (efficacy, safety, novelty, regulatory bypass)

**Why now**: the methodology note already says ranking is rank-percentile, not Kd-magnitude — a faceted UI honours that limitation by surfacing different facets of "good."

**Effort**: 1-2 days (logic + render).

**Dependencies**: 7.4 selectivity scoring.

### 8.2 Combination-screening top-up (DrugComb)

**Gap**: known cognition-enhancing combinations (donepezil + memantine; caffeine + nicotine; modafinil + caffeine) are clinically validated but absent from the v3 ranking.

**Concept**: lookup pairs of our top-50 compounds in DrugComb (open-source DOI 10.1093/nar/gkab438) for any documented synergy (delta-Bliss > 5, S-score > 5). Surface combination candidates in a separate `wet_lab_combinations.md` report. For compounds without DrugComb data, predict synergy via target-overlap heuristic (compounds hitting orthogonal panel targets = candidate synergists).

**Why now**: combination dosing is how modafinil + caffeine actually beats single-agent stimulant; the pipeline ignoring it leaves 30-40% of the cognitive-enhancement literature on the floor.

**Effort**: 3-5 days (DrugComb integration + heuristic synergy predictor).

**Dependencies**: HTTP-only; no GPU.

### 8.3 Patent / clinical-trial cross-reference

**Gap**: a compound with weak IP / many active trials is a better repurposing target than one with strong IP / no trials. The pipeline doesn't surface this.

**Concept**: for each top-25, query Lens.org (open patents) + clinicaltrials.gov (public trials API) + USPTO Patent Examination Data (PEDS, free). Annotate:
- patent_status: {strong_IP, expired, off-patent, generic-available}
- active_trials: count + most-recent phase
- cognition_trials: count of those specifically with cognitive endpoints

Surface "easily repurposable" (off-patent + active trials) vs "hard repurposing" (strong IP + no trials) as a separate annotation column.

**Why now**: wet-lab handoff needs this distinction. A donepezil follow-on is hard; a piracetam follow-on is easy.

**Effort**: 2-3 days (HTTP + rate-limit handling + name → patent normalization).

**Dependencies**: HTTP-only.

### 8.4 Literature mining per top-25 (Semantic Scholar + Europe PMC)

**Gap**: the methodology note says "this is a prioritization, not a recommendation." But humans still need the prior literature context. The pipeline currently provides none.

**Concept**: for each top-25 compound × cognition keyword, query Semantic Scholar (free, no key) + Europe PMC for the top 3-5 most-cited papers in the last 5 years. Render a `wet_lab_shortlist_v3_with_lit.md` with a "prior reports" subsection per compound. Include: paper title, authors, year, abstract (1 sentence), conclusion (1 sentence).

**Why now**: human reviewers always Google their candidates anyway; pre-rendering the literature search saves 5-15 minutes per compound.

**Effort**: 2-3 days (HTTP + summarisation).

**Dependencies**: HTTP + optional LLM for abstract summarisation.

### 8.5 Off-target liability beyond ADMET — **superseded by §8.0b**

This subsection was the author's original 6-8 protein hand-thought version. **§8.0b expands it to a 44-target PDSP-style panel** with a regulatory-grounded threshold (10× Cmax/Ki margin per Dumotier 2024). Implementation overlap is total — do §8.0b, skip this.

### 8.6 Brain region selectivity via Allen Brain Atlas expression

**Gap**: a compound that hits HCRTR1 in cerebellum has different cognitive consequences than one that hits HCRTR1 in dorsolateral PFC. The pipeline treats targets as point-objects, not regionally-distributed.

**Concept**: for each panel target, pull Allen Brain Atlas adult human expression data (ISH or RNA-seq) for 6-8 cognition-relevant regions (dlPFC, HPC, BLA, ACC, parietal cortex, cerebellum, brainstem nuclei, basal ganglia). Compute a "regional selectivity vector." Rank compounds by *cognition-region selectivity*: targets that hit dlPFC + HPC strongly, cerebellum weakly, get a regional bonus.

**Why now**: orthogonal evidence stream that's free, public, and respects the "cognition is regional" reality.

**Effort**: 3-5 days (ABA API + regional aggregation + scoring).

**Dependencies**: HTTP-only.

### 8.7 Mechanism-of-action embedding

**Gap**: MAMMAL predicts pKd without distinguishing agonist / antagonist / PAM / NAM / partial / inverse. A DRD1 *agonist* and a DRD1 *antagonist* have opposite cognitive consequences; the pipeline doesn't notice.

**Concept**: add a per-(target, compound) MoA label from ChEMBL `action_type` field (when available). For compounds without ChEMBL action_type, predict from ChEMBL-trained classifier on Morgan FP + target ESM2 embedding. Surface MoA as a column in the wet-lab shortlist.

**Why now**: trivial to add given ChEMBL SQLite is live; high value for human review.

**Effort**: 1 day for the lookup pass; 3-4 days if the classifier needs training.

**Dependencies**: ChEMBL SQLite (✅ live).

### 8.8 Cross-DTI ensemble — **superseded by §7.7**

This subsection proposed DeepPurpose as the first additional ranker. **§7.7 now lists 9 candidate alternatives** with MMAtt-DTA as the top recommendation (reports ρ > 0.72 on transporter superfamily — exactly the family MAMMAL is bad on). Start there; DeepPurpose remains a fallback for a "cheapest possible" smoke test.

### 8.9 ANI-2x neural potential validation on top-N poses

**Gap**: Boltz predicts a binding affinity from a single inferred pose. Pose accuracy is uncertain. For the top-25 calibrated, we'd want a second-opinion energy on the predicted pose.

**Concept**: for each top-25 (compound, target) Boltz pose, run ANI-2x (open-source neural network potential, ~1 s per pose on RTX 5070) to compute a refined binding energy. Compare to Boltz's affinity; if they disagree by > 2 kcal/mol, flag.

**Why now**: cheap second opinion on Boltz's structural call. ANI-2x is parameter-light, blazing fast.

**Effort**: 3-5 days (ANI-2x integration + pose extraction + comparison report).

**Dependencies**: ANI-2x package (TorchANI), GPU.

### 8.10 Reverse-engineering known nootropics (similarity ranking)

**Gap**: the panel asks "what binds these 22 targets?" The complementary question is "what binds *like* the 30 known nootropics?"

**Concept**: for each known nootropic (donepezil, modafinil, methylphenidate, galantamine, ...) compute its 22-target binding fingerprint (MAMMAL pKd vector). For every compound in the library, compute cosine similarity to each known-nootropic fingerprint. Rank compounds by *max-similarity* across the known-nootropic set. Surface as `reports/nootropic_similarity_ranking.md` — a parallel ranking to the panel-based one.

**Why now**: the calibration is uncertain at individual targets; multi-target fingerprint similarity is more robust to per-target noise (it averages it out).

**Effort**: 2-3 days.

**Dependencies**: existing MAMMAL DTI grid.

### 8.11 Open-source LLM literature agent for the methodology note

**Gap**: the methodology note v1 (`reports/methodology_v1.md`) needs to be readable by someone who's never seen the project. A future-proof version of "what is this thing" is a working LLM that can answer questions against the methodology note + the reports + the code.

**Concept**: spin up a tiny RAG over `reports/`, `design/`, `configs/`. Local Qwen 7B or similar on CPU. Bind to a Slack bot or Claude-Desktop MCP. Anyone can ask "why is SLC6A3 down-weighted?" and get the right citation back.

**Why now**: the project has matured to the point where institutional knowledge is non-trivial to onboard cold. Either we write more docs or we automate the doc-querying.

**Effort**: 5-7 days (RAG + frontend + prompt engineering).

**Dependencies**: small open-source LLM; no GPU required.

### 8.12 Sprint-style cron-driven re-calibration

**Gap**: when new Boltz data lands, the v3 calibrated weights are silently stale. There's no automation to detect this.

**Concept**: a `scripts/_v3_watch_and_recalibrate.sh` that, when `boltzina_affinity.parquet` mtime is newer than `weights_calibrated.yaml` mtime, re-runs Phase A.7 → C → D → `26_v3_wet_lab_shortlist.py`. Schedule via cron / Windows Task Scheduler. Outputs a diff vs the prior calibration to `reports/calibration_drift_log.md`.

**Why now**: closes the operational loop. The pipeline becomes self-updating as Boltz coverage grows.

**Effort**: 1 day.

**Dependencies**: shell scripting.

---

## 9. Roadmap by Priority Tier

**Sequencing principle** (synthesised from research streams):
1. *Safety risk-reduction first* — never ship a wet-lab shortlist with a 5-HT2B liability you could have caught for 1 day of effort
2. *Diagnostic before treatment* — §7.1 tells us whether the v3 INVERTED finding is manifold mismatch (→ ensemble) or rank collapse (→ fine-tune)
3. *Cheap calibration second-opinion before expensive fine-tune* — isotonic (§7.11) + conformal gating (§7.12) + residual XGBoost (§7.14) are all < 1 week each; LoRA is 5-10 days minimum and may not improve
4. *Multi-class top-N before single-list polish* — diversify outputs structurally before optimising any single ranking
5. *Publication contributions are downstream* — the calibration linchpin paper writes itself once §7.6 / §7.7 produce comparison data

### Tier 1 — DO IMMEDIATELY (gated on overnight Boltz sweep completing, ~17h)

- [ ] **Re-run Phase A.7 → C → D → wet-lab v3** once `boltzina_affinity.parquet` has ≥3 predictions on most targets — the answer to "does Boltz rescue INVERTED targets?"
- [ ] **Phase B**: Cluster C orchestrator run (`23_v3_cluster_c.py` in WSL2 `txgnn_env`); produces `kg_scores.parquet`; re-fuse with 4 real rankers
- [ ] **PDSP-style off-target liability sub-panel** (§8.0b) — **biggest risk-reduction for 1 day**
- [ ] **Per-target failure-mode diagnostic** (§7.1) — tells us which v3+ track to invest in
- [ ] **Selectivity-aware reranking + multi-class top-N** (§8.1) — uses existing DTI grid; surfaces non-HRH3 candidates
- [ ] **Add report provenance hashes** — small follow-up across all report writers (§4.9)

### Tier 2 — DO SOON (1-2 weeks each, ordered by signal-to-effort)

- [ ] **Isotonic per-target calibration** (§7.11) — cleanest sign-correction; sklearn one-liner per target. 2-3 days.
- [ ] **Cross-DTI ensemble v0 with MMAtt-DTA** (§7.7) — *MMAtt-DTA reports ρ > 0.72 on transporter superfamily — directly addresses our SLC6A3/SLC6A2 failure*. 2-3 days.
- [ ] **Residual-correction XGBoost meta-ranker** (§7.14) — cheapest meta-learner; learns per-target offset from existing parquets. 2-4 days.
- [ ] **Mechanism-of-action embedding** from ChEMBL action_type (§8.7) — 1 day; high value for wet-lab review.
- [ ] **Patent / clinical-trial cross-reference** (§8.3) — gives Pareto axis 5 (IP freedom). 2-3 days.
- [ ] **Pocket-conditioned Boltz** (§7.5) — P2Rank + PocketMiner + CryptoBench. Biggest scientific credibility win. 3-5 days.
- [ ] **Cluster D neurobiological prior** (§7.9 / §8.0c) — Mansuri 2024 + AHBA + OT Genetics. 5-7 days.
- [ ] **Combination-screening top-up** via DrugComb (§8.2). 3-5 days.
- [ ] **Reverse-engineering known nootropics** — similarity ranking (§8.10). 2-3 days.

### Tier 3 — DO AFTER (2-4 weeks each)

- [ ] **Pareto front restructure of fusion** (§8.0a) — NSGA-III over 5 axes. 5-7 days but depends on §8.3 + §8.0b.
- [ ] **Conformal prediction per-target gating** (§7.12) — automatic per-target trust. 5-7 days.
- [ ] **Scaffold-aware active learning** for calibration expansion (§7.13). 3-5 days; gated on §7.1 confirming the tropane hypothesis.
- [ ] **Sign-correction stress-test mode** (§7.2) — secondary to isotonic (§7.11).
- [ ] **GWAS-anchored panel expansion** to 40-80 targets (§7.3) — 2-3 weeks.
- [ ] **Brain region selectivity** via Allen Brain Atlas (§8.6).
- [ ] **Literature mining for top-25** (§8.4) — Semantic Scholar + Europe PMC.
- [ ] **Selectivity scoring layer** (§7.4) — quantitative metrics (Uitdehaag, Cheng, KISS-CL).
- [ ] **Cron-driven auto-recalibration** (§8.12) — operational loop.
- [ ] **ANI-2x pose validation** on top-25 (§8.9).
- [ ] **LambdaMART promotion** — now eligible (≥20 labels per target on the well-covered targets).

### Tier 4 — RESEARCH (3-10 days, may require parallelism / cloud burst)

- [ ] **LoRA fine-tune MAMMAL** on cognition-DTI corpus — BALM-style LoKr/LoHa r=8 (§7.6). Only after §7.1 diagnostic + §7.7 ensemble baseline; tells us if LoRA is worth attempting.
- [ ] **Generative top-up with REINVENT 4 / POLYGON / PILOT** conditioned on calibrated panel (§7.8).
- [ ] **Boltzina-Vina-only mode** (V3 attack plan T5) — ~11× speedup; never implemented.
- [ ] **LLM-agent prioritisation overlay** (MCP server) (§7.10) — AgentD / PharmaSwarm template.
- [ ] **LLM literature agent for the methodology note** (§8.11).

### Tier 5 — PUBLISHABLE (multi-week, methodology contributions)

- [ ] **The calibration-linchpin methodology paper** — Phase A.7's per-target Spearman ρ approach + the refined verdict matrix is itself a contribution. The pipeline that produces it is just the demonstration vehicle.
- [ ] **MAMMAL transporter-failure-mode paper** — if §7.1 confirms tropane-saturation, this is a publishable negative result with mitigation (ensemble with MMAtt-DTA, scaffold-aware AL).
- [ ] **Allosteric awareness benchmark v2** — expand the Phase 4.1 work with Boltz-rescued PAMs; n ≥ 10 per (target, binding mode) × 6+ targets.
- [ ] **Cognition virtual phenotype anchor validation paper** — TxGNN's anchor is a methodological contribution; benchmark formally against curated indication ground truth.
- [ ] **Pareto-vs-RRF benchmark in CNS repurposing** — if §8.0a delivers, we'd be among the first repurposing pipelines actually using Pareto output (Nat Med 2024 TxGNN collapses to a ranked list).

---

## 10. Risk Register (post-sprint)

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| Overnight Boltz sweep doesn't rescue INVERTED targets (most flip to DE_WEIGHT_TARGET, not BOLTZ_2X_MAMMAL) | Medium | High | The 4 transporter/NMDA targets become "no useful prediction at all"; pipeline must surface this honestly, deprecate those targets, and lean on the remaining 18. Methodology note already permits this. |
| Cross-DTI ensemble shows that ALL DTI methods are inverted at DAT/NET → the problem is the data, not MAMMAL | Medium | High | This is the most publishable negative result: foundation models trained on BindingDB share a systematic bias at transporters. Document, publish, move to LoRA fine-tune on a non-BindingDB corpus. |
| LoRA fine-tune doesn't improve per-target ρ above +0.30 even with cognition-specific corpus | Medium | Medium | Falls back on cross-DTI ensemble (multiple bad ranker × averaging > single bad ranker). |
| Cluster C (PrimeKG + TxGNN) produces low-signal scores once run (zero-shot indication for "cognition" returns noise) | Medium | Medium | Anchor the panel of 5 diseases is already weighted; consider re-weighting or replacing with a 10-disease panel including dementia subtypes, MS-cognition, post-stroke cognition. |
| GWAS-anchored panel expansion surfaces 100+ new targets, panel becomes unmanageable | Low | Low | Strict tractability filter: only proteins with known small-molecule chemistry (ChEMBL ≥ 50 records). Caps expansion at ~30-40 new targets. |
| Generative chemistry (REINVENT) produces synthesizable but commercially-blocked molecules | Medium | Low | The IP cross-reference layer (§8.3) catches this; surface as a separate "needs medicinal chemistry consult" bucket. |
| HRH3 single-target lock-in survives even multi-class top-N (most predictions still HRH3-best) | Medium | Low | Cap per-mechanism-class top-N to 5 each; explicitly diversify. |
| Roberts CA 2020 effect-size ceiling (SMD ≈ 0.21) makes any wet-lab spend hard to justify, sprint-after-sprint | High (structural) | High | Reframe the deliverable: the methodology contribution (calibration + provenance + 4-cluster fusion) is the publishable artifact; the candidate ranking is a downstream demonstration vehicle. |
| WSL2 Boltz sweep fails / corrupts the parquet | Low | Medium | Sweep is append-mode with per-row error handling; bad rows logged as NaN, recoverable. |

---

## 11. Single-Page Cheat Sheet — Commands That Work Today

```powershell
# Windows env (mammal_env, conda)
$envPython = "$env:USERPROFILE\.conda\envs\mammal_env\python.exe"
$repo = "C:\Users\Pierce Lonergan\Documents\GitHub\MAMMAL_Cognitive_Enhancement_Drug_Repurposing"
Push-Location $repo

# v3 SQLite-backed Phase A (linchpin + audit + backstop):
& $envPython scripts\_v3_sqlite_vs_rest_smoke.py         # A.5 gate
& $envPython scripts\24_v3_audit_chembl_targets_sqlite.py # A.6 clean audit
& $envPython scripts\22_v3_calibration.py                 # A.7 linchpin
& $envPython scripts\21_v3_chembl_evidence_sqlite.py --all-pairs  # A.4 backstop (99 min)

# v3 Phase C — both fusion passes:
& $envPython scripts\15_v2_fusion.py --calibrated-weights NONE --out-suffix _uncalibrated
& $envPython scripts\15_v2_fusion.py --out-suffix _calibrated

# v3 Phase D — calibrated vs uncalibrated diff:
& $envPython scripts\25_v3_fusion_diff.py

# v3 wet-lab shortlist with 4-cluster scorecards:
& $envPython scripts\26_v3_wet_lab_shortlist.py
```

```bash
# WSL2 — overnight Boltz sweep status:
wsl -d Ubuntu -- bash -c "tail -c 800 /tmp/wsl2_boltz_sweep.log | tr '\r' '\n' | tail -5"

# WSL2 — when overnight sweep done, Cluster C orchestrator:
wsl -d Ubuntu -u root -- bash -c "source /root/txgnn_env/bin/activate && python /mnt/c/Users/Pierce\ Lonergan/Documents/GitHub/MAMMAL_Cognitive_Enhancement_Drug_Repurposing/scripts/23_v3_cluster_c.py"
```

---

## 12. The Mission Statement (still true, sharpened by the calibration linchpin)

This pipeline is an honest attempt to use a foundation model (MAMMAL) plus a curated panel of cognition-relevant targets plus rigorous filtering (ADMET + structural + mechanistic) plus *ground-truth-calibrated weighting* to surface drug-repurposing candidates with mechanistic plausibility and physical feasibility for healthy adult cognitive enhancement.

The Roberts 2020 meta-analysis ceiling is real: even methylphenidate hits SMD = 0.21 on a generous metric. We are not searching for a miracle drug; we are enriching a candidate set so wet-lab cycles spend money on plausibility, not chemistry-lottery tickets.

The v2 already did this better than v1 — top 10 went from anonymous ChEMBL artifacts to actual CNS drugs. The v3 sprint added: ChEMBL SQLite (78h → 99 min), per-target Spearman ρ calibration, refined verdict matrix that honestly handles negative correlation, 4-cluster RRF with per-target weights, the methodology note, the wet-lab shortlist v3 with provenance.

What v3 exposed is uncomfortable and important: MAMMAL is anti-correlated at the most clinically validated cognition targets. The v4 work in this roadmap is the response: diagnose, ensemble, fine-tune, expand the panel, surface selectivity, add behavioural priors, generate novel candidates, and document everything so the next agent picking this up at commit `a2ff155` can pick up the right thread.

Life would be a lot better if people had better cognition. Let's keep getting on with it.

---

## Appendix A — Research citation index (consolidated)

Citations referenced inline in §7 and §8, grouped by stream.

### A.1 Foundation-model fine-tuning / sign-correction (research stream A)

**LoRA / PEFT for DTI**
- **BALM** — Gorantla et al. 2024, bioRxiv 2024.11.01.621495; github.com/meyresearch/BALM (Apache-2.0) — LoKr + LoHa PEFT on ESM-2 + ChemBERTa-2; beats full fine-tune on every BindingDB split
- Schmirler et al. 2024, *PNAS* 121:e2405840121 — "Democratizing protein language models with PEFT" — ranks ≥4 required; diminishing returns above r=8
- ESM-LoRA-Gly — Yan et al. 2025, bioRxiv 2025.08.12.669850
- SeqProFT — arXiv:2411.11530, 2024
- MAMMAL LoRA reproduction (BBBP) — github.com/Kymi808/mammal-lora-bbbp
- MAMMAL repo — github.com/BiomedSciAI/biomed-multi-alignment; Shoshan et al. 2025 *npj Drug Discovery* article s44386-026-00047-4

**Data prep + curation**
- chembl_structure_pipeline — Bento et al. 2020, *J Cheminform* 12:51
- BASE webservice — Song et al. 2024, *BMC Bioinformatics* 25 article 305

**Catastrophic forgetting**
- Kirkpatrick et al. 2017 — EWC
- Aljanabi & Khaleel 2025 — arXiv:2501.13669 — hierarchical regularization with replay
- Qiao et al. 2024 NeurIPS — "Learn more, but bother less" — adapter routing per class
- Brenndoerfer 2024 survey — mbrenndoerfer.com

**Few-shot floors**
- Stanley et al. 2021 — FS-Mol
- Lopez et al. 2024 — *J Chem Inf Model* 64 — few-shot ligand-based compound activity

**Calibration**
- Toplak et al. 2020 — *J Chem Inf Model* 60:3829 (PMID 32865408) — Venn-ABERS wins Brier on 40M pairs × 2,112 targets
- Kull et al. 2017 — beta-calibration parametric
- Neelon & Dunson 2004 — Bayesian isotonic prior

**Conformal prediction for DTI**
- Rakhshaninejad et al. 2025 — arXiv:2505.18890 — cluster-conditioned CP
- Bosc et al. 2023 — *J Cheminform* 15:79 (PMC10457664) — dynamic applicability domain

**MoE / gating**
- EnsDTI — Park et al. 2024, bioRxiv 2024.08.06.606753
- MoSE — Yan et al. 2025, arXiv:2503.15796 (AAAI 2025)

**Per-compound OOD**
- eMOSAIC — Zhang et al. 2025, *Nature Mach Intell* article s42256-025-01151-2 (bioRxiv 2024.01.05.574359)

**Alternative DTI heads**
- PSICHIC — Koh et al. 2024, *Nat Mach Intell* 6:673; github.com/huankoh/PSICHIC (MIT)
- ConPLex — Singh et al. 2023, *PNAS* 120:e2220778120; github.com/samsledje/ConPLex (MIT)
- Komet — Guichaoua et al. 2024, *J Chem Inf Model* 64:6938 (BSD-3); Zenodo 10731712 for LCIdb
- PLAPT — Rose et al. 2024, bioRxiv 2024.02.08.575577 (MIT)
- GEMS — Lemm et al. 2024, bioRxiv 2024.12.09.627482 (Apache-2.0)
- Boltz-2 — Passaro et al. 2025, bioRxiv 2025.06.14.659707; github.com/jwohlwend/boltz (MIT)
- DrugCLIP — Gao et al. 2024, arXiv:2310.06367 (MIT)
- **MMAtt-DTA** ★ — Schulman et al. 2024, *Bioinformatics* 40:btae496; github.com/AronSchulman/MMAtt-DTA (MIT) — **ρ > 0.72 on transporter superfamily**
- MGNDTI — Peng et al. 2024, *J Chem Inf Model* (PMID 39137398) (MIT)

**Active learning**
- Bailey et al. 2024 — *eLife* 12:RP89679 (bioRxiv 2023.07.26.550653) — deep batch AL for drug discovery
- Wang & Pyzer-Knapp 2024 — *J Chem Inf Model* 64 doi 10.1021/acs.jcim.4c00220
- Holzmüller `bmdal_reg` — github.com/dholzmueller/bmdal_reg
- Graff et al. 2021 — *Chem Sci* 12:7866; Fromer & Coley 2023 — arXiv:2310.10598
- Bioactivity Similarity Index — Skinnider 2025, *Front Bioinform* 5 doi 10.3389/fbinf.2025.1695353

**Chemistry context** (the tropane hypothesis)
- Carroll et al. 2004, *J Med Chem* / PMID 15566309 + PMID 14711303 — RTI-series cocaine analogs

**Residual / meta-ranker patterns**
- MLT-LE — Ivanov & Polykovskiy 2022, arXiv:2209.06274
- ResDTA — Rahman et al. 2023, arXiv:2303.11434

### A.2 Creative additions / 2026 CNS prioritization (research stream B)

**Selectivity scoring**
- Uitdehaag & Zaman 2011 — *PLOS ONE* PMC3100252 — selectivity entropy
- Cheng et al. 2010 — *J Med Chem* doi 10.1021/jm100301x — Partition Index
- Graczyk-style Gini + S(10x) — *BMC Bioinformatics* 2017 doi 10.1186/s12859-016-1413-y
- KISS-CL — *Nat Commun* 2025 article 65869-8 — contrastive selectivity learning

**Multi-objective Pareto**
- DrugEx v3 / Pareto-MOEA — Liu 2021 PMC8588612
- PMMG — *Adv Sci* 2025 article 2410640 — Pareto + MCTS, 7 objectives
- pymoo NSGA-III — pymoo.org/algorithms/moo/nsga3.html
- Pareto MCTS — PNAS Nexus 2024 PMC11368924
- TxGNN — *Nat Med* 2024 article s41591-024-03233-x (the field's current SOTA but still collapses to a list)

**Generative chemistry**
- REINVENT 4 — Loeffler et al. 2024, *J Cheminform* PMC10882833
- POLYGON — *Nat Commun* 2024 article 47120-y — dual-target ligand RL
- PILOT — *Chem Sci* 2024 article d4sc03523b — equivariant diffusion + pocket conditioning
- Pocket2Mol-RL — Deargen 2024, github.com/deargen/Pocket2Mol_RL_public

**Pocket detection / conditioning**
- P2Rank — github.com/rdk/p2rank
- PocketMiner — Meller 2023 *Nat Commun* article 36699-3
- CryptoBench — Pólya 2024 *Bioinformatics* 41 article btae745 — ESM2-based, beats P2Rank + PocketMiner
- AlphaFold3 / cofolding for cryptic discovery — Bryant 2023 *Nat Commun* PMC10373493

**Neurobiological priors**
- Mansuri et al. 2024 — PMC10941541 — 41 cortical genes vs intelligence-related fMRI
- eLife reviewed preprint 2025 — UK Biobank multimodal MRI marker of cognition (48% g-factor variance) — elifesciences.org/reviewed-preprints/108109
- Open Targets Genetics — genetics.opentargets.org
- Lit-OTAR — Sangkuhl 2024 PMC11978389 — 48M target-disease-drug associations
- Hansen 2020 AHBA + Neurosynth pipeline — bioRxiv 2020.07.16.203026

**Off-target safety**
- Dumotier 2024 — PMID 39032441 — 10× Cmax/Ki margin regulatory ask

**MoA classification**
- AiGPro 2025 — PMC11780767 — multi-task agonist/antagonist heads, AUC > 0.85 on GPCRs

**Pose validation**
- ANI-2x for virtual screening — PMC11201553 — geometry-opt + neural-potential rescoring

**Cross-DTI ensemble**
- DeepPurpose — Huang K et al. 2020, *Bioinformatics* 36:5545

**LLM-agent orchestration**
- AgentD — *J Chem Inf Model* 2025 doi 10.1021/acs.jcim.5c02454
- PharmaSwarm — arXiv:2504.17967

---

## Appendix B — How this doc was written

The deep-dive that produced §7 and §8 was assisted by two parallel sub-agent research streams (Anthropic Claude general-purpose), run during the V3 sprint:

1. **Stream A**: "Foundation-model fine-tuning + sign-correction" — focused on the MAMMAL_ONLY_INVERTED finding. Returned ~3,000 words, ~40 citations, identifying BALM-pattern LoKr/LoHa LoRA at r=8 + scaffold-aware AL + MMAtt-DTA ensemble + isotonic per-target calibration as the highest-leverage paths.
2. **Stream B**: "Creative additions / 2026 CNS prioritization" — focused on gaps in the current pipeline. Returned ~1,800 words across 6 streams (selectivity, Pareto MOO, generative, pocket conditioning, neurobiological priors, plus a creative brainstorm), with concrete platforms + wall-clock estimates.

Both streams' raw outputs are preserved in the agent transcripts. The synthesis above (priority tiers, integration into existing v3 architecture) is the author's plus this assistant. Anywhere ★ appears in §7 / §8, that's a top-recommendation from one of the research streams.

When v3+ work begins, the next assistant should re-spawn these streams (the literature moves fast — 2026 has more SOTA than 2024). Useful prompts archived in the conversation transcript at commit `a2ff155`.

---

