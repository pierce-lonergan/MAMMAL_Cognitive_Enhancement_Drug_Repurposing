# V2 Hybrid Architecture — Status, Compatibility Issues, Forward Plan

**Read this first.** This is the source-of-truth status document for the project. The architecture spec lives in [V2_HYBRID_ARCHITECTURE.md](V2_HYBRID_ARCHITECTURE.md); this doc explains *what's actually working*, *what's broken*, *what we need to fix*, and *which library-compatibility battles are worth fighting*.

---

## 1. Executive Summary

The MAMMAL Cognitive Enhancement Drug Repurposing pipeline is a 3-cluster hybrid:
- **Cluster A** (Structure): MAMMAL DTI + ESM2-650M + Boltz-2/Boltzina
- **Cluster B** (ADMET): ADMET-AI (41 endpoints + hard gates)
- **Cluster C** (Mechanism/KG): PrimeKG + TxGNN + cognition virtual anchor

Joined by Reciprocal Rank Fusion (RRF), with hard ADMET gates applied before fusion. Final outputs: ranked candidate list, provenance, disagreement diagnosis, prose narrative per candidate.

**As of now (commit `ebe938b`)**:

| Cluster | Component | Status | Empirical signal |
|---|---|---|---|
| A | MAMMAL DTI head | ✅ LIVE | 6,556 pairs scored; positive controls partly pass (4/7) |
| A | ESM2-650M embeddings | ✅ LIVE | 22 targets cached; paralog validation PASS (cos 0.997, gap 0.097) |
| A | Boltz-2 structure + affinity | ✅ LIVE (slow) | Smoke (CHRNA7, galantamine) = 149.8s/pair; cuequivariance kernel unavailable on Windows |
| A | Boltzina (Vina-pose-only) | 📐 STUBBED | No published PyPI package; would need vendored impl from Furui & Ohue 2025 |
| B | ADMET-AI 41 endpoints | ✅ LIVE | 298 compounds scored, 53 PASS / 64 FLAG / 181 CUT; positive ctrls pass BBB via regulatory bypass |
| B | Hard gates | ✅ LIVE | configs/thresholds.yaml; regulatory bypass for approved drugs |
| C | PrimeKG loader | 📐 STUBBED | Code ready, multi-GB Harvard Dataverse download not done |
| C | TxGNN scorer | 📐 STUBBED | Code ready, install blocked by PyG-on-Windows |
| C | Cognition virtual anchor | ✅ LIVE | 5-disease union (MCI+AD+ADHD+FXS+narcolepsy) loaded from configs/weights.yaml |
| Fusion | RRF | ✅ LIVE | 2-ranker fusion (MAMMAL + ADMET) working; absorbs Boltzina+TxGNN additively |
| Fusion | LambdaMART | 📐 STUBBED | Gated on ≥20 labels (we have 0) |
| Provenance | Tracker + disagreement + narrative | ✅ LIVE | Top-20 funnel narratives generated; disagreement archetypes will populate when Boltz+TxGNN land |
| Calibration | Platt vs ChEMBL ground truth | 🔄 BLOCKED | Waiting on ChEMBL v2 backstop (still grinding; ~hours remaining) |

**The empirical jump v1 → v2 is real**: top 10 went from anonymous ChEMBL structural-bias artifacts (staurosporine, fenpropimorph, chembl176261 etc.) to actual CNS drugs (modafinil, piracetam, bupropion, levetiracetam, rasagiline, rolipram, BI-409306, aniracetam, lanicemine, d-amphetamine) — 3 explicit positive controls in top 10, zero `chembl_binder` noise. This is the value of the hybrid + the ADMET gate + regulatory bypass.

**The single most important pending result**: the focused Phase 0.5 Boltz sweep (CHRNA7 + PDE4D + HRH3, ~12 pairs, in flight) will tell us whether Boltz-2 rescues the α7 nAChR PAMs (galantamine, encenicline, TC-5619) into the top 25% at CHRNA7. That single empirical answer decides whether v3 keeps investing in the structural cluster.

---

## 2. What Works Today (end-to-end runnable)

```
                Windows (mammal_env)                          Status
v1 pipeline: 02_fetch_targets ... 13_wet_lab_shortlist     ✅ committed, runnable
v2 ADMET cluster: 14_v2_cluster_b_admet → admet_gates.parquet ✅ runs in ~1 min
v2 fusion:   15_v2_fusion → rrf_ranking + provenance + narrative ✅ runs in seconds
v2 ESM2:     16_v2_esm2_embed → 22 cached .pt files        ✅ runs in ~5 min
v2 Boltz:    17_v2_boltz_smoke → single (target, compound)  ✅ runs in ~150s
             18_v2_boltzina_sweep → bulk parquet            ✅ runs (slow)
             19_v2_allosteric_gate → reports/...md          ✅ runs in seconds (depends on 18 output)
             _boltzina_focused (minimal CHRNA7+PDE4D+HRH3)  🔄 running now, ~30 min
```

**Reproducibility**: every artifact lands in `data/results/v2/`. Caches at `data/cache/{esm2,boltz_struct,boltzina,admet,txgnn}/`. Config at `configs/{thresholds,weights}.yaml`. Module API at `src/mammal_repurposing/{cluster_a,cluster_b,cluster_c,gates,fusion,provenance,pipeline}/`.

---

## 3. What's Broken (and the deep research to unblock it)

### 3.1. `cuequivariance-ops-torch-cu12` — Linux-only wheel

**Symptom**: Boltz-2's triangle attention layers try `from cuequivariance_torch.primitives.triangle import triangle_multiplicative_update`, which at runtime imports `cuequivariance_ops_torch` (the native CUDA kernel). PyPI publishes wheels for `manylinux_2_27_x86_64`, `manylinux_2_28_x86_64`, `manylinux_2_26_aarch64`, `manylinux_2_28_aarch64` — **no Windows wheel**.

**Workaround in place**: Pass `--no_kernels` to `boltz predict`. Boltz falls back to a pure-PyTorch triangle multiplicative update. Slower per pair (~150s instead of estimated ~30-50s on the same hardware with kernels).

**Deep-research findings**:
- The wheel is published by NVIDIA at the cuEquivariance repo (`NVIDIA/cuEquivariance`) v0.10.0 (released 2026-04-22), Python 3.10-3.14.
- **NVIDIA's build system does not target Windows.** No issues or PRs requesting Windows support as of investigation; this is an upstream design choice (CUDA + ninja + custom C++ build for Windows MSVC is painful).
- Building from source on Windows would require: CUDA Toolkit 12.x for Windows, Visual Studio 2022 Build Tools, ninja, the cuEquivariance source tree, ~hours of compilation and debugging.

**Recommended remediation**: **run Boltz inside WSL2 Ubuntu** with the Linux wheel. The host's Windows env stays for MAMMAL/ADMET/fetchers. See §6 for the migration plan. Estimated speedup: **2.5-5×** per Boltz call (research doc claimed 11.8× via kernels; conservative estimate of 2.5-5× accounts for the WSL2 overhead).

**Severity**: HIGH. This is the difference between a 60-hour Phase 0.4 sweep and a 6-18 hour sweep, i.e. between "infeasible in one session" and "overnight background run."

### 3.2. `tiledbsoma` — Windows wheel doesn't build

**Symptom**: `pip install tiledbsoma` fails at wheel compilation on Windows (TileDB-SOMA's native C++ extension needs CMake + a specific TileDB version that's not pre-built for Windows). PyTDC eagerly imports `tdc.multi_pred.perturboutcome` → `single_cell.py` → `cellxgene_census` → `tiledbsoma`, breaking the entire `mammal.examples.dti_bindingdb_kd.task` import chain.

**Workaround in place**: `src/mammal_repurposing/_compat.py` installs `sys.modules` stubs for `tiledbsoma`, `cellxgene_census`, `gget`, `tdc.multi_pred.perturboutcome`, `tdc.multi_pred.single_cell`, `tdc.resource.cellxgene_census`. The stubs raise loudly on actual attribute access but satisfy the import-time chain.

**Deep-research findings**:
- TileDB Inc. publishes Windows wheels for `tiledb` (the storage engine) but not for `tiledbsoma` (the SOMA spec impl).
- The tiledbsoma GitHub issue queue has ongoing Windows-support discussions; no shipped wheel as of investigation.
- conda-forge does not have a Windows tiledbsoma either (verified Phase 0.1).
- We'd need this for **Cluster D** (transcriptomic / LINCS L1000 via cellxgene-census), which the research doc deliberately defers.

**Recommended remediation**: **defer.** Tiledbsoma is only needed for Cluster D, which the v2 research doc explicitly says NOT to build until RRF is converged and a neuron-relevant L1000 dataset materializes. The stub is functional; revisit only if Cluster D becomes a priority. If WSL2 migration happens (per §6), tiledbsoma installs cleanly on Linux, unblocking Cluster D for free.

**Severity**: LOW (current scope), MEDIUM (if Cluster D enters scope).

### 3.3. `torch-geometric` (PyG) on Windows — installation gauntlet

**Symptom**: TxGNN (Cluster C) requires PyTorch Geometric. PyG's binary dependencies (`torch-scatter`, `torch-sparse`, `torch-cluster`, `torch-spline-conv`) need pre-built wheels matched to (PyTorch version, CUDA version, Python version, OS). On Windows + PyTorch nightly cu128 + Python 3.10, no matching pre-built wheels exist. Building from source needs MSVC + CUDA toolkit + ~30 min compile per package, often fails on Blackwell sm_120.

**Workaround in place**: Cluster C is stubbed at the runner level (`cluster_c/txgnn.py` raises `NotImplementedError`). Cognition anchor (cluster_c/cognition_anchor.py) is live as a dataclass.

**Deep-research findings**:
- PyG team publishes wheels at `https://data.pyg.org/whl/torch-<version>+cu<version>.html`. Latest stable: PyTorch 2.5.1, CUDA 12.1. Nothing for nightly 2.12 cu128.
- Conda-forge has `pytorch-geometric` for Linux + CUDA but no Windows-Blackwell entries.
- TxGNN's setup.py pins specific PyG versions; updating to current PyG often breaks model checkpoint loading.

**Recommended remediation**: **run TxGNN inside WSL2 Ubuntu** with the linux+cu128 PyG wheels. Same pattern as Boltz (§3.1). The cognition anchor + KG_score formulas stay portable Python (CPU-only), so the orchestration can live on Windows and call the WSL2 TxGNN inference like an RPC.

**Severity**: HIGH (Cluster C is non-trivial to enable on Windows native).

### 3.4. ChEMBL REST is the throughput bottleneck for v1.2 evidence backstop

**Symptom**: `scripts/07_chembl_evidence.py` at threshold pKd > 6.0 has 4,713 pairs to lookup. Each pair does 1 InChIKey→molecule_chembl_id lookup + 1 activity query = 2 HTTP calls per pair. ChEMBL intermittently returns 500s, which my tenacity decorator retries 4× with exponential backoff (up to ~30s per failing call). Effective throughput observed: **~1 pair/min** wall-clock. Full sweep ETA: ~78 hours. Currently running at threshold 6.5 (875 pairs, ~14 hours) in background.

**Deep-research findings**:
- ChEMBL has a bulk download (`chembl_<release>_sqlite.tar.gz`, ~4 GB) that gives instant queries against a local SQLite database. No rate limit, no 5xx flakiness.
- ChEMBL also has a SQL-over-PostgreSQL via their public mirror — works but requires Postgres client install.
- The InChIKey→molecule_chembl_id step is the slowest because for many ChEMBL-expanded compounds we already KNOW the molecule_chembl_id (it's stored in `alt_names` field of our compounds parquet from the `top_binders_for_targets` call).

**Recommended remediation (in priority order)**:
1. **Pre-extract molecule_chembl_id from compounds.parquet alt_names** — for compounds in our ChEMBL-expansion subset (198 of 298), we already have the ID; no lookup needed. Cuts ~2,000 HTTP calls.
2. **Build a local ChEMBL SQLite mirror** — one-time 4 GB download, then queries take ms not seconds. Total Phase 1.2 time drops from ~14 hours to ~5 minutes.
3. Failing both: cap the threshold to 6.7 (101 pairs, ~100 min), accept lower coverage.

**Severity**: MEDIUM. ChEMBL evidence is needed for Phase 3.1 calibration; without it, calibration is "insufficient data."

### 3.5. Boltz on Windows is slow without cuequivariance kernels

**Symptom**: 150s per (target, compound) pair on Blackwell-Windows. Full 1,500-pair Phase 0.4 sweep = ~62 hours.

**Workaround in place**: Focused sweep at 3 critical targets (CHRNA7, PDE4D, HRH3) × allosteric panel ligands only = ~12 pairs × 150s = ~30 min. Sufficient for the Phase 0.5 gate test but not for the full v2 wet-lab shortlist.

**Recommended remediation**: WSL2 migration of just the Boltz call (§6). Conservatively 2-5× speedup → 12-30 hour sweep (overnight). With cuequivariance kernels fully working (estimate ~30-50s per pair) → 8-15 hour sweep.

**Severity**: HIGH for full Phase 0.4 throughput.

### 3.6. Conda-on-Windows path-with-spaces bug

**Symptom**: `conda.exe` launcher (Scripts/conda.exe) crashes when user profile path contains spaces ("Pierce Lonergan") because the launcher constructs a command line that can't be parsed by `CreateProcess`.

**Workaround in place**: Use `_conda.exe` (the native binary, not the launcher) for all `conda` invocations. Documented in `scripts/01_setup_env.ps1` and used by all env-management code.

**Severity**: LOW (workaround is stable).

### 3.7. HuggingFace symlink-creation needs Administrator on Windows

**Symptom**: `Mammal.from_pretrained(<HF_ID>)` triggers HF cache snapshot creation, which uses symlinks (`snapshots/<sha>/file.txt → ../../blobs/<sha>`). Windows symlink creation requires `SeCreateSymbolicLinkPrivilege`, only granted to Administrator or accounts in Developer Mode. Without it, the cache populate fails partway.

**Workaround in place**: For Cluster B (`cluster_b/admet_ai_runner.py` uses HF, no issue), Cluster A MAMMAL DTI (HF, works because pre-populated cache existed), and MoleculeNet heads we use `snapshot_download(local_dir=..., local_dir_use_symlinks=False)` to flatten downloads into a non-symlinking directory. Documented in `cluster_b/admet_ai_runner.py` and the molnet sub-runners.

**Severity**: LOW once worked around. Would be cleaner with Developer Mode enabled (user setting, one click).

### 3.8. `protobuf` version conflicts between admet-ai, boltz, tensorflow

**Symptom**: admet-ai pulls a torch CPU wheel that triggers a protobuf downgrade; boltz pulls protobuf 5.x; tensorflow (in transformers pipeline) needs protobuf 6.31+. Result: silent breakage of `from transformers import AutoModel`.

**Workaround in place**: Hard-pin `protobuf>=6.31` after every install. Documented in scripts/01_setup_env.ps1 (TODO: add to that script as a post-install step).

**Severity**: LOW. Easy to fix once you know about it; would burn 30+ min if discovered fresh.

---

## 4. Throughput Bottlenecks (real numbers from current state)

| Operation | Per-call wall-clock | Total grid | ETA |
|---|---|---|---|
| MAMMAL DTI inference | 100 ms (batch 8) | 6,556 pairs | ~10 min |
| ESM2-650M embedding | 5 s (single target) | 22 targets | ~5 min |
| ADMET-AI 41 endpoints | 100 ms (CPU) | 298 compounds | ~30 s |
| Boltz-2 affinity (no_kernels) | 150 s | 1,500 pairs | ~62 h |
| Boltz-2 affinity (with kernels, projected from research doc) | 30-50 s | 1,500 pairs | ~12-20 h |
| ChEMBL evidence (current impl) | 60 s | 4,713 pairs | ~78 h |
| ChEMBL evidence (with local SQLite mirror) | 10 ms | 4,713 pairs | ~1 min |
| OpenTargets GraphQL | 200 ms | 22 targets | ~5 s |
| RRF fusion + provenance | 200 ms | one call over full grid | ~1 s |

**The two slow operations dominating the project schedule are**: Boltz (because of the missing Linux-only kernel) and ChEMBL (because of remote API latency + retries). Fixing #3.1 (WSL2) and #3.4 (local SQLite mirror) compresses the project schedule from ~weeks to ~days.

---

## 5. Component Status Matrix (one-pass reference)

### 5.1 Library installs

| Library | Windows native | WSL2 Ubuntu | Notes |
|---|---|---|---|
| torch 2.12 cu128 | ✅ (nightly) | ✅ (nightly) | Blackwell sm_120 needs nightly until stable lands |
| biomed-multi-alignment (mammal) | ✅ | ✅ | Both work; `[examples]` extra requires PyTDC stub on Windows |
| admet-ai | ✅ | ✅ | Beware: drags CPU torch as a transitive dep; pin nightly cu128 after |
| boltz | ✅ (slow path) | ✅ (fast path) | Without cuequivariance kernels on Windows |
| cuequivariance-torch | ✅ | ✅ | Python frontend; installs everywhere |
| cuequivariance-ops-torch-cu12 | ❌ Linux only | ✅ | The native CUDA kernel; this is the asymmetry |
| transformers (ESM2 host) | ✅ | ✅ | Both work |
| tiledbsoma | ❌ no wheel | ✅ | Stubbed; needed only for cellxgene-census (Cluster D, deferred) |
| cellxgene-census | ❌ depends on tiledbsoma | ✅ | Stubbed |
| torch-geometric (PyG) | ❌ no pre-built for nightly | ✅ | Needed for TxGNN (Cluster C) |
| txgnn | ❌ via PyG | ✅ | Stubbed |
| networkx | ✅ | ✅ | For PrimeKG loader |
| lightgbm | ✅ | ✅ | For LambdaMART (not needed until ≥20 labels) |

### 5.2 V2 modules

| Module | LOC | Status | Test coverage |
|---|---|---|---|
| `cluster_a/esm2_embed.py` | ~80 | LIVE | Manual via scripts/16 |
| `cluster_a/boltz_runner.py` | ~150 | LIVE | Via scripts/17 smoke |
| `cluster_a/boltzina.py` | ~210 | LIVE (full mode), STUB (Vina mode) | Via scripts/17 smoke |
| `cluster_b/admet_ai_runner.py` | ~120 | LIVE | Manual via scripts/14 |
| `cluster_c/cognition_anchor.py` | ~50 | LIVE (data only) | None |
| `cluster_c/primekg.py` | ~70 | STUB (load func raises if KG missing) | None |
| `cluster_c/txgnn.py` | ~30 | STUB (NotImplementedError) | None |
| `gates/admet_gates.py` | ~180 | LIVE | Manual via scripts/14 |
| `fusion/rrf.py` | ~140 | LIVE | Used by scripts/15 |
| `fusion/lambdamart.py` | ~90 | STUB (real fit/predict, no test) | None |
| `fusion/calibration.py` | ~70 | LIVE (small tested impl) | None |
| `provenance/tracker.py` | ~80 | LIVE | Manual via scripts/15 |
| `provenance/disagreement_report.py` | ~120 | LIVE | Manual via scripts/15 |
| `provenance/narrative.py` | ~110 | LIVE | Manual via scripts/15 |
| `_compat.py` (Windows stubs) | ~60 | LIVE | Implicit (mammal imports succeed) |
| `analysis/filters.py` (peptide/SMILES len) | ~110 | LIVE | Manual via scripts/05+09+13 |

---

## 6. WSL2 Migration Plan (recommended path for unblocking)

### 6.1 Architectural decision: hybrid Windows ↔ WSL2

**Don't migrate everything.** Keep on Windows:
- MAMMAL DTI inference (works perfectly; just needs nightly cu128)
- ADMET-AI (CPU, fast on Windows)
- Fetchers (HTTP, OS-agnostic)
- Analysis / fusion / provenance / narrative (pure pandas)

Run inside WSL2 Ubuntu:
- Boltz-2 affinity (needs cuequivariance-ops-torch-cu12)
- TxGNN inference (needs PyG cu128)
- (Future) Cluster D L1000 ops via cellxgene-census + tiledbsoma

**Data lives at `/mnt/c/Users/Pierce Lonergan/Documents/GitHub/MAMMAL_Cognitive_Enhancement_Drug_Repurposing/`** from WSL2. No data migration needed. Caches and outputs are read/written across the boundary.

### 6.2 Concrete setup steps

State on this machine (verified just now): Ubuntu distro already installed in WSL2 (stopped state). NVIDIA driver 591.59 on Windows. WSL version 2. All prerequisites met.

```bash
# Inside Windows PowerShell:
wsl -d Ubuntu                                 # boot Ubuntu
# inside Ubuntu (one-time setup):
nvidia-smi                                    # verify GPU passthrough
sudo apt update && sudo apt install -y python3.10 python3.10-venv python3.10-dev build-essential
python3.10 -m venv ~/mammal_env
source ~/mammal_env/bin/activate
pip install --pre torch torchvision --index-url https://download.pytorch.org/whl/nightly/cu128
pip install boltz cuequivariance-torch cuequivariance-ops-torch-cu12
pip install pyyaml pandas pyarrow

# Smoke test from inside Ubuntu, reading from Windows-side data:
cd "/mnt/c/Users/Pierce Lonergan/Documents/GitHub/MAMMAL_Cognitive_Enhancement_Drug_Repurposing"
python3.10 scripts/17_v2_boltz_smoke.py       # should run ~30-50s instead of 150s
```

### 6.3 Estimated effort + outcomes

- Setup time: ~30-60 min (Ubuntu + Python + pip installs)
- Disk: ~5 GB more (boltz weights re-downloaded, Linux pip cache)
- Per-pair Boltz speedup: **2.5-5× expected** (12-30 hour Phase 0.4 sweep instead of 62 hours)
- Unblocks Cluster C (PyG installs cleanly on Linux)
- Unblocks Cluster D (tiledbsoma installs cleanly on Linux)

### 6.4 Code adjustments needed

Minimal. The Windows-side `cluster_a/boltzina.py` already calls boltz via subprocess. We add a "use WSL2 boltz" code path:

```python
# In cluster_a/boltzina.py, add:
USE_WSL2_BOLTZ = os.environ.get("BOLTZ_VIA_WSL2", "0") == "1"
WSL2_DISTRO = os.environ.get("BOLTZ_WSL2_DISTRO", "Ubuntu")

def _find_boltz_executable() -> list[str]:
    if USE_WSL2_BOLTZ:
        return ["wsl", "-d", WSL2_DISTRO, "--",
                "/home/<user>/mammal_env/bin/boltz"]
    # existing Windows logic ...
```

Set `BOLTZ_VIA_WSL2=1` in the env to route Boltz calls through WSL2. The subprocess call signature stays identical; the YAML input file path needs to be `/mnt/c/...` converted from the Windows path.

---

## 7. Roadmap by Priority Tier

### Tier 1 — MUST work (gates the v2 wet-lab shortlist)

- [x] V1 MAMMAL DTI live + sanity gate (4/7 pass, peptide pollution diagnosed)
- [x] Cluster B (ADMET-AI) live with regulatory bypass
- [x] V2 RRF fusion with provenance + narrative
- [x] Boltz-2 affinity working on Windows (slow but correct)
- [ ] **Phase 0.5 CHRNA7 rescue gate** — in flight, ~30 min ETA
- [ ] WSL2 Boltz path → 2-5× speedup → full Phase 0.4 sweep feasible (~12-30 hr overnight)
- [ ] Full Phase 0.4 sweep at ~1,500 pairs
- [ ] V2 wet-lab shortlist regenerated with all 4 clusters' rankings

### Tier 2 — NICE to work (improves trustworthiness of v2)

- [ ] ChEMBL local SQLite mirror → Phase 1.2 evidence in minutes instead of hours
- [ ] Phase 3.1 calibration (per-cluster, per-target Spearman ρ)
- [ ] Phase 3.2 allosteric audit v2 (after Boltz Phase 0.4 lands)
- [ ] TxGNN running inside WSL2 against the cognition virtual anchor
- [ ] PrimeKG path scoring (compound → panel target via {drug_protein, protein_protein, pathway_protein, bioprocess_protein})
- [ ] Phase 4.1 publishable allosteric awareness benchmark (n≥10 per cell, 6 targets)

### Tier 3 — DEFER (v3 features)

- [ ] LambdaMART meta-ranker (gated on ≥20 labels)
- [ ] Cluster D (transcriptomic / L1000 via cellxgene-census) — needs WSL2 + tiledbsoma
- [ ] LLM literature-RAG layer (narrative only, not evidence)
- [ ] MCP server for interactive exploration
- [ ] Web UI / dashboard

### Tier 4 — RESEARCH (would change recommendations)

- [ ] **Cuequivariance Windows build** — if NVIDIA ever ships, drop the WSL2 detour
- [ ] **Boltzina-style Vina-pose-only mode** — implement Furui & Ohue 2025 protocol; cuts per-pair time from ~150s to ~10-15s if successful
- [ ] **TxGNN cognition-anchor validation paper** — anchor is a methodological contribution; benchmark it formally
- [ ] **MAMMAL fine-tune on cognition-relevant DTI data** — narrow MAMMAL's training distribution toward cognition targets; may resolve CHRNA7 dynamic-range collapse without needing structure

---

## 8. Risk Register

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| CHRNA7 rescue fails under Boltz | Medium | Low (publishable negative result; doc §7 explicit) | Document as "structure-aware models also fail at α7 PAM ranking"; weight CHRNA7 down in fusion |
| Full Phase 0.4 sweep too slow even with WSL2 | Low | Medium | Scope down to top-3 targets per compound; defer remainder to follow-up |
| Phase 3.1 calibration shows MAMMAL ρ < 0.3 on most targets | Medium | High (would invalidate v1's primary signal) | Restrict trustworthy panel; treat MAMMAL as a one-of-many ranker, not the anchor |
| TxGNN install on WSL2 fails (PyG version mismatch) | Medium | Medium | Pin PyG to TxGNN's tested version even at cost of slower torch |
| Boltz output schema changes in next release | Low | Low | Permissive output discovery already in place (rglob across the output tree) |
| Roberts 2020 effect-size ceiling makes any wet-lab spend hard to justify | High | High | Frame deliverable as "enrichment + provenance," not "smart drug discovery" — this is the research doc's explicit framing |

---

## 9. Single-Page Cheat Sheet — Commands That Work Today

```powershell
# Windows env (mammal_env, conda)
$envPython = "$env:USERPROFILE\.conda\envs\mammal_env\python.exe"
$repo = "C:\Users\Pierce Lonergan\Documents\GitHub\MAMMAL_Cognitive_Enhancement_Drug_Repurposing"
Push-Location $repo

# V1 pipeline (one-shot full run, all stages already done):
& $envPython scripts\02_fetch_targets.py
& $envPython scripts\03_fetch_compounds.py
& $envPython scripts\04_score_dti.py --batch-size 8 --resume
& $envPython scripts\05_sanity_check.py

# V2 ADMET + fusion (fast, ~2 min total):
& $envPython scripts\14_v2_cluster_b_admet.py
& $envPython scripts\15_v2_fusion.py

# V2 ESM2 (cached after first run, ~5 min cold):
& $envPython scripts\16_v2_esm2_embed.py

# V2 Boltz smoke (~150s; cached after first run):
$env:PYTORCH_CUDA_ALLOC_CONF = "expandable_segments:True"
& $envPython scripts\17_v2_boltz_smoke.py

# V2 Boltz focused sweep (CHRNA7+PDE4D+HRH3 panel, ~12 pairs, ~30 min):
& $envPython scripts\_boltzina_focused.py

# V2 allosteric gate (depends on boltzina_affinity.parquet):
& $envPython scripts\19_v2_allosteric_gate.py
```

```bash
# WSL2 Ubuntu side (once setup completes):
wsl -d Ubuntu
source ~/mammal_env/bin/activate
cd "/mnt/c/Users/Pierce Lonergan/Documents/GitHub/MAMMAL_Cognitive_Enhancement_Drug_Repurposing"
python3.10 scripts/17_v2_boltz_smoke.py    # expected ~30-50s with cuequivariance kernels
```
