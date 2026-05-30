# V3 Attack Plan — informed by the 4-tier deep-research findings

**Audience**: any agent (or person) picking this up at commit `d615c7d` or later. Read [V2_STATUS_AND_FORWARD_PLAN.md](V2_STATUS_AND_FORWARD_PLAN.md) first for the v2 state; this doc is the *next move* informed by [research/4-tier/](../research/4-tier/) (5 files).

---

## TL;DR (decision tree)

```
                          Where are we right now?
                                    │
        ┌───────────────────────────┼──────────────────────────┐
        ▼                           ▼                          ▼
v1 MAMMAL DTI works            ADMET-AI works            Boltz works (slow)
ESM2 cache works               (Cluster B live)          150s/pair Windows
sanity: 4/7 pos-ctrl                                      cuequivariance kernel
40 neg-ctrl hits                                          impossible on Windows
                                    │
                              Where do we go?
                                    │
        ┌───────────────┬───────────┼───────────┬─────────────────┐
        ▼               ▼           ▼           ▼                 ▼
   FIX MAMMAL      WSL2 BOLTZ   LOCAL CHEMBL   IMPLEMENT       MAMMAL
   TARGET IDs      (kernel +    SQLITE MIRROR  BOLTZINA-VINA   LORA FT
   (9 of 22 may    PyG)         (78h → sec)    PROTOCOL        ON COGNITION
    be wrong)                                  (11.8x speedup)  TARGETS
        │               │           │              │                 │
   1 hr fix,       1-3 hr       1-2 hr           1 day            3-5 days
   may resolve     setup +      setup +         (after Vina        (depends on
   CHRNA7         WSL2-VRAM     immediate       integration)      ChEMBL fix +
   collapse        risk          win                                fine-tune)
   without
   any Boltz!
```

**If only one thing can happen next**: investigate the 9 ChEMBL target ID corrections from research/4-tier/Boltzina + MAMMAL Fine-tune.md §2.3.1. If the CHRNA7 ChEMBL ID we used (`CHEMBL2492`) is wrong, the entire ChEMBL-expanded compound subset for CHRNA7 is wrong, and that's a one-hour fix vs. weeks of WSL2 + Boltz + fine-tune effort.

---

## 1. Distilled research findings (citations to /research/4-tier/)

### 1.1 cuEquivariance research → confirms WSL2, surfaces Blackwell risk

[research/4-tier/Deep Dive_ CUDA Equivariance PyTorch.md](../research/4-tier/Deep Dive_ CUDA Equivariance PyTorch.md)

| Finding | Implication |
|---|---|
| NVIDIA officially closed Windows-wheel issue #187 — no plan to support Windows | WSL2 is the *only* path |
| WSL2 GPU paravirtualization adds 0.5-1.5 GiB VRAM tax (Ampere/Ada) | Acceptable on most cards |
| **Blackwell sm_120 in WSL2: CUDA driver context = ~16 GiB** (issue microsoft/WSL #40401) | **CRITICAL — our 12 GB RTX 5070 may not fit anything else after context init** |
| Triton version conflict: PyTorch cuet path needs triton>=3.4.0 on Blackwell | Pin during setup |
| `CUEQ_TRITON_TUNING=AOT` triggers Ahead-Of-Time tuning (~hours) | Run once, cache via `CUEQUIVARIANCE_OPS_NVRTC_CACHE_DIR` |
| cuEquivariance + MACE memory: 1.4 GB → 4.6 GB with `--enable_cueq=True` (3× inflation) | Boltz likely has similar 2-3× VRAM cost with kernels |

**Decision implication**: must measure actual VRAM consumption inside WSL2 with `nvidia-smi` BEFORE committing to Phase 0.4 sweep on WSL2. If context overhead is the 16 GiB anomaly, we keep Boltz on Windows (slow) and ONLY move PyG/TxGNN to WSL2.

### 1.2 PyG Windows/WSL2 research → confirms WSL2, exact wheel commands

[research/4-tier/PyG Windows Nightly WSL2 Fixes.md](../research/4-tier/PyG Windows Nightly WSL2 Fixes.md)

- Windows native PyG on cu128 = compounding ABI mismatches + DLL pathing + MSVC errors. Unfixable in practice.
- WSL2 install (verbatim, works):
  ```bash
  pip install --pre torch torchvision --index-url https://download.pytorch.org/whl/nightly/cu128
  pip install pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv \
      -f https://data.pyg.org/whl/torch-2.7.0+cu128.html
  pip install torch-geometric
  ```
- Once PyG is installed, TxGNN install is a one-liner: `pip install git+https://github.com/mims-harvard/TxGNN.git`

### 1.3 ChEMBL local mirror → massive throughput win

[research/4-tier/ChEMBL Database Optimization Deep Dive.md](../research/4-tier/ChEMBL Database Optimization Deep Dive.md)

- The 78-hour REST backstop crisis is unnecessary. ChEMBL ships a `chembl_<release>_sqlite.tar.gz` (~4 GB) that gives instant queries against a local SQLite file.
- **Use `chembl-downloader`** (Python package): handles version-pinned downloads + extraction.
- **Optional**: install **ChemicaLite** SQLite extension for SMILES → InChIKey conversion inside SQL queries (avoids Python round-trips per row).
- Key tables to join: `activities`, `assays`, `target_dictionary`, `target_components`, `molecule_dictionary`, `chembl_id_lookup`.
- **Salt-form / virtual parent gotcha**: 30%+ of ChEMBL bioactivity records are tagged to salt forms, not the parent compound. The `molecule_dictionary.parent_chembl_id` column resolves this. Without this normalization our compound→activity joins will under-count.

**Expected throughput**: 78 hours → ~5 minutes for the full 4,713-pair lookup. Lossless.

### 1.4 Boltzina + MAMMAL fine-tune research → TWO actionable wins

[research/4-tier/Boltzina + MAMMAL Fine-tune.md](../research/4-tier/Boltzina + MAMMAL Fine-tune.md)

**§1 — Boltzina implementation (the 11.8× speedup path)**:
- Furui & Ohue 2025 protocol exact spec: Vina pose → Boltz-2 affinity head only, skip diffusion
- Suggested module: `src/mammal_repurposing/cluster_a/boltzina_vina.py`
- Requires: AutoDock Vina (apt install or precompiled), pocket DB for 22 targets (research file §1.4 lists PDB-fetched cavity centers per target)
- Wall-clock estimate on RTX 5070 + WSL2: ~10-15 s/pair (with cuequivariance kernels), vs current 150 s/pair Windows full-Boltz

**§2 — MAMMAL fine-tune (potential CHRNA7 rescue without Boltz)**:
- **§2.3.1 CRITICAL: 9 of 22 ChEMBL target mappings may need correction**
- Specifically called out: CHRNA7's correct ChEMBL ID is `CHEMBL2492` *as a single protein* but the actual nicotinic homopentamer experiments are usually tagged to `CHEMBL2107` or `CHEMBL2095192` (chicken α7 — used as model). If we expanded compound library against `CHEMBL2492`, we got the wrong subset.
- LoRA fine-tune on encoder + DTI head; ~50K parameters trainable; ~3-5 day train budget on RTX 5070
- **The cheap thing to check first**: validate the 22 ChEMBL target IDs we used in `top_binders_for_targets`. Audit each against the cognition-target intent. ~1 hour of work.

### 1.5 tiledbsoma Windows blocker → WSL2 also fixes this

[research/4-tier/tiledbsoma Windows Build Blocker.md](../research/4-tier/tiledbsoma Windows Build Blocker.md)

- Same WSL2 conclusion as everyone else
- Unblocks Cluster D (transcriptomic / cellxgene-census) which the v2 doc deferred
- No action needed unless Cluster D enters scope

---

## 2. Updated risk register (post-research)

| Risk | New evidence | Mitigation |
|---|---|---|
| **WSL2 Blackwell 16 GB context overhead** | Microsoft WSL #40401, confirmed on RTX 5090 + RTX PRO 6000 | Measure on our RTX 5070 first; fall back to Windows-Boltz-slow if real |
| **9 of 22 ChEMBL target IDs wrong** | Boltzina+MAMMAL FT research §2.3.1 | 1-hour audit script; rebuild compound library subset for affected targets |
| Boltz CPU-fallback (`--no_kernels`) on Windows produces different numeric results vs WSL2 kernel | Implicit in cuequivariance round-mode discussion (RZ vs RN precision) | Run a small overlap set on both paths; quantify ρ |
| ChEMBL local mirror is 4 GB download + 12 GB uncompressed | research §3 storage discussion | Disk OK (909 GB free in WSL2); just budget the ~10 min download |
| AutoDock Vina install failure on Windows | Common Vina issue | Vina is easy in WSL2 (`apt install autodock-vina`) |

---

## 3. The attack plan — ordered by priority × cheap-first

Tasks numbered for handoff. Each has: **What** (action), **Why** (signal value), **Files** (touch points), **Verify** (success criterion), **Owner** (any agent that picks this up).

### Tier 1 — DO NOW (cheap + high-leverage)

#### T1. Audit the 22 ChEMBL target IDs we used for compound-library expansion

**What**: Compare every `target_chembl_id` returned by `fetchers/chembl.uniprot_to_chembl_target` against literature evidence for the cognition-target intent. Cross-check against PrimeKG nodes if possible. The Boltzina + MAMMAL Fine-tune research §2.3.1 says 9 of 22 need correction; we don't know which 9.

**Why**: If CHRNA7's expansion was done against the wrong ChEMBL ID, the entire ChEMBL-expanded compound set for that target is wrong, which directly explains the dynamic-range collapse without needing structure or fine-tuning. This is *the cheapest thing that could rescue the gate*.

**Files**:
- New: `scripts/20_v3_audit_chembl_targets.py` — for each panel uniprot, fetch the full ChEMBL target list (target_type=SINGLE PROTEIN ∪ PROTEIN COMPLEX) with the matching component accession; emit a markdown report listing rejected vs accepted candidates per target.
- Update: `data/raw/targets_seed.csv` — add `chembl_target_id` column with the verified mapping.
- Reuse: `fetchers/chembl.uniprot_to_chembl_target` but with strict picking + literature cross-reference.

**Verify**: For each of the 22 targets, the chosen `chembl_target_id` either matches what we already have OR we have a written rationale for why a different ID is correct (e.g., "CHRNA7 is the human α7 nAChR; CHEMBL2107 is the homomeric pentamer functional assay used in most BindingDB records").

**Estimated time**: 1-2 hours.

**Owner**: anyone with comparative biology knowledge or willing to dig into ChEMBL.

---

#### T2. Build the local ChEMBL SQLite mirror (Phase 1.2 v3)

**What**:
1. `pip install chembl-downloader` (3-line install, MIT)
2. `chembl_downloader.download_sqlite(version="latest")` — downloads ~4 GB, extracts ~12 GB
3. Rewrite `fetchers/chembl_groundtruth.lookup_pair` to use SQLite queries instead of REST
4. Single SQL JOIN replaces the per-pair InChIKey lookup + activity lookup

**Why**: 78-hour backstop → seconds. Critical for Phase 3.1 calibration. Also enables much larger expansion sets later.

**Files**:
- New: `scripts/21_v3_chembl_mirror.py` — one-shot downloader script
- New: `src/mammal_repurposing/fetchers/chembl_sqlite.py` — replacement query module
- Update: `scripts/07_chembl_evidence.py` to use the new SQLite path with `--mirror` flag
- Schema reference for the join (research §4):
  ```sql
  SELECT md.chembl_id AS molecule_id, td.chembl_id AS target_id,
         a.standard_type, a.standard_value, a.standard_units, a.activity_comment
  FROM activities a
  JOIN assays asy ON a.assay_id = asy.assay_id
  JOIN target_dictionary td ON asy.tid = td.tid
  JOIN target_components tc ON td.tid = tc.tid
  JOIN component_sequences cs ON tc.component_id = cs.component_id
  JOIN molecule_dictionary md ON a.molregno = md.molregno
  WHERE cs.accession = ? AND md.parent_chembl_id IN (?, ?, ...);
  ```
- Salt-form gotcha: include `parent_chembl_id` resolution, not just `chembl_id`

**Verify**: Re-run the full 4,713-pair backstop. Wall-clock < 10 minutes (vs 78 hours). Label breakdown (CORROBORATED / NOVEL / CONTRADICTED) should match the REST results for the subset that already finished.

**Estimated time**: 2-3 hours including download.

**Owner**: any agent comfortable with SQL.

---

#### T3. Measure the Blackwell-WSL2 VRAM context overhead

**What**: After WSL2 setup completes, run a minimal CUDA context init and observe `nvidia-smi` consumption. Specifically:
```bash
nvidia-smi  # before
python3 -c "import torch; torch.cuda.init(); input('press enter')"  # while process is paused
nvidia-smi  # measure consumption
```
If the context eats >2 GB on our 12 GB card, we're hitting the Blackwell anomaly and Boltz inside WSL2 will OOM.

**Why**: This is the GO/NO-GO for WSL2 Boltz. If context overhead is >4 GB on the 12 GB card we don't have room for the model + activations.

**Files**:
- New: `scripts/_wsl2_vram_probe.sh` — Bash script invoking the above measurement
- Update: `design/architecture-and-plans/V2_STATUS_AND_FORWARD_PLAN.md` — record the empirical number

**Verify**: A single number in GB. If ≤2 GB, proceed with WSL2 Boltz. If 5-8 GB, scope down Boltz to smaller batch sizes. If ≥10 GB, abandon WSL2 Boltz; keep on Windows.

**Estimated time**: 5 minutes (after WSL2 setup is done).

**Owner**: anyone.

---

### Tier 2 — DO SOON (per-tool unblocks)

#### T4. WSL2 Boltz speedup smoke (assuming T3 passes)

**What**: Run `_wsl2_boltz_smoke.sh` already committed. Compare wall-clock vs the 149.8s Windows baseline.

**Why**: Quantifies the actual speedup. Expectation: 30-60s (2.5-5× per research doc).

**Verify**: pair scored, JSON parses, ETA for full 1,500-pair sweep recomputed.

**Estimated time**: 10 min (with weights pre-cached at /home/user/.boltz).

**Owner**: anyone.

---

#### T5. Implement Boltzina-Vina-only mode (the 11.8× speedup)

**What**: Per Boltzina + MAMMAL FT research §1.2-1.3, build `cluster_a/boltzina_vina.py`. Workflow:
1. AutoDock Vina (apt: `autodock-vina`) generates pose given pocket coords
2. Feed pose (as Boltz YAML with `--no_diffusion` equivalent) to Boltz-2's affinity head only
3. Skip the recycling iterations the structure pred uses

**Why**: Even on Windows-slow Boltz path, this cuts per-pair time from 150s to ~13s (research's 11.8× claim). On WSL2 with kernels could go to ~3s/pair → full 6,556-pair sweep in ~5.5 hours, finally feasible.

**Files**:
- New: `src/mammal_repurposing/cluster_a/boltzina_vina.py`
- New: `data/cache/vina_pockets/<uniprot>.json` — pocket coords for all 22 targets (research §1.4 has a script outline)
- Update: `cluster_a/boltzina.py` — add `mode="boltzina_vina"` path (it's currently NotImplementedError)

**Verify**: Same (CHRNA7, galantamine) smoke pair returns binder_prob/affinity within ±0.2 of full Boltz, in ≤30s.

**Estimated time**: 1 day (Vina integration is the bulk).

**Owner**: agent with structural biology familiarity for the pocket coord step.

---

#### T6. PrimeKG + TxGNN in WSL2 (Cluster C live)

**What**:
1. WSL2 venv (already created via T0 setup)
2. `pip install pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv -f https://data.pyg.org/whl/torch-2.7.0+cu128.html`
3. `pip install torch-geometric`
4. `pip install git+https://github.com/mims-harvard/TxGNN.git`
5. Download PrimeKG from Harvard Dataverse (DOI 10.7910/DVN/IXA7BM) to `/mnt/c/.../data/kg/primekg/`
6. Implement the real `cluster_c/{primekg,txgnn}.py` (stubs ready)

**Why**: Unblocks the entire mechanism cluster. Lets the cognition virtual anchor produce real indication scores. Fourth ranker into RRF makes the disagreement archetypes meaningful.

**Files**:
- Update: `src/mammal_repurposing/cluster_c/primekg.py` — implement real `load_primekg` + `score_compound_paths` (PPR or Katz)
- Update: `src/mammal_repurposing/cluster_c/txgnn.py` — implement real `score_compounds_against_anchor`
- New: `scripts/22_v3_kg_scores.py` — orchestrator that runs Cluster C end-to-end
- Update: `scripts/15_v2_fusion.py` — pick up txgnn parquet when present

**Verify**: TxGNN places donepezil, memantine, BPN14770 in top decile against the cognition anchor (the validation gate per the v2 research doc §3 Class C).

**Estimated time**: 1-2 days, mostly install fights + PrimeKG download.

**Owner**: agent comfortable with PyG.

---

### Tier 3 — DO LATER (Phase 0.5 + sweep)

#### T7. WSL2-Boltz Phase 0.4 sweep (full ~1,500 pairs)

Depends on T4 passing + Boltzina-Vina (T5) for the speedup. Run overnight as background.

#### T8. Phase 0.5 allosteric gate verdict

Re-runs `scripts/19_v2_allosteric_gate.py` against the WSL2-Boltz parquet.

#### T9. v2 fusion regenerated with all 4 clusters

Run `scripts/15_v2_fusion.py` once Boltz + TxGNN parquets exist. Produces the publishable wet-lab shortlist v2.

#### T10. Phase 4.1 allosteric benchmark expansion (publishable)

n≥10 per (target, binding_mode) cell × 6+ targets. Score everything through all 4 rankers. Methods note.

### Tier 4 — RESEARCH (multi-week)

#### T11. MAMMAL LoRA fine-tune on cognition-DTI corpus

Per Boltzina + MAMMAL FT research §2. Only meaningful AFTER T1 (chembl target audit) is done — otherwise the fine-tune corpus is built on the same wrong target IDs.

#### T12. Cluster D (transcriptomic / cellxgene-census)

Only after WSL2 is stable AND a neuron-relevant L1000 dataset is curated.

---

## 4. Concrete handoff state (as of commit d615c7d)

### What's committed and runnable today

```
v1 pipeline (Windows): 02-05 + 09 + 11 + 12 + 13 — full end-to-end works
v2 Cluster B (Windows): 14 + 15 — RRF over MAMMAL + ADMET, provenance + narrative
v2 Cluster A ESM2 (Windows): 16 — 22 targets cached
v2 Cluster A Boltz (Windows): 17 smoke + 18 full sweep + 19 gate; slow (150s/pair) but correct
v2 ADMET-AI: 14 — 53 PASS / 64 FLAG / 181 CUT; regulatory bypass works
WSL2: Ubuntu 24.04 + RTX 5070 passthrough confirmed; setup script in flight
```

### What's running in background right now

- WSL2 setup script (`_wsl2_setup_boltz.sh`) at Step 1: apt install
- Windows focused Boltz sweep (`_boltzina_focused.py`) at 5/10 pairs, ETA 10 min
- Windows ChEMBL v2 backstop (`07_chembl_evidence.py --threshold 6.5`) at ~target 9/22, ETA hours

### What an agent picking this up should do FIRST (sequence)

1. Read this doc (V3_ATTACK_PLAN.md) + V2_STATUS_AND_FORWARD_PLAN.md
2. Check the latest commit on `main` for any updates
3. Verify background tasks: `ls data/results/v2/boltzina_affinity.parquet` (does the focused sweep finish?), `ls data/results/chembl_evidence.parquet` (did the v2 backstop finish at 6.5 threshold?)
4. Choose ONE Tier 1 task (T1 / T2 / T3) — they're independent, all cheap, all high-leverage
5. If WSL2 setup completed: T3 (VRAM probe) is 5 minutes of work that resolves a huge architecture question

### Decision waypoints (these are the things that change the plan)

- **If T1 finds CHRNA7 ChEMBL ID was wrong**: rebuild compound library subset, re-run scripts/04 + 05, expect CHRNA7 dynamic range to improve. Maybe no Boltz needed.
- **If T3 measures >5 GB Blackwell context in WSL2**: abandon WSL2 Boltz. Keep Boltz on Windows. Use WSL2 only for PyG/TxGNN (Cluster C).
- **If T2 (local ChEMBL mirror) lands fast**: Phase 3.1 calibration suddenly becomes free; can compute ρ per target and re-weight RRF in `configs/weights_calibrated.yaml`.
- **If T5 (Boltzina-Vina) achieves <20 s/pair on Windows**: WSL2 migration is optional; keep simpler Windows-only stack.

---

## 5. The mission statement (for whoever picks this up)

This pipeline is an honest attempt to use a foundation model (MAMMAL) plus a curated panel of cognition-relevant targets plus rigorous filtering (ADMET + structural + mechanistic) to surface drug-repurposing candidates with mechanistic plausibility and physical feasibility for healthy adult cognitive enhancement.

The Roberts 2020 meta-analysis ceiling is real: even methylphenidate hits SMD=0.21 on a generous metric. We are not searching for a miracle drug; we are enriching a candidate set so wet-lab cycles spend money on plausibility, not chemistry-lottery tickets.

The v2 already does this better than v1 — the top 10 went from anonymous ChEMBL artifacts to actual CNS drugs. The v3 work in this attack plan adds:
- Mechanistic context (Cluster C live) so disagreement archetypes have meaning
- Structural feasibility (faster Boltz / Boltzina-Vina) so allosteric blindness gets a second opinion
- Calibration against ChEMBL ground truth (local SQLite) so we know which targets to trust
- Possibly: a fix to the panel that didn't need any of the above (T1: ChEMBL target ID audit)

If you're picking this up: T1 is one hour. T3 is five minutes. Either one might shrink the entire remaining roadmap. Do those first. Then go heavier.

Life would be a lot better if people had better cognition. Let's get on with it.
