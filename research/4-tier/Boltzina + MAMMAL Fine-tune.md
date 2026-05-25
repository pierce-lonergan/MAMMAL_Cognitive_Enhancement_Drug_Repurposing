# Boltzina + MAMMAL Fine-tune: A Deep Technical Roadmap for the V3 Cognition Pipeline

## TL;DR
- **Boltzina is real, vendored, and ready to plug in.** The reference implementation `ohuelab/boltzina` (MIT, 74 stars per the ohuelab GitHub organization repository listing as of May 2026) implements exactly the Furui & Ohue 2025 protocol — AutoDock Vina 1.2.7 → PDB → MMCIF (pdb-tools + MAXIT) → Boltz-2 affinity head with `recycling_steps=1`. Measured speedup on H100 + 48 CPU cores: Boltz-2 16.5 s/ligand → Boltzina 2.3 s (7.3×) → Boltzina (Cycle=1) 1.4 s (11.8×). Expected wall-clock on RTX 5070 + WSL2 for the 1,500-pair Phase 0.4 sweep: **4–6 hours**, not 62. **Critical caveat:** RTX 5070 is Blackwell sm_120; per the `boltz-blackwell` v0.1.2 PyPI README "The cuEquivariance library doesn't support sm_121 yet. Always use --no_kernels." So the expected 2–5× kernel speedup will not land on RTX 5070 today — plan for `--no_kernels` and validate on a 10-pair benchmark before scaling.
- **MAMMAL LoRA fine-tune on a cognition-DTI corpus is feasible but bounded.** The ChEMBL panel yields ~35–45k unique target-compound DTI pairs across the 22 targets (18 of 22 clear the ≥100 high-confidence threshold; GRIA3, GRIA4, GRIN2A, HCN1 are data-sparse). A LoRA rank-16 QKVO adapter on the T5 encoder-decoder is the realistic config — ~5M trainable params (≈1.2% of 458M), fits in 12 GB with bf16 + gradient checkpointing + batch 8, ~3–4 hour overnight train. **Nine of Pierce's supplied ChEMBL IDs were wrong** (CHEMBL1968 = EPHX1, not GRIN2B; CHEMBL4051 = CFTR, not KCNQ2; CHEMBL2243 = FAAH, not GRIA3; CHEMBL3231 = ROCK1, not KCNQ3; etc.) — the corrected mapping is in Topic 2 §2.3.1.
- **Order of operations:** Build Boltzina in WSL2 first (1–2 days, low risk, large throughput win). Then attempt LoRA in parallel (1 week, moderate risk, addresses the *signal-quality* problem). The fine-tune cannot conjure structural reasoning; if it fails to rescue CHRNA7 PAMs (galantamine, encenicline, TC-5619) into the top 25%, **abandon the MAMMAL allosteric ambition and rely on Boltzina/Boltz-2 structural rescoring for all allosteric-rich targets** (CHRNA7, GRIA1–4, GRIN2A/2B, PDE4D, KCNQ2/3, NTRK2).

---

## Key Findings

### Boltzina (Topic 1)
1. **Algorithm**: Vina 1.2.7 with 20 Å grid, exhaustiveness 8, best pose only (for the main MF-PCBA experiments). Pose is converted PDB → MMCIF via `pdb-tools` + RCSB `MAXIT`, then loaded as a Boltz-2 template structure that *replaces* the diffusion-sampled structure normally passed to the affinity PairFormer. The structure module is fully bypassed.
2. **Recycling**: Boltz-2 default in the affinity path is 5; Boltzina default is also 5 (yields 7.3× speedup just from skipping the diffusion structure module + batching > 1); Boltzina (Cycle=1) is 1 cycle (11.8×). Mean Average Precision: Boltz-2 0.084 → Boltzina 0.056 → Boltzina (C=1) 0.048 → Boltzina (No-Pose, ligand at origin) 0.043. **Most of Boltzina's accuracy comes from the trunk module's PairFormer**, not the pose — even no-pose mode is competitive with GNINA. This is critical: it means accuracy degradation on allosteric pockets is bounded by the trunk's protein-sequence + SMILES representation, not by Vina pose quality.
3. **Implementation surface** (`github.com/ohuelab/boltzina`): three Python entry points — `ligand_preparation.py` (SMILES → PDB), `run.py` (config-driven orchestrator), `boltzina_main.py` (the `Boltzina` class). The affinity scoring is in `boltzina/affinity/predict_affinity.py` exposing `load_boltz2_model()` and `predict_affinity()`. Two modes: full docking, and scoring-only (pre-existing poses, bypasses Vina entirely). MIT licensed — vendor freely.
4. **Pose strategy**: Top-5 Average across Vina poses gave the best ROC-AUC (0.778 vs 0.746 for Best-Pose-Only, p<0.01 Wilcoxon). The Top-N-Best-Score strategy did NOT outperform averaging — Boltz-2's binding likelihood does not discriminate individual pose quality. **Recommendation for Pierce**: use Top-5 Average for the cognition panel; the extra 4 affinity calls per ligand (~5–6 s) is cheap compared to a Vina re-run.
5. **Two-stage screening**: Boltzina top 20% → Boltz-2 rescore achieves mean AP > 0.07 at ~3× faster than Boltz-2 alone (Pareto-optimal). This is a built-in option for Phase 0.5 — use Boltzina to rank all 1,500 pairs, then re-run full Boltz-2 only on the top 300.

### MAMMAL fine-tune (Topic 2)
1. **Architecture**: `ibm/biomed.omics.bl.sm.ma-ted-458m.dti_bindingdb_pkd` is a custom `Mammal` class wrapping a 351M-parameter T5 encoder-decoder (`T5ForConditionalGeneration`) + a 106M-parameter `ClassifierMLP` scalar head (total 458M). Scalar pKd prediction goes through `model.forward_encoder_only(...)`. This is the LoRA target: 36 attention blocks × 4 projections (q/k/v/o) = 144 attention-projection modules. The repo `BiomedSciAI/biomed-multi-alignment` (Apache-2.0, 44 stars per the repo's GitHub Actions workflow page accessed May 2026) ships a working DTI fine-tune example at `mammal/examples/dti_bindingdb_kd/` invoked via `python mammal/main_finetune.py --config-name config.yaml --config-path examples/dti_bindingdb_kd`. It uses TDC's BindingDB_Kd dataset — 52,284 pairs / 10,665 drugs / 1,413 proteins per the Therapeutics Data Commons page (tdcommons.ai/multi_pred_tasks/dti/) — with Drug+Target cold split and mean/std pKd normalization.
2. **LoRA pattern is proven on this exact model**: the `Kymi808/mammal-lora-bbbp` project applies HuggingFace `peft` to MAMMAL's inner T5 attention projections (q, k, v, o) for BBBP classification on Apple Silicon. Per its SUMMARY.md: "We applied LoRA adapters via peft to MAMMAL's inner T5 attention projections (q, k, v, o — 36 attention blocks × 4 projections = 144 modules)." The same hook approach works for the DTI head.
3. **ChEMBL corpus**: across the 22 cognition targets, ~50–70k raw binding records reducing to **~35–45k unique compound-target pairs** under the standard DTI filter. The corpus is dominated by ACHE (~7k records, ~5.8k compounds), HRH3 (~6–8k), SIGMAR1 (~3.5–5.5k), DAT/SLC6A3 (~5–7k), and CHRNA7 (~3.5–5k). Eight allosteric-rich targets (CHRNA7, GRIA1–4, GRIN2A/2B, PDE4D, KCNQ2/3, NTRK2) provide the test bed for the dynamic-range repair.
4. **Catastrophic forgetting risk is real but bounded**: the BBBP MAMMAL-LoRA project saw the base model's BBBP AUROC of 0.9745 (cross-validated) hold up after LoRA touched only attention projections. For DTI regression, holding QKV+O adapters on the encoder and freezing the scalar head body until late training is a defensible default.
5. **The fine-tune addresses the *symptom* (compressed dynamic range), not the *root cause* (no pocket awareness)**. Expect the LoRA to rescue rankings on targets where ChEMBL has dense, well-labeled allosteric SAR (CHRNA7, PDE4D, GRIN2B) and to *fail to rescue* on the AMPA receptors (GRIA1–4 are subtype-sparse in ChEMBL — even after ID correction, GRIA3/4 have only ~50–150 records each). Pre-commit to falling back to Boltzina structural rescoring for the AMPA subtypes regardless of fine-tune outcome.

---

## Details

# TOPIC 1 — Boltzina Vina-Pose-Only Mode Implementation

## 1.1 Exact protocol from Furui & Ohue 2025 (arXiv:2508.17555)

The paper is 11 pages, NeurIPS 2025 AI4Mat workshop accepted (openreview.net/forum?id=OwtEQsd2hN). The Boltzina pipeline:

```
SMILES + Protein Sequence
   │
   ├─ Boltz-2 holo prediction (once per target) ────► reference complex
   │                                                       │
   │                                                       └─ centroid → Vina grid center
   │
   ├─ ligand_preparation.py: SMILES → 3D PDB (RDKit + OpenBabel)
   │
   ├─ AutoDock Vina v1.2.7, grid 20 Å, exhaustiveness 8 ─► best (or top-N) pose as PDBQT
   │
   ├─ PDB-tools + MAXIT (RCSB) ──► MMCIF complex (reusing Boltz-2's template structure processing)
   │
   └─ Boltz-2 affinity head, structure module SKIPPED,
        recycling_steps = 1 (Cycle=1 mode) or 5 (default)
        ──► affinity_pred_value (log10 IC50 μM), affinity_probability_binary
```

**Reported numbers (8 MF-PCBA assays, H100, 48 CPU cores, 1,000 ligands sampled for timing):**

| Method | Mean AP | Time / ligand | Speedup vs Boltz-2 |
|---|---|---|---|
| Boltz-2 (default) | 0.084 | 16.5 s | 1.0× |
| Boltzina (default, 5 recycles) | 0.056 | 2.3 s | **7.3×** |
| Boltzina (Cycle=1) | 0.048 | 1.4 s | **11.8×** |
| Boltzina (No Pose, ligand at origin) | 0.043 | ~1.0 s | (debug) |
| GNINA v1.3.2 | ~0.01 | — | — |
| AutoDock Vina alone | ~0.01 | 0.8 s | — |

Vina docking takes 0.8 s/ligand and is the rate-limiter in Cycle=1 mode. The 11.8× number is the ceiling.

**Direct paper quote on Cycle=1 trade-off (Furui & Ohue 2025, §3.1):** "For Boltzina (Cycle=1), the mean AP decreased from 0.056 to 0.048, but the decrease was limited, confirming that reducing the number of recycling iterations is effective for cutting computational cost without significantly compromising accuracy."

## 1.2 Implementation entry points

Recommended approach: **vendor `ohuelab/boltzina` as a git submodule under `vendor/boltzina/` and write a thin wrapper.** Do not try to reimplement — the MMCIF conversion logic alone (template-structure-as-pose injection into Boltz-2's input pipeline) is non-trivial.

```
mammal_repurposing/
├── vendor/
│   └── boltzina/          # git submodule, ohuelab/boltzina @ pinned commit
└── src/mammal_repurposing/cluster_a/
    └── boltzina_vina.py   # thin wrapper (below)
```

The vendored `boltzina/affinity/predict_affinity.py` exposes:
- `load_boltz2_model(...)` — returns a `Boltz2` instance with the affinity head wired
- `predict_affinity(...)` — batched scoring

The `boltzina_main.py` `Boltzina` class is the orchestrator (DeepWiki summary: "manages the complete workflow from receptor preparation through final result generation"; the relevant lines are 22–106 for the class declaration and 169–266 for the Boltz-2 scoring system per the deepwiki.com/ohuelab/boltzina index).

## 1.3 Concrete Python module: `src/mammal_repurposing/cluster_a/boltzina_vina.py`

```python
"""Boltzina wrapper — Vina-pose-only Boltz-2 scoring.

Vendored from ohuelab/boltzina @ <pinned-sha>. MIT licensed.
"""
from __future__ import annotations
import hashlib, json, subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import torch
from rdkit import Chem
from rdkit.Chem import AllChem

# from vendored boltzina:
from boltzina.affinity.predict_affinity import load_boltz2_model, predict_affinity
from boltzina.data.parse.mmcif import pdbqt_to_mmcif_template

CACHE_ROOT = Path("~/cache/boltzina").expanduser()
CACHE_POSE = CACHE_ROOT / "vina_poses"
CACHE_AFFY = CACHE_ROOT / "affinity"
for p in (CACHE_POSE, CACHE_AFFY):
    p.mkdir(parents=True, exist_ok=True)


@dataclass
class PocketSpec:
    """Pocket definition. Either explicit residue list OR centroid coords."""
    residues: Optional[list[int]] = None
    center: Optional[tuple[float, float, float]] = None
    box_size: tuple[float, float, float] = (20.0, 20.0, 20.0)
    exhaustiveness: int = 8


class VinaPoseGenerator:
    """SMILES + receptor (PDB) → top-N Vina poses (PDBQT). Uses Meeko 0.5+."""
    def __init__(self, vina_bin: str = "vina",
                 mk_prep: str = "mk_prepare_receptor.py",
                 mk_lig: str = "mk_prepare_ligand.py"):
        self.vina = vina_bin
        self.mk_prep = mk_prep
        self.mk_lig = mk_lig

    def prepare_receptor(self, pdb: Path, pocket: PocketSpec) -> Path:
        pdbqt = pdb.with_suffix(".pdbqt")
        cx, cy, cz = pocket.center
        sx, sy, sz = pocket.box_size
        # Meeko mk_prepare_receptor.py — current best for non-standard residues
        subprocess.run([
            self.mk_prep, "-i", str(pdb), "-o", str(pdbqt.with_suffix("")),
            "-p", "--box_center", f"{cx}", f"{cy}", f"{cz}",
            "--box_size", f"{sx}", f"{sy}", f"{sz}"
        ], check=True)
        return pdbqt

    def prepare_ligand(self, smiles: str, out: Path) -> Path:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        sdf = out.with_suffix(".sdf")
        Chem.SDWriter(str(sdf)).write(mol)
        pdbqt = out.with_suffix(".pdbqt")
        subprocess.run([self.mk_lig, "-i", str(sdf), "-o", str(pdbqt)], check=True)
        return pdbqt

    def dock(self, receptor_pdbqt: Path, ligand_pdbqt: Path,
             pocket: PocketSpec, num_modes: int = 5) -> Path:
        out = ligand_pdbqt.with_name(ligand_pdbqt.stem + "_docked.pdbqt")
        cx, cy, cz = pocket.center
        sx, sy, sz = pocket.box_size
        subprocess.run([
            self.vina, "--receptor", str(receptor_pdbqt),
            "--ligand", str(ligand_pdbqt),
            "--center_x", f"{cx}", "--center_y", f"{cy}", "--center_z", f"{cz}",
            "--size_x", f"{sx}", "--size_y", f"{sy}", "--size_z", f"{sz}",
            "--exhaustiveness", str(pocket.exhaustiveness),
            "--num_modes", str(num_modes),
            "--out", str(out)
        ], check=True, timeout=120)
        return out


class BoltzAffinityHead:
    """Loads Boltz-2 once, scores arbitrary MMCIF complexes via the affinity head."""
    def __init__(self, ckpt: Path | str = "default",
                 device: str = "cuda", recycling_steps: int = 1,
                 use_kernels: bool = False):  # False on RTX 5070 (Blackwell sm_120)
        self.model = load_boltz2_model(ckpt, device=device, use_kernels=use_kernels)
        self.recycling_steps = recycling_steps

    @torch.no_grad()
    def score(self, mmcifs: list[Path], batch_size: int = 8) -> list[dict]:
        return predict_affinity(self.model, mmcifs,
                                recycling_steps=self.recycling_steps,
                                batch_size=batch_size)


class BoltzinaPipeline:
    def __init__(self, pocket_db: dict[str, PocketSpec],
                 receptor_db: dict[str, Path],
                 affinity_head: BoltzAffinityHead,
                 vina: VinaPoseGenerator,
                 pose_strategy: str = "top5_avg"):
        self.pocket_db = pocket_db
        self.receptor_db = receptor_db
        self.head = affinity_head
        self.vina = vina
        self.pose_strategy = pose_strategy

    def _cache_key(self, target_uniprot: str, smiles: str) -> str:
        h = hashlib.sha256(f"{target_uniprot}|{smiles}".encode()).hexdigest()[:16]
        return h

    def score(self, target_uniprot: str, drug_smiles: str) -> dict:
        key = self._cache_key(target_uniprot, drug_smiles)
        aff_cache = CACHE_AFFY / f"{key}.json"
        if aff_cache.exists():
            return json.loads(aff_cache.read_text())

        receptor = self.receptor_db[target_uniprot]
        pocket = self.pocket_db[target_uniprot]

        receptor_pdbqt = self.vina.prepare_receptor(receptor, pocket)
        ligand_pdbqt = self.vina.prepare_ligand(drug_smiles, CACHE_POSE / key)
        docked = self.vina.dock(receptor_pdbqt, ligand_pdbqt, pocket, num_modes=5)

        # Split top-5 poses → 5 MMCIFs
        mmcifs = pdbqt_to_mmcif_template(docked, receptor, top_n=5)
        scores = self.head.score(mmcifs, batch_size=5)

        if self.pose_strategy == "top5_avg":
            aff = float(sum(s["affinity_pred_value"] for s in scores) / len(scores))
            prob = float(sum(s["affinity_probability_binary"] for s in scores) / len(scores))
        else:  # best
            best = min(scores, key=lambda s: s["affinity_pred_value"])
            aff, prob = best["affinity_pred_value"], best["affinity_probability_binary"]

        result = {
            "affinity_pred_value": aff,           # log10(IC50 in μM); lower = tighter
            "affinity_probability_binary": prob,  # binder probability [0,1]
            "pose_plddt": float(scores[0].get("plddt", 0.0)),
            "n_poses_scored": len(scores),
        }
        aff_cache.write_text(json.dumps(result))
        return result
```

**Caching strategy**: Unified `sha256("uniprot|smiles")[:16]` keyed JSON cache. Separate Vina pose cache by the same key for re-runs that change only Boltz-2 parameters.

## 1.4 Pocket database for the 22-target cognition panel

Hardcode the validated pocket centroids in `mammal_repurposing/data/pockets.py`. Concrete values from published structures (centroids derived from co-crystal ligand center of geometry, computable via `pymol` or by averaging chain-A CA coords of the listed residues):

| Target | PDB | Pocket | Validated residues (anchor points) |
|---|---|---|---|
| ACHE | 4EY7 | catalytic gorge (CAS+PAS) | Catalytic triad S203/H447/E334; gorge Y72, D74, W86, G121, Y124, E202, S203, N265, E268, W286, S293, V294, F295, Y337, F338, Y391, H447, G448 (per CASTp analysis of 4EY7) |
| CHRNA7 | 7EKT | TMD intersubunit PAM site (PNU-120596) | TM1–TM3 cavity. Greatest functional impact at A225 (TM1) and M253 (TM2) per Young et al. 2008 PNAS; also W148 / W54 are orthosteric anchors to avoid |
| PDE4D | 6NJJ | catalytic + UCR2 allosteric | T275 (UCR2), F271 / Y271 (subtype-selectivity switch — see Gurney et al. 2019 J Med Chem), S449, N450, Q451, H506, Mg²⁺ pocket |
| GluA2 (GRIA2) | 1LBC/3RN8 | LBD cyclothiazide dimer-interface site | S754, L751, S729 (dimer interface) |
| HRH3 | AF2 model (no cryo-EM yet) | orthosteric biogenic-amine pocket | D114(3.32), W402(6.48), Y374(6.51) |
| DRD1 | 7CKZ/7JV5 | orthosteric + extracellular vestibule (LY3154207 PAM) | D103(3.32) orthosteric + extracellular loop allosteric |
| SLC6A3 (DAT) | 4M48 (dDAT) homology | central S1 site | D79, F319, F325, Y156 |
| SLC6A2 (NET) | AF2 / homology to DAT | central S1 site | D75 (analogue to DAT D79) |
| ADRA2A | 6KUX | orthosteric | D113(3.32), F412(6.51) |
| HCRTR1 | 4ZJC | orthosteric (suvorexant pocket) | Q126(3.32), N318(6.55), Y311(6.48) |
| HCRTR2 | 4S0V/5WQC | orthosteric | Q134(3.32), H350(6.55) |
| GRIN2B | 3QEL / 5EWJ | ifenprodil ATD allosteric | F176, F114, P78, E236 (ATD dimer interface) |
| GRIN2A | 5KCJ | TCN-201 ATD allosteric | analogous to GRIN2B |
| PDE9A | 4G2L | catalytic | Y424, F441, Mg²⁺ pocket |
| NTRK2 | 4AT5 | kinase domain | hinge M635, gatekeeper F633 |
| SIGMAR1 | 5HK1 | small hydrophobic ligand pocket | D126, E172, Y103 |
| KCNQ2 | 7CR0/7CR2 | retigabine VSD-pore PAM site | W236 (S5), L243 (S5) |
| KCNQ3 | 6V00 | retigabine PAM site | W265 (analogous tryptophan anchor) |
| HCN1 | 5U6P | central cavity | C347 (analogous to ivabradine binding) |
| GRIA1/3/4 | homology to GRIA2 | LBD allosteric (assumed conserved) | conserved S754 site |

For each target, derive the centroid via:
```python
from Bio.PDB import PDBParser
parser = PDBParser(QUIET=True)
struct = parser.get_structure("X", pdb_path)
residues_of_interest = [203, 447, 334]  # ACHE example
coords = [atom.coord for res in struct[0]['A']
          if res.id[1] in residues_of_interest
          for atom in res if atom.name == 'CA']
import numpy as np
center = tuple(map(float, np.mean(coords, axis=0)))
```

## 1.5 Wall-clock estimate on RTX 5070 + WSL2

**Per pair on H100 (paper baseline):** 0.8 s Vina + 0.6 s Boltz-2 affinity (C=1) = 1.4 s.

**RTX 5070 vs H100 — concrete FLOPS data:**
- RTX 5070 (Blackwell GB205): 30.87 TFLOPS FP16 on CUDA cores (per NVIDIA Blackwell whitepaper Table 6 as compiled by gpupoet.com); 61.73 TFLOPS Tensor Core FP16 dense per WareDB.
- H100 SXM5: 989 FP16/BF16 TFLOPS Tensor Core dense (per the NVIDIA H100 datasheet as compiled at spheron.network/blog/nvidia-h100-specs/), 1,979 TFLOPS with 2:4 sparsity.
- Raw Tensor-Core ratio: ~16× on dense FP16. For transformer inference dominated by memory bandwidth the gap collapses to ~6–10× in practice; for the Boltz-2 PairFormer with its triangle attention the gap is closer to the theoretical ratio because the operation is FLOPS-bound when kernels are present (and the kernels are *not* present on sm_120, raising it further).

**Worst case (no kernels, single batch)**: 1.4 s × 8 = 11 s/pair → 1,500 pairs × 11 s = **4.6 hours**.
**Best case (kernels work hypothetically, batched 8)**: 1.4 s × 4 = 5–6 s/pair → **2 hours**.

**Bottleneck reality check**: Vina at exhaustiveness 8 on a consumer Ryzen/Intel with 16 threads runs ~3–5 s/ligand (not 0.8 s — the paper's 0.8 s used 48 parallel processes on an H100 host). So Pierce should expect Vina at 5–10 s/ligand on his hardware and Boltz-2 affinity at 5–10 s/pair. Net: **~15 s/pair × 1,500 = 6 hours overnight**. Still a >10× improvement vs the 62-hour Windows full-Boltz-2 baseline.

## 1.6 Validation protocol (positive controls)

Run these BEFORE any cognition-panel sweep. Acceptance: predicted log10(IC50 μM) within ±1.0 of experimental, and binary probability > 0.5.

| Target | PDB | Ligand | Experimental | Accept band (log10 μM) |
|---|---|---|---|---|
| ACHE | 4EY7 | donepezil | IC50 30 nM (human) | -1.5 expected; band [-2.5, -0.5] |
| CHRNA7 | 7EKT | PNU-120596 | EC50 ~200 nM | -0.7 expected; band [-1.7, +0.3] |
| PDE4D | 6NJJ | BPN14770 | Ki ~7 nM | -2.2 expected; band [-3.2, -1.2] |
| HRH3 | AF2 | pitolisant | Ki ~0.16 nM | -3.8 expected; band [-4.8, -2.8] |
| SLC6A3 | 4M48 homology | methylphenidate | Ki ~390 nM | -0.4 expected; band [-1.4, +0.6] |

If donepezil/ACHE returns log10(IC50) > 1.0 (i.e., > 10 μM predicted), the MMCIF template injection is broken — debug atom-mapping before proceeding.

## 1.7 Failure modes

1. **`pdbqt_to_mmcif_template` atom-ordering bug**: Vina output PDBQT lists atoms in a different order than RDKit canonical SMILES. MAXIT re-derives connectivity, but if the SMILES that Boltz-2 tokenizes has a different RDKit canonical order than the MMCIF residue list, the affinity head sees a mismatched ligand. **Sanity check**: log `len(rdkit_atoms)` vs `len(mmcif_ligand_atoms)` and assert equality.
2. **Boltz-2 ligand atom limit**: Per the Boltz-2 prediction docs, the affinity head is hard-limited to 128 heavy + RDKit-retained atoms and the official recommendation is "we do not recommend running the affinity module with ligands significantly larger than 56 atoms (counted as above, limit set during training)." Several of the ~298-compound library are likely > 56 atoms (any peptide-like cognition enhancer, e.g., Cerebrolysin fragments). **Pre-filter** the library before the Boltzina sweep.
3. **Allosteric pocket selection failure on AMPA**: GluA2 LBD has multiple co-existing sites (orthosteric Glu, cyclothiazide dimer-interface allosteric, perampanel TMD). Vina with a 20 Å box centered on cyclothiazide will not dock orthosteric ligands correctly. **Mitigation**: define *two* pockets per AMPA target (orthosteric + allosteric) and take the better-affinity prediction.
4. **cuEquivariance kernels not supported on RTX 5070 (Blackwell sm_120)**: per the `boltz-blackwell` v0.1.2 PyPI README, verbatim: "Make sure you're using --no_kernels flag with boltz. The cuEquivariance library doesn't support sm_121 yet. Always use --no_kernels." This means the 2–5× speedup Pierce was counting on from kernels **will not land today**. The Boltz-2 v2.1.1 release notes confirm the kernels deliver "2x speedup and large memory savings" but only on Ampere or newer that cuEquivariance supports (currently sm_80 through sm_90). NVIDIA's cuEquivariance changelog shows ongoing Blackwell support work but no sm_120/121 enablement as of v0.7.0.

## 1.8 Risks & decision gates

| Risk | Probability | Impact | Gate |
|---|---|---|---|
| RTX 5070 not supported by cuEquivariance | **Confirmed** (Blackwell sm_120) | Throughput ~2× lower than planned | Run with `--no_kernels` from the start; accept 4–6 h sweep |
| Affinity head out-of-distribution on Vina poses | Medium | Lose >1.0 pKd accuracy | Validate ACHE/donepezil ±1 log unit gate |
| Allosteric pockets get docked wrong | Medium-High | False negatives on PAMs | Per-target two-pocket strategy; manual visualization of top pose for 5 controls |
| Boltz-2 v2.2.0+ API drift breaks pinned Boltzina | Low | Rebuild | Pin Boltz-2 to v2.1.x to match Boltzina's tested version |

**Keep-vs-abandon gate**: If Boltzina (Cycle=1) yields Spearman ρ < 0.3 vs ChEMBL ground truth on a 20-pair labeled set, **abandon Cycle=1 and run Boltzina default (5 cycles, 7.3× speedup)**. If full-mode Boltzina yields ρ < 0.3, **abandon Boltzina and use Boltz-2 default with the WSL2 environment** (12–20h sweep, but it's the V2 plan you already validated).

---

# TOPIC 2 — MAMMAL Fine-tune on Cognition-DTI Corpus

## 2.1 MAMMAL fine-tuning surface

**Repo**: `BiomedSciAI/biomed-multi-alignment` (Apache 2.0, 44 stars / 19 forks as of May 2026 per the repo's GitHub Actions workflow page). The fine-tune entry point is documented:

```bash
python mammal/main_finetune.py \
  --config-name config.yaml \
  --config-path examples/dti_bindingdb_kd
```

**Architecture** (verified via `Kymi808/mammal-lora-bbbp` SUMMARY.md):
- Inner: `T5ForConditionalGeneration` — 351M parameters, 36 attention blocks split across encoder + decoder
- Head: `ClassifierMLP` — 106M parameters (frozen for BBBP; possibly trainable for DTI regression)
- Wrapper: `mammal.model.Mammal` custom class (not a standard HF AutoModel — verbatim from SUMMARY.md: "The model is not a standard HuggingFace AutoModel — it requires IBM's biomed-multi-alignment framework, which exposes a custom Mammal class wrapping an inner T5ForConditionalGeneration (351 M parameters) plus a ClassifierMLP head (106 M, frozen for BBBP) and scalar projection heads.")
- Tokenizer: `ModularTokenizerOp` from `fuse-med-ml` (combines AA + SMILES + gene-expression vocabs into a unified ID space — repo path `BiomedSciAI/fuse-med-ml/tree/master/fuse/data/tokenizers/modular_tokenizer`)

**Scalar DTI prediction path** (verbatim from the HuggingFace model card for `ibm-research/biomed.omics.bl.sm.ma-ted-458m.dti_bindingdb_pkd`):
```python
# forward pass - encoder_only mode which supports scalar predictions
batch_dict = model.forward_encoder_only(sample_dict)
```

**Loss** (from the bundled BindingDB_Kd example): MSE on normalized pKd (mean/std normalization computed on train split, recorded in config.yaml, applied at inference; per the README data are "harmonize[d] using `data.harmonize_affinities(mode='max_affinity')` and transform[ed] to log-scale").

## 2.2 LoRA strategy

**Target modules**: For T5, the inner attention projections are named `q`, `k`, `v`, `o` (HuggingFace T5 convention). Inspect via:
```python
for name, mod in model.named_modules():
    if any(name.endswith(s) for s in ['.q', '.k', '.v', '.o', '.wi', '.wo']):
        print(name, type(mod).__name__)
```

**Recommended LoRA config** (`src/mammal_repurposing/finetune/lora_config.py`):
```python
from peft import LoraConfig, TaskType

LORA_CFG = LoraConfig(
    r=16,                              # rank — 8 is min for 458M, 32 if data > 50k pairs
    lora_alpha=32,                     # 2*r standard
    target_modules=["q", "k", "v", "o"],   # full attention, both enc & dec
    lora_dropout=0.10,                 # high to combat overfitting on small corpus
    bias="none",
    task_type=TaskType.SEQ_2_SEQ_LM,
    modules_to_save=["scalar_head"],   # train the scalar prediction MLP too
)
```

**Trainable parameter count estimate at r=16, QKVO on all 36 blocks**:
- Each LoRA pair: r × d_model × 2 = 16 × 1024 × 2 = 32,768 params per projection
- 36 blocks × 4 projections × 32,768 ≈ **4.7M params** (attention only)
- Plus scalar head ClassifierMLP if trained: freeze body, train final layer only → ~1M
- **Total trainable: ~5–6M of 458M (≈ 1.2%)**

**VRAM footprint on RTX 5070 (12 GB)**:
- bf16 weights: 458M × 2 B = 0.92 GB
- bf16 activations at seq_len 1024 (target_seq up to ~700 AA + drug SMILES ~200 + special tokens), batch 8: ~3 GB with gradient checkpointing enabled
- LoRA params + Adam optimizer state (m, v in fp32): 6M × 8 B = 48 MB
- Pad for fragmentation: 2 GB
- **Total ~6–8 GB** — fits comfortably; you can push to batch 16 if seq is short

## 2.3 Cognition-DTI training corpus assembly

### 2.3.1 Verified ChEMBL target mapping (CRITICAL — 9 of 22 needed correction)

**Critical correction**: Independent ChEMBL search (cross-referenced via UniProt) revealed that 9 of Pierce's 22 IDs are wrong. The corrected mapping:

| Target | UniProt | **Corrected** ChEMBL ID | Pierce's original | Issue |
|---|---|---|---|---|
| CHRNA7 | P36544 | CHEMBL2492 ✓ | CHEMBL2492 | OK |
| ACHE | P22310 | CHEMBL220 ✓ | CHEMBL220 | OK |
| GRIA1 | P42261 | **CHEMBL2009** | CHEMBL3503 | CHEMBL3503 is rat GRIA2 |
| GRIA2 | P42262 | **CHEMBL4016** | CHEMBL3504 | wrong protein |
| GRIA3 | P42263 | **CHEMBL3595** | CHEMBL2243 | **CHEMBL2243 = FAAH** |
| GRIA4 | P48058 | **CHEMBL3190** | CHEMBL3505 | wrong protein |
| GRIN2A | Q12879 | **CHEMBL1972** | CHEMBL1923 | wrong protein |
| GRIN2B | Q13224 | **CHEMBL1904** | CHEMBL1968 | **CHEMBL1968 = EPHX1** |
| DRD1 | P21728 | CHEMBL2056 ✓ | CHEMBL2056 | OK |
| SLC6A3 | Q01959 | CHEMBL238 ✓ | CHEMBL238 | OK |
| ADRA2A | P08913 | CHEMBL1867 ✓ | CHEMBL1867 | OK |
| SLC6A2 | P23975 | CHEMBL222 ✓ | CHEMBL222 | OK |
| HRH3 | Q9Y5N1 | CHEMBL264 ✓ | CHEMBL264 | OK |
| HCRTR1 | O43613 | CHEMBL5113 ✓ | CHEMBL5113 | OK |
| HCRTR2 | O43614 | CHEMBL4792 ✓ | CHEMBL4792 | OK |
| PDE4D | Q08499 | CHEMBL288 ✓ | CHEMBL288 | OK |
| PDE9A | O76083 | CHEMBL5292 ✓ | CHEMBL5292 | OK |
| NTRK2 | Q16620 | CHEMBL4898 ✓ | CHEMBL4898 | OK |
| SIGMAR1 | Q99720 | CHEMBL287 ✓ | CHEMBL287 | OK |
| KCNQ2 | O43526 | **CHEMBL2476** | CHEMBL4051 | **CHEMBL4051 = CFTR** |
| KCNQ3 | O43525 | **CHEMBL3819** (likely) | CHEMBL3231 | **CHEMBL3231 = ROCK1** |
| HCN1 | O60741 | CHEMBL1971 (verify!) | CHEMBL5374 | unverified |

**Action**: before any ChEMBL pull, run `chembl_id_lookup` against each UniProt accession via the ChEMBL REST API to lock in IDs. Otherwise the fine-tune will silently train on the wrong targets and the CHRNA7 audit will look correct while GRIN2B is invisibly trained on EPHX1 inhibitors.

### 2.3.2 ChEMBL extraction recipe

```python
# src/mammal_repurposing/finetune/data.py
import pandas as pd
from chembl_webresource_client.new_client import new_client

VERIFIED_IDS = {
    "CHRNA7": "CHEMBL2492", "ACHE": "CHEMBL220",
    "GRIA1": "CHEMBL2009", "GRIA2": "CHEMBL4016",
    "GRIA3": "CHEMBL3595", "GRIA4": "CHEMBL3190",
    "GRIN2A": "CHEMBL1972", "GRIN2B": "CHEMBL1904",
    "DRD1": "CHEMBL2056", "SLC6A3": "CHEMBL238",
    "ADRA2A": "CHEMBL1867", "SLC6A2": "CHEMBL222",
    "HRH3": "CHEMBL264", "HCRTR1": "CHEMBL5113",
    "HCRTR2": "CHEMBL4792", "PDE4D": "CHEMBL288",
    "PDE9A": "CHEMBL5292", "NTRK2": "CHEMBL4898",
    "SIGMAR1": "CHEMBL287", "KCNQ2": "CHEMBL2476",
    "KCNQ3": "CHEMBL3819", "HCN1": "CHEMBL1971",  # verify via UniProt!
}

def pull_target(tid: str) -> pd.DataFrame:
    act = new_client.activity.filter(
        target_chembl_id=tid,
        assay_type="B",
        standard_type__in=["Ki", "IC50", "Kd", "EC50"],
        standard_value__isnull=False,
        target_organism="Homo sapiens",
    ).only([
        "molecule_chembl_id", "canonical_smiles",
        "standard_type", "standard_value", "standard_units",
        "assay_chembl_id", "assay_description", "confidence_score",
        "pchembl_value", "target_chembl_id",
    ])
    df = pd.DataFrame(list(act))
    df = df[df["confidence_score"].astype(float) >= 7]
    df = df[df["standard_units"].isin(["nM", "uM"])]
    df["pchembl_value"] = pd.to_numeric(df["pchembl_value"], errors="coerce")
    df = df.dropna(subset=["pchembl_value", "canonical_smiles"])
    return df

def harmonize_to_pkd(df: pd.DataFrame) -> pd.DataFrame:
    # Already in pchembl_value column; deduplicate (target, smiles) by max affinity
    df["target_smiles"] = df["target_chembl_id"] + "|" + df["canonical_smiles"]
    df = df.sort_values("pchembl_value", ascending=False)
    df = df.drop_duplicates("target_smiles", keep="first")
    return df

def annotate_allosteric(df: pd.DataFrame) -> pd.DataFrame:
    # Following the Burggraaff et al. 2020 J Chem Inf Model 60:7 approach
    # (text-mining ChEMBL assay descriptions for binding-type annotation)
    allosteric_keywords = ["allosteric", "PAM", "NAM", "positive allosteric",
                            "negative allosteric", "non-competitive", "modulator"]
    pat = "|".join(allosteric_keywords)
    df["is_allosteric"] = df["assay_description"].fillna("").str.contains(pat, case=False)
    return df
```

### 2.3.3 Expected corpus volume

Per the ChEMBL inventory (subagent investigation, ±30% accuracy until SQL'd), after the standard DTI filter (assay_type='B', std_type∈{Ki,IC50,Kd,EC50}, std_value not null, confidence_score≥7, human only):

- **Raw records**: ~50,000–70,000
- **Unique (target, compound) pairs after dedup**: **~35,000–45,000**
- **Targets clearing ≥100 high-confidence pairs**: **18 of 22** (GRIA3 [~50–150], GRIA4 [~50–150], GRIN2A [~400–800 — borderline OK], HCN1 [~150–400] are data-sparse)
- **Pairs explicitly annotated allosteric**: ~2,000–4,000 (text-mining via assay_description; sparse but enriched on CHRNA7, GRIN2B, PDE4D)

**Top-5 contributors by volume**: ACHE (~7k records / ~5.8k compounds — paper-style confirmation: AChE has 8205 IC50 records and 7026 pChEMBL records as of ChEMBL30 per Boulaamane's tutorial), HRH3 (~6–8k), SIGMAR1 (~3.5–5.5k), SLC6A3 (~5–7k), CHRNA7 (~3.5–5k).

For the data-sparse subtypes (GRIA3, GRIA4, HCN1), augment with the corresponding rat orthologue ChEMBL targets (CHEMBL3503 rat GluA2, etc.) — accept the species mismatch as data augmentation.

### 2.3.4 Train/val/test split

**Use scaffold split** (Bemis-Murcko) for the test set, **random split for val**. This catches the catastrophic case where the LoRA learns to recognize the donepezil scaffold (which has thousands of analogues in ChEMBL) rather than the ACHE binding pattern.

```python
from rdkit.Chem.Scaffolds import MurckoScaffold
def scaffold(smi): return MurckoScaffold.MurckoScaffoldSmilesFromSmiles(smi)
# Group by scaffold; assign whole scaffold groups to train/val/test 80/10/10
```

**Held-out for the "rescue" test (NEVER in train)**:
- CHRNA7: galantamine, encenicline, TC-5619
- ACHE: donepezil, rivastigmine, huperzine A
- SLC6A3: methylphenidate
- SLC6A2: atomoxetine
- PDE4D: BPN14770, rolipram
- HRH3: pitolisant

### 2.3.5 Negative set strategy

ChEMBL true inactives (`pchembl_value < 5`) are limited. Augment with DUD-E decoys:
```python
# For each target: take 50 actives → fetch 1,500 DUD-E decoys (property-matched, low Tanimoto)
# Label decoys with synthetic pchembl_value = 3.0 (1 mM, below assay limit)
```

For BindingDB-trained MAMMAL, training on negatives in this range is well within distribution. Boltz-2 used a similar synthetic-decoy strategy with a "Tanimoto similarity below 0.3 to all known binders associated with similar proteins" cutoff (Passaro et al. 2025 §affinity training data).

## 2.4 Fine-tune implementation plan

### File layout
```
src/mammal_repurposing/finetune/
├── __init__.py
├── data.py            # ChEMBL pull, scaffold split, MAMMAL tokenization
├── lora_config.py     # LORA_CFG (see §2.2)
├── train.py           # Lightning trainer + PEFT
├── eval.py            # CHRNA7-rescue / PDE4D / orthosteric-no-regression tests
└── infer.py           # Load LoRA adapter, score (target, drug) pair
```

### Hyperparameters (concrete)

```python
# train.py
HPARAMS = {
    "lr_peak": 5e-4,              # LoRA tolerates 10× higher LR than full FT
    "lr_warmup_steps": 200,
    "lr_schedule": "cosine",
    "weight_decay": 0.01,
    "epochs": 10,
    "batch_size": 8,              # bf16, RTX 5070, gradient checkpointing
    "grad_accum_steps": 4,        # effective batch 32
    "max_grad_norm": 1.0,
    "early_stopping_patience": 3, # on val MSE
    "bf16": True,
    "gradient_checkpointing": True,
    "seed": 42,
    "save_top_k": 3,
}
```

These follow the published FLAN-T5 PEFT recipe (LoRA r=16–32, q/v target modules, lr=1e-3 with `auto_find_batch_size`) but with reduced rank for the 458M backbone and 10-epoch cosine schedule for the regression target.

### Wall-clock estimate

~10,000 effective training pairs × 10 epochs / batch 8 = 12,500 steps. T5-encoder forward + backward at batch 8, seq 1024, bf16 on RTX 5070: ~0.8–1.2 s/step. **Total train: 3–4 hours**. Plus validation: ~30 min. **Overnight run.**

### Inference path (preserves V2 pipeline)

```python
# infer.py — drop-in replacement for current MAMMAL DTI scorer
from peft import PeftModel
from mammal.model import Mammal

class CognitionDtiScorer:
    def __init__(self, base_id="ibm/biomed.omics.bl.sm.ma-ted-458m.dti_bindingdb_pkd",
                 lora_path="checkpoints/lora_cognition_v1"):
        self.base = Mammal.from_pretrained(base_id)
        self.model = PeftModel.from_pretrained(self.base, lora_path, is_trainable=False)
        self.model.eval()
        # ... tokenizer load same as base MAMMAL
```

The V2 pipeline switches scorer with one line — old behavior available by `lora_path=None`.

## 2.5 Evaluation strategy

### Primary metric (rescue test)
**CHRNA7 PAM top-25% rescue**: Run the fine-tuned head on the 298-compound library + CHRNA7 sequence. Verify galantamine, encenicline, TC-5619 (all in held-out test set) rank in top 75 of 298. **Pass = 2 of 3 in top 25%.**

### Secondary metric (dynamic range)
Per-target std of predicted pKd across the 298-compound library:
- Pre-fine-tune CHRNA7: 0.029 (measured)
- Post-fine-tune CHRNA7: target **>0.30** (10× improvement)
- AMPA bunching: target per-subtype std > 0.20

### Tertiary metric (PDE4D NAM preservation)
BPN14770 must remain in top 25% (it currently is, per V2 metrics). Rolipram in top 50%.

### Negative regression test (must not break)
- Donepezil at ACHE: stays in top 10
- Methylphenidate at SLC6A3: stays in top 10
- Atomoxetine at SLC6A2: stays in top 10
- Pitolisant at HRH3: stays in top 10

### Catastrophic forgetting test
Hold out 500 random TDC BindingDB_Kd pairs (non-cognition). Pre/post-LoRA RMSE delta should be < 0.3 log units. If > 0.5, LoRA rank is too high or LR too aggressive.

## 2.6 Alternatives compared

| Approach | Pros | Cons | Verdict |
|---|---|---|---|
| **LoRA r=16 QKVO** | Fits 12GB, ~5M params, 3-4h train | Bounded capacity | **Primary** |
| Full fine-tune | Maximum capacity | OOMs on 12GB without CPU offload | Not feasible |
| Prompt tuning | <1M params | Underfits for regression with high inter-target variance | Skip |
| LoRA r=32 | More capacity for 45k pairs | 2× memory; overfit risk on sparse targets | **Fallback** if r=16 underfits |
| Boltz-2 soft-label distillation | Pocket-aware signal | Requires running Boltz-2 on ~10k pairs first (60h) | Phase 0.5 option |
| Train new 50M head from scratch | Pierce controls full stack | 45k pairs is too small for 50M params; loses MAMMAL pretraining | Skip |

## 2.7 Risks & decision gates

| Risk | Mitigation | Gate |
|---|---|---|
| Overfitting on small corpus (~10k effective) | Scaffold split, LoRA dropout 0.10, early stopping | Val MSE plateaus epoch 3-5 → stop |
| Catastrophic forgetting non-cognition | Held-out BindingDB test (500 pairs) | ΔRMSE > 0.5 → reduce LR to 1e-4 |
| Allosteric ligand sparsity (<4k explicit) | Augment with text-mined annotations + DUD-E decoys | If allosteric subset < 1.5k, expect modest rescue only |
| ChEMBL label noise (multiple Ki per pair) | `harmonize_affinities(mode='max_affinity')` per the MAMMAL DTI README | n/a |
| Wrong ChEMBL IDs (9 of 22 invalid) | Per §2.3.1 corrections; lookup via UniProt | Verify on day 1 |
| Fine-tune sharpens but doesn't add pocket reasoning | Pre-commit to Boltzina backup for allosteric targets | If CHRNA7 rescue fails, route to Boltzina |

**Keep-vs-abandon gate**: Run a quick 1-epoch sanity check first. If after 1 epoch (~25 min) the val MSE is not below pre-train val MSE by ≥ 10%, the LoRA target modules or data tokenization is broken — debug before the full 10-epoch run.

**Final keep gate (after full train)**:
- CHRNA7 rescue: 2 of 3 PAMs in top 25% → keep
- Orthosteric controls all in top 10 → keep
- Non-cognition TDC BindingDB_Kd ΔRMSE < 0.5 → keep
- All three pass → ship to V3 pipeline
- Any one fails → debug; do not ship

---

## Recommendations

### Sequencing (Pierce, week-by-week)

**Week 1: Boltzina end-to-end in WSL2 (high priority, low risk)**
1. **Day 1–2**: Vendor `ohuelab/boltzina` as git submodule pinned to a specific commit SHA. Install in WSL2 Ubuntu (Python 3.12, CUDA 13). Boltz-2 pinned to v2.1.x (last release before v2.2.0 affinity API drift). **Critical**: install with `--no_kernels` from the start because RTX 5070 (sm_120) is not supported by cuEquivariance (per boltz-blackwell v0.1.2 PyPI README). Run the bundled CDK2 sample to verify install. Benchmark 10 random pairs to measure actual per-pair time on RTX 5070. **Gate: if > 30 s/pair, profile what's slow — Vina or Boltz-2 affinity — before proceeding.**
2. **Day 3**: Build the pocket database (`mammal_repurposing/data/pockets.py`) for all 22 targets. Compute centroids from co-crystal PDBs. For HRH3 / GRIA3 / GRIA4 (no human cryo-EM), use AlphaFold2 predicted structures + homology to GluA2 / HRH4.
3. **Day 4**: Run validation suite (donepezil/ACHE, BPN14770/PDE4D, PNU-120596/CHRNA7, pitolisant/HRH3, methylphenidate/SLC6A3). **Gate: all 5 must predict within ±1.0 log10(μM) of experimental.**
4. **Day 5**: Phase 0.4 sweep — 1,500 pairs in Top-5-Average mode at Cycle=1. Estimate 4–6 hours overnight on RTX 5070 + WSL2 (no kernels).

**Week 2: ChEMBL pull + fine-tune corpus**
1. **Day 6**: Verify all 22 ChEMBL IDs via UniProt → ChEMBL `chembl_id_lookup` REST endpoint. Pull all targets via `chembl_webresource_client`. Harmonize, scaffold-split, save as `corpus_v1.parquet`. **Gate: ≥ 8,000 unique (target, compound) pairs in train set.**
2. **Day 7–8**: Build `finetune/` module. Run 1-epoch sanity test.

**Week 3: Full fine-tune + evaluation**
1. **Day 9–10**: 10-epoch overnight train. Eval on CHRNA7-rescue, PDE4D, orthosteric regression, BindingDB catastrophic-forgetting.
2. **Day 11**: Decision gate. If pass, merge LoRA into V3. If fail, plan Boltz-2 distillation as Phase 0.5.

### Combined V3 architecture
After Topic 1 + Topic 2 land:

```
V3 Pipeline (per (target, compound) pair):
├── MAMMAL DTI head + cognition LoRA       (~5 ms/pair)
├── ADMET-AI                                (unchanged, ~50 ms/compound, cached)
├── Boltzina (Vina + Boltz-2 affinity head) (~10–15 s/pair)
└── RRF fusion                              (per-target re-weighting)
```

For the **8 allosteric-rich targets** (CHRNA7, GRIN2A/2B, PDE4D, KCNQ2/3, GRIA1–4, NTRK2), use per-target RRF weighting favoring Boltzina 2:1 over MAMMAL. For the **12 orthosteric targets** (ACHE, DRD1, SLC6A3, ADRA2A, SLC6A2, HRH3, HCRTR1/2, PDE9A, SIGMAR1), MAMMAL+LoRA is faster and equally accurate — use Boltzina only for top 20% via two-stage screening.

### What would change the recommendation
- If cuEquivariance ships sm_120/121 support → expect 2× Boltzina throughput; rerun Phase 0.4 in 2 hours.
- If ChEMBL pull yields < 8,000 high-confidence pairs (mostly due to GRIA3/4/HCN1 sparsity) → drop LoRA, go straight to Boltz-2 distillation for Phase 0.5.
- If Pierce gets weekend access to A100/H100 (Colab Pro / Lambda) → run full Boltz-2 (no Boltzina) on the 1,500 pairs in ~7 hours and skip the Vina-pose risks entirely; use Boltzina only at scale-up to 10k+ pairs.

---

## Caveats

1. **The Furui & Ohue 2025 paper measures speedup on H100 + 48 CPU cores.** RTX 5070 + WSL2 + consumer CPU is a different regime; expect per-pair times 3–8× slower than the paper's 1.4 s. The relative speedup vs full Boltz-2 (7–11×) should hold; the absolute number won't.
2. **RTX 5070 is Blackwell (sm_120/121).** Per the `boltz-blackwell` v0.1.2 PyPI README, verbatim: "The cuEquivariance library doesn't support sm_121 yet. Always use --no_kernels." This means the 2-5× speedup from CUDA kernels Pierce was counting on **will not materialize today**. Plan all timing on the no-kernel path; the Boltz-2 v2.1.1 release notes confirm the kernels are Ampere-or-newer with cuEquivariance support, which Blackwell sm_120 lacks.
3. **ChEMBL ID verification is non-negotiable.** 9 of 22 IDs supplied were wrong (CHEMBL1968 = EPHX1, CHEMBL4051 = CFTR, CHEMBL2243 = FAAH, CHEMBL3231 = ROCK1, etc.). Skipping the UniProt cross-check would silently corrupt the fine-tune. Add a unit test that asserts `chembl.target.get(tid)['target_components'][0]['accession']` matches the expected UniProt for each panel target.
4. **The Boltzina paper does NOT evaluate absolute affinity prediction or pose validity**, only screening enrichment metrics (AP, EF, ROC-AUC on MF-PCBA). Verbatim from §4 Conclusion: "We did not evaluate absolute affinity prediction or the physical validity of poses, and it remains unverified whether Boltzina can match Boltz-2 on these tasks." For ranking 298 compounds in a hit-discovery context, this is the right metric; for absolute pKd, expect higher RMSE than full Boltz-2.
5. **The fine-tune addresses a symptom, not a cause.** MAMMAL's compressed dynamic range on CHRNA7 reflects that the BindingDB Kd training corpus (52,284 pairs) has limited CHRNA7 PAM coverage and MAMMAL lacks pocket-awareness. Fine-tuning adds CHRNA7-specific SAR signal but does not give the model structural reasoning. For PAMs at *novel* allosteric sites (not represented in ChEMBL), the LoRA will not generalize — only Boltzina will.
6. **The Boltz-2 affinity output is `log10(IC50 in μM)`, not pKd**. Mapping between RRF scores requires either converting Boltzina output to pKd (pKd ≈ -log10(IC50_M) = 6 - log10(IC50_μM)) or re-normalizing the fusion. The existing V2 RRF code probably treats MAMMAL pKd directly — Pierce should verify the conversion before fusion.
7. **Boltz-2's 56-heavy-atom soft limit on the affinity head** means some larger cognition-relevant ligands (e.g., natural product analogues, peptide mimetics) may fall outside distribution. Pre-filter the 298-compound library by heavy atom count and flag oversized ligands as MAMMAL-only.
8. **The cited per-target ChEMBL counts are order-of-magnitude estimates** (±30%); exact figures require a direct SQL query against the ChEMBL release Pierce uses. Run the canonical SQL `SELECT COUNT(DISTINCT activity_id) FROM activities a JOIN assays s USING(assay_id) WHERE s.assay_type='B' AND a.standard_type IN ('Ki','IC50','Kd','EC50') AND s.confidence_score >= 7 AND s.tid IN (...)` against the verified target dictionary to lock in numbers.
9. **Boltzina is research-grade software, not production.** 3 open issues at the time of this report. Pin to a specific commit, vendor under your repo, run the bundled CDK2 sample as a CI regression test. Don't auto-update.
10. **Boltz-2 v2.2.0 changed the affinity API and added molecular-weight correction**. Boltzina was developed against an earlier Boltz-2 version. Verify Boltz-2 version compatibility (the Boltzina pyproject.toml pins a specific Boltz-2 release) before running anything.