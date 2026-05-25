"""Self-contained WSL2-side full Boltz sweep.

Runs inside the WSL2 venv (/root/mammal_env/bin/python). Does NOT import
mammal_repurposing — uses boltz CLI directly via subprocess, reads parquets
via pandas, writes results back to /mnt/c/... so the Windows side picks them up.

Scope (conservative for overnight):
    - ADMET-PASS + ADMET-FLAG compounds from data/results/v2/admet_gates.parquet
    - Top-N targets per compound by MAMMAL pKd (default N=10)
    - Skips pairs already in boltzina_affinity.parquet (resume-safe)
    - Flushes after every pair (cheap, fault-tolerant)

Wall-clock: ~22 MSA cold-fetches × 270s + ~1100 inference × 23s ≈ 10-12 hours.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import pandas as pd
import yaml

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("wsl2_boltz_sweep")

PROJ = Path("/mnt/c/Users/Pierce Lonergan/Documents/GitHub/MAMMAL_Cognitive_Enhancement_Drug_Repurposing")
TARGETS_PARQ = PROJ / "data/interim/targets.parquet"
COMPOUNDS_PARQ = PROJ / "data/interim/compounds.parquet"
ADMET_PARQ = PROJ / "data/results/v2/admet_gates.parquet"
MAMMAL_PARQ = PROJ / "data/results/dti_scores.parquet"
OUT_PARQ = PROJ / "data/results/v2/boltzina_affinity.parquet"
CACHE_DIR = PROJ / "data/cache/boltzina"


def _find_boltz() -> str:
    """Locate boltz CLI; prefer venv-relative."""
    venv_bin = Path(sys.executable).parent / "boltz"
    if venv_bin.exists():
        return str(venv_bin)
    found = shutil.which("boltz")
    if found:
        return found
    raise FileNotFoundError("boltz CLI not on PATH and not at " + str(venv_bin))


def _pair_hash(seq: str, smiles: str) -> str:
    h = hashlib.sha1()
    h.update(seq.encode())
    h.update(b"||")
    h.update(smiles.encode())
    return h.hexdigest()


def _find_outputs(out_dir: Path):
    aff = next((p for p in out_dir.rglob("*.json") if "affinity" in p.name.lower()), None)
    conf = next((p for p in out_dir.rglob("*.json") if "confidence" in p.name.lower()), None)
    return aff, conf


def _score_pair(boltz_exe: str, target_uniprot: str, seq: str,
                compound_name: str, smiles: str,
                use_cache: bool = True) -> dict:
    h = _pair_hash(seq, smiles)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_path = CACHE_DIR / f"{h}.json"
    if use_cache and cache_path.exists():
        with open(cache_path) as f:
            return json.load(f)

    with tempfile.TemporaryDirectory(prefix="boltz_") as tmp:
        work = Path(tmp)
        run = f"aff_{h[:10]}"
        yaml_path = work / f"{run}.yaml"
        with open(yaml_path, "w") as f:
            yaml.safe_dump({
                "version": 1,
                "sequences": [
                    {"protein": {"id": "A", "sequence": seq}},
                    {"ligand": {"id": "L", "smiles": smiles}},
                ],
                "properties": [{"affinity": {"binder": "L"}}],
            }, f, sort_keys=False)
        out_dir = work / "out"
        out_dir.mkdir()

        cmd = [boltz_exe, "predict", str(yaml_path),
               "--out_dir", str(out_dir),
               "--use_msa_server",
               "--recycling_steps", "3",
               "--diffusion_samples", "1",
               "--output_format", "mmcif",
               "--accelerator", "gpu",
               "--devices", "1"]

        env = os.environ.copy()
        env.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")

        proc = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=1800)
        if proc.returncode != 0:
            logger.warning("boltz exit %d for (%s,%s):\n%s",
                           proc.returncode, target_uniprot, compound_name,
                           proc.stderr[-400:])
            raise RuntimeError(f"boltz failed exit {proc.returncode}")

        aff, conf = _find_outputs(out_dir)
        if aff is None:
            raise RuntimeError("no affinity JSON")
        with open(aff) as f:
            aff_json = json.load(f)
        pose_plddt = None
        if conf is not None:
            with open(conf) as f:
                cj = json.load(f)
            pr = cj.get("plddt") or cj.get("per_residue_plddt") or []
            if pr:
                pose_plddt = float(sum(pr) / len(pr))

    result = {
        "target_uniprot": target_uniprot,
        "compound_name": compound_name,
        "smiles": smiles,
        "affinity_pred_value": float(aff_json.get("affinity_pred_value", float("nan"))),
        "affinity_probability_binary": float(aff_json.get("affinity_probability_binary", float("nan"))),
        "pose_plddt": pose_plddt,
        "mode": "wsl2_full",
        "scored_at": dt.datetime.now(dt.timezone.utc).isoformat(),
    }
    with open(cache_path, "w") as f:
        json.dump(result, f)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--top-n", type=int, default=10,
                        help="Per surviving compound, score top-N MAMMAL targets")
    parser.add_argument("--limit", type=int, default=None,
                        help="Smoke: process only first N pairs")
    parser.add_argument("--include-flag", action="store_true", default=True)
    args = parser.parse_args()

    targets = pd.read_parquet(TARGETS_PARQ)
    compounds = pd.read_parquet(COMPOUNDS_PARQ)
    gates = pd.read_parquet(ADMET_PARQ)
    mammal = pd.read_parquet(MAMMAL_PARQ)
    logger.info("Loaded: %d targets, %d compounds, %d gated, %d mammal pairs",
                len(targets), len(compounds), len(gates), len(mammal))

    allowed = {"PASS"} | ({"FLAG"} if args.include_flag else set())
    surviving = gates[gates["gate_status"].isin(allowed)]
    surviving_lc = set(surviving["compound_name"].str.lower().str.strip())
    logger.info("Surviving ADMET: %d compounds", len(surviving_lc))

    tgt_seq = targets.set_index("uniprot")["sequence"].to_dict()
    cmp_smiles = {n.lower().strip(): (s, n)
                  for n, s in zip(compounds["name"], compounds["smiles"])}

    # Build top-N pair list
    mm = mammal[mammal["compound_name"].str.lower().str.strip().isin(surviving_lc)].copy()
    mm["lc"] = mm["compound_name"].str.lower().str.strip()
    topn = (mm.sort_values(["lc", "predicted_pkd"], ascending=[True, False])
              .groupby("lc", group_keys=False)
              .head(args.top_n))
    pairs = []
    for _, row in topn.iterrows():
        tgt = row["target_uniprot"]
        seq = tgt_seq.get(tgt)
        lc = row["compound_name"].lower().strip()
        if lc not in cmp_smiles or seq is None:
            continue
        smi, orig = cmp_smiles[lc]
        pairs.append((tgt, seq, orig, smi))

    # Resume: skip pairs already in OUT_PARQ
    existing = pd.read_parquet(OUT_PARQ) if OUT_PARQ.exists() else pd.DataFrame()
    done = set()
    if not existing.empty:
        done = set(zip(existing["target_uniprot"].astype(str),
                       existing["compound_name"].astype(str).str.lower().str.strip()))
    remaining = [p for p in pairs if (p[0], p[2].lower().strip()) not in done]
    logger.info("Pairs: built=%d, already_done=%d, remaining=%d",
                len(pairs), len(pairs) - len(remaining), len(remaining))

    if args.limit:
        remaining = remaining[:args.limit]
        logger.info("Limit: %d pairs", len(remaining))

    boltz_exe = _find_boltz()
    logger.info("Using boltz: %s", boltz_exe)

    completed = [existing] if not existing.empty else []
    new_rows = []
    last_target = None
    t0 = time.perf_counter()
    for i, (tgt, seq, name, smi) in enumerate(remaining):
        pair_start = time.perf_counter()
        try:
            r = _score_pair(boltz_exe, tgt, seq, name, smi)
            new_rows.append(r)
            ok = True
        except Exception as e:
            logger.warning("Pair (%s, %s) failed: %s; NaN", tgt, name, e)
            new_rows.append({
                "target_uniprot": tgt, "compound_name": name, "smiles": smi,
                "affinity_pred_value": float("nan"),
                "affinity_probability_binary": float("nan"),
                "pose_plddt": None, "mode": "wsl2_full",
                "scored_at": dt.datetime.now(dt.timezone.utc).isoformat(),
            })
            ok = False
        elapsed = time.perf_counter() - pair_start
        avg = (time.perf_counter() - t0) / (i + 1)
        eta = avg * (len(remaining) - i - 1) / 60
        logger.info("[%d/%d] %s + %s -> %.1fs (avg %.1fs, ETA %.0f min) ok=%s",
                    i + 1, len(remaining), tgt, name, elapsed, avg, eta, ok)
        last_target = tgt

        # Flush every pair
        df_new = pd.DataFrame(new_rows)
        combined = pd.concat(completed + [df_new], ignore_index=True, sort=False)
        tmp = OUT_PARQ.with_suffix(OUT_PARQ.suffix + ".tmp")
        combined.to_parquet(tmp, index=False)
        tmp.replace(OUT_PARQ)

    logger.info("DONE. Total rows in %s: %d", OUT_PARQ,
                len(pd.read_parquet(OUT_PARQ)) if OUT_PARQ.exists() else 0)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
