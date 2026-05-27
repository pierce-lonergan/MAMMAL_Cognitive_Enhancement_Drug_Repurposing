"""§7.17 — Pose-saving variant of `_wsl2_boltz_full_sweep.py`.

Same sweep loop as the original, but extends each parquet row with three
new columns derived from the Boltz mmCIF pose output:
    pose_centroid_x, pose_centroid_y, pose_centroid_z

These let downstream code call `pockets.classify_pose(pose_xyz, target_gene)`
to operationalise §7.5 on the full grid + §8.13 pocket-conditioned liability.

Also (optionally) saves the mmCIF pose itself to
    data/results/v2/boltz_poses/<target>_<compound_safe>.cif
when --save-poses is passed (storage ~50 KB/pair × 1165 = ~60 MB).

Operation:
    bash on WSL2:
      /root/mammal_env/bin/python \\
        /mnt/c/.../_wsl2_boltz_full_sweep_pose.py \\
        --top-n 10

For the pose-only re-run (the existing 926 pairs already have affinity in
boltzina_affinity.parquet — this script reads centroids and BACK-FILLS the
three pose columns on the existing rows without recomputing affinity):

    /root/mammal_env/bin/python \\
        /mnt/c/.../_wsl2_boltz_full_sweep_pose.py --backfill-only
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

# Make the package importable from the WSL2 mount
PROJ = Path("/mnt/c/Users/Pierce Lonergan/Documents/GitHub/MAMMAL_Cognitive_Enhancement_Drug_Repurposing")
if str(PROJ / "src") not in sys.path:
    sys.path.insert(0, str(PROJ / "src"))

from mammal_repurposing.pockets.pose_extract import (  # noqa: E402
    extract_centroid_from_file,
)

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("wsl2_boltz_sweep_pose")

TARGETS_PARQ = PROJ / "data/interim/targets.parquet"
COMPOUNDS_PARQ = PROJ / "data/interim/compounds.parquet"
ADMET_PARQ = PROJ / "data/results/v2/admet_gates.parquet"
MAMMAL_PARQ = PROJ / "data/results/dti_scores.parquet"
OUT_PARQ = PROJ / "data/results/v2/boltzina_affinity.parquet"
CACHE_DIR = PROJ / "data/cache/boltzina"
POSE_DIR = PROJ / "data/results/v2/boltz_poses"


def _safe_compound_filename(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_-]+", "_", name)[:60]


def _find_boltz() -> str:
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


def _find_outputs(out_dir: Path) -> tuple[Path | None, Path | None, Path | None]:
    """Locate (affinity_json, confidence_json, pose_cif)."""
    aff = next((p for p in out_dir.rglob("*.json") if "affinity" in p.name.lower()), None)
    conf = next((p for p in out_dir.rglob("*.json") if "confidence" in p.name.lower()), None)
    cif = next(out_dir.rglob("*.cif"), None)
    return aff, conf, cif


def _score_pair(
    boltz_exe: str,
    target_uniprot: str,
    seq: str,
    compound_name: str,
    smiles: str,
    save_pose: bool = False,
    use_cache: bool = True,
) -> dict:
    """Score one pair and additionally extract pose centroid (+ optional pose save)."""
    h = _pair_hash(seq, smiles)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_path = CACHE_DIR / f"{h}.json"
    pose_cache_path = CACHE_DIR / f"{h}_pose.json"

    cached_pose = None
    if use_cache and pose_cache_path.exists():
        try:
            cached_pose = json.loads(pose_cache_path.read_text())
        except Exception:
            cached_pose = None

    if use_cache and cache_path.exists():
        with open(cache_path) as f:
            r = json.load(f)
        # If we have a cached pose for this same hash, merge it in
        if cached_pose:
            r.update(cached_pose)
        return r

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

        aff, conf, cif = _find_outputs(out_dir)
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

        # --- §7.17 — pose centroid extraction ----------------------------
        pose_centroid = None
        pose_path_str = None
        if cif is not None:
            try:
                c = extract_centroid_from_file(cif)
                if c is not None:
                    pose_centroid = (float(c[0]), float(c[1]), float(c[2]))
            except Exception as e:
                logger.warning("Centroid extract failed for (%s,%s): %s",
                               target_uniprot, compound_name, e)

            if save_pose:
                POSE_DIR.mkdir(parents=True, exist_ok=True)
                safe = _safe_compound_filename(compound_name)
                dst = POSE_DIR / f"{target_uniprot}_{safe}.cif"
                shutil.copy(cif, dst)
                pose_path_str = str(dst.relative_to(PROJ))

    result = {
        "target_uniprot": target_uniprot,
        "compound_name": compound_name,
        "smiles": smiles,
        "affinity_pred_value": float(aff_json.get("affinity_pred_value", float("nan"))),
        "affinity_probability_binary": float(aff_json.get("affinity_probability_binary", float("nan"))),
        "pose_plddt": pose_plddt,
        "pose_centroid_x": pose_centroid[0] if pose_centroid else None,
        "pose_centroid_y": pose_centroid[1] if pose_centroid else None,
        "pose_centroid_z": pose_centroid[2] if pose_centroid else None,
        "pose_cif_path": pose_path_str,
        "mode": "wsl2_full_pose",
        "scored_at": dt.datetime.now(dt.timezone.utc).isoformat(),
    }
    # Cache both affinity (compatible with original) and pose extras separately
    aff_only = {k: v for k, v in result.items()
                if k not in ("pose_centroid_x", "pose_centroid_y", "pose_centroid_z",
                             "pose_cif_path")}
    with open(cache_path, "w") as f:
        json.dump(aff_only, f)
    pose_only = {k: result.get(k) for k in
                 ("pose_centroid_x", "pose_centroid_y", "pose_centroid_z",
                  "pose_cif_path")}
    with open(pose_cache_path, "w") as f:
        json.dump(pose_only, f)
    return result


def _backfill_pose_columns(existing_df: pd.DataFrame) -> pd.DataFrame:
    """For each row in existing_df with no pose centroid, try to read from
    cache (set by an earlier run that didn't save mmCIF, or from a pose-only
    re-run). Adds the three pose columns if missing.
    """
    out = existing_df.copy()
    for col in ("pose_centroid_x", "pose_centroid_y", "pose_centroid_z",
                "pose_cif_path"):
        if col not in out.columns:
            out[col] = None

    missing = out["pose_centroid_x"].isna()
    n_missing = int(missing.sum())
    if n_missing == 0:
        logger.info("All %d rows already have pose centroids", len(out))
        return out
    logger.info("Backfilling pose centroids for %d / %d rows", n_missing, len(out))

    filled = 0
    for i, row in out[missing].iterrows():
        h = _pair_hash(row.get("_seq", ""), row.get("smiles", ""))
        pose_cache_path = CACHE_DIR / f"{h}_pose.json"
        if pose_cache_path.exists():
            try:
                pose = json.loads(pose_cache_path.read_text())
                out.at[i, "pose_centroid_x"] = pose.get("pose_centroid_x")
                out.at[i, "pose_centroid_y"] = pose.get("pose_centroid_y")
                out.at[i, "pose_centroid_z"] = pose.get("pose_centroid_z")
                out.at[i, "pose_cif_path"] = pose.get("pose_cif_path")
                filled += 1
            except Exception as e:
                logger.warning("Backfill failed for row %d: %s", i, e)
    logger.info("Backfilled %d / %d", filled, n_missing)
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--top-n", type=int, default=10)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--include-flag", action="store_true", default=True)
    parser.add_argument("--save-poses", action="store_true",
                        help="Also copy the mmCIF pose to data/results/v2/boltz_poses/.")
    parser.add_argument("--backfill-only", action="store_true",
                        help="Don't re-run boltz; just look in the pair-hash "
                             "cache for any pose-extract results and stamp them "
                             "onto the existing boltzina_affinity.parquet.")
    args = parser.parse_args()

    targets = pd.read_parquet(TARGETS_PARQ)
    compounds = pd.read_parquet(COMPOUNDS_PARQ)
    gates = pd.read_parquet(ADMET_PARQ)
    mammal = pd.read_parquet(MAMMAL_PARQ)
    logger.info("Loaded: %d targets, %d compounds, %d gated, %d mammal pairs",
                len(targets), len(compounds), len(gates), len(mammal))

    if args.backfill_only:
        if not OUT_PARQ.exists():
            logger.error("No %s to backfill", OUT_PARQ)
            return 1
        existing = pd.read_parquet(OUT_PARQ)
        # Need the sequence to compute the pair hash — join from targets
        tgt_seq = targets.set_index("uniprot")["sequence"].to_dict()
        existing["_seq"] = existing["target_uniprot"].map(tgt_seq)
        augmented = _backfill_pose_columns(existing)
        augmented = augmented.drop(columns=["_seq"])
        tmp = OUT_PARQ.with_suffix(OUT_PARQ.suffix + ".tmp")
        augmented.to_parquet(tmp, index=False)
        tmp.replace(OUT_PARQ)
        logger.info("Wrote backfilled %s (%d rows)", OUT_PARQ, len(augmented))
        return 0

    allowed = {"PASS"} | ({"FLAG"} if args.include_flag else set())
    surviving = gates[gates["gate_status"].isin(allowed)]
    surviving_lc = set(surviving["compound_name"].str.lower().str.strip())
    logger.info("Surviving ADMET: %d compounds", len(surviving_lc))

    tgt_seq = targets.set_index("uniprot")["sequence"].to_dict()
    cmp_smiles = {n.lower().strip(): (s, n)
                  for n, s in zip(compounds["name"], compounds["smiles"])}

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
    t0 = time.perf_counter()
    for i, (tgt, seq, name, smi) in enumerate(remaining):
        pair_start = time.perf_counter()
        try:
            r = _score_pair(boltz_exe, tgt, seq, name, smi, save_pose=args.save_poses)
            new_rows.append(r)
            ok = True
        except Exception as e:
            logger.warning("Pair (%s, %s) failed: %s; NaN", tgt, name, e)
            new_rows.append({
                "target_uniprot": tgt, "compound_name": name, "smiles": smi,
                "affinity_pred_value": float("nan"),
                "affinity_probability_binary": float("nan"),
                "pose_plddt": None,
                "pose_centroid_x": None, "pose_centroid_y": None, "pose_centroid_z": None,
                "pose_cif_path": None,
                "mode": "wsl2_full_pose",
                "scored_at": dt.datetime.now(dt.timezone.utc).isoformat(),
            })
            ok = False
        elapsed = time.perf_counter() - pair_start
        avg = (time.perf_counter() - t0) / (i + 1)
        eta = avg * (len(remaining) - i - 1) / 60
        logger.info("[%d/%d] %s + %s -> %.1fs (avg %.1fs, ETA %.0f min) ok=%s",
                    i + 1, len(remaining), tgt, name, elapsed, avg, eta, ok)

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
