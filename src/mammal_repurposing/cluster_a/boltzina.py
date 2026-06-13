"""Boltzina-style affinity scoring.

The pure Boltzina protocol (Furui & Ohue, arXiv 2508.17555) — Vina pose →
Boltz-2 affinity head, skipping the diffusion structure module — would be the
ideal default. There is currently no Boltzina pip package, however, so this
module ships TWO modes:

    Mode A (default): full `boltz predict` with both structure + affinity
        prediction enabled. Slower per call (~30-60 s per pair on RTX 5070),
        but uses the official Boltz-2 CLI verbatim and produces both the
        pose and the affinity simultaneously.

    Mode B (Boltzina-style, future v2.1): pre-generate Vina poses, feed the
        pose + sequence to the Boltz-2 affinity head only. Documented as a
        TODO — implement once the boltz repo exposes a clean Python API for
        affinity-only mode, or we vendor the Furui & Ohue reference.

Output schema (consumed by fusion + provenance):
    target_uniprot | compound_name | smiles |
    affinity_pred_value (log10 IC50 in µM) |
    affinity_probability_binary (calibrated binder probability in [0,1]) |
    pose_plddt | mode | scored_at

Cache: data/cache/boltzina/<sha1(seq)+sha1(smiles)>.json.
"""

from __future__ import annotations

import datetime as dt
import gc
import hashlib
import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import yaml

from mammal_repurposing.config import DATA_DIR

logger = logging.getLogger(__name__)

_AFFINITY_CACHE_DIR = DATA_DIR / "cache" / "boltzina"


def _find_boltz_executable() -> str:
    """Locate the boltz CLI. Prefer the env's Scripts dir over PATH because
    the subprocess inherits the parent's PATH which on Windows often doesn't
    include the active conda env's Scripts directory."""
    env_scripts = Path(sys.executable).parent / "Scripts" / "boltz.exe"
    if env_scripts.exists():
        return str(env_scripts)
    env_unix = Path(sys.executable).parent / "boltz"
    if env_unix.exists():
        return str(env_unix)
    found = shutil.which("boltz")
    if found:
        return found
    raise FileNotFoundError(
        "boltz CLI not found. Tried "
        f"{env_scripts}, {env_unix}, and PATH. Run `pip install boltz`."
    )

Mode = Literal["full", "boltzina_vina"]


@dataclass
class AffinityResult:
    target_uniprot: str
    compound_name: str
    smiles: str
    affinity_pred_value: float  # log10 IC50 in µM (more negative = stronger)
    affinity_probability_binary: float  # calibrated binder probability [0, 1]
    pose_plddt: float | None
    mode: Mode
    scored_at: str

    def as_dict(self) -> dict:
        return self.__dict__.copy()


def _pair_hash(seq: str, smiles: str) -> str:
    h = hashlib.sha1()
    h.update(seq.encode("utf-8"))
    h.update(b"||")
    h.update(smiles.encode("utf-8"))
    return h.hexdigest()


def _write_affinity_yaml(
    sequence: str,
    smiles: str,
    run_name: str,
    work_dir: Path,
) -> Path:
    """Write the Boltz YAML for a (protein, ligand) affinity prediction.

    Per the current Boltz YAML schema, affinity is requested by including a
    `properties: [affinity]` block referencing the ligand id.
    """
    payload = {
        "version": 1,
        "sequences": [
            {"protein": {"id": "A", "sequence": sequence}},
            {"ligand": {"id": "L", "smiles": smiles}},
        ],
        "properties": [{"affinity": {"binder": "L"}}],
    }
    yaml_path = work_dir / f"{run_name}.yaml"
    with open(yaml_path, "w", encoding="utf-8") as f:
        yaml.safe_dump(payload, f, sort_keys=False)
    return yaml_path


def _find_affinity_output(out_dir: Path, run_name: str):
    """Walk the output dir to find affinity + confidence JSONs anywhere.

    Boltz writes outputs to several possible locations depending on version:
        <out>/boltz_results_<run>/predictions/<run>/affinity_<run>.json
        <out>/predictions/<run>/affinity_<run>.json
        <out>/<run>/affinity_<run>.json
    We just rglob for any JSON with "affinity" in its name.
    """
    affinity_json = next(
        (p for p in out_dir.rglob("*.json") if "affinity" in p.name.lower()),
        None,
    )
    confidence_json = next(
        (p for p in out_dir.rglob("*.json") if "confidence" in p.name.lower()),
        None,
    )
    return affinity_json, confidence_json


def score_affinity(
    *,
    target_uniprot: str,
    sequence: str,
    compound_name: str,
    smiles: str,
    device: str = "cuda",
    use_cache: bool = True,
    mode: Mode = "full",
    recycling_steps: int = 3,
    diffusion_samples: int = 1,
    use_msa_server: bool = True,
    timeout_sec: int = 60 * 30,  # 30 min per pair worst-case
) -> AffinityResult:
    """Score a single (target, compound) pair via Boltz-2 affinity prediction.

    On cache hit returns the stored JSON. On miss invokes `boltz predict`.
    """
    _AFFINITY_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    pair_h = _pair_hash(sequence, smiles)
    # Include non-default scoring settings in the cache key so a result computed under a different
    # mode / recycling / diffusion is not served for a different setting. Canonical defaults yield
    # an EMPTY suffix, preserving every pre-existing cache filename (backward compatible).
    settings_suffix = ("" if (mode == "full" and recycling_steps == 3 and diffusion_samples == 1)
                       else f"_{mode}_r{recycling_steps}_d{diffusion_samples}")
    cache_path = _AFFINITY_CACHE_DIR / f"{pair_h}{settings_suffix}.json"

    if use_cache and cache_path.exists():
        with open(cache_path) as f:
            data = json.load(f)
        return AffinityResult(**data)

    if mode == "boltzina_vina":
        raise NotImplementedError(
            "Boltzina-style Vina-pose-only mode is a future v2.1 optimization. "
            "Use mode='full' for now."
        )

    try:
        import boltz  # noqa: F401, PLC0415
    except ImportError as e:
        raise ImportError("boltz not installed. Run `pip install boltz`.") from e

    run_name = f"aff_{pair_h[:10]}"
    with tempfile.TemporaryDirectory(prefix="boltz_aff_") as tmp:
        work_dir = Path(tmp)
        yaml_path = _write_affinity_yaml(sequence, smiles, run_name, work_dir)
        out_dir = work_dir / "out"
        out_dir.mkdir()

        boltz_exe = _find_boltz_executable()
        cmd: list[str] = [
            boltz_exe, "predict", str(yaml_path),
            "--out_dir", str(out_dir),
            "--recycling_steps", str(recycling_steps),
            "--diffusion_samples", str(diffusion_samples),
            "--output_format", "mmcif",
            "--accelerator", "gpu" if device == "cuda" else "cpu",
            "--devices", "1",
            # Disable cuequivariance triangle-mult kernel — the native
            # cuequivariance-ops-torch-cu12 wheel is not published for Windows;
            # the PyTorch fallback is slower but functionally correct.
            "--no_kernels",
        ]
        if use_msa_server:
            cmd.append("--use_msa_server")
        env = os.environ.copy()
        env.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")

        proc = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=timeout_sec)
        if proc.returncode != 0:
            logger.error("Boltz affinity failed (exit %d). stderr tail:\n%s",
                         proc.returncode, proc.stderr[-500:])
            raise RuntimeError(f"boltz predict (affinity) failed exit {proc.returncode}")

        aff_path, conf_path = _find_affinity_output(out_dir, run_name)
        if aff_path is None:
            raise RuntimeError(f"No affinity JSON under {out_dir}. stdout tail:\n{proc.stdout[-500:]}")
        with open(aff_path) as f:
            aff = json.load(f)
        # Boltz affinity output keys (verify against installed version on first run):
        #   affinity_pred_value         — log10 IC50 (µM)
        #   affinity_probability_binary — calibrated binder probability
        pred_val = float(aff.get("affinity_pred_value", aff.get("affinity_pred_value_1", float("nan"))))
        pred_prob = float(aff.get("affinity_probability_binary",
                                   aff.get("affinity_probability_binary_1", float("nan"))))

        pose_plddt = None
        if conf_path is not None:
            with open(conf_path) as f:
                conf = json.load(f)
            per_res = conf.get("plddt") or conf.get("per_residue_plddt") or []
            if per_res:
                pose_plddt = float(sum(per_res) / len(per_res))

    result = AffinityResult(
        target_uniprot=target_uniprot,
        compound_name=compound_name,
        smiles=smiles,
        affinity_pred_value=pred_val,
        affinity_probability_binary=pred_prob,
        pose_plddt=pose_plddt,
        mode=mode,
        scored_at=dt.datetime.now(dt.timezone.utc).isoformat(),
    )
    with open(cache_path, "w") as f:
        json.dump(result.as_dict(), f)
    return result


def free_gpu() -> None:
    """Aggressive GPU cleanup. Call between target switches."""
    import torch  # noqa: PLC0415

    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()


def score_grid(
    pairs: list[tuple[str, str, str, str]],
    *,
    device: str = "cuda",
    use_cache: bool = True,
    mode: Mode = "full",
    on_pair_complete=None,
) -> list[AffinityResult]:
    """Score a list of (target_uniprot, sequence, compound_name, smiles) tuples.

    Returns the list of AffinityResult. Caller may pass ``on_pair_complete``
    callback (e.g. for incremental flush). Free GPU between target switches.
    """
    out: list[AffinityResult] = []
    last_target: str | None = None
    for tgt, seq, name, smi in pairs:
        if last_target and tgt != last_target:
            free_gpu()
        try:
            r = score_affinity(
                target_uniprot=tgt, sequence=seq,
                compound_name=name, smiles=smi,
                device=device, use_cache=use_cache, mode=mode,
            )
            out.append(r)
            if on_pair_complete is not None:
                on_pair_complete(r)
        except Exception as e:
            logger.warning("Boltzina failed for (%s, %s): %s", tgt, name, e)
            out.append(AffinityResult(
                target_uniprot=tgt, compound_name=name, smiles=smi,
                affinity_pred_value=float("nan"),
                affinity_probability_binary=float("nan"),
                pose_plddt=None, mode=mode,
                scored_at=dt.datetime.now(dt.timezone.utc).isoformat(),
            ))
        last_target = tgt
    return out
