"""Boltz-2 structure prediction runner via the `boltz predict` CLI.

Per the v2 research doc §3 Class A:
    Package: boltz (MIT license; Wohlwend/Passaro et al., bioRxiv 2025.06.14)
    Install: pip install boltz
    VRAM: ~8-10 GB structure mode; up to 1000 residues on 12 GB RTX 4070-class
    For >1000 residues: enable LMI4Boltz chunking + PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

CLI shape (verify against installed `boltz --help` if it drifts):
    boltz predict <input_dir_or_yaml> --out_dir <out>
        [--use_msa_server] [--recycling_steps 3] [--diffusion_samples 1]
        [--output_format mmcif|pdb] [--devices 1]
        [--accelerator gpu|cpu]

Input YAML shape (single-target structure prediction):
    version: 1
    sequences:
      - protein:
          id: A
          sequence: <FASTA-ish AA sequence>

Outputs land at `<out>/predictions/<run_name>/<run_name>_model_0.cif` and
`<out>/predictions/<run_name>/confidence_<run_name>_model_0.json` (pLDDT).
We cache the .cif + per-residue pLDDT under `data/cache/boltz_struct/`.
"""

from __future__ import annotations

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

import yaml

from mammal_repurposing.config import DATA_DIR

logger = logging.getLogger(__name__)

_STRUCT_CACHE_DIR = DATA_DIR / "cache" / "boltz_struct"


def _find_boltz_executable() -> str:
    """Same logic as cluster_a/boltzina._find_boltz_executable."""
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

# Sequence length threshold above which we engage LMI4Boltz-style chunking
# flags (research doc §3 Class A). The default 800 is a safe approximation
# of "approaching 1000 residues."
CHUNKING_LENGTH_THRESHOLD = 800


@dataclass
class StructureArtifact:
    """One predicted structure + its per-residue confidence."""

    cif_path: Path
    plddt_mean: float                  # mean per-residue pLDDT over the whole chain
    plddt_pocket_mean: float | None    # only set if a pocket spec was provided
    n_residues: int
    sequence_sha1: str


def _seq_hash(seq: str) -> str:
    return hashlib.sha1(seq.encode("utf-8")).hexdigest()


def _write_input_yaml(sequence: str, run_name: str, work_dir: Path) -> Path:
    """Write the minimal Boltz YAML for a single-protein structure prediction."""
    payload = {
        "version": 1,
        "sequences": [{"protein": {"id": "A", "sequence": sequence}}],
    }
    yaml_path = work_dir / f"{run_name}.yaml"
    with open(yaml_path, "w", encoding="utf-8") as f:
        yaml.safe_dump(payload, f, sort_keys=False)
    return yaml_path


def _find_output(out_dir: Path, run_name: str) -> tuple[Path | None, Path | None]:
    """Locate the Boltz prediction outputs in the run directory."""
    pred_root = out_dir / "predictions" / run_name
    if not pred_root.exists():
        # Older boltz versions may use a flat layout; try the parent
        pred_root = out_dir
    cif = next(pred_root.glob("*_model_0.cif"), None) or next(pred_root.glob("*.cif"), None)
    confidence = next(pred_root.glob("confidence_*_model_0.json"), None) or next(
        pred_root.glob("confidence_*.json"), None
    )
    return cif, confidence


def _parse_plddt(confidence_path: Path) -> tuple[float, list[float]]:
    """Read a Boltz confidence JSON, return (mean pLDDT, per-residue pLDDT)."""
    with open(confidence_path) as f:
        conf = json.load(f)
    # Boltz writes plddt as a list (per residue) and a complex_plddt summary.
    per_res = conf.get("plddt") or conf.get("per_residue_plddt") or []
    mean = float(conf.get("complex_plddt", 0.0)) or (sum(per_res) / max(1, len(per_res)))
    return mean, list(map(float, per_res))


def predict_structure(
    sequence: str,
    *,
    device: str = "cuda",
    use_cache: bool = True,
    use_msa_server: bool = True,
    recycling_steps: int = 3,
    diffusion_samples: int = 1,
    chunked: bool | None = None,
    timeout_sec: int = 60 * 60,  # 1 hour per target
) -> StructureArtifact:
    """Predict a target structure via Boltz-2. Returns a StructureArtifact.

    Caches the .cif + per-residue pLDDT at `data/cache/boltz_struct/<sha1>.{cif,plddt.json}`.
    """
    _STRUCT_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    seq_hash = _seq_hash(sequence)
    cache_cif = _STRUCT_CACHE_DIR / f"{seq_hash}.cif"
    cache_plddt = _STRUCT_CACHE_DIR / f"{seq_hash}.plddt.json"

    if use_cache and cache_cif.exists() and cache_plddt.exists():
        with open(cache_plddt) as f:
            plddt = json.load(f)
        logger.info("Boltz cache hit for sequence (%d aa, mean pLDDT %.1f).",
                    len(sequence), plddt["mean"])
        return StructureArtifact(
            cif_path=cache_cif,
            plddt_mean=float(plddt["mean"]),
            plddt_pocket_mean=plddt.get("pocket_mean"),
            n_residues=len(sequence),
            sequence_sha1=seq_hash,
        )

    if chunked is None:
        chunked = len(sequence) >= CHUNKING_LENGTH_THRESHOLD

    # Ensure boltz is installed
    try:
        import boltz  # noqa: F401, PLC0415
    except ImportError as e:
        raise ImportError(
            "boltz not installed. Run `pip install boltz` first."
        ) from e

    run_name = f"struct_{seq_hash[:10]}"
    with tempfile.TemporaryDirectory(prefix="boltz_") as tmp:
        work_dir = Path(tmp)
        yaml_path = _write_input_yaml(sequence, run_name, work_dir)
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
        ]
        if use_msa_server:
            cmd.append("--use_msa_server")

        env = os.environ.copy()
        env.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")
        if chunked:
            # LMI4Boltz chunking flags (research doc §3 Class A)
            cmd += [
                "--triangle_mult_gate_nchunks", "4",
                "--chunk_size_transition_z", "32",
                "--chunk_size_tri_attn", "64",
            ]

        logger.info("Running: %s", " ".join(cmd))
        proc = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=timeout_sec)
        if proc.returncode != 0:
            logger.error("Boltz failed (exit %d). stderr (last 500 chars):\n%s",
                         proc.returncode, proc.stderr[-500:])
            raise RuntimeError(f"boltz predict failed with exit {proc.returncode}")

        cif, conf = _find_output(out_dir, run_name)
        if cif is None or conf is None:
            raise RuntimeError(
                f"Could not locate Boltz outputs under {out_dir}. "
                f"stdout tail:\n{proc.stdout[-500:]}"
            )
        mean_plddt, per_res = _parse_plddt(conf)
        # Persist into the persistent cache
        cache_cif.write_bytes(cif.read_bytes())
        with open(cache_plddt, "w") as f:
            json.dump({"mean": mean_plddt, "per_residue": per_res}, f)

    logger.info("Boltz structure for seq (%d aa) mean pLDDT %.1f, cached at %s",
                len(sequence), mean_plddt, cache_cif)
    return StructureArtifact(
        cif_path=cache_cif,
        plddt_mean=mean_plddt,
        plddt_pocket_mean=None,
        n_residues=len(sequence),
        sequence_sha1=seq_hash,
    )
