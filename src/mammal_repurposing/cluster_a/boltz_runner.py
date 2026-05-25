"""Boltz-2 structure prediction runner (STUB).

Per the v2 research doc §3 Class A:
    Package: boltz (MIT license; Wohlwend/Passaro et al., bioRxiv 2025.06.14)
    Install: pip install boltz
    VRAM: ~8-10 GB structure mode; up to 1000 residues on 12 GB RTX 4070-class
    For >1000 residues: enable LMI4Boltz chunking + PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

This module is a STUB — the actual Boltz-2 install is a 5+ GB download. The
runner interface here documents the expected output schema so downstream
fusion + provenance can consume it as soon as the real runs land.

Workflow:
    1. For each cognition target without a usable PDB co-crystal (see research
       doc §3 Class A target inventory table), predict structure with Boltz-2
       and cache the .cif file at `data/cache/boltz_struct/<sha1(seq)>.cif`.
    2. Affinity prediction is handled by `boltzina.py` (faster mode).
"""

from __future__ import annotations

import hashlib
import logging
from pathlib import Path

from mammal_repurposing.config import DATA_DIR

logger = logging.getLogger(__name__)

_STRUCT_CACHE_DIR = DATA_DIR / "cache" / "boltz_struct"


def _seq_hash(seq: str) -> str:
    return hashlib.sha1(seq.encode("utf-8")).hexdigest()


def predict_structure(
    sequence: str,
    *,
    device: str = "cuda",
    use_cache: bool = True,
    chunked: bool = False,
) -> Path:
    """Predict a target structure via Boltz-2. Returns the .cif path.

    Args:
        sequence: protein amino acid sequence.
        device: "cuda" or "cpu".
        use_cache: read/write `data/cache/boltz_struct/<sha1>.cif`.
        chunked: enable LMI4Boltz chunking flags for sequences > 1000 residues.

    Raises:
        ImportError if `boltz` is not installed.
        NotImplementedError until the real boltz call is wired up.
    """
    _STRUCT_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_path = _STRUCT_CACHE_DIR / f"{_seq_hash(sequence)}.cif"
    if use_cache and cache_path.exists():
        logger.info("Boltz cache hit for sequence (%d aa).", len(sequence))
        return cache_path

    try:
        import boltz  # noqa: F401, PLC0415
    except ImportError as e:
        raise ImportError(
            "boltz not installed. Run `pip install boltz` first. "
            "On Windows + RTX 5070 also set "
            "PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True before invoking."
        ) from e

    # TODO: wire up actual boltz CLI/API call
    raise NotImplementedError(
        "Boltz-2 runner is a stub. Install boltz, then implement the "
        "subprocess call: `boltz predict <input.yaml> --out_dir <out>`. "
        "See https://github.com/jwohlwend/boltz for current CLI syntax."
    )
