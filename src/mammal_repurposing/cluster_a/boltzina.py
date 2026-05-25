"""Boltzina — affinity-only mode bypassing Boltz-2's diffusion structure module.

Per the v2 research doc §3 Class A:
    Source: Furui & Ohue, arXiv 2508.17555 (Aug 24 2025)
    Workflow: AutoDock Vina poses → Boltz-2 affinity head (skip diffusion)
    Speed: ~11.8× faster than full Boltz-2 (reduced recycling + batch processing)
    Accuracy: below full Boltz-2; well above Vina/GNINA on MF-PCBA
    VRAM: ~7-8 GB affinity-only on L40S; fits RTX 5070 12 GB

This module is a STUB. Real install requires:
    pip install boltz vina rdkit-pypi
    pip install git+https://github.com/Furui-Lab/Boltzina  # hypothetical; verify

Expected output schema (consumed by fusion + provenance):
    target_uniprot | compound_name | smiles | log_ic50 | binder_prob | pose_path

For the v2 hybrid, Boltzina is called only on the TOP-N compounds surviving
the ADMET hard gates (typically top 50 by MAMMAL pKd per target). This
limits the call count to ~50 × 22 = 1100 affinity predictions, ~5-10s each
in Boltzina mode = ~3-6 hours wall-clock on cold cache.
"""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path
from typing import TypedDict

import pandas as pd

from mammal_repurposing.config import DATA_DIR

logger = logging.getLogger(__name__)

_AFFINITY_CACHE_DIR = DATA_DIR / "cache" / "boltzina"


class BoltzinaResult(TypedDict):
    target_uniprot: str
    compound_name: str
    smiles: str
    log_ic50: float          # log10(IC50 in µM); more negative = stronger
    binder_prob: float       # [0, 1] calibrated binder probability
    pose_plddt: float | None # pose-confidence proxy
    pose_path: str | None    # path to saved pose (.sdf or .pdb)


def _pair_hash(seq: str, smiles: str) -> str:
    h = hashlib.sha1()
    h.update(seq.encode("utf-8"))
    h.update(b"||")
    h.update(smiles.encode("utf-8"))
    return h.hexdigest()


def score_affinity(
    target_uniprot: str,
    sequence: str,
    compound_name: str,
    smiles: str,
    *,
    structure_path: Path | None = None,
    device: str = "cuda",
    use_cache: bool = True,
) -> BoltzinaResult:
    """Score a single (target, compound) pair with Boltzina affinity-only mode.

    Args:
        target_uniprot: identifier for the protein.
        sequence: protein AA sequence.
        compound_name: identifier for the compound.
        smiles: small-molecule SMILES.
        structure_path: pre-predicted .cif from `boltz_runner.predict_structure`;
            if None, Boltzina will use Vina with a sequence-derived pocket.
        device: "cuda" or "cpu".
        use_cache: read/write cached affinity result.

    Returns:
        BoltzinaResult dict (also persisted as JSON in the cache).
    """
    _AFFINITY_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_path = _AFFINITY_CACHE_DIR / f"{_pair_hash(sequence, smiles)}.json"
    if use_cache and cache_path.exists():
        with open(cache_path) as f:
            return json.load(f)

    try:
        import boltz  # noqa: F401, PLC0415
    except ImportError as e:
        raise ImportError(
            "boltz not installed; cannot run Boltzina affinity. "
            "Run `pip install boltz`."
        ) from e

    raise NotImplementedError(
        "Boltzina runner is a stub. Implement the call to "
        "Boltz-2 affinity head with Vina-generated pose; see Furui & Ohue 2025 "
        "(arXiv 2508.17555) for the exact pipeline."
    )


def score_grid(
    pairs: list[tuple[str, str, str, str]],
    *,
    structure_paths: dict[str, Path] | None = None,
    device: str = "cuda",
    use_cache: bool = True,
) -> pd.DataFrame:
    """Score a list of (target_uniprot, sequence, compound_name, smiles) tuples.

    Returns a DataFrame with the BoltzinaResult fields. Caches every call.
    """
    structure_paths = structure_paths or {}
    rows: list[BoltzinaResult] = []
    for tgt, seq, name, smi in pairs:
        try:
            rows.append(score_affinity(
                tgt, seq, name, smi,
                structure_path=structure_paths.get(tgt),
                device=device, use_cache=use_cache,
            ))
        except (ImportError, NotImplementedError) as e:
            logger.warning("Boltzina skipped %s/%s (%s); writing NaN.", tgt, name, e)
            rows.append({
                "target_uniprot": tgt, "compound_name": name, "smiles": smi,
                "log_ic50": float("nan"), "binder_prob": float("nan"),
                "pose_plddt": None, "pose_path": None,
            })
    return pd.DataFrame(rows)
