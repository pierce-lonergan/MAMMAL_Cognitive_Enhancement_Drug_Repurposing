"""MoleculeNet head wrapper (BBBP, ClinTox toxicity, ClinTox FDA approval).

The MAMMAL repo ships a ``molnet_infer`` helper that takes a ``task_name`` and a
SMILES string, returning ``{"pred": 0|1, "score": float}`` where ``score`` is
the positive-class probability.

Supported tasks:
    - "BBBP"      -> ibm/biomed.omics.bl.sm.ma-ted-458m.moleculenet_bbbp
    - "TOXICITY"  -> ibm/biomed.omics.bl.sm.ma-ted-458m.moleculenet_clintox_tox
    - "FDA_APPR"  -> ibm/biomed.omics.bl.sm.ma-ted-458m.moleculenet_clintox_fda

VRAM management: each head loads its own ~1.8 GB model. The MCP README cautions
against >2 models loaded concurrently. Use :func:`score_task_batch` and free
between heads (see ``scripts/06_score_aux_heads.py`` for the orchestration).
"""

from __future__ import annotations

import gc
import logging
from typing import Literal, TypedDict

logger = logging.getLogger(__name__)

TaskName = Literal["BBBP", "TOXICITY", "FDA_APPR"]


class MolnetResult(TypedDict):
    smiles: str
    pred: int  # 0 or 1
    score: float  # positive-class probability in [0, 1]


def _resolve_device(device: str | None) -> str:
    if device is not None:
        return device
    import torch  # noqa: PLC0415

    return "cuda" if torch.cuda.is_available() else "cpu"


def score_task_batch(
    task_name: TaskName,
    smiles_list: list[str],
    *,
    device: str | None = None,
) -> list[MolnetResult]:
    """Score a batch of SMILES against a single MoleculeNet head.

    Loads the model once for the batch, scores serially (MAMMAL's molnet_infer
    is single-sample), then frees the model on return so the next head can load.
    """
    from mammal.examples.molnet.molnet_infer import load_model, task_infer  # noqa: PLC0415

    device = _resolve_device(device)
    logger.info("Loading MAMMAL MoleculeNet head '%s' on %s ...", task_name, device)
    task_dict = load_model(task_name=task_name, device=device)

    results: list[MolnetResult] = []
    try:
        for smiles in smiles_list:
            raw = task_infer(task_dict=task_dict, smiles_seq=smiles)
            results.append(MolnetResult(
                smiles=smiles,
                pred=int(raw.get("pred", 0)),
                score=float(raw.get("score", 0.0)),
            ))
    finally:
        _release(task_dict)

    return results


def _release(task_dict: dict) -> None:
    """Drop references to the loaded model and free CUDA cache."""
    import torch  # noqa: PLC0415

    for k in list(task_dict.keys()):
        task_dict[k] = None
    task_dict.clear()
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    logger.debug("Released MoleculeNet model and freed CUDA cache.")
