"""DTI inference wrapper around MAMMAL's released ``DtiBindingdbKdTask``.

The canonical inference pattern (verified against the HF model card for
``ibm/biomed.omics.bl.sm.ma-ted-458m.dti_bindingdb_pkd``):

    1. Build a sample_dict per pair: ``{"target_seq": AA, "drug_seq": SMILES}``.
    2. Run ``DtiBindingdbKdTask.data_preprocessing()`` on each (with
       ``norm_y_mean=None, norm_y_std=None`` because there's no ground truth
       at inference time).
    3. ``model.forward_encoder_only([sample_dict, ...])`` — accepts a list.
    4. ``DtiBindingdbKdTask.process_model_output()`` denormalizes using the
       real constants ``norm_y_mean=5.79384684128215, norm_y_std=1.33808027428196``.
    5. Read the pKd values from ``batch_dict["model.out.dti_bindingdb_kd"]``.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from mammal_repurposing.config import NORM_Y_MEAN, NORM_Y_STD

if TYPE_CHECKING:
    from fuse.data.tokenizers.modular_tokenizer.op import ModularTokenizerOp
    from mammal.model import Mammal

logger = logging.getLogger(__name__)

OUTPUT_KEY = "model.out.dti_bindingdb_kd"


def _build_sample(
    target_aa: str,
    drug_smiles: str,
    tokenizer: "ModularTokenizerOp",
    device,  # torch.device | str
) -> dict:
    """Run MAMMAL's preprocessing on a single (target, drug) pair."""
    from mammal.examples.dti_bindingdb_kd.task import DtiBindingdbKdTask  # noqa: PLC0415

    sample = {"target_seq": target_aa, "drug_seq": drug_smiles}
    return DtiBindingdbKdTask.data_preprocessing(
        sample_dict=sample,
        tokenizer_op=tokenizer,
        target_sequence_key="target_seq",
        drug_sequence_key="drug_seq",
        norm_y_mean=None,  # no ground truth at inference
        norm_y_std=None,
        device=device,
    )


def _postprocess(batch_dict: dict) -> dict:
    from mammal.examples.dti_bindingdb_kd.task import DtiBindingdbKdTask  # noqa: PLC0415

    return DtiBindingdbKdTask.process_model_output(
        batch_dict,
        scalars_preds_processed_key=OUTPUT_KEY,
        norm_y_mean=NORM_Y_MEAN,
        norm_y_std=NORM_Y_STD,
    )


def score_pair(
    model: "Mammal",
    tokenizer: "ModularTokenizerOp",
    target_aa: str,
    drug_smiles: str,
) -> float:
    """Score one (target, drug) pair and return the predicted pKd."""
    sample = _build_sample(target_aa, drug_smiles, tokenizer, device=model.device)
    batch_dict = model.forward_encoder_only([sample])
    batch_dict = _postprocess(batch_dict)
    return float(batch_dict[OUTPUT_KEY][0])


def score_batch(
    model: "Mammal",
    tokenizer: "ModularTokenizerOp",
    pairs: list[tuple[str, str]],
) -> list[float]:
    """Score a batch of (target_aa, drug_smiles) pairs.

    Returns a list of pKd values in the same order as ``pairs``. Caller is
    responsible for chunking large grids — see :mod:`mammal_repurposing.scoring.runner`.
    """
    if not pairs:
        return []
    samples = [_build_sample(t, d, tokenizer, device=model.device) for t, d in pairs]
    batch_dict = model.forward_encoder_only(samples)
    batch_dict = _postprocess(batch_dict)
    raw = batch_dict[OUTPUT_KEY]
    return [float(v) for v in raw]
