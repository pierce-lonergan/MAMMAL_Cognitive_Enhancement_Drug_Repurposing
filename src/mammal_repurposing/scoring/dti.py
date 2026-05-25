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
    sample_id: str = "unnamed",
) -> dict:
    """Run MAMMAL's preprocessing on a single (target, drug) pair.

    The ``sample_id`` is injected so MAMMAL's tokenizer error handler can
    produce a useful message instead of crashing in its own error formatter.
    """
    from mammal.examples.dti_bindingdb_kd.task import DtiBindingdbKdTask  # noqa: PLC0415

    # fuse's get_sample_id() looks for "data.sample_id" or the SAMPLE_ID key.
    # Provide both — cheap and bullet-proof against API drift in fuse.
    sample = {
        "target_seq": target_aa,
        "drug_seq": drug_smiles,
        "data.sample_id": sample_id,
    }
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
    *,
    sample_id: str = "unnamed",
) -> float:
    """Score one (target, drug) pair and return the predicted pKd.

    Raises whatever MAMMAL raises on tokenizer/forward errors. Caller in
    :func:`score_batch_safe` is responsible for catching and logging.
    """
    sample = _build_sample(target_aa, drug_smiles, tokenizer,
                           device=model.device, sample_id=sample_id)
    batch_dict = model.forward_encoder_only([sample])
    batch_dict = _postprocess(batch_dict)
    return float(batch_dict[OUTPUT_KEY][0])


def score_batch(
    model: "Mammal",
    tokenizer: "ModularTokenizerOp",
    pairs: list[tuple[str, str]],
    *,
    sample_ids: list[str] | None = None,
) -> list[float]:
    """Score a batch of (target_aa, drug_smiles) pairs.

    Returns a list of pKd values in the same order as ``pairs``. Raises on any
    pair failure — use :func:`score_batch_safe` for fault-tolerant scoring.
    """
    if not pairs:
        return []
    if sample_ids is None:
        sample_ids = [f"pair{i}" for i in range(len(pairs))]
    samples = [
        _build_sample(t, d, tokenizer, device=model.device, sample_id=sid)
        for (t, d), sid in zip(pairs, sample_ids, strict=True)
    ]
    batch_dict = model.forward_encoder_only(samples)
    batch_dict = _postprocess(batch_dict)
    raw = batch_dict[OUTPUT_KEY]
    return [float(v) for v in raw]


def score_batch_safe(
    model: "Mammal",
    tokenizer: "ModularTokenizerOp",
    pairs: list[tuple[str, str]],
    *,
    sample_ids: list[str] | None = None,
) -> list[float]:
    """Score a batch; on any batch-level failure, retry pair-by-pair so we
    can identify and skip the offending pair (returns NaN for it).

    This makes the runner robust to bad compounds (e.g. peptide SMILES that
    overflow the tokenizer max length) without losing the rest of the batch.
    """
    import math  # noqa: PLC0415

    if sample_ids is None:
        sample_ids = [f"pair{i}" for i in range(len(pairs))]
    try:
        return score_batch(model, tokenizer, pairs, sample_ids=sample_ids)
    except Exception as e:
        logger.warning(
            "Batch of %d failed at MAMMAL inference (%s: %s); retrying per-pair.",
            len(pairs), type(e).__name__, e,
        )

    out: list[float] = []
    for (t, d), sid in zip(pairs, sample_ids, strict=True):
        try:
            out.append(score_pair(model, tokenizer, t, d, sample_id=sid))
        except Exception as e:
            logger.warning(
                "Pair %s failed (%s: %s). Recording NaN. (drug_seq len=%d, target_aa len=%d)",
                sid, type(e).__name__, e, len(d), len(t),
            )
            out.append(math.nan)
    return out
