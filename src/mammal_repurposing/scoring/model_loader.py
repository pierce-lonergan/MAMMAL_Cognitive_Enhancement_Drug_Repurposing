"""Lazy singleton loader for the MAMMAL DTI head + its tokenizer.

The model is loaded once per Python process and cached. Subsequent calls return
the cached references — useful when scoring scripts dispatch through Typer
subcommands or when tests share a fixture.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from mammal_repurposing.config import MAMMAL_DTI_MODEL

if TYPE_CHECKING:  # heavy imports stay out of import time
    from fuse.data.tokenizers.modular_tokenizer.op import ModularTokenizerOp
    from mammal.model import Mammal

logger = logging.getLogger(__name__)

_model: "Mammal | None" = None
_tokenizer: "ModularTokenizerOp | None" = None


def load_dti_model(
    model_id: str = MAMMAL_DTI_MODEL,
    *,
    device: str | None = None,
    force_reload: bool = False,
) -> tuple["Mammal", "ModularTokenizerOp"]:
    """Load (or return the cached) MAMMAL DTI head + tokenizer.

    Args:
        model_id: HuggingFace model ID. Defaults to the DTI head.
        device: "cuda", "cpu", or None for auto-detect (prefers CUDA).
        force_reload: discard cache and reload from disk.

    Returns:
        (model, tokenizer_op) tuple. Model is in .eval() mode with grads disabled.
    """
    global _model, _tokenizer

    if not force_reload and _model is not None and _tokenizer is not None:
        return _model, _tokenizer

    # Heavy imports happen here, not at module top
    import torch  # noqa: PLC0415
    from fuse.data.tokenizers.modular_tokenizer.op import ModularTokenizerOp  # noqa: PLC0415
    from mammal.model import Mammal  # noqa: PLC0415

    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"
    if device == "cuda" and not torch.cuda.is_available():
        logger.warning("CUDA requested but unavailable; falling back to CPU.")
        device = "cpu"

    logger.info("Loading MAMMAL DTI model %s on %s...", model_id, device)
    model = Mammal.from_pretrained(model_id)
    model.eval()
    model.to(device)
    torch.set_grad_enabled(False)

    logger.info("Loading tokenizer for %s...", model_id)
    tokenizer = ModularTokenizerOp.from_pretrained(model_id)

    _model = model
    _tokenizer = tokenizer
    n_params = sum(p.numel() for p in model.parameters())
    logger.info("MAMMAL ready (%.1fM params, device=%s).", n_params / 1e6, device)
    return model, tokenizer


def reset_cache() -> None:
    """Drop cached references (useful for tests and memory pressure)."""
    global _model, _tokenizer
    _model = None
    _tokenizer = None
