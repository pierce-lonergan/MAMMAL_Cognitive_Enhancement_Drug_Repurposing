"""ESM2-650M target embedding cache.

Per the v2 research doc §3 Class A:
    Model: facebook/esm2_t33_650M_UR50D (MIT license)
    Params: 650M, 33 layers, 1280-dim embeddings
    VRAM at bf16: ~2.5 GB for sequences ≤1024 residues, batch 16
    Coverage: 22 targets → ~5 min one-time, cache to .pt files

Install:
    pip install fair-esm transformers
    (Or use the HF transformers wrapper directly; both work.)

Usage:
    from mammal_repurposing.cluster_a.esm2_embed import embed_target, embed_all_targets
    emb = embed_target("MNLAAA...", device="cuda")        # single target
    emb_df = embed_all_targets(targets_parquet_path)       # batch with caching
"""

from __future__ import annotations

import hashlib
import logging
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from mammal_repurposing.config import DATA_DIR

if TYPE_CHECKING:
    pass

logger = logging.getLogger(__name__)

ESM2_MODEL_ID = "facebook/esm2_t33_650M_UR50D"
_CACHE_DIR = DATA_DIR / "cache" / "esm2"

# Lazy singletons — load once per process
_model = None
_tokenizer = None


def _seq_hash(seq: str) -> str:
    return hashlib.sha1(seq.encode("utf-8")).hexdigest()


def _load_model(device: str = "cuda"):
    """Lazy-load ESM2-650M via HuggingFace transformers."""
    global _model, _tokenizer
    if _model is not None:
        return _model, _tokenizer
    try:
        from transformers import AutoModel, AutoTokenizer  # noqa: PLC0415
    except ImportError as e:
        raise ImportError(
            "transformers not installed. Run `pip install transformers`."
        ) from e

    logger.info("Loading ESM2-650M (~650M params, ~2.5 GB at bf16)...")
    _tokenizer = AutoTokenizer.from_pretrained(ESM2_MODEL_ID)
    _model = AutoModel.from_pretrained(ESM2_MODEL_ID).to(device).eval()
    return _model, _tokenizer


def embed_target(
    sequence: str,
    *,
    device: str = "cuda",
    use_cache: bool = True,
    pool: str = "mean",  # "mean" | "cls" | None (returns per-residue)
):
    """Embed a single target sequence. Cached at data/cache/esm2/<sha1>.pt."""
    import torch  # noqa: PLC0415

    _CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_path = _CACHE_DIR / f"{_seq_hash(sequence)}_{pool}.pt"
    if use_cache and cache_path.exists():
        return torch.load(cache_path, map_location=device)

    model, tokenizer = _load_model(device=device)
    inputs = tokenizer(sequence, return_tensors="pt").to(device)
    with torch.no_grad():
        outputs = model(**inputs)
    per_residue = outputs.last_hidden_state[0]  # (L, 1280)

    if pool == "mean":
        emb = per_residue.mean(dim=0)
    elif pool == "cls":
        emb = per_residue[0]
    else:
        emb = per_residue

    if use_cache:
        torch.save(emb.cpu(), cache_path)
    return emb


def embed_all_targets(
    targets_parquet_path: Path | str,
    *,
    device: str = "cuda",
    use_cache: bool = True,
) -> pd.DataFrame:
    """Embed every target in a parquet. Returns DataFrame with `embedding` column."""
    df = pd.read_parquet(targets_parquet_path)
    embs = []
    for _, row in df.iterrows():
        emb = embed_target(row["sequence"], device=device, use_cache=use_cache)
        embs.append(emb.cpu().numpy().tolist())
    df = df.copy()
    df["esm2_embedding"] = embs
    return df
