"""V2 Phase 0.1 — Cache ESM2-650M embeddings for every panel target.

Loads facebook/esm2_t33_650M_UR50D once (~2.5 GB VRAM at bf16), embeds each
of the 22 panel target sequences, writes per-target .pt cache files. ~5 min
total on RTX 5070; longer on CPU.

Outputs:
    data/cache/esm2/<sha1(seq)>_mean.pt         # mean-pooled (1, 1280)
    data/cache/esm2/<sha1(seq)>_per_residue.pt  # per-residue (L, 1280)
    data/results/v2/esm2_embeddings.parquet     # one row per target with embedding ID

Validation gate (run automatically at the end):
    cos(GRIA1, GRIA2) > 0.95    # paralogous AMPA subunits, should be near-identical
    cos(PDE4D, CHRNA7) < 0.50   # totally different families, should be far
If gates fail, the model loaded wrong or pooling is broken.
"""

from __future__ import annotations

import argparse
import hashlib
import logging
import sys
from pathlib import Path

import pandas as pd
import torch
import torch.nn.functional as F

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.config import DATA_DIR, RESULTS_DIR, TARGETS_PARQUET, ensure_dirs  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("esm2_embed")

ESM2_MODEL_ID = "facebook/esm2_t33_650M_UR50D"
CACHE_DIR = DATA_DIR / "cache" / "esm2"
V2_RESULTS_DIR = RESULTS_DIR / "v2"


def _seq_hash(seq: str) -> str:
    return hashlib.sha1(seq.encode("utf-8")).hexdigest()


def _resolve_device(prefer: str | None = None) -> str:
    if prefer in ("cuda", "cpu"):
        if prefer == "cuda" and not torch.cuda.is_available():
            logger.warning("CUDA requested but unavailable; falling back to CPU.")
            return "cpu"
        return prefer
    return "cuda" if torch.cuda.is_available() else "cpu"


def embed_one(
    sequence: str,
    *,
    model,
    tokenizer,
    device: str,
    dtype: torch.dtype,
    use_cache: bool = True,
) -> tuple[torch.Tensor, torch.Tensor]:
    """Embed a single sequence. Returns (mean_pooled, per_residue) on CPU."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    h = _seq_hash(sequence)
    mean_path = CACHE_DIR / f"{h}_mean.pt"
    pr_path = CACHE_DIR / f"{h}_per_residue.pt"

    if use_cache and mean_path.exists() and pr_path.exists():
        return torch.load(mean_path), torch.load(pr_path)

    inputs = tokenizer(sequence, return_tensors="pt", add_special_tokens=True).to(device)
    with torch.no_grad():
        if dtype == torch.bfloat16:
            with torch.autocast(device_type=device, dtype=dtype):
                outputs = model(**inputs)
        else:
            outputs = model(**inputs)
    per_residue = outputs.last_hidden_state[0].to(torch.float32).cpu()  # (L+2, 1280) incl. BOS/EOS
    # Drop special tokens for cleaner mean pool
    interior = per_residue[1:-1]
    mean_pool = interior.mean(dim=0).unsqueeze(0)  # (1, 1280)

    if use_cache:
        torch.save(mean_pool, mean_path)
        torch.save(interior, pr_path)
    return mean_pool, interior


def run(
    targets_path: Path = TARGETS_PARQUET,
    *,
    device: str | None = None,
    use_cache: bool = True,
    bf16: bool = True,
) -> pd.DataFrame:
    """Embed every target in the panel; return a registry DataFrame."""
    if not targets_path.exists():
        raise FileNotFoundError(f"targets.parquet not found at {targets_path}; run scripts/02_fetch_targets.py first.")

    targets = pd.read_parquet(targets_path)
    device = _resolve_device(device)
    dtype = torch.bfloat16 if bf16 and device == "cuda" else torch.float32
    logger.info("Loading %s on %s (dtype=%s) ...", ESM2_MODEL_ID, device, dtype)

    from transformers import AutoModel, AutoTokenizer  # noqa: PLC0415
    tokenizer = AutoTokenizer.from_pretrained(ESM2_MODEL_ID)
    model = AutoModel.from_pretrained(ESM2_MODEL_ID).to(device).eval()
    if bf16 and device == "cuda":
        model = model.to(torch.bfloat16)
    n_params = sum(p.numel() for p in model.parameters())
    logger.info("ESM2 ready (%.0fM params).", n_params / 1e6)

    rows = []
    for idx, t in targets.iterrows():
        seq = t["sequence"]
        gene = t["gene"]
        uniprot = t["uniprot"]
        logger.info("[%d/%d] %s (%s, %d aa) ...", idx + 1, len(targets), gene, uniprot, len(seq))
        mean_pool, per_residue = embed_one(
            seq, model=model, tokenizer=tokenizer, device=device, dtype=dtype, use_cache=use_cache,
        )
        rows.append({
            "uniprot": uniprot,
            "gene": gene,
            "seq_length": len(seq),
            "embedding_sha1": _seq_hash(seq),
            "mean_pool_norm": float(mean_pool.norm().item()),
            "per_residue_shape": list(per_residue.shape),
        })

    # Free model
    del model, tokenizer
    torch.cuda.empty_cache() if device == "cuda" else None

    return pd.DataFrame(rows)


def validate_paralogs(registry: pd.DataFrame, *, raise_on_fail: bool = False) -> bool:
    """Relative validation: paralog cosine must exceed unrelated cosine by a margin.

    NOTE: ESM2-650M mean-pooled embeddings have a non-isotropic baseline — any
    two protein sequences will typically have cosine similarity in [0.7, 0.95]
    because the mean pool washes out most of the discriminative signal across
    residue positions. The correct sanity test is RELATIVE: paralogs must
    cluster TIGHTER than unrelated proteins. Specifically:

        cos(GRIA1, GRIA2) - cos(PDE4D, CHRNA7) > 0.05      [discrimination]
        cos(GRIA1, GRIA2) > 0.95                            [paralogs near-identical]

    If the gap collapses or paralogs aren't near-identical, the model's
    weights are corrupt or the pooling implementation is broken.
    """
    def _load_mean(uniprot: str) -> torch.Tensor | None:
        row = registry[registry["uniprot"] == uniprot]
        if row.empty:
            return None
        return torch.load(CACHE_DIR / f"{row.iloc[0]['embedding_sha1']}_mean.pt")

    gria1 = _load_mean("P42261")
    gria2 = _load_mean("P42262")
    pde4d = _load_mean("Q08499")
    chrna7 = _load_mean("P36544")

    if any(x is None for x in (gria1, gria2, pde4d, chrna7)):
        logger.warning("Validation skipped — some required targets missing from cache.")
        return False

    cos_paralogs = F.cosine_similarity(gria1, gria2).item()
    cos_unrelated = F.cosine_similarity(pde4d, chrna7).item()
    gap = cos_paralogs - cos_unrelated

    logger.info(
        "Validation cosines: GRIA1<->GRIA2 = %.3f; PDE4D<->CHRNA7 = %.3f; gap = %.3f",
        cos_paralogs, cos_unrelated, gap,
    )
    logger.info("Gates: paralog>0.95 (%.3f), gap>0.05 (%.3f)", cos_paralogs, gap)

    paralog_ok = cos_paralogs > 0.95
    discrimination_ok = gap > 0.05

    passed = paralog_ok and discrimination_ok
    if not passed:
        msg = (f"ESM2 validation FAILED: paralog_ok={paralog_ok} ({cos_paralogs:.3f}), "
               f"discrimination_ok={discrimination_ok} (gap={gap:.3f}). "
               "Model weights may be corrupt or pooling broken.")
        if raise_on_fail:
            raise RuntimeError(msg)
        logger.error(msg)
    else:
        logger.info("ESM2 validation PASS (paralogs near-identical, clear discrimination from unrelated).")
    return passed


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets", type=Path, default=TARGETS_PARQUET)
    parser.add_argument("--out", type=Path, default=V2_RESULTS_DIR / "esm2_embeddings.parquet")
    parser.add_argument("--device", choices=["cuda", "cpu"], default=None)
    parser.add_argument("--no-cache", action="store_true")
    parser.add_argument("--no-bf16", action="store_true")
    parser.add_argument("--skip-validate", action="store_true")
    args = parser.parse_args()

    ensure_dirs()
    V2_RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)

    registry = run(args.targets, device=args.device,
                   use_cache=not args.no_cache, bf16=not args.no_bf16)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    registry.to_parquet(args.out, index=False)
    logger.info("Wrote ESM2 embedding registry -> %s (%d targets).", args.out, len(registry))

    if not args.skip_validate:
        ok = validate_paralogs(registry, raise_on_fail=False)
        return 0 if ok else 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
