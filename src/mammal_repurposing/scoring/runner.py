"""Grid runner: score every (target, compound) pair and write results to parquet.

Designed to survive interruptions:
    - Writes incrementally to a temp parquet every N batches (default 50).
    - On ``--resume``, skips any pair already present in the output parquet.

OOM strategy: the default batch size (16) is tuned for 12 GB VRAM with the
longest panel target (~1200 AA HCN1). On OOM, drop ``--batch-size`` to 8 or 4.
"""

from __future__ import annotations

import datetime as dt
import logging
from collections.abc import Iterator
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from mammal_repurposing.config import (
    DEFAULT_BATCH_SIZE,
    DTI_SCORES_PARQUET,
    MAMMAL_DTI_MODEL,
)
from mammal_repurposing.scoring.dti import score_batch_safe
from mammal_repurposing.scoring.model_loader import load_dti_model

logger = logging.getLogger(__name__)


def _build_grid(targets: pd.DataFrame, compounds: pd.DataFrame) -> pd.DataFrame:
    """Cartesian product of targets x compounds with metadata preserved."""
    t = targets[["uniprot", "gene", "sequence"]].rename(
        columns={"uniprot": "target_uniprot", "gene": "target_gene",
                 "sequence": "target_sequence"},
    )
    c = compounds[["name", "smiles"]].rename(
        columns={"name": "compound_name", "smiles": "compound_smiles"},
    )
    grid = t.merge(c, how="cross")
    return grid


def _filter_resume(grid: pd.DataFrame, existing: pd.DataFrame | None) -> pd.DataFrame:
    if existing is None or existing.empty:
        return grid
    key_cols = ["target_uniprot", "compound_name"]
    done = existing[key_cols].drop_duplicates()
    merged = grid.merge(done.assign(_done=True), on=key_cols, how="left")
    remaining = merged[merged["_done"].isna()].drop(columns="_done")
    skipped = len(grid) - len(remaining)
    if skipped:
        logger.info("Resume: skipping %d already-scored pairs.", skipped)
    return remaining


def _chunks(df: pd.DataFrame, size: int) -> Iterator[pd.DataFrame]:
    for start in range(0, len(df), size):
        yield df.iloc[start : start + size]


def score_grid(
    targets_path: Path,
    compounds_path: Path,
    out_path: Path = DTI_SCORES_PARQUET,
    *,
    batch_size: int = DEFAULT_BATCH_SIZE,
    flush_every_batches: int = 50,
    resume: bool = False,
    device: str | None = None,
) -> Path:
    """Score the full (targets x compounds) grid and persist as parquet."""
    targets = pd.read_parquet(targets_path)
    compounds = pd.read_parquet(compounds_path)
    logger.info("Loaded %d targets and %d compounds.", len(targets), len(compounds))

    grid = _build_grid(targets, compounds)
    logger.info("Built scoring grid: %d pairs.", len(grid))

    existing: pd.DataFrame | None = None
    if resume and out_path.exists():
        existing = pd.read_parquet(out_path)
        logger.info("Loaded %d existing scores from %s.", len(existing), out_path)

    grid = _filter_resume(grid, existing)
    if grid.empty:
        logger.info("Nothing to score (resume covered everything). Done.")
        return out_path

    model, tokenizer = load_dti_model(device=device)
    model_version = MAMMAL_DTI_MODEL

    buffer: list[dict] = []
    completed: list[pd.DataFrame] = [existing] if existing is not None else []
    written_total = 0

    pbar = tqdm(total=len(grid), desc="Scoring DTI", unit="pair")
    try:
        for batch_idx, chunk in enumerate(_chunks(grid, batch_size)):
            pairs = list(zip(chunk["target_sequence"], chunk["compound_smiles"]))
            sample_ids = [
                f"{row['target_uniprot']}/{row['compound_name']}"
                for _, row in chunk.iterrows()
            ]
            try:
                pkds = score_batch_safe(model, tokenizer, pairs, sample_ids=sample_ids)
            except RuntimeError as e:
                if "out of memory" in str(e).lower():
                    logger.exception(
                        "OOM at batch %d (size %d). Reduce --batch-size and use --resume.",
                        batch_idx, len(pairs),
                    )
                raise

            scored_at = dt.datetime.now(dt.timezone.utc).isoformat()
            for (_, row), pkd in zip(chunk.iterrows(), pkds, strict=True):
                buffer.append({
                    "target_uniprot": row["target_uniprot"],
                    "target_gene": row["target_gene"],
                    "compound_name": row["compound_name"],
                    "compound_smiles": row["compound_smiles"],
                    "predicted_pkd": pkd,
                    "model_version": model_version,
                    "scored_at": scored_at,
                })

            pbar.update(len(pairs))

            # Flush periodically so a crash doesn't lose everything.
            if (batch_idx + 1) % flush_every_batches == 0:
                _flush(buffer, completed, out_path)
                written_total += len(buffer)
                buffer.clear()
    finally:
        pbar.close()

    # Final flush
    if buffer:
        _flush(buffer, completed, out_path)
        written_total += len(buffer)
        buffer.clear()

    logger.info("Wrote %d new pair scores to %s.", written_total, out_path)
    return out_path


def _flush(buffer: list[dict], completed: list[pd.DataFrame], out_path: Path) -> None:
    """Append buffer to the parquet, atomically replacing the file."""
    if not buffer:
        return
    new_df = pd.DataFrame(buffer)
    completed.append(new_df)
    combined = pd.concat(completed, ignore_index=True, sort=False)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_suffix(out_path.suffix + ".tmp")
    combined.to_parquet(tmp, index=False)
    tmp.replace(out_path)
    logger.debug("Flushed %d new scores (total in file: %d).", len(buffer), len(combined))
    # Keep only one accumulator; collapse list into single frame
    completed.clear()
    completed.append(combined)
