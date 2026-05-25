"""Stage 2 - Fetch amino-acid sequences for the cognition target panel.

Reads data/raw/targets_seed.csv, queries UniProt for each accession, and writes
data/interim/targets.parquet enriched with the AA sequence, length, gene name,
and Ensembl gene ID.

Usage:
    python scripts/02_fetch_targets.py             # idempotent (skip if parquet exists)
    python scripts/02_fetch_targets.py --force     # re-fetch even if parquet exists
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

# Allow running as a script from any cwd.
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.config import (  # noqa: E402
    TARGETS_PARQUET,
    TARGETS_SEED_CSV,
    ensure_dirs,
)
from mammal_repurposing.fetchers.uniprot import fetch_many  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("fetch_targets")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true", help="Re-fetch even if parquet exists")
    parser.add_argument(
        "--seed",
        type=Path,
        default=TARGETS_SEED_CSV,
        help=f"Seed CSV path (default: {TARGETS_SEED_CSV})",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=TARGETS_PARQUET,
        help=f"Output parquet path (default: {TARGETS_PARQUET})",
    )
    args = parser.parse_args()

    ensure_dirs()

    if args.out.exists() and not args.force:
        logger.info("Output exists at %s (use --force to re-fetch). Skipping.", args.out)
        return 0

    if not args.seed.exists():
        logger.error("Seed CSV not found: %s", args.seed)
        return 1

    seed = pd.read_csv(args.seed)
    logger.info("Loaded %d targets from seed CSV.", len(seed))

    accessions = seed["uniprot"].tolist()
    logger.info("Fetching sequences from UniProt...")
    entries = fetch_many(accessions)
    logger.info("Fetched %d sequences.", len(entries))

    enrich = pd.DataFrame(entries).rename(columns={"accession": "uniprot"})
    merged = seed.merge(enrich, on="uniprot", how="left", validate="one_to_one")

    missing_seq = merged["sequence"].isna().sum()
    if missing_seq:
        logger.error("%d targets missing sequence after fetch; aborting.", missing_seq)
        return 1

    merged["seq_length"] = merged["length"].astype("int64")
    merged = merged.drop(columns=["length"])

    logger.info(
        "Sequence length stats: min=%d, median=%d, max=%d",
        int(merged["seq_length"].min()),
        int(merged["seq_length"].median()),
        int(merged["seq_length"].max()),
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    merged.to_parquet(args.out, index=False)
    logger.info("Wrote %d targets to %s.", len(merged), args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
