"""Stage 4 - Score the full target x compound grid with MAMMAL DTI head.

Reads enriched targets.parquet + compounds.parquet, runs ``score_grid``, writes
results to ``data/results/dti_scores.parquet``. Incrementally flushed so an
interruption can be resumed via ``--resume``.

Usage:
    python scripts/04_score_dti.py                       # full run
    python scripts/04_score_dti.py --resume              # skip already-scored pairs
    python scripts/04_score_dti.py --batch-size 8        # on OOM
    python scripts/04_score_dti.py --device cpu          # no-GPU fallback
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.config import (  # noqa: E402
    COMPOUNDS_PARQUET,
    DEFAULT_BATCH_SIZE,
    DTI_SCORES_PARQUET,
    TARGETS_PARQUET,
    ensure_dirs,
)
from mammal_repurposing.scoring.runner import score_grid  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("score_dti")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets", type=Path, default=TARGETS_PARQUET)
    parser.add_argument("--compounds", type=Path, default=COMPOUNDS_PARQUET)
    parser.add_argument("--out", type=Path, default=DTI_SCORES_PARQUET)
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE)
    parser.add_argument(
        "--flush-every-batches", type=int, default=50,
        help="Flush results to disk every N batches (default 50).",
    )
    parser.add_argument(
        "--resume", action="store_true",
        help="Skip pairs already present in the output parquet.",
    )
    parser.add_argument(
        "--device", choices=["cuda", "cpu"], default=None,
        help="Force device. Default: auto-detect (prefers CUDA).",
    )
    args = parser.parse_args()

    ensure_dirs()

    for required in (args.targets, args.compounds):
        if not required.exists():
            logger.error("Missing input: %s. Run earlier stages first.", required)
            return 1

    score_grid(
        targets_path=args.targets,
        compounds_path=args.compounds,
        out_path=args.out,
        batch_size=args.batch_size,
        flush_every_batches=args.flush_every_batches,
        resume=args.resume,
        device=args.device,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
