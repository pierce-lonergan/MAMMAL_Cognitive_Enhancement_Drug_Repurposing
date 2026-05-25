"""Phase 2 - Cognitive composite scoring.

Reads dti_scores.parquet + targets.parquet, computes per-compound scores for
three cognitive proxy panels (working memory, processing speed, learning rate)
plus a global composite, plus polypharm breadth/weighted score.

Output: data/results/cognitive_composites.parquet.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.analysis.composites import compute_composites  # noqa: E402
from mammal_repurposing.config import (  # noqa: E402
    DTI_SCORES_PARQUET,
    POLYPHARM_PKD_THRESHOLD,
    RESULTS_DIR,
    TARGETS_PARQUET,
    ensure_dirs,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("composites")

DEFAULT_OUT = RESULTS_DIR / "cognitive_composites.parquet"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", type=Path, default=DTI_SCORES_PARQUET)
    parser.add_argument("--targets", type=Path, default=TARGETS_PARQUET)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--threshold", type=float, default=POLYPHARM_PKD_THRESHOLD)
    args = parser.parse_args()

    ensure_dirs()
    for required in (args.scores, args.targets):
        if not required.exists():
            logger.error("Missing input: %s", required)
            return 1

    scores = pd.read_parquet(args.scores)
    targets = pd.read_parquet(args.targets)
    logger.info("Loaded %d scores spanning %d targets x %d compounds.",
                len(scores), scores["target_uniprot"].nunique(),
                scores["compound_name"].nunique())

    df = compute_composites(scores, targets, threshold=args.threshold)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(args.out, index=False)
    logger.info("Wrote %d compound rows to %s.", len(df), args.out)

    head = df.head(10)
    logger.info("Top 10 by global composite:\n%s",
                head[["compound_name", "global_composite", "working_memory",
                      "processing_speed", "learning_rate", "polypharm_breadth"]]
                .to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
