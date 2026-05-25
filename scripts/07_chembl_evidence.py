"""Phase 1.2 - ChEMBL ground-truth backstop.

For every (target, compound) pair in dti_scores.parquet where predicted pKd
exceeds the threshold (default 6.0), look up ChEMBL for reported activity at
that target and classify CORROBORATED / NOVEL / CONTRADICTED / INCONCLUSIVE.

Output: data/results/chembl_evidence.parquet.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd
from tqdm import tqdm

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.config import (  # noqa: E402
    DTI_SCORES_PARQUET,
    POLYPHARM_PKD_THRESHOLD,
    RESULTS_DIR,
    ensure_dirs,
)
from mammal_repurposing.fetchers.chembl_groundtruth import lookup_grid  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("chembl_evidence")

DEFAULT_OUT = RESULTS_DIR / "chembl_evidence.parquet"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", type=Path, default=DTI_SCORES_PARQUET)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument(
        "--threshold", type=float, default=POLYPHARM_PKD_THRESHOLD,
        help="Only look up pairs with predicted pKd above this (default 6.0).",
    )
    parser.add_argument(
        "--limit", type=int, default=None,
        help="Process only first N candidate pairs (smoke test).",
    )
    args = parser.parse_args()

    ensure_dirs()
    if not args.scores.exists():
        logger.error("Scores parquet not found: %s. Run scripts/04_score_dti.py first.", args.scores)
        return 1

    scores = pd.read_parquet(args.scores)
    hits = scores[scores["predicted_pkd"] > args.threshold].copy()
    logger.info("DTI grid has %d pairs; %d above pKd %.1f.",
                len(scores), len(hits), args.threshold)

    if args.limit:
        hits = hits.head(args.limit)

    pairs = list(zip(hits["target_uniprot"], hits["compound_name"], hits["compound_smiles"]))
    logger.info("Looking up ChEMBL evidence for %d (target, compound) pairs...", len(pairs))

    rows = []
    # Process in batches just for progress reporting; lookup_grid handles HTTP reuse.
    batch = 25
    for start in tqdm(range(0, len(pairs), batch), desc="ChEMBL lookups"):
        chunk = pairs[start : start + batch]
        rows.extend(lookup_grid(chunk))

    df = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(args.out, index=False)

    breakdown = df["label"].value_counts().to_dict()
    logger.info("Wrote %d rows to %s. Label breakdown: %s",
                len(df), args.out, breakdown)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
