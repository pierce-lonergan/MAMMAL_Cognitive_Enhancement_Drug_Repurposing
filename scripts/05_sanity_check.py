"""Stage 5 - Decision gate. Validate positive controls in dti_scores.parquet.

Reads scores + compounds, runs the sanity report + polypharmacology leaderboard,
writes data/results/sanity_report.md, and EXITS NONZERO if the gate fails.

Failure modes (exit 2):
    - Any positive-control target has no expected compound in the top 20%.
    - Any negative-control compound (peripheral drug) surfaces in any target's
      top 5%.

If this fails, STOP and debug the DTI prompt formatting / model load before
trusting any downstream output. Do NOT relax thresholds to make it pass.
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

from mammal_repurposing.analysis.filters import filter_scores_grid  # noqa: E402
from mammal_repurposing.analysis.polypharm import compute_polypharm  # noqa: E402
from mammal_repurposing.analysis.sanity import build_report, write_report  # noqa: E402
from mammal_repurposing.config import (  # noqa: E402
    COMPOUNDS_PARQUET,
    DTI_SCORES_PARQUET,
    NEGATIVE_CONTROL_FLAG_PERCENTILE,
    POSITIVE_CONTROL_TOP_PERCENTILE,
    SANITY_REPORT_MD,
    SMILES_MAX_LENGTH_FOR_RANKING,
    ensure_dirs,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("sanity_check")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", type=Path, default=DTI_SCORES_PARQUET)
    parser.add_argument("--compounds", type=Path, default=COMPOUNDS_PARQUET)
    parser.add_argument("--out", type=Path, default=SANITY_REPORT_MD)
    parser.add_argument(
        "--top-percentile", type=float, default=POSITIVE_CONTROL_TOP_PERCENTILE,
        help="Positive-control compounds must rank in top N (e.g. 0.20 = top 20%%).",
    )
    parser.add_argument(
        "--neg-flag-percentile", type=float, default=NEGATIVE_CONTROL_FLAG_PERCENTILE,
        help="Flag if any negative control ranks in top N (e.g. 0.05 = top 5%%).",
    )
    parser.add_argument(
        "--no-gate", action="store_true",
        help="Write the report but always exit 0 (don't enforce the gate).",
    )
    args = parser.parse_args()

    ensure_dirs()

    for required in (args.scores, args.compounds):
        if not required.exists():
            logger.error("Missing input: %s. Run earlier stages first.", required)
            return 1

    scores = pd.read_parquet(args.scores)
    compounds = pd.read_parquet(args.compounds)
    logger.info(
        "Loaded %d scores spanning %d targets x %d compounds (raw).",
        len(scores),
        scores["target_uniprot"].nunique(),
        scores["compound_name"].nunique(),
    )

    # Filter peptides / over-long SMILES — out of distribution for MAMMAL's
    # small-molecule DTI head. These dominate top-of-distribution rankings
    # via molecular-size bias, not real binding chemistry.
    scores_filt = filter_scores_grid(scores, compounds,
                                     max_smiles_length=SMILES_MAX_LENGTH_FOR_RANKING)
    logger.info(
        "After exclusion filter: %d scores, %d compounds remaining.",
        len(scores_filt), scores_filt["compound_name"].nunique(),
    )

    report = build_report(
        scores_filt, compounds,
        top_percentile=args.top_percentile,
        neg_flag_percentile=args.neg_flag_percentile,
    )
    polypharm = compute_polypharm(scores_filt)

    write_report(report, polypharm, args.out)

    logger.info(
        "Positive controls passing: %d/%d. Negative-control hits: %d. Overall: %s",
        report.n_targets_pass,
        len(report.target_checks),
        len(report.negative_hits),
        "PASS" if report.passed else "FAIL",
    )

    if not report.passed and not args.no_gate:
        logger.error("Sanity gate FAILED. See %s for details. Do not trust downstream outputs.",
                     args.out)
        return 2

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
