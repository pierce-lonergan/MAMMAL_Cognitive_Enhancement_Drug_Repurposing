"""V3 Phase A.4 — ChEMBL evidence backstop via local SQLite mirror.

Replaces scripts/07_chembl_evidence.py (the 78-hour REST grinder) with the
same logic running against the local SQLite mirror. Same output schema, so
downstream consumers (Phase 3.1 calibration, the wet-lab shortlist) don't
care which path produced the parquet.

Expected wall-clock for 4,713 pairs: ~5-10 minutes (vs 78 hours via REST).

Output: data/results/chembl_evidence.parquet — schema:
    target_uniprot | compound_name | smiles | inchikey |
    status (CORROBORATED/AMBIGUOUS/NOVEL/CONTRADICTED) |
    n_records | best_pchembl | best_activity_type | best_standard_value_nm
"""

from __future__ import annotations

import argparse
import datetime as dt
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
from mammal_repurposing.fetchers.chembl_sqlite import lookup_pair_evidence  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("chembl_evidence_sqlite")

DEFAULT_OUT = RESULTS_DIR / "chembl_evidence.parquet"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", type=Path, default=DTI_SCORES_PARQUET,
                        help="Source of (target, compound) pairs to look up.")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--threshold", type=float, default=POLYPHARM_PKD_THRESHOLD,
                        help="Only look up pairs with predicted pKd above this "
                             f"(default {POLYPHARM_PKD_THRESHOLD}). Use 0 for all pairs.")
    parser.add_argument("--all-pairs", action="store_true",
                        help="Override --threshold and look up EVERY pair in scores.")
    parser.add_argument("--limit", type=int, default=None,
                        help="Smoke: process only first N pairs.")
    args = parser.parse_args()

    ensure_dirs()
    if not args.scores.exists():
        logger.error("Scores parquet not found: %s. Run scripts/04_score_dti.py first.",
                     args.scores)
        return 1

    scores = pd.read_parquet(args.scores)
    if args.all_pairs:
        hits = scores.copy()
        logger.info("Looking up ALL %d pairs in scores parquet.", len(hits))
    else:
        hits = scores[scores["predicted_pkd"] > args.threshold].copy()
        logger.info("DTI grid has %d pairs; %d above pKd %.1f.",
                    len(scores), len(hits), args.threshold)

    if args.limit:
        hits = hits.head(args.limit)

    rows: list[dict] = []
    started = dt.datetime.now()
    for _, row in tqdm(hits.iterrows(), total=len(hits), desc="SQLite ChEMBL lookups"):
        try:
            ev = lookup_pair_evidence(
                target_uniprot=row["target_uniprot"],
                smiles=row["compound_smiles"],
            )
            rows.append({
                **ev,
                "compound_name": row["compound_name"],
            })
        except Exception as e:
            logger.warning("Lookup failed for (%s, %s): %s",
                           row["target_uniprot"], row["compound_name"], e)
            rows.append({
                "target_uniprot": row["target_uniprot"],
                "compound_name": row["compound_name"],
                "smiles": row["compound_smiles"],
                "inchikey": None,
                "status": "NOVEL",
                "n_records": 0,
                "best_pchembl": None,
                "best_activity_type": None,
                "best_standard_value_nm": None,
            })

    elapsed = (dt.datetime.now() - started).total_seconds()
    df = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(args.out, index=False)

    breakdown = df["status"].value_counts().to_dict()
    logger.info("Wrote %d rows to %s in %.1f s (%.1f pairs/s). Status: %s",
                len(df), args.out, elapsed, len(df) / max(elapsed, 0.001), breakdown)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
