"""Phase 1.3 - OpenTargets cross-reference for the target panel.

For each target in targets.parquet (with Ensembl ID), fetch:
    - Known drugs targeting it (with disease, mechanism, phase)
    - Count of CNS/cognition-relevant drugs vs total

Output: data/results/opentargets_context.parquet — one row per target with
flattened columns (n_known_drugs_total, n_cns_drugs, top_cns_drugs,
top_cns_diseases).
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

from mammal_repurposing.config import (  # noqa: E402
    RESULTS_DIR,
    TARGETS_PARQUET,
    ensure_dirs,
)
from mammal_repurposing.fetchers.opentargets import fetch_contexts_for_targets  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("opentargets")

DEFAULT_OUT = RESULTS_DIR / "opentargets_context.parquet"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets", type=Path, default=TARGETS_PARQUET)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument(
        "--per-target-size", type=int, default=200,
        help="Max known-drug rows to fetch per target from OpenTargets.",
    )
    args = parser.parse_args()

    ensure_dirs()
    if not args.targets.exists():
        logger.error("Targets parquet not found: %s. Run scripts/02_fetch_targets.py first.",
                     args.targets)
        return 1

    targets = pd.read_parquet(args.targets)
    if "ensembl_gene_id" not in targets.columns:
        logger.error("Targets parquet missing ensembl_gene_id column. Re-fetch with --force.")
        return 1

    pairs = list(zip(targets["ensembl_gene_id"], targets["uniprot"]))
    logger.info("Fetching OpenTargets context for %d targets...", len(pairs))
    contexts = fetch_contexts_for_targets(pairs, size=args.per_target_size)
    logger.info("Got %d context rows back.", len(contexts))

    rows = []
    for c in contexts:
        top_cns_drugs = [
            kd["drug_name"]
            for kd in c["known_drugs"]
            if kd["disease_name"] and any(
                d.lower() in kd["disease_name"].lower()
                for d in c["cns_disease_associations"][:5]
            )
        ][:10]
        rows.append({
            "target_uniprot": c["target_uniprot"],
            "ensembl_id": c["ensembl_id"],
            "approved_symbol": c["approved_symbol"],
            "n_known_drugs_total": c["n_known_drugs_total"],
            "n_cns_drugs": c["n_cns_drugs"],
            "top_cns_drugs": ";".join(top_cns_drugs),
            "top_cns_diseases": ";".join(c["cns_disease_associations"][:10]),
        })

    df = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(args.out, index=False)
    logger.info("Wrote %d rows to %s.", len(df), args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
