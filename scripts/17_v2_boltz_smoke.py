"""V2 Phase 0.2/0.3 smoke — run Boltz-2 affinity on ONE pair to validate the
CLI, GPU, output schema, and wall-clock budget before committing to the full
1,500-pair sweep.

Default pair: CHRNA7 (P36544) + galantamine — one of the v1 sanity-gate
failures (galantamine ranked at 22% percentile via MAMMAL pKd). If Boltz-2
returns a plausible affinity here in <2 minutes, the cluster is real and
the Phase 0.4 sweep is feasible.
"""

from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.cluster_a.boltzina import score_affinity  # noqa: E402
from mammal_repurposing.config import COMPOUNDS_PARQUET, TARGETS_PARQUET  # noqa: E402

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("boltz_smoke")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", default="P36544", help="UniProt accession (default CHRNA7)")
    parser.add_argument("--compound", default="galantamine")
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--no-cache", action="store_true")
    parser.add_argument("--no-msa", action="store_true",
                        help="Disable MSA server (faster, lower-quality)")
    args = parser.parse_args()

    targets = pd.read_parquet(TARGETS_PARQUET)
    compounds = pd.read_parquet(COMPOUNDS_PARQUET)
    t = targets[targets["uniprot"] == args.target]
    c = compounds[compounds["name"].str.lower() == args.compound.lower()]
    if t.empty or c.empty:
        logger.error("Could not find target=%s or compound=%s", args.target, args.compound)
        return 1

    seq = t.iloc[0]["sequence"]
    gene = t.iloc[0]["gene"]
    smiles = c.iloc[0]["smiles"]
    logger.info("Pair: %s (%s, %d aa) + %s (%s)", gene, args.target, len(seq),
                args.compound, smiles)

    t0 = time.perf_counter()
    result = score_affinity(
        target_uniprot=args.target,
        sequence=seq,
        compound_name=args.compound,
        smiles=smiles,
        device=args.device,
        use_cache=not args.no_cache,
        use_msa_server=not args.no_msa,
    )
    elapsed = time.perf_counter() - t0

    logger.info(
        "Boltz-2 result: affinity_pred_value=%.3f (log10 IC50 µM), "
        "binder_prob=%.3f, pose_plddt=%s, mode=%s, elapsed=%.1f s",
        result.affinity_pred_value,
        result.affinity_probability_binary,
        f"{result.pose_plddt:.1f}" if result.pose_plddt is not None else "n/a",
        result.mode,
        elapsed,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
