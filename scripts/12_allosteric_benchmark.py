"""Phase 4.1 - Allosteric awareness benchmark.

Scores every (target, compound) triple in data/raw/allosteric_benchmark.csv
with MAMMAL's DTI head, then analyzes binding-mode group statistics and
within-group calibration.

Output:
    data/results/allosteric_benchmark_scored.parquet
    reports/benchmarks/allosteric_benchmark.md
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

from mammal_repurposing.analysis.benchmark import analyze_benchmark, render_markdown  # noqa: E402
from mammal_repurposing.config import (  # noqa: E402
    RAW_DIR,
    RESULTS_DIR,
    TARGETS_PARQUET,
    ensure_dirs,
)
from mammal_repurposing.scoring.dti import score_batch  # noqa: E402
from mammal_repurposing.scoring.model_loader import load_dti_model  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("benchmark")

DEFAULT_BENCHMARK_CSV = RAW_DIR / "allosteric_benchmark.csv"
DEFAULT_SCORED = RESULTS_DIR / "allosteric_benchmark_scored.parquet"
DEFAULT_REPORT = ROOT / "reports" / "benchmarks" / "allosteric_benchmark.md"


def _join_sequences(benchmark: pd.DataFrame, targets: pd.DataFrame) -> pd.DataFrame:
    targets_min = targets[["uniprot", "sequence"]].rename(columns={"uniprot": "target_uniprot"})
    merged = benchmark.merge(targets_min, on="target_uniprot", how="left")
    missing = merged["sequence"].isna().sum()
    if missing:
        logger.warning("%d benchmark rows missing target sequence (run fetch_targets first).", missing)
        merged = merged.dropna(subset=["sequence"])
    return merged


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", type=Path, default=DEFAULT_BENCHMARK_CSV)
    parser.add_argument("--targets", type=Path, default=TARGETS_PARQUET)
    parser.add_argument("--scored-out", type=Path, default=DEFAULT_SCORED)
    parser.add_argument("--report-out", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--device", choices=["cuda", "cpu"], default=None)
    args = parser.parse_args()

    ensure_dirs()
    if not args.benchmark.exists():
        logger.error("Benchmark CSV not found: %s", args.benchmark)
        return 1
    if not args.targets.exists():
        logger.error("Targets parquet not found: %s. Run scripts/02_fetch_targets.py first.",
                     args.targets)
        return 1

    benchmark = pd.read_csv(args.benchmark)
    targets = pd.read_parquet(args.targets)
    merged = _join_sequences(benchmark, targets)
    logger.info("Scoring %d benchmark (target, compound) pairs...", len(merged))

    model, tokenizer = load_dti_model(device=args.device)

    pkds: list[float] = []
    rows = list(merged.itertuples(index=False))
    for start in tqdm(range(0, len(rows), args.batch_size), desc="Benchmark"):
        chunk = rows[start : start + args.batch_size]
        pairs = [(r.sequence, r.smiles) for r in chunk]
        pkds.extend(score_batch(model, tokenizer, pairs))

    merged = merged.copy()
    merged["predicted_pkd"] = pkds
    keep = [
        "target_uniprot", "target_gene", "compound_name", "smiles", "binding_mode",
        "measured_activity_nm", "activity_type", "predicted_pkd", "notes",
    ]
    out = merged[[c for c in keep if c in merged.columns]]

    args.scored_out.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(args.scored_out, index=False)
    logger.info("Wrote %d scored benchmark rows to %s.", len(out), args.scored_out)

    rows_ana = analyze_benchmark(out)
    md = render_markdown(rows_ana)
    args.report_out.parent.mkdir(parents=True, exist_ok=True)
    args.report_out.write_text(md, encoding="utf-8")
    logger.info("Wrote benchmark report to %s.", args.report_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
