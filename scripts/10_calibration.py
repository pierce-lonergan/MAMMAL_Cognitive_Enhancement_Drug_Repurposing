"""Phase 3.1 - Platt-style calibration of MAMMAL DTI vs ChEMBL ground truth.

Reads dti_scores.parquet + chembl_evidence.parquet, fits correlation + linear
regression, writes reports/calibration_report.md, exits nonzero ONLY on a
gate-failure of UNRELIABLE_RESTRICT_PANEL.
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

from mammal_repurposing.analysis.calibration import calibrate, render_markdown  # noqa: E402
from mammal_repurposing.config import (  # noqa: E402
    DTI_SCORES_PARQUET,
    RESULTS_DIR,
    TARGETS_PARQUET,
    ensure_dirs,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("calibration")

DEFAULT_CHEMBL = RESULTS_DIR / "chembl_evidence.parquet"
DEFAULT_REPORT = ROOT.parent / "reports" / "calibration_report.md"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", type=Path, default=DTI_SCORES_PARQUET)
    parser.add_argument("--chembl", type=Path, default=DEFAULT_CHEMBL)
    parser.add_argument("--targets", type=Path, default=TARGETS_PARQUET)
    parser.add_argument("--out", type=Path, default=ROOT / "reports" / "calibration_report.md")
    parser.add_argument(
        "--no-gate", action="store_true",
        help="Always exit 0 (default exits 2 if UNRELIABLE_RESTRICT_PANEL).",
    )
    args = parser.parse_args()

    ensure_dirs()
    for required in (args.scores, args.chembl, args.targets):
        if not required.exists():
            logger.error("Missing input: %s", required)
            return 1

    scores = pd.read_parquet(args.scores)
    chembl = pd.read_parquet(args.chembl)
    targets = pd.read_parquet(args.targets)

    result = calibrate(scores, chembl)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(render_markdown(result, targets), encoding="utf-8")
    logger.info("Wrote calibration report to %s.", args.out)
    logger.info("Gate: %s (n=%d, ρ=%.3f, r=%.3f, MAE=%.3f)",
                result.gate_label, result.n_pairs,
                result.spearman_rho, result.pearson_r, result.mae_log)

    if result.gate_label == "UNRELIABLE_RESTRICT_PANEL" and not args.no_gate:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
