"""Phase 3.2 - Allosteric blindness audit.

For each cognition-relevant target with known allosteric ligands in our compound
library, check whether MAMMAL ranks those allosteric ligands appropriately.
Compare against orthosteric ligands at the same target where available.

Output: reports/pipeline/allosteric_audit.md.
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

from mammal_repurposing.analysis.allosteric_audit import audit, render_markdown  # noqa: E402
from mammal_repurposing.config import DTI_SCORES_PARQUET, ensure_dirs  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("allosteric_audit")

DEFAULT_OUT = ROOT / "reports" / "pipeline" / "allosteric_audit.md"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", type=Path, default=DTI_SCORES_PARQUET)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    ensure_dirs()
    if not args.scores.exists():
        logger.error("Scores parquet not found: %s", args.scores)
        return 1

    scores = pd.read_parquet(args.scores)
    rows = audit(scores)
    md = render_markdown(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(md, encoding="utf-8")
    logger.info("Wrote allosteric audit to %s. %d target-rows; %d pass.",
                args.out, len(rows), sum(1 for r in rows if r.allosteric_passes))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
