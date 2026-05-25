"""Stage 6 (Phase 1.1) - Score every compound against BBBP + ClinTox heads.

Three serial passes, one head at a time to keep ≤1 head in VRAM:
    1. BBBP      -> p_bbb (probability of BBB-permeable)
    2. TOXICITY  -> p_tox (probability of clinical toxicity)
    3. FDA_APPR  -> p_fda (probability of looking like an approved drug)

Output: data/results/aux_scores.parquet with columns
    compound_name, smiles, p_bbb, p_tox, p_fda, pred_bbb, pred_tox, pred_fda
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
    COMPOUNDS_PARQUET,
    RESULTS_DIR,
    ensure_dirs,
)
from mammal_repurposing.scoring.bbbp import score_bbbp_batch  # noqa: E402
from mammal_repurposing.scoring.clintox import (  # noqa: E402
    score_clintox_fda_batch,
    score_clintox_tox_batch,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("score_aux_heads")

DEFAULT_OUT = RESULTS_DIR / "aux_scores.parquet"


def _run_head(label: str, fn, smiles_list: list[str]):
    logger.info("Running %s head on %d compounds...", label, len(smiles_list))
    results = fn(smiles_list)
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compounds", type=Path, default=COMPOUNDS_PARQUET)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--device", choices=["cuda", "cpu"], default=None)
    parser.add_argument("--limit", type=int, default=None,
                        help="Process only first N compounds (smoke test).")
    args = parser.parse_args()

    ensure_dirs()
    if not args.compounds.exists():
        logger.error("Compounds parquet not found: %s", args.compounds)
        return 1

    compounds = pd.read_parquet(args.compounds)
    if args.limit:
        compounds = compounds.head(args.limit)
    logger.info("Loaded %d compounds.", len(compounds))

    smiles_list = compounds["smiles"].tolist()
    name_list = compounds["name"].tolist()

    bbb = _run_head("BBBP", lambda s: score_bbbp_batch(s, device=args.device), smiles_list)
    tox = _run_head(
        "TOXICITY",
        lambda s: score_clintox_tox_batch(s, device=args.device),
        smiles_list,
    )
    fda = _run_head(
        "FDA_APPR",
        lambda s: score_clintox_fda_batch(s, device=args.device),
        smiles_list,
    )

    df = pd.DataFrame({
        "compound_name": name_list,
        "smiles": smiles_list,
        "p_bbb": [r["score"] for r in bbb],
        "pred_bbb": [r["pred"] for r in bbb],
        "p_tox": [r["score"] for r in tox],
        "pred_tox": [r["pred"] for r in tox],
        "p_fda": [r["score"] for r in fda],
        "pred_fda": [r["pred"] for r in fda],
    })

    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(args.out, index=False)
    logger.info("Wrote %d rows to %s.", len(df), args.out)
    logger.info(
        "Summary -- BBBP mean p=%.3f (pred=1: %d/%d); "
        "TOX mean p=%.3f (pred=1: %d/%d); FDA mean p=%.3f (pred=1: %d/%d)",
        df["p_bbb"].mean(), (df["pred_bbb"] == 1).sum(), len(df),
        df["p_tox"].mean(), (df["pred_tox"] == 1).sum(), len(df),
        df["p_fda"].mean(), (df["pred_fda"] == 1).sum(), len(df),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
