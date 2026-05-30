"""Score the panel's unscored targets with real MAMMAL DTI (finish panel → 31).

Runs the released MAMMAL `dti_bindingdb_pkd` head (cuda:0) over the targets that
are in the 31-target panel but were never scored — the 3 new (CHRM1/CHRM4/HTR6)
plus the 5 added after the original DTI runs (GRM2/GRM3/GRM5/GlyT1/HTR4) — each
against all 298 library compounds. Writes the same schema as `dti_scores.parquet`
so it merges cleanly into the V6.A grid.

MUST run in the MAMMAL venv (Python 3.12):
  .venv-mammal/Scripts/python.exe scripts/81_score_new_targets.py
"""

from __future__ import annotations

import argparse
import logging
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("score_new")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--panel", type=Path,
                    default=ROOT / "data" / "interim" / "targets.parquet")
    ap.add_argument("--dti", type=Path,
                    default=ROOT / "data" / "results" / "dti_scores.parquet")
    ap.add_argument("--grid", type=Path,
                    default=ROOT / "data" / "results" / "v2" / "v6a_grid_expanded.parquet")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "data" / "results" / "dti_scores_new8.parquet")
    ap.add_argument("--batch-size", type=int, default=8)
    args = ap.parse_args()

    from mammal_repurposing.scoring.model_loader import load_dti_model
    from mammal_repurposing.scoring.dti import score_batch_safe

    panel = pd.read_parquet(args.panel)
    dti = pd.read_parquet(args.dti)
    already = set(pd.read_parquet(args.grid)["target_uniprot"].astype(str)) \
        if args.grid.exists() else set(dti["target_uniprot"].astype(str))

    # targets in the panel but not yet in the V6.A grid
    todo = panel[~panel["uniprot"].astype(str).isin(already)].copy()
    gene = dict(zip(panel["uniprot"].astype(str), panel["gene"]))
    logger.info("Targets to score (%d): %s", len(todo),
                ", ".join(f"{gene[u]}" for u in todo["uniprot"].astype(str)))

    # 298 library compounds with SMILES
    cmpds = (dti.dropna(subset=["compound_smiles"])
             .drop_duplicates("compound_name")[["compound_name", "compound_smiles"]]
             .reset_index(drop=True))
    logger.info("Compounds: %d", len(cmpds))

    logger.info("Loading MAMMAL DTI head…")
    model, tok = load_dti_model()
    logger.info("Model on %s", model.device)

    stamp = datetime.now(timezone.utc).isoformat()
    rows = []
    bs = args.batch_size
    for _, t in todo.iterrows():
        uni = str(t["uniprot"]); g = str(t["gene"]); seq = str(t["sequence"])
        scores: list[float] = []
        for i in range(0, len(cmpds), bs):
            chunk = cmpds.iloc[i:i + bs]
            pairs = [(seq, s) for s in chunk["compound_smiles"]]
            sids = [f"{g}|{c}" for c in chunk["compound_name"]]
            scores.extend(score_batch_safe(model, tok, pairs, sample_ids=sids))
        n_ok = sum(1 for s in scores if s == s)  # non-NaN
        for (_, c), sc in zip(cmpds.iterrows(), scores):
            rows.append({
                "target_uniprot": uni, "target_gene": g,
                "compound_name": c["compound_name"],
                "compound_smiles": c["compound_smiles"],
                "predicted_pkd": sc, "model_version": "dti_bindingdb_pkd",
                "scored_at": stamp,
            })
        import numpy as np
        finite = [s for s in scores if s == s]
        logger.info("  %-7s %s: %d/%d scored, pkd mean=%.2f std=%.3f",
                    g, uni, n_ok, len(cmpds),
                    float(np.mean(finite)) if finite else float("nan"),
                    float(np.std(finite)) if finite else float("nan"))

    out_df = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_parquet(args.out, index=False)
    logger.info("Wrote %s (%d rows, %d targets × %d compounds)",
                args.out, len(out_df), todo.shape[0], len(cmpds))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
