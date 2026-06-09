"""Cache ADMET-AI Pgp_Broccatelli probabilities for the B3DB training set (D1).

ADMET-AI is now installed, so we can replace the Didziapetris rule-based efflux feature in
Stage-3 with a real learned P-gp-substrate probability. This caches it ONCE (the engine then
reads the cached column at train time; novel-compound inference can call admet_ai live or fall
back to the rule). CPU. Writes data/raw/logbb_b3db_pgp.csv (smiles, pgp_prob).
"""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("cache_admet_pgp")

ROOT = Path(__file__).resolve().parents[1]
B3DB = ROOT / "data" / "raw" / "logbb_b3db.csv"
OUT = ROOT / "data" / "raw" / "logbb_b3db_pgp.csv"


def main() -> int:
    df = pd.read_csv(B3DB)
    smis = df["smiles"].astype(str).tolist()
    L.info("Loading ADMET-AI and predicting Pgp for %d compounds...", len(smis))
    from admet_ai import ADMETModel
    # num_workers=0 is REQUIRED on Windows: chemprop's DataLoader spawns workers that re-import
    # the entry module; with >0 workers under spawn the probe crashed ("DataLoader worker exited
    # unexpectedly"). 0 workers runs in-process (slower but robust). CPU-only.
    try:
        model = ADMETModel(num_workers=0)
    except TypeError:                       # older/newer signature without the kwarg
        model = ADMETModel()
    preds = model.predict(smiles=smis)
    pgp_col = next((c for c in preds.columns if c.lower().startswith("pgp")
                    and "broccatelli" in c.lower()), None)
    if pgp_col is None:
        pgp_col = next((c for c in preds.columns if "pgp" in c.lower()), None)
    if pgp_col is None:
        raise RuntimeError(f"no Pgp column in ADMET-AI output: {list(preds.columns)[:30]}")
    L.info("Using ADMET-AI column '%s'", pgp_col)
    out = pd.DataFrame({"smiles": smis, "pgp_prob": preds[pgp_col].to_numpy()})
    out.to_csv(OUT, index=False)
    L.info("Wrote %s | pgp_prob range %.3f..%.3f mean %.3f",
           OUT, out["pgp_prob"].min(), out["pgp_prob"].max(), out["pgp_prob"].mean())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
