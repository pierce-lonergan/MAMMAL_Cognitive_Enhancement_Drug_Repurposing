"""Fetch the B3DB logBB regression set (CC0) as the trainable spine for PERSEUS L1 Stage-3.

HONESTY NOTE (load-bearing): B3DB's regression target is logBB = log10(total brain : total
plasma), NOT the unbound brain-to-plasma partition Kp,uu that ultimately governs free target
exposure. logBB is a passive-penetration proxy that confounds plasma-protein and tissue
binding; true rat Kp,uu data are small and license-encumbered (Friden 2009 n~41; Morales 2024
n=157; CMD-FGKpuu unlicensed). We therefore build Stage-3 as an EFFLUX-AWARE, CONFORMAL logBB
regressor (a real upgrade from the binary penetration heuristic) and keep unbound Kp,uu as the
documented residual gap. B3DB is CC0, so it is committed for reproducibility.

Source: https://github.com/theochem/B3DB (Meng et al 2021, Sci Data; CC0-1.0).
Writes data/raw/logbb_b3db.csv (smiles, logbb, name, source).
"""
from __future__ import annotations

import io
import logging
from pathlib import Path

import httpx
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("fetch_logbb")

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "raw" / "logbb_b3db.csv"
URL = "https://raw.githubusercontent.com/theochem/B3DB/main/B3DB/B3DB_regression.tsv"


def _valid(smiles: str) -> bool:
    from rdkit import Chem
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
    return bool(smiles) and Chem.MolFromSmiles(str(smiles)) is not None


def main() -> int:
    L.info("Fetching B3DB regression set (CC0) ...")
    r = httpx.get(URL, timeout=60.0, follow_redirects=True)
    r.raise_for_status()
    df = pd.read_csv(io.StringIO(r.text), sep="\t")
    smi_col = next((c for c in df.columns if c.strip().lower() == "smiles"), None)
    y_col = next((c for c in df.columns if c.strip().lower() == "logbb"), None)
    name_col = next((c for c in df.columns if "name" in c.lower()), None)
    if smi_col is None or y_col is None:
        raise RuntimeError(f"unexpected columns: {list(df.columns)}")
    out = pd.DataFrame({
        "smiles": df[smi_col].astype(str).str.strip(),
        "logbb": pd.to_numeric(df[y_col], errors="coerce"),
        "name": df[name_col].astype(str) if name_col else "",
    })
    out = out.dropna(subset=["logbb"])
    out = out[out["smiles"].map(_valid)].drop_duplicates(subset="smiles").reset_index(drop=True)
    out["source"] = "B3DB(CC0)"
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False)
    L.info("Wrote %s (%d compounds; logBB %.2f..%.2f)", OUT, len(out),
           out["logbb"].min(), out["logbb"].max())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
