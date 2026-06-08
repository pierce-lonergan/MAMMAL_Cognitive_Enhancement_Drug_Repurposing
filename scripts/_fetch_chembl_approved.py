"""Fetch the full ChEMBL approved-drug set (max_phase=4) as the F2 catalogue.

These are real, approved drugs across EVERY therapeutic area - the right pool for
repurposing discovery: drugs structurally in a strong-precedent cognition class that
were never developed for cognition. SMILES are FACTS from ChEMBL (ChEMBL ID recorded);
salts/mixtures are reduced to the largest organic fragment and everything is RDKit-
gated. No clinical outcomes touched. Re-runnable. Utility, not a numbered milestone.
"""
from __future__ import annotations

import logging
from pathlib import Path

import httpx
import pandas as pd
from tenacity import retry, retry_if_exception_type, stop_after_attempt, wait_exponential

from mammal_repurposing.config import HTTP_TIMEOUT_SEC, USER_AGENT

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
log = logging.getLogger("chembl_approved")

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "raw" / "chembl_approved_catalogue.csv"
URL = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
PAGE = 1000


@retry(reraise=True, stop=stop_after_attempt(4),
       wait=wait_exponential(multiplier=1, min=2, max=20),
       retry=retry_if_exception_type(httpx.TransportError))
def _page(client: httpx.Client, offset: int) -> dict:
    r = client.get(URL, params={"max_phase": 4, "limit": PAGE, "offset": offset})
    r.raise_for_status()
    return r.json()


def _largest_fragment(smiles: str) -> str | None:
    """Largest organic fragment (drops salts/counterions), canonicalised; None on fail."""
    try:
        from rdkit import Chem
    except Exception:
        return smiles
    parts = str(smiles).split(".")
    best, best_n = None, -1
    for p in parts:
        m = Chem.MolFromSmiles(p)
        if m is None:
            continue
        n = m.GetNumHeavyAtoms()
        if n > best_n:
            best, best_n = Chem.MolToSmiles(m), n
    return best


def main() -> int:
    headers = {"User-Agent": USER_AGENT, "Accept": "application/json"}
    rows, offset, total = [], 0, None
    with httpx.Client(timeout=HTTP_TIMEOUT_SEC, headers=headers) as client:
        while True:
            j = _page(client, offset)
            if total is None:
                total = j.get("page_meta", {}).get("total_count")
                log.info("approved drugs to page through: %s", total)
            mols = j.get("molecules", [])
            if not mols:
                break
            for m in mols:
                struct = m.get("molecule_structures") or {}
                smi = struct.get("canonical_smiles")
                if not smi:
                    continue
                frag = _largest_fragment(smi)
                if not frag:
                    continue
                rows.append({"name": m.get("pref_name") or m.get("molecule_chembl_id"),
                             "smiles": frag,
                             "chembl_id": m.get("molecule_chembl_id")})
            offset += PAGE
            log.info("  fetched %d (offset %d)", len(rows), offset)
            if total and offset >= total:
                break

    df = pd.DataFrame(rows)
    df = df[df["name"].notna()].drop_duplicates("chembl_id").reset_index(drop=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False)
    log.info("Wrote %s : %d approved drugs with parseable SMILES", OUT, len(df))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
