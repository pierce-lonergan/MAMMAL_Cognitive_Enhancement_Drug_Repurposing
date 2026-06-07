"""Expand data/raw/ledger_compound_smiles.csv to cover the research-curated
ledger compounds, so the F2 novel-compound engine has class exemplars beyond the
original base-31. SMILES are FACTS fetched from PubChem (CID recorded as
provenance) and gated through RDKit (only structures RDKit can parse are kept).
No clinical outcomes are touched. Re-runnable; only appends genuinely new,
parseable rows. Utility (underscore-prefixed), not a numbered milestone script.
"""
from __future__ import annotations

import logging
import re
from pathlib import Path

import pandas as pd

from mammal_repurposing.fetchers.pubchem import fetch_many_smiles
from mammal_repurposing.reporting.trial_watch import load_combined_ledger

logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(message)s")
log = logging.getLogger("expand_smiles")

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
SMILES_CSV = RAW / "ledger_compound_smiles.csv"
LEDGERS = [
    RAW / "clinical_outcomes_ledger.csv",
    RAW / "clinical_outcomes_ledger_EXTENSION.csv",
    RAW / "clinical_outcomes_ledger_CTGOV.csv",
    RAW / "clinical_outcomes_ledger_RESEARCH.csv",
]


def _alt_names(name: str) -> list[str]:
    """Generate fallback query names: strip a trailing parenthetical alias, and for
    combination drugs (xanomeline-trospium, drug+drug) take the leading component
    (its structure is the meaningful class exemplar)."""
    alts: list[str] = []
    stripped = re.sub(r"\s*\([^)]*\)\s*$", "", name).strip()
    if stripped and stripped != name:
        alts.append(stripped)
    # combination separators -> leading active
    for sep in ("-", "+", "/", " and "):
        if sep in stripped:
            head = stripped.split(sep)[0].strip()
            if head and head not in alts and head != name:
                alts.append(head)
    return alts


def _rdkit_ok(smiles: str) -> bool:
    try:
        from rdkit import Chem
    except Exception:
        return True  # cannot gate without rdkit; accept (engine will re-parse)
    return Chem.MolFromSmiles(str(smiles)) is not None


def main() -> int:
    led = load_combined_ledger(LEDGERS)
    have = pd.read_csv(SMILES_CSV)
    have_keys = set(have["compound"].astype(str).str.lower().str.strip())

    led["_k"] = led["compound"].astype(str).str.lower().str.strip()
    missing = led[~led["_k"].isin(have_keys)]["compound"].tolist()
    log.warning("ledger=%d  have_smiles=%d  missing=%d", len(led), len(have), len(missing))

    queries = [(m, _alt_names(m)) for m in missing]
    hits = fetch_many_smiles(queries)

    new_rows, n_unresolved, n_badparse = [], 0, 0
    for name, hit in zip(missing, hits):
        if not hit["smiles"]:
            n_unresolved += 1
            continue
        if not _rdkit_ok(hit["smiles"]):
            n_badparse += 1
            log.warning("RDKit rejected %s SMILES; skipping", name)
            continue
        new_rows.append({"compound": name, "smiles": hit["smiles"],
                         "smiles_kind": hit["smiles_kind"], "cid": hit["cid"]})

    print(f"resolved+parseable: {len(new_rows)}  unresolved: {n_unresolved}  "
          f"rdkit-rejected: {n_badparse}  (of {len(missing)} missing)")
    if not new_rows:
        print("nothing new to append.")
        return 0

    out = pd.concat([have, pd.DataFrame(new_rows)], ignore_index=True)
    out = out.drop_duplicates("compound", keep="first")
    out.to_csv(SMILES_CSV, index=False)
    print(f"wrote {SMILES_CSV} : {len(have)} -> {len(out)} rows")
    print("newly added:", sorted(r["compound"] for r in new_rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
