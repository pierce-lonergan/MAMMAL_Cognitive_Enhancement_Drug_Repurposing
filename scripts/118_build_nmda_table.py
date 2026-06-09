"""Build the curated L4b NMDA trapping-kinetics lookup table (data/raw/nmda_trapping_table.csv).

The L4b decision (reports/pipeline/l4b_second_window_synthesis.md): durability for NMDA-channel
blockers is decided by channel TRAPPING / resting-state block (Gideons 2014), NOT by structure
(ketamine and memantine are descriptor-indistinguishable). This small curated table encodes that
verified pharmacodynamic fact per compound; engine/mechanism_router.py reads it.

window_verdict rule (pre-registered): WINDOW iff (blocks_resting_NMDAR = yes AND
durable_rapid_antidepressant in {established, preclinical_only}); NEGATIVE iff established_negative
OR blocks_resting_NMDAR = no; ABSTAIN otherwise (trapper but durability not established).

PMIDs are populated from the citation-verified curation lane (reports/pipeline/
nmda_trapping_table_curation.md). Rows whose pmid is still 'pending' are NOT to be committed -
reconcile against the verified curation first. SMILES fetched from PubChem (or the ledger for N2O)
and RDKit-canonicalised so the engine matches by identity without network. CPU + network.
"""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("build_nmda_table")
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "raw" / "nmda_trapping_table.csv"

# provisional curated facts (verdicts follow the pre-registered rule above; PMIDs verified by the
# curation lane before commit). columns after compound:
#   scaffold_class, trapping_class, blocks_resting_nmdar, use_dependence,
#   durable_rapid_antidepressant, window_verdict, pmid, note
_CURATED = {
    "ketamine": ("arylcyclohexylamine", "full_trapping", "yes", "moderate", "established", "WINDOW",
                 "PMID:24912158", "blocks resting NMDARs in physiological Mg2+ -> BDNF (Gideons "
                 "2014); single-dose effect ~1 week then relapse"),
    "esketamine": ("arylcyclohexylamine", "full_trapping", "yes", "moderate", "established",
                   "WINDOW", "PMID:24912158", "S-enantiomer (Spravato); resting-block read-across "
                   "from ketamine, not a direct per-compound mEPSC measurement"),
    "arketamine": ("arylcyclohexylamine", "full_trapping", "yes", "moderate", "preclinical_only",
                   "WINDOW", "PMID:24912158", "R-enantiomer; longer-lasting preclinical. WEAKEST "
                   "WINDOW call: trapping/resting are read-across and mechanism is argued partly "
                   "NMDAR-independent - a conservative reviewer could move it to ABSTAIN"),
    "memantine": ("aminoadamantane", "partial_trapping", "no", "high", "established_negative",
                  "NEGATIVE", "PMID:24912158", "spares the resting pool (Gideons 2014); no MDD "
                  "efficacy (Zarate 2006, qualitative) - the pre-registered NEGATIVE"),
    "amantadine": ("aminoadamantane", "partial_trapping", "no", "high", "not_studied", "NEGATIVE",
                   "PMID:9120573", "aminoadamantane like memantine; resting-block read-across; "
                   "not a durable rapid antidepressant"),
    "hydroxynorketamine": ("arylcyclohexylamine", "low_trapping", "unknown", "low", "contested",
                           "ABSTAIN", "PMID:30796190", "target contested across verified papers "
                           "(Zanos NMDAR-independent / Suzuki NMDAR-dependent / Lumsden no-block); "
                           "abstain"),
    "nitrous oxide": ("gas", "low_trapping", "unknown", "moderate", "not_studied", "ABSTAIN",
                      "PMID:9546794", "fast/easily-reversible block (not ketamine-like trapping); "
                      "durability unestablished"),
    "dextromethorphan": ("morphinan", "low_trapping", "unknown", "moderate", "contested", "ABSTAIN",
                         "PMID:17848867", "chronically dosed, multi-target (sigma-1, SERT/NET); "
                         "resting-block unmeasured; not a clean single-dose durability story"),
    "phencyclidine": ("arylcyclohexylamine", "full_trapping", "unknown", "high", "not_studied",
                      "ABSTAIN", "PMID:2448800", "full trapper but not a therapeutic; resting-block "
                      "unmeasured; durability not studied"),
    "dizocilpine": ("dibenzocycloheptenamine", "full_trapping", "unknown", "high", "not_studied",
                    "ABSTAIN", "PMID:2448800", "MK-801 research tool; resting-block unmeasured; "
                    "durability not studied"),
    "lanicemine": ("other", "low_trapping", "no", "high", "established_negative", "NEGATIVE",
                   "PMID:24126931", "low-trapping (AZD6765; ket 86% vs lanicemine 54% trapping); "
                   "phase 2b failed (PMID 27681442) - a second mechanistic NEGATIVE control"),
}
_LEDGER_SMILES = {"nitrous oxide": "[N-]=[N+]=O"}


def main() -> int:
    from rdkit import Chem
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
    from mammal_repurposing.fetchers.pubchem import fetch_smiles

    rows = []
    for name, vals in _CURATED.items():
        smi = _LEDGER_SMILES.get(name)
        if smi is None:
            r = fetch_smiles(name)
            smi = r.get("smiles") if r else None
        if not smi:
            L.warning("no SMILES for %s - skipped", name); continue
        m = Chem.MolFromSmiles(str(smi))
        canon = Chem.MolToSmiles(m) if m is not None else smi
        (scaf, trap, rest, ud, dur, verdict, pmid, note) = vals
        rows.append(dict(compound=name, smiles=canon, scaffold_class=scaf, trapping_class=trap,
                         blocks_resting_nmdar=rest, use_dependence=ud,
                         durable_rapid_antidepressant=dur, window_verdict=verdict, pmid=pmid,
                         note=note))
    df = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False)
    L.info("Wrote %s (%d rows): %s", OUT, len(df),
           {v: int((df.window_verdict == v).sum()) for v in ("WINDOW", "NEGATIVE", "ABSTAIN")})
    pending = int((df.pmid == "pending").sum())
    if pending:
        L.warning("%d rows have pmid='pending' - reconcile with the verified curation BEFORE commit",
                  pending)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
