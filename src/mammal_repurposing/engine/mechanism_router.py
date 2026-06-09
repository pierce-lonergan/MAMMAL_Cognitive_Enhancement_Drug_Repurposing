"""L4b - NMDA trapping-kinetics mechanism router (curated-PD lookup, abstain-by-default).

The L4b research lanes (reports/pipeline/l4b_second_window_synthesis.md) established that for the
non-serotonergic rapid-acting antidepressants the durability discriminator is OFF both the
binding-affinity axis (what MAMMAL's DTI head sees) AND the cheap-ADMET / 2D-structure axis (what the
L4 serotonergic permeability window exploits): ketamine and the negative control memantine are
descriptor-indistinguishable (clogP 2.90/2.69, TPSA 29.1/26.0), so the L4 trick cannot transfer. The
decisive property is channel TRAPPING / resting-state, use-dependent block (Gideons 2014) - a
pharmacodynamic fact not derivable from structure.

This module therefore does NOT predict durability from structure. It:
  1. recognises the NMDA-channel-blocker SCAFFOLD class so PERSEUS can ABSTAIN-WITH-REASON rather
     than silently miss the compound, and
  2. looks up a small CURATED, citation-backed pharmacodynamic table
     (data/raw/nmda_trapping_table.csv) to return the verified durability verdict: WINDOW for
     confirmed resting-block trappers with established durability (ketamine/esketamine/arketamine),
     NEGATIVE for the established non-durable blockers (memantine/amantadine/lanicemine), ABSTAIN for
     scaffold members whose trapping/durability is unmeasured or contested (HNK/N2O/PCP/MK-801).

Pre-registered negative (mirrors the serotonin-vs-DMT proof): a STRUCTURE-ONLY model cannot separate
memantine from ketamine - only the curated PD table does. RDKit + csv only, CPU, CI-safe (the table
IS the data; absence of the table degrades gracefully to scaffold-abstain).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

_TABLE = Path(__file__).resolve().parents[3] / "data" / "raw" / "nmda_trapping_table.csv"

# arylcyclohexylamine: an sp3 carbon shared by a 6-membered carbocycle, an exocyclic amine, and an
# aromatic ring (ketamine / PCP / tiletamine / hydroxynorketamine core).
_ARYLCYCLOHEXYLAMINE = "[NX3,NX4][CX4]1([cX3])[#6][#6][#6][#6][#6]1"
# adamantane cage (aminoadamantane: memantine / amantadine), matched as a ring-system substructure.
_ADAMANTANE = "C1C2CC3CC1CC(C2)C3"


@dataclass
class MechanismRouterCall:
    mechanism_class: str | None              # "nmda_channel_blocker" | None
    verdict: str | None = None               # WINDOW | NEGATIVE | ABSTAIN | None (no class)
    window: bool = False                      # curated durability-window-positive
    scaffold: str | None = None              # arylcyclohexylamine | aminoadamantane | <curated>
    compound: str | None = None
    pmid: str | None = None
    reasons: list[str] = field(default_factory=list)


def _mol(smiles):
    try:
        from rdkit import Chem
        from rdkit import RDLogger
        RDLogger.DisableLog("rdApp.*")
    except Exception:  # pragma: no cover
        return None
    return Chem.MolFromSmiles(str(smiles))


def _canon(smiles) -> str | None:
    mol = _mol(smiles)
    if mol is None:
        return None
    from rdkit import Chem
    return Chem.MolToSmiles(mol)


@lru_cache(maxsize=1)
def _load_table() -> dict:
    """RDKit-canonical-SMILES -> curated row. Empty dict if the table or rdkit/csv is unavailable
    (engine then degrades to scaffold-only abstain - never crashes)."""
    out: dict[str, dict] = {}
    if not _TABLE.exists():
        return out
    try:
        import csv

        from rdkit import Chem
        with open(_TABLE, newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                m = Chem.MolFromSmiles(str(row.get("smiles", "")))
                key = Chem.MolToSmiles(m) if m is not None else row.get("smiles")
                if key:
                    out[key] = row
    except Exception:  # pragma: no cover
        return {}
    return out


def nmda_scaffold(smiles) -> str | None:
    """arylcyclohexylamine or aminoadamantane scaffold (NMDA-channel-blocker chemotypes), else None."""
    from rdkit import Chem
    mol = _mol(smiles)
    if mol is None:
        return None
    p = Chem.MolFromSmarts(_ARYLCYCLOHEXYLAMINE)
    if p is not None and mol.HasSubstructMatch(p):
        return "arylcyclohexylamine"
    ada = Chem.MolFromSmiles(_ADAMANTANE)
    if (ada is not None and mol.HasSubstructMatch(ada)
            and any(a.GetSymbol() == "N" for a in mol.GetAtoms())):
        return "aminoadamantane"
    return None


def nmda_router(smiles) -> MechanismRouterCall:
    """L4b durability verdict for NMDA-channel blockers. Identity-matched curated PD verdict if the
    compound is in the table (WINDOW / NEGATIVE / ABSTAIN); else, if it matches an NMDA-blocker
    scaffold, ABSTAIN-with-reason (durability is trapping kinetics, not structure); else no class."""
    canon = _canon(smiles)
    row = _load_table().get(canon) if canon else None
    if row is not None:
        verdict = row.get("window_verdict")
        pmid = row.get("pmid")
        call = MechanismRouterCall(
            mechanism_class="nmda_channel_blocker", verdict=verdict,
            window=(verdict == "WINDOW"), scaffold=row.get("scaffold_class") or "curated",
            compound=row.get("compound"),
            pmid=(None if pmid in (None, "", "pending") else pmid))
        if verdict == "WINDOW":
            call.reasons.append(
                f"curated NMDA trapping-kinetics: {row.get('compound')} blocks resting NMDARs "
                f"({row.get('trapping_class')}) with established durable plasticity -> plasticity "
                "window (off-structure, curated PD)")
        elif verdict == "NEGATIVE":
            call.reasons.append(
                f"curated NMDA negative: {row.get('compound')} ({row.get('trapping_class')}, "
                f"resting-block {row.get('blocks_resting_nmdar')}) -> not a durable plastogen")
        else:
            call.reasons.append(
                f"curated NMDA abstain: {row.get('compound')} durability not established "
                f"({row.get('durable_rapid_antidepressant')})")
        return call
    scaf = nmda_scaffold(smiles)
    if scaf is not None:
        return MechanismRouterCall(
            mechanism_class="nmda_channel_blocker", verdict="ABSTAIN", scaffold=scaf,
            reasons=[f"{scaf} NMDA-channel-blocker scaffold recognized but trapping / resting-block "
                     "kinetics are not in the curated PD table -> durability not structure-derivable, "
                     "route to the evidence layer"])
    return MechanismRouterCall(mechanism_class=None,
                               reasons=["no NMDA-channel-blocker scaffold"])
