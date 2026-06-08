"""PERSEUS L3 - mechanism reversibility classifier (state-changing vs tone-changing).

The single most defensible COMPUTABLE bridge to persistence. A drug can only produce a
durable post-cessation cognitive change if the pathway it engages changes cell STATE,
not merely neurotransmitter TONE. We place every engaged mechanism on a 5-level
persistence-SUBSTRATE ordinal (persistence-capability is monotone in the rank):

    0 transient_signaling        - reuptake / AChE / orthosteric agonism (TONE)  -> ~0 ceiling
    1 durable_transcriptional    - sustained transcriptional programs (TrkB, NRF2, ISR)
    2 structural_ecm             - perineuronal-net / ECM remodeling
    3 self_propagating_epigenetic- HDAC / DNMT chromatin marks (self-maintaining)
    4 ablative_cell_population    - senolytic / aggregate-clearing (substrate removed)

This is NECESSARY-NOT-SUFFICIENT: a high rank is required to move off a null persistence
prior, but the final call is still gated by free-brain exposure (L1), the permissive-
window firewall (L4), and the evidence-design tier (L5). A compound earns its rank from
the curated mechanism->substrate lookup (its assigned class) OR from a structural alert
(an HDACi zinc-binding group, an NRF2 electrophile) - whichever is higher. The
"reversible-signaling negative filter": if the class is tone-changing and no structural
alert fires, state_changing=False and the persistence prior is forced to ~0 (this fires
on essentially the entire current pro-cognition panel - that is the honest point).

Curation: data/raw/persistence_substrate_classes.csv + persistence_structural_alerts.csv
(every row cited). RDKit + pandas only.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

SUBSTRATE_RANK: dict[str, int] = {
    "transient_signaling": 0,
    "durable_transcriptional": 1,
    "structural_ecm": 2,
    "self_propagating_epigenetic": 3,
    "ablative_cell_population": 4,
}
RANK_TO_NAME = {v: k for k, v in SUBSTRATE_RANK.items()}

# Covalent warheads: irreversible binding DECOUPLES PK from PD but is bounded by target
# turnover - it is NOT permanence. Reported as an informational flag, never a rank-up.
_COVALENT_SMARTS = {
    "acrylamide": "[CX3](=O)[NX3][CX3]=[CX2]",
    "chloroacetamide": "[CX3](=O)[CH2][Cl]",
    "vinyl_sulfone": "[SX4](=O)(=O)[CX3]=[CX2]",
}


@dataclass
class ReversibilityCall:
    substrate_class: str
    substrate_rank: int
    state_changing: bool          # rank >= 1 (above tone-changing)
    basis: str
    source: str                   # "class:<X>" | "alert:<X>"
    alerts: list[str] = field(default_factory=list)
    covalent_flags: list[str] = field(default_factory=list)


def load_substrate_classes(csv) -> dict[str, dict]:
    out: dict[str, dict] = {}
    if not Path(csv).exists():
        return out
    df = pd.read_csv(csv)
    for _, r in df.iterrows():
        sub = str(r["substrate_class"]).strip()
        if sub not in SUBSTRATE_RANK:
            raise ValueError(f"unknown substrate_class: {sub!r}")
        out[str(r["mechanism_class"]).strip()] = {
            "substrate": sub, "rank": int(r["substrate_rank"]), "basis": str(r["basis"])}
    return out


def load_structural_alerts(csv) -> list[dict]:
    """Compile the substrate-upgrading structural alerts (HDACi ZBG, NRF2 electrophile)."""
    from rdkit import Chem
    out: list[dict] = []
    if not Path(csv).exists():
        return out
    df = pd.read_csv(csv)
    for _, r in df.iterrows():
        patt = Chem.MolFromSmarts(str(r["smarts"]))
        if patt is None:
            logger.warning("bad alert SMARTS skipped: %s", r["alert"])
            continue
        sub = str(r["substrate_class"]).strip()
        out.append({"alert": str(r["alert"]), "patt": patt, "substrate": sub,
                    "rank": int(r["substrate_rank"]), "basis": str(r["basis"])})
    return out


def _covalent_flags(mol) -> list[str]:
    from rdkit import Chem
    flags = []
    for name, sm in _COVALENT_SMARTS.items():
        p = Chem.MolFromSmarts(sm)
        if p is not None and mol.HasSubstructMatch(p):
            flags.append(name)
    return flags


def reversibility_call(smiles: str, mechanism_class: str | None,
                       classes: dict[str, dict], alerts: list[dict]) -> ReversibilityCall:
    """Resolve a compound's persistence-substrate: max of (class lookup, structural
    alerts). state_changing = rank >= 1. Tone-changing + no alert -> forced null."""
    from rdkit import Chem
    # candidate (rank, substrate, basis, source) from the assigned class
    candidates: list[tuple[int, str, str, str]] = []
    ci = classes.get(str(mechanism_class)) if mechanism_class else None
    if ci is not None:
        candidates.append((ci["rank"], ci["substrate"], ci["basis"], f"class:{mechanism_class}"))
    else:
        candidates.append((0, "transient_signaling",
                           "no curated substrate for this class - default tone-changing",
                           "class:default"))
    alert_hits: list[str] = []
    cov: list[str] = []
    mol = Chem.MolFromSmiles(str(smiles)) if smiles else None
    if mol is not None:
        for a in alerts:
            if mol.HasSubstructMatch(a["patt"]):
                alert_hits.append(a["alert"])
                candidates.append((a["rank"], a["substrate"], a["basis"], f"alert:{a['alert']}"))
        cov = _covalent_flags(mol)
    rank, sub, basis, source = max(candidates, key=lambda c: c[0])
    return ReversibilityCall(
        substrate_class=sub, substrate_rank=rank, state_changing=rank >= 1,
        basis=basis, source=source, alerts=alert_hits, covalent_flags=cov)
