"""PERSEUS L3 - mechanism reversibility (state-changing vs tone-changing), honest 3-tier.

A durable post-cessation change requires the engaged mechanism to change cell STATE in a
way that PERSISTS AFTER THE DRUG IS GONE. The earlier 5-level ordinal conflated "this
pathway can change state" with "transient engagement self-maintains after washout"; only
the top level is durable by construction. Collapsed to three honest tiers, with a separate
self-maintenance criterion:

    transient        (0) - tone / reuptake / orthosteric; reverses on washout.
    plasticity_window(1) - opens a plasticity window (TrkB / 5-HT2A / PNN); durable IFF
                           paired with experience - NOT a standalone persistent effect.
    ablative         (2) - removes the pathological substrate (senolytic, aggregate
                           clearance); self-maintaining by construction (the thing is gone).

Crucially the substrate is now keyed to MECHANISM, supplied by the curated persistence axis
(persistence.py), NOT to the L0 structural class (which misroutes). This module only adds a
STRUCTURAL signal for compounds with no curated mechanism: senolytic-likeness (Tanimoto to a
senolytic anchor set) -> ablative. Reversible epigenetic / NRF2 chemotypes (HDACi ZBG,
fumarate, isothiocyanate) are reported as CAPABILITY FLAGS only - they are state-CAPABLE but
their engagement is reversible and self-maintenance after washout is not established, so they
do NOT upgrade the substrate (appendix #6: the BDNF->TrkB consolidation loop is self-limiting,
not bistable). Final substrate = max(mechanism-axis substrate, structural senolytic). RDKit only.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)

SUBSTRATE_RANK = {"transient": 0, "plasticity_window": 1, "ablative": 2}
RANK_TO_SUBSTRATE = {v: k for k, v in SUBSTRATE_RANK.items()}

FP_RADIUS, FP_BITS = 2, 2048
SENOLYTIC_TAU = 0.55  # ECFP4 Tanimoto to a senolytic anchor -> treat as ablative

# senolytic anchors (ablative by construction - clear senescent cells). Verified SMILES.
_SENOLYTIC_ANCHORS = {
    "quercetin": "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12",
    "dasatinib": "Cc1nc(Nc2ncc(C(=O)Nc3c(C)cccc3Cl)s2)cc(N2CCN(CCO)CC2)n1",
    "piperlongumine": "COc1cc(/C=C/C(=O)N2CCC=CC2=O)cc(OC)c1OC",
}

# state-CAPABLE but reversible chemotypes -> capability FLAGS only, never a substrate upgrade
_CAPABILITY_SMARTS = {
    "hdac_hydroxamate_ZBG": "[CX3](=O)[NX3][OX2H1]",
    "hdac_ortho_aminoanilide": "O=C([#6])Nc1ccccc1N",
    "nrf2_fumarate": "[CX3](=O)[CX3]=[CX3][CX3](=O)",
    "nrf2_isothiocyanate": "[NX2]=[CX2]=[SX1]",
}
# irreversible warheads: PK-PD decoupling, bounded by target turnover - NOT permanence
_COVALENT_SMARTS = {
    "acrylamide": "[CX3](=O)[NX3][CX3]=[CX2]",
    "chloroacetamide": "[CX3](=O)[CH2][Cl]",
}


@dataclass
class ReversibilityCall:
    substrate: str                 # transient | plasticity_window | ablative
    substrate_rank: int
    self_maintaining: bool         # persists after target disengagement?
    source: str                    # axis:mechanism | structural:senolytic | ...
    capability_flags: list[str] = field(default_factory=list)  # epigenetic/NRF2 capable (reversible)
    covalent_flags: list[str] = field(default_factory=list)
    senolytic_like: bool = False
    # pulsed-HDACi self-maintenance (Nat Genet 2025: TRANSIENT/intermittent HDAC inhibition can
    # induce partially self-maintaining gene-expression + 3D-genome-folding memory after washout
    # at a SUBSET of genes). Default OFF; a HYPOTHESIS that fires only when the curated axis
    # supplies pulsed-dosing self-maintenance evidence - it does NOT auto-promote to durable.
    pulsed_self_maintaining: bool = False


def _fp(smiles: str):
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except Exception:  # pragma: no cover
        return None
    m = Chem.MolFromSmiles(str(smiles))
    return None if m is None else AllChem.GetMorganFingerprintAsBitVect(m, FP_RADIUS, nBits=FP_BITS)


def _anchor_fps():
    if not hasattr(_anchor_fps, "_cache"):
        _anchor_fps._cache = {k: _fp(v) for k, v in _SENOLYTIC_ANCHORS.items()}
    return _anchor_fps._cache


def senolytic_similarity(smiles: str) -> float:
    from rdkit import DataStructs
    q = _fp(smiles)
    if q is None:
        return 0.0
    best = 0.0
    for fp in _anchor_fps().values():
        if fp is not None:
            best = max(best, float(DataStructs.TanimotoSimilarity(q, fp)))
    return best


def _smarts_hits(smiles: str, patterns: dict[str, str]) -> list[str]:
    from rdkit import Chem
    m = Chem.MolFromSmiles(str(smiles))
    if m is None:
        return []
    out = []
    for name, sm in patterns.items():
        p = Chem.MolFromSmarts(sm)
        if p is not None and m.HasSubstructMatch(p):
            out.append(name)
    return out


def reversibility_call(smiles: str, axis_substrate: str = "transient",
                       axis_self_maintaining: bool = False, *,
                       senolytic_tau: float = SENOLYTIC_TAU) -> ReversibilityCall:
    """Resolve the persistence substrate = max(mechanism-axis substrate, structural
    senolytic). HDACi/NRF2 chemotypes set capability flags only (reversible, not durable)."""
    axis_substrate = axis_substrate if axis_substrate in SUBSTRATE_RANK else "transient"
    cands = [(SUBSTRATE_RANK[axis_substrate], axis_substrate, "axis:mechanism")]
    sim = senolytic_similarity(smiles)
    senolytic = sim >= senolytic_tau
    if senolytic:
        cands.append((SUBSTRATE_RANK["ablative"], "ablative", f"structural:senolytic({sim:.2f})"))
    cap = _smarts_hits(smiles, _CAPABILITY_SMARTS)
    cov = _smarts_hits(smiles, _COVALENT_SMARTS)
    rank, sub, src = max(cands, key=lambda c: c[0])
    # a HDACi chemotype curated as self-maintaining = the pulsed-epigenetic-memory hypothesis
    # (Nat Genet 2025); a capability flag alone stays reversible-by-default (no auto-credit).
    hdaci = any("hdac" in c for c in cap)
    pulsed = bool(hdaci and axis_self_maintaining and sub != "ablative")
    return ReversibilityCall(
        substrate=sub, substrate_rank=rank,
        self_maintaining=bool(axis_self_maintaining or sub == "ablative"),
        source=src, capability_flags=cap, covalent_flags=cov, senolytic_like=senolytic,
        pulsed_self_maintaining=pulsed)
