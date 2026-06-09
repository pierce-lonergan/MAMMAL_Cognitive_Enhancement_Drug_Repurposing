"""L4 - psychoplastogen plasticity-window detector (structure-computable, OFF the DTI axis).

The deep-research gap analysis (reports/pipeline/persistence_engineering_gaps.md) established
that durability for the psychedelic class is NOT a binding-affinity problem MAMMAL's DTI head
can see. The decisive mechanism (Vargas 2023, Science 379:700): psychedelics promote durable
structural plasticity by engaging the INTRACELLULAR 5-HT2A pool, and serotonin is
non-plastogenic ONLY because it cannot cross the membrane to reach it - serotonin and DMT are
near-isoaffine at 5-HT2A, so the discriminator is membrane PERMEABILITY (lipophilicity), an
ADMET/physicochemical property, not a DTI score. Casarotto 2021 adds the parallel TrkB-TMD
mechanism (also off the sequence-DTI axis).

This module operationalizes that as a structure-derivable plasticity-window flag:
  window = (serotonergic/monoaminergic-agonist SCAFFOLD)
           AND (intracellular ACCESS: lipophilic enough to cross the membrane, clogP gate)
           AND (CNS-penetrant)
It NEVER promotes to durable-by-construction (the L4 permissive/instructive firewall): a window
is direction-neutral and durable ONLY if paired with experience. The decisive validation is
that serotonin is window-NEGATIVE while DMT/psilocin/mescaline are window-POSITIVE despite
near-identical 5-HT2A affinity - proving the discriminator lives off the DTI axis.

RDKit-only, CPU, CI-safe.
"""
from __future__ import annotations

from dataclasses import dataclass, field

# intracellular-access (Vargas 2023): the plastogenic difference between serotonin and
# psychedelics is membrane permeability, NOT 5-HT2A affinity. The reliable passive-permeability
# determinants (Lipinski/Veber) are polar surface area + H-bond donors; serotonin is the
# canonical "too polar to cross" amine (TPSA 62, HBD 3: indole-NH + phenol-OH + primary amine),
# while DMT/psilocin/mescaline sit at TPSA <=54, HBD <=2. (Crippen clogP is unreliable here -
# it overestimates serotonin at ~1.4 - so it is only a weak floor, not the discriminator.)
ACCESS_CLOGP_MIN = 0.0
ACCESS_TPSA_MAX = 60.0
ACCESS_HBD_MAX = 2

# serotonergic / monoaminergic psychoplastogen scaffolds (5-HT2A-axis agonist chemotypes)
# indol-3-yl-ethylamine with a SMALL basic amine (primary or di-small-alkyl). The N-substituent
# vetoes (no N-benzyl, no N-aryl, no amide) separate 5-HT2A-agonist psychedelic tryptamines
# (DMT/psilocin: N is H or methyl) from bulky-amine tryptamine drugs that are NOT plastogenic
# agonists (e.g. idalopirdine, a 5-HT6 antagonist with an N-benzyl arm).
_TRYPTAMINE = "c1ccc2c(c1)c(c[nH]2)CC[NX3;!$([NX3][CH2][a]);!$([NX3][a]);!$([NX3]C=O)]"
# isoDMT (D5 recall innovation): the aminoethyl is on the indole N1, not C3 (Dunlap/Olson isoDMT
# psychoplastogens, e.g. zalsupindole/AAZ-A-154). The C3-requiring _TRYPTAMINE SMARTS misses these,
# so they are detected separately. Same small-amine guard (no N-aryl / amide) as the tryptamine.
_ISODMT = "c1ccc2c(c1)c[cH]n2[CH2][CH1,CH2][NX3;!$([NX3][a]);!$([NX3]C=O)]"
_PHENETHYL = "[cX3]1[cX3][cX3][cX3][cX3][cX3]1[CH2][CH1,CH2][NX3]"  # 2-aryl-ethyl/propyl-amine core
_ERGOLINE_SMILES = "C1CN(C)C2Cc3c[nH]c4cccc(c34)C2=C1"  # ergoline ring system (LSD/lysergamides)
_AROM_OME = "[c][OX2][CH3]"                            # aromatic methoxy (psychedelic phenethylamines)
_AROM_HALO = "[c][F,Cl,Br,I]"                          # aromatic halogen (DOI/DOB/2C-x)


@dataclass
class PsychoplastogenCall:
    window: bool                       # plasticity-window-positive?
    scaffold: str | None = None        # tryptamine | isodmt | psychedelic_phenethylamine | ergoline
    clogp: float = float("nan")
    tpsa: float = float("nan")
    hbd: int = -1
    intracellular_access: bool = False
    reasons: list[str] = field(default_factory=list)
    caveat: str = ("plasticity_window is direction-NEUTRAL and durable ONLY if paired with "
                   "experience (permissive, not instructive); never auto-durable")


def _mol(smiles: str):
    try:
        from rdkit import Chem
        from rdkit import RDLogger
        RDLogger.DisableLog("rdApp.*")
    except Exception:  # pragma: no cover
        return None
    return Chem.MolFromSmiles(str(smiles))


# Scaffold vetoes (adversarial-audit fixes): indole/ergoline-bearing drugs that share a
# psychedelic scaffold but are NOT 5-HT2A-agonist psychoplastogens. Scaffold+permeability cannot
# resolve receptor SUBTYPE or AGONISM from structure, so a few decisive non-agonist motifs are
# rejected outright:
#   - sulfonamide / cyclic carbamate (oxazolidinone) -> TRIPTAN pharmacophore (5-HT1B/1D agonists:
#     sumatriptan/zolmitriptan), D4 audit.
#   - aliphatic THIOETHER (C-S-C) -> dopaminergic CLAVINE ergolines (pergolide), D5 audit. Verified
#     absent from every compound in the persistence positive ledger, so this is collateral-free.
# DOCUMENTED RESIDUAL LIMIT (D5): the lysergamide ANTAGONIST methysergide (a 5-HT2 antagonist) is
# NOT separable here - its only structural marker vs LSD is indole-N1-methylation, which is exactly
# the feature that defines the isoDMT psychoplastogens we now DETECT as positives (zalsupindole /
# AAZ-A-154). So N1-substitution must NOT be vetoed (it would suppress the whole isoDMT class);
# agonist-vs-antagonist within the lysergamide series is left to the functional/DTI layer, not this
# structural screen. methysergide is therefore a known, accepted false-positive of the window.
_SCAFFOLD_VETO = [
    "[SX4](=O)(=O)[#7]",   # sulfonamide (sumatriptan)
    "[NX3]C(=O)[OX2]",     # carbamate / oxazolidinone (zolmitriptan)
    "[CX4][SX2][CX4]",     # aliphatic thioether (pergolide); no ledger psychoplastogen has one
]


def serotonergic_scaffold(smiles: str) -> str | None:
    """Identify a 5-HT2A-axis agonist chemotype: tryptamine (C3-aminoethyl indole), isoDMT
    (N1-aminoethyl indole), psychedelic phenethylamine (2+ aromatic OMe/halogen substituents), or
    ergoline. None if no match. Non-agonist scaffold motifs (triptan sulfonamide/carbamate;
    dopaminergic-clavine thioether) are vetoed."""
    from rdkit import Chem
    mol = _mol(smiles)
    if mol is None:
        return None
    if any((vp := Chem.MolFromSmarts(v)) is not None and mol.HasSubstructMatch(vp)
           for v in _SCAFFOLD_VETO):
        return None
    if (p := Chem.MolFromSmarts(_TRYPTAMINE)) is not None and mol.HasSubstructMatch(p):
        return "tryptamine"
    if (ip := Chem.MolFromSmarts(_ISODMT)) is not None and mol.HasSubstructMatch(ip):
        return "isodmt"
    erg = Chem.MolFromSmiles(_ERGOLINE_SMILES)
    if erg is not None and mol.HasSubstructMatch(erg):
        return "ergoline"
    pe = Chem.MolFromSmarts(_PHENETHYL)
    if pe is not None and mol.HasSubstructMatch(pe):
        ome = Chem.MolFromSmarts(_AROM_OME)
        halo = Chem.MolFromSmarts(_AROM_HALO)
        n_ome = len(mol.GetSubstructMatches(ome)) if ome is not None else 0
        n_halo = len(mol.GetSubstructMatches(halo)) if halo is not None else 0
        if n_ome + n_halo >= 2:                       # mescaline (3 OMe), DOI (2 OMe + I)
            return "psychedelic_phenethylamine"
    return None


def psychoplastogen_window(smiles: str) -> PsychoplastogenCall:
    """Structure-derivable plasticity-window verdict. Fires only when a serotonergic-agonist
    scaffold co-occurs with intracellular access (lipophilic enough to reach the intracellular
    5-HT2A pool). Encodes Vargas 2023: affinity is necessary-not-sufficient; PERMEABILITY is the
    discriminator that separates plastogenic psychedelics from non-plastogenic serotonin."""
    from mammal_repurposing.engine.cns_exposure import cns_mpo_like
    mpo = cns_mpo_like(smiles)
    if mpo is None:
        return PsychoplastogenCall(False, reasons=["unparseable SMILES"])
    scaf = serotonergic_scaffold(smiles)
    clogp, tpsa, hbd = mpo["clogp"], mpo["tpsa"], mpo["hbd"]
    access = clogp >= ACCESS_CLOGP_MIN and tpsa <= ACCESS_TPSA_MAX and hbd <= ACCESS_HBD_MAX
    call = PsychoplastogenCall(False, scaffold=scaf, clogp=round(clogp, 2),
                               tpsa=round(tpsa, 1), hbd=int(hbd), intracellular_access=access)
    if scaf is None:
        call.reasons.append("no serotonergic/monoaminergic-agonist scaffold")
        return call
    if not access:
        call.reasons.append(
            f"{scaf} scaffold but NO intracellular access (TPSA {tpsa:.0f} > {ACCESS_TPSA_MAX:.0f} "
            f"or HBD {hbd} > {ACCESS_HBD_MAX} or clogP {clogp:.2f} < {ACCESS_CLOGP_MIN}) - "
            "membrane-impermeant like serotonin")
        return call
    call.window = True
    call.reasons.append(
        f"{scaf} + intracellular access (TPSA {tpsa:.0f}, HBD {hbd}, clogP {clogp:.2f}) -> "
        "plasticity window (permissive; durable only if paired with experience)")
    return call


if __name__ == "__main__":  # quick self-test of the serotonin-vs-DMT discriminator
    cases = {
        "serotonin": "NCCc1c[nH]c2ccc(O)cc12",
        "DMT": "CN(C)CCc1c[nH]c2ccccc12",
        "psilocin": "CN(C)CCc1c[nH]c2ccc(O)cc12",
        "5-MeO-DMT": "COc1ccc2[nH]cc(CCN(C)C)c2c1",
        "mescaline": "NCCc1cc(OC)c(OC)c(OC)c1",
        "LSD": "CCN(CC)C(=O)C1CN(C)C2Cc3c[nH]c4cccc(c34)C2=C1",
        "DOI": "CC(N)Cc1cc(OC)c(I)cc1OC",
        "amphetamine": "CC(N)Cc1ccccc1",
        "donepezil": "COc1cc2c(cc1OC)C(=O)C(CC3CCN(Cc4ccccc4)CC3)C2",
        "caffeine": "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
        # idalopirdine: tryptamine-core 5-HT6 ANTAGONIST (bulky N-benzyl) - must be window-NEG
        "idalopirdine": "FC(F)C(F)(F)COc1cccc(CNCCc2c[nH]c3ccc(F)cc23)c1",
        # sumatriptan: tryptamine 5-HT1B/1D agonist, anti-migraine, very polar (sulfonamide) - NEG
        "sumatriptan": "CNS(=O)(=O)Cc1ccc2[nH]cc(CCN(C)C)c2c1",
    }
    for name, smi in cases.items():
        c = psychoplastogen_window(smi)
        print(f"{name:12} window={str(c.window):5} scaffold={c.scaffold} "
              f"clogP={c.clogp} TPSA={c.tpsa} access={c.intracellular_access}")
