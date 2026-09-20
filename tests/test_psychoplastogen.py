"""Tests for the L4 psychoplastogen plasticity-window detector.

The decisive, novel assertion (Vargas 2023): serotonin and DMT are near-ISOAFFINE at 5-HT2A,
yet serotonin is window-NEGATIVE and DMT window-POSITIVE - the discriminator is membrane
PERMEABILITY (TPSA/HBD), a structure-computable property OFF the DTI axis. Plus false-positive
guards (idalopirdine, sumatriptan) that share a tryptamine core but are not plastogenic agonists.
"""
from __future__ import annotations

import pytest

pytest.importorskip("rdkit")

from mammal_repurposing.engine.psychoplastogen import (  # noqa: E402
    L4_VALIDATED_ASSAY, WINDOW_ASSAY_FAMILIES,
    psychoplastogen_window, serotonergic_scaffold,
)

SER = "NCCc1c[nH]c2ccc(O)cc12"          # serotonin (isoaffine 5-HT2A, impermeant)
DMT = "CN(C)CCc1c[nH]c2ccccc12"          # DMT (isoaffine 5-HT2A, permeant)
PSILOCIN = "CN(C)CCc1c[nH]c2cccc(O)c12"
MESCALINE = "NCCc1cc(OC)c(OC)c(OC)c1"
LSD = "CCN(CC)C(=O)C1CN(C)C2Cc3c[nH]c4cccc(c34)C2=C1"
DOI = "CC(N)Cc1cc(OC)c(I)cc1OC"
# tryptamine-core NON-plastogens (must be window-negative)
IDALOPIRDINE = "FC(F)C(F)(F)COc1cccc(CNCCc2c[nH]c3ccc(F)cc23)c1"  # 5-HT6 antagonist, N-benzyl
SUMATRIPTAN = "CNS(=O)(=O)Cc1ccc2[nH]cc(CCN(C)C)c2c1"            # 5-HT1 agonist, sulfonamide-polar


def test_serotonin_vs_dmt_isoaffine_but_permeability_discriminates():
    # THE central thesis: same 5-HT2A binder family, opposite window call, on permeability alone
    s = psychoplastogen_window(SER, assay=L4_VALIDATED_ASSAY)
    d = psychoplastogen_window(DMT, assay=L4_VALIDATED_ASSAY)
    assert s.scaffold == "tryptamine" and d.scaffold == "tryptamine"   # both are tryptamines
    assert s.window is False and s.intracellular_access is False        # serotonin: impermeant
    assert d.window is True and d.intracellular_access is True          # DMT: permeant


def test_classic_psychedelics_are_window_positive():
    for smi in (PSILOCIN, MESCALINE, LSD, DOI):
        assert psychoplastogen_window(smi, assay=L4_VALIDATED_ASSAY).window is True


def test_tryptamine_core_nonplastogens_are_window_negative():
    # idalopirdine: bulky N-benzyl amine -> scaffold veto (not a 5-HT2A-agonist pharmacophore)
    ida = psychoplastogen_window(IDALOPIRDINE, assay=L4_VALIDATED_ASSAY)
    assert ida.scaffold is None and ida.window is False


def test_triptans_are_vetoed_not_psychoplastogens():
    # critical-audit fix: triptans are tryptamines but 5-HT1B/1D agonists, NOT 5-HT2A
    # psychoplastogens. The triptan pharmacophore (sulfonamide / cyclic carbamate) is vetoed so
    # they never reach the permeability gate (scaffold None, window False).
    sumatriptan = "CNS(=O)(=O)Cc1ccc2[nH]cc(CCN(C)C)c2c1"      # sulfonamide
    zolmitriptan = "O=C1OCC(Cc2ccc3[nH]cc(CCN(C)C)c3c2)N1"     # oxazolidinone (slips TPSA gate)
    for smi in (sumatriptan, zolmitriptan):
        c = psychoplastogen_window(smi, assay=L4_VALIDATED_ASSAY)
        assert c.scaffold is None and c.window is False


def test_permeant_endogenous_monoamine_decoys_are_window_negative():
    # specificity stress test: permeant CNS monoamines / decoys that are NOT durable
    # psychoplastogens must stay window-negative (dopamine/melatonin/methamphetamine/venlafaxine)
    for smi in ("NCCc1ccc(O)c(O)c1",                        # dopamine (polar catechol)
                "CC(=O)NCCc1c[nH]c2ccc(OC)cc12",            # melatonin (N-acetyl -> amide veto)
                "CC(NC)Cc1ccccc1",                          # methamphetamine (no OMe -> no scaffold)
                "COc1ccc(C(CN(C)C)C2(O)CCCCC2)cc1"):        # venlafaxine
        assert psychoplastogen_window(smi, assay=L4_VALIDATED_ASSAY).window is False


def test_non_serotonergic_have_no_scaffold():
    for smi in ("CC(N)Cc1ccccc1",                                  # amphetamine (no OMe/halo)
                "Cn1cnc2c1c(=O)n(C)c(=O)n2C",                      # caffeine
                "COc1cc2c(cc1OC)C(=O)C(CC3CCN(Cc4ccccc4)CC3)C2"):  # donepezil
        assert serotonergic_scaffold(smi) is None


def test_triptan_veto_covers_sulfone_and_triazole_d6():
    # D6 decoy-scan fix: the D4 veto caught only sulfonamide/carbamate triptans; the broader scan
    # found rizatriptan (1,2,4-triazolylmethyl) and eletriptan (phenyl SULFONE, not sulfonamide)
    # slipping through. The broadened sulfonyl + triazole vetoes catch them. SMILES from PubChem.
    rizatriptan = "CN(C)CCC1=CNC2=C1C=C(C=C2)CN3C=NC=N3"               # triazolylmethyl
    eletriptan = "CN1CCCC1CC2=CNC3=C2C=C(C=C3)CCS(=O)(=O)C4=CC=CC=C4"  # phenyl sulfone
    for smi in (rizatriptan, eletriptan):
        c = psychoplastogen_window(smi, assay=L4_VALIDATED_ASSAY)
        assert c.scaffold is None and c.window is False


def test_pergolide_thioether_veto_d5():
    # D5 adversarial-audit fix: pergolide is a dopaminergic CLAVINE ergoline, NOT a 5-HT2A
    # psychoplastogen. Its aliphatic thioether (CH2-S-CH3) is the veto signal - verified absent
    # from every compound in the persistence positive ledger, so the veto is collateral-free.
    # SMILES from PubChem.
    pergolide = "CCCN1CC(CC2C1CC3=CNC4=CC=CC2=C34)CSC"
    c = psychoplastogen_window(pergolide, assay=L4_VALIDATED_ASSAY)
    assert c.scaffold is None and c.window is False


def test_isodmt_recall_d5():
    # D5 recall innovation: isoDMTs put the aminoethyl on the indole N1 (C3 free), so the
    # C3-requiring tryptamine SMARTS misses them although they are genuine psychoplastogens
    # (Dunlap/Olson). zalsupindole / AAZ-A-154 is the ledger exemplar (SMILES from the ledger).
    zalsupindole = "CC(CN1C=CC2=C1C=CC(=C2)OC)N(C)C"
    c = psychoplastogen_window(zalsupindole, assay=L4_VALIDATED_ASSAY)
    assert c.scaffold == "isodmt" and c.window is True
    # N1-substitution is therefore a POSITIVE signal and must never be vetoed (it would suppress
    # the whole isoDMT class) - which is exactly why methysergide stays an accepted, documented FP.


def test_window_never_auto_durable_caveat():
    c = psychoplastogen_window(DMT, assay=L4_VALIDATED_ASSAY)
    assert "paired with experience" in c.caveat and "never auto-durable" in c.caveat


# --- B2: the assay key -------------------------------------------------------------------------
# These four FAIL on pre-B2 code (where `assay` did not exist and the verdict was a bare
# compound-level claim) and PASS after. The defect they pin: Sheynin 2019 (PMID 30766471) measured
# perceptual learning and ocular dominance in the SAME healthy adults on donepezil and got opposite
# signs, so an unscoped "opens a plasticity window" flag asserts something the data does not support.

def test_assay_is_required_and_has_no_default():
    """An unscoped call must not be expressible. A default would let every existing call site keep
    making the unscoped claim while appearing to have been fixed."""
    import pytest
    with pytest.raises(TypeError):
        psychoplastogen_window(DMT)          # no assay= -> must not be callable


def test_unknown_assay_raises_rather_than_silently_scoping():
    import pytest
    with pytest.raises(ValueError, match="unknown assay family"):
        psychoplastogen_window(DMT, assay="ocular-dominance")   # hyphen typo, not a family


def test_verdict_is_structural_but_standing_is_assay_relative():
    """Same structure in, same window out - the assay never changes the chemistry. What it changes
    is whether the verdict is VALIDATED or EXTRAPOLATED, which is the thing a reader needs."""
    val = psychoplastogen_window(DMT, assay=L4_VALIDATED_ASSAY)
    ext = psychoplastogen_window(DMT, assay="ocular_dominance")
    assert val.window is ext.window is True
    assert val.assay_evidence == "validated"
    assert ext.assay_evidence == "extrapolated"
    assert ext.assay == "ocular_dominance"


def test_extrapolated_positive_carries_the_scope_caveat_in_reasons():
    ext = psychoplastogen_window(DMT, assay="tms_ltp")
    assert any("EXTRAPOLATED" in r and "30766471" in r for r in ext.reasons), ext.reasons
    val = psychoplastogen_window(DMT, assay=L4_VALIDATED_ASSAY)
    assert not any("EXTRAPOLATED" in r for r in val.reasons)


def test_exactly_one_family_is_validated():
    """If a second family ever earns `validated`, it must be because scripts/125 measured it - not
    because someone widened a constant. This test is the tripwire on that constant."""
    assert L4_VALIDATED_ASSAY in WINDOW_ASSAY_FAMILIES
    assert L4_VALIDATED_ASSAY == "dendritic_spine"
