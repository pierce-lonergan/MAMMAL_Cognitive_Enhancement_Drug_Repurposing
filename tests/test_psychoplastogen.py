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
    s = psychoplastogen_window(SER)
    d = psychoplastogen_window(DMT)
    assert s.scaffold == "tryptamine" and d.scaffold == "tryptamine"   # both are tryptamines
    assert s.window is False and s.intracellular_access is False        # serotonin: impermeant
    assert d.window is True and d.intracellular_access is True          # DMT: permeant


def test_classic_psychedelics_are_window_positive():
    for smi in (PSILOCIN, MESCALINE, LSD, DOI):
        assert psychoplastogen_window(smi).window is True


def test_tryptamine_core_nonplastogens_are_window_negative():
    # idalopirdine: bulky N-benzyl amine -> scaffold veto (not a 5-HT2A-agonist pharmacophore)
    ida = psychoplastogen_window(IDALOPIRDINE)
    assert ida.scaffold is None and ida.window is False


def test_triptans_are_vetoed_not_psychoplastogens():
    # critical-audit fix: triptans are tryptamines but 5-HT1B/1D agonists, NOT 5-HT2A
    # psychoplastogens. The triptan pharmacophore (sulfonamide / cyclic carbamate) is vetoed so
    # they never reach the permeability gate (scaffold None, window False).
    sumatriptan = "CNS(=O)(=O)Cc1ccc2[nH]cc(CCN(C)C)c2c1"      # sulfonamide
    zolmitriptan = "O=C1OCC(Cc2ccc3[nH]cc(CCN(C)C)c3c2)N1"     # oxazolidinone (slips TPSA gate)
    for smi in (sumatriptan, zolmitriptan):
        c = psychoplastogen_window(smi)
        assert c.scaffold is None and c.window is False


def test_permeant_endogenous_monoamine_decoys_are_window_negative():
    # specificity stress test: permeant CNS monoamines / decoys that are NOT durable
    # psychoplastogens must stay window-negative (dopamine/melatonin/methamphetamine/venlafaxine)
    for smi in ("NCCc1ccc(O)c(O)c1",                        # dopamine (polar catechol)
                "CC(=O)NCCc1c[nH]c2ccc(OC)cc12",            # melatonin (N-acetyl -> amide veto)
                "CC(NC)Cc1ccccc1",                          # methamphetamine (no OMe -> no scaffold)
                "COc1ccc(C(CN(C)C)C2(O)CCCCC2)cc1"):        # venlafaxine
        assert psychoplastogen_window(smi).window is False


def test_non_serotonergic_have_no_scaffold():
    for smi in ("CC(N)Cc1ccccc1",                                  # amphetamine (no OMe/halo)
                "Cn1cnc2c1c(=O)n(C)c(=O)n2C",                      # caffeine
                "COc1cc2c(cc1OC)C(=O)C(CC3CCN(Cc4ccccc4)CC3)C2"):  # donepezil
        assert serotonergic_scaffold(smi) is None


def test_window_never_auto_durable_caveat():
    c = psychoplastogen_window(DMT)
    assert "paired with experience" in c.caveat and "never auto-durable" in c.caveat
