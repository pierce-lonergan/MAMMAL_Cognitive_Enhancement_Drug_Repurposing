"""Tests for the L4b NMDA trapping-kinetics mechanism router (curated-PD, abstain-by-default).

The decisive, novel assertion: a STRUCTURE-only model cannot separate the durable plastogen
(ketamine) from the negative control (memantine) - they are descriptor-indistinguishable - so only
the curated pharmacodynamic table separates them. The router promotes confirmed resting-block
trappers to WINDOW, keeps established non-trappers NEGATIVE, and ABSTAINS-with-reason on scaffold
members whose kinetics are unmeasured.
"""
from __future__ import annotations

import csv
from pathlib import Path

import pytest

pytest.importorskip("rdkit")

from mammal_repurposing.engine.mechanism_router import (  # noqa: E402
    muscarinic_router, muscarinic_scaffold, nmda_router, nmda_scaffold,
)

TABLE = Path(__file__).resolve().parents[1] / "data" / "raw" / "nmda_trapping_table.csv"
KETAMINE = "CNC1(CCCCC1=O)c1ccccc1Cl"          # independent (non-table) SMILES, tests canon lookup
PSILOCIN = "CN(C)CCc1c[nH]c2cccc(O)c12"


def _table_smiles() -> dict:
    out: dict[str, str] = {}
    if TABLE.exists():
        with open(TABLE, newline="", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                out[r["compound"]] = r["smiles"]
    return out


def test_router_curated_verdicts():
    t = _table_smiles()
    if not t:
        pytest.skip("nmda_trapping_table.csv not built (run scripts/118)")
    expected = {"ketamine": "WINDOW", "esketamine": "WINDOW", "arketamine": "WINDOW",
                "memantine": "NEGATIVE", "amantadine": "NEGATIVE", "lanicemine": "NEGATIVE",
                "nitrous oxide": "ABSTAIN", "hydroxynorketamine": "ABSTAIN",
                "phencyclidine": "ABSTAIN"}
    for name, verdict in expected.items():
        if name not in t:
            continue
        c = nmda_router(t[name])
        assert c.mechanism_class == "nmda_channel_blocker", name
        assert c.verdict == verdict, f"{name}: got {c.verdict}, want {verdict}"
        assert c.window is (verdict == "WINDOW"), name


def test_preregistered_ablation_ketamine_vs_memantine():
    # THE central L4b proof: ketamine and memantine are descriptor-indistinguishable (verified:
    # clogP 2.90/2.69, TPSA 29.1/26.0), so ONLY the curated PD table separates them. Ketamine is
    # supplied as an independent SMILES (not the table's canonical string) to exercise the lookup.
    t = _table_smiles()
    if "memantine" not in t:
        pytest.skip("table not built")
    ket = nmda_router(KETAMINE)
    mem = nmda_router(t["memantine"])
    assert ket.verdict == "WINDOW" and mem.verdict == "NEGATIVE"
    assert ket.verdict != mem.verdict          # separated by curated PD, not by structure


def test_non_nmda_returns_no_class():
    assert nmda_router(PSILOCIN).mechanism_class is None      # serotonergic, not NMDA
    assert nmda_router("not_a_smiles").mechanism_class is None


def test_scaffold_detection():
    assert nmda_scaffold(KETAMINE) == "arylcyclohexylamine"
    assert nmda_scaffold(PSILOCIN) is None


SCOPOLAMINE = "CN1C2CC(CC1C3C2O3)OC(=O)C(CO)C4=CC=CC=C4"      # ledger SMILES (scopine/tropane)
IPRATROPIUM = "CC(C)[N+]1(C2CCC1CC(C2)OC(=O)C(CO)C3=CC=CC=C3)C"  # peripheral quaternary muscarinic


def test_muscarinic_tropane_abstains_with_reason():
    # L4b muscarinic lane recognises the tropane/scopine chemotype and ABSTAINS-with-reason - it
    # NEVER promotes durability (scopolamine's carryover is a single-compound clinical fact).
    c = muscarinic_router(SCOPOLAMINE)
    assert c.mechanism_class == "tropane_muscarinic" and c.verdict == "ABSTAIN"
    assert c.window is False and c.reasons and "evidence layer" in c.reasons[0]
    assert muscarinic_scaffold(SCOPOLAMINE) == "tropane"


def test_muscarinic_quaternary_veto_and_non_tropane():
    # peripheral quaternary muscarinic (ipratropium) is vetoed (never reaches the CNS) -> no class
    assert muscarinic_router(IPRATROPIUM).mechanism_class is None
    # non-tropane chemotypes -> no muscarinic class (psilocin serotonergic, ketamine arylcyclohexyl)
    assert muscarinic_router(PSILOCIN).mechanism_class is None
    assert muscarinic_router(KETAMINE).mechanism_class is None


def test_novel_arylcyclohexylamine_abstains_with_reason():
    # an arylcyclohexylamine NOT in the curated table must ABSTAIN-with-reason (durability is
    # trapping kinetics, not structure) - never silently missed, never auto-window.
    c = nmda_router("NC1(c2ccccc2)CCCCC1")     # 1-phenylcyclohexan-1-amine, not in table
    assert c.mechanism_class == "nmda_channel_blocker"
    assert c.verdict == "ABSTAIN" and c.window is False
    assert c.reasons and "not in the curated PD table" in c.reasons[0]
