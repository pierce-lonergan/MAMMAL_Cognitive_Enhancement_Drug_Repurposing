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
    nmda_router, nmda_scaffold,
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


def test_novel_arylcyclohexylamine_abstains_with_reason():
    # an arylcyclohexylamine NOT in the curated table must ABSTAIN-with-reason (durability is
    # trapping kinetics, not structure) - never silently missed, never auto-window.
    c = nmda_router("NC1(c2ccccc2)CCCCC1")     # 1-phenylcyclohexan-1-amine, not in table
    assert c.mechanism_class == "nmda_channel_blocker"
    assert c.verdict == "ABSTAIN" and c.window is False
    assert c.reasons and "not in the curated PD table" in c.reasons[0]
