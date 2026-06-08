"""Tests for PERSEUS - the persistence-aware pro-cognition engine (L1 CNS gate,
L3 mechanism reversibility, two-head orchestrator).

RDKit + numpy/pandas. Deterministic unit tests for L1/L3 on literal SMILES, plus a
data-gated integration test that locks the control-panel behaviour.
"""
from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("rdkit")

from mammal_repurposing.engine.cns_exposure import (  # noqa: E402
    FAIL, PASS, cns_exposure_gate, structural_vetoes,
)
from mammal_repurposing.engine.reversibility import (  # noqa: E402
    load_structural_alerts, reversibility_call,
)

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"

# literal SMILES (canonical)
AMPHETAMINE = "CC(N)Cc1ccccc1"
DONEPEZIL = "COc1cc2c(cc1OC)C(=O)C(CC3CCN(Cc4ccccc4)CC3)C2"
TETRAMETHYLAMMONIUM = "C[N+](C)(C)C"
NEOSTIGMINE = "CN(C)C(=O)Oc1cccc([N+](C)(C)C)c1"
VORINOSTAT = "ONC(=O)CCCCCCC(=O)Nc1ccccc1"
DMF = "COC(=O)/C=C/C(=O)OC"


# --------------------------------------------------------------------------
# L1 - CNS-exposure gate
# --------------------------------------------------------------------------

def test_cns_vetoes_permanent_charge():
    assert "quaternary_ammonium" in structural_vetoes(NEOSTIGMINE)
    assert "quaternary_ammonium" in structural_vetoes(TETRAMETHYLAMMONIUM)
    assert structural_vetoes(AMPHETAMINE) == []


def test_cns_gate_three_way():
    # permanent cation -> FAIL (no passive CNS entry)
    assert cns_exposure_gate(NEOSTIGMINE).verdict == FAIL
    # classic CNS-penetrant amines -> PASS
    assert cns_exposure_gate(AMPHETAMINE).verdict == PASS
    assert cns_exposure_gate(DONEPEZIL).verdict == PASS
    # unparseable -> ABSTAIN, never a silent pass
    assert cns_exposure_gate("not_a_smiles").verdict == "ABSTAIN"


def test_cns_gate_low_tpsa_passes():
    # low TPSA is GOOD for crossing - must not be penalised into ABSTAIN
    e = cns_exposure_gate("CNC(C)Cc1ccccc1")  # methamphetamine-like
    assert e.verdict == PASS


# --------------------------------------------------------------------------
# L3 - mechanism reversibility (state vs tone)
# --------------------------------------------------------------------------

@pytest.fixture
def synth_classes():
    return {
        "catecholaminergic_ADHD": {"substrate": "transient_signaling", "rank": 0,
                                   "basis": "reuptake tone"},
        "HDAC_inhibitor": {"substrate": "self_propagating_epigenetic", "rank": 3,
                           "basis": "epigenetic"},
    }


_HAVE_ALERTS = (RAW / "persistence_structural_alerts.csv").exists()


@pytest.mark.skipif(not _HAVE_ALERTS, reason="alerts CSV not present")
def test_tone_class_is_not_state_changing(synth_classes):
    alerts = load_structural_alerts(RAW / "persistence_structural_alerts.csv")
    r = reversibility_call(AMPHETAMINE, "catecholaminergic_ADHD", synth_classes, alerts)
    assert r.substrate_class == "transient_signaling"
    assert r.state_changing is False


@pytest.mark.skipif(not _HAVE_ALERTS, reason="alerts CSV not present")
def test_hdaci_structural_alert_is_state_changing(synth_classes):
    alerts = load_structural_alerts(RAW / "persistence_structural_alerts.csv")
    # vorinostat routed to a TONE class still reads state-changing via its ZBG alert
    r = reversibility_call(VORINOSTAT, "catecholaminergic_ADHD", synth_classes, alerts)
    assert r.substrate_class == "self_propagating_epigenetic"
    assert r.state_changing is True
    assert "hydroxamate_ZBG" in r.alerts
    # NRF2 electrophile -> durable transcriptional
    r2 = reversibility_call(DMF, None, synth_classes, alerts)
    assert r2.substrate_class == "durable_transcriptional" and r2.state_changing


@pytest.mark.skipif(not _HAVE_ALERTS, reason="alerts CSV not present")
def test_alert_outranks_class(synth_classes):
    # max(class_rank, alert_rank): an HDACi in a tone class is upgraded to epigenetic
    alerts = load_structural_alerts(RAW / "persistence_structural_alerts.csv")
    r = reversibility_call(VORINOSTAT, "catecholaminergic_ADHD", synth_classes, alerts)
    assert r.substrate_rank == 3 and r.source.startswith("alert:")


# --------------------------------------------------------------------------
# integration - the two-head engine + control-panel behaviour
# --------------------------------------------------------------------------

_LEDGERS = [RAW / "clinical_outcomes_ledger.csv", RAW / "clinical_outcomes_ledger_EXTENSION.csv",
            RAW / "clinical_outcomes_ledger_CTGOV.csv", RAW / "clinical_outcomes_ledger_RESEARCH.csv"]
_SMILES = RAW / "ledger_compound_smiles.csv"
_HAVE_ENGINE = (_SMILES.exists() and all(p.exists() for p in _LEDGERS)
                and (RAW / "persistence_axis_classes.csv").exists()
                and _HAVE_ALERTS)


@pytest.fixture(scope="module")
def engine():
    from mammal_repurposing.engine.perseus import PerseusEngine
    return PerseusEngine(_LEDGERS, _SMILES, RAW / "persistence_axis_classes.csv",
                         RAW / "persistence_axis_overrides.csv",
                         RAW / "persistence_substrate_classes.csv",
                         RAW / "persistence_structural_alerts.csv")


@pytest.mark.skipif(not _HAVE_ENGINE, reason="engine data not present")
def test_perseus_two_heads_are_orthogonal(engine):
    from mammal_repurposing.engine.perseus import P_CANDIDATE, P_EXCLUDE_CNS, P_NULL
    # a permanent-cation cholinesterase inhibitor: excluded at the CNS gate
    r = engine.score("neostigmine", NEOSTIGMINE)
    assert r.cns_verdict == FAIL and r.persistence_verdict == P_EXCLUDE_CNS
    assert r.symptomatic_verdict == "EXCLUDED_NO_CNS"
    # a reversible stimulant that routes: real symptomatic tier, NULL persistence
    r2 = engine.score("methylphenidate", "COC(=O)C(C1CCCCN1)c1ccccc1")
    assert r2.cns_verdict == PASS and r2.persistence_verdict == P_NULL
    assert r2.symptomatic_verdict in ("HIGH", "MED", "LOW")
    # a state-changing HDACi: a mechanistic persistence CANDIDATE (hypothesis, not proof)
    r3 = engine.score("vorinostat", VORINOSTAT)
    assert r3.persistence_verdict == P_CANDIDATE and r3.state_changing
    assert r3.persistence_live


@pytest.mark.skipif(not _HAVE_ENGINE, reason="engine data not present")
def test_perseus_abstains_by_default_no_demonstrated(engine):
    from mammal_repurposing.engine.perseus import score_frame
    import pandas as pd
    short = ROOT / "reports" / "pipeline" / "f2_catalogue_shortlist.csv"
    if not short.exists():
        pytest.skip("F2 shortlist not present")
    cat = pd.read_csv(RAW / "chembl_approved_catalogue.csv")
    cat["k"] = cat["name"].str.lower().str.strip()
    csmi = dict(zip(cat["k"], cat["smiles"]))
    df = pd.read_csv(short)
    df["smiles"] = df["query_id"].str.lower().str.strip().map(csmi)
    df = df[df["smiles"].notna()]
    scored = score_frame(engine, df)
    # the empty gold standard: nothing is called durable enhancement in healthy people
    assert (scored["persistence_verdict"] == "DEMONSTRATED_HEALTHY").sum() == 0
    # live persistence threads are the minority; CNS gate excludes the misroutes
    assert scored["persistence_live"].sum() < len(scored) / 2
    assert (scored["persistence_verdict"] == "EXCLUDE_NO_CNS").sum() >= 1
