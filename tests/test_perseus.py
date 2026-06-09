"""Tests for PERSEUS v2 - the persistence-aware pro-cognition engine.

RDKit + numpy/pandas. Deterministic unit tests for L1 (CNS gate) and L3 (mechanism
reversibility, mechanism-fired 3-tier), plus a data-gated integration test locking the
two-head behaviour, the L0-mismatch guard, prodrug handling, and abstain-by-default.
"""
from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("rdkit")

from mammal_repurposing.engine.cns_exposure import FAIL, PASS, cns_exposure_gate, structural_vetoes  # noqa: E402
from mammal_repurposing.engine.reversibility import reversibility_call, senolytic_similarity  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"

AMPHETAMINE = "CC(N)Cc1ccccc1"
DONEPEZIL = "COc1cc2c(cc1OC)C(=O)C(CC3CCN(Cc4ccccc4)CC3)C2"
NEOSTIGMINE = "CN(C)C(=O)Oc1cccc([N+](C)(C)C)c1"
VORINOSTAT = "ONC(=O)CCCCCCC(=O)Nc1ccccc1"
FISETIN = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)ccc12"
PIPERLONGUMINE = "COc1cc(/C=C/C(=O)N2CCC=CC2=O)cc(OC)c1OC"
FLUOXETINE = "CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1"


# --------------------------------------------------------------------------
# L1 - CNS-exposure gate
# --------------------------------------------------------------------------

def test_cns_vetoes_and_three_way():
    assert "quaternary_ammonium" in structural_vetoes(NEOSTIGMINE)
    assert cns_exposure_gate(NEOSTIGMINE).verdict == FAIL
    assert cns_exposure_gate(AMPHETAMINE).verdict == PASS
    assert cns_exposure_gate(DONEPEZIL).verdict == PASS
    assert cns_exposure_gate("not_a_smiles").verdict == "ABSTAIN"


# --------------------------------------------------------------------------
# L3 - mechanism reversibility (3 honest tiers, fires from mechanism)
# --------------------------------------------------------------------------

def test_substrate_from_mechanism_axis():
    # the curated mechanism axis substrate is honoured (not the structural class)
    assert reversibility_call("CCO", "transient").substrate == "transient"
    r = reversibility_call(FLUOXETINE, "plasticity_window")
    assert r.substrate == "plasticity_window"


def test_senolytic_is_ablative():
    # a senolytic flavonol is recognised as ablative via structural similarity
    assert senolytic_similarity(FISETIN) >= 0.55
    r = reversibility_call(FISETIN, "transient")
    assert r.substrate == "ablative" and r.self_maintaining is True


def test_reversible_hdaci_is_capable_not_durable():
    # a reversible HDAC inhibitor is state-CAPABLE (flagged) but NOT auto-persistent
    r = reversibility_call(VORINOSTAT, "transient")
    assert r.substrate == "transient"          # not promoted
    assert r.self_maintaining is False
    assert any("hdac" in c for c in r.capability_flags)
    assert r.pulsed_self_maintaining is False   # default: reversible-only


def test_pulsed_hdaci_self_maintaining_hypothesis():
    # Nat Genet 2025 correction: a HDACi curated as self-maintaining (pulsed dosing) earns the
    # pulsed-epigenetic-memory HYPOTHESIS flag - but only with that curated evidence, and it
    # still does NOT auto-promote the substrate to durable.
    r = reversibility_call(VORINOSTAT, "transient", axis_self_maintaining=True)
    assert r.pulsed_self_maintaining is True
    assert r.substrate == "transient" and any("hdac" in c for c in r.capability_flags)
    # a non-HDACi self-maintaining axis call does NOT spuriously set the pulsed flag
    assert reversibility_call("CCO", "transient", axis_self_maintaining=True).pulsed_self_maintaining is False


# --------------------------------------------------------------------------
# integration - the two-head engine
# --------------------------------------------------------------------------

_LEDGERS = [RAW / "clinical_outcomes_ledger.csv", RAW / "clinical_outcomes_ledger_EXTENSION.csv",
            RAW / "clinical_outcomes_ledger_CTGOV.csv", RAW / "clinical_outcomes_ledger_RESEARCH.csv"]
_SMILES = RAW / "ledger_compound_smiles.csv"
_HAVE = (_SMILES.exists() and all(p.exists() for p in _LEDGERS)
         and (RAW / "persistence_axis_classes.csv").exists())


@pytest.fixture(scope="module")
def engine():
    from mammal_repurposing.engine.perseus import PerseusEngine
    return PerseusEngine(_LEDGERS, _SMILES, RAW / "persistence_axis_classes.csv",
                         RAW / "persistence_axis_overrides.csv")


@pytest.mark.skipif(not _HAVE, reason="engine data not present")
def test_two_heads_and_gates(engine):
    from mammal_repurposing.engine.perseus import P_CANDIDATE, P_EXCLUDE_CNS, P_NULL
    # permanent-cation cholinesterase inhibitor: excluded at the CNS gate (both heads)
    r = engine.score("neostigmine", NEOSTIGMINE)
    assert r.cns_verdict == FAIL and r.persistence_verdict == P_EXCLUDE_CNS
    assert r.symptomatic_verdict == "EXCLUDED_NO_CNS"
    # reversible stimulant that routes: real symptomatic tier, NULL persistence, transient
    r2 = engine.score("methylphenidate", "COC(=O)C(C1CCCCN1)c1ccccc1")
    assert r2.cns_verdict == PASS and r2.persistence_verdict == P_NULL
    assert r2.symptomatic_verdict in ("HIGH", "MED", "LOW")
    assert r2.substrate == "transient"
    # CNS-penetrant senolytic: the genuine ablative CANDIDATE path
    r3 = engine.score("piperlongumine", PIPERLONGUMINE)
    assert r3.persistence_verdict == P_CANDIDATE and r3.substrate == "ablative"
    assert r3.self_maintaining and r3.persistence_live


@pytest.mark.skipif(not _HAVE, reason="engine data not present")
def test_psychoplastogen_window_fires_for_permeant_psychedelic(engine):
    from mammal_repurposing.engine.perseus import P_WINDOW
    # psilocybin (phosphate prodrug) -> psilocin (active) -> serotonergic + permeant ->
    # plasticity WINDOW (off-DTI-axis, permeability-gated). The uncurated psychedelic that the
    # pre-L4 engine ABSTAINed on is now correctly a (permissive, not durable) window.
    r = engine.score("psilocybin", "CN(C)CCc1c[nH]c2ccc(OP(=O)(O)O)cc12")
    assert r.persistence_verdict == P_WINDOW and r.persistence_live
    assert any("psychoplastogen_window" in f for f in r.flags)
    # impermeant serotonergic agonist (serotonin) must NOT get the window (permeability gate)
    r2 = engine.score("serotonin", "NCCc1c[nH]c2ccc(O)cc12")
    assert r2.persistence_verdict != P_WINDOW
    assert not any("psychoplastogen_window" in f for f in r2.flags)


@pytest.mark.skipif(not _HAVE, reason="engine data not present")
def test_l4b_nmda_router_separates_ketamine_from_memantine(engine):
    from mammal_repurposing.engine.perseus import P_TESTED_NEG, P_WINDOW
    # L4b pre-registered ablation at the engine level: ketamine and memantine are
    # descriptor-indistinguishable (clogP 2.90/2.69, TPSA 29.1/26.0), so only the curated NMDA
    # trapping-kinetics table separates them. ketamine (resting-block trapper) -> plasticity
    # WINDOW; memantine (spares the resting pool) -> TESTED_NEGATIVE.
    rk = engine.score("ketamine", "CNC1(c2ccccc2Cl)CCCCC1=O")
    assert rk.persistence_verdict == P_WINDOW and rk.persistence_live
    assert any("nmda_router:WINDOW" in f for f in rk.flags)
    rm = engine.score("memantine", "CC12CC3CC(C)(C1)CC(N)(C3)C2")
    assert rm.persistence_verdict == P_TESTED_NEG and not rm.persistence_live
    assert any("nmda_router:NEGATIVE" in f for f in rm.flags)


@pytest.mark.skipif(not _HAVE, reason="engine data not present")
def test_l0_mismatch_guard_abstains_symptomatic(engine):
    from mammal_repurposing.engine.perseus import P_WINDOW
    # fluoxetine (an SSRI) is misrouted to catecholaminergic by scaffold; the symptomatic
    # head must WITHHOLD the wrong-class prior rather than emit it
    r = engine.score("fluoxetine", FLUOXETINE)
    assert r.symptomatic_verdict == "ABSTAIN"
    assert any("mismatch" in x for x in r.abstain_reasons)
    assert r.persistence_verdict == P_WINDOW   # plasticity_gated mechanism preserved


@pytest.mark.skipif(not _HAVE, reason="engine data not present")
def test_reversible_hdaci_compound_abstains(engine):
    # vorinostat: state-capable but reversible -> ABSTAIN, not a candidate
    from mammal_repurposing.engine.perseus import P_ABSTAIN
    r = engine.score("vorinostat", VORINOSTAT)
    assert r.persistence_verdict == P_ABSTAIN
    assert any("capable" in f for f in r.flags)


@pytest.mark.skipif(not _HAVE, reason="engine data not present")
def test_abstains_by_default_no_demonstrated(engine):
    import pandas as pd
    from mammal_repurposing.engine.perseus import score_frame
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
    assert (scored["persistence_verdict"] == "DEMONSTRATED_HEALTHY").sum() == 0
    assert scored["persistence_live"].sum() < len(scored) / 2
    assert (scored["persistence_verdict"] == "EXCLUDE_NO_CNS").sum() >= 1
    # evidence_design is de-broadcast: anorectics get class_extrapolation, not a borrowed tier
    anorectics = scored[scored["compound"].str.lower().isin(
        ["benzphetamine", "fenproporex", "clobenzorex", "mefenorex"])]
    assert (anorectics["evidence_design"] == "class_extrapolation").all()


@pytest.mark.skipif(not _HAVE, reason="engine data not present")
def test_negative_controls_zero_persistence_false_positives(engine):
    """Specificity: no negative control gets a DURABILITY verdict (the eval headline)."""
    import pandas as pd
    from mammal_repurposing.engine.perseus import (
        P_CANDIDATE, P_DEMONSTRATED, P_DISEASE_MOD, score_frame,
    )
    neg = RAW / "perseus_negative_controls.csv"
    if not neg.exists():
        pytest.skip("negative-control panel not present")
    scored = score_frame(engine, pd.read_csv(neg), dedup_salts=False)
    durability = {P_CANDIDATE, P_DISEASE_MOD, P_DEMONSTRATED}
    assert scored["persistence_verdict"].isin(durability).sum() == 0


@pytest.mark.skipif(not _HAVE, reason="engine data not present")
def test_ground_truth_zero_over_claims(engine):
    """Bidirectional: PERSEUS never asserts more durability than the trial-design label."""
    import pandas as pd
    from mammal_repurposing.engine.perseus import score_frame
    from mammal_repurposing.validation.persistence_eval import over_claims
    gt_path = RAW / "persistence_ground_truth.csv"
    if not gt_path.exists():
        pytest.skip("ground-truth ledger not present")
    gt = pd.read_csv(gt_path)
    gt = gt[gt["structure_scoreable"] == "yes"].copy()
    smi = pd.read_csv(RAW / "ledger_compound_smiles.csv")
    look = dict(zip(smi["compound"].str.lower().str.strip(), smi["smiles"]))
    look.setdefault("methylphenidate", "COC(=O)C(C1CCCCN1)c1ccccc1")
    gt["smiles"] = gt["compound"].str.lower().str.strip().map(look)
    gt = gt[gt["smiles"].notna()]
    scored = score_frame(engine, gt.rename(columns={"compound": "query_id"}), dedup_salts=False)
    scored = scored.merge(gt[["compound", "persistence_label"]], on="compound", how="left")
    oc = [r["compound"] for _, r in scored.iterrows()
          if over_claims(r["persistence_verdict"], r["persistence_label"])]
    assert oc == [], f"PERSEUS over-claimed durability on: {oc}"
