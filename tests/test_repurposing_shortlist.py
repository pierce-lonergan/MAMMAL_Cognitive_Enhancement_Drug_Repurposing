"""Gap 7 — Tests for the prospective repurposing shortlist (capstone).

Locks: SUCCESS-class restriction, class-prior-dominated scoring, disease-aware
novelty (schizophrenia covers CIAS), the Gap-4 engagement-reliability flag, and
the real-data headline (CIAS surfaces 5-HT1A + M1/M4 repurposing hypotheses;
FXS surfaces PDE4; xanomeline is correctly NOT novel for CIAS).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import pytest

from mammal_repurposing.reporting import repurposing_shortlist as RS
from mammal_repurposing.validation import disease_reframe as D
from mammal_repurposing.validation import retrospective as R

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "data" / "raw" / "clinical_outcomes_ledger.csv"
ANCHORS = ROOT / "data" / "raw" / "modulator_anchors_seed.csv"
GRID = ROOT / "data" / "results" / "v2" / "v6a_grid_expanded.parquet"
PANEL = ROOT / "data" / "interim" / "targets.parquet"


@dataclass
class _Prior:
    mean: float
    sd: float = 0.1
    verdict: str = "SUCCESS"


def _grid(rows):
    return pd.DataFrame(rows, columns=["compound_name", "target_uniprot",
                                       "predicted_pkd", "target_gene"])


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def test_supplement_and_known_maps_nonempty():
    assert len(RS.REPURPOSING_SUPPLEMENT) >= 8
    assert "buspirone" in RS.KNOWN_INDICATIONS
    assert "roflumilast" in RS.KNOWN_INDICATIONS


def test_success_only_filter_and_scoring():
    named = pd.DataFrame({"name": ["drugA", "drugB"],
                          "expected_top_target": ["P_succ", "P_fail"]})
    priors = {"clsS": _Prior(0.40, verdict="SUCCESS"),
              "clsF": _Prior(0.0, verdict="FAILURE")}
    tcmap = {"P_succ": "clsS", "P_fail": "clsF"}
    grid = _grid([("drugA", "P_succ", 8.0, "GA"), ("drugB", "P_fail", 8.0, "GB")])
    cands = RS.build_repurposing_shortlist("AD", named, priors, grid, tcmap)
    assert [c.compound for c in cands] == ["drugA"]      # failure class dropped
    assert cands[0].class_verdict == "SUCCESS"
    assert cands[0].score > 0


def test_disease_aware_novelty_schizophrenia_covers_cias():
    # a drug used for schizophrenia is NOT a novel CIAS repurposing
    named = pd.DataFrame({"name": ["aripiprazole"], "expected_top_target": ["P1"]})
    priors = {"D1_agonist": _Prior(0.40, verdict="SUCCESS")}
    tcmap = {"P1": "D1_agonist"}
    grid = _grid([("aripiprazole", "P1", 8.0, "DRD1")])
    cands = RS.build_repurposing_shortlist("CIAS", named, priors, grid, tcmap)
    assert cands[0].novel is False        # schizophrenia covers CIAS


def test_engagement_reliability_flag():
    named = pd.DataFrame({"name": ["x"], "expected_top_target": ["Q01959"]})  # SLC6A3 transporter
    priors = {"catecholaminergic": _Prior(0.5, verdict="SUCCESS")}
    tcmap = {"Q01959": "catecholaminergic"}
    grid = _grid([("x", "Q01959", 8.0, "SLC6A3")])
    cands = RS.build_repurposing_shortlist("ADHD", named, priors, grid, tcmap)
    assert cands[0].engagement_reliable is False    # transporter — Gap-4 flag


# ---------------------------------------------------------------------------
# Real-data headline
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not (LEDGER.exists() and GRID.exists() and PANEL.exists()),
                    reason="ledger / grid / panel absent")
def test_real_cias_surfaces_muscarinic_and_5ht1a():
    ledger = R.load_clinical_ledger(LEDGER)
    anchors = pd.read_csv(ANCHORS, comment="#")
    ev = D.load_disease_evidence(ledger, anchors)
    priors = D.build_disease_class_priors("CIAS", ev)
    comp = pd.read_parquet(ROOT / "data" / "interim" / "compounds.parquet")
    named = comp[comp.evidence_tier.notna() & comp.expected_top_target.notna()]
    named = RS.named_drugs_with_supplement(named)
    grid = pd.read_parquet(GRID)
    cands = RS.build_repurposing_shortlist(
        "CIAS", named, priors, grid, D.disease_target_class_map())
    classes = {c.mechanism_class for c in cands}
    assert "M1_M4_agonist" in classes          # xanomeline's class now scorable
    assert "5HT1A_partial_agonist" in classes
    # xanomeline (approved for CIAS) must NOT be flagged novel
    xa = [c for c in cands if c.compound.lower() == "xanomeline"]
    if xa:
        assert xa[0].novel is False
    # there is at least one genuine novel hypothesis
    assert any(c.novel for c in cands)


@pytest.mark.skipif(not (LEDGER.exists() and GRID.exists() and PANEL.exists()),
                    reason="artifacts absent")
def test_real_fxs_surfaces_pde4_repurposing():
    ledger = R.load_clinical_ledger(LEDGER)
    anchors = pd.read_csv(ANCHORS, comment="#")
    ev = D.load_disease_evidence(ledger, anchors)
    priors = D.build_disease_class_priors("FXS", ev)
    comp = pd.read_parquet(ROOT / "data" / "interim" / "compounds.parquet")
    named = comp[comp.evidence_tier.notna() & comp.expected_top_target.notna()]
    named = RS.named_drugs_with_supplement(named)
    grid = pd.read_parquet(GRID)
    cands = RS.build_repurposing_shortlist(
        "FXS", named, priors, grid, D.disease_target_class_map())
    novel_names = {c.compound.lower() for c in cands if c.novel}
    # roflumilast (approved COPD) is the marquee PDE4 repurposing hypothesis
    assert "roflumilast" in novel_names or "rolipram" in novel_names
    assert all(c.class_verdict == "SUCCESS" for c in cands)
