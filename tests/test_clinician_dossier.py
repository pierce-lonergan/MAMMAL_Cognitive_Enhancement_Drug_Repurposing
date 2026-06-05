"""Gap 5 — Tests for the clinician evidence dossier.

Locks: GRADE up/down-grading logic (incl. the rule that binding-reliability
only matters for model-predicted effects, not a drug's own trial), liability
flagging, and the real-data headline (donepezil HIGH/SUCCESS; idalopirdine
correctly flagged FAILURE with a warning).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mammal_repurposing.reporting import clinician_dossier as Dx
from mammal_repurposing.validation import retrospective as R
from mammal_repurposing.validation import disease_reframe as D

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "data" / "raw" / "clinical_outcomes_ledger.csv"
ANCHORS = ROOT / "data" / "raw" / "modulator_anchors_seed.csv"


# ---------------------------------------------------------------------------
# GRADE logic
# ---------------------------------------------------------------------------

def test_grade_high_for_strong_rct_base():
    g, reasons = Dx.grade_evidence(
        k_rcts=18, n_class_drugs=3, class_sd=0.08, own_trial=True,
        disease_direct=True, ci_width=0.17, binding_reliable=True)
    assert g == "HIGH"


def test_grade_downgrades_stack():
    # no trials + imprecision + indirectness -> low/very low
    g, reasons = Dx.grade_evidence(
        k_rcts=0, n_class_drugs=1, class_sd=0.05, own_trial=False,
        disease_direct=False, ci_width=0.6, binding_reliable=True)
    assert g in ("LOW", "VERY LOW")
    assert len(reasons) >= 3


def test_binding_reliability_ignored_for_own_trial():
    # an approved drug's own pivotal does not depend on MAMMAL binding
    g_own, _ = Dx.grade_evidence(
        k_rcts=18, n_class_drugs=3, class_sd=0.08, own_trial=True,
        disease_direct=True, ci_width=0.17, binding_reliable=False)
    assert g_own == "HIGH"
    # but a model prediction at an allosteric-blind target IS down-graded
    g_pred, reasons = Dx.grade_evidence(
        k_rcts=18, n_class_drugs=3, class_sd=0.08, own_trial=False,
        disease_direct=True, ci_width=0.17, binding_reliable=False)
    assert g_pred != "HIGH"
    assert any("allosteric-blind" in r for r in reasons)


def test_liability_concern_map_nonempty():
    assert "DRD2" in Dx.LIABILITY_CONCERNS
    assert "HTR2A" in Dx.LIABILITY_CONCERNS
    assert len(Dx.LIABILITY_CONCERNS) >= 15


def test_liability_flags_filter_offpanel_and_threshold():
    liab = pd.DataFrame({
        "compound_name": ["x", "x", "x", "x"],
        "target_uniprot": ["P_on", "DRD2u", "low", "P_panel"],
        "target_gene": ["ONPANEL", "DRD2", "HTR2A", "PANELG"],
        "predicted_pkd": [9.0, 7.5, 5.0, 9.9],
    })
    flags = Dx._liability_flags("x", liab, panel_uniprots={"P_on", "P_panel"},
                                pkd_threshold=6.5)
    genes = [f[0] for f in flags]
    assert "DRD2" in genes            # off-panel, above threshold
    assert "ONPANEL" not in genes     # on cognition panel -> not a liability
    assert "HTR2A" not in genes       # below threshold


# ---------------------------------------------------------------------------
# Real-data headline
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not (LEDGER.exists() and ANCHORS.exists()),
                    reason="ledger / anchors absent")
def test_donepezil_dossier_high_success():
    ledger = R.load_clinical_ledger(LEDGER)
    anchors = pd.read_csv(ANCHORS, comment="#")
    ev = D.load_disease_evidence(ledger, anchors)
    priors = D.build_disease_class_priors("AD", ev)
    lrow = ledger[(ledger.compound == "donepezil")].iloc[0]
    card = Dx.build_dossier("donepezil", "AD", ledger=ledger,
                            disease_priors=priors, ledger_row=lrow,
                            target_gene="ACHE", mechanism_class="AChE_inhibitor")
    assert card.grade == "HIGH"
    assert card.class_verdict == "SUCCESS"
    assert 0.30 <= card.g <= 0.42
    assert card.own_trial
    md = Dx.render_card_md(card)
    assert "donepezil" in md and ("GRADE" in md.upper() or "HIGH" in md.upper())


@pytest.mark.skipif(not (LEDGER.exists() and ANCHORS.exists()),
                    reason="ledger / anchors absent")
def test_idalopirdine_dossier_flags_failure():
    ledger = R.load_clinical_ledger(LEDGER)
    anchors = pd.read_csv(ANCHORS, comment="#")
    ev = D.load_disease_evidence(ledger, anchors)
    priors = D.build_disease_class_priors("AD", ev)
    lrow = ledger[(ledger.compound == "idalopirdine")].iloc[0]
    card = Dx.build_dossier("idalopirdine", "AD", ledger=ledger,
                            disease_priors=priors, ledger_row=lrow,
                            target_gene="HTR6", mechanism_class="5HT6_antagonist")
    assert card.class_verdict == "FAILURE"
    # the dossier must warn about the negative class track record
    assert any("NEGATIVE pivotal" in c or "failed class" in c for c in card.caveats)
