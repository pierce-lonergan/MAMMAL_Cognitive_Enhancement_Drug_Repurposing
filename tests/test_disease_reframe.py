"""Gap 2 — Tests for the disease-population reframe.

Locks: disease-bucket assignment (AD ≠ ADHD), clean target→mechanism mapping
(fixes the α7/AChE lump), k-weighted disease class priors, the fallback prior
for evidence-free classes, the disease-specific anchor override, within-disease
class-LOCO validation, and the real-data headline (each disease surfaces its
own winning mechanism; SUCCESS classes out-rank FAILURE classes).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mammal_repurposing.validation import disease_reframe as D
from mammal_repurposing.validation import retrospective as R

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "data" / "raw" / "clinical_outcomes_ledger.csv"
ANCHORS = ROOT / "data" / "raw" / "modulator_anchors_seed.csv"
V6B = ROOT / "data" / "results" / "v2" / "cluster_d_posterior_expanded_v2_mh8_ta99.parquet"


# ---------------------------------------------------------------------------
# Disease bucketing
# ---------------------------------------------------------------------------

def test_ad_does_not_match_adhd():
    assert D.disease_match("AD", "AD")
    assert D.disease_match("AD-mod-sev", "AD")
    assert D.disease_match("AD/schizophrenia", "AD")
    assert not D.disease_match("ADHD", "AD")          # the dangerous false friend
    assert D.disease_match("ADHD", "ADHD")


def test_cias_matches_schizophrenia():
    assert D.disease_match("CIAS-schizophrenia", "CIAS")
    assert D.disease_match("schizophrenia", "CIAS")
    assert D.disease_match("schizophrenia/ADHD", "CIAS")
    assert not D.disease_match("AD", "CIAS")


def test_multi_indication_buckets():
    b = D.buckets_for("AD/schizophrenia")
    assert "AD" in b and "CIAS" in b


def test_fxs_bucket():
    assert D.disease_match("FXS", "FXS")
    assert "FXS" in D.buckets_for("FXS")


# ---------------------------------------------------------------------------
# Clean target → mechanism class (fixes the v11 panel lump)
# ---------------------------------------------------------------------------

def test_alpha7_not_lumped_with_ache():
    # CHRNA7 (encenicline class) must NOT be AChE-I, or the cholinesterase
    # prior is contaminated by α7 failures.
    assert D.TARGET_TO_MECHCLASS["P36544"] == "alpha7_nAChR"
    assert D.TARGET_TO_MECHCLASS["P22303"] == "AChE_inhibitor"


def test_nmda_subunits_share_class():
    assert D.TARGET_TO_MECHCLASS["Q13224"] == "NMDA_modulator"
    assert D.TARGET_TO_MECHCLASS["Q12879"] == "NMDA_modulator"


def test_htr1a_mapped_for_expanded_grid():
    # HTR1A entered the grid via the 13->23 expansion; its class key must match
    # the CIAS prior's tandospirone class so grid scoring and prior align.
    assert D.TARGET_TO_MECHCLASS["P08908"] == "5HT1A_partial_agonist"


# ---------------------------------------------------------------------------
# V6.A grid expansion (scripts/77 output)
# ---------------------------------------------------------------------------

EXPANDED_GRID = ROOT / "data" / "results" / "v2" / "v6a_grid_expanded.parquet"


@pytest.mark.skipif(not EXPANDED_GRID.exists(),
                    reason="expanded grid not built (run scripts/77)")
def test_expanded_grid_coverage_and_no_peptides():
    g = pd.read_parquet(EXPANDED_GRID)
    # full 31-target panel once scripts/81 has scored the final 8 (CHRM1/4, HTR6,
    # GRM2/3/5, GlyT1, HTR4); >=23 if only the cached-signal expansion has run.
    assert g["target_uniprot"].nunique() >= 23
    assert {"compound_name", "target_uniprot", "predicted_pkd",
            "binding_source"} <= set(g.columns)
    # out-of-domain peptides/biologics filtered (MAMMAL DTI is small-molecule)
    cmpds = set(g["compound_name"])
    assert "semaglutide" not in cmpds
    assert "liraglutide" not in cmpds
    # targets the original 13-grid lacked
    for u in ("P36544", "Q12879", "Q99720", "P08908"):  # CHRNA7, GRIN2A, SIGMAR1, HTR1A
        assert u in set(g["target_uniprot"].astype(str))


@pytest.mark.skipif(not EXPANDED_GRID.exists(),
                    reason="expanded grid not built (run scripts/77)")
def test_full_panel_muscarinic_and_5ht6_scorable():
    """Once the panel is finished to 31 (scripts/81 MAMMAL scoring), the CIAS
    muscarinic winner (CHRM1/CHRM4) and the AD 5-HT6 failure class (HTR6) must be
    in the grid. Skips gracefully if only the 23-target cached expansion has run."""
    g = pd.read_parquet(EXPANDED_GRID)
    present = set(g["target_uniprot"].astype(str))
    if g["target_uniprot"].nunique() < 31:
        pytest.skip("panel not yet finished to 31 (run scripts/81 in the MAMMAL venv)")
    for u in ("P11229", "P08173", "P50406"):   # CHRM1, CHRM4, HTR6
        assert u in present
    assert D.TARGET_TO_MECHCLASS["P11229"] == "M1_M4_agonist"
    assert D.TARGET_TO_MECHCLASS["P50406"] == "5HT6_antagonist"


# ---------------------------------------------------------------------------
# Disease class priors
# ---------------------------------------------------------------------------

@pytest.fixture
def mini_evidence():
    return pd.DataFrame({
        "compound": ["donepezil", "galantamine", "encenicline", "idalopirdine"],
        "compound_lower": ["donepezil", "galantamine", "encenicline", "idalopirdine"],
        "target_uniprot": ["P22303", "P22303", "P36544", "P50406"],
        "mechanism_class": ["AChE_inhibitor", "AChE_inhibitor",
                            "alpha7_nAChR", "5HT6_antagonist"],
        "g": [0.36, 0.40, 0.0, -0.05],
        "k": [18, 12, 2, 3],
        "disease": ["AD", "AD", "AD", "AD"],
        "outcome": ["SUCCESS", "SUCCESS", "FAILURE", "FAILURE"],
        "source": ["ledger"] * 4,
    })


def test_k_weighted_class_mean(mini_evidence):
    pri = D.build_disease_class_priors("AD", mini_evidence)
    # AChE-I: k-weighted (18*0.36 + 12*0.40)/30 = 0.376
    assert abs(pri["AChE_inhibitor"].mean - (18 * 0.36 + 12 * 0.40) / 30) < 1e-6
    assert pri["AChE_inhibitor"].verdict == "SUCCESS"


def test_failure_classes_flagged(mini_evidence):
    pri = D.build_disease_class_priors("AD", mini_evidence)
    assert pri["alpha7_nAChR"].verdict == "FAILURE"
    assert pri["5HT6_antagonist"].verdict == "FAILURE"
    assert pri["5HT6_antagonist"].mean < 0.0


def test_success_outranks_failure(mini_evidence):
    pri = D.build_disease_class_priors("AD", mini_evidence)
    succ = [p.mean for p in pri.values() if p.verdict == "SUCCESS"]
    fail = [p.mean for p in pri.values() if p.verdict == "FAILURE"]
    assert min(succ) > max(fail)


def test_prior_table_fallback(mini_evidence):
    pri = D.build_disease_class_priors("AD", mini_evidence)
    tbl = D.disease_class_prior_table(pri, all_classes=["AChE_inhibitor", "NEVER_SEEN"])
    assert tbl["AChE_inhibitor"]["mean"] > 0.3
    assert tbl["NEVER_SEEN"] == D.FALLBACK_PRIOR     # evidence-free → weak prior


# ---------------------------------------------------------------------------
# Disease ceiling
# ---------------------------------------------------------------------------

def test_disease_ceiling_above_healthy():
    assert D.DISEASE_CEILING["healthy"] == 0.50
    assert D.DISEASE_CEILING["AD"] > 0.50
    assert D.DISEASE_CEILING["FXS"] > D.DISEASE_CEILING["AD"]   # zatolmilast 0.71


def test_diversified_shortlist_caps_per_class():
    grid = pd.DataFrame({
        "mechanism_class": ["A"] * 5 + ["B"] * 5,
        "g_predicted": [0.5, 0.4, 0.3, 0.2, 0.1, 0.45, 0.35, 0.25, 0.15, 0.05],
        "roberts_ceiling_ok": [True] * 10,
    })
    dv = D.diversified_shortlist(grid, per_class=2, n=10)
    assert (dv["mechanism_class"] == "A").sum() == 2
    assert (dv["mechanism_class"] == "B").sum() == 2


# ---------------------------------------------------------------------------
# Real-data integration
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not (LEDGER.exists() and ANCHORS.exists()),
                    reason="ledger / anchors not present")
def test_real_ad_prior_cholinergic_top():
    ledger = R.load_clinical_ledger(LEDGER)
    anchors = pd.read_csv(ANCHORS, comment="#")
    ev = D.load_disease_evidence(ledger, anchors)
    pri = D.build_disease_class_priors("AD", ev)
    # AChE-I is the AD standard of care → highest-evidenced SUCCESS class.
    assert pri["AChE_inhibitor"].verdict == "SUCCESS"
    assert pri["AChE_inhibitor"].mean > 0.30
    assert pri["AChE_inhibitor"].k_total >= 30
    # the famous AD failures are FAILURE-verdict classes
    assert pri["AMPA_PAM"].verdict == "FAILURE"
    assert pri["PDE9_PDE10"].verdict == "FAILURE"
    assert pri["5HT6_antagonist"].verdict == "FAILURE"


@pytest.mark.skipif(not (LEDGER.exists() and ANCHORS.exists()),
                    reason="ledger / anchors not present")
def test_real_cias_muscarinic_success():
    ledger = R.load_clinical_ledger(LEDGER)
    anchors = pd.read_csv(ANCHORS, comment="#")
    ev = D.load_disease_evidence(ledger, anchors)
    pri = D.build_disease_class_priors("CIAS", ev)
    # xanomeline-KarXT → muscarinic M1/M4 is the recent CIAS success
    assert "M1_M4_agonist" in pri
    assert pri["M1_M4_agonist"].verdict == "SUCCESS"
    # the α7 graveyard (encenicline et al.) is a FAILURE class
    assert pri["alpha7_nAChR"].verdict == "FAILURE"


@pytest.mark.skipif(not (LEDGER.exists() and ANCHORS.exists()),
                    reason="ledger / anchors not present")
def test_real_fxs_pde4_success():
    ledger = R.load_clinical_ledger(LEDGER)
    anchors = pd.read_csv(ANCHORS, comment="#")
    ev = D.load_disease_evidence(ledger, anchors)
    pri = D.build_disease_class_priors("FXS", ev)
    assert pri["PDE4_inhibitor"].verdict == "SUCCESS"     # zatolmilast g≈0.71
    assert pri["PDE4_inhibitor"].mean > 0.5
    assert pri["mGluR"].verdict == "FAILURE"              # basimglurant/mavoglurant


@pytest.mark.skipif(not LEDGER.exists(), reason="ledger not present")
def test_real_within_ad_class_beats_relevance():
    ledger = R.load_clinical_ledger(LEDGER)
    v6b = pd.read_parquet(V6B) if V6B.exists() else None
    wd = D.within_disease_class_loco("AD", ledger, v6b_theta=v6b)
    assert wd.n_success >= 3 and wd.n_fail >= 5
    assert wd.auroc_class >= 0.85           # class track record strongly discriminates
    assert wd.failure_recall >= 0.9         # flags the AD failures
    if np.isfinite(wd.auroc_relevance):
        assert wd.auroc_class >= wd.auroc_relevance   # class ≥ target relevance


@pytest.mark.skipif(not LEDGER.exists(), reason="ledger not present")
def test_real_ad_anchor_g_donepezil():
    ledger = R.load_clinical_ledger(LEDGER)
    a = D.disease_anchor_g("AD", ledger)
    assert "donepezil" in a
    g, g90 = a["donepezil"]
    assert 0.30 <= g <= 0.42                 # real AD ADAS-Cog g
    assert g90 > g
