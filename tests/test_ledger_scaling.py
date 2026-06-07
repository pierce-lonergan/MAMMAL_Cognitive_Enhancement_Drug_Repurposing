"""Tests for the F3 ledger-scaling / per-domain / power harness.

numpy/scipy only. Unit tests for the pure functions plus a real-ledger
regression that locks the F3 trajectory headline (class separation survives the
n=31 -> 47 expansion).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mammal_repurposing.validation.ledger_scaling import (
    load_all_ledgers, assign_domain, scaling_trajectory,
    per_domain_separation, within_class_power_roadmap, n_eff_for_rho,
)

ROOT = Path(__file__).resolve().parents[1]
PATHS = [
    ROOT / "data" / "raw" / "clinical_outcomes_ledger.csv",
    ROOT / "data" / "raw" / "clinical_outcomes_ledger_EXTENSION.csv",
    ROOT / "data" / "raw" / "clinical_outcomes_ledger_CTGOV.csv",
]
HAVE = all(p.exists() for p in PATHS)


# --------------------------------------------------------------------------
# pure functions
# --------------------------------------------------------------------------

def test_assign_domain_known_endpoints():
    assert assign_domain("ADAS-Cog") == "global_amnestic"
    assert assign_domain("SIB/ADAS-Cog") == "global_amnestic"
    assert assign_domain("DSST") == "processing_speed"
    assert assign_domain("MCCB composite") == "scz_composite_battery"
    assert assign_domain("RAVLT") == "episodic_memory"
    assert assign_domain("ADHD-RS") == "adhd_symptom"
    assert assign_domain("ESS/MWT") == "wakefulness"
    assert assign_domain("something weird") == "other"


def test_n_eff_for_rho_monotone_and_known():
    # smaller target effect -> needs more points
    assert n_eff_for_rho(0.3) > n_eff_for_rho(0.4) > n_eff_for_rho(0.5)
    # known Fisher-z values (two-sided alpha=0.05, power=0.80)
    assert n_eff_for_rho(0.5) == pytest.approx(30, abs=2)
    assert n_eff_for_rho(0.3) == pytest.approx(85, abs=2)
    assert n_eff_for_rho(0.0) == -1     # degenerate guarded


def test_load_all_ledgers_synthetic(tmp_path):
    a = tmp_path / "a.csv"
    b = tmp_path / "b.csv"
    cols = "compound,mechanism_class,target_uniprot,indication,pivotal_trial,readout_year,clinical_outcome,clinical_g,endpoint,citation"
    a.write_text(cols + "\n"
                 "d1,X,P1,AD,t,2010,SUCCESS,0.4,ADAS-Cog,c\n"
                 "d2,X,P1,AD,t,2011,FAILURE,0.0,ADAS-Cog,c\n", encoding="utf-8")
    b.write_text(cols + "\n"
                 "d2,X,P1,AD,t,2011,FAILURE,0.0,ADAS-Cog,c\n"   # dup -> dropped
                 "d3,Y,P2,SCZ,t,2012,FAILURE,0.0,MCCB,c\n"
                 "x1,Z,P3,AD,t,2013,UNADJUDICATED,0.0,ADAS,c\n", encoding="utf-8")
    led = load_all_ledgers([a, b])
    assert set(led["compound"]) == {"d1", "d2", "d3"}     # dup + non-binary dropped
    assert led["label"].tolist() == [1, 0, 0]


# --------------------------------------------------------------------------
# real-ledger regression: the F3 headline
# --------------------------------------------------------------------------

@pytest.mark.skipif(not HAVE, reason="ledger CSVs not present")
def test_scaling_preserves_separation():
    traj = scaling_trajectory(PATHS, ["base", "ext", "ctgov"], n_perm=500, seed=0)
    assert len(traj) == 3
    assert traj[0].n == 31 and traj[-1].n > 40        # ledger grows
    assert traj[-1].n_classes > traj[0].n_classes      # new classes added
    # class separation SURVIVES scaling
    assert traj[-1].auroc > 0.90
    assert traj[-1].frac_pure == pytest.approx(1.0)    # all classes outcome-pure
    assert traj[-1].frac_between > 0.90 and traj[-1].icc1 > 0.85
    for s in traj:
        assert s.perm_p < 0.01                          # always significant


@pytest.mark.skipif(not HAVE, reason="ledger CSVs not present")
def test_per_domain_and_power_roadmap():
    full = load_all_ledgers(PATHS)
    dom = per_domain_separation(full)
    assert "global_amnestic" in dom
    assert dom["global_amnestic"]["n"] >= 10
    # the largest, mixed-outcome domain shows class separation
    assert np.isfinite(dom["global_amnestic"]["auroc"])
    assert dom["global_amnestic"]["auroc"] > 0.7

    pr = within_class_power_roadmap(full)
    assert pr.cur_pooled_points > 0
    assert pr.cur_within_classes >= 1
    # needs MORE than the current ledger to reach rho=0.3 at 80% power
    assert pr.targets[0.3]["implied_total_n"] > full.shape[0]
    # easier targets need fewer points than harder ones
    assert pr.targets[0.5]["n_eff"] < pr.targets[0.3]["n_eff"]
