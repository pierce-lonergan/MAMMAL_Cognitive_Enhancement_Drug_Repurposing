"""Tests for the healthy-adult cognitive-enhancement axis (the project's actual stated goal).

The decisive, novel assertion: the mechanism-class prognostic prior that gave AUROC 1.00 on disease
pivotal trials COLLAPSES to chance on the healthy-adult ground truth, because the stimulant classes
are outcome-IMPURE - methylphenidate and d-amphetamine are the same class with the same overall SMD
(0.21) yet opposite outcomes. Only a coarse 'acute CNS stimulant' gate separates, and even that is
necessary-not-sufficient (d-amphetamine is a stimulant that does nothing).
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "data" / "raw" / "healthy_adult_cognition_ledger.csv"
SCRIPT = ROOT / "scripts" / "120_healthy_adult_axis.py"
pytest.importorskip("pandas")


def _load_script():
    spec = importlib.util.spec_from_file_location("haa", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _primary():
    import pandas as pd
    df = pd.read_csv(LEDGER)
    return df[df["evidence_tier"] == "clean_MA"].reset_index(drop=True)


def test_ledger_ground_truth_is_verified_and_stimulant_confined():
    prim = _primary()
    enh = set(prim[prim["enhances_healthy_young"] == 1]["compound"])
    nul = set(prim[prim["enhances_healthy_young"] == 0]["compound"])
    # the clean healthy-adult enhancers are exactly the four acute CNS stimulants
    assert enh == {"methylphenidate", "modafinil", "caffeine", "nicotine"}
    # every clean enhancer is in the stimulant supergroup; no non-stimulant clears the bar
    assert (prim[prim["enhances_healthy_young"] == 1]["supergroup"] == "stimulant").all()
    assert (prim[prim["supergroup"] != "stimulant"]["enhances_healthy_young"] == 0).all()
    # the decisive impurity: d-amphetamine and MPH share class AND overall SMD, opposite outcome
    mph = prim[prim["compound"] == "methylphenidate"].iloc[0]
    amp = prim[prim["compound"] == "dextroamphetamine"].iloc[0]
    assert mph["mechanism_class"] == amp["mechanism_class"] == "catecholaminergic"
    assert mph["representative_g"] == amp["representative_g"] == 0.21
    assert mph["enhances_healthy_young"] == 1 and amp["enhances_healthy_young"] == 0
    # the clean-MA null set is the 7 non-clearing compounds (stimulant nulls + every botanical/nutrient)
    assert len(nul) == 7 and {"dextroamphetamine", "guarana", "ginkgo_biloba"} <= nul


def test_class_prior_collapses_but_stimulant_gate_separates():
    from mammal_repurposing.validation.retrospective import auroc
    haa = _load_script()
    prim = _primary()
    y = prim["enhances_healthy_young"].to_numpy(float)
    # the DISEASE winner (mechanism-class prognostic prior) is at/below chance here
    au_class = auroc(haa.class_loco_score(prim), y)
    assert au_class <= 0.65, f"class prior should fail on healthy adults, got {au_class}"
    # the coarse stimulant gate is the only real separator and beats the class prior
    stim = (prim["supergroup"] == "stimulant").astype(float).to_numpy()
    au_stim = auroc(stim, y)
    assert au_stim >= 0.80 and au_stim > au_class
    # at least one multi-member mechanism class is outcome-IMPURE (the homogeneity break)
    pur = haa.class_purity(prim)
    impure_multi = [r for r in pur if (r[1] + r[2]) >= 2 and not r[3]]
    assert len(impure_multi) >= 1
