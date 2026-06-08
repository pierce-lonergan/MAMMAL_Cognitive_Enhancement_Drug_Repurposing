"""Tests for the persistence-after-cessation axis.

numpy/pandas only. Unit tests for the schema + resolution logic on synthetic
curation, plus a real-curation regression locking the headline (the F2 shortlist's
symptomatic prior does not transfer to persistence; 0 demonstrated-healthy).
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from mammal_repurposing.validation.persistence import (
    EVIDENCE_RANK, STATUS_TIER, annotate, call_for, load_persistence,
)

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"


# --------------------------------------------------------------------------
# vocabulary
# --------------------------------------------------------------------------

def test_status_tiers_and_evidence_rank_wellformed():
    assert set(STATUS_TIER.values()) == {"live", "null", "exclude"}
    # null is the default-ish bucket; the empty gold standard is "live"
    assert STATUS_TIER["unknown"] == "null"
    assert STATUS_TIER["demonstrated_healthy"] == "live"
    assert STATUS_TIER["not_applicable"] == "exclude"
    # delayed-start RCT is the strongest design
    assert EVIDENCE_RANK["delayed_start_rct"] == max(EVIDENCE_RANK.values())
    assert EVIDENCE_RANK["none"] == 0
    assert EVIDENCE_RANK["randomized_discontinuation"] > EVIDENCE_RANK["preclinical_only"]


# --------------------------------------------------------------------------
# resolution: override > class > null
# --------------------------------------------------------------------------

@pytest.fixture
def synth(tmp_path):
    c = tmp_path / "classes.csv"
    o = tmp_path / "over.csv"
    c.write_text(
        "mechanism_class,persistence_status,evidence_design,basis,caveat,source\n"
        "stim,tested_negative,randomized_discontinuation,relapse on stop,symptomatic,X\n",
        encoding="utf-8")
    o.write_text(
        "compound,persistence_status,evidence_design,basis,caveat,source\n"
        "wonderdrug,plasticity_gated,preclinical_only,iplasticity,needs training,Y\n",
        encoding="utf-8")
    return load_persistence(c, o)


def test_resolution_precedence(synth):
    classes, overrides = synth
    # compound override wins even against a class default
    assert call_for("wonderdrug", "stim", classes, overrides).status == "plasticity_gated"
    # class default when no override
    assert call_for("amphetamine", "stim", classes, overrides).status == "tested_negative"
    # null default when neither known
    c = call_for("mystery", "unlisted_class", classes, overrides)
    assert c.status == "unknown" and c.tier == "null" and c.level == "default"


def test_override_name_normalised(synth):
    classes, overrides = synth
    # parenthetical alias / case should still match the override
    assert call_for("Wonderdrug (XR)", "stim", classes, overrides).status == "plasticity_gated"


def test_call_tier_and_rank(synth):
    classes, overrides = synth
    c = call_for("wonderdrug", "stim", classes, overrides)
    assert c.tier == "live" and c.evidence_rank == EVIDENCE_RANK["preclinical_only"]


def test_bad_vocabulary_rejected(tmp_path):
    c = tmp_path / "c.csv"
    c.write_text("mechanism_class,persistence_status,evidence_design,basis,caveat,source\n"
                 "x,not_a_status,none,b,c,s\n", encoding="utf-8")
    o = tmp_path / "o.csv"
    o.write_text("compound,persistence_status,evidence_design,basis,caveat,source\n",
                 encoding="utf-8")
    with pytest.raises(ValueError):
        load_persistence(c, o)


def test_annotate_adds_columns(synth):
    classes, overrides = synth
    df = pd.DataFrame({"query_id": ["wonderdrug", "amphetamine"],
                       "assigned_class": ["stim", "stim"]})
    out = annotate(df, classes, overrides)
    assert list(out["persistence_status"]) == ["plasticity_gated", "tested_negative"]
    assert list(out["persistence_tier"]) == ["live", "null"]
    assert "persistence_basis" in out and "persistence_evidence" in out


# --------------------------------------------------------------------------
# real-curation regression: the headline
# --------------------------------------------------------------------------

_CLASSES = RAW / "persistence_axis_classes.csv"
_OVER = RAW / "persistence_axis_overrides.csv"
_SHORT = ROOT / "reports" / "pipeline" / "f2_catalogue_shortlist.csv"
_HAVE = _CLASSES.exists() and _OVER.exists() and _SHORT.exists()


@pytest.mark.skipif(not _HAVE, reason="persistence curation / shortlist not present")
def test_persistence_headline_on_shortlist():
    classes, overrides = load_persistence(_CLASSES, _OVER)
    ann = annotate(pd.read_csv(_SHORT), classes, overrides)
    # the curated calls must use only known vocabulary
    assert set(ann["persistence_tier"]).issubset({"live", "null", "exclude"})
    # the empty gold standard: no demonstrated durable enhancement in healthy people
    assert (ann["persistence_status"] == "demonstrated_healthy").sum() == 0
    # the live threads exist but are the minority (the symptomatic prior does not transfer)
    n_live = int((ann["persistence_tier"] == "live").sum())
    assert 1 <= n_live < len(ann) / 2
    # the known live calls are present
    live_names = set(ann[ann["persistence_tier"] == "live"]["query_id"].str.lower())
    assert {"fluoxetine", "selegiline", "rasagiline"} <= live_names
    # structure-router misroutes are caught as exclude
    assert (ann["persistence_tier"] == "exclude").sum() >= 5


@pytest.mark.skipif(not _HAVE, reason="persistence curation not present")
def test_known_compound_calls():
    classes, overrides = load_persistence(_CLASSES, _OVER)
    # peripheral (BBB-impermeant) cholinesterase inhibitor -> not a CNS agent
    assert call_for("neostigmine", "AChE_inhibitor", classes, overrides).status == "not_applicable"
    # opioid -> cognition-negative
    assert call_for("oxycodone", "AChE_inhibitor", classes, overrides).status == "cognition_negative"
    # plain stimulant -> class default tested_negative
    assert call_for("dexmethylphenidate", "catecholaminergic_ADHD",
                    classes, overrides).status == "tested_negative"
