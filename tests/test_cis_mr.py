"""Pin the two steps where a Mendelian randomisation silently inverts or inflates.

Neither failure mode announces itself. A harmonisation error produces a confident result pointing
the wrong way, and counting datasets instead of distinct outcome effects turns one test into an
apparent replication across four cohorts.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SCREEN = ROOT / "data" / "raw" / "cis_mr_screen.csv"


@pytest.fixture(scope="module")
def m():
    spec = importlib.util.spec_from_file_location("s140", ROOT / "scripts" / "140_cis_mr_screen.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# --- harmonisation: the sign is the whole result ----------------------------------------------

def test_effect_allele_matching_alt_keeps_the_sign(m):
    assert m.harmonise("T", "C", "C", "T", 0.30) == (1, "")


def test_effect_allele_matching_ref_flips_the_sign(m):
    """The GWAS beta refers to A1. If A1 is the eQTL REF it must be negated before the ratio."""
    assert m.harmonise("C", "T", "C", "T", 0.30) == (-1, "")


def test_strand_flipped_alleles_are_resolved_not_rejected(m):
    assert m.harmonise("A", "G", "C", "T", 0.30) == (1, "")


def test_palindromic_at_intermediate_frequency_is_dropped(m):
    """A/T at EAF 0.5 cannot be strand-resolved. Guessing is a coin flip on the direction."""
    sign, why = m.harmonise("A", "T", "A", "T", 0.50)
    assert sign == 0 and "palindromic" in why


def test_palindromic_at_extreme_frequency_is_usable(m):
    sign, _ = m.harmonise("A", "T", "A", "T", 0.12)
    assert sign in (1, -1)


def test_mismatched_alleles_are_rejected(m):
    sign, why = m.harmonise("G", "C", "A", "T", 0.20)
    assert sign == 0 and "do not match" in why


def test_indels_are_rejected(m):
    sign, why = m.harmonise("CAAA", "C", "CAAA", "C", 0.30)
    assert sign == 0 and "indel" in why


def test_harmonisation_is_symmetric_under_relabelling(m):
    """Swapping A1/A2 must flip the sign and nothing else."""
    a = m.harmonise("T", "C", "C", "T", 0.3)[0]
    b = m.harmonise("C", "T", "C", "T", 0.3)[0]
    assert a == -b


# --- instrument fallback ----------------------------------------------------------------------

def test_chooser_falls_back_to_a_measured_variant(m):
    """The top-PIP variant is often absent from the outcome GWAS; the next one tags the same
    signal, so falling back costs nothing and recovers the gene."""
    cand = pd.DataFrame([
        {"gene": "X", "dataset": "D1", "rsid": "rs_missing", "pip": 0.9},
        {"gene": "X", "dataset": "D1", "rsid": "rs_present", "pip": 0.4},
    ])
    got = m.choose_instruments(cand, {"rs_present"})
    assert len(got) == 1 and got.iloc[0]["rsid"] == "rs_present"


def test_chooser_prefers_the_highest_pip_among_measured(m):
    cand = pd.DataFrame([
        {"gene": "X", "dataset": "D1", "rsid": "rs_a", "pip": 0.3},
        {"gene": "X", "dataset": "D1", "rsid": "rs_b", "pip": 0.8},
    ])
    got = m.choose_instruments(cand, {"rs_a", "rs_b"})
    assert got.iloc[0]["rsid"] == "rs_b"


def test_chooser_returns_nothing_when_no_candidate_was_measured(m):
    cand = pd.DataFrame([{"gene": "X", "dataset": "D1", "rsid": "rs_a", "pip": 0.9}])
    assert m.choose_instruments(cand, {"rs_zzz"}).empty


# --- the replication illusion -----------------------------------------------------------------

def test_same_variant_across_datasets_is_not_independent_evidence():
    """All exposure datasets share one outcome GWAS. Same rsid means the same outcome beta, so a
    gene in four datasets can still carry exactly one piece of information about cognition."""
    d = pd.read_csv(SCREEN)
    ok = d[d["status"] == "ok"]
    multi = ok[ok["n_datasets"] > 1]
    assert len(multi) > 0
    for gene, g in multi.groupby("gene"):
        expected = g["gwas_beta_aligned"].round(6).nunique() > 1
        assert bool(g["independent_outcome_evidence"].iloc[0]) == expected, gene


def test_at_least_one_gene_is_flagged_as_repeating_an_identical_outcome_effect():
    """Measured on 2026-09-21: HTR1D and F7 select the SAME rsid in two datasets each."""
    d = pd.read_csv(SCREEN)
    ok = d[(d["status"] == "ok") & (d["n_datasets"] > 1)]
    fake = ok[ok["independent_outcome_evidence"] == False]["gene"].unique()  # noqa: E712
    assert len(fake) >= 1
    assert {"HTR1D", "F7"} <= set(fake)


# --- the live screen --------------------------------------------------------------------------

def test_the_screen_claims_no_significant_target():
    """First pass returned nothing at FDR < 0.05. Changing this must be deliberate."""
    d = pd.read_csv(SCREEN)
    ok = d[d["status"] == "ok"]
    assert len(ok) > 100
    assert (ok["mr_q"] < 0.05).sum() == 0, "a target now survives FDR; update this deliberately"


def test_dropped_instruments_carry_a_reason():
    d = pd.read_csv(SCREEN)
    bad = d[d["status"] != "ok"]
    assert bad["status"].str.len().gt(5).all()
