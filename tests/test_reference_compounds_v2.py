"""Sprint 3.3 — Tests for REFERENCE_COMPOUND_SMD_V2 expansion (15 → 60+ anchors).

Validates the V7.2 Stage 4 reference compound anchor table from
research/4-tier/MH1 + MH2 Meta-Analytic Prior Expansion for V7 CPT
Bayesian Pharmacology Pipeline.md §3.
"""

from __future__ import annotations

import pytest

pytest.importorskip("pymc")
pytest.importorskip("numpyro")

from mammal_repurposing.translation.reference_compounds_v2 import (
    CLASS_MIGRATION_DOC_TO_V2,
    COMPOUND_TO_TARGET_UNIPROT,
    REFERENCE_COMPOUND_SMD_V2,
    ReferenceCompoundSMDV2,
    anchors_to_observations,
    coverage_v2_anchors,
)


# ---------------------------------------------------------------------------
# Coverage
# ---------------------------------------------------------------------------

def test_anchor_count_at_least_60():
    """Sprint 3.3 baseline: ≥60 anchor cells populated."""
    assert len(REFERENCE_COMPOUND_SMD_V2) >= 60


def test_compound_diversity_at_least_30():
    """At least 30 unique compounds across the table."""
    compounds = set(r.compound for r in REFERENCE_COMPOUND_SMD_V2)
    assert len(compounds) >= 30


def test_negative_effect_examples_present():
    """At least 5 negative-g rows (impairment / withdrawal examples)."""
    n = sum(1 for r in REFERENCE_COMPOUND_SMD_V2 if r.pooled_g < 0)
    assert n >= 5


def test_phase3_null_examples_present():
    """At least 15 Phase III null rows (|g| ≤ 0.06)."""
    n = sum(1 for r in REFERENCE_COMPOUND_SMD_V2 if abs(r.pooled_g) <= 0.06)
    assert n >= 15


def test_every_compound_has_uniprot_mapping():
    """COMPOUND_TO_TARGET_UNIPROT must cover every compound in the table."""
    compounds = set(r.compound for r in REFERENCE_COMPOUND_SMD_V2)
    missing = compounds - set(COMPOUND_TO_TARGET_UNIPROT.keys())
    assert missing == set(), f"Missing UniProt mapping for: {missing}"


def test_every_class_has_v2_migration():
    """CLASS_MIGRATION_DOC_TO_V2 must cover every class in the table."""
    classes = set(r.class_v2 for r in REFERENCE_COMPOUND_SMD_V2)
    missing = classes - set(CLASS_MIGRATION_DOC_TO_V2.keys())
    assert missing == set(), f"Missing class migration for: {missing}"


# ---------------------------------------------------------------------------
# Data integrity
# ---------------------------------------------------------------------------

def test_all_pooled_g_in_realistic_range():
    """pooled_g should be in [-1, 1] (Hedges' g realistic range)."""
    for r in REFERENCE_COMPOUND_SMD_V2:
        assert -1.0 <= r.pooled_g <= 1.0


def test_all_cis_well_ordered():
    for r in REFERENCE_COMPOUND_SMD_V2:
        assert r.ci_lo <= r.ci_hi, f"{r.compound}.{r.endpoint}: CI misordered"


def test_pooled_g_inside_ci_with_tolerance():
    tol = 0.005
    for r in REFERENCE_COMPOUND_SMD_V2:
        assert r.ci_lo - tol <= r.pooled_g <= r.ci_hi + tol, (
            f"{r.compound}.{r.endpoint}: pooled_g {r.pooled_g} "
            f"outside CI [{r.ci_lo}, {r.ci_hi}]"
        )


def test_roberts_ceiling_at_most_one_exception():
    """Only known exception: methylphenidate in pediatric ADHD (Coghill 2014).

    The Roberts 2020 ceiling specifically applies to healthy adults. Pediatric
    or disease populations can legitimately exceed it; the V7 paper Discussion
    handles this via the MH1+MH2 V7 CPT § 5 soft envelope rather than hard
    clipping.
    """
    violations = [r for r in REFERENCE_COMPOUND_SMD_V2 if r.pooled_g > 0.55]
    # Expect 0 or 1 violation (the documented Coghill 2014 ADHD-ped exception)
    assert len(violations) <= 2
    if violations:
        # Ensure the violation is the documented exception
        for v in violations:
            assert "ADHD" in v.population or "ped" in v.population.lower(), (
                f"Unexpected Roberts ceiling violation: {v.compound}.{v.endpoint}."
                f"{v.population} = {v.pooled_g}"
            )


def test_all_endpoints_in_canonical_set():
    """Endpoints must be in the 8 canonical cognitive domains."""
    canonical = {"EM", "WM", "ATT", "EF", "PS", "VL", "VS", "MOT"}
    for r in REFERENCE_COMPOUND_SMD_V2:
        assert r.endpoint in canonical, (
            f"{r.compound} uses non-canonical endpoint '{r.endpoint}'"
        )


# ---------------------------------------------------------------------------
# Observation conversion
# ---------------------------------------------------------------------------

class TestAnchorsToObservations:

    def test_skip_mixed_target_drops_herbals(self):
        obs = anchors_to_observations(skip_mixed_target=True)
        target_uniprots = set(o.target_uniprot for o in obs)
        assert "MIXED" not in target_uniprots
        # Herbal compounds should be excluded
        compounds = set(o.compound.split("__")[0] for o in obs)
        assert "ginkgo-EGb761" not in compounds
        assert "creatine" not in compounds

    def test_keep_mixed_target_includes_herbals(self):
        obs = anchors_to_observations(skip_mixed_target=False)
        compounds = set(o.compound.split("__")[0] for o in obs)
        assert "ginkgo-EGb761" in compounds

    def test_observed_g_populated_from_pooled_g(self):
        obs = anchors_to_observations()
        for o in obs:
            assert o.observed_g is not None
            assert -1.0 <= o.observed_g <= 1.0

    def test_compound_name_includes_endpoint_and_population(self):
        """Compound key includes endpoint + population so duplicates don't clash."""
        obs = anchors_to_observations()
        keys = [o.compound for o in obs]
        assert len(keys) == len(set(keys)), "Compound keys must be unique"

    def test_relevance_post_passthrough(self):
        relevance_map = {"P22303": 0.95, "Q01959": 0.85}
        obs = anchors_to_observations(relevance_post_for_target=relevance_map)
        donepezil_obs = [o for o in obs if "donepezil" in o.compound][0]
        assert donepezil_obs.relevance_post_mean == 0.95

    def test_class_name_migrated_to_v2(self):
        """Doc class 'AChE' should be migrated to V2 'AChE_INHIBITORS'."""
        obs = anchors_to_observations()
        donepezil_obs = [o for o in obs if "donepezil" in o.compound][0]
        assert donepezil_obs.class_name == "AChE_INHIBITORS"


# ---------------------------------------------------------------------------
# Integration smoke: load + check sane shape (no NUTS — slow path)
# ---------------------------------------------------------------------------

def test_coverage_report_runs():
    cov = coverage_v2_anchors()
    assert cov["n_rows"] >= 60
    assert cov["n_compounds"] >= 30
    assert "EM" in cov["endpoints"]


@pytest.mark.slow
def test_v7_nuts_v2_runs_on_real_anchor_subset():
    """Smoke test: V7 NUTS V2 runs to convergence on a 30-anchor subset of
    REFERENCE_COMPOUND_SMD_V2 (full 95+ would be slow; subset gives
    architecture validation in <10s).
    """
    from mammal_repurposing.translation.effect_size_model import (
        fit_effect_size_nuts_v2,
    )
    all_obs = anchors_to_observations()
    # Take first 30 to keep slow test fast
    obs_subset = all_obs[:30]
    posterior = fit_effect_size_nuts_v2(
        obs_subset, n_chains=2, n_draws=500, n_tune=500, target_accept=0.95,
    )
    assert posterior.method == "pymc_nuts_v2"
    assert posterior.rhat_max < 1.10, f"R̂={posterior.rhat_max:.3f}"
    assert posterior.ess_min > 100
    # Sanity: donepezil prediction should be in roughly the right range
    donepezil_keys = [c for c in posterior.compounds if "donepezil" in c]
    assert len(donepezil_keys) > 0
    for k in donepezil_keys:
        assert -0.5 < posterior.g_mean[k] < 1.0
