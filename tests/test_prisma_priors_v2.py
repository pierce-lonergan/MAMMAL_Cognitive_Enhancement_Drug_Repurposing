"""Sprint 3.1 — Tests for V7.2 Stage 3 PRISMA prior expansion.

Validates the 12-class × 8-cognitive-domain = 96-cell PER_SUBDOMAIN_PRIORS_V2
table against:
  - Coverage >= 70 populated cells (target: 73 per current curation)
  - Schmidli 2014 τ² rule correctly applied per k
  - Roberts 2020 ceiling: no cell has pooled_g > 0.55
  - CI ordering and pooled_g containment
  - V1 → V2 class name migration
"""

from __future__ import annotations

import pytest

from mammal_repurposing.translation.prisma_priors import (
    CLASS_NAME_MIGRATION_V1_TO_V2,
    COGNITIVE_DOMAINS_V2,
    PER_SUBDOMAIN_PRIORS_V2,
    SubdomainPriorV2,
    coverage_v2,
    get_subdomain_prior_v2,
    list_class_names_v2,
    subdomain_prior_table_v2,
)


# ---------------------------------------------------------------------------
# Schema + coverage
# ---------------------------------------------------------------------------

def test_eight_canonical_cognitive_domains():
    assert COGNITIVE_DOMAINS_V2 == ["EM", "WM", "ATT", "EF", "PS", "VL", "VS", "MOT"]


def test_twelve_canonical_classes():
    classes = list_class_names_v2()
    assert len(classes) == 12
    assert "AChE_INHIBITORS" in classes
    assert "NMDA_MODULATORS" in classes
    assert "DA_STIMULANTS_MPH" in classes


def test_coverage_meets_baseline():
    """Sprint 3.1 baseline: ≥70 cells populated (out of 96)."""
    cov = coverage_v2()
    assert cov["n_classes"] == 12
    assert cov["n_domains"] == 8
    assert cov["full_grid_size"] == 96
    assert cov["populated_cells"] >= 70, (
        f"Only {cov['populated_cells']} cells populated; expected ≥70"
    )


def test_class_name_migration_keys_valid():
    """All V1 source names should map to known V2 names."""
    v2_names = set(list_class_names_v2())
    for v1, v2 in CLASS_NAME_MIGRATION_V1_TO_V2.items():
        assert v2 in v2_names, f"V1 '{v1}' maps to unknown V2 '{v2}'"


# ---------------------------------------------------------------------------
# Schmidli τ² rule
# ---------------------------------------------------------------------------

class TestSchmidliTauRule:

    def test_k_ge_5_tight_prior(self):
        sp = SubdomainPriorV2("X", "EM", 0.2, 0.1, 0.3, 5, "test")
        assert sp.tau2 == 0.02
        sp10 = SubdomainPriorV2("X", "EM", 0.2, 0.1, 0.3, 10, "test")
        assert sp10.tau2 == 0.02

    def test_k_2_4_mid_prior(self):
        for k in [2, 3, 4]:
            sp = SubdomainPriorV2("X", "EM", 0.2, 0.1, 0.3, k, "test")
            assert sp.tau2 == 0.04, f"k={k} should give tau2=0.04, got {sp.tau2}"

    def test_k_1_loose_prior(self):
        sp = SubdomainPriorV2("X", "EM", 0.2, 0.1, 0.3, 1, "test")
        assert sp.tau2 == 0.08

    def test_k_0_fallback(self):
        sp = SubdomainPriorV2("X", "EM", 0.2, 0.1, 0.3, 0, "test")
        assert sp.tau2 == 0.09    # HalfNormal(0.3)² = 0.09

    def test_tau_is_sqrt_tau2(self):
        sp = SubdomainPriorV2("X", "EM", 0.2, 0.1, 0.3, 5, "test")
        assert abs(sp.tau - 0.02 ** 0.5) < 1e-10


# ---------------------------------------------------------------------------
# Data integrity
# ---------------------------------------------------------------------------

class TestDataIntegrity:

    def test_no_cell_exceeds_roberts_ceiling(self):
        """No published pooled_g should exceed Roberts 2020 ceiling (0.55)."""
        violations = [sp for sp in PER_SUBDOMAIN_PRIORS_V2 if sp.pooled_g > 0.55]
        assert violations == [], (
            f"{len(violations)} Roberts ceiling violations: "
            f"{[(v.class_name, v.domain, v.pooled_g) for v in violations]}"
        )

    def test_all_cis_well_ordered(self):
        for sp in PER_SUBDOMAIN_PRIORS_V2:
            assert sp.ci_lo <= sp.ci_hi, (
                f"{sp.class_name}.{sp.domain}: CI_lo {sp.ci_lo} > CI_hi {sp.ci_hi}"
            )

    def test_pooled_g_inside_ci(self):
        """pooled_g should lie inside [ci_lo, ci_hi] (with small tolerance)."""
        tol = 0.005
        for sp in PER_SUBDOMAIN_PRIORS_V2:
            assert sp.ci_lo - tol <= sp.pooled_g <= sp.ci_hi + tol, (
                f"{sp.class_name}.{sp.domain}: pooled_g {sp.pooled_g} "
                f"outside CI [{sp.ci_lo}, {sp.ci_hi}]"
            )

    def test_all_pooled_g_in_reasonable_range(self):
        """pooled_g should be in [-1, 1] (Hedges' g realistic range)."""
        for sp in PER_SUBDOMAIN_PRIORS_V2:
            assert -1.0 <= sp.pooled_g <= 1.0

    def test_all_k_at_least_one(self):
        for sp in PER_SUBDOMAIN_PRIORS_V2:
            assert sp.k >= 1, f"{sp.class_name}.{sp.domain} has k={sp.k}"

    def test_all_sources_nonempty(self):
        for sp in PER_SUBDOMAIN_PRIORS_V2:
            assert sp.source.strip() != "", (
                f"{sp.class_name}.{sp.domain} has empty source"
            )

    def test_all_classes_in_canonical_domain_set(self):
        for sp in PER_SUBDOMAIN_PRIORS_V2:
            assert sp.domain in COGNITIVE_DOMAINS_V2, (
                f"{sp.class_name} uses unknown domain '{sp.domain}'"
            )


# ---------------------------------------------------------------------------
# Lookup API
# ---------------------------------------------------------------------------

class TestLookupAPI:

    def test_get_subdomain_prior_v2_canonical_name(self):
        sp = get_subdomain_prior_v2("AChE_INHIBITORS", "EM")
        assert sp is not None
        assert sp.pooled_g == 0.36
        assert sp.k == 5

    def test_get_subdomain_prior_v2_v1_name_migration(self):
        """V1 'AChE-I' should auto-migrate to V2 'AChE_INHIBITORS'."""
        sp_v1 = get_subdomain_prior_v2("AChE-I", "EM")
        sp_v2 = get_subdomain_prior_v2("AChE_INHIBITORS", "EM")
        assert sp_v1 is not None
        assert sp_v2 is not None
        assert sp_v1.pooled_g == sp_v2.pooled_g

    def test_get_subdomain_prior_v2_missing_cell(self):
        # 5HT6 only has EM populated; ATT should return None
        sp = get_subdomain_prior_v2("5HT6_ANTAGONISTS", "ATT")
        assert sp is None

    def test_get_subdomain_prior_v2_unknown_class(self):
        sp = get_subdomain_prior_v2("NONEXISTENT_CLASS", "EM")
        assert sp is None

    def test_subdomain_prior_table_v2_structure(self):
        table = subdomain_prior_table_v2()
        assert "AChE_INHIBITORS" in table
        assert "EM" in table["AChE_INHIBITORS"]
        em_cell = table["AChE_INHIBITORS"]["EM"]
        assert "pooled_g" in em_cell
        assert "ci_lo" in em_cell
        assert "ci_hi" in em_cell
        assert "k" in em_cell
        assert "tau2" in em_cell
        assert "source" in em_cell


# ---------------------------------------------------------------------------
# Production class-level cells (key sanity checks)
# ---------------------------------------------------------------------------

class TestProductionAnchors:
    """Spot-checks on individual cells that are central to the V7 paper."""

    def test_donepezil_AD_episodic_memory(self):
        sp = get_subdomain_prior_v2("AChE_INHIBITORS", "EM")
        # Birks 2018 Cochrane: AD donepezil EM g=0.36 (CI 0.27-0.44, k=5)
        assert sp.pooled_g == 0.36
        assert sp.ci_lo == 0.27
        assert sp.ci_hi == 0.44
        assert sp.k == 5
        assert "Birks2018" in sp.source

    def test_methylphenidate_healthy_episodic_memory(self):
        sp = get_subdomain_prior_v2("DA_STIMULANTS_MPH", "EM")
        # Roberts 2020: MPH HC EM g=0.43 (CI 0.21-0.65, k=24)
        assert sp.pooled_g == 0.43
        assert sp.k == 24
        assert "Roberts2020" in sp.source

    def test_modafinil_executive_function(self):
        sp = get_subdomain_prior_v2("MODAFINIL_LIKE", "EF")
        # Roberts 2020: modafinil EF g=0.28 (CI 0.03-0.53, k=14) — the only sig subdomain
        assert sp.pooled_g == 0.28
        assert sp.ci_lo == 0.03
        assert sp.k == 14

    def test_alpha7_phase3_failures_negative_or_null(self):
        """All ALPHA7_NACHR cells should have pooled_g ≤ 0.10 (Phase III nulls)."""
        for sp in PER_SUBDOMAIN_PRIORS_V2:
            if sp.class_name == "ALPHA7_NACHR":
                assert sp.pooled_g <= 0.10, (
                    f"ALPHA7_NACHR {sp.domain} has g={sp.pooled_g}; "
                    "Phase III evidence is uniformly null"
                )

    def test_5ht6_only_em_cell(self):
        """5HT6_ANTAGONISTS has only EM data; all else None per published meta-analyses."""
        cells = [sp for sp in PER_SUBDOMAIN_PRIORS_V2
                 if sp.class_name == "5HT6_ANTAGONISTS"]
        assert len(cells) == 1
        assert cells[0].domain == "EM"
