"""Sprint 1.2 — MH8 substrate-mediated AHBA-masking tests.

Validates that `build_y_obs_from_sources` correctly inflates the AHBA sigma
for substrate-degrading enzymes (ACHE, MAO-A, MAO-B, COMT), per
research/4-tier/MH8 Methods Clarity Research.md §3-§4.

Also validates that the canonical SUBSTRATE_MEDIATED_UNIPROTS frozenset
matches the docstring + that fit_cluster_d_prior_nuts surfaces divergence
count in the ClusterDPosterior.
"""

from __future__ import annotations

import numpy as np
import pytest

from mammal_repurposing.cluster_d.bayesian_prior import (
    DEFAULT_SM_SIGMA_INFLATE,
    SUBSTRATE_MEDIATED_UNIPROTS,
    build_y_obs_from_sources,
)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

def test_substrate_mediated_uniprots_is_frozenset():
    """Canonical SM set must be a frozenset (immutable, hashable)."""
    assert isinstance(SUBSTRATE_MEDIATED_UNIPROTS, frozenset)


def test_substrate_mediated_uniprots_contents():
    """ACHE, MAO-A, MAO-B, COMT — the four enzymes flagged by MH8 doc § 3."""
    expected = {
        "P22303",   # ACHE
        "P21397",   # MAOA
        "P27338",   # MAOB
        "P21964",   # COMT
    }
    assert set(SUBSTRATE_MEDIATED_UNIPROTS) == expected


def test_default_inflation_factor_is_10x():
    """MH8 doc § 4 specifies 10× sigma inflation (100× variance)."""
    assert DEFAULT_SM_SIGMA_INFLATE == 10.0


# ---------------------------------------------------------------------------
# build_y_obs_from_sources — masking behaviour
# ---------------------------------------------------------------------------

class TestBuildYObsMH8Masking:

    def setup_method(self):
        """Mini panel: 3 substrate-mediated + 3 non-SM targets."""
        self.targets = [
            "P22303",   # ACHE — SM
            "P21397",   # MAOA — SM
            "P21964",   # COMT — SM
            "P36544",   # CHRNA7 — non-SM
            "P21728",   # DRD1 — non-SM
            "Q01959",   # SLC6A3 — non-SM
        ]
        self.y_ahba = {t: 0.5 for t in self.targets}    # uniform
        self.y_l2g = {t: 0.4 for t in self.targets}
        self.y_sc = {t: 2.5 for t in self.targets}

    def test_no_mask_when_substrate_mediated_uniprots_none(self):
        """Default behaviour (backward compatibility): no inflation."""
        y_obs, sigma_obs, source_names = build_y_obs_from_sources(
            self.targets,
            y_ahba=self.y_ahba, y_l2g=self.y_l2g, y_sc=self.y_sc,
            substrate_mediated_uniprots=None,
        )
        assert source_names == ["AHBA", "L2G", "SC"]
        # All AHBA sigmas should equal the base sigma (0.30)
        ahba_row = source_names.index("AHBA")
        assert np.allclose(sigma_obs[ahba_row, :], 0.30)

    def test_mh8_mask_inflates_only_sm_targets_ahba_row(self):
        """SM targets get 10× AHBA sigma; non-SM unchanged."""
        y_obs, sigma_obs, source_names = build_y_obs_from_sources(
            self.targets,
            y_ahba=self.y_ahba, y_l2g=self.y_l2g, y_sc=self.y_sc,
            substrate_mediated_uniprots=SUBSTRATE_MEDIATED_UNIPROTS,
        )
        ahba_row = source_names.index("AHBA")
        # ACHE, MAOA, COMT (indices 0, 1, 2) → 10× = 3.0
        assert np.isclose(sigma_obs[ahba_row, 0], 3.0)
        assert np.isclose(sigma_obs[ahba_row, 1], 3.0)
        assert np.isclose(sigma_obs[ahba_row, 2], 3.0)
        # CHRNA7, DRD1, SLC6A3 (indices 3, 4, 5) → unchanged 0.30
        assert np.isclose(sigma_obs[ahba_row, 3], 0.30)
        assert np.isclose(sigma_obs[ahba_row, 4], 0.30)
        assert np.isclose(sigma_obs[ahba_row, 5], 0.30)

    def test_mh8_mask_does_not_touch_l2g_or_sc_rows(self):
        """L2G and SC contributions to θ for SM targets are preserved."""
        y_obs, sigma_obs, source_names = build_y_obs_from_sources(
            self.targets,
            y_ahba=self.y_ahba, y_l2g=self.y_l2g, y_sc=self.y_sc,
            substrate_mediated_uniprots=SUBSTRATE_MEDIATED_UNIPROTS,
        )
        l2g_row = source_names.index("L2G")
        sc_row = source_names.index("SC")
        assert np.allclose(sigma_obs[l2g_row, :], 0.20)
        assert np.allclose(sigma_obs[sc_row, :], 0.35)

    def test_mh8_mask_does_not_touch_y_obs_values(self):
        """y_obs values are unchanged — only sigma_obs is inflated."""
        y_obs, sigma_obs, source_names = build_y_obs_from_sources(
            self.targets,
            y_ahba=self.y_ahba, y_l2g=self.y_l2g, y_sc=self.y_sc,
            substrate_mediated_uniprots=SUBSTRATE_MEDIATED_UNIPROTS,
        )
        ahba_row = source_names.index("AHBA")
        # All AHBA y_obs entries should still be 0.5 (input value)
        assert np.allclose(y_obs[ahba_row, :], 0.5)

    def test_mh8_custom_inflation_factor(self):
        """Inflation factor is configurable; 100× → SM sigma = 30.0."""
        y_obs, sigma_obs, source_names = build_y_obs_from_sources(
            self.targets,
            y_ahba=self.y_ahba, y_l2g=self.y_l2g, y_sc=self.y_sc,
            substrate_mediated_uniprots={"P22303"},
            substrate_mediated_sigma_inflate=100.0,
        )
        ahba_row = source_names.index("AHBA")
        assert np.isclose(sigma_obs[ahba_row, 0], 30.0)    # ACHE × 100

    def test_mh8_mask_with_partial_sm_set(self):
        """If only ACHE is in the SM set, MAO and COMT are NOT masked."""
        y_obs, sigma_obs, source_names = build_y_obs_from_sources(
            self.targets,
            y_ahba=self.y_ahba, y_l2g=self.y_l2g, y_sc=self.y_sc,
            substrate_mediated_uniprots={"P22303"},   # ACHE only
        )
        ahba_row = source_names.index("AHBA")
        assert np.isclose(sigma_obs[ahba_row, 0], 3.0)     # ACHE inflated
        assert np.isclose(sigma_obs[ahba_row, 1], 0.30)    # MAOA NOT inflated
        assert np.isclose(sigma_obs[ahba_row, 2], 0.30)    # COMT NOT inflated

    def test_mh8_mask_handles_sm_uniprot_not_in_panel(self):
        """SM UniProts that aren't in target_uniprots are silently ignored."""
        small_panel = ["P36544", "P21728"]    # CHRNA7, DRD1 — no SM
        y_obs, sigma_obs, source_names = build_y_obs_from_sources(
            small_panel,
            y_ahba={t: 0.5 for t in small_panel},
            y_l2g={t: 0.4 for t in small_panel},
            substrate_mediated_uniprots=SUBSTRATE_MEDIATED_UNIPROTS,
        )
        ahba_row = source_names.index("AHBA")
        assert np.allclose(sigma_obs[ahba_row, :], 0.30)    # nothing inflated

    def test_mh8_mask_with_l2g_only_no_ahba(self):
        """If AHBA source is absent, masking is a no-op."""
        y_obs, sigma_obs, source_names = build_y_obs_from_sources(
            self.targets,
            y_l2g=self.y_l2g,     # only L2G
            substrate_mediated_uniprots=SUBSTRATE_MEDIATED_UNIPROTS,
        )
        assert source_names == ["L2G"]
        assert np.allclose(sigma_obs[0, :], 0.20)    # all L2G sigmas unchanged


# ---------------------------------------------------------------------------
# Posterior dataclass — divergence count + audit field
# ---------------------------------------------------------------------------

class TestClusterDPosteriorAuditFields:

    def test_default_divergence_count_is_zero(self):
        from mammal_repurposing.cluster_d.bayesian_prior import ClusterDPosterior
        p = ClusterDPosterior(targets=["A"], theta_mean={"A": 0.0})
        assert p.n_divergences == 0

    def test_default_substrate_mediated_list_is_empty(self):
        from mammal_repurposing.cluster_d.bayesian_prior import ClusterDPosterior
        p = ClusterDPosterior(targets=["A"], theta_mean={"A": 0.0})
        assert p.substrate_mediated_uniprots == []

    def test_substrate_mediated_list_round_trips(self):
        from mammal_repurposing.cluster_d.bayesian_prior import ClusterDPosterior
        p = ClusterDPosterior(
            targets=["P22303", "P36544"],
            theta_mean={"P22303": 0.5, "P36544": 0.4},
            substrate_mediated_uniprots=["P22303"],
            n_divergences=0,
        )
        assert p.substrate_mediated_uniprots == ["P22303"]
        assert p.n_divergences == 0


# ---------------------------------------------------------------------------
# Backward compatibility — pre-MH8 callers must still work
# ---------------------------------------------------------------------------

def test_pre_mh8_call_signature_still_works():
    """The pre-MH8 V6.B headline calls didn't pass substrate_mediated_*;
    those calls must continue to work identically (backward compat)."""
    targets = ["P22303", "P36544"]
    y_ahba = {"P22303": 0.5, "P36544": 0.3}
    y_l2g = {"P22303": 0.4, "P36544": 0.6}
    y_obs, sigma_obs, source_names = build_y_obs_from_sources(
        targets, y_ahba=y_ahba, y_l2g=y_l2g,
    )
    ahba_row = source_names.index("AHBA")
    # Default sigma_ahba = 0.30 — no masking applied
    assert np.allclose(sigma_obs[ahba_row, :], 0.30)
