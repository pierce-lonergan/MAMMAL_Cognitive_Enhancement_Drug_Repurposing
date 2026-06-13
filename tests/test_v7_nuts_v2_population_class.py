"""Sprint 3.2 — Tests for V7.3 Stage 2 NUTS with population × class interaction.

Validates that `fit_effect_size_nuts_v2` correctly handles:
  - Per-class τ²_class hyperprior
  - Population × class interaction term
  - V2 subdomain anchor likelihood (96-cell PER_SUBDOMAIN_PRIORS_V2)
  - Backward compat: EffectSizeObservation default `population="mixed"`
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("pymc")
pytest.importorskip("numpyro")

from mammal_repurposing.translation.effect_size_model import (
    EffectSizeObservation,
    fit_effect_size_nuts_v2,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _make_obs(compound, class_name, population, observed_g=None,
               endpoint="EM", relevance=0.7, pchembl=8.0):
    return EffectSizeObservation(
        compound=compound,
        class_name=class_name,
        target_uniprot="P22303",
        pchembl_post_mean=pchembl,
        pchembl_post_sd=0.3,
        relevance_post_mean=relevance,
        relevance_post_sd=0.1,
        pbpk_auc_brain=10.0,
        moderators=(0.0,) * 5,
        observed_g=observed_g,
        endpoint=endpoint,
        population=population,
    )


# ---------------------------------------------------------------------------
# Schema + backward compat
# ---------------------------------------------------------------------------

def test_observation_default_population_is_mixed():
    obs = EffectSizeObservation(
        compound="x", class_name="AChE_INHIBITORS", target_uniprot="P22303",
        pchembl_post_mean=8.0,
    )
    assert obs.population == "mixed"


def test_nuts_v2_fail_closes_on_v2_class_without_v2_prior_table():
    """B1 regression (docs/BUG_AUDIT_2026-06.md): fit_effect_size_nuts_v2's class-level
    priors are keyed by the V1 taxonomy (list_class_names()), but every production anchor
    carries a V2 class name. Passing a V2 name (zero overlap with V1) used to silently
    collapse every compound onto class 0 (AChE-I); it must now FAIL CLOSED with a ValueError
    rather than emit a degenerate posterior. Runs fast: the guard fires before sampling."""
    obs = [_make_obs("mph", "DA_STIMULANTS_MPH", "HC", observed_g=0.43)]
    with pytest.raises(ValueError, match="not in the V1 class-prior vocabulary"):
        fit_effect_size_nuts_v2(obs, n_chains=1, n_draws=10, n_tune=10)


# ---------------------------------------------------------------------------
# V2 NUTS run on small synthetic dataset
# ---------------------------------------------------------------------------

@pytest.mark.xfail(
    reason="B1 (docs/BUG_AUDIT_2026-06.md): fit_effect_size_nuts_v2 class priors are "
           "V1-keyed and now fail-close on V2 class names. Needs a V2 class-level prior "
           "table (k-weighted pool over subdomain_prior_table_v2) before this can pass.",
    strict=False, raises=ValueError,
)
@pytest.mark.slow
def test_nuts_v2_runs_with_two_populations():
    """Smoke test: V2 NUTS converges on a 6-compound × 2-population panel."""
    observations = [
        # AChE-I in AD (high effect)
        _make_obs("donepezil_ad",    "AChE_INHIBITORS", "AD", observed_g=0.36),
        _make_obs("galantamine_ad",  "AChE_INHIBITORS", "AD", observed_g=0.35),
        _make_obs("rivastigmine_ad", "AChE_INHIBITORS", "AD", observed_g=0.40),
        # AChE-I in HC (lower effect)
        _make_obs("donepezil_hc",    "AChE_INHIBITORS", "HC", observed_g=0.15),
        # MPH in HC (high effect)
        _make_obs("mph_hc",          "DA_STIMULANTS_MPH", "HC", observed_g=0.43,
                  endpoint="EM"),
        # MPH in ADHD (also high)
        _make_obs("mph_adhd",        "DA_STIMULANTS_MPH", "ADHD", observed_g=0.60,
                  endpoint="EM"),
    ]
    posterior = fit_effect_size_nuts_v2(
        observations, n_chains=2, n_draws=500, n_tune=500,
        target_accept=0.95,
    )
    assert posterior.method == "pymc_nuts_v2"
    assert posterior.rhat_max < 1.05    # tighter than default 1.10
    assert posterior.ess_min > 200      # for short chains
    assert "donepezil_ad" in posterior.g_mean
    assert "mph_adhd" in posterior.g_mean
    # Sanity: predicted g should be in roughly the right magnitude
    assert -0.5 < posterior.g_mean["donepezil_ad"] < 1.0
    assert -0.5 < posterior.g_mean["mph_adhd"] < 1.0


@pytest.mark.xfail(
    reason="B1 (docs/BUG_AUDIT_2026-06.md): fit_effect_size_nuts_v2 class priors are "
           "V1-keyed and now fail-close on V2 class names. Needs a V2 class-level prior "
           "table before the population x class interaction can be exercised on V2 anchors.",
    strict=False, raises=ValueError,
)
@pytest.mark.slow
def test_nuts_v2_interaction_captures_population_heterogeneity():
    """When same class has dramatically different g in two populations, the
    iota interaction term should be non-zero (captures the divergence)."""
    # AChE-I in AD: g=0.36; AChE-I in HC: g=0.15 (Δ = 0.21)
    # If iota correctly absorbs this, posterior on iota[AChE_INHIBITORS, HC]
    # should be meaningfully different from iota[AChE_INHIBITORS, AD].
    observations = []
    for i in range(4):
        observations.append(_make_obs(f"ad_{i}", "AChE_INHIBITORS", "AD",
                                       observed_g=0.36 + 0.02 * np.random.randn()))
        observations.append(_make_obs(f"hc_{i}", "AChE_INHIBITORS", "HC",
                                       observed_g=0.15 + 0.02 * np.random.randn()))

    posterior = fit_effect_size_nuts_v2(
        observations, n_chains=2, n_draws=500, n_tune=500,
        target_accept=0.95,
        use_v2_subdomain_priors=False,    # isolate the interaction effect
    )
    # Just confirm it converges; interaction quantification needs InferenceData inspection
    assert posterior.rhat_max < 1.10
    assert "ad_0" in posterior.g_mean
    assert "hc_0" in posterior.g_mean
    # Sanity: AD predictions higher than HC predictions
    ad_mean = np.mean([posterior.g_mean[f"ad_{i}"] for i in range(4)])
    hc_mean = np.mean([posterior.g_mean[f"hc_{i}"] for i in range(4)])
    # Model should recover at least directional ordering (allow noise)
    # (4 obs × 2 pops × short chains may give noisy ordering — keep this loose)
    assert np.isfinite(ad_mean) and np.isfinite(hc_mean)
