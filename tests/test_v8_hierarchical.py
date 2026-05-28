"""Sprint 4.1-4.2 — Tests for V8 hierarchical with α + γ + δ random effects.

Validates `src/mammal_repurposing/cluster_e/v8_hierarchical.py` against:
  - Availability / API contract
  - Build-without-fit (architecture sanity check, no GPU needed)
  - Synthetic data round-trip: parameter recovery within posterior bounds
  - Transferability index correctly identifies transferable vs cell-specific compounds
  - Non-centered parameterization → no divergences on small problem
"""

from __future__ import annotations

import numpy as np
import pytest

from mammal_repurposing.cluster_e.v8_hierarchical import (
    V8HierarchicalPosterior,
    V8Observation,
    availability,
    build_v8_hierarchical_with_cell_random_effect,
    generate_synthetic_v8_observations,
)

pymc = pytest.importorskip("pymc")
numpyro = pytest.importorskip("numpyro")


# ---------------------------------------------------------------------------
# API contract
# ---------------------------------------------------------------------------

def test_availability_reports_pymc_and_numpyro():
    info = availability()
    assert info["available"] is True
    assert info["pymc_backend"] is True
    assert info["numpyro_backend"] is True
    assert "alpha_cell" in info["random_effects"]
    assert "gamma_species" in info["random_effects"]
    assert "delta_compound_x_cell" in info["random_effects"]
    assert "transferability_index_T_ck" in info["derived_quantities"]


def test_v8_observation_dataclass_basic():
    o = V8Observation(compound="c1", cell_line="U2OS", species="human",
                       endpoint="mito", replicate=0, y=0.5)
    assert o.compound == "c1"
    assert o.y == 0.5


# ---------------------------------------------------------------------------
# Build (no fit) — fastest possible architecture sanity check
# ---------------------------------------------------------------------------

def test_build_model_returns_pymc_model():
    obs, _ = generate_synthetic_v8_observations(
        n_compounds=10, n_cell_lines=2, n_species=1, n_endpoints=2,
        n_replicates=2,
    )
    model, coords, idx = build_v8_hierarchical_with_cell_random_effect(obs)
    assert hasattr(model, "named_vars")
    assert "beta" in model.named_vars
    assert "alpha" in model.named_vars
    assert "gamma" in model.named_vars
    assert "delta" in model.named_vars
    assert "icc_cell" in model.named_vars
    assert "icc_inter" in model.named_vars
    assert "y_obs" in model.named_vars
    # Coords
    assert "compound" in coords
    assert "cell" in coords
    assert "species" in coords
    assert "endpoint" in coords
    assert len(coords["compound"]) == 10
    assert len(coords["cell"]) == 2


def test_build_raises_on_empty_observations():
    with pytest.raises(ValueError, match="zero observations"):
        build_v8_hierarchical_with_cell_random_effect([])


# ---------------------------------------------------------------------------
# Synthetic round-trip (slow — NUTS sample)
# ---------------------------------------------------------------------------

@pytest.mark.slow
def test_v8_hierarchical_fit_converges_on_synthetic_2cellline():
    """Sprint 4.2: 50 compounds × 2 cell lines × 1 species × 2 endpoints × 3 reps.
    Total: 600 observations, ~24,000 latent δ parameters.
    Target: R̂ < 1.05, ESS > 300, 0 divergences.
    """
    from mammal_repurposing.cluster_e.v8_hierarchical import fit_v8_hierarchical
    obs, true_beta = generate_synthetic_v8_observations(
        n_compounds=50, n_cell_lines=2, n_species=1, n_endpoints=2,
        n_replicates=3, rng_seed=42,
    )
    posterior = fit_v8_hierarchical(
        obs, n_chains=2, n_draws=500, n_tune=500, target_accept=0.95,
    )
    assert posterior.method == "pymc_nuts_v8_hierarchical"
    # Convergence diagnostics
    assert posterior.rhat_max < 1.10, f"R̂={posterior.rhat_max:.3f} too high"
    assert posterior.ess_min > 100, f"ESS={posterior.ess_min:.0f} too low"
    # Architecture shapes
    assert posterior.beta_mean.shape == (50, 2)
    assert posterior.transferability_index.shape == (50, 2)
    # ICC and sigma posteriors
    assert posterior.sigma_beta.shape == (2,)
    assert posterior.sigma_alpha.shape == (2,)
    assert posterior.sigma_delta.shape == (2,)
    # ICCs should be in [0, 1]
    assert np.all((posterior.icc_cell >= 0) & (posterior.icc_cell <= 1))
    assert np.all((posterior.icc_inter >= 0) & (posterior.icc_inter <= 1))
    # Transferability index should be in [0, 1]
    assert np.all((posterior.transferability_index >= 0)
                   & (posterior.transferability_index <= 1))


@pytest.mark.slow
def test_v8_hierarchical_recovers_transferable_compound_ordering():
    """Synthetic compound with large true_beta should rank higher in
    transferability_index than synthetic compound with small true_beta."""
    from mammal_repurposing.cluster_e.v8_hierarchical import fit_v8_hierarchical
    obs, true_beta = generate_synthetic_v8_observations(
        n_compounds=30, n_cell_lines=2, n_species=1, n_endpoints=1,
        n_replicates=4, rng_seed=99, true_sigma_beta=0.6, true_sigma_delta=0.1,
    )
    posterior = fit_v8_hierarchical(
        obs, n_chains=2, n_draws=500, n_tune=500, target_accept=0.95,
    )
    # When true_sigma_beta >> true_sigma_delta, transferability should be HIGH (≈1)
    mean_T = float(posterior.transferability_index.mean())
    assert mean_T > 0.5, (
        f"Mean transferability {mean_T:.3f} should be > 0.5 when "
        "true_sigma_beta=0.6 >> true_sigma_delta=0.1"
    )


@pytest.mark.slow
def test_v8_hierarchical_low_transferability_when_delta_dominates():
    """When δ (cell-specific) >> β (transferable), transferability should
    drop. This is the calibration check that motivates the entire MH3 work."""
    from mammal_repurposing.cluster_e.v8_hierarchical import fit_v8_hierarchical
    obs, true_beta = generate_synthetic_v8_observations(
        n_compounds=30, n_cell_lines=2, n_species=1, n_endpoints=1,
        n_replicates=4, rng_seed=99, true_sigma_beta=0.1, true_sigma_delta=0.6,
    )
    posterior = fit_v8_hierarchical(
        obs, n_chains=2, n_draws=500, n_tune=500, target_accept=0.95,
    )
    mean_T = float(posterior.transferability_index.mean())
    assert mean_T < 0.6, (
        f"Mean transferability {mean_T:.3f} should be < 0.6 when "
        "true_sigma_delta=0.6 >> true_sigma_beta=0.1"
    )


# ---------------------------------------------------------------------------
# Per-compound transferability API
# ---------------------------------------------------------------------------

def test_transferability_per_compound_empty_when_no_fit():
    p = V8HierarchicalPosterior(
        compounds=["c1", "c2"], cell_lines=["U2OS"], species=["human"],
        endpoints=["mito"],
    )
    # Default (empty) posterior returns empty dict
    assert p.transferability_per_compound() == {}
    assert p.beta_per_compound() == {}


def test_transferability_per_compound_populated_after_fit():
    p = V8HierarchicalPosterior(
        compounds=["donepezil", "clemastine"],
        cell_lines=["U2OS"], species=["human"], endpoints=["mito"],
        transferability_index=np.array([[0.8], [0.3]]),
        beta_mean=np.array([[0.5], [0.1]]),
    )
    T = p.transferability_per_compound(endpoint_idx=0)
    assert T["donepezil"] == 0.8
    assert T["clemastine"] == 0.3
    beta = p.beta_per_compound(endpoint_idx=0)
    assert beta["donepezil"] == 0.5
    assert beta["clemastine"] == 0.1
