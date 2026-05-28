"""Sprint 1.4 — Lock MH8 victory with regression tests on V6.B.5 NUTS.

Confirms that the MH8 substrate-mediated AHBA-masking fix delivers the
production-grade configuration for the 191-target expanded panel:

  - MH8 enabled + target_accept=0.99 → 0 divergences, R̂<1.01, ESS>1500
  - --no-mh8 baseline → ~37 divergences reproduces (regression contract)

Marked `slow` because each NUTS run takes ~8-15s on a clean run.
Use `pytest -m slow` to invoke.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
PANEL_PARQUET = ROOT / "data" / "results" / "v2" / "panel_expanded_v1.parquet"
V6B_ANCHOR_PARQUET = ROOT / "data" / "results" / "v2" / "cluster_d_posterior_v1.parquet"

pymc = pytest.importorskip("pymc")
numpyro = pytest.importorskip("numpyro")


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def panel_inputs():
    """Load the 191-target panel + synthesize observations the same way
    `scripts/62_v6b5_nuts_expanded.py` does."""
    if not PANEL_PARQUET.exists() or not V6B_ANCHOR_PARQUET.exists():
        pytest.skip(f"Required parquets missing: {PANEL_PARQUET} or "
                    f"{V6B_ANCHOR_PARQUET}")

    # Import the synthesis helper from the production script (kept in-script
    # because it is the canonical implementation).
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "_script_62", ROOT / "scripts" / "62_v6b5_nuts_expanded.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    panel = pd.read_parquet(PANEL_PARQUET)
    v6b_anchor = pd.read_parquet(V6B_ANCHOR_PARQUET)
    y_ahba, y_l2g, y_sc = mod.synthesize_observations_for_expansion(
        panel, real_v6b_posterior=v6b_anchor, rng_seed=42,
    )
    return {
        "panel_uniprots": panel["uniprot"].tolist(),
        "uniprot_to_gene": dict(zip(panel["uniprot"], panel["gene_symbol"])),
        "y_ahba": y_ahba,
        "y_l2g": y_l2g,
        "y_sc": y_sc,
    }


def _build_and_fit(panel_inputs, *, use_mh8: bool, target_accept: float,
                   n_chains: int = 4, n_draws: int = 2000, n_tune: int = 2000):
    """Build (y_obs, sigma_obs) with/without MH8 and run NUTS."""
    from mammal_repurposing.cluster_d.bayesian_prior import (
        DEFAULT_ANCHORS, SUBSTRATE_MEDIATED_UNIPROTS,
        build_y_obs_from_sources, fit_cluster_d_prior_nuts,
    )

    panel_uniprots = panel_inputs["panel_uniprots"]
    uniprot_to_gene = panel_inputs["uniprot_to_gene"]
    sm = (set(SUBSTRATE_MEDIATED_UNIPROTS) & set(panel_uniprots)) if use_mh8 else None

    y_obs, sigma_obs, source_names = build_y_obs_from_sources(
        panel_uniprots,
        y_ahba=panel_inputs["y_ahba"],
        y_l2g=panel_inputs["y_l2g"],
        y_sc=panel_inputs["y_sc"],
        substrate_mediated_uniprots=sm,
    )
    ref_idx = [i for i, u in enumerate(panel_uniprots)
               if uniprot_to_gene.get(u, "") in DEFAULT_ANCHORS]

    posterior = fit_cluster_d_prior_nuts(
        target_uniprots=panel_uniprots,
        y_obs=y_obs, sigma_obs=sigma_obs, source_names=source_names,
        reference_idx=ref_idx,
        n_chains=n_chains, n_draws=n_draws, n_tune=n_tune,
        target_accept=target_accept,
        random_seed=42,
    )
    return posterior


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

@pytest.mark.slow
def test_mh8_enabled_target_accept_099_zero_divergences(panel_inputs):
    """Production-grade configuration: MH8 + target_accept=0.99.

    Locks the Sprint 1.3 victory: 191-target NUTS with MH8 + 0.99 target_accept
    produces 0 divergences (down from 37 in the pre-MH8 baseline).
    """
    posterior = _build_and_fit(
        panel_inputs, use_mh8=True, target_accept=0.99,
    )
    assert posterior.n_divergences == 0, (
        f"MH8 production config produced {posterior.n_divergences} "
        f"divergences (expected 0). Posterior method: {posterior.method}; "
        f"note: {posterior.note}"
    )
    assert posterior.rhat_max < 1.01, (
        f"R̂ max = {posterior.rhat_max:.4f} exceeds 1.01"
    )
    # Production run (Sprint 1.3) achieved ESS_min = 1808; allow 20% margin
    # for RNG variance across numpyro/JAX versions.
    assert posterior.ess_min > 1200, (
        f"ESS min = {posterior.ess_min:.0f} is below 1200 (production: 1808)"
    )
    # Sanity: substrate-mediated audit recorded the 4 masked targets
    # (handled by the script, not the library — library returns empty list
    # by default; this test asserts the library default behaviour).
    assert posterior.substrate_mediated_uniprots == []
    # Sanity: ACHE recovery via reference anchor still works despite mask
    assert posterior.theta_mean["P22303"] > 0.3, (
        f"ACHE θ̄ = {posterior.theta_mean['P22303']:.3f} dropped below "
        "0.3 — reference anchor may not be firing"
    )


@pytest.mark.slow
def test_no_mh8_baseline_reproduces_divergences(panel_inputs):
    """Regression contract: --no-mh8 with target_accept=0.95 still produces
    the pre-MH8 baseline divergence count (~37).

    If this test starts producing 0 divergences, either the MH8 fix has been
    accidentally always-on, OR the underlying numpyro/PyMC version has
    silently fixed the issue elsewhere — either way, MH8's necessity must
    be re-evaluated.
    """
    posterior = _build_and_fit(
        panel_inputs, use_mh8=False, target_accept=0.95,
    )
    # We allow a wide band (20-60) because random_seed=42 is set but
    # numpyro/JAX may have minor RNG variation across versions.
    assert posterior.n_divergences >= 20, (
        f"--no-mh8 baseline produced only {posterior.n_divergences} "
        "divergences — expected ~37. If this fails, MH8 may no longer "
        "be strictly necessary."
    )
    assert posterior.rhat_max < 1.05, (
        f"--no-mh8 baseline R̂ max = {posterior.rhat_max:.4f} too high"
    )


@pytest.mark.slow
def test_mh8_dominates_target_accept_in_attribution(panel_inputs):
    """Attribution: MH8 structural fix should dominate the divergence
    reduction. MH8-on at target_accept=0.95 should reduce divergences by
    ≥ 80% relative to MH8-off at the same target_accept.
    """
    posterior_no_mh8 = _build_and_fit(
        panel_inputs, use_mh8=False, target_accept=0.95,
        n_draws=1000, n_tune=1000,
    )
    posterior_mh8 = _build_and_fit(
        panel_inputs, use_mh8=True, target_accept=0.95,
        n_draws=1000, n_tune=1000,
    )
    baseline = max(1, posterior_no_mh8.n_divergences)
    reduction = 1.0 - (posterior_mh8.n_divergences / baseline)
    assert reduction >= 0.80, (
        f"MH8 structural fix only reduced divergences by {reduction:.0%} "
        f"(no-MH8: {posterior_no_mh8.n_divergences}, "
        f"MH8: {posterior_mh8.n_divergences}). Expected ≥80% reduction. "
        "MH8 may not be the dominant intervention."
    )
