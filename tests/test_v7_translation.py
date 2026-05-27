"""V7 translation-layer pytest:
- V7.1 PBPK 9-compartment ODE: mass balance, brain-compartment uptake,
  occupancy curve, U-shape generator, PET-anchor residuals
- V7.2 PRISMA priors: 12-class coverage, robust MAP sampling, Roberts ceiling
- V7.3 Effect-size model: stub mode, Cluster D multiplicative gate, P1-P8
  prediction parser
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


# ---------------------------------------------------------------------------
# V7.1 — PBPK 9-compartment ODE
# ---------------------------------------------------------------------------
class TestPbpkOde:
    def test_availability_numpy_path_always_works(self):
        from mammal_repurposing.translation.pbpk import availability
        a = availability()
        assert a["available"] is True
        assert a["n_compartments"] == 9
        assert "cortex" in a["brain_compartments"]
        assert isinstance(a["jax_backend"], bool)

    def test_mass_balance_with_no_elimination(self):
        """Without clearance, total nmol must conserve (within RK4 numerical noise)."""
        from mammal_repurposing.translation.pbpk import (
            DrugParameters, PbpkConfig, simulate,
        )
        drug = DrugParameters(name="test", dose_mg=10.0,
                              clearance_Lph=0.0,    # no elimination
                              v_plasma_L=3.0,
                              mw_gmol=400.0)
        cfg = PbpkConfig(dt_h=0.01, t_end_h=4.0, backend="numpy")
        result = simulate(drug, cfg)
        total_nmol_start = result.state[0].sum()
        total_nmol_end = result.state[-1].sum()
        # Should match to 4 decimals (RK4 with dt=0.01)
        assert abs(total_nmol_end - total_nmol_start) / total_nmol_start < 1e-3

    def test_brain_uptake_after_oral_dose(self):
        """Cortex concentration must be non-zero by t=2h after oral dose."""
        from mammal_repurposing.translation.pbpk import (
            DrugParameters, PbpkConfig, simulate, COMPARTMENT_IDX,
        )
        drug = DrugParameters(name="test", dose_mg=10.0, ka_h=2.0,
                              clearance_Lph=10.0, mw_gmol=400.0)
        cfg = PbpkConfig(dt_h=0.05, t_end_h=6.0, backend="numpy")
        result = simulate(drug, cfg)
        t = result.t_h
        cortex_at_2h_idx = np.argmin(np.abs(t - 2.0))
        cortex_nM = result.concentration_nM[cortex_at_2h_idx,
                                             COMPARTMENT_IDX["cortex"]]
        assert cortex_nM > 1.0     # at least 1 nM after 2h

    def test_occupancy_curve_returns_valid_shape(self):
        from mammal_repurposing.translation.pbpk import (
            DrugParameters, PbpkConfig, simulate, occupancy_curve,
            OccupancyParameters,
        )
        drug = DrugParameters(name="test", dose_mg=10.0, mw_gmol=400.0)
        result = simulate(drug, PbpkConfig(t_end_h=4.0, backend="numpy"))
        occ = occupancy_curve(result, "cortex",
                              OccupancyParameters(Kd_nM=8.0, R_reserve=2.0))
        assert set(occ.keys()) == {"t_h", "C_nM", "O_obs", "O_eff", "R_avail"}
        assert len(occ["O_obs"]) == len(occ["t_h"])
        # Occupancy bounded [0, 1]
        assert np.all(occ["O_obs"] >= 0) and np.all(occ["O_obs"] <= 1)
        assert np.all(occ["O_eff"] >= 0) and np.all(occ["O_eff"] <= 1)

    def test_u_shape_generator_asymmetry(self):
        """U-shape: low O_post + zero O_auto → positive; high O_auto → negative."""
        from mammal_repurposing.translation.pbpk import u_shape_occupancy
        # Low postsynaptic + zero auto
        low = u_shape_occupancy(np.array([0.3]), np.array([0.0]),
                                 alpha_post=1.0, alpha_auto=1.5)
        # High auto + same post
        high = u_shape_occupancy(np.array([0.3]), np.array([0.8]),
                                  alpha_post=1.0, alpha_auto=1.5)
        assert low[0] > high[0]   # autoreceptor activation reduces signal
        assert high[0] < 0        # net negative at high autoreceptor

    def test_pet_anchor_residuals_have_3_drugs(self):
        from mammal_repurposing.translation.pbpk import (
            compute_pet_anchor_residuals, PET_ANCHORS,
        )
        assert len(PET_ANCHORS) == 3
        residuals = compute_pet_anchor_residuals()
        assert len(residuals) == 3
        for r in residuals:
            assert "drug_name" in r
            assert "expected" in r
            assert "predicted" in r
            assert "residual" in r
            # Residual must be a finite float
            assert np.isfinite(r["residual"])

    def test_higher_bbb_perm_increases_brain_concentration(self):
        """BBB permeability multiplier: 1.0 must give higher cortex AUC than 0.1."""
        from mammal_repurposing.translation.pbpk import (
            DrugParameters, PbpkConfig, simulate, COMPARTMENT_IDX,
        )
        cfg = PbpkConfig(t_end_h=8.0, backend="numpy")
        drug_high = DrugParameters("a", 10.0, bbb_permeability=1.0, mw_gmol=400.0)
        drug_low = DrugParameters("b", 10.0, bbb_permeability=0.1, mw_gmol=400.0)
        r_high = simulate(drug_high, cfg)
        r_low = simulate(drug_low, cfg)
        cortex_idx = COMPARTMENT_IDX["cortex"]
        auc_high = np.trapezoid(r_high.concentration_nM[:, cortex_idx], r_high.t_h)
        auc_low = np.trapezoid(r_low.concentration_nM[:, cortex_idx], r_low.t_h)
        assert auc_high > auc_low


# ---------------------------------------------------------------------------
# V7.2 — PRISMA-anchored class priors
# ---------------------------------------------------------------------------
class TestPrismaPriors:
    def test_availability_reports_12_classes(self):
        from mammal_repurposing.translation.prisma_priors import availability
        a = availability()
        assert a["available"] is True
        assert a["n_classes"] == 12
        assert "AChE-I" in a["class_names"]
        assert "NDRI" in a["class_names"]
        assert a["roberts_2020_ceiling"] == 0.50

    def test_all_12_classes_have_required_fields(self):
        from mammal_repurposing.translation.prisma_priors import PRISMA_CLASS_PRIORS
        assert len(PRISMA_CLASS_PRIORS) == 12
        for cp in PRISMA_CLASS_PRIORS:
            assert cp.prior_mean >= 0.0      # all positive or zero
            assert cp.prior_sd > 0.0
            assert cp.n_trials >= 1
            assert cp.n_subjects >= 1
            assert cp.peak_subdomain_g >= cp.prior_mean    # peak ≥ overall
            # Roberts ceiling sanity: no class should pre-claim g > 0.50
            assert cp.peak_subdomain_g <= 0.50, \
                f"{cp.class_name} peak g={cp.peak_subdomain_g} violates Roberts ceiling"

    def test_get_class_prior_raises_on_unknown(self):
        from mammal_repurposing.translation.prisma_priors import get_class_prior
        with pytest.raises(KeyError, match="not in PRISMA priors"):
            get_class_prior("UNKNOWN_CLASS")

    def test_robust_map_prior_sample_returns_ndraws(self):
        from mammal_repurposing.translation.prisma_priors import robust_map_prior_sample
        draws = robust_map_prior_sample("AChE-I", n_draws=500, w_mix=0.2, rng_seed=0)
        assert len(draws) == 500
        # Mean should be close to AChE-I prior (0.18) with small w_mix
        assert 0.05 < float(np.mean(draws)) < 0.30

    def test_assert_roberts_ceiling_flags_violations(self):
        from mammal_repurposing.translation.prisma_priors import assert_roberts_ceiling
        verdicts = assert_roberts_ceiling(
            {"safe": 0.25, "edge": 0.49, "over": 0.51, "way_over": 0.75},
            ceiling=0.50,
        )
        assert verdicts["safe"] == "PASS"
        assert verdicts["edge"] == "PASS"
        assert verdicts["over"] == "VIOLATION"
        assert verdicts["way_over"] == "VIOLATION"

    def test_class_prior_table_is_compact_dict(self):
        from mammal_repurposing.translation.prisma_priors import (
            class_prior_table, list_class_names,
        )
        t = class_prior_table()
        assert set(t.keys()) == set(list_class_names())
        for cls, entry in t.items():
            assert "mean" in entry and "sd" in entry
            assert "representative_drug" in entry


# ---------------------------------------------------------------------------
# V7.3 — Effect-size hierarchical Bayes
# ---------------------------------------------------------------------------
class TestEffectSizeModel:
    def test_availability_reports_5_moderators(self):
        from mammal_repurposing.translation.effect_size_model import availability
        a = availability()
        assert a["available"] is True
        assert a["n_moderators"] == 5
        assert a["stub_mode_works_without_pymc"] is True

    def test_fit_effect_size_stub_returns_posterior_per_compound(self):
        from mammal_repurposing.translation.effect_size_model import (
            EffectSizeObservation, fit_effect_size_stub,
        )
        obs = [
            EffectSizeObservation(
                compound="donepezil", class_name="AChE-I",
                target_uniprot="P22303",
                pchembl_post_mean=8.5, pchembl_post_sd=0.3,
                relevance_post_mean=0.90, relevance_post_sd=0.05,
                pbpk_auc_brain=12.0, moderators=(0, 0, 0, 0, 0),
            ),
            EffectSizeObservation(
                compound="methylphenidate", class_name="NDRI",
                target_uniprot="Q01959",
                pchembl_post_mean=7.0, pchembl_post_sd=0.25,
                relevance_post_mean=0.85, relevance_post_sd=0.08,
                pbpk_auc_brain=8.0, moderators=(0, 0, 0, 0, 0),
            ),
        ]
        post = fit_effect_size_stub(obs)
        assert post.method == "prisma_stub"
        assert set(post.g_mean.keys()) == {"donepezil", "methylphenidate"}
        # All Cluster D gates should be active for these high-relevance compounds
        assert post.cluster_d_gate_active["donepezil"] is True
        assert post.cluster_d_gate_active["methylphenidate"] is True

    def test_stub_respects_cluster_d_multiplicative_gate(self):
        """Low Cluster D relevance must shrink predicted g toward floor."""
        from mammal_repurposing.translation.effect_size_model import (
            EffectSizeObservation, fit_effect_size_stub,
        )
        obs_high = EffectSizeObservation(
            compound="high_relevance", class_name="AChE-I",
            target_uniprot="P22303",
            pchembl_post_mean=8.5, relevance_post_mean=0.90,
        )
        obs_low = EffectSizeObservation(
            compound="low_relevance", class_name="AChE-I",
            target_uniprot="P22303",
            pchembl_post_mean=8.5, relevance_post_mean=0.10,
        )
        post = fit_effect_size_stub([obs_high, obs_low])
        # high_relevance should predict larger |g| than low_relevance (gate scales)
        assert post.g_mean["high_relevance"] >= post.g_mean["low_relevance"]

    def test_stub_respects_moderator_penalties(self):
        """Active moderators subtract from g."""
        from mammal_repurposing.translation.effect_size_model import (
            EffectSizeObservation, fit_effect_size_stub,
        )
        obs_clean = EffectSizeObservation(
            compound="clean", class_name="NDRI",
            target_uniprot="Q01959",
            pchembl_post_mean=8.0, relevance_post_mean=0.85,
            moderators=(0, 0, 0, 0, 0),
        )
        obs_tolerant = EffectSizeObservation(
            compound="tolerant", class_name="NDRI",
            target_uniprot="Q01959",
            pchembl_post_mean=8.0, relevance_post_mean=0.85,
            moderators=(0, 0, 1, 0, 0),   # m3 = tolerance onset active
        )
        post = fit_effect_size_stub([obs_clean, obs_tolerant])
        assert post.g_mean["clean"] > post.g_mean["tolerant"]

    def test_assert_p1_through_p8_returns_8_verdicts(self):
        from mammal_repurposing.translation.effect_size_model import (
            EffectSizeObservation, fit_effect_size_stub, assert_p1_through_p8,
        )
        # Synthetic posterior with donepezil + encenicline at expected g
        obs = [
            EffectSizeObservation("donepezil", "AChE-I", "P22303",
                                  pchembl_post_mean=8.5, relevance_post_mean=0.90),
            EffectSizeObservation("encenicline", "AChE-I", "P36544",
                                  pchembl_post_mean=7.8, relevance_post_mean=0.55),
        ]
        post = fit_effect_size_stub(obs)
        verdicts = assert_p1_through_p8(post)
        assert len(verdicts) == 8
        for v in verdicts.values():
            assert v in ("PASS", "FAIL", "NO_COMPOUND")

    def test_fit_nuts_raises_without_pymc(self, monkeypatch):
        """If PYMC_AVAILABLE is False, fit_effect_size_nuts must raise ImportError."""
        from mammal_repurposing.translation import effect_size_model as esm
        monkeypatch.setattr(esm, "PYMC_AVAILABLE", False)
        with pytest.raises(ImportError, match="PyMC not installed"):
            esm.fit_effect_size_nuts(
                [esm.EffectSizeObservation("x", "AChE-I", "P22303",
                                            pchembl_post_mean=8.0)]
            )
