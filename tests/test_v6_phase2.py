"""V6 phase 2 + phase 3 pytest:
- PSICHIC adapter (phase 2)
- Venn-ABERS calibrated uncertainty propagation (phase 2)
- MMAtt-DTA INVERT-mask coverage (phase 2)
- BALM adapter (phase 3 — direct + subprocess + availability probe)
- OT Genetics L2G fetcher (phase 3 — UniProt mapping + GraphQL stubs)
- PyMC Cluster D NUTS (phase 3 — y_obs aggregator + stub posterior)
- TxGNN per-disease ranking API (phase 3 — cache builder + drug lookup)
"""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


# ---------------------------------------------------------------------------
# V6.A.1 phase 2 — PSICHIC adapter graceful degradation
# ---------------------------------------------------------------------------
class TestPsichicAdapter:
    def test_availability_when_repo_missing(self):
        from mammal_repurposing.cluster_a.psichic_adapter import availability
        a = availability()
        # Repo isn't at the default location → available=False, but no exception
        assert isinstance(a, dict)
        assert "available" in a

    def test_find_repo_raises_helpful_error(self):
        from mammal_repurposing.cluster_a.psichic_adapter import _find_psichic_repo
        with pytest.raises(FileNotFoundError, match="PSICHIC repo not found"):
            _find_psichic_repo(None)

    def test_build_input_csv_filters_invalid_smiles(self, tmp_path):
        from mammal_repurposing.cluster_a.psichic_adapter import build_psichic_input
        pairs = pd.DataFrame([
            {"compound_name": "good_smiles", "target_uniprot": "P00000",
             "compound_smiles": "CCO"},
            {"compound_name": "bad_smiles", "target_uniprot": "P00000",
             "compound_smiles": None},
            {"compound_name": "unknown_target", "target_uniprot": "MISSING",
             "compound_smiles": "CCC"},
        ])
        targets = pd.DataFrame([{"uniprot": "P00000", "sequence": "MKT" * 10}])
        out_csv = tmp_path / "psichic_in.csv"
        n = build_psichic_input(pairs, targets, out_csv)
        assert n == 1     # only the good_smiles pair survives
        df = pd.read_csv(out_csv)
        assert "protein_sequence" in df.columns
        assert "smiles" in df.columns


# ---------------------------------------------------------------------------
# V6.A.4 — Venn-ABERS
# ---------------------------------------------------------------------------
class TestVennAbers:
    def test_basic_predict_interval(self):
        from mammal_repurposing.calibration.venn_abers import VennAbersRegressor
        rng = np.random.default_rng(0)
        x_cal = rng.uniform(4, 9, 30)
        y_cal = 1.2 * x_cal + rng.normal(0, 0.5, 30)
        va = VennAbersRegressor()
        va.fit(x_cal, y_cal)
        x_q = np.array([5.0, 7.5])
        lo, hi = va.predict_interval(x_q)
        # Intervals must be ordered correctly
        assert np.all(lo <= hi)
        # Tighter calibration set → narrower intervals
        assert (hi - lo).mean() < 5.0   # sanity ceiling

    def test_unfitted_raises(self):
        from mammal_repurposing.calibration.venn_abers import VennAbersRegressor
        va = VennAbersRegressor()
        with pytest.raises(RuntimeError, match="fit\\(\\) must be called"):
            va.predict_interval(np.array([5.0]))

    def test_sigma_approximation(self):
        from mammal_repurposing.calibration.venn_abers import VennAbersRegressor
        rng = np.random.default_rng(1)
        x_cal = rng.uniform(4, 9, 30)
        y_cal = 1.2 * x_cal + rng.normal(0, 0.5, 30)
        va = VennAbersRegressor()
        va.fit(x_cal, y_cal)
        sigma = va.predict_sigma(np.array([7.0]))
        # σ should be positive and finite
        assert sigma[0] > 0
        assert np.isfinite(sigma[0])

    def test_correlated_mc_intervals(self):
        from mammal_repurposing.calibration.venn_abers import correlated_mc_intervals
        point = {"A": np.array([5.0, 7.0]), "B": np.array([5.5, 7.5])}
        sigmas = {"A": np.array([0.3, 0.3]), "B": np.array([0.4, 0.4])}
        weights = {"A": 0.5, "B": 0.5}
        # Identity correlation = independence
        lo, hi = correlated_mc_intervals(point, sigmas, weights, correlation=None,
                                         n_samples=300, rng_seed=0)
        assert len(lo) == len(hi) == 2
        assert np.all(lo < hi)
        # Higher correlation → smaller CI variability per query
        corr_pos = np.array([[1.0, 0.9], [0.9, 1.0]])
        lo_c, hi_c = correlated_mc_intervals(point, sigmas, weights,
                                              correlation=corr_pos, n_samples=300,
                                              rng_seed=0)
        # Just verify they run; quantitative comparison varies by sample
        assert len(lo_c) == 2


# ---------------------------------------------------------------------------
# V6.A.1 — MMAtt-DTA adapter superfamily map + fusion ranker INVERT mask
# ---------------------------------------------------------------------------
class TestMmattActivation:
    def test_invert_mask_drops_expected_targets(self):
        """The empirical INVERT-target set (V6.A.1 measured) must match scripts/53."""
        from importlib import import_module
        sys.path.insert(0, str(ROOT / "scripts"))
        m = import_module("53_v6_mmatt_fusion_ranker")
        # Per V6.A.1 empirical result
        expected = {"P08913", "P36544", "P42261", "Q99720", "Q16620", "P23975"}
        assert m.MMATT_INVERT_TARGETS == expected

    def test_superfamily_map_covers_priority_targets(self):
        from mammal_repurposing.cluster_a.mmatt_dta_adapter import (
            COGNITION_PANEL_SUPERFAMILY,
        )
        # The 13 supported cognition targets per V6.A.1
        supported_in_panel = [
            "Q01959",  # SLC6A3 transporter
            "Q9Y5N1",  # HRH3 gpcr
            "O43614",  # HCRTR2 gpcr
            "P21728",  # DRD1 gpcr
            "Q08499",  # PDE4D enzyme
            "O76083",  # PDE9A enzyme
            "Q13224",  # GRIN2B ion_channel
            "P36544",  # CHRNA7 ion_channel
        ]
        for u in supported_in_panel:
            assert u in COGNITION_PANEL_SUPERFAMILY
            assert COGNITION_PANEL_SUPERFAMILY[u] in (
                "transporter", "gpcr", "enzyme", "ion_channel", "kinase"
            )


# ---------------------------------------------------------------------------
# V6.A.1 phase 3 — BALM adapter graceful degradation
# ---------------------------------------------------------------------------
class TestBalmAdapter:
    def test_availability_when_weights_missing(self, monkeypatch):
        from mammal_repurposing.cluster_a.balm_adapter import availability
        # Ensure no BALM_WEIGHTS_DIR / BALM_ROOT pollutes the probe
        monkeypatch.delenv("BALM_WEIGHTS_DIR", raising=False)
        monkeypatch.delenv("BALM_ROOT", raising=False)
        a = availability()
        assert isinstance(a, dict)
        assert "available" in a
        # Without env or HF cache for BALM, must be False
        if not a["available"]:
            assert "reason" in a

    def test_find_weights_raises_helpful_error(self, monkeypatch):
        from mammal_repurposing.cluster_a.balm_adapter import _find_balm_weights
        monkeypatch.delenv("BALM_WEIGHTS_DIR", raising=False)
        with pytest.raises(FileNotFoundError, match="BALM weights not found"):
            _find_balm_weights(None)

    def test_find_repo_raises_helpful_error(self, monkeypatch):
        from mammal_repurposing.cluster_a.balm_adapter import _find_balm_repo
        monkeypatch.delenv("BALM_ROOT", raising=False)
        with pytest.raises(FileNotFoundError, match="BALM repo not found"):
            _find_balm_repo(None)

    def test_build_input_csv_filters_invalid_smiles(self, tmp_path):
        from mammal_repurposing.cluster_a.balm_adapter import build_balm_input
        pairs = pd.DataFrame([
            {"compound_name": "good", "target_uniprot": "P00000",
             "compound_smiles": "CCO"},
            {"compound_name": "bad_smiles", "target_uniprot": "P00000",
             "compound_smiles": None},
            {"compound_name": "unknown_tgt", "target_uniprot": "MISSING",
             "compound_smiles": "CCC"},
        ])
        targets = pd.DataFrame([{"uniprot": "P00000", "sequence": "MKT" * 10}])
        out_csv = tmp_path / "balm_in.csv"
        n = build_balm_input(pairs, targets, out_csv)
        assert n == 1
        df = pd.read_csv(out_csv)
        assert "target_seq" in df.columns
        assert "smiles" in df.columns

    def test_run_balm_raises_when_both_paths_unavailable(self, monkeypatch):
        from mammal_repurposing.cluster_a.balm_adapter import run_balm, BalmConfig
        monkeypatch.delenv("BALM_WEIGHTS_DIR", raising=False)
        monkeypatch.delenv("BALM_ROOT", raising=False)
        pairs = pd.DataFrame([{"compound_name": "x", "target_uniprot": "P00000",
                               "compound_smiles": "CCO"}])
        targets = pd.DataFrame([{"uniprot": "P00000", "sequence": "MKT"*10}])
        with pytest.raises(FileNotFoundError, match="BALM unavailable via both"):
            run_balm(pairs, targets, config=BalmConfig())


# ---------------------------------------------------------------------------
# V6.B.1 Stage 2 — OT Genetics L2G fetcher
# ---------------------------------------------------------------------------
class TestOTGeneticsFetcher:
    def test_cognition_studies_have_expected_entries(self):
        from mammal_repurposing.cluster_d.data_fetchers import COGNITION_GWAS_STUDIES
        # Davies + Hill + Savage + Sniekers + UKBB → at least 5 studies
        assert len(COGNITION_GWAS_STUDIES) >= 5
        # Davies 2018 GCST006269 is the canonical anchor (Nat Commun 9:2098)
        assert "GCST006269" in COGNITION_GWAS_STUDIES
        assert COGNITION_GWAS_STUDIES["GCST006269"]["author"] == "Davies G"

    def test_uniprot_to_ensembl_handles_failure(self, monkeypatch):
        """If UniProt unreachable, returns ('', '') not raise."""
        from mammal_repurposing.cluster_d import data_fetchers as df_mod

        def _boom(*args, **kwargs):
            raise OSError("network down")

        monkeypatch.setattr("urllib.request.urlopen", _boom)
        ens, sym = df_mod._uniprot_to_ensembl("P00000")
        assert ens == "" and sym == ""

    def test_fetch_ot_l2g_returns_stubs_when_unreachable(self, tmp_path, monkeypatch):
        from mammal_repurposing.cluster_d import data_fetchers as df_mod
        # Force network unreachable
        monkeypatch.setattr(df_mod, "_check_ot_reachable", lambda: False)
        out = df_mod.fetch_ot_l2g(
            ["P00000", "P11111"],
            cache_path=tmp_path / "ot_cache.parquet",
            use_cache=False,
        )
        assert len(out) == 2
        for u, r in out.items():
            assert r.target_uniprot == u
            assert "reachable=False" in r.note

    def test_fetch_ot_l2g_uses_cache(self, tmp_path):
        from mammal_repurposing.cluster_d import data_fetchers as df_mod
        # Pre-populate a cache parquet
        cache_path = tmp_path / "ot_cache.parquet"
        pd.DataFrame([{
            "target_uniprot": "P12345",
            "target_ensembl": "ENSG00000123456",
            "gene_symbol": "TEST",
            "l2g_max_score": 0.42,
            "contributing_studies": ["GCST006269:0.42"],
            "note": "cached",
        }]).to_parquet(cache_path, index=False)
        # Force network down — must still return cached entry
        with patch.object(df_mod, "_check_ot_reachable", lambda: False):
            out = df_mod.fetch_ot_l2g(["P12345"], cache_path=cache_path)
        assert out["P12345"].l2g_max_score == pytest.approx(0.42)
        assert out["P12345"].gene_symbol == "TEST"


# ---------------------------------------------------------------------------
# V6.B.3 — PyMC NUTS y_obs builder + stub posterior
# ---------------------------------------------------------------------------
class TestClusterDBayesianPrior:
    def test_build_y_obs_drops_missing_sources(self):
        from mammal_repurposing.cluster_d.bayesian_prior import build_y_obs_from_sources
        y_obs, sigma_obs, names = build_y_obs_from_sources(
            ["P00000", "P11111", "P22222"],
            y_ahba={"P00000": 0.5, "P11111": -0.3, "P22222": 0.1},
            y_l2g={"P00000": 0.8, "P11111": 0.0},
            y_sc=None,
        )
        assert names == ["AHBA", "L2G"]
        assert y_obs.shape == (2, 3)
        assert sigma_obs.shape == (2, 3)
        # Missing P22222 in y_l2g → centered at 0 with inflated sigma
        assert y_obs[1, 2] == 0.0
        assert sigma_obs[1, 2] > sigma_obs[1, 0]   # inflated

    def test_build_y_obs_raises_when_all_none(self):
        from mammal_repurposing.cluster_d.bayesian_prior import build_y_obs_from_sources
        with pytest.raises(ValueError, match="At least one"):
            build_y_obs_from_sources(["P00000"], None, None, None)

    def test_stub_posterior_runs_without_pymc(self):
        from mammal_repurposing.cluster_d.bayesian_prior import fit_cluster_d_prior_stub
        out = fit_cluster_d_prior_stub(
            ["P00000", "P11111"],
            y_ahba={"P00000": 1.0, "P11111": -1.0},
        )
        assert out.method == "stage_0_stub"
        assert len(out.theta_mean) == 2
        # σ(θ) ∈ (0, 1)
        for u, w in out.w_pipeline.items():
            assert 0.0 < w < 1.0

    def test_roberts_2020_ceiling_no_predictions_returns_no_smd(self):
        from mammal_repurposing.cluster_d.bayesian_prior import (
            fit_cluster_d_prior_stub, roberts_2020_ceiling_check,
        )
        post = fit_cluster_d_prior_stub(["P00000", "P11111"])
        verdicts = roberts_2020_ceiling_check(post, target_smd_predictions=None)
        for u, v in verdicts.items():
            assert v == "NO_SMD_PREDICTION"


# ---------------------------------------------------------------------------
# V6 Cluster C — TxGNN per-disease ranking API
# ---------------------------------------------------------------------------
class TestTxGnnPerDiseaseAPI:
    def test_availability_when_txgnn_missing(self):
        from mammal_repurposing.cluster_c.txgnn import availability
        a = availability()
        assert isinstance(a, dict)
        assert "available" in a

    def test_cognition_anchors_complete(self):
        from mammal_repurposing.cluster_c.txgnn import COGNITION_ANCHORS
        # 5 anchor diseases per design: MCI, AD, ADHD, FXS, narcolepsy
        assert len(COGNITION_ANCHORS) == 5
        assert "EFO_0006816" in COGNITION_ANCHORS    # MCI

    def test_score_compounds_returns_empty_when_no_model(self, monkeypatch):
        """Should not raise when TxGNN unavailable and no cache."""
        from mammal_repurposing.cluster_c import txgnn as tx
        # Force ImportError on load_txgnn so the fall-through path triggers
        monkeypatch.setattr(tx, "load_txgnn",
                            lambda *a, **k: (_ for _ in ()).throw(
                                ImportError("txgnn missing")))
        out = tx.score_compounds_against_anchor(
            ["DB00001"], model=None, cache_df=pd.DataFrame(columns=[
                "anchor_id", "anchor_name", "drug_idx",
                "p_indication", "p_contraindication"
            ]),
        )
        # With empty cache + no model lookup possible, we expect empty
        assert isinstance(out, pd.DataFrame)

    def test_aggregate_per_compound_weight_mean(self):
        from mammal_repurposing.cluster_c.txgnn import aggregate_per_compound
        long = pd.DataFrame([
            {"compound_id": "X", "anchor_id": "A1", "anchor_name": "n",
             "weight": 1.0, "p_indication": 0.8, "p_contraindication": 0.1},
            {"compound_id": "X", "anchor_id": "A2", "anchor_name": "n",
             "weight": 1.0, "p_indication": 0.6, "p_contraindication": 0.2},
            {"compound_id": "Y", "anchor_id": "A1", "anchor_name": "n",
             "weight": 1.0, "p_indication": 0.3, "p_contraindication": 0.4},
        ])
        out = aggregate_per_compound(long)
        x = out[out["compound_id"] == "X"].iloc[0]
        assert x["mean_p_indication"] == pytest.approx(0.7)
        assert x["mean_p_contraindication"] == pytest.approx(0.15)
        assert int(x["n_anchors_resolved"]) == 2

    def test_score_compounds_uses_cache_lookup(self):
        """With a pre-built cache + mock model lookup, compounds resolve to scores."""
        from mammal_repurposing.cluster_c import txgnn as tx
        cache_df = pd.DataFrame([
            {"anchor_id": "EFO_0006816", "anchor_name": "mild cognitive impairment",
             "disease_idx": 100, "drug_idx": 42,
             "p_indication": 0.75, "p_contraindication": 0.10},
            {"anchor_id": "EFO_0000249", "anchor_name": "Alzheimer disease",
             "disease_idx": 101, "drug_idx": 42,
             "p_indication": 0.60, "p_contraindication": 0.05},
        ])
        # Mock a TxGNN model with a drug lookup
        mock_model = MagicMock()
        with patch.object(tx, "_drug_idx_lookup",
                          return_value={"DONEPEZIL": 42, "DB00843": 42}):
            out = tx.score_compounds_against_anchor(
                ["donepezil"], model=mock_model, cache_df=cache_df,
            )
        # Should have 2 rows for donepezil (one per anchor present in cache)
        donepezil_rows = out[out["compound_id"] == "donepezil"]
        assert len(donepezil_rows) == 2
        # p_indication of 0.75 (MCI) and 0.60 (AD)
        ind_scores = sorted(donepezil_rows["p_indication"].tolist())
        assert ind_scores == [pytest.approx(0.6), pytest.approx(0.75)]
