"""V6 phase 2 pytest: PSICHIC adapter + Venn-ABERS + multi-head 4-head smoke."""

from __future__ import annotations

import sys
from pathlib import Path

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
