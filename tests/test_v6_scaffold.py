"""Pytest coverage for V6 scaffold modules + §7.16 detector ensemble:

  - §7.16 pockets/detector_ensemble.py — stub mode + consensus voting
  - V6 diagnostics/per_head_bias.py — PC/SN/OOD signature + trust matrix
  - V6 fusion/bayesian_router.py — OOD + confidence gates + routing
  - V6 cluster_d/bayesian_prior.py — stub posterior + reference anchors
  - V6 cluster_d/data_fetchers.py — availability probes + AHBA stub
"""

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
# §7.16 Detector ensemble
# ---------------------------------------------------------------------------
class TestDetectorEnsemble:
    def test_availability_returns_3_keys(self):
        from mammal_repurposing.pockets.detector_ensemble import detector_availability
        a = detector_availability()
        assert set(a.keys()) == {"p2rank", "pocketminer", "cryptobench"}
        for v in a.values():
            assert isinstance(v, bool)

    def test_stub_detectors_return_pockets(self, tmp_path):
        from mammal_repurposing.pockets.detector_ensemble import run_ensemble
        # Any path works in stub mode
        fake_pdb = tmp_path / "fake.pdb"
        fake_pdb.write_text("HEADER\n", encoding="utf-8")
        out = run_ensemble(fake_pdb)
        # At least p2rank should return its STUB pocket
        assert "p2rank" in out
        assert len(out["p2rank"]) >= 1

    def test_consensus_vote_distance_threshold(self):
        from mammal_repurposing.pockets.detector_ensemble import (
            DetectedPocket, consensus_vote,
        )
        # Three detectors all reporting one pocket at ~origin
        det = {
            "p2rank": [DetectedPocket("p2rank", "p1", np.zeros(3), 0.8)],
            "pocketminer": [DetectedPocket("pocketminer", "pm1",
                                            np.array([1.0, 0, 0]), 0.7,
                                            is_cryptic_predicted=True)],
            "cryptobench": [DetectedPocket("cryptobench", "cb1",
                                            np.array([2.0, 0, 0]), 0.6,
                                            is_cryptic_predicted=True)],
        }
        verdicts = consensus_vote(det, distance_threshold=8.0)
        assert len(verdicts) == 1
        v = verdicts[0]
        assert v.n_voters == 3
        assert v.consensus_class == "cryptic_consensus"
        assert v.is_cryptic

    def test_consensus_vote_far_pockets_are_separate(self):
        from mammal_repurposing.pockets.detector_ensemble import (
            DetectedPocket, consensus_vote,
        )
        det = {
            "p2rank": [DetectedPocket("p2rank", "p1", np.zeros(3), 0.8)],
            "pocketminer": [DetectedPocket("pocketminer", "pm1",
                                            np.array([50.0, 0, 0]), 0.7)],
        }
        verdicts = consensus_vote(det, distance_threshold=8.0)
        assert len(verdicts) == 2
        for v in verdicts:
            assert v.n_voters == 1
            assert v.consensus_class == "single_detector"


# ---------------------------------------------------------------------------
# V6 per-head bias signature + trust matrix
# ---------------------------------------------------------------------------
class TestPerHeadBias:
    def test_pc_ratio_detects_collapse(self):
        from mammal_repurposing.diagnostics.per_head_bias import compute_pc_ratio
        # MAMMAL-like prior collapse: σ ≈ 0.1 against training SD 1.34
        rng = np.random.default_rng(0)
        collapsed = rng.normal(5.8, 0.1, 1000)
        ratio = compute_pc_ratio(collapsed, training_label_std=1.34)
        assert ratio < 0.3      # SEVERE
        # Healthy predictions
        healthy = rng.normal(5.8, 1.2, 1000)
        ratio_h = compute_pc_ratio(healthy, training_label_std=1.34)
        assert ratio_h > 0.5    # ACCEPTABLE

    def test_pc_severity_buckets(self):
        from mammal_repurposing.diagnostics.per_head_bias import HeadBiasSignature
        s_sev = HeadBiasSignature("X", "T", 10, pc_ratio=0.2, sn_rho=0, ood_fraction=0, calibration_tier="C")
        s_mod = HeadBiasSignature("X", "T", 10, pc_ratio=0.4, sn_rho=0, ood_fraction=0, calibration_tier="C")
        s_ok = HeadBiasSignature("X", "T", 10, pc_ratio=0.8, sn_rho=0, ood_fraction=0, calibration_tier="C")
        assert s_sev.pc_severity == "SEVERE"
        assert s_mod.pc_severity == "MODERATE"
        assert s_ok.pc_severity == "ACCEPTABLE"

    def test_sn_bias_high_for_correlated(self):
        from mammal_repurposing.diagnostics.per_head_bias import compute_sn_bias
        rng = np.random.default_rng(7)
        x = rng.uniform(0, 1, 20)
        # Predictions strongly correlated with Tanimoto (similarity searcher)
        sim = x + rng.normal(0, 0.05, 20)
        sn = compute_sn_bias(x, sim)
        assert sn > 0.85

    def test_ood_fraction_extremes(self):
        from mammal_repurposing.diagnostics.per_head_bias import compute_ood_fraction
        rng = np.random.default_rng(2)
        T = rng.normal(0, 1, (50, 4))
        # Query near training mean → OOD ~0
        Q_in = rng.normal(0, 0.5, (20, 4))
        ood_in = compute_ood_fraction(Q_in, T)
        assert ood_in < 0.5
        # Query far from training → OOD high
        Q_out = rng.normal(8, 0.1, (20, 4))
        ood_out = compute_ood_fraction(Q_out, T)
        assert ood_out >= ood_in

    def test_trust_matrix_softmax_renormalises(self):
        from mammal_repurposing.diagnostics.per_head_bias import (
            HeadBiasSignature, build_trust_matrix,
        )
        sigs = [
            HeadBiasSignature("MAMMAL", "T1", 10, 0.1, 0.0, 0.5, "D"),
            HeadBiasSignature("Tanimoto", "T1", 10, 1.0, 0.9, 0.0, "A"),
            HeadBiasSignature("MMAtt-DTA", "T1", 10, 0.8, 0.4, 0.2, "A"),
        ]
        T = build_trust_matrix(sigs)
        assert "T1" in T.index
        # Rows sum to ~1
        assert abs(T.loc["T1"].sum() - 1.0) < 1e-6
        # MAMMAL gets least weight (D tier + collapsed)
        assert T.loc["T1"]["MAMMAL"] < T.loc["T1"]["Tanimoto"]
        assert T.loc["T1"]["MAMMAL"] < T.loc["T1"]["MMAtt-DTA"]


# ---------------------------------------------------------------------------
# V6 Bayesian router
# ---------------------------------------------------------------------------
class TestBayesianRouter:
    def test_gate_ood_decays_beyond_threshold(self):
        from mammal_repurposing.fusion.bayesian_router import gate_ood
        # Inside threshold → ~1
        assert abs(gate_ood(5.0, 10.0) - 1.0) < 1e-6
        # Far beyond → ~0
        assert gate_ood(20.0, 10.0) < 0.01

    def test_gate_confidence_inverse_var(self):
        from mammal_repurposing.fusion.bayesian_router import gate_confidence
        # Tight σ → high weight; loose σ → low weight
        h_tight = gate_confidence(0.1)
        h_loose = gate_confidence(5.0)
        assert h_tight > h_loose

    def test_route_one_collapses_with_zero_trust(self):
        from mammal_repurposing.fusion.bayesian_router import route_one
        r = route_one("c1", "T1",
                      head_predictions={"M": 7.0, "T": 8.0},
                      head_sigmas={"M": 0.5, "T": 0.5},
                      trust_row={"M": 0.0, "T": 0.0})
        assert np.isnan(r.y_hat)
        assert r.n_active_heads == 0
        assert "all heads gated to zero" in r.note

    def test_route_one_picks_higher_trust(self):
        from mammal_repurposing.fusion.bayesian_router import route_one
        r = route_one("c1", "T1",
                      head_predictions={"M": 4.0, "T": 8.0},
                      head_sigmas={"M": 0.5, "T": 0.5},
                      trust_row={"M": 0.1, "T": 0.9})
        # Ensemble pulled toward T (high trust)
        assert r.y_hat > 7.0
        # Weights normalised
        assert abs(sum(r.weights.values()) - 1.0) < 1e-6

    def test_identifiability_theorem_n_star(self):
        from mammal_repurposing.fusion.bayesian_router import identifiability_diagnostic
        import pandas as pd
        T = pd.DataFrame([[0.5, 0.5]],
                         index=["TARGET_A"], columns=["MAMMAL", "Tanimoto"])
        diag = identifiability_diagnostic(T, {"TARGET_A": 20})
        # n* should be ~720 per the theorem (4·σ²/0.01)
        assert diag["n_star_for_id"].iloc[0] >= 700
        # 20 << 720 → not identifiable
        assert not diag["identifiable"].iloc[0]


# ---------------------------------------------------------------------------
# V6 Cluster D Bayesian prior (stub path)
# ---------------------------------------------------------------------------
class TestClusterDPrior:
    def test_stub_returns_uniform_when_no_inputs(self):
        from mammal_repurposing.cluster_d.bayesian_prior import (
            fit_cluster_d_prior_stub,
        )
        post = fit_cluster_d_prior_stub(["T1", "T2", "T3"])
        # No streams → theta=0 → sigmoid(0) = 0.5
        for t in ["T1", "T2", "T3"]:
            assert abs(post.w_pipeline[t] - 0.5) < 1e-6
        assert post.method == "stage_0_stub"

    def test_stub_respects_reference_anchors(self):
        from mammal_repurposing.cluster_d.bayesian_prior import (
            fit_cluster_d_prior_stub, DEFAULT_ANCHORS,
        )
        # BDNF should be in anchors; passing it shifts θ toward 0.7
        post = fit_cluster_d_prior_stub(
            target_uniprots=["BDNF"],
            y_ahba={"BDNF": 0.65},
        )
        # The reference-anchor blend gives θ = 0.5*z_score + 0.5*0.7
        # With single-target, z-score is 0 (no variance), so θ ≈ 0.35
        assert "BDNF" in post.theta_mean

    def test_roberts_ceiling_check(self):
        from mammal_repurposing.cluster_d.bayesian_prior import (
            ClusterDPosterior, roberts_2020_ceiling_check,
        )
        post = ClusterDPosterior(
            targets=["T1", "T2"],
            theta_mean={"T1": 0.3, "T2": 0.5},
            w_pipeline={"T1": 0.57, "T2": 0.62},
        )
        # No SMD predictions → all NO_SMD_PREDICTION
        gate = roberts_2020_ceiling_check(post)
        assert gate["T1"] == "NO_SMD_PREDICTION"
        # With SMD predictions, one inside / one outside ceiling
        gate2 = roberts_2020_ceiling_check(post,
                                            target_smd_predictions={"T1": 0.3, "T2": 0.7},
                                            smd_ceiling=0.5)
        assert gate2["T1"] == "REGIME_OK"
        assert gate2["T2"] == "REGIME_VIOLATION"


# ---------------------------------------------------------------------------
# V6 Cluster D data fetchers
# ---------------------------------------------------------------------------
class TestClusterDFetchers:
    def test_availability_returns_3_keys(self):
        from mammal_repurposing.cluster_d.data_fetchers import availability
        a = availability()
        assert set(a.keys()) == {"abagen", "ot_genetics_reachable", "tiledbsoma"}

    def test_ahba_stub_falls_back_to_brain_region(self):
        from mammal_repurposing.cluster_d.data_fetchers import fetch_ahba_expression
        # When abagen not installed, returns stub from §8.6 BRAIN_REGION_BIAS
        out = fetch_ahba_expression(["Q01959", "P36544"])
        assert "Q01959" in out
        # SLC6A3 brainstem bias → cortical_r should be -0.2
        assert out["Q01959"].gene_symbol == "SLC6A3"

    def test_unknown_target_returns_empty_signature(self):
        from mammal_repurposing.cluster_d.data_fetchers import fetch_ahba_expression
        out = fetch_ahba_expression(["UNKNOWN_TARGET"])
        assert "UNKNOWN_TARGET" in out
        # Empty gene_symbol when not in BRAIN_REGION_BIAS
        assert out["UNKNOWN_TARGET"].gene_symbol == ""
