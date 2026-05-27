"""V8 advanced pytest:
- V8.2 chemCPA training scaffold + Tanimoto-stub imputation + availability
- V8.3 MOFA+ joint embedding + SVD fallback + variance attribution
- V8.5 V7+V8 joint posterior + 8-cell classification + I_novel score
- V8.6 v10 wet-lab composer + Roberts ceiling filter + 4-axis annotation
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
# V8.2 — chemCPA training scaffold
# ---------------------------------------------------------------------------
class TestChemcpaTrain:
    def test_availability_reports_dependencies(self):
        from mammal_repurposing.cluster_e.chemcpa_train import availability
        a = availability()
        assert isinstance(a, dict)
        assert "torch_backend" in a
        assert "rdkit_backend" in a
        assert "stub_mode_works_without_torch" in a

    def test_morgan_fingerprint_returns_n_bits(self):
        from mammal_repurposing.cluster_e.chemcpa_train import (
            morgan_fingerprint, DEFAULT_MORGAN_NBITS, RDKIT_AVAILABLE,
        )
        if not RDKIT_AVAILABLE:
            pytest.skip("rdkit not installed")
        fp = morgan_fingerprint("CCO")  # ethanol
        assert fp is not None
        assert len(fp) == DEFAULT_MORGAN_NBITS

    def test_morgan_fingerprint_returns_none_on_invalid(self):
        from mammal_repurposing.cluster_e.chemcpa_train import (
            morgan_fingerprint, RDKIT_AVAILABLE,
        )
        if not RDKIT_AVAILABLE:
            pytest.skip("rdkit not installed")
        fp = morgan_fingerprint("not_a_valid_smiles")
        assert fp is None

    def test_tanimoto_self_is_one(self):
        from mammal_repurposing.cluster_e.chemcpa_train import (
            tanimoto_similarity, morgan_fingerprint, RDKIT_AVAILABLE,
        )
        if not RDKIT_AVAILABLE:
            pytest.skip("rdkit not installed")
        fp = morgan_fingerprint("CCO")
        assert tanimoto_similarity(fp, fp) == 1.0

    def test_tanimoto_disjoint_is_zero(self):
        from mammal_repurposing.cluster_e.chemcpa_train import (
            tanimoto_similarity,
        )
        a = np.array([1, 0, 1, 0], dtype=np.uint8)
        b = np.array([0, 1, 0, 1], dtype=np.uint8)
        assert tanimoto_similarity(a, b) == 0.0

    def test_train_chemcpa_raises_without_torch(self, monkeypatch):
        from mammal_repurposing.cluster_e import chemcpa_train as cc
        monkeypatch.setattr(cc, "TORCH_AVAILABLE", False)
        with pytest.raises(ImportError, match="torch required"):
            cc.train_chemcpa(lincs_data=None)

    def test_impute_stub_tanimoto_high_confidence(self):
        from mammal_repurposing.cluster_e.chemcpa_train import (
            impute_signature_tanimoto_stub, RDKIT_AVAILABLE,
        )
        if not RDKIT_AVAILABLE:
            pytest.skip("rdkit not installed")
        n_landmark = 977
        rng = np.random.default_rng(0)
        # Reference: ethanol-like signature
        reference_sigs = {"BRD-001": rng.normal(0, 1, n_landmark),
                          "BRD-002": rng.normal(0, 1, n_landmark)}
        reference_smiles = {"BRD-001": "CCO", "BRD-002": "CCC"}
        # Query is ethanol itself → max Tanimoto = 1.0 → high confidence
        result = impute_signature_tanimoto_stub(
            "CCO", reference_sigs, reference_smiles,
        )
        assert result.max_tanimoto_to_train == pytest.approx(1.0)
        assert result.tau_uncertainty == pytest.approx(1.0)
        assert result.flag == "chemCPA.imputed.high_confidence"

    def test_impute_stub_tanimoto_low_confidence_inflates_tau(self):
        """A novel scaffold (low Tanimoto to refs) should inflate τ_chemCPA."""
        from mammal_repurposing.cluster_e.chemcpa_train import (
            impute_signature_tanimoto_stub, RDKIT_AVAILABLE,
        )
        if not RDKIT_AVAILABLE:
            pytest.skip("rdkit not installed")
        n_landmark = 977
        rng = np.random.default_rng(1)
        # Reference: simple aliphatic
        reference_sigs = {"BRD-001": rng.normal(0, 1, n_landmark)}
        reference_smiles = {"BRD-001": "CCO"}
        # Novel-scaffold query: heteroaromatic
        result = impute_signature_tanimoto_stub(
            "c1ccc2nc3ccccc3nc2c1", reference_sigs, reference_smiles,
        )
        # Heteroaromatic vs ethanol Morgan-FP should be very low Tanimoto
        assert result.max_tanimoto_to_train < 0.3
        assert result.tau_uncertainty > 1.0
        assert "low_confidence" in result.flag

    def test_impute_stub_handles_invalid_query(self):
        from mammal_repurposing.cluster_e.chemcpa_train import (
            impute_signature_tanimoto_stub, RDKIT_AVAILABLE,
        )
        if not RDKIT_AVAILABLE:
            pytest.skip("rdkit not installed")
        result = impute_signature_tanimoto_stub(
            "not_a_smiles", {"X": np.zeros(977)}, {"X": "CCO"},
        )
        assert result.flag == "chemCPA.imputed.invalid_smiles"
        assert result.tau_uncertainty >= 1.0


# ---------------------------------------------------------------------------
# V8.3 — MOFA+ joint embedding
# ---------------------------------------------------------------------------
class TestMofaEmbed:
    def test_availability_reports_7_default_views(self):
        from mammal_repurposing.cluster_e.mofa_embed import availability
        a = availability()
        assert a["available"] is True
        assert a["n_views_default"] == 7
        assert "L1000_zscore" in a["default_view_dims"]

    def test_svd_fallback_on_random_views(self):
        from mammal_repurposing.cluster_e.mofa_embed import fit_mofa_plus, MofaConfig
        rng = np.random.default_rng(0)
        n_compounds = 50
        views = {
            "L1000": rng.normal(0, 1, (n_compounds, 30)),
            "CP_DeepProfiler": rng.normal(0, 1, (n_compounds, 20)),
            "MEA": rng.normal(0, 1, (n_compounds, 10)),
        }
        cfg = MofaConfig(K=5, backend="numpy_svd")
        result = fit_mofa_plus(views, cfg=cfg)
        assert result.method == "numpy_svd_fallback"
        assert result.factor_matrix.shape == (n_compounds, 5)
        assert "L1000" in result.per_view_variance
        assert len(result.per_view_variance["L1000"]) == 5

    def test_view_dimension_mismatch_raises(self):
        from mammal_repurposing.cluster_e.mofa_embed import fit_mofa_plus
        views = {
            "A": np.zeros((10, 5)),
            "B": np.zeros((20, 5)),
        }
        with pytest.raises(ValueError, match="share n_compounds"):
            fit_mofa_plus(views)

    def test_empty_views_raises(self):
        from mammal_repurposing.cluster_e.mofa_embed import fit_mofa_plus
        with pytest.raises(ValueError, match="empty views"):
            fit_mofa_plus({})

    def test_variance_attribution_table_structure(self):
        from mammal_repurposing.cluster_e.mofa_embed import (
            fit_mofa_plus, variance_attribution_table, MofaConfig,
        )
        rng = np.random.default_rng(1)
        n_compounds = 30
        views = {
            "L1000": rng.normal(0, 1, (n_compounds, 20)),
            "CP_DeepProfiler": rng.normal(0, 1, (n_compounds, 15)),
        }
        cfg = MofaConfig(K=4, backend="numpy_svd")
        result = fit_mofa_plus(views, cfg=cfg)
        table = variance_attribution_table(result, top_k_factors=3)
        assert table["K"] == 4
        assert table["n_views"] == 2
        assert "L1000" in table["top_factors_per_view"]
        assert len(table["top_factors_per_view"]["L1000"]) == 3

    def test_cosine_similarity_diagonal_is_one(self):
        from mammal_repurposing.cluster_e.mofa_embed import cosine_similarity_matrix
        rng = np.random.default_rng(2)
        F = rng.normal(0, 1, (10, 5))
        sim = cosine_similarity_matrix(F)
        assert sim.shape == (10, 10)
        # Diagonal must be 1.0 (cosine with self)
        np.testing.assert_allclose(np.diag(sim), 1.0, atol=1e-6)

    def test_cognition_centroid_similarity_zero_when_no_centroid(self):
        from mammal_repurposing.cluster_e.mofa_embed import cognition_centroid_similarity
        rng = np.random.default_rng(3)
        F = rng.normal(0, 1, (10, 5))
        sim = cognition_centroid_similarity(F, {"empty_centroid": []})
        assert "empty_centroid" in sim
        assert np.all(sim["empty_centroid"] == 0.0)


# ---------------------------------------------------------------------------
# V8.5 — V7+V8 joint posterior
# ---------------------------------------------------------------------------
class TestJointPosterior:
    def test_availability_reports_5_centroids(self):
        from mammal_repurposing.cluster_e.joint_phenotype import availability
        a = availability()
        assert a["available"] is True
        assert len(a["5_moa_centroids"]) == 5
        assert a["n_axes"] == 4
        assert a["n_8cell_tags"] == 8

    def test_compute_joint_posterior_basic(self):
        from mammal_repurposing.cluster_e.joint_phenotype import compute_joint_posterior
        post = compute_joint_posterior(
            v6a_pchembl_post={"donepezil": (8.5, 0.3), "encenicline": (7.8, 0.4)},
            v6b_theta_post={"P22303": (0.5, 0.15), "P36544": (0.2, 0.20)},
            v7_g_post={"donepezil": (0.18, 0.10), "encenicline": (0.00, 0.12)},
            v8_phen_cosine={"donepezil": 0.82, "encenicline": 0.05},
            compound_to_target={"donepezil": "P22303", "encenicline": "P36544"},
        )
        by_c = post.by_compound()
        # Both compounds present
        assert "donepezil" in by_c
        assert "encenicline" in by_c
        # Each has 8-cell tag + I_novel + JSD assigned
        for c in ("donepezil", "encenicline"):
            assert by_c[c].eight_cell_tag != ""
            assert np.isfinite(by_c[c].i_novel_score)
            assert np.isfinite(by_c[c].three_way_jsd)

    def test_8cell_tag_canonical_positive(self):
        """High pchembl + high theta + high phen → agreement.all_high"""
        from mammal_repurposing.cluster_e.joint_phenotype import compute_joint_posterior
        post = compute_joint_posterior(
            v6a_pchembl_post={"x": (9.0, 0.2)},      # very high pchembl
            v6b_theta_post={"P22303": (1.0, 0.1)},   # very high theta
            v7_g_post={"x": (0.25, 0.05)},
            v8_phen_cosine={"x": 0.85},              # very high phen
            compound_to_target={"x": "P22303"},
        )
        e = post.by_compound()["x"]
        assert e.eight_cell_tag == "agreement.all_high"

    def test_8cell_tag_target_true_phenotype_failed(self):
        """Encenicline pattern: high pchembl + high theta + low phen"""
        from mammal_repurposing.cluster_e.joint_phenotype import compute_joint_posterior
        post = compute_joint_posterior(
            v6a_pchembl_post={"encenicline": (8.0, 0.3)},
            v6b_theta_post={"P36544": (0.7, 0.10)},
            v7_g_post={"encenicline": (0.0, 0.15)},
            v8_phen_cosine={"encenicline": -0.2},    # negative phen
            compound_to_target={"encenicline": "P36544"},
        )
        e = post.by_compound()["encenicline"]
        assert e.eight_cell_tag == "target_true.phenotype_failed"

    def test_8cell_tag_novel_mechanism_LLH(self):
        """Clemastine pattern: low pchembl + low theta + high phen → (L, L, H)"""
        from mammal_repurposing.cluster_e.joint_phenotype import compute_joint_posterior
        post = compute_joint_posterior(
            v6a_pchembl_post={"clemastine": (5.0, 0.4)},    # weak binding
            v6b_theta_post={"P22303": (-0.5, 0.15)},        # weak Cluster D
            v7_g_post={"clemastine": (0.0, 0.20)},
            v8_phen_cosine={"clemastine": 0.78},            # strong phenotype
            compound_to_target={"clemastine": "P22303"},
        )
        e = post.by_compound()["clemastine"]
        assert e.eight_cell_tag == "phenotype_only.novel_mechanism"

    def test_three_way_jsd_in_valid_range(self):
        from mammal_repurposing.cluster_e.joint_phenotype import (
            compute_joint_posterior, three_way_jsd,
        )
        post = compute_joint_posterior(
            v6a_pchembl_post={"x": (8.0, 0.3), "y": (5.0, 0.5)},
            v6b_theta_post={"P22303": (0.5, 0.15)},
            v7_g_post={"x": (0.2, 0.10)},
            v8_phen_cosine={"x": 0.5, "y": -0.5},
            compound_to_target={"x": "P22303", "y": "P22303"},
        )
        js3 = three_way_jsd(post)
        for c, v in js3.items():
            assert 0.0 <= v <= np.log(3) + 0.1   # bounded by log 3

    def test_i_novel_score_in_unit_interval(self):
        from mammal_repurposing.cluster_e.joint_phenotype import (
            compute_joint_posterior, i_novel_score,
        )
        post = compute_joint_posterior(
            v6a_pchembl_post={f"c{i}": (5.0 + i * 0.1, 0.3) for i in range(8)},
            v6b_theta_post={"P22303": (0.5, 0.15)},
            v7_g_post={f"c{i}": (0.1, 0.10) for i in range(8)},
            v8_phen_cosine={f"c{i}": -0.5 + i * 0.2 for i in range(8)},
            compound_to_target={f"c{i}": "P22303" for i in range(8)},
        )
        novel = i_novel_score(post)
        for c, v in novel.items():
            assert 0.0 <= v <= 1.0

    def test_roberts_ceiling_violation_flag(self):
        from mammal_repurposing.cluster_e.joint_phenotype import (
            compute_joint_posterior, violates_roberts_ceiling,
        )
        # Force a violation by giving very high g + low sd
        post = compute_joint_posterior(
            v6a_pchembl_post={"big_g": (10.0, 0.05)},
            v6b_theta_post={"P22303": (2.0, 0.05)},
            v7_g_post={"big_g": (0.80, 0.05)},     # well above 0.50 ceiling
            v8_phen_cosine={"big_g": 0.95},
            compound_to_target={"big_g": "P22303"},
            roberts_ceiling_g=0.50,
        )
        violations = violates_roberts_ceiling(post, ceiling_g=0.50)
        assert violations["big_g"] is True

    def test_wet_lab_priority_respects_ceiling(self):
        from mammal_repurposing.cluster_e.joint_phenotype import (
            compute_joint_posterior, wet_lab_priority,
        )
        post = compute_joint_posterior(
            v6a_pchembl_post={"safe": (8.0, 0.3), "violator": (10.0, 0.05)},
            v6b_theta_post={"P22303": (0.5, 0.15)},
            v7_g_post={"safe": (0.20, 0.10), "violator": (0.80, 0.05)},
            v8_phen_cosine={"safe": 0.6, "violator": 0.9},
            compound_to_target={"safe": "P22303", "violator": "P22303"},
        )
        prio = wet_lab_priority(post, enforce_ceiling=True)
        # Violator should be zeroed
        assert prio["violator"] == 0.0
        # Safe should get a positive priority
        assert prio["safe"] != 0.0


# ---------------------------------------------------------------------------
# V8.6 — v10 wet-lab composer
# ---------------------------------------------------------------------------
class TestWetLabComposerV10:
    def test_availability_reports_18_columns(self):
        from mammal_repurposing.fusion.joint_composition import availability
        a = availability()
        assert a["available"] is True
        assert a["n_columns_output"] == 18
        assert a["roberts_ceiling_default"] == 0.50

    def test_compose_with_no_inputs_returns_empty(self):
        from mammal_repurposing.fusion.joint_composition import (
            compose_wet_lab_shortlist_v10,
        )
        df = compose_wet_lab_shortlist_v10()
        assert df.empty
        assert "rank" in df.columns

    def test_compose_with_minimal_v6a_only(self):
        from mammal_repurposing.fusion.joint_composition import (
            compose_wet_lab_shortlist_v10,
        )
        v6a = pd.DataFrame([
            {"compound_name": "donepezil", "target_uniprot": "P22303",
             "predicted_pkd": 8.5},
            {"compound_name": "MPH", "target_uniprot": "Q01959",
             "predicted_pkd": 7.0},
        ])
        df = compose_wet_lab_shortlist_v10(v6a_pchembl=v6a)
        assert len(df) == 2
        assert "rank" in df.columns
        # Both rows should have evidence_axes='v6a'
        assert (df["evidence_axes"] == "v6a").all()

    def test_compose_respects_roberts_ceiling_filter(self):
        """A compound whose V7 g_90_upper > 0.50 should rank lower (priority=0)."""
        from mammal_repurposing.fusion.joint_composition import (
            compose_wet_lab_shortlist_v10, JointCompositionConfig,
        )
        from mammal_repurposing.cluster_e.joint_phenotype import (
            compute_joint_posterior,
        )
        v6a = pd.DataFrame([
            {"compound_name": "safe", "target_uniprot": "P22303",
             "predicted_pkd": 8.0},
            {"compound_name": "violator", "target_uniprot": "P22303",
             "predicted_pkd": 9.5},
        ])
        v8_post = compute_joint_posterior(
            v6a_pchembl_post={"safe": (8.0, 0.3), "violator": (9.5, 0.05)},
            v6b_theta_post={"P22303": (0.5, 0.10)},
            v7_g_post={"safe": (0.20, 0.10), "violator": (0.75, 0.05)},
            v8_phen_cosine={"safe": 0.6, "violator": 0.95},
            compound_to_target={"safe": "P22303", "violator": "P22303"},
        )
        cfg = JointCompositionConfig(enforce_roberts_ceiling=True, top_n=10)
        df = compose_wet_lab_shortlist_v10(
            v6a_pchembl=v6a, v8_posterior=v8_post, cfg=cfg,
        )
        # Safe should outrank violator
        safe_row = df[df["compound"] == "safe"].iloc[0]
        viol_row = df[df["compound"] == "violator"].iloc[0]
        assert safe_row["wet_lab_priority"] > viol_row["wet_lab_priority"]
        assert viol_row["roberts_ceiling_ok"] is False or viol_row["roberts_ceiling_ok"] is np.False_
        assert viol_row["wet_lab_priority"] == 0.0

    def test_compose_includes_8cell_tag_when_v8_provided(self):
        from mammal_repurposing.fusion.joint_composition import (
            compose_wet_lab_shortlist_v10,
        )
        from mammal_repurposing.cluster_e.joint_phenotype import (
            compute_joint_posterior,
        )
        v8_post = compute_joint_posterior(
            v6a_pchembl_post={"x": (8.5, 0.3)},
            v6b_theta_post={"P22303": (0.7, 0.10)},
            v7_g_post={"x": (0.20, 0.10)},
            v8_phen_cosine={"x": 0.80},
            compound_to_target={"x": "P22303"},
        )
        df = compose_wet_lab_shortlist_v10(v8_posterior=v8_post)
        assert "eight_cell_tag" in df.columns
        # x should be tagged agreement.all_high
        assert df.iloc[0]["eight_cell_tag"] == "agreement.all_high"

    def test_compose_top_n_truncation(self):
        from mammal_repurposing.fusion.joint_composition import (
            compose_wet_lab_shortlist_v10, JointCompositionConfig,
        )
        # 100 compounds; top-N=10
        v6a = pd.DataFrame([
            {"compound_name": f"c{i}", "target_uniprot": "P22303",
             "predicted_pkd": 5.0 + (i * 0.05)}
            for i in range(100)
        ])
        cfg = JointCompositionConfig(top_n=10, enforce_roberts_ceiling=False)
        df = compose_wet_lab_shortlist_v10(v6a_pchembl=v6a, cfg=cfg)
        assert len(df) == 10

    def test_render_v10_markdown_report_writes_file(self, tmp_path):
        from mammal_repurposing.fusion.joint_composition import (
            compose_wet_lab_shortlist_v10, render_v10_markdown_report,
        )
        v6a = pd.DataFrame([
            {"compound_name": "donepezil", "target_uniprot": "P22303",
             "predicted_pkd": 8.5},
        ])
        df = compose_wet_lab_shortlist_v10(v6a_pchembl=v6a)
        out_path = tmp_path / "v10_report.md"
        render_v10_markdown_report(df, out_path)
        assert out_path.exists()
        body = out_path.read_text(encoding="utf-8")
        assert "Wet-Lab Shortlist v10" in body
        assert "8-cell disagreement legend" in body
        assert "donepezil" in body or "Empty shortlist" in body
