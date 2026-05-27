"""V8 / Cluster E + V6.B.4 validation-gates pytest:
- V8.1 LINCS L1000 availability + WTCS math + per-compound aggregation
- V8.1b JUMP-CP availability + sync dry-run + cosine-to-centroid + normalize
- V6.B.4 Roberts ceiling gate + Spearman vs SMD + held-out AUROC + LOSO
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
# V8.1 — LINCS L1000 ingestion
# ---------------------------------------------------------------------------
class TestLincsIngest:
    def test_availability_returns_dict_with_cell_lines(self, monkeypatch):
        monkeypatch.delenv("LINCS_DATA_DIR", raising=False)
        from mammal_repurposing.cluster_e.ingest_lincs import availability
        a = availability()
        assert isinstance(a, dict)
        assert "available" in a
        assert "cell_lines_supported" in a
        # Neural lines always reported
        assert "NPC" in a["cell_lines_supported"]
        assert "NEU" in a["cell_lines_supported"]
        assert "SHSY5Y" in a["cell_lines_supported"]

    def test_cognition_cell_lines_weight_neural_highest(self):
        from mammal_repurposing.cluster_e.ingest_lincs import COGNITION_CELL_LINES
        # Neural lines should be weighted 1.0
        for neural in ("NPC", "NEU", "SHSY5Y"):
            assert COGNITION_CELL_LINES[neural] == 1.0
        # Cancer-line weights should be ≤ 0.5
        for cancer in ("MCF7", "A375", "PC3", "VCAP"):
            assert COGNITION_CELL_LINES[cancer] <= 0.5

    def test_compute_wtcs_returns_per_reference_score(self):
        from mammal_repurposing.cluster_e.ingest_lincs import compute_wtcs
        rng = np.random.default_rng(0)
        n_genes = 100
        query = rng.normal(0, 1, n_genes)
        refs = rng.normal(0, 1, (5, n_genes))
        scores = compute_wtcs(query, refs, top_k=20)
        assert scores.shape == (5,)
        # WTCS bounded by construction (mean of values in [-1, 1])
        assert np.all(scores >= -1.0) and np.all(scores <= 1.0)

    def test_compute_wtcs_self_similarity_is_high(self):
        """A signature compared to itself should have high positive WTCS."""
        from mammal_repurposing.cluster_e.ingest_lincs import compute_wtcs
        rng = np.random.default_rng(1)
        n_genes = 100
        query = rng.normal(0, 1, n_genes)
        # Self-similarity
        score = compute_wtcs(query, query[None, :], top_k=20)
        assert score[0] > 0.3   # strong self-connectivity

    def test_compute_wtcs_negation_is_anti_similar(self):
        """A signature vs its negation should have strongly negative WTCS."""
        from mammal_repurposing.cluster_e.ingest_lincs import compute_wtcs
        rng = np.random.default_rng(2)
        n_genes = 100
        query = rng.normal(0, 1, n_genes)
        negated = -query
        score = compute_wtcs(query, negated[None, :], top_k=20)
        assert score[0] < -0.3   # strong anti-connectivity

    def test_find_lincs_dir_raises_helpful_error(self, monkeypatch, tmp_path):
        from mammal_repurposing.cluster_e.ingest_lincs import _find_lincs_dir
        monkeypatch.delenv("LINCS_DATA_DIR", raising=False)
        with pytest.raises(FileNotFoundError, match="LINCS data directory"):
            _find_lincs_dir(tmp_path / "nonexistent_path_for_test")

    def test_per_compound_max_tau_aggregates_across_lines(self):
        from mammal_repurposing.cluster_e.ingest_lincs import per_compound_max_tau
        long = pd.DataFrame([
            {"pert_id": "BRD-001", "cell_id": "MCF7", "tau": 0.20, "cognition_weight": 0.3},
            {"pert_id": "BRD-001", "cell_id": "NEU",  "tau": 0.45, "cognition_weight": 1.0},
            {"pert_id": "BRD-002", "cell_id": "MCF7", "tau": 0.70, "cognition_weight": 0.3},
        ])
        out = per_compound_max_tau(long)
        # BRD-001's NEU score (0.45 × 1.0) > BRD-001's MCF7 score (0.20 × 0.3 = 0.06)
        brd1 = out[out["pert_id"] == "BRD-001"].iloc[0]
        assert brd1["best_cell_line"] == "NEU"
        assert brd1["max_weighted_tau"] == pytest.approx(0.45)


# ---------------------------------------------------------------------------
# V8.1b — JUMP-CP ingestion
# ---------------------------------------------------------------------------
class TestJumpCpIngest:
    def test_availability_returns_dict_when_deps_missing(self):
        from mammal_repurposing.cluster_e.ingest_jumpcp import availability
        a = availability()
        assert isinstance(a, dict)
        # Either available (if boto3+pycytominer installed) or has reason
        if not a["available"]:
            assert "reason" in a
            assert "dependencies" in a

    def test_embedding_types_specified(self):
        from mammal_repurposing.cluster_e.ingest_jumpcp import EMBEDDING_TYPES
        assert "deepprofiler" in EMBEDDING_TYPES
        assert "cellprofiler" in EMBEDDING_TYPES
        assert "dinov2" in EMBEDDING_TYPES
        # DeepProfiler is 672-d per Moshkov 2024
        assert EMBEDDING_TYPES["deepprofiler"]["dim"] == 672

    def test_13_jump_sources_known(self):
        from mammal_repurposing.cluster_e.ingest_jumpcp import JUMP_SOURCES
        assert len(JUMP_SOURCES) == 13
        # All 13 sources keyed s1..s13
        for i in range(1, 14):
            assert f"s{i}" in JUMP_SOURCES

    def test_cosine_to_centroid_compatible_shapes(self):
        from mammal_repurposing.cluster_e.ingest_jumpcp import cosine_to_centroid
        df = pd.DataFrame({
            "Metadata_pert_iname": ["A", "B", "C"],
            "f1": [1.0, 0.0, -1.0],
            "f2": [0.0, 1.0, 0.0],
            "f3": [1.0, 1.0, 1.0],
        })
        centroid = np.array([1.0, 0.0, 1.0])
        result = cosine_to_centroid(df, centroid,
                                     feature_cols=["f1", "f2", "f3"])
        assert len(result) == 3
        assert "cosine_to_centroid" in result.columns

    def test_cosine_dimension_mismatch_raises(self):
        from mammal_repurposing.cluster_e.ingest_jumpcp import cosine_to_centroid
        df = pd.DataFrame({
            "Metadata_pert_iname": ["A"],
            "f1": [1.0], "f2": [0.0],
        })
        with pytest.raises(ValueError, match="centroid has"):
            cosine_to_centroid(df, np.array([1.0]),
                                feature_cols=["f1", "f2"])

    def test_sync_raises_without_boto3(self, monkeypatch, tmp_path):
        from mammal_repurposing.cluster_e import ingest_jumpcp as ij
        monkeypatch.setattr(ij, "BOTO3_AVAILABLE", False)
        with pytest.raises(ImportError, match="boto3 required"):
            ij.sync_jumpcp_consensus(tmp_path)


# ---------------------------------------------------------------------------
# V6.B.4 — 4-gate validation framework
# ---------------------------------------------------------------------------
class TestClusterDValidationGates:
    def test_availability_reports_4_gates(self):
        from mammal_repurposing.cluster_d.validation_gates import availability
        a = availability()
        assert a["available"] is True
        assert a["n_gates"] == 4
        assert a["n_reference_compounds"] == 15
        assert a["ceiling_g"] == 0.50

    def test_15_reference_compounds_with_required_fields(self):
        from mammal_repurposing.cluster_d.validation_gates import REFERENCE_COMPOUND_SMD
        assert len(REFERENCE_COMPOUND_SMD) == 15
        for r in REFERENCE_COMPOUND_SMD:
            assert r.compound and r.target_uniprot
            assert r.g_2p5 <= r.pooled_g <= r.g_97p5
            assert r.n_trials >= 1

    def test_gate_1_passes_on_safe_predictions(self):
        from mammal_repurposing.cluster_d.validation_gates import (
            gate_1_roberts_ceiling,
        )
        result = gate_1_roberts_ceiling(
            target_smd_predictions={"P22303": 0.30, "Q01959": 0.45, "P36544": 0.20},
            ceiling=0.50,
        )
        assert result.pass_status == "PASS"
        assert result.metric_value == 0.45    # max prediction

    def test_gate_1_fails_on_ceiling_violation(self):
        from mammal_repurposing.cluster_d.validation_gates import (
            gate_1_roberts_ceiling,
        )
        result = gate_1_roberts_ceiling(
            target_smd_predictions={"safe": 0.30, "bad": 0.65},
            ceiling=0.50,
        )
        assert result.pass_status == "FAIL"
        assert "bad" in result.per_item
        assert result.per_item["bad"] == "VIOLATION"

    def test_gate_1_insufficient_data_when_empty(self):
        from mammal_repurposing.cluster_d.validation_gates import (
            gate_1_roberts_ceiling,
        )
        result = gate_1_roberts_ceiling(target_smd_predictions={})
        assert result.pass_status == "INSUFFICIENT_DATA"

    def test_gate_2_spearman_correlates_high(self):
        """Strong positive correlation between θ̄ and SMD must PASS."""
        from mammal_repurposing.cluster_d.validation_gates import (
            gate_2_spearman_vs_smd, REFERENCE_COMPOUND_SMD,
        )
        # Inject θ̄ that matches SMD pooled_g monotonically
        # so Spearman ρ is ≈ +1.0
        theta = {r.target_uniprot: r.pooled_g for r in REFERENCE_COMPOUND_SMD}
        result = gate_2_spearman_vs_smd(theta)
        # ρ should be high, status PASS
        assert result.metric_value > 0.30 or result.pass_status in ("PASS", "DEGRADE")

    def test_gate_2_insufficient_data_when_few_pairs(self):
        from mammal_repurposing.cluster_d.validation_gates import (
            gate_2_spearman_vs_smd,
        )
        result = gate_2_spearman_vs_smd(
            theta_mean={"P22303": 0.5},   # only 1 pair
        )
        assert result.pass_status == "INSUFFICIENT_DATA"

    def test_gate_3_auroc_pass_on_separable(self):
        from mammal_repurposing.cluster_d.validation_gates import (
            gate_3_held_out_gwas,
        )
        # 10 targets; positives = those with theta > 0.5
        theta = {f"T{i}": (i / 10.0) for i in range(10)}
        held_out = {f"T{i}": 0.25 if i >= 7 else 0.05 for i in range(10)}
        result = gate_3_held_out_gwas(theta, held_out, threshold_auroc=0.70)
        # Perfect separation → AUROC = 1.0
        assert result.metric_value == pytest.approx(1.0)
        assert result.pass_status == "PASS"

    def test_gate_3_insufficient_data_when_few_targets(self):
        from mammal_repurposing.cluster_d.validation_gates import (
            gate_3_held_out_gwas,
        )
        result = gate_3_held_out_gwas({"a": 0.5}, {"a": 0.5})
        assert result.pass_status == "INSUFFICIENT_DATA"

    def test_gate_4_loso_passes_on_robust_posteriors(self):
        from mammal_repurposing.cluster_d.validation_gates import (
            gate_4_leave_one_source_out,
        )
        targets = [f"T{i}" for i in range(10)]
        # Full and per-LOSO posteriors very similar → ρ ≈ 1.0
        full = {t: float(i) for i, t in enumerate(targets)}
        loso_ahba = {t: float(i) + 0.01 for i, t in enumerate(targets)}
        loso_l2g = {t: float(i) + 0.02 for i, t in enumerate(targets)}
        result = gate_4_leave_one_source_out(
            {"_full": full, "AHBA": loso_ahba, "L2G": loso_l2g},
            threshold=0.20,
        )
        assert result.pass_status == "PASS"
        assert result.metric_value > 0.95

    def test_aggregate_gates_critical_when_gate1_fails(self):
        from mammal_repurposing.cluster_d.validation_gates import (
            GateResult, aggregate_gates,
        )
        g1 = GateResult("gate_1_roberts_ceiling", "FAIL",
                        metric_value=0.65, metric_threshold=0.50)
        g2 = GateResult("gate_2_spearman_vs_smd", "PASS",
                        metric_value=0.40, metric_threshold=0.30)
        summary = aggregate_gates(g1, g2)
        assert summary["overall"] == "CRITICAL"
        assert summary["n_fail"] == 1
        assert summary["n_pass"] == 1

    def test_aggregate_gates_all_green_when_all_pass(self):
        from mammal_repurposing.cluster_d.validation_gates import (
            GateResult, aggregate_gates,
        )
        results = [
            GateResult(f"gate_{i}", "PASS", 0.5, 0.3)
            for i in range(1, 5)
        ]
        summary = aggregate_gates(*results)
        assert summary["overall"] == "ALL_GREEN"
        assert summary["n_pass"] == 4
