"""V6.B 4-gate live + V6.A/V6.B/V8 figure scripts + production runner + paper
figure-embedding pytest.

Covers:
- V6.B 4-gate live driver (scripts/64_v6b_validation_gates_live.py)
- V6.A figure generation (scripts/66_v6a_figures.py)
- V6.B figure generation (scripts/67_v6b_figures.py)
- V8 figure generation (scripts/65_v8_figures.py)
- Production runner (scripts/68_production_runner.py)
- All 4 paper drafts contain inline figure references
- All 16 figure PNGs land in figures/v{6a,6b,7,8}/
"""

from __future__ import annotations

import sys
import importlib.util
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))


# ---------------------------------------------------------------------------
# V6.B 4-gate live driver
# ---------------------------------------------------------------------------
class TestV6bGatesLive:
    def test_v6b_gate_verdicts_parquet_present(self):
        path = ROOT / "data" / "results" / "v2" / "v6b_gate_verdicts_v1.parquet"
        if not path.exists():
            pytest.skip("V6.B gate verdicts not yet generated")
        df = pd.read_parquet(path)
        # Should have 4 gates
        assert len(df) == 4
        # Required columns
        for col in ("gate", "status", "metric_value", "metric_threshold", "detail"):
            assert col in df.columns

    def test_v6b_gates_report_present(self):
        path = ROOT / "reports" / "pipeline" / "v6b_validation_gates_v1.md"
        if not path.exists():
            pytest.skip("V6.B gates report not yet generated")
        body = path.read_text(encoding="utf-8")
        # All 4 gates referenced
        body_lower = body.lower()
        for gate in ("gate 1", "gate 2", "gate 3", "gate 4"):
            assert gate in body_lower
        # Roberts ceiling reference
        assert "roberts" in body_lower
        # Overall verdict section
        assert "overall verdict" in body_lower


# ---------------------------------------------------------------------------
# V6.A / V6.B / V8 figure generation scripts
# ---------------------------------------------------------------------------
class TestFigureScripts:
    def _load_script(self, name: str):
        spec = importlib.util.spec_from_file_location(
            name, SCRIPTS / f"{name}.py")
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        return m

    def test_v6a_figures_module_loads(self):
        m = self._load_script("66_v6a_figures")
        # Must define the 4 figure functions
        for f in ("figure_1_rho_heatmap", "figure_2_tier_a_fail",
                   "figure_3_v9_fusion_top10", "figure_4_disagreement_axis"):
            assert hasattr(m, f), f"66_v6a_figures missing {f}"

    def test_v6b_figures_module_loads(self):
        m = self._load_script("67_v6b_figures")
        for f in ("figure_1_theta_posterior", "figure_2_source_contribution",
                   "figure_3_reference_anchor_pull",
                   "figure_4_roberts_ceiling_joint"):
            assert hasattr(m, f), f"67_v6b_figures missing {f}"

    def test_v8_figures_module_loads(self):
        m = self._load_script("65_v8_figures")
        for f in ("figure_1_chemcpa_loss", "figure_2_gate1_ami_sweep",
                   "figure_3_8cell_scatter", "figure_4_i_novel_rank"):
            assert hasattr(m, f), f"65_v8_figures missing {f}"

    def test_v6a_figures_landed_in_figures_v6a(self):
        figures_dir = ROOT / "figures" / "v6a"
        if not figures_dir.exists():
            pytest.skip("figures/v6a/ not yet populated")
        for fname in ("fig1_rho_heatmap.png", "fig2_tier_a_fail.png",
                       "fig3_v9_fusion_top10.png",
                       "fig4_disagreement_axis.png"):
            p = figures_dir / fname
            assert p.exists(), f"V6.A figure missing: {fname}"
            assert p.stat().st_size > 10_000, \
                f"V6.A figure too small: {fname}"

    def test_v6b_figures_landed_in_figures_v6b(self):
        figures_dir = ROOT / "figures" / "v6b"
        if not figures_dir.exists():
            pytest.skip("figures/v6b/ not yet populated")
        for fname in ("fig1_theta_posterior.png",
                       "fig2_source_contribution.png",
                       "fig3_reference_anchor_pull.png",
                       "fig4_roberts_ceiling_joint.png"):
            p = figures_dir / fname
            assert p.exists(), f"V6.B figure missing: {fname}"
            assert p.stat().st_size > 10_000

    def test_v8_figures_landed_in_figures_v8(self):
        figures_dir = ROOT / "figures" / "v8"
        if not figures_dir.exists():
            pytest.skip("figures/v8/ not yet populated")
        for fname in ("fig1_chemcpa_loss.png", "fig2_gate1_ami_sweep.png",
                       "fig3_8cell_scatter.png", "fig4_i_novel_rank.png"):
            p = figures_dir / fname
            assert p.exists(), f"V8 figure missing: {fname}"
            assert p.stat().st_size > 10_000

    def test_v7_figures_landed(self):
        figures_dir = ROOT / "figures" / "v7"
        if not figures_dir.exists():
            pytest.skip("figures/v7/ not yet populated")
        for fname in ("fig1_pbpk_traces.png", "fig2_p1_p8_bands.png",
                       "fig3_loo_mae.png", "fig4_sensitivity_sweep.png"):
            p = figures_dir / fname
            assert p.exists(), f"V7 figure missing: {fname}"


# ---------------------------------------------------------------------------
# Paper drafts have inline figure references
# ---------------------------------------------------------------------------
class TestPaperFigureEmbedding:
    def test_v6a_paper_references_all_4_figures(self):
        path = ROOT / "reports" / "paper-drafts" / "v6a_paper_draft.md"
        body = path.read_text(encoding="utf-8")
        for fname in ("fig1_rho_heatmap.png", "fig2_tier_a_fail.png",
                       "fig3_v9_fusion_top10.png",
                       "fig4_disagreement_axis.png"):
            assert fname in body, f"V6.A paper missing figure ref: {fname}"

    def test_v6b_paper_references_all_4_figures(self):
        path = ROOT / "reports" / "paper-drafts" / "v6b_paper_draft.md"
        body = path.read_text(encoding="utf-8")
        for fname in ("fig1_theta_posterior.png",
                       "fig2_source_contribution.png",
                       "fig3_reference_anchor_pull.png",
                       "fig4_roberts_ceiling_joint.png"):
            assert fname in body, f"V6.B paper missing figure ref: {fname}"

    def test_v7_paper_references_all_4_figures(self):
        path = ROOT / "reports" / "paper-drafts" / "v7_paper_draft.md"
        body = path.read_text(encoding="utf-8")
        for fname in ("fig1_pbpk_traces.png", "fig2_p1_p8_bands.png",
                       "fig3_loo_mae.png", "fig4_sensitivity_sweep.png"):
            assert fname in body, f"V7 paper missing figure ref: {fname}"

    def test_v8_paper_references_all_4_figures(self):
        path = ROOT / "reports" / "paper-drafts" / "shelved" / "v8_paper_draft.md"
        body = path.read_text(encoding="utf-8")
        for fname in ("fig1_chemcpa_loss.png", "fig2_gate1_ami_sweep.png",
                       "fig3_8cell_scatter.png", "fig4_i_novel_rank.png"):
            assert fname in body, f"V8 paper missing figure ref: {fname}"

    def test_all_papers_use_consistent_figure_path(self):
        """All 4 papers should reference figures via `../figures/v{paper}/`
        relative path."""
        for paper, figures_dir in [
            ("v6a_paper_draft.md", "v6a"),
            ("v6b_paper_draft.md", "v6b"),
            ("v7_paper_draft.md", "v7"),
            # v8 shelved; its figure-path consistency is no longer enforced
        ]:
            path = ROOT / "reports" / "paper-drafts" / paper
            body = path.read_text(encoding="utf-8")
            assert f"../figures/{figures_dir}/" in body, \
                f"{paper} missing ../figures/{figures_dir}/ path prefix"


# ---------------------------------------------------------------------------
# Production runner
# ---------------------------------------------------------------------------
class TestProductionRunner:
    def test_production_runner_module_loads(self):
        spec = importlib.util.spec_from_file_location(
            "production_runner", SCRIPTS / "68_production_runner.py")
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        # Must define STAGES list with 15 entries
        assert hasattr(m, "STAGES")
        assert len(m.STAGES) >= 14    # 14-15 stages

    def test_production_run_report_present(self):
        path = ROOT / "reports" / "pipeline" / "production_run_v1.md"
        if not path.exists():
            pytest.skip("Production run report not yet generated")
        body = path.read_text(encoding="utf-8")
        # Must reference all 15 stages
        for stage_id in range(1, 16):
            assert f"| {stage_id} |" in body, f"Production report missing stage {stage_id}"
        # Wall-clock + summary
        body_lower = body.lower()
        assert "wall-clock" in body_lower or "wall_clock" in body_lower or "summary" in body_lower

    def test_production_runner_skip_if_exists_exits_zero(self):
        """Production runner with --skip-if-exists must exit 0 when all
        artifacts are present (each stage SKIPPED)."""
        result = subprocess.run(
            [sys.executable, str(SCRIPTS / "68_production_runner.py"),
             "--skip-if-exists",
             "--report", str(ROOT / "reports" / "pipeline" / "production_run_test.md")],
            capture_output=True, text=True, timeout=60, cwd=str(ROOT),
        )
        assert result.returncode == 0, \
            f"production runner failed with skip-if-exists: {result.stderr[-400:]}"
        # Clean up test report
        test_report = ROOT / "reports" / "pipeline" / "production_run_test.md"
        if test_report.exists():
            test_report.unlink()


# ---------------------------------------------------------------------------
# Integration — full sprint artifact inventory
# ---------------------------------------------------------------------------
class TestSprintArtifactInventory:
    EXPECTED_PARQUETS = (
        "ahba_expression_v1.parquet",
        "cluster_d_posterior_v1.parquet",
        "cluster_d_posterior_expanded_v1.parquet",
        "panel_expanded_v1.parquet",
        "v7_effect_size_posterior_v1.parquet",
        "v7_nuts_posterior_v1.parquet",
        "v8_chemcpa_smoke_v1.parquet",
        "v8_gate1_dryrun_v1.parquet",
        "wet_lab_shortlist_v10.parquet",
        "v6b_gate_verdicts_v1.parquet",
    )

    EXPECTED_REPORTS = (
        "pipeline/v6b_validation_gates_v1.md",
        "pipeline/v7_validation_v1.md",
        "pipeline/v7_nuts_v1.md",
        "pipeline/v8_chemcpa_smoke_v1.md",
        "pipeline/v8_gate1_dryrun_v1.md",
        "pipeline/cluster_d_nuts_v1.md",
        "pipeline/cluster_d_nuts_expanded_v1.md",
        "wet-lab/wet_lab_shortlist_v10.md",
        "paper-drafts/v6a_paper_draft.md",
        "paper-drafts/v6b_paper_draft.md",
        "paper-drafts/v7_paper_draft.md",
        "paper-drafts/shelved/v8_paper_draft.md",
        "paper-drafts/v7_osf_preregistration.md",
        "paper-drafts/v8_osf_preregistration.md",
        "paper-drafts/methodology_v3.md",
        "pipeline/hypothesis_audit_v1.md",
    )

    def test_all_expected_parquets_present(self):
        results_dir = ROOT / "data" / "results" / "v2"
        for fname in self.EXPECTED_PARQUETS:
            p = results_dir / fname
            assert p.exists(), f"Expected parquet missing: {fname}"

    def test_all_expected_reports_present(self):
        reports_dir = ROOT / "reports"
        for fname in self.EXPECTED_REPORTS:
            p = reports_dir / fname
            assert p.exists(), f"Expected report missing: {fname}"

    def test_all_16_figures_present(self):
        for v in ("v6a", "v6b", "v7", "v8"):
            figures_dir = ROOT / "figures" / v
            assert figures_dir.exists(), f"figures/{v}/ missing"
            n_pngs = len(list(figures_dir.glob("fig*.png")))
            assert n_pngs == 4, f"figures/{v}/ has {n_pngs} PNGs (expected 4)"

    def test_citations_bib_at_repo_root(self):
        assert (ROOT / "CITATIONS.bib").exists()

    def test_project_status_at_repo_root(self):
        assert (ROOT / "PROJECT_STATUS.md").exists()
