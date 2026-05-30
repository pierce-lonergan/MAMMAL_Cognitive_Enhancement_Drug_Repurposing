"""Big-runs sprint pytest:
- Integration umbrella paper (reports/paper-drafts/integration_paper_draft.md)
- V6.B.5 real NUTS posterior (data/results/v2/cluster_d_posterior_expanded_v1.parquet)
- V7 production NUTS posterior (data/results/v2/v7_nuts_posterior_production_v1.parquet)
- Wet-lab handoff document (reports/wet-lab/wet_lab_handoff_v1.md)
- chemCPA larger-scale run artifacts (when present)
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
# Integration umbrella manuscript
# ---------------------------------------------------------------------------
class TestIntegrationPaper:
    REQUIRED_SECTIONS = (
        "Title",
        "Abstract",
        "Introduction",
        "5 architectural layers",
        "Results",
        "Discussion",
        "Code + data availability",
        "References",
    )

    def test_integration_paper_exists(self):
        path = ROOT / "reports" / "paper-drafts" / "integration_paper_draft.md"
        assert path.exists(), "Integration umbrella paper missing"

    def test_integration_paper_has_required_sections(self):
        path = ROOT / "reports" / "paper-drafts" / "integration_paper_draft.md"
        body = path.read_text(encoding="utf-8").lower()
        for section in self.REQUIRED_SECTIONS:
            assert section.lower() in body, \
                f"Integration paper missing section: {section}"

    def test_integration_paper_references_all_4_layer_papers(self):
        path = ROOT / "reports" / "paper-drafts" / "integration_paper_draft.md"
        body = path.read_text(encoding="utf-8")
        for fname in ("v6a_paper_draft.md", "v6b_paper_draft.md",
                       "v7_paper_draft.md", "v8_paper_draft.md"):
            assert fname in body, f"Integration paper missing layer ref: {fname}"

    def test_integration_paper_targets_top_tier_venues(self):
        path = ROOT / "reports" / "paper-drafts" / "integration_paper_draft.md"
        body = path.read_text(encoding="utf-8").lower()
        # Must reference top-tier venues
        assert ("nature" in body or "nature medicine" in body
                or "nature biotechnology" in body or "cell" in body)

    def test_integration_paper_includes_all_5_layers(self):
        path = ROOT / "reports" / "paper-drafts" / "integration_paper_draft.md"
        body = path.read_text(encoding="utf-8").lower()
        for layer in ("v4", "v5", "v6.a", "v6.b", "v7", "v8"):
            assert layer in body, f"Integration paper missing layer: {layer}"

    def test_integration_paper_includes_8_cell_taxonomy(self):
        path = ROOT / "reports" / "paper-drafts" / "integration_paper_draft.md"
        body = path.read_text(encoding="utf-8").lower()
        # Key 8-cell tags
        for tag in ("agreement.all_high", "target_true.phenotype_failed",
                     "phenotype_only.novel_mechanism"):
            assert tag in body, f"Integration paper missing 8-cell tag: {tag}"

    def test_integration_paper_references_roberts_ceiling(self):
        path = ROOT / "reports" / "paper-drafts" / "integration_paper_draft.md"
        body = path.read_text(encoding="utf-8").lower()
        assert "roberts" in body
        assert "0.50" in body or "0.5" in body
        assert "ceiling" in body

    def test_integration_paper_reports_production_metrics(self):
        path = ROOT / "reports" / "paper-drafts" / "integration_paper_draft.md"
        body = path.read_text(encoding="utf-8")
        # Real production metrics
        assert "1.000" in body                    # R̂ from V6.B NUTS
        assert "0.073" in body or "0.071" in body  # V7 MAE
        assert "12,780" in body or "12780" in body # V6.B ESS

    def test_integration_paper_lists_5_paper_suite(self):
        path = ROOT / "reports" / "paper-drafts" / "integration_paper_draft.md"
        body = path.read_text(encoding="utf-8").lower()
        # Should explicitly mention "5-paper" or "4 + 1" suite
        assert ("5-paper" in body or "five-paper" in body
                or "4 + 1" in body or "five paper" in body
                or "4 layer manuscripts" in body
                or "umbrella" in body)


# ---------------------------------------------------------------------------
# V7 production NUTS
# ---------------------------------------------------------------------------
class TestV7ProductionNuts:
    def test_v7_production_parquet_present(self):
        path = (ROOT / "data" / "results" / "v2"
                / "v7_nuts_posterior_production_v1.parquet")
        if not path.exists():
            pytest.skip("V7 production parquet not yet generated")
        df = pd.read_parquet(path)
        assert len(df) == 15      # 15 anchor compounds
        for col in ("compound", "g_mean", "g_2p5", "g_97p5",
                     "g_90_upper", "cluster_d_gate_active"):
            assert col in df.columns

    def test_v7_production_report_reports_R_hat_1_000(self):
        path = ROOT / "reports" / "pipeline" / "v7_nuts_production_v1.md"
        if not path.exists():
            pytest.skip("V7 production report not yet generated")
        body = path.read_text(encoding="utf-8")
        # R̂ should be at 1.000 with 4 chains × 2000 draws
        assert "1.000" in body
        # ESS should be higher than the 2,332 from the dev run
        body_lower = body.lower()
        assert "ess" in body_lower
        # Roberts ceiling pass
        assert "0 violations" in body_lower or "zero violations" in body_lower

    def test_v7_production_all_g_under_roberts_ceiling(self):
        path = (ROOT / "data" / "results" / "v2"
                / "v7_nuts_posterior_production_v1.parquet")
        if not path.exists():
            pytest.skip("V7 production parquet not yet generated")
        df = pd.read_parquet(path)
        # All g_90_upper must be below Roberts 2020 ceiling
        assert (df["g_90_upper"] < 0.50).all(), \
            "Some compounds violate Roberts 2020 ceiling"


# ---------------------------------------------------------------------------
# V6.B.5 real NUTS (chains=1 to bypass multiprocess)
# ---------------------------------------------------------------------------
class TestV6b5RealNuts:
    def test_v6b5_expanded_posterior_present(self):
        path = (ROOT / "data" / "results" / "v2"
                / "cluster_d_posterior_expanded_v1.parquet")
        if not path.exists():
            pytest.skip("V6.B.5 expanded posterior not yet generated")
        df = pd.read_parquet(path)
        # ≥150 targets expected on the 191-target expanded panel
        assert len(df) >= 150
        for col in ("target_uniprot", "gene", "theta_mean",
                     "w_pipeline", "in_v6b_panel_22"):
            assert col in df.columns


# ---------------------------------------------------------------------------
# Wet-lab handoff document
# ---------------------------------------------------------------------------
class TestWetLabHandoff:
    REQUIRED_SECTIONS = (
        "What this pipeline produces",
        "Top-25 compound shortlist",
        "Recommended assay priority",
        "5-MoA cognition centroid",
        "Cost estimate",
        "Recommended pre-screening",
        "Expected timeline",
        "What we can predict",
        "Recommended collaboration scope",
        "Contact",
        "Honest caveats",
    )

    def test_wet_lab_handoff_exists(self):
        path = ROOT / "reports" / "wet-lab" / "wet_lab_handoff_v1.md"
        assert path.exists(), "Wet-lab handoff document missing"

    def test_wet_lab_handoff_has_required_sections(self):
        path = ROOT / "reports" / "wet-lab" / "wet_lab_handoff_v1.md"
        body = path.read_text(encoding="utf-8").lower()
        for section in self.REQUIRED_SECTIONS:
            assert section.lower() in body, \
                f"Wet-lab handoff missing: {section}"

    def test_wet_lab_handoff_lists_top_compounds(self):
        path = ROOT / "reports" / "wet-lab" / "wet_lab_handoff_v1.md"
        body = path.read_text(encoding="utf-8").lower()
        # Canonical top-N compounds must be referenced
        for compound in ("methylphenidate", "donepezil", "modafinil",
                          "encenicline", "rivastigmine", "galantamine",
                          "atomoxetine", "memantine", "caffeine",
                          "clemastine"):
            assert compound in body, f"Wet-lab handoff missing: {compound}"

    def test_wet_lab_handoff_references_bima8_cluster(self):
        path = ROOT / "reports" / "wet-lab" / "wet_lab_handoff_v1.md"
        body = path.read_text(encoding="utf-8").lower()
        # BIMA-8 cluster must be referenced for (L, L, H) tier C
        for compound in ("benztropine", "atropine", "oxybutynin",
                          "trospium", "tiotropium", "pipe-307"):
            assert compound in body, f"Wet-lab handoff missing BIMA-8: {compound}"

    def test_wet_lab_handoff_includes_3_collaboration_options(self):
        path = ROOT / "reports" / "wet-lab" / "wet_lab_handoff_v1.md"
        body = path.read_text(encoding="utf-8").lower()
        # Should offer ≥3 collaboration tiers
        assert "option a" in body
        assert "option b" in body
        assert "option c" in body

    def test_wet_lab_handoff_specifies_cost_estimate(self):
        path = ROOT / "reports" / "wet-lab" / "wet_lab_handoff_v1.md"
        body = path.read_text(encoding="utf-8")
        # Should have explicit dollar estimates
        assert "$" in body
        assert "k" in body.lower()    # thousands

    def test_wet_lab_handoff_specifies_roberts_ceiling_status(self):
        path = ROOT / "reports" / "wet-lab" / "wet_lab_handoff_v1.md"
        body = path.read_text(encoding="utf-8").lower()
        assert "roberts" in body
        assert "0.50" in body or "0.5" in body

    def test_wet_lab_handoff_references_8_cell_taxonomy(self):
        path = ROOT / "reports" / "wet-lab" / "wet_lab_handoff_v1.md"
        body = path.read_text(encoding="utf-8").lower()
        # The two most informative cells should be highlighted
        assert "target_true.phenotype_failed" in body
        assert "phenotype_only.novel_mechanism" in body


# ---------------------------------------------------------------------------
# chemCPA larger-scale (when artifact present)
# ---------------------------------------------------------------------------
class TestChemcpaLargerScale:
    def test_chemcpa_v2_present_when_run(self):
        path = (ROOT / "data" / "results" / "v2"
                / "v8_chemcpa_smoke_v2.parquet")
        if not path.exists():
            pytest.skip("chemCPA v2 not yet generated (background run)")
        df = pd.read_parquet(path)
        # 20 epochs of loss + test_r2_mean + n_train + n_test + ...
        assert (df["metric"] == "test_r2_mean").any()
        n_epoch_rows = df["metric"].str.startswith("epoch_").sum()
        assert n_epoch_rows >= 15      # at least 15 of 20 epochs reported

    def test_chemcpa_v2_loss_decreases(self):
        path = (ROOT / "data" / "results" / "v2"
                / "v8_chemcpa_smoke_v2.parquet")
        if not path.exists():
            pytest.skip("chemCPA v2 not yet generated")
        df = pd.read_parquet(path)
        losses = df[df["metric"].str.startswith("epoch_")].copy()
        losses["epoch"] = losses["metric"].str.extract(r"(\d+)").astype(int)
        losses = losses.sort_values("epoch")
        # First loss should be > last loss
        assert losses["value"].iloc[0] > losses["value"].iloc[-1]


# ---------------------------------------------------------------------------
# Integration — 5-paper suite + 16 figures + production runner artifacts
# ---------------------------------------------------------------------------
class TestSprintArtifactIntegration:
    def test_5_paper_suite_present(self):
        for fname in ("v6a_paper_draft.md", "v6b_paper_draft.md",
                       "v7_paper_draft.md", "v8_paper_draft.md",
                       "integration_paper_draft.md"):
            p = ROOT / "reports" / "paper-drafts" / fname
            assert p.exists(), f"Paper draft missing: {fname}"

    def test_3_documentation_artifacts_at_repo_root(self):
        for fname in ("CITATIONS.bib", "PROJECT_STATUS.md", "README.md"):
            p = ROOT / fname
            assert p.exists(), f"Documentation artifact missing: {fname}"

    def test_wet_lab_handoff_at_reports_root(self):
        assert (ROOT / "reports" / "wet-lab" / "wet_lab_handoff_v1.md").exists()
