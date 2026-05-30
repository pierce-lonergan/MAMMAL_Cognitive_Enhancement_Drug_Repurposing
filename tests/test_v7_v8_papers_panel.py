"""V7.2 Stage 2 + V8.4 Gate 1 dry-run + V6.A/V7 paper drafts + V6.B.5 panel.

Covers:
- V7.2 Stage 2 per-(class, endpoint) subdomain priors + NUTS integration
- V8.4 Gate 1 dry-run driver (synthetic phenotype generation, clustering,
  AMI/ARI verdict thresholds)
- V6.A + V7 paper draft structural validators
- V6.B.5 panel expansion (~210 targets, subset constraints, dedup)
"""

from __future__ import annotations

import sys
import importlib.util
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
# V7.2 Stage 2 — per-(class, endpoint) subdomain priors
# ---------------------------------------------------------------------------
class TestV7Stage2SubdomainPriors:
    def test_per_subdomain_priors_table_populated(self):
        from mammal_repurposing.translation.prisma_priors import (
            PER_SUBDOMAIN_PRIORS,
        )
        # At least 25 (class, endpoint) cells per V7.2 Stage 2 spec
        assert len(PER_SUBDOMAIN_PRIORS) >= 25

    def test_subdomain_priors_respect_roberts_ceiling(self):
        """No (class, endpoint) prior may exceed Roberts 2020 g=0.50 ceiling."""
        from mammal_repurposing.translation.prisma_priors import (
            PER_SUBDOMAIN_PRIORS,
        )
        for key, (g_mean, g_sd) in PER_SUBDOMAIN_PRIORS.items():
            assert g_mean <= 0.50, \
                f"{key} prior_g={g_mean} violates Roberts ceiling"

    def test_get_subdomain_prior_returns_known_pair(self):
        from mammal_repurposing.translation.prisma_priors import (
            get_subdomain_prior,
        )
        # AChE-I × delayed_recall = (0.31, 0.14) per Birks 2018
        g_mean, g_sd = get_subdomain_prior("AChE-I", "delayed_recall")
        assert g_mean == pytest.approx(0.31)
        assert g_sd == pytest.approx(0.14)

    def test_get_subdomain_prior_falls_back_to_class(self):
        from mammal_repurposing.translation.prisma_priors import (
            get_subdomain_prior,
        )
        # (AChE-I, unknown_endpoint) → fall back to class prior (0.18, 0.15)
        g_mean, g_sd = get_subdomain_prior("AChE-I", "unknown_endpoint_xyz")
        assert g_mean == pytest.approx(0.18)
        assert g_sd == pytest.approx(0.15)

    def test_list_subdomain_endpoints(self):
        from mammal_repurposing.translation.prisma_priors import (
            list_subdomain_endpoints,
        )
        endpoints = list_subdomain_endpoints("AChE-I")
        assert "delayed_recall" in endpoints
        assert "ADAS-Cog" in endpoints

    def test_subdomain_prior_table_nested_structure(self):
        from mammal_repurposing.translation.prisma_priors import (
            subdomain_prior_table,
        )
        t = subdomain_prior_table()
        assert "AChE-I" in t
        assert "delayed_recall" in t["AChE-I"]
        assert "mean" in t["AChE-I"]["delayed_recall"]
        assert "sd" in t["AChE-I"]["delayed_recall"]

    def test_stub_uses_subdomain_when_available(self):
        """fit_effect_size_stub method='prisma_stub_subdomain' when matching."""
        from mammal_repurposing.translation.effect_size_model import (
            EffectSizeObservation, fit_effect_size_stub,
        )
        obs = EffectSizeObservation(
            compound="donepezil", class_name="AChE-I", target_uniprot="P22303",
            pchembl_post_mean=8.5, relevance_post_mean=0.90,
            endpoint="delayed_recall",   # has subdomain prior
        )
        post = fit_effect_size_stub([obs])
        assert post.method == "prisma_stub_subdomain"

    def test_stub_falls_back_to_class_when_no_subdomain(self):
        from mammal_repurposing.translation.effect_size_model import (
            EffectSizeObservation, fit_effect_size_stub,
        )
        obs = EffectSizeObservation(
            compound="x", class_name="AChE-I", target_uniprot="P22303",
            pchembl_post_mean=8.5, relevance_post_mean=0.90,
            endpoint="never_seen_endpoint_xyz",
        )
        post = fit_effect_size_stub([obs])
        # No subdomain hit → method stays "prisma_stub"
        assert post.method == "prisma_stub"


# ---------------------------------------------------------------------------
# V8.4 Gate 1 dry-run (scripts/60_v8_gate1_dryrun.py)
# ---------------------------------------------------------------------------
class TestV8Gate1Dryrun:
    def _load_driver(self):
        spec = importlib.util.spec_from_file_location(
            "v8_gate1_dryrun", SCRIPTS / "60_v8_gate1_dryrun.py")
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        return m

    def test_synthetic_phenotype_shape(self):
        driver = self._load_driver()
        F, labels, names = driver.generate_synthetic_phenotype(
            n_classes=5, n_per_class=20, K=10,
        )
        assert F.shape == (100, 10)
        assert labels.shape == (100,)
        assert len(names) == 5
        # Each class has exactly n_per_class compounds
        for c in range(5):
            assert (labels == c).sum() == 20

    def test_synthetic_centroids_separable(self):
        """Orthogonal centroids → Agglomerative AMI should be ≈ 1.0."""
        driver = self._load_driver()
        F, labels, _ = driver.generate_synthetic_phenotype(
            n_classes=5, n_per_class=20, K=30, noise_sigma=0.3,
        )
        result = driver.cluster_and_score(
            F, labels, method="agglomerative", n_clusters=5,
        )
        # AMI close to 1.0 (perfect recovery on synthetic separable centroids)
        assert result["ami"] > 0.90

    def test_gate1_verdict_thresholds(self):
        driver = self._load_driver()
        # PASS band
        assert driver.gate1_verdict(ami=0.60, ari=0.50) == "PASS"
        assert driver.gate1_verdict(ami=0.50, ari=0.40) == "PASS"
        # DEGRADE band
        assert driver.gate1_verdict(ami=0.35, ari=0.30) == "DEGRADE"
        # FAIL band
        assert driver.gate1_verdict(ami=0.20, ari=0.15) == "FAIL"

    def test_5_cognition_centroid_names(self):
        driver = self._load_driver()
        names = driver.COGNITION_CENTROID_NAMES
        assert "cholinergic" in names
        assert "catecholaminergic" in names
        assert "glutamatergic" in names
        assert "trophic_ISR" in names
        assert "remyelination" in names


# ---------------------------------------------------------------------------
# V6.A + V7 paper drafts structural validators
# ---------------------------------------------------------------------------
class TestPaperDrafts:
    REQUIRED_SECTIONS = (
        "Title",
        "Abstract",
        "Methods",
        "Results",
        "Discussion",
        "Code + data availability",
        "References",
    )

    def test_v6a_paper_draft_exists_with_required_sections(self):
        path = ROOT / "reports" / "paper-drafts" / "v6a_paper_draft.md"
        assert path.exists(), "V6.A paper draft missing"
        body = path.read_text(encoding="utf-8").lower()
        for section in self.REQUIRED_SECTIONS:
            assert section.lower() in body, \
                f"V6.A paper missing section: {section}"

    def test_v6a_paper_reports_tier_a_fail(self):
        path = ROOT / "reports" / "paper-drafts" / "v6a_paper_draft.md"
        body = path.read_text(encoding="utf-8").lower()
        # Tier-A FAIL at SLC6A3 must be the headline finding
        assert "tier-a fail" in body or "tier-a fails" in body or "tier-a failed" in body or "tier-a criterion: fail" in body or "fails the tier-a criterion" in body
        # MMAtt-DTA ρ +0.65 must be reported
        assert "+0.65" in body or "0.65" in body
        # Tanimoto +0.90 baseline must be referenced
        assert "+0.90" in body or "0.90" in body
        # INVERT-mask architecture
        assert "invert-mask" in body or "invert mask" in body

    def test_v6a_paper_cites_key_dti_methods(self):
        path = ROOT / "reports" / "paper-drafts" / "v6a_paper_draft.md"
        body = path.read_text(encoding="utf-8").lower()
        for citation in ("shoshan", "schulman", "koh", "gorantla", "mervin",
                          "park"):
            assert citation in body, f"V6.A paper missing citation: {citation}"

    def test_v7_paper_draft_exists_with_required_sections(self):
        path = ROOT / "reports" / "paper-drafts" / "v7_paper_draft.md"
        assert path.exists(), "V7 paper draft missing"
        body = path.read_text(encoding="utf-8").lower()
        for section in self.REQUIRED_SECTIONS:
            assert section.lower() in body, \
                f"V7 paper missing section: {section}"

    def test_v7_paper_reports_real_nuts_metrics(self):
        path = ROOT / "reports" / "paper-drafts" / "v7_paper_draft.md"
        body = path.read_text(encoding="utf-8")
        # R̂ = 1.000 must appear
        assert "1.000" in body
        # MAE = 0.073 (Gate 3 PASS metric)
        assert "0.073" in body
        # 0 Roberts ceiling violations
        body_lower = body.lower()
        assert "0 violations" in body_lower or "zero" in body_lower

    def test_v7_paper_cites_pet_anchors_and_pbpk(self):
        path = ROOT / "reports" / "paper-drafts" / "v7_paper_draft.md"
        body = path.read_text(encoding="utf-8").lower()
        for citation in ("bohnen", "volkow", "kapur", "watson",
                          "schmidli", "roberts"):
            assert citation in body, f"V7 paper missing citation: {citation}"

    def test_v7_paper_specifies_cpt_venue(self):
        path = ROOT / "reports" / "paper-drafts" / "v7_paper_draft.md"
        body = path.read_text(encoding="utf-8").lower()
        assert "clinical pharmacology" in body
        # CPT:PSP fallback
        assert "cpt:psp" in body or "cpt: pharmacometrics" in body


# ---------------------------------------------------------------------------
# V6.B.5 panel expansion
# ---------------------------------------------------------------------------
class TestV6b5PanelExpansion:
    def test_availability_reports_panel_size(self):
        from mammal_repurposing.cluster_d.panel_expansion import availability
        a = availability()
        assert a["available"] is True
        assert a["stub_mode"] is True
        # Should produce ~100-250 panel targets
        assert 50 <= a["n_panel_targets_total"] <= 300

    def test_22_panel_strict_subset(self):
        from mammal_repurposing.cluster_d.panel_expansion import (
            build_expanded_panel, validate_panel, PANEL_22_TARGETS,
        )
        df = build_expanded_panel()
        val = validate_panel(df)
        # All 22 panel UniProts in expanded panel
        assert val["v6b_panel_22_subset_ok"] is True
        assert val["n_in_v6b_panel_22"] == 22

    def test_all_uniprots_unique(self):
        from mammal_repurposing.cluster_d.panel_expansion import (
            build_expanded_panel,
        )
        df = build_expanded_panel()
        # No duplicate UniProts
        assert df["uniprot"].nunique() == len(df)

    def test_each_target_has_at_least_one_rule(self):
        from mammal_repurposing.cluster_d.panel_expansion import (
            build_expanded_panel,
        )
        df = build_expanded_panel()
        # Every row must have ≥1 rule fired
        assert (df["n_rules"] >= 1).all()

    def test_rules_distribution_diverse(self):
        from mammal_repurposing.cluster_d.panel_expansion import (
            build_expanded_panel, validate_panel,
        )
        df = build_expanded_panel()
        val = validate_panel(df)
        # At least 5 distinct inclusion rules used (l2g/magma/ahba/sc/lit_otar
        # + v6b/liability anchors)
        assert len(val["rules_distribution"]) >= 5

    def test_v6b5_driver_exit_code(self):
        """Smoke: run the driver and check it exits 0."""
        import subprocess
        result = subprocess.run(
            [sys.executable, str(SCRIPTS / "61_v6b5_panel_expand.py"),
             "--out", str(ROOT / "data" / "results" / "v2"
                          / "panel_expanded_test.parquet"),
             "--report", str(ROOT / "reports" / "pipeline" / "panel_expansion_test.md")],
            capture_output=True, text=True, timeout=60,
        )
        # Driver must exit 0 (all sanity checks pass)
        assert result.returncode == 0, \
            f"V6.B.5 driver failed: stderr={result.stderr[-500:]}"
        # Clean up
        for p in (ROOT / "data" / "results" / "v2"
                  / "panel_expanded_test.parquet",
                  ROOT / "reports" / "pipeline" / "panel_expansion_test.md"):
            if p.exists():
                p.unlink()


# ---------------------------------------------------------------------------
# Integration: V7 NUTS report + V8 Gate 1 + V6.B.5 panel artifacts
# ---------------------------------------------------------------------------
class TestSprintArtifacts:
    def test_v7_nuts_v1_report_present(self):
        path = ROOT / "reports" / "pipeline" / "v7_nuts_v1.md"
        assert path.exists()
        body = path.read_text(encoding="utf-8").lower()
        assert "gate 1" in body
        assert "gate 2" in body
        assert "gate 3" in body

    def test_v8_gate1_dryrun_report_present(self):
        path = ROOT / "reports" / "pipeline" / "v8_gate1_dryrun_v1.md"
        assert path.exists()
        body = path.read_text(encoding="utf-8").lower()
        assert "ami" in body
        assert "ari" in body
        assert "agglomerative" in body

    def test_panel_expansion_report_present(self):
        path = ROOT / "reports" / "pipeline" / "panel_expansion_v1.md"
        assert path.exists()
        body = path.read_text(encoding="utf-8").lower()
        assert "22-target" in body or "22 panel" in body or "22/22" in body
