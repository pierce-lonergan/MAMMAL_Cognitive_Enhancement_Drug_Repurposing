"""V7.4 validation driver pytest + OSF pre-registration structural validator.

Covers:
- scripts/57_v7_validation_gates.py — build_observations correctness;
  PANEL_TARGET_CLASS_MAP completeness; COMPOUND_CLASS/TARGET_OVERRIDE coverage
- P1-P8 compound-name parser handles dose-suffixed names (P3_methylphenidate_20mg)
- OSF pre-registration documents have required sections per OSF + AsPredicted
  templates (hypothesis, model spec, predictions, validation gates,
  sensitivity, falsifiability)
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

# Add scripts/ for the V7 driver import
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))


# ---------------------------------------------------------------------------
# V7.4 driver — scripts/57_v7_validation_gates.py
# ---------------------------------------------------------------------------
class TestV7ValidationDriver:
    def _load_driver(self):
        """Dynamic-import the V7 driver script (digit-prefixed filename)."""
        spec = importlib.util.spec_from_file_location(
            "v7_driver", SCRIPTS / "57_v7_validation_gates.py")
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        return m

    def test_panel_target_class_map_covers_22_targets(self):
        driver = self._load_driver()
        assert len(driver.PANEL_TARGET_CLASS_MAP) >= 22
        # All 22 panel UniProts present
        canonical_panel = {
            "P22303", "P36544", "P42261", "P42262", "P42263", "P48058",
            "Q12879", "Q13224", "P21728", "Q01959", "P08913", "P23975",
            "Q9Y5N1", "O43613", "O43614", "Q08499", "O76083", "Q16620",
            "Q99720", "O43526", "O43525", "O60741",
        }
        assert canonical_panel.issubset(set(driver.PANEL_TARGET_CLASS_MAP.keys()))

    def test_compound_class_override_covers_p1_through_p8_anchors(self):
        driver = self._load_driver()
        # The 8 P1-P8 anchor compound base names
        anchor_compounds = [
            "donepezil", "encenicline", "methylphenidate", "modafinil",
            "memantine", "intepirdine", "pridopidine", "lecanemab",
        ]
        for c in anchor_compounds:
            assert c in driver.COMPOUND_CLASS_OVERRIDE
            assert c in driver.COMPOUND_TARGET_OVERRIDE

    def test_compound_target_override_resolves_to_panel_targets(self):
        driver = self._load_driver()
        # Every target in the override map must be in PANEL_TARGET_CLASS_MAP
        # (which is built from the 22-target cognition panel)
        for compound, target in driver.COMPOUND_TARGET_OVERRIDE.items():
            assert target in driver.PANEL_TARGET_CLASS_MAP, \
                f"{compound} maps to {target} which is not in panel"

    def test_build_observations_handles_empty_v6b(self):
        driver = self._load_driver()
        # V6.A with 2 compounds, V6.B empty (no Cluster D yet)
        v6a = pd.DataFrame([
            {"compound_name": "donepezil", "target_uniprot": "P22303",
             "predicted_pkd": 8.5},
            {"compound_name": "memantine", "target_uniprot": "Q13224",
             "predicted_pkd": 6.0},
        ])
        v6b = pd.DataFrame(columns=["target_uniprot", "theta_mean",
                                     "theta_2p5", "theta_97p5"])
        obs = driver.build_observations(v6a, v6b, driver.PANEL_TARGET_CLASS_MAP)
        assert len(obs) == 2
        # Both observations should have relevance_post_mean = 0.5 (sigmoid(0)
        # because theta defaults to 0 when V6.B has no data)
        for o in obs:
            assert o.relevance_post_mean == pytest.approx(0.5, abs=1e-6)

    def test_build_observations_uses_v6b_theta_when_present(self):
        driver = self._load_driver()
        v6a = pd.DataFrame([
            {"compound_name": "donepezil", "target_uniprot": "P22303",
             "predicted_pkd": 8.5},
        ])
        v6b = pd.DataFrame([
            {"target_uniprot": "P22303", "theta_mean": 1.5,
             "theta_2p5": 1.0, "theta_97p5": 2.0},
        ])
        obs = driver.build_observations(v6a, v6b, driver.PANEL_TARGET_CLASS_MAP)
        assert len(obs) == 1
        # sigmoid(1.5) ≈ 0.818
        assert 0.80 < obs[0].relevance_post_mean < 0.85

    def test_p1_through_p8_parser_handles_dose_suffix(self):
        """The fixed parser must extract 'methylphenidate' from 'P3_methylphenidate_20mg'."""
        from mammal_repurposing.translation.effect_size_model import (
            EffectSizeObservation, fit_effect_size_stub, assert_p1_through_p8,
        )
        obs = [
            EffectSizeObservation("methylphenidate", "NDRI", "Q01959",
                                  pchembl_post_mean=7.0,
                                  relevance_post_mean=0.85),
        ]
        post = fit_effect_size_stub(obs)
        verdicts = assert_p1_through_p8(post)
        # P3_methylphenidate_20mg must resolve to "methylphenidate" → not NO_COMPOUND
        assert verdicts["P3_methylphenidate_20mg"] in ("PASS", "FAIL")
        assert verdicts["P3_methylphenidate_20mg"] != "NO_COMPOUND"

    def test_p1_through_p8_parser_handles_no_dose_suffix(self):
        """For predictions without _Nmg suffix (P1, P6, P7, P8), parser still works."""
        from mammal_repurposing.translation.effect_size_model import (
            EffectSizeObservation, fit_effect_size_stub, assert_p1_through_p8,
        )
        obs = [
            EffectSizeObservation("donepezil", "AChE-I", "P22303",
                                  pchembl_post_mean=8.5,
                                  relevance_post_mean=0.90),
            EffectSizeObservation("intepirdine", "multimodal_5HT", "O43614",
                                  pchembl_post_mean=8.0,
                                  relevance_post_mean=0.50),
        ]
        post = fit_effect_size_stub(obs)
        verdicts = assert_p1_through_p8(post)
        assert verdicts["P1_donepezil"] in ("PASS", "FAIL")
        assert verdicts["P6_intepirdine"] in ("PASS", "FAIL")


# ---------------------------------------------------------------------------
# OSF pre-registration structural validators
# ---------------------------------------------------------------------------
class TestOsfPreRegistrationStructure:
    """OSF pre-reg docs must include required sections per template."""

    REQUIRED_SECTIONS = (
        "Hypothesis",
        "model specification",
        "Sensitivity analyses",
        "Falsifiability",
        "Publication plan",
        "Caveats",
        "OSF + AsPredicted template fields",
    )

    def test_v7_pre_reg_has_all_required_sections(self):
        path = ROOT / "reports" / "paper-drafts" / "v7_osf_preregistration.md"
        assert path.exists(), "V7 OSF pre-reg doc missing"
        body = path.read_text(encoding="utf-8")
        # Case-insensitive matches
        body_lower = body.lower()
        for section in self.REQUIRED_SECTIONS:
            assert section.lower() in body_lower, \
                f"V7 OSF pre-reg missing required section: '{section}'"

    def test_v7_pre_reg_locks_p1_through_p8(self):
        path = ROOT / "reports" / "paper-drafts" / "v7_osf_preregistration.md"
        body = path.read_text(encoding="utf-8").lower()
        # All 8 pre-registered predictions must appear
        for pid in ("p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8"):
            assert pid in body, f"V7 pre-reg missing prediction {pid.upper()}"
        # The 4 validation gates must be enumerated
        assert "gate 1" in body
        assert "gate 2" in body
        assert "gate 3" in body
        assert "gate 4" in body
        # Roberts ceiling must be locked at 0.50
        assert "0.50" in body
        # CPT target venue must be specified
        assert "clinical pharmacology" in body

    def test_v8_pre_reg_has_all_required_sections(self):
        path = ROOT / "reports" / "paper-drafts" / "v8_osf_preregistration.md"
        assert path.exists(), "V8 OSF pre-reg doc missing"
        body = path.read_text(encoding="utf-8")
        body_lower = body.lower()
        for section in self.REQUIRED_SECTIONS:
            assert section.lower() in body_lower, \
                f"V8 OSF pre-reg missing required section: '{section}'"

    def test_v8_pre_reg_locks_clustering_thresholds(self):
        path = ROOT / "reports" / "paper-drafts" / "v8_osf_preregistration.md"
        body = path.read_text(encoding="utf-8").lower()
        # MOFA+ K=30 must be locked
        assert "k = 30" in body or "k=30" in body
        # Leiden γ range
        assert "leiden" in body
        # AMI thresholds
        assert "ami" in body
        assert "0.50" in body or "0.5" in body
        # 9+1 nootropic anchors
        assert "encenicline" in body
        assert "donepezil" in body
        # I_novel novel-mechanism score
        assert "i_novel" in body
        # Nat Mach Intell target venue
        assert ("nature machine intelligence" in body
                or "nat mach intell" in body)

    def test_pre_reg_docs_reference_companion_design_docs(self):
        """Both pre-reg docs must point at design/ companions."""
        for fname in ("v7_osf_preregistration.md", "v8_osf_preregistration.md"):
            body = (ROOT / "reports" / "paper-drafts" / fname).read_text(encoding="utf-8").lower()
            assert "v4_status_and_forward_plan.md" in body or "v4_status" in body
            assert "v6_architecture_plan.md" in body or "v6_architecture" in body

    def test_pre_reg_docs_include_falsifiability_action(self):
        """Each pre-reg must explicitly state what happens if gates fail."""
        v7_body = (ROOT / "reports" / "paper-drafts" / "v7_osf_preregistration.md").read_text(encoding="utf-8").lower()
        # V7 must mention CPT:PSP negative-result fallback
        assert "negative-result" in v7_body or "negative result" in v7_body
        v8_body = (ROOT / "reports" / "paper-drafts" / "v8_osf_preregistration.md").read_text(encoding="utf-8").lower()
        # V8 must mention degrade band or negative result
        assert ("negative" in v8_body or "degrade" in v8_body
                or "fallback" in v8_body)


# ---------------------------------------------------------------------------
# Methodology v3 + README structural sanity
# ---------------------------------------------------------------------------
class TestMethodologyV3:
    def test_methodology_v3_exists_and_covers_all_5_layers(self):
        path = ROOT / "reports" / "paper-drafts" / "methodology_v3.md"
        assert path.exists(), "methodology_v3.md missing"
        body = path.read_text(encoding="utf-8").lower()
        for layer in ("v4", "v5", "v6.a", "v6.b", "v7", "v8"):
            assert layer in body, f"methodology_v3 missing layer: {layer}"

    def test_methodology_v3_includes_falsifier_table(self):
        path = ROOT / "reports" / "paper-drafts" / "methodology_v3.md"
        body = path.read_text(encoding="utf-8").lower()
        assert "falsif" in body
        # Roberts ceiling reference
        assert "roberts" in body
        # Hypothesis audit count
        assert "hypothesis audit" in body or "hypothesis-audit" in body

    def test_readme_reflects_v4_through_v8(self):
        path = ROOT / "README.md"
        body = path.read_text(encoding="utf-8").lower()
        for layer in ("v4", "v5", "v6.a", "v6.b", "v7", "v8"):
            assert layer in body, f"README missing layer: {layer}"
        # Production NUTS metric
        assert "r̂" in body or "rhat" in body or "r-hat" in body or "r=1.000" in body or "1.000" in body
        # Wet-lab shortlist v10
        assert "v10" in body
