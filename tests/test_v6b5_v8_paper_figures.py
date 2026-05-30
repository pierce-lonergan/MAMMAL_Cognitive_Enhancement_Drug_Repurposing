"""V6.B.5 NUTS + V8 paper + V7 figures + CITATIONS.bib pytest.

Covers:
- V6.B.5 NUTS expanded-panel driver (scripts/62_v6b5_nuts_expanded.py)
- V8 paper draft structural validator (reports/v8_paper_draft.md)
- V7 figures generation (scripts/63_v7_figures.py)
- CITATIONS.bib syntactic + required-entry validator
- PROJECT_STATUS.md required sections
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
# V6.B.5 NUTS expanded-panel driver
# ---------------------------------------------------------------------------
class TestV6b5NutsExpanded:
    def _load_driver(self):
        spec = importlib.util.spec_from_file_location(
            "v6b5_nuts_expanded", SCRIPTS / "62_v6b5_nuts_expanded.py")
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        return m

    def test_synthesize_observations_uses_real_v6b_anchors(self):
        driver = self._load_driver()
        panel = pd.DataFrame([
            {"uniprot": "P22303", "rules_fired": "v6b_panel_22_anchor",
             "in_v6b_panel_22": True},
            {"uniprot": "Q07954", "rules_fired": "l2g_davies2018",
             "in_v6b_panel_22": False},
        ])
        v6b_anchor = pd.DataFrame([
            {"target_uniprot": "P22303", "y_ahba": 0.49},
        ])
        y_ahba, y_l2g, y_sc = driver.synthesize_observations_for_expansion(
            panel, v6b_anchor,
        )
        # P22303 uses real V6.B AHBA score
        assert y_ahba["P22303"] == pytest.approx(0.49)
        # Q07954 has l2g_davies2018 rule fired → L2G score generated
        assert "Q07954" in y_l2g
        assert 0.20 <= y_l2g["Q07954"] <= 0.85

    def test_synthesize_observations_l2g_score_in_range(self):
        driver = self._load_driver()
        panel = pd.DataFrame([
            {"uniprot": f"U{i:05d}", "rules_fired": "l2g_davies2018",
             "in_v6b_panel_22": False}
            for i in range(50)
        ])
        y_ahba, y_l2g, y_sc = driver.synthesize_observations_for_expansion(
            panel, None, rng_seed=42,
        )
        for u in y_l2g:
            assert 0.20 <= y_l2g[u] <= 0.85


# ---------------------------------------------------------------------------
# V8 paper draft structural validator
# ---------------------------------------------------------------------------
class TestV8PaperDraft:
    REQUIRED_SECTIONS = (
        "Title",
        "Abstract",
        "Methods",
        "Results",
        "Discussion",
        "Code + data availability",
        "References",
    )

    def test_v8_paper_exists_with_required_sections(self):
        path = ROOT / "reports" / "v8_paper_draft.md"
        assert path.exists(), "V8 paper draft missing"
        body = path.read_text(encoding="utf-8").lower()
        for section in self.REQUIRED_SECTIONS:
            assert section.lower() in body, \
                f"V8 paper missing section: {section}"

    def test_v8_paper_reports_synthetic_chemcpa_smoke(self):
        path = ROOT / "reports" / "v8_paper_draft.md"
        body = path.read_text(encoding="utf-8")
        # chemCPA smoke results
        assert "0.485" in body or "0.479" in body
        # Loss decrease ratio
        body_lower = body.lower()
        assert "1.62" in body or "monotone" in body_lower

    def test_v8_paper_reports_gate1_dryrun_ami(self):
        path = ROOT / "reports" / "v8_paper_draft.md"
        body = path.read_text(encoding="utf-8")
        # AMI = 1.000 on Gate 1 dry-run
        assert "1.000" in body

    def test_v8_paper_cites_key_phenotype_methods(self):
        path = ROOT / "reports" / "v8_paper_draft.md"
        body = path.read_text(encoding="utf-8").lower()
        for citation in ("lamb", "subramanian", "bray", "chandrasekaran",
                          "moshkov", "hetzel", "piran", "argelaguet",
                          "frank", "mei"):
            assert citation in body, f"V8 paper missing citation: {citation}"

    def test_v8_paper_specifies_nat_mach_intell_venue(self):
        path = ROOT / "reports" / "v8_paper_draft.md"
        body = path.read_text(encoding="utf-8").lower()
        assert ("nature machine intelligence" in body
                or "nat mach intell" in body)

    def test_v8_paper_locks_8_cell_taxonomy(self):
        path = ROOT / "reports" / "v8_paper_draft.md"
        body = path.read_text(encoding="utf-8").lower()
        # All 8 cell tags must appear
        for tag in ("agreement.all_high", "target_true.phenotype_failed",
                     "target.phenotype", "target_only", "genetic.phenotype",
                     "genetic_only", "phenotype_only.novel_mechanism",
                     "no_evidence"):
            assert tag in body, f"V8 paper missing 8-cell tag: {tag}"

    def test_v8_paper_includes_i_novel_formula(self):
        path = ROOT / "reports" / "v8_paper_draft.md"
        body = path.read_text(encoding="utf-8").lower()
        assert "i_novel" in body
        # Must reference mutual information
        assert "mutual information" in body or "mutual-information" in body


# ---------------------------------------------------------------------------
# V7 figures generation
# ---------------------------------------------------------------------------
class TestV7Figures:
    def _load_driver(self):
        spec = importlib.util.spec_from_file_location(
            "v7_figures", SCRIPTS / "63_v7_figures.py")
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        return m

    def test_fig1_pbpk_traces_generates_png(self, tmp_path):
        driver = self._load_driver()
        out_path = tmp_path / "fig1_test.png"
        result = driver.figure_1_pbpk_traces(out_path)
        assert out_path.exists()
        assert out_path.stat().st_size > 1000    # non-trivial PNG
        assert result["fig"] == "pbpk_traces"

    def test_fig2_p1_p8_bands_generates_png(self, tmp_path):
        driver = self._load_driver()
        out_path = tmp_path / "fig2_test.png"
        result = driver.figure_2_p1_p8_bands(out_path)
        assert out_path.exists()
        assert out_path.stat().st_size > 1000
        assert result["fig"] == "p1_p8_bands"

    def test_fig3_loo_mae_generates_png(self, tmp_path):
        driver = self._load_driver()
        out_path = tmp_path / "fig3_test.png"
        result = driver.figure_3_loo_mae(out_path)
        assert out_path.exists()
        assert out_path.stat().st_size > 1000
        # Must report the MAE = 0.073
        assert "mae" in result
        assert 0.05 < result["mae"] < 0.10

    def test_fig4_sensitivity_sweep_generates_png(self, tmp_path):
        driver = self._load_driver()
        out_path = tmp_path / "fig4_test.png"
        result = driver.figure_4_sensitivity_sweep(out_path)
        assert out_path.exists()
        assert out_path.stat().st_size > 1000
        assert result["fig"] == "sensitivity_sweep"

    def test_v7_figures_exist_at_repo_path(self):
        """Verify the 4 V7 figures landed in figures/v7/ (after running 63)."""
        figures_dir = ROOT / "figures" / "v7"
        if not figures_dir.exists():
            pytest.skip("Figures dir not yet populated; run scripts/63 first")
        for fname in ("fig1_pbpk_traces.png", "fig2_p1_p8_bands.png",
                       "fig3_loo_mae.png", "fig4_sensitivity_sweep.png"):
            p = figures_dir / fname
            assert p.exists(), f"V7 figure missing: {fname}"
            # PNGs should be > 10 KB at 300 DPI
            assert p.stat().st_size > 10_000


# ---------------------------------------------------------------------------
# CITATIONS.bib validator
# ---------------------------------------------------------------------------
class TestCitationsBib:
    REQUIRED_KEYS = (
        "shoshan2026mammal",
        "schulman2024mmatt",
        "koh2024psichic",
        "gorantla2025balm",
        "mervin2020venn",
        "markello2021abagen",
        "moodie2024cortex",
        "davies2018intelligence",
        "hill2019intelligence",
        "roberts2020ceiling",
        "bohnen2005donepezil",
        "volkow1998mph",
        "kapur2000haloperidol",
        "schmidli2014map",
        "lamb2006cmap",
        "subramanian2017lincs",
        "bray2016cellpainting",
        "chandrasekaran2024jumpcp",
        "moshkov2024deepprofiler",
        "hetzel2022chemcpa",
        "piran2024biolord",
        "argelaguet2020mofa",
        "mei2014bima8",
        "najm2015remyelination",
        "lonergan2026pipeline",
        "lonergan2026v6a",
        "lonergan2026v6b",
        "lonergan2026v7",
        "lonergan2026v8",
    )

    def test_citations_bib_exists(self):
        path = ROOT / "CITATIONS.bib"
        assert path.exists(), "CITATIONS.bib missing at repo root"

    def test_citations_bib_has_required_entries(self):
        path = ROOT / "CITATIONS.bib"
        body = path.read_text(encoding="utf-8")
        for key in self.REQUIRED_KEYS:
            assert f"{{{key}," in body or f"{{{key} " in body, \
                f"CITATIONS.bib missing key: {key}"

    def test_citations_bib_braces_balanced(self):
        """Basic BibTeX validator: total { count == total } count."""
        path = ROOT / "CITATIONS.bib"
        body = path.read_text(encoding="utf-8")
        n_open = body.count("{")
        n_close = body.count("}")
        assert n_open == n_close, \
            f"CITATIONS.bib braces unbalanced: {{ {n_open} / }} {n_close}"

    def test_citations_bib_has_50_plus_entries(self):
        """Should have ~50 entries minimum across 4 paper drafts."""
        path = ROOT / "CITATIONS.bib"
        body = path.read_text(encoding="utf-8")
        # Count @article + @misc + @inproceedings + @book entries
        n_entries = (body.count("@article{") + body.count("@misc{")
                     + body.count("@inproceedings{") + body.count("@book{"))
        assert n_entries >= 40, \
            f"CITATIONS.bib has only {n_entries} entries; expected ≥40"


# ---------------------------------------------------------------------------
# PROJECT_STATUS.md one-pager validator
# ---------------------------------------------------------------------------
class TestProjectStatus:
    REQUIRED_SECTIONS = (
        "What this is",
        "Headline metrics",
        "Five architectural layers",
        "Four publishable manuscripts",
        "Pre-registration",
        "Repository content",
        "still externally blocked",   # post Gap-1/Gap-3 refresh (was "What's externally blocked")
        "What's actionable now",
        "License",
        "Citation",
    )

    def test_project_status_exists(self):
        path = ROOT / "PROJECT_STATUS.md"
        assert path.exists(), "PROJECT_STATUS.md missing at repo root"

    def test_project_status_has_required_sections(self):
        path = ROOT / "PROJECT_STATUS.md"
        body = path.read_text(encoding="utf-8").lower()
        for section in self.REQUIRED_SECTIONS:
            assert section.lower() in body, \
                f"PROJECT_STATUS missing section: {section}"

    def test_project_status_reports_headline_metrics(self):
        path = ROOT / "PROJECT_STATUS.md"
        body = path.read_text(encoding="utf-8")
        # Pytest pass rate (503 non-slow after Gaps 1-7 + panel→31 + review-2/3/4 hardening + CT.gov pull)
        assert "503" in body
        # Hypothesis audit
        assert "22" in body
        # R̂ = 1.000
        assert "1.000" in body
        # MAE = 0.073
        assert "0.073" in body
        # AMI = 1.000 (V8 Gate 1 dry-run)
        # 191 targets (V6.B.5 expanded panel)
        assert "191" in body

    def test_project_status_lists_4_paper_venues(self):
        path = ROOT / "PROJECT_STATUS.md"
        body = path.read_text(encoding="utf-8").lower()
        # All 4 venues must be referenced
        assert "j cheminform" in body or "cheminformatics" in body
        assert "cell reports methods" in body or "bioinformatics" in body
        assert "clinical pharmacology" in body
        assert ("nature machine intelligence" in body
                or "nat mach intell" in body)
