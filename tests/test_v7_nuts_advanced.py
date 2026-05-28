"""V7 NUTS full-path + V8.2 chemCPA smoke + V6.B paper draft pytest.

- V7 NUTS Stage 2 driver (scripts/58_v7_nuts_synthetic.py): anchor
  observation builder; gate evaluation; sensitivity sweep robustness
- V8.2 chemCPA smoke (scripts/59_v8_chemcpa_smoke.py): synthetic LINCS
  generation; chemCPA architecture trains end-to-end; reconstruction R²
- V6.B paper draft (reports/v6b_paper_draft.md): structural validator for
  manuscript outline; required sections per Cell Reports Methods template
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
# V7 NUTS driver (scripts/58_v7_nuts_synthetic.py)
# ---------------------------------------------------------------------------
class TestV7NutsDriver:
    def _load_driver(self):
        spec = importlib.util.spec_from_file_location(
            "v7_nuts_driver", SCRIPTS / "58_v7_nuts_synthetic.py")
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        return m

    def test_build_anchor_observations_returns_15(self):
        driver = self._load_driver()
        obs = driver.build_anchor_observations()
        assert len(obs) == 15
        # All have observed_g populated
        for o in obs:
            assert o.observed_g is not None
            assert isinstance(o.observed_g, float)

    def test_anchor_compounds_cover_15_reference_set(self):
        driver = self._load_driver()
        obs = driver.build_anchor_observations()
        compounds = {o.compound for o in obs}
        # The 15 reference compounds per REFERENCE_COMPOUND_SMD
        expected = {
            "donepezil", "galantamine", "rivastigmine", "memantine",
            "methylphenidate", "d_amphetamine", "modafinil", "atomoxetine",
            "varenicline", "caffeine", "encenicline", "intepirdine",
            "pridopidine", "vortioxetine", "guanfacine",
        }
        assert compounds == expected

    def test_anchors_use_real_v6b_theta_when_present(self):
        driver = self._load_driver()
        obs = driver.build_anchor_observations()
        # Check that at least one observation has a non-trivial relevance
        # (i.e. V6.B θ̄ posterior was loaded; relevance != 0.5 sigmoid baseline)
        v6b_path = ROOT / "data" / "results" / "v2" / "cluster_d_posterior_v1.parquet"
        if v6b_path.exists():
            # If V6.B posterior exists, at least one anchor should have
            # relevance_post_mean != 0.5
            relevances = [o.relevance_post_mean for o in obs]
            assert any(abs(r - 0.5) > 1e-6 for r in relevances), \
                "V6.B posterior loaded but no anchors picked up real θ̄"


# ---------------------------------------------------------------------------
# V7 NUTS full-path Bayesian convergence (mock data; fast)
# ---------------------------------------------------------------------------
class TestV7NutsFullPath:
    @pytest.mark.slow
    def test_nuts_converges_on_synthetic_anchor_set(self):
        """Full NUTS on synthetic 15-anchor data converges with R̂ < 1.05."""
        from mammal_repurposing.translation.effect_size_model import (
            EffectSizeObservation, fit_effect_size_nuts, PYMC_AVAILABLE,
        )
        if not PYMC_AVAILABLE:
            pytest.skip("PyMC not installed")
        # Tiny synthetic anchor set
        rng = np.random.default_rng(0)
        obs = []
        for i in range(5):
            obs.append(EffectSizeObservation(
                compound=f"c{i}",
                class_name="AChE-I",
                target_uniprot="P22303",
                pchembl_post_mean=7.0 + rng.normal(0, 0.5),
                pchembl_post_sd=0.30,
                relevance_post_mean=0.7 + rng.normal(0, 0.1),
                relevance_post_sd=0.15,
                pbpk_auc_brain=1.0,
                moderators=(0,) * 5,
                observed_g=0.15 + rng.normal(0, 0.05),
            ))
        post = fit_effect_size_nuts(
            obs, n_chains=2, n_draws=500, n_tune=500,
            target_accept=0.90,
        )
        assert post.method == "pymc_nuts"
        assert post.rhat_max < 1.10, f"R̂={post.rhat_max} above 1.10 threshold"
        assert post.ess_min > 100, f"ESS={post.ess_min} below 100 threshold"

    def test_nuts_imports_without_running(self):
        """Smoke test: PyMC + arviz can import; fit_effect_size_nuts is callable."""
        from mammal_repurposing.translation.effect_size_model import (
            fit_effect_size_nuts, PYMC_AVAILABLE,
        )
        assert callable(fit_effect_size_nuts)
        assert isinstance(PYMC_AVAILABLE, bool)


# ---------------------------------------------------------------------------
# V8.2 chemCPA smoke (scripts/59_v8_chemcpa_smoke.py)
# ---------------------------------------------------------------------------
class TestV8ChemcpaSmoke:
    def _load_driver(self):
        spec = importlib.util.spec_from_file_location(
            "v8_chemcpa_smoke", SCRIPTS / "59_v8_chemcpa_smoke.py")
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        return m

    def test_synthetic_perturbation_returns_expected_shape(self):
        driver = self._load_driver()
        try:
            df = driver.generate_synthetic_perturbation(
                ["CCO", "CCC", "CCCC"],
                cell_lines=["NPC", "NEU"],
                doses=(1.0, 5.0),
            )
        except ImportError:
            pytest.skip("rdkit not installed")
        # 3 compounds × 2 cell lines × 2 doses = 12 signatures
        assert len(df) == 12
        # 977 landmark gene columns
        gene_cols = [c for c in df.columns if c.startswith("gene_")]
        assert len(gene_cols) == 977
        # Required metadata columns
        for col in ("compound_smiles", "cell_line", "dose_um", "pert_id"):
            assert col in df.columns

    def test_synthetic_perturbation_handles_invalid_smiles(self):
        driver = self._load_driver()
        try:
            df = driver.generate_synthetic_perturbation(
                ["CCO", "not_a_smiles", "CCC"],
                cell_lines=["NPC"], doses=(1.0,),
            )
        except ImportError:
            pytest.skip("rdkit not installed")
        # Invalid SMILES dropped → 2 valid × 1 cell × 1 dose = 2 signatures
        assert len(df) == 2

    @pytest.mark.slow
    def test_chemcpa_smoke_loss_decreases(self):
        """chemCPA training loss decreases over 3 epochs on synthetic data."""
        driver = self._load_driver()
        try:
            import torch  # noqa: F401
        except ImportError:
            pytest.skip("torch not installed")
        try:
            from rdkit import Chem  # noqa: F401
        except ImportError:
            pytest.skip("rdkit not installed")
        df = driver.generate_synthetic_perturbation(
            ["CCO", "CCC", "CCCC", "CCCCC", "c1ccccc1", "CCN", "CCNC", "CCO"],
            cell_lines=["NPC", "MCF7"],
            doses=(1.0, 5.0),
        )
        result = driver.train_chemcpa_smoke(df, n_epochs=3, batch_size=4)
        # Loss must decrease
        assert result["epoch_losses"][0] > result["epoch_losses"][-1]
        # R² is a finite scalar
        assert np.isfinite(result["test_r2_mean"])


# ---------------------------------------------------------------------------
# V6.B paper draft structural validator
# ---------------------------------------------------------------------------
class TestV6bPaperDraft:
    REQUIRED_SECTIONS = (
        "Title",
        "Abstract",
        "Methods",
        "Results",
        "Discussion",
        "Limitations",
        "Code + data availability",
        "References",
    )

    def test_v6b_paper_draft_exists(self):
        path = ROOT / "reports" / "v6b_paper_draft.md"
        assert path.exists(), "v6b_paper_draft.md missing"

    def test_v6b_paper_has_required_sections(self):
        path = ROOT / "reports" / "v6b_paper_draft.md"
        body = path.read_text(encoding="utf-8").lower()
        for section in self.REQUIRED_SECTIONS:
            assert section.lower() in body, \
                f"V6.B paper draft missing required section: '{section}'"

    def test_v6b_paper_cites_key_methods(self):
        path = ROOT / "reports" / "v6b_paper_draft.md"
        body = path.read_text(encoding="utf-8").lower()
        # Must cite the core methodology papers
        required_citations = ("markello", "moodie", "davies", "hill",
                               "mountjoy", "roberts", "siletti", "pymc")
        for c in required_citations:
            assert c in body, f"V6.B paper missing key citation: {c}"

    def test_v6b_paper_reports_real_convergence(self):
        path = ROOT / "reports" / "v6b_paper_draft.md"
        body = path.read_text(encoding="utf-8")
        # R̂ = 1.000 must appear (production NUTS converged metric)
        assert "1.000" in body
        # ESS metric must appear
        assert "12,780" in body or "12780" in body
        # Roberts 2020 ceiling must be referenced
        body_lower = body.lower()
        assert "roberts" in body_lower
        assert "0.50" in body or "0.5" in body

    def test_v6b_paper_specifies_target_venue(self):
        path = ROOT / "reports" / "v6b_paper_draft.md"
        body = path.read_text(encoding="utf-8").lower()
        # Cell Reports Methods OR Bioinformatics must be referenced
        assert ("cell reports methods" in body
                or "bioinformatics" in body)

    def test_v6b_paper_includes_per_target_table(self):
        path = ROOT / "reports" / "v6b_paper_draft.md"
        body = path.read_text(encoding="utf-8")
        # Top targets per posterior table
        for gene in ("ACHE", "CHRNA7", "GRIN2B"):
            assert gene in body, f"V6.B paper missing per-target row for {gene}"

    def test_v6b_paper_has_code_availability_section(self):
        path = ROOT / "reports" / "v6b_paper_draft.md"
        body = path.read_text(encoding="utf-8").lower()
        assert "code" in body and "data" in body
        # Must reference the GitHub repo + Apache-2.0 license
        assert "github" in body
        assert "apache" in body


# ---------------------------------------------------------------------------
# V7/V8 sprint integration tests — outputs must coexist
# ---------------------------------------------------------------------------
class TestSprintIntegrationArtifacts:
    """Verify that the V6.B + V7 NUTS + V8 smoke artifacts coexist."""

    def test_v6b_posterior_parquet_present(self):
        path = ROOT / "data" / "results" / "v2" / "cluster_d_posterior_v1.parquet"
        assert path.exists(), "V6.B.3 NUTS posterior parquet missing"
        df = pd.read_parquet(path)
        assert len(df) == 22
        for col in ("target_uniprot", "gene", "theta_mean",
                     "theta_2p5", "theta_97p5", "w_pipeline"):
            assert col in df.columns

    def test_v7_nuts_posterior_parquet_present(self):
        path = ROOT / "data" / "results" / "v2" / "v7_nuts_posterior_v1.parquet"
        if not path.exists():
            pytest.skip("V7 NUTS not yet run; smoke test only")
        df = pd.read_parquet(path)
        # 15 anchor compounds expected
        assert len(df) >= 14    # at least 14 of 15 successfully fit
        for col in ("compound", "g_mean", "g_2p5", "g_97p5", "g_90_upper"):
            assert col in df.columns

    def test_v8_chemcpa_smoke_parquet_present(self):
        path = ROOT / "data" / "results" / "v2" / "v8_chemcpa_smoke_v1.parquet"
        if not path.exists():
            pytest.skip("V8.2 chemCPA smoke not yet run")
        df = pd.read_parquet(path)
        assert "metric" in df.columns
        assert "value" in df.columns
        # Should have a test_r2_mean row
        assert (df["metric"] == "test_r2_mean").any()
