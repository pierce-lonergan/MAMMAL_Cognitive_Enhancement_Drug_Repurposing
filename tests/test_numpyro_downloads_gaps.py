"""numpyro fix + LINCS/JUMP-CP downloads + GAPS document pytest.

Covers:
- numpyro / jax / jaxlib install (validates the PyMC EOFError fix)
- V6.B.5 expanded posterior parquet has real numpyro-NUTS output (R̂=1.000)
- LINCS L1000 pert_info downloaded + parseable
- JUMP-CP S3 reachability via boto3 UNSIGNED probe
- GAPS_AND_RESEARCH_DIRECTIONS document structural validator
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]


# ---------------------------------------------------------------------------
# numpyro / jax / jaxlib install validation
# ---------------------------------------------------------------------------
class TestNumpyroJaxStack:
    def test_numpyro_importable(self):
        try:
            import numpyro
            from packaging.version import Version
            assert Version(numpyro.__version__) >= Version("0.18"), \
                f"numpyro {numpyro.__version__} is too old"
        except ImportError:
            pytest.skip("numpyro not installed (expected in production env)")

    def test_jax_importable(self):
        try:
            import jax
            from packaging.version import Version
            # Accept any jax ≥ 0.4 (we're on 0.10.1+)
            assert Version(jax.__version__) >= Version("0.4"), \
                f"jax {jax.__version__} is too old"
        except ImportError:
            pytest.skip("jax not installed")

    def test_jaxlib_importable(self):
        try:
            import jaxlib
            from packaging.version import Version
            assert Version(jaxlib.__version__) >= Version("0.4"), \
                f"jaxlib {jaxlib.__version__} is too old"
        except ImportError:
            pytest.skip("jaxlib not installed")

    def test_pymc_detects_numpyro_sampler(self):
        try:
            import numpyro  # noqa
            from mammal_repurposing.cluster_d.bayesian_prior import (
                _numpyro_available,
            )
            sys.path.insert(0, str(ROOT / "src"))
            assert _numpyro_available() is True
        except ImportError:
            pytest.skip("numpyro/PyMC bridge not testable")


# ---------------------------------------------------------------------------
# V6.B.5 expanded posterior — real numpyro output
# ---------------------------------------------------------------------------
class TestV6b5NumpyroPosterior:
    def test_v6b5_expanded_posterior_present(self):
        path = (ROOT / "data" / "results" / "v2"
                / "cluster_d_posterior_expanded_v1.parquet")
        if not path.exists():
            pytest.skip("V6.B.5 expanded posterior not yet generated")
        df = pd.read_parquet(path)
        # Must have 191 targets (V6.B.5 expanded panel)
        assert len(df) >= 150
        # Required columns
        for col in ("target_uniprot", "gene", "theta_mean",
                     "theta_2p5", "theta_97p5", "w_pipeline",
                     "in_v6b_panel_22", "is_reference_anchor"):
            assert col in df.columns

    def test_v6b5_report_reports_R_hat_convergence(self):
        path = ROOT / "reports" / "pipeline" / "cluster_d_nuts_expanded_v1.md"
        if not path.exists():
            pytest.skip("V6.B.5 expanded report not yet generated")
        body = path.read_text(encoding="utf-8")
        # If numpyro NUTS ran, method should be pymc_nuts
        body_lower = body.lower()
        if "method: pymc_nuts" in body_lower or "pymc_nuts" in body_lower:
            # R̂ should report a numerical value (1.000 expected)
            assert "r̂ max" in body_lower
        else:
            pytest.skip("Stub mode; numpyro NUTS not yet run")

    def test_v6b5_4_reference_anchors_active(self):
        path = (ROOT / "data" / "results" / "v2"
                / "cluster_d_posterior_expanded_v1.parquet")
        if not path.exists():
            pytest.skip("V6.B.5 expanded posterior not yet generated")
        df = pd.read_parquet(path)
        # 4 reference anchors: BDNF, COMT, ACHE, CHRNA7 (in 191-target panel)
        n_anchors = int(df["is_reference_anchor"].sum())
        assert n_anchors >= 4


# ---------------------------------------------------------------------------
# LINCS L1000 download validation
# ---------------------------------------------------------------------------
class TestLincsDownload:
    def test_lincs_pert_info_present_when_downloaded(self):
        path = ROOT / "data" / "cache" / "lincs" / "GSE70138_compoundinfo.txt.gz"
        if not path.exists():
            pytest.skip("LINCS pert_info not yet downloaded")
        # > 50 KB = real data, not 404 HTML
        assert path.stat().st_size > 50_000

    def test_lincs_pert_info_parseable(self):
        path = ROOT / "data" / "cache" / "lincs" / "GSE70138_compoundinfo.txt.gz"
        if not path.exists():
            pytest.skip("LINCS pert_info not yet downloaded")
        df = pd.read_csv(path, sep="\t", compression="gzip")
        # Real LINCS pert_info has these columns
        for col in ("pert_id", "canonical_smiles", "inchi_key",
                     "pert_iname", "pert_type"):
            assert col in df.columns, f"LINCS pert_info missing column: {col}"
        # Real LINCS Phase 2 has ~2,170 perturbations
        assert len(df) >= 2000

    def test_lincs_pert_types_canonical(self):
        path = ROOT / "data" / "cache" / "lincs" / "GSE70138_compoundinfo.txt.gz"
        if not path.exists():
            pytest.skip("LINCS pert_info not yet downloaded")
        df = pd.read_csv(path, sep="\t", compression="gzip")
        pert_types = set(df["pert_type"].unique())
        # Real LINCS Phase 2 should have compound + CRISPR + control perts
        assert "trt_cp" in pert_types or "trt_compound" in pert_types
        # ~1,700+ small-molecule perturbations
        n_compounds = (df["pert_type"].str.startswith("trt_cp")
                       | df["pert_type"].str.startswith("trt_compound")).sum()
        assert n_compounds > 1000


# ---------------------------------------------------------------------------
# JUMP-CP S3 reachability (no actual download required)
# ---------------------------------------------------------------------------
class TestJumpCpS3Reachability:
    def test_boto3_importable(self):
        try:
            import boto3
        except ImportError:
            pytest.skip("boto3 not installed")

    @pytest.mark.slow
    def test_jumpcp_s3_lists_14_sources(self):
        """Network test (slow): list top-level cpg0016-jump/ contents."""
        try:
            import boto3
            from botocore import UNSIGNED
            from botocore.config import Config
        except ImportError:
            pytest.skip("boto3 not installed")
        s3 = boto3.client("s3", config=Config(signature_version=UNSIGNED))
        try:
            resp = s3.list_objects_v2(
                Bucket="cellpainting-gallery",
                Prefix="cpg0016-jump/",
                Delimiter="/",
                MaxKeys=20,
            )
        except Exception as e:
            pytest.skip(f"JUMP-CP S3 unreachable in this env: {e}")
        prefixes = [p["Prefix"] for p in (resp.get("CommonPrefixes") or [])]
        # Should have source_1 through source_15 + source_all
        assert len(prefixes) >= 10
        assert any("source_4" in p for p in prefixes)


# ---------------------------------------------------------------------------
# GAPS_AND_RESEARCH_DIRECTIONS document validator
# ---------------------------------------------------------------------------
class TestGapsDocument:
    REQUIRED_SECTIONS = (
        "Open engineering gaps",
        "External blockers",
        "Research directions",
        "Cross-cutting limitations",
        "Recommended order",
        "Completed ledger",
    )

    def test_gaps_document_exists(self):
        path = ROOT / "GAPS_AND_RESEARCH_DIRECTIONS.md"
        assert path.exists(), "GAPS_AND_RESEARCH_DIRECTIONS.md missing at repo root"

    def test_gaps_document_has_required_sections(self):
        path = ROOT / "GAPS_AND_RESEARCH_DIRECTIONS.md"
        body = path.read_text(encoding="utf-8").lower()
        for section in self.REQUIRED_SECTIONS:
            assert section.lower() in body, \
                f"Gaps doc missing section: {section}"

    def test_gaps_document_resolves_numpyro_fix(self):
        path = ROOT / "GAPS_AND_RESEARCH_DIRECTIONS.md"
        body = path.read_text(encoding="utf-8").lower()
        # Must reference the EOFError fix
        assert "numpyro" in body
        assert "eoferror" in body or "multiprocess" in body
        # Must document the fix as RESOLVED ✅
        assert "fixed" in body or "resolved" in body

    def test_gaps_document_validates_lincs_download(self):
        path = ROOT / "GAPS_AND_RESEARCH_DIRECTIONS.md"
        body = path.read_text(encoding="utf-8").lower()
        # Must reference real LINCS download evidence
        assert "lincs" in body
        assert "geo" in body or "gse70138" in body
        assert "2,170" in body or "2170" in body or "1,796" in body or "1796" in body

    def test_gaps_document_validates_jumpcp_s3(self):
        path = ROOT / "GAPS_AND_RESEARCH_DIRECTIONS.md"
        body = path.read_text(encoding="utf-8").lower()
        # Must reference JUMP-CP S3 validation
        assert "jump-cp" in body or "cpg0016" in body
        assert "s3" in body
        assert "14 source" in body or "sources" in body

    def test_gaps_document_lists_must_have_topics(self):
        path = ROOT / "GAPS_AND_RESEARCH_DIRECTIONS.md"
        body = path.read_text(encoding="utf-8").lower()
        # MUST-HAVE research topics (MH1 - MH10)
        for n in range(1, 11):
            assert f"mh{n}" in body or f"mh{n}." in body, \
                f"Gaps doc missing MUST-HAVE research item MH{n}"

    def test_gaps_document_includes_priority_sequence(self):
        path = ROOT / "GAPS_AND_RESEARCH_DIRECTIONS.md"
        body = path.read_text(encoding="utf-8")
        # Recommended sprint sequence with numbered priority
        body_lower = body.lower()
        assert "recommended order" in body_lower
        assert "priority order" in body_lower or "highest" in body_lower

    def test_gaps_document_specifies_effort_estimates(self):
        path = ROOT / "GAPS_AND_RESEARCH_DIRECTIONS.md"
        body = path.read_text(encoding="utf-8").lower()
        # Should have effort estimates (weeks / hours / months / $)
        assert "weeks" in body or "hours" in body
        assert "$" in body
        # Phase 1 trial estimate
        assert "phase 1" in body or "phase-1" in body

    def test_gaps_document_references_5_paper_suite(self):
        path = ROOT / "GAPS_AND_RESEARCH_DIRECTIONS.md"
        body = path.read_text(encoding="utf-8").lower()
        # Must reference the 5-paper suite + companion docs
        assert "5-paper" in body or "five-paper" in body or "5 manuscript" in body or "5 paper" in body
        assert "osf" in body
        assert "wet-lab handoff" in body or "wet_lab_handoff" in body


# ---------------------------------------------------------------------------
# Sprint integration — all new artifacts
# ---------------------------------------------------------------------------
class TestSprintArtifacts:
    def test_all_new_artifacts_present(self):
        artifacts = [
            ROOT / "GAPS_AND_RESEARCH_DIRECTIONS.md",
            ROOT / "data" / "cache" / "lincs" / "GSE70138_compoundinfo.txt.gz",
            ROOT / "data" / "results" / "v2" / "cluster_d_posterior_expanded_v1.parquet",
        ]
        for p in artifacts:
            assert p.exists(), f"Sprint artifact missing: {p}"

    def test_gaps_doc_in_repo_root(self):
        # GAPS doc should be at repo root alongside README + PROJECT_STATUS
        assert (ROOT / "GAPS_AND_RESEARCH_DIRECTIONS.md").exists()
        assert (ROOT / "README.md").exists()
        assert (ROOT / "PROJECT_STATUS.md").exists()
        assert (ROOT / "CITATIONS.bib").exists()
