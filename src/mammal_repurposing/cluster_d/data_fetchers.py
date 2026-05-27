"""V6 §13.2 — Cluster D data-fetch scaffolds for AHBA / OT Genetics / cellxgene.

Each fetcher is code-complete with graceful degradation. When the
corresponding heavy dependency lands (abagen / tiledbsoma / brainsmash), the
real-path branch activates.

The static §8.6 brain_region.py annotations remain as the V5 fallback when
none of these are installed.

References per fetcher:
  - AHBA via abagen — Markello et al. 2021 eLife 10:e72129
  - OT Genetics L2G GraphQL — Mountjoy 2021 Nat Genet 53:1527 (V4 endpoint at
    https://api.genetics.opentargets.org/graphql); Suzuki et al. 2024 Cell
    Genomics is the newer L2G release
  - cellxgene-census tiledbsoma — Siletti 2023 Science add7046; CZ CELLxGENE
    Discover Census 2025-11-08 LTS
  - BrainSMASH — Burt et al. 2020 NeuroImage 220:117038
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)


@dataclass
class AHBAExpressionData:
    target_uniprot: str
    gene_symbol: str
    region_expression: dict[str, float]    # DK68 parcel → expression z-score
    cortical_r: float = float("nan")       # corr with cognition reference map
    brainsmash_p: float = float("nan")     # spatial null
    stability_r: float = float("nan")      # AHBA donor stability (Hawrylycz 2015)
    note: str = "stub"


@dataclass
class OTGeneticsL2G:
    target_ensembl: str
    target_uniprot: str
    gene_symbol: str
    l2g_max_score: float = float("nan")
    contributing_studies: list[str] = field(default_factory=list)
    note: str = "stub"


@dataclass
class CellTypeEnrichment:
    target_uniprot: str
    gene_symbol: str
    enrichment_by_celltype: dict[str, float] = field(default_factory=dict)
    cognition_salient_z: float = float("nan")  # weighted-sum over §A.8 cognition C* set
    note: str = "stub"


# ===========================================================================
# AHBA
# ===========================================================================
def _check_abagen_available() -> bool:
    try:
        import abagen   # noqa: F401
        return True
    except ImportError:
        return False


def fetch_ahba_expression(target_uniprots: list[str]) -> dict[str, AHBAExpressionData]:
    """Fetch AHBA expression for a target list. Stub returns static §8.6 biases."""
    if not _check_abagen_available():
        logger.info("abagen not installed; returning §8.6 static-bias stub")
        from mammal_repurposing.analysis.brain_region import BRAIN_REGION_BIAS
        out: dict[str, AHBAExpressionData] = {}
        for u in target_uniprots:
            e = BRAIN_REGION_BIAS.get(u, {})
            out[u] = AHBAExpressionData(
                target_uniprot=u,
                gene_symbol=e.get("gene", ""),
                region_expression={e.get("primary_region", "unknown"): 1.0},
                cortical_r={"cortex-biased": 0.5, "subcortical": 0.0,
                            "brainstem": -0.2, "hippocampal": 0.3,
                            "mixed": 0.1}.get(e.get("bias", ""), float("nan")),
                note="stub from §8.6 static brain_region map",
            )
        return out
    raise NotImplementedError(
        "abagen real-mode wiring TBD — see Markello 2021 + §13.2 V6 plan"
    )


# ===========================================================================
# OT Genetics L2G
# ===========================================================================
def _check_ot_reachable() -> bool:
    """Quick HEAD check for OT Genetics GraphQL endpoint."""
    import urllib.request
    try:
        req = urllib.request.Request(
            "https://api.genetics.opentargets.org/graphql",
            headers={"User-Agent": "mammal-repurposing/0.1"},
        )
        with urllib.request.urlopen(req, timeout=3.0) as _:
            return True
    except Exception:
        return False


def fetch_ot_l2g(target_uniprots: list[str]) -> dict[str, OTGeneticsL2G]:
    """Pull max L2G score per target across cognition GWAS studies.

    Studies: Davies 2018 (GCST006269), Hill 2019 (GCST006716), Savage 2018
    (GCST006250), Sniekers 2017, UKBB fluid intelligence.
    """
    if not _check_ot_reachable():
        logger.info("OT Genetics GraphQL unreachable; returning stub")
        return {u: OTGeneticsL2G(
            target_ensembl="", target_uniprot=u, gene_symbol="",
            note="OT reachable=False; stub returned",
        ) for u in target_uniprots}
    raise NotImplementedError(
        "OT Genetics real-mode GraphQL wiring TBD — see V6 plan §13.2 Stage 1"
    )


# ===========================================================================
# cellxgene-census
# ===========================================================================
def _check_tiledbsoma_available() -> bool:
    try:
        import tiledbsoma    # noqa: F401
        import cellxgene_census    # noqa: F401
        return True
    except ImportError:
        return False


# Per Cluster D §E.1 cognition-salient cell-type set
COGNITION_SALIENT_CELLTYPES: list[str] = [
    "L2/3 IT pyramidal",
    "L5 ET pyramidal",
    "L6 CT",
    "CA1 hippocampal",
    "CA3 hippocampal",
    "DG hippocampal",
    "basal forebrain cholinergic (CHAT+)",
    "LC noradrenergic (DBH+)",
    "raphe serotonergic (TPH2+)",
    "VTA dopaminergic (TH+)",
    "PV+ interneurons",
    "SST+ interneurons",
    "VIP+ interneurons",
]


def fetch_celltype_enrichment(
    target_uniprots: list[str],
) -> dict[str, CellTypeEnrichment]:
    """Compute per-target cognition-salient cell-type z-score.

    Stub: returns 0.0 for every target. Real path uses tiledbsoma to slice
    the cellxgene-census brain atlas (Siletti 2023 + Allen + Mathys 2019).
    """
    if not _check_tiledbsoma_available():
        logger.info("tiledbsoma not installed; returning enrichment stub")
        return {u: CellTypeEnrichment(
            target_uniprot=u, gene_symbol="",
            note="tiledbsoma not installed; stub returned",
        ) for u in target_uniprots}
    raise NotImplementedError(
        "cellxgene-census real-mode TileDB-SOMA wiring TBD — see V6 plan §13.2 Stage 1"
    )


def availability() -> dict[str, bool]:
    """Reports which Cluster D data sources are operational at import time."""
    return {
        "abagen": _check_abagen_available(),
        "ot_genetics_reachable": _check_ot_reachable(),
        "tiledbsoma": _check_tiledbsoma_available(),
    }
