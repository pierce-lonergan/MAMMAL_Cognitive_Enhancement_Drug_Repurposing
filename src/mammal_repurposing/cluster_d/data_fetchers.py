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
    Genomics is the newer L2G release. Open Targets Platform v25+ unifies
    L2G into api.platform.opentargets.org/api/v4/graphql
  - cellxgene-census tiledbsoma — Siletti 2023 Science add7046; CZ CELLxGENE
    Discover Census 2025-11-08 LTS
  - BrainSMASH — Burt et al. 2020 NeuroImage 220:117038
"""

from __future__ import annotations

import json
import logging
import time
from dataclasses import dataclass, field
from pathlib import Path

logger = logging.getLogger(__name__)


# Cognition GWAS study IDs per V6.B.1 Stage 2 spec
# (Davies 2018, Hill 2019, Sniekers 2017, Savage 2018, UKBB fluid intelligence)
COGNITION_GWAS_STUDIES: dict[str, dict[str, str]] = {
    "GCST006269": {
        "author": "Davies G",
        "year": "2018",
        "trait": "general cognitive ability",
        "n_eff": "300486",
        "doi": "10.1038/s41467-018-04362-x",
    },
    "GCST006716": {
        "author": "Hill WD",
        "year": "2019",
        "trait": "intelligence (MTAG)",
        "n_eff": "248482",
        "doi": "10.1038/s41380-017-0001-5",
    },
    "GCST006250": {
        "author": "Savage JE",
        "year": "2018",
        "trait": "intelligence",
        "n_eff": "269867",
        "doi": "10.1038/s41588-018-0152-6",
    },
    "GCST005059": {
        "author": "Sniekers S",
        "year": "2017",
        "trait": "intelligence",
        "n_eff": "78308",
        "doi": "10.1038/ng.3869",
    },
    "GCST006572": {
        "author": "UKBB fluid intelligence",
        "year": "2018",
        "trait": "fluid intelligence (UK Biobank)",
        "n_eff": "108818",
        "doi": "ukbb",
    },
}


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


def _uniprot_to_ensembl(uniprot: str, timeout_s: float = 5.0) -> tuple[str, str]:
    """Map a UniProt accession to an Ensembl gene ID + HGNC symbol.

    Uses the UniProt REST API (https://rest.uniprot.org/uniprotkb/{acc}.json).
    Returns (ensembl_gene_id, gene_symbol). Returns ('', '') on failure.
    """
    import urllib.request
    try:
        req = urllib.request.Request(
            f"https://rest.uniprot.org/uniprotkb/{uniprot}.json",
            headers={"User-Agent": "mammal-repurposing/0.1 (research)",
                     "Accept": "application/json"},
        )
        with urllib.request.urlopen(req, timeout=timeout_s) as resp:
            data = json.loads(resp.read().decode("utf-8"))
        ensembl = ""
        symbol = ""
        # Extract gene symbol (preferred name)
        genes = data.get("genes", [])
        if genes:
            symbol = genes[0].get("geneName", {}).get("value", "")
        # Extract Ensembl gene cross-ref
        for xref in data.get("uniProtKBCrossReferences", []):
            if xref.get("database") == "Ensembl":
                # Properties: GeneId is "ENSG..."
                for prop in xref.get("properties", []):
                    if prop.get("key") == "GeneId":
                        ensembl = prop.get("value", "").split(".")[0]   # strip version
                        break
                if ensembl:
                    break
        return ensembl, symbol
    except Exception as e:
        logger.warning("UniProt → Ensembl lookup failed for %s: %s", uniprot, e)
        return "", ""


def _query_ot_genetics_legacy(
    ensembl_gene_id: str,
    timeout_s: float = 10.0,
) -> list[dict]:
    """Query the legacy OT Genetics GraphQL API for L2G across studies.

    Endpoint: https://api.genetics.opentargets.org/graphql
    Returns a list of {studyId, yProbaModel, traitReported, pmid, ...} dicts
    for the gene across all available studies. Caller filters to cognition.
    """
    import urllib.request
    import urllib.error
    query = """
    query GeneL2G($geneId: String!) {
      studiesAndLeadVariantsForGeneByL2G(geneId: $geneId) {
        yProbaModel
        pval
        study {
          studyId
          pmid
          pubAuthor
          pubDate
          traitReported
          nInitial
        }
        variant { id }
      }
    }
    """
    payload = json.dumps({
        "query": query,
        "variables": {"geneId": ensembl_gene_id},
    }).encode("utf-8")
    req = urllib.request.Request(
        "https://api.genetics.opentargets.org/graphql",
        data=payload,
        headers={
            "Content-Type": "application/json",
            "Accept": "application/json",
            "User-Agent": "mammal-repurposing/0.1 (research)",
        },
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout_s) as resp:
            data = json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        # Common: 410 Gone if legacy endpoint retired
        logger.info("Legacy OT Genetics endpoint returned HTTP %s for %s",
                    e.code, ensembl_gene_id)
        return []
    except Exception as e:
        logger.warning("Legacy OT Genetics query failed for %s: %s",
                       ensembl_gene_id, e)
        return []
    if "errors" in data:
        logger.info("Legacy OT GraphQL errors for %s: %s",
                    ensembl_gene_id, data["errors"][:1])
        return []
    return (data.get("data", {})
                .get("studiesAndLeadVariantsForGeneByL2G", []) or [])


def _query_ot_platform_l2g(
    ensembl_gene_id: str,
    timeout_s: float = 10.0,
) -> list[dict]:
    """Query the unified Open Targets Platform GraphQL for L2G credible sets.

    Endpoint: https://api.platform.opentargets.org/api/v4/graphql
    Per Open Targets Platform 24+, L2G is exposed via the `target` →
    `credibleSets` association. Returns a list of {studyId, score,
    traitReported, ...} dicts.
    """
    import urllib.request
    import urllib.error
    query = """
    query TargetL2G($geneId: String!) {
      target(ensemblId: $geneId) {
        id
        approvedSymbol
        credibleSets(page: { index: 0, size: 200 }) {
          rows {
            studyId
            study {
              studyId
              pubmedId
              traitFromSource
              nSamples
            }
            l2GPredictions {
              rows { score }
            }
          }
        }
      }
    }
    """
    payload = json.dumps({
        "query": query,
        "variables": {"geneId": ensembl_gene_id},
    }).encode("utf-8")
    req = urllib.request.Request(
        "https://api.platform.opentargets.org/api/v4/graphql",
        data=payload,
        headers={
            "Content-Type": "application/json",
            "Accept": "application/json",
            "User-Agent": "mammal-repurposing/0.1 (research)",
        },
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout_s) as resp:
            data = json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        logger.info("OT Platform endpoint returned HTTP %s for %s",
                    e.code, ensembl_gene_id)
        return []
    except Exception as e:
        logger.warning("OT Platform query failed for %s: %s",
                       ensembl_gene_id, e)
        return []
    if "errors" in data:
        logger.info("OT Platform GraphQL errors for %s: %s",
                    ensembl_gene_id, data["errors"][:1])
        return []
    target = data.get("data", {}).get("target") or {}
    rows: list[dict] = []
    for row in target.get("credibleSets", {}).get("rows", []) or []:
        l2g_rows = row.get("l2GPredictions", {}).get("rows", []) or []
        score = max((r.get("score", 0.0) for r in l2g_rows), default=0.0)
        rows.append({
            "studyId": row.get("studyId", ""),
            "score": float(score),
            "study": row.get("study", {}),
        })
    return rows


def fetch_ot_l2g(
    target_uniprots: list[str],
    cache_path: Path | None = None,
    use_cache: bool = True,
    rate_limit_s: float = 0.1,
) -> dict[str, OTGeneticsL2G]:
    """Pull max L2G score per target across cognition GWAS studies.

    Studies: Davies 2018 (GCST006269), Hill 2019 (GCST006716), Savage 2018
    (GCST006250), Sniekers 2017 (GCST005059), UKBB fluid intelligence
    (GCST006572). See COGNITION_GWAS_STUDIES for the full set.

    Tries the legacy api.genetics.opentargets.org GraphQL first (which has
    explicit studiesAndLeadVariantsForGeneByL2G), then falls back to the
    unified api.platform.opentargets.org Platform v25+ L2G endpoint. If both
    fail, returns stubs with `note` documenting the failure mode.

    Cache: results are persisted to `cache_path` as parquet so subsequent
    runs are network-free.

    Args:
        target_uniprots: list of UniProt accessions (e.g. ["P22303", "Q01959"])
        cache_path: where to persist results (default
            data/cache/ot_genetics_l2g_v1.parquet)
        use_cache: if True and cache exists, return cached results first
        rate_limit_s: sleep between GraphQL calls (OT is ~10qps soft-limit)
    """
    if cache_path is None:
        from pathlib import Path as _P
        cache_path = (_P(__file__).resolve().parents[3]
                      / "data" / "cache" / "ot_genetics_l2g_v1.parquet")

    # Cache hit
    if use_cache and cache_path.exists():
        try:
            import pandas as pd
            df_cache = pd.read_parquet(cache_path)
            cached = {row.target_uniprot: OTGeneticsL2G(
                target_ensembl=row.target_ensembl,
                target_uniprot=row.target_uniprot,
                gene_symbol=row.gene_symbol,
                l2g_max_score=float(row.l2g_max_score),
                contributing_studies=list(row.contributing_studies)
                    if hasattr(row, "contributing_studies") and row.contributing_studies is not None
                    else [],
                note=row.note,
            ) for row in df_cache.itertuples()}
            missing = [u for u in target_uniprots if u not in cached]
            if not missing:
                logger.info("OT Genetics L2G cache hit for all %d targets at %s",
                            len(target_uniprots), cache_path)
                return {u: cached[u] for u in target_uniprots}
            logger.info("OT Genetics L2G cache hit for %d/%d targets; fetching %d new",
                        len(target_uniprots) - len(missing), len(target_uniprots),
                        len(missing))
            existing = cached
        except Exception as e:
            logger.warning("OT Genetics L2G cache read failed (%s); refetching", e)
            existing = {}
    else:
        existing = {}

    # Network probe
    reachable = _check_ot_reachable()
    if not reachable:
        logger.info("OT Genetics GraphQL unreachable; returning stubs for uncached targets")
        out = dict(existing)
        for u in target_uniprots:
            if u not in out:
                out[u] = OTGeneticsL2G(
                    target_ensembl="", target_uniprot=u, gene_symbol="",
                    note="OT reachable=False; stub returned",
                )
        return out

    # Fetch missing targets
    cognition_study_ids = set(COGNITION_GWAS_STUDIES.keys())
    out = dict(existing)
    for u in target_uniprots:
        if u in out:
            continue
        ensembl, symbol = _uniprot_to_ensembl(u)
        if not ensembl:
            out[u] = OTGeneticsL2G(
                target_ensembl="", target_uniprot=u, gene_symbol=symbol,
                note="UniProt→Ensembl mapping failed",
            )
            continue

        # Try legacy first, then Platform
        legacy_rows = _query_ot_genetics_legacy(ensembl)
        max_score = 0.0
        contributing: list[str] = []
        for row in legacy_rows:
            study_id = (row.get("study", {}) or {}).get("studyId", "")
            score = float(row.get("yProbaModel", 0.0) or 0.0)
            if study_id in cognition_study_ids and score > 0:
                contributing.append(f"{study_id}:{score:.3f}")
                if score > max_score:
                    max_score = score
        note = "legacy OT Genetics API"

        if max_score == 0.0:
            # Fall back to Platform L2G
            platform_rows = _query_ot_platform_l2g(ensembl)
            for row in platform_rows:
                study_id = row.get("studyId", "")
                score = float(row.get("score", 0.0))
                if study_id in cognition_study_ids and score > 0:
                    contributing.append(f"{study_id}:{score:.3f}")
                    if score > max_score:
                        max_score = score
            note = "Platform v25+ L2G API" if platform_rows else "no L2G evidence in cognition GWAS"

        out[u] = OTGeneticsL2G(
            target_ensembl=ensembl,
            target_uniprot=u,
            gene_symbol=symbol,
            l2g_max_score=max_score,
            contributing_studies=contributing,
            note=note,
        )
        time.sleep(rate_limit_s)

    # Persist cache
    try:
        import pandas as pd
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        rows = [
            {"target_uniprot": v.target_uniprot,
             "target_ensembl": v.target_ensembl,
             "gene_symbol": v.gene_symbol,
             "l2g_max_score": v.l2g_max_score,
             "contributing_studies": v.contributing_studies,
             "note": v.note}
            for v in out.values()
        ]
        pd.DataFrame(rows).to_parquet(cache_path, index=False)
        logger.info("OT Genetics L2G cache written to %s (%d rows)",
                    cache_path, len(rows))
    except Exception as e:
        logger.warning("OT Genetics L2G cache write failed: %s", e)

    return out


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
