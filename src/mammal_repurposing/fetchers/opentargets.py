"""OpenTargets Platform GraphQL fetcher for target/drug context.

Endpoint: https://api.platform.opentargets.org/api/v4/graphql
No auth required. Use a polite User-Agent.

OpenTargets identifies targets by Ensembl gene ID (e.g., ENSG00000175344 for
CHRNA7). We get Ensembl IDs from UniProt via :mod:`fetchers.uniprot`.

Three lookups per target:
    1. Known drugs targeting it (drug name + max trial phase)
    2. Tractability scores (small-molecule / antibody / etc.)
    3. CNS/cognition-adjacent disease associations (filtered by EFO ID)
"""

from __future__ import annotations

import logging
from typing import TypedDict

import httpx
from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential,
)

from mammal_repurposing.config import HTTP_MAX_RETRIES, HTTP_TIMEOUT_SEC, USER_AGENT

logger = logging.getLogger(__name__)

OT_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"

# EFO descendants under "nervous system disease" (EFO_0000618) and "cognitive disorder"
# (EFO_0003900). The full descendant list is large; we instead filter post-hoc on
# therapeutic area names client-side (more robust to ontology drift).
COGNITIVE_THERAPEUTIC_AREAS = {
    "neurological disorder",
    "psychiatric disorder",
    "nervous system disease",
    "cognitive disorder",
    "alzheimer disease",
    "parkinson disease",
    "schizophrenia",
    "depression",
    "attention deficit hyperactivity disorder",
    "autism spectrum disorder",
    "huntington disease",
    "amyotrophic lateral sclerosis",
}


class KnownDrug(TypedDict):
    drug_name: str
    drug_chembl_id: str | None
    max_phase: float | None
    mechanism_of_action: str | None
    disease_name: str | None


class TargetContext(TypedDict):
    target_uniprot: str
    ensembl_id: str
    approved_symbol: str | None
    known_drugs: list[KnownDrug]
    cns_disease_associations: list[str]  # disease names with cognitive/CNS therapeutic area
    n_known_drugs_total: int
    n_cns_drugs: int


_QUERY_TARGET_CONTEXT = """
query TargetContext($ensemblId: String!, $size: Int!) {
  target(ensemblId: $ensemblId) {
    id
    approvedSymbol
    biotype
    associatedDiseases(page: {index: 0, size: $size}) {
      count
      rows {
        score
        disease {
          name
          id
          therapeuticAreas {
            name
          }
        }
      }
    }
  }
}
"""


@retry(
    reraise=True,
    stop=stop_after_attempt(HTTP_MAX_RETRIES),
    wait=wait_exponential(multiplier=1, min=1, max=10),
    retry=retry_if_exception_type((httpx.TransportError, httpx.HTTPStatusError)),
)
def _post_graphql(
    client: httpx.Client,
    query: str,
    variables: dict,
) -> dict:
    resp = client.post(
        OT_GRAPHQL,
        json={"query": query, "variables": variables},
        headers={"Content-Type": "application/json", "User-Agent": USER_AGENT},
    )
    if resp.status_code >= 500:
        resp.raise_for_status()
    if resp.status_code != 200:
        raise httpx.HTTPStatusError(
            f"OpenTargets HTTP {resp.status_code}",
            request=resp.request, response=resp,
        )
    body = resp.json()
    if "errors" in body:
        logger.warning("OpenTargets GraphQL errors: %s", body["errors"])
    return body.get("data", {}) or {}


def _is_cns_disease(disease: dict) -> bool:
    """Heuristic: disease's therapeutic areas include a CNS/cognition keyword."""
    areas = disease.get("therapeuticAreas") or []
    for a in areas:
        name = (a.get("name") or "").lower()
        if any(kw in name for kw in COGNITIVE_THERAPEUTIC_AREAS):
            return True
    name = (disease.get("name") or "").lower()
    return any(kw in name for kw in COGNITIVE_THERAPEUTIC_AREAS)


def fetch_target_context(
    ensembl_id: str,
    target_uniprot: str,
    *,
    size: int = 200,
    client: httpx.Client | None = None,
) -> TargetContext:
    """Pull known drugs + CNS-disease filter for a single target."""
    own_client = client is None
    if own_client:
        client = httpx.Client(timeout=HTTP_TIMEOUT_SEC)

    try:
        data = _post_graphql(client, _QUERY_TARGET_CONTEXT,
                              {"ensemblId": ensembl_id, "size": size})
        target = data.get("target") or {}
        ad_block = target.get("associatedDiseases") or {}
        rows = ad_block.get("rows") or []

        known_drugs: list[KnownDrug] = []  # left empty in v4 schema (knownDrugs field removed)
        cns_diseases: set[str] = set()
        n_cns_total = 0

        for row in rows:
            disease = row.get("disease") or {}
            if _is_cns_disease(disease):
                n_cns_total += 1
                if disease.get("name"):
                    cns_diseases.add(disease["name"])

        return TargetContext(
            target_uniprot=target_uniprot,
            ensembl_id=ensembl_id,
            approved_symbol=target.get("approvedSymbol"),
            known_drugs=known_drugs,
            cns_disease_associations=sorted(cns_diseases),
            n_known_drugs_total=ad_block.get("count") or 0,
            n_cns_drugs=n_cns_total,
        )
    finally:
        if own_client:
            client.close()


def fetch_contexts_for_targets(
    ensembl_uniprot_pairs: list[tuple[str, str]],
    *,
    size: int = 200,
) -> list[TargetContext]:
    """Batch fetcher, reuses one HTTP client. Order preserved."""
    out: list[TargetContext] = []
    with httpx.Client(timeout=HTTP_TIMEOUT_SEC) as client:
        for ensembl, uniprot in ensembl_uniprot_pairs:
            if not ensembl:
                logger.warning("Missing Ensembl ID for %s; skipping OpenTargets lookup.", uniprot)
                continue
            try:
                out.append(fetch_target_context(ensembl, uniprot, size=size, client=client))
            except Exception:
                logger.exception("OpenTargets fetch failed for %s/%s", uniprot, ensembl)
    return out
