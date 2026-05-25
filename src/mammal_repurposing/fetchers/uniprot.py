"""UniProt REST fetcher for amino-acid sequences.

Endpoint: https://rest.uniprot.org/uniprotkb/{accession}.json
No authentication required. Polite User-Agent set per :mod:`mammal_repurposing.config`.
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

UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"


class UniprotEntry(TypedDict):
    accession: str
    sequence: str
    length: int
    gene_name: str | None
    ensembl_gene_id: str | None


class UniprotFetchError(RuntimeError):
    """Raised when a UniProt fetch fails after retries."""


def _extract_ensembl_gene_id(entry: dict) -> str | None:
    """Pull the first Ensembl gene ID from a UniProt JSON xrefs block, if any."""
    for xref in entry.get("uniProtKBCrossReferences", []):
        if xref.get("database") == "Ensembl":
            for prop in xref.get("properties", []):
                if prop.get("key") == "GeneId":
                    return prop.get("value")
    return None


def _extract_gene_name(entry: dict) -> str | None:
    genes = entry.get("genes", [])
    if not genes:
        return None
    first = genes[0].get("geneName") or {}
    return first.get("value")


@retry(
    reraise=True,
    stop=stop_after_attempt(HTTP_MAX_RETRIES),
    wait=wait_exponential(multiplier=1, min=1, max=10),
    retry=retry_if_exception_type((httpx.TransportError, httpx.HTTPStatusError)),
)
def fetch_sequence(
    accession: str,
    *,
    client: httpx.Client | None = None,
) -> UniprotEntry:
    """Fetch a UniProt entry and return the AA sequence plus metadata.

    Args:
        accession: UniProt accession (e.g. "P36544").
        client: optional pre-configured httpx.Client (for testing / connection reuse).

    Returns:
        UniprotEntry with sequence, length, gene_name, ensembl_gene_id.

    Raises:
        UniprotFetchError: on non-recoverable HTTP errors or malformed responses.
    """
    url = f"{UNIPROT_BASE}/{accession}.json"
    headers = {"User-Agent": USER_AGENT, "Accept": "application/json"}

    own_client = client is None
    if own_client:
        client = httpx.Client(timeout=HTTP_TIMEOUT_SEC, headers=headers)

    try:
        response = client.get(url, headers=headers)
        if response.status_code >= 500:
            response.raise_for_status()  # triggers tenacity retry
        if response.status_code == 404:
            raise UniprotFetchError(f"UniProt accession not found: {accession}")
        if response.status_code != 200:
            raise UniprotFetchError(
                f"UniProt fetch failed for {accession}: HTTP {response.status_code}"
            )

        entry = response.json()
        seq_block = entry.get("sequence") or {}
        sequence = seq_block.get("value")
        length = seq_block.get("length")
        if not sequence or not length:
            raise UniprotFetchError(f"UniProt response missing sequence for {accession}")

        return UniprotEntry(
            accession=accession,
            sequence=sequence,
            length=int(length),
            gene_name=_extract_gene_name(entry),
            ensembl_gene_id=_extract_ensembl_gene_id(entry),
        )
    finally:
        if own_client:
            client.close()


def fetch_many(
    accessions: list[str],
) -> list[UniprotEntry]:
    """Fetch sequences for a list of accessions, reusing one HTTP client.

    Errors on individual accessions propagate — we want a hard failure on
    target-panel construction rather than a silently incomplete panel.
    """
    headers = {"User-Agent": USER_AGENT, "Accept": "application/json"}
    with httpx.Client(timeout=HTTP_TIMEOUT_SEC, headers=headers) as client:
        return [fetch_sequence(acc, client=client) for acc in accessions]
