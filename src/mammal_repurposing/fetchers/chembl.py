"""ChEMBL REST fetcher for top binders per target.

Endpoint root: https://www.ebi.ac.uk/chembl/api/data/
No auth. Used to expand the curated compound seed with known high-affinity
binders for each panel target.

Workflow:
    1. Map UniProt accession -> ChEMBL target_chembl_id via /target.json
    2. Pull activities with type IN (Ki, IC50, Kd, EC50) and standard_value below
       a Ki/IC50 cutoff (default 1000 nM = 1 µM) for that target_chembl_id.
    3. For each activity, return molecule_chembl_id + preferred name + canonical SMILES.
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

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"


class ChemblBinder(TypedDict):
    molecule_chembl_id: str
    pref_name: str | None
    smiles: str | None
    activity_type: str
    standard_value_nm: float | None
    target_chembl_id: str


@retry(
    reraise=True,
    stop=stop_after_attempt(HTTP_MAX_RETRIES),
    wait=wait_exponential(multiplier=1, min=1, max=10),
    retry=retry_if_exception_type((httpx.TransportError, httpx.HTTPStatusError)),
)
def _get_json(client: httpx.Client, url: str, params: dict | None = None) -> dict:
    """GET with retry; raises on 5xx (triggering retry)."""
    resp = client.get(url, params=params)
    if resp.status_code >= 500:
        resp.raise_for_status()
    if resp.status_code != 200:
        raise httpx.HTTPStatusError(
            f"ChEMBL HTTP {resp.status_code} for {resp.request.url}",
            request=resp.request,
            response=resp,
        )
    return resp.json()


def uniprot_to_chembl_target(
    accession: str,
    *,
    client: httpx.Client | None = None,
) -> str | None:
    """Resolve a UniProt accession to a single ChEMBL target_chembl_id.

    Picks the first SINGLE PROTEIN target with the matching component accession.
    Returns None if no mapping found.
    """
    own_client = client is None
    if own_client:
        headers = {"User-Agent": USER_AGENT, "Accept": "application/json"}
        client = httpx.Client(timeout=HTTP_TIMEOUT_SEC, headers=headers)

    try:
        url = f"{CHEMBL_BASE}/target.json"
        params = {
            "target_components__accession": accession,
            "target_type": "SINGLE PROTEIN",
            "limit": 5,
        }
        try:
            payload = _get_json(client, url, params=params)
        except httpx.HTTPStatusError as e:
            logger.warning("ChEMBL target lookup failed for %s: %s", accession, e)
            return None

        targets = payload.get("targets", [])
        if not targets:
            logger.warning("ChEMBL no SINGLE PROTEIN target for UniProt %s", accession)
            return None
        return targets[0]["target_chembl_id"]
    finally:
        if own_client:
            client.close()


def top_binders(
    target_chembl_id: str,
    *,
    n: int = 15,
    max_standard_nm: float = 1000.0,
    client: httpx.Client | None = None,
) -> list[ChemblBinder]:
    """Return the top-N most-potent binders for a ChEMBL target.

    Pulls activities with standard_type in {Ki, IC50, Kd, EC50}, standard_units = nM,
    standard_value <= max_standard_nm. Sorts by standard_value ascending (most potent
    first). Deduplicates by molecule_chembl_id.
    """
    own_client = client is None
    if own_client:
        headers = {"User-Agent": USER_AGENT, "Accept": "application/json"}
        client = httpx.Client(timeout=HTTP_TIMEOUT_SEC, headers=headers)

    try:
        # Pull enough activities to dedupe down to N molecules. 4x is a safe buffer.
        params = {
            "target_chembl_id": target_chembl_id,
            "standard_type__in": "Ki,IC50,Kd,EC50",
            "standard_units": "nM",
            "standard_value__lte": max_standard_nm,
            "order_by": "standard_value",
            "limit": n * 4,
        }
        try:
            payload = _get_json(client, f"{CHEMBL_BASE}/activity.json", params=params)
        except httpx.HTTPStatusError as e:
            logger.warning("ChEMBL activity query failed for %s: %s", target_chembl_id, e)
            return []

        activities = payload.get("activities", [])
        seen: set[str] = set()
        binders: list[ChemblBinder] = []
        for act in activities:
            mol_id = act.get("molecule_chembl_id")
            if not mol_id or mol_id in seen:
                continue
            seen.add(mol_id)
            binders.append(
                ChemblBinder(
                    molecule_chembl_id=mol_id,
                    pref_name=act.get("molecule_pref_name"),
                    smiles=act.get("canonical_smiles"),
                    activity_type=act.get("standard_type") or "?",
                    standard_value_nm=(
                        float(act["standard_value"]) if act.get("standard_value") else None
                    ),
                    target_chembl_id=target_chembl_id,
                )
            )
            if len(binders) >= n:
                break
        return binders
    finally:
        if own_client:
            client.close()


def top_binders_for_targets(
    uniprots: list[str],
    *,
    per_target: int = 15,
    max_standard_nm: float = 1000.0,
) -> list[ChemblBinder]:
    """One-shot helper: resolve each UniProt, fetch top binders, return flat list.

    Reuses a single HTTP client. Per-target failures (5xx, timeouts) are logged
    and the target is skipped — we'd rather lose 10-15 binders from one target
    than abort the whole compound-library build because ChEMBL is flaky.
    """
    out: list[ChemblBinder] = []
    headers = {"User-Agent": USER_AGENT, "Accept": "application/json"}
    with httpx.Client(timeout=HTTP_TIMEOUT_SEC, headers=headers) as client:
        for acc in uniprots:
            try:
                chembl_id = uniprot_to_chembl_target(acc, client=client)
            except Exception as e:
                logger.warning("ChEMBL target lookup raised for %s: %s; skipping.", acc, e)
                continue
            if chembl_id is None:
                continue
            try:
                binders = top_binders(
                    chembl_id,
                    n=per_target,
                    max_standard_nm=max_standard_nm,
                    client=client,
                )
            except Exception as e:
                logger.warning("ChEMBL binders fetch raised for %s/%s: %s; skipping.",
                               chembl_id, acc, e)
                continue
            logger.info("ChEMBL %s (UniProt %s): %d binders", chembl_id, acc, len(binders))
            out.extend(binders)
    return out
