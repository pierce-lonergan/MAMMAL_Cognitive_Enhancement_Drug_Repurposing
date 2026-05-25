"""PubChem PUG-REST fetcher for canonical SMILES.

Endpoint: https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/JSON
No auth. Rate limit: 5 requests per second (we use ~4.8 to stay safe).

Important gotcha: the response may return ``ConnectivitySMILES`` under the
``Properties`` key instead of ``CanonicalSMILES`` for some compounds. Both are
valid for our purposes — extract whichever is present.
"""

from __future__ import annotations

import logging
import time
from typing import TypedDict
from urllib.parse import quote

import httpx
from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential,
)

from mammal_repurposing.config import (
    HTTP_MAX_RETRIES,
    HTTP_TIMEOUT_SEC,
    PUBCHEM_RATE_LIMIT_SEC,
    USER_AGENT,
)

logger = logging.getLogger(__name__)

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


class PubChemHit(TypedDict):
    name_queried: str
    cid: int | None
    smiles: str | None
    smiles_kind: str | None  # "canonical" | "connectivity" | None


class _RateLimiter:
    """Simple monotonic-clock spacing throttle. Not thread-safe (we run sync)."""

    def __init__(self, min_interval_sec: float) -> None:
        self._min_interval = min_interval_sec
        self._last_call: float = 0.0

    def wait(self) -> None:
        now = time.monotonic()
        elapsed = now - self._last_call
        if elapsed < self._min_interval:
            time.sleep(self._min_interval - elapsed)
        self._last_call = time.monotonic()


_throttle = _RateLimiter(PUBCHEM_RATE_LIMIT_SEC)


def _extract_smiles(payload: dict) -> tuple[str | None, str | None, int | None]:
    """Return (smiles, smiles_kind, cid) from a PubChem property response."""
    try:
        props = payload["PropertyTable"]["Properties"][0]
    except (KeyError, IndexError):
        return None, None, None
    cid = props.get("CID")
    if "CanonicalSMILES" in props and props["CanonicalSMILES"]:
        return props["CanonicalSMILES"], "canonical", cid
    if "ConnectivitySMILES" in props and props["ConnectivitySMILES"]:
        return props["ConnectivitySMILES"], "connectivity", cid
    if "IsomericSMILES" in props and props["IsomericSMILES"]:
        return props["IsomericSMILES"], "isomeric", cid
    return None, None, cid


@retry(
    reraise=True,
    stop=stop_after_attempt(HTTP_MAX_RETRIES),
    wait=wait_exponential(multiplier=1, min=1, max=10),
    retry=retry_if_exception_type(httpx.TransportError),
)
def _get_property(client: httpx.Client, name: str) -> dict | None:
    """Single-name PubChem property lookup. Returns parsed JSON or None on 404."""
    # PubChem's path component requires URL-encoding for names with commas,
    # spaces, parens, etc. (e.g. "7,8-dihydroxyflavone", "(R)-modafinil").
    safe_name = quote(name, safe="")
    url = (
        f"{PUBCHEM_BASE}/compound/name/{safe_name}"
        "/property/CanonicalSMILES,IsomericSMILES,ConnectivitySMILES/JSON"
    )
    _throttle.wait()
    resp = client.get(url)
    if resp.status_code == 404:
        return None
    if resp.status_code >= 500:
        resp.raise_for_status()  # retry
    if resp.status_code != 200:
        logger.warning("PubChem %s returned HTTP %d", name, resp.status_code)
        return None
    try:
        return resp.json()
    except ValueError:
        logger.warning("PubChem %s returned non-JSON body", name)
        return None


def fetch_smiles(
    name: str,
    *,
    alt_names: list[str] | None = None,
    client: httpx.Client | None = None,
) -> PubChemHit:
    """Resolve a compound name to SMILES via PubChem.

    Tries the primary name first, then each alt_name in order. Returns
    a PubChemHit with smiles=None if nothing resolves.
    """
    own_client = client is None
    if own_client:
        headers = {"User-Agent": USER_AGENT, "Accept": "application/json"}
        client = httpx.Client(timeout=HTTP_TIMEOUT_SEC, headers=headers)

    try:
        candidates = [name] + list(alt_names or [])
        for candidate in candidates:
            if not candidate:
                continue
            payload = _get_property(client, candidate.strip())
            if payload is None:
                continue
            smiles, kind, cid = _extract_smiles(payload)
            if smiles:
                if candidate != name:
                    logger.info("PubChem resolved %r via alt name %r", name, candidate)
                return PubChemHit(
                    name_queried=candidate,
                    cid=cid,
                    smiles=smiles,
                    smiles_kind=kind,
                )
        logger.warning("PubChem could not resolve any name in %r", candidates)
        return PubChemHit(name_queried=name, cid=None, smiles=None, smiles_kind=None)
    finally:
        if own_client:
            client.close()


def fetch_many_smiles(
    names_with_alts: list[tuple[str, list[str]]],
) -> list[PubChemHit]:
    """Resolve a batch of (name, alt_names) tuples, reusing one HTTP client."""
    headers = {"User-Agent": USER_AGENT, "Accept": "application/json"}
    out: list[PubChemHit] = []
    with httpx.Client(timeout=HTTP_TIMEOUT_SEC, headers=headers) as client:
        for name, alts in names_with_alts:
            out.append(fetch_smiles(name, alt_names=alts, client=client))
    return out
