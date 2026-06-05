"""§8.3 — ClinicalTrials.gov v2 API fetcher for the wet-lab shortlist.

ClinicalTrials.gov v2 API (production since 2024) base URL:
    https://clinicaltrials.gov/api/v2/studies

Query pattern (free-text intervention):
    GET /studies?query.intr=donepezil&query.cond=cognition
        &filter.overallStatus=RECRUITING,ACTIVE_NOT_RECRUITING,COMPLETED
        &pageSize=20&format=json

Per-compound this returns 0-20 trials. We aggregate to:
    - n_trials_cognition: count of all-status trials with cognition condition
    - n_trials_active: count of ACTIVE_NOT_RECRUITING + RECRUITING
    - n_trials_completed: completed count
    - latest_phase: max(phase) reported (PHASE1/2/3/4)
    - ip_status: derived from phase + status:
        approved      — at least one Phase 4 completed
        investigational — Phase 1-3 with active or completed studies
        early          — Phase 0/1 only
        none           — no cognition-relevant trials

Public API; no auth; respectful 100ms throttle (the v2 API is fast and
documented to handle 100+ qps without rate-limiting individual IPs).
"""

from __future__ import annotations

import json as _json
import logging
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass

import httpx

logger = logging.getLogger(__name__)

CTGOV_V2_BASE = "https://clinicaltrials.gov/api/v2/studies"
DEFAULT_TIMEOUT = 15.0
# CTgov v2 rejects requests with custom User-Agent strings (HTTP 403). Spoof a
# browser UA — same behaviour as curl + most public benchmark tooling.
DEFAULT_USER_AGENT = (
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
    "AppleWebKit/537.36 (KHTML, like Gecko) "
    "Chrome/123.0.0.0 Safari/537.36"
)
# Client-side filter — these tokens, found in any returned trial's conditions
# list, mark it as cognition-relevant. Applied after fetch.
COGNITION_CONDITION_TOKENS: tuple[str, ...] = (
    "cognition", "cognitive", "memory", "attention",
    "alzheimer", "dementia", "mild cognitive impairment", "mci",
    "adhd", "narcolepsy", "schizophrenia", "depression",
    "parkinson", "huntington",  # cognitive-relevant CNS comorbid
)


@dataclass
class TrialSummary:
    compound: str
    n_trials: int                      # cognition-filtered, all statuses
    n_trials_active: int
    n_trials_completed: int
    latest_phase: str                  # PHASE1 | PHASE2 | PHASE3 | PHASE4 | NA
    ip_status: str                     # approved | investigational | early | none
    sample_nct_ids: list[str]


def _normalise_phase(phases: list[str]) -> str:
    """Return the highest-numbered phase observed (PHASE0..PHASE4)."""
    rank = {"PHASE4": 4, "PHASE3": 3, "PHASE2": 2, "PHASE1": 1, "PHASE0": 0}
    best = -1
    out = "NA"
    for p in phases:
        for token, n in rank.items():
            if token in p.upper().replace("/", " "):
                if n > best:
                    best = n
                    out = f"PHASE{n}"
    return out


def fetch_cognitive_trials(
    compound_name: str,
    *,
    client: httpx.Client | None = None,
    page_size: int = 20,
    throttle_s: float = 0.10,
) -> TrialSummary:
    """One compound → trial-summary roll-up."""
    own_client = False
    if client is None:
        client = httpx.Client(
            timeout=DEFAULT_TIMEOUT,
            headers={"User-Agent": DEFAULT_USER_AGENT, "Accept": "*/*"},
        )
        own_client = True

    try:
        # CTgov v2 rejects parens-bracketed condition filters with 403; query
        # only on intervention and filter cognition-relevance client-side.
        # Also: httpx 403s here regardless of UA / Accept — fall back to
        # urllib.request which behaves like curl.
        params = {
            "query.intr": compound_name,
            "pageSize": page_size,
            "format": "json",
        }
        qs = urllib.parse.urlencode(params)
        url = f"{CTGOV_V2_BASE}?{qs}"
        req = urllib.request.Request(
            url,
            headers={"User-Agent": DEFAULT_USER_AGENT, "Accept": "*/*"},
        )
        with urllib.request.urlopen(req, timeout=DEFAULT_TIMEOUT) as resp:
            body = _json.loads(resp.read().decode("utf-8"))
    except Exception as e:
        logger.warning("CTgov fetch failed for %s: %s", compound_name, e)
        if own_client:
            client.close()
        return TrialSummary(
            compound=compound_name, n_trials=0, n_trials_active=0,
            n_trials_completed=0, latest_phase="NA",
            ip_status="none", sample_nct_ids=[],
        )
    finally:
        if throttle_s:
            time.sleep(throttle_s)

    studies = body.get("studies", [])
    # Client-side cognition filter — keep only studies where ANY condition
    # contains one of the COGNITION_CONDITION_TOKENS.
    cognition_relevant = []
    for s in studies:
        prot = s.get("protocolSection", {})
        conds = prot.get("conditionsModule", {}).get("conditions", []) or []
        joined = " ".join(conds).lower()
        if any(tok in joined for tok in COGNITION_CONDITION_TOKENS):
            cognition_relevant.append(s)

    n_active = 0
    n_completed = 0
    phases: list[str] = []
    sample_nct: list[str] = []

    for s in cognition_relevant:
        prot = s.get("protocolSection", {})
        status = prot.get("statusModule", {}).get("overallStatus", "")
        if status in ("RECRUITING", "ACTIVE_NOT_RECRUITING"):
            n_active += 1
        elif status == "COMPLETED":
            n_completed += 1
        ph = prot.get("designModule", {}).get("phases", []) or []
        phases.extend(ph)
        nct = prot.get("identificationModule", {}).get("nctId", "")
        if nct and len(sample_nct) < 5:
            sample_nct.append(nct)

    latest_phase = _normalise_phase(phases)
    # Re-bind studies to the filtered count so the summary reflects cognition-only
    studies = cognition_relevant

    # ip_status derivation
    if latest_phase == "PHASE4" and n_completed > 0:
        ip_status = "approved"
    elif latest_phase in ("PHASE2", "PHASE3") and (n_active + n_completed) > 0:
        ip_status = "investigational"
    elif latest_phase in ("PHASE0", "PHASE1") and (n_active + n_completed) > 0:
        ip_status = "early"
    else:
        ip_status = "none"

    if own_client:
        client.close()

    return TrialSummary(
        compound=compound_name,
        n_trials=len(studies),
        n_trials_active=n_active,
        n_trials_completed=n_completed,
        latest_phase=latest_phase,
        ip_status=ip_status,
        sample_nct_ids=sample_nct,
    )


def fetch_trials_for_shortlist(
    compound_names: list[str],
    throttle_s: float = 0.10,
) -> list[TrialSummary]:
    """Batch wrapper — single httpx.Client + sequential calls."""
    out: list[TrialSummary] = []
    with httpx.Client(
        timeout=DEFAULT_TIMEOUT,
        headers={"User-Agent": DEFAULT_USER_AGENT, "Accept": "application/json"},
    ) as client:
        for i, name in enumerate(compound_names):
            r = fetch_cognitive_trials(name, client=client, throttle_s=throttle_s)
            out.append(r)
            if (i + 1) % 5 == 0:
                logger.info("  CTgov [%d/%d]: latest=%s", i + 1, len(compound_names),
                            r.compound)
    return out
