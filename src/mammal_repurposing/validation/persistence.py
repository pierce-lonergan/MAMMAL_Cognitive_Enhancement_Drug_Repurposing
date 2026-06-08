"""Persistence-after-cessation axis.

The cognition ledger and the F2 screen score a SYMPTOMATIC class prior: the
clinical-g a drug delivers WHILE it occupies its target, which reverses on washout.
That is a different axis from a DISEASE-MODIFYING / structurally-persistent effect:
a durable change in trajectory so the patient is better off after STOPPING than they
would otherwise be. "Permanent gain that persists after cessation" is the second
category, and across psychopharmacology it is nearly empty in healthy people - almost
every cognition drug delivers a state-dependent, reversible effect.

This module makes that distinction first-class. It is null-by-default: a class or
compound with no persistence evidence is `unknown`, not assumed persistent. It also
encodes the EVIDENCE-DESIGN hierarchy, because a persistence claim is only as good as
the design that tested it - a randomized delayed-start trial (the ADAGIO template) is
the gold standard, a randomized-discontinuation/relapse study is next, and mechanistic
inference ("doesn't cross the BBB") is weakest.

Curation lives in `data/raw/persistence_axis_classes.csv` (mechanism-class defaults)
and `persistence_axis_overrides.csv` (compound-level, mainly the structure-router
misroutes). Every row carries a cited basis; nothing is fabricated and the rare
non-null calls (MAO-B contested, fluoxetine plasticity-gated) are hedged to match the
evidence. numpy/pandas only, fully testable in CI.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from mammal_repurposing.reporting.trial_watch import _norm_drug

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Vocabulary
# ---------------------------------------------------------------------------
# persistence_status: the verdict on a durable post-cessation cognitive gain.
STATUS_TIER: dict[str, str] = {
    # LIVE - persistence is plausible / formally tested and worth tracking
    "demonstrated_healthy": "live",       # durable enhancement in healthy people (EMPTY today)
    "disease_modifying_patients": "live",  # durable slowing in patients, not enhancement
    "contested": "live",                   # delayed-start tested, equivocal (MAO-B/ADAGIO)
    "plasticity_gated": "live",            # durable IF paired with training (SSRI iPlasticity)
    # NULL - no persistence; symptomatic or untested
    "symptomatic": "null",                 # reverses on washout
    "tested_negative": "null",             # discontinuation/washout showed loss
    "unknown": "null",                     # no persistence data (default)
    # EXCLUDE - not a valid central cognition agent at all
    "not_applicable": "exclude",           # no CNS exposure / wrong mechanism (routing artifact)
    "cognition_negative": "exclude",       # chronic use impairs cognition
}

# evidence_design: strength of the design behind a persistence claim (high = stronger).
EVIDENCE_RANK: dict[str, int] = {
    "delayed_start_rct": 6,            # ADAGIO template - gold standard for disease-modification
    "randomized_discontinuation": 5,  # randomized withdrawal / relapse
    "longitudinal_followup": 4,       # long-term cohort / RCT follow-up
    "washout_observation": 3,         # observational washout
    "preclinical_only": 2,            # animal / mechanistic
    "mechanistic_inference": 1,       # pharmacology-based reasoning (e.g. BBB exclusion)
    "none": 0,
}


@dataclass
class PersistenceCall:
    status: str
    evidence_design: str
    basis: str
    caveat: str
    source: str
    level: str  # "compound" | "class" | "default"

    @property
    def tier(self) -> str:
        return STATUS_TIER.get(self.status, "null")

    @property
    def evidence_rank(self) -> int:
        return EVIDENCE_RANK.get(self.evidence_design, 0)


_UNKNOWN = PersistenceCall("unknown", "none", "no persistence-specific evidence assessed",
                           "", "", "default")


# ---------------------------------------------------------------------------
# Load curation
# ---------------------------------------------------------------------------

def _row_to_call(row, level: str) -> PersistenceCall:
    return PersistenceCall(
        status=str(row["persistence_status"]).strip(),
        evidence_design=str(row["evidence_design"]).strip(),
        basis=str(row.get("basis", "")), caveat=str(row.get("caveat", "")),
        source=str(row.get("source", "")), level=level)


def load_persistence(classes_csv, overrides_csv
                     ) -> tuple[dict[str, PersistenceCall], dict[str, PersistenceCall]]:
    """Return (class_defaults_by_mechanism_class, overrides_by_normalised_compound)."""
    classes: dict[str, PersistenceCall] = {}
    if Path(classes_csv).exists():
        cdf = pd.read_csv(classes_csv)
        for _, r in cdf.iterrows():
            classes[str(r["mechanism_class"]).strip()] = _row_to_call(r, "class")
    overrides: dict[str, PersistenceCall] = {}
    if Path(overrides_csv).exists():
        odf = pd.read_csv(overrides_csv)
        for _, r in odf.iterrows():
            overrides[_norm_drug(str(r["compound"]))] = _row_to_call(r, "compound")
    # validate vocabulary so a typo never silently becomes a tier
    for d in (classes, overrides):
        for c in d.values():
            if c.status not in STATUS_TIER:
                raise ValueError(f"unknown persistence_status: {c.status!r}")
            if c.evidence_design not in EVIDENCE_RANK:
                raise ValueError(f"unknown evidence_design: {c.evidence_design!r}")
    return classes, overrides


def call_for(compound: str, mechanism_class: str,
             classes: dict[str, PersistenceCall],
             overrides: dict[str, PersistenceCall]) -> PersistenceCall:
    """Resolve the persistence call: compound override wins, else class default,
    else null (`unknown`)."""
    ov = overrides.get(_norm_drug(str(compound)))
    if ov is not None:
        return ov
    cl = classes.get(str(mechanism_class))
    if cl is not None:
        return cl
    return _UNKNOWN


def annotate(df: pd.DataFrame, classes: dict[str, PersistenceCall],
             overrides: dict[str, PersistenceCall], *,
             name_col: str = "query_id",
             class_col: str = "assigned_class") -> pd.DataFrame:
    """Add persistence columns to a routed shortlist (compound + assigned class)."""
    out = df.copy()
    calls = [call_for(r[name_col], r[class_col], classes, overrides)
             for _, r in out.iterrows()]
    out["persistence_status"] = [c.status for c in calls]
    out["persistence_tier"] = [c.tier for c in calls]
    out["persistence_evidence"] = [c.evidence_design for c in calls]
    out["persistence_basis"] = [c.basis for c in calls]
    out["persistence_caveat"] = [c.caveat for c in calls]
    out["persistence_source"] = [c.source for c in calls]
    out["persistence_call_level"] = [c.level for c in calls]
    return out
