"""Persistence-after-cessation axis (the single mechanism-based source of truth).

The cognition ledger and the F2 screen score a SYMPTOMATIC class prior: the clinical-g
a drug delivers WHILE it occupies its target, which reverses on washout. That is a
different axis from a DISEASE-MODIFYING / structurally-persistent effect: a durable
change in trajectory so the patient is better off after STOPPING. "Persists after
cessation" is the second category, and in healthy people it is nearly empty.

This module is the curated, mechanism-keyed source for every persistence-relevant fact
about a compound or class:
  - persistence_status  (the verdict: symptomatic / tested_negative / contested /
                         plasticity_gated / disease_modifying_patients / not_applicable /
                         cognition_negative / demonstrated_healthy / unknown);
  - evidence_design     (the DESIGN behind any persistence claim - a claim is only as
                         good as its design; delayed-start RCT > discontinuation > ... >
                         class_extrapolation, the new weakest tier for a class fact
                         borrowed onto an individual member with no per-compound study);
  - substrate           (the mechanism's persistence CAPABILITY on 3 honest tiers:
                         transient < plasticity_window < ablative - see reversibility.py);
  - self_maintaining    (does the altered state persist after target disengagement? only
                         ablative is True by construction);
  - known_mechanism_class (set ONLY where the L0 structural router is known to misroute,
                         e.g. fluoxetine=SSRI - lets PERSEUS abstain the symptomatic head
                         on a structural/mechanism mismatch instead of emitting a wrong-
                         class prior).

Null by default: a class/compound with no curated entry is `unknown` / `transient`, never
assumed persistent. Curation: data/raw/persistence_axis_{classes,overrides}.csv (every
non-null row cited). numpy/pandas only.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from mammal_repurposing.reporting.trial_watch import _norm_drug

logger = logging.getLogger(__name__)

# persistence_status -> coarse tier
STATUS_TIER: dict[str, str] = {
    "demonstrated_healthy": "live",        # durable enhancement in healthy people (EMPTY)
    "disease_modifying_patients": "live",  # durable slowing in patients, not enhancement
    "contested": "live",                   # delayed-start tested, equivocal (MAO-B/ADAGIO)
    "plasticity_gated": "live",            # durable IF paired with training (SSRI iPlasticity)
    "symptomatic": "null",                 # real on-drug effect, reverses on washout
    "tested_negative": "null",             # formally tested and found negative
    "unknown": "null",                     # no persistence data (default)
    "not_applicable": "exclude",           # no CNS exposure / not a cognition agent
    "cognition_negative": "exclude",       # chronic use impairs cognition
}

# evidence_design strength (high = stronger support for a persistence claim).
# class_extrapolation is the weakest non-null tier: a class fact borrowed onto a member
# with no per-compound study (de-broadcasts the F2-style class-prior artifact).
EVIDENCE_RANK: dict[str, int] = {
    "delayed_start_rct": 7,            # ADAGIO template - gold standard
    "randomized_discontinuation": 6,  # randomized withdrawal / relapse
    "longitudinal_followup": 5,       # cohort / open-label extension / RCT follow-up
    "washout_observation": 4,         # observational washout
    "preclinical_only": 3,            # animal / mechanistic
    "mechanistic_inference": 2,       # pharmacology reasoning (e.g. BBB exclusion)
    "class_extrapolation": 1,         # class fact borrowed onto a member, no per-compound study
    "none": 0,
}

# persistence-substrate capability, 3 honest tiers (see reversibility.py for the rationale)
VALID_SUBSTRATE = ("transient", "plasticity_window", "ablative")


@dataclass
class PersistenceCall:
    status: str
    evidence_design: str
    basis: str
    caveat: str
    source: str
    level: str                              # "compound" | "class" | "default"
    substrate: str = "transient"
    self_maintaining: bool = False
    known_mechanism_class: str | None = None

    @property
    def tier(self) -> str:
        return STATUS_TIER.get(self.status, "null")

    @property
    def evidence_rank(self) -> int:
        return EVIDENCE_RANK.get(self.evidence_design, 0)


_UNKNOWN = PersistenceCall("unknown", "none", "no persistence-specific evidence assessed",
                           "", "", "default")


def _row_to_call(row, level: str) -> PersistenceCall:
    sm = str(row.get("self_maintaining", "no")).strip().lower() in ("yes", "true", "1")
    kmc = row.get("known_mechanism_class", None)
    if kmc is None or (isinstance(kmc, float) and pd.isna(kmc)) or str(kmc).strip() == "":
        kmc = None
    else:
        kmc = str(kmc).strip()
    return PersistenceCall(
        status=str(row["persistence_status"]).strip(),
        evidence_design=str(row["evidence_design"]).strip(),
        basis=str(row.get("basis", "")), caveat=str(row.get("caveat", "")),
        source=str(row.get("source", "")), level=level,
        substrate=(str(row.get("substrate", "transient")).strip() or "transient"),
        self_maintaining=sm, known_mechanism_class=kmc)


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
            if c.substrate not in VALID_SUBSTRATE:
                raise ValueError(f"unknown substrate: {c.substrate!r}")
    return classes, overrides


def call_for(compound: str, mechanism_class: str,
             classes: dict[str, PersistenceCall],
             overrides: dict[str, PersistenceCall]) -> PersistenceCall:
    """Resolve the persistence call: compound override wins, else class default, else
    null (`unknown`)."""
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
    out["persistence_substrate"] = [c.substrate for c in calls]
    out["persistence_self_maintaining"] = [c.self_maintaining for c in calls]
    out["persistence_basis"] = [c.basis for c in calls]
    out["persistence_caveat"] = [c.caveat for c in calls]
    out["persistence_source"] = [c.source for c in calls]
    out["persistence_call_level"] = [c.level for c in calls]
    return out
