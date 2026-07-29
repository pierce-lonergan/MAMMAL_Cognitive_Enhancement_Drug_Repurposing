"""Integrity guard for the healthy-adult cognition ledger.

The 2026-06 robustness audit (scripts/121) found two structural weaknesses in this ledger that were
invisible because nothing enforced them:

  1. **The stated inclusion rule was not checkable.** The ledger declares
     "enhances_healthy_young = 1 iff a clean healthy-adult MA has a CI excluding 0", but nothing
     verified that the labels obeyed it. Exactly one compound (l-theanine: g=+0.35, CI [+0.10,+0.61],
     labelled 0) silently violated it -- and the project's entire headline turned out to hinge on
     that single row.
  2. **Provenance was optional.** Rows could carry an effect size with no PMID/DOI, and 5 of 11
     clean-MA rows had no CI at all, which is what let statistical power masquerade as efficacy.

This module turns both into enforced, testable contracts. It NEVER edits data -- it reports
violations so a human decides. `LEDGER_RULE_DOC` is the single place the rule is written down.
"""
from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

LEDGER_RULE_DOC = (
    "enhances_healthy_young = 1 iff a CLEAN healthy-adult meta-analysis reports a pooled effect "
    "whose confidence interval excludes 0 in a cognitive domain. Any row whose label departs from "
    "this rule MUST carry an explicit justification in `robustness` (e.g. only one sub-domain "
    "significant), so the departure is a recorded judgement rather than an invisible one."
)

# Detecting a PATIENT sample by blocklist, not a whitelist of "healthy" spellings. The ledger uses a
# rich descriptive vocabulary for healthy samples (healthy_nonSD, rested_healthy, nonsmokers,
# healthy_unimpaired, mixed_stressed, ...), and a whitelist flagged 7 of those as suspicious --
# noise that would dilute the real signals. What actually matters is that the sample is not a
# clinical population, so we look for clinical markers instead.
PATIENT_POPULATION_MARKERS = (
    "dementia", "alzheimer", "mci", "mild_cognitive_impairment", "adhd", "schizophreni",
    "parkinson", "depress", "bipolar", "stroke", "tbi", "epilep", "patient", "clinical",
    "impaired_", "diseased",
)
VALID_TIERS = {"clean_MA", "mixed_pop", "contested", "absent"}
VALID_SUPERGROUPS = {"stimulant", "nonstimulant", "other"}


@dataclass
class Violation:
    compound: str
    severity: str          # "error" (blocks) | "warn" (records)
    rule: str
    detail: str

    def __str__(self) -> str:      # pragma: no cover - display only
        return f"[{self.severity.upper()}] {self.compound}: {self.rule} -- {self.detail}"


def _has_ci(row) -> bool:
    return pd.notna(row.get("ci_lo")) and pd.notna(row.get("ci_hi"))


def validate_ledger(df: pd.DataFrame) -> list[Violation]:
    """Check the healthy-adult ledger against its own stated contract.

    Returns a list of Violations (empty = clean). Errors are contract breaches that would corrupt an
    analysis; warnings are recorded weaknesses (e.g. a missing CI) that are permitted but must be
    visible, because a missing CI is exactly what lets study volume masquerade as efficacy.
    """
    v: list[Violation] = []

    dupes = df["compound"][df["compound"].duplicated()].tolist()
    for c in sorted(set(dupes)):
        v.append(Violation(c, "error", "duplicate_compound",
                           "compound appears more than once; the ledger must be one row per compound"))

    for _, r in df.iterrows():
        c = str(r.get("compound", "<missing>"))
        tier = r.get("evidence_tier")

        if tier not in VALID_TIERS:
            v.append(Violation(c, "error", "bad_evidence_tier", f"{tier!r} not in {sorted(VALID_TIERS)}"))
        if r.get("supergroup") not in VALID_SUPERGROUPS:
            v.append(Violation(c, "error", "bad_supergroup", f"{r.get('supergroup')!r}"))

        if tier == "absent":
            # An "absent" row asserts that NO healthy-adult MA exists; it must not carry an effect.
            if pd.notna(r.get("representative_g")):
                v.append(Violation(c, "error", "absent_row_has_effect",
                                   "evidence_tier='absent' means no MA exists, so it cannot report g"))
            if pd.notna(r.get("enhances_healthy_young")):
                v.append(Violation(c, "error", "absent_row_has_label",
                                   "an absent row cannot carry an enhances_healthy_young label"))
            continue

        # --- rows that DO claim evidence -------------------------------------------------------
        if pd.isna(r.get("representative_g")):
            v.append(Violation(c, "error", "missing_effect_size",
                               f"evidence_tier={tier!r} claims evidence but representative_g is null"))
        if not str(r.get("pmid_doi") or "").strip():
            v.append(Violation(c, "error", "missing_provenance",
                               "every row claiming evidence needs a real PMID/DOI; unsourced effect "
                               "sizes are unciteable and unverifiable"))
        pop = str(r.get("population") or "").strip().lower()
        hit = next((m for m in PATIENT_POPULATION_MARKERS if m in pop), None)
        if hit:
            v.append(Violation(c, "error", "patient_population",
                               f"population={pop!r} matches clinical marker {hit!r}; this ledger is "
                               "healthy-adult ground truth and must not contain patient samples "
                               "(mixed samples belong in evidence_tier='mixed_pop')"))
        if not _has_ci(r):
            v.append(Violation(c, "warn", "missing_ci",
                               "no CI recorded: this row cannot distinguish a true null from an "
                               "under-powered one, and inflates the power confound (audit R1/R3)"))

        # --- THE inclusion rule ----------------------------------------------------------------
        if tier == "clean_MA" and _has_ci(r) and pd.notna(r.get("enhances_healthy_young")):
            excludes_zero = bool(r["ci_lo"] > 0)
            label = bool(r["enhances_healthy_young"] == 1)
            if excludes_zero != label:
                justified = len(str(r.get("robustness") or "").strip()) >= 20
                v.append(Violation(
                    c, "warn" if justified else "error", "label_rule_conflict",
                    f"CI [{r['ci_lo']:+.2f}, {r['ci_hi']:+.2f}] excludes_zero={excludes_zero} but "
                    f"label={int(label)}. {LEDGER_RULE_DOC} "
                    + ("A justification IS recorded in `robustness`, so this is a recorded judgement "
                       "-- but the headline's sensitivity to it must be reported (see scripts/121)."
                       if justified else "NO justification is recorded in `robustness`.")))
    return v


def assert_ledger_valid(df: pd.DataFrame) -> None:
    """Raise on any error-severity violation. Warnings are returned by validate_ledger, not raised."""
    errs = [x for x in validate_ledger(df) if x.severity == "error"]
    if errs:
        raise ValueError("healthy-adult ledger failed validation:\n" + "\n".join(str(e) for e in errs))
