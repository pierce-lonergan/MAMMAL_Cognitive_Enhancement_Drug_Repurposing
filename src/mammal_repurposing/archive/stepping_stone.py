"""The stepping-stone archive: what this project killed, and what would bring it back.

WHY THIS EXISTS. Before 2026-09-20 this repository had thirty-eight ledgers and every one of them
held positives, anchors or controls. Nothing recorded a killed hypothesis. The B1 drug-by-training
contrast, the B3 prediction that was retracted, fourteen compounds labelled non-enhancing, and the
L4 plasticity window demoted this week all died into prose inside markdown files, where no program
can reach them. So a kill was permanent by default, not by evidence: when a capability landed that
would have changed a verdict, nothing re-opened it, because nothing knew the verdict existed.

That is not hypothetical. On 2026-09-15 the B2 harness scored the L4 window on parent-compound
SMILES and recorded psilocybin as a false negative. PERSEUS had resolved psilocybin to psilocin
through `PRODRUG_TO_ACTIVE` for months. The correction was found by hand, and agreement moved
0.33 -> 0.50. A machine-checkable record of what each verdict DEPENDED ON would have caught it.

THE DISCRIMINATION PROBLEM. An archive of failures is only useful if it separates two things that
look identical in a results table:

  - nature refused the mechanism            -> the hypothesis is closed, and reviving it is crankery
  - we asked the wrong question             -> the hypothesis is untested, and the kill was ours

This is the whole difficulty. A discard pile mined without that distinction produces pseudoscience,
and a discard pile never mined at all throws away the cases where the experiment, not the world,
was at fault. `FAILURE_MODES` below is that separation made explicit and machine-readable, and it is
deliberately conservative: a mode is revivable only when the recorded evidence CANNOT distinguish
"no effect" from "not measurable by the test that was run".

The first thing it found was in this project's own primary ledger. Of fourteen compounds labelled
`enhances_healthy_young = 0`, only seven have a confidence interval whose upper bound excludes a
target-sized effect (g = 0.25, the project's own replication-power target). Three have intervals
that ADMIT one (dextroamphetamine reaches +0.47), and four carry no interval at all. A binary label
cannot hold that difference, so the ledger asserts fourteen refutations it has evidence for seven of.

WHAT THIS MODULE DOES NOT DO. It never revives anything. It reports which dead hypotheses have a
keystone predicate that is now satisfied, and a human decides. This follows `validation/ledger_guard`,
which reports violations rather than editing data, for the same reason: an automatic revival is an
automatic claim, and claims in this project are made by people.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd


@dataclass(frozen=True)
class FailureMode:
    """One way a hypothesis can die, and whether death is the end of it."""
    name: str
    revivable: bool
    definition: str


#: THE DISCRIMINATION LAYER. Ordered nature-refused first, then ours-to-fix.
#:
#: `revivable=False` is a strong claim and is reserved for evidence that positively excludes the
#: effect, not merely evidence that failed to find it. The distinction is the entire point: "the
#: interval excludes a target-sized effect" closes a question, "the study found p > 0.05" does not.
FAILURE_MODES: dict[str, FailureMode] = {m.name: m for m in [
    # ---- nature refused. Do not re-open without NEW primary evidence. --------------------------
    FailureMode("physics_refused", False,
                "violates a physical, thermodynamic or biophysical constraint (insoluble at the "
                "required concentration, cannot cross the barrier at any tolerable dose, the "
                "molecule does not exist in the claimed form)"),
    FailureMode("measured_null", False,
                "adequately powered test in the RIGHT population whose interval EXCLUDES an effect "
                "of practical size. Not 'p > 0.05'. The interval must do the work."),
    FailureMode("measured_harm", False,
                "the interval excludes zero in the harmful direction"),
    FailureMode("mechanism_contradicted", False,
                "the asserted causal mechanism was directly tested and found false (knockout "
                "retains the effect, the active species is something else, the target is absent)"),

    # ---- we asked wrong. These are the stepping stones. -----------------------------------------
    FailureMode("wrong_population", True,
                "tested in patients or an impaired sample when the claim is about healthy adults, "
                "or vice versa"),
    FailureMode("wrong_endpoint", True,
                "measured something other than the claim (an on-drug readout for a durability "
                "claim, a biomarker for a cognitive claim, one assay family for another)"),
    FailureMode("wrong_dose_regime", True,
                "a monotonic or single-dose design applied to a compound with a non-monotonic "
                "dose-response, so the effective window was never sampled"),
    FailureMode("wrong_formulation", True,
                "delivery failed rather than the mechanism (no CNS exposure, prodrug not cleaved, "
                "peak-concentration toxicity forced a sub-therapeutic dose)"),
    FailureMode("underpowered", True,
                "the interval ADMITS an effect of practical size; the test could not have detected "
                "what was being looked for"),
    FailureMode("unknown_precision", True,
                "no interval was recorded at all, so a null cannot be distinguished from silence"),
    FailureMode("instrument_blind", True,
                "OUR measurement could not see it: the model class is structurally blind to the "
                "mechanism, the screen has no opinion about this chemotype, the corpus does not "
                "index this kind of result"),
    FailureMode("confounded", True,
                "a design confound makes the result uninterpretable in either direction"),
    FailureMode("provenance_failed", True,
                "the supporting source could not be verified, so the claim is unresolved rather "
                "than refuted"),
    FailureMode("superseded", True,
                "replaced by a better-specified version of the same idea"),
]}

REVIVABLE_MODES = {k for k, v in FAILURE_MODES.items() if v.revivable}
CLOSED_MODES = {k for k, v in FAILURE_MODES.items() if not v.revivable}

VALID_STATUS = {"DEAD", "REVIVABLE", "REVIVED", "PERMANENTLY_CLOSED"}

REQUIRED_COLUMNS = (
    "hypothesis_id",        # stable key, never reused
    "claim",                # what was asserted, in one sentence
    "domain",               # which arc/lane it belonged to
    "verdict",              # KILLED | RETRACTED | REFUTED | DEMOTED | ABSTAINED
    "died_on",              # ISO date
    "killed_by",            # the artifact that did it: script, report, or commit
    "evidence",             # the MEASURED number that killed it, not a summary
    "failure_mode",         # a key of FAILURE_MODES -- the discrimination decision
    "keystone",             # prose: what was missing
    "keystone_predicate",   # MACHINE-CHECKABLE condition for re-opening (see predicates.py)
    "revival_test",         # what to run if the predicate fires
    "status",
)


def _blank(value) -> bool:
    """Is this cell empty?

    pandas reads an empty CSV cell as NaN, and `str(nan)` is the four-character string "nan", which
    is truthy. Every emptiness check in this module routes through here because the first version
    did not: it reported a permanently-closed entry as carrying a revival predicate, and, far worse,
    it would have PASSED a revivable entry that named no keystone at all. That is the one rule the
    archive exists to enforce, and a truthy NaN had silently disabled it.
    """
    if value is None:
        return True
    if isinstance(value, float) and pd.isna(value):
        return True
    return str(value).strip() == ""


@dataclass
class ArchiveViolation:
    hypothesis_id: str
    severity: str           # "error" | "warn"
    rule: str
    detail: str

    def __str__(self) -> str:      # pragma: no cover - display only
        return f"[{self.severity.upper()}] {self.hypothesis_id}: {self.rule} -- {self.detail}"


def validate_archive(df: pd.DataFrame) -> list[ArchiveViolation]:
    """Check the archive against its own contract. Reports; never edits.

    The load-bearing rule is the last one: a hypothesis whose failure mode is REVIVABLE must carry a
    keystone predicate. Without it the entry is a tombstone wearing a stepping stone's clothes, and
    the retro-validation engine cannot ever re-open it. That is precisely the failure this archive
    was built to end, so it is an error rather than a warning.
    """
    v: list[ArchiveViolation] = []

    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        return [ArchiveViolation("<schema>", "error", "missing_columns", str(missing))]

    dupes = df["hypothesis_id"][df["hypothesis_id"].duplicated()].tolist()
    for h in sorted(set(dupes)):
        v.append(ArchiveViolation(h, "error", "duplicate_id",
                                  "hypothesis_id must be unique and is never reused"))

    for _, r in df.iterrows():
        h = str(r.get("hypothesis_id", "<missing>"))
        mode = str(r.get("failure_mode", "")).strip()

        if mode not in FAILURE_MODES:
            v.append(ArchiveViolation(h, "error", "bad_failure_mode",
                                      f"{mode!r} is not one of {sorted(FAILURE_MODES)}"))
            continue

        status = str(r.get("status", "")).strip()
        if status not in VALID_STATUS:
            v.append(ArchiveViolation(h, "error", "bad_status",
                                      f"{status!r} not in {sorted(VALID_STATUS)}"))

        if _blank(r.get("evidence")):
            v.append(ArchiveViolation(h, "error", "no_evidence",
                                      "a kill must record the MEASURED result that caused it; "
                                      "'it did not work' is not an archive entry"))
        if _blank(r.get("killed_by")):
            v.append(ArchiveViolation(h, "error", "no_provenance",
                                      "record the artifact that killed it (script, report, commit) "
                                      "so the verdict can be re-run"))

        revivable = FAILURE_MODES[mode].revivable
        predicate = ("" if _blank(r.get("keystone_predicate"))
                     else str(r.get("keystone_predicate")).strip())

        if revivable and not predicate:
            v.append(ArchiveViolation(
                h, "error", "revivable_without_keystone",
                f"failure_mode={mode!r} is revivable, so this entry MUST name a machine-checkable "
                "keystone_predicate. An entry that can come back but cannot say what would bring "
                "it back is a tombstone, and the archive exists to stop producing those."))
        if revivable and _blank(r.get("revival_test")):
            v.append(ArchiveViolation(h, "warn", "no_revival_test",
                                      "name what to RUN when the keystone fires, or the flag will "
                                      "arrive with nothing to do about it"))

        if not revivable:
            if status != "PERMANENTLY_CLOSED":
                v.append(ArchiveViolation(
                    h, "warn", "closed_mode_not_marked_closed",
                    f"failure_mode={mode!r} positively excludes the effect; status should be "
                    "PERMANENTLY_CLOSED so no sweep re-opens it"))
            if predicate:
                v.append(ArchiveViolation(
                    h, "warn", "closed_mode_carries_predicate",
                    "a permanently closed entry naming a revival condition invites exactly the "
                    "crank re-opening the discrimination layer is meant to prevent"))

    return v


def assert_archive_valid(df: pd.DataFrame) -> None:
    """Raise on any error-severity violation."""
    errs = [x for x in validate_archive(df) if x.severity == "error"]
    if errs:
        raise ValueError("stepping-stone archive failed validation:\n"
                         + "\n".join(str(e) for e in errs))


# -------------------------------------------------------------------------------------------------
# The discrimination helper that produced this archive's first finding.
# -------------------------------------------------------------------------------------------------

#: The smallest effect this project treats as practically meaningful. A "null" whose interval
#: reaches past this has NOT excluded the effect being looked for.
#:
#: CORRECTED 2026-09-20. This was introduced as 0.25 and described as "the project's replication-
#: power target", citing scripts/126. That was wrong on both counts. scripts/126 defines 0.25 as
#: META_AUGMENTATION_CEILING, the d-cycloserine IPD meta-analytic PEAK across 21 RCTs (n = 1047),
#: which is an empirical upper bound on what drug augmentation achieves, not a target anyone set.
#: The project's actual meaningfulness floor is MEANINGFUL_G = 0.20 in scripts/121, which predates
#: this module. Two constants for one concept is how a ledger ends up disagreeing with itself, so
#: there is now one, defined here and imported by scripts/121.
#:
#: The direction of the correction is worth stating: LOWERING the bar means FEWER nulls count as
#: refuted, because excluding a smaller effect is a harder thing for an interval to do.
TARGET_EFFECT_G = 0.20


def classify_null(ci_lo, ci_hi, target: float = TARGET_EFFECT_G) -> str:
    """Is a non-significant result EVIDENCE OF NO EFFECT, or an absence of evidence?

    Returns a `FAILURE_MODES` key. This is the single most reusable piece of the discrimination
    layer, because the mistake it prevents is the commonest one in the whole field: reading
    'the interval included zero' as 'the compound does nothing'.

      no interval            -> unknown_precision  (a null cannot be distinguished from silence)
      interval reaches target-> underpowered       (the data ADMIT the effect being sought)
      interval excludes it   -> measured_null      (the question is closed)
      excludes zero, harmful -> measured_harm

    `target` is the smallest effect the project cares about, not a statistical threshold. A result
    can be highly significant and still be a measured_null for a target above its interval.
    """
    if pd.isna(ci_lo) or pd.isna(ci_hi):
        return "unknown_precision"
    lo, hi = float(ci_lo), float(ci_hi)
    if hi < 0:
        return "measured_harm"
    if hi < target:
        return "measured_null"
    return "underpowered"


@dataclass
class ArchiveSummary:
    total: int = 0
    revivable: int = 0
    closed: int = 0
    by_mode: dict = field(default_factory=dict)


def summarise(df: pd.DataFrame) -> ArchiveSummary:
    s = ArchiveSummary(total=len(df))
    s.by_mode = dict(df["failure_mode"].value_counts()) if len(df) else {}
    s.revivable = int(df["failure_mode"].isin(REVIVABLE_MODES).sum()) if len(df) else 0
    s.closed = int(df["failure_mode"].isin(CLOSED_MODES).sum()) if len(df) else 0
    return s
