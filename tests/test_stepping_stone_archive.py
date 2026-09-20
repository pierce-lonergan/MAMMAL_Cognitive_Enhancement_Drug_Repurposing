"""The stepping-stone archive must enforce the discrimination it exists to make.

These tests fail on pre-archive code (nothing to import) and pass after. The ones that matter are
the two that pin actual defects found while building it:

  - `test_blank_cell_is_blank_even_as_nan` pins the NaN bug. pandas reads an empty CSV cell as NaN
    and `str(nan)` is truthy, so the first version of `validate_archive` would have PASSED a
    revivable entry carrying no keystone at all. That is the single rule the archive exists to
    enforce, silently disabled.
  - `test_durability_predicate_rejects_the_known_confounded_precedent` pins the loose-predicate bug.
    The first `durable_healthy_rows` matched "positive" as a substring, swept in the authors' own
    hedges, ignored replication, and flagged a revival on the strength of the three studies this
    project had already examined and excluded.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from mammal_repurposing.archive.predicates import (
    UnknownPredicate, evaluate, parse, registered,
)
from mammal_repurposing.archive.stepping_stone import (
    FAILURE_MODES, REQUIRED_COLUMNS, TARGET_EFFECT_G, _blank, classify_null, validate_archive,
)

ARCHIVE = Path(__file__).resolve().parents[1] / "data" / "raw" / "stepping_stone_archive.csv"


def _row(**kw):
    r = {"hypothesis_id": "H-x", "claim": "c", "domain": "d", "verdict": "KILLED",
         "died_on": "2026-01-01", "killed_by": "scripts/1.py", "evidence": "g = 0.0, CI [-1, 1]",
         "failure_mode": "underpowered", "keystone": "k",
         "keystone_predicate": "dose_axis_exists()", "revival_test": "re-run",
         "status": "DEAD"}
    r.update(kw)
    return pd.DataFrame([r])


# --- the discrimination layer --------------------------------------------------------------------

def test_the_only_thing_that_closes_a_question_is_an_interval_that_excludes_the_target():
    """'p > 0.05' is not a refutation. The interval has to do the work."""
    assert classify_null(-0.05, 0.06) == "measured_null"          # excludes a target-sized effect
    assert classify_null(-0.06, 0.47) == "underpowered"           # ADMITS one: dextroamphetamine
    assert classify_null(None, None) == "unknown_precision"       # cannot tell null from silence
    assert classify_null(-0.62, -0.06) == "measured_harm"
    # exactly at the target is NOT closed: the data still admit the effect being sought
    assert classify_null(-0.1, TARGET_EFFECT_G) == "underpowered"


def test_closed_modes_are_only_those_that_positively_exclude_the_effect():
    closed = {k for k, v in FAILURE_MODES.items() if not v.revivable}
    assert closed == {"physics_refused", "measured_null", "measured_harm", "mechanism_contradicted"}
    # everything that is about OUR test rather than the world must stay revivable
    for m in ("underpowered", "unknown_precision", "instrument_blind", "wrong_population",
              "wrong_endpoint", "wrong_dose_regime", "wrong_formulation", "confounded",
              "provenance_failed"):
        assert FAILURE_MODES[m].revivable, m


# --- the NaN defect ------------------------------------------------------------------------------

def test_blank_cell_is_blank_even_as_nan():
    """`str(float('nan'))` is 'nan', which is truthy. That silently disabled the archive's one rule."""
    assert _blank(float("nan")) and _blank(None) and _blank("") and _blank("   ")
    assert not _blank("dose_axis_exists()")


def test_revivable_entry_with_no_keystone_is_an_error_even_when_read_from_csv(tmp_path):
    """The regression: round-trip through CSV so the empty cell arrives as NaN, as it does live."""
    p = tmp_path / "a.csv"
    _row(keystone_predicate="", revival_test="").to_csv(p, index=False)
    v = validate_archive(pd.read_csv(p))
    assert any(x.rule == "revivable_without_keystone" and x.severity == "error" for x in v), \
        "a revivable entry naming no keystone must be an ERROR; NaN must not read as a predicate"


def test_permanently_closed_entry_should_not_advertise_a_revival_route():
    v = validate_archive(_row(failure_mode="measured_null", status="PERMANENTLY_CLOSED"))
    assert any(x.rule == "closed_mode_carries_predicate" for x in v)


def test_a_kill_must_record_the_number_that_killed_it():
    assert any(x.rule == "no_evidence" for x in validate_archive(_row(evidence="")))
    assert any(x.rule == "no_provenance" for x in validate_archive(_row(killed_by="")))


def test_ids_are_never_reused():
    df = pd.concat([_row(), _row()], ignore_index=True)
    assert any(x.rule == "duplicate_id" for x in validate_archive(df))


# --- predicates ----------------------------------------------------------------------------------

def test_an_unimplemented_predicate_raises_rather_than_reading_as_still_dead():
    """A predicate the engine cannot evaluate must never silently mean False."""
    with pytest.raises(UnknownPredicate):
        evaluate("no_such_predicate()")
    with pytest.raises(UnknownPredicate):
        evaluate("not a call at all")
    with pytest.raises(UnknownPredicate):
        evaluate("ci_recorded()")            # wrong arity


def test_unresolvable_predicate_is_distinct_from_false():
    """cluster_d_v2g is deliberately unresolvable: returning False would assert an unmade finding."""
    with pytest.raises(UnknownPredicate, match="no source"):
        evaluate("cluster_d_v2g_attribution_measured()")


def test_parse_handles_zero_and_multi_argument_forms():
    assert parse("dose_axis_exists()") == ("dose_axis_exists", [])
    assert parse("ledger_rows_at_least(x.csv, 5)") == ("ledger_rows_at_least", ["x.csv", "5"])


def test_the_psilocybin_dependency_is_detectable():
    """The worked example: the exact staleness that was caught by hand on 2026-09-15."""
    ok, detail = evaluate("prodrug_resolution_covers(psilocybin)")
    assert ok and "PRESENT" in detail
    ok2, detail2 = evaluate("prodrug_resolution_covers(ibogaine)")
    assert not ok2 and "ABSENT" in detail2, "ibogaine is still a live false positive"


@pytest.mark.skipif(not ARCHIVE.exists(), reason="archive not built")
def test_durability_predicate_rejects_the_known_confounded_precedent():
    """The loose-predicate regression.

    Three off-drug healthy rows exist: Rokem & Silver 2013 (n=8, unreplicated, placebo arm also
    retained), Shellshear 2015 (effect unobtainable) and Chamoun 2017 (n=9, preliminary). A
    substring match on "positive" promoted all three and manufactured a revival. The predicate must
    agree with the project's headline that ZERO verified durable rows exist.
    """
    ok, detail = evaluate("durable_healthy_rows(1)")
    assert not ok, f"predicate must not fire on the known-confounded precedent: {detail}"
    assert "0 row(s) meeting the full durability standard" in detail
    # and it must show its work rather than just saying no
    assert "off-drug and healthy at all" in detail


@pytest.mark.skipif(not ARCHIVE.exists(), reason="archive not built")
def test_the_shipped_archive_satisfies_its_own_contract():
    df = pd.read_csv(ARCHIVE)
    errs = [x for x in validate_archive(df) if x.severity == "error"]
    assert not errs, "archive has contract breaches:\n" + "\n".join(str(e) for e in errs)
    assert set(REQUIRED_COLUMNS) <= set(df.columns)


@pytest.mark.skipif(not ARCHIVE.exists(), reason="archive not built")
def test_every_keystone_in_the_shipped_archive_parses_to_a_real_predicate():
    """A typo in a keystone would report as 'still dead' forever. Catch it at test time."""
    df = pd.read_csv(ARCHIVE)
    known = set(registered())
    for _, r in df.iterrows():
        if not FAILURE_MODES[r["failure_mode"]].revivable:
            continue
        name, _args = parse(str(r["keystone_predicate"]))
        assert name in known, f"{r['hypothesis_id']} names unknown predicate {name!r}"


@pytest.mark.skipif(not ARCHIVE.exists(), reason="archive not built")
def test_every_null_in_the_ledger_is_now_classifiable():
    """The resolving-power gain, pinned.

    ON 2026-09-20 MORNING the derived null set was 14 entries: 7 measured_null/harm, 3
    under-powered, and 4 `unknown_precision` -- rows carrying a point estimate and NO interval, where
    a true null cannot be distinguished from an undetected effect. The ledger asserted fourteen
    refutations and had evidence for seven.

    That afternoon the intervals were recovered from primary sources, and the recovery exposed two
    extraction defects rather than merely filling blanks. bacopa_monnieri's recorded g = 0.00 was
    read from the wrong cell of a network meta-analysis league table (BME vs anthocyanin, not BME vs
    placebo); the real estimate is 0.17 with an interval that ADMITS a target-sized effect. omega_3
    pooled cognitively normal older adults together with MCI participants, which fails this ledger's
    own clean-healthy rule, so it left the primary set entirely.

    The headline is `unknown_precision == 0`: every remaining null is now either genuinely refuted
    or honestly open. If future curation adds a row without an interval this fails, which is the
    point.
    """
    df = pd.read_csv(ARCHIVE)
    nulls = df[df["hypothesis_id"].str.startswith("H-null-")]
    modes = dict(nulls["failure_mode"].value_counts())
    assert modes.get("unknown_precision", 0) == 0, (
        "a null with no interval is an unclassifiable assertion; recover the interval or re-tier "
        f"the row. Classification now: {modes}")
    closed = modes.get("measured_null", 0) + modes.get("measured_harm", 0)
    assert closed + modes.get("underpowered", 0) == len(nulls),         f"every null must be classified as closed or open: {modes}"
    # the two defects the recovery exposed, pinned so a silent revert is visible
    assert not nulls["hypothesis_id"].eq("H-null-omega_3").any(),         "omega_3 pools MCI participants and must not be in the clean healthy-adult null set"
    bac = nulls[nulls["hypothesis_id"] == "H-null-bacopa_monnieri"]
    assert len(bac) == 1 and bac["failure_mode"].iloc[0] == "underpowered",         "bacopa_monnieri's interval admits a target-sized effect; it is open, not refuted"


def test_the_meaningfulness_threshold_is_defined_exactly_once():
    """Two constants for one concept is how a ledger ends up disagreeing with itself.

    Found 2026-09-20. `archive/stepping_stone.TARGET_EFFECT_G` was introduced as 0.25 and described
    as "the project's replication-power target", citing scripts/126. Both halves were wrong.
    scripts/126 defines 0.25 as META_AUGMENTATION_CEILING, the d-cycloserine IPD meta-analytic PEAK
    across 21 RCTs, which is an empirical upper bound on what drug augmentation achieves rather than
    a target anyone chose. Meanwhile scripts/121 had carried its own MEANINGFUL_G = 0.20 for the
    same concept since well before this module existed.

    The consequence was that the robustness report and the stepping-stone archive classified the
    same rows differently, silently. There is now one definition, here, imported by scripts/121.
    """
    import importlib.util
    from mammal_repurposing.archive.stepping_stone import TARGET_EFFECT_G

    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location(
        "_rob", root / "scripts" / "121_healthy_adult_robustness.py")
    rob = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(rob)

    assert rob.MEANINGFUL_G == TARGET_EFFECT_G, (
        f"scripts/121 uses {rob.MEANINGFUL_G} while the archive uses {TARGET_EFFECT_G}. One "
        "concept, one constant: import it rather than redeclaring it.")
    assert TARGET_EFFECT_G == 0.20, (
        "the meaningfulness floor is 0.20 (scripts/121, predating the archive). 0.25 is the "
        "d-cycloserine augmentation CEILING in scripts/126 and is a different quantity.")

    # scripts/121 must not have re-introduced a local definition. Checked on CODE lines only:
    # the comment recording this history mentions the old declaration, and the first version of
    # this test matched its own explanation.
    src = (root / "scripts" / "121_healthy_adult_robustness.py").read_text(encoding="utf-8")
    code = [ln for ln in src.split("\n") if not ln.lstrip().startswith("#")]
    offenders = [ln for ln in code if ln.lstrip().startswith("MEANINGFUL_G")
                 and "=" in ln and "import" not in ln]
    assert not offenders, (
        "scripts/121 has redeclared the threshold instead of importing it: " + str(offenders))
