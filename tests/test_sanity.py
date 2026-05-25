"""Tests for the sanity-gate logic (no MAMMAL model required)."""

from __future__ import annotations

from mammal_repurposing.analysis.polypharm import compute_polypharm
from mammal_repurposing.analysis.sanity import (
    build_report,
    check_negative_controls,
    check_positive_controls,
    render_markdown,
)


def test_positive_controls_pass(mini_scores_pass):
    checks = check_positive_controls(mini_scores_pass)
    assert len(checks) > 0
    # Every check we ran should pass — donepezil/methylphenidate/pitolisant are
    # at pkd=8.0 against their targets, others ~4.0.
    for c in checks:
        if not c.found_compounds:
            continue  # target had no expected compound in our mini library
        assert c.passed, f"{c.target_gene} unexpectedly failed: {c.found_compounds}"


def test_positive_controls_fail(mini_scores_fail):
    checks = check_positive_controls(mini_scores_fail)
    failing = [c for c in checks if c.found_compounds and not c.passed]
    assert failing, "Expected at least one positive-control target to fail."


def test_negative_controls_flagged_when_top(mini_scores_fail, mini_compounds):
    hits = check_negative_controls(mini_scores_fail, mini_compounds)
    # loratadine has pkd=8.5 against every target in the fail fixture; should be flagged.
    assert hits
    assert all(h.compound_name == "loratadine" for h in hits)


def test_negative_controls_clean(mini_scores_pass, mini_compounds):
    hits = check_negative_controls(mini_scores_pass, mini_compounds)
    assert hits == []


def test_build_report_passes(mini_scores_pass, mini_compounds):
    report = build_report(mini_scores_pass, mini_compounds)
    assert report.passed
    assert report.n_targets_pass > 0
    assert report.total_pairs == len(mini_scores_pass)


def test_build_report_fails(mini_scores_fail, mini_compounds):
    report = build_report(mini_scores_fail, mini_compounds)
    assert not report.passed
    assert report.negative_hits


def test_render_markdown_includes_outcome(mini_scores_pass, mini_compounds):
    report = build_report(mini_scores_pass, mini_compounds)
    polypharm = compute_polypharm(mini_scores_pass)
    md = render_markdown(report, polypharm)
    assert "PASS" in md
    assert "Positive-Control Checks" in md
    assert "pKd Distribution Summary" in md


def test_render_markdown_fail_path(mini_scores_fail, mini_compounds):
    report = build_report(mini_scores_fail, mini_compounds)
    md = render_markdown(report, compute_polypharm(mini_scores_fail))
    assert "FAIL" in md
    assert "Negative-Control Hits" in md


def test_polypharm_basic(mini_scores_pass):
    pp = compute_polypharm(mini_scores_pass, threshold=6.0)
    # Only the 3 positive controls cross 6.0 (at pKd 8.0 against their targets).
    assert len(pp) == 3
    assert (pp["n_hits"] == 1).all()
    assert (pp["max_pkd"] == 8.0).all()


def test_polypharm_no_hits(mini_scores_fail):
    pp = compute_polypharm(mini_scores_fail, threshold=10.0)  # nothing above 10
    assert pp.empty
