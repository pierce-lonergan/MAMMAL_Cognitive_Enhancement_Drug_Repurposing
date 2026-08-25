"""The Tanimoto ranker, and specifically that it no longer reads its own label.

Until 2026-08-24 `tanimoto_score` was the maximum similarity to the target's
ChEMBL actives at pChEMBL >= 8.0, and the query compound was a member of that
set whenever its own affinity cleared the threshold. 143 of 289 rows in the
Gap-4 evaluation set scored exactly 1.000. Every case below is a way that leak
can come back.
"""

from __future__ import annotations

import math

import pytest

from mammal_repurposing.cluster_a.tanimoto_ranker import (
    TanimotoRankerConfig,
    _skeleton_key,
    score_library_against_target,
    score_library_against_target_audited,
)

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
ASPIRIN_NA = "CC(=O)Oc1ccccc1C(=O)O.[Na+]"
S_AMPHETAMINE = "C[C@H](N)Cc1ccccc1"
R_AMPHETAMINE = "C[C@@H](N)Cc1ccccc1"
IBUPROFEN = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
CAFFEINE = "Cn1cnc2c1c(=O)n(C)c(=O)n2C"


def score(lib, actives, **cfg):
    return score_library_against_target(lib, actives, TanimotoRankerConfig(**cfg))


# -- the leak itself ---------------------------------------------------------

def test_the_query_does_not_score_against_its_own_record() -> None:
    assert score([ASPIRIN], [ASPIRIN, IBUPROFEN])[0] < 1.0


def test_a_salt_of_the_query_is_still_the_query() -> None:
    """Fingerprint equality would MISS this: aspirin and its sodium salt score
    0.96 against each other, not 1.000, so an exclusion keyed on identical
    fingerprints leaves the self-read in place at 0.96."""
    with_salt = score([ASPIRIN], [ASPIRIN_NA, IBUPROFEN])[0]
    without = score([ASPIRIN], [IBUPROFEN])[0]
    assert with_salt == pytest.approx(without)


def test_the_other_enantiomer_is_still_the_query() -> None:
    """ECFP4 is built here without useChirality, so enantiomers have Tanimoto
    exactly 1.000 and this feature cannot tell them apart. Excluding one and
    keeping the other would leave the self-read intact."""
    assert score([S_AMPHETAMINE], [R_AMPHETAMINE, IBUPROFEN])[0] < 1.0


def test_a_compound_that_is_genuinely_absent_is_not_penalised() -> None:
    """The exclusion must not quietly remove real neighbours."""
    _, dropped = score_library_against_target_audited([CAFFEINE], [ASPIRIN, IBUPROFEN])
    assert dropped == [0]


def test_actives_consisting_only_of_the_query_give_nan_not_one() -> None:
    """The purest form of the leak: nothing left to compare against, and the
    old code would have returned 1.000."""
    assert math.isnan(score([ASPIRIN], [ASPIRIN])[0])


def test_exclusion_is_per_compound_not_global() -> None:
    """Excluding aspirin for the aspirin row must not remove it from the
    ibuprofen row, where it is a legitimate neighbour."""
    actives = [ASPIRIN, CAFFEINE]
    a_self = score([ASPIRIN], actives)[0]
    ibu = score([IBUPROFEN], actives)[0]
    ibu_without_aspirin = score([IBUPROFEN], [CAFFEINE])[0]
    assert a_self < 1.0
    assert ibu > ibu_without_aspirin


def test_audit_counts_what_was_dropped() -> None:
    scores, dropped = score_library_against_target_audited(
        [ASPIRIN, CAFFEINE], [ASPIRIN, ASPIRIN_NA, IBUPROFEN])
    assert dropped == [2, 0]
    assert scores[0] < 1.0


# -- the escape hatch, which must be loud and explicit -----------------------

def test_the_old_behaviour_is_reachable_only_on_request() -> None:
    assert score([ASPIRIN], [ASPIRIN, IBUPROFEN], exclude_self=False) == [1.0]


def test_exclusion_is_the_default() -> None:
    assert TanimotoRankerConfig().exclude_self is True


# -- the key ----------------------------------------------------------------

def test_skeleton_key_collapses_salt_and_stereochemistry() -> None:
    k = _skeleton_key(ASPIRIN)
    assert k is not None and len(k) == 14
    assert _skeleton_key(ASPIRIN_NA) == k
    assert _skeleton_key(S_AMPHETAMINE) == _skeleton_key(R_AMPHETAMINE)


def test_skeleton_key_separates_different_compounds() -> None:
    assert _skeleton_key(ASPIRIN) != _skeleton_key(IBUPROFEN)


def test_unparseable_input_is_none_and_never_matches() -> None:
    assert _skeleton_key("not a molecule") is None
    assert _skeleton_key("") is None
    assert _skeleton_key(None) is None  # type: ignore[arg-type]
    # An unparseable ACTIVE must not be excluded from every query by matching
    # None against None.
    _, dropped = score_library_against_target_audited(
        ["not a molecule either"], [ASPIRIN])
    assert dropped == [0]


def test_an_unparseable_query_still_scores_nan() -> None:
    assert math.isnan(score(["not a molecule"], [ASPIRIN])[0])


# -- the wrapper keeps its shape --------------------------------------------

def test_plain_scorer_returns_only_scores() -> None:
    out = score_library_against_target([ASPIRIN, IBUPROFEN], [CAFFEINE])
    assert isinstance(out, list) and len(out) == 2
    assert all(isinstance(v, float) for v in out)
