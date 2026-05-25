"""Sanity tests for the config module — catches drift between POSITIVE_CONTROLS
and the targets seed CSV without needing the network."""

from __future__ import annotations

import pandas as pd

from mammal_repurposing.config import POSITIVE_CONTROLS, TARGETS_SEED_CSV


def test_positive_control_uniprots_are_in_target_panel():
    """Every UniProt in POSITIVE_CONTROLS must appear in targets_seed.csv."""
    seed = pd.read_csv(TARGETS_SEED_CSV)
    panel_uniprots = set(seed["uniprot"].tolist())
    for uniprot in POSITIVE_CONTROLS:
        assert uniprot in panel_uniprots, (
            f"POSITIVE_CONTROLS references {uniprot} but it's not in {TARGETS_SEED_CSV}. "
            f"Update one or the other."
        )


def test_targets_seed_has_no_duplicate_uniprots():
    seed = pd.read_csv(TARGETS_SEED_CSV)
    dup = seed[seed["uniprot"].duplicated()]
    assert dup.empty, f"Duplicate uniprots in seed: {dup['uniprot'].tolist()}"


def test_targets_seed_has_no_duplicate_genes():
    seed = pd.read_csv(TARGETS_SEED_CSV)
    dup = seed[seed["gene"].duplicated()]
    assert dup.empty, f"Duplicate genes in seed: {dup['gene'].tolist()}"
