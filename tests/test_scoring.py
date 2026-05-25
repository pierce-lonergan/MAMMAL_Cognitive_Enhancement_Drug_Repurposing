"""Tests for the MAMMAL scoring path.

These tests are marked ``slow`` because they load the real MAMMAL DTI model
(~1.8 GB download on first run, ~5 s warm load). Run with::

    pytest tests/test_scoring.py -m slow

Skip with the default ``pytest`` invocation.
"""

from __future__ import annotations

import math

import pytest


@pytest.mark.slow
def test_model_loads_and_scores_reference_pair():
    """End-to-end smoke: load the released DTI head, score the HF reference pair,
    expect a finite float roughly in [3, 10] (pKd range observed in BindingDB).
    """
    from mammal_repurposing.scoring.dti import score_pair
    from mammal_repurposing.scoring.model_loader import load_dti_model

    model, tokenizer = load_dti_model()

    # Reference inputs from the HF model card README
    target_seq = "NLMKRCTRGFRKLGKCTTLEEEKCKTLYPRGQCTCSDSKMNTHSCDCKSC"
    drug_seq = "CC(=O)NCCC1=CNc2c1cc(OC)cc2"

    pkd = score_pair(model, tokenizer, target_seq, drug_seq)
    assert math.isfinite(pkd), f"Got non-finite pKd: {pkd}"
    assert 0.0 < pkd < 15.0, f"pKd {pkd} outside plausible BindingDB range"


@pytest.mark.slow
def test_batch_matches_single():
    """Scoring N pairs in a batch should match scoring each individually."""
    from mammal_repurposing.scoring.dti import score_batch, score_pair
    from mammal_repurposing.scoring.model_loader import load_dti_model

    model, tokenizer = load_dti_model()

    pairs = [
        ("NLMKRCTRGFRKLGKCTTLEEEKCKTLYPRGQCTCSDSKMNTHSCDCKSC",
         "CC(=O)NCCC1=CNc2c1cc(OC)cc2"),
        ("MNLAAAMNLAAAMNLAAAMNLAAAMNLAAAMNLAAA",
         "CCO"),
    ]

    singles = [score_pair(model, tokenizer, t, d) for t, d in pairs]
    batched = score_batch(model, tokenizer, pairs)

    assert len(singles) == len(batched)
    for s, b in zip(singles, batched):
        assert abs(s - b) < 1e-3, f"single={s} batched={b} diverge"
