"""TxGNN zero-shot indication scoring (STUB).

Per the v2 research doc §3 Class C:
    Source: Huang, Chandak, Wang et al. "A foundation model for clinician-centered
            drug repurposing." Nature Medicine 2024 Dec;30(12):3601–3613.
            DOI: 10.1038/s41591-024-03233-x
    GitHub: mims-harvard/TxGNN
    Training: PrimeKG (129K nodes, 4M edges, 17,080 disease coverage)
    Reports: +49.2% indication, +35.1% contraindication vs 8 baselines (zero-shot)
    VRAM: <2 GB at inference

Install:
    pip install git+https://github.com/mims-harvard/TxGNN.git
    (Or follow repo's setup instructions; depends on PyG.)

Usage (planned):
    from mammal_repurposing.cluster_c.txgnn import score_compounds_against_anchor
    df = score_compounds_against_anchor(
        compound_drugbank_ids=[...],
        anchor=CognitionAnchor.from_config(),
    )
    # Returns: compound_name, indication_score, contraindication_score per row
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from mammal_repurposing.cluster_c.cognition_anchor import CognitionAnchor

logger = logging.getLogger(__name__)


def score_compounds_against_anchor(
    compound_drugbank_ids: list[str],
    *,
    anchor: "CognitionAnchor",
    model_checkpoint: Path | str | None = None,
) -> pd.DataFrame:
    """Run TxGNN against each compound for each disease in the anchor.

    Returns DataFrame with columns:
        compound_drugbank_id | compound_name |
        indication_score (mean across anchor diseases, weighted) |
        contraindication_score (mean across anchor diseases, weighted) |
        per_disease_scores (dict)
    """
    try:
        import txgnn  # noqa: F401, PLC0415
    except ImportError as e:
        raise ImportError(
            "txgnn not installed. Run "
            "`pip install git+https://github.com/mims-harvard/TxGNN.git`. "
            "Note: requires PyG (torch-geometric) which can be tricky on Windows."
        ) from e

    raise NotImplementedError(
        "TxGNN runner is a stub. Load the pretrained checkpoint, query "
        "drug↔disease scores for each (compound, anchor_disease) pair, "
        "aggregate weighted-mean per anchor.disease_weights."
    )
