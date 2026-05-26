"""§7.5 Pocket-conditioned Boltz-2 — curated centroid DB + geometric classifier.

Per research/4-tier/Pocket-Conditioned-Boltz2.md.

MVP scope (this v1):
  pocket_database.py  - YAML loader + RCSB PDB fetch + biopython centroid extractor
  pocket_classifier.py - geometric assignment of Boltz pose centroid → pocket_class

Deferred to Sprint 2 (research doc §1):
  p2rank, pocketminer, cryptobench  - detector ensemble (Java + PyTorch deps)
  boltz2_cofold_oracle             - second-opinion runner on disagreement pairs
"""

from .pocket_database import (
    PocketSpec,
    PocketDatabase,
    load_pocket_database,
    fetch_pdb,
    compute_centroid_for_pocket,
)
from .pocket_classifier import (
    classify_pose,
    PocketClassification,
    RANK_MULTIPLIER,
)

__all__ = [
    "PocketSpec",
    "PocketDatabase",
    "load_pocket_database",
    "fetch_pdb",
    "compute_centroid_for_pocket",
    "classify_pose",
    "PocketClassification",
    "RANK_MULTIPLIER",
]
