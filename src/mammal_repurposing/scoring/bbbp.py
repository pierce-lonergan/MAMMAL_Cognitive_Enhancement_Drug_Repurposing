"""BBBP (blood-brain-barrier permeability) wrapper around the MAMMAL MoleculeNet head."""

from __future__ import annotations

from mammal_repurposing.scoring.molnet import MolnetResult, score_task_batch


def score_bbbp_batch(smiles_list: list[str], *, device: str | None = None) -> list[MolnetResult]:
    """Predict BBB permeability per SMILES. ``score`` is P(BBB-permeable)."""
    return score_task_batch("BBBP", smiles_list, device=device)
