"""ClinTox wrappers (toxicity + FDA-approval probability) around MAMMAL MoleculeNet heads."""

from __future__ import annotations

from mammal_repurposing.scoring.molnet import MolnetResult, score_task_batch


def score_clintox_tox_batch(
    smiles_list: list[str],
    *,
    device: str | None = None,
) -> list[MolnetResult]:
    """Predict clinical toxicity per SMILES. ``score`` is P(toxic-in-trials)."""
    return score_task_batch("TOXICITY", smiles_list, device=device)


def score_clintox_fda_batch(
    smiles_list: list[str],
    *,
    device: str | None = None,
) -> list[MolnetResult]:
    """Predict FDA-approval-similarity per SMILES. ``score`` is P(looks-like-approved-drug)."""
    return score_task_batch("FDA_APPR", smiles_list, device=device)
