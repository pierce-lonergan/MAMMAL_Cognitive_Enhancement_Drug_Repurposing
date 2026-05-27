"""Shared pytest fixtures."""

from __future__ import annotations

import sys
from pathlib import Path

# Make src/ importable for all test modules without requiring pip install -e .
_SRC = Path(__file__).resolve().parents[1] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import pandas as pd
import pytest


@pytest.fixture
def mini_targets() -> pd.DataFrame:
    """3-target panel for unit tests."""
    return pd.DataFrame([
        {"uniprot": "P22303", "gene": "ACHE",
         "mechanism_class": "cholinergic", "cognitive_domain": "memory",
         "notes": "donepezil target", "sequence": "MNLAAA" * 20, "seq_length": 120,
         "gene_name": "ACHE", "ensembl_gene_id": "ENSG00000087085"},
        {"uniprot": "Q01959", "gene": "SLC6A3",
         "mechanism_class": "dopaminergic", "cognitive_domain": "attention",
         "notes": "DAT", "sequence": "MSKSAA" * 20, "seq_length": 120,
         "gene_name": "SLC6A3", "ensembl_gene_id": "ENSG00000142319"},
        {"uniprot": "Q9Y5N1", "gene": "HRH3",
         "mechanism_class": "histaminergic", "cognitive_domain": "arousal",
         "notes": "pitolisant", "sequence": "MERPAA" * 20, "seq_length": 120,
         "gene_name": "HRH3", "ensembl_gene_id": "ENSG00000101180"},
    ])


@pytest.fixture
def mini_compounds() -> pd.DataFrame:
    """5-compound library for unit tests: 3 positive ctrls, 1 named, 1 negative."""
    return pd.DataFrame([
        {"name": "donepezil", "smiles": "O=C1CC2=CC=C(OC)C=C2CC1CC1CCN(CC2=CC=CC=C2)CC1",
         "smiles_kind": "canonical", "cid": 3152, "source": "seed",
         "mechanism_class": "cholinergic", "evidence_tier": "positive_control",
         "expected_top_target": "P22303", "notes": "AChE inhibitor"},
        {"name": "methylphenidate", "smiles": "COC(=O)C(C1CCCCN1)C1=CC=CC=C1",
         "smiles_kind": "canonical", "cid": 4158, "source": "seed",
         "mechanism_class": "dopaminergic", "evidence_tier": "positive_control",
         "expected_top_target": "Q01959", "notes": "DAT inhibitor"},
        {"name": "pitolisant", "smiles": "C(CCN1CCCCC1)(CCC1=CC=C(OCCC2=CC=CC=C2)C=C1)Cl",
         "smiles_kind": "canonical", "cid": 9948102, "source": "seed",
         "mechanism_class": "histaminergic", "evidence_tier": "positive_control",
         "expected_top_target": "Q9Y5N1", "notes": "H3 inverse agonist"},
        {"name": "aniracetam", "smiles": "O=C1N(CCO)CCC1C(=O)C1=CC=C(OC)C=C1",
         "smiles_kind": "canonical", "cid": 2196, "source": "seed",
         "mechanism_class": "glutamatergic_ampa", "evidence_tier": "named_in_research",
         "expected_top_target": "", "notes": "racetam"},
        {"name": "loratadine", "smiles": "CCOC(=O)N1CCC(=C2C3=NC=CC=C3CCC3=CC(Cl)=CC=C23)CC1",
         "smiles_kind": "canonical", "cid": 3957, "source": "negative_control",
         "mechanism_class": "peripheral_antihistamine", "evidence_tier": "negative_control",
         "expected_top_target": "", "notes": "peripheral H1"},
    ])


@pytest.fixture
def mini_scores_pass(mini_targets, mini_compounds) -> pd.DataFrame:
    """Score grid where positive controls rank at top of their targets — should PASS."""
    rows = []
    for _, t in mini_targets.iterrows():
        for _, c in mini_compounds.iterrows():
            # Give expected pair a pKd of 8.0, everyone else 4.0 (clear winner).
            if c["expected_top_target"] == t["uniprot"]:
                pkd = 8.0
            elif c["name"] == "loratadine":
                pkd = 3.5
            else:
                pkd = 4.0 + (hash(t["uniprot"] + c["name"]) % 100) / 200.0  # 4.0-4.5 noise
            rows.append({
                "target_uniprot": t["uniprot"],
                "target_gene": t["gene"],
                "compound_name": c["name"],
                "compound_smiles": c["smiles"],
                "predicted_pkd": pkd,
                "model_version": "test",
                "scored_at": "2026-05-24T00:00:00+00:00",
            })
    return pd.DataFrame(rows)


@pytest.fixture
def mini_scores_fail(mini_targets, mini_compounds) -> pd.DataFrame:
    """Score grid where positive controls rank LOW — should FAIL the gate."""
    rows = []
    for _, t in mini_targets.iterrows():
        for _, c in mini_compounds.iterrows():
            if c["expected_top_target"] == t["uniprot"]:
                pkd = 3.0  # buried at the bottom
            elif c["name"] == "loratadine":
                pkd = 8.5  # negative control near the top — should trigger neg flag too
            else:
                pkd = 5.0
            rows.append({
                "target_uniprot": t["uniprot"],
                "target_gene": t["gene"],
                "compound_name": c["name"],
                "compound_smiles": c["smiles"],
                "predicted_pkd": pkd,
                "model_version": "test",
                "scored_at": "2026-05-24T00:00:00+00:00",
            })
    return pd.DataFrame(rows)
