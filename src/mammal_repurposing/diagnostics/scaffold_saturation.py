"""Diagnostic A — Murcko scaffold saturation analysis.

For each target, compute what fraction of the library's compounds share
the most-common generic Bemis-Murcko scaffold of high-affinity ChEMBL
binders. Routing thresholds (literature-anchored, Zhang ICLR 2025 +
Dablander 2023):

  > 60%  → rank-resolution loss in saturated cluster (Scenario 2; LoRA worth)
  25-60% → ambiguous → defer to Diagnostic D
  < 25%  → manifold mismatch (Scenario 1; cross-DTI ensemble worth)
"""

from __future__ import annotations

import logging
from collections import Counter
from dataclasses import dataclass

import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class ScaffoldResult:
    target_uniprot: str
    target_gene: str
    n_library: int
    n_chembl_actives: int
    top_chembl_scaffold_smiles: str
    top_chembl_scaffold_count: int
    library_in_top_scaffold_pct: float
    library_unique_scaffolds: int
    decision: str   # rank_resolution | manifold_mismatch | ambiguous


def _generic_murcko(smiles: str) -> str | None:
    """Bemis-Murcko generic scaffold (twice-applied per RDKit Discussion #6844)."""
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    s = MurckoScaffold.GetScaffoldForMol(mol)
    if s is None or s.GetNumAtoms() == 0:
        return ""   # acyclic compound; valid but distinct from "failed to parse"
    g = MurckoScaffold.MakeScaffoldGeneric(s)
    s2 = MurckoScaffold.GetScaffoldForMol(g) if g else None
    if s2 is None:
        return ""
    return Chem.MolToSmiles(s2, canonical=True)


def scaffold_overlap(
    target_uniprot: str,
    target_gene: str,
    library_smiles: list[str],
    chembl_active_smiles: list[str],
) -> ScaffoldResult:
    """Compute the headline Diagnostic A metric."""
    lib_scaffolds = [_generic_murcko(s) for s in library_smiles]
    chembl_scaffolds = [_generic_murcko(s) for s in chembl_active_smiles]
    chembl_scaffolds = [s for s in chembl_scaffolds if s]  # drop None/empty

    if not chembl_scaffolds:
        return ScaffoldResult(
            target_uniprot=target_uniprot, target_gene=target_gene,
            n_library=len(library_smiles), n_chembl_actives=0,
            top_chembl_scaffold_smiles="", top_chembl_scaffold_count=0,
            library_in_top_scaffold_pct=float("nan"),
            library_unique_scaffolds=0,
            decision="insufficient_chembl_actives",
        )

    counter = Counter(chembl_scaffolds)
    top_scaffold, top_count = counter.most_common(1)[0]
    n_lib_in_top = sum(1 for s in lib_scaffolds if s == top_scaffold and s)
    pct = 100.0 * n_lib_in_top / max(len(library_smiles), 1)
    n_unique_lib = len({s for s in lib_scaffolds if s})

    if pct > 60:
        decision = "rank_resolution_loss"        # Scenario 2 → LoRA worth
    elif pct < 25:
        decision = "manifold_mismatch"            # Scenario 1 → ensemble worth
    else:
        decision = "ambiguous_defer_to_diag_D"

    return ScaffoldResult(
        target_uniprot=target_uniprot, target_gene=target_gene,
        n_library=len(library_smiles), n_chembl_actives=len(chembl_scaffolds),
        top_chembl_scaffold_smiles=top_scaffold,
        top_chembl_scaffold_count=top_count,
        library_in_top_scaffold_pct=pct,
        library_unique_scaffolds=n_unique_lib,
        decision=decision,
    )
