"""§7.17 — Heavy-atom-mean centroid extraction from Boltz mmCIF poses.

Stand-alone, dependency-light: parses the _atom_site loop directly from mmCIF
text without needing biopython. Suitable for inclusion in the WSL2 sweep
script, where adding biopython would balloon the venv.

The Boltz mmCIF output emits the ligand under HETATM with a non-standard
residue label_comp_id (e.g. 'LIG'). We collect every HETATM line, filter to
heavy atoms (type_symbol != 'H'), and take the per-axis mean.

This file is also exported from the package for use by downstream scripts
that already use biopython.
"""

from __future__ import annotations

import io
import logging
from pathlib import Path
from typing import Iterable

import numpy as np

logger = logging.getLogger(__name__)


def _split_loop_columns(loop_header: list[str]) -> dict[str, int]:
    """Map mmCIF column name → 0-based index in the data lines."""
    return {col: i for i, col in enumerate(loop_header)}


def _iter_atom_lines(cif_text: str) -> Iterable[list[str]]:
    """Yield split-token lists for every _atom_site data line in the mmCIF.

    Handles the standard loop_ block:
        loop_
        _atom_site.<col1>
        _atom_site.<col2>
        ...
        ATOM   ...
        HETATM ...

    Stops at the next # or loop_ block boundary.
    """
    lines = cif_text.splitlines()
    in_loop = False
    header: list[str] = []
    yielding_data = False

    for raw in lines:
        line = raw.strip()
        if not in_loop and line == "loop_":
            in_loop = True
            header = []
            yielding_data = False
            continue
        if in_loop and line.startswith("_atom_site."):
            header.append(line.split(".", 1)[1])
            continue
        if in_loop and header and (line.startswith("ATOM") or line.startswith("HETATM")):
            yielding_data = True
            yield header, line.split()
            continue
        if yielding_data and (not line or line == "#" or line.startswith("loop_")):
            # End of this atom_site block
            in_loop = False
            header = []
            yielding_data = False


def extract_ligand_centroid(
    cif_text: str,
    ligand_label_comp_id: str | None = None,
) -> np.ndarray | None:
    """Return the heavy-atom-mean (x, y, z) for the ligand record set, or None.

    Args:
        cif_text: full mmCIF file contents as a string.
        ligand_label_comp_id: optional filter on the residue label (e.g. 'LIG');
            when omitted, ALL HETATM records are pooled (typical for Boltz
            single-ligand poses where there's exactly one HETATM block).

    Returns:
        np.ndarray shape (3,) [x, y, z] in Å, or None if no heavy HETATM found.
    """
    coords: list[tuple[float, float, float]] = []
    for header, tokens in _iter_atom_lines(cif_text):
        cols = _split_loop_columns(header)
        try:
            group = tokens[cols["group_PDB"]]
        except (KeyError, IndexError):
            continue
        if group != "HETATM":
            continue
        if ligand_label_comp_id is not None:
            try:
                if tokens[cols["label_comp_id"]] != ligand_label_comp_id:
                    continue
            except (KeyError, IndexError):
                continue
        try:
            elem = tokens[cols["type_symbol"]]
        except (KeyError, IndexError):
            elem = ""
        if elem == "H":
            continue
        try:
            x = float(tokens[cols["Cartn_x"]])
            y = float(tokens[cols["Cartn_y"]])
            z = float(tokens[cols["Cartn_z"]])
        except (KeyError, ValueError, IndexError):
            continue
        coords.append((x, y, z))
    if not coords:
        return None
    arr = np.asarray(coords, dtype=float)
    return arr.mean(axis=0)


def extract_centroid_from_file(
    cif_path: Path | str,
    ligand_label_comp_id: str | None = None,
) -> np.ndarray | None:
    """Convenience wrapper around extract_ligand_centroid that reads the file."""
    p = Path(cif_path)
    if not p.exists():
        return None
    return extract_ligand_centroid(p.read_text(encoding="utf-8"),
                                   ligand_label_comp_id=ligand_label_comp_id)
