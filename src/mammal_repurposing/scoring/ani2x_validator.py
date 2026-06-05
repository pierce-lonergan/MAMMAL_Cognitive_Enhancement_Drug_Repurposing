"""§8.9 — ANI-2x neural-potential validation for Boltz-2 poses.

Concept: after §7.5 + §7.17 give us pose centroids and (eventually) full
ligand+pocket coordinates from Boltz-2, validate each top-25 pose by computing
its ANI-2x energy and comparing to a baseline conformation. Large positive
ΔE (relative to a known good ligand pose) flags strain — the pose is high-
energy and probably wrong.

ANI-2x reference:
  Devereux et al. 2020 J Chem Theory Comput 16(7):4192 — extended ANI-1x to
  H/C/N/O/F/Cl/S. Trained on >25M DFT energies (ωB97X / 6-31G*).
  TorchANI: github.com/aiqm/torchani.

OPERATIONAL STATE (commit time):
  - TorchANI is NOT yet installed in mammal_env. Install: `pip install torchani`.
  - Pose mmCIF data not yet on the live grid (depends on §7.17 WSL2 re-run).
  - This module ships as a STUB with the correct API so when both
    dependencies arrive, it activates by setting `ani_available=True` and
    swapping `_compute_energy_stub` for a real `torchani.models.ANI2x()` call.

API:
    from mammal_repurposing.scoring.ani2x_validator import (
        validate_pose_set, summarise_validation,
    )
    results = validate_pose_set(pose_cif_paths, baseline_pose_path)
    summary = summarise_validation(results)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)


# Detect whether the heavy deps are installed; ship gracefully if not.
try:
    import torch                  # noqa: F401
    import torchani               # noqa: F401
    ANI_AVAILABLE = True
except ImportError:
    ANI_AVAILABLE = False


@dataclass
class PoseValidation:
    pose_path: str
    energy_kcal_mol: float | None = None
    energy_baseline: float | None = None
    delta_e_kcal_mol: float | None = None
    status: str = "STUB"          # STUB | OK | HIGH_STRAIN | ERROR
    n_atoms: int = 0
    note: str = ""


HIGH_STRAIN_THRESHOLD_KCAL = 50.0     # Δ above baseline → flag pose


def _parse_xyz_from_cif(cif_text: str) -> tuple[list[str], np.ndarray]:
    """Minimal mmCIF parser: pull (element_symbol_list, coord_array) for
    HETATM records (the docked ligand). Reuses the parser pattern from
    `pockets.pose_extract` but returns per-atom info instead of centroid.
    """
    from mammal_repurposing.pockets.pose_extract import _iter_atom_lines
    elements: list[str] = []
    coords: list[tuple[float, float, float]] = []
    for header, tokens in _iter_atom_lines(cif_text):
        try:
            cols = {c: i for i, c in enumerate(header)}
            group = tokens[cols["group_PDB"]]
            if group != "HETATM":
                continue
            elem = tokens[cols["type_symbol"]]
            if elem == "H":
                continue       # ANI-2x supports H but we strip for centroid consistency
            x = float(tokens[cols["Cartn_x"]])
            y = float(tokens[cols["Cartn_y"]])
            z = float(tokens[cols["Cartn_z"]])
            elements.append(elem)
            coords.append((x, y, z))
        except (KeyError, IndexError, ValueError):
            continue
    return elements, np.array(coords, dtype=float) if coords else np.zeros((0, 3))


def _compute_energy_stub(elements: list[str], coords: np.ndarray) -> float:
    """Stub used when TorchANI isn't installed. Returns a deterministic
    fake energy that scales with atom count so downstream code can still
    flow. Real implementation: torchani.models.ANI2x().forward((Z, xyz))."""
    if len(elements) == 0:
        return float("nan")
    # Synthetic: -50 kcal/mol per heavy atom (typical scale for small molecules)
    return -50.0 * len(elements)


def _compute_energy_real(elements: list[str], coords: np.ndarray) -> float:
    """Real ANI-2x computation. Only call when ANI_AVAILABLE."""
    import torch
    import torchani
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = torchani.models.ANI2x().to(device)
    # ANI species map: H=0, C=1, N=2, O=3, F=4, Cl=5, S=6
    species_map = {"H": 0, "C": 1, "N": 2, "O": 3, "F": 4, "Cl": 5, "S": 6}
    z = torch.tensor([[species_map.get(e, -1) for e in elements]],
                     dtype=torch.long, device=device)
    if (z < 0).any():
        raise ValueError(f"Element not supported by ANI-2x: {elements}")
    xyz = torch.tensor([coords], dtype=torch.float32, device=device)
    with torch.no_grad():
        _, e_hartree = model((z, xyz))
    return float(e_hartree.item() * 627.509)   # Hartree → kcal/mol


def compute_pose_energy(cif_path: Path | str) -> tuple[float, int]:
    """Read an mmCIF, extract ligand heavy atoms, return (energy_kcal, n_atoms)."""
    cif_text = Path(cif_path).read_text(encoding="utf-8")
    elements, coords = _parse_xyz_from_cif(cif_text)
    if len(coords) == 0:
        return float("nan"), 0
    if ANI_AVAILABLE:
        return _compute_energy_real(elements, coords), len(elements)
    return _compute_energy_stub(elements, coords), len(elements)


def validate_pose_set(
    pose_paths: list[Path | str],
    baseline_path: Path | str | None = None,
    high_strain_threshold: float = HIGH_STRAIN_THRESHOLD_KCAL,
) -> list[PoseValidation]:
    """Validate a list of poses; flag those with ΔE > threshold vs baseline."""
    baseline_e: float | None = None
    if baseline_path is not None:
        try:
            baseline_e, _ = compute_pose_energy(baseline_path)
        except Exception as e:
            logger.warning("Baseline energy failed: %s", e)
    out: list[PoseValidation] = []
    for p in pose_paths:
        try:
            e, n = compute_pose_energy(p)
            de = (e - baseline_e) if baseline_e is not None else None
            if not ANI_AVAILABLE:
                status = "STUB"
            elif de is not None and de > high_strain_threshold:
                status = "HIGH_STRAIN"
            else:
                status = "OK"
            out.append(PoseValidation(
                pose_path=str(p),
                energy_kcal_mol=e,
                energy_baseline=baseline_e,
                delta_e_kcal_mol=de,
                status=status,
                n_atoms=n,
                note=("ANI-2x not installed; stub energy used"
                      if not ANI_AVAILABLE else ""),
            ))
        except Exception as e:
            out.append(PoseValidation(
                pose_path=str(p),
                status="ERROR",
                note=str(e),
            ))
    return out


def summarise_validation(results: list[PoseValidation]) -> dict:
    """Roll-up stats for the report."""
    counts: dict[str, int] = {}
    energies: list[float] = []
    deltas: list[float] = []
    for r in results:
        counts[r.status] = counts.get(r.status, 0) + 1
        if r.energy_kcal_mol is not None and not np.isnan(r.energy_kcal_mol):
            energies.append(r.energy_kcal_mol)
        if r.delta_e_kcal_mol is not None and not np.isnan(r.delta_e_kcal_mol):
            deltas.append(r.delta_e_kcal_mol)
    return {
        "n_total": len(results),
        "status_counts": counts,
        "energy_mean": float(np.mean(energies)) if energies else None,
        "energy_std": float(np.std(energies)) if energies else None,
        "delta_e_mean": float(np.mean(deltas)) if deltas else None,
        "delta_e_max": float(np.max(deltas)) if deltas else None,
        "ani_available": ANI_AVAILABLE,
    }
