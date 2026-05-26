"""Geometric pocket classifier — Boltz-2 pose centroid → pocket_class.

Per research/4-tier/Pocket-Conditioned-Boltz2.md §1 decision tree:

  d_known     = min distance to any orthosteric / allosteric_known centroid
  d_buried    = secondary threshold for no_pocket_match vs surface_artifact

  IF d_ortho ≤ 8 Å                  → orthosteric
  ELIF d_allo  ≤ 8 Å                → allosteric_known (or allosteric_putative)
  ELIF d_any ≤ 15 Å (no known match) → no_pocket_match
  ELSE                              → surface_artifact

Plus a 'dual_site' detection: if pose contacts BOTH centroids of a
multi-pocket target (e.g. ACHE CAS + PAS), the class is the special
'orthosteric' variant for donepezil-mode.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np

from .pocket_database import PocketDatabase, PocketSpec

logger = logging.getLogger(__name__)

# Rank multipliers per research doc §3.2 (applied after RRF, before per-target ρ
# calibration).
RANK_MULTIPLIER = {
    "orthosteric":         1.00,
    "allosteric_known":    1.00,
    "allosteric_putative": 1.00,    # tagged manual_review=True
    "cryptic_predicted":   1.00,    # tagged manual_review=True
    "no_pocket_match":     0.50,
    "surface_artifact":    0.50,
    "NA_no_pose":          0.30,
}


@dataclass
class PocketClassification:
    pocket_class: str
    matched_pocket_tag: str        # which pocket spec matched (or empty for no_pocket_match)
    matched_pdb: str
    distance_to_match: float       # Å
    detector_votes: dict[str, float]
    manual_review: bool
    pose_centroid: tuple[float, float, float]
    rank_multiplier: float
    is_dual_site: bool = False     # both CAS + PAS or similar dual-pocket contact


def classify_pose(
    pose_centroid: np.ndarray,
    target_gene: str,
    db: PocketDatabase,
) -> PocketClassification:
    """Classify a single Boltz-2 pose centroid against the per-target pocket DB.

    Args:
        pose_centroid: (x, y, z) heavy-atom mean of the docked ligand
        target_gene: gene symbol (e.g. 'CHRNA7')
        db: loaded PocketDatabase

    Returns:
        PocketClassification with class label, matched pocket, distance, etc.
    """
    pose_centroid = np.asarray(pose_centroid, dtype=float).flatten()
    if pose_centroid.shape != (3,):
        return PocketClassification(
            pocket_class="NA_no_pose", matched_pocket_tag="",
            matched_pdb="", distance_to_match=float("nan"),
            detector_votes={}, manual_review=False,
            pose_centroid=(float("nan"),) * 3,
            rank_multiplier=RANK_MULTIPLIER["NA_no_pose"],
        )

    specs = db.pockets_by_target.get(target_gene, [])
    specs_with_centroid = [s for s in specs if s.centroid is not None]

    if not specs_with_centroid:
        return PocketClassification(
            pocket_class="NA_no_pose",
            matched_pocket_tag="",
            matched_pdb="",
            distance_to_match=float("nan"),
            detector_votes={},
            manual_review=False,
            pose_centroid=tuple(pose_centroid.tolist()),
            rank_multiplier=RANK_MULTIPLIER["NA_no_pose"],
        )

    # Compute distance to every pocket centroid in this target
    distances = []
    for spec in specs_with_centroid:
        d = float(np.linalg.norm(pose_centroid - spec.centroid))
        distances.append((spec, d))
    distances.sort(key=lambda t: t[1])
    closest_spec, closest_d = distances[0]

    # Check for dual-site (e.g. donepezil CAS+PAS): if TWO centroids both ≤
    # known threshold, label as orthosteric-dual. Per research doc §2G:
    # "a 'dual' pose is one whose ligand heavy-atom span makes contact within
    # 5 Å of BOTH centroids" — donepezil-mode CAS+PAS span.
    within_known = [(s, d) for s, d in distances if d <= db.distance_threshold_known]
    is_dual = len(within_known) >= 2 and any(
        s.pocket_class == "orthosteric" for s, _ in within_known
    ) and any(
        s.pocket_class == "allosteric_known" for s, _ in within_known
    )

    # Decision tree
    if closest_d <= db.distance_threshold_known:
        if is_dual:
            # Donepezil-mode: pose spans both subsites — call it `orthosteric`
            # per research doc §2G; the orthosteric subsite (CAS for ACHE,
            # OBP for DRD1) is the dominant pharmacology.
            pclass = "orthosteric"
            manual = False
            # Re-select the matched_pocket_tag to the orthosteric one
            ortho_specs = [(s, d) for s, d in within_known if s.pocket_class == "orthosteric"]
            if ortho_specs:
                closest_spec, closest_d = ortho_specs[0]
        else:
            pclass = closest_spec.pocket_class
            manual = pclass == "allosteric_putative" or pclass == "cryptic_predicted"
    elif closest_d <= db.distance_threshold_buried:
        pclass = "no_pocket_match"
        manual = False
    else:
        pclass = "surface_artifact"
        manual = False

    return PocketClassification(
        pocket_class=pclass,
        matched_pocket_tag=closest_spec.tag if closest_d <= db.distance_threshold_buried else "",
        matched_pdb=closest_spec.pdb if closest_d <= db.distance_threshold_buried else "",
        distance_to_match=closest_d,
        detector_votes={"geometric": 1.0},
        manual_review=manual,
        pose_centroid=tuple(pose_centroid.tolist()),
        rank_multiplier=RANK_MULTIPLIER.get(pclass, 1.0),
        is_dual_site=is_dual,
    )
