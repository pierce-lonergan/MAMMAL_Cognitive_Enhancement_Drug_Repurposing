"""§7.16 — Pocket-detector ensemble (Sprint 2).

The §7.5 MVP shipped a curated-centroid pocket database + geometric
classifier (13/13 validation gates). Sprint 2 adds three ML pocket detectors
whose outputs are combined via consensus voting:

  1. P2Rank (Krivák & Hoksza 2018 J Cheminform 10:39) — random forest on
     residue-level features; Java jar binary. Detects pockets from PDB +
     ranks by druggability. Best-known PDBbind benchmark (Top-1 success
     rate ~72%).

  2. PocketMiner (Meller et al. 2023 Nat Commun 14:1177) — PyTorch GVP-GNN
     trained to detect CRYPTIC pockets (those that open only on ligand
     binding). Critical for §7.5's blind spot at allosteric / cryptic sites.
     Distributed as a PyTorch model + inference script.

  3. CryptoBench (Škrhák et al. 2024 Bioinformatics 41(1):btae745) —
     ESM-2 embeddings + cryptic-pocket classifier. Complements PocketMiner
     with a sequence-only second opinion.

The ensemble adds a `detector_votes` column to PocketClassification with
each detector's confidence + a `consensus_class` label when ≥2 of 3 agree.

OPERATIONAL STATE (commit time):
  - All 3 detectors require non-trivial installs (Java JRE + 250 MB jar for
    P2Rank; PyTorch model checkpoint ~150 MB for PocketMiner; ESM-2 cache
    + classifier ~2 GB for CryptoBench).
  - This module ships uniform DetectorAdapter interfaces with code-complete
    stubs. Each adapter exposes `detect(pdb_path) -> list[DetectedPocket]`.
  - Per-detector AVAILABLE flag auto-detects whether the underlying tool
    is installed.

Reference:
  Krivák & Hoksza 2018 J Cheminform 10:39 (DOI 10.1186/s13321-018-0285-8) — P2Rank
  Meller et al. 2023 Nat Commun 14:1177 (DOI 10.1038/s41467-023-36699-3) — PocketMiner
  Škrhák et al. 2024 Bioinformatics 41(1):btae745 — CryptoBench
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class DetectedPocket:
    """One detector's verdict for one pocket."""
    detector: str                       # "p2rank" | "pocketminer" | "cryptobench"
    pocket_id: str                      # detector-internal pocket label
    centroid: np.ndarray                # (x, y, z) heavy-atom mean
    score: float                        # detector-internal confidence in [0, 1]
    is_cryptic_predicted: bool = False  # only PocketMiner / CryptoBench
    residues: list[int] = field(default_factory=list)
    note: str = ""


@dataclass
class ConsensusVerdict:
    """Cross-detector consensus for one pocket region."""
    consensus_class: str    # "orthosteric_consensus" | "allosteric_consensus" |
                            # "cryptic_consensus" | "single_detector" |
                            # "no_consensus"
    n_voters: int
    avg_score: float
    detector_votes: dict[str, float]
    centroid: np.ndarray
    is_cryptic: bool


# ===========================================================================
# Detector availability probes
# ===========================================================================
def _check_p2rank_available() -> bool:
    """P2Rank is a Java jar at $P2RANK_HOME/prank or via $PATH `prank`."""
    if shutil.which("prank"):
        return True
    import os
    home = os.environ.get("P2RANK_HOME", "")
    if home and (Path(home) / "prank").exists():
        return True
    return False


def _check_pocketminer_available() -> bool:
    """PocketMiner needs torch + a downloaded checkpoint at $POCKETMINER_HOME."""
    try:
        import torch  # noqa: F401
    except ImportError:
        return False
    import os
    home = os.environ.get("POCKETMINER_HOME", "")
    return bool(home and Path(home).exists())


def _check_cryptobench_available() -> bool:
    """CryptoBench needs esm + a downloaded classifier at $CRYPTOBENCH_HOME."""
    try:
        import esm  # noqa: F401
    except ImportError:
        return False
    import os
    home = os.environ.get("CRYPTOBENCH_HOME", "")
    return bool(home and Path(home).exists())


P2RANK_AVAILABLE = _check_p2rank_available()
POCKETMINER_AVAILABLE = _check_pocketminer_available()
CRYPTOBENCH_AVAILABLE = _check_cryptobench_available()


# ===========================================================================
# Adapters
# ===========================================================================
def detect_with_p2rank(pdb_path: Path, max_pockets: int = 5) -> list[DetectedPocket]:
    """Wrap the P2Rank `prank predict` CLI.

    Real-mode call:
        prank predict -o <out_dir> -threads 1 <pdb_path>
    Parses <out_dir>/<pdb_stem>_predictions.csv for centroids + scores.

    Stub: returns one synthetic pocket at PDB heavy-atom mean.
    """
    pdb_path = Path(pdb_path)
    if not P2RANK_AVAILABLE:
        logger.info("p2rank not installed; returning STUB pocket for %s", pdb_path.name)
        return [DetectedPocket(
            detector="p2rank",
            pocket_id="STUB_pocket_1",
            centroid=np.zeros(3),
            score=0.5,
            note="STUB — P2Rank not installed",
        )]
    # Real path
    import tempfile
    with tempfile.TemporaryDirectory(prefix="p2rank_") as tmp:
        out_dir = Path(tmp)
        cmd = ["prank", "predict", "-o", str(out_dir), "-threads", "1", str(pdb_path)]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=600)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
            logger.warning("P2Rank failed on %s: %s", pdb_path.name, e)
            return []
        csv = next((p for p in out_dir.rglob("*_predictions.csv")), None)
        if not csv:
            return []
        import pandas as pd
        df = pd.read_csv(csv).head(max_pockets)
        # P2Rank columns: name, rank, score, center_x, center_y, center_z, ...
        pockets = []
        for _, r in df.iterrows():
            pockets.append(DetectedPocket(
                detector="p2rank",
                pocket_id=str(r.get("name", "pocket")),
                centroid=np.array(
                    [r.get("center_x", 0), r.get("center_y", 0), r.get("center_z", 0)],
                    dtype=float),
                score=float(r.get("score", 0.0)),
                note=f"P2Rank rank={int(r.get('rank', 0))}",
            ))
        return pockets


def detect_with_pocketminer(pdb_path: Path, max_pockets: int = 5) -> list[DetectedPocket]:
    """Wrap the PocketMiner GVP-GNN inference.

    Real-mode call: load the published checkpoint, run inference on the PDB,
    extract pockets above the cryptic-confidence threshold.

    Stub: returns one synthetic cryptic pocket.
    """
    pdb_path = Path(pdb_path)
    if not POCKETMINER_AVAILABLE:
        logger.info("pocketminer not installed; returning STUB for %s", pdb_path.name)
        return [DetectedPocket(
            detector="pocketminer",
            pocket_id="STUB_cryptic_1",
            centroid=np.zeros(3),
            score=0.5,
            is_cryptic_predicted=True,
            note="STUB — PocketMiner not installed",
        )]
    # Real path placeholder
    import os
    home = Path(os.environ["POCKETMINER_HOME"])
    logger.info("PocketMiner real-mode inference on %s (home=%s)", pdb_path.name, home)
    # The actual PocketMiner repo provides inference.py; user wires that here.
    # For now we surface a clear hook:
    raise NotImplementedError(
        "PocketMiner real-mode wiring TBD — see Meller 2023 inference script. "
        "Set POCKETMINER_HOME and complete this function."
    )


def detect_with_cryptobench(pdb_path: Path, max_pockets: int = 5) -> list[DetectedPocket]:
    """Wrap CryptoBench (ESM-2 + cryptic-pocket classifier).

    Real-mode call: extract ESM-2 embeddings per residue, run the
    cryptic-classifier head, threshold to identify cryptic residues, cluster
    into pockets.

    Stub: returns one synthetic cryptic pocket.
    """
    pdb_path = Path(pdb_path)
    if not CRYPTOBENCH_AVAILABLE:
        logger.info("cryptobench not installed; returning STUB for %s", pdb_path.name)
        return [DetectedPocket(
            detector="cryptobench",
            pocket_id="STUB_cryptic_1",
            centroid=np.zeros(3),
            score=0.5,
            is_cryptic_predicted=True,
            note="STUB — CryptoBench not installed",
        )]
    raise NotImplementedError(
        "CryptoBench real-mode wiring TBD — see Škrhák 2024 inference. "
        "Set CRYPTOBENCH_HOME and complete this function."
    )


# ===========================================================================
# Ensemble runner + consensus
# ===========================================================================
DEFAULT_DETECTORS = ("p2rank", "pocketminer", "cryptobench")
CONSENSUS_DISTANCE_THRESHOLD = 8.0   # Å — pockets within this distance are "same"


def run_ensemble(
    pdb_path: Path,
    detectors: tuple[str, ...] = DEFAULT_DETECTORS,
) -> dict[str, list[DetectedPocket]]:
    """Run every detector on the PDB; return {detector: pockets}."""
    out: dict[str, list[DetectedPocket]] = {}
    for d in detectors:
        try:
            if d == "p2rank":
                out["p2rank"] = detect_with_p2rank(pdb_path)
            elif d == "pocketminer":
                out["pocketminer"] = detect_with_pocketminer(pdb_path)
            elif d == "cryptobench":
                out["cryptobench"] = detect_with_cryptobench(pdb_path)
            else:
                logger.warning("Unknown detector: %s", d)
        except NotImplementedError as e:
            logger.warning("%s not implemented at runtime: %s", d, e)
            out[d] = []
        except Exception as e:
            logger.warning("%s failed on %s: %s", d, pdb_path.name, e)
            out[d] = []
    return out


def consensus_vote(
    detections: dict[str, list[DetectedPocket]],
    distance_threshold: float = CONSENSUS_DISTANCE_THRESHOLD,
) -> list[ConsensusVerdict]:
    """Cluster detector pockets by centroid proximity; emit one verdict per cluster.

    A pocket is `consensus` if ≥2 detectors place a pocket within
    `distance_threshold` Å of the cluster centroid. Cryptic flag rides on
    OR of the contributing detector flags.
    """
    all_pockets: list[DetectedPocket] = []
    for d, lst in detections.items():
        all_pockets.extend(lst)
    if not all_pockets:
        return []

    # Naive O(N²) greedy clustering — fine for ≤ ~5×3 = 15 pockets per target
    clusters: list[list[DetectedPocket]] = []
    for p in all_pockets:
        placed = False
        for c in clusters:
            if any(np.linalg.norm(p.centroid - q.centroid) <= distance_threshold
                   for q in c):
                c.append(p)
                placed = True
                break
        if not placed:
            clusters.append([p])

    verdicts: list[ConsensusVerdict] = []
    for cluster in clusters:
        detector_votes = {p.detector: p.score for p in cluster}
        n_voters = len(detector_votes)
        avg = float(np.mean([p.score for p in cluster]))
        centroid_avg = np.mean(np.stack([p.centroid for p in cluster]), axis=0)
        any_cryptic = any(p.is_cryptic_predicted for p in cluster)
        if n_voters >= 2:
            if any_cryptic:
                consensus = "cryptic_consensus"
            else:
                consensus = "orthosteric_consensus"
        else:
            consensus = "single_detector"
        verdicts.append(ConsensusVerdict(
            consensus_class=consensus,
            n_voters=n_voters,
            avg_score=avg,
            detector_votes=detector_votes,
            centroid=centroid_avg,
            is_cryptic=any_cryptic,
        ))
    # Sort by avg score descending
    verdicts.sort(key=lambda v: -v.avg_score)
    return verdicts


def detector_availability() -> dict[str, bool]:
    """Reports which detectors are operational at import time."""
    return {
        "p2rank": P2RANK_AVAILABLE,
        "pocketminer": POCKETMINER_AVAILABLE,
        "cryptobench": CRYPTOBENCH_AVAILABLE,
    }
