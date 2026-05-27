"""§8.0a — Multi-axis Pareto NSGA-III non-dominated sort.

Restructures the wet-lab shortlist as a Pareto frontier across 5 axes:

  axis 1 — efficacy:    RRF score (HIGHER is better)
  axis 2 — safety:      −(n_tier_1 + 0.5 × n_tier_2) liability hits      (HIGHER is better)
  axis 3 — selectivity: top-target Partition Index                       (HIGHER is better)
  axis 4 — IP freedom:  CTgov ip_status mapped to numeric novelty score  (HIGHER is better)
                            approved        → 0.0
                            investigational → 0.4
                            early           → 0.7
                            none (novel)    → 1.0
  axis 5 — novelty:     nootropic-similarity tag                         (HIGHER is better)
                            analog          → 0.0
                            intermediate    → 0.5
                            novel_scaffold  → 1.0

A compound `a` dominates `b` if `a` is at least as good as `b` on every axis
AND strictly better on at least one. The Pareto front (frontier-rank 0) is the
set of non-dominated compounds. Subsequent frontiers (rank 1, 2, ...) are
peeled off in turn.

Hypervolume + spread metrics are computed on the front for the audit.

Reference: Deb & Jain 2014 IEEE Trans Evol Comput 18(4):577 — NSGA-III.
Hand-rolled (no external dep) for the 5-axis × ≤300-compound problem size.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Iterable

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


IP_NOVELTY_SCORE: dict[str, float] = {
    "approved": 0.0,
    "investigational": 0.4,
    "early": 0.7,
    "none": 1.0,
    "": 0.5,                # missing → neutral
}

NOVELTY_TAG_SCORE: dict[str, float] = {
    "analog": 0.0,
    "intermediate": 0.5,
    "novel_scaffold": 1.0,
    "unknown": 0.3,
    "": 0.3,
}


@dataclass
class ParetoConfig:
    rrf_col: str = "rrf_score"
    n_tier_1_col: str = "n_tier_1"
    n_tier_2_col: str = "n_tier_2"
    pi_col: str = "partition_index_top"
    ip_status_col: str = "ip_status"
    novelty_tag_col: str = "nootropic_novelty_tag"


def build_axis_matrix(
    df: pd.DataFrame,
    config: ParetoConfig | None = None,
) -> tuple[np.ndarray, list[str]]:
    """Return (N × 5) axis matrix (every axis: higher=better) + axis names."""
    cfg = config or ParetoConfig()
    n = len(df)
    A = np.zeros((n, 5), dtype=float)

    A[:, 0] = df[cfg.rrf_col].to_numpy(dtype=float)              # efficacy
    # safety = -(n_tier_1 + 0.5*n_tier_2). Default 0 if missing.
    n_t1 = df.get(cfg.n_tier_1_col, pd.Series(np.zeros(n)))
    n_t2 = df.get(cfg.n_tier_2_col, pd.Series(np.zeros(n)))
    A[:, 1] = -(n_t1.to_numpy(dtype=float) + 0.5 * n_t2.to_numpy(dtype=float))
    # selectivity = top-target Partition Index (higher = more selective)
    pi = df.get(cfg.pi_col, pd.Series(np.zeros(n)))
    A[:, 2] = pi.to_numpy(dtype=float)
    # IP novelty score
    ip = df.get(cfg.ip_status_col, pd.Series([""] * n)).astype(str)
    A[:, 3] = ip.map(IP_NOVELTY_SCORE).fillna(0.5).to_numpy(dtype=float)
    # Nootropic novelty score
    nv = df.get(cfg.novelty_tag_col, pd.Series([""] * n)).astype(str)
    A[:, 4] = nv.map(NOVELTY_TAG_SCORE).fillna(0.3).to_numpy(dtype=float)

    return A, ["efficacy_rrf", "safety_neg_liability",
               "selectivity_pi", "ip_novelty", "scaffold_novelty"]


def _dominates(a: np.ndarray, b: np.ndarray, tol: float = 1e-9) -> bool:
    """a dominates b iff a[i] ≥ b[i] for all i AND a[i] > b[i] for some i."""
    return np.all(a >= b - tol) and np.any(a > b + tol)


def non_dominated_sort(A: np.ndarray) -> list[int]:
    """Return a list of frontier-rank assignments (0 = Pareto front, 1 = next, ...).

    Naive O(N²×D) but fine for N ≤ 300, D = 5.
    """
    n = A.shape[0]
    ranks = np.full(n, -1, dtype=int)
    domination_counts = np.zeros(n, dtype=int)
    dominated_by = [[] for _ in range(n)]

    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if _dominates(A[i], A[j]):
                dominated_by[i].append(j)
            elif _dominates(A[j], A[i]):
                domination_counts[i] += 1

    current_rank = 0
    current_front = [i for i in range(n) if domination_counts[i] == 0]
    while current_front:
        for i in current_front:
            ranks[i] = current_rank
        next_front: list[int] = []
        for i in current_front:
            for j in dominated_by[i]:
                domination_counts[j] -= 1
                if domination_counts[j] == 0:
                    next_front.append(j)
        current_rank += 1
        current_front = next_front

    return ranks.tolist()


def hypervolume_2d(front: np.ndarray, reference_point: np.ndarray) -> float:
    """Hypervolume of a 2D Pareto front w.r.t. a reference point (lower-left).

    For >2D the exact computation is exponential; we use a Monte-Carlo
    estimator instead (see `hypervolume_mc`).
    """
    if front.shape[1] != 2 or len(front) == 0:
        return 0.0
    # Sort by first axis descending
    f = front[np.argsort(-front[:, 0])]
    hv = 0.0
    prev_y = reference_point[1]
    for x, y in f:
        if y > prev_y:
            hv += (x - reference_point[0]) * (y - prev_y)
            prev_y = y
    return float(hv)


def hypervolume_mc(
    front: np.ndarray,
    reference_point: np.ndarray,
    n_samples: int = 100_000,
    seed: int = 42,
) -> float:
    """Monte-Carlo hypervolume estimator for arbitrary dimension."""
    if len(front) == 0:
        return 0.0
    d = front.shape[1]
    # Bounding box: [ref, ideal]
    ideal = front.max(axis=0)
    rng = np.random.default_rng(seed)
    samples = rng.uniform(reference_point, ideal, size=(n_samples, d))
    # A sample is dominated by the front if ANY front point dominates it
    # (component-wise ≥).
    dominated = np.zeros(n_samples, dtype=bool)
    for p in front:
        dominated |= np.all(samples <= p, axis=1)
    volume_box = float(np.prod(ideal - reference_point))
    return float(volume_box * dominated.mean())


def crowding_distance(front_points: np.ndarray) -> np.ndarray:
    """Per-point crowding distance — diversity within a single Pareto front.

    Used as the secondary sort key for picking representative compounds when
    the front contains too many candidates to ship.
    """
    n, d = front_points.shape
    cd = np.zeros(n, dtype=float)
    if n <= 2:
        cd[:] = np.inf
        return cd
    for axis in range(d):
        order = np.argsort(front_points[:, axis])
        cd[order[0]] = np.inf
        cd[order[-1]] = np.inf
        vmin, vmax = front_points[order[0], axis], front_points[order[-1], axis]
        span = vmax - vmin
        if span == 0:
            continue
        for k in range(1, n - 1):
            cd[order[k]] += (
                front_points[order[k + 1], axis] - front_points[order[k - 1], axis]
            ) / span
    return cd


def rank_pareto(
    df: pd.DataFrame,
    config: ParetoConfig | None = None,
) -> pd.DataFrame:
    """Augment df with `pareto_rank` (frontier number) + `crowding_distance`
    + per-axis scaled value columns. Returns a copy."""
    cfg = config or ParetoConfig()
    A, axis_names = build_axis_matrix(df, cfg)
    ranks = non_dominated_sort(A)

    out = df.copy()
    out["pareto_rank"] = ranks
    for i, name in enumerate(axis_names):
        out[f"_axis_{name}"] = A[:, i]

    # Crowding distance computed per-front
    out["crowding_distance"] = 0.0
    for r in sorted(set(ranks)):
        idx = np.array([i for i, x in enumerate(ranks) if x == r])
        cd = crowding_distance(A[idx])
        for k, i in enumerate(idx):
            out.iat[i, out.columns.get_loc("crowding_distance")] = float(cd[k])

    return out
