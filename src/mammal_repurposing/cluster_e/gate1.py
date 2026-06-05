"""V8 Gate 1: mechanism-class recovery from perturbational signatures.

The pre-registered V8 primary gate asks whether a compound's phenotypic
signature recovers its mechanism class. The class labels are pharmacology-
grounded (target / mechanism), NOT derived from the signatures themselves, so
the test is not circular. Clusters are scored against those labels with
Adjusted Mutual Information; the verdict bands follow the V8 OSF pre-registration
section 5.1: PASS at AMI >= 0.50 and ARI >= 0.40, DEGRADE in [0.30, 0.50),
otherwise FAIL (a publishable negative result).

This module holds the reusable, testable gate. The synthetic dry-run
(scripts/60) and the real-data runner (scripts/92) both compute the same
verdict so the pipeline is identical on synthetic and real inputs.
"""
from __future__ import annotations

import numpy as np

try:
    from sklearn.cluster import AgglomerativeClustering
    from sklearn.metrics import (
        adjusted_mutual_info_score, adjusted_rand_score,
        v_measure_score, fowlkes_mallows_score,
    )
    SKLEARN_AVAILABLE = True
except ImportError:                                            # pragma: no cover
    SKLEARN_AVAILABLE = False


def compound_consensus(X: np.ndarray, group_keys) -> tuple[np.ndarray, list]:
    """Collapse multiple signatures per compound to one consensus row (mean over
    that compound's signatures), so clustering is over compounds (the unit the
    mechanism-class label applies to) rather than over cell-line/dose replicates.

    Returns (X_consensus [n_compounds, n_features], ordered unique keys)."""
    keys = list(group_keys)
    uniq = sorted(set(keys))
    rows = []
    for k in uniq:
        mask = np.array([g == k for g in keys])
        rows.append(X[mask].mean(axis=0))
    return np.vstack(rows), uniq


def cluster_and_score(X: np.ndarray, labels, *, method: str = "agglomerative",
                      n_clusters: int | None = None,
                      hdbscan_min_size: int = 5,
                      leiden_gamma: float = 1.0, knn: int = 10) -> dict:
    """Cluster X and score the partition against integer/string `labels` with
    AMI / ARI / V-measure / Fowlkes-Mallows. `method` in
    {agglomerative, hdbscan, leiden}."""
    if not SKLEARN_AVAILABLE:                                  # pragma: no cover
        raise ImportError("scikit-learn required for clustering + metrics")
    truth = np.asarray(labels)
    n_truth = len(set(truth.tolist()))

    if method == "agglomerative":
        k = n_clusters or n_truth
        pred = AgglomerativeClustering(n_clusters=k).fit_predict(X)
    elif method == "hdbscan":
        import hdbscan
        pred = hdbscan.HDBSCAN(min_cluster_size=hdbscan_min_size).fit_predict(X)
    elif method == "leiden":
        import igraph as ig
        import leidenalg as la
        from sklearn.neighbors import NearestNeighbors
        nn = NearestNeighbors(n_neighbors=min(knn, len(X) - 1),
                              metric="cosine").fit(X)
        _, ind = nn.kneighbors(X)
        edges = [(int(i), int(j)) for i, nb in enumerate(ind) for j in nb[1:]]
        g = ig.Graph(edges=edges, directed=False)
        part = la.find_partition(g, la.RBConfigurationVertexPartition,
                                 resolution_parameter=leiden_gamma)
        pred = np.array(part.membership)
    else:
        raise ValueError(f"unknown method: {method}")

    return {
        "method": method,
        "n_clusters_predicted": int(len(set(pred.tolist()))),
        "n_clusters_truth": int(n_truth),
        "ami": float(adjusted_mutual_info_score(truth, pred)),
        "ari": float(adjusted_rand_score(truth, pred)),
        "v_measure": float(v_measure_score(truth, pred)),
        "fm": float(fowlkes_mallows_score(truth, pred)),
        "hdbscan_min_size": hdbscan_min_size if method == "hdbscan" else None,
        "leiden_gamma": leiden_gamma if method == "leiden" else None,
    }


def gate1_verdict(ami: float, ari: float) -> str:
    """V8.4 Gate 1 bands (V8 OSF pre-reg section 5.1)."""
    if ami >= 0.50 and ari >= 0.40:
        return "PASS"
    if ami >= 0.30 and ari >= 0.25:
        return "DEGRADE"
    return "FAIL"
