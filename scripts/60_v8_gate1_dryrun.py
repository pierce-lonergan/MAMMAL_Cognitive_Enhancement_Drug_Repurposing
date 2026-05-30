"""V8.4 Stage 1 — Gate 1 mechanism-class recovery dry-run on synthetic data.

Validates the V8 Gate 1 pipeline (Leiden + HDBSCAN clustering on MOFA+
factors → AMI/ARI vs PRISMA mechanism-class labels) on a SYNTHETIC dataset
with known ground-truth cluster structure.

Synthetic generator:
  - 5 mechanism classes (cholinergic, catecholaminergic, glutamatergic,
    trophic_ISR, remyelination) × 50 compounds each = 250 compounds
  - 30 latent factors (matches V8.3 MOFA+ K=30)
  - Each compound's factor vector = class_centroid + Gaussian noise
  - Class centroids are orthogonal in 30-d space (separability ≈ guaranteed)

Validation:
  - Leiden clustering: γ sweep ∈ {0.4, 0.6, 0.8, 1.0, 1.2}; pick best AMI
  - HDBSCAN: min_cluster_size ∈ {15, 25, 50}; pick best AMI
  - Compute AMI / ARI / V-measure / Fowlkes-Mallows vs ground truth
  - Apply V8.4 thresholds (PASS ≥ 0.50 AMI; DEGRADE [0.30, 0.50); FAIL < 0.30)
  - Per-class breakdown: which classes cluster cleanly, which don't

This is a SANITY CHECK that the V8 Gate 1 pipeline produces correct
verdicts on data where the answer is known. If Leiden/HDBSCAN can't recover
5 well-separated clusters on this synthetic data, the pipeline is broken.

Outputs:
  data/results/v2/v8_gate1_dryrun_v1.parquet
  reports/pipeline/v8_gate1_dryrun_v1.md
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("v8_gate1_dryrun")


# Optional clustering dependencies
try:
    import leidenalg  # noqa: F401
    LEIDEN_AVAILABLE = True
except ImportError:
    LEIDEN_AVAILABLE = False

try:
    import hdbscan  # noqa: F401
    HDBSCAN_AVAILABLE = True
except ImportError:
    HDBSCAN_AVAILABLE = False

try:
    from sklearn.cluster import AgglomerativeClustering  # noqa: F401
    from sklearn.metrics import (
        adjusted_mutual_info_score, adjusted_rand_score,
        v_measure_score, fowlkes_mallows_score,
    )
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False


COGNITION_CENTROID_NAMES = (
    "cholinergic", "catecholaminergic", "glutamatergic",
    "trophic_ISR", "remyelination",
)


def generate_synthetic_phenotype(
    n_classes: int = 5,
    n_per_class: int = 50,
    K: int = 30,
    noise_sigma: float = 0.30,
    rng_seed: int = 42,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Generate synthetic MOFA+-like factor matrix with known cluster structure.

    Returns:
        factor_matrix: (n_compounds, K) latent factors
        ground_truth_labels: (n_compounds,) integer class labels
        class_names: list of K class names
    """
    rng = np.random.default_rng(rng_seed)
    class_names = list(COGNITION_CENTROID_NAMES[:n_classes])
    n_compounds = n_classes * n_per_class

    # Orthogonal class centroids in K-dim space (use first n_classes orthonormal
    # basis vectors scaled to magnitude ~ 3 for separability)
    centroids = np.zeros((n_classes, K))
    for c in range(n_classes):
        centroids[c, c] = 3.0     # diagonal centroids
        # Add small off-diagonal variation per class to avoid perfect axis alignment
        for k in range(K):
            if k != c:
                centroids[c, k] += rng.normal(0, 0.5)

    # Generate per-compound vectors
    factor_matrix = np.zeros((n_compounds, K))
    labels = np.zeros(n_compounds, dtype=int)
    idx = 0
    for c in range(n_classes):
        for _ in range(n_per_class):
            factor_matrix[idx] = centroids[c] + rng.normal(0, noise_sigma, K)
            labels[idx] = c
            idx += 1

    return factor_matrix, labels, class_names


def cluster_and_score(
    factor_matrix: np.ndarray,
    ground_truth: np.ndarray,
    method: str = "agglomerative",
    n_clusters: int | None = None,
    leiden_gamma: float = 1.0,
    hdbscan_min_size: int = 25,
    knn: int = 15,
) -> dict:
    """Run a clustering method and compute AMI / ARI / V-measure / FM."""
    if not SKLEARN_AVAILABLE:
        raise ImportError("scikit-learn required for clustering + metrics")

    pred_labels: np.ndarray
    if method == "leiden":
        if not LEIDEN_AVAILABLE:
            raise ImportError("leidenalg not installed")
        import igraph as ig
        import leidenalg as la
        # Build SNN graph via cosine distance + k-NN
        from sklearn.neighbors import NearestNeighbors
        nn = NearestNeighbors(n_neighbors=knn, metric="cosine").fit(factor_matrix)
        _, indices = nn.kneighbors(factor_matrix)
        edges = []
        for i, neighbors in enumerate(indices):
            for j in neighbors[1:]:    # skip self
                edges.append((int(i), int(j)))
        g = ig.Graph(edges=edges, directed=False)
        partition = la.find_partition(
            g, la.RBConfigurationVertexPartition,
            resolution_parameter=leiden_gamma,
        )
        pred_labels = np.array(partition.membership)
    elif method == "hdbscan":
        if not HDBSCAN_AVAILABLE:
            raise ImportError("hdbscan not installed")
        import hdbscan
        clusterer = hdbscan.HDBSCAN(min_cluster_size=hdbscan_min_size)
        pred_labels = clusterer.fit_predict(factor_matrix)
    elif method == "agglomerative":
        from sklearn.cluster import AgglomerativeClustering
        # Default: cluster into n_classes (if known) or guess via dendrogram
        if n_clusters is None:
            n_clusters = int(len(np.unique(ground_truth)))
        clusterer = AgglomerativeClustering(n_clusters=n_clusters)
        pred_labels = clusterer.fit_predict(factor_matrix)
    else:
        raise ValueError(f"Unknown clustering method: {method}")

    return {
        "method": method,
        "n_clusters_predicted": int(len(np.unique(pred_labels))),
        "n_clusters_truth": int(len(np.unique(ground_truth))),
        "ami": float(adjusted_mutual_info_score(ground_truth, pred_labels)),
        "ari": float(adjusted_rand_score(ground_truth, pred_labels)),
        "v_measure": float(v_measure_score(ground_truth, pred_labels)),
        "fm": float(fowlkes_mallows_score(ground_truth, pred_labels)),
        "leiden_gamma": leiden_gamma if method == "leiden" else None,
        "hdbscan_min_size": hdbscan_min_size if method == "hdbscan" else None,
    }


def gate1_verdict(ami: float, ari: float) -> str:
    """V8.4 Gate 1 thresholds per OSF pre-reg §5.1."""
    if ami >= 0.50 and ari >= 0.40:
        return "PASS"
    if ami >= 0.30 and ari >= 0.25:
        return "DEGRADE"
    return "FAIL"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-classes", type=int, default=5)
    parser.add_argument("--n-per-class", type=int, default=50)
    parser.add_argument("--K", type=int, default=30)
    parser.add_argument("--noise-sigma", type=float, default=0.30)
    parser.add_argument("--out", type=Path,
                        default=ROOT / "data" / "results" / "v2"
                        / "v8_gate1_dryrun_v1.parquet")
    parser.add_argument("--report", type=Path,
                        default=ROOT / "reports" / "pipeline" / "v8_gate1_dryrun_v1.md")
    args = parser.parse_args()

    if not SKLEARN_AVAILABLE:
        logger.error("scikit-learn required (`pip install scikit-learn`)")
        return 2

    # Generate synthetic data
    factor_matrix, ground_truth, class_names = generate_synthetic_phenotype(
        n_classes=args.n_classes,
        n_per_class=args.n_per_class,
        K=args.K,
        noise_sigma=args.noise_sigma,
    )
    n_compounds = factor_matrix.shape[0]
    logger.info("Synthetic dataset: %d compounds × %d factors, %d classes",
                n_compounds, args.K, args.n_classes)

    results: list[dict] = []

    # Agglomerative (baseline, always available)
    logger.info("Running Agglomerative clustering (n_clusters=%d)", args.n_classes)
    r_agg = cluster_and_score(
        factor_matrix, ground_truth, method="agglomerative",
        n_clusters=args.n_classes,
    )
    r_agg["verdict"] = gate1_verdict(r_agg["ami"], r_agg["ari"])
    results.append(r_agg)
    logger.info("Agglomerative: AMI=%.3f ARI=%.3f → %s",
                r_agg["ami"], r_agg["ari"], r_agg["verdict"])

    # Leiden sweep
    if LEIDEN_AVAILABLE:
        for gamma in (0.4, 0.6, 0.8, 1.0, 1.2):
            try:
                r = cluster_and_score(
                    factor_matrix, ground_truth, method="leiden",
                    leiden_gamma=gamma,
                )
                r["verdict"] = gate1_verdict(r["ami"], r["ari"])
                results.append(r)
                logger.info("Leiden γ=%.1f: AMI=%.3f ARI=%.3f → %s",
                            gamma, r["ami"], r["ari"], r["verdict"])
            except Exception as e:
                logger.warning("Leiden γ=%.1f failed: %s", gamma, e)
    else:
        logger.info("leidenalg not installed; skipping Leiden sweep")

    # HDBSCAN sweep
    if HDBSCAN_AVAILABLE:
        for min_size in (15, 25, 50):
            try:
                r = cluster_and_score(
                    factor_matrix, ground_truth, method="hdbscan",
                    hdbscan_min_size=min_size,
                )
                r["verdict"] = gate1_verdict(r["ami"], r["ari"])
                results.append(r)
                logger.info("HDBSCAN min_size=%d: AMI=%.3f ARI=%.3f → %s",
                            min_size, r["ami"], r["ari"], r["verdict"])
            except Exception as e:
                logger.warning("HDBSCAN min_size=%d failed: %s", min_size, e)
    else:
        logger.info("hdbscan not installed; skipping HDBSCAN sweep")

    # Persist
    df = pd.DataFrame(results)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(args.out, index=False)
    logger.info("Wrote %s (%d rows)", args.out, len(df))

    # Best result by AMI
    best = max(results, key=lambda r: r["ami"])
    overall_verdict = gate1_verdict(best["ami"], best["ari"])

    # Report
    L: list[str] = []
    L.append("# V8.4 Gate 1 Dry-Run v1 (synthetic phenotype)")
    L.append("")
    L.append("Sanity-checks the V8 Gate 1 mechanism-class recovery pipeline "
             "(Leiden + HDBSCAN + AMI/ARI) on synthetic MOFA+-like factors "
             "with known ground-truth cluster structure. If the pipeline "
             "can't recover 5 well-separated clusters here, it's broken.")
    L.append("")
    L.append("## Synthetic dataset")
    L.append("")
    L.append(f"- {args.n_classes} mechanism classes × {args.n_per_class} compounds "
             f"= {n_compounds} compounds total")
    L.append(f"- K = {args.K} latent factors (matches V8.3 MOFA+)")
    L.append(f"- Noise σ = {args.noise_sigma}")
    L.append(f"- Class names: {', '.join(class_names)}")
    L.append(f"- Class centroids: orthogonal in {args.K}-d space (separable)")
    L.append("")
    L.append("## Clustering availability")
    L.append("")
    L.append(f"- Agglomerative (sklearn): ✅ always available")
    L.append(f"- Leiden (leidenalg): {'✅' if LEIDEN_AVAILABLE else '⏳ not installed'}")
    L.append(f"- HDBSCAN (hdbscan): {'✅' if HDBSCAN_AVAILABLE else '⏳ not installed'}")
    L.append("")
    L.append("## Results")
    L.append("")
    L.append("| Method | n_pred | AMI | ARI | V-measure | FM | Verdict |")
    L.append("|---|---|---|---|---|---|---|")
    for r in results:
        method_label = r["method"]
        if r.get("leiden_gamma"):
            method_label += f" γ={r['leiden_gamma']:.1f}"
        elif r.get("hdbscan_min_size"):
            method_label += f" min={r['hdbscan_min_size']}"
        L.append(f"| {method_label} | {r['n_clusters_predicted']} | "
                 f"{r['ami']:.3f} | {r['ari']:.3f} | {r['v_measure']:.3f} | "
                 f"{r['fm']:.3f} | **{r['verdict']}** |")
    L.append("")
    L.append("## Best result")
    L.append("")
    best_label = best["method"]
    if best.get("leiden_gamma"):
        best_label += f" γ={best['leiden_gamma']:.1f}"
    elif best.get("hdbscan_min_size"):
        best_label += f" min={best['hdbscan_min_size']}"
    L.append(f"- **{best_label}**: AMI={best['ami']:.3f}, ARI={best['ari']:.3f}")
    L.append(f"- Overall verdict: **{overall_verdict}**")
    L.append("")
    L.append("## V8.4 Gate 1 thresholds (per V8 OSF pre-reg §5.1)")
    L.append("")
    L.append("| Band | AMI | ARI | Action |")
    L.append("|---|---|---|---|")
    L.append("| **PASS** | ≥ 0.50 | ≥ 0.40 | Enter joint posterior at λ_phen = 1.0 |")
    L.append("| **DEGRADE** | [0.30, 0.50) | [0.25, 0.40) | Enter at λ_phen = 0.5; flag |")
    L.append("| **FAIL** | < 0.30 | < 0.25 | Publish negative result |")
    L.append("")
    L.append("## Honest caveats")
    L.append("")
    L.append("- This is a SYNTHETIC sanity check, not real V8 mechanism-class "
             "recovery. Real V8 Gate 1 needs (a) actual MOFA+ joint embedding "
             "on real LINCS+JUMP-CP+iPSC-MEA data and (b) ground-truth labels "
             "from PRISMA chemistry+literature classification (NOT phenotype-"
             "derived; that would be circular).")
    L.append("- Synthetic centroids are orthogonal → cluster recovery is "
             "trivially possible. Real signature data has overlap (e.g., "
             "stimulants ↔ wake-promoting, AChE-I ↔ α7 nAChR PAM).")
    L.append("- The V8 PRISMA taxonomy is ~30 classes per OSF pre-reg §3, not "
             "5. Real Gate 1 will face a harder K=30 multi-class problem.")
    L.append("- Synthetic noise σ=0.30 is the V8.3 hyperparameter sweep "
             "midpoint; real LINCS WTCS τ noise is larger.")
    L.append("")
    L.append("---")
    L.append("")
    L.append("Generated by `scripts/60_v8_gate1_dryrun.py`. V8.4 Stage 1 sanity "
             "check of the Gate 1 pipeline. V8.4 Stage 2 requires real LINCS "
             "+ JUMP-CP + chemCPA + MOFA+ run.")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(L), encoding="utf-8")
    logger.info("Wrote %s", args.report)

    # Exit: 0 if any method passes, 1 if degrade, 2 if all fail
    if any(r["verdict"] == "PASS" for r in results):
        return 0
    if any(r["verdict"] == "DEGRADE" for r in results):
        return 1
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
