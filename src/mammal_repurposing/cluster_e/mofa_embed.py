"""V8.3 — MOFA+ joint Bayesian factor embedding across 7 views.

Per Argelaguet R, Arnol D, Bredikhin D, et al. 2020 *Genome Biol* 21:111
"MOFA+: a statistical framework for comprehensive integration of multi-modal
single-cell data" + V8 plan §B.3.

Architecture:
    K = 30 latent factors with ARD sparsity priors
    Views = {L1000_zscore, CP_CellProfiler, CP_DeepProfiler, CP_DINO,
             MEA_features, snRNA_pseudobulk, chemCPA_latent}
    Groups = {neural_lineage, non_neural_lineage, imputed}

Why MOFA+ over Deep CCA / Biolord (per V8 plan §B.3):
    (i) Bayesian sparse factor model with explicit per-modality variance
        decomposition — defends U2OS-to-brain transfer (factor 7 may
        explain 75% of JUMP-CP variance but ~0% of iPSC-MEA variance →
        flagged morphology-only, downweighted for cognition).
    (ii) Native missing-modality handling — most compounds have L1000 but
         not iPSC-MEA, and chemCPA-imputed signatures should be a clearly-
         labelled view.
    (iii) Scales to ~200K compound × ~50 factor on RTX 5070 + 32 GB RAM.

Graceful degradation:
    - mofapy2 missing → numpy SVD fallback on concatenated views for
      Stage 1 sanity (no per-view variance attribution; degraded utility
      but the pipeline doesn't break)
    - All views must be (n_compounds × n_features); per-view scaling
      handled internally
    - Missing values handled via mask propagation in mofapy2; SVD fallback
      uses mean-imputation

API:
    views = {
        "L1000": np.array (n_compounds, 977),
        "CP_DeepProfiler": np.array (n_compounds, 672),
        "MEA": np.array (n_compounds, 25),
        ...
    }
    groups = {"neural_lineage": [...uniprots], "non_neural_lineage": [...]}
    result = fit_mofa_plus(views, groups, K=30)
    # result.factor_matrix → (n_compounds, K)
    # result.per_view_variance → {view_name: array of K floats}
    # result.identifiability_diagnostics → dict
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)


# Optional MOFA+ backend -------------------------------------------------
try:
    import mofapy2  # noqa: F401
    MOFAPY2_AVAILABLE = True
except ImportError:
    MOFAPY2_AVAILABLE = False


# Default 7-view configuration per V8 plan §B.3
DEFAULT_VIEW_DIMS: dict[str, int] = {
    "L1000_zscore": 977,
    "CP_CellProfiler": 700,        # post feature-selection
    "CP_DeepProfiler": 672,        # CellPainting_CNN per Moshkov 2024
    "CP_DINO": 384,                # DINOv2 ViT-S/16
    "MEA_features": 25,            # MFR, WMFR, BR, NBF, NBD, ISI CV, IBI, sync...
    "snRNA_pseudobulk": 1000,      # scVI latent (typical Lopez 2018)
    "chemCPA_latent": 128,         # Hetzel 2022 latent dim
}


@dataclass
class MofaConfig:
    """Configuration for MOFA+ joint factor model."""
    K: int = 30                              # number of latent factors
    n_iter: int = 1000
    convergence_threshold: float = 0.0001
    ard_per_factor: bool = True              # ARD sparsity prior
    spikeslab_weights: bool = True
    backend: str = "auto"                    # 'auto' picks mofapy2 if available
    sparsify_groups: bool = True             # per-group ARD
    seed: int = 42


@dataclass
class MofaResult:
    """Output of fit_mofa_plus()."""
    factor_matrix: np.ndarray                # (n_compounds, K)
    per_view_variance: dict[str, np.ndarray] = field(default_factory=dict)
    per_factor_per_view_variance: dict[str, dict[str, float]] = field(default_factory=dict)
    weights_per_view: dict[str, np.ndarray] = field(default_factory=dict)
    K: int = 30
    n_compounds: int = 0
    n_views: int = 0
    n_iter_used: int = 0
    converged: bool = False
    method: str = "mofapy2"
    note: str = ""


def fit_mofa_plus(
    views: dict[str, np.ndarray],
    groups: dict[str, list[int]] | None = None,
    cfg: MofaConfig | None = None,
) -> MofaResult:
    """Fit MOFA+ K-factor joint model across N views.

    Args:
        views: {view_name: (n_compounds, n_features)} dict; all matrices must
            share the n_compounds axis (NaN allowed for missing observations)
        groups: optional {group_name: list_of_compound_indices} for per-group
            ARD; default = single group
        cfg: MofaConfig

    Returns MofaResult with factor matrix + per-view variance attribution.
    """
    cfg = cfg or MofaConfig()
    if not views:
        raise ValueError("Cannot fit MOFA+ on empty views dict")

    # Sanity: all views share n_compounds
    n_compounds_list = [v.shape[0] for v in views.values()]
    if len(set(n_compounds_list)) != 1:
        raise ValueError(
            f"All views must share n_compounds; got {dict(zip(views.keys(), n_compounds_list))}"
        )
    n_compounds = n_compounds_list[0]

    backend = cfg.backend
    if backend == "auto":
        backend = "mofapy2" if MOFAPY2_AVAILABLE else "numpy_svd"
    if backend == "mofapy2" and not MOFAPY2_AVAILABLE:
        logger.warning("MOFA+ requested but mofapy2 missing; falling back to numpy SVD")
        backend = "numpy_svd"

    if backend == "mofapy2":
        return _fit_mofapy2(views, groups, cfg, n_compounds)
    return _fit_svd_fallback(views, cfg, n_compounds)


def _fit_svd_fallback(
    views: dict[str, np.ndarray],
    cfg: MofaConfig,
    n_compounds: int,
) -> MofaResult:
    """Numpy SVD fallback: concatenate views (mean-imputed) + truncated SVD.

    Loses per-view variance attribution (everything attributed to "joint").
    Used when mofapy2 is unavailable so the V8 pipeline doesn't break.
    """
    # Mean-impute missing values per view, then z-score per feature
    stacked_cols: list[np.ndarray] = []
    view_widths: dict[str, tuple[int, int]] = {}    # {view: (start, end)}
    col_idx = 0
    for name, X in views.items():
        Xv = X.copy()
        # Mean-impute NaNs per column
        col_means = np.nanmean(Xv, axis=0)
        col_means = np.where(np.isfinite(col_means), col_means, 0.0)
        for c in range(Xv.shape[1]):
            mask = ~np.isfinite(Xv[:, c])
            Xv[mask, c] = col_means[c]
        # Z-score per feature
        std = Xv.std(axis=0, keepdims=True)
        std[std == 0] = 1.0
        Xv = (Xv - Xv.mean(axis=0, keepdims=True)) / std
        stacked_cols.append(Xv)
        view_widths[name] = (col_idx, col_idx + Xv.shape[1])
        col_idx += Xv.shape[1]
    X_full = np.hstack(stacked_cols)    # (n_compounds, total_features)

    # Truncated SVD via numpy
    K = min(cfg.K, n_compounds, X_full.shape[1])
    Xc = X_full - X_full.mean(axis=0, keepdims=True)
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    factor_matrix = U[:, :K] * S[:K]      # (n_compounds, K) factor scores
    V = Vt[:K, :]                          # (K, total_features) loadings

    # Per-view variance attribution: for each factor k, sum squared loadings
    # within each view's column slice
    per_view_var: dict[str, np.ndarray] = {}
    per_factor_per_view: dict[str, dict[str, float]] = {}
    factor_total_var = (S[:K] ** 2) / Xc.shape[0]
    for name, (start, end) in view_widths.items():
        view_loadings = V[:, start:end]   # (K, view_features)
        # Variance explained per factor in this view = sum(loadings²) × eigenval / total_var
        view_var_per_factor = np.sum(view_loadings ** 2, axis=1) * factor_total_var
        per_view_var[name] = view_var_per_factor
        per_factor_per_view[name] = {f"factor_{k}": float(view_var_per_factor[k])
                                     for k in range(K)}

    return MofaResult(
        factor_matrix=factor_matrix,
        per_view_variance=per_view_var,
        per_factor_per_view_variance=per_factor_per_view,
        weights_per_view={name: V[:, start:end].T   # (view_features, K) loadings
                          for name, (start, end) in view_widths.items()},
        K=K, n_compounds=n_compounds, n_views=len(views),
        n_iter_used=1,
        converged=True,
        method="numpy_svd_fallback",
        note=("mofapy2 unavailable. Used truncated-SVD on z-scored concatenated "
              "views with mean-imputation. Lacks the ARD sparsity + Bayesian "
              "per-group decomposition; for production V8.3 use "
              "`pip install mofapy2` to activate the full model."),
    )


def _fit_mofapy2(
    views: dict[str, np.ndarray],
    groups: dict[str, list[int]] | None,
    cfg: MofaConfig,
    n_compounds: int,
) -> MofaResult:
    """Full MOFA+ fit via mofapy2 (Argelaguet 2020).

    Requires mofapy2 installed. Returns MofaResult with Bayesian variance
    attribution per view per factor (the load-bearing diagnostic for V8.3).
    """
    if not MOFAPY2_AVAILABLE:
        raise ImportError("mofapy2 not available; use SVD fallback.")
    import mofapy2.run.entry_point as me

    ent = me.entry_point()
    ent.set_data_options(scale_views=True, center_groups=True)
    # Build per-view data matrices in mofapy2's expected shape: list of lists
    # [view_idx][group_idx] = (n_compounds_in_group, n_features) matrix
    if groups is None:
        groups = {"single_group": list(range(n_compounds))}
    group_names = list(groups.keys())
    view_names = list(views.keys())

    data_lists: list[list[np.ndarray]] = []
    for vn in view_names:
        per_group: list[np.ndarray] = []
        for gn in group_names:
            idx = groups[gn]
            per_group.append(views[vn][idx])
        data_lists.append(per_group)
    ent.set_data_matrix(data_lists,
                        groups_names=group_names,
                        views_names=view_names)
    ent.set_model_options(factors=cfg.K,
                          ard_factors=cfg.ard_per_factor,
                          spikeslab_weights=cfg.spikeslab_weights,
                          ard_weights=cfg.sparsify_groups)
    ent.set_train_options(iter=cfg.n_iter,
                          convergence_mode="medium",
                          dropR2=0.0,
                          seed=cfg.seed)
    ent.build()
    ent.run()

    # Extract factor matrix
    Z = ent.model.getExpectations("Z", expand=True)
    if isinstance(Z, dict):    # per-group dict
        factor_matrix = np.vstack(list(Z.values()))
    else:
        factor_matrix = np.asarray(Z)

    # Per-view per-factor variance via getVarianceExplained
    var_decomp = ent.model.calculate_variance_explained()
    per_view_var: dict[str, np.ndarray] = {}
    per_factor_per_view: dict[str, dict[str, float]] = {}
    for vn in view_names:
        # var_decomp shape: (n_groups, n_factors) per view
        try:
            r2 = var_decomp["r2_per_factor"][vn]
            arr = np.asarray(r2).mean(axis=0)    # mean across groups
            per_view_var[vn] = arr
            per_factor_per_view[vn] = {f"factor_{k}": float(arr[k])
                                       for k in range(len(arr))}
        except Exception as e:
            logger.warning("Could not extract var for view %s: %s", vn, e)

    # Weights per view
    W = ent.model.getExpectations("W")
    weights_per_view: dict[str, np.ndarray] = {}
    for i, vn in enumerate(view_names):
        try:
            weights_per_view[vn] = np.asarray(W[i])
        except Exception:
            pass

    converged = ent.model.getTrainingStats().get("converged", False)

    return MofaResult(
        factor_matrix=factor_matrix,
        per_view_variance=per_view_var,
        per_factor_per_view_variance=per_factor_per_view,
        weights_per_view=weights_per_view,
        K=cfg.K, n_compounds=n_compounds, n_views=len(views),
        n_iter_used=ent.model.getTrainingStats().get("iter", cfg.n_iter),
        converged=converged,
        method="mofapy2",
        note=f"MOFA+ converged={converged} iter={cfg.n_iter}",
    )


def variance_attribution_table(
    result: MofaResult,
    top_k_factors: int = 10,
) -> dict:
    """Compact dict for V8.3 per-factor per-view variance attribution.

    Useful for the V8.3 paper figure + the V8.6 Discussion section defending
    U2OS-to-brain transfer ("factor 7 = morphology-only; downweighted").
    """
    out: dict = {
        "method": result.method,
        "K": result.K,
        "n_views": result.n_views,
        "top_factors_per_view": {},
    }
    for view_name, var_per_factor in result.per_view_variance.items():
        ordering = np.argsort(-var_per_factor)[:top_k_factors]
        out["top_factors_per_view"][view_name] = [
            {"factor": int(f), "variance_explained": float(var_per_factor[f])}
            for f in ordering
        ]
    return out


def cosine_similarity_matrix(
    factor_matrix: np.ndarray,
) -> np.ndarray:
    """Compound × compound cosine similarity on MOFA+ factor vectors.

    Used by V8.5 joint posterior as the phenotypic-match feature φ_c.
    """
    norms = np.linalg.norm(factor_matrix, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    normed = factor_matrix / norms
    return normed @ normed.T


def cognition_centroid_similarity(
    factor_matrix: np.ndarray,
    centroid_indices: dict[str, list[int]],
) -> dict[str, np.ndarray]:
    """For each cognition centroid (cholinergic / catecholaminergic / etc.),
    compute per-compound cosine to the centroid in MOFA+ factor space.

    Returns {centroid_name: (n_compounds,) cosine array}.
    """
    out: dict[str, np.ndarray] = {}
    norms = np.linalg.norm(factor_matrix, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    normed = factor_matrix / norms
    for name, idx_list in centroid_indices.items():
        if not idx_list:
            out[name] = np.zeros(factor_matrix.shape[0])
            continue
        centroid = factor_matrix[idx_list].mean(axis=0)
        c_norm = centroid / max(np.linalg.norm(centroid), 1e-12)
        out[name] = normed @ c_norm
    return out


def availability() -> dict[str, object]:
    """Probe MOFA+ availability."""
    return {
        "available": True,            # SVD fallback always works
        "mofapy2_backend": MOFAPY2_AVAILABLE,
        "stub_mode_works_without_mofapy2": True,
        "default_K": 30,
        "default_view_dims": dict(DEFAULT_VIEW_DIMS),
        "n_views_default": len(DEFAULT_VIEW_DIMS),
    }
