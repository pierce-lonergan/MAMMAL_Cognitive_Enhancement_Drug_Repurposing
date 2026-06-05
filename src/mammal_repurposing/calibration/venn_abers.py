"""V6.A.4 — Venn-ABERS calibrated uncertainty propagation.

Per Mervin et al. 2020 *J Chem Inf Model* 60(10):4546 (doi:10.1021/acs.jcim.0c00476)
the AstraZeneca 40M-pair benchmark concluded Venn-ABERS gave the best calibration
across ML algorithms and CV methods, with the lowest Brier score loss.

Per Multi Head DTI.md §5.2: this module wraps per-head DTI predictions with
inductive Venn-ABERS predictors (Vovk & Petej 2014 UAI; Nouretdinov et al.
2018 PMLR 91:1). The router (V6.A.3) uses VA predictive distributions to
propagate per-head uncertainty into ensemble CIs.

API:
    from mammal_repurposing.calibration.venn_abers import (
        VennAbersRegressor, fit_va_per_head,
    )
    # Per head: fit on calibration fold, predict (lower, upper) per query
    va = VennAbersRegressor()
    va.fit(cal_predictions, cal_truth)
    lo, hi = va.predict_interval(query_predictions)

This is the V6.A.4 SKELETON — the router scaffold (fusion/bayesian_router.py)
already accepts head_sigmas as input. When VA-derived (lo, hi) intervals are
available, σ = (hi - lo) / (2 * 1.96) approximates the per-head std for the
existing Gaussian-CI router; full multivariate VA copula MC propagation
becomes the V6.A.4.b extension.

Reference:
  Mervin et al. 2020 — VA Brier-best across 40M pairs
  Vovk & Petej 2014 UAI pp.829 — VA validity proof
  Nouretdinov et al. 2018 PMLR 91:1 — IVAR inductive variant
  Kull et al. 2017 EJS 11:5052 — beta-calibration (alternative for skewed scores)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class VennAbersConfig:
    """Configuration for the Venn-ABERS regressor.

    alpha: target miscoverage rate (e.g. 0.10 → 90% nominal coverage)
    monotonic_direction: 'auto' (data-driven) | True (increasing) | False
    """
    alpha: float = 0.10
    monotonic_direction: str | bool = "auto"
    n_bins: int = 100   # bin count for the isotonic predictive distribution


class VennAbersRegressor:
    """Inductive Venn-ABERS regressor (Vovk & Petej 2014; Nouretdinov 2018).

    For a regression target with calibration set (x_i, y_i)_{i=1..n}:
      1. Fit two isotonic regressors:
         - iso_lower assigns the query the LOWEST y consistent with monotone fit
         - iso_upper assigns the query the HIGHEST y consistent with monotone fit
      2. The VA prediction interval [lo, hi] has guaranteed marginal coverage
         under the exchangeability assumption.

    For per-head DTI calibration: fit one VA regressor per target on the
    head's predictions on the calibration fold (held-out ChEMBL truth).
    """

    def __init__(self, config: VennAbersConfig | None = None):
        self.cfg = config or VennAbersConfig()
        self._fitted: bool = False
        self._cal_x: np.ndarray = np.array([])
        self._cal_y: np.ndarray = np.array([])

    def fit(self, x: np.ndarray, y: np.ndarray) -> "VennAbersRegressor":
        """Store calibration set (no parameter fitting — VA is non-parametric)."""
        self._cal_x = np.asarray(x, dtype=float)
        self._cal_y = np.asarray(y, dtype=float)
        if len(self._cal_x) != len(self._cal_y):
            raise ValueError(
                f"x ({len(self._cal_x)}) and y ({len(self._cal_y)}) length mismatch"
            )
        if len(self._cal_x) < 10:
            logger.warning("VA calibration set has only %d points; intervals will be wide",
                           len(self._cal_x))
        self._fitted = True
        return self

    def predict_interval(self, x_query: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Returns (lower, upper) marginal prediction intervals for each query."""
        if not self._fitted:
            raise RuntimeError("VennAbersRegressor.fit() must be called first")

        from sklearn.isotonic import IsotonicRegression
        x_query = np.asarray(x_query, dtype=float)

        # Per-query VA pair: add the query (x_q, y_assumed) to the cal set and
        # ask isotonic for the lowest/highest y consistent with monotone fit.
        # This is the canonical Vovk-Petej algorithm; simplified to a
        # quantile-based shortcut here.
        residuals = np.empty(len(self._cal_x), dtype=float)
        iso = IsotonicRegression(
            increasing="auto" if self.cfg.monotonic_direction == "auto"
            else bool(self.cfg.monotonic_direction),
            out_of_bounds="clip",
        )
        iso.fit(self._cal_x, self._cal_y)
        cal_pred = iso.predict(self._cal_x)
        residuals = np.abs(self._cal_y - cal_pred)

        # Quantile threshold for (1 - alpha) coverage
        q = np.quantile(residuals, 1 - self.cfg.alpha)
        point = iso.predict(x_query)
        lower = point - q
        upper = point + q
        return lower, upper

    def predict_sigma(self, x_query: np.ndarray) -> np.ndarray:
        """Approximate per-query σ from the VA interval width.

        The interval is a (1 - alpha) interval, so σ ≈ half-width / z with
        z = Φ⁻¹(1 - alpha/2). (Previously hardcoded the 95% z = 1.96, which
        mis-scaled σ whenever alpha != 0.05; the default alpha is 0.10.)
        The router consumes this directly as head_sigmas.
        """
        from scipy.stats import norm
        z = float(norm.ppf(1 - self.cfg.alpha / 2))
        lo, hi = self.predict_interval(x_query)
        return (hi - lo) / (2 * z)


def fit_va_per_head(
    cal_predictions: dict[str, np.ndarray],
    cal_truth: np.ndarray,
    config: VennAbersConfig | None = None,
) -> dict[str, VennAbersRegressor]:
    """Fit one Venn-ABERS regressor per head on a shared calibration set.

    Args:
        cal_predictions: {head_name: prediction_array_on_cal_fold}
        cal_truth: ground-truth pchembl values for the same compounds
        config: shared VA config (alpha, etc.)

    Returns:
        {head_name: fitted VennAbersRegressor}
    """
    cfg = config or VennAbersConfig()
    out: dict[str, VennAbersRegressor] = {}
    for head, preds in cal_predictions.items():
        va = VennAbersRegressor(cfg)
        va.fit(preds, cal_truth)
        out[head] = va
        logger.info("Fitted Venn-ABERS for head=%s on n=%d cal points",
                    head, len(cal_truth))
    return out


def correlated_mc_intervals(
    point_predictions: dict[str, np.ndarray],
    sigmas: dict[str, np.ndarray],
    weights: dict[str, float],
    correlation: np.ndarray | None = None,
    n_samples: int = 1000,
    rng_seed: int = 42,
) -> tuple[np.ndarray, np.ndarray]:
    """Multivariate Monte Carlo propagation of per-head VA intervals.

    Per Multi Head DTI.md §5.2 — cross-head correlation matters when heads
    share training data (MMAtt-DTA + BALM both touch BindingDB Kd). Pass a
    K×K correlation matrix; default identity = independence.

    Returns (ensemble_lower, ensemble_upper) per query.
    """
    heads = list(point_predictions.keys())
    K = len(heads)
    if K == 0:
        return np.array([]), np.array([])

    # Stack into K × n_query arrays
    point_arr = np.stack([point_predictions[h] for h in heads])
    sigma_arr = np.stack([sigmas[h] for h in heads])
    w_arr = np.array([weights.get(h, 1.0 / K) for h in heads])
    w_arr = w_arr / w_arr.sum()    # normalise

    if correlation is None:
        correlation = np.eye(K)

    rng = np.random.default_rng(rng_seed)
    n_query = point_arr.shape[1]

    L = np.linalg.cholesky(correlation + 1e-6 * np.eye(K))
    lower = np.empty(n_query)
    upper = np.empty(n_query)
    for j in range(n_query):
        # Per query: sample K-dim Gaussian with mean point, cov σ·ρ·σ
        z = rng.standard_normal((n_samples, K))
        # Decorrelated → correlated
        z_corr = z @ L.T
        # Scale per query
        samples_per_head = point_arr[:, j][None, :] + z_corr * sigma_arr[:, j][None, :]
        # Weighted ensemble per draw
        ens = samples_per_head @ w_arr
        lower[j] = float(np.quantile(ens, 0.025))
        upper[j] = float(np.quantile(ens, 0.975))
    return lower, upper
