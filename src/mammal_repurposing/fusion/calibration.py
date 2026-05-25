"""Optional Platt scaling — used only if CombSUM-style fusion is requested.

RRF is rank-based and immune to score-distribution heterogeneity, so for the
default v2 path this module is unused. Keep it for two reasons:
    1. CombSUM is an obvious "what about a simple weighted sum?" comparison
       that reviewers will ask about. We can answer "we calibrated first."
    2. LambdaMART's tree models are scale-invariant, but Platt-scaled features
       are useful for downstream linear models / regression heads.

Reference: Platt JC. "Probabilistic Outputs for Support Vector Machines and
Comparisons to Regularized Likelihood Methods." Adv. Large Margin Classifiers, 1999.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from dataclasses import dataclass


@dataclass
class PlattScaler:
    """Two-parameter sigmoid: P(positive | score) = 1 / (1 + exp(A*score + B))."""

    A: float
    B: float

    def transform(self, scores: np.ndarray | pd.Series) -> np.ndarray | pd.Series:
        if isinstance(scores, pd.Series):
            return 1.0 / (1.0 + np.exp(self.A * scores.values + self.B))
        return 1.0 / (1.0 + np.exp(self.A * np.asarray(scores) + self.B))


def fit_platt(scores: np.ndarray, labels: np.ndarray) -> PlattScaler:
    """Fit a Platt scaler via Newton-Raphson on the log-likelihood.

    Tiny implementation; sklearn would also work but adds a heavy dep.
    """
    scores = np.asarray(scores, dtype=float)
    labels = np.asarray(labels, dtype=int)

    n_pos = (labels == 1).sum()
    n_neg = (labels == 0).sum()
    if n_pos == 0 or n_neg == 0:
        raise ValueError("Platt fit requires both positive and negative labels.")

    # Targets per Platt (1999) to avoid overfitting at boundaries
    hi = (n_pos + 1.0) / (n_pos + 2.0)
    lo = 1.0 / (n_neg + 2.0)
    t = np.where(labels == 1, hi, lo)

    A, B = 0.0, np.log((n_neg + 1.0) / (n_pos + 1.0))
    for _ in range(100):
        fApB = scores * A + B
        p = 1.0 / (1.0 + np.exp(fApB))
        d_A = float(((p - t) * scores).sum())
        d_B = float((p - t).sum())
        h_AA = float((p * (1 - p) * scores * scores).sum()) + 1e-12
        h_BB = float((p * (1 - p)).sum()) + 1e-12
        h_AB = float((p * (1 - p) * scores).sum())
        det = h_AA * h_BB - h_AB * h_AB + 1e-12
        dA = (h_BB * d_A - h_AB * d_B) / det
        dB = (-h_AB * d_A + h_AA * d_B) / det
        A -= dA
        B -= dB
        if abs(dA) < 1e-6 and abs(dB) < 1e-6:
            break

    return PlattScaler(A=A, B=B)
