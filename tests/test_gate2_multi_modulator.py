"""Sprint 2.2 — Tests for gate_2_multi_modulator_spearman.

Validates the multi-modulator Gate 2 evaluator against:
  - Synthetic positive control (perfect rank agreement → ρ ≈ 1.0)
  - Synthetic negative control (anti-correlated → ρ ≈ -1.0 → FAIL)
  - Synthetic noise (independent → ρ ≈ 0.0 → FAIL or DEGRADE)
  - Real modulator_anchors.parquet (Sprint 2.1) when present
  - Edge cases: too few targets, single modulator per target, etc.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mammal_repurposing.cluster_d.validation_gates import (
    gate_2_multi_modulator_spearman,
)


# ---------------------------------------------------------------------------
# Synthetic fixtures
# ---------------------------------------------------------------------------

def _synth_rows(
    n_targets: int = 20,
    n_modulators_per_target: int = 2,
    relationship: str = "positive",
    noise: float = 0.05,
    rng_seed: int = 42,
) -> tuple[dict[str, float], list[dict]]:
    """Build a (theta_mean, rows) pair where each target has n_modulators.

    relationship:
      - "positive": pooled_g ≈ θ̄ + small noise (perfect rank agreement)
      - "negative": pooled_g ≈ -θ̄ + small noise
      - "noise":    pooled_g is random, independent of θ̄
    """
    rng = np.random.default_rng(rng_seed)
    # Generate target θ̄ spread evenly in [-0.5, +0.5]
    thetas = np.linspace(-0.5, 0.5, n_targets)
    rng.shuffle(thetas)
    theta_mean = {f"U{i:03d}": float(thetas[i]) for i in range(n_targets)}

    rows = []
    for i in range(n_targets):
        u = f"U{i:03d}"
        for m in range(n_modulators_per_target):
            if relationship == "positive":
                g = thetas[i] * 0.6 + rng.normal(0, noise)
            elif relationship == "negative":
                g = -thetas[i] * 0.6 + rng.normal(0, noise)
            elif relationship == "noise":
                g = rng.normal(0, 0.3)
            else:
                raise ValueError(relationship)
            g = float(np.clip(g, -1.0, 1.0))
            rows.append({
                "target_uniprot": u,
                "target_gene": f"G{i}",
                "compound": f"cmpd_{i}_{m}",
                "mechanism": "test",
                "pooled_g": g,
                "CI_lo": g - 0.10,
                "CI_hi": g + 0.10,
                "k": 1,
                "endpoint": "synth",
                "citation_doi": "n/a",
                "population": "synth",
                "notes": "",
            })
    return theta_mean, rows


# ---------------------------------------------------------------------------
# Positive / negative / noise sanity checks
# ---------------------------------------------------------------------------

def test_positive_control_passes():
    """Perfect rank agreement → ρ ≈ +1 → PASS."""
    theta_mean, rows = _synth_rows(relationship="positive", n_targets=30,
                                    n_modulators_per_target=2)
    res = gate_2_multi_modulator_spearman(theta_mean, rows, threshold=0.30)
    assert res.pass_status == "PASS", res.detail
    assert res.metric_value > 0.80, f"ρ should be > 0.80, got {res.metric_value:.3f}"


def test_negative_control_fails():
    """Anti-correlated → ρ ≈ -1 → FAIL."""
    theta_mean, rows = _synth_rows(relationship="negative", n_targets=30,
                                    n_modulators_per_target=2)
    res = gate_2_multi_modulator_spearman(theta_mean, rows, threshold=0.30)
    assert res.pass_status == "FAIL", res.detail
    assert res.metric_value < -0.50, f"ρ should be < -0.50, got {res.metric_value:.3f}"


def test_noise_does_not_pass():
    """Independent noise → ρ ≈ 0 → DEGRADE or FAIL."""
    theta_mean, rows = _synth_rows(relationship="noise", n_targets=30,
                                    n_modulators_per_target=2)
    res = gate_2_multi_modulator_spearman(theta_mean, rows, threshold=0.30)
    assert res.pass_status in {"DEGRADE", "FAIL"}, res.detail


# ---------------------------------------------------------------------------
# Aggregation strategies
# ---------------------------------------------------------------------------

class TestAggregationStrategies:

    def setup_method(self):
        self.theta_mean = {"U001": 0.5, "U002": -0.3, "U003": 0.1}
        self.rows = [
            # U001: high effect compounds
            {"target_uniprot": "U001", "compound": "c1", "pooled_g": 0.4,
             "CI_lo": 0.3, "CI_hi": 0.5, "k": 5},
            {"target_uniprot": "U001", "compound": "c2", "pooled_g": 0.6,
             "CI_lo": 0.5, "CI_hi": 0.7, "k": 3},
            # U002: low-effect compounds
            {"target_uniprot": "U002", "compound": "c3", "pooled_g": -0.2,
             "CI_lo": -0.3, "CI_hi": -0.1, "k": 2},
            {"target_uniprot": "U002", "compound": "c4", "pooled_g": -0.4,
             "CI_lo": -0.5, "CI_hi": -0.3, "k": 4},
            # U003: mixed
            {"target_uniprot": "U003", "compound": "c5", "pooled_g": 0.0,
             "CI_lo": -0.1, "CI_hi": 0.1, "k": 1},
            {"target_uniprot": "U003", "compound": "c6", "pooled_g": 0.2,
             "CI_lo": 0.1, "CI_hi": 0.3, "k": 1},
        ]

    def test_mean_aggregation(self):
        res = gate_2_multi_modulator_spearman(
            self.theta_mean, self.rows, aggregation="mean", threshold=0.30,
        )
        # 3 pairs: (0.5, 0.5), (-0.3, -0.3), (0.1, 0.1) → perfect rank
        # Note: too few pairs for INSUFFICIENT_DATA branch (need ≥5)
        # so this should return INSUFFICIENT_DATA
        assert res.pass_status == "INSUFFICIENT_DATA"

    def test_max_aggregation_picks_best_modulator(self):
        # Add 2 more targets to get past the n≥5 floor
        for i in range(4, 8):
            self.rows.append({"target_uniprot": f"U00{i}", "compound": f"c{i}",
                              "pooled_g": 0.1 * i, "CI_lo": 0.1 * i - 0.05,
                              "CI_hi": 0.1 * i + 0.05, "k": 1})
            self.theta_mean[f"U00{i}"] = 0.1 * i
        res = gate_2_multi_modulator_spearman(
            self.theta_mean, self.rows, aggregation="max", threshold=0.30,
        )
        assert res.pass_status in {"PASS", "DEGRADE"}, res.detail

    def test_unknown_aggregation_raises(self):
        with pytest.raises(ValueError, match="Unknown aggregation"):
            gate_2_multi_modulator_spearman(
                self.theta_mean, self.rows, aggregation="bogus",
            )


# ---------------------------------------------------------------------------
# DataFrame input acceptance
# ---------------------------------------------------------------------------

def test_dataframe_input_accepted():
    theta_mean, rows = _synth_rows(relationship="positive", n_targets=20)
    df = pd.DataFrame(rows)
    res = gate_2_multi_modulator_spearman(theta_mean, df, threshold=0.30)
    assert res.pass_status == "PASS"


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

def test_insufficient_data_returned_when_too_few_targets():
    theta_mean = {"U001": 0.5, "U002": -0.5, "U003": 0.0}
    rows = [
        {"target_uniprot": "U001", "compound": "c1", "pooled_g": 0.4,
         "CI_lo": 0.3, "CI_hi": 0.5, "k": 1},
        {"target_uniprot": "U002", "compound": "c2", "pooled_g": -0.4,
         "CI_lo": -0.5, "CI_hi": -0.3, "k": 1},
    ]
    res = gate_2_multi_modulator_spearman(theta_mean, rows)
    assert res.pass_status == "INSUFFICIENT_DATA"


def test_targets_not_in_theta_are_dropped():
    """If a modulator's target isn't in theta_mean, that row is silently dropped."""
    theta_mean, rows = _synth_rows(relationship="positive", n_targets=10)
    rows.append({"target_uniprot": "MISSING_UNIPROT", "compound": "x",
                 "pooled_g": 0.99, "CI_lo": 0.9, "CI_hi": 1.0, "k": 1})
    res = gate_2_multi_modulator_spearman(theta_mean, rows, threshold=0.30)
    # MISSING is silently dropped — should still PASS on the 10 real targets
    assert res.pass_status == "PASS"


# ---------------------------------------------------------------------------
# Real modulator anchors integration (Sprint 2.1 output)
# ---------------------------------------------------------------------------

ROOT = Path(__file__).resolve().parents[1]
MODULATOR_PARQUET = ROOT / "data" / "interim" / "modulator_anchors.parquet"
V6B5_POSTERIOR = ROOT / "data" / "results" / "v2" / "cluster_d_posterior_expanded_v2_mh8_ta99.parquet"
V6B_HEADLINE_POSTERIOR = ROOT / "data" / "results" / "v2" / "cluster_d_posterior_v1.parquet"


@pytest.mark.skipif(
    not MODULATOR_PARQUET.exists(),
    reason=f"Modulator parquet missing: {MODULATOR_PARQUET}",
)
def test_real_modulator_anchors_load_clean():
    df = pd.read_parquet(MODULATOR_PARQUET)
    assert len(df) >= 60, f"Expected >=60 rows, got {len(df)}"
    assert df["target_uniprot"].nunique() >= 25
    # Phase III nulls present
    n_nulls = (df["pooled_g"].abs() <= 0.06).sum()
    assert n_nulls >= 15, f"Expected >=15 Phase III nulls, got {n_nulls}"


@pytest.mark.skipif(
    not MODULATOR_PARQUET.exists() or not V6B5_POSTERIOR.exists(),
    reason="Requires both modulator parquet AND V6.B.5 posterior",
)
def test_gate2_runs_on_real_v6b5_posterior(capsys):
    """End-to-end smoke: load real V6.B.5 posterior + real modulator anchors
    + run Gate 2. Records the actual rho as a sanity check — doesn't enforce
    PASS because the small n (~32 targets surviving panel intersection) and
    noisy synthetic AHBA on 169/191 targets may legitimately yield FAIL.

    This test documents the production Gate 2 number for the V6.B paper's
    Results section. Sprint 2.2 deliverable.
    """
    posterior_df = pd.read_parquet(V6B5_POSTERIOR)
    theta_mean = dict(zip(posterior_df["target_uniprot"], posterior_df["theta_mean"]))
    modulators = pd.read_parquet(MODULATOR_PARQUET)

    results: dict[str, tuple[str, float]] = {}
    for agg in ["mean", "median", "max", "weighted_mean"]:
        res = gate_2_multi_modulator_spearman(theta_mean, modulators,
                                               aggregation=agg, threshold=0.30)
        assert res.pass_status in {"PASS", "DEGRADE", "FAIL", "INSUFFICIENT_DATA"}
        assert np.isfinite(res.metric_value) or res.pass_status == "INSUFFICIENT_DATA"
        results[agg] = (res.pass_status, res.metric_value)

    # Sanity: at least one aggregation should produce a finite rho.
    assert any(np.isfinite(v) for _, v in results.values())
    # Write production-number summary to capsys for human reading.
    print("V6.B.5 production Gate 2 results:")
    for agg, (status, rho) in results.items():
        print(f"  agg={agg}: status={status}, rho={rho:+.3f}")


@pytest.mark.skipif(
    not MODULATOR_PARQUET.exists() or not V6B_HEADLINE_POSTERIOR.exists(),
    reason="Requires both modulator parquet AND V6.B headline posterior",
)
def test_gate2_runs_on_v6b_headline_posterior(capsys):
    """Gate 2 on the V6.B 22-target headline panel (real AHBA throughout).

    Hypothesis: Gate 2 should perform BETTER on the headline panel than on
    V6.B.5 because all 22 anchor AHBA values are real, not synthetic. If
    Gate 2 is FAIL on both panels, the issue is NOT the synthetic AHBA — it
    is genuine model-vs-clinic disagreement (likely driven by Phase III
    nulls anchored at clinically-validated high-affinity targets).
    """
    posterior_df = pd.read_parquet(V6B_HEADLINE_POSTERIOR)
    theta_mean = dict(zip(posterior_df["target_uniprot"], posterior_df["theta_mean"]))
    modulators = pd.read_parquet(MODULATOR_PARQUET)

    results: dict[str, tuple[str, float]] = {}
    for agg in ["mean", "median", "max", "weighted_mean"]:
        res = gate_2_multi_modulator_spearman(theta_mean, modulators,
                                               aggregation=agg, threshold=0.30,
                                               min_modulators_per_target=1)
        results[agg] = (res.pass_status, res.metric_value)

    assert any(np.isfinite(v) for _, v in results.values())
    print("V6.B HEADLINE Gate 2 results (real AHBA on all 22 targets):")
    for agg, (status, rho) in results.items():
        print(f"  agg={agg}: status={status}, rho={rho:+.3f}")
