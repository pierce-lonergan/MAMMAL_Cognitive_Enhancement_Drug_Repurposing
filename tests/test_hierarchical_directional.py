"""C2 regression: sign-agnostic hierarchical slope + inverted-family routing (BUG_AUDIT_2026-06.md).

The NUTS slope prior was HalfNormal (slope >= 0), so the model could not fit the negative-rho
families it exists to rescue, and the slope SIGN was not estimable. The fix makes the slope
sign-agnostic (Normal) and routes confidently-negative-slope ("inverted but informative") families
to a separate exploratory output, out of the primary positive-direction shortlist.
"""
from __future__ import annotations

import numpy as np
import pytest


def test_classify_and_route_isolates_negative_slope():
    """Fast/deterministic: a confidently-negative-slope target is removed from the primary pooled
    rho and recorded in exploratory_negative_rho. (classify_and_route does not exist pre-fix.)"""
    from mammal_repurposing.calibration.hierarchical_bayes import classify_and_route
    pooled = {"POS": 0.6, "NEG": 0.7, "AMB": 0.4}
    sign_prob = {"POS": 0.99, "NEG": 0.02, "AMB": 0.6}   # NEG: P(beta>0)=0.02 -> inverted
    primary, exploratory, direction = classify_and_route(pooled, sign_prob)
    assert "NEG" not in primary and exploratory["NEG"] == 0.7
    assert direction["NEG"] == "negative"
    assert primary["POS"] == 0.6 and direction["POS"] == "positive"
    assert "AMB" in primary and direction["AMB"] == "ambiguous"   # unresolved sign stays primary


@pytest.mark.slow
def test_nuts_sign_agnostic_recovers_negative_slope_and_routes_to_exploratory():
    """The directional-rescue verification test. On a genuinely NEGATIVE-slope family the
    sign-agnostic model must estimate P(beta>0) < 0.5, label the direction "negative", and route
    the family to exploratory_negative_rho (excluded from pooled_rho). Under the old HalfNormal
    prior every beta draw was >= 0, so slope_sign_prob ~ 1.0 and direction would be "positive" ->
    these assertions FAIL on the pre-fix code and PASS post-fix."""
    pytest.importorskip("pymc")
    pytest.importorskip("numpyro")
    from mammal_repurposing.calibration.hierarchical_bayes import hierarchical_bayesian_nuts
    rng = np.random.default_rng(0)
    x1 = rng.uniform(5, 9, 25); y1 = -1.0 * x1 + 14 + rng.normal(0, 0.3, 25)   # slope ~ -1
    x2 = rng.uniform(5, 9, 25); y2 = -0.8 * x2 + 13 + rng.normal(0, 0.3, 25)   # slope ~ -0.8
    data = {"T_NEG1": (x1, y1), "T_NEG2": (x2, y2)}
    res = hierarchical_bayesian_nuts(
        "NEGFAM", data, n_chains=2, n_tune=400, n_draws=400, random_seed=0,
    )
    assert res.slope_sign_prob["T_NEG1"] < 0.5, "sign-agnostic prior must estimate a negative slope"
    assert res.direction["T_NEG1"] == "negative"
    assert "T_NEG1" in res.exploratory_negative_rho      # surfaced as an exploratory inverted family
    assert "T_NEG1" not in res.pooled_rho                # and excluded from the primary shortlist
