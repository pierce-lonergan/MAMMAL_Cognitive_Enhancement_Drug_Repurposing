"""B4 regression: the clinician-dossier fallback credible interval must be a TRUE two-sided 90% CrI
(z = 1.6449), matching its "90% CrI" label and the V7 pre-registered 90% commitment — not the
two-sided 80% interval (z = 1.2816) it previously used (BUG_AUDIT_2026-06.md, B4).
"""
from __future__ import annotations

import pandas as pd

from mammal_repurposing.reporting.clinician_dossier import build_dossier
from mammal_repurposing.validation.disease_reframe import DiseaseClassPrior


def test_dossier_fallback_ci_is_true_90pct():
    prior = DiseaseClassPrior(
        mechanism_class="testclass", mean=0.30, sd=0.10, n_drugs=3, k_total=5,
        n_success=2, n_fail=1, success_rate=0.667, drugs=["a", "b", "c"],
    )
    card = build_dossier(
        "tc_compound", "AD", ledger=pd.DataFrame(),
        disease_priors={"testclass": prior}, mechanism_class="testclass",
    )
    half = card.g_ci_hi - card.g
    # true two-sided 90% half-width = 1.6449 * sd = 0.16449 (fails-before: old code gave 0.12816)
    assert abs(half - 1.6449 * 0.10) < 1e-3, f"expected true-90% half-width 0.16449, got {half:.5f}"
    assert abs(half - 1.2816 * 0.10) > 1e-2, "must NOT be the old two-sided 80% interval"
