"""C4 convention-lock: analysis/benchmark.py reports a DESCRIPTIVE SD with ddof=0 (population),
the go-forward project convention (BUG_AUDIT_2026-06.md, C4). Also gives this previously
zero-coverage module a smoke test.

C4 is the documented immaterial-by-design item (the value does not change), so this is a forward
guard: it LOCKS ddof=0 so a future flip to ddof=1 fails, and it confirms the input distinguishes the
two conventions.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from mammal_repurposing.analysis.benchmark import analyze_benchmark


def _frame(preds):
    n = len(preds)
    return pd.DataFrame({
        "target_uniprot": ["P1"] * n,
        "target_gene": ["G1"] * n,
        "compound_name": [f"c{i}" for i in range(n)],
        "binding_mode": ["orthosteric"] * n,
        "measured_activity_nm": [np.nan] * n,
        "activity_type": ["Ki"] * n,
        "predicted_pkd": list(preds),
    })


def test_benchmark_pkd_std_uses_population_ddof0():
    preds = [5.0, 6.0, 7.0, 8.0]
    grp = analyze_benchmark(_frame(preds))[0].groups["orthosteric"]
    # population SD (ddof=0) is the convention; a flip to sample SD (ddof=1) would fail the first.
    assert abs(grp.pkd_std - float(np.std(preds, ddof=0))) < 1e-9
    assert abs(grp.pkd_std - float(np.std(preds, ddof=1))) > 1e-6   # input distinguishes the two


def test_benchmark_single_member_group_std_is_zero_not_nan():
    # ddof=0 degrades gracefully on n=1 (ddof=1 would be nan) — the reason for the convention.
    grp = analyze_benchmark(_frame([6.5]))[0].groups["orthosteric"]
    assert grp.pkd_std == 0.0
    assert grp.n == 1
