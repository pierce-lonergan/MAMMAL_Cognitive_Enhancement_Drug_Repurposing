"""Pin the Lo-Hi band logic in scripts/133.

The whole point of that script is that a band comparison is only meaningful if the two bands are
scored against the SAME negatives and split on a similarity axis that behaves correctly. Both of
those are easy to break silently, so both are pinned here.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "133_dti_scale_lohi.py"


@pytest.fixture(scope="module")
def s133():
    spec = importlib.util.spec_from_file_location("s133", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_self_similarity_is_one(s133):
    """A molecule compared against itself must score 1.0, or the Hi band is mislabelled."""
    smi = ["CN1C=NC2=C1C(=O)N(C)C(=O)N2C"]  # caffeine
    fps = s133.fingerprints(smi)
    assert s133.max_sim_to(fps, fps) == [1.0]


def test_unparseable_smiles_is_nan_not_zero(s133):
    """An unparseable structure must not silently land in the Lo band as similarity 0."""
    fps = s133.fingerprints(["not_a_molecule"])
    ref = s133.fingerprints(["CN1C=NC2=C1C(=O)N(C)C(=O)N2C"])
    out = s133.max_sim_to(fps, ref)
    assert len(out) == 1 and out[0] != out[0], "expected NaN, got a usable number"


def test_distinct_chemotypes_fall_below_boundary(s133):
    """Caffeine against a steroid should sit well inside the Lo regime."""
    ref = s133.fingerprints(["CN1C=NC2=C1C(=O)N(C)C(=O)N2C"])
    other = s133.fingerprints(["CC12CCC3C(CCC4=CC(=O)CCC34C)C1CCC2O"])  # testosterone
    assert s133.max_sim_to(other, ref)[0] < s133.LOHI_BOUNDARY


def test_both_bands_share_one_negative_pool(s133):
    """Hi and Lo must be scored against identical negatives; only the actives may differ.

    Constructed so the Hi actives outscore the negatives and the Lo actives do not. If the function
    ever started drawing per-band negatives, these two AUROCs would stop being comparable and this
    test would no longer hold its shape.
    """
    rows = []
    for i in range(20):
        rows.append(dict(gene="X", role="active", predicted_pkd=9.0 + i * 0.01, max_sim_ref=0.80))
    for i in range(20):
        rows.append(dict(gene="X", role="active", predicted_pkd=5.0 + i * 0.01, max_sim_ref=0.10))
    for i in range(20):
        rows.append(dict(gene="X", role="negative", predicted_pkd=7.0 + i * 0.01,
                         max_sim_ref=float("nan")))
    s1, s2, _ = s133.analyse_frame(pd.DataFrame(rows))

    bands = {(r.axis, r.band): r for r in s2.itertuples()}
    assert set(a for a, _ in bands) == {"reference"}, "only the reference axis was supplied"
    bands = {b: r for (_, b), r in bands.items()}
    assert bands["Hi"].n == 20 and bands["Lo"].n == 20
    assert bands["Hi"].interpretable and bands["Lo"].interpretable
    assert bands["Hi"].auroc == pytest.approx(1.0)
    assert bands["Lo"].auroc == pytest.approx(0.0)
    # Stage 1 pools both bands, so it must land between them rather than tracking either.
    assert 0.0 < float(s1.iloc[0].auroc) < 1.0


def test_tiny_band_is_flagged_uninterpretable(s133):
    """A band below MIN_BAND must be reported but marked so it cannot be quoted as a result."""
    rows = [dict(gene="X", role="active", predicted_pkd=9.0 + i * 0.01, max_sim_ref=0.80)
            for i in range(20)]
    rows += [dict(gene="X", role="active", predicted_pkd=5.0 + i * 0.01, max_sim_ref=0.10)
             for i in range(4)]
    rows += [dict(gene="X", role="negative", predicted_pkd=7.0 + i * 0.01,
                  max_sim_ref=float("nan")) for i in range(20)]
    _, s2, _ = s133.analyse_frame(pd.DataFrame(rows))
    bands = {r.band: r for r in s2.itertuples()}
    assert bands["Lo"].n == 4
    assert not bands["Lo"].interpretable
    assert bands["Hi"].interpretable


def test_panel_targets_are_unique(s133):
    """A duplicated gene would double-count a target in the stage-1 pass rate."""
    assert len(s133.PANEL) == len(set(s133.PANEL))


def test_both_axes_reported_when_both_columns_present(s133):
    """The pre-registered axis must survive alongside the deviation, not be replaced by it."""
    rows = [dict(gene="X", role="active", predicted_pkd=9.0 + i * 0.01,
                 max_sim_ref=0.80, max_sim_pool=0.10) for i in range(20)]
    rows += [dict(gene="X", role="active", predicted_pkd=5.0 + i * 0.01,
                  max_sim_ref=0.10, max_sim_pool=0.80) for i in range(20)]
    rows += [dict(gene="X", role="negative", predicted_pkd=7.0 + i * 0.01,
                  max_sim_ref=float("nan"), max_sim_pool=float("nan")) for i in range(20)]
    _, s2, _ = s133.analyse_frame(pd.DataFrame(rows))
    assert set(s2.axis) == {"reference", "neighbourhood"}
    got = {(r.axis, r.band): r.auroc for r in s2.itertuples()}
    # The two axes disagree by construction, which is exactly why both are reported.
    assert got[("reference", "Hi")] == pytest.approx(1.0)
    assert got[("neighbourhood", "Hi")] == pytest.approx(0.0)


def test_gradient_recovers_a_monotone_similarity_effect(s133):
    """A planted near-neighbour effect must show up as a rising Q1 to Q4 AUROC."""
    rows = []
    for i in range(80):
        sim = 0.05 + i * 0.01           # 0.05 .. 0.84, uniform across quartiles
        rows.append(dict(gene="X", role="active", predicted_pkd=6.0 + sim * 4,
                         max_sim_ref=sim, max_sim_pool=sim))
    for i in range(40):
        rows.append(dict(gene="X", role="negative", predicted_pkd=7.0 + i * 0.001,
                         max_sim_ref=float("nan"), max_sim_pool=float("nan")))
    g = s133.gradient_frame(pd.DataFrame(rows))
    ref = g[g.axis == "reference"].set_index("quartile")["auroc"]
    assert list(ref.index) == ["Q1", "Q2", "Q3", "Q4"]
    assert ref["Q1"] < ref["Q4"], "planted similarity effect not recovered"
    assert all(ref[f"Q{i}"] <= ref[f"Q{i+1}"] for i in (1, 2, 3))


def test_gradient_is_flat_when_similarity_is_irrelevant(s133):
    """If the score does not depend on similarity, no quartile trend may be manufactured."""
    import numpy as np
    rng = np.random.default_rng(0)
    rows = [dict(gene="X", role="active", predicted_pkd=8.0 + rng.normal(0, 0.1),
                 max_sim_ref=0.05 + i * 0.01, max_sim_pool=0.05 + i * 0.01) for i in range(80)]
    rows += [dict(gene="X", role="negative", predicted_pkd=6.0 + rng.normal(0, 0.1),
                  max_sim_ref=float("nan"), max_sim_pool=float("nan")) for _ in range(40)]
    g = s133.gradient_frame(pd.DataFrame(rows))
    ref = g[g.axis == "reference"].set_index("quartile")["auroc"]
    assert ref.max() - ref.min() < 0.05, "flat input produced a spurious gradient"
