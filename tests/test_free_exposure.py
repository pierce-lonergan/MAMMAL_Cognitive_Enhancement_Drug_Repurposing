"""Tests for L1 Stage-3 - efflux-aware conformal free-exposure (logBB) regressor.

Pure logic (gate, conformal quantile, P-gp rule, scaffold split) is RDKit/numpy and always
runs; the model-dependent path skips if the fitted joblib (scripts/111) or lightgbm is absent.
"""
from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("rdkit")

from mammal_repurposing.engine.free_exposure import (  # noqa: E402
    ABSTAIN, FAIL, PASS, FreeExposurePrediction, _conf_q, featurize, free_exposure_gate,
    mondrian_quantiles, pgp_substrate, scaffold_split,
)

ROOT = Path(__file__).resolve().parents[1]
MODEL = ROOT / "data" / "interim" / "free_exposure_model.joblib"


def test_featurize_and_pgp_rule():
    import numpy as np
    from rdkit import Chem
    feat = featurize("CN(C)CCc1c[nH]c2ccccc12")        # DMT: small, few heteroatoms
    assert feat is not None
    vec, cat = feat
    assert isinstance(vec, np.ndarray) and len(vec) == 15 and cat == "nonsubstrate"
    assert featurize("not_a_smiles") is None
    # a big, H-bond-rich molecule reads as a likely P-gp substrate (Didziapetris)
    big = Chem.MolFromSmiles("OC(=O)C1N2C(=O)C(NC(=O)C(N)c3ccccc3)C2SC1(C)C")  # ampicillin-like
    assert pgp_substrate(big)[1] in ("substrate", "uncertain")
    small = Chem.MolFromSmiles("CC(N)Cc1ccccc1")        # amphetamine
    assert pgp_substrate(small)[1] == "nonsubstrate"


def test_conformal_quantile_finite_sample():
    import numpy as np
    res = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    q = _conf_q(res, alpha=0.1)               # ceil(11*0.9)=10 -> 10th sorted value
    assert q == 1.0
    md = mondrian_quantiles(res, ["a"] * 10, alpha=0.1)
    assert "_pooled" in md and md["_pooled"] == 1.0


def test_scaffold_split_is_disjoint_by_scaffold():
    smis = ["c1ccccc1CCN", "c1ccccc1CCO", "C1CCCCC1CN", "C1CCCCC1CO",
            "c1ccncc1C", "c1ccncc1CC", "CCO", "CCC"]
    tr, ca, te = scaffold_split(smis, fracs=(0.6, 0.2, 0.2))
    assert len(tr) + len(ca) + len(te) == len(smis)
    assert not (set(tr) & set(ca)) and not (set(ca) & set(te)) and not (set(tr) & set(te))


def test_gate_logic_pure():
    # confidently penetrant (lower bound >= -1) -> PASS
    p = FreeExposurePrediction(logbb=0.5, lo=-0.2, hi=1.2, pgp_category="nonsubstrate",
                               pgp_efflux_score=0.0, in_domain=True)
    assert free_exposure_gate(p).verdict == PASS
    # confident P-gp substrate below threshold -> FAIL (efflux)
    e = FreeExposurePrediction(logbb=-1.5, lo=-2.2, hi=-0.8, pgp_category="substrate",
                               pgp_efflux_score=1.0, in_domain=True)
    assert free_exposure_gate(e).verdict == FAIL
    # confidently excluded (upper bound < -1) -> FAIL
    x = FreeExposurePrediction(logbb=-2.0, lo=-2.8, hi=-1.3, pgp_category="nonsubstrate",
                               pgp_efflux_score=0.0, in_domain=True)
    assert free_exposure_gate(x).verdict == FAIL
    # out of applicability domain -> ABSTAIN
    o = FreeExposurePrediction(logbb=0.0, lo=-1.0, hi=1.0, pgp_category="uncertain",
                               pgp_efflux_score=0.5, in_domain=False, ad_distance=9.9)
    assert free_exposure_gate(o).verdict == ABSTAIN
    # band straddles threshold -> ABSTAIN
    s = FreeExposurePrediction(logbb=-0.5, lo=-1.4, hi=0.4, pgp_category="uncertain",
                               pgp_efflux_score=0.5, in_domain=True)
    assert free_exposure_gate(s).verdict == ABSTAIN


@pytest.mark.skipif(not MODEL.exists(), reason="Stage-3 model not trained (run scripts/111)")
def test_trained_model_passes_a_clean_cns_drug():
    pytest.importorskip("lightgbm")
    from mammal_repurposing.engine.free_exposure import FreeExposureModel
    m = FreeExposureModel.load(MODEL)
    pred = m.predict("COc1cc2c(cc1OC)C(=O)C(CC3CCN(Cc4ccccc4)CC3)C2")  # donepezil
    assert pred is not None and pred.in_domain
    assert free_exposure_gate(pred).verdict == PASS


@pytest.mark.skipif(not MODEL.exists(), reason="Stage-3 model not trained")
def test_stage3_preserves_cns_penetrant_psychedelic():
    # psilocin must keep CNS PASS (Stage-3 only downgrades confident P-gp substrates), so the
    # L4 psychoplastogen window is not broken by the new gate
    from mammal_repurposing.engine.cns_exposure import PASS as CPASS, cns_exposure_gate
    assert cns_exposure_gate("CN(C)CCc1c[nH]c2cccc(O)c12").verdict == CPASS