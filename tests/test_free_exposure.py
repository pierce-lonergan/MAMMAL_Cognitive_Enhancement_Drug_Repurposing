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


def test_pluggable_external_pgp_overrides_rule():
    from rdkit import Chem
    from mammal_repurposing.engine.free_exposure import pgp_substrate, try_admet_ai_pgp
    m = Chem.MolFromSmiles("CC(N)Cc1ccccc1")          # amphetamine -> rule says nonsubstrate
    assert pgp_substrate(m)[1] == "nonsubstrate"
    # an external (e.g. ADMET-AI) probability overrides the rule with 0.3/0.7 cut-points
    assert pgp_substrate(m, external_prob=0.9)[1] == "substrate"
    assert pgp_substrate(m, external_prob=0.05)[1] == "nonsubstrate"
    assert pgp_substrate(m, external_prob=0.5)[1] == "uncertain"
    # the ADMET-AI hook degrades gracefully when admet_ai is not installed
    assert try_admet_ai_pgp("CCO") in (None, ) or isinstance(try_admet_ai_pgp("CCO"), float)


def test_featurize_pgp_override_threads_through():
    # D1: a cached probability (e.g. ADMET-AI Pgp_Broccatelli) passed as pgp_override becomes the
    # last feature AND sets the Mondrian category via the 0.3/0.7 cut-points - without it the
    # Didziapetris rule applies. amphetamine: rule -> nonsubstrate (efflux 0.0).
    rule = featurize("CC(N)Cc1ccccc1")
    ovr = featurize("CC(N)Cc1ccccc1", pgp_override=0.9)
    assert rule is not None and ovr is not None
    assert rule[1] == "nonsubstrate" and rule[0][-1] == 0.0
    assert ovr[1] == "substrate" and abs(ovr[0][-1] - 0.9) < 1e-9


def test_admet_featurize_abstains_when_live_call_unavailable(monkeypatch):
    # D1-scrutiny fix: featurize(use_admet_ai=True) must REFUSE (return None) when the live ADMET
    # call is unavailable, NOT silently feed a rule-derived efflux value into a model trained on
    # ADMET features (that train/inference contract mismatch would be a quiet wrong prediction).
    import mammal_repurposing.engine.free_exposure as fe
    monkeypatch.setattr(fe, "try_admet_ai_pgp", lambda smi: None)
    assert fe.featurize("CCO", use_admet_ai=True) is None       # ADMET requested but down -> abstain
    assert fe.featurize("CCO") is not None                       # default rule path unaffected
    assert fe.featurize("CCO", pgp_override=0.9) is not None     # explicit cached override still works


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


@pytest.mark.skipif(not MODEL.exists(), reason="Stage-3 model not trained")
def test_predict_from_feature_matches_predict():
    # D1: the batch feature-level path must equal the smiles path for the (rule) shipped model,
    # so training-set coverage eval is byte-identical to live inference.
    pytest.importorskip("lightgbm")
    from mammal_repurposing.engine.free_exposure import FreeExposureModel
    m = FreeExposureModel.load(MODEL)
    smi = "COc1cc2c(cc1OC)C(=O)C(CC3CCN(Cc4ccccc4)CC3)C2"   # donepezil
    p1 = m.predict(smi)
    feat = featurize(smi, use_admet_ai=getattr(m, "use_admet_ai", False))
    assert p1 is not None and feat is not None
    p2 = m.predict_from_feature(feat[0], feat[1])
    assert abs(p1.logbb - p2.logbb) < 1e-9 and p1.pgp_category == p2.pgp_category
    assert abs(p1.lo - p2.lo) < 1e-9 and abs(p1.hi - p2.hi) < 1e-9


def test_stage3_default_loads_rule_model_not_admet(monkeypatch):
    # D1: with PERSEUS_STAGE3_ADMET unset the loader must pick the CI-safe rule model
    # (use_admet_ai False) even when the heavier ADMET variant joblib sits alongside it.
    import mammal_repurposing.engine.cns_exposure as cx
    monkeypatch.delenv("PERSEUS_STAGE3_ADMET", raising=False)
    cx._S3_MODEL = "__unset__"                  # reset the lazy singleton
    try:
        m = cx._load_free_exposure_model()
        if m is not None:                       # only meaningful if a model is present
            assert getattr(m, "use_admet_ai", False) is False
    finally:
        cx._S3_MODEL = "__unset__"              # leave the cache clean for other tests


def test_admet_variant_carries_contract_flag():
    # D1: the saved ADMET variant must record use_admet_ai=True so predict() featurizes the same
    # way it was trained (no silent rule/ADMET feature mismatch). Loading the joblib does NOT call
    # admet_ai (only predict() would), so this stays CI-safe.
    admet = ROOT / "data" / "interim" / "free_exposure_model_admet.joblib"
    if not admet.exists():
        pytest.skip("ADMET Stage-3 variant not trained (run scripts/115 + 111)")
    pytest.importorskip("lightgbm")
    from mammal_repurposing.engine.free_exposure import FreeExposureModel
    assert getattr(FreeExposureModel.load(admet), "use_admet_ai", False) is True