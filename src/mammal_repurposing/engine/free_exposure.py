"""PERSEUS L1 Stage-3 - efflux-aware, conformal free-brain-penetration regressor.

Replaces the ABSTAIN stub in cns_exposure.py with a quantitative, uncertainty-banded CNS-
penetration call. Design follows the deep-research Gap-2 recommendation:

  * TARGET = logBB (log10 total brain:plasma) from B3DB (CC0). HONEST SCOPE: logBB is a
    passive-penetration proxy, NOT the unbound Kp,uu that finally governs free target exposure
    (Kp,uu is efflux-dominated and its public data are tiny/license-encumbered). So Stage-3 is
    an efflux-AWARE logBB predictor with a calibrated abstain-wide band - a real upgrade over
    the binary heuristic - and true Kp,uu stays the documented residual gap.
  * EFFLUX FEATURE = a cited rule-based P-gp-substrate descriptor (Didziapetris 2003): high
    (N+O) count and MW favour efflux. This is the "model-within-a-model" lever that makes the
    regressor efflux-aware rather than a passive-permeability rehash (ADMET-AI's Pgp head can
    be swapped in where available; the rule keeps it CI-safe and dependency-free).
  * UNCERTAINTY = Mondrian split-conformal (numpy): per-P-gp-category absolute-residual
    quantiles give a distribution-free prediction interval with valid per-subpopulation
    coverage. A kNN applicability-domain check ABSTAINS out-of-manifold queries.
  * The field ceiling is real (public rat Kp,uu R2 ~0.3-0.6); the deliverable is the calibrated
    interval + abstain rule, not a high R2.

  TWO documented upgrade paths (both currently data/dependency-blocked, NOT faked):
    1. EFFLUX: swap the Didziapetris rule for ADMET-AI's Pgp_Broccatelli probability - the
       featurizer threads an external probability through `try_admet_ai_pgp`; `pip install
       admet-ai` then retrain (scripts/111). Falls back to the rule when admet_ai is absent.
    2. TARGET: replace logBB with true unbound rat Kp,uu by adding a Friden-2009 / Morales-2024
       Kp,uu spine (small, partly license-encumbered - request reuse before vendoring). Same
       regressor + conformal pipeline; only the training table changes. Until that data is in
       hand the honest target is logBB and Kp,uu stays the residual gap.

featurize / P-gp / conformal / AD math are RDKit+numpy (CI-safe); training needs lightgbm.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

# fixed, ordered RDKit 2D descriptor set (kept small + robust for a ~1k-row regime)
_DESC = [
    "MolWt", "MolLogP", "TPSA", "NumHDonors", "NumHAcceptors", "NumRotatableBonds",
    "NumAromaticRings", "RingCount", "FractionCSP3", "NumHeteroatoms", "NHOHCount",
    "NOCount", "LabuteASA", "qed",
]
PENETRATION_LOGBB = -1.0   # logBB >= -1 is the established CNS+ cutoff (Clark 2003)


def _descriptors(mol) -> list[float] | None:
    from rdkit.Chem import Crippen, Descriptors, QED, rdMolDescriptors
    try:
        vals = {
            "MolWt": Descriptors.MolWt(mol), "MolLogP": Crippen.MolLogP(mol),
            "TPSA": rdMolDescriptors.CalcTPSA(mol), "NumHDonors": rdMolDescriptors.CalcNumHBD(mol),
            "NumHAcceptors": rdMolDescriptors.CalcNumHBA(mol),
            "NumRotatableBonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
            "NumAromaticRings": rdMolDescriptors.CalcNumAromaticRings(mol),
            "RingCount": rdMolDescriptors.CalcNumRings(mol),
            "FractionCSP3": rdMolDescriptors.CalcFractionCSP3(mol),
            "NumHeteroatoms": rdMolDescriptors.CalcNumHeteroatoms(mol),
            "NHOHCount": Descriptors.NHOHCount(mol), "NOCount": Descriptors.NOCount(mol),
            "LabuteASA": rdMolDescriptors.CalcLabuteASA(mol), "qed": QED.qed(mol),
        }
    except Exception:
        return None
    return [float(vals[k]) for k in _DESC]


def try_admet_ai_pgp(smiles: str) -> float | None:
    """Optional upgrade: ADMET-AI's Pgp_Broccatelli probability as the efflux feature (the
    deep-research Gap-2 recommendation). Returns None if admet_ai is not installed, so the
    rule-based fallback below keeps the pipeline CI-safe and dependency-free. To ADOPT it,
    install admet-ai and retrain with scripts/111 (the featurizer threads it through cleanly)."""
    try:
        from admet_ai import ADMETModel  # noqa: PLC0415
    except Exception:
        return None
    try:
        model = _admet_model()
        pred = model.predict(smiles=smiles)
        for k in ("Pgp_Broccatelli", "Pgp_Inhibition", "Pgp_substrate"):
            if k in pred:
                return float(pred[k])
    except Exception:  # pragma: no cover
        return None
    return None


def _admet_model():
    if not hasattr(_admet_model, "_m"):
        from admet_ai import ADMETModel  # noqa: PLC0415
        _admet_model._m = ADMETModel()
    return _admet_model._m


def pgp_substrate(mol, *, external_prob: float | None = None) -> tuple[float, str]:
    """P-gp-substrate efflux likelihood -> (score in [0,1], category).

    If ``external_prob`` is supplied (e.g. ADMET-AI Pgp_Broccatelli), it is used directly with
    0.3/0.7 cut-points; otherwise the cited rule-based fallback (Didziapetris 2003): efflux is
    favoured by many H-bonding heteroatoms (N+O) and high MW. The rule is the default so the
    fitted Stage-3 model and CI stay dependency-free; the external hook is the documented
    upgrade path."""
    if external_prob is not None:
        if external_prob >= 0.7:
            return float(external_prob), "substrate"
        if external_prob <= 0.3:
            return float(external_prob), "nonsubstrate"
        return float(external_prob), "uncertain"
    from rdkit.Chem import Descriptors, rdMolDescriptors
    n_no = Descriptors.NOCount(mol)            # N + O count
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if n_no >= 8 and mw > 400:
        return 1.0, "substrate"
    if n_no <= 4 and mw < 400:
        return 0.0, "nonsubstrate"
    return 0.5, "uncertain"


def featurize(smiles: str, *, use_admet_ai: bool = False) -> tuple[np.ndarray, str] | None:
    """Return (feature_vector, pgp_category) or None if unparseable. The efflux score is
    appended as the final feature (the model-within-a-model lever). With use_admet_ai=True the
    ADMET-AI Pgp probability replaces the rule when admet_ai is installed (else it falls back
    silently) - keep the train (scripts/111) and inference settings consistent."""
    from rdkit import Chem
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    desc = _descriptors(mol)
    if desc is None:
        return None
    ext = try_admet_ai_pgp(smiles) if use_admet_ai else None
    efflux, cat = pgp_substrate(mol, external_prob=ext)
    return np.asarray(desc + [efflux], dtype=float), cat


FEATURE_NAMES = _DESC + ["pgp_efflux_score"]


# --------------------------------------------------------------------------
# scaffold split + Mondrian split-conformal + applicability domain (numpy)
# --------------------------------------------------------------------------

def bemis_murcko(smiles: str) -> str:
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return ""
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
    except Exception:
        return ""


def scaffold_split(smiles, fracs=(0.7, 0.15, 0.15), seed: int = 0):
    """Group by Bemis-Murcko scaffold, then assign whole scaffolds to train/cal/test so no
    scaffold leaks across splits (blocks analog memorization). Deterministic."""
    groups: dict[str, list[int]] = {}
    for i, s in enumerate(smiles):
        groups.setdefault(bemis_murcko(s) or f"_singleton{i}", []).append(i)
    order = sorted(groups, key=lambda g: (-len(groups[g]), g))   # big scaffolds first, stable
    n = len(smiles)
    cut1, cut2 = fracs[0] * n, (fracs[0] + fracs[1]) * n
    tr, ca, te, seen = [], [], [], 0
    for g in order:
        idx = groups[g]
        bucket = tr if seen < cut1 else (ca if seen < cut2 else te)
        bucket.extend(idx)
        seen += len(idx)
    return sorted(tr), sorted(ca), sorted(te)


def mondrian_quantiles(residuals: np.ndarray, cats: list[str], alpha: float = 0.1) -> dict:
    """Per-category (1-alpha) conformal quantile of absolute residuals, with the finite-sample
    correction ceil((n+1)(1-alpha))/n. Falls back to the pooled quantile for thin categories."""
    res = np.asarray(residuals, dtype=float)
    pooled = _conf_q(res, alpha)
    out = {"_pooled": pooled}
    for c in set(cats):
        r = res[np.asarray(cats) == c]
        out[c] = _conf_q(r, alpha) if len(r) >= 20 else pooled
    return out


def _conf_q(res: np.ndarray, alpha: float) -> float:
    n = len(res)
    if n == 0:
        return float("nan")
    k = math.ceil((n + 1) * (1 - alpha))
    if k > n:
        return float(np.max(res))
    return float(np.sort(res)[k - 1])


@dataclass
class FreeExposurePrediction:
    logbb: float
    lo: float
    hi: float
    pgp_category: str
    pgp_efflux_score: float
    in_domain: bool
    ad_distance: float = float("nan")


@dataclass
class FreeExposureModel:
    model: object                     # fitted LightGBM regressor
    scaler_mean: np.ndarray
    scaler_std: np.ndarray
    train_scaled: np.ndarray          # standardized train descriptors (for kNN AD)
    conformal: dict                   # Mondrian quantiles + "_pooled"
    ad_threshold: float
    ad_k: int = 5
    alpha: float = 0.1
    metrics: dict = field(default_factory=dict)

    def _ad_distance(self, x_scaled: np.ndarray) -> float:
        d = np.sqrt(((self.train_scaled - x_scaled) ** 2).sum(axis=1))
        kth = np.partition(d, min(self.ad_k, len(d) - 1))[:self.ad_k]
        return float(kth.mean())

    def predict(self, smiles: str) -> FreeExposurePrediction | None:
        feat = featurize(smiles)
        if feat is None:
            return None
        x, cat = feat
        x_scaled = (x - self.scaler_mean) / self.scaler_std
        yhat = float(self.model.predict(x.reshape(1, -1))[0])
        q = self.conformal.get(cat, self.conformal["_pooled"])
        adist = self._ad_distance(x_scaled)
        return FreeExposurePrediction(
            logbb=yhat, lo=yhat - q, hi=yhat + q, pgp_category=cat,
            pgp_efflux_score=float(x[-1]), in_domain=adist <= self.ad_threshold, ad_distance=adist)

    def save(self, path) -> None:
        import joblib
        joblib.dump(self, path)

    @staticmethod
    def load(path) -> "FreeExposureModel":
        import joblib
        return joblib.load(path)


# --------------------------------------------------------------------------
# the Stage-3 gate
# --------------------------------------------------------------------------

FAIL, PASS, ABSTAIN = "FAIL", "PASS", "ABSTAIN"


@dataclass
class FreeExposureCall:
    verdict: str
    logbb: float = float("nan")
    lo: float = float("nan")
    hi: float = float("nan")
    pgp_category: str = "?"
    reasons: list[str] = field(default_factory=list)


def free_exposure_gate(pred: FreeExposurePrediction | None,
                       threshold: float = PENETRATION_LOGBB) -> FreeExposureCall:
    """Efflux-aware 3-way Stage-3 verdict on conformal logBB:
      PASS    - lower conformal bound >= threshold (confidently penetrant).
      FAIL    - upper conformal bound < threshold (confidently excluded), OR a predicted P-gp
                substrate whose point logBB is below threshold (efflux likely kills free exposure).
      ABSTAIN - out of the conformal applicability domain, or the band straddles the threshold.
    """
    if pred is None:
        return FreeExposureCall(ABSTAIN, reasons=["unparseable SMILES"])
    c = FreeExposureCall(ABSTAIN, logbb=round(pred.logbb, 2), lo=round(pred.lo, 2),
                         hi=round(pred.hi, 2), pgp_category=pred.pgp_category)
    if not pred.in_domain:
        c.reasons.append(f"out of conformal applicability domain (kNN dist {pred.ad_distance:.2f})")
        return c
    if pred.pgp_category == "substrate" and pred.logbb < threshold:
        c.verdict = FAIL
        c.reasons.append(f"predicted P-gp substrate with logBB {pred.logbb:.2f} < {threshold} "
                         "- efflux likely prevents free brain exposure")
        return c
    if pred.lo >= threshold:
        c.verdict = PASS
        c.reasons.append(f"conformal logBB lower bound {pred.lo:.2f} >= {threshold} "
                         "(confidently CNS-penetrant)")
        return c
    if pred.hi < threshold:
        c.verdict = FAIL
        c.reasons.append(f"conformal logBB upper bound {pred.hi:.2f} < {threshold} "
                         "(confidently non-penetrant)")
        return c
    c.reasons.append(f"logBB band [{pred.lo:.2f}, {pred.hi:.2f}] straddles {threshold} "
                     "- penetration uncertain")
    return c
