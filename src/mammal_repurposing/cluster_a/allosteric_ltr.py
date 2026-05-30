"""Gap 4 — Allosteric learn-to-rank head (MAMMAL's structural blindness fix).

The honestly-disclosed core weakness: MAMMAL's `dti_bindingdb_pkd` head is a
sequence-only model and is **structurally blind to allosteric / transporter
pharmacology**. Quantified on the cited 21-compound allosteric benchmark, its
within-target predicted-pKd has near-zero variance (std 0.01-0.05) across
ligands spanning three log-units of measured affinity — it literally cannot
rank binders within a target. The Boltz-2 / Boltzina patch (Tier-A criterion)
failed at the transporters.

This module is the real ML contribution: a **learn-to-rank head** that fuses
the heterogeneous binding evidence already in the pipeline —

    [ MAMMAL pKd  ⊕  Tanimoto-to-actives  ⊕  Boltz affinity  ⊕
      RDKit physicochemical descriptors ]

— trained on real ChEMBL pChEMBL affinity labels and evaluated, held-out, on
the allosteric benchmark. The question it answers: *does fusing 3D-structure
(Boltz) + ligand-similarity (Tanimoto) + physicochemistry onto the
sequence-only MAMMAL score recover the within-target affinity ranking that
MAMMAL alone cannot produce — including for the allosteric PAMs/NAMs
(galantamine, zatolmilast)?*

Honest scope: the allosteric benchmark is small (n=21, 5 targets); this is a
proof-of-concept that the fusion direction works + a quantified negative result
on MAMMAL's within-target ranking, not a production affinity predictor. All
metrics use within-target Spearman (the only meaningful axis — any
per-target-constant feature is useless for ranking *within* a target).

sklearn + RDKit + numpy/pandas.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

PHYSCHEM_COLS = [
    "mw", "mollogp", "tpsa", "n_hdonor", "n_hacceptor",
    "n_rotatable", "n_aromatic_rings", "fraction_csp3", "n_heavy",
]
FUSION_FEATURES = ["mammal_pkd", "tanimoto", "boltz_affinity",
                   "boltz_prob", "has_boltz"] + PHYSCHEM_COLS


# ---------------------------------------------------------------------------
# Featurisation
# ---------------------------------------------------------------------------

def physchem_descriptors(smiles: str) -> dict[str, float]:
    """RDKit physicochemical descriptors for one SMILES (NaNs if unparsable)."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
    except Exception:
        return {c: float("nan") for c in PHYSCHEM_COLS}
    m = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if m is None:
        return {c: float("nan") for c in PHYSCHEM_COLS}
    return {
        "mw": float(Descriptors.MolWt(m)),
        "mollogp": float(Descriptors.MolLogP(m)),
        "tpsa": float(rdMolDescriptors.CalcTPSA(m)),
        "n_hdonor": float(Lipinski.NumHDonors(m)),
        "n_hacceptor": float(Lipinski.NumHAcceptors(m)),
        "n_rotatable": float(Descriptors.NumRotatableBonds(m)),
        "n_aromatic_rings": float(rdMolDescriptors.CalcNumAromaticRings(m)),
        "fraction_csp3": float(rdMolDescriptors.CalcFractionCSP3(m)),
        "n_heavy": float(m.GetNumHeavyAtoms()),
    }


def build_feature_table(pairs: pd.DataFrame, *,
                        mammal: pd.DataFrame | None = None,
                        tanimoto: pd.DataFrame | None = None,
                        boltz: pd.DataFrame | None = None) -> pd.DataFrame:
    """Assemble the fusion feature table for a set of (compound, target) pairs.

    `pairs` must have columns: compound_name, target_uniprot, smiles.
    Optional join sources (long, keyed on compound_name + target_uniprot):
      mammal:   predicted_pkd      -> mammal_pkd
      tanimoto: tanimoto_score     -> tanimoto
      boltz:    affinity_pred_value, affinity_probability_binary
    Missing numeric joins are imputed to the per-target mean (then global mean);
    Boltz also carries a has_boltz indicator.
    """
    df = pairs.copy()
    df["target_uniprot"] = df["target_uniprot"].astype(str)
    df["_ck"] = df["compound_name"].str.lower() + "|" + df["target_uniprot"]

    def _join(src, col, newname):
        if src is None or col not in src.columns:
            df[newname] = np.nan
            return
        s = src.copy()
        s["target_uniprot"] = s["target_uniprot"].astype(str)
        s["_ck"] = s["compound_name"].str.lower() + "|" + s["target_uniprot"]
        m = s.drop_duplicates("_ck").set_index("_ck")[col]
        df[newname] = df["_ck"].map(m)

    _join(mammal, "predicted_pkd", "mammal_pkd")
    _join(tanimoto, "tanimoto_score", "tanimoto")
    _join(boltz, "affinity_pred_value", "boltz_affinity")
    _join(boltz, "affinity_probability_binary", "boltz_prob")
    df["has_boltz"] = df["boltz_affinity"].notna().astype(float)

    # physicochemical
    phys = df["smiles"].apply(physchem_descriptors).apply(pd.Series)
    df = pd.concat([df, phys], axis=1)

    # impute numeric features: per-target mean -> global mean -> 0
    for c in [c for c in FUSION_FEATURES if c not in ("has_boltz",)]:
        if c not in df.columns:
            df[c] = np.nan
        df[c] = df.groupby("target_uniprot")[c].transform(
            lambda s: s.fillna(s.mean()))
        df[c] = df[c].fillna(df[c].mean()).fillna(0.0)
    return df.drop(columns=["_ck"])


# ---------------------------------------------------------------------------
# Metrics — within-target Spearman (the only meaningful ranking axis)
# ---------------------------------------------------------------------------

def _spearman(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 3:
        return float("nan")
    def rank(x):
        o = x.argsort(kind="mergesort"); r = np.empty_like(o, dtype=float)
        r[o] = np.arange(len(x)); return r
    ra, rb = rank(a), rank(b)
    ra -= ra.mean(); rb -= rb.mean()
    d = np.sqrt((ra**2).sum() * (rb**2).sum())
    return float((ra * rb).sum() / d) if d > 0 else float("nan")


def within_target_spearman(df: pd.DataFrame, score_col: str, label_col: str,
                           group: str = "target_uniprot",
                           min_n: int = 3) -> tuple[float, dict[str, float]]:
    """Sample-size-weighted mean of per-target Spearman(score, label).
    Returns (pooled_rho, {target: rho})."""
    per: dict[str, float] = {}
    weights, vals = [], []
    for t, g in df.groupby(group):
        if len(g) < min_n:
            continue
        rho = _spearman(g[score_col].to_numpy(float), g[label_col].to_numpy(float))
        if np.isfinite(rho):
            per[str(t)] = rho
            weights.append(len(g)); vals.append(rho)
    pooled = float(np.average(vals, weights=weights)) if vals else float("nan")
    return pooled, per


def mammal_flatness(df: pd.DataFrame, score_col: str = "mammal_pkd",
                    group: str = "target_uniprot") -> dict[str, float]:
    """Within-target std of the MAMMAL score per target — the quantified
    structural-blindness finding (near-zero => cannot rank within target)."""
    return {str(t): float(g[score_col].std())
            for t, g in df.groupby(group) if len(g) >= 2}


# ---------------------------------------------------------------------------
# Learn-to-rank head
# ---------------------------------------------------------------------------

@dataclass
class LTRResult:
    pooled_rho: dict[str, float]                  # condition -> pooled within-target rho
    per_target_rho: dict[str, dict[str, float]]   # condition -> {target: rho}
    flatness: dict[str, float] = field(default_factory=dict)
    n_train: int = 0
    n_eval: int = 0
    allosteric_correct: dict[str, bool] = field(default_factory=dict)
    feature_importance: dict[str, float] = field(default_factory=dict)


def train_fusion_ranker(train: pd.DataFrame, label_col: str,
                        features: list[str] | None = None,
                        *, seed: int = 0):
    """Gradient-boosted regressor mapping fusion features -> affinity (pAct).
    Small + shallow to suit the modest label set."""
    from sklearn.ensemble import GradientBoostingRegressor
    features = features or FUSION_FEATURES
    X = train[features].to_numpy(float)
    y = train[label_col].to_numpy(float)
    model = GradientBoostingRegressor(
        n_estimators=200, max_depth=2, learning_rate=0.05,
        subsample=0.8, random_state=seed)
    model.fit(X, y)
    return model, features


def evaluate(train: pd.DataFrame, benchmark: pd.DataFrame, *,
             label_col: str = "pact", features: list[str] | None = None,
             seed: int = 0) -> LTRResult:
    """Train the fusion head on `train` (ChEMBL pChEMBL) and evaluate within-
    target ranking on the held-out `benchmark`, against single-feature
    baselines (MAMMAL alone, Tanimoto alone, physchem-only model)."""
    features = features or FUSION_FEATURES
    model, feats = train_fusion_ranker(train, label_col, features, seed=seed)
    bench = benchmark.copy()
    bench["fused_score"] = model.predict(bench[feats].to_numpy(float))

    # physchem-only model (no MAMMAL/Tanimoto/Boltz) for ablation
    phys_model, _ = train_fusion_ranker(train, label_col, PHYSCHEM_COLS, seed=seed)
    bench["physchem_score"] = phys_model.predict(bench[PHYSCHEM_COLS].to_numpy(float))

    conditions = {
        "mammal_only": "mammal_pkd",
        "tanimoto_only": "tanimoto",
        "physchem_only": "physchem_score",
        "fused_ltr": "fused_score",
    }
    pooled, per = {}, {}
    for name, col in conditions.items():
        p, pt = within_target_spearman(bench, col, label_col)
        pooled[name] = p
        per[name] = pt

    fi = dict(sorted(zip(feats, model.feature_importances_),
                     key=lambda kv: -kv[1]))
    return LTRResult(
        pooled_rho=pooled, per_target_rho=per,
        flatness=mammal_flatness(bench),
        n_train=int(len(train)), n_eval=int(len(bench)),
        feature_importance=fi,
    )


def availability() -> dict:
    try:
        import sklearn  # noqa
        sk = True
    except Exception:
        sk = False
    try:
        import rdkit  # noqa
        rd = True
    except Exception:
        rd = False
    return {"available": sk and rd, "sklearn": sk, "rdkit": rd,
            "features": FUSION_FEATURES,
            "metric": "within_target_spearman"}
