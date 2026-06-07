"""F1 — the compound-level resolution test.

The headline class-prognostic predictor (``retrospective.class_loco_g``) assigns
every member of a mechanism class the same predicted clinical g (the class
mean). The deepest open question is whether that is the *resolution limit* of
in-silico cognition-drug prognosis, or whether a compound-LEVEL feature can rank
drugs WITHIN a class. This module answers it, pre-registered and leakage-safe.

Two complementary read-outs:

1. **Variance decomposition** (the ceiling). One-way decomposition of clinical_g
   by mechanism_class: how much of the total g variance is BETWEEN classes vs
   WITHIN classes. If the within-class variance is ~0, no compound feature can
   help by construction -- "class is the resolution limit" is then a property of
   the data, not of any particular feature.

2. **Within-class association** (the test). For each candidate compound feature,
   the pooled within-class Spearman correlation with clinical_g (class removed),
   a within-class permutation p-value, a class-cluster bootstrap CI, and a
   leave-one-compound-out test of whether augmenting the class mean with the
   feature lowers held-out g error.

The statistical core (this file) is numpy/scipy only, so it runs in CI and is
fully tested on synthetic data. Molecular-feature construction (RDKit
descriptors, structural typicality) is optional and degrades gracefully.

Honest expectation: on the n=31 ledger the failure classes carry ~0 within-class
g variance (all g≈0, the outcome-pure finding) and the success classes are small
(n<=5), so the test is low-powered. Distinguishing "class is the true ceiling"
from "we lack power" is exactly what the F3 ledger expansion would resolve; this
module reports the power explicitly so the negative is interpretable.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# 1. Variance decomposition (the ceiling on compound-level resolution)
# ---------------------------------------------------------------------------

@dataclass
class VarianceDecomposition:
    n: int
    n_classes: int
    total_var: float
    between_var: float           # SS_between / N  (population convention)
    within_var: float            # SS_within  / N
    frac_between: float          # eta^2 = SS_between / SS_total
    frac_within: float           # 1 - eta^2  (the ceiling for compound features)
    icc1: float                  # one-way random-effects ICC(1)
    per_class: dict[str, dict]   # class -> {n, mean_g, within_sd}


def variance_decomposition(
    ledger: pd.DataFrame,
    value: str = "clinical_g",
    group: str = "mechanism_class",
) -> VarianceDecomposition:
    """Decompose ``value`` variance into between-class and within-class parts.

    ``frac_within`` is the fraction of total variance that lives inside classes:
    it is the hard upper bound on the variance any compound-level feature could
    ever explain. ``icc1`` is the standard one-way random-effects intraclass
    correlation (how much a drug's g is determined by its class identity).
    """
    x = ledger[value].to_numpy(dtype=float)
    g = ledger[group].to_numpy()
    n = len(x)
    if n == 0:
        raise ValueError("empty ledger")
    grand = float(x.mean())
    ss_total = float(((x - grand) ** 2).sum())

    classes = pd.unique(g)
    k = len(classes)
    ss_between = 0.0
    ss_within = 0.0
    per_class: dict[str, dict] = {}
    group_sizes = []
    for c in classes:
        xc = x[g == c]
        nc = len(xc)
        mc = float(xc.mean())
        ss_between += nc * (mc - grand) ** 2
        ss_within += float(((xc - mc) ** 2).sum())
        group_sizes.append(nc)
        per_class[str(c)] = {
            "n": int(nc),
            "mean_g": mc,
            "within_sd": float(xc.std(ddof=1)) if nc > 1 else 0.0,
        }

    frac_between = ss_between / ss_total if ss_total > 0 else 0.0
    # One-way random-effects ICC(1): (MSB - MSW) / (MSB + (n0-1) MSW)
    if k > 1 and n - k > 0:
        msb = ss_between / (k - 1)
        msw = ss_within / (n - k)
        # n0: size-corrected mean group size
        n0 = (n - (sum(s ** 2 for s in group_sizes) / n)) / (k - 1)
        denom = msb + (n0 - 1) * msw
        icc1 = float((msb - msw) / denom) if denom > 0 else 0.0
    else:
        msw = 0.0
        icc1 = float("nan")
    icc1 = max(-1.0, min(1.0, icc1)) if np.isfinite(icc1) else float("nan")

    return VarianceDecomposition(
        n=n, n_classes=k,
        total_var=ss_total / n,
        between_var=ss_between / n,
        within_var=ss_within / n,
        frac_between=float(frac_between),
        frac_within=float(1.0 - frac_between),
        icc1=icc1,
        per_class=per_class,
    )


# ---------------------------------------------------------------------------
# 2. Pooled within-class association (the test)
# ---------------------------------------------------------------------------

def _avg_ranks(a: np.ndarray) -> np.ndarray:
    """1-based average (mid) ranks, matching scipy.stats.rankdata('average')."""
    a = np.asarray(a, dtype=float)
    order = np.argsort(a, kind="mergesort")
    ranks = np.empty(len(a), dtype=float)
    ranks[order] = np.arange(1, len(a) + 1, dtype=float)
    # average ties
    sa = a[order]
    i = 0
    while i < len(sa):
        j = i
        while j + 1 < len(sa) and sa[j + 1] == sa[i]:
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = (i + 1 + j + 1) / 2.0
        i = j + 1
    return ranks


def _pooled_within_rho(feat: np.ndarray, val: np.ndarray, grp: np.ndarray,
                       min_class: int = 2) -> tuple[float, int, int]:
    """Within-class partial Spearman: rank feature and value WITHIN each class,
    center each within its class, then Pearson-correlate the pooled centered
    ranks. Removes the between-class signal, isolating compound-level ordering.

    Returns (rho, n_pooled_points, n_classes_used). Classes that are too small
    or constant in either feature or value contribute nothing.
    """
    fc_all: list[np.ndarray] = []
    vc_all: list[np.ndarray] = []
    n_classes = 0
    for c in pd.unique(grp):
        m = grp == c
        if m.sum() < min_class:
            continue
        f = feat[m].astype(float)
        v = val[m].astype(float)
        ok = np.isfinite(f) & np.isfinite(v)
        if ok.sum() < min_class:
            continue
        f, v = f[ok], v[ok]
        if np.ptp(f) == 0 or np.ptp(v) == 0:   # constant -> no within-class info
            continue
        rf = _avg_ranks(f)
        rv = _avg_ranks(v)
        fc_all.append(rf - rf.mean())
        vc_all.append(rv - rv.mean())
        n_classes += 1
    if n_classes == 0:
        return float("nan"), 0, 0
    fcat = np.concatenate(fc_all)
    vcat = np.concatenate(vc_all)
    denom = np.sqrt((fcat ** 2).sum() * (vcat ** 2).sum())
    if denom == 0:
        return float("nan"), len(fcat), n_classes
    return float((fcat * vcat).sum() / denom), len(fcat), n_classes


@dataclass
class WithinClassResult:
    feature: str
    rho: float                    # pooled within-class partial Spearman
    n_points: int
    n_classes: int
    perm_p: float                 # within-class permutation p (two-sided)
    ci_lo: float                  # class-cluster bootstrap 90% CI
    ci_hi: float
    per_class: dict[str, dict] = field(default_factory=dict)


def within_class_spearman(
    df: pd.DataFrame,
    feature: str,
    value: str = "clinical_g",
    group: str = "mechanism_class",
    min_class: int = 2,
    n_perm: int = 5000,
    n_boot: int = 2000,
    seed: int = 0,
) -> WithinClassResult:
    """Pooled within-class Spearman of ``feature`` vs ``value``, with a
    within-class permutation p-value and a class-cluster bootstrap CI.

    - Permutation null: shuffle ``value`` WITHIN each class (preserving class
      means and sizes), recompute the pooled rho. p = (1+#|null|>=|obs|)/(1+B).
    - Bootstrap CI: resample whole CLASSES with replacement (the unit of
      independence is the class, matching the headline class-cluster CI), 90%.
    """
    rng = np.random.default_rng(seed)
    sub = df[[feature, value, group]].copy()
    feat = sub[feature].to_numpy(dtype=float)
    val = sub[value].to_numpy(dtype=float)
    grp = sub[group].to_numpy()

    rho, n_pts, n_cls = _pooled_within_rho(feat, val, grp, min_class)

    # per-class detail (Spearman within each eligible class)
    per_class: dict[str, dict] = {}
    for c in pd.unique(grp):
        m = grp == c
        if m.sum() < min_class:
            continue
        f, v = feat[m], val[m]
        ok = np.isfinite(f) & np.isfinite(v)
        if ok.sum() < min_class or np.ptp(f[ok]) == 0 or np.ptp(v[ok]) == 0:
            per_class[str(c)] = {"n": int(m.sum()), "rho": float("nan")}
            continue
        r, _, _ = _pooled_within_rho(f[ok], v[ok], np.zeros(ok.sum()), min_class)
        per_class[str(c)] = {"n": int(ok.sum()), "rho": r}

    # within-class permutation p
    if np.isfinite(rho) and n_cls > 0:
        ge = 0
        for _ in range(n_perm):
            vperm = val.copy()
            for c in pd.unique(grp):
                idx = np.where(grp == c)[0]
                vperm[idx] = rng.permutation(vperm[idx])
            r, _, _ = _pooled_within_rho(feat, vperm, grp, min_class)
            if np.isfinite(r) and abs(r) >= abs(rho):
                ge += 1
        perm_p = (1 + ge) / (1 + n_perm)
    else:
        perm_p = float("nan")

    # class-cluster bootstrap CI
    classes = pd.unique(grp)
    if np.isfinite(rho) and len(classes) > 1:
        boots = []
        for _ in range(n_boot):
            pick = rng.choice(classes, size=len(classes), replace=True)
            f_parts, v_parts, g_parts = [], [], []
            for j, c in enumerate(pick):
                m = grp == c
                f_parts.append(feat[m]); v_parts.append(val[m])
                # relabel each drawn class uniquely so duplicates stay separate
                g_parts.append(np.full(m.sum(), f"{c}__{j}"))
            r, _, _ = _pooled_within_rho(
                np.concatenate(f_parts), np.concatenate(v_parts),
                np.concatenate(g_parts), min_class)
            if np.isfinite(r):
                boots.append(r)
        if boots:
            ci_lo = float(np.percentile(boots, 5))
            ci_hi = float(np.percentile(boots, 95))
        else:
            ci_lo = ci_hi = float("nan")
    else:
        ci_lo = ci_hi = float("nan")

    return WithinClassResult(
        feature=feature, rho=rho, n_points=n_pts, n_classes=n_cls,
        perm_p=perm_p, ci_lo=ci_lo, ci_hi=ci_hi, per_class=per_class,
    )


# ---------------------------------------------------------------------------
# 3. Leave-one-compound-out: does the feature beat the class mean?
# ---------------------------------------------------------------------------

@dataclass
class LocoGain:
    feature: str
    mae_classmean: float          # baseline: class mean (excl self), strict LOO
    mae_augmented: float          # class mean + within-class linear feature adj
    delta_mae: float              # baseline - augmented (positive = improvement)
    n: int
    n_adjusted: int               # how many drugs actually received an adjustment


def loco_within_class_mae(
    df: pd.DataFrame,
    feature: str,
    value: str = "clinical_g",
    group: str = "mechanism_class",
    min_train: int = 3,
) -> LocoGain:
    """Leave-one-compound-out MAE: class-mean baseline vs class-mean + a
    within-class OLS slope on the (centered) feature.

    For held-out drug d in class c, the baseline is the mean g of c's *other*
    members (strict LOO; singleton -> global mean excl self). The augmented
    predictor adds a slope*centered-feature term fit on c's training members
    (only when >= ``min_train`` non-constant training points exist; otherwise it
    falls back to the baseline). If the feature carries no within-class signal,
    ``delta_mae`` <= 0 and the class mean is unbeaten.
    """
    g_all = df[value].to_numpy(dtype=float)
    f_all = df[feature].to_numpy(dtype=float)
    grp = df[group].to_numpy()
    n = len(df)
    err_base = np.zeros(n)
    err_aug = np.zeros(n)
    n_adj = 0
    idx_all = np.arange(n)
    for i in range(n):
        mask = idx_all != i
        c = grp[i]
        in_c = mask & (grp == c)
        if in_c.sum() == 0:
            base = float(g_all[mask].mean())   # singleton -> global mean excl self
        else:
            base = float(g_all[in_c].mean())
        err_base[i] = abs(g_all[i] - base)

        aug = base
        if in_c.sum() >= min_train:
            ft = f_all[in_c]; gt = g_all[in_c]
            ok = np.isfinite(ft) & np.isfinite(gt)
            if ok.sum() >= min_train and np.ptp(ft[ok]) > 0:
                ft, gt = ft[ok], gt[ok]
                fmean = ft.mean()
                fc = ft - fmean
                slope = float((fc * (gt - gt.mean())).sum() / (fc ** 2).sum())
                if np.isfinite(f_all[i]):
                    aug = base + slope * (f_all[i] - fmean)
                    n_adj += 1
        err_aug[i] = abs(g_all[i] - aug)

    mae_base = float(err_base.mean())
    mae_aug = float(err_aug.mean())
    return LocoGain(
        feature=feature, mae_classmean=mae_base, mae_augmented=mae_aug,
        delta_mae=mae_base - mae_aug, n=n, n_adjusted=n_adj,
    )


# ---------------------------------------------------------------------------
# 4. Optional molecular features (RDKit-gated; degrade gracefully)
# ---------------------------------------------------------------------------

def rdkit_descriptors(smiles: str) -> dict[str, float]:
    """Physicochemical / CNS-drug-likeness descriptors for one SMILES.

    A CNS-permeable, drug-like molecule is more likely to reach brain exposure
    sufficient for target engagement -- a (loose) proxy for the dose-adequacy
    axis. Returns an empty dict if RDKit is unavailable or the SMILES is
    unparseable. Requires RDKit (not installed in CI).
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Crippen, Descriptors, QED, rdMolDescriptors
    except ImportError:
        return {}
    if not smiles:
        return {}
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {}
    mw = float(Descriptors.MolWt(mol))
    logp = float(Crippen.MolLogP(mol))
    tpsa = float(rdMolDescriptors.CalcTPSA(mol))
    hbd = int(rdMolDescriptors.CalcNumHBD(mol))
    hba = int(rdMolDescriptors.CalcNumHBA(mol))
    # A simple CNS-MPO-style desirability: reward MW<=360, logP in [1,3],
    # TPSA in [40,90], HBD<=1 (higher = more CNS-drug-like). Bounded [0, 5].
    def _hump(v, lo, hi):
        return 1.0 if lo <= v <= hi else max(0.0, 1.0 - min(abs(v - lo), abs(v - hi)) / (hi - lo + 1e-9))
    cns_mpo = (
        _hump(mw, 0, 360) + _hump(logp, 1.0, 3.0) + _hump(tpsa, 40.0, 90.0)
        + (1.0 if hbd <= 1 else max(0.0, 1.0 - (hbd - 1) / 3.0))
        + _hump(float(rdMolDescriptors.CalcNumRotatableBonds(mol)), 0, 7)
    )
    return {
        "mw": mw, "logp": logp, "tpsa": tpsa, "hbd": float(hbd), "hba": float(hba),
        "rotb": float(rdMolDescriptors.CalcNumRotatableBonds(mol)),
        "arom_rings": float(rdMolDescriptors.CalcNumAromaticRings(mol)),
        "fcsp3": float(rdMolDescriptors.CalcFractionCSP3(mol)),
        "qed": float(QED.qed(mol)),
        "cns_mpo": float(cns_mpo),
    }


def class_centroid_tanimoto(
    smiles_by_compound: dict[str, str],
    class_by_compound: dict[str, str],
) -> dict[str, float]:
    """Mean ECFP4 Tanimoto of each compound to its OTHER class members
    (structural typicality within its mechanism class). Empty dict if RDKit is
    unavailable. NaN for a compound with no same-class peer or no fingerprint.
    """
    try:
        from rdkit import Chem, DataStructs
        from rdkit.Chem import rdMolDescriptors
    except ImportError:
        return {}
    fps: dict[str, object] = {}
    for c, smi in smiles_by_compound.items():
        if not smi:
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            fps[c] = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    out: dict[str, float] = {}
    for c in smiles_by_compound:
        cls = class_by_compound.get(c)
        peers = [d for d in fps
                 if d != c and class_by_compound.get(d) == cls]
        if c not in fps or not peers:
            out[c] = float("nan")
            continue
        sims = [DataStructs.TanimotoSimilarity(fps[c], fps[d]) for d in peers]
        out[c] = float(np.mean(sims))
    return out
