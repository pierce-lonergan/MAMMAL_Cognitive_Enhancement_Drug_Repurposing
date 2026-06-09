"""Distill ADMET-AI Pgp into the RDKit descriptors (dependency-free efflux feature).

INNOVATION (post-D1 scrutiny): the ADMET-AI Pgp Stage-3 model is better than the Didziapetris
rule (logBB R2 0.276 vs 0.214) but needs a ~40-model chemprop ensemble LIVE at inference. Can we
keep most of that gain with ZERO runtime dependency by DISTILLING ADMET-AI's Pgp_Broccatelli into
the 14 RDKit 2D descriptors we already compute? P-gp substrate-ness is largely driven by MW / TPSA
/ H-bonding / N+O count, which are in the descriptor set, so a small student model should capture
much of the teacher's signal.

Leakage-controlled measurement (one Bemis-Murcko scaffold split, shared with Stage-3):
  1. DISTILL FIDELITY: train descriptor->Pgp (student) on the TRAIN split only; report test-split
     R2 / Spearman vs the ADMET-AI teacher labels (logbb_b3db_pgp.csv).
  2. STAGE-3 with distilled efflux: generate the student's Pgp for train/cal/test (student never
     saw cal/test), use it as the Stage-3 efflux feature, and measure logBB R2 + conformal
     coverage vs the rule (0.214) and the live-ADMET (0.276) baselines.

Verdict is mechanical from the numbers. If distilled Stage-3 lands at/above the rule, it is a free
win (better efflux feature, no admet_ai at inference) and should become the default. CPU.
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from mammal_repurposing.engine.free_exposure import (
    _descriptors, featurize, mondrian_quantiles, scaffold_split,
)

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("distill_pgp")

ROOT = Path(__file__).resolve().parents[1]
B3DB = ROOT / "data" / "raw" / "logbb_b3db.csv"
PGP = ROOT / "data" / "raw" / "logbb_b3db_pgp.csv"
ALPHA = 0.1
_CATS = ("substrate", "uncertain", "nonsubstrate")
RULE_R2, ADMET_R2 = 0.214, 0.276    # D1 baselines (scripts/111), same split


def _mol(s):
    from rdkit import Chem
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
    return Chem.MolFromSmiles(str(s))


def main() -> int:
    import lightgbm as lgb
    from scipy.stats import spearmanr

    df = pd.read_csv(B3DB)
    pgp = pd.read_csv(PGP)
    pgp_map = dict(zip(pgp["smiles"].astype(str), pgp["pgp_prob"].astype(float)))

    # descriptors + teacher Pgp label, aligned, parseable only
    desc, ys_logbb, ys_pgp, smis = [], [], [], []
    for _, r in df.iterrows():
        s = str(r["smiles"])
        m = _mol(s)
        if m is None or s not in pgp_map:
            continue
        d = _descriptors(m)
        if d is None:
            continue
        desc.append(d); ys_logbb.append(float(r["logbb"])); ys_pgp.append(pgp_map[s]); smis.append(s)
    D = np.vstack(desc); y_pgp = np.asarray(ys_pgp); y_logbb = np.asarray(ys_logbb)
    L.info("Aligned %d compounds (descriptors + teacher Pgp + logBB)", len(D))

    tr, ca, te = scaffold_split(smis, fracs=(0.7, 0.15, 0.15))
    L.info("Scaffold split: train %d / cal %d / test %d", len(tr), len(ca), len(te))

    # 1. DISTILL FIDELITY: descriptor -> ADMET-Pgp, student trained on TRAIN only
    student = lgb.LGBMRegressor(n_estimators=500, learning_rate=0.05, num_leaves=31,
                                min_child_samples=20, subsample=0.8, colsample_bytree=0.8,
                                random_state=0, verbose=-1)
    student.fit(D[tr], y_pgp[tr])
    pgp_hat_te = np.clip(student.predict(D[te]), 0.0, 1.0)
    ss_res = float(((y_pgp[te] - pgp_hat_te) ** 2).sum())
    ss_tot = float(((y_pgp[te] - y_pgp[te].mean()) ** 2).sum())
    fid_r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    fid_rho = float(spearmanr(y_pgp[te], pgp_hat_te).statistic)
    fid_mae = float(np.abs(y_pgp[te] - pgp_hat_te).mean())
    L.info("[distill fidelity] desc->Pgp test R2 %.3f | Spearman %.3f | MAE %.3f",
           fid_r2, fid_rho, fid_mae)

    # student Pgp for EVERY compound (student never saw cal/test -> no leakage into Stage-3 eval)
    pgp_hat_all = np.clip(student.predict(D), 0.0, 1.0)

    # 2. STAGE-3 logBB with the DISTILLED efflux feature
    feats, cats = [], []
    for i, s in enumerate(smis):
        f = featurize(s, pgp_override=float(pgp_hat_all[i]))
        feats.append(f[0]); cats.append(f[1])
    X = np.vstack(feats); cats = np.asarray(cats)

    model = lgb.LGBMRegressor(n_estimators=400, learning_rate=0.05, num_leaves=31,
                              min_child_samples=20, subsample=0.8, colsample_bytree=0.8,
                              random_state=0, verbose=-1)
    model.fit(X[tr], y_logbb[tr])
    yhat_te = model.predict(X[te])
    ss_res = float(((y_logbb[te] - yhat_te) ** 2).sum())
    ss_tot = float(((y_logbb[te] - y_logbb[te].mean()) ** 2).sum())
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    rmse = float(np.sqrt(((y_logbb[te] - yhat_te) ** 2).mean()))

    res_ca = np.abs(y_logbb[ca] - model.predict(X[ca]))
    conformal = mondrian_quantiles(res_ca, list(cats[ca]), alpha=ALPHA)
    cov = 0
    for i in te:
        q = conformal.get(cats[i], conformal["_pooled"])
        yh = float(model.predict(X[i].reshape(1, -1))[0])
        if yh - q <= y_logbb[i] <= yh + q:
            cov += 1
    coverage = cov / len(te)

    L.info("[stage3 distilled] logBB R2 %.3f RMSE %.2f | coverage %.2f", r2, rmse, coverage)
    L.info("[baselines]        rule R2 %.3f | live-ADMET R2 %.3f", RULE_R2, ADMET_R2)

    # mechanical verdict
    gain_vs_rule = r2 - RULE_R2
    recovered = (r2 - RULE_R2) / (ADMET_R2 - RULE_R2) if ADMET_R2 != RULE_R2 else float("nan")
    if gain_vs_rule >= 0.02:
        verdict = (f"DISTILLED efflux BEATS the rule by {gain_vs_rule:+.3f} R2 and recovers "
                   f"{recovered:.0%} of the live-ADMET gain with ZERO inference dependency -> "
                   "WIRE as the dependency-free default efflux feature.")
    elif gain_vs_rule >= -0.01:
        verdict = (f"DISTILLED efflux ~matches the rule ({gain_vs_rule:+.3f} R2). No free win; the "
                   "3-level rule already captures the descriptor-accessible efflux signal. Keep "
                   "the rule default; live-ADMET stays the (heavier) ceiling.")
    else:
        verdict = (f"DISTILLED efflux is WORSE than the rule ({gain_vs_rule:+.3f} R2) - do not "
                   "adopt; the student adds noise. Keep the rule default.")
    L.info("[VERDICT] %s", verdict)
    print("\nVERDICT:", verdict)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
