"""Gap 3 — Leakage-audited retrospective clinical-outcome validation.

The question: *does the pipeline's evidence anticipate the real pivotal-trial
outcomes of cognition drugs it was never told the outcome of?* — including
the famous Phase III failures (encenicline, idalopirdine, intepirdine,
bitopertin, pomaglumetad …).

This module defines three predictors of clinical SUCCESS vs FAILURE, ordered
by how much information each is allowed to use, and a leakage audit for each:

  P1  TARGET-FIRST (σ(θ̄) at the drug's known target, from V6.B Cluster D)
      — and optionally V6.A within-target binding percentile.
      LEAKAGE: none. V6.B θ̄ is built from GWAS/AHBA/single-cell brain data;
      V6.A binding from ChEMBL bioactivity. Neither saw any cognition trial.
      HYPOTHESIS: target relevance + binding affinity do NOT discriminate
      clinical winners from losers (encenicline binds α7 well, α7 is
      cognition-relevant — yet it failed). Expected AUROC ≈ 0.5.

  P2  CLASS-STRUCTURE LEAVE-ONE-COMPOUND-OUT (the pipeline's class-hierarchical
      prior). For each held-out drug, predict its effect size from the
      meta-analytic track record of its mechanism-class SIBLINGS (empirical-
      Bayes shrinkage toward the global mean), with the drug's OWN outcome
      removed.
      LEAKAGE: the predictor uses siblings' clinical g but NEVER the held-out
      drug's own outcome. This is legitimate inductive generalization
      ("other α7 cognition drugs failed → predict this α7 drug fails").
      HYPOTHESIS: mechanism-class track record is strongly prognostic.

  P3  LEAVE-ONE-CLASS-OUT (the hard extrapolation bound). Predict each drug
      from the mean g of all OTHER mechanism classes — its own class entirely
      removed.
      LEAKAGE: none w.r.t. the class. Tests whether the pipeline can
      extrapolate clinical outcome to an unseen mechanism. Expected weak.

The honest scientific story is the CONTRAST: P1 ≈ chance, P2 ≫ chance, P3 ≈
chance. Picking the right mechanism class (not target affinity) is what
predicts cognition-drug success — exactly what the class-hierarchical
Bayesian pipeline encodes, and consistent with the V6.B Gate-2 falsification.

All metrics (AUROC, bootstrap CI, permutation p, Spearman, failure-recall)
are computed with numpy only — no sklearn dependency.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Ledger
# ---------------------------------------------------------------------------

def load_clinical_ledger(path) -> pd.DataFrame:
    """Load the curated clinical-outcome ledger (binary primary analysis)."""
    df = pd.read_csv(path, comment="#")
    # Restrict to the binary primary analysis
    df = df[df["clinical_outcome"].isin(["SUCCESS", "FAILURE"])].copy()
    df["label"] = (df["clinical_outcome"] == "SUCCESS").astype(int)
    df["compound_lower"] = df["compound"].str.lower()
    return df.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Predictors
# ---------------------------------------------------------------------------

def class_loco_g(ledger: pd.DataFrame, shrinkage_k0: float = 1.0) -> dict[str, float]:
    """P2 — leave-one-COMPOUND-out class-structure predicted g.

    For each compound, predicted g = empirical-Bayes shrinkage of its
    mechanism-class siblings' mean clinical_g toward the global mean, with
    the compound itself excluded:

        ĝ_c = (n_sib · mean_sib + k0 · global_mean) / (n_sib + k0)

    Singleton classes (no siblings) fall back to the global mean.
    """
    global_mean = float(ledger["clinical_g"].mean())
    out: dict[str, float] = {}
    for i, row in ledger.iterrows():
        cls = row["mechanism_class"]
        siblings = ledger[(ledger["mechanism_class"] == cls)
                          & (ledger.index != i)]
        n_sib = len(siblings)
        if n_sib == 0:
            out[row["compound"]] = global_mean
        else:
            sib_mean = float(siblings["clinical_g"].mean())
            out[row["compound"]] = ((n_sib * sib_mean + shrinkage_k0 * global_mean)
                                    / (n_sib + shrinkage_k0))
    return out


def leave_one_class_out_g(ledger: pd.DataFrame) -> dict[str, float]:
    """P3 — leave-one-CLASS-out predicted g (hard extrapolation bound).

    For each compound, predicted g = mean clinical_g of all compounds NOT in
    its mechanism class (its own class entirely removed).
    """
    out: dict[str, float] = {}
    for i, row in ledger.iterrows():
        cls = row["mechanism_class"]
        others = ledger[ledger["mechanism_class"] != cls]
        out[row["compound"]] = (float(others["clinical_g"].mean())
                                if len(others) else float(ledger["clinical_g"].mean()))
    return out


def target_relevance_score(ledger: pd.DataFrame,
                           v6b_theta: pd.DataFrame) -> dict[str, float]:
    """P1a — target cognition-relevance σ(θ̄) at the drug's known target,
    from the V6.B Cluster D posterior. Leakage-free (θ̄ never saw trials)."""
    theta_by_t = {}
    for _, r in v6b_theta.iterrows():
        u = str(r["target_uniprot"])
        w = r.get("w_pipeline", float("nan"))
        if not (isinstance(w, float) and np.isfinite(w)):
            tm = float(r.get("theta_mean", 0.0))
            w = 1.0 / (1.0 + np.exp(-tm))
        theta_by_t[u] = float(w)
    out: dict[str, float] = {}
    for _, row in ledger.iterrows():
        u = str(row["target_uniprot"])
        if u in theta_by_t:
            out[row["compound"]] = theta_by_t[u]
    return out


def target_popularity_score(ledger: pd.DataFrame,
                            chembl_evidence: pd.DataFrame) -> dict[str, float]:
    """Gap-6 baseline — 'well-studied target' / knowledge-hubness proxy: the
    log10 total ChEMBL bioactivity records at the drug's target. Stands in for
    the knowledge-graph repurposing paradigm (target popularity).

    CAVEAT (hindsight confound — NOT leakage-free): a target accrues ChEMBL
    records partly *because* a drug succeeded there (ACHE is saturated with
    cholinesterase-inhibitor records because donepezil worked). The count is a
    contemporaneous snapshot, so high popularity is partly a CONSEQUENCE of the
    clinical outcome, not an a-priori predictor. Reported as an instructive
    confound, not a clean baseline — contrast with the genuinely leakage-free
    affinity / genetics predictors."""
    col = "n_records" if "n_records" in chembl_evidence.columns else None
    if col is None:
        return {}
    pop = (chembl_evidence.assign(_u=chembl_evidence["target_uniprot"].astype(str))
           .groupby("_u")[col].sum().to_dict())
    out: dict[str, float] = {}
    for _, row in ledger.iterrows():
        u = str(row["target_uniprot"])
        if u in pop:
            out[row["compound"]] = float(np.log10(pop[u] + 1.0))
    return out


def binding_score(ledger: pd.DataFrame, v6a_grid: pd.DataFrame) -> dict[str, float]:
    """P1b — V6.A within-target binding percentile at the drug's known target.
    Leakage-free (ChEMBL bioactivity never saw cognition trials). Only defined
    for ledger drugs present in the V6.A grid at their known target."""
    pcol = next((c for c in ("predicted_pkd", "pchembl_mean") if c in v6a_grid.columns), None)
    if pcol is None:
        return {}
    g = v6a_grid.copy()
    g["pct"] = g.groupby("target_uniprot")[pcol].rank(pct=True)
    key = {(str(r["compound_name"]).lower(), str(r["target_uniprot"])): float(r["pct"])
           for _, r in g.iterrows()}
    out: dict[str, float] = {}
    for _, row in ledger.iterrows():
        k = (row["compound_lower"], str(row["target_uniprot"]))
        if k in key:
            out[row["compound"]] = key[k]
    return out


# ---------------------------------------------------------------------------
# Metrics (numpy-only)
# ---------------------------------------------------------------------------

def auroc(scores: np.ndarray, labels: np.ndarray) -> float:
    """AUROC via the Mann–Whitney U relationship. labels ∈ {0,1}."""
    pos = scores[labels == 1]
    neg = scores[labels == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    # rank-based U
    allv = np.concatenate([pos, neg])
    order = allv.argsort(kind="mergesort")
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, len(allv) + 1)
    # average ties
    _, inv, counts = np.unique(allv, return_inverse=True, return_counts=True)
    tie_mean = np.zeros(len(counts))
    cum = 0
    sums = np.zeros(len(counts))
    np.add.at(sums, inv, ranks)
    avg = sums / counts
    ranks = avg[inv]
    r_pos = ranks[:len(pos)].sum()
    u = r_pos - len(pos) * (len(pos) + 1) / 2.0
    return float(u / (len(pos) * len(neg)))


def bootstrap_auroc_ci(scores: np.ndarray, labels: np.ndarray,
                       n_boot: int = 2000, seed: int = 42,
                       alpha: float = 0.10) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    n = len(scores)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        s, l = scores[idx], labels[idx]
        if l.sum() == 0 or l.sum() == len(l):
            continue
        vals.append(auroc(s, l))
    if not vals:
        return (float("nan"), float("nan"))
    lo = float(np.percentile(vals, 100 * alpha / 2))
    hi = float(np.percentile(vals, 100 * (1 - alpha / 2)))
    return (lo, hi)


def permutation_p(scores: np.ndarray, labels: np.ndarray,
                  n_perm: int = 5000, seed: int = 42) -> float:
    """One-sided permutation p-value: P(AUROC_perm >= AUROC_obs) under label
    shuffling."""
    rng = np.random.default_rng(seed)
    obs = auroc(scores, labels)
    if not np.isfinite(obs):
        return float("nan")
    ge = 0
    for _ in range(n_perm):
        perm = rng.permutation(labels)
        if perm.sum() == 0 or perm.sum() == len(perm):
            ge += 1
            continue
        if auroc(scores, perm) >= obs:
            ge += 1
    return (ge + 1) / (n_perm + 1)


def paired_auroc_bootstrap(scores_a: np.ndarray, scores_b: np.ndarray,
                           labels: np.ndarray, n_boot: int = 2000,
                           seed: int = 42) -> dict:
    """Gap-6 — paired bootstrap of AUROC_a − AUROC_b on the SAME resampled rows
    (so the two predictors are compared on identical bootstrap draws). Returns
    the observed delta, a 90% CI, and P(AUROC_a > AUROC_b)."""
    rng = np.random.default_rng(seed)
    n = len(labels)
    obs = auroc(scores_a, labels) - auroc(scores_b, labels)
    deltas, wins = [], 0
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        l = labels[idx]
        if l.sum() == 0 or l.sum() == len(l):
            continue
        da = auroc(scores_a[idx], l) - auroc(scores_b[idx], l)
        if np.isfinite(da):
            deltas.append(da)
            wins += int(da > 0)
    if not deltas:
        return {"delta": obs, "ci_lo": float("nan"), "ci_hi": float("nan"),
                "p_a_gt_b": float("nan")}
    return {
        "delta": float(obs),
        "ci_lo": float(np.percentile(deltas, 5)),
        "ci_hi": float(np.percentile(deltas, 95)),
        "p_a_gt_b": float(wins / len(deltas)),
    }


def spearman(pred: np.ndarray, obs: np.ndarray) -> float:
    if len(pred) < 3:
        return float("nan")
    def rank(a):
        order = a.argsort(kind="mergesort")
        r = np.empty_like(order, dtype=float)
        r[order] = np.arange(len(a))
        return r
    rp, ro = rank(pred), rank(obs)
    rp = (rp - rp.mean()); ro = (ro - ro.mean())
    denom = np.sqrt((rp**2).sum() * (ro**2).sum())
    return float((rp * ro).sum() / denom) if denom > 0 else float("nan")


def failure_recall(scores: dict[str, float], ledger: pd.DataFrame,
                   threshold: float) -> tuple[float, list[str]]:
    """Fraction of true FAILUREs correctly flagged (score < threshold).
    Returns (recall, list of correctly-flagged failure compounds)."""
    fails = ledger[ledger["label"] == 0]
    flagged = []
    n_have = 0
    for _, row in fails.iterrows():
        c = row["compound"]
        if c in scores:
            n_have += 1
            if scores[c] < threshold:
                flagged.append(c)
    recall = (len(flagged) / n_have) if n_have else float("nan")
    return recall, flagged


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

@dataclass
class PredictorResult:
    name: str
    leakage_note: str
    n: int
    auroc: float
    ci_lo: float
    ci_hi: float
    perm_p: float
    spearman_g: float = float("nan")
    extra: dict = field(default_factory=dict)


def evaluate_predictor(name: str, leakage_note: str,
                       scores: dict[str, float], ledger: pd.DataFrame,
                       seed: int = 42) -> PredictorResult:
    rows = ledger[ledger["compound"].isin(scores.keys())]
    s = np.array([scores[c] for c in rows["compound"]])
    y = rows["label"].to_numpy()
    obs_g = rows["clinical_g"].to_numpy()
    au = auroc(s, y)
    lo, hi = bootstrap_auroc_ci(s, y, seed=seed)
    p = permutation_p(s, y, seed=seed)
    rho = spearman(s, obs_g)
    return PredictorResult(
        name=name, leakage_note=leakage_note, n=int(len(rows)),
        auroc=au, ci_lo=lo, ci_hi=hi, perm_p=p, spearman_g=rho,
    )


def availability() -> dict:
    return {
        "available": True,
        "predictors": ["target_relevance(P1a)", "binding(P1b)",
                       "class_loco(P2)", "leave_one_class_out(P3)"],
        "metrics": ["auroc", "bootstrap_ci", "permutation_p", "spearman", "failure_recall"],
        "sklearn_free": True,
    }
