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


# ---------------------------------------------------------------------------
# Class-taxonomy sensitivity — is the result an artifact of the grouping?
# ---------------------------------------------------------------------------

# Coarser grouping: collapse the 11 mechanism classes into ~4 neurotransmitter
# systems. This deliberately lumps success and failure classes (AChE-inhibitors
# with α7 agonists under 'cholinergic'), so a drop in AUROC under this map shows
# the signal is specific to the mechanism-class granularity, not any grouping.
COARSE_SYSTEM_MAP = {
    "AChE_inhibitor": "cholinergic", "alpha7_nAChR": "cholinergic",
    "catecholaminergic_ADHD": "monoaminergic", "wake_promoting": "monoaminergic",
    "H3_cognition": "monoaminergic", "multimodal_5HT": "monoaminergic",
    "5HT6_antagonist": "monoaminergic",
    "AMPA_PAM": "glutamatergic", "mGluR": "glutamatergic",
    "NMDA_modulator": "glutamatergic", "PDE9_PDE10": "phosphodiesterase",
}


def auroc_under_taxonomy(ledger: pd.DataFrame, taxonomy, *,
                         shrinkage_k0: float = 1.0) -> float:
    """Class-LOCO AUROC after remapping mechanism_class through `taxonomy`
    (a dict old→new, or a same-length array of new labels)."""
    l2 = ledger.copy().reset_index(drop=True)
    if isinstance(taxonomy, dict):
        l2["mechanism_class"] = l2["mechanism_class"].map(taxonomy)
    else:
        l2["mechanism_class"] = list(taxonomy)
    pred = class_loco_g(l2, shrinkage_k0=shrinkage_k0)
    s = np.array([pred[c] for c in l2["compound"]], float)
    return auroc(s, l2["label"].to_numpy())


def taxonomy_perturbation_test(ledger: pd.DataFrame, *, n_perm: int = 2000,
                               seed: int = 0, shrinkage_k0: float = 1.0) -> dict:
    """Permute the class labels across drugs (preserving class sizes) and
    recompute the class-LOCO AUROC. If the real taxonomy carries genuine signal,
    the observed AUROC sits far above this null. Returns the null mean/SD, 95%
    interval, the fraction of permutations reaching the observed AUROC, and a
    permutation p-value."""
    rng = np.random.default_rng(seed)
    obs = auroc_under_taxonomy(ledger, dict(zip(ledger["mechanism_class"],
                                                 ledger["mechanism_class"])),
                               shrinkage_k0=shrinkage_k0)
    labels = ledger["mechanism_class"].to_numpy()
    vals = []
    for _ in range(n_perm):
        perm = labels.copy(); rng.shuffle(perm)
        vals.append(auroc_under_taxonomy(ledger, perm, shrinkage_k0=shrinkage_k0))
    vals = np.array(vals, float)
    ge = int(np.sum(vals >= obs - 1e-9))
    return {"observed": float(obs), "null_mean": float(vals.mean()),
            "null_sd": float(vals.std()), "null_lo": float(np.percentile(vals, 2.5)),
            "null_hi": float(np.percentile(vals, 97.5)),
            "frac_reaching_observed": float(ge / n_perm),
            "perm_p": float((ge + 1) / (n_perm + 1))}


# ---------------------------------------------------------------------------
# Temporal validation (pseudo-prospective) — the curation-proof test
# ---------------------------------------------------------------------------

def temporal_holdout_auroc(ledger: pd.DataFrame, cutoff_year: int, *,
                           shrinkage_k0: float = 1.0) -> dict:
    """Train the class prior ONLY on drugs that read out ≤ cutoff_year; predict
    the strictly-later drugs. This is a pseudo-prospective test: every prediction
    uses only information available before the test drug's readout.

    A held-out drug is scored by the empirical-Bayes-shrunken mean clinical_g of
    its mechanism-class siblings *in the training window*; drugs whose class is
    unseen before the cutoff fall back to the training-window global mean (an
    honest 'no class history' prediction). Returns the test AUROC plus coverage
    (how many test drugs had a same-class precedent)."""
    train = ledger[ledger["readout_year"] <= cutoff_year]
    test = ledger[ledger["readout_year"] > cutoff_year]
    if len(train) == 0 or test["label"].nunique() < 2:
        # AUROC undefined (test set one-class); still report coverage + recall
        gmean = float(train["clinical_g"].mean()) if len(train) else 0.0
        cls_mean = train.groupby("mechanism_class")["clinical_g"].mean().to_dict()
        seen = {c for c in train["mechanism_class"].unique()}
        scores, labels, covered = [], [], 0
        for _, r in test.iterrows():
            c = r["mechanism_class"]
            covered += int(c in seen)
            n = int((train["mechanism_class"] == c).sum())
            s = ((n * cls_mean[c] + shrinkage_k0 * gmean) / (n + shrinkage_k0)
                 if c in seen else gmean)
            scores.append(s); labels.append(int(r["label"]))
        scores = np.array(scores, float); labels = np.array(labels)
        # failure-recall: of test FAILURES in a class already failing by cutoff,
        # how many are predicted below the train global mean
        fail_mask = labels == 0
        recall = (float(np.mean(scores[fail_mask] < gmean)) if fail_mask.any()
                  else float("nan"))
        return {"cutoff": int(cutoff_year), "auroc": float("nan"),
                "n_train": int(len(train)), "n_test": int(len(test)),
                "test_pos": int(labels.sum()), "test_neg": int((1 - labels).sum()),
                "coverage": covered / max(1, len(test)),
                "failure_recall_vs_trainmean": recall}
    gmean = float(train["clinical_g"].mean())
    cls_mean = train.groupby("mechanism_class")["clinical_g"].mean().to_dict()
    seen = set(train["mechanism_class"].unique())
    scores, labels, covered = [], [], 0
    for _, r in test.iterrows():
        c = r["mechanism_class"]; covered += int(c in seen)
        n = int((train["mechanism_class"] == c).sum())
        s = ((n * cls_mean[c] + shrinkage_k0 * gmean) / (n + shrinkage_k0)
             if c in seen else gmean)
        scores.append(s); labels.append(int(r["label"]))
    scores = np.array(scores, float); labels = np.array(labels)
    return {"cutoff": int(cutoff_year), "auroc": float(auroc(scores, labels)),
            "n_train": int(len(train)), "n_test": int(len(test)),
            "test_pos": int(labels.sum()), "test_neg": int((1 - labels).sum()),
            "coverage": covered / len(test),
            "failure_recall_vs_trainmean": float("nan")}


def prequential_class_loco(ledger: pd.DataFrame, *,
                           shrinkage_k0: float = 1.0) -> dict:
    """As-of (prequential) evaluation: predict EACH drug using only drugs that
    read out STRICTLY BEFORE it. The gold-standard temporal design — no fixed
    cutoff, every drug judged against the knowledge available at its own readout.

    For drug d (readout year y_d), siblings = same-class drugs with
    readout_year < y_d. If ≥1, predicted score = EB-shrunken sibling mean g vs
    the strictly-earlier global mean; otherwise the drug is 'uninformed'
    (no same-class precedent) and excluded from the informed AUROC (its
    fallback = earlier global mean is still recorded for the full-coverage AUROC).

    Returns AUROC over informed predictions, AUROC over all (with fallback),
    coverage, and the per-drug table."""
    rows = []
    for _, r in ledger.iterrows():
        y = r["readout_year"]; c = r["mechanism_class"]
        earlier = ledger[ledger["readout_year"] < y]
        sibs = earlier[earlier["mechanism_class"] == c]
        gmean = float(earlier["clinical_g"].mean()) if len(earlier) else float("nan")
        informed = len(sibs) > 0 and np.isfinite(gmean)
        if informed:
            n = len(sibs)
            score = (n * float(sibs["clinical_g"].mean()) + shrinkage_k0 * gmean) / (n + shrinkage_k0)
        else:
            score = gmean
        rows.append({"compound": r["compound"], "year": int(y),
                     "mechanism_class": c, "label": int(r["label"]),
                     "informed": bool(informed),
                     "n_prior_sibs": int(len(sibs)), "score": score})
    tab = pd.DataFrame(rows)
    inf = tab[tab["informed"] & np.isfinite(tab["score"])]
    full = tab[np.isfinite(tab["score"])]
    au_inf = (float(auroc(inf["score"].to_numpy(), inf["label"].to_numpy()))
              if inf["label"].nunique() == 2 else float("nan"))
    au_full = (float(auroc(full["score"].to_numpy(), full["label"].to_numpy()))
               if full["label"].nunique() == 2 else float("nan"))
    return {"auroc_informed": au_inf, "n_informed": int(len(inf)),
            "informed_pos": int(inf["label"].sum()),
            "informed_neg": int((1 - inf["label"]).sum()),
            "auroc_full_with_fallback": au_full, "n_full": int(len(full)),
            "coverage_informed": len(inf) / len(tab), "table": tab}


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


def kg_network_score(ledger: pd.DataFrame, kg_scores: pd.DataFrame,
                     col: str = "kg_ppr_sum") -> dict[str, float]:
    """Gap-6 baseline — network-propagation comparator (the KG/GNN paradigm).

    Personalised-PageRank mass that diffuses from the drug node to the
    cognition target panel over PrimeKG (`kg_ppr_sum`), i.e. a random-walk /
    network-propagation repurposing score. This is the concrete stand-in for
    the connectivity-/network-propagation methods reviewers expect as a
    comparator.

    CAVEAT (hindsight confound — NOT leakage-free, like target popularity): a
    drug accrues PrimeKG edges — and therefore PageRank mass to the panel —
    partly *because* it was studied and succeeded. High network connectivity is
    thus partly a CONSEQUENCE of the outcome. Reported as an instructive
    confounded baseline, not a clean a-priori predictor."""
    m = {str(r["compound_name"]).lower(): float(r[col])
         for _, r in kg_scores.iterrows()
         if col in kg_scores.columns and np.isfinite(float(r[col]))}
    return {row["compound"]: m[row["compound_lower"]]
            for _, row in ledger.iterrows() if row["compound_lower"] in m}


def structure_nn_success_score(ledger: pd.DataFrame, compounds: pd.DataFrame,
                               *, radius: int = 2, n_bits: int = 2048) -> dict[str, float]:
    """Structure-similarity comparator — leave-one-out nearest-successful-drug.

    For each ledger drug, the maximum Morgan-fingerprint Tanimoto to any OTHER
    ledger drug that historically SUCCEEDED. This is the chemical-structure
    analogue of the class-track-record predictor: both use the historical
    outcomes of the *other* drugs, but one aggregates by chemical similarity
    and the other by mechanism class. The contrast (structure ≪ class) isolates
    *which kind* of historical aggregation carries the signal.

    Requires RDKit; returns {} (gracefully) if unavailable or no SMILES join.
    Only defined for ledger drugs whose SMILES are in the compound library."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, DataStructs
    except Exception:
        return {}
    sm = {str(r["name"]).lower(): r["smiles"] for _, r in compounds.iterrows()}
    fps, lab = {}, {}
    for _, row in ledger.iterrows():
        s = sm.get(row["compound_lower"])
        if not s:
            continue
        mol = Chem.MolFromSmiles(s)
        if mol is None:
            continue
        fps[row["compound"]] = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
        lab[row["compound"]] = int(row["label"])
    out: dict[str, float] = {}
    for c in fps:
        sims = [DataStructs.TanimotoSimilarity(fps[c], fps[o])
                for o in fps if o != c and lab[o] == 1]  # LOO: other SUCCESSES only
        out[c] = float(max(sims)) if sims else 0.0
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


def class_cluster_bootstrap_auroc(ledger: pd.DataFrame, *, n_boot: int = 2000,
                                  seed: int = 42, alpha: float = 0.10,
                                  shrinkage_k0: float = 1.0) -> dict:
    """Class-LEVEL (cluster) bootstrap CI for the class-LOCO AUROC.

    The drug-level bootstrap underestimates uncertainty when the predictor is
    class-aggregated: it cannot vary the *class* composition. Here the
    resampling unit is the **mechanism class** — we sample the classes with
    replacement (renaming replicates so each is an independent cluster with its
    own LOCO siblings), rebuild the ledger, recompute the class-LOCO predictor
    and its AUROC. This propagates the genuine between-class variance (with only
    ~11 classes, the resulting CI is much wider than the drug-level one).

    Returns observed AUROC, class-level 90% CI, and the fraction of resamples
    that were degenerate (all-success or all-failure class draws → AUROC
    undefined, excluded)."""
    rng = np.random.default_rng(seed)
    classes = ledger["mechanism_class"].unique()
    # observed (point estimate) on the real ledger
    pred0 = class_loco_g(ledger, shrinkage_k0=shrinkage_k0)
    s0 = np.array([pred0[c] for c in ledger["compound"]], float)
    obs = auroc(s0, ledger["label"].to_numpy())

    vals, degenerate = [], 0
    for _ in range(n_boot):
        draw = rng.choice(classes, size=len(classes), replace=True)
        parts = []
        for i, c in enumerate(draw):
            sub = ledger[ledger["mechanism_class"] == c].copy()
            sub["mechanism_class"] = f"{c}__rep{i}"        # independent cluster
            sub["compound"] = sub["compound"].astype(str) + f"__rep{i}"
            parts.append(sub)
        rl = pd.concat(parts, ignore_index=True)
        if rl["label"].nunique() < 2:
            degenerate += 1
            continue
        pred = class_loco_g(rl, shrinkage_k0=shrinkage_k0)
        s = np.array([pred[c] for c in rl["compound"]], float)
        a = auroc(s, rl["label"].to_numpy())
        if np.isfinite(a):
            vals.append(a)
    if not vals:
        return {"auroc": obs, "ci_lo": float("nan"), "ci_hi": float("nan"),
                "n_classes": int(len(classes)), "frac_degenerate": float("nan")}
    return {
        "auroc": float(obs),
        "ci_lo": float(np.percentile(vals, 100 * alpha / 2)),
        "ci_hi": float(np.percentile(vals, 100 * (1 - alpha / 2))),
        "median": float(np.median(vals)),
        "n_classes": int(len(classes)),
        "frac_degenerate": float(degenerate / n_boot),
    }


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
