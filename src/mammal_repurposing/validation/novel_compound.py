"""F2 - novel-compound onboarding engine.

The retrospective validator (retrospective.py) and the scaling study (F3) showed
that a compound's mechanism *class* carries the prognostic signal for cognition
(class-LOCO AUROC ~0.92 at n=125), while no within-class compound feature beats
the class mean (F1, ICC 0.95). This module turns that finding into a *prospective*
screen: given an arbitrary novel SMILES, route it to a known cognition mechanism
class and return that class's calibrated clinical-g prior - or ABSTAIN.

The decision path (GAPS F2):

    novel SMILES
      -> structural class assignment:
           (a) max Tanimoto (ECFP4) to each class's ledger exemplars
           (b) Murcko generic-scaffold match to class exemplars
           (c) [pluggable] nearest class in multi-head DTI-profile space
      -> class prior: EB-shrunk class mean clinical_g + 90% bootstrap CrI
      -> confidence tier with two non-negotiable guardrails:
           1. ABSTAIN on genuinely novel mechanisms. A novel compound that is not
              structurally near ANY known cognition class (max Tanimoto < TAU_OOD)
              is out-of-manifold; the engine re-ranks KNOWN mechanisms, it does not
              invent new ones (the leave-one-class-out=0.00 result forbids it).
           2. The V6.A allosteric blindness makes profile/structure class
              assignment unreliable for allosteric chemotypes, so allosteric-flagged
              classes are downgraded (HIGH -> MED) and carry a note.

The statistical/structural core here is numpy/pandas + RDKit; it runs in CI and is
fully testable. The thresholds are not guessed: scripts/95 calibrates them from
leave-one-compound-out class recovery on the exemplar base, then locks the values
below. Nothing in this module fabricates a clinical outcome - the returned g is a
model PREDICTION from the class prior, always labelled as such.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Decision thresholds (calibrated by scripts/95 leave-one-compound-out, then
# locked here). Tanimoto is ECFP4 / 2048-bit, "max-to-class-exemplar".
# ---------------------------------------------------------------------------
# Calibrated by scripts/95 leave-one-compound-out on the exemplar base: every
# observed mis-route sat in Tanimoto [0.26, 0.34], while correct routes were at
# 0.25 or >= 0.40. A 0.35 out-of-manifold floor (also the conventional ECFP4
# weak-similarity cutoff) excludes every observed error -> 100% class recovery on
# the routed set, abstaining on the rest. Precision is prioritised over recall:
# in a prospective screen a confidently wrong route is worse than an abstention.
TAU_OOD = 0.35       # max sim below this -> out-of-manifold -> ABSTAIN
TAU_HIGH = 0.45      # max sim at/above this (+ margin + populated prior) -> HIGH
TAU_MARGIN = 0.05    # (top1 - top2) below this with opposite-sign priors -> ambiguous
MIN_CLASS_N = 2      # a class needs >= this many ledger members to carry a prior

FP_RADIUS = 2
FP_BITS = 2048


# ---------------------------------------------------------------------------
# Allosteric-blind classes (V6.A). Structure/DTI-profile class assignment is
# unreliable for allosteric chemotypes, so these are downgraded, not trusted at
# HIGH. Matched by name pattern so new class names are covered automatically.
# ---------------------------------------------------------------------------
_ALLOSTERIC_PATTERNS = re.compile(
    r"(pam|nam|allosteric|alpha7|nachr|ampa|mglur|m1_m4|m4_pam|muscarinic|sigma)",
    re.IGNORECASE,
)


def is_allosteric_class(mechanism_class: str) -> bool:
    """True if the class is an allosteric/modulator chemotype where DTI-profile
    and structural class assignment are unreliable (the V6.A blindness)."""
    return bool(_ALLOSTERIC_PATTERNS.search(str(mechanism_class)))


# ---------------------------------------------------------------------------
# RDKit helpers (lazy; degrade gracefully if RDKit is unavailable)
# ---------------------------------------------------------------------------

def _morgan_fp(smiles: str):
    """ECFP4 (radius 2, 2048-bit) Morgan fingerprint, or None on parse failure."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except Exception:  # pragma: no cover - environment without rdkit
        return None
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, FP_RADIUS, nBits=FP_BITS)


def _tanimoto(fp_a, fp_b) -> float:
    from rdkit import DataStructs
    return float(DataStructs.TanimotoSimilarity(fp_a, fp_b))


def _murcko_generic(smiles: str) -> str | None:
    """Generic (atom/bond-agnostic) Murcko scaffold SMILES, or None on failure.
    Generic so that, e.g., two donepezil-like benzylpiperidines match on skeleton
    even with different substituents."""
    try:
        from rdkit import Chem
        from rdkit.Chem.Scaffolds import MurckoScaffold
    except Exception:  # pragma: no cover
        return None
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    try:
        scaf = MurckoScaffold.GetScaffoldForMol(mol)
        generic = MurckoScaffold.MakeScaffoldGeneric(scaf)
        return Chem.MolToSmiles(generic) if generic is not None else None
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Class priors (clinical-g prior + 90% CrI per mechanism class)
# ---------------------------------------------------------------------------

@dataclass
class ClassPrior:
    mechanism_class: str
    n: int                  # ledger members in the class
    prior_g: float          # EB-shrunk class mean clinical_g
    g_ci_lo: float          # 5th pctl, class-member bootstrap
    g_ci_hi: float          # 95th pctl
    p_success: float        # EB-shrunk class success rate
    allosteric: bool


def build_class_priors(ledger: pd.DataFrame, *, k0: float = 1.0,
                       n_boot: int = 2000, seed: int = 0) -> dict[str, ClassPrior]:
    """Per mechanism class: EB-shrunk clinical-g prior + 90% class-member bootstrap
    CrI + EB-shrunk success rate. Shrinkage is toward the global mean with pseudo-
    count k0, matching trial_watch / retrospective.

        prior_g = (n * class_mean + k0 * global_mean) / (n + k0)

    The CrI resamples the class's own members (with replacement) and recomputes the
    shrunk mean, so a thin class gets an honestly wide interval."""
    global_g = float(ledger["clinical_g"].mean())
    base_p = float(ledger["label"].mean())
    rng = np.random.default_rng(seed)
    out: dict[str, ClassPrior] = {}
    for cls, g in ledger.groupby("mechanism_class"):
        gv = g["clinical_g"].to_numpy(dtype=float)
        n = len(gv)
        shrunk = (n * gv.mean() + k0 * global_g) / (n + k0)
        if n >= 2:
            boots = np.empty(n_boot, dtype=float)
            for b in range(n_boot):
                samp = rng.choice(gv, size=n, replace=True)
                boots[b] = (n * samp.mean() + k0 * global_g) / (n + k0)
            lo, hi = float(np.percentile(boots, 5)), float(np.percentile(boots, 95))
        else:
            lo = hi = float(shrunk)
        k = int(g["label"].sum())
        p = (k + k0 * base_p) / (n + k0)
        out[str(cls)] = ClassPrior(str(cls), n, float(shrunk), lo, hi, float(p),
                                   is_allosteric_class(cls))
    return out


# ---------------------------------------------------------------------------
# Class exemplar library (fingerprints + scaffolds keyed by class)
# ---------------------------------------------------------------------------

@dataclass
class _Exemplar:
    compound: str
    fp: object
    scaffold: str | None


def build_exemplars(ledger: pd.DataFrame,
                    smiles_df: pd.DataFrame) -> dict[str, list[_Exemplar]]:
    """Group ledger compounds (that have a parseable SMILES) into a per-class
    exemplar library of (compound, fingerprint, generic-scaffold). ``smiles_df``
    has columns [compound, smiles]; matched case-insensitively to the ledger."""
    smi = smiles_df.copy()
    smi["_k"] = smi["compound"].astype(str).str.lower().str.strip()
    smap = dict(zip(smi["_k"], smi["smiles"]))
    lib: dict[str, list[_Exemplar]] = {}
    for _, row in ledger.iterrows():
        k = str(row["compound"]).lower().strip()
        s = smap.get(k)
        if s is None or (isinstance(s, float) and np.isnan(s)):
            continue
        fp = _morgan_fp(s)
        if fp is None:
            continue
        lib.setdefault(str(row["mechanism_class"]), []).append(
            _Exemplar(str(row["compound"]), fp, _murcko_generic(s)))
    return lib


# ---------------------------------------------------------------------------
# Assignment + scoring
# ---------------------------------------------------------------------------

@dataclass
class NovelScore:
    query_id: str
    smiles: str
    assigned_class: str | None
    similarity: float          # max Tanimoto to assigned class
    runner_up_class: str | None
    margin: float              # similarity - runner-up similarity
    scaffold_hit: bool
    n_class_members: int
    prior_g: float
    g_ci_lo: float
    g_ci_hi: float
    p_success: float
    predicted_outcome: str     # SUCCESS | FAILURE | n/a
    allosteric_flag: bool
    tier: str                  # HIGH | MED | LOW | ABSTAIN
    reason: str


def assign_class(query_fp, exemplars: dict[str, list[_Exemplar]],
                 query_scaffold: str | None = None,
                 *, external_class_scores: dict[str, float] | None = None):
    """Rank classes for a query fingerprint by max-Tanimoto to class exemplars.

    Returns (ranked, scaffold_hit) where ranked is a list of (class, sim) sorted
    descending. ``external_class_scores`` (e.g. a DTI-profile-nearest-class score
    in [0,1] per class) is averaged 50/50 with the Tanimoto score when provided -
    the pluggable hook for the GPU multi-head profile signal."""
    sims: dict[str, float] = {}
    scaffold_hit_class: set[str] = set()
    for cls, exs in exemplars.items():
        best = 0.0
        for e in exs:
            t = _tanimoto(query_fp, e.fp)
            if t > best:
                best = t
            if query_scaffold is not None and e.scaffold == query_scaffold:
                scaffold_hit_class.add(cls)
        sims[cls] = best
    if external_class_scores:
        for cls in sims:
            if cls in external_class_scores:
                sims[cls] = 0.5 * sims[cls] + 0.5 * float(external_class_scores[cls])
    ranked = sorted(sims.items(), key=lambda kv: kv[1], reverse=True)
    top_class = ranked[0][0] if ranked else None
    return ranked, (top_class in scaffold_hit_class if top_class else False)


def score_compound(query_id: str, smiles: str,
                   exemplars: dict[str, list[_Exemplar]],
                   priors: dict[str, ClassPrior], *,
                   tau_ood: float = TAU_OOD, tau_high: float = TAU_HIGH,
                   tau_margin: float = TAU_MARGIN, min_class_n: int = MIN_CLASS_N,
                   external_class_scores: dict[str, float] | None = None
                   ) -> NovelScore:
    """Route one novel compound and return its scored onboarding record.

    Tier logic (after the structural assignment):
      ABSTAIN  - SMILES unparseable, no exemplar library, or max sim < tau_ood
                 (out-of-manifold / novel mechanism), or ambiguous top-2 with
                 opposite-sign priors.
      LOW      - assigned, but class prior is thin (n < min_class_n).
      MED      - assigned with a populated prior; moderate similarity, or an
                 allosteric-blind class (downgraded from HIGH).
      HIGH     - clean single-class map: sim >= tau_high, margin >= tau_margin,
                 populated non-allosteric prior.
    """
    na = float("nan")
    fp = _morgan_fp(smiles)
    if fp is None:
        return NovelScore(query_id, smiles, None, na, None, na, False, 0,
                          na, na, na, na, "n/a", False, "ABSTAIN",
                          "unparseable SMILES")
    if not exemplars:
        return NovelScore(query_id, smiles, None, na, None, na, False, 0,
                          na, na, na, na, "n/a", False, "ABSTAIN",
                          "no exemplar library")
    scaf = _murcko_generic(smiles)
    ranked, scaffold_hit = assign_class(fp, exemplars, scaf,
                                        external_class_scores=external_class_scores)
    top_class, top_sim = ranked[0]
    run_class, run_sim = (ranked[1] if len(ranked) > 1 else (None, 0.0))
    margin = top_sim - run_sim

    # guardrail 1: out-of-manifold -> abstain (do not invent a mechanism)
    if top_sim < tau_ood:
        return NovelScore(query_id, smiles, None, top_sim, run_class, margin,
                          scaffold_hit, 0, na, na, na, na, "n/a", False,
                          "ABSTAIN",
                          f"out-of-manifold: max Tanimoto {top_sim:.2f} < {tau_ood:.2f} "
                          f"(nearest known class '{top_class}')")

    pr = priors.get(top_class)
    run_pr = priors.get(run_class) if run_class else None
    # ambiguity: near-tie whose two classes predict OPPOSITE outcomes
    if (run_pr is not None and pr is not None and margin < tau_margin
            and (pr.p_success >= 0.5) != (run_pr.p_success >= 0.5)):
        return NovelScore(query_id, smiles, None, top_sim, run_class, margin,
                          scaffold_hit, pr.n, na, na, na, na, "n/a", False,
                          "ABSTAIN",
                          f"ambiguous: '{top_class}' vs '{run_class}' tie "
                          f"(margin {margin:.2f}) with opposite-sign priors")

    if pr is None:
        return NovelScore(query_id, smiles, top_class, top_sim, run_class, margin,
                          scaffold_hit, 0, na, na, na, na, "n/a", False,
                          "ABSTAIN", f"assigned '{top_class}' has no class prior")

    allo = pr.allosteric
    outcome = "SUCCESS" if pr.p_success >= 0.5 else "FAILURE"

    if pr.n < min_class_n:
        tier, reason = "LOW", (f"thin prior: class '{top_class}' has n={pr.n} "
                               f"(< {min_class_n}); single-precedent routing")
    elif top_sim >= tau_high and margin >= tau_margin and not allo:
        tier, reason = "HIGH", (f"clean map to '{top_class}' (sim {top_sim:.2f}, "
                                f"margin {margin:.2f}, n={pr.n})")
    elif allo:
        tier, reason = "MED", (f"assigned '{top_class}' (sim {top_sim:.2f}) but "
                               f"allosteric-blind class: structural routing "
                               f"downgraded (V6.A)")
    else:
        tier, reason = "MED", (f"moderate map to '{top_class}' (sim {top_sim:.2f}, "
                               f"margin {margin:.2f}, n={pr.n})")

    return NovelScore(query_id, smiles, top_class, top_sim, run_class, margin,
                      bool(scaffold_hit), pr.n, pr.prior_g, pr.g_ci_lo, pr.g_ci_hi,
                      pr.p_success, outcome, bool(allo), tier, reason)


def score_catalogue(catalogue: pd.DataFrame,
                    exemplars: dict[str, list[_Exemplar]],
                    priors: dict[str, ClassPrior], **kw) -> pd.DataFrame:
    """Score a catalogue DataFrame [id, smiles] and return one row per compound,
    sorted by (non-abstained first, then predicted prior_g desc). The actionable
    output: in-precedent compounds ranked by predicted clinical g."""
    id_col = "id" if "id" in catalogue.columns else catalogue.columns[0]
    smi_col = "smiles" if "smiles" in catalogue.columns else catalogue.columns[1]
    recs = [score_compound(str(r[id_col]), str(r[smi_col]), exemplars, priors, **kw)
            for _, r in catalogue.iterrows()]
    df = pd.DataFrame([r.__dict__ for r in recs])
    df["_abstain"] = (df["tier"] == "ABSTAIN").astype(int)
    df["_g_sort"] = df["prior_g"].fillna(-1e9)
    return (df.sort_values(["_abstain", "_g_sort"], ascending=[True, False])
            .drop(columns=["_abstain", "_g_sort"]).reset_index(drop=True))


# ---------------------------------------------------------------------------
# Leave-one-compound-out class recovery (validation + threshold calibration)
# ---------------------------------------------------------------------------

def loco_class_recovery(ledger: pd.DataFrame, smiles_df: pd.DataFrame, *,
                        tau_ood: float = TAU_OOD,
                        min_class_n: int = MIN_CLASS_N) -> dict:
    """Hold out each SMILES-backed compound, rebuild the exemplar library + priors
    from the rest, and re-route the held-out compound. The honest test of the
    assignment step: does structure recover the TRUE mechanism class?

    Only compounds whose class keeps >= 1 sibling after holdout are evaluable
    (otherwise the true class has no exemplar to match). Returns top-1 accuracy
    among ROUTED (non-abstained) evaluable compounds, the abstention rate, and a
    per-similarity-bin breakdown to calibrate tau_high.
    """
    smi = smiles_df.copy()
    smi["_k"] = smi["compound"].astype(str).str.lower().str.strip()
    have = set(smi["_k"])
    led = ledger[ledger["compound"].astype(str).str.lower().str.strip().isin(have)].copy()

    rows = []
    for i, row in led.iterrows():
        true_cls = str(row["mechanism_class"])
        rest = ledger.drop(index=i)
        # need a sibling left in the true class, else unevaluable
        if (rest["mechanism_class"] == true_cls).sum() < 1:
            continue
        ex = build_exemplars(rest, smi.rename(columns={"_k": "_drop"}))
        pr = build_class_priors(rest)
        k = str(row["compound"]).lower().strip()
        s = dict(zip(smi["_k"], smi["smiles"])).get(k)
        sc = score_compound(str(row["compound"]), str(s), ex, pr,
                            tau_ood=tau_ood, min_class_n=min_class_n)
        rows.append({
            "compound": row["compound"], "true_class": true_cls,
            "assigned": sc.assigned_class, "similarity": sc.similarity,
            "tier": sc.tier, "correct": (sc.assigned_class == true_cls),
            "abstained": (sc.tier == "ABSTAIN"),
        })
    res = pd.DataFrame(rows)
    if res.empty:
        return {"n_evaluable": 0, "n_routed": 0, "top1_acc": float("nan"),
                "abstain_rate": float("nan"), "detail": res}
    routed = res[~res["abstained"]]
    acc = float(routed["correct"].mean()) if len(routed) else float("nan")
    return {
        "n_evaluable": int(len(res)),
        "n_routed": int(len(routed)),
        "top1_acc": acc,
        "abstain_rate": float(res["abstained"].mean()),
        "acc_when_sim_ge_high": float(
            res[(res["similarity"] >= TAU_HIGH) & (~res["abstained"])]["correct"].mean()
        ) if (res["similarity"] >= TAU_HIGH).any() else float("nan"),
        "detail": res.sort_values("similarity", ascending=False).reset_index(drop=True),
    }
