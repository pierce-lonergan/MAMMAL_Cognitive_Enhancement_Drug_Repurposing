"""Gap 4 robustness — scaffold-split evaluation of the allosteric fusion ranker.

Implements `docs/PREREG_ALLOSTERIC_ROBUSTNESS.md` (LOCKED 2026-08-24).

The published Gap-4 number (`reports/pipeline/allosteric_ltr_v1.md`: MAMMAL alone
-0.12 -> fused +0.61 under leave-one-TARGET-out) is leakage-aware ACROSS targets
but not WITHIN one. ChEMBL affinity data at a single target is dominated by
analogue series, and two channels can carry a within-target ranking without any
generalisation:

  1. **Analogue-series leakage between targets.** 39 of 120 generic Murcko
     scaffolds in the labelled set span more than one target, covering 184 of 299
     rows. Leave-one-target-out leaves every one of those in the training set.
  2. **Tanimoto self-match.** The published `tanimoto` feature is max-Tanimoto to
     the target's ChEMBL actives at pChEMBL >= 8.0 — and the query compound is
     itself in that actives set whenever it is an active, so the feature reads
     exactly 1.000 on 143 of the 289 joined rows. No train/test split can remove
     this, because the feature is derived from the test row's own activity
     record; it has to be recomputed.

This module builds four arms over identical evaluation rows so every contrast is
paired:

  A  reference leave-one-target-out (the published configuration)
  B  scaffold-disjoint training (channel 1 removed)
  C  scaffold-disjoint training + scaffold-clean Tanimoto (both removed) [PRIMARY]
  D  scaffold-blocked grouped CV, K=5, targets ignored (secondary, not gated)

and scores six pre-registered feature blocks on each.

This module does NOT modify `allosteric_ltr`; it imports and reuses its
featurisation, its GBM and its `within_target_spearman`, so the reference arm is
the published code path rather than a re-implementation of it.

numpy / pandas / sklearn / RDKit.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from . import allosteric_ltr as A

logger = logging.getLogger(__name__)

# --- Pre-registered constants (section 5 and 3 of the pre-registration) -------
BOOTSTRAP_SEED = 20260824          # fixed in the pre-registration, written as a literal
BOOTSTRAP_B = 2000                 # resamples
BOOTSTRAP_ALPHA = 0.05             # 95% percentile interval
MODEL_SEED = 0                     # train_fusion_ranker seed, unchanged from published
BLOCK_MARGIN = 0.05                # section 3.1 "one block matches the fusion" margin
ACTIVE_PCHEMBL = 8.0               # published TanimotoRankerConfig threshold
FP_RADIUS = 2                      # ECFP4, as published
FP_BITS = 2048

# --- Pre-registered run-time gates (section 6) --------------------------------
MIN_ROWS_PER_FOLD = 4              # the existing min_n
MIN_TEST_SCAFFOLDS_PER_FOLD = 3
MIN_TRAIN_ROWS_PER_FOLD = 100
MIN_LABEL_IQR = 0.5                # log units
MIN_FOLDS_FOR_HEADLINE = 12
MIN_ROWS_FOR_HEADLINE = 150

# --- Feature blocks (section 3) -----------------------------------------------
BOLTZ_COLS = ["boltz_affinity", "boltz_prob", "has_boltz"]

#: block name -> (kind, columns). kind "raw" scores by the feature value itself
#: (matching the published `mammal_only` / `tanimoto_only` conditions, so that a
#: monotone model wrapper cannot change their rank order); kind "gbm" fits
#: `train_fusion_ranker` on that column subset.
FEATURE_BLOCKS: dict[str, tuple[str, list[str]]] = {
    "M_mammal":     ("raw", ["mammal_pkd"]),
    "T_tanimoto":   ("raw", ["tanimoto"]),
    "B_boltz":      ("gbm", list(BOLTZ_COLS)),
    "P_physchem":   ("gbm", list(A.PHYSCHEM_COLS)),
    "F_fusion":     ("gbm", list(A.FUSION_FEATURES)),
    "F_minus_T":    ("gbm", [c for c in A.FUSION_FEATURES if c != "tanimoto"]),
}

BLOCK_LABELS = {
    "M_mammal":   "M — MAMMAL pKd alone (raw)",
    "T_tanimoto": "T — Tanimoto-to-actives alone (raw)",
    "B_boltz":    "B — Boltz block (GBM)",
    "P_physchem": "P — physicochemistry block (GBM)",
    "F_fusion":   "F — full fusion (GBM, 14 features)",
    "F_minus_T":  "F-minus-T — fusion without Tanimoto (GBM, 13 features)",
}

ARM_LABELS = {
    "A_loto":            "Arm A — leave-one-target-out (reference, as published)",
    "B_scaffold_train":  "Arm B — scaffold-disjoint training",
    "C_scaffold_clean":  "Arm C — scaffold-disjoint training + scaffold-clean Tanimoto (PRIMARY)",
    "D_blocked_cv":      "Arm D — scaffold-blocked grouped CV, K=5 (secondary, not gated)",
}


# ---------------------------------------------------------------------------
# Scaffold keying (section 2.1 / 2.2)
# ---------------------------------------------------------------------------

def scaffold_key(smiles: str, row_id: str) -> str:
    """Generic Bemis-Murcko scaffold SMILES, or a UNIQUE sentinel key.

    The scaffold itself is computed exactly as `_murcko_generic` in
    `validation/novel_compound.py` — the repo's existing convention, chosen there
    so two donepezil-like benzylpiperidines match on skeleton. It is deliberately
    coarse (it merges rings differing only in heteroatom identity), and coarse in
    the CONSERVATIVE direction: it removes more training rows than a strict
    analogue-series definition would, so it cannot inflate the surviving score.

    Two edge cases get a *unique* sentinel key, per section 2.2. Unique matters:
    a unique key matches nothing, so the row is never removed from training for
    scaffold reasons and never removes anything from training itself. Bucketing
    all failures under one shared key would over-remove — it would make
    acetylcholine and glutamate "the same series".

      * RDKit cannot parse the SMILES        -> ``UNPARSEABLE::<row_id>``
      * acyclic molecule, empty generic scaffold -> ``ACYCLIC::<row_id>``
    """
    try:
        from rdkit import Chem
        from rdkit.Chem.Scaffolds import MurckoScaffold
    except Exception:  # pragma: no cover - environment without rdkit
        raise RuntimeError("scaffold_key requires RDKit")

    mol = Chem.MolFromSmiles(str(smiles)) if isinstance(smiles, str) else None
    if mol is None:
        return f"UNPARSEABLE::{row_id}"
    try:
        scaf = MurckoScaffold.GetScaffoldForMol(mol)
        generic = MurckoScaffold.MakeScaffoldGeneric(scaf)
        key = Chem.MolToSmiles(generic) if generic is not None else ""
    except Exception:
        return f"UNPARSEABLE::{row_id}"
    if not key:
        return f"ACYCLIC::{row_id}"
    return key


def is_sentinel(key: str) -> bool:
    """True for the unique per-row sentinel keys of section 2.2."""
    return str(key).startswith(("UNPARSEABLE::", "ACYCLIC::"))


def assign_scaffold_keys(df: pd.DataFrame, *, smiles_col: str = "smiles",
                         id_col: str = "inchikey",
                         out_col: str = "scaffold_key") -> pd.DataFrame:
    """Add the scaffold key column, using `id_col` to make sentinels unique."""
    out = df.copy()
    ids = out[id_col] if id_col in out.columns else out.index.astype(str)
    out[out_col] = [scaffold_key(s, str(i)) for s, i in zip(out[smiles_col], ids)]
    return out


# ---------------------------------------------------------------------------
# Folds and training-set construction (section 2.3)
# ---------------------------------------------------------------------------

def target_folds(df: pd.DataFrame, *, min_n: int = MIN_ROWS_PER_FOLD,
                 group: str = "target_uniprot") -> list[str]:
    """The leave-one-TARGET-out folds: one per target with >= min_n labelled rows."""
    counts = df.groupby(group).size()
    return sorted(str(t) for t, n in counts.items() if n >= min_n)


def scaffold_disjoint_train_mask(df: pd.DataFrame, test_idx: pd.Index, *,
                                 group: str = "target_uniprot",
                                 scaffold_col: str = "scaffold_key",
                                 id_col: str = "inchikey") -> pd.Series:
    """Boolean mask over `df` selecting Arm B/C training rows for one fold.

    A row trains iff it is not at the held-out target AND its scaffold key is not
    in the held-out fold's scaffold keys AND its InChIKey is not in the held-out
    fold's InChIKeys.

    The InChIKey clause is stated separately from the scaffold clause because the
    scaffold rule subsumes it only for REAL keys: two rows of the same compound at
    two targets share a scaffold key, unless that key is a per-row sentinel, in
    which case they do not. 59 of 198 compounds appear at more than one target.
    """
    test = df.loc[test_idx]
    held_targets = set(test[group].astype(str))
    held_scaffolds = {k for k in test[scaffold_col] if not is_sentinel(k)}
    held_ids = set(test[id_col].astype(str))

    mask = ~df[group].astype(str).isin(held_targets)
    mask &= ~df[scaffold_col].isin(held_scaffolds)
    mask &= ~df[id_col].astype(str).isin(held_ids)
    return mask


def blocked_cv_assignment(df: pd.DataFrame, *, k: int = 5,
                          scaffold_col: str = "scaffold_key") -> pd.Series:
    """Arm D: partition rows into K scaffold-disjoint folds, deterministically.

    Groups are assigned largest-first to the currently smallest fold (ties broken
    by fold index, then by scaffold key for stable ordering). No RNG, so the
    partition is reproducible without a seed.
    """
    sizes = df[scaffold_col].value_counts()
    order = sorted(sizes.index, key=lambda s: (-int(sizes[s]), str(s)))
    loads = [0] * k
    assign: dict[str, int] = {}
    for key in order:
        f = min(range(k), key=lambda i: (loads[i], i))
        assign[key] = f
        loads[f] += int(sizes[key])
    return df[scaffold_col].map(assign).astype(int)


# ---------------------------------------------------------------------------
# Arm C — scaffold-clean Tanimoto (section 2.3)
# ---------------------------------------------------------------------------

def scaffold_clean_tanimoto(rows: pd.DataFrame, actives: pd.DataFrame, *,
                            group: str = "target_uniprot",
                            smiles_col: str = "smiles",
                            id_col: str = "inchikey",
                            scaffold_col: str = "scaffold_key",
                            active_smiles_col: str = "canonical_smiles",
                            active_id_col: str = "inchikey") -> pd.Series:
    """Recompute max-Tanimoto-to-actives with the self-match channel removed.

    Same fingerprint (ECFP4, Morgan radius 2, 2048 bits), same `max` aggregator,
    same pChEMBL >= 8.0 actives set as the published feature — the ONLY change is
    which actives are eligible. For a query row at target T, the actives set is
    filtered to drop

      (a) the query molecule's own InChIKey — this is what makes the published
          feature read 1.000 on 143 of 289 rows, and
      (b) every active whose generic Murcko scaffold key equals the query's —
          the analogue-series channel inside the target's own actives.

    Both filters depend only on the query row, not on the fold, so the result is a
    fold-independent per-row column. Returns NaN where no eligible active remains
    (imputed downstream on training statistics, like any other missing feature).
    """
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import AllChem
    RDLogger.DisableLog("rdApp.*")

    def _fp(smi):
        if not isinstance(smi, str) or not smi:
            return None
        m = Chem.MolFromSmiles(smi)
        if m is None:
            return None
        return AllChem.GetMorganFingerprintAsBitVect(m, radius=FP_RADIUS, nBits=FP_BITS)

    act = actives.copy()
    act[group] = act[group].astype(str)
    # Scaffold-key the actives once. Sentinels are made unique on the active's own
    # id so an unparseable active can never collide with a query's scaffold key.
    act_keys: list[str] = []
    act_fps: list = []
    for smi, aid in zip(act[active_smiles_col], act[active_id_col].astype(str)):
        act_fps.append(_fp(smi))
        act_keys.append(scaffold_key(smi, f"ACTIVE::{aid}"))
    act = act.assign(_fp=act_fps, _scaf=act_keys)
    act = act[act["_fp"].notna()]

    by_target = {t: g for t, g in act.groupby(group)}

    out = np.full(len(rows), np.nan, dtype=float)
    n_no_actives = 0
    for pos, (_, r) in enumerate(rows.iterrows()):
        g = by_target.get(str(r[group]))
        if g is None or g.empty:
            n_no_actives += 1
            continue
        q_fp = _fp(r[smiles_col])
        if q_fp is None:
            continue
        q_id = str(r[id_col])
        q_scaf = str(r[scaffold_col])
        keep = (g[active_id_col].astype(str) != q_id)
        if not is_sentinel(q_scaf):
            keep &= (g["_scaf"] != q_scaf)
        fps = g.loc[keep, "_fp"].tolist()
        if not fps:
            continue
        sims = DataStructs.BulkTanimotoSimilarity(q_fp, fps)
        out[pos] = float(max(sims)) if sims else np.nan
    if n_no_actives:
        logger.warning("scaffold_clean_tanimoto: %d rows had no actives at their target",
                       n_no_actives)
    return pd.Series(out, index=rows.index, name="tanimoto_clean")


# ---------------------------------------------------------------------------
# Scoring one fold
# ---------------------------------------------------------------------------

def _impute_train_test(train: pd.DataFrame, test: pd.DataFrame,
                       cols: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per-fold imputation on TRAIN statistics only — identical to `loto_evaluate`.

    Train: per-target mean, then train global mean, then 0.
    Test: the train global mean (never the held-out target's own statistics).
    """
    train = train.copy()
    test = test.copy()
    for c in cols:
        if c not in train.columns:
            continue
        train[c] = train.groupby("target_uniprot")[c].transform(lambda s: s.fillna(s.mean()))
        gmean = train[c].mean()
        train[c] = train[c].fillna(gmean).fillna(0.0)
        test[c] = test[c].fillna(gmean).fillna(0.0)
    return train, test


def score_fold(train: pd.DataFrame, test: pd.DataFrame, *, label_col: str = "pact",
               blocks: dict[str, tuple[str, list[str]]] | None = None,
               seed: int = MODEL_SEED) -> pd.DataFrame:
    """Score every feature block on one fold's held-out rows.

    Returns `test` with one `score__<block>` column per block. Imputation happens
    once per fold on train statistics only, then every block reads the same
    imputed frame, so blocks differ only by which columns they see.
    """
    blocks = blocks or FEATURE_BLOCKS
    impute_cols = [c for c in A.FUSION_FEATURES if c != "has_boltz"]
    train, test = _impute_train_test(train, test, impute_cols)

    for name, (kind, cols) in blocks.items():
        if kind == "raw":
            test[f"score__{name}"] = test[cols[0]].astype(float)
        else:
            model, feats = A.train_fusion_ranker(train, label_col, cols, seed=seed)
            test[f"score__{name}"] = model.predict(test[feats].to_numpy(float))
    return test


# ---------------------------------------------------------------------------
# Run-time gates (section 6)
# ---------------------------------------------------------------------------

@dataclass
class FoldGate:
    fold: str
    n_test: int
    n_test_scaffolds: int
    n_train_scaffold_disjoint: int
    label_iqr: float
    passed: bool
    reasons: list[str] = field(default_factory=list)


def evaluate_fold_gates(df: pd.DataFrame, folds: list[str], *,
                        label_col: str = "pact",
                        group: str = "target_uniprot") -> list[FoldGate]:
    """Apply the pre-registered per-fold gates. A fold failing any clause is
    dropped from the primary metric and reported by name with its reason."""
    gates: list[FoldGate] = []
    for t in folds:
        test_idx = df.index[df[group].astype(str) == t]
        test = df.loc[test_idx]
        n_test = len(test)
        n_scaf = test["scaffold_key"].nunique()
        n_train = int(scaffold_disjoint_train_mask(df, test_idx).sum())
        q75, q25 = np.percentile(test[label_col].to_numpy(float), [75, 25])
        iqr = float(q75 - q25)

        reasons = []
        if n_test < MIN_ROWS_PER_FOLD:
            reasons.append(f"n_test {n_test} < {MIN_ROWS_PER_FOLD}")
        if n_scaf < MIN_TEST_SCAFFOLDS_PER_FOLD:
            reasons.append(f"test scaffolds {n_scaf} < {MIN_TEST_SCAFFOLDS_PER_FOLD}")
        if n_train < MIN_TRAIN_ROWS_PER_FOLD:
            reasons.append(f"scaffold-disjoint train rows {n_train} < {MIN_TRAIN_ROWS_PER_FOLD}")
        if iqr < MIN_LABEL_IQR:
            reasons.append(f"label IQR {iqr:.2f} < {MIN_LABEL_IQR}")
        gates.append(FoldGate(fold=t, n_test=n_test, n_test_scaffolds=n_scaf,
                              n_train_scaffold_disjoint=n_train, label_iqr=iqr,
                              passed=not reasons, reasons=reasons))
    return gates


# ---------------------------------------------------------------------------
# Arms
# ---------------------------------------------------------------------------

def run_arm(df: pd.DataFrame, folds: list[str], arm: str, *,
            label_col: str = "pact", tanimoto_col: str | None = None,
            blocks: dict[str, tuple[str, list[str]]] | None = None,
            seed: int = MODEL_SEED, k_blocked: int = 5) -> pd.DataFrame:
    """Run one arm and return the pooled held-out rows with per-block scores.

    `df` is the FULL labelled pool that training may draw on; `folds` are the
    targets to evaluate. These are deliberately separate: the section-6 gates drop
    a fold from the primary METRIC, not from other folds' training sets, and the
    published `loto_evaluate` likewise trains on every row whose target differs —
    including rows at targets too small to form a fold of their own. Held-out rows
    are then restricted to `folds`, so every arm reports on identical evaluation
    rows.

    `tanimoto_col`, when given, replaces the `tanimoto` feature column before any
    model is fit (Arm C's scaffold-clean recomputation). Both train and test rows
    get the replacement, per section 2.3.
    """
    work = df.copy()
    if tanimoto_col is not None:
        work["tanimoto"] = work[tanimoto_col].astype(float)
    keep = set(str(t) for t in folds)

    held: list[pd.DataFrame] = []
    if arm == "D_blocked_cv":
        assign = blocked_cv_assignment(work, k=k_blocked)
        for f in range(k_blocked):
            test_idx = work.index[assign == f]
            train_idx = work.index[assign != f]
            if len(test_idx) == 0 or len(train_idx) == 0:
                continue
            out = score_fold(work.loc[train_idx], work.loc[test_idx],
                             label_col=label_col, blocks=blocks, seed=seed)
            out["cv_fold"] = f
            # Score only the evaluation rows the gated arms use.
            held.append(out[out["target_uniprot"].astype(str).isin(keep)])
    else:
        for t in folds:
            test_idx = work.index[work["target_uniprot"].astype(str) == t]
            if arm == "A_loto":
                train_idx = work.index[work["target_uniprot"].astype(str) != t]
            elif arm in ("B_scaffold_train", "C_scaffold_clean"):
                train_idx = work.index[scaffold_disjoint_train_mask(work, test_idx)]
            else:
                raise ValueError(f"unknown arm: {arm}")
            out = score_fold(work.loc[train_idx], work.loc[test_idx],
                             label_col=label_col, blocks=blocks, seed=seed)
            out["cv_fold"] = t
            out["n_train_rows"] = len(train_idx)
            held.append(out)

    if not held:
        return work.iloc[0:0]
    allheld = pd.concat(held, ignore_index=True)
    allheld["arm"] = arm
    return allheld


def per_fold_rho(held: pd.DataFrame, block: str, *, label_col: str = "pact",
                 group: str = "target_uniprot",
                 folds: list[str] | None = None) -> pd.DataFrame:
    """Per-target Spearman and fold size for one block, via the published
    `within_target_spearman` (mid-ranks for ties).

    When `folds` is given, the table is ALIGNED to exactly that fold list, with
    `rho = NaN` on folds where the block's score is constant and Spearman is
    therefore undefined. Aligning rather than dropping matters: it keeps the fold
    universe fixed at the section-6 gate survivors and lets each block pool over
    the folds where it is defined, instead of letting the least-informative block
    (Boltz, degenerate wherever a target has no Boltz coverage) silently shrink
    the fold set for every other block and for the primary metric.
    """
    _, per = A.within_target_spearman(held, f"score__{block}", label_col, group=group)
    n = held.groupby(group).size()
    keys = list(folds) if folds is not None else sorted(per)
    return pd.DataFrame({
        "fold": keys,
        "rho": [per.get(str(t), float("nan")) for t in keys],
        "n": [int(n[t]) if t in n.index else 0 for t in keys],
    }).reset_index(drop=True)


def pooled_rho(fold_table: pd.DataFrame) -> float:
    """Sample-size-weighted pooled rho over the folds where it is DEFINED.

    The primary metric of section 4. NaN folds (undefined Spearman) are excluded
    from both the numerator and the weights; `n_contributing_folds` reports how
    many actually contributed.
    """
    if fold_table.empty:
        return float("nan")
    rho = fold_table["rho"].to_numpy(float)
    n = fold_table["n"].to_numpy(float)
    ok = np.isfinite(rho) & (n > 0)
    if not ok.any():
        return float("nan")
    return float(np.average(rho[ok], weights=n[ok]))


def n_contributing_folds(fold_table: pd.DataFrame) -> int:
    """Folds on which this block's rho is defined — reported per cell so a block
    pooled over fewer folds is never mistaken for one pooled over all of them."""
    if fold_table.empty:
        return 0
    return int((np.isfinite(fold_table["rho"].to_numpy(float))
                & (fold_table["n"].to_numpy(float) > 0)).sum())


# ---------------------------------------------------------------------------
# Cluster bootstrap (section 5)
# ---------------------------------------------------------------------------

def bootstrap_indices(n_folds: int, *, b: int = BOOTSTRAP_B,
                      seed: int = BOOTSTRAP_SEED) -> np.ndarray:
    """The resampled fold-index sets, drawn ONCE and reused across every arm and
    every feature block so all contrasts are paired on the same draws."""
    rng = np.random.default_rng(seed)
    return rng.integers(0, n_folds, size=(b, n_folds))


def bootstrap_pooled(fold_table: pd.DataFrame, idx: np.ndarray) -> np.ndarray:
    """Pooled rho under each resample of the fold-level (rho, n) pairs.

    `idx` indexes the SAME aligned fold list for every arm and every block, so
    contrasts are paired on identical draws (section 5). Folds where this block's
    rho is undefined contribute zero weight within a resample rather than voiding
    it; a resample in which no fold is defined yields NaN.

    Declared limitation (section 5): models are NOT refit inside the resample, so
    this reflects between-target variability in the held-out scores, not
    model-refitting variability.
    """
    rho = fold_table["rho"].to_numpy(float)
    n = fold_table["n"].to_numpy(float)
    r = rho[idx]
    w = n[idx].astype(float)
    ok = np.isfinite(r) & (w > 0)
    w = np.where(ok, w, 0.0)
    r = np.where(ok, r, 0.0)
    tot = w.sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = (r * w).sum(axis=1) / tot
    return np.where(tot > 0, out, np.nan)


def percentile_ci(draws: np.ndarray, alpha: float = BOOTSTRAP_ALPHA) -> tuple[float, float]:
    """95% percentile interval — the only interval the pre-registration reports."""
    if draws is None or not np.isfinite(np.asarray(draws, dtype=float)).any():
        return float("nan"), float("nan")
    lo, hi = np.nanpercentile(draws, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return float(lo), float(hi)


# ---------------------------------------------------------------------------
# Band verdict (section 4)
# ---------------------------------------------------------------------------

def band_verdict(rho_primary: float, ci_lo: float, rho_reference: float) -> tuple[str, str]:
    """The section-4 bands, evaluated in their own words.

    SURVIVES  rho >= 0.45 AND bootstrap 95% LB >= 0.25 AND rho >= rho_ref - 0.15
    DEGRADES  0.25 <= rho < 0.45, OR rho < rho_ref - 0.15 while rho >= 0.25
    FAILS     rho < 0.25, OR bootstrap 95% LB <= 0
    """
    drop = rho_reference - rho_primary
    if rho_primary < 0.25 or ci_lo <= 0.0:
        why = []
        if rho_primary < 0.25:
            why.append(f"rho_C(F) = {rho_primary:+.3f} < 0.25")
        if ci_lo <= 0.0:
            why.append(f"bootstrap 95% lower bound = {ci_lo:+.3f} <= 0")
        return "FAILS", " and ".join(why)
    if rho_primary >= 0.45 and ci_lo >= 0.25 and drop <= 0.15:
        return "SURVIVES", (f"rho_C(F) = {rho_primary:+.3f} >= 0.45, 95% LB = {ci_lo:+.3f} >= 0.25, "
                            f"drop from Arm A = {drop:+.3f} <= 0.15")
    why = []
    if 0.25 <= rho_primary < 0.45:
        why.append(f"0.25 <= rho_C(F) = {rho_primary:+.3f} < 0.45")
    if drop > 0.15:
        why.append(f"drop from Arm A = {drop:+.3f} > 0.15")
    if rho_primary >= 0.45 and ci_lo < 0.25:
        why.append(f"95% LB = {ci_lo:+.3f} < 0.25")
    return "DEGRADES", " and ".join(why) if why else f"rho_C(F) = {rho_primary:+.3f}"


def block_interpretation(pooled: dict[str, float], *, margin: float = BLOCK_MARGIN) -> list[str]:
    """Section 3.1 — the pre-registered reading if a single block matches the fusion."""
    notes: list[str] = []
    f = pooled.get("F_fusion", float("nan"))
    names = {"T_tanimoto": "structural similarity to known actives",
             "P_physchem": "physicochemistry alone",
             "B_boltz": "the Boltz block alone"}
    for blk, phrase in names.items():
        v = pooled.get(blk, float("nan"))
        if not np.isfinite(v) or not np.isfinite(f):
            continue
        if v > f:
            notes.append(f"{blk} ({v:+.3f}) EXCEEDS the fusion ({f:+.3f}) — "
                         f"the fusion is a net negative against {phrase}.")
        elif v >= f - margin:
            notes.append(f"{blk} ({v:+.3f}) matches the fusion ({f:+.3f}) within the "
                         f"pre-registered margin of {margin:.2f} — "
                         f"{phrase} ranks within-target affinity as well as the fusion.")
    fm = pooled.get("F_minus_T", float("nan"))
    if np.isfinite(f) and np.isfinite(fm) and (f - fm) > 0.15:
        notes.append(f"rho(F) - rho(F-minus-T) = {f - fm:+.3f} > 0.15 — the fusion's "
                     f"skill depends on the single Tanimoto feature by that much.")
    return notes


def _spearman_ordinal(a: np.ndarray, b: np.ndarray) -> float:
    """The PRE-2026-06-06 Spearman: ordinal rank (argsort-of-argsort), which
    breaks ties by array order instead of averaging them.

    Kept only to quantify what that bug did to the published numbers. It was
    replaced by mid-ranks in commit 615bb1c, AFTER
    `reports/pipeline/allosteric_ltr_v1.md` was generated (2026-05-30), and that
    report was never regenerated. Not used anywhere in the analysis itself.
    """
    if len(a) < 3:
        return float("nan")

    def rank(x):
        o = x.argsort(kind="mergesort")
        r = np.empty_like(o, dtype=float)
        r[o] = np.arange(len(x))
        return r

    ra, rb = rank(np.asarray(a, float)), rank(np.asarray(b, float))
    ra -= ra.mean()
    rb -= rb.mean()
    d = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / d) if d > 0 else float("nan")


def tie_convention_contrast(held: pd.DataFrame, blocks: list[str], *,
                            label_col: str = "pact",
                            group: str = "target_uniprot") -> dict[str, dict[str, float]]:
    """Pooled rho for each block under mid-rank vs ordinal-rank Spearman.

    Quantifies how much of a published number came from tie handling. It matters
    here because the `tanimoto` feature is extremely tie-heavy — the self-match
    puts 143 of 289 rows at exactly 1.000 — and ordinal ranking scatters those
    ties by array order, destroying rank information the mid-rank convention
    keeps.

    A second, subtler defect of the ordinal convention shows up here too: on a
    fold whose score column is CONSTANT (e.g. a feature that is missing at that
    target and has been fully imputed to one training mean), argsort still hands
    out n distinct ranks, so the old code returned a finite correlation
    manufactured purely from row order where the correct answer is "undefined".
    Both poolings are therefore reported: `*_own_folds` pools each convention over
    the folds IT calls defined — reproducing how the published number was
    computed — while `*_paired` restricts both to folds where both are defined.
    """
    out: dict[str, dict[str, float]] = {}
    for blk in blocks:
        col = f"score__{blk}"
        if col not in held.columns:
            continue
        mid: list = []
        ordi: list = []
        pair: list = []
        for _, g in held.groupby(group):
            if len(g) < 3:
                continue
            s = g[col].to_numpy(float)
            y = g[label_col].to_numpy(float)
            m, o = A._spearman(s, y), _spearman_ordinal(s, y)
            n = float(len(g))
            if np.isfinite(m):
                mid.append((m, n))
            if np.isfinite(o):
                ordi.append((o, n))
            if np.isfinite(m) and np.isfinite(o):
                pair.append((m, o, n))

        def _avg(pairs, i=0):
            if not pairs:
                return float("nan")
            v = np.array([p[i] for p in pairs], float)
            w = np.array([p[-1] for p in pairs], float)
            return float(np.average(v, weights=w))

        out[blk] = {
            "mid_rank_own_folds": _avg(mid),
            "ordinal_rank_own_folds": _avg(ordi),
            "mid_rank_paired": _avg(pair, 0),
            "ordinal_rank_paired": _avg(pair, 1),
            "n_folds_mid": len(mid),
            "n_folds_ordinal": len(ordi),
            "n_folds_paired": len(pair),
        }
    return out


def availability() -> dict:
    try:
        import sklearn  # noqa: F401
        sk = True
    except Exception:
        sk = False
    try:
        import rdkit  # noqa: F401
        rd = True
    except Exception:
        rd = False
    return {"available": sk and rd, "sklearn": sk, "rdkit": rd,
            "arms": list(ARM_LABELS), "blocks": list(FEATURE_BLOCKS),
            "bootstrap_seed": BOOTSTRAP_SEED, "bootstrap_b": BOOTSTRAP_B}
