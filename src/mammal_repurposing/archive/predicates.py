"""Machine-checkable keystone predicates: the conditions under which a dead hypothesis re-opens.

A keystone written as prose ("it failed because we had no durability data") cannot be checked, so
nothing ever checks it, so the hypothesis stays dead no matter what the repository learns later.
This module is the difference between an archive and a graveyard: every revivable entry names a
predicate here, and `scripts/129_retro_validation.py` evaluates all of them against CURRENT repo
state on demand.

The predicate language is deliberately tiny: `name(arg, arg)`, resolved against the registry below.
There is no eval, no expression parser, and no way to write a predicate that is not implemented,
because a predicate the engine silently cannot evaluate is worse than no predicate at all -- it
reads as "checked and still dead".

Every checker returns `(satisfied, detail)` and the detail is mandatory. A bare True is not useful
six months later; "PRODRUG_TO_ACTIVE covers 3 compounds, ibogaine is not among them" is.
"""
from __future__ import annotations

import re
from collections.abc import Callable
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
RAW = ROOT / "data" / "raw"

#: name -> (arity, checker). Checkers take string args and return (satisfied, detail).
_REGISTRY: dict[str, tuple[int, Callable[..., tuple[bool, str]]]] = {}


def predicate(name: str, arity: int):
    def deco(fn):
        _REGISTRY[name] = (arity, fn)
        return fn
    return deco


class UnknownPredicate(ValueError):
    """A keystone naming a predicate that does not exist. Fails loudly, by design."""


_CALL = re.compile(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*\((.*)\)\s*$", re.S)


def parse(expr: str) -> tuple[str, list[str]]:
    m = _CALL.match(expr or "")
    if not m:
        raise UnknownPredicate(
            f"keystone predicate {expr!r} is not of the form name(arg, ...). "
            f"Known predicates: {sorted(_REGISTRY)}")
    name, raw = m.group(1), m.group(2).strip()
    args = [a.strip().strip("'\"") for a in raw.split(",")] if raw else []
    return name, [a for a in args if a != ""]


def evaluate(expr: str) -> tuple[bool, str]:
    """Evaluate one keystone predicate against current repository state."""
    name, args = parse(expr)
    if name not in _REGISTRY:
        raise UnknownPredicate(
            f"no checker registered for {name!r}. A predicate the engine cannot evaluate would "
            f"report as 'still dead' forever, so this is an error. Known: {sorted(_REGISTRY)}")
    arity, fn = _REGISTRY[name]
    if arity >= 0 and len(args) != arity:
        raise UnknownPredicate(f"{name} takes {arity} argument(s), got {len(args)}: {args}")
    try:
        return fn(*args)
    except UnknownPredicate:
        raise
    except Exception as e:                       # noqa: BLE001 - a broken checker must not look FALSE
        raise UnknownPredicate(f"{name}({', '.join(args)}) raised {type(e).__name__}: {e}") from e


def _read(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        return None
    try:
        return pd.read_csv(path)
    except Exception:                            # noqa: BLE001
        return None


# -------------------------------------------------------------------------------------------------
# Capability predicates -- "does the pipeline now have a thing it lacked?"
# -------------------------------------------------------------------------------------------------

@predicate("prodrug_resolution_covers", 1)
def _prodrug_covers(compound: str) -> tuple[bool, str]:
    """Is this compound resolved to its active species before structural scoring?

    The keystone for every verdict that scored a parent molecule when the biology is done by a
    metabolite. Psilocybin made this concrete: the L4 window called it negative on the phosphate
    ester's polarity while PERSEUS had been resolving it to psilocin all along.
    """
    from mammal_repurposing.engine.perseus import PRODRUG_TO_ACTIVE
    covered = sorted(PRODRUG_TO_ACTIVE)
    hit = compound.strip().lower() in {c.lower() for c in covered}
    return hit, (f"PRODRUG_TO_ACTIVE covers {len(covered)} compound(s) {covered}; "
                 f"{compound!r} is {'PRESENT' if hit else 'ABSENT'}")


@predicate("metabolite_prediction_available", 0)
def _metabolite_pred() -> tuple[bool, str]:
    """Is there a GENERAL metabolite predictor, as opposed to a hand-written map?

    A hand-curated map protects exactly the compounds someone thought of. The B2 measurement found
    the hole the moment it looked: psilocybin was in the map and ibogaine was not, and ibogaine is
    the one L4 got wrong.
    """
    mod = ROOT / "src" / "mammal_repurposing" / "engine" / "metabolite.py"
    return mod.exists(), (f"{mod.relative_to(ROOT)} "
                          f"{'exists' if mod.exists() else 'does not exist'}; until it does, "
                          "metabolite coverage is whatever PRODRUG_TO_ACTIVE happens to list")


@predicate("dose_axis_exists", 0)
def _dose_axis() -> tuple[bool, str]:
    """Can any ledger represent a dose-response, or only a single effect per compound?

    A compound with an inverted-U response cannot be represented by one effect size. Until some
    ledger carries dose, every hypothesis killed by a single-dose or monotonic design is untestable
    rather than refuted.
    """
    wanted = {"dose", "dose_mg", "dose_mg_kg", "dose_level", "concentration"}
    found = []
    for p in sorted(RAW.glob("*.csv")):
        df = _read(p)
        if df is not None and wanted & set(df.columns):
            found.append(f"{p.name}:{sorted(wanted & set(df.columns))}")
    return bool(found), (f"ledgers carrying a dose column: {found or 'NONE'}")


@predicate("baseline_ability_axis_exists", 0)
def _baseline_axis() -> tuple[bool, str]:
    """Can any ledger represent baseline-dependent direction reversal?

    The inverted-U in healthy adults is usually baseline-dependent: the same dose helps low-baseline
    and harms high-baseline participants. A ledger with one pooled effect per compound averages that
    to approximately zero and reports it as a null.
    """
    wanted = {"baseline_ability", "baseline_strata", "baseline_group", "baseline_percentile"}
    found = [p.name for p in sorted(RAW.glob("*.csv"))
             if (df := _read(p)) is not None and wanted & set(df.columns)]
    return bool(found), f"ledgers carrying a baseline-ability column: {found or 'NONE'}"


# -------------------------------------------------------------------------------------------------
# Evidence predicates -- "has the world (or our curation of it) produced the missing row?"
# -------------------------------------------------------------------------------------------------

@predicate("ci_recorded", 1)
def _ci_recorded(compound: str) -> tuple[bool, str]:
    """Does the healthy-adult ledger now carry an interval for this compound?

    The keystone for every `unknown_precision` kill. Without an interval a null is indistinguishable
    from silence, and four of this ledger's fourteen non-enhancing labels are in exactly that state.
    """
    df = _read(RAW / "healthy_adult_cognition_ledger.csv")
    if df is None:
        return False, "healthy_adult_cognition_ledger.csv not readable"
    row = df[df["compound"].astype(str).str.lower() == compound.strip().lower()]
    if row.empty:
        return False, f"{compound!r} is not in the healthy-adult ledger"
    r = row.iloc[0]
    has = pd.notna(r.get("ci_lo")) and pd.notna(r.get("ci_hi"))
    return has, (f"{compound}: ci_lo={r.get('ci_lo')}, ci_hi={r.get('ci_hi')} -> "
                 f"{'interval recorded' if has else 'STILL no interval'}")


@predicate("interval_excludes_target", 1)
def _interval_excludes(compound: str) -> tuple[bool, str]:
    """Has the evidence become precise enough to actually CLOSE the question?

    The keystone for every `underpowered` kill. It fires only when the interval's upper bound falls
    below the project's target effect, which is the condition for promoting `underpowered` to
    `measured_null`. Note the direction: this predicate firing means the hypothesis is now DEAD for
    good, not alive. Some keystones resolve downward, and recording that is the point.
    """
    from mammal_repurposing.archive.stepping_stone import TARGET_EFFECT_G, classify_null
    df = _read(RAW / "healthy_adult_cognition_ledger.csv")
    if df is None:
        return False, "healthy_adult_cognition_ledger.csv not readable"
    row = df[df["compound"].astype(str).str.lower() == compound.strip().lower()]
    if row.empty:
        return False, f"{compound!r} is not in the healthy-adult ledger"
    r = row.iloc[0]
    cls = classify_null(r.get("ci_lo"), r.get("ci_hi"))
    closed = cls in {"measured_null", "measured_harm"}
    return closed, (f"{compound}: CI [{r.get('ci_lo')}, {r.get('ci_hi')}] vs target "
                    f"g={TARGET_EFFECT_G} -> {cls}")


@predicate("durable_healthy_rows", 1)
def _durable_rows(min_n: str) -> tuple[bool, str]:
    """Does a VERIFIED durable, post-washout, healthy-adult cognitive row exist yet?

    This is the project's binding constraint expressed as a predicate, and the first version of it
    was WRONG in a way worth recording. It accepted any direction matching "positive", which swept
    in `positive_conditional` and `positive_preliminary`, ignored `replicated`, and ignored sample
    size. It therefore returned 3 rows and flagged a hypothesis for revival on the strength of
    Rokem & Silver 2013 (n = 8, unreplicated), Shellshear 2015 (effect unobtainable) and Chamoun
    2017 (n = 9, preliminary) -- the three studies this project had already examined and excluded.
    `ledger_guard.DURABILITY_REQUIRED_FIELDS` exists precisely because that n=8 precedent has a
    placebo arm that also retained its learning.

    A predicate looser than the standard it is checking against manufactures false revivals, and a
    few of those would discredit the whole archive faster than any single wrong verdict. So this
    now applies the project's own criteria: unqualified positive direction, replicated, and clean
    under `validate_durability_claim`.
    """
    from mammal_repurposing.validation.ledger_guard import validate_durability_claim
    n = int(min_n)
    df = _read(RAW / "paired_experience_ledger.csv")
    if df is None:
        return False, "paired_experience_ledger.csv not readable"
    if not {"off_drug_at_readout", "population", "direction"} <= set(df.columns):
        return False, f"ledger lacks durability columns; has {sorted(df.columns)[:8]}..."

    off = pd.to_numeric(df["off_drug_at_readout"], errors="coerce").fillna(0).astype(int) == 1
    healthy = df["population"].astype(str).str.contains("healthy", case=False, na=False)
    # EXACT "positive" only. `positive_conditional` and `positive_preliminary` are the authors'
    # own hedges and must not be promoted to a clean positive by a substring match.
    clean_pos = df["direction"].astype(str).str.strip().str.lower() == "positive"
    cand = df[off & healthy & clean_pos]

    kept, rejected = [], []
    for _, r in cand.iterrows():
        errs = [v for v in validate_durability_claim(r) if v.severity == "error"]
        warns = [v for v in validate_durability_claim(r) if v.severity == "warn"]
        if errs:
            rejected.append(f"{r['compound']}(underspecified)")
        elif warns:
            rejected.append(f"{r['compound']}(unreplicated)")
        else:
            kept.append(str(r["compound"]))
    near = int((off & healthy).sum())
    return len(kept) >= n, (
        f"{len(kept)} row(s) meeting the full durability standard vs required {n}"
        + (f": {sorted(kept)}" if kept else "")
        + f". {near} row(s) are off-drug and healthy at all; {len(cand)} carry an unqualified "
          f"positive direction; rejected on the standard: {rejected or 'none'}")


@predicate("assay_family_signed_compounds", 1)
def _signed_multifamily(min_n: str) -> tuple[bool, str]:
    """Are there enough compounds with a SIGNED direction in >=2 assay families to test sign flips?

    The B2 sign-consistency test ran on four compounds and could not estimate a flip rate. Its kill
    came from the permutation gate instead. This is the keystone that would let the flip question
    actually be answered.
    """
    n = int(min_n)
    df = _read(RAW / "plasticity_window_assays.csv")
    if df is None:
        return False, "plasticity_window_assays.csv not readable"
    signed = df[df["direction"].isin(["opens_window", "closes_window"])]
    per = signed.groupby("compound")["assay_family"].nunique()
    k = int((per >= 2).sum())
    return k >= n, f"{k} compound(s) with a signed direction in >=2 assay families vs required {n}"


@predicate("paired_studies_at_n", 2)
def _paired_at_n(min_n_per_study: str, min_count: str) -> tuple[bool, str]:
    """Are there enough ADEQUATELY SIZED paired drug-plus-training studies to resolve the contrast?

    The B1 kill rested on six studies at n >= 25, and the decisive observation was that the largest
    n among any positive-direction study was 21. A contrast of -0.076 at p = 0.33 on six studies
    does not refute the interaction hypothesis; it fails to resolve it. This is the keystone.
    """
    per_study, want = int(min_n_per_study), int(min_count)
    df = _read(RAW / "paired_experience_ledger.csv")
    if df is None:
        return False, "paired_experience_ledger.csv not readable"
    if not {"paired_experience", "n_randomised"} <= set(df.columns):
        return False, f"ledger lacks paired_experience/n_randomised; has {sorted(df.columns)[:8]}..."
    paired = pd.to_numeric(df["paired_experience"], errors="coerce").fillna(0).astype(int) == 1
    big = pd.to_numeric(df["n_randomised"], errors="coerce").fillna(0) >= per_study
    k = int((paired & big).sum())
    return k >= want, f"{k} paired study/studies at n >= {per_study} vs required {want}"


@predicate("ledger_rows_at_least", 2)
def _ledger_rows(filename: str, min_n: str) -> tuple[bool, str]:
    """Generic size gate on any raw ledger."""
    n = int(min_n)
    df = _read(RAW / filename)
    if df is None:
        return False, f"{filename} not readable"
    return len(df) >= n, f"{filename} has {len(df)} row(s) vs required {n}"


@predicate("tanimoto_excludes_self", 0)
def _tanimoto_self() -> tuple[bool, str]:
    """Does the Tanimoto ranker still read each compound's own ChEMBL record?

    The keystone for the retracted per-target correlations. Unlike most entries in the archive this
    one is already satisfied: the fix landed at 2026-08-24. It is kept, and swept, so the retracted
    numbers stay findable by anyone who meets them in an old document.
    """
    from mammal_repurposing.cluster_a.tanimoto_ranker import TanimotoRankerConfig
    on = bool(TanimotoRankerConfig().exclude_self)
    return on, (f"TanimotoRankerConfig.exclude_self defaults to {on}; the query compound is "
                f"{'excluded from' if on else 'STILL INSIDE'} its own actives set")


@predicate("cluster_d_v2g_attribution_measured", 0)
def _cluster_d_v2g() -> tuple[bool, str]:
    """Is Cluster D's failure attributable to its variant-to-gene step specifically?

    Deliberately UNRESOLVABLE. Cluster D is measured as at chance for predicting trial success, but
    nothing decomposes that failure across its stages, so the question of whether a better
    variant-to-gene mapping would help is open and unmeasured. Returning False here would assert a
    negative nobody established; the sweep reports it as unresolvable and exits non-zero instead.
    """
    raise UnknownPredicate(
        "cluster_d_v2g_attribution_measured() has no source. Cluster D's at-chance performance is "
        "not decomposed by stage anywhere in the repository, so whether variant-to-gene mapping is "
        "the binding constraint within it is unmeasured. Build the stage-wise ablation before any "
        "keystone depends on this.")


@predicate("allosteric_head_beats", 2)
def _allosteric_beats(metric: str, threshold: str) -> tuple[bool, str]:
    """Has the DTI head stopped being blind at allosteric sites?

    Measured state at archive creation: AMPA-PAM AUROC 0.26, AMPA orthosteric 0.09, MMP9 0.42, all
    permutation p > 0.7. This predicate has no automated source yet and reports UNRESOLVED rather
    than False, because reporting False would quietly assert a measurement nobody made.
    """
    raise UnknownPredicate(
        f"allosteric_head_beats({metric}, {threshold}) has no automated source. The measurement "
        "lives in reports/pipeline/allosteric_ltr_v1.md and is not machine-readable. Wire a parser "
        "or emit a metrics JSON before relying on this keystone; until then it must not report a "
        "verdict it has not made.")


def registered() -> list[str]:
    return sorted(_REGISTRY)
