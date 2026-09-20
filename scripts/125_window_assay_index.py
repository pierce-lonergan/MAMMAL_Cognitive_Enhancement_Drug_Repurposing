"""B2 - is "opens a plasticity window" a property of a COMPOUND, or of a (compound, assay) pair?

The L4 window (engine/psychoplastogen.py) scores compounds. Scoring compounds presumes the window
is a compound-level property. Donepezil is the reason to doubt that: in healthy adults it
AUGMENTS motion-direction perceptual learning (Rokem & Silver 2010, PMID 20850321) and REDUCES the
ocular-dominance shift after monocular deprivation (Sheynin 2019, PMID 30766471, t(11) = -4.9,
p < 0.001). One drug, opposite signs by readout.

CORRECTED 2026-09-20. This file previously attributed BOTH readouts to PMID 30766471 and said they
were measured "in the same people". Verified against the primary source: PMID 30766471 reports ONLY
ocular dominance and contains no learning task. The four donepezil studies come from two
laboratories, and the only pair sharing participants is Rokem & Silver 2010 <-> 2013 (8 of the
original 12). The contrast is also partly confounded with DOSING: both positives used 8-day
steady-state dosing, while the ocular-dominance reduction and the texture-discrimination null both
used a single dose, which Sheynin notes explicitly. See docs/PREREG_DEVIATIONS_2026-06.md.

This script measures it against `data/raw/plasticity_window_assays.csv`, a citation-verified index
of (compound, assay-family, direction) outcomes. Two pre-registered criteria, both decided here:

  SUCCESS  >= 40 rows across >= 3 assay families, >= 8 compounds in >= 2 families, AND L4 agreement
           beats a size-matched permutation gate (p < 0.05) in at least one family.
  KILL     among multi-assay compounds, window-label sign inconsistent across families at >= 30%
           -> DEMOTE L4 from a ranking signal to a per-assay annotation, and stop scoring compounds
           on it. Also kill if L4 fails the permutation gate in EVERY family.

The script never invents a row. If the literature does not contain 40 verified rows, it reports the
number it found and declares SUCCESS unmet; a short honest table is the result, not a padded one.

Writes reports/pipeline/window_assay_index_v1.md. CPU; network only for SMILES (PubChem, cached).
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("b2_assay_index")
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

INDEX_CSV = ROOT / "data" / "raw" / "plasticity_window_assays.csv"
SMILES_CACHE = ROOT / "data" / "raw" / "plasticity_window_smiles.csv"
REPORT = ROOT / "reports" / "pipeline" / "window_assay_index_v1.md"

REQUIRED_COLS = ("compound", "assay_family", "assay", "species", "direction",
                 "effect_detail", "citation_short", "pmid_doi", "verification")
#: only these two are DIRECTIONAL. no_effect and mixed carry no sign and cannot vote on consistency.
SIGNED = {"opens_window": 1, "closes_window": -1}

MIN_ROWS, MIN_FAMILIES = 40, 3
MIN_MULTI_COMPOUNDS, MULTI_FAMILY_MIN = 8, 2
INCONSISTENCY_KILL = 0.30
SEED, N_PERM = 42, 20000


def load_index() -> pd.DataFrame:
    if not INDEX_CSV.exists():
        raise SystemExit(f"missing {INDEX_CSV}; B2 has no data to measure and must not guess")
    df = pd.read_csv(INDEX_CSV, keep_default_na=False, na_values=[""])
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise SystemExit(f"{INDEX_CSV} missing required columns: {missing}")
    # FAIL CLOSED on provenance. An unsourced plasticity row is exactly the failure mode this
    # project treats as fatal, and a silent drop would hide it.
    unsourced = df[df["pmid_doi"].fillna("").astype(str).str.strip() == ""]
    if len(unsourced):
        raise SystemExit("rows without a PMID/DOI:\n" + unsourced[["compound", "assay"]].to_string())
    bad = set(df["direction"]) - (set(SIGNED) | {"no_effect", "mixed"})
    if bad:
        raise SystemExit(f"unknown direction values: {sorted(bad)}")
    return df


def sign_consistency(df: pd.DataFrame) -> tuple[pd.DataFrame, float]:
    """Per compound: does its window SIGN agree across the assay families that measured it?

    A compound counts as INCONSISTENT when its signed directions disagree across families. Within a
    family, disagreement between studies is ordinary replication noise and is resolved by majority;
    ACROSS families it is the Sheynin phenomenon, which is what B2 is testing for.
    """
    rows = []
    for cmpd, g in df[df["direction"].isin(SIGNED)].groupby("compound"):
        per_family = {}
        for fam, gf in g.groupby("assay_family"):
            tot = int(gf["direction"].map(SIGNED).sum())
            if tot == 0:          # family is internally tied -> it casts no vote
                continue
            per_family[fam] = 1 if tot > 0 else -1
        if len(per_family) < MULTI_FAMILY_MIN:
            continue
        signs = set(per_family.values())
        detail = "; ".join(f"{k}:{'+' if v > 0 else '-'}" for k, v in sorted(per_family.items()))
        rows.append(dict(compound=cmpd, n_families=len(per_family), families=detail,
                         inconsistent=int(len(signs) > 1)))
    out = pd.DataFrame(rows)
    rate = float(out["inconsistent"].mean()) if len(out) else float("nan")
    return out, rate


def _smiles_for(names: list[str]) -> dict[str, str]:
    cache: dict[str, str] = {}
    if SMILES_CACHE.exists():
        c = pd.read_csv(SMILES_CACHE)
        cache = dict(zip(c["compound"], c["smiles"]))
    todo = [n for n in names if n not in cache]
    if todo:
        from mammal_repurposing.fetchers.pubchem import fetch_smiles
        for n in todo:
            try:
                r = fetch_smiles(n.replace("_", " "))
                if r and r.get("smiles"):
                    cache[n] = r["smiles"]
            except Exception as e:                      # noqa: BLE001 - network is best-effort
                L.warning("SMILES lookup failed for %s: %s", n, e)
        SMILES_CACHE.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"compound": list(cache), "smiles": [cache[k] for k in cache]}
                     ).to_csv(SMILES_CACHE, index=False)
    return cache


def l4_agreement(df: pd.DataFrame) -> pd.DataFrame:
    """Per assay family: does the L4 structural window predict the EMPIRICAL direction?

    Returns (per-family summary, per-compound call sheet). The call sheet is what makes a
    KILL actionable: a verdict that only says "it does not work" cannot be acted on, and the
    individual misses turn out to share one mechanism.

    Truth = opens_window (1) vs closes_window/no_effect (0). Prediction = the structural call,
    scoped to that family. The permutation gate shuffles the truth labels within the family, so it
    asks whether L4 beats chance AT THE FAMILY'S OWN BASE RATE.
    """
    from mammal_repurposing.engine.perseus import PRODRUG_TO_ACTIVE
    from mammal_repurposing.engine.psychoplastogen import (
        WINDOW_ASSAY_FAMILIES, psychoplastogen_window)
    smi = _smiles_for(sorted(df["compound"].unique()))
    # Score the engine AS DEPLOYED. PERSEUS resolves a curated set of prodrugs to their active
    # species before calling L4, so measuring L4 on the PubChem parent SMILES would invent a false
    # negative the engine does not actually make (psilocybin is in that map). The map is also the
    # measurement: whatever it does NOT cover is a real miss, and ibogaine is not in it.
    prodrug = {k: v[1] for k, v in PRODRUG_TO_ACTIVE.items()}
    resolved = {c: prodrug.get(c, s) for c, s in smi.items()}
    rng = np.random.default_rng(SEED)
    out, calls = [], []
    for fam, g in df.groupby("assay_family"):
        if fam not in WINDOW_ASSAY_FAMILIES:
            L.warning("assay family %r is outside the engine vocabulary; skipped", fam)
            continue
        # one row per compound: a compound that opened the window in ANY study of this family
        # counts as an opener (the screen's job is to FIND openers, so recall is the live question)
        agg = (g.assign(opened=(g["direction"] == "opens_window").astype(int))
                .groupby("compound")["opened"].max())
        y, yhat = [], []
        for cmpd, truth in agg.items():
            s = resolved.get(cmpd)
            if not s:
                calls.append(dict(assay_family=fam, compound=cmpd, truth=int(truth),
                                  l4=None, scaffold=None, tpsa=None, hbd=None,
                                  prodrug_resolved=0))
                continue
            call = psychoplastogen_window(s, assay=fam)
            calls.append(dict(assay_family=fam, compound=cmpd, truth=int(truth),
                              l4=int(call.window), scaffold=call.scaffold,
                              tpsa=call.tpsa, hbd=call.hbd,
                              prodrug_resolved=int(cmpd in prodrug)))
            y.append(int(truth))
            yhat.append(int(call.window))
        y, yhat = np.asarray(y), np.asarray(yhat)
        n = len(y)
        rec = dict(assay_family=fam, n_compounds=n, n_openers=int(y.sum()),
                   n_l4_positive=int(yhat.sum()))
        constant_pred = bool(n) and int(yhat.sum()) in (0, n)
        constant_truth = bool(n) and int(y.sum()) in (0, n)
        if n < 4 or constant_pred or constant_truth:
            # DEGENERATE. A constant prediction or a constant truth makes agreement uninformative
            # and the permutation p meaningless. Say so rather than emit a number.
            note = ("L4 predicts a CONSTANT in this family" if constant_pred
                    else "truth is constant" if constant_truth else "too few compounds")
            rec.update(agreement=float((y == yhat).mean()) if n else float("nan"),
                       perm_p=float("nan"), note=note)
        else:
            obs = float((y == yhat).mean())
            null = np.array([(rng.permutation(y) == yhat).mean() for _ in range(N_PERM)])
            rec.update(agreement=obs,
                       perm_p=float(((null >= obs).sum() + 1) / (N_PERM + 1)), note="")
        out.append(rec)
    return pd.DataFrame(out), pd.DataFrame(calls)


def main() -> int:
    df = load_index()
    verified = df[df["verification"].astype(str).str.upper().str.startswith("CONFIRM")]
    L.info("index: %d rows | %d verified | %d families | %d compounds",
           len(df), len(verified), df["assay_family"].nunique(), df["compound"].nunique())

    cons, rate = sign_consistency(df)
    agree, calls = l4_agreement(df)

    # TWO different counts, and conflating them would overstate the evidence. `n_in2` is how
    # many compounds appear in >= 2 families at all (the pre-registered size bar). `n_multi`
    # is how many have a SIGNED direction in >= 2 families - only those can vote on a sign
    # flip, and it is the smaller, binding number.
    fam_per_compound = df.groupby('compound')['assay_family'].nunique()
    n_in2 = int((fam_per_compound >= MULTI_FAMILY_MIN).sum())
    n_multi = len(cons)
    gate_size = bool(len(df) >= MIN_ROWS and df["assay_family"].nunique() >= MIN_FAMILIES
                     and n_in2 >= MIN_MULTI_COMPOUNDS)
    testable = agree[agree["perm_p"].notna()] if len(agree) else agree
    gate_perm = bool(len(testable) and (testable["perm_p"] < 0.05).any())
    kill_inconsistent = bool(n_multi and rate >= INCONSISTENCY_KILL)
    kill_perm = bool(len(testable) and not gate_perm)

    verdict = ("KILL - demote L4 to a per-assay annotation" if (kill_inconsistent or kill_perm)
               else "SUCCESS" if (gate_size and gate_perm) else "UNMET - insufficient evidence")
    L.info("multi-family compounds %d | inconsistency %.2f | verdict %s", n_multi, rate, verdict)
    write_report(df, verified, cons, rate, agree, calls, testable, n_in2, n_multi,
                 gate_size, gate_perm, kill_inconsistent, kill_perm, verdict)
    return 0


def write_report(df, verified, cons, rate, agree, calls, testable, n_in2, n_multi,
                 gate_size, gate_perm, kill_inconsistent, kill_perm, verdict) -> None:
    fams = df["assay_family"].nunique()
    out = [
        "# B2 - the plasticity window is measured, not owned", "",
        f"**Verdict: {verdict}**", "",
        "The L4 psychoplastogen window scores COMPOUNDS. That presumes window-opening is a property "
        "a compound has. Donepezil is the reason to doubt it: in healthy adults it AUGMENTS "
        "motion-direction perceptual learning (Rokem & Silver 2010, PMID 20850321, n = 12, 8-day "
        "steady-state dosing) and REDUCES the ocular-dominance shift after monocular deprivation "
        "(Sheynin 2019, PMID 30766471, t(11) = -4.9, p < 0.001). It is NULL for letter "
        "identification (Levi 2020, PMID 32347910, no placebo arm) and NULL for texture "
        "discrimination (Byrne 2020, PMID 32511666, placebo-controlled, single dose).", "",
        "Two qualifications, both material. No two of those studies share participants except "
        "Rokem & Silver 2010 and 2013, and Sheynin is a different laboratory; an earlier version of "
        "this report wrongly said all readouts came from the same people. And the contrast is "
        "partly confounded with DOSING, since both positives used 8-day steady-state dosing while "
        "the ocular-dominance reduction and the texture null used a single dose. What strengthens "
        "it is that one author (Silver) is on the positive and on both nulls, so the assay contrast "
        "is partly within-laboratory.", "",
        "## Coverage", "",
        f"- rows: **{len(df)}** (pre-registered success bar {MIN_ROWS})",
        f"- assay families: **{fams}** (bar {MIN_FAMILIES})",
        f"- compounds: **{df['compound'].nunique()}**",
        f"- compounds measured in >= {MULTI_FAMILY_MIN} families: **{n_in2}** "
        f"(bar {MIN_MULTI_COMPOUNDS})",
        f"- ...of which carry a SIGNED direction in >= {MULTI_FAMILY_MIN} families, and can "
        f"therefore vote on a sign flip: **{n_multi}**. A `no_effect` or `mixed` result is "
        "real evidence but has no sign, so the flip test runs on the smaller set.",
        f"- rows CONFIRMED by independent verification: **{len(verified)}/{len(df)}**", "",
        "Rows per family:", "",
        "| assay family | rows | compounds |", "| --- | ---: | ---: |",
    ]
    for fam, g in df.groupby("assay_family"):
        out.append(f"| {fam} | {len(g)} | {g['compound'].nunique()} |")
    out += ["", "## Sign consistency across assay families", ""]
    if n_multi:
        out += [f"Among the **{n_multi}** compounds with a directional result in two or more assay "
                f"families, **{int(cons['inconsistent'].sum())}** flip sign between families: "
                f"**{rate:.0%}** (kill threshold {INCONSISTENCY_KILL:.0%}).", "",
                "| compound | families | sign per family | inconsistent |",
                "| --- | ---: | --- | :---: |"]
        ordered = cons.sort_values(["inconsistent", "compound"], ascending=[False, True])
        for _, r in ordered.iterrows():
            flag = "**YES**" if r["inconsistent"] else "no"
            out.append(f"| {r['compound']} | {r['n_families']} | {r['families']} | {flag} |")
    else:
        out.append("No compound has a directional result in two or more assay families. The "
                   "cross-assay question is therefore **not answerable from the published "
                   "literature as indexed here** - which is itself the finding: the field almost "
                   "never measures two plasticity readouts in the same subjects, so nearly every "
                   "apparent sign flip is a between-study comparison carrying between-study "
                   "confounds, dose regimen chief among them.")
    out += ["", "## Does the L4 structural window predict the empirical direction?", "",
            "Truth = the compound opened the window in at least one study of that family. "
            "Prediction = the structural call scoped to that family. The permutation gate shuffles "
            "truth within the family, so it tests L4 against that family's own base rate.", "",
            "| assay family | compounds | openers | L4-positive | agreement | perm p | note |",
            "| --- | ---: | ---: | ---: | ---: | ---: | --- |"]
    for _, r in agree.iterrows():
        pp = "n/a" if pd.isna(r["perm_p"]) else f"{r['perm_p']:.3f}"
        ag = "n/a" if pd.isna(r["agreement"]) else f"{r['agreement']:.2f}"
        out.append(f"| {r['assay_family']} | {r['n_compounds']} | {r['n_openers']} | "
                   f"{r['n_l4_positive']} | {ag} | {pp} | {r['note']} |")
    # The call sheet for the ONE family the rule was derived from. A KILL that cannot say WHY is
    # not actionable, and here the individual misses turn out to share a single mechanism.
    from mammal_repurposing.engine.psychoplastogen import L4_VALIDATED_ASSAY
    home = calls[calls["assay_family"] == L4_VALIDATED_ASSAY] if len(calls) else calls
    if len(home):
        out += ["", f"## Call sheet on the home family (`{L4_VALIDATED_ASSAY}`)", "",
                "This is the family the structural rule was derived from and the only one that was "
                "statistically testable. Every row is a compound the rule had an opinion about.", "",
                "| compound | empirically opens | L4 says | scaffold | TPSA | HBD | |",
                "| --- | :---: | :---: | --- | ---: | ---: | --- |"]
        for _, r in home.sort_values("compound").iterrows():
            if pd.isna(r["l4"]):
                verdict_cell, mark = "no SMILES", "excluded"
            else:
                verdict_cell = str(int(r["l4"]))
                mark = "HIT" if int(r["l4"]) == int(r["truth"]) else "**MISS**"
            tpsa = "" if pd.isna(r["tpsa"]) else f"{r['tpsa']:.1f}"
            hbd = "" if pd.isna(r["hbd"]) else f"{int(r['hbd'])}"
            out.append(f"| {r['compound']} | {int(r['truth'])} | {verdict_cell} | "
                       f"{r['scaffold'] or '-'} | {tpsa} | {hbd} | {mark} |")
        scored = home[home["l4"].notna()]
        miss = scored[scored["l4"].astype(int) != scored["truth"].astype(int)]
        blind = miss[miss["scaffold"].isna()]
        inscope = miss[miss["scaffold"].notna()]
        n_pro = int(scored["prodrug_resolved"].sum()) if "prodrug_resolved" in scored else 0
        out += ["", f"Of the **{len(miss)}** misses, **{len(blind)}** are out-of-scope blindness "
                f"({', '.join(sorted(blind['compound'])) or 'none'}): no serotonergic scaffold, so "
                "the screen was never built to see them, and counting those as errors would be "
                "scoring it on a job it does not claim. The remaining "
                f"**{len(inscope)}** ({', '.join(sorted(inscope['compound'])) or 'none'}) "
                f"{'is' if len(inscope) == 1 else 'are'} inside its own scaffold class.", "",
                "**The parent-compound assumption, and the hand-curated patch over it.** L4 reads "
                "the structure it is handed, and the plasticity is produced by whatever the liver "
                "hands the brain. PERSEUS already knows this: it resolves a curated prodrug map "
                f"before calling L4, and this table is scored THROUGH that map ({n_pro} of "
                f"{len(scored)} compounds here were resolved). Psilocybin is why the map exists - "
                "the phosphate ester fails the permeability gate at TPSA 86 while psilocin passes "
                "at TPSA 39, and without the map L4 would call a real opener negative.", "",
                "The map is three compounds long and written by hand, so its coverage is the whole "
                "of the protection. Ibogaine is not in it. Ly 2018 measured NO structural effect "
                "for ibogaine itself and identified noribogaine as the active species, but "
                "ibogaine clears the gate on its own descriptors (TPSA 28, HBD 1), so L4 calls a "
                "measured null positive. That is the single in-scope error in the family, and it "
                "is not a threshold that can be moved: the fix is either another hand-curated map "
                "entry, which does not generalise, or metabolite prediction, which the pipeline "
                "does not have.", ""]
    out += ["", "## Gate arithmetic", "",
            f"- size bar met: **{gate_size}**",
            f"- beats permutation in >= 1 family: **{gate_perm}** "
            f"({len(testable)} of {len(agree)} families were statistically testable)",
            f"- kill by sign inconsistency (>= {INCONSISTENCY_KILL:.0%}): **{kill_inconsistent}**",
            f"- kill by failing permutation in EVERY testable family: **{kill_perm}**", "",
            "## Honest scope", "",
            "The L4 window is a SEROTONERGIC screen: scaffold plus membrane permeability. Most of "
            "the plasticity literature indexed here is about compounds it was never built to see "
            "(cholinergic, SSRI, GABAergic, ECM-degrading). Where the table above shows L4 "
            "predicting a constant, that is not a tie - it is the screen having no opinion about "
            "the drugs the field actually studies, and it bounds how much of this literature the "
            "compound-level score can ever be validated against.", "",
            "Generated by `scripts/125_window_assay_index.py`.", ""]
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(out), encoding="utf-8")
    L.info("wrote %s", REPORT)


if __name__ == "__main__":
    raise SystemExit(main())
