"""Retro-validation: walk the graveyard and ask which deaths the repository has since outgrown.

Every revivable entry in `data/raw/stepping_stone_archive.csv` names a machine-checkable keystone
predicate: the specific condition whose absence killed it. This script evaluates all of them against
the repository AS IT IS NOW and reports which have become satisfied.

WHY IT EXISTS, CONCRETELY. On 2026-09-15 the B2 harness scored the L4 plasticity window on
parent-compound SMILES and recorded psilocybin as a false negative. PERSEUS had been resolving
psilocybin to psilocin through `PRODRUG_TO_ACTIVE` for months; the harness simply did not consult
it. The error was caught by hand, and the measured agreement moved 0.33 -> 0.50. Nothing in the
repository would have caught it on its own, because nothing recorded that the psilocybin verdict
DEPENDED ON prodrug resolution. This script is that missing check, generalised: it runs
`prodrug_resolution_covers(psilocybin)` and gets True, which is the signal that a verdict resting
on that dependency is stale.

THREE OUTCOMES, AND THE THIRD MATTERS MOST.

  SATISFIED     the keystone is now present. The hypothesis is flagged for human re-adjudication.
                It is NOT revived. A revival is a scientific claim and this script does not make
                claims; it makes appointments.
  UNSATISFIED   still missing, with the detail of what was actually looked at. "6 paired studies at
                n >= 25 vs required 15" is an agenda item; "still dead" is not.
  UNRESOLVABLE  the predicate names a measurement nothing in the repository can currently produce.
                This is reported LOUDLY and never silently folded into UNSATISFIED, because a
                keystone the engine cannot evaluate would otherwise read as "checked, still dead"
                forever, which is the exact failure the archive was built to end.

Writes reports/pipeline/retro_validation_v1.md. CPU only, no network.
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("retro_validation")
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.archive.predicates import (  # noqa: E402
    UnknownPredicate, evaluate, registered,
)
from mammal_repurposing.archive.stepping_stone import (  # noqa: E402
    FAILURE_MODES, assert_archive_valid,
)

ARCHIVE = ROOT / "data" / "raw" / "stepping_stone_archive.csv"
REPORT = ROOT / "reports" / "pipeline" / "retro_validation_v1.md"

#: The worked example. `prodrug_resolution_covers(psilocybin)` is True today, which is exactly the
#: dependency whose staleness was caught by hand on 2026-09-15. Running it here proves the engine
#: detects the class of error rather than merely asserting that it would.
SELF_TEST = [
    ("prodrug_resolution_covers(psilocybin)", True,
     "the dependency that was stale on 2026-09-15 and was corrected by hand"),
    ("prodrug_resolution_covers(ibogaine)", False,
     "the same dependency, still missing, and still producing a live false positive"),
]


def sweep(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in df.iterrows():
        mode = str(r["failure_mode"])
        if not FAILURE_MODES[mode].revivable:
            rows.append(dict(hypothesis_id=r["hypothesis_id"], claim=r["claim"],
                             failure_mode=mode, outcome="CLOSED",
                             detail="permanently closed: the evidence positively excludes the "
                                    "effect, so no sweep re-opens it", predicate=""))
            continue
        pred = str(r["keystone_predicate"]).strip()
        try:
            ok, detail = evaluate(pred)
            outcome = "SATISFIED" if ok else "UNSATISFIED"
        except UnknownPredicate as e:
            outcome, detail = "UNRESOLVABLE", str(e)
        rows.append(dict(hypothesis_id=r["hypothesis_id"], claim=r["claim"], failure_mode=mode,
                         outcome=outcome, detail=detail, predicate=pred))
    return pd.DataFrame(rows)


def main() -> int:
    if not ARCHIVE.exists():
        raise SystemExit(f"missing {ARCHIVE}; run scripts/129_stepping_stone_archive.py first")
    df = pd.read_csv(ARCHIVE)
    assert_archive_valid(df)          # fail closed: never sweep an archive that breaks its contract

    selftest = []
    for expr, expected, why in SELF_TEST:
        got, detail = evaluate(expr)
        selftest.append(dict(predicate=expr, expected=expected, got=got,
                             passed=(got == expected), why=why, detail=detail))
    for s in selftest:
        (L.info if s["passed"] else L.error)(
            "self-test %s: %s -> %s (expected %s)",
            "PASS" if s["passed"] else "FAIL", s["predicate"], s["got"], s["expected"])

    res = sweep(df)
    counts = dict(res["outcome"].value_counts())
    L.info("archive %d entries | %s", len(df), counts)
    for _, r in res[res.outcome == "SATISFIED"].iterrows():
        L.warning("REVIVABLE NOW: %s -- %s", r["hypothesis_id"], r["detail"])
    for _, r in res[res.outcome == "UNRESOLVABLE"].iterrows():
        L.error("UNRESOLVABLE keystone: %s -- %s", r["hypothesis_id"], r["detail"][:160])

    write_report(df, res, counts, selftest)
    # An unresolvable keystone is a defect in the archive, not a neutral state: it silently
    # converts to "still dead" on every future sweep. Exit non-zero so CI can see it.
    return 1 if (res.outcome == "UNRESOLVABLE").any() else 0


def write_report(df, res, counts, selftest) -> None:
    sat = res[res.outcome == "SATISFIED"]
    unsat = res[res.outcome == "UNSATISFIED"]
    unres = res[res.outcome == "UNRESOLVABLE"]
    closed = res[res.outcome == "CLOSED"]

    out = [
        "# Retro-validation sweep", "",
        f"**{len(sat)} of {len(df)} archived hypotheses have a keystone that is now satisfied.**", "",
        "Every revivable entry in the stepping-stone archive names the specific condition whose "
        "absence killed it. This sweep evaluates those conditions against the repository as it "
        "stands and reports which have been met. It revives nothing: a revival is a scientific "
        "claim, and this script makes appointments rather than claims.", "",
        "## Does the mechanism actually work?", "",
        "The engine exists because of a specific error. On 2026-09-15 the B2 harness scored the L4 "
        "window on parent-compound SMILES and recorded psilocybin as a false negative, while "
        "PERSEUS had been resolving psilocybin to psilocin for months. The harness did not consult "
        "the map. That was caught by hand and moved the measured agreement from 0.33 to 0.50. "
        "These two probes are that error, and its still-live twin, stated as predicates:", "",
        "| predicate | expected | got | pass | what it is |", "| --- | :---: | :---: | :---: | --- |",
    ]
    for s in selftest:
        out.append(f"| `{s['predicate']}` | {s['expected']} | {s['got']} | "
                   f"{'yes' if s['passed'] else '**NO**'} | {s['why']} |")
    out += ["", "The first returns True: a verdict resting on that dependency is stale, and the "
            "engine says so without anyone remembering to look. The second returns False, and it "
            "is not historical. Ibogaine is still absent from the map, and L4 still calls a "
            "measured null positive because of it.", ""]

    out += ["## Sweep", "",
            f"- satisfied (flag for re-adjudication): **{len(sat)}**",
            f"- unsatisfied (still dead, with the gap named): **{len(unsat)}**",
            f"- unresolvable (the engine cannot check this): **{len(unres)}**",
            f"- permanently closed (never swept): **{len(closed)}**", ""]

    if len(sat):
        out += ["### Satisfied now", "", "| hypothesis | mode | keystone | what the engine saw |",
                "| --- | --- | --- | --- |"]
        for _, r in sat.iterrows():
            out.append(f"| {r['hypothesis_id']} | {r['failure_mode']} | `{r['predicate']}` | "
                       f"{r['detail']} |")
        out.append("")

    if len(unres):
        out += ["### Unresolvable keystones (defects in the archive)", "",
                "A keystone the engine cannot evaluate reports as 'still dead' on every future "
                "sweep, which is indistinguishable from a verdict nobody ever revisits. These are "
                "failures of the archive, not of the hypotheses, and the sweep exits non-zero "
                "while any remain.", "",
                "| hypothesis | keystone | why it cannot be checked |", "| --- | --- | --- |"]
        for _, r in unres.iterrows():
            out.append(f"| {r['hypothesis_id']} | `{r['predicate']}` | {r['detail'][:300]} |")
        out.append("")

    out += ["### Still dead, and what is missing", "",
            "| hypothesis | mode | keystone | the gap, measured |", "| --- | --- | --- | --- |"]
    for _, r in unsat.iterrows():
        out.append(f"| {r['hypothesis_id']} | {r['failure_mode']} | `{r['predicate']}` | "
                   f"{r['detail']} |")

    modes = dict(df["failure_mode"].value_counts())
    out += ["", "## What the graveyard is made of", "",
            "The discrimination decision is the whole value of the archive: a kill that nature "
            "delivered is closed, and a kill that our own instrument or sample size delivered is an "
            "open question wearing a verdict's clothes.", "",
            "| failure mode | entries | revivable |", "| --- | ---: | :---: |"]
    for m, n in sorted(modes.items(), key=lambda kv: -kv[1]):
        out.append(f"| {m} | {n} | {'yes' if FAILURE_MODES[m].revivable else 'no'} |")

    n_rev = int(df["failure_mode"].map(lambda m: FAILURE_MODES[m].revivable).sum())
    out += ["", f"**{n_rev} of {len(df)} entries are revivable.** That is not a claim that "
            f"{n_rev} hypotheses are alive. It is a claim that for {n_rev} of them the recorded "
            "evidence cannot distinguish 'no effect' from 'not measurable by the test that was "
            "run', which is a different statement and a much weaker one than the ledgers currently "
            "make.", "",
            "The clearest case is in the primary healthy-adult ledger. Fourteen compounds carry "
            "`enhances_healthy_young = 0`. Seven have an interval whose upper bound excludes a "
            "target-sized effect and are properly closed. Three have intervals that admit one "
            "(dextroamphetamine reaches +0.47, against a target of 0.25). Four carry no interval "
            "at all. The ledger asserts fourteen refutations and has evidence for seven.", "",
            "That finding is not new here. `healthy_adult_robustness_v1.md` section R3 said it "
            "plainly and the ledger did not change, because a finding written into a report changes "
            "no data structure and therefore changes nothing downstream. The archive is the data "
            "structure.", "",
            f"Registered predicates: {', '.join('`' + p + '`' for p in registered())}.", "",
            "Generated by `scripts/130_retro_validation.py`.", ""]

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(out), encoding="utf-8")
    L.info("wrote %s", REPORT.relative_to(ROOT))


if __name__ == "__main__":
    raise SystemExit(main())
