"""Gate — no published report may be older than the code that generates it.

    python scripts/128_report_freshness.py                 # check, exit 1 on failures
    python scripts/128_report_freshness.py --list           # every stale pair, no exit code
    python scripts/128_report_freshness.py --reproduce reports/pipeline/x_v1.md
    python scripts/128_report_freshness.py --record reports/pipeline/x_v1.md --outcome unverified --why "..."

`--reproduce` is the only route to a clearing attestation, and it earns it: it
re-runs the generator into a temporary path and diffs the output against the
committed report. `reproduced` means those bytes matched, here, today.
`regenerated` is written for you when they did not match and you committed the
new output. Neither can be typed by hand.

`--record --outcome unverified` exists for generators that cannot run in this
environment (GPU, licensed data, an external service). It does NOT clear the
report. It puts it on a list, against the code commit that outran it, so the
backlog is countable and so that the next commit to that generator makes it
fail again rather than sinking further out of sight.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.provenance.report_freshness import (  # noqa: E402
    MANIFEST, VERIFIED_OUTCOMES, check, generator_args, generator_for,
    load_manifest, scan,
)


def _write_manifest(data: dict) -> None:
    p = ROOT / MANIFEST
    p.parent.mkdir(parents=True, exist_ok=True)
    data["attestations"] = dict(sorted(data.get("attestations", {}).items()))
    p.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n",
                 encoding="utf-8")


def _stale_pairs_for(report: str):
    return [s for s in scan(ROOT)[0] if s.report == report]


def cmd_reproduce(report: str, timeout: float = 900.0) -> int:
    pairs = _stale_pairs_for(report)
    if not pairs:
        print(f"{report}: not stale; nothing to attest")
        return 0
    script = generator_for(ROOT, report)
    if script is None:
        print(f"{report}: no runnable generator in its trailer; cannot re-run")
        return 2

    tmpdir = Path(tempfile.mkdtemp(prefix="freshness-"))
    try:
        out = tmpdir / Path(report).name
        try:
            proc = subprocess.run(
                [sys.executable, script, *generator_args(ROOT, report),
                 "--report", str(out)],
                cwd=ROOT, capture_output=True, text=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            print(f"{report}: `{script}` exceeded {timeout}s")
            print("    use --record --outcome unverified --why '...' instead")
            return 2
        if proc.returncode != 0 or not out.exists():
            tail = (proc.stderr or proc.stdout or "").strip().splitlines()[-4:]
            print(f"{report}: `{script} --report` did not produce a report "
                  f"(exit {proc.returncode})")
            for line in tail:
                print(f"    {line}")
            print("    use --record --outcome unverified --why '...' instead")
            return 2

        published = (ROOT / report).read_bytes()
        fresh = out.read_bytes()
        matched = published == fresh
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    if not matched:
        print(f"{report}: RE-RUN DIFFERS from the committed report.")
        print( "    The published numbers are not what this code produces.")
        print(f"    Regenerate in place (`python {script}`), review the diff, commit it,")
        print( "    and this report clears on its own commit date -- no attestation needed.")
        return 1

    data = load_manifest(ROOT)
    atts = data.setdefault("attestations", {})
    for s in pairs:
        atts[s.key] = {
            "outcome": "reproduced",
            "code_commit": s.code_commit,
            "verified_on": date.today().isoformat(),
            "evidence": " ".join([script, *generator_args(ROOT, report),
                                  "--report <tmp>"])
                        + " reproduced the committed report byte for byte",
        }
    _write_manifest(data)
    print(f"{report}: reproduced byte for byte; recorded {len(pairs)} attestation(s)")
    return 0


def cmd_record(report: str, outcome: str, why: str) -> int:
    if outcome in VERIFIED_OUTCOMES:
        print(f"--record cannot write '{outcome}'. Use --reproduce, which checks.")
        return 2
    pairs = _stale_pairs_for(report)
    if not pairs:
        print(f"{report}: not stale; nothing to record")
        return 0
    data = load_manifest(ROOT)
    atts = data.setdefault("attestations", {})
    for s in pairs:
        atts[s.key] = {
            "outcome": outcome,
            "code_commit": s.code_commit,
            "recorded_on": date.today().isoformat(),
            "why": why,
        }
    _write_manifest(data)
    print(f"{report}: recorded {len(pairs)} pair(s) as {outcome}")
    return 0


def cmd_record_backlog(why: str) -> int:
    """Record every currently-failing NON-CITED pair as unverified, in one scan.

    This is a bulk write, which is the shape of an exemption, so it is fenced:
    it refuses cited reports outright, it can only write the outcome that
    clears nothing, and every entry it writes pins the code SHA that outran the
    report -- so the next commit to any of those generators drops the report
    straight back to failing.
    """
    findings, _ = check(ROOT)
    todo = [f for f in findings if f.fails and not f.cited]
    refused = [f for f in findings if f.fails and f.cited]
    if refused:
        print(f"refusing: {len({f.stale.report for f in refused})} CITED report(s) "
              f"are stale. Those are published outside this repo and must be "
              f"re-run, not recorded:")
        for f in refused:
            print(f"    {f.stale.report}  (outran by {f.stale.code})")
        return 2
    if not todo:
        print("no unrecorded backlog")
        return 0
    data = load_manifest(ROOT)
    atts = data.setdefault("attestations", {})
    for f in todo:
        atts[f.stale.key] = {
            "outcome": "unverified",
            "code_commit": f.stale.code_commit,
            "recorded_on": date.today().isoformat(),
            "why": why,
        }
    _write_manifest(data)
    print(f"recorded {len(todo)} pair(s) across "
          f"{len({f.stale.report for f in todo})} report(s) as unverified.")
    print("None of them is cleared. Each needs --reproduce, or a re-run and a commit.")
    return 0


def cmd_check(list_only: bool) -> int:
    findings, missing = check(ROOT)
    failing = [f for f in findings if f.fails]
    backlog = [f for f in findings if not f.fails and not f.cleared]
    cleared = [f for f in findings if f.cleared]

    reports = {f.stale.report for f in findings}
    print(f"report freshness: {len(findings)} stale pair(s) across "
          f"{len(reports)} report(s)")
    print(f"  cleared by re-run   {len(cleared)}")
    print(f"  recorded backlog    {len(backlog)}")
    print(f"  failing             {len(failing)}")
    if missing:
        print(f"  trailer names a file that does not exist   {len(missing)}")

    if list_only:
        for f in sorted(findings, key=lambda f: (not f.fails, f.stale.key)):
            state = "FAIL" if f.fails else ("ok  " if f.cleared else "todo")
            print(f"  [{state}] {f.describe()}")
        for m in missing:
            print(f"  [warn] {m}")
        return 0

    for m in missing:
        print(f"\n  WARN  {m}")
    for f in failing:
        print(f"\n  FAIL  {f.describe()}")
    if failing:
        print(f"\n{len(failing)} stale pair(s) not accounted for. Either re-run the "
              f"generator and commit,\nor run --reproduce, or record the backlog "
              f"with --record --outcome unverified --why '...'.")
        return 1
    print("\nno unaccounted staleness")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", action="store_true",
                    help="print every stale pair and its state; always exits 0")
    ap.add_argument("--reproduce", metavar="REPORT",
                    help="re-run the generator and attest only if bytes match")
    ap.add_argument("--record", metavar="REPORT",
                    help="record a non-clearing outcome for a report")
    ap.add_argument("--record-backlog", action="store_true",
                    help="record every failing non-cited pair as unverified")
    ap.add_argument("--timeout", type=float, default=900.0,
                    help="seconds a generator may run under --reproduce")
    ap.add_argument("--outcome", default="unverified")
    ap.add_argument("--why", default="")
    a = ap.parse_args()

    if a.reproduce:
        return cmd_reproduce(a.reproduce, a.timeout)
    if a.record_backlog:
        if not a.why:
            print("--record-backlog requires --why")
            return 2
        return cmd_record_backlog(a.why)
    if a.record:
        if not a.why:
            print("--record requires --why: an unexplained backlog entry is an "
                  "exemption with extra steps")
            return 2
        return cmd_record(a.record, a.outcome, a.why)
    return cmd_check(a.list)


if __name__ == "__main__":
    raise SystemExit(main())
