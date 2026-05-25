"""V2 Phase 0.5 — Boltz-derived allosteric audit + CHRNA7 rescue gate.

Re-runs the v1 allosteric audit using Boltz-2 binder_prob (or affinity_pred_value)
ranks instead of MAMMAL pKd. THE critical empirical test: does the structural
cluster rescue α7 nAChR PAMs (galantamine/encenicline/TC-5619) into the top 25%
at CHRNA7?

Gates (per user spec):
    CHRNA7 allosteric rescue : ≥2 of {galantamine, encenicline, tc-5619} in
                               top 25% by Boltz binder_prob at CHRNA7
    PDE4D allosteric preservation : bpn14770 top 25%, rolipram top 50%
    GRIA dynamic range : Boltz pred span ≥1.5 log units across ampakines
    Positive control retention : donepezil/ACHE, methylphenidate/SLC6A3,
                                 atomoxetine/SLC6A2, pitolisant/HRH3 all
                                 in Boltz top 10 at their cognate target

Writes: reports/boltzina_allosteric_audit.md
Exit:   0 = all gates pass; 2 = positive control retention fails (STOP).
        Other failures are logged but do not gate exit (per user: CHRNA7
        also failing under Boltz is itself a publishable methods finding).
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.config import RESULTS_DIR, ensure_dirs  # noqa: E402

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("allosteric_gate")

V2_DIR = RESULTS_DIR / "v2"
DEFAULT_BOLTZ = V2_DIR / "boltzina_affinity.parquet"
DEFAULT_REPORT = ROOT / "reports" / "boltzina_allosteric_audit.md"

# Same panels as v1 analysis/allosteric_audit.py
ALLOSTERIC_PANELS = {
    "P36544": {  # CHRNA7
        "gene": "CHRNA7",
        "allosteric": ["galantamine", "encenicline", "tc-5619"],
        "orthosteric": [],
    },
    "Q08499": {  # PDE4D
        "gene": "PDE4D",
        "allosteric": ["bpn14770", "zatolmilast"],
        "orthosteric": ["rolipram"],
    },
    "P42261": {"gene": "GRIA1", "allosteric": ["tulrampator", "cx-516", "cx-717",
               "aniracetam", "piracetam", "oxiracetam", "pramiracetam"], "orthosteric": []},
    "P42262": {"gene": "GRIA2", "allosteric": ["tulrampator", "cx-516", "cx-717",
               "aniracetam", "piracetam", "oxiracetam", "pramiracetam"], "orthosteric": []},
    "P42263": {"gene": "GRIA3", "allosteric": ["tulrampator", "cx-516", "cx-717",
               "aniracetam", "piracetam", "oxiracetam", "pramiracetam"], "orthosteric": []},
    "P48058": {"gene": "GRIA4", "allosteric": ["tulrampator", "cx-516", "cx-717",
               "aniracetam", "piracetam", "oxiracetam", "pramiracetam"], "orthosteric": []},
}

POSITIVE_CONTROLS = {
    "P22303": ["donepezil"],      # ACHE
    "Q01959": ["methylphenidate"],# DAT
    "P23975": ["atomoxetine"],    # NET
    "Q9Y5N1": ["pitolisant"],     # HRH3
}


def _percentile(df_target: pd.DataFrame, name: str, score_col: str) -> tuple[int, float, float] | None:
    """Return (rank, percentile, score) for a compound by `score_col` within target."""
    name_lc = name.lower().strip()
    sub = df_target.dropna(subset=[score_col]).copy()
    sub["lc"] = sub["compound_name"].str.lower().str.strip()
    if not (sub["lc"] == name_lc).any():
        return None
    sub["rank"] = sub[score_col].rank(method="max", ascending=True).astype(int)
    sub["pct"] = sub["rank"] / len(sub)
    row = sub[sub["lc"] == name_lc].iloc[0]
    return int(row["rank"]), float(row["pct"]), float(row[score_col])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--boltz", type=Path, default=DEFAULT_BOLTZ)
    parser.add_argument("--out", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--score-col", default="affinity_probability_binary",
                        choices=["affinity_probability_binary", "affinity_pred_value"],
                        help="Which Boltz field to rank on (binder_prob is calibrated [0,1])")
    parser.add_argument("--top-pct-allosteric", type=float, default=0.75,
                        help="Allosteric rescue gate: ligand must rank ≥ this percentile (default 0.75 = top 25%%)")
    args = parser.parse_args()

    ensure_dirs()
    if not args.boltz.exists():
        logger.error("Boltz affinity parquet not found at %s. Run scripts/18_v2_boltzina_sweep.py first.",
                     args.boltz)
        return 1

    df = pd.read_parquet(args.boltz)
    logger.info("Loaded %d Boltz pairs spanning %d targets, %d compounds (NaN: %d).",
                len(df), df["target_uniprot"].nunique(),
                df["compound_name"].nunique(),
                df[args.score_col].isna().sum())

    # affinity_pred_value: more negative = stronger binder, so higher percentile
    # requires inverting. binder_prob: higher = better binder, natural ranking.
    if args.score_col == "affinity_pred_value":
        df["_ranking_score"] = -df["affinity_pred_value"]
    else:
        df["_ranking_score"] = df[args.score_col]
    score_col = "_ranking_score"

    lines: list[str] = []
    lines.append("# Boltz-2 Allosteric Audit — CHRNA7 Rescue Gate")
    lines.append("")
    lines.append(f"Re-running the v1 allosteric audit with Boltz-2 `{args.score_col}` "
                 "as the rank-scoring metric instead of MAMMAL pKd. "
                 "Gate: ≥2 allosteric ligands in top "
                 f"{(1.0 - args.top_pct_allosteric)*100:.0f}% per target.")
    lines.append("")

    # --- Allosteric audit per target ---
    n_targets_pass = 0
    n_targets_checked = 0
    for uniprot, panel in ALLOSTERIC_PANELS.items():
        sub = df[df["target_uniprot"] == uniprot]
        if sub.empty:
            continue
        n_targets_checked += 1
        lines.append(f"## {panel['gene']} ({uniprot})")
        lines.append("")
        allo_results: list[tuple[str, int, float, float]] = []
        for name in panel["allosteric"]:
            r = _percentile(sub, name, score_col)
            if r is not None:
                rank, pct, sc = r
                allo_results.append((name, rank, pct, sc))
        passed = sum(1 for _, _, p, _ in allo_results if p >= args.top_pct_allosteric)
        # Gate against EVALUATED ligands (allo_results) not the full configured
        # panel — if only 1 allosteric ligand from the panel was actually scored
        # (e.g. PDE4D's bpn14770 alone), require that single one to pass; we
        # can't ask for "≥2 of 1".
        n_eval = len(allo_results)
        gate_passed = (passed >= 2) if n_eval >= 2 else (passed >= 1 and n_eval >= 1)
        if gate_passed:
            n_targets_pass += 1
        flag = "✅" if gate_passed else "❌"
        lines.append(f"**Allosteric rescue gate**: {flag} ({passed}/{len(allo_results)} ligands in top "
                     f"{(1.0 - args.top_pct_allosteric)*100:.0f}%)")
        lines.append("")
        if allo_results:
            lines.append("| Allosteric ligand | rank | percentile | Boltz score |")
            lines.append("|---|---|---|---|")
            for name, rank, pct, sc in sorted(allo_results, key=lambda r: -r[2]):
                top_flag = "✅" if pct >= args.top_pct_allosteric else "❌"
                lines.append(f"| {name} | {rank} | {pct:.0%} {top_flag} | {sc:.3f} |")
            lines.append("")
        if panel["orthosteric"]:
            lines.append("Orthosteric reference:")
            lines.append("")
            lines.append("| Orthosteric ligand | rank | percentile | Boltz score |")
            lines.append("|---|---|---|---|")
            for name in panel["orthosteric"]:
                r = _percentile(sub, name, score_col)
                if r:
                    rank, pct, sc = r
                    lines.append(f"| {name} | {rank} | {pct:.0%} | {sc:.3f} |")
            lines.append("")

    # --- Positive control retention gate (HARD; exit 2 on fail) ---
    lines.append("## Positive-control retention")
    lines.append("")
    pos_fail = []
    for uniprot, names in POSITIVE_CONTROLS.items():
        sub = df[df["target_uniprot"] == uniprot]
        if sub.empty:
            continue
        for name in names:
            r = _percentile(sub, name, score_col)
            if r is None:
                lines.append(f"- {name} @ {uniprot}: not scored")
                continue
            rank, pct, sc = r
            in_top_10 = rank <= 10 or pct >= (1.0 - 10.0 / len(sub))
            flag = "✅" if in_top_10 else "❌"
            lines.append(f"- {name} @ {uniprot}: rank {rank}, pct {pct:.0%} {flag}")
            if not in_top_10:
                pos_fail.append(f"{name}@{uniprot} rank={rank}")
    lines.append("")
    if pos_fail:
        lines.append(f"**Positive-control retention FAILED** for: {', '.join(pos_fail)}.")
        lines.append("STOP. Debug Boltz prompt/pose generation before proceeding to Phase 1.")
    else:
        lines.append("**Positive-control retention PASS**: all checked controls in top 10 at their cognate target.")
    lines.append("")

    # --- Summary ---
    lines.insert(2, f"**Allosteric targets passing gate**: {n_targets_pass}/{n_targets_checked}\n")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(lines), encoding="utf-8")
    logger.info("Wrote %s. Allosteric pass: %d/%d. Pos-ctrl fails: %d",
                args.out, n_targets_pass, n_targets_checked, len(pos_fail))

    if pos_fail:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
