"""V3 §8.0b — 44-target off-target liability panel orchestrator.

Three stages:
  1. (CPU/HTTP) Enrich targets_liability_seed.csv with UniProt sequences →
     data/interim/targets_liability.parquet
  2. (GPU, ~1 hr) Run MAMMAL DTI on 298 compounds × 44 liability targets →
     data/results/liability_dti.parquet  -- skipped if already exists
  3. (CPU) Apply tier-stratified gates → reports/pipeline/liability_audit_v1.md +
     data/results/v2/liability_gates.parquet + combined-with-ADMET parquet

Flags:
  --skip-dti     skip the GPU re-run if liability_dti.parquet already exists
  --dti-only     stop after stage 2 (no gating)
  --gates-only   skip stage 1+2 (assumes liability_dti.parquet exists)

Wall-clock estimate on RTX 5070: ~1 hr for the 13,112-pair DTI grid.
"""

from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.analysis.filters import filter_scores_grid  # noqa: E402
from mammal_repurposing.config import (  # noqa: E402
    COMPOUNDS_PARQUET,
    DEFAULT_BATCH_SIZE,
    INTERIM_DIR,
    RESULTS_DIR,
    ensure_dirs,
)
from mammal_repurposing.fetchers.uniprot import fetch_sequence  # noqa: E402
from mammal_repurposing.gates.liability_panel import (  # noqa: E402
    apply_liability_gates,
    combine_admet_and_liability,
    load_panel,
)

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("v3_liability")

DEFAULT_SEED = ROOT / "data" / "raw" / "targets_liability_seed.csv"
DEFAULT_INTERIM = INTERIM_DIR / "targets_liability.parquet"
DEFAULT_DTI = RESULTS_DIR / "liability_dti.parquet"
DEFAULT_GATES = RESULTS_DIR / "v2" / "liability_gates.parquet"
DEFAULT_COMBINED = RESULTS_DIR / "v2" / "combined_gates.parquet"
DEFAULT_REPORT = ROOT / "reports" / "pipeline" / "liability_audit_v1.md"


def enrich_targets(seed_csv: Path, out_parquet: Path) -> pd.DataFrame:
    """Stage 1 — fetch UniProt sequences for each liability target.

    The output parquet shape matches what scripts/04_score_dti.py expects:
    columns include uniprot, gene, sequence (plus liability metadata).
    """
    seed = pd.read_csv(seed_csv)
    logger.info("Enriching %d liability targets with UniProt sequences ...", len(seed))

    if out_parquet.exists():
        existing = pd.read_parquet(out_parquet)
        already_done = set(existing["uniprot"])
        logger.info("  Already cached: %d", len(already_done))
        seed_to_fetch = seed[~seed["uniprot_accession"].isin(already_done)]
    else:
        existing = pd.DataFrame()
        seed_to_fetch = seed

    fetched = []
    for _, row in seed_to_fetch.iterrows():
        acc = row["uniprot_accession"]
        try:
            payload = fetch_sequence(acc)
            seq = payload.get("sequence", "")
            length = payload.get("length", 0)
            logger.info("  %s (%s): %d AA", row["gene_symbol"], acc, length)
        except Exception as e:
            logger.warning("  %s (%s): UniProt fetch failed (%s) — skipping",
                           row["gene_symbol"], acc, e)
            continue
        fetched.append({
            "uniprot": acc,
            "gene": row["gene_symbol"],
            "sequence": seq,
            "length": length,
            "panel_type": "liability",
            "severity_tier": int(row["severity_tier"]),
            "liability_category": row["liability_category"],
        })
        time.sleep(0.05)    # polite throttle

    if fetched:
        new_df = pd.DataFrame(fetched)
        full = pd.concat([existing, new_df], ignore_index=True) if not existing.empty else new_df
        out_parquet.parent.mkdir(parents=True, exist_ok=True)
        full.to_parquet(out_parquet, index=False)
        logger.info("  Wrote %d enriched targets to %s", len(full), out_parquet)
        return full
    return existing


def run_mammal_dti(
    targets_parquet: Path,
    compounds_parquet: Path,
    out_parquet: Path,
    batch_size: int,
    resume: bool,
) -> None:
    """Stage 2 — MAMMAL DTI on the liability panel × library."""
    from mammal_repurposing.scoring.runner import score_grid  # noqa: PLC0415
    score_grid(
        targets_path=targets_parquet,
        compounds_path=compounds_parquet,
        out_path=out_parquet,
        batch_size=batch_size,
        resume=resume,
    )


def apply_and_render(
    liability_dti_parquet: Path,
    seed_csv: Path,
    compounds_parquet: Path,
    admet_gates_parquet: Path,
    out_gates: Path,
    out_combined: Path,
    out_report: Path,
    znorm: bool = False,
    z_cut_tier1: float = 2.0,
    z_flag_tier2: float = 1.5,
    z_flag_tier3: float = 1.0,
) -> None:
    """Stage 3 — apply gates + render markdown audit."""
    panel = load_panel(seed_csv)
    liability_dti = pd.read_parquet(liability_dti_parquet)
    logger.info("Loaded %d liability DTI predictions (%d compounds × %d targets)",
                len(liability_dti),
                liability_dti["compound_name"].nunique(),
                liability_dti["target_uniprot"].nunique())

    # Filter peptides + ADMET CUT before gating
    compounds = pd.read_parquet(compounds_parquet)
    liability_dti = filter_scores_grid(liability_dti, compounds)

    if admet_gates_parquet.exists():
        admet_gates = pd.read_parquet(admet_gates_parquet)
        cut = set(admet_gates[admet_gates["gate_status"] == "CUT"]
                  ["compound_name"].str.lower().str.strip())
        before = liability_dti["compound_name"].nunique()
        liability_dti = liability_dti[
            ~liability_dti["compound_name"].str.lower().str.strip().isin(cut)
        ]
        logger.info("After ADMET CUT filter: %d → %d compounds",
                    before, liability_dti["compound_name"].nunique())

    # Make sure target_gene is set
    if "target_gene" not in liability_dti.columns:
        seed = pd.read_csv(seed_csv)
        uni_to_gene = dict(zip(seed["uniprot_accession"], seed["gene_symbol"]))
        liability_dti = liability_dti.assign(
            target_gene=liability_dti["target_uniprot"].map(uni_to_gene),
        )

    gates = apply_liability_gates(
        liability_dti, panel,
        znorm=znorm,
        z_cut_tier1=z_cut_tier1,
        z_flag_tier2=z_flag_tier2,
        z_flag_tier3=z_flag_tier3,
    )
    out_gates.parent.mkdir(parents=True, exist_ok=True)
    gates.to_parquet(out_gates, index=False)
    logger.info("Wrote %s (%d rows).", out_gates, len(gates))

    # Combine with ADMET
    if admet_gates_parquet.exists():
        admet = pd.read_parquet(admet_gates_parquet)
        combined = combine_admet_and_liability(admet, gates)
        combined.to_parquet(out_combined, index=False)
        logger.info("Wrote %s (%d rows).", out_combined, len(combined))
    else:
        combined = gates

    # Render markdown audit
    L: list[str] = []
    L.append("# Off-Target Liability Audit v1 (§8.0b)")
    L.append("")
    L.append("Per research/4-tier/Cognition-44Target-Liability-Panel.md. "
             "Bowes-44 + Brennan-77 − peripheral irrelevancies, tier-stratified "
             "with cognition-context thresholds.")
    L.append("")
    L.append(f"Compounds evaluated: **{len(gates)}**.")
    L.append("")

    cnt = gates["liability_status"].value_counts().to_dict()
    L.append(f"**Status breakdown**: CUT={cnt.get('CUT', 0)} | "
             f"FLAG={cnt.get('FLAG', 0)} | PASS={cnt.get('PASS', 0)}")
    L.append("")

    # Tier 1 hits (the headline)
    cuts = gates[gates["liability_status"] == "CUT"].copy()
    L.append(f"## Tier 1 CUTs ({len(cuts)} compounds)")
    L.append("")
    if len(cuts):
        L.append("| Compound | Note | Top 3 liabilities |")
        L.append("|---|---|---|")
        for _, r in cuts.iterrows():
            L.append(f"| {r['compound_name']} | {r['liability_note']} | "
                     f"{r['top_3_liabilities']} |")
    else:
        L.append("_None — no compounds hit Tier 1 thresholds._")
    L.append("")

    # Tier 2 FLAGs
    flags = gates[gates["liability_status"] == "FLAG"].copy()
    L.append(f"## Tier 2 FLAGs ({len(flags)} compounds)")
    L.append("")
    if len(flags):
        L.append("| Compound | n_T2 | Note | Top 3 liabilities |")
        L.append("|---|---|---|---|")
        for _, r in flags.iterrows():
            L.append(f"| {r['compound_name']} | {int(r['n_tier_2'])} | "
                     f"{r['liability_note']} | {r['top_3_liabilities']} |")
    L.append("")

    # PASSes (just count by mechanism class)
    passes = gates[gates["liability_status"] == "PASS"].copy()
    L.append(f"## PASS ({len(passes)} compounds — clean panel)")
    L.append("")
    L.append("Compounds passed all tier-1 thresholds and <2 tier-2 thresholds. "
             "These are the v4 wet-lab-eligible set after liability gating.")
    L.append("")

    # If combined with ADMET, show the final_status distribution
    if admet_gates_parquet.exists():
        L.append("## Combined with ADMET (final_status)")
        L.append("")
        cnt_final = combined["final_status"].value_counts().to_dict()
        L.append(f"**Final breakdown**: CUT={cnt_final.get('CUT', 0)} | "
                 f"FLAG={cnt_final.get('FLAG', 0)} | PASS={cnt_final.get('PASS', 0)}")
        L.append("")
        # How many ADMET-clean compounds get newly CUT by liability?
        new_cut = combined[(combined["admet_status"] != "CUT")
                           & (combined["liability_status"] == "CUT")]
        L.append(f"**ADMET-clean but liability-CUT**: {len(new_cut)} compounds")
        L.append("")
        if len(new_cut):
            L.append("| Compound | ADMET | Liability | Note |")
            L.append("|---|---|---|---|")
            for _, r in new_cut.iterrows():
                L.append(f"| {r['compound_name']} | {r['admet_status']} | "
                         f"{r['liability_status']} | {r['liability_note']} |")
        L.append("")

    L.append("---")
    L.append("")
    L.append("Generated by `scripts/29_v3_liability_panel.py`.")
    out_report.parent.mkdir(parents=True, exist_ok=True)
    out_report.write_text("\n".join(L), encoding="utf-8")
    logger.info("Wrote %s.", out_report)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=Path, default=DEFAULT_SEED)
    parser.add_argument("--targets-out", type=Path, default=DEFAULT_INTERIM)
    parser.add_argument("--dti-out", type=Path, default=DEFAULT_DTI)
    parser.add_argument("--gates-out", type=Path, default=DEFAULT_GATES)
    parser.add_argument("--combined-out", type=Path, default=DEFAULT_COMBINED)
    parser.add_argument("--report-out", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--compounds", type=Path, default=COMPOUNDS_PARQUET)
    parser.add_argument("--admet-gates", type=Path,
                        default=RESULTS_DIR / "v2" / "admet_gates.parquet")
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--skip-dti", action="store_true",
                        help="Don't re-run MAMMAL; assume liability_dti.parquet exists.")
    parser.add_argument("--dti-only", action="store_true",
                        help="Run only stages 1-2; skip gating.")
    parser.add_argument("--gates-only", action="store_true",
                        help="Skip stages 1-2; assume parquets exist.")
    parser.add_argument("--enrich-only", action="store_true",
                        help="Stage 1 only — fetch UniProt sequences, no DTI / gating.")
    parser.add_argument("--znorm", action="store_true",
                        help="Apply per-target Z-norm to liability_dti before "
                             "gating (the §8.0b-zn fix for MAMMAL prior collapse).")
    parser.add_argument("--z-cut-tier1", type=float, default=2.0,
                        help="Z-score threshold for Tier 1 CUT verdict (default 2.0σ).")
    parser.add_argument("--z-flag-tier2", type=float, default=1.5,
                        help="Z-score threshold for Tier 2 FLAG verdict (default 1.5σ).")
    parser.add_argument("--z-flag-tier3", type=float, default=1.0,
                        help="Z-score threshold for Tier 3 informational hit (default 1.0σ).")
    args = parser.parse_args()

    ensure_dirs()

    # Stage 1
    if not args.gates_only:
        enrich_targets(args.seed, args.targets_out)
        if args.enrich_only:
            return 0

    # Stage 2
    if not args.gates_only and not args.skip_dti:
        logger.info("Stage 2: running MAMMAL DTI on liability panel "
                    "(this takes ~1 hr on RTX 5070; ensure GPU is free)")
        run_mammal_dti(
            args.targets_out, args.compounds, args.dti_out,
            args.batch_size, args.resume,
        )
        if args.dti_only:
            return 0
    elif args.skip_dti:
        logger.info("Skipping stage 2 (--skip-dti); using existing %s", args.dti_out)

    if args.dti_only:
        return 0

    # Stage 3
    if not args.dti_out.exists():
        logger.error("No %s — run with --skip-dti=false first (needs GPU).",
                     args.dti_out)
        return 1
    apply_and_render(
        args.dti_out, args.seed, args.compounds, args.admet_gates,
        args.gates_out, args.combined_out, args.report_out,
        znorm=args.znorm,
        z_cut_tier1=args.z_cut_tier1,
        z_flag_tier2=args.z_flag_tier2,
        z_flag_tier3=args.z_flag_tier3,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
