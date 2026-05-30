"""V5 wet-lab shortlist v6 — compose calibrated+znorm v6 fusion, faceted v5,
ADMET-AI gates, and §8.0b-zn liability gates into one wet-lab handoff document.

The v4 faceted shortlist surfaced compounds by mechanism class but didn't see
the off-target liability panel (because §8.0b infrastructure was queued at the
time). This v6 deliverable adds:

  - liability_status column (PASS / FLAG / CUT from §8.0b-zn) on every row
  - final_status column (CUT > FLAG > PASS precedence over ADMET + liability)
  - a section ranking the calibrated+znorm v6 fusion top-25 directly
  - explicit wet-lab-eligible set (final_status == PASS only)
  - explicit FLAG list (manual review needed)
  - explicit CUT list (excluded with reason)

Inputs:
    data/results/v2/final_ranking_v6_calibrated_znorm.parquet  (v6 fusion)
    data/results/v2/faceted_shortlist_v5.parquet               (v5 faceted)
    data/results/v2/combined_gates.parquet                     (ADMET + liability)

Output:
    reports/wet-lab/wet_lab_shortlist_v6_full.md
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

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("v5_shortlist")

DEFAULT_V6_RANKING = ROOT / "data" / "results" / "v2" / "final_ranking_v6_calibrated_znorm.parquet"
DEFAULT_FACETED = ROOT / "data" / "results" / "v2" / "faceted_shortlist_v5.parquet"
DEFAULT_COMBINED_GATES = ROOT / "data" / "results" / "v2" / "combined_gates.parquet"
DEFAULT_REPORT = ROOT / "reports" / "wet-lab" / "wet_lab_shortlist_v6_full.md"


def _ascii(s: object) -> str:
    """ASCII-safe rendering for Windows console + markdown tables."""
    if s is None or (isinstance(s, float) and pd.isna(s)):
        return ""
    return (str(s).encode("ascii", "replace").decode("ascii")
            .replace("?", "≥"))  # heuristic — only fixes the common one


def _annotate(ranking_df: pd.DataFrame, gates_df: pd.DataFrame) -> pd.DataFrame:
    """Join ADMET+liability gate verdicts onto a ranked compound table."""
    g = gates_df[[
        "compound_name", "admet_status", "liability_status",
        "liability_note", "final_status",
    ]].copy()
    # Normalise compound names for join (lower + strip)
    g["_join"] = g["compound_name"].str.lower().str.strip()
    out = ranking_df.copy()
    out["_join"] = out["compound_name"].str.lower().str.strip()
    out = out.merge(g.drop(columns=["compound_name"]), on="_join", how="left")
    out = out.drop(columns=["_join"])
    # Fill NaN for compounds not in gates (shouldn't happen, but defensive)
    for col in ("admet_status", "liability_status", "final_status"):
        if col in out.columns:
            out[col] = out[col].fillna("UNKNOWN")
    out["liability_note"] = out["liability_note"].fillna("")
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v6-ranking", type=Path, default=DEFAULT_V6_RANKING)
    parser.add_argument("--faceted", type=Path, default=DEFAULT_FACETED)
    parser.add_argument("--combined-gates", type=Path, default=DEFAULT_COMBINED_GATES)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--top-n", type=int, default=25)
    args = parser.parse_args()

    rk = pd.read_parquet(args.v6_ranking)
    fac = pd.read_parquet(args.faceted)
    gates = pd.read_parquet(args.combined_gates)
    logger.info("Loaded: v6 ranking %d rows; faceted %d rows; combined gates %d rows.",
                len(rk), len(fac), len(gates))

    rk_ann = _annotate(rk, gates)
    fac_ann = _annotate(fac, gates)

    # --- Render markdown ----------------------------------------------------
    L: list[str] = []
    L.append("# Wet-lab Shortlist v6 — Full Composition")
    L.append("")
    L.append("**Pipeline**: V5 calibrated+znorm RRF fusion (V4 §4.4 + §4.8 wired) "
             "+ V4 faceted top-5-per-mechanism-class + §8.0b-zn within-target "
             "Z-norm liability gating.")
    L.append("")
    L.append("This document is the V5 wet-lab handoff. Every compound carries "
             "two independent gate verdicts:")
    L.append("- `admet_status` from `gates/admet_gates.py` (BBBP / DILI / hERG / "
             "P-gp / CYP / Ames / clearance + regulatory bypass for approved drugs)")
    L.append("- `liability_status` from `gates/liability_panel.py` in **z-norm "
             "mode** (§8.0b-zn): Tier 1 CUT @ z≥+2σ within target on KCNH2 / "
             "OPRM1 / HTR2B / CNR1 / HRH1 / CHRM1 / MAOA; Tier 2 FLAG @ ≥+1.5σ; "
             "Tier 3 informational @ ≥+1.0σ.")
    L.append("")
    L.append("**`final_status` precedence**: CUT > FLAG > PASS. A compound is "
             "wet-lab-eligible iff `final_status == PASS`.")
    L.append("")

    final_counts = rk_ann["final_status"].value_counts().to_dict()
    L.append("## Pipeline-wide composition (298 compounds total)")
    L.append("")
    L.append(f"- **PASS** (wet-lab eligible): **{final_counts.get('PASS', 0)}**")
    L.append(f"- **FLAG** (manual review): **{final_counts.get('FLAG', 0)}**")
    L.append(f"- **CUT** (excluded): **{final_counts.get('CUT', 0)}**")
    L.append("")

    # === Section 1: top-25 calibrated+znorm v6 fusion =======================
    L.append(f"## §1. Top {args.top_n} (calibrated MAMMAL + Z-norm + Tanimoto, "
             "ALL gate states)")
    L.append("")
    L.append("This is the raw v6 ranking with gate annotations — useful for "
             "sanity-checking against the v3/v4 reports. **Note**: FLAG/CUT "
             "rows are present but should not be acted on without resolving "
             "their liability flag.")
    L.append("")
    L.append("| # | Compound | Tier | RRF | MAMMAL best (target) | ADMET | Liability | Final |")
    L.append("|---|---|---|---|---|---|---|---|")
    for i, r in rk_ann.head(args.top_n).iterrows():
        L.append(f"| {i+1} | {_ascii(r['compound_name'])} | "
                 f"{_ascii(r.get('evidence_tier', '?'))} | "
                 f"{r['rrf_score']:.3f} | "
                 f"{r['mammal_best_pkd']:.2f} ({_ascii(r.get('mammal_best_target', '?'))}) | "
                 f"{_ascii(r['admet_status'])} | "
                 f"{_ascii(r['liability_status'])} | "
                 f"**{_ascii(r['final_status'])}** |")
    L.append("")

    # === Section 2: PASS-only top-25 =======================================
    pass_only = rk_ann[rk_ann["final_status"] == "PASS"].reset_index(drop=True)
    L.append(f"## §2. Top {args.top_n} (PASS-only — wet-lab-eligible set)")
    L.append("")
    L.append("After applying both ADMET-AI and §8.0b-zn liability gates with "
             "CUT > FLAG > PASS precedence. This is the **production wet-lab "
             "handoff**.")
    L.append("")
    L.append(f"_{len(pass_only)} compounds eligible total; showing top "
             f"{min(args.top_n, len(pass_only))}._")
    L.append("")
    L.append("| # | Compound | Tier | RRF | MAMMAL best (target) | Mech class hint |")
    L.append("|---|---|---|---|---|---|")
    for i, r in pass_only.head(args.top_n).iterrows():
        L.append(f"| {i+1} | {_ascii(r['compound_name'])} | "
                 f"{_ascii(r.get('evidence_tier', '?'))} | "
                 f"{r['rrf_score']:.3f} | "
                 f"{r['mammal_best_pkd']:.2f} ({_ascii(r.get('mammal_best_target', '?'))}) | "
                 f"— |")
    L.append("")

    # === Section 3: faceted shortlist + liability ==========================
    L.append("## §3. Faceted shortlist v5 + §8.0b-zn liability annotation")
    L.append("")
    L.append("The V4 faceted shortlist (top-5 per mechanism class + 9 "
             "targeted-pair facets) now annotated with §8.0b-zn liability. "
             "A FLAG/CUT compound in a facet means the mechanism class is "
             "useful but THIS specific compound has off-target concerns.")
    L.append("")

    if "facet_type" in fac_ann.columns:
        for facet_type in sorted(fac_ann["facet_type"].dropna().unique()):
            L.append(f"### Facet type: `{_ascii(facet_type)}`")
            L.append("")
            L.append("| Facet | Rank | Compound | Top target | Score | Gini | Liability | Final |")
            L.append("|---|---|---|---|---|---|---|---|")
            sub = fac_ann[fac_ann["facet_type"] == facet_type].sort_values(
                ["facet_name", "facet_rank"])
            for _, r in sub.iterrows():
                L.append(f"| {_ascii(r['facet_name'])} | "
                         f"{int(r['facet_rank'])} | "
                         f"{_ascii(r['compound_name'])} | "
                         f"{_ascii(r.get('top_target', '?'))} | "
                         f"{r.get('composite_score', float('nan')):.3f} | "
                         f"{r.get('gini', float('nan')):.2f} | "
                         f"{_ascii(r['liability_status'])} | "
                         f"**{_ascii(r['final_status'])}** |")
            L.append("")

    # === Section 4: FLAG compounds in v6 top-50 ============================
    flag_top = rk_ann[rk_ann["final_status"] == "FLAG"].head(20)
    L.append(f"## §4. FLAG compounds (top {len(flag_top)} by v6 RRF) — manual review")
    L.append("")
    L.append("These compounds rank well on efficacy signals but flag on "
             "ADMET or §8.0b-zn liability. Each should be triaged against the "
             "specific liability before any wet-lab spend.")
    L.append("")
    L.append("| # | Compound | RRF | MAMMAL best (target) | ADMET | Liability | Note |")
    L.append("|---|---|---|---|---|---|---|")
    for i, r in flag_top.iterrows():
        L.append(f"| {i+1} | {_ascii(r['compound_name'])} | "
                 f"{r['rrf_score']:.3f} | "
                 f"{r['mammal_best_pkd']:.2f} ({_ascii(r.get('mammal_best_target', '?'))}) | "
                 f"{_ascii(r['admet_status'])} | "
                 f"{_ascii(r['liability_status'])} | "
                 f"{_ascii(r['liability_note'])[:80]} |")
    L.append("")

    # === Section 5: CUT compounds in v6 top-100 ============================
    cut_top = rk_ann[rk_ann["final_status"] == "CUT"].head(25)
    L.append(f"## §5. CUT compounds (top {len(cut_top)} by v6 RRF) — excluded with reason")
    L.append("")
    L.append("These compounds rank well on efficacy but are CUT by ADMET, "
             "§8.0b-zn liability, or both. Listed for transparency — they "
             "would otherwise dominate the top-N if gates were absent.")
    L.append("")
    L.append("| # | Compound | RRF | ADMET | Liability | Cut reason |")
    L.append("|---|---|---|---|---|---|")
    for i, r in cut_top.iterrows():
        L.append(f"| {i+1} | {_ascii(r['compound_name'])} | "
                 f"{r['rrf_score']:.3f} | "
                 f"{_ascii(r['admet_status'])} | "
                 f"{_ascii(r['liability_status'])} | "
                 f"{_ascii(r['liability_note'])[:80]} |")
    L.append("")

    # === Section 6: methodology + caveats ===================================
    L.append("## §6. Methodology & honest caveats")
    L.append("")
    L.append("- **Calibrated MAMMAL** (V4 §7.11): isotonic per-target post-hoc "
             "calibration was applied to the raw DTI grid. 18 targets received "
             "isotonic calibrators; 4 pass through. Tier-A targets (SLC6A3, "
             "Spearman ρ_post = +0.62) are well-calibrated; Tier-D targets are "
             "down-weighted in the per-target weights file.")
    L.append("- **Z-norm within target** (V4 §4.8 / §7.18): per-target Z-norm "
             "applied before RRF to prevent cross-target scale heterogeneity "
             "(PDE9A calibrated mean 10.4 vs SLC6A3 calibrated mean 6.6 would "
             "otherwise bias panel-wide ranking).")
    L.append("- **Tanimoto-to-actives baseline** (V4 §2.1): the cluster_a_tanimoto "
             "ranker beats MAMMAL ρ at every audited target (SLC6A3: +0.90 vs "
             "-0.70; ACHE: +0.81 vs +0.24). It is weighted 1.5× MAMMAL in the "
             "global config and dominates the fusion signal at most targets.")
    L.append("- **§8.0b-zn liability gating** (V4 §8.0b, this sprint): "
             "absolute-mode gates CUT 115/115 because MAMMAL pKd on liability "
             "targets has per-target std 0.02-0.17 (prior collapse). Z-norm "
             "mode gates on within-target outlier rank (≥+2σ Tier 1 CUT, "
             "≥+1.5σ Tier 2 FLAG, ≥+1.0σ Tier 3 info). Result: 80 PASS / 21 "
             "FLAG / 14 CUT — pharmacology-consistent (hydroxyzine→hERG; "
             "aripiprazole/risperidone/lurasidone→broad polypharmacology; "
             "donepezil/methylphenidate→clean).")
    L.append("- **Faceted shortlist limitation**: facets use the v5 Tanimoto-"
             "vector selectivity (Graczyk Gini + S(10×)). Re-rendering facets "
             "under calibrated+znorm MAMMAL is V5 followup (§7.18 production "
             "swap).")
    L.append("- **Cluster C absent**: PrimeKG + TxGNN still queued (separate "
             "txgnn_env venv). When live, becomes a 5th ranker in fusion.")
    L.append("- **GRIN2A / GRIN2B**: confirmed Scenario 3 structural blindness "
             "(NMDA ATD heterodimer interface invisible to single-chain "
             "inference). Calibration ceiling ~+0.2. Methodology note flags "
             "these targets for `INVERTED_TARGET_TOP` provenance.")
    L.append("- **Roberts 2020 ceiling**: even methylphenidate hits SMD = 0.21 "
             "on a generous metric. This pipeline is not searching for a "
             "miracle drug; it is enriching a candidate set so wet-lab cycles "
             "spend money on plausibility, not chemistry lottery tickets.")
    L.append("")
    L.append("---")
    L.append("")
    def _rel(p: Path) -> str:
        try:
            return str(p.resolve().relative_to(ROOT))
        except ValueError:
            return str(p)

    L.append("Generated by `scripts/36_v5_wet_lab_shortlist.py` from:")
    L.append(f"- {_rel(args.v6_ranking)}")
    L.append(f"- {_rel(args.faceted)}")
    L.append(f"- {_rel(args.combined_gates)}")

    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(L), encoding="utf-8")
    logger.info("Wrote %s (%d lines).", args.report, len(L))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
