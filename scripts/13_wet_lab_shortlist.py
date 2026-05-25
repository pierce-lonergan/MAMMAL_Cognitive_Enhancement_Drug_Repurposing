"""Phase 5 - Final wet-lab handoff shortlist.

Joins every artifact produced by the pipeline (DTI scores, aux head scores,
ChEMBL evidence, OpenTargets context, cognitive composites) and ranks the
top-N compounds by a composite priority score:

    priority = cognitive_composite * p_bbb * (1 - p_tox) * novelty_score

Where novelty_score is high if the compound is NOT already in the well-studied
healthy-cognition RCT literature (rewards under-explored chemistry) and low if
it's already a canonical nootropic.

Output:
    data/results/final_ranked.parquet  (full sorted table)
    reports/wet_lab_shortlist.md       (top-N with rich per-compound profile)
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

from mammal_repurposing.config import (  # noqa: E402
    COMPOUNDS_PARQUET,
    DTI_SCORES_PARQUET,
    RESULTS_DIR,
    ensure_dirs,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("shortlist")

DEFAULT_AUX = RESULTS_DIR / "aux_scores.parquet"
DEFAULT_CHEMBL = RESULTS_DIR / "chembl_evidence.parquet"
DEFAULT_OT = RESULTS_DIR / "opentargets_context.parquet"
DEFAULT_COMP = RESULTS_DIR / "cognitive_composites.parquet"
DEFAULT_OUT = RESULTS_DIR / "final_ranked.parquet"
DEFAULT_REPORT = ROOT / "reports" / "wet_lab_shortlist.md"

# Compounds in the canonical healthy-cognition RCT literature (lower novelty).
# Penalize these because we're trying to surface NEW chemistry, not recapitulate
# methylphenidate/modafinil for the Nth time.
KNOWN_HEALTHY_RCT = {
    "methylphenidate", "modafinil", "d-amphetamine", "atomoxetine",
    "donepezil", "rivastigmine", "galantamine", "memantine",
    "armodafinil", "lisdexamfetamine",
}


def _novelty_score(name: str, evidence_tier: str, source: str) -> float:
    """Novelty in [0, 1]. 1 = totally novel chemistry, 0 = canonical healthy-RCT drug."""
    name_lc = name.lower().strip()
    if name_lc in KNOWN_HEALTHY_RCT:
        return 0.10
    if evidence_tier == "positive_control":
        return 0.20
    if evidence_tier == "named_in_research":
        return 0.50
    if source == "chembl":
        return 0.70
    if evidence_tier == "negative_control":
        return 0.0   # explicitly exclude
    return 0.60  # default for extended_cns + anything else


def _top_targets_per_compound(scores: pd.DataFrame, top_n: int = 3) -> pd.DataFrame:
    """For each compound, return the top-N (target, pkd) hits as a flattened summary."""
    ranked = scores.sort_values(["compound_name", "predicted_pkd"], ascending=[True, False])
    grouped = ranked.groupby("compound_name").head(top_n)

    def _summarize(g):
        items = [
            f"{r.target_gene} ({r.target_uniprot}): {r.predicted_pkd:.2f}"
            for r in g.itertuples(index=False)
        ]
        return pd.Series({
            "top_targets": "; ".join(items),
            "top_target_uniprots": ";".join(r.target_uniprot for r in g.itertuples(index=False)),
            "top_pkd_max": g["predicted_pkd"].max(),
        })

    return grouped.groupby("compound_name").apply(
        _summarize, include_groups=False
    ).reset_index()


def _aggregate_chembl(ev: pd.DataFrame) -> pd.DataFrame:
    """One row per compound summarizing ChEMBL evidence across its predicted hits."""
    if ev.empty:
        return pd.DataFrame(columns=["compound_name", "chembl_label_summary",
                                     "n_corroborated", "n_novel", "n_contradicted"])
    g = ev.groupby("compound_name")
    out = pd.DataFrame({
        "n_corroborated": g.apply(lambda d: (d["label"] == "CORROBORATED").sum(), include_groups=False),
        "n_novel":        g.apply(lambda d: (d["label"] == "NOVEL").sum(), include_groups=False),
        "n_contradicted": g.apply(lambda d: (d["label"] == "CONTRADICTED").sum(), include_groups=False),
    }).reset_index()
    out["chembl_label_summary"] = out.apply(
        lambda r: f"C:{r['n_corroborated']} N:{r['n_novel']} X:{r['n_contradicted']}",
        axis=1,
    )
    return out


def _recommended_experiment(top_targets: str, n_corroborated: int) -> str:
    """Suggest next wet-lab step based on prediction profile."""
    if not top_targets:
        return "skip (no panel target)"
    if n_corroborated > 0:
        return "literature triage (already corroborated) — pursue only if NEW (compound, target) pair surfaces"
    # primary novel hit = radioligand binding at the top target
    return "Eurofins/Cerep radioligand binding assay at top predicted target (~$500-2000)"


def build_shortlist(
    *,
    scores: pd.DataFrame,
    compounds: pd.DataFrame,
    aux: pd.DataFrame | None,
    chembl: pd.DataFrame | None,
    opentargets: pd.DataFrame | None,  # noqa: ARG001 (reserved for future per-target context joins)
    composites: pd.DataFrame | None,
    top_n: int = 20,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (full_ranked, top_n_rows)."""
    top_targets = _top_targets_per_compound(scores)
    base = compounds[["name", "smiles", "evidence_tier", "source", "mechanism_class"]].rename(
        columns={"name": "compound_name"}
    )
    df = base.merge(top_targets, on="compound_name", how="left")

    if aux is not None:
        df = df.merge(
            aux[["compound_name", "p_bbb", "p_tox", "p_fda"]],
            on="compound_name", how="left",
        )
    else:
        df["p_bbb"] = None
        df["p_tox"] = None
        df["p_fda"] = None

    if chembl is not None and not chembl.empty:
        ev = _aggregate_chembl(chembl)
        df = df.merge(ev, on="compound_name", how="left")
    else:
        df["n_corroborated"] = 0
        df["n_novel"] = 0
        df["n_contradicted"] = 0
        df["chembl_label_summary"] = "no_data"

    if composites is not None:
        df = df.merge(
            composites[["compound_name", "global_composite", "working_memory",
                        "processing_speed", "learning_rate",
                        "polypharm_breadth", "polypharm_weighted_score"]],
            on="compound_name", how="left",
        )
    else:
        df["global_composite"] = None
        df["polypharm_breadth"] = 0

    df["novelty_score"] = df.apply(
        lambda r: _novelty_score(r["compound_name"], r["evidence_tier"], r["source"]),
        axis=1,
    )

    # Build the priority score. Missing aux scores default to neutral (0.5).
    p_bbb = df["p_bbb"].fillna(0.5)
    p_tox = df["p_tox"].fillna(0.5)
    composite = df["global_composite"].fillna(0.0)
    # Shift composite to a non-negative range for multiplicative weighting.
    composite_pos = composite - composite.min() + 1e-3 if composite.notna().any() else 1.0
    df["priority"] = composite_pos * p_bbb * (1.0 - p_tox) * df["novelty_score"]

    df["recommended_experiment"] = df.apply(
        lambda r: _recommended_experiment(r.get("top_targets") or "", r.get("n_corroborated") or 0),
        axis=1,
    )

    # Drop negative controls from the final ranking
    df = df[df["evidence_tier"] != "negative_control"].copy()
    df = df.sort_values("priority", ascending=False).reset_index(drop=True)
    return df, df.head(top_n)


def render_markdown(top: pd.DataFrame, full_n: int) -> str:
    lines: list[str] = []
    lines.append("# Wet-Lab Handoff Shortlist")
    lines.append("")
    lines.append(f"Top {len(top)} candidates from a pool of {full_n} compounds. "
                 "Ranked by `priority = global_composite_pos × p_BBB × (1 − p_tox) × novelty_score`.")
    lines.append("")
    lines.append("**Important caveats**:")
    lines.append("- All `predicted_pkd` values are MAMMAL DTI head outputs. "
                 "Trust rank-order, not absolute Kd — see `reports/calibration_report.md`.")
    lines.append("- These are computational predictions, NOT proven cognitive enhancers. "
                 "Wet-lab confirmation required before any further investment.")
    lines.append("- Novelty score down-weights canonical healthy-RCT drugs (methylphenidate, modafinil, etc.) "
                 "to surface less-explored chemistry. If you want the canonical list instead, sort by `global_composite`.")
    lines.append("")
    lines.append("## Ranked Shortlist")
    lines.append("")

    for idx, row in top.iterrows():
        rank = idx + 1
        lines.append(f"### {rank}. {row['compound_name']}")
        lines.append("")
        lines.append(f"- **Priority score**: {row['priority']:.4f}")
        lines.append(f"- **Mechanism class**: {row.get('mechanism_class', '?')}  |  "
                     f"**Source**: {row['source']}  |  **Evidence tier**: {row['evidence_tier']}")
        lines.append(f"- **SMILES**: `{row['smiles']}`")
        lines.append(f"- **Top predicted targets**: {row.get('top_targets', '—')}")
        comp_cell = (
            f"global {row['global_composite']:.2f} "
            f"(WM {row['working_memory']:.2f} / PS {row['processing_speed']:.2f} / LR {row['learning_rate']:.2f})"
            if pd.notna(row.get("global_composite"))
            else "—"
        )
        lines.append(f"- **Cognitive composite**: {comp_cell}")
        lines.append(f"- **Polypharm breadth** (panel targets at pKd > 6): {int(row['polypharm_breadth'])}")
        if pd.notna(row.get("p_bbb")):
            lines.append(f"- **BBB**: P(permeable) = {row['p_bbb']:.2f}  |  "
                         f"**Toxicity**: P(toxic) = {row['p_tox']:.2f}  |  "
                         f"**FDA-similarity**: P = {row['p_fda']:.2f}")
        lines.append(f"- **ChEMBL evidence**: {row.get('chembl_label_summary', 'no_data')}")
        lines.append(f"- **Novelty**: {row['novelty_score']:.2f}")
        lines.append(f"- **Recommended next experiment**: {row['recommended_experiment']}")
        lines.append("")

    lines.append("---")
    lines.append("")
    lines.append("## Scoring Formula")
    lines.append("")
    lines.append("```")
    lines.append("priority = (global_composite - min + ε) × p_BBB × (1 − p_tox) × novelty_score")
    lines.append("```")
    lines.append("")
    lines.append("Novelty score priors:")
    lines.append("- canonical healthy-RCT drug (methylphenidate, modafinil, donepezil, etc.): **0.10**")
    lines.append("- positive control in seed: 0.20")
    lines.append("- named in research deep-dive: 0.50")
    lines.append("- ChEMBL-expanded compound: 0.70")
    lines.append("- extended CNS / unclassified: 0.60")
    lines.append("- negative control: 0 (excluded)")
    lines.append("")
    lines.append("Generated by `scripts/13_wet_lab_shortlist.py`.")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", type=Path, default=DTI_SCORES_PARQUET)
    parser.add_argument("--compounds", type=Path, default=COMPOUNDS_PARQUET)
    parser.add_argument("--aux", type=Path, default=DEFAULT_AUX)
    parser.add_argument("--chembl", type=Path, default=DEFAULT_CHEMBL)
    parser.add_argument("--opentargets", type=Path, default=DEFAULT_OT)
    parser.add_argument("--composites", type=Path, default=DEFAULT_COMP)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--top-n", type=int, default=20)
    args = parser.parse_args()

    ensure_dirs()
    if not args.scores.exists() or not args.compounds.exists():
        logger.error("Missing required inputs (scores or compounds). Run earlier stages first.")
        return 1

    scores = pd.read_parquet(args.scores)
    compounds = pd.read_parquet(args.compounds)
    aux = pd.read_parquet(args.aux) if args.aux.exists() else None
    chembl = pd.read_parquet(args.chembl) if args.chembl.exists() else None
    ot = pd.read_parquet(args.opentargets) if args.opentargets.exists() else None
    comp = pd.read_parquet(args.composites) if args.composites.exists() else None

    logger.info(
        "Joining: scores=%d, compounds=%d, aux=%s, chembl=%s, opentargets=%s, composites=%s",
        len(scores), len(compounds),
        len(aux) if aux is not None else "missing",
        len(chembl) if chembl is not None else "missing",
        len(ot) if ot is not None else "missing",
        len(comp) if comp is not None else "missing",
    )

    full, top = build_shortlist(
        scores=scores, compounds=compounds, aux=aux,
        chembl=chembl, opentargets=ot, composites=comp,
        top_n=args.top_n,
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    full.to_parquet(args.out, index=False)
    logger.info("Wrote full ranked table (%d rows) to %s.", len(full), args.out)

    md = render_markdown(top, len(full))
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(md, encoding="utf-8")
    logger.info("Wrote shortlist report to %s.", args.report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
