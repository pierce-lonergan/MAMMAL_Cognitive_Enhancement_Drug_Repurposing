"""V3 Phase C completion — wet-lab shortlist with 4-cluster scorecards.

Successor to scripts/13_wet_lab_shortlist.py. Reads the calibrated v2 fusion
output (`data/results/v2/final_ranking_calibrated.parquet`) plus the per-
cluster raw parquets, and renders a top-N markdown table where every entry
gets a 4-cluster scorecard:

    Cluster A.1 (MAMMAL DTI):   best target + pKd; polypharmacology (n targets)
    Cluster A.2 (Boltzina):     best target + binder_prob (when covered)
    Cluster B   (ADMET):        score + gate verdict + regulatory bypass
    Cluster C   (KG/TxGNN):     placeholder until kg_scores.parquet exists
    ChEMBL ground truth:        per-compound verdict counts (CORROB/AMBIG/NOVEL)
                                when chembl_evidence.parquet exists (Phase A.4)

Output: reports/wet-lab/wet_lab_shortlist_v3.md
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
    DTI_SCORES_PARQUET,
    RESULTS_DIR,
    TARGETS_PARQUET,
)

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("v3_shortlist")

V2_DIR = RESULTS_DIR / "v2"
DEFAULT_FINAL = V2_DIR / "final_ranking_calibrated.parquet"
DEFAULT_BOLTZINA = V2_DIR / "boltzina_affinity.parquet"
DEFAULT_KG = V2_DIR / "kg_scores.parquet"
DEFAULT_CHEMBL_EV = RESULTS_DIR / "chembl_evidence.parquet"
REPORT_OUT = ROOT / "reports" / "wet-lab" / "wet_lab_shortlist_v3.md"


def _top_targets_for_compound(mammal: pd.DataFrame, compound: str,
                              gene_map: dict, n: int = 3) -> list[tuple[str, str, float]]:
    """Return up to N (gene, uniprot, pKd) tuples sorted by pKd desc."""
    hits = mammal[mammal["compound_name"] == compound]
    if hits.empty:
        return []
    top = hits.sort_values("predicted_pkd", ascending=False).head(n)
    return [(gene_map.get(r["target_uniprot"], "?"),
             r["target_uniprot"],
             float(r["predicted_pkd"])) for _, r in top.iterrows()]


def _boltz_hits_for_compound(boltzina: pd.DataFrame | None, compound: str,
                             gene_map: dict, n: int = 3) -> list[tuple[str, str, float, float]]:
    """Return up to N (gene, uniprot, binder_prob, log_ic50)."""
    if boltzina is None or boltzina.empty:
        return []
    hits = boltzina[boltzina["compound_name"].str.lower() == compound.lower()]
    if hits.empty:
        return []
    prob_col = ("affinity_probability_binary"
                if "affinity_probability_binary" in hits.columns
                else "binder_prob")
    aff_col = ("affinity_pred_value"
               if "affinity_pred_value" in hits.columns
               else "log_ic50")
    top = hits.sort_values(prob_col, ascending=False).head(n)
    return [(gene_map.get(r["target_uniprot"], "?"),
             r["target_uniprot"],
             float(r[prob_col]),
             float(r[aff_col])) for _, r in top.iterrows()]


def _chembl_verdicts_for_compound(chembl_ev: pd.DataFrame | None,
                                  compound: str) -> dict[str, int]:
    """Count CORROBORATED / AMBIGUOUS / NOVEL / CONTRADICTED across all targets."""
    if chembl_ev is None or chembl_ev.empty:
        return {}
    hits = chembl_ev[chembl_ev["compound_name"].str.lower() == compound.lower()]
    if hits.empty:
        return {}
    return hits["status"].value_counts().to_dict()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--final", type=Path, default=DEFAULT_FINAL)
    parser.add_argument("--mammal", type=Path, default=DTI_SCORES_PARQUET)
    parser.add_argument("--boltzina", type=Path, default=DEFAULT_BOLTZINA)
    parser.add_argument("--kg", type=Path, default=DEFAULT_KG)
    parser.add_argument("--chembl-evidence", type=Path, default=DEFAULT_CHEMBL_EV)
    parser.add_argument("--targets", type=Path, default=TARGETS_PARQUET)
    parser.add_argument("--out", type=Path, default=REPORT_OUT)
    parser.add_argument("--top-n", type=int, default=25)
    args = parser.parse_args()

    if not args.final.exists():
        logger.error("Need %s — run scripts/15_v2_fusion.py --out-suffix _calibrated first.",
                     args.final)
        return 1

    final = pd.read_parquet(args.final)
    mammal = pd.read_parquet(args.mammal)
    targets = pd.read_parquet(args.targets)
    gene_map = dict(zip(targets["uniprot"], targets["gene"]))

    boltzina = pd.read_parquet(args.boltzina) if args.boltzina.exists() else None
    kg = pd.read_parquet(args.kg) if args.kg.exists() else None
    chembl_ev = pd.read_parquet(args.chembl_evidence) if args.chembl_evidence.exists() else None

    top = final.head(args.top_n)
    logger.info("Rendering shortlist for top %d compounds.", len(top))

    md: list[str] = []
    md.append("# Wet-Lab Shortlist v3 (4-cluster, calibrated)")
    md.append("")
    md.append("Source: `data/results/v2/final_ranking_calibrated.parquet` "
              "(produced by `scripts/15_v2_fusion.py --out-suffix _calibrated`).")
    md.append("Per-target weights from `configs/weights_calibrated.yaml` "
              "(see `reports/pipeline/calibration_report.md` for verdicts).")
    md.append("")
    md.append("Coverage at this snapshot:")
    md.append(f"  - MAMMAL DTI: {len(mammal):,} (target, compound) pairs")
    md.append(f"  - Boltzina:   {len(boltzina) if boltzina is not None else 0:,} pairs "
              f"({'overnight WSL2 sweep in progress' if boltzina is None or len(boltzina) < 500 else 'sweep complete'})")
    md.append(f"  - PrimeKG/TxGNN: {'absent' if kg is None else f'{len(kg):,} compounds'}")
    md.append(f"  - ChEMBL evidence backstop: "
              f"{'absent (Phase A.4 still running)' if chembl_ev is None else f'{len(chembl_ev):,} rows'}")
    md.append("")
    md.append("**This shortlist is a PRIORITISATION, not a wet-lab recommendation.** "
              "Read `reports/paper-drafts/methodology_v1.md` for the known failure modes "
              "(4 MAMMAL_ONLY_INVERTED targets including DAT/NET, Boltz coverage "
              "still partial). Each compound's calibrated rank reflects which "
              "targets it scores well on AFTER down-weighting WEAK / INVERTED targets.")
    md.append("")

    # Summary table first
    md.append("## Top compounds (summary)")
    md.append("")
    md.append("| Rank | Compound | RRF | ADMET | Gate | Tier | MAMMAL pKd@target |")
    md.append("|---|---|---|---|---|---|---|")
    for i, (_, r) in enumerate(top.iterrows(), 1):
        best_tgt_id = r.get("mammal_best_target", "")
        best_gene = gene_map.get(best_tgt_id, "?") if best_tgt_id else "—"
        pkd = r.get("mammal_best_pkd")
        md.append(f"| {i} | {r['compound_name']} | "
                  f"{r.get('rrf_score', float('nan')):.4f} | "
                  f"{r.get('admet_score', float('nan')):.3f} | "
                  f"{r.get('gate_status', '?')} | "
                  f"{r.get('evidence_tier', '?')} | "
                  f"{best_gene}={pkd:.2f}" + ("" if pd.isna(pkd) else "") + " |")
    md.append("")

    # Per-compound detail cards
    md.append("## Per-compound 4-cluster scorecards")
    md.append("")
    for i, (_, r) in enumerate(top.iterrows(), 1):
        name = r["compound_name"]
        md.append(f"### {i}. {name}")
        md.append("")
        md.append(f"- **Calibrated RRF**: {r.get('rrf_score', float('nan')):.4f} "
                  f"(rank {int(r.get('rrf_rank', i))})")
        md.append(f"- **Evidence tier**: `{r.get('evidence_tier', '?')}` "
                  f"  |  **Gate**: `{r.get('gate_status', '?')}`"
                  f"{'  (regulatory bypass)' if r.get('regulatory_bypass') else ''}")
        md.append("")

        # Cluster A.1 — MAMMAL
        top_tgts = _top_targets_for_compound(mammal, name, gene_map, n=3)
        md.append("**Cluster A.1 — MAMMAL DTI:**")
        if top_tgts:
            for g, u, p in top_tgts:
                md.append(f"  - {g} ({u}): pKd = {p:.2f}")
            md.append(f"  - Polypharmacology (n MAMMAL targets): "
                      f"{int(r.get('mammal_polypharm_n', 0))}")
        else:
            md.append("  - _no MAMMAL hits_")
        md.append("")

        # Cluster A.2 — Boltzina
        bz_hits = _boltz_hits_for_compound(boltzina, name, gene_map, n=3)
        md.append("**Cluster A.2 — Boltzina (structure-aware):**")
        if bz_hits:
            for g, u, prob, aff in bz_hits:
                md.append(f"  - {g} ({u}): binder_prob = {prob:.3f}, "
                          f"affinity log10(IC50 µM) = {aff:.3f}")
        else:
            md.append("  - _not yet covered by Boltz sweep_")
        md.append("")

        # Cluster B — ADMET
        admet = r.get("admet_score", float('nan'))
        md.append("**Cluster B — ADMET-AI:**")
        md.append(f"  - admet_score = {admet:.3f}  |  gate = "
                  f"`{r.get('gate_status', '?')}`")
        if r.get("gates_failed") and not pd.isna(r.get("gates_failed")):
            md.append(f"  - gates_failed: `{r.get('gates_failed')}`")
        if r.get("gates_flagged") and not pd.isna(r.get("gates_flagged")):
            md.append(f"  - gates_flagged: `{r.get('gates_flagged')}`")
        md.append("")

        # Cluster C — KG/TxGNN (placeholder until Phase B runs)
        md.append("**Cluster C — PrimeKG + TxGNN:**")
        if kg is not None:
            hit = kg[kg["compound_name"].str.lower() == name.lower()]
            if not hit.empty:
                row = hit.iloc[0]
                kg_summary: list[str] = []
                for col in ("kg_ppr_sum", "kg_shortest_path_min",
                            "kg_n_targets_reachable", "txgnn_mean_p_indication",
                            "txgnn_mean_p_contraindication"):
                    if col in row.index and not pd.isna(row[col]):
                        kg_summary.append(f"`{col}`={row[col]:.3f}")
                md.append("  - " + ", ".join(kg_summary) if kg_summary else "  - _no KG signal_")
            else:
                md.append("  - _not present in KG scores_")
        else:
            md.append("  - _Cluster C not yet run (run `scripts/23_v3_cluster_c.py` in WSL2)_")
        md.append("")

        # ChEMBL ground-truth evidence
        verdicts = _chembl_verdicts_for_compound(chembl_ev, name)
        md.append("**ChEMBL ground truth (across panel):**")
        if verdicts:
            parts = [f"{v}: {c}" for v, c in sorted(verdicts.items())]
            md.append("  - " + " | ".join(parts))
        else:
            md.append("  - _Phase A.4 backstop still running; check `data/results/chembl_evidence.parquet`_")
        md.append("")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(md), encoding="utf-8")
    logger.info("Wrote %s (%d compounds).", args.out, len(top))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
