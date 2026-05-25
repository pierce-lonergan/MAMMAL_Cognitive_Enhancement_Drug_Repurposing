"""V2 Fusion — combine MAMMAL DTI + ADMET gates into RRF + provenance + narrative.

This is the v2 equivalent of v1's wet_lab_shortlist, but built on the hybrid
architecture. Currently combines Clusters A (MAMMAL DTI) and B (ADMET-AI).
Clusters Boltzina (A) and TxGNN (C) join later when their data is available;
the RRF layer absorbs them additively via the same RankerInput interface.

Inputs:
    data/results/dti_scores.parquet                (v1 MAMMAL DTI grid)
    data/interim/compounds.parquet                 (compound metadata)
    data/results/v2/admet_gates.parquet            (v2 ADMET-AI + gates)
    data/results/v2/boltzina_affinity.parquet      (optional, Cluster A v2)
    data/results/v2/kg_scores.parquet              (optional, Cluster C v2)

Outputs:
    data/results/v2/rrf_ranking.parquet
    data/results/v2/provenance.parquet
    data/results/v2/disagreement_report.md
    data/results/v2/funnel_narrative.md
    data/results/v2/final_ranking.parquet          (joined view sorted by RRF)
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

from mammal_repurposing.analysis.filters import filter_scores_grid  # noqa: E402
from mammal_repurposing.config import (  # noqa: E402
    COMPOUNDS_PARQUET,
    DTI_SCORES_PARQUET,
    RESULTS_DIR,
    ensure_dirs,
)
from mammal_repurposing.fusion.rrf import (  # noqa: E402
    RankerInput,
    rrf,
    rrf_per_target_then_compound,
)
from mammal_repurposing.provenance.disagreement_report import (  # noqa: E402
    render_markdown as render_disagreement,
)
from mammal_repurposing.provenance.narrative import render_narrative  # noqa: E402
from mammal_repurposing.provenance.tracker import build_provenance  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("v2_fusion")

V2_DIR = RESULTS_DIR / "v2"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mammal", type=Path, default=DTI_SCORES_PARQUET)
    parser.add_argument("--compounds", type=Path, default=COMPOUNDS_PARQUET)
    parser.add_argument("--admet-gates", type=Path, default=V2_DIR / "admet_gates.parquet")
    parser.add_argument("--boltzina", type=Path, default=V2_DIR / "boltzina_affinity.parquet")
    parser.add_argument("--txgnn", type=Path, default=V2_DIR / "txgnn_scores.parquet")
    parser.add_argument("--kg", type=Path, default=V2_DIR / "kg_scores.parquet")
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--k-const", type=int, default=60, help="RRF k_const (Cormack default 60)")
    parser.add_argument("--include-cut", action="store_true",
                        help="Include CUT compounds (default: drop them).")
    args = parser.parse_args()

    ensure_dirs()
    V2_DIR.mkdir(parents=True, exist_ok=True)

    # --- Load inputs --------------------------------------------------------
    mammal = pd.read_parquet(args.mammal)
    compounds = pd.read_parquet(args.compounds)
    gates = pd.read_parquet(args.admet_gates)
    logger.info("Loaded: %d MAMMAL pairs, %d compounds, %d gated compounds.",
                len(mammal), len(compounds), len(gates))

    # --- Apply v1's compound-exclusion filter (peptides etc.) --------------
    mammal = filter_scores_grid(mammal, compounds)
    logger.info("After v1 exclusion filter: %d MAMMAL pairs.", len(mammal))

    # --- Apply v2 ADMET gates ----------------------------------------------
    if not args.include_cut:
        cut_compounds = set(
            gates[gates["gate_status"] == "CUT"]["compound_name"].str.lower().str.strip()
        )
        before_n = len(mammal)
        mammal = mammal[~mammal["compound_name"].str.lower().str.strip().isin(cut_compounds)]
        logger.info("After ADMET CUT removal: %d MAMMAL pairs (dropped %d compounds).",
                    len(mammal), before_n - len(mammal))

    # --- Build RRF input from MAMMAL (per-target) + ADMET (global score) ---
    mammal_long = mammal[["target_uniprot", "compound_name", "predicted_pkd"]].copy()
    mammal_long["ranker_name"] = "cluster_a_mammal"

    # ADMET as a global per-compound ranker. Broadcast it as a per-target ranker
    # by replicating across every target. The ADMET score has no per-target
    # variation, but RRF still benefits from its constant signal because it
    # tilts the fusion toward physically clean compounds.
    admet_long_rows: list[pd.DataFrame] = []
    targets = mammal_long["target_uniprot"].unique()
    admet_by_compound = gates.set_index(
        gates["compound_name"].str.lower().str.strip()
    )["admet_score"]
    for t in targets:
        df = pd.DataFrame({
            "target_uniprot": t,
            "compound_name": admet_by_compound.index,
            "predicted_pkd": admet_by_compound.values,
            "ranker_name": "cluster_b_admet",
        })
        admet_long_rows.append(df)
    admet_long = pd.concat(admet_long_rows, ignore_index=True)

    # Optional clusters
    additional_long: list[pd.DataFrame] = []
    boltzina = None
    if args.boltzina.exists():
        boltzina = pd.read_parquet(args.boltzina)
        logger.info("Loaded Boltzina affinity: %d pairs.", len(boltzina))
        b_long = boltzina[["target_uniprot", "compound_name", "binder_prob"]].rename(
            columns={"binder_prob": "predicted_pkd"}
        )
        b_long["ranker_name"] = "cluster_a_boltzina"
        additional_long.append(b_long)

    txgnn = None
    if args.txgnn.exists():
        txgnn = pd.read_parquet(args.txgnn)
        logger.info("Loaded TxGNN scores: %d compounds.", len(txgnn))
        # TxGNN is per-compound, broadcast like ADMET
        tx_rows = []
        tx_by_compound = txgnn.set_index(
            txgnn["compound_name"].str.lower().str.strip()
        )["indication_score"]
        for t in targets:
            df = pd.DataFrame({
                "target_uniprot": t,
                "compound_name": tx_by_compound.index,
                "predicted_pkd": tx_by_compound.values,
                "ranker_name": "cluster_c_txgnn",
            })
            tx_rows.append(df)
        additional_long.append(pd.concat(tx_rows, ignore_index=True))

    long_scores = pd.concat([mammal_long, admet_long, *additional_long], ignore_index=True)
    logger.info("Fusion input: %d rows across %d rankers.",
                len(long_scores), long_scores["ranker_name"].nunique())

    # --- Run per-target then aggregate ---------------------------------------
    rrf_ranking = rrf_per_target_then_compound(
        long_scores,
        item_col="compound_name",
        target_col="target_uniprot",
        score_col="predicted_pkd",
        ranker_col="ranker_name",
        ascending=False,
        k_const=args.k_const,
    )
    rrf_ranking.to_parquet(V2_DIR / "rrf_ranking.parquet", index=False)
    logger.info("Wrote rrf_ranking.parquet (%d compounds).", len(rrf_ranking))

    # --- Build provenance ----------------------------------------------------
    prov = build_provenance(
        compounds=compounds,
        gates=gates,
        mammal_scores=mammal,
        boltzina_scores=boltzina,
        txgnn_scores=txgnn,
        rrf_ranking=rrf_ranking,
    )
    prov.to_parquet(V2_DIR / "provenance.parquet", index=False)

    # --- Disagreement diagnosis report --------------------------------------
    disagreement_md = render_disagreement(prov)
    (V2_DIR / "disagreement_report.md").write_text(disagreement_md, encoding="utf-8")
    logger.info("Wrote disagreement_report.md.")

    # --- Narrative for top-N -------------------------------------------------
    narrative_md = render_narrative(prov, top_n=args.top_n)
    (V2_DIR / "funnel_narrative.md").write_text(narrative_md, encoding="utf-8")
    logger.info("Wrote funnel_narrative.md (top %d).", args.top_n)

    # --- Final ranking parquet (joined view) --------------------------------
    final = prov.sort_values("rrf_score", ascending=False).reset_index(drop=True)
    final.to_parquet(V2_DIR / "final_ranking.parquet", index=False)
    logger.info("Wrote final_ranking.parquet (%d compounds).", len(final))

    # --- Console summary ----------------------------------------------------
    logger.info("Top 10 v2 candidates:")
    cols = ["compound_name", "evidence_tier", "rrf_score", "admet_score",
            "mammal_best_pkd", "mammal_best_target", "gate_status"]
    cols = [c for c in cols if c in final.columns]
    print(final.head(10)[cols].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
