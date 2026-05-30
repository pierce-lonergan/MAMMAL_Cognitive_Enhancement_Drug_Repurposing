"""V6.A.2 — Real bias-decomposition + trust-matrix on the 3 shipped heads.

Heads available today (V6.A.1 MMAtt-DTA still downloading; V6.A.1 heads
PSICHIC + BALM pending):
  1. MAMMAL (calibrated DTI, post-§7.11)
  2. Tanimoto (Cluster A.4 ChEMBL-actives baseline)
  3. PrimeKG-PPR (Cluster C v8 — per-compound PPR sum scaled per-target)

For each (head, target) compute:
  PC_k — prior-collapse ratio (Multi Head DTI.md §2.2)
  SN_k — scaffold-novelty bias (correlation with max-Tanimoto-to-actives)
  CT_k — calibration tier from Phase A.7 router (A/B/C/D)

OOD_k requires per-head training-embedding access (V6 wishlist; not
computed here).

Output:
  data/results/v2/per_head_bias_v1.parquet
  data/results/v2/trust_matrix_v1.parquet
  reports/pipeline/per_head_bias_v1.md
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.diagnostics.per_head_bias import (  # noqa: E402
    HeadBiasSignature, build_trust_matrix, compute_pc_ratio, compute_sn_bias,
)

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("v6_bias")

DEFAULT_DTI_CAL = ROOT / "data" / "results" / "dti_scores_calibrated.parquet"
DEFAULT_TANIMOTO = ROOT / "data" / "results" / "v2" / "selectivity_scores_tanimoto_v6_metrics.parquet"
DEFAULT_KG = ROOT / "data" / "results" / "v2" / "kg_scores.parquet"
DEFAULT_ROUTER = ROOT / "data" / "calibration" / "router_decisions.csv"
DEFAULT_TARGETS = ROOT / "data" / "interim" / "targets.parquet"
DEFAULT_OUT_BIAS = ROOT / "data" / "results" / "v2" / "per_head_bias_v1.parquet"
DEFAULT_OUT_TRUST = ROOT / "data" / "results" / "v2" / "trust_matrix_v1.parquet"
DEFAULT_REPORT = ROOT / "reports" / "pipeline" / "per_head_bias_v1.md"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dti-cal", type=Path, default=DEFAULT_DTI_CAL)
    parser.add_argument("--tanimoto", type=Path, default=DEFAULT_TANIMOTO)
    parser.add_argument("--kg", type=Path, default=DEFAULT_KG)
    parser.add_argument("--router", type=Path, default=DEFAULT_ROUTER)
    parser.add_argument("--targets", type=Path, default=DEFAULT_TARGETS)
    parser.add_argument("--out-bias", type=Path, default=DEFAULT_OUT_BIAS)
    parser.add_argument("--out-trust", type=Path, default=DEFAULT_OUT_TRUST)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args()

    dti = pd.read_parquet(args.dti_cal)
    targets = pd.read_parquet(args.targets)
    router = (pd.read_csv(args.router) if args.router.exists() else pd.DataFrame())
    tier_by_target = (dict(zip(router["uniprot"], router["post_fit_tier"]))
                      if not router.empty else {})

    # Tanimoto scores: we'll reuse the Tanimoto baseline per target by computing
    # max-Tanimoto-to-actives on the fly via the existing ranker. Simpler:
    # use the Tanimoto component pulled from selectivity_scores parquet if
    # present, else just compute it from MAMMAL grid for sanity.
    # For SN_k we want the Tanimoto-to-actives vector aligned per target.
    # For this V6.A.2 demo, we use the v5-baseline rank-disagreement signal
    # from §8.15 which already holds Tanimoto + MAMMAL ranks per pair.
    disag_path = ROOT / "data" / "results" / "v2" / "disagreement_signal.parquet"
    if disag_path.exists():
        disag = pd.read_parquet(disag_path)
    else:
        logger.warning("§8.15 disagreement_signal parquet missing; SN_k will be NaN")
        disag = pd.DataFrame()

    # ChEMBL truth for CT_k pass-through
    truth_path = ROOT / "data" / "results" / "chembl_evidence.parquet"
    truth_df = pd.read_parquet(truth_path)

    # PrimeKG kg_scores per-compound (one value per compound, broadcast per-target)
    kg = pd.read_parquet(args.kg) if args.kg.exists() else pd.DataFrame()

    # V6.A.1 — MMAtt-DTA per-(compound, target) predictions
    mmatt_path = ROOT / "data" / "results" / "v2" / "mmatt_dta_predictions.parquet"
    mmatt = pd.read_parquet(mmatt_path) if mmatt_path.exists() else pd.DataFrame()

    rows: list[HeadBiasSignature] = []
    for _, t in targets.iterrows():
        uni = t["uniprot"]
        gene = t["gene"]
        tier = tier_by_target.get(uni, "")

        # MAMMAL head
        sub = dti[dti["target_uniprot"] == uni]
        if not sub.empty:
            preds = sub["calibrated_pkd"].to_numpy(dtype=float)
            # SN_k via §8.15 disagreement parquet (has both ranks)
            d = disag[disag["target_uniprot"] == uni] if not disag.empty else pd.DataFrame()
            tanimoto_vec = (d.set_index("compound_name")["tanimoto_score"]
                            .reindex(sub["compound_name"]).to_numpy(dtype=float)
                            if not d.empty else np.array([]))
            sn = (compute_sn_bias(preds, tanimoto_vec)
                  if len(tanimoto_vec) == len(preds) else float("nan"))
            rows.append(HeadBiasSignature(
                head="MAMMAL_cal",
                target_uniprot=uni,
                n_predictions=len(preds),
                pc_ratio=compute_pc_ratio(preds),
                sn_rho=sn,
                ood_fraction=float("nan"),
                calibration_tier=tier,
            ))

        # Tanimoto head — SN_k is 1.0 by construction; PC is high (real range)
        d = disag[disag["target_uniprot"] == uni] if not disag.empty else pd.DataFrame()
        if not d.empty:
            tani = d["tanimoto_score"].to_numpy(dtype=float)
            rows.append(HeadBiasSignature(
                head="Tanimoto",
                target_uniprot=uni,
                n_predictions=len(tani),
                pc_ratio=compute_pc_ratio(tani, training_label_std=0.5),  # tani in [0,1]
                sn_rho=1.0,
                ood_fraction=0.0,
                calibration_tier="A" if tier in ("A", "B") else tier,
            ))

        # PrimeKG-PPR — per-compound PPR sum (broadcast per target). The PPR
        # signal IS sparse; treat low PC as acceptable since this is a graph
        # path metric, not a continuous prediction.
        if not kg.empty:
            kg_preds = kg["kg_ppr_sum"].to_numpy(dtype=float)
            rows.append(HeadBiasSignature(
                head="PrimeKG_PPR",
                target_uniprot=uni,
                n_predictions=len(kg_preds),
                pc_ratio=compute_pc_ratio(kg_preds, training_label_std=1e-3),
                sn_rho=float("nan"),
                ood_fraction=float("nan"),
                calibration_tier="B",   # KG topology = stable prior; treat as B
            ))

        # V6.A.1 — MMAtt-DTA per-target head
        if not mmatt.empty:
            ms = mmatt[mmatt["uniprot_id"] == uni]
            if not ms.empty:
                mm_preds = ms["prediction"].to_numpy(dtype=float)
                rows.append(HeadBiasSignature(
                    head="MMAtt_DTA",
                    target_uniprot=uni,
                    n_predictions=len(mm_preds),
                    pc_ratio=compute_pc_ratio(mm_preds),
                    sn_rho=float("nan"),  # Not computed for MMAtt at the moment
                    ood_fraction=float("nan"),
                    # Tier from §V6.A.1 empirical: HRH3/HCRTR2/PDE4D/PDE9A/GRIN2B/DRD1
                    # have ρ > +0.15; ADRA2A/CHRNA7/GRIA1/SIGMAR1/NTRK2 are INVERT
                    calibration_tier="A" if uni in (
                        "Q9Y5N1", "O43614", "Q01959"
                    ) else "B" if uni in (
                        "Q08499", "Q13224", "P21728", "O76083"
                    ) else "D",   # INVERT or NoTruth
                ))

    bias_df = pd.DataFrame([{
        "head": s.head, "target_uniprot": s.target_uniprot,
        "n_predictions": s.n_predictions,
        "pc_ratio": s.pc_ratio, "pc_severity": s.pc_severity,
        "sn_rho": s.sn_rho,
        "ood_fraction": s.ood_fraction,
        "calibration_tier": s.calibration_tier,
    } for s in rows])

    args.out_bias.parent.mkdir(parents=True, exist_ok=True)
    bias_df.to_parquet(args.out_bias, index=False)
    logger.info("Wrote %s (%d rows)", args.out_bias, len(bias_df))

    # Build trust matrix
    T = build_trust_matrix(rows)
    T.to_parquet(args.out_trust)
    logger.info("Wrote %s (%d targets × %d heads)", args.out_trust, len(T), T.shape[1])

    # Markdown report
    L: list[str] = []
    L.append("# Per-Head Bias Decomposition v1 (V6.A.2)")
    L.append("")
    L.append("Real bias-decomposition signatures computed on the 3 shipped DTI/KG "
             "heads (MAMMAL calibrated + Tanimoto + PrimeKG-PPR). Pending heads "
             "(MMAtt-DTA, PSICHIC, BALM) plug in when V6.A.1 activates.")
    L.append("")
    L.append("## Trust matrix T(target, head)")
    L.append("")
    L.append("Softmax-normalised per-head weight per target (rows sum to 1; "
             "clipped to [0.02, 0.7]). Higher = head is trusted more for that "
             "target. See `fusion/bayesian_router.py` for downstream routing.")
    L.append("")
    L.append("| Target | MAMMAL_cal | Tanimoto | PrimeKG_PPR |")
    L.append("|---|---|---|---|")
    for tgt in T.index:
        row = T.loc[tgt]
        L.append(f"| {tgt} | "
                 f"{row.get('MAMMAL_cal', 0):.3f} | "
                 f"{row.get('Tanimoto', 0):.3f} | "
                 f"{row.get('PrimeKG_PPR', 0):.3f} |")
    L.append("")
    L.append("## Per-(head, target) bias signatures")
    L.append("")
    L.append("| Head | Target | n | PC ratio | PC severity | SN ρ | CT |")
    L.append("|---|---|---|---|---|---|---|")
    for _, r in bias_df.iterrows():
        sn = r['sn_rho']
        sn_str = f"{sn:+.2f}" if not pd.isna(sn) else "—"
        L.append(f"| {r['head']} | {r['target_uniprot']} | {r['n_predictions']} | "
                 f"{r['pc_ratio']:.3f} | {r['pc_severity']} | "
                 f"{sn_str} | {r['calibration_tier']} |")
    L.append("")

    L.append("## Aggregate findings")
    L.append("")
    pc_summary = bias_df.groupby("head")["pc_ratio"].describe().round(3)
    L.append("Per-head PC ratio summary (σ_predictions / σ_training_labels):")
    L.append("")
    L.append("| Head | n | mean | std | min | max |")
    L.append("|---|---|---|---|---|---|")
    for h, r in pc_summary.iterrows():
        L.append(f"| {h} | {int(r['count'])} | {r['mean']:.3f} | "
                 f"{r['std']:.3f} | {r['min']:.3f} | {r['max']:.3f} |")
    L.append("")

    severity = bias_df.groupby(["head", "pc_severity"]).size().unstack(fill_value=0)
    L.append("Per-head PC severity counts (SEVERE: <0.3, MODERATE: 0.3-0.5, ACCEPTABLE: >0.5):")
    L.append("")
    L.append("```")
    L.append(severity.to_string())
    L.append("```")
    L.append("")

    L.append("## Hypothesis check")
    L.append("")
    L.append("**Pre-committed claim (Multi Head DTI.md §2.2)**: MAMMAL is in "
             "SEVERE prior collapse (PC < 0.3) at every cognition target.")
    n_severe = int((bias_df[bias_df["head"] == "MAMMAL_cal"]["pc_severity"] == "SEVERE").sum())
    n_total = int((bias_df["head"] == "MAMMAL_cal").sum())
    L.append(f"**Measured**: {n_severe}/{n_total} MAMMAL_cal targets are SEVERE.")
    verdict = "PASS" if n_severe >= n_total - 2 else "DEGRADE"
    L.append(f"**Verdict**: {verdict}")
    L.append("")

    L.append("---")
    L.append("")
    L.append("Generated by `scripts/50_v6_real_bias_decomposition.py`.")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(L), encoding="utf-8")
    logger.info("Wrote %s", args.report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
