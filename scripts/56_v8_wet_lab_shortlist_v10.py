"""V8.6 — Wet-lab shortlist v10 (three-factor joint composition).

Composes the full V6.A × V6.B × V7 × V8 stack into the v10 wet-lab handoff.

Inputs (any can be missing; the composer collapses to neutral defaults):
  - V6.A: data/results/v2/mmatt_for_fusion.parquet (or any 4-head pchembl
    output with [compound_name, target_uniprot, predicted_pkd])
  - V6.B: data/results/v2/cluster_d_posterior_v1.parquet (from scripts/55)
  - V7:   constructed inline from PRISMA priors via fit_effect_size_stub
  - V8:   constructed inline from a smoke MOFA+ + heuristic 5-MoA centroids

Outputs:
  data/results/v2/wet_lab_shortlist_v10.parquet
  reports/wet_lab_shortlist_v10.md

Honest caveat: this is the **architectural composition** demonstration.
Real-data v10 awaits all 4 posteriors flowing (V6.A.4 Venn-ABERS shipped,
V6.B.3 PyMC NUTS scaffold shipped, V7.3 effect-size scaffold shipped,
V8.3 MOFA+ scaffold shipped). The composer always produces SOMETHING
usable; missing axes are visibly flagged in the `evidence_axes` column.
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

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("v8_wet_lab_shortlist_v10")


def construct_v7_posterior_from_v6a(v6a_pchembl: pd.DataFrame,
                                     v6b_theta: pd.DataFrame,
                                     panel_class_map: dict[str, str]):
    """Build a V7 posterior via PRISMA stub from V6.A pchembl + V6.B θ̄."""
    from mammal_repurposing.translation.effect_size_model import (
        EffectSizeObservation, fit_effect_size_stub,
    )
    # Cluster D θ̄ per target
    theta_by_target = {row["target_uniprot"]: float(row["theta_mean"])
                       for _, row in v6b_theta.iterrows()
                       if "theta_mean" in row.index and np.isfinite(row["theta_mean"])}
    theta_sd_by_target = {row["target_uniprot"]:
                          float((row.get("theta_97p5", 0.5) -
                                 row.get("theta_2p5", -0.5)) / 4.0)
                          for _, row in v6b_theta.iterrows()
                          if "theta_mean" in row.index}

    # Per-compound (mean pchembl across V6.A 4 heads) + target + class
    obs: list[EffectSizeObservation] = []
    for c, g in v6a_pchembl.groupby("compound_name"):
        target_u = str(g["target_uniprot"].iloc[0])
        pchembl_col = next((cc for cc in ("predicted_pkd", "pchembl_mean")
                            if cc in g.columns), None)
        if pchembl_col is None:
            continue
        pchembl_m = float(g[pchembl_col].mean())
        pchembl_sd = float(g[pchembl_col].std()) if len(g) > 1 else 0.3
        relevance_m = theta_by_target.get(target_u, 0.5)
        relevance_sd = theta_sd_by_target.get(target_u, 0.15)
        class_name = panel_class_map.get(target_u, "AChE-I")    # default class
        obs.append(EffectSizeObservation(
            compound=str(c),
            class_name=class_name,
            target_uniprot=target_u,
            pchembl_post_mean=pchembl_m,
            pchembl_post_sd=pchembl_sd,
            relevance_post_mean=float(1.0 / (1.0 + np.exp(-relevance_m))),
            relevance_post_sd=relevance_sd,
            pbpk_auc_brain=1.0,
            moderators=(0, 0, 0, 0, 0),
        ))
    return fit_effect_size_stub(obs)


def construct_v8_posterior_from_pchembl(v6a_pchembl: pd.DataFrame,
                                         v6b_theta: pd.DataFrame,
                                         v7_post,
                                         panel_compound_class_map: dict[str, str]):
    """Build a V8 JointPosterior using V6.A pchembl + V6.B θ̄ + V7 g_mean +
    a heuristic phenotype-cosine = sigmoid(pchembl - 6) * sigmoid(theta)."""
    from mammal_repurposing.cluster_e.joint_phenotype import (
        compute_joint_posterior,
    )
    theta_by_target = {row["target_uniprot"]: (float(row["theta_mean"]),
                                                 float((row.get("theta_97p5", 0.5) -
                                                        row.get("theta_2p5", -0.5)) / 4.0))
                       for _, row in v6b_theta.iterrows()
                       if "theta_mean" in row.index}

    pchembl_by_compound: dict[str, tuple[float, float]] = {}
    compound_to_target: dict[str, str] = {}
    phen_by_compound: dict[str, float] = {}
    phen_centroid_by_compound: dict[str, str] = {}
    for c, g in v6a_pchembl.groupby("compound_name"):
        cstr = str(c)
        target_u = str(g["target_uniprot"].iloc[0])
        pchembl_col = next((cc for cc in ("predicted_pkd", "pchembl_mean")
                            if cc in g.columns), None)
        if pchembl_col is None:
            continue
        pchembl_m = float(g[pchembl_col].mean())
        pchembl_sd = float(g[pchembl_col].std()) if len(g) > 1 else 0.3
        pchembl_by_compound[cstr] = (pchembl_m, pchembl_sd)
        compound_to_target[cstr] = target_u
        # Heuristic phenotype cosine: high when pchembl + theta both high
        theta_m, _ = theta_by_target.get(target_u, (0.0, 0.5))
        phen = float(np.tanh((pchembl_m - 6.0) / 2.0)
                     * (1.0 / (1.0 + np.exp(-theta_m))))
        phen_by_compound[cstr] = phen
        cls = panel_compound_class_map.get(target_u, "cholinergic")
        phen_centroid_by_compound[cstr] = cls

    v7_g_post: dict[str, tuple[float, float]] = {}
    if hasattr(v7_post, "compounds"):
        for c in v7_post.compounds:
            g_m = v7_post.g_mean.get(c, 0.0)
            g_sd = float((v7_post.g_97p5.get(c, 0.0)
                          - v7_post.g_2p5.get(c, 0.0)) / 4.0)
            v7_g_post[c] = (g_m, g_sd)

    return compute_joint_posterior(
        v6a_pchembl_post=pchembl_by_compound,
        v6b_theta_post=theta_by_target,
        v7_g_post=v7_g_post,
        v8_phen_cosine=phen_by_compound,
        compound_to_target=compound_to_target,
        v8_phen_centroid=phen_centroid_by_compound,
    )


# Per V4 §13.Y + PRISMA priors — coarse mapping from UniProt → PRISMA class
PANEL_TARGET_CLASS_MAP: dict[str, str] = {
    "P22303": "AChE-I",            # ACHE
    "Q01959": "NDRI",              # SLC6A3 DAT
    "P23975": "NRI",               # SLC6A2 NET
    "Q9Y5N1": "wake_promoting",    # HRH3
    "O43614": "wake_promoting",    # HCRTR2
    "O43613": "wake_promoting",    # HCRTR1
    "P21728": "wake_promoting",    # DRD1
    "P08913": "alpha2A_agonist",   # ADRA2A
    "Q13224": "NMDA_antagonist",   # GRIN2B
    "Q12879": "NMDA_antagonist",   # GRIN2A
    "P36544": "AChE-I",            # CHRNA7 — closest non-cholinergic class
    "Q08499": "AMPA_pos_mod",      # PDE4D
    "O76083": "AMPA_pos_mod",      # PDE9A
    "Q99720": "multimodal_5HT",    # SIGMAR1
    "Q16620": "creatine",          # NTRK2
    "P42261": "AMPA_pos_mod",      # GRIA1
    "P42262": "AMPA_pos_mod",      # GRIA2
    "P42263": "AMPA_pos_mod",      # GRIA3
    "P48058": "AMPA_pos_mod",      # GRIA4
    "O43526": "alpha2A_agonist",   # KCNQ2 — channel modulator placeholder
    "O43525": "alpha2A_agonist",   # KCNQ3
    "O60741": "alpha2A_agonist",   # HCN1
}

# 5-MoA centroid mapping per V8 plan
PANEL_TARGET_CENTROID_MAP: dict[str, str] = {
    "P22303": "cholinergic",
    "P36544": "cholinergic",
    "Q01959": "catecholaminergic",
    "P23975": "catecholaminergic",
    "P21728": "catecholaminergic",
    "P08913": "catecholaminergic",
    "Q9Y5N1": "catecholaminergic",
    "O43614": "catecholaminergic",
    "O43613": "catecholaminergic",
    "Q13224": "glutamatergic",
    "Q12879": "glutamatergic",
    "P42261": "glutamatergic",
    "P42262": "glutamatergic",
    "P42263": "glutamatergic",
    "P48058": "glutamatergic",
    "Q08499": "glutamatergic",
    "O76083": "glutamatergic",
    "Q16620": "trophic_ISR",
    "Q99720": "trophic_ISR",
    "O43526": "remyelination",   # KCNQ2/3 are not remyelination but placeholder
    "O43525": "remyelination",
    "O60741": "remyelination",
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v6a", type=Path,
                        default=ROOT / "data" / "results" / "v2"
                        / "mmatt_for_fusion.parquet",
                        help="V6.A 4-head pchembl input")
    parser.add_argument("--v6b", type=Path,
                        default=ROOT / "data" / "results" / "v2"
                        / "cluster_d_posterior_v1.parquet",
                        help="V6.B Cluster D θ̄ posterior")
    parser.add_argument("--admet", type=Path,
                        default=ROOT / "data" / "results" / "v2"
                        / "admet_gates.parquet",
                        help="V4 ADMET PASS/FLAG/CUT")
    parser.add_argument("--out", type=Path,
                        default=ROOT / "data" / "results" / "v2"
                        / "wet_lab_shortlist_v10.parquet")
    parser.add_argument("--report", type=Path,
                        default=ROOT / "reports" / "wet_lab_shortlist_v10.md")
    parser.add_argument("--top-n", type=int, default=50)
    parser.add_argument("--no-roberts", action="store_true",
                        help="Skip the Roberts 2020 SMD ceiling pre-filter")
    args = parser.parse_args()

    from mammal_repurposing.fusion.joint_composition import (
        JointCompositionConfig, compose_wet_lab_shortlist_v10,
        render_v10_markdown_report,
    )

    # Load inputs
    v6a = pd.read_parquet(args.v6a) if args.v6a.exists() else None
    v6b = pd.read_parquet(args.v6b) if args.v6b.exists() else None
    admet = pd.read_parquet(args.admet) if args.admet.exists() else None

    if v6a is None:
        logger.error("V6.A pchembl input missing at %s", args.v6a)
        return 2
    logger.info("V6.A: %d (compound, target) rows; %d unique compounds",
                len(v6a), v6a["compound_name"].nunique())
    if v6b is None:
        logger.warning("V6.B Cluster D posterior missing at %s; using stubs", args.v6b)
        v6b = pd.DataFrame(columns=["target_uniprot", "gene", "theta_mean",
                                     "theta_2p5", "theta_97p5", "w_pipeline"])
    else:
        logger.info("V6.B: %d targets with posterior θ̄", len(v6b))
    if admet is None:
        logger.info("ADMET gates not loaded (optional)")

    # Construct V7 from V6.A + V6.B via PRISMA stub
    v7_post = construct_v7_posterior_from_v6a(v6a, v6b, PANEL_TARGET_CLASS_MAP)
    logger.info("V7 stub posterior: %d compounds (PRISMA class-mean × "
                "Cluster D gate × moderator debit)", len(v7_post.compounds))

    # Construct V8 joint posterior
    v8_post = construct_v8_posterior_from_pchembl(
        v6a, v6b, v7_post, PANEL_TARGET_CENTROID_MAP,
    )
    logger.info("V8 joint posterior: %d entries with 8-cell tags",
                len(v8_post.entries))

    # Compose
    cfg = JointCompositionConfig(
        enforce_roberts_ceiling=not args.no_roberts,
        top_n=args.top_n,
    )
    df = compose_wet_lab_shortlist_v10(
        v6a_pchembl=v6a,
        v6b_theta=v6b,
        v7_posterior=v7_post,
        v8_posterior=v8_post,
        admet_gates=admet,
        cfg=cfg,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(args.out, index=False)
    logger.info("Wrote %s (%d rows, top_n=%d)", args.out, len(df), cfg.top_n)

    # Render report
    render_v10_markdown_report(df, args.report, cfg)
    logger.info("Wrote %s", args.report)

    # Headline
    if not df.empty:
        n_ceiling_ok = int(df["roberts_ceiling_ok"].sum())
        n_novel = int((df["eight_cell_tag"] == "phenotype_only.novel_mechanism").sum())
        n_agreement = int((df["eight_cell_tag"] == "agreement.all_high").sum())
        logger.info("Headline: %d rows; %d ceiling PASS; %d novel-mechanism cell "
                    "(L, L, H); %d full-agreement (H, H, H)",
                    len(df), n_ceiling_ok, n_novel, n_agreement)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
