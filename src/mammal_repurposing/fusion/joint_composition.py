"""V8.6 — Joint-composition wet-lab shortlist v10.

Composes V6.A 4-head pchembl + V6.B Cluster D θ̄ + V7 effect-size + V8
phenotype into the v10 wet-lab handoff DataFrame, ranked by predicted
Hedges' g with 95% CrI, pre-filtered by Roberts 2020 SMD ceiling,
annotated with 4-axis disagreement + 8-cell classification + I_novel.

Inputs are intentionally loose (any of the 4 axes can be missing or stub
mode) so the composer always produces SOMETHING usable. Missing axes
collapse to neutral defaults; the per-compound output flags which axes
contributed.

API:
    df = compose_wet_lab_shortlist_v10(
        v6a_parquet=...,
        v6b_parquet=...,
        v7_posterior=...,    # optional EffectSizePosterior
        v8_posterior=...,    # optional JointPosterior
        admet_gates=...,     # optional V4 PASS/FLAG/CUT
        enforce_roberts_ceiling=True,
        top_n=50,
    )

Output columns:
    rank, compound, target_uniprot, target_gene,
    pchembl_mean, pchembl_sd, theta_mean, theta_sd,
    g_predicted, g_90_upper, phen_cosine, phen_centroid,
    three_way_jsd, i_novel_score, eight_cell_tag,
    admet_status, roberts_ceiling_ok, wet_lab_priority, evidence_axes
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class JointCompositionConfig:
    """Configuration for the v10 wet-lab composer."""
    roberts_ceiling_g: float = 0.50
    enforce_roberts_ceiling: bool = True
    top_n: int = 50
    weight_g: float = 1.0
    weight_jsd: float = 0.25
    weight_novel: float = 0.60
    weight_ci_penalty: float = 0.30
    cluster_d_floor: float = 0.10


def _safe_first(s: pd.Series, default=float("nan")):
    try:
        return s.iloc[0]
    except Exception:
        return default


def compose_wet_lab_shortlist_v10(
    v6a_pchembl: pd.DataFrame | None = None,
    v6b_theta: pd.DataFrame | None = None,
    v7_posterior=None,
    v8_posterior=None,
    admet_gates: pd.DataFrame | None = None,
    cfg: JointCompositionConfig | None = None,
    compound_to_target: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Compose the v10 wet-lab shortlist DataFrame.

    Args:
        v6a_pchembl: long-format DataFrame [compound_name, target_uniprot,
            predicted_pkd, predicted_pkd_sd] from V6.A 4-head Bayesian router
        v6b_theta: long-format DataFrame [target_uniprot, gene, theta_mean,
            theta_2p5, theta_97p5, w_pipeline] from V6.B.3 NUTS posterior
        v7_posterior: optional EffectSizePosterior from V7.3
        v8_posterior: optional JointPosterior from V8.5
        admet_gates: optional [compound_name, status] {PASS|FLAG|CUT}
        cfg: JointCompositionConfig
        compound_to_target: optional {compound: target_uniprot} override
    """
    cfg = cfg or JointCompositionConfig()

    # Collect all compounds across axes
    compounds: set[str] = set()
    if v6a_pchembl is not None and "compound_name" in v6a_pchembl.columns:
        compounds.update(v6a_pchembl["compound_name"].astype(str).tolist())
    if v7_posterior is not None and hasattr(v7_posterior, "compounds"):
        compounds.update(v7_posterior.compounds)
    if v8_posterior is not None and hasattr(v8_posterior, "entries"):
        compounds.update(e.compound for e in v8_posterior.entries)
    if compound_to_target:
        compounds.update(compound_to_target.keys())

    if not compounds:
        logger.warning("compose_wet_lab_shortlist_v10: no compounds in any input")
        return pd.DataFrame(columns=[
            "rank", "compound", "target_uniprot", "target_gene",
            "pchembl_mean", "pchembl_sd", "theta_mean", "theta_sd",
            "g_predicted", "g_90_upper", "phen_cosine", "phen_centroid",
            "three_way_jsd", "i_novel_score", "eight_cell_tag",
            "admet_status", "roberts_ceiling_ok", "wet_lab_priority",
            "evidence_axes",
        ])

    # Build per-compound rows
    v6a_by_compound: dict[str, pd.DataFrame] = {}
    if v6a_pchembl is not None and "compound_name" in v6a_pchembl.columns:
        for c, g in v6a_pchembl.groupby("compound_name"):
            v6a_by_compound[str(c)] = g

    v6b_by_target: dict[str, pd.Series] = {}
    target_to_gene: dict[str, str] = {}
    if v6b_theta is not None and "target_uniprot" in v6b_theta.columns:
        for _, row in v6b_theta.iterrows():
            v6b_by_target[str(row["target_uniprot"])] = row
            if "gene" in row.index:
                target_to_gene[str(row["target_uniprot"])] = str(row.get("gene", ""))

    v7_by_compound: dict[str, tuple[float, float, float]] = {}
    if v7_posterior is not None:
        for c in getattr(v7_posterior, "compounds", []):
            v7_by_compound[c] = (
                v7_posterior.g_mean.get(c, float("nan")),
                float(v7_posterior.g_97p5.get(c, float("nan"))
                       - v7_posterior.g_2p5.get(c, float("nan"))) / 4.0,
                v7_posterior.g_90_upper.get(c, float("nan")),
            )

    v8_by_compound: dict[str, object] = {}
    if v8_posterior is not None:
        v8_by_compound = v8_posterior.by_compound()

    admet_by_compound: dict[str, str] = {}
    if admet_gates is not None and "compound_name" in admet_gates.columns:
        status_col = next((c for c in ("status", "admet_status", "verdict")
                           if c in admet_gates.columns), None)
        if status_col:
            for _, row in admet_gates.iterrows():
                admet_by_compound[str(row["compound_name"])] = str(row[status_col])

    rows: list[dict] = []
    compound_to_target = compound_to_target or {}
    for c in sorted(compounds):
        # Resolve target
        target_u = compound_to_target.get(c, "")
        if not target_u and c in v6a_by_compound:
            target_u = str(_safe_first(v6a_by_compound[c]["target_uniprot"], ""))
        gene = target_to_gene.get(target_u, "")

        # V6.A pchembl
        pchembl_m, pchembl_sd = float("nan"), float("nan")
        evidence: list[str] = []
        if c in v6a_by_compound:
            g6 = v6a_by_compound[c]
            pchembl_col = next((cc for cc in ("predicted_pkd", "pchembl_mean",
                                               "predicted_pchembl")
                                if cc in g6.columns), None)
            sd_col = next((cc for cc in ("predicted_pkd_sd", "pchembl_sd")
                           if cc in g6.columns), None)
            if pchembl_col:
                pchembl_m = float(g6[pchembl_col].mean())
                evidence.append("v6a")
            if sd_col:
                pchembl_sd = float(g6[sd_col].mean())

        # V6.B Cluster D
        theta_m, theta_sd = float("nan"), float("nan")
        if target_u in v6b_by_target:
            tr = v6b_by_target[target_u]
            theta_m = float(tr.get("theta_mean", float("nan")))
            if "theta_97p5" in tr.index and "theta_2p5" in tr.index:
                theta_sd = float(tr["theta_97p5"] - tr["theta_2p5"]) / 4.0
            evidence.append("v6b")

        # V7
        g_pred, g_sd, g_90_upper = float("nan"), float("nan"), float("nan")
        if c in v7_by_compound:
            g_pred, g_sd, g_90_upper = v7_by_compound[c]
            evidence.append("v7")

        # V8 phenotype
        phen_cos, phen_cent, jsd, novel, tag = (
            float("nan"), "", float("nan"), float("nan"), "no_evidence"
        )
        if c in v8_by_compound:
            v8e = v8_by_compound[c]
            phen_cos = v8e.phen_cosine
            phen_cent = v8e.phen_centroid
            jsd = v8e.three_way_jsd
            novel = v8e.i_novel_score
            tag = v8e.eight_cell_tag
            if np.isfinite(phen_cos):
                evidence.append("v8")
            # If V8 has joint_g_mean, prefer it over V7 standalone
            if np.isfinite(v8e.joint_g_mean):
                g_pred = v8e.joint_g_mean
                g_sd = v8e.joint_g_sd
                g_90_upper = v8e.joint_g_90_upper

        # Roberts ceiling check
        ceiling_ok = True
        if np.isfinite(g_90_upper) and g_90_upper > cfg.roberts_ceiling_g:
            ceiling_ok = False

        # Wet-lab priority
        if cfg.enforce_roberts_ceiling and not ceiling_ok:
            priority = 0.0
        else:
            g_val = g_pred if np.isfinite(g_pred) else 0.0
            jsd_val = jsd if np.isfinite(jsd) else 0.0
            ci_val = 4 * g_sd if np.isfinite(g_sd) else 1.0
            novel_val = novel if np.isfinite(novel) else 0.0
            priority = (cfg.weight_g * g_val
                        + cfg.weight_jsd * np.log1p(jsd_val)
                        - cfg.weight_ci_penalty * ci_val
                        + cfg.weight_novel * novel_val)

        rows.append({
            "compound": c,
            "target_uniprot": target_u,
            "target_gene": gene,
            "pchembl_mean": pchembl_m,
            "pchembl_sd": pchembl_sd,
            "theta_mean": theta_m,
            "theta_sd": theta_sd,
            "g_predicted": g_pred,
            "g_90_upper": g_90_upper,
            "phen_cosine": phen_cos,
            "phen_centroid": phen_cent,
            "three_way_jsd": jsd,
            "i_novel_score": novel,
            "eight_cell_tag": tag,
            "admet_status": admet_by_compound.get(c, "UNKNOWN"),
            "roberts_ceiling_ok": ceiling_ok,
            "wet_lab_priority": priority,
            "evidence_axes": ",".join(evidence) if evidence else "none",
        })

    df = pd.DataFrame(rows)
    df = df.sort_values("wet_lab_priority", ascending=False).reset_index(drop=True)
    df.insert(0, "rank", df.index + 1)
    if cfg.top_n and cfg.top_n > 0:
        df = df.head(cfg.top_n)
    return df


def render_v10_markdown_report(
    df: pd.DataFrame,
    output_path: Path | str,
    cfg: JointCompositionConfig | None = None,
) -> None:
    """Render the v10 wet-lab shortlist as a markdown report."""
    cfg = cfg or JointCompositionConfig()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    L: list[str] = []
    L.append("# Wet-Lab Shortlist v10 — Three-Factor Joint Composition (V6 × V7 × V8)")
    L.append("")
    L.append("Composes V6.A 4-head pchembl + V6.B Cluster D θ̄ + V7 effect-size + V8 πphen "
             "phenotype into a single ranked handoff. Pre-filtered by Roberts 2020 SMD "
             f"ceiling (no g > {cfg.roberts_ceiling_g} at 90% upper CrI). Annotated "
             "with 4-axis disagreement (V6.A multi-head + V6.B D_i + V8 three-way JSD + "
             "I_novel novel-mechanism score) and 8-cell classification.")
    L.append("")
    L.append("## Composition")
    L.append("")
    L.append("- **V6.A** — Multi-head DTI ensemble (MAMMAL + Tanimoto + MMAtt + PSICHIC + BALM) "
             "via Bayesian router with Venn-ABERS calibration")
    L.append("- **V6.B** — Cluster D Bayesian θ̄ (AHBA + OT Genetics L2G + cellxgene single-cell)")
    L.append("- **V7** — Clinical Effect-Size Translation: PBPK + PRISMA-anchored hierarchical "
             "Bayes; β_target gated multiplicatively by θ̄")
    L.append("- **V8** — πphen Perturbational Evidence: LINCS L1000 + JUMP-CP + iPSC-MEA + "
             "chemCPA imputation; cosine to 5-MoA cognition centroids")
    L.append("")
    L.append("## 8-cell disagreement legend")
    L.append("")
    L.append("| Tag | (T, G, P) | Interpretation |")
    L.append("|---|---|---|")
    L.append("| `agreement.all_high` | (H, H, H) | canonical positive (donepezil, MPH) |")
    L.append("| `target_true.phenotype_failed` | (H, H, L) | encenicline/intepirdine/pridopidine |")
    L.append("| `target.phenotype` | (H, L, H) | binding + functional, no genetics |")
    L.append("| `target_only` | (H, L, L) | binding artifact / off-pathway |")
    L.append("| `genetic.phenotype` | (L, H, H) | genetic + functional, no good binder |")
    L.append("| `genetic_only` | (L, H, L) | GWAS but no actionable binder |")
    L.append("| **`phenotype_only.novel_mechanism`** | **(L, L, H)** | **clemastine territory** |")
    L.append("| `no_evidence` | (L, L, L) | nothing across axes |")
    L.append("")

    if df.empty:
        L.append("⚠️ Empty shortlist — no V6/V7/V8 posteriors fed into composer.")
        output_path.write_text("\n".join(L), encoding="utf-8")
        return

    n_total = len(df)
    n_ceiling_ok = int(df["roberts_ceiling_ok"].sum())
    n_violations = n_total - n_ceiling_ok
    L.append("## Headline")
    L.append("")
    L.append(f"- Total ranked: **{n_total}** compounds")
    L.append(f"- Roberts ceiling PASS: **{n_ceiling_ok}** ({n_ceiling_ok/n_total:.0%})")
    L.append(f"- Roberts ceiling violations: **{n_violations}**")
    # 8-cell distribution
    tag_counts = df["eight_cell_tag"].value_counts()
    L.append("- 8-cell distribution:")
    for tag, n in tag_counts.items():
        marker = " ★" if tag == "phenotype_only.novel_mechanism" else ""
        L.append(f"  - `{tag}`: {n}{marker}")
    L.append("")
    L.append("## Top-25 by wet-lab priority")
    L.append("")
    L.append("| Rank | Compound | Target | g | g₉₀ | I_novel | 8-cell | AXES |")
    L.append("|---|---|---|---|---|---|---|---|")
    for _, r in df.head(25).iterrows():
        g_pred = f"{r['g_predicted']:+.2f}" if np.isfinite(r['g_predicted']) else "—"
        g_upper = f"{r['g_90_upper']:+.2f}" if np.isfinite(r['g_90_upper']) else "—"
        i_novel = f"{r['i_novel_score']:.2f}" if np.isfinite(r['i_novel_score']) else "—"
        L.append(f"| {r['rank']} | {r['compound']} | "
                 f"{r['target_gene']}/{r['target_uniprot']} | "
                 f"{g_pred} | {g_upper} | {i_novel} | "
                 f"{r['eight_cell_tag']} | {r['evidence_axes']} |")
    L.append("")
    L.append("## Honest caveats")
    L.append("")
    L.append("- v10 is the **architectural composition**. Real-data v10 awaits all 4 "
             "posteriors flowing (V6.A.4 Venn-ABERS shipped, V6.B.3 PyMC NUTS scaffold "
             "shipped, V7.3 effect-size scaffold shipped, V8.3 MOFA+ scaffold shipped). ")
    L.append("- The Roberts 2020 ceiling (g = 0.50 at 90% upper) is the HARD pre-filter "
             "per V4 §13.Y + V7.4 Gate 2.")
    L.append("- I_novel highlighting (L, L, H) cell compounds is the V8 publishable "
             "contribution — `phenotype_only.novel_mechanism` tag.")
    L.append("- Bayesian-copula correlation correction is a closed-form Gaussian "
             "approximation in V8.5 Stage 1; full PyMC NUTS in Stage 2.")
    L.append("")
    L.append("---")
    L.append("")
    L.append("Generated by `src/mammal_repurposing/fusion/joint_composition.py` via "
             "`scripts/56_v8_wet_lab_shortlist_v10.py`.")
    output_path.write_text("\n".join(L), encoding="utf-8")


def availability() -> dict[str, object]:
    """Probe v10 composer availability."""
    return {
        "available": True,
        "n_columns_output": 18,
        "roberts_ceiling_default": 0.50,
        "weights_default": {
            "g_mean": 1.0,
            "jsd": 0.25,
            "novelty": 0.60,
            "ci_penalty": 0.30,
        },
    }


# ===========================================================================
# V11 — Grid-based composition (compound × target repurposing hypotheses)
# ===========================================================================
#
# WHY V11 EXISTS: the v10 composer (above) collapsed every compound to a
# single target via `g["target_uniprot"].iloc[0]`, which — because the V6.A
# parquet is ordered with ACHE (P22303) first — assigned ACHE to all 298
# compounds. That destroyed the target dimension and produced a degenerate
# shortlist (everything → ACHE, all g ≈ +0.07, all violating the ceiling).
#
# The correct scientific unit is the **(compound, target) repurposing
# hypothesis**. V11 scores the FULL grid and never collapses the target axis.
#
# Per-pair score (transparent, fully real-signal — no stubs):
#   B(c,t)  = within-target percentile of V6.A predicted_pkd  ∈ [0,1]
#             ("is c a strong binder at t, relative to the library?")
#   R(t)    = σ(θ̄_t) cognition relevance from V6.B NUTS posterior ∈ (0,1)
#   R̃(t)   = min(1, R(t) / relevance_anchor)  (rescaled discriminative gate)
#   μ_m(t)  = PRISMA meta-analytic class-prior mean for t's mechanism class
#   g(c,t)  = μ_m(t) · B(c,t) · R̃(t)            predicted Hedges' g
#   For ANCHOR compounds (donepezil, MPH, ...) with a real V7 NUTS posterior,
#   g is OVERRIDDEN by the real g_mean at the compound's best-binding target.
#
# g_90_upper = g + 1.2816 · σ_g  (σ_g blends class-prior σ_m + engagement unc.)
# Roberts ceiling filters pairs with g_90_upper > 0.50.


@dataclass
class GridCompositionConfig:
    """Configuration for the V11 grid composer."""
    roberts_ceiling_g: float = 0.50
    enforce_roberts_ceiling: bool = True
    relevance_anchor: float = 0.62          # σ(θ̄) of the most cognition-relevant target
    binding_high_threshold: float = 0.60    # B percentile for 8-cell "T high"
    relevance_high_theta: float = 0.30      # θ̄ for 8-cell "G high"
    phenotype_high_threshold: float = 0.60  # transferability/cosine for 8-cell "P high"
    ci_z90: float = 1.2816                  # one-sided 90% normal quantile
    novelty_weight: float = 0.15            # small bonus for (L,L,H) novelty in priority
    top_n_pairs: int = 50


def _eight_cell_tag(b_high: bool, g_high: bool, p_high: bool) -> str:
    key = (b_high, g_high, p_high)
    return {
        (True, True, True):    "agreement.all_high",
        (True, True, False):   "target_true.phenotype_failed",
        (True, False, True):   "target.phenotype",
        (True, False, False):  "target_only",
        (False, True, True):   "genetic.phenotype",
        (False, True, False):  "genetic_only",
        (False, False, True):  "phenotype_only.novel_mechanism",
        (False, False, False): "no_evidence",
    }[key]


def compose_grid_shortlist_v11(
    v6a_pchembl: pd.DataFrame,
    v6b_theta: pd.DataFrame,
    target_class_map: dict[str, str],
    class_prior_table: dict[str, dict],
    *,
    v7_anchor_g: dict[str, tuple[float, float]] | None = None,
    phenotype_by_compound: dict[str, float] | None = None,
    target_gene_map: dict[str, str] | None = None,
    cfg: GridCompositionConfig | None = None,
    anchor_compound_target: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Compose the V11 wet-lab shortlist over the FULL (compound × target) grid.

    Args:
        v6a_pchembl: long DataFrame [compound_name, target_uniprot, predicted_pkd]
            — the real V6.A binding grid (no collapse).
        v6b_theta: DataFrame [target_uniprot, gene, theta_mean, theta_2p5,
            theta_97p5, w_pipeline] from the V6.B NUTS posterior.
        target_class_map: {uniprot: PRISMA mechanism class}.
        class_prior_table: {class: {"mean": μ_m, "sd": σ_m, ...}} from prisma_priors.
        v7_anchor_g: optional {base_compound_lower: (g_mean, g_90_upper)} real
            V7 NUTS posterior — overrides the class-prior estimate at the
            compound's KNOWN mechanism target (see anchor_compound_target).
        phenotype_by_compound: optional {compound_lower: phenotype score ∈ [0,1]}
            (e.g. chemCPA transferability or cosine-to-centroid). Drives the
            P axis of the 8-cell tag. Missing → P treated as low (unknown).
        target_gene_map: optional {uniprot: gene} for display.
        cfg: GridCompositionConfig.
        anchor_compound_target: authoritative {compound_lower: uniprot} mechanism
            map. The V7 clinical g is a COMPOUND-LEVEL effect size and must be
            placed at the drug's KNOWN mechanism target — NOT at MAMMAL's
            best-binding target, because MAMMAL is structurally blind to
            allosteric/transporter pharmacology and its binding argmax is
            biologically unreliable for exactly these drugs (e.g. it would
            rank rivastigmine's strongest binding at orexin, not ACHE).
            If a compound's known target is absent from the grid, its anchor
            g is simply not placed (the class-prior pathway runs instead).

    Returns:
        Per-(compound, target) DataFrame with binding/relevance/g/ceiling/8-cell
        columns, ranked by wet_lab_priority. NOT collapsed — one row per pair.
    """
    cfg = cfg or GridCompositionConfig()
    v7_anchor_g = v7_anchor_g or {}
    phenotype_by_compound = phenotype_by_compound or {}
    target_gene_map = target_gene_map or {}
    anchor_compound_target = {k.lower(): v for k, v in
                              (anchor_compound_target or {}).items()}

    pchembl_col = next((c for c in ("predicted_pkd", "pchembl_mean",
                                    "predicted_pchembl")
                        if c in v6a_pchembl.columns), None)
    if pchembl_col is None:
        raise ValueError("v6a_pchembl must have a predicted_pkd/pchembl column")

    # --- Within-target binding percentile B(c,t) ---
    df = v6a_pchembl[["compound_name", "target_uniprot", pchembl_col]].copy()
    df["B"] = (
        df.groupby("target_uniprot")[pchembl_col]
        .rank(method="average", pct=True)
    )

    # --- Relevance R(t) from V6.B ---
    theta_by_t: dict[str, float] = {}
    w_by_t: dict[str, float] = {}
    for _, r in v6b_theta.iterrows():
        u = str(r["target_uniprot"])
        theta_by_t[u] = float(r.get("theta_mean", float("nan")))
        w = r.get("w_pipeline", float("nan"))
        if not np.isfinite(w):
            # derive σ(θ̄) if w_pipeline absent
            tm = theta_by_t[u]
            w = 1.0 / (1.0 + np.exp(-tm)) if np.isfinite(tm) else 0.5
        w_by_t[u] = float(w)

    # --- Determine each compound's best-binding target (for V7 anchor override) ---
    best_target_by_compound: dict[str, str] = {}
    for c, g in df.groupby("compound_name"):
        best_target_by_compound[str(c)] = str(
            g.sort_values("B", ascending=False)["target_uniprot"].iloc[0]
        )

    rows: list[dict] = []
    for _, r in df.iterrows():
        c = str(r["compound_name"])
        t = str(r["target_uniprot"])
        B = float(r["B"])
        pkd = float(r[pchembl_col])
        R = w_by_t.get(t, 0.5)
        theta = theta_by_t.get(t, 0.0)
        R_tilde = min(1.0, R / cfg.relevance_anchor)

        cls = target_class_map.get(t)
        if cls is None:
            logger.warning("compose_grid_shortlist_v11: target %s has no class "
                           "mapping; using the weak FAILURE-class fallback", t)
            cls = "AMPA_pos_mod"
        prior = class_prior_table.get(cls)
        if prior is None:
            logger.warning("compose_grid_shortlist_v11: class %s has no prior; "
                           "using the weak high-variance fallback", cls)
            prior = {"mean": 0.05, "sd": 0.20}
        mu_m = float(prior["mean"])
        sd_m = float(prior.get("sd", 0.15))

        # Predicted g from class-prior pathway
        g_pred = mu_m * B * R_tilde
        # Uncertainty: class-prior σ scaled by engagement + a floor
        engagement = B * R_tilde
        g_sd = sd_m * (0.5 + 0.5 * engagement)
        source = "class_prior"

        # V7 anchor override (real meta-analytic g) at the compound's KNOWN
        # mechanism target. We use the authoritative compound→target map, NOT
        # MAMMAL's best-binding argmax (which is structurally unreliable for
        # the allosteric/transporter drugs that dominate this panel).
        c_lower = c.lower()
        known_t = anchor_compound_target.get(c_lower)
        anchor_target = known_t if known_t else best_target_by_compound.get(c)
        if c_lower in v7_anchor_g and anchor_target == t:
            g_anchor, g90_anchor = v7_anchor_g[c_lower]
            if np.isfinite(g_anchor):
                g_pred = float(g_anchor)
                g_sd = max(1e-3, (float(g90_anchor) - g_pred) / cfg.ci_z90) \
                    if np.isfinite(g90_anchor) else g_sd
                source = ("v7_nuts_anchor" if known_t
                          else "v7_nuts_anchor_bindingfallback")

        g_90 = g_pred + cfg.ci_z90 * g_sd
        ceiling_ok = not (np.isfinite(g_90) and g_90 > cfg.roberts_ceiling_g)

        # Phenotype axis
        phen = phenotype_by_compound.get(c_lower, float("nan"))
        p_high = bool(np.isfinite(phen) and phen >= cfg.phenotype_high_threshold)

        b_high = B >= cfg.binding_high_threshold
        g_high = theta >= cfg.relevance_high_theta
        tag = _eight_cell_tag(b_high, g_high, p_high)

        # Priority: clinically-meaningful g among ceiling-passing pairs,
        # with a small novelty bonus for the (L,L,H) cell.
        if cfg.enforce_roberts_ceiling and not ceiling_ok:
            priority = -1.0
        else:
            novelty = cfg.novelty_weight if tag == "phenotype_only.novel_mechanism" else 0.0
            priority = g_pred + novelty

        rows.append({
            "compound": c,
            "target_uniprot": t,
            "target_gene": target_gene_map.get(t, ""),
            "mechanism_class": cls,
            "predicted_pkd": pkd,
            "binding_percentile": B,
            "theta_mean": theta,
            "relevance_w": R,
            "class_prior_g": mu_m,
            "g_predicted": g_pred,
            "g_sd": g_sd,
            "g_90_upper": g_90,
            "g_source": source,
            "phenotype_score": phen,
            "eight_cell_tag": tag,
            "roberts_ceiling_ok": ceiling_ok,
            "wet_lab_priority": priority,
        })

    out = pd.DataFrame(rows)
    out = out.sort_values("wet_lab_priority", ascending=False).reset_index(drop=True)
    out.insert(0, "rank", out.index + 1)
    return out


def best_target_per_compound(grid_df: pd.DataFrame) -> pd.DataFrame:
    """Collapse the V11 grid to each compound's single best (highest-priority)
    target hypothesis. This is the clinician-facing 'which target should we
    test this drug against' view."""
    idx = grid_df.groupby("compound")["wet_lab_priority"].idxmax()
    best = grid_df.loc[idx].sort_values("wet_lab_priority", ascending=False)
    best = best.reset_index(drop=True)
    best["rank"] = best.index + 1
    return best
