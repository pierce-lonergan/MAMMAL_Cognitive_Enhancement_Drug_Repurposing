"""V3 Diagnostic Protocol — diagnose the MAMMAL_ONLY_INVERTED finding.

Runs the cheap, CPU-only diagnostics from
research/4-tier/Diagnosing MAMMAL DTI Anti-Correlation.md while the Boltz
sweep continues in the background. Does NOT touch the GPU.

Diagnostics run:
  0. Prior-collapse sanity check (panel-wide; all 22 targets)
  1. Power analysis (Bonett-Wright + permutation CI)
  2. Diagnostic A — Murcko scaffold saturation
  3. Diagnostic B — pchembl distribution K-S + Wasserstein
  4. Diagnostic D — Tanimoto-to-known-actives vs MAMMAL pKd (THE most important)
  5. Lateral 6.1 — binding-mode mix (ChEMBL action_type)
  6. Lateral 6.2 — temporal stratification of ChEMBL truth

Targets:
  4 INVERTED (Q01959 SLC6A3, P23975 SLC6A2, Q12879 GRIN2A, Q13224 GRIN2B)
  2 STRONG controls (P21728 DRD1, O43613 HCRTR1)
  1 WEAK control (P22303 ACHE — ρ=+0.20 n=10)

Output: reports/pipeline/diagnostics_v1.md + reports/data/diagnostics_v1.parquet
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections import Counter
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.config import (  # noqa: E402
    DTI_SCORES_PARQUET,
    TARGETS_PARQUET,
)
from mammal_repurposing.diagnostics import (  # noqa: E402
    binding_mode_mix, distribution_shift, power_analysis, prior_collapse,
    scaffold_saturation, tanimoto_correlation, temporal_strat,
)
from mammal_repurposing.fetchers.chembl_sqlite import (  # noqa: E402
    chembl_actives_with_smiles_for_target,
    get_conn,
    per_target_pchembl_records,
)

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("v3_diagnose")

# Target panel split (from Phase A.7 calibration)
INVERTED = ["Q01959", "P23975", "Q12879", "Q13224"]
STRONG_CONTROLS = ["P21728", "O43613"]
WEAK_CONTROL = ["P22303"]
ALL_TARGETS = INVERTED + STRONG_CONTROLS + WEAK_CONTROL


def _smiles_to_inchikey(smi: str) -> str | None:
    if not isinstance(smi, str) or not smi:
        return None
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    return Chem.MolToInchiKey(m)


def join_library_to_truth(
    dti_grid: pd.DataFrame,
    target_uniprot: str,
    per_target_truth: pd.DataFrame,
) -> pd.DataFrame:
    """Return (smiles, predicted_pkd, inchikey, best_pchembl) for library compounds
    at this target that have a ChEMBL truth record."""
    lib = dti_grid[dti_grid["target_uniprot"] == target_uniprot].copy()
    lib["inchikey"] = lib["compound_smiles"].map(_smiles_to_inchikey)
    lib = lib.dropna(subset=["inchikey"])
    joined = lib.merge(
        per_target_truth[["inchikey", "best_pchembl"]],
        on="inchikey", how="left",
    )
    return joined.rename(columns={"compound_smiles": "smiles"})


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", type=Path, default=DTI_SCORES_PARQUET)
    parser.add_argument("--targets", type=Path, default=TARGETS_PARQUET)
    parser.add_argument("--out", type=Path, default=ROOT / "reports" / "pipeline" / "diagnostics_v1.md")
    parser.add_argument("--parquet-out", type=Path,
                        default=ROOT / "reports" / "data" / "diagnostics_v1.parquet")
    parser.add_argument("--n-perm", type=int, default=10000)
    parser.add_argument("--active-pchembl", type=float, default=8.0,
                        help="Active threshold for Diagnostic D (≥10nM = 8.0)")
    args = parser.parse_args()

    dti_grid = pd.read_parquet(args.scores)
    targets = pd.read_parquet(args.targets)
    gene_map = dict(zip(targets["uniprot"], targets["gene"]))

    # --- 0. Prior-collapse panel-wide -------------------------------------------
    logger.info("=== 0. Prior-collapse sanity check (panel-wide) ===")
    prior_df = prior_collapse.evaluate_panel(dti_grid)
    logger.info("Prior collapse verdicts: %s",
                Counter(prior_df["verdict"]).most_common())

    # --- 1. Power analysis on the joined (pred, truth) per target --------------
    logger.info("=== 1. Power analysis on the joined (pred, truth) per target ===")
    conn = get_conn()
    truth_cache: dict[str, pd.DataFrame] = {}
    joined_cache: dict[str, pd.DataFrame] = {}
    for u in ALL_TARGETS:
        truth_cache[u] = per_target_pchembl_records(u)
        joined_cache[u] = join_library_to_truth(dti_grid, u, truth_cache[u])
    power_input = {
        u: {
            "gene": gene_map.get(u, "?"),
            "pred": joined_cache[u].dropna(subset=["best_pchembl"])["predicted_pkd"].to_numpy(),
            "truth": joined_cache[u].dropna(subset=["best_pchembl"])["best_pchembl"].to_numpy(),
        }
        for u in ALL_TARGETS
    }
    power_df = power_analysis.power_panel(power_input)
    logger.info("Power verdicts: %s", dict(Counter(power_df["verdict"])))

    # --- 2. Diagnostic A — scaffold saturation ---------------------------------
    logger.info("=== 2. Diagnostic A — Murcko scaffold saturation ===")
    scaffold_rows = []
    actives_cache: dict[str, pd.DataFrame] = {}
    for u in ALL_TARGETS:
        gene = gene_map.get(u, "?")
        actives = chembl_actives_with_smiles_for_target(u, min_pchembl=args.active_pchembl)
        actives_cache[u] = actives
        lib_smi = joined_cache[u]["smiles"].dropna().tolist()
        active_smi = actives["canonical_smiles"].dropna().tolist()
        logger.info("  %s (%s): library=%d, ChEMBL pchembl≥%.1f actives=%d",
                    gene, u, len(lib_smi), args.active_pchembl, len(active_smi))
        if not lib_smi or not active_smi:
            continue
        r = scaffold_saturation.scaffold_overlap(u, gene, lib_smi, active_smi)
        scaffold_rows.append({
            "target_uniprot": r.target_uniprot, "gene": r.target_gene,
            "n_library": r.n_library, "n_chembl_actives": r.n_chembl_actives,
            "library_in_top_scaffold_pct": r.library_in_top_scaffold_pct,
            "library_unique_scaffolds": r.library_unique_scaffolds,
            "top_scaffold_smiles": r.top_chembl_scaffold_smiles,
            "scaffold_decision": r.decision,
        })
    scaffold_df = pd.DataFrame(scaffold_rows)

    # --- 3. Diagnostic B — pchembl distribution shift --------------------------
    logger.info("=== 3. Diagnostic B — pchembl distribution shift ===")
    dist_rows = []
    for u in ALL_TARGETS:
        gene = gene_map.get(u, "?")
        truth = truth_cache[u]
        lib_pchembl = joined_cache[u].dropna(subset=["best_pchembl"])["best_pchembl"].to_numpy()
        chembl_all = truth["best_pchembl"].to_numpy()
        if len(lib_pchembl) < 3 or len(chembl_all) < 3:
            continue
        r = distribution_shift.diagnose(u, gene, lib_pchembl, chembl_all)
        dist_rows.append({
            "target_uniprot": r.target_uniprot, "gene": r.target_gene,
            "n_library_with_truth": r.n_library, "n_chembl_all": r.n_chembl_all,
            "ks_stat": r.ks_stat, "ks_p": r.ks_pvalue, "wasserstein": r.wasserstein,
            "lib_median_pchembl": r.library_pchembl_median,
            "chembl_median_pchembl": r.chembl_pchembl_median,
            "dist_decision": r.decision,
        })
    dist_df = pd.DataFrame(dist_rows)

    # --- 4. Diagnostic D — Tanimoto to known actives ---------------------------
    logger.info("=== 4. Diagnostic D — Tanimoto-to-known-actives vs MAMMAL pKd ===")
    tani_rows = []
    for u in ALL_TARGETS:
        gene = gene_map.get(u, "?")
        actives_smi = actives_cache[u]["canonical_smiles"].dropna().tolist()
        lib_df = joined_cache[u][["smiles", "predicted_pkd", "best_pchembl"]].copy()
        r = tanimoto_correlation.diagnose(u, gene, lib_df, actives_smi)
        tani_rows.append({
            "target_uniprot": r.target_uniprot, "gene": r.target_gene,
            "n_library": r.n_library, "n_known_actives": r.n_known_actives,
            "spearman_mammal_vs_tanimoto": r.spearman_mammal_vs_tanimoto,
            "pearson_mammal_vs_tanimoto": r.pearson_mammal_vs_tanimoto,
            "spearman_mammal_vs_truth": r.spearman_mammal_vs_truth,
            "mean_max_tanimoto": r.mean_max_tanimoto,
            "tanimoto_decision": r.decision,
        })
    tani_df = pd.DataFrame(tani_rows)

    # --- 5. Lateral 6.1 — binding-mode mix -------------------------------------
    logger.info("=== 5. Lateral 6.1 — binding-mode mix (ChEMBL action_type) ===")
    mode_rows = []
    for u in ALL_TARGETS:
        gene = gene_map.get(u, "?")
        try:
            at_df = binding_mode_mix.fetch_action_types_for_target(conn, u)
        except Exception as e:
            logger.warning("  %s: action_type fetch failed: %s", u, e)
            continue
        lib_keys = joined_cache[u]["inchikey"].dropna().tolist()
        r = binding_mode_mix.diagnose(u, gene, lib_keys, at_df)
        mode_rows.append({
            "target_uniprot": r.target_uniprot, "gene": r.target_gene,
            "n_lib_with_action": r.n_library_with_action_type,
            "top_action": r.top_action_type, "top_action_pct": r.top_action_type_pct,
            "n_allosteric": r.n_allosteric, "pct_allosteric": r.pct_allosteric,
            "binding_mode_verdict": r.verdict,
        })
    mode_df = pd.DataFrame(mode_rows)

    # --- 6. Lateral 6.2 — temporal stratification ------------------------------
    logger.info("=== 6. Lateral 6.2 — temporal stratification (split year=2015) ===")
    temp_rows = []
    for u in ALL_TARGETS:
        gene = gene_map.get(u, "?")
        try:
            year_df = temporal_strat.fetch_chembl_year_pchembl(conn, u)
        except Exception as e:
            logger.warning("  %s: year fetch failed: %s", u, e)
            continue
        lib = joined_cache[u][["smiles", "predicted_pkd", "inchikey"]].dropna()
        r = temporal_strat.temporal_split(u, gene, lib, year_df, split_year=2015)
        temp_rows.append({
            "target_uniprot": r.target_uniprot, "gene": r.target_gene,
            "split_year": r.split_year, "n_pre": r.n_pre, "n_post": r.n_post,
            "pre_rho": r.pre_rho, "post_rho": r.post_rho,
            "pre_median_pchembl": r.pre_median_pchembl,
            "post_median_pchembl": r.post_median_pchembl,
            "temporal_verdict": r.verdict,
        })
    temp_df = pd.DataFrame(temp_rows)

    # --- Merge into one results parquet ----------------------------------------
    aggregate = (
        power_df.rename(columns={"verdict": "power_verdict"})
        .merge(scaffold_df.drop(columns=["gene"]), on="target_uniprot", how="left")
        .merge(dist_df.drop(columns=["gene"]), on="target_uniprot", how="left")
        .merge(tani_df.drop(columns=["gene"]), on="target_uniprot", how="left")
        .merge(mode_df.drop(columns=["gene"]), on="target_uniprot", how="left")
        .merge(temp_df.drop(columns=["gene"]), on="target_uniprot", how="left")
    )
    args.parquet_out.parent.mkdir(parents=True, exist_ok=True)
    aggregate.to_parquet(args.parquet_out, index=False)
    prior_df.to_parquet(args.parquet_out.with_name("diagnostics_v1_prior_collapse.parquet"),
                       index=False)

    # --- Render markdown report ------------------------------------------------
    md = _render_markdown(prior_df, power_df, scaffold_df, dist_df,
                          tani_df, mode_df, temp_df, gene_map)
    args.out.write_text(md, encoding="utf-8")
    logger.info("Wrote %s and %s", args.out, args.parquet_out)
    return 0


def _render_markdown(prior_df, power_df, scaffold_df, dist_df,
                     tani_df, mode_df, temp_df, gene_map) -> str:
    L: list[str] = []
    L.append("# V3 Diagnostic Protocol — Investigating MAMMAL_ONLY_INVERTED")
    L.append("")
    L.append("**Source**: `research/4-tier/Diagnosing MAMMAL DTI Anti-Correlation.md`.")
    L.append("")
    L.append("All diagnostics CPU-only; GPU stays free for the WSL2 Boltz sweep.")
    L.append("")

    # 0. Prior collapse — the headline finding
    L.append("## 0. Prior-collapse sanity check (panel-wide)")
    L.append("")
    L.append("MAMMAL's training prior: `norm_y_mean = 5.794`, `norm_y_std = 1.338`.")
    L.append("If predictions cluster tightly around the prior mean with std << training "
             "std, the \"ranking\" within that target is noise, not learned signal.")
    L.append("")
    L.append("| Target | Gene | n | pred_mean | pred_std | range | IQR | collapse vs SD=1.34 | Verdict |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    for _, r in prior_df.iterrows():
        L.append(f"| {r['target_uniprot']} | {r['target_gene']} | {int(r['n'])} | "
                 f"{r['pred_mean']:.3f} | {r['pred_std']:.4f} | "
                 f"{r['pred_range']:.3f} | {r['iqr']:.3f} | "
                 f"**{r['collapse_ratio']:.1f}×** | `{r['verdict']}` |")
    L.append("")
    severe = (prior_df["verdict"] == "SEVERE").sum()
    moderate = (prior_df["verdict"] == "MODERATE").sum()
    L.append(f"**Headline**: {severe}/{len(prior_df)} targets show **SEVERE** prior collapse "
             f"(pred std < 1/10 of training std). {moderate} are MODERATE. ")
    L.append("This reframes the entire calibration: per-target Spearman ρ values are "
             "computed over predictions that span <0.5 log unit, so the rank order is "
             "noise-driven for most compounds. Even the STRONG-control targets are collapsed; "
             "DRD1's ρ = +0.31 is a real signal extracted from a narrow band.")
    L.append("")

    # 1. Power analysis
    L.append("## 1. Power analysis (Bonett-Wright Fisher-z + permutation)")
    L.append("")
    L.append("Per Bonett & Wright 2000 *Psychometrika* — is the observed ρ distinguishable "
             "from zero at the joined sample size?")
    L.append("")
    L.append("| Target | Gene | n | ρ | Fisher-z 95% CI | perm p | perm 95% null CI | Distinguishable? | Verdict |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    for _, r in power_df.iterrows():
        rho = f"{r['rho']:+.2f}" if not pd.isna(r['rho']) else "—"
        L.append(f"| {r['target_uniprot']} | {r['gene']} | {int(r['n'])} | {rho} | "
                 f"{r['fisher_ci']} | {r['perm_p']:.3f} | {r['perm_null_ci']} | "
                 f"{'✅' if r['distinguishable_from_zero'] else '❌'} | `{r['verdict']}` |")
    L.append("")

    # 2. Scaffold saturation
    L.append("## 2. Diagnostic A — Murcko scaffold saturation")
    L.append("")
    L.append("For each target, what fraction of the library matches the most-common "
             "generic Bemis-Murcko scaffold among ChEMBL high-affinity binders (pchembl ≥ 8.0)?")
    L.append("")
    L.append("| Target | Gene | n_lib | n_chembl_actives | lib-in-top-scaffold % | unique_lib_scaffolds | Decision |")
    L.append("|---|---|---|---|---|---|---|")
    for _, r in scaffold_df.iterrows():
        pct = f"{r['library_in_top_scaffold_pct']:.1f}%" if not pd.isna(r['library_in_top_scaffold_pct']) else "—"
        L.append(f"| {r['target_uniprot']} | {r['gene']} | {int(r['n_library'])} | "
                 f"{int(r['n_chembl_actives'])} | **{pct}** | "
                 f"{int(r['library_unique_scaffolds'])} | `{r['scaffold_decision']}` |")
    L.append("")
    L.append("Routing per research doc: >60% = `rank_resolution_loss` (Scenario 2, LoRA worth); "
             "<25% = `manifold_mismatch` (Scenario 1, ensemble worth); else = ambiguous.")
    L.append("")

    # 3. Distribution shift
    L.append("## 3. Diagnostic B — pchembl distribution shift (K-S + Wasserstein)")
    L.append("")
    L.append("Is the library's pchembl distribution at this target consistent with ChEMBL's full distribution?")
    L.append("")
    L.append("| Target | Gene | n_lib_w_truth | n_chembl_all | K-S | K-S p | Wasserstein | lib_med | chembl_med | Decision |")
    L.append("|---|---|---|---|---|---|---|---|---|---|")
    for _, r in dist_df.iterrows():
        L.append(f"| {r['target_uniprot']} | {r['gene']} | {int(r['n_library_with_truth'])} | "
                 f"{int(r['n_chembl_all'])} | {r['ks_stat']:.3f} | {r['ks_p']:.3g} | "
                 f"{r['wasserstein']:.3f} | {r['lib_median_pchembl']:.2f} | "
                 f"{r['chembl_median_pchembl']:.2f} | `{r['dist_decision']}` |")
    L.append("")

    # 4. Tanimoto correlation — THE most important
    L.append("## 4. Diagnostic D — Tanimoto-to-known-actives vs MAMMAL pKd ★")
    L.append("")
    L.append("**Highest-value diagnostic.** ρ(MAMMAL pred, max Tanimoto-to-pchembl≥8 actives):")
    L.append("")
    L.append("- ρ > +0.30 → model rewards structural similarity correctly. Inversion vs ChEMBL is "
             "  activity-cliff driven within the cluster → **Scenario 2 (LoRA worth)**")
    L.append("- −0.20 < ρ < +0.30 → model has no usable signal → **Scenario 1 (manifold mismatch)**")
    L.append("- ρ < −0.20 → **ACTIVE INVERSION**: model penalises the right structural class → "
             "  **Scenario 4 (label-sign error)** — audit BindingDB rows BEFORE LoRA.")
    L.append("")
    L.append("| Target | Gene | n_lib | n_actives | ρ(MAMMAL, Tanimoto) | ρ(MAMMAL, truth) for cross-ref | mean max Tanimoto | Decision |")
    L.append("|---|---|---|---|---|---|---|---|")
    for _, r in tani_df.iterrows():
        rt = f"{r['spearman_mammal_vs_tanimoto']:+.2f}" if not pd.isna(r['spearman_mammal_vs_tanimoto']) else "—"
        rtruth = f"{r['spearman_mammal_vs_truth']:+.2f}" if not pd.isna(r['spearman_mammal_vs_truth']) else "—"
        mean_tan = f"{r['mean_max_tanimoto']:.3f}" if not pd.isna(r['mean_max_tanimoto']) else "—"
        L.append(f"| {r['target_uniprot']} | {r['gene']} | {int(r['n_library'])} | "
                 f"{int(r['n_known_actives'])} | **{rt}** | {rtruth} | "
                 f"{mean_tan} | `{r['tanimoto_decision']}` |")
    L.append("")

    # 5. Binding mode mix
    L.append("## 5. Lateral 6.1 — binding-mode mix (ChEMBL action_type)")
    L.append("")
    L.append("Is the target dominated by a single binding mode (orthosteric inhibitor), "
             "or is it a mix that includes allosteric pharmacology MAMMAL's single-chain "
             "sequence cannot represent (e.g., ifenprodil-class at GluN1/GluN2B interface)?")
    L.append("")
    L.append("| Target | Gene | n_with_action | top action | top % | n_allosteric | % allosteric | Verdict |")
    L.append("|---|---|---|---|---|---|---|---|")
    for _, r in mode_df.iterrows():
        ta = f"{r['top_action_pct']:.1f}%" if not pd.isna(r['top_action_pct']) else "—"
        pa = f"{r['pct_allosteric']:.1f}%" if not pd.isna(r['pct_allosteric']) else "—"
        L.append(f"| {r['target_uniprot']} | {r['gene']} | {int(r['n_lib_with_action'])} | "
                 f"{r['top_action']} | {ta} | {int(r['n_allosteric'])} | {pa} | "
                 f"`{r['binding_mode_verdict']}` |")
    L.append("")

    # 6. Temporal stratification
    L.append("## 6. Lateral 6.2 — temporal stratification (split year 2015)")
    L.append("")
    L.append("Are ChEMBL records driving the inversion clustered in post-2015 chemistry "
             "MAMMAL's pre-2018 BindingDB training never saw?")
    L.append("")
    L.append("| Target | Gene | n_pre | n_post | pre ρ | post ρ | pre_med pchembl | post_med pchembl | Verdict |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    for _, r in temp_df.iterrows():
        pre = f"{r['pre_rho']:+.2f}" if not pd.isna(r['pre_rho']) else "—"
        post = f"{r['post_rho']:+.2f}" if not pd.isna(r['post_rho']) else "—"
        L.append(f"| {r['target_uniprot']} | {r['gene']} | {int(r['n_pre'])} | {int(r['n_post'])} | "
                 f"**{pre}** | **{post}** | {r['pre_median_pchembl']:.2f} | "
                 f"{r['post_median_pchembl']:.2f} | `{r['temporal_verdict']}` |")
    L.append("")

    L.append("## Summary verdict table")
    L.append("")
    L.append("| Target | Gene | n | ρ | power | scaffold | distribution | Tanimoto | mode | temporal | Final routing |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|")
    merged = (power_df[["target_uniprot", "gene", "n", "rho", "verdict"]]
              .rename(columns={"verdict": "power_v"}))
    if not scaffold_df.empty:
        merged = merged.merge(scaffold_df[["target_uniprot", "scaffold_decision"]],
                              on="target_uniprot", how="left")
    if not dist_df.empty:
        merged = merged.merge(dist_df[["target_uniprot", "dist_decision"]],
                              on="target_uniprot", how="left")
    if not tani_df.empty:
        merged = merged.merge(tani_df[["target_uniprot", "tanimoto_decision"]],
                              on="target_uniprot", how="left")
    if not mode_df.empty:
        merged = merged.merge(mode_df[["target_uniprot", "binding_mode_verdict"]],
                              on="target_uniprot", how="left")
    if not temp_df.empty:
        merged = merged.merge(temp_df[["target_uniprot", "temporal_verdict"]],
                              on="target_uniprot", how="left")
    for _, r in merged.iterrows():
        rho = f"{r['rho']:+.2f}" if not pd.isna(r['rho']) else "—"
        # Heuristic final routing
        if r.get("power_v") == "NOISE":
            routing = "🟡 EXPAND n FIRST (S5)"
        elif r.get("binding_mode_verdict") == "allosteric_dominant":
            routing = "🔴 DEPRECATE / Boltz dimer (S3)"
        elif r.get("tanimoto_decision") == "correctly_correlated" and r.get("scaffold_decision") == "rank_resolution_loss":
            routing = "🟢 LoRA worth (S2)"
        elif r.get("tanimoto_decision") == "systematic_inversion":
            routing = "🟠 AUDIT BindingDB (S4)"
        elif r.get("tanimoto_decision") == "pure_noise":
            routing = "🔵 ENSEMBLE worth (S1)"
        else:
            routing = "❔ See diagnostics above"
        L.append(f"| {r['target_uniprot']} | {r['gene']} | {int(r['n'])} | {rho} | "
                 f"`{r.get('power_v', '?')}` | `{r.get('scaffold_decision', '—')}` | "
                 f"`{r.get('dist_decision', '—')}` | `{r.get('tanimoto_decision', '—')}` | "
                 f"`{r.get('binding_mode_verdict', '—')}` | `{r.get('temporal_verdict', '—')}` | "
                 f"**{routing}** |")
    L.append("")
    L.append("---")
    L.append("")
    L.append("_Scenario legend per research doc §3: S1=manifold mismatch, S2=rank-resolution loss "
             "(LoRA), S3=representational gap (deprecate), S4=active inversion (label bug), "
             "S5=insufficient n._")
    L.append("")
    L.append("Generated by `scripts/30_v3_diagnose_inverted.py`.")
    return "\n".join(L)


if __name__ == "__main__":
    raise SystemExit(main())
