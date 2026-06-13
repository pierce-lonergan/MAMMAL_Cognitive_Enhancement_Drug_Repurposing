"""Model-disagreement diagnosis report.

Renders the v2 research doc's §7 model-disagreement protocol as a visible
artifact. For each archetype, surface the compounds that fall into it so
reviewers can see WHERE clusters disagree and HOW the pipeline resolves.

Archetypes:
    1. mammal_strong_boltzina_weak  — sequence-only false positive; downrank
    2. boltzina_strong_mammal_weak  — structure-aware win; trust Boltzina if pLDDT≥70
    3. txgnn_strong_no_panel_hit    — off-panel mechanism; flag for review
    4. admet_gate_failed_other_strong — physical kill; do NOT bypass
    5. all_clusters_agree           — high-confidence shortlist members
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import yaml

from mammal_repurposing.config import PROJECT_ROOT

logger = logging.getLogger(__name__)

DEFAULT_WEIGHTS_PATH = PROJECT_ROOT / "configs" / "weights.yaml"


def classify_disagreements(
    provenance: pd.DataFrame,
    *,
    weights_path: Path | str = DEFAULT_WEIGHTS_PATH,
) -> pd.DataFrame:
    """Add a `disagreement_archetype` column based on per-cluster contributions."""
    with open(weights_path, encoding="utf-8") as _fh:
        cfg = yaml.safe_load(_fh)
    d = cfg.get("disagreement", {})
    pkd_strong = float(d.get("mammal_pkd_strong", 7.0))
    pkd_weak = float(d.get("mammal_pkd_weak", 5.5))
    bp_high = float(d.get("boltzina_binder_prob_high", 0.70))
    bp_low = float(d.get("boltzina_binder_prob_low", 0.30))
    tx_strong = float(d.get("txgnn_indication_strong", 0.50))

    df = provenance.copy()
    has_mammal = "mammal_best_pkd" in df.columns
    has_boltzina = "boltzina_best_binder_prob" in df.columns
    has_txgnn = "txgnn_max_indication" in df.columns

    def _arch(row) -> str:
        m_str = has_mammal and pd.notna(row.get("mammal_best_pkd")) and row["mammal_best_pkd"] >= pkd_strong
        m_weak = has_mammal and pd.notna(row.get("mammal_best_pkd")) and row["mammal_best_pkd"] <= pkd_weak
        b_str = has_boltzina and pd.notna(row.get("boltzina_best_binder_prob")) and row["boltzina_best_binder_prob"] >= bp_high
        b_weak = has_boltzina and pd.notna(row.get("boltzina_best_binder_prob")) and row["boltzina_best_binder_prob"] <= bp_low
        t_str = has_txgnn and pd.notna(row.get("txgnn_max_indication")) and row["txgnn_max_indication"] >= tx_strong
        panel_hit = has_mammal and pd.notna(row.get("mammal_best_pkd")) and row["mammal_best_pkd"] >= 6.0
        gate_cut = row.get("gate_status") == "CUT"

        if gate_cut:
            return "admet_gate_failed"
        if m_str and b_weak:
            return "mammal_strong_boltzina_weak"
        if b_str and m_weak:
            return "boltzina_strong_mammal_weak"
        if t_str and not panel_hit:
            return "txgnn_strong_no_panel_hit"
        if m_str and (b_str or not has_boltzina) and (t_str or not has_txgnn):
            return "all_clusters_agree"
        return "neutral"

    df["disagreement_archetype"] = df.apply(_arch, axis=1)
    return df


def render_markdown(provenance: pd.DataFrame, top_per_archetype: int = 10) -> str:
    classified = classify_disagreements(provenance)
    lines: list[str] = []
    lines.append("# Model-Disagreement Diagnosis Report")
    lines.append("")
    lines.append("Concrete realization of the v2 research doc's §7 model-disagreement protocol. "
                 "Each archetype below is a class of (compound) where clusters disagree; the "
                 "**Action** column states how the v2 pipeline resolves the disagreement.")
    lines.append("")
    counts = classified["disagreement_archetype"].value_counts()
    lines.append("## Archetype counts")
    lines.append("")
    lines.append("| Archetype | n compounds | Action |")
    lines.append("|---|---|---|")
    action_map = {
        "mammal_strong_boltzina_weak": "Downrank — sequence-only false positive",
        "boltzina_strong_mammal_weak": "Trust Boltzina if pose pLDDT ≥ 70",
        "txgnn_strong_no_panel_hit": "Flag for off-panel mechanism review",
        "admet_gate_failed": "Physical kill — do NOT bypass",
        "all_clusters_agree": "High-confidence shortlist member",
        "neutral": "No strong cross-cluster signal",
    }
    for arch in [
        "all_clusters_agree", "boltzina_strong_mammal_weak",
        "mammal_strong_boltzina_weak", "txgnn_strong_no_panel_hit",
        "admet_gate_failed", "neutral",
    ]:
        n = int(counts.get(arch, 0))
        lines.append(f"| `{arch}` | {n} | {action_map.get(arch, '—')} |")
    lines.append("")

    for arch in [
        "all_clusters_agree", "boltzina_strong_mammal_weak",
        "mammal_strong_boltzina_weak", "txgnn_strong_no_panel_hit",
    ]:
        subset = classified[classified["disagreement_archetype"] == arch]
        if subset.empty:
            continue
        lines.append(f"## {arch}")
        lines.append("")
        lines.append(f"_{action_map.get(arch, '')}_")
        lines.append("")
        cols_to_show = [
            c for c in [
                "compound_name", "evidence_tier",
                "mammal_best_target", "mammal_best_pkd",
                "boltzina_best_target", "boltzina_best_binder_prob",
                "txgnn_max_indication", "admet_score", "gate_status",
                "rrf_score",
            ]
            if c in subset.columns
        ]
        # Show the genuinely top-ranked members, not the input-order-first rows (the
        # provenance frame arrives in compounds-seed order; rrf_score is a merged column).
        ranked = (subset.sort_values("rrf_score", ascending=False)
                  if "rrf_score" in subset.columns else subset)
        head = ranked.head(top_per_archetype)[cols_to_show]
        lines.append(head.to_markdown(index=False))
        lines.append("")

    lines.append("---")
    lines.append("Generated by `provenance/disagreement_report.py`.")
    return "\n".join(lines)
