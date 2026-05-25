"""Per-candidate prose narratives explaining the funnel journey.

Template-based (no LLM). Generates a short paragraph per top-N candidate
describing WHICH clusters voted strongly, WHICH gates flagged, and WHAT
next experiment is recommended. Optional LLM rewrite layer for Phase 5+.
"""

from __future__ import annotations

import pandas as pd


def _format_value(name: str, val) -> str:
    if pd.isna(val):
        return "n/a"
    if isinstance(val, float):
        if name.endswith("_prob") or name.startswith("admet_") or "indication" in name:
            return f"{val:.2f}"
        return f"{val:.2f}"
    return str(val)


def _cluster_summary(row: pd.Series) -> list[str]:
    bits: list[str] = []
    if pd.notna(row.get("mammal_best_pkd")):
        bits.append(
            f"MAMMAL predicts pKd {row['mammal_best_pkd']:.2f} at "
            f"{row.get('mammal_best_target', '?')} (polypharm hits: "
            f"{int(row.get('mammal_polypharm_n', 0))})"
        )
    if pd.notna(row.get("boltzina_best_binder_prob")):
        bits.append(
            f"Boltzina pose-conditioned binder probability "
            f"{row['boltzina_best_binder_prob']:.2f} at "
            f"{row.get('boltzina_best_target', '?')}"
        )
    if pd.notna(row.get("admet_score")):
        bits.append(f"ADMET composite score {row['admet_score']:.2f}")
    if pd.notna(row.get("txgnn_max_indication")):
        bits.append(f"TxGNN indication probability {row['txgnn_max_indication']:.2f}")
    return bits


def _gate_summary(row: pd.Series) -> str:
    status = row.get("gate_status", "unknown")
    if status == "PASS":
        return "All ADMET hard gates passed."
    if status == "FLAG":
        flags = row.get("gates_flagged") or ""
        bypass = row.get("regulatory_bypass", False)
        suffix = " (regulatory bypass: approved drug)" if bypass else ""
        return f"Flagged but kept{suffix}; flags: {flags}."
    if status == "CUT":
        return f"CUT by gates: {row.get('gates_failed', '?')} — excluded from final ranking."
    return f"Gate status: {status}"


def _recommended_experiment(row: pd.Series) -> str:
    target = row.get("mammal_best_target") or row.get("boltzina_best_target")
    if not target:
        return "Recommended next: literature triage; no clear panel target."
    if row.get("evidence_tier") == "positive_control":
        return f"Recommended next: literature confirmation (already a positive control for {target})."
    if (row.get("mammal_best_pkd") or 0) >= 6.5:
        return (f"Recommended next: radioligand binding assay at {target} "
                f"(~$500-2000 via Eurofins/Cerep) to confirm predicted affinity.")
    return f"Recommended next: literature triage at {target} before further investment."


def render_narrative(
    provenance: pd.DataFrame,
    *,
    top_n: int = 20,
) -> str:
    if "rrf_score" not in provenance.columns:
        ordered = provenance.copy()
    else:
        ordered = provenance.sort_values("rrf_score", ascending=False)
    top = ordered.head(top_n)

    lines: list[str] = []
    lines.append("# Funnel Narrative — Top Candidates")
    lines.append("")
    lines.append(f"Per-compound prose explanation for the top {len(top)} candidates. "
                 "Each entry describes which clusters voted strongly, which gates flagged, "
                 "and what the recommended next experimental step is.")
    lines.append("")

    for idx, row in top.iterrows():
        rank = ordered.index.get_loc(idx) + 1
        lines.append(f"## #{rank}. {row['compound_name']}")
        lines.append("")
        lines.append(f"**Evidence tier**: `{row.get('evidence_tier', 'unknown')}`")
        lines.append("")
        lines.append("**Cluster contributions**:")
        for bit in _cluster_summary(row):
            lines.append(f"- {bit}")
        lines.append("")
        lines.append(f"**Gate status**: {_gate_summary(row)}")
        lines.append("")
        if pd.notna(row.get("rrf_score")):
            lines.append(f"**Fusion**: RRF score {row['rrf_score']:.4f}, "
                         f"{int(row.get('n_clusters_contributing', 0))} clusters contributed.")
            lines.append("")
        lines.append(_recommended_experiment(row))
        lines.append("")
        lines.append("---")
        lines.append("")
    return "\n".join(lines)
