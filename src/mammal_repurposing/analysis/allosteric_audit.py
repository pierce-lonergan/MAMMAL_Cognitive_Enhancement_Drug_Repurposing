"""Allosteric blindness audit for the MAMMAL DTI head.

MAMMAL's training data (BindingDB pKd labels) conflates orthosteric and
allosteric binding, with strong skew toward orthosteric ligands. We expect the
DTI head to systematically underrank allosteric ligands at their targets.

This audit defines, for each cognition-relevant allosteric target:
    (a) Known orthosteric ligands (expected to rank near the top)
    (b) Known allosteric ligands (the prediction we want to interrogate)
    (c) Optional: known non-binders or off-class ligands (control)

For each target, compute rank percentiles within MAMMAL's score distribution
for both classes. The audit passes a target if at least one allosteric ligand
ranks in the top quartile (top 25%) — anything worse suggests MAMMAL is
genuinely blind to that allosteric chemotype.

Targets covered (matched to our compound library):
    - CHRNA7 (alpha-7 nAChR PAMs): galantamine, encenicline, TC-5619
    - PDE4D (NAM): bpn14770 / zatolmilast (vs orthosteric rolipram)
    - GRIA1-4 (AMPA PAMs / ampakines): cx-717, cx-516, tulrampator, aniracetam
"""

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd


# Curated allosteric ligand panels. Names must match `compound_name` in dti_scores
# (case-insensitive). Orthosteric panel is used as a within-target reference.
ALLOSTERIC_AUDIT_PANELS: dict[str, dict[str, list[str]]] = {
    "P36544": {  # CHRNA7
        "gene": ["CHRNA7"],
        "orthosteric": [],  # full orthosteric ligands are not well-represented in our lib
        "allosteric": ["galantamine", "encenicline", "tc-5619"],
        "notes": "alpha-7 nAChR PAMs. Galantamine is the only approved member.",
    },
    "Q08499": {  # PDE4D
        "gene": ["PDE4D"],
        "orthosteric": ["rolipram"],
        "allosteric": ["bpn14770", "zatolmilast"],
        "notes": "PDE4D NAM (BPN14770) vs orthosteric (rolipram). "
                 "Direct head-to-head comparison.",
    },
    "P42261": {  # GRIA1
        "gene": ["GRIA1", "GRIA2", "GRIA3", "GRIA4"],
        "orthosteric": [],  # AMPA orthosteric ligands not in our library
        "allosteric": ["cx-717", "cx-516", "tulrampator", "aniracetam",
                       "piracetam", "oxiracetam", "pramiracetam"],
        "notes": "AMPA PAMs (ampakines). Aggregate across all 4 AMPA subunits in render.",
    },
}


@dataclass
class CompoundRanking:
    compound_name: str
    rank: int  # 1-indexed within target's distribution
    percentile: float  # 0-1; higher = better pKd
    pkd: float


@dataclass
class TargetAuditRow:
    target_uniprot: str
    target_gene: str
    notes: str
    orthosteric: list[CompoundRanking] = field(default_factory=list)
    allosteric: list[CompoundRanking] = field(default_factory=list)
    n_compounds_at_target: int = 0

    @property
    def allosteric_passes(self) -> bool:
        """Audit passes if ≥1 allosteric ligand is in top 25% of target's distribution."""
        return any(c.percentile >= 0.75 for c in self.allosteric)

    @property
    def gap_to_orthosteric(self) -> float | None:
        """Mean percentile gap (orthosteric - allosteric). Positive = MAMMAL prefers orthosteric."""
        if not self.allosteric or not self.orthosteric:
            return None
        ortho_mean = sum(c.percentile for c in self.orthosteric) / len(self.orthosteric)
        allo_mean = sum(c.percentile for c in self.allosteric) / len(self.allosteric)
        return ortho_mean - allo_mean


def _rank_within_target(
    scores: pd.DataFrame,
    target_uniprot: str,
    compound_names: list[str],
) -> list[CompoundRanking]:
    sub = scores[scores["target_uniprot"] == target_uniprot].copy()
    sub["rank"] = sub["predicted_pkd"].rank(method="max", ascending=True).astype(int)
    sub["percentile"] = sub["rank"] / len(sub)
    sub["name_lc"] = sub["compound_name"].str.lower().str.strip()

    found: list[CompoundRanking] = []
    for name in compound_names:
        match = sub[sub["name_lc"] == name.lower().strip()]
        if match.empty:
            continue
        row = match.iloc[0]
        found.append(CompoundRanking(
            compound_name=row["compound_name"],
            rank=int(row["rank"]),
            percentile=float(row["percentile"]),
            pkd=float(row["predicted_pkd"]),
        ))
    return found


def audit(scores: pd.DataFrame) -> list[TargetAuditRow]:
    """Run the allosteric audit across the configured panels."""
    rows: list[TargetAuditRow] = []
    for uniprot, panel in ALLOSTERIC_AUDIT_PANELS.items():
        # For multi-subunit targets (GRIA), audit each subunit separately
        for gene in panel["gene"]:
            # Resolve uniprot for THIS subunit gene; default to the panel's uniprot
            sub_uniprot = uniprot
            sub = scores[scores["target_gene"].str.upper() == gene.upper()]
            if sub.empty:
                continue
            sub_uniprot = sub["target_uniprot"].iloc[0]

            ortho = _rank_within_target(scores, sub_uniprot, panel["orthosteric"])
            allo = _rank_within_target(scores, sub_uniprot, panel["allosteric"])
            n_at_target = len(scores[scores["target_uniprot"] == sub_uniprot])

            if not ortho and not allo:
                continue

            rows.append(TargetAuditRow(
                target_uniprot=sub_uniprot,
                target_gene=gene,
                notes=panel["notes"],
                orthosteric=ortho,
                allosteric=allo,
                n_compounds_at_target=n_at_target,
            ))
    return rows


def render_markdown(rows: list[TargetAuditRow]) -> str:
    lines: list[str] = []
    n_pass = sum(1 for r in rows if r.allosteric_passes)
    overall = "PASS" if rows and all(r.allosteric_passes for r in rows) else "PARTIAL"
    if not rows:
        overall = "INSUFFICIENT_DATA"

    lines.append("# Allosteric Blindness Audit")
    lines.append("")
    lines.append(f"**Overall**: `{overall}` — {n_pass}/{len(rows)} targets pass allosteric gate (≥1 allosteric ligand in top 25%).")
    lines.append("")
    lines.append("Audit interrogates MAMMAL's known weakness: BindingDB-trained DTI heads systematically "
                 "underrank allosteric ligands relative to orthosteric. If allosteric ligands cluster in the "
                 "bottom half of their target's distribution, downstream candidate ranking should down-weight "
                 "those targets or apply target-specific recalibration.")
    lines.append("")

    for row in rows:
        verdict = "✅ pass" if row.allosteric_passes else "❌ fail"
        gap = row.gap_to_orthosteric
        gap_str = f"{gap:+.2f}" if gap is not None else "N/A (no orthosteric)"
        lines.append(f"## {row.target_gene} ({row.target_uniprot}) — {verdict}")
        lines.append("")
        lines.append(f"_{row.notes}_")
        lines.append("")
        lines.append(f"- Compounds scored at this target: **{row.n_compounds_at_target}**")
        lines.append(f"- Orthosteric vs allosteric percentile gap: **{gap_str}**")
        lines.append("")

        if row.orthosteric:
            lines.append("**Orthosteric ligands** (expected high pKd):")
            lines.append("")
            lines.append("| Compound | pKd | Rank | Percentile |")
            lines.append("|---|---|---|---|")
            for c in sorted(row.orthosteric, key=lambda x: -x.percentile):
                lines.append(f"| {c.compound_name} | {c.pkd:.2f} | {c.rank} | {c.percentile:.0%} |")
            lines.append("")

        if row.allosteric:
            lines.append("**Allosteric ligands** (interrogation target):")
            lines.append("")
            lines.append("| Compound | pKd | Rank | Percentile | Top-25%? |")
            lines.append("|---|---|---|---|---|")
            for c in sorted(row.allosteric, key=lambda x: -x.percentile):
                flag = "✅" if c.percentile >= 0.75 else "❌"
                lines.append(f"| {c.compound_name} | {c.pkd:.2f} | {c.rank} | {c.percentile:.0%} | {flag} |")
            lines.append("")

    lines.append("---")
    lines.append("Generated by `scripts/11_allosteric_audit.py`.")
    return "\n".join(lines)
