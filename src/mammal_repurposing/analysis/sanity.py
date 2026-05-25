"""Positive-control sanity gate for the MAMMAL DTI score grid.

For each (target, expected_compound_list) in :data:`config.POSITIVE_CONTROLS`,
compute the rank percentile of the named compounds among that target's score
distribution. A target PASSES if at least one named compound ranks in the top
:data:`config.POSITIVE_CONTROL_TOP_PERCENTILE` of the distribution.

Also checks that negative-control compounds (peripheral-only drugs) don't
surface in any target's top 5% — that would suggest the model is producing
noise rather than signal.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from mammal_repurposing.config import (
    NEGATIVE_CONTROL_FLAG_PERCENTILE,
    POSITIVE_CONTROL_TOP_PERCENTILE,
    POSITIVE_CONTROLS,
)

logger = logging.getLogger(__name__)


@dataclass
class TargetCheck:
    target_uniprot: str
    target_gene: str
    expected_compounds: list[str]
    found_compounds: list[tuple[str, float, float]]  # (name, pkd, percentile)
    passed: bool


@dataclass
class NegativeControlHit:
    target_uniprot: str
    target_gene: str
    compound_name: str
    pkd: float
    percentile: float


@dataclass
class SanityReport:
    target_checks: list[TargetCheck] = field(default_factory=list)
    negative_hits: list[NegativeControlHit] = field(default_factory=list)
    pkd_summary: dict[str, float] = field(default_factory=dict)
    total_pairs: int = 0

    @property
    def passed(self) -> bool:
        """Overall gate: all positive-control targets pass AND no negative-control hits."""
        return all(t.passed for t in self.target_checks) and not self.negative_hits

    @property
    def n_targets_pass(self) -> int:
        return sum(1 for t in self.target_checks if t.passed)


def _rank_percentile(target_df: pd.DataFrame, compound_name: str) -> tuple[float, float] | None:
    """Return (predicted_pkd, percentile) for a compound within a target's distribution.

    Percentile is in [0, 1] where 1.0 == highest pKd. Returns None if compound
    not present in target_df.
    """
    name_lc = compound_name.lower().strip()
    matches = target_df[target_df["compound_name"].str.lower().str.strip() == name_lc]
    if matches.empty:
        return None
    pkd = float(matches["predicted_pkd"].iloc[0])
    # Rank ascending then convert to percentile [0, 1] — higher pkd = higher percentile.
    ranks = target_df["predicted_pkd"].rank(method="max", ascending=True)
    rank_of_compound = float(ranks[matches.index[0]])
    percentile = rank_of_compound / len(target_df)
    return pkd, percentile


def check_positive_controls(
    scores: pd.DataFrame,
    *,
    top_percentile: float = POSITIVE_CONTROL_TOP_PERCENTILE,
    controls: dict[str, list[str]] | None = None,
) -> list[TargetCheck]:
    """For each entry in POSITIVE_CONTROLS, check if any expected compound is
    in the top ``top_percentile`` (default 0.20 = top 20%) of that target's scores."""
    controls = controls if controls is not None else POSITIVE_CONTROLS
    threshold = 1.0 - top_percentile

    results: list[TargetCheck] = []
    for uniprot, expected in controls.items():
        sub = scores[scores["target_uniprot"] == uniprot]
        if sub.empty:
            logger.warning("Sanity: no scores for target %s in input.", uniprot)
            continue
        gene = sub["target_gene"].iloc[0]
        found: list[tuple[str, float, float]] = []
        for compound in expected:
            r = _rank_percentile(sub, compound)
            if r is None:
                continue
            pkd, pct = r
            found.append((compound, pkd, pct))
        passed = any(pct >= threshold for _, _, pct in found)
        results.append(TargetCheck(uniprot, gene, expected, found, passed))
    return results


def check_negative_controls(
    scores: pd.DataFrame,
    compounds: pd.DataFrame,
    *,
    flag_percentile: float = NEGATIVE_CONTROL_FLAG_PERCENTILE,
) -> list[NegativeControlHit]:
    """Flag any (target, negative_control_compound) pair that ranks in top
    ``flag_percentile`` (default 5%) of that target's distribution."""
    neg_names = (
        compounds[compounds["evidence_tier"] == "negative_control"]["name"]
        .str.lower().str.strip().tolist()
    )
    if not neg_names:
        logger.info("Sanity: no negative-control compounds present in compounds parquet.")
        return []

    threshold = 1.0 - flag_percentile
    scores = scores.copy()
    scores["_name_lc"] = scores["compound_name"].str.lower().str.strip()
    scores["_pct"] = scores.groupby("target_uniprot")["predicted_pkd"].rank(
        method="max", pct=True, ascending=True
    )

    hits: list[NegativeControlHit] = []
    flagged = scores[(scores["_name_lc"].isin(neg_names)) & (scores["_pct"] >= threshold)]
    for _, row in flagged.iterrows():
        hits.append(
            NegativeControlHit(
                target_uniprot=row["target_uniprot"],
                target_gene=row["target_gene"],
                compound_name=row["compound_name"],
                pkd=float(row["predicted_pkd"]),
                percentile=float(row["_pct"]),
            )
        )
    return hits


def build_report(
    scores: pd.DataFrame,
    compounds: pd.DataFrame,
    *,
    top_percentile: float = POSITIVE_CONTROL_TOP_PERCENTILE,
    neg_flag_percentile: float = NEGATIVE_CONTROL_FLAG_PERCENTILE,
) -> SanityReport:
    """Compute the full sanity report (pos + neg + pKd summary)."""
    target_checks = check_positive_controls(
        scores, top_percentile=top_percentile
    )
    negative_hits = check_negative_controls(
        scores, compounds, flag_percentile=neg_flag_percentile
    )
    pkd_summary = {
        "min": float(scores["predicted_pkd"].min()),
        "p25": float(scores["predicted_pkd"].quantile(0.25)),
        "median": float(scores["predicted_pkd"].median()),
        "p75": float(scores["predicted_pkd"].quantile(0.75)),
        "max": float(scores["predicted_pkd"].max()),
        "mean": float(scores["predicted_pkd"].mean()),
        "std": float(scores["predicted_pkd"].std()),
    }
    return SanityReport(
        target_checks=target_checks,
        negative_hits=negative_hits,
        pkd_summary=pkd_summary,
        total_pairs=len(scores),
    )


def render_markdown(
    report: SanityReport,
    polypharm: pd.DataFrame | None = None,
    *,
    top_percentile: float = POSITIVE_CONTROL_TOP_PERCENTILE,
    neg_flag_percentile: float = NEGATIVE_CONTROL_FLAG_PERCENTILE,
) -> str:
    """Render the report as Markdown for ``sanity_report.md``."""
    lines: list[str] = []
    overall = "PASS" if report.passed else "FAIL"
    lines.append(f"# Sanity Report — {overall}")
    lines.append("")
    lines.append(f"- Total pairs scored: **{report.total_pairs:,}**")
    lines.append(
        f"- Positive-control targets passing: **{report.n_targets_pass}/{len(report.target_checks)}** "
        f"(threshold: any named compound in top {top_percentile:.0%})"
    )
    lines.append(
        f"- Negative-control flagged hits in top {neg_flag_percentile:.0%}: "
        f"**{len(report.negative_hits)}**"
    )
    lines.append("")

    lines.append("## pKd Distribution Summary")
    lines.append("")
    s = report.pkd_summary
    lines.append("| Statistic | Value |")
    lines.append("|---|---|")
    lines.append(f"| min    | {s['min']:.3f} |")
    lines.append(f"| p25    | {s['p25']:.3f} |")
    lines.append(f"| median | {s['median']:.3f} |")
    lines.append(f"| p75    | {s['p75']:.3f} |")
    lines.append(f"| max    | {s['max']:.3f} |")
    lines.append(f"| mean   | {s['mean']:.3f} |")
    lines.append(f"| std    | {s['std']:.3f} |")
    lines.append("")

    lines.append("## Positive-Control Checks")
    lines.append("")
    lines.append("| Target (gene) | UniProt | Expected | Found (pKd / percentile) | Passed |")
    lines.append("|---|---|---|---|---|")
    for c in report.target_checks:
        found_str = (
            "; ".join(f"{n} ({pkd:.2f} / {pct:.0%})" for n, pkd, pct in c.found_compounds)
            if c.found_compounds
            else "_(none of the named compounds present)_"
        )
        expected_str = ", ".join(c.expected_compounds)
        flag = "✅" if c.passed else "❌"
        lines.append(
            f"| {c.target_gene} | {c.target_uniprot} | {expected_str} | {found_str} | {flag} |"
        )
    lines.append("")

    if report.negative_hits:
        lines.append("## ⚠️ Negative-Control Hits (peripheral drugs ranking too high)")
        lines.append("")
        lines.append("| Target | UniProt | Compound | pKd | Percentile |")
        lines.append("|---|---|---|---|---|")
        for hit in report.negative_hits:
            lines.append(
                f"| {hit.target_gene} | {hit.target_uniprot} | {hit.compound_name} | "
                f"{hit.pkd:.2f} | {hit.percentile:.0%} |"
            )
        lines.append("")

    if polypharm is not None and not polypharm.empty:
        lines.append("## Polypharmacology Leaderboard (top compounds by hit count)")
        lines.append("")
        lines.append("| Compound | Targets hit (pKd > threshold) | Mean pKd across hits |")
        lines.append("|---|---|---|")
        for _, row in polypharm.head(20).iterrows():
            lines.append(
                f"| {row['compound_name']} | {int(row['n_hits'])} | "
                f"{row['mean_pkd_hits']:.2f} |"
            )
        lines.append("")

    lines.append("---")
    lines.append("")
    lines.append("Generated by `scripts/05_sanity_check.py`.")
    return "\n".join(lines)


def write_report(report: SanityReport, polypharm: pd.DataFrame | None, out_path: Path) -> None:
    """Render and write the markdown report."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(render_markdown(report, polypharm), encoding="utf-8")
    logger.info("Wrote sanity report to %s.", out_path)
