"""Allosteric Awareness Benchmark - scoring + analysis.

Loads the curated benchmark CSV (data/raw/allosteric_benchmark.csv), scores
each (target, compound) pair with MAMMAL's DTI head, and produces a benchmark
report quantifying:

    1. Does MAMMAL differentiate allosteric from orthosteric ligands at each
       target? (effect-size between binding-mode groups)
    2. Does predicted pKd correlate with measured activity within each binding
       mode? (Spearman rho per group)

This benchmark is the contribution: no published benchmark exists for "does a
DTI foundation model correctly handle allosteric vs orthosteric binding."
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
import pandas as pd


@dataclass
class BenchmarkGroup:
    binding_mode: str
    n: int
    pkd_mean: float
    pkd_std: float
    spearman_predicted_vs_measured: float | None


@dataclass
class BenchmarkTargetRow:
    target_uniprot: str
    target_gene: str
    groups: dict[str, BenchmarkGroup] = field(default_factory=dict)

    def gap(self, mode_a: str, mode_b: str) -> float | None:
        """Mean pKd gap between two binding-mode groups (a - b)."""
        ga = self.groups.get(mode_a)
        gb = self.groups.get(mode_b)
        if not ga or not gb:
            return None
        return ga.pkd_mean - gb.pkd_mean


def _spearman(x: np.ndarray, y: np.ndarray) -> float | None:
    if len(x) < 3:
        return None
    rx = pd.Series(x).rank().to_numpy()
    ry = pd.Series(y).rank().to_numpy()
    c = np.corrcoef(rx, ry)
    return float(c[0, 1]) if not math.isnan(c[0, 1]) else None


def analyze_benchmark(scored: pd.DataFrame) -> list[BenchmarkTargetRow]:
    """scored DataFrame columns: target_uniprot, target_gene, compound_name,
    binding_mode, measured_activity_nm, activity_type, predicted_pkd."""
    rows: list[BenchmarkTargetRow] = []
    for (uniprot, gene), grp in scored.groupby(["target_uniprot", "target_gene"]):
        target_row = BenchmarkTargetRow(target_uniprot=uniprot, target_gene=gene)
        for mode, mode_grp in grp.groupby("binding_mode"):
            preds = mode_grp["predicted_pkd"].to_numpy(dtype=float)
            measured = mode_grp.dropna(subset=["measured_activity_nm"])
            if len(measured) >= 3:
                measured_pkd = 9.0 - np.log10(measured["measured_activity_nm"].astype(float))
                rho = _spearman(
                    measured["predicted_pkd"].to_numpy(dtype=float),
                    measured_pkd.to_numpy(dtype=float),
                )
            else:
                rho = None
            target_row.groups[mode] = BenchmarkGroup(
                binding_mode=mode,
                n=len(mode_grp),
                pkd_mean=float(np.mean(preds)),
                pkd_std=float(np.std(preds)),
                spearman_predicted_vs_measured=rho,
            )
        rows.append(target_row)
    return rows


def render_markdown(rows: list[BenchmarkTargetRow]) -> str:
    lines: list[str] = []
    lines.append("# Allosteric Awareness Benchmark — MAMMAL DTI Head")
    lines.append("")
    lines.append("Benchmark dataset: `data/raw/allosteric_benchmark.csv`. "
                 "Each (target, compound) triple is labeled by binding mode "
                 "(orthosteric_agonist/antagonist/inhibitor, allosteric_pam/nam/partial_agonist) "
                 "with measured affinity from literature.")
    lines.append("")
    lines.append("Two questions per target:")
    lines.append("")
    lines.append("1. **Mode-discrimination**: does MAMMAL's predicted pKd differ between binding-mode groups? "
                 "A blind model produces similar pKd for orthosteric and allosteric ligands.")
    lines.append("2. **Within-group calibration**: within a single binding mode, does predicted pKd correlate "
                 "with measured affinity? (Spearman ρ.)")
    lines.append("")

    for row in rows:
        lines.append(f"## {row.target_gene} ({row.target_uniprot})")
        lines.append("")
        lines.append("| Binding mode | n | predicted pKd (mean ± std) | Spearman ρ (pred vs measured) |")
        lines.append("|---|---|---|---|")
        for mode_name in sorted(row.groups.keys()):
            g = row.groups[mode_name]
            rho_s = f"{g.spearman_predicted_vs_measured:+.2f}" if g.spearman_predicted_vs_measured is not None else "n/a (n<3)"
            lines.append(f"| `{mode_name}` | {g.n} | {g.pkd_mean:.2f} ± {g.pkd_std:.2f} | {rho_s} |")

        # Surface key contrasts if both modes are present
        modes = list(row.groups.keys())
        ortho = [m for m in modes if m.startswith("orthosteric")]
        allo = [m for m in modes if m.startswith("allosteric")]
        if ortho and allo:
            o_mean = np.mean([row.groups[m].pkd_mean for m in ortho])
            a_mean = np.mean([row.groups[m].pkd_mean for m in allo])
            gap = o_mean - a_mean
            verdict = "✅ comparable" if abs(gap) < 0.5 else (
                f"⚠️ orthosteric > allosteric by {gap:.2f} pKd" if gap > 0
                else f"⚠️ allosteric > orthosteric by {-gap:.2f} pKd"
            )
            lines.append("")
            lines.append(f"**Mode gap**: orthosteric mean {o_mean:.2f} vs allosteric mean {a_mean:.2f} → {verdict}")
        lines.append("")

    lines.append("---")
    lines.append("## Interpretation")
    lines.append("")
    lines.append("- **Mode gap close to zero (|<0.5|)** at a target = MAMMAL is allosteric-aware at that target. "
                 "Allosteric ligands ranked alongside orthosteric.")
    lines.append("- **Mode gap >1 pKd in favor of orthosteric** = MAMMAL underranks allosteric mechanism. "
                 "Downstream candidate ranking should down-weight allosteric-heavy targets or apply "
                 "target-specific recalibration. This is the failure mode the benchmark was built to detect.")
    lines.append("- **Within-group ρ low (<0.3)** = MAMMAL does not rank-order even within one binding mode at this target. "
                 "Strong signal that the training data for this target is sparse or noisy.")
    lines.append("")
    lines.append("This benchmark dataset and code are released as a methods contribution. "
                 "Cite as: MAMMAL Allosteric Awareness Benchmark, 2026.")
    lines.append("")
    lines.append("Generated by `scripts/12_allosteric_benchmark.py`.")
    return "\n".join(lines)
