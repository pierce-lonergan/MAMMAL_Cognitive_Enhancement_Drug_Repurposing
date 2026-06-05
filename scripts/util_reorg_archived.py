"""One-off housekeeping: reorganize research/4-tier/archived for open-source.

Sanitizes broken filenames (spaces, '+', ',', em-dash) into clean kebab-case and
groups the 24 archived research docs into thematic subdirs. Already-clean names
(Title-Case-Kebab) are kept verbatim, only relocated (per the "only fix broken
names" preference). Uses `git mv` to preserve history.

Usage: python scripts/util_reorg_archived.py
"""
from __future__ import annotations

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = "research/4-tier/archived"

# old basename -> new relative path under BASE
MOVES = {
    # methodology / formal research
    "Clinical Effect-Size Translation Function.md":
        "methodology/clinical-effect-size-translation.md",
    "Clinical Effect-Size Translation Function A Methodology Pre-Registration for Bayesian Cognition-Enhancement Drug Repurposing.md":
        "methodology/clinical-effect-size-translation-prereg.md",
    "Cluster D Methodology Report — Gate 3 (Held-out Cognition GWAS L2G) and Gate 2 (Multi-Modulator Curation).md":
        "methodology/cluster-d-methodology-gate3-gate2.md",
    "Multi-Source Neurobiological Prior for Cognition Target Prioritization.md":
        "methodology/multi-source-neurobiological-prior.md",
    "Perturbational Evidence Axis.md":
        "methodology/perturbational-evidence-axis.md",
    "MH1 + MH2 Meta-Analytic Prior Expansion for V7 CPT Bayesian Pharmacology Pipeline.md":
        "methodology/mh1-mh2-meta-analytic-prior-expansion.md",
    "MH1+MH2_ PRISMA, Anchor Expansion.md":
        "methodology/mh1-mh2-prisma-anchor-expansion.md",
    "MH3_per_cell_line_random_effect_deep_research.md":
        "methodology/mh3-per-cell-line-random-effect.md",
    "MH8 Methods Clarity Research.md":
        "methodology/mh8-methods-clarity.md",
    "Hybrid Architecture for MAMMAL-Based Cognitive-Enhancement Drug Repurposing.md":
        "methodology/hybrid-architecture-mammal.md",
    "MAMMAL for Healthy Cognitive Enhancement Drug Repurposing.md":
        "methodology/mammal-for-cognitive-enhancement.md",
    "Multi Head DTI.md":
        "methodology/multi-head-dti.md",
    # deep dives
    "Deep Dive_ CUDA Equivariance PyTorch.md":
        "deep-dives/cuda-equivariance-pytorch.md",
    "Deep Dive_ LINCS L1000 chemCPA Training.md":
        "deep-dives/lincs-l1000-chemcpa-training.md",
    "Technical Feasibility Deep-Dive Adding a Phenotypic.md":
        "deep-dives/phenotypic-axis-feasibility.md",
    "Diagnosing MAMMAL DTI Anti-Correlation.md":
        "deep-dives/diagnosing-mammal-dti-anticorrelation.md",
    "ChEMBL Database Optimization Deep Dive.md":
        "deep-dives/chembl-database-optimization.md",
    # analysis notes (clean names kept verbatim; broken ones sanitized)
    "Cognition-44Target-Liability-Panel.md":
        "analysis-notes/Cognition-44Target-Liability-Panel.md",
    "Isotonic-PerTarget-Calibration.md":
        "analysis-notes/Isotonic-PerTarget-Calibration.md",
    "Graczyk-Selectivity-Faceted-Shortlist.md":
        "analysis-notes/Graczyk-Selectivity-Faceted-Shortlist.md",
    "Pocket-Conditioned-Boltz2.md":
        "analysis-notes/Pocket-Conditioned-Boltz2.md",
    "Boltzina + MAMMAL Fine-tune.md":
        "analysis-notes/boltzina-mammal-finetune.md",
    # infrastructure / build notes
    "PyG Windows Nightly WSL2 Fixes.md":
        "infrastructure/pyg-windows-wsl2-fixes.md",
    "tiledbsoma Windows Build Blocker.md":
        "infrastructure/tiledbsoma-windows-build-blocker.md",
}


def main() -> int:
    moved = 0
    for old, new in MOVES.items():
        src = ROOT / BASE / old
        dst = ROOT / BASE / new
        if not src.exists():
            print(f"  SKIP (missing): {old}")
            continue
        dst.parent.mkdir(parents=True, exist_ok=True)
        r = subprocess.run(["git", "mv", str(src), str(dst)], cwd=ROOT,
                           capture_output=True, text=True)
        if r.returncode != 0:
            print(f"  ERR git mv {old}: {r.stderr.strip()}")
        else:
            moved += 1
    print(f"moved {moved}/{len(MOVES)}")
    # update the 2 provenance comments that name an archived doc
    for rel, old_ref, new_ref in [
        ("scripts/39_v5_pocket_conditional_liability.py",
         "research/4-tier/archived/", "research/4-tier/archived/analysis-notes/"),
        ("src/mammal_repurposing/gates/liability_panel.py",
         "research/4-tier/archived/Pocket-Conditioned-Boltz2.md",
         "research/4-tier/archived/analysis-notes/Pocket-Conditioned-Boltz2.md"),
    ]:
        p = ROOT / rel
        t = p.read_text(encoding="utf-8")
        if old_ref in t:
            p.write_text(t.replace(old_ref, new_ref), encoding="utf-8")
            print(f"  updated ref in {rel}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
