"""V8 / Cluster E — πphen Perturbational Evidence Axis.

Target-agnostic phenotypic evidence axis parallel to V6.A (target-binding)
and V6.B (target-relevance). Per `research/4-tier/Perturbational Evidence
Axis.md` + `Technical Feasibility Deep-Dive Adding a Phenotypic.md`.

Subsystems:
  - ingest_lincs.py   — LINCS L1000 cmapPy WTCS index + cell-line metadata
  - ingest_jumpcp.py  — JUMP-CP cpg0016 S3 sync + pycytominer consensus
  - chemcpa_train.py  — chemCPA generative imputation training (V8.2)
  - mofa_embed.py     — MOFA+ joint K=30 embedding across 7 views (V8.3)
  - joint_phenotype.py — V7+V8 joint posterior PyMC (V8.5)

Target venue: Nature Machine Intelligence (A realistic) or Nature Methods
(A+ stretch at Gate 1 AMI ≥ 0.6).
"""

from __future__ import annotations
