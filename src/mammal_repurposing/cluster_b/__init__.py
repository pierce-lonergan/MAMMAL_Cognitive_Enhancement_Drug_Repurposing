"""Cluster B: ADMET / safety predictions via ADMET-AI.

Per the v2 research:
    "ADMET-AI has the best average rank among all models that have been
     evaluated on the 22 datasets in the TDC ADMET Leaderboard"
     (Swanson, Walters & Zou, Bioinformatics 2024, btae416).

The cluster outputs 41 endpoints per SMILES and feeds two downstream layers:
    1. Hard gates (gates/admet_gates.py) — physical kill criteria
    2. ADMET_score — soft signal for the fusion layer

Cluster B runs on CPU. It does NOT contend with the GPU for MAMMAL / Boltz-2 /
ESM2 — schedule it in parallel.
"""
