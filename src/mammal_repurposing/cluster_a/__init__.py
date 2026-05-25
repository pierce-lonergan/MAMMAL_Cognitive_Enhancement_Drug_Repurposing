"""Cluster A — Structure / pocket-aware affinity.

Components (per v2 research doc §3 Class A):
    esm2_embed.py     — ESM2-650M target embeddings (cached, one-time)
    boltz_runner.py   — Boltz-2 structure prediction (heavy, one-time per target)
    boltzina.py       — Boltzina affinity-only mode (Furui & Ohue 2025)

VRAM sequencing: never load Boltz-2 and ESM2 simultaneously. Use
`torch.cuda.empty_cache()` + explicit `del` between phases. Affinity-only
Boltzina mode fits ~7-8 GB on the RTX 5070 per Nebius's L40S benchmark.
"""
