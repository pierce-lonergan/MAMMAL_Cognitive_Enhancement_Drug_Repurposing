"""Cluster C — Mechanism / knowledge graph evidence.

Components (per v2 research doc §3 Class C):
    primekg.py            — PrimeKG loader + subgraph extraction
    txgnn.py              — TxGNN zero-shot indication scoring
    cognition_anchor.py   — Virtual cognition phenotype anchor (union of 2-hop
                            neighborhoods around MCI/AD/ADHD/FXS/narcolepsy)

PrimeKG: Chandak, Huang & Zitnik, Scientific Data 2023 (DOI 10.1038/s41597-023-01960-3).
         129,375 nodes × 4,050,249 edges across 30 relation types.
         Harvard Dataverse DOI 10.7910/DVN/IXA7BM. Multi-GB download.

TxGNN:   Huang, Chandak, Wang et al., Nature Medicine 2024 (DOI 10.1038/s41591-024-03233-x).
         Foundation model for clinician-centered drug repurposing.
         Reports +49.2% indication / +35.1% contraindication vs baselines, zero-shot.
"""
