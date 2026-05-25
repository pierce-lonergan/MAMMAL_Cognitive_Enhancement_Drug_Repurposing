"""Virtual cognition phenotype anchor for KG-based repurposing.

Because "healthy cognitive enhancement" is not a disease node in PrimeKG,
we anchor the KG query to the UNION of 2-hop neighborhoods around these
disease nodes (per v2 research doc §3 Class C):

    - Mild cognitive impairment (EFO_0006816)
    - Alzheimer's disease (EFO_0000249) [mechanistic proxy; weight 0.8]
    - Attention-deficit/hyperactivity disorder (EFO_0003888)
    - Fragile X syndrome (EFO_0004247) [PDE4D anchor]
    - Narcolepsy (EFO_0003781) [HCRTR + HRH3 anchor]

Each compound gets four sub-scores per the doc's KG_score formula:
    indication        — mean TxGNN indication probability across the 5 anchors
    contraindication  — mean TxGNN contraindication probability (subtracted)
    path_count        — PrimeKG path count via {drug_protein, protein_protein,
                        pathway_protein, bioprocess_protein} to any panel target
    side_effect       — overlap with cognition-degrading side-effect terms
                        (sedation, somnolence, cognitive impairment, anticholinergic)

Composite: KG_score = w_ind*ind - w_con*con + w_path*log(1+paths) - w_se*se
Default weights from configs/weights.yaml: (0.4, 0.3, 0.2, 0.1).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml

from mammal_repurposing.config import PROJECT_ROOT


@dataclass
class CognitionAnchor:
    """The 5-disease anchor + side-effect penalty terms."""

    disease_efo_ids: list[str]
    disease_weights: dict[str, float]
    k_hops: int
    cognition_degrading_side_effects: list[str]

    @classmethod
    def from_config(cls, path: Path | str = PROJECT_ROOT / "configs" / "weights.yaml"):
        cfg = yaml.safe_load(open(path, encoding="utf-8"))["cognition_virtual_anchor"]
        ids = [d["id"] for d in cfg["disease_nodes"]]
        weights = {d["id"]: float(d["weight"]) for d in cfg["disease_nodes"]}
        return cls(
            disease_efo_ids=ids,
            disease_weights=weights,
            k_hops=int(cfg.get("k_hops", 2)),
            cognition_degrading_side_effects=cfg["cognition_degrading_side_effects"],
        )
