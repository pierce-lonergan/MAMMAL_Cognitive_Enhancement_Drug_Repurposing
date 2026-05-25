"""PrimeKG loader + cognition-relevant subgraph extractor (STUB).

Per the v2 research doc §3 Class C:
    Source: Chandak, Huang & Zitnik, Scientific Data 2023
    DOI: 10.1038/s41597-023-01960-3
    Dataverse: 10.7910/DVN/IXA7BM (Harvard)
    Scale: 129,375 nodes × 4,050,249 edges, 30 relation types, 10 node types

Install / download:
    1. Download kg.csv / nodes.csv / edges.csv from Harvard Dataverse
       (≈400 MB compressed; >2 GB uncompressed)
    2. Place at data/kg/primekg/{kg.csv, nodes.csv, edges.csv}
    3. Optional: convert to parquet for faster loading

Usage (planned):
    from mammal_repurposing.cluster_c.primekg import load_primekg, score_compound_paths
    kg = load_primekg()
    scores = score_compound_paths(kg, compound_chembl_ids=[...], target_uniprots=[...])
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from mammal_repurposing.config import DATA_DIR

if TYPE_CHECKING:
    import networkx as nx

logger = logging.getLogger(__name__)

PRIMEKG_DIR = DATA_DIR / "kg" / "primekg"


def load_primekg(path: Path | str = PRIMEKG_DIR) -> "nx.Graph":
    """Load PrimeKG into a NetworkX graph.

    Raises:
        FileNotFoundError if the KG files haven't been downloaded.
        ImportError if networkx isn't installed.
    """
    try:
        import networkx as nx  # noqa: PLC0415
    except ImportError as e:
        raise ImportError("networkx not installed. `pip install networkx`.") from e

    kg_path = Path(path) / "kg.csv"
    if not kg_path.exists():
        raise FileNotFoundError(
            f"PrimeKG not found at {path}. Download from "
            f"https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/IXA7BM"
        )

    logger.info("Loading PrimeKG from %s (this can take ~30s)...", kg_path)
    df = pd.read_csv(kg_path)
    g = nx.from_pandas_edgelist(
        df, source="x_id", target="y_id", edge_attr=["relation", "display_relation"],
    )
    logger.info("PrimeKG loaded: %d nodes, %d edges.",
                g.number_of_nodes(), g.number_of_edges())
    return g


def score_compound_paths(
    g,
    *,
    compound_chembl_ids: list[str],
    target_uniprots: list[str],
    relations: list[str] | None = None,
    max_path_len: int = 3,
) -> pd.DataFrame:
    """Per-compound count of paths to any panel target via cognition-relevant relations.

    STUB — wire up after PrimeKG is downloaded. Suggested implementation: use
    Personalized PageRank from each compound restricted to the cognition subgraph
    (anchor at panel targets) for a fast continuous score.
    """
    raise NotImplementedError(
        "PrimeKG path scoring is a stub. Download PrimeKG first and implement "
        "Personalized PageRank (networkx.pagerank) or a Katz-style decayed walk."
    )
