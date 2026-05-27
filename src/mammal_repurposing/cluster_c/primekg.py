"""PrimeKG loader + cognition-relevant path scoring (real implementation).

Source: Chandak, Huang & Zitnik, Scientific Data 2023, DOI 10.1038/s41597-023-01960-3
Dataverse: 10.7910/DVN/IXA7BM. Scale: 129,375 nodes × 4,050,249 edges.

igraph is preferred over networkx for path queries — ~10× faster on this
edge count.
"""

from __future__ import annotations

import logging
from functools import lru_cache
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from mammal_repurposing.config import DATA_DIR

if TYPE_CHECKING:
    import igraph as ig

logger = logging.getLogger(__name__)

PRIMEKG_DIR = DATA_DIR / "kg" / "primekg"
PRIMEKG_CSV = PRIMEKG_DIR / "kg.csv"

# Relations we traverse for compound→target path scoring (drug-relevant).
DRUG_RELEVANT_RELATIONS = {
    "drug_protein",
    "drug_drug",
    "protein_protein",
    "pathway_protein",
    "bioprocess_protein",
    "molfunc_protein",
    "indication",
    "contraindication",
}


@lru_cache(maxsize=1)
def load_primekg(path: Path | str = PRIMEKG_CSV) -> "ig.Graph":
    """Load PrimeKG into an igraph.Graph. Cached process-wide.

    PrimeKG kg.csv columns: relation, display_relation, x_id, x_type, x_name,
    x_source, y_id, y_type, y_name, y_source.
    """
    try:
        import igraph as ig  # noqa: PLC0415
    except ImportError as e:
        raise ImportError("igraph not installed. `pip install igraph`.") from e

    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(
            f"PrimeKG not found at {p}. Download from "
            "https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/IXA7BM "
            "or run scripts/_wsl2_download_primekg.sh"
        )

    logger.info("Loading PrimeKG from %s (this can take ~30s) ...", p)
    df = pd.read_csv(p, dtype=str, low_memory=False)
    logger.info("  raw rows: %d", len(df))

    # Build deduplicated node table from both endpoints.
    x = df[["x_id", "x_type", "x_name"]].rename(columns=lambda c: c[2:])
    y = df[["y_id", "y_type", "y_name"]].rename(columns=lambda c: c[2:])
    nodes = pd.concat([x, y], ignore_index=True).drop_duplicates("id").reset_index(drop=True)
    id_to_idx: dict[str, int] = {nid: i for i, nid in enumerate(nodes["id"])}
    logger.info("  unique nodes: %d", len(nodes))

    edges = [(id_to_idx[x_id], id_to_idx[y_id])
             for x_id, y_id in zip(df["x_id"], df["y_id"])]
    g = ig.Graph(n=len(nodes), edges=edges, directed=False)
    g.vs["id"] = nodes["id"].tolist()
    g.vs["type"] = nodes["type"].tolist()
    g.vs["name"] = nodes["name"].tolist()
    g.es["relation"] = df["relation"].tolist()
    g.es["display_relation"] = df["display_relation"].tolist()
    logger.info("PrimeKG loaded: %d nodes, %d edges.",
                g.vcount(), g.ecount())
    return g


def resolve_uniprot_to_node(g: "ig.Graph", uniprot: str,
                            gene_symbol: str | None = None) -> int | None:
    """Find the PrimeKG gene/protein node for a UniProt accession.

    PrimeKG gene/protein nodes are keyed by NCBI gene ID with `name` = canonical
    gene symbol (e.g. id='1813', name='DRD2'). The cheapest resolver path is
    therefore by gene symbol against the `name` attribute.
    """
    # Cheap path — gene symbol against name
    if gene_symbol:
        matches = g.vs.select(name_eq=gene_symbol, type_eq="gene/protein")
        if len(matches) > 0:
            return matches[0].index
    # Exact id match (rare — only when caller already has NCBI gene id)
    matches = g.vs.select(id_eq=uniprot)
    if len(matches) > 0:
        return matches[0].index
    return None


def resolve_compound_to_node(g: "ig.Graph", chembl_id: str | None,
                             drugbank_id: str | None = None,
                             name: str | None = None) -> int | None:
    """Find the PrimeKG drug node for a compound by DrugBank ID, ChEMBL ID,
    or lowercased name. PrimeKG drug nodes are keyed by `id` = DrugBank ID
    with `name` = lowercased drug name (e.g. 'donepezil')."""
    if drugbank_id:
        matches = g.vs.select(id_eq=drugbank_id)
        if len(matches) > 0:
            return matches[0].index
    if chembl_id:
        for prefix in (chembl_id, f"CHEMBL.{chembl_id}", f"CHEMBL:{chembl_id}"):
            matches = g.vs.select(id_eq=prefix)
            if len(matches) > 0:
                return matches[0].index
    if name:
        # PrimeKG drug nodes store the lowercased drug name in `name` attr.
        # First exact match (cheap, igraph indexes name).
        matches = g.vs.select(name_eq=name.lower().strip(), type_eq="drug")
        if len(matches) > 0:
            return matches[0].index
        # Then try title-case (PrimeKG x_name "Donepezil")
        matches = g.vs.select(name_eq=name.strip().capitalize(), type_eq="drug")
        if len(matches) > 0:
            return matches[0].index
    return None


def score_compound_paths_ppr(
    g: "ig.Graph",
    compound_node: int,
    target_nodes: list[int],
    damping: float = 0.85,
) -> dict[int, float]:
    """Personalized PageRank from one compound; return PPR mass per target node.

    Faster + better-resolution than enumerating all simple paths.
    """
    n = g.vcount()
    reset = [0.0] * n
    reset[compound_node] = 1.0
    ppr = g.personalized_pagerank(damping=damping, reset=reset)
    return {t: float(ppr[t]) for t in target_nodes}


def shortest_path_length(g: "ig.Graph", src: int, dst: int) -> int:
    """Shortest-path edges between two nodes; -1 if disconnected."""
    d = g.shortest_paths_dijkstra(source=src, target=dst)[0][0]
    return int(d) if d != float("inf") else -1


def score_compound_against_panel(
    g: "ig.Graph",
    *,
    compound_chembl_id: str | None,
    compound_drugbank_id: str | None,
    target_uniprots: list[str],
    compound_name: str | None = None,
    target_gene_symbols: list[str] | None = None,
) -> dict:
    """End-to-end: resolve compound + panel targets to nodes, compute PPR + shortest path.

    target_gene_symbols (if provided, aligned with target_uniprots) enables
    the cheap name-based protein resolver against PrimeKG's gene/protein
    nodes (`name` = canonical gene symbol).
    """
    src = resolve_compound_to_node(g, compound_chembl_id, compound_drugbank_id, compound_name)
    if src is None:
        return {
            "compound_node_found": False,
            "target_nodes_found": 0,
            "ppr_sum": 0.0,
            "shortest_path_min": -1,
            "n_targets_reachable": 0,
        }

    target_indices: list[int] = []
    gene_map = (dict(zip(target_uniprots, target_gene_symbols))
                if target_gene_symbols else {})
    for u in target_uniprots:
        n = resolve_uniprot_to_node(g, u, gene_symbol=gene_map.get(u))
        if n is not None:
            target_indices.append(n)

    if not target_indices:
        return {
            "compound_node_found": True,
            "target_nodes_found": 0,
            "ppr_sum": 0.0,
            "shortest_path_min": -1,
            "n_targets_reachable": 0,
        }

    ppr = score_compound_paths_ppr(g, src, target_indices)
    ppr_sum = float(sum(ppr.values()))
    dists = [shortest_path_length(g, src, t) for t in target_indices]
    reachable = [d for d in dists if d > 0]

    return {
        "compound_node_found": True,
        "target_nodes_found": len(target_indices),
        "ppr_sum": ppr_sum,
        "shortest_path_min": min(reachable) if reachable else -1,
        "n_targets_reachable": len(reachable),
    }
