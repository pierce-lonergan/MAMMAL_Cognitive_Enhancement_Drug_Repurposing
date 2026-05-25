"""Local ChEMBL SQLite mirror — replaces REST API for ground-truth lookups.

The v1 REST backstop (`fetchers/chembl_groundtruth.py`) takes ~60s per pair
because of round-trip latency and intermittent 5xx retries — 78 hours wall-clock
for the full 4,713-pair Phase 1.2 sweep. The local SQLite mirror does the same
work as a SQL JOIN in milliseconds.

Built on `chembl-downloader` (0.5.2+), which handles version-pinned download
of the official `chembl_<release>_sqlite.tar.gz` (~4 GB) and extraction to a
local `.db` file at `~/.data/chembl/`.

Critical gotcha: 30%+ of ChEMBL bioactivity records are tagged to SALT FORMS,
not the parent compound. The `molecule_dictionary.parent_molregno` column
resolves this — we always query by parent. Without this normalization the
audit silently under-counts.

Quality filters applied (cap noise from low-confidence assays):
- `assays.assay_type = 'B'` (Binding assays, not Functional or ADMET)
- `activities.standard_type IN ('Ki', 'IC50', 'Kd', 'EC50')`
- `assays.confidence_score >= 7` (5-9 scale; 7+ = direct single protein target)
"""

from __future__ import annotations

import logging
import sqlite3
from functools import lru_cache
from pathlib import Path
from typing import Literal, Optional, TypedDict

import pandas as pd

logger = logging.getLogger(__name__)

# Singleton connection — opening the 12 GB DB takes a few hundred ms and we
# want every script in the project to share one handle.
_CONN: Optional[sqlite3.Connection] = None
_DB_PATH: Optional[Path] = None

EvidenceLabel = Literal["CORROBORATED", "AMBIGUOUS", "NOVEL", "CONTRADICTED"]


class PairEvidence(TypedDict):
    target_uniprot: str
    smiles: str
    inchikey: Optional[str]
    status: EvidenceLabel
    n_records: int
    best_pchembl: Optional[float]
    best_activity_type: Optional[str]
    best_standard_value_nm: Optional[float]


def get_conn() -> sqlite3.Connection:
    """Lazy-open the local ChEMBL SQLite mirror via chembl-downloader.

    Note on chembl-downloader 0.5.2 API quirks:
        - `latest()` returns a version string like "36" (NOT a path).
        - `download_extract_sqlite()` returns a `Path` to the .db file and is
          idempotent (no re-download if cached). We use it for path resolution.
        - `connect()` is a context manager; not useful for a singleton pattern.
    """
    global _CONN, _DB_PATH
    if _CONN is not None:
        return _CONN

    import chembl_downloader  # noqa: PLC0415

    try:
        _DB_PATH = Path(chembl_downloader.download_extract_sqlite())
    except Exception as e:
        raise FileNotFoundError(
            "ChEMBL SQLite not available. Run:\n"
            "  python -c \"import chembl_downloader; "
            "chembl_downloader.download_extract_sqlite()\"\n"
            "(~4 GB download + 12 GB extracted; ~15 min on a fast connection.)"
        ) from e
    if not _DB_PATH.exists():
        raise FileNotFoundError(f"ChEMBL SQLite expected at {_DB_PATH} but missing.")

    logger.info("Opening ChEMBL SQLite at %s", _DB_PATH)
    _CONN = sqlite3.connect(f"file:{_DB_PATH}?mode=ro", uri=True, check_same_thread=False)
    _CONN.row_factory = sqlite3.Row
    return _CONN


def db_path() -> Path:
    """Return the SQLite file path (after lazy open)."""
    get_conn()
    assert _DB_PATH is not None
    return _DB_PATH


@lru_cache(maxsize=1)
def chembl_release() -> str:
    """Return the ChEMBL release number from the DB path filename."""
    p = db_path()
    # filename is like chembl_35.db
    parts = p.stem.split("_")
    return parts[-1] if parts and parts[-1].isdigit() else "unknown"


# --- SMILES → InChIKey (Python-side; ChemicaLite is optional) ----------------

def smiles_to_inchikey(smiles: str) -> Optional[str]:
    """Compute InChIKey from SMILES via RDKit. Returns None on parse failure."""
    try:
        from rdkit import Chem  # noqa: PLC0415
        from rdkit import RDLogger  # noqa: PLC0415
        RDLogger.DisableLog("rdApp.*")  # silence rdkit's warning spam
    except ImportError as e:
        raise ImportError("rdkit required for SMILES->InChIKey. `pip install rdkit`.") from e
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToInchiKey(mol)


# --- Lookup queries ----------------------------------------------------------

# ChEMBL 36 schema (also verified for 33-36):
#   activities.molregno → molecule_dictionary.molregno (the salt form)
#   molecule_hierarchy.molregno → molecule_hierarchy.parent_molregno
#                                 (NULL row if the molecule has no parent record)
#   molecule_dictionary.molregno ↔ compound_structures.molregno → standard_inchi_key
#   activities.assay_id → assays.assay_id → assays.tid → target_dictionary.tid
#   target_dictionary.tid → target_components.tid → component_sequences.accession (UniProt)
#
# IMPORTANT: parent_molregno lives in molecule_hierarchy (a SEPARATE table),
# NOT as a column on molecule_dictionary. Earlier ChEMBL docs sometimes show
# it on molecule_dictionary; that's a different (older / aggregated) schema.


_LOOKUP_BY_INCHIKEY_SQL = """
WITH canonical AS (
    -- Resolve InChIKey → set of canonical parent molregnos. Use molecule_hierarchy
    -- to walk salt → parent; fall back to the row itself if no hierarchy entry.
    SELECT DISTINCT
        COALESCE(mh.parent_molregno, md.molregno) AS parent_molregno
    FROM compound_structures cs
    JOIN molecule_dictionary md ON cs.molregno = md.molregno
    LEFT JOIN molecule_hierarchy mh ON md.molregno = mh.molregno
    WHERE cs.standard_inchi_key = ?
),
expanded AS (
    -- All molregnos that resolve to one of these canonical parents —
    -- includes the parent itself and every salt form.
    SELECT md.molregno
    FROM molecule_dictionary md
    LEFT JOIN molecule_hierarchy mh ON md.molregno = mh.molregno
    WHERE COALESCE(mh.parent_molregno, md.molregno) IN (SELECT parent_molregno FROM canonical)
)
SELECT
    md.chembl_id           AS molecule_chembl_id,
    td.chembl_id           AS target_chembl_id,
    td.pref_name           AS target_pref_name,
    a.standard_type        AS standard_type,
    a.standard_value       AS standard_value,
    a.standard_units       AS standard_units,
    a.pchembl_value        AS pchembl_value,
    a.activity_comment     AS activity_comment,
    s.confidence_score     AS confidence_score,
    s.assay_type           AS assay_type,
    s.description          AS assay_description
FROM activities a
JOIN assays s                ON a.assay_id      = s.assay_id
JOIN target_dictionary td    ON s.tid           = td.tid
JOIN target_components tc    ON td.tid          = tc.tid
JOIN component_sequences cseq ON tc.component_id = cseq.component_id
JOIN molecule_dictionary md  ON a.molregno      = md.molregno
WHERE cseq.accession = ?
  AND md.molregno IN (SELECT molregno FROM expanded)
  AND s.assay_type = 'B'
  AND a.standard_type IN ('Ki', 'IC50', 'Kd', 'EC50')
  AND s.confidence_score >= 7
ORDER BY a.pchembl_value DESC NULLS LAST
"""


def lookup_activity(
    target_uniprot: str,
    *,
    smiles: Optional[str] = None,
    inchikey: Optional[str] = None,
) -> pd.DataFrame:
    """Return all bioactivity records for a (target, compound) pair.

    Either ``smiles`` or ``inchikey`` must be provided (smiles is converted
    via RDKit). Handles the parent_molregno salt-form gotcha by expanding to
    all molregnos sharing the same canonical parent.
    """
    if not (smiles or inchikey):
        raise ValueError("provide smiles or inchikey")

    if inchikey is None:
        inchikey = smiles_to_inchikey(smiles)
        if inchikey is None:
            return pd.DataFrame()  # unparseable SMILES

    conn = get_conn()
    return pd.read_sql_query(_LOOKUP_BY_INCHIKEY_SQL, conn,
                              params=(inchikey, target_uniprot))


def lookup_pair_evidence(
    target_uniprot: str,
    smiles: str,
) -> PairEvidence:
    """Return CORROBORATED / AMBIGUOUS / NOVEL / CONTRADICTED for a pair.

    pchembl thresholds (Ki/IC50/Kd in M, pchembl = -log10(M)):
        >= 6.0  → ≤ 1 µM     → CORROBORATED
        5.0-6.0 →   1-10 µM  → AMBIGUOUS
        < 5.0  → > 10 µM     → CONTRADICTED
        no records           → NOVEL
    """
    inchikey = smiles_to_inchikey(smiles)
    df = lookup_activity(target_uniprot, smiles=smiles, inchikey=inchikey)
    if df.empty:
        return PairEvidence(
            target_uniprot=target_uniprot, smiles=smiles, inchikey=inchikey,
            status="NOVEL", n_records=0, best_pchembl=None,
            best_activity_type=None, best_standard_value_nm=None,
        )

    n_records = len(df)
    pchembls = df["pchembl_value"].dropna()
    if pchembls.empty:
        # We have records but no harmonized pchembl — too noisy to classify
        return PairEvidence(
            target_uniprot=target_uniprot, smiles=smiles, inchikey=inchikey,
            status="NOVEL", n_records=n_records, best_pchembl=None,
            best_activity_type=None, best_standard_value_nm=None,
        )

    best_idx = pchembls.idxmax()
    best = float(pchembls.loc[best_idx])
    best_row = df.loc[best_idx]

    if best >= 6.0:
        status: EvidenceLabel = "CORROBORATED"
    elif best < 5.0:
        status = "CONTRADICTED"
    else:
        status = "AMBIGUOUS"

    return PairEvidence(
        target_uniprot=target_uniprot,
        smiles=smiles,
        inchikey=inchikey,
        status=status,
        n_records=n_records,
        best_pchembl=best,
        best_activity_type=str(best_row.get("standard_type") or ""),
        best_standard_value_nm=(
            float(best_row["standard_value"])
            if pd.notna(best_row.get("standard_value")) else None
        ),
    )


# --- Target ID audit helpers (replaces the REST audit in scripts/20) --------

def all_chembl_targets_for_uniprot(uniprot: str) -> pd.DataFrame:
    """Return every ChEMBL target whose components include this UniProt accession,
    with activity counts via a sub-query (no remote API).
    """
    sql = """
    SELECT
        td.chembl_id            AS target_chembl_id,
        td.pref_name            AS pref_name,
        td.target_type          AS target_type,
        td.organism             AS organism,
        (SELECT COUNT(*) FROM activities a
            JOIN assays s ON a.assay_id = s.assay_id
            WHERE s.tid = td.tid)  AS n_activities
    FROM target_dictionary td
    JOIN target_components tc ON td.tid = tc.tid
    JOIN component_sequences cseq ON tc.component_id = cseq.component_id
    WHERE cseq.accession = ?
    ORDER BY n_activities DESC
    """
    conn = get_conn()
    return pd.read_sql_query(sql, conn, params=(uniprot,))


def per_target_pchembl_records(target_uniprot: str) -> pd.DataFrame:
    """All distinct (molecule_chembl_id, inchikey, best_pchembl) tuples for a
    target, used by Phase 3.1 calibration to compute Spearman ρ vs MAMMAL/Boltz.
    """
    sql = """
    SELECT
        md.chembl_id                                AS molecule_chembl_id,
        cs.standard_inchi_key                       AS inchikey,
        MAX(a.pchembl_value)                        AS best_pchembl,
        COUNT(*)                                    AS n_records
    FROM activities a
    JOIN assays s                ON a.assay_id = s.assay_id
    JOIN target_dictionary td    ON s.tid = td.tid
    JOIN target_components tc    ON td.tid = tc.tid
    JOIN component_sequences cseq ON tc.component_id = cseq.component_id
    JOIN molecule_dictionary md  ON a.molregno = md.molregno
    LEFT JOIN molecule_hierarchy mh ON md.molregno = mh.molregno
    JOIN compound_structures cs  ON COALESCE(mh.parent_molregno, md.molregno) = cs.molregno
    WHERE cseq.accession = ?
      AND s.assay_type = 'B'
      AND a.standard_type IN ('Ki', 'IC50', 'Kd', 'EC50')
      AND s.confidence_score >= 7
      AND a.pchembl_value IS NOT NULL
    GROUP BY md.chembl_id, cs.standard_inchi_key
    HAVING MAX(a.pchembl_value) IS NOT NULL
    """
    conn = get_conn()
    return pd.read_sql_query(sql, conn, params=(target_uniprot,))


def chembl_actives_with_smiles_for_target(
    target_uniprot: str,
    min_pchembl: float = 6.0,
) -> pd.DataFrame:
    """All compounds at this target with pchembl ≥ threshold, returned with
    canonical SMILES from compound_structures.canonical_smiles.

    Cols: molecule_chembl_id, inchikey, canonical_smiles, best_pchembl, n_records.
    Used by diagnostics A (scaffold saturation) and D (Tanimoto to known actives).
    """
    sql = """
    SELECT
        md.chembl_id              AS molecule_chembl_id,
        cs.standard_inchi_key     AS inchikey,
        cs.canonical_smiles       AS canonical_smiles,
        MAX(a.pchembl_value)      AS best_pchembl,
        COUNT(*)                  AS n_records
    FROM activities a
    JOIN assays s                ON a.assay_id = s.assay_id
    JOIN target_dictionary td    ON s.tid = td.tid
    JOIN target_components tc    ON td.tid = tc.tid
    JOIN component_sequences cseq ON tc.component_id = cseq.component_id
    JOIN molecule_dictionary md  ON a.molregno = md.molregno
    LEFT JOIN molecule_hierarchy mh ON md.molregno = mh.molregno
    JOIN compound_structures cs  ON COALESCE(mh.parent_molregno, md.molregno) = cs.molregno
    WHERE cseq.accession = ?
      AND s.assay_type = 'B'
      AND a.standard_type IN ('Ki', 'IC50', 'Kd', 'EC50')
      AND s.confidence_score >= 7
      AND a.pchembl_value IS NOT NULL
      AND a.pchembl_value >= ?
      AND cs.canonical_smiles IS NOT NULL
    GROUP BY md.chembl_id, cs.standard_inchi_key, cs.canonical_smiles
    HAVING MAX(a.pchembl_value) IS NOT NULL
    """
    conn = get_conn()
    return pd.read_sql_query(sql, conn, params=(target_uniprot, min_pchembl))


def close_conn() -> None:
    """Close the singleton connection (tests, teardown)."""
    global _CONN
    if _CONN is not None:
        _CONN.close()
        _CONN = None
