"""Lateral 6.1 — binding-mode mix per target.

For each (target, compound) joined to ChEMBL, parse mechanism_of_action +
action_type to find out what TYPE of pharmacology dominates the library.

Special case for GRIN2A/GRIN2B:
  - If >50% of library compounds at GRIN2B are ifenprodil-class (i.e.,
    listed mechanism_of_action contains 'NEGATIVE ALLOSTERIC MODULATOR'
    or 'antagonist' at the ATD-allosteric site), the target's pharmacology
    is dominated by a DIMER-INTERFACE binding mode that MAMMAL's single-chain
    sequence input CANNOT see.
  - Action: deprecate GRIN2B from the MAMMAL panel; rely on Boltz-2 with
    explicit heterodimeric template (Karakas et al. Nature 2011).

For the SLC6 transporters:
  - Most compounds are 'INHIBITOR' at the substrate site. Predominantly
    one mode → if MAMMAL still inverts, the failure is rank-resolution
    within that mode, not mode-confusion.

Reference: ChEMBL action_type taxonomy (see ChEMBL schema docs);
  Karakas et al. 2011 Nature on GluN2B ATD ifenprodil site;
  Mony et al. 2012 Mol Pharmacol on the UL-UL dimer interface mapping.
"""

from __future__ import annotations

import logging
import sqlite3
from collections import Counter
from dataclasses import dataclass

import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class BindingModeMixResult:
    target_uniprot: str
    target_gene: str
    n_library_with_action_type: int
    top_action_type: str
    top_action_type_pct: float
    action_type_counts: dict[str, int]
    n_allosteric: int
    pct_allosteric: float
    verdict: str   # single_mode | mixed_modes | allosteric_dominant | unknown


_ACTION_TYPE_SQL = """
SELECT
    md.chembl_id            AS molecule_chembl_id,
    cs.standard_inchi_key   AS inchikey,
    mr.action_type          AS action_type,
    mech.mechanism_of_action AS mechanism_of_action,
    mech.binding_site_comment AS binding_site_comment
FROM molecule_dictionary md
LEFT JOIN molecule_hierarchy mh ON md.molregno = mh.molregno
JOIN compound_structures cs  ON COALESCE(mh.parent_molregno, md.molregno) = cs.molregno
LEFT JOIN drug_mechanism mech ON md.molregno = mech.molregno
LEFT JOIN action_type mr ON mech.action_type = mr.action_type
LEFT JOIN target_dictionary td ON mech.tid = td.tid
LEFT JOIN target_components tc ON td.tid = tc.tid
LEFT JOIN component_sequences cseq ON tc.component_id = cseq.component_id
WHERE cseq.accession = ? OR cseq.accession IS NULL
"""

# Simpler version: just look up action_type for any (molecule, target) row
_TARGETED_ACTION_TYPE_SQL = """
SELECT DISTINCT
    md.chembl_id            AS molecule_chembl_id,
    cs.standard_inchi_key   AS inchikey,
    mech.action_type        AS action_type,
    mech.mechanism_of_action AS mechanism_of_action,
    mech.binding_site_comment AS binding_site_comment
FROM drug_mechanism mech
JOIN target_dictionary td ON mech.tid = td.tid
JOIN target_components tc ON td.tid = tc.tid
JOIN component_sequences cseq ON tc.component_id = cseq.component_id
JOIN molecule_dictionary md ON mech.molregno = md.molregno
LEFT JOIN molecule_hierarchy mh ON md.molregno = mh.molregno
JOIN compound_structures cs ON COALESCE(mh.parent_molregno, md.molregno) = cs.molregno
WHERE cseq.accession = ?
"""


def fetch_action_types_for_target(
    conn: sqlite3.Connection,
    target_uniprot: str,
) -> pd.DataFrame:
    """Returns DataFrame: molecule_chembl_id, inchikey, action_type,
    mechanism_of_action, binding_site_comment."""
    return pd.read_sql_query(_TARGETED_ACTION_TYPE_SQL, conn, params=(target_uniprot,))


def _is_allosteric(row: pd.Series) -> bool:
    """Detect allosteric binding from mechanism_of_action / binding_site_comment / action_type."""
    haystacks = []
    for col in ("action_type", "mechanism_of_action", "binding_site_comment"):
        v = row.get(col)
        if isinstance(v, str):
            haystacks.append(v.lower())
    text = " ".join(haystacks)
    keywords = [
        "allosteric", "negative allosteric modulator", "positive allosteric modulator",
        "modulator", "ifenprodil", "atd", "amino-terminal",
    ]
    return any(k in text for k in keywords)


def diagnose(
    target_uniprot: str,
    target_gene: str,
    library_inchikeys: list[str],
    chembl_action_types: pd.DataFrame,
) -> BindingModeMixResult:
    """Join library compounds to ChEMBL action_type records; classify dominant mode."""
    sub = chembl_action_types[
        chembl_action_types["inchikey"].isin(library_inchikeys)
    ].copy()

    if sub.empty:
        return BindingModeMixResult(
            target_uniprot=target_uniprot, target_gene=target_gene,
            n_library_with_action_type=0,
            top_action_type="UNKNOWN", top_action_type_pct=float("nan"),
            action_type_counts={},
            n_allosteric=0, pct_allosteric=float("nan"),
            verdict="no_action_type_annotations",
        )

    sub["is_allosteric"] = sub.apply(_is_allosteric, axis=1)
    action_counts = Counter(sub["action_type"].dropna().astype(str))
    top_action, top_count = action_counts.most_common(1)[0] if action_counts else ("UNKNOWN", 0)
    top_pct = 100.0 * top_count / max(len(sub), 1)
    n_alo = int(sub["is_allosteric"].sum())
    pct_alo = 100.0 * n_alo / max(len(sub), 1)

    if pct_alo >= 50:
        verdict = "allosteric_dominant"
    elif top_pct >= 70:
        verdict = "single_mode_dominant"
    elif len(action_counts) >= 3:
        verdict = "mixed_modes"
    else:
        verdict = "two_modes_or_uncertain"

    return BindingModeMixResult(
        target_uniprot=target_uniprot, target_gene=target_gene,
        n_library_with_action_type=len(sub),
        top_action_type=top_action, top_action_type_pct=top_pct,
        action_type_counts=dict(action_counts),
        n_allosteric=n_alo, pct_allosteric=pct_alo,
        verdict=verdict,
    )
