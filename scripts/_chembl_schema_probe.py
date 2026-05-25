"""Quick ChEMBL 36 schema introspection — find the right table for parent compound linkage."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.fetchers.chembl_sqlite import get_conn

conn = get_conn()

print("=== molecule_dictionary columns ===")
for row in conn.execute("PRAGMA table_info(molecule_dictionary)").fetchall():
    print(f"  {row[1]} ({row[2]})")

print("\n=== molecule_* tables ===")
for row in conn.execute(
    "SELECT name FROM sqlite_master WHERE type='table' AND name LIKE 'molecule%'"
).fetchall():
    print(f"  {row[0]}")

print("\n=== molecule_hierarchy columns ===")
for row in conn.execute("PRAGMA table_info(molecule_hierarchy)").fetchall():
    print(f"  {row[1]} ({row[2]})")

# Sample a parent relationship
print("\n=== molecule_hierarchy sample (where molregno != parent_molregno) ===")
sql = """
SELECT mh.molregno, mh.parent_molregno, mh.active_molregno,
       md.chembl_id AS molecule_id,
       (SELECT chembl_id FROM molecule_dictionary WHERE molregno = mh.parent_molregno) AS parent_id
FROM molecule_hierarchy mh
JOIN molecule_dictionary md ON mh.molregno = md.molregno
WHERE mh.molregno != mh.parent_molregno
LIMIT 5
"""
for row in conn.execute(sql).fetchall():
    print(f"  {dict(row)}")

print("\n=== sample lookup: pitolisant (CHEMBL2106101 expected) ===")
sql2 = """
SELECT md.chembl_id, md.molregno, mh.parent_molregno,
       cs.standard_inchi_key
FROM molecule_dictionary md
JOIN compound_structures cs ON md.molregno = cs.molregno
LEFT JOIN molecule_hierarchy mh ON md.molregno = mh.molregno
WHERE md.chembl_id LIKE 'CHEMBL2106101'
LIMIT 5
"""
for row in conn.execute(sql2).fetchall():
    print(f"  {dict(row)}")
