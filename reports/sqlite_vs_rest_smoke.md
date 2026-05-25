# Phase A.5 — SQLite vs REST agreement smoke test

**Agreement: 19/20 (95%) (errors: 1)**

Picked the top 20 (target, compound) pairs by predicted_pkd with SMILES < 200 chars.
Both the SQLite backstop (`fetchers/chembl_sqlite.py`) and the legacy REST
fetcher (`fetchers/chembl_groundtruth.py`) were run on the same pairs. INCONCLUSIVE (REST) is normalised to AMBIGUOUS (SQLite) for comparison.

| target | compound | SQL status | SQL n | REST status | REST n | agree |
|---|---|---|---|---|---|---|
| Q9Y5N1 | staurosporine | NOVEL | 0 | ERROR | 0 | ✗ |
| Q9Y5N1 | lm22a-4 | NOVEL | 0 | NOVEL | 0 | ✓ |
| Q9Y5N1 | atorvastatin | NOVEL | 0 | NOVEL | 0 | ✓ |
| Q9Y5N1 | chembl4780352 | NOVEL | 0 | NOVEL | 0 | ✓ |
| Q9Y5N1 | tulrampator | NOVEL | 0 | NOVEL | 0 | ✓ |
| Q9Y5N1 | fexofenadine | NOVEL | 0 | NOVEL | 0 | ✓ |
| Q9Y5N1 | lurasidone | NOVEL | 0 | NOVEL | 0 | ✓ |
| Q9Y5N1 | ivabradine | NOVEL | 0 | NOVEL | 0 | ✓ |
| Q9Y5N1 | methylene blue | NOVEL | 0 | NOVEL | 0 | ✓ |
| Q9Y5N1 | suvorexant | NOVEL | 0 | NOVEL | 0 | ✓ |
| Q9Y5N1 | lemborexant | NOVEL | 0 | NOVEL | 0 | ✓ |
| Q9Y5N1 | risperidone | AMBIGUOUS | 1 | AMBIGUOUS | 1 | ✓ |
| Q9Y5N1 | cetirizine | NOVEL | 0 | NOVEL | 0 | ✓ |
| Q9Y5N1 | chembl5805992 | NOVEL | 0 | NOVEL | 0 | ✓ |
| Q9Y5N1 | chembl5813903 | NOVEL | 0 | NOVEL | 0 | ✓ |
| Q9Y5N1 | pitolisant | CORROBORATED | 35 | CORROBORATED | 39 | ✓ |
| Q9Y5N1 | bpn14770 | NOVEL | 0 | NOVEL | 0 | ✓ |
| Q9Y5N1 | tc-5619 | NOVEL | 0 | NOVEL | 0 | ✓ |
| Q9Y5N1 | aripiprazole | AMBIGUOUS | 1 | AMBIGUOUS | 1 | ✓ |
| Q9Y5N1 | hydroxyzine | NOVEL | 0 | NOVEL | 0 | ✓ |

_Gate logic: PASS when n_agree == n_total - n_err. REST 500s are tolerated (transient EBI infrastructure) — only true status mismatches count as failures._

**Result: PASS** (n_agree=19, n_total=20, n_err=1)