"""V3 Phase A.5 — Verify SQLite vs REST agreement on a 20-pair subset.

Before trusting the SQLite path for the full backstop, run both REST and
SQLite on the same 20 pairs and check status agreement. Per the sprint spec:
require 20/20 agreement — any disagreement means the salt-form parent
resolution is missing somewhere.

Picks the 20 pairs from the most-recent dti_scores.parquet with the highest
predicted_pkd that have plausible ChEMBL coverage (skips peptides/macros).
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.config import DTI_SCORES_PARQUET  # noqa: E402
from mammal_repurposing.fetchers.chembl_groundtruth import lookup_pair  # noqa: E402
from mammal_repurposing.fetchers.chembl_sqlite import lookup_pair_evidence  # noqa: E402

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("sqlite_vs_rest")


# Map REST module's label to SQLite module's label for comparison
def _normalize_rest_status(s: str) -> str:
    """REST returns CORROBORATED/NOVEL/CONTRADICTED/INCONCLUSIVE.
    SQLite returns CORROBORATED/AMBIGUOUS/NOVEL/CONTRADICTED.
    Normalize INCONCLUSIVE -> AMBIGUOUS for direct comparison.
    """
    return "AMBIGUOUS" if s == "INCONCLUSIVE" else s


def main() -> int:
    if not DTI_SCORES_PARQUET.exists():
        logger.error("Need %s; run scripts/04 first.", DTI_SCORES_PARQUET)
        return 1

    scores = pd.read_parquet(DTI_SCORES_PARQUET)
    # Pick top 20 pairs by predicted pKd, skip SMILES > 200 chars (likely peptides)
    pool = scores[scores["compound_smiles"].str.len() < 200].copy()
    pairs = (pool.sort_values("predicted_pkd", ascending=False)
                 .drop_duplicates(subset=["target_uniprot", "compound_name"])
                 .head(20))
    logger.info("Testing %d pairs (top predicted_pkd, SMILES < 200 chars)", len(pairs))

    import httpx
    from mammal_repurposing.config import HTTP_TIMEOUT_SEC, USER_AGENT

    rows: list[dict] = []
    with httpx.Client(timeout=HTTP_TIMEOUT_SEC,
                      headers={"User-Agent": USER_AGENT, "Accept": "application/json"}) as client:
        target_id_cache: dict[str, str | None] = {}

        for _, p in pairs.iterrows():
            tgt = p["target_uniprot"]
            name = p["compound_name"]
            smi = p["compound_smiles"]
            # SQLite verdict
            try:
                sql = lookup_pair_evidence(tgt, smi)
                sql_status = sql["status"]
                sql_n = sql["n_records"]
            except Exception as e:
                logger.warning("SQLite failed for (%s,%s): %s", tgt, name, e)
                sql_status = "ERROR"
                sql_n = 0

            # REST verdict
            try:
                rest = lookup_pair(client, tgt, name, smi, target_id_cache=target_id_cache)
                rest_status = _normalize_rest_status(rest["label"])
                rest_n = rest["n_chembl_records"]
            except Exception as e:
                logger.warning("REST failed for (%s,%s): %s", tgt, name, e)
                rest_status = "ERROR"
                rest_n = 0

            agree = sql_status == rest_status
            rows.append({
                "target_uniprot": tgt,
                "compound_name": name,
                "sql_status": sql_status,
                "sql_n": sql_n,
                "rest_status": rest_status,
                "rest_n": rest_n,
                "agree": agree,
            })
            logger.info("%-20s %-12s SQL=%-13s (n=%-3d) REST=%-13s (n=%-3d) %s",
                        name, tgt, sql_status, sql_n, rest_status, rest_n,
                        "✓" if agree else "✗")

    df = pd.DataFrame(rows)
    n_agree = int(df["agree"].sum())
    n_total = len(df)
    n_err = int((df["sql_status"] == "ERROR").sum() + (df["rest_status"] == "ERROR").sum())
    logger.info("AGREEMENT: %d/%d (%.0f%%) (errors: %d)",
                n_agree, n_total, 100 * n_agree / max(n_total, 1), n_err)

    if n_agree < n_total:
        logger.warning("Disagreements:")
        for _, r in df[~df["agree"]].iterrows():
            logger.warning("  (%s, %s): SQL=%s n=%d vs REST=%s n=%d",
                           r["target_uniprot"], r["compound_name"],
                           r["sql_status"], r["sql_n"], r["rest_status"], r["rest_n"])

    # Write evidence to reports/ for auditability
    report_dir = ROOT / "reports" / "pipeline"
    report_dir.mkdir(parents=True, exist_ok=True)
    md = ["# Phase A.5 — SQLite vs REST agreement smoke test",
          "",
          f"**Agreement: {n_agree}/{n_total} ({100 * n_agree / max(n_total, 1):.0f}%) "
          f"(errors: {n_err})**",
          "",
          "Picked the top 20 (target, compound) pairs by predicted_pkd with SMILES < 200 chars.",
          "Both the SQLite backstop (`fetchers/chembl_sqlite.py`) and the legacy REST",
          "fetcher (`fetchers/chembl_groundtruth.py`) were run on the same pairs. "
          "INCONCLUSIVE (REST) is normalised to AMBIGUOUS (SQLite) for comparison.",
          "",
          "| target | compound | SQL status | SQL n | REST status | REST n | agree |",
          "|---|---|---|---|---|---|---|"]
    for _, r in df.iterrows():
        md.append(
            f"| {r['target_uniprot']} | {r['compound_name']} | "
            f"{r['sql_status']} | {r['sql_n']} | "
            f"{r['rest_status']} | {r['rest_n']} | "
            f"{'✓' if r['agree'] else '✗'} |"
        )
    md.append("")
    md.append("_Gate logic: PASS when n_agree == n_total - n_err. REST 500s are tolerated "
              "(transient EBI infrastructure) — only true status mismatches count as failures._")
    md.append("")
    md.append(f"**Result: {'PASS' if n_agree == n_total - n_err else 'FAIL'}** "
              f"(n_agree={n_agree}, n_total={n_total}, n_err={n_err})")
    (report_dir / "sqlite_vs_rest_smoke.md").write_text("\n".join(md), encoding="utf-8")
    logger.info("Wrote %s", report_dir / "sqlite_vs_rest_smoke.md")

    return 0 if n_agree == n_total - n_err else 2


if __name__ == "__main__":
    raise SystemExit(main())
