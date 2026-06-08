"""PERSEUS SENSITIVITY against a verified positive-persistence ledger.

Until now PERSEUS could only report SPECIFICITY (0 over-claims / 0 false positives) because the
positive class - compounds with durable post-cessation cognitive/neuroplastic change - was
empty (label budget ~381). This script consumes `data/raw/persistence_positive_ledger.csv`
(verified, cited positives assembled by the deep-research workflow; psychedelics /
dissociatives / psychoplastogens / neurotrophics) and measures, for the first time, what
fraction PERSEUS does NOT dismiss as null/abstain - i.e. recall.

This is expected to expose a GAP: PERSEUS's persistence head is built around senolytic-ablative
+ curated-mechanism substrate, so it will likely MISS psychoplastogen-mediated persistence
(5-HT2A / TrkB-cascade). Quantifying that miss is the deliverable - it scopes exactly where the
L4 TrkB/5-HT2A head and the psychoplastogen-target module must go.

CPU. Resolves missing SMILES via PubChem. Writes reports/pipeline/perseus_sensitivity_v1.md.
"""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from mammal_repurposing.engine.perseus import PerseusEngine, score_frame
from mammal_repurposing.validation.persistence_eval import sensitivity

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("perseus_sensitivity")

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
LEDGERS = [RAW / "clinical_outcomes_ledger.csv", RAW / "clinical_outcomes_ledger_EXTENSION.csv",
           RAW / "clinical_outcomes_ledger_CTGOV.csv", RAW / "clinical_outcomes_ledger_RESEARCH.csv"]
SMILES = RAW / "ledger_compound_smiles.csv"
POS = RAW / "persistence_positive_ledger.csv"
REPORT = ROOT / "reports" / "pipeline" / "perseus_sensitivity_v1.md"


# Curated PubChem query names: the salvaged ledger compound strings carry parentheticals /
# codes ("Psilocybin (psilocin)", "DOI (2,5-dimethoxy-4-iodoamphetamine)") that do not resolve
# by raw name. Map to a clean chemical name (ordered: more specific keys first).
_QUERY_TABLE = [
    ("psilocyb", "psilocybin"), ("5-meo-dmt", "5-methoxy-N,N-dimethyltryptamine"),
    ("n,n-dmt", "N,N-dimethyltryptamine"), ("lsd", "lysergic acid diethylamide"),
    ("lysergic", "lysergic acid diethylamide"), ("doi", "2,5-dimethoxy-4-iodoamphetamine"),
    ("mescaline", "mescaline"), ("hydroxynorketamine", "(2R,6R)-hydroxynorketamine"),
    ("ketamine", "ketamine"), ("ibogaine", "ibogaine"), ("scopolamine", "scopolamine"),
    ("nitrous", "nitrous oxide"), ("zuranolone", "zuranolone"), ("nsi-189", "NSI-189"),
    ("aaz-a-154", "AAZ-A-154"), ("zalsupindole", "AAZ-A-154"), ("gm-2505", "GM-2505"),
    ("bretisilocin", "GM-2505"), ("tabernanthalog", "tabernanthalog"),
    ("dmt", "N,N-dimethyltryptamine"),
]


def _query_name(compound: str) -> str:
    c = compound.lower()
    for key, q in _QUERY_TABLE:
        if key in c:
            return q
    return compound.split("(")[0].split("/")[0].strip()


def _resolve_smiles(df: pd.DataFrame) -> pd.DataFrame:
    if "smiles" in df and df["smiles"].notna().all():
        return df
    from mammal_repurposing.fetchers.pubchem import fetch_many_smiles
    df = df.copy()
    if "smiles" not in df:
        df["smiles"] = None
    pairs = [(c, _query_name(c)) for c in df.loc[df["smiles"].isna(), "compound"].unique()]
    L.info("Resolving %d positive-ledger SMILES via PubChem...", len(pairs))
    qhits = dict(zip([q for _, q in pairs],
                     fetch_many_smiles([(q, []) for _, q in pairs])))
    smap = {c: (qhits.get(q) or {}).get("smiles") for c, q in pairs}
    df["smiles"] = df.apply(
        lambda r: r["smiles"] if pd.notna(r.get("smiles")) else smap.get(r["compound"]), axis=1)
    # persist resolved SMILES back into the ledger (durable + auditable)
    if smap and POS.exists():
        full = pd.read_csv(POS)
        if "smiles" not in full:
            full["smiles"] = None
        full["smiles"] = full.apply(
            lambda r: r["smiles"] if pd.notna(r.get("smiles")) else smap.get(r["compound"]), axis=1)
        full.to_csv(POS, index=False)
        L.info("Persisted SMILES into %s", POS)
    return df


def main() -> int:
    if not POS.exists():
        L.warning("Positive ledger %s not present yet - run the deep-research workflow first.", POS)
        return 0
    pos = pd.read_csv(POS)
    # keep verified, structure-scoreable small molecules
    keep = pos[(pos.get("verification_status", "verified") == "verified")
               & (pos.get("is_small_molecule", True).astype(str).str.lower().isin(["true", "1", "yes"]))].copy()
    keep = _resolve_smiles(keep)
    n_no_smiles = int(keep["smiles"].isna().sum())
    keep = keep[keep["smiles"].notna()]

    eng = PerseusEngine(LEDGERS, SMILES, RAW / "persistence_axis_classes.csv",
                        RAW / "persistence_axis_overrides.csv")
    scored = score_frame(eng, keep.rename(columns={"compound": "query_id"}), dedup_salts=False)
    meta_cols = [c for c in ("compound", "domain", "drug_class") if c in keep.columns]
    scored = scored.merge(keep[meta_cols], on="compound", how="left")
    if "drug_class" not in scored:
        scored["drug_class"] = "?"

    records = [{"compound": r["compound"], "persistence_verdict": r["persistence_verdict"],
                "domain": r.get("domain", "?")} for _, r in scored.iterrows()]
    s = sensitivity(records)

    Ls = ["# PERSEUS sensitivity - against the verified positive-persistence ledger", "",
          "First recall measurement for PERSEUS, enabled by a non-empty verified-positive "
          "ledger (`data/raw/persistence_positive_ledger.csv`; psychedelics / dissociatives / "
          "psychoplastogens / neurotrophics with cited durable post-cessation effects). A "
          "positive is FLAGGED if PERSEUS asserts any durability (not null/abstain/excluded). "
          "Reproduced by `scripts/107_perseus_sensitivity.py`.", "",
          f"## Recall: **{s['n_flagged']} / {s['n']}** = "
          f"{s['sensitivity']:.2f} " + ("(n/a)" if s["n"] == 0 else "") +
          f"  (+{n_no_smiles} verified positives had no resolvable SMILES)", "",
          "Per domain:", ""]
    for d, v in sorted(s["by_domain"].items()):
        Ls.append(f"- {d}: {v['flagged']}/{v['n']} flagged")
    Ls += ["", "| compound | class | domain | PERSEUS verdict | flagged? |", "|---|---|---|---|---|"]
    for _, r in scored.sort_values("persistence_verdict").iterrows():
        from mammal_repurposing.validation.persistence_eval import VERDICT_DURABILITY
        fl = VERDICT_DURABILITY.get(r["persistence_verdict"], 0) >= 1
        Ls.append(f"| {r['compound']} | {r.get('drug_class', '?')} | {r.get('domain', '?')} | "
                  f"{r['persistence_verdict']} | {'yes' if fl else 'no'} |")
    Ls += ["", "## Reading", "",
           "A LOW recall here is the expected, honest finding: PERSEUS's persistence head is "
           "built around senolytic-ablative + curated-mechanism substrate and the size-matched "
           "DTI channel (BCL2 confirmed, NTRK2 plasticity), none of which captures "
           "psychoplastogen-mediated (5-HT2A / TrkB-cascade) persistence. The missed compounds "
           "scope exactly where the L4 TrkB/5-HT2A access head and a psychoplastogen-target "
           "module must extend the engine. This is the recall baseline those additions must beat.", ""]
    REPORT.write_text("\n".join(Ls), encoding="utf-8")
    L.info("Wrote %s", REPORT)
    L.info("SENSITIVITY: %d/%d = %.2f flagged; by domain %s",
           s["n_flagged"], s["n"], s["sensitivity"] if s["n"] else float("nan"), s["by_domain"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
