"""Does the DTI head rank anything, and does it rank anything NEW?

Every existing measurement of the MAMMAL DTI head in this repository is taken on three to five
anchor compounds. scripts/114 reported AUROC 0.09 for AMPA-orthosteric on THREE actives, and a
three-active AUROC has a 95% null interval that spans nearly the whole unit line. Those numbers
have been load-bearing for the allosteric-blindness finding (G2) and were about to become
load-bearing for closing de novo design. They deserve a real n.

This runs two stages against ChEMBL 36, pre-registered before the model was loaded.

STAGE 1, CAN IT RANK AT ALL. Twelve cognition-relevant human targets, up to N_PER actives each
(pChEMBL >= 7, binding assays only), against a size-matched hard negative pool drawn from actives
at OTHER panel targets with no recorded activity at this one. Hard negatives, not random chemistry:
every negative is a real drug-like bioactive. Per-target AUROC plus a permutation test.

STAGE 2, THE Lo-Hi SPLIT. For every target that clears stage 1, split its actives by maximum ECFP4
Tanimoto to that target's ten REFERENCE ligands, defined as the ten earliest-published actives,
which is the subset most likely to sit in any pretraining corpus and the subset a chemist would
name unprompted. Hi band is >= 0.4, Lo band is < 0.4, the boundary the Lo-Hi benchmark uses for the
hard regime. Both bands are scored against the identical negative pool, so the only thing that
changes between them is distance from the familiar.

PRE-REGISTERED PREDICTIONS, written before running:
  P1. Stage 1 passes on at least some targets. A head that scored at chance everywhere would not
      have been shipped, and the existing 0.09/0.26 are small-sample noise rather than measurement.
  P2. The Lo band is materially worse than the Hi band. If it is, the head is a near-neighbour
      lookup and target-directed generation is refuted for the same reason structure-based
      generation already was.
  P3. GRIA1 specifically will rank better than the 0.09/0.26 that scripts/114 reported, because
      n = 3 and n = 4 cannot resolve anything.

Whatever comes back, it is a measurement this project has never made. The interesting outcomes are
symmetric: if P2 holds, the de novo door closes on evidence rather than on n = 3. If P2 fails and
the head generalises, then G2 is overstated and the allosteric finding needs re-examination at
honest n, which is a results-changing correction that must be declared.

Writes reports/pipeline/dti_scale_lohi_v1.md.
"""
from __future__ import annotations

import json
import logging
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd

from mammal_repurposing.engine.persistence_dti import calibrate_target
from mammal_repurposing.provenance.trailer import stamp

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("dti_lohi")

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "reports" / "pipeline" / "dti_scale_lohi_v1.md"
RAWOUT = ROOT / "data" / "interim" / "dti_scale_lohi_scores.csv"
CHEMBL = Path.home() / ".data" / "chembl" / "36" / "chembl_36.db"

# Pre-registered panel. Chosen for relevance to this project's mechanism classes and for having
# enough ChEMBL data to resolve an AUROC, NOT for expected performance.
PANEL = [
    "SLC6A3",   # DAT, methylphenidate
    "SLC6A2",   # NET, atomoxetine
    "ADORA2A",  # caffeine
    "ACHE",     # donepezil, the assay-dependence compound
    "CHRNA7",   # nicotinic, the v1 allosteric-blindness target
    "GRIA1",    # AMPA, the G2 target
    "PDE4D",    # zatolmilast / rolipram class
    "MMP9",     # PNN/ECM
    "DRD2", "HTR2A", "HRH3", "SLC6A4",
]
N_PER = 120           # actives per target, and an equal number of hard negatives
MIN_ACTIVES = 40      # below this a target is reported but not interpreted
N_REFERENCE = 10      # earliest-published actives that define the familiar region
LOHI_BOUNDARY = 0.40  # the Lo-Hi hard-regime boundary
MIN_BAND = 15         # a band below this is reported but not interpreted
PASS_AUROC = 0.70     # the project's existing channel gate
SEED = 0


def fetch(con: sqlite3.Connection) -> tuple[pd.DataFrame, dict[str, str]]:
    """Pull actives + sequences for the panel. One row per (gene, molecule)."""
    marks = ",".join("?" * len(PANEL))
    q = f"""
    SELECT cs.component_synonym AS gene,
           act.molregno           AS molregno,
           MAX(act.pchembl_value) AS pchembl,
           MIN(d.year)            AS first_year,
           cst.canonical_smiles   AS smiles
    FROM component_synonyms cs
    JOIN target_components tc ON tc.component_id = cs.component_id
    JOIN target_dictionary td ON td.tid = tc.tid
    JOIN assays a              ON a.tid = td.tid
    JOIN activities act        ON act.assay_id = a.assay_id
    JOIN compound_structures cst ON cst.molregno = act.molregno
    LEFT JOIN docs d           ON d.doc_id = act.doc_id
    WHERE cs.syn_type = 'GENE_SYMBOL' AND cs.component_synonym IN ({marks})
      AND td.organism = 'Homo sapiens' AND td.target_type = 'SINGLE PROTEIN'
      AND a.assay_type = 'B' AND act.pchembl_value >= 7
      AND cst.canonical_smiles IS NOT NULL
    GROUP BY cs.component_synonym, act.molregno
    """
    acts = pd.read_sql_query(q, con, params=PANEL)

    qs = f"""
    SELECT cs.component_synonym AS gene, MAX(seq.sequence) AS sequence
    FROM component_synonyms cs
    JOIN component_sequences seq ON seq.component_id = cs.component_id
    JOIN target_components tc ON tc.component_id = cs.component_id
    JOIN target_dictionary td ON td.tid = tc.tid
    WHERE cs.syn_type = 'GENE_SYMBOL' AND cs.component_synonym IN ({marks})
      AND td.organism = 'Homo sapiens' AND td.target_type = 'SINGLE PROTEIN'
    GROUP BY cs.component_synonym
    """
    seqs = dict(con.execute(qs, PANEL).fetchall())
    return acts, seqs


def any_activity(con: sqlite3.Connection, gene: str) -> set[int]:
    """Every molregno with ANY recorded activity at this gene, at any potency.

    Used to keep the negative pool clean: a compound that is merely weak at the target is not a
    negative, it is a weak active, and scoring it as a negative would understate the head.
    """
    q = """
    SELECT DISTINCT act.molregno
    FROM component_synonyms cs
    JOIN target_components tc ON tc.component_id = cs.component_id
    JOIN target_dictionary td ON td.tid = tc.tid
    JOIN assays a ON a.tid = td.tid
    JOIN activities act ON act.assay_id = a.assay_id
    WHERE cs.syn_type='GENE_SYMBOL' AND cs.component_synonym = ?
      AND td.organism='Homo sapiens' AND td.target_type='SINGLE PROTEIN'
    """
    return {r[0] for r in con.execute(q, (gene,))}


def fingerprints(smiles: list[str]):
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdFingerprintGenerator
    RDLogger.DisableLog("rdApp.*")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    out = []
    for s in smiles:
        m = Chem.MolFromSmiles(s)
        out.append(None if m is None else gen.GetFingerprint(m))
    return out


def max_sim_to(fps, ref_fps) -> list[float]:
    from rdkit import DataStructs
    ref = [f for f in ref_fps if f is not None]
    out = []
    for f in fps:
        if f is None or not ref:
            out.append(float("nan"))
        else:
            out.append(max(DataStructs.BulkTanimotoSimilarity(f, ref)))
    return out


def build_rows(con: sqlite3.Connection) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    acts, seqs = fetch(con)
    L.info("ChEMBL returned %d (gene, molecule) actives across %d genes",
           len(acts), acts["gene"].nunique())

    known = {g: any_activity(con, g) for g in PANEL}
    rows = []
    for gene in PANEL:
        if gene not in seqs:
            L.warning("%s: no sequence, skipped", gene)
            continue
        sub = acts[acts["gene"] == gene]
        if len(sub) < MIN_ACTIVES:
            L.warning("%s: only %d actives", gene, len(sub))
        # Deterministic sample, but keep the earliest-published ones so references are available.
        sub = sub.sort_values(["first_year", "molregno"], na_position="last")
        refs = sub.head(N_REFERENCE)
        rest = sub.iloc[N_REFERENCE:]
        take = min(N_PER, len(rest))
        pick = rest.iloc[rng.permutation(len(rest))[:take]] if take else rest
        for _, r in pick.iterrows():
            rows.append(dict(gene=gene, role="active", molregno=int(r.molregno),
                             smiles=r.smiles, pchembl=float(r.pchembl), seq=seqs[gene]))
        for _, r in refs.iterrows():
            rows.append(dict(gene=gene, role="reference", molregno=int(r.molregno),
                             smiles=r.smiles, pchembl=float(r.pchembl), seq=seqs[gene]))

        # Hard negatives: actives at OTHER panel targets, with no activity at all here.
        pool = acts[(acts["gene"] != gene) & (~acts["molregno"].isin(known[gene]))]
        pool = pool.drop_duplicates("molregno")
        n_neg = min(take, len(pool))
        negp = pool.iloc[rng.permutation(len(pool))[:n_neg]]
        for _, r in negp.iterrows():
            rows.append(dict(gene=gene, role="negative", molregno=int(r.molregno),
                             smiles=r.smiles, pchembl=float("nan"), seq=seqs[gene]))
        L.info("%s: %d actives + %d references + %d hard negatives",
               gene, take, len(refs), n_neg)
    return pd.DataFrame(rows)


def score(df: pd.DataFrame) -> pd.DataFrame:
    from mammal_repurposing.scoring.dti import score_batch_safe
    from mammal_repurposing.scoring.model_loader import load_dti_model
    model, tok = load_dti_model()
    pkds: list[float] = []
    pairs = list(zip(df["seq"], df["smiles"], strict=True))
    ids = [f"{g}|{m}" for g, m in zip(df["gene"], df["molregno"], strict=True)]
    step = 8
    for i in range(0, len(pairs), step):
        pkds.extend(score_batch_safe(model, tok, pairs[i:i + step], sample_ids=ids[i:i + step]))
        if (i // step) % 25 == 0:
            L.info("  scored %d / %d", min(i + step, len(pairs)), len(pairs))
    df = df.copy()
    df["predicted_pkd"] = pkds
    return df


def analyse(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    stage1, stage2 = [], []
    for gene, g in df.groupby("gene"):
        neg = g[g.role == "negative"]["predicted_pkd"].dropna().tolist()
        pos = g[g.role == "active"]["predicted_pkd"].dropna().tolist()
        if len(pos) < 3 or len(neg) < 3:
            continue
        c = calibrate_target(pos, neg, min_pos=3, min_auroc=PASS_AUROC)
        stage1.append(dict(gene=gene, n_pos=len(pos), n_neg=len(neg),
                           auroc=c["auroc"], perm_p=c["perm_p"],
                           ranks=bool(c["auroc"] >= PASS_AUROC and c["perm_p"] < 0.05)))

        act = g[g.role == "active"].dropna(subset=["predicted_pkd"])
        ref = g[g.role == "reference"]
        fps_a = fingerprints(act["smiles"].tolist())
        fps_r = fingerprints(ref["smiles"].tolist())
        sims = max_sim_to(fps_a, fps_r)
        act = act.assign(max_sim_ref=sims)
        for band, mask in (("Hi", act.max_sim_ref >= LOHI_BOUNDARY),
                           ("Lo", act.max_sim_ref < LOHI_BOUNDARY)):
            p = act[mask]["predicted_pkd"].tolist()
            if len(p) < 3:
                stage2.append(dict(gene=gene, band=band, n=len(p), auroc=float("nan"),
                                   perm_p=float("nan"), interpretable=False))
                continue
            cb = calibrate_target(p, neg, min_pos=3, min_auroc=PASS_AUROC)
            stage2.append(dict(gene=gene, band=band, n=len(p), auroc=cb["auroc"],
                               perm_p=cb["perm_p"], interpretable=len(p) >= MIN_BAND))
    return pd.DataFrame(stage1), pd.DataFrame(stage2)


def write_report(s1: pd.DataFrame, s2: pd.DataFrame, df: pd.DataFrame) -> None:
    lines = [
        "# Does the DTI head rank anything, and does it rank anything new?",
        "",
        "Pre-registered two-stage test of the MAMMAL DTI head at ChEMBL scale. Every previous",
        "measurement in this repository used three to five anchor compounds; this uses up to",
        f"{N_PER} actives per target against an equal number of hard negatives, where a hard",
        "negative is a compound that is a confirmed active at a different panel target and has no",
        "recorded activity at this one. Predictions P1 to P3 were written before the model loaded;",
        "see the docstring of `scripts/133_dti_scale_lohi.py`.",
        "",
        "## Stage 1: can it rank at all?",
        "",
        "| gene | n actives | n negatives | AUROC | perm-p | ranks |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for _, r in s1.sort_values("auroc", ascending=False).iterrows():
        lines.append(f"| {r.gene} | {r.n_pos} | {r.n_neg} | {r.auroc:.2f} | "
                     f"{r.perm_p:.3f} | {'**yes**' if r.ranks else 'no'} |")
    n_rank = int(s1["ranks"].sum())
    lines += [
        "",
        f"{n_rank} of {len(s1)} targets clear the project's existing channel gate "
        f"(AUROC >= {PASS_AUROC}, permutation p < 0.05).",
        "",
        "## Stage 2: the Lo-Hi split",
        "",
        "Actives split by maximum ECFP4 Tanimoto to the ten earliest-published actives for the same",
        f"target. Hi is >= {LOHI_BOUNDARY}, Lo is below it. Both bands are scored against the",
        "identical negative pool, so the only variable is distance from the familiar region. A band",
        f"under n = {MIN_BAND} is shown but not interpreted.",
        "",
        "| gene | band | n | AUROC | perm-p |",
        "|---|---|---:|---:|---:|",
    ]
    for _, r in s2.sort_values(["gene", "band"]).iterrows():
        au = "n/a" if pd.isna(r.auroc) else f"{r.auroc:.2f}"
        pp = "n/a" if pd.isna(r.perm_p) else f"{r.perm_p:.3f}"
        flag = "" if r.interpretable else " *"
        lines.append(f"| {r.gene} | {r.band}{flag} | {r.n} | {au} | {pp} |")
    lines += ["", "`*` band too small to interpret.", ""]

    piv = s2[s2.interpretable].pivot_table(index="gene", columns="band", values="auroc")
    if {"Hi", "Lo"}.issubset(piv.columns):
        paired = piv.dropna()
        if len(paired):
            delta = (paired["Hi"] - paired["Lo"]).mean()
            lines += [
                f"Across the {len(paired)} targets where both bands are interpretable, mean AUROC "
                f"in the Hi band exceeds the Lo band by {delta:+.3f}.",
                "",
            ]
    lines += [
        "## Reading",
        "",
        "Fill this in by hand after reading the tables. Do not let the script assert a conclusion.",
        "",
        f"Scored {len(df)} (target, compound) pairs. Raw scores in "
        "`data/interim/dti_scale_lohi_scores.csv`.",
        "",
    ]
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(stamp("\n".join(lines), "scripts/133_dti_scale_lohi.py"), encoding="utf-8")


def main() -> int:
    if not CHEMBL.exists():
        L.error("ChEMBL not found at %s", CHEMBL)
        return 2
    con = sqlite3.connect(f"file:{CHEMBL}?mode=ro", uri=True)
    df = build_rows(con)
    con.close()
    L.info("Scoring %d (target, compound) pairs", len(df))
    df = score(df)
    RAWOUT.parent.mkdir(parents=True, exist_ok=True)
    df.drop(columns=["seq"]).to_csv(RAWOUT, index=False)
    s1, s2 = analyse(df)
    write_report(s1, s2, df)
    print(json.dumps({"stage1": s1.to_dict("records"), "stage2": s2.to_dict("records")},
                     indent=1, default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
