"""cis-Mendelian-randomisation screen: druggable genes against human cognitive performance.

THE PIVOT. Three boundaries close the previous approach, and all three were measured here rather
than assumed: the DTI head ranks at chance (pooled AUROC 0.468 over 12 targets at n = 120 a side),
the healthy-adult literature yields four modest ACUTE enhancers and zero durable ones, and the
prospective registry cannot be resolved in useful time (25 of 26 due rows unpublished, and the one
that resolved was retracted on attribution).

So the source of ground truth changes. Human genetics supplies causal evidence that does not depend
on anybody running a trial, publishing it, or reporting it honestly. A cis-eQTL is a lifetime,
randomised perturbation of one gene's expression; regressing cognitive performance on it asks
whether moving that protein moves cognition, in humans, at zero wet-lab cost.

WHAT THIS IS AND IS NOT. It is a SCREEN. A Wald-ratio cis-MR on a single fine-mapped instrument is
a hypothesis generator, not a causal claim. Two things must follow before any target is believed:
colocalisation (to rule out that the eQTL and the GWAS signal are distinct variants in linkage
disequilibrium, which is the dominant false-positive mode), and a pleiotropy screen. Neither is in
this file. Its output is a ranked list to take into those steps, and it says so in its own report.

DATA, all openly downloadable, nothing behind an application:
  exposure  eQTL Catalogue SuSiE credible sets, human cortex:
            ROSMAP DLPFC n=560, CommonMind DLPFC n=586, BrainSeq DLPFC n=479, GTEx cortex n=205
  outcome   Lee et al. 2018 Cognitive Performance GWAS, GCST006572, n=257,841
  targets   genes carrying a ChEMBL drug-mechanism annotation (891 genes, 447 with an approved drug)

TWO TRAPS THIS HANDLES EXPLICITLY, because both silently invert results rather than failing loudly:

  ALLELE HARMONISATION. The eQTL variant id encodes chr_pos_REF_ALT and its beta is for ALT. The
  GWAS reports beta for its own effect allele A1. If A1 is the eQTL REF the GWAS beta must be
  sign-flipped before the ratio is taken. Getting this wrong does not produce noise, it produces a
  confident result with the direction reversed, which is worse than no result at all.

  PALINDROMIC VARIANTS. A/T and C/G SNPs cannot be strand-resolved from alleles alone. Where the
  allele frequency is near one half, even frequency cannot break the tie. Those are DROPPED rather
  than guessed, because a coin-flip on strand is a coin-flip on direction.

A CONFOUND THAT IS NOT HANDLED HERE AND MUST BE, stated so it is not discovered later. Cognitive
Performance is genetically correlated with educational attainment at roughly 0.9, and EA carries
social and demographic pathways that have nothing to do with neurobiology. A hit in this screen may
act on schooling rather than on cognitive machinery. The mitigation is multivariable MR conditioning
on EA (GCST006442, n=1,131,881), and it is the next thing to build, not an optional refinement.

Writes data/interim/cis_mr_instruments.csv, data/raw/cis_mr_screen.csv and
reports/pipeline/cis_mr_screen_v1.md.
"""
from __future__ import annotations

import csv
import gzip
import logging
import math
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.provenance.trailer import stamp  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("cis_mr")

HOME = Path.home()
GWAS = HOME / ".data" / "gwas" / "GWAS_CP_all.txt"
EQTL_DIR = HOME / ".data" / "eqtl"
HGNC = HOME / ".data" / "ref" / "hgnc.txt"
DRUGGABLE = ROOT / "data" / "interim" / "druggable_genes_chembl.csv"
INSTR = ROOT / "data" / "interim" / "cis_mr_instruments.csv"
OUT = ROOT / "data" / "raw" / "cis_mr_screen.csv"
REPORT = ROOT / "reports" / "pipeline" / "cis_mr_screen_v1.md"

#: dataset id -> (study, tissue, n). Cortical gene-expression QTLs only.
DATASETS = {
    "QTD000434": ("ROSMAP", "brain (DLPFC)", 560),
    "QTD000075": ("CommonMind", "brain (DLPFC)", 586),
    "QTD000051": ("BrainSeq", "brain (DLPFC)", 479),
    "QTD000171": ("GTEx", "brain (cortex)", 205),
}

MIN_PIP = 0.20          # fine-mapping confidence for the chosen instrument
MAX_EQTL_P = 5e-8       # instrument must be a genome-wide-significant cis-eQTL
PALINDROMIC_EAF = 0.42  # below this distance from 0.5 a palindrome cannot be strand-resolved
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def ensembl_to_symbol() -> dict[str, str]:
    with HGNC.open(encoding="utf-8") as f:
        return {r["ensembl_gene_id"]: r["symbol"]
                for r in csv.DictReader(f, delimiter="\t") if r.get("ensembl_gene_id")}


def candidate_instruments(druggable: set[str], sym: dict[str, str]) -> pd.DataFrame:
    """ALL usable credible-set variants per (gene, dataset), PIP-ranked.

    The first version kept only the single highest-PIP variant per gene, and 74 of 248 instruments
    were then lost because that one variant is absent from the GWAS. Keeping the whole credible set
    lets the chooser fall back to the next-best variant that the outcome study actually measured,
    which costs nothing: every variant in a credible set is instrumenting the same signal, so the
    fallback is a different tag for the same thing rather than a different hypothesis.
    """
    rows: list[dict] = []
    for ds, (study, tissue, n) in DATASETS.items():
        path = EQTL_DIR / f"{ds}.credible_sets.tsv.gz"
        if not path.exists():
            L.warning("%s missing, skipped", path.name)
            continue
        with gzip.open(path, "rt", encoding="utf-8") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                gene = sym.get(r["gene_id"])
                if gene is None or gene not in druggable:
                    continue
                try:
                    pip, p = float(r["pip"]), float(r["pvalue"])
                    beta, se = float(r["beta"]), float(r["se"])
                except (TypeError, ValueError):
                    continue
                if pip < MIN_PIP or p > MAX_EQTL_P or se <= 0:
                    continue
                parts = r["variant"].split("_")
                if len(parts) < 4:
                    continue
                rows.append({"gene": gene, "dataset": ds, "study": study, "tissue": tissue,
                             "eqtl_n": n, "rsid": r["rsid"], "variant": r["variant"],
                             "cs_id": r.get("cs_id", ""),
                             "ref": parts[2].upper(), "alt": parts[3].upper(),
                             "pip": pip, "eqtl_p": p, "eqtl_beta": beta, "eqtl_se": se})
        L.info("%s (%s, %s): %d candidate variants over %d druggable genes",
               ds, study, tissue,
               sum(1 for x in rows if x["dataset"] == ds),
               len({x["gene"] for x in rows if x["dataset"] == ds}))
    return pd.DataFrame(rows)


def choose_instruments(cand: pd.DataFrame, in_gwas: set[str]) -> pd.DataFrame:
    """Per (gene, dataset), the highest-PIP candidate that the GWAS actually measured."""
    cand = cand[cand["rsid"].isin(in_gwas)]
    if cand.empty:
        return cand
    idx = cand.groupby(["gene", "dataset"])["pip"].idxmax()
    return cand.loc[idx].reset_index(drop=True)


def harmonise(gwas_a1: str, gwas_a2: str, ref: str, alt: str, eaf: float) -> tuple[int, str]:
    """Sign to apply to the GWAS beta so it refers to the eQTL ALT allele, or a rejection.

    Returns (+1 | -1, "") on success, or (0, reason). Rejecting is always preferable to guessing:
    a wrong sign here yields a confident result pointing the wrong way.
    """
    a1, a2 = gwas_a1.upper(), gwas_a2.upper()
    if len(ref) != 1 or len(alt) != 1 or len(a1) != 1 or len(a2) != 1:
        return 0, "indel or multi-base allele; not harmonised"
    palindromic = COMPLEMENT.get(ref) == alt
    if palindromic:
        if eaf is None or math.isnan(eaf) or abs(eaf - 0.5) < (0.5 - PALINDROMIC_EAF):
            return 0, "palindromic variant with ambiguous frequency; strand unresolvable"
    if {a1, a2} == {ref, alt}:
        return (1, "") if a1 == alt else (-1, "")
    flipped = {COMPLEMENT.get(a1, "?"), COMPLEMENT.get(a2, "?")}
    if flipped == {ref, alt}:
        return (1, "") if COMPLEMENT.get(a1) == alt else (-1, "")
    return 0, f"alleles do not match ({a1}/{a2} vs {ref}/{alt})"


def load_gwas_for(rsids: set[str]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    with GWAS.open(encoding="utf-8") as f:
        rdr = csv.DictReader(f, delimiter="\t")
        for i, r in enumerate(rdr):
            if r["MarkerName"] in rsids:
                try:
                    out[r["MarkerName"]] = {"a1": r["A1"], "a2": r["A2"],
                                            "eaf": float(r["EAF"]), "beta": float(r["Beta"]),
                                            "se": float(r["SE"]), "p": float(r["Pval"])}
                except (TypeError, ValueError):
                    pass
            if i and i % 2_000_000 == 0:
                L.info("  scanned %dM GWAS rows, matched %d", i // 1_000_000, len(out))
    return out


def bh_fdr(pvals: list[float]) -> list[float]:
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [1.0] * n
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        idx = n - rank + 1
        prev = min(prev, pvals[i] * n / idx)
        q[i] = prev
    return q


def main() -> int:
    for p in (GWAS, HGNC, DRUGGABLE):
        if not p.exists():
            L.error("missing input: %s", p)
            return 2
    sym = ensembl_to_symbol()
    druggable = set(pd.read_csv(DRUGGABLE)["gene"].astype(str))
    L.info("druggable genes: %d | ensembl->symbol: %d", len(druggable), len(sym))

    cand = candidate_instruments(druggable, sym)
    if cand.empty:
        L.error("no candidate instruments; are the credible-set files present in %s?", EQTL_DIR)
        return 2
    L.info("candidates: %d variants over %d unique genes",
           len(cand), cand["gene"].nunique())

    L.info("scanning the CP GWAS for %d candidate rsids ...", cand["rsid"].nunique())
    g = load_gwas_for(set(cand["rsid"]))
    L.info("GWAS measured %d of %d candidate rsids", len(g), cand["rsid"].nunique())

    ins = choose_instruments(cand, set(g))
    if ins.empty:
        L.error("no candidate instrument is present in the GWAS")
        return 2
    INSTR.parent.mkdir(parents=True, exist_ok=True)
    ins.to_csv(INSTR, index=False, encoding="utf-8")
    L.info("instruments chosen: %d over %d unique genes (fallback to next-best PIP where the "
           "top variant was unmeasured)", len(ins), ins["gene"].nunique())

    rows = []
    for _, r in ins.iterrows():
        hit = g.get(r["rsid"])
        if hit is None:
            rows.append({**r.to_dict(), "status": "rsid not in GWAS"})
            continue
        sign, why = harmonise(hit["a1"], hit["a2"], r["ref"], r["alt"], hit["eaf"])
        if sign == 0:
            rows.append({**r.to_dict(), "status": why})
            continue
        b_out, se_out = sign * hit["beta"], hit["se"]
        b_exp, se_exp = r["eqtl_beta"], r["eqtl_se"]
        wald = b_out / b_exp
        # delta-method SE for a ratio, first order (standard for a single-instrument Wald ratio)
        se_wald = abs(wald) * math.sqrt((se_out / b_out) ** 2 + (se_exp / b_exp) ** 2) \
            if b_out != 0 else float("nan")
        z = wald / se_wald if se_wald and not math.isnan(se_wald) else float("nan")
        rows.append({**r.to_dict(), "status": "ok", "gwas_a1": hit["a1"], "gwas_a2": hit["a2"],
                     "gwas_eaf": hit["eaf"], "harmonise_sign": sign,
                     "gwas_beta_aligned": b_out, "gwas_se": se_out, "gwas_p": hit["p"],
                     "wald_ratio": wald, "wald_se": se_wald, "wald_z": z})

    df = pd.DataFrame(rows)
    ok = df[df["status"] == "ok"].copy()
    # A gene appearing in several eQTL datasets LOOKS replicated and usually is not. All four
    # datasets are tested against the SAME outcome GWAS, so if they select the same variant, or
    # variants in perfect LD, the outcome beta is literally the same number and only the eQTL
    # denominator changes. Count distinct outcome effects, not distinct datasets.
    if len(ok):
        n_distinct = ok.groupby("gene")["gwas_beta_aligned"].transform(lambda s: s.round(6).nunique())
        ok["n_datasets"] = ok.groupby("gene")["dataset"].transform("nunique")
        ok["n_distinct_outcome_betas"] = n_distinct
        ok["independent_outcome_evidence"] = n_distinct > 1
    if len(ok):
        from scipy.stats import norm
        ok["mr_p"] = 2 * norm.sf(ok["wald_z"].abs())
        ok["mr_q"] = bh_fdr(ok["mr_p"].tolist())
        df = df.merge(ok[["gene", "dataset", "mr_p", "mr_q", "n_datasets",
                          "n_distinct_outcome_betas", "independent_outcome_evidence"]],
                      on=["gene", "dataset"], how="left")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False, encoding="utf-8")
    write_report(df, ok)
    L.info("usable tests: %d | FDR<0.05: %d", len(ok),
           int((ok["mr_q"] < 0.05).sum()) if len(ok) else 0)
    return 0


def write_report(df: pd.DataFrame, ok: pd.DataFrame) -> None:
    dropped = df[df["status"] != "ok"]["status"].value_counts()
    lines = [
        "# cis-MR screen: druggable genes against human cognitive performance",
        "",
        "**A SCREEN, not a result.** A Wald-ratio cis-MR on one fine-mapped instrument generates",
        "hypotheses. Two things must follow before any target here is believed: colocalisation, to",
        "rule out that the eQTL and the cognition signal are distinct variants in linkage",
        "disequilibrium, which is the dominant false-positive mode; and a pleiotropy screen.",
        "Neither is in this file.",
        "",
        "## Why genetics replaced the clinical literature",
        "",
        "Three measured boundaries, not a change of taste. The DTI head ranks at chance (pooled",
        "AUROC 0.468). The healthy-adult literature yields four modest ACUTE enhancers and zero",
        "durable ones. The prospective registry cannot be resolved in useful time: 25 of 26 due rows",
        "unpublished, and the single row that did resolve was retracted on attribution. A cis-eQTL",
        "is a lifetime randomised perturbation of one gene, and it needs nobody to run, publish or",
        "honestly report a trial.",
        "",
        "## Data",
        "",
        "| layer | source | n |",
        "|---|---|---|",
        "| outcome | Lee 2018 Cognitive Performance, GCST006572 | 257,841 |",
    ]
    for ds, (study, tissue, n) in DATASETS.items():
        lines.append(f"| exposure | {study} {tissue} (`{ds}`, eQTL Catalogue SuSiE) | {n} |")
    lines += [
        "| targets | ChEMBL drug-mechanism genes | 891 (447 with an approved drug) |",
        "",
        f"Instrument rule: highest-PIP variant in the gene's credible set, PIP >= {MIN_PIP}, "
        f"cis-eQTL p <= {MAX_EQTL_P:g}.",
        "",
        "## Harmonisation, which is where MR silently inverts",
        "",
        "The eQTL beta refers to the ALT allele in `chr_pos_REF_ALT`; the GWAS beta refers to its",
        "own A1. Where A1 is the eQTL REF the GWAS beta is sign-flipped before the ratio is taken.",
        "Getting this wrong yields a confident result with the direction reversed, which is worse",
        "than no result. Palindromic A/T and C/G variants whose frequency is near one half cannot",
        "be strand-resolved and are DROPPED rather than guessed.",
        "",
        "| outcome | tests |", "|---|---:|",
        f"| usable | {len(ok)} |",
    ]
    for reason, k in dropped.items():
        lines.append(f"| dropped: {reason} | {k} |")

    if len(ok):
        top = ok.sort_values("mr_p").head(30)
        lines += [
            "",
            f"**{int((ok['mr_q'] < 0.05).sum())} gene-tissue tests at FDR < 0.05.** Ranked by "
            "p-value; `wald_ratio` is the change in cognitive performance (SD units) per unit "
            "increase in normalised brain expression, so its SIGN is the directional requirement: "
            "positive means a drug should RAISE the target's activity.",
            "",
            "| gene | tissue | rsid | PIP | eQTL beta | CP beta (aligned) | Wald | p | q |",
            "|---|---|---|---:|---:|---:|---:|---:|---:|",
        ]
        for _, r in top.iterrows():
            lines.append(
                f"| **{r['gene']}** | {r['study']} | {r['rsid']} | {r['pip']:.2f} | "
                f"{r['eqtl_beta']:+.3f} | {r['gwas_beta_aligned']:+.4f} | "
                f"{r['wald_ratio']:+.4f} | {r['mr_p']:.2e} | {r['mr_q']:.3f} |")
    if len(ok):
        multi = ok[ok["n_datasets"] > 1]
        fake = multi[~multi["independent_outcome_evidence"]]["gene"].nunique()
        real = multi[multi["independent_outcome_evidence"]]["gene"].nunique()
        lines += [
            "",
            "## Why \"replicated across cohorts\" is mostly an illusion here",
            "",
            "A gene appearing in several eQTL datasets reads as replication and usually is not. All",
            "four exposure datasets are tested against the SAME outcome GWAS. Where they select the",
            "same variant, or variants in perfect LD, the outcome beta is literally the same number",
            "and only the eQTL denominator changes, so the second test carries no new information",
            "about cognition.",
            "",
            f"Of the genes appearing in more than one dataset, **{fake}** repeat an identical",
            f"outcome effect and **{real}** carry more than one distinct outcome estimate. Even the",
            "latter are cis-variants at one locus and therefore LD-correlated, so they are not",
            "independent either. Nothing in this screen is replicated in the sense that word",
            "normally carries; genuine replication requires a second, non-overlapping cognition",
            "GWAS.",
            "",
        ]
    lines += [
        "",
        "## The confound this does not handle",
        "",
        "Cognitive Performance is genetically correlated with educational attainment at roughly",
        "0.9, and EA carries social and demographic pathways with no neurobiological content. A hit",
        "above may act on schooling rather than on cognitive machinery. The mitigation is",
        "multivariable MR conditioning on EA (GCST006442, n = 1,131,881). That is the next thing to",
        "build, not an optional refinement, and until it exists every row here is provisional.",
        "",
        "## One tension to decide deliberately",
        "",
        "Screening candidates against adverse psychiatric phenotypes is sensible, but cognition and",
        "schizophrenia share substantial genetic architecture, so a blanket rejection of",
        "schizophrenia-associated genes would discard some of the strongest cognition signals. That",
        "should be a recorded per-target safety judgement, not a silent filter.",
        "",
    ]
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(stamp("\n".join(lines), "scripts/140_cis_mr_screen.py"), encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
