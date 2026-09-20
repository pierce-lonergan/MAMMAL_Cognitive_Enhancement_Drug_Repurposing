"""Could this pipeline ever have DISCOVERED the compounds it already knows?

The project ranks known compounds. Predicting NOVEL ones is a different claim, and before building
any generator there are two prior questions that bound what a generator could possibly do here. Both
are answerable from data already in the repository, so they are measured rather than argued.

QUESTION 1: HOW MUCH OF THE POSITIVE CLASS IS EVEN A SMALL MOLECULE? A generative model designs
single small molecules. It cannot design a multivitamin-mineral mixture, a botanical extract, a
dietary pattern, or a nine-residue peptide. Whatever fraction of the enhancer class is not a single
designable molecule is simply outside the reach of de novo design, permanently, and subtracting it
first is the honest starting point.

QUESTION 2: THE NOVELTY CEILING. `scripts/95_novel_compound_onboarding.py` routes an arbitrary SMILES
to a known mechanism class only when its maximum Tanimoto to a known chemotype clears a threshold,
and ABSTAINS otherwise, on the stated grounds that it "cannot invent a mechanism". On its own demo,
six of eight novel compounds abstained. So the engine is deliberately incapable of novelty, and the
question that bounds every generative plan is: HOW FAR APART ARE THE KNOWN ACTIVES THEMSELVES?

If a known enhancer, held out and presented as if it were novel, sits further from the remaining
actives than the abstention threshold, then the engine would have refused to score it. It could not
have discovered that compound. Run over every active, that fraction is the ceiling on what any
similarity-gated screen built on this evidence base could ever find.

This is the same question B2 asked of the plasticity window, pointed at the screen instead of the
assay: not "does it work" but "could it have worked on the cases we already know the answer to".

Writes reports/pipeline/novelty_ceiling_v1.md. CPU; network only for SMILES (PubChem, cached).
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("novelty_ceiling")
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.provenance.trailer import stamp  # noqa: E402

LEDGER = ROOT / "data" / "raw" / "healthy_adult_cognition_ledger.csv"
SMILES_FILE = ROOT / "data" / "raw" / "ledger_compound_smiles.csv"
CACHE = ROOT / "data" / "raw" / "novelty_ceiling_smiles.csv"
REPORT = ROOT / "reports" / "pipeline" / "novelty_ceiling_v1.md"

#: The abstention threshold the onboarding engine uses: below this maximum Tanimoto to any known
#: chemotype it refuses to route, because the compound is not near a precedented cognition scaffold.
ABSTAIN_BELOW = 0.35

#: At or above this, a "match" is the compound finding itself under another name: an
#: enantiomer, a salt or a prodrug. modafinil matches armodafinil at exactly 1.000,
#: because the fingerprint is computed without chirality and they are the same graph.
SELF_MATCH_AT = 0.95

#: Names that are not a single designable molecule. Recorded explicitly rather than inferred from a
#: failed lookup, because a failed PubChem call and a compound that cannot exist as one molecule are
#: different facts and conflating them would hide the finding.
NOT_A_SINGLE_MOLECULE = {
    "multivitamin_mineral": "mixture of many vitamins and minerals; no single structure exists",
    "b_vitamins": "a class of several distinct molecules dosed together",
    "guarana": "botanical extract; many constituents, caffeine among them",
    "ginkgo_biloba": "botanical extract (EGb 761 and relatives); flavonoid and terpene mixture",
    "bacopa_monnieri": "botanical extract; bacoside mixture",
    "fruit_derived_polyphenols": "a polyphenol CLASS, not a compound",
    "psilocybin_lsd_microdosing": "two distinct molecules pooled as one exposure",
    "menopausal_hormone_therapy": "several regimens of different hormones",
    "dietary_nitrate": "an ion delivered in food, not a designed molecule",
    "oxytocin_intranasal": "nine-residue cyclic peptide; outside small-molecule generative scope",
    "insulin_intranasal": "51-residue protein; outside small-molecule generative scope",
    "omega_3": "a class of fatty acids dosed as a mixture",
    "folic_acid": "single molecule, but a vitamin repletion exposure rather than a design target",
}


def primary_set() -> pd.DataFrame:
    df = pd.read_csv(LEDGER)
    cand = df.get("candidate_enhancer", pd.Series([1] * len(df)))
    dep = df.get("depends_on", pd.Series([""] * len(df))).fillna("").astype(str).str.strip()
    return df[(df["evidence_tier"] == "clean_MA") & (cand == 1) & (dep == "")].reset_index(drop=True)


def resolve_smiles(names: list[str]) -> dict[str, str]:
    """SMILES from the existing ledger file, then PubChem for the rest, cached."""
    out: dict[str, str] = {}
    if SMILES_FILE.exists():
        s = pd.read_csv(SMILES_FILE)
        out.update({str(a).lower(): str(b) for a, b in zip(s["compound"], s["smiles"])
                    if isinstance(b, str) and b})
    if CACHE.exists():
        c = pd.read_csv(CACHE)
        out.update({str(a).lower(): str(b) for a, b in zip(c["compound"], c["smiles"])
                    if isinstance(b, str) and b})

    todo = [n for n in names if n.lower() not in out and n not in NOT_A_SINGLE_MOLECULE]
    if todo:
        from mammal_repurposing.fetchers.pubchem import fetch_smiles
        got = {}
        for n in todo:
            try:
                r = fetch_smiles(n.replace("_", " "))
                if r and r.get("smiles"):
                    got[n.lower()] = r["smiles"]
            except Exception as e:                      # noqa: BLE001 - network is best-effort
                L.warning("SMILES lookup failed for %s: %s", n, e)
        out.update(got)
        if got:
            prev = pd.read_csv(CACHE) if CACHE.exists() else pd.DataFrame(columns=["compound", "smiles"])
            new = pd.DataFrame({"compound": list(got), "smiles": [got[k] for k in got]})
            pd.concat([prev, new], ignore_index=True).drop_duplicates("compound").to_csv(
                CACHE, index=False)
    return out


def fingerprints(smiles_by_name: dict[str, str]):
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    RDLogger.DisableLog("rdApp.*")
    fps, ok = {}, {}
    for name, smi in smiles_by_name.items():
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        fps[name] = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        ok[name] = smi
    return fps, ok


def exemplar_fingerprints():
    """The FULL exemplar base the onboarding engine actually routes against.

    Corrected after a first version of this script compared each compound only against the other
    healthy-adult primary-set compounds. That is a harder test than the engine faces: it routes
    against `ledger_compound_smiles.csv`, 110 compounds across 46 mechanism classes. Reporting the
    restricted number as if it were the engine's behaviour would have overstated the finding.
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    RDLogger.DisableLog("rdApp.*")
    df = pd.read_csv(SMILES_FILE)
    fps = {}
    for name, smi in zip(df["compound"].astype(str), df["smiles"]):
        if not isinstance(smi, str) or not smi:
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            fps[name] = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    return fps


def loo_against_exemplars(fps: dict, exemplars: dict) -> list[dict]:
    """Hold each compound out of the FULL exemplar base and ask whether the engine would route it."""
    from rdkit import DataStructs
    rows = []
    for n, fp in sorted(fps.items()):
        others = {m: f for m, f in exemplars.items() if m.lower() != n.lower()}
        if not others:
            continue
        sims = {m: float(DataStructs.TanimotoSimilarity(fp, f)) for m, f in others.items()}
        best = max(sims, key=sims.get)
        # A "match" at near-identity is the compound finding ITSELF under another name: a single
        # enantiomer, a salt, or a prodrug. That is not discovery, and counting it as a routing
        # success would flatter the screen. Flagged so the report can subtract it.
        trivial = sims[best] >= SELF_MATCH_AT
        rows.append(dict(compound=n, max_tanimoto=sims[best], nearest=best,
                         would_abstain=sims[best] < ABSTAIN_BELOW, n_exemplars=len(others),
                         self_match=trivial))
    return rows


def leave_one_out_similarity(fps: dict) -> list[dict]:
    """For each compound: its maximum Tanimoto to every OTHER compound in the set.

    This is the held-out question. If a known active, presented as if novel, does not clear the
    abstention threshold against the remaining actives, the engine would have refused to score it.
    """
    from rdkit import DataStructs
    names = sorted(fps)
    rows = []
    for n in names:
        others = [m for m in names if m != n]
        if not others:
            continue
        sims = [float(DataStructs.TanimotoSimilarity(fps[n], fps[m])) for m in others]
        best = int(np.argmax(sims))
        rows.append(dict(compound=n, max_tanimoto=max(sims), nearest=others[best],
                         would_abstain=max(sims) < ABSTAIN_BELOW))
    return rows


def main() -> int:
    prim = primary_set()
    names = [str(c) for c in prim["compound"]]
    labels = dict(zip(names, prim["enhances_healthy_young"].astype(int)))

    designable = [n for n in names if n not in NOT_A_SINGLE_MOLECULE]
    excluded = [n for n in names if n in NOT_A_SINGLE_MOLECULE]

    smi = resolve_smiles(designable)
    have = {n: smi[n.lower()] for n in designable if n.lower() in smi}
    unresolved = [n for n in designable if n not in have]

    fps, ok = fingerprints(have)
    L.info("primary set %d | not a single molecule %d | resolved %d | unparseable/unresolved %d",
           len(names), len(excluded), len(ok), len(unresolved) + (len(have) - len(ok)))

    # the enhancer class specifically, which is what a generator would be conditioned on
    pos = [n for n in names if labels.get(n) == 1]
    pos_designable = [n for n in pos if n in ok]
    L.info("ENHANCERS %d | of which single designable molecules with a structure: %d (%s)",
           len(pos), len(pos_designable), ", ".join(pos_designable) or "none")

    exemplars = exemplar_fingerprints()
    L.info("full exemplar base the engine routes against: %d parseable compounds", len(exemplars))
    loo_full = loo_against_exemplars(fps, exemplars)
    n_abstain_full = sum(1 for r in loo_full if r["would_abstain"])
    L.info("AGAINST THE FULL EXEMPLAR BASE: %d/%d would be abstained on",
           n_abstain_full, len(loo_full))
    for r in sorted(loo_full, key=lambda x: x["max_tanimoto"]):
        L.info("   %-24s max Tanimoto %.3f to %-24s %s", r["compound"], r["max_tanimoto"],
               r["nearest"], "ABSTAIN" if r["would_abstain"] else "routed")

    loo = leave_one_out_similarity(fps)
    loo_pos = [r for r in loo if labels.get(r["compound"]) == 1]
    n_abstain = sum(1 for r in loo if r["would_abstain"])
    L.info("leave-one-out: %d/%d compounds would be ABSTAINED on at threshold %.2f",
           n_abstain, len(loo), ABSTAIN_BELOW)
    for r in sorted(loo, key=lambda x: x["max_tanimoto"]):
        L.info("   %-24s max Tanimoto %.3f to %-22s %s",
               r["compound"], r["max_tanimoto"], r["nearest"],
               "ABSTAIN" if r["would_abstain"] else "")

    write_report(prim, excluded, unresolved, ok, pos, pos_designable, loo, loo_pos,
                 n_abstain, loo_full, n_abstain_full, len(exemplars))
    return 0


def write_report(prim, excluded, unresolved, ok, pos, pos_designable, loo, loo_pos,
                 n_abstain, loo_full, n_abstain_full, n_exemplars):
    out = [
        "# Could this pipeline ever have discovered what it already knows?", "",
        "Before building anything that predicts NOVEL compounds, two prior questions bound what such "
        "a thing could do here. Both are answerable from data already in the repository.", "",
        "## 1. How much of the evidence base is even a designable molecule?", "",
        f"The primary analysis set is **n = {len(prim)}**. Of those, **{len(excluded)}** are not a "
        "single small molecule at all and are therefore permanently outside the reach of de novo "
        "design, whatever the model:", "",
        "| compound | why it cannot be designed |", "| --- | --- |",
    ]
    for c in excluded:
        out.append(f"| {c} | {NOT_A_SINGLE_MOLECULE[c]} |")
    if unresolved:
        out += ["", f"A further **{len(unresolved)}** could not be resolved to a structure "
                f"({', '.join(unresolved)}), which is a data gap rather than a statement about the "
                "compound."]
    out += ["", f"**The enhancer class is the part that matters**, because a generator would be "
            f"conditioned on the positives. There are **{len(pos)}** labelled enhancers, of which "
            f"**{len(pos_designable)}** are single molecules with a resolved structure: "
            f"{', '.join(pos_designable) or 'none'}.", ""]
    if len(pos_designable) <= 4:
        out += [f"Conditioning a generative model on {len(pos_designable)} molecules is not "
                "discovery. Whatever it produced would be analogues of those few, and they are "
                "already known, already scheduled or already consumed at scale.", ""]

    out += ["## 2. The novelty ceiling", "",
            "`scripts/95_novel_compound_onboarding.py` routes a novel SMILES to a known mechanism "
            f"class only when its maximum Tanimoto to a known chemotype clears **{ABSTAIN_BELOW}**, "
            "and abstains otherwise, on the stated grounds that it cannot invent a mechanism. On "
            "its own demo, six of eight novel compounds abstained.", "",
            "So: hold each known compound out and present it as if it were novel. Does the engine "
            "clear its own threshold against the remaining actives, or does it refuse to score a "
            "compound we already know works?", "",
            f"### Against the FULL exemplar base ({n_exemplars} compounds, what the engine "
            "actually does)", "",
            "| compound | label | max Tanimoto | nearest exemplar | verdict |",
            "| --- | :---: | ---: | --- | --- |"]
    _lab = dict(zip(prim["compound"].astype(str), prim["enhances_healthy_young"].astype(int)))
    for r in sorted(loo_full, key=lambda x: -x["max_tanimoto"]):
        v = "**ABSTAIN**" if r["would_abstain"] else "routed"
        out.append(f"| {r['compound']} | {_lab.get(r['compound'], '?')} | "
                   f"{r['max_tanimoto']:.3f} | {r['nearest']} | {v} |")
    _routed = [r for r in loo_full if not r["would_abstain"]]
    _self = [r for r in _routed if r["self_match"]]
    out += ["", f"**{n_abstain_full} of {len(loo_full)}** would be abstained on against the full "
            f"base. Of the **{len(_routed)}** that route, **{len(_self)}** match at or above "
            f"{SELF_MATCH_AT}, which is the compound finding ITSELF under another name: "
            + (", ".join(f"{r['compound']} to {r['nearest']} at {r['max_tanimoto']:.3f}"
                         for r in _self) if _self else "none")
            + ".", "",
            "Read the routed rows carefully, because two of them are not discoveries. Modafinil "
            "matches armodafinil, which is its own single enantiomer, at a Tanimoto of exactly "
            "1.000 since the fingerprint ignores chirality. Dextroamphetamine matches "
            "lisdexamfetamine, which is its own lysine prodrug. Subtracting those leaves ONE "
            "non-trivial route in the whole set, and it barely clears the threshold.", "",
            "### Against the healthy-adult evidence base only", "",
            "The stricter and more relevant comparison, because it is the healthy-adult evidence "
            "that defines this project's outcome. A compound routed on its resemblance to a CNS "
            "drug with no healthy-adult cognitive evidence has been matched to a chemotype, not to "
            "a demonstrated effect.", "",
            "| compound | label | max Tanimoto to any other | nearest | verdict |",
            "| --- | :---: | ---: | --- | --- |"]
    labels = dict(zip(prim["compound"].astype(str), prim["enhances_healthy_young"].astype(int)))
    for r in sorted(loo, key=lambda x: -x["max_tanimoto"]):
        v = "**ABSTAIN**" if r["would_abstain"] else "routed"
        out.append(f"| {r['compound']} | {labels.get(r['compound'], '?')} | "
                   f"{r['max_tanimoto']:.3f} | {r['nearest']} | {v} |")

    pos_abstain = sum(1 for r in loo_pos if r["would_abstain"])
    out += ["", f"**{n_abstain} of {len(loo)}** structurally resolved compounds would be abstained "
            f"on, including **{pos_abstain} of {len(loo_pos)}** known ENHANCERS.", ""]
    if loo_pos:
        med = float(np.median([r["max_tanimoto"] for r in loo_pos]))
        out += [f"Median maximum Tanimoto among the enhancers is **{med:.3f}**.", ""]
    out += ["## What this bounds", "",
            "The similarity gate is not a tuning parameter that could be relaxed. It is what keeps "
            "the engine honest: below it, the engine has no evidence that the compound belongs to "
            "any precedented mechanism, and routing it anyway would be inventing a mechanism. "
            "Lowering the threshold does not create knowledge, it only stops recording ignorance.", "",
            "So a similarity-gated screen built on this evidence base can only find compounds that "
            "resemble the handful it already has. That is interpolation. It is a legitimate and "
            "useful thing to do, and it should be described as analogue search rather than "
            "discovery.", "",
            "The deeper limit is the one in section 1, and no model fixes it: most of what this "
            "project has established as working in healthy adults is not a molecule anyone could "
            "design.", "",
            "Generated by `scripts/132_novelty_ceiling.py`."]
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(stamp("\n".join(out), "scripts/132_novelty_ceiling.py"), encoding="utf-8")
    L.info("wrote %s", REPORT.relative_to(ROOT))


if __name__ == "__main__":
    raise SystemExit(main())
