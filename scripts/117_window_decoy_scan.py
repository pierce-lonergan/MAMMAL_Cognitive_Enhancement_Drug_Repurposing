"""L4 psychoplastogen-window precision/FPR scan against a broad scaffold-matched decoy panel.

Systematic hardening of the window (#1 after the D5 pergolide/isoDMT audit): a larger panel of
marketed CNS drugs that SHARE the window's scaffold families (tryptamine/indole, ergoline,
di-substituted phenethylamine, N1-substituted indole) but are NOT 5-HT2A-agonist psychoplastogens.
Every such decoy MUST be window-NEGATIVE; every known psychedelic MUST stay window-POSITIVE. Any
new false positive is reported (then fixed only with a ledger-checked, non-overfit veto).

Drug-class labels are textbook-established pharmacology; SMILES are fetched from PubChem (verifiable,
never hand-written). Writes data/raw/l4_decoy_panel.csv (reproducible cache) +
reports/pipeline/l4_window_precision_v1.md. CPU + network.
"""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("l4_decoy_scan")
ROOT = Path(__file__).resolve().parents[1]
OUT_CSV = ROOT / "data" / "raw" / "l4_decoy_panel.csv"
REPORT = ROOT / "reports" / "pipeline" / "l4_window_precision_v1.md"

# (name, scaffold family, pharmacology class) - all NON-psychoplastogens -> must be window-NEGATIVE
DECOYS = [
    # triptans (5-HT1B/1D agonists; tryptamine core)
    ("sumatriptan", "tryptamine", "triptan 5-HT1B/1D"),
    ("rizatriptan", "tryptamine", "triptan 5-HT1B/1D"),
    ("eletriptan", "tryptamine", "triptan 5-HT1B/1D"),
    ("almotriptan", "tryptamine", "triptan 5-HT1B/1D"),
    ("naratriptan", "tryptamine", "triptan 5-HT1B/1D"),
    ("frovatriptan", "tryptamine", "triptan 5-HT1B/1D"),
    # 5-HT3 antagonists (indole/carbazole)
    ("ondansetron", "indole", "5-HT3 antagonist"),
    ("granisetron", "indole", "5-HT3 antagonist"),
    # melatonergics (tryptamine/indole)
    ("melatonin", "tryptamine", "MT1/MT2 agonist"),
    ("ramelteon", "indole", "MT1/MT2 agonist"),
    ("agomelatine", "naphthalene", "MT agonist/5-HT2C antag"),
    ("tasimelteon", "benzofuran", "MT1/MT2 agonist"),
    # indole alkaloids / nootropics (non-psychedelic)
    ("yohimbine", "indole", "alpha2 antagonist"),
    ("reserpine", "indole", "VMAT inhibitor"),
    ("vinpocetine", "indole", "PDE1/vasodilator nootropic"),
    ("tadalafil", "indole", "PDE5 inhibitor"),
    # precursors / endogenous
    ("tryptophan", "tryptamine", "amino acid"),
    ("oxitriptan", "tryptamine", "5-HTP precursor"),
    ("serotonin", "tryptamine", "endogenous (impermeant)"),
    # ergolines (dopaminergic / antimigraine / antagonist)
    ("bromocriptine", "ergoline", "D2 agonist"),
    ("cabergoline", "ergoline", "D2 agonist"),
    ("pergolide", "ergoline", "D2 agonist (clavine)"),
    ("methysergide", "ergoline", "5-HT2 antagonist"),
    ("metergoline", "ergoline", "5-HT antagonist"),
    ("ergotamine", "ergoline", "5-HT1/alpha agonist"),
    ("dihydroergotamine", "ergoline", "antimigraine"),
    ("nicergoline", "ergoline", "alpha blocker"),
    # di-substituted phenethylamines that are NOT psychedelic (2+ aromatic OMe)
    ("verapamil", "phenethylamine", "Ca channel blocker (3,4-diOMe)"),
    ("trimetazidine", "phenethylamine", "anti-anginal (3,4,5-triOMe)"),
    ("methoxamine", "phenethylamine", "alpha1 agonist (2,5-diOMe)"),
    # N1-substituted indoles (stress the new isoDMT detector's FPR)
    ("indomethacin", "indole", "NSAID (N1-aroyl indole)"),
    ("frovatriptan", "tryptamine", "triptan (dup-guard)"),
]
# known psychoplastogens - MUST stay window-POSITIVE (recall regression)
POSITIVES = [
    "lysergide", "psilocin", "dimethyltryptamine", "5-methoxy-N,N-dimethyltryptamine",
    "mescaline", "ibogaine",
]


def _smiles(name: str) -> str | None:
    from mammal_repurposing.fetchers.pubchem import fetch_smiles
    r = fetch_smiles(name)
    return r.get("smiles") if r else None


def main() -> int:
    from mammal_repurposing.engine.psychoplastogen import psychoplastogen_window
    seen, rows = set(), []
    for name, fam, klass in DECOYS:
        if name in seen:
            continue
        seen.add(name)
        smi = _smiles(name)
        if not smi:
            L.warning("no SMILES for %s", name); continue
        c = psychoplastogen_window(smi)
        rows.append(dict(name=name, role="decoy", family=fam, klass=klass, smiles=smi,
                         scaffold=c.scaffold, window=bool(c.window), clogp=c.clogp,
                         tpsa=c.tpsa, hbd=c.hbd, expected=False))
    for name in POSITIVES:
        smi = _smiles(name)
        if not smi:
            L.warning("no SMILES for positive %s", name); continue
        c = psychoplastogen_window(smi)
        rows.append(dict(name=name, role="positive", family="psychedelic", klass="psychoplastogen",
                         smiles=smi, scaffold=c.scaffold, window=bool(c.window), clogp=c.clogp,
                         tpsa=c.tpsa, hbd=c.hbd, expected=True))
    df = pd.DataFrame(rows)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False)

    dec = df[df.role == "decoy"]; pos = df[df.role == "positive"]
    fps = dec[dec.window]; fns = pos[~pos.window]
    fpr = len(fps) / len(dec) if len(dec) else float("nan")
    recall = pos.window.mean() if len(pos) else float("nan")
    L.info("Decoys %d | FALSE POSITIVES %d (FPR %.2f)", len(dec), len(fps), fpr)
    if len(fps):
        for _, r in fps.iterrows():
            L.info("  FP: %s (%s, %s) scaffold=%s TPSA=%s", r["name"], r.family, r.klass,
                   r.scaffold, r.tpsa)
    L.info("Positives %d | recall %.2f | MISSES %d", len(pos), recall, len(fns))
    if len(fns):
        for _, r in fns.iterrows():
            L.info("  MISS: %s scaffold=%s", r["name"], r.scaffold)

    write_report(df, dec, pos, fps, fns, fpr, recall)
    return 0


def write_report(df, dec, pos, fps, fns, fpr, recall) -> None:
    Ls = ["# L4 psychoplastogen window - precision / FPR scan", "",
          "Systematic hardening of the L4 window against a broad scaffold-matched decoy panel "
          "(marketed CNS drugs sharing the window's scaffold families but NOT 5-HT2A-agonist "
          "psychoplastogens). Reproduced by `scripts/117_window_decoy_scan.py`; SMILES cached in "
          "`data/raw/l4_decoy_panel.csv` (PubChem-sourced, never hand-written).", "",
          f"## Result: FPR **{len(fps)}/{len(dec)} = {fpr:.2f}** | positive recall "
          f"**{int(pos.window.sum())}/{len(pos)} = {recall:.2f}**", ""]
    if len(fps) == 0:
        Ls.append("No false positives: every scaffold-matched non-psychoplastogen decoy is "
                  "window-NEGATIVE. The serotonergic-scaffold + permeability gate, with the "
                  "triptan-sulfonamide/carbamate and clavine-thioether vetoes, cleanly excludes "
                  "the look-alike drug classes.")
    else:
        Ls.append("Remaining false positives (candidates for a ledger-checked veto):")
        for _, r in fps.iterrows():
            Ls.append(f"- **{r['name']}** ({r.klass}, {r.family}) - scaffold {r.scaffold}, "
                      f"TPSA {r.tpsa}, clogP {r.clogp}")
    Ls += ["", "## Decoy families covered", ""]
    for fam in sorted(dec.family.unique()):
        sub = dec[dec.family == fam]
        Ls.append(f"- {fam}: {len(sub)} decoys, {int(sub.window.sum())} FP")
    Ls += ["", "## Positive recall regression (must stay window-positive)", ""]
    for _, r in pos.iterrows():
        Ls.append(f"- {r['name']}: window={r.window} (scaffold {r.scaffold})")
    if len(fns):
        Ls += ["", "MISSES (recall gaps):"] + [f"- {r['name']}" for _, r in fns.iterrows()]
    Ls += ["", "## Honest scope", "",
           "This scan stresses PRECISION on serotonergic look-alikes only. The window remains "
           "deliberately blind to NON-serotonergic durable-plasticity classes (NMDA/dissociative, "
           "GABA-A neurosteroid, muscarinic) - that is a scope boundary, not a precision failure, "
           "and is the subject of the separate L4b research lane.", ""]
    REPORT.write_text("\n".join(Ls), encoding="utf-8")
    L.info("Wrote %s and %s", OUT_CSV, REPORT)


if __name__ == "__main__":
    raise SystemExit(main())
