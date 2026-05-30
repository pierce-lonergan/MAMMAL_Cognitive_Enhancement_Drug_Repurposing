"""Gap 7 — Prospective, mechanism-justified repurposing shortlist (capstone).

This is the deliverable the whole pipeline was building toward: a ranked set of
**approved drugs as repurposing hypotheses** for a chosen cognitive-impairment
disease, justified by the one signal the project proved is prognostic — the
drug's mechanism-class clinical track record IN THAT DISEASE (Gap 3 / Gap 6) —
and engaging a disease-relevant target.

It integrates every prior piece:
  * Gap 2 disease-conditioned class priors — is this drug's mechanism class a
    SUCCESS-track-record class in this disease?
  * the 31-target binding grid — does the drug engage a relevant target?
  * Gap 4 — flags engagement reliability (MAMMAL is allosteric-blind at
    transporter / allosteric targets, so binding there is uncertain).
  * Gap 5 — top picks get full GRADE evidence dossiers.

The headline output is the **NOVEL** column: an approved drug, in a mechanism
class with a positive pivotal-trial record in disease D, that engages a
D-relevant target but is **not currently indicated for D** — i.e. a
mechanism-justified repurposing hypothesis (roflumilast→cognition, buspirone→
CIAS, …). These are hypotheses worth evaluation, not predicted cures — the
Roberts-2020 ceiling and the binding-reliability caveats still hold.

numpy/pandas only.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# Curated current/approved indication(s) per named drug (lower-case). Used ONLY
# to flag novelty: a drug is a NOVEL repurposing hypothesis for disease D if D
# is not already in its indication set. Conservative — when unsure, list the
# known use so we DON'T over-claim novelty.
KNOWN_INDICATIONS: dict[str, set[str]] = {
    "donepezil": {"AD"}, "rivastigmine": {"AD", "PD"}, "galantamine": {"AD"},
    "huperzine a": set(),                      # supplement; not an approved AD drug
    "memantine": {"AD"}, "tacrine": {"AD"},
    "methylphenidate": {"ADHD", "narcolepsy"}, "d-amphetamine": {"ADHD"},
    "lisdexamfetamine": {"ADHD"}, "dextroamphetamine": {"ADHD"},
    "modafinil": {"narcolepsy"}, "armodafinil": {"narcolepsy"},
    "atomoxetine": {"ADHD"}, "guanfacine": {"ADHD"}, "clonidine": {"ADHD"},
    "pitolisant": {"narcolepsy"},
    "encenicline": set(), "tc-5619": set(), "abt-126": set(),  # failed clinical agents
    "suvorexant": {"insomnia"}, "lemborexant": {"insomnia"},
    "roflumilast": {"COPD"},                   # PDE4 — repurposing target for cognition
    "rolipram": set(), "bpn14770": {"FXS"}, "zatolmilast": {"FXS"},
    "buspirone": {"anxiety"},                  # 5-HT1A — repurposing target for CIAS
    "tandospirone": {"anxiety"},
    "vortioxetine": {"MDD"}, "fluvoxamine": {"OCD", "MDD"},
    "fluoxetine": {"MDD"}, "sertraline": {"MDD"}, "paroxetine": {"MDD"},
    "citalopram": {"MDD"}, "venlafaxine": {"MDD"}, "duloxetine": {"MDD"},
    "aripiprazole": {"schizophrenia", "MDD"}, "lurasidone": {"schizophrenia"},
    "xanomeline": {"CIAS"}, "pimavanserin": {"PD"},
    "blarcamesine": set(), "pridopidine": set(),
    "lamotrigine": {"epilepsy"}, "levetiracetam": {"epilepsy"},
    "topiramate": {"epilepsy"}, "retigabine": {"epilepsy"},
    "riluzole": {"ALS"}, "troriluzole": set(), "memantine_xr": {"AD"},
    "ivabradine": {"angina"}, "semaglutide": {"diabetes"},
    "liraglutide": {"diabetes"}, "clemastine": {"allergy"},
    "metformin": {"diabetes"}, "minocycline": {"infection"},
}


# Curated supplement of APPROVED drugs at cognition-relevant targets — especially
# the newly-added panel targets (CHRM1/CHRM4 = M1/M4, HTR6, HTR4, GlyT1) that the
# original 298-compound screening library predates. Each maps a real, recognisable
# approved drug to its AUTHORITATIVE primary cognition target (UniProt), so the
# repurposing shortlist can surface the muscarinic / 5-HT4 / GlyT1 hypotheses the
# 31-target panel expansion enabled. (name, target_uniprot)
REPURPOSING_SUPPLEMENT: list[tuple[str, str]] = [
    ("cevimeline", "P11229"),     # M1/M3 agonist, approved Sjögren's dry mouth
    ("pilocarpine", "P11229"),    # muscarinic agonist, approved glaucoma/xerostomia
    ("xanomeline", "P11229"),     # M1/M4 — now CIAS (known, for contrast)
    ("emraclidine", "P08173"),    # M4 PAM — CIAS (investigational)
    ("buspirone", "P08908"),      # 5-HT1A partial agonist, approved anxiety
    ("tandospirone", "P08908"),   # 5-HT1A partial agonist (anxiety, Asia)
    ("roflumilast", "Q08499"),    # PDE4 inhibitor, approved COPD
    ("prucalopride", "Q13639"),   # 5-HT4 agonist, approved constipation
    ("sarcosine", "P48067"),      # GlyT1 inhibitor (investigational adjunct)
    ("idalopirdine", "P50406"),   # 5-HT6 antagonist (failed AD — for contrast)
]

# add their current indications to the novelty map
KNOWN_INDICATIONS.update({
    "cevimeline": {"Sjogren"}, "pilocarpine": {"glaucoma"},
    "emraclidine": set(), "prucalopride": {"constipation"},
    "sarcosine": set(), "idalopirdine": set(),
})


@dataclass
class RepurposingCandidate:
    compound: str
    indication: str
    mechanism_class: str
    target_gene: str
    target_uniprot: str
    class_g: float                 # disease-conditioned class prior mean
    class_verdict: str
    engagement: float              # within-target binding percentile [0,1]
    engagement_reliable: bool      # False at allosteric/transporter targets (Gap 4)
    score: float
    novel: bool                    # not currently indicated for this disease
    current_use: str
    liabilities: list[str] = field(default_factory=list)

    @property
    def rationale(self) -> str:
        nov = "NOVEL repurposing" if self.novel else "known/standard-of-care"
        rel = "" if self.engagement_reliable else " (engagement uncertain — allosteric/transporter)"
        return (f"{self.mechanism_class} is a SUCCESS-track-record class in "
                f"{self.indication} (class g≈{self.class_g:+.2f}); {self.compound} "
                f"engages {self.target_gene} at {self.engagement:.0%} percentile{rel}. "
                f"Currently used for: {self.current_use or 'n/a'} → {nov}.")


# For novelty, a disease's "already used" set includes its parent indication:
# CIAS (cognitive impairment in schizophrenia) is covered by any schizophrenia
# use, AD-mod-sev by AD, etc. A drug used for the parent is NOT novel.
DISEASE_NOVELTY_ALIASES: dict[str, set[str]] = {
    "CIAS": {"CIAS", "schizophrenia"},
    "AD": {"AD", "AD-mod-sev", "dementia"},
    "FXS": {"FXS"},
    "ADHD": {"ADHD"},
}


# Targets where MAMMAL's sequence-only binding is unreliable (Gap 4: allosteric /
# transporter). Engagement at these is flagged uncertain.
ALLOSTERIC_BLIND_UNIPROTS = {
    "P36544", "Q08499", "O76083", "P42261", "P42262", "P42263", "P48058",
    "Q01959", "P23975", "P50406", "P08908", "Q13639", "P48067",
}


def named_drugs_with_supplement(named: pd.DataFrame) -> pd.DataFrame:
    """Append the curated approved-drug supplement (mapped to the new panel
    targets) to the screening library's named drugs, de-duplicated by name."""
    sup = pd.DataFrame(REPURPOSING_SUPPLEMENT, columns=["name", "expected_top_target"])
    have = set(named["name"].str.lower())
    sup = sup[~sup["name"].str.lower().isin(have)]
    return pd.concat([named[["name", "expected_top_target"]], sup],
                     ignore_index=True)


def build_repurposing_shortlist(
    disease: str,
    named_drugs: pd.DataFrame,
    disease_priors: dict,
    grid: pd.DataFrame,
    target_class_map: dict[str, str],
    *,
    liability_df: pd.DataFrame | None = None,
    panel_uniprots: set[str] | None = None,
    success_only: bool = True,
) -> list[RepurposingCandidate]:
    """Rank named/approved drugs as repurposing hypotheses for `disease`.

    `named_drugs`: DataFrame [name, expected_top_target] (the curated drug list).
    `disease_priors`: {clean_class: DiseaseClassPrior} from disease_reframe.
    `grid`: the (compound × target) binding grid with within-target percentile.
    `target_class_map`: {uniprot: clean mechanism class}.

    Score = class_prior_g × (0.5 + 0.5·engagement): the disease-class prognostic
    prior dominates (the validated signal), engagement is a secondary modifier.
    """
    panel_uniprots = panel_uniprots or set()
    # within-target binding percentile per (compound, target)
    g = grid.copy()
    g["target_uniprot"] = g["target_uniprot"].astype(str)
    if "binding_percentile" in g.columns:
        pct = g.set_index([g["compound"].str.lower() if "compound" in g else
                           g["compound_name"].str.lower(), "target_uniprot"])["binding_percentile"]
        bind = {(c, t): float(v) for (c, t), v in pct.items()}
    else:
        col = "predicted_pkd"
        g["pct"] = g.groupby("target_uniprot")[col].rank(pct=True)
        ck = g["compound_name"].str.lower() if "compound_name" in g else g["compound"].str.lower()
        bind = {(c, t): float(v) for c, t, v in zip(ck, g["target_uniprot"], g["pct"])}

    out: list[RepurposingCandidate] = []
    for _, r in named_drugs.iterrows():
        name = str(r["name"])
        uni = str(r.get("expected_top_target", "") or "")
        if not uni or uni not in target_class_map:
            continue
        cls = target_class_map[uni]
        prior = disease_priors.get(cls)
        if prior is None:
            continue
        if success_only and prior.verdict != "SUCCESS":
            continue
        eng = bind.get((name.lower(), uni), float("nan"))
        if not np.isfinite(eng):
            eng = 0.5  # unknown engagement → neutral
        reliable = uni not in ALLOSTERIC_BLIND_UNIPROTS
        score = prior.mean * (0.5 + 0.5 * eng)
        known = KNOWN_INDICATIONS.get(name.lower(), None)
        # novel if we know its uses and none of them cover this disease (incl.
        # its parent indication — schizophrenia covers CIAS, etc.)
        disease_cover = DISEASE_NOVELTY_ALIASES.get(disease, {disease})
        novel = (known is not None) and not (known & disease_cover)
        current = ", ".join(sorted(known)) if known else ("unknown" if known is None else "—")
        liab = _liab(name, liability_df, panel_uniprots)
        gene = uni
        if "target_gene" in grid.columns:
            gg = grid[grid["target_uniprot"].astype(str) == uni]["target_gene"]
            if len(gg):
                gene = str(gg.iloc[0])
        out.append(RepurposingCandidate(
            compound=name, indication=disease, mechanism_class=cls,
            target_gene=gene, target_uniprot=uni, class_g=prior.mean,
            class_verdict=prior.verdict, engagement=eng,
            engagement_reliable=reliable, score=score, novel=novel,
            current_use=current, liabilities=liab,
        ))
    out.sort(key=lambda c: c.score, reverse=True)
    return out


def _liab(name: str, liability_df, panel_uniprots, *, thresh: float = 6.5,
          top: int = 3) -> list[str]:
    if liability_df is None or not len(liability_df):
        return []
    sub = liability_df[liability_df["compound_name"].str.lower() == name.lower()]
    sub = sub[~sub["target_uniprot"].astype(str).isin(panel_uniprots)]
    sub = sub[sub["predicted_pkd"] >= thresh].sort_values("predicted_pkd", ascending=False).head(top)
    return [str(r.get("target_gene", r["target_uniprot"])) for _, r in sub.iterrows()]


def availability() -> dict:
    return {"available": True, "n_known_indications": len(KNOWN_INDICATIONS),
            "scores": "class_prior × (0.5 + 0.5·engagement); SUCCESS classes only"}
