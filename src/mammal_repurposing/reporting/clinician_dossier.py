"""Gap 5 — Clinician-legible evidence dossier.

The pipeline already speaks the right currency — Hedges' g with credible
intervals — but a clinician does not read 90 markdown reports. This module
distils, per (compound, indication), a **one-page evidence card** in the
language a doctor actually uses:

  * predicted cognition effect size (Hedges' g) + credible interval
  * a **GRADE-style** evidence-quality rating (HIGH / MODERATE / LOW / VERY LOW)
    with the explicit reasons it was up/down-graded (the Cochrane framework a
    Cambridge clinician will recognise on sight)
  * the **mechanism-class clinical track record** (the Gap-3 prognostic signal:
    did this class succeed or fail in pivotal trials, and on how many?)
  * predicted **off-target liability flags** mapped to clinical concerns
    (the 44-target off-target panel)
  * the **provenance trail** (which pipeline layers contributed)
  * explicit **failure-mode caveats** (MAMMAL allosteric blindness, the Roberts
    2020 ceiling, healthy-vs-disease indirectness, small-n)

GRADE is applied honestly: meta-analytic RCT evidence starts HIGH and is
down-graded for imprecision, inconsistency, indirectness, and model-reliability
risk — exactly as Cochrane would. The output is meant to be the single artifact
a doctor would call useful.

numpy/pandas only.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# Off-target gene -> recognisable clinical concern (the 44-target liability panel)
LIABILITY_CONCERNS: dict[str, str] = {
    "KCNH2": "hERG: QT prolongation / torsades risk",
    "KCNA5": "Kv1.5: atrial electrophysiology",
    "DRD2": "D2: extrapyramidal symptoms / hyperprolactinaemia",
    "DRD4": "D4: psychiatric off-target",
    "HTR2A": "5-HT2A: sedation / metabolic",
    "HTR2C": "5-HT2C: weight gain / appetite",
    "HTR7": "5-HT7: thermoregulation / circadian",
    "HRH1": "H1: sedation / weight gain",
    "ADRA1A": "α1A: orthostatic hypotension",
    "ADRA1B": "α1B: orthostatic hypotension",
    "ADRA1D": "α1D: orthostatic hypotension",
    "CHRM1": "M1: anticholinergic (cognition/confusion)",
    "CHRM3": "M3: anticholinergic (dry mouth, GI, urinary)",
    "OPRM1": "µ-opioid: dependence / respiratory depression",
    "CNR1": "CB1: psychiatric / dependence",
    "MAOA": "MAO-A: hypertensive / serotonergic interaction",
    "MAOB": "MAO-B: dopaminergic interaction",
    "GABRA1": "GABA-A α1: sedation / dependence",
    "GABRA5": "GABA-A α5: cognition / sedation",
    "TACR1": "NK1: emetic / CNS off-target",
    "NR3C1": "glucocorticoid receptor: endocrine",
    "ESR1": "oestrogen receptor: endocrine",
    "NTRK1": "TrkA: pain / off-target kinase",
}

GRADE_LEVELS = ["VERY LOW", "LOW", "MODERATE", "HIGH"]


@dataclass
class DossierCard:
    compound: str
    indication: str
    mechanism_class: str
    target_gene: str
    g: float
    g_ci_lo: float
    g_ci_hi: float
    grade: str
    grade_reasons: list[str] = field(default_factory=list)
    class_verdict: str = "UNKNOWN"          # SUCCESS / FAILURE / MIXED / UNKNOWN
    class_n_drugs: int = 0
    class_k_rcts: int = 0
    class_siblings: list[str] = field(default_factory=list)
    liabilities: list[tuple[str, str, float]] = field(default_factory=list)
    provenance: list[str] = field(default_factory=list)
    caveats: list[str] = field(default_factory=list)
    own_trial: bool = False


def grade_evidence(*, k_rcts: int, n_class_drugs: int, class_sd: float,
                   own_trial: bool, disease_direct: bool, ci_width: float,
                   binding_reliable: bool) -> tuple[str, list[str]]:
    """GRADE-style quality rating. RCT evidence starts HIGH; down-grade for the
    standard Cochrane domains. Returns (level, reasons)."""
    score = 3  # HIGH
    reasons: list[str] = []
    if k_rcts >= 3:
        reasons.append(f"RCT meta-analytic base (k={k_rcts}) → start HIGH")
    elif k_rcts >= 1:
        reasons.append(f"limited trial base (k={k_rcts})")
    else:
        score -= 1
        reasons.append("no direct trials for this class in this indication (−1)")
    # imprecision
    if not np.isfinite(ci_width) or ci_width > 0.40 or k_rcts < 2:
        score -= 1
        reasons.append("imprecision: wide CI / few RCTs (−1)")
    # inconsistency
    if np.isfinite(class_sd) and class_sd > 0.15:
        score -= 1
        reasons.append(f"inconsistency: heterogeneous class effects (sd={class_sd:.2f}) (−1)")
    # indirectness — predicted from class siblings, not the drug's own trial
    if not own_trial:
        score -= 1
        reasons.append("indirectness: effect predicted from mechanism-class siblings, "
                       "not this compound's own pivotal trial (−1)")
    if not disease_direct:
        score -= 1
        reasons.append("indirectness: healthy-adult extrapolation, not the disease "
                       "population (−1)")
    # model-reliability risk of bias — only relevant when the effect is a
    # model PREDICTION relying on binding; an approved drug's own pivotal trial
    # does not depend on MAMMAL's (un)reliable target-engagement score.
    if not binding_reliable and not own_trial:
        score -= 1
        reasons.append("risk of bias: effect relies on MAMMAL sequence-only "
                       "target-engagement, which is allosteric-blind here (−1)")
    score = max(0, min(3, score))
    return GRADE_LEVELS[score], reasons


def _liability_flags(compound: str, liability_df: pd.DataFrame,
                     panel_uniprots: set[str], *, pkd_threshold: float = 6.5,
                     top_n: int = 4) -> list[tuple[str, str, float]]:
    if liability_df is None or not len(liability_df):
        return []
    sub = liability_df[liability_df["compound_name"].str.lower() == compound.lower()]
    sub = sub[~sub["target_uniprot"].astype(str).isin(panel_uniprots)]
    sub = sub[sub["predicted_pkd"] >= pkd_threshold]
    sub = sub.sort_values("predicted_pkd", ascending=False).head(top_n)
    out = []
    for _, r in sub.iterrows():
        gene = str(r.get("target_gene", r["target_uniprot"]))
        concern = LIABILITY_CONCERNS.get(gene, f"{gene}: off-target engagement")
        out.append((gene, concern, float(r["predicted_pkd"])))
    return out


def build_dossier(compound: str, indication: str, *,
                  ledger: pd.DataFrame,
                  disease_priors: dict,
                  anchor_row: pd.Series | None = None,
                  ledger_row: pd.Series | None = None,
                  liability_df: pd.DataFrame | None = None,
                  panel_uniprots: set[str] | None = None,
                  binding_reliable: bool = True,
                  target_gene: str = "",
                  mechanism_class: str = "") -> DossierCard:
    """Assemble one (compound, indication) evidence card from the pipeline's
    real outputs. `disease_priors` is a {class: DiseaseClassPrior} dict from
    `validation.disease_reframe.build_disease_class_priors(indication, ...)`.
    Effect-size source priority: the compound's own indication-matched pivotal
    outcome (`ledger_row`) > a modulator-anchor trial (`anchor_row`, gives k +
    CI) > the mechanism-class prediction (`disease_priors`)."""
    panel_uniprots = panel_uniprots or set()
    prior = disease_priors.get(mechanism_class)
    class_sd = prior.sd if prior is not None else 0.15

    # effect size + CI, by source priority
    own_trial = (ledger_row is not None) or (anchor_row is not None)
    if ledger_row is not None:
        # indication-matched real pivotal outcome (point estimate)
        g = float(ledger_row.get("clinical_g", np.nan))
        # k + CI from the anchor table if the same drug is there, else derive
        if anchor_row is not None:
            ci_lo = float(anchor_row.get("CI_lo", g - 1.2816 * class_sd))
            ci_hi = float(anchor_row.get("CI_hi", g + 1.2816 * class_sd))
            k_rcts = int(anchor_row.get("k", 3))
        else:
            ci_lo, ci_hi = g - 1.2816 * class_sd, g + 1.2816 * class_sd
            k_rcts = prior.k_total if prior is not None else 3
    elif anchor_row is not None:
        g = float(anchor_row.get("pooled_g", np.nan))
        ci_lo = float(anchor_row.get("CI_lo", np.nan))
        ci_hi = float(anchor_row.get("CI_hi", np.nan))
        k_rcts = int(anchor_row.get("k", 1))
    elif prior is not None:
        g = prior.mean
        ci_lo = prior.mean - 1.2816 * prior.sd
        ci_hi = prior.mean + 1.2816 * prior.sd
        k_rcts = prior.k_total
    else:
        g, ci_lo, ci_hi, k_rcts = float("nan"), float("nan"), float("nan"), 0

    class_sd = prior.sd if prior is not None else float("nan")
    n_class = prior.n_drugs if prior is not None else 0
    verdict = prior.verdict if prior is not None else "UNKNOWN"
    siblings = [d.split("_")[0] for d in (prior.drugs if prior is not None else [])][:5]

    ci_width = (ci_hi - ci_lo) if (np.isfinite(ci_hi) and np.isfinite(ci_lo)) else float("nan")
    disease_direct = indication.lower() not in ("healthy", "hc")
    grade, reasons = grade_evidence(
        k_rcts=k_rcts, n_class_drugs=n_class, class_sd=class_sd,
        own_trial=own_trial, disease_direct=disease_direct,
        ci_width=ci_width, binding_reliable=binding_reliable)

    liabilities = _liability_flags(compound, liability_df, panel_uniprots)

    provenance = ["V6.A multi-head DTI (target engagement)",
                  "V6.B Cluster D θ̄ (cognition relevance)",
                  "V7 / disease-conditioned class prior (effect size)",
                  "Gap-3 mechanism-class track record (prognosis)"]
    if liabilities:
        provenance.append("44-target off-target liability panel")
    if not binding_reliable:
        provenance.append("Gap-4 allosteric LTR (binding-reliability flag)")

    caveats = []
    caveats.append("Predicted cognition effect is bounded by the Roberts 2020 ceiling "
                   "(healthy-adult SMD ≈ 0.2-0.5); disease effects can be larger.")
    if not own_trial:
        caveats.append("This is a mechanism-class prediction, NOT a per-compound "
                       "clinical result — the effect size is the class ceiling scaled "
                       "by predicted target engagement.")
    if not binding_reliable:
        caveats.append("MAMMAL's sequence-only binding score is structurally blind at "
                       "this (allosteric/transporter) target; engagement is uncertain "
                       "(see Gap-4 allosteric learn-to-rank).")
    if verdict == "FAILURE":
        caveats.append("⚠ This mechanism class has a NEGATIVE pivotal-trial track "
                       "record in this indication — proceed only with a hypothesis "
                       "that distinguishes this compound from its failed class-mates.")
    if liabilities:
        caveats.append("Off-target liability flags are MODEL-PREDICTED (MAMMAL DTI), "
                       "unvalidated — confirm against the drug's known safety profile.")

    return DossierCard(
        compound=compound, indication=indication, mechanism_class=mechanism_class,
        target_gene=target_gene, g=g, g_ci_lo=ci_lo, g_ci_hi=ci_hi,
        grade=grade, grade_reasons=reasons, class_verdict=verdict,
        class_n_drugs=n_class, class_k_rcts=k_rcts, class_siblings=siblings,
        liabilities=liabilities, provenance=provenance, caveats=caveats,
        own_trial=own_trial,
    )


def render_card_md(card: DossierCard) -> str:
    """One-page markdown dossier for a single (compound, indication)."""
    L: list[str] = []
    L.append(f"### {card.compound} — for {card.indication}")
    L.append("")
    grade_badge = {"HIGH": "🟢 HIGH", "MODERATE": "🟡 MODERATE",
                   "LOW": "🟠 LOW", "VERY LOW": "🔴 VERY LOW"}.get(card.grade, card.grade)
    g_str = (f"**{card.g:+.2f}** (90% CrI {card.g_ci_lo:+.2f} to {card.g_ci_hi:+.2f})"
             if np.isfinite(card.g) else "not estimable")
    L.append(f"| | |")
    L.append(f"|---|---|")
    L.append(f"| **Mechanism** | {card.mechanism_class}"
             f"{' at ' + card.target_gene if card.target_gene else ''} |")
    L.append(f"| **Predicted cognition effect** | Hedges' g = {g_str} |")
    L.append(f"| **Evidence quality (GRADE)** | {grade_badge} |")
    L.append(f"| **Mechanism-class track record** | {card.class_verdict} "
             f"({card.class_n_drugs} class drugs, {card.class_k_rcts} pooled RCTs) |")
    if card.class_siblings:
        L.append(f"| **Class exemplars** | {', '.join(card.class_siblings)} |")
    L.append(f"| **Effect basis** | "
             f"{'this compound’s own pivotal trial' if card.own_trial else 'mechanism-class prediction'} |")
    L.append("")
    L.append(f"**Why this GRADE rating:** " + "; ".join(card.grade_reasons) + ".")
    L.append("")
    if card.liabilities:
        L.append("**Predicted off-target liability flags** (model-based, unvalidated):")
        for gene, concern, pkd in card.liabilities:
            L.append(f"- `{gene}` (pKd≈{pkd:.1f}) — {concern}")
        L.append("")
    L.append("**Provenance:** " + " · ".join(card.provenance) + ".")
    L.append("")
    L.append("**Caveats & failure modes:**")
    for c in card.caveats:
        L.append(f"- {c}")
    L.append("")
    return "\n".join(L)


def availability() -> dict:
    return {"available": True, "grade_levels": GRADE_LEVELS,
            "n_liability_concerns": len(LIABILITY_CONCERNS)}
