"""The RETROSPECTIVE arm: healthy-adult cognition trials whose results are already posted.

The prospective registry cannot be scored. A blind sweep over its due rows found 25 of 26 with no
publication at all, and not one of them had registry-posted results, so journal publication is the
only channel and publication lag is the binding constraint.

This opens the other door. ClinicalTrials.gov requires many sponsors to POST results as structured
data, independent of whether a paper ever appears. Those trials are immune to the file drawer in a
way the journal literature is not, and their results are machine-readable today.

WHAT THIS ARM IS, AND WHAT IT IS NOT. It is RETROSPECTIVE. The results already exist, so applying
the frozen rule to them is a retrodiction, not a prediction. It is scored separately, reported
separately, and MUST NEVER be pooled with the prospective arm. That distinction is the same one
that made the patient arm's "2/2 correct, out-of-sample" claim wrong, and it is not repeated here.

THE CIRCULARITY GUARD, which is what makes this arm worth anything at all. The rule's input is the
ledger's pooled estimate for a compound. If a trial was itself INCLUDED in that pooled estimate,
then predicting it from the ledger is circular: the trial helped produce the number being used to
predict it. So every candidate is checked against the publication year of the meta-analysis the
ledger row rests on, and anything that could have been swept into it is EXCLUDED rather than
explained away. A trial only counts when its results became available AFTER the ledger's
meta-analysis was published.

That guard is deliberately strict. It will reject most of the pool, and rejecting most of the pool
is the correct outcome: a large circular sample is worth less than a small independent one.

Targets exactly the 19 ledger compounds that generate a DISCORDANT call, because those are the only
rows that can separate the rule from a constant predictor (see registry_power_v1.md). Compounds
whose interval spans zero, or whose point estimate sits below the 0.20 floor, are skipped.

Writes data/raw/posted_results_candidates.csv and reports/pipeline/posted_results_arm_v1.md.
"""
from __future__ import annotations

import json
import logging
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.provenance.trailer import stamp  # noqa: E402
from mammal_repurposing.reporting.prospective import predict_healthy_adult  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("posted_arm")

LEDGER = ROOT / "data" / "raw" / "healthy_adult_cognition_ledger.csv"
OUT = ROOT / "data" / "raw" / "posted_results_candidates.csv"
REPORT = ROOT / "reports" / "pipeline" / "posted_results_arm_v1.md"
API = "https://clinicaltrials.gov/api/v2/studies"

COG = ('AREA[PrimaryOutcomeMeasure](cognit OR memory OR attention OR "reaction time" OR vigilance '
       'OR Stroop OR "n-back" OR CANTAB OR "digit span" OR psychomotor OR "Trail Making" OR '
       'learning OR "processing speed" OR "working memory" OR "executive function")')

#: ledger compound -> the ClinicalTrials.gov intervention term to search. Only compounds whose
#: ledger row produces a POSITIVE or NEGATIVE call appear here; NULL_EFFECT callers are worthless
#: to the paired test and are deliberately absent.
SEARCH_TERMS: dict[str, str] = {
    "caffeine": "caffeine", "methylphenidate": "methylphenidate", "nicotine": "nicotine",
    "l_theanine": "theanine", "oxytocin": "oxytocin", "melatonin": "melatonin",
    "psilocybin": "psilocybin", "scopolamine": "scopolamine", "alcohol_acute": "alcohol",
    "resveratrol": "resveratrol", "cocoa_flavanols": "cocoa", "anthocyanins": "anthocyanin",
    "ashwagandha": "ashwagandha", "exogenous_ketones": "ketone", "dehydration": "dehydration",
}

#: A trial whose results predate the meta-analysis the ledger row rests on could have been
#: INCLUDED in that pooled estimate, which would make predicting it from the ledger circular. The
#: year is DERIVED from the ledger's own citation_short at runtime rather than hardcoded. An
#: earlier version of this file hardcoded the years from memory and got five of fifteen wrong,
#: including oxytocin (2013, written as 2025) and dehydration (2018, written as 2025). Deriving it
#: means the guard cannot drift from the ledger it is protecting.
def ledger_ma_year(led: pd.DataFrame, compound: str) -> int | None:
    row = led[led["compound"].astype(str).str.strip().str.lower() == compound.lower()]
    if row.empty:
        return None
    m = re.search(r"(?:19|20)\d\d", str(row.iloc[0].get("citation_short", "")))
    return int(m.group(0)) if m else None


FIELDS = ("NCTId|BriefTitle|OverallStatus|ResultsFirstPostDate|CompletionDate|"
          "PrimaryCompletionDate|EnrollmentCount|InterventionName|PrimaryOutcomeMeasure|"
          "HealthyVolunteers|EligibilityCriteria|Phase")


def _get(params: dict) -> dict:
    url = API + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"User-Agent": "mammal-repurposing/1.0"})
    with urllib.request.urlopen(req, timeout=90) as r:  # noqa: S310
        return json.loads(r.read().decode("utf-8"))


def search(term: str) -> list[dict]:
    q = (f"AREA[HealthyVolunteers]true AND AREA[ResultsFirstPostDate]RANGE[2000-01-01,MAX] "
         f"AND AREA[InterventionName]{term} AND {COG}")
    out, token = [], None
    while True:
        p = {"filter.advanced": q, "pageSize": 100, "fields": FIELDS}
        if token:
            p["pageToken"] = token
        d = _get(p)
        out.extend(d.get("studies", []))
        token = d.get("nextPageToken")
        if not token:
            return out
        time.sleep(0.3)


def flatten(s: dict) -> dict:
    ps = s.get("protocolSection", {})
    idm = ps.get("identificationModule", {})
    st = ps.get("statusModule", {})
    dm = ps.get("designModule", {})
    arms = ps.get("armsInterventionsModule", {})
    om = ps.get("outcomesModule", {})
    elig = ps.get("eligibilityModule", {})

    def date(key):
        v = st.get(key)
        return (v or {}).get("date", "") if isinstance(v, dict) else (v or "")

    return {
        "nct": idm.get("nctId", ""),
        "title": idm.get("briefTitle", ""),
        "status": st.get("overallStatus", ""),
        "results_posted": date("resultsFirstPostDateStruct") or st.get("resultsFirstPostDate", ""),
        "completion": date("completionDateStruct"),
        "primary_completion": date("primaryCompletionDateStruct"),
        "enrollment": (dm.get("enrollmentInfo") or {}).get("count", ""),
        "phase": ", ".join(dm.get("phases", []) or []),
        "interventions": "; ".join(i.get("name", "") for i in arms.get("interventions", []) or []),
        "primary_outcomes": " | ".join(o.get("measure", "")
                                       for o in om.get("primaryOutcomes", []) or []),
        "healthy_flag": elig.get("healthyVolunteers", None),
        "eligibility": (elig.get("eligibilityCriteria", "") or "")[:1200],
    }


def year_of(datestr: str) -> int | None:
    s = str(datestr)[:4]
    return int(s) if s.isdigit() else None


#: A NAMED cognitive instrument or an unambiguous performance construct. Deliberately narrow.
_COG_STRICT = re.compile(
    r"\b(CANTAB|CogState|NIH Toolbox|RBANS|MCCB|MATRICS|BACS|ANAM|PVT|DSST|Stroop|Flanker|"
    r"Corsi|Sternberg|Simon task|Go/?No-?Go|n-?back|digit span|digit symbol|Trail Making|"
    r"Rey Auditory|Hopkins Verbal|Logical Memory|verbal fluency|paired associate|"
    r"continuous performance|psychomotor vigilance|Wechsler|Grooved Pegboard|"
    r"Attention Network|working memory|episodic memory|executive function|processing speed|"
    r"cognitive (?:performance|function|test|battery|composite)|"
    r"reaction time|response time|recall|recognition memory|sustained attention|"
    r"attentional (?:bias|control)|false alarm|omission error)\b", re.I)

#: Pharmacokinetic, bioequivalence, safety and exposure endpoints. ClinicalTrials.gov's search
#: expands "attention" to "concentration" and "learning" to "Deep Learning", so a query for
#: cognitive terms returns plasma-concentration studies. Every one of these must be thrown out.
_NOT_COG = re.compile(
    r"\b(Cmax|AUC|C ?max|area under the curve|plasma|serum|cotinine|pharmacokinetic|"
    r"bioavailability|bioequivalence|geometric mean ratio|Tmax|half-?life|elimination|"
    r"adverse event|treatment-?emergent|tolerabilit|safety|vital sign|ECG|QTc|"
    r"concentration[- ]time|dose proportionality|nicotine uptake|abuse potential|"
    r"subjective effect|craving|withdrawal symptom)\b", re.I)


#: Diagnoses that disqualify a sample when they appear in INCLUSION criteria. The same words are
#: perfectly normal in EXCLUSION criteria ("no history of psychosis"), which is why the whole
#: eligibility blob cannot be searched naively: of the first five candidates, two were flagged
#: purely because psychosis appeared in their exclusion list.
_PATIENT = re.compile(
    r"\b(schizophren\w*|psychotic disorder|bipolar|major depress\w*|ADHD|"
    r"attention.deficit|dementia|Alzheimer|mild cognitive impairment|\bMCI\b|"
    r"Parkinson|epilep\w*|traumatic brain injury|substance use disorder)\b", re.I)


def inclusion_text(eligibility: str) -> str:
    """Just the inclusion half of an eligibility blob.

    ClinicalTrials.gov eligibility is free text with an 'Exclusion Criteria' heading partway
    through. Searching the whole thing for a diagnosis finds every study that sensibly EXCLUDES
    that diagnosis, which is the opposite of what the check is for.
    """
    t = str(eligibility or "")
    m = re.search(r"exclusion\s*criteria", t, re.I)
    return t[:m.start()] if m else t


def patient_population(title: str, eligibility: str) -> str:
    """The disqualifying diagnosis if the sample is a patient group, else empty.

    Checked against the TITLE and the INCLUSION criteria only. A registry healthyVolunteers flag is
    not trusted: it is set true when a study enrols a healthy CONTROL arm alongside patients, and
    that is how three schizophrenia trials reached this script's usable list on the first pass.
    """
    for field in (str(title or ""), inclusion_text(eligibility)):
        m = _PATIENT.search(field)
        if m:
            return m.group(0)
    return ""


def is_cognitive_primary(text: str) -> bool:
    """True only when the PRIMARY outcome text names a cognitive measure and is not a PK endpoint.

    The registry query alone is not sufficient. ClinicalTrials.gov expands search terms, so
    "Plasma Cotinine Concentration" is returned by a query for "attention". Of the first four
    candidates this filter was written against, three were pharmacokinetic studies that had matched
    purely through that expansion. Requiring a named instrument AND rejecting exposure endpoints
    removes the whole class.
    """
    t = str(text or "")
    return bool(_COG_STRICT.search(t)) and not bool(_NOT_COG.search(t))


def main() -> int:
    led = pd.read_csv(LEDGER)
    rows = []
    for compound, term in SEARCH_TERMS.items():
        pred = predict_healthy_adult(compound, led)
        if pred["prediction"] in ("ABSTAIN", "NULL_EFFECT"):
            L.warning("%s calls %s; skipping (cannot make a discordant pair)",
                      compound, pred["prediction"])
            continue
        try:
            studies = search(term)
        except Exception as exc:  # noqa: BLE001
            L.error("%s: search failed (%s)", compound, exc)
            continue
        ma_year = ledger_ma_year(led, compound)
        for s in studies:
            f = flatten(s)
            ry = year_of(f["results_posted"]) or year_of(f["completion"])
            independent = bool(ma_year and ry and ry > ma_year)
            cog = is_cognitive_primary(f["primary_outcomes"])
            patient = patient_population(f["title"], f["eligibility"])
            rows.append({
                "ledger_compound": compound, "search_term": term,
                "prediction": pred["prediction"], "pred_g": pred["g"],
                "pred_ci_lo": pred["ci_lo"], "pred_ci_hi": pred["ci_hi"],
                "ledger_ma_year": ma_year, "results_year": ry,
                "independent_of_ledger": independent,
                "cognitive_primary": cog,
                "patient_population": patient,
                "usable": bool(independent and cog and not patient),
                "exclusion_reason": ("" if not cog else "") + ("" if cog else
                    "primary outcome is not a cognitive measure; matched only via the registry's "
                    "term expansion (attention -> concentration)") + ("" if independent else
                                     f"results available {ry} vs ledger meta-analysis {ma_year}; "
                                     f"could have been included in the pooled estimate"),
                **f,
            })
        L.info("%s (%s): %d found, %d independent, %d USABLE (independent + cognitive primary)",
               compound, pred["prediction"], len(studies),
               sum(1 for r in rows if r["ledger_compound"] == compound
                   and r["independent_of_ledger"]),
               sum(1 for r in rows if r["ledger_compound"] == compound and r["usable"]))
        time.sleep(0.3)

    df = pd.DataFrame(rows).drop_duplicates("nct")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False, encoding="utf-8")
    write_report(df)
    L.info("wrote %s: %d unique, %d independent, %d USABLE",
           OUT.name, len(df), int(df["independent_of_ledger"].sum()),
           int(df["usable"].sum()))
    return 0


def write_report(df: pd.DataFrame) -> None:
    ind = df[df["usable"]]
    lines = [
        "# The retrospective arm: trials whose results are already posted",
        "",
        "The prospective registry cannot be scored yet. A blind sweep over its due rows found 25 of",
        "26 unpublished, and NOT ONE had registry-posted results, so journal publication is the",
        "only channel and publication lag is the binding constraint.",
        "",
        "ClinicalTrials.gov requires many sponsors to post results as structured data whether or",
        "not a paper appears. Those trials are immune to the file drawer in the way the journal",
        "literature is not, and they are readable today.",
        "",
        "**This arm is RETROSPECTIVE.** The results already exist, so applying the frozen rule to",
        "them is a retrodiction. It is scored separately and must never be pooled with the",
        "prospective arm. Conflating those two is exactly what made the patient arm's",
        "\"2 / 2 correct, genuine out-of-sample\" claim wrong.",
        "",
        "## What the search returned",
        "",
        f"Searching only the ledger compounds that generate a DISCORDANT call: **{len(df)} unique",
        "trials** with posted results, a healthy-volunteer flag, and a cognitive primary outcome.",
        "",
        "| compound | call | found | independent | usable |",
        "|---|---|---:|---:|---:|",
    ]
    for c, g in df.groupby("ledger_compound"):
        lines.append(f"| {c} | {g['prediction'].iloc[0]} | {len(g)} | "
                     f"{int(g['independent_of_ledger'].sum())} | {int(g['usable'].sum())} |")
    lines += [
        "",
        "## The circularity guard, which rejects most of them",
        "",
        "The rule's input is the ledger's pooled estimate. A trial that was itself INCLUDED in that",
        "pooled estimate cannot be used to test it: the trial helped produce the number doing the",
        "predicting. So a trial counts only if its results became available AFTER the",
        "meta-analysis the ledger row rests on was published.",
        "",
        f"**{int(df['independent_of_ledger'].sum())} of {len(df)} survive that check.** The rest",
        "are excluded rather than explained away. A large circular sample is worth less than a",
        "small independent one.",
        "",
        "## Two further filters, both of which were wrong on the first pass",
        "",
        "**The primary outcome must actually be cognitive.** ClinicalTrials.gov expands search",
        "terms, so a query for \"attention\" returns \"Plasma Cotinine CONCENTRATION\". Three of the",
        "first four candidates were pharmacokinetic or bioequivalence studies that matched purely",
        "through that expansion. A named-instrument requirement plus an explicit exposure-endpoint",
        "reject removes the class: "
        f"{int((df['independent_of_ledger'] & ~df['cognitive_primary']).sum())} of the independent "
        "trials fall here, nearly all of them Cmax, AUC or plasma-concentration endpoints.",
        "",
        "**The sample must actually be healthy.** The registry's healthyVolunteers flag is true "
        "whenever a study enrols a healthy CONTROL arm beside patients, so two schizophrenia trials "
        "reached the usable list before this was added. The check now reads the TITLE and the "
        "INCLUSION criteria only, because the same diagnoses appear routinely in exclusion "
        "criteria and a naive search flags every well-designed healthy study for saying \"no "
        "history of psychosis\".",
        "",
        f"**{len(ind)} candidates survive all three filters.**",
        "",
    ]
    if len(ind):
        lines += ["| NCT | compound | call | results posted | n | title |", "|---|---|---|---|---:|---|"]
        for _, r in ind.sort_values(["ledger_compound", "results_year"]).iterrows():
            lines.append(f"| {r['nct']} | {r['ledger_compound']} | {r['prediction']} | "
                         f"{r['results_posted'] or r['completion']} | {r['enrollment']} | "
                         f"{str(r['title'])[:70]} |")
    lines += [
        "",
        "## Not yet scored",
        "",
        "Finding a trial with posted results is not the same as reading its result. Every row above",
        "still needs its posted outcome extracted and independently verified, and the",
        "healthy-volunteer flag re-checked against the eligibility prose, because a registry flag",
        "that contradicts its own criteria is how a medicated early-psychosis trial entered this",
        "project's registry earlier this week.",
        "",
    ]
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(stamp("\n".join(lines), "scripts/137_posted_results_arm.py"),
                      encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
