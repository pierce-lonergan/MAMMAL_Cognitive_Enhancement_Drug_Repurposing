"""The healthy-adult arm of the forward registry, which the project did not have.

Every row of `data/raw/prospective_predictions.csv` is a PATIENT population: CIAS, Fragile X,
Alzheimer's, schizophrenia. Nothing in the project's forward test has ever concerned healthy adults,
which is the population the whole question is about. This builds that arm.

BE PRECISE ABOUT WHAT IT DOES AND DOES NOT TEST. It tests whether the healthy-adult ledger
PREDICTS new healthy-adult readouts, which is a real test of the ledger's forward calibration and
one the project has never run. It does NOT test G1. G1 is the claim that no DURABLE post-washout
gain has been verified, and the ledger's representative_g values are acute, on-drug pooled
estimates, so this arm deliberately ABSTAINS on post-acute endpoints rather than predicting them
from the wrong quantity. Of 64 confirmed readouts only a handful study persistence at all, and the
one with a genuine post-acute cognitive endpoint is excluded for exactly this reason. Testing G1
forward would need post-washout designs that mostly do not exist.

INPUT is `data/raw/pending_readouts_2026-09.csv`: 83 candidate registrations found by an
eight-lane blind sweep of ClinicalTrials.gov, PROSPERO, ISRCTN, UMIN and OSF, each one then handed
to an independent skeptic who had to re-find it and could reject it as not-found, wrong-population,
not-cognitive, already-published or misdescribed. 69 survived. Zero came back
REJECT_ID_NOT_FOUND, which is the fabrication signature, across 83 checks.

THE RULE IS NOT CHOSEN HERE. `predict_healthy_adult` was written and committed at 2baedf8, before
this file existed and before the sweep's compound list had been read. This script only applies it.

THE MAPPING IS THE ONLY JUDGEMENT IN THIS FILE, so it is written out in full below rather than
inferred by string matching. A registration maps to a ledger compound only when its intervention is
unambiguously that ONE compound. Everything else abstains, including:

  - multi-compound reviews. CRD420261324802 pools methylphenidate, dextroamphetamine,
    lisdexamfetamine, atomoxetine and modafinil, so its pooled estimate is not a statement about
    any single ledger row.
  - blends. NCT06780774 is 120 mg caffeine plus vitamins, minerals and botanicals, so a positive
    result there is not a caffeine result.
  - compounds with no ledger row at all: alpha-GPC, magnesium, klotho, rapamycin, fisetin,
    nicotinamide riboside, tocotrienol, haskap, carotenoids, antipsychotics.

Abstaining is expected to be the common outcome and is the honest one. A registry that forces a
call on every row manufactures a track record out of guesses.

Writes data/raw/prospective_healthy_adult.csv and reports/pipeline/prospective_healthy_adult_v1.md.
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.provenance.trailer import stamp  # noqa: E402
from mammal_repurposing.reporting.prospective import (  # noqa: E402
    PRED_ABSTAIN,
    healthy_adult_baseline,
    predict_healthy_adult,
)

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("ha_forward")

PENDING = ROOT / "data" / "raw" / "pending_readouts_2026-09.csv"
LEDGER = ROOT / "data" / "raw" / "healthy_adult_cognition_ledger.csv"
OUT = ROOT / "data" / "raw" / "prospective_healthy_adult.csv"
REPORT = ROOT / "reports" / "pipeline" / "prospective_healthy_adult_v1.md"

FROZEN_ON = "2026-09-20"

#: registration identifier -> the ONE ledger compound it tests. Anything absent abstains.
#: Read the docstring before adding to this: an entry here is a claim that the registration's
#: readout will be an estimate FOR that compound, not merely that the compound is involved.
COMPOUND_MAP: dict[str, str] = {
    # --- PROSPERO / OSF registered reviews. These matter most: the ledger's inclusion rule is a
    # --- clean healthy-adult meta-analysis, so each of these IS a future ledger row.
    "CRD420251140305": "creatine",
    "CRD420251160876": "l_theanine",
    "CRD420251163859": "caffeine",
    "CRD420251488583": "caffeine",
    "CRD420250650923": "dextroamphetamine",
    "CRD420261386551": "cocoa_flavanols",
    "CRD420251266075": "anthocyanins",
    "CRD420261499408": "caffeine_plus_taurine",
    "CRD420261280163": "oxytocin",
    "CRD420261431354": "menopausal_hormone_therapy",
    "kw54p": "tyrosine",
    # --- single-compound trials
    "NCT05325502": "caffeine",
    "NCT05588934": "caffeine",
    "NCT06763172": "caffeine",
    "NCT07469852": "caffeine",
    "NCT06041048": "modafinil",
    "NCT07231497": "modafinil",
    "NCT06906848": "citicoline",
    "NCT07408180": "nicotine",
    "NCT05301608": "psilocybin",
    "NCT06367738": "psilocybin",
    "NCT06692192": "psilocybin",
    "NCT06768944": "psilocybin",
    "NCT07079852": "psilocybin",
    "NCT07449351": "psilocybin_lsd_microdosing",  # MICRODOSING: g -0.34, not the -1.13 full-dose row
    "NCT07732166": "psilocybin",
    "NCT07818564": "psilocybin",
    "ISRCTN10543404": "scopolamine",
    "ISRCTN16663180": "creatine",
    "ISRCTN29176549": "cocoa_flavanols",
    "NCT05519644": "exogenous_ketones",
    "NCT07051655": "exogenous_ketones",
    "NCT07109245": "insulin_intranasal",
    "NCT06530459": "menopausal_hormone_therapy",
    "NCT07666685": "dietary_nitrate",
    "NCT07814300": "omega_3",
    "NCT07242430": "multivitamin_mineral",
}

#: Mapped registrations whose PRIMARY outcome is not an objective cognitive measure. The sweep's
#: own verifiers flagged every one of these in prose and the first version of this script ignored
#: them. A prediction about a cognitive effect cannot be graded on an imaging, self-report or BOLD
#: primary; that is precisely the category error found in the patient arm, where a cognition claim
#: was being scored against a PANSS negative-symptom primary.
ENDPOINT_EXCLUDE: dict[str, str] = {
    "NCT06367738": "primary outcome is imaging (mean diffusivity), not cognition",
    "NCT06768944": "primary outcome is self-report (Persisting Effects Questionnaire)",
    "NCT06041048": "sole primary outcome is BOLD brain activity, not a cognitive test",
    "NCT07732166": "verifier flagged the primary as borderline on domain; excluded rather than "
                   "resolved by preference",
}

#: Mapped registrations measuring a POST-ACUTE or persisting effect. The ledger's
#: representative_g values are ACUTE, on-drug pooled estimates, so they cannot predict what
#: remains after washout. Predicting a post-washout readout from an on-drug estimate would be
#: exactly the durability confusion (G1) this project exists to keep straight.
TIMING_EXCLUDE: dict[str, str] = {
    "NCT06692192": "follow-up study of persisting effects; the ledger estimate is acute/on-drug",
}

#: Mapped registrations where cognition is a SECONDARY endpoint. These keep their call, because a
#: secondary objective cognitive measure is still a real readout, but they are tagged so any
#: future scoring can stratify rather than pooling them with primary-endpoint predictions.
SECONDARY_ENDPOINT: set[str] = {
    "NCT06530459", "ISRCTN16663180", "ISRCTN29176549", "CRD420261280163",
}

#: Why each non-mapped confirmed registration abstains. Stated so the abstentions are auditable
#: rather than looking like an oversight.
ABSTAIN_REASON: dict[str, str] = {
    "CRD420261307872": "review pools methylphenidate AND modafinil; not a single-compound estimate",
    "CRD420261324802": "review pools five stimulants; not a single-compound estimate",
    "CRD42024549774": "'any natural nootropic', individual and combined; scope is not a compound",
    "CRD420250627818": "carotenoids; no ledger row",
    "CRD420261356864": "antipsychotics; no ledger row and not an enhancement question",
    "NCT06780774": "caffeine blended with vitamins, minerals and botanicals; caffeine not isolable",
    "NCT07530185": "multi-arm tyrosine/caffeine/combination; registration does not state doses",
    "NCT07135232": "sarmentosin plus L-theanine, factorially crossed with resistance training",
    "NCT07267845": "alpha-GPC; no ledger row",
    "NCT07235878": "magnesium; no ledger row",
    "NCT07216781": "klotho gene therapy; no ledger row",
    "NCT07285629": "klotho plus follistatin gene therapy; no ledger row",
    "NCT06682767": "rapamycin; no ledger row",
    "NCT07475546": "multi-drug geroprotective stack; no ledger row and not isolable",
    "NCT06431932": "fisetin; no ledger row",
    "NCT05483465": "nicotinamide riboside; no ledger row",
    "NCT07637487": "tocotrienol-rich fraction; no ledger row",
    "NCT07119788": "haskap/honeyberry powder; no ledger row",
    "NCT06690892": "beef portions; a whole-food exposure, not a compound",
    "NCT07388576": "AP-Brain proprietary supplement; composition not isolable",
    "NCT07451496": "personalised supplement stack plus supervised exercise; not isolable",
    "NCT07319117": "multi-ingredient 'Think Tank' formulation; citicoline not isolable",
    "ISRCTN85608856": "omega-3 plus curcumin plus citicoline; no single compound isolable",
    "ISRCTN66268301": "all arms on an exercise and diet intervention; compound not isolable",
    "ISRCTN75484092": "four-arm design; registration does not isolate a single compound",
    "UMIN000057063": "nucleotide-containing capsule; no ledger row",
}


def build() -> pd.DataFrame:
    pend = pd.read_csv(PENDING)
    conf = pend[pend["verdict"] == "CONFIRMED"].drop_duplicates("identifier").copy()
    led = pd.read_csv(LEDGER)
    L.info("%d confirmed registrations (%d rows before de-duplicating cross-lane hits)",
           len(conf), int((pend["verdict"] == "CONFIRMED").sum()))

    rows = []
    for _, r in conf.iterrows():
        ident = str(r["identifier"]).strip()
        compound = COMPOUND_MAP.get(ident)
        if compound is not None and ident in ENDPOINT_EXCLUDE:
            pred = {"prediction": PRED_ABSTAIN, "g": None, "ci_lo": None, "ci_hi": None,
                    "basis": "endpoint_not_cognitive",
                    "reason": ENDPOINT_EXCLUDE[ident]}
            compound = None
        elif compound is not None and ident in TIMING_EXCLUDE:
            pred = {"prediction": PRED_ABSTAIN, "g": None, "ci_lo": None, "ci_hi": None,
                    "basis": "post_acute_endpoint",
                    "reason": TIMING_EXCLUDE[ident]}
            compound = None
        elif compound is None:
            pred = {"prediction": PRED_ABSTAIN, "g": None, "ci_lo": None, "ci_hi": None,
                    "basis": "not_mapped",
                    "reason": ABSTAIN_REASON.get(ident, "no mapping recorded for this "
                                                        "registration")}
        else:
            pred = predict_healthy_adult(compound, led)
        rows.append({
            "identifier": ident, "registry": r["registry"], "lane": r["lane"],
            "ledger_compound": compound or "", "intervention": r["intervention"],
            "population": r["population"], "n": r["n"],
            "primary_cognitive_outcome": r["primary_cognitive_outcome"],
            "registry_status": r["status"], "expected_readout": r["expected_readout"],
            "prediction": pred["prediction"], "pred_g": pred["g"],
            "pred_ci_lo": pred["ci_lo"], "pred_ci_hi": pred["ci_hi"],
            "prediction_basis": pred["basis"], "prediction_reason": pred["reason"],
            "endpoint_tier": "secondary" if ident in SECONDARY_ENDPOINT else "primary",
            "prediction_date": FROZEN_ON, "rule_committed_at": "2baedf8",
            "status": "PENDING", "actual_outcome": "", "url": r["url"],
        })
    return pd.DataFrame(rows).sort_values(["prediction", "ledger_compound", "identifier"])


def write_report(df: pd.DataFrame, base: dict) -> None:
    called = df[df["prediction"] != PRED_ABSTAIN]
    counts = df["prediction"].value_counts()
    lines = [
        "# The healthy-adult forward registry",
        "",
        "The project's existing forward registry is entirely patient populations, so nothing in it",
        "concerns the population the question is actually about. This is that arm.",
        "",
        "**What it tests, precisely.** Whether the healthy-adult ledger predicts new healthy-adult",
        "readouts. That is a test of the ledger's forward calibration, and the project has never",
        "run one. It is NOT a test of G1. G1 concerns DURABLE post-washout gain, and the ledger's",
        "estimates are acute and on-drug, so this arm abstains on post-acute endpoints rather than",
        "predicting them from the wrong quantity. Testing G1 forward would need post-washout",
        "designs, and the sweep found that they mostly do not exist.",
        "",
        f"**{len(df)} confirmed pending readouts.** Found by an eight-lane blind sweep of",
        "ClinicalTrials.gov, PROSPERO, ISRCTN, UMIN and OSF; every registration was then handed to",
        "an independent skeptic who had to re-find it and could reject it as not-found,",
        "wrong-population, not-cognitive, already-published or misdescribed. 69 of 83 candidates",
        "survived, and **zero** came back as the fabrication signature (identifier does not",
        "resolve). Five identifiers were found independently by two lanes, which is the only",
        "cross-check the design gives for free.",
        "",
        "## The predictions",
        "",
        "The rule was committed at `2baedf8`, BEFORE this file existed and before the sweep's",
        "compound list was read. This script applies it and does not choose it. Check the commit",
        "order in git rather than taking that sentence for it.",
        "",
        "| call | n |",
        "|---|---:|",
    ]
    for k, v in counts.items():
        lines.append(f"| {k} | {v} |")
    lines += [
        "",
        f"**{len(called)} of {len(df)} registrations get a call; {counts.get(PRED_ABSTAIN, 0)} "
        f"abstain.** Abstention is the expected common outcome, not a failure. Every abstention "
        f"carries its reason in the CSV, so the set can be audited rather than assumed to be an "
        f"oversight.",
        "",
        "### Every call, with the evidence behind it",
        "",
        "| identifier | registry | compound | call | ledger g [CI] | readout |",
        "|---|---|---|---|---|---|",
    ]
    for _, r in called.iterrows():
        g = "" if pd.isna(r["pred_g"]) else f"{r['pred_g']:.2f}"
        ci = ("" if pd.isna(r["pred_ci_lo"])
              else f"[{r['pred_ci_lo']:.2f}, {r['pred_ci_hi']:.2f}]")
        lines.append(f"| `{r['identifier']}` | {r['registry']} | {r['ledger_compound']} | "
                     f"**{r['prediction']}** | {g} {ci} | {str(r['expected_readout'])[:42]} |")

    lines += [
        "",
        "## The opponent",
        "",
        "Any accuracy this registry eventually reports must be read against a constant predictor,",
        "because that is the mistake the patient arm made: it reported 100% accuracy that a fixed",
        "answer also achieves.",
        "",
        f"- Ledger rows carrying an interval: **{base['n_with_interval']}**",
        f"- Positive: {base['n_positive']}, null: {base['n_null']}, negative: {base['n_negative']}",
        f"- So always answering {base['constant_call']} scores "
        f"**{base['constant_accuracy']:.3f}**, and that is the bar.",
        "",
        "One caveat on that number, stated rather than buried: it is the fraction of LEDGER ROWS",
        "that are null, which is a proxy for the fraction of FUTURE READOUTS that will be null. The",
        "two are not the same quantity. It is the best available opponent, not an exact one.",
        "",
        "## What makes this arm worth waiting on",
        "",
        "Of the calls above, the PROSPERO and OSF registrations are the ones that matter most. The",
        "ledger's inclusion rule is specifically a clean healthy-adult meta-analysis whose interval",
        "lies entirely above zero, so a registered-but-unpublished systematic review on a ledger",
        "compound is not merely evidence about the prediction: it is a future ledger row. When it",
        "publishes it will either confirm or move the very estimate the prediction was made from.",
        "",
        "## Nothing here is scored yet",
        "",
        "Every row is PENDING and `actual_outcome` is empty. That is the point. The registry",
        "becomes evidence only as readouts land, and the rows must not be edited when they do:",
        "a frozen prediction rewritten to match what happened is not a prediction. Re-run this",
        "script to refresh registry status; fill `actual_outcome` and nothing else.",
        "",
    ]
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(stamp("\n".join(lines), "scripts/135_prospective_healthy_adult.py"),
                      encoding="utf-8")


def main() -> int:
    df = build()
    led = pd.read_csv(LEDGER)
    base = healthy_adult_baseline(led)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False, encoding="utf-8")
    write_report(df, base)
    L.info("wrote %s (%d rows)", OUT.name, len(df))
    L.info("calls: %s", df["prediction"].value_counts().to_dict())
    L.info("constant-%s baseline to beat: %.3f", base["constant_call"],
           base["constant_accuracy"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
