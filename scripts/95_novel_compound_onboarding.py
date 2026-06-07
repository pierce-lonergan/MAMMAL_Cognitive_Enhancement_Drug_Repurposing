"""F2 - novel-compound onboarding engine (end-to-end run + validation).

Turns the class-prognostic retrospective result into a PROSPECTIVE screen: given an
arbitrary SMILES, route it to a known cognition mechanism class and return that
class's calibrated clinical-g prior, or ABSTAIN. Two stages:

  1. VALIDATION - leave-one-compound-out class recovery on the exemplar base. The
     honest internal metric: does structure recover the TRUE mechanism class? Held-
     out drugs are re-routed using only the remaining exemplars + priors.
  2. DEMO - score a small set of novel compounds (real CNS drugs not in the ledger,
     plus peripheral out-of-manifold negatives) to show routing + abstention.

Outputs reports/pipeline/novel_compound_onboarding_v1.md + a scored CSV. CPU only
(RDKit + numpy/pandas), so it runs in CI. Nothing here fabricates a clinical
outcome: the returned g is a PREDICTION from the class prior, always labelled.
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from mammal_repurposing.reporting.trial_watch import load_combined_ledger
from mammal_repurposing.validation.novel_compound import (
    MIN_CLASS_N, TAU_HIGH, TAU_MARGIN, TAU_OOD,
    build_class_priors, build_exemplars, loco_class_recovery, score_catalogue,
)

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("f2")

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
LEDGERS = [
    RAW / "clinical_outcomes_ledger.csv",
    RAW / "clinical_outcomes_ledger_EXTENSION.csv",
    RAW / "clinical_outcomes_ledger_CTGOV.csv",
    RAW / "clinical_outcomes_ledger_RESEARCH.csv",
]
SMILES = RAW / "ledger_compound_smiles.csv"
DEMO = RAW / "novel_demo_compounds.csv"
REPORT = ROOT / "reports" / "pipeline" / "novel_compound_onboarding_v1.md"
SCORED = ROOT / "reports" / "pipeline" / "novel_compound_scored.csv"


def _fmt_g(row) -> str:
    if not np.isfinite(row["prior_g"]):
        return "-"
    return f"{row['prior_g']:+.2f} [{row['g_ci_lo']:+.2f}, {row['g_ci_hi']:+.2f}]"


def main() -> int:
    led = load_combined_ledger(LEDGERS)
    smi = pd.read_csv(SMILES)[["compound", "smiles"]]
    ex = build_exemplars(led, smi)
    priors = build_class_priors(led)
    n_ex = sum(len(v) for v in ex.values())
    n_multi = sum(1 for v in ex.values() if len(v) >= 2)

    # ---- 1. validation: leave-one-compound-out class recovery ----
    rec = loco_class_recovery(led, smi)
    det = rec["detail"]
    misroutes = det[(~det["correct"]) & (~det["abstained"])]

    # ---- 2. demo scoring ----
    scored = pd.DataFrame()
    if DEMO.exists():
        demo = pd.read_csv(DEMO).rename(columns={"name": "id"})
        scored = score_catalogue(demo[["id", "smiles"]], ex, priors)
        scored = scored.merge(
            demo[["id", "rationale"]].rename(columns={"id": "query_id"}),
            on="query_id", how="left")
        SCORED.parent.mkdir(parents=True, exist_ok=True)
        scored.to_csv(SCORED, index=False)

    # ---- report ----
    Ls: list[str] = []
    Ls.append("# F2 - Novel-compound onboarding engine")
    Ls.append("")
    Ls.append("**Question.** Can the class-prognostic result (F3: class-LOCO AUROC "
              "~0.92; F1: class is the resolution limit) be turned into a *prospective* "
              "screen? Given an arbitrary novel SMILES, route it to a known cognition "
              "mechanism class and return that class's calibrated clinical-*g* prior - "
              "or ABSTAIN. Reproduced by `scripts/95_novel_compound_onboarding.py`.")
    Ls.append("")
    Ls.append("The engine re-ranks KNOWN mechanisms for cognition; it does not invent "
              "new ones. The leave-one-class-out=0.00 result (no signal for genuinely "
              "novel mechanisms) is enforced as a hard ABSTAIN guardrail.")
    Ls.append("")
    Ls.append("## Exemplar base")
    Ls.append("")
    Ls.append(f"- Ledger compounds with a parseable SMILES: **{n_ex}** "
              f"(`data/raw/ledger_compound_smiles.csv`).")
    Ls.append(f"- Mechanism classes with >= 1 exemplar: **{len(ex)}** of "
              f"{led['mechanism_class'].nunique()}; with >= 2 exemplars: **{n_multi}**.")
    Ls.append("")
    Ls.append("## 1. Validation - leave-one-compound-out class recovery")
    Ls.append("")
    Ls.append("Each SMILES-backed ledger compound is held out, the exemplar library + "
              "class priors are rebuilt from the rest, and the held-out compound is "
              "re-routed. Only compounds whose class keeps a sibling after holdout are "
              "evaluable. The test asks: does structure recover the TRUE class?")
    Ls.append("")
    Ls.append("| metric | value |")
    Ls.append("|---|---|")
    Ls.append(f"| evaluable held-out compounds | {rec['n_evaluable']} |")
    Ls.append(f"| routed (not abstained) | {rec['n_routed']} |")
    Ls.append(f"| **top-1 class recovery (routed)** | **{rec['top1_acc']:.3f}** |")
    Ls.append(f"| abstention rate | {rec['abstain_rate']:.3f} |")
    Ls.append(f"| accuracy when sim >= {TAU_HIGH:.2f} | {rec['acc_when_sim_ge_high']:.3f} |")
    Ls.append("")
    Ls.append(f"Decision thresholds (calibrated here, locked in `novel_compound.py`): "
              f"out-of-manifold floor TAU_OOD={TAU_OOD}, HIGH-tier TAU_HIGH={TAU_HIGH}, "
              f"ambiguity margin TAU_MARGIN={TAU_MARGIN}, min class members "
              f"MIN_CLASS_N={MIN_CLASS_N}.")
    Ls.append("")
    if len(misroutes):
        Ls.append(f"**Mis-routes ({len(misroutes)}):**")
        Ls.append("")
        Ls.append("| compound | true class | assigned | similarity |")
        Ls.append("|---|---|---|---|")
        for _, r in misroutes.iterrows():
            Ls.append(f"| {r['compound']} | {r['true_class']} | {r['assigned']} "
                      f"| {r['similarity']:.2f} |")
        Ls.append("")
        Ls.append("Note: the residual mis-route is an *enantiomer* blind spot - the "
                  "2D ECFP4 fingerprint cannot separate stereoisomers whose mechanisms "
                  "differ (e.g. (-)-phenserine, an AChE inhibitor, vs its (+)-enantiomer "
                  "posiphen/buntanetap, an APP-translation inhibitor; identical 2D "
                  "structure -> Tanimoto 1.0).")
    else:
        Ls.append("**No mis-routes.**")
    Ls.append("")
    Ls.append("High abstention is the guardrail working: where a held-out drug has no "
              "close structural analog among its class siblings, the engine refuses "
              "rather than guess. Coverage (exemplar SMILES per class) is the lever to "
              "lower abstention; `scripts/_expand_ledger_smiles.py` grew it 31 -> "
              f"{n_ex}.")
    Ls.append("")

    if len(scored):
        n_ab = int((scored["tier"] == "ABSTAIN").sum())
        Ls.append("## 2. Demo - novel compounds (not in the ledger)")
        Ls.append("")
        Ls.append(f"{len(scored)} compounds: real CNS drugs + peripheral out-of-manifold "
                  f"negatives. {len(scored) - n_ab} routed, {n_ab} abstained. The "
                  "predicted *g* is the assigned class's prior (a model output), not a "
                  "measured outcome.")
        Ls.append("")
        Ls.append("| compound | tier | assigned class | sim | predicted g [90% CrI] | "
                  "P(success) | basis |")
        Ls.append("|---|---|---|---|---|---|---|")
        for _, r in scored.iterrows():
            ac = r["assigned_class"] or "-"
            sim = f"{r['similarity']:.2f}" if np.isfinite(r["similarity"]) else "-"
            ps = f"{r['p_success']:.2f}" if np.isfinite(r["p_success"]) else "-"
            Ls.append(f"| {r['query_id']} | {r['tier']} | {ac} | {sim} | {_fmt_g(r)} "
                      f"| {ps} | {r['reason']} |")
        Ls.append("")

    Ls.append("## Guardrails (non-negotiable)")
    Ls.append("")
    Ls.append("1. **Out-of-manifold -> ABSTAIN.** Max Tanimoto to any known class "
              f"< {TAU_OOD} means the compound is not near any precedented cognition "
              "chemotype; the engine abstains (it cannot invent a mechanism).")
    Ls.append("2. **Allosteric downgrade (V6.A).** Structural/DTI-profile class "
              "assignment is unreliable for allosteric chemotypes, so allosteric-flagged "
              "classes are capped at MED with a note.")
    Ls.append("3. **Thin prior -> LOW; ambiguous opposite-sign tie -> ABSTAIN.**")
    Ls.append("")
    Ls.append("## Limitations")
    Ls.append("")
    Ls.append("- 2D-structural routing (ECFP4 + Murcko); the pluggable multi-head "
              "DTI-profile signal (MAMMAL/MMAtt-DTA/PSICHIC/BALM nearest-class) is a "
              "GPU upgrade wired via `external_class_scores` but not run here.")
    Ls.append("- Enantiomer mechanism-switches are a known blind spot (see mis-routes).")
    Ls.append("- The prior is only as good as the ledger (n=125); singleton classes "
              "carry no usable prior and route LOW.")
    Ls.append("")

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(Ls), encoding="utf-8")
    L.info("Wrote %s", REPORT)
    if len(scored):
        L.info("Wrote %s", SCORED)
    L.info("F2: LOCO top1 %.3f on %d routed (%.0f%% abstain); exemplars %d in %d classes",
           rec["top1_acc"], rec["n_routed"], 100 * rec["abstain_rate"], n_ex, len(ex))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
