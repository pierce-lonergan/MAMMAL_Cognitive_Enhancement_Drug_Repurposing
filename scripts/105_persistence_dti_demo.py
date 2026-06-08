"""Stage 3 of the persistence-target DTI module (PERSEUS roadmap #2): a HELD-OUT
demonstration that the substrate-hypothesis head fires correctly end-to-end.

None of these compounds were calibration anchors, so this is a generalization test of the
calibrated channels:
  - ablative_heldout (obatoclax/gossypol/sabutoclax): BH3-mimetics NOT used to calibrate
    BCL2/BCL-xL -> should engage the ablative channel -> promotes_durable=True.
  - senolytic_flavonoid (fisetin/quercetin): senolytic but NOT direct BH3-mimetics -> the
    DTI-BCL2 channel should stay SILENT (precision); these are instead caught by the L3
    Tanimoto-to-flavonoid senolytic detector, so the two channels are complementary.
  - capability_heldout (chidamide/givinostat/tazemetostat): chromatin writers -> capability
    flag at most, never durable.
  - non_persistence (galantamine/rivastigmine/atomoxetine/citalopram/lamotrigine): off-
    substrate -> abstain.

GPU. Run with the MAMMAL venv:
  .venv-mammal/Scripts/python.exe scripts/105_persistence_dti_demo.py
"""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from mammal_repurposing.engine.persistence_dti import (
    load_calibration, load_panel, score_compound_against_panel, substrate_hypothesis,
)

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("persistence_dti_demo")

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
INTERIM = ROOT / "data" / "interim"
RESULTS = ROOT / "data" / "results"
PANEL = INTERIM / "persistence_targets.csv"
CALIB = RESULTS / "persistence_dti_calibration.json"
DEMO = RAW / "persistence_dti_demo.csv"
SCORED_OUT = RESULTS / "persistence_dti_demo_scored.csv"
REPORT = ROOT / "reports" / "pipeline" / "persistence_dti_demo.md"


def main() -> int:
    panel = load_panel(PANEL)
    calib = load_calibration(CALIB)
    demo = pd.read_csv(DEMO)

    from mammal_repurposing.scoring.model_loader import load_dti_model
    model, tok = load_dti_model()

    rows = []
    for _, r in demo.iterrows():
        scores = score_compound_against_panel(model, tok, r["smiles"], panel)
        h = substrate_hypothesis(scores, panel, calib)
        engaged = "; ".join(f"{e['gene']}={e['pkd']:.2f}" for e in h.engaged) or "-"
        rows.append({
            "compound": r["compound"], "expected_class": r["expected_class"],
            "substrate_hypothesis": h.substrate_hypothesis or "ABSTAIN",
            "promotes_durable": h.promotes_durable,
            "engaged": engaged, "capability_flags": ",".join(h.capability_flags) or "-",
            "abstained_targets": ",".join(h.abstained_targets) or "-",
        })
        L.info("%-14s [%s] -> %s durable=%s | %s", r["compound"], r["expected_class"],
               h.substrate_hypothesis or "ABSTAIN", h.promotes_durable, engaged)

    out = pd.DataFrame(rows)
    RESULTS.mkdir(parents=True, exist_ok=True)
    out.to_csv(SCORED_OUT, index=False)

    # honest scorecard against the expectation per class
    durable = out[out["promotes_durable"]]["compound"].tolist()
    ablative_hits = out[(out["expected_class"] == "ablative_heldout") & out["promotes_durable"]]
    flav_false = out[(out["expected_class"] == "senolytic_flavonoid") & out["promotes_durable"]]
    nonpersist_false = out[(out["expected_class"] == "non_persistence") & (out["substrate_hypothesis"] != "ABSTAIN")]

    write_report(out, durable, ablative_hits, flav_false, nonpersist_false)
    L.info("DEMO: durable=%s | held-out ablative recovered %d/%d | flavonoid false-durable %d | "
           "non-persistence leaks %d", durable, len(ablative_hits),
           int((out["expected_class"] == "ablative_heldout").sum()),
           len(flav_false), len(nonpersist_false))
    return 0


def write_report(out, durable, ablative_hits, flav_false, nonpersist_false) -> None:
    n_abl = int((out["expected_class"] == "ablative_heldout").sum())
    Ls = ["# Persistence-target DTI module - held-out demonstration", "",
          "End-to-end test of `engine/persistence_dti.py:substrate_hypothesis` on compounds "
          "NOT used to calibrate the panel. Reproduced by "
          "`scripts/105_persistence_dti_demo.py`.", "",
          "## Scorecard", "",
          f"- **Held-out ablative (BH3-mimetic) recovered: {len(ablative_hits)}/{n_abl}** "
          f"({', '.join(ablative_hits['compound']) or 'none'}) - generalization of the senolytic channel.",
          f"- **Flavonoid-senolytic false-durable: {len(flav_false)}** "
          "(expected 0 - flavonoid senolytics are not BH3-mimetics; the DTI channel should be "
          "silent and the L3 Tanimoto detector handles them).",
          f"- **Non-persistence substrate leaks: {len(nonpersist_false)}** "
          "(expected 0 - off-substrate drugs should ABSTAIN).", "",
          "## Per-compound", "",
          "| compound | expected | substrate hypothesis | durable? | engaged (calibrated) | capability flags | abstained (un-calibrated) |",
          "|---|---|---|---|---|---|---|"]
    order = {"ablative_heldout": 0, "senolytic_flavonoid": 1, "capability_heldout": 2, "non_persistence": 3}
    for _, r in out.sort_values("expected_class", key=lambda s: s.map(order)).iterrows():
        Ls.append(f"| {r['compound']} | {r['expected_class']} | {r['substrate_hypothesis']} | "
                  f"{'**yes**' if r['promotes_durable'] else 'no'} | {r['engaged']} | "
                  f"{r['capability_flags']} | {r['abstained_targets']} |")
    Ls += ["", "## Reading", "",
           "Only an ABLATIVE (senolytic) engagement on a calibration-passing target promotes "
           "to durable. capability/window engagements are reported as hypotheses (flags) and "
           "never auto-promoted. Engagements on un-calibrated targets are listed but IGNORED "
           "(the `abstained` column) - the engine refuses to trust a channel MAMMAL cannot "
           "route. This is the structure-computable persistence prior the design doc deferred, "
           "now gated on measured per-target calibration.", ""]
    REPORT.write_text("\n".join(Ls), encoding="utf-8")
    L.info("Wrote %s", REPORT)


if __name__ == "__main__":
    raise SystemExit(main())
