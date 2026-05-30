"""V5 Hypothesis audit — falsifiable re-test of every pre-committed claim
from the V3/V4/V5 design docs against the live production artifacts.

The point is brutal honesty. Each claim either PASSES, DEGRADES, or
FAILS against the current shipped data. No claim is grandfathered.

Claims tested (per V4_STATUS_AND_FORWARD_PLAN.md):
  H1. **Tanimoto beats MAMMAL 7/0** at every audited cognition target
      (reports/pipeline/tanimoto_baseline_v1.md)
  H2. **§7.11 SLC6A3 isotonic post-cal ρ ∈ [+0.45, +0.65]** at fit time
      (paper-commitment from Isotonic-PerTarget-Calibration.md)
  H3. **§7.11 SLC6A2 isotonic post-cal ρ ∈ [+0.30, +0.55]**
  H4. **Positive-control sanity gate**: at least 5 of 7 named compounds
      rank in the top-20% at their expected target
      (donepezil@ACHE, methylphenidate@SLC6A3, atomoxetine@SLC6A2,
       pitolisant@HRH3, encenicline@CHRNA7, rolipram@PDE4D,
       aripiprazole partial@DRD1)
  H5. **Negative-control suppression**: peripheral-only compounds
      (loratadine, naproxen, simvastatin, insulin, warfarin, metformin,
       enalapril, omeprazole, atenolol, ibuprofen) average BELOW the
       library median rrf_score in the V5/V7 fusion
  H6. **Top-25 PASS mechanism diversity**: ≥5 distinct mechanism classes
  H7. **§8.0b-zn pharmacology consistency**: hydroxyzine → KCNH2 CUT;
      aripiprazole → broad polypharmacology CUT; donepezil → PASS;
      methylphenidate → PASS; bupropion → PASS
  H8. **§7.5 pocket DB validation**: 13/13 gates pass (P1+P2+P3)
  H9. **§8.15 disagreement signal real**: at least one
      novel_scaffold_suspect AND one activity_cliff_suspect surface
  H10. **CHRNA7 rescue**: TC-5619 + encenicline both in CHRNA7 top-25
       under Boltz-rescued ranking
  H11. **Calibrator QC drift**: SLC6A3 audit Δρ should be within ±0.20
       of the fitted +0.62 (§8.16 audit; today's run was -0.19 — at boundary)
  H12. **MoA ranker preserves top-10**: adding cluster_b_moa shouldn't
       displace any of d-amph/methylphenidate/bupropion from top-3

Output: reports/pipeline/hypothesis_audit_v1.md + JSON per-claim ledger at
data/results/v2/hypothesis_audit_v1.json.

Exit code: 0 if no claim FAILED, 1 if any FAILED, 2 if any DEGRADED.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("v5_hypothesis")


@dataclass
class HypothesisVerdict:
    id: str
    claim: str
    status: str            # PASS | DEGRADE | FAIL | INSUFFICIENT_DATA
    measured: Any = None
    expected: Any = None
    note: str = ""
    raw: dict = field(default_factory=dict)


# ---------------------------------------------------------------------------
# H1 — Tanimoto-beats-MAMMAL 7/0
# ---------------------------------------------------------------------------
def h1_tanimoto_beats_mammal() -> HypothesisVerdict:
    path = ROOT / "reports" / "pipeline" / "tanimoto_baseline_v1.md"
    if not path.exists():
        return HypothesisVerdict(
            id="H1", claim="Tanimoto ρ beats MAMMAL ρ at every audited cognition target (7/0)",
            status="INSUFFICIENT_DATA",
            note=f"No tanimoto_baseline_v1.md at {path}",
        )
    text = path.read_text(encoding="utf-8")
    # Look for the Δρ table; count rows where Tanimoto wins.
    # Actual table shape:
    #   | uniprot | gene | n | n_actives | mammal_rho | tanimoto_rho | delta | verdict |
    # So cells[4]=mammal_rho, cells[5]=tanimoto_rho.
    wins = 0
    losses = 0
    details: list[tuple[str, float, float]] = []
    for line in text.splitlines():
        if line.startswith("|") and "|" in line[1:]:
            cells = [c.strip() for c in line.split("|")[1:-1]]
            if len(cells) >= 6 and cells[0] not in ("Target", "UniProt", "---", "")\
                    and not cells[0].startswith("-"):
                try:
                    mammal = float(cells[4])
                    tani = float(cells[5].replace("*", ""))
                    if tani > mammal:
                        wins += 1
                    elif tani < mammal:
                        losses += 1
                    details.append((cells[1], mammal, tani))
                except (ValueError, IndexError):
                    continue
    return HypothesisVerdict(
        id="H1",
        claim="Tanimoto ρ beats MAMMAL ρ at every audited cognition target",
        status="PASS" if (wins >= 7 and losses == 0) else
               ("DEGRADE" if wins >= 5 else "FAIL"),
        measured=f"{wins} wins, {losses} losses",
        expected="≥7 wins, 0 losses",
        note="Parsed from reports/pipeline/tanimoto_baseline_v1.md table",
    )


# ---------------------------------------------------------------------------
# H2 / H3 — SLC6A3 / SLC6A2 isotonic ρ predictions
# ---------------------------------------------------------------------------
def _read_router_csv() -> pd.DataFrame | None:
    p = ROOT / "data" / "calibration" / "router_decisions.csv"
    if not p.exists():
        return None
    return pd.read_csv(p)


def h2_slc6a3_post_cal() -> HypothesisVerdict:
    rd = _read_router_csv()
    if rd is None:
        return HypothesisVerdict(
            id="H2", claim="SLC6A3 post-cal ρ ∈ [+0.45, +0.65]",
            status="INSUFFICIENT_DATA", note="No router_decisions.csv")
    row = rd[rd["uniprot"] == "Q01959"]
    if row.empty:
        return HypothesisVerdict(
            id="H2", claim="SLC6A3 post-cal ρ ∈ [+0.45, +0.65]",
            status="INSUFFICIENT_DATA", note="SLC6A3 not in router")
    rho = float(row["post_fit_loco_rho"].iloc[0])
    in_band = 0.45 <= rho <= 0.65
    return HypothesisVerdict(
        id="H2",
        claim="SLC6A3 (DAT) isotonic post-cal ρ ∈ [+0.45, +0.65] (fit time)",
        status="PASS" if in_band else ("DEGRADE" if rho >= 0.30 else "FAIL"),
        measured=f"ρ={rho:+.3f}",
        expected="[+0.45, +0.65]",
        note="From router_decisions.csv at fit time (NOT audit). "
             "H11 tests calibrator drift since fit.",
    )


def h3_slc6a2_post_cal() -> HypothesisVerdict:
    rd = _read_router_csv()
    if rd is None:
        return HypothesisVerdict(
            id="H3", claim="SLC6A2 post-cal ρ ∈ [+0.30, +0.55]",
            status="INSUFFICIENT_DATA")
    row = rd[rd["uniprot"] == "P23975"]
    if row.empty:
        return HypothesisVerdict(
            id="H3", claim="SLC6A2 post-cal ρ ∈ [+0.30, +0.55]",
            status="INSUFFICIENT_DATA")
    rho = float(row["post_fit_loco_rho"].iloc[0])
    in_band = 0.30 <= rho <= 0.55
    return HypothesisVerdict(
        id="H3",
        claim="SLC6A2 (NET) isotonic post-cal ρ ∈ [+0.30, +0.55] (fit time)",
        status="PASS" if in_band else ("DEGRADE" if rho >= 0.20 else "FAIL"),
        measured=f"ρ={rho:+.3f}",
        expected="[+0.30, +0.55]",
    )


# ---------------------------------------------------------------------------
# H4 — Positive control sanity (top-20% at expected target)
# ---------------------------------------------------------------------------
POSITIVE_CONTROLS: dict[str, list[str]] = {
    "P22303": ["donepezil", "rivastigmine", "galantamine"],      # ACHE
    "Q01959": ["methylphenidate", "d-amphetamine", "modafinil"], # SLC6A3
    "P23975": ["atomoxetine", "methylphenidate"],                # SLC6A2
    "Q9Y5N1": ["pitolisant"],                                    # HRH3
    "P36544": ["encenicline", "galantamine", "tc-5619"],         # CHRNA7
    "Q08499": ["rolipram", "bpn14770"],                          # PDE4D
    "P21728": ["aripiprazole"],                                  # DRD1 (partial)
}


def h4_positive_control_sanity() -> HypothesisVerdict:
    # Use the calibrated MAMMAL DTI grid (more relevant than v6 RRF for
    # per-target sanity).
    dti_path = ROOT / "data" / "results" / "dti_scores_calibrated.parquet"
    src_col = "calibrated_pkd"
    if not dti_path.exists():
        dti_path = ROOT / "data" / "results" / "dti_scores.parquet"
        src_col = "predicted_pkd"
    dti = pd.read_parquet(dti_path)

    passes: dict[str, list[str]] = {}
    fails: dict[str, list[str]] = {}
    for uniprot, expected in POSITIVE_CONTROLS.items():
        sub = dti[dti["target_uniprot"] == uniprot].copy()
        if sub.empty:
            continue
        sub = sub.sort_values(src_col, ascending=False).reset_index(drop=True)
        sub["rank"] = sub.index + 1
        sub["top_pct"] = sub["rank"] / len(sub)
        sub["compound_lc"] = sub["compound_name"].str.lower().str.strip()

        target_passes = []
        target_fails = []
        for c in expected:
            row = sub[sub["compound_lc"] == c.lower().strip()]
            if row.empty:
                continue
            top_pct = float(row["top_pct"].iloc[0])
            if top_pct <= 0.20:
                target_passes.append(c)
            else:
                target_fails.append(f"{c}({top_pct:.1%})")
        if target_passes:
            passes[uniprot] = target_passes
        if target_fails:
            fails[uniprot] = target_fails

    n_pass = len(passes)
    n_required = 5    # ≥5 of 7
    status = "PASS" if n_pass >= n_required else \
             ("DEGRADE" if n_pass >= 3 else "FAIL")
    return HypothesisVerdict(
        id="H4",
        claim="≥5 of 7 positive controls in top-20% at expected target "
              f"(via {src_col})",
        status=status,
        measured=f"{n_pass} of {len(POSITIVE_CONTROLS)} targets PASS",
        expected="≥5 / 7",
        note=f"DTI source: {dti_path.name}",
        raw={"passes": passes, "fails": fails},
    )


# ---------------------------------------------------------------------------
# H5 — Negative control suppression
# ---------------------------------------------------------------------------
NEGATIVE_CONTROLS = [
    "loratadine", "naproxen", "simvastatin", "insulin", "warfarin",
    "metformin", "enalapril", "omeprazole", "atenolol", "ibuprofen",
    "cetirizine", "ranitidine",
]


def h5_negative_control_suppression() -> HypothesisVerdict:
    rank_path = ROOT / "data" / "results" / "v2" / "final_ranking_v7_moa.parquet"
    if not rank_path.exists():
        rank_path = ROOT / "data" / "results" / "v2" / "final_ranking_v6_calibrated_znorm.parquet"
    rk = pd.read_parquet(rank_path)
    rk = rk.copy()
    rk["compound_lc"] = rk["compound_name"].str.lower().str.strip()
    median = float(rk["rrf_score"].median())

    in_lib = rk[rk["compound_lc"].isin([n.lower() for n in NEGATIVE_CONTROLS])]
    if in_lib.empty:
        return HypothesisVerdict(
            id="H5", claim="Negative controls average below median rrf_score",
            status="INSUFFICIENT_DATA",
            note=f"No negative controls in {rank_path.name}",
        )
    neg_mean = float(in_lib["rrf_score"].mean())
    suppressed = neg_mean < median
    return HypothesisVerdict(
        id="H5",
        claim="Negative controls (peripheral-only) average BELOW library "
              "median rrf_score",
        status="PASS" if suppressed else "DEGRADE",
        measured=f"neg_mean={neg_mean:.3f}; lib_median={median:.3f}; "
                 f"n_neg_in_lib={len(in_lib)}",
        expected=f"< {median:.3f}",
        note=f"From {rank_path.name}",
        raw={"neg_compound_scores":
             dict(zip(in_lib["compound_name"], in_lib["rrf_score"].round(3)))},
    )


# ---------------------------------------------------------------------------
# H6 — Top-25 PASS mechanism diversity
# ---------------------------------------------------------------------------
def h6_top25_diversity() -> HypothesisVerdict:
    # Use the v7 PASS-only top-25 from combined_gates
    rank_path = ROOT / "data" / "results" / "v2" / "final_ranking_v7_moa.parquet"
    if not rank_path.exists():
        rank_path = ROOT / "data" / "results" / "v2" / "final_ranking_v6_calibrated_znorm.parquet"
    rk = pd.read_parquet(rank_path)
    gates = pd.read_parquet(ROOT / "data" / "results" / "v2" / "combined_gates.parquet")
    g = gates[["compound_name", "final_status"]]
    g["_join"] = g["compound_name"].str.lower().str.strip()
    rk = rk.copy()
    rk["_join"] = rk["compound_name"].str.lower().str.strip()
    j = rk.merge(g.drop(columns=["compound_name"]), on="_join", how="left")
    pass_only = j[j["final_status"] == "PASS"].head(25)

    # Find mechanism class for each compound's best target
    compounds_meta = pd.read_parquet(ROOT / "data" / "interim" / "compounds.parquet")
    meta_map = dict(zip(compounds_meta["name"].str.lower(), compounds_meta["mechanism_class"]))
    classes = set()
    for _, r in pass_only.iterrows():
        cl = meta_map.get(r["compound_name"].lower(), "")
        if cl:
            classes.add(cl)
    n_classes = len(classes)
    return HypothesisVerdict(
        id="H6",
        claim="Top-25 PASS-only set spans ≥5 distinct mechanism classes",
        status="PASS" if n_classes >= 5 else
               ("DEGRADE" if n_classes >= 3 else "FAIL"),
        measured=f"{n_classes} classes: {sorted(classes)}",
        expected="≥5",
    )


# ---------------------------------------------------------------------------
# H7 — §8.0b-zn pharmacology consistency
# ---------------------------------------------------------------------------
EXPECTED_LIABILITY = {
    "hydroxyzine":    ("CUT",  ["KCNH2"]),
    "aripiprazole":   ("CUT",  ["CHRM1", "MAOA"]),
    "risperidone":    ("CUT",  ["HTR2B"]),
    "donepezil":      ("PASS", None),
    "methylphenidate": ("PASS", None),
    "bupropion":      ("PASS", None),
    "rivastigmine":   ("PASS", None),
    "galantamine":    ("PASS", None),
}


def h7_liability_pharmacology() -> HypothesisVerdict:
    gates_path = ROOT / "data" / "results" / "v2" / "liability_gates.parquet"
    if not gates_path.exists():
        return HypothesisVerdict(
            id="H7", claim="§8.0b-zn pharmacology consistent",
            status="INSUFFICIENT_DATA", note="No liability_gates.parquet")
    g = pd.read_parquet(gates_path)
    correct = 0
    wrong: list[str] = []
    for c, (expected_status, expected_genes) in EXPECTED_LIABILITY.items():
        sub = g[g["compound_name"].str.lower().str.strip() == c.lower()]
        if sub.empty:
            wrong.append(f"{c}:MISSING")
            continue
        actual = str(sub["liability_status"].iloc[0])
        if actual != expected_status:
            wrong.append(f"{c}:expected={expected_status} got={actual}")
            continue
        if expected_genes is None:
            correct += 1
            continue
        # For CUT compounds, also verify at least one expected gene is in the
        # tier_1_hits string
        hits = str(sub["tier_1_hits"].iloc[0] or "")
        if any(g_ in hits for g_ in expected_genes):
            correct += 1
        else:
            wrong.append(f"{c}:expected_genes_not_in_T1={hits}")

    return HypothesisVerdict(
        id="H7",
        claim="§8.0b-zn assigns expected status + tier-1 hits for 8 reference compounds",
        status="PASS" if correct >= 7 else
               ("DEGRADE" if correct >= 5 else "FAIL"),
        measured=f"{correct} of {len(EXPECTED_LIABILITY)} correct",
        expected="≥7 of 8",
        raw={"misclassified": wrong},
    )


# ---------------------------------------------------------------------------
# H8 — Pocket DB 13/13 validation
# ---------------------------------------------------------------------------
def h8_pocket_db_gates() -> HypothesisVerdict:
    report = ROOT / "reports" / "pipeline" / "pocket_database_v1.md"
    if not report.exists():
        return HypothesisVerdict(
            id="H8", claim="§7.5 13/13 validation gates pass",
            status="INSUFFICIENT_DATA", note="No pocket_database_v1.md")
    text = report.read_text(encoding="utf-8")
    # Parse per-gate counters. The report uses markdown like "**Gate P1**: 4/4 passed"
    import re
    matches = re.findall(r"Gate\s+P\d\**:\s*\**\s*(\d+)/(\d+)", text)
    if not matches:
        matches = re.findall(r"P\d:\s*(\d+)/(\d+)", text)
    if not matches:
        return HypothesisVerdict(
            id="H8", claim="§7.5 13/13 validation gates pass",
            status="INSUFFICIENT_DATA",
            note="Couldn't parse P1/P2/P3 numbers from report")
    total_pass = sum(int(p) for p, _ in matches)
    total = sum(int(t) for _, t in matches)
    return HypothesisVerdict(
        id="H8",
        claim="§7.5 pocket DB validation: 13/13 gates pass",
        status="PASS" if total_pass == total and total >= 13 else
               ("DEGRADE" if total_pass >= total - 1 else "FAIL"),
        measured=f"{total_pass}/{total}",
        expected="13/13",
    )


# ---------------------------------------------------------------------------
# H9 — §8.15 disagreement signal real
# ---------------------------------------------------------------------------
def h9_disagreement_signal() -> HypothesisVerdict:
    path = ROOT / "data" / "results" / "v2" / "disagreement_signal.parquet"
    if not path.exists():
        return HypothesisVerdict(
            id="H9", claim="§8.15 surfaces novel_scaffold + activity_cliff",
            status="INSUFFICIENT_DATA")
    df = pd.read_parquet(path)
    n_novel = int((df["disagreement_tag"] == "novel_scaffold_suspect").sum())
    n_cliff = int((df["disagreement_tag"] == "activity_cliff_suspect").sum())
    return HypothesisVerdict(
        id="H9",
        claim="§8.15 disagreement signal surfaces ≥1 novel_scaffold AND ≥1 activity_cliff",
        status="PASS" if (n_novel >= 1 and n_cliff >= 1) else
               ("DEGRADE" if (n_novel >= 1 or n_cliff >= 1) else "FAIL"),
        measured=f"novel_scaffold={n_novel}, activity_cliff={n_cliff}",
        expected="≥1 each",
    )


# ---------------------------------------------------------------------------
# H10 — CHRNA7 PAM rescue
# ---------------------------------------------------------------------------
def h10_chrna7_rescue() -> HypothesisVerdict:
    """TC-5619 + encenicline should rank in CHRNA7 top-25 by Boltzina."""
    bp = ROOT / "data" / "results" / "v2" / "boltzina_affinity.parquet"
    if not bp.exists():
        return HypothesisVerdict(
            id="H10", claim="CHRNA7 rescue: TC-5619 + encenicline in Boltz top-25",
            status="INSUFFICIENT_DATA")
    df = pd.read_parquet(bp)
    sub = df[df["target_uniprot"] == "P36544"].copy()
    if sub.empty:
        return HypothesisVerdict(
            id="H10", claim="CHRNA7 rescue",
            status="INSUFFICIENT_DATA",
            note="No CHRNA7 Boltzina rows")
    sub = sub.sort_values("affinity_probability_binary", ascending=False).reset_index(drop=True)
    sub["rank"] = sub.index + 1
    sub["compound_lc"] = sub["compound_name"].str.lower().str.strip()
    pos_ctrls = ["tc-5619", "encenicline"]
    found = {}
    for c in pos_ctrls:
        rows = sub[sub["compound_lc"] == c]
        if not rows.empty:
            found[c] = int(rows["rank"].iloc[0])
    in_top25 = sum(1 for r in found.values() if r <= 25)
    return HypothesisVerdict(
        id="H10",
        claim="CHRNA7 rescue: TC-5619 + encenicline in top-25 by Boltzina affinity",
        status="PASS" if in_top25 >= 2 else
               ("DEGRADE" if in_top25 == 1 else "FAIL"),
        measured=f"found={found}; in_top25={in_top25}/2",
        expected="2/2 in top-25",
    )


# ---------------------------------------------------------------------------
# H11 — Calibrator drift (SLC6A3 audit)
# ---------------------------------------------------------------------------
def h11_calibrator_drift() -> HypothesisVerdict:
    qc = ROOT / "data" / "calibration" / "qc" / "Q01959.json"
    if not qc.exists():
        return HypothesisVerdict(
            id="H11", claim="SLC6A3 calibrator drift |Δρ| ≤ 0.20",
            status="INSUFFICIENT_DATA")
    obj = json.loads(qc.read_text())
    delta = obj.get("delta_rho")
    if delta is None:
        return HypothesisVerdict(
            id="H11", claim="SLC6A3 calibrator drift",
            status="INSUFFICIENT_DATA",
            note="delta_rho is null (no reported_post_cal_rho)")
    abs_d = abs(float(delta))
    return HypothesisVerdict(
        id="H11",
        claim="SLC6A3 (Tier A) calibrator drift |Δρ| ≤ 0.20 between fit and audit",
        status="PASS" if abs_d <= 0.10 else
               ("DEGRADE" if abs_d <= 0.20 else "FAIL"),
        measured=f"Δρ = {delta:+.3f}",
        expected="|Δρ| ≤ 0.10 (PASS) / ≤ 0.20 (DEGRADE) / > 0.20 (FAIL)",
        note=f"Audit n={obj.get('audit_n')} vs fit n={obj.get('calibrator_n_reported')}",
    )


# ---------------------------------------------------------------------------
# H12 — MoA ranker preserves top-3
# ---------------------------------------------------------------------------
def h12_moa_preserves_top3() -> HypothesisVerdict:
    v6 = ROOT / "data" / "results" / "v2" / "final_ranking_v6_calibrated_znorm.parquet"
    v7 = ROOT / "data" / "results" / "v2" / "final_ranking_v7_moa.parquet"
    if not (v6.exists() and v7.exists()):
        return HypothesisVerdict(
            id="H12", claim="MoA ranker preserves top-3",
            status="INSUFFICIENT_DATA")
    v6_top3 = pd.read_parquet(v6).head(3)["compound_name"].str.lower().tolist()
    v7_top3 = pd.read_parquet(v7).head(3)["compound_name"].str.lower().tolist()
    preserved = set(v6_top3) == set(v7_top3)
    return HypothesisVerdict(
        id="H12",
        claim="Adding §8.7 MoA ranker preserves v6 top-3 in v7",
        status="PASS" if preserved else "DEGRADE",
        measured=f"v6={v6_top3}; v7={v7_top3}",
        expected="set equality",
    )


# ---------------------------------------------------------------------------
# V6 phase 3 + V7 + V8 + V6.B.4 hypotheses (this sprint additions)
# ---------------------------------------------------------------------------
def h13_balm_availability_probe() -> HypothesisVerdict:
    """BALM adapter graceful degradation: availability() must never throw,
    must return {'available': bool}; if weights missing, must have 'reason'."""
    try:
        from mammal_repurposing.cluster_a.balm_adapter import availability
        a = availability()
        ok = isinstance(a, dict) and "available" in a
        if not a.get("available") and "reason" not in a:
            return HypothesisVerdict(
                id="H13", claim="BALM adapter availability probe is well-formed",
                status="FAIL",
                measured=str(a)[:80],
                expected="dict with 'available' + (if False, 'reason')",
            )
        return HypothesisVerdict(
            id="H13",
            claim="BALM adapter availability probe is well-formed",
            status="PASS" if ok else "FAIL",
            measured=str(a)[:80],
            expected="dict with 'available' boolean",
        )
    except Exception as e:
        return HypothesisVerdict(
            id="H13", claim="BALM adapter availability probe", status="FAIL",
            measured=f"raised {type(e).__name__}: {e}",
        )


def h14_ot_genetics_cognition_studies() -> HypothesisVerdict:
    """OT Genetics fetcher includes the 5 cognition GWAS studies per
    V6.B.1 Stage 2 spec (Davies + Hill + Savage + Sniekers + UKBB)."""
    try:
        from mammal_repurposing.cluster_d.data_fetchers import COGNITION_GWAS_STUDIES
        canonical = {"GCST006269", "GCST006716", "GCST006250", "GCST005059"}
        present = set(COGNITION_GWAS_STUDIES.keys())
        missing = canonical - present
        if missing:
            return HypothesisVerdict(
                id="H14",
                claim="OT Genetics fetcher includes ≥4 canonical cognition GWAS",
                status="FAIL",
                measured=f"present={sorted(present)}",
                expected=f"superset of {sorted(canonical)}",
            )
        return HypothesisVerdict(
            id="H14",
            claim="OT Genetics fetcher includes ≥4 canonical cognition GWAS",
            status="PASS" if len(COGNITION_GWAS_STUDIES) >= 4 else "DEGRADE",
            measured=f"{len(COGNITION_GWAS_STUDIES)} studies registered",
            expected="≥4 cognition GWAS",
        )
    except Exception as e:
        return HypothesisVerdict(
            id="H14", claim="OT Genetics cognition studies", status="FAIL",
            measured=f"raised {type(e).__name__}: {e}",
        )


def h15_pymc_nuts_stub_path() -> HypothesisVerdict:
    """Stage 0 stub posterior runs without PyMC and produces σ(θ) ∈ (0, 1)
    per target. Falsifier: stub returns nan or values outside (0, 1)."""
    try:
        from mammal_repurposing.cluster_d.bayesian_prior import (
            fit_cluster_d_prior_stub,
        )
        post = fit_cluster_d_prior_stub(
            ["P22303", "Q01959", "P36544"],
            y_ahba={"P22303": 0.5, "Q01959": -0.2, "P36544": 0.1},
            y_l2g={"P22303": 0.7, "Q01959": 0.4},
        )
        all_valid = all(0 < w < 1 for w in post.w_pipeline.values())
        return HypothesisVerdict(
            id="H15",
            claim="Cluster D PyMC stub path produces σ(θ) ∈ (0, 1) per target",
            status="PASS" if all_valid else "FAIL",
            measured=f"w_pipeline={post.w_pipeline}",
            expected="all w ∈ (0, 1)",
        )
    except Exception as e:
        return HypothesisVerdict(
            id="H15", claim="Cluster D NUTS stub path", status="FAIL",
            measured=f"raised {type(e).__name__}: {e}",
        )


def h16_txgnn_per_disease_api() -> HypothesisVerdict:
    """TxGNN per-disease ranking API: 5 cognition anchors registered,
    availability probe returns dict, aggregate_per_compound preserves
    weighted mean."""
    try:
        from mammal_repurposing.cluster_c.txgnn import (
            COGNITION_ANCHORS, availability, aggregate_per_compound,
        )
        anchors_ok = len(COGNITION_ANCHORS) == 5
        avail = availability()
        avail_ok = isinstance(avail, dict) and "available" in avail
        # Smoke aggregate_per_compound
        long = pd.DataFrame([
            {"compound_id": "X", "anchor_id": "A1", "anchor_name": "n",
             "weight": 1.0, "p_indication": 0.8, "p_contraindication": 0.1},
            {"compound_id": "X", "anchor_id": "A2", "anchor_name": "n",
             "weight": 1.0, "p_indication": 0.6, "p_contraindication": 0.2},
        ])
        agg = aggregate_per_compound(long)
        wmean_ok = abs(agg.iloc[0]["mean_p_indication"] - 0.7) < 1e-6
        if anchors_ok and avail_ok and wmean_ok:
            return HypothesisVerdict(
                id="H16",
                claim="TxGNN per-disease API: 5 anchors + availability + weighted mean",
                status="PASS",
                measured=f"5 anchors, avail={avail.get('available')}, wmean=0.70",
            )
        return HypothesisVerdict(
            id="H16",
            claim="TxGNN per-disease API: 5 anchors + availability + weighted mean",
            status="FAIL",
            measured=f"anchors_ok={anchors_ok}, avail_ok={avail_ok}, wmean_ok={wmean_ok}",
        )
    except Exception as e:
        return HypothesisVerdict(
            id="H16", claim="TxGNN per-disease API", status="FAIL",
            measured=f"raised {type(e).__name__}: {e}",
        )


def h17_pbpk_pet_anchors_finite() -> HypothesisVerdict:
    """V7.1 PBPK PET-anchor residuals are finite + bounded for all 3 anchors
    (donepezil cortical AChE / MPH DAT / haloperidol D2)."""
    try:
        from mammal_repurposing.translation.pbpk import (
            compute_pet_anchor_residuals, PET_ANCHORS,
        )
        if len(PET_ANCHORS) < 3:
            return HypothesisVerdict(
                id="H17", claim="V7.1 PBPK PET anchors finite",
                status="FAIL",
                measured=f"only {len(PET_ANCHORS)} anchors",
                expected="≥3 PET anchors",
            )
        residuals = compute_pet_anchor_residuals()
        all_finite = all(np.isfinite(r["residual"]) for r in residuals)
        bounded = all(abs(r["predicted"]) <= 1.0 for r in residuals)
        if all_finite and bounded:
            return HypothesisVerdict(
                id="H17",
                claim="V7.1 PBPK PET anchors produce finite bounded residuals",
                status="PASS",
                measured=f"residuals={[round(r['residual'], 3) for r in residuals]}",
            )
        return HypothesisVerdict(
            id="H17",
            claim="V7.1 PBPK PET anchors produce finite bounded residuals",
            status="FAIL",
            measured=f"finite={all_finite}, bounded={bounded}",
        )
    except Exception as e:
        return HypothesisVerdict(
            id="H17", claim="V7.1 PBPK PET anchors", status="FAIL",
            measured=f"raised {type(e).__name__}: {e}",
        )


def h18_prisma_priors_12_classes() -> HypothesisVerdict:
    """V7.2 PRISMA priors: exactly 12 mechanism classes, all with
    peak_subdomain_g ≤ Roberts 2020 ceiling (0.50)."""
    try:
        from mammal_repurposing.translation.prisma_priors import (
            PRISMA_CLASS_PRIORS,
        )
        n_classes = len(PRISMA_CLASS_PRIORS)
        n_violations = sum(1 for cp in PRISMA_CLASS_PRIORS
                           if cp.peak_subdomain_g > 0.50)
        if n_classes == 12 and n_violations == 0:
            return HypothesisVerdict(
                id="H18",
                claim="V7.2 PRISMA: 12 classes, all peak_g ≤ Roberts ceiling 0.50",
                status="PASS",
                measured=f"n_classes={n_classes}, n_violations=0",
            )
        return HypothesisVerdict(
            id="H18",
            claim="V7.2 PRISMA: 12 classes, all peak_g ≤ Roberts ceiling 0.50",
            status="FAIL",
            measured=f"n_classes={n_classes}, n_violations={n_violations}",
            expected="n_classes=12, n_violations=0",
        )
    except Exception as e:
        return HypothesisVerdict(
            id="H18", claim="V7.2 PRISMA 12 classes", status="FAIL",
            measured=f"raised {type(e).__name__}: {e}",
        )


def h19_lincs_availability_well_formed() -> HypothesisVerdict:
    """V8.1 LINCS availability probe returns dict with cell_lines_supported
    including the 3 neural lines (NPC, NEU, SHSY5Y)."""
    try:
        from mammal_repurposing.cluster_e.ingest_lincs import (
            availability, COGNITION_CELL_LINES,
        )
        a = availability()
        ok = (isinstance(a, dict)
              and "cell_lines_supported" in a
              and {"NPC", "NEU", "SHSY5Y"}.issubset(
                  set(a["cell_lines_supported"])))
        neural_weighted = all(COGNITION_CELL_LINES[c] == 1.0
                              for c in ("NPC", "NEU", "SHSY5Y"))
        if ok and neural_weighted:
            return HypothesisVerdict(
                id="H19",
                claim="V8.1 LINCS: 3 neural cell lines weighted 1.0; probe well-formed",
                status="PASS",
                measured="NPC/NEU/SHSY5Y all weight=1.0",
            )
        return HypothesisVerdict(
            id="H19",
            claim="V8.1 LINCS: 3 neural cell lines weighted 1.0; probe well-formed",
            status="FAIL",
            measured=f"probe_ok={ok}, weighted={neural_weighted}",
        )
    except Exception as e:
        return HypothesisVerdict(
            id="H19", claim="V8.1 LINCS availability", status="FAIL",
            measured=f"raised {type(e).__name__}: {e}",
        )


def h20_jumpcp_13_sources_known() -> HypothesisVerdict:
    """V8.1b JUMP-CP: 13 source codes registered (s1-s13) per
    Chandrasekaran 2024 design + 3 embedding types."""
    try:
        from mammal_repurposing.cluster_e.ingest_jumpcp import (
            JUMP_SOURCES, EMBEDDING_TYPES,
        )
        n_sources = len(JUMP_SOURCES)
        n_embed = len(EMBEDDING_TYPES)
        sources_ok = all(f"s{i}" in JUMP_SOURCES for i in range(1, 14))
        embed_ok = {"deepprofiler", "cellprofiler", "dinov2"}.issubset(
            set(EMBEDDING_TYPES.keys()))
        if n_sources == 13 and n_embed == 3 and sources_ok and embed_ok:
            return HypothesisVerdict(
                id="H20",
                claim="V8.1b JUMP-CP: 13 sources × 3 embedding types registered",
                status="PASS",
                measured=f"sources={n_sources}, embeddings={n_embed}",
            )
        return HypothesisVerdict(
            id="H20",
            claim="V8.1b JUMP-CP: 13 sources × 3 embedding types registered",
            status="FAIL",
            measured=f"sources={n_sources}, embeddings={n_embed}, "
                     f"sources_ok={sources_ok}, embed_ok={embed_ok}",
        )
    except Exception as e:
        return HypothesisVerdict(
            id="H20", claim="V8.1b JUMP-CP coverage", status="FAIL",
            measured=f"raised {type(e).__name__}: {e}",
        )


def h21_validation_gates_framework() -> HypothesisVerdict:
    """V6.B.4 validation gates: 15 reference compounds + 4 gates + Gate 1
    correctly identifies a Roberts-ceiling violation."""
    try:
        from mammal_repurposing.cluster_d.validation_gates import (
            REFERENCE_COMPOUND_SMD, gate_1_roberts_ceiling, availability,
        )
        avail = availability()
        n_compounds = avail["n_reference_compounds"]
        n_gates = avail["n_gates"]
        # Synthetic violation: predict 0.65 at one target
        result = gate_1_roberts_ceiling(
            target_smd_predictions={"safe": 0.30, "bad": 0.65}, ceiling=0.50,
        )
        gate_fires = result.pass_status == "FAIL"
        if n_compounds == 15 and n_gates == 4 and gate_fires:
            return HypothesisVerdict(
                id="H21",
                claim="V6.B.4: 15 ref compounds + 4 gates + Gate 1 fires on g > 0.50",
                status="PASS",
                measured=f"n_compounds={n_compounds}, gate1='{result.pass_status}'",
            )
        return HypothesisVerdict(
            id="H21",
            claim="V6.B.4: 15 ref compounds + 4 gates + Gate 1 fires on g > 0.50",
            status="FAIL",
            measured=f"compounds={n_compounds}, gates={n_gates}, fires={gate_fires}",
        )
    except Exception as e:
        return HypothesisVerdict(
            id="H21", claim="V6.B.4 validation gates", status="FAIL",
            measured=f"raised {type(e).__name__}: {e}",
        )


def h22_effect_size_cluster_d_gate() -> HypothesisVerdict:
    """V7.3 effect-size stub: Cluster D multiplicative gate produces
    monotonic dependence on relevance_post_mean."""
    try:
        from mammal_repurposing.translation.effect_size_model import (
            EffectSizeObservation, fit_effect_size_stub,
        )
        obs = [
            EffectSizeObservation("high", "AChE-I", "P22303",
                                  pchembl_post_mean=8.5,
                                  relevance_post_mean=0.90),
            EffectSizeObservation("low", "AChE-I", "P22303",
                                  pchembl_post_mean=8.5,
                                  relevance_post_mean=0.10),
        ]
        post = fit_effect_size_stub(obs)
        monotonic = post.g_mean["high"] >= post.g_mean["low"]
        return HypothesisVerdict(
            id="H22",
            claim="V7.3: Cluster D multiplicative gate monotonic in relevance_post",
            status="PASS" if monotonic else "FAIL",
            measured=f"g(high)={post.g_mean['high']:.3f}, "
                     f"g(low)={post.g_mean['low']:.3f}",
        )
    except Exception as e:
        return HypothesisVerdict(
            id="H22", claim="V7.3 Cluster D gate monotonicity", status="FAIL",
            measured=f"raised {type(e).__name__}: {e}",
        )


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------
def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path,
                        default=ROOT / "data" / "results" / "v2" / "hypothesis_audit_v1.json")
    parser.add_argument("--report", type=Path,
                        default=ROOT / "reports" / "pipeline" / "hypothesis_audit_v1.md")
    args = parser.parse_args()

    tests = [
        h1_tanimoto_beats_mammal,
        h2_slc6a3_post_cal,
        h3_slc6a2_post_cal,
        h4_positive_control_sanity,
        h5_negative_control_suppression,
        h6_top25_diversity,
        h7_liability_pharmacology,
        h8_pocket_db_gates,
        h9_disagreement_signal,
        h10_chrna7_rescue,
        h11_calibrator_drift,
        h12_moa_preserves_top3,
        # V6 phase 3 + V7 + V8 + V6.B.4 falsifiers (this sprint)
        h13_balm_availability_probe,
        h14_ot_genetics_cognition_studies,
        h15_pymc_nuts_stub_path,
        h16_txgnn_per_disease_api,
        h17_pbpk_pet_anchors_finite,
        h18_prisma_priors_12_classes,
        h19_lincs_availability_well_formed,
        h20_jumpcp_13_sources_known,
        h21_validation_gates_framework,
        h22_effect_size_cluster_d_gate,
    ]

    verdicts = [t() for t in tests]
    counts = {"PASS": 0, "DEGRADE": 0, "FAIL": 0, "INSUFFICIENT_DATA": 0}
    for v in verdicts:
        counts[v.status] = counts.get(v.status, 0) + 1
        logger.info("  %s [%s]: %s", v.id, v.status, v.claim[:90])

    # JSON ledger
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps([asdict(v) for v in verdicts], indent=2, default=str),
        encoding="utf-8",
    )

    # Markdown report
    L: list[str] = []
    L.append("# Hypothesis Audit v1 — V5 Pre-Commitment Validation")
    L.append("")
    L.append("Falsifiable re-test of every pre-committed claim from V3/V4/V5 design "
             "docs against the live production artifacts. **Brutal honesty mode** — "
             "no claim is grandfathered.")
    L.append("")
    L.append(f"**Summary**: PASS={counts['PASS']} | DEGRADE={counts['DEGRADE']} | "
             f"FAIL={counts['FAIL']} | INSUFFICIENT_DATA={counts['INSUFFICIENT_DATA']} "
             f"(of {len(verdicts)} hypotheses)")
    L.append("")
    L.append("## Verdicts")
    L.append("")
    L.append("| ID | Claim | Status | Measured | Expected |")
    L.append("|---|---|---|---|---|")
    for v in verdicts:
        L.append(f"| {v.id} | {v.claim[:80]} | **{v.status}** | {str(v.measured)[:60]} | {str(v.expected)[:40]} |")
    L.append("")

    L.append("## Detail per hypothesis")
    L.append("")
    for v in verdicts:
        L.append(f"### {v.id} — {v.status}")
        L.append("")
        L.append(f"**Claim**: {v.claim}")
        L.append("")
        if v.expected is not None:
            L.append(f"**Expected**: {v.expected}")
        L.append(f"**Measured**: {v.measured}")
        L.append("")
        if v.note:
            L.append(f"_Note_: {v.note}")
            L.append("")
        if v.raw:
            L.append("_Raw_:")
            L.append("```json")
            L.append(json.dumps(v.raw, indent=2, default=str))
            L.append("```")
            L.append("")

    L.append("---")
    L.append("")
    L.append("Generated by `scripts/41_v5_hypothesis_audit.py`. JSON ledger at "
             f"`{args.out.relative_to(ROOT)}`.")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(L), encoding="utf-8")
    logger.info("Wrote %s and %s", args.out, args.report)

    if counts["FAIL"] > 0:
        return 1
    if counts["DEGRADE"] > 0:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
