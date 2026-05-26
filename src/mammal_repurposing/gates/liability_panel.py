"""§8.0b — Cognition-context-tailored 44-target off-target liability gate.

Per research/4-tier/Cognition-44Target-Liability-Panel.md. Loads the 44-row
data/raw/targets_liability_seed.csv with tier-stratified pKi thresholds,
applies them to MAMMAL DTI predictions on the liability panel, and emits
{PASS, FLAG, CUT} verdicts with mechanistic notes.

Runs AFTER ADMET-AI gates, BEFORE RRF fusion. Hard gate precedence:
    final_status = CUT  if (admet_status == CUT or liability_status == CUT)
                   else FLAG if any FLAG
                   else PASS

The regulatory bypass for approved drugs (currently in ADMET) does NOT
extend to liability — aripiprazole is approved for schizophrenia but its
5-HT2B + α1 + D2 profile makes it inappropriate for chronic healthy-adult
cognitive enhancement. The §8.0b panel is the discrimination that ADMET
alone cannot deliver.

References:
  Bowes 2012 Nat Rev Drug Discov 11:909
  Brennan 2024 Nat Rev Drug Discov 23:525
  Dumotier & Urban 2024 J Pharmacol Toxicol Methods 128:107542 (>10x Cmax/Ki)
  Connolly 1997 NEJM 337:581 (fen-phen valvulopathy)
  Schade 2007 NEJM 356:29 (pergolide valvulopathy)
  Roth 2007 NEJM 356:6 (5-HT2B class warning)
  Gray 2015 JAMA Intern Med 175:401 (anticholinergic dementia HR 1.54)
  Coupland 2019 JAMA Intern Med 179:1084 (anticholinergic ~50% dementia risk)
  Topol 2010 Lancet 376:517 (CRESCENDO rimonabant 32% neuropsych AEs)
  Sharretts 2020 NEJM 383:1000 (lorcaserin cancer withdrawal)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Result type
# ---------------------------------------------------------------------------
@dataclass
class LiabilityResult:
    compound_name: str
    status: str                          # PASS | FLAG | CUT
    tier_1_hits: list[str]
    tier_2_hits: list[str]
    tier_3_hits: list[str]
    top_3_liabilities: list[tuple[str, float, int]]    # (gene, pKi, tier)
    liability_note: str
    liability_summary: str               # e.g. "α1A=6.8(T2), 5-HT2B=5.2(clean), hERG=6.4(T1)"


# ---------------------------------------------------------------------------
# Loader
# ---------------------------------------------------------------------------
def load_panel(panel_csv: Path | str) -> pd.DataFrame:
    """Load the 44-target liability panel spec from CSV."""
    df = pd.read_csv(panel_csv)
    # Coerce thresholds to float; treat empty as NaN
    for col in ("cut_threshold_pki", "flag_threshold_pki"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


# ---------------------------------------------------------------------------
# Per-compound gating logic
# ---------------------------------------------------------------------------
def evaluate_compound(
    compound_name: str,
    per_target_pki: dict[str, float],          # {gene_symbol: predicted_pki}
    panel: pd.DataFrame,
) -> LiabilityResult:
    """Apply tier-stratified gates for one compound."""
    tier_1_hits: list[str] = []
    tier_2_hits: list[str] = []
    tier_3_hits: list[str] = []
    top_with_tier: list[tuple[str, float, int]] = []

    for _, row in panel.iterrows():
        gene = row["gene_symbol"]
        pki = per_target_pki.get(gene)
        if pki is None or pd.isna(pki):
            continue
        tier = int(row["severity_tier"])
        cut_t = row["cut_threshold_pki"]
        flag_t = row["flag_threshold_pki"]
        top_with_tier.append((gene, float(pki), tier))

        if tier == 1 and not pd.isna(cut_t) and pki >= cut_t:
            tier_1_hits.append(gene)
        elif tier == 2 and not pd.isna(flag_t) and pki >= flag_t:
            tier_2_hits.append(gene)
        elif tier == 3 and not pd.isna(flag_t) and pki >= flag_t:
            tier_3_hits.append(gene)

    # Verdict
    if tier_1_hits:
        status = "CUT"
        note = f"Tier 1 hard fail: {', '.join(tier_1_hits)}"
    elif len(tier_2_hits) >= 2:
        status = "FLAG"
        note = f"Composite Tier 2 ({len(tier_2_hits)} hits): {', '.join(tier_2_hits)}"
    elif len(tier_2_hits) == 1:
        status = "FLAG"
        note = f"Tier 2: {tier_2_hits[0]}"
    elif tier_3_hits:
        status = "PASS"
        note = f"Informational Tier 3: {', '.join(tier_3_hits)}"
    else:
        status = "PASS"
        note = "clean"

    # Build a one-line summary of the top-3 liabilities by pKi
    top_with_tier.sort(key=lambda t: -t[1])
    top3 = top_with_tier[:3]
    summary = ", ".join(f"{g}={p:.2f}(T{t})" for g, p, t in top3)

    return LiabilityResult(
        compound_name=compound_name,
        status=status,
        tier_1_hits=tier_1_hits,
        tier_2_hits=tier_2_hits,
        tier_3_hits=tier_3_hits,
        top_3_liabilities=top3,
        liability_note=note,
        liability_summary=summary,
    )


def apply_liability_gates(
    liability_dti: pd.DataFrame,
    panel: pd.DataFrame,
    compound_col: str = "compound_name",
    gene_col: str = "target_gene",
    pki_col: str = "predicted_pkd",         # MAMMAL output we treat as pKi
) -> pd.DataFrame:
    """Vectorised gate application.

    Args:
        liability_dti: long-format DTI predictions on the 44-target liability panel
        panel: loaded panel spec from load_panel()
        compound_col / gene_col / pki_col: column names in liability_dti

    Returns:
        DataFrame indexed by compound_name with columns:
            liability_status, liability_note, liability_summary,
            tier_1_hits, tier_2_hits, tier_3_hits, top_3_liabilities
    """
    # Pivot to {compound: {gene: pki}}
    pivot = liability_dti.pivot_table(
        index=compound_col, columns=gene_col, values=pki_col, aggfunc="first",
    )

    rows = []
    for compound, row in pivot.iterrows():
        per_target = row.dropna().to_dict()
        r = evaluate_compound(str(compound), per_target, panel)
        rows.append({
            "compound_name": r.compound_name,
            "liability_status": r.status,
            "liability_note": r.liability_note,
            "liability_summary": r.liability_summary,
            "tier_1_hits": ";".join(r.tier_1_hits),
            "tier_2_hits": ";".join(r.tier_2_hits),
            "tier_3_hits": ";".join(r.tier_3_hits),
            "n_tier_1": len(r.tier_1_hits),
            "n_tier_2": len(r.tier_2_hits),
            "n_tier_3": len(r.tier_3_hits),
            "top_3_liabilities": "; ".join(
                f"{g}={p:.2f}(T{t})" for g, p, t in r.top_3_liabilities
            ),
        })

    return pd.DataFrame(rows)


def combine_admet_and_liability(
    admet_gates: pd.DataFrame,
    liability_gates: pd.DataFrame,
    compound_col: str = "compound_name",
) -> pd.DataFrame:
    """Merge ADMET and liability verdicts into a final_status per compound.

    Precedence: CUT > FLAG > PASS. Compound is CUT if EITHER cluster CUT.
    """
    a = admet_gates[[compound_col, "gate_status"]].rename(
        columns={"gate_status": "admet_status"}
    )
    l = liability_gates[
        [compound_col, "liability_status", "liability_note",
         "liability_summary", "tier_1_hits", "tier_2_hits", "tier_3_hits",
         "n_tier_1", "n_tier_2", "n_tier_3"]
    ]
    merged = a.merge(l, on=compound_col, how="outer")

    def _final(row):
        a_s = (row.get("admet_status") or "").upper()
        l_s = (row.get("liability_status") or "").upper()
        if a_s == "CUT" or l_s == "CUT":
            return "CUT"
        if a_s == "FLAG" or l_s == "FLAG":
            return "FLAG"
        return "PASS"

    merged["final_status"] = merged.apply(_final, axis=1)
    return merged
