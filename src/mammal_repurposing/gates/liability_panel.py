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


def evaluate_compound_znorm(
    compound_name: str,
    per_target_z: dict[str, float],            # {gene_symbol: within-target z-score}
    panel: pd.DataFrame,
    z_cut_tier1: float = 2.0,
    z_flag_tier2: float = 1.5,
    z_flag_tier3: float = 1.0,
) -> LiabilityResult:
    """Z-norm gate (the §8.0b-zn fix for MAMMAL prior collapse).

    Why: MAMMAL DTI on liability targets has per-target std 0.02-0.17 (noise
    floor) — the absolute pKi thresholds in the panel CSV all sit inside the
    prior mean and trip on every compound. Instead, threshold on within-target
    Z-score: a compound is a "Tier 1 hit" iff it ranks more than +2σ above the
    library mean for that target — i.e. MAMMAL flags it as a within-target
    outlier, not just a prior-mean draw.

    Mirrors §7.18 selectivity rescue logic.
    """
    tier_1_hits: list[str] = []
    tier_2_hits: list[str] = []
    tier_3_hits: list[str] = []
    top_with_tier: list[tuple[str, float, int]] = []

    for _, row in panel.iterrows():
        gene = row["gene_symbol"]
        z = per_target_z.get(gene)
        if z is None or pd.isna(z):
            continue
        tier = int(row["severity_tier"])
        top_with_tier.append((gene, float(z), tier))

        if tier == 1 and z >= z_cut_tier1:
            tier_1_hits.append(gene)
        elif tier == 2 and z >= z_flag_tier2:
            tier_2_hits.append(gene)
        elif tier == 3 and z >= z_flag_tier3:
            tier_3_hits.append(gene)

    if tier_1_hits:
        status = "CUT"
        note = f"Tier 1 within-target z≥{z_cut_tier1}σ: {', '.join(tier_1_hits)}"
    elif len(tier_2_hits) >= 2:
        status = "FLAG"
        note = (f"Composite Tier 2 z≥{z_flag_tier2}σ "
                f"({len(tier_2_hits)} hits): {', '.join(tier_2_hits)}")
    elif len(tier_2_hits) == 1:
        status = "FLAG"
        note = f"Tier 2 z≥{z_flag_tier2}σ: {tier_2_hits[0]}"
    elif tier_3_hits:
        status = "PASS"
        note = f"Informational Tier 3 z≥{z_flag_tier3}σ: {', '.join(tier_3_hits)}"
    else:
        status = "PASS"
        note = "clean"

    top_with_tier.sort(key=lambda t: -t[1])
    top3 = top_with_tier[:3]
    summary = ", ".join(f"{g}=z{p:+.2f}(T{t})" for g, p, t in top3)

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
    znorm: bool = False,
    z_cut_tier1: float = 2.0,
    z_flag_tier2: float = 1.5,
    z_flag_tier3: float = 1.0,
) -> pd.DataFrame:
    """Vectorised gate application.

    Args:
        liability_dti: long-format DTI predictions on the 44-target liability panel
        panel: loaded panel spec from load_panel()
        compound_col / gene_col / pki_col: column names in liability_dti
        znorm: if True, z-score within target and apply z-threshold gating
            (the §8.0b-zn fix for MAMMAL prior collapse). When False uses raw
            predicted_pkd against the CSV's absolute thresholds.

    Returns:
        DataFrame indexed by compound_name with columns:
            liability_status, liability_note, liability_summary,
            tier_1_hits, tier_2_hits, tier_3_hits, top_3_liabilities
    """
    df = liability_dti
    value_col = pki_col

    if znorm:
        # Per-target z-normalisation. Targets with zero/NaN std → all-zero
        # (no signal contribution, no hits).
        # Vectorised z-norm avoids the deprecated groupby.apply(include_groups=)
        # signature entirely.
        df = liability_dti.copy()
        grp = df.groupby(gene_col)[pki_col]
        mu = grp.transform("mean")
        sigma = grp.transform("std")
        # Guard against zero/NaN std (single-point or constant targets).
        z = (df[pki_col] - mu) / sigma
        z = z.where(sigma.notna() & (sigma != 0), 0.0)
        df["_pkd_z"] = z
        value_col = "_pkd_z"

    pivot = df.pivot_table(
        index=compound_col, columns=gene_col, values=value_col, aggfunc="first",
    )

    rows = []
    for compound, row in pivot.iterrows():
        per_target = row.dropna().to_dict()
        if znorm:
            r = evaluate_compound_znorm(
                str(compound), per_target, panel,
                z_cut_tier1=z_cut_tier1,
                z_flag_tier2=z_flag_tier2,
                z_flag_tier3=z_flag_tier3,
            )
        else:
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
                (f"{g}=z{p:+.2f}(T{t})" if znorm else f"{g}={p:.2f}(T{t})")
                for g, p, t in r.top_3_liabilities
            ),
        })

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# §8.13 — Pocket-class-conditioned liability gate
# ---------------------------------------------------------------------------
# Per research/4-tier/archived/Pocket-Conditioned-Boltz2.md §3.3 + V4 plan §8.13:
# the absolute-mode 5-HT2B/hERG/HRH1/CB1 CUT is too aggressive when the
# predicted pose binds OUTSIDE the orthosteric pocket. The literature precedent
# for the demotion logic:
#   - 5-HT2B: only Roth-2007-class agonists driving valvulopathy bind orthosteric;
#     allosteric NAMs don't trigger the warning (Connolly 1997, Schade 2007).
#   - hERG: central-pore Y652/F656/T623 binding is the classical block; allosteric /
#     auxiliary-subunit / vestibule binding is much weaker risk
#     (Dumotier & Urban 2024).
#   - HRH1: peripheral antihistamines bind orthosteric; allosteric H1 ligands have
#     no documented cognition-impairing precedent (Gray 2015).
#   - CB1: Topol 2010 CRESCENDO neuropsych AEs were rimonabant orthosteric only;
#     allosteric NAMs have different clinical risk profile.

POCKET_AWARE_DEMOTABLE: dict[str, set[str]] = {
    # gene -> set of pocket_classes for which CUT can be demoted to FLAG
    "HTR2B": {"allosteric_known", "allosteric_putative"},
    "KCNH2": {"allosteric_known", "allosteric_putative", "no_pocket_match"},
    "HRH1":  {"allosteric_known", "allosteric_putative"},
    "CNR1":  {"allosteric_known", "allosteric_putative"},
    "CHRM1": {"allosteric_known", "allosteric_putative"},
    "OPRM1": {"allosteric_known", "allosteric_putative"},
    "MAOA":  {"allosteric_known", "allosteric_putative"},
}


def pocket_aware_liability_gate(
    liability_gates_df: pd.DataFrame,
    pose_classifications_df: pd.DataFrame,
    compound_col: str = "compound_name",
    gene_col: str = "target_gene",
    pocket_class_col: str = "pocket_class",
) -> pd.DataFrame:
    """Compose §7.5 pocket classifier output with §8.0b liability gates.

    For each compound, if a Tier-1 hard-fail liability hit can be demoted (the
    target's `gene` appears in POCKET_AWARE_DEMOTABLE AND the pose classifies
    into a demotable `pocket_class`), the gate verdict gets reconsidered:

        Tier 1 CUT @ demotable-pocket   →   Tier 2 FLAG @ pocket-aware-demoted
        Tier 1 CUT @ non-demotable      →   Tier 1 CUT (unchanged)
        Tier 2 FLAG                     →   unchanged
        Tier 3 / PASS                   →   unchanged

    The output appends three new columns:
        pocket_demotions  — semicolon-separated list of (gene, original_tier, new_tier)
        n_pocket_demoted  — count of demotions
        liability_status_pocket_aware — new status after demotion

    Args:
        liability_gates_df: output of apply_liability_gates (one row per compound,
            with tier_1_hits / tier_2_hits / tier_3_hits semicolon strings).
        pose_classifications_df: per-(compound, target_gene) classifier output,
            columns: compound_name, target_gene, pocket_class.

    Returns:
        Augmented copy of liability_gates_df.
    """
    out = liability_gates_df.copy()

    # Build (compound, gene) -> pocket_class map (lowercased compound names)
    pose_map: dict[tuple[str, str], str] = {}
    for _, r in pose_classifications_df.iterrows():
        cname = str(r[compound_col]).lower().strip()
        gene = str(r[gene_col])
        pose_map[(cname, gene)] = str(r[pocket_class_col])

    new_status: list[str] = []
    demotions: list[str] = []
    n_demoted: list[int] = []

    for _, r in out.iterrows():
        cname = str(r["compound_name"]).lower().strip()
        original_status = str(r["liability_status"])
        tier1_hits = [g for g in str(r.get("tier_1_hits", "") or "").split(";") if g]
        if original_status != "CUT" or not tier1_hits:
            new_status.append(original_status)
            demotions.append("")
            n_demoted.append(0)
            continue

        demoted: list[tuple[str, str, str]] = []
        remaining_t1: list[str] = []
        for gene in tier1_hits:
            demotable_classes = POCKET_AWARE_DEMOTABLE.get(gene, set())
            pose_class = pose_map.get((cname, gene))
            if pose_class and pose_class in demotable_classes:
                demoted.append((gene, "T1_CUT", "T2_FLAG_pocket_aware"))
            else:
                remaining_t1.append(gene)

        if remaining_t1:
            # Some Tier 1 hits couldn't be demoted → still CUT
            new = "CUT"
        elif demoted:
            # All Tier 1 hits demoted → upgrade to FLAG (pocket-aware)
            new = "FLAG"
        else:
            new = original_status

        new_status.append(new)
        demotions.append("; ".join(f"{g}:{a}→{b}" for g, a, b in demoted))
        n_demoted.append(len(demoted))

    out["liability_status_pocket_aware"] = new_status
    out["pocket_demotions"] = demotions
    out["n_pocket_demoted"] = n_demoted
    return out


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

    def _norm(v):
        # Robust against NaN (float) appearing from outer-merge gaps.
        if v is None:
            return ""
        if isinstance(v, float) and pd.isna(v):
            return ""
        return str(v).upper()

    def _final(row):
        a_s = _norm(row.get("admet_status"))
        l_s = _norm(row.get("liability_status"))
        if a_s == "CUT" or l_s == "CUT":
            return "CUT"
        if a_s == "FLAG" or l_s == "FLAG":
            return "FLAG"
        return "PASS"

    merged["final_status"] = merged.apply(_final, axis=1)
    return merged
