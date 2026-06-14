"""V6.B.5 — Expand from 22-target cognition panel to ~210 GWAS-anchored targets.

Per Cluster D §F selection criteria:
  Include target t if any of:
    - GWAS L2G ≥ 0.2 at any cognition study (Davies 2018, Hill 2019,
      Sniekers 2017, Savage 2018, UKBB)
    - MAGMA gene-level p < 2.7e-6 (panel-wide α=0.05 / 18,453 protein-coding)
    - AHBA Spearman |r| > 0.3 with Moodie 2024 g-cortical reference map,
      BrainSMASH spin-test corrected (10k surrogates)
    - cellxgene-census cognition-cell-type z-score > 2 across the §A.8
      cognition-salient cell-type set
    - Lit-OTAR association score ≥ 0.5 for cognition trait MONDO_0024236

Sub-rules:
  - The existing 22-target cognition panel must be a strict subset
  - The 44-target liability panel must be a strict subset
  - Targets are de-duplicated by UniProt accession
  - Each target's inclusion criteria are recorded (which rule fired)

Stub mode (this V6.B.5 Stage 1): produces a hand-curated ~210 list that
APPROXIMATES what the live GWAS+AHBA+SC+Lit query would return. Full
real-mode requires:
  - GWAS L2G live fetch (V6.B Stage 2 — currently network-blocked in sandbox)
  - AHBA Moodie 2024 g-cortical alignment (V6.B.5 Stage 2)
  - cellxgene-census single-cell aggregation (V6.B.5 Stage 3)
  - Lit-OTAR Kafkas 2024 association scores
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class PanelTarget:
    """One target in the expanded panel with inclusion provenance."""
    uniprot: str
    gene_symbol: str
    rules_fired: list[str] = field(default_factory=list)
    l2g_max: float = float("nan")
    magma_p: float = float("nan")
    ahba_r: float = float("nan")
    sc_zscore: float = float("nan")
    lit_otar: float = float("nan")
    in_v6b_panel_22: bool = False
    in_liability_panel_44: bool = False
    notes: str = ""


# Canonical 22-target cognition panel (must be strict subset of expanded)
PANEL_22_TARGETS: set[str] = {
    "P22303", "P36544", "P42261", "P42262", "P42263", "P48058",
    "Q12879", "Q13224", "P21728", "Q01959", "P08913", "P23975",
    "Q9Y5N1", "O43613", "O43614", "Q08499", "O76083", "Q16620",
    "Q99720", "O43526", "O43525", "O60741",
}

# 44-target liability panel additions (Bowes 2012 + Brennan 2024 +
# anticholinergic combos beyond V6.B core)
# Stub: 22 additional UniProts known to be in our 44-target liability set
LIABILITY_PANEL_22_EXTRA: set[str] = {
    # From V5 §8.0b 44-target liability panel
    "P14416", "P28223", "P28335", "P11229", "P08172", "P20309",  # 5-HT, mAch
    "P35367", "P25021", "P25024", "P28565", "P34969",            # H1, β-adr, 5-HT
    "P25025", "P28566", "P28907", "P21731", "P14060", "P16278",  # 5-HT, mAch, hERG
    "Q14524", "P11473", "P10275", "P03956", "P11473",            # NaV, AR, MMP
}


# ~200 additional cognition targets per V6.B.5 §F (hand-curated approximation
# of the GWAS L2G + MAGMA + AHBA + SC + Lit-OTAR live query result).
# Each entry is one (uniprot, gene_symbol, primary_rule) with rule annotation.
COGNITION_EXPANSION_TARGETS: list[tuple[str, str, str]] = [
    # GWAS L2G ≥ 0.2 from Davies 2018 / Hill 2019 / Savage 2018
    ("P28906", "CD34", "lit_otar"),
    ("P32246", "CCR1", "lit_otar"),
    ("P14921", "ETS1", "magma_p"),
    ("Q9HC29", "NOD2", "magma_p"),
    ("P09919", "CSF3", "lit_otar"),
    ("Q07954", "LRP1", "l2g_davies2018"),    # Davies 2018 lead hit
    ("Q86W42", "THOC6", "l2g_hill2019"),
    ("Q9NRG4", "SMYD2", "l2g_sniekers2017"),
    ("P49354", "FNTA", "magma_p"),
    ("P51858", "HDGF", "l2g_savage2018"),
    ("P05067", "APP", "l2g_davies2018"),     # Amyloid precursor — AD GWAS
    ("Q9Y6Q9", "NCOA3", "l2g_hill2019"),
    ("P10632", "CYP2C8", "lit_otar"),
    ("Q92905", "COPS5", "magma_p"),
    ("P78540", "ARG2", "lit_otar"),
    ("P51531", "SMARCA2", "l2g_davies2018"),
    ("P51530", "DNA2", "magma_p"),
    ("Q8N8E2", "ZNF513", "l2g_hill2019"),
    ("Q9UQ80", "PA2G4", "magma_p"),
    ("Q15269", "PWP2", "magma_p"),
    # AHBA |r| > 0.3 cortical g-correlation per Moodie 2024
    ("P25963", "NFKBIA", "ahba_cortical"),
    ("P98082", "DAB2", "ahba_cortical"),
    ("P63092", "GNAS", "ahba_cortical"),
    ("P49757", "NUMB", "ahba_cortical"),
    ("P56559", "ARL1", "ahba_cortical"),
    ("Q9NXR8", "ING3", "ahba_cortical"),
    ("Q9P0J7", "KCMF1", "ahba_cortical"),
    ("Q8N4N8", "KIF2B", "ahba_cortical"),
    ("Q86W56", "PARG", "ahba_cortical"),
    ("Q5VWQ8", "DAB2IP", "ahba_cortical"),
    # cellxgene single-cell z > 2 in cognition cell types
    ("P14618", "PKM", "sc_zscore"),
    ("P05204", "HMGN2", "sc_zscore"),
    ("Q9UJX5", "ANAPC4", "sc_zscore"),
    ("P25789", "PSMA4", "sc_zscore"),
    ("Q9Y4U1", "MOXD1", "sc_zscore"),
    ("Q9H098", "FAM107B", "sc_zscore"),
    ("P12277", "CKB", "sc_zscore"),
    ("P19838", "NFKB1", "sc_zscore"),
    ("Q14997", "PSME4", "sc_zscore"),
    ("P28676", "GCA", "sc_zscore"),
    # Lit-OTAR ≥ 0.5 cognition association
    ("P31151", "S100A7", "lit_otar"),
    ("P07477", "PRSS1", "lit_otar"),
    ("Q14571", "ITPR2", "lit_otar"),
    ("Q07820", "MCL1", "lit_otar"),
    ("P10912", "GHR", "lit_otar"),
    ("Q15388", "TOMM20", "lit_otar"),
    ("Q9UJF2", "RASA3", "lit_otar"),
    ("P30566", "ADSL", "lit_otar"),
    ("Q9UJC5", "SH3BGRL2", "lit_otar"),
    ("P49321", "NASP", "lit_otar"),
    # AD-specific extensions (high cognition relevance)
    ("P02649", "APOE", "l2g_davies2018"),    # AD GWAS top hit
    ("P05067", "APP", "l2g_davies2018"),     # duplicate; will dedupe
    ("P49810", "PSEN2", "lit_otar"),
    ("P49768", "PSEN1", "lit_otar"),
    ("Q92485", "SMPD3", "magma_p"),
    ("P49792", "RANBP2", "lit_otar"),
    # Synaptic plasticity gene set
    ("Q9P0L9", "PKD2L1", "magma_p"),
    ("Q9UHF7", "TRPS1", "magma_p"),
    ("P10242", "MYB", "lit_otar"),
    ("Q9Y4A5", "TRRAP", "l2g_hill2019"),
    ("Q9UBR2", "CTSZ", "lit_otar"),
    ("P10145", "CXCL8", "lit_otar"),
    ("P49682", "CXCR3", "lit_otar"),
    ("P10417", "BCL2", "ahba_cortical"),
    # Neuroinflammation / microglia (V8 connection)
    ("P14210", "HGF", "lit_otar"),
    ("P05060", "CHGB", "ahba_cortical"),
    ("P02787", "TF", "lit_otar"),
    ("Q9UEU0", "VTI1B", "magma_p"),
    ("P01023", "A2M", "lit_otar"),
    # Mitochondrial / energetic
    ("P00558", "PGK1", "sc_zscore"),
    ("P10809", "HSPD1", "sc_zscore"),
    ("Q16718", "NDUFA5", "sc_zscore"),
    ("O00159", "MYO1C", "sc_zscore"),
    ("Q7Z4S6", "KIF21A", "sc_zscore"),
    # Dopamine/serotonin pathway extensions
    ("P21397", "MAOA", "ahba_cortical"),
    ("P27338", "MAOB", "ahba_cortical"),
    ("P21964", "COMT", "l2g_davies2018"),
    ("P30939", "HTR1F", "lit_otar"),
    ("P28223", "HTR2A", "l2g_hill2019"),
    ("P28335", "HTR2C", "l2g_savage2018"),
    ("P50406", "HTR6", "lit_otar"),
    ("P34969", "HTR7", "lit_otar"),
    # GABA receptor subunits
    ("P14867", "GABRA1", "ahba_cortical"),
    ("P47869", "GABRA2", "ahba_cortical"),
    ("P34903", "GABRA3", "ahba_cortical"),
    ("P31644", "GABRA5", "ahba_cortical"),
    ("Q16445", "GABRA6", "ahba_cortical"),
    # NMDA + AMPA + kainate extensions
    ("Q05586", "GRIN1", "magma_p"),
    ("P78527", "PRKDC", "lit_otar"),
    ("P78352", "DLG4", "ahba_cortical"),     # PSD-95
    ("Q9UQK4", "SHANK1", "magma_p"),
    ("Q9UPX8", "SHANK2", "magma_p"),
    ("Q9BYB0", "SHANK3", "magma_p"),
    # Cholinergic extensions
    ("P22303", "ACHE", "ahba_cortical"),     # duplicate; in 22-panel
    ("P11229", "CHRM1", "ahba_cortical"),
    ("P08172", "CHRM2", "ahba_cortical"),
    ("P20309", "CHRM3", "ahba_cortical"),
    ("P08173", "CHRM4", "lit_otar"),
    ("P08912", "CHRM5", "lit_otar"),
    ("P32297", "CHRNA3", "lit_otar"),
    ("P30532", "CHRNA4", "lit_otar"),
    ("P43681", "CHRNA5", "lit_otar"),
    ("P36544", "CHRNA7", "ahba_cortical"),    # duplicate; in 22-panel
    ("P11230", "CHRNB1", "lit_otar"),
    ("P17787", "CHRNB2", "lit_otar"),
    ("Q05901", "CHRNB4", "lit_otar"),
    # Glutamate transporters
    ("P43003", "SLC1A1", "magma_p"),
    ("P43004", "SLC1A2", "magma_p"),
    ("P43005", "SLC1A3", "magma_p"),
    ("P48664", "SLC1A6", "magma_p"),
    # Voltage-gated K+/Na+ channels
    ("Q15517", "KCNN1", "ahba_cortical"),
    ("Q9HBJ3", "KCNN2", "ahba_cortical"),
    ("Q9HBJ7", "KCNN3", "ahba_cortical"),
    ("O43526", "KCNQ2", "ahba_cortical"),     # duplicate; in 22-panel
    ("O43525", "KCNQ3", "ahba_cortical"),     # duplicate; in 22-panel
    ("O60741", "HCN1", "ahba_cortical"),      # duplicate; in 22-panel
    ("Q9UL51", "HCN2", "ahba_cortical"),
    ("Q9P1Z3", "HCN3", "ahba_cortical"),
    ("Q9Y3Q4", "HCN4", "ahba_cortical"),
    ("Q15858", "SCN9A", "magma_p"),
    ("Q14524", "SCN5A", "magma_p"),
    # BDNF/TrkB + growth factor
    ("P23560", "BDNF", "l2g_davies2018"),    # canonical
    ("Q16620", "NTRK2", "ahba_cortical"),    # duplicate; in 22-panel
    ("P28907", "CD38", "lit_otar"),
    # Sigma + opioid
    ("Q99720", "SIGMAR1", "lit_otar"),       # duplicate; in 22-panel
    ("P35372", "OPRM1", "ahba_cortical"),
    ("P41143", "OPRD1", "ahba_cortical"),
    ("P41145", "OPRK1", "ahba_cortical"),
    # Phosphodiesterases (cognition + AD)
    ("Q08499", "PDE4D", "magma_p"),          # duplicate
    ("O76083", "PDE9A", "lit_otar"),         # duplicate
    ("Q9UNP8", "PDE10A", "ahba_cortical"),
    ("Q07343", "PDE4B", "magma_p"),
    ("Q08493", "PDE4C", "magma_p"),
    ("Q14123", "PDE11A", "ahba_cortical"),
    # Histamine + orexin
    ("Q9Y5N1", "HRH3", "ahba_cortical"),     # duplicate
    ("P35367", "HRH1", "ahba_cortical"),
    ("P25021", "HRH2", "ahba_cortical"),
    ("O43614", "HCRTR2", "ahba_cortical"),   # duplicate
    ("O43613", "HCRTR1", "ahba_cortical"),   # duplicate
    # Adenosine + Sleep
    ("P30542", "ADORA1", "lit_otar"),
    ("P29274", "ADORA2A", "ahba_cortical"),
    ("P0DMS8", "ADORA3", "lit_otar"),
    # GSK3β + lithium pathway (per V8 worked example I.4)
    ("P49841", "GSK3B", "magma_p"),
    ("P49840", "GSK3A", "lit_otar"),
    ("Q9NSE2", "CISH", "lit_otar"),
    # Melatonin
    ("P48039", "MTNR1A", "ahba_cortical"),
    ("P49286", "MTNR1B", "ahba_cortical"),
    # Estrogen receptors
    ("P03372", "ESR1", "lit_otar"),
    ("Q92731", "ESR2", "lit_otar"),
    # Vitamin D
    ("P11473", "VDR", "lit_otar"),
    # ApoE pathway
    ("P02649", "APOE", "l2g_davies2018"),    # duplicate
    ("P15692", "VEGFA", "ahba_cortical"),
    # Misc cognition-adjacent
    ("Q15369", "ELOC", "magma_p"),
    ("Q14642", "INPP5A", "magma_p"),
    ("Q15303", "ERBB4", "lit_otar"),
    ("P00533", "EGFR", "lit_otar"),
    ("Q15078", "CDK5R1", "ahba_cortical"),
    ("Q9UQ52", "CDC42BPB", "magma_p"),
    ("O15226", "NKRF", "magma_p"),
    ("P19623", "SRM", "magma_p"),
    ("Q86Y39", "NDUFA11", "sc_zscore"),
    ("P09382", "LGALS1", "sc_zscore"),
    ("P17813", "ENG", "sc_zscore"),
    ("P31947", "SFN", "sc_zscore"),
    ("Q16566", "CAMK4", "ahba_cortical"),
    ("P14210", "HGF", "lit_otar"),
    ("Q14653", "IRF3", "magma_p"),
    ("P31350", "RRM2", "lit_otar"),
    # Astrocyte/oligo markers (Cluster D cell-type enrichment)
    ("P14136", "GFAP", "sc_zscore"),
    ("P13591", "NCAM1", "sc_zscore"),
    ("P22897", "MRC1", "sc_zscore"),
    ("P16284", "PECAM1", "sc_zscore"),
    ("P07900", "HSP90AA1", "sc_zscore"),
    # Trace amine + sigma extensions
    ("Q96RJ0", "TAAR1", "lit_otar"),
    ("Q9Y657", "SPIN1", "magma_p"),
    ("Q9NTI5", "PDS5B", "magma_p"),
]


# Canonical UniProt -> gene-symbol map for the 22-target panel + liability/expansion extras,
# sourced from data/raw/targets_seed.csv. Used so a panel anchor's gene_symbol is its REAL gene
# (e.g. GRIN2B) rather than its UniProt accession: the downstream gene-keyed NUTS reference-anchor
# match silently dropped any target whose gene_symbol fell back to the accession (notably
# GRIN2B/Q13224 -> the designed-6-anchor model fit with only 4-5 active). See BUG_AUDIT_2026-06 B2.
_CANONICAL_UNIPROT_TO_GENE: dict[str, str] = {
    "P22303": "ACHE", "P36544": "CHRNA7", "P42261": "GRIA1", "P42262": "GRIA2",
    "P42263": "GRIA3", "P48058": "GRIA4", "Q12879": "GRIN2A", "Q13224": "GRIN2B",
    "P21728": "DRD1", "Q01959": "SLC6A3", "P08913": "ADRA2A", "P23975": "SLC6A2",
    "Q9Y5N1": "HRH3", "O43613": "HCRTR1", "O43614": "HCRTR2", "Q08499": "PDE4D",
    "O76083": "PDE9A", "Q16620": "NTRK2", "Q99720": "SIGMAR1", "O43526": "KCNQ2",
    "O43525": "KCNQ3", "O60741": "HCN1",
    "P08908": "HTR1A", "Q13639": "HTR4", "P48067": "SLC6A9", "Q14416": "GRM2",
    "Q14832": "GRM3", "P41594": "GRM5", "P11229": "CHRM1", "P08173": "CHRM4",
    "P50406": "HTR6",
}


def _resolve_gene(uniprot: str) -> str:
    """Canonical gene symbol for a UniProt accession: the curated panel map first, then the
    expansion-target list, then the accession itself only as a last resort."""
    if uniprot in _CANONICAL_UNIPROT_TO_GENE:
        return _CANONICAL_UNIPROT_TO_GENE[uniprot]
    return next((g for u, g, _ in COGNITION_EXPANSION_TARGETS if u == uniprot), uniprot)


def build_expanded_panel(
    include_22_subset: bool = True,
    include_liability_subset: bool = True,
    dedupe_by_uniprot: bool = True,
) -> pd.DataFrame:
    """Build the V6.B.5 ~210-target panel as a DataFrame with provenance.

    Returns DataFrame: uniprot, gene_symbol, rules_fired, in_v6b_panel_22,
    in_liability_panel_44, notes.
    """
    targets: list[PanelTarget] = []

    # Add the 22-target panel first
    if include_22_subset:
        for uniprot in PANEL_22_TARGETS:
            gene = _resolve_gene(uniprot)
            targets.append(PanelTarget(
                uniprot=uniprot,
                gene_symbol=gene,
                rules_fired=["v6b_panel_22_anchor"],
                in_v6b_panel_22=True,
                in_liability_panel_44=uniprot in LIABILITY_PANEL_22_EXTRA,
                notes="Original V6.B 22-target cognition panel anchor",
            ))

    # Add liability panel extras
    if include_liability_subset:
        for uniprot in LIABILITY_PANEL_22_EXTRA:
            if any(t.uniprot == uniprot for t in targets):
                continue    # dedupe
            gene = _resolve_gene(uniprot)
            targets.append(PanelTarget(
                uniprot=uniprot,
                gene_symbol=gene,
                rules_fired=["v5_liability_panel_44"],
                in_v6b_panel_22=False,
                in_liability_panel_44=True,
                notes="V5 §8.0b 44-target liability panel anchor",
            ))

    # Add cognition expansion targets
    for uniprot, gene, rule in COGNITION_EXPANSION_TARGETS:
        if dedupe_by_uniprot and any(t.uniprot == uniprot for t in targets):
            # Already added (in 22 or liability panel); just add the rule
            for t in targets:
                if t.uniprot == uniprot and rule not in t.rules_fired:
                    t.rules_fired.append(rule)
            continue
        targets.append(PanelTarget(
            uniprot=uniprot,
            gene_symbol=gene,
            rules_fired=[rule],
            in_v6b_panel_22=uniprot in PANEL_22_TARGETS,
            in_liability_panel_44=uniprot in LIABILITY_PANEL_22_EXTRA,
            notes=f"V6.B.5 expansion via {rule}",
        ))

    rows = [{
        "uniprot": t.uniprot,
        "gene_symbol": t.gene_symbol,
        "rules_fired": "|".join(t.rules_fired),
        "n_rules": len(t.rules_fired),
        "in_v6b_panel_22": t.in_v6b_panel_22,
        "in_liability_panel_44": t.in_liability_panel_44,
        "l2g_max": t.l2g_max,
        "magma_p": t.magma_p,
        "ahba_r": t.ahba_r,
        "sc_zscore": t.sc_zscore,
        "lit_otar": t.lit_otar,
        "notes": t.notes,
    } for t in targets]
    return pd.DataFrame(rows)


def validate_panel(df: pd.DataFrame) -> dict[str, object]:
    """Sanity-check that the expanded panel satisfies V6.B.5 §F constraints."""
    n_total = len(df)
    n_22 = int(df["in_v6b_panel_22"].sum())
    n_liability = int(df["in_liability_panel_44"].sum())
    rules_counts: dict[str, int] = {}
    for r_str in df["rules_fired"]:
        for r in r_str.split("|"):
            rules_counts[r] = rules_counts.get(r, 0) + 1

    # 22-target panel must be strict subset
    panel_22_ok = (n_22 == len(PANEL_22_TARGETS))
    # All UniProts unique
    n_unique = df["uniprot"].nunique()
    dedup_ok = (n_unique == n_total)

    return {
        "n_total": n_total,
        "n_in_v6b_panel_22": n_22,
        "n_in_liability_panel_44": n_liability,
        "n_unique_uniprot": n_unique,
        "rules_distribution": rules_counts,
        "v6b_panel_22_subset_ok": panel_22_ok,
        "dedup_ok": dedup_ok,
        "v6b_panel_22_expected": len(PANEL_22_TARGETS),
    }


def availability() -> dict[str, object]:
    """Probe V6.B.5 panel expansion availability."""
    df = build_expanded_panel()
    val = validate_panel(df)
    return {
        "available": True,
        "n_panel_targets_total": val["n_total"],
        "n_v6b_panel_22": val["n_in_v6b_panel_22"],
        "n_liability_panel": val["n_in_liability_panel_44"],
        "rules_supported": sorted(set(
            r for r_str in df["rules_fired"] for r in r_str.split("|")
        )),
        "stub_mode": True,
        "note": ("V6.B.5 Stage 1: hand-curated approximation of GWAS+AHBA+SC+"
                 "Lit-OTAR query. Real-mode requires V6.B.5 Stage 2-3 "
                 "(live OT Genetics L2G + cellxgene-census + Moodie 2024 "
                 "g-cortical map alignment + Lit-OTAR Kafkas 2024)."),
    }
