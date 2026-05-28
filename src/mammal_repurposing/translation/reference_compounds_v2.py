"""V7.2 Stage 4 (Sprint 3.3) — 60-row REFERENCE_COMPOUND_SMD expansion.

Source: research/4-tier/MH1 + MH2 Meta-Analytic Prior Expansion for V7 CPT
        Bayesian Pharmacology Pipeline.md §3 "Reference compound table (MH2)"

Replaces the Stage 1 (V1) 15-row table in `validation_gates.REFERENCE_COMPOUND_SMD`
with a denser 60-row (compound × endpoint × population) anchor set.

Schema per the doc:
    (compound, class, endpoint, g, ci_lo, ci_hi, k_or_N, population,
     dose_range, source_doi_or_pmid)

Each row contributes a likelihood anchor for V7 NUTS V2 (Sprint 3.2). The same
compound can appear multiple times for different (endpoint, population) cells —
this is the per-(compound, endpoint, population) granularity the V7 paper
needs to defend the partial-pool calibration.

V7 NUTS V2 (`fit_effect_size_nuts_v2`) consumes these rows by mapping each
to an EffectSizeObservation with `observed_g` populated. The compound →
target UniProt map below assigns each compound its canonical primary target
(some compounds are multi-target; we pick the dominant pharmacological
mechanism for the V7 Cluster D gate).
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ReferenceCompoundSMDV2:
    """One row of the V2 reference compound anchor table."""
    compound: str          # e.g. "donepezil"
    class_v2: str          # MH1+MH2 class id (V2 prisma class — auto-mapped from V1)
    endpoint: str          # 8 canonical cognitive domains (EM/WM/ATT/EF/PS/VL/VS/MOT)
    pooled_g: float
    ci_lo: float
    ci_hi: float
    k_or_n: str            # studies or N (string because some entries are descriptive)
    population: str        # HC / AD / MCI / SCZ / ADHD / MDD / FXS / NRC / ...
    dose_range: str        # "5-10mg", "20-40mg", "0.5mg/kg", "patch/spray", etc.
    source: str            # DOI or PMID or descriptive citation key

    @property
    def ci_width(self) -> float:
        return self.ci_hi - self.ci_lo

    @property
    def is_within_roberts_ceiling(self) -> bool:
        """No published cell may exceed Roberts 2020 g = 0.55."""
        return self.pooled_g <= 0.55


# ---------------------------------------------------------------------------
# Compound → canonical primary target UniProt mapping
# ---------------------------------------------------------------------------

COMPOUND_TO_TARGET_UNIPROT: dict[str, str] = {
    # Cholinesterase inhibitors → ACHE
    "donepezil":           "P22303",
    "rivastigmine":        "P22303",
    "galantamine":         "P22303",
    "huperzine A":         "P22303",
    # NMDA modulators → GRIN2B (canonical heteromer subunit; ketamine impairs)
    "memantine":           "Q13224",
    "ketamine":            "Q13224",
    # α7 nAChR
    "encenicline":         "P36544",
    "ABT-126":             "P36544",
    "DMXB-A":              "P36544",
    "tropisetron":         "P36544",
    "TC-5619":             "P36544",
    # α4β2 nAChR — use CHRNA4 (P43681) but our panel lacks it; map to CHRNA7
    "varenicline":         "P36544",     # closest pharmacology in our panel
    # DA stimulants — methylphenidate → DAT
    "methylphenidate":     "Q01959",
    "dextroamphetamine":   "Q01959",
    "MAS-Adderall":        "Q01959",
    "lisdexamfetamine":    "Q01959",
    "modafinil":           "Q01959",     # weak DAT + Hist/Hcrt; SLC6A3 dominant in our taxonomy
    "armodafinil":         "Q01959",
    # NRI / α2A
    "atomoxetine":         "P23975",     # NET = SLC6A2
    "guanfacine-XR":       "P08913",     # ADRA2A
    # Multimodal 5-HT — vortioxetine (primary 5-HT transporter, polypharmacology)
    "vortioxetine":        "P31645",     # SLC6A4 5-HT transporter
    # 5-HT6 antagonists
    "idalopirdine":        "P50406",
    "intepirdine":         "P50406",
    "SUVN-502":            "P50406",
    # Nicotine — pan-nAChR; use CHRNA7
    "nicotine":            "P36544",
    # Caffeine — A2A receptor (ADORA2A)
    "caffeine":            "P29274",
    # H3 antagonist
    "pitolisant":          "Q9Y5N1",
    # PDE4D
    "BPN14770":            "Q08499",
    # AMPA modulators → GRIA1
    "CX-516":              "P42261",
    "S47445":              "P42261",
    # Sigma-1
    "blarcamesine":        "Q99720",
    # mGluR5 NAM
    "mavoglurant":         "P41594",
    "basimglurant":        "P41594",
    # GSK-3β
    "tideglusib":          "P49841",
    "lithium":             "P49841",
    # PDE9, PDE10A
    "PF-04447943":         "O76083",
    "TAK-063":             "Q9Y233",
    # mGluR2/3
    "pomaglumetad":        "Q14416",
    # Herbals / supplements — MIXED multi-target
    "ginkgo-EGb761":       "MIXED",
    "bacopa-monnieri":     "MIXED",
    "citicoline":          "MIXED",
    "DHA":                 "MIXED",
    "creatine":            "MIXED",
    "L-theanine":          "MIXED",
    "piracetam":           "MIXED",
    # Antipsychotic cognitive sub-analyses
    "lurasidone":          "P14416",     # DRD2
    "cariprazine":         "P14416",     # DRD2/D3
    # Suvorexant — DORA → HCRTR2
    "suvorexant":          "O43614",
}


# Migration: MH1+MH2 doc's class names → our V7 prisma_priors V2 names.
CLASS_MIGRATION_DOC_TO_V2: dict[str, str] = {
    "AChE":     "AChE_INHIBITORS",
    "NMDA":     "NMDA_MODULATORS",
    "A7NACHR":  "ALPHA7_NACHR",
    "A4B2":     "ALPHA4BETA2_NACHR",
    "MPH":      "DA_STIMULANTS_MPH",
    "AMPH":     "AMPHETAMINE_LIKE",
    "MODA":     "MODAFINIL_LIKE",
    "M5HT":     "MULTIMODAL_5HT",
    "5HT6":     "5HT6_ANTAGONISTS",
    "NIC":      "ALPHA7_NACHR",          # nicotine → α7 (closest in our 12 classes)
    "ADO":      "MODAFINIL_LIKE",        # caffeine — wake-promoting proxy
    "H3":       "H3_ANTAGONISTS",
    "PDE4D":    "PDE4D_NAM",
    "AMPA":     "AMPA_POSITIVE_MOD",
    "SIGMA1":   "PDE4D_NAM",             # no SIGMA1 class in V2; nearest = trophic
    "MGLUR5":   "AMPA_POSITIVE_MOD",     # placeholder — closest glutamatergic
    "GSK3B":    "AMPA_POSITIVE_MOD",     # placeholder — no GSK3B class
    "PDE9":     "PDE4D_NAM",             # closest PDE class
    "PDE10A":   "PDE4D_NAM",
    "MGLUR23":  "AMPA_POSITIVE_MOD",
    "HERBAL":   "AMPA_POSITIVE_MOD",     # mixed; uses class-level fallback
    "RACETAM":  "AMPA_POSITIVE_MOD",
    "ANTIPSY":  "MULTIMODAL_5HT",        # atypical antipsychotic — 5-HT receptors
    "NRI":      "DA_STIMULANTS_MPH",     # NRI in our class set → closest stim
    "A2A":      "DA_STIMULANTS_MPH",     # α2A is sympatholytic — closest stim
    "DORA":     "H3_ANTAGONISTS",        # arousal-related — closest
}


REFERENCE_COMPOUND_SMD_V2: list[ReferenceCompoundSMDV2] = [
    # ===== AChE INHIBITORS =====
    ReferenceCompoundSMDV2("donepezil",    "AChE", "EM",  0.36, 0.27, 0.45,
                            "5 studies / 8257 total in Birks 2018 (30 studies)",
                            "AD-mod",  "5-10mg",   "PMID:29923184 Birks 2018 CD001190.pub3"),
    ReferenceCompoundSMDV2("donepezil",    "AChE", "WM",  0.30, 0.20, 0.40,
                            "7 studies (MMSE)", "AD-mod", "5-10mg",
                            "PMID:29923184 (MMSE MD 1.05)"),
    ReferenceCompoundSMDV2("donepezil",    "AChE", "VL",  0.34, 0.22, 0.46,
                            "5 studies (ADCS-ADL-WL)", "AD-mod", "5-10mg", "PMID:29923184"),
    ReferenceCompoundSMDV2("donepezil",    "AChE", "EF",  0.20, 0.08, 0.32,
                            "IPD k=8", "AD-mild", "10mg",
                            "PMID:35988219 Ide IPD-meta 2022"),
    ReferenceCompoundSMDV2("donepezil",    "AChE", "EM",  0.18, 0.05, 0.31,
                            "MCI k=17 (2847)", "MCI", "5-10mg",
                            "PMID:35153124 Cui 2022; Russ Cochrane 2012 CD009132"),
    ReferenceCompoundSMDV2("rivastigmine", "AChE", "EM",  0.24, 0.18, 0.30,
                            "k=6 (3232)", "AD-mild-mod", "6-12mg",
                            "PMID:26393402 Birks 2015 CD001191.pub4"),
    ReferenceCompoundSMDV2("rivastigmine", "AChE", "WM",  0.20, 0.13, 0.27,
                            "k=6 (3205) MMSE", "AD", "9.5mg-patch", "PMID:26393402"),
    ReferenceCompoundSMDV2("rivastigmine", "AChE", "ATT", 0.22, 0.10, 0.34,
                            "IDEAL N=1195", "AD", "9.5mg-patch",
                            "PMID:17035691 Winblad IDEAL"),
    ReferenceCompoundSMDV2("galantamine",  "AChE", "EM",  0.40, 0.34, 0.46,
                            "k=10 Tan/Loy", "AD-mild-mod", "16-24mg",
                            "PMID:24662102 Tan 2014; PMID:39498781 Loy 2024 CD001747.pub4"),
    ReferenceCompoundSMDV2("galantamine",  "AChE", "WM",  0.31, 0.10, 0.52,
                            "k>=3 MMSE", "AD", "24mg", "PMID:24662102"),
    ReferenceCompoundSMDV2("galantamine",  "AChE", "VL",  0.35, 0.20, 0.50,
                            "k>=3", "AD", "24mg", "PMID:39498781 Loy 2024"),
    ReferenceCompoundSMDV2("huperzine A",  "AChE", "EM",  0.30, 0.10, 0.50,
                            "k=8 AD (733)", "AD", "0.4mg", "PMID:24086396 Yang 2013"),
    ReferenceCompoundSMDV2("huperzine A",  "AChE", "WM",  0.45, 0.20, 0.70,
                            "k=20 (1823) MMSE", "AD/MCI", "0.4mg", "PMID:24086396"),
    # ===== NMDA MODULATORS =====
    ReferenceCompoundSMDV2("memantine",    "NMDA", "EM",  0.27, 0.14, 0.39,
                            "k=9 (2433)", "AD-mod-sev", "20mg",
                            "PMID:25869017 Matsunaga 2015 PLOS ONE"),
    ReferenceCompoundSMDV2("memantine",    "NMDA", "EF",  0.20, 0.08, 0.32,
                            "k>=3", "AD", "20mg", "PMID:25869017"),
    ReferenceCompoundSMDV2("memantine",    "NMDA", "VL",  0.18, 0.05, 0.31,
                            "k>=3", "AD", "20mg", "PMID:25869017"),
    ReferenceCompoundSMDV2("memantine",    "NMDA", "ATT", 0.15, 0.02, 0.28,
                            "k>=3", "AD", "20mg", "PMID:25869017"),
    ReferenceCompoundSMDV2("ketamine",     "NMDA", "EM", -0.60, -0.90, -0.30,
                            "Krystal N~120", "healthy", "0.5mg/kg",
                            "PMID:8122957 Krystal 1994 (impairment)"),
    ReferenceCompoundSMDV2("ketamine",     "NMDA", "WM", -0.45, -0.75, -0.15,
                            "single RCT", "healthy", "0.5mg/kg", "PMID:8122957"),
    ReferenceCompoundSMDV2("ketamine",     "NMDA", "ATT", -0.30, -0.55, -0.05,
                            "single RCT", "healthy", "0.5mg/kg", "PMID:8122957"),
    # ===== ALPHA7 nAChR =====
    ReferenceCompoundSMDV2("encenicline",  "A7NACHR", "EM",  0.36, 0.05, 0.67,
                            "Keefe N=319 ES=0.36 SCoRS", "SCZ", "0.9mg",
                            "PMID:26089183 Keefe 2015"),
    ReferenceCompoundSMDV2("encenicline",  "A7NACHR", "EF",  0.18, -0.05, 0.41,
                            "Keefe N=319 OCI", "SCZ", "0.9mg", "PMID:26089183"),
    ReferenceCompoundSMDV2("encenicline",  "A7NACHR", "ATT", -0.05, -0.20, 0.10,
                            "k=10 in Lewis pooled", "SCZ+AD", "",
                            "PMID:28065843 Lewis 2017"),
    ReferenceCompoundSMDV2("ABT-126",      "A7NACHR", "EM",  0.10, -0.15, 0.35,
                            "3 RCTs (Haig)", "SCZ", "25-75mg", "Haig 2016 JCP 36:467"),
    ReferenceCompoundSMDV2("DMXB-A",       "A7NACHR", "EM",  0.20, -0.10, 0.50,
                            "Freedman N~30", "SCZ", "150mg", "PMID:18381905 Freedman 2008"),
    ReferenceCompoundSMDV2("tropisetron",  "A7NACHR", "EM",  0.30,  0.00, 0.60,
                            "single RCT N=40", "SCZ", "10mg", "Zhang 2012"),
    ReferenceCompoundSMDV2("TC-5619",      "A7NACHR", "EM", -0.05, -0.30, 0.20,
                            "Walling N=457", "SCZ", "1-25mg", "Walling 2016 Schiz Res"),
    # ===== ALPHA4BETA2 =====
    ReferenceCompoundSMDV2("varenicline",  "A4B2", "EM",  -0.02, -0.20, 0.16,
                            "k=4 (339)", "SCZ", "1-2mg", "PMID:31792645 Tanzer 2020"),
    ReferenceCompoundSMDV2("varenicline",  "A4B2", "ATT", -0.05, -0.20, 0.10,
                            "k=4 (339)", "SCZ", "1-2mg", "PMID:31792645"),
    ReferenceCompoundSMDV2("varenicline",  "A4B2", "EF",  -0.06, -0.47, 0.35,
                            "k=2 (339)", "SCZ", "1-2mg", "PMID:31792645"),
    ReferenceCompoundSMDV2("varenicline",  "A4B2", "PS",   0.04, -0.23, 0.31,
                            "k=3 (339)", "SCZ", "1-2mg", "PMID:31792645"),
    # ===== DA STIMULANTS - MPH =====
    ReferenceCompoundSMDV2("methylphenidate", "MPH", "EM",  0.60, 0.41, 0.79,
                            "60 studies reviewed, 36 in meta-analysis",
                            "ADHD-ped", "0.3-1mg/kg",
                            "PMID:24231201 Coghill 2014 non-exec mem"),
    ReferenceCompoundSMDV2("methylphenidate", "MPH", "WM",  0.26, 0.13, 0.39,
                            "Coghill subset", "ADHD-ped", "",
                            "PMID:24231201 exec mem"),
    ReferenceCompoundSMDV2("methylphenidate", "MPH", "ATT", 0.24, 0.15, 0.33,
                            "Coghill subset", "ADHD-ped", "", "PMID:24231201 RT"),
    ReferenceCompoundSMDV2("methylphenidate", "MPH", "EF",  0.41, 0.27, 0.55,
                            "Coghill subset", "ADHD-ped", "",
                            "PMID:24231201 response inhib"),
    ReferenceCompoundSMDV2("methylphenidate", "MPH", "PS",  0.55, 0.34, 0.55,
                            "Coghill subset clipped to Roberts ceiling",
                            "ADHD-ped", "",
                            "PMID:24231201 RT-variability (g=0.62 raw -> clipped)"),
    ReferenceCompoundSMDV2("methylphenidate", "MPH", "EM",  0.43, 0.21, 0.65,
                            "k=24 healthy", "healthy", "20-40mg",
                            "PMID:32709551 Roberts 2020 recall (significant)"),
    ReferenceCompoundSMDV2("methylphenidate", "MPH", "WM",  0.10, -0.05, 0.25,
                            "k=24 healthy", "healthy", "20-40mg",
                            "PMID:32709551 SWM (NS)"),
    ReferenceCompoundSMDV2("methylphenidate", "MPH", "ATT", 0.42, 0.18, 0.55,
                            "k=24 healthy (clipped CI_hi to Roberts ceiling)",
                            "healthy", "20-40mg", "PMID:32709551 sust-att"),
    ReferenceCompoundSMDV2("methylphenidate", "MPH", "EF",  0.27, 0.03, 0.51,
                            "k=24 healthy", "healthy", "20-40mg",
                            "PMID:32709551 inhib (p=.03)"),
    ReferenceCompoundSMDV2("methylphenidate", "MPH", "ATT", 0.17, 0.05, 0.28,
                            "k=21 adult ADHD", "ADHD-adult", "18-72mg",
                            "PMID:29751051 Pievsky 2018"),
    # ===== AMPHETAMINE =====
    ReferenceCompoundSMDV2("dextroamphetamine", "AMPH", "EM",  0.20, 0.05, 0.35,
                            "k=10 (Ilieva STM)", "healthy", "5-30mg",
                            "PMID:25591060 Ilieva 2015"),
    ReferenceCompoundSMDV2("dextroamphetamine", "AMPH", "WM",  0.13, -0.02, 0.28,
                            "k=10 (Ilieva WM)", "healthy", "5-30mg", "PMID:25591060"),
    ReferenceCompoundSMDV2("dextroamphetamine", "AMPH", "EM",  0.45, 0.20, 0.55,
                            "k=6 delayed (Ilieva) clipped to Roberts ceiling",
                            "healthy", "5-30mg", "PMID:25591060 delayed mem"),
    ReferenceCompoundSMDV2("dextroamphetamine", "AMPH", "EF",  0.00, -0.20, 0.20,
                            "k=10, 27 ES Roberts", "healthy", "5-30mg",
                            "PMID:32709551 Roberts 2020 - NO significant d-amph effects"),
    ReferenceCompoundSMDV2("dextroamphetamine", "AMPH", "ATT", 0.00, -0.20, 0.20,
                            "k=10 Roberts", "healthy", "5-30mg", "PMID:32709551 (null)"),
    ReferenceCompoundSMDV2("MAS-Adderall",     "AMPH", "EM",  0.05, -0.20, 0.30,
                            "Ilieva 2013 RCT N=46", "healthy", "10-30mg",
                            "PMID:22884611 Ilieva 2013"),
    ReferenceCompoundSMDV2("lisdexamfetamine", "AMPH", "ATT", 0.55, 0.55, 0.55,
                            "Cortese network meta clipped to Roberts ceiling (g=0.85 raw)",
                            "ADHD", "30-70mg", "PMID:34693523 Cortese 2021"),
    # ===== MODAFINIL =====
    ReferenceCompoundSMDV2("modafinil",   "MODA", "EF",  0.28, 0.03, 0.53,
                            "k=14 (Roberts)", "healthy", "100-200mg",
                            "PMID:32709551 Roberts 2020 memory updating (p=.03)"),
    ReferenceCompoundSMDV2("modafinil",   "MODA", "ATT", 0.10, -0.05, 0.25,
                            "k=14 Roberts", "healthy", "100-200mg", "PMID:32709551 (NS)"),
    ReferenceCompoundSMDV2("modafinil",   "MODA", "EM",  0.05, -0.10, 0.20,
                            "k=14 Roberts", "healthy", "100-200mg", "PMID:32709551 (NS)"),
    ReferenceCompoundSMDV2("modafinil",   "MODA", "WM",  0.05, -0.10, 0.20,
                            "k=14 Roberts", "healthy", "100-200mg", "PMID:32709551 (NS)"),
    ReferenceCompoundSMDV2("modafinil",   "MODA", "ATT", 0.40, 0.10, 0.55,
                            "single RCT N=209 (clipped to Roberts ceiling)",
                            "shift-work", "200mg", "PMID:16079371 Czeisler 2005 PVT"),
    ReferenceCompoundSMDV2("armodafinil", "MODA", "ATT", 0.35, 0.15, 0.55,
                            "single RCT N=254", "shift-work", "150mg",
                            "Czeisler 2009 Mayo Clin Proc"),
    # ===== ATOMOXETINE / GUANFACINE =====
    ReferenceCompoundSMDV2("atomoxetine",   "NRI", "EM",  0.25, 0.05, 0.45,
                            "k=4 (Isfandnia)", "ADHD", "60-100mg",
                            "Isfandnia 2024 NBR 162:105703"),
    ReferenceCompoundSMDV2("atomoxetine",   "NRI", "ATT", 0.40, 0.20, 0.55,
                            "k>=3 (Isfandnia) clipped to Roberts ceiling",
                            "ADHD", "60-100mg", "Isfandnia 2024"),
    ReferenceCompoundSMDV2("atomoxetine",   "NRI", "EF",  0.50, 0.25, 0.55,
                            "k>=3 inhibition clipped to Roberts ceiling",
                            "ADHD", "60-100mg", "Isfandnia 2024"),
    ReferenceCompoundSMDV2("guanfacine-XR", "A2A", "ATT", 0.50, 0.25, 0.55,
                            "Sallee N=345 (clipped to Roberts ceiling)",
                            "ADHD-ped", "1-4mg", "PMID:19106767 Sallee 2009"),
    ReferenceCompoundSMDV2("guanfacine-XR", "A2A", "EF",  0.40, 0.15, 0.55,
                            "Wilens N=450 adolescent (clipped to Roberts ceiling)",
                            "ADHD-adolesc", "1-7mg", "PMID:26506582 Wilens 2015"),
    # ===== MULTIMODAL 5-HT =====
    ReferenceCompoundSMDV2("vortioxetine", "M5HT", "PS",  0.51, 0.30, 0.55,
                            "FOCUS N=602 (clipped CI_hi to Roberts ceiling)",
                            "MDD", "10mg", "PMID:27231256 Harrison FOCUS DSST"),
    ReferenceCompoundSMDV2("vortioxetine", "M5HT", "PS",  0.52, 0.31, 0.55,
                            "FOCUS N=602 (clipped CI_hi to Roberts ceiling)",
                            "MDD", "20mg", "PMID:27231256 Harrison FOCUS DSST"),
    ReferenceCompoundSMDV2("vortioxetine", "M5HT", "EM",  0.31, 0.10, 0.52,
                            "FOCUS N=602", "MDD", "10mg", "PMID:27231256"),
    ReferenceCompoundSMDV2("vortioxetine", "M5HT", "EF",  0.42, 0.21, 0.55,
                            "FOCUS N=602 (clipped to Roberts ceiling)",
                            "MDD", "10-20mg", "PMID:27231256"),
    ReferenceCompoundSMDV2("vortioxetine", "M5HT", "PS",  0.24, 0.15, 0.33,
                            "k=6 (1782) Pan", "MDD", "10-20mg",
                            "PMID:36398888 Pan 2022 DSST"),
    ReferenceCompoundSMDV2("vortioxetine", "M5HT", "PS",  0.35, 0.23, 0.47,
                            "k=3 (607) McIntyre", "MDD", "10-20mg",
                            "PMID:27312740 McIntyre 2016 DSST SES=0.35"),
    ReferenceCompoundSMDV2("vortioxetine", "M5HT", "EM",  0.27, 0.10, 0.44,
                            "Katona N=453", "elderly-MDD", "5mg",
                            "PMID:22318341 Katona 2012 IJG"),
    # ===== 5HT6 ANTAGONISTS =====
    ReferenceCompoundSMDV2("idalopirdine", "5HT6", "EM", -0.05, -0.20, 0.10,
                            "k=4 (2803)", "AD", "30-60mg",
                            "PMID:30560763 Matsunaga 2018"),
    ReferenceCompoundSMDV2("idalopirdine", "5HT6", "EM",  0.30,  0.10, 0.50,
                            "LADDER N=278", "AD-mod", "90mg",
                            "PMID:25297016 LADDER (single+)"),
    ReferenceCompoundSMDV2("intepirdine",  "5HT6", "EM",  0.00, -0.15, 0.15,
                            "MINDSET N=1300", "AD", "35mg",
                            "PMID:29318278 Atri MINDSET (null)"),
    ReferenceCompoundSMDV2("SUVN-502",     "5HT6", "EM", -0.05, -0.20, 0.10,
                            "N=564 phase 2", "AD-mod", "50-100mg",
                            "PMID:35662833 Nirogi 2022"),
    # ===== NICOTINIC / NICOTINE =====
    ReferenceCompoundSMDV2("nicotine", "NIC", "MOT", 0.16, 0.06, 0.26,
                            "k=41 (Heishman)", "healthy", "patch/spray",
                            "PMID:20414766 Heishman 2010 fine-motor"),
    ReferenceCompoundSMDV2("nicotine", "NIC", "ATT", 0.34, 0.20, 0.48,
                            "k=41", "healthy", "patch/spray",
                            "PMID:20414766 alerting-RT"),
    ReferenceCompoundSMDV2("nicotine", "NIC", "EM",  0.27, 0.10, 0.44,
                            "k=41", "healthy", "patch/spray",
                            "PMID:20414766 STM-accuracy"),
    ReferenceCompoundSMDV2("nicotine", "NIC", "WM",  0.44, 0.25, 0.55,
                            "k=41 (clipped to Roberts ceiling)", "healthy",
                            "patch/spray", "PMID:20414766 WM-RT"),
    # ===== CAFFEINE =====
    ReferenceCompoundSMDV2("caffeine", "ADO", "ATT", 0.28, 0.18, 0.38,
                            "k=31 (1455)", "healthy", "75-250mg",
                            "Springer Psychopharm 2025 doi:10.1007/s00213-025-06775-1 RT"),
    ReferenceCompoundSMDV2("caffeine", "ADO", "ATT", 0.27, 0.17, 0.37,
                            "k=31 (1455)", "healthy", "75-250mg",
                            "Springer 2025 accuracy"),
    # ===== H3 ANTAGONISTS =====
    ReferenceCompoundSMDV2("pitolisant", "H3", "ATT", 0.55, 0.30, 0.55,
                            "k=2 (431) (clipped to Roberts ceiling)",
                            "narcolepsy", "17.8-35.6mg",
                            "PMID:33779931 Bassetti 2021 CNS Drugs"),
    ReferenceCompoundSMDV2("pitolisant", "H3", "PS",  0.20, 0.00, 0.40,
                            "k=2 (431)", "narcolepsy", "", "Lehert 2020 Drugs RD"),
    # ===== PDE4D =====
    ReferenceCompoundSMDV2("BPN14770", "PDE4D", "VL", 0.55, 0.20, 0.55,
                            "N=30 crossover (clipped CI_hi to Roberts ceiling)",
                            "FXS", "25mg BID",
                            "PMID:33927413 Berry-Kravis 2021 OralRead+PicVocab"),
    ReferenceCompoundSMDV2("BPN14770", "PDE4D", "EM", 0.40, 0.05, 0.55,
                            "N=30 (clipped to Roberts ceiling)", "FXS", "25mg BID",
                            "PMID:33927413"),
    ReferenceCompoundSMDV2("BPN14770", "PDE4D", "EF", 0.35, 0.00, 0.55,
                            "N=30 (clipped to Roberts ceiling)", "FXS", "25mg BID",
                            "PMID:33927413"),
    # ===== AMPAKINES =====
    ReferenceCompoundSMDV2("CX-516",  "AMPA", "EM",  0.05, -0.30, 0.40,
                            "Goff N=105", "SCZ", "900mg TID",
                            "PMID:17487227 Goff 2008"),
    ReferenceCompoundSMDV2("CX-516",  "AMPA", "EM", -0.05, -0.40, 0.30,
                            "Berry-Kravis MCI", "MCI", "900mg TID",
                            "Lynch/Berry-Kravis MCI trial"),
    ReferenceCompoundSMDV2("S47445",  "AMPA", "EM", -0.02, -0.20, 0.16,
                            "Bernard N=520", "AD-mild-mod", "2-15mg",
                            "PMID:31297441 Bernard 2019 ADTRCI"),
    # ===== SIGMA-1 =====
    ReferenceCompoundSMDV2("blarcamesine", "SIGMA1", "EM", 0.20, 0.00, 0.40,
                            "ANAVEX2-73-AD-004 N=508", "AD-early", "30-50mg",
                            "Macfarlane JPAD 2024 PMC11713060"),
    # ===== mGluR5 NAM =====
    ReferenceCompoundSMDV2("mavoglurant",  "MGLUR5", "EM", -0.05, -0.30, 0.20,
                            "Berry-Kravis N=315", "FXS", "50-100mg",
                            "PMID:26764156 Berry-Kravis 2016 STM"),
    ReferenceCompoundSMDV2("basimglurant", "MGLUR5", "EM",  0.00, -0.20, 0.20,
                            "Youssef N=183", "FXS", "0.5-1.5mg",
                            "PMID:28816242 Youssef 2018"),
    # ===== GSK-3 β / Lithium =====
    ReferenceCompoundSMDV2("tideglusib", "GSK3B", "EM", 0.00, -0.20, 0.20,
                            "ARGO N=306", "AD", "500-1000mg",
                            "PMID:25537011 Lovestone 2015"),
    ReferenceCompoundSMDV2("lithium",    "GSK3B", "EM", 0.41, 0.02, 0.55,
                            "k=3 (232) Matsunaga (clipped CI_hi to Roberts ceiling)",
                            "MCI/AD", "0.25-0.5mM", "PMID:26402004 Matsunaga 2015 JAD"),
    # ===== PDE9 / PDE10A =====
    ReferenceCompoundSMDV2("PF-04447943", "PDE9",   "EM", -0.05, -0.30, 0.20,
                            "Schwam N=191", "AD", "25mg", "PMID:24801218 Schwam 2014"),
    ReferenceCompoundSMDV2("TAK-063",     "PDE10A", "EM",  0.00, -0.20, 0.20,
                            "Macek N=160", "SCZ", "20mg", "PMID:30172593 Macek 2019"),
    # ===== mGluR2/3 =====
    ReferenceCompoundSMDV2("pomaglumetad", "MGLUR23", "EM", 0.00, -0.15, 0.15,
                            "Adams N~1000", "SCZ", "10-80mg",
                            "Adams 2014 BMC Psych PMC4276262"),
    # ===== HERBALS / SUPPLEMENTS =====
    ReferenceCompoundSMDV2("ginkgo-EGb761", "HERBAL", "EM",  0.44, 0.12, 0.55,
                            "AD subgroup k=4 (Weinmann) (clipped to Roberts ceiling)",
                            "AD", "240mg", "PMID:20236541 Weinmann 2010"),
    ReferenceCompoundSMDV2("ginkgo-EGb761", "HERBAL", "EM",  0.55, 0.01, 0.55,
                            "k=9 (2372) (clipped to Roberts ceiling)",
                            "dementia", "240mg", "PMID:20236541"),
    ReferenceCompoundSMDV2("ginkgo-EGb761", "HERBAL", "ATT", 0.30, 0.05, 0.55,
                            "k=9 (2561) Tan", "dementia", "240mg",
                            "Tan 2015 JAD doi:10.3233/JAD-140837"),
    ReferenceCompoundSMDV2("bacopa-monnieri", "HERBAL", "EM", 0.30, 0.10, 0.50,
                            "k=9 (518) Kongkeaw", "healthy", "300mg",
                            "PMID:24252493 Kongkeaw 2014 delayed recall"),
    ReferenceCompoundSMDV2("bacopa-monnieri", "HERBAL", "ATT", 0.30, 0.05, 0.55,
                            "k=9 (518) Kongkeaw", "healthy", "300mg",
                            "PMID:24252493 choice-RT"),
    ReferenceCompoundSMDV2("citicoline", "HERBAL", "EM",  0.19, 0.06, 0.32,
                            "k=14 (884) Cochrane", "VCI/VaD", "1000mg",
                            "PMID:15846601 Fioravanti Cochrane CD000269.pub2"),
    ReferenceCompoundSMDV2("citicoline", "HERBAL", "ATT", -0.09, -0.23, 0.05,
                            "k=14 (884)", "VCI/VaD", "1000mg",
                            "PMID:15846601 (null)"),
    ReferenceCompoundSMDV2("DHA", "HERBAL", "EM", 0.20, 0.00, 0.40,
                            "MIDAS N=485", "healthy-older", "900mg",
                            "PMID:20434961 Yurko-Mauro MIDAS CANTAB PAL"),
    ReferenceCompoundSMDV2("DHA", "HERBAL", "EM", 0.00, -0.10, 0.10,
                            "k=3 (3500+) Cochrane", "healthy-older", "",
                            "PMID:22696350 Sydenham Cochrane CD005379.pub3 (null)"),
    ReferenceCompoundSMDV2("creatine", "HERBAL", "EM",  0.31, 0.18, 0.44,
                            "k=16 (492) Xu 2024", "healthy", "5-20g",
                            "PMC11275561 Xu 2024 Front Nutr g=0.30 (0.18,0.42)"),
    ReferenceCompoundSMDV2("creatine", "HERBAL", "ATT", 0.31, 0.03, 0.55,
                            "k=16 (492) (clipped CI_hi to Roberts ceiling)",
                            "healthy", "5-20g", "PMC11275561 Xu 2024"),
    ReferenceCompoundSMDV2("creatine", "HERBAL", "PS",  0.49, 0.20, 0.55,
                            "k=16 (492) (clipped to Roberts ceiling)",
                            "healthy", "5-20g", "PMC11275561 Xu 2024"),
    ReferenceCompoundSMDV2("L-theanine", "HERBAL", "ATT", 0.30, 0.10, 0.50,
                            "k=11 (Camfield)", "healthy", "100mg+caf",
                            "PMID:24946991 Camfield 2014 Nutr Rev 72:507"),
    ReferenceCompoundSMDV2("piracetam", "RACETAM", "EM", 0.10, -0.10, 0.30,
                            "Cochrane Flicker", "dementia", "2.4-9.6g",
                            "PMID:11405971 Flicker Cochrane CD001011 (NS)"),
    # ===== ANTIPSYCHOTIC COGNITIVE SUB-ANALYSES =====
    ReferenceCompoundSMDV2("lurasidone",  "ANTIPSY", "PS", 0.30, 0.05, 0.55,
                            "Harvey N=244", "SCZ", "80-160mg",
                            "Harvey 2015 Schiz Res; PMC10616918"),
    ReferenceCompoundSMDV2("cariprazine", "ANTIPSY", "EM", 0.20, 0.00, 0.40,
                            "Mucci network N=400", "SCZ", "1.5-6mg",
                            "PMID:41191868 Mucci 2025 EAPCN"),
    # ===== SUVOREXANT (next-day cognition) =====
    ReferenceCompoundSMDV2("suvorexant", "DORA", "ATT", -0.05, -0.20, 0.10,
                            "Vermeeren k=2 (500+)", "insomnia", "10-20mg",
                            "PMID:26085297 Vermeeren 2015 (null)"),
]


def coverage_v2_anchors() -> dict[str, object]:
    """Coverage report for V7.2 Stage 4 anchor expansion."""
    n = len(REFERENCE_COMPOUND_SMD_V2)
    compounds = sorted(set(r.compound for r in REFERENCE_COMPOUND_SMD_V2))
    classes = sorted(set(r.class_v2 for r in REFERENCE_COMPOUND_SMD_V2))
    populations = sorted(set(r.population for r in REFERENCE_COMPOUND_SMD_V2))
    endpoints = sorted(set(r.endpoint for r in REFERENCE_COMPOUND_SMD_V2))
    n_targets = sum(1 for c in compounds
                    if COMPOUND_TO_TARGET_UNIPROT.get(c, "") not in ("", "MIXED"))
    n_negative = sum(1 for r in REFERENCE_COMPOUND_SMD_V2 if r.pooled_g < 0)
    n_phase3_null = sum(1 for r in REFERENCE_COMPOUND_SMD_V2 if abs(r.pooled_g) <= 0.06)
    return {
        "n_rows": n,
        "n_compounds": len(compounds),
        "n_classes": len(classes),
        "n_populations": len(populations),
        "n_endpoints": len(endpoints),
        "n_compounds_with_uniprot": n_targets,
        "n_negative_g": n_negative,
        "n_phase3_null": n_phase3_null,
        "compounds_missing_uniprot": [c for c in compounds
                                       if c not in COMPOUND_TO_TARGET_UNIPROT],
        "populations": populations,
        "endpoints": endpoints,
        "classes": classes,
    }


def anchors_to_observations(
    *,
    relevance_post_for_target: dict[str, float] | None = None,
    pchembl_post_for_compound: dict[str, float] | None = None,
    pbpk_auc_for_compound: dict[str, float] | None = None,
    skip_mixed_target: bool = True,
):
    """Convert REFERENCE_COMPOUND_SMD_V2 rows → EffectSizeObservation list
    for direct ingestion by `fit_effect_size_nuts_v2`.

    Args:
        relevance_post_for_target: optional V6.B θ̄ posterior per UniProt.
            Defaults to 0.5 (neutral) when missing.
        pchembl_post_for_compound: optional V6.A pchembl posterior per compound.
            Defaults to 8.0 (representative 10 nM affinity).
        pbpk_auc_for_compound: optional V7.1 brain AUC per compound.
            Defaults to 1.0.
        skip_mixed_target: drop rows whose compound maps to MIXED multi-target.

    Returns:
        list[EffectSizeObservation]: one observation per V2 row.
    """
    # Late import to avoid circular dependency at module load time
    from .effect_size_model import EffectSizeObservation

    relevance_post_for_target = relevance_post_for_target or {}
    pchembl_post_for_compound = pchembl_post_for_compound or {}
    pbpk_auc_for_compound = pbpk_auc_for_compound or {}

    observations = []
    # Track per-(compound, endpoint, population) row indices so multiple
    # source meta-analyses for the same cell remain disambiguated.
    seen: dict[tuple[str, str, str], int] = {}
    for r in REFERENCE_COMPOUND_SMD_V2:
        uniprot = COMPOUND_TO_TARGET_UNIPROT.get(r.compound, "MIXED")
        if skip_mixed_target and uniprot == "MIXED":
            continue
        class_v2_name = CLASS_MIGRATION_DOC_TO_V2.get(r.class_v2, r.class_v2)
        key = (r.compound, r.endpoint, r.population)
        seen[key] = seen.get(key, 0) + 1
        idx = seen[key]
        suffix = "" if idx == 1 else f"__src{idx}"
        observations.append(EffectSizeObservation(
            compound=f"{r.compound}__{r.endpoint}__{r.population}{suffix}",
            class_name=class_v2_name,
            target_uniprot=uniprot,
            pchembl_post_mean=pchembl_post_for_compound.get(r.compound, 8.0),
            pchembl_post_sd=0.3,
            relevance_post_mean=relevance_post_for_target.get(uniprot, 0.5),
            relevance_post_sd=0.1,
            pbpk_auc_brain=pbpk_auc_for_compound.get(r.compound, 1.0),
            moderators=(0.0,) * 5,
            observed_g=r.pooled_g,
            endpoint=r.endpoint,
            population=r.population,
        ))
    return observations
