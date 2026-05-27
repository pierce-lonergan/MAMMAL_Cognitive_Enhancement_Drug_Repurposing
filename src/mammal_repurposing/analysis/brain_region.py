"""§8.6 — Brain-region selectivity preview (V6 Cluster D entry point).

Maps the 22-target cognition panel to coarse anatomical biases via a
hand-curated lookup synthesised from:
  - Allen Human Brain Atlas (AHBA, Hawrylycz 2012 Nature 489:391) regional
    expression patterns
  - GTEx brain subregion bulk RNA-seq
  - Siletti et al. 2023 Science single-cell brain atlas
  - Literature review of receptor binding-site distribution

Per V6 Cluster D plan (research/4-tier/Multi-Source Neurobiological Prior...md),
the full Bayesian model is ~16 weeks. This is a STATIC FALLBACK so the V5
wet-lab shortlist can still surface "this compound's top target sits in the
hippocampus" annotations without the full PyMC pipeline.

The categories follow Hansen 2022 (Nat Neurosci 25:1569) cognitive-gradient
PC1:
    cortex-biased    — cortical pyramidal / interneuron enriched
    subcortical      — striatum / thalamus / hypothalamus
    brainstem        — LC, raphe, VTA
    hippocampal      — CA1/CA3/DG/subiculum
    mixed            — multi-region expression
"""

from __future__ import annotations

import pandas as pd

# Per-target brain-region bias — hand-curated, V6-Cluster-D-style preview.
# Source notes per target:
#   ACHE         — broadly expressed; mixed (cortex + striatum + brainstem)
#   CHRNA7       — cortical L2/3 + hippocampal CA1 + basal forebrain (mixed)
#   GRIA1-4      — cortical glutamate; pyramidal AMPA (cortex-biased)
#   GRIN2A/2B    — cortical pyramidal + hippocampus (cortex+hippocampal)
#   DRD1         — striatal medium spiny neurons (subcortical)
#   DRD2         — striatal MSNs + pituitary (subcortical)
#   SLC6A3 (DAT) — VTA / substantia nigra dopamine neurons (brainstem)
#   SLC6A2 (NET) — LC noradrenergic neurons (brainstem)
#   ADRA2A       — LC + prefrontal cortex (mixed)
#   HRH3         — broad cortical + striatal H3 autoreceptors (mixed)
#   HCRTR1/2     — lateral hypothalamus + downstream cortical projections (mixed)
#   PDE4D        — broadly expressed; cortex + hippocampus (cortex+hippocampal)
#   PDE9A        — striatum + cortex (subcortical+cortex)
#   NTRK2 (TrkB) — cortical L5 + hippocampal CA3 + cerebellum (mixed)
#   SIGMAR1      — broad ER chaperone (mixed)
#   KCNQ2/3      — cortical L5/6 axon-initial-segment (cortex-biased)
#   HCN1         — hippocampal CA1 distal dendrites (hippocampal)
BRAIN_REGION_BIAS: dict[str, dict] = {
    "P22303":  {"gene": "ACHE",    "bias": "mixed",
                "primary_region": "cortex+striatum+brainstem",
                "note": "Broadly expressed cholinesterase"},
    "P36544":  {"gene": "CHRNA7",  "bias": "mixed",
                "primary_region": "cortical L2/3 + hippocampal CA1 + basal forebrain",
                "note": "Nicotinic α7; CA1 enrichment per Siletti 2023"},
    "P42261":  {"gene": "GRIA1",   "bias": "cortex-biased",
                "primary_region": "cortical pyramidal",
                "note": "AMPA GluA1; classic cortex marker"},
    "P42262":  {"gene": "GRIA2",   "bias": "cortex-biased",
                "primary_region": "cortical pyramidal + cerebellum",
                "note": "AMPA GluA2; broadest AMPA subunit"},
    "P42263":  {"gene": "GRIA3",   "bias": "cortex-biased",
                "primary_region": "cortical pyramidal",
                "note": "AMPA GluA3"},
    "P48058":  {"gene": "GRIA4",   "bias": "cortex-biased",
                "primary_region": "cerebellar Purkinje + cortical interneurons",
                "note": "AMPA GluA4; cerebellum-enriched"},
    "Q12879":  {"gene": "GRIN2A",  "bias": "cortex+hippocampal",
                "primary_region": "cortical pyramidal + hippocampus",
                "note": "NMDA NR2A; mature postsynaptic"},
    "Q13224":  {"gene": "GRIN2B",  "bias": "cortex+hippocampal",
                "primary_region": "cortical pyramidal + hippocampus + striatum",
                "note": "NMDA NR2B; ifenprodil-class antagonists at ATD interface"},
    "P21728":  {"gene": "DRD1",    "bias": "subcortical",
                "primary_region": "striatum (direct pathway MSNs) + cortical L5",
                "note": "D1 dopamine receptor"},
    "Q01959":  {"gene": "SLC6A3",  "bias": "brainstem",
                "primary_region": "VTA / substantia nigra dopaminergic neurons",
                "note": "DAT; presynaptic on midbrain DA terminals"},
    "P08913":  {"gene": "ADRA2A",  "bias": "mixed",
                "primary_region": "locus coeruleus + prefrontal cortex",
                "note": "α2A adrenergic; LC autoreceptor + PFC"},
    "P23975":  {"gene": "SLC6A2",  "bias": "brainstem",
                "primary_region": "locus coeruleus noradrenergic neurons",
                "note": "NET; presynaptic on LC terminals"},
    "Q9Y5N1":  {"gene": "HRH3",    "bias": "mixed",
                "primary_region": "cortical + striatal H3 autoreceptors",
                "note": "Histamine H3"},
    "O43613":  {"gene": "HCRTR1",  "bias": "mixed",
                "primary_region": "lateral hypothalamus + LC + cortical projections",
                "note": "Orexin receptor 1"},
    "O43614":  {"gene": "HCRTR2",  "bias": "mixed",
                "primary_region": "lateral hypothalamus + tuberomammillary",
                "note": "Orexin receptor 2"},
    "Q08499":  {"gene": "PDE4D",   "bias": "cortex+hippocampal",
                "primary_region": "cortical pyramidal + hippocampus",
                "note": "cAMP PDE; broadly expressed"},
    "O76083":  {"gene": "PDE9A",   "bias": "subcortical+cortex",
                "primary_region": "striatum + cortical interneurons",
                "note": "cGMP PDE"},
    "Q16620":  {"gene": "NTRK2",   "bias": "mixed",
                "primary_region": "cortical L5 + hippocampal CA3 + cerebellum",
                "note": "TrkB BDNF receptor"},
    "Q99720":  {"gene": "SIGMAR1", "bias": "mixed",
                "primary_region": "ER chaperone — broad",
                "note": "Sigma-1 receptor; ER chaperone-receptor"},
    "O43526":  {"gene": "KCNQ2",   "bias": "cortex-biased",
                "primary_region": "cortical L5/6 axon-initial-segment",
                "note": "Kv7.2 M-current"},
    "O43525":  {"gene": "KCNQ3",   "bias": "cortex-biased",
                "primary_region": "cortical L5/6 axon-initial-segment",
                "note": "Kv7.3 M-current; KCNQ2 partner"},
    "O60741":  {"gene": "HCN1",    "bias": "hippocampal",
                "primary_region": "hippocampal CA1 distal dendrites",
                "note": "Hyperpolarization-activated cation channel"},
}


def annotate(
    df: pd.DataFrame,
    target_col: str = "mammal_best_target",
) -> pd.DataFrame:
    """Add `brain_bias`, `brain_primary_region`, `brain_note` columns to df.

    Looks up target via the uniprot column (preferred) or gene_symbol column.
    """
    out = df.copy()
    biases: list[str] = []
    regions: list[str] = []
    notes: list[str] = []
    for _, r in out.iterrows():
        t = r.get(target_col)
        if t and t in BRAIN_REGION_BIAS:
            e = BRAIN_REGION_BIAS[t]
            biases.append(e["bias"])
            regions.append(e["primary_region"])
            notes.append(e["note"])
        else:
            biases.append("")
            regions.append("")
            notes.append("")
    out["brain_bias"] = biases
    out["brain_primary_region"] = regions
    out["brain_note"] = notes
    return out


def summary_counts() -> dict[str, int]:
    """Tally of bias categories across the 22-target panel."""
    from collections import Counter
    return dict(Counter(e["bias"] for e in BRAIN_REGION_BIAS.values()))
