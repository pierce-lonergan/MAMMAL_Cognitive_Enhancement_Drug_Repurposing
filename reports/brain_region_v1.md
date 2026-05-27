# Brain-Region Selectivity v1 (§8.6 preview)

Hand-curated mapping of each 22-target cognition panel member to a coarse anatomical bias category, synthesised from AHBA (Hawrylycz 2012), GTEx, and Siletti 2023 single-cell atlas. **This is the V6 Cluster D entry point** — the full Bayesian AHBA+OT-Genetics+cellxgene model is documented at `research/4-tier/Multi-Source Neurobiological Prior...md` (~16 weeks).

## Per-target panel bias distribution (22 targets)

| Bias category | N targets |
|---|---|
| mixed | 8 |
| cortex-biased | 6 |
| cortex+hippocampal | 3 |
| brainstem | 2 |
| subcortical | 1 |
| subcortical+cortex | 1 |
| hippocampal | 1 |

## Panel map

| UniProt | Gene | Bias | Primary region | Note |
|---|---|---|---|---|
| P23975 | SLC6A2 | brainstem | locus coeruleus noradrenergic neurons | NET; presynaptic on LC terminals |
| Q01959 | SLC6A3 | brainstem | VTA / substantia nigra dopaminergic neurons | DAT; presynaptic on midbrain DA terminals |
| Q12879 | GRIN2A | cortex+hippocampal | cortical pyramidal + hippocampus | NMDA NR2A; mature postsynaptic |
| Q13224 | GRIN2B | cortex+hippocampal | cortical pyramidal + hippocampus + striatum | NMDA NR2B; ifenprodil-class antagonists at ATD interface |
| Q08499 | PDE4D | cortex+hippocampal | cortical pyramidal + hippocampus | cAMP PDE; broadly expressed |
| P42261 | GRIA1 | cortex-biased | cortical pyramidal | AMPA GluA1; classic cortex marker |
| P42262 | GRIA2 | cortex-biased | cortical pyramidal + cerebellum | AMPA GluA2; broadest AMPA subunit |
| P42263 | GRIA3 | cortex-biased | cortical pyramidal | AMPA GluA3 |
| P48058 | GRIA4 | cortex-biased | cerebellar Purkinje + cortical interneurons | AMPA GluA4; cerebellum-enriched |
| O43526 | KCNQ2 | cortex-biased | cortical L5/6 axon-initial-segment | Kv7.2 M-current |
| O43525 | KCNQ3 | cortex-biased | cortical L5/6 axon-initial-segment | Kv7.3 M-current; KCNQ2 partner |
| O60741 | HCN1 | hippocampal | hippocampal CA1 distal dendrites | Hyperpolarization-activated cation channel |
| P22303 | ACHE | mixed | cortex+striatum+brainstem | Broadly expressed cholinesterase |
| P08913 | ADRA2A | mixed | locus coeruleus + prefrontal cortex | α2A adrenergic; LC autoreceptor + PFC |
| P36544 | CHRNA7 | mixed | cortical L2/3 + hippocampal CA1 + basal forebrain | Nicotinic α7; CA1 enrichment per Siletti 2023 |
| O43613 | HCRTR1 | mixed | lateral hypothalamus + LC + cortical projections | Orexin receptor 1 |
| O43614 | HCRTR2 | mixed | lateral hypothalamus + tuberomammillary | Orexin receptor 2 |
| Q9Y5N1 | HRH3 | mixed | cortical + striatal H3 autoreceptors | Histamine H3 |
| Q16620 | NTRK2 | mixed | cortical L5 + hippocampal CA3 + cerebellum | TrkB BDNF receptor |
| Q99720 | SIGMAR1 | mixed | ER chaperone — broad | Sigma-1 receptor; ER chaperone-receptor |
| P21728 | DRD1 | subcortical | striatum (direct pathway MSNs) + cortical L5 | D1 dopamine receptor |
| O76083 | PDE9A | subcortical+cortex | striatum + cortical interneurons | cGMP PDE |

## Top-30 v7 ranking with brain-region annotation

| # | Compound | Tier | RRF | MAMMAL best | Brain bias | Primary region |
|---|---|---|---|---|---|---|
| 1 | d-amphetamine | positive_control | 1.185 | P23975 | brainstem | locus coeruleus noradrenergic neurons |
| 2 | methylphenidate | positive_control | 1.163 | Q01959 | brainstem | VTA / substantia nigra dopaminergic neurons |
| 3 | bupropion | extended_cns | 1.160 | Q16620 | mixed | cortical L5 + hippocampal CA3 + cerebellum |
| 4 | aniracetam | named_in_research | 1.122 | Q13224 | cortex+hippocampal | cortical pyramidal + hippocampus + striatum |
| 5 | rasagiline | extended_cns | 1.092 | Q01959 | brainstem | VTA / substantia nigra dopaminergic neurons |
| 6 | pridopidine | named_in_research | 1.085 | P08913 | mixed | locus coeruleus + prefrontal cortex |
| 7 | levetiracetam | extended_cns | 1.084 | P23975 | brainstem | locus coeruleus noradrenergic neurons |
| 8 | lisdexamfetamine | extended_cns | 1.082 | Q08499 | cortex+hippocampal | cortical pyramidal + hippocampus |
| 9 | lanicemine | named_in_research | 1.069 | Q01959 | brainstem | VTA / substantia nigra dopaminergic neurons |
| 10 | pramiracetam | extended_cns | 1.067 | P48058 | cortex-biased | cerebellar Purkinje + cortical interneurons |
| 11 | cx-717 | named_in_research | 1.066 | Q99720 | mixed | ER chaperone — broad |
| 12 | rivastigmine | positive_control | 1.064 | P48058 | cortex-biased | cerebellar Purkinje + cortical interneurons |
| 13 | selegiline | extended_cns | 1.064 | Q16620 | mixed | cortical L5 + hippocampal CA3 + cerebellum |
| 14 | cx-516 | named_in_research | 1.057 | Q01959 | brainstem | VTA / substantia nigra dopaminergic neurons |
| 15 | guanfacine | extended_cns | 1.050 | Q16620 | mixed | cortical L5 + hippocampal CA3 + cerebellum |
| 16 | piracetam | extended_cns | 1.040 | P48058 | cortex-biased | cerebellar Purkinje + cortical interneurons |
| 17 | modafinil | positive_control | 1.040 | P42261 | cortex-biased | cortical pyramidal |
| 18 | atomoxetine | positive_control | 1.035 | P48058 | cortex-biased | cerebellar Purkinje + cortical interneurons |
| 19 | donepezil | positive_control | 1.033 | P42261 | cortex-biased | cortical pyramidal |
| 20 | chembl1255723 | chembl_expanded | 1.023 | P23975 | brainstem | locus coeruleus noradrenergic neurons |
| 21 | pramipexole | extended_cns | 1.023 | Q01959 | brainstem | VTA / substantia nigra dopaminergic neurons |
| 22 | chembl1256414 | chembl_expanded | 1.020 | Q99720 | mixed | ER chaperone — broad |
| 23 | cep-26401 | named_in_research | 1.020 | P08913 | mixed | locus coeruleus + prefrontal cortex |
| 24 | chembl302231 | chembl_expanded | 1.018 | P23975 | brainstem | locus coeruleus noradrenergic neurons |
| 25 | chembl1256378 | chembl_expanded | 1.018 | P48058 | cortex-biased | cerebellar Purkinje + cortical interneurons |
| 26 | atenolol | negative_control | 1.013 | Q08499 | cortex+hippocampal | cortical pyramidal + hippocampus |
| 27 | isrib | named_in_research | 1.012 | O60741 | hippocampal | hippocampal CA1 distal dendrites |
| 28 | levodopa | extended_cns | 1.011 | Q16620 | mixed | cortical L5 + hippocampal CA3 + cerebellum |
| 29 | enalapril | negative_control | 1.011 | P42262 | cortex-biased | cortical pyramidal + cerebellum |
| 30 | rolipram | positive_control | 1.010 | P42262 | cortex-biased | cortical pyramidal + cerebellum |

## Hypothesis test

**Claim**: the V5 wet-lab top-25 spans ≥3 distinct brain-region bias categories. This is a sanity check that the pipeline isn't lopsided toward a single anatomy.

**Measured**: 4 distinct bias categories in top-25: ['brainstem', 'cortex+hippocampal', 'cortex-biased', 'mixed']
**Verdict**: PASS

---

Generated by `scripts/45_v5_brain_region.py`. Full V6 Cluster D Bayesian pipeline pending.