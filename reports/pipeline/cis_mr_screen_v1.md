# cis-MR screen: druggable genes against human cognitive performance

**A SCREEN, not a result.** A Wald-ratio cis-MR on one fine-mapped instrument generates
hypotheses. Two things must follow before any target here is believed: colocalisation, to
rule out that the eQTL and the cognition signal are distinct variants in linkage
disequilibrium, which is the dominant false-positive mode; and a pleiotropy screen.
Neither is in this file.

## Why genetics replaced the clinical literature

Three measured boundaries, not a change of taste. The DTI head ranks at chance (pooled
AUROC 0.468). The healthy-adult literature yields four modest ACUTE enhancers and zero
durable ones. The prospective registry cannot be resolved in useful time: 25 of 26 due rows
unpublished, and the single row that did resolve was retracted on attribution. A cis-eQTL
is a lifetime randomised perturbation of one gene, and it needs nobody to run, publish or
honestly report a trial.

## Data

| layer | source | n |
|---|---|---|
| outcome | Lee 2018 Cognitive Performance, GCST006572 | 257,841 |
| exposure | ROSMAP brain (DLPFC) (`QTD000434`, eQTL Catalogue SuSiE) | 560 |
| exposure | CommonMind brain (DLPFC) (`QTD000075`, eQTL Catalogue SuSiE) | 586 |
| exposure | BrainSeq brain (DLPFC) (`QTD000051`, eQTL Catalogue SuSiE) | 479 |
| exposure | GTEx brain (cortex) (`QTD000171`, eQTL Catalogue SuSiE) | 205 |
| targets | ChEMBL drug-mechanism genes | 891 (447 with an approved drug) |

Instrument rule: highest-PIP variant in the gene's credible set, PIP >= 0.2, cis-eQTL p <= 5e-08.

## Harmonisation, which is where MR silently inverts

The eQTL beta refers to the ALT allele in `chr_pos_REF_ALT`; the GWAS beta refers to its
own A1. Where A1 is the eQTL REF the GWAS beta is sign-flipped before the ratio is taken.
Getting this wrong yields a confident result with the direction reversed, which is worse
than no result. Palindromic A/T and C/G variants whose frequency is near one half cannot
be strand-resolved and are DROPPED rather than guessed.

| outcome | tests |
|---|---:|
| usable | 226 |
| dropped: palindromic variant with ambiguous frequency; strand unresolvable | 5 |

**0 gene-tissue tests at FDR < 0.05.** Ranked by p-value; `wald_ratio` is the change in cognitive performance (SD units) per unit increase in normalised brain expression, so its SIGN is the directional requirement: positive means a drug should RAISE the target's activity.

| gene | tissue | rsid | PIP | eQTL beta | CP beta (aligned) | Wald | p | q |
|---|---|---|---:|---:|---:|---:|---:|---:|
| **TIE1** | CommonMind | rs3120047 | 0.35 | +0.317 | -0.0123 | -0.0389 | 2.72e-04 | 0.061 |
| **CD40** | ROSMAP | rs4239702 | 0.88 | +0.357 | +0.0112 | +0.0314 | 1.18e-03 | 0.133 |
| **MC1R** | CommonMind | rs72813442 | 0.61 | -0.311 | +0.0105 | -0.0339 | 3.04e-03 | 0.159 |
| **CD5** | BrainSeq | rs519296 | 0.33 | -0.435 | +0.0091 | -0.0210 | 4.12e-03 | 0.159 |
| **SOST** | BrainSeq | rs1877633 | 1.00 | -0.609 | -0.0084 | +0.0138 | 5.98e-03 | 0.159 |
| **SOST** | CommonMind | rs9910625 | 0.40 | -0.359 | -0.0087 | +0.0242 | 6.38e-03 | 0.159 |
| **HTR1D** | ROSMAP | rs61777102 | 0.36 | +0.426 | +0.0085 | +0.0199 | 6.42e-03 | 0.159 |
| **F7** | CommonMind | rs2774033 | 0.74 | -0.466 | +0.0132 | -0.0283 | 6.68e-03 | 0.159 |
| **HTR1D** | GTEx | rs61777102 | 0.72 | +0.541 | +0.0085 | +0.0157 | 7.54e-03 | 0.159 |
| **F7** | BrainSeq | rs2774033 | 0.61 | -0.426 | +0.0132 | -0.0310 | 7.57e-03 | 0.159 |
| **PAH** | BrainSeq | rs9782 | 1.00 | +0.599 | -0.0083 | -0.0138 | 7.75e-03 | 0.159 |
| **GPA33** | BrainSeq | rs10918612 | 0.91 | -0.784 | -0.0113 | +0.0144 | 8.51e-03 | 0.159 |
| **PAH** | ROSMAP | rs9782 | 0.50 | +0.508 | -0.0083 | -0.0163 | 9.82e-03 | 0.159 |
| **PAH** | GTEx | rs9782 | 0.54 | +0.607 | -0.0083 | -0.0137 | 1.05e-02 | 0.159 |
| **BCHE** | CommonMind | rs2668196 | 0.51 | -0.375 | +0.0097 | -0.0258 | 1.16e-02 | 0.159 |
| **BCHE** | ROSMAP | rs9881048 | 0.35 | -0.410 | +0.0097 | -0.0235 | 1.20e-02 | 0.159 |
| **TNFRSF13C** | ROSMAP | rs6002545 | 0.53 | +0.845 | +0.0115 | +0.0136 | 1.22e-02 | 0.159 |
| **SOST** | GTEx | rs12600549 | 0.30 | -0.545 | -0.0079 | +0.0145 | 1.26e-02 | 0.159 |
| **CDK6** | CommonMind | rs42035 | 0.25 | +0.302 | -0.0091 | -0.0300 | 1.34e-02 | 0.159 |
| **CA4** | CommonMind | rs75506524 | 0.47 | -0.644 | -0.0169 | +0.0262 | 1.53e-02 | 0.173 |
| **TSHR** | CommonMind | rs9646167 | 0.87 | +0.297 | +0.0070 | +0.0237 | 1.79e-02 | 0.190 |
| **CA1** | CommonMind | rs1845891 | 0.28 | +0.240 | +0.0073 | +0.0303 | 1.85e-02 | 0.190 |
| **SSTR5** | BrainSeq | rs112392455 | 0.55 | +0.581 | +0.0223 | +0.0384 | 2.18e-02 | 0.215 |
| **SSTR5** | GTEx | rs112392455 | 0.50 | +1.305 | +0.0223 | +0.0171 | 2.36e-02 | 0.220 |
| **PSCA** | ROSMAP | rs4736371 | 0.99 | +0.780 | +0.0079 | +0.0101 | 2.43e-02 | 0.220 |
| **RGMA** | BrainSeq | rs17706224 | 0.71 | -0.472 | -0.0091 | +0.0192 | 2.66e-02 | 0.227 |
| **TPH1** | CommonMind | rs3958112 | 0.46 | -0.709 | -0.0064 | +0.0090 | 2.80e-02 | 0.227 |
| **CACNA1B** | ROSMAP | rs12001881 | 0.25 | +0.170 | -0.0067 | -0.0394 | 2.90e-02 | 0.227 |
| **SCN9A** | GTEx | rs16851931 | 0.87 | +0.756 | -0.0085 | -0.0112 | 2.91e-02 | 0.227 |
| **FKBP1A** | BrainSeq | rs6041749 | 1.00 | -0.210 | -0.0069 | +0.0327 | 3.19e-02 | 0.240 |

## Why "replicated across cohorts" is mostly an illusion here

A gene appearing in several eQTL datasets reads as replication and usually is not. All
four exposure datasets are tested against the SAME outcome GWAS. Where they select the
same variant, or variants in perfect LD, the outcome beta is literally the same number
and only the eQTL denominator changes, so the second test carries no new information
about cognition.

Of the genes appearing in more than one dataset, **9** repeat an identical
outcome effect and **40** carry more than one distinct outcome estimate. Even the
latter are cis-variants at one locus and therefore LD-correlated, so they are not
independent either. Nothing in this screen is replicated in the sense that word
normally carries; genuine replication requires a second, non-overlapping cognition
GWAS.


## The confound this does not handle

Cognitive Performance is genetically correlated with educational attainment at roughly
0.9, and EA carries social and demographic pathways with no neurobiological content. A hit
above may act on schooling rather than on cognitive machinery. The mitigation is
multivariable MR conditioning on EA (GCST006442, n = 1,131,881). That is the next thing to
build, not an optional refinement, and until it exists every row here is provisional.

## One tension to decide deliberately

Screening candidates against adverse psychiatric phenotypes is sensible, but cognition and
schizophrenia share substantial genetic architecture, so a blanket rejection of
schizophrenia-associated genes would discard some of the strongest cognition signals. That
should be a recorded per-target safety judgement, not a silent filter.

Generated by `scripts/140_cis_mr_screen.py`.
