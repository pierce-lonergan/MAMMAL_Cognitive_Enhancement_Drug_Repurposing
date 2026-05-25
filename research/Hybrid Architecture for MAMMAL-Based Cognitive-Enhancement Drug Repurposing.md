# Hybrid Architecture for MAMMAL-Based Cognitive-Enhancement Drug Repurposing on a Single RTX 5070

## TL;DR
- **Build a HYBRID architecture**: tight-couple MAMMAL's DTI head with Boltz-2 (affinity-only, "Boltzina"-style) inside a structure-aware cluster, run ADMET-AI as an independent ADMET cluster, and overlay PrimeKG/TxGNN as a knowledge-graph cluster — then fuse the three clusters at the top with Reciprocal Rank Fusion (RRF), promoted to a LightGBM LambdaMART meta-ranker once you have labeled positives.
- **The top 3 complementary classes for cognitive enhancement specifically** are: (A) structure/pocket-aware (Boltz-2 + ESM2-650M), (B) ADMET-AI, (C) PrimeKG + TxGNN. Class D (Geneformer/L1000) is recommended as a Phase 4 add-on, not a core member — the LINCS L1000 evidence base for healthy cognition is too thin to bear weight on the primary rank.
- **You cannot beat the ceiling.** Roberts CA, Jones A, Sumnall H, Gage SH, Montgomery C. "How effective are pharmaceuticals for cognitive enhancement in healthy adults? A series of meta-analyses..." European Neuropsychopharmacology 2020;38:40–62 (PMID 32709551, k=47 studies) reports pooled SMDs of methylphenidate +0.21 (recall +0.43, sustained attention +0.42, inhibitory control +0.27), modafinil +0.12 (memory updating +0.28), and null for d-amphetamine. A perfect in-silico pipeline still bottoms out on a noisy, small-effect endpoint. The win is enrichment of plausible candidates, not discovery of a "smart drug."

---

## Key Findings

1. **MAMMAL's structural blindness is the single biggest gap for cognition targets.** Roughly half the 22-target panel is dominated by allosteric pharmacology (α7 nAChR type-I/II PAMs at the TMD subunit interface; PDE4D UCR2-directed partial inhibitors; AMPA-receptor LBD modulators at the cyclothiazide site; GluN2B ifenprodil-site negative allosteric modulators; HCN1 ivabradine-like blockers). A sequence-only DTI head cannot resolve these. **Class A is non-negotiable.**

2. **ADMET-AI dominates ADMET coverage in the open-source space.** Per Swanson, Walters & Zou, Bioinformatics 2024 (doi:10.1093/bioinformatics/btae416): "ADMET-AI has the best average rank among all models that have been evaluated on the 22 datasets in the TDC ADMET Leaderboard." It exposes 41 endpoints (10 regression + 31 classification) and runs in seconds on CPU for hundreds of SMILES. It is the only single open package that simultaneously covers BBB, P-gp, CYP1A2/2C9/2C19/2D6/3A4, hERG, DILI, Caco-2, PAMPA, VDss, t½ — all first-order kill criteria for chronic healthy-adult dosing. **Class B is non-negotiable.**

3. **PrimeKG + TxGNN provide the mechanism layer MAMMAL completely lacks.** PrimeKG (Chandak, Huang & Zitnik, *Scientific Data* 2023) has 129,375 nodes and 4,050,249 edges across 30 relation types. TxGNN (Huang K, Chandak P, Wang Q, et al., *Nature Medicine* 2024 Dec;30(12):3601–3613, doi:10.1038/s41591-024-03233-x) is trained on this graph and "rank[s] drugs as potential indications and contraindications across 17,080 diseases," beating 8 baselines by 49.2% on indication and 35.1% on contraindication accuracy under zero-shot evaluation. The cognition use case — "rank compounds whose target/pathway neighborhood overlaps with cognition-relevant phenotypes" — is exactly what this class was built for.

4. **Hybrid > tight > loose.** Pure loose coupling wastes information (ESM2 protein embeddings should *augment* MAMMAL's protein side, not just feed a separate ranker; Boltz-2 affinity prediction already uses ESM-style protein features internally). Pure tight coupling forces premature embedding alignment between modalities (sequence ↔ graph ↔ ADMET) that have no joint training corpus on a single RTX 5070. The right answer is *cluster-tight, system-loose*: tight inside each of (A/B/C), loose at the cross-cluster rank fusion layer.

5. **The 5070 is sufficient if you sequence carefully.** Peak per-stage VRAM with bfloat16 inference: MAMMAL (~2–3 GB), ESM2-650M (~2.5 GB), Boltz-2 affinity-only mode (~7–8 GB on L40S per Nebius's published benchmark, fits 12 GB), Boltz-2 structure mode (≤1,000 residues on RTX 4070-class 12 GB per UCSF ChimeraX benchmarks Dec 2025), ADMET-AI (CPU-friendly, <1 GB), TxGNN inference (<2 GB). No two large models run concurrently; load/unload between phases.

6. **Throughput estimate for 22 × 300 = 6,600 pairs**: MAMMAL DTI ~10 min total; ESM2 embeddings cached once for 22 targets (<5 min); ADMET-AI for 300 compounds <2 min; Boltz-2 structure prediction for the 7 targets without good co-crystals (~30–60 min, cached); Boltz-2/Boltzina affinity for top-50 compounds × 22 targets ≈ 1,100 calls × ~20 s = ~6 hours; PrimeKG/TxGNN scoring <10 min. **Full cold-cache pipeline: ~8 hours wall-clock; warm-cache refresh: ~30–60 minutes per 50 new compounds.**

7. **RRF is the right default; promote to LambdaMART once you have ≥20 labeled positives.** Cormack, Clarke & Buettcher (SIGIR 2009) showed "Reciprocal Rank Fusion Outperforms Condorcet and Individual Rank Learning Methods" when relevance labels are sparse and rankers are heterogeneous — exactly your regime. Seal et al. (J Cheminformatics 2013, "Enhanced ranking of PknB inhibitors using data fusion methods") replicated this for virtual screening with structure + ligand-based methods.

---

## Details

### 1. Architectural framework — decision and rationale

**Loose coupling (rejected as primary).** Run MAMMAL, Boltz-2, ADMET-AI, TxGNN as black boxes, collect per-compound rankings, fuse. Pros: fault-tolerant, modular, easy to version. Cons: ignores the fact that ESM2 protein embeddings are a strict superset of MAMMAL's sequence representation for protein-side reasoning, and that Boltz-2 affinity prediction already uses ESM-style protein features internally. You'd duplicate work and lose cross-model signal.

**Tight coupling (rejected as primary).** Force everything into one embedding space — contrastively align ESM2 protein embeddings, MoLFormer molecular embeddings, and PrimeKG node embeddings, then train a joint head. This is the BioBridge / multimodal-alignment direction. Cons on a 5070: alignment training requires triplets and gradient memory you don't have, and the alignment quality on a 22-target panel will be statistically meaningless. You're optimizing engineering complexity for a problem you don't have.

**Hybrid (RECOMMENDED).** Three model-class **clusters**, each internally tight-coupled where it makes physical sense, joined at the top by rank fusion.

```
                        ┌─────────────────────────────┐
                        │   INPUT LAYER               │
                        │  SMILES list (~300)         │
                        │  Target panel (22, FASTA)   │
                        └────────────┬────────────────┘
                                     │
        ┌────────────────────────────┼─────────────────────────────┐
        │                            │                             │
        ▼                            ▼                             ▼
┌───────────────────┐     ┌────────────────────┐     ┌──────────────────────┐
│ CLUSTER A         │     │ CLUSTER B          │     │ CLUSTER C            │
│ Structure/pocket  │     │ ADMET / safety     │     │ Mechanism / KG       │
│                   │     │                    │     │                      │
│ ESM2-650M ──┐     │     │ ADMET-AI (41 EP)   │     │ PrimeKG subgraph     │
│             ├──►  │     │  ├─ BBB            │     │  ├─ drug-target-     │
│ MAMMAL DTI ─┘     │     │  ├─ P-gp           │     │  │   pathway-        │
│  (pKd grid)       │     │  ├─ CYP3A4/2D6     │     │  │   phenotype       │
│      │            │     │  ├─ hERG           │     │  └─ side-effect      │
│      ▼            │     │  ├─ DILI           │     │      neighborhood    │
│ Boltz-2 affinity- │     │  └─ Caco-2/PAMPA   │     │ TxGNN zero-shot      │
│  only (Boltzina)  │     │                    │     │  indication score    │
│  for top-50 cmpds │     │                    │     │                      │
└─────────┬─────────┘     └─────────┬──────────┘     └──────────┬───────────┘
          │                         │                            │
          │ rank_A(c,t)             │ rank_B(c)                  │ rank_C(c,t)
          │ + pose feasibility      │ + hard gates              │ + provenance
          │                         │                            │
          └────────────┬────────────┴───────────────┬────────────┘
                       ▼                            ▼
            ┌─────────────────────┐      ┌─────────────────────┐
            │  HARD GATES         │      │  RANK FUSION        │
            │  - hERG > 0.7 → cut │      │  Stage 1: RRF       │
            │  - BBB < 0.3 → cut  │      │  Stage 2: LambdaMART│
            │  - DILI > 0.8 → cut │      │  (when labels avail)│
            └──────────┬──────────┘      └──────────┬──────────┘
                       └─────────────┬──────────────┘
                                     ▼
                       ┌──────────────────────────┐
                       │  RANKED CANDIDATE LIST   │
                       │  + per-model provenance  │
                       │  + Sanity-gate flag      │
                       │  + Polypharmacology summ │
                       └──────────────────────────┘
```

**Why this works on a 5070.** Cluster B is CPU-friendly; runs in parallel to GPU work. Cluster C inference is graph lookup + a small GNN forward; <2 GB VRAM. Cluster A has the only heavy GPU spend, and you sequence it: MAMMAL DTI on all 6,600 pairs first (~10 min), then load ESM2 to embed the 22 targets once (cached), unload everything, then load Boltz-2 for affinity on the top-N after hard ADMET gating.

**Decision tree.**
- <100 compounds, no labels, no time → **pure loose RRF**. (Ship in 2 days.)
- Curated panel + ≥20 positives + structural targets dominate → **hybrid as above**. (Your situation.)
- Single target, ample labels, want SOTA on one number → **tight coupling within one cluster**, skip fusion. (Not your case.)

### 2. Top 3 model classes for cognitive enhancement (decision)

**Recommended: A (Structure) + B (ADMET) + C (KG).**

The complementarity argument for the cognition use case specifically:

- **MAMMAL gives you binding-affinity ranking.** Its 458M-parameter multimodal architecture, trained on >2 billion biological samples across proteins/small-molecules/single-cell gene expression, achieves SOTA on 9/11 downstream benchmarks (Shoshan et al., arXiv 2410.22367; npj Drug Discovery 2026). On DTI it gives you a pKd estimate. That's it. It cannot tell you:
  1. Whether the predicted binder fits a real pocket — half your cognition targets care about *allosteric* sites, not orthosteric. **→ Class A required.**
  2. Whether the compound crosses the BBB, gets effluxed by P-gp, blocks hERG, induces CYP3A4. For a chronic, healthy-adult dosing scenario these are 80% of attrition. **→ Class B required.**
  3. Whether the target sits in a cognition-relevant pathway in the first place vs. a spurious hit. **→ Class C required.**

- **Why structure beats transcriptomics here.** Cognition targets are mostly ligand-gated ion channels, GPCRs, transporters, and a phosphodiesterase — *not* transcriptionally driven phenotypes. The LINCS L1000 panel is dominated by cancer cell lines (MCF7, PC3, A375, HA1E), with very limited representation of human neurons, no cognitive-task readout, and only ~1,000 landmark transcripts. CMap reproducibility for drug repositioning was independently shown to be limited (Lin et al., Nature Sci Rep 2021, "Evaluation of connectivity map shows limited reproducibility in drug repositioning"). **Class D is high-effort, low-yield for this specific endpoint; defer it.**

- **Why KG beats LLM-RAG here.** A literature-RAG layer (Class E) would re-extract claims already in PrimeKG (which integrates DrugBank, ChEMBL, OpenTargets, MeSH, DisGeNET, etc.). For a 22-target curated panel, the KG approach gives you structured, queryable, provenance-tracked evidence with no hallucination surface. Add the LLM later as an *explanation* layer over the KG paths, not as an evidence source.

**Trade-off acknowledged for what's left out.** Skipping Class D means no phenotypic complement to target-affinity — if a compound works via a polypharmacological mechanism hitting a target *not* in your panel, you won't see it. Mitigation: in Phase 4 add a CMap/LINCS L1000 query against a curated "cognition-relevant gene set" (CREB targets, BDNF/TrkB downstream, LTP-associated transcripts) using `signatureSearcher` or `metaLINCS` and bolt it on as a 4th cluster.

### 3. Per-class technical detail

#### Class A — Structure / pocket-aware: ESM2-650M + Boltz-2 (affinity-only)

- **ESM2-650M** (Meta, MIT license). HuggingFace `facebook/esm2_t33_650M_UR50D`. 650M params, 33 layers, 1280-dim embeddings. Peak VRAM at bf16 inference: ~2.5 GB for sequences ≤1024 residues at batch 16; the original Meta benchmark on A100 establishes the upper bound. FlashAttention-equipped ESME (Hallee et al., Cell Reports Methods 2025) runs lighter. On the 5070 you will batch one target at a time and cache embeddings to parquet — only 22 targets to embed.

- **Boltz-2** (MIT license, Wohlwend/Passaro et al., bioRxiv 2025.06.14.659707, June 18 2025; MIT + Recursion). GitHub `jwohlwend/boltz`. PyPI `boltz`. Approaches FEP accuracy at ~1000× speed. **Sized for 12 GB**: per UCSF ChimeraX Boltz documentation (T. Goddard, updated Dec 15 2025), Boltz 2.1.1 on an RTX 4070 (12 GB, Windows) handles up to ~1,000 residues in full structure mode; an RTX 3070 (8 GB) caps at ~700 residues. UCSF doc states verbatim: *"Consumer Nvidia GPUs with 8 or 12 GB of memory (e.g. RTX 3070) only handle 300–500 residues before using CPU memory on Windows that slows the prediction speed by 10–20 fold."* For the 22 cognition targets, all single chains fit (largest are GluN1/GluN2A ~880-residue dimer chain, just under the cap). **Nebius's L40S benchmark** reports ~11 GB for structure prediction and ~7–8 GB for affinity prediction, so affinity-only mode comfortably fits 12 GB.

- **Boltzina mode** (Furui & Ohue, arXiv 2508.17555, Aug 24 2025): skip Boltz-2's diffusion structure module and feed AutoDock Vina poses straight to Boltz-2's affinity head. Boltzina paper reports it "achieved up to 11.8× faster through reduced recycling iterations and batch processing." Accuracy is below full Boltz-2 but well above Vina/GNINA on the MF-PCBA benchmark. This is your production mode after Phase 1.

- **Throughput**: Furui & Ohue 2025 cite Boltz-2 as "requiring approximately 20 seconds per compound per GPU" (their motivating quote, referring to H100); expect 2–3× slower on the 5070 (~40–60 s in full mode, ~5–10 s in Boltzina mode after the 11.8× speedup). LMI4Boltz (Litfin et al., bioRxiv 2025.10.29.684571) reports the low-memory branch reduces a 911-residue prediction from 28 min to 10 min on RTX 3070 — proof that the inference is tractable on consumer hardware with chunking enabled. The biorxiv preprint quotes: *"Boltz-2 supports predictions with size up to ∼1600 tokens using a GPU with 24 GB VRAM"*; LMI4Boltz extends this by 66.7%.

**Role in the pipeline.** ESM2-650M produces per-residue and per-protein embeddings used as augmenting features for the meta-ranker and for ChEMBL nearest-neighbor protein search. Boltz-2 (or Boltzina) produces a pose-conditioned affinity (`affinity_pred_value` = log10(IC50 in µM)) and a binary binder probability (`affinity_probability_binary` ∈ [0,1]) — the latter is the calibrated "is this even plausibly a binder" signal.

**Cognition target structural inventory (verified PDB IDs):**

| Target | Representative PDB | Ligand class |
|---|---|---|
| CHRNA7 (α7 nAChR) | 7EKT (EVP-6124 + PNU-120596 PAM), 7KOX, 7KOQ; recent ago-PAM GAT107 (Cell Discovery 2025) | Type-II PAM at TMD interface — critical |
| ACHE | 4EY7 (donepezil) | Orthosteric + peripheral |
| GRIA1 | 3SAJ (NTD only) | Limited; use GluA2 homolog |
| GRIA2 | 1LBC, 3TKD (LBD + cyclothiazide); 4U5B (full receptor + CTZ) | CTZ allosteric site — key |
| GRIA3 | 3DLN; PMC12422969 (GluA3-G–TARPγ2 cryo-EM with CTZ, Nature 2025) | Allosteric |
| GRIA4 | 3KFM (LBD + kainate); cryo-EM GluA4/TARPγ2 + glutamate + CTZ (bioRxiv 2025.06.12.659357) | Allosteric |
| GRIN2A | 7EU7 (S-ketamine), 5VII (GluN2A-selective antagonist) | Mixed |
| GRIN2B | 4PE5 (ifenprodil, full receptor), 3QEL (ATD + ifenprodil), 5IOV (Ro25-6981) | NAM at ATD — key |
| DRD1 | 7CKW (fenoldopam), 7LJD (dopamine + LY3154207 PAM), 7LJC (SKF81297 + LY3154207) | PAM site important |
| SLC6A3 (DAT) | dDAT 4XP1 (cocaine), 4XPG (β-CFT); hDAT cryo-EM (Wang/Gouaux Nature 2024) | Orthosteric + S2 |
| ADRA2A | 6KUY (RES), 7W6P (Gαo + biased agonist) | Mixed |
| SLC6A2 (NET) | 8ZP1 (reboxetine), 8ZP2 (atomoxetine), 8ZOY (NE) | Orthosteric |
| HRH3 | 7F61 (PF-03654746) | Orthosteric antagonist |
| HCRTR1 | 4ZJ8 (suvorexant), 6V9S (JH112) | Orthosteric antagonist |
| HCRTR2 | 4S0V (suvorexant), 5WQC (EMPA) | Orthosteric antagonist |
| PDE4D | **6NJJ (BPN14770), 6NJH, 6NJI** (Burgin/Gurney 2019); 3G45 PDE4B (Burgin 2010 seminal) | UCR2 allosteric — key |
| PDE9A | 3K3E (BAY73-6691), 4GH6, 8BPY | Orthosteric |
| NTRK2 | 4AT3 (TrkB kinase + CPD5N), 4ASZ apo | Kinase ATP-site |
| SIGMAR1 | 5HK1 (PD144418), 5HK2 (4-IBP); haloperidol/NE-100/pentazocine series | Orthosteric |
| KCNQ2 | 7CR2 (retigabine), 7CR0 apo, 8XO1 (QO-83) | Pore opener |
| KCNQ3 | No homomeric ligand-bound PDB; use 7BYM (KCNQ4+retigabine) homolog | — fold via Boltz-2 |
| HCN1 | 8Y60 (ivabradine), 5U6O apo, 6UQF (S4-down) | Pore blocker |

For KCNQ3 and any pocket without a co-crystal, predict structure with Boltz-2 (single-chain runs ≤500 residues, well within 12 GB), then run Boltzina-style affinity-only inference against the predicted pocket.

#### Class B — ADMET: ADMET-AI

- **ADMET-AI** (Swanson K, Walters WP, Zou J. "ADMET-AI: a machine learning ADMET platform for evaluation of large-scale chemical libraries." *Bioinformatics* 2024;40(7):btae416, doi:10.1093/bioinformatics/btae416, MIT license). GitHub `swansonk14/admet_ai`. CLI `admet_predict`. Trained on 41 TDC datasets (10 regression, 31 classification). Per the paper: "ADMET-AI has the best average rank among all models that have been evaluated on the 22 datasets in the TDC ADMET Leaderboard." Inference is CPU-fast — full DrugBank-percentile context for 300 SMILES runs in <2 minutes.
- **VRAM: negligible** (CPU model); leaves GPU for MAMMAL/Boltz.

**Cognition-specific hard gates** (kill filters applied before rank fusion):

| Endpoint | Threshold | Rationale |
|---|---|---|
| BBB penetration (prob) | <0.30 → cut | Cognition target ⇒ must reach CNS |
| P-gp substrate (prob) | >0.85 → cut | P-gp efflux kills CNS exposure for many cognitive enhancers |
| hERG inhibition (prob) | >0.70 → cut | Cardiotox risk; chronic dosing context |
| DILI (prob) | >0.80 → cut | Hepatotox risk; chronic dosing |
| CYP3A4 inhibition (prob) | >0.85 → flag (not cut) | DDI risk; flag for review |
| Ames mutagenicity | >0.85 → cut | Hard kill |
| Caco-2 permeability | <-5.5 logPapp → flag | Oral bioavailability concern |

**Cognition-specific weighted ADMET score** (used as a feature in the meta-ranker after gating):
```
ADMET_score(c) = 0.35*BBB + 0.20*(1 - hERG) + 0.15*(1 - P_gp_substrate)
               + 0.10*(1 - DILI) + 0.10*Caco2_norm + 0.10*(1 - CYP3A4_inhib)
```

#### Class C — Knowledge graph: PrimeKG + TxGNN

- **PrimeKG** (Chandak P, Huang K, Zitnik M. *Scientific Data* 2023, doi:10.1038/s41597-023-01960-3). GitHub `mims-harvard/PrimeKG`. 129,375 nodes (10 types: drug, gene/protein, disease, phenotype, pathway, biological process, molecular function, cellular component, anatomy, exposure) × 4,050,249 edges across 30 relation types. Harvard Dataverse persistent DOI 10.7910/DVN/IXA7BM.
- **TxGNN** (Huang K, Chandak P, Wang Q, et al. "A foundation model for clinician-centered drug repurposing." *Nature Medicine* 2024 Dec;30(12):3601–3613, doi:10.1038/s41591-024-03233-x, Epub 2024 Sep 25). GitHub `mims-harvard/TxGNN`. Trained on PrimeKG. Zero-shot drug-disease scoring with metric-learning module; from the abstract: "rank drugs as potential indications and contraindications across 17,080 diseases." Reports +49.2% on indication and +35.1% on contraindication accuracy vs. 8 baselines under zero-shot evaluation. Inference VRAM <2 GB.

**Cognition subgraph extraction strategy.** Because "healthy cognitive enhancement" is not a disease node in PrimeKG, build a *virtual phenotype anchor* by taking the union of these disease nodes' k-hop neighborhoods (k=2):
- Mild cognitive impairment (MCI)
- Alzheimer's disease (used cautiously as a mechanistic proxy, not a clinical proxy)
- Attention-deficit/hyperactivity disorder (for processing-speed / working-memory targets)
- Fragile X syndrome (anchor for PDE4D)
- Narcolepsy (anchor for HCRTR1/2 and HRH3)

For each candidate compound:
- Score 1: TxGNN indication probability against the virtual anchor (mean of the 5 disease scores).
- Score 2: TxGNN contraindication probability — subtract.
- Score 3: PrimeKG path count between compound and any of the 22 targets via relations in {`drug_protein`, `protein_protein`, `pathway_protein`, `bioprocess_protein`}. Use Personalized PageRank or a Katz-style decayed walk.
- Score 4: Side-effect overlap with cognition-degrading side-effect terms (sedation, somnolence, cognitive impairment, anticholinergic side effects) — *penalty*.

**KG cluster output** is a composite `KG_score(c) = w1*indication − w2*contraindication + w3*log(1+path_count) − w4*side_effect_penalty`, with weights initialized (0.4, 0.3, 0.2, 0.1) and tuned on the positive-control panel.

### 4. Concrete hybrid architecture spec

**Sequencing (load/unload on 12 GB)**:

```
Phase 0 (one-time, cached):
  - ESM2-650M → embeddings for 22 targets → cache .pt files          [~5 min, 2.5 GB VRAM]
  - Boltz-2 → predicted structures for 7 targets w/o co-crystal      [~30-60 min, ~8-10 GB VRAM]
  - PrimeKG load + TxGNN model load (CPU + RAM)                       [3 min, no GPU]

Phase 1 (per-compound, fast):
  - MAMMAL DTI 22×300 grid                                            [~10 min, 3 GB VRAM]
  - ADMET-AI 300 SMILES                                               [~2 min, CPU]
  - TxGNN scoring 300 compounds                                       [~5 min, 1 GB VRAM]
  - Apply ADMET hard gates → typically removes 30-50% of compounds

Phase 2 (per-surviving-compound × target, expensive):
  - Boltzina affinity for top-N pairs (e.g., MAMMAL top-50 surviving) [~6 hr, 7-8 GB VRAM]

Phase 3 (fusion):
  - RRF over (MAMMAL, Boltzina, TxGNN, ADMET_score) per (c, t)        [seconds, CPU]
  - Aggregate to per-compound score across targets                     [seconds]
  - Polypharmacology view (target-count-weighted)                      [seconds]

Phase 4 (validation):
  - Sanity-gate positive controls (donepezil, galantamine, modafinil,
    methylphenidate, atomoxetine, BPN14770, varenicline, vortioxetine)
    must rank in top decile or pipeline fails                          [seconds]
```

**Caching strategy.** Memoize every expensive computation as parquet keyed by content hash:
- `cache/esm2/<sha1(seq)>.pt` — ESM2 embedding per target
- `cache/boltz_struct/<sha1(seq)>.cif` — predicted structure
- `cache/boltzina/<sha1(seq)>__<sha1(smiles)>.json` — affinity + binder prob
- `cache/admet/<sha1(smiles)>.parquet` — 41-EP ADMET vector
- `cache/txgnn/<sha1(smiles)>.json` — KG scores

**VRAM management.** PyTorch `torch.cuda.empty_cache()` between phases; explicit `del model; gc.collect()`. Never load Boltz-2 and ESM2 simultaneously. Use `PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True` (LMI4Boltz default) to avoid fragmentation. Enable LMI4Boltz chunking flags (`triangle_mult_gate_nchunks=4`, `chunk_size_transition_z=32`, `chunk_size_tri_attn=64`) for any target sequence approaching 1,000 residues.

**Estimated wall-clock for 22×300 problem on RTX 5070**:
- Cold cache (first run): ~8–10 hours, dominated by Boltz-2 structure prediction (one-time per target) and Boltzina affinity calls.
- Warm cache (subsequent runs with new compounds): ~30–60 minutes per batch of 50 new compounds.

### 5. Rank fusion methodology

**Stage 1: Reciprocal Rank Fusion (RRF).** No labels needed. For each (compound c, target t) pair, given K ranker outputs:

```
RRF_score(c, t) = Σ_k  1 / (k_const + rank_k(c, t))
```

with `k_const = 60` (Cormack et al. 2009 default; the value is famously insensitive in the 30–80 range). This is your **default for Phase 1–2** and what you ship.

Python pseudocode:

```python
import pandas as pd

def rrf(rankings: dict[str, pd.Series], k: int = 60) -> pd.Series:
    """
    rankings: {model_name: Series indexed by compound_id, values = rank (1 = best)}
    returns: Series of RRF scores, higher = better
    """
    score = pd.Series(0.0, index=next(iter(rankings.values())).index)
    for name, r in rankings.items():
        score = score.add(1.0 / (k + r), fill_value=0.0)
    return score.sort_values(ascending=False)
```

Why RRF: Cormack, Clarke & Buettcher SIGIR 2009 showed RRF beats Condorcet and individual learning-to-rank baselines when relevance labels are sparse and rankers are heterogeneous. Drug-discovery follow-ups (Seal et al., J Cheminformatics 2013) replicated this on virtual screening with structure + ligand-based methods. Willett's analysis (Mol Inf 2013) attributes RRF's effectiveness to "the close relationship between the reciprocal rank of a database structure and its probability of activity."

**Stage 2: LambdaMART (LightGBM/XGBoost ranker).** Once you have ≥20 positive controls + a comparable number of "expected-inactive" decoys, train a supervised meta-ranker:

```python
import lightgbm as lgb

features = [
    "mammal_pkd", "mammal_polyphar_n",
    "boltzina_logIC50", "boltzina_binder_prob",
    "admet_bbb", "admet_pgp", "admet_herg", "admet_dili",
    "admet_cyp3a4", "admet_cyp2d6", "admet_caco2",
    "txgnn_indication", "txgnn_contraindication",
    "kg_path_count", "kg_side_effect_penalty",
    "esm2_cos_to_canonical_binder",
    "morgan_cos_to_canonical_binder"
]
ranker = lgb.LGBMRanker(
    objective="lambdarank",
    metric="ndcg",
    label_gain=[0, 1, 3],
    n_estimators=300, num_leaves=31, learning_rate=0.05
)
ranker.fit(X_train, y_train, group=groups_train)
```

Use leave-one-positive-out CV (few positives expected). Don't overfit hyperparameters — 200–400 trees, num_leaves 15–31, learning_rate 0.05.

**Calibration.** Score distributions across rankers are wildly heterogeneous (Boltz-2 logIC50 is signed real-valued; MAMMAL pKd is real-valued; TxGNN indication is [0,1]; ADMET-AI is [0,1]). RRF is rank-based and immune; LambdaMART handles raw features fine because trees are scale-invariant. Avoid CombSUM-style score-sum methods unless you Platt-scale every input to [0,1] first.

**Why NOT Borda / Condorcet.** Borda is sensitive to ties in long tails; Condorcet is expensive and gives no improvement over RRF in the Cormack benchmark. CombMNZ is competitive but requires score normalization that RRF avoids.

### 6. Integration with existing MAMMAL pipeline

**Existing artifacts (v1)**:
- `dti_grid.parquet`: 22 × 300 predicted pKd, scaled scores
- `sanity_gate.parquet`: positive-control retrieval check
- `polypharm.parquet`: target-count weighted aggregation
- planned: BBBP/ClinTox MAMMAL heads (auxiliary), ChEMBL backstop, OpenTargets cross-reference

**Augmentation map (per new cluster)**:

| New cluster | Augments / replaces | New artifacts | New decision gate |
|---|---|---|---|
| ESM2-650M | Augments MAMMAL protein side; feeds ChEMBL NN | `cache/esm2/*.pt`, `target_esm_embeddings.parquet` | None (feature only) |
| Boltz-2 / Boltzina | Augments DTI grid with structural feasibility | `boltz_structs/*.cif`, `boltzina_affinity.parquet` | `pose_feasible` boolean (binder_prob > 0.3) |
| ADMET-AI | **Replaces** v1's BBBP/ClinTox heads (drop them); provides 41 EPs | `admet_predictions.parquet`, `admet_gates_passed.parquet` | Hard gates (§3 Class B) |
| PrimeKG + TxGNN | Provides mechanism complement; subsumes OpenTargets v2 plan | `kg_scores.parquet`, `kg_paths.parquet` (provenance) | `kg_evidence_present` boolean |
| Fusion layer | New | `rrf_ranking.parquet`, `lambdamart_ranking.parquet`, `final_ranking.parquet` | Sanity-gate must place ≥6/8 positive controls in top 30% |

**Migration**: each cluster is additive; v1 MVP continues to ship `dti_grid.parquet` unchanged. The fusion stage *reads* dti_grid + new artifacts and *writes* `final_ranking.parquet`. v1 results remain comparable across versions because the underlying dti_grid is byte-for-byte identical.

**File-tree additions** for `C:\Users\Pierce Lonergan\Documents\GitHub\MAMMAL_Cognitive_Enhancement_Drug_Repurposing\`:

```
src/
  models/
    mammal_dti.py          (existing)
    esm2_embed.py          (NEW)
    boltz_runner.py        (NEW)
    boltzina_affinity.py   (NEW)
    admet_ai_runner.py     (NEW)
    txgnn_runner.py        (NEW)
    kg_score.py            (NEW; PrimeKG path/PPR scoring)
  fusion/
    rrf.py                 (NEW)
    lambdamart.py          (NEW)
    calibration.py         (NEW)
  gates/
    admet_gates.py         (NEW)
    sanity_gate.py         (existing → extended)
  pipeline/
    run_phase0_cache.py    (NEW; ESM + Boltz struct cache)
    run_phase1_fast.py     (NEW; MAMMAL + ADMET + TxGNN)
    run_phase2_boltzina.py (NEW; expensive affinity)
    run_phase3_fusion.py   (NEW; rank fusion)
data/
  cache/                   (NEW; memoized expensive ops)
  kg/primekg/              (NEW)
artifacts/
  v1/                      (existing)
  v2/                      (NEW)
configs/
  thresholds.yaml          (NEW; ADMET gate thresholds)
  weights.yaml             (NEW; fusion weights)
```

### 7. Critical assessment and failure modes

**Where the hybrid genuinely adds value:**
1. **Allosteric pharmacology coverage.** MAMMAL alone will completely miss the α7-PAM and PDE4D-UCR2 allosteric chemotypes — they don't look like orthosteric binders. Boltz-2 with pocket-conditioning fixes this. **High-value addition.**
2. **CNS-safety triage.** ADMET-AI hard gates remove ~30–50% of compounds before expensive Boltz-2 calls. **High-value, near-zero cost.**
3. **Mechanism explanations**. TxGNN/PrimeKG paths give you a publishable "why" alongside each candidate. **High-value, low cost.**

**Where it's complexity for complexity's sake (be honest):**
1. **Adding LLM literature reasoning before you have ranked candidates.** Speculative. Defer to Phase 5.
2. **Training a meta-ranker before ≥20 positives.** RRF is near-ceiling at low label counts. Don't build LambdaMART in Phase 2 just because you can.
3. **Embedding alignment between modalities (BioBridge-style).** Sounds elegant; will not produce better candidates on a 22-target panel with 5070 compute.

**Model disagreement protocol (concrete rules):**
- MAMMAL says high pKd (≤7 nM predicted), Boltzina says no plausible binding pose (`affinity_probability_binary` < 0.3): **trust Boltzina, downrank**. Sequence-only DTI is more failure-prone for allosteric pockets.
- Boltzina says strong binder, MAMMAL says no: **trust Boltzina** if structure is high-confidence (pLDDT ≥ 70 in pocket); otherwise flag.
- TxGNN says strong indication but no target in the panel is hit: **flag for "off-panel mechanism" review** — interesting, not a failure.
- ADMET gate fails but everything else is strong: **do not bypass**. ADMET failures are physical, not statistical.

**Calibration challenges.**
- Boltz-2 affinity is calibrated on PDBBind/MF-PCBA, both biased toward orthosteric pockets. Allosteric calibration is unproven — treat `affinity_pred_value` as ordinal, not as a real IC50.
- MAMMAL pKd is calibrated on BindingDB. Watch for the ~1-log over-prediction for novel scaffolds (a known DTI failure mode across the literature).
- TxGNN was trained for *disease* indications; using it for the cognition virtual anchor is an off-label use. Validate by injecting known cognitive enhancers (donepezil, memantine, BPN14770) and confirming they rank in TxGNN's top decile against the anchor.

**Compute bottlenecks.**
- Phase 2 (Boltzina) is the wall-clock dominator. Mitigation: only call it on the top-50–80 surviving compounds after ADMET gates.
- ESM2-650M embedding for 22 targets is one-time; not a bottleneck.
- LMI4Boltz chunking is your safety valve if any target sequence pushes past 1,000 residues — enable by default.

**The fundamental ceiling — honest take.** Roberts CA, Jones A, Sumnall H, Gage SH, Montgomery C, *European Neuropsychopharmacology* 2020;38:40–62 (k=47 studies, n≈3,000 healthy adults) reports pooled SMDs of **methylphenidate +0.21** (driven by recall +0.43, sustained attention +0.42, inhibitory control +0.27), **modafinil +0.12** (memory updating +0.28), and **null for d-amphetamine**. From the abstract: *"There is a user perception that these drugs are effective cognitive enhancers, but this is not supported by the evidence so far."* Battleday & Brem's earlier modafinil review reached the same conclusion. Even with a flawless in-silico pipeline, you are searching a space where the *best known agents* deliver d≈0.2. Your pipeline will not deliver a "smart drug"; it can credibly enrich for plausible candidates with mechanistic justification, defensible ADMET, and structural feasibility. **Frame the deliverable as enrichment + provenance, not discovery.**

### 8. Execution plan

**Phase 1 (Week 1) — ADMET first.** Integrate ADMET-AI as Cluster B. Replace v1's BBBP/ClinTox heads. Implement the hard-gate stage. **Why first**: highest signal-to-effort ratio, no GPU contention with MAMMAL, immediately filters ~30–50% of the 300-compound list, makes everything downstream cheaper. **Validation gate**: positive controls (donepezil, modafinil, methylphenidate, atomoxetine, galantamine, memantine, BPN14770) must all pass BBB gate and ≤1 may fail hERG.

**Phase 2 (Week 2–3) — Structure cluster.** Integrate ESM2-650M (target-side embeddings, cached). Integrate Boltz-2 structure mode for the 7 targets without good co-crystals (KCNQ3, GRIA1, GRIA3, NTRK2 ECD if needed, etc.). Add Boltzina affinity-only inference for top-50 MAMMAL hits per target. **Validation gate**: re-score known α7-PAMs (PNU-120596, EVP-6124), PDE4D inhibitors (BPN14770, D159687/MK-8189), AMPA potentiators (CX-516, IDRA-21) and confirm Boltzina places them in top 20% for their respective targets.

**Phase 3 (Week 3–4) — KG cluster.** Load PrimeKG, install TxGNN. Build the virtual cognition phenotype anchor. Compute compound-level KG scores. Add the RRF fusion layer over (MAMMAL, Boltzina, TxGNN, ADMET_score). **Validation gate**: top-30 final-rank compounds must include ≥6 of 8 positive controls; if not, debug weights/anchor.

**Phase 4 (Week 5) — Promote to LambdaMART** if you have ≥20 labels by then. Otherwise stay on RRF and ship. **STOP adding complexity here unless empirical validation demands more.**

**Phase 5 (optional, Week 6+) — Transcriptomic and LLM layers.** Add a curated LINCS L1000 query against a cognition-relevant gene set as Cluster D *only if* RRF is converged and the team has cycles. Add Llama-3.1-8B (4-bit quantization, ~6 GB VRAM) over PubMedBERT-retrieved abstracts as an *explanation* layer over KG paths — not as a new evidence source.

**Stop conditions (when to stop building, start validating empirically):**
- Sanity gate stable: ≥6/8 positive controls in top 30, ≥4/8 in top 10
- RRF ranking reproducible across seeds (Spearman ρ > 0.95)
- ADMET hard gates remove >25% of candidates
- Top-20 list reviewed by a pharmacologist (or by Pierce with explicit mechanistic justification per candidate)

**Publication-worthy methodology contributions** (realistic):
1. **A reproducible 458M-MAMMAL + Boltz-2(Boltzina) + ADMET-AI + TxGNN hybrid for CNS repurposing on a single consumer GPU** — no open-source replication of this stack exists in the published literature. The single-RTX-5070 budget is itself the novelty.
2. **A virtual "healthy cognition" phenotype anchor for KG-based repurposing**, with explicit overlap of MCI/ADHD/Fragile-X/narcolepsy neighborhoods. Clean methodological contribution; benchmarkable.
3. **A calibration analysis of MAMMAL pKd vs. Boltzina affinity vs. TxGNN scores** on a cognition target panel. Useful as a benchmarking paper even if no candidate ever enters in vitro testing.

Frame any publication around methodology + benchmark, not around the candidate list. The candidate list will not survive peer review without wet-lab confirmation, which is outside scope.

---

## Recommendations

1. **Build the hybrid in this order: ADMET-AI → ESM2+Boltzina → PrimeKG/TxGNN → RRF fusion → (later) LambdaMART.** Do not start with the structure cluster; the ADMET gate makes everything downstream cheaper.
2. **Use Boltzina mode by default**, not full Boltz-2. The up-to-11.8× speedup (Furui & Ohue 2025) is worth the small accuracy loss on the 5070's compute budget.
3. **Ship RRF, not LambdaMART, for v2.** Promote to LambdaMART only after ≥20 labels exist.
4. **Apply ADMET hard gates BEFORE rank fusion**, not after. Gates are physical kill criteria, not relative rankings.
5. **Cache aggressively to parquet** keyed by content hash. The pipeline must be cheap to re-run; assume Pierce will re-run it dozens of times.
6. **Defer transcriptomic (Class D) and LLM (Class E) clusters to Phase 5+.** Both are high-effort, low-yield for this specific endpoint.
7. **Frame the publication as methodology + benchmark, not discovery.** The healthy-adult cognitive-enhancement effect-size ceiling (SMD ≤ 0.43 for the best known agents on specific sub-domains, ≤ 0.21 pooled) makes "we found a smart drug" claims unsupportable.

**Thresholds that would change these recommendations:**
- **24 GB GPU available**: switch from Boltzina back to full Boltz-2 in structure+affinity mode; consider tight-coupling MAMMAL+ESM2 via fine-tuning.
- **Neuron-relevant L1000 dataset** (iPSC-derived neurons, NPC panels): elevate Class D from Phase 5 to Phase 3.
- **≥100 labeled positives**: skip RRF and start with LambdaMART/XGBoost ranker.
- **Experimental validation budget secured**: drop Phase 5 entirely and pour cycles into hit-to-lead rather than more in silico layers.

---

## Caveats

1. **Boltz-2 / Boltzina on 12 GB is at the edge of supported configuration.** NVIDIA NIM requires ≥48 GB VRAM for production Boltz-2 deployment (NIM documentation). The CLI/Python distribution works on 12 GB up to ~1,000 residues for structure prediction (UCSF ChimeraX Boltz documentation, T. Goddard, updated Dec 15 2025), but expect to use LMI4Boltz chunking and `expandable_segments:True` allocator for the largest targets. Affinity-only mode is more comfortable (~7–8 GB on L40S per Nebius's published benchmark).
2. **Boltzina-style affinity-only inference accuracy is below full Boltz-2.** Furui & Ohue 2025 (arXiv 2508.17555) report Boltzina as "below Boltz-2" on MF-PCBA but well above AutoDock Vina/GNINA. For your relative ranking purpose this is fine; for absolute IC50 quotation it is not.
3. **TxGNN was not trained for cognitive enhancement.** Its 17,080-disease coverage does not include a "healthy cognitive enhancement" node. The virtual anchor approach is reasonable but not validated. Treat TxGNN scores as informative-but-noisy.
4. **PrimeKG side-effect data is biased toward post-market reports** (DrugBank, SIDER ingestion). Cognition-degrading side effects of well-characterized drugs are well covered; novel scaffolds will have empty side-effect neighborhoods and look artificially safe.
5. **The 300-compound list is curated by Pierce.** Curation bias is the largest, unmeasured source of variance in the eventual top-20. Document curation criteria explicitly; consider augmenting with DrugBank approved-drugs + ChEMBL Phase-3 layer to dilute curation effects.
6. **5070 Blackwell sm_120 PyTorch compatibility**: as of mid-2025, the Blackwell consumer architecture requires PyTorch 2.5+ with CUDA 12.6+ for stable inference; verify before installing Boltz/cuEquivariance kernels (cuEquivariance v0.5+ recommended for triangle attention acceleration on Blackwell).
7. **Effect-size ceiling**: Roberts et al., *European Neuropsychopharmacology* 2020;38:40–62 (k=47 studies, n≈3,000): methylphenidate SMD=0.21, modafinil SMD=0.12, d-amphetamine null. Your pipeline cannot deliver candidates exceeding this ceiling in healthy adults; it can deliver candidates with credibly higher *prior* probability of joining this small-effect club.
8. **PDB IDs for SLC6A3/DAT human cocaine cryo-EM structures** (Wang/Gouaux Nature 2024) should be confirmed at rcsb.org by Pierce before pipeline integration; the dDAT entries (4XP1, 4XPG, etc., Penmatsa & Gouaux) are well-established backups. KCNQ3 has no homomeric ligand-bound PDB entry as of pre-2026 retrieval; predict with Boltz-2 or use KCNQ4 (7BYM) as homolog template.