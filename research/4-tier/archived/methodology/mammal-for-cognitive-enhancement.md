# MAMMAL for Healthy Cognitive Enhancement Drug Repurposing — A Technical Deep-Dive

## TL;DR

- **MAMMAL is a 458M-parameter T5-style multi-modal encoder-decoder (Shoshan et al., arXiv 2410.22367; *npj Drug Discovery* 3:14, 2026) that hits state-of-the-art on 9 of 11 drug-discovery benchmarks, but none of those benchmarks measure cognition, behavior, or anything closer to "working memory" than IC50 in a tumor cell line — so for healthy cognitive enhancement, MAMMAL is a target-affinity oracle and a transcriptomic response simulator, not a cognition predictor.**
- **Realistically, Pierce can run MAMMAL inference and LoRA-style fine-tuning on a single RTX 5070 (12 GB) — the base model is ~1.8 GB at fp32 / ~0.9 GB at bf16 and IBM ships a working Colab T4 notebook — but the value will come from using MAMMAL to score known nootropic chemotypes against a curated panel of ~18 cognition-relevant human targets (CHRNA7, GRIA1-4, HRH3, PDE4D, PDE9A, DRD1, ADRA2A, SIGMAR1, NTRK2/TrkB, KCNQ2/3, HCN1, EIF2B subunits, etc.) using the released `dti_bindingdb_pkd` head, then triangulating with LINCS L1000 connectivity-map signatures — not by querying MAMMAL for "working memory enhancers" as a label, because that label does not exist in its training data.**
- **The most defensible MAMMAL-surfaceable cognition-relevant repurposing candidates are the ones whose targets are well-represented in BindingDB/PubChem/STRING: galantamine (CHRNA7 PAM + AChEi), donepezil (AChE), methylphenidate (DAT/NET), atomoxetine (NET), modafinil (DAT-adjacent), pitolisant (HRH3/SIGMAR1), zatolmilast/BPN14770 (PDE4D), and clemastine (M1/M3 + oligodendrocyte differentiation). The genuinely interesting but speculative class — ISRIB / 2BAct / DNL343 (EIF2B activators) — will likely score poorly on a vanilla MAMMAL DTI query because eIF2B is not a classic ligand-binding protein in the training distribution; that's a model limitation, not biological irrelevance.**

## Key Findings

1. **MAMMAL is publicly available, weight-released, Apache-2.0 licensed, and Pierce-runnable today.** GitHub: `BiomedSciAI/biomed-multi-alignment` (44 stars, 19 forks); weights: `ibm/biomed.omics.bl.sm.ma-ted-458m` (HF). Install: `pip install biomed-multi-alignment[examples]`. Tested on Python ≥3.10, PyTorch ≥2.0. The team also ships an MCP server (`mammal_mcp`) that exposes 4 task tools — PPI, protein solubility, DTI (BindingDB pKd), and TCR-epitope binding — over STDIO or Streamable-HTTP on port 8001, with Claude Desktop and Ollama/MCPHost configs documented.

2. **458M parameter T5-derived encoder-decoder, hidden-dim 768, trained on 2B samples sourced from 6 datasets — UniRef90 (proteins), OAS (antibodies), ZINC + PubChem (small molecules), CELLxGENE Census (single-cell expression), and STRING (PPI pairs)** — using **7 pretraining tasks** (span-masking LM on each modality, antibody denoise, PPI classification, PPI generation). Pretraining used 16 × A100-80GB across 2 nodes for ~3 months under PyTorch FSDP.

3. **None of MAMMAL's 11 benchmark tasks measure cognition, behavior, or central nervous system endpoints.** Its closest CNS-relevant capability is BBBP (blood-brain-barrier permeability classification, AUROC 0.957 ± 0.006 vs MoLFormer 0.937), which only tells you whether the molecule enters the brain, not what it does once there. MAMMAL has **no DDI benchmark**, no behavioral endpoint, no neurological disease label, and was not trained on any cognition-stratified data.

4. **MAMMAL has demonstrated drug repurposing capability in cancer.** In the npj Drug Discovery paper, MAMMAL correctly predicted the relative potency ranking of four oncology drugs not in its GDSC training data (carfilzomib > nintedanib > infigratinib > vemurafenib, three with Tanimoto < 0.7 to training compounds), prospectively validated with CellTiter-Glo. This is the only published prospective repurposing prediction. **There is no published cognition-relevant repurposing prediction from the MAMMAL group.**

5. **The fundamental scoping reality:** "Healthy working memory enhancement" is not a label in any database MAMMAL was trained on. Pretraining sources contain (i) protein sequences, (ii) molecule SMILES, (iii) PPI edges, (iv) single-cell expression — none encode "improves digit span by N%". The legitimate MAMMAL-shaped questions are therefore narrower: *Does compound X bind target Y? Does compound X reverse a transcriptomic signature Z? Does this protein-protein edge exist?* You map cognition to MAMMAL queries by curating the target list yourself.

## Details

### 1. MAMMAL architecture and capabilities

**Citation.** Shoshan Y, Raboh M, Ozery-Flato M, Ratner V, Golts A, Weber JK, Barkan E, Rabinovici-Cohen S, Polaczek S, Amos I, Shapira B, Hazan L, Ninio M, Ravid S, Danziger MM, Shamay Y, Kurant S, Morrone JA, Suryanarayanan P, Rosen-Zvi M, Hexter E. "MAMMAL — Molecular Aligned Multi-Modal Architecture and Language for biomedical discovery." *npj Drug Discovery* 3:14 (May 2026). Preprint: arXiv:2410.22367 (v1 Oct 2024, v3 May 2025). Corresponding author: yoels@il.ibm.com. IBM Research–Israel (Haifa) + TJ Watson + Technion.

**Architecture.** Transformer encoder-decoder explicitly inspired by T5 (Raffel et al.) and Vaswani et al. Three architectural innovations vs. vanilla T5:

- **Dual-mode encoder.** Shared-weight encoder stack can be used (a) encoder-only with token + scalar prediction heads on the encoder's final hidden states (classification/regression), or (b) encoder-decoder autoregressive with residual injection of encoder final-layer hidden states into each decoder layer (generation, variable-length output).
- **Modular tokenizer.** Distinct sub-tokenizers for amino acids, SMILES atoms, gene-name vocabulary, and special meta-tokens are stitched into one unified ID space. Meta-tokens like `<@TOKENIZER-TYPE=AA>`, `<@TOKENIZER-TYPE=SMILES@MAX-LEN=10>`, `<MOLECULAR_ENTITY_GENERAL_PROTEIN>`, `<SEQUENCE_NATURAL_START>` partition the prompt across modalities. Code path: `fuse/data/tokenizers/modular_tokenizer/`.
- **Continuous scalar embeddings.** Numerical values (e.g., expression levels, binding affinities) are projected through a learned linear layer into the 768-dim embedding space and added to the token embedding sequence — no quantization buckets, no float-as-string. There's a dedicated encoder scalar prediction head for regression outputs.

Hidden dimension is 768. Exact encoder/decoder layer count, attention heads, and FFN dimensions are **not stated in the paper or model card** — Pierce will need to inspect `config.json` in the HF safetensors. The "ted" in `ma-ted-458m` is "T5-Encoder-Decoder."

**Pretraining corpus (6 datasets, ~2B samples, 3 modalities).**

| Dataset | Modality | Pretraining task | Samples |
|---|---|---|---|
| UniRef90 | General proteins | Span-masking LM | 180M |
| OAS (Observed Antibody Space) | Antibody sequences | LM + denoise | 650M (each task) |
| ZINC + PubChem | Small molecule SMILES | Span-masking LM | 200M |
| CELLxGENE Census | scRNA-seq | Cell-genes LM | 30M |
| STRING | PPI edges | Classification + generation | 780M + 390M |

Optimizer: AdamW (β1=0.9, β2=0.999), weight decay 0.01, gradient clip 1.0, 2K warmup, cosine decay to 10% peak LR.

**Benchmarks (11 tasks; SOTA on 9, comparable on 2).** All numbers below come directly from the published benchmark table.

| # | Benchmark | Metric | Prior SOTA | MAMMAL |
|---|---|---|---|---|
| 1 | Cell type (Zheng68k PBMC) | F1 | 0.710 (scBERT) | **0.763 ± 0.012** |
| 2 | BBBP (MoleculeNet) | AUROC | 0.937 (MoLFormer) | **0.957 ± 0.006** |
| 3 | ClinTox (avg) | AUROC | 0.948 (MoLFormer) | **0.986 ± 0.007** |
| 4 | Cancer-Drug Resp 1 (GDSC1) | Pearson | 0.887 (TxLLM) | **0.917 ± 0.001** |
| 5 | Cancer-Drug Resp 2 (GDSC2) | Pearson | 0.900 (TxLLM) | **0.931 ± 0.002** |
| 6 | Cancer-Drug Resp 3 (DeepCDR) | Pearson | 0.923 (DeepCDR) | 0.928 ± 0.000 (comparable) |
| 7 | Antibody infilling (SAbDab/dyMEAN) | CDRH3-AAR | 0.375 | **0.446 ± 0.002** |
| 8 | Ab-Ag binding (HER2/Mason) | AUROC | 0.924 (AbLang) | 0.928 ± 0.002 (comparable) |
| 9 | TCR β-epitope binding (Weber/TITAN) | AUROC | 0.862 (TITAN) | **0.879 ± 0.003** |
| 10 | PPI ΔΔG (SKEMPI S1131) | Pearson | 0.663 (seq-only AttABseq) | **0.852 ± 0.041** |
| 11 | DTI (PEER/BindingDB pKd, 4 held-out target classes) | NRMSE↓ | 0.942 ± 0.028 (PEER) | **0.906 ± 0.011** |

The cognition-relevant items here are #2 (BBBP), #11 (DTI), and arguably #4-6 (cell-line response — if you treat a neural cell line as your "target tissue" by fine-tuning). **There is no behavioral, organ-level, or in vivo benchmark.**

**Drug repurposing handling.** MAMMAL doesn't have a dedicated "repurposing" head. The released repurposing demonstration uses the DTI (pKd) and Cancer-Drug-Response (IC50) heads in two modes: (a) score a fixed compound against a library of proteins (what new targets does this approved drug hit?), or (b) score a fixed protein against a compound library (what approved drugs hit my target?). This is *target-affinity ranking* repurposing, not *signature-reversal* repurposing — for the latter, you'd combine MAMMAL's cell-genes LM head with LINCS L1000 externally.

**Availability and licensing.**
- Code: https://github.com/BiomedSciAI/biomed-multi-alignment (Apache 2.0, 44 stars, 19 forks, PyPI: `biomed-multi-alignment`)
- Base weights: https://huggingface.co/ibm-research/biomed.omics.bl.sm.ma-ted-458m (also mirrored at `ibm/...`)
- Fine-tuned heads on HF (`ibm-research/` namespace):
  - `…ma-ted-458m.dti_bindingdb_pkd` — drug-target pKd regression
  - `…ma-ted-458m.protein_solubility` — solubility classifier
  - `…ma-ted-458m.moleculenet_bbbp` — BBB permeability
  - `…ma-ted-458m.moleculenet_clintox_tox` — clinical toxicity
  - `…ma-ted-458m.moleculenet_clintox_fda` — FDA approval probability
  - `…ma-ted-458m.tcr_epitope_bind` — TCR-epitope binding
- MCP server: `mammal_mcp/` subdirectory; toggles 4 tools via env flags; STDIO or HTTP transport.
- Colab quickstart: `tutorials/begginer_inference.ipynb` (works on free T4 / 16 GB).

**Compute footprint and your RTX 5070.** Base model is **~1.8 GB at fp32, ~0.9 GB at bf16**. With KV cache + activations for typical DTI prompts (target AA sequence ≤ 1500 tokens + drug SMILES ≤ 200 tokens), realistic inference VRAM is comfortably under 4 GB. **A 12 GB RTX 5070 will run base inference, batched inference of ~32-64 (protein, drug) pairs at once, and LoRA fine-tuning of the encoder/decoder on a modest dataset (~10k-100k pairs) without offload.** Full-parameter fine-tuning is more marginal — likely possible at small batch (1-4) with gradient checkpointing and bf16; if VRAM is tight, use PEFT/LoRA on the encoder. IBM does not publish an official VRAM figure, so I'm extrapolating from the parameter count and the fact that their Colab tutorial runs on a free 16 GB T4. The MCP README explicitly cautions "not using anymore than 2 models at one time," corroborating a per-model footprint of a few GB.

### 2. Practical workflow for cognitive-enhancement repurposing

#### Step 1 — Environment setup (~15 min on a clean Ubuntu)

```bash
conda create -n mammal_env python=3.10 -y
conda activate mammal_env
conda install pytorch pytorch-cuda=12.1 -c pytorch -c nvidia
pip install biomed-multi-alignment[examples]
# Optional: pip install fastmcp if you want the MCP server
```

Then download weights lazily on first use (`Mammal.from_pretrained("ibm/biomed.omics.bl.sm.ma-ted-458m")` triggers the HF download — ~1.8 GB).

#### Step 2 — Build the cognition target panel

This is the most important design decision and it lives **outside MAMMAL**. You're committing to a hypothesis about which proteins, when modulated, plausibly improve working memory, processing speed, or learning rate in healthy humans. Pierce's prior list is solid; here's the curated panel I'd actually feed into MAMMAL queries, with UniProt accessions for direct prompt construction:

| Mechanism class | Gene | UniProt | Rationale |
|---|---|---|---|
| Cholinergic | CHRNA7 (α7 nAChR) | P36544 | Encenicline, galantamine PAM; attention/WM |
| Cholinergic | ACHE | P22303 | Donepezil, rivastigmine, galantamine |
| Glutamatergic | GRIA1-4 (AMPA) | P42261-3, P48058 | Ampakines, CX-717, tulrampator |
| Glutamatergic | GRIN2A/2B (NMDA) | Q12879, Q13224 | Memantine, rapastinel |
| Dopaminergic | DRD1 | P21728 | D1 agonists (DAR-0100A), processing speed |
| Dopaminergic | SLC6A3 (DAT) | Q01959 | Methylphenidate, modafinil-adjacent |
| Noradrenergic | ADRA2A | P08913 | Guanfacine, prefrontal WM |
| Noradrenergic | SLC6A2 (NET) | P23975 | Atomoxetine |
| Histaminergic | HRH3 | Q9Y5N1 | Pitolisant; arousal + WM |
| Orexinergic | HCRTR1/2 | O43613, O43614 | Wake-promoting; suvorexant antagonism reverse |
| Phosphodiesterases | PDE4D | Q08499 | BPN14770/zatolmilast; cAMP/memory |
| Phosphodiesterases | PDE9A | O76083 | BI 409306, PF-04447943 |
| BDNF/TrkB | NTRK2 | Q16620 | 7,8-DHF, LM22A-4 |
| Sigma-1 | SIGMAR1 | Q99720 | Fluvoxamine, blarcamesine, pridopidine |
| Channels | KCNQ2/3 | O43526, O43525 | Retigabine, XEN-1101 |
| Channels | HCN1 | O60741 | Ivabradine analogs |
| ISR | EIF2B1-5 | Q14232, Q14202, O00303, Q9NR50, Q13144 | ISRIB; **caveat: not a classical druggable receptor — MAMMAL DTI will likely underscore.** |
| Mitochondrial | NDUFA-family / mtETC | various | Methylene blue, MitoQ |

#### Step 3 — Build the compound library

Curate a ~200-500 compound library: all FDA-approved CNS-active drugs, well-known nootropic chemotypes, ampakines, ISR modulators, sigma ligands, GLP-1 agonists (for the emerging cognitive-data angle), and a negative-control set (peripheral-only drugs that should not score on CNS targets). Pull SMILES from DrugBank or PubChem.

#### Step 4 — Run DTI scoring (the main MAMMAL call)

The `dti_bindingdb_pkd` head outputs predicted pKd (higher = stronger affinity). Reference code from the HF model card, adapted:

```python
from mammal.model import Mammal
from fuse.data.tokenizers.modular_tokenizer.op import ModularTokenizerOp
from mammal.keys import *
import torch

model = Mammal.from_pretrained("ibm-research/biomed.omics.bl.sm.ma-ted-458m.dti_bindingdb_pkd").eval().cuda()
tok = ModularTokenizerOp.from_pretrained("ibm-research/biomed.omics.bl.sm.ma-ted-458m.dti_bindingdb_pkd")

def score(target_aa: str, drug_smiles: str) -> float:
    sample = {ENCODER_INPUTS_STR:
        f"<@TOKENIZER-TYPE=AA><BINDING_AFFINITY_VALUE><SENTINEL_ID_0>"
        f"<MOLECULAR_ENTITY><MOLECULAR_ENTITY_GENERAL_PROTEIN>"
        f"<SEQUENCE_NATURAL_START>{target_aa}<SEQUENCE_NATURAL_END>"
        f"<@TOKENIZER-TYPE=SMILES>"
        f"<MOLECULAR_ENTITY><MOLECULAR_ENTITY_SMALL_MOLECULE>"
        f"<SEQUENCE_NATURAL_START>{drug_smiles}<SEQUENCE_NATURAL_END><EOS>"}
    tok(sample_dict=sample, key_in=ENCODER_INPUTS_STR,
        key_out_tokens_ids=ENCODER_INPUTS_TOKENS,
        key_out_attention_mask=ENCODER_INPUTS_ATTENTION_MASK)
    # batch, run forward, decode scalar pKd from regression head
    # (Exact tensor plumbing in mammal/examples/dti_bindingdb_pkd.py)
    ...
```

Loop this over your 18 × ~300 ≈ 5,400 (target, drug) pairs. At ~50 ms/pair on a 5070 with batch-32, this is single-digit minutes wall-clock. Output: ~5,400 predicted pKd values.

#### Step 5 — Rank and triangulate

- **Sanity check:** the top hits for ACHE should include donepezil, rivastigmine, galantamine (positive controls). If not, the prompt formatting is wrong.
- **Cross-target promiscuity score:** count how many panel targets a compound binds with pKd > 6. Polypharmacology is generally *good* for cognition (modafinil, galantamine, atomoxetine all hit ≥2 panel targets), and bad for safety. Track both.
- **BBBP gate:** run the `moleculenet_bbbp` head on every compound, filter for predicted BBB+. Anything not crossing the BBB is filtered out unless you have a prodrug strategy.
- **ClinTox/FDA gate:** run the two ClinTox heads to get a basic safety prior. A drug that scores high on ClinTox-FDA-approval is "looks like an approved drug" — useful for repurposing prioritization.
- **External triangulation (essential):**
  - **LINCS L1000 / Connectivity Map (CMap):** pull the L1000 signatures of your compound list (clue.io), compute signature reversal against published "youthful working memory" or "long-term-potentiation" signatures. MAMMAL cannot do this end-to-end yet (no signature-distance head), so this is a separate step.
  - **OpenTargets / DrugBank:** check whether existing literature already implicates the predicted (compound, target) pair. Strong corroboration = green light.
  - **ChEMBL bioactivity:** for any compound where MAMMAL predicts pKd > 6 at a target with no ChEMBL evidence, that's either a novel insight or a hallucination — treat with skepticism.

#### Step 6 — Constructing "working memory / processing speed / learning rate" queries

You cannot directly ask MAMMAL "what improves working memory." What you can do:

- **Proxy 1 (target-based):** define each cognitive endpoint as a weighted target panel.
  - *Working memory* ≈ HRH3 + ADRA2A + DRD1 + CHRNA7 + GRIN2B (prefrontal/working memory circuits)
  - *Processing speed* ≈ DAT/SLC6A3 + NET/SLC6A2 + HRH3 + HCRTR (catecholamine/arousal)
  - *Learning rate / LTP* ≈ GRIA1-4 + GRIN2A + NTRK2 + EIF2B + PDE4D (synaptic plasticity)
  - Score each compound on the weighted panel; rank.
- **Proxy 2 (signature-based, via cell-genes LM):** fine-tune MAMMAL's cell-genes head on a small dataset of brain-region pseudo-bulk profiles from young vs older healthy donors, or on hippocampal expression before/after LTP induction in mice. Then use the model to predict which compound transcriptionally shifts a brain pseudo-bulk profile toward the "young/LTP-induced" state. Requires fine-tuning data Pierce would need to assemble (CELLxGENE has the source data) and is non-trivial.
- **Proxy 3 (cancer-cell-line drug-response repurposing):** misuse the GDSC head — feed it expression profiles from neural cell lines (SH-SY5Y, iPSC-derived cortical neurons) instead of cancer cells. The model has never seen this, so outputs are out-of-distribution and untrustworthy without empirical validation. Don't lead with this.

#### Step 7 — Output interpretation

The pKd head is a regression head, so output is a continuous number. Calibration is the open question — NRMSE 0.906 on the held-out PEER test is decent but not amazing on absolute scale. A pKd prediction of 7.5 means "MAMMAL thinks Kd ≈ 30 nM"; you should not believe the absolute number to better than ±0.5-1 log units. **Rank order is more trustworthy than absolute values.** Always interpret in z-score within your panel, not as standalone Ki estimates.

### 3. Compound survey — what MAMMAL would likely flag

This is a hybrid analysis: where MAMMAL's training data lets it score these compounds confidently vs. where it cannot. Effect sizes in healthy adults are from the published RCT literature.

#### Cholinergic (well-represented in MAMMAL training; high confidence in DTI scores)

- **Donepezil** (AChE, Kd ~6 nM). Healthy adults: mixed; one famous Yesavage 2002 pilot study in pilots showed flight-simulator training retention improvement, but subsequent attempts did not replicate. Expected MAMMAL behavior: high pKd at ACHE, low at CHRNA7. Will surface but evidence in healthy is weak.
- **Galantamine** (AChE Ki ~360 nM + CHRNA7 PAM). MAMMAL should pick up the AChE binding, will likely *miss* the α7 allosteric modulation (PAMs are notoriously hard for orthosteric-trained DTI models — the binding site is different and the assay endpoint isn't simple Kd). In healthy elderly (Hänsel et al., crossover), an effect on attention was reported only in a subset stratified by basal forebrain volume.
- **Rivastigmine** (AChE + BuChE). Will surface for AChE/BCHE. Healthy-adult cognitive data effectively nonexistent.

#### Alpha-7 nicotinic PAMs / partial agonists (under-represented; medium confidence)

- **Encenicline (EVP-6124)** — α7 partial agonist + 5-HT3 antagonist, half-life >50h. In schizophrenia (Keefe et al. 2015, *Neuropsychopharmacology* 40:3053-3060), CogState OCI improvements at 0.27 and 0.9 mg QD; *clinical development halted ~2015 due to GI adverse events*. Single-dose studies in healthy male volunteers showed beneficial effects on cognition with mild side-effect profile. **Galantamine remains the only approved compound in this mechanistic family.**
- **TC-5619** — exploratory schizophrenia trial showed working-memory effect *only in tobacco users*. Healthy data limited to Phase 1 safety.
- MAMMAL handling: α7 nAChR is in STRING/UniProt; small-molecule partial agonists are in PubChem. DTI head will produce predictions. Allosteric modulation is structurally distinct from orthosteric binding and BindingDB's pKd labels conflate the two — expect noise.

#### AMPA modulators / ampakines (medium-low confidence)

- **CX-717, CX-516, tulrampator (S-47445)** — RespireRx / Servier ampakines. Cortex Pharmaceuticals' CX-717 showed positive effects on sleep-deprivation cognition in healthy adults (DARPA-funded studies, Porrino et al. 2005), then development stalled. Tulrampator is in MDD trials. **No durable healthy-adult RCT signal of clinically meaningful magnitude.**
- **Aniracetam** — racetam, no compelling healthy-adult RCT data despite popularity.
- MAMMAL handling: GRIA1-4 receptors are large multi-subunit channels; ampakines are positive allosteric modulators. Same caveat as α7 PAMs.

#### NMDA modulators

- **Memantine** (low-affinity NMDA antagonist). In healthy: not a cognitive enhancer; mild impairment in some tasks.
- **Lanicemine / rapastinel** — failed in MDD; no healthy data.
- MAMMAL handling: will surface NMDA binding; whether the direction (antagonism vs partial agonism) maps to cognition is target-label-not-effect.

#### Dopaminergic / arousal stimulants (best-established healthy-adult cognition class)

The most rigorous meta-analyses are Roberts CA, Jones A, Sumnall H, Gage SH, Montgomery C, "How effective are pharmaceuticals for cognitive enhancement in healthy adults? A series of meta-analyses…" *Eur Neuropsychopharmacol* 2020 Sep;38:40-62 (PMID 32709551):
- **Modafinil** — overall SMD = 0.12, p = 0.01; memory updating SMD = 0.28, p = 0.03; small effect, mostly on executive subdomains.
- **Methylphenidate** — overall SMD = 0.21, p = 0.0004; recall SMD = 0.43, p = 0.0002; sustained attention SMD = 0.42, p = 0.0004. **The largest, most consistent healthy-adult signal in the field.**
- **D-amphetamine** — heterogeneous, smaller effect sizes.
- **Atomoxetine** — selective NET; some healthy-adult evidence for response inhibition.
- MAMMAL handling: DAT (SLC6A3), NET (SLC6A2), DRD1, DRD2 are well-represented in BindingDB. Methylphenidate/atomoxetine binding should be predicted accurately. Modafinil's mechanism is partly DAT-mediated but with substantial uncharacterized polypharmacology (orexin/histamine indirect effects) — MAMMAL will catch the DAT piece, miss the rest.

#### Histamine H3 antagonists

- **Pitolisant** — FDA-approved for narcolepsy 2019. The PEACE study (NCT05849675) is actively testing acute pitolisant 36 mg on n-back working memory + Affective Go/NoGo in healthy adults; results not yet reported in my searches. Mechanism includes σ-1 agonism, which MAMMAL may or may not register.
- **CEP-26401** (Cephalon) — in healthy volunteers (Baakman et al. 2019, *Br J Clin Pharmacol*): subtle improvement in spatial working memory at the *lowest* dose (20 µg), inverted-U dose-response. Compared favorably to modafinil and donepezil in arousal/alertness measures.
- MAMMAL handling: HRH3 is a small molecule druggable GPCR, well represented. Expect reliable affinity prediction.

#### Orexin agents
- **Suvorexant / lemborexant** — antagonists, marketed for insomnia. Cognitive effects in healthy are mostly negative (sedation residue).
- Orexin agonists (TAK-925/danavorexton) are emerging for narcolepsy. Healthy-adult cognition data scant.

#### Integrated stress response inhibitors — the speculative interesting case

- **ISRIB** — Sidrauski/Walter discovery; reverses age-related spatial memory deficit and ameliorates working memory in old mice (Krukowski et al. 2020, *eLife* 9:e62048); also rescues neuronal/dendritic spine changes after TBI (Frias et al. 2022, *PNAS*). **No human cognition RCT to date** — Calico/AbbVie's clinical candidate (likely 2BAct/DNL343 class) is in early development. Drug-likeness historically poor; ISRIB itself is a tool compound.
- **DNL343** (Denali, EIF2B activator) — in Phase 2-3 for ALS; cognitive endpoints not primary.
- MAMMAL handling: **EIF2B subunits are not classic druggable receptors and the ISRIB binding site at the EIF2B decameric interface is poorly represented in BindingDB.** Vanilla DTI scoring will likely underrank ISRIB — this is a model blind spot, not a biological one. To capture this with MAMMAL you'd need to fine-tune on a custom EIF2B-ISRIB-analog dataset, which doesn't really exist publicly.

#### BDNF/TrkB modulators
- **7,8-Dihydroxyflavone (7,8-DHF)** and **LM22A-4** — both controversial TrkB "agonists" with ongoing debate about whether they truly activate TrkB or work through off-target mechanisms (Boltaev et al. 2017, *Sci Signal*). No healthy-adult RCT cognition data.
- MAMMAL handling: NTRK2 ECD binding pocket isn't well annotated; expect unreliable scoring.

#### Sigma-1 agonists
- **Fluvoxamine** — SSRI with σ-1 agonism. Healthy cognition effects neutral/sedating.
- **Pridopidine** — Huntington's failed Phase 3, σ-1 mechanism. No healthy data.
- **Blarcamesine (ANAVEX 2-73)** — Anavex; mixed Alzheimer's results. No healthy data.
- MAMMAL handling: SIGMAR1 is a small membrane protein, ligands in BindingDB. Affinity predictions should be okay; cognition mapping is the issue.

#### KCNQ / HCN channel modulators
- **Retigabine** (Kv7 opener) — withdrawn for AED side effects; XEN-1101 in development.
- **Ivabradine** (HCN) — cardiac use; cognitive effects minimal in healthy.
- MAMMAL handling: ion channels are tricky — large multi-subunit complexes, often allosteric binding. Limited training signal.

#### PDE inhibitors
- **Rolipram** — selective PDE4, intolerable nausea, abandoned.
- **BPN14770 / zatolmilast** — PDE4D-selective negative allosteric modulator. Phase 1 in healthy elderly suggested working/immediate memory benefit (Alzheimer's Drug Discovery Foundation Cognitive Vitality summary). Phase 2 in Fragile X — Berry-Kravis EM et al., "Inhibition of phosphodiesterase-4D in adults with fragile X syndrome: a randomized, placebo-controlled, phase 2 clinical trial." *Nat Med* 2021;27(5):862–870 (doi: 10.1038/s41591-021-01321-w) — was a 24-week randomized, placebo-controlled, two-way crossover trial in 30 adult males aged 18–41 with FXS, 25 mg BID; primary endpoint was safety/tolerability (met) with significant cognitive improvement on NIH-Toolbox Oral Reading Recognition (LSMean +2.80, p=0.0157), Picture Vocabulary (+5.79, p=0.0342), and Cognition Crystallized Composite (+5.29, p=0.0018). FDA Orphan + Fast Track for FXS; Shionogi acquired.
- **PF-04447943** (PDE9) — Pfizer, Alzheimer's failed.
- MAMMAL handling: PDE4D / PDE9A active sites are well-characterized in PDB/BindingDB. **PDE4D-NAM (allosteric) is harder than PDE4D-orthosteric.** Likely to surface rolipram and orthosteric PDE4 inhibitors; may underscore BPN14770 due to allosteric mechanism.

#### Mitochondrial / metabolic
- **Low-dose methylene blue** (~1-4 mg/kg). Healthy-adult RCT — Rodriguez P, Zhou W, Barrett DW, Altmeyer W, Gutierrez JE, Li J, Lancaster JL, Gonzalez-Lima F, Duong TQ. "Multimodal Randomized Functional MR Imaging of the Effects of Methylene Blue in the Human Brain." *Radiology* 2016 Nov;281(2):516–526 (doi:10.1148/radiol.2016152893): 280 mg single oral dose, n=26 healthy adults aged 22–62, double-blind placebo-controlled, **7% increase in correct responses during memory retrieval (P = .01)** on delayed match-to-sample, with increased bilateral insular cortex fMRI BOLD during sustained attention (Z=2.9–3.4, P=.01–.008). **Biphasic dose-response — pro-oxidant above ~4 mg/kg.** This is one of the only nootropic mechanisms with a published, positive healthy-adult fMRI + behavior RCT.
- **MitoQ, urolithin A** — limited cognition data in healthy.
- MAMMAL handling: methylene blue's primary mechanism is mitochondrial electron-cycling, not target-receptor binding — MAMMAL DTI won't capture it. You'd need to use the cell-genes/L1000 path.

#### Remyelination
- **Clemastine** — first-gen H1 antihistamine, also M1/M3 muscarinic antagonist, drives oligodendrocyte differentiation. ReBUILD trial (Green et al. 2017, *Lancet* 390:2481-2489) in MS: shortened P100 VEP latency by 1.7 ms/eye (95% CI 0.5–2.9; p=0.0048). No healthy-adult cognition data per se. Of interest because remyelination is a plausible substrate for processing-speed improvement.
- MAMMAL handling: muscarinic targets well-represented; M1/M3 binding will surface. The oligodendrocyte-differentiation phenotypic mechanism will not be captured by DTI alone.

#### Glutamate modulators
- **Riluzole / troriluzole** — Na-channel / glutamate release modulator; ALS approval; cognition in healthy weak.

#### GLP-1 agonists (emerging)
- **Semaglutide** — Phase 3 EVOKE + EVOKE+ in early Alzheimer's reported 1-year follow-up at CTAD 2025: did **not** beat placebo on primary cognitive endpoint, though safety confirmed. In MDD (Mansur et al. 2025) semaglutide failed primary executive-function endpoint but improved global cognition in secondary analyses. **OxSENSE** is testing single-dose 0.5 mg semaglutide in healthy volunteers on cognition + reward sensitivity.
- **Liraglutide** — earlier MDD pilot showed executive function improvements (McIntyre lab).
- MAMMAL handling: GLP1R is a class B GPCR, well-represented for liraglutide/semaglutide peptide ligands. But these are peptides, not small molecules — MAMMAL's small-molecule DTI head may handle them suboptimally. Expect noisy predictions.

### 4. Critical methodology evaluation

**Training-data bias.** MAMMAL's substrate is target-binding labels (BindingDB Kd/IC50/Ki), antibody-antigen binding, single-cell expression. None of this encodes behavioral phenotype. Every "cognitive enhancement" inference Pierce makes through MAMMAL is necessarily *target-mediated* — it's a triangle: compound → target affinity (MAMMAL's strength) → cognitive effect (your mapping, not MAMMAL's). **The model has no opinion about whether α7 nicotinic stimulation improves working memory.** It just tells you whether your compound binds CHRNA7.

**Absent cognitive labels.** Even in CELLxGENE, the gene expression data are tissue/disease-stratified, not cognition-stratified. The closest you can get is to query "what compound shifts hippocampal pyramidal neuron expression profile X toward profile Y" — which assumes (i) the relevant cell type is in CELLxGENE, (ii) the "young/learning-enhanced" reference profile is well-defined, both of which are non-trivial.

**Polypharmacology handling.** MAMMAL scores one (compound, protein) pair at a time. It produces no integrated polypharmacology summary. For cognition, polypharmacology is the rule (modafinil = DAT + adenosine + orexin indirect + histamine indirect). You can manually aggregate per-target scores, but MAMMAL's pretraining doesn't include any joint reasoning over a target set.

**Healthy-vs-disease translation.** This is field-level, not MAMMAL-specific. Drugs that work in disease may show ceiling effects in healthy populations (cf. galantamine: dramatic in AD, modest at best in healthy elderly subgroup). The Roberts 2020 meta-analyses show even the "best" healthy enhancers (MPH, modafinil) have SMD < 0.5 across most domains. MAMMAL has no notion of this — its predictions are about ligand binding, not therapeutic context.

**Out-of-distribution failures to expect.**
- **Allosteric modulators** (α7 PAM, PDE4D NAM, ampakines) — under-represented in BindingDB; expect score deflation.
- **Peptide drugs hitting peptide-binding pockets** (GLP-1, orexin agonists) — handled poorly by small-molecule DTI.
- **Metabolic/mitochondrial agents** (methylene blue, MitoQ) — no clean target → DTI head will be uninformative; use cell-genes head.
- **Novel mechanism, no training analog** (ISR pathway, eIF2B oligomer stabilization) — MAMMAL will both underrank true positives and potentially hallucinate spurious targets.
- **Multi-target small molecules where the cognitive mechanism is the minor binding** (clemastine — H1 binding dominates the chemistry, oligodendrocyte M1/M3 effect is the cognitive piece) — DTI ranking will surface the wrong target.

**Confidence calibration.** MAMMAL outputs continuous regression values without uncertainty estimates. There is no published calibration curve, no held-out reliability diagram. **Treat predictions as rank-ordering signals, not as Kd estimates.** If you need uncertainty, ensemble across multiple seeds / fine-tunes, or use Monte Carlo dropout (not native).

**Empirical validation pipeline you'd actually use:**
1. **In silico** → MAMMAL DTI + BBBP + ClinTox + LINCS L1000 signature triangulation.
2. **Biochemical** → ChEMBL bioactivity confirmation for top hits; for novel predictions, a single-point radioligand binding assay at Eurofins/Cerep ($500–$2000 per compound-target pair) is cheap and decisive.
3. **iPSC neuron / LTP** → hippocampal slice LTP induction assays, or iPSC-derived cortical neuron MEA (multi-electrode array) firing patterns under compound exposure.
4. **Animal** → 5-CSRTT (5-choice serial reaction time, attention/processing speed), Delayed Non-Match to Sample (DNMS, working memory), Morris water maze (spatial learning), Barnes maze. Use young healthy rodents specifically, not aged or disease-model.
5. **Human RCT** → small (n=20-40) crossover, within-subject, healthy non-sleep-deprived adults, Cambridge Cognition CANTAB battery (PAL, SWM, RVP) or NIH Toolbox Cognition Battery. This is where most candidates die.

**Where MAMMAL adds value vs. alternatives.**
- **vs. CMap/L1000 alone:** MAMMAL provides target-resolution that L1000 lacks; L1000 gives you global transcriptomic response that MAMMAL's DTI head lacks. They are complementary.
- **vs. classic docking (AutoDock, Glide):** MAMMAL is sequence-based, doesn't need a crystal structure, ~100-1000x faster, but loses spatial reasoning. For targets without good structures, MAMMAL wins; for known crystallized pockets, docking is more interpretable.
- **vs. phenotypic screening:** Phenotypic screens *would* directly measure your endpoint (e.g., LTP in slices) but cost orders of magnitude more. MAMMAL is the in silico prefilter to phenotypic validation.
- **vs. specialized DTI models (DeepPurpose, etc.):** MAMMAL is competitive on benchmarks (NRMSE 0.906 on PEER/BindingDB vs PEER's 0.942) and has unified prompt syntax + multimodal capability that specialized DTI models lack.

**Failure modes / risk register.**
1. *Hallucinated affinity at understudied targets.* MAMMAL extrapolates; predicted high pKd at a target with little training coverage may be unreliable. Sanity check via Tanimoto similarity to training compounds — IBM did this in their carfilzomib validation.
2. *Off-target liabilities not flagged.* MAMMAL's ClinTox head is a coarse signal; it will not flag hERG, CYP3A4 interaction, idiosyncratic hepatotox. Add a separate ADMET pass (e.g., ADMET-AI, SwissADME).
3. *Allosteric blindness* as noted above.
4. *No PK/PD model.* MAMMAL doesn't know about brain penetration kinetics, half-life, plasma protein binding. BBBP gives a binary in/out, not exposure.

### 5. Concrete deliverables for Pierce

#### A turnkey workflow

```
project/
├── 01_setup/
│   └── conda env (mammal_env), pip install, weight download verification
├── 02_targets/
│   ├── targets.csv  # 18 panel proteins above, UniProt + AA sequence
│   └── fetch_aa.py  # UniProt REST → AA sequences
├── 03_compounds/
│   ├── compounds.csv  # ~300 compounds: SMILES, name, mechanism class, RCT evidence flag
│   └── (DrugBank XML parser if you have a license; otherwise PubChem REST)
├── 04_score/
│   ├── score_dti.py  # 18 × 300 = 5400 (target, drug) pairs, batched
│   ├── score_bbbp.py  # 300 compounds
│   ├── score_clintox.py  # 300 compounds × 2 heads
│   └── results.parquet
├── 05_triangulate/
│   ├── lincs_l1000.py  # pull signatures from clue.io API, compute reversal score
│   ├── opentargets.py  # cross-reference predicted (compound, target) pairs
│   └── ranked_candidates.csv
└── 06_followup/
    └── (handoff to wet-lab partner for radioligand binding or LTP slice work)
```

#### Models / databases to triangulate against (brief)

- **LINCS L1000 / CMap (clue.io)** — gold standard transcriptomic perturbagen library; complementary to MAMMAL's target view.
- **OpenTargets Platform** — target-disease evidence; cross-reference predicted (compound, target) hits for prior literature.
- **ChEMBL bioactivity** — ground truth for sanity-checking DTI predictions.
- **BBB Predictor / B3DB** — independent BBB permeability for cross-check.
- **PubChem BioAssay** — high-throughput screening results for many CNS targets.
- *(per scoping constraint: not surveying competing foundation models)*

#### Open research questions / gaps where novel work is publishable

1. **Cognition-stratified L1000 signatures.** No public reference signature for "young hippocampus" or "LTP-induced cortical neuron." Building one (from BrainSpan + Allen Brain + relevant published RNA-seq) and using MAMMAL's cell-genes head against it would be a contribution.
2. **MAMMAL + ISR co-fine-tune.** Fine-tune MAMMAL on the eIF2B-ISRIB-analog literature (PDB structures, ~20-30 compounds) to give it a representation of an under-druggable target — proof-of-concept for using MAMMAL as a few-shot learner on a novel mechanism.
3. **Allosteric awareness benchmark.** No benchmark exists for "does this DTI model correctly handle allosteric PAMs vs orthosteric agonists." Building such a benchmark using α7 nAChR (PAM literature) and PDE4D (NAM literature) and evaluating MAMMAL would be valuable.
4. **MAMMAL polypharmacology prompt.** The current prompt syntax is one-target-at-a-time. Designing a multi-target prompt (a panel "fingerprint" embedding) and validating it on known polypharmacologic enhancers (modafinil, galantamine) would extend the model's expressiveness.
5. **Prospective wet-lab validation of cognition-relevant MAMMAL predictions.** The Carfilzomib repurposing demonstration is the template — pick one MAMMAL-flagged compound-target pair with no prior literature, run a binding assay + (if positive) a primary neuron LTP assay. Publication-worthy regardless of outcome.

## Recommendations

**Stage 1 — Week 1-2.** Set up the environment, run the Colab tutorial end-to-end, then reproduce one published MAMMAL result (DTI on BindingDB pKd or BBBP) locally to confirm your 5070 environment is working and you trust the outputs.

**Stage 2 — Week 3-4.** Build `targets.csv` (the 18-target panel) and `compounds.csv` (300 compounds). Score all pairs. Sanity-check positive controls (donepezil/AChE, methylphenidate/DAT, pitolisant/HRH3) score in the top 5% of their respective target columns. **Decision gate:** if positive controls fail to rank in the top 5%, stop — your prompt formatting or model setup is wrong.

**Stage 3 — Week 5-6.** Triangulate top-50 candidates against ChEMBL (known evidence?), OpenTargets (literature?), LINCS L1000 (transcriptomic reversal?), and BBBP/ClinTox MAMMAL heads. Produce a final ranked candidate list with explicit columns: MAMMAL pKd, target panel overlap count, BBBP score, ClinTox score, prior healthy-adult RCT evidence (Y/N/unclear), suggested follow-up.

**Stage 4 — Optional Month 2-3.** Pick 2-3 candidates where MAMMAL surfaces a non-obvious (compound, target) prediction (e.g., a known drug at a cognition-relevant target with no prior literature) and outsource a radioligand binding assay (Eurofins/Cerep, $500-2000 each) for empirical confirmation.

**Stage 5 — Month 3+.** Pursue one of the open research questions above as a publishable mini-project, ideally the "MAMMAL allosteric awareness benchmark" (lowest experimental cost, highest reusability).

**Benchmarks/thresholds that change the recommendation:**
- *If MAMMAL's positive control ranking fails:* switch to a docking/CMap approach; MAMMAL isn't the right tool.
- *If your top hits are all already-approved CNS drugs with known healthy-adult RCT data:* MAMMAL is recapitulating literature, not generating insight — consider fine-tuning on a more targeted dataset (e.g., specific to ampakines, or to allosteric PAMs) to push it into less-explored chemistry.
- *If you find a novel high-pKd prediction at an under-studied cognition target:* radioligand-assay-confirm before any downstream investment.
- *If wet-lab confirmation triple-fails on top predictions:* MAMMAL's calibration is worse than the paper suggests for cognition-relevant targets, and the right move is to incorporate a calibration step (Platt scaling against ChEMBL ground truth).

## Caveats

1. **Speculative inference about MAMMAL on cognition.** The MAMMAL paper, supplementary materials, IBM Research blog posts, and the team's published collaborations (cancer cell-line response, influenza HA antibodies) **do not contain a single cognition-focused application or prediction.** Everything in Section 3 above is my mapping of MAMMAL's published capabilities onto the cognition problem — it is not a recapitulation of published MAMMAL output. Where I write "MAMMAL handling: …", I am predicting based on the model's training distribution and architecture, not citing a published result.

2. **Layer-count / architectural specifics.** Hidden dimension is 768 (paper explicitly states). Exact encoder depth, decoder depth, attention heads, FFN dim, vocab size, and precise FLOPs are **not published**. Pierce will need to load `config.json` from the HF safetensors to confirm.

3. **VRAM extrapolation.** "Runs on a 12 GB 5070" is inferred from parameter count and IBM's confirmed Colab T4 (16 GB) compatibility, not from a measured benchmark on a 5070. Worst case for full-parameter fine-tuning with sequence length 2048 may push against 12 GB — fall back to LoRA / gradient checkpointing if so.

4. **Healthy-adult RCT effect sizes throughout Section 3 are small.** Methylphenidate's recall SMD = 0.43 is one of the field's strongest signals — for context, this is "small-to-moderate" by Cohen's conventions. Any "novel cognitive enhancer" surfaced by MAMMAL in healthy adults should be expected to produce effects in the same SMD = 0.1-0.4 range, not the dramatic effects of Alzheimer's drug trials in disease populations. The field's ceiling appears low.

5. **GLP-1 narrative shift.** Initial observational data (Xu et al., target trial emulation) suggested 40-70% reduced AD diagnosis with semaglutide; the prospective Phase 3 EVOKE/EVOKE+ trials read out at CTAD 2025 as **negative on primary cognitive endpoint** in early symptomatic AD. The healthy-adult OxSENSE study has not yet read out. Treat the GLP-1-cognition story as currently unsettled.

6. **Allosteric and PAM bias.** Many of the most interesting cognitive-enhancement mechanisms (α7 PAM, PDE4D NAM, AMPA PAM, eIF2B stabilization) involve allosteric or interface-disrupting binding sites that are under-represented in BindingDB. MAMMAL's DTI head will systematically underweight these mechanisms relative to their actual biological importance. **This is the single most important methodological caveat for this use case.**

7. **arXiv vs. journal version discrepancy.** The arXiv v3 (May 2025) reports MAMMAL beats AlphaFold3 on antibody-antigen binding for "3 of 4 targets"; the npj Drug Discovery journal version (May 2026) expanded this to 7 targets with MAMMAL better on 5 of 7. Cite the journal version for any external use.

8. **No regulatory or therapeutic guidance is implied here.** This report is a methodology and engineering analysis for an in silico computational workflow. None of the compounds listed are recommended for self-administration; many are unapproved or unavailable.