# Diagnosing MAMMAL DTI Anti-Correlation at SLC6A3 / SLC6A2 / GRIN2A / GRIN2B: A 1–2 Day Diagnostic Protocol

## TL;DR
- The four inverted targets fail in two qualitatively different regimes — DAT/NET because BindingDB's monoamine-transporter coverage collapses to tropane/phenethylamine scaffolds (Vaughan/Newman RTI series; Carroll-school cocaine analogs), and GRIN2A/2B because their pharmacology is dominated by allosteric phenylethanolamine NAMs (ifenprodil class) bound at an ATD dimer interface that MAMMAL's single-chain protein-sequence tokens never touch. Diagnostic D (Tanimoto-vs-prediction inversion) plus Diagnostic E (Boltz-2 positive controls) are the two highest-information cheap tests; run those first.
- Pierce's three hypotheses are not mutually exclusive. The literature predicts SLC6A3/SLC6A2 will route to Scenario 2 (rank-resolution loss within a saturated tropane cluster — §7.6 LoRA) and GRIN2A/GRIN2B will route to Scenario 3 (representational gap, allosteric site — deprecate + Boltz-only). Pre-commit thresholds: scaffold overlap >60% AND K-S shift <0.2 → rank-resolution; scaffold overlap <25% AND K-S shift >0.5 → manifold mismatch; ESM2 pocket-attention z-score <0 vs random-residue null → representational gap.
- Statistical power is the blocking concern Pierce hasn't called out: ρ = −0.71 at n=26 has a 95% Fisher-z CI of roughly [−0.86, −0.45]; ρ = −0.30 at n=14 has CI [−0.71, +0.27], and ρ = −0.35 at n=8 has CI [−0.83, +0.43]. **GRIN2A and GRIN2B are not statistically distinguishable from zero with current n.** Before spending a week on LoRA, run the diagnostics on DAT/NET (where the signal is real) and treat GRIN2A/2B as "insufficient evidence, expand n first."

---

## Section 1 — Literature on DTI Foundation-Model Failure Modes

### 1.1 What MAMMAL actually is, and what it was actually evaluated on

MAMMAL (Shoshan et al., arXiv:2410.22367; npj Drug Discovery 2026, Nature) is a 458M-parameter T5-style encoder–decoder pretrained on **2 billion biological samples** across proteins, small molecules, and single-cell modalities, with a modular tokenizer that gives amino-acid residues, SMILES atoms, and gene names disjoint token IDs in a unified vocabulary. The DTI head (`ibm-research/biomed.omics.bl.sm.ma-ted-458m.dti_bindingdb_pkd`) is a regression head producing scalar pKd. Inputs are normalized with `norm_y_mean=5.79384684128215, norm_y_std=1.33808027428196` — i.e., the BindingDB training pKd distribution centers around 5.79 with SD 1.34. **Any predictions clustering tightly around 5.79 are the model returning its prior, not learning a target-specific signal.** This is the first sanity check Pierce should run before any diagnostics.

The npj Drug Discovery paper reports the DTI BindingDB pKd benchmark using **NRMSE = 0.906 ± 0.011**, against PEER-benchmark SOTA NRMSE = 0.942 ± 0.028 — a 3.8% improvement. Verbatim: *"Our model achieves an average NRMSE of 0.906 (Table 1), demonstrating a solid improvement of 3.8% over the SOTA."* NRMSE near 0.9 means the model explains only ~18% of the variance (1 − 0.9²) on the held-out class — weak in absolute terms. Note also that the paper reports **only NRMSE**: no Pearson r, no RMSE, no MAE, no Spearman ρ. The benchmark uses PEER's split, which **holds out four entire protein classes** (estrogen receptor, GPCRs, ion channels, RTKs) as the test set. SLC6 transporters (ion channels family in some taxonomies) and ionotropic glutamate receptors are exactly the held-out classes — meaning the benchmark numbers do not tell us anything about in-distribution interpolation behaviour at the targets Pierce cares about.

The HuggingFace fine-tuned head uses the TDC `Drug+Target` cold-split instead — distinct from the PEER 4-class holdout used in the paper. **Neither evaluation regime mirrors Pierce's use case**: he's asking the model to rank 25–26 compounds against a single SLC6/GRIN target where similar binders demonstrably exist in training. That's not a cold-split scenario — that's an interpolation scenario where the canonical failure mode is rank-resolution loss within a known scaffold cluster.

### 1.2 Imbalanced training / scaffold-saturation failures

**Sundar & Colwell, J. Chem. Inf. Model. 60(1):56–62, 2020 (online Dec 11, 2019), doi:10.1021/acs.jcim.9b00415** — "The Effect of Debiasing Protein–Ligand Binding Data on Generalization" — showed that aggressive debias splits (asymmetric validation embedding, Wallach & Heifets AVE) can *destroy* generalizability rather than improve it, because they remove the very compounds the model needs to interpolate. Lesson: when BindingDB coverage at a target is dense around one chemotype, models learn "scaffold present → high pKd" and lose the within-scaffold gradient. This is the canonical scaffold-saturation regime.

**Zhang, Hu, Jiang, Chen, Xu, Zhang 2025 ("Rethinking the generalization of drug target affinity prediction algorithms via similarity aware evaluation," ICLR 2025 Oral, arXiv:2504.09481)** — demonstrates that *"the canonical randomized split of a test set in conventional evaluation leaves the test set dominated by samples with high similarity to the training set. The performance of models is severely degraded on samples with lower similarity to the training set but the drawback is highly overlooked."* This directly predicts that compounds in our 298-library whose nearest BindingDB neighbour has Tanimoto < ~0.4 will be predicted near the prior mean (5.79), and that the rank order over those compounds will be effectively random with respect to true pKd. Combined with the fact that several true high-pKd binders in ChEMBL may be the *most* dissimilar to the BindingDB cocaine-analog cluster, you get anti-correlation by construction.

**Graber, Stockinger, Meyer, Mishra, Horn, Buller, Nat. Mach. Intell. 7:1713–1725 (Oct 21, 2025), doi:10.1038/s42256-025-01124-5** ("Resolving data bias improves generalization in binding affinity prediction," PDBbind CleanSplit) — verbatim: *"Retraining current top-performing models on CleanSplit caused their benchmark performance to drop substantially, indicating that the performance of existing models is largely driven by data leakage."* Confirms that benchmark numbers like NRMSE 0.906 are inflated by within-cluster nearest-neighbour effects.

**FS-Mol (Stanley et al. NeurIPS 2021)** — established that QSAR tasks with a **median of 94 compounds per task** ("the median number of compounds per task is 94, far below alternative datasets") are the regime where pretrained backbones routinely fail to beat a Random Forest on Morgan fingerprints. Pierce's n=8–26 per inverted target is far below that floor: a frozen MAMMAL head with no per-target adaptation is operating below the support where its representations were validated.

**Dablander et al. 2023 (J. Cheminform., "Exploring QSAR models for activity-cliff prediction")** — used dopamine receptor D2 as one of three case studies; showed that *"QSAR models frequently fail to predict ACs"* and that AC-sensitivity collapses when both compounds in a matched pair are unseen. For DAT/NET this is the dominant noise channel because tropane SAR is famously cliff-rich (Carroll RTI series literature; the RTI-31/RTI-55/RTI-32 sweep over X = Cl, I, Me on the phenyl gives DAT pKd differences of ~0.4–0.5 log unit per atom swap).

### 1.3 Tokenisation / sequence-representation failures

MAMMAL uses a per-residue AA tokenizer (`<@TOKENIZER-TYPE=AA>`) inherited from the modular tokenizer (`fuse.data.tokenizers.modular_tokenizer`); each natural amino acid maps to a single token, so there is no BPE-style subword fragmentation across residues. **There is no truncation artifact at the substrate pocket per se** — the AA tokenizer is faithful at the residue level. What CAN fail is whether MAMMAL's encoder attention is concentrated on the pocket residues that determine ligand discrimination, vs distributed across the whole 600-residue SLC6 protein.

The interpretability literature is converging on this:
- **Lin et al. 2023 (ESM-2/ESMFold, Science)** — per-residue attention in ESM-2 layers 25–32 is highly concentrated on contact-forming residues; coevolutionary "categorical Jacobian" signal (Zhang et al. 2024) cleanly recovers known contacts.
- **Luo et al. 2024 (Brief. Bioinform. 25(2):bbad534)** — *"Interpretable feature extraction and dimensionality reduction in ESM2 for protein localization prediction"* — confirms that the CLS-token representation captures the functional features that downstream heads predict from.
- **ESM2_AMP (2025, BMC Bioinformatics)** explicitly showed *"strong correlations between segments with high attention weights and known functional regions of amino acid sequences."*

This means: ESM2-650M attention extracted at the known substrate residues *should* be measurably elevated at DAT/NET — if it isn't, that's a representational gap. **The hDAT cryo-EM (Srivastava et al., Nature 632:672–677, 2024, doi:10.1038/s41586-024-07739-9) gives exact residues to test**: subsite A — F76, A77, V78, D79, A81, F320, G323; subsite B — V152, G153, Y156, F326, S422, A423, G426; gate H-bond pair Y156–D79. Verbatim from the paper: *"residues F76, A77, D79, A81, F320 and G323 from subsite A frequently interact with the tropane moiety, and residues V152, G153, Y156, F326, S422, A423 and G426 interact with the fluorophenyl moiety."*

For hNET (Yuan et al., Cell Research 2024, doi:10.1038/s41422-024-01024-0, PDB 8ZP1): F72, D75, A77, N78 (subsite A); Y152, N153, G423 (subsite B); V148, F323, F329, M424 (subsite C); gate Y152–D75. Verbatim: *"They are observed in the same binding pocket as NE, mainly through hydrophobic interactions with F72 (SA), Y152/G423 (SB), and V148/F323 (SC). In this outward-open state, Y152 (SB) forms an H-bond with the side chain of D75 (SA)... The secondary amine group of atomoxetine and the morpholine nitrogen of reboxetine form an H-bond with the carbonyl oxygen of F72."*

For GRIN2B the situation is qualitatively different. Karakas et al. 2011 (Nature) and the Mony / Stroebel / Paoletti reviews establish that **the dominant pharmacology in ChEMBL CHEMBL1904 is the ifenprodil-class allosteric site, which sits at the GluN1/GluN2B amino-terminal-domain DIMER INTERFACE — not on the GluN2B monomer at all**. Verbatim from Mony et al. (the GluN2B-selective NAM mapping paper, Mol. Pharmacol. 2012): *"Ifenprodil and related compounds inhibit receptors incompletely (e.g., 90%) at saturating concentrations… ifenprodil makes direct interactions with the UL-UL mainly through hydrophobic interactions."* MAMMAL is fed a single GluN2B sequence. There is no token that represents "the GluN1-facing interface of GluN2B ATD when dimerized with GluN1." A T5 encoder operating on a single chain physically cannot learn ifenprodil-site SAR. This isn't a tokenization bug; it's a representational impossibility within the inputs.

### 1.4 Why a regression model goes ANTI-correlated

For a model to be anti-correlated rather than uncorrelated, *something* in the input must actively predict the wrong direction. Three published mechanisms:

1. **Label-encoding sign flips during fine-tuning.** If a non-trivial fraction of BindingDB rows at SLC6A3 misencoded Ki vs IC50 vs Kd (BindingDB does not enforce assay-type consistency at the row level), then `harmonize_affinities(mode='max_affinity')` (the published MAMMAL preprocessing) will systematically pull *inhibitor* IC50 values into the same scalar bucket as *substrate* Kd values. For a transporter, substrate Km is typically larger (worse-binding) than tightest inhibitor Ki — the harmonization can invert the ordering for transporter-substrate pairs.

2. **Activity-cliff inversion** (Dablander 2023). When the test compound is on the wrong side of a cliff but nearest-neighbour interpolation picks up the cliff's "low" partner, the prediction lands below the mean while truth is above it.

3. **Dataset poisoning / silent leakage** — Graber et al. 2025 demonstrated in PDBbind that removing leakage didn't help unless you also fixed underlying label errors. If BindingDB labels for some DAT compounds were entered as `1/Kd` rather than `Kd` (a common LIMS-to-paper transcription error), the model learns the wrong sign on those rows and produces anti-correlated predictions on chemically related compounds.

The most informative single experiment that distinguishes (1)/(3) from manifold-mismatch is **Diagnostic D below**.

### 1.5 Statistical power: what ρ = −0.71 at n=26 actually means

Using the Fisher-z approximation (Bonett & Wright 2000, *Psychometrika* 65(1):23–28) for Spearman ρ:

| Target | n | ρ | 95% CI (Fisher-z) | Distinguishable from 0? | Distinguishable from ρ = −0.71? |
|---|---|---|---|---|---|
| SLC6A3 | 26 | −0.71 | [−0.86, −0.45] | yes (p < 0.001) | reference |
| SLC6A2 | 25 | −0.53 | [−0.76, −0.18] | yes (p ≈ 0.007) | no (CIs overlap) |
| GRIN2B | 14 | −0.30 | [−0.71, +0.27] | **no** | no (CIs overlap) |
| GRIN2A | 8 | −0.35 | [−0.83, +0.43] | **no** | no (CIs overlap) |

**The DAT and NET inversions are statistically real. The GRIN2A and GRIN2B inversions are not.** Pierce should not spend a week on LoRA for GRIN2A based on 8 compounds. The right action for GRIN2A/2B is to expand n first (cheap: pull more compounds from the cognition-relevant ChEMBL set into the 298-library) and re-measure.

---

## Section 2 — The Five Diagnostics

All diagnostics live under `src/mammal_repurposing/diagnostics/`. Wall-clock estimates assume RTX 5070 (12 GB VRAM), WSL2 Ubuntu, with cached ESM2-650M embeddings and the existing DTI grid parquet already on disk.

### Diagnostic A — Scaffold saturation analysis (Murcko clustering)

```python
# diagnostics/scaffold_saturation.py
from pathlib import Path
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, DataStructs

def murcko_scaffold_overlap(
    target_uniprot: str,
    library_compounds_parquet: Path,
    chembl_db_path: Path,
    bindingdb_tsv: Path | None = None,
    tanimoto_threshold: float = 0.4,
    pchembl_active_threshold: float = 6.0,
) -> dict:
    """
    Returns:
        n_library_compounds, n_chembl_compounds, n_bindingdb_compounds,
        library_scaffold_clusters (list of (generic_murcko_smiles, count)),
        bindingdb_top_cluster_smiles,
        library_in_bindingdb_top_cluster_pct,
        decision: 'manifold_mismatch' | 'rank_resolution' | 'ambiguous'
    """
```

**Protocol (per target):**
1. Pull the 25–26 SMILES from `dti_grid.parquet` filtered on `target == uniprot`.
2. Pull all ChEMBL canonical SMILES joined on activities with `target_dictionary.chembl_id == 'CHEMBL238'` (DAT), `'CHEMBL222'` (NET), `'CHEMBL1972'` (GRIN2A), `'CHEMBL1904'` (GRIN2B) and `pchembl_value IS NOT NULL`.
3. Pull BindingDB rows for the same UniProt via `tdc.multi_pred.DTI(name='BindingDB_Kd').get_data()` filtered on `Target_ID`.
4. For each SMILES: `s = MurckoScaffold.GetScaffoldForMol(mol)`; then `g = MurckoScaffold.MakeScaffoldGeneric(s)`; canonicalize, then run `GetScaffoldForMol` once more (per the RDKit Discussion #6844 — this matches the original Bemis–Murcko 1996 definition; otherwise exocyclic atoms persist).
5. Cluster on Morgan-2 fingerprints (ECFP4, 2048 bits) with single-linkage Tanimoto > 0.4. This matches the widely used scaffold-cluster definition (e.g., the ChemDiv workflow, RDKit MaxMin).
6. Compute `library_in_bindingdb_top_cluster_pct` = fraction of library SMILES whose generic-Murcko scaffold matches the BindingDB top-cluster scaffold OR has Tanimoto ≥ 0.4 to its centroid.

**Decision thresholds (literature-anchored, refined from Pierce's guess):**

| `library_in_bindingdb_top_cluster_pct` | Route | Justification |
|---|---|---|
| < 25% | manifold mismatch | Zhang et al. 2025 (ICLR Oral) show severe degradation begins at Tanimoto-to-training < 0.4 |
| 25–60% | ambiguous → check Diagnostic D | mixed regime |
| > 60% | rank-resolution loss in saturated cluster | Dablander 2023 AC-failure regime |

Pierce's 30%/70% guess is in the right ballpark; literature pulls the upper bound down to ~60% because activity-cliff failure manifests well before saturation reaches 70%.

**Wall-clock:** ~3–8 min per target on CPU. RDKit Murcko + Morgan FP for 3,000–4,000 ChEMBL compounds and 1,000–10,000 BindingDB compounds is I/O-dominated.

**Positive control: DRD1 (CHEMBL2056) and HCRTR1 (CHEMBL5113).** Expect `library_in_bindingdb_top_cluster_pct` ≈ 30–60%. The control checks that we observe an INTERMEDIATE saturation at successful targets — not zero, not 100%.

**Negative control: pick one MAMMAL_ONLY_WEAK target with ρ near 0 and n ≥ 15.** Expect saturation in the 30–60% intermediate range with no clear routing signal — if your diagnostic confidently routes the negative control to a scenario, your thresholds are over-calling.

### Diagnostic B — Pchembl distribution comparison (K-S + Wasserstein)

```python
# diagnostics/distribution_shift.py
from scipy.stats import ks_2samp, wasserstein_distance

def pchembl_distribution_shift(
    target_uniprot: str,
    library_pred_parquet: Path,
    chembl_db_path: Path,
    bindingdb_tsv: Path,
) -> dict:
    """
    Returns: ks_stat_lib_vs_chembl, ks_pvalue_lib_vs_chembl,
             ks_stat_lib_vs_bindingdb, wasserstein_lib_vs_chembl,
             wasserstein_lib_vs_bindingdb,
             mean_pchembl_library, mean_pchembl_chembl, mean_pkd_bindingdb,
             decision: 'panel_revision' | 'scaffold_aware_AL' | 'exclude'
    """
```

**Protocol (per target):**
1. Distribution L: ChEMBL pchembl_value for the 25–26 library compounds AT target T.
2. Distribution C: all ChEMBL pchembl_value for target T (the 3,391 / 3,635 / 3,228 / 239 records).
3. Distribution B: BindingDB pKd for target T (via TDC).
4. Run `ks_2samp(L, C)` and `ks_2samp(L, B)`. Compute Wasserstein distances.

**Decision thresholds:**

| K-S statistic (L vs C) | Route |
|---|---|
| > 0.5 | §7.3 panel revision — library is not a representative ChEMBL sample at this target |
| 0.2 – 0.5 | §7.13 scaffold-aware active learning |
| < 0.2 | exclude distribution mismatch as a cause |

A K-S of 0.5 corresponds roughly to "the median library pchembl falls outside the inter-quartile range of the ChEMBL distribution" — a useful operational definition of "the panel does not represent the target." Note: SIMPD (Landrum et al., J. Cheminform. 15:119, 2023) provides validated machinery for generating simulated time-splits using a multi-objective GA over Novartis project data, but does not itself prescribe K-S = 0.5 as a canonical cutoff — the 0.5 threshold here is engineering judgment, not literature precedent.

**Wall-clock:** ~2 min per target. Pure pandas + scipy.

**Positive control: DRD1 should give K-S < 0.2** (library was constructed to match ChEMBL distribution at the strong targets). If it doesn't, the panel has a sampling bug and that must be fixed before any of the other diagnostics' conclusions are valid.

### Diagnostic C — ESM2 attention at substrate-binding residues

```python
# diagnostics/esm2_attention.py
import torch
from transformers import AutoTokenizer, EsmModel

def attention_at_pocket(
    target_uniprot: str,
    pocket_residues_human_numbering: list[int],
    esm2_embeddings_cache: Path,
    layer_range: tuple[int, int] = (15, 30),
) -> dict:
    """
    Returns: per_layer_pocket_attention_mass,
             mean_pocket_attention,
             null_distribution_mean,
             null_distribution_std,
             z_score_vs_null,
             decision: 'pocket_seen' | 'pocket_invisible'
    """
```

**Protocol:**
1. Load ESM2-650M with `output_attentions=True`. If only embeddings were cached, re-run forward pass for the 7 proteins (one-off ~30 s GPU operation per ~600-residue protein on the 5070).
2. For each layer in 15–30 (the functional-site-rich middle-late band per the ESM-2 interpretability literature), compute the column-sum attention received at each residue, averaged across heads.
3. Sum column-sum attention over the target's pocket residues = `pocket_attention_mass`.
4. Build a null: 1,000 random residue sets of the same cardinality, compute mass each time. Compute z-score of true pocket mass vs null.
5. Repeat for DRD1, HCRTR1 (positive controls), and one MAMMAL_ONLY_WEAK target (negative control).

**Residue lists (human numbering throughout):**

- **SLC6A3 (hDAT, UniProt Q01959)**: Subsite A — F76, A77, V78, D79, A81, F320, G323. Subsite B — V152, G153, Y156, F326, S422, A423, G426. Gate H-bond pair — Y156, D79. *Source: Srivastava et al. Nature 632:672–677, 2024 — "the tropane moiety of β-CFT is positioned toward subsite A, facing D79 and A81 on TM1b, F76 on TM1a and G323 on TM6… V152, S422 and Y156 participate in van der Waals contacts, and F326 on the TM6a–TM6b linker forms an edge-to-face contact with the phenyl ring of the fluorophenyl group."*
- **SLC6A2 (hNET, UniProt P23975)**: Subsite A — F72, D75, A77, N78. Subsite B — Y152, N153, G423. Subsite C — V148, F323, F329, M424. *Source: Yuan et al. Cell Research 2024, PDB 8ZP1.* Cross-validated with mutagenesis (Jha, Ragnarsson & Lewis 2020, Front. Pharmacol.): *"the functional binding pocket for NE comprised residues A73, A77, N78, V148, N153, I156, G320, F329, N350, S420, G423, and M424."*
- **GRIN2B (UniProt Q13224)**: Ifenprodil-site residues on GluN2B-ATD lobe — P78, F114, Q110, F176, I150, D101, E236, M207 (Karakas et al. Nature 2011; rat numbering — apply +6 offset for human numbering check). **NOTE: this is a dimer-interface site; expect low pocket attention even in a well-trained model because the interface partner (GluN1) isn't in the input.** This is the *intended* finding for Scenario 3.
- **GRIN2A (UniProt Q12879)**: ATD residues analogous to GluN2B at the equivalent dimer interface — but GluN2A lacks the ifenprodil site; competitive/glycine-site ligands dominate. Test the glutamate ligand-binding-domain residues (S511, T513, R518, T690 in human numbering, derived from the agonist-binding cleft).

**Decision thresholds:**

| z-score vs null | Pocket-attention regime |
|---|---|
| z > 2 | model "sees" the pocket; tokenization not the cause |
| 0 < z < 2 | weak signal; ambiguous |
| z < 0 | pocket invisible; rule out LoRA for that target |

**Cross-check**: if DAT/NET show z > 2 but GRIN2A/B show z ≈ 0 or negative, that is the cleanest possible diagnostic split for the "tokenization fine for transporters, broken for NMDA allosteric" hypothesis.

**Wall-clock:** ~5 min per target.

### Diagnostic D — MAMMAL prediction vs Tanimoto to known high-affinity binder

**This is the highest-value diagnostic.** It distinguishes "model is noise" from "model is actively inverted."

```python
# diagnostics/tanimoto_correlation.py
from scipy.stats import spearmanr, pearsonr

def mammal_vs_tanimoto_to_known_actives(
    target_uniprot: str,
    library_pred_parquet: Path,
    chembl_db_path: Path,
    pchembl_active_threshold: float = 8.0,
    fp_radius: int = 2,
    fp_bits: int = 2048,
) -> dict:
    """
    Returns: spearman_mammal_vs_tanimoto, pearson_mammal_vs_tanimoto,
             n_known_actives, mean_max_tanimoto,
             decision: 'systematic_inversion' | 'pure_noise' | 'correctly_correlated'
    """
```

**Protocol:**
1. From ChEMBL, get all SMILES at target T with `pchembl_value >= 8.0` (≥10 nM). Call this set K_T.
2. For each of the 26 library compounds at target T, compute max Tanimoto (ECFP4, 2048 bits) to any compound in K_T. Call this `t_i`.
3. Correlate `t_i` with MAMMAL's predicted pKd for compound i.
4. Compare to the same correlation on DRD1/HCRTR1 (positive control — should be ρ ≈ +0.3 to +0.6) and a MAMMAL_ONLY_WEAK target (negative control — should be ρ ≈ 0).

**Decision thresholds:**

| Spearman ρ(MAMMAL_pKd, max_Tanimoto_to_K_T) | Interpretation |
|---|---|
| ρ > +0.3 | Model rewards structural similarity to known binders correctly. Anti-correlation with ChEMBL pchembl is then driven by activity cliffs WITHIN the cluster → Scenario 2 (rank-resolution loss; LoRA worth it). |
| −0.2 < ρ < +0.3 | Model has no usable signal at this target → Scenario 1 (pure manifold mismatch / model returning prior). Verify predictions cluster around `norm_y_mean = 5.79`. |
| ρ < −0.2 | **Active inversion**: model penalizes the right structural class → Scenario 4 (label-sign error in training subset, or systematic `harmonize_affinities` damage). Highest investigative priority; routes to a BindingDB row-level audit, not LoRA. |

This diagnostic isolates *absence of signal* (Scenario 1) vs *presence of perverse signal* (Scenario 4) — which route to completely different remedies. It is the diagnostic Pierce was missing the literature for.

**Wall-clock:** ~3 min per target.

### Diagnostic E — Boltz-2 sanity check on canonical binders

```python
# diagnostics/boltz_positive_control.py
def boltz_positive_controls(
    target_uniprot: str,
    canonical_binder_smiles: list[str],
    boltz_predictions_parquet: Path,
    library_pred_parquet: Path,
) -> dict:
    """
    Returns: n_controls, n_in_top25_pct, mean_affinity_probability_binary,
             mean_predicted_pIC50,
             decision: 'cluster_rescue' | 'shared_failure' | 'partial_rescue'
    """
```

**Canonical positive controls (textbook DAT/NET pharmacology):**
- **DAT**: cocaine (CHEMBL370805), methylphenidate (CHEMBL796), GBR-12909 / vanoxerine (CHEMBL12713), mazindol (CHEMBL1525), bupropion (CHEMBL894). All have published pKd > 7 at hDAT.
- **NET**: atomoxetine (CHEMBL641), reboxetine (CHEMBL471), desipramine (CHEMBL72), nisoxetine, nortriptyline (CHEMBL84). All have pKd > 7.
- **GRIN2B**: ifenprodil (CHEMBL34833), Ro 25-6981 (CHEMBL285520), EVT-101 (CHEMBL2347571), CP-101,606/traxoprodil (CHEMBL45816). Published EC50 ranking from the Pápai et al. 2018/2019 GluN2B literature: EVT-101 (22 ± 8 nM) > Ro 25-6981 (60 ± 30 nM) > ifenprodil (100 ± 40 nM) > eliprodil (1300 ± 700 nM). A model that doesn't reproduce roughly this rank order on Boltz-2 has not learned GluN2B ATD pharmacology.
- **GRIN2A**: GNE-6901 / GNE-8324 (CHEMBL3645525), TCN-201, MPX-004.

**Decision thresholds:**

| Controls in Boltz-2 top 25% of library | Routing |
|---|---|
| 4–5 of 5 | Cluster rescue: Boltz-2 sees what MAMMAL misses → §7.7 cross-DTI ensemble; no LoRA needed |
| 2–3 of 5 | Partial rescue → ensemble + targeted LoRA |
| 0–1 of 5 | Shared failure: representational gap in BOTH sequence and structure models → §7.6 LoRA on cognition corpus IS justified, or deprecate |

**Wall-clock:** Boltz-2 takes ~3–8 min per protein–ligand pair (full structure + affinity) on RTX 5070 / 12 GB at batch size 1. For 5 controls × 4 targets = 20 pairs, the overnight sweep is realistic (3–5 h total). If the overnight sweep is already running, this diagnostic is downstream parsing only — ~5 min.

---

## Section 3 — Decision Framework

Routing is a hierarchical decision tree. **Run Diagnostic D first**; it's the cheapest and most discriminating. Then A for the targets where D didn't fire. Then C and E to refine. B is a sanity check on the panel itself.

```
═══════════════════════════════════════════════════════════════════════════
Scenario 1: Manifold mismatch from training under-coverage
  Diagnostic D: ρ(MAMMAL, Tanimoto) ∈ [-0.2, +0.3]
  Diagnostic A: library_in_bindingdb_top_cluster_pct < 25%
  Diagnostic B: K-S(library_pchembl, chembl_pchembl) < 0.3
  Diagnostic E: Boltz-2 rescues 4-5/5 controls
  → ACTION: §7.7 cross-DTI ensemble with MMAtt-DTA-style baseline (2-3 days)
  → CONFIDENCE: high
  → MOST-LIKELY-TARGETS: none of the 4 inverted; this fits MAMMAL_ONLY_WEAK targets
  → JUSTIFICATION: Zhang ICLR 2025 (low-Tanimoto degradation);
                   Sundar & Colwell 2020 (debias-improved coverage)

═══════════════════════════════════════════════════════════════════════════
Scenario 2: Rank-resolution loss within saturated tropane/phenethylamine cluster
  Diagnostic D: ρ(MAMMAL, Tanimoto) > +0.3 (model correctly identifies cluster)
  Diagnostic A: library_in_bindingdb_top_cluster_pct > 60%
  Diagnostic B: K-S < 0.2
  Diagnostic C: z-score of pocket attention > 2
  Diagnostic E: Boltz-2 rescues 2-3/5 controls
  → ACTION: §7.6 LoRA fine-tune on cognition-specific corpus, rank-stratified loss
            (Dablander 2023 twin-network training or pairwise ranking loss)
  → CONFIDENCE: medium
  → MOST-LIKELY-TARGETS: SLC6A3, SLC6A2
  → JUSTIFICATION: Dablander 2023 (AC failure within saturated clusters);
                   Carroll RTI series literature (intra-tropane SAR cliffs)

═══════════════════════════════════════════════════════════════════════════
Scenario 3: Representational gap (sequence input cannot see the binding site)
  Diagnostic C: z-score of pocket attention < 0
  Diagnostic E: Boltz-2 also fails (0-1/5 controls in top 25%)
  → ACTION: Deprecate target from MAMMAL panel; rely on Boltz-2 with explicit
            dimer/heterodimer template; apply §7.9 neurobiological prior if the
            target stays in the cognition panel
  → CONFIDENCE: high
  → MOST-LIKELY-TARGETS: GRIN2B (ifenprodil site = GluN1/GluN2B INTERFACE,
                                 MAMMAL gets only GluN2B chain)
  → JUSTIFICATION: Karakas et al. Nature 2011 (ifenprodil binds dimer interface);
                   Mony et al. Mol. Pharmacol. 2012 (UL-UL interaction mapping)

═══════════════════════════════════════════════════════════════════════════
Scenario 4: Active inversion (label-sign error or harmonize_affinities damage)
  Diagnostic D: ρ(MAMMAL, Tanimoto) < -0.2  (PERVERSE)
  → ACTION: BindingDB row-level audit at target T:
            (i) IC50/Ki/Kd type confusion in rows where harmonize picked max
            (ii) substrate-vs-inhibitor encoding (substrate Km as Kd)
            (iii) inverse-Kd (1/Kd) encoding errors
            If audit fails, file an issue at BiomedSciAI/biomed-multi-alignment.
            Stopgap: FLIP MAMMAL's sign at that target.
  → CONFIDENCE: medium-high (sign-flip is empirically rescuable)
  → MOST-LIKELY-TARGETS: SLC6A3 (n=26, ρ=-0.71 is too clean to be random)
  → JUSTIFICATION: Graber et al. Nat. Mach. Intell. 2025 (label leakage drives
                   benchmark numbers; only label-fix recovers performance)

═══════════════════════════════════════════════════════════════════════════
Scenario 5: Insufficient evidence — n too small
  Bonett-Wright Fisher-z CI for ρ at observed n includes 0
  → ACTION: Expand n 5-10× via §7.13 scaffold-aware sampling from ChEMBL
            before committing engineering time. For GRIN2A: pull 50+ compounds
            from CHEMBL1972 with pchembl ≥ 7, score, re-evaluate ρ.
  → CONFIDENCE: high (literature standard for n < 15)
  → MOST-LIKELY-TARGETS: GRIN2A (n=8), GRIN2B (n=14)
  → JUSTIFICATION: Bonett & Wright Psychometrika 2000

═══════════════════════════════════════════════════════════════════════════
Scenario 6: Mixed signal (diagnostics disagree)
  Default rules:
    • If D and A disagree, trust D (Tanimoto-vs-prediction is harder to fool
      than scaffold-percentage choice).
    • If C and E disagree (attention says pocket visible but Boltz-2 also fails),
      trust E (Boltz-2 has structure access; its failure is an independent signal
      that features aren't the bottleneck).
  → ACTION: 10,000-iteration permutation CI on ρ at the inverted target
            using compounds bootstrap-resampled from MAMMAL_ONLY_WEAK targets
            matched on pchembl distribution. If ρ at inverted target is outside
            the bootstrap 95% CI, inversion is real.
```

---

## Section 4 — Implementation Specification

Repository layout under `src/mammal_repurposing/diagnostics/`:

```
diagnostics/
├── __init__.py
├── scaffold_saturation.py        # Diagnostic A
├── distribution_shift.py         # Diagnostic B
├── esm2_attention.py             # Diagnostic C
├── tanimoto_correlation.py       # Diagnostic D
├── boltz_positive_control.py     # Diagnostic E
├── decision_tree.py              # Section 3 routing logic
├── controls.py                   # DRD1, HCRTR1 + one MAMMAL_ONLY_WEAK negative control
└── run_all.py                    # CLI: runs A-E × 4 targets + 3 controls
```

**Shared expected data shapes (pandas + pyarrow):**

```python
# dti_grid.parquet schema (existing artifact)
# columns: compound_chembl_id (str), smiles (str), target_uniprot (str),
#          target_chembl_id (str), mammal_pkd_pred (float),
#          chembl_pchembl_value (float, nullable),
#          assay_chembl_id (str, nullable),
#          mechanism_of_action (str, nullable)

# Aggregate results parquet:
#   target_uniprot, target_chembl_id, n_library, observed_rho_ci_low,
#   observed_rho, observed_rho_ci_high,
#   diagnostic_A_overlap_pct, diagnostic_B_ks_stat,
#   diagnostic_C_pocket_z, diagnostic_D_tanimoto_rho,
#   diagnostic_E_controls_in_top25, decision, confidence, notes
```

**Library version pinning (validated for WSL2 Ubuntu + RTX 5070 + CUDA 12.4):**
- `rdkit==2024.09.5`
- `pyarrow==17.0.0`, `pandas==2.2.3`
- `scipy==1.14.1`
- `tdc==1.0.6`
- `transformers==4.45.0` (for ESM2 attention extraction)
- `torch==2.4.0+cu124`
- Boltz-2: MIT/Recursion `boltz` repo; ~5 GB weights; batch size 1 fits 12 GB VRAM

**Validation gates per module:**
- Each diagnostic must produce identical outputs on DRD1 across two runs (determinism check; seed numpy and torch to 42).
- Each diagnostic must complete the full panel (4 inverted + 2 positive + 1 negative = 7 targets) in under 90 min total wall-clock.
- Decision routing must be unanimous on synthetic inputs in unit tests.

**Wall-clock budget for the full 1-2 day protocol:**
- Day 1 AM: implement A, B, D (~6h coding) → run on 7 targets (~30 min compute).
- Day 1 PM: implement C (~3h coding, attention extraction from cached embeddings) → run (~20 min compute).
- Day 1 evening: launch Boltz-2 overnight on 20 control pairs (~3–5 h).
- Day 2 AM: implement E parser (~2h), run decision_tree.py to aggregate.
- Day 2 PM: lateral diagnostic from §6 (binding-mode mix for GRIN2B). Decide v4 direction.

---

## Section 5 — Validation Gates and Publication Angle

### 5.1 Internal-consistency checks per diagnostic

| Diagnostic | Internal-consistency check on DRD1 / HCRTR1 (positive controls) | What success looks like |
|---|---|---|
| A | DRD1 should show 30–60% library overlap with BindingDB top cluster | Saturation isn't 0 at successful targets |
| B | DRD1 K-S(library, ChEMBL) should be < 0.2 | Panel construction is sound |
| C | DRD1 pocket-attention z should be > 2 (orthosteric pocket well-characterized) | ESM2 attention sees known pockets |
| D | DRD1 ρ(MAMMAL, Tanimoto) should be > +0.3 (model rewards known scaffold) | Tanimoto-vs-prediction direction is meaningful |
| E | Boltz-2 should rank dopamine, SKF-38393, dihydrexidine in DRD1 top 25% | Boltz-2 baseline works |

### 5.2 Negative-control target

Pick a MAMMAL_ONLY_WEAK target with |ρ| < 0.1 and n ≥ 15. Diagnostic outputs should be AMBIGUOUS — none of Scenarios 1/2/3/4 should cleanly fire. If your diagnostics confidently route the negative control, the thresholds need raising.

### 5.3 Publication angle

If the diagnostic produces a clean per-target attribution — *e.g., "SLC6A3 and SLC6A2 = Scenario 2 (saturated tropane cluster, rank-resolution loss confirmed by D and A); GRIN2B = Scenario 3 (ifenprodil-site allosteric pharmacology invisible to single-chain protein input)"* — the paper is:

**Title**: *"Diagnosing failure modes of foundation-model DTI heads at clinically-validated cognition targets: when scaffold saturation, label-noise inversion, and quaternary-state invisibility masquerade as the same metric."*

**Venue ranking**:
1. *J. Cheminform.* (BMC) — best fit for the diagnostic-protocol framing; ~3-month review cycle; open-access.
2. *Briefings in Bioinformatics* — better impact factor, slower review.
3. arXiv preprint first (q-bio.QM). Templates: Graber et al. 2025 (PDBbind CleanSplit, Nat. Mach. Intell.); Dablander et al. 2023 (J. Cheminform.); Zhang et al. 2025 (ICLR Oral).

**The publication-worthy contribution is the protocol + the 4-target attribution table.** The downstream LoRA/ensemble work is engineering, not the paper.

### 5.4 Risks and failure modes of the protocol itself

1. **Cached ESM2 embeddings may not include attention weights.** If only final-layer hidden states cached, re-run with `output_attentions=True` for the 7 proteins — ~3 min per protein on the 5070; not blocking.
2. **BindingDB TSV is ~3 GB**; use TDC wrapper (~70k-row Kd subset) unless patent-year stratification needed.
3. **`harmonize_affinities('max_affinity')` is the published MAMMAL preprocessing** — replicate it exactly when comparing distributions, otherwise Diagnostic B is apples-to-oranges.
4. **GRIN2B residue numbering is rat-vs-human inconsistent in the literature** (Karakas 2011 used rat GluN2B; human UniProt Q13224 differs by ~6 residues in ATD). Cross-check residue numbers against human UniProt before computing attention masses.
5. **Boltz-2's affinity head was trained predominantly on kinases / nuclear receptors / GPCRs.** Transporters and NMDA receptors are out-of-distribution; a Boltz-2 failure at DAT does NOT automatically mean MAMMAL is correctly diagnosing the protein-side problem — both models may share a transporter-coverage gap. Note this caveat in the paper.

---

## Section 6 — What Pierce hasn't asked but should have

Five lateral diagnostic angles, ordered by expected information content:

### 6.1 Binding-mode mix per target (highest priority)

For GRIN2B specifically, ChEMBL CHEMBL1904 contains a *mix* of pharmacologies:
- ifenprodil-class NAMs (ATD dimer interface)
- glycine-site competitive (LBD, on GluN1)
- glutamate-site competitive (LBD, on GluN2B)
- channel-blocking (TMD)
- intracellular C-terminal modulators

MAMMAL gets one GluN2B sequence and outputs one scalar — it cannot represent ANY of those mechanisms distinctly. If 80% of the 14 GRIN2B library compounds are ifenprodil-class, the model is being asked to predict allosteric affinity from a sequence input that doesn't contain the dimerization partner — not a model failure, a task-specification failure.

**Action**: Parse ChEMBL `mechanism_of_action` and `binding_site` fields for the 14 GRIN2B compounds. If >50% ATD allosteric, downweight or remove GRIN2B from the MAMMAL panel — replace with Boltz-2 evaluation on a heterodimeric GluN1/GluN2B template (e.g., 5IOV or 4PE5).

### 6.2 Temporal stratification of ChEMBL ground truth

ChEMBL records added pre-2010 vs post-2010 sample different SAR space:
- Pre-2010 DAT: Carroll RTI series + Kuhar lab tropanes; pKd dynamic range wide (5–10).
- Post-2015 DAT: Newman lab atypical inhibitors (benztropines, modafinil analogs); pKd range narrower (6.5–8), pharmacology shifted toward atypical binding modes.

If MAMMAL's BindingDB training set was assembled pre-2018 (likely — paper published Oct 2024), the inversion may largely be "MAMMAL learned classical DAT pharmacology; ChEMBL records driving the inversion are post-2018 atypicals."

**Action**: Split ChEMBL records at target T into pre-2015 vs post-2015 cohorts. Re-compute Spearman ρ separately. If pre-2015 ρ is high-positive and post-2015 ρ strongly negative, the failure is temporal drift, not architecture. **This is the single most likely confound and Pierce has not enumerated it.** ChEMBL 33+ added `CHEMBL_RELEASE_ID` precisely to support this analysis (Zdrazil et al., Nucleic Acids Res. 52(D1):D1180, 2024). TDC's `BindingDB_Patent` benchmark group confirms the magnitude: *"Out-of-distribution Pearson's correlation degrades from ~0.70 in-distribution to 0.42–0.43 with sophisticated DG methods (ERM, MMD, CORAL)."*

### 6.3 Family-level signature (SLC6 vs GRIN)

The 4 inverted targets cluster into 2 families. Test whether MAMMAL produces anti-correlation at OTHER family members not in the cognition panel:
- SLC6 family: SERT (SLC6A4, CHEMBL228), GAT-1 (SLC6A1, CHEMBL2074), GlyT-1 (SLC6A9, CHEMBL2068).
- GRIN family: GluN2C (CHEMBL3084), GluN2D (CHEMBL3110), GluN1 (CHEMBL1981).

If SERT shows ρ ≈ +0.5 but DAT shows ρ = −0.71 *within the same family*, the failure is target-specific not family-specific — pointing to label noise. If SERT also inverts, it's a family-level failure pointing to SLC6 representation rather than data.

### 6.4 Loss-function inspection: is MAMMAL's training loss INVERTED on a subset?

ρ = −0.71 at SLC6A3 (n=26, p < 0.001) is too large to be activity-cliff alone. There is non-trivial probability that during MAMMAL's BindingDB pretraining, a subset of SLC6 rows had pKd encoded with the wrong sign (1/Kd rather than -log10(Kd), or pKi concatenated with pKd losing the prefix). A 458M T5 will dutifully learn what it's shown.

**Diagnostic**: pull all BindingDB rows used in `dti_bindingdb_pkd` fine-tuning at SLC6A3 and SLC6A2. Check `Y ∈ [3, 12]` (sensible pKd range) and that values correlate with corresponding ChEMBL records at compound-id match. If BindingDB Y is anti-correlated with ChEMBL pchembl at the row level, the bug is in the training labels, not the model. This is the highest-impact diagnostic IF Diagnostic D returns Scenario 4. Pre-commit: if D-ρ < −0.2 at any target, file a GitHub issue on `BiomedSciAI/biomed-multi-alignment` with the row-level audit attached.

### 6.5 Statistical-power resampling (operationalize Pierce's stated concern)

For each inverted target, run 10,000 permutation tests:

```python
from scipy.stats import spearmanr
import numpy as np

def permutation_ci(pred, truth, n_iter=10000, seed=42):
    rng = np.random.default_rng(seed)
    null = np.array([spearmanr(rng.permutation(pred), truth)[0]
                     for _ in range(n_iter)])
    return np.percentile(null, [2.5, 97.5])
```

Empirically expect: SLC6A3 (n=26) significant; SLC6A2 (n=25) significant; GRIN2B (n=14) marginal; GRIN2A (n=8) NOT significant. Don't invest engineering time on hypotheses the data can't support.

---

## Section 7 — Bottom line for v4

Run Diagnostics **D and E first**; they're cheap and discriminate Scenario 2 (LoRA-worth) from Scenario 4 (training-label-bug-worth) from Scenario 5 (n-too-small).

**Most likely outcome for the four targets, with pre-commitments:**

- **SLC6A3 (DAT)**: Scenario 2 + possibly Scenario 4. Tropane saturation will be high; if Diagnostic D returns ρ > +0.3 with truth-anti-correlation, it's saturated rank-loss (LoRA worth it). If D returns ρ < −0.2, audit BindingDB labels FIRST.

- **SLC6A2 (NET)**: Scenario 2 with high confidence. NET BindingDB coverage is dominated by NRI/TCA scaffolds (reboxetine, atomoxetine, desipramine analogs); same mechanism as DAT, slightly less saturated.

- **GRIN2B**: Scenario 3 with high confidence — dominant ChEMBL pharmacology is ifenprodil-site allosteric; MAMMAL's protein-sequence input physically cannot represent the GluN1/GluN2B dimer interface. Deprecate from MAMMAL panel; rely on Boltz-2 with explicit dimer template.

- **GRIN2A**: Scenario 5 (insufficient n) blocking attribution to 1 or 3. Expand n to ≥30 before deciding anything.

The 1–2 day diagnostic produces these per-target attributions concretely. The paper then writes itself: protocol + 4-target attribution table + explicit per-scenario action map. That is the publishable scientific deliverable, distinct from the engineering follow-up.

---

## Caveats

1. **The reported MAMMAL DTI metric is NRMSE only.** No Pearson r, RMSE, MAE, or Spearman ρ appears in the paper or model card. Direct rank-correlation comparisons against the published benchmark are not possible.
2. **The MAMMAL paper's DTI benchmark (PEER 4-class holdout) is not the same evaluation regime as the HuggingFace fine-tuned head (TDC Drug+Target cold-split), and neither matches Pierce's interpolation use case.** Quoting paper performance as evidence of expected DAT/NET behavior is invalid.
3. **The K-S = 0.5 threshold in Diagnostic B is engineering judgment, not literature precedent.** SIMPD (Landrum et al. 2023) provides the multi-objective GA machinery for simulated time-splits but does not prescribe this cutoff.
4. **GRIN2A and GRIN2B inversions are not statistically distinguishable from zero** at observed n. Diagnostic outputs at these targets must be reported with explicit Bonett-Wright CIs.
5. **Boltz-2 affinity-head training distribution skews kinase/GPCR/NR.** Transporters and NMDA receptors are OOD for Boltz-2 too; shared failure between MAMMAL and Boltz-2 at DAT does not uniquely implicate sequence representation.
6. **Karakas et al. 2011 GluN2B residue numbering is rat**; human (UniProt Q13224) requires offset cross-check before Diagnostic C is meaningful.
7. **The diagnostic protocol assumes Pierce can reproduce MAMMAL preprocessing exactly** (`harmonize_affinities(mode='max_affinity')`, the published normalization constants). If the existing DTI grid parquet used different preprocessing, Diagnostic B comparisons will be confounded — verify before running.