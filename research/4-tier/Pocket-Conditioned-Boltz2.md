# §7.5 — Pocket-Conditioned Boltz-2 Inference with Per-Pose Provenance Flagging for the MAMMAL Cognitive-Enhancement Repurposing Pipeline

## TL;DR

- **Build a 4-detector consensus (P2Rank + PocketMiner + CryptoBench + Boltz-2 cofold) and curate per-target orthosteric/allosteric centroid databases for the 7 priority targets (CHRNA7, GRIN2B, PDE4D, SIGMAR1, HRH3, DRD1, ACHE) first — this is the minimum viable §7.5.** Total compute cost is dominated by PocketMiner (~12 min for the 22-target panel on an RTX 5070), not by training; the gating value comes from a curated centroid table, not from more ML.
- **The textbook payoff is PDE4D: BPN14770 binds the UCR2 allosteric helix (PDB 6NJJ) while rolipram/piclamilast bind the catalytic site (PDB 1XOR / 1XOQ family). Boltz-2 already preserved this ordering (BPN14770 0.963 > rolipram 0.907 in Phase 0.5) but cannot currently label it.** Once labelled, you can promote BPN14770-like allosteric NAMs and demote orthosteric inhibitors via the §8.0b emesis-liability gate.
- **Auto-demote `surface_artifact` and `no_pocket_match` (×0.5 rank), flag `cryptic_predicted` and `allosteric_putative` for review, and ship `orthosteric` / `allosteric_known` without penalty.** Wire the same pocket_class column into §8.0b (5-HT2B agonism, hERG pore, HRH1, CB1 → hard-cut only on orthosteric matches) and §8.1 (PDE4D-allosteric, CHRNA7-allosteric, DRD1-PAM as new facets).

---

## Key Findings

1. **Pocket-detector landscape (May 2026)**: Three open-source detectors cover complementary regimes. P2Rank (Krivák & Hoksza, *J. Cheminform.* 10:39, 2018; DOI 10.1186/s13321-018-0285-8) is a random-forest on solvent-accessible-surface points; the published paper states only that "P2Rank belongs to the fastest available tools (requires under 1 s for prediction on one protein)" — informal benchmarks at "a few seconds per target" are widely cited but not from the original paper. It is the established baseline for orthosteric pockets. PocketMiner (Meller et al., *Nat. Commun.* 14:1177, 2023; DOI 10.1038/s41467-023-36699-3) is a GVP graph neural network trained on MD trajectories; the paper reports **"ROC-AUC: 0.87 …> 1,000-fold faster than existing methods"** on 39 experimentally confirmed cryptic pockets. CryptoBench (Škrhák et al., *Bioinformatics* 41(1):btae745, December 2024; DOI 10.1093/bioinformatics/btae745) uses ESM-2-3B residue embeddings + an MLP head and **"significantly outperforms PocketMiner in key metrics such as AUC and AUPRC"** on the CryptoBench test set. They are not redundant — they fail differently — which is exactly what you want in a consensus.

2. **Boltz-2 has a documented pocket-quality dependency**: from the Boltz-2 preprint (Passaro et al., bioRxiv 2025.06.14.659707): **"If the model fails to identify the correct pocket or inaccurately reconstructs the binding interface or conformational state of the protein, downstream affinity predictions are unlikely to be reliable."** This is the structural reason that ungated Boltz-2 affinity rankings need pocket provenance — it is a known failure mode of the model itself, not a hypothetical.

3. **The PDE4D allosteric site is structurally well-defined and is the highest-value test case in this panel.** The mechanism was established by Burgin et al. (*Nat. Biotechnol.* 28:63, 2010) on UCR2 closure; BPN14770 was co-crystallized in PDB **6NJJ** ("Crystal Structure of the PDE4D Catalytic Domain and UCR2 Regulatory Helix with BPN14770") and BPN5004 in **6BOJ**. The mechanism (Wilson et al., *Nat. Commun.* 9:3334, 2018): **"PDE4D allosteric inhibitors bind in the catalytic site and complete a hydrophobic surface that allows closure of the amphipathic UCR2 regulatory helix … Closure of UCR2 inhibits the access of cAMP to the catalytic site."** The pocket centroid that distinguishes BPN14770 from rolipram is anchored by **PDE4D Phe271 (UCR2) — F instead of Y, the single residue that makes the subtype selectivity work in primate PDE4D**. Operationally: a Boltz-2 pose within 5 Å of Phe271 Cα with the UCR2 helix closed = allosteric_known; a pose deep in the catalytic gorge near Gln610/Phe613 with no UCR2 contact = orthosteric.

4. **CHRNA7 has TWO mechanistically distinct allosteric sites that must be tracked separately.** PDB **7KOX** (Noviello et al., *Cell* 184:2121, 2021) shows the human α7 nAChR bound to epibatidine (orthosteric) AND PNU-120596 (type-II PAM) at the **TMD intersubunit interface** in the activated state; sister structures **7KOO** (resting/α-bungarotoxin) and **7KOQ** (desensitized/epibatidine) provide the conformational ensemble. Galantamine is **NOT** at the type-II site: Ludwig et al. (*J. Recept. Signal Transduct.* 30:469, 2010) localized the galantamine binding site to **β-strand 10 residues T197, I196, F198 in the extracellular domain**: "residues T197, I196, and F198 of β-strand 10 represent major elements of the galantamine binding site." Subsequent photo-affinity work (Hamouda, Kimm & Cohen, "Physostigmine and galanthamine bind in the presence of agonist at the canonical and noncanonical subunit interfaces of a nicotinic acetylcholine receptor," *J. Neurosci.* 33(2):485–494, 2013; DOI 10.1523/JNEUROSCI.3483-12.2013) flagged additional ECD subunit-interface sites that did not overlap β-strand 10. Note the literature dispute: Kowal et al. (*Br. J. Pharmacol.* 2018) titled their paper **"Galantamine is not a positive allosteric modulator of human α4β2 or α7 nicotinic acetylcholine receptors"** — a real disagreement Pierce should flag honestly. For §7.5 we adopt the Maelicke/Ludwig model and label galantamine as `allosteric_known (type-I, β-strand 10 / ECD)` with a note.

5. **GRIN2B/NR2B ifenprodil site is a heterodimer interface invisible to single-chain inference.** Karakas et al. (*Nature* 475:249, 2011; PDB **3QEL**) showed **"the GluN1 and GluN2B ATDs form heterodimer and that phenylethanolamine binds at the GluN1-GluN2B subunit interface rather than within the GluN2B cleft"**, with the essential residues **GluN1-A75, GluN2B-I82, GluN2B-F114**. EVT-101 was later shown (Stroebel et al., PMC4859819) to occupy a partially overlapping but distinct cavity. **Operational implication**: MAMMAL receives only the GluN2B (Q13224) sequence — the ifenprodil pocket physically does not exist in single-chain mode. Pocket classification at GRIN2B must run Boltz-2 cofolding on the GluN1+GluN2B ATD heterodimer (templated on 3QEL or the full-receptor 5FXG/5FXH/5FXI ensemble of Tajima et al.) or all ifenprodil-class poses will appear as `surface_artifact` on the single chain.

6. **DRD1 has a clean two-pocket dichotomy with PDB-quality ground truth.** The orthosteric pocket is the canonical TM3/TM5/TM6 aminergic cleft. The PAM site for LY3154207 (mevidalen) is in a **membrane-embedded pocket formed by ICL2 / TM3 / TM4 about 12 Å from the orthosteric site**. PDB codes from cryo-EM:
   - **7CKZ** (Xiao P, Yan W, Gou L et al., "Ligand recognition and allosteric regulation of DRD1-Gs signaling complexes," *Cell* 184(4):943–956.e18, 2021; DOI 10.1016/j.cell.2021.01.028) — DRD1 + dopamine + LY3154207 + Gs.
   - **7LJD** (Zhuang et al., *Cell Res.* 31:593, 2021; DOI 10.1038/s41422-021-00482-0) — DRD1 + dopamine + LY3154207 + miniGs.
   - **7LJC** (same paper) — DRD1 + SKF-81297 + LY3154207 + miniGs. This is the cleanest reference for the PAM cavity itself.
   - **7X2F** — additional LY3154207 complex catalogued by later virtual-screening studies.
   The downstream MD-derived mechanism (LY3154207 stabilizing ICL2 helical conformation and tightening the Na⁺ binding cluster) is from Teng et al., *Nat. Commun.* 13:3186, 2022 (DOI 10.1038/s41467-022-30929-w), not the 2021 cryo-EM papers.

7. **HRH3 has at least two usable structural references.** Crystal structure **7F61** (deposited 2022, released 2022-10-26) shows human H3R + PF-03654746 (antagonist, orthosteric). The 2025 cryo-EM paper (Wang et al., "Decoding ligand recognition and constitutive activation of histamine H3 and H4 receptors," *Acta Pharmacol. Sin.* 2025; DOI 10.1038/s41401-025-01633-4) reports **"two Gi-coupled structures of H3R and H4R in complex with histamine"** with the active-state H3R-histamine-Gi complex deposited as **PDB 8YUU**. Pitolisant is documented as Ki = 0.16 nM (competitive antagonist) and EC50 = 1.5 nM (inverse agonist) at recombinant human H3R per Ligneau et al. 2007 / the Wakix EPAR (EMA/CHMP/4151/2015) and is expected to occupy the same canonical aminergic pocket; use 7F61 as the inactive-state orthosteric template and 8YUU as the active-state reference.

8. **ACHE has a single gorge but two functional sub-sites that must be coded as separate centroids.** PDB **4EY7** (Cheung et al., 2012) shows donepezil binding **both** the catalytic anionic site (CAS; Trp86, Tyr133, Glu202, Ser203, Tyr337, Phe338, His447) and the peripheral anionic site (PAS; Tyr72, Asp74, Tyr124, Trp286, Phe295, Tyr341) — the dual-site binding mechanism that distinguishes donepezil from tacrine (CAS-only; PDB 1ACJ). Galantamine sits primarily in CAS with weak PAS contact (PDB 4M0E / 4EY6 family). For §7.5, code CAS and PAS as separate centroids; a "dual" pose is one whose ligand heavy-atom span makes contact within 5 Å of **both** centroids.

9. **SIGMAR1 has a single occluded β-barrel cavity and ligand class (agonist vs antagonist) is the discriminator, not pocket.** PDB **5HK1** (PD144418, antagonist) and **5HK2** (4-IBP, agonist) from Schmidt et al. (*Nature* 532:527, 2016) show the same cavity used by both; Glu172 forms the obligatory salt bridge ("the highly conserved Glu172 located near the center of the cavity, forming a salt bridge with the ligands"). Newer Xenopus σ1R structures (PDB 7W2B–7W2H, from Zhou X et al., "An open-like conformation of the sigma-1 receptor reveals its ligand entry pathway," *Nat. Commun.* 13:1267, 2022; DOI 10.1038/s41467-022-28946-w) reveal closed and open conformations. For §7.5, σ1R is a one-pocket target with conformation conditioning, not a multi-pocket target.

10. **DAT (SLC6A3) poses a special problem and is a key lateral lever.** Nielsen et al. (*Nature* 2024, DOI 10.1038/s41586-024-07804-3) deposited a cryo-EM structure of hDAT-cocaine showing cocaine at the central S1 site. Pedersen et al. (*J. Neurochem.* 168:2043, 2024) showed AC-4-248 at the orthosteric site of dDAT. A separate S2 vestibule site has been characterized computationally (Cheng et al., *J. Chem. Inf. Model.* 2020, DOI 10.1021/acs.jcim.0c00346): **"cocaine binding to this new site allosterically reduces the binding of DA/cocaine to the central binding pocket."** This is the right lateral lever for §7.11 isotonic + §7.13 active learning: if MAMMAL saturates on tropanes at S1, AL should prefer non-tropane scaffolds predicted to dock at S2.

---

## Details

### 1. Comparative Methodology Table

| Method | Input | Output | Wall-clock (22 targets) | Strength | Weakness | Recommended role |
|---|---|---|---|---|---|---|
| **P2Rank 2.4** (Krivák & Hoksza 2018) | PDB / mmCIF | Ranked pocket centroids + ligandability score + per-residue assignments | ~1 min total (<1 s/target per published claim; informal benchmarks ~3 s cold) | Mature, single Java jar, integrates with PrankWeb 4 (Polak et al., *NAR* 53:W466, 2025) and PDBe-KB; rescores Fpocket, DeepSite, PUResNetV2.0 outputs | Trained on holo structures; biased toward orthosteric concave sites; cryptic recall poor on apo | **Baseline triage on all 22 targets**; rescore step for any other detector |
| **PocketMiner** (Meller 2023) | Single apo PDB | Per-residue cryptic-opening probability | ~12 min (~30 s/target on RTX 5070) | Specifically trained on MD-derived cryptic-opening labels; ROC-AUC 0.87 on 39 experimentally confirmed cryptics; >1000× faster than running MD yourself | Output is broad/smooth (independent reviewers note non-trivial false-positive rate); training labels depend on LIGSITE pocket volume during MD, so pockets that open on slow ringflip/SS timescales are under-represented | **Targeted cryptic flag** on CHRNA7 apo, HRH3 apo, KCNQ2/3 retigabine site, GRIN2B ATD |
| **CryptoBench-pLM-NN** (Škrhák 2024) | Sequence (UniProt) | Per-residue cryptic-binding-residue probability | ~5 min one-time head training + ~5 s/target inference | Sequence-only (reuses your ESM-2 cache); outperforms PocketMiner on the CryptoBench test set in AUC, AUPRC, MCC, F1 | New (Dec 2024 Bioinformatics), limited independent validation; uses ESM-2-3B embeddings (your cache is 650M — check dimension compatibility before claiming "free") | **Sequence-side cross-check** on every target; primary detector for unstable/AlphaFold-only targets |
| **Boltz-2 cofolding** (Passaro et al. 2025) | Sequence + SMILES | Per-pose structure + binder/affinity head | Free (already in overnight sweep) | Reveals genuine cryptic openings on the specific compound-target pair; same model whose affinity we already trust | Pocket "by inspection"; compute-bound; affinity head explicitly warned to be unreliable when pocket is wrong (per Boltz-2 "Limitations" section) | **Second opinion** on any compound where Boltz-2 affinity disagrees with MAMMAL pKd by >1.5 log units |

**Deployment decision tree (per target × per compound)**:
```
1. ALWAYS run P2Rank on the curated PDB → pocket-list with centroids
2. For 7 priority targets: ALSO run PocketMiner (apo state) + CryptoBench
3. For Boltz-2 pose P at target T with predicted affinity A:
     d_ortho     = min_a∈ortho_centroids(T)    ||P_centroid − a||
     d_allo      = min_a∈allo_centroids(T)     ||P_centroid − a||
     d_predicted = min_a∈P2Rank∪PM∪CB(T)       ||P_centroid − a||
     
     IF d_ortho ≤ 8 Å                          → orthosteric
     ELIF d_allo  ≤ 8 Å                        → allosteric_known
     ELIF d_predicted ≤ 8 Å AND ≥2 detectors agree
                                                → allosteric_putative or cryptic_predicted
                                                  (allosteric_putative if site is buried;
                                                   cryptic_predicted if PocketMiner score >0.5
                                                   and CryptoBench score >0.5)
     ELIF d_any ≤ 15 Å                         → no_pocket_match
     ELSE (solvent-exposed, SASA > 80 Å²)      → surface_artifact
4. IF MAMMAL pKd vs Boltz-2 affinity disagree by >1.5 log → trigger Boltz-2 cofold
   second pass and re-classify
```

### 2. Per-Target Pocket Curation (7 priority targets)

**Format**: target → UniProt → reference PDBs → orthosteric centroid + residues → allosteric centroid(s) + residues → decision-rule special cases.

**2A. CHRNA7 — α7 nAChR (UniProt P36544)**
- Orthosteric: ACh-binding site at α7/α7 ECD interface, 5 sites per pentamer; principal-face residues Tyr92, Trp148, Tyr187, Tyr194 (loop C); complementary-face residues Trp53, Leu108, Gln114. Centroid = heavy-atom mean of Tyr92, Trp148, Tyr187, Tyr194, Trp53 across interfaces. Reference: PDB **7KOQ** (epibatidine-bound desensitized).
- Allosteric type-II (TMD intersubunit): PNU-120596 site between M1/M2/M3 helices of adjacent subunits; key residues Met253, Ser223 (M1), Leu247 (M2), Ile281 (M3). Reference: PDB **7KOX** (epibatidine + PNU-120596 activated). Wang et al. 2025 (*Cell Discov.* DOI s41421-025-00788-y) note GAT-107 occupies a position spanning "from the middle to the top of the transmembrane helix" in the same pocket but with vertical orientation differing from PNU.
- Allosteric type-I (ECD β-strand 10 / vestibule): galantamine APL site; key residues T197, I196, F198 (Ludwig et al. 2010); additional ECD subunit-interface sites flagged by Hamouda, Kimm & Cohen 2013 (*J. Neurosci.* 33:485). Centroid is on the principal face of the ECD, above the orthosteric C-loop. A separate vestibule allosteric pocket near the β8–β9 loop was identified in α7-AChBP (Spurny et al., PMC5818190, 2018).
- Conformational ensemble fallback: 7KOO (resting), 7KOQ (desensitized), 7KOX (activated). Run PocketMiner on all three states; cryptic opening at the TMD intersubunit site is state-dependent.
- **Disagreement honestly noted**: Kowal et al. 2018 (*Br. J. Pharmacol.*) report no PAM effect from galantamine on human α7. Label galantamine `allosteric_known` only if Phase 0.5 evidence holds; otherwise `allosteric_putative` with the Kowal caveat.

**2B. GRIN2B — NR2B (UniProt Q13224)**
- LBD orthosteric (glutamate binding cleft): principal residues Tyr754, Ser690, Thr691, Ser688; bi-lobed clamshell. Heterodimer not required for this site. Reference: PDB 4PE5 family.
- ATD allosteric (ifenprodil / Ro 25-6981 / EVT-101 / traxoprodil): **at the GluN1-GluN2B ATD heterodimer interface**, NOT within the GluN2B cleft. Key residues GluN1-Ala75, GluN2B-Ile82, GluN2B-Phe114, GluN1-Ser132 (H-bond to ligand hydroxyl), GluN2B-Gln110 (Karakas et al. 2011 mutagenesis: Q110A reduces ifenprodil potency 10-fold). EVT-101 occupies a partially overlapping cavity (PMC4859819). References: PDB **3QEL** (GluN1b-GluN2B ATD + ifenprodil), 5IOV family, 5EWJ/5EWL family (EVT-101, MK-22), 5FXG/5FXH/5FXI/5FXJ (full-tetramer states from Tajima et al.).
- **Critical engineering note**: Boltz-2 cofold for ifenprodil-class compounds MUST include both GluN1 (UniProt Q05586) and GluN2B chains, with template constraint to 3QEL. Single-chain GluN2B cofolds will systematically misclassify ifenprodil as `surface_artifact` because the pocket physically does not exist.

**2C. PDE4D (UniProt Q08499)**
- Orthosteric (catalytic): cAMP-binding pocket; conserved PDE active-site residues Gln610 (in 6NJJ numbering), Ile577, Phe613, two divalent metals (Zn²⁺ + Mg²⁺). Rolipram, piclamilast, roflumilast, GEBR-7B bind here. Reference: PDB 1XOR (rolipram), 1XOQ.
- Allosteric (UCR2 closure): UCR2 helix folds over the catalytic site; **Phe271 (UCR2, primate-specific Y→F switch)** anchors the π–π edge-on stack with BPN14770's pyrimidine core. Mechanism (Wilson et al. 2018): the inhibitor "completes a hydrophobic surface that allows closure of the amphipathic UCR2 regulatory helix." Centroid = midpoint between catalytic Q610 and UCR2 F271 with UCR2 in closed conformation. References: **6NJJ** (BPN14770), **6NJH** (T-48), **6NJI** (T-49), **6BOJ** (BPN5004).
- Decision rule: pose with heavy atoms within 5 Å of BOTH F271 (UCR2) and Q610 (catalytic) = `allosteric_known` (BPN14770-like); pose within 5 Å of Q610/F613 catalytic only, no UCR2 contact = `orthosteric` (rolipram-like). The two classes route to different §8.0b emetic-liability rules.

**2D. SIGMAR1 (UniProt Q99720)**
- Single ligand-binding cavity in cupin-like β-barrel; obligatory residues **Glu172** (salt bridge to ligand amine), **Tyr103**, **Asp126** (H-bond network), Met93, Leu95, Tyr120. Antagonist (PD144418, 5HK1) and agonist (4-IBP, 5HK2; PRE-084, 6DK1) occupy the same cavity. Discriminator is ligand chemotype + conformational state (closed vs open, Xenopus structures 7W2B–7W2H from Zhou et al. 2022 *Nat. Commun.* 13:1267), not pocket location.
- Pocket-class assignment: `orthosteric` for any pose within 6 Å of Glu172 forming a salt bridge to a basic nitrogen; `surface_artifact` otherwise. Agonist/antagonist split is delegated to a downstream functional classifier (out of scope for §7.5).

**2E. HRH3 (UniProt Q9Y5N1)**
- Single canonical aminergic GPCR pocket: TM3 Asp114³·³² (salt bridge to amine), TM5/TM6 aromatic stack, ECL2 cap. Reference inactive: **PDB 7F61** (H3R + PF-03654746 antagonist, deposited 2022). Reference active: **PDB 8YUU** (H3R-histamine-Gi cryo-EM, Wang et al. 2025 *Acta Pharmacol. Sin.*). Pitolisant (Ki = 0.16 nM antagonist, EC50 = 1.5 nM inverse agonist at recombinant human H3R; Wakix EPAR EMA/CHMP/4151/2015) is expected to occupy the 7F61-style inactive pose.
- No widely accepted secondary allosteric site for H3R yet — literature is preliminary. Anything Boltz-2 docks far from D114 should be flagged `cryptic_predicted` and reviewed before promotion.

**2F. DRD1 (UniProt P21728)**
- Orthosteric (OBP + EBP extended): TM3 Asp103³·³², TM5 Ser198/Ser202, TM6 W321, TM7. Fenoldopam binds with two molecules (one at OBP, one at EBP); tavapadon spans both. References: PDB 7CRH (apo-like), 7JOZ (Gs complex), 7JV5 (orthosteric agonist).
- Allosteric PAM (LY3154207 / mevidalen / DETQ): **membrane-embedded pocket formed by ICL2 / TM3 / TM4**, ~12 Å from OBP, predominantly hydrophobic. References: **7LJD** (dopamine + LY3154207 + miniGs, Zhuang et al. *Cell Res.* 2021), **7LJC** (SKF-81297 + LY3154207 + miniGs — the cleanest reference for the PAM cavity), **7CKZ** (Xiao P, Yan W, Gou L et al., *Cell* 184:943, 2021; senior co-author Jin-Peng Sun), **7X2F**. Downstream MD (Teng et al., *Nat. Commun.* 13:3186, 2022) shows the PAM stabilizes ICL2 helical conformation and tightens the Na⁺ binding cluster.
- Pocket-class assignment: pose centroid within 6 Å of Asp103 = `orthosteric`; centroid within 6 Å of the ICL2 helical region (~residues 130–145) on the lipid-facing side = `allosteric_known (PAM)`.

**2G. ACHE (UniProt P22303)**
- Two centroids in a single deep gorge:
  - **CAS**: Trp86, Tyr133, Glu202, Ser203 (catalytic triad nucleophile), Tyr337, Phe338, His447.
  - **PAS**: Tyr72, Asp74, Tyr124, Trp286, Phe295, Tyr341 — at the gorge mouth.
  - Mid-gorge linker: Asp74, Leu76, Phe297, Phe338.
- Reference: **PDB 4EY7** (donepezil dual-site CAS+PAS), 4EY6 (galantamine), 4EY5 (huperzine A, CAS), 1ACJ (tacrine, CAS-only).
- Pocket-class rule: `orthosteric` if pose contacts CAS within 5 Å OR dual-CAS+PAS within 5 Å of each centroid (donepezil-mode); `allosteric_known (PAS)` if pose contacts only PAS within 5 Å with no CAS contact (rare in approved drugs, common in propidium-class probes); `surface_artifact` if outside the gorge entirely.

### 3. Provenance Schema and Integration

**3.1 Column specification**
```python
# v5 shortlist columns
pocket_class: Literal[
    "orthosteric",
    "allosteric_known",
    "allosteric_putative",
    "cryptic_predicted",
    "surface_artifact",
    "no_pocket_match",
    "NA_no_pose",
]
pocket_confidence: float          # ∈ [0, 1], weighted across detectors
pocket_reference_pdb: str         # the PDB anchor used (e.g., "6NJJ")
pocket_centroid_xyz: tuple[float, float, float]
detector_votes: dict[str, float]  # {"p2rank":0.91, "pocketminer":0.42, "cryptobench":0.55, "boltz2_cofold":1.0}
```

**3.2 Auto-demotion logic in fusion**
```python
RANK_MULTIPLIER = {
    "orthosteric":         1.00,
    "allosteric_known":    1.00,
    "allosteric_putative": 1.00,   # but tagged manual_review=True
    "cryptic_predicted":   1.00,   # tagged manual_review=True
    "no_pocket_match":     0.50,
    "surface_artifact":    0.50,
    "NA_no_pose":          0.30,
}
# applied AFTER RRF fusion, BEFORE the per-target Spearman ρ calibration
```

**3.3 §8.0b liability-panel integration (the pharmacological-granularity upgrade)**

| Off-target | Orthosteric pose | Allosteric pose | Rationale |
|---|---|---|---|
| 5-HT2B | **HARD CUT** if predicted agonist | **FLAG** | Fenfluramine/norfenfluramine valvulopathy precedent — Rothman et al. (*Circulation* 102:2836, 2000) showed direct 5-HT2B agonism drives valve interstitial cell mitogenesis. |
| KCNH2 (hERG) | **HARD CUT** if pose in central pore cavity (Y652, F656, T623) | **FLAG** | Auxiliary-subunit / membrane-vestibule binding has less clean torsadogenic precedent. |
| HRH1 | **HARD CUT** if predicted antagonist at orthosteric site | **FLAG** | Anticholinergic-burden + dementia: Gray et al. *JAMA Intern. Med.* 175:401, 2015 (HR 1.54 [95% CI 1.21–1.96] at >1095 TSDD); Coupland et al. *JAMA Intern. Med.* 179(8):1084–1093, 2019; DOI 10.1001/jamainternmed.2019.0677 (nested case-control of 58,769 dementia cases vs 225,574 matched controls; **adjusted odds ratio 1.49 [95% CI 1.44–1.54]** for the highest cumulative exposure category, >1095 TSDD over 1–11 years pre-index). Allosteric H1 binding has no comparable dementia signal. |
| CB1 | **HARD CUT** at orthosteric pose (rimonabant precedent — withdrawn for suicidality) | **FLAG** | Allosteric CB1 NAMs (Org-27569 class) have a different clinical risk profile. |

Without pocket conditioning, every 5-HT2B / hERG / HRH1 / CB1 hit looks identical to the gate. With pocket conditioning, only the mechanistically-loaded pose triggers the hard cut — preserving compounds that touch these proteins at non-canonical sites for further triage.

**3.4 §8.1 selectivity-faceting integration**

New facets enabled:
- **PDE4D-allosteric facet**: BPN14770, GEBR-7B family, GEBR-32a, ApremiPDE4D-selective leads.
- **CHRNA7-allosteric (type-II) facet**: PNU-120596, GAT107, A-867744.
- **CHRNA7-allosteric (type-I) facet**: galantamine, NS-1738, 5-hydroxyindole.
- **DRD1-PAM facet**: LY3154207 (mevidalen), DETQ, Lilly compound B series.
- **ACHE dual-site facet**: donepezil-class (CAS+PAS span), distinct from CAS-only (tacrine, rivastigmine).

Within-facet bonus: +20% selectivity score if the compound binds the cognition-relevant pocket and is *selective* against the canonical-pocket sibling (e.g., PDE4D-UCR2 selective vs PDE4D-catalytic). Cross-facet conflict resolution: assigned to the facet closer to a clinical/Phase-2 lead.

**3.5 Module layout**
```
src/mammal_repurposing/pockets/
├── pocket_database.py        # curated centroid table for 22 targets
├── pocket_detector.py        # P2Rank / PocketMiner / CryptoBench routers
├── boltz2_cofold_oracle.py   # second-opinion runner on disagreement pairs
├── pocket_classifier.py      # geometric assignment to pocket_class
├── consensus.py              # 4-detector vote aggregation
└── provenance.py             # joins pocket_class onto the v5 shortlist
```

### 4. Expected Impact on the v3 Top-25

| Compound | Top target | Predicted pocket_class | Action |
|---|---|---|---|
| BPN14770 | PDE4D | allosteric_known (UCR2, F271 anchor) | ship; anchor of new PDE4D-allosteric facet |
| rolipram | PDE4D | orthosteric (catalytic Q610) | flag emetic-liability via §8.0b; do NOT demote, but mark as orthosteric class |
| piclamilast | PDE4D | orthosteric | same as rolipram |
| TC-5619 | CHRNA7 | orthosteric (α7/α7 ECD interface) | ship |
| encenicline | CHRNA7 | orthosteric (partial agonist) | ship |
| galantamine | CHRNA7 | allosteric_known (β-strand 10 / type-I) WITH Kowal caveat | ship to type-I PAM facet, manual_review=True |
| galantamine | ACHE | orthosteric (CAS + minor PAS) | ship |
| donepezil (if added) | ACHE | orthosteric (dual CAS+PAS) | anchor of dual-site facet |
| pitolisant | HRH3 | orthosteric (D114-anchored aminergic pocket) | ship |
| 7,8-DHF | NTRK2 | allosteric_putative (TrkB-D5 binding; not the BDNF-ECD primary site, modeled at K312/L315/I334; literature is 2D-docking-only) | manual_review; pocket database fragile |
| blarcamesine (ANAVEX 2-73) | SIGMAR1 | orthosteric (single cavity) | ship |
| 2BAct | PDE4D | surface_artifact (true target eIF2B, off-panel) | demote ×0.5 |
| methylphenidate | SLC6A3 | orthosteric (S1 central) | ship after §7.11 sign-flip |
| solriamfetol | SLC6A3 + SLC6A2 | orthosteric S1 | ship |
| aripiprazole | DRD1 | orthosteric | already CUT by §8.0b via 5-HT2B + α1 + D2 |
| amitriptyline | HRH1 (off-target) | orthosteric | CUT by §8.0b (HRH1 + anticholinergic) |
| (R,S)-AMPA | GRIA1-4 | orthosteric (LBD glutamate cleft) | ship (research tool) |
| lemborexant | HCRTR2 | orthosteric | ship |
| lithium carbonate | none in panel directly | NA_no_pose (Boltz-2 returns no confident pose for Li⁺ as small inorganic) | demote heavily; expect manual handling |

Predicted overall distribution across v3 top-25: ~55% orthosteric, ~20% allosteric_known, ~10% allosteric_putative, ~5% cryptic_predicted, ~10% surface_artifact / no_pocket_match. The 10% surface_artifact band is where §7.5 earns its keep — those are the silent failures of the current pipeline.

### 5. Validation Gates

- **Gate P1 (orthosteric positive controls)**: pitolisant @ HRH3 (D114 contact); donepezil @ ACHE (dual-site CAS+PAS); methylphenidate @ SLC6A3 (S1); (R,S)-AMPA @ GRIA1-4 (LBD cleft). **All must classify `orthosteric`.**
- **Gate P2 (allosteric_known positive controls)**: BPN14770 @ PDE4D (UCR2 F271 contact); PNU-120596 @ CHRNA7 (TMD intersubunit, type-II); LY3154207 @ DRD1 (ICL2/TM3/TM4 PAM site); ifenprodil @ GluN1+GluN2B heterodimer (ATD interface; **requires Boltz-2 cofold of heterodimer**). **All must classify `allosteric_known`.**
- **Gate P3 (negative controls)**: 2BAct @ PDE4D classified `surface_artifact` (true target is eIF2B); aspirin @ HRH3 classified `surface_artifact`; a random promiscuous binder (e.g., suramin) at PDE4D should not classify `allosteric_known`.
- **Gate P4 (detector consensus)**: For the 7 priority targets, P2Rank + PocketMiner + CryptoBench must agree on ≥70% of compound-target pairs at the pocket-class level. Disagreement → manual review queue.
- **Gate P5 (no regression)**: §8.1 v5 top-5 within each existing facet must still contain TC-5619 + encenicline (CHRNA7), donepezil (ACHE), pitolisant (HRH3), BPN14770 (PDE4D-allosteric facet specifically), LY3154207 (DRD1-PAM facet).
- **Gate P6 (galantamine controversy)**: explicitly check what happens when Kowal et al. 2018 is taken as ground truth and galantamine is relabeled `allosteric_putative`. The pipeline must not crash and the type-I facet must remain meaningful.

### 6. Lateral Angles

**6.1 Conformational-ensemble pocket prediction.** CHRNA7 has 7KOO (resting) + 7KOQ (desensitized) + 7KOX (activated) — three states where the TMD intersubunit pocket has different geometry. Run P2Rank and PocketMiner separately on each state and merge centroids. A pose that classifies `allosteric_known (type-II)` in the activated state but `surface_artifact` in the resting state has a known mechanistic interpretation (state-dependent binding) rather than being a contradiction. Do the same with the SIGMAR1 closed/open ensemble (7W2B–7W2H, Zhou et al. 2022) and the GluN1-GluN2B ATD active/inactive ensemble (5FXG/H/I/J).

**6.2 ChEMBL `binding_site` stratified recalibration.** ChEMBL records carry a `binding_site` field when reported. Re-run §7's Phase A.7 Spearman-ρ-vs-pchembl calibration STRATIFIED by binding_site. If MAMMAL's ρ at PDE4D differs between orthosteric ChEMBL records and allosteric ChEMBL records, you have a quantitative measure of how much of the rescue Boltz-2 is actually doing. This becomes a pre-registered Phase 0.6 analysis.

**6.3 §7.13 active-learning routing.** If §7.1 confirms tropane saturation at SLC6A3, AL should preferentially sample non-tropane scaffolds predicted to bind DAT S2 vestibule (per Cheng et al. 2020 J. Chem. Inf. Model.: "cocaine binding to this new site allosterically reduces the binding of DA/cocaine to the central binding pocket"). The S2 centroid lives near Trp51, Phe471, His472, Asp475, Tyr547.

**6.4 §7.5 + §7.11 monotone-isotonic break.** The §7.11 isotonic calibrator at SLC6A3/SLC6A2 assumed a single monotone MAMMAL-pKd-vs-ChEMBL-pchembl relationship. If compounds bind S1 vs S2 at DAT, the assumption is violated and the isotonic fit will mix two regimes. Solution: route by `pocket_class` BEFORE isotonic; fit separate calibrators for S1 (orthosteric) and S2 (cryptic_predicted).

**6.5 Cross-cluster disagreement operationalized.** MAMMAL is pocket-blind. Boltz-2 is pocket-aware (implicitly). When MAMMAL ranks compound C high at target T but the 4-detector consensus classifies the Boltz-2 pose as `surface_artifact`, the disagreement is mechanistically meaningful and should be exposed as an explicit demotion signal (already in `provenance/disagreement_report.py` but not currently typed; add `pocket_class_disagreement` field).

### 7. Publication Angle and Reproducibility

- **Title**: "Pocket-conditioned structure-aware affinity prediction for cognitive-enhancement drug repurposing: orthosteric / allosteric / cryptic / surface-artifact discrimination via a four-detector ensemble on a curated 22-target panel."
- **Venue priority**: *J. Cheminform.* (best home — same journal as P2Rank); *Bioinformatics* (CryptoBench home); *Drug Discovery Today: Technologies* (methods venue).
- **Novelty claims to make explicitly**:
  1. First repurposing pipeline that systematically classifies docked poses into orthosteric / allosteric_known / allosteric_putative / cryptic_predicted / surface_artifact, applied to a clinical-cognition target panel.
  2. First operational PDE4D-allosteric-vs-orthosteric split inside a foundation-model DTI pipeline (BPN14770 vs rolipram class).
  3. First systematic CHRNA7 type-I (galantamine) vs type-II (PNU-120596) PAM-site discrimination in a virtual-screening workflow.
  4. Pocket-class-conditioned liability gating (§8.0b) — converting "hERG hit" from a flat cut into "hERG-pore hit = cut, hERG-vestibule hit = flag."
- **Reproducibility deliverables**: curated pocket DB (YAML) for 22 targets; consensus code; gate definitions; positive/negative control list; honest disagreement notes (Kowal 2018 galantamine controversy, Boltz-2 OOD warnings at transporters/NMDA).

### 8. Risks and Failure Modes

1. **PocketMiner false-positive halo**. Independent reviewers flag broad, smooth output. Mitigation: require ≥2-detector agreement before promoting to `cryptic_predicted`.
2. **CryptoBench dimensional mismatch**. CryptoBench was trained on ESM-2-3B embeddings (2560-dim per residue); the cache is ESM-2-650M (1280-dim). Verify before claiming "free re-use" — if a re-embedding pass is required, budget ~30 min/target for embedding instead of "free."
3. **Boltz-2 affinity head OOD at transporters and NMDA**. The Boltz-2 paper explicitly warns affinity is structure-dependent; SLC6A3/SLC6A2/GRIN2A/GRIN2B are kinase/GPCR-poor in training distribution. Pocket classification at these targets must carry an OOD flag.
4. **GRIN2B heterodimer requirement**. Single-chain Boltz-2 cofold cannot reveal the ifenprodil pocket. Pipeline must template on 3QEL with GluN1 chain; otherwise every ifenprodil-class hit is silently misclassified.
5. **Galantamine literature dispute**. The β-strand 10 site (Ludwig et al. 2010) and the no-PAM-effect finding (Kowal et al. 2018) are both peer-reviewed. Pierce should label galantamine `allosteric_known` with an explicit citation pair and flag for manual triage.
6. **Curated centroid drift**. The 7 priority targets are curated as of the May 2026 PDB; newer cryo-EM structures may add or shift centroids. Schedule a 6-month re-curation pass.
7. **State-dependence cost**. Conformational ensembles multiply pocket-database entries; budget storage and disambiguation overhead.

---

## Recommendations

**Sprint 1 (this week — minimum viable §7.5)**:
1. Curate `pocket_database.yml` for the 7 priority targets using the PDB IDs above (CHRNA7: 7KOO/7KOQ/7KOX; GRIN2B: 3QEL+5IOV family; PDE4D: 6NJJ+1XOR; SIGMAR1: 5HK1/5HK2; HRH3: 7F61+8YUU; DRD1: 7LJC/7LJD/7CKZ; ACHE: 4EY7). Pull centroids via `MDAnalysis.select_atoms("around 5.0 resname <ligand>")` and store heavy-atom mean.
2. Install P2Rank 2.4 in WSL2; run on all 22 AlphaFold or PDB structures; cache to `data/cache/p2rank/`.
3. Build `pocket_classifier.py` with the decision tree above; require it to pass Gate P1 + Gate P2 against the curated set.
4. Wire `pocket_class` column into v5 shortlist; add to `provenance/disagreement_report.py`.

**Sprint 2 (next 2 weeks — full ensemble)**:
5. Install PocketMiner from `https://github.com/Mickdub/gvp` (branch `pocket_pred`); run on apo and ligand-bound states of the 7 priority targets.
6. Validate ESM-2 cache dimension; either retrain CryptoBench-pLM-NN head on ESM-2-650M (cheaper, ~10 min) or pay the re-embed cost for ESM-2-3B.
7. Implement consensus voting; pass Gate P4 (>70% agreement).
8. Wire pocket_class into §8.0b gate table; demonstrate that the 5-HT2B / hERG / HRH1 / CB1 hard-cuts only fire on orthosteric matches.

**Sprint 3 (3-4 weeks — facet upgrade + paper draft)**:
9. Add PDE4D-allosteric, CHRNA7-type-I, CHRNA7-type-II, DRD1-PAM, ACHE-dual-site facets to §8.1.
10. Run the stratified-ChEMBL Phase A.7 recalibration (lateral 6.2).
11. Wire pocket-routed isotonic at SLC6A3/SLC6A2 (lateral 6.4).
12. Draft *J. Cheminform.* manuscript with the four novelty claims; supplementary = the pocket DB YAML + control list + positive/negative gate logs.

**Thresholds that change these recommendations**:
- If Gate P4 (detector agreement) <50% on the 7 priority targets → drop CryptoBench, keep P2Rank+PocketMiner only.
- If Gate P2 fails on BPN14770 @ PDE4D → the §7.5 thesis is broken; investigate before scaling.
- If GRIN2B heterodimer cofold doesn't converge → defer GRIN2B from the priority-7 to a special-handling tier.

---

## Caveats

- **Boltz-2 affinity-head OOD warning** is explicit in the source paper: "If the model fails to identify the correct pocket … downstream affinity predictions are unlikely to be reliable" (Passaro et al. 2025). Pocket classification does not fix this — it only makes the failure mode visible.
- **PocketMiner is new (2023) and has independent reviewer concerns** about false-positive density. ROC-AUC 0.87 is on a curated set of 39 confirmed cryptics; real-world precision at decision thresholds matters more than ROC-AUC for pipeline use.
- **CryptoBench is newer still (Dec 2024)**; almost no independent validation outside the authors' benchmark.
- **The galantamine α7 PAM literature is genuinely contested** (Ludwig 2010 + Hamouda 2013 vs Kowal 2018). Do not assert a clean type-I-PAM mechanism without flagging the controversy.
- **GRIN2B and DAT analyses depend on heterodimer / multi-conformation cofolds** that the current pipeline does not run by default. The §7.5 promise at these targets is conditional on the cofold infrastructure being extended.
- **The pocket database is a snapshot.** Cryo-EM and crystallography move; expect HRH3, DRD1-PAM, and DAT structural ground truth to shift within 12–24 months.
- **Centroid coordinates in this report are residue-list-derived, not precomputed XYZ**: Pierce must compute heavy-atom means from the cited PDBs (script in `pocket_database.py`); we do not bake in numerical XYZ here because they depend on which PDB chain/alternate-location is taken.
- **Phase 0.5 numbers (TC-5619 100%, encenicline 80%, galantamine 40%, BPN14770 0.963 > rolipram 0.907, pitolisant 100%) are the pipeline's internal metrics** and have not yet been wet-validated; §7.5 sharpens these claims by tying each rescue to a named pocket, but does not replace wet-lab confirmation.
- **Coupland 2019 effect estimate is an adjusted odds ratio (AOR 1.49 [95% CI 1.44–1.54]) from a nested case-control design**, not a cohort hazard ratio as some secondary citations describe it. Cite carefully.
- **Coordinates and identities of the LY3154207-bound DRD1 PDB structures are dual-sourced**: Xiao et al. *Cell* 2021 deposited 7CKZ; Zhuang et al. *Cell Res.* 2021 deposited 7LJC and 7LJD. Use 7LJC as the cleanest PAM-cavity reference; the Na⁺-cluster compactness and ICL2-helicity mechanism findings come from later MD work (Teng et al. *Nat. Commun.* 13:3186, 2022), not the cryo-EM papers themselves.