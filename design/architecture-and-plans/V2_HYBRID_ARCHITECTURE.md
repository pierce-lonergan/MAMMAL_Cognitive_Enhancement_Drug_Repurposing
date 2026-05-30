# V2 — Hybrid Cluster Architecture

**Source**: [research/Hybrid Architecture for MAMMAL-Based Cognitive-Enhancement Drug Repurposing.md](../research/Hybrid%20Architecture%20for%20MAMMAL-Based%20Cognitive-Enhancement%20Drug%20Repurposing.md)

V1 (already shipped) was a single-cluster pipeline: MAMMAL DTI head → sanity gate → ranked list. The empirical failure modes we measured — peptide pollution, structural-bias inflation, allosteric blindness at CHRNA7, dynamic-range collapse — are exactly what a single-model approach predicts. V2 is the response: a 3-cluster hybrid that joins **target-affinity (MAMMAL)**, **physical safety (ADMET-AI)**, and **mechanism/network plausibility (PrimeKG + TxGNN)**, with optional pose-conditioned affinity (Boltz-2/Boltzina) as the structure-aware refinement layer.

## Architecture

```
                    ┌────────────────────────────────┐
                    │  V1 ARTIFACTS (frozen)         │
                    │  data/results/dti_scores.parq  │
                    │  data/results/aux_scores.parq  │  (BBBP/ClinTox - retained as audit)
                    └──────────────┬─────────────────┘
                                   │
        ┌──────────────────────────┼──────────────────────────┐
        ▼                          ▼                          ▼
┌──────────────────┐    ┌──────────────────┐      ┌──────────────────────┐
│ CLUSTER A        │    │ CLUSTER B        │      │ CLUSTER C            │
│ Structure/pocket │    │ ADMET / safety   │      │ Mechanism / KG       │
│                  │    │                  │      │                      │
│ ESM2-650M cache  │    │ ADMET-AI (41 EP) │      │ PrimeKG subgraph     │
│ Boltz-2 struct   │    │   - BBB          │      │ TxGNN zero-shot      │
│ Boltzina aff     │    │   - hERG, DILI   │      │ Cognition virtual    │
│ (top-N pairs)    │    │   - P-gp, CYPs   │      │ phenotype anchor     │
└────────┬─────────┘    └────────┬─────────┘      └──────────┬───────────┘
         │                       │                            │
         │              ┌────────▼────────┐                   │
         │              │ HARD GATES      │                   │
         │              │ Physical kill   │                   │
         │              │ criteria        │                   │
         │              └────────┬────────┘                   │
         │                       │ (drop ~30-50%)             │
         └────────┬──────────────┼────────────────────────────┘
                  ▼              ▼
            ┌───────────────────────────────┐
            │  FUSION LAYER                 │
            │  Stage 1: RRF (k=60)          │  ← ship this for v2
            │  Stage 2: LambdaMART          │  ← gate on ≥20 labels
            └────────────────┬──────────────┘
                             ▼
            ┌───────────────────────────────┐
            │  V2 OUTPUTS                   │
            │  - ranked candidates parquet  │
            │  - provenance (per-cluster)   │
            │  - disagreement diagnosis     │
            │  - funnel narrative           │
            └───────────────────────────────┘
```

## Non-negotiable design rules (from the research doc)

1. **ADMET hard gates apply BEFORE rank fusion**. They are physical kill criteria, not relative rankings. Position: between Cluster B and Fusion.
2. **Never load Boltz-2 and ESM2 simultaneously** on the 12 GB RTX 5070. Sequence loads with `torch.cuda.empty_cache()` + explicit `del model; gc.collect()` between.
3. **Boltzina mode by default** (Furui & Ohue 2025, arXiv 2508.17555) — pose from AutoDock Vina + Boltz-2 affinity head. ~11.8× faster than full Boltz-2 with accuracy below full-Boltz but above Vina/GNINA.
4. **RRF first, LambdaMART later**. Cormack et al. 2009 SIGIR showed RRF beats Condorcet + learning-to-rank under sparse-label regimes. Gate the LambdaMART promotion on ≥20 labeled positives.
5. **Score distributions are heterogeneous** (signed real, [0,1] prob, ordinal). RRF is rank-based and immune. Avoid CombSUM without Platt scaling.
6. **Cache every expensive computation** as parquet keyed by content hash. The pipeline must be cheap to re-run.
7. **Defer Cluster D (transcriptomic)**. LINCS L1000 is dominated by cancer cell lines; the cognition evidence base is too thin to bear primary-rank weight. Only add if RRF is converged AND a neuron-relevant L1000 dataset exists.

## V2 module layout (additive, v1 preserved)

```
src/mammal_repurposing/
├── config.py                    (existing, unchanged)
├── _compat.py                   (existing, unchanged)
├── fetchers/                    (existing, unchanged)
├── scoring/                     (existing, unchanged — used by Cluster A as MAMMAL DTI)
├── analysis/                    (existing, unchanged)
├── cli.py                       (existing, extended with v2 subcommands)
│
├── cluster_a/                   (NEW)
│   ├── __init__.py
│   ├── esm2_embed.py            (ESM2-650M target embedding cacher)
│   ├── boltz_runner.py          (Boltz-2 structure prediction)
│   └── boltzina.py              (Boltzina affinity-only mode)
│
├── cluster_b/                   (NEW)
│   ├── __init__.py
│   └── admet_ai_runner.py       (41-endpoint ADMET-AI prediction)
│
├── cluster_c/                   (NEW)
│   ├── __init__.py
│   ├── primekg.py               (PrimeKG loader + subgraph extractor)
│   ├── txgnn.py                 (TxGNN zero-shot indication scoring)
│   └── cognition_anchor.py      (Virtual phenotype anchor: MCI ∪ AD ∪ ADHD ∪ FXS ∪ narcolepsy)
│
├── gates/                       (NEW)
│   ├── __init__.py
│   ├── admet_gates.py           (BBB / hERG / P-gp / DILI / Ames thresholds)
│   └── disagreement.py          (MAMMAL vs Boltzina disagreement protocol)
│
├── fusion/                      (NEW)
│   ├── __init__.py
│   ├── rrf.py                   (Reciprocal Rank Fusion — ship this)
│   ├── lambdamart.py            (LightGBM LambdaMART — promotion path)
│   └── calibration.py           (Platt scaling for CombSUM if ever needed)
│
├── provenance/                  (NEW — creative addition)
│   ├── __init__.py
│   ├── tracker.py               (per-candidate cluster contribution log)
│   ├── disagreement_report.py   (renders model-disagreement diagnosis)
│   └── narrative.py             (prose funnel explanation for top-N)
│
└── pipeline/                    (NEW)
    ├── __init__.py
    ├── run_phase0_cache.py      (ESM2 embed + Boltz struct, one-time)
    ├── run_phase1_fast.py       (MAMMAL DTI + ADMET-AI + TxGNN parallel)
    ├── run_phase2_boltzina.py   (Boltzina affinity, top-N expensive)
    └── run_phase3_fusion.py     (RRF + gates + provenance + narrative)

configs/                         (NEW)
├── thresholds.yaml              (ADMET hard-gate cutoffs)
└── weights.yaml                 (RRF k, KG_score component weights, ADMET_score weights)

data/
├── cache/                       (NEW — content-hash-keyed memoization)
│   ├── esm2/<sha1(seq)>.pt
│   ├── boltz_struct/<sha1(seq)>.cif
│   ├── boltzina/<sha1(seq)+sha1(smi)>.json
│   ├── admet/<sha1(smi)>.parquet
│   └── txgnn/<sha1(smi)>.json
├── kg/primekg/                  (NEW — large download, gitignored)
└── results/
    ├── v1/                      (existing artifacts mirrored under v1/)
    └── v2/                      (NEW)
        ├── admet_predictions.parquet
        ├── admet_gates_passed.parquet
        ├── boltzina_affinity.parquet
        ├── kg_scores.parquet
        ├── kg_paths.parquet
        ├── rrf_ranking.parquet
        ├── final_ranking.parquet
        ├── provenance.parquet
        ├── disagreement_report.md
        └── funnel_narrative.md
```

## Cognition-specific configuration

### `configs/thresholds.yaml` — hard ADMET gates (cut compound below/above)

| Endpoint | Threshold | Direction | Action |
|---|---|---|---|
| BBB penetration | 0.30 | < | CUT (must reach CNS) |
| P-gp substrate  | 0.85 | > | CUT (efflux kills CNS exposure) |
| hERG inhibition | 0.70 | > | CUT (cardiotox; chronic dosing) |
| DILI            | 0.80 | > | CUT (hepatotox; chronic) |
| Ames mutagen.   | 0.85 | > | CUT (hard kill) |
| CYP3A4 inhib.   | 0.85 | > | FLAG (DDI risk) |
| Caco-2 logPapp  | -5.5 | < | FLAG (oral bioavailability) |

### `configs/weights.yaml`

```yaml
fusion:
  rrf_k: 60                       # Cormack et al. 2009 default
  lambdamart_min_labels: 20       # promote only when this many positives exist
  cluster_rrf_weights:            # equal by default; sweep if calibration data exists
    cluster_a_mammal: 1.0
    cluster_a_boltzina: 1.0
    cluster_b_admet: 0.5          # ADMET enters fusion as a soft signal too (gates run separately)
    cluster_c_txgnn: 1.0
    cluster_c_kg_paths: 0.5

admet_score:
  bbb: 0.35
  one_minus_herg: 0.20
  one_minus_pgp: 0.15
  one_minus_dili: 0.10
  caco2_norm: 0.10
  one_minus_cyp3a4: 0.10

kg_score:
  indication: 0.4
  contraindication: 0.3
  path_count_log: 0.2
  side_effect_penalty: 0.1

disagreement_thresholds:
  mammal_pkd_strong: 7.0          # log10 Kd, ~100 nM
  boltzina_binder_prob_low: 0.30
  boltzina_binder_prob_high: 0.70
  structure_plddt_high_confidence: 70
```

## Phase plan with VRAM sequencing (per research doc §4)

| Phase | When | Operations | Peak VRAM | Wall-clock |
|---|---|---|---|---|
| 0 | One-time, cached | ESM2-650M → 22 target embeddings → cache; Boltz-2 → 7 missing target structures → cache | 2.5 GB (ESM) → 8-10 GB (Boltz) | 5 min + 30-60 min |
| 1 | Per-batch | MAMMAL DTI grid + ADMET-AI + TxGNN (CPU/light GPU) + apply ADMET hard gates | 3 GB (MAMMAL) + ~2 GB (TxGNN) | 10 min + 2 min + 5 min |
| 2 | Per top-N surviving compounds | Boltzina affinity for top-50 surviving × 22 targets ≈ 1100 calls | 7-8 GB | ~6 hours (cold cache) |
| 3 | Final | RRF over (MAMMAL, Boltzina, TxGNN, ADMET_score) per pair → final ranking + provenance + sanity gate | <1 GB | seconds |

## Cognition-specific virtual phenotype anchor (Cluster C)

Because "healthy cognitive enhancement" is not a node in PrimeKG, we anchor the KG query to the *union* of 2-hop neighborhoods around these disease nodes:

| Disease | EFO / MeSH | Rationale for inclusion |
|---|---|---|
| Mild cognitive impairment | EFO_0006816 | Direct cognitive-decline proxy |
| Alzheimer's disease | EFO_0000249 | Mechanistic proxy (cholinergic, AMPA, etc.); not clinical |
| ADHD | EFO_0003888 | Processing-speed / working-memory targets |
| Fragile X syndrome | EFO_0004247 | Anchor for PDE4D (BPN14770 evidence) |
| Narcolepsy | EFO_0003781 | Anchor for HCRTR1/2 + HRH3 |

`KG_score(compound) = w_ind * mean(TxGNN_indication over 5 anchors)
                    - w_con * mean(TxGNN_contraindication)
                    + w_path * log(1 + PrimeKG_path_count(compound → any panel target))
                    - w_se * side_effect_overlap(compound, {sedation, somnolence, cognitive impairment, anticholinergic})`

## What this addresses vs. the v1 measured failures

| V1 failure mode | V2 mitigation |
|---|---|
| Peptide pollution (liraglutide/semaglutide@top) | ADMET gate: peptides fail BBB filter (already empirically true — semaglutide p_bbb=1.2e-7 in v1 aux scores) |
| Structural-bias inflation (atorvastatin/fexofenadine@top) | Boltzina pose feasibility; ADMET hERG/DILI gates remove some; TxGNN contraindication signal for known-peripheral drugs |
| CHRNA7 allosteric blindness (std=0.029) | Boltz-2 + Boltzina at α7 nAChR's known PAM site (PDB 7EKT) gives a structure-aware second opinion |
| Compressed dynamic range at most targets | RRF over multiple rankers naturally widens spread; Boltzina contributes signed values from a different distribution |
| Calibration unknown without ChEMBL ground truth | Per-cluster, per-target rho computed by `provenance/disagreement_report.py` |

## Creative additions (mine, not in the research doc)

These complement the doc's prescription without contradicting it.

### 1. Provenance tracker (`provenance/tracker.py`)
For every compound surviving to the final list, store a structured record of *which cluster placed it where*: MAMMAL rank, ADMET gate status (passed/flagged/cut), TxGNN rank, Boltzina rank if available. Drives the disagreement-diagnosis report + funnel narrative. Output: `provenance.parquet`.

### 2. Disagreement-diagnosis report (`provenance/disagreement_report.py`)
The research doc specifies a per-pair model-disagreement protocol. I render it as a *visible artifact* — one section per disagreement archetype:
- `mammal_strong_boltzina_weak`: MAMMAL pKd > 7, Boltzina binder_prob < 0.3 → "downrank candidate, likely sequence-only false positive"
- `boltzina_strong_mammal_weak`: symmetric, structure-aware win
- `txgnn_strong_no_panel_hit`: off-panel mechanism, flag for follow-up
- `admet_gate_failed_other_strong`: do not bypass; physical kill criterion

This is a publishable methodology contribution per the doc's recommendation §8.3.

### 3. Funnel narrative (`provenance/narrative.py`)
For each top-N candidate, generate a short prose explanation (template, not LLM):
> "Compound X reached the shortlist via cluster A (MAMMAL pKd 6.8 at HRH3, rank 12/271); cluster B passed all hard gates (BBB 0.99, hERG 0.03, DILI 0.12); cluster C TxGNN indication probability 0.42 against the narcolepsy anchor; Boltzina affinity_pred -1.2 (binder_prob 0.71). Three of four clusters voted strong; recommended next: radioligand binding at HRH3."

Optional LLM-rewrite layer using local model in Phase 5 if cycles allow.

### 4. Dry-run mode
Each pipeline phase has a `--dry-run` flag that emits stubbed cluster outputs with realistic shapes. Lets you wire-test the full pipeline before installing Boltz-2 (5+ GB) or PrimeKG (multi-GB).

## Execution order (matches research doc §8 with concrete check)

1. **Week 1 (NOW)**: Cluster B (ADMET-AI) + RRF + gates. Validation: positive controls pass BBB. ← _highest signal/effort, no GPU contention_
2. **Week 2-3**: Cluster A (ESM2 + Boltzina). Validation: known PAMs (PNU-120596, EVP-6124) rank top-20% at their targets.
3. **Week 3-4**: Cluster C (PrimeKG + TxGNN). Validation: donepezil/memantine/BPN14770 top-decile against cognition anchor.
4. **Week 5**: LambdaMART promotion *only if* ≥20 labels exist. Else stay on RRF and ship.
5. **Week 6+**: Optional Cluster D / LLM explanation layer.

## What does NOT enter v2

- Cluster D (transcriptomic / L1000) — defer per doc §2
- LLM literature reasoning as evidence source — defer per doc §3, §7
- Embedding alignment between modalities (BioBridge-style) — rejected per doc §1
- Train LambdaMART before ≥20 labels — rejected per doc §5
- "Smart drug discovery" framing — the publication framing is methodology + benchmark + enrichment, not discovery (doc §8)
