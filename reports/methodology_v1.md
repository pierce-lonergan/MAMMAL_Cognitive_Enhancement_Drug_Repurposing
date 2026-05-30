# Methodology Note v1 — MAMMAL Cognitive Enhancement Drug Repurposing

**Repository:** `MAMMAL_Cognitive_Enhancement_Drug_Repurposing` (private)
**Pipeline snapshot:** main @ `530dc40` (diagnostics + Tanimoto ranker)
**Author:** Pierce Lonergan, with Claude-assisted scaffolding
**Date:** 2026-05-25

This is the v1 methodology note. Its purpose is to describe what the pipeline does, what it does NOT do, and where to find the artifacts. It is not a publication and is not a wet-lab handoff. The point of the note is to make the system's failure modes legible so that downstream consumers (us, future agents, collaborators) do not over-interpret its rankings.

> **⚠️ POST-SHIP UPDATE (commit `530dc40`)**: A diagnostic protocol against
> `research/4-tier/Diagnosing MAMMAL DTI Anti-Correlation.md` revealed
> **MAMMAL prior collapse is panel-wide** (19/22 targets at >10× collapse vs
> training SD) and a **1996-vintage Tanimoto-on-Morgan-FP baseline beats
> MAMMAL at every audited target** (SLC6A3 +0.90 vs -0.70; DRD1 +0.85 vs +0.29).
> The Phase A.7 per-target ρ values are statistics on noise for most compounds.
> A new Cluster A.4 (Tanimoto-to-ChEMBL-actives ranker) has been added to fusion
> as the immediate v4 fix. See `reports/diagnostics_v1.md`,
> `reports/tanimoto_baseline_v1.md`, `reports/fusion_tanimoto_addition_diff.md`.

---

## 1. Question

> Given a panel of 22 cognition-relevant human protein targets and a curated library of ~300 nootropic or putatively-pro-cognitive compounds, can we produce a defensible *rank ordering* over compounds that integrates (a) sequence-only DTI prediction, (b) structure-aware affinity prediction, (c) drug-likeness/ADMET, and (d) knowledge-graph indication signal — and *honestly weight each cluster against ground truth at each target*?

The pipeline does **not** answer: "is compound X likely to enhance working memory in humans?" That requires behavioural data MAMMAL has never seen. We answer the binding-affinity-and-context proxy question only.

---

## 2. Four clusters

| Cluster | Identity | Output per (target, compound) | Source |
|---|---|---|---|
| **A.1** | MAMMAL DTI (IBM, 458M-param T5-style, BindingDB-trained) | predicted pKd | `scripts/04_score_dti.py` → `dti_scores.parquet` |
| **A.2** | Boltzina (Boltz-2 affinity-only, structure-aware) | log10 IC50 (µM) + binder probability | `scripts/_wsl2_boltz_full_sweep.py` → `boltzina_affinity.parquet` |
| **B** | ADMET-AI (41 endpoints, Caco-2, BBB, hERG, DILI, etc.) | per-compound `admet_score` + gate verdict | `scripts/14_v2_cluster_b_admet.py` → `admet_gates.parquet` |
| **C** | PrimeKG (~129K nodes, ~4M edges, igraph) + TxGNN (zero-shot indication) | per-compound knowledge-graph score | `scripts/23_v3_cluster_c.py` → `kg_scores.parquet` *(WSL2 venv)* |

Each cluster is run independently and stored in a versioned parquet. Fusion is downstream and reversible.

### Why these four

* MAMMAL alone is a single-cluster prior — fast, but has a documented dynamic-range collapse on allosteric pockets ([T1 audit](../design/decision-records/T1_CHEMBL_AUDIT_VERDICT.md) confirmed the target IDs were correct, ruling out the easy fix; [Phase 0.5](../design/decision-records/PHASE_0_5_DECISION_RECORD.md) confirmed Boltz-2 rescues CHRNA7 PAMs).
* Boltzina gives a structure-conditioned second opinion at the cost of ~80 s per pair on RTX 5070 in WSL2 (cuequivariance kernels). It cannot replace MAMMAL for throughput, only complement it.
* ADMET-AI is a regulatory filter, not a ranker — `gate_status ∈ {PASS, FLAG, CUT}` plus a composite `admet_score`. Approved drugs get a `regulatory_bypass` flag so an FDA history doesn't drown in toxicity gates the molecule has empirically already cleared.
* PrimeKG + TxGNN provide an orthogonal "does this compound show up in published indication graphs for cognition-adjacent diseases" signal. The cognition virtual phenotype anchor is a weighted union of 5 EFO disease IDs (MCI, AD, ADHD, FXS, narcolepsy — weights in `configs/weights.yaml`).

---

## 3. Calibration — the linchpin

This is the thing that distinguishes v3 from "vote everyone equally." The script lives at `scripts/22_v3_calibration.py`.

### What it does

For each panel target:

1. Query the local ChEMBL 36 SQLite mirror (`fetchers/chembl_sqlite.py`) for every (molecule, target) record with `assay_type='B'`, `confidence_score ≥ 7`, `standard_type ∈ {Ki, IC50, Kd, EC50}`, `pchembl_value IS NOT NULL`. Salt forms are walked to parents via `molecule_hierarchy`.
2. Join MAMMAL's predicted pKd at that target to the truth set by InChIKey (RDKit). Compute Spearman ρ.
3. Same for Boltzina's `affinity_pred_value` (negated to align direction with pchembl).
4. Verdict — refined from the original sprint spec to handle negative ρ and single-cluster cases:

   | Verdict | Condition | Weight effect |
   |---|---|---|
   | `MAMMAL_2X_BOLTZ`     | both ≥3 preds; ρ_M ≥ 0.30 AND ρ_M ≥ ρ_B + 0.1 | MAMMAL ×2.0, Boltz ×1.0 |
   | `BOLTZ_2X_MAMMAL`     | both ≥3 preds; ρ_B ≥ 0.30 AND ρ_B ≥ ρ_M + 0.1 | Boltz ×2.0, MAMMAL ×1.0 |
   | `EQUAL_WEIGHTS`       | both ≥0.30, within 0.1 | no override |
   | `DE_WEIGHT_TARGET`    | both <0.30 (including negative) | both ×0.3 |
   | `MAMMAL_ONLY_STRONG`  | only MAMMAL has data; ρ_M ≥ 0.30 | no override |
   | `MAMMAL_ONLY_WEAK`    | only MAMMAL has data; 0 ≤ ρ_M < 0.30 | MAMMAL ×0.6 |
   | `MAMMAL_ONLY_INVERTED`| only MAMMAL has data; ρ_M ≤ -0.30 | MAMMAL ×0.3 + manual-review flag |
   | `BOLTZ_ONLY_*`        | symmetric (rare; will become common after overnight sweep) | as above |
   | `INSUFFICIENT_DATA` / `NO_CLUSTER_DATA` | truth < 5 or no joined preds | no override |

   We deliberately do NOT auto-invert MAMMAL when ρ < -0.30. The sign could flip with more data and we don't want a non-monotonic policy.
5. Write `reports/calibration_report.md` (legible) and `configs/weights_calibrated.yaml` (machine-readable per-target overrides).

### What we found at v1 (this snapshot)

`reports/calibration_report.md` summary:

* **1 BOLTZ_2X_MAMMAL** — HRH3 (ρ_B = +0.87, n=3 vs ρ_M = -0.14, n=12). Small Boltz sample; fragile but directionally encouraging.
* **2 MAMMAL_ONLY_STRONG** — DRD1 (ρ=+0.31, n=21), HCRTR1 (ρ=+0.37, n=6). The only two targets where MAMMAL's ranking is materially positively correlated with experimental pchembl.
* **14 MAMMAL_ONLY_WEAK** — including ACHE (ρ=+0.20, n=10) — the donepezil/rivastigmine target — and the GRIA1-4 family. MAMMAL ranks at these targets are weakly informative.
* **4 MAMMAL_ONLY_INVERTED** — **SLC6A3 (DAT, ρ=-0.71, n=26)**, **SLC6A2 (NET, ρ=-0.53, n=25)**, **GRIN2A (ρ=-0.35, n=8)**, **GRIN2B (ρ=-0.30, n=14)**. MAMMAL's top-ranked compounds at these targets are systematically the *worst* binders, not the best.
* **1 NO_CLUSTER_DATA** — KCNQ3 (only 2 ChEMBL InChIKeys joined to our 298 compounds).

### What this means in practice

The pipeline still ships rankings, but the calibrated weights down-weight 18 of 22 targets and de-weight to 0.3 for the 4 inverted ones. That means:

* In the calibrated 4-cluster fusion, HRH3 (the only un-down-weighted target with a positive Boltz signal) dominates many compound rankings — top gainers under calibration are all HRH3-best-target compounds.
* The previously-passing positive-control sanity gate (donepezil top-5% at ACHE etc.) is now legible as "donepezil's high rank at ACHE was a *rank percentile* property of the curated library, not a property of MAMMAL ranking ACHE binders correctly across all of ChEMBL." Both statements can be true; they answer different questions.
* The inversion at the monoamine transporters is the headline negative finding. The likely mechanism is that BindingDB has dense coverage of high-affinity DAT/NET inhibitors (the cocaine/methylphenidate chemical space) but MAMMAL's per-protein representation may not separate them by potency — it ranks chemically similar molecules without distinguishing magnitudes. Re-running with Boltz at these targets (post overnight sweep) is the test of whether a structure-conditioned model can recover the correct ordering.

---

## 4. Fusion — RRF with per-target weights

`src/mammal_repurposing/fusion/rrf.py` implements Reciprocal Rank Fusion (Cormack, Clarke & Buettcher, SIGIR 2009):

$$
\text{RRF}(c) = \sum_{r \in \text{rankers}} \frac{w_r}{k + \text{rank}_r(c)}, \quad k = 60
$$

where `rank_r(c)` is the 1-indexed rank of compound `c` by ranker `r`. Per-target overrides from `configs/weights_calibrated.yaml` replace `w_r` at that target.

### Why RRF, not CombSUM

Boltzina log-IC50 (signed real), MAMMAL pKd (positive real), ADMET-AI [0,1], TxGNN [0,1] live on incompatible scales. CombSUM would require Platt scaling to make them comparable; RRF needs none. The cost is that absolute score magnitudes lose meaning — RRF outputs are unit-less and only the within-list ordering is preserved.

### Per-target then aggregate

`rrf_per_target_then_compound` runs RRF inside each target (so each target votes independently with its own per-cluster weights), then sums per-target RRF scores per compound. A compound that ranks well across many targets gets credit (polypharmacology). A compound that wins one target by a huge margin doesn't dominate.

### Output

`scripts/15_v2_fusion.py` writes both passes (uncalibrated + calibrated) under `data/results/v2/`:

```
rrf_ranking_uncalibrated.parquet      rrf_ranking_calibrated.parquet
final_ranking_uncalibrated.parquet    final_ranking_calibrated.parquet
disagreement_report_uncalibrated.md   disagreement_report_calibrated.md
funnel_narrative_uncalibrated.md      funnel_narrative_calibrated.md
```

Phase D's diff (`scripts/25_v3_fusion_diff.py` → `reports/fusion_calibration_diff.md`) shows: Spearman ρ between the two orderings is **+0.994**. Calibration preserves macro-order but materially re-shuffles compounds whose top-target was inverted or weak.

---

## 5. Provenance + reproducibility

Every parquet has timestamps. Every report cites the script that produced it. Major artifacts:

* `data/results/dti_scores.parquet` (6,556 rows = 22 targets × 298 compounds) — MAMMAL
* `data/results/v2/boltzina_affinity.parquet` (growing — overnight sweep in progress; 92 rows at this snapshot)
* `data/results/v2/admet_gates.parquet` (298 rows)
* `data/results/v2/kg_scores.parquet` — pending; produced by `scripts/23_v3_cluster_c.py` in the WSL2 `txgnn_env` venv (separate from `mammal_env` due to PyG ABI incompatibility with PyTorch 2.12 nightly — see `scripts/_wsl2_setup_cluster_c.sh` header for the why)
* `reports/calibration_report.md` — Phase A.7 verdicts per target
* `reports/fusion_calibration_diff.md` — Phase D before/after
* `configs/weights.yaml` + `configs/weights_calibrated.yaml` — global + per-target overrides
* `reports/chembl_target_id_audit_sqlite.md` — Phase A.6: 21/22 ChEMBL target picks ALIGNED (and the original 13 NO_RECORDS were all REST timeouts, not real misses)
* `reports/sqlite_vs_rest_smoke.md` — Phase A.5 gate (19/20 agreement; 1 was a REST 500, no real disagreement)

The local ChEMBL 36 SQLite mirror (`~/.data/chembl/36/chembl_36.db`, ~13 GB extracted) is regenerated via:
```bash
python -c "import chembl_downloader; chembl_downloader.download_extract_sqlite()"
```
This is idempotent. The salt-form parent resolution lives in `_LOOKUP_BY_INCHIKEY_SQL`'s `canonical`/`expanded` CTEs and uses `LEFT JOIN molecule_hierarchy mh ON md.molregno = mh.molregno`. Skipping this loses ~30% of records on common salt-form compounds.

---

## 6. Known limitations

1. **Boltz coverage is currently 1/22 targets with n≥3.** The overnight WSL2 sweep (1,165 pairs queued) is at ~7% at this snapshot. Almost every Phase A.7 verdict is `MAMMAL_ONLY_*` rather than the 2-cluster verdicts the calibration was designed to deliver. Re-run Phase A.7 → Phase C → Phase D when the sweep completes. Expect the inverted-MAMMAL targets to either flip to `DE_WEIGHT_TARGET` (both clusters poor) or `BOLTZ_2X_MAMMAL` (structure rescues sequence) — that distinction is the most important downstream signal we don't yet have.

2. **The truth set is ChEMBL only.** We use pchembl_value (a unified log of Ki/IC50/Kd/EC50 in mol). Assay heterogeneity is partially controlled by `confidence_score ≥ 7` and `assay_type = 'B'` (binding), but inter-assay variance at a single target can still be ±0.5–1.0 log units. Spearman ρ is robust to monotone scale issues but not to true assay disagreement.

3. **MAMMAL_ONLY_INVERTED at DAT/NET is not yet explained.** It could be (a) BindingDB sampling bias toward a narrow chemical space at these transporters, (b) tokenisation artefacts in the protein sequence near the substrate-binding pocket, (c) a real model bug, or (d) ChEMBL pchembl having different chemical-space coverage than BindingDB (so the comparison isn't apples-to-apples). Diagnosing this is on the v2 roadmap.

4. **No wet-lab confirmation.** The wet-lab shortlist in `reports/wet_lab_shortlist.md` is a *prioritisation* output, not a recommendation. Every entry has the verdict, the contributing clusters, the gate status, and the ChEMBL backstop verdict so a human can decide whether to send it to assay.

5. **PrimeKG + TxGNN (Cluster C) is implemented but not yet run on this snapshot.** WSL2 setup script is `scripts/_wsl2_setup_cluster_c.sh`; download is `scripts/_wsl2_download_primekg.sh`. The fusion will accept the parquet the moment Cluster C runs; no upstream change required.

6. **RRF k_const = 60 is the Cormack default.** Sensitivity in [30, 80] is documented as small but we have not swept it. Per-cluster weights are static within a target; we have not tried compound-class-conditioned weights.

7. **LambdaMART is not yet eligible.** `src/mammal_repurposing/fusion/lambdamart.py` exists as a stub. Per the v2 design doc it should not be promoted past RRF until we have ≥20 labelled (compound, target, true-binder?) tuples. We have ~2,500 ChEMBL truth rows on the surviving 298 compounds, but the labelling is heterogeneous (pchembl is a regression target, not a binary). A first cut would binarise at pchembl ≥ 6.0 and try LambdaRanker on the 22 per-target lists; the v1 calibrated RRF is sufficient until then.

---

## 7. Sprint history (the path to this snapshot)

| Commit | What |
|---|---|
| `9f800f8` | Phase D: calibrated vs uncalibrated diff (this snapshot) |
| `b716ec7` | Phase C: 4-cluster RRF wiring + per-target overrides |
| `7c1d55e` | Phase A.6 + A.7: SQLite audit and calibration linchpin |
| `1551a3e` | Phase A.5: SQLite vs REST agreement (19/20) + schema fix |
| `9d40e19` | Phase A.1-A.4: chembl_sqlite module + calibration scripts |
| `42b6597` | Cluster C: real PrimeKG + TxGNN + WSL2 isolated venv |
| `7467ac1` | Phase 0.5: Boltz-2 rescues CHRNA7 PAMs (gate passed) |
| `ed4761c` | T1 audit verdict (hypothesis falsified — picks were correct) |
| `3baf422` | V3 attack plan |
| `6d9aa91` | v2 hybrid architecture skeleton |
| `627c528` | v1 MVP (Phases 0-5) |

---

## 8. What the next agent should do first

1. Check the WSL2 Boltz sweep is alive (`bash scripts/_wsl2_sweep_status.sh` from a WSL2 shell, or `wsl -d Ubuntu -- tail /tmp/wsl2_boltz_sweep.log`).
2. When the sweep completes (currently ETA ~24 h from this snapshot): re-run Phase A.7 (`scripts/22_v3_calibration.py`), then Phase C (both fusion passes), then Phase D (`scripts/25_v3_fusion_diff.py`). Diff the new `calibration_report.md` against this snapshot's — that diff is the answer to "did Boltz actually rescue the inverted targets?"
3. Run Cluster C in the WSL2 `txgnn_env` venv. Once `kg_scores.parquet` exists, the fusion absorbs it without code changes.
4. If LambdaMART becomes eligible (≥20 labels per target with sufficient class balance), promote from RRF.
5. Do not touch the salt-form parent resolution in `chembl_sqlite.py`. Re-confirmed in this sprint via the `molecule_hierarchy` JOIN fix — bypassing it loses ~30% of records.

---

_End of methodology note v1._
