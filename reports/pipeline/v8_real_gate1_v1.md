# V8 Gate 1 on real LINCS L1000 (v1)

The pre-registered V8 primary gate, run on REAL LINCS Level-5 signatures of cognition compounds (not synthetic). Each compound is labelled by its pharmacology (target / mechanism), never by its signature, so the test is not circular. Signatures are collapsed to one consensus per compound, PCA-reduced, clustered, and scored against the pharmacology labels with Adjusted Mutual Information.

- Compounds with real LINCS signatures: **16** of 31 labelled
- Mechanism classes represented: **10**
- Signatures loaded (pre-consensus): see log; consensus = 1 per compound

## Result

| Set | Method | n compounds | n classes | AMI | ARI | V-measure | Verdict |
|---|---|---|---|---|---|---|---|
| all_classes | agglomerative | 16 | 10 | 0.090 | 0.078 | 0.781 | **FAIL** |
| all_classes | hdbscan min=2 | 16 | 10 | 0.000 | 0.000 | 0.000 | **FAIL** |
| multi_member_classes | agglomerative | 9 | 3 | 0.126 | 0.071 | 0.393 | **FAIL** |
| multi_member_classes | hdbscan min=2 | 9 | 3 | 0.000 | 0.000 | 0.000 | **FAIL** |

**Best AMI = 0.126 (multi_member_classes, agglomerative) -> FAIL**

## Gate 1 bands (V8 OSF pre-reg section 5.1)

| Band | AMI | ARI | Action |
|---|---|---|---|
| PASS | >= 0.50 | >= 0.40 | enter joint posterior at lambda_phen = 1.0 |
| DEGRADE | [0.30, 0.50) | [0.25, 0.40) | enter at lambda_phen = 0.5; flag |
| FAIL | < 0.30 | < 0.25 | publish the negative result |

## Honest reading

- Real L1000 signatures do NOT recover the cognition mechanism classes at the pre-registered bar. This is a pre-registered negative result for the phenotype axis on this compound set, reported honestly rather than buried. It is consistent with the known difficulty of recovering fine pharmacology from bulk transcriptional response (cell-line and assay variance compete with mechanism).
- Labels are pharmacology-grounded, so the score is not inflated by circular signature-derived labels.
- GSE70138 predates several newer compounds (zatolmilast, encenicline, blarcamesine), so the labelled set is whatever is actually present in the 2017 release; the count above is the real coverage, not a target. With 16 compounds across 10 classes the test is underpowered, which is itself part of the honest finding.
- This real result SUPERSEDES the synthetic dry-run (`v8_gate1_dryrun_v1.md`, AMI = 1.00). That dry-run used orthogonal centroids that are trivially separable and was only ever a pipeline sanity check; it does not represent real performance. The number to cite for V8 mechanism-class recovery is the one on this page.
- The FAIL is consistent with this project's central finding: the target-centric and phenotype-centric paradigms (binding affinity, target genetics, and now bulk perturbational signature) all underperform the mechanism-class clinical track record. V8 is honestly a weak axis on the available real data, not a validated one. This is a single-modality (LINCS-only) test; the full V8 MOFA+ multi-modal stack is unrun and would need its own gate.
