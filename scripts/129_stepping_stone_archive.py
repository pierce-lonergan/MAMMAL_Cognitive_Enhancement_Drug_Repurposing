"""Build the stepping-stone archive: every hypothesis this project killed, and what would revive it.

Before this script the repository had thirty-eight ledgers and not one of them held a failure. Every
kill lived as prose in a markdown file, which is to say it lived somewhere no program could reach,
which is to say it was permanent by accident rather than by evidence.

The archive has two sources and they are deliberately different:

  DERIVED    the fourteen compounds labelled `enhances_healthy_young = 0` are classified straight
             from the live ledger by `classify_null`, so the archive cannot drift from the data it
             describes. Re-running this script re-derives them.
  CURATED    the architectural kills (a pre-registered test, a retracted prediction, a demoted
             screen, an unverifiable source) are hand-written here, each carrying the measured
             number that killed it and the artifact that recorded it.

On the derived half, the classification reproduces a finding the project already made and then did
nothing with. `healthy_adult_robustness_v1.md` section R3 says in as many words that most of these
nulls are not refuted but under-powered. That was written, published, and left there: the ledger
still carries a binary label, and every downstream consumer still reads fourteen refutations. A
finding that changes no data structure is a finding that dies. This archive is the data structure.

Writes data/raw/stepping_stone_archive.csv. CPU only, no network.
"""
from __future__ import annotations

import csv
import logging
import sys
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
L = logging.getLogger("stepping_stones")
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mammal_repurposing.archive.stepping_stone import (  # noqa: E402
    FAILURE_MODES, REQUIRED_COLUMNS, TARGET_EFFECT_G, classify_null, validate_archive,
)

LEDGER = ROOT / "data" / "raw" / "healthy_adult_cognition_ledger.csv"
OUT = ROOT / "data" / "raw" / "stepping_stone_archive.csv"


# -------------------------------------------------------------------------------------------------
# CURATED: the architectural kills. Every `evidence` field is a number this repository measured.
# -------------------------------------------------------------------------------------------------
CURATED = [
    dict(
        hypothesis_id="H-B1-pairing",
        claim="Durable cognitive gain is a drug-by-experience interaction: pairing a compound with "
              "training produces post-washout gain that the compound alone does not.",
        domain="PERSEUS/durability",
        verdict="KILLED", died_on="2026-07-28",
        killed_by="scripts/124_paired_experience_contrast.py -> paired_experience_contrast_v1.md",
        evidence="n-weighted contrast (paired - unpaired) = -0.071, permutation p = 0.2692, 15 "
                 "studies; at n >= 25 only: -0.076, p = 0.3345, 6 studies. Point estimate points "
                 "the WRONG WAY. Largest n among any positive-direction study = 21.",
        failure_mode="underpowered",
        keystone="Six studies at n >= 25 cannot resolve a between-study contrast. The pattern that "
                 "decided it (every positive came from a small study) is the signature of a "
                 "small-study effect, which is a reason to distrust the positives AND a reason the "
                 "contrast itself is unresolved. The hypothesis was not refuted; the test could "
                 "not see.",
        keystone_predicate="paired_studies_at_n(25, 15)",
        revival_test="Re-run scripts/124 once 15 or more paired studies at n >= 25 are curated. "
                     "Pre-register the contrast direction again before looking.",
        status="DEAD",
    ),
    dict(
        hypothesis_id="H-L4-window-verdict",
        claim="The L4 structural psychoplastogen screen identifies compounds that open a plasticity "
              "window, and that call is strong enough to carry a persistence verdict on its own.",
        domain="PERSEUS/L4",
        verdict="DEMOTED", died_on="2026-09-15",
        killed_by="scripts/125_window_assay_index.py -> window_assay_index_v1.md",
        evidence="On dendritic_spine, the family the rule was derived from and the only testable "
                 "one: agreement 0.50 at permutation p = 1.000 against 48 curated (compound, "
                 "assay) rows. The family base rate is 5 openers of 6, so a constant 'always "
                 "opener' scores 0.83. In ocular_dominance, perceptual_learning and tms_ltp the "
                 "screen predicts a constant (zero positives).",
        failure_mode="instrument_blind",
        keystone="Two different blindnesses, and only one is a defect. The screen is serotonergic "
                 "by construction, so having no opinion about cholinergic, SSRI, GABAergic and "
                 "ECM-degrading compounds is scope, not error. The defect is that it scores the "
                 "administered molecule when the biology is done by the metabolite, and the "
                 "protection against that is a three-entry hand-written map.",
        keystone_predicate="metabolite_prediction_available()",
        revival_test="Re-run scripts/125 once a general metabolite predictor replaces "
                     "PRODRUG_TO_ACTIVE. The screen may be re-promoted only if it beats the "
                     "permutation gate in at least one family on that footing.",
        status="DEAD",
    ),
    dict(
        hypothesis_id="H-L4-ibogaine",
        claim="Ibogaine is plasticity-window positive (serotonergic scaffold, TPSA 28, HBD 1, "
              "clears the intracellular-access gate).",
        domain="PERSEUS/L4",
        verdict="REFUTED", died_on="2026-09-15",
        killed_by="scripts/125_window_assay_index.py; primary source Ly 2018 Cell Rep",
        evidence="Ly 2018 (PMID 29898390) measured NO structural-plasticity effect for ibogaine "
                 "itself and identified noribogaine as the active species. L4 calls it positive on "
                 "the parent's descriptors. This is L4's only in-scope error in its home family.",
        failure_mode="instrument_blind",
        keystone="Ibogaine is absent from PRODRUG_TO_ACTIVE, so the parent is scored. Psilocybin is "
                 "present, which is exactly why psilocybin is a hit and ibogaine is not. The map's "
                 "coverage IS the protection, and it was written by hand.",
        keystone_predicate="prodrug_resolution_covers(ibogaine)",
        revival_test="Re-score ibogaine through L4 once noribogaine resolution exists; the window "
                     "call should flip to negative for the parent.",
        status="DEAD",
    ),
    dict(
        hypothesis_id="H-B2-signflip-rate",
        claim="The rate at which a compound's plasticity-window sign flips between assay families "
              "is estimable, and is at or above the 30% threshold that would condemn compound-level "
              "window scoring outright.",
        domain="PERSEUS/L4",
        verdict="ABSTAINED", died_on="2026-09-15",
        killed_by="scripts/125_window_assay_index.py -> window_assay_index_v1.md",
        evidence="Measured 25% (1 of 4). Donepezil does flip (perceptual_learning positive, "
                 "Rokem & Silver 2010 PMID 20850321; ocular_dominance negative, Sheynin 2019 "
                 "PMID 30766471, t(11) = -4.9, p < 0.001), so the phenomenon is real - though the "
                 "two results come from DIFFERENT laboratories and different participants, and "
                 "differ in dose regimen (8-day steady state vs single dose), so it is a "
                 "between-study contrast. But only 4 compounds carry a SIGNED "
                 "direction in two or more families, so the RATE is not estimable. The L4 kill came "
                 "from the permutation gate instead, not from this test.",
        failure_mode="underpowered",
        keystone="The field almost never measures two plasticity readouts in the same subjects, "
                 "so the cross-assay question is mostly unanswerable from the published literature "
                 "as indexed, and nearly every apparent flip is a between-study comparison carrying "
                 "between-study confounds. Four compounds is not a rate.",
        keystone_predicate="assay_family_signed_compounds(12)",
        revival_test="Re-run the sign-consistency arm of scripts/125 at 12 or more signed "
                     "multi-family compounds and report the rate with an interval.",
        status="DEAD",
    ),
    dict(
        hypothesis_id="H-B3-ndcg-prediction",
        claim="Adding paired_training_washout_retest to EVIDENCE_RANK would raise ranking NDCG@25 "
              "from 0.8912 to approximately 0.9117.",
        domain="meta/self-prediction",
        verdict="RETRACTED", died_on="2026-07-28",
        killed_by="scripts/47 re-run; recorded in docs/PREREG_DEVIATIONS_2026-06.md",
        evidence="Measured NDCG@25 = 0.8716. The number went DOWN, not up. The hypothesis under "
                 "test still passes; the PREDICTION ABOUT IT was wrong and was retracted in the "
                 "deviations ledger rather than quietly restated.",
        failure_mode="mechanism_contradicted",
        keystone="",
        keystone_predicate="",
        revival_test="",
        status="PERMANENTLY_CLOSED",
    ),
    dict(
        hypothesis_id="H-R1-power-confound",
        claim="The healthy-adult stimulant gate is a power artifact: a predictor knowing only how "
              "heavily a compound was studied matches or beats it.",
        domain="healthy-adult axis",
        verdict="REFUTED", died_on="2026-07-28",
        killed_by="scripts/121_healthy_adult_robustness.py on the expanded n=19 ledger",
        evidence="At n = 11 the n_studies power proxy reached AUROC 0.88. On the expanded n = 19 "
                 "primary set it collapsed to 0.59 (permutation p = 0.2835), while the stimulant "
                 "gate strengthened to p = 0.018 and survived the l-theanine sensitivity "
                 "(p = 0.046, previously 0.176). The confound was real at n = 11 and is not real "
                 "at n = 19.",
        failure_mode="underpowered",
        keystone="This entry is the archive's own cautionary tale: a finding generated at n = 11, "
                 "published, and then refuted by more data. It is recorded as a kill of the "
                 "FINDING, not of the gate.",
        keystone_predicate="ledger_rows_at_least(healthy_adult_cognition_ledger.csv, 60)",
        revival_test="Re-run scripts/121 at a materially larger ledger. A confound that vanished "
                     "between n=11 and n=19 may return; it is not settled either way.",
        status="DEAD",
    ),
    dict(
        hypothesis_id="H-gervain-absolute-pitch",
        claim="Valproate reopened a critical period for absolute-pitch acquisition in healthy adult "
              "men, and is the field's precedent for drug-induced durable human learning.",
        domain="durability precedent",
        verdict="REFUTED", died_on="2026-07-28",
        killed_by="research verification pass -> REJECT_MISDESCRIBED",
        evidence="Gervain 2013 (PMID 24348349): on-drug readout only, arm-1-only analysis, blind "
                 "penetrated in 17 of 18 participants, and DURABILITY WAS NEVER ASSESSED. The "
                 "most-cited precedent for the project's central claim does not test the claim.",
        failure_mode="wrong_endpoint",
        keystone="Nobody has measured whether the learning persists off drug. The study cannot be "
                 "repaired by re-reading it; the measurement was not taken.",
        keystone_predicate="durable_healthy_rows(1)",
        revival_test="If any verified off-drug healthy-adult positive row appears, re-examine "
                     "whether a valproate-plus-training design is the right replication target.",
        status="DEAD",
    ),
    dict(
        hypothesis_id="H-knecht-2004",
        claim="Knecht 2004 provides a levodopa-plus-training human result bearing directly on the "
              "durability verdict, and is the cheapest available falsification of it.",
        domain="durability precedent",
        verdict="ABSTAINED", died_on="2026-07-28",
        killed_by="access failure (publisher returned HTTP 403); recorded in the open-items list",
        evidence="No verdict was reached. The full text could not be retrieved, so the row was "
                 "never adjudicated and never entered any ledger. This is an unresolved question "
                 "recorded as unresolved, which is the correct handling and is also why it must be "
                 "in the archive rather than in a to-do list.",
        failure_mode="provenance_failed",
        keystone="The PDF. Nothing else.",
        keystone_predicate="ledger_rows_at_least(paired_experience_ledger.csv, 16)",
        revival_test="Obtain the full text, adjudicate the row, and re-run scripts/124.",
        status="DEAD",
    ),
    dict(
        hypothesis_id="H-tanimoto-selfread",
        claim="Tanimoto-to-known-actives predicts within-target affinity strongly (SLC6A3 +0.90, "
              "SLC6A2 +0.91, DRD1 +0.85).",
        domain="cluster_a/ranker",
        verdict="RETRACTED", died_on="2026-08-24",
        killed_by="commit fd28421 area / tanimoto_ranker.py exclude_self fix",
        evidence="The maximum was taken over an actives set CONTAINING THE QUERY, so the feature "
                 "read each row's own ChEMBL record whenever its affinity cleared the threshold. "
                 "On the Gap-4 evaluation set 143 of 289 rows scored exactly 1.000 and membership "
                 "predicted that 1.000 with no errors in either direction. No train/test split can "
                 "close it: the leak is inside the feature.",
        failure_mode="confounded",
        keystone="The quoted correlations are uninterpretable, not merely optimistic. Current "
                 "numbers live in tanimoto_baseline_v1.md and are deliberately not restated here, "
                 "because a figure typed into prose is how the previous ones outlived their code.",
        keystone_predicate="tanimoto_excludes_self()",
        revival_test="Nothing to revive; the corrected feature is live. The entry exists so the "
                     "retracted numbers stay findable by anyone who meets them in an old document.",
        status="REVIVED",
    ),
    dict(
        hypothesis_id="H-alphagenome-fit",
        claim="The AlphaGenome family (AlphaGenome, AlphaMissense, AlphaFold3) addresses one or "
              "more of this project's measured gaps.",
        domain="external systems",
        verdict="KILLED", died_on="2026-09-15",
        killed_by="reports/pipeline/external_systems_review_2026-09.md",
        evidence="Addresses NONE of the four measured gaps. G1, G3 and G4 are human-phenotype "
                 "MEASUREMENT problems; no variant-effect predictor at any accuracy produces a row "
                 "of durable-cognition ground truth. G2 is the one gap where a structural model is "
                 "the right category, and AlphaGenome is not a structural model.",
        failure_mode="wrong_endpoint",
        keystone="One narrow legitimate use survives: Cluster D's variant-to-gene mapping, which "
                 "the project has already measured as at chance for predicting trial success. The "
                 "payoff is bounded by how much of Cluster D's failure is attributable to that "
                 "step, and nothing measures that yet.",
        keystone_predicate="cluster_d_v2g_attribution_measured()",
        revival_test="Only if Cluster D's variant-to-gene step is independently shown to be the "
                     "binding constraint WITHIN Cluster D.",
        status="DEAD",
    ),
    # --- 2026-09-20: compounds supplied as "overlooked" candidates, fact-checked before entry.
    # Working in reports/pipeline/trash_can_candidates_2026-09.md: 40 of 73 claims did not survive.
    dict(
        hypothesis_id="H-dihexa-healthy",
        claim="Dihexa enhances cognition in healthy adults via HGF/c-Met-driven synaptogenesis, at "
              "far lower concentrations than BDNF.",
        domain="trash-can candidate",
        verdict="REFUTED", died_on="2026-09-20",
        killed_by="reports/pipeline/trash_can_candidates_2026-09.md",
        evidence="The mechanism does not operate on an undamaged substrate. ADDF Cognitive "
                 "Vitality, verbatim: Dihexa did not improve cognitive functions in rats with "
                 "normal cognition, and the HGF/c-Met system is not engaged during normal learning "
                 "under healthy conditions. Separately the founding mechanism paper (Benoist 2014, "
                 "PMID 25187433) was RETRACTED 2025-04-29 after a Washington State University "
                 "investigation found image manipulation, and four papers from that lab carry "
                 "Expressions of Concern (JPET 2021;378(3):311-314). The seven-orders-of-magnitude "
                 "potency figure comes from a 2012 press release; the two were never directly "
                 "compared. Zero human data exist. The mechanism WAS tested at scale via the "
                 "successor fosgonimeton: 5 registered trials, about 1,090 participants.",
        failure_mode="mechanism_contradicted",
        keystone="", keystone_predicate="", revival_test="",
        status="PERMANENTLY_CLOSED",
    ),
    dict(
        hypothesis_id="H-4mu-durability",
        claim="4-methylumbelliferone produces DURABLE cognitive gain by thinning perineuronal nets.",
        domain="trash-can candidate",
        verdict="REFUTED", died_on="2026-09-20",
        killed_by="reports/pipeline/trash_can_candidates_2026-09.md; source Dubisova 2022",
        evidence="Dubisova 2022 (PMID 35066096) is the only study that measured the washout, and it "
                 "dissociates the biomarker from the behaviour: the memory-enhancing effect did not "
                 "persist 1 month after treatment ended and the SOR score returned to pre-treatment "
                 "values, WHILE THE PNN REDUCTION DID PERSIST. That is this project's binding "
                 "constraint failing on the mechanism's own home assay. Regimen was 5 percent w/w "
                 "chow, about 6.7 mg/g/day, for six months in mice. The mechanism is also "
                 "misdescribed as PNN degradation: 4-MU inhibits hyaluronan SYNTHESIS. No study "
                 "measures 4-MU in brain or CSF, and no registered hymecromone trial carries a "
                 "cognitive or CNS endpoint.",
        failure_mode="unknown_precision",
        keystone="The washout result is a point comparison with no interval reported, so it cannot "
                 "formally exclude a smaller durable effect. Recorded as imprecise rather than as a "
                 "refutation, which is the conservative reading. The dissociation it shows is the "
                 "load-bearing part and does not depend on the interval.",
        keystone_predicate="durable_healthy_rows(1)",
        revival_test="Only if a verified durable healthy-adult row appears anywhere, at which point "
                     "re-examine whether an ECM mechanism is implicated.",
        status="DEAD",
    ),
    dict(
        hypothesis_id="H-tak071-ampa",
        claim="TAK-071 is a low-impact AMPA receptor positive allosteric modulator, and evidence "
              "that the AMPA-PAM class avoids first-generation ampakine liabilities.",
        domain="trash-can candidate",
        verdict="REFUTED", died_on="2026-09-20",
        killed_by="reports/pipeline/trash_can_candidates_2026-09.md",
        evidence="TAK-071 is a muscarinic M1 receptor positive allosteric modulator with no AMPA "
                 "mechanism at all (PMID 34240455 first-in-human, 177 healthy volunteers, "
                 "safety/PK/qEEG only, no cognitive endpoints; PMID 39761063 Phase 2 in Parkinson "
                 "disease). The low-impact-versus-first-generation framing is additionally a "
                 "category error: CX516 is both the first ampakine into human trials (PMID 9270067) "
                 "and the prototypical low-impact compound.",
        failure_mode="mechanism_contradicted",
        keystone="", keystone_predicate="", revival_test="",
        status="PERMANENTLY_CLOSED",
    ),
    dict(
        hypothesis_id="H-nsi189-healthy",
        claim="NSI-189 improves objective cognition, and the improvement persists after washout.",
        domain="trash-can candidate",
        verdict="REFUTED", died_on="2026-09-20",
        killed_by="reports/pipeline/trash_can_candidates_2026-09.md",
        evidence="The Phase 2 signal is SELF-REPORT only and dose-discordant: at 40 mg, SDQ -8.2 "
                 "(p = 0.04) and MGH-CPFQ -1.9 (p = 0.035); at 80 mg all null. The objective "
                 "Cogstate battery showed NO significant difference at either dose (PMID 30626911). "
                 "The trial had no washout or follow-up period, so the persistence claim asserts "
                 "something never measured. The molecule is now ALTO-100 / amdiglurax, and Alto's "
                 "Phase 2b (NCT05712187, n = 301) enriched for a memory-linked biomarker and still "
                 "missed its primary endpoint in October 2024. The one healthy-volunteer study "
                 "(NCT01310881) administered no cognitive battery.",
        failure_mode="wrong_population",
        keystone="Every cognitive datum is in depressed patients and is subjective. No healthy adult "
                 "has been given NSI-189 and tested. The claim is untested in the population it is "
                 "about, rather than refuted in it.",
        keystone_predicate="ci_recorded(nsi_189)",
        revival_test="Re-adjudicate if a healthy-adult NSI-189 cognitive row with an interval is "
                     "ever curated into the healthy-adult ledger.",
        status="DEAD",
    ),
    dict(
        hypothesis_id="H-klotho-glun2b",
        claim="Peripherally administered klotho enhances cognition by driving GluN2B-containing "
              "NMDA receptors into the postsynaptic density.",
        domain="trash-can candidate",
        verdict="REFUTED", died_on="2026-09-20",
        killed_by="reports/pipeline/trash_can_candidates_2026-09.md",
        evidence="The GluN2B result belongs to LIFELONG TRANSGENIC OVEREXPRESSION (Dubal 2014, PMID "
                 "24813892) and is explicitly contradicted in the acute peripheral-protein paradigm "
                 "that the macaque protocol actually models: Leon 2017 (PMID 28793260) reports that "
                 "synaptic GluN2B did NOT differ between vehicle and drug despite cognitive "
                 "enhancement. The cognitive effect is real; this explanation of it is not.",
        failure_mode="mechanism_contradicted",
        keystone="", keystone_predicate="", revival_test="",
        status="PERMANENTLY_CLOSED",
    ),
    dict(
        hypothesis_id="H-roflumilast-durability",
        claim="Low-dose roflumilast (100 ug) produces durable cognitive gain in healthy adults.",
        domain="trash-can candidate",
        verdict="ABSTAINED", died_on="2026-09-20",
        killed_by="reports/pipeline/trash_can_candidates_2026-09.md",
        evidence="NOT TESTED, rather than failed. Roflumilast 100 ug is the only compound in the "
                 "supplied list with genuine healthy-adult cognitive evidence: PMID 29241652 "
                 "(healthy young, immediate recall up 2-3 words on a 30-word verbal learning task, "
                 "enhanced P600) and PMID 30776650 (healthy 60-80, n = 20, delayed recall d = 0.69). "
                 "EVERY such study is a single acute dose about 1 h before testing, so durability "
                 "has never been measured. Caveats on the acute claim: the two studies disagree "
                 "about which endpoint moves (immediate in the young, delayed in the old, spatial "
                 "in neither, which is G3 assay-dependence inside one drug and one group); both "
                 "come from the same group with industry co-authors; NCT02051335 (n = 27) found "
                 "monotherapy did NOT improve verbal recall; and ROSTMEMA (PMID 42032844, n = 100) "
                 "missed its pre-specified primary analysis.",
        failure_mode="wrong_endpoint",
        keystone="Nobody has dosed a healthy adult with roflumilast and retested off drug. This is "
                 "the most testable durability question the supplied list produced, and it is a "
                 "curation and study-design task rather than a modelling one.",
        keystone_predicate="durable_healthy_rows(1)",
        revival_test="If any verified durable healthy-adult row appears, roflumilast is the cheapest "
                     "replication target in the list: approved, oral, sub-emetic at 100 ug, with two "
                     "acute positives to anchor a washout design.",
        status="DEAD",
    ),
    # --- 2026-09-20: recommendations from the frontier-systems survey, checked and killed.
    # These exist so the next survey cannot re-propose them blind. The survey that produced them
    # re-proposed Rokem & Silver 2013 as a novel finding precisely because this project's reasoning
    # about it was not machine-findable. Full working: reports/pipeline/frontier_systems_review_2026-09.md
    dict(
        hypothesis_id="H-chembl-allosteric-g2",
        claim="Expanding the ChEMBL directional-allosteric label layer makes G2 measurable, giving a "
              "roughly 500x bigger labelled test set than the handful of AMPA rows the AUROC 0.26 "
              "verdict rests on.",
        domain="external systems / G2",
        verdict="REFUTED", died_on="2026-09-20",
        killed_by="direct query of the local ChEMBL 36 mirror; "
                  "reports/pipeline/frontier_systems_review_2026-09.md",
        evidence="The aggregate claim holds: 6,863 PAM + 1,986 NAM = 8,849 directional allosteric "
                 "activity rows. The breakdown refutes the use. Top targets are muscarinic M4 (691), "
                 "adenosine A3 (613), GLP-1 (297), FFA4 (294), alpha-7 nAChR (278), overwhelmingly "
                 "GPCRs. ALL ionotropic glutamate receptors together are 396 rows of 8,849. At AMPA, "
                 "the target whose AUROC 0.26 the G2 verdict rests on: 75 PAM rows and ZERO NAM "
                 "rows. A PAM-versus-NAM classifier cannot be tested at a target with no negatives.",
        failure_mode="instrument_blind",
        keystone="The corpus does not index the thing. ChEMBL records no AMPA negative allosteric "
                 "modulators at all, so the label layer can only build a cross-target, "
                 "GPCR-dominated benchmark. That is a different question, and transfer from GPCR "
                 "allostery to ionotropic channel gating is the assumption this project should not "
                 "make for free.",
        keystone_predicate="chembl_allosteric_negatives_at(AMPA, 40)",
        revival_test="If ChEMBL ever indexes 40 or more AMPA NAM rows, re-run the G2 permutation "
                     "test on an AMPA-specific PAM-vs-NAM set. Until then the honest next step is "
                     "AMPA-specific labelled data, which this proposal does not supply.",
        status="DEAD",
    ),
    dict(
        hypothesis_id="H-knowledge-graph-repurposing",
        claim="A knowledge-graph repurposing platform (TxGNN, PrimeKG, Open Targets, Every Cure "
              "MATRIX) can extend or beat this project's healthy-adult cognition ledger.",
        domain="external systems / G1",
        verdict="REFUTED", died_on="2026-09-20",
        killed_by="reports/pipeline/frontier_systems_review_2026-09.md (live API checks)",
        evidence="Structural, not a performance question. Every platform is keyed to a disease "
                 "vocabulary, MONDO or EFO. An OLS4 query of MONDO for 'cognitive enhancement' "
                 "returns numFound = 0; every cognition-adjacent MONDO term is a DEFICIT. Open "
                 "Targets EFO_0008354 'cognitive function measurement' has 1,945 associated targets "
                 "and drugAndClinicalCandidates = 0. The parent term 'cognition' has 4,040 targets "
                 "and 11 drug candidates in total, of which only caffeine and creatine are "
                 "plausible enhancers and lorazepam is an impairer. Roughly two usable positive "
                 "edges, against six labelled enhancers the ledger already has.",
        failure_mode="instrument_blind",
        keystone="The outcome is not nameable in the vocabulary these systems are built on. Healthy "
                 "cognitive enhancement is not a disease, and a disease ontology has no node for "
                 "it, so the platforms cannot reproduce the existing ledger let alone extend it.",
        keystone_predicate="ontology_names_the_outcome()",
        revival_test="Re-check only if MONDO or EFO gains a non-deficit cognitive-enhancement term "
                     "with drug edges attached.",
        status="DEAD",
    ),
    dict(
        hypothesis_id="H-cohort-prescription-emulation",
        claim="A prescription-exposure analysis or target-trial emulation in a population cohort "
              "(UK Biobank and similar) can supply durable healthy-adult cognitive rows.",
        domain="external systems / G1",
        verdict="REFUTED", died_on="2026-09-20",
        killed_by="reports/pipeline/frontier_systems_review_2026-09.md",
        evidence="UK Biobank has baseline cognitive data on 480,416 participants and CANNOT support "
                 "the analysis: its prescription records stop in 2016-2017 while the later "
                 "cognitive wave is 2021, so no exposure-then-cognition contrast can be formed. No "
                 "cohort in the lane measures cognition after a documented discontinuation. "
                 "Reliability compounds it: UK Biobank cognitive tests have 4-week retest r = "
                 "0.41-0.61 and longitudinal ICC 0.16-0.65 (PMID 32310977, PMID 27110937), "
                 "overlapping the ICC 0.16-0.53 for which this project already demoted cTBS.",
        failure_mode="wrong_endpoint",
        keystone="The cohorts can only manufacture NULLS, and nulls are the row type worth least: "
                 "measured here, +100 nulls buys +0.11% power per row against +4.60% for a "
                 "positive. A route that supplies only the cheap row type does not relieve a "
                 "constraint set by the expensive one.",
        keystone_predicate="durable_healthy_rows(1)",
        revival_test="Re-assess only if a cohort links prescription or exposure records to a "
                     "post-discontinuation cognitive wave.",
        status="DEAD",
    ),
    dict(
        hypothesis_id="H-cloudlab-experimental-arm",
        claim="This project can acquire an experimental arm by buying assays from a self-driving or "
              "cloud lab, without building a laboratory.",
        domain="external systems",
        verdict="REFUTED", died_on="2026-09-20",
        killed_by="reports/pipeline/frontier_systems_review_2026-09.md",
        evidence="No cloud lab sells the assay. Emerald Cloud Lab advertises over 200 instrument "
                 "models and names HPLC, MS, NMR, PCR, ELISA, SPR and flow cytometry, but lists no "
                 "multielectrode arrays, no patch clamp, no high-content spine imaging and no "
                 "primary or iPSC neuron culture. Arctoris runs biochemistry, biophysics, cell "
                 "biology and structural biology with no electrophysiology. Where prices are "
                 "published they are high and stale: ECL ~$25k/month on a one-year minimum, ~$300k "
                 "entry; Strateos ~$130k entry for a single automated method (DOI "
                 "10.33552/OJRAT.2022.01.000511). Neither publishes a current rate card, so a real "
                 "quote is NOT FOUND.",
        failure_mode="wrong_endpoint",
        keystone="Even if the assay were purchasable it would add in-vitro rows, and G1 is a "
                 "HUMAN-cognition constraint. At n = 21 with a 0.73 critical AUROC, no assay "
                 "purchase changes what the ledger can credit. The one free item in the lane that "
                 "touches a measured gap is NIMH PDSP screening, which is recorded separately as "
                 "open rather than killed.",
        keystone_predicate="durable_healthy_rows(1)",
        revival_test="Only relevant once the ledger can adjudicate predictions at all.",
        status="DEAD",
    ),
]


def derive_compound_entries() -> list[dict]:
    """Classify every compound the ledger labels non-enhancing, straight from its interval.

    Derived rather than typed so the archive cannot drift from the ledger. The classification is
    `classify_null`, which asks the only question that matters: does the interval EXCLUDE an effect
    of the size this project is looking for, or merely fail to find one?
    """
    df = pd.read_csv(LEDGER)
    cand = df.get("candidate_enhancer", pd.Series([1] * len(df)))
    # Mirror the primary-set rule in scripts/121 and 131: a row whose effect is not independent of
    # another already here (a caffeine combination alongside caffeine) is never pooled. No such row
    # is currently a NULL, but the filter belongs here so a future one cannot leak in unnoticed.
    dep = df.get("depends_on", pd.Series([""] * len(df))).fillna("").astype(str).str.strip()
    sel = df[(df["evidence_tier"] == "clean_MA") & (cand == 1) & (dep == "")
             & (df["enhances_healthy_young"] == 0)]
    rows = []
    for _, r in sel.iterrows():
        c = str(r["compound"])
        mode = classify_null(r.get("ci_lo"), r.get("ci_hi"))
        closed = not FAILURE_MODES[mode].revivable
        ci = (f"CI [{r['ci_lo']:+.3f}, {r['ci_hi']:+.3f}]"
              if pd.notna(r.get("ci_lo")) and pd.notna(r.get("ci_hi")) else "NO INTERVAL RECORDED")
        if mode == "unknown_precision":
            keystone = ("No interval was recorded, so this compound's null cannot be distinguished "
                        "from silence. It is not evidence of absence; it is absence of evidence.")
            pred, test = f"ci_recorded({c})", (
                "Re-extract the meta-analytic interval from the source, then re-classify. If the "
                "upper bound falls below the target this becomes a measured_null and closes.")
        elif mode == "underpowered":
            keystone = (f"The interval ADMITS an effect of practical size: its upper bound reaches "
                        f"{r['ci_hi']:+.3f} against a target of {TARGET_EFFECT_G}. A label of 0 "
                        "here asserts a refutation the data do not contain.")
            pred, test = f"interval_excludes_target({c})", (
                "Await or curate more precise evidence. This keystone resolves DOWNWARD: when it "
                "fires, the hypothesis is properly closed rather than revived.")
        else:
            keystone, pred, test = "", "", ""
        rows.append(dict(
            hypothesis_id=f"H-null-{c}",
            claim=f"{c} enhances cognition in healthy adults.",
            domain="healthy-adult ledger",
            verdict="REFUTED", died_on="2026-07-28",
            killed_by="data/raw/healthy_adult_cognition_ledger.csv (enhances_healthy_young = 0); "
                      "classification derived by scripts/129_stepping_stone_archive.py",
            evidence=f"g = {r['representative_g']:+.3f}, {ci}, n_studies = {r.get('n_studies')}, "
                     f"source {r.get('pmid_doi')}. Classified {mode} against target "
                     f"g = {TARGET_EFFECT_G}.",
            failure_mode=mode,
            keystone=keystone, keystone_predicate=pred, revival_test=test,
            status="PERMANENTLY_CLOSED" if closed else "DEAD",
        ))
    return rows


def main() -> int:
    rows = derive_compound_entries() + CURATED
    ids = [r["hypothesis_id"] for r in rows]
    assert len(ids) == len(set(ids)), "duplicate hypothesis_id"

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(REQUIRED_COLUMNS))
        w.writeheader()
        for r in sorted(rows, key=lambda x: x["hypothesis_id"]):
            w.writerow({k: r.get(k, "") for k in REQUIRED_COLUMNS})

    df = pd.read_csv(OUT)
    viol = validate_archive(df)
    errs = [v for v in viol if v.severity == "error"]
    for v in viol:
        (L.error if v.severity == "error" else L.warning)("%s", v)
    L.info("wrote %s: %d entries (%d derived, %d curated)",
           OUT.relative_to(ROOT), len(df), len(rows) - len(CURATED), len(CURATED))
    L.info("failure modes: %s", dict(df["failure_mode"].value_counts()))
    revivable = df[df["failure_mode"].map(lambda m: FAILURE_MODES[m].revivable)]
    L.info("revivable: %d | permanently closed: %d", len(revivable), len(df) - len(revivable))
    if errs:
        L.error("archive FAILED its own contract (%d error-severity violations)", len(errs))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
