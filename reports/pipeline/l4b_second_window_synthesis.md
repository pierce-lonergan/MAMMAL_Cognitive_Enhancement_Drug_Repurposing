# L4b second-window research synthesis - non-serotonergic durable plasticity

Synthesis of two parallel Opus-4.8 research lanes (full reports:
`research_nmda_durability_window.md`, `research_neurosteroid_convergent_window.md`) plus two
in-repo verifications. The question: can PERSEUS add a SECOND off-axis, structure-computable
durability window for the non-serotonergic rapid-acting antidepressants (the L4 serotonergic window
is blind to them, which is the bulk of the 7/16 ledger recall gap)?

## Bottom line (decision)

**No structural second window is justified.** For these classes the durability discriminator lives
off BOTH the binding-affinity axis (what MAMMAL's DTI head sees) AND the cheap-ADMET / 2D-descriptor
axis (what the serotonergic L4 permeability trick exploits). The honest architecture is:
mechanism-class scaffold ROUTERS that ABSTAIN-with-reason on durability (not predict it), a small
curated pharmacodynamic lookup for the one class where a real discriminator exists but is not
structural (NMDA trapping kinetics), and the existing L1 free-exposure gate + evidence-design
governor as the genuine cross-class objects. This is a publishable negative, parallel to the
TrkB-TMD and AMPA-PAM off-axis results.

## Lane A - NMDA / dissociative class

- MECHANISM (established, cited in the lane-A report): a brief uncompetitive NMDA-channel block
  disinhibits pyramidal neurons -> prefrontal glutamate surge -> AMPA/GluA1 + halted eEF2K ->
  de-suppressed BDNF -> TrkB -> mTORC1 -> dendritic-spine synthesis. The antagonist is the transient
  TRIGGER; the BDNF-TrkB-mTORC1-spine program is the self-sustaining substrate that outlasts the
  ~2-3 h half-life (Li/Duman 2010; Autry/Monteggia 2011; Moda-Sava 2019).
- HONESTY (per compound, from the lane-A report): ketamine/esketamine single-dose benefit is real
  but ~1 week then relapse (maintenance re-dosing needed); HNK is genuinely contested
  (NMDAR-independence disputed, replication inconsistent, higher plasma HNK predicted WORSE human
  response); N2O durability beyond ~1 week is unestablished; DXM/AXS-05 is chronically dosed and
  multi-target. memantine is the established NEGATIVE (no MDD efficacy; no BDNF/eEF2 response).
- THE DISCRIMINATOR is channel TRAPPING / resting-state, use-dependent block (Gideons 2014): in
  physiological Mg2+ only ketamine blocks resting/spontaneous NMDARs and triggers BDNF; memantine is
  a partial trapper that spares the resting pool. Subtype (GluN2B) preference is SHARED and does not
  discriminate.
- VERIFIED IN-REPO (the linchpin): RDKit descriptors cannot separate the positive from the negative.
  ketamine clogP 2.90 / TPSA 29.1 vs memantine clogP 2.69 / TPSA 26.0 (delta 0.21 / 3.1); both are
  L4 window-NEGATIVE (no serotonergic scaffold). So the L4 permeability strategy provably does not
  transfer - the discriminator is a pharmacodynamic channel-state property a structure model cannot
  derive. (amantadine 1.91/26.0, esketamine 2.90/29.1, dextromethorphan 3.38/12.5 also computed.)
- PROPOSED (lane A): an L4b NMDA router = cheap scaffold flag (arylcyclohexylamine / aminoadamantane
  / morphinan + L1 CNS gate, necessary-not-sufficient) GATED BY a small curated PD table
  (`trapping_class`, `blocks_resting_NMDAR`, `use_dependence` with source PMIDs); ABSTAIN by default
  for any blocker whose trapping/resting-block is unmeasured (HNK, N2O, DXM all abstain).
  Pre-registered falsifiers: memantine + amantadine window-NEGATIVE; ketamine/esketamine
  window-POSITIVE; a structure-only ablation must FAIL to separate memantine from ketamine (already
  confirmed above).

## Lane B - GABA-A neurosteroids + the convergent pathway

- NEUROSTEROID DURABILITY (honest, from lane-B report): brexanolone sustained within-group to Day 30
  (controlled drug-vs-placebo significance at Day 30 NOT verified by the agent). Zuranolone is the
  crux: a fixed 14-day course; WATERFALL met Day-15 primary with an ~6-week within-responder tail,
  but the companion MOUNTAIN trial MISSED its Day-15 primary, and neither MDD trial showed Day-42
  placebo separation; SHORELINE shows episodic repeat-as-needed use. Lane-B verdict: rapid with a
  multi-week within-responder tail, NOT post-cessation disease-modification.
- CONVERGENT PATHWAY (real, well-cited: Zanos/Duman/Gould 2018; Duman 2019; Kavalali & Monteggia
  2012; Anacker 2018): chemically diverse rapid-acting antidepressants converge on
  glutamate -> AMPA/GluA1 -> BDNF -> TrkB -> mTORC1 -> synaptogenesis. BUT this convergence is purely
  downstream/intracellular (effectors, not drug targets), so it is INVISIBLE from ligand structure.
  Confirmed by the repo's own grouped-LOMO audit (serotonergic window = 0.00 recall on
  gaba_neurosteroid / muscarinic / nmda).
- VERDICT (lane B): no single cross-class structure-computable signature. scopolamine -> a narrow
  tropane/tertiary-amine/CNS-PASS class router is buildable (with a quaternary-tropane veto mirroring
  the permanent-charge logic) but detects MECHANISM, not durability, so it routes to the evidence
  layer only. neurosteroids -> ABSTAIN at the structural layer (the pregnane scaffold maps to a
  shallow/episodic effect, so a window flag would be a false positive).

## Ledger integrity check (verification CAUGHT a research overstatement)

Lane B recommended re-tagging two ledger ground-truth positives (zuranolone, scopolamine) as
"washout-observation, not durable". I checked this against the ledger's ACTUAL cited basis before
acting and found the recommendation OVERSTATED - no re-tag is warranted:
- **Zuranolone** (ledger row 12) is grounded in the POSTPARTUM trial (Deligiannidis 2021, JAMA
  Psychiatry, PMID 34190962), which DID show a post-cessation Day-45 separation (HAMD-17 difference
  -4.1, 95% CI -6.7 to -1.4, P=.003; Day-45 remission 53% vs 30%). Lane B critiqued the weaker MDD
  program (WATERFALL/MOUNTAIN), which is NOT the ledger's basis. The ledger row already states
  "durability beyond Day 45 was not assessed" - the scope is honestly bounded, not over-credited.
- **Scopolamine** (ledger row 10) cites a REPLICATED 12-16 day post-cessation carryover (Furey &
  Drevets 2006 + the 2010 replication), where improvement persisted into a subsequent placebo block.
  That is a genuine after-cessation persistence, accurately described.
Conclusion: both rows are correctly cited and stay as-is. The verification prevented a wrong
downgrade of two true positives - the same "check the hypothesis against the whole ground-truth set"
discipline that produced the D5/D6 wins. The legitimate residual point (these are shorter-window /
course-based positives than the single-dose-with-long-followup psychedelics) is a candidate for a
future STRUCTURED durability-strength tier field in the ledger; it is already captured narratively in
each row's `persistence_design`.

## Carried-forward uncertainties (do NOT promote to any ground-truth/tradable claim)

From the agents, flagged honestly: paywalled publishers (PNAS/Science/Nature) returned 403, so some
volume/page numbers came from indexer metadata and need spot-checking before formal citation; HNK
replication specifics are "contested" not numeric; brexanolone Day-30 CONTROLLED significance
unverified; the zuranolone 86.1% figure is a company within-responder number, not a controlled
effect size. None of these entered an engine decision or a ground-truth label.

## Next concrete build (designed, deferred for greenlight)

Not built in this lane (it is curated-data + engine work, not a structural innovation, and the
honest headline is the negative): (1) an abstain-with-reason mechanism-class router module
(NMDA arylcyclohexylamine/aminoadamantane, neurosteroid pregnane, muscarinic tropane) that emits
"mechanism recognized; durability not structure-derivable -> route to curated-PD / evidence layer"
instead of a silent miss, plus the pre-registered structural-ablation negative (memantine vs
ketamine) as a regression test; (2) the small curated NMDA trapping-kinetics table with PMIDs. Both
reuse the `psychoplastogen.py` pattern and keep abstain-by-default.
