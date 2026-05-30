# Prospective Repurposing Shortlist (Gap 7 — capstone)

The deliverable the pipeline was built for: **approved drugs as mechanism-justified repurposing hypotheses** for cognitive-impairment diseases. Each candidate is ranked by the one signal the project proved is prognostic — its mechanism class's pivotal-trial track record in that disease (Gap 3 / Gap 6) — times how well it engages a disease-relevant target. Restricted to **SUCCESS-track-record mechanism classes**.

> These are hypotheses worth evaluation, **not** predicted cures. Effect sizes remain bounded by the Roberts 2020 ceiling; engagement flagged *uncertain* is at an allosteric/transporter target where MAMMAL's sequence-only binding is unreliable (Gap 4). Novelty is flagged against a curated indication map — a **NOVEL** row is an approved drug whose mechanism class succeeds in this disease but which is not currently used for it.

---

## AD

10 drugs in SUCCESS-track-record classes engage a AD-relevant target; **4 are NOVEL repurposing hypotheses** (approved elsewhere, not used for AD).

| Rank | Drug | Currently used for | Mechanism class | Target | Class g | Engagement | NOVEL? |
|---|---|---|---|---|---|---|---|
| 1 | donepezil | AD | AChE_inhibitor | P22303 | +0.37 | 96% | standard |
| 2 | huperzine A | — | AChE_inhibitor | P22303 | +0.37 | 90% | 🔵 **NOVEL** |
| 3 | rivastigmine | AD, PD | AChE_inhibitor | P22303 | +0.37 | 62% | standard |
| 4 | rapastinel | unknown | NMDA_modulator | Q12879 | +0.29 | 77% | standard |
| 5 | galantamine | AD | AChE_inhibitor | P22303 | +0.37 | 12% | standard |
| 6 | lanicemine | unknown | NMDA_modulator | Q12879 | +0.29 | 42% | standard |
| 7 | blarcamesine | — | sigma1 | Q99720 | +0.24 | 58% | 🔵 **NOVEL** |
| 8 | fluvoxamine | MDD, OCD | sigma1 | Q99720 | +0.24 | 45% | 🔵 **NOVEL** |
| 9 | pridopidine | — | sigma1 | Q99720 | +0.24 | 41% | 🔵 **NOVEL** |
| 10 | memantine | AD | NMDA_modulator | Q13224 | +0.29 | 10% | standard |

### Top NOVEL AD repurposing hypotheses

- **huperzine A** (AChE_inhibitor, currently —): AChE_inhibitor is a SUCCESS-track-record class in AD (class g≈+0.37); huperzine A engages P22303 at 90% percentile. Currently used for: — → NOVEL repurposing.
  - liability flags: CHRM3, TACR1, OPRM1
- **blarcamesine** (sigma1, currently —): sigma1 is a SUCCESS-track-record class in AD (class g≈+0.24); blarcamesine engages Q99720 at 58% percentile. Currently used for: — → NOVEL repurposing.
  - liability flags: CHRM3, TACR1, OPRM1
- **fluvoxamine** (sigma1, currently MDD, OCD): sigma1 is a SUCCESS-track-record class in AD (class g≈+0.24); fluvoxamine engages Q99720 at 45% percentile. Currently used for: MDD, OCD → NOVEL repurposing.
  - liability flags: CHRM3, TACR1, OPRM1
- **pridopidine** (sigma1, currently —): sigma1 is a SUCCESS-track-record class in AD (class g≈+0.24); pridopidine engages Q99720 at 41% percentile. Currently used for: — → NOVEL repurposing.
  - liability flags: CHRM3, TACR1, OPRM1

---

## CIAS

7 drugs in SUCCESS-track-record classes engage a CIAS-relevant target; **5 are NOVEL repurposing hypotheses** (approved elsewhere, not used for CIAS).

| Rank | Drug | Currently used for | Mechanism class | Target | Class g | Engagement | NOVEL? |
|---|---|---|---|---|---|---|---|
| 1 | buspirone | anxiety | 5HT1A_partial_agonist | P08908 | +0.40 | 73% ⚠ | 🔵 **NOVEL** |
| 2 | aripiprazole | MDD, schizophrenia | D1_agonist | P21728 | +0.40 | 68% | standard |
| 3 | tandospirone | anxiety | 5HT1A_partial_agonist | P08908 | +0.40 | 50% ⚠ | 🔵 **NOVEL** |
| 4 | cevimeline | Sjogren | M1_M4_agonist | P11229 | +0.38 | 50% | 🔵 **NOVEL** |
| 5 | pilocarpine | glaucoma | M1_M4_agonist | P11229 | +0.38 | 50% | 🔵 **NOVEL** |
| 6 | xanomeline | CIAS | M1_M4_agonist | P11229 | +0.38 | 50% | standard |
| 7 | emraclidine | — | M1_M4_agonist | P08173 | +0.38 | 50% | 🔵 **NOVEL** |

### Top NOVEL CIAS repurposing hypotheses

- **buspirone** (5HT1A_partial_agonist, currently anxiety): 5HT1A_partial_agonist is a SUCCESS-track-record class in CIAS (class g≈+0.40); buspirone engages P08908 at 73% percentile (engagement uncertain — allosteric/transporter). Currently used for: anxiety → NOVEL repurposing.
  - liability flags: CHRM3, TACR1, OPRM1
- **tandospirone** (5HT1A_partial_agonist, currently anxiety): 5HT1A_partial_agonist is a SUCCESS-track-record class in CIAS (class g≈+0.40); tandospirone engages P08908 at 50% percentile (engagement uncertain — allosteric/transporter). Currently used for: anxiety → NOVEL repurposing.
- **cevimeline** (M1_M4_agonist, currently Sjogren): M1_M4_agonist is a SUCCESS-track-record class in CIAS (class g≈+0.38); cevimeline engages P11229 at 50% percentile. Currently used for: Sjogren → NOVEL repurposing.
- **pilocarpine** (M1_M4_agonist, currently glaucoma): M1_M4_agonist is a SUCCESS-track-record class in CIAS (class g≈+0.38); pilocarpine engages P11229 at 50% percentile. Currently used for: glaucoma → NOVEL repurposing.

---

## FXS

3 drugs in SUCCESS-track-record classes engage a FXS-relevant target; **2 are NOVEL repurposing hypotheses** (approved elsewhere, not used for FXS).

| Rank | Drug | Currently used for | Mechanism class | Target | Class g | Engagement | NOVEL? |
|---|---|---|---|---|---|---|---|
| 1 | bpn14770 | FXS | PDE4_inhibitor | Q08499 | +0.71 | 80% ⚠ | standard |
| 2 | rolipram | — | PDE4_inhibitor | Q08499 | +0.71 | 76% ⚠ | 🔵 **NOVEL** |
| 3 | roflumilast | COPD | PDE4_inhibitor | Q08499 | +0.71 | 50% ⚠ | 🔵 **NOVEL** |

### Top NOVEL FXS repurposing hypotheses

- **rolipram** (PDE4_inhibitor, currently —): PDE4_inhibitor is a SUCCESS-track-record class in FXS (class g≈+0.71); rolipram engages Q08499 at 76% percentile (engagement uncertain — allosteric/transporter). Currently used for: — → NOVEL repurposing.
  - liability flags: CHRM3, TACR1, OPRM1
- **roflumilast** (PDE4_inhibitor, currently COPD): PDE4_inhibitor is a SUCCESS-track-record class in FXS (class g≈+0.71); roflumilast engages Q08499 at 50% percentile (engagement uncertain — allosteric/transporter). Currently used for: COPD → NOVEL repurposing.

---

## Evidence dossiers — top novel repurposing picks

### huperzine A — for AD

| | |
|---|---|
| **Mechanism** | AChE_inhibitor at P22303 |
| **Predicted cognition effect** | Hedges' g = **+0.37** (90% CrI +0.27 to +0.47) |
| **Evidence quality (GRADE)** | 🟡 MODERATE |
| **Mechanism-class track record** | SUCCESS (3 class drugs, 40 pooled RCTs) |
| **Class exemplars** | Donepezil, Galantamine, Rivastigmine |
| **Effect basis** | mechanism-class prediction |

**Why this GRADE rating:** RCT meta-analytic base (k=40) → start HIGH; indirectness: effect predicted from mechanism-class siblings, not this compound's own pivotal trial (−1).

**Predicted off-target liability flags** (model-based, unvalidated):
- `CHRM3` (pKd≈7.1) — M3: anticholinergic (dry mouth, GI, urinary)
- `TACR1` (pKd≈6.9) — NK1: emetic / CNS off-target
- `OPRM1` (pKd≈6.8) — µ-opioid: dependence / respiratory depression
- `GABRA1` (pKd≈6.7) — GABA-A α1: sedation / dependence

**Provenance:** V6.A multi-head DTI (target engagement) · V6.B Cluster D θ̄ (cognition relevance) · V7 / disease-conditioned class prior (effect size) · Gap-3 mechanism-class track record (prognosis) · 44-target off-target liability panel.

**Caveats & failure modes:**
- Predicted cognition effect is bounded by the Roberts 2020 ceiling (healthy-adult SMD ≈ 0.2-0.5); disease effects can be larger.
- This is a mechanism-class prediction, NOT a per-compound clinical result — the effect size is the class ceiling scaled by predicted target engagement.
- Off-target liability flags are MODEL-PREDICTED (MAMMAL DTI), unvalidated — confirm against the drug's known safety profile.

---

### buspirone — for CIAS

| | |
|---|---|
| **Mechanism** | 5HT1A_partial_agonist at P08908 |
| **Predicted cognition effect** | Hedges' g = **+0.40** (90% CrI +0.21 to +0.59) |
| **Evidence quality (GRADE)** | 🟠 LOW |
| **Mechanism-class track record** | SUCCESS (1 class drugs, 2 pooled RCTs) |
| **Class exemplars** | Tandospirone |
| **Effect basis** | mechanism-class prediction |

**Why this GRADE rating:** limited trial base (k=2); indirectness: effect predicted from mechanism-class siblings, not this compound's own pivotal trial (−1); risk of bias: effect relies on MAMMAL sequence-only target-engagement, which is allosteric-blind here (−1).

**Predicted off-target liability flags** (model-based, unvalidated):
- `CHRM3` (pKd≈7.2) — M3: anticholinergic (dry mouth, GI, urinary)
- `TACR1` (pKd≈6.9) — NK1: emetic / CNS off-target
- `OPRM1` (pKd≈6.9) — µ-opioid: dependence / respiratory depression
- `GABRA2` (pKd≈6.7) — GABRA2: off-target engagement

**Provenance:** V6.A multi-head DTI (target engagement) · V6.B Cluster D θ̄ (cognition relevance) · V7 / disease-conditioned class prior (effect size) · Gap-3 mechanism-class track record (prognosis) · 44-target off-target liability panel · Gap-4 allosteric LTR (binding-reliability flag).

**Caveats & failure modes:**
- Predicted cognition effect is bounded by the Roberts 2020 ceiling (healthy-adult SMD ≈ 0.2-0.5); disease effects can be larger.
- This is a mechanism-class prediction, NOT a per-compound clinical result — the effect size is the class ceiling scaled by predicted target engagement.
- MAMMAL's sequence-only binding score is structurally blind at this (allosteric/transporter) target; engagement is uncertain (see Gap-4 allosteric learn-to-rank).
- Off-target liability flags are MODEL-PREDICTED (MAMMAL DTI), unvalidated — confirm against the drug's known safety profile.

---

### rolipram — for FXS

| | |
|---|---|
| **Mechanism** | PDE4_inhibitor at Q08499 |
| **Predicted cognition effect** | Hedges' g = **+0.71** (90% CrI +0.52 to +0.90) |
| **Evidence quality (GRADE)** | 🔴 VERY LOW |
| **Mechanism-class track record** | SUCCESS (1 class drugs, 1 pooled RCTs) |
| **Class exemplars** | Zatolmilast |
| **Effect basis** | mechanism-class prediction |

**Why this GRADE rating:** limited trial base (k=1); imprecision: wide CI / few RCTs (−1); indirectness: effect predicted from mechanism-class siblings, not this compound's own pivotal trial (−1); risk of bias: effect relies on MAMMAL sequence-only target-engagement, which is allosteric-blind here (−1).

**Predicted off-target liability flags** (model-based, unvalidated):
- `CHRM3` (pKd≈7.2) — M3: anticholinergic (dry mouth, GI, urinary)
- `TACR1` (pKd≈7.0) — NK1: emetic / CNS off-target
- `OPRM1` (pKd≈6.9) — µ-opioid: dependence / respiratory depression
- `ADRA1D` (pKd≈6.8) — α1D: orthostatic hypotension

**Provenance:** V6.A multi-head DTI (target engagement) · V6.B Cluster D θ̄ (cognition relevance) · V7 / disease-conditioned class prior (effect size) · Gap-3 mechanism-class track record (prognosis) · 44-target off-target liability panel · Gap-4 allosteric LTR (binding-reliability flag).

**Caveats & failure modes:**
- Predicted cognition effect is bounded by the Roberts 2020 ceiling (healthy-adult SMD ≈ 0.2-0.5); disease effects can be larger.
- This is a mechanism-class prediction, NOT a per-compound clinical result — the effect size is the class ceiling scaled by predicted target engagement.
- MAMMAL's sequence-only binding score is structurally blind at this (allosteric/transporter) target; engagement is uncertain (see Gap-4 allosteric learn-to-rank).
- Off-target liability flags are MODEL-PREDICTED (MAMMAL DTI), unvalidated — confirm against the drug's known safety profile.

---

Generated by `scripts/82_repurposing_shortlist.py` via `reporting/repurposing_shortlist.py` (integrates Gap 2 disease priors + Gap 4 engagement-reliability + Gap 5 dossiers).