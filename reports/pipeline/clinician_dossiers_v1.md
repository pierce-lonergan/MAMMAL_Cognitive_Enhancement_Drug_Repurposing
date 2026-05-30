# Clinician Evidence Dossiers (Gap 5)

One-page, GRADE-style evidence cards — the single artifact a clinician reads instead of the full report suite. Each distils the pipeline's real outputs into: predicted cognition effect size + credible interval, Cochrane-GRADE evidence quality with explicit up/down-grade reasons, the mechanism-class pivotal-trial track record (the Gap-3 prognostic signal), predicted off-target liability flags, the provenance trail, and explicit failure-mode caveats.

> Effect sizes are predicted cognition Hedges' g bounded by the Roberts 2020 ceiling; off-target flags are model-predicted (MAMMAL DTI) and unvalidated. This is a triage aid, not a prescribing guide.

---

### donepezil — for AD

| | |
|---|---|
| **Mechanism** | AChE_inhibitor at ACHE |
| **Predicted cognition effect** | Hedges' g = **+0.36** (90% CrI +0.27 to +0.44) |
| **Evidence quality (GRADE)** | 🟢 HIGH |
| **Mechanism-class track record** | SUCCESS (3 class drugs, 18 pooled RCTs) |
| **Class exemplars** | Donepezil, Galantamine, Rivastigmine |
| **Effect basis** | this compound’s own pivotal trial |

**Why this GRADE rating:** RCT meta-analytic base (k=18) → start HIGH.

**Predicted off-target liability flags** (model-based, unvalidated):
- `CHRM3` (pKd≈7.3) — M3: anticholinergic (dry mouth, GI, urinary)
- `TACR1` (pKd≈7.0) — NK1: emetic / CNS off-target
- `OPRM1` (pKd≈7.0) — µ-opioid: dependence / respiratory depression
- `ADRA1D` (pKd≈6.8) — α1D: orthostatic hypotension

**Provenance:** V6.A multi-head DTI (target engagement) · V6.B Cluster D θ̄ (cognition relevance) · V7 / disease-conditioned class prior (effect size) · Gap-3 mechanism-class track record (prognosis) · 44-target off-target liability panel.

**Caveats & failure modes:**
- Predicted cognition effect is bounded by the Roberts 2020 ceiling (healthy-adult SMD ≈ 0.2-0.5); disease effects can be larger.
- Off-target liability flags are MODEL-PREDICTED (MAMMAL DTI), unvalidated — confirm against the drug's known safety profile.

---

### memantine — for AD

| | |
|---|---|
| **Mechanism** | NMDA_modulator at GRIN2B |
| **Predicted cognition effect** | Hedges' g = **+0.29** (90% CrI +0.14 to +0.43) |
| **Evidence quality (GRADE)** | 🟢 HIGH |
| **Mechanism-class track record** | SUCCESS (1 class drugs, 4 pooled RCTs) |
| **Class exemplars** | Memantine |
| **Effect basis** | this compound’s own pivotal trial |

**Why this GRADE rating:** RCT meta-analytic base (k=4) → start HIGH.

**Predicted off-target liability flags** (model-based, unvalidated):
- `CHRM3` (pKd≈7.1) — M3: anticholinergic (dry mouth, GI, urinary)
- `TACR1` (pKd≈6.9) — NK1: emetic / CNS off-target
- `OPRM1` (pKd≈6.8) — µ-opioid: dependence / respiratory depression
- `GABRA1` (pKd≈6.8) — GABA-A α1: sedation / dependence

**Provenance:** V6.A multi-head DTI (target engagement) · V6.B Cluster D θ̄ (cognition relevance) · V7 / disease-conditioned class prior (effect size) · Gap-3 mechanism-class track record (prognosis) · 44-target off-target liability panel.

**Caveats & failure modes:**
- Predicted cognition effect is bounded by the Roberts 2020 ceiling (healthy-adult SMD ≈ 0.2-0.5); disease effects can be larger.
- Off-target liability flags are MODEL-PREDICTED (MAMMAL DTI), unvalidated — confirm against the drug's known safety profile.

---

### methylphenidate — for ADHD

| | |
|---|---|
| **Mechanism** | catecholaminergic at SLC6A3 |
| **Predicted cognition effect** | Hedges' g = **+0.50** (90% CrI +0.10 to +0.32) |
| **Evidence quality (GRADE)** | 🟡 MODERATE |
| **Mechanism-class track record** | SUCCESS (3 class drugs, 1 pooled RCTs) |
| **Class exemplars** | Lisdexamfetamine, dextroamphetamine, methylphenidate |
| **Effect basis** | this compound’s own pivotal trial |

**Why this GRADE rating:** limited trial base (k=1); imprecision: wide CI / few RCTs (−1).

**Predicted off-target liability flags** (model-based, unvalidated):
- `CHRM3` (pKd≈7.2) — M3: anticholinergic (dry mouth, GI, urinary)
- `TACR1` (pKd≈6.9) — NK1: emetic / CNS off-target
- `OPRM1` (pKd≈6.8) — µ-opioid: dependence / respiratory depression
- `GABRA1` (pKd≈6.8) — GABA-A α1: sedation / dependence

**Provenance:** V6.A multi-head DTI (target engagement) · V6.B Cluster D θ̄ (cognition relevance) · V7 / disease-conditioned class prior (effect size) · Gap-3 mechanism-class track record (prognosis) · 44-target off-target liability panel · Gap-4 allosteric LTR (binding-reliability flag).

**Caveats & failure modes:**
- Predicted cognition effect is bounded by the Roberts 2020 ceiling (healthy-adult SMD ≈ 0.2-0.5); disease effects can be larger.
- MAMMAL's sequence-only binding score is structurally blind at this (allosteric/transporter) target; engagement is uncertain (see Gap-4 allosteric learn-to-rank).
- Off-target liability flags are MODEL-PREDICTED (MAMMAL DTI), unvalidated — confirm against the drug's known safety profile.

---

### pitolisant — for narcolepsy

| | |
|---|---|
| **Mechanism** | H3_cognition at HRH3 |
| **Predicted cognition effect** | Hedges' g = **+0.61** (90% CrI +0.30 to +0.92) |
| **Evidence quality (GRADE)** | 🟡 MODERATE |
| **Mechanism-class track record** | SUCCESS (1 class drugs, 2 pooled RCTs) |
| **Class exemplars** | Pitolisant |
| **Effect basis** | this compound’s own pivotal trial |

**Why this GRADE rating:** limited trial base (k=2); imprecision: wide CI / few RCTs (−1).

**Predicted off-target liability flags** (model-based, unvalidated):
- `CHRM3` (pKd≈7.4) — M3: anticholinergic (dry mouth, GI, urinary)
- `TACR1` (pKd≈7.0) — NK1: emetic / CNS off-target
- `OPRM1` (pKd≈7.0) — µ-opioid: dependence / respiratory depression
- `GABRA1` (pKd≈6.8) — GABA-A α1: sedation / dependence

**Provenance:** V6.A multi-head DTI (target engagement) · V6.B Cluster D θ̄ (cognition relevance) · V7 / disease-conditioned class prior (effect size) · Gap-3 mechanism-class track record (prognosis) · 44-target off-target liability panel.

**Caveats & failure modes:**
- Predicted cognition effect is bounded by the Roberts 2020 ceiling (healthy-adult SMD ≈ 0.2-0.5); disease effects can be larger.
- Off-target liability flags are MODEL-PREDICTED (MAMMAL DTI), unvalidated — confirm against the drug's known safety profile.

---

### idalopirdine — for AD

| | |
|---|---|
| **Mechanism** | 5HT6_antagonist at P50406 |
| **Predicted cognition effect** | Hedges' g = **-0.05** (90% CrI -0.16 to +0.05) |
| **Evidence quality (GRADE)** | 🟢 HIGH |
| **Mechanism-class track record** | FAILURE (3 class drugs, 3 pooled RCTs) |
| **Class exemplars** | Idalopirdine, SUVN-502, intepirdine |
| **Effect basis** | this compound’s own pivotal trial |

**Why this GRADE rating:** RCT meta-analytic base (k=3) → start HIGH.

**Provenance:** V6.A multi-head DTI (target engagement) · V6.B Cluster D θ̄ (cognition relevance) · V7 / disease-conditioned class prior (effect size) · Gap-3 mechanism-class track record (prognosis) · Gap-4 allosteric LTR (binding-reliability flag).

**Caveats & failure modes:**
- Predicted cognition effect is bounded by the Roberts 2020 ceiling (healthy-adult SMD ≈ 0.2-0.5); disease effects can be larger.
- MAMMAL's sequence-only binding score is structurally blind at this (allosteric/transporter) target; engagement is uncertain (see Gap-4 allosteric learn-to-rank).
- ⚠ This mechanism class has a NEGATIVE pivotal-trial track record in this indication — proceed only with a hypothesis that distinguishes this compound from its failed class-mates.

---

### encenicline — for CIAS

| | |
|---|---|
| **Mechanism** | alpha7_nAChR at CHRNA7 |
| **Predicted cognition effect** | Hedges' g = **+0.00** (90% CrI -0.10 to +0.10) |
| **Evidence quality (GRADE)** | 🟢 HIGH |
| **Mechanism-class track record** | FAILURE (5 class drugs, 2 pooled RCTs) |
| **Class exemplars** | ABT-126, DMXB-A, Encenicline, TC-5619, Varenicline |
| **Effect basis** | this compound’s own pivotal trial |

**Why this GRADE rating:** limited trial base (k=2).

**Predicted off-target liability flags** (model-based, unvalidated):
- `CHRM3` (pKd≈7.2) — M3: anticholinergic (dry mouth, GI, urinary)
- `TACR1` (pKd≈7.0) — NK1: emetic / CNS off-target
- `OPRM1` (pKd≈6.9) — µ-opioid: dependence / respiratory depression
- `ADRA1D` (pKd≈6.9) — α1D: orthostatic hypotension

**Provenance:** V6.A multi-head DTI (target engagement) · V6.B Cluster D θ̄ (cognition relevance) · V7 / disease-conditioned class prior (effect size) · Gap-3 mechanism-class track record (prognosis) · 44-target off-target liability panel · Gap-4 allosteric LTR (binding-reliability flag).

**Caveats & failure modes:**
- Predicted cognition effect is bounded by the Roberts 2020 ceiling (healthy-adult SMD ≈ 0.2-0.5); disease effects can be larger.
- MAMMAL's sequence-only binding score is structurally blind at this (allosteric/transporter) target; engagement is uncertain (see Gap-4 allosteric learn-to-rank).
- ⚠ This mechanism class has a NEGATIVE pivotal-trial track record in this indication — proceed only with a hypothesis that distinguishes this compound from its failed class-mates.
- Off-target liability flags are MODEL-PREDICTED (MAMMAL DTI), unvalidated — confirm against the drug's known safety profile.

---

### galantamine — for AD

| | |
|---|---|
| **Mechanism** | AChE_inhibitor at ACHE |
| **Predicted cognition effect** | Hedges' g = **+0.37** (90% CrI +0.31 to +0.42) |
| **Evidence quality (GRADE)** | 🟢 HIGH |
| **Mechanism-class track record** | SUCCESS (3 class drugs, 12 pooled RCTs) |
| **Class exemplars** | Donepezil, Galantamine, Rivastigmine |
| **Effect basis** | this compound’s own pivotal trial |

**Why this GRADE rating:** RCT meta-analytic base (k=12) → start HIGH.

**Predicted off-target liability flags** (model-based, unvalidated):
- `CHRM3` (pKd≈7.2) — M3: anticholinergic (dry mouth, GI, urinary)
- `TACR1` (pKd≈7.0) — NK1: emetic / CNS off-target
- `OPRM1` (pKd≈6.9) — µ-opioid: dependence / respiratory depression
- `GABRA1` (pKd≈6.8) — GABA-A α1: sedation / dependence

**Provenance:** V6.A multi-head DTI (target engagement) · V6.B Cluster D θ̄ (cognition relevance) · V7 / disease-conditioned class prior (effect size) · Gap-3 mechanism-class track record (prognosis) · 44-target off-target liability panel.

**Caveats & failure modes:**
- Predicted cognition effect is bounded by the Roberts 2020 ceiling (healthy-adult SMD ≈ 0.2-0.5); disease effects can be larger.
- Off-target liability flags are MODEL-PREDICTED (MAMMAL DTI), unvalidated — confirm against the drug's known safety profile.

---

Generated by `scripts/80_clinician_dossier.py` via `reporting/clinician_dossier.py`.