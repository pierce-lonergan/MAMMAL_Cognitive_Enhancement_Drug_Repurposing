# Pareto-Restructured Shortlist v1 (§8.0a)

Five-axis non-dominated sort over the V7 wet-lab handoff. Every compound on **frontier rank 0** is non-dominated — there is no compound that beats it on every axis.

**Axes** (all directions = HIGHER is better):
- `efficacy_rrf` — v7 RRF score from 5-cluster fusion (MAMMAL+Tanimoto+Boltzina+ADMET+MoA)
- `safety_neg_liability` — −(n_tier_1 + 0.5×n_tier_2) from §8.0b-zn liability panel
- `selectivity_pi` — top-target Cheng 2010 Partition Index from §7.4 v2
- `ip_novelty` — CTgov §8.3 IP-status mapped to [0,1]: approved=0.0, investigational=0.4, early=0.7, none=1.0
- `scaffold_novelty` — §8.10 nootropic-similarity tag: analog=0.0, intermediate=0.5, novel_scaffold=1.0

**Pareto front (rank 0)**: 12 compounds
**Hypervolume (MC, 200k samples)**: 7.7138
**Frontier count**: 193 (max rank = 192)

## Pareto front (rank 0) — top 30 by RRF

| # | Compound | Tier | RRF | Safety | PI | IP | Novelty |
|---|---|---|---|---|---|---|---|
| 1 | d-amphetamine | positive_control | 1.185 | -0.0 | 0.47 | investigational | `novel_scaffold` |
| 3 | bupropion | extended_cns | 1.160 | -0.0 | 0.71 | approved | `novel_scaffold` |
| 4 | aniracetam | named_in_research | 1.122 | -0.0 | 0.50 | none | `intermediate` |
| 6 | pridopidine | named_in_research | 1.085 | -0.0 | 0.88 | investigational | `novel_scaffold` |
| 7 | levetiracetam | extended_cns | 1.084 | -0.0 | 0.97 | approved | `intermediate` |
| 14 | cx-516 | named_in_research | 1.057 | -0.0 | 0.88 | investigational | `novel_scaffold` |
| 20 | chembl1255723 | chembl_expanded | 1.023 | -0.0 | 0.92 | none | `novel_scaffold` |
| 21 | pramipexole | extended_cns | 1.023 | -0.0 | 0.97 | approved | `novel_scaffold` |
| 23 | cep-26401 | named_in_research | 1.020 | -0.0 | 0.98 | early | `novel_scaffold` |
| 24 | chembl302231 | chembl_expanded | 1.018 | -0.0 | 0.96 | none | `novel_scaffold` |
| 46 | riluzole | named_in_research | 0.997 | -0.0 | 0.99 | none | `novel_scaffold` |
| 53 | chembl294061 | chembl_expanded | 0.982 | -0.0 | 1.00 |  | `novel_scaffold` |

## Rank-1 frontier (one peel below Pareto) — top 15

| # | Compound | Tier | RRF | Safety | PI | IP | Novelty |
|---|---|---|---|---|---|---|---|
| 2 | methylphenidate | positive_control | 1.163 | -0.0 | 0.45 | approved | `novel_scaffold` |
| 5 | rasagiline | extended_cns | 1.092 | -0.0 | 0.71 | approved | `novel_scaffold` |
| 9 | lanicemine | named_in_research | 1.069 | -0.0 | 0.83 | investigational | `intermediate` |
| 10 | pramiracetam | extended_cns | 1.067 | -0.5 | 0.41 | none | `intermediate` |
| 11 | cx-717 | named_in_research | 1.066 | -0.0 | 0.78 | investigational | `novel_scaffold` |
| 22 | chembl1256414 | chembl_expanded | 1.020 | -0.0 | 0.92 | none | `novel_scaffold` |
| 25 | chembl1256378 | chembl_expanded | 1.018 | -0.0 | 0.93 | none | `novel_scaffold` |
| 38 | risperidone | extended_cns | 1.005 | -11.0 | 0.95 | approved | `novel_scaffold` |
| 39 | chembl292558 | chembl_expanded | 1.004 | -0.0 | 0.96 | none | `novel_scaffold` |
| 59 | aripiprazole | positive_control | 0.975 | -9.0 | 0.98 |  | `novel_scaffold` |
| 66 | buspirone | extended_cns | 0.971 | -0.0 | 0.99 |  | `novel_scaffold` |
| 114 | huperzine A | extended_cns | 0.801 | -0.0 | 1.00 |  | `novel_scaffold` |

## Frontier-size distribution

| Rank | N compounds |
|---|---|
| 0 | 12 |
| 1 | 12 |
| 2 | 11 |
| 3 | 6 |
| 4 | 12 |
| 5 | 10 |
| 6 | 9 |
| 7 | 9 |
| 8 | 5 |
| 9 | 5 |
| 10 | 6 |
| 11 | 4 |
| 12 | 5 |
| 13 | 3 |
| 14 | 4 |
| 15 | 2 |
| 16 | 1 |
| 17 | 1 |
| 18 | 1 |
| 19 | 1 |
| 20 | 1 |
| 21 | 1 |
| 22 | 1 |
| 23 | 1 |
| 24 | 1 |
| 25 | 1 |
| 26 | 1 |
| 27 | 1 |
| 28 | 1 |
| 29 | 1 |
| 30 | 1 |
| 31 | 1 |
| 32 | 1 |
| 33 | 1 |
| 34 | 1 |
| 35 | 1 |
| 36 | 1 |
| 37 | 2 |
| 38 | 1 |
| 39 | 1 |
| 40 | 1 |
| 41 | 1 |
| 42 | 1 |
| 43 | 1 |
| 44 | 1 |
| 45 | 1 |
| 46 | 1 |
| 47 | 1 |
| 48 | 1 |
| 49 | 1 |
| 50 | 1 |
| 51 | 1 |
| 52 | 1 |
| 53 | 2 |
| 54 | 1 |
| 55 | 1 |
| 56 | 1 |
| 57 | 1 |
| 58 | 2 |
| 59 | 1 |
| 60 | 1 |
| 61 | 1 |
| 62 | 2 |
| 63 | 1 |
| 64 | 1 |
| 65 | 1 |
| 66 | 1 |
| 67 | 1 |
| 68 | 2 |
| 69 | 1 |
| 70 | 1 |
| 71 | 1 |
| 72 | 1 |
| 73 | 1 |
| 74 | 1 |
| 75 | 1 |
| 76 | 1 |
| 77 | 1 |
| 78 | 1 |
| 79 | 1 |
| 80 | 1 |
| 81 | 1 |
| 82 | 1 |
| 83 | 1 |
| 84 | 1 |
| 85 | 2 |
| 86 | 1 |
| 87 | 1 |
| 88 | 1 |
| 89 | 1 |
| 90 | 1 |
| 91 | 1 |
| 92 | 1 |
| 93 | 1 |
| 94 | 1 |
| 95 | 1 |
| 96 | 1 |
| 97 | 1 |
| 98 | 1 |
| 99 | 1 |
| 100 | 1 |
| 101 | 1 |
| 102 | 1 |
| 103 | 1 |
| 104 | 1 |
| 105 | 1 |
| 106 | 1 |
| 107 | 1 |
| 108 | 1 |
| 109 | 1 |
| 110 | 1 |
| 111 | 1 |
| 112 | 1 |
| 113 | 1 |
| 114 | 1 |
| 115 | 1 |
| 116 | 1 |
| 117 | 1 |
| 118 | 1 |
| 119 | 1 |
| 120 | 1 |
| 121 | 1 |
| 122 | 1 |
| 123 | 1 |
| 124 | 1 |
| 125 | 1 |
| 126 | 1 |
| 127 | 1 |
| 128 | 1 |
| 129 | 1 |
| 130 | 1 |
| 131 | 1 |
| 132 | 1 |
| 133 | 1 |
| 134 | 1 |
| 135 | 1 |
| 136 | 1 |
| 137 | 1 |
| 138 | 1 |
| 139 | 1 |
| 140 | 1 |
| 141 | 1 |
| 142 | 1 |
| 143 | 1 |
| 144 | 1 |
| 145 | 1 |
| 146 | 1 |
| 147 | 1 |
| 148 | 1 |
| 149 | 1 |
| 150 | 1 |
| 151 | 1 |
| 152 | 1 |
| 153 | 1 |
| 154 | 1 |
| 155 | 1 |
| 156 | 1 |
| 157 | 1 |
| 158 | 1 |
| 159 | 1 |
| 160 | 1 |
| 161 | 1 |
| 162 | 1 |
| 163 | 1 |
| 164 | 1 |
| 165 | 1 |
| 166 | 1 |
| 167 | 1 |
| 168 | 1 |
| 169 | 1 |
| 170 | 1 |
| 171 | 1 |
| 172 | 1 |
| 173 | 1 |
| 174 | 1 |
| 175 | 1 |
| 176 | 1 |
| 177 | 1 |
| 178 | 1 |
| 179 | 1 |
| 180 | 1 |
| 181 | 1 |
| 182 | 1 |
| 183 | 1 |
| 184 | 1 |
| 185 | 1 |
| 186 | 1 |
| 187 | 1 |
| 188 | 1 |
| 189 | 1 |
| 190 | 1 |
| 191 | 1 |
| 192 | 1 |

## Interpretation

- **Front (rank 0)** is the *production* wet-lab handoff for constrained-budget runs. Every compound here strictly dominates nothing on the front.
- **Crowding distance** (saved as `crowding_distance` column in the parquet) within rank 0 ranks compounds by frontier diversity — pick the top-k by crowding when the front is too large.
- **Rank 1-2** are useful when expanding the candidate pool: these are compounds beaten on exactly one or two axes.
- **Hypervolume** measures the volume of objective space dominated by the front; higher = better coverage. Use as a regression metric when re-running with new evidence.

---

Generated by `scripts/42_v5_pareto_shortlist.py`. Input ranking: `final_ranking_v7_moa.parquet`. PASS-only=False.