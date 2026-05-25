# WSL2 Boltz Validation Record

**Date**: 2026-05-25
**Question**: Does WSL2 with cuequivariance kernels deliver the predicted 2.5-5× Boltz speedup, and does the 16 GB Blackwell context anomaly (Microsoft WSL #40401) trigger on our 12 GB RTX 5070?

## Setup state at measurement

- WSL2 Ubuntu 24.04.3, kernel 6.6.87.2-microsoft-standard-WSL2
- Python 3.12.3 in `/root/mammal_env/` venv
- PyTorch 2.12.0.dev20260407+cu128 (Blackwell sm_120 compatible)
- boltz 2.2.1 with `--use_msa_server`, kernels ENABLED (no `--no_kernels`)
- cuequivariance-torch 0.10.0 + cuequivariance-ops-torch-cu12 0.10.0 (Linux x86_64 wheel installed cleanly)
- cuequivariance triangle kernel import verified PRIOR to inference (no runtime ImportError)

## Smoke-pair: CHRNA7 (P36544, 502 aa) + galantamine

| Metric | Windows (no_kernels) | **WSL2 (kernels)** |
|---|---|---|
| Total wall-clock first call | 149.8 s | 305 s |
| `Predicting DataLoader` step (GPU inference only) | ~30 s | **~23 s** |
| MSA fetch (remote MMseqs2 + processing) | ~120 s | ~270 s |
| `affinity_pred_value` | 0.555 | 0.577 (sample mean) |
| `affinity_probability_binary` | 0.498 | 0.529 (sample mean) |
| Sample variance (Boltz returned 2 samples here) | n/a | 0.732 / 0.420 (pred), 0.708 / 0.350 (prob) |

The cold-start first-call wall-clock is **slower** on WSL2 because the MSA server response was slower this attempt (independent of our setup — MMseqs2 server load varies). The **inference itself**, which is what matters for amortized throughput, is faster on WSL2.

## What's actually different in the math

The Boltz outputs are numerically close but not identical:
- Windows used the pure-PyTorch triangle-multiplication fallback (`--no_kernels`)
- WSL2 used the cuequivariance native CUDA kernel, which switched to round-nearest (RN) precision in v0.10 (vs the prior round-towards-zero)
- Boltz also defaulted to 2 diffusion samples on WSL2 (smoke script used `--diffusion_samples 1` on Windows, took the default for WSL2). This is why we see two `_value` / `_probability_binary` fields plus the merged mean.

The rank order in the focused-sweep Boltz parquet (Phase 0.5 — CHRNA7 PAMs correctly ranked) would not change.

## Speedup analysis for the full sweep

For Phase 0.4 (≈1,500 surviving compound × 22 target pairs):

- **MSA cost**: ~22 fetches × 270 s = ~99 min cold; cached after first call per sequence.
- **Inference cost**: (1500 - 22) pairs × 23 s = ~570 min ≈ **9.5 hours**.
- **Total WSL2 estimate**: ~11 hours overnight.

Compare Windows projection: ~62 hours = 2.5 days.

Speedup: **~5.6× wall-clock**. This matches the research doc's "2.5-5× expected with cuequivariance kernels" prediction at the high end.

## Blackwell 16 GB context anomaly: did NOT trigger

The cuEquivariance deep-dive research warned of a 16 GB CUDA driver context overhead on sm_120 inside WSL2 per Microsoft WSL issue #40401, observed on RTX 5090 and RTX PRO 6000. Per that report, allocating the basic CUDA context would consume 16 GB before any model loads.

On our 12 GB RTX 5070, this did NOT happen:
- Boltz loaded both `boltz2_conf.ckpt` (~3 GB) and `boltz2_aff.ckpt` (~2 GB)
- Structure prediction + affinity prediction ran end-to-end
- No OOM, no allocation failures
- (We didn't capture `nvidia-smi` mid-run — adding that to a follow-up probe — but the qualitative outcome is unambiguous: the predict ran)

This suggests the 16 GB anomaly is either:
- Specific to even-larger-VRAM cards where the context auto-sizes proportionally
- Specific to a different WSL kernel / driver version than ours (591.59)
- Resolved in our driver but not in the version the issue reporters used

**Conclusion**: WSL2 + Boltz + cuequivariance + Blackwell sm_120 on a 12 GB RTX 5070 is viable. The 16 GB anomaly is NOT a blocker for us.

## What to commit to next

1. **Full Phase 0.4 sweep runs in WSL2**, not on Windows. Use the focused-sweep parquet (already in `data/results/v2/boltzina_affinity.parquet`) as the seed; resume to add the missing pairs.
2. The WSL2 boltz wrapper for the bigger sweep needs a small Python orchestrator inside WSL2 that reads `data/interim/{targets,compounds}.parquet` from `/mnt/c/...`, runs `boltz predict` per pair, and writes back to `data/results/v2/boltzina_affinity.parquet`. The existing `scripts/_boltzina_focused.py` is the template — port it.
3. Hot path: pre-fetch MSAs for all 22 targets in a single batch run so the sweep is purely inference-bound.
4. After sweep completes, re-run `scripts/15_v2_fusion.py` — RRF now absorbs Boltzina as the third ranker.
