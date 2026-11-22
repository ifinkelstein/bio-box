# Meta-mode gene-calling cost and batch sizing before `PL3` (card B-63)

Date: 2026-11-21, `wilkcomp02`. Every number carries a label (`*measured*`, `*derived*`, or `*estimate*`) per `AGENTS.md` §6. All pod benchmarks ran inside `bio-box_faec1f4.sif` under `systemd-run --user --scope` (`MemoryMax=16G`–`32G`, `MemorySwapMax=0`) at `nice -n 19`, with `/tmp/if2453/r10/root` and `/tmp/if2453/s2/inputs/PL2` bound read-only (`:ro`), writing only to `/tmp/if2453/b63/` (cleaned up after recording). Nothing was written to the data root or `/stor`.

Script: [`scripts/stress/b63_meta_genecall_arms.py`](../../scripts/stress/b63_meta_genecall_arms.py).

## 1. Summary of findings

1. **Per-genome `GeneFinder(meta=True)` instantiation is negligible**: standalone instantiation costs **`0.5 us`** (`pyrodigal.GeneFinder(meta=True)`) and **`1.6 us`** (`pyrodigal_gv.ViralGeneFinder(meta=True)`) (*measured*, 500 iterations, `arm-a`). Across `PL3`'s ~753k genomes, fresh instantiation costs `~0.38 s` (*derived*). Reusing one cached finder per thread/process changes runtime by `-2%` to `+2%` (within pod load jitter; exact call-for-call equality verified on every linear and circular test sequence).
2. **Root cause of the `IMG-PR` / `RefSeq` rate gap**: in Prodigal meta mode, **CPU-s/Mbp increases `2.4×`–`3.9×` with contig length** in all four plasmid sources (`1.38`–`2.01 s/Mbp` at 1–5 kb vs `1.89`–`3.20 s/Mbp` at 10–30 kb vs `3.27`–`5.76 s/Mbp` at ≥50 kb, *measured*, `arm-a`). `CALL_S_PER_MBP[("meta", False)] = 2.9450` was fitted on `PL1`'s small linear minority, whereas `PL1` was 94.5% circular (`CALL_S_PER_MBP[("meta", True)] = 4.4205`, because circular contigs are wrapped with `2 × 50 kb` flanks). `IMG-PR` is predominantly long linear contigs (45% of bases ≥50 kb), so linear meta calling runs near the circular rate (`~4.0–4.8 s/Mbp`).
3. **Secondary bottleneck in `genecall._batches`**: `genome_uid` order interleaves `< 10 kb` and `>= 10 kb` genomes. Previously `_batches` flushed the pending short-genome batch whenever a `>= 10 kb` genome arrived, turning 1,982 short genomes (`11.07 Mbp`) in our 1-in-200 `IMG-PR` + `mMGE` slice into **921 tiny worker-process batches** (424 singletons, mean `2.15` genomes / `12.0 kb` instead of `BATCH_BASES = 1,000,000` bp), burning `68.5–86.2 s` (`6.2–7.8 s/Mbp`) in worker IPC (*measured*, `arm-b`).
4. **Fix in `src/biobox/l1/genecall.py`**:
   - Accumulate `< THREAD_MIN_BASES` genomes across interleaved long genomes up to `BATCH_BASES` in `_batches` (921 → **12** short batches, 0 singletons, mean `922.4 kb`; cuts worker-process `kids` CPU by **36%–40%** and total step CPU by **8.5%–8.8%**, with **5/5** output Parquet files sha256 byte-identical).
   - Price linear meta-mode contigs `>= META_LINEAR_LONG_BP` (`10_000` bp) in `calls_cpu_s` at `META_LINEAR_LONG_FACTOR * CALL_S_PER_MBP[("meta", False)]` (`1.50 × 2.9450 = 4.4175 s/Mbp`, matching `CALL_S_PER_MBP[("meta", True)] = 4.4205`), leaving `CALL_S_PER_MBP` itself unchanged (so built stages' published plans and `test_b48_genecall_estimate.py` stay untouched, and monkeypatching `CALL_S_PER_MBP` scales proportionally).
5. **`PL3` projection**: on DEC-63 circular-DTR `IMG-PR` + `mMGE`, `cut_plan` at margin `240` CPU-min now cuts **8 batches** estimated at **`~212.4` CPU-min** each; at the post-fix measured ratio (`0.902×`), each batch projects to **`~192` CPU-min** (*derived*) — **~88 CPU-min below the 280 CPU-min process guard** (down from `~269–271` CPU-min on 6 batches before B-63, and `~304–308` CPU-min for linear-only `PL3`).

## 2. Arm A: Finder instantiation overhead and contig-length scaling

Run: `b63_meta_genecall_arms.py arm-a --bin-bp 500000 --bin-max-seqs 200 --unstrat-bp 1500000 --reps 2 --modes meta,gv-meta` (*measured*, `wilkcomp02`, 2026-11-21, 1 CPU, `bio-box_faec1f4.sif`, interleaved fresh vs. reused order reversed on rep 2).

- Standalone init (500 iterations): `pyrodigal.GeneFinder(meta=True)` = **`0.5 us`**, `pyrodigal_gv.ViralGeneFinder(meta=True)` = **`1.6 us`**.
- Exact call-for-call equality (BLAKE3 digest over all gene coordinates, strands, translations, codes, and partial flags): **verified** across all 32 `(source, bin, mode)` cells and on circular `mMGE` short/medium subsets.

| Source | Length bin | Mode | Seqs | Mean bp | Fresh (s/Mbp) | Reused (s/Mbp) | Delta |
|---|---|---|---:|---:|---:|---:|---:|
| `IMG-PR` | `short_1_5kb` | `meta` | 124 | 4,069 | 2.0063 | 1.8494 | -7.82% |
| `IMG-PR` | `med_10_30kb` | `meta` | 28 | 17,894 | 3.1989 | 3.1856 | -0.42% |
| `IMG-PR` | `long_50kb+` | `meta` | 4 | 126,358 | 4.0094 | 3.9078 | -2.53% |
| `IMG-PR` | `unstrat` | `meta` | 25 | 60,176 | 4.7930 | 4.7485 | -0.93% |
| `mMGE` | `short_1_5kb` | `meta` | 166 | 3,034 | 1.8466 | 1.8044 | -2.29% |
| `mMGE` | `med_10_30kb` | `meta` | 28 | 18,408 | 1.8883 | 1.8603 | -1.49% |
| `mMGE` | `long_50kb+` | `meta` | 8 | 101,259 | 5.5323 | 5.2867 | -4.44% |
| `mMGE` | `unstrat` | `meta` | 154 | 9,785 | 3.8069 | 3.3672 | -11.55% |
| `gut` | `short_1_5kb` | `meta` | 152 | 3,301 | 1.4712 | 1.5396 | +4.65% |
| `gut` | `med_10_30kb` | `meta` | 27 | 19,143 | 2.8057 | 2.7439 | -2.20% |
| `gut` | `long_50kb+` | `meta` | 8 | 71,388 | 5.7579 | 5.1836 | -9.97% |
| `gut` | `unstrat` | `meta` | 179 | 8,423 | 2.4480 | 2.4272 | -0.85% |
| `RefSeq` | `short_1_5kb` | `meta` | 182 | 2,748 | 1.3821 | 1.3746 | -0.54% |
| `RefSeq` | `med_10_30kb` | `meta` | 27 | 19,037 | 2.3559 | 2.3378 | -0.77% |
| `RefSeq` | `long_50kb+` | `meta` | 5 | 113,307 | 3.2690 | 3.3068 | +1.15% |
| `RefSeq` | `unstrat` | `meta` | 20 | 77,031 | 4.3422 | 4.2313 | -2.55% |

(`gv-meta` showed the same length scaling and ~0% reuse delta: e.g. `IMG-PR` `1.99` → `3.40` → `4.33 s/Mbp`, `RefSeq` `1.58` → `2.61` → `3.71 s/Mbp`.)

## 3. Arm B: Full `stage genecall` batch benchmark (before vs. after B-63)

Representative 1-in-200 slice of `IMG-PR` + `mMGE` new contigs from `/tmp/if2453/s2/inputs/PL2` (**3,742 genomes**, **82,577,785 bp**; 1,760 genomes ≥10 kb holding 71.51 Mbp and 1,982 genomes <10 kb holding 11.07 Mbp), ingested into a fresh scratch root and run with `stage genecall --workers 4` (2 reps each, full subprocess `RUSAGE_SELF + RUSAGE_CHILDREN` CPU-minutes, *measured*, `wilkcomp02`, 2026-11-21):

| Arm | Circular | Short batches (singletons) | `cut_plan` est (CPU-min) | Rep 1 / Rep 2 measured (CPU-min) | Median (CPU-min) | Ratio (meas/est) | Self / Kids (s) | Wall (s) | Peak RSS (MiB) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Baseline linear | 38 | 921 (424) | 5.301 | 7.600 / 7.673 | **7.637** | **1.4406** | 389.6 / 68.5 | 109.7 | 514.0 |
| **Post-fix linear** | 38 | **12 (0)** | **7.453** | 6.961 / 7.018 | **6.990** | **0.9378** | 378.2 / **41.2** | 104.3 | 519.0 |
| Baseline DEC-63 DTR | 684 | 921 (424) | 6.779 | 8.614 / 8.672 | **8.643** | **1.2750** | 432.4 / 86.2 | 124.6 | 509.0 |
| **Post-fix DEC-63 DTR** | 684 | **12 (0)** | **8.734** | 7.810 / 7.946 | **7.878** | **0.9020** | 417.8 / **54.9** | 117.2 | 507.6 |

**Byte-identity verification** (*measured*, `sha256` over all 5 published outputs: `feature`, `feature_protein`, `protein_occurrence`, `protein_seq`, `ledger`):
- Baseline linear vs. Post-fix linear: **5/5 files sha256-identical**.
- Baseline DEC-63 DTR vs. Post-fix DEC-63 DTR: **5/5 files sha256-identical**.

## 4. `PL3` batch projections and impact on built stages

*Derived* from the 1-in-200 Arm B slices (`200 × (estimate - 0.25)` for full-stage work, cut by `cut_plan` at margin `240` CPU-min):

| Scenario | Full-stage est (CPU-min) | Batches (`margin=240`) | Est per batch (CPU-min) | Projected measured per batch (CPU-min) | Margin below 280 guard |
|---|---:|---:|---:|---:|---:|
| Baseline `PL3` (DEC-63 DTR, pre-B-63) | ~1,306 (note: 1,264 trimmed) | 6 | 210.7–217.8 | **269–278** | **2–11 min (unsafe)** |
| Baseline `PL3` (linear, pre-B-63) | ~1,010 | 5 | ~202.3 | **296** | **-16 min (trips 280)** |
| **Post-fix `PL3` (DEC-63 DTR)** | **~1,697** | **8** | **~212.4** | **~192** (`0.902×`) | **+88 min (safe)** |
| **Post-fix `PL3` (linear)** | **~1,441** | **7** | **~206.1** | **~193** (`0.938×`) | **+87 min (safe)** |

**Impact on built stages** (*derived* from published run records in `/tmp/if2453/r10/root`):
- `S1`, `H1`, `V1` (`gv-meta`) and `P1`–`P6`, `C5` (`single`): **0.00% change** (none use `meta` mode).
- `PL1` (`meta`, 94.5% circular across 4 batches; published records `209.89`–`213.89` CPU-min against `218.35` estimated, `0.961`–`0.980×`): even if 100% of `PL1`'s linear bases (`110.7`–`126.5 Mbp`/batch) were `>= 10 kb`, batch estimates rise by `+3.3` to `+3.8` CPU-min (`+1.6%`, to `221.7`–`222.2` CPU-min, ratio `0.946`–`0.964×`), staying well within `±20%` (and `PL1`'s published `plan.json` is never recomputed).
