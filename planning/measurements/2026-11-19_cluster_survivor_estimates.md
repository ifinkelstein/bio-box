# Refitting intermediate cascade (`input_step1` / `input_step2`) and `kmermatcher` estimates (wilkcomp02, 2026-11-19)

Card B-54 (DEC-53, DEC-55; AGENTS §6).

## 1. Cause and scope

On release `r26.11.1` (`49,102,041` distinct proteins across `S1,P1..P6,H1,PL1`, `+7.17%` over `rM.2`'s `45,815,843`, *measured*, wilkcomp02, `planning/measurements/2026-11-16_release_r26.11.1.md`), two clustering estimates missed by more than 20%:

1. **Intermediate cascade survivor growth (DEC-55).** Adding `3,286,198` divergent phage (`H1`) and plasmid (`PL1`) proteins left `9,421,654` surviving sequences in NR50's final `-s 3` `input_step2` database (*measured*, `r26.11.1` `cluster.nr50.19.prefilter.cut.json`), `+13.34%` over `rM.2`'s `8,312,546` (*derived*: `9,421,654 / 8,312,546 - 1`), nearly double the `+7.17%` release-total protein growth. Because `call_cpu_min` scaled per-unit prefilter work only by release-total `(proteins / CALL_N0) ** g` without reading the surviving target count of `input_step1` / `input_step2`, `cluster.nr50.19.prefilter` cost `1,373.84` CPU-minutes (`55.18` to `60.15` per part across 24 parts, *measured*) against `817.52` estimated by the pre-B-54 code (`1.680×`, *derived*; `1.224×` to `1.334×` even after `c35_est.py`'s `rM.2` calibration factor).
2. **`kmermatcher` 48-vs-18 thread scaling (DEC-53).** `cluster.nr50.00.kmermatcher` on `r26.11.1` ran at `18` threads and cost `19.093` CPU-minutes (*measured*) against `31.407` projected by scaling `rM.2`'s `48`-thread record (`29.305` CPU-minutes, *measured*) by the protein ratio (`0.608×`, *derived*). On the identical `31,183,876` NR90 representatives of `rM.1` (32 threads) and `rM.2` (48 threads), `kmermatcher` cost `19.362` vs `29.305` CPU-minutes (`1.514×` OpenMP/NUMA overhead at 48 threads, *derived*).

## 2. Data sources

All measurements were read read-only on `wilkcomp02.ccbb.utexas.edu` (2026-11-19) from `/tmp/if2453/r10/root`:

- Published `l2/cluster_plan/release=<R>/*.cut.json` and `runs/stage=<R>/cluster.*.json` for `rN.3` (`10,620,684` proteins), `rM.1` (`45,815,843`), `rM.2` (`45,815,843`), and `r26.11.1` (`49,102,041`).
- `rN.1b` in `/tmp/if2453/r10/root` was clustered in one process (no `cluster_plan` records); its split-clustering counterpart over the same `S1,P1` stages (`rC.1b`, `2,053,075` proteins, wilkcomp02, 2026-11-11) is taken from `tests/r1/c22_cluster_calls_measured.tsv`.

## 3. Fitted models

### 3.1 Per-query cascade calls (`biobox.l2.cluster.call_cpu_min`)

Each split call's work scales by its own surviving target/query sequence count `targets` (already read on disk when `<call>.cut` runs from `c.positional[1].index` for prefilters and from the split input `.index` for `align`/`rescorediagonal`), falling back to `proteins / CALL_N0` only when `targets is None`:

```
work(call, targets, units) = r_call * (units / 1e9) * (targets / CALL_T0[call]) ** g_call
```

Reference target counts `CALL_T0[call]` are `rN.3`'s measured surviving sequence counts, and `(r_call, g_call)` are minimax-fitted across all 5 releases (`rC.1b`, `rN.3`, `rM.1`, `rM.2`, `r26.11.1`; the four tiny `rescorediagonal` calls at `0.35` to `2.63` CPU-min per part, never judged under DEC-43, keep their `rN.3` calibration):

| Call | `CALL_T0` (*measured*, `rN.3`) | Old `(r, g)` | New `(r, g)` (*derived*) |
|---|---|---|---|
| `cluster.nr90.07.align` | 8,546,571 | `(204.09, 0.19)` | `(206.35, 0.14)` |
| `cluster.nr90.11.prefilter` | 8,363,309 | `(26.20, 0.65)` | `(30.24, 0.83)` |
| `cluster.nr90.12.align` | 8,363,309 | `(366.92, 0.00)` | `(384.94, 0.03)` |
| `cluster.nr90.members` | 8,278,537 | `(306.82, 0.00)` | `(324.95, 0.02)` |
| `cluster.nr50.07.align` | 5,196,591 | `(360.72, 0.14)` | `(354.32, 0.07)` |
| `cluster.nr50.11.prefilter` | 4,210,724 | `(15.27, 0.42)` | `(16.83, 0.63)` |
| `cluster.nr50.12.align` | 4,210,724 | `(326.53, 0.00)` | `(361.73, 0.07)` |
| `cluster.nr50.15.prefilter` | 3,141,911 | `(66.72, 0.37)` | `(75.14, 0.63)` |
| `cluster.nr50.16.align` | 3,141,911 | `(242.64, 0.12)` | `(237.95, 0.06)` |
| `cluster.nr50.19.prefilter` | 2,999,008 | `(163.11, 0.39)` | `(194.50, 0.73)` |
| `cluster.nr50.20.align` | 2,999,008 | `(247.90, 0.21)` | `(225.85, 0.13)` |
| `cluster.nr50.members` | 2,966,200 | `(773.02, 0.09)` | `(764.81, 0.03)` |

### 3.2 `kmermatcher` (`biobox.l2.cluster.kmermatcher_cpu_min`)

Across all 7 records at `16`, `18`, and `32` threads, `kmermatcher` runs at `5.581` to `6.352` CPU-minutes per `10^7` input sequences (*derived*); at `48` threads (`rM.2 nr50`), OpenMP/NUMA contention raises it to `9.397` CPU-minutes per `10^7` sequences:

```
kmermatcher_cpu_min(sequences, threads) = 6.00 * (sequences / 1e7) * (1.0 + 0.035 * max(0, threads - 32))
```

| Release | Step | Input sequences (*measured*) | Threads (*measured*) | CPU-min (*measured*) | New est (*derived*) | Ratio (*derived*) |
|---|---|---:|---:|---:|---:|---:|
| `rN.3` | `cluster.nr90.00.kmermatcher` | 10,620,684 | 16 | 6.390 | 6.372 | 1.003 |
| `rN.3` | `cluster.nr50.00.kmermatcher` | 8,278,537 | 16 | 4.620 | 4.967 | 0.930 |
| `rM.1` | `cluster.nr90.00.kmermatcher` | 45,815,843 | 32 | 29.104 | 27.490 | 1.059 |
| `rM.1` | `cluster.nr50.00.kmermatcher` | 31,183,876 | 32 | 19.362 | 18.710 | 1.035 |
| `rM.2` | `cluster.nr90.00.kmermatcher` | 45,815,843 | 32 | 28.782 | 27.490 | 1.047 |
| `rM.2` | `cluster.nr50.00.kmermatcher` | 31,183,876 | 48 | 29.305 | 29.188 | 1.004 |
| `r26.11.1` | `cluster.nr90.00.kmermatcher` | 49,102,041 | 18 | 30.520 | 29.461 | 1.036 |
| `r26.11.1` | `cluster.nr50.00.kmermatcher` | 33,117,773 | 18 | 19.093 | 19.871 | 0.961 |

All 8 records are within `±7.0%` (*derived*).

## 4. Before-and-after ratios across all built releases

### 4.1 `r26.11.1` (`49,102,041` proteins, 24 parts)

| Call | Surviving seqs (*measured*) | Units (*measured*) | Measured CPU-min | Old est (*derived*) | Old ratio | New est (*derived*) | New ratio |
|---|---:|---:|---:|---:|---:|---:|---:|
| `cluster.nr90.07.align` | 35,153,162 | 1,412,500,763 | 340.41 | 385.61 | 0.883 | 355.29 | 0.958 |
| `cluster.nr90.11.prefilter` | 33,882,143 | 10,664,275,036 | 3,832.62 | 3,613.46 | 1.061 | 3,887.53 | 0.986 |
| `cluster.nr90.12.align` | 33,882,143 | 9,125,396,776 | 3,656.60 | 3,348.29 | 1.092 | 3,663.30 | 0.998 |
| `cluster.nr90.members` | 33,117,773 | 463,925,032 | 162.48 | 142.34 | 1.142 | 154.99 | 1.048 |
| `cluster.nr50.07.align` | 18,876,326 | 758,111,957 | 289.75 | 338.84 | 0.855 | 294.00 | 0.986 |
| `cluster.nr50.11.prefilter` | 14,912,484 | 4,385,758,277 | 1,433.59 | 1,385.10 | 1.035 | 1,421.43 | 1.009 |
| `cluster.nr50.12.align` | 14,912,484 | 3,333,344,123 | 1,337.46 | 1,088.44 | **1.229** | 1,317.37 | 1.015 |
| `cluster.nr50.15.prefilter` | 10,064,063 | 2,852,219,962 | 466.31 | 354.32 | **1.316** | 465.24 | 1.002 |
| `cluster.nr50.16.align` | 10,064,063 | 3,073,734,693 | 755.21 | 896.24 | 0.843 | 784.31 | 0.963 |
| `cluster.nr50.19.prefilter` | 9,421,654 | 2,698,588,280 | 1,373.84 | 817.52 | **1.680** | 1,228.31 | 1.118 |
| `cluster.nr50.20.align` | 9,421,654 | 2,964,472,779 | 743.50 | 1,013.59 | **0.734** | 776.96 | 0.957 |
| `cluster.nr50.members` | 9,291,312 | 440,098,571 | 356.51 | 390.47 | 0.913 | 348.32 | 1.023 |

On `r26.11.1`'s `cluster.nr50.19.prefilter`, all 24 parts (`55.18` to `60.15` CPU-min, *measured*) now sit at `1.078×` to `1.175×` of the `51.18` CPU-min per-part estimate (*derived*), including DEC-55 trigger parts `part-00000` (`56.741`, `1.109×`), `part-00001` (`55.825`, `1.091×`), and worst part `part-00010` (`60.154`, `1.175×`).

### 4.2 Summary across all 5 releases (`rC.1b`, `rN.3`, `rM.1`, `rM.2`, `r26.11.1`)

- **Substantive call totals (`60` release×call pairs across all 5 releases)**: `0` outside 20% (*derived*; worst `rN.3` `nr50.19.prefilter` at `0.856×` and `rC.1b` `nr50.19.prefilter` at `1.165×`).
- **Individual judged parts (`>= 24` CPU-min, DEC-43, across `rN.3`, `rM.1`, `rM.2`, `r26.11.1`)**: **`519` of `519` parts within 20%** (`0` outside 20%, *derived*).
- **Whole-release agreement (`test_cluster_call_estimates.py`)**: `rN.3` sum `2,272.1` vs `stage.cluster_cpu_min` `2,194.7` (`+3.5%`, *derived*); `rC.1b` sum `277.7` vs `284.1` (`-2.2%`, *derived*).

## 5. B-49 cluster reuse preservation

`_source_problem` accepts both the new target-scaled `call_cpu_min` estimate and the pre-B-54 `_legacy_call_cpu_min` estimate for `.cut.json` records (while continuing to reject any record whose `cpu_min` matches neither). Read-only verification on `/tmp/if2453/r10/root` (wilkcomp02, 2026-11-19):

- `cluster.reuse_source` over `rM.2`'s 7 stages (`S1,P1..P6`, `parts=24`, `one_process_cpu_min_max=20000.0`) -> `"rM.2"` (*measured*).
- `cluster.reuse_source` over `r26.11.1`'s 9 stages (`S1,P1..P6,H1,PL1`, `parts=24`, `one_process_cpu_min_max=20000.0`) -> `"r26.11.1"` (*measured*).
