# Card B-64: Release Clustering Calibration Factors Across `rM.2`, `r26.11.1`, and `r26.11.2`

**Date**: 2026-11-21
**Host**: `wilkcomp02` (read-only inspection of `/tmp/if2453/r10/root`; zero writes to root or `/stor`)
**Scope**: Card `B-64` (DEC-61; depends on `B-54`, `C-46`, `r26.11.2` clustering published)

## 1. Root Cause Summary

Inspecting all published cluster plan records (`l2/cluster_plan/release=<R>/*.cut.json`) and run records (`runs/stage=<R>/*.json`) across `rM.2` (`45,815,843` distinct proteins, *measured*), `r26.11.1` (`49,102,041` distinct proteins, *measured*), and `r26.11.2` (`57,780,176` distinct proteins, *measured*) shows that **card B-54's `biobox.l2.cluster` cost model (`CALL_FIT` + `CALL_T0` + `kmermatcher_cpu_min`) already predicts every single judged clustering step and part of `r26.11.2` within ±13% (`0.870` to `1.080`)** — no change to `src/biobox/l2/cluster.py` is needed, and no published cluster plan hash is touched.

Instead, all four `r26.11.2` clustering misses noted in `DEC-61` (`nr50.20.align` at `1.186`–`1.265`, `nr50.19.prefilter` at `0.629`–`0.643`, `nr50.15.prefilter` at `0.725`–`0.768`, `nr50.12.align` at `0.759`–`0.829`) plus three hidden `rescorediagonal` under-estimates (`0.300`–`0.457`) and `rM.2 -> r26.11.1` misses (`nr50.00.kmermatcher` at `0.608`, `nr50.19.prefilter` up to `1.252`) came from **four bugs in `scripts/stress/build_stage_pod/train_est.py`**:

1. **Pre-B-54 vs post-B-54 calibration mismatch**: `r26.11.1` was clustered *before* card B-54, so its `*.cut.json` files store legacy (`_LEGACY_CALL_FIT`) `part["cpu_min"]` estimates. `r26.11.2` was clustered *after* card B-54, so its `*.cut.json` files already store B-54 target-scaled `part["cpu_min"]` estimates. `train_est.py` divided `r26.11.1`'s measured totals by `r26.11.1`'s *stored legacy* `estimate_sum` (producing `0.7335` for `nr50.20.align`, `1.6805` for `nr50.19.prefilter`, `1.3161` for `nr50.15.prefilter`, `1.2288` for `nr50.12.align`) and multiplied `r26.11.2`'s *post-B-54* cut estimates by those legacy factors.
2. **Unjudged `rescorediagonal` 24-part -> 1-part factor inflation**: In `r26.11.1` (`parts = 24`), `nr50.01.rescorediagonal` (`15.06` CPU-min total B-54 estimate) and `nr90.01.rescorediagonal` (`18.31` CPU-min) paid 24× process-startup overhead (`44.22` and `50.85` CPU-min measured, factors `2.936` and `2.777`). In `r26.11.2` (`part_cpu_min = 240.0`, `k = 1` part), multiplying their `16.13` and `22.10` CPU-min cut estimates by `2.936` and `2.777` inflated them to `47.36` and `61.37` CPU-min (turning unjudged `< 24` CPU-min steps into judged failures at `0.365` and `0.300`).
3. **Single-release noise amplification inside the B-54 noise band**: Across all three releases (`rM.2`, `r26.11.1`, `r26.11.2`), every substantive call's raw `measured_sum / b54_sum` lies within `0.902` to `1.119` (within ±12% of `1.0`). Multiplying an already unbiased B-54 model by a single earlier release's residual doubles variance: on `rM.2 -> r26.11.1`, `rM.2`'s `nr50.19.prefilter` residual (`0.9387`) amplified `r26.11.1`'s worst part (`1.175` raw) to `1.252`. Keeping `factor = 1.0` inside a `±15%` deadband (`CALIBRATION_DEADBAND = 0.15`) eliminates noise amplification while still tracking any genuine systematic drift (`> 15%`) at larger scales.
4. **Pre-launch projection & `kmermatcher` thread scaling**: `train_est.py before` called `l2.call_cpu_min(name, n, u)` without `targets=...`, kept `k = len(rows_)` (`24` fixed parts) instead of computing `k` from `part_cpu_min = 240.0` (projecting `nr50.11.prefilter` at `1,711.3` CPU-min vs `232.9` measured), and scaled `kmermatcher` steps by protein ratio `r` without adjusting for `l2.kmermatcher_cpu_min` thread-contention overhead (`48t -> 18t` on `rM.2 -> r26.11.1` gave `0.608`).

## 2. Measured Clustering Calls Across `rM.2`, `r26.11.1`, and `r26.11.2`

All `targets`, `units`, `stored_sum`, and `meas_sum` values are *measured* from `/tmp/if2453/r10/root` on `wilkcomp02` (2026-11-21); `b54_sum` values are *derived* via `biobox.l2.cluster.call_cpu_min(call, proteins, units, targets=targets) + k * overhead_cpu_min`.

### 2.1 `rM.2` (`45,815,843` proteins, `parts = 24`, legacy pre-B-54 cuts)

| Call | k | Targets (*measured*) | Units (*measured*) | Stored sum (*measured*) | B-54 sum (*derived*) | Measured sum (*measured*) | `meas / B-54` (*derived*) | Part `meas / B-54` range (*derived*) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `cluster.nr90.07.align` | 24 | 33,085,999 | 1,356,134,142 | 365.38 | 338.23 | 340.23 | 1.0059 | 0.941–1.094 |
| `cluster.nr90.11.prefilter` | 24 | 31,920,514 | 10,194,781,853 | 3,382.92 | 3,629.20 | 3,598.11 | 0.9914 | 0.964–1.034 |
| `cluster.nr90.12.align` | 24 | 31,920,514 | 8,616,200,048 | 3,161.46 | 3,452.71 | 3,500.65 | 1.0139 | 0.983–1.036 |
| `cluster.nr90.members` | 24 | 31,183,876 | 432,415,353 | 132.67 | 144.29 | 148.46 | 1.0289 | 0.809–1.555 (unjudged, 6.0/part) |
| `cluster.nr50.07.align` | 24 | 17,348,804 | 698,116,275 | 309.01 | 269.14 | 279.40 | 1.0381 | 0.934–1.107 |
| `cluster.nr50.11.prefilter` | 24 | 13,520,716 | 4,070,572,296 | 1,255.18 | 1,283.19 | 1,307.65 | 1.0191 | 0.996–1.063 |
| `cluster.nr50.12.align` | 24 | 13,520,716 | 3,038,129,206 | 992.04 | 1,192.49 | 1,259.10 | 1.0559 | 0.979–1.104 |
| `cluster.nr50.15.prefilter` | 24 | 8,915,865 | 2,601,223,781 | 314.91 | 393.90 | 421.93 | 1.0712 | 1.040–1.112 |
| `cluster.nr50.16.align` | 24 | 8,915,865 | 2,717,279,156 | 785.74 | 688.33 | 707.66 | 1.0281 | 0.953–1.095 |
| `cluster.nr50.19.prefilter` | 24 | 8,312,546 | 2,454,500,276 | 723.70 | 1,020.53 | 957.96 | 0.9387 | 0.889–0.978 |
| `cluster.nr50.20.align` | 24 | 8,312,546 | 2,611,996,035 | 880.18 | 673.52 | 685.88 | 1.0183 | 0.951–1.056 |
| `cluster.nr50.members` | 24 | 8,192,065 | 409,423,542 | 361.00 | 322.82 | 332.32 | 1.0294 | 0.890–1.108 |

### 2.2 `r26.11.1` (`49,102,041` proteins, `parts = 24`, legacy pre-B-54 cuts)

| Call | k | Targets (*measured*) | Units (*measured*) | Stored sum (*measured*) | B-54 sum (*derived*) | Measured sum (*measured*) | Old factor `meas/stored` (*derived*) | Raw `meas/B-54` (*derived*) | New factor (*derived*) | Part `meas / B-54` range (*derived*) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `cluster.nr90.07.align` | 24 | 35,153,162 | 1,412,500,763 | 385.61 | 355.29 | 340.41 | 0.8828 | 0.9581 | 1.0000 | 0.907–1.021 |
| `cluster.nr90.11.prefilter` | 24 | 33,882,143 | 10,664,275,036 | 3,613.46 | 3,887.53 | 3,832.62 | 1.0607 | 0.9859 | 1.0000 | 0.963–1.024 |
| `cluster.nr90.12.align` | 24 | 33,882,143 | 9,125,396,776 | 3,348.29 | 3,663.30 | 3,656.60 | 1.0921 | 0.9982 | 1.0000 | 0.971–1.026 |
| `cluster.nr90.members` | 24 | 33,117,773 | 463,925,032 | 142.34 | 154.99 | 162.48 | 1.1415 | 1.0483 | 1.0000 | 0.808–1.580 (unjudged, 6.5/part) |
| `cluster.nr50.07.align` | 24 | 18,876,326 | 758,111,957 | 338.84 | 294.00 | 289.75 | 0.8551 | 0.9855 | 1.0000 | 0.892–1.060 |
| `cluster.nr50.11.prefilter` | 24 | 14,912,484 | 4,385,758,277 | 1,385.10 | 1,421.43 | 1,433.59 | 1.0350 | 1.0086 | 1.0000 | 0.979–1.043 |
| `cluster.nr50.12.align` | 24 | 14,912,484 | 3,333,344,123 | 1,088.44 | 1,317.37 | 1,337.46 | 1.2288 | 1.0153 | 1.0000 | 0.937–1.053 |
| `cluster.nr50.15.prefilter` | 24 | 10,064,063 | 2,852,219,962 | 354.32 | 465.24 | 466.31 | 1.3161 | 1.0023 | 1.0000 | 0.966–1.076 |
| `cluster.nr50.16.align` | 24 | 10,064,063 | 3,073,734,693 | 896.24 | 784.31 | 755.21 | 0.8426 | 0.9629 | 1.0000 | 0.905–1.021 |
| `cluster.nr50.19.prefilter` | 24 | 9,421,654 | 2,698,588,280 | 817.52 | 1,228.31 | 1,373.84 | 1.6805 | 1.1185 | 1.0000 | 1.078–1.175 |
| `cluster.nr50.20.align` | 24 | 9,421,654 | 2,964,472,779 | 1,013.59 | 776.96 | 743.50 | 0.7335 | 0.9569 | 1.0000 | 0.903–1.007 |
| `cluster.nr50.members` | 24 | 9,291,312 | 440,098,571 | 390.47 | 348.32 | 356.51 | 0.9130 | 1.0235 | 1.0000 | 0.919–1.094 |

### 2.3 `r26.11.2` (`57,780,176` proteins, `part_cpu_min = 240.0`, post-B-54 cuts): Old vs New `train_est.py`

| Call | k | Targets (*measured*) | Units (*measured*) | Cut/B-54 sum (*measured*) | Measured sum (*measured*) | Old `train_est` part ratios (*derived*) | New `train_est` part ratios (*derived*) |
|---|---:|---:|---:|---:|---:|---:|---:|
| `cluster.nr90.07.align` | 2 | 37,975,640 | 1,519,158,769 | 386.27 | 351.67 | 1.009–1.054 | **0.891–0.930** |
| `cluster.nr90.11.prefilter` | 11 | 36,380,963 | 11,254,104,592 | 2,559.36 | 2,314.89 | 0.830–0.899 | **0.880–0.953** |
| `cluster.nr90.12.align` | 17 | 36,380,963 | 9,811,385,661 | 3,947.10 | 3,977.82 | 0.903–0.935 | **0.986–1.022** |
| `cluster.nr90.members` | 1 | 35,475,280 | 544,385,754 | 182.12 | 180.48 | 0.868 | **0.991** |
| `cluster.nr50.07.align` | 2 | 20,249,969 | 821,248,164 | 320.05 | 311.92 | 1.109–1.171 | **0.948–1.001** |
| `cluster.nr50.11.prefilter` | 1 | 16,085,986 | 4,633,835,937 | 237.97 | 232.88 | 0.946 | **0.979** |
| `cluster.nr50.12.align` | 6 | 16,085,986 | 3,597,818,614 | 1,429.45 | 1,409.13 | 0.759–0.829 (OUT) | **0.933–1.018** |
| `cluster.nr50.15.prefilter` | 3 | 10,927,842 | 3,029,785,935 | 501.84 | 495.16 | 0.725–0.768 (OUT) | **0.954–1.010** |
| `cluster.nr50.16.align` | 4 | 10,927,842 | 3,343,819,070 | 857.45 | 796.61 | 1.063–1.130 | **0.896–0.952** |
| `cluster.nr50.19.prefilter` | 6 | 10,235,527 | 2,868,261,896 | 1,371.69 | 1,463.42 | 0.629–0.643 (OUT) | **1.057–1.080** |
| `cluster.nr50.20.align` | 4 | 10,235,527 | 3,225,662,877 | 854.57 | 770.78 | 1.186–1.265 (OUT, DEC-61) | **0.870–0.928** |
| `cluster.nr50.members` | 2 | 10,091,519 | 519,001,993 | 411.79 | 405.57 | 1.036–1.121 | **0.946–1.024** |
| `cluster.nr50.01.rescorediagonal` | 1 | 35,475,280 | 5,488,078,645 | 16.13 | 17.27 | 0.365 (OUT, est 47.36) | **unjudged (est 16.13, 1.071)** |
| `cluster.nr50.06.rescorediagonal` | 1 | 20,249,969 | 2,078,959,494 | 8.73 | 8.69 | unjudged (est 17.46) | **unjudged (est 8.73, 0.996)** |
| `cluster.nr90.01.rescorediagonal` | 1 | 57,780,176 | 4,040,497,318 | 22.10 | 18.42 | 0.300 (OUT, est 61.37) | **unjudged (est 22.10, 0.834)** |
| `cluster.nr90.06.rescorediagonal` | 1 | 37,975,640 | 1,892,665,885 | 15.18 | 11.17 | 0.457 (OUT, est 24.43) | **unjudged (est 15.18, 0.736)** |

## 3. Non-Split Clustering Steps (`kmermatcher`, `.cut`, `cluster`, `cluster.plan`)

| Step | `rM.2` (*measured*) | `r26.11.1` (*measured*) | `r26.11.2` (*measured*) | `rM.2 -> r26.11.1` new ratio (*derived*) | `r26.11.1 -> r26.11.2` new ratio (*derived*) |
|---|---:|---:|---:|---:|---:|
| `cluster.plan` | 11.41 | 10.92 | 13.04 | 0.893 (unjudged) | 1.015 (unjudged) |
| `cluster.nr90.00.kmermatcher` | 28.78 (32t) | 30.52 (18t) | 35.56 (18t) | **0.989** | **0.990** |
| `cluster.nr50.00.kmermatcher` | 29.30 (48t) | 19.09 (18t) | 20.35 (18t) | **0.948** (unjudged; was 0.608) | **0.906** (unjudged) |
| `cluster.nr50.15.prefilter.cut` | 32.21 (48t) | 34.46 (18t) | 36.45 (18t) | **0.998** | **0.899** |
| `cluster.nr50.19.prefilter.cut` | 39.47 (48t) | 42.16 (18t) | 44.63 (18t) | **0.997** | **0.900** |
| `cluster` (finish) | 80.44 (48t) | 75.08 (18t) | 94.65 (18t) | **0.871** | **1.071** |

`biobox.l2.cluster.kmermatcher_cpu_min` independently matches all six `kmermatcher` records within ±5% (`0.9562` to `1.0470`, *derived*).

## 4. Verification

All judged clustering steps and parts (`>= 24` CPU-min, `DEC-43`) across both `r26.11.1` (recomputed under B-54, calibrated from `rM.2`) and `r26.11.2` (calibrated from legacy `r26.11.1`) lie within 20% (`0.80 <= ratio <= 1.20`), with `outside_20pct == []`, `over_20pct == []`, `recorded_over_20pct == []`, and `under_20pct == []`.

Live read-only verification on `wilkcomp02` (2026-11-21) running `train_est.py before --earlier r26.11.1 --release r26.11.2 --cpus 18` followed by `python3 -I train_est.py rates --root /tmp/if2453/r10/root --release r26.11.2` across all `122` published `r26.11.2` run records (*measured*):
- `records = 122`, `judged = 65`, `cpu_min = 13,235.8`, `estimate_of_recorded = 13,612.08` (`-2.76%`), `worst_over = 1.080` (`cluster.nr50.19.prefilter.part-00002`, `246.92` vs `228.615`), `outside_20pct = []`, `over_20pct = []`, `recorded_over_20pct = []`, `under_20pct = []`, exit `0`.

