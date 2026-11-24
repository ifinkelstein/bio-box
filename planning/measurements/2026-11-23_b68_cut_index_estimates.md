# Card B-68: Shared-Index `.cut` (`mmseqs createindex`) Estimates in `train_est.py`

**Date**: 2026-11-23
**Host**: `wilkcomp02` (read-only inspection of `/tmp/if2453/r10/root`, `/tmp/if2453/c46/train_r26_11_3/estimate/before.json`, and `/tmp/if2453/c46/train/estimate/before.json`; zero writes to data root, train directories, or `/stor`)
**Scope**: Card `B-68` (`DEC-53`, `DEC-64`; depends on `B-43`, `B-54`, `B-64`, `C-46`)

## 1. Root Cause Summary

On live `r26.11.3` (`102,337,875` distinct proteins across 13 stages `S1..SPX1`, *measured* from `/tmp/if2453/r10/root` on `wilkcomp02`, 2026-11-23), `cluster.nr50.15.prefilter.cut` (`mmseqs createindex` `-s 2`, 18 threads) completed in `651.443 s` wall (`113.827` CPU-min, `55.98 GB` peak RSS, *measured*), well under the `280` CPU-min process cap, but `train_est.py rates` stopped the train at `ratio = 1.763` (`113.827 / 64.562`, *derived*) > `STOP_ABOVE = 1.50` (`ACCEPT_CEILING = 1.60`).

Two effects drove the `1.763×` miss:
1. `train_est.py` placed shared-index `.cut` steps (`cluster.nr50.15.prefilter.cut`, `cluster.nr50.19.prefilter.cut`) in `calib["others"]` and scaled `r26.11.2`'s record linearly by release-total distinct proteins `r = 102,337,875 / 57,780,176 = 1.771×` (`36.452 × 1.771 = 64.562` CPU-min, *derived*), whereas environmental metagenomes + archaea (`SPX1 + G1`) collapsed less at NR90/NR50 step 0 than strain-rich `C5`, so `input_step1` survivors jumped `2.224×` in target count (`24,298,242` vs `10,927,842`, *measured*) and `2.138×` in residues (`6,476,355,769` vs `3,029,785,935`, *measured*).
2. `mmseqs createindex` k-mer table construction scales super-linearly with target count: `12.08`, `12.03`, and `17.58` CPU-min per `10^9` residues at `10.06M`, `10.93M`, and `24.30M` targets for `-s 2` (*derived*: `measured_cpu_min / (residues * 1e-9)`), and `15.62` and `15.56` CPU-min per `10^9` residues for `-s 3` on `r26.11.1` and `r26.11.2` (*derived*).

Because `A-15` resumes keep `estimate/before.json` untouched, `train_est.py` `check_rates` (stdlib-only via `python3 -I` on the host) now prices shared-index `.cut` steps directly from `<root>/l2/cluster_plan/release=<R>/<call>.cut.json` (`parts` `units`/`keys` and `index.argv` `-s`), and `train_est.py before` projects shared-index `.cut` steps with the same model whenever `k > 1` (even when `k_old == 1` in the earlier release). `src/biobox/l2/cluster.py` was left untouched so published `cluster_plan` files and `B-49` reuse hashes never change.

## 2. Model and Fit Across All Five Measured 18-Thread `createindex` `.cut` Records

Formula (*derived*, `scripts/stress/build_stage_pod/train_est.py::createindex_cpu_min`):

```text
createindex_cpu_min(residues, targets, s) = rate(s) * 1e-9 * residues * (targets / 10_000_000) ** 0.48
```

with `CREATEINDEX_T0 = 10_000_000`, `CREATEINDEX_GROWTH = 0.48`, and `CREATEINDEX_RATE_BY_SENS = {"2": 11.80, "3": 15.75}` CPU-min per `10^9` residues at `10M` targets.

| Release | Step | `-s` | Targets (*measured*) | Residues (*measured*) | Measured CPU-min (*measured*) | `createindex_cpu_min` (*derived*) | Ratio `meas / est` (*derived*) |
|---|---|---|---:|---:|---:|---:|---:|
| `r26.11.1` | `cluster.nr50.15.prefilter.cut` | `2` | 10,064,063 | 2,852,219,962 | 34.457 | 33.759 | 1.021 |
| `r26.11.1` | `cluster.nr50.19.prefilter.cut` | `3` | 9,421,654 | 2,698,588,280 | 42.162 | 41.303 | 1.021 |
| `r26.11.2` | `cluster.nr50.15.prefilter.cut` | `2` | 10,927,842 | 3,029,785,935 | 36.452 | 37.307 | 0.977 |
| `r26.11.2` | `cluster.nr50.19.prefilter.cut` | `3` | 10,235,527 | 2,868,261,896 | 44.629 | 45.683 | 0.977 |
| `r26.11.3` | `cluster.nr50.15.prefilter.cut` | `2` | 24,298,242 | 6,476,355,769 | 113.827 | 117.028 | 0.973 |

All five measured 18-thread `createindex` `.cut` records sit within `±2.7%` (`0.973`–`1.021×`, *derived*), well inside the `±10%` acceptance band and `±20%` guard band. (Note: `r26.11.2` and `r26.11.3` `createindex_cpu_min` values above reflect the exact per-part sum of `keys` and `units` in the live `<call>.cut.json` files on `wilkcomp02`.)

## 3. Read-Only Live Verification on `wilkcomp02` (2026-11-23)

All commands ran read-only with `python3 -I` on `wilkcomp02` against `/tmp/if2453/r10/root`.

### 3.1 `r26.11.3` (109 published records, launch 1 `before.json` unchanged)

Command:
```sh
python3 -I /tmp/if2453/ci/B-68/scripts/stress/build_stage_pod/train_est.py rates \
  --root /tmp/if2453/r10/root --release r26.11.3 \
  --est /tmp/if2453/c46/train_r26_11_3/estimate/before.json
```

Result (*measured*, `wilkcomp02`, 2026-11-23):
- Exit code: `0`
- `records = 109`, `judged = 74`, `cpu_min = 15,419.35`, `estimate_of_recorded = 20,013.25`
- `over_20pct = []`, `recorded_over_20pct = []`
- `worst_over`: `cluster.nr50.00.kmermatcher` (`42.53` vs `36.05` CPU-min, `ratio = 1.18`)
- `cluster.nr50.15.prefilter.cut`: `113.827` CPU-min *measured* vs `117.028` CPU-min *derived* (`ratio = 0.973`, `wall_s = 651.443`, `peak_rss_gb = 55.98`, basis `cut createindex (-s 2, 24,298,242 targets)`)

### 3.2 `r26.11.2` (122 published records, launch 1 `before.json` unchanged)

Command:
```sh
python3 -I /tmp/if2453/ci/B-68/scripts/stress/build_stage_pod/train_est.py rates \
  --root /tmp/if2453/r10/root --release r26.11.2 \
  --est /tmp/if2453/c46/train/estimate/before.json
```

Result (*measured*, `wilkcomp02`, 2026-11-23):
- Exit code: `0`
- `records = 122`, `judged = 68`, `cpu_min = 13,235.8`, `estimate_of_recorded = 15,175.21`
- `over_20pct = []`
- `cluster.nr50.15.prefilter.cut`: `36.452` CPU-min *measured* vs `37.307` CPU-min *derived* (`ratio = 0.977`, `wall_s = 210.593`, `peak_rss_gb = 19.42`, basis `cut createindex (-s 2, 10,927,842 targets)`)
- `cluster.nr50.19.prefilter.cut`: `44.629` CPU-min *measured* vs `45.683` CPU-min *derived* (`ratio = 0.977`, `wall_s = 235.001`, `peak_rss_gb = 22.65`, basis `cut createindex (-s 3, 10,235,527 targets)`)
