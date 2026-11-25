# Hybrid `linclust NR90 + cascade NR50` Clustering & `r26.11.4`/`r26.11.5+` Scaling Brief (`DEC-75` Follow-Up)

- **Date:** 2026-11-24
- **Decisions informed:** `DEC-24`, `DEC-29`, `DEC-38`, `DEC-53`, `DEC-54` (Q9), `DEC-72`, `DEC-73`, `DEC-74`, `DEC-75` (3)
- **Measured sources:** [`2026-11-16_release_r26.11.1.md`](../measurements/2026-11-16_release_r26.11.1.md), [`2026-11-19_cluster_survivor_estimates.md`](../measurements/2026-11-19_cluster_survivor_estimates.md), [`2026-11-21_release_r26.11.2.md`](../measurements/2026-11-21_release_r26.11.2.md), [`2026-11-23_b68_cut_index_estimates.md`](../measurements/2026-11-23_b68_cut_index_estimates.md), [`2026-11-24_a26_r26_11_4_confs.md`](../measurements/2026-11-24_a26_r26_11_4_confs.md), [`2026-11-24_c38c_linclust_vs_cascade.md`](../measurements/2026-11-24_c38c_linclust_vs_cascade.md)
- **Code audited:** `src/biobox/l2/cluster.py`, `src/biobox/stage.py`, `src/biobox/registry/schemas.py`, `profiles/wilkcomp02.yaml`, `scripts/stress/build_stage_pod/train_est.py`, `release_train.sh`, `r50_watch.py`

Every number below carries a label per `AGENTS.md` §6: *measured* (host, date, record), *derived* (formula and inputs), or *estimate*.

---

## 1. Executive Summary & Headline Forecast

**Current full cascade (`src/biobox/l2/cluster.py` as-is) will NOT complete `r26.11.4` (`32` stages, `287.03 Gbp`, `250–300M` distinct proteins *estimate*, central `265M`) on `wilkcomp02`.** It hits **four independent hard walls** before finishing:

1. **Hard code refusal (`cluster_cpu_min_max: 20000`):** `cluster.nr90.11.prefilter` (`--diag-score 0`, unsplittable above `68.30M` targets under `DEC-29`) is projected at **`25,681` CPU-min** (`23,200–31,000` across `250–300M`, *derived*), exceeding `profiles/wilkcomp02.yaml`'s `20,000` cap -> refused by `_one_process_doc` (`l2/cluster.py:1539`) and `one_process_refusal` (`stage.py:1014`).
2. **Physical `/tmp` exhaustion (`~450 GiB` free vs `~666 GiB` needed):** Unpruned full-cascade scratch scales at `2,696 B/prot` (`r26.11.3`: `275.9 GB` at `102.34M`, *measured*) -> **`715 GB` (`666 GiB`)** at `265M` (`629–754 GiB` across `250–300M`, *derived*), tripping `release_train.sh`'s `SCRATCH_MAX_GB=250` guard and physically exceeding `/tmp`'s `~450 GiB` free space.
3. **Process cap kill (`CPU_MAX=280` / `AGENTS.md` §7.2 `300` CPU-min) on `createindex`:** `cluster.nr50.15.prefilter.cut` (`-s 2`) and `cluster.nr50.19.prefilter.cut` (`-s 3`) scale super-linearly (`targets^0.48`, card `B-68`, *measured*) -> **`560` CPU-min** and **`685` CPU-min** at `265M` (*derived*), with no `DEC-29` exemption in `r50_watch.py`.
4. **Process cap kill (`CPU_MAX=280`) on `cluster` (`finish`):** Single-process table export scales at `1.64–1.85` CPU-min/M (*measured* `r26.11.2`: `94.65` CPU-min at `57.78M`) -> **`435–490` CPU-min** at `265M` (`410–555` across `250–300M`, *derived*), with no `DEC-29` exemption.

**Recommended path (zero NR50 neighbour-fragmentation loss for `r26.11.4`, clear runway through `r26.11.5+`):**
- **Adopt Hybrid A (`linclust NR90 + 3-step cascade NR50`) — or Hybrid B1 (`linclust NR90 + 2-step -s 1,2 cascade NR50`) — PLUS dead-intermediate scratch pruning + `createindex`/`finish` budget fixes before cutting `r26.11.4`.**
  - **Hybrid A** eliminates Wall #1 (`nr90.11.prefilter` `25,681` CPU-min -> `0`), cuts total `r26.11.4` clustering CPU by **`35.7%`** (`119,400 -> 76,830` CPU-min = `16.0` h ideal / `~19–21` h effective at `80` night cores, *derived*), lowers single-process peak RAM from `322 GB` to `204 GB` (`ru_maxrss`), and with dead-intermediate pruning fits `/tmp` scratch at **`265–318 GB` (`247–296 GiB`)** (*estimate*) inside `~450 GiB` free, while keeping **100% of cascade NR50 sensitivity** (`+2.31%` NR90 clusters, `0.00%` discordant NR90 joins, `C-38c` *measured*).
  - **Hybrid B1** (dropping only NR50's final `-s 3` step `nr50.18..21`) cuts `r26.11.4` clustering CPU by **`69.4%`** (`119,400 -> 36,570` CPU-min = **`7.6` h ideal / `~9–11` h effective** at `80` night cores, *derived* — one overnight run!) and eliminates the `685` CPU-min `-s 3` `createindex` step, at a measured cluster-count delta of only **`+1.40%` (`r26.11.1`) / `+1.43%` (`r26.11.2`)** (*measured*) versus `+60.84%` for pure linclust (`C-38c`).
  - **Hybrid B2 (full two-level `linclust`)** (`~7,670–8,085` CPU-min at `265M`, `~10,100–10,700` at `350M`, `~14,500–15,800` at `500M`, *derived*) should wait until **`400–500M+` proteins (`r26.11.6+`)** when even 2-step NR50 cascade outgrows `/tmp` scratch or overnight wall budgets, because `C-38c` measured `18.35×` mean fragmentation of `size >= 100` families and `20–42%` top-cluster locus dilution on divergent partners.

---

## 2. Measured Baseline (`r26.11.1`, `r26.11.2`, `r26.11.3`, and `C-38c`)

| Corpus / Run | Stages | Stored Gbp | Distinct Proteins | NR90 Clusters (`%` prot) | NR50 Clusters (`%` prot) | Clustering CPU-min | Peak Single-Proc RSS / Unit MemPeak | Peak Scratch | Label |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `r26.11.1` (cascade, 24 parts) | 9 | `59.88` | `49,102,041` | `33,117,773` (`67.45%`) | `9,291,312` (`18.92%`) | `15,255.59` (`14,644.04` core) | `58.61 GB` (`nr90.11`) / `172.8 GB` | `158.73 GB` (`3,233 B/prot`) | *measured* (`wilkcomp02`, 2026-11-18) |
| `r26.11.1` (`C-38c` linclust-only) | 9 | `59.88` | `49,102,041` | `33,882,143` (`+2.31%`) | `14,943,715` (`+60.84%`) | `767.74` core (`409.65` NR90 + `352.55` NR50) | `37.77 GB` `ru_maxrss` (`769.1 B/prot`) / `26.71 GB` cgroup (`544.0 B/prot`) | — | *measured* (`wilkcomp02`, 2026-11-24) |
| `r26.11.2` (cascade, `B-54` dynamic cut) | 11 | `74.16` | `57,780,176` | `35,475,280` (`61.40%`) | `10,091,519` (`17.47%`) | `13,168.1` (`13,235.80` all steps) | `61.42 GB` (`nr90.11`) / `182.83 GB` | `181.6 GB` (`3,143 B/prot`) | *measured* (`wilkcomp02`, 2026-11-21) |
| `r26.11.3` (cascade, `DEC-29` on `nr90.11`) | 13 | `119.89` | `102,337,875` | `71.21M` (`69.58%`) | `22.16M` (`21.65%`) | `25,451.77` | `124.41 GB` (`nr90.11`) / `170.74 GiB` (`183.33 GB`) | `275.9 GB` (`2,696 B/prot`) | *measured* (`wilkcomp02`, 2026-11-23/24) |

Key measured scaling rates used in all derivations below:
- **NR90 linclust CPU (`C-38c`):** `409.65 / 49.102041M = 8.343` CPU-min/M proteins (*derived* from *measured* `C-38c`; `5.11%` of `r26.11.1` NR90 cascade `8,014.42` CPU-min).
- **Both-level linclust core CPU (`C-38c`):** `767.74 / 49.102041M = 15.636` CPU-min/M proteins (*derived*; `19.1×` faster than `r26.11.1` cascade).
- **Linclust peak RAM (`C-38c`, `nr90.00.kmermatcher`):** `769.1 B/prot` `ru_maxrss` (`544.0 B/prot` cgroup `MemoryPeak`, *measured*; matches `r26.11.1` `771.7 B/prot` and `r26.11.2` `764.4 B/prot`, *derived*).
- **Cascade `nr90.11.prefilter` peak RAM (`r26.11.3`):** `124.41 GB / 102.338M = 1,216 B/prot` (`~1,709 B/target`, *derived* from *measured* `DEC-73`).
- **Shared-index `createindex` CPU (`B-68`, `train_est.py:61`):** `rate(s) * 1e-9 * residues * (targets / 10M)^0.48`, with `rate("2") = 11.80`, `rate("3") = 15.75` (`±2.7%` across all 5 measured records; `r26.11.3` `-s 2` measured `113.827` CPU-min and `55.98 GB` RSS = `2,304 B/target` on `24,298,242` targets / `6.476B` residues, *measured*).
- **Cascade `nr50` `-s 3` marginal yield (`r26.11.1` / `r26.11.2`, *measured*):** `r26.11.1` `input_step2` `9,421,654` -> final `9,291,312` (`+1.40%` clusters without `-s 3`); `r26.11.2` `10,235,527` -> `10,091,519` (`+1.43%` clusters without `-s 3`).
- **`cluster` (`finish`) CPU (`r26.11.1` / `r26.11.2`, *measured*):** `75.08` CPU-min at `49.10M` (`1.529` CPU-min/M); `94.65` CPU-min at `57.78M` (`1.638` CPU-min/M).
- **Cascade `CALL_FIT` (`B-54`, `src/biobox/l2/cluster.py:645`):** `work = r * 1e-9 * units * (targets / CALL_T0)^g`, with `nr90.11.prefilter = (30.24, 0.83)`, `nr50.11.prefilter = (16.83, 0.63)`, `nr50.15.prefilter = (75.14, 0.63)`, `nr50.19.prefilter = (194.50, 0.73)`.
- **Whole-release cascade power law `N17` (`stage.py:950`):** `284.10 * (n / 2,053,075)^1.244` CPU-min (*derived*).

---

## 3. Corpus Size Projection for `r26.11.4` (`32` Stages) and `r26.11.5+`

From [`2026-11-24_a26_r26_11_4_confs.md`](../measurements/2026-11-24_a26_11_4_confs.md) §3 (*measured*):
- `r26.11.3` (`13` built stages `S1..SPX1`): `119.89 Gbp` stored, `102.34M` distinct proteins (*measured*).
- `19` queued `r26.11.4` stages (`V2 -> PL2 -> PL3 -> CG1 -> MGE1 -> SAT1 -> ECE1 -> O1 -> O2 -> LPL1 -> SYM1 -> LRV1 -> LGRT1 -> MGE2 -> MV1 -> V3 -> LR1 -> GART1 -> GART2`): `10,134,020` net new contigs, **`167.14 Gbp`** (`167,137,613,539` bp) net new bases (*measured*).
- **`r26.11.4` total stored sequence:** `119.89 + 167.14` = **`287.03 Gbp`** (*derived*).
- **Distinct protein yield:** `r26.11.3`'s environmental/non-GTDB wave (`G1 + SPX1`) added `44.56M` distinct proteins over `45.73 Gbp` = **`0.974 M prot/Gbp`** (*derived*). Viral stages (`V1`: `1.044 M prot/Gbp`, *derived*) and `GART1`/`GART2` (`1.0–1.1 M prot/Gbp`, *estimate*) span `0.90–1.10 M new distinct prot/Gbp` -> `150.4M–183.9M` new distinct proteins -> **`253M–286M` total distinct proteins** (*derived*; bounded as **`250M` low / `265M` central / `300M` stress** *estimate*).

| Target Release | Stored Gbp | Distinct Proteins | Cascade NR90 Reps (`~72%`) | Linclust NR90 Reps (`+2.31%`) | Cascade NR50 Reps (`~24%`) | NR50 `input_step1` (`-s 2`) / `input_step2` (`-s 3`) Targets | Label |
|---|---:|---:|---:|---:|---:|---:|---|
| `r26.11.3` (actual) | `119.89` | `102.34M` | `71.21M` (`69.6%`) | `~72.85M` (*derived*) | `22.16M` (`21.7%`) | `24.30M` (*meas*) / `~22.8M` (*est*) | *measured* / *derived* |
| **`r26.11.4` low** | `287.03` | `250.0M` | `180.0M` | `184.2M` | `58.0M` | `64.0M` / `60.0M` | *estimate* / *derived* |
| **`r26.11.4` central** | **`287.03`** | **`265.0M`** | **`190.8M`** | **`195.2M`** | **`64.0M`** | **`70.0M` / `66.0M`** | *estimate* / *derived* |
| **`r26.11.4` stress** | `287.03` | `300.0M` | `216.0M` | `221.0M` | `74.0M` | `81.0M` / `76.0M` | *estimate* / *derived* |
| **`r26.11.5` (`LR2+NX1+SAG1+GV1+SPX2`)** | `~380–400` | **`350.0M`** | `252.0M` | `257.8M` | `85.0M` | `93.0M` / `88.0M` | *estimate* / *derived* |
| **`r26.11.6+` (`G2..G10+SPX3..`)** | `~550–600` | **`500.0M`** | `360.0M` | `368.3M` | `120.0M` | `132.0M` / `124.0M` | *estimate* / *derived* |

---

## 4. Architecture Comparison on `r26.11.4` (`265M` Central; `250M–300M` Range) and `r26.11.5+`

Wall-clock hours show **ideal** (`CPU-h / cores`) and **effective** (in parentheses, accounting for `J=4 × 18 = 72` night / `J=2 × 18 = 36` day concurrency and serial `cut`/`join`/`finish` steps, `~80%` efficiency). Peak RAM is compared against `wilkcomp02`'s **`750G`** (`750 GiB` = `805.3 GB`, `DEC-73`) scope; peak `/tmp` scratch is compared against **`~450 GiB` (`~483 GB`)** free disk.

### 4.1 `r26.11.4` (`265M` proteins central; `250M–300M` range in brackets)

| Architecture | Total CPU-min (*derived*) | CPU-hours (*derived*) | Wall @ `80` night cores | Wall @ `40` day cores | Peak RAM vs `750G` (*derived*) | Peak `/tmp` Scratch vs `450 GiB` (*derived*/*estimate*) | NR50 Quality vs Cascade | Hard Walls Hit (unpatched) |
|---|---:|---:|---:|---:|---:|---:|---|---|
| **(a) Current Full Cascade** (`l2/cluster.py` as-is) | **`119,400`** `[108k–144k]` | `1,990` `[1,800–2,400]` | `24.9 h` (`~31 h` eff) | `49.8 h` (`~62 h` eff) | `322 GB` (`300 GiB`) `nr90.11` / `322 GB` `J=2` `nr50.19` (**OK**) | **`715 GB` (`666 GiB`)** `[629–754 GiB]` (**FAILS `/tmp`!**) | Baseline (`1.00×` NR50) | **4 walls:** `nr90.11` `25.7k > 20k` refusal; `/tmp` `ENOSPC`; `createindex` `560/685 > 280`; `finish` `460 > 280` |
| **(b) Hybrid A** (`linclust NR90` + 3-step `-s 1,2,3` cascade NR50) | **`76,830`** `[69.8k–92.2k]` | `1,281` `[1,163–1,537]` | **`16.0 h`** (`~20 h` eff) | `32.0 h` (`~39 h` eff) | `204 GB` (`190 GiB`) `kmer` / `322 GB` `J=2` `nr50.19` (**OK**) | `477–504 GB` unpruned (**TIGHT**); **`265–318 GB` (`247–296 GiB`) pruned (OK)** | **Lossless NR50** (`+2.31%` NR90 only) | **2 walls** (eliminates `nr90.11` refusal & cuts scratch): needs pruning + `createindex` & `finish` cap fixes |
| **(c1) Hybrid B1** (`linclust NR90` + 2-step `-s 1,2` cascade NR50, drop `-s 3`) | **`36,570`** `[33.5k–43.5k]` | `610` `[558–725]` | **`7.6 h`** (`~10 h` eff) | `15.2 h` (`~19 h` eff) | `204 GB` (`190 GiB`) `kmer` / `322 GB` `J=2` `nr50.15` (**OK**) | **`230–280 GB` (`214–261 GiB`) pruned (OK)** | **`+1.4%` NR50 clusters** (`r26.11.1/2` *meas*) | Needs pruning + `-s 2` `createindex` (`560`) & `finish` (`460`) fixes (drops `685` `-s 3` `createindex`!) |
| **(c2) Hybrid B2** (full 2-level `linclust`, `C-38c` mode) | **`7,670–8,085`** `[7.2k–9.2k]` | `128–135` `[121–153]` | **`1.6–1.7 h`** (`~2.5 h` eff) | `3.2–3.4 h` (`~4.5 h` eff) | `204 GB` (`190 GiB`) `ru_maxrss` / `144 GB` cgroup (**OK**) | `320–370 GB` unpruned; **`210–250 GB` (`195–233 GiB`) pruned (OK)** | **`+60.84%` NR50**, `18.35×` split on `size>=100` (`C-38c` *meas*) | Needs `finish` (`460`) cap fix only |

*Derivation notes for `265M` central:*
- Bottom-up `CALL_FIT` sum for (a) Full Cascade (`119,439` CPU-min: NR90 linclust `2,211` + NR90 cascade `44,200` [`nr90.11.prefilter` `25,681` + `nr90.12.align` `18,500`] + NR50 linclust `1,903` + NR50 3-step cascade `67,600` [`11.pf` `3,991` + `12.aln` `8,800` + `15.cut` `560` + `15.pf` `9,890` + `16.aln` `5,000` + `19.cut` `685` + `19.pf` `32,670` + `20.aln` `6,000`] + overhead `3,525`) agrees within **`0.7%`** of `N17` (`284.10 × (265M / 2.053M)^1.244 = 120,288` CPU-min).
- Hybrid A removes NR90 cascade (`-44,200` CPU-min) and adds a conservative `+2.31%` (`+1,590` CPU-min) to NR50 for the `+2.31%` extra NR90 reps -> `76,830` CPU-min.
- Hybrid B1 additionally removes NR50 step `-s 3` (`19.cut` `685` + `19.pf` `32,670` + `20.aln` `6,000`, `× 1.023` = `-40,260` CPU-min) -> `36,570` CPU-min.
- Hybrid B2 scales `C-38c`'s `15.64` CPU-min/M (`4,145` linear to `4,560` `n log n`) plus `3,525` CPU-min for `plan` + `members` + `finish` -> `7,670–8,085` CPU-min.

### 4.2 `r26.11.5+` Scaling Horizon (`350M`, `500M`, and `1B` Proteins, *derived* / *extrapolation*)

| Scale | Architecture | Total CPU-min | Wall @ `80` cores (ideal / eff) | Peak RAM (`ru_maxrss` / `J=2`) vs `750G` | Pruned `/tmp` Scratch vs `450 GiB` | Feasibility on `wilkcomp02` |
|---|---|---:|---:|---:|---:|---|
| **`350M`** (`r26.11.5`) | (a) Full Cascade | `169,809` (`N17`) | `35.4 h` / `~44 h` | `426 GB` single (`nr90.11`) | `~945 GB` unpruned (**FAIL**) | **Infeasible** (`nr90.11` `~44k` CPU-min; `/tmp` overflow) |
| **`350M`** (`r26.11.5`) | (b) Hybrid A | `~113,000` | `23.5 h` / `~29 h` | `269 GB` `kmer` / `430 GB` `J=2` (**OK**) | `~350–420 GB` (**TIGHT**) | Marginal (`nr50.15/19.cut` `93M` targets hits `l2/cluster.py:1466` `90.96M` refusal; `createindex` `~890/1,090` CPU-min) |
| **`350M`** (`r26.11.5`) | (c1) Hybrid B1 (`-s 1,2`) | **`~52,000`** | **`10.8 h` / `~14 h`** | `269 GB` `kmer` / `430 GB` `J=2` (**OK**) | **`~310–370 GB` (OK)** | **Sweet spot for `r26.11.5`** (with `createindex` cap fix + pruning) |
| **`350M`** (`r26.11.5`) | (c2) Hybrid B2 (linclust) | `~10,100–10,700` | `2.1–2.2 h` / `~3.5 h` | `269 GB` (`251 GiB`) (**OK**) | `~280–330 GB` (**OK**) | Fast fallback (`+60.8%` NR50 fragmentation) |
| **`500M`** (`r26.11.6+`) | (b) Hybrid A | `~188,000` | `39.2 h` / `~49 h` | `385 GB` `kmer` / `580 GB` `J=2` | `~500–600 GB` (**FAIL `/tmp`**) | Exceeds `/tmp` & wall budget |
| **`500M`** (`r26.11.6+`) | (c1) Hybrid B1 (`-s 1,2`) | `~83,000` | `17.3 h` / `~22 h` | `385 GB` `kmer` (`J=1` for `nr50.15`) | `~440–520 GB` (**AT `/tmp` LIMIT**) | Borderline on `/tmp`; `nr50.11` `DEC-29` `~10k` CPU-min |
| **`500M`** (`r26.11.6+`) | (c2) Hybrid B2 (linclust) | **`~14,500–15,800`** | **`3.1 h` / `~4.5 h`** | **`385 GB` (`358 GiB`) (OK)** | **`~400–475 GB` (TIGHT)** | **Recommended transition point (`>= 450–500M`)**; needs `/tmp` relief or streaming `finish` |
| **`1.0B`** (ceiling) | (c2) Hybrid B2 (linclust) | `~29,000–31,600` (`15.6k–18.3k` core) | `~6.3 h` / `~9 h` | **`769 GB` (`716 GiB`) `ru_maxrss` / `544 GB` (`507 GiB`) cgroup** | `~800–950 GB` (**FAIL `/tmp`**) | At `750G` RAM ceiling (`C-38c` §4); requires **2-pass hash-partitioned `linclust` pre-collapse** (`95%`/`90%`) + external scratch or TACC |

*Note on chunked / two-pass `linclust` pre-collapse (`95%` -> `90%`):* Partitioning proteins deterministically by `protein_uid` hash prefix into `K` chunks (e.g. `K=4` at `1B` proteins -> `250M` per chunk, `~192 GB` peak RAM each), running `linclust 95%` or `90%` per chunk, and merging chunk representatives globally is **not needed below `~900M` proteins** (since single-pass `kmermatcher` fits `750G` up to `~950M` proteins) and should be deferred until `>750M` proteins because chunk boundaries slightly reduce first-pass collapse efficiency while adding orchestration complexity (though hash-prefix partitioning still satisfies `DEC-24` order-independence).

---

## 5. Imminent `r26.11.4` Guard & Resource Risk Audit

| Risk / Guard | Where Enforced | Threshold | `r26.11.4` (`265M`) Value (*derived*) | Status under Full Cascade | Status under Hybrid A / B1 (patched) |
|---|---|---|---|---|---|
| `nr90.11.prefilter` `DEC-29` cap | `l2/cluster.py:1539`, `stage.py:1014`, `wilkcomp02.yaml:40` | `cluster_cpu_min_max = 20,000` | `25,681` CPU-min `[23.2k–31.0k]` | **HARD REFUSAL** | **ELIMINATED** (`0` CPU-min) |
| `/tmp` clustering scratch | `release_train.sh:63` (`SCRATCH_MAX_GB=250`, `TMP_MIN_GB=150`) + `~450 GiB` physical free | `250 GB` guard / `~483 GB` physical | `715 GB` (cascade) / `265–318 GB` (Hybrid A pruned) / `230–280 GB` (B1 pruned) | **PHYSICAL `ENOSPC`** | **PASS** with dead-intermediate pruning + `SCRATCH_MAX_GB=350` |
| `createindex` (`nr50.15.cut`, `nr50.19.cut`) | `r50_watch.py` (`CPU_MAX=280`), `AGENTS.md` §7.2 (`300`), `l2/cluster.py:1466` (`240` linear = `90.96M` targets) | `280` CPU-min kill (no `DEC-29` exemption) | `560` (`-s 2`, `70M` targets) & `685` (`-s 3`, `66M` targets) CPU-min | **KILLED AT 280m** | **PASS** after extending `DEC-29` one-process exemption to shared-index `createindex` (B1 drops `-s 3` altogether) |
| `cluster` (`finish`) single process | `r50_watch.py` (`CPU_MAX=280`), `AGENTS.md` §7.2 (`300`) | `280` CPU-min kill | `435–490` CPU-min `[410–555]` | **KILLED AT 280m** | **PASS** after splitting `finish` by level (`nr90` + `nr50`, each `~220–245` CPU-min) or exempting under `cluster_cpu_min_max` |
| `nr50.11.prefilter` `DEC-29` | `l2/cluster.py:1539`, `r50_watch.py` (`A-22`) | `20,000` cap (triggers `>68.3M` targets) | `3,991` CPU-min (`~100M` targets) | **PASS** (`A-22` handles `DEC-29`) | **PASS** (`3,991 << 20,000`, `~50` min wall @ 80 cores) |
| `kmermatcher` (`nr90.00`, `nr50.00`) | `r50_watch.py` (`CPU_MAX=280`) | `280` CPU-min | `159` (`nr90`) / `117` (`nr50`) CPU-min | **PASS** (hits `280` at `~450M` prot) | **PASS** |
| Peak RAM vs `750G` scope | `profiles/wilkcomp02.yaml:8` (`MemoryMax=750G`, `MemorySwapMax=0`) | `750 GiB` (`805.3 GB`) | `322 GB` (cascade) / `204–322 GB` (Hybrid A/B1) | **PASS** | **PASS** (`>400 GiB` margin) |
| `train_est.py` rate guard | `train_est.py:456` (`DEC-64`: stop `>1.50`) | `ratio <= 1.50` | `check_rates` uses live `.cut.json` (`B-54`/`B-68`); `before` scales `r26.11.3` | **PASS** (within `1.50`) | **PASS** once `train_est.py before` filters `calib["calls"]` to active `l2.calls(mode)` |

---

## 6. Exact Code Specification for Hybrid A (and Optional Hybrid B1)

### 6.1 `src/biobox/l2/cluster.py`
1. **Clustering mode parameter (`CLUSTERING_MODES = ("cascade", "hybrid-a", "hybrid-b1", "linclust")`):**
   - `workflow(level: str, mode: str = "cascade")`:
     - When `level == "nr90"` and `mode in ("hybrid-a", "hybrid-b1", "linclust")`: return only the 10 linclust calls (`00.kmermatcher` through `09.mergeclusters`), with `09.mergeclusters` writing directly to `clu90` (omitting `10.createsubdb`, `11.prefilter`, `12.align`, `13.clust`, `14.mergeclusters`).
     - When `level == "nr50"`: `"cascade"` and `"hybrid-a"` keep all 23 calls (`steps = [s1, s2, s3]`); `"hybrid-b1"` uses `steps = [s1, s2]` (19 calls, dropping `18..21`); `"linclust"` uses the 10 linclust calls only.
   - `calls(mode: str = "cascade")`: passes `mode` to `workflow()`, preserving exact step names (`cluster.nr90.00..09`, `cluster.rep90`, `cluster.nr50.00..22`, `cluster.merge`, `cluster.nr90.members`, `cluster.nr50.members`) so `train_est.py` calibration factors map 1-to-1 with no renamed steps!
   - `run(..., mode: str = "cascade")`: calls `mm("linclust", "db", "clu90", "tmp90", ...)` instead of `mm("cluster", ...)` when `mode != "cascade"`, keeping one-process `run()` and split `run_split()` **byte-identical** (`AGENTS.md` §4.3, card `B-41`).
2. **`PLAN_FORMAT` and `B-49` cluster reuse:**
   - Bump `PLAN_FORMAT = 2` and add `"method": mode` to `_INPUT_KEYS` (lines 576, 1214).
   - Retain read-compatibility for `PLAN_FORMAT == 1` (treating absent `"method"` as `"cascade"`) in `_ctx`, `read_record`, and `_source_problem`, so published `rM.2`/`r26.11.1`/`r26.11.2`/`r26.11.3` plans stay valid and verifiable while never cross-reusing between `"cascade"` and `"hybrid-a"`.
3. **Determinism & byte-identity (`AGENTS.md` §4.3):**
   - `mmseqs linclust` with `--shuffle 0`, sorted `protein_uid` FASTA, explicit `--threads`, and `B-41` key-range cut on `nr90.01/06.rescorediagonal` and `nr90.07.align` is deterministic and thread-count-independent. Pin a fixture test verifying `run(threads=1) == run(threads=4) == run_split(parts=2, threads=2)` byte-for-byte.
4. **`cluster_member_aln` schema impact:**
   - **Zero schema change.** `cluster_member`, `cluster_member_aln`, and `member_cluster` schemas in `src/biobox/registry/schemas.py` (and `tests/registry/schema_lock.json`) are completely untouched; `cluster.nr90.members` and `cluster.nr50.members` still run `mmseqs align db db clu aln_<level> -a -e inf --min-seq-id 0 -c 0 --alignment-mode 3` with `B-70` empty-alignment normalization.
5. **`cluster_uid` preservation (PI Decision Q3 below):**
   - `cluster_uid = blake3(sorted(member_uids) + params_hash(level))`.
   - **Recommended:** Keep `params_hash("nr90")` semantic (unchanged: threshold/coverage/mode/tool) and record `"method": "hybrid-a"` in `cluster_plan` and the clustering ledger (`tool_versions_json` / `reason_code`), so **~97.7% of unchanged NR90 clusters preserve their exact `cluster_uid`** across the `r26.11.3 -> r26.11.4` boundary and `c35_check.py`'s `100%` unchanged-member-set `cluster_uid` invariant holds. (Alternative: include `"method"` in `params_hash("nr90")`, which rotates 100% of NR90 `cluster_uid`s once while keeping NR50 `cluster_uid`s unchanged.)
6. **Dead-intermediate scratch pruning (`clean_consumed`):**
   - Add `clean_consumed(root, release, step_name, tmp, out=None)` called from `stage.py::_cluster_in_parts`:
     - After `cluster.rep90` publishes -> unlink `tmp90/w/`.
     - After `cluster.nr50.09.mergeclusters` publishes -> unlink `tmp50/w/linclust/`.
     - After `cluster.nr50.13.clust` / `17.clust` / `21.clust` publish -> unlink `tmp50/w/pref_step{k}`, `aln_step{k}`, and consumed `input_step*`.
   - In `_write_tables` (`finish`): process `nr90` end-to-end and delete `finish/*nr90*` temp files before `nr50` (halving `finish/` scratch).
   - Update `lost()` (line 1780) to skip pruned intermediate steps once their downstream milestone (`clu90` + `rep90`, `clu50`, etc.) is verified on disk — preserving `redo()` safety.
7. **`createindex` (`_cut_doc:1466`) & `finish` process-cap handling:**
   - In `_cut_doc`: compute `createindex` cost using `B-68`'s super-linear formula (moving `createindex_cpu_min` into `l2/cluster.py` so `cluster.py` and `train_est.py` share one source of truth). When `index_cpu_min > budget`, allow it under `one_process_cpu_min_max` (recording `"one_process": {"decision": "DEC-76", "cpu_min": ..., "cpu_min_max": cap}` on `<call>.cut.json`, refusing only above `cluster_cpu_min_max`).
   - Split `finish` into two `run-step` calls (`cluster.finish.nr90` + `cluster.finish.nr50`, plus lightweight `cluster` publish/ledger) OR exempt `finish` under `cluster_cpu_min_max` (see PI Decision Q4).

### 6.2 `src/biobox/stage.py`, `train_est.py`, `r50_watch.py`, `release_train.sh`
- **`src/biobox/stage.py`**: pass `mode` (`--cluster-method` / profile `cluster_method`) to `l2cluster.plan` / `calls` / `reuse`; call `l2cluster.clean_consumed` after milestone steps; pass `budget_reason="exempt: DEC-76"` on `createindex` `.cut` steps estimated above `part_cpu_min` and under `cluster_cpu_min_max`.
- **`scripts/stress/build_stage_pod/train_est.py`**: filter `calib["calls"]` in `before()` to active `l2.calls(mode)` so dropped `nr90.11..12` records from `r26.11.3` do not inflate `r26.11.4`'s projected CPU; import/share `createindex_cpu_min`; record `"cluster_method"` in `before.json`.
- **`scripts/stress/build_stage_pod/r50_watch.py`**: extend `_dec29_cap` (line 39) to recognize `one_process` records on `.cut` (`createindex`) steps and (if chosen) `finish`.
- **`scripts/stress/build_stage_pod/release_train.sh`**: raise default `SCRATCH_MAX_GB` from `250` to `350` (with `TMP_MIN_GB=100`) once dead-intermediate pruning is merged, and wire `CLUSTER_METHOD=${CLUSTER_METHOD:-hybrid-a}`.

---

## 7. PI Decision Table (5 Yes/No Questions to Unblock `r26.11.4`)

| # | Question | Recommendation | Impact if Yes | Impact if No |
|---|---|---|---|---|
| **Q1** | Switch `r26.11.4` from full cascade to **Hybrid A** (`linclust NR90 + 3-step cascade NR50`)? | **YES** (or choose Q2 for Hybrid B1) | Saves `~44,200` CPU-min (`-35.7%`, `~16–20 h` wall @ 80 cores), eliminates `nr90.11` `25.7k > 20k` refusal, lossless NR50 (`+2.31%` NR90). | Must raise `cluster_cpu_min_max` to `>= 35,000` AND solve `~715 GB` `/tmp` overflow (`~31 h` wall @ 80 cores). |
| **Q2** | Also drop NR50's final `-s 3` step (**Hybrid B1**: `linclust NR90 + 2-step -s 1,2 cascade NR50`) for `r26.11.4` & `r26.11.5`? | **YES** (optional fast-track; **NO** if strict 3-step NR50 parity is desired on `r26.11.4`) | Cuts `r26.11.4` clustering to `~36.6k` CPU-min (**`~7.6–10 h` wall @ 80 cores — one overnight cut!**), drops `685`-min `-s 3` `createindex`, at only **`+1.4%` NR50 clusters** (*measured* `r26.11.1/2`). | Keep 3-step NR50 (Hybrid A: `~16–20 h` wall @ 80 cores on `r26.11.4`; revisit at `r26.11.5` `350M`). |
| **Q3** | Keep `params_hash("nr90")` semantic (unchanged) so ~97.7% of unchanged NR90 clusters preserve their `cluster_uid`, recording `"method": "hybrid-a"` in `cluster_plan` (`PLAN_FORMAT=2`) and ledger? | **YES** | Preserves `100%` `cluster_uid` stability on unchanged member sets across `r26.11.3 -> r26.11.4` (`c35_check.py` passes cleanly). | Rotates 100% of NR90 `cluster_uid`s once on `r26.11.4` (NR50 `cluster_uid`s unaffected). |
| **Q4** | Extend `DEC-29`'s `cluster_cpu_min_max` one-process exemption to shared-index `createindex` (`.cut`, `560–685` CPU-min) AND either split `cluster finish` by level (`nr90`/`nr50`, `<250` CPU-min each) or exempt `finish` (`435–490` CPU-min)? | **YES** (exempt `createindex` under `cluster_cpu_min_max`; **split `finish` by level** so `finish` stays `<280` CPU-min) | Unblocks both `280` CPU-min process-cap walls cleanly with `<= 15` min wall per `createindex` and zero `finish` cap breach through `~400M` proteins. | `r26.11.4` will be killed by `r50_watch.py` at `280` CPU-min on `nr50.15.prefilter.cut` and `cluster`. |
| **Q5** | Enable **dead-intermediate scratch pruning** in `l2/cluster.py` (`clean_consumed` + level-sequential `finish`) and raise `release_train.sh` `SCRATCH_MAX_GB` from `250` to `350`? | **YES** | Brings `r26.11.4` peak `/tmp` scratch to **`230–318 GB` (`214–296 GiB`)**, fitting safely inside `~450 GiB` free on `wilkcomp02`. | `r26.11.4` hits `SCRATCH_MAX_GB=250` and physically exhausts `/tmp` (`ENOSPC`). |

---

## 8. Follow-Up Implementation Card Specs (Ready Once PI Answers Q1–Q5)

1. **Card `B-71` (Tier M): Hybrid clustering mode (`hybrid-a` / `hybrid-b1`), `PLAN_FORMAT = 2`, and `createindex` one-process exemption.**
   - *Files:* `src/biobox/l2/cluster.py`, `src/biobox/stage.py`, `src/biobox/cli.py`, `tests/r1/test_cluster_hybrid.py`, `tests/r2/test_cluster_split.py`.
   - *Acceptance:* Test-first (`AGENTS.md` §5.1); byte-identity `run(mode="hybrid-a", threads=1) == run(threads=4) == run_split(parts=2, threads=2)`; `PLAN_FORMAT = 1` read/reuse compatibility preserved; `createindex` `> budget` emits `one_process` under `cluster_cpu_min_max` and refuses above cap; `schema_lock.json` unchanged.
2. **Card `B-72` (Tier M): Dead-intermediate clustering scratch pruning (`clean_consumed`) and level-sequential/split `cluster finish`.**
   - *Files:* `src/biobox/l2/cluster.py`, `src/biobox/stage.py`, `tests/r1/test_cluster_pruning.py`.
   - *Acceptance:* Fixture test verifies byte-identical output tables (`cluster_member`, `cluster_member_aln`, `member_cluster`) with pruning enabled, `lost()`/`redo()` recovery after mid-run scratch loss, and `>= 40%` reduction in peak scratch directory bytes on the fixture; `finish` split keeps each level's `run-step` `< 280` CPU-min.
3. **Card `A-27` (Tier L): Release-train driver updates for hybrid clustering (`train_est.py`, `r50_watch.py`, `release_train.sh`).**
   - *Files:* `scripts/stress/build_stage_pod/train_est.py`, `scripts/stress/build_stage_pod/r50_watch.py`, `scripts/stress/build_stage_pod/release_train.sh`, `tests/r1/test_train_est.py`, `tests/r1/test_r50_watch.py`.
   - *Acceptance:* `train_est.py before` filters `calib["calls"]` to active `l2.calls(mode)` without inflating CPU when transitioning from a cascade `--earlier` release; `r50_watch.py` honors `one_process` on `.cut` (`createindex`) steps; `release_train.sh` sets `SCRATCH_MAX_GB=${SCRATCH_MAX_GB:-350}` and passes `--cluster-method`; rehearsal passes on pod.
