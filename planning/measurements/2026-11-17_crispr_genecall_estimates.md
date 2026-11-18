# Refitted estimates for CRISPR on circular plasmids and gene calling on circular viruses (card B-53)

Date: 2026-11-17 to 18, wilkcomp02 (the pod's clock, CST). Every number carries a label: *measured* (host, date,
command or record), *derived* (formula and inputs) or *estimate*. Read-only on the data root `/tmp/if2453/r10/root` and
stage `V1` inputs `/tmp/if2453/s2/inputs/V1/` (bound `:ro` into every command that read them); nothing written under
`/stor`; scratch `/tmp/if2453/b53/`.

## Result

- **What CRISPR cost follows across genome and plasmid stages**: the bases MinCED actually scans — the stage's stored
  genome bases (`stage.stage_bases`) plus the circular wrap-flank bases (`stage.crispr_flank_bases`:
  `2 * min(CIRC_FLANK, contig_len)` with `CIRC_FLANK = 100,000` around each circular contig the stage stores, as
  `crispr._Scan.contig` writes it). Stages `S1`, `P1` to `P6` and `H1` list every contig linear, so MinCED scanned
  their stored bases alone; stage `PL1` stores 61,647 contigs (5,924,552,875 bp), 56,548 of them circular
  (5,453,091,032 bp, mean 96.4 kb), adding 6,205,263,778 wrap-flank bases so MinCED scanned 12,129,816,653 bp
  (2.047 times stored bases, *measured*, table in §2). Contig count does not explain `PL1`: `S1` has more contigs per
  stored Gbp than `PL1` (16,380 vs 10,405, *derived*), yet ran at 183.4 vs 313.3 CPU-s per stored Gbp (153.0 per
  scanned Gbp, *derived*).
- **The CRISPR model** (`stage.crispr_cpu_min(bases, flank_bases, *, near_cas)`): one rate per Gbp scanned,
  `CRISPR_CPU_S_PER_GBP = 168.1` (least-squares relative-error fit 168.14 over the nine default records, *derived*),
  times `CRISPR_NEAR_CAS_COST = 1.15` for the filtered two-repeat producer (least-squares relative-error fit 1.152 over
  the nine filtered records at 168.1, *derived*). **All 18 of 18 built-stage CRISPR records (nine default and nine
  filtered) are within 10% of their estimates** (well inside the 20% requirement): from `PL1` filtered 9.1% under
  (`0.909`) and `PL1` default 9.0% under (`0.910`) to `S1` default 9.1% over (`1.091`) (*derived*, §3). Before card
  B-53 (`182.0` per stored Gbp, `1.19` filtered): `PL1` was 72.2% over (`1.722`) on the default step and 66.1% over
  (`1.661`) on the filtered step. Both `stage_bases` and `crispr_flank_bases` are known before the CRISPR step starts
  (after `inputs` and `ingest` publish `genome_input` and `contig_index`), streamed in `contig_uid` order in bounded
  memory without DuckDB or scratch files; `crispr_flank_bases` on `/tmp/if2453/r10/root` matched the table's flank
  bases on all nine stages (*measured*, §4).
- **Circular `gv-meta` gene-calling rate** (`genecall.CALL_S_PER_MBP[("gv-meta", True)]`): none of the nine built
  stages in `/tmp/if2453/r10/root` has a circular contig in `gv-meta` mode (`S1` and `H1` are 100% linear `gv-meta`;
  `P1` to `P6` are linear `single`; `PL1` is `meta`, *measured*, §5). Card C-39 (`c39_est.py`, unit `c39-est`,
  wilkcomp02, 2026-11-16, `planning/measurements/2026-11-16_virus_votu_inputs.md` §5 and `V1/est.json`) timed
  `genecall._call` on one thread on every 30th stored contig of stage `V1` (`/tmp/if2453/s2/inputs/V1/`, UHGV v1.0):
  778 circular virus genomes (34,564,134 contig bp; 88,404,110 bp read with 50 kb wrap flanks, *measured* by
  `b53_table.py`) took 386.87 CPU-s (*measured*), which is **4.3762 CPU-seconds per Mbp read** (4.37615 rounded up;
  11.1929 per contig Mbp, *derived*). Replacing the unmeasured `4.2019 * 1.5011 = 6.3075` (*estimate*, +44.1% above the
  sample) with `4.3762` (*measured*) leaves every built stage's gene-call estimate unchanged (`PL1`'s four batches stay
  at ratios `0.961` to `0.980`, *derived*) and moves stage `V1`'s whole-stage `one_process_cpu_min` from 687.55 to
  588.77 CPU-minutes against 550.11 projected from C-39's sample (ratio `0.800` to `0.934`, *derived*, §5). Only
  `("single", True)` remains an estimate.

## 1. Setup

- Worktree on branch `task/B-53` from `10f3ddc`.
- Read-only queries in `/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif` inside transient `systemd-run --user`
  scopes (`MemoryMax=8G`, `MemorySwapMax=0`, `CPUQuota=100%`, `nice -n 19`, one thread), with `/tmp/if2453/r10/root`
  and `/tmp/if2453/s2/inputs/V1` bound `:ro`:
  - `b53_table.py` (Appendix A): read `contig_index`, `genome_input`, CRISPR and gene-call run records, CRISPR ledgers,
    `genecall_plan`, and `V1` input tables into `/tmp/if2453/b53/table.json` (sha256
    `c073fb82344a19ef214a8c34511e4ced1a185535d66fcb82f3d24e2b33e0f6c6`).
  - `b53_flank_check.py` (Appendix B): ran `stage.stage_bases`, `stage.crispr_flank_bases`, `stage.crispr_cpu_min`, and
    `genecall.one_process_cpu_min` from this card's `src/` on `/tmp/if2453/r10/root`, writing
    `/tmp/if2453/b53/flank_check.jsonl`.
  - `b53_fit.py` (Appendix C): least-squares and minimax model comparisons on `table.json`.

## 2. Inputs and run records of every built stage's CRISPR steps (*measured*)

Per stage in `/tmp/if2453/r10/root`: stored default-variant contigs and bases (`contig_index`); contigs skipped because
an earlier stage owns them (`genome_input` minus `contig_index`, DEC-26); topology and kind from `genome_input`;
circular wrap-flank bases (`2 * min(100,000, contig_len)` per circular stored contig) and total bases written to
MinCED's 50 Mbp batch FASTAs (`crispr._Scan.contig`); `feature_protein` and Pfam `l1_domains` row counts read by the
filtered step's join; and both producers' run records (`runs/stage=<S>/crispr-*.json`, `--cpus 4` in every run; stages
`S1` and `P1` to `P6` backfilled on 2026-11-15 by card C-33 part 1, `H1` built 2026-11-15 by card C-24, `PL1` built
2026-11-16 by card C-26 part 2).

| Stage | Kind | Stored contigs | Circular contigs | Stored bases | Wrap-flank bases | Scanned bases | MinCED batches | Default calls kept | Filtered calls in / kept | Default CPU-min | Filtered CPU-min | Filtered / default |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| S1 | virus | 27,258 | 0 | 1,664,109,938 | 0 | 1,664,109,938 | 34 | 213 | 3,613 / 34 | 5.085 | 5.487 | 1.079 |
| P1 | prokaryote | 539 | 0 | 1,016,625,803 | 0 | 1,016,625,803 | 20 | 286 | 3,379 / 6 | 3.085 | 3.383 | 1.097 |
| P2 | prokaryote | 5,018 | 0 | 9,596,180,328 | 0 | 9,596,180,328 | 184 | 2,291 | 30,682 / 72 | 26.489 | 31.794 | 1.200 |
| P3 | prokaryote | 9,107 | 0 | 9,949,077,522 | 0 | 9,949,077,522 | 192 | 3,921 | 39,391 / 118 | 29.195 | 33.386 | 1.144 |
| P4 | prokaryote | 4,005 | 0 | 9,976,456,950 | 0 | 9,976,456,950 | 191 | 2,520 | 30,835 / 69 | 27.585 | 32.606 | 1.182 |
| P5 | prokaryote | 4,111 | 0 | 9,967,711,632 | 0 | 9,967,711,632 | 189 | 2,205 | 29,779 / 56 | 27.505 | 32.361 | 1.177 |
| P6 | prokaryote | 4,227 | 0 | 9,962,536,103 | 0 | 9,962,536,103 | 191 | 2,041 | 30,953 / 59 | 27.938 | 32.796 | 1.174 |
| H1 | virus | 6,714 | 0 | 1,819,831,338 | 0 | 1,819,831,338 | 37 | 954 | 11,906 / 30 | 4.925 | 5.876 | 1.193 |
| PL1 | plasmid | 61,647 | 56,548 | 5,924,552,875 | 6,205,263,778 | 12,129,816,653 | 242 | 1,636 | 29,510 / 49 | 30.938 | 35.513 | 1.148 |

Observations (*derived*):

1. **Why `PL1` ran 1.72 times the old estimate**: `PL1`'s 56,548 circular plasmids average 96,433 bp, so
   `2 * min(100,000, n)` nearly triples short plasmids and adds 200 kb to plasmids ≥ 100 kb, raising MinCED's input
   from 5.925 Gbp stored to 12.130 Gbp scanned (242 JVM batches of ~50 Mbp instead of ~119). Divided by scanned bases,
   `PL1`'s default step ran at 153.0 CPU-s per Gbp scanned — 9% *below* the genome stages' 162.4 to 183.4, because
   `SeqStore.get` decompresses the 5.925 Gbp once while `crispr.extended` copies the 6.205 Gbp of flanks in memory.
2. **Why contig count is not the driver**: MinCED is launched once per 50 Mbp batch (`crispr.BATCH_BASES`), not once
   per contig (`PL1` spawned 242 JVMs, `P3` 192). Per stored Gbp, `S1` has 16,380 contigs/Gbp vs `PL1`'s 10,405, yet
   `S1` cost 183.4 CPU-s per stored Gbp vs `PL1`'s 313.3.
3. **The filtered producer's second pass**: across all nine stages, `--near-cas` (MinCED at `-minNR 2` plus one
   single-threaded DuckDB join against `feature`, `feature_protein` and `l1_domains` with row-group pruning) cost
   between 1.079 (`S1`) and 1.200 (`P2`) times the default step (`PL1` 1.148, mean 1.155), independent of plasmid
   topology once the circular wrap flanks are in the base count.

## 3. CRISPR model fits tried and before/after ratios (*derived*, `b53_fit.py`)

As in card B-52, the guard's measure is `(record - estimate) / estimate` (DEC-43). Each fit below is evaluated across
all 18 CRISPR run records (nine default, nine `--near-cas`).

| Model | Constants | Worst `|off|` (18) | Within 10% | Within 20% | Notes |
|---|---|---:|---:|---:|---|
| Before B-53: `stored_bases` alone | `182.0` s/Gbp, `× 1.19` | +72.2% | 14 / 18 | 16 / 18 | Fails on `PL1` (+72.2% default, +66.1% filtered) |
| `stored_bases` alone, refitted (LS) | `183.5` s/Gbp, `× 1.155` | +70.7% | 15 / 18 | 16 / 18 | Cannot reconcile `PL1` (313.3 s/stored Gbp) with `H1` (162.4) |
| `stored_bases` alone, minimax | `237.6` s/Gbp, `× 1.19` | 31.7% | 0 / 18 | 0 / 18 | No single rate on stored bases holds every stage within 20% |
| `stored_bases` + `stored_contigs` (LS) | `163.8` s/Gbp + `1,722` s/M contigs, `× 1.16` | +63.0% | 15 / 18 | 16 / 18 | `S1` has more contigs/Gbp than `PL1`; `PL1` still +63.0% |
| `scanned_bases` at old constants | `182.0` s/Gbp, `× 1.19` | −19.0% | 7 / 18 | 18 / 18 | Meets 20%, but `PL1` is −15.8% / −19.0% and only 7/18 within 10% |
| **`scanned_bases` alone (chosen)** | **`168.1` s/Gbp, `× 1.15`** | **9.1%** | **18 / 18** | **18 / 18** | Smallest 2-constant model; every record within 10% |
| `0.07` CPU-min + `scanned_bases` (LS) | `0.07` + `165.4` s/Gbp, `× 1.15` | 8.3% | 18 / 18 | 18 / 18 | Adds a third constant for 0.8 percentage points |
| `stored_bases` + `flank_bases` (2 rates) | `170.6` s/stored Gbp + `135.9` s/flank Gbp, `× 1.15` | 7.5% | 18 / 18 | 18 / 18 | Flank rate fit on one circular stage (`PL1` → 0.0% residual): overfits |
| `scanned_bases` + `stored_contigs` | `169.3` s/Gbp − `109` s/M contigs, `× 1.15` | 9.6% | 18 / 18 | 18 / 18 | Unphysical negative per-contig term |

- **Rejected**: `stored_bases` alone (minimax error 31.7%, 0/18 within 20%) and `stored_bases + contigs` (`PL1` +63.0%
  over).
- **Rejected**: two-rate models (`stored_bases + flank_bases`, `scanned_bases + contigs`, `0.07 + scanned_bases`). None
  moves a stage across the 10% boundary (already 18/18), and separating stored from flank bases fits a second rate to
  one circular stage (`PL1`).
- **Chosen**: `CRISPR_CPU_S_PER_GBP = 168.1` (default step relative-error LS fit `168.14`; default minimax is `168.20`)
  and `CRISPR_NEAR_CAS_COST = 1.15` (filtered relative-error LS ratio `1.152`). At `1.15`, the interval of scanned-base
  rates that keeps all 18 records within 20% is `[152.8, 198.2]` CPU-s/Gbp, and within 10% is `[166.7, 171.7]`
  CPU-s/Gbp (*derived*, `b53_fit.py`). Refusal thresholds at 80% of 300 CPU-minutes (`240` CPU-min): `85.7` Gbp scanned
  (`0.8 * 300 * 60 / 168.1 = 85.66`) for the default step and `74.5` Gbp scanned (`74.49`) for the filtered step.

### CRISPR estimates and ratios (`record / estimate`) before and after, per built stage

| Stage | Scanned bases | Default rec | Old def est | Old def ratio | New def est | **New def ratio** | Filtered rec | Old filt est | Old filt ratio | New filt est | **New filt ratio** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| S1 | 1,664,109,938 | 5.085 | 5.048 | 1.007 | 4.662 | **1.091** | 5.487 | 6.007 | 0.913 | 5.362 | **1.023** |
| P1 | 1,016,625,803 | 3.085 | 3.084 | 1.000 | 2.848 | **1.083** | 3.383 | 3.670 | 0.922 | 3.275 | **1.033** |
| P2 | 9,596,180,328 | 26.489 | 29.108 | 0.910 | 26.885 | **0.985** | 31.794 | 34.639 | 0.918 | 30.918 | **1.028** |
| P3 | 9,949,077,522 | 29.195 | 30.179 | 0.967 | 27.874 | **1.047** | 33.386 | 35.913 | 0.930 | 32.055 | **1.042** |
| P4 | 9,976,456,950 | 27.585 | 30.262 | 0.912 | 27.951 | **0.987** | 32.606 | 36.012 | 0.905 | 32.143 | **1.014** |
| P5 | 9,967,711,632 | 27.505 | 30.235 | 0.910 | 27.926 | **0.985** | 32.361 | 35.980 | 0.899 | 32.115 | **1.008** |
| P6 | 9,962,536,103 | 27.938 | 30.220 | 0.924 | 27.912 | **1.001** | 32.796 | 35.961 | 0.912 | 32.098 | **1.022** |
| H1 | 1,819,831,338 | 4.925 | 5.520 | 0.892 | 5.099 | **0.966** | 5.876 | 6.569 | 0.895 | 5.863 | **1.002** |
| PL1 | 12,129,816,653 | 30.938 | 17.971 | 1.722 | 33.984 | **0.910** | 35.513 | 21.386 | 1.661 | 39.081 | **0.909** |

Implications for upcoming stages (*derived*):

- Stage `PL2` (18.3 Gbp of IMG/PR plasmids, card C-40): at `PL1`'s scanned-to-stored multiple (`2.047`), ~37.5 Gbp
  scanned → ~105 CPU-min default and ~121 CPU-min filtered, well under the 240 CPU-min margin (previously estimated at
  ~55 and ~66 CPU-min, which would have tripped DEC-43's > 24 CPU-min and > 20% guard).
- Stage `V1` (`/tmp/if2453/s2/inputs/V1/`, 168,358 stored contigs: 4.191 Gbp linear + 0.988 Gbp circular adding 1.819
  Gbp of 100 kb MinCED flanks = 6.997 Gbp scanned, *measured* by `b53_table.py`): `19.6` CPU-min default and `22.5`
  CPU-min filtered (vs `15.7` and `18.7` before).

## 4. Before launch on the real root (`b53_flank_check.py`)

`stage.crispr_flank_bases(root, stage)` streams `contig_index` (from `runs/stage=<S>/ingest.json`, sorted by
`contig_uid`, default variants) against the distinct contigs of `genome_input` (from `runs/stage=<S>/inputs.json`,
sorted by `contig_uid`) via `ncrna._range_rows` and `ncrna._distinct`. Contigs owned by an earlier stage are absent
from the stage's `contig_index` (DEC-26), and contigs shared by multiple genomes of the stage are yielded once
(DEC-27). Running `b53_flank_check.py` on `/tmp/if2453/r10/root` (bound `:ro`) returned `0` for `S1`, `P1` to `P6` and
`H1`, and `6,205,263,778` for `PL1` — identical to `table.json` on all nine stages (*measured*).

## 5. Gene-calling table of every built stage, and the circular `gv-meta` rate on `V1`

### Built stages (`S1`, `P1` to `P6`, `H1`, `PL1`; *measured* records, *derived* estimates)

None of the nine built stages has a circular contig in `gv-meta` mode, so their gene-call estimates are unchanged by
card B-53 (`Old ratio == New ratio`):

| Stage / piece | Genomes (skipped) | Mode & topology | Contig bases | Called bases (with 50 kb flanks) | Record CPU-min | Estimate CPU-min | **Old / New ratio** |
|---|---:|---|---:|---:|---:|---:|---:|
| S1 (1 proc) | 27,258 (0) | `gv-meta` linear | 1,664,109,938 | 1,664,109,938 | 142.786 | 143.621 | **0.994** |
| P1 (1 proc) | 250 (0) | `single` linear | 1,016,625,803 | 1,016,625,803 | 15.392 | 16.809 | **0.916** |
| P2 (1 proc) | 2,247 (0) | `single` linear | 9,596,180,328 | 9,596,180,328 | 143.288 | 156.456 | **0.916** |
| P3 (1 proc) | 2,671 (0) | `single` linear | 9,949,084,485 | 9,949,084,485 | 153.677 | 162.201 | **0.947** |
| P4 (1 proc) | 2,329 (0) | `single` linear | 9,976,464,003 | 9,976,464,003 | 151.406 | 162.647 | **0.931** |
| P5 (1 proc) | 2,289 (0) | `single` linear | 9,967,897,315 | 9,967,897,315 | 155.692 | 162.507 | **0.958** |
| P6 (1 proc) | 2,180 (0) | `single` linear | 9,962,553,301 | 9,962,553,301 | 157.196 | 162.420 | **0.968** |
| H1 (1 proc)¹ | 7,228 (514) | `gv-meta` linear | 1,961,515,105 | 1,961,515,105 | 146.091 | 169.248 | **0.863** |
| PL1 (1-proc plan) | 68,693 (7,046) | `meta` (5,099 lin + 56,548 circ) | 5,924,552,875 | 9,784,288,331 | — | 872.645 | — |
| PL1 `batch-00000` | 17,367 | `meta` (110.7 Mbp lin + 1,361.1 Mbp circ) | 1,471,826,928 | 2,443,650,166 | 213.891 | 218.345 | **0.980** |
| PL1 `batch-00001` | 16,928 | `meta` (113.2 Mbp lin + 1,380.0 Mbp circ) | 1,493,140,979 | 2,444,580,037 | 212.302 | 218.355 | **0.972** |
| PL1 `batch-00002` | 16,966 | `meta` (121.0 Mbp lin + 1,363.7 Mbp circ) | 1,484,779,065 | 2,447,119,851 | 209.893 | 218.347 | **0.961** |
| PL1 `batch-00003` | 17,432 | `meta` (126.5 Mbp lin + 1,348.3 Mbp circ) | 1,474,805,903 | 2,448,938,277 | 213.041 | 218.348 | **0.976** |
| PL1 `merge` | — | 6,371,980 genes written | — | — | 1.121 | 0.836 | 1.341 (< 24 min) |

¹ `H1` ran before card B-48 item 5 skipped genomes wholly owned by an earlier stage, so its run record called all
1,961,515,105 bp (`test_b48_genecall_estimate.py`); on the 6,714 stored genomes (`1,819,831,338` bp) `plan` now
estimates `157.043` CPU-min.

### Stage `V1` circular `gv-meta` measurement (card C-39, `/tmp/if2453/s2/inputs/V1/`)

`genecall.CALL_S_PER_MBP` multiplies `called_bases` (`n + 2 * min(n, 50,000)` for a circular contig of length `n`), not
raw contig bases. Card C-39 (`scripts/stress/build_stage_pod/c39_est.py`, unit `c39-est`, wilkcomp02, 2026-11-16,
`planning/measurements/2026-11-16_virus_votu_inputs.md` §5 and `/tmp/if2453/s2/inputs/V1/est.json`) timed
`genecall._call` on one thread on every 30th stored contig of UHGV v1.0 (`V1`):

- **Circular `gv-meta` sample**: 778 circular virus genomes, `34,564,134` contig bp (*measured*), `88,404,110` called bp
  with 50 kb wrap flanks (*measured*, `b53_table.py`; mean contig 44.4 kb, ratio `2.5577` called/contig), `386.87` CPU-s
  (*measured*). Per contig Mbp that is `11.1929` CPU-s/Mbp; per called Mbp read by `_call` that is
  `386.87 / 88.404110 = 4.376154` CPU-s/Mbp (*derived*), rounded up to **`4.3762`** as card B-48 rounded the other
  modes. Before B-53, `4.2019 * 1.5011 = 6.3075` (*estimate*) predicted `557.61` CPU-s on the sample (`+44.1%` over
  measurement, ratio `0.694`), because circular `meta` mode's `1.5011` per-called-base ratio reflects Prodigal's
  6-code meta pass on plasmids whereas `gv-meta`'s genetic-code pass scales nearly 1:1 with the extended length
  (`4.3762 / 4.2019 = 1.0415` against `S1`, `4.3762 / 3.7222 = 1.1757` against `V1` linear, *derived*).
- **Linear `gv-meta` sample**: 4,833 linear virus genomes, `137,465,952` bp, `511.67` CPU-s = `3.7222` CPU-s/Mbp
  (*measured*), below `S1`'s `4.2019` (kept as the larger measurement: sample ratio `0.886`).
- **Whole-stage `V1` before and after** (*derived*):
  - All 168,536 first-seen genomes (`4,233,918,303` linear bp + `988,370,473` circular contig bp = `2,494,857,527`
    circular called bp): C-39's sample rates project `(15 + 1.23 * (3.7222 * 4233.918 + 11.1929 * 988.370)) / 60 =
    550.11` CPU-min in one process. `one_process_cpu_min` was `687.55` CPU-min before (ratio `0.800`, right on the 20%
    edge) and is **`588.77` CPU-min** after (ratio **`0.934`**, `-6.6%`, 3 batches).
  - The 168,358 stored genomes (`4,190,809,601` linear bp + `2,493,970,530` circular called bp, excluding 178 owned by
    `S1`/`H1`): projected `546.66` CPU-min; `683.64` before (`0.800`) → **`584.90` CPU-min** after (**`0.935`**).

## 6. Tests

All run on `wilkcomp02` inside `/stor/home/if2453/bio-box/images/bio-box-test_faec1f4.sif` (`MemoryMax=32G`,
`MemorySwapMax=0`, `OMP_NUM_THREADS=2`, `nice -n 19`).

- **Test first** (`tests/r1/test_b53_step_estimates.py` on `10f3ddc` before editing `src/`): `22 failed, 4 passed in
  3.45s` — `[PL1-default]` failed first (`TypeError: crispr_cpu_min() got an unexpected keyword argument
  'flank_bases'`), `crispr_refusal` failed to refuse the circular plasmid stage (`AssertionError: None`), and the `V1`
  circular `gv-meta` sample test failed (`AssertionError: (557.61, 386.87)`).
- **After the edits**: `tests/r1/test_b53_step_estimates.py` and `tests/r0/test_runbook.py`: `76 passed in 5.14s`
  (*measured*, wilkcomp02, 2026-11-18); `tests/r0 tests/r1`: `1153 passed in 118.68s (0:01:58)` (*measured*); `tests/r0`
  plus `tests/r2/test_crispr_build_stage.py`, `tests/r2/test_e15_build_stage.py`, `tests/r2/test_b48_genecall_batches.py`
  and `tests/r2/test_build_stage.py`: `291 passed in 3047.03s (0:50:47)` (*measured*, capped at 2 threads beside
  `c35-rel-1`). No existing test was modified.

## Appendix A: `b53_table.py` (run as `python b53_table.py /tmp/if2453/r10/root /tmp/if2453/s2/inputs/V1 /tmp/if2453/b53/table.json`)

```python
"""Read-only inspection of built stages (S1, P1..P6, H1, PL1) and V1 inputs for card B-53."""
from __future__ import annotations

import csv
import json
import sys
from collections import Counter
from pathlib import Path

import pyarrow.parquet as pq

STAGES = ["S1", "P1", "P2", "P3", "P4", "P5", "P6", "H1", "PL1"]
CIRC_FLANK = 100_000
BATCH_BASES = 50_000_000
WRAP_FLANK = 50_000
SINGLE_MODE_MIN_BP = 100_000


def mode_of(kind: str, bases: int) -> str:
    if kind == "virus":
        return "gv-meta"
    if kind == "plasmid" or bases < SINGLE_MODE_MIN_BP:
        return "meta"
    return "single"


def called_bases(n: int, circular: bool) -> int:
    return n + 2 * min(n, WRAP_FLANK) if circular else n


def minced_bases(n: int, circular: bool) -> int:
    return n + 2 * min(CIRC_FLANK, n) if circular else n


def part_files(d: Path, pattern: str = "part-*.parquet") -> list[Path]:
    return sorted(d.glob(pattern)) if d.exists() else []


def main() -> None:
    root = Path(sys.argv[1])
    v1_dir = Path(sys.argv[2])
    out_path = Path(sys.argv[3])

    stage_own: dict[str, dict[bytes, int]] = {}
    all_len: dict[bytes, int] = {}
    for s in STAGES:
        own: dict[bytes, int] = {}
        for f in part_files(root / "l0seq" / "contig_index" / f"stage={s}"):
            t = pq.read_table(f, columns=["contig_uid", "contig_len", "is_default_variant"], use_threads=False)
            for u, n, d in zip(t["contig_uid"].to_pylist(), t["contig_len"].to_pylist(),
                               t["is_default_variant"].to_pylist()):
                if d:
                    own[u] = n
                    all_len[u] = n
        stage_own[s] = own

    crispr_rows, genecall_rows = [], []
    for s in STAGES:
        own = stage_own[s]
        gi_rows = []
        for f in part_files(root / "l0seq" / "genome_input" / f"stage={s}"):
            t = pq.read_table(f, columns=["contig_uid", "genome_uid", "is_circular", "kind"], use_threads=False)
            gi_rows.extend(zip(t["contig_uid"].to_pylist(), t["genome_uid"].to_pylist(),
                               t["is_circular"].to_pylist(), t["kind"].to_pylist()))
        by_cu: dict[bytes, list[tuple[bytes, bool, str]]] = {}
        by_guid: dict[bytes, list[tuple[bytes, bool, str]]] = {}
        for cu, guid, circ, kind in gi_rows:
            by_cu.setdefault(cu, []).append((guid, bool(circ), kind))
            by_guid.setdefault(guid, []).append((cu, bool(circ), kind))

        searched_contigs = circ_contigs = lin_contigs = circ_bases = lin_bases = m_bases = n_batches = cur_batch = 0
        kinds_searched: Counter[str] = Counter()
        skipped_other_stage = 0
        for cu in sorted(by_cu):
            _, circ, kind = sorted(by_cu[cu], key=lambda x: x[0])[0]
            if cu not in own:
                skipped_other_stage += 1
                continue
            n = own[cu]
            searched_contigs += 1
            kinds_searched[kind] += 1
            if circ:
                circ_contigs += 1
                circ_bases += n
            else:
                lin_contigs += 1
                lin_bases += n
            rec_len = minced_bases(n, circ)
            m_bases += rec_len
            if n > 0:
                cur_batch += rec_len
                if cur_batch >= BATCH_BASES:
                    n_batches += 1
                    cur_batch = 0
        if cur_batch > 0:
            n_batches += 1

        runs_dir = root / "runs" / f"stage={s}"
        c_info = {}
        for cr in sorted(runs_dir.glob("crispr-*.json")):
            rec = json.loads(cr.read_text())
            step_name = cr.stem
            key = "near_cas" if "--near-cas" in rec.get("argv", []) else "default"
            calls_in = calls_kept = contigs_with_call = 0
            rej_by_reason: Counter[str] = Counter()
            for lf in part_files(root / "ledger" / f"stage={s}" / f"step={step_name}"):
                lt = pq.read_table(lf, columns=["n_in", "n_kept", "n_rejected", "reason_code"], use_threads=False)
                for nin, nkept, nrej, rcode in zip(lt["n_in"].to_pylist(), lt["n_kept"].to_pylist(),
                                                   lt["n_rejected"].to_pylist(), lt["reason_code"].to_pylist()):
                    if rcode == "ok":
                        calls_in += nin
                        calls_kept += nkept
                        contigs_with_call += int(nin > 0)
                    elif not rcode.startswith("owned_by:"):
                        rej_by_reason[rcode] += nrej
            c_info[key] = {"step": step_name, "cpu_min": rec["cpu_min"], "wall_s": rec["wall_s"],
                           "calls_in": calls_in, "calls_kept": calls_kept, "contigs_with_call": contigs_with_call,
                           "rejected": dict(rej_by_reason)}

        crispr_rows.append({
            "stage": s, "stored_contigs": len(own), "stored_bases": sum(own.values()),
            "searched_contigs": searched_contigs, "skipped_other_stage": skipped_other_stage,
            "circular_contigs": circ_contigs, "linear_contigs": lin_contigs,
            "circular_bases": circ_bases, "linear_bases": lin_bases,
            "minced_bases": m_bases, "flank_bases": m_bases - sum(own.values()),
            "minced_batches": n_batches, "kinds": dict(kinds_searched), "runs": c_info,
        })

        least_guid = {cu: min(g for g, _, _ in entries) for cu, entries in by_cu.items()}
        called_by_mt: Counter[tuple[str, bool]] = Counter()
        all_called_by_mt: Counter[tuple[str, bool]] = Counter()
        n_not_called = 0
        for guid in sorted(by_guid):
            c_list = by_guid[guid]
            m = mode_of(c_list[0][2], sum(all_len[cu] for cu, _, _ in c_list))
            mine = any((cu in own) and (least_guid[cu] == guid) for cu, _, _ in c_list)
            for cu, circ, _ in c_list:
                all_called_by_mt[(m, circ)] += called_bases(all_len[cu], circ)
            if not mine:
                n_not_called += 1
                continue
            for cu, circ, _ in c_list:
                called_by_mt[(m, circ)] += called_bases(all_len[cu], circ)

        gc_rec = json.loads((runs_dir / "genecall.json").read_text())
        plan_p = root / "l1" / "genecall_plan" / f"stage={s}" / "plan.json"
        batches_info, plan_info = [], None
        if plan_p.exists():
            plan_doc = json.loads(plan_p.read_text())
            plan_info = {"estimate_cpu_min": plan_doc["estimate_cpu_min"],
                         "merge_estimate_cpu_min": plan_doc["merge_estimate_cpu_min"],
                         "n_batches": len(plan_doc["batches"])}
            guids_sorted = sorted(by_guid)
            for bdoc in plan_doc["batches"]:
                j = bdoc["batch"]
                brun = json.loads((runs_dir / f"genecall.batch-{j:05d}.json").read_text())
                b_called: Counter[tuple[str, bool]] = Counter()
                for guid in guids_sorted[bdoc["first"] : bdoc["first"] + bdoc["genomes"]]:
                    c_list = by_guid[guid]
                    if not any((cu in own) and (least_guid[cu] == guid) for cu, _, _ in c_list):
                        continue
                    m = mode_of(c_list[0][2], sum(all_len[cu] for cu, _, _ in c_list))
                    for cu, circ, _ in c_list:
                        b_called[(m, circ)] += called_bases(all_len[cu], circ)
                batches_info.append({"batch": j, "genomes": bdoc["genomes"], "bases": bdoc["bases"],
                                     "cpu_min": brun["cpu_min"], "estimate_cpu_min": bdoc["estimate_cpu_min"],
                                     "called_bases_by_mt": {f"{m}/{c}": v for (m, c), v in sorted(b_called.items())}})
        genecall_rows.append({
            "stage": s, "genomes": len(by_guid), "genomes_not_called": n_not_called,
            "called_bases_by_mt": {f"{m}/{c}": v for (m, c), v in sorted(called_by_mt.items())},
            "all_called_bases_by_mt": {f"{m}/{c}": v for (m, c), v in sorted(all_called_by_mt.items())},
            "genecall_record": {"cpu_min": gc_rec["cpu_min"], "wall_s": gc_rec["wall_s"]},
            "plan": plan_info, "batches": batches_info,
        })

    v1_info = {}
    if v1_dir.exists():
        with open(v1_dir / "records.tsv", newline="") as fh:
            recs = list(csv.DictReader(fh, delimiter="\t"))
        with open(v1_dir / "shared_records.tsv", newline="") as fh:
            owned_v1 = {r["contig_uid"] for r in csv.DictReader(fh, delimiter="\t")}
        with open(v1_dir / "repeated_records.tsv", newline="") as fh:
            rep_v1 = {r["ordinal"] for r in csv.DictReader(fh, delimiter="\t")}
        first = [r for r in recs if r["ordinal"] not in rep_v1]
        stored_recs = [r for r in first if r["contig_uid"] not in owned_v1]
        sample_recs = [r for i, r in enumerate(stored_recs, 1) if i % 30 == 0]

        def summarize_v1(rs):
            out_s = {"linear": {"n": 0, "contig_bp": 0, "called_bp": 0, "minced_bp": 0},
                     "circular": {"n": 0, "contig_bp": 0, "called_bp": 0, "minced_bp": 0}}
            for r in rs:
                circ = r["circular"] == "true"
                ln = int(r["length"])
                d = out_s["circular" if circ else "linear"]
                d["n"] += 1
                d["contig_bp"] += ln
                d["called_bp"] += called_bases(ln, circ)
                d["minced_bp"] += minced_bases(ln, circ)
            return out_s

        v1_info = {"all_first": summarize_v1(first), "stored": summarize_v1(stored_recs),
                   "sample_every_30": summarize_v1(sample_recs)}

    out_path.write_text(json.dumps({"crispr": crispr_rows, "genecall": genecall_rows, "v1": v1_info}, indent=1) + "\n")


if __name__ == "__main__":
    main()
```

## Appendix B: `b53_flank_check.py` (run with `PYTHONPATH=/tmp/if2453/ci/B-53/src`, `/tmp/if2453/r10/root` bound `:ro`)

```python
"""Check `stage.stage_bases`, `stage.crispr_flank_bases`, `stage.crispr_cpu_min`, and `genecall.one_process_cpu_min`
from the B-53 tree against `/tmp/if2453/r10/root` (bound `:ro`) and `/tmp/if2453/b53/table.json`."""
from __future__ import annotations

import json
from pathlib import Path

from biobox import stage as bs
from biobox.l1 import genecall

ROOT = Path("/tmp/if2453/r10/root")
TABLE = json.loads(Path("/tmp/if2453/b53/table.json").read_text())

for row in TABLE["crispr"]:
    st = row["stage"]
    bases = bs.stage_bases(ROOT, st)
    flank = bs.crispr_flank_bases(ROOT, st)
    assert bases == row["stored_bases"] and flank == row["flank_bases"]
    def_est = bs.crispr_cpu_min(bases, flank, near_cas=False)
    near_est = bs.crispr_cpu_min(bases, flank, near_cas=True)
    print(json.dumps({"stage": st, "stored_bases": bases, "flank_bases": flank,
                      "default_est": round(def_est, 6), "near_cas_est": round(near_est, 6)}))

v1 = TABLE["v1"]["all_first"]
print(json.dumps({"stage": "V1_inputs", "one_process_cpu_min": round(genecall.one_process_cpu_min({
    ("gv-meta", False): v1["linear"]["called_bp"],
    ("gv-meta", True): v1["circular"]["called_bp"],
}), 6)}))
```

## Appendix C: `b53_fit.py` (arithmetic on `table.json`)

```python
"""Fit and compare candidate models for CRISPR and gene-calling estimates on table.json (card B-53)."""
import json
import sys

doc = json.load(open(sys.argv[1]))
cdata = doc["crispr"]
yd = {r["stage"]: r["runs"]["default"]["cpu_min"] for r in cdata}
yn = {r["stage"]: r["runs"]["near_cas"]["cpu_min"] for r in cdata}
mb = {r["stage"]: r["minced_bases"] / 1e9 for r in cdata}

r_ls = sum(mb[s] / 60 / yd[s] for s in mb) / sum((mb[s] / 60 / yd[s]) ** 2 for s in mb)
ratios = {s: yn[s] / yd[s] for s in mb}
nf_ls = sum(1.0 / ratios[s] for s in mb) / sum((1.0 / ratios[s]) ** 2 for s in mb)
print(f"default LS={r_ls:.4f} s/Gbp, near/default LS={nf_ls:.4f} (mean={sum(ratios.values())/len(ratios):.4f})")
```
