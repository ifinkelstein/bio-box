# Card C-38c: Linclust alone against our cascade on `r26.11.1`'s proteins

- **Date:** 2026-11-24
- **Host:** `wilkcomp02.ccbb.utexas.edu` (`MemoryMax=64G`, `MemorySwapMax=0`, `CPUQuota=1200%`, `--threads 12`, `nice 19`, detached `systemd-run --user --unit=c38c-*`)
- **Image:** `/stor/home/if2453/bio-box/images/bio-box-test_faec1f4.sif` (`MMseqs2 18.8cc5c`)
- **Corpus (read-only):** `r26.11.1` in `/tmp/if2453/r10/root` — `49,102,041` distinct proteins (`33,117,773` cascade NR90 clusters, `9,291,312` cascade NR50 clusters, *measured*).
- **Supporting files:** [`linclust_calls.tsv`](2026-11-24_c38c_linclust_vs_cascade/linclust_calls.tsv), [`compare.json`](2026-11-24_c38c_linclust_vs_cascade/compare.json), [`nbr_comparison.json`](2026-11-24_c38c_linclust_vs_cascade/nbr_comparison.json), [`n19_n20_checks.json`](2026-11-24_c38c_linclust_vs_cascade/n19_n20_checks.json), [`extrapolations.json`](2026-11-24_c38c_linclust_vs_cascade/extrapolations.json).
- **Decisions informed:** `DEC-54` (Q9), `DEC-75` (3).

---

## Executive summary & `AGENTS.md` §6 (20% rule) flags

1. **Speed:** Linclust alone across both levels (`NR90` + `rep90` + `NR50` + `mergeclusters`) took **`767.74` CPU-minutes** (*measured*, `15.64` CPU-min/M proteins; `782.87` CPU-min including `prepare` + `createtsv`) versus **`14,644.04` CPU-minutes** for `r26.11.1`'s published cascade clustering (`15,255.59` CPU-min total including `plan`, `members`, `finish`, *measured* from `runs/stage=r26.11.1/cluster*.json`) — a **`19.1×` speedup** (*derived*). Every individual process stayed `< 86` CPU-minutes (with `nr90.07.align` and `nr50.07.align` split into 4 parts via `B-41`'s key-range cut).
2. **>20% deviations from design note `N19` and `N20` (`AGENTS.md` §6):**
   - **FLAG — `N20` (linclust memory per protein): EXCEEDS 20% RULE (`+128.9%` by `ru_maxrss`, `+61.9%` by cgroup `MemoryPeak`).** Peak memory occurred in `linclust-nr90-00-kmermatcher`: `ru_maxrss = 37.77 GB` (`35.17 GiB`, **`769.1 B/protein`**, *measured*) vs `N20`'s `336 B/seq` (*estimate*, `2.29×` higher). Even cgroup `MemoryPeak` (`26.71 GB` = `544.0 B/protein`, *measured*, excluding pre-cached `mmap` pages of `db`) is `1.62×` (`+61.9%`) `N20`. This matches `r26.11.1`'s own `cluster.nr90.00.kmermatcher.json` (`37.89 GB`, `771.7 B/protein`, *measured*).
   - **FLAG — `N19` (NR90 linclust share of NR90 cascade): EXCEEDS 20% RULE (`-30.9%`).** NR90 linclust (`409.65` CPU-min, *measured*) is **`5.11%`** of `r26.11.1`'s NR90 cascade (`8,014.42` CPU-min, *measured*), versus `N19`'s `7.4%` (*derived* at `10.6M` proteins), because the cascade's prefilter grows superlinearly (`n^1.244`) while linclust grows near-linearly.
   - *Within 20%:* Both-levels linclust share (`5.24%` of `14,644.04` CPU-min vs `N19`'s `6.0%`: `-12.6%`) and both-levels total CPU (`767.74` CPU-min vs `N19`'s `14.0 × 49.102M = 687.43` CPU-min: `+11.7%`).
3. **Clustering quality:**
   - **NR90 is nearly lossless:** `+2.31%` more clusters (`33,882,143` vs `33,117,773`) and `+3.00%` more singletons (`28,387,341` vs `27,561,851`, *measured*). `0` of `15,219,898` linclust NR90 joins disagree with the cascade (`0.00%`).
   - **NR50 fragments substantially below 70% identity:** `+60.84%` more clusters (`14,943,715` vs `9,291,312`, `1.608×`) and `+71.68%` more singletons (`10,186,340` vs `5,933,293`, `1.717×`, *measured*). Among cascade NR50 clusters of `size >= 100` (`64,654` clusters holding `38.0%` of all proteins), **`98.88%` split**, into a mean of **`18.35`** linclust clusters (median `15`, p95 `45`, max `303`, *measured*). Linclust misses `71.99%` of cascade joins in the `50-60%` identity band and `47.57%` in `60-70%`, versus only `6.36%` in `80-90%` and `0.71%` at `>=90%`.
4. **Worked examples (`nbr summarize --by nr50`):**
   - On **TniQ** (`4,545` loci, `10g` / `20kb`) and **RNA_pol** (`2,667` loci, `10g` / `5g`), total distinct NR50 neighbour keys grow by **`+10.9%` to `+20.0%`** (`1.109×`–`1.200×`, *measured*), while Pfam keys are 100% identical (`5,613`, `7,381`, `1,068`, `567`).
   - Conserved top partners retain their top ranks (`TnsB_C` #1 in TniQ: `637 -> 575` loci, `-9.7%`; `DNA_ligase_A_M` #1 in RNA_pol: `665 -> 654`, `-1.7%`; `Phage_T7_1_1` #2: `493 -> 491`, `-0.4%`), but divergent families suffer sharp single-cluster locus dilution (`DDE_Tnp_IS240` `363 -> 231`, `-36.4%`; `Response_reg` `214 -> 123`, `-42.5%`; `Phage_T7_kinase` `492 -> 298`, `-39.4%`; `DUF3310` `367 -> 229`, `-37.6%`; `S-AdoMet_lyase` `253 -> 166`, `-34.4%`).
5. **Extrapolation (`350M` and `1B` proteins, *derived*):**
   - At **350M proteins**, linclust needs `5,472–6,069` CPU-min (`91–101` CPU-h, `28–31×` faster than `N17`'s `169,809` CPU-min) and **`250.7 GiB`** peak RSS (`ru_maxrss` linear; `177.3 GiB` cgroup) — **above** `N20`'s `118 GB` estimate and a 200 GB scope, but well within `wilkcomp02`'s `750G` scope (`DEC-73`).
   - At **1B proteins**, linclust needs `15,636–18,285` CPU-min (`261–305` CPU-h, `34–40×` faster than `N17`'s `626,881` CPU-min) and **`716.3 GiB`** peak RSS (`ru_maxrss` linear; `506.7 GiB` cgroup) — right at the ceiling of `wilkcomp02`'s `750G` scope (`DEC-73`).

---

## 1. Per-call runtime and peak memory on `r26.11.1` (`49,102,041` proteins)

All calls ran `l2/cluster.workflow`'s exact 10 redundancy-phase calls (`kmermatcher`, `rescorediagonal`, `clust`, `createsubdb`, `filterdb`, `align`, `mergeclusters`; `--shuffle 0`, `-c 0.8`, `--cov-mode 1`, `--cluster-mode 2`, `--min-seq-id 0.9` then `0.5` over `rep90`, followed by `mergeclusters db clu50 clu90 clu50r`). Both `07.align` steps were cut into 4 key-range parts (`cluster.cut`) to stay strictly below the 300 CPU-minute single-process cap (`AGENTS.md` §7.2).

| Step | Wall (s) | Cgroup CPU (min) | Cgroup MemPeak (GB) | `ru_maxrss` (GB) | Output keys | Label |
|---|---:|---:|---:|---:|---:|---|
| `prepare` (FASTA + `createdb --shuffle 0`) | 816 | 8.28 | 48.56 | — | 49,102,041 | *measured* |
| `linclust-nr90-00-kmermatcher` | 216 | 30.24 | **26.71** | **37.77** | 49,102,041 | *measured* |
| `linclust-nr90-01-rescorediagonal` | 138 | 13.50 | 6.77 | 22.19 | 49,102,041 | *measured* |
| `linclust-nr90-02-clust` | 88 | 4.30 | 4.83 | 4.73 | 35,153,162 | *measured* |
| `linclust-nr90-03-createsubdb` | 18 | 0.34 | 2.23 | 2.36 | 35,153,162 | *measured* |
| `linclust-nr90-04-createsubdb` | 32 | 0.91 | 2.20 | 2.30 | 35,153,162 | *measured* |
| `linclust-nr90-05-filterdb` | 87 | 9.72 | 7.16 | 7.65 | 35,153,162 | *measured* |
| `linclust-nr90-06-rescorediagonal` | 81 | 9.59 | 5.62 | 19.84 | 35,153,162 | *measured* |
| `nr90-07-align` (cut + 4 parts + join) | 2,008 | 330.66 | 3.94 | 18.27 | — | *measured* (max part: `85.41` CPU-min) |
| `linclust-nr90-08-clust` | 40 | 2.62 | 3.72 | 4.34 | 33,882,143 | *measured* |
| `linclust-nr90-09-mergeclusters` (`clu90`) | 78 | 7.76 | 7.31 | 5.84 | 33,882,143 | *measured* |
| **NR90 linclust subtotal** | **2,786** | **409.65** | **26.71** | **37.77** | **33,882,143** | *derived* (sum/max) |
| `cluster-rep90` (`createsubdb`) | 50 | 0.80 | 13.30 | 16.94 | 33,882,143 | *measured* |
| `linclust-nr50-00-kmermatcher` | 161 | 19.78 | 21.60 | 26.02 | 33,882,143 | *measured* |
| `linclust-nr50-01-rescorediagonal` | 142 | 15.79 | 5.22 | 18.30 | 33,882,143 | *measured* |
| `linclust-nr50-02-clust` | 86 | 3.42 | 3.31 | 3.63 | 18,970,453 | *measured* |
| `linclust-nr50-03-createsubdb` | 14 | 0.25 | 1.46 | 1.62 | 18,970,453 | *measured* |
| `linclust-nr50-04-createsubdb` | 27 | 0.87 | 1.46 | 1.60 | 18,970,453 | *measured* |
| `linclust-nr50-05-filterdb` | 77 | 9.39 | 5.18 | 7.81 | 18,970,453 | *measured* |
| `linclust-nr50-06-rescorediagonal` | 63 | 7.68 | 3.30 | 14.11 | 18,970,453 | *measured* |
| `nr50-07-align` (cut + 4 parts + join) | 1,666 | 289.80 | 2.59 | 12.75 | — | *measured* (max part: `73.23` CPU-min) |
| `linclust-nr50-08-clust` | 27 | 1.50 | 2.04 | 2.73 | 14,943,715 | *measured* |
| `linclust-nr50-09-mergeclusters` (`clu50r`) | 43 | 4.07 | 4.58 | 3.70 | 14,943,715 | *measured* |
| **NR50 linclust subtotal** | **2,306** | **352.55** | **21.60** | **26.02** | **14,943,715** | *derived* (sum/max) |
| `cluster-merge` (`mergeclusters` -> `clu50`) | 54 | 4.74 | 6.05 | 5.46 | 14,943,715 | *measured* |
| **Linclust core total (NR90 + rep90 + NR50 + merge)** | **5,196** | **767.74** | **26.71** | **37.77** | **14,943,715** | *derived* |
| `createtsv` (`clu90.tsv` + `clu50.tsv`) | 116 | 6.85 | 12.61 | — | — | *measured* |

*Note on `ru_maxrss` vs cgroup `MemoryPeak`:* `kmermatcher` mmaps the `db` sequence database (`~17 GB`). Because `prepare` had just written `db` into Linux page cache in an earlier cgroup, `kmermatcher`'s own cgroup was charged `26.71 GB` (`544.0 B/prot`), whereas `ru_maxrss` (`37.77 GB`, `769.1 B/prot`) includes resident mapped pages of `db` and is the conservative figure for cold-cache sizing (`r26.11.1`'s own `cluster.nr90.00.kmermatcher.json` recorded `37.89 GB`).

---

## 2. Comparison against `r26.11.1`'s published cascade

### 2.1 Clusters and singletons per level (*measured*, bounded DuckDB `memory_limit='8GB'`)

| Level | Metric | Cascade (`r26.11.1`) | Linclust alone | Delta (`lin - cas`) | Ratio (`lin / cas`) |
|---|---|---:|---:|---:|---:|
| **NR90** | Total clusters | 33,117,773 | 33,882,143 | +764,370 | 1.0231× (+2.31%) |
| **NR90** | Singletons | 27,561,851 | 28,387,341 | +825,490 | 1.0300× (+3.00%) |
| **NR90** | Singleton share of clusters | 83.22% | 83.78% | +0.56 pp | — |
| **NR90** | Singleton share of proteins | 56.13% | 57.81% | +1.68 pp | — |
| **NR50** | Total clusters | 9,291,312 | 14,943,715 | +5,652,403 | **1.6084× (+60.84%)** |
| **NR50** | Singletons | 5,933,293 | 10,186,340 | +4,253,047 | **1.7168× (+71.68%)** |
| **NR50** | Singleton share of clusters | 63.86% | 68.16% | +4.31 pp | — |
| **NR50** | Singleton share of proteins | 12.08% | 20.75% | +8.66 pp | — |

Reverse-join sanity check (*measured*): at NR90, `0 / 15,219,898` (`0.000%`) linclust non-representative joins fell in different cascade clusters; at NR50, `127,013 / 34,158,326` (`0.372%`) differed (due to the 2.3% larger `rep90` input feeding NR50 linclust).

### 2.2 Fragmentation of cascade NR50 clusters into linclust NR50 clusters (*measured*)

| Cascade NR50 stratum | Cascade clusters | Proteins in stratum | Unsplit (`k=1`) | Split (`k>=2`) | Split share | Mean `k` | Median `k` | p95 `k` | Max `k` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `all` (`size >= 1`) | 9,291,312 | 49,102,041 | 7,231,626 | 2,059,686 | 22.17% | 1.61 | 1 | 4 | 303 |
| `size >= 2` | 3,358,019 | 43,168,748 | 1,298,333 | 2,059,686 | **61.34%** | 2.70 | 2 | 7 | 303 |
| `size >= 10` | 652,439 | 33,828,126 | 61,470 | 590,969 | **90.58%** | 6.45 | 5 | 18 | 303 |
| `size >= 100` | 64,654 | 18,660,161 | 725 | 63,929 | **98.88%** | **18.35** | **15** | **45** | 303 |

Histogram of linclust cluster count `k` per cascade NR50 cluster (`all`, *measured*):
- `k = 1`: `7,231,626` clusters (`10,769,200` proteins, `21.9%` of corpus)
- `k = 2`: `1,104,398` clusters (`5,213,045` proteins, `10.6%`)
- `k = 3–5`: `680,506` clusters (`8,620,049` proteins, `17.6%`)
- `k = 6–10`: `181,996` clusters (`7,493,417` proteins, `15.3%`)
- `k = 11–50`: `90,461` clusters (`14,796,132` proteins, `30.1%`)
- `k > 50`: `2,325` clusters (`2,210,198` proteins, `4.5%`)

### 2.3 Cascade joins missed by linclust, by `cluster_member_aln.pident` identity band (*measured*)

A cascade non-representative member (`member_uid != rep_uid`) is counted as **missed by linclust** when it is not in the same linclust cluster as its cascade representative (`lin_rep(member) != lin_rep(cascade_rep)`).

**NR50 (`39,810,729` cascade non-rep members):**

| Identity band | Cascade non-rep | Joined by linclust | Missed by linclust | Missed share | Missed -> singleton | Singleton share of missed |
|---|---:|---:|---:|---:|---:|---:|
| `<50%` | 3,122,229 | 301,802 | 2,820,427 | **90.33%** | 218,415 | 7.74% |
| `50-60%` | 16,679,576 | 4,672,353 | 12,007,223 | **71.99%** | 2,141,303 | 17.83% |
| `60-70%` | 9,550,468 | 5,007,305 | 4,543,163 | **47.57%** | 701,172 | 15.43% |
| `70-80%` | 5,265,449 | 4,114,922 | 1,150,527 | 21.85% | 208,588 | 18.13% |
| `80-90%` | 3,034,061 | 2,841,185 | 192,876 | 6.36% | 45,953 | 23.83% |
| `>=90%` | 2,158,946 | 2,143,667 | 15,279 | 0.71% | 4,566 | 29.88% |
| **Total** | **39,810,729** | **19,081,234** | **20,729,495** | **52.07%** | **3,319,997** | **16.02%** |

**NR90 (`15,984,268` cascade non-rep members):**

| Identity band | Cascade non-rep | Joined by linclust | Missed by linclust | Missed share | Missed -> singleton |
|---|---:|---:|---:|---:|---:|
| `<80%` (merged/transitive) | 2,748 | 2,381 | 367 | 13.36% | 19 |
| `80-90%` | 378,635 | 232,409 | 146,226 | 38.62% | 2,007 |
| `>=90%` | 15,602,885 | 14,452,474 | 1,150,411 | **7.37%** | 578,931 |
| **Total** | **15,984,268** | **14,687,264** | **1,297,004** | **8.11%** | **580,957** |

### 2.4 Explicit `AGENTS.md` §6 (20% rule) audit against `N19` and `N20`

| Design-note item | Estimate (`C-38_import_depths.md`) | Measured (`r26.11.1`, `49.10M`) | Deviation | >20% flag |
|---|---:|---:|---:|---|
| `N19`: NR90 linclust share of NR90 cascade CPU | `7.40%` | `5.11%` (`409.65 / 8,014.42`) | **`-30.93%`** | **YES (>20%)** |
| `N19`: Both-levels linclust share of cascade CPU | `6.00%` | `5.24%` (`767.74 / 14,644.04`) | `-12.62%` | No |
| `N19`: Both-levels CPU-min at `14.0`/M | `687.43` | `767.74` (`15.64`/M) | `+11.68%` | No |
| `N20`: Peak memory per protein (`ru_maxrss`) | `336.0 B/seq` | **`769.1 B/prot`** (`37.77 GB`) | **`+128.91%`** | **YES (>20%)** |
| `N20`: Peak memory per protein (cgroup `MemoryPeak`) | `336.0 B/seq` | **`544.0 B/prot`** (`26.71 GB`) | **`+61.91%`** | **YES (>20%)** |

---

## 3. Worked examples: `nbr summarize --by nr50` under cascade vs linclust

Both worked-example searches ran read-only against `r26.11.1` (`TniQ` PF06527 at GA `27.00` bits: `4,545` loci of `3,164` proteins; `RNA_pol` PF00940 at `30` bits: `2,667` loci of `2,037` proteins — matching `examples/*/recorded/r26.11.1.tsv` exactly). `bio-box nbr summarize --by pfam,nr50 --top 0` then ran on `/tmp/if2453/r10/root` (cascade) and on a scratch release view symlinking `r26.11.1`'s read-only stage/dim files with linclust's `cluster_member` + `member_cluster`. Pfam key counts were 100% identical across both views in all four windows.

### 3.1 Summary of NR50 key expansion (*measured*)

| Search & window | Anchor loci | Neighbours | Pfam keys (cas / lin) | NR50 keys (cascade) | NR50 keys (linclust) | Key ratio |
|---|---:|---:|---:|---:|---:|---:|
| TniQ `10g` | 4,545 | 90,098 | 5,613 / 5,613 | 39,161 | 43,415 | **1.109× (+10.9%)** |
| TniQ `20kb` | 4,545 | 174,549 | 7,381 / 7,381 | 70,995 | 79,600 | **1.121× (+12.1%)** |
| RNA_pol `10g` | 2,667 | 47,600 | 1,068 / 1,068 | 7,270 | 8,724 | **1.200× (+20.0%)** |
| RNA_pol `5g` | 2,667 | 25,321 | 567 / 567 | 4,008 | 4,768 | **1.190× (+19.0%)** |

### 3.2 Partner fragmentation in TniQ (`10g`, top neighbour families by cascade NR50 loci, *measured*)

| Dominant Pfam | Cascade NR50 clusters | Linclust NR50 clusters | Top cluster loci (cas -> lin) | Top share (cas -> lin) | Top-cluster loci delta |
|---|---:|---:|---:|---:|---:|
| `TnsB_C` | 327 | 386 | 637 -> 575 | 14.02% -> 12.65% | -9.7% |
| `DDE_Tnp_IS240` | 33 | 40 | 363 -> 231 | 7.99% -> 5.08% | **-36.4%** (rank #2 -> #6) |
| `Resolvase` | 115 | 152 | 353 -> 307 | 7.77% -> 6.75% | -13.0% |
| `MerR_1` | 26 | 17 | 311 -> 310 | 6.84% -> 6.82% | -0.3% (rank #4 -> #2) |
| `DEDD_Tnp_IS110` | 53 | 76 | 237 -> 237 | 5.21% -> 5.21% | 0.0% |
| `Response_reg` | 116 | 133 | 214 -> 123 | 4.71% -> 2.71% | **-42.5%** (drops out of top 15) |
| `Tn7_TnsA-like_N` | 295 | 374 | 208 -> 208 | 4.58% -> 4.58% | 0.0% |
| `AAA_22` | 622 | 844 | 207 -> 207 | 4.55% -> 4.55% | 0.0% |
| `(unannotated)` | 12,595 | 13,967 | 207 -> 207 | 4.55% -> 4.55% | 0.0% |

### 3.3 Partner fragmentation in RNA_pol (`10g` and `5g`, *measured*)

| Dominant Pfam (`10g`) | Cascade NR50 clusters | Linclust NR50 clusters | Top cluster loci (cas -> lin) | Top share (cas -> lin) | Top-cluster loci delta |
|---|---:|---:|---:|---:|---:|
| `DNA_ligase_A_M` | 28 | 40 | 665 -> 654 | 24.93% -> 24.52% | -1.7% |
| `Phage_T7_1_1` | 21 | 34 | 493 -> 491 | 18.49% -> 18.41% | -0.4% |
| `Phage_T7_kinase` | 22 | 31 | 492 -> 298 | 18.45% -> 11.17% | **-39.4%** (rank #3 -> #6; splits into 298 + 194) |
| `Amidase_2` | 20 | 31 | 414 -> 338 | 15.52% -> 12.67% | **-18.4%** |
| `BC10` | 3 | 12 | 386 -> 301 | 14.47% -> 11.29% | **-22.0%** |
| `SBB_T7` | 38 | 52 | 369 -> 339 | 13.84% -> 12.71% | -8.1% |
| `(unannotated)` | 5,002 | 5,862 | 368 -> 253 | 13.80% -> 9.49% | **-31.3%** |
| `DUF3310` | 34 | 54 | 367 -> 229 | 13.76% -> 8.59% | **-37.6%** |
| `Terminase-T7_RNaseH-like` | 6 | 10 | 355 -> 250 | 13.31% -> 9.37% | **-29.6%** |
| `S-AdoMet_lyase` | 20 | 29 | 253 -> 166 | 9.49% -> 6.22% | **-34.4%** |

At `5g` (`RNA_pol` standard window), the exact same pattern holds: `DNA_ligase_A_M` (`663 -> 652`) and `Phage_T7_1_1` (`493 -> 491`) stay intact at #1–#2, while `Phage_T7_kinase` splits `492 -> 298` (`+194` in a second cluster at rank #12), `BC10` drops `321 -> 251` (`3 -> 11` clusters), and `S-AdoMet_lyase` drops `230 -> 146` (`16 -> 25` clusters).

---

## 4. Extrapolation to `350,000,000` and `1,000,000,000` proteins (*derived* / *extrapolation*)

Scaled from `r26.11.1`'s measured linclust core (`767.74` CPU-min at `n = 49,102,041`; peak `ru_maxrss = 37.77 GB` = `769.1 B/prot`; cgroup `MemoryPeak = 26.71 GB` = `544.0 B/prot`). Cascade references use `N17` (`284.10 × (n / 2,053,075)^1.244` CPU-min) and `N18` (`178 GB × n / 45,815,843`).

| Target proteins | Scale vs `r26.11.1` | Linclust CPU-min (linear / `n log n`) | Linclust peak RSS (`ru_maxrss` / cgroup, linear) | Min parts (<240 CPU-min) | Cascade `N17` CPU-min | Cascade `N18` peak RSS | Speedup (lin / `n log n`) | Label |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `49,102,041` (`r26.11.1`) | 1.00× | `767.7` | `35.2 GiB` / `24.9 GiB` | 4 | `14,644` (*meas*) | `178 GB` (`rM.2`, *meas*) | `19.1×` | *measured* |
| **`350,000,000`** | 7.13× | **`5,472` / `6,069`** (`91–101` CPU-h) | **`250.7 GiB`** / `177.3 GiB` | 23–26 | `169,809` (`2,830` CPU-h) | `1,266 GiB` (`1.36 TB`) | `31.0×` / `28.0×` | *derived (extrapolation)* |
| **`1,000,000,000`** | 20.37× | **`15,636` / `18,285`** (`261–305` CPU-h) | **`716.3 GiB`** / `506.7 GiB` | 66–77 | `626,881` (`10,448` CPU-h) | `3,618 GiB` (`3.88 TB`) | `40.1×` / `34.3×` | *derived (extrapolation)* |

---

## 5. Practical takeaway for `DEC-54` (Q9) and `DEC-75` (3)

- **Keep the cascade while `n <= ~200M` (`r26.11.4` under `DEC-75`):** Cascade NR50 preserves `20.7M` distant (`50–70%` identity) joins that linclust misses, avoiding `18×` average fragmentation of `size >= 100` families and `20–42%` locus-count dilution in divergent neighbour families (`DDE_Tnp_IS240`, `Phage_T7_kinase`, `DUF3310`, `Response_reg`).
- **Switch to linclust when crossing `250–350M+` proteins (or consider hybrid NR90-linclust + NR50-cascade below 200M):**
  - **NR90 linclust is essentially free** (`+2.3%` clusters, `5.1%` of NR90 cascade CPU, saving `~7,600` CPU-min per release at 49M proteins!).
  - Full two-level linclust scales to `350M` (`~95` CPU-h, `~251 GiB` peak RSS) and up to `~1B` proteins (`~280` CPU-h, `~716 GiB` peak RSS, right at `wilkcomp02`'s `750G` scope) while keeping top-ranked conserved biological neighbours (`TnsB_C`, `DNA_ligase_A_M`, `Phage_T7_1_1`, `Tn7_TnsA-like_N`, `MerR_1`) clearly discoverable — at the documented cost of `1.61×` more NR50 clusters and `10–42%` top-cluster dilution on divergent families.
