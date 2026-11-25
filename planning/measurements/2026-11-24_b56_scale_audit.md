# Pure-Code Scale Audit: `B-56` vs `B-66` Across `r26.11.4` (`4.35M -> 14.48M` Contigs)

- **Date**: 2026-11-24
- **Scope**: Pure-code static audit of every `KeyLookup`, `SeqStore`, and `ContigOwners` lookup path across `ingest`, `genecall`, `crispr`, and `release` (`tapes`, `tokens`, `postings`, `dim`) as the root grows from `r26.11.3` (`13` stages, `4,350,219` stored contigs, `46` row groups) across the `19`-stage `r26.11.4` linear chain (`+10,134,020` net new contigs -> `14,484,239` cumulative contigs, `32` stages, `141` cumulative row groups).
- **Compute**: Zero pod compute, zero SSH, zero local `pytest`/`python` runs (`AGENTS.md` §5.5).
- **Labels (`AGENTS.md` §6)**: Every number is labelled *measured* (with source note, host, and date) or *derived* (with explicit formula and inputs).

---

## 0. Executive Verdict: **NO-GO for `GART1`/`GART2` without a 3-line fix in `genecall._genome_bases`; GO immediately after**

1. **`B-66` completely protects `crispr` and owned-contig `SeqStore.get` in `genecall`**: `KeyLookup.find_first` (`src/biobox/l0seq/lookup.py:117-131`) + single-locate `SeqStore.get` (`src/biobox/l0seq/store.py:1646-1657, 1697-1718`) stop at `stages[0] == stage`, so owned-contig lookups only touch the current stage's own `contig_index` and `seq_locator` (`<= 43` row groups on `GART1` `<= 64` LRU capacity).
2. **`B-66` missed one `KeyLookup.get` call in `genecall._genome_bases` (`src/biobox/l1/genecall.py:562, 868`)**: `lengths = KeyLookup(store.index_files, "contig_uid", ["contig_len", "is_default_variant"])` still calls **`lengths.get(cu)`** instead of `lengths.find_first(cu, lambda r: r["is_default_variant"])`. On **every contig of every owned genome**, `lengths.get(cu)` scans **all `S` stages** in the root and touches **all `R_cum` cumulative `contig_index` row groups** through a single `64`-entry LRU cache.
3. **Why this stops `r26.11.4` at `GART1` and `GART2`**:
   - Cumulative `contig_index` row groups cross `64` at stage 25 (`LRV1`, `65` RGs), reach `88` at `LR1`, **`131` at `GART1`**, and **`141` at `GART2`** (`22` single-RG stages permanently pin `22` of the `64` LRU slots, leaving only **`42` slots** for `109–119` multi-RG stage row groups).
   - Because `GART1` (`9,234` MAGs, mean `565` contigs/MAG) and `GART2` (`15,320` MAGs, mean `79` contigs/MAG) iterate each MAG's `contig_uids` in sorted `contig_uid` order (`inputs.py:139`), every MAG sweeps `109–111.5` multi-RG row groups `> 42` free slots, causing **100% LRU eviction between MAGs**:
     - **`GART1`** (`3` batches, `213.8` CPU-min/batch est.): **`~1.01M` RG decodes** + `30`-stage scan adds **`+121 to +238` CPU-min total (`+40 to +79` CPU-min/batch)** (*derived* at `3–10 ms/RG`), risking the **`280` CPU-min hard process cap** (`DEC-64`) at `>= 7.5 ms/RG`.
     - **`GART2`** (`1` batch, `192.6` CPU-min est.): **`~1.71M` RG decodes** + `31`-stage scan adds **`+102 to +302` CPU-min** (*derived* at `3–10 ms/RG`), pushing the single process to **`295–494` CPU-min (`1.53×–2.57×`)** — tripping **both** the `280` CPU-min process cap and the `1.50×` `DEC-64` rate guard even at the conservative `3 ms/RG` lower bound.
4. **Minimal fix (`3` lines in `src/biobox/l1/genecall.py:562-565`)**: replace `lengths.get(cu)` with `lengths.find_first(cu, lambda r: r["is_default_variant"])`. Output tables, ledger, plan hashes, and batch digests are **100% bit-for-bit identical**. Once landed, `B-56` only bites when a **single** stage exceeds `64 × 122,880 = 7,864,320` stored contigs (`GART1` is `5.19M` = `43` RGs, leaving `21` RGs of headroom).

---

## 1. Cumulative Row-Group Growth Across `r26.11.3` + `r26.11.4` (`ROW_GROUP_ROWS = 122,880`)

Per `src/biobox/io.py:35`, `ROW_GROUP_ROWS = 122_880`. Every stage with `N >= 1` stored contigs writes at least `ceil(N / 122,880)` row groups in both `contig_index` and `seq_locator`. Base counts for `r26.11.3` (`13` stages) are *measured* from [`2026-11-22_b66_seqstore_and_crispr_contigs.md`](2026-11-22_b66_seqstore_and_crispr_contigs.md) §4–§5 (`wilkcomp02`, 2026-11-22); counts for the `19` `r26.11.4` stages are *measured* from [`2026-11-24_a26_r26_11_4_confs.md`](2026-11-24_a26_r26_11_4_confs.md) §3 (`wilkcomp02`, 2026-11-24). Row-group counts and cumulative totals are *derived* (`ceil(new_contigs / 122_880)`).

| # | Stage | Net New Contigs (*measured*) | Stage RGs (*derived*) | Cum. Contigs (*derived*) | Cum. RGs (*derived*) | 1-RG Stages Pinned | Multi-RG Pool vs Free LRU (`64 - pinned`) |
|---:|---|---:|---:|---:|---:|---:|---|
| 1–10 | `S1, P1..P6, H1, PL1, C5` | `125,279` | `10` (`1` each) | `125,279` | `10` | `10` | `0` vs `54` |
| 11 | `V1` | `168,358` | `2` | `293,637` | `12` | `10` | `2` vs `54` |
| 12 | `G1` | `1,227,214` | `10` | `1,520,851` | `22` | `10` | `12` vs `54` |
| 13 | `SPX1` (`r26.11.3` end) | `2,829,368` | `24` | `4,350,219` | **`46`** | `10` | `36` vs `54` (fits) |
| 14 | `V2` | `243,466` | `2` | `4,593,685` | `48` | `10` | `38` vs `54` (fits) |
| 15 | `PL2` | `112,562` | `1` | `4,706,247` | `49` | `11` | `38` vs `53` (fits) |
| 16 | `PL3` | `745,873` | `7` | `5,452,120` | `56` | `11` | `45` vs `53` (fits) |
| 17–24 | `CG1..SYM1` (`8` stages) | `195,352` | `8` (`1` each) | `5,647,472` | **`64`** | `19` | `45` vs `45` (**at cap**) |
| 25 | `LRV1` | `2,535` | `1` | `5,650,007` | **`65`** | `20` | `45` vs `44` (**crosses 64**) |
| 26–27 | `LGRT1, MGE2` (`2` stages) | `22,725` | `2` (`1` each) | `5,672,732` | `67` | `22` | `45` vs `42` |
| 28 | `MV1` | `871,018` | `8` | `6,543,750` | **`75`** | `22` | **`53` vs `42`** |
| 29 | `V3` | `686,211` | `6` | `7,229,961` | **`81`** | `22` | **`59` vs `42`** |
| 30 | `LR1` | `856,624` | `7` | `8,086,585` | **`88`** | `22` | **`66` vs `42`** |
| 31 | `GART1` | `5,194,697` | **`43`** | `13,281,282` | **`131`** | `22` | **`109` vs `42`** |
| 32 | `GART2` (`r26.11.4` end) | `1,202,957` | `10` | `14,484,239` | **`141`** | `22` | **`119` vs `42`** |

**Key structural observation**: Because every small stage occupies a distinct Parquet file with `1` row group whose `[min, max]` covers nearly `[0x00.., 0xff..]` of the 16-byte BLAKE3 `contig_uid` space, the cumulative `64`-RG threshold is crossed at **`5.65M` contigs (`LRV1`, stage 25)** — well before `7.86M` — and `22` single-RG stages permanently pin `22` of the `64` LRU slots in any lookup that scans all stages.

---

## 2. Exhaustive Trace of All Lookup Paths

Every stage driver (`src/biobox/stage.py:1359`) and builder (`genecall.py:825, 1190`; `crispr.py:736`) orders `stages` with the current stage **first**: `stages = [stage, *(s for s in stages_on_disk(root) if s != stage)]`, where `stages_on_disk` (`store.py:138`) returns alphabetical order.

| Step / Component | Exact Code Location | Lookup Mechanism | Stage Order / Stop Condition | RGs Touched per Owned Contig | Thrashing Risk at `14.48M` (`32` stages)? |
|---|---|---|---|---|---|
| **`ingest`** (`store.ingest` & batch ingest) | `src/biobox/l0seq/store.py:696-724` | DuckDB `ANTI JOIN old_loc` / `old_idx` (`memory_limit="1GB"`, spills) | N/A (no `KeyLookup`) | `0` | **None** |
| **`genecall` plan** (`genecall.stage_items`) | `src/biobox/l1/genecall.py:1030-1082` | DuckDB `LEFT JOIN ci` (`memory_limit="2GB"`, spills) | N/A (no `KeyLookup`) | `0` | **None** |
| **`genecall` ownership** (`ContigOwners.owner`) | `src/biobox/l0seq/store.py:168-183`, `genecall.py:579` | One `KeyLookup(cache=64)` **per stage** in `self._index[s]` | `self.stages[0] == stage`; returns on first stage with `is_default_variant` | `1` RG of `stage` (`<= 43 <= 64`) | **None** (each stage has its own 64-RG LRU; max single stage = `43` RGs) |
| **`genecall` bases** (`_genome_bases` -> `lengths.get`) | **`src/biobox/l1/genecall.py:556-566, 868`** | Single `KeyLookup(store.index_files, ..., cache=64)` calling **`lengths.get(cu)`** | **BUG: scans ALL `S` stages (`lookup.py:107-115`), never stops early** | **`S` RGs (all `32` stages; `141` cumulative RGs > `64`)** | **CRITICAL — trips `280` CPU-min cap and `1.50×` `DEC-64` guard on `GART2` (and risks `280` cap on `GART1`)** |
| **`genecall` sequences** (`_sequences` -> `SeqStore.get`) | `src/biobox/l1/genecall.py:379-380`, `store.py:1641-1657, 1697-1718` | `store.index.find_first` + single `store.locator.find_first` (`B-66`) | Stops at `stages[0] == stage` for owned contigs | `1` index RG + `1` locator RG of `stage` (`<= 43 <= 64`) | **Low** (owned contigs 100% in-cache; minor spill on `GART1`'s `23,461` held contigs, §3.2) |
| **`crispr` ownership** (`ContigOwners.owner`) | `src/biobox/l1/crispr.py:672, 678`, `store.py:168-183` | Per-stage `KeyLookup(cache=64)` | Stops at `stage`; iterates in sorted `contig_uid` order (`_distinct(_range_rows(...))`) | `1` RG of `stage` sequentially | **None** |
| **`crispr` sequences** (`SeqStore.get`) | `src/biobox/l1/crispr.py:671, 679-686`, `store.py:1697-1718` | `find_first` + single `_where` (`B-66`); held contigs skipped at line `679-682` before `store.get` | **Only owned contigs** reach `store.get`; stops at `stage` (`<= 43` RGs) | `1` index RG (sequential) + `1` locator RG of `stage` (`<= 43 <= 64`) | **None** |
| **`release tapes`** (`tapes.run`) | `src/biobox/l3/tapes.py:208-221` | DuckDB single-stage sort/join (`2GB`, spills) | N/A | `0` | **None** |
| **`release tokens`** (`tokens.run`) | `src/biobox/l3/tokens.py:46-68` | DuckDB + hash-partitioned range routing (`B-46`) | N/A | `0` | **None** |
| **`release postings`** (`postings.run`) | `src/biobox/l3/postings.py:34-46, 78-81` | Arrow row-group streaming + routed `_family_run` + coalesced runs (`B-62`) | N/A | `0` | **None** |
| **`release dim`** (`dim.run`) | `src/biobox/dim.py:120-197` | DuckDB multi-stage join (`2GB`, spills); `stored_by` only on error (`dim.py:154`) | N/A | `0` | **None** |

---

## 3. Quantitative Impact Analysis (`*derived*` Formulas & `DEC-64` Guard Check)

### 3.1 Primary Defect: `genecall._genome_bases` (`src/biobox/l1/genecall.py:562`)

Two additive CPU penalties hit every owned contig in `_genome_bases`:

1. **In-cache `KeyLookup` stage-scan penalty**: From `B-66`'s measurement on `G1` (`2026-11-22_b66_seqstore_and_crispr_contigs.md` §2, `wilkcomp02`, *measured*), eliminating `34` extra in-cache `KeyLookup` file probes cut `SeqStore.get` by `0.9254 ms/contig`, i.e., **`27.2 us` per extra stage probe** (*derived*: `0.9254 ms / 34`). Scanning `S - 1` extra stages on `N` contigs costs:
   $$\text{CPU}_{\text{scan}} = N \times (S - 1) \times 27.2\ \mu\text{s}.$$
2. **LRU cache-miss row-group decode penalty**: Each `contig_index` row group holds `122,880` rows × `(16 + 8 + 1) = 25 B` = `3.072 MB` uncompressed (`~2.3–2.5 MB` zstd). Opening `pq.ParquetFile(path, metadata=md).read_row_group(g, columns=3)` + `combine_chunks()` (`lookup.py:94-96`) costs `t_rg ∈ [3, 10] ms` CPU (conservative lower `3 ms`, central `5 ms`, upper `10 ms`; *estimate* bounded by `C-42`'s `>3×` CPU inflation at `18%` miss rate in `2026-11-18_spire_pilot.md` §4.4 #2).
   - **1-contig genomes (`MV1`, `V3`, `LR1`, random `contig_uid` order)**: Using Che's independent-reference LRU approximation, the characteristic time $T$ (in calls) solves $\sum_{s \in \text{multi}} R_s (1 - e^{-T/R_s}) = 64 - 22 = 42$, and expected misses per contig call are $M = \sum_{s \in \text{multi}} e^{-T/R_s}$.
   - **Multi-contig MAGs (`GART1` `9,234` MAGs mean `565` contigs; `GART2` `15,320` MAGs mean `79` contigs)**: `iter_genome_input` (`inputs.py:139`) orders contigs ascending by `contig_uid` within each genome. One MAG of `m` contigs sweeps $K(m) = \sum_{s \in \text{multi}} R_s \bigl(1 - (1 - 1/R_s)^m\bigr)$ distinct multi-RG row groups monotonically from `0x00..` to `0xff..`. Because $K(m) \approx 109$ (`GART1`) and $111.5$ (`GART2`) both vastly exceed the `42` free LRU slots, **100% of multi-RG row groups are evicted between consecutive MAGs**, yielding $\text{decodes} \approx N_{\text{MAG}} \times K(m)$.

| Stage | Genomes / Contigs (*measured*) | `S` / `R_cum` | Miss Model (*derived*) | Total RG Decodes (*derived*) | Scan CPU-min (*derived*) | Decode CPU-min (`3 / 5 / 10 ms`, *derived*) | Total Overhead CPU-min (*derived*) | Baseline `genecall` Est. (batches, *measured* `A-26`) | Projected Batch CPU-min (`3 / 5 / 10 ms`) | `DEC-64` (`1.50×`) or `280` Cap Risk? |
|---|---|---:|---|---:|---:|---:|---:|---|---|---|
| `MV1` | `~750k` / `871k` owned | `28` / `75` | IRM $T=22.8$, $M=0.585$/contig | `~510k` | `10.7` | `25.5 / 42.5 / 85.0` | `36–96` | `~2,630` (`16` batches, `~164`/b) | `166 / 167 / 170` | Safe (`< 4%`) |
| `V3` | `686k` / `686k` (1-contig) | `29` / `81` | IRM $T=15.6$, $M=1.057$/contig | `~725k` | `8.7` | `36.3 / 60.4 / 120.8` | `45–130` | `3,286` (`14` batches, `235`/b) | `238 / 240 / 244` | Safe (`< 4%`) |
| `LR1` | `857k` / `857k` (1-contig) | `30` / `88` | IRM $T=11.2$, $M=1.766$/contig | `~1.51M` | `11.3` | `75.7 / 126.1 / 252.2` | `87–264` | `1,775` (`8` batches, `222`/b) | `233 / 239 / 255` | Recorded (`1.05–1.15×`), under `280` |
| **`GART1`** | `9,234` MAGs / `5.218M` (`m=565`) | `31` / `131` | MAG sweep $K=109.0$/MAG | **`~1.01M`** | **`70.4`** | `50.3 / 83.9 / 167.8` | **`121–238`** | **`641.4` (`3` batches, `213.8`/b)** | **`254 / 265 / 293`** | **HIGH RISK: breaches `280` cap at `>= 7.5 ms/RG`; `1.19–1.37×`** |
| **`GART2`** | `15,320` MAGs / `1.213M` (`m=79`) | `32` / `141` | MAG sweep $K=111.5$/MAG | **`~1.71M`** | **`16.9`** | **`85.4 / 142.3 / 284.7`** | **`102–302`** | **`192.6` (`1` batch, `192.6`/b)** | **`295 / 352 / 494`** | **CERTAIN STOP: `1.53×–2.57×` (breaches BOTH `280` cap and `1.50×` guard even at `3 ms/RG`)** |

*Notes*:
- `GART1` (`9,234` MAGs, `0` genomes `< 100 kb`) and `GART2` (`15,320` MAGs, `0` genomes `< 100 kb`) run `100%` in single mode (`genecall.py:283, 614`) via `ThreadPoolExecutor` (`genecall.py:669-674`), so `_genome_bases` runs **single-threaded in the main thread**, serializing the entire `1.01M–1.71M` row-group decode stream and starving worker threads.
- **RSS impact**: `lengths` exists only in the main process (`genecall.py:868`), holding at most `64` decoded row groups × `3.07 MB` ≈ `197 MB` Arrow buffers (`~250–350 MB` RSS). Worker subprocesses (used only for `< 10 kb` genomes in `MV1`/`V3`/`LR1`) do not create `lengths`, and their `SeqStore` caches stay bounded by `<= 8` stage RGs (`~100 MB/worker`). Thus **RSS is NOT an OOM risk** in `r26.11.4` (unlike `C-42`'s pre-`B-66` `>100 GB` across `35` workers).

### 3.2 Secondary Effect After Fixing `lengths.find_first`: Held Contigs in Mixed MAGs (`GART1` / `GART2`)

In `genecall`, a MAG that holds both owned and earlier-stage contigs is called whole (`genecall.py:41-44`), so `lengths.find_first`, `store.index.find_first`, and `store.locator.find_first` must fall past `stages[0]` for its held contigs. On `GART1`, `23,277` of `23,461` held contigs belong to `SPX1` (`26th` alphabetically out of `30` earlier stages, preceded by `25` stages / `49` RGs + `SPX1`'s `24` RGs). Because `43` (`GART1`) + `26` (one RG per probed stage) = `69 > 64`, each `SPX1`-held contig can evict a small number (`~1–5`) of `GART1`'s cached RGs. Across `3` batches (`~7,820` held contigs/batch), this adds at most `~1–5` CPU-min/batch (`< 2.5%` of `213.8` CPU-min/batch) — negligible and nowhere near `DEC-64` or `280` CPU-min. On `GART2`, `10` (`GART2`) + `26` = `36 <= 64`, so **zero** `GART2` RGs are evicted by `SPX1`-held contigs.

### 3.3 Exact Threshold Where `B-56` Actually Bites After Fixing `genecall.py:562`

Once `genecall._genome_bases` uses `lengths.find_first`, every owned-contig lookup in `genecall` (`ContigOwners`, `lengths`, `SeqStore.index`, `SeqStore.locator`) and `crispr` (`ContigOwners`, `SeqStore.index`, `SeqStore.locator`) stops at `stages[0] == stage`. The `64`-RG LRU cliff therefore shifts from **cumulative across the root** to **per single stage**:
- **Exact threshold**: **`> 64 × 122,880 = 7,864,320` stored contigs in ONE stage**.
- **Headroom in `r26.11.4`**: largest single stage is `GART1` at `5,194,697` new contigs = **`43` RGs** (`67.2%` of the `64`-RG cap; `21` RGs / `2.67M` contigs of headroom).
- **When `B-56` must be revisited**: before any future stage (e.g., an unsplit `SPIRE`, `IMG-UMAG`, or `GlobDB` tranche) exceeds `7.86M` stored contigs in a single stage (or keep every stage `<= 6M` contigs / `<= 49` RGs as `C-42` did).

---

## 4. Ready-to-Launch Minimal Fix Card (`B-70` / `B-56a`)

```markdown
### B-70 Use `KeyLookup.find_first` in `genecall._genome_bases`
**L** · Pod: no (pure R1 unit test; run `pod_pytest.sh` per AGENTS §5.5) · Depends: B-66 · uncovered 2026-11-24 by `2026-11-24_b56_scale_audit.md` before `r26.11.4` (`GART1`/`GART2`)
Files: `src/biobox/l1/genecall.py` (lines 556-566 only), `tests/r1/test_b66_seqstore_and_crispr_contigs.py` (or new `tests/r1/test_b70_genecall_genome_bases.py`)
Why: Card B-66 switched `SeqStore.default_digest` and `SeqStore._where` to `KeyLookup.find_first`, but left `genecall._genome_bases` (`src/biobox/l1/genecall.py:562`) calling `lengths.get(cu)`. On `r26.11.4` (`31–32` stages, `131–141` cumulative `contig_index` row groups vs `64`-RG LRU), every owned contig scans all `31–32` stages, causing `~1.01M` RG decodes on `GART1` (`+121–238` CPU-min, risking the `280` CPU-min cap) and `~1.71M` RG decodes on `GART2` (`+102–302` CPU-min on a `192.6` CPU-min single batch -> `1.53×–2.57×`, guaranteed `DEC-64` and `280` CPU-min stop).
Do, test first:
1. Test first: on a synthetic multi-stage store (`>= 10` stages, `stages[0]` owning the queried contig), call `genecall._genome_bases(store, lengths, genome)` and assert `lengths.row_groups_read == 1` (fails on current code which reads all `10` stages' row groups), while still returning the exact base count for contigs owned by an earlier stage (`stages[-1]`) and raising `KeyError` naming the owning stage on disk when absent.
2. In `src/biobox/l1/genecall.py:562-565`, replace:
   `rows = [r for r in lengths.get(cu) if r["is_default_variant"]]`
   `if not rows: ...`
   `n += rows[0]["contig_len"]`
   with:
   `row = lengths.find_first(cu, lambda r: r["is_default_variant"])`
   `if row is None: ...`
   `n += row["contig_len"]`
3. Verify bit-for-bit identical `genecall.run` output (feature, feature_protein, protein_occurrence, protein_seq, ledger) and `micro_s0` manifest hash.
Accept: new test fails before edit (`row_groups_read == 10`) and passes after (`row_groups_read == 1`); full pod pytest suite passes; zero byte-identity drift.
```
