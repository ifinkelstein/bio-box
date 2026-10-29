# Audit Report 7: Database, Storage-Engine & Query-Path Audit of `all-the-data` (v2.0)

**Auditor:** Database / storage-engine reviewer (columnar formats, inverted indexes, compressed bitmaps, DuckDB/Parquet internals, LSM stores, query planning).
**Targets:** `00-spec.md`, `PLAN.md` (v2.0), `DEPLOYMENT-AND-MALLEABILITY.md`, `PRIOR-ART.md`, `audits/critique_2_systems.md`, `audits/SYNTHESIS.md`.
**Verdict:** The two-projection idea (contig-centric tapes plus anchor-centric postings) and the use of compressed bitmaps for co-occurrence are sound. The numbers attached to them are not. The plan has no single cardinality model, so row counts, byte sizes and latencies disagree with each other by 2–10×. Three index structures (16-slot windows, `PackedSlotDescriptor`, Stage-0 dedup) return wrong answers as specified, not just slow ones. Several latency guarantees assume DuckDB/Parquet behaviours that DuckDB does not document. None of this sinks the project; all of it has to be fixed before code is written, and the fix for most items is "measure at 1% scale against a brute-force oracle".

---

## Step A: the deeper version of this task

**(1) What a top professional would do that the request does not ask for.**
Write the workload first (query classes, output sizes, concurrency, p50/p99 warm and cold), then derive the physical design from it. Keep one sizing ledger in which every count (contigs, slots, postings, pairs, tokens per level) is derived from per-record measurements, so totals cannot drift between documents. Insist that every index has a correctness oracle (brute-force answer on a sample) before its speed is discussed. Treat identifier stability, schema evolution, concurrent access, determinism of the build and crash recovery as first-class design requirements, not deployment details.

**(2) The rawest material to start from.**
Not the plan's derived totals. Start from: a few hundred real Logan contig files and an IMG/VR and GTDB-rep sample; the actual `pyrodigal`/`pyhmmer` output on them (slots per contig, domains per ORF, intergenic gap distribution); the real token-frequency distribution per level (Pfam/COG hit counts are public via AnnoTree/GTDB annotations); the source and docs of the libraries (DuckDB Parquet reader, CRoaring, RocksDB); and `fio` numbers from the actual target machine.

**(3) Methods and tools.**
Explicit back-of-envelope arithmetic with named cardinalities; Zipf/skew modelling of token frequency; CRoaring `roaring_bitmap_statistics` for container-type mix; DuckDB `EXPLAIN ANALYZE` and `parquet_metadata()`; `fio`, `perf`, `bpftrace`, `vmtouch` for I/O and page-cache behaviour; RocksDB `db_bench`; external sort-merge group-by instead of streaming hash state; mergeable sketches (bitmasks, histograms, HLL) for delta-epoch aggregates; property-based tests (`proptest`/`hypothesis`) on bit-packed encodings; versioned-catalog patterns (Iceberg/DuckLake-style manifests, atomic rename).

**(4) What other than you or me would show the result is right.**
Exact agreement between each index and a brute-force scan on a 1% sample (zero missed and zero extra hits); measured p50/p99 on hardware whose baseline has been characterised with `fio`; bit-identical rebuilds (same inputs give the same BLAKE3 manifest); the libraries' own documentation and source; and a third party re-running the harness and getting the same numbers.

**Rewritten prompt.**
> Audit the storage and query layer of `all-the-data` as a database engineer. (a) Build one sizing ledger: derive contigs, slots, windows, postings per token level, pair rows, unique pairs and tokens per level from stated per-record assumptions, and flag every place in the five documents that contradicts it. (b) For each index (Index 1 windows plus descriptors, `contig_tapes`, `window_postings`, `directon_pair_summary`, Index 3 bitvectors, Stage-0 dedup), state whether it returns *correct* answers as specified, with a concrete counterexample if not. (c) Check every latency claim against what DuckDB, Parquet, CRoaring and RocksDB actually do, citing docs, and recompute it with explicit I/O and decode arithmetic, separating warm from cold and bounding by output size. (d) Audit bit widths, integer ranges, ID namespaces and ID stability across epochs and releases. (e) Audit concurrency and malleability semantics (who writes the catalog, which aggregates survive delta overlays, what a remap or compaction does to keys). (f) Propose a workload model with SLOs and a 1%-scale benchmark harness whose pass criterion is exact agreement with a brute-force oracle. Show all arithmetic. Rank findings by severity with a concrete fix and a validation test for each.

---

## Step B: findings

Severity: **S1** returns wrong results or cannot be built as written. **S2** misses a stated SLO or size by >3×. **S3** correctness edge case or moderate cost. **S4** hygiene.

### DB-01 (S1) There is no consistent cardinality model
*Where:* `00-spec.md` §1, §2.6; `PLAN.md` §1.2, §2.2 table, §4.1, §4.2; `critique_2` §2–3; `DEPLOYMENT` §4.1.

Contradictions:
- Contigs: `PLAN` says 2.5×10⁹ post-sieve contigs; `critique_2` says 5×10⁸. With 3.6×10⁹ slots, `PLAN` implies **1.44 slots per contig**. The Stage-0 gate rejects contigs with <2 features unless the longest ORF is ≥250 aa, so the mean must be ≥~2 slots per contig. 2.5×10⁹ contigs would need ≥5×10⁹ slots. `critique_2`'s 5×10⁸ gives 7.2 slots per contig, which fits the gate.
- Tier 1: "1.8M assemblies" yields "~480M dereplicated slots". 1.8×10⁶ genomes × ~3,500 CDS × ~1.4 segments ≈ **8.8×10⁹ slots**. 480M is roughly the 113,104 species reps alone (113,104 × 3,500 × 1.2 ≈ 4.8×10⁸). Yet the `core` edition advertises "full Tier 1 (1.8M HQ genomes)" in 420 GB, and Index 3's `pan_any_bits` needs presence across all 596,859 GTDB genomes. Either the tapes hold species reps only, and the pangenome planes come from a separate per-genome presence table that is not in the schema, or the slot count is ~18× too low.
- Bytes per slot: the `contig_tapes` schema in `PLAN` §4.2 adds up to about **81 B/slot uncompressed**. Summing: regimes 1, vehicle 4, orf 2, seg 1, parent 1, topo flags 1, kind 1, strand 1, start/end 8, env start/end 4, delta 2, junction 1, directon 2, regulon 2, edges 4, clan 4, family 4, subfamily 4, motif 4, intact 1, fold 4, super-cluster 4, badge 8, nr90 8. 3.6×10⁹ × 81 B = **292 GB raw**, which compresses to roughly 100–150 GB. `PLAN` claims 1.05 TB compressed, and the §2.2 tier table claims 2.4 TB for the same 3.6×10⁹ slots (667 B/slot). Compressed output cannot be 3.6–8× larger than raw input. So either slot counts are low by about 5–10× or sizes are high by that much.
- Windows: `PLAN` §4.1 gives `window_id = global_slot_idx >> 4 ∈ [0, 2.8×10⁹)`. But 3.6×10⁹ >> 4 = **2.25×10⁸**. 2.8×10⁹ is only reachable through the separate rule that every contig starts a new window (≥1 window per contig × 2.5×10⁹ contigs). At 1.44 slots per contig, windows are ~9% full and the "16 slots ≈ 10–12 ORFs" description is meaningless. At 5×10⁸ contigs, windows ≈ Σ⌈nᵢ/16⌉ ≈ 5×10⁸ + ~3×10⁷ (long Tier-1 chromosomes) ≈ **5.3×10⁸**.

*Fix:* Add a `SIZING.md` ledger (or a table in `PLAN` §1.2) with one row per entity: per-record bytes (raw and measured compressed), count derivation, and total. Make every other number in the suite cite it. Seed it from the 1% build (see harness below), not from estimates.
*Validation:* The harness prints the ledger from the measured 1% build × 100, alongside the plan's values. CI fails if any document number differs from the ledger by >20%.

### DB-02 (S1) 16-slot aligned windows silently drop 22% of true neighbour pairs
*Where:* `00-spec.md` §2.6; `PLAN` §4.1 Layer 1A; `critique_2` §2.2.

`roaring_and(win(A), win(B))` only proves co-occurrence when A and B share a window. With aligned, non-overlapping 16-slot windows, a pair at slot distance d straddles a boundary with probability d/16. Example: slots 15 and 16 sit in different windows. At d = 6 the loss is 37.5%; averaged uniformly over d = 1..6 it is (1+…+6)/(6·16) = 3.5/16 = **21.9% false negatives**. Separately, the window is in *slots* while the query predicate is in *ORFs* (±6 ORFs). A neighbourhood containing one 10-domain NRPS/PKS spans >16 slots inside 6 ORFs and is missed entirely. Repeat collapse (BIO-04) does not cover modular megasynthases.

The formula `window_id = global_slot_idx >> 4` also contradicts the rule that contigs start on a window boundary. With padding, the ID must be `contig_window_base[c] + (slot_in_contig >> 4)`.

*Fix (pick one, then benchmark):*
1. Index at slot (or ORF) granularity: `token_to_slots` bitmaps over a 32-bit slot ID space, partitioned per tier/epoch (see DB-08). Proximity becomes `A ∧ (⋁_{δ=-k..k} B ≪ δ)` using `roaring_bitmap_add_offset`. This is exact, and a contig-boundary check comes from a `contig_start_slots` bitmap.
2. Keep windows but make them overlapping (width 16, stride 8, so each slot is in 2 windows; 2× postings), and define windows on `orf_idx`, not slot index.
3. Keep aligned windows but intersect `win(A)` with `win(B) ∪ (win(B)+1) ∪ (win(B)−1)` (about 3× bitmap work), then verify in Layer 2.

*Validation:* On the 1% build, enumerate every (A, B, d ≤ 6) pair by brute-force scan of the unnested tape and require **exact set equality** with Index 1 output, across 10⁴ random queries stratified by token frequency.

### DB-03 (S1) `PackedSlotDescriptor` cannot answer the question it exists for
*Where:* `PLAN` §4.1 Layer 2; `critique_2` §2.2 (different layout).

`PLAN` layout: bits 0–11 `contig_offset_in_block` (4,096 contigs per block), bits 12–14 `orf_mod_8`, bit 15 strand. Problems:
- No slot or ORF position within the contig is stored, only ORF mod 8. Inside a window the true ORF distance ranges 0..15, so mod 8 aliases d with d ± 8. A stored difference of 7 can mean −1 or +7, so neither order nor the ±6 test can be decided. Example: A at ORF 3 and B at ORF 10 (d = 7, outside ±6) gives (10−3) mod 8 = 7. A at ORF 3 and B at ORF 2 (d = −1, inside ±6) also gives (2−3) mod 8 = 7.
- `contig_offset_in_block` is redundant with the window ID, which already determines the contig.
- The slices are keyed by `(token_id, block_id = dense_contig_id >> 12)`, while Roaring containers are keyed by `window_id >> 16`. The two arrays cannot be aligned by container rank.
- When a token occurs twice in one window, the bitmap holds the window once but there are two descriptors. A variable-length mapping with offsets is needed and not specified. `critique_2` handles this with a 2-bit "overflow index" pointing to an unspecified list.
- `PLAN` and `critique_2` specify incompatible bit layouts. Neither says which is authoritative.

*Fix:* If DB-02 option 1 is adopted, the descriptor disappears (position is the bitmap value) and per-slot attributes are read from a fixed-stride slot attribute array (`slot_attrs.u32[global_slot]`: strand, kind, junction, directon-local ID, edge flag), mmapped at 3.6×10⁹ × 4 B = 14.4 GB. If windows are kept: per posting `uint16 = slot_in_window(4) | orf_delta_from_window_first_orf(4) | strand(1) | junction_prev(3) | regulon_local(3) | edge(1)`, stored in a per-container posting array with a `uint16` prefix-count array for repeated tokens. Key slices by the same `window_id >> 16` as the Roaring container.
*Validation:* Property tests over random contig layouts (random domain counts, strands, repeats, contig edges): encode, decode, and check that the distance predicate matches brute force.

### DB-04 (S1) Stage-0 dedup is logically broken and nondeterministic
*Where:* `00-spec.md` §2.2; `PLAN` §2.3 Filter 4; `critique_2` §4.2.

1. **First sightings never reach RocksDB.** On "NOT PRESENT" the record is only added to the Cuckoo filter and appended to `tier3_unique_contigs.bin.zst`. On the second sighting, RocksDB is looked up for `(hash, biome_bucket)` and holds nothing, so there is no representative sketch to compare against. The branch is undefined.
2. **Key mismatch.** The Cuckoo filter is keyed by `stage0_arch_hash`, RocksDB by `(hash, biome_bucket)`. A hash first seen in biome X and then in biome Y passes the filter and misses RocksDB.
3. **One representative per key.** If the ANI check fails (same architecture, divergent sequence), the value can hold only one rep, and what happens to the second contig is not specified.
4. **Order dependence.** The representative is whichever contig the S3 listing order delivers first. Builds are not reproducible, and the rep is not the best one (longest, circular, highest-quality sample).
5. **Cuckoo filter is mis-sized by 8×.** 10⁹ × 16-bit fingerprints = 2.0 GB at 100% load, 2.1 GB at the usual 95% load. 16 GB holds 8×10⁹ fingerprints. False-positive rate with 4-slot buckets is about 2b/2^f = 8/65,536 = 1.2×10⁻⁴, about 1.2×10⁵ false "present" results over 10⁹ inserts, each of which then hits bug 1. Cuckoo filters also cannot be resized when the 6×10⁸ estimate turns out low.
6. **The ANI sketch threshold is wrong by ~4 points of ANI.** For a 64-register 1-bit (b = 1) one-permutation MinHash, P(bit equal) ≈ ½ + ½J. Hamming ≤ 6/64 means P_eq ≥ 0.906, so J ≥ 0.81. Under the Mash relation D = −(1/k)·ln(2J/(1+J)) with k = 15, that gives D = 0.0074, i.e. **≈99.3% ANI, not 95%**. 95% ANI (D = 0.05) corresponds to J ≈ 0.31, P_eq ≈ 0.65, Hamming ≈ 22. The estimator's standard error at 64 bits is √(0.65·0.35/64) ≈ 0.06 on P_eq, ≈0.12 on J, too noisy to threshold. Jaccard is also the wrong measure for fragments of unequal length; a 2 kb fragment inside a 20 kb contig needs *containment*.
7. **Irreversible loss.** Collapsed contigs keep only a `(rep_contig_id, sample_id, biome_id)` posting. If Stage 2 shows two collapsed contigs differ (BIO-01/BIO-11), the stage-2 "canonicalization" can only merge, never split, and the 40.8 TB stream would have to be re-run.

*Fix:* Replace the streaming probabilistic design with a deterministic **two-pass external sort-merge**, which is the textbook DB approach and fits in the same RAM:
- Pass 1 (streaming, stateless): for each contig that passes Filters 1–3, emit `(stage0_arch_hash, biome_bucket, quality_rank, contig_ref, sketch[ ≥1,024-bit FracMinHash or 128 × 16-bit bottom-k ])` into 256 hash partitions on NVMe.
- Pass 2: sort each partition (RAM-sized), group by key, and within each group run greedy clustering by containment ≥ threshold, using `quality_rank` to choose reps deterministically. Write `contig_membership.parquet (member_ref → rep_ref, containment, ani_est)` so dedup can be re-done at a different threshold without re-streaming.
- Drop the Cuckoo filter. If an online variant is kept, insert every first sighting into RocksDB (5×10⁸ × ~40 B ≈ 20 GB SST, bloom 10 bits/key ≈ 0.6 GB) and let RocksDB's bloom filters do the negative lookups.

*Validation:* Same input in shuffled order gives an identical `contig_membership` BLAKE3 hash. On 10⁴ sampled pairs, sketch-estimated ANI is checked against `skani`/`fastANI`, requiring |error| < 1 ANI point at 95%.

### DB-05 (S2) `contig_tapes.parquet` cannot serve 1,000 random contigs in <12 ms
*Where:* `00-spec.md` §2.6; `PLAN` §1.2, §4.2 "Retrieval Guarantee"; `critique_2` §3.2.

- DuckDB's `rowid` pseudo-column belongs to native tables. For Parquet the equivalent is the virtual `file_row_number` column, which DuckDB documents as an output column. I found no documentation that a predicate on it is pushed into the scan ([DuckDB Parquet docs](https://duckdb.org/docs/2.0/data/parquet/overview.html)). Pushdown is documented via zone maps (row-group min/max statistics). The Parquet project lists DuckDB's page-index support as read-only ([Parquet implementation status](https://parquet.apache.org/docs/file-format/implementationstatus/)). The minimum decode unit is therefore a page per column, and without a page index it is a row group per column.
- Arithmetic: 1,000 contigs selected by a biological query are scattered across ~1,000 distinct row groups. Each needed column (≈25 LIST columns) costs at least one page per row group. With 64 KB pages that is 25 × 1,000 × 64 KB = **1.6 GB** read and decompressed: ~115 ms at 14 GB/s, ~25,000 random reads ≈ 12.5 ms at 2 M IOPS (only at very high queue depth), plus zstd decode. With DuckDB defaults (122,880-row row groups and larger pages) it is worse by 10× or more. Realistic warm p50: **50–300 ms**; cold: seconds.
- The sort key `(tier, tax_order_id, stage2_arch_hash)` (`PLAN`) disagrees with `(tier, tax_order_id, biome_id, dense_contig_id)` (`critique_2`).

*Fix:* Serve neighbourhood fetches from a purpose-built **offset-indexed tape store**: `tape_offsets.u64[dense_contig_id]` (5×10⁸ × 8 B = 4 GB, mmapped) pointing into a row-oriented, per-contig, block-compressed heap (~24–32 B/slot → ~100 GB). 1,000 contigs ≈ 1,000 × 4 KB reads ≈ 4 MB, **<2 ms warm, ~5–10 ms cold** at queue depth 64. Keep `contig_tapes.parquet` as the analytic and export projection. If a single format is wanted, benchmark Lance (designed for random row access) against the custom heap at 1% scale.
*Validation:* Harness Q2 (below): p50/p99 for 1,000 random and 1,000 query-derived contigs, cold (after `echo 3 > /proc/sys/vm/drop_caches`) and warm, each mode measured separately.

### DB-06 (S2) `window_postings` size is a 10× arithmetic error, and hub anchors are 100× bigger than claimed
*Where:* `PLAN` §1.2, §4.2 Projection 2B; `critique_2` §3.2; `DEPLOYMENT` §3.

- 1.8×10¹⁰ rows × 7.5 B = 1.35×10¹¹ B = **135 GB, not 1.35 TB**.
- 7.5 B/row is not credible either. The row has ~47 B raw: anchor 4, anchor_subfam 4, partner 4, partner_subfam 4, fold 4, offset 1, topology 1, flags 1, delta 2, tax 2, vehicle 4, badge 8, rep contig 4, rep slot 2, count 2. The two leading sort keys compress to ≈0. The remaining ~39 B are mostly high-entropy (badge, vehicle, rep IDs, subfamilies), so a realistic figure is **15–25 B/row**, giving 2×10¹⁰ × 20 B ≈ **0.3–0.5 TB**.
- Row count: with 7.2 slots per contig, each slot sees ~6 partners after edge truncation, so 3.6×10⁹ × 6 ≈ 2.2×10¹⁰ at one token level. Keying `anchor_token_id` at more than one level (family and subfamily, as the queries need) multiplies this.
- **Hub skew.** For `PF00005` (ABC_tran) the claim is "1–4 row groups, 5–80 MB, 3–10 ms". Estimate ≥5×10⁷ occurrences × 6–12 partners = 3–6×10⁸ rows × 20 B = **6–12 GB for that one anchor**. That is ~0.5–1 s of I/O plus decode, and the specific pair `(ABC_tran, ABC_membrane)` alone holds tens of millions of rows.
- Footer cost: about 2×10¹⁰ / 122,880 ≈ 1.6×10⁵ row groups × ~15 column stats. A monolithic footer runs to tens or hundreds of MB and is parsed per query unless the daemon caches it.

*Fix:* Correct the size to a ledger-derived value. Route anchor-profile queries for any anchor with df > 10⁵ to `directon_pair_summary` (the precomputed answer), and use `window_postings` only for (anchor, partner) drill-down with `LIMIT`/stratified sampling by `tax_order_id`. Partition into files by `anchor_token_id` bucket (e.g. 4,096 files) so a query opens one small footer. Keep an `(anchor, partner) → (file, first_row_group, n_rows)` offset table in the daemon.
*Validation:* Harness Q3/Q4 split by df stratum: top-100 hubs, df 10⁴–10⁶, tail.

### DB-07 (S2) The hub definition is arithmetically impossible, and hub latency ignores output size
*Where:* `PLAN` §4.1 Layer 1B; `critique_2` §2.2.

"Top H = 2,048 hub tokens with >2×10⁶ occurrences" requires ≥2,048 × 2×10⁶ = 4.1×10⁹ postings, which exceeds the 3.6×10⁹ total slots. At most 1,800 tokens per level can clear 2×10⁶, and under any Zipf-like distribution it will be a few hundred. Separately, materialising `hub_pair_windows[(ABC_tran, HTH)]` in 0.05 ms does not make the *query* fast: that bitmap may hold ~10⁷ windows, and enumerating, verifying and returning them is output-bound (10⁷ × ~10 ns ≈ 100 ms just to iterate, before fetching any tape).

*Fix:* Define hubs by a df threshold measured on the 1% build (×100). State SLOs separately for `COUNT` (bitmap cardinality, O(#containers), sub-ms), top-k (sampled), and full enumeration (output-bound, streamed with a cursor).

### DB-08 (S2) Index 1 RAM is underestimated, and 32-bit headroom is thin
*Where:* `PLAN` §1.2, §4.1; `critique_2` §2.2.

- Roaring array containers cost 2 B per value; only bitmap containers (density > 4,096/65,536 = 6.25%) and run containers go below that. After architecture dedup, adjacent contigs carry *different* architectures, so runs are rare. Hub density example: 2×10⁶ / 5.3×10⁸ windows ≈ 0.4%, i.e. ~246 values per container, all array containers at 2 B each. So almost every posting costs ~2 B. `critique_2`'s 1.15 B/posting and `PLAN`'s 3.2 GB for ≥3 token levels (≈0.3 B/posting) are not reachable.
- Estimate: ~3.6×10⁹ postings × 3 levels (family, subfamily, fold; clan is cheap) ≈ 1.1×10¹⁰ × 2 B ≈ **22 GB**, plus 1 GB hub bigrams and 10–14 GB position data. Index 1 is **~35 GB, not 14.2 GB**. That still fits in 256 GB, so this changes the RAM budget, not feasibility.
- Headroom: `global_slot_idx` at 3.6×10⁹ uses 84% of 2³² (4.29×10⁹). If the slot count is in fact higher (DB-01), it overflows. CRoaring has a native 64-bit API (`roaring64_bitmap_*`; the Rust `croaring` crate exposes `Bitmap64`/`Treemap`, [docs.rs/croaring](https://docs.rs/croaring), [croaring-sys](https://docs.rs/croaring-sys)). The 64-bit variants are slower and their serialization is less standardised.

*Fix:* Use per-(tier, epoch) 32-bit bitmaps over local slot IDs (each ≤2³²), evaluate queries per partition and union the results. Size from `roaring_bitmap_statistics` on the 1% build.

### DB-09 (S2) Index 3 dense bitvectors do not scale to subfamily tokens
*Where:* `PLAN` §4.4; `DEPLOYMENT` §4.2 (`pangenome_3plane.bits` 4.5 GB).

Each token costs 3 planes × 14,144 B = 42,432 B. 4.5 GB therefore covers only **~106,000 tokens**. `PLAN` indexes Index 3 by `subfamily_50_id` *and* `seq_family_id`; subfamilies number in the millions (from ~85M NR50 reps): 10⁶ → 42 GB, 10⁷ → 424 GB. The `lite` edition (85 GB) cannot contain them. Most subfamilies are present in few species, so dense rows are mostly zeros.

*Fix:* Store planes as Roaring bitmaps over species leaf IDs (113,104 < 2¹⁷, i.e. two containers); keep dense padded rows only for the top ~10⁵ tokens by species count, where AVX-512 popcount helps. `transition_edges` at <400 B/token × 10⁷ = 4 GB is fine. Also check that GTDB r220 order count < 2,048 before fixing the order bitmask width in `directon_pair_summary`; otherwise use a Roaring bitmap.

### DB-10 (S2) The "O(1)-memory, 4.5-minute" pair summary covers only the easy columns
*Where:* `PLAN` §4.3 Table 2C; `critique_2` §3.3; `DEPLOYMENT` §3.1 (delta overlay).

A sorted streaming group-by gives O(1) state per group for counts, bitmasks and histograms; that part is correct. I/O check: 0.4 TB at 14 GB/s ≈ 30 s, so minutes is plausible. But the table also promises:
- `n_veh_indep_events`: single-linkage over ±15-ORF flank Jaccard between all loci of a pair. That is O(n²) per pair: a hub pair with n = 10⁷ loci means 5×10¹³ comparisons.
- `gene_tree_homoplasy_b_on_a`: needs a phylogeny of family A, with families of up to 10⁶–10⁷ members.
- `empirical_fdr_qval`: needs decoy datasets processed through the same pipeline.

None of these is computable in a single streaming pass. Several aggregates are also **not mergeable**: median, IQR (merge the histogram instead), entropy (needs the neighbour histogram), LLR and q-values (global). The `DEPLOYMENT` §3.1 overlay that "adds delta pair counts on the fly" is wrong for every non-additive column.

*Fix:* Store sufficient statistics (count vectors, 32-bin Δbp histograms, phylum/order bitmasks, HLL for genera, top-K neighbour histograms) and compute derived statistics at read time or at compaction. Specify the vehicle algorithm as MinHash-LSH over flank sets (O(n) per pair with banding) and give its budget. Mark LLR and q columns `stale_since_epoch` until compaction.

### DB-11 (S2) Malleability semantics break keys and sort orders
*Where:* `00-spec.md` §2.12; `DEPLOYMENT` §3.1–3.4; `PLAN` §6.2.

- **Renumbering.** `bio-box compact` "merges delta epochs into a new base". If the base is re-sorted by `(tier, tax_order, arch_hash)`, every `dense_contig_id` and slot index changes, and sidecars keyed by `(dense_contig_id, slot_idx)` silently point at wrong rows.
- **Unstable cluster IDs.** `subfamily_50_id` and `nr50_domain_id` come from MMseqs2 clustering, which changes when re-run on new input. Sidecars keyed by them do not survive a release.
- **Remap.** A token *merge* can be done by OR-ing bitmaps and scanning several anchor ranges. A *split* "based on a new active-site alignment" needs per-slot reassignment by a criterion not stored on the slot. That is a 3.6×10⁹-row slot sidecar, not a "<1 s overlay". Remapped tokens also break `window_postings`' sort order: the new token's rows are scattered across old anchor ranges.
- **Delta downloads.** Content-addressed 64 MB chunks only help when files are byte-identical. After a rebuild that renumbers IDs nearly every chunk changes, so "download only changed chunks" ≈ a full download.

*Fix:* (a) Give every contig a stable external key (xxh128 of sequence plus source accession) and every cluster a stable key (hash of the representative sequence), with per-release crosswalk tables (`id_map_vN_to_vN+1.parquet`). Dense IDs become release-internal and sidecars must use stable keys. (b) Never renumber inside a release; compaction writes a new release. (c) Document remap as merge-only; a split is a slot-level sidecar plus a planner rewrite. (d) Order files so stable content stays in stable chunks (per-tier, per-epoch files).

### DB-12 (S2) Multi-process writes to `catalog.duckdb` are not supported
*Where:* `DEPLOYMENT` §2.3, §3.1–3.2, §4.2 step 3.

`serve` attaches `catalog.duckdb` read-write while `append`, `register-sidecar` and `remap-token` run as separate `docker run` processes that also write it. DuckDB allows either one read-write process or many read-only processes; writing the native format from multiple processes is not supported, and DuckLake with a PostgreSQL catalog is the documented route for multi-writer setups ([DuckDB concurrency](https://duckdb.org/docs/current/connect/concurrency.html)). The second writer will fail to open the file, or, if the daemon opened it read-only, `register-sidecar` cannot update views.

*Fix:* The daemon owns all writes; the CLI subcommands call the daemon API. Alternatively, make `catalog.duckdb` a derived, read-only artifact regenerated from a small `catalog.json` manifest committed by atomic rename, and have each process build its views in memory at startup. Readers open with `access_mode = READ_ONLY`.

### DB-13 (S3) `catalytic_motif_hash` cannot support Hamming distance
*Where:* `00-spec.md` §2.10; `PLAN` §4.3 Table 2D (`exp(−1.4·d_Hamming)`); `critique_3` §3.2.

A hash destroys distance: two motifs differing in one residue give unrelated `xxhash32` values, so C_iso is not computable from stored data. xxhash32 also collides: with 10⁶ distinct motifs the birthday estimate is N²/2³³ ≈ 116 collisions.

*Fix:* Store `catalytic_motif_packed: uint64` as up to 12 residues × 5 bits = 60 bits, plus a 4-bit length, scoped by `seq_family_id` (columns are family-specific; motifs from different families are not comparable). Hamming distance is a per-lane XOR plus popcount on 5-bit lanes. Keep the hash only as a convenience join key, if at all.

### DB-14 (S3) `biophysical_badge` bit layout overflows
*Where:* `PLAN` §3.3; `critique_1` BIO-06 fix (bits 40/41).

Bits 9–14 (6 bits) hold `sorting_tag`, but 8 flags are listed: PEP-CTERM, LPXTG, GlyGly-CTERM, PGF-CTERM, Archaeosortase PGF, BMC-EP, Encapsulin CLP, RiPP leader. `critique_1` assigns bits 40/41 to frameshift and intron flags, which collide with `repeat_unit_count` (bits 40–47). Other fields fit: `pi_quantized` 5 bits for 30 levels; `tm_count` 5 bits.

*Fix:* Re-lay out the badge: tm 0–4, signal 5–8, sorting 9–16 (8 bits), pI 17–21, motif 22–29, length 30–41, repeat 42–49, edit/flags 50–53, reserved 54–63. Add a `badge_version` in the file metadata, and macros `badge_tm(x)` etc. in `catalog.duckdb`.
*Validation:* Round-trip property test over all fields at their maxima.

### DB-15 (S3) Integer widths overflow on real biology
*Where:* `PLAN` §4.2 tables.

- `delta_bp_prev_slots: int16` (±32,767). Gaps >32 kb occur after large CRISPR arrays, unannotated islands, and stitched intron spans. In `window_postings.delta_bp: int16`, ±6 ORFs routinely exceeds 32 kb in BGCs: six 5–10 kb NRPS/PKS genes give 30–60 kb. Fix: `int32`, or a saturating int16 plus an overflow flag.
- `parent_seg_idxs_in_orf: int8` (max 127) versus `seg_idxs_in_orf: uint8` (max 255). Use `int16` or the same domain.
- `global_slot_idx` in uint32: 84% used (DB-08).
- `dense_contig_id: uint32` with 2.5×10⁹ uses 58% of the range; fine at 5×10⁸.
- `orf_idx: uint16` is fine (the largest bacterial chromosomes have ~13k genes). `n_slots: uint16` is fine.

### DB-16 (S3) Token namespace and many-to-many annotation are unspecified
*Where:* `PLAN` §4.1 Layer 1A ("maps each `token_id` (`seq_family_id`, `subfamily_50_id`, or `struct_fold_id`)"), §4.2.

`token_to_windows.roar` mixes levels and sources (Pfam, COG, NCBIfam, Rfam, DARK, ECOD) in one ID space, but each level is a separate `uint32` column with no stated namespacing, so clan 5, family 5 and fold 5 would collide. A domain segment usually hits Pfam, COG and NCBIfam at once, yet the tape has one `seq_family_id` per slot, and the precedence rule is not stated.

*Fix:* A global `token_dictionary` with `token_id = level(3 bits) | source(5 bits) | local(24 bits)`, or 64-bit if counts require it. A side table `slot_annotations(slot_key, token_id, evalue, env)` holds the many-to-many hits; the tape column is the canonical pick under a documented precedence.

### DB-17 (S3) Cold-start and "zero-copy mmap" claims overstate what happens
*Where:* `00-spec.md` §2.12; `DEPLOYMENT` §1, §4.2 step 3.

- `MADV_WILLNEED` is an asynchronous readahead hint, not a residency guarantee. 8.7 GB at 14 GB/s = 0.62 s only if readahead saturates both drives; a single-threaded page-cache fill typically runs at 2–4 GB/s, about 2–4 s. Everything outside those 8.7 GB (tape store, postings, summary) is cold, so p99 for first queries is NVMe-bound.
- Parquet is never "zero-copy": DuckDB reads pages with I/O into its buffer manager and decompresses them. Roaring can be zero-copy only through CRoaring's frozen format (`roaring_bitmap_frozen_view`, with alignment requirements). The portable format needs deserialization.
- Download times: 3.45 TB at 10 Gbps = 3.45×10¹²×8/10¹⁰ = **2,760 s = 46 min** at line rate, so 52 min requires ~89% sustained efficiency. Hetzner dedicated servers default to a 1 Gbps uplink, at which the planetary pull takes **7.7 h**. State the uplink assumption.
- `lite` on macOS Docker runs inside a VM; bind-mounted mmap I/O through the VM file-sharing layer is much slower than native. Recommend native or Apptainer-on-Linux for `lite`.

*Fix:* Restate as "Index 1 resident in ≈3–5 s; full warm-up via a `vmtouch`-style prefetch of hot files; cold p99 SLO stated separately".

### DB-18 (S3) Sidecar joins are not "<1 ms"
*Where:* `DEPLOYMENT` §3.2.

DuckDB keeps no persistent indexes over Parquet. A view joining a 2 GB `(token_a, token_b)` sidecar builds a hash table per query unless the query filters first: ~10⁸ rows at ~10⁸ rows/s/core is seconds. Slot-level sidecars keyed by `(dense_contig_id, slot_idx)` are up to 3.6×10⁹ rows, not "10 MB–2 GB".

*Fix:* Require sidecars to be sorted by key and joined *after* the selective predicate (the planner applies the query filter first). The daemon keeps a key-sorted mmap and probes it by binary search for small result sets. Fold heavily used sidecars into the base at compaction.

### DB-19 (S4) File, row-group and block geometry are inconsistent
`PLAN`: Layer-2 blocks of 4,096 contigs. `critique_2`: blocks of 65,536 windows, and "1 Row Group == 1 Index-1 Window/Contig Block" with 65,536-row row groups. These units differ. Specify the target file size (1–4 GB), rows per row group, page size (64–256 KB for random-access columns) and page-index writing. Note that DuckDB is listed as *reading* page indexes only ([Parquet status](https://parquet.apache.org/docs/file-format/implementationstatus/)), so writing them needs another writer such as arrow-rs or parquet-cpp. Cache footers in the daemon.

### DB-20 (S3) Build determinism, provenance and durability are absent
- No table records source accession, source release, tool versions or parameters per contig or slot. Add `provenance` columns or a `build_manifest.json` per epoch, keyed by stable contig key.
- RAID-0 holds `epoch_002_user_private`, the only non-reproducible data, with no backup. Snapshot user epochs and sidecars to the Storage Box.
- The atomic symlink flip must keep old files until the daemon re-mmaps (reference-counted release directories).

### DB-21 (S4) Prior-art claim to correct
`PRIOR-ART.md` §5 still describes "Z-order/Hilbert-sorted Parquet" and an "inverted positional n-gram index", which v2.0 has replaced. Align it with the dual projection and the offset-indexed tape store.

---

## Workload model and SLOs (proposed for `PLAN.md` §1.2)

Concurrency: 8 simultaneous clients (MCP agents plus a UI). Token draws are stratified by df: hubs (df > 10⁶), mid (10³–10⁶), tail (<10³). Report p50/p99 per stratum, warm (after a prefetch) and cold (after `drop_caches`).

| Class | Example | Served by | Output bound | Warm p50 / p99 | Cold p99 |
|---|---|---|---|---|---|
| Q1 count | `COUNT` loci with A, B within ±k ORFs, filtered by tier/phylum | Index 1 | scalar | 5 / 50 ms | 500 ms |
| Q2 neighbourhoods | fetch 1,000 contig tapes | tape store | 1,000 contigs | 5 / 30 ms | 300 ms |
| Q3 anchor profile | partners of A ranked by LLR | pair summary | top 500 | 5 / 50 ms | 300 ms |
| Q4 pair drill-down | loci of (A, B) with geometry | window_postings | ≤10⁴ rows, sampled | 20 / 300 ms | 2 s |
| Q5 NOGD hole | `core(A)&core(B)&~pan(C)&masks` | Index 3 | species set | 5 / 50 ms | 300 ms |
| Q6 co-transmission screen | top-k X for C over 10⁵ candidates | transition_edges | top 100 | 20 / 100 ms | 500 ms |
| Q7 slot prediction | gLM-Domain on one locus, CPU | Index 4 | K = 16 beams | 50 / 200 ms | 1 s |
| Q8 ad-hoc SQL | analytic scans | DuckDB | — | no SLO; memory-capped, cancellable | — |

Full enumeration of hub results is streamed with a cursor and has no latency SLO.

## Benchmark harness at 1% scale (proposed for `PLAN.md` §5 as Benchmark 0)

1. **Sample:** 1% of contigs stratified by tier × biome × phylum, plus the full token dictionary, so the hub/tail token distribution is preserved and not truncated. Repeat at 3% and 10% to fit scaling exponents per build stage. Flag any stage whose exponent is >1.1, such as vehicle clustering or Fitch.
2. **Build** every artifact with the production code. Record per-stage wall time, peak RSS and bytes written, and emit the sizing ledger (DB-01) from measured bytes per record.
3. **Oracle:** an unnested slot table in DuckDB, brute-force evaluated. For 10⁴ random Q1/Q4/Q5 queries, require exact set equality between the index and the oracle (catches DB-02, DB-03, DB-15). For dedup, require shuffled-input determinism and an ANI-error check against `skani` (DB-04).
4. **Performance:** characterise the hardware with `fio` (4 KB random QD1/QD64, 1 MB sequential). Run Q1–Q8 with the stratified token draw under 8-way concurrency, warm and cold. Record `EXPLAIN ANALYZE` for every DuckDB query and `roaring_bitmap_statistics` for every bitmap file. Use `vmtouch` for residency and `perf`/`bpftrace` for page faults.
5. **Pass criteria:** oracle agreement 100%; ledger ×100 within 20% of the documents; SLO table met at 1% and extrapolated to 100% with fitted exponents; rebuild BLAKE3 identical across two runs.
6. **External check:** publish the harness and sample manifest so a second machine or person reproduces the ledger and latency table.

---

## Sources
- DuckDB concurrency: https://duckdb.org/docs/current/connect/concurrency.html
- DuckDB Parquet reader (filter pushdown via zone maps; `file_row_number`): https://duckdb.org/docs/2.0/data/parquet/overview.html
- Parquet implementation status (DuckDB page index: read-only): https://parquet.apache.org/docs/file-format/implementationstatus/
- Parquet page index spec: https://parquet.apache.org/docs/file-format/pageindex/
- CRoaring 64-bit API via Rust bindings: https://docs.rs/croaring and https://docs.rs/croaring-sys
- Cuckoo filter sizing and FPR: Fan, Andersen, Kaminsky, Mitzenmacher, "Cuckoo Filter: Practically Better Than Bloom", CoNEXT 2014.
- Mash distance: Ondov et al., Genome Biology 2016 (D = −(1/k)·ln(2J/(1+J))).
