# bio-box task cards

Read [`AGENTS.md`](../../AGENTS.md) first. One card per session. Evidence for every design choice
here is in [`planning/measurements/2026-11-09_design_stress_test.md`](../measurements/2026-11-09_design_stress_test.md).

Tiers: **L** a small model does it alone. **M** a small model does it and a stronger model or the PI
reviews before merge. **H** design or decision; not for a small model. Model assignment (PI,
2026-11-09): Sonnet for **L** cards, Opus for **M** cards; the orchestrating session audits every result.

Card fields: *Depends* (cards that must be merged first), *Files* (the only files you may edit;
`+` marks a new file), *Do*, *Accept* (commands that must pass), *Stop if*.

Order of phases: 0 → A → B → C and D in parallel → E → F → G → H. Nothing in phase B starts before
card 0-01 is done.

---

## 1. Verdicts on the four hard-to-reverse decisions

| Decision | Verdict | Measured evidence (2026-11-09) | Cost to change now | Cost after `rS1.1` |
|---|---|---|---|---|
| Slot key `(contig_uid, slot_ordinal)` | **Change** to `feature_uid` | One added ncRNA on *E. coli* K-12 rewrote 1,214 of 4,318 slot rows, 1,006 of 3,562 IGR rows and 24,342 of 86,360 pair rows; a gene split rewrote 2,372, 1,946 and 47,430. With a coordinate key: 0–1, 1–2 and 110–130. Worse, 86,340 pair primary keys stayed valid while pointing at different genes | 8 tables, 12 modules (about 95 references); no release exists, so the lock may be regenerated | Every released table is rewritten under new names, `id_history` rows for every slot, and every stored search result that holds an ordinal is invalid |
| Materialized `slot_pairs` | **Change** to `token_postings`, pairs computed at query time | 54 rows and 422 bytes per gene at three token levels on the golden set (40 and 308 at two), against 334 bytes per gene for the sequence shards; the builder holds 544 bytes of heap per row. Projected with GTDB R226 counts: S5 16–23 × 10⁹ rows, 130–460 GB against a 25–60 GB budget line | One builder, four macros, one table; `pairs_of_contig` is kept as the oracle | A released table and its queries are replaced; hundreds of GB already written |
| Release-local token ids | **Keep**, with the API rule in `AGENTS.md` §3.4 | Not stress-tested beyond reading the code: ids are recomputed per release by design and never leave a `release=` partition | none | low: token columns are a release-time join |
| Digest-sorted bgzip shards | **Keep the format**, change what surrounds it | samtools 1.19.2 returns identical bytes for 33,713 of 33,713 gene slices and builds byte-identical `.fai` and `.gzi`. At 1,000× the fixture with 12 shards: 1.15 blocks per ±10-gene neighborhood; a 500-contig export reads 0.36 bytes per bp (contiguous cost 0.31). But: the reader loads the whole index (2 s and 414 MB per 10⁶ contigs), joined lookups scan the whole index (9 of 9 row groups), ingest holds a stage in RAM, and one tape row group holds 8,192 contigs | Reader, ingest, two tables' partitioning, tape writer | The same code, plus a rewrite of released index and tape files |

Also measured: the laptop (arm64, Python 3.13) and the pod image (x86-64, Python 3.14) produce the
same manifest hash for the stage-A build, 15 of 15 files byte-identical. Gene calls equal the
Prodigal 2.6.3 binary for 32,195 of 32,198 genes on ten golden prokaryotes (the 3 differences are
origin-wrapping genes, by design). Against RefSeq CDS for K-12: 95.7% matched by stop codon, 94.4%
of those with the same start.

---

## 2. Target data model

This section is the specification for phase B. Types use the registry's names (`UID` = 16-byte
binary, `COORD` = int64, `STRAND` = int8).

**Identifiers** (`src/biobox/ids.py`)

- `KIND = {"cds": 1, "ncrna": 2, "trna": 3, "crispr_array": 4, "other": 15}` (uint8).
- `feature_uid(contig_uid, start, end, strand, kind)` = BLAKE3-128 of `contig_uid ‖ start (8 bytes, big-endian) ‖ end (8 bytes, big-endian) ‖ strand (1 byte, signed) ‖ kind (1 byte)`.
- A domain segment is not a feature. It is a child of a CDS, identified by `(feature_uid, seg_idx)`, where `seg_idx` is the rank of the domain by `env_start` on the protein, from 0.

**Tables**

| Table | Location, partitions | Columns | Primary key; sort key |
|---|---|---|---|
| `genome_input` (new) | `l0seq/genome_input`, `stage` | `contig_uid` UID, `genome_uid` UID, `is_circular` bool, `kind` string (`prokaryote`, `virus`, `plasmid`), `genetic_code` uint8, `assembly_quality` string | `(contig_uid, stage)`; `contig_uid` |
| `feature` (new) | `l1/feature`, `producer`, `stage` | `feature_uid` UID, `contig_uid` UID, `start` COORD, `end` COORD, `strand`, `kind` uint8 | `(feature_uid, producer)`; `(contig_uid, start, end, strand, kind)` |
| `feature_protein` (replaces `slot_protein`) | `l1/feature_protein`, `gene_call_version`, `stage` | `contig_uid`, `feature_uid`, `start`, `end`, `strand`, `protein_uid`, `genetic_code`, `edit_type` | `(feature_uid, gene_call_version)`; `(contig_uid, start, end)` |
| `protein_occurrence` | unchanged location and partitions | `protein_uid`, `contig_uid`, `feature_uid`, `start`, `end`, `strand`, `dist_to_contig_end_bp` | `(feature_uid, gene_call_version)`; `(protein_uid, contig_uid, start)` |
| `igr` | unchanged location; `igr_version`, `stage` | `contig_uid`, `start`, `end`, `strand`, `left_feature_uid` UID null, `right_feature_uid` UID null | unchanged |
| `contig_dim` (new) | `dim/contig_dim`, `release` | `contig_uid`, `genome_uid`, `contig_len`, `is_circular`, `tier` uint8 null, `assembly_quality` string null, `tax_order` uint32 null, `is_species_rep` bool | `(contig_uid, release)`; `contig_uid` |
| `tape_version_def` (new, `small_ok`) | `l3/tape_version_def.parquet` | `tape_version` string, `producer` string | `(tape_version, producer)`; same |
| `contig_tapes` | `l3/tapes`, `tape_version`, `stage` | `contig_uid`, `contig_len`, `is_circular`, `n_slots`, then lists of length `n_slots`: `feature_uids`, `feature_kinds`, `strands`, `start_bps`, `end_bps`, `delta_bp_prev_slots`, `junction_topology_prevs`, `orfs_to_contig_starts`, `orfs_to_contig_ends`, `slot_flags`, `protein_uids` (null entry for a non-CDS slot) | `(contig_uid, tape_version)`; `contig_uid` |
| `tape_tokens` | `l3/tape_tokens`, `tape_version`, `release` | unchanged | `(contig_uid, tape_version, release)`; `contig_uid` |
| `token_postings` (replaces `slot_pairs`) | `l3/token_postings`, `release`, `tape_version`, `token_bucket` | `token_id` uint64, `contig_uid`, `ordinal` uint32, `seg_idx` uint16, `feature_uid`, `start`, `end`, `strand`, `n_slots` uint32, `is_circular` bool, `contig_len` COORD | `(token_id, contig_uid, ordinal, seg_idx, release, tape_version)`; `(token_id, contig_uid, ordinal, seg_idx)` |
| `seq_locator` | `l0seq/seq_locator`, `stage` (a directory, no longer one file; not `small_ok`) | unchanged | `raw_digest`; `raw_digest` |
| `neighborhoods` | unchanged | add `anchor_feature_uid`, `feature_uid`, `tape_version`; ordinals stay as plain columns | `(contig_uid, anchor_feature_uid, feature_uid, policy_hash)` |
| `search_region` | unchanged | `left_feature_uid`, `right_feature_uid` replace `left_slot`, `right_slot` | unchanged |

Rules that go with the tables:

- `tape_version` = `"t-" +` first 12 hex digits of sha256 over the sorted producer names joined by `\n`. The slots of a tape are the union of the `feature` rows of those producers, ordered by `(start, end, strand, kind)`. `ordinal` is the index in that order.
- `token_postings` holds, per slot: one row per NR90 token and one per NR50 token (`seg_idx` 0), and one row per included domain for family tokens (`seg_idx` as defined above). `token_bucket` is the low byte of the token as hex at the width `biobox.layout.prefix_width` picks, as `anchor_bucket` is today.
- Sequence shards: one ingest batch writes one shard. A batch is a run of source objects, taken in order of object sha256, that holds at most `shard_bases` bases (a single object larger than that is its own batch). Records are sorted by `raw_digest` inside a shard. There is no order across shards.
- Tape and tape-token files: a row group holds whole contigs and at most 131,072 slots.
- `tier`, `tax_order`, `assembly_quality`, `is_species_rep`, `genome_uid` are removed from `protein_occurrence`, `contig_tapes` and pairs. Queries that filter on them join `contig_dim`.

**Query macros** (`src/biobox/catalog.py`)

- `slot_dist(i, j, n, circ)` = `CASE WHEN circ THEN least(abs(i - j), n - abs(i - j)) ELSE abs(i - j) END`.
- `q1_count(a, b, k, tier)`: self-join of `token_postings` on `contig_uid` with `token_id = a` on one side and `b` on the other, `slot_dist BETWEEN 1 AND k`, counting distinct `(contig_uid, ordinal)` of the `a` side; `tier` filters through `contig_dim`.
- `pairs_of(a, b, k)`: the same join returning `contig_uid, anchor_feature_uid, partner_feature_uid, anchor_ordinal, partner_ordinal, orf_offset, junction_topology, delta_bp, flags`. Offsets, topology, gap and flags must equal `biobox.l3.pairs.pairs_of_contig` on the same tape, including its rule for a partner exactly half way round a circle.
- `fusions_of(a, b)`: rows where both tokens sit on the same slot with different `seg_idx`.
- `partners_of(a, k)`: for the loci of `a`, read `tape_tokens` for those contigs and count the tokens within `k` slots. For a token with more than 100,000 loci, refuse with an error that names `pair_stats` (hub routing, PLAN §4.2).
- `bb_seq(contig_uid, start, end, strand)`: one Python function that resolves the contig, slices and reverse-complements. `seq_get`, `seq_slice` and `seq_get_circular` call it and contain no subquery on `contig_index` or `seq_locator`.

---

## Phase 0. Decisions (PI)

### 0-01 Approve the phase-B decisions
**H** · Depends: none
Do: the PI reads §1 and §2 and adds rows to `planning/DECISIONS.md`:
DEC-18 feature key replaces the slot ordinal key; DEC-19 `token_postings` replaces `slot_pairs`;
DEC-20 dimension tables hold taxonomy, tier and quality; DEC-21 one shard per ingest batch, sorted
inside the shard only, and `seq_locator` partitioned by stage (amends DEC-17); DEC-22 tape row groups
are bounded by slots; DEC-23 token ids stay release-local and never cross an API.
Accept: the six rows exist. Phase B cards cite them.

### 0-02 Registry credential and HEAD image on the pod
**H** (human) · Depends: none · Pod: yes
Do: `gh auth refresh -h github.com -s read:packages` on the laptop; on the pod pull
`ghcr.io/ifinkelstein/bio-box-test@sha256:48296e3a6b611b2ef6ba643fbba1d1f759251dd51cbe079481870ca907c6944c`
and `ghcr.io/ifinkelstein/bio-box@sha256:50143ef2df2035e285e9b4def994d5495357d3d51946baa202fff4c5f901702f`
into `~/bio-box/images/`. Create `/stor/work/Wilke/if2453/all-the-data` (setgid) with a README like the scratch one (DEC-11).
Accept: `singularity exec <test sif> pytest -q` is green on the pod; the directory exists.

### 0-03 First-10 questions
**H** · Depends: none
Do: the PI signs off or replaces the ten draft questions in `planning/preflight_gaps.md` §3 (DEC-05 said to ask again before the S1 gate). Card G-09 needs Q1.

---

## Phase A. Repair and record

### A-01 Fix the runner-death race
**M** · Depends: none
Files: `src/biobox/search/runner.py`, `tests/r2/test_micro_s0_search.py`
Background: CI run 41096564184 failed `test_runner_death_is_reported_and_the_worker_does_not_outlive_it`: `status` returned `running` for a killed runner. `runner._alive` uses `waitpid` then `kill(pid, 0)`; a killed runner that is not our child, and not yet reaped by its parent, is a zombie and still answers `kill`. A pid can also be reused.
Do: 1. Reproduce: run that test 200 times in the test image (`pytest --count` is not installed; loop in a small script under `tests/`-free scratch). 2. Replace pid liveness with a lock: the runner opens `<work>/<run_id>/runner.lock` (local disk, never under the data root) and holds `fcntl.flock(LOCK_EX)` for its lifetime; the pending JSON records the lock path; `status` tries `LOCK_EX | LOCK_NB`, and if it gets the lock the runner is dead. 3. Keep the pid only for `cancel`. 4. Remove `_alive` from `runner.py`.
Accept: the test passes 200 of 200 in the image; `uv run pytest -q tests/r2/test_micro_s0_search.py` green.
Stop if: the failure does not reproduce in 200 runs; report the loop output.

### A-02 Correct recorded numbers
**L** · Depends: none
Files: `STAGES.md`, `planning/sources_and_sizing.md`, `PLAN.md`, `planning/measurements/2026-10-30_wilkcomp02.md` (this card is the exception to AGENTS §2 for these four)
Do: replace each value below wherever it appears, labelled *measured, GTDB R226 metadata, 2026-11-09* or *measured, golden set, 2026-11-09*:
- golden set size: 36,475,789 bp in 21 contigs and 15 genomes (not 43 Mbp); gene calling 1.0–1.1 Mbp/s per Haswell thread (36.5 Mbp in 31.9–34.1 s), 2.7 Mbp/s on Apple Silicon.
- S2C (complete genomes, ≤20 per species): 25,458 genomes, 106.4 Gbp, 56,387 contigs, 97.7 M proteins.
- S3 (species representatives under 2 Mbp): 38,750 genomes, 52.4 Gbp, 4.17 M contigs, 53.2 M proteins (the plan said 9–12k genomes and 10–13 Gbp).
- S5 (all species representatives): 143,614 genomes, 465.5 Gbp, 24.95 M contigs, 437 M proteins.
- S6 (all genomes): 732,475 genomes, 2,588 Gbp, 111.5 M contigs, 2,436 M proteins.
- measured storage per unit (golden set): shards 0.309 bytes/bp; `protein_seq` 227 bytes per distinct protein (0.73 bytes/aa); `protein_occurrence` 32 bytes/gene; tapes 30 bytes/slot; IGR 12.7 bytes/row; `cluster_member` 36.5 and `cluster_member_aln` 26.2 bytes per row per level.
Accept: `rg -n "43 Mbp|9–12k|10–13" STAGES.md planning/sources_and_sizing.md PLAN.md` returns no stale value; R0 doc tests green.

### A-03 Stress scripts as a maintained tool
**L** · Depends: none
Files: `scripts/stress/README.md`+, `scripts/stress/*.py`
Do: write the README: what each script measures, its command line, and which card uses it. Do not change script behavior.
Accept: every script runs with `--help` or prints its usage.

### A-04 `doctor` reads the cgroup cap and `/stor` mounts
**L** · Depends: none
Files: `src/biobox/doctor.py`, `tests/r1/test_doctor.py`
Do: report `memory.max` of the current cgroup (v2 path from `/proc/self/cgroup`), and whether each path in a `--paths` argument is visible and writable.
Accept: new unit tests with a fake `/proc` tree pass.

### A-05 Spill-watcher thread race in the K-02 test
**L** · Depends: none
Files: `tests/r2/test_micro_s0_duckdb.py`
Background: `test_k02_k03_oversized_group_by_spills_to_temp_directory` starts a `watch` thread that sums `p.stat().st_size` over the spill directory; DuckDB deletes a temp file between `rglob` and `stat`, the thread dies with `FileNotFoundError`, and pytest reports only a warning (seen 2026-11-09). A dead watcher can under-report the spill and let the test pass or fail for the wrong reason.
Do: make the watcher tolerate a file that vanishes (skip it) and fail the test if the watcher thread raised anything else.
Accept: 50 consecutive runs of that test with `-W error::pytest.PytestUnhandledThreadExceptionWarning` pass.

---

## Phase B. Data model before any release

Every card here cites DEC-18..DEC-23. After each card: `uv run python scripts/update_schema_lock.py`
only if the card says so.

### B-01 One Parquet writer
**M** · Depends: 0-01
Files: `src/biobox/io.py`+, `tests/r1/test_io.py`+
Do: `write_table(batches, spec, path, row_group_rows=122_880, max_slots=None, slots_column=None)` takes an iterator of `pyarrow.RecordBatch`, casts to `spec.schema`, writes zstd Parquet under a temporary name and renames it, with statistics on `spec.sort_key`. It refuses input that is not sorted by `spec.sort_key` (check across batch boundaries). With `max_slots`, it closes a row group before the sum of `slots_column` would pass `max_slots`. Add `duck_batches(con, sql, spec)` that streams a DuckDB result as batches cast to the spec (DuckDB returns `BLOB`; cast to fixed 16- and 32-byte binary).
Accept: tests for the cast, the sort refusal, the slot-bounded row groups, atomic rename (kill between write and rename leaves no final file), and identical bytes for one batch versus many.

### B-02 `feature_uid` and kinds
**L** · Depends: 0-01
Files: `src/biobox/ids.py`, `tests/r1/test_conformance.py`
Do: add `KIND` and `feature_uid` exactly as §2. Hypothesis test: distinct inputs give distinct uids; the same feature on the reverse strand differs.
Accept: tests pass.

### B-03 Registry: new and changed tables
**M** · Depends: B-02
Files: `src/biobox/registry/schemas.py`, `tests/registry/schema_lock.json` (by script), `tests/registry/pending_tables.yaml`, `tests/registry/operations.yaml`, `tests/r0/test_registries.py`
Do: add `genome_input`, `feature`, `feature_protein`, `contig_dim`, `tape_version_def`, `token_postings`; change `protein_occurrence`, `igr`, `contig_tapes`, `tape_tokens`, `seq_locator`, `neighborhoods`, `search_region` as §2; remove `slot_protein` and `slot_pairs`. Regenerate the lock with the script. Add a registry test: no primary key contains a column named `ordinal`, `slot_ordinal`, `anchor_ordinal` or `partner_ordinal`.
Accept: R0 green. Other rings will fail until B-04..B-13 land; list the failing test files in the report.

### B-04 Persist genome inputs
**L** · Depends: B-01, B-03
Files: `src/biobox/l0seq/inputs.py`+, `tests/r2/test_micro_s0_inputs.py`+
Do: `write_genome_input(root, stage, genomes)` and `read_genome_input(root, stage)` for the `genome_input` table. `GenomeInput` objects are built only from this table afterwards.
Accept: round trip on the micro-S0 fixture; conformance check passes.

### B-05 Gene calling writes features
**M** · Depends: B-04
Files: `src/biobox/l1/genecall.py`, `tests/r2/test_micro_s0_genecall.py`
Do: read genomes from `genome_input`; process one genome at a time; write `feature` (producer = gene-call version, kind 1), `feature_protein` and `protein_occurrence` through `io.write_table`; no `slot_ordinal` anywhere. Distinct proteins: write per-genome candidates to a temporary Parquet, then deduplicate by `protein_uid` with DuckDB (`duck_batches`) into the bucketed `protein_seq` files.
Accept: migrated tests pass; gene set equals the old one on the fixture (same starts, ends, strands, proteins); peak RSS at 100× synthetic under 1.5 GB.

### B-06 IGRs reference features
**L** · Depends: B-05
Files: `src/biobox/l1/igr.py`, `tests/r1/test_igr.py`, `tests/r2/test_micro_s0_igr.py`
Do: flanks are `left_feature_uid` and `right_feature_uid`; input is every `feature` row of the tape's producers, not only genes; one contig at a time.
Accept: migrated tests pass; the IGR coordinates on the fixture are unchanged.

### B-07 Tapes by tape version
**M** · Depends: B-05
Files: `src/biobox/l3/tapes.py`, `tests/r1/test_l3.py`, `tests/r2/test_micro_s0_l3.py`
Do: build from `feature` rows of the listed producers; write `tape_version_def`; `feature_uids` list; row groups bounded at 131,072 slots (`io.write_table(max_slots=...)`); drop the dimension columns.
Accept: migrated tests; new test: a tape read for one contig scans one row group whose slot total is ≤ 131,072 (`profile_query`).

### B-08 The stability test
**L** · Depends: B-07
Files: `tests/r2/test_micro_s0_feature_stability.py`+
Do: build the fixture; add one kind-2 feature in an IGR of the M. genitalium contig as a second producer and rebuild features, IGRs and tapes. Assert: every `feature`, `feature_protein` and `protein_occurrence` row of the first build is present unchanged; exactly one new `feature` row; exactly one IGR row replaced by two. Mark `covers("MR-20")` and add MR-20 to `tests/registry/relations.yaml` (allowed for this card).
Accept: the test passes, and fails if `feature_uid` is replaced by an ordinal (show this by a temporary local edit; do not commit it).

### B-09 Tokens
**L** · Depends: B-07
Files: `src/biobox/l3/tokens.py`, tests for tokens in `tests/r1/test_l3.py`
Do: partition `tape_tokens` by `(tape_version, release)`; stream contigs; token ids as today.
Accept: migrated tests pass.

### B-10 Postings builder
**M** · Depends: B-09
Files: `src/biobox/l3/postings.py`+, `tests/r2/test_micro_s0_postings.py`+
Do: one pass over tapes and tape tokens, one row per slot per cluster token, plus one row per included domain from `l1_domains` for family tokens (`seg_idx` per §2). Sort and bucket with DuckDB, write with `io.write_table`. Remove `pairs.run`; move `pairs_of_contig`, `window` and the flags to `tests/oracles/pairs_oracle.py`+ unchanged.
Accept: row count equals slots × cluster levels + included domains on the fixture; LY-02, LY-03, LY-05 hold; peak RSS at 100× synthetic under 1.5 GB.

### B-11 Pair macros over postings
**M** · Depends: B-10
Files: `src/biobox/catalog.py`, `tests/r2/test_micro_s0_l3.py`
Do: implement `slot_dist`, `q1_count`, `pairs_of`, `fusions_of`, `partners_of` as §2. Test `pairs_of` and `q1_count` against the oracle for every token pair on the fixture, on linear and circular contigs, including a circular contig with an even slot count (the half-way rule).
Accept: oracle equality; `profile_query` shows Q1 and Q4 on a non-hub token scan at most 2 row groups per side at 1,000× synthetic.

### B-12 `contig_dim`
**L** · Depends: B-04
Files: `src/biobox/dim.py`+, `src/biobox/catalog.py`, `tests/r2/test_micro_s0_dim.py`+
Do: build `contig_dim` for a release from `genome_input` and `contig_index` (taxonomy columns null until a taxonomy step exists). Catalog view `contig_dim` reads only the catalog's release.
Accept: one row per default contig; a tier filter in `q1_count` works through the join.

### B-13 Sequence index that does not load
**M** · Depends: B-03
Files: `src/biobox/l0seq/store.py`, `src/biobox/l0seq/lookup.py`+, `tests/r2/test_micro_s0_l0seq.py`, `tests/r1/test_lookup.py`+
Do: `KeyLookup(files, key_column)` finds a key by reading Parquet footers once, binary-searching row-group min/max, and reading one row group; LRU of 64 row groups. `SeqStore` uses it for `seq_locator` and `contig_index` and holds no dict of all contigs. `seq_locator` is written per stage and never rewritten.
Accept: on `scripts/stress/builder_limits.py`'s 10⁶-contig index, `SeqStore(root)` opens in under 0.2 s with under 50 MB of heap, and one lookup reads one row group.

### B-14 Sequence macros without index scans
**M** · Depends: B-13
Files: `src/biobox/catalog.py`, `tests/r2/test_micro_s0_catalog.py`, `tests/r2/test_micro_s0_duckdb.py`
Do: register `bb_seq(contig_uid, start, end, strand)`; rewrite `seq_slice`, `seq_get`, `seq_get_circular` on it; bounds and strand errors keep their messages.
Accept: `scripts/stress/macro_scan.py` adapted to `bb_seq`: a 1,000-contig joined fetch shows no scan of `contig_index` or `seq_locator`; existing K and R tests pass.

### B-15 Batch ingest
**M** · Depends: B-13
Files: `src/biobox/l0seq/store.py`, `src/biobox/contracts.py` (`check_shards` only), `tests/r2/test_micro_s0_l0seq.py`, `tests/r2/test_micro_s0_contracts.py`
Do: pass 1 streams every source object once and writes `contig_occurrence`; DuckDB picks, per `raw_digest`, the first occurrence by `(object sha256, ordinal)` that no earlier stage's `seq_locator` holds. Pass 2 writes one shard per batch (§2). `check_shards` checks order inside each shard only.
Accept: MR-01 (input order does not change any output byte); RSS at 1,000× synthetic minus RSS at 100× is under 200 MB with `--shard-bases 16777216`; samtools test (R-07) passes in the image.

### B-16 Search code on features
**M** · Depends: B-07, B-11
Files: `src/biobox/search/nbr.py`, `src/biobox/search/pack.py`, `src/biobox/search/export.py`, `src/biobox/search/orfs.py`, `tests/r1/test_search_units.py`, `tests/r2/test_micro_s0_search.py`, `tests/registry/json_schema_lock.json` (bump the schema version as S-13 requires)
Do: every output row that names a slot carries `feature_uid` and `tape_version`; region flanks are feature uids.
Accept: S-09 and S-13 tests pass after migration.

### B-17 Contracts and foreign keys
**L** · Depends: B-05..B-12
Files: `src/biobox/contracts.py`, `src/biobox/provenance.py`, `tests/r2/test_micro_s0_contracts.py`, `tests/r3/test_golden.py`
Do: replace the foreign-key list with the new tables; `VERSION_PARTITION` for the new tables; golden tests use the new builders.
Accept: R0–R2 green on the laptop; R3 green with `BIOBOX_GOLDEN_DIR`.

### B-18 Phase-B exit measurement
**M** · Depends: B-01..B-17 · Pod: yes (after C-02; until then the PI runs it)
Files: `scripts/stress/*.py`, `planning/measurements/<date>_phase_b_exit.md`+
Do: update the stress scripts to the new model and rerun all of them: cross-host manifest hash, 1,000× synthetic queries, feature stability, independent-tool checks, GTDB projection of every table.
Accept: laptop and pod manifest hashes equal; no query scans more than 2 row groups per lookup; no projected table exceeds its budget in `sources_and_sizing.md` (update the budget lines with measured values, labelled).

---

## Phase C. Running on the pod

### C-01 Host profile
**L** · Depends: none
Files: `profiles/wilkcomp02.yaml`+, `src/biobox/profile.py`+, `tests/r1/test_profile.py`+
Do: the caps of `DEPLOYMENT-AND-MALLEABILITY.md` §0.3 as data: day and night thread caps with their hours, `mem_scope`, `tmp`, `bulk`, `durable`, nice values, the 300 CPU-minute limit, the 64 MB small-file floor.
Accept: loader validates types and refuses unknown keys.

### C-02 `bio-box run-step`
**M** · Depends: C-01
Files: `src/biobox/runstep.py`+, `src/biobox/cli.py`, `tests/r2/test_runstep.py`+
Do: run one step command: check `STOP`; create a staging directory under the profile's `tmp`; start the command in `systemd-run --user --scope -p MemoryMax=<mem_scope>` at the profile's nice (plain subprocess where `systemd-run` is absent); set every thread variable (`OMP_NUM_THREADS`, `POLARS_MAX_THREADS`, DuckDB `threads`, MMseqs2 `--threads`) from the cap for the current hour; on success copy outputs with one `rsync` and rename into place; refuse outputs that are files under 64 MB unless the table is `small_ok`; write a run record `runs/stage=<S>/<step>.json` with wall time, CPU-minutes, peak RSS and the sha256 of each output (outside the release manifest, because times and host names would break cross-host manifest equality); on failure leave no output in the data root. The night thread cap applies only when the load average shows no other compute.
Accept: tests with a fake step: STOP honored; failure leaves nothing; small-file refusal; thread cap by hour, weekday and load; resume skips a step whose run record and outputs exist.

### C-03 Kill and resume
**M** · Depends: C-02
Files: `tests/r2/test_runstep_resume.py`+
Files also: `src/biobox/runstep.py`
Background (audit of C-02, 2026-11-09): two gaps. (a) If `run-step` dies after renaming some outputs into the root and before writing the run record, the next attempt reruns the command and then refuses to publish because the files exist; the step is stuck. (b) A SIGTERM or SIGINT to `run-step` leaves the step command running.
Do: D-03. 1. Before the first rename, write an intent file `runs/stage=<S>/<step>.publishing.json` listing every output with bytes and sha256. On start, if an intent file exists: an output present in the root with the listed sha256 is kept; one still in `.incoming` is renamed in; if every output is then in place, write the run record and exit 0; otherwise delete only the root files that match the intent's sha256, delete the intent and rerun. Never delete a file the intent does not list. 2. Forward SIGTERM and SIGINT to the child's process group and exit after it does, publishing nothing. 3. Kill a step at each phase (running, during rsync, between renames, before the record), also with the staging directory deleted; rerun; outputs equal an uninterrupted run byte for byte.
Accept: the test passes for every kill point; no test sleeps to hit a kill point (inject a hook).

### C-04 `bio-box build-stage`
**M** · Depends: C-02, B-17
Files: `src/biobox/stage.py`+, `src/biobox/cli.py`, `tests/r2/test_build_stage.py`+
Do: a fixed list of steps for a stage, each through `run-step`: inputs → ingest → gene calling → other feature producers → IGRs → tapes → domain scans → clustering → tokens → postings → `contig_dim` → manifest → catalog. No workflow engine; a list and the ledger. A scan step is split into shard processes that each stay under 300 CPU-minutes.
Accept: the fixture builds end to end with one command; a second run does nothing.

### C-05 Search runner under systemd scopes
**M** · Depends: A-01, C-02
Files: `src/biobox/search/runner.py`, `tests/r2/test_micro_s0_search.py`
Do: start workers through the same scope helper as `run-step`; CPU and memory limits come from the scope, not from polling.
Accept: S-05 and S-06 tests pass in the image.

### C-06 Runbook
**L** · Depends: C-04
Files: `RUNBOOK.md`+
Do: launch, pause (`STOP`), resume, cancel a search, clean `/tmp`, who to email, how to read the ledger. Commands only, each one tested.

---

## Phase D. Test strength

### D-01 Mutation testing in the repo
**L** · Depends: none
Files: `pyproject.toml` (add `mutmut` to the dev group), `.github/workflows/mutation.yml`+, `planning/tests/TEST_ARCHITECTURE.md` §11 item 3 (status line only)
Do: configure mutmut for `src/biobox/l0seq`, `l1`, `l3`, with tests `tests/r1` and `tests/r2`; a weekly workflow that runs it in the test image and uploads the survivor list.
Accept: the workflow runs to completion on a manual dispatch.

### D-02 Kill the surviving mutants
**L**, one module per session · Depends: D-01
Files: tests only
Do: for each survivor listed for your module in the measurements file §6 and classified *real gap*, add a test that fails on the mutant. Check a strand sign, a coordinate boundary or a wrap case with an independent oracle (AGENTS §5.2).
Accept: the listed mutants are killed; module score (killed ÷ (killed + survived)) is at least 0.80.

### D-03 Independent tools in ring R3
**L** · Depends: none
Files: `tests/r3/test_golden_independent.py`+, `container/Dockerfile` (add `prodigal=2.6.3` only if the PI adds it to `container/environment.yml`)
Do: port `scripts/stress/step6_export.py` and `step6_compare.py` into tests: every gene slice equals `samtools faidx`; `.fai` and `.gzi` equal samtools' own; gene calls equal Prodigal 2.6.3 except origin-wrapping genes; RefSeq CDS for K-12 (pinned as an L0 truth source with sha256) matched by stop at ≥ 0.95.
Accept: green in the image.

### D-04 Cross-host equality as a gate row
**L** · Depends: C-04
Files: `tests/r3/test_cross_host.py`+, `tests/registry/gates.yaml`
Do: the test reads a manifest hash recorded by another host (path from `BIOBOX_PEER_MANIFEST`) and compares it with its own build; writes a `test_run` row covering D-05.
Accept: passes with the laptop's manifest on the pod.

### D-05 Scale test in CI
**L** · Depends: B-18
Files: `tests/r2/test_scale_layout.py`+
Do: `synth_scale.py` at 100× as a test: assert the row-group and block counts of AGENTS §4.4 and §4.6.
Accept: under 5 minutes in CI.

---

## Phase E. Close S0 and cut `rS0.1`

### E-01 Non-coding RNA producer
**M** · Depends: B-07, C-02
Files: `src/biobox/l1/ncrna.py`+, `src/biobox/registry/schemas.py` (`l1_ncrna`), `tests/r2/test_micro_s0_ncrna.py`+
Do: `cmsearch` with Rfam 15.1 at `Z_nt` (DEC-04) and each model's GA cutoff over contigs; hits become `feature` rows (kind 2, producer `infernal-1.1.5-rfam-15.1`) and `l1_ncrna` rows keyed by `feature_uid`.
Accept: oracle: the `cmsearch` binary's tabular output on the fixture; tapes with both producers build.

### E-02 Rfam spike-in panel (G-04)
**L** · Depends: E-01
Do: `cmemit` 5 sequences from each of 20 fixed families into a synthetic contig at recorded coordinates; recovery reported with an exact binomial interval. Files: `tests/r3/test_g04_spikeins.py`+.

### E-03 Truth tables (G-01, G-02)
**H** then **L** · Depends: D-03
Do (H): the PI chooses licensed truth sources (EcoCyc licensing is unverified; RefSeq and RegulonDB are the fallback) and pins them in `sources.yaml`. Do (L): score gene starts and the lambda gpG-T and T4 intron cases against them.

### E-04 Recipes on paper (S-12), Domainator concordance (S-02d), replays (S-01)
**M**, one per session · Depends: B-16
Do: as specified in `planning/data_layer_and_tests.md` §4.10. S-02d needs Domainator in the test image (PI adds the pin).

### E-05 NFS versus local (K-10)
**M** · Depends: C-04 · Pod: yes, at night
Do: the K-10 test of `data_layer_and_tests.md` §4.7 on the golden release.

### E-06 Build and gate `rS0.1`
**M** · Depends: all of B, C, D-03, D-04, E-01..E-05 · Pod: yes
Do: `bio-box build-stage S0 --release rS0.1` in the HEAD image; run R0–R3 in the image with `--image` and `--release rS0.1`; `python -m biobox.gate S0`.
Accept: the gate query reports every blocking test passing on that release and image, or the report lists exactly which do not. Q5 and Q6 need taxonomy; if the PI has not moved them to the S3 gate, report them as the blockers.

---

## Phase F. Statistics

### F-01 Design `pair_stats`
**H** · Depends: 0-03
Do: a design note fixing the unit of counting (locus, vehicle, species representative), the score test and its sufficient statistics, the permutation and empirical nulls (00-spec invariant 11), and the table schema. It is computed by streaming tapes and `tape_tokens`, never from materialized pairs.

### F-02 `synthworld` generator
**M**, one spec section per session · Depends: B-07
Do: implement `planning/tests/synthworld_spec.md` on ZOMBI2 0.54.0 under the DEC-15 conditions; outputs are FASTA plus truth tables, fed through the normal builders.

### F-03 Implement `pair_stats`, then calibrate
**M** · Depends: F-01, F-02
Accept: MR-18 on null worlds (p-values uniform by KS; FDR at nominal within the binomial bound) and planted pairs recovered at the powered effect size.

---

## Phase G. S1: complete phages

| Card | Tier | Do | Accept |
|---|---|---|---|
| G-01 | H | `sources.yaml` contract for INPHARED (version, URL, license, redistribution flag); the PI confirms the scripts are called, not vendored (AGPL) | file merged |
| G-02 | L | fetch and checksum into `l0/` through `run-step`; lock file | checksums recorded; refetch verifies |
| G-03 | L | `genome_input` for S1 from the INPHARED metadata table: kind `virus`, topology from the source field | row count equals the metadata |
| G-04 | H then L | pin PHROGs in `refdb.yaml`; scan as a library (`O-SCAN`) | hmmsearch-binary oracle on a sample |
| G-05 | M | `bio-box build-stage S1 --release rS1.1` on the pod | ledger reconciles (L-01..L-03) |
| G-06 | M | accretion equivalence A-01..A-07: S0+S1 built incrementally equals built from scratch | content hashes equal |
| G-07 | L | S-03 Z invariance and S-08 projection on S1 | pass |
| G-08 | L | measured rates and bytes per unit for every table, first point of the scaling fit | measurement file |
| G-09 | M | cookbook notebook answering question Q1 from the CLI with sequences | a lab member runs it unaided (S-10) |
| G-10 | H | S1 gate review; write the S2 cards from the template below using G-08's numbers | cards merged |

---

## Phase H. S2 and later

Each stage starts with an **H** planning card: read the previous gate's measurement file, refit
`sizing.yaml`, and write that stage's cards from this template. Small models do not plan a stage.

Template (every stage): source contract → fetch and lock → `genome_input` → stage-specific
producers → `build-stage` → accretion equivalence → Z invariance and projection → measurements →
question cookbook → gate review.

Stage-specific work the planner must card:

| Stage | Adds |
|---|---|
| S2 plasmids (PLSDB) | license check before any redistribution; contig identity and deduplication (`genome_identity`, `contig_member`, skani); plasmid gene calling in meta mode; the DuckLake, K-07 and K-12 decisions |
| S2C complete genomes | selection of 25,458 genomes from GTDB metadata as a reproducible query; host retraining of plasmid gene calls (DEC-13); taxonomy into `contig_dim`; search tier ST2; null-world calibration at stage scale (kill criterion K7) |
| P1 | Hetzner rehearsal: same image digest, manifest hash equal to the pod's; costs into `planning/costs.md` |
| S3 small genomes | 38,750 genomes and 52.4 Gbp, not 9–12k (re-budget first); genetic-code detection; NR90-representative scans with projection (S-08); Q5 and Q6; edge censoring becomes mandatory |
| S4 IMG/VR | local-only source (JGI policy); statistics only in any pack; L2 drift test |
| S5 species representatives | 465.5 Gbp and 24.95 M contigs; go or no-go from the scaling fit (kill criterion K1); memory exception request if needed |
| S6 all genomes | 111.5 M contigs; per-genome work on the pod, global steps on Hetzner; pangenome planes |
