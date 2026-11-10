# bio-box task cards

Read [`AGENTS.md`](../../AGENTS.md) first. One card per session. What is done, in progress and
next is in [`STATUS.md`](STATUS.md), updated at every merge; a gap found while doing or auditing a
card becomes a new card here in the same commit that records the finding. Evidence for every design choice
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
| `tape_version_def` (new, `small_ok`) | `l3/tape_version_def`, `tape_version` (one small file per tape version, so a new version never rewrites a file) | `producer` string | `(tape_version, producer)`; `producer` |
| `contig_tapes` | `l3/tapes`, `tape_version`, `stage` | `contig_uid`, `contig_len`, `is_circular`, `n_slots`, then lists of length `n_slots`: `feature_uids`, `feature_kinds`, `strands`, `start_bps`, `end_bps`, `delta_bp_prev_slots`, `junction_topology_prevs`, `orfs_to_contig_starts`, `orfs_to_contig_ends`, `slot_flags`, `protein_uids` (null entry for a non-CDS slot) | `(contig_uid, tape_version)`; `contig_uid` |
| `tape_tokens` | `l3/tape_tokens`, `tape_version`, `release` | unchanged | `(contig_uid, tape_version, release)`; `contig_uid` |
| `token_postings` (replaces `slot_pairs`) | `l3/token_postings`, `release`, `tape_version`, `token_bucket` | `token_id` uint64, `contig_uid`, `ordinal` uint32, `seg_idx` uint16, `feature_uid`, `start`, `end`, `strand`, `n_slots` uint32, `is_circular` bool, `contig_len` COORD | `(token_id, contig_uid, feature_uid, seg_idx, release, tape_version)`; `(token_id, contig_uid, ordinal, seg_idx)` (`ordinal` orders the rows but is not part of the key) |
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
**Done 2026-11-09:** the PI approved DEC-18 to DEC-23; the rows are in `planning/DECISIONS.md`.

### 0-02 Registry credential and HEAD image on the pod
**H** (human) · Depends: none · Pod: yes
Do: `gh auth refresh -h github.com -s read:packages` on the laptop; on the pod pull
`ghcr.io/ifinkelstein/bio-box-test@sha256:48296e3a6b611b2ef6ba643fbba1d1f759251dd51cbe079481870ca907c6944c`
and `ghcr.io/ifinkelstein/bio-box@sha256:50143ef2df2035e285e9b4def994d5495357d3d51946baa202fff4c5f901702f`
into `~/bio-box/images/`. Create `/stor/work/Wilke/if2453/all-the-data` (setgid) with a README like the scratch one (DEC-11).
Accept: `singularity exec <test sif> pytest -q` is green on the pod; the directory exists.

### 0-04 Slim `token_postings` or keep it self-contained
**H** (PI) · Depends: B-26 · Decide before `rS0.1` is signed (after that the columns are locked, AGENTS §3.8)
Evidence: `planning/measurements/2026-11-10_postings_bytes.md`. A posting is 35.7 bytes (*measured*, `rL.2`): `feature_uid` 17.0, `contig_uid` 4.8, `start` and `end` 8.2, `contig_len` 2.2, `ordinal` 1.4, `n_slots` 1.2, `token_id` 0.9. At GTDB scale that is about 36 GB (*derived* there).
Options: (a) keep every column: a neighbour query needs only postings, no join; (b) drop `start`, `end`, `strand`, `n_slots`, `is_circular`, `contig_len`, which the tape holds: 32.6% smaller (*derived*), and every coordinate or distance question joins the tape; (c) drop only `contig_len`, `n_slots`, `is_circular` (contig attributes, arguably against AGENTS §3.6): 9.4% smaller (*derived*: 3.35 of 35.7), with a join to `contig_dim` or the tape for wrap-around distances only.
Recommendation of the orchestrating session: (a) for now. 36 GB is small beside Pfam scan cost and sequence shards, and the measured queries are fast because they do not join. Revisit only if a budget is missed.

### 0-05 Two S5 storage budgets are exceeded
**H** (PI) · Depends: B-18 · Evidence: `planning/measurements/2026-11-10_phase_b_exit.md` (table "Misses")
At S5 (*derived* from bytes per unit *measured* on the pod): per-genome annotations project to 107 GB against a budget of about 50 GB; tapes plus postings to 75 GB against 25 to 60 GB. The whole S5 store projects to about 497 GB and S6 to about 2.8 TB. The cause is the phase-B decision itself: a 16-byte `feature_uid` beside every gene in `feature`, `feature_protein`, `protein_occurrence`, tapes and postings costs about 144 bytes per gene, about 63 GB at S5.
Options: (a) raise the two budget lines to the measured values: disk is not the binding constraint on the pod, and stable keys were bought on purpose (DEC-18); (b) carry `feature_uid` in fewer tables and join for it; (c) with card 0-04, slim postings.
Recommendation of the orchestrating session: (a). Decide with 0-04, before `rS0.1`.

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

**Expand, migrate, contract.** Every card in this phase leaves the whole suite green. New tables and
columns are added first (B-03 to B-16) while the legacy ones (`slot_protein`, `slot_ordinal`,
`left_slot`, `right_slot`, the dimension columns on fact tables) stay written and readable; B-17
removes them once nothing reads them. A card that adds a column to an existing table edits
`schemas.py` and regenerates the lock in the same commit as the builder that writes it.
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
Do: add the six new tables only: `genome_input`, `feature`, `feature_protein`, `contig_dim`, `tape_version_def`, `token_postings`, exactly as §2 (no location clash with existing tables). Do not change or remove any existing table. Regenerate the lock with the script. Add a registry test: no primary key contains a column named `ordinal`, `slot_ordinal`, `anchor_ordinal` or `partner_ordinal`, except for the tables in a module-level allowlist `LEGACY_ORDINAL_KEYS` that holds exactly the tables that have one today; B-17 empties it.
Accept: the whole suite green.

### B-04 Persist genome inputs
**L** · Depends: B-01, B-03
Files: `src/biobox/l0seq/inputs.py`+, `tests/r2/test_micro_s0_inputs.py`+
Do: `write_genome_input(root, stage, genomes)` and `read_genome_input(root, stage)` for the `genome_input` table. `GenomeInput` objects are built only from this table afterwards.
Accept: round trip on the micro-S0 fixture; conformance check passes.

### B-05 Gene calling writes features
**M** · Depends: B-04
Files: `src/biobox/l1/genecall.py`, `tests/r2/test_micro_s0_genecall.py`
Files also: `src/biobox/registry/schemas.py` and the lock (add non-null `feature_uid` to `protein_occurrence`; its primary key becomes `(feature_uid, gene_call_version)`).
Do: read genomes from `genome_input` when no list is passed; process one genome at a time; write `feature` (producer = gene-call version, kind 1), `feature_protein` and `protein_occurrence` through `io.write_table`. Legacy stays for now: `slot_ordinal` remains a plain column of `protein_occurrence` and `slot_protein` is still written (B-17 removes both). Distinct proteins: write per-genome candidates to a temporary Parquet, then deduplicate by `protein_uid` with DuckDB (`duck_batches`) into the bucketed `protein_seq` files.
Accept: migrated tests pass; gene set equals the old one on the fixture (same starts, ends, strands, proteins); peak RSS at 100× synthetic under 1.5 GB.

### B-06 IGRs reference features
**L** · Depends: B-05
Files: `src/biobox/l1/igr.py`, `tests/r1/test_igr.py`, `tests/r2/test_micro_s0_igr.py`
Files also: `src/biobox/registry/schemas.py` and the lock (add nullable `left_feature_uid`, `right_feature_uid`).
Do: flanks are `left_feature_uid` and `right_feature_uid` (`left_slot` and `right_slot` stay written until B-17); input is every `feature` row of the tape's producers, not only genes; one contig at a time.
Accept: migrated tests pass; the IGR coordinates on the fixture are unchanged.

### B-07 Tapes by tape version
**M** · Depends: B-05
Files: `src/biobox/l3/tapes.py`, `tests/r1/test_l3.py`, `tests/r2/test_micro_s0_l3.py`
Files also: `src/biobox/registry/schemas.py` and the lock; `src/biobox/l3/tokens.py`, `src/biobox/l3/pairs.py`, `src/biobox/catalog.py`, `src/biobox/provenance.py`, `src/biobox/search/*.py` and their tests, only where they name the tape partition.
Do: build from `feature` rows of the listed producers; write `tape_version_def`; the partition key `tape_version` replaces `gene_call_version` for `contig_tapes` and `tape_tokens` (update every reader of those paths in the same commit); add the `feature_uids` list; row groups bounded at 131,072 slots (`io.write_table(max_slots=...)`). The dimension columns stay, nullable, until B-17.
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
Do: one pass over tapes and tape tokens, one row per slot per cluster token, plus one row per included domain from `l1_domains` for family tokens (`seg_idx` per §2). Sort and bucket with DuckDB, write with `io.write_table`. `pairs.run` and `slot_pairs` stay until B-11.
Accept: row count equals slots × cluster levels + included domains on the fixture; LY-02, LY-03, LY-05 hold; peak RSS at 100× synthetic under 1.5 GB.

### B-11 Pair macros over postings
**M** · Depends: B-10
Files: `src/biobox/catalog.py`, `tests/r2/test_micro_s0_l3.py`
Files also: `src/biobox/l3/pairs.py` (removed), `tests/oracles/pairs_oracle.py`+, `src/biobox/registry/schemas.py` and the lock (remove `slot_pairs`), `src/biobox/contracts.py`, `src/biobox/provenance.py`, `scripts/stress/*.py`, and every test that builds or reads `slot_pairs`.
Do: implement `slot_dist`, `q1_count`, `pairs_of`, `fusions_of`, `partners_of` as §2. Then remove `pairs.run` and the `slot_pairs` table; move `pairs_of_contig`, `window` and the flags unchanged to `tests/oracles/pairs_oracle.py`. Test `pairs_of` and `q1_count` against the oracle for every token pair on the fixture, on linear and circular contigs, including a circular contig with an even slot count (the half-way rule).
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

### B-17 Contract: remove the legacy keys
**M** · Depends: B-05..B-16
Files: `src/biobox/registry/schemas.py` and the lock, every builder and reader that still writes or reads a legacy item, `src/biobox/contracts.py`, `src/biobox/provenance.py`, `tests/r0/test_registries.py`, `tests/r2/*`, `tests/r3/test_golden.py`
Do: remove the `slot_protein` table; remove `slot_ordinal` from `protein_occurrence`, `left_slot` and `right_slot` from `igr` and `search_region`, and `tier`, `tax_order`, `tax_order_bucket`, `assembly_quality`, `is_species_rep`, `genome_uid` from `protein_occurrence` and `contig_tapes`; make `neighborhoods` key on feature uids; empty `LEGACY_ORDINAL_KEYS`; replace the foreign-key list and `VERSION_PARTITION` with the final tables. `rg -n "slot_ordinal|left_slot|right_slot|slot_protein" src` must return nothing.
Accept: R0–R2 green on the laptop; R3 green with `BIOBOX_GOLDEN_DIR`.

### B-18 Phase-B exit measurement
**M** · Depends: B-01..B-17 · Pod: yes (after C-02; until then the PI runs it)
Files: `scripts/stress/*.py`, `planning/measurements/<date>_phase_b_exit.md`+
Do: update the stress scripts to the new model and rerun all of them: cross-host manifest hash, 1,000× synthetic queries, feature stability, independent-tool checks, GTDB projection of every table.
Amended 2026-11-09: `scripts/stress/README.md` (card A-03) lists what is stale in each script: `build_s0.py --cluster`, `pilot.py` and `pilot_resume.py` import the removed `biobox.l3.pairs`; `step6_export.py` and `step4_keys.py` key on `slot_ordinal`; `gtdb_pairs.py` is obsolete (delete it and project `token_postings` instead); `resume.sh` applies no memory scope. No compute runs on the laptop any more, so the cross-host check compares the pod with a GitHub Actions runner (a workflow step that builds the fixture store and prints the manifest hash), not with the laptop. Also add the CAMI II pilot driver (measurements §6c; its store is at `/tmp/if2453/pilot_cami` on the pod) to `scripts/stress/`.
Accept: CI-runner and pod manifest hashes equal; no query scans more than 2 row groups per lookup; no projected table exceeds its budget in `sources_and_sizing.md` (update the budget lines with measured values, labelled).

### B-19 Streaming ledger and multi-part tables
**M** · Depends: B-07 (found while auditing B-05 and B-07)
Files: `src/biobox/ledger.py`, the builders that call `ledger.write`, every reader that opens `part-0.parquet` by name, their tests
Do: `ledger.write` accepts an iterator of rows and writes through `io.write_table` (today each builder holds one ledger dict per contig of the stage; S6 has 111.5 M contigs). Tables that are still one `part-0.parquet` per stage (`feature`, `feature_protein`, `protein_occurrence`, `igr`, `contig_index`, `contig_occurrence`) are cut into `part-<k>.parquet` at `layout.MAX_FILE_BYTES` while streaming, as `tapes.run` and `genecall._write_protein_seq` already do; every reader globs `part-*.parquet`.
Accept: with `max_file_bytes` forced small on the fixture, each table has several parts whose concatenation equals the single-file output row for row; suite green.

### B-21 Per-domain family postings
**M** · Depends: B-11 (split out of B-10 on 2026-11-09)
Files: `src/biobox/l3/tokens.py`, `src/biobox/l3/postings.py`, tests
Do: B-10 writes one posting per slot and level from `tape_tokens`, so a multi-domain protein is indexed only by its canonical family. Extend `token_dict` to every model with an included domain on a slot (not only canonical ones), and write one family posting per included domain with `seg_idx` = the domain's rank by `env_start`; `tape_tokens.family_tokens` keeps the canonical family. Then `fusions_of(a, b)` finds two families on one slot.
Accept: on a fixture protein with two Pfam domains, both families have a posting on the same slot with `seg_idx` 0 and 1; token ids of cluster levels are unchanged.

### B-22 `l0_manifest` per stage
**L** · Depends: B-15 (found in its audit)
Files: `src/biobox/l0seq/store.py`, `src/biobox/registry/schemas.py`, tests
Do: `l0_manifest` is one file, `l0/l0_manifest.parquet`, re-read and rewritten by every stage, against AGENTS §3.10. Make it a directory partitioned by stage (`l0/l0_manifest/stage=<S>/part-0.parquet`), written once per stage; readers take all stages. Remove the exemption for that file from the R-09 test.
Amended 2026-11-10 from the B-23 report: `l0seq/store_params.json` is claimed by no table or artifact, so `run-step` refuses to publish the first stage of a new store (`runstep._claim`, LY-10). Register it as an artifact in `schemas.py` in this card, written once by the first stage and verified, never rewritten, by later ones; remove the workaround in `tests/r2/test_ingest_restart.py` that pre-seeds the file. With both fixes a first and a second stage can be ingested under `run-step`; add that as a test. Also remove the `l0_manifest` allowance from `tests/r2/test_micro_s0_release_def.py` (B-25 left it for this card). Raised to **M**.
Accept: after a second stage is ingested, every file of the first stage is byte-identical, with no exception.

### B-23 Ingest: recover the second-read cost and make a crashed stage restartable
**M** · Depends: B-15 (found in its audit)
Files: `src/biobox/l0seq/store.py`, tests
Do: batch ingest reads each source twice and re-hashes kept records in the second pass, which made it 28% slower at 1,000× the fixture (measured). Check length only in the second pass and hash once, unless a measurement shows the hash is not the cost. Separately, an ingest that dies mid-stage leaves files the writer then refuses to overwrite: under `run-step` the stage is published from staging, so make `ingest` safe to rerun into a clean staging directory and document that a direct call needs a clean stage directory.
Amended 2026-11-10 after the agent's measurement (*measured*, wilkcomp02, 3 interleaved runs each, 1,000× synthetic, one thread): batch ingest is 0.9% faster than the pre-B-15 code at the default shard size and 12.4% slower at 16 MiB shards, not 28% slower (that was a laptop figure). A profile shows zlib compression is 182 to 188 s of about 240 s, the second FASTA parse about 25 s, and the second hash about 3 s. The premise of this card was therefore wrong: hashing is not the cost, and the pass-2 hash check stays, because it is what ties the stored bytes to their digest. What remains of the card is the restart safety. Recovering the second parse (pass 1 records byte offsets, pass 2 seeks) is not worth a card until ingest is a measured bottleneck: compression dominates.
Accept: a killed and rerun ingest under `run-step` gives byte-identical output; a source changed between the passes is refused.

### B-24 Streaming rewrites of clustering and domain scans
**M** · Depends: B-19
Files: `src/biobox/l2/cluster.py`, `src/biobox/l1/domains.py`, tests
Do: these two builders still hold every protein, every alignment row and every domain row in Python lists (B-19 only changes their final write). Stream proteins to MMseqs2 from Parquet, parse its output tables with DuckDB, and write through `io.write_parts`; scan domains shard by shard, each shard its own process under the 300 CPU-minute limit (the sharded scan run by hand on 2026-11-09 is the model: 17 shards, 121 CPU-minutes each), merging sorted shard outputs.
Amended 2026-11-10: needed before S2, not after. *Measured* on wilkcomp02: clustering a whole release took 473.5 s and 5.35 GB at 1.13 M proteins and 1,320 s and 9.57 GB at 2.05 M (about n^1.73 and 4.7 KB per protein, *derived* from two points). Also decide, with evidence from the 10.6 M protein run (`rL.3`), whether each release must recluster everything or can use `mmseqs clusterupdate`; if that changes cluster identity rules it is an **H** question: stop and report.
Accept: outputs content-identical on the fixture; peak memory does not grow with protein count (K-copy measurement as in B-05).

### B-25 A release names its stages
**M** · Depends: B-17 (found by the import ladder, 2026-11-10: `planning/measurements/2026-11-10_import_ladder/README.md` finding 1)
Files: `src/biobox/l3/tokens.py`, `src/biobox/l3/postings.py`, `src/biobox/dim.py`, `src/biobox/l2/cluster.py`, `src/biobox/catalog.py` (`build_manifest`), `src/biobox/cli.py`, `src/biobox/registry/schemas.py` (a `release_def` table: `release`, ordered `stages`), tests
Why: every release-level builder reads whatever `stage=*` directories exist. A release built while a later stage was being written took that stage's files in (release `rL.2p` on the pod), against AGENTS §3.9 (explicit lists, never globs).
Do: a release is defined once, by a `release_def` row listing its stages. Every release-level builder and `build_manifest` takes that list and reads only those stages; a stage directory on disk that the release does not list is ignored, and a listed stage that is missing is an error.
Accept: test first: build release R1 over stage S0; add a second stage's files to the store; rebuild R1 into a fresh root from the same `release_def`: every file is byte-identical to the first build. `rg -n "stage=\*" src/biobox` returns nothing in a release-level builder.

### B-28 Query-time readers take the release's stages
**M** · Depends: B-25 (found in its report)
Files: `src/biobox/search/pack.py`, `src/biobox/search/nbr.py`, `src/biobox/search/worker.py`, `src/biobox/l0seq/store.py` (`SeqStore` only), `src/biobox/ledger.py` or `src/biobox/catalog.py` (`ledger_stage`), tests
Why: after B-25 the release builders read only the release's stages, but `pack build`, `nbr`, the search worker and `SeqStore` still read every `stage=*` directory on disk, so a search against release R can return genes of a stage R does not contain. Separately, two releases whose names share a stage prefix (`rS0.1`, `rS0.2`) collide on `ledger/stage=S0/step=cluster`, and the second clustering is refused.
Do: each of these readers takes the file list of the release (from the manifest when one is given, else from `catalog.release_stages`), never a glob (AGENTS §3.9). The clustering ledger is partitioned by release (`ledger/release=<R>/step=cluster/`), a new location in the registry; other ledgers stay per stage.
Accept: test first: a store with stages S0 and S1 and release R1 over S0 only: `pack build`, a search run and `nbr` against R1 return nothing from S1; two releases over the same stage both cluster. `rg -n "stage=\*" src/biobox/search src/biobox/l0seq/store.py` returns only per-stage builders.

### B-29 Postings builder runs out of memory at 21 M slots
Correction, 2026-11-10: the title and the evidence line below say "about 21 M slots"; that was the orchestrating session's error. Release `rL.3` has 12,048,866 slots and 25,822,709 postings (*measured* by the B-29 agent, `planning/measurements/2026-11-10_postings_memory.md`).
**M** · Depends: B-25 · Pod: yes (found by the import ladder, stage C, 2026-11-10 02:53 CST)
Files: `src/biobox/l3/postings.py`, `scripts/stress/synth_scale.py` (a `--slots-per-contig` option), `tests/r2/test_micro_s0_postings.py`, `planning/measurements/<date>_postings_memory.md`+
Evidence (*measured*, wilkcomp02): `postings.run` for release `rL.3` (32,815 contigs, about 21 M slots, of which 5,018 prokaryote contigs hold thousands of slots each) failed after 40 s with `OutOfMemoryException: could not allocate block of size 256.0 KiB (1.8 GiB/1.8 GiB used)` at `con.sql(sql).to_arrow_reader(...)`. It passed at 7.6 M and 8.4 M postings. The query joins whole tape rows to whole token rows, unnests every list by position, and sorts all levels at once inside one 2 GB connection; a spill directory is set, so the operator that holds the memory is one that cannot spill. The B-10 memory test used the 100x synthetic store, whose contigs are fixture-sized, so it never saw long lists.
Do: (1) reproduce first, small: a synthetic store with prokaryote-like contigs (`--slots-per-contig 5000`) on which the current builder fails under the 2 GB limit; that is the failing test (AGENTS §5.1), at a size that runs in under two minutes. (2) Find the operator with `EXPLAIN ANALYZE` or by bisecting the query; do not guess. (3) Fix without raising the memory limit: for example explode tapes and tokens one tape row group at a time (a row group holds at most 131,072 slots, AGENTS §4.6) into scratch Parquet runs, then let DuckDB sort the flat rows, which does spill. Output bytes must not change on the fixture.
Accept: the new test passes; the fixture's postings files are byte-identical to before; on the pod, `postings.run` for `rL.3` in `/tmp/if2453/ladder/root` completes inside the 2 GB DuckDB limit, with wall time and peak RSS recorded (the ladder's driver can resume: `planning/measurements/2026-11-10_import_ladder/README.md`, but call the builder through the image with your tree first on the path, not the ladder's old checkout).

### B-31 Catalog macros that do not prune
**M** · Depends: B-28 (same file) · found by B-18 (`planning/measurements/2026-11-10_phase_b_exit.md`; AGENTS §4.4)
Files: `src/biobox/catalog.py`, `tests/r2/test_micro_s0_duckdb.py`, `tests/r2/test_scale_layout.py`+
Evidence (*measured*, wilkcomp02, 2026-11-10, `scripts/stress/synth_scale.py` at 1,000× with `catalog.profile_query`): a key passed through `token_of(...)` (a scalar subquery) is not pushed into the `token_postings` scan, so `token_rows`, `q1_count`, `pairs_of` and `partners_of` read 20 of 20 row groups where the same key through a plain join reads 1; the near side of `partners_of` filters `token_postings` by `contig_uid`, which cannot prune because the table is sorted by `token_id` (TASKS §2 says to read `tape_tokens` there); `default_digest`, `digest_len` and `contig_len_of` read 9 of 9 row groups for a joined key at 10⁶ contigs (`scripts/stress/macro_scan.py`).
Do: make each of these resolve its key so that the scan prunes: join the key in, or resolve it first and pass a literal, as B-14 did for the sequence macros with a Python UDF. `partners_of` takes its near side from `tape_tokens` by `contig_uid`. Results must not change.
Accept: test first, on the 100× synthetic store built in the test: each macro above, with a literal key and with keys from a joined table, reads at most 2 row groups per key (`catalog.profile_query`); the existing Q-tests give the same rows.

### B-32 Gene calling of short contigs uses its workers
**L** · Depends: none · found by B-18
Files: `src/biobox/l1/genecall.py`, tests, `planning/measurements/<date>_genecall_short_contigs.md`+ · Pod: yes
Evidence (*measured*, wilkcomp02): the CAMI II sample (1.48 M contigs) took 1,989 s with 8 workers against 2,322 s on one thread, at about 2 CPUs busy; complete prokaryotes run at 12 Mbp/s on 16 workers.
Do: profile first on a 100,000-contig slice (where does the parent spend its time: reading sequences, pickling, the scratch writer?). Fix what the profile shows, for example by sending workers batches of genomes sized by bases, not one genome at a time. Output bytes must not change.
Accept: on the slice, 8 workers are at least 4 times faster than 1; the fixture's gene-call files are byte-identical; numbers recorded with labels.

### B-26 Re-measure postings size and correct its budget
**L** · Depends: B-17 · Pod: yes
Files: `src/biobox/l3/postings.py` (the estimate in its docstring or constant), `planning/sources_and_sizing.md` (the postings line only), `planning/measurements/<date>_postings_bytes.md`+
Why: `token_postings` measured 34.6 to 36.8 bytes per row on the pod against the 20 the code states (73 to 84% over; AGENTS §6 threshold is 20%).
Do: on the `rL.2` store on the pod, report bytes per row by column (`pq` metadata: compressed size per column chunk), find which columns carry the excess (`feature_uid` is 16 random bytes and will not compress; check the others' encodings and the sort key). Change no schema. Replace the estimate with the measured value, labelled, and recompute the GTDB-scale projection.
Accept: the measurement file has the per-column table and the command; no unlabelled number remains for postings.

### B-20 IGR version names its producer set
**L** · Depends: B-07
Files: `src/biobox/l1/igr.py`, `src/biobox/provenance.py`, tests
Do: `igr_version` is `<tape_version>-igr1`, so IGRs computed for a second producer set land beside the first instead of at the same path. Today the version names only the gene caller.
Accept: the stability test (B-08) runs both producer sets in one store without removing files.

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
Amended 2026-11-10 (twice). An agent started this card and stopped with evidence (three probe tests, `sh scripts/pod_pytest.sh C-04`): no builder except ingest can write into `$BIOBOX_STAGING`, and `run-step` refuses `store_params.json`, the manifest and the catalog because nothing in the registry claims them. So: **Depends: C-02, B-17, B-25, B-22, C-08.** Rulings for this card:
- Each step is `run-step` running a builder with `root=$BIOBOX_ROOT, out=$BIOBOX_STAGING` (card C-08 gives every builder `out=`). No symlink views.
- `scripts/pod_run.sh` launches `build-stage` once on the host; every step's `run-step` inherits that unit's memory cap. `build-stage` never calls `pod_run.sh`.
- Ingest is one process in this card. Before starting it, estimate its CPU-minutes from the input bytes (*measured* rate: 42.9 CPU-min for 9.6 Gbp on wilkcomp02, `planning/measurements/2026-11-10_import_ladder/`) and refuse a stage that would exceed the profile's `cpu_min_per_process`, naming card B-30 (ingest by batch). The domain scan is one process the same way until B-24 gives a sharded scan; same refusal.
- The one-shard-per-batch decision is DEC-21 (the earlier text said DEC-22).
Accept: the fixture builds end to end with one command; a second run does nothing.

### C-08 Every builder can write to a separate output root
**M** · Depends: B-23, B-29 (they edit two of the same files) · found by the C-04 agent, 2026-11-10
Files: `src/biobox/l0seq/inputs.py`, `src/biobox/l1/genecall.py`, `src/biobox/l1/igr.py`, `src/biobox/l3/tapes.py`, `src/biobox/l1/domains.py`, `src/biobox/l2/cluster.py`, `src/biobox/l3/tokens.py`, `src/biobox/l3/postings.py`, `src/biobox/dim.py`, `src/biobox/ledger.py`, `src/biobox/catalog.py` (`write_manifest`, `generate_catalog`, `define_release`), `src/biobox/registry/schemas.py` (artifacts), `src/biobox/cli.py` (`release` and a new `stage` command group), `tests/r2/test_builders_out.py`+
Why: `run-step` publishes what a command wrote into `$BIOBOX_STAGING` (its module docstring is the contract). Only `ingest` has an `out=` argument (card B-23); every other builder reads and writes the one `root` it is given, so none can run as a step. And the registry claims neither `manifests/release=<R>.json` nor `catalog/release=<R>.duckdb`, so `run-step` refuses them (LY-10).
Do: (1) every builder's `run` takes `out: Path | None = None`: inputs are read from `root` (earlier steps and stages), every output including its ledger and any definition row (`tape_version_def`, `release_def`) is written under `out` (default: `root`, so direct calls and all existing tests are unchanged). A definition row that already exists in `root` is verified, not rewritten. Follow `ingest`'s pattern exactly, including its refusal when the target partition already holds files. (2) Register the manifest and the catalog as artifacts so `run-step` can publish them. (3) CLI: `bio-box stage inputs|ingest|genecall|igr|tapes|domains --root R --stage S [--out DIR]` beside the existing `bio-box release ...`, each defaulting `--out` to `$BIOBOX_STAGING` when that is set.
Accept: tests first: for each builder, run as a `run-step` command on the fixture, the published root is file-for-file identical (sha256) to the root built by direct calls; the whole chain inputs → catalog run step by step under `run-step` gives the same manifest hash as today's fixture build.

### B-30 Ingest by batch processes
**M** · Depends: B-22, C-08
Files: `src/biobox/l0seq/store.py`, `src/biobox/stage.py` (if C-04 has merged), tests
Why: one ingest process runs at about 3.7 Mbp/s (*measured*, wilkcomp02), so a stage over about 60 Gbp exceeds 300 CPU-minutes (AGENTS §7.2). `ingest` decides new sequences and `contig_index` for the whole stage in one call, so it cannot be split today.
Do: split a stage's ingest into a plan step (pass 1 over all sources, itself sharded by source object, writing a scratch decision table: which records are new, their batch, their shard) and one process per batch that writes that batch's shard and rows (DEC-21: one shard per batch), then a merge step for the per-stage index tables. Batches are a deterministic function of the sorted inputs and the batch size. Output must be byte-identical to the single-process ingest.
Accept: with the batch size forced small on the fixture, the batch build equals the single-process build file for file; each process's CPU-minutes are in its run record and under the budget on the 1,000× synthetic stage.

### C-05 Search runner under systemd scopes
**M** · Depends: A-01, C-02
Files: `src/biobox/search/runner.py`, `tests/r2/test_micro_s0_search.py`
Do: start workers through the same scope helper as `run-step`; CPU and memory limits come from the scope, not from polling.
Accept: S-05 and S-06 tests pass in the image.

### C-06 Runbook
**L** · Depends: C-04
Files: `RUNBOOK.md`+
Do: launch, pause (`STOP`), resume, cancel a search, clean `/tmp`, who to email, how to read the ledger. Commands only, each one tested.

### C-07 Scope on the host, around the container
**M** · Depends: C-05 · Pod: yes (the acceptance probe only)
Files: `src/biobox/runstep.py`, `src/biobox/search/runner.py`, `scripts/pod_run.sh`+, `scripts/pod_pytest.sh`, `src/biobox/profile.py`, `profiles/*.yaml`, `tests/r1/test_profile.py`, `tests/r1/test_runstep_units.py`, the two scope-argv assertions in `tests/r2/test_runstep.py` and `tests/r1/test_search_units.py`, `RUNBOOK.md` (if C-06 has merged)
Amended 2026-11-09 after the agent's first stop: inside the image `/sys/fs/cgroup` is empty unless bound (`--bind /sys/fs/cgroup:/sys/fs/cgroup:ro`; *measured*, wilkcomp02, 2026-11-09, probe units `c07-cgprobe*`), so both launch scripts add that bind; "the profile requires one" is a new optional profile key `require_mem_cap` (default false, true for wilkcomp02); the pod has 1117 GB of swap (*measured*, same probe), which is why `MemorySwapMax=0` is needed for a cap to kill.
Why: `systemd-run` is not visible inside the `.sif` (found by the import ladder on wilkcomp02, 2026-11-09), so a `run-step` or a search worker started inside the image cannot open its own scope and today runs with no memory cap.
Do: (1) `scripts/pod_run.sh <unit> <bio-box args>`: the one host-side launcher. It starts `systemd-run --user --unit=<unit> --collect -p MemoryMax=<profile> -p MemorySwapMax=0 nice -n 19 singularity exec <sif by digest> bio-box ...`, detached. (2) Inside the container `run-step` and the search runner detect that they are already in a scope (read `/proc/self/cgroup` and the cgroup's `memory.max`), do not call `systemd-run`, and write the detected cap into the run record. If there is no `systemd-run` and no cap, and the profile requires one, they refuse to start with a message that names `scripts/pod_run.sh`. (3) Add `MemorySwapMax=0` wherever a scope is built.
Accept: unit tests for the three cases (scope available; already capped; neither) using a fake cgroup tree, written first and seen to fail. On the pod: a step launched through `scripts/pod_run.sh` records the cap in its run record, and a deliberate 1 GB over-allocation in a 512 MB unit is killed by the kernel, not by polling.

---

## Phase W. One-command discovery workflow

From the walk-through in [`planning/measurements/2026-11-09_ssrnap_workflow/`](../measurements/2026-11-09_ssrnap_workflow/README.md):
mining single-subunit RNA polymerases with a custom HMM worked in about 30 s of compute but needed
150 lines of ad hoc code, and the recipe runner stopped in six places. These cards close that gap.
They do not depend on phase B, except where noted; cards that touch slot keys wait for B-16.

### W-01 Pack build survives genetic code 15
**L** · Depends: none
Files: `src/biobox/search/orfs.py`, `tests/r1/test_search_units.py`, `tests/r2/test_micro_s0_search.py`
Background: `orfs._translate` calls `pyhmmer.easel.GeneticCode(15)`, which does not exist; 61 phage contigs carry code 15 (TAG read as Gln) and the dark-ORF stage crashes.
Do: translate code 15 as table 11 with TAG → Q (a small explicit table, tested against a hand-written oracle); any other code pyhmmer lacks is a clear error naming the code and the contig, raised before the build starts.
Accept: a fixture contig with code 15 builds a pack; the translation test passes.

### W-02 A user's HMM file as a seed
**M** · Depends: none
Files: `src/biobox/search/recipe.py`, `src/biobox/search/runner.py`, `src/biobox/search/worker.py`, `src/biobox/search/engine.py`, tests
Do: a seed kind `hmm_file` (path; the file's sha256 and the model's name, accession and length are recorded in `search_run`); the recipe hash covers the file's sha256, not its path. A seed kind `sequences` (FASTA of seeds) that builds a profile by aligning the seeds to a single-sequence profile of the first seed and running the builder, as the walk-through did.
Accept: the Pfam `RNA_pol` model given as a file returns the same hits as the same model named from a library; a moved copy of the file gives the same recipe hash.

### W-03 Pack and runs outside the data root
**M** · Depends: none
Files: `src/biobox/search/pack.py`, `src/biobox/search/runner.py`, `src/biobox/cli.py`, `src/biobox/catalog.py` (`attach_search`), tests
Do: the pack is built once per release as a step of `build-stage` (card C-04) into `search_pack/release=<R>/` and is read-only afterwards; a run's output directory is `--out DIR` (default: a per-user directory from the host profile, never the data root), and `attach_search` takes a list of run directories.
Accept: a full run against a read-only store root succeeds and writes nothing under it.

### W-04 Worker threads, exit and kill labelling
**M** · Depends: C-05
Files: `src/biobox/search/worker.py`, `src/biobox/search/runner.py`, tests
Do: the worker sets the Arrow CPU and I/O pools, OpenMP and DuckDB threads to the run's CPU limit before importing anything that starts a pool; any uncaught exception ends the process with a non-zero exit (no hang); a SIGKILL that the runner did not send and that the cgroup does not report as an out-of-memory kill is recorded as `failed`, not `oom`.
Accept: the thread count of a worker with `cpus=2` on a many-core host stays under 16 (measured in a test through `/proc` or `psutil`); an injected exception in round 2 of an iterative run ends the run as `failed` within seconds.

### W-05 Exports a person can open
**M** · Depends: B-16
Files: `src/biobox/search/export.py`, `src/biobox/cli.py`, `tests/registry/json_schema_lock.json` (version bump), tests
Do: `bio-box search hits <run> --format tsv|faa` with one row per gene: accession, header, contig length, start, end, strand, protein length, bitscore, E-value at the reference Z, NR90 and NR50 cluster ids, the gene's Pfam architecture. `bio-box nbr <run> --format tsv` writes one table per locus. Scores rounded to one decimal. `bio-box export proteins --release R` writes the distinct proteins as FASTA named by protein uid.
Amended 2026-11-10: card B-17 changed the nested occurrence objects of a search-hit record (`slot_ordinal` became `ordinal`; `genome_uid` and `tier` were dropped, they live in `contig_dim`) without a version bump, because the lock covers top-level keys only. This card's bump to `search_hit/2` records that change too; add `genome_uid` and `tier` back to the export by joining `contig_dim` at export time.
Accept: the walk-through's `hits.tsv` is reproduced by one command; S-13 passes with the bumped schema version.

### W-06 Neighbor summary
**M** · Depends: B-11
Files: `src/biobox/search/nbr.py`, `src/biobox/cli.py`, tests
Do: `bio-box nbr summarize <run> --window 5g --by pfam,nr50`: per neighbor family or cluster, the number and share of loci, the offsets it occurs at, the number of genomes, and for clusters with no domain the representative's length. Computed from postings and tapes.
Accept: on the fixture, counts equal a brute-force loop over tapes.

### W-07 One command
**M** · Depends: W-02, W-03, W-05, W-06
Files: `src/biobox/cli.py`, `RUNBOOK.md`, tests
Do: `bio-box search hmm --hmm FILE --release R --bits 30 --neighbors 5 --out DIR` writes the recipe it ran, the hit table, the FASTA, the neighbor summary and per-locus tables into `DIR`. A positive control may be given as `--positive ACCESSION:START-END` or a header substring, resolved to a protein uid; without one the command refuses, as the runner does.
Accept: the walk-through's result (1,931 proteins for Pfam `RNA_pol` at 30 bits on `rS1.pilot1`) is reproduced by this one command on the pod (Pod: yes).

### W-08 Small recipe and alarm fixes
**L** · Depends: none
Files: `src/biobox/search/recipe.py`, `src/biobox/search/worker.py`, tests
Do: the shorthand `{protein_uid: X, label: ...}` for a positive control is accepted; the drift alarm writes a warning into `warnings_json` when no domain table exists instead of evaluating against nothing; a catalog view `contig_label` gives one accession and header per contig (`contig_occurrence` has several rows per contig).
Accept: tests for each.

### W-09 Every pyhmmer call takes its thread count from the caller
**L** · Depends: W-04
Files: `src/biobox/search/engine.py`, `tests/r1/test_search_units.py`
Why: `engine.families_of` calls `pyhmmer.hmmsearch` without `cpus`, so pyhmmer sizes its pool from the host (found by the import ladder on the 64-core pod, 2026-11-09). On a shared machine that breaks the thread cap of AGENTS.md §7.2.
Do: `families_of` and every other pyhmmer entry point in `src/biobox` take a `threads` argument with no default and pass it as `cpus`. Callers pass the worker's capped thread count.
Accept: a test that replaces `pyhmmer.hmmsearch` with a recorder and fails when any call lacks `cpus` or passes 0, seen to fail first; `rg "hmmsearch\(|hmmscan\(|phmmer\(" src/biobox` shows `cpus=` on every call.

### W-11 Export follow-ups from the W-05 audit
**L** · Depends: W-05, B-25
Files: `src/biobox/search/export.py`, `src/biobox/cli.py` (`export`, `search hits`), `tests/r2/test_micro_s0_search.py`
Do: (1) `export proteins` without `manifests/release=<R>.json` refuses with exit 2 instead of reading every `protein_seq` file it finds (AGENTS §3.9: readers use the manifest's file list, never globs); migrate the W-05 test that compares "with and without manifest" to assert the refusal. (2) FASTA written for a person (`--format faa`, `export proteins`) has no `*`: a trailing `*` is dropped, and a protein with an internal `*` is written as stored and counted in a note on stderr. First check on the pod whether `protein_seq` holds any `*` at all (`rL.2`); if it holds none, add only a test that says so and change nothing.
Accept: tests for each, written first.

### W-12 Hit and locus exports read only the manifest's files
**L** · Depends: B-28, W-11 (found in the W-11 report; AGENTS §3.9)
Files: `src/biobox/search/export.py`, `src/biobox/cli.py` (`search hits`, `nbr`), `tests/r2/test_micro_s0_search.py`
Do: `hits_tsv`, `hits_faa` and `locus_tables` still fall back to a registry glob when the release has no manifest. Make them refuse as `export proteins` does since W-11 (exit 2, naming the manifest and `bio-box release manifest`), or take the release's stage list through the readers B-28 changed; no path may glob over stages.
Accept: tests first: each command without a manifest exits 2 with that message; with a manifest and an extra stage on disk, output is unchanged.

### B-27 Lookups that cannot prune: contig labels and cluster membership
**M** · Depends: B-25 (found in the W-05 report; AGENTS §4.4)
Files: `src/biobox/registry/schemas.py` and lock (two tables), their builders (`src/biobox/l0seq/store.py`, `src/biobox/l2/cluster.py`), `src/biobox/catalog.py` (`contig_label`, macros), `src/biobox/search/export.py`, `src/biobox/search/nbr.py`, tests
Why: `contig_occurrence` is looked up by `contig_uid` but sorted by (`l0_object_sha256`, `record_ordinal`); `cluster_member` is looked up by `member_uid` but sorted by (`level`, `cluster_uid`, `member_uid`). Neither lookup can use min/max statistics, so exports and `nbr summarize` scan the whole table: 10.5 M `cluster_member` rows at `rL.2p`, 190 M at GTDB scale (*derived*: 2 levels x 95 M proteins).
Do: add two lookup tables written by the same builders, sorted by the lookup key: `contig_label` (per stage: `contig_uid`, first source accession and header by (`l0_object_sha256`, `record_ordinal`), `n_records`) and `member_cluster` (per release: `member_uid`, `level`, `cluster_uid`, `rep_uid`). Point the `contig_label` view, the exports and `nbr.summarize` at them. Do not re-sort the existing tables.
Accept: `catalog.profile_query` on the 100x synthetic store: a label lookup and a membership lookup by 100 keys from a joined table each read at most 2 row groups per key; outputs of the exports unchanged byte for byte.

### W-10 `bio-box search run` takes the host profile
**L** · Depends: C-07 (found in its report)
Files: `src/biobox/cli.py`, `tests/r1/test_cli_search.py`+
Why: `Limits.require_cap`, `mem_gb` and `cpus` cannot be set from the command line, so a search started by hand on the pod never refuses to run uncapped and ignores the profile's thread cap.
Do: `bio-box search run` and `search submit` take `--profile PATH` (default: `$BIOBOX_PROFILE` when set). With a profile, `Limits.require_cap` is the profile's `require_mem_cap`, `cpus` is at most `runstep.choose_threads` for the current hour, and `nice` is the profile's. Explicit `--cpus` and `--mem-gb` flags may lower, never raise, what the profile allows.
Accept: tests, written first: with a profile that requires a cap and neither `systemd-run` nor a cgroup cap (injected), the command exits 2 and names `scripts/pod_run.sh`; `--cpus` above the profile's cap is refused.

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

### D-06 The append-only check covers small tables
**L** · Depends: B-22 (found in its report)
Files: `src/biobox/contracts.py`, `tests/r2/test_micro_s0_contracts.py`, `tests/registry/operations.yaml` (add `store_params` to the outputs of O-INGEST)
Why: `contracts.check_append_only` skips every `small_ok` table, with a comment that `l0_manifest` is rewritten per stage by design. After B-22 nothing is rewritten by design, and the skip would hide a rewrite of `seq_locator`, `ledger`, `release_def` or `tape_version_def`.
Do: remove the skip. If a table is then reported on the two-stage fixture, do not re-add the skip: report which table and where it is rewritten (`id_history` is one file and is a candidate).
Accept: a fault-injection test, written first: rewrite one byte of a small table of an earlier stage and the check reports it; the check is clean on the two-stage fixture.

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
