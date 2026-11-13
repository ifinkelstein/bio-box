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

### 0-08 Non-coding RNA: which hits to keep, and whether to scan everything
**H** (PI) · Depends: E-01 · Evidence: `planning/measurements/2026-11-10_ncrna_rate.md` (*measured*, wilkcomp02, 2026-11-10)
Two questions. (1) Overlapping hits. The producer keeps every `cmsearch` hit at or above its model's GA cutoff, because neither the spec nor a decision row says otherwise. On one small genome, four small-subunit rRNA models and three large-subunit models all hit the same operon, and tRNA and tRNA-Sec models hit overlapping loci. Each distinct locus becomes one feature and one tape slot, so overlapping models put several slots where there is one gene. Rfam's own practice is clan competition: within a clan, keep the best-scoring hit per locus. Options: (a) keep all hits (today); (b) clan competition before features are written, keeping the losing hits in `l1_ncrna` flagged but off the tape. (2) Cost. Full Rfam 15.1 scans at 1,319 bases per CPU-second on the golden set (*derived* from measured CPU time): about 2,100 CPU-hours for 10 Gbp and 10,500 for 50 Gbp (*derived*), five times the Pfam scan, with up to 12.3 GiB per `cmsearch` process. Options: (a) full Rfam everywhere; (b) a curated subset (rRNA, tRNA, tmRNA, RNase P, SRP, riboswitches and other models with bacterial, archaeal or phage hits) with full Rfam on the golden set only; (c) full Rfam on species representatives only.
Recommendation of the orchestrating session: (1b) and (2b): the tape should carry one slot per gene, and a five-fold cost over Pfam is not justified for models that cannot hit prokaryotes. Both need the PI's word because they define what an ncRNA feature is.
**Decided 2026-11-11 (DEC-28)**: the configuration of `planning/design/0-08_ncrna_models.md` §6 as ticked there: 588 models by rule, tRNA and rRNA on Rfam's models without the four eukaryotic rRNA models, intergenic regions with a 200-base margin and `--notrunc`, 23 in-CDS models on whole genomes, clan competition. Cards E-07 (amended), E-08, 0-08a, 0-08b.

### 0-07 What names a protein cluster from one release to the next
**H** (PI) · Depends: DEC-24 · Evidence: `planning/measurements/2026-11-10_cluster_id_churn.md` (*measured*, wilkcomp02, 2026-11-10, releases `rL.1`, `rL.2`, `rL.3`)
`cluster_uid` is a hash of the cluster's whole member set, so a cluster that gains one protein gets a new id. From `rL.2` to `rL.3` (2.05 M to 10.6 M proteins), 12.4% of NR90 clusters and 39.9% of NR50 clusters lost their id, and 44.0% of carried-over proteins have a new NR50 id; yet 94.4% of those proteins still sit with the largest fragment of their old cluster (99.6% at NR90), and three quarters of the NR50 clusters that lost their id only gained members. Representatives are not stable either: 30.9% of carried-over proteins have a new NR50 representative. Among NR50 clusters of at least 10 members, 29% were split. So today nobody can follow "the same family" across releases by `cluster_uid`, although the members mostly stay together.
Options: (a) keep `cluster_uid` as it is (it names an exact member set, which is useful as a checksum) and add a lineage table per release that maps each cluster to its best-overlap predecessor with the Jaccard value: best overlap is at least 0.5 for 91.9% of the changed NR50 clusters (98.7% at NR90), and such a match is always mutual; (b) name a cluster by its representative: simple, but the representative changes for a third of proteins at this growth step; (c) incremental clustering, ruled out by DEC-24.
Not measured: no release was clustered twice, so run-to-run variation of MMseqs2 is not separated from real change; one growth path only.
Recommendation of the orchestrating session: (a). It keeps from-scratch clustering and the exact-set id, and gives users a documented way to follow a cluster, with the overlap stated so weak matches are visible.

### 0-06 How a release is clustered as the corpus grows
**Decided 2026-11-10 by the PI: from scratch, not incremental (DEC-24).** The text below is the record of the question. Still to measure in the planned builds: MMseqs2's own peak memory at 10 M and 45 M proteins, and how many cluster ids survive from one release to the next.
**H** (PI) · Depends: B-24 · Evidence: `planning/measurements/2026-11-10_cluster_memory.md`, `planning/measurements/2026-11-10_import_ladder/rL3.md`
After B-24 the Python side of clustering holds about 1 GB at 2.05 M proteins (was 8.9 GB), so the 410 GB projection for 95 M proteins no longer stands: it was our own process's memory, not MMseqs2's. What remains is MMseqs2 itself: at least 5.3 GB at 2.05 M proteins (*measured*, a sampled lower bound); how it grows is not measured. Time: 473 s, 1,320 s, 9,464 s at 1.13 M, 2.05 M and 10.63 M proteins (*measured*), about 36 h projected at 95 M (*derived*), and every release reclusters everything.
Question: keep "each release clusters all its proteins from scratch" (simple; cluster ids depend only on the release's proteins; cost grows faster than linearly), or move to `mmseqs clusterupdate` (incremental; cluster identity then depends on history, which touches DEC rules on cluster uids)? Not needed for the 50 Gbp build; needed before S6.
Recommendation of the orchestrating session: keep from-scratch through S2 and measure MMseqs2's memory at 10 M and 45 M proteins in the planned builds before deciding.

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

### A-07 Two search tests fail intermittently on GitHub Actions
**M** · Depends: C-07 · blocks the merge to `main`
Files: `src/biobox/search/runner.py`, `src/biobox/search/worker.py`, `tests/r2/test_micro_s0_search.py`, `.github/workflows/ci.yml` (only if the cause is in the runner's environment), `planning/measurements/<date>_ci_flakes.md`+
Evidence (*measured*, GitHub Actions `ci` on `plan/task-cards`, 2026-11-10): 7 of the last 12 runs failed (run ids 42823715863, 42822182025, 42819334492, 42807319670, 42798125650, 42798006548, 42795631143); the pod suite passed on the same commits. Two tests: (1) `test_s05_injected_failures_are_states_never_zero_hits[oom-oom]` ends `timeout` ("wall-clock timeout") instead of `oom`, with `Limits(cpus=1, wall_s=3, mem_gb=4)`: the injected allocation does not reach the memory limit within 3 s on the runner; (2) `test_runner_death_is_reported_and_the_worker_does_not_outlive_it` reads `running` after the runner pid is gone (seen in run 42823715863), the race card A-01 was meant to close.
Do: read all seven failed logs first (`gh run view <id> --log-failed`) and tabulate which test failed in each. Find each root cause; do not add a retry, a sleep or a longer timeout as the fix (AGENTS §2, §5.3). For (1): find which limit path the runner takes on GitHub's Ubuntu image (scope, inherited cap or `RLIMIT_AS`) and why the `oom` fault takes over 3 s there; the fault should hit its limit by construction (for example allocate past the limit in one step, or have the test's wall limit not race the fault at all), not by being fast. For (2): after the runner is killed, something still holds the lock or the status is read before the kernel releases it; show which with evidence (an fd inherited by the worker is the first suspect) and fix the liveness rule, not the test's wait. Reproduce each on the pod if you can (stress the machine with a CPU-bound loop under the same scope) or in a draft workflow run; say which.
Accept: each cause stated with evidence; a test for each that fails deterministically on the old code under an injected condition (not by chance); 20 consecutive green `ci` runs on the branch after the fix (trigger them with `workflow_dispatch` if the workflow has it; the orchestrating session pushes).

### A-06 Stress scripts follow the B-28 reader signatures
**L** · Depends: B-28 · Pod: yes (to run them)
Files: `scripts/stress/{builder_limits.py,synth_scale.py,step6_export.py}`, `scripts/stress/pilot_s1/{pilot.py,store_probe.py}`, `scripts/stress/pilot_cami/{pilot_cami.py,probe_cami.py}`, `scripts/stress/README.md`, `tests/r0/test_stress_scripts.py`+
Do: card B-28 made `SeqStore(root, stages)` require its stages and `pack.build` require a defined release; eight call sites in these scripts still call `SeqStore(root)` and now raise `TypeError`. Pass the script's own stage list (or `catalog.ReleaseFiles.of(root, release)` where the script has a release). Add an R0 test that imports nothing heavy: it parses every script under `scripts/stress` with `ast` and fails on a call to `SeqStore` with fewer than two positional arguments, and on any builder call that B-25 gave a required `stages` argument without it, so the next signature change is caught without running the scripts.
Accept: the R0 test, written first and seen to fail on the eight sites; `synth_scale.py` at 10× and `builder_limits.py` run to completion on the pod.

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

### B-33 A release's manifest covers only that release's files
**M** · Depends: B-28 (found in its report)
Files: `src/biobox/catalog.py` (`build_manifest`, `verify_manifest`, `ReleaseFiles`), `tests/r2/test_release_readers.py`, `tests/r2/test_micro_s0_dim.py` (it relies on the old behaviour: migrate it)
Why: `build_manifest` and `verify_manifest` take in every file under any `release=` partition (`token_dict`, `cluster_member`, postings, cluster ledgers of other releases). So release R1's manifest stops verifying as soon as release R2 writes its files, and R2's manifest lists R1's. On the pod this would make every earlier release fail verification after each new one.
Do: a file with a `release=` partition belongs to that release only. A release's manifest lists the stage files of its `release_def` stages, its own `release=` files, and unpartitioned artifacts; `verify_manifest` ignores other releases' files instead of reporting them as strays.
Accept: test first: build R1 and its manifest, then build R2 over the same stage: R1's manifest still verifies with no problems and is byte-identical when rebuilt; R2's manifest lists no R1 file.

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
Added 2026-11-10 from the `rL.3` release (`planning/measurements/2026-11-10_import_ladder/rL3.md`, *measured*): in every neighbour query the partner side reads the whole token level: 99 of 211 row groups for NR50 neighbours of one cluster (1.3 to 1.4 s at 12 M slots, against 0.14 s at 2.4 M), so that read grows linearly with the release. Postings are sorted by `token_id`, and partners are found by position on the anchors' contigs, so they cannot be pruned from postings: take the partner slots from `tape_tokens` of the anchors' contigs (sorted by `contig_uid`; one row group per contig), as TASKS §2 specifies, and resolve token ids to stable keys through `token_dict` afterwards. This is the main item of the card.
Do: make each of these resolve its key so that the scan prunes: join the key in, or resolve it first and pass a literal, as B-14 did for the sequence macros with a Python UDF. `partners_of` takes its near side from `tape_tokens` by `contig_uid`. Results must not change.
Since card D-05 the failing tests already exist: `tests/r2/test_scale_layout.py` holds eight `xfail(strict=True, reason="card B-31")` tests, one per defect, with the row groups each reads today (32 of 32 for the token macros, 10 of 10 for the digest macros). Remove the markers as you fix each; a marker left on a fixed lookup fails the suite (XPASS).
Accept: test first, on the 100× synthetic store built in the test: each macro above, with a literal key and with keys from a joined table, reads at most 2 row groups per key (`catalog.profile_query`); the existing Q-tests give the same rows.

### B-34 Catalog follow-ups from B-31
**M** · Depends: B-31 (found in its report)
Files: `src/biobox/l3/tokens.py` (how `token_dict` is written), `src/biobox/registry/schemas.py` and lock if a table is added, `src/biobox/catalog.py` (`token_of`, `open_catalog`), tests
Do: (1) `token_of(level, key)` looks `token_dict` up by (level, stable key), but the table is sorted by `token_id`, so the lookup cannot prune; it is 1 row group on the fixtures and 11.3 M rows at `rL.3` (*measured* row count). Add a lookup table `token_key` per release sorted by (level, key) written by the same builder, or sort `token_dict` so that both directions prune: decide by measuring both on the 1,000× store, and do not re-sort a table other readers rely on without checking them. (2) `open_catalog` cannot open a second connection to the same catalog file in one process ("Scalar Function with name bb_seq already exists"): register the UDFs per connection so it can. (3) `partners_of` and the sequence macros need `open_catalog`; a plain DuckDB connection to the catalog file fails with "bb_slot_tokens does not exist". Make that failure say what to do: a stored macro of the same name that raises a clear error naming `biobox.catalog.open_catalog`, shadowed by the real one on an opened connection.
Accept: tests first: `token_of` for 100 keys from a joined table reads at most 2 row groups per key at 100× with at least 8 row groups in the lookup table; two connections to one catalog in one process both answer a sequence query; the plain-connection error names `open_catalog`.

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

### C-10 First real build with `build-stage` on the pod
**M** · Depends: C-04, C-09, C-07 · Pod: yes
Files: `planning/measurements/<date>_build_stage_pod.md`+, `scripts/stress/` (input-list helpers only, if needed)
Why: `build-stage` is tested on the fixture only. Nothing has yet been built on the pod through `scripts/pod_run.sh` → `build-stage` → `run-step` with the real profile (200 G cap, 64 MiB floor, thread caps by hour, `/tmp` staging then publish).
Do: in a fresh data root under `/tmp/if2453/c10/root` (local disk; nothing on `/stor` in this card), build two stages and two releases with the current code: stage `S1a` = the first 2,000 INPHARED phage genomes and stage `P1` = the 250 complete prokaryotes of the import ladder (inputs are already on the pod: `planning/measurements/2026-11-10_import_ladder/README.md`, `selection/`), Pfam as the one domain library; release `rC.1` over `S1a`, then `rC.2` over both. Launch each `build-stage` once with `scripts/pod_run.sh`, detached. Exercise, on real data: a `STOP` file mid-build and resume; a second identical call (must start nothing); a kill of the unit mid-step and resume. Compare stage `P1`'s per-table row counts with the ladder's `rL.2` (the code differs by cards B-17 onward, so compare counts and gene coordinates, not bytes).
Accept: both releases verify; every step has a run record with its memory cap recorded as inherited; the three exercises behave as the fixture tests say; per-step wall, CPU-minutes and peak RSS recorded with labels; every refusal or surprise reported with the command and its output. If `build-stage` refuses the Pfam scan for CPU-minutes (294.5 M residues took 2,150 CPU-min on the ladder, over the 300-minute budget), that is the expected stop until card B-24: record it, then finish that step by the sharded scan the ladder used, and say exactly what was done outside `build-stage`.

### C-11 Domain-scan shards run several at a time, with a margin
**M** · Depends: B-24 · needed before any real Pfam scan through `build-stage`
Files: `src/biobox/stage.py`, `src/biobox/l1/domains.py` (only if a shard needs a change to run beside others), `src/biobox/cli.py` (`build-stage` options), `tests/r2/test_build_stage.py`, `RUNBOOK.md`, `tests/r0/test_runbook.py`
Why: card B-24 made `build-stage` scan in shards, each sized exactly to the CPU-minute budget, one after another. On real data that is slow by construction: the 10 Gbp set has 2,777 M new residues, about 68 shards of 41 M residues, each about 300 CPU-minutes; in sequence at 8 threads that is about 42 hours, where 40 threads' worth of shards at once would take about 8.5 (*derived* from 7.30 CPU-min per million residues, *measured* on wilkcomp02). And a shard sized with no margin exceeds the budget whenever a library is slower than Pfam or the host is loaded.
Do: `build-stage` runs up to `--scan-jobs J` shards at once, each as its own `run-step` process with `--scan-cpus C` threads, with J × C at most the thread count `run-step` allows for the hour (`runstep.choose_threads`); defaults from the profile so that one command uses the allowance and no more. Shards are planned at 80% of the budget (a named constant with that reason). Deterministic output whatever J is; a failed shard stops new launches, lets running ones finish, and returns its code; `STOP` is checked before each launch; a second call starts nothing; resume runs only the shards without a run record. No polling tighter than a second and no sleeps in tests (wait on processes).
Accept: tests first: with J of 1, 2 and 4 the merged output is byte-identical to the unsharded scan; a shard made to fail stops the build with the others' records intact, and the next call runs only the missing shard; the number of shard processes alive at once never exceeds J (count through an injected runner); J × C above the hour's allowance is refused.

### C-12 Long builds respect the hour, survive lost scratch, and publish side by side
**M** · Depends: C-11 (found in its report) · needed before a build that runs past 08:00 or longer than a day
Files: `src/biobox/stage.py`, `src/biobox/runstep.py`, `src/biobox/l1/domains.py`, `scripts/pod_run.sh`, `RUNBOOK.md`, tests
Do: (1) `build-stage` decides its thread allowance once, at the start: a build begun at night with 56 threads keeps 56 after the day window opens, against AGENTS §7.2. Re-evaluate `runstep.choose_threads` before each shard launch and before each step, and never start more work than the current hour allows (running shards finish). `pod_run.sh` sets the unit's `CPUQuota` from one thread count at launch: set it to the night cap and let `build-stage` hold the real limit, or document why not. (2) `runstep._remove_leftovers` removes `<root>/.incoming` when it finds it empty, which races with another step publishing; C-11 works around it from `stage.py` with a marker directory. Fix it in `run-step` itself (create the step's incoming directory with parents; tolerate a non-empty `.incoming`) and remove the marker. (3) If scratch is cleared between the scan shards and their merge (a reboot clears `/tmp`), the merge refuses and the shards are never rescanned, because their run records exist: the build is stuck. Make `build-stage` detect a shard whose run record exists but whose scratch run is gone or fails its recorded checksum, and rescan exactly those shards (the shard record holds what is needed to verify a rescan equals the first scan).
Accept: tests first: a build whose clock moves from night to day between launches (inject the clock) never has more threads in flight than the hour allows after the change; two steps publishing at once with no marker both succeed, in both orders, on the `copytree` path and on the `rsync` path where `rsync` exists; a build with one shard's scratch run deleted before the merge rescans that shard only and ends byte-identical.

### C-13 Defects from the first real build (C-10)
**M** · Depends: C-12 (same files) · Evidence: `planning/measurements/2026-11-10_build_stage_pod.md`
Files: `src/biobox/stage.py`, `src/biobox/runstep.py`, `src/biobox/cli.py` (`release verify`), `scripts/pod_run.sh`, `RUNBOOK.md`, `tests/r0/test_runbook.py`, tests
Do, each with a test written first: (1) a killed step leaves its scratch behind (1.2 GB after one killed clustering step, *measured*): give each step a scratch directory named for the step under the profile's scratch, and have the next attempt of that step remove the previous attempt's directory before it starts. (2) Unit logs are out of time order because `run-step`'s and `build-stage`'s own lines are block-buffered: write them line-buffered to stderr with a timestamp, so the last line of a log shows progress. (3) `pod_run.sh` writes no final line after `systemctl --user stop`: trap termination in the wrapper (also in card W-13; do it once, here, if W-13 has not merged). (4) Add `bio-box release verify --root R --release REL`: runs `catalog.verify_manifest` and the append-only check against the previous release when one is named; exit 0 or 1; put it in the runbook's verify section. (5) The Pfam rate constant is 7.30 CPU-min per million residues; the second real scan measured 7.96 (19 shards, 2,343.74 CPU-min for 294.5 M residues, *measured*). Set the constant to the larger measured value, labelled, and note that the 80% shard margin of C-11 covers a further 25%. (6) The ingest estimate is 29% low on a 63 MB stage (fixed start-up cost) and within 4% on a 1 Gbp stage: add a measured constant term, labelled. (7) A `STOP` that appears during a step discards that step's finished work and reruns it (57.7 s of gene calling lost in the test): this is the documented rule; add the fixture test for STOP mid-step that C-10 found missing, and state the cost in the runbook. Do not change the rule.
Accept: the tests; the runbook doc test passes.

### C-14 Runbook errors found by the 10 Gbp rebuild
**L** · Depends: W-13 (same file) · found 2026-11-10 by the rebuild agent on wilkcomp02
Files: `RUNBOOK.md`, `tests/r0/test_runbook.py`
Do: (1) §2 says a finished build prints `all 13 steps done`; without `--library` it is 12 (*measured*, units `r10-S1-1`, `r10-P1-1`). Say that `--library` is optional and what the line is in each case. (2) §2 recommends `--scan-jobs 10 --scan-cpus 4` by day, which is 40 threads: the account's whole cap, which is shared by everything of ours on the pod. Recommend a count that leaves room and say the cap is a total. (3) §2 says shards running at a `STOP` "finish and publish"; §3 says a step that ends while `STOP` exists publishes nothing. Read `stage.py` and `runstep.py`, state the true rule once, and make both sections agree. (4) §0's nested-quote `python -c` one-liners fail when sent through `ssh '…' bash -c "…"` (`unexpected EOF`); replace them with a script file or a form that survives one level of quoting.
Accept: doc tests first for (1) to (3); the runbook's parser test passes.

### C-15 A run record's peak memory is the launcher's, not the step's
**M** · Depends: C-13 · found 2026-11-10 by the B-30 agent
Files: `src/biobox/runstep.py`, `tests/r2/test_runstep.py`
Evidence (*measured*, wilkcomp02, the B-30 probe on the 1,000× stage): every step's run record, the oracle's included, shows a peak RSS of 344.0 MB. The likely cause is that `ru_maxrss` from the wait call counts the forked image of the parent before `exec`. Not yet confirmed.
Do: reproduce first (a step that allocates a known amount, and one that allocates almost nothing, launched from a parent holding more memory than either; the two records must differ and the small one must be below the parent's). Then record the step's own peak: the scope's `memory.peak` where a cgroup is present, else a measure that excludes the pre-`exec` image. Say in the record which source was used.
Accept: the reproduction fails first and passes after; AGENTS §4.2 memory tests that read run records still pass.

### W-14 A search run's peak memory is the runner's, not the worker's
**L** · Depends: C-15 (found in its report)
Files: `src/biobox/search/runner.py`, `tests/r2/test_micro_s0_search.py`
Why: `search_run.peak_rss_bytes` comes from `wait4` on a forked worker, as run records did before C-15, so it is at least the runner's own peak.
Do: reproduce first, as C-15's test does (a runner holding more memory than its worker). Then measure through the same helper as `run-step` (`runstep.measured` and `runstep.step_peak`); do not write a second helper. Record the source beside the number if the table has room without a schema change; if a schema change is needed, stop and report.
Accept: the reproduction fails first and passes after.

### B-35 `partners_of` reads every `tape_tokens` row group
**M** · Depends: B-31, B-34 · found 2026-11-10 by the 10 Gbp rebuild (AGENTS §4.4)
Files: `src/biobox/catalog.py` (the neighbour macros), `src/biobox/registry/schemas.py` only if the sort key must change (then the lock, by `scripts/update_schema_lock.py`), tests
Evidence (*measured*, wilkcomp02, releases `rN.1a` and `rN.1b`, `scripts/stress/build_stage_pod/r10_queries.py`, entry `partners_prok_top_nr50`): the partner side read 19 of 19 and 26 of 26 `tape_tokens` row groups and took 0.65 to 0.85 s, against 0.13 s for the ladder's hand-written query; results equal (12,778 and 2,065 pairs). Contig uids are random, so anchors on more contigs than there are row groups touch all of them. On the full `rN.1` (10 Gbp, 3,509 anchors) `partners_of` took 9.5 s and made 556 row-group reads of the 93-row-group `tape_tokens`, against 0.47 to 0.80 s for the ladder's query; both give 21,034 pairs (*measured*, same script). Its 32-row-group lookup cache is smaller than the table, so row groups are read again. Also: the lookup's row-group cache hides reads from a later `profile_query` on the same connection.
Do: reproduce on the 100× synthetic store with `catalog.profile_query`. Then make the partner lookup read only the row groups of the anchors' contigs when those are few, and state the bound when they are many (the honest bound may be "all row groups above N contigs"; then say so in the macro's doc and test it). Make `profile_query` count reads that the cache served.
Accept: the reproduction fails first; a 10-anchor query reads at most 10 row groups; results unchanged.

### B-36 Clustering is one process far over the CPU budget
**H** (design) · found 2026-11-10 by the 10 Gbp rebuild (AGENTS §7.2)
Evidence (*measured*, wilkcomp02, unit `r10-P2-1`, run record `rN.1.cluster`): clustering 10,620,684 proteins used 2,195.9 CPU-minutes in one process, 7.3 times the 300 CPU-minute rule; the record says `over_cpu_budget: true` and `build-stage` gates ingest and scans but not clustering. The same record shows a peak of 11.5 GiB against the ladder's 46.1 GB and the unit's 45.4 GiB `MemoryPeak` (which counts page cache); which is right is not settled (see card C-15).
Decided 2026-11-10 (DEC-25): an exemption with its own cap now (this card), and the split (card B-37); the exemption ends when B-37 merges.
**M** to implement · Depends: C-13, B-30 (same file) · Files: `src/biobox/stage.py`, `src/biobox/profile.py` or wherever the host profile is read (one new key), `RUNBOOK.md`, `tests/r1/test_stage_estimates.py`, `tests/r2/test_build_stage.py`
Do: `build-stage` estimates the clustering step's CPU-minutes before launching it, from measured points only (473 s, 1,320 s and 9,464 s of wall time at 1.13, 2.05 and 10.63 M proteins on 16 threads, and 195.08 and 2,195.9 CPU-minutes at 2.05 and 10.62 M, all *measured*, wilkcomp02, 2026-11-10; fit and label the fit *derived*, and say that it is extrapolated above 10.6 M). It refuses when the estimate exceeds the profile's `cluster_cpu_min_max`, naming the key. Set that key in the wilkcomp02 profile to 12,000 (*estimate*: about five times the measured 10 Gbp run; the PI raised it to 20,000 on 2026-11-10 for the 50 Gbp build, which allows about 62.7 M proteins, *derived* from the fit). The run record keeps `over_cpu_budget` and gains the reason `exempt: DEC-25`. Other steps keep the 300-minute rule. State the exemption and the cap in the runbook.
Amended 2026-11-10 (the B-36 agent found the card's 195.08 point mislabelled: it is 982,049 proteins at 8 threads, `rC.2.cluster`). The fit's points are the rebuild's three, all at 16 threads (*measured*, wilkcomp02, 2026-11-10, branch `build/r10`, `planning/measurements/2026-11-10_rebuild_10gbp.md`): 84.75 CPU-minutes at 1,134,550 proteins (`rN.1a`, phages only), 284.10 at 2,053,075 (`rN.1b`), 2,195.91 at 10,620,684 (`rN.1`). Fit a power law through the upper two; it overestimates the phage-only point by about 60% (*derived*), which is the safe side for a refusal: say so in the comment. Ignore the wall-time points and the 195.08 above. `src/biobox/runstep.py` is added to Files for one keyword on `run_step` that writes the reason into the record's meta when given, and nothing else. The estimate needs the release's protein count: refuse before any step when every stage of the release has gene-call records, else before the cluster step's command, as the scan estimate does. The cluster check comes after the ingest check in `_refusal`.
Accept: tests first: a release whose estimate is over the cap is refused before the cluster step does any work; one under it runs and its record names DEC-25; a profile without the key refuses clustering over 300 CPU-minutes; the estimate is within 20% of the two upper measured points and below none of the three.

### B-37 Clustering as several processes, each under the budget
**H** (design first, then **M**) · Depends: B-36 · DEC-25
Why: DEC-24 clusters every release from scratch; at 95 M proteins the single process is about 36 h (*derived*, DEC-24).
Design question for a stronger model or the PI before any code: which MMseqs2 modules replace the `linclust`/`cluster` workflow call (`prefilter`, `align`, `clust`, with `--split` or database splits), how the splits are made a pure function of the sorted input, and whether the result is byte-identical to the workflow's at any split count. Write the answer as a short design note with a measurement on the 2.05 M-protein stage; if the result is not identical, stop: that needs a new decision.
Accept (for the implementing card): the split build equals the single-process build row for row on the fixture and on the 2.05 M stage; every process's run record is under 300 CPU-minutes; then remove the DEC-25 exemption.

### B-38 `partners_of` is slow with many anchors
**M** · Depends: B-35 (which fixed the counting, not the time)
Files: `src/biobox/catalog.py` (the neighbour lookups), tests
Evidence (*measured*, wilkcomp02, 2026-11-10): on `rN.1` (10 Gbp, 93 `tape_tokens` row groups) `partners_of` for 3,509 anchors took 9.5 s and made 556 row-group reads, against 0.47 to 0.80 s for the ladder's hand-written query with the same 21,034 pairs (`r10_queries.py`). On `rN.1b`, 1,146 anchor contigs took 0.804 s for 26 row groups (B-35 probe). B-35 showed the lookup reads the right row groups (at most min(contigs, row groups)); the cost is that its 32-row-group cache is smaller than the table, so row groups are read again, and that whole row groups are decoded in Python.
Do: reproduce on a synthetic store with more `tape_tokens` row groups than the cache holds, with `profile_query` (`row_groups_read` now counts re-reads). Then make a many-anchor query read each row group at most once (visit anchors in row-group order, or size the cache to the query within a stated memory bound), and compare the time with the plain SQL join on the same store; if the join is faster above some anchor count, use it there and say where the switch is. Bounded memory (AGENTS §4.1): state the bound.
Accept: the reproduction fails first; no row group is read twice in one query; on the 100× store the many-anchor query is no slower than the join; results equal the B-35 oracle.

### C-16 A killed step or search has no peak memory
**M** · Depends: C-15, W-14 (found in the W-14 report)
Files: `src/biobox/runstep.py`, `src/biobox/search/runner.py`, tests
Why: the measuring helper writes its file after the command ends. When the helper itself is killed (wall timeout, the OOM killer, SIGKILL of the group) nothing is written, so `peak_rss_bytes` is null (`unmeasured`) for exactly the runs where the peak matters most. Before C-15 and W-14 those runs had a number, a wrong one.
Do: reproduce first (a step and a search killed for memory under a cap, and for wall time). Then record a peak for them: where the run has a cgroup of its own, read its `memory.peak` or `memory.events` from the parent after the kill; else have the helper sample and write as it goes, at a stated interval, so the last sample survives. Say in the record which source was used. No schema change in `search_run`; if one is needed, stop and report.
Accept: the reproduction fails first; an OOM-killed run records a peak within 20% of its cap.

### B-39 `partners_of` reports only canonical families
**M** · Depends: B-21, B-38 (same code)
Files: `src/biobox/catalog.py`, tests
Why: since B-21 the postings hold one family posting per included domain, so anchors can be any domain of a protein, but `partners_of` takes partner tokens from `tape_tokens.family_tokens`, the canonical family only. A neighbour's second domain is invisible to it, while `fusions_of` and the Q1/Q4 macros see it.
Do: test first on the B-21 fixture (a neighbour with two domains). Then take family-level partners from the postings of the neighbour slots, keeping the lookup's pruning (AGENTS §4.4), and say in the macro's doc what a partner is at each level. Cluster levels are unchanged.
Amended 2026-11-11 (the B-39 agent stopped with evidence): the postings cannot serve a lookup by slot. `token_postings` is sorted by (token, contig, ordinal), so one slot's family postings touch every family-level row group (9 of 9 on the 100× store, *measured*, `test_b39_family_postings_of_neighbour_slots_cannot_prune`). Ruling: add one nullable column to `tape_tokens`, appended last, holding for each slot the tokens of every family with an included domain on it, in the order of the postings' `seg_idx` (`family_tokens` keeps the canonical family; an empty list where there is none). `l3/tokens.py` writes it from the same domain rows B-21's postings use, and the two must agree (test it). `partners_of` reads it at the family level through `TapeTokenLookup`, with the pruning it has today. Files gain `src/biobox/registry/schemas.py`, `src/biobox/l3/tokens.py` and the lock, regenerated only by `scripts/update_schema_lock.py` (allowed before `rS0.1`, AGENTS §3.8). Row groups stay bounded by slots (DEC-22). Releases built before this card lack the column: `partners_of` then answers from `family_tokens` and says so once on stderr. Report the size of `tape_tokens` before and after on the fixture and the 100× store (*measured*). Remove the `xfail` from the oracle test in the same commit as the fix.
Accept: the test fails first; cluster-level results are byte-identical to before; the family-level result equals a brute-force oracle over the domain table; a 10-anchor family query reads the same row groups as before the card.

### B-40 Measure the per-domain postings step at 10 Gbp
**M** · Depends: B-21 · Pod: yes, when `rN.2` (Pfam on all three stages) exists
Why: B-21's family run ranks every domain row of the release with a DuckDB window and joins it to the unnested slots; it was tested at about 2,000 domain rows only (AGENTS §4.1, §4.2).
Do: on a copy of the release definition in a scratch root (never the rebuild's root), rebuild tokens and postings of `rN.2` with the current code under a 32G cap and 8 threads; record peak memory (run record, after C-15), wall time, CPU-minutes, family postings and total postings, against B-21's derived 3.5 M family postings without P2 and about 14.6 M more with it. If the peak exceeds 16 GB or a number misses its estimate by more than 20%, stop and report.
Accept: a measurement file with every number labelled.

### W-15 Find a positive control by accession without scanning
**M** · Depends: W-07, B-27 (found in the W-07 report; AGENTS §4.4)
Files: `src/biobox/registry/schemas.py` (one lookup table, then the lock by `scripts/update_schema_lock.py`), the builder that writes `contig_label`, `src/biobox/cli.py` (`search hmm --positive`), tests
Why: `search hmm --positive ACCESSION:START-END` scans every row group of `contig_label`, which is sorted by contig, not by accession. Not measured on a large release. A header substring cannot prune at all, and on real releases nearly always names a contig with many genes, so it is refused with a list.
Do: measure first on `rN.1` (read-only) how long the positive lookup takes and how many row groups it reads. If it is over a second, add a per-stage lookup table sorted by accession (accession, `contig_uid`), written by the step that writes `contig_label`, and resolve `--positive` through it with O(1) row groups; keep the header-substring form as the documented slow path. Also accept a protein accession where the release records one. Stop and report if the measurement shows the scan is cheap: then only document it.
Accept: the measurement, labelled; if the table is added, `profile_query` shows one row group per stage for an accession lookup.

### W-16 `pack build` holds every region of the release in memory
**M** · blocks every search on a release of 10 Gbp or more · found 2026-11-10 by the trial search on `rN.1` (AGENTS §4.1, §4.2)
Files: `src/biobox/search/pack.py`, tests
Evidence (*measured*, wilkcomp02, unit `trial-pack`, `planning/measurements/2026-11-10_trial_search_rN1.md`): on `rN.1` (10.6 M proteins) the `unique` and `igr` parts finished; the dark-ORF part reached the unit's 32 GiB cap, then thrashed at about one CPU for 41 minutes without being killed; stopped by hand after 9,528 s and 155.1 CPU-minutes, with 4.66 GB on disk and no manifest. `_build_igr` and `_build_dark_orf` hold every region of the release in Python lists before writing. No search could run.
Do: reproduce first with a memory test on the 100× synthetic store (peak RSS of the pack build under a bound that the current code breaks; state the bound). Then stream: read regions one row group or one contig at a time, write shards as they fill, and keep the output byte-identical to today's on the fixture (the pack id must not change for the same inputs). Print one progress line per part and per shard to stderr (a build takes hours and says nothing today). A build stopped by a signal removes its `.building-*` directory or reuses it on the next call; say which. The refusal message of a search without a pack, and the runbook, must not suggest `--model` is required when the release has no domain library (add `RUNBOOK.md` and `src/biobox/cli.py` to Files for those lines only).
Accept: the memory test fails first and passes after; fixture pack byte-identical; then on the pod (Pod: yes) the pack of `rN.1` builds under a 32G cap, read-only root, at most 4 threads, with peak memory, wall time, CPU-minutes and size recorded.

### W-17 `nbr summarize` reports every neighbour as unannotated when the release has no `member_cluster`
**L** · found 2026-11-10 by the trial search on `rN.1` · a wrong answer with exit 0
Files: `src/biobox/search/nbr.py` or wherever `nbr summarize` resolves cluster levels, `src/biobox/cli.py` (`release verify`, one check), tests
Evidence (*measured*, wilkcomp02): on `rN.1`, built before card B-27, `nbr summarize <1,957 uids> --window 5g --by nr50` exits 0 with `n_keys {"nr50": 0}` and all 24,473 neighbours `unannotated`; an SQL check finds no neighbour without an NR50 cluster. The missing `member_cluster` becomes an empty view.
Do: test first on a fixture release with its `member_cluster` files removed from the manifest. `nbr summarize` (and any export that reads `member_cluster`) then either answers correctly from `cluster_member`, with a note on stderr that it is the slow path, or refuses with exit 2 naming the missing table and the command that builds it; choose the fallback if it is a few lines, and say which. `release verify` reports, as a warning line and not a failure, each lookup table that current readers use and the release lacks. Do not edit `catalog.py` (card B-38 is in it): if the fix must be there, stop and report.
Accept: the test fails first; with the table present, output is unchanged.

### W-18 Small defects from the trial search
**L** · Depends: W-16, B-38 (same files) · found 2026-11-10
Files: `src/biobox/catalog.py`, `src/biobox/cli.py`, `src/biobox/search/recipe.py`, `scripts/pod_run.sh`, `RUNBOOK.md`, tests
Do, each with a test first: (1) `open_catalog` puts DuckDB's temporary directory under the data root by default: put it under the profile's `tmp` or the system temporary directory, never under the root (W-03's rule). (2) A recipe's `protein_accession` positive control is accepted by the recipe parser and refused at submit: make the two agree. (3) `nbr` takes anchors only as one comma-separated argument (64 KB for 1,957 uids): accept `--anchors-file`. (4) `pod_run.sh --collect` loses the unit's CPU and memory totals when it ends: write them to the log before the exit line. (5) Add `bio-box find --release R (--accession ACC[:START-END] | --protein UID)`, printing contig, coordinates, strand, `feature_uid` and protein uid; `search hmm --positive` uses the same function (see card W-15 for its speed).
Accept: the tests.

### C-17 Give `rN.1` and `rN.2` the lookup tables current readers need
**M** · Depends: the rebuild's phase 2 · Pod: yes
Why: the 10 Gbp rebuild runs on the tree at `740d049`, before cards B-27 (`contig_label`, `member_cluster`) and B-21 (per-domain postings). Its releases lack those tables, so `nbr summarize` and the exports are wrong or slow on them (card W-17).
Do: when `rN.2` is built and verified, define a new release `rN.3` over the same three stages and build only its release steps (cluster, tokens, postings, dim, manifest, catalog) with the current code, in the same root, after checking on a copy of the release definition in a scratch root that the current `build-stage` accepts stages built at `740d049` (if it refuses, stop and report what it refuses). No stage file is rewritten (AGENTS §3.10). Record whether the clusters equal `rN.2`'s row for row, the peak memory of the per-domain postings step (card B-40 is then done by this run), and the time. `contig_label` is per stage and stays the fallback view.
Accept: `rN.3` verifies; `nbr summarize --by nr50,pfam` on it answers with clusters and families; the measurement file.

### 0-09 Annotation is too slow for larger imports: choose how to cut it
**H** (PI decision, after three measurements) · raised by the PI, 2026-11-10
Why: annotation is nearly all the cost of an import. Pfam: 443 CPU-hours per 10 Gbp (*derived* from 7.96 CPU-min per M residues, *measured*); Infernal with full Rfam: about 2,100 CPU-hours per 10 Gbp (*derived* from 1,319 bases per CPU-second, *measured*); ingest and gene calling together about 2 (*measured*, `rN.1`). At 50 Gbp Pfam is about 2,200 CPU-hours, two days of the pod's whole allowance; at 1 Tbp about 44,000 (*derived*, linear).
Options, not exclusive: (a) run scan shards on a larger system (TACC): they are sharded and deterministic already; needs an allocation and the image there; no science trade-off. (b) Scan cluster representatives only and infer the rest: at 10 Gbp of diverse genomes NR90 saves 22% and NR50 72% of proteins (*derived* from `rN.1`: 8,278,537 and 2,966,200 of 10,620,684); more on redundant sets; domain calls on members become inferred. (c) Annotate on demand: scan only the neighbours of a search's hits, cached by protein uid; loses the corpus-wide family index for unscanned tiers. (d) A faster first pass for Pfam (MMseqs2 profile search) for the bulk tier or as a prefilter. (e) Infernal on intergenic regions plus a margin, with a curated model subset (ties to card 0-08), dedicated tools for tRNA and rRNA, and a check of the producer's `cmsearch` filter settings (card E-02 measured that default filters at our Z lose 3 of 100 detectable spike-ins).
Recommendation: a tiered corpus (complete genomes fully scanned; larger sets by representatives or on demand), Infernal on intergenic regions with a curated subset, and a TACC allocation in any case.
Measured 2026-11-11 (all three done; each number is *measured* on wilkcomp02 in the file named):
- 0-09a, `planning/measurements/2026-11-11_mmseqs_profile_vs_hmmsearch.md`: MMseqs2 profile search on stage P1 is 27 to 95 times cheaper than `hmmsearch` and recovers 50 to 63% of (protein, model) pairs at 99% precision, 18 to 35% of the pairs within 5 bits of the cutoff. As a prefilter by protein it saves about 13%. Caveat: profiles were sampled from the HMMs, since the seed alignments are not on the pod; repeating with `Pfam-A.seed` needs it added to `refdb.yaml` (PI). Not a replacement for the family index.
- 0-09b, `planning/measurements/2026-11-11_representatives_only.md`: scanning NR90 representatives plus members whose length differs from their representative's by more than 20% scans 71% of residues on `rN.1b` with 99.7% of pairs right and 0.4% of reported pairs wrong; NR50 the same way scans 45% with 3.1% wrong. The saving is redundancy: 58% on the phage set (many isolates per species), 6% on one-genome-per-species prokaryotes. Envelopes transfer through the stored alignment within 2 residues (NR90) and about 10 (NR50) at the 95th percentile.
- 0-09c, `planning/measurements/2026-11-11_infernal_intergenic.md`: on 8 genomes, whole-genome full Rfam is 1,910 CPU-hours per 10 Gbp (*derived* from the measured rate). Intergenic regions with a 200-base margin and `--notrunc` cost 796 and keep 98.0% of hits (lost: RNAs inside coding sequence); without `--notrunc` the truncated-hit re-search of short records eats most of the saving (1,231). With no margin 9 to 13% of hits are lost (spurious ORFs across tRNAs, gene starts over riboswitch leaders). `--rfam` changes nothing at our Z; `--hmmonly` loses half the hits. Cost is additive per model: the 3,969 models with no hit carry 84% of the CPU; the 258 with a hit would cost 314 whole-genome, about 176 intergenic (*derived*, not run as one scan).
Recommendation after the measurements: (1) Pfam on redundant sets: NR90 representatives plus length-different members, labelled as inferred; full scans for one-genome-per-species sets. (2) Infernal: intergenic regions with a 200-base margin, `--notrunc`, and a curated model set drawn from a larger sample than 8 genomes (card 0-08 decides the list). (3) MMseqs2 profile search only as a labelled lower-recall tier, if at all. (4) TACC for whatever remains (card C-18).
The measurement cards, as they were written (cards 0-09a to 0-09c, **M**, Pod: yes): (a) recall and precision of MMseqs2 profile search against `hmmsearch` at the gathering cutoff on stage P1 (918,765 proteins), with CPU-minutes for both; (b) on a redundant set (pick one from `planning/sources_and_sizing.md`), the fraction of proteins that are not NR90 or NR50 representatives, and on P1 how often a member's `hmmsearch` domain set differs from its representative's; (c) Infernal on intergenic regions plus margins of 0, 50 and 200 bases against whole-genome scanning on the golden set and the E-02 panel: speed, and which hits are lost.

### C-18 Run scan shards on TACC
**H** (design first) · Depends: 0-09 · the PI said on 2026-11-10 that scaling to TACC will eventually be needed; not scheduled
Why: the scans are sharded, deterministic and independent, and the pod's 40 to 56 threads are the limit (card 0-09).
Design questions, to answer in a short note before any code: (1) which system and allocation, and its rules for wall time, node sharing and file systems (the PI supplies these; nothing is fetched or submitted before that). (2) A host profile for it: `run-step` needs a memory cap it can detect; under SLURM the job's cgroup should read as an inherited cap (`detect_cap`), to be checked there, not assumed. `scripts/pod_run.sh` is systemd-specific: its counterpart is a batch script per shard or a job array, with the same last-line contract. (3) The image by digest under Apptainer there; reference libraries (Pfam, Rfam) staged once, checked by sha256 against `refdb.yaml`. (4) What moves: a shard's input is a slice of a stage's `protein_seq` (or contigs), its output a scratch run and a record; define a self-contained shard bundle so only inputs go out and runs come back, each checked by sha256, and the merge stays on the pod. (5) Equality: a shard run there must reproduce the bytes of the same shard on the pod; test it on one real shard before trusting a batch (the cross-host check of card B-18, at shard scale). (6) AGENTS §7 and §8 are written for the pod: propose the amendments, for the PI to approve.
Accept (design): the note, and one real Pfam shard run on both hosts with equal sha256.

### B-41 Implement split clustering (the M half of B-37)
**M** · Depends: B-36, B-30 · Design: `planning/design/B-37_split_clustering.md` (read all of it first; §6 is this card)
Files: `src/biobox/l2/cluster.py`, `src/biobox/stage.py`, `src/biobox/cli.py` (`release cluster plan|step|finish|redo`), `src/biobox/registry/schemas.py` (one `ArtifactSpec`, `cluster_plan`, `small_ok`, not in manifests), `RUNBOOK.md`, tests
Established by the design (*measured*, wilkcomp02, MMseqs2 18.8cc5c, 2026-11-11): the two `cluster` workflow calls expand to 38 module calls; with `prefilter`, `rescorediagonal` and `align` run as N processes over consecutive ranges of sorted query keys, joined by `concatdbs --preserve-keys 1`, and the prefilter parts run with `--add-self-matches 1`, every intermediate database and the final clusters equal the workflow's at 1, 2 and 4 parts on the fixture, on 100,000 S1 proteins and on a constructed set with dense families and ties. MMseqs2's own target split (`--split-mode 0`) changes results and must never be used; a part whose log lacks `Process prefiltering step 1 of 1` is refused.
Do first: find out why `prefilter` rebuilt its target index although `createindex` had built one (each part rebuilding it would cost about 196 CPU-minutes per part at 95 M proteins, *derived*); if a shared index cannot be used, say what part size keeps every process under the budget at 10 M and at 50 M proteins. Then implement as §6 of the note: parts cut from the sorted keys at 80% of the budget, never from thread counts or the hour; one `run-step` per module call or part; each step's record holds a sha256 over (key, entry), not over files; what `cluster` publishes is unchanged; lost scratch reruns from the earliest missing step and must reproduce its digest.
Accept: tests first: the hard-coded module list equals the image's `-v 3` log; builds at 1, 2 and 4 parts equal `cluster.run` file for file on the fixture and the constructed set; without `--add-self-matches 1` the intermediates differ; a deleted scratch database is rebuilt to its digest. Then on the pod (Pod: yes) the 2.05 M-protein release `rN.1b` reclustered in parts equals its published clusters row for row, with every run record under 300 CPU-minutes. Then remove the DEC-25 exemption in a separate commit.

### W-19 A search run loads whole pack shards into memory
**M** · Depends: W-16 (found in its report; AGENTS §4.1) · measure first
Files: `src/biobox/search/runner.py`, `src/biobox/search/worker.py`, `src/biobox/search/fastz.py`, tests, `RUNBOOK.md` (§9: leftovers of a killed build or verify)
Why: the runner builds a set of every target name and the worker and `fastz.read_records` read a whole shard; on `rN.1` the `unique` part's first shard holds about 2 GiB of text (*measured* sizes in `planning/measurements/2026-11-11_pack_build_rN1.md`), so a search will probably hold several GB and grow with the release.
Measured by the second trial (*measured*, wilkcomp02, 2026-11-11, unit `trial2-hmm`, `planning/measurements/2026-11-11_trial_search_rN1_second.md`): one `search hmm` on `rN.1` (10.6 M proteins) gave the right answer (1,957 proteins, equal to `hmmsearch` on the export) in 534 s, but the worker peaked at 26.27 GiB and the runner at 5.65 GiB under a 32 GiB cap. Of the 534 s: about 110 s re-verifying the whole 5.5 GB pack; about 50 s building a set of all 10.6 M protein ids to check the positive control; about 200 s loading the whole target into Python (a dict, a text block, a digitized copy); about 120 s of `hmmsearch` itself; 32 s of exports. The verify also writes a `.verify-*` scratch directory into the shared pack directory, so a user who cannot write there cannot search; and no run id or status exists for the first 160 s. A 50 Gbp release would not fit in 32 GiB.
Do, test first for each: (1) the worker hands pack shards to the search tool by path (or streams records) and keeps only hits: peak memory must not grow with the pack. (2) The positive control is checked by a lookup in the pack's metadata, not a set of every id. (3) A search does not re-verify the whole pack: verify once at build, record the result in the pack's manifest, and check only the manifest and file sizes at search time (`pack verify` stays the full check). (4) Nothing is written into the pack directory by a search. (5) The pending record, with the run id, is written before any long work. Add to the runbook's cleaning section the `.building-*`, `.scratch-*` and `.verify-*` directories a killed pack build or verify leaves behind.
Accept: a memory test on the 100× store whose bound does not grow with the pack; a search against a read-only pack directory succeeds; hits unchanged on the fixture; then on the pod (Pod: yes) the same `search hmm` on `rN.1` gives the same 1,957 proteins, with peak memory, wall time and CPU-minutes recorded against the numbers above.

### W-20 Output wording and small defects from the second trial
**L** · Depends: W-19 (same files) · found 2026-11-11
Files: `src/biobox/search/export.py`, `src/biobox/search/nbr.py`, `src/biobox/cli.py`, `scripts/pod_run.sh`, `RUNBOOK.md`, tests
Do, each with a test first: (1) `pod_run.sh`'s usage line (W-18) is always empty in a real unit, because systemd expands `${...}` in the unit's command: escape it or compute the totals outside the unit's command line; check on one real unit. (2) `nbr summarize --by pfam,...` on a release with no domain library exits 0 with every neighbour `unannotated` at the pfam level: say on stderr that the release has no domain library and drop the level, as `search hmm` does. (3) The slow-path note for a missing `member_cluster` prints nine times in one command: once. (4) The `search hits` table stops at 50 rows without saying so: say how many of how many. (5) The pending record stays after a run closes: remove it. (6) Outputs: state in each header or a `README` written beside the files that coordinates are 0-based half-open and how strand and the neighbour offsets are written (one convention for strand: `+1`/`-1`, AGENTS §3.1; offsets relative to the anchor's strand, said so); `tier` empty in every row and `dist_to_contig_end_bp` 0 for genes that run off an edge: fix or drop; the FASTA's `example=` should name the hit's own contig; document `n_records`.
Accept: the tests; the fixture outputs change only as listed, and the W-07 byte-identity test is updated in the same commit.

### C-20 A contig shared by two genomes in different stages stops the build
**H** (PI decision, then **M**) · blocks the 50 Gbp build (unit `r50-P3-1`) · found 2026-11-11
Evidence (*measured*, wilkcomp02, 2026-11-11, log `/tmp/if2453/r50/tmp/logs/r50-P3-1.log`; probes in `/tmp/if2453/c20/`): stage P3's `inputs` and `ingest` were published; `genecall` failed with `KeyError: "contig_uid ee816c57… is in none of the stages ['P3'] this store reads"`. The store is content-addressed: a sequence is stored once, by the first stage that sees it, and gene calling opens the store with its own stage only. 4 of the 54,265 contigs of the 50 Gbp set belong to two genomes in different stages: two plasmids shared by related isolates (P1/P3, P2/P4) and two phage genomes of S1 that are also contigs of prokaryote assemblies (P5, P6). The same bytes get different gene calls in different genomes, even in the same mode and genetic code, because single mode trains on the whole genome: the *Roseovarius* plasmid has 8 genes in P2's genome and 6 in P4's; phage c-st has 197 in S1 (gv-meta) and 198 in P5 (single), 11 of them different. `feature_uid`, `contig_tapes` and `contig_dim` hold one call set and one genome per contig, and `dim.run` refuses a contig in two genomes (an existing test asserts it), so fixing gene calling alone moves the failure to the release's `dim` step.
Options: (A1) the stage that first stored a contig owns its genes, IGRs and tape; a later stage's gene calling reads earlier stages' sequences for training but writes nothing for that contig and records the skip; `contig_dim` keeps the first genome and `dim` no longer refuses. No schema change; the later genome is one contig short in `contig_dim`; the result depends on stage order; needs a decision row amending DEC-13 and a change to the `dim` refusal test. (A2) as A1 for the calls, with genome membership many-to-many (`contig_dim` keyed by contig and genome): a schema and lock change; right in the long run, since identical plasmids across isolates will be common. (B) refuse such a genome at `stage inputs` and drop the 4 genomes (or contigs) from the selection. (C) key features by genome and contig: breaks content addressing; not recommended.
Recommendation: A1 now, A2 as its own card before the first large strain collection.
Decided 2026-11-11 (DEC-26): option A1. Now **M**.
Files: `src/biobox/l1/genecall.py`, `src/biobox/l1/ncrna.py`, `src/biobox/l1/igr.py` or wherever IGRs are built, `src/biobox/l3/tapes.py`, `src/biobox/dim.py`, `src/biobox/stage.py`, tests
Do, test first, with a fixture whose second stage's genome repeats a first-stage contig (one case in the same mode and code, one where a phage contig of the first stage is a contig of a prokaryote genome in the second): `build-stage` passes the earlier stages of the root to gene calling and the ncRNA producer for reading sequences; they train on the whole genome but write rows only for contigs in their own stage's `contig_index`, and record each skipped contig (contig uid, the stage that owns it) in the step's result and ledger; IGRs and tapes already follow `contig_index`; `dim` keeps the owning stage's genome for a shared contig, lists the later row as skipped, and still refuses a real conflict within one stage; the existing refusal test is migrated to say so. A release that lists a later stage without the owning stage is refused with a message that names the owning stage. Check that no key is duplicated in any table a release over both stages reads, and state in the report which results depend on stage order.
Accept: the fixture build fails first with the KeyError above and passes after; on the pod, stage P3 resumes at `genecall` (its published `inputs` and `ingest` stay) and completes.

### C-21 A contig can belong to several genomes
**H** (design) then **M** · Depends: C-20 · DEC-26 defers this; do it before the first strain collection
Why: under DEC-26 a genome that shares a contig with an earlier stage is one contig short in `contig_dim`, and a shared contig carries the gene calls of whichever genome was built first. With one genome per species that is 4 contigs in 54,265 (*measured*); identical plasmids across isolates of one species will make it common.
Design questions: `contig_dim` keyed by (contig, genome), or a separate membership table, and what every join by `contig_uid` then means (a hit on a shared plasmid belongs to several genomes: counts by genome fan out, counts by contig do not); whether gene calls stay one set per contig (content-addressed, the first genome's) or become per genome context, and what that does to `feature_uid`; how releases stay independent of build order. Write the note with the numbers from a real redundant set before any code.

### C-21a `genome_input` lists a shared contig in every genome
**M** · Depends: C-20 · DEC-27 · Spec: `planning/design/C-21_shared_contigs.md` §3, §6
Files: `src/biobox/l0seq/inputs.py`, `src/biobox/registry/schemas.py` (`genome_input` key and sort key only), `src/biobox/l0seq/store.py` (ownership helper only), `src/biobox/l1/genecall.py`, `src/biobox/l1/ncrna.py` (minimal), `src/biobox/dim.py`, `src/biobox/l3/tapes.py` (one line, added 2026-11-12: its join on raw `genome_input` rows put each feature of a contig in N genomes on the tape N times and duplicated posting keys; found by this card's agent), `scripts/stress/build_stage_pod/make_inputs.py`, tests
Do, test first: `genome_input` may list a contig in more than one genome, within a stage and across stages; key (`contig_uid`, `genome_uid`, stage), sort key (`contig_uid`, `genome_uid`). `inputs` stops refusing a contig in two genomes of one stage; it still refuses a contig twice in one genome, a genome listed twice, and genomes that disagree on a shared contig's kind, genetic code or topology. Within a stage the genome with the least `genome_uid` writes the shared contig's calls; across stages the owning stage does (DEC-26). Every genome is still called whole; the other genomes' rows are dropped before writing and each skip goes to the result and the ledger. `contig_dim` keeps one row per contig with the owning genome; the others are listed as skipped (`n_genomes`, `genome_dim` and the membership view are card C-21b). `make_inputs.py` lists every record in its genome; `--drop-repeats` keeps the old behaviour byte for byte, for the stages already made with it. The same stage built with its genomes in two orders gives identical files.
Accept: the tests; no duplicated key in any table a release reads.

### C-21b Genome membership: the view, `genome_dim`, and `contig_dim.n_genomes`
**M** · Depends: C-21a · DEC-27 · Spec: the design note §3 and §6
Files: `src/biobox/registry/schemas.py` (`genome_dim`, one nullable appended column on `contig_dim`), the locks by `scripts/update_schema_lock.py`, `src/biobox/dim.py`, `src/biobox/catalog.py` (a `genome_member` view over the release's stages' `genome_input`, and the macros that report a genome), tests
Do: as the design note specifies: membership read through a view, so it depends on the release's set of stages and not their order, and no earlier file is rewritten; `contig_dim` keeps one row per contig with a deterministic representative genome, `n_genomes`, and genome-level attributes null where its genomes disagree; genome-level attributes per release in `genome_dim`. A lookup of a contig's genomes prunes (AGENTS §4.4; check with `profile_query`, including on a store split into 20 stages). Releases built before this card answer as today.
Accept: tests first on the C-21a fixture; the pruning check.

### C-21c Counts by genome fan out; counts by locus do not
**M** · Depends: C-21b · DEC-27
Files: `src/biobox/catalog.py`, `src/biobox/search/nbr.py`, `src/biobox/search/export.py`, tests
Do: every count that is by genome takes distinct genomes through `genome_member`; every count by contig, locus or neighbour never joins it except as a semi-join. The design note found `nbr summarize`'s `_aggregate` would inflate neighbour and locus counts once a contig has two genomes: test first with a genome that repeats an existing contig, across `q1_count`'s tier filter, a new `q1_genomes`, `partners_of`, `nbr summarize`, `hits.tsv`, the JSON output and the locus tables. A hit on a shared contig lists every genome that carries it.
Added 2026-11-12 (W-23 merge): `bio-box releases` counts genomes as the distinct `genome_uid` of `contig_dim`, which under C-21a holds only each contig's owner genome, so a genome whose every contig another genome owns is not counted. Count them through `genome_member` (or `genome_dim`, card C-21b); add `src/biobox/cli.py` to Files for that count only.
Ruled 2026-11-12 (byte-identity against the design note's new fields): every output keeps one schema on every release; the new columns and keys (`hits.tsv` `n_genomes`, the JSON occurrences' `n_genomes` and `genome_uids`, the locus facts' `n_genomes=` and `genome_uids=`) are appended always, never only when a contig is shared. On a release without shared contigs every output equals today's except those appended fields and the README lines that describe them (there 1 and the contig's one genome). `nbr summarize`'s `genome_basis` does not switch from `contig_uid` to `genome_uid` on a release that has `genome_input` but no `contig_dim`.
Accept: the tests; outputs on a release without shared contigs are byte-identical to before.

### C-21f Measure context-free gene calling
**M** (measurement) · Pod: yes · Depends: C-21a · DEC-27 defers the decision on gene calls to this
Do: on stage P1's 250 genomes, call genes under the design note's policy p3 (viral contigs as today; other contigs of 100 kb or more trained on the contig alone; shorter ones in meta mode) and compare with the stored calls: per contig class, the fraction of genes with identical coordinates, start-only differences, genes gained and lost, and the CPU cost (meta mode was about 5.8 times single mode per base, *derived* in the note). About 1 CPU-hour (*estimate*). An annotation truth set would need a source added to `sources.yaml` by the PI; without it, report agreement with today's calls only and say so.
Accept: a measurement file; the PI then decides whether calls become a function of the contig alone.

### C-22 Split-clustering estimates run 25% low at 10 M proteins
**L** · Depends: B-41 · found 2026-11-11 building `rN.3`
Evidence (*measured*, wilkcomp02, unit `r50-rN3-1`): the three parts of the NR90 `align` call used 263.5, 260.3 and 260.6 CPU-minutes against an estimate of 208.5 each: under the 300 budget, 25% over the estimate. The per-call shares in `cluster.CALL_CPU_S` were measured on 100,000 proteins.
More evidence from the finished `rN.3` (*measured*, `planning/measurements/2026-11-11_build_r50.md`): `cluster.nr50.12.align` ran as one part at 313.84 CPU-minutes against an estimate of 230.5 (+36%, `over_cpu_budget: true`: the one process over 300); `nr90.07.align` +24%, `nr50.16.align` +10%, `nr90.11.prefilter` −33%, `nr50.19.prefilter` −41%. The whole clustering was 2,174.3 CPU-minutes in 64 steps against 2,195.9 in one process. The 50 Gbp build was launched with `--cluster-parts 24` and a guard at 280 CPU-minutes as a stopgap; its run records (`runs/stage=rM.1/`) are a second, larger sample when they exist.
Do: refit the per-call constants from `rN.3`'s run records (every call and part, in `runs/stage=rN.3/`), labelled *measured* at 10.6 M proteins; keep the 100,000-protein figures in the comment; cut parts so the estimate with the new constants stays at 80% of the budget; restate the table of parts at 10 M, 50 M and 95 M.
Accept: a test that the estimate of each `rN.3` part is within 20% of its record.

### W-22 One place for a release's search pack, and no writes into it
**M** · found 2026-11-11 writing `GUIDE.md` and `AGENT_GUIDE.md`
Files: `src/biobox/search/pack.py`, `src/biobox/cli.py`, `profiles/wilkcomp02.yaml`, `RUNBOOK.md`, `GUIDE.md`, `AGENT_GUIDE.md`, tests
Why: the only pack of `rN.1` is in an agent's scratch directory (`/tmp/if2453/w16/packs`, and a copy with its verification record in `/tmp/if2453/w19/packs`); every search must name it with `--pack-dir`. `pack verify`, and `pack build` on an existing pack, still write a `.verify-*` scratch directory into the pack directory (card W-19 moved only the search's).
Do: the profile names one pack directory per host, and `search` and `nbr` find a release's pack there without `--pack-dir`; `pack verify` and a no-op `pack build` write their scratch under the profile's `tmp`; a pack directory may be read-only for everyone but its builder. Update both guides and the runbook.
Accept: tests first; a search with only `--profile` finds the pack; `pack verify` against a read-only pack directory succeeds.

### W-23 Commands a user needs and does not have
**M** · found 2026-11-11 writing the guides · each item small; test first for each
Files: `src/biobox/cli.py`, `src/biobox/search/export.py`, `src/biobox/catalog.py` (reading only, unless an item needs a helper), `GUIDE.md`, `AGENT_GUIDE.md`, `RUNBOOK.md`, tests
Do: (1) `bio-box releases --root R`: one line per release with its stages, proteins, domain libraries, whether it has the lookup tables, and its pack if the profile's pack directory holds one. (2) `bio-box seq --root R --release REL --protein UID` (and `--contig ACC:START-END`): one sequence as FASTA, by key, without exporting the release (today a cluster representative's sequence means exporting 10.6 M proteins and searching the file). (3) `bio-box sql --root R --release REL FILE`: runs a query through `open_catalog`, read-only, temporary files under the profile's `tmp`, and prints the row groups read (`profile_query`). (4) `hits.tsv` gains `stage` and the genome's `kind` beside each gene, so a split by source needs no SQL. (5) `find` and `--positive` accept 1-based inclusive coordinates as `ACC:START..END` (the GenBank form) beside the 0-based `ACC:START-END`, and print both. (6) After a search, one line per output file with its row count.
Accept: the tests; both guides use the new commands in place of their work-arounds.

### W-24 Silent fallbacks in the search readers
**M** · found 2026-11-11 writing `AGENT_GUIDE.md` · each is a wrong or incomplete answer with exit 0
Files: `src/biobox/cli.py`, `src/biobox/search/nbr.py`, `src/biobox/search/export.py`, `src/biobox/catalog.py`, tests
Do, test first for each: (1) `search hits` and `nbr` on a run fall back to the gene-call tape version when the run's pack cannot be found (`_run_tape_version`): refuse, naming the pack. (2) `nbr summarize RUN --release OTHER` mixes one release's run with another release's files: refuse unless the releases are equal. (3) On a release built before card B-21, `fusions_of` and the family postings hold only canonical families and say nothing: detect it (the release has domains but no non-canonical family posting where the domain table has one, or a marker the token step records) and print one note, as `partners_of` does for card B-39's column. (4) `neighbors.tsv` gives `key` and `cluster_uid` the same value and an empty `model_name` at cluster levels: drop the duplicate or fill `model_name` with the representative's uid and length, and say which in the README.
Accept: the tests.

### C-23 A release does not record the code that built it
**L** · found 2026-11-11 writing `AGENT_GUIDE.md`
Files: `src/biobox/stage.py`, `src/biobox/catalog.py` (`catalog_info`), tests
Why: only `STATUS.md` says `rN.1` was built by the tree at `740d049`; readers behave differently on releases built before cards B-21, B-27 and B-39, and nothing in the release says which it is. Run records hold the image digest, but the source was bound in from outside the image.
Do: each run record gains the source tree's identity (the git commit when the tree is a checkout, else a hash of the `biobox` package's files, computed once per process); `catalog_info` gains `built_by` (the distinct identities of the release's run records) without changing the manifest's bytes for equal inputs (the manifest must stay a function of the data, not of who built it: say how you keep that). `bio-box releases` (card W-23) prints it.
Accept: tests first; the manifest hash of the fixture is unchanged.

### 0-10 Where a user's results live
**Decided 2026-11-12: option (b), one archive per run on `/stor` (DEC-33; card W-28).** · raised 2026-11-11 by the agent guide
A search's outputs are small files (tables of kilobytes to megabytes). The pod's local `/tmp` is cleared on reboot; `/stor` refuses files under 64 MB except registry tables marked `small_ok` (AGENTS §7.3). So a user's results have no durable home on the pod today. Options: (a) results are copied off the pod by the user and nothing durable is kept there; (b) a results area on `/stor` where each run is one archive (a tar per run or per day), written by a command; (c) an exemption for a results directory, agreed with the storage's owner. Recommendation: (b).

### C-25 S2C phase 2: the rest of the ≤20-per-species complete genomes
**M** (selection and fetch) then **H** (how to annotate) · Depends: C-21a, C-21b, C-21c (DEC-27), and the phase-1 build being done · confirmed by the PI 2026-11-11
Why: `STAGES.md` defines S2C as at most 20 complete genomes per GTDB species (25,458 genomes, 106.4 Gbp, *measured* from GTDB R226 metadata). Phase 1 built one per species (11,973 selected). Phase 2 is the remainder: about 13,485 genomes and 55.9 Gbp (*derived*).
Do: (1) Reproduce the plan's selection from the GTDB R226 metadata on the pod (complete genomes, at most 20 per species, by CheckM2 quality then ANI diversity): the ANI-diversity step needs skani, which is not in the image; say what the selection is without it (quality order only) and how many genomes differ, and stop for the PI if the tool is needed. Check that every phase-1 genome is in the set; report those that are not. (2) Fetch the remainder through the listed NCBI Datasets source (`sources.yaml`), verified and archived as phase 1 was, in stages of about 10 Gbp. (3) Before building: report, by content digest, how many contigs each new stage shares with stages already stored and within itself (the many-genomes case), and the fraction of proteins that would be new at NR90, estimated on one stage. (4) The annotation choice is the PI's (card 0-09): full Pfam on about 56 Gbp is about 2,500 CPU-hours (*derived*: 443 per 10 Gbp); representatives-only would cut it by the redundancy measured in (3).
Step (1) done 2026-11-12 (`planning/measurements/2026-11-12_s2c_selection.md`): the quality-only selection equals the plan's counts exactly (25,458 genomes, 106.4 Gbp; remainder 13,485 and 55.9 Gbp); the ANI step could change at most 3,368 selected genomes in 243 species, 485 of which were chosen by accession inside quality ties. Step (1b) done 2026-11-12: GTDB publishes no pairwise ANI, so DEC-31's fallback applies (card C-27: skani). Decided 2026-11-12 (DEC-31, DEC-32): step (1b), before any fetch: find out whether GTDB R226 publishes ANI values that cover every within-species pair among the complete genomes of the 243 species with more than 20 (name the file, its size and what it holds; propose the `sources.yaml` entry with its sha256; fetch it only after the entry is in the file); if it covers them, select 20 per species for ANI diversity starting from the best-quality genome (state the algorithm: greedy farthest-first on ANI, ties by accession), report how many genomes change against the quality-only list, and check phase 1 stays a subset; if it does not, stop and report what skani would need. Step (2) then fetches that list, retrying the 7 genomes NCBI did not serve once. Step (3) measures on one fetched stage the saving and the wrong-pair rate of representatives-only Pfam (DEC-32; stop above 1% wrong). Step (4), the annotation choice, is decided: DEC-32 (card B-47 implements it).
Accept: the selection file and its comparison with phase 1; the fetch record; the overlap table; then a build card once the PI chooses the annotation.

### C-24 Build stage H1, the huge phages
**M** · Pod: yes · Depends: C-20, and the 50 Gbp build's stage steps being done (one build at a time in the root)
Inputs: `/tmp/if2453/s2/inputs/H1/` on the pod (copy in `all-the-data/provisional/s2_inputs/s2_inputs.tar`); record `planning/measurements/2026-11-11_s2_inputs.md`. 7,228 genomes after the input script drops 2 byte-identical repeats; 516 records are already contigs of stage S1, so under DEC-26 their genes stay there and H1 adds 6,714 contigs and 1.82 Gbp (*measured* and *derived* in the record).
The PI confirmed the inputs as written on 2026-11-11: topology linear (the file records none; a gene across the start of a circular genome is split), assembly quality `complete`, and the 43 records a later preprint calls non-viral stay in (no list of them is available).
Do: `build-stage` for H1 into the root that holds S1 and P1 to P6, without `--library` first, then with Pfam (about 0.6 M proteins by the phage stage's density, *estimate*: a few hours of scan); then a release that adds H1. Check the gene-call result's skipped contigs against the 516 (514 distinct sequences); report genes, proteins, how many proteins are new at NR90 and NR50 against the existing release, and the genetic-code question: these are called in the phage mode with code 11, and jumbo phages with alternative codes will show as short broken genes: report the fraction of genomes whose mean gene length is under half the stage's median, as a flag for a later card.
Accept: the release verifies; the measurement file.

### C-19 An older release fails verification once a library is scanned on its stages
**M** · found 2026-11-11 at the end of the 10 Gbp rebuild · check first whether current code still does it
Evidence (*measured*, wilkcomp02, tree `740d049`, `/tmp/if2453/r10/phase2/verify.json`): after Pfam was scanned on stages S1, P1 and P2, releases `rN.1a`, `rN.1b` and `rN.1` (built before the scan, so without domains) each report `not in manifest: l1/domains/library=pfam-38.2/stage=<S>/part-0.parquet` and the matching ledger file, and `verify` exits 1. `rN.2`, built after, verifies with no problems. Nothing in the old releases changed: a later step added files to their stages, which the rules allow (AGENTS §3.10: a stage never rewrites an earlier stage's file; adding a library's files is not a rewrite).
Do: reproduce on the fixture with the current code (build a release without a library, scan the library on its stage, verify the first release). If it fails: a release's manifest scope is the files of the producers and libraries the release uses; files under its stages from a library or producer it does not list are not its concern and must not fail it (card B-33 defined the scope per release; extend it). Keep failing on a file the release does list that is missing or changed, and on an unlisted file under a path the release does own.
Accept: the reproduction fails first and passes after; the existing scope tests of B-33 pass unchanged.

### B-42 Finish split clustering: the real release, then end the exemption
**M** · Depends: B-41, the rebuild's phase 2 (nothing may write in the root meanwhile) · Pod: yes
Do: (1) on the pod, recluster `rN.1b` (2,053,075 proteins) in parts with the current code into a scratch root that reads the stage files without changing them (never the rebuild's root as the output), and compare with its published clusters row for row; every run record under 300 CPU-minutes; record wall time, CPU-minutes, the work directory's peak size on disk and each step's peak memory. (2) If equal: remove the DEC-25 exemption from `stage.py`, the profile key and the runbook in a separate commit, and propose the text of a decision row saying it has ended (the PI adds the row). If not equal: stop and report which members move.
Amended 2026-11-12 (DEC-29): step (2) narrows the exemption instead of removing it. The one-process clustering step of DEC-25 ends; `cluster_cpu_min_max` and its gate stay for the unsplittable calls of DEC-29 (card B-45). Do (2) after B-45 is merged, and change the runbook's text on DEC-25 to DEC-29.
Accept: the measurement file; the equality; then the removal.

### B-47 Pfam by representatives, with inferred rows (DEC-32)
**H** (design first), then **M** · Depends: C-25 step (3) · DEC-32
Design questions, in a note before code: where inferred rows live (the same `l1_domains` table with a nullable appended column that marks them and names the representative, or a table of their own), so that readers that need measured domains can filter and family tokens can include them; how envelopes are mapped through the stored alignment and what is written when a member's alignment does not cover the domain; how the scan plans its shards from the cluster table (clustering must run before the scan of that stage's proteins: say what that does to `build-stage`'s step order and to stages added later, whose new members join old clusters only at the next release); how `release verify` and the ledger count inferred rows; what a later full scan of the same proteins does to them. Measurement 0-09b is the evidence; card C-25 step (3) gives the rates on a strain-rich stage.

### C-27 Select S2C phase 2 for ANI diversity with skani (DEC-31's fallback)
**M** · Pod: yes · Depends: C-25 steps (1) and (1b) · DEC-31
Files: `container/environment.yml` (one line: skani, pinned; DEC-31 names it), `scripts/stress/build_stage_pod/c25_select.py` (`--ani`), a fetch and a skani driver beside it, `sources.yaml` if the fetch needs a new entry, `tests/r1/test_c25_select.py`, a measurement note
Why (card C-25 step 1b, `planning/measurements/2026-11-12_s2c_ani.md`): GTDB R226 publishes no pairwise ANI, so the 20 genomes kept in each of 243 species need ANI computed here. 28,483 genomes (118.8 Gbp) must be downloaded to choose at most 4,617; about 34.5 GB compressed and 120.3 GB of FASTA (*derived*); greedy farthest-first needs 499,624 pairs.
Do: (1) pin skani in `container/environment.yml` (the current bioconda release, exact version; say which and why) and nothing else there; the image is rebuilt when `main` next takes `plan/task-cards`, and skani runs only from that image (AGENTS §7.1). (2) Selection code first, with tests: `c25_select.py --ani <table>` takes pairwise ANI (and alignment fraction) per species and picks 20 by greedy farthest-first from the best-quality genome, ties by quality order then accession; pairs skani does not report (below its screen) count as the largest distance; a brute-force oracle written separately; input-order independence. (3) Fetch the 28,483 genomes through the listed NCBI Datasets source, verified as phase 1 was (retry the phase-1 genomes NCBI did not serve once, DEC-31), into scratch on the pod's local disk, in batches so that FASTA is removed once sketched (state the peak disk; stop if `/tmp` would fall under 200 GB free); archive nothing on `/stor` except, as one file over 64 MB, the genomes phase 2 will keep (after the selection). (4) When the image with skani is on the pod: measure skani on one species first (sketch time per genome, time per pair, memory) against the 4 to 17 CPU-hour *estimate*; then sketch all, compute the pairs farthest-first needs (or all pairs within a species where that is cheaper to state: say which), write the table, select. Report: genomes that differ from the quality-only list (and by species for the ten largest), the distribution of the minimum pairwise ANI among the 20 picked per species under each rule (the evidence that the rule matters), phase 1 a subset, the new remainder's genomes and bases, list hashes.
Accept: the tests; the measurement note; the phase-2 list ready for the fetch of the rest.

### C-26 Build the PLSDB plasmid stage
**M** · Pod: yes · Depends: C-24 (one build at a time in the root), C-21a to C-21c (merged) · approved by the PI 2026-11-12
Inputs: PLSDB FASTA and metadata, fetched and archived (`sources.yaml`; `planning/measurements/2026-11-11_s2_inputs.md`: 72,556 records, 7,139 of them already contigs of P1 to P6, *measured*).
Do: first the inputs and a short record, before any build: one genome per PLSDB record, kind `plasmid`, topology from PLSDB's metadata (say which field, and how many records it calls circular or linear or leaves empty), genetic code 11, quality from the metadata; no record dropped for repeating a sequence (DEC-27: a record whose bytes are already a contig of P1 to P6 becomes a second genome of that contig, and the earlier stage keeps its calls); report how many records share a contig with an earlier stage and with each other, and the bases and expected proteins the stage adds. Then `build-stage` without a library, then with Pfam, then a release that adds the stage: check that counts by genome include the shared plasmids' second genomes (`q1_genomes`, `bio-box releases`) and that counts by locus do not fan out; report genes, proteins, how many are new at NR90 and NR50, and time and peak memory per step.
Accept: the release verifies; the measurement file.

### W-28 Archive a run's results on `/stor` (DEC-33)
**M** · Depends: W-03, W-22 · DEC-33
Files: `src/biobox/cli.py` (`bio-box results archive RUN`, `bio-box results list`), `src/biobox/search/` (a small module), `profiles/wilkcomp02.yaml` (a results area, if `profile.py` allows a key without change; else say what it needs), `RUNBOOK.md`, `GUIDE.md`, `AGENT_GUIDE.md`, tests
Do, test first: one archive per run (tar, deterministic: sorted names, fixed mtimes and owners, so the same run gives the same bytes), written under a temporary name in the results area and renamed into place; a manifest line per archive (run id, release, sha256, size, date of the run from its record, not of the archiving); never overwrites: an existing archive with the same sha256 is a no-op, with another sha256 a refusal; refuses a results area that is not on the profile's durable storage; `results list` reads the manifests. The 64 MB rule (AGENTS §7.3) applies to the archive too: say what is done with a run smaller than that (the PI's rule covers files under 64 MB on `/stor`; propose bundling by day, or ask) and stop for the PI if a small archive cannot be avoided.
Accept: the tests; one real run archived on the pod and read back equal.

### E-10 A dedicated tRNA caller for phage and archaeal tRNAs
**H** (the PI adds `trnascan-se` to `container/environment.yml`; approved in principle 2026-11-12, as a later card), then **M** · Depends: 0-08b
Why (card 0-08b): Infernal's filters at our fixed Z miss 6 of 239 tRNA loci on six genomes, 5 of them in the jumbo phage LAK-A1, and cut three more short (two 91-base archaeal tRNAs); turning the filters off costs about 7,900 CPU-hours per 10 Gbp (*derived*).
Do (M, after the pin): measure tRNAscan-SE on the same eight genomes against the two Infernal scans (which of the 6 loci it finds, what it adds, isotype and intron calls, CPU per 10 Gbp); then, if it earns its place, a producer beside the ncRNA producer whose tRNA calls replace RF00005's where they overlap, by a stated rule.

### W-27 The synteny figure of `nbr flanks`
**M** · Depends: W-25 (merged) · approved by the PI 2026-11-12: add matplotlib (`uv add matplotlib`, and one line in `container/environment.yml`, which the PI named for this card; the image is rebuilt at the next merge to `main`)
Files: `src/biobox/search/flanks.py`, `src/biobox/cli.py` (`nbr flanks --figure`), `pyproject.toml` and `uv.lock` (`uv add matplotlib`, once the PI approves; `container/environment.yml` is the PI's), `GUIDE.md`, `AGENT_GUIDE.md`, tests
Why: card W-25 writes every table the figure needs (each leaf's order in the tree, the tokens per position, context types, metadata) and no figure, because the plotting library is not a locked dependency.
Do, test first: one figure per group (SVG and PDF, fixed metadata so the bytes carry no date): genes as arrows coloured by family, families below a set share in grey, rows in tree order with the tree beside them, columns for context type, genome kind, stage, taxonomy and the anchor's NR90; a bound on rows per page (say what happens above it). Byte-identical output across runs, thread counts and hash seeds. Also check that the GFF export loads in LoVis4u and clinker if either can be run off the pod, and say so if not.
Accept: the tests; the figure of the T7-like group of `rN.3`'s RNA polymerase search in `planning/measurements/`.

### E-09 Freeze `prok-2` as the ncRNA model set (DEC-30)
**M** · Depends: E-07, E-08, 0-08a (merged) · DEC-28, DEC-30
Files: `src/biobox/l1/rfam_sets/prok-2.tsv` (new), `src/biobox/l1/ncrna.py` (`DEFAULT_SET`, `FROZEN_SETS`, docstrings), `src/biobox/stage.py` and `src/biobox/cli.py` only where they name the default set or its shard rate, `RUNBOOK.md`, tests
Do, test first: `prok-2.tsv` is `planning/design/0-08a_prok-2_proposed.tsv` with a header that says frozen (770 models, 34 whole genome; check every accession is in Rfam 15.1 and that `prok-1`'s 588 rows are in it with their flags unchanged, by a test that reads both files); its sha256 goes into `FROZEN_SETS`; `DEFAULT_SET` becomes `prok-2`, so the producer name carries `prok2`; `prok-1` stays a frozen set that still loads and gives its old producer name (a stage scanned with it is still read). The provenance check accepts the new producer name. The shard planner's rate: the set is 18% dearer (*derived*, 299.7 against 254.5 CPU-hours per 10 Gbp), so shards cut for `prok-1` would run 18% longer: scale the planning rate by the set (state the formula) so a shard stays inside 80% of the budget; do not change `prok-1`'s. Then measure on the pod (Pod: yes, read-only on the root, into scratch `/tmp/if2453/e09`): the producer with `prok-2` on the 8 genomes of measurement 0-09c, as cards E-07 and E-08 did: hits, families, clan losers, CPU-hours per 10 Gbp against the derived 299.7 (stop and report if it misses by more than 20%), and that every `prok-1` feature is still a feature (or say which clan competition now takes).
Corrected 2026-11-12 (the E-09 agent stopped with evidence): "flags unchanged" above was wrong. DEC-30's set is the proposed file as it is: `prok-1`'s 588 models are all in it, and 4 of them (RF01670 P17, RF02111 IS009, RF02820 V_AS9, RF02908 rli30) change from intergenic to whole genome by the sample rule, with 7 new models (34 whole genome = 23 + 4 + 7). The test asserts that: every `prok-1` row is in `prok-2`, with its flag unchanged except those 4, which are 1. Those 4 are then searched on whole contigs without `--notrunc`: the measurement reports, for them, each `prok-1` feature that moves, changes or disappears. A set file that is not a frozen set is planned at the dearest frozen set's rate.
Accept: the tests; the measurement file.

### W-26 An accession lookup table, before the first stage of a million contigs
**M** · Depends: W-15 (measured) · not needed for `rM.2`; do it with the first draft-assembly or metagenome stage
Files: `src/biobox/registry/schemas.py` (one lookup table, then the lock by `scripts/update_schema_lock.py`), the builder that writes `contig_label`, `src/biobox/cli.py` (`find`, `search hmm --positive`, `seq --contig`), `GUIDE.md`, tests
Why (card W-15, `planning/measurements/2026-11-12_positive_lookup_rN3.md`): resolving an accession scans every label row of the release. On `rN.3` (39,633 rows, one row group per stage) that is 0.034 to 0.047 s (*measured*); it crosses one second at about 1 M label rows and would be about 110 s at S6's 111.5 M contigs (*derived*, linear, about 1 µs per row at 2 threads). Also: since card B-27 `contig_label` holds only each contig's first record accession, while `contig_occurrence` holds every record's; on `rN.3` a second accession of a contig (for example `MK327939` beside `NC_054918`) resolves only because its stages predate B-27 and the lookup falls back to `contig_occurrence` (*inferred* from the B-27 definition, not run on a B-27 release).
Do, test first: reproduce the second point on a fixture built by the current code (a contig with two record accessions: the second must resolve). Then a per-stage table sorted by accession holding every record's accession and its `contig_uid`, written by the step that writes `contig_label`; `find --accession`, `--positive` and `seq --contig ACCESSION` resolve through it with O(1) row groups per stage (`profile_query`, on the 100× store split into 20 stages); the header-substring form stays the documented slow path; stages built before the card answer as today with one note. Accept a protein accession where the release records one.
Accept: the tests; the pruning check; `GUIDE.md` §9 updated.

### W-25 Group a search's loci by their flanks, with a tree and a synteny figure
**M** · approved by the PI 2026-11-12; merged; its defaults (`--level nr50`, `--group-by anchor-nr90`, `--window 5g`, `--cut 0.5`, `--distance j`) and the rule that an unannotated gene is a token of its own were accepted by the PI on 2026-11-12 · Depends: W-23, C-21b, C-21c, W-24 (merged) · source: Thorn, Aarestrup and Munk, "Flankophile: a bioinformatic pipeline for prokaryotic genomic synteny analysis", Microbiol Spectr 12(1):e02413-23 (2023), PMC10783016, CC BY 4.0
Why: Flankophile finds a reference gene set in assemblies (ABRicate, 99% identity and coverage), keeps hits with at least 1,500 bp of flank on each side, groups hits by their closest reference gene, computes distances between flanks (KMA k-mer Jaccard on nucleotides), builds UPGMA trees, annotates flanks (Prokka) and draws synteny plots with metadata. A search here already gives the anchors (`search hmm`), their gene neighbourhoods (`nbr`, `loci/*.tsv`, tapes) and family and cluster names (Pfam, NR50, NR90), but no way to say which loci share a context. Doing it on tokens instead of nucleotides needs no new tool and reads only the tapes the release has. Method: distance on flank tokens as in Domainator (Jaccard and adjacency index, from BiG-SCAPE) and UPGMA as in Flankophile; review in `planning/design/W-25_flank_tools_review.md` (amended 2026-11-12 from it).
Do (design note first, half a page, then test first): `bio-box nbr flanks RUN --window Ng --level nr50|nr90|pfam [--group-by anchor-nr90|anchor-nr50|none] [--max-loci K] [--cut D] [--export gff]`. Each locus is the anchor's flanks as two ordered token lists (upstream, downstream) oriented to the anchor's strand. Each token is (level id, side, strand relative to the anchor). At `pfam` the id is the protein's Pfam architecture. A gene without an id at the level gets a token unique to it (never shared). Censoring by a linear contig end is stated per locus. Distances, both reported: J = 1 − Jaccard of the two loci's token sets (default), and A = 1 − adjacency index (the Jaccard of adjacent oriented token pairs, BiG-SCAPE). The combined distance 1 − (w·J_sim + (1 − w)·A_sim) is offered with w stated in the note and chosen on the planted test, not copied from BiG-SCAPE. Where a contig end cuts a flank, a pair is compared only over the positions both loci observe on that side, and the number of positions compared is written per pair. Loci are grouped by the anchor's own cluster. Within a group, loci with identical token lists collapse to one leaf with a count before any pairwise step. Per group: the distance matrices (TSV, leaves named by `feature_uid`), a UPGMA tree (Newick, inputs sorted by `feature_uid`), context types from cutting the tree at D (one row per locus: `feature_uid`, group, `context_type`), and a conservation table per context type and neighbour token (share of loci, modal offset along the anchor's strand, share on the same strand, distinct genomes through `genome_member`). Also one figure (SVG and PDF, matplotlib; no R): genes as arrows coloured by family, families below a set share in grey, rows in tree order, columns for context type, genome kind, stage, taxonomy and the anchor's NR90. Optional `--export gff`: one GFF3 per locus with `protein_uid` and the token as attributes, plus a protein-to-group TSV, readable by LoVis4u, clinker and Domainator. Bounded: at most K leaves per group after collapsing (default 2,000; above it, a deterministic sample by `feature_uid`, said in the output). Deterministic output (AGENTS §4.3). Every row names `feature_uid` beside any ordinal (AGENTS §3.3). The method is reimplemented from the papers; no code is copied (BiG-SCAPE is AGPL-3.0; the others state no code licence).
Accept: tests on a constructed release with planted contexts. (i) Two flank types around one anchor family, one shared by phages and prokaryotes, are recovered as two context types. (ii) A pair with the same families in a rearranged order is separated by A and not by J. (iii) One inserted gene changes A by at most the pairs it touches. (iv) A censored locus is compared only over the positions it observes. (v) Unannotated neighbours do not make two loci similar. The oracle for J and A is a brute-force loop written separately, and Domainator `compare_contigs` on the export once the image has it. On `rN.3`, the RNA polymerase run (2,581 loci for 1,957 proteins, `GUIDE.md` §6) is grouped with time and peak memory recorded (Pod: yes, read-only).
Not proposed: running Flankophile itself (its tools, Prokka, ABRicate, KMA, CD-HIT and R, are not in the image; `container/environment.yml` is the PI's), importing its evaluation data (2,006 Danish clinical isolates, PRJEB37711, and 273 pig gut metagenomes: strain-redundant clinical species, and the metagenomes are not public as assemblies), or ResFinder as a nucleotide reference set (a source entry for the PI if wanted).

### B-46 The tokens step runs out of memory at 46 M proteins
**M** · blocks `rM.2` (the 50 Gbp release) · found 2026-11-12 by the 50 Gbp build · Pod: yes (read-only on the root)
Files: `src/biobox/l3/tokens.py`, tests; and, only if the same failure is found there, `src/biobox/l3/postings.py` and `src/biobox/dim.py` (say so in the report)
Evidence (*measured*, wilkcomp02, 2026-11-12, unit `r50-P6-5`, log `/tmp/if2453/r50/tmp/logs/r50-P6-5.log`): release `rM.1` (stages S1, P1 to P6; 45,815,843 proteins; 31,183,876 NR90 and 8,192,065 NR50 clusters) published its clustering at 10:53 CST; `rM.1.tokens` failed 4 minutes later in the first `COPY` of `tokens.run` (the slot join: tapes unnested, left-joined to `fam` and to `cl` twice) with `_duckdb.OutOfMemoryException: failed to allocate data of size 32.0 MiB (1.8 GiB/1.8 GiB used)`. The connection has `threads 1`, `memory_limit DUCK_MEMORY_LIMIT`, `preserve_insertion_order True`. The same code built `rN.3` (10.6 M proteins). Nothing was published. The build now scans Pfam on P3 to P6 and waits for this card before `rM.2` (`scripts/stress/build_stage_pod/r50_driver.sh`, `SKIP_RM1`).
Do: reproduce first, at a size that fails with the current code: on the 100× synthetic store with a lowered memory limit if that reproduces the same operator's failure, and on the pod against `rM.1`'s real inputs, read-only, writing to a scratch `out=` under `/tmp/if2453/b46` (never into `/tmp/if2453/r10/root`; bind it read-only). Find which operator cannot spill (the `arg_max` aggregates of `cl` and `fam` over BLOB keys, the order-preserving `COPY`, or the joins) and change the query so memory stays bounded by DuckDB's limit at any size (AGENTS §4.1): for example materialise `cl` and `fam` as sorted Parquet first, or join in ranges of `protein_uid`. Output byte-identical to today's on the fixture and on `rN.3`'s inputs (compare with `rN.3`'s published `tape_tokens` and `token_dict`, read-only). Then run the release steps after tokens (postings, dim, manifest, catalog) of `rM.1` into the same scratch `out=`, to find the next failure before the build meets it; report each step's peak memory, wall time and CPU-minutes. Pod threads: at most 8, inside the account's total (a Pfam scan of 32 threads is running); `MemoryMax=32G` scopes, nice 19.
Accept: the reproduction fails first; byte-identity on the fixture and `rN.3`; `rM.1`'s tokens step completes into scratch with its peak memory recorded; the later steps' results reported.

### B-45 Run an unsplittable clustering call as one process (DEC-29)
**M** · Depends: B-41, B-43, C-22 (merged) · DEC-29 · before any release above about 74 M proteins
Files: `src/biobox/l2/cluster.py`, `src/biobox/stage.py`, `RUNBOOK.md`, tests
Why: the cut refuses a prefilter call whose per-part table rebuild alone exceeds the part budget (`cluster.py`, the C-22 table: the NR90 step's prefilter from about 74 M proteins, NR50 step 0 from about 108 M, *derived*). DEC-29 runs that call unsplit instead.
Do, test first: when the cut would refuse a call for that reason and the profile sets `cluster_cpu_min_max`, plan the call as one part (its own table, `--add-self-matches` as today, the log check for `Process prefiltering step 1 of 1` unchanged), with its estimate (work plus table) checked against `cluster_cpu_min_max` before launch and refused above it, and run it with `budget_reason` `exempt: DEC-29`; without `cluster_cpu_min_max` keep today's refusal and name DEC-29 in it. Every other call is cut as today. The one-part call's output equals the workflow's (the 1-part equality test of B-41 covers it; add the case to it with a forced overhead over the budget on the fixture and the constructed set). `release cluster plan` shows which calls run as one process and why. Restate the C-22 table's 95 M row with the call as one process, and the total.
Amended 2026-11-12 (the B-45 agent stopped with evidence: the cut runs in its own process under `run-step`, which carries no profile, and the plan's `cut` record holds only `parts` or `part_cpu_min`): Files gain `src/biobox/cli.py`, for one option on `release cluster plan` only (the one-process cap), which `stage.py` passes when the profile sets `cluster_cpu_min_max` and the plan records in `cut`, so the cut stays a function of the plan; a forced overhead for tests is likewise an explicit plan parameter recorded in `cut`. Also: when clustering runs split, the release-level DEC-25 check (`stage.cluster_refusal`, which compares the one-process estimate T(n) with the cap and refuses above about 62.7 M proteins) checks instead the estimate of each call the plan runs as one process against the cap; the unsplit one-process mode keeps today's check until card B-42 (2).
Accept: the tests fail first; equality at 1, 2 and 4 parts unchanged; with the profile key absent the refusal is unchanged except its text.

### B-43 A shared prefilter index for split clustering
**M** · Depends: B-41 · needed above about 68 M proteins
Evidence (*measured*, wilkcomp02, MMseqs2 18.8cc5c, 2026-11-11, B-41 report, `/tmp/if2453/b41/idx.sh`): `prefilter` uses a prebuilt index only when the index itself is given as the target (`db.idx`); it then equals the workflow's output on the sensitive NR50 calls (`-s 2`, `-s 3`, with `--comp-bias-corr 1 --diag-score 1`) at 0.5 to 0.8 CPU-s against 5.2, but aborts with `munmap_chunk(): invalid pointer` on the `-s 1 --comp-bias-corr 0 --diag-score 0` calls (the NR90 step and NR50 step 0). Without a shared index each part rebuilds it: at 95 M proteins that alone is about 333 CPU-minutes per part (*derived*), over the budget, and at 50 M the largest prefilter call costs about 3.7 times its unsplit CPU (*derived*).
Do: use the shared index on the calls where it works, with a test that the result equals the rebuilt one; for the calls where it crashes, find the cause (the parameters the index was built with against those the call asks for; try an index built with the call's own settings) or, if it is a defect of this MMseqs2 version, say which version fixes it and propose the pin change for the PI (`container/environment.yml` is the PI's to change). State the part sizes and total CPU at 50 M and 95 M with whatever works.
Accept: equality tests on the fixture and the constructed set; the derived budget table.

### B-44 Clustering above about 87 M proteins needs an MMseqs2 without the index abort
**H** (PI) · deferred by the PI on 2026-11-11: do nothing now · Depends: B-43
Evidence (*measured*, wilkcomp02, MMseqs2 18.8cc5c, card B-43): `prefilter` aborts (`munmap_chunk(): invalid pointer`, exit 134) whenever its target is a prebuilt index and `--diag-score` is 0, at every index and run setting tried (sixteen variations). Reproducer: `createdb in.faa db; createindex db tmp; prefilter db db.idx pf --diag-score 0`. The fast clustering calls (the NR90 step and NR50 step 0) use that setting, so their parts each rebuild the index; the cut refuses the NR90 step from about 87.2 M proteins and NR50 step 0 from about 108.1 M (*derived*, card B-43). 18-8cc5c (2025-07-27) is the latest release on GitHub and bioconda as of 2026-11-11; no upstream issue describes the abort.
When it is needed (a release over about 87 M proteins): (a) build the current upstream commit in a scratch image and run the reproducer; (b) report the abort upstream with the reproducer (from the PI's account, text approved by the PI); (c) pin a fixed version only if it passes the reproducer, the module-list test, the 1, 2 and 4-part tests with those calls added to `INDEX_SETTINGS`, and gives clusters equal row for row to today's on the fixture and on `rN.1b` (DEC-24: otherwise cluster ids change between releases and that is its own decision). `container/environment.yml` is the PI's to change.

### W-21 Read the next target block while searching
**L** · Depends: W-19
Why: after W-19 one thread reads each block of the target while the search threads wait: about 2 of 4 CPUs were busy on `rN.1` and the search passes took about 234 s against 121 s before (*measured*, W-19 report).
Do: prefetch the next block in a background thread (one block ahead; memory bound: two blocks); results identical. Also make `worker._families` stream `l1_domains` instead of loading every row when `foreign_domain_alarm` is on (it would grow with a Pfam-annotated release).
Accept: the W-19 memory test still passes with its bound restated; on the fixture, hits unchanged; on the pod the `rN.1` search time recorded against 294 s.

### C-09 The file-size rule does not refuse small release records
**L** · Depends: C-08 (found in its report) · blocks the second release under C-04
Files: `src/biobox/contracts.py` (`check_file_sizes`), `src/biobox/registry/schemas.py` (a `small_ok` flag on `ArtifactSpec` if it lacks one), `tests/r2/test_micro_s0_contracts.py`, `tests/registry/operations.yaml` (drop `catalog` and `release_manifest` from `external` if R0 accepts it)
Why: `check_file_sizes` groups artifacts by (artifact, stage) and refuses a group with more than one file under the profile's floor (64 MiB on wilkcomp02). The release manifest and the catalog have no partitions, so every release's manifest falls in one group: the second release's manifest and catalog would be refused by `run-step` (`_check_sizes`). The same grouping would refuse two small domain-library files in one stage. Derived from the code by the C-08 agent; not yet run at the pod's floor.
Do: reproduce first at the real floor (a test profile with `small_file_floor_bytes: 67108864` and two releases built through `run-step`). Then: artifacts and tables that are small by nature (manifest, catalog, `store_params`, and every `small_ok` table) are exempt by a registry flag, not by name in the check; for the rest, group by the full partition tuple, not by stage alone. Do not lower or bypass the floor for data tables.
Accept: the reproduction fails first and passes after; a deliberately tiny second file in a data table's partition is still refused.

### B-30 Ingest by batch processes
**M** · Depends: B-22, C-08
Files: `src/biobox/l0seq/store.py`, `src/biobox/stage.py`, `src/biobox/registry/schemas.py` (one `ArtifactSpec`, `ingest_plan`, and nothing else), `src/biobox/cli.py` (`stage ingest plan|pass1|batch|merge`), tests
Amended 2026-11-10 (found by the B-30 agent before editing): `run-step` refuses an empty or unclaimed staging, so the plan and pass-1 steps need a registered record to publish. Add one artifact `ingest_plan` at `l0seq/ingest_plan/stage=*`, `in_manifest=False`, `small_ok=True`, as B-24 added `domain_scan`. The commands are `bio-box stage ingest plan|pass1|batch|merge`, matching `stage domains`. Pass 1 stays sharded by source object as written below; one shard is a valid plan for a small stage. Lost scratch is rescanned and checked against the published sha256, as C-12 does for domain scans.
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
Amended 2026-11-10 after the agent's stop (probe on the pod: the fixture pack is 7 files, 2.14 MB, byte-identical across two builds; the size rule refuses it at the 64 MiB floor, and would at any scale because the pack manifest is a few KB): the pack does **not** go into the data root and is not a `build-stage` step or a registered artifact. It is a rebuildable cache of a release in a pack directory outside the root (`--pack-dir`, default under the profile's `tmp`), at `<pack-dir>/release=<R>/<pack_id>/`, written once and never overwritten. A search whose pack is missing exits 2 and names the command that builds it. The sentence above about `search_pack/release=<R>/` under the root is superseded.
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

### W-13 A search launched on the pod can be watched and cancelled
**M** · Depends: C-06, A-07 (same files) · found while testing the runbook on the pod, 2026-11-10
Files: `src/biobox/search/runner.py`, `src/biobox/cli.py` (`search`), `scripts/pod_run.sh`, `RUNBOOK.md`, `tests/r0/test_runbook.py`, tests
Evidence (*measured*, wilkcomp02, units `c06-*`): (1) `search run` without `--wait`, launched through `pod_run.sh`, prints a run id; the unit then ends with its main process and systemd kills the detached runner, so the run ends `failed: runner exited without closing the run`. (2) A run started with `--wait` has no pending record while it runs, so `search status` finds nothing and `search cancel` cannot be used; the only cancel is `systemctl --user stop <unit>`. (3) After that stop the unit's log has no `pod_run: <unit> exited` line (the wrapper shell is killed too), and `run-step`/`build-stage` output reaches the log only when the unit exits (buffered). (4) The default search work directory is `/tmp/atd/search`, not under the profile's `tmp`.
Do: a run started with `--wait` writes the same pending record as a submitted one, so `status` and `cancel` work on it from another shell; `search run` without `--wait` under a unit either keeps the unit alive until the run closes or refuses with a message that says to use `--wait` there (decide by what `systemd-run --collect` allows; say which); `pod_run.sh` traps termination so the last log line is always written, and runs Python unbuffered; with a profile, the work directory defaults to a per-user directory under the profile's `tmp`. Update the runbook sections and its appendix rows for each.
Amended 2026-11-10 after W-03 merged: (a) the per-user default work directory is done (W-03: `<tmp>/search/<user>/work`, runs in `<tmp>/search/<user>/runs`, packs in `<tmp>/packs`); do not redo it. (b) The termination trap and unbuffered Python in `pod_run.sh` belong to card C-13; do them here only if C-13 has not merged. (c) W-03 moved the pack and the runs out of the data root and left `RUNBOOK.md` untouched: bring the runbook up to date here. `pack build` and `pack verify` need `--pack-dir` or `--profile` and a sentence that the pack lives outside the root; `search run`, `status`, `cancel`, `hits` and `nbr` need `--runs` and `--pack-dir`, or `--profile`; the work-directory note, the locations and `/tmp` cleaning notes (add `<tmp>/packs` and `<tmp>/search/<user>/`), the known-issue notes and the appendix rows change with them. (d) `catalog.EXCLUDED_TREES` still lists `search_pack` and `search` for old roots; correct its comment, which names `search/run_index.json` (add `src/biobox/catalog.py` to Files for that comment only).
Accept: tests first for the pending record and the default work directory; on the pod, a search launched through `pod_run.sh` is visible to `search status` while running and ends `stopped` after `search cancel` from another shell; the log of a stopped unit ends with the exit line.

### W-12 Hit and locus exports read only the manifest's files
**L** · Depends: B-28, W-11 (found in the W-11 report; AGENTS §3.9)
Files: `src/biobox/search/export.py`, `src/biobox/cli.py` (`search hits`, `nbr`), `tests/r2/test_micro_s0_search.py`
Do: `hits_tsv`, `hits_faa` and `locus_tables` still fall back to a registry glob when the release has no manifest. Make them refuse as `export proteins` does since W-11 (exit 2, naming the manifest and `bio-box release manifest`), or take the release's stage list through the readers B-28 changed; no path may glob over stages.
Accept: tests first: each command without a manifest exits 2 with that message; with a manifest and an extra stage on disk, output is unchanged.

### B-27 Lookups that cannot prune: contig labels and cluster membership
**M** · Depends: B-25 (found in the W-05 report; AGENTS §4.4)
Files: `src/biobox/registry/schemas.py` and lock (two tables), their builders (`src/biobox/l0seq/store.py`, `src/biobox/l2/cluster.py`), `src/biobox/catalog.py` (`contig_label`, macros), `src/biobox/search/export.py`, `src/biobox/search/nbr.py`, tests
Why: `contig_occurrence` is looked up by `contig_uid` but sorted by (`l0_object_sha256`, `record_ordinal`); `cluster_member` is looked up by `member_uid` but sorted by (`level`, `cluster_uid`, `member_uid`). Neither lookup can use min/max statistics, so exports and `nbr summarize` scan the whole table: 10.5 M `cluster_member` rows at `rL.2p`, 190 M at GTDB scale (*derived*: 2 levels x 95 M proteins).
Also (from the W-12 report): remove the dead `manifest=None` glob branch of `nbr._store_files`; nothing reaches it since W-12, and with it gone `rg -n "rglob|glob\(" src/biobox/search` returns nothing.
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

### E-07 ncRNA scan in shards, and in `build-stage`
**M** · Depends: E-01, C-12, PI card 0-08 (which hits, which models)
Files: `src/biobox/l1/ncrna.py`, `src/biobox/stage.py`, `src/biobox/cli.py` (`stage ncrna plan|shard|merge`, `build-stage --ncrna-cm`), tests, `RUNBOOK.md`
Why: one `cmsearch` process scans about 23.7 Mbp inside the 300 CPU-minute budget (*derived* from 1,319 bases per CPU-second, *measured*), so any real stage needs shards; and the producer is not in `build-stage`'s step list.
Do: shard the ncRNA scan as card B-24 did the domain scan (plan by cumulative bases over contigs sorted by `contig_uid`, one process per shard, merge in sort order, byte-identical to the unsharded run), reuse C-11's shard pool with a memory term (each `cmsearch` peaked at 12.3 GiB: jobs × 13 GiB must fit the unit's cap); fill the `features` slot of `build-stage` when `--ncrna-cm` is given; pass both producers to the IGR and tape steps, which changes the tape version by design.
Accept: sharded equals unsharded on the fixture; a `build-stage` run with and without `--ncrna-cm` gives two tape versions and the gene-call tables are byte-identical between them.
**Amended 2026-11-11 after DEC-28** (the design note `planning/design/0-08_ncrna_models.md` is the specification; read §2, §4 and §6 whole). Files gain `src/biobox/l1/rfam_sets/prok-1.tsv` (the frozen model set: accession, name, and whether the model is scanned on whole genomes; generated once from the note's rule script `planning/design/0-08_model_rules.py` and its tables, then never regenerated: a new set is a new file with a new version) and new test files. The producer scans (a) intergenic regions with a 200-base margin into the flanking genes, overlaps merged, with a map back to contig coordinates and strand and the origin of circular contigs handled, using `cmsearch --notrunc`, against the set's intergenic models; and (b) whole contigs against the set's 23 in-CDS models; never the same model on the same base twice. Hits are mapped back to contig coordinates (0-based half-open; a feature across the origin has `end > contig_len`). The producer name carries the set's version and the region rule, so a change to either is a new producer and a new tape version. A contig another stage owns is skipped (DEC-26, card C-20). Shards are cut by genome bases at 80% of the CPU budget from the measured rate (about 202 Mbp per process at 247.4 CPU-hours per 10 Gbp, *derived* in the note); memory per `cmsearch --cpu 2` was at most 6.4 GiB on intergenic records (*measured*, 0-09c), which replaces the 12.3 GiB above for the pool's memory term. Clan competition is card E-08, not this one: keep every hit, as today.
Tests first (each against an independent oracle): the regions of a fixture's gene calls equal a separately written per-base loop (margin, merging, linear ends, origin wrap); hits map back to the coordinates the `cmsearch` binary reports on the whole contig; the subset's hits equal the full fixture library's hits restricted to the subset; in-CDS models are scanned whole and nothing is scanned twice; sharded output is byte-identical to unsharded at 1, 2 and 4 shards.
Accept: the tests; then on the pod (Pod: yes) on the 8 genomes of measurement 0-09c the step reproduces that measurement's margin-200 `--notrunc` hits restricted to the set plus the whole-genome hits of the in-CDS class, at a CPU cost within 20% of 247.4 CPU-hours per 10 Gbp; and the two-tape-version check above.

### E-08 Clan competition for ncRNA hits
**M** · Depends: E-07 (same file) · DEC-28
Files: `src/biobox/l1/ncrna.py`, `src/biobox/l1/rfam_sets/Rfam-15.1.clanin` (in the repository; pinned by sha256 in `refdb.yaml`, entry `rfam_clanin`: check the hash at load), `src/biobox/registry/schemas.py` (`l1_ncrna`: append nullable `clan` and `clan_winner`), the schema locks by `scripts/update_schema_lock.py` only, tests
Do: as §5 of `planning/design/0-08_ncrna_models.md`: among overlapping hits of models in one clan, on the same strand, the best bit score wins; state the overlap rule exactly and where it differs from Rfam's own (`cmscan --oclan`), and test against that binary's output on a constructed case if the image's Infernal offers it. Every hit stays in `l1_ncrna` with its clan and whether it won; only winners (and hits of models in no clan) become `feature` rows and go on the tape. Hits of different clans that overlap are all kept (nested and antisense features are real).
Also in this card, from the E-07 report (2026-11-11): (a) `provenance.NCRNA_RE` does not accept the new producer name (`infernal-1.1.5-rfam-15.1-prok1-igr200nt`), so `check_provenance` on a release with ncRNA tapes would reject it: accept an optional set and region tail, test first (add `src/biobox/provenance.py` to Files for that). (b) Shard plans and records were put in the `domain_scan` artifact for lack of their own: since this card already changes `schemas.py`, add an `ncrna_scan` `ArtifactSpec` and move them, with releases built before it still readable. (c) The producer name does not say which gene-call version's intergenic regions were scanned: add it to the ledger row and the run record's result (not the name), and say so in the docstring. (d) Per-base cost differed 1.5 times between the two measured shards (206.6 and 304.0 CPU-hours per 10 Gbp, *derived*): cut shards from the slower rate, so an 80% shard stays under the budget, and restate the shard size.
Tests first, on constructed truth (AGENTS §5.6): the E-02 spike-in panel extended with a planted rRNA operon (several SSU and LSU models must give one kept feature each) and an antisense pair (both kept); a unit test against a separately written brute-force loop.
Accept: the tests; on the 8 genomes of measurement 0-09c, one kept SSU and one kept LSU feature per operon, and every other truth hit kept.

### 0-08a Which long-tail Rfam models hit on a sample of the corpus
**M** (measurement) · Pod: yes · DEC-28 · M1 of the design note
Do: scan the full Rfam 15.1 library on whole genomes of a 100-genome sample of the 2,497 prokaryotes of release `rN.2` (stratified by phylum from GTDB's metadata already on the pod, archaea included) plus 300 phages of stage S1, with the producer's settings. Report: families with hits as genomes are added; which long-tail models (rule d) hit and should join the set; any excluded model (rule c) that hits, which would be a wrong exclusion; any rule-(b) model named by description that never hits; hits lying inside genes, by family. About 81 CPU-hours (*derived* in the note): run it in processes under 300 CPU-minutes, at a thread count that fits beside the 50 Gbp build (check the account's use before each launch), as detached units.
Accept: `planning/measurements/<date>_rfam_sample.md` with every number labelled, and a proposed `prok-2.tsv` for the PI.

### 0-08b The ncRNA set as one scan, tRNAs lost to filters, and short contigs
**M** (measurement) · Pod: yes · Depends: E-07 for the set file (or build the list from the note's tables) · M2, M4, M5 of the design note
Do: (M4) run the 588 models at margin 200 with `--notrunc` and the 23 in-CDS models on whole genomes, on the 8 genomes of 0-09c, as one configuration, and compare cost and hits with the derived 247.4 CPU-hours per 10 Gbp and 99.4% recall. (M2) scan the tRNA and tRNA-Sec models alone with `--max` on the same genomes and count loci at or above GA that the filtered scan misses. (M5) the same configuration on about 50 Mbp of phage contigs from stage S1, since intergenic cost follows the number of records.
Accept: `planning/measurements/<date>_ncrna_set.md` with every number labelled; if the measured cost misses 247.4 by more than 20%, say so first.

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
Approved by the PI 2026-11-12: add `zombi2==0.54.0` with `uv add` (the one new dependency; `pyproject.toml` and `uv.lock` are in Files for that line only; if the resolution would change any existing pin, stop and report). `container/environment.yml` is not changed: until an image holds the package, pod tests import it from a directory of its own (never the shared `/tmp/if2453/pilot_s1/pydeps`), and a test that needs it skips where it is absent and runs on GitHub Actions.
Sessions: (1) the backbone (spec §2.1, §2.2, §6.1, the `micro` preset's backbone parameters of §5.1, the three DEC-15 conditions as tests, and the corrections of the spike's §6): gene order, strands, families, copies and groups with their truth, no sequences; later sessions take §2.3 to §2.5 (sequences, IGRs, RNA motifs), §3 (truth tables), §4 (null worlds), §6.2 to §6.4 and §7, one section each, each ending with the list of what remains.
Session 3 merged 2026-11-12 (families and coding sequences: profile-F81 on the gene trees, synthetic panel fitted to Pfam, GC3 codon model, code 11 only; output keyed by `copy_label`). Chosen values for the PI to confirm: F81 with site-specific frequencies; gamma alpha 1 in 4 categories; uncorrelated log-normal clock (sigma 0.3); kappa 3 for synonymous changes; tail and linker indel rate 0.05; GC3 by Brownian logit (sigma 0.25); start codons ATG/GTG/TTG 0.90/0.08/0.02, stops TAA/TGA/TAG 0.64/0.29/0.07; difficulty weights 0.5/0.35/0.15; domain counts 0.6/0.3/0.1. Open before session 6: detectability does not fall with divergence under this model (PI). Left over from §2.3: genetic codes 4 and 25, `real_seed_member` roots, named promiscuous domains, `ref_bitscore_median`, a Pfam panel usable in CI, a bounded `iter_genome`, and a way to bind the Pfam file in `pod_pytest.sh` for the parity tests. Session 2 merged 2026-11-12 (taxonomy, plasmids by excision, vehicles, hubs, decoys, distance decay, lineage-varying rates, typed events and gene trees; format `synthworld-backbone/2`); what remains starts at item (3) below, and the truth-table session must decide whether plasmid and vehicle pairs are excluded from null evaluation; N0 must drop decoys. Session 1 merged 2026-11-12 (backbone; output `out/backbone/{genomes,genes,families,groups,truth_units,truth_adjacency}.tsv` and `backbone.json`, documented in `backbone.py`). Remaining, in a workable order: (2) the rest of the backbone of §2.1 (plasmids as further replicons, vehicles, hubs, decoy recipients, distance-decay transfer, lineage-varying rates for the §5.1 size range, typed `gene_trees` and `genome_events` tables) and the synthetic taxonomy; (3) §2.3 sequences; (4) §2.4 and §2.5 IGRs and RNA motifs; (5) §2.1 S6 to S11 (phages and CRISPR, fusions, call errors, assembly, contamination, emit); (6) §3 truth tables; (7) §4 null worlds; (8) §6.2 to §6.4 and §7. Each session also applies the spike corrections (its §7, ten items) that touch its sections: items 1, 2, 3, 4, 6, 8 and 9 are still open. The `micro` strength grid (member transfer rate {0.05, 0.15, 0.30} by together fraction {0.3, 0.7, 1.0}) was accepted by the PI on 2026-11-12.

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
