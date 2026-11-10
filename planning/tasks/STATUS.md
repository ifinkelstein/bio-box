# Task status

The running record of the cards in [`TASKS.md`](TASKS.md). Updated at every merge. Branch
`plan/task-cards` holds all merged work; nothing is on `main` yet. Each card was done by an agent on
its own branch `task/<ID>`, audited by the orchestrating session, and merged after the full suite
passed on the integration branch.

Last update: 2026-11-09, 23:00 CST. Suite on the integration branch after the B-17 merge: 861 passed, 10 skipped; golden ring
(`--golden tests/r3`): 32 passed, 5 skipped
(*measured*, wilkcomp02, 2026-11-10 about 03:00 CST, `sh scripts/pod_pytest.sh integration --shards 8`, in the
pinned image; 3 min 40 s wall for the slowest shard, against 13.5 min in one process). The skips
need samtools, bgzip, Prodigal or the golden set. Tests no longer run on the laptop (AGENTS.md §5.5).

## Merged (49 cards)

| Card | Merge | What it did | Audit notes |
|---|---|---|---|
| A-01 | `14acede` | Search-runner liveness by a file lock instead of pid checks | Audit made the probe a shared lock and stopped `cancel` signalling a reused pid. First Linux CI run since: green |
| C-01 | `3bdacfa` | Host profile for wilkcomp02 | The 08:00 to 20:00 day window was my assumption, now labelled so; weekends use the night cap |
| C-02 | `85d4da8` | `bio-box run-step`: scope, nice, thread caps, staging, publish, run record | Audit found two gaps, fixed in C-03 |
| C-03 | `914320d` | Kill and resume for `run-step`; signals forwarded | Checked on the pod that a systemd scope keeps the pid and process group |
| A-02 | `386ed06` | Corrected stage sizes and rates in the plan documents | Audit removed the stored-IGR term from the store model (0.528 GB per Gbp) |
| A-05 | `fc14de0` | Spill-watcher thread race in the K-02 test | |
| B-01 | `106f387` | One streaming Parquet writer (`io.write_table`) | Bytes are independent of input batching |
| B-02 | `2646ace` | `feature_uid` and feature kinds | |
| B-03 | `2f6f942` | Six new registry tables | |
| B-04 | `315de8e` | Persisted genome inputs | |
| B-05 | `3b55cc0` | Gene calling: streaming, parallel, writes feature tables | 7× faster at 8 workers; memory grows 87 MB from 30 to 300 genomes |
| B-06 | `d12fcf7` | IGRs keyed on features | |
| B-07 | `452d573` | Tapes keyed by tape version; slot-bounded row groups; pruning lookup | Audit: one ledger per tape version |
| B-08 | `65854d1` | Stability relation MR-20 | 0 rows change by `feature_uid`; 514 of 527 by ordinal |
| B-09 | `1bac11e` | Token building through DuckDB | Stopped once for usage limits and resumed; its memory measurement was never run |
| B-12 | `c0019e5` | `contig_dim` | |
| B-20 | `f6a1175` | IGR version names its producer set | |
| B-10 | `8b27396` | `token_postings` builder | Narrowed to one posting per slot and level; per-domain postings are card B-21 |
| B-13 | `6e741b4` | Sequence index looked up, not loaded; `.fai` bisected | Open at 10⁶ contigs: 0.001 s and 0.008 MB |
| W-01 | `e7ed141` | Pack build survives genetic code 15 | |
| W-08 | `6b5c0f2` | Recipe shorthand with label; drift-alarm warning; `contig_label` view | |
| B-14 | `258dc11` | Sequence macros without index scans | 1,000-contig fetch at 10⁶ contigs: 3.2 s to 0.25 s |
| B-15 | `f9af5b6` | Batch ingest with bounded memory | 28% slower (sources read twice); found the old ledger depended on input order |
| B-11 | `f331c20` | Neighbor queries over postings; `slot_pairs` removed | Equal to the pair oracle on all 7,735 fixture pairs |
| W-02 | `db24106` | A user's HMM file or seed FASTA as a search seed | The HMMER-binary oracle test has not run (image only) |
| D-01 | `f682baf` | Mutation testing with a two-pass verdict | Workflow written, not executed |
| A-04 | `4ad0e82` | `doctor`: cgroup memory cap and path visibility | |
| B-16 | `c25699d` | Search layer names slots by `feature_uid` | |
| D-03 | `b108be6` | samtools, bgzip, Prodigal and RefSeq checks as golden-ring tests | Only the RefSeq test has executed; Prodigal tolerance set from the measured 4 |
| B-19 | `4972f71` | Streaming ledger; builder tables cut into parts | |
| W-06 | `5f31050` | `nbr summarize`: neighbor families and clusters per anchor window | Same commit adds `scripts/pod_pytest.sh` and fixes a cluster-ledger test that only fails where MMseqs2 exists |
| C-05, W-04 | `0a4ae15` | Search workers start in scopes; worker thread caps; exit and kill labelling | The scope cannot be opened from inside the image: card C-07 |
| W-09 | this merge | Every pyhmmer search in `search/` takes a required thread count | Agent's full pod run on its branch: 830 passed, 10 skipped. `l1/domains.py` `scan` still defaults `cpus=1` (bounded, so left) |
| A-03 | this merge | README for `scripts/stress`; every script prints usage | 11 of 17 scripts had no usage path; README lists which scripts still use the old data model (folded into card B-18) |
| C-07 | this merge | `scripts/pod_run.sh` opens the memory-capped unit on the host; `run-step` and the search runner inherit a cap they find, or refuse when the profile requires one | Probed on the pod: a 1 GiB allocation in a 512 MiB unit was killed by the kernel (`oom_kill 1`); the pod has 1117 GB of swap, so `MemorySwapMax=0` is set everywhere. The search CLI cannot yet require a cap: card W-10 |
| W-10 | this merge | `bio-box search run --profile` (or `$BIOBOX_PROFILE`): thread cap, nice and required memory cap from the host profile; flags may only lower them | Agent's full pod run on its branch: 866 passed, 10 skipped. With a profile, `mem_gb` defaults to the profile's `mem_scope` |
| B-17 | this merge | Legacy keys removed: `slot_protein` table; ordinal, flank-slot and dimension columns; foreign keys and lineage now run through `feature` and `feature_protein`; `LEGACY_ORDINAL_KEYS` is empty | Schema lock regenerated on the pod with `--breaking`. Nested occurrence objects in search-hit exports changed (`ordinal`, no `genome_uid` or `tier`) without a schema-version bump: folded into card W-05 |
| W-05 | this merge | `search hits --format tsv\|faa`, `nbr <run> --format tsv`, `export proteins`; search-hit JSON is `search_hit/2` with `genome_uid` and `tier` joined from `contig_dim` | TSV checked against a brute-force join written in the test. Agent's full pod run on its branch: 866 passed, 10 skipped. Two lookups cannot prune (card B-27); `export proteins` globs when there is no manifest and FASTA keeps `*` (card W-11) |
| B-26 | this merge | Postings measured by column: 35.7 bytes per row, of which `feature_uid` is 17.0; estimate in the code corrected from 20 to 36 | The audit migrated one test constant that mirrors the estimate. Whether to slim the table is a PI decision: card 0-04 |
| B-25 | this merge | A release is defined once by `release_def` (ordered stages); clustering, tokens, postings, `contig_dim` and the manifest read only those stages; new `bio-box release define\|cluster\|tokens\|postings\|dim\|manifest` | The agent first reproduced the pod bug on the old code (all six release files changed when a later stage appeared). Query-time readers (`pack`, `nbr`, worker, `SeqStore`) still read every stage: card B-28 |
| B-23 | this merge | Ingest takes `out=` and can be rerun after a kill under `run-step` (SIGKILL at nine points, byte-identical result); refuses a stage directory that is not clean before reading anything | The card's premise was wrong: on the pod the second hash is about 3 s of 240 s and compression is 182 to 188 s. The audit had the agent restore the pass-2 digest check and test that a same-length edit between passes is refused. Agent's full pod run on the merged state: 901 passed, 10 skipped |
| B-22 | this merge | `l0_manifest` is partitioned by stage; `store_params.json` is a registered artifact, written once; a first and a second stage both publish under `run-step` and equal two direct ingests file for file | 32 other fixture files byte-identical before and after. Agent's full pod run: 903 passed, 10 skipped; golden 32 passed, 5 skipped. `contracts.check_append_only` still skips all small tables: card D-06 |
| B-18 | this merge | Phase-B exit measurements (`planning/measurements/2026-11-10_phase_b_exit.md`); stress scripts moved to the new model; CAMI II driver added; `.github/workflows/cross_host.yml` | Two of three acceptance criteria are **not met**: some catalog macros read every row group (card B-31), and two S5 storage budgets are exceeded (PI card 0-05). Cross-host equality holds on its first run: GitHub Actions and the pod both give manifest hash `ccd621b1...` for commit `f37c366` (fixture store; no clustering or domain step in it). Short-contig results: store open 0.03 s (was 7.2 s), tape lookup 1 of 15 row groups (was 181 of 181) |
| B-29 | this merge | Postings builder fixed: tokens are exploded one row group at a time into flat scratch runs, then sorted, so the sort spills; memory limit unchanged | Cause found by measurement: the `ORDER BY` over a union fed from a materialized CTE did not spill (1.9 GB) when family postings are sparse. Reproduced small (690,000 slots at a scaled limit). On real `rL.3`: old builder fails at 38 s; new one writes 25,822,709 rows in 105 s at 2.4 GB peak RSS, byte-identical to the old builder run with an 8 GB cap. 771 fixture files unchanged |
| W-11 | this merge | `export proteins` refuses without the release manifest; a check of 2,053,315 stored proteins in `rL.2` found no `*` anywhere, so FASTA output needed no change (a fixture test pins it) | Agent's full pod run: 906 passed, 10 skipped. Hit and locus exports still glob without a manifest: card W-12 |
| D-06 | this merge | The append-only check compares every file and artifact, small tables included | Fault injection on `l0_manifest`, `ledger`, `release_def`, `tape_version_def` and `store_params.json`: each rewrite is reported; clean on the two-stage, two-release fixture. Agent's full pod run: 906 passed, 10 skipped |
| B-28 | this merge | `pack build`, `nbr`, the search worker and `SeqStore` read only the release's stages (no default that globs); the clustering ledger is per release | The agent showed the bug on the old code: proteins planted in a stage outside the release were in its pack and its search hits. Other fixture files byte-identical. Stress scripts still call `SeqStore(root)`: card A-06. One release's manifest takes in other releases' files: card B-33 |
| A-06 | this merge | Ten stale `SeqStore(root)` calls in eight stress scripts fixed; a static R0 test parses the scripts and `src` signatures and fails on a stale call | `synth_scale.py` at 10× and `builder_limits.py` ran to completion on the pod; the pilot scripts were checked statically only. Agent's full pod run: 917 passed, 10 skipped |

Decisions DEC-18 to DEC-23 were approved by the PI on 2026-11-09 and are in `planning/DECISIONS.md`.

## In progress

| Card | Model | State |
|---|---|---|
| B-32 gene calling of short contigs uses its workers | Opus | started 2026-11-10 |
| C-08 every builder writes to an output root; manifest and catalog registered | Opus | started 2026-11-10 |
| B-33 a release's manifest covers only that release's files | Opus | started 2026-11-10 || C-04 `bio-box build-stage` | Opus | stopped 2026-11-10 with evidence, nothing written: builders cannot write into `run-step` staging. Now depends on B-22 and the new card C-08; card amended |
| Pod import ladder (not a card): phages, then 1 Gbp and 10 Gbp of complete prokaryotes with the new code | Opus | detached jobs on the pod; agent waiting on them |

## Next, in order

1. The path to a one-command stage build, in order: B-23 (in rework) → B-22 (`l0_manifest` per stage; register `store_params.json`) → C-08 (every builder writes to `out=`; after B-29) → C-04 (`build-stage`) → C-06 (runbook). Until C-04 merges, stages are built by a hand-written driver, as the import ladder was.
2. B-28 (query-time readers take the release's stages), B-24 (streaming clustering and sharded scans; memory, not time, is the limit), B-30 (ingest by batch), B-27 and W-11.
3. Pull request `plan/task-cards` → `main`, merged when the suite and CI are green (the PI delegated the timing). The merge builds the image and runs the image-only tests for the first time.
4. C-04 `build-stage`; C-06 runbook.
5. W-03 pack and runs outside the data root; W-05 readable exports; W-07 the one command.
6. B-21 per-domain postings.
7. Phase E (close S0), phase F (statistics; F-01 is a design card for a stronger model or the PI), phase G (S1 as a release).

## Not yet executed anywhere

These exist and pass their laptop-runnable parts, but their real run needs the image:
the clustering tests; the samtools slice and index tests; the HMMER-binary oracles (two tests);
the golden pairs test over postings; the Prodigal test (also needs `prodigal=2.6.3` added to
`container/environment.yml`, a PI decision); the mutation workflow; the systemd path of `run-step`
and of the search worker (only their command lines are unit-tested; the scope behaviour itself was
probed by hand on the pod).

## Open items for the PI

- Done 2026-11-09: the PI approved `prodigal=2.6.3` in `container/environment.yml` (test oracle only). The Prodigal test first runs in the image built from `main`.
- A `read:packages` credential, so the image built from `main` can be pulled on the pod.
- Confirm the day window for the pod's thread caps (profile assumes 08:00 to 20:00 on weekdays).
- First-10 questions (card 0-03), before the S1 gate.
- Card 0-05: two S5 storage budgets are exceeded (annotations 107 GB against 50; tapes plus postings 75 GB against 25 to 60), because every table now carries a 16-byte `feature_uid`. Recommendation: raise the budgets.
- Card 0-04: keep `token_postings` self-contained (about 36 GB at GTDB scale) or drop the columns the tape already holds (about a third smaller, at the price of a join). Recommendation: keep.
- `l0_manifest` is still one file rewritten per stage, against the append-only rule (found in B-15; no card yet: see B-22 in `TASKS.md`).

## On the pod (wilkcomp02)

- `~/bio-box`: `main`, unchanged. `~/bio-box-next`: `plan/task-cards` at `f9af5b6`, used by the import ladder.
- Provisional phage database with Pfam domains: `/stor/scratch/Wilke/if2453/all-the-data/provisional/rS1.pilot1` (old-code format; verified).
- Old-code pilots kept on local disk: `/tmp/if2453/pilot_s1` (phages), `/tmp/if2453/pilot_cami` (1.48 M short contigs).
- Import ladder with the new code: work in `/tmp/if2453/ladder`, results in `/stor/scratch/Wilke/if2453/all-the-data/provisional/rL`. Stage A (phages, `rL.1`), stage B (1 Gbp of complete prokaryotes, `rL.2`, finished 2026-11-09 23:12 CST) are done and clean. `rL.2p` (meant to be B plus Pfam) is **contaminated**: it was built while stage C's files were being written and took them in (its 12,048,866 features include stage C's genes, with clusters computed without them); use it only to try queries, not for counts. Stage C (10 Gbp, `rL.3`) failed at 02:53 CST on 2026-11-10 in the postings step (card B-29), and was finished at 04:10 CST with the fixed builder's postings: **`rL.3` is complete and verified on `/stor`** (32,815 contigs, 12,048,866 slots, 10,620,684 proteins, 25,822,709 postings, 13.5 GB; `planning/measurements/2026-11-10_import_ladder/rL3.md`). It has no Pfam for its newest stage (about 338 CPU-hours, not started), and it was built by the pre-B-17 code, so it is in the old layout. On `/stor`, `rL.1` and `rL.2` each fail verification on one file, `l0/l0_manifest.parquet`, which the old ingest rewrote per stage (fixed in the code by B-22). Results, scripts and findings: `planning/measurements/2026-11-10_import_ladder/`. New cards from it: B-25 (a release names its stages), B-26 (postings bytes); B-24 and C-04 amended.
- Workflow walk-through outputs: `/tmp/if2453/ssrnap`.
