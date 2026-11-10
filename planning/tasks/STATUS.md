# Task status

The running record of the cards in [`TASKS.md`](TASKS.md). Updated at every merge. Branch
`plan/task-cards` holds all merged work; nothing is on `main` yet. Each card was done by an agent on
its own branch `task/<ID>`, audited by the orchestrating session, and merged after the full suite
passed on the integration branch.

Last update: 2026-11-10, after midnight. Suite on the integration branch: 829 passed, 10 skipped
(*measured*, wilkcomp02, 2026-11-10, `sh scripts/pod_pytest.sh integration --shards 8`, in the
pinned image; 3 min 40 s wall for the slowest shard, against 13.5 min in one process). The skips
need samtools, bgzip, Prodigal or the golden set. Tests no longer run on the laptop (AGENTS.md §5.5).

## Merged (33 cards)

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

Decisions DEC-18 to DEC-23 were approved by the PI on 2026-11-09 and are in `planning/DECISIONS.md`.

## In progress

| Card | Model | State |
|---|---|---|
| B-17 remove the legacy ordinal columns | Opus | started 2026-11-10 |
| A-03 README for `scripts/stress` | Sonnet | started 2026-11-10 |
| C-07 scope on the host around the container | Opus | started 2026-11-10 |
| W-09 pyhmmer thread count from the caller | Sonnet | started 2026-11-10 |
| Pod import ladder (not a card): phages, then 1 Gbp and 10 Gbp of complete prokaryotes with the new code | Opus | detached jobs on the pod; agent waiting on them |

## Next, in order

1. B-18 phase-B exit measurements (after B-17).
2. B-22 `l0_manifest` per stage and W-05 readable exports (held until B-17 merges: they share `schemas.py`, the lock and `export.py` with it).
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

- Add `prodigal=2.6.3` to the test image, or accept that the Prodigal test stays skipped.
- A `read:packages` credential, so the image built from `main` can be pulled on the pod.
- Confirm the day window for the pod's thread caps (profile assumes 08:00 to 20:00 on weekdays).
- First-10 questions (card 0-03), before the S1 gate.
- `l0_manifest` is still one file rewritten per stage, against the append-only rule (found in B-15; no card yet: see B-22 in `TASKS.md`).

## On the pod (wilkcomp02)

- `~/bio-box`: `main`, unchanged. `~/bio-box-next`: `plan/task-cards` at `f9af5b6`, used by the import ladder.
- Provisional phage database with Pfam domains: `/stor/scratch/Wilke/if2453/all-the-data/provisional/rS1.pilot1` (old-code format; verified).
- Old-code pilots kept on local disk: `/tmp/if2453/pilot_s1` (phages), `/tmp/if2453/pilot_cami` (1.48 M short contigs).
- Import ladder with the new code: `/tmp/if2453/ladder` (in progress).
- Workflow walk-through outputs: `/tmp/if2453/ssrnap`.
