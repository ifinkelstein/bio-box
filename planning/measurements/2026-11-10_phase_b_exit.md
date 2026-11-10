# Phase-B exit measurement, 2026-11-10 (card B-18)

Purpose: rerun the 2026-11-09 design stress test
([`2026-11-09_design_stress_test.md`](2026-11-09_design_stress_test.md), "the 2026-11-09 file") on the
phase-B data model: feature keys (DEC-18), `token_postings` instead of `slot_pairs` (DEC-19), dimension
tables (DEC-20), batch shards (DEC-21), slot-bounded tape row groups (DEC-22). The per-table bytes and
rates of the import ladder ([`2026-11-10_import_ladder/README.md`](2026-11-10_import_ladder/README.md),
"the ladder") are reused, not re-measured, except where card B-17 changed a table after the ladder ran
(§6 corrects those from column-chunk sizes).

Host: wilkcomp02, x86-64, test image `bio-box-test_74b3e05.sif` (Python 3.14.10, pyarrow 26.0.0,
duckdb 2.0.0, pyrodigal 3.7.1, pyrodigal-gv 0.3.2, blake3 1.0.11, numpy 2.5.4; MMseqs2 and HMMER in
the image; no samtools, bgzip or Prodigal) with the checkout's `src/` first on `PYTHONPATH`, as
`scripts/pod_pytest.sh` sets it. This card changes no `src/`. Code, two rounds:

- **first round**, `src/` at `7a2a098`, tree `/tmp/if2453/ci/B-18`, launcher `run.sh`;
- **final round**, after `git merge plan/task-cards` at `c285538` (card B-25: release-level builders
  and `build_manifest` take the release's stage list, `release_def`), tree `/tmp/if2453/ci/B-18m`
  (`sh scripts/pod_pytest.sh B-18m tests/r0 -k nothing_matches_this`), launcher `run2.sh`. The
  golden build, the 1,000x run, the feature-stability run, the fixture hash, the projection, the
  golden ring and the suite were rerun here; the numbers below are from this round unless a row
  says "first round". Query counts at 1,000x and every feature-stability count were identical in
  both rounds; so were all 27 golden files the two rounds share (§6).

Every run below was a detached unit,

    systemd-run --user --unit=b18-<name> --collect -p MemoryMax=<M> -p MemorySwapMax=0 -p CPUQuota=<Q> \
      sh /tmp/if2453/b18/run2.sh <name> <script> <args>

where `run2.sh` (`run.sh` in the first round, with the tree `/tmp/if2453/ci/B-18`) checks the `STOP`
file and runs `nice -n 19 singularity exec --cleanenv --env
PYTHONPATH=/tmp/if2453/ci/B-18m/src:/tmp/if2453/pilot_s1/pydeps --env OMP_NUM_THREADS=2 --bind
/tmp/if2453 --bind /stor/scratch/Wilke/if2453/all-the-data/l0/refdb:/refdb:ro --bind
/sys/fs/cgroup:/sys/fs/cgroup:ro --pwd /tmp/if2453/ci/B-18m <image> python <script> <args>` with stdout
and stderr in `/tmp/if2453/b18/logs/<name>.out|.err` and a `.done` or `.failed` marker. Work in
`/tmp/if2453/b18`; nothing was written under `/stor`. Each unit had at most 16 threads
(`CPUQuota=1600%`), but the quotas of overlapping units summed to more at times (first round: the
golden build at 16 and the 1,000x run at 16, which mostly used one core; final round: CAMI at 8, the
golden build at 8 and one-core jobs); the
1,000x and CAMI runs had `MemoryMax=64G`, the others 4 to 32G. Times are 2026-11-10, 02:29 to 03:35
CST, while `ladder-C` (until 02:53) and other agents' jobs ran.

Labels: *measured*, *derived*, *estimate* as in `AGENTS.md` §6. Every number below is *measured* on
wilkcomp02 on 2026-11-10 by the command of its section unless it says otherwise.

## Acceptance at a glance

| Card criterion | Result |
|---|---|
| CI-runner and pod manifest hashes equal | **pending first CI run**: pod `b865e9e4d5beb6124c24d4c2f4995e67` (twice, final round); the CI half needs a push, which agents may not do (§5) |
| no query scans more than 2 row groups per lookup | **not met**: literal keys and keys from a joined table read 1 row group per lookup (sequence, tapes, postings, Q1, Q4), but a key passed through `token_of(...)` reads every postings row group, and `partners_of` reads every postings row group on its near side (§2; Open questions) |
| no projected table exceeds its budget | **not met**: S5 annotations 107 GB against ~50 GB, tapes plus postings 75 GB against 25–60 GB (§6, Misses); budget lines updated with the measured values |

## 1. Stress scripts on the phase-B model

| Script | Change |
|---|---|
| `build_s0.py` | writes `genome_input`; `--cluster` runs clustering, tokens, `postings.run`, `dim.run` (was `pairs.run`); `--workers` |
| `pilot_s1/pilot.py`, `pilot_resume.py` | postings and `contig_dim` steps replace `pairs_count` / `pairs_sample`; part files; `--workers`; `pilot_resume.py` is now `pilot.py --resume` |
| `pilot_s1/resume.sh` | runs in a `systemd-run --scope` with `MemoryMax` and `MemorySwapMax=0`; paths from the environment |
| `step6_export.py`, `step6_compare.py` | genes from `feature_protein`, ordered by coordinates and named by `feature_uid`; every shard's indexes |
| `step4_keys.py` | adds one ncRNA feature with the repo's builders and compares every table by `feature_uid` and by ordinal (§3) |
| `synth_scale.py` | features, synthetic clusters and families, then the real tape, token, postings and `contig_dim` builders; key lookups scored per key (§2) |
| `gtdb_pairs.py` | deleted; `gtdb_postings.py` projects every table (§6) |
| new | `build_fixture.py` (§5), `column_bytes.py` (§6), `gtdb_units_2026-11-10.json`, `pilot_cami/` (the CAMI II driver from the pod, card amendment) |
| `builder_limits.py` | also counts the point lookups `bb_seq` makes at 10^6 contigs (§7) |
| after the B-25 merge | `build_s0.py`, `build_fixture.py`, `synth_scale.py`, `macro_scan.py`, both pilots pass the release's stages to `cluster.run`, `tokens.run`, `postings.run`, `dim.run` and `build_manifest`; `step4_keys.py` reads them from the release's `release_def`; `store_probe.py` takes `--release` or `--stages` |

Every script's `--help` path exits 0 in the image (22 of 22, `sh /tmp/if2453/b18/helpcheck.sh`, a loop
of `python <script> --help` in the image; first round). Usage runs (first round): `pilot.py` on the first 30 INPHARED records
(`--limit 30 --threads 4 --workers 4`, all 11 steps without error in 28 s) and `pilot_resume.py`
on the same work directory (every step skipped); `step6_export.py` on the golden store (33,716 genes,
33,713 slices, 3 origin-wrapping genes skipped, as on 2026-11-09). Final round: `pilot.py` on 30
records again (`b18-mpilot30`, all 11 steps without error with the stage lists) and `store_probe.py
--release rS0.b18` on the final golden store (stages read from `release_def`: `['S0']`).

## 2. Layout at 1,000x the fixture

`synth_scale.py /tmp/if2453/ci/B-18m /tmp/if2453/b18/m/synth 1000 --shard-bases 67108864` (unit
`b18-msynth`, 64G, 4 CPUs; first round: unit `b18-synth1000`, 16 CPUs, the same counts in every row
below): 3,000 contigs, 797.5 Mbp, 854,000 slots, 20 shards. Sequences go
through `ingest`; gene features are the fixture's calls replicated per copy; clusters (NR90 per 10
copies, NR50 per 100) and families (664 families, 20% of proteins without one, 3% in one hub) are
synthetic, so bytes per row here are not representative (§6 uses real stores). Everything else is the
repo's builders: `tapes.run` (854,000 slots, 7 row groups), `tokens.run` (85,400 NR90, 8,540 NR50, 664
family tokens), `postings.run` (2,397,000 rows, 1 file, bucket width 0, 20 row groups), `dim.run`.
Build and queries 277 s in all (ingest 194 s, 4.1 Mbp/s; first round 309 s and 213 s); peak RSS of
the process 1,037 MB after the builders.

The `postings.run` out-of-memory failure that stopped ladder stage C (DuckDB `OutOfMemoryException` at
about 21 M slots with prokaryote-sized contigs; ladder README, "Stage C outcome"; card B-29) did
**not** reproduce here: this store has 854,000 slots (25 times fewer) in contigs of 48 to 580 kb, and
`postings.run` took 5.1 s. The 1,000x fixture is too small to test B-29; `synth_scale.py` was left
without a `--slots-per-contig` option for B-29 to add.

Row groups read, from `catalog.profile_query`, of the scan's total; for sequence, the row groups the
`bb_seq` point lookups read (`profile_query` `index`: `seq_locator` + `contig_index`, each 1 row group
in all at 3,000 contigs). "Per key" = row groups read / keys.

| Lookup | Keys | Read | Per key |
|---|---|---|---|
| sequence by contig, literal (`seq_get(<uid>, 0, 1000, 1)`) | 1 | index 1 + 1; 1 bgzf block | 2 |
| sequence, keys from a registered table (`FROM pick<k>`) | 1 / 10 / 1,000 | index 1 + 1 each; 1 / 10 / 996 blocks | 2 / 0.2 / 0.002 |
| sequence, 1,000 keys from `contig_dim_now()` joined to the table | 1,000 | `contig_dim` 1 of 1; index 1 + 1 | 0.002 |
| tape of a contig, literal (`tapes_of([<uid>])`) | 1 | `contig_tapes` 1 of 7, `tape_tokens` 1 of 7 | 1 |
| tape, keys from a table (`tapes_of((SELECT list(contig_uid) FROM pick<k>))`) | 1 / 10 / 1,000 | 1 / 6 / 7 of 7 each | 1 / 0.6 / 0.007 |
| tape, 1 key, plain join to `contig_tapes` | 1 | 1 of 7 | 1 |
| postings of a token, literal (`token_rows(?)`), family token with 1,000 loci; NR90 token with 10 | 1 | `token_postings` 1 of 20 | 1 |
| postings, key from `token_dict` through `token_of('nr90', ?)` | 1 | **20 of 20** (+ `token_dict` 1 of 1, three scans) | **20** |
| postings, key from `token_dict`, plain join on `token_id` | 1 | 1 of 20 | 1 |
| postings, 100 NR90 keys from a table, plain join | 100 | 8 of 20 | 0.08 |
| neighbours Q1 `q1_count(a, b, 3, NULL)`, literal family pair; literal NR90 pair | 2 | 1 + 1 of 20 | 1 per side |
| neighbours Q1 NR90 pair through `token_of` | 2 | **20 + 20 of 20** | **20 per side** |
| neighbours Q4 `pairs_of(a, b, 10)`, literal family pair | 2 | 1 + 1 of 20 | 1 per side |
| neighbours Q4 NR90 pair through `token_of` | 2 | **20 + 20 of 20** | **20 per side** |
| neighbours Q3 `partners_of(a, 10)`, literal family | 1 | anchor 1 + 1 of 20; near side **20 of 20** | **20** |
| neighbours Q3 one NR90 through `token_of` | 1 | **20 of 20** on all three scans | **20** |
| hub family (26,000 loci): Q1 hub x family; Q3 partners of the hub (36,222 partner tokens) | | 1 + 1; near side 20 of 20 | not a key lookup |

Bulk reads, against 2026-11-09 §5 (12 shards then, 20 now):

| Read | 2026-11-09 (laptop, pre-phase-B) | 2026-11-10 |
|---|---|---|
| 1,000 neighbourhoods of +-10 genes, sequence | 1,146 blocks, 22.8 MB; `contig_index` scanned 3 times, `seq_locator` twice | 17.0 Mbp: 1,153 blocks, 23.0 MB; index row groups 1 + 1, no index scan (B-13, B-14) |
| tapes and tokens of those 1,000 neighbourhoods | the only row group, all 3,000 contigs | 7 of 7 row groups (at most 131,072 slots each, DEC-22) |
| export of 500 whole contigs | 135 Mbp, 2,463 blocks, 48.7 MB (0.36 B/bp), 12 of 12 shards | 131.1 Mbp, 2,216 blocks, 43.8 MB (0.334 B/bp, *derived*), 20 of 20 shards |

So the layout prunes: every lookup whose key is a literal or comes from a joined table reads at most 2
row groups per key. Two macro paths do not, and both are defects in `src/biobox/catalog.py` (Open
questions): a key given as `token_of(...)`, a scalar subquery, is not pushed into the postings scan
(the same key through a plain join is); and `partners_of` finds the near-side postings by `contig_uid
IN (...)` on a table sorted by `token_id`, so it reads every row group, while TASKS.md §2 specifies
that it read `tape_tokens` for those contigs (sorted by `contig_uid`, which prunes).

## 3. Feature stability (the 2026-11-09 §4 measurement)

`step4_keys.py /tmp/if2453/b18/m/golden/root /tmp/if2453/b18/m/step4 rS0.b18` (unit `b18-mstep4`) on
the golden store of §6 (built with Pfam and clusters); the first round (`b18-step4b`) gave the same
counts in every row. It copies the store, adds one ncRNA (kind 2, 377 nt,
minus strand, 3,318,066 to 3,318,443 on the *E. coli* K-12 chromosome, 60 nt into the first IGR of at
least 500 nt at or after 3.27 Mb: the same place as on 2026-11-09) as a second producer `b18-ncrna`,
and runs `igr.run`, `tapes.run`, `tokens.run` and `postings.run` for the tape version of both producers
(`t-cce4e4c57a9f`, against `t-ca596dd85040`). The new feature is slot 3,104 of 4,319 on that contig.
"Rewritten" = key gone or any other column changed.

| Table (key) | Rows before | Rewritten | New | 2026-11-09, ordinal key |
|---|---|---|---|---|
| `feature` (`feature_uid`) | 33,716 | 0 | 1 | |
| `feature_protein`, `protein_occurrence` (`feature_uid`) | 33,716 each | 0 | 0 | 1,214 slot rows of 4,318 |
| `igr` (coordinates; flanks by `feature_uid`) | 28,161 | 1 | 2 | 1,006 of 3,562 |
| `token_dict` (`token_id`) | 68,799 | 0 | 0 | |
| tape slot entries keyed by `feature_uid`: coordinates, kind, protein | 33,716 | 0 | 1 | |
| tape slot entries keyed by `feature_uid`: all list columns | 33,716 | 1 (the slot after the ncRNA: its junction and gap to the previous slot) | 1 | |
| `tape_tokens` slot entries keyed by `feature_uid` | 33,716 | 0 | 1 | |
| `token_postings` (`token_id, contig_uid, feature_uid, seg_idx`) without `ordinal` and `n_slots` | 94,827 | 0 | 0 | 24,342 of 86,360 pair rows |
| `token_postings`, primary key, all columns | 94,827 | 12,755 (every posting on that contig: `n_slots`; 3,585 of them also `ordinal`) | 0 | |
| the same slot entries keyed by ordinal | 33,716 | 1,214 | 1 | |
| `token_postings` keyed by `(token_id, contig_uid, ordinal, seg_idx)`, without `n_slots` | 94,827 | 3,585 | 3,542 | |

Under the feature key nothing that names a feature changes. The 1,214 slots after the insertion keep
their `feature_uid` and move one ordinal (as on 2026-11-09, where the same 1,214 rows were rewritten
under the ordinal key). Every posting row carries its contig's `n_slots`, so a feature added to a
contig changes that column in all of the contig's postings of the new tape version (12,755 rows here);
they are written as new files of the new tape version, never rewritten in place (AGENTS §3.10). The
contig rows of `contig_tapes` and `tape_tokens` (one row per contig, lists inside) change for that one
contig.

## 4. Independent-tool checks

`sh scripts/pod_pytest.sh B-18g --golden tests/r3 -rs` (golden set at `/tmp/if2453/golden`), final
round: **32 passed, 5 skipped in 227 s** (first round: the same, 221 s). The skips:

| Test (`tests/r3/test_golden_independent.py`) | Reason |
|---|---|
| line 143, 171: samtools slices and `.fai` / `.gzi` | samtools is not in the `74b3e05` image on the pod |
| line 187: `bgzip -t` | bgzip is not in that image |
| line 197: Prodigal 2.6.3 | Prodigal is not in that image (added to `container/environment.yml` on 2026-11-09; the image built from `main` has it) |
| line 239: RefSeq CDS of *E. coli* K-12 | the truth file `truth/GCF_000005845.2_ASM584v2_genomic.gff.gz` is not under the pod's golden directory |

The four tool skips wait for the image built from `main` (card 0-02: the registry credential). For
the fifth, a golden view with the truth file (`/tmp/if2453/b18/golden_view`: symlinks to every golden
genome directory plus the laptop's copy of the GFF, sha256 `afdf03dc…f495e`):
`BIOBOX_POD_GOLDEN=/tmp/if2453/b18/golden_view sh scripts/pod_pytest.sh B-18g --golden
tests/r3/test_golden_independent.py -rs` gave **11 passed, 4 skipped in 136 s** (first round): the
RefSeq test passes on the phase-B store.

## 5. Cross-host equality

The laptop no longer computes (AGENTS §5.5), so the second host is a GitHub Actions runner.
`.github/workflows/cross_host.yml` (manual `workflow_dispatch`, and pushes to `plan/task-cards` and
`main`) sets up the environment exactly as `ci.yml` does (`actions/checkout@v7`,
`astral-sh/setup-uv@v10.3.3`, `uv sync --python 3.14`; no new action or dependency) and runs
`uv run python scripts/stress/build_fixture.py . "$RUNNER_TEMP/fixture_root" rX.fixture`, then prints
the report (`show_report.py`, whose first line has the manifest hash, and the JSON).

Pod, in the image:

    systemd-run --user --unit=b18-mfixture --collect -p MemoryMax=8G -p MemorySwapMax=0 -p CPUQuota=200% \
      sh /tmp/if2453/b18/run2.sh m_fixture scripts/stress/build_fixture.py /tmp/if2453/ci/B-18m \
      /tmp/if2453/b18/m/fixture/root rX.fixture

| | Pod run 1 | Pod run 2 (`m_fixture2`, another root) | CI runner |
|---|---|---|---|
| manifest hash, final round (code `c285538` merged) | **`b865e9e4d5beb6124c24d4c2f4995e67`** | `b865e9e4d5beb6124c24d4c2f4995e67` | **pending first run** |
| manifest hash, first round (code `7a2a098`) | `36664c57ecbc9aac01c61c3de3909a0e` | `36664c57ecbc9aac01c61c3de3909a0e`, 19 of 19 files byte-identical | |
| build | 2.6 s: 3 genomes, 797,481 bp, 854 genes, 525 IGRs, 3 tapes, 3 `contig_dim` rows | | |

The hash changed between the rounds because card B-25 added the `release_def` file and the release's
stages to the manifest (`MANIFEST_VERSION` 3); compare a CI run with the hash of the same code.

The CI half is **pending its first run**: an agent may not push (AGENTS §1.3, §10), so the workflow has
never run. When it has, compare its hash with the pod's above (or save its JSON and run
`diff_reports.py`). The report names the package versions on both hosts; the pod's are in the header
of this file. A mismatch in pyarrow's version alone changes every Parquet footer (`created_by`) and
the ledger's tool versions.

## 6. Projection to GTDB R226 scale

Unit counts: `gtdb_counts.py /tmp/if2453/b18/gtdb` (GTDB R226 `bac120` and `ar53` metadata, the copy
made for the ladder) reproduced the 2026-11-09 counts exactly (S2C 25,458 genomes, 106.4 Gbp, 56,387
contigs, 97.65 M proteins; S3 38,750, 52.4 Gbp, 4.17 M, 53.2 M; S5 143,614, 465.5 Gbp, 24.95 M,
437.2 M; S6 732,475, 2,587.6 Gbp, 111.5 M, 2,436 M).

Bytes per unit, current schema. Sources: **ladder P1** = stage P1 of the import ladder (250 complete
GTDB genomes, 1.02 Gbp, 539 contigs, 927,920 genes), read with `column_bytes.py
/tmp/if2453/ladder/root/manifests/release=rL.2.json --where stage=P1` (and `--where release=rL.2` for
release tables), code `f9af5b6`; **golden** = the B-18 golden build,
`build_s0.py /tmp/if2453/ci/B-18m /tmp/if2453/golden /tmp/if2453/b18/m/golden/root rS0.b18 --cluster 8
--pfam /refdb/pfam/38.2/Pfam-A.hmm.gz --cpus 8 --workers 8` (unit `b18-mgolden`, manifest
`c3db3865eb834e17fec837576d9bd137`: ingest 10.1 s, gene calling 11.7 s at 8 workers, Pfam 587 s at 8
threads, clustering 41.5 s at 8 threads, tokens 0.6 s, postings 0.3 s; peak RSS 1,288 MB). The first
round built the same store at 16 workers and 16 threads (unit `b18-golden`, manifest
`5e09e031cbbaf7005068f63dc9e2d39c`, Pfam 336 s, clustering 28 s): `diff_reports.py` finds all 27
files the two builds share byte-identical (gene calls, Pfam domains, MMseqs2 clusters and alignments,
tokens and postings included); the final build has one more file, its `release_def`. So the outputs do
not depend on the thread count between 8 and 16 (AGENTS §4.3). Card B-17 removed columns after the
ladder ran; column chunks compress independently, so the ladder's bytes minus those columns' chunks
are *derived* current-schema values.

| Table | Unit | Bytes per unit | Source | 2026-11-09 (old schema) |
|---|---|---|---|---|
| shards + `.gzi` | bp | 0.3056 | ladder P1 | 0.309 with `.fai` (golden) |
| `.fai` | contig | 88.0 | ladder P1 | |
| `protein_seq` | distinct protein | 228.5 | ladder P1 (golden 226.6) | 227 |
| `feature` | gene | 26.0 | ladder P1 (golden 26.7) | |
| `feature_protein` | gene | 43.0 | ladder P1 (golden 44.4) | 19.5 (`slot_protein`) |
| `protein_occurrence` | gene | 51.6 | *derived*: ladder P1 54.44 − `slot_ordinal` 1.82 − `genome_uid` 1.02 − four dimension columns 0.01 (golden 47.6) | 32.2 |
| `igr` | row; gene | 43.0; 35.5 | *derived*: ladder P1 46.45 − `left_slot` 1.71 − `right_slot` 1.72; 0.824 rows per gene (golden 44.6 per row) | 12.7 per row |
| `contig_tapes` | slot | 47.9 | *derived*: ladder P1 47.91 − four dimension columns 0.01 (golden 48.4) | 30.4 |
| `tape_tokens` | slot | 11.0 | golden, three levels (ladder rL.2, mostly two levels: 9.2) | 9.0, two levels |
| `token_postings` | row; gene | 35.74; 100.5 | bytes: the rL.2 store, 13,330,867 rows (card B-26, [`2026-11-10_postings_bytes.md`](2026-11-10_postings_bytes.md); the rL.2 release file alone, 7,593,431 rows, is 36.59 by `column_bytes.py`; golden 21.1, 21 contigs); rows: golden 2.81 per gene (NR90 and NR50 on every gene, a Pfam family on 81.3%; B-26's 2.329 per slot is the S1+P1 mix, mostly phages) | `slot_pairs` 7.8 per row, 422 per gene |
| `token_dict` | gene | 12.7 | ladder rL.2: 0.650 tokens per gene × 19.5 B | |
| `cluster_member`, `_aln` | protein and level | 36.2; 34.2 | ladder rL.2 (golden 36.5; 26.2) | 36.5; 26.2 |
| `l1_domains`, Pfam 38.2 | row; gene | 52.5; 89.1 | golden: 56,786 domains on 33,442 proteins | 53.8; 120 with NCBIfam |
| `contig_index`, `contig_occurrence`, `seq_locator`, `genome_input`, `contig_dim` | contig | 57.7, 95.8, 40.0, 28.7, 39.7 | ladder P1 (`contig_dim`: rL.2) | |
| ledger; `l0_manifest` | contig; genome | 77.6; 235 | ladder P1 (2.46 rows per contig × 31.5 B); golden | |

Projection, `gtdb_postings.py /tmp/if2453/b18/gtdb scripts/stress/gtdb_units_2026-11-10.json` (unit
`b18-mgtdbproj`), *derived*: units × bytes per unit; `protein_seq` and clusters assume no
deduplication (upper bound). GB = 10⁹ bytes.

| Table | S2C | S3 | S5 | S6 |
|---|---|---|---|---|
| shards + `.gzi` + `.fai` | 32.5 | 16.4 | 144.5 | 800.7 |
| `protein_seq` | 22.3 | 12.2 | 99.9 | 556.7 |
| `feature` | 2.5 | 1.4 | 11.4 | 63.3 |
| `feature_protein` | 4.2 | 2.3 | 18.8 | 104.8 |
| `protein_occurrence` | 5.0 | 2.7 | 22.6 | 125.7 |
| `igr` | 3.5 | 1.9 | 15.5 | 86.4 |
| `l1_domains` (Pfam only) | 8.7 | 4.7 | 39.0 | 217.1 |
| `cluster_member` | 7.1 | 3.9 | 31.6 | 176.3 |
| `cluster_member_aln` | 6.7 | 3.6 | 30.0 | 166.9 |
| `contig_tapes` | 4.7 | 2.5 | 20.9 | 116.7 |
| `tape_tokens` | 1.1 | 0.6 | 4.8 | 26.7 |
| `token_postings` | 9.8 | 5.3 | 44.0 | 244.9 |
| `token_dict` | 1.2 | 0.7 | 5.6 | 31.0 |
| contig tables (`contig_index`, `contig_occurrence`, `seq_locator`, `genome_input`, `contig_dim`), ledger, `l0_manifest` | 0.0 | 1.4 | 8.5 | 38.0 |
| **total** | **109.4** | **59.7** | **496.9** | **2,755.2** |

For comparison, the 2026-11-09 projection of `slot_pairs` at S5 was 128–178 GB at 7.7 B per row and
330–465 GB at 20; `token_postings` is 44.0 GB, and its 2026-11-09 estimate was 24 GB. B-26's sizing
line gives ~36 GB from the same bytes per row with 2.329 rows per slot (the phage-heavy S1+P1 mix);
the golden set's 2.81 rows per gene is the prokaryote mix, with a Pfam family on 81% of genes.

Against the budget lines of `planning/sources_and_sizing.md` (S5 disk per layer, §4), now updated with
these values and labelled:

| Budget line (S5) | Budget | Projected | Within? |
|---|---|---|---|
| sequence store (nt + aa, bgzip + indexes) | ~207 GB | 244.4 GB (+ 8.5 GB contig tables) | yes, +18% |
| per-genome annotations (features, IGRs, Pfam hits) | ~50 GB | 107.2 GB | **no** |
| cluster membership | ~32 GB | 31.6 GB (+ 30.0 GB alignments, no budget line) | yes |
| tapes + pairs (now postings) + `pair_stats` | 25–60 GB | 75.3 GB without `pair_stats` | **no** |
| DuckDB catalog | <1 GB | 1.06 MB for rL.2 (46 files) | yes |

Two columns will grow at scale beyond these measured values (*derived*, not in the table above): in
`protein_occurrence` (sorted by `protein_uid`) `contig_uid` took 1.2 B per row with 539 contigs, but
every other unclustered 16-byte UID column took 17.0; at S5's 24.95 M contigs it would approach that,
giving 67.4 B per gene and 29.5 GB at S5. In `token_postings`, `contig_uid` took 4.76 B per row at
27,797 contigs (B-26); at 17.03 B the row is 48.0 B and S5 59.0 GB. Card B-26 owns the postings bytes.

### Misses

Measured or derived values more than 20% from the estimate a plan document or the code gives. Card
B-26 handles `token_postings`; nothing in `src/` was changed here.

| Item | Estimate (source) | Measured / derived (source) | Miss |
|---|---|---|---|
| `token_postings` bytes per row | 20 (`postings.py` `EST_POSTING_ROW` before B-26; PLAN.md §4.2 size model) | 35.74 (rL.2 store, card B-26; 36.59 for the rL.2 release file alone) | +79% (card B-26 corrected the code estimate to 36) |
| `token_postings` at S5 | 24 GB (2026-11-09 §3, *estimate*) | 44.0 GB (*derived*, §6) | +83% |
| S5 tapes + postings (+ `pair_stats`) | 25–60 GB (`sources_and_sizing.md` §4) | 75.3 GB without `pair_stats` (*derived*) | +25% above the upper end |
| S5 per-genome annotations | ~50 GB (`sources_and_sizing.md` §4, from a superseded estimate) | 107.2 GB (*derived*) | +114% |
| tape bytes per slot | 30.4 (S5 tapes line; golden 2026-11-09) and 30 (`tapes.py` `EST_SLOT_BYTES`) | 47.9 (ladder P1, *derived* for the B-17 columns): `feature_uids` adds 17.1 B per slot | +58%, +60% |
| `.fai` per record | ≈ 50 B (`sources_and_sizing.md` §3) | 88.0 B (ladder P1) | +76% |
| `tape_tokens` bytes per slot | 30 (`tokens.py` `EST_TOKEN_SLOT`) | 11.0 (golden), 9.2 (ladder rL.2) | −63% |
| `token_dict` bytes per row | 60 (`tokens.py` `EST_DICT_ROW`) | 19.5 (ladder rL.2) | −67% |
| `contig_occurrence` bytes per row | 120 (`store.py` `EST_OCCURRENCE_ROW`) | 95.8 (ladder P1) | −20.1% |
| `genome_input` bytes per row | 36 (`inputs.py` `EST_ROW_BYTES`) | 28.7 (ladder P1) | −20.3% |

Within 20% of their estimates: shards 0.3056 against 0.309 B/bp; protein 0.713 against 0.73 B/aa;
`feature` 26.0 against 24; `feature_protein` 43.0 against 40; `protein_occurrence` 51.6 against 48;
IGR 43.0 against 45 B per row; `contig_index` 57.7 against 56; `seq_locator` 40.0 against 40;
`contig_dim` 39.7 against 44; `l1_domains` 52.5 against 54 B per row and 89 against 90 B per protein
(§4 annotations line); `cluster_member` 36.2 against 37 and `_aln` 34.2 against 31. The code estimates
only cut parts and pick bucket widths; none decides correctness.

## 7. Other reruns

First round, except the CAMI row (`builder_limits.py`, `macro_scan.py` and the `SeqStore` paths they
measure are untouched by card B-25; `store_probe.py` gave the same open time and blocks per fetch in
the final round).

| Script and command | Result | 2026-11-09 |
|---|---|---|
| `builder_limits.py /tmp/if2453/b18/golden/root /tmp/if2453/b18/limits` | pair oracle `pairs_of_contig` on the largest tape (4,903 slots, 98,060 rows, one level): 352 B of heap per row (nothing builds pairs now). `SeqStore` open on the synthetic 10^6-contig index: 0.0 s, 0 MB held. `bb_seq` point lookups: 1 key reads 1 of 9 `contig_index` and 1 of 9 `seq_locator` row groups; 1,000 random keys with no row-group cache read 1,000 + 1,000 (2.0 per key) in 24.7 s | 544 B per row (builder); open 2.0 s and 414 MB |
| `macro_scan.py /tmp/if2453/b18/limits/locator_1e6` | user macros `default_digest` / `digest_len`: literal key 1 of 9 row groups each; 1 joined key and 1,000 joined keys 9 of 9 each (unchanged) | the same |
| `pilot_s1/store_probe.py /tmp/if2453/b18/golden/root ...` | open 0.033 s, 2 footers, 0 row groups; 1,000 +-10-gene neighbourhoods 0.52 ms each, 1.17 blocks each (1.34 without the block cache) | |
| `pilot_cami/pilot_cami.py` and `probe_cami.py` | below | §6c |

**CAMI II marine sample 0, rerun** (the 2026-11-09 §6c pilot; `pilot_cami.py
/tmp/if2453/pilot_cami/src/x/simulation_short_read/2018.08.15_09.49.32_sample_0/contigs/anonymous_gsa.fasta.gz
/tmp/if2453/b18/cami rCAMI.b18 --threads 8 --workers 8`, then `probe_cami.py /tmp/if2453/b18/cami
rCAMI.b18 --threads 8`; 64G, 8 CPUs). Ingest to tapes ran in the first round (unit `b18-cami`, tree
`B-18`); the unit was stopped when clustering started and resumed in the final round (unit `b18-cami2`,
tree `B-18m`), which skipped the finished steps (the resume path of the driver). Every step ran without
error; manifest `32001f752ecd7fa6cf62268a8e8e6f0e`.

| Item | 2026-11-09 (pre-phase-B) | 2026-11-10 |
|---|---|---|
| contigs; bases | 1,475,972 records; 909.5 Mbp | the same; 1,475,957 distinct |
| ingest | 296 s, 3.9 GB | 327 s (2.8 Mbp/s), 2.4 GB |
| gene calling | 2,322 s on 1 thread, 7.8 GB | 1,989 s with 8 workers, 2.7 GB; 1,963,951 genes, 98.8% distinct proteins |
| IGRs; tapes | ; 824 s | 15 s; 54 s, 15 row groups of at most 100,672 contigs |
| clustering | | 702 s at 8 threads, 8.9 GB; NR90 89.4% singletons (89%) |
| tokens; postings | | 33 s; 14 s, 3,927,902 rows, **41.6 B per row** |
| tapes, bytes per slot; per contig | 52.5; 70 | 55.0; 73.2 |
| `.fai` per contig | 85 | 85.0 |
| `SeqStore` open | 7.2 s, 1.15 GB | **0.03 s, 7.3 MB**, 2 footers, 0 row groups |
| one-contig `tapes_of` | 181 of 181 row groups, 53 MB | **1 of 15** row groups (tapes and tape tokens), 3.6 MB |
| 1,000 contigs, `tapes_of` | | 15 of 15 (0.015 per key) |
| 1,000 gene fetches (gene +-1 kb), direct | 8.5 s, 1.01 blocks each | 1.19 s, 1.01 blocks each |
| 1,000 gene fetches, catalog macro | 15.2 s; `contig_index` scanned whole 3 times, `seq_locator` twice, 326 MB of Parquet | **0.97 s; no Parquet scan**, 974 blocks; 0 mismatches against the direct fetch |

With 1.48 M contigs the postings' `contig_uid` no longer dictionary-encodes: 41.6 B per row against
35.7 on rL.2, in line with the *derived* 48 B at S5 scale in §6. Postings did not run out of memory
here (1.96 M slots, short contigs; B-29's failure needs about 21 M slots in long contigs).

The user macros `default_digest`, `digest_len` and `contig_len_of` scan the index views by design
(`catalog.macro_sql` docstring, card B-14), and the sequence macros no longer call them; but AGENTS
§4.4 binds every SQL macro that resolves a key, so they are listed under Open questions.

## 8. Not rerun, and why

- The S1 pilot (2026-11-09 §6a): the ladder's stage A is the same input built with the phase-B code
  (release `rL.1`); its numbers are reused. The driver was run on 30 records to check its usage path.
- NCBIfam on the golden set: its HMM file is a `.tgz` that needs unpacking; Pfam alone gives the
  family coverage the postings projection needs.
- Mutation testing (2026-11-09 §6) belongs to cards D-01 and D-02.

## 9. Defects and open questions found here (no `src/` was changed)

1. **A key given through `token_of(...)` is not pruned.** `token_rows(token_of('nr90', ?))`, and
   `q1_count`, `pairs_of`, `partners_of` called with `token_of(...)` arguments, read 20 of 20
   `token_postings` row groups per side at 1,000x; the same key through a plain join from
   `token_dict` reads 1 of 20 (§2; `synth_scale.py`, queries "... through `token_of`"). DuckDB does not
   push a scalar subquery's value into the Parquet scan's zone-map filter, while it does push a join
   key. AGENTS §4.4 binds every macro that resolves a key "including when the keys come from a joined
   table". A fix would join the macros' key in (`... JOIN (SELECT _tok AS t) k ON p.token_id = k.t`)
   rather than compare to it; it belongs to a card on `src/biobox/catalog.py` with a
   `profile_query` test that passes `token_of(...)`.
2. **`partners_of` reads every postings row group.** Its near side selects `token_postings` by
   `contig_uid IN (anchor contigs)` on a table sorted by `token_id`, so nothing prunes (20 of 20 at
   1,000x for a family with 1,000 loci, and for one NR90 token with 10). TASKS.md §2 specifies that
   `partners_of` read `tape_tokens` for those contigs, which is sorted by `contig_uid`; at S5
   (about 1.2 x 10⁹ postings, *derived* from §6) every call would scan the whole table.
3. **The user lookup macros still scan the index.** `default_digest`, `digest_len` and
   `contig_len_of` read 9 of 9 row groups for one key from a joined table at 10^6 contigs (§7,
   `macro_scan.py`), as on 2026-11-09. The sequence macros no longer use them (card B-14), but they
   are macros that resolve a key (AGENTS §4.4), and `contig_len_of` is what an export query reaches
   for first.
4. **Budgets.** Two S5 budget lines are exceeded (§6): the per-genome annotations (107 GB against
   ~50) and tapes plus postings (75 GB against 25–60). Most of the growth is the 16-byte random UIDs
   that phase B added (`feature_uid` costs 17.0 B per row in each of `feature`, `feature_protein`,
   `protein_occurrence` and the tape lists, in each of the 2.81 postings per gene, and twice in each
   of the 0.824 IGRs per gene): 8.46 copies x 17.0 B = 144 B per gene, 63 GB at S5 (*derived*). Whether to raise the budget lines or to store fewer UID
   copies is a PI decision.
5. **Gene calling of short contigs does not use its workers.** CAMI II gene calling with
   `--workers 8` took 1,989 s for 1.48 M contigs, at 194% CPU in one `top` sample 13 minutes in, against
   2,322 s on one thread on 2026-11-09. The likely bound (not measured) is per-contig Python work in the
   calling thread (sequence reads, scratch rows), not the finder.
6. The cross-host CI half is pending its first run (§5), and the four tool-dependent golden-ring tests
   wait for the image built from `main` (§4).
