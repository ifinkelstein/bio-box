# Stress scripts

These scripts produced the measurements in
[`planning/measurements/2026-11-09_design_stress_test.md`](../../planning/measurements/2026-11-09_design_stress_test.md)
("the 2026-11-09 file") and, after card B-18 moved them to the phase-B data model, in
[`planning/measurements/2026-11-10_phase_b_exit.md`](../../planning/measurements/2026-11-10_phase_b_exit.md)
("the exit file"). Section numbers below are those files' sections. Card A-03 documented the scripts;
card B-18 updated them and reran them. The numbers live in the measurement files, not here.

Rules that apply to every run (AGENTS.md sections 4, 6, 7):

- No compute on the laptop. Run on the pod (wilkcomp02) in the test image, with the checkout's `src/`
  first on `PYTHONPATH`. `sh scripts/pod_pytest.sh <name> tests/r0 -k nothing_matches_this` mirrors the
  working tree to `/tmp/if2453/ci/<name>`; the exact command prefix is in `scripts/pod_pytest.sh`
  (`systemd-run --user ... -p MemoryMax=... -p MemorySwapMax=0 nice -n 19 singularity exec --cleanenv
  --env PYTHONPATH=<tree>/src:<pydeps> --bind /tmp/if2453 --pwd <tree> <image> python ...`). A
  script's `REPO` argument is that tree. Anything longer than a few minutes is launched detached
  (`systemd-run --user --unit=<name> --collect ...`); the exit file lists the commands B-18 used.
- Every script prints usage and exits with `--help` (or `-h`) before doing any work. With too few
  arguments the stdlib-style scripts print the same usage and exit 2.
- Scripts print JSON or text to stdout; redirect it to a file under `/tmp/if2453/`. Nothing here
  writes outside the output directory it is given, and none writes under its input store.

## Index

| Script | Host | Measures | Section | Status |
|---|---|---|---|---|
| `build_s0.py` | pod | build of the golden set, per-table size, manifest hash | 2026-11-09 §1, §3; exit §3, §6 | current (B-18) |
| `build_fixture.py` | pod and CI | stage-A build of the micro-S0 fixture and its manifest hash: the cross-host check | exit §5 | new (B-18) |
| `show_report.py` | any (stdlib) | compact view of a `build_s0.py` or `build_fixture.py` report | | current |
| `diff_reports.py` | any (stdlib) | two reports: manifest hash and per-file sha256 | exit §5 | current |
| `column_bytes.py` | pod | bytes per table and per column from the Parquet footers of a manifest | exit §6 | new (B-18) |
| `step6_export.py` | pod | exports gene slices from a built store for the independent check | 2026-11-09 §2 | current (B-18); its checks now also run as `tests/r3/test_golden_independent.py` |
| `hash_fa.py` | pod (stdlib) | sha256 per record of a `samtools faidx` FASTA | 2026-11-09 §2 | current |
| `prodigal_run.py` | pod (stdlib + `prodigal`) | Prodigal single mode on the ten golden prokaryotes | 2026-11-09 §2 | current |
| `step6_compare.py` | any | compares bio-box with samtools, Prodigal and RefSeq CDS | 2026-11-09 §2 | current (B-18) |
| `gtdb_counts.py` | pod (needs `duckdb`) | unit counts per stage from GTDB R226 metadata | 2026-11-09 §3 | current |
| `gtdb_postings.py` | pod (needs `duckdb`) | every table projected to GTDB R226 scale from measured bytes per unit | exit §6 | new (B-18); replaces `gtdb_pairs.py` |
| `step4_keys.py` | pod | one added ncRNA feature: rows rewritten in every table, by `feature_uid` and by ordinal | 2026-11-09 §4; exit §3 | current (B-18) |
| `synth_scale.py` | pod | 1,000x layout and row groups read per key lookup | 2026-11-09 §5; exit §2 | current (B-18) |
| `builder_limits.py` | pod | heap per pair row of the pair oracle and `SeqStore` open cost at 10^6 contigs | 2026-11-09 §5 | runs; the pair half measures the test oracle only (pairs are not built, AGENTS §3.5) |
| `macro_scan.py` | pod | row groups read by the user lookup macros at 10^6 contigs | 2026-11-09 §5 | runs; the sequence macros no longer call these (card B-14) |
| `postings_bytes.py` | pod | `token_postings` bytes per column from Parquet footers (`python scripts/stress/postings_bytes.py ROOT RELEASE [--json]`; read only, card B-26) | `planning/measurements/2026-11-10_postings_bytes.md` | current |
| `protein_stars.py` | pod | sequences of `protein_seq` that end in `*` or hold an internal `*`, files from the release manifest (`python scripts/stress/protein_stars.py ROOT RELEASE`; read only, one thread, card W-11) | `planning/measurements/2026-11-10_protein_stars.md` | new (W-11) |
| `build_stage_pod/` | pod | card C-10: `make_inputs.py` (the `--genomes`/`--objects` files of `build-stage` from the import ladder's inputs), `pfam_plan.py`, `pfam_shards.py` (host stdlib launcher, at most 8 shard processes), `pfam_merge.py` (the shards as `bio-box stage domains` writes them, run as a `run-step` command), `watch_then.sh` (creates `STOP` or stops a unit when a path appears); usage in each file's header | `planning/measurements/2026-11-10_build_stage_pod.md` | new (C-10) |
| `pilot_s1/pilot.py` | pod | S1 rate pilot on INPHARED phages | 2026-11-09 §6a, §6b | current (B-18): postings replace the pair steps |
| `pilot_s1/pilot_resume.py` | pod | `pilot.py --resume` | 2026-11-09 §6a | current (B-18) |
| `pilot_s1/resume.sh` | pod | waits for an earlier process, then runs `pilot_resume.py` in a memory scope | 2026-11-09 §6a | current (B-18): scope added, paths from the environment |
| `pilot_s1/store_probe.py` | pod | `SeqStore` open cost and 1,000 neighbourhood fetches in a fresh process | 2026-11-09 §6a, §6c | current |
| `pilot_cami/pilot_cami.py` | pod | short-contig pilot on CAMI II marine sample 0 | 2026-11-09 §6c | added by B-18 from the pod, moved to the phase-B model |
| `pilot_cami/probe_cami.py` | pod | short-contig measurements on the CAMI store | 2026-11-09 §6c | as above |
| `pilot_cami/run.sh` | pod | driver then probe, each in a memory scope | 2026-11-09 §6c | as above |

## Per script

Paths are relative to the repository root. `REPO` is the checkout (its `src/` is on `PYTHONPATH`).
`python` means the image's Python inside the prefix above.

### `build_s0.py`

```
python scripts/stress/build_s0.py REPO GOLDEN_DIR ROOT RELEASE [--cluster THREADS] [--pfam PATH] [--ncbifam PATH] [--cpus N] [--workers N] > report.json
```

Builds S0 from the golden set with the repo's own builders and prints a JSON report: host, Python,
manifest hash, per-step timings and results, peak RSS, per-table files/rows/bytes, sha256 of every
file. Stage A (ingest, `genome_input`, gene calling with `--workers` genomes at once, IGRs, tapes)
always runs. `--cluster` adds MMseqs2 clustering, tokens, `token_postings` and `contig_dim`; `--pfam`
/ `--ncbifam` add domain scans (HMM files, `--cpus` threads) before tokens. Card B-18 replaced the
removed pairs step by the postings builder. Writes `ROOT/manifests/release=<R>.json`.

### `build_fixture.py`

```
python scripts/stress/build_fixture.py REPO ROOT RELEASE > report.json
```

Stage A of the micro-S0 fixture (`tests/fixtures/micro_s0`: three genomes, 0.8 Mbp, no network, no
MMseqs2) plus `contig_dim`, then `catalog.build_manifest`. ROOT must not exist. The report has
`build_s0.py`'s shape plus the versions of the packages whose bytes reach the files. The cross-host
check of card B-18: `.github/workflows/cross_host.yml` runs it on a GitHub runner set up as `ci.yml`
is, and the pod runs it in the image; `diff_reports.py` compares the two reports.

### `show_report.py`, `diff_reports.py`

```
python scripts/stress/show_report.py REPORT.json
python scripts/stress/diff_reports.py A.json B.json
```

One line per table (rows, bytes, bytes per row, per gene, per Mbp); and the manifest hash of two
reports, equal or different, with the files that differ.

### `column_bytes.py`

```
python scripts/stress/column_bytes.py MANIFEST_JSON [--root ROOT] [--where KEY=VALUE ...] > bytes.json
```

Reads only Parquet footers. Per table: files, rows, bytes; per column the compressed bytes of its
column chunks. `--where stage=P1` measures one stage of a multi-stage store. Column chunks compress
independently, so a table's size without some columns is its size minus their chunks: B-18 uses this
to correct the import ladder's per-unit sizes, measured before card B-17 removed the legacy columns.

### `step6_export.py`, `hash_fa.py`, `prodigal_run.py`, `step6_compare.py`

```
python scripts/stress/step6_export.py ROOT OUT_DIR
samtools faidx SHARD.fa.gz -r regions_plus.txt -o plus.fa     # outside these scripts
python3 scripts/stress/hash_fa.py plus.fa > plus.tsv
python3 scripts/stress/prodigal_run.py GOLDEN_DIR OUT_DIR
python scripts/stress/step6_compare.py STEP6_DIR ROOT > compare.json
```

The hand-run independent check of the 2026-11-09 file §2. `step6_export.py` reads genes from
`feature_protein`, orders them by (accession, start, end, strand) and names each by its `feature_uid`
(sixth column of `genes.tsv`). `step6_compare.py` compares each shard's indexes with
`pod/<shard name>.fai` / `.gzi`. The same checks run as golden-ring tests
(`tests/r3/test_golden_independent.py`, card D-03); B-18 used those (exit §4).

### `gtdb_counts.py`, `gtdb_postings.py`

```
python scripts/stress/gtdb_counts.py GTDB_DIR > counts.json
python scripts/stress/gtdb_postings.py GTDB_DIR UNITS_JSON > projection.json
```

`GTDB_DIR` has `bac120_metadata_r226.tsv.gz` and `ar53_metadata_r226.tsv.gz` (AGENTS.md section 8
allows fetching only URLs listed in `golden.yaml`, `refdb.yaml` or `sources.yaml`; B-18 used the copy
already made for the import ladder). `gtdb_counts.py` prints genomes, bp, contigs, proteins and
contig-count quantiles for S6, S5, S3, S2C. `gtdb_postings.py` multiplies those unit counts by the
bytes per unit in `UNITS_JSON` (`{"table": {"unit": "bp"|"contig"|"gene", "bytes_per_unit": B,
"source": "..."}}`, or `rows_per_unit` and `bytes_per_row`) for every table, `token_postings` included.
It replaces `gtdb_pairs.py`, which bounded `slot_pairs` rows: pairs are no longer materialized
(AGENTS.md section 3.5), so card B-18 deleted it.

### `step4_keys.py`

```
python scripts/stress/step4_keys.py ROOT SCRATCH RELEASE > keys.json
```

ROOT is a golden store from `build_s0.py ... --cluster N`. The script copies it to `SCRATCH/after`,
adds one ncRNA feature (kind 2, 377 nt, minus strand, 60 nt into the first IGR of at least 500 nt at or
after 3.27 Mb of *E. coli* K-12, as on 2026-11-09) as a second producer, and runs the repo's IGR,
tape, token and postings builders for the tape version of both producers. Every table is then
compared with the original under its registry key (`feature_uid`) and under the slot ordinal, the key
phase B removed. This is the evidence behind AGENTS.md section 3.3.

### `synth_scale.py`

```
python scripts/stress/synth_scale.py REPO OUT COPIES [--shard-bases N] [--postings-width {0,1,2}] [--seed S] > synth.json
```

Copies each micro-S0 contig `COPIES` times with 3% substitutions and ingests them; replicates the
fixture's gene features per copy; writes synthetic clusters (NR90 per 10 copies, NR50 per 100) and a
synthetic family library (one family per fixture protein, 20% without, 3% in one hub); then runs the
repo's tape, token, postings and `contig_dim` builders and builds the catalog. Each key lookup
(sequence by contig, tape of a contig, postings of a token, neighbours of a token) runs with a literal
key and with keys from a joined table, scored by `catalog.profile_query`: row groups read of the total
per scan, the sequence index row groups, and `rg_per_key`. `key_lookups_over_2_rg_per_key` lists any
lookup above AGENTS.md section 4.4's bound. `--postings-width` replaces the old `--pair-width`.
AGENTS.md section 4.2 names this script for the 100x builder memory tests.

### `builder_limits.py`, `macro_scan.py`

```
python scripts/stress/builder_limits.py ROOT SCRATCH > limits.json
python scripts/stress/macro_scan.py LOCATOR_ROOT > macro.json
```

`builder_limits.py`: heap per pair row of `pairs_of_contig` (the test oracle; nothing builds pairs
since B-11) on the largest tape, and `SeqStore` open time and heap on a synthetic 10^6-contig index
written under `SCRATCH` (about 100 MB). `macro_scan.py`: on that index, row groups read by the user
macros `default_digest` / `digest_len` for a literal key, one joined key and 1,000 joined keys. These
macros scan the index views by design (`catalog.macro_sql`); the sequence macros resolve keys through
`bb_seq` instead (card B-14), which `synth_scale.py` measures.

### `pilot_s1/pilot.py`, `pilot_s1/pilot_resume.py`, `pilot_s1/resume.sh`

```
python scripts/stress/pilot_s1/pilot.py SRC_FASTA WORKDIR RELEASE [--limit N] [--threads 16] [--workers 1] [--step-timeout 10800] [--overhead-n 200] [--resume]
python scripts/stress/pilot_s1/pilot_resume.py SRC_FASTA WORKDIR RELEASE [same options]
sh scripts/stress/pilot_s1/resume.sh [PID_OF_EARLIER_PROCESS]
```

Throwaway S1 rate pilot on the INPHARED FASTA (`14Apr2025_genomes.fa.gz`): ingest, genome inputs,
gene calling, IGRs, tapes, MMseqs2 clustering, tokens, postings, `contig_dim`, manifest. `WORKDIR`
gets `root/`, `mmseqs_tmp/`, `manifest.json` and `result.json` (rewritten after each step). Pod only:
long, keep to AGENTS.md section 7.2 and check for the `STOP` file. `pilot_resume.py` is `pilot.py
--resume`. `resume.sh` takes its paths from `PILOT_*` environment variables (see the file) and runs
the resume inside a `systemd-run --scope` with `MemoryMax` and `MemorySwapMax=0`. The 2026-11-09 run
used the pre-phase-B versions; its store is pre-phase-B format.

### `pilot_s1/store_probe.py`

```
python scripts/stress/pilot_s1/store_probe.py ROOT OUT_JSON [--n 1000] [--k 10] [--release R | --stages S1,P1]
```

In a fresh process: `SeqStore` open time and memory, index footers and row groups read on open, and
`--n` random neighbourhoods of +-`--k` genes fetched with and without the block cache. Reads tapes
through `biobox.l3.tokens.read_tapes` of the gene-call tape version, for the stages of `--release`
(its `release_def`) or `--stages`; with neither, every stage directory of that tape version.

Since card B-25 every release-level builder (`cluster.run`, `tokens.run`, `postings.run`, `dim.run`)
and `catalog.build_manifest` takes the release's stage list; the scripts here pass the one stage they
build (`S0`, `S1` or the CAMI stage), and `step4_keys.py` reads it from the release's `release_def`.
Since card B-28 `SeqStore(root, stages)` has no default stage list: `builder_limits.py`, `synth_scale.py` and
`step6_export.py` open it with `["S0"]`, `pilot.py` with `["S1"]`, the CAMI scripts with the pilot's stage,
and `store_probe.py` with the stages of `--release` or `--stages`. `tests/r0/test_stress_scripts.py` parses
every script here with `ast` and fails on a `SeqStore` call without stages or a builder or reader call
with fewer arguments than the function requires, so the next signature change is caught statically.

### `pilot_cami/pilot_cami.py`, `pilot_cami/probe_cami.py`, `pilot_cami/run.sh`

```
sh scripts/stress/pilot_cami/run.sh WORKDIR SRC_FASTA RELEASE
python scripts/stress/pilot_cami/pilot_cami.py SRC_FASTA WORKDIR RELEASE [--threads 8] [--workers 1] [...]
python scripts/stress/pilot_cami/probe_cami.py WORKDIR RELEASE [--threads 8] [--n 1000] [--flank 1000]
```

The short-contig pilot of the 2026-11-09 file §6c, copied into the repo by card B-18 from
`/tmp/if2453/pilot_cami/scripts` on the pod and moved to the phase-B model (postings instead of the
pair steps, part files, the per-stage `seq_locator`, `tape_version` paths, a `SeqStore` that loads no
index). The source is the gold-standard assembly of CAMI II marine sample 0, already on the pod under
`/tmp/if2453/pilot_cami/src/x/`; its download script stays on the pod because its URL is not listed in
`sources.yaml` (AGENTS.md section 8). The 2026-11-09 store at `/tmp/if2453/pilot_cami/work/root` is
pre-phase-B format; B-18 reran the pilot into `/tmp/if2453/b18/cami` (exit file §7).
