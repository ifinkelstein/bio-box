# Stress scripts

These scripts produced the measurements in
[`planning/measurements/2026-11-09_design_stress_test.md`](../../planning/measurements/2026-11-09_design_stress_test.md)
(called "the measurements file" below; section numbers are its sections). Card A-03 documents them;
card B-18 ("Phase-B exit measurement") updates them to the new data model and reruns them. The numbers
live in the measurements file, not here.

Rules that apply to every run (AGENTS.md sections 4, 6, 7):

- No compute on the laptop. Run on the pod (wilkcomp02) through the test image, with the checkout's
  `src/` first on `PYTHONPATH`. `scripts/pod_pytest.sh <name>` mirrors a tree to
  `/tmp/if2453/ci/<name>` and shows the exact prefix (`systemd-run --user --scope -p MemoryMax=...
  nice -n 19 singularity exec --cleanenv --env PYTHONPATH=<tree>/src:<pydeps> --bind /tmp/if2453
  --pwd <tree> <image> python ...`). A script's `REPO` argument is that tree.
- Every script prints usage and exits with `--help` (or `-h`) before doing any work. With too few
  arguments the stdlib-only scripts print the same usage and exit 2.
- Scripts print JSON or text to stdout; redirect it to a file under `/tmp/if2453/`. Nothing here
  writes outside the output directory it is given.
- The CAMI II pilot of section 6c has no driver in this directory (its store is on the pod under
  `/tmp/if2453/pilot_cami`). Section 1 compares a laptop build with a pod build; with the rule above
  the second host must now be another pod run or a CI run, decided at B-18.

Status column: "pre-B" means the script was written for the data model before phase B and card B-18
has to update it. This was found by reading the code, and each reason is stated.

## Index

| Script | Host | Measures | Measurements file | Status |
|---|---|---|---|---|
| `build_s0.py` | pod | build of the golden set, per-table size, cross-host manifest hash | 1, 6b | stage A current; stage B (`--cluster`) broken |
| `show_report.py` | any (stdlib) | compact view of a `build_s0.py` report | 1, 3 | current |
| `diff_reports.py` | any (stdlib) | two `build_s0.py` reports: manifest hash and per-file sha256 | 1 | current |
| `step6_export.py` | pod | exports gene slices from a built store for the independent check | 2 | pre-B |
| `hash_fa.py` | pod (stdlib) | sha256 per record of a `samtools faidx` FASTA | 2 | current |
| `prodigal_run.py` | pod (stdlib + `prodigal`) | Prodigal single mode on the ten golden prokaryotes | 2 | current |
| `step6_compare.py` | any | compares bio-box with samtools, Prodigal and RefSeq CDS | 2 | pre-B |
| `gtdb_counts.py` | any (needs `duckdb`) | unit counts per stage from GTDB R226 metadata | 3 | current |
| `gtdb_pairs.py` | any (needs `duckdb`) | lower bound of `slot_pairs` rows per level | 3 | obsolete by design (see below) |
| `step4_keys.py` | pod | slot-key stability: ordinal key against coordinate key | 4 | pre-B |
| `synth_scale.py` | pod | 1,000x layout and query read counts | 5 | pre-B (slot_pairs section) |
| `builder_limits.py` | pod | heap per pair row and `SeqStore` open cost at 10^6 contigs | 5 (limits table) | pre-B |
| `macro_scan.py` | pod | row groups read by the sequence macros at 10^6 contigs | 5 (limits table) | needs rerun for B-14 |
| `pilot_s1/pilot.py` | pod | S1 rate pilot on INPHARED phages | 6a, 6b | pre-B, `biobox.l3.pairs` is gone |
| `pilot_s1/pilot_resume.py` | pod | same, skipping steps that already finished | 6a | pre-B, as `pilot.py` |
| `pilot_s1/resume.sh` | pod | wrapper that waits for an earlier process, then runs `pilot_resume.py` | 6a | pod paths hard-coded |
| `pilot_s1/store_probe.py` | pod | `SeqStore` open cost and 1,000 neighborhood fetches in a fresh process | 6a, 6c | check at B-18 |

## Per script

Paths are relative to the repository root. `REPO` is the checkout (its `src/` is on `PYTHONPATH`).
`python` means the image's Python inside the prefix above.

### `build_s0.py`

Builds S0 from the golden set with the repo's own builders and prints a JSON report: host, Python,
manifest hash, per-step timings, peak RSS, per-table files/rows/bytes, sha256 of every file.

```
python scripts/stress/build_s0.py REPO GOLDEN_DIR ROOT RELEASE [--cluster THREADS] [--pfam PATH] [--ncbifam PATH] [--cpus N] > report.json
```

Stage A (ingest, gene calling, IGRs, tapes) always runs. `--cluster` adds clustering, tokens and pairs;
`--pfam` / `--ncbifam` add domain scans (HMM files; they need `hmmsearch`, so the pod). Used by
measurements section 1 (cross-host hash, stage B, full build) and section 3 (bytes per unit). Cards:
B-18 (cross-host manifest hash).

Pre-phase-B code: with `--cluster` it imports `biobox.l3.pairs` and calls `pairs.run`, a module card
B-11 removed (`src/biobox/l3` now has `postings.py` and `geometry.py`). Stage A and `--pfam` /
`--ncbifam` do not touch it. B-18 must replace the `pairs` step with the postings builder.

### `show_report.py`

```
python scripts/stress/show_report.py REPORT.json
```

Prints one line per table (rows, bytes, bytes per row, per gene, per Mbp) from a `build_s0.py` report.
Needs only the report's keys, so it works on any host. Used for the per-unit costs of section 3.

### `diff_reports.py`

```
python scripts/stress/diff_reports.py A.json B.json
```

Compares two reports from `build_s0.py` (for example laptop and pod): manifest hash equal or
different, and which files differ. This is the cross-host determinism check of section 1 and the
acceptance of B-18.

### `step6_export.py`

Exports every gene slice of the built golden store for an independent check by `samtools faidx`.

```
python scripts/stress/step6_export.py ROOT OUT_DIR
```

Writes `regions_{plus,minus}.txt`, `expected_{plus,minus}.tsv` (sha256 of the bytes bio-box returns)
and `genes.tsv`. Run on the store built by `build_s0.py`; then copy `OUT_DIR` and the shard to the pod
(section 2 used a laptop-built shard; build on the pod now). Pre-phase-B: reads
`protein_occurrence` and orders by `slot_ordinal`, with no `feature_uid` (AGENTS.md section 3.3).

### `hash_fa.py`

```
samtools faidx SHARD.fa.gz -r regions_plus.txt -o plus.fa     # outside this script
python3 scripts/stress/hash_fa.py plus.fa > plus.tsv
```

Writes `name<TAB>sha256` per record (the part of the name before `/`), stdlib only.
Output goes to `step6_compare.py` as `pod/plus.tsv` and `pod/minus.tsv` under its `STEP6_DIR`.

### `prodigal_run.py`

```
python3 scripts/stress/prodigal_run.py GOLDEN_DIR OUT_DIR
```

Runs the Prodigal binary on the ten golden prokaryotes (genetic code per genome is hard-coded) at
`nice -n 19`; writes `OUT_DIR/<id>.fna` and `.gff` and prints the sha256 and CDS count per genome.
`prodigal` must be on `PATH` (an oracle only; AGENTS.md section 7.1 forbids using pod tools to produce
outputs). Used for the Prodigal row of section 2.

### `step6_compare.py`

```
python scripts/stress/step6_compare.py STEP6_DIR ROOT > compare.json
```

`STEP6_DIR` holds the `step6_export.py` output plus a `pod/` directory with `plus.tsv`, `minus.tsv`,
`shard-0000.fa.gz.fai`, `shard-0000.fa.gz.gzi` (made by samtools) and the Prodigal `*.gff` files, plus
`GCF_000005845.2.gff.gz`. Compares slices, indexes, Prodigal calls and RefSeq CDS (section 2). Reads
`ROOT/l0seq/seq_shards/stage=S0/shard-0000.fa.gz`. Pre-phase-B only through its input
(`step6_export.py` output).

### `gtdb_counts.py`

```
python scripts/stress/gtdb_counts.py GTDB_DIR > counts.json
```

`GTDB_DIR` has `bac120_metadata_r226.tsv.gz` and `ar53_metadata_r226.tsv.gz` (listed in
AGENTS.md section 8 allows fetching only URLs listed in `golden.yaml`, `refdb.yaml` or
`sources.yaml`). Prints genomes, bp, contigs, proteins and contig-count quantiles for S6, S5, S3, S2C
and complete genomes. Used by the unit counts of section 3. It is a metadata query; run it on the
pod like the others.

### `gtdb_pairs.py`

```
python scripts/stress/gtdb_pairs.py GTDB_DIR > pairs.json
```

Lower bound (and the upper bound of 20) on `slot_pairs` rows per level for S2C, S5 and S6. Section 3
used it for the pair-row projection. AGENTS.md section 3.5 now forbids materializing pairs, so B-18
should replace it by a `token_postings` projection (one row per slot and token) rather than rerun it.

### `step4_keys.py`

```
python scripts/stress/step4_keys.py ROOT > keys.json
```

On the *E. coli* K-12 chromosome of a built golden store: adds one ncRNA slot, splits the longest gene,
and counts which slot, IGR and pair rows survive under the ordinal key and under the coordinate key.
This is the evidence behind AGENTS.md section 3.3. Pre-phase-B: the "current key" arm is
`(contig_uid, slot_ordinal)`; it imports the pair enumeration from `tests/oracles/pairs_oracle.py`
(already moved there by B-11) and reads `slot_ordinal` from `protein_occurrence`. B-18 changes it to
check `feature_uid` stability (expected: 0 rows changed) and drops the pair arm.

### `synth_scale.py`

```
python scripts/stress/synth_scale.py REPO OUT COPIES [--shard-bases N] [--pair-width {0,1,2}] [--seed S] > synth.json
```

Copies each micro-S0 fixture contig `COPIES` times with 3% substitutions, ingests them, writes
synthetic gene, token and tape tables in the repo's layout, and scores canonical queries with
`catalog.profile_query` (row groups and bytes read). Section 5 used `COPIES=1000`. It reports its
own peak RSS (`peak_rss_mb_*` keys); run it on the pod under the memory scope. AGENTS.md section 4.2
names it (as `synth_scale.py`) for the 100x builder memory tests.

Pre-phase-B: the file says so in a comment. The `slot_pairs` section (lines from "slot_pairs" down)
writes a table the registry no longer has and calls `tests/oracles/pairs_oracle.py`; the query list
still uses `pairs_of`, `partners_of`, `q1_count` through the catalog. Gene tables are keyed by
`slot_ordinal`. B-18 moves it to `token_postings` and `feature_uid`.

### `builder_limits.py`

```
python scripts/stress/builder_limits.py ROOT SCRATCH > limits.json
```

Python heap per pair row from `pairs_of_contig` (the oracle) on the largest tape, and `SeqStore`
open time and heap on a synthetic 10^6-contig locator written under `SCRATCH`. These two numbers
produced the "pair builder" and "`SeqStore` open" rows of the limits table in section 5, and the
evidence for AGENTS.md section 4.1. The pair half is moot after B-11; the open-time half is the
baseline for cards B-13 and B-18. `SCRATCH` gets about 100 MB.

### `macro_scan.py`

```
python scripts/stress/macro_scan.py LOCATOR_ROOT > macro.json
```

`LOCATOR_ROOT` is the synthetic store `builder_limits.py` wrote (`SCRATCH/locator_1e6`). Builds a
manifest and catalog (written next to it as `locator.duckdb`, in the parent directory of
`LOCATOR_ROOT`) and counts row groups scanned by `default_digest(...)` / `digest_len(...)` for a literal key,
one joined key and 1,000 joined keys (AGENTS.md section 4.4). Section 5 limits table. B-14 and
B-18 rerun it to show the joined case now reads O(1) row groups.

### `pilot_s1/pilot.py`

```
python scripts/stress/pilot_s1/pilot.py SRC_FASTA WORKDIR RELEASE [--limit N] [--threads 16] [--step-timeout 10800] [--overhead-n 200] [--pair-sample 0.02]
```

Throwaway S1 rate pilot on the INPHARED FASTA (`14Apr2025_genomes.fa.gz`): ingest, gene calling with
`pyrodigal-gv`, IGRs, tapes, MMseqs2 clustering, tokens, pair-row counting plus a 2% pair sample,
manifest. `WORKDIR` gets `root/`, `mmseqs_tmp/`, `pairs_sample/`, `manifest.json` and `result.json`
(rewritten after each step). Needs `pyrodigal_gv` and `mmseqs`. Pod only: long (hours), uses
`--threads` threads, so keep to the thread and CPU-minute limits of AGENTS.md section 7.2 and check
for the `STOP` file. Sections 6a and 6b. Pre-phase-B: the `pairs_count` and `pairs_sample` steps import
`biobox.l3.pairs` (`K_MAX`, `window`, `BYTES_PER_ROW`, `pairs_of_contig`), which no longer exists, and
write the removed `slot_pairs` table. Everything up to tokens still matches the code. Section 6b says
its store is "pre-phase-B format".

### `pilot_s1/pilot_resume.py`

Same arguments as `pilot.py`; identical except that it reads an existing `WORKDIR/result.json` and
skips steps that finished without error. Same pre-phase-B status.

```
python scripts/stress/pilot_s1/pilot_resume.py SRC_FASTA WORKDIR RELEASE [--limit N] [--threads 16] [--step-timeout 10800]
```

### `pilot_s1/resume.sh`

```
sh scripts/stress/pilot_s1/resume.sh [PID_OF_EARLIER_PROCESS]
```

Waits for `PID` to exit, refuses to start if the `STOP` file exists in the data root, then runs
`pilot_resume.py` in the image at `nice -n 19` with 16 threads and a 12,600 s step timeout. Paths
(`/tmp/if2453/pilot_s1`, the image digest, `/stor/scratch/Wilke/if2453/all-the-data/STOP`) are
hard-coded for the 2026-11-09 run; it does not use `systemd-run`, so the memory scope of AGENTS.md
section 7.2 is not applied. Edit paths and add the scope before any reuse. Section 6a.

### `pilot_s1/store_probe.py`

```
python scripts/stress/pilot_s1/store_probe.py ROOT OUT_JSON [--n 1000] [--k 10]
```

In a fresh process: `SeqStore` open time and memory, index footers and row groups read on open, and
`--n` random neighborhoods of +-`--k` genes fetched with and without the block cache. Runs on a store
built by the pilots (sections 6a and 6c). It reads tapes through `biobox.l3.tokens.read_tapes`,
which still exists; check its keys against the new tape schema at B-18.
