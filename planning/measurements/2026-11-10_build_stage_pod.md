# First real build with `build-stage` on the pod (card C-10, wilkcomp02, 2026-11-10)

Two stages and two releases built in a fresh data root on local disk, each `build-stage` launched
detached as `scripts/pod_run.sh` → `bio-box build-stage` → one `run-step` per step, with the real
wilkcomp02 profile. Every number is *measured* on wilkcomp02 on 2026-11-10 (CST) unless labelled
*derived* or *estimate*. Sources: a run record is `/tmp/if2453/c10/root/runs/stage=<S>/<step>.json`
(wall, CPU-minutes and peak RSS from `os.wait4` in `run-step`); a unit log is
`/tmp/if2453/c10/tmp/logs/<unit>.log`; the scan shards' numbers are in `/tmp/if2453/c10/pfam_P1/shards.json`
(`os.wait4` in `pfam_shards.py`) and `scan/shard-<k>.json`. Nothing was written under `/stor`.

## Code that ran

`plan/task-cards` at `121503d` (C-04 merged, C-09 merged) plus this card's helper scripts, no change
under `src/` or `tests/`, mirrored once to `/tmp/if2453/ci/C-10` (`sh scripts/pod_pytest.sh C-10 tests/r0
-k nothing_matches_this`) and used for the whole build through
`--env PYTHONPATH=/tmp/if2453/ci/C-10/src:/tmp/if2453/pilot_s1/pydeps`. While this build ran,
`plan/task-cards` moved on: cards B-24 and C-11 make `build-stage` shard the domain scan itself
(`--scan-jobs`, `--scan-cpus`). This build kept the mirrored tree to the end, so **the sharded-scan
path of B-24/C-11 was not exercised here**: the scan of P1 was refused by the one-process estimate and
run by hand (below), as the card foresaw.

## Setup

- Image `/stor/home/if2453/bio-box/images/bio-box-test_74b3e05.sif`,
  `sha256:c30be58557d756095346308c58cc3098af4bbe66404f67ab87be3f64f147c898` (computed by `pod_run.sh`
  at each launch; every run record carries it in `image_digest`).
- Profile `/tmp/if2453/c10/profile.yaml`: a copy of `profiles/wilkcomp02.yaml` with one key changed,
  `tmp: /tmp/if2453/bio-box` → `tmp: /tmp/if2453/c10/tmp` (staging, scratch and unit logs). `durable`
  and `bulk` unchanged: nothing `build-stage` runs reads them (only `doctor`). `mem_scope: 200G`,
  `small_file_floor_bytes: 67108864`, thread caps 40/56, `cpu_min_per_process: 300` and
  `require_mem_cap: true` unchanged.
- Data root `/tmp/if2453/c10/root` (did not exist; the first step created it). Pfam 38.2 read in place
  from `/stor/scratch/Wilke/if2453/all-the-data/l0/refdb/pfam/38.2/Pfam-A.hmm.gz`, bound read-only.
- Pod before the first launch (07:15 CST; `top -b -n 1`, `ps -eo user,pcpu`): load 27.4, every process
  above 12% CPU was our own account (other agents' pytest scopes and unit `b31-after-1000`, nice 19);
  no other user computing. `df -h /tmp`: 882 G free (863 G at the end). 1007 G RAM, 1117 G swap.
- Every `build-stage` launch had this form (unit, stage, release, `--stages`, inputs and thread count
  vary; `--threads` was 16 for `c10-S1a-*` and `c10-P1-1`, 8 from `c10-P1-2` on, see "Thread budget"):

```
cd /tmp/if2453/c10 && bash /tmp/if2453/ci/C-10/scripts/pod_run.sh --profile /tmp/if2453/c10/profile.yaml \
  --image /stor/home/if2453/bio-box/images/bio-box-test_74b3e05.sif --threads 16 --bind /tmp/if2453 \
  --bind /stor/scratch/Wilke/if2453/all-the-data/l0/refdb/pfam/38.2:/stor/scratch/Wilke/if2453/all-the-data/l0/refdb/pfam/38.2:ro \
  --env PYTHONPATH=/tmp/if2453/ci/C-10/src:/tmp/if2453/pilot_s1/pydeps \
  c10-S1a-1 build-stage --profile /tmp/if2453/c10/profile.yaml --root /tmp/if2453/c10/root \
  --stage S1a --release rC.1 --stages S1a \
  --genomes /tmp/if2453/c10/inputs/S1a/genomes.json --objects /tmp/if2453/c10/inputs/S1a/objects.json \
  --library pfam 38.2 /stor/scratch/Wilke/if2453/all-the-data/l0/refdb/pfam/38.2/Pfam-A.hmm.gz gathering \
  --threads 16
```

  P1: `--stage P1 --release rC.2 --stages S1a,P1` and `inputs/P1/*.json`.

## Inputs (`scripts/stress/build_stage_pod/make_inputs.py`, run in the image in a 8G scope)

| Stage | Source | Records | Genomes | Contigs | Description |
|---|---|---|---|---|---|
| `S1a` | first 2,000 records of INPHARED `/tmp/if2453/pilot_s1/src/14Apr2025_genomes.fa.gz` (the ladder's S1 source, in file order), copied byte for byte to `inputs/S1a/phage_first2000.fa` (63,357,113 bytes) | 2,000 | 1,947 | 1,947 | each record its own genome, virus, code 11, linear (as the ladder's S1) |
| `P1` | the ladder's 250 GTDB R226 genomes (`/tmp/if2453/ladder/genomes/P1/lock.json`; files read in place) | 539 | 250 | 539 | one genome per accession, prokaryote, complete, linear, code from `selection/stage_P1.tsv` (11: 244, 4: 6) |

**Surprise 1.** 53 of the first 2,000 INPHARED records repeat an earlier record's sequence (equal
`ids.contig_uid`), so S1a has 1,947 genomes, not 2,000. `genome_input` refuses a contig listed in two
genomes, so `make_inputs.py` lists each sequence once, under its first record. The ladder did the same
implicitly (one genome per default contig of the index: 27,258 genomes from 34,076 records, the latter
*measured* by `zcat ... | grep -c '^>'`). Ingest read all 2,000 records (`records: 2000,
new_sequences: 1947`, log `c10-S1a-1`). P1's genome inputs equal the ladder's
(`/tmp/if2453/ladder/results/B/genome_input.json`: 250 genomes, 539 contigs, no duplicate, same codes).

## What ran, in order (CST)

| Time | Unit | What happened | Exit |
|---|---|---|---|
| 07:21 | `c10-S1a-1` | inputs, ingest published; `STOP` created at 07:21:57 during genecall (exercise a) | 75 |
| 07:25 | `c10-S1a-2` | `STOP` still present: nothing run | 75 |
| 07:26 | `c10-S1a-3` | `STOP` removed; inputs, ingest skipped; genecall to catalog of `rC.1` | 0 |
| 07:38 | `c10-S1a-4` | finished build relaunched (exercise b): every step skipped | 0 |
| 07:39 | `c10-P1-1` | inputs to tapes published; refused before the Pfam scan | 2 |
| 07:46 to 09:58 | `c10-pfam-P1`, `c10-pfam-P1b` | the P1 scan in 19 shard processes, outside `build-stage` (below) | 0 |
| 09:58 | `c10-P1-pfam-merge` | `run-step --stage P1 --step domains-pfam-38.2 -- python pfam_merge.py ...` published the scan | 0 |
| 09:59 | `c10-P1-2` | P1 steps skipped; unit stopped 120 s into `rC.2.cluster` (exercise c) | (143, see defect 1) |
| 10:03 | `c10-P1-3` | resumed at `rC.2.cluster`; cluster to catalog of `rC.2` | 0 |
| 10:33 | `c10-P1-4` | finished build relaunched (exercise b): every step skipped | 0 |

## Acceptance

**Both releases verify.** Command (in the image, 8G scope):
`python -c "... m = catalog.read_manifest(catalog.manifest_path(r, rel)); catalog.verify_manifest(m, r)"`
(there is no `bio-box release verify`):

```
rC.1 ['S1a'] 25 files 4 artifacts efbcf11464c93657b0d5ea0c4806ea8f problems: []
rC.2 ['S1a', 'P1'] 41 files 7 artifacts 1fa3a73a827833946cbb8210bf47cc9a problems: []
```

**Every step has a run record whose `mem_cap` is inherited from its unit.** All 24 records say
`"how": "inherited"`, `memory_max_bytes: 214748364800` (200 GiB, the profile's `mem_scope`), with
`limited_at` the cgroup of the unit that ran it (`.../app.slice/c10-<unit>.service`; table below).
No record says `opened` or `none`. C-09 held at the real 64 MiB floor: `rC.2`'s manifest and catalog
were published beside `rC.1`'s, and no step was refused for small files.

### (a) `STOP` while a build runs

`watch_then.sh` (unit `c10-watch-stop`) created `/tmp/if2453/c10/root/STOP` at 07:21:57, 20 s after
`staging/S1a.genecall` appeared. Log of `c10-S1a-1`:

```
run-step: /tmp/if2453/c10/root/STOP appeared; not publishing S1a.genecall, staging left at /tmp/if2453/c10/tmp/staging/S1a.genecall
build-stage: S1a.genecall exited 75; stopping here; the next call resumes at this step
pod_run: c10-S1a-1 exited 75
```

Relaunched with `STOP` present (`c10-S1a-2`): `build-stage: /tmp/if2453/c10/root/STOP exists; not
starting stage S1a`, exit 75, root unchanged. `STOP` removed, relaunched (`c10-S1a-3`): inputs and
ingest skipped, genecall run again, the build finished (exit 0).

**Difference from the fixture test.** `test_a_stop_file_stops_before_the_next_step_and_a_later_call_continues`
creates `STOP` *between* steps (hook `after:genecall`), so genecall is published and the build stops
before igr. On real data a `STOP` almost always lands *inside* a step; `run-step` then lets the command
finish and publishes nothing (its documented rule, C-02 step 6), so the step's whole work is lost and
repeated: here genecall's 57.7 s and 4.98 CPU-min (the rerun's record). The build does stop before
the next step with exit 75 and resumes correctly; what is discarded is the step in flight.

### (b) A finished build relaunched

`c10-S1a-4` and `c10-P1-4`: every one of the 12 commands printed `run-step: <stage>.<step> skipped:
run record and outputs present`, `features` was passed over, exit 0, and `sha256sum` of every file
under the root, run records included (43 files after S1a, 84 after P1), was identical before and after (`cmp` of
`/tmp/if2453/c10/logs/tree_{before,after}_{S1a,P1}-4.sha256`). No staging directory was created.

### (c) Unit stopped in the middle of a long step

`watch_then.sh` (unit `c10-watch-kill`) ran `systemctl --user stop c10-P1-2` at 10:02:03, 120 s after
`staging/rC.2.cluster` appeared (MMseqs2 clustering of 982,049 proteins). The unit stopped within one
second. Log:

```
run-step: SIGTERM received; forwarded to rC.2.cluster, which exited -15; nothing published, staging left at /tmp/if2453/c10/tmp/staging/rC.2.cluster
build-stage: rC.2.cluster exited 143; stopping here; the next call resumes at this step
```

No process of the build was left (`pgrep -af c10/root` empty), no run record or `.incoming` for
`rC.2`. `c10-P1-3` skipped the six P1 steps, removed the stale staging, ran `rC.2.cluster` from the
start and finished every release step (exit 0). Two things differ from what the fixture tests and
`pod_run.sh` promise: defects 1 and 2 below.

## Expected stop: the P1 scan, and what was done outside `build-stage`

`c10-P1-1` exited 2 after tapes, verbatim:

```
build-stage: stopped before P1.domains-pfam-38.2: the scan of 294,516,869 residues is estimated at 2,150.0 CPU-minutes in one process (7.3 CPU-min per million residues, measured for Pfam on wilkcomp02, 2026-11-10), over the profile's cpu_min_per_process of 300 (AGENTS §7.2). The earlier steps are published; a scan in shard processes is card B-24
```

S1a's scan was not refused: 14,266,576 residues, estimate 104.1 CPU-min (*derived*, `stage.scan_cpu_min`);
it ran in one process with 16 threads, 114.14 CPU-min measured (+9.6% on the estimate).

Done outside `build-stage` for P1, the ladder's procedure, and nothing else:

1. `pfam_plan.py /tmp/if2453/c10/root P1 pyrodigal-3.7.1-gv-0.3.2-p2 /tmp/if2453/c10/pfam_P1` (image,
   16G scope): the 918,765 proteins of the `protein_seq` files listed in `runs/stage=P1/genecall.json`
   (all of the stage's proteins, which `bio-box stage domains` would scan; the ladder scanned only the
   918,525 not in S1), sorted by uid, cut into 19 shards of at most 16,000,000 residues (18 of about
   16.0 M, one of 6.5 M) under `/tmp/if2453/c10/pfam_P1/input/`. Nothing written in the root.
2. `pfam_shards.py` in units `c10-pfam-P1` (MemoryMax=200G, MemorySwapMax=0, CPUQuota=3200%, nice 19)
   and then `c10-pfam-P1b` (CPUQuota=1600%): each shard its own process, `nice -n 19 singularity exec
   ... python scan_shard.py shard-<k>.parquet Pfam-A.hmm.gz gathering scan/shard-<k>.parquet --cpus 4`
   (the ladder's `scan_shard.py`, unchanged, calling `domains.scan`), outputs under `pfam_P1/scan/`.
   At most 8 shards (32 threads) at a time until 07:50, then at most 4 (16 threads); see "Thread budget".
3. `run-step --profile /tmp/if2453/c10/profile.yaml --root /tmp/if2453/c10/root --stage P1 --step
   domains-pfam-38.2 --threads 1 -- python .../pfam_merge.py /tmp/if2453/c10/pfam_P1 pfam 38.2`,
   launched by `pod_run.sh` as unit `c10-P1-pfam-merge`. `pfam_merge.py` checks the plan against the
   gene-call record and each shard's input sha256 and protein count, checks that the shards
   concatenated in order are in `domains.scan`'s row order, and writes into `$BIOBOX_STAGING` exactly
   what `domains.run` writes: `l1/domains/library=pfam-38.2/stage=P1/part-0.parquet` through
   `io.write_parts` (1,552,653 rows, 80,777,664 bytes) and the step's two ledger rows (n_in 918,765,
   750,383 with a hit). `run-step` published both (`run-step: P1.domains-pfam-38.2 published 2 files`)
   and wrote the step's run record, so the step's layout allowed publishing; nothing was written into
   the root by hand. That record's `argv` names `pfam_merge.py`, and its wall and CPU (7.2 s, 0.10
   CPU-min) are the merge's alone; the scan's CPU is in `shards.json`.
4. `build-stage` relaunched (`c10-P1-2`, `c10-P1-3`): it found the scan step's record and skipped it.

Shards (*measured*, `shards.json`): 19 ok, 2,343.74 CPU-min in all (wall of each shard 685 s to
1,980 s, 4 threads), largest 131.46 CPU-min, peak RSS at most 636 MB; 7.96 CPU-min per million
residues (*derived*: 2,343.74 / 294.52), 9% above the ladder's 7.30. Every shard stayed under the 300
CPU-minute budget.

**Surprise 2.** The first shard launch (07:46:51) failed: all 19 shards exited 2 within 2 s with
`usage: scan_shard.py [-h] [--seed SEED] [--b B] [--c C] gtdb out` — the usage of the ladder's
`select.py`. Python puts the script's directory first on `sys.path`, and
`planning/measurements/2026-11-10_import_ladder/scripts/select.py` shadows the standard module
`select`, which `subprocess` imports. A copy of `scan_shard.py` (sha256 `328037bf...`, identical) in
`/tmp/if2453/c10/bin/` ran (defect 3).

**Pfam rows against the ladder.** Every one of the ladder's 1,552,451 P1 domain rows (`rL.2p`'s file,
proteins new in P1) is in this build's file with equal values in every column (`EXCEPT ALL`: 0 only in
the ladder); the 202 rows only here belong to 129 proteins that the ladder did not scan (they were
also in its S1). So the sharded scan and merge reproduce the ladder's sharded scan row for row.

## Thread budget

S1a (`c10-S1a-1` to `-4`) and `c10-P1-1` ran with `--threads 16` (pod_run `--threads 16`, CPUQuota
1600%). At 07:50 CST the coordinator lowered our share to at most 16 threads until 20:00 (the day
window; our account's total was about 43 CPUs with four other agents' pytest). The eight running
shards (32 threads, started 07:47:53) were left to finish; no new shard was started from 07:50:32 (a
`STOP` file in the root, text "C-10 throttle ...", which `pfam_shards.py` honours; nothing else ran in
the root then) until the first eight ended at 08:22:35; the STOP was removed and the remaining 11
shards ran 4 at a time from 08:22:48 to 09:57:58. **The scan's wall time (07:47:53 to 09:57:58) is
therefore not a rate.** `c10-P1-2` to `-4` ran with `--threads 8` (CPUQuota 800%), so `rC.2.cluster`'s
wall time is at 8 threads.

## Per-step numbers (measured, wilkcomp02, 2026-11-10, from the run records `runs/stage=<S>/<step>.json`)

Threads is the `run-step` thread count. Every `mem_cap` is `inherited`, 214,748,364,800 bytes, at the
named unit.

| stage | step | started (UTC) | wall s | CPU-min | peak RSS MiB | threads | unit |
|---|---|---|---|---|---|---|---|
| S1a | inputs | 13:21:11 | 2.4 | 0.07 | 323 | 16 | c10-S1a-1 |
| S1a | ingest | 13:21:14 | 19.8 | 0.36 | 322 | 16 | c10-S1a-1 |
| S1a | genecall | 13:26:29 | 57.7 | 4.98 | 483 | 16 | c10-S1a-3 |
| S1a | igr | 13:27:27 | 3.4 | 0.08 | 323 | 16 | c10-S1a-3 |
| S1a | tapes | 13:27:30 | 3.7 | 0.09 | 322 | 16 | c10-S1a-3 |
| S1a | domains-pfam-38.2 | 13:27:34 | 432.2 | 114.14 | 1,015 | 16 | c10-S1a-3 |
| rC.1 | cluster | 13:34:46 | 34.0 | 5.15 | 919 | 16 | c10-S1a-3 |
| rC.1 | tokens | 13:35:20 | 3.7 | 0.09 | 348 | 16 | c10-S1a-3 |
| rC.1 | postings | 13:35:24 | 2.7 | 0.08 | 432 | 16 | c10-S1a-3 |
| rC.1 | dim | 13:35:27 | 2.6 | 0.08 | 323 | 16 | c10-S1a-3 |
| rC.1 | manifest | 13:35:29 | 2.2 | 0.07 | 321 | 16 | c10-S1a-3 |
| rC.1 | catalog | 13:35:32 | 5.1 | 0.07 | 322 | 16 | c10-S1a-3 |
| P1 | inputs | 13:39:29 | 2.6 | 0.07 | 315 | 16 | c10-P1-1 |
| P1 | ingest | 13:39:32 | 268.3 | 4.44 | 1,216 | 16 | c10-P1-1 |
| P1 | genecall | 13:44:03 | 83.2 | 15.58 | 2,539 | 16 | c10-P1-1 |
| P1 | igr | 13:45:28 | 7.6 | 0.16 | 557 | 16 | c10-P1-1 |
| P1 | tapes | 13:45:35 | 7.7 | 0.16 | 737 | 16 | c10-P1-1 |
| P1 | domains-pfam-38.2 (merge only) | 15:58:18 | 7.2 | 0.10 | 736 | 1 | c10-P1-pfam-merge |
| P1 | (scan shards, `shards.json`) | 13:47:53 | see above | 2,343.74 | 636 (max per shard) | 19 x 4 | c10-pfam-P1, -P1b |
| rC.2 | cluster | 16:03:28 | 1,575.2 | 195.08 | 4,625 | 8 | c10-P1-3 |
| rC.2 | tokens | 16:29:45 | 19.4 | 0.33 | 1,061 | 8 | c10-P1-3 |
| rC.2 | postings | 16:30:04 | 13.3 | 0.24 | 999 | 8 | c10-P1-3 |
| rC.2 | dim | 16:30:18 | 2.8 | 0.06 | 321 | 8 | c10-P1-3 |
| rC.2 | manifest | 16:30:21 | 7.1 | 0.13 | 632 | 8 | c10-P1-3 |
| rC.2 | catalog | 16:30:28 | 7.4 | 0.06 | 322 | 8 | c10-P1-3 |

Not in any record: the discarded S1a genecall of `c10-S1a-1` (STOP) and the first 120 s of
`rC.2.cluster` (killed); `run-step` records only published steps. No step was over the 300 CPU-minute
budget (`over_cpu_budget: false` in all 24). Printed by `python3
scripts/stress/build_stage_pod/runs_table.py /tmp/if2453/c10/root`.

Estimates in `stage.py` against these records: ingest of P1 4.60 CPU-min (*derived*, 1,029,384,352
bytes) against 4.44 measured (−3.5%); ingest of S1a 0.28 CPU-min (63,357,113 bytes) against 0.36
measured (**+29%, over the 20% of AGENTS §6**; the fixed start-up cost of a step, 0.07 CPU-min for
`inputs` here, dominates a 63 MB stage; harmless against a 300-minute gate, reported as finding 6);
scan of S1a 104.1 against 114.14 (+9.6%).

Against the ladder (`planning/measurements/2026-11-10_import_ladder/README.md`, *measured* there, 16
threads both): P1 ingest 268.3 s here against 253.5 s; P1 gene calling 83.2 s against 86.8 s. The
cluster runs are not comparable (different protein sets and thread counts).

## Stage P1 against the ladder's `rL.2`

`compare_ladder.py /tmp/if2453/c10/root rC.2 /tmp/if2453/ladder/root rL.2 P1` (image, files from the
two manifests; the ladder root was only read). Rows per table of stage P1:

| table | `rC.2` | `rL.2` |
|---|---|---|
| contig_index, contig_occurrence, contig_tapes, genome_input, seq_locator | 539 each | 539 each |
| feature, feature_protein, protein_occurrence | 927,920 each | 927,920 each |
| protein_seq (distinct proteins) | 918,765 | 918,765 |
| igr | 764,935 | 764,935 |
| l1_domains | 1,552,653 | (none in `rL.2`; 1,552,451 in `rL.2p`, above) |
| l0_manifest | 250 | (none: per-stage since card B-22) |
| slot_protein | (none: removed by card B-17) | 927,920 |
| ledger | 1,330 | 1,328 |

Not a sample but every gene: the (contig_uid, start, end, strand, kind) rows of `feature` are the same
set on both sides (0 only here, 0 only in the ladder) and every matched gene has the same
`feature_uid`; every IGR (contig_uid, start, end, strand) matches; the protein uids are the same set
and every sequence is equal. The first five genes by sort key are identical on both sides (contig
`00a9c53b...`, [0, 1416), [1900, 3025), [3089, 3986), [3996, 5286), [5278, 5905), all + strand, kind
1). Output: `/tmp/if2453/c10/logs/compare_P1_rL2.json`.

`rC.2.cluster` aligned 982,049 proteins = S1a's 63,284 + P1's 918,765 (no protein shared between the
two stages); 505,045 NR50 and 909,655 NR90 clusters (log `c10-P1-3`).

## Defects and surprises

1. **`pod_run.sh`'s last log line is missing after `systemctl --user stop`.** Its header says the log's
   last line is `pod_run: <unit> exited <code>`. Stopping the unit sends SIGTERM to every process in
   the unit's cgroup, the wrapping `/bin/sh -c` included, so the line is never written: `c10-P1-2.log`
   ends with `run-step: P1.domains-pfam-38.2 skipped: ...`, and `build-stage`'s exit code (143) is only
   in its own stderr line. Reproduce: launch any build with `pod_run.sh`, `systemctl --user stop <unit>`
   during a step, `tail -n 1 <log>`. A poller that waits for that line waits forever.
2. **A killed step leaves its scratch behind.** After the stop in (c), `run-step`'s rerun removed the
   stale staging, but the clustering scratch `/tmp/if2453/c10/tmp/scratch/tmp60hywwi5` (1.2 G) is still
   there; nothing removes scratch of a killed step. Reproduce: stop a unit during `rC.2.cluster`,
   `du -sh <profile tmp>/scratch/*` after the resumed build. Left in place as evidence.
3. **The ladder's `scan_shard.py` cannot run from its directory**: `select.py` beside it shadows the
   standard `select` module (Surprise 2). Reproduce: `python
   planning/measurements/2026-11-10_import_ladder/scripts/scan_shard.py a b c d` prints `select.py`'s
   usage and exits 2.
4. **The unit log is not a timeline.** `run-step`'s `published` / `skipped` lines go to `build-stage`'s
   stdout, block-buffered when redirected to the log, so they all appear at the end, after every step's
   JSON and after `build-stage`'s stderr lines (see the `c10-S1a-1` log above: `S1a.inputs published`
   is printed after `S1a.genecall exited 75`). Watching "the log's last line" shows no progress during
   a build. Reproduce: any `c10-*` log.
5. **A `STOP` inside a step discards that step's finished work** (difference (a) above). Documented
   `run-step` behaviour, not covered by a fixture test with `STOP` mid-step.
6. **Ingest estimate misses by 29% on a small stage** (0.28 against 0.36 CPU-min, S1a); within 4% on
   P1. Not a gate problem; reported as AGENTS §6 asks.
7. No `bio-box` command verifies a release; verification used `catalog.verify_manifest` directly.
8. The `S1a` stage has 1,947 genomes, not 2,000 (Surprise 1).

## Helper scripts (`scripts/stress/build_stage_pod/`)

`make_inputs.py` (the `--genomes`/`--objects` files), `pfam_plan.py`, `pfam_shards.py` (host python3,
shard launcher), `pfam_merge.py` (run as the `run-step` command), `watch_then.sh` (the `STOP` and
unit-stop triggers), `runs_table.py` (the table above), `compare_ladder.py` (the comparison above).
Work products on the pod: `/tmp/if2453/c10/{inputs,pfam_P1,logs,tmp/logs}`.
