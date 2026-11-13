# Split clustering reproduces `rN.1b` (card B-42 step 1)

Card B-42 (1): recluster release `rN.1b` (stages S1 and P1; 2,053,075 proteins) in parts with the current code,
through `build-stage` and `run-step` as a release build does, into a scratch root that reads the stage files
without changing them, and compare with the published clusters row for row. Every number below is *measured* on
wilkcomp02 (`wilkcomp02.ccbb.utexas.edu`), 2026-11-12 (CST; run records in UTC, 2026-11-13), unless labelled
*derived* or *estimate*.

**Result: equal, and byte-identical.** Clustered in parts with the default cut (run A, one part per call) and with
every per-query call forced into three parts (run B), `cluster_member` and `cluster_member_aln` of release `rN.1b`
are the published files byte for byte (same sha256), so equal row for row in file order; no member moves at either
level. Every step's run record is under 300 CPU-minutes (largest 92.72, `cluster.nr90.12.align`); none says
`over_cpu_budget`.

## Setup

- Code: branch `task/B-42` at `be9c198` (the tip of `plan/task-cards` on 2026-11-12), synced to
  `/tmp/if2453/ci/B-42` by `sh scripts/pod_pytest.sh B-42 tests/r0/test_runbook.py` (38 passed). Run records name
  it as `code_id bf0eedd6…` (card C-23; the synced tree is not a git checkout, so `code_commit` is null).
- Image `/stor/home/if2453/bio-box/images/bio-box-test_74b3e05.sif`, `sha256:c30be585…c898` (the digest in every
  run record), MMseqs2 `18.8cc5c` (the plan record's `mmseqs2`), `PYTHONPATH=/tmp/if2453/ci/B-42/src:/tmp/if2453/pilot_s1/pydeps`.
- The data root `/tmp/if2453/r10/root` was bound read-only into every unit
  (`--bind /tmp/if2453/r10/root:/tmp/if2453/r10/root:ro`). The scratch roots `/tmp/if2453/b42/view-A` and
  `view-B` were made by `/tmp/if2453/b42/bin/mkview.py SRC DST rN.1b inputs,ingest,genecall,igr,tapes` (the pattern
  of card B-46's `mkview.py`): real directories holding symlinks to exactly the files and artifacts that
  `manifests/release=rN.1b.json` lists outside its own `release=` partitions (36: the stages' `protein_seq`,
  `contig_index`, store shards, features, IGRs, tapes, ledgers, ...), plus the ten stage run records of S1 and P1
  for those five steps (the script checks that every output those records list is among the 36). Nothing of the
  release (its cluster tables, `release_def`, runs, manifest, catalog), no domain file of the running Pfam scan and
  no `STOP` file was linked. After both runs the real root's `runs/stage=rN.1b/` still holds its six records and
  `l2/cluster_plan/` has no `release=rN.1b`.
- Profile: a copy of `profiles/wilkcomp02.yaml` (`/tmp/if2453/b42/profile.yaml`) whose only change is
  `tmp: /tmp/if2453/b42/tmp`, so run-step's staging is not beside the running build's.
- Launch (`/tmp/if2453/b42/bin/b42run.sh`): `scripts/pod_run.sh` of the synced tree, `--mem 64G` (MemoryMax, no
  swap), `--threads N` (CPUQuota N×100%), nice 19, then
  `build-stage --profile /tmp/if2453/b42/profile.yaml --root <view> --stage P1 --release rN.1b --stages S1,P1
  --genomes /tmp/if2453/r10/inputs/P1/genomes.json --objects /tmp/if2453/r10/inputs/P1/objects.json
  --tmp /tmp/if2453/b42/scratch-<run> --threads N [...]`. The five P1 stage steps were passed over
  (`run-step: P1.<step> skipped: run record and outputs present`); then clustering, tokens, postings, dim, manifest,
  catalog of `rN.1b` in the scratch root.
  - Run A, unit `b42-A`, N = 12, no cluster option: `build-stage` chose parts by itself (`clustering 2,053,315
    proteins is estimated at 284.1 CPU-minutes in one process, over 80% of the profile's cpu_min_per_process of
    300 ...; it runs in parts`; 2,053,315 is the stages' `protein_seq` rows, 240 shared proteins counted twice).
    Started 20:05:16 CST at an account total of about 33 CPUs (`ps -u if2453 -o pcpu=` summed), night cap 56.
  - Run B, unit `b42-B`, N = 8, `--cluster-parts 3 --scan-jobs 2 --scan-cpus 4`: every per-query call in three
    parts (two at a time, four threads each). Started 20:42:22 CST; the account was at 46 CPUs (a 5-second `top`
    sample), so 8 threads, not 12, kept it inside 56. At 20:46 the synced tree was replaced by this card's step-2
    code (a test sync, by mistake while B ran): the 41 part steps started after that ran with it (their records
    carry `code_id 68a67760…`; the rest of A and B carry `bf0eedd6…`). That code changed only `src/biobox/stage.py`'s
    release-level check and the reason given to the one-process cluster step, neither of which a part step runs;
    `src/biobox/l2/cluster.py` was the same file in both, and the outputs are byte-identical anyway.
- Disk: a sampler unit (`/tmp/if2453/b42/bin/dusample.sh`, `du -sB1` every 15 s) recorded the allocated bytes of the
  work directory `<scratch>/biobox-cluster/`, of the whole scratch directory and of run-step's staging.
- Comparison (`/tmp/if2453/b42/bin/compare.py`, in the image, unit `b42-cmp-*`, 2 CPUs).

## What was compared, and how

For each of `l2/cluster_member/release=rN.1b/` and `l2/cluster_member_aln/release=rN.1b/`, published
(`/tmp/if2453/r10/root`) against reclustered (the view):

- the files: names, bytes, sha256, row groups;
- the Arrow schemas (`Schema.equals`) and the tables row for row in file order (pyarrow `Table.equals`: every column,
  every row, the same order);
- rows present on one side only, as multisets (DuckDB `EXCEPT ALL`, both ways);
- members whose cluster or representative differs (`cluster_member`: `cluster_uid`, `rep_uid`), or whose alignment
  row differs (`cluster_member_aln`: `rep_uid`, coordinates, `pident`, `alnlen`, `cigar`), per level, by a full join
  on (`member_uid`, `level`);
- per level: rows, distinct members, distinct representatives, distinct clusters.

`member_cluster` exists only in the reclustering (card B-27 came after `rN.1b` was built): it was checked to hold
exactly the new `cluster_member` rows (`EXCEPT ALL` both ways on its four columns) and to be sorted by
(`level`, `member_uid`) in file order.

| | published `rN.1b` | run A (default cut) | run B (3 parts per call) |
|---|---|---|---|
| `cluster_member` file | 148,567,302 bytes, sha256 `8bc4461d…c709`, 34 row groups | the same bytes | the same bytes |
| `cluster_member_aln` file | 140,592,805 bytes, sha256 `709d9ba1…4fd1`, 34 row groups | the same bytes | the same bytes |
| rows (each table) | 4,106,150 | 4,106,150 | 4,106,150 |
| row for row equal / schema equal | | yes / yes | yes / yes |
| rows on one side only | | 0 and 0 | 0 and 0 |
| members that move (either level) | | 0 | 0 |
| NR90: members, clusters (= representatives) | 2,053,075; 1,395,712 | the same | the same |
| NR50: members, clusters (= representatives) | 2,053,075; 761,833 | the same | the same |
| `member_cluster` | none | 198,251,815 bytes, `8d1c2ed8…109f`; equal to `cluster_member` as a set; 0 rows out of order | the same bytes as run A |

## CPU, wall time, memory, disk

Peak memory is each run record's `peak_rss_bytes`, source `child_maxrss` in every record (the step command's
process tree, card C-15/C-16). Wall time of the whole clustering is from the plan step's start to the end of
`cluster` (the tables); CPU-minutes are the sum of the clustering steps' records (`cluster.plan`, every module call,
cut, part and join, `cluster`).

| | `rN.1b` published (one process, 2026-11-10) | run A, 12 threads | run B, 8 threads (parts 2 × 4) |
|---|---|---|---|
| clustering steps (run records) | 1 | 60 | 108 (16 calls in 3 parts each, with their cuts and joins) |
| CPU-minutes, all clustering steps | 284.10 | 291.68 | 313.05 |
| largest step, CPU-minutes | 284.10 | 92.72 (`cluster.nr90.12.align`) | 33.03 (`cluster.nr90.12.align.part-00000`) |
| wall time of clustering | 1,287.9 s (16 threads) | 1,978.7 s (20:05:16 to 20:38:15 CST) | 3,792.6 s (20:42:22 to 21:45:35 CST) |
| largest peak RSS of a step | 4,234,657,792 bytes (source not recorded then) | 3,683,426,304 (`cluster.nr50.16.align`) | 2,786,729,984 (`cluster.nr90.12.align.part-00002`) |
| work directory peak on disk | not measured | 5,648,736,256 bytes allocated (5.26 GiB), at 20:38:12, during `cluster` | 6,262,779,904 bytes allocated (5.83 GiB), at 21:32:36, during the parts of `cluster.nr50.20.align` |
| run-step staging peak | | 441,909,248 bytes | 231,763,968 bytes |
| whole unit (`CPUUsageNSec`; `MemoryPeak`, page cache included) | | 293.60 CPU-min; 7,779,356,672 bytes | 317.86 CPU-min; 8,928,849,920 bytes |
| tokens, postings, dim, manifest, catalog | | 1.69 CPU-min in all, peak 2,209,050,624 bytes (`postings`) | 1.65 CPU-min in all, peak 2,165,448,704 bytes (`postings`) |

The disk figures are samples every 15 s, so a peak shorter than that can be missed; the scratch directory's peak
(5,648,744,448 bytes in run A) is the work directory plus a few kilobytes.

## Against the estimates

**The release's one-process estimate** (`stage.cluster_cpu_min`, card B-36), 284.1 CPU-minutes for 2,053,315 rows:
run A's clustering used 291.68, 2.7% over (*derived*: 291.68 / 284.10 − 1); run B used 313.05, 10.2% over
(*derived*), the cost of three parts per call: 32 more processes (48 parts and 16 joins against 16 call steps), each
paying a start-up of about 0.1 to 0.25 CPU-minutes, every NR90-step and NR50-step-0 prefilter part rebuilding its
table, and the shared indexes built by the two sensitive prefilters' cut steps (2.31 and 3.02 CPU-minutes).

**C-22's estimate per call** (the cut records' `parts[].cpu_min`, the model `cluster.call_cpu_min` plus the part
overhead) against each run record's `cpu_min`. Run A, every call one part, 12 threads:

| call | estimate | measured | measured / estimate |
|---|---|---|---|
| `cluster.nr90.01.rescorediagonal` | 0.359 | 0.399 | 1.11 |
| `cluster.nr90.06.rescorediagonal` | 0.250 | 0.344 | **1.38** |
| `cluster.nr90.07.align` | 4.505 | 4.659 | 1.03 |
| `cluster.nr90.11.prefilter` | 8.489 | 8.746 | 1.03 |
| `cluster.nr90.12.align` | 87.545 | 92.715 | 1.06 |
| `cluster.nr90.members` | 5.127 | 5.577 | 1.09 |
| `cluster.nr50.01.rescorediagonal` | 0.171 | 0.298 | **1.74** |
| `cluster.nr50.06.rescorediagonal` | 0.138 | 0.268 | **1.95** |
| `cluster.nr50.07.align` | 6.817 | 6.778 | 0.99 |
| `cluster.nr50.11.prefilter` | 5.147 | 5.363 | 1.04 |
| `cluster.nr50.12.align` | 35.771 | 40.522 | 1.13 |
| `cluster.nr50.15.prefilter` | 10.390 | 11.085 | 1.07 |
| `cluster.nr50.16.align` | 31.543 | 32.521 | 1.03 |
| `cluster.nr50.19.prefilter` | 20.464 | 22.321 | 1.09 |
| `cluster.nr50.20.align` | 35.611 | 36.192 | 1.02 |
| `cluster.nr50.members` | 10.718 | 10.887 | 1.02 |
| all 16 calls | 263.05 | 278.68 | 1.06 (*derived*) |

Run B, every call in three parts at 4 threads each (48 parts; the estimate of each part is the same within a call,
since the cut divides the call's units equally; measured / estimate per part, lowest to highest):

| call | estimate per part | measured parts | measured / estimate |
|---|---|---|---|
| `cluster.nr90.01.rescorediagonal` | 0.120 | 0.243, 0.263, 0.263 | **2.03 to 2.20** |
| `cluster.nr90.06.rescorediagonal` | 0.083 | 0.139, 0.135, 0.141 | **1.62 to 1.69** |
| `cluster.nr90.07.align` | 1.502 | 1.974, 1.696, 1.495 | 1.00 to **1.31** |
| `cluster.nr90.11.prefilter` | 6.126 | 5.662, 5.004, 5.101 | 0.82 to 0.92 |
| `cluster.nr90.12.align` | 29.182 | 33.031, 31.658, 31.499 | 1.08 to 1.13 |
| `cluster.nr90.members` | 1.709 | 2.751, 1.688, 1.418 | 0.83 to **1.61** |
| `cluster.nr50.01.rescorediagonal` | 0.057 | 0.208, 0.211, 0.210 | **3.65 to 3.69** |
| `cluster.nr50.06.rescorediagonal` | 0.046 | 0.119, 0.123, 0.123 | **2.60 to 2.69** |
| `cluster.nr50.07.align` | 2.272 | 2.428, 2.380, 2.385 | 1.05 to 1.07 |
| `cluster.nr50.11.prefilter` | 3.875 | 2.857, 2.905, 3.069 | **0.74 to 0.79** |
| `cluster.nr50.12.align` | 11.924 | 14.103, 13.510, 14.260 | 1.13 to 1.20 |
| `cluster.nr50.15.prefilter` (shared index) | 2.606 | 2.974, 2.660, 2.648 | 1.02 to 1.14 |
| `cluster.nr50.16.align` | 10.514 | 11.192, 10.945, 11.783 | 1.04 to 1.12 |
| `cluster.nr50.19.prefilter` (shared index) | 5.983 | 6.592, 6.166, 6.451 | 1.03 to 1.10 |
| `cluster.nr50.20.align` | 11.870 | 12.556, 12.111, 13.013 | 1.02 to 1.10 |
| `cluster.nr50.members` | 3.573 | 3.649, 4.084, 3.382 | 0.95 to 1.14 |
| all 48 parts | 274.32 | 293.26 | 1.07 (*derived*) |

The shared index (card B-43) is built in the cut step of each sensitive prefilter when it has more than one part: its
estimate (`index_cpu_min`) was 2.071 and 2.025 CPU-minutes, and the cut steps used 2.307 (`cluster.nr50.15`, +11%) and
**3.021** (`cluster.nr50.19`, +49%), the latter including the cut's own start-up of about 0.1.

Off by more than 20% (AGENTS §6), in bold above:

- Every `rescorediagonal` call and part (each under 0.4 CPU-minutes): a step's fixed start-up, about 0.1 to 0.2
  CPU-minutes (each `.cut` step, which does almost no work, used 0.09 to 0.20), is as large as the work. C-22's own
  sample shows the same (`rC.1b`: 0.148 estimated, 0.357 measured for `cluster.nr50.01.rescorediagonal`,
  `tests/r1/c22_cluster_calls_measured.tsv`); C-22 claims 10% only for calls of 2 CPU-minutes or more.
- In run B, parts of 1.4 to 3.9 CPU-minutes: `cluster.nr90.07.align` part 0 (+31%), `cluster.nr90.members` part 0
  (+61%; its other parts −1% and −17%: an equal share of entry bytes is not an equal share of work at this size),
  and the three parts of `cluster.nr50.11.prefilter` (−21% to −26%; the same call in one part, run A, was +4%: a
  part's own table is estimated at the conservative rebuild rate of card B-41, 2.1e-4 CPU-s per target, which card
  B-43 measured as 1.5e-4 to 1.8e-4 at larger targets).
- The shared index of `cluster.nr50.19.prefilter`, +49% (above).

Every call or part of 4 CPU-minutes or more is within 20% of its estimate (largest gap `cluster.nr50.12.align` part
2, +20%; *derived* 14.260 / 11.924 = 1.196), and in both runs the per-query work as a whole is 6% to 7% over the
estimate, so at this size the cut errs slightly low. Nothing here is near the part budget of 240 CPU-minutes: the
largest process of either run is 92.72.

## Earlier evidence

The 50 Gbp build's compatibility check (`planning/measurements/2026-11-11_build_r50.md`, A1) built a release `rC.1b`
over the same stages with the tree of 2026-11-11 (before cards B-43, C-22 and B-45), also one part per call, and found
its clusters equal to `rN.1b`'s row for row (by count of rows on each side only). This card repeats that with the
current code, compares the bytes, and adds run B: real multi-part cuts, with the parts of the NR50 steps 1 and 2
searching the shared prefilter index of card B-43, which the one-part runs never use.

## Files on the pod

`/tmp/if2453/b42/`: `bin/` (the scripts named above, and `records.py`, which prints the per-step and per-call tables
from the records), `logs/b42-A.log`, `logs/b42-B.log`, `logs/cmp-*.log` (the comparison JSON), `logs/du-A.tsv`,
`logs/du-B.tsv`, and the scratch roots `view-A`, `view-B` (their own files: the release's tables, run records and
`cluster_plan` records). They can be removed.
