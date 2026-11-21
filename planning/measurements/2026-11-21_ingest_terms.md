# The ingest estimate prices gzip input and bases read but not stored (card B-61)

Question: `stage.ingest_cpu_min` charges every byte of FASTA at the stored rate (4.11 CPU-minutes per 10^9 bytes, card
B-52), less the bases an earlier stage holds. Card B-58 measured three things it does not know: gzip input, bases read
and not stored, and a dearer stored rate on many small contigs. For `SPX1` (304 `.fa.gz`, 256.9 GB of FASTA, built
with `--min-contig-len 5000`, 85.23% of its bases under it) the code says about 1,056 CPU-minutes where card B-58's
arms derive about 264.

Labels (AGENTS §6): *measured* (host, date, command or file), *derived* (formula and inputs), *estimate*.

## 1. Design: how the estimate learns, before ingest, what will be stored

What is known before the ingest step starts (*measured* by reading the code, `src/biobox/stage.py`,
`src/biobox/l0seq/store.py`, and `/tmp/if2453/s2/inputs/SPX1`, 2026-11-21):

- `--objects` gives each object's path; `object_sizes` already streams every `.gz` object through zlib to count its
  bytes of FASTA (`fasta_bytes`), in `build-stage`'s own process, for the refusal and for the batch decision. For
  `SPX1` that is a full decompression of 75.3 GB of gzip into 256.9 GB of FASTA each time `build-stage` plans the
  ingest.
- The `inputs` step publishes `genome_input` before ingest: the `contig_uid` of every contig the stage keeps (ingest
  refuses a stage whose `genome_input` lists a contig under the least length or misses one over it,
  `store._check_genome_input`). It holds no lengths, so the kept bases are not in the root before ingest.
- Ingest's own pass 1 counts every object's records, bases and bases under the least length; in batches
  (`ingest.pass1-<k>`) these counts are published in the pass-1 records, but only after the pass-1 processes ran,
  and those processes are the ones the estimate must size.
- The staging of `SPX1` wrote `samples.tsv` (per sample: bases at or over 5 kb and under it) and, in scratch,
  `/tmp/if2453/c48e1/work/contigs_ge5000.parquet` (the kept contigs and their lengths). Neither is part of what
  `build-stage` is given (`--objects`, `--genomes`): reading them would make the C-48e staging's own files an input
  of the code, a new column or file of the inputs.

The three options of the card, judged:

1. **A table of record lengths the staging already wrote**: `samples.tsv` has exactly the figure, but only for
   `SPX1` and in a form the code does not read; making it an input is a change to the inputs' form (`make_inputs`,
   `c48e_inputs.py`, the objects' fields), outside this card's files. Not taken.
2. **A first pass that already runs**: two exist. (a) Ingest's pass 1 knows the dropped bases exactly, but only
   after it has run: it cannot size its own shards, nor the one-process decision. (b) `object_sizes` already reads
   every `.gz` object whole, before any step, in `build-stage`. Counting the records' lengths in that same stream
   (`stage.fasta_counts`: the bytes of FASTA of the records under the least length) needs no new pass over a gzip
   stage, holds nothing but counters (AGENTS §4.1), and is known before the first step. For a stage of plain FASTA
   with a least length it is one read of each file (none built so far). **Taken.** Its cost in `build-stage`'s
   process is measured in section 5.
3. **A figure the stage's conf states**: a number typed into a build's arguments with a note is a figure no test
   checks against the files; and `build-stage` has no such argument (a change to `cli.py` and the drivers). Not
   taken, since 2(b) gives the same figure from the files themselves.

What the estimate assumes when it cannot know: without `--min-contig-len` nothing is dropped and every byte not held
elsewhere is stored, as before (card B-52; it errs high for repeated records, as for `S1`). The counter counts a line
of only spaces or tabs as bases where `fasta.iter_records` drops it, so a record can only look longer than ingest
finds it: kept, not dropped, and the estimate errs high.

Where the terms enter: `ingest_cpu_min(n_bytes, held_bases, *, gz_bytes=0, below_bytes=0)`, called with the counts by
`_batched` (the one-process decision and its message) and `pass1_bytes` (the first-pass cut: a pass-1 shard is
priced at the stage's average cost per byte, still as though every byte not dropped were stored, since a pass-1
shard stores nothing). With both keywords 0 the function and the cut are card B-52's, byte for byte.

What is not priced: held bases stay at 0 in the sum (card B-52's stored rate was fitted with them at 0, and `PL1`,
`H1`, `P3` to `P6` and `V1` hold some; a held term would move their estimates, which this card may not do). Held
bytes cost little beyond reading: card B-58's arm `held B` (every contig held) is within 1% of the new estimate with
the gzip term alone (section 3), because the bytes of a FASTA that are not bases (headers, newlines) are still
counted at the stored rate.

## 2. The arms

Every arm is a `build-stage --depth stored` in a fresh scratch root under `/tmp/if2453/b61/roots`, a detached
`pod_run.sh` unit of 2 threads at `nice 19` in `bio-box_faec1f4.sif` (digest `sha256:d67300aa...`), with the tree of
integration `f60ed52` as it was before this card's code (`/tmp/if2453/b61/tree`), one arm at a time, interleaved
(rep 1 of every arm, then rep 2, then rep 3), `/tmp/if2453/r10/root/STOP` checked before each, the root removed after
its record was read. The figure is `cpu_min` of the ingest step's run record. Script: `scripts/stress/b58_ingest_arms.py`
(card B-58's, extended: modes `groups-inputs`, `g1-inputs`, `c5-counts`, `fit`, `spx1-plan`; units `b61-inputs`,
`b61-arms`, `b61-arms2`, `b61-g1-inputs`, `b61-c5-counts`, `b61-spx1-plan`). Logs and records:
`/tmp/if2453/b61/logs/`, `/tmp/if2453/b61/results.jsonl`, `/tmp/if2453/b61/report.json` (written as
`inputs/report.json`; the samples' inputs were removed after the arms), `/tmp/if2453/b61/plan.json`
(*measured*, wilkcomp02, 2026-11-21 06:04 to about 11:05 CST).

**The 12 files** (rule in the script's docstring): `SPX1`'s 304 objects in the order of (file size, path), less card
B-58's 4 middle files, are 300; the files at ranks floor((i + 0.5) x 300 / 12), i = 0 to 11, an even spread by rank,
in 6 groups of two neighbours. Groups by size, with what they read (*measured* by the script's own loop over the
decompressed lines, `line_counts`, not the code's counter; the code's counter gave the same bytes under 5 kb for
every one of the 12 files, section 5):

| Group | FASTA bytes | Records | Bytes per record | Bytes in records under 5 kb | Share of bytes under 5 kb |
|---|---:|---:|---:|---:|---:|
| spxg1 | 3,129,093 | 6,386 | 490 | 2,932,981 | 93.7% |
| spxg2 | 134,563,454 | 144,360 | 932 | 105,088,104 | 78.1% |
| spxg3 | 450,359,726 | 711,802 | 633 | 427,508,895 | 94.9% |
| spxg4 | 963,086,613 | 1,208,452 | 797 | 841,176,755 | 87.3% |
| spxg5 | 2,124,973,510 | 2,118,659 | 1,003 | 1,366,543,629 | 64.3% |
| spxg6 | 5,563,406,950 | 6,899,266 | 806 | 4,901,131,525 | 88.1% |
| card B-58's `spx` (4 files) | 1,297,869,659 | | | 1,023,387,001 | 78.9% |
| `SPX1` (304 files) | 256,948,810,775 | 325,805,489 | 789 | 220,918,826,783 | 86.0% |

Kept shares by group run from 5.1% to 35.7% of the bytes, around `SPX1`'s 14.0% (card B-58's 4 files: 21.1%).

**The G1 sample** (at the coordinator's request): every 10th object of `G1`'s `objects.json` (714 GlobDB MAGs as
`.fa.gz`), 996,211,156 bytes of FASTA, 121,666 contigs (8,188 bytes each), one genome per object, every contig kept.

**C5's samples counted again** (card B-58's note gave bytes only for the first): `objects[::10]` 944,466,060 bytes of
FASTA, 932,761,780 bases, 510 records; `objects[::3]` 3,159,336,047 bytes, 1,706 records.

The arms (CPU-minutes of the ingest record; `B` is the same objects again as stage B of the same root, every contig
held; `m` is `--min-contig-len 5000`; the estimate is the new code's, with B-52's two constants and the two below):

| Arm | Runs | Mean | Spread (min to max) | New estimate | Mean / estimate | Card B-52's estimate | Mean / it |
|---|---|---:|---|---:|---:|---:|---:|
| gz (B-58) | 6 | 4.3904 | 4.2027 to 4.5305 | 4.1973 | 1.046 | 3.9517 | 1.111 |
| plain (B-58) | 3 | 4.1443 | 3.9943 to 4.2284 | 3.9517 | 1.049 | 3.9517 | 1.049 |
| few (B-58) | 3 | 4.0746 | 3.8344 to 4.2043 | 3.9517 | 1.031 | 3.9517 | 1.031 |
| view (B-58) | 3 | 4.1124 | 3.9525 to 4.3381 | 4.1973 | 0.980 | 3.9517 | 1.041 |
| held A (B-58) | 3 | 4.2385 | 4.0956 to 4.4781 | 4.1973 | 1.010 | 3.9517 | 1.073 |
| held B (B-58) | 3 | 0.3620 | 0.3333 to 0.3803 | 0.3637 | 0.995 | 0.1181 | 3.065 |
| scale (B-58) | 3 | 14.6723 | 14.5195 to 14.8125 | 13.8763 | 1.057 | 13.0549 | 1.124 |
| spx (B-58) | 3 | 6.4752 | 6.0604 to 6.7426 | 5.7417 | 1.128 | 5.4042 | 1.198 |
| spx5000 (B-58) | 3 | 1.7196 | 1.6454 to 1.7576 | 1.6788 | 1.024 | 5.4042 | 0.318 |
| spxg1 | 3 | 0.0677 | 0.0624 to 0.0735 | 0.0837 | 0.809 | 0.0829 | 0.817 |
| spxg1 B | 3 | 0.0532 | 0.0491 to 0.0574 | 0.0719 | 0.740 | 0.0710 | 0.749 |
| spxg1m | 3 | 0.0543 | 0.0520 to 0.0578 | 0.0720 | 0.755 | 0.0829 | 0.655 |
| spxg2 | 3 | 0.6999 | 0.6627 to 0.7428 | 0.6580 | 1.064 | 0.6231 | 1.123 |
| spxg2 B | 3 | 0.1196 | 0.1028 to 0.1323 | 0.1294 | 0.924 | 0.0944 | 1.266 |
| spxg2m | 3 | 0.2283 | 0.2154 to 0.2440 | 0.2408 | 0.948 | 0.6231 | 0.366 |
| spxg3 | 3 | 2.2747 | 2.1650 to 2.4663 | 2.0381 | 1.116 | 1.9210 | 1.184 |
| spxg3 B | 3 | 0.3603 | 0.3156 to 0.4177 | 0.3108 | 1.159 | 0.1937 | 1.860 |
| spxg3m | 3 | 0.3481 | 0.3326 to 0.3636 | 0.3409 | 1.021 | 1.9210 | 0.181 |
| spxg4 | 3 | 4.7097 | 4.5318 to 5.0553 | 4.2787 | 1.101 | 4.0283 | 1.169 |
| spxg4 B | 3 | 0.5497 | 0.5361 to 0.5645 | 0.5289 | 1.039 | 0.2785 | 1.974 |
| spxg4m | 6 | 0.9078 | 0.8757 to 0.9833 | 0.9392 | 0.967 | 4.0283 | 0.225 |
| spxg5 | 3 | 10.1675 | 9.5035 to 10.6522 | 9.3561 | 1.087 | 8.8036 | 1.155 |
| spxg5 B | 3 | 1.0505 | 1.0221 to 1.1025 | 0.9919 | 1.059 | 0.4394 | 2.391 |
| spxg5m | 3 | 3.7635 | 3.6765 to 3.8459 | 3.9310 | 0.957 | 8.8036 | 0.427 |
| spxg6 | 3 | 27.6782 | 26.5625 to 29.8981 | 24.3821 | 1.135 | 22.9356 | 1.207 |
| spxg6 B | 3 | 3.3747 | 3.1879 to 3.5407 | 2.7310 | 1.236 | 1.2845 | 2.627 |
| spxg6m | 3 | 4.8609 | 4.7145 to 4.9451 | 4.9246 | 0.987 | 22.9356 | 0.212 |
| g1 | 3 | 4.2400 | 4.1936 to 4.2737 | 4.4234 | 0.959 | 4.1644 | 1.018 |

(Estimates and ratios *derived* from the means with `stage.ingest_cpu_min` as committed; held bases of a `B` arm are
the bases of every contig, bytes as counted above.)

**Load** (the 1-minute load average at each arm's start is in `results.jsonl`): rep 1 of groups 1 to 4 ran at 47 to
55 (the release train and two rehearsals beside it), every later run at 17 to 25 (the train had stopped at 06:39 and
the integration suite ran in 8 shards). The same arm cost 7% to 13% more at about 50 than at about 20: spxg2 0.743
against 0.678, spxg3 2.466 against 2.179, spxg4 5.055 against 4.537, spxg2m 0.244 against 0.220, spxg3m 0.364
against 0.340 (*derived*, rep 1 against the mean of reps 2 and 3). Card B-58 saw 6.7% between about 30 and 55. The
means mix both loads as card B-58's did; a build at night at a load near 55 runs some 5% over a mean of this kind.

## 3. The fit, term by term

Model (bytes of FASTA throughout; held bases as card B-52 counts them):

    estimate = 0.07 + 4.11 x stored + D x below + G x gz      (per 10^9 bytes; stored = read - held - below)

The start-up (0.07) and the stored rate (4.11) are card B-52's and are not refitted: every stage built from plain
FASTA without a least length keeps its estimate to the byte only so. G and D were fitted together by least squares of
the relative error (`b58_ingest_arms.py fit`, *derived*) over the 14 arms whose contigs are long or which drop the
short ones and which are not start-up bound: gz, plain, few, view, held A, held B, scale, spx5000, spxg2m to spxg6m,
g1. Two constants from 14 means: 12 degrees of freedom left.

| Constant | Mechanism | Fit | With one arm left out | Largest error of the arm left out | Sources | In the code |
|---|---|---:|---|---:|---|---:|
| G, gzip | zlib inflates every byte of a `.gz` object, in both reading passes of `store.ingest` | 0.258 | 0.252 to 0.259 | 5.8% (spxg3m) | C5's genomes (gz against plain, held B), SPX1's unbinned contigs (the `m` arms), G1's MAGs | 0.26 |
| D, under the least length | a record under `--min-contig-len` is parsed in both passes and hashed once (`contig_uid`), never written | 0.131 | 0.113 to 0.138 | same | SPX1, 6 groups (B-58's 4 files and 5 of card B-61's groups) | 0.14 |

The arms alone, each with the other constant at 0.26: D from spx5000 0.179, spxg2m 0.020, spxg3m 0.156, spxg4m 0.102,
spxg5m 0.017, spxg6m 0.126 (*derived*). The spread is wide because D is a small part of each arm: in spxg5m the
dropped bytes are 64% of what it reads and cost 0.19 of its 3.76 CPU-minutes, so 5% of noise in the arm moves D from 0
to 0.3. The term is small for the same reason: at 0.14, `SPX1`'s 220.9 GB under 5 kb cost 30.9 CPU-minutes of 245.9
(section 4). The data leave a choice between 0.113 and 0.138; the code takes 0.14 (the top, rounded up: an estimate
too low stops a build at night, DEC-53).

G is not fitted on C5. C5's own gz-against-plain arms give 0.261 unpaired and 0.164 paired by run (card B-58); at
0.164 the arm `held B` (which reads and inflates every byte and stores nothing) would be 1.33 times its estimate, and
the fit's 0.258 puts it at 1.00. Card B-58's reading tests predicted 0.17 per GB for SPX1's files and 0.44 for C5's
(two passes): 0.26 lies between.

**The third term (the stored rate of many short contigs) is not written.** The data say it differs: every arm that
stores SPX's short contigs runs over the estimate, spxg2 1.064, spxg5 1.087, spxg4 1.101, spxg3 1.116, spx 1.128,
spxg6 1.135, by about 18 to 29 microseconds per stored record (*derived*: the arm's excess over its estimate divided
by its records, spxg2 17.5, spxg3 19.9, spxg4 21.4, spxg5 23.0, spxg6 28.7); the arms that hold them (`B`) by about 4
to 6 microseconds per record (spxg3 B 4.2, spxg6 B 5.6). A cost per record is a mechanism (Python's work per record in
both passes, the index and occurrence rows), from two sources (B-58's `spx`, card B-61's groups). But the code does not
know the number of records before ingest: it counts records' lengths only with a least length (section 1), and then
the records stored are the long ones; a term on the contigs `genome_input` lists would move every plain stage's
estimate (S1 lists 27,258 contigs; at 25 microseconds, 0.011 CPU-minutes), which this card may not do. Every
all-kept short-contig arm is within 15% anyway (largest 1.135). No stage to be built stores contigs this short: SPIRE
stages are built with `--min-contig-len 5000` (DEC-54, DEC-60), and G1's contigs (8.2 kb on average) run at 0.960. If
a stage of unbinned contigs is ever built without a least length, expect about 1.1 to 1.14 of the estimate at a load
of 20, more at 55: a card then.

**Arms outside 15%, named** (the test holds them separately):

- spxg1, spxg1 B, spxg1m (3 MB): the start-up of a process; the records are 0.049 to 0.074 CPU-minutes against
  estimates of 0.072 to 0.084, 0.016 to 0.019 CPU-minutes over (DEC-43: a step that small has no rate). The estimate
  errs high.
- spxg3 B (1.159) and spxg6 B (1.236): 0.7 and 6.9 million short contigs read again and every one held, the per-record
  cost above. No stage holds millions of short contigs of an earlier stage (`SPX1` holds none of the root's; card
  C-48e).

**Before and after, the built stages** (*derived*, `ingest_cpu_min`; records from `runs/stage=<S>/ingest.json`):

| Stage | Bytes of FASTA | Held bases | Before (B-52) | After (B-61) | Record | Record / after |
|---|---:|---:|---:|---:|---:|---:|
| S1 | 2,121,124,779 | 0 | 8.787822841690001 | the same | 7.540384 | 0.858 |
| P1 | 1,029,384,352 | 0 | 4.300769686720001 | the same | 4.611177 | 1.072 |
| P2 | 9,716,599,273 | 0 | 40.00522301203001 | the same | 38.33493 | 0.958 |
| P3 | 10,074,413,081 | 6,963 | 41.475809144980005 | the same | 43.889413 | 1.058 |
| P4 | 10,101,497,437 | 7,053 | 41.587125478240004 | the same | 42.979174 | 1.034 |
| P5 | 10,092,827,773 | 185,683 | 41.550758989900004 | the same | 43.854939 | 1.055 |
| P6 | 10,087,447,605 | 17,198 | 41.529338972770006 | the same | 41.485983 | 0.999 |
| H1 | 1,962,116,403 | 141,683,767 | 7.551978133960001 | the same | 7.644566 | 1.012 |
| PL1 | 7,425,668,755 | 1,177,838,866 | 25.748580843790002 | the same | 26.61184 | 1.034 |
| V1 | 5,228,187,536 | 43,795,699 | 21.377850450070003 | the same | 21.816238 | 1.021 |
| C5 | 9,509,853,070 (all `.gz`) | 267,273,180 | 38.0570033479 | 40.5296 | 45.803255 | **1.130** |

`C5` moves by the gzip term alone (0.26 x 9.510 = 2.473 CPU-minutes): 1.2035 before, 1.130 after. That is under the
1.15 the card expected: the gzip rate comes from the arms (0.258 fitted without C5's record), not from C5; with the
paired 0.164 it would be 1.156 and `held B` would miss its estimate by 33%. C5's record stays the sum of three small
effects (card B-58): gzip, C5-type genomes at 4.39 per GB in plain copies, and load.

**Plans of built stages**: none changes. No built stage was ingested in batches (no `l0seq/ingest_plan` in the root,
*measured*, 2026-11-21), so no `pass1_bytes` was ever recorded; a published pass-1 record is read back by
`pass1_bytes` before any estimate (unchanged); none of the eleven used `--min-contig-len` (no run record names it,
*measured*); only `C5` read `.gz` objects (the root's `l0_manifest`: 2,032 of 2,032 `.gz` for C5, 0 for the other
ten, *measured*), and C5's ingest is published, so `build-stage` never estimates it again (`_refusal` and
`_batched` run only while the ingest record is missing). The estimate is not written into any published file.

## 4. `SPX1` and `G1` by the new code

**`SPX1`** (`b58_ingest_arms.py spx1-plan`, unit `b61-spx1-plan`, read-only on `/tmp/if2453/s2/inputs/SPX1`, the tree
with this card's code, *measured* counts, *derived* estimates; `/tmp/if2453/b61/plan.json`): 304 objects, 75,267,693,104
bytes on disk, 256,948,810,775 bytes of FASTA (all `.gz`), 220,918,826,783 of them in records under 5,000 bp; held 0.

- One-process estimate: 0.07 + 4.11 x 36.030 (stored) + 0.14 x 220.919 (under 5 kb) + 0.26 x 256.949 (gzip) = 0.07 +
  148.08 + 30.93 + 66.81 = **245.9 CPU-minutes**, against card B-58's 264 (6.9% under it) and card B-52's 1,056.1.
  It is over `SHARD_MARGIN` of 300 (240), so ingest runs in batches, as before.
- First pass: `pass1_bytes` 73,464,662,516 bytes on disk per shard (17,100,339,219 before, card C-48e): **2 shards**,
  296 objects (73.28 GB on disk, 250.15 GB of FASTA, estimated at 239.1 CPU-minutes as though stored) and 8 objects
  (1.99 GB, 6.80 GB of FASTA, 6.8). Before: 5 shards. A pass-1 shard stores nothing, so its own cost is far under
  that estimate (the `B` arms, which read, inflate and hash everything and store nothing, cost 0.6 to 0.9
  CPU-minutes per GB of SPX FASTA, *derived*: about 150 to 225 for the large shard).

**`G1`** (planning/measurements/2026-11-18_globdb_inputs/G1/est.json: 7,139 `.fa.gz`, 10,087,932,724 bytes of FASTA,
5,699,697 bases held when staged, no least length): card B-52's estimate 41.508, the new one 41.508 + 0.26 x 10.088 =
**44.131**. Carried from its own sample (the g1 arm, 4.186 CPU-minutes per GB after start-up, at a load of 20 to 23),
the stage would cost about 42.27 (*derived*): 1.018 of the old estimate and 0.958 of the new. `G1` behaves like C5's
genomes, not like SPX's short contigs (8.2 kb per contig against 0.8 kb), and would not have stopped under either
estimate at that load; at a load near 55 add 7% to 13% (1.09 to 1.15 of the old, 1.03 to 1.08 of the new). If
`SPX1` is built first, `G1` holds 12,052 more contigs (196.6 Mbp, card C-48e) and its estimate falls by 0.81.

## 5. What counting costs, and the counter checked

On card B-61's 12 files (`spx1-plan`, *measured*, getrusage in one process): `fasta_bytes` (decompression only) 39.5
CPU-seconds, `fasta_counts` 57.8, for 9.24 GB of FASTA: counting adds 46% to the read `object_sizes` already makes
(*derived*). On the whole of `SPX1`, `object_sizes` with the counts took 1,702.5 CPU-seconds (28.4 CPU-minutes) in
`build-stage`'s own process, where decompression alone would be about 19.4 (*derived*: 28.4 / 1.46). The counter's
bytes under 5 kb equal the script's independent line loop on every one of the 12 files (e.g. spxg1: 99,519 +
2,833,462 = 2,932,981; spxg2: 37,353,407 + 67,734,697 = 105,088,104), and the tests hold it against a loop over the
lines and `fasta.iter_records` at chunk sizes of 1 to 2^20 bytes.

## 6. Open points

- The per-record cost of storing short contigs (18 to 29 microseconds per record) is not in the estimate (section 3).
  It matters only for a stage storing millions of contigs under about 1 kb.
- Load moves an arm by 7% to 13% between about 20 and 50; the constants are means over both. The guard's 20% leaves
  room for it at the fitted stages, not for a further 15% error.
- `build-stage` spends about 28 CPU-minutes reading `SPX1`'s objects before the first step (19 of them before this
  card), outside any run-step, once per invocation.
- The batch warning (`ingest_cpu_min(x["bases_read"])` per batch of the plan) still prices every base read as stored
  (it errs high and never fires for `SPX1`: its largest object holds 15.6 Gbp, 64 CPU-minutes so priced).
