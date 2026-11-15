# Gene calling in batch steps under the CPU budget (card B-48)

Date: 2026-11-15, wilkcomp02 (the pod's clock, CST). Every number carries a label: *measured* (host, date, command or
file), *derived* (formula and inputs) or *estimate*. Read-only on the data root `/tmp/if2453/r10/root` (bound
read-only into every command that read it); nothing under `/stor` or `/tmp/if2453/s2/inputs` was written. Stage H1's
half-built state (card C-24) and its units were not touched; only its published `genome_input` and `contig_index` and
its gene-call run record were read.

Small results are committed beside this file in [`2026-11-15_genecall_batches/`](2026-11-15_genecall_batches/).

## Summary

- **Design** (`src/biobox/l1/genecall.py`, module docstring): a batch is a run of whole genomes, consecutive in
  `genome_uid` order, called in its own process exactly as the one-process step calls them, writing its unsorted
  scratch rows, its skips and a record into the plan's work directory in scratch (never the root). The merge is the
  run-step `genecall`: `_publish`, the one-process step's own writer, once over every batch's rows, so it writes the
  same files, ledger and result. Nothing in the output can depend on the number of batches: training is per genome
  (single mode trains on one genome's contigs, meta modes train nothing), every table is sorted on a key unique per
  gene, file cuts and the protein bucket width are computed from the sorted stream and the whole stage's distinct
  proteins in the merge, the ledger's first-of-protein window runs over every gene in the merge, and a contig's
  owner (DEC-26, DEC-27) is read from the whole stage's inputs by every batch. Card W-35's parts are the model;
  the batches are not `run-step` steps (below, deviations).
- **Byte identity on the fixture**: at 1, 2, 3 and 5 batches, 1 and 3 batches at once, 1 to 3 workers each, through
  build-stage at 1 and 3 jobs, after `STOP` between batches and after SIGTERM and SIGKILL of a batch: every gene-call
  file, the ledger and the printed result equal the one-process call's (tests); a release built on top (IGRs, tapes,
  clusters, tokens, postings, `contig_dim`) is the same file for file.
- **Rates** (*measured*, one thread, the calls alone, every 30th genome): single 0.7941 s per Mbp read (P6; P2 0.7838),
  gv-meta 4.2019 (S1; H1 3.7795), meta 2.9450 linear and 4.4205 circular (card C-26's PL1 sample). The step is the
  calls times 1.23 (*derived*: the largest of the eight built stages' record over its calls) plus 15 s a process
  (*estimate*). Every built stage's record is at its estimate or under it by at most 15.8% (H1).
- **PL1** on its real inputs (*derived*, `b48_pl1_plan.py`, the code's own `cut_plan`): 68,693 genomes, 7,046 not
  called (every contig owned by an earlier stage), 5,924,552,875 bases called; **872.6 CPU-minutes in one process**
  (the card's figure was about 860, calls only and every genome); **4 batches of 218.3 CPU-minutes each**, all under
  the 240 of the margin, and a merge of 0.8.
- **Byte identity on S1 in scratch**: S1 gene-called again by build-stage in 4 batches at 2 at once (4 workers each),
  in a view of the root, the runtime image of `faec1f4`: all five gene-call files (the ledger among them) have the
  sha256 of the root's own, made on 2026-11-10 by the code and image of then, and the merge's run record lists the
  same outputs (path, bytes, sha256). Each batch cost 30.70 to 31.87 CPU-minutes against 36.09 estimated, the merge
  0.43 (*measured*); the root was only read.

## 1. What the built stages asked of gene calling, and what they cost

`b48_stage_inputs.py` (unit `b48-stats`, 11:19:19 to 11:19:23): for each gene-called stage, its genomes from
`genome_input`, each contig's length from the `contig_index` that holds its default variant, the mode the code
calls each genome in (`genecall.mode_of`, the rule of `_finder`), and the bases the caller reads (a circular contig
with its wrap flanks). Every stage built so far is linear and of one mode (*measured*, `inputs.json`):

| Stage | mode | genomes | of which every contig owned elsewhere | bases read | genes written | record `cpu_min` | threads |
|---|---|---|---|---|---|---|---|
| S1 | gv-meta | 27,258 | 0 | 1,664,109,938 | 2,398,087 | 142.786 | 16 |
| H1 | gv-meta | 7,228 | 514 (141,683,767 bases) | 1,961,515,105 | 2,276,451 | 146.091 | 36 |
| P1 | single | 250 | 0 | 1,016,625,803 | 927,920 | 15.392 | 16 |
| P2 | single | 2,247 | 0 | 9,596,180,328 | 8,722,859 | 143.288 | 16 |
| P3 | single | 2,671 | 0 | 9,949,084,485 | 9,101,793 | 153.677 | 16 |
| P4 | single | 2,329 | 0 | 9,976,464,003 | 9,078,865 | 151.406 | 16 |
| P5 | single | 2,289 | 0 | 9,967,897,315 | 9,144,001 | 155.692 | 16 |
| P6 | single | 2,180 | 0 | 9,962,553,301 | 9,064,930 | 157.196 | 32 |

H1's 514 genomes all of whose contigs S1 owns were called by the code of its build and their rows dropped (DEC-26);
their bases are in its "bases read". Records: `runs/stage=<S>/genecall.json`, 2026-11-10 (S1, P1, P2), 2026-11-12 (P3
to P6), 2026-11-15 (H1).

## 2. The calls alone, by mode and topology

`b48_call_rate.py` (units `b48-rate-<S>`, 11:22:10 to 11:26:57, one thread each): every 30th genome of the stage in
`genome_uid` order, read from the store, then called with the step's own `genecall._call`, the CPU time of the calls
alone (`getrusage` around each call); the PL1 rows are card C-26's `c26_rate.py`, the same method on every 30th new
PLSDB record. The bases read of the PL1 sample (with the wrap flanks) were counted again here by `b48_pl1_lengths.py`
(unit `b48-pl1`), which found the same 2,054 records and bases as `rate.json`.

| Sample | genomes | bases | bases read | genes | CPU-s | s per Mbp read |
|---|---|---|---|---|---|---|
| S1, gv-meta | 908 | 55,446,997 | 55,446,997 | 79,300 | 232.980 | 4.2019 |
| H1, gv-meta | 240 | 65,344,071 | 65,344,071 | 81,208 | 246.970 | 3.7795 |
| P2, single | 74 | 304,636,972 | 304,636,972 | 276,162 | 238.786 | 0.7838 |
| P6, single | 72 | 351,376,483 | 351,376,483 | 316,956 | 279.027 | 0.7941 |
| PL1 circular, meta | 1,874 | 169,713,014 | 295,462,404 | 183,128 | 1306.076 | 4.4204 |
| PL1 linear, meta | 180 | 17,141,185 | 17,141,185 | 17,688 | 50.480 | 2.9449 |

All *measured*, wilkcomp02, 2026-11-15, image `bio-box_faec1f4.sif`. One rate per mode does not fit meta mode's two
topologies even per base read (4.42 against 2.94): the circular records are mostly small plasmids, so the population
and not only the flanks differs. The estimate therefore takes a rate per mode and topology (the card's "bases by
kind, topology and mode"), each the larger sample, rounded up: `CALL_S_PER_MBP`. Circular contigs in single and
gv-meta mode have not been measured (no such stage is built); they take the linear rate times meta mode's circular to
linear ratio, 1.5011 (*derived*: 4.4205 / 2.9449; an *estimate*).

## 3. The step over its calls, and the estimate against every record

Each record over its calls at those rates (*derived*: `cpu_min` × 60 / (rate × bases read / 10^6)): S1 1.2252, H1
1.0635, P1 1.144, P2 1.128, P3 1.167, P4 1.147, P5 1.180, P6 1.192. The rest of the step (reading the sequences,
building the rows, hashing, the scratch files, the sort and the writing, and the workers' contention) is 6% to 23%
of the calls, with no clear dependence on genes or bases (per gene it is 0.12 to 0.17 ms on the P stages against
0.56 to 0.66 on S1 and H1, *derived*). One factor fits every stage: `STEP_OVERHEAD` = 1.23, the largest, rounded up.
With `PROCESS_S` = 15 s a process (*estimate*), `genecall.one_process_cpu_min` gives (*derived*; the test
`tests/r1/test_b48_genecall_estimate.py` holds each within [record, 1.2 × record]):

| Stage | estimate | record | estimate over record |
|---|---|---|---|
| S1 | 143.59 | 142.79 | +0.6% |
| H1 | 169.21 | 146.09 | +15.8% |
| P1 | 16.80 | 15.39 | +9.2% |
| P2 | 156.47 | 143.29 | +9.2% |
| P3 | 162.21 | 153.68 | +5.6% |
| P4 | 162.66 | 151.41 | +7.4% |
| P5 | 162.52 | 155.69 | +4.4% |
| P6 | 162.43 | 157.20 | +3.3% |

and the calls-only samples: S1 232.98 (+0.0%), H1 274.57 (+11.2%), P2 241.91 (+1.3%), P6 279.03 (+0.0%), PL1
circular 1306.09 (+0.0%), PL1 linear 50.48 (+0.0%), PL1 as a whole 1356.57 against 1356.556 (*derived*).

## 4. Stage PL1: the estimate and the plan on its real inputs

`b48_pl1_lengths.py` (unit `b48-pl1`, 11:19:19 to 11:20:19): every record of the PLSDB FASTA by `ids.contig_uid`,
68,693 contigs, their lengths into `lengths.parquet` (PL1 is not ingested, so no `contig_index` holds them).
`b48_pl1_plan.py` (unit `b48-pl1plan`, 11:48:50 to 11:48:58): the genomes of `genomes.json` in `genome_uid` order, a
contig owned by an earlier stage when some `contig_index` of the root holds its default variant (DEC-26), each genome's
estimate by `genecall.genome_cost_us` (0 when not called), cut by `genecall.cut_plan` at 80% of 300 CPU-minutes, as
build-stage cuts it when the step starts (*derived*, `plan.json`):

| | value |
|---|---|
| genomes; of which not called (every contig owned by S1 or P1 to P6) | 68,693; 7,046 (C-26 counted 7,046 such contigs, *measured* there) |
| bases of all genomes; called | 7,102,391,741; 5,924,552,875 (C-26's "bases the stage adds", to the base) |
| one process | **872.6 CPU-minutes**: over 240, so batches |
| the largest genome | 1.08 CPU-minutes (genome `5ea2e633…`) |
| batch 0 | 17,367 genomes, 1,471,826,928 bases, **218.3** CPU-minutes |
| batch 1 | 16,928 genomes, 1,493,140,979 bases, **218.4** |
| batch 2 | 16,966 genomes, 1,484,779,065 bases, **218.3** |
| batch 3 | 17,432 genomes, 1,474,805,903 bases, **218.3** |
| the merge | **0.8** CPU-minutes (6,393,185 genes expected at 1,079.1 per Mbp; 15 s and 5.5 µs a gene, section 5) |
| refusal | none |

(The first run of this plan, 11:48, gave the merge 10.9 CPU-minutes at a provisional 0.1 ms a gene; the constant was
set from S1's merge, section 5, and the plan made again at 12:08:43 with the final code: the same batches, to the
genome.)

The card's 860 was the calls of every genome at the sample's 7.26 s per Mbp (*derived* in card C-26); this estimate
leaves out the 7,046 genomes not called (about 143 of the 860, *derived* there) and adds the step's overhead (1.23):
717 × 1.23 + 0.25 ≈ 882 (*derived*), here 872.6 from each genome's own lengths and topology. Within 20% of the card's.
At 4 batches × 4 workers the batches would take about an hour of wall time (*estimate*: 218 / 4 × 1.1).

## 5. Byte identity on stage S1 in scratch

The view (`b48_view.py`, `/tmp/if2453/b48/s1/view`): real directories holding links to S1's `genome_input`,
`contig_index`, `contig_occurrence`, `seq_locator`, `seq_shards`, `l0_manifest`, ingest ledger and `store_params.json`
in the root, and copies of S1's `inputs` and `ingest` run records, so build-stage passes over those steps and reaches
gene calling with the stage looking ingested and not gene-called. No gene-call directory was linked; they and the
gene-call ledger were written into the view's own directories, the root bound read-only.

The launch (`b48_driver.sh batches`, unit `b48-s1-batches`, 11:49:22 to 12:07:05, `pod_run.sh --digest
sha256:d67300aa…f098a`, `PYTHONPATH` the branch's `src`): `build-stage --stage S1 --release rB48.s1 --stages S1
--threads 8 --scan-jobs 2 --scan-cpus 4` with a scratch profile whose `cpu_min_per_process` is 60, so that S1's
one-process estimate (143.6) is over its 80% (48) and build-stage cuts batches; every other key the pod's. A watcher
wrote `STOP` in the view once the gene-call record appeared, and the launch ended there (`STOP exists; stopping before
S1.features`, exit 75), as card C-26's stage-only launch does.

| Process | genomes | bases read | genes | CPU-minutes | estimate | wall s | peak RSS |
|---|---|---|---|---|---|---|---|
| batch 0 | 6,818 | 416,023,553 | 598,871 | 30.86 | 36.09 | 464.0 | 419.2 MiB |
| batch 1 | 6,912 | 416,029,336 | 599,041 | 30.70 | 36.09 | 462.2 | 397.8 MiB |
| batch 2 | 6,710 | 416,075,372 | 601,366 | 31.87 | 36.09 | 479.5 | 405.2 MiB |
| batch 3 | 6,818 | 415,981,677 | 598,809 | 30.78 | 36.08 | 463.1 | 413.9 MiB |
| merge (run-step `genecall`) | 27,258 | | 2,398,087 | 0.432 | 4.25 (provisional), 0.470 (final) | 22.7 | 1,249.9 MiB |

All *measured* (each batch's `done.json`, the merge's run record). Batches 0 and 1 ran at once, then 2 and 3; batch 3
was held a minute because a test suite of ours (ours elsewhere, 1 CPU) brought the allowance to 7 (C-28's rule). The
unit's totals: 124.85 CPU-minutes (`CPUUsageNSec` 7,490,743,709,000), peak 2.45 GB (*measured*, `pod_run` usage line);
the one-process step of 2026-11-10 recorded 142.79 at 16 workers. Each batch's estimate is 13% to 18% above its
record (within 20%; the test holds them). The merge's estimate was a provisional constant (0.1 ms a gene) set before
it was measured; it was 10 times too high, and `PUBLISH_S_PER_GENE` is now derived from this record ((25.95 s − 15
s) / 2,398,087 genes = 4.57 µs, raised by a fifth: 5.5 µs).

The comparison (`b48_compare.py`, unit `b48-compare-root`, 12:07:34 to 12:07:42; `compare-root-batches.json`):

| File | sha256 (first 16) | root = batches |
|---|---|---|
| `l0seq/protein_seq/…/stage=S1/protein_prefix=all/part-0.parquet` | `2e6172311575a8b3` | equal |
| `l1/feature/…/stage=S1/part-0.parquet` | `d97de7210f0e1c9b` | equal |
| `l1/feature_protein/…/stage=S1/part-0.parquet` | `0da65b4a24609243` | equal |
| `l1/protein_occurrence/…/stage=S1/part-0.parquet` | `78fc23427aeddd1c` | equal |
| `ledger/stage=S1/step=genecall/part-0.parquet` | `29423c6b21348a88` | equal |

and the two run records' `outputs` (path, bytes, sha256) are equal. The root's files were made by the image
`bio-box-test_74b3e05.sif` (its record's digest `sha256:c30be585…`) and the code of 2026-11-10, in one process at 16
workers; these by the image of `faec1f4` and this branch, in batches. The bytes agree, so the stronger check holds and
the one-process call in the current image (`b48_driver.sh one`) was not needed and not run. The result printed by the
merge: 2,398,087 genes, 1,134,550 distinct proteins, modes `["gv-meta"]`, no wrapping gene, no skipped contig.

## 6. Memory

`tests/r2/test_b48_genecall_memory.py` (the form of card W-16's test; *measured*, wilkcomp02, 2026-11-15, each
process through `runstep.measured`, `ru_maxrss` of the largest process of its tree): the synthetic store of
`synth_scale.py` at 10 copies (30 genomes) gene-called as one batch, and at 100 copies (300 genomes) as 10 batches,
2 workers each:

| | 10 copies | 100 copies |
|---|---|---|
| a batch | 325.0 MiB | 323.8 to 325.4 MiB (each of 10) |
| the merge | 325.0 MiB (11,888 genes) | 434.8 MiB (119,452 genes) |

A batch's peak does not grow with the stage. The merge is the one-process step's own writer: DuckDB's sort grows with
the rows up to `DUCK_MEMORY_LIMIT` (2 GB) and spills past it, as in the one-process step. The test's first run held the
merge to the same 32 MiB growth bound as a batch and failed on it (434.8 against 325.0 MiB): that bound was mine and
wrong for a sort under a cap; the test now holds every process under 2,560 MiB (DuckDB's limit and 512 MiB) and the
merge's growth to twice the measured 1.07 KB a gene. On S1 the merge peaked at 1.25 GB for 2.4 M genes and a batch at
398 to 419 MiB (section 5).

## 7. Pod use and scratch

| Unit | CST | What | Threads |
|---|---|---|---|
| `b48-stats` | 11:19:19 to 11:19:23 | `b48_stage_inputs.py` | 2 |
| `b48-pl1` | 11:19:19 to 11:20:19 | `b48_pl1_lengths.py` | 1 |
| `b48-rate-S1`, `-H1`, `-P2`, `-P6` | 11:22:10 to 11:26:57 | `b48_call_rate.py`, one thread each | 4 |
| `b48-pl1plan`, `b48-pl1plan2` | 11:48:50 to 11:48:58, 12:08:3x to 12:08:43 | `b48_pl1_plan.py` | 2 |
| `b48-driver-batches`, `b48-s1-batches` | 11:49:22 to 12:07:05 | S1 in batches (section 5) | 8 |
| `b48-compare-root` | 12:07:34 to 12:07:42 | `b48_compare.py` | 2 |

Beside them, `pod_pytest.sh` runs (32 GB scopes, 2 threads and their subprocesses). At most about 12 threads of ours at
once; nothing else of ours computed (`top` at 11:17 and 11:49: load 0.27 and 1.51, only our own processes). Scratch:
`/tmp/if2453/b48` was 330 MB at 12:01 with the batches' rows (293 MB) and 483 MB at 12:10 with the published files in
the view (`du -sh`, *measured*); its peak, while the merge's files were in staging beside the batches' rows, was under
1.3 GB (*derived*: 330 MB + twice S1's 467 MB of gene-call files, staging and the view). `/tmp` had 225 GB free at
11:17 and 224 GB at 12:10 (`df -h`). Nothing was written in the data root, under `/stor` or under
`/tmp/if2453/s2/inputs`.

## 8. Deviations from the card, and open questions

- **The batches are not `run-step` steps.** The card says one `run-step` process per batch. `run-step` publishes only
  files that a registered table or artifact claims (it refuses an empty or unclaimed staging), so a batch as a
  run-step would need an `ArtifactSpec` in `src/biobox/registry/schemas.py` (as card B-30 added `ingest_plan`),
  which is not among this card's files. A batch is instead a process like card W-35's parts: launched by build-stage's
  `_launch` within the allowance of the hour, `STOP` checked before each launch, under the memory cap and nice a
  run-step gives (`_batch_process`), writing only into the plan's work directory in scratch, with its CPU-minutes,
  wall time and peak in its `done.json` and the log. The merge is the run-step `genecall`. So the root holds
  nothing that says the step was cut, and lost scratch only makes batches run again (no rescan machinery). If the
  PI wants each batch's record in `runs/`, that is one artifact in `schemas.py` and a small change here.
- **The merge's provisional constant** was 10 times too high (section 5); it is now derived from S1's merge. It is my
  own constant, not a card's estimate.
- **My memory test's first bound for the merge was wrong** (section 6); corrected in the same session, before any
  merge of the test.
- Circular contigs in single and gv-meta mode: rates are an estimate (meta mode's ratio); none is built yet.
- Open: the plan is cut at the rates of this note; a stage of other genomes (for example many small circular
  prokaryote plasmids in single mode) is estimated, not measured. Each batch's `done.json` line in the log is the
  check: a batch more than 20% under its estimate in a real build is a finding (AGENTS §6).
