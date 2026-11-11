# A search on `rN.1` after card W-19: memory, time and the first `search status` (2026-11-11)

Question: card W-19 changed how a search reads its target (a block at a time, not held), checks its positive
control (a key lookup, not a set of every name), checks its pack (the build's record and the files' sizes, not a
full verify), where it writes (nothing into the pack directory) and when its pending record appears (before the
long work). Does the second trial's `search hmm` on `rN.1` still give the same answer, and what does it cost now?
Every number is *measured* on wilkcomp02, 2026-11-11, between 02:37 and 02:54 CST, by the command named beside it,
unless it is labelled *derived*. The baseline is the second trial
(`planning/measurements/2026-11-11_trial_search_rN1_second.md`, same host, same day, 01:39 to 01:49 CST).

## Result in one paragraph

Same answer, a twentieth of the memory, about half the wall time. The search worker peaked at 490 MiB resident
in the 15 s samples while it searched (26.27 GiB in the trial) and at 1.34 GiB in its run record, which is the
whole worker including the neighbourhood step after the search; the runner held 377 MiB until the exports
started (5.65 GiB in the trial). The command took 294 s (534 s), the run 258 s (342 s), and `search status`
first answered 10.4 s after launch (160 s). All 1,984 hit rows equal the trial's without the run id (every
score and E-value compared exactly), the 1,957 proteins are the same, and `hits.tsv`, `hits.faa`,
`neighbors.tsv` and the 2,581 locus tables are byte-identical to the trial's. The pack id is unchanged
(`71917a4167c7aebe1b7ff44d18a99306`): the build's verification record sits beside the manifest. The same search
against the pack built before W-19 (no record, bound read-only) verified it in full first and gave the same
bytes in 400 s.

## Setup

- Code: branch `task/W-19` at `fb6466f` (`origin/plan/task-cards` `950ffc2` with `origin/task/W-16` and
  `origin/trial/search-rN1-2` merged, then the W-19 commits), mirrored to `/tmp/if2453/ci/w19acc` by
  `sh scripts/pod_pytest.sh w19acc tests/r0/test_runbook.py` (`32 passed in 3.74s`).
- Data, read-only: root `/tmp/if2453/r10/root`, release `rN.1`, bound `:ro` in every command. Packs:
  `/tmp/if2453/w19/packs/release=rN.1/71917a4167c7aebe1b7ff44d18a99306/`, a copy (`cp -a`) of the W-16 pack
  `/tmp/if2453/w16/packs/...` (5,510,724,073 B, `du -sb`), and the W-16 pack itself; each bound `:ro` into the
  search that used it.
- Image `bio-box-test_74b3e05.sif`, sha256 `c30be585...c898`; profile `/tmp/if2453/w19/profile.yaml`
  (`profiles/wilkcomp02.yaml` with `tmp: /tmp/if2453/w19/tmp`).
- The trial's query and reference: `/tmp/if2453/trial2/work/RNA_pol.hmm`, the trial's run
  `/tmp/if2453/trial2/tmp/search/if2453/runs/run=87021d87e0f41e453c14d35763d51d17` and outputs
  `/tmp/if2453/trial2/out/rnap/`, read only.
- Load: `load.sh` (sum of `ps -u if2453 -o pcpu=` without the `ps` line, plus 400%, against 5,600% at night):
  4,174.5% (02:37), 4,690.3% (02:40), 5,025.9% (02:46), each with 400% under the cap. The pod suite of this card
  (`pod_pytest.sh W-19 --shards 2`) ran at the same time.
- Searches through `/tmp/if2453/w19/go.sh <unit> <pack dir> <out>`: `scripts/pod_run.sh --mem 32G --threads 4`
  with the trial's command (section 2 of the trial file, with `--pack-dir`, `--out` and `--profile` changed),
  a sampler unit `<unit>-sample` (`sample.sh`, the trial's: `VmHWM`/`VmRSS` of every process of the unit,
  `MemoryPeak`, `CPUUsageNSec`, every 15 s), then, as soon as the log names the run id, one `search status` in
  a capped scope (`bb.sh`), timed from just before `pod_run.sh`.

## 1. The pack's record, without a rebuild

The W-16 pack has no `pack_verified.json` (built before W-19). `pack build` with the same inputs on the copy:

```
bash /tmp/if2453/ci/w19acc/scripts/pod_run.sh --profile /tmp/if2453/w19/profile.yaml --image ... --digest ... \
  --mem 16G --threads 1 --bind /tmp/if2453 --bind /tmp/if2453/r10:/tmp/if2453/r10:ro --env PYTHONPATH=... \
  --log /tmp/if2453/w19/logs/w19-pack.log w19-pack pack build --root /tmp/if2453/r10/root --release rN.1 \
  --pack-dir /tmp/if2453/w19/packs
```

Log: `02:38:15 hashing the input files`, then `02:40:30 pack 71917a4167c7aebe1b7ff44d18a99306 is there, built
from the same inputs: nothing to do` and `pack ... verified; its record pack_verified.json added (the pack's files
are unchanged)`; `pod_run: w19-pack exited 0`. 154 s from the first log line (02:37:56) to the last (*derived*),
one CPU. The record (324 B) names `pack_id` `71917a41...` and `manifest_blake3` `71917a41...a2ea`: the manifest,
and so the pack id, are those of the W-16 build. No rebuild was needed.

## 2. The search, against the recorded pack

`bash /tmp/if2453/w19/go.sh w19-hmm /tmp/if2453/w19/packs /tmp/if2453/w19/out/rnap`, launched 02:40:50.6. The
log, without the nine slow-path notes (rough edge 18 of the trial, card W-20):

```
search run: afe0cce5b53d1086d2744ce31ad3f738 running; pending record /tmp/if2453/w19/tmp/search/if2453/runs/pending/run=afe0cce5b53d1086d2744ce31ad3f738.json
afe0cce5b53d1086d2744ce31ad3f738	ok
hits	/tmp/if2453/w19/out/rnap/hits.tsv	2581 genes of 1957 proteins
loci	/tmp/if2453/w19/out/rnap/loci	2581 tables
pod_run: w19-hmm exited 0
```

| Search | W-19 | Trial | Source |
|---|---|---|---|
| Wall time, whole command | 293.6 s | 534 s | *derived*: launch 02:40:50.6 (`go.sh`) to the log's last write 02:45:44.2 (`stat`) |
| Wall time, the run | 258.1 s | 341.7 s | `run.parquet` `wall_s` |
| Run id in the log | 6.4 s after launch | about 160 s | `go.sh` (polls the log every 0.5 s) |
| First `search status` answer | 10.4 s after launch, `{"run_id": "afe0cce5...", "status": "running"}` (the command itself 3.92 s, 328,248 kB) | 160 s or more | `go.sh` |
| CPU, whole unit | 9.08 CPU-minutes at least (545,058,720,000 ns) | 12.64 at least | `CPUUsageNSec` at the last sample, 02:45:40, 4 s before the end |
| CPU, the worker | 7.92 CPU-minutes | 9.00 | `run.parquet` `cpu_min_max_pid` |
| Peak memory, worker, while it searched | 501,952 kB = 490.2 MiB | 27,542,276 kB = 26.27 GiB | largest sampled `VmHWM` of the worker (`python -c import json, os, sys`), 02:42:23 |
| Peak memory, worker, whole process | 1,435,045,888 B = 1.34 GiB | 28,203,290,624 B = 26.27 GiB | `run.parquet` `peak_rss_bytes` (the measuring helper, card W-14) |
| Peak memory, runner, until the exports | 386,060 kB = 377.0 MiB | 5,920,924 kB = 5.65 GiB | sampled `VmHWM` of the `search hmm` process, constant from 02:41:08 to 02:45:09 |
| Peak memory, runner, with the exports | 1,407,396 kB = 1.34 GiB | (the trial's runner peak, 5.65 GiB, came before them) | sampled `VmHWM`, 02:45:25 |
| Unit memory | `MemoryPeak` 4,074,512,384 B (page cache included) | 31,688,663,040 B | sampler |

The worker's run-record peak (1.34 GiB) is above every sample of it while it searched (at most 490.2 MiB,
02:41:23 to 02:45:09, 16 samples). Its last sample is at 02:45:09 and the run closed after it; between them the
worker writes its hits and fetches the neighbourhoods of the 1,957 hit proteins (`nbr.fetch`). So the 1.34 GiB is
probably that step, as the separate `nbr` and `search hits` exports of the trial each reached about 1.3 GB
(trial file, section 2). Not measured separately; it does not depend on the target's size (card W-19 is about
the target), but it does depend on the number of hits and on the release's tapes.

Phases, from the samples (*derived* from the times and `CPUUsageNSec` deltas): 02:40:50 to 02:40:57 the
positive lookup of `search hmm` (about 4 s, the `contig_label` slow path) and the pending record; to about
02:41:15 the runner's check (record, sizes, key lookups in the two `unique` shards) and the staging of the two
shards (2,147,246,641 B, copied and hashed); 02:41:15 to 02:45:09 the worker: two passes over 3.34 G residues
(the model, then the shuffled negative control), at about 2.0 CPUs (`CPUUsageNSec` rises about 30 s per 15 s);
02:45:10 to 02:45:44 the four exports (in the runner, as before). The search itself took longer than in the
trial (about 234 s against 121 s, *derived*) because each block is read by one thread while the others wait;
but the 200 s of loading and the 160 s of checks are gone.

## 3. Against the trial's answer

`bash /tmp/if2453/w19/py.sh check1 /tmp/if2453/w19/check.py <trial run> <new run> <trial out> <new out>` (in the
image, 1.91 s, 102,212 kB):

| Check | Result |
|---|---|
| Hit rows (`hits.parquet`), every column but `run_id`, compared exactly | trial 1,984, new 1,984, equal: True |
| Run record: status, target size, rounds, final round, stop reason, pack id, shard sha256, Z | all equal (`ok`, 10,620,684, 1, 1, `max_rounds`, `71917a41...`, the two shards' sha256, 1e8) |
| Positive control | equal: T7 gp1 `2eef240a...` recovered |
| Negative control | 0 hits, expectation 0.000183, pass (its shuffle is seeded by the run id, so it differs from the trial's draw) |
| Proteins in `hits.tsv` | 1,957 and 1,957; 0 only in the trial, 0 only in the new run |
| `hits.tsv`, `hits.faa`, `neighbors.tsv` | byte-identical to the trial's |
| `loci/` | 2,581 files each, all byte-identical |

The trial's own check found its 1,957 proteins equal to `hmmsearch` (HMMER 3.4) on the export, so this run's are
too.

## 4. The same search against the pack built before W-19, read-only

`bash /tmp/if2453/w19/go.sh w19-old /tmp/if2453/w16/packs /tmp/if2453/w19/out/rnap_old`, the W-16 pack bound
`:ro`, launched 02:46:43.1. The log has, after the pending record line, `note: search pack 71917a41... of release
rN.1 was built before card W-19 and holds no record that its build verified it: verifying it in full for this
run (C-12); bio-box pack build with the same inputs adds the record without building again`, then `ok` and the
same output lines.

| Search | Value | Source |
|---|---|---|
| Wall time, whole command | 399.6 s | *derived*: 02:46:43.1 to the log's last write 02:53:22.6 |
| Wall time, the run | 363.1 s | `run.parquet` |
| First `search status` answer | 11.2 s after launch (`running`) | `go.sh` |
| The full verify | about 130 s (the worker's first sample is at 02:49:00) | *derived* from the samples |
| CPU, whole unit / worker | 10.92 CPU-minutes at least (02:53:17) / 7.97 | sampler / `run.parquet` |
| Peak memory, worker | 498,236 kB sampled (02:50:01); 1,577,553,920 B in the run record | sampler / `run.parquet` |
| Peak memory, runner, with the full verify and until the exports | 329,988 kB = 322.3 MiB | sampled `VmHWM` before 02:53 |
| Answer | `check.py`: 1,984 hit rows equal, 1,957 proteins, the four outputs and 2,581 locus tables byte-identical to the trial's | `check2` |

The verify's scratch went to the run's work directory, which the run removed.

## 5. What was written where

- W-16 pack: a listing of every path under `/tmp/if2453/w16/packs` with its type, size and modification time
  (`find -printf '%y %s %T@ %p'`, 15 lines) at 02:34 before anything ran and at 02:53 after the last run:
  identical (`diff` prints nothing). The trial had found `release=rN.1`'s time changed by the verify's scratch.
- Root: the top two levels (`find -mindepth 1 -maxdepth 2 -printf '%y %T@ %p'`, 56 entries) at 02:40, before the searches, and at 02:53:
  identical. The root was bound read-only in every command.
- Written, all under `/tmp/if2453/w19`: the pack copy (5.5 GB) with its `pack_verified.json`, the two runs
  under `tmp/search/if2453/runs/` (their work directories removed by the runs), `out/rnap` and `out/rnap_old`,
  `logs/`, the helper scripts and listings. Nothing under `/stor`. No `r10-*` unit was touched.

## Commands and last lines

| Time | Command | Last line |
|---|---|---|
| 02:36 | (laptop) `sh scripts/pod_pytest.sh w19acc tests/r0/test_runbook.py` | `32 passed in 3.74s` |
| 02:37 | `cp -a /tmp/if2453/w16/packs/release=rN.1/71917a41... /tmp/if2453/w19/packs/release=rN.1/`; `du -sb` | `5510724073 /tmp/if2453/w19/packs` |
| 02:37 | `pod_run.sh ... w19-pack pack build ...` (section 1) | `pod_run: w19-pack exited 0` |
| 02:40 | `bash go.sh w19-hmm /tmp/if2453/w19/packs /tmp/if2453/w19/out/rnap` | `first status answer 02:41:00.99 after 10.4 s: {"run_id": "afe0cce5...", "status": "running"}`; unit log `pod_run: w19-hmm exited 0` |
| 02:46 | `bash py.sh check1 check.py ...` | `loci: trial 2581 files, new 2581, all byte-identical: True` |
| 02:46 | `bash go.sh w19-old /tmp/if2453/w16/packs /tmp/if2453/w19/out/rnap_old` | `first status answer ... after 11.2 s`; unit log `pod_run: w19-old exited 0` |
| 02:53 | `bash py.sh check2 check.py ...` | `loci: trial 2581 files, new 2581, all byte-identical: True` |
| 02:53 | W-16 pack and root listings, `diff` | no output for either |

## What is left

1. The worker reads each block in one thread while the search threads wait: about 2 of 4 CPUs were busy, and
   the two passes took about 234 s against the trial's 121 s of searching (*derived*). Reading the next block
   while the current one is searched would bring the search back to about the trial's time.
2. The neighbourhood step at the end of the worker (probably the 1.34 GiB of its run record) and the exports in
   the runner (1.34 GiB) hold more than the search does; neither grows with the target, but both grow with the
   hits and read the release's tapes.
3. `worker._families` still reads every `l1_domains` row of the release into a dict when the target is
   `unique` and the recipe's `curation.foreign_domain_alarm` is on; `rN.1` has no domain table, so this run did
   not show it.
   On a release with Pfam hits it would grow with the release.

## Addendum, card W-21: the next block is read while the current one is searched (2026-11-11)

Question: with a reader thread one block ahead (`engine.prefetch`, two blocks in memory), does the same search on
`rN.1` still give the same answer, and how much of the 294 s and the 234 s of searching is left? Every number is
*measured* on wilkcomp02, 2026-11-11, 06:03 to 06:08 CST, unless it is labelled *derived*.

Setup. Code: branch `task/W-21` (`origin/plan/task-cards` `fed123f`, then the W-21 commit), mirrored to
`/tmp/if2453/ci/W-21` by `scripts/pod_pytest.sh`. Same image, profile (`tmp` moved to `/tmp/if2453/w21/tmp`),
root (`/tmp/if2453/r10/root`, bound `:ro`), pack (`/tmp/if2453/w19/packs`, the copy with its record, bound `:ro`),
query (`/tmp/if2453/trial2/work/RNA_pol.hmm`) and command as section 2, through `/tmp/if2453/w21/go.sh` (the W-19
`go.sh` with the paths changed): `pod_run.sh` unit `w21-hmm`, `--mem 32G --threads 4`, sampler every 15 s. Written
only under `/tmp/if2453/w21`; nothing under `/tmp/if2453/r10` or `/stor`; no `r10-*` unit touched.
Load before launch (06:03:47, Wednesday night, cap 5,600%): the sum of `ps -u if2453 -o pcpu=` without `ps` was
3,987.1% (load average 41.3), plus 400% is 4,387% (the helper's own cap line says 4,000%: it treats 05:00 to 20:00
as day, the card says 08:00 to 20:00 on weekdays). The CPUs were that busy with this account's other pod runs
(the pod suites of other cards), so the search ran on a loaded pod and its times are, if anything, long.

| Search | W-21 | W-19 | Source |
|---|---|---|---|
| Wall time, whole command | 244.6 s | 293.6 s (trial 534 s) | *derived*: launch 06:03:47.5 to the log's last write 06:07:52.1 (`stat`) |
| Wall time, the run | 210.2 s | 258.1 s | `run.parquet` `wall_s` |
| Search passes (staging done to the worker's neighbourhood step) | about 190 s | about 234 s (trial about 121 s) | *derived* from the samples: about 06:04:12 (launch plus the 25 s the W-19 run took to stage) to about 06:07:20, where the worker's resident size jumps (the `nbr.fetch` step) |
| CPUs busy while searching | about 2.6 of 4 | about 2.0 of 4 | *derived*: `CPUUsageNSec` rises 387.9 s over the 151 s from 06:04:19 to 06:06:50 |
| CPU, whole unit / the worker | 9.53 CPU-minutes at least (571,858,183,000 ns at 06:07:50) / 8.27 | 9.08 / 7.92 | sampler / `run.parquet` `cpu_min_max_pid` |
| Run id in the log / first `search status` | 6.4 s / 10.0 s after launch | 6.4 s / 10.4 s | `go.sh` |
| Peak memory, worker, while it searched | 737,692 kB = 720.4 MiB | 490.2 MiB | largest sampled `VmHWM` of the worker (06:06:05 on, the second pass) |
| Peak memory, worker, whole process | 1,916,227,584 B = 1.78 GiB | 1.34 GiB | `run.parquet` `peak_rss_bytes` |
| Peak memory, runner, until the exports / with them | 382,280 kB = 373.3 MiB / 1,286,708 kB = 1.23 GiB (06:07:50) | 377.0 MiB / 1.34 GiB | sampled `VmHWM` of `search hmm` |
| Unit memory | `MemoryPeak` 4,570,214,400 B (page cache included) | 4,074,512,384 B | sampler / `pod_run` usage line |

Answer. `check.py` (the W-19 run `afe0cce5...` against this run `b7d67df5...`): the 1,984 hit rows of
`hits.parquet` are equal in every column but the run id (every score and E-value compared exactly, order too);
status, target size, rounds, final round, stop reason, pack id, both shards' sha256 and Z equal; the positive
control (T7 gp1 `2eef240a...`) is recovered; the negative control has 0 hits; the 1,957 proteins of `hits.tsv` are
the same ones (0 only in either run); `neighbors.tsv` is byte-identical. `hits.tsv`, `hits.faa` and the locus
tables are not byte-identical to the W-19 run's, because the code under test includes card W-20 (merged into the
base before this card): its `hits.tsv` has no `tier` column, `hits.faa` has `strand=+1` in `example=`, and there
is a `README.txt` (card W-20's wording changes). Which other cells of `hits.tsv` and of the locus tables differ
was not itemised, and no pre-W-21 run on the W-20 code was made, so the byte-identity of the exports to the same
code without W-21 rests on `hits.parquet` above and on the unit tests that pin the bytes (W-07, W-19), which pass
unchanged. No hit row, score or order differs.

The W-19 memory test with its bound unchanged (the growth of the 100x search over the 10x search: at most
8 MiB): worker 338,198,528 B = 322.5 MiB at 10x (8,701 dark ORFs) and 339,566,592 B = 323.8 MiB at 100x (86,427),
1.3 MiB more (W-19 had 322.7 and 324.0 MiB, 1.2 MiB more), runner heap 16,798,805 and 16,798,848 B
(`pod_pytest.sh W-21 -s tests/r2/test_w19_search_memory.py`, 2026-11-11 05:12 to 05:14 CST). Two blocks of 2**16
residues are too small to show in the allocator noise; the bound needed no restating.

Reading. The search is faster but not back to the trial's 121 s: 190 s against 234 s on a loaded pod, with about
2.6 of 4 CPUs busy against about 2. The wall time of the command fell by 49 s (17%), the run by 48 s (19%). The
memory rose: while searching by about 230 MiB, which is about the second block (a block is about 2**25 residues
digitized, *estimate* 100 to 150 MB; the bound is two blocks), and in the run record by 0.44 GiB (1.34 to 1.78
GiB, +33%, above the card's 20% line, so reported here). The run record's peak is reached after the search, in
the neighbourhood step; that the reader's blocks (freed by then) leave it higher is probably the allocator keeping
the reader thread's arena, not tested. Neither number grows with the target (card W-19's test) and both are far
under the 32 GiB cap.
