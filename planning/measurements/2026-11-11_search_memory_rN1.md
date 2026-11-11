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
