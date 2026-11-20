# Who kills the peak helper: card A-10's diagnosis and the evidence card A-10a keeps (2026-11-20)

Numbers are *measured* unless labelled. GitHub numbers come from the logs of workflow `ci` (`gh run view <id> --log`,
attempt 1 unless said, read 2026-11-20); the logs are not copied here. Pod numbers are from wilkcomp02 on 2026-11-20.

## 1. Card A-10's diagnosis

Test `tests/r2/test_micro_s0_search.py::test_c16_an_oom_killed_search_records_a_peak_within_20_percent_of_its_cap`
failed once on Actions (run 44761392180, attempt 1, commit `b17935e`), with `recorded 596.7 MiB
(scope_memory_peak_sampled); cap 953.7 MiB; status oom`. Attempt 2 of the same run passed.

**The sampling interval is not the cause.** The search runner and run-step measure the same way: the same helper
(`runstep.measured`) runs inside the scope, and the same reader (`runstep.step_peak`) reads its file. When the
command ends, the helper reads the scope's `memory.peak` a final time (source `scope_memory_peak`). A `*_sampled`
source means the helper itself died before that final read, so all that is left is its last sample, taken
every 0.2 s.

Across the 31 job attempts read, the search OOM test recorded:

| Recorded | Source | Job attempts |
|---|---|---|
| 953.7 MiB, exactly the cap | `scope_memory_peak` (final) | 30 |
| 596.7 MiB | `scope_memory_peak_sampled` | 1 (44761392180 attempt 1) |

The step OOM test (`tests/r2/test_runstep.py::test_c16_an_oom_killed_step_records_a_peak_within_20_percent_of_its_cap`)
recorded `256.0 MiB (scope_memory_peak)`, the cap, in 31 of 31.

The 31 job attempts are 44761392180 (attempts 1 and 2) and these 29 runs:
44514423589, 44514832450, 44520684745, 44528317634, 44536752564, 44537776188,
44544970222, 44545556597, 44547221491, 44553797820, 44554083872, 44560235771, 44564776363, 44579413866,
44588896397, 44597860277, 44605025119, 44607147466, 44617991500, 44621676786, 44626292775, 44685837120,
44694435396, 44753536440, 44763961704, 44773315133, 44773693350, 44774354892, 44775879527.

**What the failing run shows.**
- The status was `oom` and the helper wrote no final file, so the runner saw exit -9: an exit of 86 or of 1 would
  have given a final file or the status `failed`. So the helper was SIGKILLed before its final read.
- Nothing in the runner sent that kill. A wall timeout would have ended `timeout`, a cancel `stopped`, and the
  runner was alive.
- The helper took no sample after 596.7 MiB, although one was due every 0.2 s and `memory.peak` never falls. So
  either the helper died at about the moment the worker was OOM-killed, or its sampler was blocked until then.

**Who killed it is not known.** The three candidates:
- the kernel's memcg OOM killer, taking the helper as a second victim;
- systemd escalating the scope's stop (`OOMPolicy=stop` by default) to SIGKILL;
- systemd-oomd.

The record cannot tell them apart. `cgroup_oom` is true either from `memory.events` or from the unit's
`Result=oom-kill`, and an oomd kill also sets that Result. The worker's `stderr.txt` was deleted with the work
directory.

**Why ten green runs prove nothing.** At 1 failure in 31, ten consecutive runs pass by chance with probability
(30/31)^10 = 0.72 (*derived*). A fix that only sometimes worked could not be told from no fix. The design choice
is with the PI (STATUS open item 20).

**Records on the pod.** The data root's `/tmp/if2453/r10/root/runs` holds 1,908 `"peak_rss_source": "child_maxrss"`
and no `*_sampled` (wilkcomp02, `grep -rhoE` over `runs/`, read-only). By the code it cannot hold one: run-step
publishes a record only when the step exits 0, and a helper that exits 0 has written its final read (runstep.py,
step 6). A sampled source appears only in two places:
- a failed attempt's record, `<tmp>/staging/<stage>.<step>.failed.json`, outside the root;
- a search run's `limits_json`, for runs that ended `oom`, `timeout` or killed `failed`.

## 2. The evidence card A-10a keeps

`runstep.peak_and_evidence` reads the helper's file as `step_peak` does. Only when the source is a sampled one
(`scope_memory_peak_sampled`, `child_hwm_sampled`) does it also return `runstep.peak_evidence`; otherwise it reads
nothing else and runs nothing. It is called once the last process of the run has been reaped: by run-step right
after `_execute`, and by the search runner right after `_wait`. The search runner calls it before
`scope_oom_killed`, which may `reset-failed` the scope unit and with that discard its Result.

The helper now names its cgroup in every write of its file (`"cgroup"`, None outside a scope). That file is
temporary, read once and removed, and no record changes for it.

Where the evidence goes:

| Run | Where | When |
|---|---|---|
| run-step | field `peak_evidence` of the failed attempt's record `<tmp>/staging/<stage>.<step>.failed.json` | sampled source only; a sampled source never reaches a published record (above) |
| search | key `peak_evidence` in the run row's `limits_json` | sampled source only; otherwise `limits_json` is byte for byte what it was |

The fields are read as they are and interpreted nowhere in the code:

| Field | What it is |
|---|---|
| `cgroup` | the scope's cgroup path, as the helper found it in `/proc/self/cgroup` (None outside a scope) |
| `cgroup_present` | whether that directory still existed when read (systemd removes a scope's cgroup soon after it empties) |
| `memory_events` | the scope's `memory.events` as `{key: count}`, or None when gone |
| `memory_peak`, `memory_max` | the scope's `memory.peak` (int) and `memory.max` (as written), or None |
| `unit`, `unit_show` | the scope unit, and `systemctl --user show -p LoadState -p ActiveState -p Result -p OOMPolicy -p MemoryPeak <unit>` as a dict (None without systemctl) |
| `systemd_oomd` | stdout of `systemctl is-active systemd-oomd` (None without systemctl) |
| `stderr_tail` | search only: the last 4,096 characters of the worker's `stderr.txt`, where the helper's own stderr also goes |

## 3. How to read a sampled case

Among the kernel's `memory.events` keys:
- `oom` counts the scope's own charges that failed at `memory.max`;
- `oom_kill` counts the processes of the scope that any OOM killer killed;
- `oom_group_kill` counts kills of the whole group (memory.oom.group).

The helper and the worker are the only processes in the scope. Read a sampled case as follows (*inferred* from the
kernel and systemd documentation, not yet observed):

| What the evidence shows | Reading |
|---|---|
| `oom_kill` 2 or more, `oom_group_kill` 0 | the kernel's OOM killer killed the helper too, as a second victim; option (c) of card A-10 (the helper charges no memory) or (a) |
| `oom_group_kill` 1 or more | the scope was killed as a group (memory.oom.group set, so `OOMPolicy=kill`): look at `unit_show.OOMPolicy` |
| `oom_kill` 1 and `systemd_oomd` `active` | the helper was not killed by the kernel; oomd is a candidate (it kills a cgroup without the kernel's count) |
| `oom_kill` 1, oomd not active, `unit_show.OOMPolicy` `stop` | the helper was killed outside the kernel's OOM killer: systemd's stop of the scope is the candidate; option (b) or `OOMPolicy=continue` |
| `cgroup_present` false | the counts are lost; only `unit_show.Result` remains, which an oomd kill sets to `oom-kill` too |
| a traceback in `stderr_tail` from the helper (`Exception in thread ... (sample)`) | the helper's sampler failed while the scope was full |

**Known killers, for comparison.** Two tests kill the helper themselves and are sampled on every Actions run:
- `test_c16_a_search_killed_at_its_wall_time_records_its_peak`: the runner's wall timeout SIGKILLs the process
  group;
- `test_c16_a_killed_step_records_its_peak_and_where_it_came_from[group]`: the test SIGKILLs the process group.

Their evidence lines (`C16PeakEvidence`) now print `peak_evidence` too. That gives the usual values of a run whose
helper was killed by someone other than the kernel: no OOM kill counted. In the two OOM tests the evidence line
prints `peak_evidence null` whenever the final read was taken, which is about 30 jobs in 31.

## 4. Tests and where they ran

- New, `tests/r1/test_a10a_peak_evidence.py`, on the pod (wilkcomp02, image `bio-box-test_faec1f4.sif`). It uses a
  fake cgroup tree, a stand-in `systemctl` that logs its calls, and run-step end to end with a stand-in
  `systemd-run`. Before the change: 10 failed, 1 passed (`FFFFFFFFFF.`). The one that passed is the step killed
  alone, where no evidence is expected and none was kept. After the change all 11 pass.
- The two OOM-killer tests skip on the pod (no `systemd-run` inside the image) and run only on Actions. Card A-10a's
  acceptance is the branch run there until the sampled case is seen (about one job in 31) or 60 jobs; this
  section is to be completed with what it showed.

## 5. The first run on Actions: 60 iterations of the four tests alone (the orchestrator, 2026-11-20)

Workflow `peak-evidence` (card A-10b), run 44790766798 on `ubuntu-24.04` (runner image 20261116.376.1), commit
`plan/task-cards`' head pushed to `peak-evidence/run-1`, 18:13 to 18:41 UTC: 60 iterations run, 0 failed (*measured*,
the job's summary and `gh run view --log`).

- **The sampled case did not occur**: 0 of 60 iterations recorded `scope_memory_peak_sampled` for the OOM-killed
  search, and 0 of 60 for the OOM-killed step; every one printed `recorded 953.7 MiB (scope_memory_peak)` and
  `peak_evidence null`. At the rate seen under the whole suite (1 job in 31) the chance of none in 60 is 14%
  (*derived*: (30/31)^60), so this neither shows nor excludes that the case needs the suite's load. The evidence for
  it now waits for an ordinary `ci` job, whose log prints the same lines.
- **`systemd-oomd` is not the killer on this runner image**: `"systemd_oomd": "inactive"` in all 120 evidence
  records of the run (*measured*: the two tests that kill their own helper, 60 iterations each). One of card A-10's
  three candidates is out, as long as the image keeps it inactive.
- **The step's scope sets `OOMPolicy=stop`** (*measured*, the 60 records of the killed step: unit `active`,
  `LoadState` `loaded`, `Result` `success`, `memory.events` all zero, `memory.max` 1073741824, `MemoryPeak` equal to
  the cgroup's `memory.peak`). So a kernel OOM kill in the scope makes systemd stop the unit; whether that stop
  reaches the helper before its final write is what a sampled OOM case would show.
- **A limit of the evidence for searches, found here**: in all 60 records of the search killed at its wall time the
  scope was already gone when the evidence was read (`"cgroup_present": false`, `memory_events` null, unit
  `LoadState` `not-found`, `Result` `success`, `stderr_tail` empty), where the step's scope was still there in all
  60 (`"cgroup_present": true`). The search runner reads after its worker's scope has emptied and been collected. A
  sampled OOM case of a **search**, the one that failed on 2026-11-20, would therefore keep only `systemd_oomd`, the
  worker's `stderr_tail`, and the unit's `Result` if a failed unit outlives its cgroup (not observed yet: the unit of
  a scope that ended without failure is gone). Reading the counters for a search needs them taken while the scope
  still has a process in it, which is card A-10's design question (a keeper process, option (a)), not a change of
  order.
- The workflow's own parts that could not be run before the merge worked: the input's default on a push, the loop,
  the per-iteration evidence lines, the totals, the upload and the final check ("all 60 iterations passed").
