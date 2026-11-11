# Runbook: bio-box on the pod

Commands only. Each command in a `bash` block below was run on the pod on 2026-11-10; Appendix A
lists the last line each one printed. The rules behind the limits are in `AGENTS.md` section 7. If
a step here and a rule there disagree, the rule wins: stop and tell the PI.

Where a command runs:

- Everything runs on the pod host, in `bash`, except section 8 (`scripts/pod_pytest.sh`), which
  runs on the laptop and reaches the pod by `ssh`.
- `bio-box ...` inside the image is the function defined in section 0. A long job is started with
  `scripts/pod_run.sh`, which opens the memory-capped systemd unit on the host and runs
  `python -m biobox.cli ...` inside the image.
- Values like `S0`, `S1`, `rC6.1`, `rC6.2`, `c06-build` and `w13-search` are the fixture names used
  when the commands were tested. Replace them with your own stage, release and unit names.

The commands were tested on a small fixture build in a scratch data root with a scratch copy of
the host profile (sections 1 to 6 and 8 in `/tmp/if2453/c06`, card C-06; section 7 and the search
part of section 9 in `/tmp/if2453/w13`, card W-13). Where the real profile differs, the text says so. The real profile is
`profiles/wilkcomp02.yaml`; section 10 copies its values.

## 0. Set up the shell

Paste this once per shell. Lines marked `real:` give the value for a real build.

```bash
export TREE=/tmp/if2453/ci/C-06                      # real: the checkout that `git push pod` updated
export IMAGE="$HOME/bio-box/images/bio-box-test_74b3e05.sif"
export DIGEST=sha256:c30be58557d756095346308c58cc3098af4bbe66404f67ab87be3f64f147c898
export PYDEPS=/tmp/if2453/pilot_s1/pydeps
export PYTHONPATH_POD="$TREE/src:$PYDEPS"            # the image is old: the source is bound in (section 11)
export PROFILE=/tmp/if2453/c06/profile.yaml          # real: "$TREE/profiles/wilkcomp02.yaml"
export MEM=8G                                        # real: 200G (the profile's mem_scope)
export WORK=/tmp/if2453/c06                          # your scratch directory for inputs and outputs
export BBTMP=/tmp/if2453/c06/tmp                     # the profile's tmp (real: /tmp/if2453/bio-box)
export ROOT=/tmp/if2453/c06/root                     # the data root you build into
SCOPE=(systemd-run --user --scope --quiet -p MemoryMax="$MEM" -p MemorySwapMax=0 nice -n 19)
SING=(singularity exec --cleanenv --env PYTHONPATH="$PYTHONPATH_POD" --env OMP_NUM_THREADS=2 --bind /tmp/if2453 --bind /sys/fs/cgroup:/sys/fs/cgroup:ro)
bio-box() { "${SCOPE[@]}" "${SING[@]}" "$IMAGE" python -m biobox.cli "$@"; }
biopy() { "${SCOPE[@]}" "${SING[@]}" "$IMAGE" python "$@"; }
```

`bio-box` and `biopy` run a short command in its own capped scope at nice 19 (AGENTS 7.2). Do not
use them for a long job: use `scripts/pod_run.sh` (sections 2 and 7).

## 1. Before you start

1. Look for a `STOP` file. If this prints a path, do not start anything.

```bash
ls -d "$ROOT/STOP"
```

`No such file or directory` is the good answer.

2. Look for other compute. A load average above 8 means someone is computing: the night cap does
not apply, and you should wait or take the day cap (`run-step` does this by itself).

```bash
uptime
top -b -n 1 -o %CPU -w 120 > "$WORK/top.txt"
head -n 15 "$WORK/top.txt"
systemctl --user list-units --type=service --type=scope --no-pager
```

The unit list shows the units of this account. Other agents and the PI use the same account, so
a unit you did not start is not yours. Processes of other users are in `top`.

3. Check free space where the work goes.

```bash
df -h /tmp
```

4. Record the host and image. The `cgroup` block must show your cap (`8 GiB` here, `200 GiB` for
the real profile). `paths` shows the profile's `tmp`, `durable` and `bulk`; `missing` is normal
before the first build.

```bash
biopy -m biobox.doctor --profile "$PROFILE"
```

The card calls this `bio-box doctor`. There is no such subcommand; the command is the module
`biobox.doctor`.

## 2. Build a stage and a release with one command

`build-stage` runs inputs, ingest, gene calling, IGRs, tapes, one domain scan per `--library`,
clustering, tokens, postings, `contig_dim`, manifest and catalog. Each step is a `run-step`. The
two input files are JSON lists: `--genomes` (fields of `GenomeInput`, contig uids as hex) and
`--objects` (fields of `SourceObject`). `--library` takes a name, a version, an HMM file and a
cutoff, and repeats. Appendix B makes fixture inputs.

Launch it once, detached. Use a new unit name for every launch.

```bash
bash "$TREE/scripts/pod_run.sh" --profile "$PROFILE" --image "$IMAGE" --digest "$DIGEST" --bind /tmp/if2453 --env PYTHONPATH="$PYTHONPATH_POD" c06-build build-stage --profile "$PROFILE" --root "$ROOT" --stage S0 --release rC6.1 --stages S0 --genomes "$WORK/inputs/genomes-S0.json" --objects "$WORK/inputs/objects-S0.json" --library test 1 "$WORK/inputs/lib.hmm" gathering
```

It returns at once and prints the log path. The memory cap, the thread cap and nice come from the
profile. `--digest` refuses to start if the image is not the one you named. Pass one more `--bind`
for every directory outside `/tmp/if2453` that the build reads.

Where the log is: `$BBTMP/logs/<unit>.log`, the profile's `tmp` plus `logs`.

Watch it while it runs:

```bash
systemctl --user status c06-build --no-pager
tail -n 3 "$BBTMP/logs/c06-build.log"
```

Each step prints one JSON line to the log as it ends. Every `run-step` and `build-stage` line starts
with a timestamp (local time and its UTC offset, for example `2026-11-10T07:21:57-06:00`) and is
written as it happens, so the last line of the log shows where the build is: the step that was last
`published` or `skipped`, or the shard that last ended.

How to tell it finished:

```bash
tail -n 3 "$BBTMP/logs/c06-build.log"
grep "steps done" "$BBTMP/logs/c06-build.log"
ls "$ROOT/runs/stage=S0"
```

- The last line is `pod_run: c06-build exited 0`. The number is the exit code. The line before it,
  `pod_run: c06-build usage CPUUsageNSec=<ns> MemoryPeak=<bytes>`, is what the unit cost (card W-18):
  CPU time in nanoseconds and the peak memory of the whole unit, read from `systemctl --user show` just
  before the unit ends, because `--collect` discards them afterwards. It is absent when systemd keeps no
  accounting for the unit.
- The log has a line ending in `build-stage: stage S0 and release rC6.1: all N steps done`, and a
  `run-step: ... published N files` line for each step that runs a command (the `features` slot is
  empty). `--library` is optional and repeats: N is 12 plus one per `--library` (derived from
  `plan` in `src/biobox/stage.py`, which prints `len(plan)`; measured, wilkcomp02, 2026-11-10: units
  `r10-S1-1` and `r10-P1-1`, which ran without `--library`, printed `all 12 steps done`).
  With one library, as in the commands here, the line is `all 13 steps done` and 12 steps run a command;
  with none it is `all 12 steps done` and 11 run a command. A scan that runs in shards adds plan, shard and
  merge steps, but the count in the line stays the same.
- `systemctl --user status` says `could not be found`: the unit has ended.
- `pod_run: c06-build received SIGTERM` just before the last line: the unit was stopped
  (`systemctl --user stop`, section 3). Relaunch to resume (section 4).
- A log whose last line is not `pod_run: ... exited`: the unit was killed without a chance to
  write it (SIGKILL, or a reboot). Go to section 4.

Exit codes: `0` done. `75` a `STOP` file (section 3). `2` refused before or at a step: read the
`build-stage:` or `run-step:` line in the log. Any other number is the failing step's own code, and
the log names the step. `128` plus a number means a signal was forwarded.

A second identical call starts no step command. Every step prints `skipped: run record and outputs
present`:

```bash
bash "$TREE/scripts/pod_run.sh" --profile "$PROFILE" --image "$IMAGE" --digest "$DIGEST" --bind /tmp/if2453 --env PYTHONPATH="$PYTHONPATH_POD" c06-build-again build-stage --profile "$PROFILE" --root "$ROOT" --stage S0 --release rC6.1 --stages S0 --genomes "$WORK/inputs/genomes-S0.json" --objects "$WORK/inputs/objects-S0.json" --library test 1 "$WORK/inputs/lib.hmm" gathering
```

A second stage and a release over both stages: change `--stage` to `S1`, `--release` to `rC6.2`,
`--stages` to `S0,S1`, and the input files to the second stage's. The first release stays valid.

**Domain scans in shards.** A scan estimated above 80% of the profile's `cpu_min_per_process` runs
as a plan step, one step per shard and a merge step. The plan cuts the stage's proteins into shards
estimated at 80% of the budget at most; it depends only on the inputs and the budget. Several
shards run at once, each its own `run-step`:

- `--scan-cpus C`: threads per shard. Default 4.
- `--scan-jobs J`: shards at once. Default: the allowance when `build-stage` starts, divided by C,
  rounded down.
- The allowance is the thread count `run-step` allows: 40 by day; 56 at night and at weekends, only
  while the load average is at most 8.
- J × C above the allowance when `build-stage` starts is refused before any step (exit 2). The
  message names J, C and the allowance.
- The allowance is read again before each step and before each shard launch. A build started at
  night keeps to 40 threads once the day window opens: shards already running finish, and new ones
  start only while the shards in flight fit the allowance, so fewer run at once. If one shard's C
  is above the allowance, the launch waits (the log says `holding it`) and starts when the hour or
  the load allows. The same holds for a step whose `--threads` is above the allowance; the default
  `--threads` (the smaller cap) never waits.
- The output is the same, byte for byte, for any J and C, and whenever the allowance changes.

The unit's CPU quota is `pod_run.sh --threads`, which defaults to the profile's night cap
(`max_threads_night`, 56), the most any hour allows. It is only a ceiling, fixed at launch:
`build-stage` holds the limit of the hour itself, before each step and before each shard launch.
This fixture launch runs two shards of one thread at once; the scratch profile allows 2 threads:

```bash
bash "$TREE/scripts/pod_run.sh" --profile "$PROFILE" --image "$IMAGE" --digest "$DIGEST" --bind /tmp/if2453 --env PYTHONPATH="$PYTHONPATH_POD" --threads 2 c11-build build-stage --profile "$PROFILE" --root "$ROOT" --stage S0 --release rC6.1 --stages S0 --genomes "$WORK/inputs/genomes-S0.json" --objects "$WORK/inputs/objects-S0.json" --library test 1 "$WORK/inputs/lib.hmm" gathering --scan-jobs 2 --scan-cpus 1
```

For a real Pfam scan, with the real profile, do not take the whole cap. The cap (40 threads by day, 56 at
night and at weekends, section 10) is a total over everything of ours on the pod: your build, other
agents' builds and tests, searches, and the PI's work, all together (AGENTS 7.2). `build-stage` counts
only its own threads, so it cannot see the others; you must. Count the threads the other jobs of ours use
now (`systemctl --user list-units`, `top`, section 1) and choose J so that all together stay at or under the cap:

```text
J = floor((cap - other threads) / C)       C = --scan-cpus, default 4; J at least 1
```

`cap` is the one for the hour (the allowance). `other threads` is what else of ours runs, and what you
will start during the build. A worked example (derived, from the formula; the 40 is the profile's
`max_threads_day`): by day, with cap 40, one other build that runs 16 threads and C = 4, J =
floor((40 - 16) / 4) = 6, so `--scan-jobs 6 --scan-cpus 4`, 24 threads, 40 in all. With nothing else
running, the same formula gives J = 10, which is the whole cap: leave a margin if anything may start
during the hours the scan runs. A smaller J only makes the scan take longer; the output is the same.

**Clustering's own cap (DEC-25).** A release's clustering is one process, and on a large release it
passes 300 CPU-minutes: 2,195.91 CPU-minutes for 10,620,684 proteins (measured, wilkcomp02,
2026-11-10, run record `rN.1.cluster`). DEC-25 exempts this one step from the 300 CPU-minute rule
under its own cap, the profile's `cluster_cpu_min_max` (20000 in `profiles/wilkcomp02.yaml`, an
estimate: about five times that run; the PI may change it). Every other step keeps the 300-minute rule.

- `build-stage` estimates clustering's CPU-minutes before it starts it, from the release's protein
  count: a power law through 284.10 CPU-minutes at 2,053,075 proteins and 2,195.91 at 10,620,684
  (measured, wilkcomp02, 2026-11-10, `rN.1b` and `rN.1`; the fit is derived and extrapolated above
  10.6 M proteins). At the cap of 20000 it refuses a release of more than about 62.7 M proteins
  (derived: 10,620,684 × (20000 / 2195.91)^(1 / 1.244)).
- When every stage of the release has its gene-call run record, the check is made before any step;
  otherwise just before the cluster step's command, once this call's stage is gene-called.
- An estimate over `cluster_cpu_min_max` is refused, exit 2. The log line starts `build-stage: refused`
  and names `cluster_cpu_min_max`.
- A profile without `cluster_cpu_min_max` (the scratch profile of Appendix B) gives no exemption:
  clustering estimated over `cpu_min_per_process` is refused, exit 2, and the log line says so.
- With the key, the cluster step's run record holds `cpu_budget_reason: exempt: DEC-25`, and
  `over_cpu_budget` is still true when the step used more than 300 CPU-minutes (section 5).
- The exemption ends when card B-37 splits clustering into processes under the budget.

**Clustering in parts (card B-41).** When clustering in one process is estimated above 80% of the
profile's `cpu_min_per_process` (240 CPU-minutes of 300), `build-stage` clusters the release in
parts, each step its own `run-step` under `runs/stage=<release>/`:

- `cluster.plan` (the proteins into MMseqs2's database), then one step per MMseqs2 module call that
  the two `cluster` workflow calls make (`cluster.<level>.<kk>.<module>`, 42 in all, with the step
  between and after them), then `cluster` (the tables) last, under the cluster step's own name.
- Each per-query call (`prefilter`, `align`, `rescorediagonal`) has a `.cut` step and, when cut into
  more than one part, one `.part-<p>` step per part before the call's own step, which joins them.
  Parts run several at once like scan shards, with `--scan-jobs` and `--scan-cpus`, within the hour's
  allowance. A part is estimated at 80% of the budget at most, from the call's measured share.
- `--cluster-parts N` cuts every per-query call into N parts whatever the estimate (for trials).
- The files are those of clustering in one process, byte for byte (tested on the fixture and a
  constructed set; not yet on the pod). The steps' records are in `l2/cluster_plan/release=<R>/`,
  never in a manifest.
- MMseqs2's databases are under `$BBTMP/scratch/biobox-cluster/<release>-<key>/`, and are removed when
  `cluster` publishes. Their size is not yet measured.
- The DEC-25 estimate and its refusal above still come first, until the exemption is removed.

A shard that fails stops new launches. Shards already running go on to their own end and publish if
no `STOP` file exists then, and `build-stage` exits with the failed shard's code. A `STOP` file is
different: a shard whose command ends while `STOP` exists publishes nothing (section 3). The next call
runs only the shards that have no run record, then the merge. While shards run, the log holds one
`run-step:` line per shard as it ends, in the order they end.

## 3. Pause with STOP, and resume

Pause: create the file in the data root.

```bash
touch "$ROOT/STOP"
```

- A new launch exits 75 and runs nothing. The log says `STOP exists; not starting`.
- A running build checks before each step and exits 75 at the next step boundary. During a scan
  in shards it checks before each shard launch, launches no more, and waits for the shards that run
  (`_launch` in `src/biobox/stage.py`: `launching no more shards, waiting for the N running`).
- A step that is running when the file appears finishes its command. `run-step` then publishes
  nothing, leaves its staging directory and exits 75. The log says `STOP appeared; not publishing`.
  The next call runs that step again from the start. A scan shard is such a step (each shard is its
  own `run-step`), so this is the one rule for shards too: a shard that ended before the file appeared
  is published and kept; a shard still running when the file appears runs its command to the end and
  publishes nothing (the check is after the command, in `run_step` of `src/biobox/runstep.py`), and
  `build-stage` launches no more and exits 75 once they have ended. Those shards run again from the
  start in the next call; the shards with a run record do not. A shard is cut to at most 80% of 300 CPU-minutes (240, derived: 300 × 0.8), so
  a `STOP` costs at most that much per running shard, as an upper bound from the plan, not measured.
- So a `STOP` inside a step costs that step's whole work, done twice. On the pod a `STOP` during
  stage S1a's gene calling cost 57.7 s and 4.98 CPU-minutes (measured, wilkcomp02, 2026-11-10,
  card C-10, the rerun's run record); during a release's clustering it would cost up to 1,575.2 s
  and 195.08 CPU-minutes (measured the same day, `rC.2.cluster` at 8 threads). To lose less, create
  `STOP` just after the log shows a `run-step: ... published` line. This is the rule, not a fault.

To stop at once, without waiting for the step, stop the unit. `run-step` forwards the SIGTERM to
the step and publishes nothing from it.

```bash
systemctl --user stop c06-build
```

The log of a stopped unit ends with `pod_run: c06-build received SIGTERM` and then
`pod_run: c06-build exited 143`. The code is the command's own: a search that is stopped records the
run as `stopped` and exits 1 (section 7).

Resume: remove the file, then launch the same build command again (section 2), with a new unit name.

```bash
rm "$ROOT/STOP"
```

Steps with a run record are skipped. The step that stopped runs again.

## 4. After a kill or a reboot

Relaunch the same command from section 2 with a new unit name. That is the whole recovery.

First check that the data root and `tmp` are still there. After a reboot `/tmp` may have been
cleared. If `$ROOT` was under `/tmp` and is gone, the build is gone with it. The search packs and
runs are under `tmp` too (`$BBTMP/packs`, `$BBTMP/search/$USER`, section 7): if they are gone,
build the pack again; the runs are lost. A search whose runner was killed reads as `failed` with
`runner exited without closing the run`; start it again.

```bash
ls -d "$ROOT" "$BBTMP"
```

What `run-step` leaves, and what it does at the next start:

| Left behind | When | What `run-step` does |
|---|---|---|
| `$BBTMP/staging/<stage>.<step>/` | the step was killed, failed or stopped | removes it and runs the step again |
| `$BBTMP/scratch/<stage>.<step>/` (a step's own scratch) | the step was killed, failed or stopped | `build-stage` removes it before that step's next attempt (log: `removing ..., left by an earlier attempt`), and removes it once the step succeeds. In C-10 a stopped clustering step left 1.2 GB there (measured, wilkcomp02, 2026-11-10) |
| a scan shard's run missing from `$BBTMP/scratch/biobox-domain-scan/<plan>/` while its run record exists | scratch was cleared (a reboot) or damaged between the shards and the merge | `build-stage`, before the merge, scans exactly those shards again into scratch (log: `scanning it again`); the new run must have the sha256 the shard's published record holds, else the build stops naming the shard (tell the PI) |
| a clustering step's database missing from `$BBTMP/scratch/biobox-cluster/<release>-<key>/` while its record exists (card B-41) | scratch was cleared (a reboot) or damaged during clustering in parts | `build-stage`, before it goes on, makes exactly those steps again, in order (log: `making them again in order`); each must reproduce its published record, else the build stops naming the step (tell the PI) |
| `$ROOT/runs/stage=<S>/<step>.json` (run record) | the step finished | skips the step if every listed output is in the root with its recorded size |
| `$ROOT/runs/stage=<S>/<step>.publishing.json` (intent file) | killed while publishing | finishes the publish, or undoes it, then reruns (see below) |
| `$ROOT/.incoming/<stage>.<step>/` | killed while publishing | uses it to finish the publish, then removes it |
| `$ROOT/.incoming/` (empty) | a step was killed after removing its own directory there, or the root's filesystem cannot lock a directory (NFS) | nothing; the next step to publish removes it when it is empty and no other step is creating into it |

An intent file lists every output with its sha256. At the next start, an output that is in the
root with the listed sha256 is kept, and one still in `.incoming` is moved in. If every output is
then in place, `run-step` writes the run record and exits 0 without running the step. Otherwise it
deletes only the root files that match the intent's sha256, deletes the intent and runs the step
again. A root file at a listed path with another sha256 is not this step's: `run-step` exits 2,
names it and changes nothing. Do not touch such a file; tell the PI.

List the intent files:

```bash
find "$ROOT/runs" -name '*.publishing.json'
```

When a staging directory may be removed by hand: only when no unit of yours is running, only under
`$BBTMP/staging/`, and only for a step that has no intent file. `run-step` removes it itself at
the next start, so removing it is only to free space (section 9).

When it must not be removed: the staging directory, the `.incoming` directory and every root file
of a step that has an intent file and no run record. They are what the recovery reads. Leave them
and relaunch.

If a step failed with its own error (not a kill), the log names the staging directory it left
for inspection. Read what the step wrote, fix the cause, and relaunch.

## 5. Read a run record and the ledger

A run record is one JSON file per step: `$ROOT/runs/stage=<S>/<step>.json`. A release step is
under `stage=<release>`. `runs/` is not in the release manifest, because it holds times and host
names. Fields:

| Field | Meaning |
|---|---|
| `stage`, `step` | what ran |
| `argv` | the command, as started (after `nice`) |
| `threads`, `nice` | the thread count and nice value |
| `started_at` | UTC start time |
| `wall_s`, `cpu_min`, `peak_rss_bytes` | wall seconds, CPU-minutes, peak resident bytes of the step command |
| `mem_cap` | `how` is `inherited` (the unit's cap, with `memory_max_bytes`), `opened` or `none` |
| `exit_code`, `over_cpu_budget` | `over_cpu_budget` is true above the profile's `cpu_min_per_process` |
| `cpu_budget_reason` | only in a release's `cluster` record, when the profile sets `cluster_cpu_min_max`: `exempt: DEC-25` (section 2) |
| `image_digest`, `host` | the image sha256 that `pod_run.sh --digest` passed in, and the host |
| `outputs` | each published file: `path`, `bytes`, `sha256` |

These commands run `$WORK/rb.py`, a script file saved from Appendix B. It is a Python file and not a
`python -c` string, so a command survives one level of quoting, such as `ssh host 'bash -c "..."'`; a
`python -c` string with quotes inside it does not (`unexpected EOF`, found in card C-14).

One record:

```bash
biopy "$WORK/rb.py" record "$ROOT/runs/stage=S0/genecall.json"
```

Every record in the root, in start order:

```bash
biopy "$WORK/rb.py" records "$ROOT"
```

The ledger is a Parquet table, one file per step: `ledger/stage=<S>/step=<name>/part-0.parquet`,
and for clustering `ledger/release=<R>/step=cluster/part-0.parquet`. One row per (object, reason):
`object_key`, `object_kind`, `n_in`, `n_kept`, `n_merged`, `n_rejected`, `reason_code`,
`tool_versions_json` (the versions of the tools that shaped the output). `n_in`, `n_kept` and
`n_merged` repeat on every reason row of an object. `n_rejected` counts that row's reason only.
A release's catalog has the view `ledger`. Totals per step, with each object counted once:

```bash
biopy "$WORK/rb.py" ledger "$ROOT" rC6.2
```

Two checks. Balance (L-01): per object, `n_in = n_kept + n_merged + sum(n_rejected)`. Chain (L-03):
along each edge, the downstream `n_in` equals the upstream `n_kept`. Each returns a list of
problems; an empty list is a pass. The chain check needs a release whose stages are one stage.

```bash
biopy "$WORK/rb.py" checks "$ROOT" S0 rC6.1 ingest,genecall,domains-test-1,cluster
```

## 6. Verify a release

`release verify` re-reads every file the release's manifest lists and compares size and sha256. It
also reports any unlisted file in the release's scope. It writes nothing. It prints each problem,
then `release rC6.2: verified (...)` with the stages, the file count and the manifest hash. Exit 0
is a pass, 1 a failure (a missing manifest is one).

```bash
bio-box release verify --root "$ROOT" --release rC6.2
```

With `--previous`, it also checks that the later release kept the earlier one's files: every file
of the earlier manifest, except the earlier release's own `release=` files, is in the later manifest
with the same sha256 (append-only, R-09). It then prints `append-only against rC6.1: ok`, or each
file as `removed` or `rewritten`.

```bash
bio-box release verify --root "$ROOT" --release rC6.2 --previous rC6.1
```

Open the catalog and read from it. The catalog is a DuckDB file that names the root it was built
from; it reads the root's Parquet files by the manifest's list.

```bash
biopy "$WORK/rb.py" catalog "$ROOT" rC6.2
```

An earlier release must still verify after a later one is built. Check `rC6.1` the same way:

```bash
bio-box release verify --root "$ROOT" --release rC6.1
```

## 7. Search

The search path writes nothing under the data root (card W-03): a user who can read a release may
not be able to write its root. The search pack and the runs live outside the data root, under the
profile's `tmp`:

- the pack of release `R`: `$BBTMP/packs/release=<R>/<pack id>/`, one per host (a pack is a cache
  of a release, shared by every user);
- the runs: `$BBTMP/search/$USER/runs/run=<run id>/`, one directory per user; each run's pending
  record is `$BBTMP/search/$USER/runs/pending/run=<run id>.json`;
- the runs' work directories: `$BBTMP/search/$USER/work/<run id>/`, removed when the run closes.

Every `pack`, `search` and `nbr` command below passes `--profile "$PROFILE"`, which gives these
directories. `scripts/pod_run.sh` does not export `BIOBOX_PROFILE`, and the `bio-box` function of
section 0 does not either, so a command without `--profile` falls back to `/tmp/atd/packs` and
`/tmp/atd/search`, which are not ours. In place of `--profile`, name the directories: `--pack-dir`
(pack, run, hits, nbr), `--runs` (run, status, cancel, hits, nbr) and `--work` (run). A directory
under the data root is refused before anything is read.

Search commands that run for more than a moment go through `scripts/pod_run.sh`, so they run in a
capped unit.

**Pack.** Build the search pack of a release once. `--model` is `<library>-<version>=<HMM file>`,
once for each domain library of the release; it is optional, and a release with no domain library
(the 10 Gbp release `rN.1` has none) is built without it. A search seeded by an HMM file or by
sequences needs no model in the pack. A second build from the same inputs finds the pack
and changes nothing; a pack is never overwritten. A search whose pack is missing exits 2 and names
the `pack build` command.

```bash
bash "$TREE/scripts/pod_run.sh" --profile "$PROFILE" --image "$IMAGE" --digest "$DIGEST" --bind /tmp/if2453 --env PYTHONPATH="$PYTHONPATH_POD" w13-pack pack build --root "$ROOT" --release rC6.1 --model test-1="$WORK/inputs/lib.hmm" --profile "$PROFILE"
bio-box pack verify --root "$ROOT" --release rC6.1 --profile "$PROFILE"
```

**Recipe.** A YAML file; its format is in the docstring of `src/biobox/search/recipe.py`. A recipe
without a positive control is refused. Appendix B writes a fixture recipe.

**Run.** Pass `--profile` so the run takes the profile's thread cap, nice and memory cap, and
refuses to start uncapped. `--cpus` and `--mem-gb` may only lower what the profile allows. Pass
`--wait`: the run then runs in the unit's main process, which keeps the unit alive until the run
closes.

```bash
bash "$TREE/scripts/pod_run.sh" --profile "$PROFILE" --image "$IMAGE" --digest "$DIGEST" --bind /tmp/if2453 --env PYTHONPATH="$PYTHONPATH_POD" w13-search search run "$WORK/recipe.yaml" --root "$ROOT" --release rC6.1 --profile "$PROFILE" --wait
```

As soon as the run is pending, the log `$BBTMP/logs/w13-search.log` has the line
`search run: <run id> running; pending record <path>`. Set `RUN` to that run id. The log ends with
the run id and its status, then `pod_run: w13-search exited 0` (`exited 1` when the status is not
`ok`). To choose the id yourself, add `--run-id` and 32 hex digits.

`search run` without `--wait` inside a `pod_run.sh` unit is refused (exit 2); the log says to add
`--wait`. The reason is what `systemd-run --collect` starts: a transient service, which ends when
its main process exits, and systemd then stops every process left in the unit's cgroup
(`KillMode=control-group`, from `systemctl --user show -p KillMode <unit>`). A detached runner
cannot outlive the command that started it: before card W-13 such a run ended
`failed: runner exited without closing the run` (measured, wilkcomp02, 2026-11-10, card C-06).
`pod_run.sh` sets `BIOBOX_POD_UNIT` to the unit's name, and that is how `bio-box` knows it is in a
unit.

**Status.** `search status` reads the run's pending record while it runs and its record once it has
closed. Run it from any shell of this account:

```bash
systemctl --user status w13-search --no-pager
RUN=77e72f990d18af9b0e4420807244e2e3
bio-box search status "$RUN" --root "$ROOT" --profile "$PROFILE"
```

The status line is JSON. While the run runs it is `{"run_id": "...", "status": "running"}`. Once it
has closed, `status` is `ok`, `failed`, `stopped`, `oom`, `timeout` or `output_cap`, with `error`,
`rounds`, `final_round` and `stop_reason`. A runner that was killed without closing its run (SIGKILL,
a reboot) reads as `failed` with `runner exited without closing the run`.

**Cancel.** From any shell of this account:

```bash
bio-box search cancel "$RUN" --root "$ROOT" --profile "$PROFILE"
```

It sends SIGTERM to the run's runner (the process its pending record names, and only while that
process holds the run's lock). The runner ends the worker and records the run as `stopped` with the
error `cancelled`. It prints nothing and exits 0, also when the run has already closed. The unit's
log then ends with `<run id>` and `stopped`, then `pod_run: <unit> exited 1`.

Stopping the unit does the same: `systemctl --user stop <unit>`. Its log then ends with
`pod_run: <unit> received SIGTERM` and `pod_run: <unit> exited 1`.

**Hits and neighbourhoods.** All of these read the release the run ran against.

```bash
bio-box search hits "$RUN" --root "$ROOT" --profile "$PROFILE"
bio-box search hits "$RUN" --root "$ROOT" --profile "$PROFILE" --format tsv --out "$WORK/hits.tsv"
bio-box search hits "$RUN" --root "$ROOT" --profile "$PROFILE" --format faa --out "$WORK/hits.faa"
bio-box nbr "$RUN" --root "$ROOT" --profile "$PROFILE" --window 5g --format tsv --out "$WORK/nbr"
bio-box nbr summarize "$RUN" --root "$ROOT" --profile "$PROFILE" --window 5g --by pfam,nr50
bio-box export proteins --root "$ROOT" --release rC6.1 --out "$WORK/proteins.faa"
```

`search hits` without `--format` prints the kept hits as a table (`--json` for JSON). `--format tsv`
is one row per gene; `--format faa` is one record per protein, named by protein uid. `nbr ...
--format tsv --out DIR` writes one file per anchor locus. `nbr summarize` counts families and
clusters around the anchors. `export proteins` needs the release's manifest; it reads the root only
and takes no `--profile`.

**Anchors from a file, and finding a gene** (card W-18). `nbr` takes a long list of protein uids from a
file, one uid per line (blank lines and `#` lines are skipped), in place of the comma-separated
argument; `find` prints the genes an accession or a protein uid names, one row each: accession,
contig uid, start, end (0-based, end exclusive), strand, feature uid, protein uid. `ACC:START-END`
names one gene, as `search hmm --positive` does (it uses the same lookup). `--accession ACC` alone lists
every gene of that contig; `--protein UID` every gene that carries the protein.

```bash
bio-box nbr summarize --anchors-file "$WORK/anchors.txt" --root "$ROOT" --profile "$PROFILE" --release rC6.1 --window 5g --by pfam,nr50
bio-box find --root "$ROOT" --release rC6.1 --accession NC_001416.1:200-1000
bio-box find --root "$ROOT" --release rC6.1 --protein 0123456789abcdef0123456789abcdef
```

**One command.** `search hmm` (card W-07) does all of the above for an HMM file: it writes the recipe,
runs it, and writes the exports into `--out DIR`. `--bits` is the inclusion threshold;
`--neighbors N` the genes on each side (window `Ng`, default 5); `--model` picks one model of a file
that holds several. `--positive` is required: a protein the search must find, given as
`ACCESSION:START-END` (one gene, with the accession, start and end that `search hits --format tsv`
prints: 0-based, end exclusive) or as a substring of its contig's header (case-sensitive). It must
name exactly one protein of the release; otherwise the command exits 2 and lists the genes that
matched, each as `ACCESSION:START-END`, so one can be given instead. `$POSITIVE` below is such a
gene. The run runs in the foreground, as `search run --wait`, so it may run in a `pod_run.sh` unit;
`search status` and `search cancel` work on it from another shell.

```bash
bash "$TREE/scripts/pod_run.sh" --profile "$PROFILE" --image "$IMAGE" --digest "$DIGEST" --bind /tmp/if2453 --env PYTHONPATH="$PYTHONPATH_POD" w07-hmm search hmm --root "$ROOT" --release rC6.1 --profile "$PROFILE" --hmm "$WORK/inputs/lib.hmm" --model lam0 --bits 25 --neighbors 5 --positive "$POSITIVE" --out "$WORK/hmm"
```

`DIR` must be new or empty and outside the data root. It then holds `recipe.yaml` (the recipe the
run ran; `search run DIR/recipe.yaml` runs it again), `hits.tsv` and `hits.faa` (as `search hits
--format tsv|faa`), `neighbors.tsv` (as `nbr summarize --format tsv`) and `loci/` (one table per
locus, as `nbr --format tsv`); each is written by that command, so the files are the same. The log
ends with the run id and `ok`, then one line per file. The run itself is in the runs directory, as
for `search run`. A release with no domain library (the 10 Gbp release `rN.1` has none) is
searched all the same: the neighbour summary is by `nr50,nr90` instead of `pfam,nr50`, and the log
says so in a line that begins `note`. A run that ends other than `ok` exits 1 and leaves only
`recipe.yaml` in `DIR`. When the release has no pack, the command exits 2 and prints the
`pack build` command to run; it does not build the pack.

## 8. Run the test suite on the pod

Run these on the laptop, in the repository checkout. The script copies the working tree to the pod
(committed or not) and runs pytest there, in the image, in a 32 GB capped scope. `<name>` is a
label, usually the task ID.

```bash
sh scripts/pod_pytest.sh C-06 tests/r0/test_runbook.py
sh scripts/pod_pytest.sh C-06 --shards 8
```

The first runs one file. The second runs every test file in 8 processes and prints one summary line
per shard; it fails if any shard fails. The whole suite takes about 13 minutes in one process and
about 3 with 8 shards (the script's header, not re-measured here). An exit code 5 with `no tests
ran` only means the selection matched nothing.

## 9. Clean up /tmp safely

What is ours, all under `/tmp/if2453/`:

- `$BBTMP/staging/<stage>.<step>/`: leftovers of killed or failed steps.
- `$BBTMP/scratch/`: builder scratch: one directory per step, `<stage>.<step>/`, which
  `build-stage` removes itself (section 4), and the scan shards' runs in `biobox-domain-scan/`,
  removed after the merge.
- `$BBTMP/logs/`: unit logs.
- `$BBTMP/packs/release=<R>/<pack id>/`: the search pack of a release (section 7). It is a cache of
  the release, shared by every user of the host: `pack build` makes it again from the release.
- `$BBTMP/search/$USER/`: your searches (section 7). `runs/run=<run id>/` is a closed run's record
  and hits; `runs/pending/run=<run id>.json` is a run's pending record (and `.stderr`, the runner's
  stderr of a run submitted outside a unit); `work/<run id>/` is a run's work directory, which the
  run removes when it closes. One is left only when the runner was killed (SIGKILL, a reboot).
- `$TREE`, `$WORK` and `/tmp/if2453/ci-tmp/<name>-*`: your mirrored trees and pytest directories.

Never delete:

- anything under `/stor`;
- anything that is not yours: another user's files, another user's `$BBTMP/search/<user>/`, and
  another agent's directories under `/tmp/if2453/` (`c10`, `b24`, `ci/<other name>`, and so on);
- a staging directory, `.incoming` directory or root file of a step that has an intent file
  (`*.publishing.json`) and no run record (section 4);
- the pending record or work directory of a search whose `bio-box search status` says `running`,
  and a pack while any search runs;
- any file of a data root: a root is append-only. If a root is wrong, tell the PI.

Before you remove anything, check that none of the units you launched is still running, and look
for intent files.

```bash
systemctl --user list-units --type=service --no-pager
find "$ROOT/runs" -name '*.publishing.json'
ls -l "$BBTMP/staging"
ls "$BBTMP/packs" "$BBTMP/search/$USER/runs" "$BBTMP/search/$USER/work"
```

Then remove leftovers one directory at a time, by exact path. Never use a glob and never `rm -rf`.

```bash
rm -r "$BBTMP/staging/S9.demo"
rm -r "$BBTMP/scratch/S9.demo"
rm -r "$BBTMP/search/$USER/work/3d96a70ba10c7f5cd895de3d838153ae"
```

Old logs: list first.

```bash
find "$BBTMP/logs" -name '*.log' -mtime +30
```

## 10. Limits

Copied from `profiles/wilkcomp02.yaml` (the profile) and `AGENTS.md` section 7 (AGENTS), each
with its source. `scripts/pod_run.sh` reads `mem_scope`, `max_threads_night` (the unit's CPU
quota), `nice_monolithic` and `tmp` from the profile; `run-step` and `build-stage` read the rest. The test `tests/r0/test_runbook.py` fails when
a profile value is missing from this table.

| Limit | Value | Source |
|---|---|---|
| Threads by day | 40 | profile `max_threads_day`; AGENTS 7.2 |
| Threads at night and weekends | 56, only when no other compute runs | profile `max_threads_night`; AGENTS 7.2 |
| Day window | 08:00 to 20:00, Monday to Friday | profile `day_start_hour` 8, `day_end_hour` 20. An estimate until the PI confirms it with the pod owner (AGENTS 7.2) |
| "Other compute" | 1-minute load average above 8 | `BUSY_LOAD` in `src/biobox/runstep.py`. An estimate, not measured |
| Memory per process or unit | 200G, swap 0 | profile `mem_scope`; AGENTS 7.2; `MemorySwapMax=0` in `scripts/pod_run.sh` |
| Memory cap required | true: `run-step` and the search runner refuse to start uncapped | profile `require_mem_cap` |
| Nice | 19 (one big process), 15 (many small ones) | profile `nice_monolithic`, `nice_fanout`; AGENTS 7.2 says 19 |
| CPU-minutes per process | 300 | profile `cpu_min_per_process`; AGENTS 7.2 |
| CPU-minutes of a release's clustering | 20000, estimated before launch; refused above it. The one exemption from the 300 (DEC-25), until card B-37. Set by the PI on 2026-11-10 for the 50 Gbp build: about nine times the 2,195.91 measured for 10,620,684 proteins (wilkcomp02, 2026-11-10). Without the key, clustering over 300 is refused | profile `cluster_cpu_min_max`; DEC-25; `cluster_refusal` in `src/biobox/stage.py`; card B-36 |
| Scan shard size | 30,150,753 residues: 80% of `cpu_min_per_process` (300) at 7.96 CPU-minutes per million residues (derived: 300 × 0.8 / 7.96 × 10^6; the rate of the Pfam scan of stage P1 in card C-10, measured on wilkcomp02, 2026-11-10: 2,343.74 CPU-minutes for 294.5 M residues). The 80% leaves room for a scan 25% slower per residue (derived: 1 / 0.8) | `SHARD_MARGIN` and `SCAN_CPU_MIN_PER_M` in `src/biobox/stage.py`; cards C-11, C-13 |
| Scan shards at once | `--scan-jobs` J × `--scan-cpus` C at most the threads `run-step` allows when `build-stage` starts (40 by day, 56 at night when idle); defaults C = 4, J = allowance // C (10 by day, 14 at night). The cap is a total over everything of ours on the pod, and `build-stage` counts only its own threads: choose J = floor((cap - other threads) / C) (section 2). The allowance is read again before each step and before each shard launch; a launch that would pass it waits | `scan_slots` and `_launch` in `src/biobox/stage.py`; cards C-11, C-12 |
| Unit CPU quota | the night cap, 56: a ceiling fixed at launch; `build-stage` holds the limit of the hour | `scripts/pod_run.sh --threads`, profile `max_threads_night`; card C-12 |
| Smallest data file | 64 MiB (67108864 bytes), except tables the registry marks `small_ok` | profile `small_file_floor_bytes`; AGENTS 7.3 says 64 MB |
| Work directory | `/tmp/if2453/bio-box` | profile `tmp`; AGENTS 7.3 |
| Search pack directory | `<tmp>/packs` (`/tmp/if2453/bio-box/packs`), outside the data root; with `--profile`, else `--pack-dir` | profile `tmp`; `default_pack_base` in `src/biobox/search/runner.py`; card W-03 |
| Search runs and work directories | `<tmp>/search/<user>/runs` and `<tmp>/search/<user>/work`, one per user, outside the data root; with `--profile`, else `--runs` and `--work` | profile `tmp`; `default_runs_dir`, `default_work_base` in `src/biobox/search/runner.py`; card W-03 |
| `search run` in a `pod_run.sh` unit | only with `--wait`; without it, refused (exit 2) | `BIOBOX_POD_UNIT` in `scripts/pod_run.sh` and `src/biobox/cli.py`; card W-13 |
| Durable results | `/stor/work/Wilke/if2453/all-the-data` | profile `durable` |
| Bulk results | `/stor/scratch/Wilke/if2453/all-the-data` | profile `bulk` |
| Results reach the data root | one `rsync`, then an atomic rename | AGENTS 7.3 |
| `STOP` file in the data root | exit cleanly before each step | AGENTS 7.4 |
| Under `/stor` | never delete or overwrite | AGENTS 7.5 |
| Docker, network services | never | AGENTS 7.5 |
| Where code comes from | `git push pod` from the PI's laptop; never edit `~/bio-box` on the pod | AGENTS 7.6 |
| How to run | only through `bio-box run-step`, in the image pulled by digest | AGENTS 7.1 |
| Which pod runs are allowed | only a card marked **Pod: yes** | AGENTS 7 |

The scratch profile used for the tests in this file differs from the real one in these keys:
`max_threads_day` 2, `max_threads_night` 2, `mem_scope` 8G, `tmp`, `durable` and `bulk` under
`/tmp/if2453/c06`, `small_file_floor_bytes` 1000, and no `cluster_cpu_min_max` (so no clustering
exemption; the fixture's clustering is far under 300 CPU-minutes). Appendix B lists it.

## 11. Known limits today

| Limit | What happens | Card |
|---|---|---|
| Clustering's DEC-25 cap | above 80% of the budget a release's clustering runs in parts (card B-41, section 2), but the DEC-25 estimate still refuses one over the profile's `cluster_cpu_min_max` (20000; about 62.7 M proteins, derived, section 2), exit 2, until the split has passed on the pod and the exemption is removed. The estimate is a fit by protein count, extrapolated above 10.6 M proteins. At 95 M proteins a prefilter part's own k-mer index of the whole target would take about 333 CPU-minutes (derived: 12.44 CPU-seconds per 59,000 targets, B-37 design §3), so no part fits the budget there; the cut then refuses | B-36, B-37, B-41 |
| Ingest of a stage over about 60 Gbp | `build-stage` refuses before any step (exit 2): one ingest process runs at about 3.7 Mbp/s (measured, wilkcomp02), so it would pass 300 CPU-minutes. Split the stage into smaller stages | B-30 |
| Shard size is Pfam's | shards are cut at 80% of the budget at Pfam's rate (7.96 CPU-minutes per million residues, measured in C-10). A library more than 1.25 times slower per residue than that (derived: 1 / 0.8; above 9.95 CPU-minutes per million residues, derived) passes 300 CPU-minutes per shard, and the run record says `over_cpu_budget`. Lower `--scan-cpus` does not help, since the budget counts CPU-minutes | C-11, C-13 |
| `pod_run.sh` CPU quota | the unit's quota is `--threads` of `pod_run.sh` (default the night cap, 56), fixed at launch, whatever `--scan-jobs` × `--scan-cpus` is. By day it is above the day cap: the day cap is held by `build-stage` and `run-step`, not by the quota. A `run-step` launched alone keeps the thread count it chose at its start until it ends (at most 300 CPU-minutes) | C-07, C-11, C-12 |
| The image on the pod is old | the pod has the test image `bio-box-test_74b3e05.sif`, not an image built from the current source. The current source is bound in with `PYTHONPATH` (section 0) | 0-02 |
| `rL.*` releases on `/stor` | they are in the pre-B-17 layout: current readers and `verify_manifest` may not read them. `rL.1` and `rL.2` also fail verification on one file, `l0/l0_manifest.parquet` | B-17, B-22 |
| Search packs and runs inside an old root | a root written before card W-03 may hold `search_pack/` and `search/` (the pack and runs were then written into the root). Nothing reads or writes them now; the manifest leaves them out (`catalog.EXCLUDED_TREES`). Build the pack again into the pack directory (section 7) | W-03, W-13 |
| No detached search in a unit | `search run` without `--wait` in a `pod_run.sh` unit is refused (exit 2): the unit would stop the detached runner when the command returns (section 7). Use `--wait`; `search status` and `search cancel` work on a `--wait` run from another shell | W-13 |
| `pod_run.sh` and the profile | `pod_run.sh --profile` sets the unit's limits only; it does not pass the profile to `bio-box` (no `BIOBOX_PROFILE` in the image). Give `bio-box` its own `--profile`, or the pack, runs and work directories (section 7) | W-03, W-13 |
| `search hmm --positive` reads all contig labels | the accession or header of the positive control is looked up in every row group of `contig_label`: no table is sorted by accession, and a header substring cannot prune. The contig's genes are then read from its tape by contig uid, which prunes. Not measured on a large release | W-07 |

Found while writing this runbook, no card yet (reported to the PI):

- The card text says `bio-box doctor`; the command is `python -m biobox.doctor`.

## 12. Contact

Contact the PI.

## Appendix A. Commands run

Run on 2026-11-10 between 07:25 and 07:45 (pod local time), in `/tmp/if2453/c06`, with units named
`c06-*`, the scratch profile of Appendix B, and the image of section 0. "Last line" is the last line
the command printed. A command that prints nothing has `(no output)`. Variables were set as in
section 0.

The rows of section 7, and the rows marked W-13, were run on 2026-11-10 between 15:35 and 15:42 (pod
local time) by card W-13, in `/tmp/if2453/w13` (`WORK`; `ROOT=$WORK/root`, `BBTMP=$WORK/tmp`,
`PROFILE=$WORK/profile.yaml`), with units named `w13-*`, the tree of the W-13 branch
(`TREE=/tmp/if2453/ci/W-13`) and the same image. Its scratch profile is Appendix B's with `tmp`,
`durable` and `bulk` under `/tmp/if2453/w13` and `max_threads_day` and `max_threads_night` 4. The
fixture of Appendix B was built there first, with the section 2 command as unit `w13-build` and
`--threads 4` (`pod_run: w13-build exited 0`). The two searches that were cancelled and stopped
used the runner's test hook `--fault timeout --wall-s 900` (its worker waits until it is stopped),
so that they ran long enough to be watched; `--fault` is not for real searches.

| Section | Command | Date run | Last line |
|---|---|---|---|
| 1 | `ls -d "$ROOT/STOP"` | 2026-11-10 | `ls: cannot access '/tmp/if2453/c06/root/STOP': No such file or directory` |
| 1 | `uptime` | 2026-11-10 | ` 07:26:38 up 11 days, 12:56,  7 users,  load average: 23.99, 23.10, 17.69` |
| 1 | `top -b -n 1 -o %CPU -w 120 > "$WORK/top.txt"` | 2026-11-10 | `(no output)` |
| 1 | `head -n 15 "$WORK/top.txt"` | 2026-11-10 | `2375504 if2453      39  19 8039140   1.2g  97688 R  18.8   0.1   1:50.74 python` |
| 1 | `systemctl --user list-units --type=service --type=scope --no-pager` | 2026-11-10 | `To show all installed unit files use 'systemctl list-unit-files'.` |
| 1 | `df -h /tmp` | 2026-11-10 | `/dev/sdc2       1.1T  165G  879G  16% /` |
| 1 | `biopy -m biobox.doctor --profile "$PROFILE"` | 2026-11-10 | `}` (the JSON shows `"cgroup_memory_max": "8 GiB"`) |
| 2 | `pod_run.sh ... c06-build build-stage ...` | 2026-11-10 | `pod_run: log /tmp/if2453/c06/tmp/logs/c06-build.log; watch: systemctl --user status c06-build` (the first launch ran without `--digest`; the same command with `--digest` ran as `c06-rebuild`, below) |
| 2 | `systemctl --user status c06-build --no-pager` | 2026-11-10 | `└─2081436 /opt/bio-box/.venv/bin/python -m biobox.cli stage inputs --root /tmp/if2453/c06/root --gene-call-version pyrodigal-3.7.1-gv-0.3.2-p2 --stage S0 --genomes /tmp/if2453/c06/inputs/genomes-S0.json` |
| 2 | `tail -n 3 "$BBTMP/logs/c06-build.log"` | 2026-11-10 | `pod_run: c06-build exited 0` |
| 2 | `grep "steps done" "$BBTMP/logs/c06-build.log"` | 2026-11-10 | `build-stage: stage S0 and release rC6.1: all 13 steps done` (run on the `c06-rebuild` log) |
| 2 | `ls "$ROOT/runs/stage=S0"` | 2026-11-10 | `tapes.json` |
| 2 | second call, `pod_run.sh ... --digest "$DIGEST" ... c06-build-again build-stage ...` (and `c06-rebuild`) | 2026-11-10 | `pod_run: c06-rebuild exited 0` (every step `skipped: run record and outputs present`) |
| 2 | second stage and release, `... c06-s1g build-stage ... --stage S1 --release rC6.2 --stages S0,S1 ...` | 2026-11-10 | `pod_run: c06-s1g exited 0` |
| 2 | `pod_run.sh ... --threads 2 c11-build build-stage ... --scan-jobs 2 --scan-cpus 1` | not run with `pod_run.sh`; tested on the fixture | `tests/r2/test_build_stage.py` (card C-11) runs `build-stage --scan-jobs 2 --scan-cpus 2` on the fixture with the scan in 7 shards (2026-11-10). No real Pfam scan was run |
| 3 | `touch "$ROOT/STOP"`, then a launch | 2026-11-10 | `pod_run: c06-s1 exited 75` (log: `STOP exists; not starting stage S1`) |
| 3 | `touch "$ROOT/STOP"` while `S1.igr` ran (unit `c06-s1d`) | 2026-11-10 | `pod_run: c06-s1d exited 75` (log: `STOP appeared; not publishing S1.igr, staging left at ...`) |
| 3 | `rm "$ROOT/STOP"` | 2026-11-10 | `(no output)` |
| 3 | `systemctl --user stop c06-s1e` while `rC6.2.cluster` ran | 2026-11-10 | `(no output)`; log: `build-stage: rC6.2.cluster exited 143; stopping here; the next call resumes at this step`, the log's last line then (before card C-13 the wrapper wrote no exit line after a stop; its trap is tested by `tests/r1/test_pod_run.py`, and on the pod by W-13, next row) |
| 3 | W-13: `systemctl --user stop w13-stop` while a search ran in it (section 7) | 2026-11-10 | `(no output)`, exit 0; the log's last two lines: `pod_run: w13-stop received SIGTERM`, then `pod_run: w13-stop exited 1` |
| 4 | relaunch after the stop (unit `c06-s1g`) | 2026-11-10 | `pod_run: c06-s1g exited 0` |
| 4 | `ls -d "$ROOT" "$BBTMP"` | 2026-11-10 | `/tmp/if2453/c06/tmp` |
| 4 | `find "$ROOT/runs" -name '*.publishing.json'` | 2026-11-10 | `(no output)` |
| 5 | record one-liner (before card C-14: now `rb.py record`, same output; see Not yet tested) | 2026-11-10 | `{'stage': 'S0', 'step': 'genecall', 'wall_s': 3.446, 'cpu_min': 0.066808, 'peak_rss_bytes': 335257600, 'threads': 2, 'exit_code': 0, 'over_cpu_budget': False} inherited 5 outputs` |
| 5 | all records, DuckDB | 2026-11-10 | `  24 rows                                                         6 columns` |
| 5 | ledger totals, DuckDB | 2026-11-10 | `│ NULL    │ rC6.2   │ cluster             │       2 │   1707 │   1684 │       23 │          0 │` (then the table's bottom border) |
| 5 | balance and chain checks | 2026-11-10 | `chain (L-03): []` |
| 6 | verify `rC6.2` | 2026-11-10 | `rC6.2 stages ['S0', 'S1'] files 41 hash 5d2522bbd240fc3c9f40c397405ad4cb problems 0` (run as a python one-liner calling `catalog.verify_manifest`, before `bio-box release verify` existed) |
| 6 | `bio-box release verify --root "$ROOT" --release rC6.2`, and with `--previous rC6.1` | not run on the pod; tested on the fixture | `tests/r2/test_build_stage.py` (card C-13): exit 0 and `release ...: verified` on a built release, exit 1 naming the file on a damaged copy, and `append-only against ...: ok` or `... rewritten` |
| 6 | open the catalog of `rC6.2` | 2026-11-10 | `contigs 3` |
| 6 | verify `rC6.1` | 2026-11-10 | `rC6.1 stages ['S0'] files 25 hash 98fbb4855e3d7e997db13336c32b7864 problems 0` (the python one-liner, as above) |
| 7 | `pod_run.sh ... w13-pack pack build ... --profile "$PROFILE"` | 2026-11-10 | `pod_run: w13-pack exited 0` (the pack in `$BBTMP/packs/release=rC6.1/16f032d9f9785d4852e2192fd63996be`; nothing under `$ROOT`) |
| 7 | `bio-box pack verify --root "$ROOT" --release rC6.1 --profile "$PROFILE"` | 2026-11-10 | `pack rC6.1: verified` |
| 7 | `pod_run.sh ... w13-search search run ... --profile "$PROFILE" --wait` | 2026-11-10 | `pod_run: w13-search exited 0` (log, before it: `search run: 77e72f990d18af9b0e4420807244e2e3 running; pending record /tmp/if2453/w13/tmp/search/if2453/runs/pending/run=77e72f990d18af9b0e4420807244e2e3.json`, then the run id and `ok`) |
| 7 | `systemctl --user status w13-cancel --no-pager` while a search ran in it | 2026-11-10 | the unit's process tree (`sh`, `Singularity runtime parent`, the runner `python -m biobox.cli search run ...`, its worker); `systemctl --user show -p Type -p KillMode w13-cancel`: `Type=simple`, `KillMode=control-group` |
| 7 | `bio-box search status "$RUN" --root "$ROOT" --profile "$PROFILE"`, from another shell while unit `w13-cancel` ran the search `131313131313131313131313131313c1` | 2026-11-10 | `{"run_id": "131313131313131313131313131313c1", "status": "running"}` |
| 7 | `bio-box search status "$RUN" --root "$ROOT" --profile "$PROFILE"` of the `ok` run | 2026-11-10 | `{"run_id": "77e72f990d18af9b0e4420807244e2e3", "status": "ok", "error": null, "rounds": 1, "final_round": 1, "stop_reason": "max_rounds"}` |
| 7 | `bio-box search cancel "$RUN" --root "$ROOT" --profile "$PROFILE"`, from another shell, of the running `131313131313131313131313131313c1` | 2026-11-10 | `(no output)`, exit 0. The unit's log then ended `131313131313131313131313131313c1	stopped` and `pod_run: w13-cancel exited 1`; `search status` then: `{"run_id": "131313131313131313131313131313c1", "status": "stopped", "error": "cancelled", "rounds": null, "final_round": null, "stop_reason": null}` |
| 7 | `systemctl --user stop w13-stop` while it ran the search `13131313131313131313131313131351`, then `search status` | 2026-11-10 | `{"run_id": "13131313131313131313131313131351", "status": "stopped", "error": "cancelled", "rounds": null, "final_round": null, "stop_reason": null}`; the log ended `pod_run: w13-stop received SIGTERM`, `pod_run: w13-stop exited 1` |
| 7 | `pod_run.sh ... w13-nowait search run ... --profile "$PROFILE"` (no `--wait`) | 2026-11-10 | `pod_run: w13-nowait exited 2`; log, before it: ``refused: search run without --wait in the pod_run.sh unit w13-nowait: the unit ends when this command returns, and systemd stops the detached runner with it. Add --wait: the run then keeps the unit until it closes, and `bio-box search status` and `search cancel` work on it from another shell`` |
| 7 | `bio-box search hits "$RUN" --root "$ROOT" --profile "$PROFILE"` | 2026-11-10 | `protein	01563e192ab22a5931d55e224a2dbb07	2598.22509765625	2598.05859375	0.0	0	1132	1	keep` |
| 7 | `bio-box search hits ... --format tsv --out "$WORK/hits.tsv"` | 2026-11-10 | `(no output)`, exit 0 |
| 7 | `bio-box search hits ... --format faa --out "$WORK/hits.faa"` | 2026-11-10 | `(no output)`, exit 0 |
| 7 | `bio-box nbr "$RUN" ... --format tsv --out "$WORK/nbr"` | 2026-11-10 | `(no output)`, exit 0; one file `locus_NC_001416.1_1afb879accb6db3171514adc98a253cb.tsv` |
| 7 | `bio-box nbr summarize "$RUN" ... --by pfam,nr50` | 2026-11-10 | `nr50		true					0	0.0	0	0	genome_uid			` |
| 7 | `bio-box export proteins ... --out "$WORK/proteins.faa"` | 2026-11-10 | `(no output)`, exit 0 |
| 7 | `bio-box search cancel "$RUN" ...` of a run that has closed | not run on the pod by W-13; tested | `tests/r2/test_micro_s0_search.py` (card W-13): exit 0, and nothing is signalled. Card C-06 ran it on the pod before W-03 (`--root` alone, runner already gone): `(no output)`, exit 0 |
| 7 | `ls -d "$ROOT/search_pack" "$ROOT/search"` after every command above | 2026-11-10 | `ls: cannot access '/tmp/if2453/w13/root/search': No such file or directory` (and the same for `search_pack`): the search wrote nothing under the root |
| 8 | `sh scripts/pod_pytest.sh C-06 tests/r0 -k nothing_matches_this` (laptop) | 2026-11-10 | `140 deselected in 1.09s` (exit 5: nothing selected) |
| 8 | `sh scripts/pod_pytest.sh C-06 tests/r0/test_runbook.py` (laptop) | 2026-11-10 | `20 passed in 3.48s` |
| 8 | `sh scripts/pod_pytest.sh C-06 --shards 8` (laptop) | 2026-11-10 | `shard 7: 121 passed, 8 xfailed in 509.90s (0:08:29)` (eight shard lines, no failure; 1023 passed, 10 skipped, 8 xfailed in all) |
| 9 | `rm -r` of a staging directory and a scratch directory | 2026-11-10 | `(no output)`, exit 0 each. The staging directory was an empty directory made for this test |
| 9 | W-13: `ls "$BBTMP/packs" "$BBTMP/search/$USER/runs" "$BBTMP/search/$USER/work"` | 2026-11-10 | `/tmp/if2453/w13/tmp/search/if2453/work:` (empty: every run had removed its work directory; `runs` held `pending` and the three runs of section 7, `packs` held `release=rC6.1`) |
| 9 | W-13: `rm -r "$BBTMP/search/$USER/work/3d96a70ba10c7f5cd895de3d838153ae"` | 2026-11-10 | `(no output)`, exit 0. An empty directory made for this test |
| 9 | `find "$BBTMP/logs" -name '*.log' -mtime +30` | 2026-11-10 | `(no output)` |

## Appendix B. Rehearse on the fixture

The scratch profile `$WORK/profile.yaml`. It differs from the real profile only in the values that
section 10 lists.

```yaml
max_threads_day: 2
max_threads_night: 2
mem_scope: 8G
tmp: /tmp/if2453/c06/tmp
durable: /tmp/if2453/c06/durable
bulk: /tmp/if2453/c06/bulk
nice_monolithic: 19
nice_fanout: 15
day_start_hour: 8
day_end_hour: 20
cpu_min_per_process: 300
small_file_floor_bytes: 1000
require_mem_cap: true
```

Helper script `$WORK/rb.py`, used by sections 5 and 6 and below. Save it with an editor or `scp`; the
commands that use it then have no quotes inside quotes.

```python
"""Read run records, the ledger and the catalog: rb.py record|records|ledger|checks|catalog|recipe ARGS."""
import json
import sys
from pathlib import Path

from biobox import catalog, ledger


def open_cat(root, release):
    return catalog.open_catalog(Path(root) / "catalog" / f"release={release}.duckdb", threads=2)


def record(path):
    r = json.load(open(path))
    keys = ("stage", "step", "wall_s", "cpu_min", "peak_rss_bytes", "threads", "exit_code", "over_cpu_budget")
    print({k: r[k] for k in keys}, r["mem_cap"]["how"], len(r["outputs"]), "outputs")


def records(root):
    import duckdb
    con = duckdb.connect(config={"threads": 2})
    con.sql("SELECT stage, step, wall_s, cpu_min, peak_rss_bytes, mem_cap.how AS cap "
            "FROM read_json_auto(?, filename = true, union_by_name = true) "
            "WHERE filename NOT LIKE ? ORDER BY started_at",
            params=[root + "/runs/*/*.json", "%.publishing.json"]).show(max_width=200, max_rows=60)


def totals(root, release):
    open_cat(root, release).sql(
        "SELECT stage, release, step, count(*) AS objects, sum(n_in) AS n_in, sum(n_kept) AS n_kept, "
        "sum(n_merged) AS n_merged, sum(n_rejected) AS n_rejected FROM (SELECT stage, release, step, object_key, "
        "min(n_in) AS n_in, min(n_kept) AS n_kept, min(n_merged) AS n_merged, sum(n_rejected) AS n_rejected "
        "FROM ledger GROUP BY ALL) GROUP BY ALL ORDER BY ALL").show(max_width=250)


def checks(root, stage, release, steps):
    print("balance (L-01):", ledger.check_balance(Path(root), stage, steps.split(","), release))
    print("chain (L-03):", ledger.check_chain(Path(root), stage, release))


def info(root, release):
    con = open_cat(root, release)
    print(con.execute("SELECT key, value FROM catalog_info ORDER BY key").fetchall())
    print("contigs", con.execute("SELECT count(*) FROM contig_index").fetchone()[0])


def recipe(root, out):
    con = open_cat(root, "rC6.1")
    uid = con.execute("SELECT lower(hex(protein_uid)) FROM l1_domains ORDER BY 1 LIMIT 1").fetchone()[0]
    Path(out).write_text(f"recipe: fixture_phmmer\nwritten_against: rC6.1\nseed: [{{protein_uid: {uid}}}]\n"
                         f"target: unique\ntool: phmmer\niterations: {{max: 1, inclusion_bits: 25}}\n"
                         f"controls: {{positive: [{{protein_uid: {uid}}}]}}\n")
    print(uid)


COMMANDS = {"record": record, "records": records, "ledger": totals, "checks": checks, "catalog": info,
            "recipe": recipe}
COMMANDS[sys.argv[1]](*sys.argv[2:])
```

Fixture inputs: the genome and object lists for stages `S0` and `S1` and an HMM library, made with
the helpers of `tests/r2/test_builders_out.py`. Save this as `$WORK/prep.py`.

```python
"""Make the fixture inputs: genomes-S0.json, objects-S0.json and lib.hmm, in $WORK/inputs."""
import shutil
import sys
from pathlib import Path

tree, work = Path(sys.argv[1]), Path(sys.argv[2])
sys.path.insert(0, str(tree / "tests" / "r2"))
import test_builders_out as c08  # the helpers the tests use

w = c08.Work(work / "inputs")
w.prepare()
scratch_root = work / "libroot"
for step in c08.stage_steps("S0", w)[:3]:  # inputs, ingest, genecall: the proteins the HMMs are built from
    step.direct(scratch_root)
c08.write_library(scratch_root, w.lib)
shutil.rmtree(scratch_root)
print("inputs written under", w.base)
```

```bash
biopy "$WORK/prep.py" "$TREE" "$WORK"
```

Last line: `inputs written under /tmp/if2453/c06/inputs`. After the build of `rC6.1`, this writes a
fixture search recipe, `$WORK/recipe.yaml`, that searches for one library-hit protein and takes it
as its own positive control. It prints the protein uid.

```bash
biopy "$WORK/rb.py" recipe "$ROOT" "$WORK/recipe.yaml"
```

Last line: `01563e192ab22a5931d55e224a2dbb07`.

## Not yet tested

- `$WORK/rb.py` (Appendix B), which replaced the `python -c` one-liners of sections 5 and 6 and Appendix B in card C-14. It holds the same code, but was not run on the pod: C-14 may run only the test suite there. Run each of its commands once, and through `ssh`, and correct this line.
- Recovery from an intent file (a kill during publishing). It cannot be hit by hand at fixture
  size; `tests/r2/test_runstep_resume.py` kills a step at each publish point with an injected hook.
- A reboot. A kill by `systemctl --user stop` was tested; the machine was not rebooted.
- A hard kill of a unit (SIGKILL). `systemctl --user kill --signal=SIGKILL` printed
  `Failed to send signal SIGKILL to auxiliary processes: Invalid argument`, though the unit ended;
  it is not documented here.
- The night thread cap. The pod was busy (load above 8) at every run, so `run-step` used the day
  cap.
- `sh scripts/pod_pytest.sh <name> --golden` (the golden ring needs the golden set on the pod).
- `git push pod`. Agents never push.
- Single release steps by hand (`bio-box release cluster ...`, `bio-box stage genecall ...`).
  Only `build-stage` ran them.
- Clustering in parts (card B-41, section 2) in a build unit, and on a real release: it was tested in
  the test suite on the pod (the fixture and a constructed set, through `build-stage` on the fixture),
  not in a build launched with `pod_run.sh`. The reclustering of `rN.1b` in parts (2.05 M proteins,
  the card's acceptance) and the removal of the DEC-25 exemption have not been done.
- A build at the real profile's 64 MiB floor and 200G cap, and a real Pfam scan (card C-10).
- Inputs outside `/tmp/if2453` (they need an extra `--bind`, section 2), and `--digest` refusing a
  wrong image (the script's `die` path was read, not run).
- Card C-13's changes in a real build unit: the timestamped log lines, the per-step scratch
  directories and `bio-box release verify` were tested in the test suite on the pod
  (`scripts/pod_pytest.sh`), not in a build launched with `pod_run.sh`. The exit line after
  `systemctl --user stop` was seen in a real unit by card W-13, in a search unit (section 3), not in
  a build unit.
- `search run` without `--wait` outside a `pod_run.sh` unit (a submitted, detached runner, for
  example from the `bio-box` function of section 0). It is tested in the test suite, not on the pod.
- `search hmm` (section 7) in a `pod_run.sh` unit, and on a real release. Card W-07 tested it in the
  test suite on the pod, on the fixture: its files equal those of the separate commands byte for
  byte, its recipe run again through `search run` gives the same hits, a fixture release without a
  domain library gives the summary by `nr50,nr90`, and its refusals (no positive control, one that
  matches no protein or several, no pack, `DIR` under the root or not empty) exit 2.
