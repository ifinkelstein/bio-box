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
- Values like `S0`, `S1`, `rC6.1`, `rC6.2` and `c06-build` are the fixture names used when the
  commands were tested. Replace them with your own stage, release and unit names.

The commands were tested on a small fixture build in a scratch data root with a scratch copy of
the host profile. Where the real profile differs, the text says so. The real profile is
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

Each step prints one JSON line to the log as it ends. `run-step` and `build-stage` lines come at
the end of the log, because the unit's output is buffered until it exits.

How to tell it finished:

```bash
tail -n 3 "$BBTMP/logs/c06-build.log"
grep "all 13 steps done" "$BBTMP/logs/c06-build.log"
ls "$ROOT/runs/stage=S0"
```

- The last line is `pod_run: c06-build exited 0`. The number is the exit code.
- The log has `build-stage: stage S0 and release rC6.1: all 13 steps done`, and a
  `run-step: ... published N files` line for each of the 12 steps that run a command (the
  `features` slot is empty).
- `systemctl --user status` says `could not be found`: the unit has ended.
- A log whose last line is not `pod_run: ... exited`: the unit was stopped from outside
  (`systemctl --user stop`, or a reboot). Go to section 4.

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
- `--scan-jobs J`: shards at once. Default: the allowance divided by C, rounded down.
- The allowance is the thread count `run-step` allows when `build-stage` starts: 40 by day; 56 at
  night and at weekends, only while the load average is at most 8.
- J × C above the allowance is refused before any step (exit 2). The message names J, C and the
  allowance.
- The output is the same, byte for byte, for any J and C.

The unit's CPU quota is `pod_run.sh --threads`, which defaults to the day cap (40). For a night build
at the night cap, pass `--threads 56` to `pod_run.sh` too. This fixture launch runs two shards of one
thread at once; the scratch profile allows 2 threads:

```bash
bash "$TREE/scripts/pod_run.sh" --profile "$PROFILE" --image "$IMAGE" --digest "$DIGEST" --bind /tmp/if2453 --env PYTHONPATH="$PYTHONPATH_POD" --threads 2 c11-build build-stage --profile "$PROFILE" --root "$ROOT" --stage S0 --release rC6.1 --stages S0 --genomes "$WORK/inputs/genomes-S0.json" --objects "$WORK/inputs/objects-S0.json" --library test 1 "$WORK/inputs/lib.hmm" gathering --scan-jobs 2 --scan-cpus 1
```

For a real Pfam scan by day, with the real profile: `--scan-jobs 10 --scan-cpus 4`.

A shard that fails stops new launches. Shards already running finish and publish, then
`build-stage` exits with the failed shard's code. The next call runs only the shards that have no
run record, then the merge. While shards run, the log holds one `run-step:` line per shard as it
ends, in the order they end.

## 3. Pause with STOP, and resume

Pause: create the file in the data root.

```bash
touch "$ROOT/STOP"
```

- A new launch exits 75 and runs nothing. The log says `STOP exists; not starting`.
- A running build checks before each step and exits 75 at the next step boundary. During a scan
  in shards it checks before each shard launch, launches no more, and waits for the shards that run.
- A step that is running when the file appears finishes its command. `run-step` then publishes
  nothing, leaves its staging directory and exits 75. The log says `STOP appeared; not publishing`.
  The next call runs that step again from the start.

To stop at once, without waiting for the step, stop the unit. `run-step` forwards the SIGTERM to
the step and publishes nothing from it.

```bash
systemctl --user stop c06-build
```

A stopped unit leaves no `pod_run: ... exited` line in its log.

Resume: remove the file, then launch the same build command again (section 2), with a new unit name.

```bash
rm "$ROOT/STOP"
```

Steps with a run record are skipped. The step that stopped runs again.

## 4. After a kill or a reboot

Relaunch the same command from section 2 with a new unit name. That is the whole recovery.

First check that the data root and `tmp` are still there. After a reboot `/tmp` may have been
cleared. If `$ROOT` was under `/tmp` and is gone, the build is gone with it.

```bash
ls -d "$ROOT" "$BBTMP"
```

What `run-step` leaves, and what it does at the next start:

| Left behind | When | What `run-step` does |
|---|---|---|
| `$BBTMP/staging/<stage>.<step>/` | the step was killed, failed or stopped | removes it and runs the step again |
| `$BBTMP/scratch/...` | a builder was killed | nothing; a builder makes a new scratch directory |
| `$ROOT/runs/stage=<S>/<step>.json` (run record) | the step finished | skips the step if every listed output is in the root with its recorded size |
| `$ROOT/runs/stage=<S>/<step>.publishing.json` (intent file) | killed while publishing | finishes the publish, or undoes it, then reruns (see below) |
| `$ROOT/.incoming/<stage>.<step>/` | killed while publishing | uses it to finish the publish, then removes it |
| `$ROOT/.incoming/_build-stage.<stage>.<scan step>/` (empty) | `build-stage` was killed while scan shards ran | nothing; `build-stage` removes it when the shards of its next call end |

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
| `image_digest`, `host` | the image sha256 that `pod_run.sh --digest` passed in, and the host |
| `outputs` | each published file: `path`, `bytes`, `sha256` |

One record:

```bash
biopy -c "import json, sys; r = json.load(open(sys.argv[1])); print({k: r[k] for k in ('stage', 'step', 'wall_s', 'cpu_min', 'peak_rss_bytes', 'threads', 'exit_code', 'over_cpu_budget')}, r['mem_cap']['how'], len(r['outputs']), 'outputs')" "$ROOT/runs/stage=S0/genecall.json"
```

Every record in the root, in start order:

```bash
biopy -c 'import duckdb, sys; con = duckdb.connect(config={"threads": 2}); con.sql("SELECT stage, step, wall_s, cpu_min, peak_rss_bytes, mem_cap.how AS cap FROM read_json_auto(?, filename = true, union_by_name = true) WHERE filename NOT LIKE ? ORDER BY started_at", params=[sys.argv[1] + "/runs/*/*.json", "%.publishing.json"]).show(max_width=200, max_rows=60)' "$ROOT"
```

The ledger is a Parquet table, one file per step: `ledger/stage=<S>/step=<name>/part-0.parquet`,
and for clustering `ledger/release=<R>/step=cluster/part-0.parquet`. One row per (object, reason):
`object_key`, `object_kind`, `n_in`, `n_kept`, `n_merged`, `n_rejected`, `reason_code`,
`tool_versions_json` (the versions of the tools that shaped the output). `n_in`, `n_kept` and
`n_merged` repeat on every reason row of an object. `n_rejected` counts that row's reason only.
A release's catalog has the view `ledger`. Totals per step, with each object counted once:

```bash
biopy -c 'import sys; from pathlib import Path; from biobox import catalog; con = catalog.open_catalog(Path(sys.argv[1]) / "catalog" / f"release={sys.argv[2]}.duckdb", threads=2); con.sql("SELECT stage, release, step, count(*) AS objects, sum(n_in) AS n_in, sum(n_kept) AS n_kept, sum(n_merged) AS n_merged, sum(n_rejected) AS n_rejected FROM (SELECT stage, release, step, object_key, min(n_in) AS n_in, min(n_kept) AS n_kept, min(n_merged) AS n_merged, sum(n_rejected) AS n_rejected FROM ledger GROUP BY ALL) GROUP BY ALL ORDER BY ALL").show(max_width=250)' "$ROOT" rC6.2
```

Two checks. Balance (L-01): per object, `n_in = n_kept + n_merged + sum(n_rejected)`. Chain (L-03):
along each edge, the downstream `n_in` equals the upstream `n_kept`. Each returns a list of
problems; an empty list is a pass. The chain check needs a release whose stages are one stage.

```bash
biopy -c 'import sys; from pathlib import Path; from biobox import ledger; root = Path(sys.argv[1]); print("balance (L-01):", ledger.check_balance(root, "S0", ["ingest", "genecall", "domains-test-1", "cluster"], "rC6.1")); print("chain (L-03):", ledger.check_chain(root, "S0", "rC6.1"))' "$ROOT"
```

## 6. Verify a release

Manifest verification re-reads every listed file and compares size and sha256. It also reports any
unlisted file in the release's scope. `problems 0` is a pass.

```bash
biopy -c 'import sys; from pathlib import Path; from biobox import catalog; m = catalog.read_manifest(catalog.manifest_path(Path(sys.argv[1]), sys.argv[2])); problems = catalog.verify_manifest(m); print(m["release_id"], "stages", m["stages"], "files", len(m["files"]), "hash", catalog.manifest_hash(m), "problems", len(problems)); print(*problems, sep="\n")' "$ROOT" rC6.2
```

Open the catalog and read from it. The catalog is a DuckDB file that names the root it was built
from; it reads the root's Parquet files by the manifest's list.

```bash
biopy -c 'import sys; from pathlib import Path; from biobox import catalog; con = catalog.open_catalog(Path(sys.argv[1]) / "catalog" / f"release={sys.argv[2]}.duckdb", threads=2); print(con.execute("SELECT key, value FROM catalog_info ORDER BY key").fetchall()); print("contigs", con.execute("SELECT count(*) FROM contig_index").fetchone()[0])' "$ROOT" rC6.2
```

An earlier release must still verify after a later one is built. Check `rC6.1` the same way:

```bash
biopy -c 'import sys; from pathlib import Path; from biobox import catalog; m = catalog.read_manifest(catalog.manifest_path(Path(sys.argv[1]), sys.argv[2])); problems = catalog.verify_manifest(m); print(m["release_id"], "stages", m["stages"], "files", len(m["files"]), "hash", catalog.manifest_hash(m), "problems", len(problems)); print(*problems, sep="\n")' "$ROOT" rC6.1
```

## 7. Search

Search commands that run for more than a moment go through `scripts/pod_run.sh`, so they run in a
capped unit. Always give `--work` a directory under your own scratch: the default is
`/tmp/atd/search`, which is not ours.

**Pack.** Build the search pack of a release once. `--model` is `<library>-<version>=<HMM file>`
and must match a domain library of the release.

```bash
bash "$TREE/scripts/pod_run.sh" --profile "$PROFILE" --image "$IMAGE" --bind /tmp/if2453 --env PYTHONPATH="$PYTHONPATH_POD" c06-pack pack build --root "$ROOT" --release rC6.1 --model test-1="$WORK/inputs/lib.hmm"
bio-box pack verify --root "$ROOT" --release rC6.1
```

**Recipe.** A YAML file; its format is in the docstring of `src/biobox/search/recipe.py`. A recipe
without a positive control is refused. Appendix B writes a fixture recipe.

**Run.** Pass `--profile` so the run takes the profile's thread cap, nice and memory cap, and
refuses to start uncapped. `--cpus` and `--mem-gb` may only lower what the profile allows. Pass
`--wait`: the runner then stays in the unit until the run ends. Without `--wait` the runner is a
detached child that dies when the unit's main process exits, and the run is recorded as `failed`
(section 11).

```bash
bash "$TREE/scripts/pod_run.sh" --profile "$PROFILE" --image "$IMAGE" --bind /tmp/if2453 --env PYTHONPATH="$PYTHONPATH_POD" c06-search search run "$WORK/recipe.yaml" --root "$ROOT" --release rC6.1 --work "$WORK/search-work" --profile "$PROFILE" --wait
```

The log `$BBTMP/logs/c06-search.log` ends with the run id and its status, then
`pod_run: c06-search exited 0`. Set `RUN` to that run id. To know the id before the run starts, add
`--run-id` and a 32-digit hex number of your own.

**Status.** While the unit runs, use the unit and the log. `search status` finds a run only when it
has ended.

```bash
systemctl --user status c06-search --no-pager
RUN=6592f8dad89ce024c8140f2228c1ab24
bio-box search status "$RUN" --root "$ROOT"
```

The status line is JSON: `status` is `ok`, `failed`, `stopped`, `oom`, `timeout` or `output_cap`.

**Cancel.** Stop the unit. The runner takes the SIGTERM, ends the worker and records the run as
`stopped` with the error `cancelled`.

```bash
systemctl --user stop c06-search3
bio-box search status a1b2c3d4e5f60718293a4b5c6d7e8f90 --root "$ROOT"
```

`bio-box search cancel <run id> --root R` signals the runner of a run started without `--wait`. It
prints nothing and exits 0, also when the runner is already gone. Under `pod_run.sh` there is no
such runner, so use `systemctl --user stop`.

**Hits and neighbourhoods.** All of these read the release the run ran against.

```bash
RUN=6592f8dad89ce024c8140f2228c1ab24
bio-box search hits "$RUN" --root "$ROOT"
bio-box search hits "$RUN" --root "$ROOT" --format tsv --out "$WORK/hits.tsv"
bio-box search hits "$RUN" --root "$ROOT" --format faa --out "$WORK/hits.faa"
bio-box nbr "$RUN" --root "$ROOT" --window 5g --format tsv --out "$WORK/nbr"
bio-box nbr summarize "$RUN" --root "$ROOT" --window 5g --by pfam,nr50
bio-box export proteins --root "$ROOT" --release rC6.1 --out "$WORK/proteins.faa"
```

`search hits` without `--format` prints the kept hits as a table (`--json` for JSON). `--format tsv`
is one row per gene; `--format faa` is one record per protein, named by protein uid. `nbr ...
--format tsv --out DIR` writes one file per anchor locus. `nbr summarize` counts families and
clusters around the anchors. `export proteins` needs the release's manifest.

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
- `$BBTMP/scratch/`: builder scratch.
- `$BBTMP/logs/`: unit logs.
- `$WORK/search-work/<run id>/`: a search run's work directory.
- `$TREE`, `$WORK` and `/tmp/if2453/ci-tmp/<name>-*`: your mirrored trees and pytest directories.

Never delete:

- anything under `/stor`;
- anything that is not yours: another user's files, and another agent's directories under
  `/tmp/if2453/` (`c10`, `b24`, `ci/<other name>`, and so on);
- a staging directory, `.incoming` directory or root file of a step that has an intent file
  (`*.publishing.json`) and no run record (section 4);
- any file of a data root: a root is append-only. If a root is wrong, tell the PI.

Before you remove anything, check that none of the units you launched is still running, and look
for intent files.

```bash
systemctl --user list-units --type=service --no-pager
find "$ROOT/runs" -name '*.publishing.json'
ls -l "$BBTMP/staging"
```

Then remove leftovers one directory at a time, by exact path. Never use a glob and never `rm -rf`.

```bash
rm -r "$BBTMP/staging/S9.demo"
rm -r "$BBTMP/scratch/tmpflnpepuv"
rm -r "$WORK/search-work/3d96a70ba10c7f5cd895de3d838153ae"
```

Old logs: list first.

```bash
find "$BBTMP/logs" -name '*.log' -mtime +30
```

## 10. Limits

Copied from `profiles/wilkcomp02.yaml` (the profile) and `AGENTS.md` section 7 (AGENTS), each
with its source. `scripts/pod_run.sh` reads `mem_scope`, `max_threads_day`, `nice_monolithic` and
`tmp` from the profile; `run-step` reads the rest. The test `tests/r0/test_runbook.py` fails when
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
| Scan shard size | 32,876,712 residues: 80% of `cpu_min_per_process` (300) at 7.30 CPU-minutes per million residues (derived: 300 × 0.8 / 7.30 × 10^6; the rate measured for Pfam on wilkcomp02, 2026-11-10) | `SHARD_MARGIN` and `SCAN_CPU_MIN_PER_M` in `src/biobox/stage.py`; card C-11 |
| Scan shards at once | `--scan-jobs` J × `--scan-cpus` C at most the threads `run-step` allows when `build-stage` starts (40 by day, 56 at night when idle); defaults C = 4, J = allowance // C (10 by day, 14 at night) | `scan_slots` in `src/biobox/stage.py`; card C-11 |
| Smallest data file | 64 MiB (67108864 bytes), except tables the registry marks `small_ok` | profile `small_file_floor_bytes`; AGENTS 7.3 says 64 MB |
| Work directory | `/tmp/if2453/bio-box` | profile `tmp`; AGENTS 7.3 |
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
`/tmp/if2453/c06`, `small_file_floor_bytes` 1000. Appendix B lists it.

## 11. Known limits today

| Limit | What happens | Card |
|---|---|---|
| Ingest of a stage over about 60 Gbp | `build-stage` refuses before any step (exit 2): one ingest process runs at about 3.7 Mbp/s (measured, wilkcomp02), so it would pass 300 CPU-minutes. Split the stage into smaller stages | B-30 |
| The scan allowance is decided once | `build-stage` checks `--scan-jobs` × `--scan-cpus` against the allowance when it starts. A build started at night with 56 threads of shards keeps 56 after 08:00 on a weekday. If a build will run into the day, start it with the day allowance (`--scan-jobs 10 --scan-cpus 4`) | C-11 |
| Shard size is Pfam's | shards are cut at 80% of the budget at Pfam's rate. A library more than 1.25 times slower per residue than Pfam (derived: 1 / 0.8) passes 300 CPU-minutes per shard, and the run record says `over_cpu_budget`. Lower `--scan-cpus` does not help, since the budget counts CPU-minutes | C-11 |
| `pod_run.sh` CPU quota | the unit's quota is `--threads` of `pod_run.sh` (default 40), whatever `--scan-jobs` × `--scan-cpus` is. At night pass `--threads 56` to use the night cap | C-07, C-11 |
| The image on the pod is old | the pod has the test image `bio-box-test_74b3e05.sif`, not an image built from the current source. The current source is bound in with `PYTHONPATH` (section 0) | 0-02 |
| `rL.*` releases on `/stor` | they are in the pre-B-17 layout: current readers and `verify_manifest` may not read them. `rL.1` and `rL.2` also fail verification on one file, `l0/l0_manifest.parquet` | B-17, B-22 |

Found while writing this runbook, no card yet (reported to the PI):

- `bio-box search run` without `--wait`, started with `scripts/pod_run.sh`, records the run as
  `failed` with `runner exited without closing the run`, because the unit ends with its main
  process and takes the detached runner with it. Use `--wait`.
- `search status` of a run started with `--wait` finds nothing until the run has ended.
- After `systemctl --user stop`, the unit's log has no `pod_run: ... exited` line.
- The card text says `bio-box doctor`; the command is `python -m biobox.doctor`.
- A scan's shard runs live in scratch (`$BBTMP/scratch/biobox-domain-scan/`) until the merge. If
  scratch is cleared after some shards have published (a reboot clears `/tmp`), the merge refuses,
  because those runs are gone, and their shards are not scanned again, because their run records
  exist. Found by card C-11, from the code of `domains.merge` and `domains.scan_shard`; not run.

## 12. Contact

Contact the PI.

## Appendix A. Commands run

Run on 2026-11-10 between 07:25 and 07:45 (pod local time), in `/tmp/if2453/c06`, with units named
`c06-*`, the scratch profile of Appendix B, and the image of section 0. "Last line" is the last line
the command printed. A command that prints nothing has `(no output)`. Variables were set as in
section 0.

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
| 2 | `grep "all 13 steps done" "$BBTMP/logs/c06-build.log"` | 2026-11-10 | `build-stage: stage S0 and release rC6.1: all 13 steps done` (run on the `c06-rebuild` log) |
| 2 | `ls "$ROOT/runs/stage=S0"` | 2026-11-10 | `tapes.json` |
| 2 | second call, `pod_run.sh ... --digest "$DIGEST" ... c06-build-again build-stage ...` (and `c06-rebuild`) | 2026-11-10 | `pod_run: c06-rebuild exited 0` (every step `skipped: run record and outputs present`) |
| 2 | second stage and release, `... c06-s1g build-stage ... --stage S1 --release rC6.2 --stages S0,S1 ...` | 2026-11-10 | `pod_run: c06-s1g exited 0` |
| 2 | `pod_run.sh ... --threads 2 c11-build build-stage ... --scan-jobs 2 --scan-cpus 1` | not run with `pod_run.sh`; tested on the fixture | `tests/r2/test_build_stage.py` (card C-11) runs `build-stage --scan-jobs 2 --scan-cpus 2` on the fixture with the scan in 7 shards (2026-11-10). No real Pfam scan was run |
| 3 | `touch "$ROOT/STOP"`, then a launch | 2026-11-10 | `pod_run: c06-s1 exited 75` (log: `STOP exists; not starting stage S1`) |
| 3 | `touch "$ROOT/STOP"` while `S1.igr` ran (unit `c06-s1d`) | 2026-11-10 | `pod_run: c06-s1d exited 75` (log: `STOP appeared; not publishing S1.igr, staging left at ...`) |
| 3 | `rm "$ROOT/STOP"` | 2026-11-10 | `(no output)` |
| 3 | `systemctl --user stop c06-s1e` while `rC6.2.cluster` ran | 2026-11-10 | `(no output)`; log: `build-stage: rC6.2.cluster exited 143; stopping here; the next call resumes at this step` |
| 4 | relaunch after the stop (unit `c06-s1g`) | 2026-11-10 | `pod_run: c06-s1g exited 0` |
| 4 | `ls -d "$ROOT" "$BBTMP"` | 2026-11-10 | `/tmp/if2453/c06/tmp` |
| 4 | `find "$ROOT/runs" -name '*.publishing.json'` | 2026-11-10 | `(no output)` |
| 5 | record one-liner | 2026-11-10 | `{'stage': 'S0', 'step': 'genecall', 'wall_s': 3.446, 'cpu_min': 0.066808, 'peak_rss_bytes': 335257600, 'threads': 2, 'exit_code': 0, 'over_cpu_budget': False} inherited 5 outputs` |
| 5 | all records, DuckDB | 2026-11-10 | `  24 rows                                                         6 columns` |
| 5 | ledger totals, DuckDB | 2026-11-10 | `│ NULL    │ rC6.2   │ cluster             │       2 │   1707 │   1684 │       23 │          0 │` (then the table's bottom border) |
| 5 | balance and chain checks | 2026-11-10 | `chain (L-03): []` |
| 6 | verify `rC6.2` | 2026-11-10 | `rC6.2 stages ['S0', 'S1'] files 41 hash 5d2522bbd240fc3c9f40c397405ad4cb problems 0` |
| 6 | open the catalog of `rC6.2` | 2026-11-10 | `contigs 3` |
| 6 | verify `rC6.1` | 2026-11-10 | `rC6.1 stages ['S0'] files 25 hash 98fbb4855e3d7e997db13336c32b7864 problems 0` |
| 7 | `pod_run.sh ... c06-pack pack build ...` | 2026-11-10 | `pod_run: c06-pack exited 0` |
| 7 | `bio-box pack verify --root "$ROOT" --release rC6.1` | 2026-11-10 | `pack rC6.1: verified` |
| 7 | `pod_run.sh ... c06-search search run ... --wait` | 2026-11-10 | `pod_run: c06-search exited 0` (log: run id and `ok`) |
| 7 | `bio-box search status "$RUN" --root "$ROOT"` | 2026-11-10 | `{"run_id": "6592f8dad89ce024c8140f2228c1ab24", "status": "ok", "error": null, "rounds": 1, "final_round": 1, "stop_reason": "max_rounds"}` |
| 7 | `systemctl --user stop c06-search3` (run with `--run-id a1b2c3d4e5f60718293a4b5c6d7e8f90`), then `search status` | 2026-11-10 | `{"run_id": "a1b2c3d4e5f60718293a4b5c6d7e8f90", "status": "stopped", "error": "cancelled", "rounds": null, "final_round": null, "stop_reason": null}` |
| 7 | `bio-box search hits "$RUN" --root "$ROOT"` | 2026-11-10 | `protein	01563e192ab22a5931d55e224a2dbb07	2598.22509765625	2598.05859375	0.0	0	1132	1	keep` |
| 7 | `bio-box search hits ... --format tsv --out "$WORK/hits.tsv"` | 2026-11-10 | `(no output)`, exit 0 |
| 7 | `bio-box search hits ... --format faa --out "$WORK/hits.faa"` | 2026-11-10 | `(no output)`, exit 0 |
| 7 | `bio-box nbr "$RUN" ... --format tsv --out "$WORK/nbr"` | 2026-11-10 | `(no output)`, exit 0; one file `locus_NC_001416.1_1afb879accb6db3171514adc98a253cb.tsv` |
| 7 | `bio-box nbr summarize "$RUN" ... --by pfam,nr50` | 2026-11-10 | `nr50		true					0	0.0	0	0	genome_uid			` |
| 7 | `bio-box export proteins ... --out "$WORK/proteins.faa"` | 2026-11-10 | `(no output)`, exit 0 |
| 7 | `bio-box search cancel 3d96a70ba10c7f5cd895de3d838153ae --root "$ROOT"` (runner already gone) | 2026-11-10 | `(no output)`, exit 0 |
| 8 | `sh scripts/pod_pytest.sh C-06 tests/r0 -k nothing_matches_this` (laptop) | 2026-11-10 | `140 deselected in 1.09s` (exit 5: nothing selected) |
| 8 | `sh scripts/pod_pytest.sh C-06 tests/r0/test_runbook.py` (laptop) | 2026-11-10 | `20 passed in 3.48s` |
| 8 | `sh scripts/pod_pytest.sh C-06 --shards 8` (laptop) | 2026-11-10 | `shard 7: 121 passed, 8 xfailed in 509.90s (0:08:29)` (eight shard lines, no failure; 1023 passed, 10 skipped, 8 xfailed in all) |
| 9 | `rm -r` of a staging directory, a scratch directory and a search work directory | 2026-11-10 | `(no output)`, exit 0 each. The staging directory was an empty directory made for this test |
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
biopy -c 'import sys; from pathlib import Path; from biobox import catalog; con = catalog.open_catalog(Path(sys.argv[1]) / "catalog" / "release=rC6.1.duckdb", threads=2); uid = con.execute("SELECT lower(hex(protein_uid)) FROM l1_domains ORDER BY 1 LIMIT 1").fetchone()[0]; Path(sys.argv[2]).write_text(f"recipe: fixture_phmmer\nwritten_against: rC6.1\nseed: [{{protein_uid: {uid}}}]\ntarget: unique\ntool: phmmer\niterations: {{max: 1, inclusion_bits: 25}}\ncontrols: {{positive: [{{protein_uid: {uid}}}]}}\n"); print(uid)' "$ROOT" "$WORK/recipe.yaml"
```

Last line: `01563e192ab22a5931d55e224a2dbb07`.

## Not yet tested

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
- A build at the real profile's 64 MiB floor and 200G cap, and a real Pfam scan (card C-10).
- Inputs outside `/tmp/if2453` (they need an extra `--bind`, section 2), and `--digest` refusing a
  wrong image (the script's `die` path was read, not run).
