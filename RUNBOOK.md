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
export IMAGE="$HOME/bio-box/images/bio-box_64f1bf0.sif"
export DIGEST=sha256:cb7353f85316981be5268c2212100322b5fbf1086c32b07e2f9b47802b6afba3
export PYTHONPATH_POD="$TREE/src"                    # the tree is newer than the image: its source is bound in (section 11)
export PROFILE=/tmp/if2453/c06/profile.yaml          # real: "$TREE/profiles/wilkcomp02.yaml"
export MEM=8G                                        # real: 200G (the profile's mem_scope)
export WORK=/tmp/if2453/c06                          # your scratch directory for inputs and outputs
export BBTMP=/tmp/if2453/c06/tmp                     # the profile's tmp (real: /tmp/if2453/bio-box)
export ROOT=/tmp/if2453/c06/root                     # the data root you build into
export DURABLE=/tmp/if2453/w28/one/durable           # the profile's durable (real: /stor/work/Wilke/if2453/all-the-data)
SCOPE=(systemd-run --user --scope --quiet -p MemoryMax="$MEM" -p MemorySwapMax=0 nice -n 19)
SING=(singularity exec --cleanenv --env PYTHONPATH="$PYTHONPATH_POD" --env OMP_NUM_THREADS=2 --bind /tmp/if2453 --bind /sys/fs/cgroup:/sys/fs/cgroup:ro)
bio-box() { "${SCOPE[@]}" "${SING[@]}" "$IMAGE" python -m biobox.cli "$@"; }
biopy() { "${SCOPE[@]}" "${SING[@]}" "$IMAGE" python "$@"; }
```

`bio-box` and `biopy` run a short command in its own capped scope at nice 19 (AGENTS 7.2). Do not
use them for a long job: use `scripts/pod_run.sh` (sections 2 and 7).

The image (card C-29). `bio-box_64f1bf0.sif` is the runtime image built from `main` at `64f1bf0`,
pulled by digest from `ghcr.io/ifinkelstein/bio-box@sha256:2faf4b8de020430d3c5c9dace5569445fdb0eb0057088972793874a799e25bea`.
`DIGEST` is the sha256 of the `.sif` file, which is what `pod_run.sh --digest` compares (measured with
`sha256sum` on the pod, 2026-11-13). Its test image, `bio-box-test_64f1bf0.sif`, runs the test suite
(section 8). The image holds the locked Python environment (`uv.lock`: numpy, matplotlib, zombi2 and
the rest), MMseqs2 18.8cc5c, HMMER 3.4, Infernal 1.1.5 and skani 0.4.0. Only the tree's `src` goes on
the path: a dependency directory put before it would shadow the locked packages. Both `.sif` files
are in `$HOME/bio-box/images/` and in `/tmp/if2453/images/`.

Until card C-29 the pod ran everything in the test image `bio-box-test_74b3e05.sif` (built
2026-10-30; sha256 `c30be58557d756095346308c58cc3098af4bbe66404f67ab87be3f64f147c898`), which has
no numpy: the dependency directory `/tmp/if2453/pilot_s1/pydeps` supplied numpy 2.5.4, the
version `uv.lock` pins. The 50 Gbp build (unit `r50-driver`) runs in that image until it is done: a
release must not mix image digests in its run records, so do not change its driver, its tree or its
image while it runs. The next build takes the new image in its driver
(`scripts/stress/build_stage_pod/r50_driver.sh` or its copy): set `IMAGE` and `DIGEST` to the two
values above, and `PYPATH` to `$TREE/src` alone, with no `pydeps` directory after it.

## 1. Before you start

1. Look for a `STOP` file. If this prints a path, do not start anything.

```bash
ls -d "$ROOT/STOP"
```

`No such file or directory` is the good answer.

2. Look for other compute. The night cap holds only while no other user computes (AGENTS 7.2).
`run-step` and `build-stage` check this by themselves (card C-28): other users' processes using more
than half a core together over 5 seconds mean the day cap at night. The load average does not tell:
our own build holds it near 42. This prints what they read now, and the threads a 4-thread launch
would get (read-only; it starts nothing):

```bash
biopy -m biobox.runstep --profile "$PROFILE" --threads 4
uptime
top -b -n 1 -o %CPU -w 120 > "$WORK/top.txt"
head -n 15 "$WORK/top.txt"
systemctl --user list-units --type=service --type=scope --no-pager
```

The first command prints a number (the allowance) and one line: the cap and why, other users' CPU,
system accounts' CPU (root and uids up to 999: daemons such as the virus scanner, never counted),
ours elsewhere (subtracted from the cap) and ours in this unit, over the 5 s window, and the load
average beside them. The unit list shows the units of this account. Other agents and the PI use the
same account, so a unit you did not start is not yours. Processes of other users are in `top`.

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
  `pod_run: c06-build usage MemoryPeak=<bytes> CPUUsageNSec=<ns>` (the two in the order systemd prints
  them; before card W-20 the line was always empty, and a real unit printed `pod_run: w20-usage usage
  MemoryPeak=426106880 CPUUsageNSec=4151944000`, measured on wilkcomp02, 2026-11-11), is what the unit cost (card W-18):
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
- The cap is the thread count `run-step` allows one process: 40 by day; 56 at night and at weekends,
  only while no other user computes (card C-28: other users' processes, not root's or a system
  account's, use more than half a core together over a 5 s window; a reading that fails gives 40).
- J × C above the cap when `build-stage` starts is refused before any step (exit 2). The
  message names J, C and the cap (`run-step allows`).
- The allowance of a launch is the cap less what ours elsewhere use now: every process of this
  account outside the build's own unit and process tree (a test suite, another agent's build, a
  search), to the nearest CPU. The build's own shards are not subtracted; `build-stage` counts them.
- The allowance is read again before each step and before each shard launch. A build started at
  night keeps to 40 threads once the day window opens: new launches start only while the shards in
  flight fit the allowance, so fewer run at once, and a launch at night whose budget could still run
  when the day window opens (300 CPU-minutes over C threads: 75 minutes for C = 4, so from 06:45 on a
  weekday, derived) already counts against 40. If one shard's C is above the allowance, the launch waits (the
  log says `holding it`) and starts when the hour, other users or ours elsewhere allow. The same
  holds for a step whose `--threads` is above the allowance, the default `--threads` (the smaller
  cap) included: it waits while ours elsewhere leave less than it.
- `build-stage` writes the allowance and the readings behind it in one line at its start and again
  whenever it changes; every `run-step` writes one such line before its command (`S0.genecall:
  allowance ...`). `run-step` refuses only `--threads` above the cap; above the allowance it runs
  and says so, since a fixed-size job of ours is not squeezed out (`build-stage` waits instead).
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
agents' builds and tests, searches, and the PI's work, all together (AGENTS 7.2). Since card C-28
`build-stage` subtracts what ours elsewhere use at each launch, but it cannot see what will start
after its shards are launched, and a shard runs about an hour (3,291 to 3,727 s at 4 threads,
measured, `planning/measurements/2026-11-11_mmseqs_profile_vs_hmmsearch.md`). Count the threads the other jobs of
ours use now (the first command of section 1, `systemctl --user list-units`, `top`) and what you will
start, and choose J so that all together stay at or under the cap:

```text
J = floor((cap - other threads) / C)       C = --scan-cpus, default 4; J at least 1
```

`cap` is the one for the hour (the allowance). `other threads` is what else of ours runs, and what you
will start during the build. A worked example (derived, from the formula; the 40 is the profile's
`max_threads_day`): by day, with cap 40, one other build that runs 16 threads and C = 4, J =
floor((40 - 16) / 4) = 6, so `--scan-jobs 6 --scan-cpus 4`, 24 threads, 40 in all. With nothing else
running, the same formula gives J = 10, which is the whole cap: leave a margin if anything may start
during the hours the scan runs. A smaller J only makes the scan take longer; the output is the same.

**Non-coding RNA (card E-07, DEC-28).** `--ncrna-cm <Rfam.cm.gz>` puts the ncRNA producer in the
`features` slot, as a step named `ncrna-<producer>`. The producer scans the model set `--ncrna-set`
(default `prok-2`, DEC-30, the frozen file `src/biobox/l1/rfam_sets/prok-2.tsv`: 770 Rfam 15.1 models,
`prok-1`'s 588 plus 182 long-tail models with a hit on measurement 0-08a's sample; `--ncrna-set prok-1`
still scans card E-07's set, and a stage scanned with it keeps its producer name): its
intergenic models on the stage's intergenic regions widened by 200 bases into the flanking genes, with
`cmsearch --notrunc`, and its whole-genome models (RNA inside coding sequence; 34 in `prok-2`, 23 in
`prok-1`) on whole contigs. It
fetches the set's models from the library with `cmfetch` into scratch. Overlapping hits of models in one
Rfam clan on one strand compete (card E-08): the best bit score wins, every hit stays in `l1_ncrna` with its
`clan` and `clan_winner`, and only winners become features on the tape (clan file
`src/biobox/l1/rfam_sets/Rfam-15.1.clanin`, pinned in `refdb.yaml` as `rfam_clanin`). The producer's name
carries the library, the set, the region rule and the clan rule, for example
`infernal-1.1.5-rfam-15.1-prok2-igr200nt-clan1` (`-prok1-` with `prok-1`), so a new set or rule is a new
step and a new tape version
(card E-07's producer, without `-clan1`, kept every hit on the tape). The IGR, tape, token and posting steps then get `--producers
<gene caller>,<ncRNA producer>`: the release's tapes are those of both producers. The gene-call files are
the same with or without `--ncrna-cm`. The step counts as one step in `all N steps done`.

- Threads: each `cmsearch` runs with `--scan-cpus` threads (default 4, at most `--threads`).
- Shards: a scan estimated above 80% of `cpu_min_per_process` runs as `ncrna-<p>.plan`, one
  `ncrna-<p>.shard-<k>` per shard and the merge, as a domain scan does. The estimate is the stage's
  genome bases at the set's rate (card E-09): 304.0 CPU-hours per 10 Gbp for `prok-1`, the slower of the
  two shards measured on the pod (derived in `planning/measurements/2026-11-11_ncrna_producer.md`; card
  E-08), times the set's cost relative to `prok-1`. `prok-2` costs 299.7 against 254.5 (derived in
  `planning/measurements/2026-11-12_rfam_sample.md`), so its rate is 358.0 (derived: 304.0 × 299.7 /
  254.5). A shard holds at most 111,734,541 genome bases with `prok-2` and 131,578,947 with `prok-1`, at
  300 CPU-minutes (derived: 0.8 × 300 × 1e10 / (rate × 60)). A set file that is not a frozen set is
  planned at the dearest frozen set's rate.
- Memory: each `cmsearch` is counted at 6.4 GiB per two threads (measured peak 6.4 GiB at `--cpu 2` on
  intergenic records, measurement 0-09c). Shards at once are at most the unit's `MemoryMax` divided by
  that, so a 16G unit runs two shards of two threads. A unit too small for one `cmsearch` is refused
  (exit 2, `above the memory cap`) before the scan; raise `pod_run.sh --mem` or lower `--scan-cpus`.
- A contig another stage owns is skipped and named in the ledger (`owned_by:<stage>`), as gene calling
  skips it (DEC-26).
- Its plan and shard records are in `l1/ncrna_scan/producer=<producer>/stage=<S>/` (card E-08; a stage
  planned by card E-07's code has them in `l1/domain_scan/library=<producer>/stage=<S>/`, where they are
  still read), the shards' runs in `$BBTMP/scratch/biobox-ncrna-scan/<plan>/`; a run lost from scratch is
  scanned again before the merge (`python -m biobox.l1.ncrna rescan`, log: `scanning it again`), as for a
  domain scan (section 4).
- The ledger's tool versions and the step's result name the gene-call version whose intergenic regions were
  scanned (`gene_call_version`).

The library must be the `refdb.yaml` file for the name to say `rfam-15.1`; another file is named by its
sha256, which the release's C-08 check reports as differing from the lock.

**Clustering's CPU budget (DEC-29).** Clustering in one process has no exemption from the 300
CPU-minute rule any more: DEC-25 exempted it until clustering was split, and card B-42 ended that once
the split clustering of `rN.1b` (2,053,075 proteins) equalled its published clusters row for row
(`planning/measurements/2026-11-12_split_clustering_rN1b.md`). A release whose clustering would pass
the budget clusters in parts (below), every process under 300 CPU-minutes. The one exemption left is
DEC-29's: a single clustering call that no number of parts fits runs as one process under the
profile's `cluster_cpu_min_max` (20000 in `profiles/wilkcomp02.yaml`, set by the PI; below). Every
other step keeps the 300-minute rule.

- `build-stage` estimates clustering's CPU-minutes in one process before it starts it, from the
  release's protein count: a power law through 284.10 CPU-minutes at 2,053,075 proteins and 2,195.91
  at 10,620,684 (measured, wilkcomp02, 2026-11-10, `rN.1b` and `rN.1`; the fit is derived and
  extrapolated above 10.6 M proteins). Above 80% of `cpu_min_per_process` (240 of 300; about 1.79 M
  proteins, derived: 2,053,075 × (240 / 284.10)^(1 / 1.244)) the release clusters in parts.
- When every stage of the release has its gene-call run record, the check is made before any step;
  otherwise just before the cluster step's command, once this call's stage is gene-called.
- In one process (only below 80% of the budget, so it is never refused in practice) an estimate over
  `cpu_min_per_process` is refused, exit 2, whatever `cluster_cpu_min_max` says; the log line starts
  `build-stage: refused`, names DEC-29 and says to cluster in parts (`--cluster-parts N`, or let
  `build-stage` choose). The cluster step's run record holds no `cpu_budget_reason`.
- In parts, with or without the key, the check does not use the release's one-process estimate: it
  checks each call the published cuts run as one process (DEC-29, below), against `cluster_cpu_min_max`,
  and refuses such a call when the profile has no key (cards B-45, B-42). A release over the budget
  whose calls all fit in parts is not refused, with or without the key.
- A release clustered in one process by the code of cards B-36 to B-41 with the key has
  `cpu_budget_reason: exempt: DEC-25` in its `cluster` record; nothing reads it.

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
- The files are those of clustering in one process, byte for byte: tested on the fixture and a
  constructed set, and on the pod, where `rN.1b` (2,053,075 proteins) reclustered in parts by
  `build-stage` gave its published `cluster_member` and `cluster_member_aln` files with the same sha256,
  with the default cut and with three parts per call (card B-42,
  `planning/measurements/2026-11-12_split_clustering_rN1b.md`). The steps' records are in
  `l2/cluster_plan/release=<R>/`, never in a manifest.
- MMseqs2's databases are under `$BBTMP/scratch/biobox-cluster/<release>-<key>/`, and are removed when
  `cluster` publishes. At 2.05 M proteins they peaked at about 5.6 GB (measured, card B-42).
- The check above still comes first (in parts, call by call: above).

**A call that cannot be split runs as one process (DEC-29, card B-45).** The NR90 step's prefilter
(`cluster.nr90.11.prefilter`) and NR50 step 0 (`cluster.nr50.11.prefilter`) cannot search a shared index
in MMseqs2 18.8cc5c (card B-44), so each of their parts rebuilds the k-mer table of the whole target. From
about 68.3 M targets that table alone is 240 CPU-minutes a part (derived: 240 × 60 / 2.108e-4 CPU-seconds
per target), and no number of parts fits the budget: about 74 M proteins for the NR90 step and 108 M for
NR50 step 0 (derived, DEC-29). Then:

- With `cluster_cpu_min_max` in the profile, `build-stage` gives it to the plan
  (`release cluster plan --one-process-cpu-min-max`), which records it in its `cut`. The call's `.cut`
  step makes it one part, the workflow's own call, and its record holds `one_process` with the reason and
  the estimate (its work plus its own table). Above the cap the cut is refused, exit 2.
- Just before the call's step `build-stage` checks the estimate against the profile's cap again and
  refuses (exit 2, `refused before <release>.<call>`) above it or when the profile has no key. Its run
  record holds `cpu_budget_reason: exempt: DEC-29`, with `over_cpu_budget` true above 300 CPU-minutes.
- Without the key the cut refuses the call as before: `no number of parts fits it`, naming DEC-29.
- At 95 M proteins the NR90 step's prefilter as one process is about 3,945 CPU-minutes (derived: 3,629 of
  work and 316.4 for the table, the extrapolated model in `src/biobox/l2/cluster.py`), and the whole
  clustering in parts about 37,963 (derived, same place). Nothing that large has been measured.
- `release cluster plan` prints, under `one_process`, which calls may run as one process, from how many
  targets, and whether this release can reach that. `release cluster step --name <call>.cut` prints the
  cut's `one_process` when it is one.
- This ends when an MMseqs2 without the abort is pinned (card B-44).

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
| a scan shard's run missing from `$BBTMP/scratch/biobox-domain-scan/<plan>/` (or, for the ncRNA scan, `biobox-ncrna-scan/<plan>/`) while its run record exists | scratch was cleared (a reboot) or damaged between the shards and the merge | `build-stage`, before the merge, scans exactly those shards again into scratch (log: `scanning it again`); the new run must have the sha256 the shard's published record holds, else the build stops naming the shard (tell the PI) |
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
| `cpu_budget_reason` | only in the record of a clustering call run as one process in parts: `exempt: DEC-29` (section 2). Before card B-42 a release's `cluster` record clustered in one process with `cluster_cpu_min_max` in the profile held `exempt: DEC-25` |
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

Files that a later step adds to the release's stages for a library, producer or step the release
does not use there are not its concern (card C-19): a Pfam scan of a stage after the release was
built does not fail it. A second part beside a file the release lists still does. If the release
has a library's domain files for some of its stages only, `verify` prints a `warning:` line naming
the stages on each side; the exit code is unchanged. That release's family tokens cover only the
first stages' genes; a release built after the other stages are scanned covers them all.

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

**Before a release is announced: the worked examples** (AGENTS 5.7, card E-11). A release that verifies is not yet a
milestone. Before it is announced (in `STATUS.md`, `README.md` or `GUIDE.md` as the one to use, or copied to durable
storage as such), run every worked example under `examples/` on it with one command, read-only on the root, each
launch inside the account's thread cap (`examples/README.md` says what an example is and holds). It needs the
release's search pack in the profile's pack directory, and the examples' environment: `TREE`, `PROFILE`, `IMAGE` (with
matplotlib: the runtime image of section 0 has it), `DIGEST`, `ROOT`, and `OUT`, a new directory outside the root.
Look first, then start it as a unit of its own, since it runs for about an hour with the account near its cap:

```bash
export OUT="$WORK/examples-rN.3"
bash "$TREE/scripts/release_examples.sh" --dry-run rN.3
systemd-run --user --unit ex-rN3 --collect --quiet -p MemoryMax=2G -p MemorySwapMax=0 -p StandardOutput=append:"$WORK/examples-rN3.log" -p StandardError=append:"$WORK/examples-rN3.log" --setenv=TREE="$TREE" --setenv=PROFILE="$PROFILE" --setenv=IMAGE="$IMAGE" --setenv=DIGEST="$DIGEST" --setenv=ROOT="$ROOT" --setenv=OUT="$OUT" --setenv=UNIT=ex-rN3 --setenv=USER="$USER" /usr/bin/nice -n 19 bash "$TREE/scripts/release_examples.sh" rN.3
tail -n 1 "$WORK/examples-rN3.log"
```

The last line says whether every invariant holds; `$OUT/summary.tsv` has one row per check (example, check, kind,
expected, observed, verdict), and the script exits 1 when any invariant fails or is missing. **A failed or missing
invariant stops the announcement**: report it with the example's row, its README's table of what that failure would
mean, and the step's log under `$OUT/<example>/logs/`. When every invariant holds, finish `$OUT/note.md` (a
measurement-note skeleton: the run's facts, per-example counts and CPU, every check, `TODO` where the reading goes)
into `planning/measurements/`, and copy each `$OUT/<example>/recorded_new.tsv` to
`examples/<example>/recorded/<release>.tsv` in the same commit. A recorded value that `differs` on a new release is
expected; on a rerun of the same release it means the code or the release changed.

## 7. Search

The search path writes nothing under the data root (card W-03): a user who can read a release may
not be able to write its root. The search pack and the runs live outside the data root, under the
profile's `tmp`:

- the pack of release `R`: `$BBTMP/packs/release=<R>/<pack id>/`, one per host (a pack is a cache
  of a release, shared by every user). `$BBTMP/packs` (real: `/tmp/if2453/bio-box/packs`) is the
  host's one pack directory (card W-22): `pack`, `search`, `nbr` and `releases` find a release's pack
  there with only `--profile`. Only a `pack build` that makes a pack writes into it, so it may be
  read-only for everyone but its builder: `pack verify`, a `pack build` that finds the pack already
  built, and every search write nothing there (their scratch goes to `$BBTMP/search/$USER/work`), and
  a `pack build` that cannot write it exits 2;
- the runs: `$BBTMP/search/$USER/runs/run=<run id>/`, one directory per user; each running run's
  pending record is `$BBTMP/search/$USER/runs/pending/run=<run id>.json`, removed when the run closes
  (card W-20; the closed run's own record, `run.parquet`, then answers `status` and `cancel`);
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
the `pack build` command. `pack verify` re-reads every file of the pack and writes nothing into the
pack directory; it sorts the release's protein uids in `$BBTMP/search/$USER/work/.verify-*`, removed
when it ends (card W-22).

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
error `cancelled`. It prints nothing and exits 0. For a run that has already closed it exits 0 and
says on stderr that the run is closed; a run id that never existed is refused (exit 2). The unit's
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

`search hits` without `--format` prints the kept hits as a table (`--json` for JSON) and says on stderr
how many rows it shows of how many (`--top_n 0` shows all). `nbr summarize --by pfam,...` on a
release with no domain library says so on stderr and leaves the `pfam` level out. Strand is `+1` or
`-1` in every table, coordinates are 0-based and half-open, and the offsets of `neighbors.tsv` are
relative to the anchor's strand (card W-20; `search hmm` writes these and every column's meaning into
`README.txt`; `tier` is a column of `hits.tsv` only when the release has tiers, none has yet; its last two columns,
`stage` and `kind`, give each gene's build stage and its genome's kind, card W-23).
`--format tsv`
is one row per gene; `--format faa` is one record per protein, named by protein uid. `nbr ...
--format tsv --out DIR` writes one file per anchor locus. `nbr summarize` counts families and
clusters around the anchors. `export proteins` needs the release's manifest; it reads the root only
and takes no `--profile`. Every one of these given `--out` prints one stderr line per file it wrote,
`wrote PATH: N rows` (`sequences` for FASTA, `records` for JSON lines, `N tables, M rows` for a
directory of locus tables; card W-23).

**What a root holds, one sequence, and SQL** (card W-23). `releases` prints one line per release
of the root (a manifest or only a `release_def`): stages, genomes, contigs, distinct proteins, domain
libraries, ncRNA producers, whether it has the lookup tables of cards B-27 and B-39, and with
`--profile` (or `--pack-dir`) its pack. `seq` prints proteins by uid (`--protein`, several, or
`--proteins-file`) or one stretch of a contig (`--contig` by accession or contig uid, `--strand -1`)
as FASTA, read by key, never by exporting the release. `sql` runs SELECT statements on the release's
catalog through `open_catalog` on a read-only connection (external access off, only the release's and
the runs' files readable; `COPY ... TO`, `CREATE`, `ATTACH`, `SET` and the like exit 2), with DuckDB's
temporary files in a scratch directory under the profile's `tmp`, and prints on stderr the row groups
each statement read; `--runs` (or the profile) attaches runs of the release.

```bash
bio-box releases --root "$ROOT" --profile "$PROFILE"
bio-box seq --root "$ROOT" --release rC6.1 --protein 0123456789abcdef0123456789abcdef
bio-box seq --root "$ROOT" --release rC6.1 --contig NC_001416.1:201..1000 --strand -1
bio-box sql --root "$ROOT" --release rC6.1 --profile "$PROFILE" -e "SELECT key, value FROM catalog_info"
```

**Anchors from a file, and finding a gene** (card W-18). `nbr` takes a long list of protein uids from a
file, one uid per line (blank lines and `#` lines are skipped), in place of the comma-separated
argument; `find` prints the genes an accession or a protein uid names, one row each: accession,
contig uid, start, end (0-based, end exclusive), strand, feature uid, protein uid, and the gene in both
coordinate forms (`span`, `genbank_span`). `ACC:START-END` (0-based, end exclusive) or `ACC:START..END`
(GenBank: 1-based, both ends included; card W-23) names one gene, as `search hmm --positive` does (it
uses the same lookup); any other spelling exits 2 and shows both forms. `--accession ACC` alone lists
every gene of that contig; `--protein UID` every gene that carries the protein.

```bash
bio-box nbr summarize --anchors-file "$WORK/anchors.txt" --root "$ROOT" --profile "$PROFILE" --release rC6.1 --window 5g --by pfam,nr50
bio-box find --root "$ROOT" --release rC6.1 --accession NC_001416.1:200-1000
bio-box find --root "$ROOT" --release rC6.1 --accession NC_001416.1:201..1000
bio-box find --root "$ROOT" --release rC6.1 --protein 0123456789abcdef0123456789abcdef
```

**One command.** `search hmm` (card W-07) does all of the above for an HMM file: it writes the recipe,
runs it, and writes the exports into `--out DIR`. `--bits` is the inclusion threshold;
`--neighbors N` the genes on each side (window `Ng`, default 5); `--model` picks one model of a file
that holds several. `--positive` is required: a protein the search must find, given as
`ACCESSION:START-END` (one gene, with the accession, start and end that `search hits --format tsv`
prints: 0-based, end exclusive), as `ACCESSION:START..END` (the same gene in GenBank's 1-based form) or
as a substring of its contig's header (case-sensitive). It must
name exactly one protein of the release; otherwise the command exits 2 and lists the genes that
matched, each as `ACCESSION:START-END`, so one can be given instead. `$POSITIVE` below is such a
gene. The run runs in the foreground, as `search run --wait`, so it may run in a `pod_run.sh` unit;
`search status` and `search cancel` work on it from another shell.

```bash
bash "$TREE/scripts/pod_run.sh" --profile "$PROFILE" --image "$IMAGE" --digest "$DIGEST" --bind /tmp/if2453 --env PYTHONPATH="$PYTHONPATH_POD" w07-hmm search hmm --root "$ROOT" --release rC6.1 --profile "$PROFILE" --hmm "$WORK/inputs/lib.hmm" --model lam0 --bits 25 --neighbors 5 --positive "$POSITIVE" --out "$WORK/hmm"
```

`DIR` must be new or empty and outside the data root. It then holds `README.txt` (what each file and
column means), `recipe.yaml` (the recipe the
run ran; `search run DIR/recipe.yaml` runs it again), `hits.tsv` and `hits.faa` (as `search hits
--format tsv|faa`), `neighbors.tsv` (as `nbr summarize --format tsv`) and `loci/` (one table per
locus, as `nbr --format tsv`); each is written by that command, so the files are the same. The log
has a `wrote PATH: N ...` line for each of `hits.tsv`, `hits.faa`, `neighbors.tsv` and `loci/`, then
ends with the run id and `ok`, then one line per file. The run itself is in the runs directory, as
for `search run`. A release with no domain library (the 10 Gbp release `rN.1` has none) is
searched all the same: the neighbour summary is by `nr50,nr90` instead of `pfam,nr50`, and the log
says so in a line that begins `note`. A run that ends other than `ok` exits 1 and leaves only
`recipe.yaml` in `DIR`. When the release has no pack, the command exits 2 and prints the
`pack build` command to run; it does not build the pack.

**Keep the results** (card W-28, DEC-33). The runs directory and every `--out` are on `/tmp`, which
a reboot may clear. `results` keeps a closed run on the profile's durable storage, in the results
area `$DURABLE/results` (real: `/stor/work/Wilke/if2453/all-the-data/results`), with no file under
64 MB on `/stor` (AGENTS 7.3) and nothing there overwritten or deleted (AGENTS 7.5):

- `results archive RUN --with DIR ...` writes the run's directory (its record and tables) and each
  output named with `--with` (the `search hmm` DIR, the `nbr --format tsv` DIR, an `nbr flanks` DIR,
  a lone `hits.tsv`) as one tar, with one manifest line, into the spool `$BBTMP/results/spool/` on
  local disk. Since card W-29 an output directory written from a run holds a `RUN` file (the run id and
  the release), and `--under DIR` finds every such directory up to four levels below DIR, so naming each
  one is needed only for single files and for directories written before that card; a DIR whose `recipe.yaml` is not the run's recipe is another
  run's and is refused. The same run archived again is a no-op (exit 0, said on stderr); with other
  files or outputs it is refused (exit 2), naming both sha256. The tar has sorted names, mtime 0, no
  owner, no host or path, and the run's own date (`submitted_at`) in its manifest line: the same run
  gives the same bytes on any host.
- `results flush` packs every spooled run into one bundle, `results-<first day>-<last day>-<hash>.tar`,
  and renames it into the results area only when it is at least `small_file_floor_bytes` (64 MiB);
  below that it writes nothing there, the runs stay spooled, and it says so on stderr (exit 0). Most
  runs are far smaller (the `rN.3` RNA polymerase run with its `search hmm` DIR is 12,492,800 bytes,
  measured below), so they wait in the spool until enough gather. `--below-floor` writes a smaller
  bundle; it breaks AGENTS 7.3 and needs the PI's agreement. A bundle is copied to a dot-named
  `.results-*.part-*` file there, read back, then renamed; a flush killed before the rename leaves
  that file, which no command reads and `results list` names: only the PI removes it.
- `results list` prints one line per kept run, flushed or spooled, with where it is, from the
  manifests only (a bundle's first member). `results get RUN --out DIR` restores the run into a new
  or empty DIR (`FILES.tsv`, `run/`, `outputs/<name>/`) and checks every file's sha256; a mismatch
  exits 1 and writes no DIR.

The image sees the durable storage only when it is bound, so these run through `pod_run.sh` with
`--bind "$DURABLE"`. Until the PI has made `/stor/work/Wilke/if2453/all-the-data` (it is not on the
pod as of 2026-11-12), every `results` command with the real profile fails: without the bind it exits
2 with `the profile's durable storage ... does not exist`, and with it the image cannot start, since
it cannot bind a missing directory. The first flush onto `/stor` is the PI's step.

```bash
bash "$TREE/scripts/pod_run.sh" --profile "$PROFILE" --image "$IMAGE" --digest "$DIGEST" --mem 8G --threads 1 --bind /tmp/if2453 --bind "$DURABLE" --env PYTHONPATH="$PYTHONPATH_POD" w28-archive results archive "$RUN" --profile "$PROFILE" --with "$WORK/hmm" --with "$WORK/nbr"
bash "$TREE/scripts/pod_run.sh" --profile "$PROFILE" --image "$IMAGE" --digest "$DIGEST" --mem 8G --threads 1 --bind /tmp/if2453 --bind "$DURABLE" --env PYTHONPATH="$PYTHONPATH_POD" w28-flush results flush --profile "$PROFILE"
bash "$TREE/scripts/pod_run.sh" --profile "$PROFILE" --image "$IMAGE" --digest "$DIGEST" --mem 8G --threads 1 --bind /tmp/if2453 --bind "$DURABLE" --env PYTHONPATH="$PYTHONPATH_POD" w28-list results list --profile "$PROFILE"
bash "$TREE/scripts/pod_run.sh" --profile "$PROFILE" --image "$IMAGE" --digest "$DIGEST" --mem 8G --threads 1 --bind /tmp/if2453 --bind "$DURABLE" --env PYTHONPATH="$PYTHONPATH_POD" w28-get results get "$RUN" --profile "$PROFILE" --out "$WORK/back"
tail -n 4 "$BBTMP/logs/w28-flush.log"
```

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

The image is the test image built from `main` at `64f1bf0`,
`$HOME/bio-box/images/bio-box-test_64f1bf0.sif` (sha256 of the file
`3d77df6b7f0b894b5371e42eaf077dfda09ebc97b1f78f2ed1d2cfc5fa9fbda4`, pulled from
`ghcr.io/ifinkelstein/bio-box-test@sha256:cb7cd0f1e7146350576e4658de39290647d1cbba39494ce81ba072f3041707bb`),
with the tree's `src` alone on `PYTHONPATH` (section 0). `BIOBOX_POD_IMAGE` names another image and
`BIOBOX_POD_PYDEPS` adds directories after `src`, for an image that lacks a package; give neither
for the image above. A suite that must survive a dropped connection runs as detached units on the
pod: sync the tree with the first command above, then, on the pod,
`bash "$TREE/scripts/orchestrator/detached_suite.sh" <name> 2` (the script's header says how to
poll it).

The image check. Before a new image is used, run its own suite on the pod: the source and tests
built into the image (`/opt/bio-box`), no `PYTHONPATH`, as the `image` workflow runs it on Actions
(`.github/workflows/image.yml`). Run on the pod, one suite at a time, after checking the thread cap
(section 1):

```bash
BIOBOX_POD_IMAGE="$HOME/bio-box/images/bio-box-test_64f1bf0.sif" BIOBOX_PEER_MANIFEST=/tmp/if2453/ci-peer/right.json bash "$TREE/scripts/orchestrator/detached_suite.sh" --own image-64f1bf0 2
systemctl --user is-active ci-image-64f1bf0-0 ci-image-64f1bf0-1
tail -n 1 /tmp/if2453/ci-tmp/image-64f1bf0-detached/log.0 /tmp/if2453/ci-tmp/image-64f1bf0-detached/log.1
```

Add the two summary lines and compare with the Actions run that built the image. For `64f1bf0`:
Actions ran 1,720 passed, 22 skipped (measured, run 43380434543); the pod ran 1,722 passed, 20
skipped, no failure (measured, wilkcomp02, 2026-11-13, card C-29: shards of 894 passed, 14 skipped
in 53 min and 828 passed, 6 skipped in 80 min, beside the 50 Gbp build). The two that pass only
here are `test_d05_cross_host_manifest_hash_equal` (`tests/r3/test_cross_host.py`: the peer file
was given, and the fixture's manifest hash was `f6126b00d45e9284d2799c3a699c992c` on both hosts)
and `test_c07_detects_the_real_cap_it_runs_under` (`tests/r1/test_runstep_units.py`: it needs a
cgroup memory cap, which a unit on the pod has and `docker run` on Actions has not). That pair is
derived: the Actions log prints counts, not skip reasons.
`BIOBOX_PEER_MANIFEST` is a file holding the `cross-host` workflow's manifest hash for the commit
(card D-04, `tests/r3/test_cross_host.py` says how to make it); without it that test skips.

The tree's suite in the new image. `plan/task-cards` with card C-29, in `bio-box-test_64f1bf0.sif`
with the tree's `src` alone on the path: 1,740 passed, 19 skipped, no failure; the same tree in
`bio-box-test_74b3e05.sif` with `/tmp/if2453/pilot_s1/pydeps` after `src`: 1,696 passed, 63
skipped (both measured, wilkcomp02, 2026-11-13, card C-29, `detached_suite.sh` with 2 shards). The
44 that skipped in the old image and pass in the new are the tests that need zombi2 (32), matplotlib
(8), samtools (3) or polars (1), none of which the old image has; no test changed outcome otherwise.
The old image's only Python package from the dependency directory was numpy 2.5.4, the version the
new image holds; every package in both has the same version in both.

## 9. Clean up /tmp safely

What is ours, all under `/tmp/if2453/`:

- `$BBTMP/staging/<stage>.<step>/`: leftovers of killed or failed steps.
- `$BBTMP/scratch/`: builder scratch: one directory per step, `<stage>.<step>/`, which
  `build-stage` removes itself (section 4), and the scan shards' runs in `biobox-domain-scan/` and
  `biobox-ncrna-scan/`, removed after the merge.
- `$BBTMP/logs/`: unit logs.
- `$BBTMP/packs/release=<R>/<pack id>/`: the search pack of a release (section 7). It is a cache of
  the release, shared by every user of the host: `pack build` makes it again from the release. Its
  `pack_verified.json` is the build's record that the pack verified (card W-19); a search reads it and
  writes nothing into the pack directory. A pack built before W-19 has none, and every search of it
  verifies it in full first (its log says so); `pack build` with the same inputs adds the record
  without building again.
- `$BBTMP/packs/release=<R>/.building-<pid>-<hex>/` and `.scratch-<pid>-<hex>/`: a pack build's
  temporary pack and its sort's scratch (the build's check of its new pack sorts there too). Each is
  removed by the command that made it, also when it is stopped by SIGTERM, SIGHUP or Ctrl-C; one is
  left only when that command was killed (SIGKILL, the OOM killer, a reboot). They are dot-named, so
  `ls` shows them only with `-a`. Remove one only when no `pack build` of that release is running
  (`ps -u "$USER" -o pid,args`), by exact path. Nothing else writes into `$BBTMP/packs` (card W-22).
- `$BBTMP/search/$USER/`: your searches (section 7). `runs/run=<run id>/` is a closed run's record
  and hits; `runs/pending/run=<run id>.json` is a running run's pending record, which the run
  removes when it closes (`.stderr` beside it is the runner's stderr of a run submitted outside a unit,
  and stays); `work/<run id>/` is a run's work directory, which the
  run removes when it closes. One is left only when the runner was killed (SIGKILL, a reboot).
- `$BBTMP/search/$USER/work/.verify-<pid>-<hex>/`: the scratch of a `pack verify`, or of the verify a
  `pack build` runs on a pack it finds already built (card W-22; before W-22 it was in the pack
  directory, `$BBTMP/packs/release=<R>/`, where one left by an older tree may still be). Removed by
  the command that made it; one is left only when that command was killed. Remove one only when no
  `pack verify` or `pack build` is running, by exact path.
- `$BBTMP/results/spool/` (card W-28, section 7): runs archived and not yet flushed, `<run id>.tar` and
  `<run id>.tsv` each, which `results flush` removes once their bundle is in the results area and read
  back. `.archive.tmp-*` and `.bundle.tmp-*` there are what a killed `results archive` or `results flush`
  left; remove one by exact path when no `results` command runs. `.lock` stays.
- `$TREE`, `$WORK` and `/tmp/if2453/ci-tmp/<name>-*`: your mirrored trees and pytest directories.

Never delete:

- anything under `/stor`, the results area's `.results-*.part-*` files included (ask the PI);
- a spooled run, `$BBTMP/results/spool/<run id>.tar` or `.tsv`, that `results list` shows only in the
  spool: it is the run's only archive;
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
| "Other compute" | other users' processes (uid above 999 and not ours; root's and system accounts' are daemons, never counted) above 0.5 of a core together, over a 5 s window (two reads of every `/proc/<pid>/stat`); a failed reading counts as busy. Estimates: other login users together used 0.01 cores and root's daemons about 0.3 (ps, wilkcomp02, 2026-11-13, measured); one single-threaded job uses 1. The load average is reported, not used: our own scan holds it near 42 | `OTHERS_BUSY_CPUS`, `CPU_WINDOW_S`, `SYSTEM_UID_MAX`, `measure_cpu` in `src/biobox/runstep.py`; card C-28 |
| Allowance of a launch | the cap (above) less ours elsewhere: this account's processes outside the build's unit (its cgroup) and process tree, to the nearest CPU, at least 0. At night a launch of T threads counts against the day cap when 300 / T minutes from now reach the day window | `allowance`, `launch_cap` in `src/biobox/runstep.py`; card C-28 |
| Memory per process or unit | 200G, swap 0 | profile `mem_scope`; AGENTS 7.2; `MemorySwapMax=0` in `scripts/pod_run.sh` |
| Memory cap required | true: `run-step` and the search runner refuse to start uncapped | profile `require_mem_cap` |
| Nice | 19 (one big process), 15 (many small ones) | profile `nice_monolithic`, `nice_fanout`; AGENTS 7.2 says 19 |
| CPU-minutes per process | 300 | profile `cpu_min_per_process`; AGENTS 7.2 |
| CPU-minutes of one clustering call run as one process | 20000, estimated before launch; refused above it. The one exemption from the 300 (DEC-29): a call that no number of parts fits (the NR90 step's prefilter from about 74 M proteins, NR50 step 0 from about 108 M, derived). Set by the PI on 2026-11-10 for the 50 Gbp build: about nine times the 2,195.91 measured for 10,620,684 proteins (wilkcomp02, 2026-11-10). Without the key such a call is refused. Clustering in one process has no exemption since card B-42 (DEC-25's ended): above 80% of the budget a release clusters in parts | profile `cluster_cpu_min_max`; DEC-29 (narrowing DEC-25); `cluster_refusal`, `one_process_refusal` in `src/biobox/stage.py`; cards B-36, B-45, B-42 |
| ncRNA scan shard size | `prok-2` (the default): 111,734,541 genome bases, 80% of `cpu_min_per_process` (300) at 358.0 CPU-hours per 10 Gbp (derived: 304.0 × 299.7 / 254.5, the set's cost relative to `prok-1` from `planning/measurements/2026-11-12_rfam_sample.md`). `prok-1`: 131,578,947 at 304.0, the slower of the two shards of card E-07's pod check (derived in `planning/measurements/2026-11-11_ncrna_producer.md`; the design note's 247.4 gave 161,681,487). Shards at once: at most the unit's `MemoryMax` / (6.4 GiB per two `cmsearch` threads) | `NCRNA_CPU_H_PER_10GBP`, `NCRNA_SET_COST`, `NCRNA_JOB_GIB` in `src/biobox/stage.py`; cards E-07, E-08, E-09 |
| Scan shard size | 30,150,753 residues: 80% of `cpu_min_per_process` (300) at 7.96 CPU-minutes per million residues (derived: 300 × 0.8 / 7.96 × 10^6; the rate of the Pfam scan of stage P1 in card C-10, measured on wilkcomp02, 2026-11-10: 2,343.74 CPU-minutes for 294.5 M residues). The 80% leaves room for a scan 25% slower per residue (derived: 1 / 0.8) | `SHARD_MARGIN` and `SCAN_CPU_MIN_PER_M` in `src/biobox/stage.py`; cards C-11, C-13 |
| Scan shards at once | `--scan-jobs` J × `--scan-cpus` C at most the cap `run-step` allows when `build-stage` starts (40 by day, 56 at night when no other user computes); defaults C = 4, J = cap // C (10 by day, 14 at night). The cap is a total over everything of ours on the pod: `build-stage` subtracts what ours elsewhere use at each launch, not what starts later, so choose J = floor((cap - other threads) / C) (section 2). The allowance is read again before each step and before each shard launch; a launch that would pass it waits | `scan_slots` and `_launch` in `src/biobox/stage.py`; cards C-11, C-12, C-28 |
| Unit CPU quota | the night cap, 56: a ceiling fixed at launch; `build-stage` holds the limit of the hour | `scripts/pod_run.sh --threads`, profile `max_threads_night`; card C-12 |
| Smallest data file | 64 MiB (67108864 bytes), except tables the registry marks `small_ok` | profile `small_file_floor_bytes`; AGENTS 7.3 says 64 MB |
| Work directory | `/tmp/if2453/bio-box` | profile `tmp`; AGENTS 7.3 |
| Search pack directory | `<tmp>/packs` (`/tmp/if2453/bio-box/packs`), the host's one, outside the data root; with `--profile`, else `--pack-dir`. Only a `pack build` that makes a pack writes into it: it may be read-only for all but its builder | profile `tmp`; `default_pack_base` in `src/biobox/search/runner.py`; cards W-03, W-22 |
| Search runs and work directories | `<tmp>/search/<user>/runs` and `<tmp>/search/<user>/work`, one per user, outside the data root; with `--profile`, else `--runs` and `--work` | profile `tmp`; `default_runs_dir`, `default_work_base` in `src/biobox/search/runner.py`; card W-03 |
| `search run` in a `pod_run.sh` unit | only with `--wait`; without it, refused (exit 2) | `BIOBOX_POD_UNIT` in `scripts/pod_run.sh` and `src/biobox/cli.py`; card W-13 |
| Durable results | `/stor/work/Wilke/if2453/all-the-data` | profile `durable` |
| Results area | `<durable>/results` (`/stor/work/Wilke/if2453/all-the-data/results`), or `--results` inside `durable` (decided on the resolved path and the file system's device, never the spelling); only bundles of at least `small_file_floor_bytes`, never overwritten or deleted; `--below-floor` only with the PI's agreement | profile `durable`, `small_file_floor_bytes`; DEC-33; `biobox.search.results`; card W-28 |
| Results spool | `<tmp>/results/spool` (`/tmp/if2453/bio-box/results/spool`), local disk: archived runs wait there until a flush reaches the floor | profile `tmp`; card W-28 |
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
`/tmp/if2453/c06`, `small_file_floor_bytes` 1000, and no `cluster_cpu_min_max` (so no DEC-29
exemption; the fixture's clustering is far under 300 CPU-minutes). Appendix B lists it.

## 11. Known limits today

| Limit | What happens | Card |
|---|---|---|
| Clustering's budget | above 80% of the budget a release's clustering runs in parts (card B-41, section 2); in one process there is no exemption (card B-42). In parts only a call run as one process is checked: refused over the profile's `cluster_cpu_min_max` (20000) or without the key, exit 2 (cards B-45, B-42). The one-process estimate that chooses parts is a fit by protein count, extrapolated above 10.6 M proteins. From about 74 M proteins a part of the NR90 step's prefilter would spend 240 CPU-minutes or more on its own k-mer table of the whole target (derived, section 2), so no part fits the budget; that call then runs as one process under `cluster_cpu_min_max` (DEC-29), or is refused without the key | B-36, B-37, B-41, B-45, B-42 |
| Ingest of a stage over about 60 Gbp | `build-stage` refuses before any step (exit 2): one ingest process runs at about 3.7 Mbp/s (measured, wilkcomp02), so it would pass 300 CPU-minutes. Split the stage into smaller stages | B-30 |
| Shard size is Pfam's | shards are cut at 80% of the budget at Pfam's rate (7.96 CPU-minutes per million residues, measured in C-10). A library more than 1.25 times slower per residue than that (derived: 1 / 0.8; above 9.95 CPU-minutes per million residues, derived) passes 300 CPU-minutes per shard, and the run record says `over_cpu_budget`. Lower `--scan-cpus` does not help, since the budget counts CPU-minutes | C-11, C-13 |
| `pod_run.sh` CPU quota | the unit's quota is `--threads` of `pod_run.sh` (default the night cap, 56), fixed at launch, whatever `--scan-jobs` × `--scan-cpus` is. By day it is above the day cap: the day cap is held by `build-stage` and `run-step`, not by the quota. A `run-step` launched alone keeps the thread count it chose at its start until it ends (at most 300 CPU-minutes) | C-07, C-11, C-12 |
| What the allowance sees | each reading is 5 s of CPU by owner, taken once per launch: work of ours that starts after a shard is launched is subtracted only at the next launch (a shard runs about an hour), a process born and ended inside the window is not seen, and `run-step` does not refuse `--threads` above the allowance (only above the cap). So the total can pass the cap until the build's next launch; leave a margin (section 2). The reading costs 5 s of wall time per launch and per `run-step` | C-28 |
| The tree is newer than the image | the images on the pod are built from `main` at `64f1bf0` (2026-11-13), which holds the locked environment; `plan/task-cards` has cards since, so the tree's source is bound in with `PYTHONPATH`, with no dependency directory (section 0). A new image is checked with its own suite before use (section 8). The 50 Gbp build (`r50-driver`) keeps the old image `bio-box-test_74b3e05.sif` until it is done (section 0) | 0-02, C-29 |
| `rL.*` releases on `/stor` | they are in the pre-B-17 layout: current readers and `verify_manifest` may not read them. `rL.1` and `rL.2` also fail verification on one file, `l0/l0_manifest.parquet` | B-17, B-22 |
| Search packs and runs inside an old root | a root written before card W-03 may hold `search_pack/` and `search/` (the pack and runs were then written into the root). Nothing reads or writes them now; the manifest leaves them out (`catalog.EXCLUDED_TREES`). Build the pack again into the pack directory (section 7) | W-03, W-13 |
| No detached search in a unit | `search run` without `--wait` in a `pod_run.sh` unit is refused (exit 2): the unit would stop the detached runner when the command returns (section 7). Use `--wait`; `search status` and `search cancel` work on a `--wait` run from another shell | W-13 |
| `pod_run.sh` and the profile | `pod_run.sh --profile` sets the unit's limits only; it does not pass the profile to `bio-box` (no `BIOBOX_PROFILE` in the image). Give `bio-box` its own `--profile`, or the pack, runs and work directories (section 7) | W-03, W-13 |
| Older outputs are named by hand | An output directory written before card W-29 has no `RUN` file, and neither has a lone file such as `hits.tsv`: give those to `results archive` with `--with` (a DIR whose `recipe.yaml` is not the run's is refused). Directories written since are found with `--under`, and the run's record holds the user's command line (`limits_json`, key `command`) | W-28, W-29 |
| The durable storage is not there yet | `/stor/work/Wilke/if2453/all-the-data` (profile `durable`) does not exist on the pod (2026-11-12); `results` with the real profile exits 2 until the PI makes it. The image sees it only with `--bind` | W-28 |
| Spooled runs are on `/tmp` | until a flush reaches the 64 MiB floor, a run's archive is only in the spool on local disk, which a reboot may clear, as it clears the run itself | W-28 |
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

Rows dated before 2026-11-13 ran in the image that section 0 named then, `bio-box-test_74b3e05.sif`,
with the dependency directory `/tmp/if2453/pilot_s1/pydeps` after the tree's `src`. The rows marked
C-29 were run on 2026-11-13 by card C-29 with the images of sections 0 and 8, the tree of the C-29
branch (`TREE=/tmp/if2453/ci/C-29`) and the real profile, in units named `c29-*` and `ci-*`.

The rows of section 7, and the rows marked W-13, were run on 2026-11-10 between 15:35 and 15:42 (pod
local time) by card W-13, in `/tmp/if2453/w13` (`WORK`; `ROOT=$WORK/root`, `BBTMP=$WORK/tmp`,
`PROFILE=$WORK/profile.yaml`), with units named `w13-*`, the tree of the W-13 branch
(`TREE=/tmp/if2453/ci/W-13`) and the same image. Its scratch profile is Appendix B's with `tmp`,
`durable` and `bulk` under `/tmp/if2453/w13` and `max_threads_day` and `max_threads_night` 4. The
fixture of Appendix B was built there first, with the section 2 command as unit `w13-build` and
`--threads 4` (`pod_run: w13-build exited 0`). The two searches that were cancelled and stopped
used the runner's test hook `--fault timeout --wall-s 900` (its worker waits until it is stopped),
so that they ran long enough to be watched; `--fault` is not for real searches.

The rows marked W-28 were run on 2026-11-12 between 23:44 and 23:53 (pod local time) by card W-28, in
`/tmp/if2453/w28`, on the `rN.3` RNA polymerase run `599bea5d471c84759a0eec02e2ca9225` of card W-25 and
its outputs in `/tmp/if2453/w25` (bound read-only), with two scratch profiles: the real profile with
`tmp`, `durable` and `bulk` under `/tmp/if2453/w28/one` or `/tmp/if2453/w28/two`, and the real floor
(67108864). `DURABLE` was the scratch `durable`, never `/stor`. Each ran in a unit `w28-*` with
`MemoryMax=8G`, `CPUQuota=100%` and nice 19, through the script `/tmp/if2453/w28/run.sh` (the same
`singularity exec` as `pod_run.sh`, waiting for the unit), and `results list` once through
`pod_run.sh` as written in section 7.

| Section | Command | Date run | Last line |
|---|---|---|---|
| 1 | `ls -d "$ROOT/STOP"` | 2026-11-10 | `ls: cannot access '/tmp/if2453/c06/root/STOP': No such file or directory` |
| 1 | `uptime` | 2026-11-10 | ` 07:26:38 up 11 days, 12:56,  7 users,  load average: 23.99, 23.10, 17.69` |
| 1 | `top -b -n 1 -o %CPU -w 120 > "$WORK/top.txt"` | 2026-11-10 | `(no output)` |
| 1 | `head -n 15 "$WORK/top.txt"` | 2026-11-10 | `2375504 if2453      39  19 8039140   1.2g  97688 R  18.8   0.1   1:50.74 python` |
| 1 | `systemctl --user list-units --type=service --type=scope --no-pager` | 2026-11-10 | `To show all installed unit files use 'systemctl list-unit-files'.` |
| 1 | `df -h /tmp` | 2026-11-10 | `/dev/sdc2       1.1T  165G  879G  16% /` |
| 1 | `biopy -m biobox.doctor --profile "$PROFILE"` | 2026-11-10 | `}` (the JSON shows `"cgroup_memory_max": "8 GiB"`) |
| 1 | C-29: `biopy -m biobox.doctor --profile "$PROFILE"` in `bio-box_64f1bf0.sif` | 2026-11-13 | `}` (the JSON shows `"cgroup_memory_max": "8 GiB"`, `mmseqs` `18.8cc5c`, `HMMER 3.4`, `INFERNAL 1.1.5`, `pyrodigal` `3.7.1`) |
| 2 | `pod_run.sh ... c06-build build-stage ...` | 2026-11-10 | `pod_run: log /tmp/if2453/c06/tmp/logs/c06-build.log; watch: systemctl --user status c06-build` (the first launch ran without `--digest`; the same command with `--digest` ran as `c06-rebuild`, below) |
| 2 | C-29: `pod_run.sh --profile "$PROFILE" --image "$IMAGE" --digest "$DIGEST" --mem 2G --threads 1 ... c29-digest --help` (section 0's new image and digest) | 2026-11-13 | `pod_run: c29-digest exited 0` (before it, `pod_run: image ... sha256:cb7353f8...afba3` at launch) |
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
| 7 | W-28: `results list --profile` the real profile | 2026-11-12 | `refused: the profile's durable storage /stor/work/Wilke/if2453/all-the-data does not exist (it is made by the PI, not by bio-box)`, exit 2 |
| 7 | W-28: `results archive "$RUN" ... --with` the `search hmm` DIR alone (profile `two`) | 2026-11-12 | `spool /tmp/if2453/w28/two/tmp/results/spool: 1 runs, 12492800 bytes; ...` (2,590 files) |
| 7 | W-28: `results flush` of that spool | 2026-11-12 | `results flush: the spool's 1 runs make a bundle of 12503040 bytes, under the floor of 67108864 bytes (...); nothing written to /tmp/if2453/w28/two/durable/results, the runs stay spooled. ...`, exit 0; the results area stayed empty |
| 7 | W-28: `results archive "$RUN" ...` with the `search hmm` DIR and the five `nbr flanks` DIRs (profile `one`) | 2026-11-12 | `spool ...: 1 runs, 164618240 bytes; ...` (14,264 files, sha256 `64cdbd66...7b77`, 9.6 s) |
| 7 | W-28: the same command again | 2026-11-12 | the same line; before it `results archive: run 599bea5d471c84759a0eec02e2ca9225 is already kept, sha256 64cdbd66...7b77, in .../spool/599bea5d471c84759a0eec02e2ca9225.tar; nothing written`, exit 0 |
| 7 | W-28: `results archive "$RUN" ...` with other outputs into profile `two`, which keeps the run already | 2026-11-12 | `refused: run 599bea5d471c84759a0eec02e2ca9225 is already archived in ... with sha256 9c305398...d7cf; this archive would be sha256 481bbb76...6fef (...). An archive is never overwritten`, exit 2 |
| 7 | W-28: `results flush --profile "$PROFILE"` (profile `one`, unit `w28-flush`) | 2026-11-12 | `/tmp/if2453/w28/one/durable/results/results-20261112-20261112-493d4af66eeb856b.tar	1 runs	164628480 bytes	sha256 5ae58462ba605c4289e08804eb03a605ad34f42a0a86db54d6606785026b379c` (7.8 s; `sha256sum` of the file printed the same) |
| 7 | W-28: `pod_run.sh ... w28-podrun-list results list --profile "$PROFILE"` (with `--bind "$DURABLE"`) | 2026-11-12 | `pod_run: w28-podrun-list exited 0` (before it, the run's line naming the bundle, and `spool ...: 0 runs, 0 bytes; ...`) |
| 7 | W-28: `results get "$RUN" --profile "$PROFILE" --out "$WORK/back"` (unit `w28-get`) | 2026-11-12 | the run's manifest line; before it `wrote /tmp/if2453/w28/back: 14264 files` (8.7 s). `diff -r` of `back/run` and each `back/outputs/<name>` against the originals: no difference |
| 8 | `sh scripts/pod_pytest.sh C-06 tests/r0 -k nothing_matches_this` (laptop) | 2026-11-10 | `140 deselected in 1.09s` (exit 5: nothing selected) |
| 8 | `sh scripts/pod_pytest.sh C-06 tests/r0/test_runbook.py` (laptop) | 2026-11-10 | `20 passed in 3.48s` |
| 8 | `sh scripts/pod_pytest.sh C-06 --shards 8` (laptop) | 2026-11-10 | `shard 7: 121 passed, 8 xfailed in 509.90s (0:08:29)` (eight shard lines, no failure; 1023 passed, 10 skipped, 8 xfailed in all) |
| 8 | C-29: the image check, `detached_suite.sh --own image-64f1bf0 2` (pod) | 2026-11-13 | the full suite ran through a script with the same units, binds and test split before `--own` was added to `detached_suite.sh`: shard lines `894 passed, 14 skipped, 59 warnings in 3195.94s (0:53:15)` and `828 passed, 6 skipped, 24 warnings in 4777.45s (1:19:37)`. `--own` itself, on a two-test selection (`-k 'test_the_runbook_has_commands or test_d05_cross_host_manifest_hash_equal'`): `1 passed, 907 deselected in 6.98s` and `1 passed, 1 skipped, 832 deselected, 1 warning in 12.34s` |
| 8 | C-29: `sh scripts/pod_pytest.sh C-29 tests/r0/test_runbook.py -k c29` (laptop) | 2026-11-13 | `5 failed, 2 passed, 38 deselected in 3.97s` (before the change, in the old image: the card's tests failing first) |
| 8 | C-29: `sh scripts/pod_pytest.sh C-29 tests/r0` (laptop; the new default image) | 2026-11-13 | `203 passed in 7.72s` |
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

- `results flush` onto `/stor` (card W-28): the PI's step. It was run only into a scratch results
  area under `/tmp/if2453/w28`, at the real floor.

- `$WORK/rb.py` (Appendix B), which replaced the `python -c` one-liners of sections 5 and 6 and Appendix B in card C-14. It holds the same code, but was not run on the pod: C-14 may run only the test suite there. Run each of its commands once, and through `ssh`, and correct this line.
- Recovery from an intent file (a kill during publishing). It cannot be hit by hand at fixture
  size; `tests/r2/test_runstep_resume.py` kills a step at each publish point with an injected hook.
- A reboot. A kill by `systemctl --user stop` was tested; the machine was not rebooted.
- A hard kill of a unit (SIGKILL). `systemctl --user kill --signal=SIGKILL` printed
  `Failed to send signal SIGKILL to auxiliary processes: Invalid argument`, though the unit ended;
  it is not documented here.
- The night thread cap in a build. Until card C-28 the load average decided it, and our own scan
  held it above 8, so `run-step` used the day cap at every run. Card C-28 read the new measure on
  the pod once (`biopy -m biobox.runstep --profile "$PROFILE" --threads 4`, section 1; its numbers are
  in the card's report); no build has run with it yet.
- `sh scripts/pod_pytest.sh <name> --golden` (the golden ring needs the golden set on the pod).
- `git push pod`. Agents never push.
- Single release steps by hand (`bio-box release cluster ...`, `bio-box stage genecall ...`).
  Only `build-stage` ran them.
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
