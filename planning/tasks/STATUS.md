# Task status

The running record of the cards in [`TASKS.md`](TASKS.md). Updated at every merge. Branch
`plan/task-cards` holds all merged work; `main` holds the work up to pull request 2. Each card was done
by an agent on its own branch `task/<ID>`, audited by the orchestrating session, and merged after its
tests passed on the merged tree.

## Handoff, 2026-11-12 14:30 CST (current; the dated sections under "History" are superseded where they differ)

Written so a new session can continue without the conversation that produced this file. The orchestrating session
was restarted at 13:13 CST (Fable 5.2 again from about 14:20). **Agents survived that restart**: the four listed
under "Agents running" kept working in their worktrees. Before restarting any card, list the live agents and look
for a running `pod_pytest.sh <ID>` process; a second agent started in C-21c's worktree by mistake and stopped
itself when it saw the first one's files. Pod jobs (systemd units) keep running in any case.

**For a new session.**
1. `git fetch origin`. The integration worktree is `.claude/worktrees/integration` in the main checkout, on local
   branch `integration`, which tracks `origin/plan/task-cards`: merge there and push with
   `git push origin HEAD:plan/task-cards`. The main checkout itself is on an old branch; do not merge there.
2. Re-arm watchers on the pod jobs below with `scripts/orchestrator/watch_until.sh '<remote test>' '<remote report>'
   [interval]` (one short ssh call per interval; long-held ssh sessions drop with exit 255 when the laptop sleeps
   or the pod closes them). For a full suite that must survive drops: `pod_pytest.sh <name> tests/r0/test_runbook.py`
   to sync, then on the pod `bash /tmp/if2453/ci-tmp/detached_suite.sh <name> 2`, poll `systemctl --user is-active
   'ci-<name>-*'` and `tail -n 1 /tmp/if2453/ci-tmp/<name>-detached/log.*`.
3. Pass `watch_ci.sh` the full sha from `git rev-parse`, never a typed one.

**How the work runs.** One card per agent, in its own worktree (`.claude/worktrees/agent-*`) on branch
`task/<ID>` (measurements: `measure/<ID>`); tests only on the pod (`sh scripts/pod_pytest.sh <ID> [--shards 2]
[paths]`). The orchestrating session audits each report, merges `--no-ff --no-commit`
into `plan/task-cards` in an integration worktree, runs the affected tests on the merged tree with
`env BIOBOX_TREE=<worktree> sh scripts/pod_pytest.sh integration <paths>` (a merged tree identical to the agent's
tested tree needs no rerun), commits, adds the card's row to the table below, pushes, and checks the GitHub
Actions run of the pushed commit (`scripts/orchestrator/watch_ci.sh <worktree> <full sha>`). Sonnet for tier-L
cards, Opus for the rest, at most four agents at once. Long pod jobs run as detached systemd units with
DONE/FAILED markers, because ssh from a sleeping laptop drops (exit 255 from a watcher means the connection,
not the job).

**Branches.** `plan/task-cards` on GitHub (repository `ifinkelstein/bio-box`, public): 105 cards (and the first session of F-02) and the
measurements in the table below; also the pptx deck (`presentations/all-the-data-deck.pptx`, 51 slides) and the
flank-tools review behind card W-25 (`planning/design/W-25_flank_tools_review.md`; W-25 approved). Last
`ci` and `cross-host` confirmed green: `cb6d3b2` (2026-11-12, the C-21c merge). **`ci` is red at `7d00cc1`** (run 43327988144): one C-16 test, `test_c16_a_search_killed_at_its_wall_time_records_its_peak`, compares a scope's `memory.peak` (269 MiB) with the worker's VmHWM (357 MiB), a wrong oracle that only runs where a scope opens (Actions, not the pod) and is flaky there; the C-16 agent is correcting it on `task/C-16` and must show 5 consecutive green runs. `cross-host` is green. Manifest hash
`f6126b00d45e9284d2799c3a699c992c` (since C-21b's new schema; `5b7fceb8…` after C-21a). The
manifest hash includes a hash of every registered schema, so it changes with any schema change (found by C-23). `main` has pull
requests 1 and 2 (87 cards) and is behind; its image passes its own suite (`ghcr.io/ifinkelstein/bio-box-test@sha256:7f97bc7b…`,
`ghcr.io/ifinkelstein/bio-box@sha256:927cf824…`, on the pod as `~/bio-box/images/bio-box{,-test}_80d4c55.sif`),
but the runbook and `scripts/pod_pytest.sh` still use `bio-box-test_74b3e05.sif`.

**Agents running (14:30 CST).**

| Card | Branch | State | To continue |
|---|---|---|---|
| W-25 group a search's loci by their flanks | `task/W-25` | started 2026-11-12 evening; design note first (`planning/design/W-25_flank_tools_review.md` is the method review) | audit, merge; run on `rN.3`'s RNA polymerase search (2,581 loci) |

**On the pod now** (wilkcomp02; account thread cap 40 by day, 56 at night, a total over everything we run).
- The 50 Gbp build, unit `r50-driver` (log `/tmp/if2453/r50/driver.log`, step log `/tmp/if2453/r50/tmp/logs/r50-P6-5.log`;
  state `/tmp/if2453/r50/driver/{DONE,FAILED}`; pause `touch /tmp/if2453/r50/driver/HALT`; `FAILED.1409-C20` there is
  the old C-20 failure, renamed; `FAILED.1057-tokens` the one below). Stages P3 to P6 are ingested and gene-called.
  Release `rM.1` published its clustering at 10:53 CST (*measured*: 45,815,843 proteins, 31,183,876 NR90 and
  8,192,065 NR50 clusters), then **`rM.1.tokens` failed at 10:57 (DuckDB out of memory in the slot join; card
  B-46)**. The driver was relaunched at 11:00 with `SKIP_RM1=1`: it skips `rM.1`, scans Pfam on P3 to P6 now (about
  midday Saturday, *estimate*; each scan launch ends with the same tokens failure after the scan is published,
  which the driver accepts), then waits for `/tmp/if2453/r50/driver/RELEASE_OK`. **B-46 is merged (2026-11-12 14:45): the fix is in `plan/task-cards`.** At 13:07 CST it was scanning P3
  (unit `r50-P3-6`, launched 11:00). **To finish:** sync the fixed tree of `plan/task-cards` to
  `/tmp/if2453/ci/r50` (rsync as `pod_pytest.sh` does) only while the driver is waiting (its log says "scans done;
  waiting for .../RELEASE_OK"), never during a scan, then `touch /tmp/if2453/r50/driver/RELEASE_OK`; the driver
  builds `rM.2` (with C-21b's `dim` and B-45's plan option; at 46 M proteins B-45 changes no cut). `rM.1` can then
  be finished by relaunching its release steps with the fixed code. The driver and its `SKIP_RM1` mode are in
  `scripts/stress/build_stage_pod/r50_driver.sh` (the pod copy in `/tmp/if2453/r50/bin` equals it). Watch: remote
  test `[ -e /tmp/if2453/r50/driver/DONE ] || [ -e /tmp/if2453/r50/driver/FAILED ] || [ -e /tmp/if2453/r10/root/STOP ]
  || ! systemctl --user is-active --quiet r50-driver || grep -q 'waiting for /tmp/if2453/r50/driver/RELEASE_OK'
  /tmp/if2453/r50/driver.log`.
- Measurement 0-08a, the Rfam sample scan (`/tmp/if2453/m08a`, units `m08a-*`, 20 shards done at 13:07 CST; branch
  `measure/0-08a` pushed at `935dcac`). When `DONE` exists: copy `report/2026-11-12_rfam_sample.md` to
  `planning/measurements/` and `0-08a_prok-2_proposed.tsv` to `planning/design/`, write the interpretation, merge.
- Search pack directory (W-22): `/tmp/if2453/bio-box/packs` holds the `rN.3` pack (copied 2026-11-12 from
  `/tmp/if2453/r50/packs/release=rN.3/a334c89ccf048904f297de3b2597af22`; equal `du -s`, 7,573,492 KB; with its
  `pack_verified.json`; `chmod -R go-w` done). The `rN.1` pack (`/tmp/if2453/w19/packs/release=rN.1/71917a4167c7aebe1b7ff44d18a99306`, with its
  verification record) may be copied the same way. Nothing was moved or deleted.
- An orphaned search worker from W-14's `test_runner_death_is_reported` sat idle in a test scope from 2026-11-10
  until the PI stopped it on 2026-11-12: a killed runner leaves its worker alive. Put that in card C-16's brief.
- Data root `/tmp/if2453/r10/root` (46 GB; stages S1, P1 to P6; releases `rN.1a`, `rN.1b`, `rN.1`, `rN.2`, `rN.3`)
  is on local disk, which a reboot clears, and **not yet copied to `/stor`**. Our files on the local disk total
  about 350 to 390 GB of its 1.1 TB (543 GB free), mostly scratch.

**Next, in order.**
1. Merge B-46, C-23, C-21c and W-24 as they report; then C-16 (after C-23), W-15 (after W-24) and W-25 (after
   C-21c and W-24).
2. When `rM.2` is built and verified: copy the data root to `/stor/scratch/Wilke/if2453/all-the-data/` (one
   `rsync`, atomic rename; files under 64 MB only as `small_ok` tables allow, AGENTS §7.3), then clear obsolete
   scratch (old pack copies under `w16`, `w19`; finished measurement directories).
3. C-24 (stage H1, huge phages) once the Pfam scan of P3 to P6 is done (one build in the root at a time); B-42
   after B-45; then PLSDB (no card yet; needs C-21b and C-21c).
4. A pull request bringing `main` up to `plan/task-cards`; switch the runbook and `pod_pytest.sh` to the new image.
5. Finish 0-08a (above). S2C phase 2 (card C-25) waits for the PI.

**Waiting for the PI.** Listed under "Open items for the PI" below.

## History

### State on 2026-11-10 evening (after 17:10 CST)

- Merged since the section below: B-27, E-01, C-12. 66 cards. The pod suite with C-12 merged:
  1,148 passed, 12 skipped. `ci` and `cross-host` were green on the last three pushes
  (`740d049`, `dd6103a`, `94b6f8d`); the fixture's manifest hash is equal on Actions and the pod.
- All agents were cut off by the usage limit in the afternoon and resumed after 17:10 CST.
- Running: the 10 Gbp rebuild (`build/r10`; phages and the 250-genome stage built, the 2,247-genome
  stage in ingest when the agent was cut off), W-03, C-13, B-30 (card amended: one `ingest_plan`
  artifact and `stage ingest` commands). C-13 and B-30 are built on C-12.
- Held: W-13, until W-03 merges (same search files). E-07 waits for PI card 0-08.
- Then: the merge to `main`; then pull the image by digest; then the 50 Gbp build.

### Trial search on `rN.1`, 22:50 CST: a search cannot run on it yet

Record: `planning/measurements/2026-11-10_trial_search_rN1.md` (all *measured*, wilkcomp02, read-only root).
- Current code opens and verifies `rN.1` (60 s), finds T7 by accession, and exports all 10.6 M proteins (88 s).
- **Blocker:** `pack build` holds every region in memory; it hit a 32 GiB cap in the dark-ORF part and
  thrashed, so no pack exists and `search run` refuses. Card W-16.
- **Wrong answer:** `nbr summarize --by nr50` exits 0 with every neighbour `unannotated`, because `rN.1`
  predates `member_cluster` (B-27). Card W-17. The rebuild's releases need the newer tables: card C-17.
- The answer by hand (`hmmsearch` on the export, Pfam `RNA_pol`, 30 bits): 1,957 proteins, 1,931 in phages
  (equal to the 2026-11-09 walk-through) and 26 in prokaryotes; T7 gp1 scores 367.7 bits.
- Ten more rough edges: cards W-18 and W-15.
- Night threads: the PI raised the scan to 12 shards (48 threads) from the P2 stage, leaving 8 for tests;
  driver unit `r10-p2-driver2` replaced `r10-p2-driver` at 22:09 CST.

### State at 01:40 CST, 2026-11-11

- `plan/task-cards`: 81 cards merged; pod suite 1,255 passed, 13 skipped.
- Rebuild: `rN.1a` and `rN.1b` have Pfam. The large stage is scanning under driver `r10-p2-driver3` at 10 night
  shards (40 threads). 12 shards (48 threads) was refused at 23:15: `run-step` allows only the day cap, 40,
  whenever the load average is above 8, which our own test runs cause. Marker `phase2/FAILED.2315-refused-48`.
- Design B-37 (`planning/design/B-37_split_clustering.md`): split clustering equals the workflow row for row
  with one flag; implementing card B-41. No PI decision needed.
- Measurement 0-09a (`planning/measurements/2026-11-11_mmseqs_profile_vs_hmmsearch.md`): MMseqs2 profile
  search is 27 to 95 times cheaper than `hmmsearch` and recovers 50 to 63% of protein-family pairs at 99%
  precision, with profiles sampled from the HMMs (no seed alignments on the pod).
- 50 Gbp inputs (`planning/measurements/2026-11-11_r50_inputs.md`): 11,973 genomes selected, 2,497 present
  and verified, 9,476 to fetch in four stages P3 to P6. **The fetch is not started: no `sources.yaml` exists
  and the NCBI Datasets v2 source is on no approved list (AGENTS §8). Waiting for the PI.**
- Waiting for the PI: the download source; best-complete-per-species (11,973) or GTDB representatives only
  (9,635); adding `Pfam-A.seed` to `refdb.yaml` for a repeat of 0-09a.
- **A real search ran end to end on `rN.1`** (second trial, `planning/measurements/2026-11-11_trial_search_rN1_second.md`,
  branch `trial/search-rN1-2`): `bio-box search hmm` with Pfam `RNA_pol` gave the same 1,957 proteins as
  `hmmsearch` on the export, in 534 s, but the worker peaked at 26.27 GiB. Cards W-19 (memory and start-up,
  running) and W-20 (output wording). Outputs: `/tmp/if2453/trial2/out/rnap` on the pod.
- Running: W-19, B-41 (split clustering, fixture only), B-39.

### State at 20:45 CST, 2026-11-11

- `main` has the task-card work (pull requests 1 and 2); the image built from it passes its own suite on
  Actions and on the pod (1,322 passed, 9 skipped) and is at `~/bio-box/images/bio-box{,-test}_80d4c55.sif`.
  The repository and both container packages are public since 2026-11-11.
- `plan/task-cards`: 89 cards merged. Decisions DEC-25 to DEC-28 recorded today.
- **`rN.3` is built and verified** (`planning/measurements/2026-11-11_build_r50.md`): the current-code release
  over the 10 Gbp stages, with clusters equal to `rN.2` row for row, built by clustering in parts. A search on
  it gives the same 1,957 RNA polymerase hits, now with Pfam names on the neighbours; worker peak 2.51 GiB.
  Its pack is at `/tmp/if2453/r50/packs`.
- **The 50 Gbp build is running** under unit `r50-driver`: P3 is gene-called (one shared contig skipped,
  DEC-26), P4 to P6 follow tonight; then release `rM.1` (about 46.6 M proteins, clustering in 24 parts,
  expected about midday Thursday), then the Pfam scan (about midday Saturday), then `rM.2` (early Sunday).
  All projections are *derived* or *estimate* in that file. State: `/tmp/if2453/r50/driver/{DONE,FAILED}`,
  log `/tmp/if2453/r50/driver.log`; pause with `touch /tmp/if2453/r50/driver/HALT`. A guard writes `STOP` in
  the root if a clustering step passes 280 CPU-minutes.
- Split-clustering estimates miss: NR90 `align` parts +25%, NR50 `align` +36% and over 300 CPU-minutes once
  (313.84). Card C-22 refits them.
- Fetched, verified and archived in `all-the-data/provisional/`: the 50 Gbp genomes, HPGC huge phages (stage
  H1 inputs ready, card C-24), PLSDB (waits for DEC-27's cards), INPHARED and the GTDB metadata.
  `sources.yaml` and `refdb.yaml` were checked against the files on disk.
- Not yet on permanent storage: the data root itself (`/tmp/if2453/r10/root`, growing).
- Running: E-07 (ncRNA producer for DEC-28). Guides: `GUIDE.md`, `AGENT_GUIDE.md`.

### 10 Gbp rebuild, 18:20 CST

- Phase 1 is done: releases `rN.1a`, `rN.1b`, `rN.1` in `/tmp/if2453/r10/root` on the pod, all verified
  (problems `[]`). `rN.1` equals `rL.3` in every gene, protein, IGR and cluster row (12,048,866 genes,
  10,620,684 proteins, 2,966,200 NR50, 8,278,537 NR90; *measured*, branch `build/r10`,
  `planning/measurements/2026-11-10_rebuild_10gbp.md`). No Pfam yet. Built with the tree at `740d049`.
- Phase 2 (Pfam 38.2, about 103 shards, about 443 CPU-hours, *derived*) started at 18:14 CST under unit
  `r10-p2-driver`: 16 threads by day, 32 from 20:00 to 06:30; expected to end about Wednesday midday,
  then about 2.5 h for the `rN.2` clustering. State is under `/tmp/if2453/r10/phase2/`: `driver.log`,
  then `DONE` or `FAILED`; `touch HALT` there stops it before its next launch; launching the driver
  again resumes. Its STOP-based switch of the shard count first runs at 20:00 and is untested.
- 21:24 CST: the 20:00 switch to 8 shards worked (`r10-p2-S1-2` launched at 20:22 with 32 threads); S1's
  scan and release `rN.1a` finished at 21:18 (exit 0, 9 shards); P1 is scanning (`r10-p2-P1-3`).
- `plan/task-cards` at `7a34965`: 76 cards; pod suite 1,227 passed, 12 skipped; `ci` and `cross-host`
  green again after the one red run on `f7c511e` (a C-15 test, corrected).
- Not yet synced to `/stor`. `build/r10` is pushed but not merged into `plan/task-cards`.
- Found: cards B-35 (`partners_of` 12 times slower than the query it replaced), B-36 (clustering
  2,195.9 CPU-minutes in one process; PI decision), C-14 (runbook), C-15 (peak memory in run records).

### State on 2026-11-10 at about 11:30 CST

The session continued after the handoff below was written (its clock time was wrong: it was
written at about 10:05 CST). Since then:
- The suite on `plan/task-cards` with C-11 merged passed on the pod: 1,064 passed, 10 skipped.
- The source review of the paper is merged: `planning/sources_from_s41467-026-79223-z.md`; eleven
  candidate databases, none added to the import plan.
- **A-07 is merged**: its branch had 20 consecutive green `ci` runs. `ci` on `plan/task-cards`
  should now be reliably green; a red run from here on is a new problem, not the known one.
- Running: C-10 (real build on the pod), B-27, E-01, C-12. All four were resumed after the session
  restart; their branches are `task/<ID>`.
- C-10 is merged. The 10 Gbp rebuild with the current code started on the pod (not a card; root
  `/tmp/if2453/r10/root`, units `r10-*`): first without Pfam, as release `rN.1`, so a queryable store
  exists in hours; then the Pfam scan through `build-stage`'s sharded path and release `rN.2`.
- Next: W-13 and W-03 (unblocked by A-07), B-30, then the merge to `main` once several consecutive
  `ci` runs on `plan/task-cards` are green.

### Handoff, 2026-11-10 (written about 10:05 CST)

Written so a new session can continue without the conversation that produced this file.

**Where things are.** All merged work is on branch `plan/task-cards`, pushed to GitHub. Nothing is on
`main`. 61 cards are merged (table below). Rules for agents are in `AGENTS.md`; cards in `TASKS.md`.
How the work is run: one card per agent in its own git worktree on branch `task/<ID>`; tests only on
the pod (`sh scripts/pod_pytest.sh <ID> --shards 8`); the orchestrating session audits each report,
merges `--no-ff` into `plan/task-cards`, runs the suite on the merged tree, commits, updates this
file, pushes, and checks the GitHub Actions result for the pushed commit. Sonnet for tier-L cards,
Opus for the rest, at most four agents at once. Helper scripts for watching CI and pod units are in
`scripts/orchestrator/`.

**The last merge (C-11) was committed before its merged-tree suite run finished.** The run
(`pod_pytest.sh integration --shards 8`) was in flight at handoff. C-11's own branch passed the full
suite (1,054 passed, 10 skipped). Rerun the suite on `plan/task-cards` first.

**Cards whose agents were running at handoff.** Each has a branch; a branch with commits holds real
work, one without holds none yet. Worktrees are under `.claude/worktrees/` in the main checkout
(`git worktree list`).

| Card | Branch | State at handoff | To continue |
|---|---|---|---|
| A-07 CI flakes | `task/A-07`, pushed, complete | Fix done and audited. 16 of the 20 required consecutive `ci` runs on that branch were green, none red; the loop (`scripts/orchestrator/ci_repeat.sh ifinkelstein/bio-box task/A-07 20`) was still going | Check the remaining runs (`gh run list --workflow ci --branch task/A-07`); if all green, merge into `plan/task-cards`. Until then `ci` on `plan/task-cards` can fail on `test_s05_...[oom-oom]`: that is this known cause, not a new break |
| E-01 ncRNA producer | `task/E-01`, has commits | On its pod test runs; no report yet | Read the branch; rerun its tests; finish the card's acceptance and report |
| B-27 lookup tables | `task/B-27`, has commits | Working; no report yet | Same |
| C-10 first real `build-stage` on the pod | `task/C-10`, no commits yet | Building on the pod in `/tmp/if2453/c10` (units `c10-*`); the Pfam scan of the 250-genome stage was running, throttled to 16 threads | Read `/tmp/if2453/c10/logs` and the run records under `/tmp/if2453/c10/root/runs`; the card says what to measure and write |
| Source review of a paper | `docs/natcomm-sources`, no commits yet | Reading `s41467-026-79223-z` to list databases not in `planning/sources_and_sizing.md` | Rerun from the brief in card form below if nothing was committed |

Brief for the source review, since it has no card: read the paper, list every database it uses or
mentions, mark which are not in `planning/sources_and_sizing.md`, `STAGES.md`, `refdb.yaml` or
`planning/search/05_external_modern.md` and `06_infra_inventory.md`; write
`planning/sources_from_s41467-026-79223-z.md` and a "candidates, not yet approved" subsection in
`sources_and_sizing.md`; add nothing to the import plan itself; no network access; sizes as the
paper states them, with page references.

**Next cards, in order.** C-12 (builds respect the hour; lost scratch; publish race in `run-step`),
then W-13 and W-03 (searches under a pod unit; runs outside the data root: both after A-07 merges),
B-30 (ingest by batch), then phase E after E-01. Then the merge of `plan/task-cards` into `main`,
which the PI delegated: it waits for A-07 so that `ci` is reliably green; the merge builds the new
image, and the `gh` token now has `read:packages` to pull it.

**Larger imports** (plan approved by the PI, section below): the 10 Gbp rebuild waits for C-10's
report; launch it with `build-stage --scan-jobs/--scan-cpus` at a total of 40 threads so it is legal
at any hour (C-12 is not merged). The 50 Gbp build follows.

**Decisions waiting for the PI:** cards 0-04, 0-05, 0-07 and the catalog byte-identity exception,
all under "Open items for the PI" below. Decided today: DEC-24 (cluster from scratch).

**On the pod now:** `rL.1`, `rL.2`, `rL.3` on `/stor` (old layout; current code cannot open them);
C-10's build on local disk; test mirrors under `/tmp/if2453/ci`. The account's thread cap is a total
across everything we run: 40 by day, 56 at night.

Last full suite on the integration branch before the C-11 merge: 1,053 passed, 10 skipped
(*measured*, wilkcomp02, 2026-11-10, commit `6225c32`). Older note, kept for the record: Suite on the integration branch after the B-17 merge: 861 passed, 10 skipped; golden ring
(`--golden tests/r3`): 32 passed, 5 skipped
(*measured*, wilkcomp02, 2026-11-10 about 03:00 CST, `sh scripts/pod_pytest.sh integration --shards 8`, in the
pinned image; 3 min 40 s wall for the slowest shard, against 13.5 min in one process). The skips
need samtools, bgzip, Prodigal or the golden set. Tests no longer run on the laptop (AGENTS.md §5.5).

## Merged (105 cards, and measurements)

| Card | Merge | What it did | Audit notes |
|---|---|---|---|
| A-01 | `14acede` | Search-runner liveness by a file lock instead of pid checks | Audit made the probe a shared lock and stopped `cancel` signalling a reused pid. First Linux CI run since: green |
| C-01 | `3bdacfa` | Host profile for wilkcomp02 | The 08:00 to 20:00 day window was my assumption, now labelled so; weekends use the night cap |
| C-02 | `85d4da8` | `bio-box run-step`: scope, nice, thread caps, staging, publish, run record | Audit found two gaps, fixed in C-03 |
| C-03 | `914320d` | Kill and resume for `run-step`; signals forwarded | Checked on the pod that a systemd scope keeps the pid and process group |
| A-02 | `386ed06` | Corrected stage sizes and rates in the plan documents | Audit removed the stored-IGR term from the store model (0.528 GB per Gbp) |
| A-05 | `fc14de0` | Spill-watcher thread race in the K-02 test | |
| B-01 | `106f387` | One streaming Parquet writer (`io.write_table`) | Bytes are independent of input batching |
| B-02 | `2646ace` | `feature_uid` and feature kinds | |
| B-03 | `2f6f942` | Six new registry tables | |
| B-04 | `315de8e` | Persisted genome inputs | |
| B-05 | `3b55cc0` | Gene calling: streaming, parallel, writes feature tables | 7× faster at 8 workers; memory grows 87 MB from 30 to 300 genomes |
| B-06 | `d12fcf7` | IGRs keyed on features | |
| B-07 | `452d573` | Tapes keyed by tape version; slot-bounded row groups; pruning lookup | Audit: one ledger per tape version |
| B-08 | `65854d1` | Stability relation MR-20 | 0 rows change by `feature_uid`; 514 of 527 by ordinal |
| B-09 | `1bac11e` | Token building through DuckDB | Stopped once for usage limits and resumed; its memory measurement was never run |
| B-12 | `c0019e5` | `contig_dim` | |
| B-20 | `f6a1175` | IGR version names its producer set | |
| B-10 | `8b27396` | `token_postings` builder | Narrowed to one posting per slot and level; per-domain postings are card B-21 |
| B-13 | `6e741b4` | Sequence index looked up, not loaded; `.fai` bisected | Open at 10⁶ contigs: 0.001 s and 0.008 MB |
| W-01 | `e7ed141` | Pack build survives genetic code 15 | |
| W-08 | `6b5c0f2` | Recipe shorthand with label; drift-alarm warning; `contig_label` view | |
| B-14 | `258dc11` | Sequence macros without index scans | 1,000-contig fetch at 10⁶ contigs: 3.2 s to 0.25 s |
| B-15 | `f9af5b6` | Batch ingest with bounded memory | 28% slower (sources read twice); found the old ledger depended on input order |
| B-11 | `f331c20` | Neighbor queries over postings; `slot_pairs` removed | Equal to the pair oracle on all 7,735 fixture pairs |
| W-02 | `db24106` | A user's HMM file or seed FASTA as a search seed | The HMMER-binary oracle test has not run (image only) |
| D-01 | `f682baf` | Mutation testing with a two-pass verdict | Workflow written, not executed |
| A-04 | `4ad0e82` | `doctor`: cgroup memory cap and path visibility | |
| B-16 | `c25699d` | Search layer names slots by `feature_uid` | |
| D-03 | `b108be6` | samtools, bgzip, Prodigal and RefSeq checks as golden-ring tests | Only the RefSeq test has executed; Prodigal tolerance set from the measured 4 |
| B-19 | `4972f71` | Streaming ledger; builder tables cut into parts | |
| W-06 | `5f31050` | `nbr summarize`: neighbor families and clusters per anchor window | Same commit adds `scripts/pod_pytest.sh` and fixes a cluster-ledger test that only fails where MMseqs2 exists |
| C-05, W-04 | `0a4ae15` | Search workers start in scopes; worker thread caps; exit and kill labelling | The scope cannot be opened from inside the image: card C-07 |
| W-09 | this merge | Every pyhmmer search in `search/` takes a required thread count | Agent's full pod run on its branch: 830 passed, 10 skipped. `l1/domains.py` `scan` still defaults `cpus=1` (bounded, so left) |
| A-03 | this merge | README for `scripts/stress`; every script prints usage | 11 of 17 scripts had no usage path; README lists which scripts still use the old data model (folded into card B-18) |
| C-07 | this merge | `scripts/pod_run.sh` opens the memory-capped unit on the host; `run-step` and the search runner inherit a cap they find, or refuse when the profile requires one | Probed on the pod: a 1 GiB allocation in a 512 MiB unit was killed by the kernel (`oom_kill 1`); the pod has 1117 GB of swap, so `MemorySwapMax=0` is set everywhere. The search CLI cannot yet require a cap: card W-10 |
| W-10 | this merge | `bio-box search run --profile` (or `$BIOBOX_PROFILE`): thread cap, nice and required memory cap from the host profile; flags may only lower them | Agent's full pod run on its branch: 866 passed, 10 skipped. With a profile, `mem_gb` defaults to the profile's `mem_scope` |
| B-17 | this merge | Legacy keys removed: `slot_protein` table; ordinal, flank-slot and dimension columns; foreign keys and lineage now run through `feature` and `feature_protein`; `LEGACY_ORDINAL_KEYS` is empty | Schema lock regenerated on the pod with `--breaking`. Nested occurrence objects in search-hit exports changed (`ordinal`, no `genome_uid` or `tier`) without a schema-version bump: folded into card W-05 |
| W-05 | this merge | `search hits --format tsv\|faa`, `nbr <run> --format tsv`, `export proteins`; search-hit JSON is `search_hit/2` with `genome_uid` and `tier` joined from `contig_dim` | TSV checked against a brute-force join written in the test. Agent's full pod run on its branch: 866 passed, 10 skipped. Two lookups cannot prune (card B-27); `export proteins` globs when there is no manifest and FASTA keeps `*` (card W-11) |
| B-26 | this merge | Postings measured by column: 35.7 bytes per row, of which `feature_uid` is 17.0; estimate in the code corrected from 20 to 36 | The audit migrated one test constant that mirrors the estimate. Whether to slim the table is a PI decision: card 0-04 |
| B-25 | this merge | A release is defined once by `release_def` (ordered stages); clustering, tokens, postings, `contig_dim` and the manifest read only those stages; new `bio-box release define\|cluster\|tokens\|postings\|dim\|manifest` | The agent first reproduced the pod bug on the old code (all six release files changed when a later stage appeared). Query-time readers (`pack`, `nbr`, worker, `SeqStore`) still read every stage: card B-28 |
| B-23 | this merge | Ingest takes `out=` and can be rerun after a kill under `run-step` (SIGKILL at nine points, byte-identical result); refuses a stage directory that is not clean before reading anything | The card's premise was wrong: on the pod the second hash is about 3 s of 240 s and compression is 182 to 188 s. The audit had the agent restore the pass-2 digest check and test that a same-length edit between passes is refused. Agent's full pod run on the merged state: 901 passed, 10 skipped |
| B-22 | this merge | `l0_manifest` is partitioned by stage; `store_params.json` is a registered artifact, written once; a first and a second stage both publish under `run-step` and equal two direct ingests file for file | 32 other fixture files byte-identical before and after. Agent's full pod run: 903 passed, 10 skipped; golden 32 passed, 5 skipped. `contracts.check_append_only` still skips all small tables: card D-06 |
| B-18 | this merge | Phase-B exit measurements (`planning/measurements/2026-11-10_phase_b_exit.md`); stress scripts moved to the new model; CAMI II driver added; `.github/workflows/cross_host.yml` | Two of three acceptance criteria are **not met**: some catalog macros read every row group (card B-31), and two S5 storage budgets are exceeded (PI card 0-05). Cross-host equality holds on its first run: GitHub Actions and the pod both give manifest hash `ccd621b1...` for commit `f37c366` (fixture store; no clustering or domain step in it). Short-contig results: store open 0.03 s (was 7.2 s), tape lookup 1 of 15 row groups (was 181 of 181) |
| B-29 | this merge | Postings builder fixed: tokens are exploded one row group at a time into flat scratch runs, then sorted, so the sort spills; memory limit unchanged | Cause found by measurement: the `ORDER BY` over a union fed from a materialized CTE did not spill (1.9 GB) when family postings are sparse. Reproduced small (690,000 slots at a scaled limit). On real `rL.3`: old builder fails at 38 s; new one writes 25,822,709 rows in 105 s at 2.4 GB peak RSS, byte-identical to the old builder run with an 8 GB cap. 771 fixture files unchanged |
| W-11 | this merge | `export proteins` refuses without the release manifest; a check of 2,053,315 stored proteins in `rL.2` found no `*` anywhere, so FASTA output needed no change (a fixture test pins it) | Agent's full pod run: 906 passed, 10 skipped. Hit and locus exports still glob without a manifest: card W-12 |
| D-06 | this merge | The append-only check compares every file and artifact, small tables included | Fault injection on `l0_manifest`, `ledger`, `release_def`, `tape_version_def` and `store_params.json`: each rewrite is reported; clean on the two-stage, two-release fixture. Agent's full pod run: 906 passed, 10 skipped |
| B-28 | this merge | `pack build`, `nbr`, the search worker and `SeqStore` read only the release's stages (no default that globs); the clustering ledger is per release | The agent showed the bug on the old code: proteins planted in a stage outside the release were in its pack and its search hits. Other fixture files byte-identical. Stress scripts still call `SeqStore(root)`: card A-06. One release's manifest takes in other releases' files: card B-33 |
| A-06 | this merge | Ten stale `SeqStore(root)` calls in eight stress scripts fixed; a static R0 test parses the scripts and `src` signatures and fails on a stale call | `synth_scale.py` at 10× and `builder_limits.py` ran to completion on the pod; the pilot scripts were checked statically only. Agent's full pod run: 917 passed, 10 skipped |
| D-05 | this merge | `tests/r2/test_scale_layout.py`: 100× synthetic store; row-group limits; lookups that prune asserted at ≤ 2 row groups per key; eight lookups that do not prune are strict expected-failures tied to card B-31 | Runs in 40 to 48 s. At the real row-group limits a 100× store has 1 or 2 row groups, so the test scales the writer's caps by 1/16 while building (test side only) to get ≥ 8; the real limits are checked at 1,000× by `scripts/stress/synth_scale.py`, not in CI. B-31 must remove the eight markers |
| B-33 | this merge | A release's manifest lists its stages' files, its own `release=` files, the tape versions it uses and unpartitioned artifacts; `verify_manifest` ignores other releases | An earlier release's manifest still verifies and rebuilds byte-identical after a later release. Fixture manifest hash unchanged. Stopped once on two tests outside the card: one assertion migrated with authorisation; a latent bug fixed (a manifest with no Parquet files could not generate a catalog) |
| C-08 | this merge | Every builder takes `out=`; the release manifest and the catalog are registered artifacts; `bio-box stage ...` and `bio-box release ...` default `--out` to `$BIOBOX_STAGING` | Each builder as a real `run-step` subprocess equals the direct build file for file, MMseqs2 and HMMER steps included; the 12-step chain gives the same manifest hash; a catalog built in staging embeds the final root. Merged tree on the pod: 983 passed, 10 skipped, 8 xfailed; golden 32 passed, 5 skipped. The size rule will refuse a second release's manifest at the pod's floor: card C-09 |
| C-09 | this merge | Size rule: small-by-nature artifacts (manifest, catalog, `store_params`) are exempt by a registry flag; other files are grouped by their partition tuple without bucket keys | Reproduced first at the pod's 64 MiB floor: the second release's manifest was refused by `run-step`. The audit had bucket keys taken out of the grouping, so tiny files spread across buckets are still refused. Agent's full pod run: 988 passed, 10 skipped, 8 xfailed |
| C-04 | this merge | `bio-box build-stage`: one command runs inputs → ingest → gene calling → IGRs → tapes → domain scans → clustering → tokens → postings → `contig_dim` → manifest → catalog, each as a `run-step` | Equals the direct build file for file on the fixture; a second call starts nothing; resumes after a kill; honours `STOP`; refuses a stage whose ingest or scan would exceed the CPU-minute budget (cards B-30, B-24). Branch `task/C-04-attempt2`. Tested on the fixture only: not yet run on real data on the pod (card C-10) |
| W-12 | this merge | Hit and locus exports read only the release's files (manifest, else its `release_def` stages); an undefined release exits 2 | Shown first: without a manifest an S1 gene appeared in a release over S0. Agent's full pod run: 994 passed, 10 skipped, 8 xfailed. A dead glob branch remains in `nbr._store_files`; card B-27 removes it |
| B-32 | this merge | Gene calling sends short genomes (under 10 kb) to worker processes in batches bounded by bases; long genomes stay in threads | On 100,000 CAMI II contigs: 8 workers 32.4 s against 1 worker about 160 s (5.0×; was 1.12× before). Complete genomes within 3% (fixture) and 11% faster at 8 workers (golden). All five gene-call outputs byte-identical across 96 runs. The profile showed sequence reads under the GIL were the bound, not pickling; two follow-ups in `SeqStore` noted in the measurement file |
| C-06 | this merge | `RUNBOOK.md`: build, pause, resume, verify, search, export, clean up, limits; an R0 test checks every command in it against the real CLI parser and the limits table against the profile | Every command was run on the pod on a fixture-sized root (appendix A lists each with its last line). Not exercised: a reboot, the real 64 MiB floor and 200G cap, real Pfam (card C-10). Found four defects in running searches under a unit: card W-13 |
| B-31 | this merge | Catalog lookups prune: a key through `token_of` reads 1 row group; digest and length macros read no scan; the neighbour query's partner side reads the tape tokens of the anchors' contigs, not the whole token level | 1,000× synthetic, cold: `token_of` postings 0.32 s to 0.09 s. Partner queries are slower cold at that size (family 0.25 s to 0.37 s) but no longer grow with the release; crossover about 1,200× to 1,480× (*derived*, not measured). The two D-05 tests with an impossible bound were replaced by tests of the per-anchor-contig bound; no xfail remains. `partners_of` now needs `open_catalog`. No real-data number: current code cannot open `rL.3`. Follow-ups: card B-34 |
| B-24 | this merge | Clustering streams: Python peak 997 MB at 2.05 M proteins, was 8,921 MB; outputs byte-identical. Domain scan can be planned, scanned in shards and merged, and `build-stage` does that instead of refusing | The 46 GB at 10.6 M proteins and the 410 GB projection were our Python process, not MMseqs2; MMseqs2 itself is at least 5.3 GB at 2.05 M and its growth is unmeasured (PI card 0-06). Shards run one after another and are sized with no margin: card C-11, needed before a real Pfam scan |
| B-34 | this merge | `token_of` prunes (13 row groups to 1 for a literal or one joined key at 1,000×) through footer statistics on `token_dict.key`, no new table; a catalog can be opened several times in one process; a plain DuckDB connection gets an error that names `open_catalog` | Found: the catalog `.duckdb` file was never byte-deterministic: two builds differ in about 65 bytes of DuckDB internals (file identifier, checksums, padding). Its logical content is now tested instead. The manifest, not the catalog, is the record, so this is an exception to AGENTS §4.3 to state, not a defect to fix: listed under open items. A cold `token_of` over 100 spread keys is slower than before at this size (0.059 s against 0.018 s) |
| C-11 | this merge | `build-stage --scan-jobs J --scan-cpus C`: scan shards run J at a time, J × C within the hour's thread allowance; shards planned at 80% of the CPU-minute budget (32.9 M residues for Pfam, *derived*) | Byte-identical to the unsharded scan at J of 1, 2 and 4; a failed shard stops new launches and the next call runs only what is missing. Found and worked around a race between two steps publishing at once. Open: the allowance is fixed at the start of a build, so a night build keeps the night cap by day; lost scratch between shards and merge leaves the build stuck: card C-12 |
| A-07 | this merge | Two intermittent `ci` failures fixed at their causes: the injected-oom test no longer races a 3 s wall clock (a memory-capped scope takes 4.5 to 5.2 s to kill at 4 GB, *measured*); a runner counts as live only if it holds the lock and its process is not exiting (the lock outlives the main thread by 25 to 59 ms, *measured* over 50 kills) | 20 consecutive green `ci` runs on its branch. Also fixed: an oom fault that was not stopped was reported as `oom`; it now ends `failed: that limit is not enforced`. Not established: that GitHub's runner takes the scope path (inferred) |
| C-10 | `2b8e421`+1 | First real `build-stage` run on the pod with the real profile: 1,947 phages and 250 prokaryotes, two releases, both verify; all 24 run records show the inherited 200 G cap | The 250-genome stage equals the earlier ladder build gene for gene (927,920 genes, 918,765 proteins, 764,935 IGRs; all compared). STOP, a repeat call and a mid-step kill behave as on the fixture. Ran the code from before B-24/C-11, so the sharded scan inside `build-stage` was not exercised. Seven defects: card C-13 |
| B-27 | this merge | Lookup tables `contig_label` (per stage, by contig) and `member_cluster` (per release, by level then protein); exports and `nbr summarize` read them | On the real 10 Gbp release a protein's cluster lookup read 1 of 173 row groups in 0.04 s, against 173 of 173 in 0.51 s (*measured*, with the first sort order; the final order reads up to 3). Every export byte-identical. The audit changed the sort to (level, member) so a many-neighbour summary reads one level's half. `member_cluster` is about 860 MB at 10 Gbp. Joined-key lookups on the multi-stage label view still scan |
| E-01 | this merge | ncRNA producer: Infernal `cmsearch` with Rfam at the fixed Z, hits at or above GA become `feature` rows (producer `infernal-1.1.5-rfam-15.1`) and `l1_ncrna` rows; `bio-box stage ncrna` | Oracle: the `cmsearch` binary's own table, parsed separately. A hit across the origin of a circular contig is found once. Byte-identical at 1 and 4 threads. **Cost: 1,319 bases per CPU-second with full Rfam: about 2,100 CPU-hours for 10 Gbp, five times Pfam** (*derived*). Keeps every overlapping hit (seven rRNA models on one operon): PI card 0-08. Not in `build-stage` and not sharded: card E-07 |
| C-12 | this merge | `build-stage` rechecks the hour's thread allowance every 60 s while shards run; two steps publishing at once are serialized by a lock in `run-step`; a shard whose scratch run is lost is rescanned alone (`python -m biobox.l1.domains rescan`); `pod_run.sh` sets the unit's CPU quota to the night cap | Tests inject the clock: a build that crosses from night to day never has more threads in flight than the hour allows. Both publish orders pass on the `copytree` and `rsync` paths. A build with one shard's scratch deleted ends byte-identical. `tests/r2/test_build_stage.py` now takes about 15 minutes on the pod (*measured*, agent's run) |
| W-03 | this merge | The search pack and search runs live outside the data root: packs at `<pack-dir>/release=<R>/<pack_id>/`, runs in `--runs`; with a profile the defaults are `<tmp>/packs` and `<tmp>/search/<user>/{runs,work}`. A pack is never overwritten; a second build is a no-op. `run_index.json` is gone | A full search (pack verify, run, status, hits, `nbr`) against a read-only root leaves every file's sha256 unchanged. A pack, run or work directory under the root is refused. Two pack builds are byte-identical. The fixture pack is 2.1 MB (*measured*); about 3.4 GB at 10 Gbp (*derived*, `planning/measurements/2026-11-10_search_pack_size.md`). `RUNBOOK.md` is not yet updated: card W-13. Old runs and packs under `<root>/search` and `<root>/search_pack` on the pod are no longer read; removing them is the PI's call |
| C-13 | `b23839b` | Seven defects from the first real build: each step has its own scratch directory, removed before a new attempt; `run-step` and `build-stage` lines are timestamped and line-buffered; `pod_run.sh` writes its exit line after a stop and runs Python unbuffered; `bio-box release verify`; Pfam rate constant 7.96 CPU-min per M residues (*measured*, C-10); a fixed ingest term 0.07 CPU-min (*measured*, C-10); a test for `STOP` inside a step | Pod suite 1,172 passed, 12 skipped. Shards are now planned at 30,150,753 residues (*derived*: 300 × 0.8 / 7.96). The `STOP` test passed on first run (it tests the existing rule). Not run in a real pod unit by this card; W-13 confirmed the exit line under systemd |
| W-13 | this merge | A `search run --wait` writes the pending record, so `search status` and `search cancel` work from another shell; `search run` without `--wait` inside a `pod_run.sh` unit is refused (exit 2), because the unit is `Type=simple`, `KillMode=control-group` and systemd kills a detached runner when the main process exits (*measured*, wilkcomp02); `RUNBOOK.md` brought up to date with W-03 | Pod suite 1,183 passed, 12 skipped. On the pod (units `w13-*`): status showed `running`, cancel ended the run `stopped`, and a unit stopped with `systemctl` ended its log with `received SIGTERM` and the exit line. The two long runs used the test hook `--fault timeout` |
| C-15 | this merge | A run record's `peak_rss_bytes` is the step's own, measured by a small helper process that spawns the step, with `peak_rss_source` beside it: `scope_memory_peak` when `run-step` opened a scope for the step alone, else `child_maxrss` | Cause confirmed: Linux carries the launcher's image into the child's `ru_maxrss`. Under a 512 MiB launcher a 192 MiB step and an empty step both recorded about 606 MiB before, 203.0 and 10.5 MiB after (*measured*, wilkcomp02). On the pod (inherited cap) the number is the largest single process of the step, not the sum of processes running at once. **Every peak in a run record written before this merge is the launcher's**, including the 10 Gbp rebuild's. Pod suite 1,186 passed, 12 skipped. Same defect in the search runner: card W-14 |
| C-14 | merge before this one | Runbook: the done line is 12 steps plus one per `--library`; the day shard count is a formula of what else of ours is running under the 40-thread total; the `STOP` rule stated once (a shard still running when `STOP` appears runs to the end, publishes nothing and reruns; shards that ended before are kept); six quoted `python -c` one-liners replaced by a script `rb.py` printed in Appendix B | The one-liners were in §5, §6 and Appendix B, not §0. `rb.py` has not been run on the pod; the runbook says so |
| B-30 | this merge | Ingest as batch processes when the one-process estimate exceeds 80% of the CPU budget: pass-1 shards by source object, a plan, one process per batch (one shard each, DEC-21), a merge; `bio-box stage ingest pass1\|plan\|batch\|merge`; one new artifact `ingest_plan`, not in manifests | Byte-identical to the single-process ingest on the fixture and on the 1,000× stage, where the largest process used 0.248 CPU-min against 3.761 for one process, and the batch path cost 33% more CPU in total (*measured*, wilkcomp02, before the merge with C-13). Small stages build as before. Pod suite with C-14: 1,205 passed, 12 skipped. Open: pass-1 shards are cut at the whole-ingest rate (about 13 times too many); a batch reads whole objects, so it can read far more than `shard_bases`; `stage.py` imports a private helper from `cli.py` |
| B-35 | `17ce0ed` | `profile_query` counts row groups served by the lookup cache (`row_groups_read`, `row_groups_cached`); the many-anchor bound of `partners_of` is stated and tested | The lookup itself was already right: a 10-anchor query reads 7 of 26 `tape_tokens` row groups on `rN.1b` (*measured*); with many anchor contigs it reads all of them, by design of the sort key. Not fixed: 9.5 s for 3,509 anchors on `rN.1`: card B-38. Pod suite 1,207 passed, 12 skipped. This commit also corrects a C-15 test that failed on Actions (`f7c511e`), where `run-step` opens a scope |
| B-21 | merge after `17ce0ed` | One family posting per included domain (`seg_idx` = rank by `env_start`); `token_dict` holds every model on a slot; `tape_tokens.family_tokens` keeps the canonical family; no schema change | `fusions_of` equals a brute-force oracle over the domain table. Cluster tokens and postings unchanged. Family postings about double at 10 Gbp, about 7% more postings in all (*derived*). Not measured at scale: card B-40. `partners_of` still reports canonical families only: card B-39 |
| W-14 | merge after B-21 | A search run's peak memory is the worker's own, through `run-step`'s measuring helper | With a 768 MiB runner: 1,014 MiB before, 322 MiB after (*measured*, wilkcomp02). The helper dies with the runner (PR_SET_PDEATHSIG, Linux), so the worker still notices. A run killed for timeout or memory now has no peak: card C-16 |
| B-36 | this merge | DEC-25: `build-stage` estimates clustering CPU-minutes before launch and refuses above the profile's `cluster_cpu_min_max` (20,000 on wilkcomp02, set by the PI on 2026-11-10; it was 12,000 at the merge); without the key there is no exemption; the cluster record carries `cpu_budget_reason: exempt: DEC-25` | Fit: 284.10 × (n / 2,053,075)^1.244 CPU-minutes, exact at the two upper measured points, 60% high for the phage-only set (*derived*). At 20,000 the cap is reached at about 62.7 M proteins (*derived*), above the 50 Gbp build's roughly 50 M (*estimate*). A stage with no gene-call run record cannot be estimated and is refused |
| W-07 | this merge | `bio-box search hmm --hmm FILE --release R --bits 30 --neighbors 5 --positive ACCESSION:START-END --out DIR`: writes the recipe, the hit table, the FASTA, the neighbour summary and the per-locus tables; always waits; needs the release's pack; `nbr summarize` gains `--out` | On the fixture its files are byte-identical to the separate commands, and its recipe rerun gives the same hits. Works on a release with no domain library (cluster levels only, with a note). The card's old acceptance (`rS1.pilot1`, 1,931 proteins) was replaced: that release cannot be opened by current code. Pod suite 1,235 passed, 12 skipped. Not yet run on a real release. The positive lookup scans all of `contig_label`: card W-15 |
| E-02 | merge after W-07 | G-04 spike-in panel, test only: 5 `cmemit` sequences from each of 20 Rfam families, each one detectable alone at its family's GA, planted in a 68 kb random contig; recovery with an exact binomial interval; skips without Rfam | **97 of 100 recovered, interval (0.915, 0.994); no hit on background** (*measured*, wilkcomp02, Infernal 1.1.5); asserted bound 0.90. The producer's hits equal `cmsearch`'s own. The 3 losses (two tRNA, one group II intron) score above GA alone and are found with `--max`: Infernal's filter pipeline at Z = 1e5 Mb drops them, not the producer. The first panel, with unscreened draws, recovered 67 of 100, all but 3 for weak draws. Family recovery only; subtype assignment waits for PI card 0-08. Rfam for the test is a copy at `/tmp/if2453/golden/rfam/` on the pod |
| W-17 | `79632dc` | On a release without `member_cluster` (built before B-27), `nbr summarize` and the exports read `cluster_member` and say on stderr that it is the slow path; the exports' label columns fall back to `contig_occurrence` the same way; `release verify` warns (exit unchanged) for each lookup table the release lacks | Before: `nr50` keys 0 where the brute-force oracle has 62, with exit 0. The catalog's own `contig_label` view is still empty on such a release (card W-18). Branch suite 1,240 passed, 12 skipped |
| B-38 | this merge | `partners_of` reads each `tape_tokens` row group once per query, in batches of at most 16 row groups on up to 8 threads; memory bound about 300 MB (*derived*) | On `rN.1` (*measured*, wilkcomp02, read-only, 4 CPUs): 3,509 anchors 12.35 s before, 0.52 s after, 0.78 s for the plain join; 10 anchors 0.33 s, 0.14 s, 0.47 s. Answers equal by sha256. No switch to a join: it ties in a middle band and loses elsewhere. Branch suite 1,237 passed, 12 skipped |
| W-18 | merge after B-38 | Six small defects from the trial search: `open_catalog`'s DuckDB temp directory is under the system temp directory, never the data root; a `protein_accession` positive control is refused by the recipe parser (the release records none); `nbr --anchors-file`; `pod_run.sh` logs the unit's CPU and memory totals before the exit line; `bio-box find --accession ACC[:START-END] | --protein UID`, shared with `search hmm --positive`; the catalog's `contig_label` view falls back to `contig_occurrence` on a release without label files | One existing test migrated (the `protein_accession` case of a W-08 test, now a refusal). The `usage` line of `pod_run.sh` is tested with fakes only: check one real unit log. Merge conflict with W-17 in one test file, both sides appended; both kept. Pod suite 1,255 passed, 13 skipped |
| W-16 | this merge | `pack build` and `pack verify` stream: regions are written through a stable external merge sort and shards as they fill; progress lines on stderr; a signal removes the build's temporary directories; the missing-pack message no longer implies `--model` is required | Peak memory no longer grows with regions: 330 MiB at 10× and at 200× (it grew 175 MiB before). Fixture pack byte-identical, same pack id. **The pack of `rN.1` built under a 32 GiB cap: 750 MiB peak, 5,282 s wall, 86.0 CPU-minutes, 5.51 GB, verified** (*measured*, wilkcomp02, unit `w16-pack`); it is at `/tmp/if2453/w16/packs`. The pack is 32% over the earlier size estimate (dark ORFs twice as dense as on the fixture). Pod suite 1,264 passed, 13 skipped |
| B-41 | `012b4f0` | Clustering in parts (design B-37): `prefilter`, `rescorediagonal` and `align` run as processes over ranges of sorted query keys, joined in key order; one `run-step` per call or part; `bio-box release cluster plan|step|finish|redo`; one artifact `cluster_plan`; used only when the estimate exceeds 80% of the budget | Builds at 1, 2 and 4 parts equal the single process byte for byte on the fixture and a constructed set; without `--add-self-matches 1` the intermediates differ; a lost scratch database is rebuilt to its digest. A shared prefilter index works only on the sensitive NR50 calls and crashes MMseqs2 18.8cc5c on the fast ones, so each part builds its own: every process fits up to about 68 M proteins (*derived*); card B-43. **Not yet run on a real release, and the DEC-25 exemption is still in place: card B-42.** Pod suite 1,289 passed, 13 skipped |
| W-19 | this merge | A search streams the target in blocks and keeps only hits; the positive control is looked up by key; the pack is verified once at build (a record beside its manifest, pack id unchanged); a search writes nothing into the pack directory; the pending record is written before the long work | **The same search on `rN.1`: worker peak 26.27 GiB before, 1.34 GiB after (0.49 GiB while searching); runner 5.65 GiB, 1.34 GiB; wall 534 s, 294 s; first `search status` 160 s, 10 s; the same 1,957 proteins, every output byte-identical** (*measured*, wilkcomp02, `planning/measurements/2026-11-11_search_memory_rN1.md`). Memory no longer grows with the target (+1.2 MiB from 10× to 100×). The search passes are slower (234 s against 121 s; one thread reads each block): card W-21 |
| W-20 | merge after W-19 | Output wording and small defects from the second trial: strand is `+1`/`-1` in every exported file; `search hmm` writes a `README.txt` (0-based half-open coordinates, every column, offsets relative to the anchor's strand); `tier` is a column only when some gene has one (no tiering step exists yet); the missing-table note prints once per command; the hit table says how many rows it shows; `nbr summarize` drops the pfam level with a note on a release without a domain library; a closed run's pending record is removed; `pod_run.sh`'s usage line works | The usage line was checked on a real unit (`w20-usage`: MemoryPeak and CPUUsageNSec printed; systemd had been expanding `${...}` in the unit's command). `dist_to_contig_end_bp` of 0 is not a defect: it marks a partial gene at a linear contig's end; documented. Five existing tests migrated for the new wording. Its suite ran on this same tree: 1,306 passed, 13 skipped |
| W-21 | this merge | The search reads the next target block in a background thread (at most two blocks alive); `l1_domains` is read one row group at a time for the proteins asked for, not whole | The same search on `rN.1`, on a loaded pod: wall 244.6 s (294 s after W-19), search passes about 190 s (234 s), about 2.6 of 4 CPUs busy (2); same 1,957 proteins, scores and order (*measured*, wilkcomp02, unit `w21-hmm`). **The worker's recorded peak rose from 1.34 to 1.78 GiB (+33%)**: about 230 MiB is the second block; the rest appears after the search, cause not established. It does not grow with the target (memory test: +1.3 MiB from 10× to 100×). Still short of the 121 s the search passes took before W-19. Its suite ran on this same tree: 1,313 passed, 13 skipped |
| B-39 | merge after W-21 | `tape_tokens` gains a nullable last column, `domain_family_tokens`: for each slot the tokens of every family with an included domain, in the order of the postings' `seg_idx`. `partners_of` reads it at the family level, so a neighbour's non-canonical domains count; cluster levels unchanged. Schema lock regenerated by `scripts/update_schema_lock.py` on the pod (one appended entry) | The postings could not serve this: one slot's family postings touch every family-level row group (9 of 9, *measured*). Equal to a brute-force oracle over the domain table; a 10-anchor query still reads 8 of 11 row groups; read-once and the timing test hold. `tape_tokens` grows 7.0% on the fixture and 14.3% on the 100× store (*measured*). A release without the column answers from the canonical family with one note. **The cross-host fixture hash changes with this card.** Pod suite, run detached: 1,318 passed, 13 skipped |
| C-20 | `7a4406a` | DEC-26: a contig shared across stages belongs to the stage that stored it first. Gene calling and the ncRNA producer read earlier stages' sequences (training still sees the whole genome) and write nothing for a contig another stage owns, recording the skip in the ledger; `dim` keeps the owning stage's genome and lists the later row as skipped; a release without the owning stage is refused, naming it | Found when stage P3 of the 50 Gbp build failed at gene calling. 4 of 54,265 contigs are shared across stages (*measured*); the same bytes get different calls in different genome contexts (8 against 6 genes; 197 against 198). No duplicated key in any table a release over both stages reads. Which stage holds a shared contig's calls depends on build order. Branch suite 1,333 passed, 13 skipped. Membership in several genomes: DEC-27 and cards C-21a to C-21f |
| B-43 | this merge | Split clustering shares one prebuilt prefilter index per call where MMseqs2 allows it (the NR50 steps 1 and 2); the choice is a table of call settings, tested against the module list | MMseqs2 18.8cc5c aborts whenever the target is an index and `--diag-score` is 0, at every setting tried, so the NR90 step and NR50 step 0 still rebuild per part. Results with the index equal the rebuilt ones at 1, 2 and 4 parts and up to 1 M targets. At 50 M proteins: 18,910 CPU-minutes rebuilt, 16,918 shared (*derived*). Every process fits the budget up to about 87 M proteins (*derived*; the 68 M in card B-41 took every protein as a target). The index's bytes differ between builds; parts are checked by their output digests. 18.8cc5c is the latest release; the PI deferred any change (card B-44). Branch suite 1,315 passed, 13 skipped |
| E-07 | merge after B-43 | DEC-28's ncRNA producer: the frozen set `prok-1` (588 Rfam 15.1 models, 23 scanned on whole contigs), intergenic regions with a 200-base margin and `--notrunc`, hits mapped back to contig coordinates, the scan in shards, and `build-stage --ncrna-cm` | On the 8 genomes of measurement 0-09c: 1,027 of 1,027 expected hits, none extra, 254 families; 254.5 CPU-hours per 10 Gbp (*derived* from 27.65 measured CPU-minutes), +2.9% against the projected 247.4; about 2 GB per scan process. Spike-in panel unchanged at 97 of 100. Sharded equals unsharded at 1, 2 and 4 shards; with and without `--ncrna-cm` the gene-call tables are byte-identical. Three existing ncRNA tests migrated to pass a set. Follow-ups in card E-08: the provenance check rejects the new producer name; shard records sit in the `domain_scan` artifact; cost per base differed 1.5 times between shards. Branch suite 1,380 passed, 11 skipped |
| C-22 | `538c87d` | Split-clustering CPU estimates refitted: one rate and one growth exponent per call, from the two real runs (2.05 M and 10.6 M proteins), on each call's own input, not a share of the release total | Every part of `rN.3` within 1% of its record (a calibration); at 2.05 M every call of 2 CPU-minutes or more within 10%. Whole-release and per-call estimates agree within 3%. Above 10.6 M everything is extrapolated from two samples: at 50 M the NR90 align call needs about 26 parts and the total is about 18,100 CPU-minutes (*derived*); the NR90 prefilter cannot be cut under the budget from about 74 M proteins (card B-44). The 50 Gbp build was launched before this, with 24 parts and a guard at 280 CPU-minutes. Branch suite 1,352 passed, 13 skipped |
| C-19 | `aa5d1a4` | `release verify` covers what the release uses: a file a later scan added to its stages, for a library, producer or version it does not list, no longer fails it; a warning (exit unchanged) when a library has domain files for some of the release's stages only | Decided from the registry's partition keys, per stage; no manifest change. Still reported: a listed file missing or changed, a second part beside a listed file, any file no table claims. Rerun on the real root after the merge (read-only, 2026-11-11): `rN.1a`, `rN.1b`, `rN.1`, `rN.2` and `rN.3` all verify, exit 0, with the expected warnings about lookup tables the older releases lack. Branch suite 1,348 passed, 13 skipped |
| E-08 | `4322400` | Clan competition for ncRNA hits (DEC-28): among hits of one Rfam clan on the same strand that share a base, the best bit score is kept, greedily; every hit stays in `l1_ncrna` with `clan` and `clan_winner`; only winners become features and go on the tape. Producer name gains `-clan1`. Also: the provenance check accepts the new producer names; shard records move to an `ncrna_scan` artifact (E-07's are still read); the ledger records the gene-call version behind the regions; shards are cut from the slower measured rate (131.6 Mbp) | On the 8 genomes of 0-09c: 1,027 of 1,027 hits found, 52 lose, 975 features; all 25 SSU and all 25 LSU rRNA loci end with one kept feature. Spike-in recovery unchanged at 97 of 100. Tested against `cmscan --oclan` on constructed cases; the greedy rule differs from `cmscan`'s only in chains of four or more overlapping hits (none seen). Two real features lose by Rfam's own clans (small SRP to large SRP in *B. subtilis*; GlmY to GlmZ in *E. coli*): for the PI to confirm. 261.0 CPU-hours per 10 Gbp (*derived*). Branch suite 1,427 passed, 11 skipped |
| C-21a | `dcc726d` | DEC-27, first part: `genome_input` lists a contig once per genome that holds it, within a stage and across stages (key `contig_uid, genome_uid, stage`). Within a stage the genome with the least `genome_uid` writes a shared contig's gene calls and ncRNA hits; across stages the owning stage does (DEC-26). Every genome is still called whole; the other genomes' rows are dropped before writing, and each skip is in the result and the ledger (`owned_by:<stage>/<genome_uid>`). `contig_dim` keeps one row per contig, the owner's. `make_inputs.py` lists every record; `--drop-repeats` reproduces the stages already made. Also fixed: the tape join on raw `genome_input` rows put a shared contig's features on its tape once per genome (26 extra slots and duplicated posting keys on the fixture) | Same stage built with its genomes in two orders gives identical files; no key duplicated in any table a release reads. Schema lock unchanged (the lock holds columns only). Membership view, `genome_dim` and `n_genomes` are card C-21b; counts by genome card C-21c. Branch suite 1,421 passed, 17 skipped |
| W-23 | `3ff6bc4` | Commands a user lacked: `bio-box releases` (stages, genomes, contigs, proteins, libraries, lookup tables, pack), `bio-box seq --protein UID` and `--contig ACC:START-END` (one sequence by key), `bio-box sql FILE` (read-only, SELECT only; DuckDB itself refuses writes, ATTACH, INSTALL and other files: external access off, allowed paths the catalog's files only), `hits.tsv` gains `stage` and `kind` (appended), `find` and `--positive` take the GenBank form `ACC:START..END` and print both forms, one line per output file with its row count. Both guides use the new commands | `profile_query` gains `profile_dir` and returns `columns`. `find` gains `span` and `genbank_span` (W-18 test migrated); the W-05 oracle gains the two columns. `releases` counts genomes from `contig_dim`, which under C-21a misses a genome whose every contig another genome owns: card C-21c. Merged-tree tests 230 passed; branch suite 1,415 passed, 13 skipped |
| 0-08b | `3ff6bc4` | Measurement: the ncRNA producer on phage contigs, and real tRNAs lost to Infernal's filters (`planning/measurements/2026-11-12_ncrna_set.md`) | On 826 S1 phage contigs (47.0 Mbp) the producer costs 165.2 CPU-hours per 10 Gbp, 35% below the prokaryote rate: same intergenic records per base, and rRNA models cost little where there are no operons; cutting shards at 304.0 is safe for phage stages (about 1.8 times conservative). RF00005 with filters off finds 6 tRNA loci of 239 the default misses (2.5%, 95% CI 0.9 to 5.4%), 5 of them in the jumbo phage LAK-A1; it costs about 10,000 times more per base, so the remedy is a dedicated tRNA caller (M3, needs the PI's environment change), not `--max`. M2 cut short on the cost rule (estimate missed 2.3 times); M4 already done by E-07 and E-08; metagenome fragments not measured (no such stage) |
| C-21f | `fb443cc` | Measurement: context-free gene calling (policy p3) against today's calls on stage P1 (`planning/measurements/2026-11-12_context_free_genecall.md`); comparison with today's calls only, no truth set | Control: re-calling P1 under today's policy reproduced 927,920 of 927,920 calls. Under p3: chromosomes 99.78% identical; other contigs of 100 kb or more trained alone 96.39% (median 91% at 100 to 200 kb); contigs under 100 kb in meta mode 86.14%; all of P1 99.47%, 0.13% lost, 0.19% gained. Lost genes are short (median 228 nt) and 19.3% carry a Pfam hit (237 genes, 0.026% of the stage). CPU +1.06%. P1 has no viral, plasmid-kind or circular contigs, so those cases are untested. For the PI: whether calls become a function of the contig alone (DEC-27); judging plasmids needs a reference-CDS source in `sources.yaml` |
| B-45 | `eb58e99` | DEC-29 in code: a prefilter whose parts would each rebuild a k-mer table at or over the part budget runs as one part, the workflow's own call, when the plan holds the profile's `cluster_cpu_min_max` (new option `release cluster plan --one-process-cpu-min-max`, passed by `build-stage` only when the profile sets the key); its estimate (work plus table) is checked by the cut and again before launch (`one_process_refusal`), and its run record says `exempt: DEC-29`. In parts and with the key, the release-level check compares those calls with the cap instead of the one-process estimate T(n). Without the key, refusals are unchanged except that they name DEC-29 | Plans without the key are byte for byte as before, so builds below about 68 M proteins are unchanged. At 95 M the NR90 prefilter becomes one process of about 3,945 CPU-minutes, total about 37,963 (*derived*, extrapolated from C-22's model). Open: without the key the T(n) check still applies in parts (left for B-42 (2)). Merged-tree r0 tests 186 passed; branch suite 1,495 passed, 17 skipped |
| C-21b | `770f886` | DEC-27, second part: a `genome_member` view over the release's stages' `genome_input` (membership depends on the set of stages, not their order) with a `genomes_of(list)` macro; a new table `genome_dim` (one row per genome per release: contigs, bases, kind, genetic code, quality, tier, taxonomy; an attribute null where the genome's stages disagree); `contig_dim` gains `n_genomes` (appended, nullable; null in older files) and its `genome_uid` is the least of the contig's genomes in the release. Schema lock regenerated by `scripts/update_schema_lock.py` | Pruning (*measured*, a 20-stage store of 4,059 row groups): one key reads 19 row groups (one per stage file), 10 joined keys 183, 100 joined keys 1,463: the footer bound. A table macro called once per joined row cannot prune in this DuckDB (read 4,059 of 4,059), so `genomes_of` takes a list, as `tapes_of` does. Open for the PI: a lookup reads one row group per stage; a per-release copy sorted by contig would make it one. `dim`'s result reports `shared_contigs`, no longer `skipped`. Changes the fixture's manifest hash (a new schema). Merged-tree suite 1,519 passed, 17 skipped |
| W-22 | `8664cf2` | One pack directory per host and no writes into it: `pack verify`, and the full verify a no-op `pack build` runs, put their sort scratch under the profile's `<tmp>/search/<user>/work` (system temp without a profile), never beside the pack; only a build that makes a new pack writes into the pack directory, and one that cannot exits 2 with a refusal naming the directory. The guides and the runbook find packs through `--profile` alone, at `<tmp>/packs` (`/tmp/if2453/bio-box/packs` on wilkcomp02, which now holds `rN.3`'s pack) | No new profile key: `profile.load` refuses unknown keys and `profile.py` was not in the card. The profile-only search already worked since W-03, so that test pins it. The read-only tests lower `pack.RUN_ROWS` so the fixture writes a scratch run at all. Left for W-24: the silent fallback to the gene-call tape version. Merged-tree tests 401 passed; branch suite 1,470 passed, 17 skipped |
| B-46 | `228d507` | The tokens step stays inside its DuckDB limit at any size: slots, cluster rows and domain rows are split into protein ranges by hash (the count from the input sizes and the limit), keys and tokens are joined per range, cluster keys are numbered from a sorted stream, and the per-tape lists are built from one sorted stream. The cause: one LEFT JOIN of 48.4 M slots to 45.8 M NR50 cluster rows cannot spill in DuckDB 2.0.0 (hash table 1,860 MiB at a 2 GB limit), while aggregates and sorts do | Output unchanged: byte-identical on the fixture and on `rN.3`'s published `token_dict` and `tape_tokens`; the bytes do not depend on the number of ranges (tested). On `rM.1`'s real inputs into scratch (*measured*, wilkcomp02, 2026-11-12, process peak RSS; the cgroup peaks of 20 to 26 GB are page cache): tokens 1,333.8 s, 19.2 CPU-min, 2,258 MiB (17,841 family, 8,192,065 NR50, 31,183,876 NR90 tokens); with 8 ranges instead of the 64 chosen it still fits (2,332 MiB): at least an 8-fold margin (*derived*). The later steps then ran clean: postings 781.2 s, 11.7 CPU-min, 2,337 MiB, 115,248,867 rows; dim 1.9 s (54,265 contigs, 39,224 genomes, 4 shared contigs); manifest 743.3 s, 2,327 MiB; catalog 24.2 s. `dim` reports 0 circular contigs for `rM.1`: to check. Scratch `/tmp/if2453/b46` can be cleared. Branch suite 1,522 passed, 17 skipped; merged-tree tests 240 passed |
| W-24 | `dda23ce` | Silent fallbacks in the search readers become refusals or notes: (1) `search hits` and `nbr` on a run whose pack is not found exit 2 naming the run, its pack id, its release and the directory looked in, instead of falling back to the gene-call tape version; (2) `nbr` on a run with a `--release` other than the one it ran against exits 2; (3) a family lookup on a release built before card B-21 prints one note that only canonical families are indexed (recognised without a new marker: `tape_tokens.domain_family_tokens` present, else a postings row group with `seg_idx` above 0, else the domain table decides); (4) `neighbors.tsv` drops `cluster_uid` at the cluster levels, where `key` is the cluster uid (the JSON keeps it) | The W-07 byte-identity test needed no change; one W-06 unit test migrated. `GUIDE.md` §6 and `AGENT_GUIDE.md` §4 updated with the merge. Not run on the real `rN.2` (it should now print the B-21 note). Merged-tree tests 339 passed; branch suite 1,535 passed, 17 skipped |
| C-23 | `f13634d` | A release records the code that built it: every run record gains `code_id` (sha256 over the `biobox` package's files), `code_version` and `code_commit` (the git commit when the tree is a checkout); `catalog_info.built_by` lists the distinct identities of the release's run records, with the count of steps built before the card; `release verify` prints one `built by:` line | The manifest stays a function of the data: run records live under `runs/`, outside every manifest. Tested without a pinned hash (the manifest hash covers every registered schema's hash, so a pin breaks at each schema change): the fixture is built from three package copies (the same code at two paths; one with an extra file and another commit), each through a real `run-step`, and the run records' identities differ while the manifests and listed files are byte-identical. The pod image has no `git`: the commit path is tested with a stand-in. Releases built before the card print `built by: unknown`. Merged-tree tests 296 passed, 3 skipped; branch suite 1,532 passed, 18 skipped |
| W-15 | `132230f` | Measurement only (the card's stop condition): resolving a positive control by accession is cheap on today's releases, so no lookup table was added (`planning/measurements/2026-11-12_positive_lookup_rN3.md`) | On `rN.3` (*measured*, wilkcomp02, 2026-11-12): the lookup reads 3 of 3 row groups (one per stage), 832,975 bytes, in 0.034 to 0.047 s; `search hmm --positive`'s resolution 0.15 to 0.31 s; `find --accession` 2.06 to 2.67 s end to end, 1.98 s of it Python start-up. The scan crosses 1 s at about 1 M label rows (*derived*); `rM.2` has about 0.2 M. Found: `rN.3` still uses the contig-name fallback (its stages were ingested before B-27), and a B-27 release would resolve only a contig's first record accession: card W-26. `GUIDE.md` §1 and §9 corrected with the merge. Branch suite 1,538 passed, 17 skipped; merged-tree r0 tests 195 passed |
| D-04 | `f938b13` | Cross-host equality as a test (`tests/r3/test_cross_host.py`, plan test D-05, in `gates.yaml` on the S0 warning list): it builds the fixture, and compares its manifest hash with the one another host recorded in the file `BIOBOX_PEER_MANIFEST` names; without the variable it skips, naming it | The peer is the GitHub Actions `cross-host` job, since the laptop no longer runs tests: `gh run view <run> --log`, the `"manifest_hash"` line into a file on the pod. A wrong peer hash fails; the right one passes (pod and Actions both `f6126b00d45e9284d2799c3a699c992c` at `132230f`, *measured*). `scripts/pod_pytest.sh` now passes the variable into the image (added with the merge). Merged-tree tests 197 passed, none skipped, with the peer set |
| C-21c | `cb6d3b2` | DEC-27, third part: counts by genome include every genome that carries a contig; counts by contig, locus and neighbour do not fan out. A new macro `q1_genomes` counts genomes through `genome_member`; `q1_count`'s tier filter counts a locus once in each tier its genomes have; `nbr summarize`'s `n_genomes` and `bio-box releases`' genome count go through the view; `hits.tsv`, the JSON occurrences and each locus table's facts line gain `n_genomes` and `genome_uids` (appended, on every release, as ruled) | Tested on two copies of the fixture that differ by one genome repeating two existing contigs: postings, `q1_count`, `partners_of` and the locus and neighbour counts are unchanged, every genome count is up by one. On releases without shared contigs every output equals the old code's once the appended fields are dropped (three roots, old code against new in one process). `q1_genomes` reads at most 2 `genome_input` row groups per contig of 188 or more (*measured* in the pruning test). `genome_basis` keeps its old rule. Not added: an `n_genomes` column on `partners_of` (an existing test selects its three columns). `open_catalog` re-creates the four macros older catalogs hold in an older form in one place (resolved with W-24's). Guides updated with the merge. Merged-tree tests 248 passed; branch suite 1,561 passed, 18 skipped |
| C-16 | `316d724` | A killed step or search keeps a peak: the measuring helper writes the peak as the command runs (every 0.2 s, when it has grown: the scope's `memory.peak`, else the largest VmHWM of the command's process tree) and once more at the end, so a kill of the helper leaves its last sample; the record names the source (`scope_memory_peak`, `child_maxrss`, or the `_sampled` forms). A failed, killed or signalled step writes `<tmp>/staging/<stage>.<step>.failed.json` (not under `runs/`, which `built_by` reads), removed by the next attempt; a search's source goes into `limits_json` (no schema change). The helper now starts the command with fork and exec under PR_SET_PDEATHSIG, so a command dies with its helper | The orphan: a worker still importing when its runner died took its new parent as the one to watch and lived on (reproduced with a start-up gate); fixed by the helper change. The 36-hour orphan seen on the pod most likely predates W-14's own fix (*inferred*). *Measured* on the pod: a step killed with its group records 265.7 MiB against 256 MiB written; a search killed at its wall time 323.9 MiB, equal to the worker's VmHWM. The two real OOM-killer tests need a scope and skip inside the image; on Actions they are inferred to have passed from the progress characters, with no printed peak. Left: a step still outlives a SIGKILLed `run-step` (on the pod the unit's cgroup ends it). Merged-tree tests 451 passed, 5 skipped; branch suite 1,558 passed, 20 skipped |
| C-25 (1) | `8b563f2` | Measurement: the S2C selection reproduced from the GTDB R226 metadata without an ANI step (`planning/measurements/2026-11-12_s2c_selection.md`; script `scripts/stress/build_stage_pod/c25_select.py`; lists on the pod under `/tmp/if2453/c25/out`, rebuilt byte for byte by the script) | *Measured*: 49,324 complete genomes in 11,973 species; S2C 25,458 genomes and 106,409,469,843 bp, the plan's numbers exactly; every phase-1 genome is rank 1 of its species; remainder 13,485 genomes, 55,910,526,156 bp, 51,597,124 GTDB proteins (Pseudomonadota 30.23 Gbp). The ANI step can act only in the 243 species with more than 20 complete genomes: at most 3,368 selected genomes could differ (*derived*; 25% of the remainder); in 102 of them the cut falls inside a quality tie and 485 genomes were chosen by accession alone. Five species depend on genomes NCBI did not serve. For the PI: skani or not |
| F-02 (1) | this merge | `synthworld` backbone (`src/biobox/synthworld/`): from a preset and a seed, the species tree and per genome the realized gene order, strands, families, copies and groups with their truth, on ZOMBI2 0.54.0 (`zombi2==0.54.0` added with `uv add`, approved by the PI; the lock gains only `zombi2` and `tqdm`). DEC-15's conditions are enforced and tested: grouped fates for scored units, replacement refused where groups exist, the root scrambled on a stem of synthworld's own and checked at every build | *Measured* (wilkcomp02, 2026-11-12, `micro`, seed 20261112): 46 genomes, 1,509 genes each on average (1,393 to 1,708), 69,411 adjacency rows, 17.5 to 19.6 s against a 60 s limit; grouped fates break 0 of 2,943 units, per-family fates 264 of 3,842. The image has no ZOMBI2: pod tests import it from `/tmp/if2453/pydeps_f02` (the locked wheel, hash-checked) and skip without it; Actions installs it from the lock. For the PI: the `micro` strength grid uses member transfer rate {0.05, 0.15, 0.30} by together {0.3, 0.7, 1.0}, since a per-family loss axis is not available for adjacent groups. Remaining sessions and seven spec corrections are listed on the card. Merged-tree r0, r1 and backbone tests pass |
| DEC-29 | `c3dc901` | The PI chose a CPU exemption for clustering calls that cannot be split under the budget (above about 74 M proteins), over waiting for an MMseqs2 fix or incremental clustering. Card B-45 implements it; card B-42 narrows DEC-25 to it instead of removing it | |

Decisions DEC-18 to DEC-23 were approved by the PI on 2026-11-09 and are in `planning/DECISIONS.md`.

## CI was not reliably green (found 2026-11-10; resolved by A-07, merged 2026-11-10)

GitHub Actions `ci` failed on 7 of the last 12 pushes to `plan/task-cards` while the pod suite passed on
the same commits: two search tests fail intermittently there (card A-07). The orchestrating session
checked CI only after the first push and missed this for a day. The merge to `main` waits for A-07.
The `gh` token has `read:packages` since 2026-11-10, so the image can be pulled once `main` builds it.

## Plan for larger imports (approved by the PI, 2026-11-10)

1. When card C-10 reports a clean real run of `build-stage`: rebuild the 10 Gbp set (the phages and
   2,497 prokaryotes of `rL.3`) with the current code, under the night cap, on the pod's local disk,
   then sync to the provisional area on `/stor`. It replaces `rL.3`, which current code cannot open.
2. When card B-24 merges: build all 11,973 complete GTDB species representatives (about 50 Gbp,
   *derived*: 4.25 Mbp per genome from the 2,497 already fetched) as several stages of about 10 Gbp.
   Expected wall time about 2.5 days, of which Pfam about 1.5 (*derived* from rates in
   `planning/measurements/2026-11-10_import_ladder/`).
The thread cap is the sum over all our processes: 40 by day, 56 at night.

## In progress, and next

See the handoff at the top. Cards written and not started, besides those it names: C-21c, C-16, W-15, W-24,
C-24, B-42, C-25 (S2C phase 2; waits for the PI), B-44 and C-18 (deferred by the PI), E-06 (`rS0.1`, the first
signed release: it locks the schemas, AGENTS §3.8), phase F (statistics), phase G. B-40 is done in effect by the
`rN.3` build.

## Not yet executed anywhere

These exist and pass their laptop-runnable parts, but their real run needs the image:
the clustering tests; the samtools slice and index tests; the HMMER-binary oracles (two tests);
the golden pairs test over postings; the Prodigal test (also needs `prodigal=2.6.3` added to
`container/environment.yml`, a PI decision); the mutation workflow; the systemd path of `run-step`
and of the search worker (only their command lines are unit-tested; the scope behaviour itself was
probed by hand on the pod).

## Open items for the PI

Added 2026-11-12, recommendation first:
- Card C-25: skani for the S2C selection. Without an ANI step the 20 genomes kept in each of 243 heavily sequenced species are
  the best by CheckM2 and then by accession, which can be near-identical isolates of one project (485 genomes chosen by
  accession alone; at most 3,368 could differ, *derived*). Recommendation: add skani to `container/environment.yml` (yours)
  and select for diversity, since strain diversity is why S2C takes 20 per species. Also: retry the 7 genomes NCBI did not
  serve (5 species have no other complete genome)?
- Card F-02: accept the `micro` preset's strength grid (member transfer rate by together fraction) in place of the spec's
  per-family loss axis. Recommendation: accept; it is what ZOMBI2 can vary for adjacent groups.
- Done 2026-11-12: the PI approved card W-25 (group a search's loci by their flanks; review in
  `planning/design/W-25_flank_tools_review.md`). It starts after C-21c and W-24 merge (same files).
- Card C-21b: a contig's genomes are looked up at one row group per stage file (19 at 20 stages, *measured*); a
  per-release copy of the membership sorted by contig would make it one. Recommendation: not yet, until a search
  needs it.
- Card W-22: the host's pack directory is the profile's `<tmp>/packs`; an explicit `packs:` key in the profile needs
  a change to `profile.py`. Recommendation: not needed now.
- Card C-21f: whether gene calls become a function of the contig alone (DEC-27 deferred it to this measurement).
  Under policy p3, 99.47% of P1's calls are unchanged (chromosomes 99.78%; contigs under 100 kb 86.14%), at +1.06%
  CPU. Recommendation: keep today's calls for now and decide before S2C phase 2, when shared plasmids become common;
  judging which calls are better needs a reference-CDS source in `sources.yaml`.
- Card E-08: two real features lose by Rfam's own clans (small SRP to large SRP in *B. subtilis*; GlmY to GlmZ in
  *E. coli*), and the greedy rule differs from `cmscan`'s only in chains of four or more overlapping hits.
  Recommendation: accept.
- Card 0-08b: Infernal's filters miss 6 of 239 tRNA loci on six genomes, 5 in the jumbo phage LAK-A1; turning them
  off costs about 7,900 CPU-hours per 10 Gbp. Recommendation: a dedicated tRNA caller (tRNAscan-SE) for phage and
  archaeal tRNAs, which needs a change to `container/environment.yml` (yours).
- Card C-25: S2C phase 2 needs skani (not in the image) for its ANI-diversity selection, and a choice of annotation
  (full Pfam about 2,500 CPU-hours, *derived*; or representatives only).
- Card 0-10: where users' search results live (recommendation: one archive per run on `/stor`).
- Whether to build PLSDB, and when (inputs fetched; 7,139 of its 72,556 records are already contigs of P1 to P6).
- The data root's README on the pod still says the repository is private.
- Seven genomes of the 50 Gbp selection that NCBI does not serve (recommendation: leave them out).
- Old search runs and packs under old `/stor` roots: keep or delete (nothing is deleted without you).

Earlier items:
- Done 2026-11-09: the PI approved `prodigal=2.6.3` in `container/environment.yml` (test oracle only). The Prodigal test first runs in the image built from `main`.
- Done 2026-11-10: the `read:packages` credential; the packages are now public.
- Confirm the day window for the pod's thread caps (profile assumes 08:00 to 20:00 on weekdays).
- First-10 questions (card 0-03), before the S1 gate.
- The generated catalog file is not byte-identical between two builds of the same manifest (about 65 bytes of DuckDB internals differ; *measured*, card B-34). AGENTS §4.3 says every output is byte-identical. Proposed wording for the PI to accept: "the catalog is generated from the manifest and is excluded from byte identity; its logical content is tested".
- Done 2026-11-11: card 0-08 (ncRNA) is decided by DEC-28.
- Card 0-07: a cluster's id is a hash of its member set, so it changes whenever the cluster gains a protein: 39.9% of NR50 clusters lost their id from `rL.2` to `rL.3` while 94.4% of proteins stayed with their old cluster's main fragment (*measured*). Recommendation: keep the id and add a lineage table that maps each cluster to its best-overlap predecessor.
- Card 0-05: two S5 storage budgets are exceeded (annotations 107 GB against 50; tapes plus postings 75 GB against 25 to 60), because every table now carries a 16-byte `feature_uid`. Recommendation: raise the budgets.
- Card 0-04: keep `token_postings` self-contained (about 36 GB at GTDB scale) or drop the columns the tape already holds (about a third smaller, at the price of a join). Recommendation: keep.
- `l0_manifest` is still one file rewritten per stage, against the append-only rule (found in B-15; no card yet: see B-22 in `TASKS.md`).

## On the pod (wilkcomp02)

Current (2026-11-12): the data root `/tmp/if2453/r10/root` and the 50 Gbp build are described in the handoff at
the top. Release `rN.3`'s search pack: `/tmp/if2453/r50/packs`; `rN.1`'s: `/tmp/if2453/w16/packs` (a copy with its
verification record in `/tmp/if2453/w19/packs`). Durable copies of every fetched input, with hashes in
`sources.yaml`: `/stor/scratch/Wilke/if2453/all-the-data/provisional/` (`r50_inputs/`, `s2_inputs/` with H1, HPGC and
PLSDB, `s1_inputs/`). Pfam seed: `/stor/scratch/Wilke/if2453/all-the-data/l0/refdb/pfam/38.2/Pfam-A.seed.gz`. Test
mirrors: `/tmp/if2453/ci/<ID>`; measurement scratch: `/tmp/if2453/{m08a,m08b,c21f,...}`. The lines below are the
state on 2026-11-10 and older.

- `~/bio-box`: `main`, unchanged. `~/bio-box-next`: `plan/task-cards` at `f9af5b6`, used by the import ladder.
- Provisional phage database with Pfam domains: `/stor/scratch/Wilke/if2453/all-the-data/provisional/rS1.pilot1` (old-code format; verified).
- Old-code pilots kept on local disk: `/tmp/if2453/pilot_s1` (phages), `/tmp/if2453/pilot_cami` (1.48 M short contigs).
- Import ladder with the new code: work in `/tmp/if2453/ladder`, results in `/stor/scratch/Wilke/if2453/all-the-data/provisional/rL`. Stage A (phages, `rL.1`), stage B (1 Gbp of complete prokaryotes, `rL.2`, finished 2026-11-09 23:12 CST) are done and clean. `rL.2p` (meant to be B plus Pfam) is **contaminated**: it was built while stage C's files were being written and took them in (its 12,048,866 features include stage C's genes, with clusters computed without them); use it only to try queries, not for counts. Stage C (10 Gbp, `rL.3`) failed at 02:53 CST on 2026-11-10 in the postings step (card B-29), and was finished at 04:10 CST with the fixed builder's postings: **`rL.3` is complete and verified on `/stor`** (32,815 contigs, 12,048,866 slots, 10,620,684 proteins, 25,822,709 postings, 13.5 GB; `planning/measurements/2026-11-10_import_ladder/rL3.md`). It has no Pfam for its newest stage (about 338 CPU-hours, not started), and it was built by the pre-B-17 code, so it is in the old layout. On `/stor`, `rL.1` and `rL.2` each fail verification on one file, `l0/l0_manifest.parquet`, which the old ingest rewrote per stage (fixed in the code by B-22). Results, scripts and findings: `planning/measurements/2026-11-10_import_ladder/`. New cards from it: B-25 (a release names its stages), B-26 (postings bytes); B-24 and C-04 amended.
- Workflow walk-through outputs: `/tmp/if2453/ssrnap`.
