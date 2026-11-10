# Task status

The running record of the cards in [`TASKS.md`](TASKS.md). Updated at every merge. Branch
`plan/task-cards` holds all merged work; nothing is on `main` yet. Each card was done by an agent on
its own branch `task/<ID>`, audited by the orchestrating session, and merged after the full suite
passed on the integration branch.

## Current state, 2026-11-10 evening (after 17:10 CST; supersedes the sections below where they differ)

- Merged since the section below: B-27, E-01, C-12. 66 cards. The pod suite with C-12 merged:
  1,148 passed, 12 skipped. `ci` and `cross-host` were green on the last three pushes
  (`740d049`, `dd6103a`, `94b6f8d`); the fixture's manifest hash is equal on Actions and the pod.
- All agents were cut off by the usage limit in the afternoon and resumed after 17:10 CST.
- Running: the 10 Gbp rebuild (`build/r10`; phages and the 250-genome stage built, the 2,247-genome
  stage in ingest when the agent was cut off), W-03, C-13, B-30 (card amended: one `ingest_plan`
  artifact and `stage ingest` commands). C-13 and B-30 are built on C-12.
- Held: W-13, until W-03 merges (same search files). E-07 waits for PI card 0-08.
- Then: the merge to `main`; then pull the image by digest; then the 50 Gbp build.

## State at about 11:30 CST

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

## Handoff, 2026-11-10 (written about 10:05 CST)

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

## Merged (67 cards)

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

Decisions DEC-18 to DEC-23 were approved by the PI on 2026-11-09 and are in `planning/DECISIONS.md`.

## CI is not reliably green (found 2026-11-10)

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

## In progress

| Card | Model | State |
|---|---|---|
| C-13 defects from the first real build | Opus | started 2026-11-10 |
| B-30 ingest by batch processes | Opus | started 2026-11-10 |
| 10 Gbp rebuild with `build-stage` (not a card; branch `build/r10`, pod root `/tmp/if2453/r10`) | Opus | started 2026-11-10 |

## Next, in order

1. The path to a one-command stage build, in order: B-23 (in rework) → B-22 (`l0_manifest` per stage; register `store_params.json`) → C-08 (every builder writes to `out=`; after B-29) → C-04 (`build-stage`) → C-06 (runbook). Until C-04 merges, stages are built by a hand-written driver, as the import ladder was.
2. B-28 (query-time readers take the release's stages), B-24 (streaming clustering and sharded scans; memory, not time, is the limit), B-30 (ingest by batch), B-27 and W-11.
3. Pull request `plan/task-cards` → `main`, merged when the suite and CI are green (the PI delegated the timing). The merge builds the image and runs the image-only tests for the first time.
4. C-04 `build-stage`; C-06 runbook.
5. W-03 pack and runs outside the data root; W-05 readable exports; W-07 the one command.
6. B-21 per-domain postings.
7. Phase E (close S0), phase F (statistics; F-01 is a design card for a stronger model or the PI), phase G (S1 as a release).

## Not yet executed anywhere

These exist and pass their laptop-runnable parts, but their real run needs the image:
the clustering tests; the samtools slice and index tests; the HMMER-binary oracles (two tests);
the golden pairs test over postings; the Prodigal test (also needs `prodigal=2.6.3` added to
`container/environment.yml`, a PI decision); the mutation workflow; the systemd path of `run-step`
and of the search worker (only their command lines are unit-tested; the scope behaviour itself was
probed by hand on the pod).

## Open items for the PI

- Done 2026-11-09: the PI approved `prodigal=2.6.3` in `container/environment.yml` (test oracle only). The Prodigal test first runs in the image built from `main`.
- A `read:packages` credential, so the image built from `main` can be pulled on the pod.
- Confirm the day window for the pod's thread caps (profile assumes 08:00 to 20:00 on weekdays).
- First-10 questions (card 0-03), before the S1 gate.
- The generated catalog file is not byte-identical between two builds of the same manifest (about 65 bytes of DuckDB internals differ; *measured*, card B-34). AGENTS §4.3 says every output is byte-identical. Proposed wording for the PI to accept: "the catalog is generated from the manifest and is excluded from byte identity; its logical content is tested".
- Card 0-08: ncRNA. Keep every overlapping Rfam hit or apply clan competition (recommended: clan competition, one tape slot per gene); scan full Rfam everywhere (about 2,100 CPU-hours per 10 Gbp, five times Pfam) or a curated prokaryote and phage subset (recommended).
- Card 0-07: a cluster's id is a hash of its member set, so it changes whenever the cluster gains a protein: 39.9% of NR50 clusters lost their id from `rL.2` to `rL.3` while 94.4% of proteins stayed with their old cluster's main fragment (*measured*). Recommendation: keep the id and add a lineage table that maps each cluster to its best-overlap predecessor.
- Card 0-05: two S5 storage budgets are exceeded (annotations 107 GB against 50; tapes plus postings 75 GB against 25 to 60), because every table now carries a 16-byte `feature_uid`. Recommendation: raise the budgets.
- Card 0-04: keep `token_postings` self-contained (about 36 GB at GTDB scale) or drop the columns the tape already holds (about a third smaller, at the price of a join). Recommendation: keep.
- `l0_manifest` is still one file rewritten per stage, against the append-only rule (found in B-15; no card yet: see B-22 in `TASKS.md`).

## On the pod (wilkcomp02)

- `~/bio-box`: `main`, unchanged. `~/bio-box-next`: `plan/task-cards` at `f9af5b6`, used by the import ladder.
- Provisional phage database with Pfam domains: `/stor/scratch/Wilke/if2453/all-the-data/provisional/rS1.pilot1` (old-code format; verified).
- Old-code pilots kept on local disk: `/tmp/if2453/pilot_s1` (phages), `/tmp/if2453/pilot_cami` (1.48 M short contigs).
- Import ladder with the new code: work in `/tmp/if2453/ladder`, results in `/stor/scratch/Wilke/if2453/all-the-data/provisional/rL`. Stage A (phages, `rL.1`), stage B (1 Gbp of complete prokaryotes, `rL.2`, finished 2026-11-09 23:12 CST) are done and clean. `rL.2p` (meant to be B plus Pfam) is **contaminated**: it was built while stage C's files were being written and took them in (its 12,048,866 features include stage C's genes, with clusters computed without them); use it only to try queries, not for counts. Stage C (10 Gbp, `rL.3`) failed at 02:53 CST on 2026-11-10 in the postings step (card B-29), and was finished at 04:10 CST with the fixed builder's postings: **`rL.3` is complete and verified on `/stor`** (32,815 contigs, 12,048,866 slots, 10,620,684 proteins, 25,822,709 postings, 13.5 GB; `planning/measurements/2026-11-10_import_ladder/rL3.md`). It has no Pfam for its newest stage (about 338 CPU-hours, not started), and it was built by the pre-B-17 code, so it is in the old layout. On `/stor`, `rL.1` and `rL.2` each fail verification on one file, `l0/l0_manifest.parquet`, which the old ingest rewrote per stage (fixed in the code by B-22). Results, scripts and findings: `planning/measurements/2026-11-10_import_ladder/`. New cards from it: B-25 (a release names its stages), B-26 (postings bytes); B-24 and C-04 amended.
- Workflow walk-through outputs: `/tmp/if2453/ssrnap`.
