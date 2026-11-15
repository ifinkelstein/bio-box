# CRISPR arrays of both producers on the seven built stages of the data root (card C-33 part 1, DEC-40)

Date: 2026-11-15, wilkcomp02, 01:29 to 02:08 CST (the pod's clock). Every number carries a label: *measured* (host,
date, command or record), *derived* (formula and inputs) or *estimate*. Part 1 of the card as amended on 2026-11-15:
the two producers onto stages S1, P1 to P6 of `/tmp/if2453/r10/root`, nothing else. No release was defined and no
release step was run; nothing under any `release=` directory was written; `rM.1` was not touched.

## Result

- **56 files added to the root, 26,958,140 bytes, and nothing else changed** (*measured*, `c33_check.py`, below):
  for each of the 7 stages and 2 producers one `feature` part, one `l1_crispr` part, one ledger part and one run
  record. The list with sizes and sha256 is [`root_added.tsv`](2026-11-15_crispr_arrays_on_built_stages/root_added.tsv)
  (file sha256 `a53473d9c4ed8883b60cdc65b8ecb556233b0c59798fe370590a15a20b020e06`, the same as the pod's
  `/tmp/if2453/c33/root_added.tsv`). By directory: `l1/crispr` 14 files, 24,259,014 bytes; `l1/feature` 14, 515,769;
  ledgers 14, 2,157,457; run records 14, 25,900 (*derived*: sums of the list).
- **No existing file changed**: the 3,347 files of the root before are all there after with the same size and mtime,
  and the 90 files of stages P1 and S1 have the same sha256; the 56 new files are exactly the expected ones.
- **Oracle: byte-identical on S1, P1 and P2.** The new `feature`, `l1_crispr` and ledger files of both producers
  equal those of `/tmp/if2453/e15/view`, and the default producer's those of `/tmp/if2453/e12/view`, by sha256:
  27 comparisons of files, all equal ([`oracle.tsv`](2026-11-15_crispr_arrays_on_built_stages/oracle.tsv)).
- **Both releases still verify** with the manifest hashes they had (below).
- Arrays: 13,477 of the default producer and 414 of the filtered one over 52,132,698,276 bases; 146.9 and 171.8
  CPU-minutes, 7% and 9% under the estimates; no process over 33.4 CPU-minutes.

## Setup

- Code: the tree of branch `task/C-33` (`origin/plan/task-cards` at `0d82b05` plus this card's scripts; no `src/`
  change), synced to `/tmp/if2453/ci/C-33` by `scripts/pod_pytest.sh C-33 tests/r0` before any run. Every run record
  of this card names `code_id` `6d3c4606bd136e8a66d3a5804819dbbb5d0b5494c8228c0a00cbdba70e8070f4` and `code_commit`
  null (the synced tree has no `.git`); the views were made by other trees (E-15's `code_id` `d05b496e...`), with the
  same producer code as the bytes show.
- Image: `/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif`, the runtime image of `faec1f4` (RUNBOOK §0), checked
  at every launch by `pod_run.sh --digest`; every run record of this card names `image_digest`
  `sha256:d67300aa1321b571fa8da79e74dfccdd87a8dca2e56518f324902efa6b9f098a`. The earlier steps of the seven stages name
  the old image `sha256:c30be58557d7...` (the `rM.2` build, STATUS); `release verify` reads no image digest, and both
  releases verify with the mixed records (below).
- Profile: `profiles/wilkcomp02.yaml` of the tree (staging under `/tmp/if2453/bio-box/staging/`, removed by run-step
  at each publish). Scratch (`--tmp`): `/tmp/if2453/c33/scratch/<stage>.<step>`. Logs: `/tmp/if2453/c33/logs/<unit>.log`;
  driver log `/tmp/if2453/c33/driver.log`.
- Driver: [`scripts/stress/build_stage_pod/c33_driver.sh`](../../scripts/stress/build_stage_pod/c33_driver.sh), with
  `c33_view.py` (the rehearsal root) and `c33_check.py` (the checks), run from copies in `/tmp/if2453/c33/bin`, one
  mode per unit: `systemd-run --user --unit=c33-driver-<mode> --collect -p MemoryMax=4G -p MemorySwapMax=0
  -p StandardOutput=append:/tmp/if2453/c33/driver.out -p StandardError=append:/tmp/if2453/c33/driver.out nice -n 19
  bash /tmp/if2453/c33/bin/c33_driver.sh <mode>`, modes `rehearse`, `before`, `run`, `after`, `check` in that order.
- Each step: one unit through `pod_run.sh --profile <profile> --image <sif> --digest <digest> --mem 16G --threads 4
  --bind /tmp/if2453 --env PYTHONPATH=/tmp/if2453/ci/C-33/src --log <log> <unit> run-step --profile <profile> --root
  /tmp/if2453/r10/root --stage S --step crispr-<producer> --threads 4 -- python -m biobox.cli stage crispr
  [--near-cas --library pfam --library-version 38.2] --root /tmp/if2453/r10/root --stage S --cpus 4 --tmp <scratch>
  --seq-stages <S, then the other six>`. Units `c33-nr3-<s>-<n>` (default producer) and `c33-near-<s>-<n>` (filtered).
- The argv is the scratch views' (cards E-12, E-15; their run records), with two differences that do not reach the
  Parquet bytes: `--cpus 4` (E-12 ran 4, E-15 ran 3; the producer's bytes do not depend on it, tested in E-12 and
  shown again by the oracle), and `--seq-stages`, which is what `build-stage` gives (`with_seq_stages` in
  `src/biobox/stage.py`: the stage itself, then every other stage with a `seq_locator` in the root, sorted by name),
  so `S1,P1,P2,P3,P4,P5,P6` for S1 and for example `P3,P1,P2,P4,P5,P6,S1` for P3. Every stage of the root had to be
  given: P3 to P6 each hold one contig an earlier stage stored (below), and a contig owned by no listed stage is an
  error (`ContigOwners`). The views listed only S1, P1 and P2, which own all their contigs; one stage owns each
  contig (DEC-26), so the longer list cannot change their files, as the oracle shows. No `--gene-call-version` was
  given, as in the views: the default `pyrodigal-3.7.1-gv-0.3.2-p2` is the root's one version.
- Threads: at most 5 units at once of 4 threads, 20 threads, within the 24 of the brief; before each launch the
  driver read the account's allowance (`python -m biobox.runstep --threads 4`, card C-28): 41 to 56, other users 0.01
  to 0.11 CPUs, the night cap of 56 (*measured*, driver log). Another agent's work did not show (ours elsewhere 0 to
  14.7 CPUs, all of it this run's own units).

## The rehearsal (*measured*, wilkcomp02, 2026-11-15, 01:29 to 01:34, units `c33-rh-nr3-s1`, `c33-rh-near-s1`)

`c33_view.py /tmp/if2453/r10/root /tmp/if2453/c33/rehearsal S1 S1,P1,P2,P3,P4,P5,P6`: real directories holding 30
symlinks (the `contig_index`, `seq_locator`, `seq_shards`, `genome_input` of all seven stages, S1's gene calls and
S1's Pfam 38.2 scan) and `store_params.json`, with its own `runs/`, `ledger/`, `l1/feature`, `l1/crispr`; a scratch
copy of the profile with `tmp: /tmp/if2453/c33/rehearsal-tmp`. Then the exact commands of the real run for S1 with
`--root /tmp/if2453/c33/rehearsal`: the default producer, then the filtered one. The six new files had the sha256 of
the views' (driver log: `rehearsal: S1 nr3: 3 files byte-identical to /tmp/if2453/e15/view`, the same against
`/tmp/if2453/e12/view`, and `S1 near: 3 files byte-identical to /tmp/if2453/e15/view`;
[`rehearsal-cmp/`](2026-11-15_crispr_arrays_on_built_stages/rehearsal-cmp)). Cost 5.282 and 5.499 CPU-minutes
(run records).

What differs between the rehearsal (and the views) and the real run is only where the records go: the views and
the rehearsal were made with `--root` the scratch root and run-step's staging, not with `--out` (the views' run
records show `--root /tmp/if2453/e15/view` and no `--out`), so their run records and ledgers went under the scratch
root; in the real run they go under the data root at the same relative paths. The run record's `argv`, `started_at`,
`mem_cap` (its unit), times and `code_id` differ; its `outputs` (paths, bytes, sha256) are the same. The Parquet
bytes, ledgers included, are the same.

## Before the first write (*measured*, 01:34:59 to 01:39:11, unit `c33-driver-before`)

- No `STOP` in the root. Active units of the account: `c33-driver-before` and two `w27-fl-nr50*` units in state
  `active exited` (finished `nbr flanks` runs of 2026-11-13, no process); no `r50-*`, `c24-*`, `c25-*`, `c26-*`,
  `b47-*` unit.
- `release verify` (root bound read-only, [`verify-logs/before/`](2026-11-15_crispr_arrays_on_built_stages/verify-logs)):
  - `release rM.2: verified (stages S1,P1,P2,P3,P4,P5,P6, 158 files, 154 artifacts, manifest hash 2239df7b92ec775058fb9c2ddae3609e)`
  - `release rN.3: verified (stages S1,P1,P2, 60 files, 34 artifacts, manifest hash eee598acd16324d1dd73349cfd6d931f)`
- Snapshot `/tmp/if2453/c33/before/`: `listing.tsv` (relative path, bytes, mtime of every file: 3,347 files,
  80,969,212,936 bytes) and `p1s1.sha256` (the 90 files under a `stage=P1` or `stage=S1` directory: tables, run
  records, ledgers).

## The run (*measured*, 01:40:31 to 02:02:38, unit `c33-driver-run`)

The default producer on S1, P1, ..., P6 (01:40:39 to 01:50:34), then the filtered one (01:51:12 to 02:02:08): 14
units, each `exited 0` with `published 3 files`; 22 min 7 s in all (*derived* from the driver log). No STOP, no
refusal, no guard action (the guard stops a unit of this run above 280 CPU-minutes; the largest used 33.4).

Per stage and producer ([`per_stage.tsv`](2026-11-15_crispr_arrays_on_built_stages/per_stage.tsv); arrays are the
`feature` rows, equal to the ledger's kept; calls in and rejections from the ledger, counted once per contig; cost
from the run record; unit peak from `pod_run: <unit> usage`):

| stage | bases | default: arrays | default: CPU-min (estimate) | filtered: calls at 2 | rejected (a) overlap / (b) far / (c) inside | filtered: kept | filtered: CPU-min (estimate) |
|---|---|---|---|---|---|---|---|
| S1 | 1,664,109,938 | 213 | 5.085 (5.048) | 3,613 | 219 / 3,360 / 0 | 34 | 5.487 (6.007) |
| P1 | 1,016,625,803 | 286 | 3.085 (3.084) | 3,379 | 294 / 3,076 / 3 | 6 | 3.383 (3.670) |
| P2 | 9,596,180,328 | 2,291 | 26.489 (29.108) | 30,682 | 2,354 / 28,243 / 13 | 72 | 31.794 (34.639) |
| P3 | 9,949,077,522 | 3,921 | 29.195 (30.179) | 39,391 | 4,039 / 35,212 / 22 | 118 | 33.386 (35.913) |
| P4 | 9,976,456,950 | 2,520 | 27.585 (30.262) | 30,835 | 2,567 / 28,194 / 5 | 69 | 32.606 (36.012) |
| P5 | 9,967,711,632 | 2,205 | 27.505 (30.235) | 29,779 | 2,278 / 27,437 / 8 | 56 | 32.361 (35.980) |
| P6 | 9,962,536,103 | 2,041 | 27.938 (30.220) | 30,953 | 2,090 / 28,794 / 10 | 59 | 32.796 (35.961) |
| all | 52,132,698,276 | 13,477 | 146.882 (158.136) | 168,632 | 13,841 / 154,316 / 61 | 414 | 171.813 (188.182) |

- Bases: `stage_bases` (the default variants' `contig_len` of each stage's `contig_index`, `src/biobox/stage.py`),
  *measured*; their sum is the brief's 52,132,698,276. The default producer's calls in equal its arrays (no array
  longer than the circle on any stage); no call of either producer was rejected as longer than the circle.
- Estimates (*derived*): bases × 182.0 CPU-seconds per Gbp / 60 (`CRISPR_CPU_S_PER_GBP`), times 1.19 for the filtered
  producer (`CRISPR_NEAR_CAS_COST`). Measured over estimate: 0.910 to 1.007 for the default producer, 0.899 to 0.930
  for the filtered one; all within 20% (AGENTS §6). In all 146.882 against 158.136 (0.929) and 171.813 against
  188.182 (0.913), *derived*. The default producer ran at 165.6 to 183.4 CPU-seconds per Gbp, 169.0 over the seven
  stages (*derived*: 146.882 × 60 / 52.133).
- The filtered step over the default on the same stage (*derived*, CPU-minutes): S1 1.079, P1 1.097, P2 1.200, P3
  1.144, P4 1.182, P5 1.177, P6 1.174. P2's 1.200 is 0.9% above the factor 1.19 that E-15 took from its own P2 run
  (31.790 / 26.723); the filtered steps' estimates still hold (the default steps ran 7% to 9% under 182.0), and no
  refusal threshold is near (66.5 Gbp a stage against 10.0 at most).
- Contigs searched (ledger `ok` rows): S1 27,258; P1 539; P2 5,018; P3 9,107; P4 4,005; P5 4,111; P6 4,227. Skipped
  as another stage's (DEC-26, ledger reason `owned_by:<stage>`, both producers alike): one contig each of P3
  (owned by P1), P4 (P2), P5 (S1) and P6 (S1), the four contigs DEC-26 counts as shared across stages. No contig of
  any stage is in several of its genomes (no `owned_by:<stage>/<genome>` reason).
- S1, P1 and P2 equal cards E-12 and E-15 (arrays 213, 286, 2,291; kept 34, 6, 72; the rejections by rule; E-15's
  table): as the byte identity implies. The CPU-minutes are within 1% of E-15's on P1, P2 (3.109, 26.723, 31.790)
  and 10% and 15% under on S1 (5.670, 6.479), *derived*.
- Kept by the filtered producer per Gbp (*derived*, kept / bases): P3 11.86, P4 6.92, P5 5.62, P6 5.92; P1 to P6
  together 7.53 (380 / 50.469 Gbp), against 7.35 for P1 and P2 in card E-13 and 20.43 for S1 (card E-16 is the open
  measurement of phages). P3's rate is the highest of the prokaryote stages; nothing here says whether it is chance.
- Wall seconds 50.1 to 516.0; peak RSS of one process (`peak_rss_bytes`, `child_maxrss`) 355.1 to 561.3 MB for the
  default producer and 454.6 to 1,199.3 MB for the filtered one (P3; the filter's DuckDB connection, capped at 2 GB);
  unit MemoryPeak 1.46 to 1.96 GB (systemd, page cache included), under the units' 16 GiB. Threads 4 in every record;
  no `over_cpu_budget`.

## After the last step (*measured*, 02:02:53 to 02:07:16, units `c33-driver-after`, `c33-driver-check`)

- Snapshot `/tmp/if2453/c33/after/`: 3,403 files, 80,996,171,076 bytes; 106 files under `stage=P1` or `stage=S1`.
- `c33_check.py` (in the image, root read-only; [`report.txt`](2026-11-15_crispr_arrays_on_built_stages/report.txt)):
  `before: 3347 files; after: 3403 files; new 56; gone 0; size or mtime changed 0; P1 and S1 files hashed before 90,
  after 106, sha256 changed 0`; the 56 new files are the 14 run records and the 42 Parquet parts of the producers'
  directories, none else (16 of them among P1 and S1); `RESULT: all checks pass`. 56 = 3,403 - 3,347 (*derived*).
- Top-level entries of the root newer than the snapshot: the root itself and `l1` (the new `l1/crispr` directory and
  `.incoming`, which run-step creates and removes at each publish); no `.incoming` or intent file is left.
- `release verify` again ([`verify-logs/after/`](2026-11-15_crispr_arrays_on_built_stages/verify-logs)):
  - `release rM.2: verified (stages S1,P1,P2,P3,P4,P5,P6, 158 files, 154 artifacts, manifest hash 2239df7b92ec775058fb9c2ddae3609e)`
  - `release rN.3: verified (stages S1,P1,P2, 60 files, 34 artifacts, manifest hash eee598acd16324d1dd73349cfd6d931f)`

  The same lines as before. The second line, `built by:`, is not the same: it reads the run records of the release's
  stages, and now counts the new records too (`6d3c4606bd13 (version 0.0.1, commit unknown, 14 steps)` for `rM.2`,
  `6 steps` for `rN.3`, beside `ec1272b840c6 ... 449 steps` and the steps from before card C-23). The releases'
  manifests do not list these steps; a reader of that line could take it that they built the release.

## For the durable copy and part 2

- The durable copy needs a tar of exactly the 56 files of `root_added.tsv` (the orchestrator's; the four tars of
  2026-11-14 stay valid, since no file of the root changed).
- Which release holds the arrays is open for the PI (STATUS open item 6); `rM.2` and `rN.3` do not list them, and
  `nbr` on those releases does not show them.
- Not changed here: the guides' worked command that names `/tmp/if2453/e14/view`, and the example's `recorded/`
  files (part 2).
