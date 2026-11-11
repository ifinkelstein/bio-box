# B-37 design: clustering as several processes, each under the budget

Card B-37 (design half, tier H). DEC-24 clusters every release from scratch; DEC-25 exempts the one-process
clustering from the 300 CPU-minute rule until it is split. This note says how to split it and whether the
split gives the same clusters. Measurements: wilkcomp02, 2026-11-10/11 (CST), image
`bio-box-test_74b3e05.sif`, **MMseqs2 `18.8cc5c`** (`mmseqs version`; the same string in `bio-box_ed2dceb8.sif`),
2 threads per process, in `systemd-run --user --scope -p MemoryMax=16G -p MemorySwapMax=0 nice -n 19`. CPU is
user + system seconds of each child from `os.wait4`. Script: `scripts/stress/b37_split_cluster.py`; commands
in the appendix. Inputs: the micro-S0 fixture (854 proteins); the first 100,000 proteins of stage S1 by
`protein_uid` (read only from `/tmp/if2453/r10/root`); a constructed set of 29,998 proteins with dense
families and many length ties (`synth`, seed 37).

**Answer to the key question: yes.** With every per-query module split by query key range and the prefilter
parts run with `--add-self-matches 1`, every MMseqs2 database of the pipeline (43 per run) is equal entry for
entry to the workflow's, and the final (member, representative) pairs and member alignments are equal row for
row, at 1, 2 and 4 splits, on all three inputs. Without `--add-self-matches 1` the intermediate databases
differ (below) although the final clusters happened to be equal. No new PI decision is needed.

## 1. What the `cluster` call expands to (from its `-v 3` log)

`biobox.l2.cluster.run` calls `createdb`, then `cluster` at NR90 (`--min-seq-id 0.9 -c 0.8 --cov-mode 1
--cluster-mode 2`), `createsubdb` of the representatives, `cluster` at NR50 over them, `mergeclusters`, and per
level `createtsv`, `align` of members to representatives (`-a -e inf --min-seq-id 0 -c 0 --alignment-mode 3`)
and `convertalis`. The log shows each `cluster` call is the cascaded workflow, without reassignment:

| | NR90 (`Set cluster sensitivity to -s 1`, `iterations to 1`) | NR50 (`-s 3`, `iterations to 3`) |
|---|---|---|
| linclust (redundancy) | `kmermatcher` (`--kmer-per-seq 21 --alph-size aa:13`) → `rescorediagonal` (`--rescore-mode 0`) → `clust` → `order_redundancy` (`awk` over `pre_clust.index`, in `linclust.sh`) → `createsubdb` ×2 → `filterdb` → `rescorediagonal` (`--rescore-mode 1 --filter-hits 1`) → `align` (`-e 0.001`) → `clust` → `mergeclusters` | the same, with `--comp-bias-corr 1` |
| cascade | `createsubdb`; step 0: `prefilter -s 1 --max-seqs 20 --split 0 --split-mode 2` → `align` → `clust`; `mergeclusters` | steps 0, 1, 2: `prefilter -s 1`, `-s 2`, `-s 3` (the last two `--diag-score 1 --min-ungapped-score 15 --comp-bias-corr 1`), each → `align` → `clust`, `createsubdb` between steps; `mergeclusters` of redundancy and three steps |

`clust` is `--cluster-mode 2 --similarity-type 2 --cluster-weight-threshold 0.9` throughout. Replaying these
38 module calls (15 at NR90, 23 at NR50) one process each, with the logged arguments (`split --splits 1`), gives every database equal to the
workflow's and 489.8 CPU-s against the workflow's 485.6 on the S1 sample (*measured*).

## 2. Where the CPU goes (*measured*, S1 sample, 100,000 proteins, `s1_100k_s1`)

| Module calls | CPU-s | Share |
|---|---|---|
| `align`, NR90 cascade step 0 (one call) | 139.5 | 28.5% |
| `align`, NR50 (linclust 14.6, steps 51.4, 42.6, 51.0) | 159.7 | 32.6% |
| `prefilter`, NR50 steps (19.9, 28.6, 58.3) | 106.8 | 21.8% |
| `prefilter`, NR90 step 0 | 26.1 | 5.3% |
| member `align` (builder), NR50 + NR90 | 37.1 | 7.6% |
| `align`, NR90 linclust | 9.5 | 1.9% |
| `kmermatcher` NR90 + NR50 | 4.7 | 1.0% |
| `clust` (8 calls) | 1.0 | 0.2% |
| everything else (`rescorediagonal` 1.2, `mergeclusters`, `createsubdb`, `filterdb`, `createdb`, `createtsv`, `convertalis`) | 5.4 | 1.1% |
| total | 489.8 | |

Per-query modules (`prefilter`, `align`, `rescorediagonal`) are 98.0% (*measured*). The one-process modules'
share was the same at 25,000, 50,000 and 100,000 proteins (`kmermatcher` NR90 0.6% each time); the total grew
as n^1.03 there (*derived*: ln(489.8/118.0)/ln 4). Where the time goes at 10.6 M proteins is **not measured**
(the release builds log at `-v 1`).

## 3. How each step splits, and what stays one process

- **`prefilter`**: N processes, each over a consecutive key range of the query database (`createsubdb
  --subdb-mode 1`), target the whole database; parts joined in key order (`concatdbs --preserve-keys 1`, `mvdb`).
  A part's query db is not the target db by path, so MMseqs2 does not treat it as a self-search: the self hit
  scores by k-mers (100 instead of 255 on the fixture) and short or masked proteins lose it (fixture keys 413
  and 550; 8 to 26 `align` entries per step on the S1 sample). `--add-self-matches 1` on the parts restores the
  workflow's entries exactly (*measured*, all inputs). MMseqs2's own query split in one process
  (`--split 3 --split-mode 1`) is also exact; its **target split (`--split 3 --split-mode 0`) is not**: 3,081 of
  78,316 entries (NR90 step 0) and 3,452 of 58,787 (NR50 step 2) have another target set (*measured*, `pfcmp`).
  So the target must never be split, and every part's log must say `Process prefiltering step 1 of 1` (with
  `--split 0`, MMseqs2 chooses a split from the memory it sees; the implementing card checks the line).
- **`rescorediagonal`, `align`**: N processes over consecutive key ranges of the input result db
  (`createsubdb --subdb-mode 0`), query and target dbs whole (so the same-db path is kept). Exact.
- **One process**: `kmermatcher` (k-mers of the whole db), `clust` (the whole graph), `mergeclusters`,
  `createsubdb`, `filterdb`, and the cheap `createdb`, `createtsv`, `convertalis`.

Cost of splitting: 4 splits used 646.7 CPU-s against 489.8 (*measured*), almost all in `prefilter`: 282.2 CPU-s
in parts against 132.9 unsplit, because each part builds the whole target's k-mer index: 12.4 CPU-s per extra
part at 59,000 to 78,000 targets (*derived*: (282.2 − 132.9) / 12). Splitting `align` costs little (NR90
cascade step: 137.3 in parts against 139.5; all NR50 `align`: 167.8 against 159.7). A shared precomputed index (`createindex`) gave identical output but `prefilter` rebuilt its table
anyway, so the saving is not shown (B-41 found why, and card B-43 uses it where it works: section 7).

Budget (*derived*). Total for one process, the B-36 fit `284.10 × (n / 2,053,075)^1.244`: 2,195.9 CPU-min at
10.62 M (measured), about 33,600 at 95 M (extrapolated). The one-process modules, at their measured share
(largest: `kmermatcher` NR90, 0.6%): 13 CPU-min at 10.6 M and about 200 at 95 M; at the measured n^1.03 from 100,000
instead: about 46 at 95 M (2.9 s × 950). Both under 300; at 95 M the share bound uses about 83% of the 240 CPU-min
margin of C-11. The largest per-query step at its 100,000-protein share (28.5%): 626 CPU-min at 10.6 M (3 parts
at 240) and about 9,560 at 95 M (40 parts). The index rebuilt per `prefilter` part, linear in target size: about
22 CPU-min per part at 10.6 M and about 196 at 95 M (12.4 s × 106 and × 950). **At 10.6 M every process fits
easily. At 95 M the `prefilter` parts do not fit unless the index is built once and shared**; that is open
(section 6), as is MMseqs2's peak memory per part (each part holds the whole target index; MMseqs2 estimated
1 GB at 100,000, not measured above).

## 4. Equality (*measured*, `compare`)

| Input | Clusters NR90 / NR50 (with ≥ 2 members) | 1 split | 2 splits | 4 splits | 4 splits, no `--add-self-matches 1` |
|---|---|---|---|---|---|
| fixture, 854 | 853 (1) / 831 (6) | all dbs equal | all dbs equal | all dbs equal | 7 dbs differ (every `pref_step*`, 3 `aln_step*`); pairs equal |
| S1 sample, 100,000 | 78,037 (11,021) / 58,529 (14,759) | all equal | all equal | all equal | 8 dbs differ (`pref` 58,762 to 72,268 entries, `aln` 8 to 26); pairs equal |
| synth, 29,998 | 10,038 (7,619) / 2,353 (2,350) | all equal | all equal | all equal | 4 `pref_step*` dbs differ; pairs equal |

"All equal": every database (key and bytes per entry; the data files' layout differs) and, per level, the
member–representative pairs (0 moved) and the sorted member-alignment rows. The workflow at 1 and 2 threads
gives equal pairs and alignment rows on the S1 sample and synth. The TSVs' line order differs between runs even
when every db is equal; the builder sorts them, so its tables do not depend on it. Because the per-query modules
are exact under any key partition and `clust` sees the same alignment db, the result cannot depend on N; the
naive split shows the risk is the self hit, not `clust` order.

## 5. Determinism, records, publishing, lost scratch

- **Pure plan.** The FASTA is the release's proteins sorted by `protein_uid` (as today); `createdb --shuffle 0`
  keys them by rank. Each per-query step's parts are consecutive ranges of its input's sorted keys, cut by
  cumulative residues (sequence step) or entry bytes (result step) at a size from a measured rate and 80% of the
  budget; never from J, C or the hour. Output is equal at any N, so N is scheduling only.
- **Steps.** `cluster.prep` (FASTA, `createdb`), then one `run-step` per module call in the fixed list of
  section 1, per-query calls as `cluster.<level>.<kk>.<module>.part-<p>` plus a join, then the final step
  under the name `cluster` (TSVs, member align in parts before it, the Python side as today). The list is
  hard-coded and a test checks it against the `-v 3` log of the image's MMseqs2, so an upgrade fails loudly.
- **Records.** Each step writes a small record to a new `small_ok`, `in_manifest=False` artifact
  `cluster_plan` (`l2/cluster_plan/release=<R>/<step>.json`, as `domain_scan`): its input digests, its key range
  and the digest of its output db, sha256 over (key, entry) in key order, not over files (the `.0/.1` layout
  depends on threads). Parts run up to J at once within the allowance (C-11, C-12 machinery).
- **Published.** Only what `cluster` publishes today (`cluster_member`, `cluster_member_aln`, `member_cluster`,
  the ledger, `release_def`) plus the records. MMseqs2 dbs stay in `<scratch>/biobox-cluster/<plan digest>/` and
  are removed after `cluster` publishes.
- **Lost scratch.** The chain is sequential, so lost scratch means rerunning from the earliest step whose output
  is missing or fails its digest; each step with a record must reproduce it, and a mismatch stops the build
  naming the step (as C-12 and B-30). Worst case one more full pass.

## 6. Proposed implementing card (the M half of B-37)

**M** · Depends: B-36, C-13 · Files: `src/biobox/l2/cluster.py`, `src/biobox/stage.py`, `src/biobox/cli.py`
(`release cluster plan|step|finish|redo`), `src/biobox/registry/schemas.py` (one `ArtifactSpec`,
`cluster_plan`, nothing else), `RUNBOOK.md`, `tests/r2/test_micro_s0_cluster.py`, `tests/r2/test_build_stage.py`.
Do: section 5; parts run `prefilter` with `--add-self-matches 1`, refuse a part log without `Process
prefiltering step 1 of 1`, and never split the target. Before the step list is fixed, find out why `prefilter`
did not use a `createindex` index; if a shared index gives identical output, build it once per cascade step.
Tests first: the module list equals the workflow's `-v 3` log; split builds at 1, 2 and 4 parts equal
`cluster.run` file for file on the fixture and on `synth`; a part without `--add-self-matches 1` makes the
intermediate dbs differ (the reason the flag is there); a deleted scratch db is rebuilt to its digest and the
build ends byte-identical; a part's estimate over 80% of the budget is cut smaller.
Accept: as card B-37 (equal row for row on the fixture and on the 2.05 M stage; every record under 300
CPU-minutes); then remove the DEC-25 exemption from `stage.py` and the runbook. Stop if any equality fails, or if
at 2.05 M a one-process module passes 240 CPU-minutes.

## 7. A shared prefilter index (card B-43)

Measured on wilkcomp02, 2026-11-11, MMseqs2 `18.8cc5c` in `bio-box-test_74b3e05.sif`, at most 2 threads, in
`systemd-run --user --scope -p MemoryMax=16G -p MemorySwapMax=0 nice -n 19`; scripts and outputs in
`/tmp/if2453/b43` (`exp1.sh` to `exp6.sh`, `e1/<input>/table.tsv`, `e*.out`). `prefilter` uses a prebuilt index
only when the index itself (`<db>.idx`) is its target (B-41).

**Why `-s 1 --comp-bias-corr 0 --diag-score 0` aborts.** Every attempt below aborted with `munmap_chunk(): invalid
pointer` (exit -6) after its scoring pass, on the fixture (854 proteins) and on the constructed set (29,998): an
index built with the call's own `-s 1 --comp-bias-corr 0 -k 0 --mask 1 --spaced-kmer-mode 1 --alph-size aa:21,nucl:5
--split 0`; the same with the whole query, half the keys as query, `--add-self-matches 0` or `1`; prefilter at 2
threads; `--db-load-mode 1`, `2`, `3`; the index built at 2 threads; with `--split 1`; `--index-subset 1` and `4`
(`2` and `3` are refused: "Index was not built with `prefilter` support"); `--check-compatible 1` and `2` over an
old index of other settings and over a fresh one; `-k 6` and `-k 7` in both commands. One setting changed at a time,
with the index rebuilt to match:

| `-s` | `--comp-bias-corr` | `--diag-score` | `--min-ungapped-score` | against the index |
|---|---|---|---|---|
| 1, 2, 3 | 0 or 1 | 0 | 0 or 15 | aborts (6 combinations) |
| 1, 2, 3 | 0 or 1 | 1 | 0 or 15 | equal to the rebuilt call, whole query and half (6 combinations, both inputs) |

So the abort follows `--diag-score 0` alone, whatever the index was built with: it is a defect of this MMseqs2 in
prefilter's path without diagonal scoring when the target is an index, and no index setting avoids it. Smallest
reproducer (`exp4.sh`, fixture, one thread): `mmseqs createdb in.faa db; mmseqs createindex db tmp; mmseqs prefilter
db db.idx pf --diag-score 0` exits 134 (`Aborted`); the same with `--diag-score 1` exits 0. Changing the calls'
`--diag-score` is not an option: it changes the workflow's results. Whether a newer MMseqs2 differs is **unknown**
(no changelog offline; nothing was downloaded).

**Where it works it is exact.** With `--diag-score 1` (the NR50 cascade's `-s 2` and `-s 3` steps) a part against
the index, asking for the self hit, equals the rebuilt call entry for entry: on the fixture and the constructed set
(whole query and half); on the 100,000-protein S1 sample's step-1 and step-2 targets (59,338 and 58,787) whole and in
4 parts; for a 1,000-query part at 250,000, 500,000 and 1,000,000 S1 targets. Every part log says `Process
prefiltering step 1 of 1`. Tests (`tests/r2/test_cluster_split.py`, `test_b43_*`): at 2 and 4 parts the joined result
of both calls equals the workflow's call and B-41's rebuilt parts; the whole query against the index (1 part) equals
it too; the parts' logs show the index read and no table built; a lost index is made again and the build ends
byte-identical.

**The index's bytes are not reproducible.** The same index built twice from the split path's target at one thread
gave other bytes in its DBR index entries (2 of 22 entries, 40 of 12,077 bytes, on the constructed set; on the
fixture another entry at 2 threads); built from `createdb`'s own database it was equal at 1 and 2 threads and run to
run (`exp3.sh`, `exp5.sh`). The parts' outputs were the same every time. So the cut's record holds what the index is
made of (path, `createindex` argv, target digest, and their sha256 `id`), not a digest of its bytes; a part's record
names that `id`; a part made again over a new index must reproduce its own output digest, which is the check.

**Cost** (one thread unless said). `createindex`: 39.58, 78.05, 136.40 CPU-s at 250,000, 500,000, 1,000,000 targets
(15.0 to 16.8 at 58,787 to 78,316, fixed cost dominating); the code takes 1.583e-4 CPU-s per target (*derived*,
the largest rate at 250,000 or more). A part's load of the index: 0.28 CPU-s per extra part at 59,338 targets (2
threads) and at most 1.29 CPU-s at 1,000,000; the code takes 4.72e-6 per target (*derived*). A part rebuilding the
table: 46.0, 72.9, 147.7 CPU-s at the three sizes (1.5e-4 to 1.8e-4 per target), so B-41's 2.108e-4 is on the safe
side and is kept. Index size: 1.51, 1.83, 2.45 GB at the three sizes, about 1.26 KB per target over 1.2 GB
(*derived*); a part's peak RSS is the same as a rebuilding part's (2.12 against 2.11 GB at 1,000,000 targets).

Budget (*derived*, the code's model; the formulas are in `biobox/l2/cluster.py` beside
`PREFILTER_INDEX_CPU_S_PER_TARGET`):

| n proteins | one process | in parts, rebuilt | in parts, shared index | NR50 steps 1, 2 parts (rebuilt / shared) | largest single process |
|---|---|---|---|---|---|
| 10 M | 2,036 | 2,148 | 2,103 | 1, 2 / 1, 2 | kmermatcher 12.2; index 15.5 |
| 50 M | 15,079 | 18,910 | 16,918 | 7, 14 / 4, 8 | kmermatcher 90.5; index 78.3 |
| 95 M | 33,506 | refused | refused | 47, 92 / 9, 17 | NR90-step part: 261.5 of overhead alone |

At 95 M the shared index turns the NR50 steps' rebuilt tables (27,371 CPU-min) into 410, but the NR90 step (78% of
the proteins as targets) cannot be cut: each part's own table is estimated at 261.5 CPU-min, over 240. The cut
refuses the NR90 step from 87.2 M proteins and NR50 step 0 from 108.1 M; the first process over 300 CPU-minutes would
be an NR90-step part at 109.0 M; `kmermatcher` (0.6% of the one-process total) passes 300 at 131.1 M and
`createindex` at 191.6 M (over 240 from 153.3 M). With the measured rebuild rate at 1,000,000 targets (1.46e-4)
instead of 2.108e-4, an NR90-step part at 95 M would carry about 181 CPU-min of table and fit, at about 31 parts and
5,600 CPU-min of tables (*estimate*, extrapolated 74 times).

### Proposal: an MMseqs2 for the `--diag-score 0` calls (for the PI; `container/environment.yml` is not changed)

A version change is the only way found to give the NR90 step and NR50 step 0 a shared index. Proposal: pin a newer
MMseqs2 release in `container/environment.yml` only after it passes, in the image built from it: (1) `exp4.sh`'s
reproducer exits 0 with `--diag-score 0`; (2) the B-41 module list test (`-v 3` log) passes, or the list is updated in
the same card; (3) `tests/r2/test_cluster_split.py` at 1, 2 and 4 parts on the fixture and the constructed set, with
those two calls added to `INDEX_SETTINGS` (their parts against the index equal the rebuilt calls entry for entry);
(4) `cluster.run` with the new version equals today's clusters row for row on the fixture and on `rN.1b` (2,053,075
proteins): DEC-24 clusters every release from scratch, so a version whose clusters differ would move members between
releases without any input changing. If (4) fails, the change needs a decision row, not a pin. Which version fixes
the abort is unknown here; nothing was downloaded.

## Appendix: commands (on wilkcomp02, in `/tmp/if2453/b37`)

`x.sh U SRC ARGS` runs `python tree/scripts/stress/b37_split_cluster.py ARGS` in the image in the capped scope
`b37-U`, after checking the account's `ps -u if2453 -o pcpu=` total plus 200% against 5,600% (20:00 to 05:00)
or 4,000%.

```
python tree/scripts/stress/build_fixture.py tree /tmp/if2453/b37/fixture rF        # tree = this branch
x.sh export-fx tree/src export fixture S0 in/fixture.faa
x.sh export-s1 /tmp/if2453/ci/r10/src export /tmp/if2453/r10/root S1 in/s1_100k.faa --limit 100000
x.sh synth tree/src synth in/synth.faa
x.sh wf-I  tree/src workflow in/I.faa runs/I_wf --threads 2          # I = fixture, s1_100k, synth; also --threads 1
x.sh s1-I  tree/src split in/I.faa runs/I_s1 --splits 1 --threads 2 --plan-from runs/I_wf
x.sh x4-I  tree/src split in/I.faa runs/I_x4 --splits 4 --threads 2 --plan-from runs/I_wf --self-matches   # and 2
x.sh s4-I  tree/src split in/I.faa runs/I_s4 --splits 4 --threads 2 --plan-from runs/I_wf                  # naive
x.sh cmp-I tree/src compare runs/I_wf runs/I_x4 > runs/cmp_I_x4.json
head -n 50000 in/s1_100k.faa > in/s1_25k.faa         # and 100000 lines for 50k; then split --splits 1
sh pf.sh      # in the image: prefilter on runs/s1_100k_wf's NR90 step-0 and NR50 step-2 inputs, as logged,
              # with --split 0, --split 3 --split-mode 0, --split 3 --split-mode 1, and after createindex
python3 pfcmp.py runs/pf A90 T90 Q90 X90              # and A50 T50 Q50 X50
```
