# Clustering memory before and after card B-24 (wilkcomp02, 2026-11-10)

Card B-24, part 1. Where does `biobox.l2.cluster.run` hold its memory, and what is left after the
streaming rewrite? Every number is *measured* on wilkcomp02 on 2026-11-10 (CST) unless it is
labelled *derived* (formula given).

## Method

- **Input.** The import ladder's release `rL.2` inputs: stages S1 and P1 of
  `/tmp/if2453/ladder/root` (read only; its `protein_seq` layout is the current one, so no
  synthetic stand-in was needed), 2,053,075 distinct proteins. The release was named `rB24` and
  written to `/tmp/if2453/b24/<run>/out` (`out=`), so nothing in the ladder root changed.
- **Code.** "Before" is `plan/task-cards` at `121503d`; "after" is `task/B-24` at the part-1
  commit. Both ran in `bio-box-test_74b3e05.sif` with MMseqs2 at 16 threads (`threads=16`), as the
  ladder did.
- **Command.** `systemd-run --user --unit=b24-<run> --collect -p MemoryMax=64G -p MemorySwapMax=0
  nice -n 19 sh /tmp/if2453/b24/launch_prof.sh <run> <tree> 16`, which runs
  `singularity exec ... python /tmp/if2453/b24/profile_cluster.py /tmp/if2453/ladder/root rB24
  S1,P1 <out> <tmp> 16 <result.json>`. The profiler (scratch, on the pod only) calls
  `cluster.run` in its own process and records:
  - **Python:** `VmHWM` of `/proc/self/status` (the process's own high-water RSS; a child's
    memory is never in it), plus a 0.5 s `VmRSS` timeline with marks at the start and end of each
    MMseqs2 command and each builder phase.
  - **MMseqs2:** each command is reaped with `os.wait4` (wall, CPU). Its peak is the largest
    `VmHWM` of any process in its tree, sampled every 0.5 s from `/proc`. `ru_maxrss` from `wait4`
    is not used for it: Linux carries the forking parent's RSS into a child's `ru_maxrss` at exec,
    so in the first "before" run every short MMseqs2 command reported the Python process's
    1,996.9 MB. A sampled peak can miss a spike shorter than 0.5 s.

## Results at 2,053,075 proteins

Three runs: "before" twice (units `b24-prof-before`, `b24-prof-before2`) and "after" once
(`b24-prof-after`).

| | Before, run 1 | Before, run 2 | After |
|---|---|---|---|
| Python process, peak RSS (`VmHWM`) | **8,920.9 MB** | **8,915.1 MB** | **997.2 MB** |
| Python process, CPU | 136.4 s | 202.7 s | 121.6 s |
| wall | 1,317.1 s | 1,368.1 s | 1,302.4 s |
| MMseqs2 `cluster` NR90, peak | 4,620.9 MB (`wait4`) | 3,412.8 MB (sampled) | 3,274.5 MB (sampled) |
| MMseqs2 `cluster` NR50, peak | 3,016.3 MB (`wait4`) | 3,165.7 MB (sampled) | 5,312.0 MB (sampled) |
| MMseqs2, CPU (all commands) | 16,985.5 s | 16,846.7 s | 17,222.3 s |
| sha256 of `cluster_member`, `cluster_member_aln`, ledger | `8bc4461d…`, `709d9ba1…`, none | the same three | the same three |

The ledger's sha256 is `8a02d3fe…`. Run 1 failed after writing both tables, in the profiler's own
wrapper around the tool-version call, so it has no ledger; its Python peak was already reached
(the timeline below). Run 1 used the `wait4` peaks; its two `cluster` values lie above the floor
the forking parent sets (1,996.9 MB), so they are MMseqs2's own. Runs 2 and 3 used the sampled
`VmHWM`. The same commands gave 3,016 to 5,312 MB per `cluster` across the runs, so a sampled value
is a lower bound and MMseqs2's peak at this size is at least 5.3 GB. The "before" Python peak
agrees with the ladder's 9,570.9 MB for `rL.2` (`RUSAGE_SELF` of the step process,
`results/B/cluster.json`) to 7%.

MMseqs2's other commands, after (sampled `VmHWM`, MB): `createdb` 184.7, `createsubdb` 399.5,
`mergeclusters` 232.9; NR90 / NR50 `createtsv` 351.9 / 293.4, `align` 911.4 / 1,162.2,
`convertalis` 539.8 / 518.4. Run 2 gave the same to within 31%.

### Where the memory went before

Python `VmRSS` and `VmHWM` at the marks of the "before" run:

| Time s | Phase | VmRSS MB | VmHWM MB |
|---|---|---|---|
| 14.8 | proteins read into a dict, FASTA written | 1,388.0 | 1,996.9 |
| 24.5 to 1,124.0 | MMseqs2 `createdb`, both `cluster`s, `mergeclusters` | 1,388 | 1,996.9 |
| 1,127.7 | NR90 TSV parsed into a dict | 1,919.7 | 2,060.1 |
| 1,158.1 | NR90 alignment lines parsed into lists | 3,025.2 | 3,175.3 |
| 1,209.2 | NR50 TSV and alignment lines parsed | 3,509.9 | 3,563.9 |
| 1,216 to 1,290 | row dicts for both tables built and sorted | 4,674 → 8,612 | |
| 1,309.5 | tables written | 8,613.9 | **8,920.9** |

So the Python side holds the memory: the protein dict (about 1.4 GB), the parsed TSV lines (about
2.1 GB more) and above all the row dicts (about 5.1 GB more), against at most 5.3 GB for MMseqs2.
That is 4.3 KB of Python peak per protein (*derived*: 8,920.9 MB / 2,053,075).

### What is left after

| Time s | Phase | VmRSS MB at its end | VmHWM MB |
|---|---|---|---|
| 2.8 | imports | 203.4 | 321.5 |
| 15.1 | FASTA written by merging the stages' sorted files | 637.5 | 637.5 |
| 1,085 to 1,206 | MMseqs2; the four TSVs converted to flat Parquet | 316 to 357 | 650.1 |
| 1,231.6, 1,256.7 | per level: sorts by member and by (target, line), merge with the proteins | 921.2, 994.1 | 994.1 |
| 1,269.8, 1,282.2 | per level: sort by (rep, member), cluster uids hashed | 752.0, 697.1 | 997.2 |
| 1,291.1, 1,300.6 | both tables sorted and written | 876.7, 851.1 | 997.2 |

The rest is fixed in size: DuckDB's buffers (capped at `DUCK_MEMORY_LIMIT` = 2 GB; its sorts
spill), one input row group of `protein_seq` per stage (at most 122,880 rows, AGENTS §4.6), a
writer row group, and batches of 16,384 rows. None of it holds a list of all proteins, members or
alignments.

**The remaining peak of the clustering step is MMseqs2's own**: at least 5,312.0 MB (sampled, the
NR50 `cluster` command) at 2.05 M proteins, five times the Python process's 997.2 MB. This card does not
change MMseqs2 (its command lines, parameters and thread counts are those of `121503d`; the
outputs above are byte-identical). How MMseqs2's peak grows with the release is not measured
here: the ladder's 46.1 GB at 10.63 M proteins and the 410 GB projection for 95 M (README,
`rL3.md`) were the Python process's `RUSAGE_SELF`, not MMseqs2's. Bounding MMseqs2 itself
(`--split-memory-limit`, splitting the database, or `mmseqs clusterupdate`) changes its
parameters or the cluster identity rules and is a question for the PI, not for this card.

## The memory test on synthetic proteins

`tests/r2/test_micro_s0_cluster.py::test_b24_python_peak_does_not_grow_with_protein_count`
builds a K-copy protein set (5,000 random base proteins of 60 to 149 residues, each copied K times
with 4% of its residues substituted, seed 24) and runs `cluster.run` in a child process at
K = 30 (150,000 proteins) and K = 80 (400,000), MMseqs2 at 4 threads, DuckDB capped at 64 MB.
It compares the child's own `VmHWM`. Run with `sh scripts/pod_pytest.sh B-24
tests/r2/test_micro_s0_cluster.py -k b24` (after) and the same with `BIOBOX_TREE` set to a
`121503d` checkout carrying the new test (before):

| | 150,000 proteins | 400,000 proteins | Growth |
|---|---|---|---|
| Before | 749.3 MB | 1,742.8 MB | 993.5 MB (4.0 KB per added protein, *derived*) |
| After | 308.8 MB | 350.0 MB | 41.1 MB |

The test's bound is 96 MB of growth (DuckDB's 64 MB cap plus 32 MB). Before the test read
`VmHWM`, it read `ru_maxrss` of the child, which carried the pytest process's RSS (the protein
set it had just built) and so grew with K for the streaming code too (340 to 517 MB).

## A DuckDB limit found on the way

DuckDB 2.0.0 (the image's) sorts out of core, but reading a sorted result needs memory of its
own. A probe (`/tmp/if2453/b24/duck_sort.py`: 1,000,000 rows of two 32-character strings, one
thread, `preserve_insertion_order`) completed `ORDER BY` and failed to drain the result under a
32 MB limit (`failed to pin block of size 256.0 KiB (30.2 MiB/30.5 MiB used)`); it drained under
64, 128 and 256 MB (process peaks 206, 270 and 344 MB). Two such sorted results read at once
failed at 64 MB and passed at 128 MB. So the builder never reads two sorts at once: each sort is a
`COPY (... ORDER BY ...) TO` a flat file, one at a time, and the merges read those files with
pyarrow. The first rewrite, which joined the alignment lines to the members with a window
function, ran out of a 32 MB cap at 120,000 proteins in the test.

## Fixture identity

On the micro-S0 fixture plus a mutated lambda (`tests/r2/test_micro_s0_cluster.py::_build`), the
four files `cluster.run` writes (`cluster_member`, `cluster_member_aln`, the release's clustering
ledger and its `release_def`) have the same sha256 before and after, at 1 and at 4 threads
(`/tmp/if2453/b24/fixture_hash.sh`, unit `b24-fixture`).
