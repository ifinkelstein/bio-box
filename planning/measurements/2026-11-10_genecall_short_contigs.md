# Gene calling of short contigs, 2026-11-10 (card B-32)

Every number here is *measured* on wilkcomp02 (112 cores, about 96% idle during the runs; load
average 7 to 9 from other jobs) on 2026-11-10, in the pinned image
`bio-box-test_74b3e05.sif` (Python 3.14.10, pyrodigal 3.7.1, pyrodigal-gv 0.3.2, pyarrow 26.0.0,
duckdb 2.0.0), unless it is labelled *derived* or *estimate*. Each run is its own
`systemd-run --user --unit=b32-<name> -p MemoryMax=32G -p MemorySwapMax=0` unit, at `nice -n 19`, with
`singularity exec --cleanenv --env PYTHONPATH=<tree>/src:/tmp/if2453/pilot_s1/pydeps --env OMP_NUM_THREADS=2
--bind /tmp/if2453 --bind /sys/fs/cgroup:/sys/fs/cgroup:ro` (the prefix `scripts/pod_pytest.sh` uses), at
most 8 workers plus the parent. The scripts are in `/tmp/if2453/b32/scripts/` on the pod and are
described in the appendix. They are not in the repository (the card lists no script file).

## Inputs

- **Slice**: 100,000 contigs, 61,953,033 bp, a seeded sample (`random.Random(20261110)`) of the
  1,475,957 default contigs of the B-18 CAMI II rerun store (`/tmp/if2453/b18/cami/root`, stage `SM`,
  read only). They were written to a FASTA file and ingested into a fresh root,
  `/tmp/if2453/b32/slice/root` (one shard), each contig its own prokaryote genome, as
  `pilot_cami.py` does. By length: 91,594 contigs under 1 kb (52.4% of bases), 8,102 from 1 to 10 kb
  (27.5%), 252 from 10 to 40 kb (6.9%), 28 from 40 to 100 kb (2.5%) and 24 of 100 kb or more (10.6%,
  single mode). Gene calling finds 133,754 genes and 133,555 distinct proteins.
- **Fixture**: the micro-S0 fixture (lambda, T4, *M. genitalium*; 854 genes).
- **Golden set**: `/tmp/if2453/golden` as `scripts/stress/build_s0.py` builds its genome list
  (33,712 genes).

Each run copies the base root, builds the genome list from `genome_input` outside the timed region,
and times `genecall.run(root, "SM", genomes, workers=W)`. "Wall" is the whole `run` call. "Publish" is
its `_publish` (the DuckDB sort and the writes; wrapped by the script, not changed). "Calling" is
wall minus publish.

## Before: the profile (tree `plan/task-cards` at fc3f0ec)

| Slice | wall | publish | CPU (cgroup) | peak RSS |
|---|---|---|---|---|
| 1 worker | 145.9 s | 2.0 s | 146.0 s | 482 MB |
| 8 workers (threads) | 130.5 s | 2.5 s | 281.6 s | 572 MB |

That is 1.12x with 8 workers (*derived*), the B-18 finding (1,989 s against 2,322 s on all of CAMI).
Per-thread CPU at 8 workers (from `/proc/<pid>/task/*/stat`, sampled every 5 s): the calling thread
used 89.4 s and the eight worker threads 22.8 to 26.3 s each, 192 s in all. One thread did all the
work in 146 s, so the threaded run spent about twice the CPU for the same work (*derived*).

cProfile of the calling thread, 1 worker (189.6 s under the profiler):

| Where | tottime or cumtime | calls |
|---|---|---|
| `_genes` (`find_genes`, C, releases the GIL) | 93.5 s | 233,754 |
| `_sequences` (`SeqStore.get`), in all | 81.9 s cum | 100,000 |
| of which `KeyLookup.get` (index, locator, locator again) | 33.0 s cum | 300,000 |
| of which `_Values.__getitem__` (one pyarrow scalar per bisect probe) | 15.2 s | 10.0 M |
| of which `zlib.decompress` (BGZF blocks) | 30.6 s | 94,060 |
| of which the `.fai` bisect in `bgzf.py` | 10.2 s cum | 100,000 |
| `_finder` (a meta-mode `GeneFinder` per genome) | 3.0 s | 100,000 |
| `ids.feature_uid` + `ids.protein_uid` | 2.3 s | 133,754 each |
| scratch rows (`_Scratch.add`) | 2.1 s | 367,508 |
| `_publish` (sort and write) | 3.0 s | 1 |
| pickling | none (threads) | |

With 8 worker threads, cProfile (3.14 profiles every thread) puts the calling thread's sequence
reads at the same cost, and `SimpleQueue.get` / `Future.result` waits at 153 s cumulative. **What
bounds it**: about 0.5 ms per contig of Python under the GIL. Most of it is reading the sequence
(`KeyLookup` and the BGZF reader, *derived*: 1.46 ms per contig in all unprofiled, less 0.93 ms in
`find_genes`). It ran in the one calling thread, and the per-gene Python work ran in the workers
under the same GIL. Model selection (`_finder`, 3%), uid hashing (1%), the scratch writer (1%) and
the sort (2%) are not the bound.

## The fix (`src/biobox/l1/genecall.py`)

1. Consecutive genomes shorter than `THREAD_MIN_BASES` (10 kb) are grouped into batches of at most
   `batch_bases` bases (default `BATCH_BASES` = 1 Mbp, chosen). Contig lengths come from
   `contig_index`, one key lookup per contig in the parent. A batch goes to a worker *process* that
   opens its own `SeqStore`, reads its own sequences and returns two Arrow record batches (genes,
   protein candidates) instead of per-gene dicts. The parent appends them to the scratch files.
2. A genome of 10 kb or more is a task alone in a worker *thread*, as before the card. Worker
   processes start only when the first batch of short genomes comes. Each costs about 2.8 s of
   imports: `python -X importtime -c "import biobox.l1.genecall"` gives 2.81 s, of which
   `pyrodigal_gv.meta` 1.19 s, duckdb 0.51 s and pyarrow 0.39 s.
3. With workers, results are taken in the order they finish, at most `workers` batches at once. The
   output does not depend on that order: `_publish` sorts every table on a key that is unique when
   no contig is in two genomes (MR-01 pins this).
4. The workers are `python -c` subprocesses with pickles on stdin and stdout, and a reader thread
   per worker. `multiprocessing` is not used: its spawn and forkserver methods re-run the caller's
   main script in each worker. The first try at this card crashed that way (`prof.py`, which has no
   main guard), and `build_s0.py`, `pilot.py` and `pilot_cami.py` have none either. Fork is unsafe in
   a process with threads.

Two intermediate versions were measured and dropped:

| Version | slice 8 workers | fixture 8 workers (median of 5) | golden 8 workers (median of 3) |
|---|---|---|---|
| fc3f0ec (threads) | 130.5 s | 1.07 s | 7.88 s |
| processes for every batch (3741420) | 30.2 s | 3.88 s | 9.81 s |
| plus long genomes in threads, results in input order | 53.7 s | | |
| plus results as they finish (cfc76e9; final) | see below | | |

Processes for everything made the complete genomes slower: startup imports dominate the fixture,
and the golden set is 25% slower. With long genomes in threads but results handed back in input
order, a 100 kb+ single-mode contig at the head of the queue held back every batch behind it, and
the processes idled.

## After: the profile (merged tree, task/B-32 at f379ef2; `plan/task-cards` 9edb1e8 as the base)

Interleaved runs on the slice (`cmp2.sh fslice ... 2`: base w1, base w8, new w1, new w8, twice):

| Slice | base 9edb1e8 (threads) | new f379ef2 | 
|---|---|---|
| 1 worker, wall (2 runs) | 172.4, 150.8 s | 165.9, 159.8 s |
| 8 workers, wall (2 runs) | 129.8, 129.2 s | 32.42, 32.37 s |
| 8 workers, publish | 2.4, 2.5 s | 2.6, 2.8 s |
| 8 workers, CPU (cgroup) | 280.6, 280.1 s | 225.6, 224.9 s |
| 8 workers, peak RSS parent / largest worker | | 551 / 322 MB |

**Speedup, 8 workers against 1 (*derived*, medians): 162.88 / 32.39 = 5.03x on the new tree;
161.59 / 32.39 = 4.99x against the base tree's 1 worker.** The card asks for at least 4x. One-thread
runs vary from 144 to 178 s across the day on both trees (fc3f0ec 145.9 and 143.8 s, this tree
150.4 to 178.2 s, each run measured). That spread is the host, not the code: the profiled
one-thread runs of the new code (169.2 to 180.5 s, three runs) were faster than the base's (189.6 s). The 8-worker runs
vary by under 1%.

cProfile, 8 workers (f379ef2, `prof.py ... 8 --profile`; 34.9 s under the profiler, 1 worker 176.0 s;
parent process only, all its threads):

| Where (parent) | time | calls |
|---|---|---|
| `_genome_bases` (the length lookups), cum | 22.6 s | 100,000 |
| of which `KeyLookup.get` | 14.8 s cum | 100,912 |
| reader threads blocked in `pickle.load` (waiting, not CPU) | 31.5 s cum | 311 |
| `pickle.dump` of batches | 0.23 s | 303 |
| long genomes in threads (`_genes`) | 7.3 s | 12,377 |
| `_publish` | 3.1 s cum | 1 |

(The same run on cfc76e9, before the merge, gave 35.3 s and 28.6 s in `_genome_bases`.)

Unprofiled (fslice r1 and r2 at 8 workers): the parent's calling thread used 15.5 to 15.8 CPU-s,
the three busiest long-genome threads 8.3 to 9.1 s each, and the eight worker processes 15.4 to
21.9 CPU-s each (181.8 and 181.9 s in all). The whole cgroup used 225 CPU-s against 160 to 165 s for
the same work on one thread. **What bounds it now**: the worker processes are about 65% busy
(*derived*: 19.5 CPU-s mean of 29.7 s calling). The rest is process start (imports, about 2.8 s
each, 22 CPU-s for eight), the tail of the last batches, and the parent's length lookups at about
0.1 ms per key in `KeyLookup.get`. Those lookups would cap the parent at roughly 10,000 contigs/s,
about 15 workers' worth at 667 contigs/s each (*estimate*, from the profiled 22.6 to 28.6 s per 100,000,
discounted for the profiler).

## Complete genomes do not regress

Interleaved runs (`cmp2.sh`; base = `plan/task-cards` 9edb1e8, new = task/B-32 f379ef2):

| | base | new | change (*derived*) |
|---|---|---|---|
| fixture, 1 worker, median of 5 | 1.92 s (1.46 to 2.39) | 1.95 s (1.51 to 2.15) | +2% |
| fixture, 8 workers, median of 5 | 0.91 s (0.85 to 1.08) | 0.94 s (0.82 to 1.58) | +3% |
| golden, 1 worker, median of 3 | 31.42 s (31.40 to 31.88) | 31.77 s (31.72 to 32.26) | +1% |
| golden, 8 workers, median of 3 | 9.06 s (8.04 to 9.14) | 8.02 s (7.92 to 8.87) | -11% |

## Output bytes

Every run records the sha256 of each of the five gene-call output files (`feature`,
`feature_protein`, `protein_occurrence`, `protein_seq`, the genecall ledger). Every run of fc3f0ec,
9edb1e8, cfc76e9 and f379ef2, at 1 and at 8 workers, has the same five hashes
as the first fc3f0ec 1-worker run of its input: 50 fixture runs, 30 golden-set runs and 16 slice
runs (`/tmp/if2453/b32/runs/*/result.json`; checked with a one-line comparison of the `hashes`
fields).

## What is not fixed here (outside `genecall.py`)

- **BGZF read amplification on short contigs.** One block decompression per contig (94,060 for
  100,000 contigs, 30.6 s of the 189.6 s profiled single-thread run, 16%). A 620 bp contig costs a
  64 KiB block inflate, because genomes come in `genome_uid` order, which is random against shard
  order. Sorting short genomes by storage position, or a larger block cache, belongs to
  `SeqStore`/`bgzf.py`.
- **`KeyLookup.get` costs about 0.1 ms per key** (profiled), half of it one pyarrow scalar `as_py()` per
  bisect probe (`lookup.py` `_Values.__getitem__`). It is three lookups per contig in each worker's
  sequence read, and the parent's length lookup (above). A vectorized search on the fixed-size key
  column would cut both; that is `src/biobox/l0seq/lookup.py`, which card B-28's agent and this card
  do not own.

## Appendix: the scripts (on the pod, `/tmp/if2453/b32/scripts/`)

- `make_slice.py SRC_ROOT OUT_DIR 100000`: the seeded sample, FASTA, `ingest`, `write_genome_input`.
- `make_complete.py fixture|golden REPO OUT_DIR [GOLDEN_DIR]`: the fixture and golden base roots
  (stage `SM`).
- `prof.py BASE WORK WORKERS [--profile]`: copies `BASE/root`, times `genecall.run`, and writes
  `WORK/result.json` (wall, publish, cgroup `cpu.stat` CPU, getrusage, per-thread and
  per-descendant CPU from `/proc`, output sha256) and, with `--profile`, `WORK/prof.txt`.
- `run.sh NAME TREE ARGS...`: the image command above, with STOP-file check and done/failed markers.
- `seq.sh`, `cmp.sh`, `cmp2.sh`, `final.sh`: run sequences (interleaved base/new; reps as stated).
- `summ.py TAG...`: medians and the hash check per tag.
