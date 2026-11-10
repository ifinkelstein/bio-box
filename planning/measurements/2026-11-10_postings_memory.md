# Postings builder memory (card B-29, wilkcomp02, 2026-11-10)

Every number is *measured* on wilkcomp02 on 2026-11-10 (CST) unless it is marked *derived*. DuckDB
2.0.0 and pyarrow 26.0.0 inside `bio-box-test_74b3e05.sif`; every run in a `systemd-run --user` unit
named `b29-*` with `MemoryMax=16G` (32G for the byte checks), `MemorySwapMax=0`, nice 19, and the
builder's own DuckDB settings (one thread, `memory_limit` 2 GB, spill directory, insertion order
preserved). Scripts: `/tmp/if2453/b29/probe/` on the pod (not in the repository).

## Input

Release `rL.3` in `/tmp/if2453/ladder/root` (stages S1, P1, P2; tape version `t-ca596dd85040`),
built by the ladder's older checkout. The store was only read: `/tmp/if2453/b29/root` holds symlinks
to its `l3/tapes/tape_version=t-ca596dd85040` and `l3/tape_tokens/.../release=rL.3` directories, and
the builder wrote `release_def` and `token_postings` there. The current code reads the older layout
(its tapes have extra legacy columns; postings reads only the columns it needs).

| Item | Value |
|---|---|
| contigs | 32,815 (S1 27,258; P1 539; P2 5,018) |
| slots | 12,048,866 |
| postings | 25,822,709 (family 1,724,977; nr50 12,048,866; nr90 12,048,866) |
| tape row groups | 19 + 8 + 68, at most 131,068 slots each |
| tape-token row groups | 93, at most 131,071 slots each |

Card B-29 says "about 21 M slots"; the store has 12,048,866 slots, which equals the gene count of
the ladder README's input table (2,398,087 + 927,920 + 8,722,859). The card's number is off by 74%.

## Failure before the fix

`postings.run` from `plan/task-cards` (c285538) on that store, 2 GB DuckDB limit:
`OutOfMemoryException: could not allocate block of size 256.0 KiB (1.8 GiB/1.8 GiB used)` after
38.0 s wall, peak RSS 2,189.5 MB (`os.wait4`, `measure.py`).

## Which operator holds the memory

`probe.py` runs parts of the old statement on rL.3 and polls `duckdb_memory()` every 0.1 s from a
second cursor; the table gives the peak per tag. `EXPLAIN` shows the plan: a hash join of the tape
rows to the token rows (whole lists), an unnest, a materialized CTE `s` (ColumnDataCollection) scanned
by three `UNION ALL` branches (family, nr50, nr90), and one `ORDER BY` over the union.

| Statement (rL.3, 2 GB unless stated) | Result | ORDER_BY MB | COLUMN_DATA MB | HASH_TABLE MB | spilled MB |
|---|---|---|---|---|---|
| whole statement | OOM at 28.9 s | 1,899 | 1,292 | 340 | 981 (COLUMN_DATA only) |
| whole statement, 12 GB | ok, 90.3 s | 3,502 | 1,292 | 340 | 0 |
| join, unnest, CTE, union; no ORDER BY | ok, 6.5 s | | 1,292 | 340 | 0 |
| the CTE not materialized | OOM at 10.8 s | 1,193 | | 1,021 (one per branch) | 0 |
| only P2's tapes | OOM at 23.5 s | 1,884 | 901 | 320 | 672 (COLUMN_DATA only) |
| only P2's tapes, family branch last | ok, 46.2 s | 1,900 | 901 | 320 | 688 ORDER_BY, 666 COLUMN_DATA |
| the union written to one flat Parquet file, then ORDER BY over it | ok, 62.7 s for the sort | 1,895 | | | 1,515 ORDER_BY |
| `_check` (full outer join of whole list rows) | ok, 2.8 s | | | 340 | 0 |

The operator that runs out is the `ORDER BY`. Fed by the union over the materialized CTE, it grows to
the limit and never spills: zero ORDER_BY bytes reached the spill directory in each failing run.
It spills, and the statement completes, when the union's first branch is large. It also spills when
its input is a flat Parquet scan. In rL.3 the first branch (family) is small: 1.7 M of the 25.8 M
rows, and 26,485 of P2's 8.7 M slots. Synthetic stores with families on 75% of slots passed at
the same and at smaller scaled limits (13.1 M postings at 300 MB; 31.3 M postings at 2 GB, sorts
spilling 1.1 and 2.2 GB). The test store therefore puts families on 3 slots in 1,000. Why the sort's
spilling depends on the first branch is inside DuckDB and was not traced further. The hash join's
build side (the whole token lists, 340 MB here) and the CTE (1.3 GB) both grow with the release, so
the fix removes them as well as the sort-over-union.

## The fix

Postings are exploded with Arrow one tape-token row group at a time (at most 131,072 slots), beside
a cursor on each tape file in `contig_uid` order, into flat zstd Parquet runs in the private scratch
directory. One DuckDB statement then sorts the runs by (bucket, sort key), and the existing writer
streams it into the bucket and part files. `_check` joins list lengths, not lists.

| Run (rL.3, 2 GB DuckDB limit) | Wall s | Peak RSS MB | CPU s (user + sys) | Result |
|---|---|---|---|---|
| `plan/task-cards` builder | 38.0 | 2,189.5 | 30.8 + 7.1 | OutOfMemoryException |
| B-29 builder (commit c6e1272) | 105.4 | 2,417.4 | 90.4 + 17.3 | 25,822,709 rows, one bucket, one file |

Peak RSS is the process's, measured with `os.wait4` on a forked child (`measure.py`). It exceeds the
2 GB DuckDB cap by the Python, Arrow and Parquet-writer memory outside DuckDB. Rate: 245 k rows/s
(*derived*: 25,822,709 / 105.4 s), against 272 to 287 k rows/s for the old builder on `rL.1` and
`rL.2` (ladder README). The output file is 1,001,038,504 bytes, 38.8 bytes per row (*derived*;
`EST_POSTING_ROW` is 36, 7.7% under).

The processes show about 260 threads (`ps -o nlwp`) but use one core (100% CPU in `ps`; CPU time
about equal to wall time).

## Bytes are unchanged

- rL.3: the B-29 builder at 2 GB and the `plan/task-cards` builder with its DuckDB cap patched to 8 GB
  (`old_bigmem.py`) wrote the same file: sha256 `92d69f44695fcb4611da09701918dd030b0ed1ace428201d1ef36f5d82b9c66a`.
- Fixture and synthetic stores (`digests.py`): micro-S0 at bucket widths 0, 1, 2 and with parts cut at
  7 rows (widths 0 and 1), the B-10 synthetic store (246,967 postings) and the B-29 store
  (1,312,511 postings). Inputs were built once and postings were written with each tree. 771 files,
  and every sha256 is equal. The two JSON reports are byte-identical (sha256
  `cd3aae8cf24346bead040ef2bcb1672dcc70a6552eb374acafd249175044a9e3`).

## tokens.run

Not changed. It unnests only `protein_uids` before any join. Its joins have narrow build sides (the
canonical families, the cluster members, the ranked keys), and its list aggregate is the unordered
`list()` that spills. The old postings shape (whole list rows joined, then a union sorted from a
CTE) does not occur in it. On rL.3 it took 216 s and 2.3 GB peak RSS (ladder README, stage C).

## Commands

    sh probe.sh <tree> <mode> 2GB [root release stages producers]   # modes as in the table
    sh launch.sh rl3-old   /tmp/if2453/b29/base  /tmp/if2453/b29/root-old
    sh launch.sh rl3-fixed /tmp/if2453/b29/fixed /tmp/if2453/b29/root
    sh runimg.sh dig-before /tmp/if2453/b29/base  digests.py post <inputs> <dest>
    sh runimg.sh dig-after  /tmp/if2453/ci/B-29   digests.py post <inputs> <dest>
