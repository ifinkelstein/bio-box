# Contig labels and cluster membership by their lookup keys (card B-27), 2026-11-10

Before this card two lookups could not prune. `contig_occurrence` is looked up by `contig_uid` (a contig's
accession and header) but sorted by (`l0_object_sha256`, `record_ordinal`). `cluster_member` is looked up by
`member_uid` (a protein's NR90 and NR50 cluster) but sorted by (`level`, `cluster_uid`, `member_uid`). The hit
and locus exports, `nbr summarize` and the catalog's `contig_label` view therefore read every row group of
those tables for any number of keys.

What changed:

- Two lookup tables, written by the builders of their sources in the same step (`build-stage` has no new step):
  - `contig_label` (`l0seq/contig_label/stage=<S>`), written by `ingest`. One row per contig of the stage,
    sorted by `contig_uid`: the accession and header of the contig's first record by (`l0_object_sha256`,
    `record_ordinal`), that record's `l0_object_sha256` and `record_ordinal`, and `n_records`. It is built
    from one DuckDB sort of the stage's records that spills, read a batch at a time.
  - `member_cluster` (`l2/member_cluster/release=<R>`), written by `cluster.run` from the `cluster_member`
    files it has just written (`cluster.write_member_cluster`). It holds `member_uid`, `level`,
    `cluster_uid`, `rep_uid`, sorted by (`member_uid`, `level`), from one DuckDB sort that spills.
- The catalog's `contig_label` view, the exports (`hits_tsv`, `hits_faa`, `locus_tables`) and `nbr.summarize`
  read the new tables. `contig_dims` reads `contig_dim` only and is unchanged.
- `open_catalog`, `export.open_store` and `nbr.summarize` set DuckDB's `dynamic_or_filter_threshold` to the row
  groups of their largest table, at least 50 and at most 4,096 (`catalog.set_join_key_filter`). Section 2 says
  why.
- The dead `manifest=None` glob branch of `nbr._store_files` is gone.

All numbers are *measured* on wilkcomp02 on 2026-11-10 in `bio-box-test_74b3e05.sif` (DuckDB 2.0.0), under
`systemd-run --user -p MemoryMax=32G -p MemorySwapMax=0 -p CPUQuota=400% nice -n 19` (16G for the fixture and
the ladder), unless labelled otherwise. The pod was shared with other agents' builds during every run, so wall
times are noisy. Row-group counts are not.

## 1. Synthetic stores at 100x and 1,000x, before and after

Driver: `/tmp/if2453/b27/exp/scale.py <tree> <out> <copies>`, run by `/tmp/if2453/b27/exp/run.sh` for the tree
at `plan/task-cards` (6225c32, units `b27-before`, `b27-before-fx`) and this branch (711134d, unit
`b27-after`). It runs `scripts/stress/synth_scale.py` unchanged, with the row-group caps of
`tests/r2/test_scale_layout.py` scaled down:

- row tables: 7,680 rows;
- the contig-level tables, now including `contig_occurrence` and `contig_label`: 32 rows;
- `token_dict`: 960 rows.

With the branch it then writes `member_cluster` from the synthetic `cluster_member` with the builder's function
and rebuilds the manifest and catalog. Each catalog lookup runs on a fresh `open_catalog` connection
(4 threads), "cold", then again, "warm". "rg" means row groups read of the table's total
(`catalog.profile_query`).

| Table | 100x rows, rg | 1,000x rows, rg | Bytes per row at 1,000x |
|---|---|---|---|
| `contig_occurrence` | 300, 10 | 3,000, 94 | 102.6 |
| `contig_label` | 300, 10 | 3,000, 94 | 75.0 (synthetic headers; section 3 has real ones) |
| `cluster_member` | 170,800, 23 | 1,708,000, 223 | 19.7 |
| `member_cluster` | 170,800, 23 | 1,708,000, 223 | 40.6 |

`member_cluster` is about twice the bytes of `cluster_member` here. The synthetic `cluster_member` is sorted by
cluster, so its member uids come in runs that compress. On real clusters the two are close (section 3).

### Catalog lookups

"Before" reads the old `contig_label` view (over `contig_occurrence`) and `cluster_member`. "After" reads the
new `contig_label` view and `member_cluster`. The keys come from a joined table. "Adjacent" means 100 keys next
to each other in key order. The 10,000 contig keys are all 3,000 contigs at 1,000x and all 300 at 100x.

| Lookup, 1,000x | Before rg | Before s (cold/warm) | After rg | After s (cold/warm) |
|---|---|---|---|---|
| label, literal key | 94/94 | 0.038/0.009 | 1/94 | 0.044/0.006 |
| label, 1 joined key | 94/94 | 0.035/0.012 | 1/94 | 0.043/0.008 |
| label, 100 random keys | 94/94 | 0.045/0.017 | 63/94 | 0.028/0.009 |
| label, 100 adjacent keys | 94/94 | 0.049/0.012 | 4/94 | 0.019/0.006 |
| label, all 3,000 contigs | 94/94 | 0.059/0.012 | 94/94 | 0.049/0.013 |
| membership, literal key | 223/223 | 0.058/0.032 | 1/223 | 0.041/0.010 |
| membership, 1 joined key | 223/223 | 0.058/0.031 | 1/223 | 0.022/0.011 |
| membership, 100 random keys | 223/223 | 0.084/0.044 | 84/223 | 0.048/0.035 |
| membership, 100 adjacent keys | 223/223 | 0.076/0.039 | 1/223 | 0.027/0.012 |
| membership, 10,000 random keys | 223/223 | 0.082/0.067 | 223/223 | 0.107/0.077 |

At 100x the pattern is the same:

- label: 10/10 before for every key set; after, 1/10 for a literal or 1 key, 4/10 for 100 adjacent keys, and
  10/10 for 100 random keys and for all 300 contigs;
- membership: 23/23 before; after, 1/23 for a literal, 1 key or 100 adjacent keys, and 23/23 for 100 random
  keys and for 10,000 keys.

100 random keys in 223 row groups lie in about 223 x (1 - e^(-100/223)) = 80 of them (*derived*). The 84 read
are the row groups that hold a key. At 100x the 100 keys outnumber the 23 row groups, so reading 23 is less than
one row group per key.

### The exports' own lookups

These are `export._annotate_contigs` (`_lab`) and `export._annotate_proteins` (`_cl`), on `open_store`'s
connection with every statement profiled, for k contigs and k proteins.

| k, 1,000x | `_lab` before | `_lab` after | `_cl` before | `_cl` after |
|---|---|---|---|---|
| 1 | 94/94, 0.041 s | 1/94, 0.020 s | 223/223, 0.046 s | 1/223, 0.017 s |
| 100 | 94/94, 0.042 s | 63/94, 0.032 s | 223/223, 0.081 s | 84/223, 0.051 s |
| 10,000 (3,000 contigs) | 94/94, 0.128 s | 94/94, 0.035 s | 223/223, 0.094 s | 223/223, 0.120 s |

The plan for 10,000 keys is a `HASH_JOIN` of type `SEMI` whose probe side is the `READ_PARQUET` of the lookup
table (the `EXPLAIN` in `ladder.json`, section 3). The table has fewer row groups than there are keys, so DuckDB
pushes the keys' min and max and reads the table front to back in key order. That is the sorted merge the card
asks for: every row group is read once, and fewer than one per key.

### `nbr summarize`, 100 anchor proteins, window 5g, levels pfam, nr50, nr90

| | Before | After |
|---|---|---|
| neighbors' clusters (`_k`, by member) | `cluster_member` 112/223 | `member_cluster` 223/223 |
| clusters' members and representatives (by cluster) | `cluster_member` 112/223 per scan | `cluster_member` 112/223 per scan, unchanged |
| wall, profiled | 1.75 s | 1.94 s |

The neighbor lookup by member now reads `member_cluster`. Here it has about 1,400 neighbor keys, more than the
223 row groups, so it reads the whole table. The old scan read half of `cluster_member`, because `level` leads
that table's sort key and the query names one level. Both reads are under one row group per key. For fewer
neighbors than row groups the new lookup reads only the row groups that hold them. The trade-off for wide
summaries is under "Open" below.

### End to end: same bytes before and after

| Export, 1,000x | sha256 before = after | Lines | Wall before | Wall after |
|---|---|---|---|---|
| `hits_tsv`, 10,000 hit proteins | 7add6af1... | 10,001 | 0.79 s | 0.98 s |
| `locus_tables`, 100 anchor proteins, 5g | e2be038f... | 1,396 | 0.32 s | 0.39 s |
| `nbr.summarize`, 3 levels, `top_n` 0 | 4e5a4d84... | (JSON) | 1.24 s | 1.44 s |

The outputs are also identical at 100x (`hits_tsv` 070c9ac8..., `locus_tables` 61971937..., `summarize`
007d7f2e...).

At this size the end-to-end times are dominated by setup, not by these lookups. The table above shows each
lookup getting faster. The synthetic build itself took 338 s before and 316 s after; writing `member_cluster`
added 3.3 s at 1,000x.

## 2. Why the join filter is raised: two probes

### `probe.py`

`/tmp/if2453/b27/exp/probe.py`: 16,384 random 16-byte keys, sorted, in 16 row groups of 1,024.

| Query | Plain view over the file | View with `GROUP BY k` |
|---|---|---|
| `k = literal` | 1/16 | 1/16 |
| join or `IN` with 1 key | 1/16 | 16/16 |
| 10 keys in 10 row groups | 10/16 | 16/16 |
| 100 adjacent keys | 1/16 | 16/16 |
| 100 random keys | 16/16 | 16/16 |
| 10,000 keys | 16/16 | 16/16 |

The same holds when the `GROUP BY` sits over a `UNION ALL` of two files (8/8 and 12/12 for every joined form).

DuckDB 2.0 pushes no join filter through an aggregate. So:

- The catalog's `contig_label` view is a plain projection of the `contig_label` files when the release has one
  stage. With several stages it must merge rows with a `GROUP BY`; there a literal key still prunes, but keys
  from a joined table read the whole label table (still smaller than `contig_occurrence`).
- The exports filter the label rows by key before they merge them, so they prune either way.
  `test_b27_contig_label_view_merges_the_release_stages` checks the merged answer.

### `probe2.py`

`/tmp/if2453/b27/exp/probe2.py`: 2^20 random keys, sorted, in 1,024 row groups. Each column gives row groups
read (seconds, `IN` form).

| Random keys | Row groups holding them | Default threshold (50) | Threshold 100,000 |
|---|---|---|---|
| 1 | 1 | 1 (0.06 s) | 1 (0.04 s) |
| 10 | 10 | 10 (0.08 s) | 10 (0.06 s) |
| 60 | 59 | 1,015 (0.17 s) | 59 (0.14 s) |
| 100 | 97 at default, 94 raised (different draws) | 981 (0.13 s) | 94 (0.07 s) |
| 1,000 | 634, 631 | 1,022 (0.14 s) | 631 (0.26 s) |
| 10,000 | 1,024 | 1,024 (0.14 s) | 1,024 (0.48 s) |

Up to `dynamic_or_filter_threshold` keys, the hash join hands the scan the list of keys, and the scan reads only
the row groups whose statistics hold one. Above the threshold it hands only the keys' min and max, so random
keys read nearly every row group.

The list costs time on every row read: 0.48 s against 0.14 s at 10,000 keys, when every row group is read
anyway. So the threshold is set to the row groups of the connection's largest table, capped at 4,096:

- With no more keys than that, a lookup reads the row groups that hold a key, at most two per key (AGENTS §4.4).
- With more keys, the table has fewer row groups than there are keys, and the cheap min/max filter is used.

The cap covers the lookup tables at GTDB scale: 190 M / 122,880 = 1,547 and 111.5 M / 122,880 = 908 row groups
(*derived*).

## 3. The real store `rL.3` (import ladder), read only

Driver: `/tmp/if2453/b27/exp/ladder.py /tmp/if2453/ladder/root /tmp/if2453/b27/ladder_out` (unit `b27-ladder`,
this branch). It builds both tables from `rL.3`'s files at the real caps (122,880 rows per row group), into
`/tmp/if2453/b27/ladder_out`.

| | Rows | Row groups | Bytes | Bytes per row |
|---|---|---|---|---|
| `contig_occurrence`, stages S1, P1, P2 | 39,633 | 3 | 2,611,395 | 65.9 |
| `contig_label`, the same stages | 32,815 | 3 | 1,571,095 | 47.9 |
| `cluster_member`, rL.3 | 21,241,368 | 173 | 762,906,723 | 35.9 |
| `member_cluster`, rL.3 | 21,241,368 | 173 | 859,640,737 | 40.5 |

Build cost:

- `member_cluster`: 56.3 s on one thread, peak RSS 2,171 MB. That is the process with DuckDB's 2 GB cap; the
  RSS was 196 MB before the sort.
- `contig_label`: 0.04 to 0.15 s per stage.

The rows equal the card's figure of 21 M `cluster_member` rows at the 10 Gbp release. The code's size
estimates (`EST_LABEL_ROW`, `EST_MEMBER_CLUSTER_ROW`) were set to these measurements. They only cut parts.

Lookups use the `IN (SELECT ...)` form of the exports, on a plain DuckDB connection with 4 threads. The table
gives row groups read (cold seconds).

| Lookup, rL.3 | `cluster_member`, default | `member_cluster`, default | `member_cluster`, threshold 173 |
|---|---|---|---|
| literal key | 173/173 (0.51 s) | 1/173 (0.04 s) | 1/173 (0.13 s) |
| 1 joined key | 173/173 (0.63 s) | 1/173 (0.03 s) | 1/173 (0.05 s) |
| 100 random keys | 173/173 (0.85 s) | 171/173 (0.30 s) | 74/173 (0.15 s) |
| 1,000 random keys | 173/173 (0.38 s) | 173/173 (0.32 s) | 173/173 (0.35 s) |
| 10,000 random keys | 173/173 (0.39 s) | 173/173 (1.40 s cold, 0.33 warm) | 173/173 (0.34 s) |

100 random keys in 173 row groups lie in about 173 x (1 - e^(-100/173)) = 76 of them (*derived*); 74 were read.

The `contig_label` tables of `rL.3` are one row group per stage (3 in all), as is `contig_occurrence`, so at
this size neither label table can show pruning. The synthetic stores in section 1 show it.

## 4. Projection (*derived*)

Bytes = rows x the bytes per row measured on `rL.3` (section 3).

| Table | 10 Gbp release (`rL.3`) | GTDB scale |
|---|---|---|
| `member_cluster`, 40.5 B/row | 21.24 M rows: 860 MB (*measured*) | 190 M rows (card: 2 levels x 95 M proteins) x 40.5 = 7.7 GB |
| `contig_label`, 47.9 B/row | 32,815 contigs: 1.6 MB (*measured*) | S5, 24.95 M contigs: 1.2 GB; S6 cumulative, 111.5 M contigs: 5.3 GB |

Notes on the projection:

- Contig counts are from `sources_and_sizing.md`, GTDB R226 metadata (*measured*, 2026-11-09).
- A contig that occurs in several stages has one row in each.
- `member_cluster` is about 1.13 times the bytes of `cluster_member` (40.5 / 35.9).
- `contig_label` is about 0.60 times the bytes of `contig_occurrence` (1,571,095 / 2,611,395); label bytes
  grow with header length.

## 5. Fixture: exports byte-identical, files changed

Driver: `/tmp/if2453/b27/exp/fixture.py <tree> <out>`, run by `runfx.sh` for both trees. It builds micro-S0 as
`tests/r2/test_micro_s0_l3.build` does, with real MMseqs2 clusters and the test HMM library, plus
`genome_input` and `contig_dim`, then writes the release manifest. It runs these exports and the
`contig_label` view, 49 outputs in all:

- `hits_tsv` and `hits_faa` for 128 hit proteins: every protein with a domain hit, and every 7th other;
- `locus_tables` for the first 40 of them (window 5g);
- `nbr.summarize` at windows 5g, 10g and 3kb, as NDJSON and TSV.

Results:

- All 49 outputs are byte-identical before and after.
- Two files are new: `l0seq/contig_label/stage=S0/part-0.parquet` and
  `l2/member_cluster/release=rS0.l3/part-0.parquet`.
- One file changed: `manifests/release=rS0.l3.json`, because it lists the two new files. Its hash went from
  `04b91a217f2ffe2958659ab13c61a6b3` to `2dd511dbe58d1a9d179627c4b454ddc6`.
- Every other file under the data root (27 of 29) is byte-identical.

## 6. Tests (`tests/r2/test_scale_layout.py`, 100x)

These failed first on the pod against `plan/task-cards` with the new tests:

- the label lookups read `contig_occurrence` 10 of 10;
- the membership tests errored: `module 'biobox.l2.cluster' has no attribute 'write_member_cluster'`;
- the two-stage test found no `contig_label` files.

They pass on this branch (45 passed). They check:

- **The new tables against brute-force oracles.** Each contig's label is its first `contig_occurrence` record
  by (sha256, ordinal). Each (member, level) of `cluster_member` appears in `member_cluster` exactly once, with
  the same cluster and representative.
- **Reads for a literal key and for 40 or 100 joined keys,** for labels, membership and the exports' own
  lookups:
  - nothing is read from the old tables;
  - at most 2 row groups per key and no more than the keys' spread;
  - exactly the row groups that hold a key when the keys fit the join filter list;
  - 40 keys from both ends of the key order read only the row groups that hold them (2 for members), fewer
    than the spread, which is the whole table.
- **The merged view of a two-stage release.**

## Open

- **`nbr summarize` with many neighbors.** It reads all of `member_cluster` where it used to read the half of
  `cluster_member` that holds one level (112/223 against 223/223 at 1,000x, about 1,400 neighbor keys). With
  fewer neighbors than row groups it reads less than before. Sorting `member_cluster` by (`level`,
  `member_uid`) would keep both properties at 2 row groups per member for the exports, but the card specifies
  `member_uid` order. It is left as the card says, and reported.
- **Multi-stage catalogs.** Joined keys on the catalog's `contig_label` view scan the label table (section 2).
  The exports are not affected.
