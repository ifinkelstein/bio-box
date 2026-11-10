# `token_of` by stable key, and catalog connections (card B-34), 2026-11-10

Card B-34 has three parts, all in `src/biobox/catalog.py` (and one line of `src/biobox/l3/tokens.py`):

1. `token_of(level, key)` resolved a stable key with a scalar subquery on `token_dict`, which is sorted by
   `token_id`, so it read every row group of `token_dict` for a literal key and for joined keys.
2. `open_catalog` could not open a second connection to the same catalog file in one process.
3. On a plain `duckdb` connection, the sequence macros and `partners_of` failed with "Scalar Function with
   name bb_seq does not exist!" (`bb_slot_tokens` for `partners_of`), which does not say what to do.

What changed:

1. `token_dict` keeps its order (by `token_id`; the reverse lookup token -> key prunes on it, below) and
   gains min/max statistics on `key`. `open_catalog` shadows the stored `token_of` with a temporary macro over
   `bb_token_of_udf`, a Python point lookup (`TokenKeyLookup`) over the release's `token_dict` footers.
   Option (b) of the card; no table was added and the schema lock is unchanged (the choice is below).
2. `open_catalog` attaches the catalog file read-only to a DuckDB instance of the connection's own
   (`duckdb.connect(":memory:")`, then `ATTACH ... (READ_ONLY)` and `USE`).
3. The stored catalog holds `bb_seq`, `bb_seq_circular`, `bb_seq_raw` and `bb_slot_tokens` as stand-in
   macros that raise "<macros>: open this catalog with biobox.catalog.open_catalog(path); <name> is a
   Python function of biobox that a plain duckdb connection does not have". `open_catalog` registers the
   real functions as `<name>_udf` and shadows each stand-in with a temporary macro of the same name.

All numbers below are *measured* on wilkcomp02 on 2026-11-10 in `bio-box-test_74b3e05.sif` (DuckDB 2.0.0),
under `systemd-run --user -p MemoryMax=32G -p MemorySwapMax=0 nice -n 19`, unless labelled otherwise.

## 1. `token_of`: before, the two options, after

### The store

`scripts/stress/synth_scale.py <tree> <out> 1000 --seed 1`, run unchanged by a driver,
`/tmp/if2453/b34/exp/tokq.py <tree> <out> 1000 7680` (units `b34-before-1000` with the tree at `951e3f8`, and
`b34-after-1000` with this branch). The driver wraps `biobox.io.write_table` as
`tests/r2/test_scale_layout.py` does, but only for `token_dict`: its row groups hold at most 7,680 rows
(1/16 of `ROW_GROUP_ROWS`); every other table keeps the real caps. Then it scores each lookup on a fresh
`open_catalog` connection (cold caches) with `catalog.profile_query`. A point lookup's reads are listed as
a `token_dict` scan with `"lookup": "bb_token_of"`.

- `token_dict`: 94,604 rows in 13 row groups, 2,054,252 bytes before and 2,056,150 bytes after. The
  difference, 1,898 bytes (*derived*), is the new `key` statistics in the footer.
- The manifest hash changes, from `052b10ff755847f8a9d69a9d982420dc` to `846855f1232ebe507428501d12f3dfc8`,
  because the `token_dict` file bytes change. It is the only file of the 23 under the two data roots that
  differs (`/tmp/if2453/b34/exp/treediff.py`, sha256 of every file).
- All 29 queries of `synth_scale.py`'s report give the same result before and after.
- Script totals: `total_s` 311.6 s before and 297.4 s after; `peak_rss_mb_queries` 1,093 MB and 1,117 MB.

### Row groups read

"rg" means `token_dict` row groups read, of 13. Results are equal before and after in every row.

| Lookup | Before | After |
|---|---|---|
| `token_of('nr90', k)`, literal key | 13 | 1 |
| `token_of('nr90', key)`, 1 key from a joined table | 13 | 1 |
| `token_of('nr90', key)`, 100 keys from a joined table (keys spread by hash; called twice in the query) | 13 + 13 | 12 (0.12 per key) |
| `token_rows(token_of('nr90', k))` (card B-31's query) | 13 + 13, and postings 1 of 20 | 1, and postings 1 of 20 |
| reverse, token -> key, literal token (`token_dict` by `token_id`) | 1 | 1 |
| reverse, 100 tokens from a joined table (`JOIN token_dict ON token_id`) | 12 | 12 |
| plain SQL `token_dict WHERE level = 'nr90' AND key = ?` (not `token_of`) | 13 | 3 |

The 100 joined keys lie in 12 of the 13 row groups, since the nr90 keys fill 12 of them. So reading 12 is
reading exactly the row groups that hold a key. `tests/r2/test_scale_layout.py` checks this exactly at
100x: 100 keys taken from 2 row groups of 11 read those 2.

The reverse lookup is a join on `token_id`, the sort key. Nothing in this card changes it. It reads the row
groups the 100 tokens fall in, and one for a literal token. This is why `token_dict` keeps its order.

The last row is ad-hoc SQL on the view, not `token_of`. It now prunes on the new `key` statistics, to 3
row groups and not 13. Three is the key's own row group plus the two row groups that straddle a level
boundary, whose `key` ranges span every key. Zone maps cannot do better here; `token_of` can.

### The choice: option (b), no new table

Option (a) was prototyped in the same driver. The same rows were re-sorted by `(level, key)` into one
Parquet file at the same cap, with DuckDB `COPY` and statistics on every column (12 row groups), and queried
with plain SQL:

| Lookup on the (a) prototype | rg read |
|---|---|
| literal key | 2 of 12 |
| 1 key from a joined table | 2 of 12 |
| 100 keys from a joined table, as a join | 11 of 12 |
| 100 keys from a joined table, as `token_of`'s correlated subquery | 11 of 12 |

- A table sorted by `(level, key)` prunes a literal key in SQL, to 2 row groups (one is a level-boundary
  row group).
- It does not prune joined keys. DuckDB pushes only a min/max filter of the 100 build-side keys into the
  scan, and that range spans the table. So joined keys need a point-lookup function with either layout.
- The prototype file is 1,751,900 bytes, against 2,056,150 for `token_dict`. So option (a) adds about 85%
  of `token_dict`'s bytes to every release (*derived*: 1,751,900 / 2,056,150). At `rL.3` (11,258,843
  `token_dict` rows, `2026-11-10_import_ladder/rL3.md`) that is about 0.85 x 11.3 M x 19.5 bytes per row
  = 187 MB (*derived*; 19.5 bytes per row is the measured `rL.2` value in
  `scripts/stress/gtdb_units_2026-11-10.json`).

With a point lookup in either case, `token_dict`'s own order is enough:

- `token_id` is level | source | local, and local is the key's rank within (level, source). So inside each
  (level, source) run, the rows are in key order.
- A row group's `token_id` min and max are its first and last rows, so they name the runs it starts and
  ends in.
- With `key` statistics added, `TokenKeyLookup` reads one row group for a key of a run that fills whole row
  groups. The class docstring states the rule for row groups that hold parts of several runs.
- The new table would cost about 85% more token bytes and buy nothing a footer lookup does not already
  give. So the choice is (b). `tokens.run` writes `key` statistics, and `catalog.py` does the rest.
- Option (a) would also have needed a new registry table and a schema-lock change; (b) needs neither.

### Wall time

Wall time, *measured*, `/tmp/if2453/b34/exp/tokt.py <tree> <catalog>` on the two 1,000x stores above.
It is the median of 7 runs, each on a fresh connection: "cold" is the first execution, "warm" the second
on the same connection. The pod's load average was about 55 of 112 cores during these runs.

| Query | Before, cold | After, cold | Before, warm | After, warm |
|---|---|---|---|---|
| `token_of`, literal key | 0.0113 s | 0.0066 s | 0.0092 s | 0.0016 s |
| `token_of`, 100 joined keys | 0.0179 s | 0.0590 s | 0.0162 s | 0.0159 s |
| `token_rows(token_of(...))` | 0.0256 s | 0.0184 s | 0.0240 s | 0.0120 s |
| reverse, literal token | 0.0060 s | 0.0065 s | 0.0051 s | 0.0054 s |
| reverse, 100 joined tokens | 0.0120 s | 0.0076 s | 0.0114 s | 0.0067 s |

- 100 joined keys, cold, take longer after: 0.018 s to 0.059 s. Before, DuckDB scanned 13 row groups
  natively. After, Python decodes the 12 row groups that hold the keys. Warm, both take 0.016 s.
- The old time grows with the size of `token_dict`, and the new one with the row groups that hold the keys.
  At `rL.3`'s 11.3 M rows (92 row groups at the real cap, *derived*: 11,258,843 / 122,880), the old cost
  is about 92 / 13 = 7 times the scan above. That figure is *derived*, not measured.
- An earlier version wrapped the function in a scalar subquery, `(SELECT bb_token_of_udf(...))`. For the
  100 joined keys it took 0.059 s warm, in 8 calls of 20 to 30 rows. The direct call took 0.016 s, in 2
  calls of 100 rows (*measured*, `/tmp/if2453/b34/exp/tokt2.py`). For a literal or a parameter, both forms
  call the function equally often and prune `token_postings` to 1 of 20 in `token_rows` and `q1_count`.
  So the macro calls the function directly. The `b34-after-1000` run above used the subquery form. With
  the direct form, `tokt2.py` profiles the same `token_postings` reads, and the 100x tests below were run
  with the direct form. The wall-time table is for the direct form. Which row groups `TokenKeyLookup`
  reads does not depend on the form: it gets the same keys.

### 100x test (scaled caps)

`tests/r2/test_scale_layout.py`, with `token_dict` capped at 960 rows, gives 11 row groups at 100x. The
tests failed first on the pod against the tree at `951e3f8` with the new tests copied in. Each read
11 of 11 row groups, per call:

- literal key, family, nr50 and nr90: 11 each;
- 100 joined keys: 22 = 11 + 11, where the oracle says 2 row groups hold the keys;
- 1 joined key: 11.

The existing `test_token_of_key_prunes_postings`, `test_q1_token_of_keys_prune` and
`test_q4_token_of_keys_prune` failed too, at 11 row groups per key. Before this card they passed only
because `token_dict` had 1 or 2 row groups. All of these pass after. The property `token_dict_reads` in the
junit report (`sh scripts/pod_pytest.sh B-34 tests/r2/test_scale_layout.py ... --junitxml=...`) records:

| 100x lookup | rg read after, of 11 |
|---|---|
| literal family key | 1 |
| literal nr50 key | 2 (about 850 nr50 tokens at a 960-row cap, *derived* from the synthetic layout: the level can fill no row group of its own, and the class docstring's rule then reads both row groups it touches) |
| literal nr90 key | 1 |
| 100 joined nr90 keys, held by row groups 1 and 10 | 2 |

## 2. Two connections to one catalog file

*Measured* with `/tmp/if2453/b34/exp/probe.py` (DuckDB 2.0.0):

- A Python function belongs to a DuckDB database instance: it lives in that instance's system catalog. It
  is not scoped to one connection or to the process.
- `duckdb.connect(path, read_only=True)`, called twice in one process, returns two connections to one
  cached instance. `create_function` on the second then fails with `Catalog Error: Scalar Function with name
  "bb_py" already exists!`.
- A second `connect` with other settings, such as another `threads`, fails earlier, with `Connection Error:
  Can't open a connection to same database file with a different configuration than existing connections`.
  This is the failing-first output of the new test, against the old `open_catalog`.
- Two `duckdb.connect(":memory:")` instances can each `ATTACH` the same file `(READ_ONLY)` and each
  register functions with the same name. Closing one leaves the other working.
- A plain `duckdb.connect(path)` can be open beside an attached one.

So `open_catalog` now attaches the file to its own instance, which keeps its `threads`, `memory_limit` and
`temp_directory` per connection.

- `temp_directory` defaults to `<catalog path>.tmp`, as for a file connection. An in-memory instance would
  otherwise spill into `./.tmp`.
- The database keeps the name a plain connection gives it, the file name without `.duckdb`.

`test_b34_two_connections_to_one_catalog_in_one_process` covers this. It opens two connections to one
catalog at 1 and 2 threads, plus one to a second catalog, and checks:

- each answers `seq_get`, `partners_of` and `token_of`;
- after one is closed, the rest give the same answers;
- a fourth connection opened after the close does too.

## 3. A plain connection says what to do

*Measured* with `/tmp/if2453/b34/exp/probe.py` and `probe2.py`:

- A stored macro `bb_x(_a) AS error('...')` raises when it is called directly.
- Called inside another function (`bb_x(_a) || '!'`), the same macro returns NULL with no error. DuckDB's
  binder folds a function whose argument has the NULL type to NULL.
- `CAST(error('...') AS VARCHAR)` raises in both places, and also inside `unnest(...)` of a list-typed
  stand-in. So each stand-in is cast to its function's return type (`VARCHAR`, or `UBIGINT[][]` for
  `bb_slot_tokens`).
- Over zero rows the cast stand-in raises nothing, because nothing is evaluated.
- A Python function cannot be registered under the name of a stored macro: `create_function('bb_x', ...)`
  fails with `bb_x is not an scalar function`. So the functions are registered as `<name>_udf`.
- A temporary macro of the stand-in's name shadows the stored macro, also inside stored macro bodies
  (`seq_get` calls `bb_seq`, `partners_of` calls `bb_slot_tokens`).

`test_b34_plain_connection_error_names_open_catalog` checks that the message names
`biobox.catalog.open_catalog` for `seq_get`, `seq_slice`, `seq_get_circular`, `seq_raw`,
`length(seq_get(...)) + 1` and two `partners_of` queries. It also checks that `token_of` and the views
still work on that connection.

It failed first with `Catalog Error: Scalar Function with name bb_seq does not exist! Did you mean "abs"?`.

## Catalog bytes (AGENTS §4.3)

*Measured*, `/tmp/if2453/b34/exp/catbytes.py <tree> <out>`: one manifest of the micro-S0 fixture, given
twice to `generate_catalog`. The catalog file was never byte-identical across builds, and still is not.

| Tree | File size, both builds | Bytes that differ |
|---|---|---|
| before (`951e3f8`) | 1,060,864 | 67 |
| after | 1,060,864 | 65 |

The differences are in the same places in both trees:

- the file header checksum (bytes 0 to 7) and a 16-byte random database identifier (124 to 139), each also
  in the second header copy (8192 to 8199 and 8316 to 8331);
- the checksums of later blocks;
- a few bytes inside blocks that are uninitialised memory. In one build bytes 537,418 to 537,423 hold
  `a770fd7f0000`, a stack address; in the other they are zeros.

So the catalog cannot be compared by its sha256. `test_b34_catalog_content_is_the_same_on_two_builds`
compares what the file holds instead: views with their SQL, macros with their definitions, `catalog_info`,
`catalog_files` and `catalog_artifacts`. It records the count of differing bytes.

The manifest, which carries the release's record, is unaffected (`catalog.manifest_hash`, tested by
`test_mr01_manifest_hash_ignores_input_order_and_root`). The cross-host check builds the stage-A fixture
store without tokens (`scripts/stress/build_fixture.py`), so this card does not change its manifest hash.

## What is left

- Ad-hoc SQL on `token_dict` by key reads up to 3 row groups at 1,000x with the scaled cap (above). Only
  `token_of` prunes to one.
- Family-level lookups with many small libraries can read more than 2 row groups. Runs that start and end
  inside one row group have keys that no footer bounds. At the real cap every family token of a release
  fits in one or two row groups (`rL.3`: tens of thousands of families, against 122,880 rows per row
  group, *estimate*), so this is not measured.
- `_READERS` and `_LOOKUPS` are keyed by `id(connection)` and never cleared on close; this is unchanged
  from card B-31.
