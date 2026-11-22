# Card B-62: Bounded memory and open files in `biobox.l3.postings`

Date: 2026-11-21
Host: `wilkcomp02` (inside `/stor/home/if2453/bio-box/images/bio-box-test_faec1f4.sif`, `systemd-run --user --scope -p MemoryMax=16G -p MemorySwapMax=0 nice -n 19`, data root `/tmp/if2453/r10/root` bound `:ro`).

## 1. Root causes

1. **Fault 1 (`_family_run` OOM at `memory_limit="2GB"`)**:
   - `_family_run` ran a single DuckDB statement over every `l1_domains` row across all release stages (`SELECT DISTINCT`, `row_number() OVER (PARTITION BY protein_uid ...)`, and `s JOIN g ON g.protein_uid = s.pu`). Because DuckDB's hash join build side does not spill cleanly once `g` exceeds the memory limit (card B-46), `r26.11.2` failed with `_duckdb.OutOfMemoryException` at `1.8 GiB / 1.8 GiB used` after `229 s` / `3.77 CPU-min` (*measured*, wilkcomp02, 2026-11-21, `c46-train-rel-2.log`).
   - Why `tests/r2/test_micro_s0_postings.py::test_b29_long_contigs_stay_inside_the_memory_limit` did not catch it: its synthetic fixture (`_long_contigs`) puts a family domain on only 3 CDSs in 1,000 (~1,967 domain rows out of 690,000 slots, *derived*) at `memory_limit="115MB"`, so `g` stayed tiny (~200 KB) even without range partitioning.

2. **Fault 2 (`Too many open files` on run/canon files)**:
   - `_explode_runs` wrote one `run-XXXXXX.parquet` and one `canon-XXXXXX.parquet` per `tape_tokens` row group (at most 131,072 slots per row group, AGENTS §4.6). On `r26.11.2`, that produced >545 `canon-*.parquet` and >545 `run-*.parquet` files, and `_check_canonical` / the final sort passed all of them to a single `read_parquet([...])`. Under a `systemd-run --user` unit's default soft limit `RLIMIT_NOFILE = 1024` (*measured*, wilkcomp02, 2026-11-21, `probe8.log`), `_check_canonical` failed on `canon-000544.parquet` with `IOException: Too many open files`.

## 2. Fix

- **`_family_run` range partitioning (`src/biobox/l3/postings.py`)**:
  - When `k = tokens.n_ranges(max(n_dom, n_slots), memory_limit) > 1`, domain rows and unnested non-null protein slots are routed by `hash(protein_uid) % k` into `k` independent scratch Parquet pairs via `_route` (at most `_ROUTE_OPEN_FILES = 64` open writers and `_ROUTE_HELD_ROWS = 262_144` buffered rows in Python).
  - Each range `r` runs the exact `DISTINCT` + `row_number()` window + protein join independently (since both the window and the join are strictly per `protein_uid`), streams into the single `family.parquet` writer (preserving the missing-token `ValueError` check identically), and unlinks its range files immediately.
- **Coalesced runs and bounded `read_parquet` (`src/biobox/l3/postings.py`)**:
  - `_explode_runs` coalesces exploded cluster batches and canonical batches into multi-row-group Parquet files of up to `_RUN_MAX_ROWS = 50_000_000` rows (~400 row groups per file), cutting `r26.11.2` from >545 files each to 3 cluster run files and 2 canon files.
  - `_coalesce_parquet` guarantees that no `read_parquet([...])` call in `_check_canonical` or the final sort ever receives more than `_MAX_READ_FILES = 64` files (streaming with at most 2 open descriptors at a time if a future corpus exceeds `_MAX_READ_FILES * _RUN_MAX_ROWS`).

## 3. Small reproductions (test first)

- **Fault 1** (`tests/r2/test_b62_postings_bounds.py::test_b62_family_run_stays_inside_memory_limit`, 400,000 slots, 1,600,000 domain rows):
  - Before fix (and with unsplit `k = 1`): failed with `_duckdb.OutOfMemoryException` in `_family_run` at `45MB` (`42.9 MiB/42.9 MiB used`), `50MB` (`47.7 MiB/47.6 MiB used`), and `55MB` (`51.6 MiB/52.4 MiB used`) (*measured*, wilkcomp02, 2026-11-21, `pod_pytest.sh`).
  - After fix: passes at `55MB` and produces postings byte-identical (`sha256`) to a `2GB` run.
- **Fault 2** (`tests/r2/test_b62_postings_bounds.py::test_b62_many_row_groups_stay_under_nofile_limit`, 180 `tape_tokens` row groups under soft `RLIMIT_NOFILE = 128`):
  - Before fix: failed in `_check_canonical` with `_duckdb.IOException: Cannot open file ".../runs/canon-000016.parquet": Too many open files` (*measured*, wilkcomp02, 2026-11-21, `pod_pytest.sh`).
  - After fix: passes and produces postings byte-identical (`sha256`) to an unconstrained run; also passes when `_RUN_MAX_ROWS = 1` and `_MAX_READ_FILES = 16` force `_coalesce_parquet` under `RLIMIT_NOFILE = 128`.

## 4. Real-release verification on `wilkcomp02`

Both runs used `memory_limit="2GB"`, soft `RLIMIT_NOFILE = 1024`, data root `/tmp/if2453/r10/root` bound `:ro`, and scratch under `/tmp/if2453/b62/` (removed after recording).

| Release | Stages | Status | Byte-identical vs published | Rows (total) | Family | NR50 | NR90 | Buckets / Files | Output bytes | Wall (s) | CPU (s / min) | Peak RSS (KB / GiB) |
|---|---|---|---|---:|---:|---:|---:|---|---:|---:|---:|---:|
| `r26.11.1` | `S1,P1..P6,H1,PL1` | ok (*measured*) | **yes** (16/16 `sha256` match) | 201,822,651 | 87,648,879 | 57,086,886 | 57,086,886 | 16 / 16 | — | 2,974.1 | 1,685.6 / 28.09 | 2,814,696 / 2.68 GiB |
| `r26.11.2` | `S1,P1..P6,H1,PL1,C5,V1` | ok (*measured*) | n/a (unpublished) | 252,132,516 | 106,276,604 | 72,927,956 | 72,927,956 | 16 / 16 | 10,093,179,143 | 3,679.1 | 2,165.1 / 36.09 | 2,938,376 / 2.80 GiB |

Notes:
- `r26.11.1` previously measured `27.5 CPU-min` and `2.72 GB` (`2.53 GiB`) peak RSS (*measured*, `TASKS.md` card B-62); the ranged build took `28.09 CPU-min` (+2.1%, *derived*) and `2.68 GiB` (`2.88 GB`) peak RSS, with exact byte-for-byte `sha256` equality on all 16 files.
- `r26.11.2` peak RSS is `2.80 GiB` (`3.01 GB`, i.e., `2 GB` DuckDB limit plus `~1.0 GB` Python/Arrow/OS margin), well below `probe8.log`'s `8.41 GB` at `8GB`, and `36.09 CPU-min` is well under the 280/300 CPU-minute process cap.
