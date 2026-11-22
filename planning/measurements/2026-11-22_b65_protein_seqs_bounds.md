# Card B-65: Bounded Row-Group-Pruning `export.protein_seqs` and `hits_faa` at Multi-Stage Scale

**Date**: 2026-11-22
**Host**: `wilkcomp02` (read-only verification on `/tmp/if2453/r10/root`, release `r26.11.2`; outputs written to scratch `/tmp/if2453/b65` and cleaned afterward; zero writes to root or `/stor`)
**Image**: `bio-box_faec1f4.sif` (`sha256:d67300aa1321b571fa8da79e74dfccdd87a8dca2e56518f324902efa6b9f098a`)
**Scope**: Card `B-65` (depends on `W-23`, `B-27`)

## 1. Root Cause of `r26.11.2` `seq-oracle` OOM

On release `r26.11.2` (`11` stages, `16` `protein_seq` Parquet files, `494` row groups, `59,418,086` `protein_seq` rows, *measured* `wilkcomp02` 2026-11-22), both worked examples (`crispr_cas_neighbourhoods`, `4,264` unique uids in `4,265` lines of `oracle_uids.txt`; `rnap`, `2,751` unique uids in `2,752` lines of `oracle_uids.txt`) failed in `bio-box seq` (`src/biobox/search/export.py:996`) with `_duckdb.OutOfMemoryException: Out of Memory Error: failed to allocate data of size 2.0 MiB (1.8 GiB/1.8 GiB used)` and cgroup `MemoryPeak = 2.36 GB` (*measured*, `release_examples.sh r26.11.2`, 2026-11-22):

1. `export.protein_seqs` opened DuckDB over `prot` (a `UNION ALL` of 11 per-stage `read_parquet` scans) and capped `dynamic_or_filter_threshold` at the row-group count (`494`). Above that threshold, DuckDB dropped the key-list filter and pushed only `[min(uid), max(uid)]`, scanning all `59.4M` `seq` strings (~18 GB uncompressed, *estimate*) across all 11 branches concurrently.
2. `open_store` left DuckDB's `preserve_insertion_order = True` enabled, buffering unordered Parquet chunks in RAM.
3. `examples/crispr_cas_neighbourhoods/run.sh` (line 128) and `examples/rnap/run.sh` (line 85) did not pass `--threads "$THREADS"` to `bb seq-oracle seq`, leaving DuckDB at `EXPORT_THREADS = 4` inside a `2`-CPU unit.
4. `export.hits_faa` used the same `SELECT protein_uid, min(seq) AS seq FROM prot WHERE protein_uid IN (SELECT protein_uid FROM _p) GROUP BY protein_uid` pattern.

## 2. Changes Implemented

- `src/biobox/search/export.py`:
  - `open_store.__enter__`: set `"preserve_insertion_order": False` in `duckdb.connect(config={...})`.
  - `protein_seqs`: validate `files = release_files(root, release, manifest)` first, resolve explicit `protein_seq` files in deterministic sorted path order, read only Parquet footers (`pq.read_metadata`), enforce `LY-02`/`LY-03` sort order and min/max statistics on `protein_uid`, bisect `wanted = sorted(set(uids))` against each row group's `[min, max]`, open at most one `pq.ParquetFile(..., pre_buffer=False)` at a time, check `protein_uid` hits first (`use_threads=False`, avoiding PyArrow's `112`-thread host pool on `wilkcomp02`, *measured*), stream matching row groups through `pf.iter_batches(batch_size=_BATCH_ROWS, ...)` (`_BATCH_ROWS = 8,192`) so only one batch of `seq` strings is decoded at a time, aggregate `(max(length), min(seq))` across stages, and call `pa.default_memory_pool().release_unused()` after each file.
  - `hits_faa`: removed `CREATE OR REPLACE TEMP TABLE _seq` from DuckDB; streams metadata rows from `_h`/`_genes`/`_plen` in chunks of `_BATCH_ROWS` (`8,192`) and resolves each chunk's sequences via `protein_seqs`.
- `examples/crispr_cas_neighbourhoods/run.sh` (line 128) and `examples/rnap/run.sh` (line 85): added `--threads "$THREADS"` to `bb seq-oracle seq`.

## 3. Live `r26.11.2` `seq-oracle` Verification (`MemoryMax=1G`, `--threads 2`)

Commands run on `wilkcomp02` (2026-11-22) via `scripts/pod_run.sh --mem 1G --threads 2 --bind /tmp/if2453/r10/root:/tmp/if2453/r10/root:ro --bind /tmp/if2453/c46/train:/tmp/if2453/c46/train:ro`:

```bash
bio-box seq --root /tmp/if2453/r10/root --release r26.11.2 --threads 2 \
  --proteins-file /tmp/if2453/c46/train/examples/crispr_cas_neighbourhoods/oracle_uids.txt \
  --out /tmp/if2453/b65/crispr_oracle.faa

bio-box seq --root /tmp/if2453/r10/root --release r26.11.2 --threads 2 \
  --proteins-file /tmp/if2453/c46/train/examples/rnap/oracle_uids.txt \
  --out /tmp/if2453/b65/rnap_oracle.faa
```

| Example | Input uids (*measured*) | Output FASTA records (*measured*) | Output bytes (*measured*) | Exit (*measured*) | Wall time (*measured*) | CPU time (*measured*) | Cgroup `MemoryPeak` (*measured*) | Process `ru_maxrss` (*measured*) | Reduction vs 2.36 GB OOM (*derived*) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `crispr_cas_neighbourhoods` | 4,264 (`4,265` lines incl. `#protein_uid`) | 4,264 | 2,116,470 | 0 | 44.30 s | 41.23 s (`41233873000` ns) | `455,041,024` B (`433.96` MiB) | `318.2` MiB | `-81.6%` cgroup peak (`433.96` MiB vs `2,360` MiB) |
| `rnap` | 2,751 (`2,752` lines incl. `#protein_uid`) | 2,751 | 2,304,732 | 0 | 21.82 s | 18.47 s (`18470572000` ns) | `448,155,648` B (`427.39` MiB) | `318.2` MiB | `-81.9%` cgroup peak (`427.39` MiB vs `2,360` MiB) |

Note on `use_threads=False` + `iter_batches`: on `wilkcomp02`, PyArrow's default thread pool reports `pa.cpu_count() = 112` (*measured*, 2026-11-22). An initial prototype calling `pf.read_row_group(..., columns=["protein_uid", "length", "seq"])` with default `use_threads=True` exited `0` under `MemoryMax=1G` at `MemoryPeak = 648,110,080` B (`618.1` MiB, *measured*) and `wall = 57.39 s`. Switching to `pre_buffer=False`, `use_threads=False`, a `protein_uid`-first presence check, `pf.iter_batches(batch_size=8192, ...)`, and per-file `pa.default_memory_pool().release_unused()` cut cgroup `MemoryPeak` to `433.96` MiB (`crispr`) and `427.39` MiB (`rnap`), process `ru_maxrss` to `318.2` MiB, and wall time to `44.30 s` and `21.82 s` — both well below the `512` MiB acceptance bound.
