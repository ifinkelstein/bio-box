# Release `rL.3` finished and synced (wilkcomp02, 2026-11-10)

Stage C of the import ladder (see `README.md`) stopped in its postings step at 02:53 CST. This note
records how `rL.3` was finished, verified and synced, and the three-point scaling of `rL.1`, `rL.2`
and `rL.3`.

Every number is *measured* on wilkcomp02 on 2026-11-10 (CST) unless it is marked *derived*
(formula given). The steps were timed by `scripts/driver.py` (`os.wait4` on each step's process
tree: wall, CPU, peak RSS in MiB, reported as "MB" in the records), run by
`systemd-run --user --unit=ladder-C2 --collect -p MemoryMax=200G -p MemorySwapMax=0 nice -n 19 sh
/tmp/if2453/ladder/scripts/stage.sh C`, from 04:05:07 to 04:10:51 CST. The code was the old
checkout `~/bio-box-next` (`f9af5b6`) inside `bio-box-test_74b3e05.sif`. Raw records:
`results/C/`. The new records are `dim.json`, `catalog.json`, `appendonly.json`, `queries.json`,
`scratch_catalog.json`, `extra_queries.json`, `verify_scratch.json` and the updated `result.json`.
`cluster.json` and `tokens.json` are from the first run and are copied here for the first time.
`postings.json` is the failed first attempt.

## What was done

1. **STOP check.** `/stor/scratch/Wilke/if2453/all-the-data/STOP` was absent at 04:02, and again
   before each unit.
2. **Postings placed.** The B-29 builder had already written `rL.3`'s postings from this store's
   own tapes and tape tokens, to
   `/tmp/if2453/b29/root/l3/token_postings/release=rL.3/tape_version=t-ca596dd85040/token_bucket=all/part-0.parquet`.
   That partition layout equals the old layout (`rL.2` uses
   `release=rL.2/tape_version=t-ca596dd85040/token_bucket=all/part-0.parquet`). The old builder
   would choose the same path:
   - bucket width 0, because `prefix_width(20 × 25,822,709 B)` is about 516 MB, and divided by 16
     that is under the 256 MiB target;
   - one part, because the per-part row cap is 2^30 / 20 = 53.7 M rows, more than 25.8 M (*derived*,
     from `layout.py` and `postings.py` at `f9af5b6`).

   The old builder writes no ledger rows and no `release_def`, so only the parquet file was placed
   (`release_def` was left in the b29 tree). It went in as a hard link under the temporary name
   `.part-0.parquet.tmp`, then `mv -n` into place. Its sha256 after placing is
   `92d69f44695fcb4611da09701918dd030b0ed1ace428201d1ef36f5d82b9c66a` (1,001,038,504 bytes). That
   equals the sha256 that `../2026-11-10_postings_memory.md` records for both the B-29 builder and
   the `f9af5b6` builder with an 8 GB DuckDB cap.
3. **Postings marked done.** On start the driver drops every `result.json` step entry that has an
   `error`, then reruns it. A rerun of the old `postings.run` would now raise `FileExistsError`,
   which is a hard step failure. `/tmp/if2453/ladder/c2/mark_postings.py` (host python3) did three
   things:
   - saved the failed entry under the new top-level key `postings_failed_attempts`;
   - wrote `steps.postings` with `error: null`, the B-29 run's wall 105.4 s, CPU 107.7 s and peak
     2,417.4 MB (from `../2026-11-10_postings_memory.md`), its result (`rows`, `rows_by_level`,
     `buckets` 1, `files` 1, `tape_version`), the file's sha256 and bytes, and an `external` note
     saying the driver did not run this step;
   - removed the top-level `finished`, `failed_hard_step` and `soft_failures`, which the driver
     rewrites.

   The pre-edit copy is `/tmp/if2453/ladder/c2/result.json.before-C2`.
4. **Driver resumed.** It skipped `wait_inputs` through `postings` and ran `dim`, `catalog`,
   `appendonly`, `queries` and `sync`, all without error. The marker is `logs/C.done`, and the
   wrapper exited with rc 0.
5. **Pre-sync dry run.** The driver's own rsync command was run with `-n`, recorded in
   `logs/sync-C.dryrun-pre.log`. It listed 5 new files and no change to any existing file.
   Stage P2's files were already on `/stor`: the `Bp` sync of 00:43 had copied them (README
   finding 1).
6. **After the sync**, three more units ran:
   - `verify_scratch.py` on the `/stor` copy for `rL.1`, `rL.2` and `rL.3` (unit `ladder-C2-verify`);
   - three extra queries and the coverage checks through the `/stor` catalog
     (`/tmp/if2453/ladder/c2/myq.py`, unit `ladder-C2-q2`, with `MemoryMax=32G` and
     `MemorySwapMax=0`);
   - a check of the contigs that have no postings (`nopost.py`).
7. **NOTES.** An addendum was appended to `provisional/rL/NOTES`. Its first 741 bytes are
   unchanged (sha256 `8d2efe8f…` before and after).

The Pfam scan of P2 was not started (about 338 CPU-h, README).

## What is on `/stor` now

`/stor/scratch/Wilke/if2453/all-the-data/provisional/rL`:

- **Data files.** The rsync added 7 data files of `rL.3`: `cluster_member`,
  `cluster_member_aln`, `contig_dim`, `tape_tokens`, `token_dict`, `token_postings`, and
  `ledger/stage=rL.3/step=cluster`. That is 2,847,504,036 bytes in 10.5 s, and no existing file
  was updated (`updated_files: []`).
- **Manifest and catalog.** The scratch catalog build then added
  `manifests/release=rL.3.json` (43,932 B) and `catalog/release=rL.3.duckdb` (1,060,864 B).
- **Rewritten by the driver.** The driver's sync rewrites `README` ("Releases synced so far: rL.1,
  rL.2, rL.2p, rL.3") and `.ladder_releases.json` on every sync, as it did for the earlier stages.
  These are the only existing files whose contents changed apart from the `NOTES` append.
- **Small files.** The run added 4 files under 64 MB: `contig_dim` (1,264,516 B), the `rL.3`
  ledger (2,529 B), the manifest and the catalog. The other 5 are 121 MB to 1.0 GB.

## Verification

| Check | Result |
|---|---|
| `rL.3` manifest on `/stor` (`verify_scratch.py`) | 77 files and 33 artifacts, 13,508,242,096 bytes. `missing_or_changed: []`, 0 files outside the manifest. Its (path, sha256) list equals the local manifest's. Hash `e08feda18e1e13d901472b76c4913cb6` in both. 45.4 s |
| `rL.3` catalog built on `/stor` (`catalog_build.py` in the sync) | `verify_problems: []`. `view_vs_manifest_row_mismatch: {}`: every view's row count equals the manifest's rows |
| `rL.3` rows through the catalog (`queries` step) | `cluster_member` 21,241,368; `token_dict` 11,258,843; `tape_tokens` 32,815; `token_postings` 25,822,709; `contig_dim` 32,815. Each equals the manifest's new rows for that table |
| Append-only (`appendonly` step, local) | 103 files of the `rL.1`, `rL.2` and `rL.2p` manifests rehashed (10,660,738,060 B in 41.1 s). 102 identical; `l0/l0_manifest.parquet` is the expected per-stage rewrite (now 2,498 rows); nothing missing or different |
| `rL.1` and `rL.2` manifests on `/stor` | one problem each: `size differs: l0/l0_manifest.parquet`. This is the per-stage rewrite of README finding 5. The `B` sync already overwrote that file on `/stor` (`updated_files` in `results/B/result.json`), and the `rL.3` sync changed it no further. Every other file verifies. The (path, sha256) lists equal the local manifests' |
| Stages in the `rL.3` manifest | paths name `stage=` S1 (19 files), P1 (19) and P2 (42), plus the ledger stages `rL.1`, `rL.2`, `rL.2p` and `rL.3` (1 file each). No other stage is on disk: `l0seq/contig_index` holds only `stage=S1`, `P1` and `P2` |
| Default-variant contigs per stage (`contig_index`, `/stor` catalog) | S1 27,258; P1 539; P2 5,018; sum 32,815 |
| `contig_dim` of `rL.3` | 32,815 rows and 32,815 distinct contigs. 0 contigs not in `contig_index`, 0 index contigs missing from it |
| Postings contigs | 32,805 of 32,815. The 10 without postings are P2 contigs with 0 slots (no gene called; 423 to 33,111 bp). `contig_tapes` has exactly 10 zero-slot tapes |
| Slots | 12,048,866 = 2,398,087 + 927,920 + 8,722,859 (`contig_tapes` by stage) |

The `rL.3` manifest also lists the release-partitioned files of `rL.1`, `rL.2` and `rL.2p` (6
each). Its catalog views therefore hold every release on disk, `rL.2p` included. A query must
filter on `release` or use `contig_dim_now()`. That is the old layout, and the earlier releases
work the same way.

## What `rL.3` contains and lacks

- Clusters over 10,620,684 proteins of S1, P1 and P2: 2,966,200 NR50 and 8,278,537 NR90 clusters
  (`cluster.json`).
- Tokens: 14,106 family, 2,966,200 NR50 and 8,278,537 NR90 (`tokens.json`). Every slot carries an
  NR50 and an NR90 token.
- **No Pfam for P2.** `l1_domains` has S1 (605,564 rows) and P1 (1,552,451 rows) only. Family
  postings by stage:

  | Stage | Family postings | Share of the stage's slots (*derived*) |
  |---|---|---|
  | S1 | 941,262 | 39.3% |
  | P1 | 757,230 | 81.6% |
  | P2 | 26,485 | 0.30% |

  P2 slots get a family token only where their protein is identical to a scanned S1 or P1 protein.
  Family-level questions on prokaryotes are not meaningful in `rL.3`. The representative of the
  largest prokaryote NR50 cluster (query b) has no Pfam row.
- Postings were written by the B-29 builder, not by the `f9af5b6` builder that ran the other
  steps. The bytes are identical, but the timing comes from a different builder.

## Steps run in this session

| Step | Wall s | CPU s | Peak RSS MB | Notes |
|---|---|---|---|---|
| dim | 4.0 | 4.3 | 326.4 | 32,815 contigs, 29,755 genomes |
| catalog (local) | 142.0 | 121.0 | 2,057.8 | build_manifest 78.44 s, generate_catalog 10.52 s, verify_manifest 49.63 s |
| appendonly | 44.6 | 44.1 | 312.1 | 10.66 GB hashed |
| queries | 16.1 | 34.9 | 1,597.9 | DuckDB 4 threads |
| sync | 137.6 | | | rsync 10.5 s, 14.0 MB peak. Scratch catalog 127.1 s, 2,127.5 MB peak (build_manifest 71.39 s, verify 47.47 s, over NFS) |
| extra queries (`myq.py`) | 17.17 | | 2,579 | `/usr/bin/time -v`, peak includes the singularity wrapper |
| verify_scratch (rL.1, rL.2, rL.3) | 60.6 (sum of `verify_s`) | | | not timed by wait4 |

## Bytes (from the manifests)

| Table, new files of the release | `rL.1` | `rL.2` | `rL.3` |
|---|---|---|---|
| `cluster_member` rows / bytes | 2,269,100 / 70,337,946 | 4,106,150 / 148,588,878 | 21,241,368 / 762,906,723 |
| `cluster_member`, bytes per row (*derived*) | 31.0 | 36.2 | 35.9 |
| `cluster_member_aln`, bytes per row (*derived*) | 31.0 | 34.2 | 35.0 |
| `token_dict` rows / bytes | 836,585 / 16,605,981 | 2,162,947 / 42,277,387 | 11,258,843 / 216,384,222 |
| `token_dict`, bytes per token (*derived*) | 19.8 | 19.5 | 19.2 |
| `tape_tokens` bytes, bytes per slot (*derived*) | 18,329,399, 7.64 | 30,523,869, 9.18 | 121,417,336, 10.08 |
| `token_postings` rows / bytes | 5,737,436 / 198,642,343 | 7,593,431 / 277,833,790 | 25,822,709 / 1,001,038,504 |
| `token_postings`, bytes per row (*derived*) | 34.6 | 36.6 | 38.8 |
| `token_postings`, rows per slot (*derived*) | 2.39 | 2.28 | 2.14 |
| `contig_dim`, bytes per contig (*derived*) | 39.6 | 39.7 | 38.5 |
| whole manifest: files + artifacts, bytes | 25 + 3, 1,623,256,748 | 46 + 6, 2,996,804,176 | 77 + 33, 13,508,242,096 |
| release-only bytes added (*derived*, sum of new files) | | | 2,847,504,036 (= rsync transfer) |

Postings bytes per row keep rising (34.6, 36.6, 38.8) because rows per slot fall: the family share
drops as unscanned P2 slots dominate. `EST_POSTING_ROW` is 20 in the old code and 36 in the B-29 builder
(`../2026-11-10_postings_memory.md`).

## Queries

The `queries` step runs locally on `/tmp/if2453/ladder/root` and its catalog. The `rL.1` and
`rL.2` columns come from `results/A` and `results/B`, and the `rL.2p` column from the README.
"rg" is row groups scanned / row groups the scan could read (`catalog.profile_query`).

| Query (`queries` step) | `rL.1` | `rL.2` | `rL.2p` | `rL.3` |
|---|---|---|---|---|
| open catalog | 0.38 s | 0.40 s | 0.39 s | 0.36 s |
| 1,000 random genes ±500 bp | 0.589 s | 0.656 s | 0.746 s | 0.787 s |
| the same, warm | 0.477 s | 0.532 s | | 0.613 s |
| index row groups read for it (locator, index) | 1, 1 | 2, 2 | 3 | 3, 3 |
| tapes of one contig, rg | 1/19 | 1/26 | 1/93 | 1/93 |
| NR50 postings lookup (most abundant NR50), rg | 1/47 | 1/62 | 1/69 | 1/211 |
| NR50 postings lookup, wall | 0.045 s | 0.034 s | | 0.052 s |
| terminase family PF03237.22 → family partners | 0.089 s, rg 1/47 + 8/47 | 0.194 s, 1/62 + 8/62 | 0.21 s | 0.225 s, 1/211 + 15/211 |
| terminase family → NR50 partners | 0.124 s, 1/47 + 21/47 | 0.293 s, 1/62 + 28/62 | | 0.412 s, 1/211 + 99/211 |
| most abundant NR50 → NR50 partners | 0.144 s, 1/47 + 21/47 | 0.267 s, 1/62 + 28/62 | | 1.312 s, 1/211 + 99/211 |
| most abundant prokaryote NR50 → NR50 partners | | 0.301 s, 1/62 + 28/62 (345 anchors) | | 1.210 s, 1/211 + 99/211 (3,509 anchors) |

The terminase family has 6,200 anchors in `rL.1` and `rL.2` and 6,261 in `rL.3`. The most
abundant NR50 is a different token in `rL.3`: `3a04ab75…`, 3,509 slots, an EF-Tu representative,
against `346053ce…` with 2,139 slots before.

Three extra queries ran through the **`/stor`** catalog (`myq.py`; `results/C/extra_queries.json`).
The `.duckdb` was copied to local disk, and the store was read over NFS.

| Query (`rL.3`, `/stor`) | Wall | Read |
|---|---|---|
| open catalog | 0.691 s | |
| (a) 1,000 random genes ±500 bp (`ORDER BY hash(feature_uid, 'rL.3-c2')`; 813 prokaryote, 187 phage, 881 contigs) | 0.904 s; warm 0.670 s; picking 0.51 s | locator and index row groups 3 and 3 (0 and 0 warm). 1,033 BGZF blocks, 20,402,604 compressed bytes. 1,856,765 bp returned = the expected bp |
| (b) NR50 neighbours (±3 slots) of the largest prokaryote NR50 cluster (`f66315ba…`: 2,899 distinct member proteins on prokaryote contigs, 2,908 in all; 2,950 slots on 2,331 contigs) | 1.433 s; top-10 partners 1.958 s; finding the cluster 3.68 s | anchors rg 1/211. Partners rg 99/211 (10,052,622 rows scanned, 86.7 MB read). 17,678 pair rows, every anchor has a partner |
| (c) `token_postings` ⋈ `contig_dim_now()`, count by level and `assembly_quality` | 1.466 s | postings rg 211/211 (a full aggregate), `contig_dim` 1/1. Total 25,822,709 = manifest. Family 1,724,977 on 28,965 contigs; NR50 and NR90 12,048,866 each on 32,805 contigs; every row has `assembly_quality` "complete" |

**Partner lookups do not prune by contig.** In every neighbour query the anchor side reads 1 row
group, but the partner side reads the whole token-level range of the release:

| Partner level | Row groups read (`rL.1`, `rL.2`, `rL.3`) | `rL.3` fraction of row groups (*derived*) | Level's share of `rL.3` postings (*derived*) |
|---|---|---|---|
| NR50 | 21/47, 28/62, 99/211 | 99/211 = 0.469 | 12,048,866 / 25,822,709 = 0.467 |
| family | 8/47, 8/62, 15/211 | 15/211 = 0.071 | 1,724,977 / 25,822,709 = 0.067 |

The filter `contig_uid IN (SELECT … FROM a)` is a joined key on a file sorted by `token_id`. So
this read grows linearly with the release, against AGENTS §4.4. The README's "lookups do constant
work" holds for the anchor and tape lookups, not for partners.

## Scaling, three points

The size columns give each release's size as a multiple of `rL.1`'s. The "exponent" columns are
*derived* as ln(t₂/t₁) / ln(n₂/n₁) between consecutive points, with n the unit in the second
column.

| Step | Unit n (`rL.1` / `rL.2` / `rL.3`) | Size ×`rL.1` | Wall s | Peak RSS MB | Exponent `rL.1`→`rL.2` | Exponent `rL.2`→`rL.3` | Rate (*derived*) |
|---|---|---|---|---|---|---|---|
| clustering, 16 threads | proteins 1,134,550 / 2,053,075 / 10,620,684 | 1 / 1.81 / 9.36 | 473.5 / 1,319.9 / 9,463.5 | 5,351.6 / 9,570.9 / 46,073.1 | 1.73 | 1.20 | 4,946 / 4,888 / 4,549 B of peak RSS per protein. The memory exponent `rL.2`→`rL.3` is 0.96 |
| tokens | slots 2,398,087 / 3,326,007 / 12,048,866 | 1 / 1.39 / 5.02 | 21.2 / 38.1 / 219.7 | 1,299.2 / 1,762.5 / 2,336.6 | 1.79 | 1.36 | 113 / 87 / 55 k slots/s. Per `cluster_member` row of the release: 107 / 108 / 97 k rows/s |
| postings | rows 5,737,436 / 7,593,431 / 25,822,709 | 1 / 1.32 / 4.50 | 20.0 / 27.9 / 105.4 (`rL.3`: B-29 builder) | 1,353.3 / 1,748.5 / 2,417.4 | 1.19 | 1.09 | 287 / 272 / 245 k rows/s |
| catalog (local) | manifest bytes 1.62 / 3.00 / 13.51 GB | 1 / 1.85 / 8.32 | 17.6 / 33.5 / 142.0 | 421.0 / 745.7 / 2,057.8 | 1.05 | 0.96 | 10.8 / 11.2 / 10.5 s per GB. build_manifest 255 / 231 / 172 MB/s; verify 278 / 250 / 272 MB/s |
| catalog on `/stor` (in sync) | same | | 18.1 / 30.3 / 127.1 | 346.1 / 661.2 / 2,127.5 | | | |
| dim | contigs 27,258 / 27,797 / 32,815 | | 3.9 / 4.0 / 4.0 | 326.4 / 325.1 / 326.4 | | | |
| queries step, whole | | | 6.1 / 7.8 / 16.1 | 505.0 / 778.0 / 1,597.9 | | | |
| 1,000 genes ±500 bp | stages 1 / 2 / 3 | | 0.589 / 0.656 / 0.787 | | | | one locator and one index row group per stage |
| NR50 postings lookup | postings row groups 47 / 62 / 211 | | 0.045 / 0.034 / 0.052 | | | | 1 row group at every size |
| neighbours, terminase → family | | | 0.089 / 0.194 / 0.225 | | | | partner rg 8 / 8 / 15 |
| neighbours, terminase → NR50 | | | 0.124 / 0.293 / 0.412 | | | | partner rg 21 / 28 / 99 |
| neighbours, top NR50 → NR50 | | | 0.144 / 0.267 / 1.312 (other anchor token in `rL.3`) | | | | partner rg 21 / 28 / 99 |

Readings:

- **Clustering** is about n^1.2 in time and linear in memory, at about 4.5 to 4.9 KB per protein,
  between 2.05 M and 10.6 M proteins. This confirms the README's stage C outcome.
- **Tokens** slow per slot (113 to 55 k slots/s). Their time tracks the release's
  `cluster_member` rows, which grow faster than slots (21.2 M members for 12.0 M slots in `rL.3`),
  at about 100 k rows/s.
- **Postings** are near-linear, but the `rL.3` point is the B-29 builder.
- **The catalog** is linear in manifest bytes, at about 11 s per GB. Its peak RSS grew 4.9 times
  for 8.3 times the bytes, which was not investigated. In `build_manifest`, hashing fell to
  172 MB/s on `rL.3` from 231 to 255 MB/s earlier. That is 31% below the README's "about 250 MB/s
  per hashing pass". The cause was not investigated.
- **Point lookups** stay at 1 row group, and the sequence fetch adds one index row group per stage.
  Partner scans grow with the release (above).

## Open points

- Partner-side pruning of neighbour queries over `token_postings` (above). I know of no card that
  covers it.
- `rL.3`'s family tokens cover only 0.30% of P2 slots until P2 is scanned with Pfam (about
  338 CPU-h, README).
- On `/stor`, the `rL.1` and `rL.2` manifests fail verification on `l0/l0_manifest.parquet` alone
  (README finding 5, card B-22).
