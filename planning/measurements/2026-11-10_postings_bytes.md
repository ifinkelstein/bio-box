# `token_postings` bytes per column (card B-26)

Card B-26. The import ladder measured `token_postings` at 34.6 to 36.8 bytes per row; `postings.py`
stated 20. This note says which columns carry the bytes.

## Method

- Host: wilkcomp02 (`if2453@wilkcomp02.ccbb.utexas.edu`), date 2026-11-10, one thread, image
  `bio-box-test_74b3e05.sif`.
- Store: `/stor/scratch/Wilke/if2453/all-the-data/provisional/rL`, bound read only; release `rL.2`
  (stages S1 + P1). File list from `manifests/release=rL.2.json` (2 files, 109 row groups).
- Script: `scripts/stress/postings_bytes.py`. It sums `ColumnChunkMetaData.total_compressed_size` and
  `total_uncompressed_size` over every row group of every file. Rows per level: `token_id >> 61`
  from the row groups' min/max statistics; 4 of 109 row groups spanned two levels, and only their
  `token_id` column was read. No other data page was read.
- Command (prefix as built by `scripts/pod_pytest.sh`, with the store bound read only):

```
systemd-run --user --scope --quiet -p MemoryMax=32G -p MemorySwapMax=0 nice -n 19 \
  singularity exec --cleanenv --env PYTHONPATH=/tmp/if2453/ci/B-26/src:/tmp/if2453/pilot_s1/pydeps \
  --bind /tmp/if2453 \
  --bind /stor/scratch/Wilke/if2453/all-the-data/provisional/rL:/stor/scratch/Wilke/if2453/all-the-data/provisional/rL:ro \
  --pwd /tmp/if2453/ci/B-26 /stor/home/if2453/bio-box/images/bio-box-test_74b3e05.sif \
  python scripts/stress/postings_bytes.py /stor/scratch/Wilke/if2453/all-the-data/provisional/rL rL.2
```

## Result (*measured*)

13,330,867 rows; 476,476,133 file bytes (35.74 per row, footer overhead included).

Rows per token level: nr90 5,724,094; nr50 5,724,094; family 1,882,679. That is 2.329 rows per CDS
slot (13,330,867 / 5,724,094), with a family token on 0.329 of the slots.

| Column | Compressed B | Uncompressed B | Compressed B/row | Share | Encodings (codec zstd) |
|---|---|---|---|---|---|
| `token_id` | 12,479,401 | 34,187,272 | 0.936 | 2.6% | dictionary + RLE, PLAIN |
| `contig_uid` | 63,499,118 | 63,622,885 | 4.763 | 13.3% | dictionary + RLE, PLAIN |
| `ordinal` | 18,296,957 | 20,190,337 | 1.373 | 3.8% | dictionary + RLE, PLAIN |
| `seg_idx` | 41,436 | 33,624 | 0.003 | 0.0% | dictionary + RLE, PLAIN |
| `feature_uid` | 226,995,112 | 227,553,530 | 17.028 | 47.6% | dictionary + RLE, PLAIN |
| `start` | 54,318,291 | 97,183,172 | 4.075 | 11.4% | dictionary + RLE, PLAIN |
| `end` | 54,397,081 | 97,453,892 | 4.081 | 11.4% | dictionary + RLE, PLAIN |
| `strand` | 1,654,572 | 1,735,208 | 0.124 | 0.3% | dictionary + RLE, PLAIN |
| `n_slots` | 15,783,215 | 16,101,905 | 1.184 | 3.3% | dictionary + RLE, PLAIN |
| `is_circular` | 31,011 | 1,682,949 | 0.002 | 0.0% | RLE, PLAIN (no dictionary) |
| `contig_len` | 28,888,373 | 39,396,872 | 2.167 | 6.1% | dictionary + RLE, PLAIN |
| total | 476,384,567 | 599,141,646 | 35.735 | 100% | |

The writer (`io.write_table`) uses pyarrow defaults with `compression="zstd"`: dictionary encoding
on, no delta encodings. Every column but `is_circular` has a dictionary page in all 109 row groups,
and the PLAIN entry is the fallback once a dictionary page passes pyarrow's 1 MiB default limit
(*inferred* from both encodings being listed; page-level encoding stats were not read).

## Where the excess is

- `feature_uid`: 17.03 B per row, 47.6% of the file. It is 16 random bytes; zstd takes 227.6 MB
  to 227.0 MB. This alone is 85% of the old estimate of 20 B.
- `contig_uid`: 4.76 B per row (13.3%), 99.8% of its uncompressed size. The sort key puts `token_id`
  before `contig_uid`, so within a row group each contig appears only a few times
  (4.76 / 16 = 0.30 dictionary entries of 16 B per row, *derived*, so about 3.4 rows per distinct
  contig per row group), and the 16-byte dictionary entry is paid again in every row group.
- `start`, `end`: 4.08 B each (22.8% together). They are 8-byte integers (7.3 B per row before zstd);
  zstd halves them. Neighbouring rows are different tokens on different contigs, so the values are
  unrelated and neither delta encoding nor the sort key helps them.
- `contig_len` 2.17 B, `n_slots` 1.18 B, `ordinal` 1.37 B, `token_id` 0.94 B: the sort key does help
  `token_id` (runs; 34.2 MB to 12.5 MB) and `strand`, `is_circular` and `seg_idx` are near zero.
- The 20 B estimate (PLAN.md §4.2, 15 to 25 B) covers `feature_uid` plus 3 B. The row also carries
  two UIDs and five coordinate-like integers.

The three excess groups together: the two UIDs and the two coordinates are 30.0 B of 35.7 B (84%).

## Consequences

- `EST_POSTING_ROW` in `src/biobox/l3/postings.py` is 36 (35.74, rounded up to an integer).
  It cuts parts at `max_file_bytes` and picks the bucket width through `layout.prefix_width`.
- Projection at GTDB scale (*derived*): 437 M proteins (`planning/sources_and_sizing.md`, *measured*
  from GTDB R226 metadata, 2026-11-09) × 2.329 rows per slot × 35.74 B = 1.018e9 rows and about
  36.4 GB. The 2.329 comes from S1 + P1 (phages and a few prokaryote sets); the family share
  (0.329 of slots) will differ with Pfam coverage on GTDB, so this is an *estimate* in that factor.
  With the old 20 B the same formula gives 20.4 GB.

## Open questions (a design decision, not made here)

Dropping the columns that the tape already holds under `feature_uid` and `contig_uid` (`start`,
`end`, `strand`, `n_slots`, `is_circular`, `contig_len`) would remove 4.075 + 4.081 + 0.124 + 1.184
+ 2.167 + 0.002 = 11.633 B of 35.735 B per row (*derived* from the table above), a 32.6% saving,
above the 20% line of AGENTS.md §3 and §6. The cost is a join to the tapes at query time for
anything that needs coordinates, and a schema change. Replacing `contig_uid` by a release-local
integer would save at most about 4.7 B (13%), but AGENTS.md §3.4 limits release-local ids to a
`release=` partition and this is a design question. Neither was made; both belong to an **H** card.
