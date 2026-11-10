# Neighbour and lookup macros that prune (card B-31), 2026-11-10

Card B-31 lists three defects in `src/biobox/catalog.py`, all found by B-18 and D-05:

1. a token key given through `token_of(...)` (a scalar subquery) was not pushed into the `token_postings` scan;
2. the partner side of `partners_of` filtered `token_postings` by `contig_uid`, which cannot prune a table
   sorted by `token_id`;
3. `default_digest`, `digest_len` and `contig_len_of` read every index row group for a joined key.

What changed:

1. `token_rows` joins the token in (`(SELECT _tok) _k JOIN token_postings ON token_id = _k._t`), so DuckDB
   pushes the one build-side value into the postings scan as a join filter. With a literal, the bucket
   filter still prunes files. `token_bucket_ok` now names its token once, so a `token_of` key adds
   one `token_dict` scan instead of three.
2. `partners_of` still takes its anchors from the postings by `token_id`. The tokens near them come from
   `tape_tokens` of the anchors' contigs through `bb_slot_tokens(contig_uid, slots)`, a Python function
   that `open_catalog` registers. It is a footer-statistics point lookup (`KeyLookup`, as B-13 and B-14),
   called once per anchor contig. `profile_query` reports its reads as one more `tape_tokens` scan.
3. `open_catalog` shadows the three lookup macros, on its own connection, with temporary macros over
   point lookups in the connection's `SeqStore`. These are counted in `profile_query`'s `index`. A plain
   `duckdb` connection keeps the stored SQL macros, which give the same answers without pruning.

No table, schema or builder changed. No pair is materialized, and `partners_of` returns the same
columns as before.

## `rL.3`: the current catalog code cannot serve it

*Measured*, wilkcomp02, 2026-11-10. This was run read only: `/tmp/if2453/b31/exp/rl3.py` under
`systemd-run --user --scope -p MemoryMax=32G -p MemorySwapMax=0 nice -n 19` in
`bio-box-test_74b3e05.sif`, and it wrote only under `/tmp/if2453/b31/rl3`.

- `manifests/release=rL.3.json` is manifest version 2. It has no `stages` and lists 77 files, which
  include `slot_protein`, a table the registry no longer has (card B-17).
- `catalog.generate_catalog(manifest, root=<rL>)` raises `KeyError: 'slot_protein'`. The current code
  therefore cannot build a catalog for `rL.3`, and has nothing to open.
- The stored `catalog/release=rL.3.duckdb` (a copy, opened read only) holds the old macros. In it,
  `partners_of` returns 0 rows for the top NR50 token (`3a04ab75…`, 3,509 slots) and for the terminase
  family `PF03237.22`, and `token_rows` does not exist. So the old catalog cannot give a "before" for these macros.

Following the card, before and after are measured on the 1,000x synthetic store instead.

## 1,000x synthetic store, before and after

*Measured*, wilkcomp02, 2026-11-10. `scripts/stress/synth_scale.py <tree> <out> 1000 --seed 1` was run
in the same image and scope (`systemd-run --user --unit=b31-{before,after2}-1000`), once with the tree
at `22b55c4` (before) and once with this branch (after). The builders are unchanged, so both builds
give the same manifest (`manifest_hash` `44f284301445751ad9ca2e4e418bb96c`). The store has 3,000 contigs
and 854,000 slots, at the real row-group caps. `tape_tokens` has 7 row groups, `token_postings` 20,
`contig_index` and `seq_locator` 1 each.

Each lookup below was then run alone on a fresh connection (cold lookup caches) with
`/tmp/if2453/b31/exp/q.py <out> v3.py`. "rg" is row groups read of the scan's total
(`catalog.profile_query`). For a point lookup, rg is the lookup's own count: `index` for the sequence
store, the `tape_tokens` entry for `bb_slot_tokens`.

| Lookup | Before: rg read | After: rg read | Result (the same before and after) |
|---|---|---|---|
| `token_rows(token_of('nr90', k))` | postings 20 of 20; `token_dict` 1 of 1, 3 scans | postings 1 of 20; `token_dict` 1 of 1, 2 scans | (10, 560) |
| `q1_count(token_of(a), token_of(b), 3, NULL)` | postings 20 + 20 of 20 | postings 1 + 1 of 20 | 10 |
| `pairs_of(token_of(a), token_of(b), 10)` | postings 20 + 20 of 20 | postings 1 + 1 of 20 | (10, -220) |
| `partners_of(<family>, 10)`, literal; 1,000 loci on 1,000 contigs | anchors 1 + 1 of 20; partner side postings 20 of 20 (1,482,520 rows) | anchors 1 + 1 of 20; partner side `tape_tokens` 7 of 7 | (2,214 tokens, 57,000 pairs) |
| `partners_of(token_of('nr90', k), 10)`; 10 loci on 10 contigs | postings 20 + 20 of 20 (anchors); partner side 20 of 20 | anchors 1 + 1 of 20; partner side `tape_tokens` 5 of 7 | (54, 570) |
| `partners_of(<nr90>, 10)`, literal | anchors 1 + 1 of 20; partner side 20 of 20 | anchors 1 + 1 of 20; `tape_tokens` 5 of 7 | (54, 570) |
| `partners_of(<hub family>, 10)`; 26,000 loci (not a key lookup) | anchors 1 + 1; partner side 20 of 20 | anchors 1 + 1; `tape_tokens` 7 of 7 | (36,222, 1,469,000) |
| `default_digest(contig_uid) FROM pick1` | `contig_index` scan 1 of 1 (all 3,000 rows) | no scan; index lookups 0 + 1 | 1 |
| `digest_len(default_digest(contig_uid)) FROM pick1` | `seq_locator` 1 of 1 + `contig_index` 1 of 1 scans | no scan; lookups 1 + 1 | 168,903 |
| `contig_len_of(contig_uid) FROM pick1` / `pick1000` | `contig_index` scan 1 of 1 | no scan; lookups 1 + 1 | 168,903 / 275,448,318 |

At 1,000x each index has a single row group, so the digest macros' "1 of 1" before already reads
everything. The 100x test with scaled caps shows the difference: 10 of 10 before
(`tests/r2/test_scale_layout.py`, D-05 markers), and at most 2 lookup row groups per key after.

Wall time, *measured*, in the same runs. "Script" means `synth_scale.py`'s own report, run in its order
on one connection, so a later partner query can find its contigs already cached. "Cold" means the
`q.py` runs above, one fresh connection per query. Both use 4 DuckDB threads.

| Query | Before, script | After, script | Before, cold | After, cold |
|---|---|---|---|---|
| postings, key through `token_of` | 0.275 s | 0.044 s | 0.321 s | 0.088 s |
| Q1, keys through `token_of` | 0.593 s | 0.099 s | 0.594 s | 0.094 s |
| Q4, keys through `token_of` | 0.648 s | 0.141 s | 0.610 s | 0.112 s |
| Q3 partners of one family, literal key | 0.259 s | 0.317 s | 0.252 s | 0.373 s |
| Q3 partners of one NR90 through `token_of` | 0.586 s | 0.101 s (warm) | 0.552 s | 0.209 s |
| Q3 partners of the hub (26,000 loci) | 0.846 s | 0.794 s | 0.848 s | 1.015 s |
| export proxy, 500 whole contigs (`contig_len_of` inside) | 1.130 s | 1.316 s | | |
| whole script, `total_s` | 310.7 s | 324.2 s | | |
| `peak_rss_mb_queries` | 1,114 MB | 1,078 MB | | |

On this small store, a partners query whose anchors cover the whole release (family, hub) costs 0.1 to
0.17 s more than before: the partner side is Python work per anchor contig, against DuckDB's scan of
1.5 to 2.2 M postings rows. That holds for these timings only (*measured*, 3,000 contigs); what
scales is the row-group count.

The partner side now costs Python work per anchor contig. A first version that looked up one
(contig, slot) per row took 286.1 s for the hub and 4.4 s for the family (the same script, *measured*;
`cProfile`: 12,015 pyarrow calls per family query). The current version decodes a row group once to
numpy and calls the function once per contig with all of that contig's slots: 0.79 s and 0.32 s.

## The partner side grows with the anchors, not with the release

*Measured*, wilkcomp02, 2026-11-10, `sh scripts/pod_pytest.sh B-31j tests/r2/test_scale_layout.py -k
grow_with_the_anchors -o junit_family=xunit1 --junitxml=...` (the property `partner_reads`). The test
builds the D-05 store at 100x and at 200x with the scaled caps. Copies 0 to 99 are the same contigs at
both sizes, and the test runs `partners_of` of the NR90 token of copies 0 to 9 (10 anchor contigs) on
a fresh connection.

| Store | `tape_tokens` rg read | `token_postings` rg read (anchors, hub check) | Before B-31: partner-side postings rg |
|---|---|---|---|
| 100x | 6 of 11 | 1 of 32, 1 of 32 | 32 of 32 (*measured*: this test's failing run before the fix) |
| 200x | 9 of 22 | 1 of 63, 1 of 63 | 63 of 63 (*derived*: the old near side read every postings row group at 100x and 1,000x) |

At both sizes the row groups read equal the row groups whose footer `contig_uid` range holds an anchor
contig (the test's oracle), at most one per anchor contig. Results resolved to stable keys through
`token_dict` are equal at both sizes. The count rose from 6 to 9 because, at 100x, several anchors
shared a row group; doubling the store moved the row-group boundaries. It is bounded by the 10 anchor
contigs, not by the release's 11 or 22.

## What still does not meet 2 row groups per key

- `partners_of` reads one `tape_tokens` row group per anchor contig, at most. For a token on many
  contigs this is more than 2 per token: 7 of 7 for the family token with 1,000 loci at 1,000x, and 6
  of 11 for the NR90 token at 100x. The two D-05 tests that score `partners_of` per token
  (`test_partners_of_near_side_prunes_literal_key`, `test_partners_of_token_of_key_prunes`) keep their
  `xfail` markers, and they now measure on a fresh connection. In the script's order the NR90 query ran
  warm, after the family query had cached its contigs, and read 0 row groups.
- `token_of(level, key)` itself scans `token_dict`, which is sorted by `token_id` and not by
  `(level, key)`. It read 1 of 1 row group here, but at `rL.3` `token_dict` holds 11,258,843 rows
  (`2026-11-10_import_ladder/rL3.md`). The card does not cover this; it is an open question.
