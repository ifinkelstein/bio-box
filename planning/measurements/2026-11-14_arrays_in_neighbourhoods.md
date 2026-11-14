# CRISPR arrays beside a locus (card E-14)

Date: 2026-11-14, wilkcomp02. Every number carries a label: *measured* (host, date, command or record), *derived*
(formula and inputs) or *estimate*.

## What changed

- `bio-box nbr --format tsv` (so also `search hmm`'s `loci/`) and `nbr summarize --format tsv` show the CRISPR arrays
  of card E-12 beside each locus, on a release that holds `l1_crispr` files (`src/biobox/search/arrays.py`). Arrays
  are never tape slots and never carry an ordinal; they are rows of their own after the gene rows (a section that
  opens with `# crispr_arrays<TAB>window_start=..<TAB>window_end=..`), and one more table at the end of the summary.
  `search hmm`'s README.txt gains a paragraph on them on such a release. The JSON forms are unchanged (the record
  schemas are locked, `tests/registry/json_schema_lock.json`, which this card may not change).
- The C-08 provenance check (`src/biobox/provenance.py`) knows the producer `minced-<v>-w<n>-nr<n>-r<a>x<b>-s<a>x<b>`
  and checks the `crispr-*` step: `minced` and `pyarrow` against the lock, `java` recorded and non-empty (the JVM is
  not pinned), `minced_args` equal to the arguments the producer's name encodes.
- The CRISPR-Cas example gains array checks marked `[requires l1_crispr]`, which the shared checker reports `n/a` on a
  release without arrays (`examples/_shared/check.py`).

## The window in bases

Per locus (one anchor gene): a gene window (`--window Ng`) spans from the least start to the greatest end of the
locus's rows (the anchor and the genes the window shows, after `--strand` and `--operon-gap`; a row reached across a
circular origin moved by one contig length to the anchor's side); a base window (`--window W`) spans
`[anchor start - W, anchor end + W)`. An array is within the window when it shares at least one base with the span,
`partial` when part of it lies outside; on a circular contig (tape `is_circular`, `--circular auto`) it is placed at
the shift (-n, 0, +n) that puts most of it inside (ties: the smaller gap to the anchor, then no shift). A gene in a
base window is shown when its gap to the anchor is at most W, so a gene starting exactly W bases after the anchor's
end is shown and an array starting there is not (the one-base difference of a half-open span).

## Real data: a scratch view of rN.3 with card E-12's arrays

Built by `/tmp/if2453/e14/bin/mkview_e14.py /tmp/if2453/r10/root /tmp/if2453/e12/view /tmp/if2453/e14/view` in the
runtime image of `faec1f4` (`/tmp/if2453/images/bio-box_faec1f4.sif`, sha256 `d67300aa...f098a`, checked by the
driver), inside the user unit `e14-real` (MemoryMax 16G, CPUQuota 200%, nice 19), with the data root and E-12's view
bound read-only. The view holds symlinks only: every top-level directory of the root but `l1`, `ledger`, `manifests`,
`catalog` and `runs`; per table of `l1` (and per producer of `l1/feature`, E-12's `minced-` producer added);
`l1/crispr` -> E-12's; per stage of `ledger` (per step for P1, P2, S1, E-12's crispr steps added). Its own
`manifests/release=rN.3.json` is the root's manifest of rN.3 with `root` set to the view and 9 entries added (the
`feature` and `l1_crispr` files and the crispr ledgers of P1, P2, S1, made as `catalog.build_manifest` makes them);
its own `catalog/release=rN.3.duckdb` is `catalog.generate_catalog` of it. Nothing was written to the data root or
to `/tmp/if2453/e12`. Cost: 10.9 s wall, 6.5 CPU-s, 408 MB (*measured*, `/usr/bin/time -v`).

**C-08 on that manifest**: no problem over its 69 files, the three real crispr ledgers included (*measured*,
`mkview.json`).

**The control loci** (*measured*, `/tmp/if2453/e14/bin/e14_real.py`, `bio-box nbr <uids> --release rN.3 --format
tsv` on the view, uids from `bio-box find`: TniQ `28e32e6533f7f55e5bca64f44c89b426`, Cas7f
`cf65ec9336feb980ebf09b593078bf7d`, the complete system's Cas7f `07ba31e0f0e5ffcc30388139b5ec30e4`; each one locus,
all on `NZ_AP026781.1`):

| Locus | Window 10g (bases) | Arrays within it |
|---|---|---|
| TniQ 426340-427525 | 413722-444027 | none |
| Cas7f 423377-424433 | 411868-440339 | none |
| complete I-F Cas7f 687857-688841 | 676009-697441 | 679543-681134 (27 repeats, 26 spacers, `side_rel -1`, delta 6,723) and 689534-690762 (21 repeats, 20 spacers, `side_rel 1`, delta 693), repeat `CTTCACTGCCGAGTAGGCAGCTTAGAAA`, neither partial |

The same at 20 kb and 10 kb windows: the complete system's two arrays, none for the Tn7-like element (E-12: its
array has 2 repeats, below MinCED's default of 3). `nbr summarize` over the three: 1 of 3 loci with an array, 2
arrays. Write lines: `3 tables, 65 rows` at 10g (63 gene rows and 2 array rows).

**The example's two runs** (`/tmp/if2453/ex2/rN.3-b/runs`, read-only; `nbr <run> --window 10g` on the view):
13 of the 635 TniQ loci and 138 of the 157 Cas7f loci have an array within the window; 173 distinct arrays near the
Cas7f loci (*measured*, `nbr summarize` array table). The example's `sql/checks_arrays.sql` on the view's catalog:
`arrays_pdam_complete_if_window_has_array` true, `arrays_pdam_complete_if_window` 2, control TniQ and Cas7f windows
0, nearest array to the control Cas6f 256,168 bases, 7 of the 20 I-F3 TniQ loci, 6 of the 19 Cas1-less Cas7f loci
with a TniQ, 78 of 78 complete prokaryote I-F loci, 45 of 46 phage Cas7f loci, ICP1's Cas7f 1 (*measured*).
Cross-check: the (locus, array) pairs `nbr` wrote into the 792 locus tables equal those the SQL definition lists
from the runs' windows, 190 and 190, none on one side only (*measured*, `/tmp/if2453/e14/bin/e14_cross.py`).

**Reads** (`catalog.profile_query` scans of `arrays.read`, *measured*, `e14_cross.py`). On the view, `feature` of
the minced producer is 3 files of 1 row group (286, 2,291, 213 rows), `l1_crispr` 3 files of 4 row groups,
`contig_tapes` 95 row groups. "Holding" is the row groups whose footer `contig_uid` range holds a locus contig.

| Loci | Contigs | feature scanned / holding / total | l1_crispr | contig_tapes | CPU-s |
|---|---|---|---|---|---|
| 1 (complete I-F Cas7f) | 1 | 3 / 3 / 3 | 3 / 3 / 4 | 3 / 3 / 95 | 0.30 |
| 157 (Cas7f run) | 149 | 3 / 3 / 3 | 4 / 4 / 4 | 95 / 87 / 95 | 0.35 |
| 635 (TniQ run) | 406 | 3 / 3 / 3 | 4 / 4 / 4 | 95 / 95 / 95 | 0.31 |

One query per table for all loci. A file of one row group spans every contig uid, so it is always read (the array
tables are small: 2,790 arrays and 145,148 elements in all, *measured*). For the Cas7f run's 149 contigs the
`contig_tapes` scan read 95 row groups where 87 hold a contig: above as many keys as the largest table has row groups
DuckDB pushes only the keys' min and max, the rule of card B-27 (`catalog.set_join_key_filter`; fewer row groups
than keys, so under one row group per key). The array step adds about 0.3 CPU-s to an `nbr` of 24 CPU-s (*derived*:
0.35 / 24.05 = 1.5%, the Cas7f run's `nbr --format tsv`).

**Cost of the real-data pass** (`e14_real.py`, all commands above in one process): 152 s wall, 275 CPU-s, 4.25 GB
peak (*measured*, `/usr/bin/time -v`), the largest step `nbr summarize` of the TniQ run (57 CPU-s, 4.15 GB, as
before the card: 43.5 CPU-s and 3.8 GB unit peak at 20 kb in the example's run).

## Fixture

`tests/r2/test_e14_arrays_nbr.py`: the micro-S0 search world with planted arrays (feature row groups of 4 rows,
l1_crispr of 8). Reads (*measured*, wilkcomp02, 2026-11-14, `test_arrays_are_read_by_contig_whatever_the_number_of_loci`,
the test's `e14 scans` line), scanned of total row groups, summed over the two producers' files:

| Loci | Contigs | feature | l1_crispr | contig_tapes |
|---|---|---|---|---|
| 1 (a lambda gene) | 1 | 3 of 27 | 21 of 165 | 1 of 1 |
| 62 (every lambda gene) | 1 | 3 of 27 | 21 of 165 | 1 of 1 |
| 854 (every gene of the three contigs) | 3 | 27 of 27 | 165 of 165 | 1 of 1 |

The test asserts that the scans equal the row groups whose footer range holds the loci's contigs, and that 62 loci
read what one locus reads.

## Byte identity on a release without arrays

On the fixture without arrays: the array reader is never called, and every output (locus tables on stdout and in
`--out`, the `wrote` line, summary TSV and JSON) equals the output with the array path cut off, both in one process
(`test_on_a_release_without_arrays_nothing_is_read_and_every_output_is_the_unflagged_one`); `search hmm`'s README is
the old text. On the real release rN.3 (no arrays), old code against new: this branch's `nbr <run> --window 10g
--format tsv` and `nbr summarize <run> --window 10g --by pfam,nr50 --format tsv` of the example's two runs equal,
byte for byte, the `loci/` tables (635 and 157) and `neighbors.tsv` files `search hmm` wrote with the code before the
card on 2026-11-13 (`/tmp/if2453/ex2/rN.3-b/{tniq,cas7f}`; only the `RUN` file card W-29 adds to an `--out`
directory left out) (*measured*, `/tmp/if2453/e14/bin/e14_ident.py`, unit `e14-ident`, 93 s wall).
