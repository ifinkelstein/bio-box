# Resolving a positive control by accession on `rN.3` (card W-15), 2026-11-12

Card W-15 asks first how long `search hmm --positive ACCESSION:START-END` and `bio-box find --accession` take to
resolve an accession on a real release, and how many row groups of the label relation they read. If the lookup
takes over a second, a per-stage table sorted by accession is added; if not, the card ends with this note. **It
takes 0.03 to 0.05 s and reads one row group per stage, so no table was added.**

Every number is *measured* on wilkcomp02, 2026-11-12, image `bio-box-test_74b3e05.sif`, tree `task/W-15`
(= `plan/task-cards` at `dda23ce`, no code change), unit `w15-measure` under `systemd-run --user -p MemoryMax=16G
-p MemorySwapMax=0 -p CPUQuota=200% nice -n 19`, DuckDB at 2 threads, the root bound read-only, unless labelled
otherwise. The pod's load average was about 39 during the run (other agents' builds), so wall times are noisy;
row-group counts are not. Driver: [`2026-11-12_positive_lookup_rN3/w15_measure.py`](2026-11-12_positive_lookup_rN3/w15_measure.py)
`/tmp/if2453/r10/root rN.3 <out>`; raw output: [`w15_rN.3.json`](2026-11-12_positive_lookup_rN3/w15_rN.3.json).

## 1. What the lookup reads on `rN.3`

The card names `rN.1`; this used `rN.3`, the release over the same stages (S1, P1, P2) built by current code
(card C-17). Its stages were ingested at `740d049`, before card B-27, so the manifest lists **no `contig_label`
file**: `export.open_store` reads `contig_occurrence` instead and prints the slow-path note. `rN.1` reads the same
three files, so the numbers hold for it too (not run separately).

| stage | `contig_occurrence` rows | row groups |
|---|---:|---:|
| S1 | 34,076 | 1 |
| P1 | 539 | 1 |
| P2 | 5,018 | 1 |

`rN.3` has 32,815 contigs (`contig_dim` rows); S1 has 34,076 records for 27,258 contigs (several accessions per
contig, e.g. RefSeq and GenBank copies of one phage).

## 2. Timings

The lookup itself, `SELECT DISTINCT contig_uid FROM lab WHERE accession = ...` as `cli._contig_labels` runs it,
under `catalog.profile_query`:

| accession | seconds | row groups read | bytes read |
|---|---:|---:|---:|
| PQ006425 (S1) | 0.047 | 3 of 3 | 832,975 |
| NC_054918 (S1) | 0.034 | 3 of 3 | 832,975 |
| MK327939 (S1, same contig as NC_054918) | 0.045 | 3 of 3 | 832,975 |
| CP089929.1 (P2, the longest contig, 10,737 genes) | 0.040 | 3 of 3 | 832,975 |
| header substring `Escherichia coli` (38 contigs) | 0.042 | 3 of 3 | 1,158,313 |

The whole resolution in process (open the store, label lookup, read the tapes, pick the gene):

| accession | `_contig_labels` s | `find_genes` s | `_positive_uid` s (what `search hmm --positive` calls) | genes |
|---|---:|---:|---:|---:|
| PQ006425 | 0.194 | 0.274 | 0.160 | 4 |
| NC_054918 | 0.111 | 0.145 | 0.152 | 266 |
| MK327939 | 0.101 | 0.147 | 0.155 | 266 |
| CP089929.1 | 0.076 | 0.299 | 0.310 | 10,737 |

`bio-box find --accession ACC:START-END --threads 2` as a user runs it, a fresh process each: 2.06 to 2.67 s
(4 runs, each exit 0 with one gene). `python -c 'import biobox.cli'` alone: 1.98 s. Nearly all of a `find`'s wall
time is Python start-up and imports, not the lookup. `GUIDE.md` §9 says "about 4 s on `rN.1`"; that figure has no
label and was not reproduced here.

With 3 row groups in all and one per stage, an index sorted by accession would read the same 3 row groups: on
this release it can save nothing.

## 3. Where it stops being cheap

- `rM.2` (the 50 Gbp release) is about five times `rN.3` in bases (*estimate*, from the brief). Its genomes are
  complete or near-complete, so its contig count grows far less than its bases: the 50 Gbp set has 54,265
  contigs (*measured*, card C-20). At five times `rN.3`'s 39,633 label rows, about 200,000 rows (*derived*, upper
  bound), the scan reads about 4 MB and takes about 0.2 s (*derived*, linear from 0.04 s and 832,975 bytes per
  39,633 rows), and each stage file still holds one or two row groups of 122,880 rows.
- The scan crosses one second at about 1 M label rows (*derived*: 0.04 s per 39,633 rows, about 1 µs per row at 2
  threads). Stage S6 has 111.5 M contigs (*measured*, AGENTS §4.1); its `contig_label` alone would be about 908 row
  groups (*derived*: 111.5 M / 122,880) and a scan about 110 s and 2.3 GB read (*derived*, linear). The table this
  card describes is needed before a draft-assembly (MAG) stage, not before `rM.2`.

## 4. Found while measuring

- On a release built after card B-27, `contig_label` holds one accession per contig (its first record's), while
  `contig_occurrence`, read on `rN.3`, holds every record's. Both NC_054918 and MK327939 resolve to the same
  contig on `rN.3` because the fallback reads `contig_occurrence`; on a B-27 release only the first record's
  accession would, and the other would match nothing (from the B-27 definition; not run, no such release over
  these stages exists). A lookup table sorted by accession should hold every record's accession.
- `GUIDE.md` §9 says `rN.3` will not need the slow path for contig names (card C-17). It does: `contig_label` is
  written per stage at ingest, and `rN.3`'s stages were ingested before card B-27.
