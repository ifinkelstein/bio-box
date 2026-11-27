# Bytes that are not bases: the inputs of the built stages scanned, and `GV1` re-staged without its trailing gaps (card R-23)

- **Date and host:** 2026-11-27, `wilkcomp02`. Every number carries a label (AGENTS §6): *measured* (on `wilkcomp02`,
  that day, by the command or file named), *derived* (formula and inputs), or *estimate*.
- **Scripts:** card R-19's `scripts/stress/build_stage_pod/r19_fasta_scan.py`, used as it is (its copy on the pod,
  `/tmp/if2453/r23/bin/r19_fasta_scan.py`, sha256 `dd2a0f5f…`, is the tree's file); for objects that are no longer plain
  files on disk (tar members, a file the archiver compressed) a small pod-side wrapper,
  `/tmp/if2453/r23/bin/r23_stream_scan.py`, that imports the scanner's `_Scan`, `_HashingRaw` and `_compression`
  unchanged and feeds them the object's bytes (§1.2). The re-stage is
  `scripts/stress/build_stage_pod/r23_gv1_trailing_gaps.py` (tests `tests/r1/test_r23_gv1_trailing_gaps.py`); template
  `scripts/stress/build_stage_pod/stages/templates/GV1.conf`.
- **Read only:** the root `/tmp/if2453/r10/root` (its `runs/`, `l0/l0_manifest/`, `l0seq/contig_index/` and the sizes of
  `l0seq/seq_shards/`), `/tmp/if2453/s2/inputs/` (`V2`, `PL2`, `GV1`, the chain), the copies of earlier inputs under
  `/tmp/if2453/s3/` and `/stor/scratch/Wilke/if2453/all-the-data/provisional/`, and the reports of cards R-19
  (`/tmp/if2453/r19/`) and R-20 (`/tmp/if2453/r20/`). Nothing in the root, under `/tmp/if2453/s2/inputs/` or under
  `/tmp/if2453/stages/` was written; nothing under `/stor` was overwritten or deleted.
- **Written:** under `/tmp/if2453/r23/`: `scan/` (one JSON per object, `<STAGE>__<object>.json`, and one summary per
  stage and source, `<STAGE>[.<source>].stage.json`), `built_inputs.json`, `l0_manifest.json`, `stored.json`, `bin/`
  (the helpers), `logs/`.

## 0. Summary

@@SUMMARY@@

## 1. The built stages and where their inputs are

### 1.1 Which stages, which objects

The root stores 16 stages (`ls /tmp/if2453/r10/root/l0seq/contig_index`: `S1`, `P1` to `P6`, `H1`, `PL1`, `C5`, `V1`,
`G1`, `SPX1`, `V2`, `PL2`, `PL3`; *measured*); `PL3` is the wave's stage now building (past ingest; its Pfam scan
runs), and no other wave stage has started. For each stage the objects ingest read are named twice in the root: the
`--objects` list of its ingest run record (`runs/stage=<S>/ingest.json`, `argv`) and the rows of its
`l0/l0_manifest/stage=<S>/` (object name, bytes, `object_sha256`) (*measured*, `/tmp/if2453/r23/built_inputs.json` and
`l0_manifest.json`, the pod's Python and the image's `pyarrow`).

| Stage | Ingest started (UTC) | `--objects` of the run record | Still there? | Objects | Bytes | Copy scanned (§2) |
|---|---|---|---|---:|---:|---|
| `S1` | 2026-11-10 16:57 | `/tmp/if2453/r10/inputs/S1/objects.json` | no | 1 | 2,121,124,779 | member `inputs/S1/phage_first34076.fa` of `provisional/r50_inputs/r10_inputs.tar` |
| `P1` | 2026-11-10 18:13 | `/tmp/if2453/r10/inputs/P1/objects.json` | no | 250 | 1,029,384,352 | members of `provisional/r50_inputs/ladder_genomes_P1.tar` |
| `P2` | 2026-11-10 18:43 | `/tmp/if2453/r10/inputs/P2/objects.json` | no | 2,247 | 9,716,599,273 | members of `provisional/r50_inputs/ladder_genomes_P2.tar` |
| `P3` | 2026-11-11 19:22 | `/tmp/if2453/r50/inputs/P3/objects.json` | no | 2,671 | 10,074,413,081 | members of `provisional/r50_inputs/genomes_P3.tar` |
| `P4` | 2026-11-12 02:26 | `/tmp/if2453/r50/inputs/P4/objects.json` | no | 2,329 | 10,101,497,437 | members of `provisional/r50_inputs/genomes_P4.tar` |
| `P5` | 2026-11-12 03:27 | `/tmp/if2453/r50/inputs/P5/objects.json` | no | 2,289 | 10,092,827,773 | members of `provisional/r50_inputs/genomes_P5.tar` |
| `P6` | 2026-11-12 04:29 | `/tmp/if2453/r50/inputs/P6/objects.json` | no | 2,180 | 10,087,447,605 | members of `provisional/r50_inputs/genomes_P6.tar` |
| `H1` | 2026-11-15 13:44 | `/tmp/if2453/s2/inputs/H1/objects.json` | no | 1 | 1,962,116,403 | `provisional/s2_inputs/HPGC.7230.genomes.final.fasta.gz`, decompressed |
| `PL1` | 2026-11-16 03:06 | `/tmp/if2453/s2/inputs/PL1/objects.json` | no | 1 | 7,425,668,755 | `provisional/s2_inputs/plsdb_2024_05_31_v2_sequences.fasta` |
| `C5` | 2026-11-19 17:24 | `/tmp/if2453/s2/inputs/C5/objects.json` | no | 2,032 | 2,781,209,046 | @@C5COPY@@ |
| `V1` | 2026-11-20 15:28 | `/tmp/if2453/s2/inputs/V1/objects.json` | no | 1 | 5,228,187,536 | @@V1COPY@@ |
| `G1` | 2026-11-22 13:20 | `/tmp/if2453/s2/inputs/G1/objects.json` | no | 7,139 | 2,903,786,890 | `/tmp/if2453/s3/globdb/genomes/G1/files/` (card C-41's download, every object by name and bytes) |
| `SPX1` | 2026-11-23 10:19 | `/tmp/if2453/s2/inputs/SPX1/objects.json` | no | 304 | 75,267,693,104 | **none on the pod** (§3) |
| `V2` | 2026-11-26 02:01 | `/tmp/if2453/s2/inputs/V2/objects.json` | yes | 2 | 9,508,728,867 | the objects themselves |
| `PL2` | 2026-11-27 02:34 | `/tmp/if2453/s2/inputs/PL2/objects.json` | yes | 8 | 16,580,610,887 | the objects themselves |
| `PL3` | 2026-11-27 05:16 | `/tmp/if2453/s2/inputs/PL3/objects.json` | yes | 2 | 18,166,241,620 | scanned by card R-19 (its note §3; not again) |

`provisional/` is `/stor/scratch/Wilke/if2453/all-the-data/provisional/` (the archive of the first stages' inputs, made
2026-11-11 to 2026-11-13). Objects and bytes are the manifest's (*measured*).

So the `objects.json` that ingest read is still on the pod for `V2`, `PL2` and `PL3` only. For the others a copy of
every object was looked for by name and size (a `find` over `/tmp/if2453` and the `/stor` areas of this project,
stopped once the archived copies were found, and the listings of `/tmp/if2453/s3/` and
`/stor/scratch/Wilke/if2453/all-the-data/provisional/`), and each copy scanned was then checked against the
manifest's `object_sha256` (the scanner hashes the object's bytes as it reads them): **every object scanned is,
byte for byte, the object ingest read** (§2).

### 1.2 How objects that are not plain files were scanned

`r19_fasta_scan.py` takes a path and reads a plain or gzip file. Four kinds of copy are not that: the members of a
tar archive (`P1` to `P6`, `S1`, `C5`), a gzip file whose decompressed bytes are the object (`H1`), and a zstd file
whose decompressed bytes are the object (`V1`). `r23_stream_scan.py tar` reads a tar archive as a stream and, for each
member whose base name is an object of the stage, runs the scanner's own `_Scan` over the member's bytes through its
`_HashingRaw` (so the sha256 is the member's, comparable with the manifest) and `_compression` (a member that is itself
gzip, as `G1`'s and `C5`'s objects are, is decompressed by the scanner's rule); `r23_stream_scan.py file --decompress`
does the same over a decompressed stream. The report has `scan_file`'s fields. A self-test on the pod
(`r23_stream_selftest.py`, `r23_zstd_selftest.py`: a clean file, one with trailing and inner gaps, a gzip member with a
header that is not UTF-8, a byte `0xe5` and an empty record, a gzip and a zstd file) gave **every field equal to
`r19_fasta_scan.scan_file` on the same bytes** (`ALL EQUAL`, `ZSTD EQUAL`; *measured*). Memory stays bounded: one
object at a time, read in pieces.

## 2. The scan of the built stages

Detached units at `nice -n 19`, `CPUQuota=100%`, `MemoryMax=8G`, `MemorySwapMax=0`, `IOSchedulingClass=idle`, one
process each, at most four of this card at a time (two while a unit of card R-21 ran), each started after a look at
the active units and with no `STOP` file in the root; the wave was in `PL3`'s Pfam scan throughout, and no `MV1` or
`V3` unit ran. Units `r23-scan-1` (`V2`), `-2` (`PL2`), `-3` (`G1`, `scan-stage` over the 7,139 files), `-4` (`PL1`,
`H1`, `S1`, `P1`, `P2`), `-5` (`P3` to `P6`; its `V1` step stopped on a zstd window, below), `-6` (`V1`), `-7` (`C5`),
05:16 to 06:00 CST.

| Stage | Objects | Records | Bases | Not UTF-8 | Not ASCII | Non-IUPAC records | Non-IUPAC bytes | Empty | Before header | Read errors | Byte values | sha256 = manifest | Status |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---|
| `S1` | 1 | 34,076 | 2,084,350,440 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | none | 1 of 1 | clean |
| `P1` | 250 | 539 | 1,016,625,803 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | none | 250 of 250 | clean |
| `P2` | 2,247 | 5,018 | 9,596,180,328 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | none | 2,247 of 2,247 | clean |
| `P3` | 2,671 | 9,108 | 9,949,084,485 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | none | 2,671 of 2,671 | clean |
| `P4` | 2,329 | 4,006 | 9,976,464,003 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | none | 2,329 of 2,329 | clean |
| `P5` | 2,289 | 4,112 | 9,967,897,315 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | none | 2,289 of 2,289 | clean |
| `P6` | 2,180 | 4,228 | 9,962,553,301 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | none | 2,180 of 2,180 | clean |
| `H1` | 1 | 7,230 | 1,961,987,211 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | none | 1 of 1 | clean |
| `PL1` | 1 | 72,556 | 7,297,669,596 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | none | 1 of 1 | clean |
| `C5` | 2,032 | 5,228 | 9,392,003,904 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | none | 2,032 of 2,032 | clean |
| `V1` | 1 | 168,536 | 5,222,288,776 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | none | 1 of 1 | clean |
| `G1` | 7,139 | 1,227,287 | 10,051,402,938 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | none | 7,139 of 7,139 | clean |
| `V2` | 2 | 262,013 | 9,336,787,044 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | none | 2 of 2 | clean |
| `PL2` | 8 | 270,652 | 16,319,048,429 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | none | 8 of 8 | clean |
| `PL3` (card R-19, cited) | 2 | 792,463 | 17,842,081,525 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | none | 2 of 2 | clean |
| `SPX1` | inputs gone (§3) | | | | | | | | | | | | stored shards scanned (§3) |

All *measured* (`/tmp/if2453/r23/scan/<STAGE>__<object>.json`, one per object, summed by
`/tmp/if2453/r23/bin/r23_table.py`; `PL3` from `/tmp/if2453/r19/r26_11_5/scan/`; the sha256 column by
`r23_match.py check` against `l0_manifest.json`). **In the 21,151 objects of the 15 stages built from inputs still on
the pod there is no header that is not UTF-8 or not printable ASCII, no byte outside the IUPAC letters, no empty
record, no line before a header and no read error.** Every object was read to its end (its sha256 is that of its
whole bytes, equal to the manifest's). Records and bases
are of the objects as read; the root stores fewer where a record repeats another or another stage owns it (DEC-26,
DEC-27): for example `P3` 9,108 records, 9,107 stored contigs. Since ingest stores a record's bytes as read (card
R-20 §8.3) and the staging of these stages wrote no byte into an object afterwards, a clean object means a clean
stored contig; so **the catalog holds no byte that is not a base in these 15 stages**.

Per object (one JSON each; too many to list here; the stages of one to eight objects): `S1`
`phage_first34076.fa` (sha256 `e62d5658…`); `H1` `HPGC.7230.genomes.final.fasta`; `PL1`
`plsdb_2024_05_31_v2_sequences.fasta`; `V1` `uhgv_v1.0_votus_full.fna` (sha256 `9c595d26…`); `V2`
`imgvr_v4.1_hc_hq_votu_reps.fna` 256,049 records, 9,287,551,193 bases (`68795441…`), `inovirus_roux2019.fna` 5,964,
49,235,851 (`209f5a17…`); `PL2` `plasmidscope_RefSeq.fasta` 85,995 / 6,948,506,952, `plasmidscope_GenBank.fasta`
91,774 / 6,840,061,592, `plasmidscope_COMPASS.fasta` 12,084 / 1,260,922,758, `plasmidscope_DDBJ.fasta` 5,305 /
446,998,007, `gut_plasmids_yu2024.fasta` 68,350 / 397,149,786, `plasmidscope_ENA.fasta` 6,239 / 375,652,221,
`plasmidscope_Kraken2.fasta` 898 / 49,062,410, `plasmidscope_TPA.fasta` 7 / 694,703; every count of every one 0.

The `V1` copy is zstd with a 2^31 window (`--long`), which Python's default decoder refuses ("Frame requires too much
memory"); the first try (in `r23-scan-5`) stopped there and wrote nothing, and `r23-scan-6` read it with
`window_log_max=31` (`r23_stream_scan_zstdlong.py`, the same wrapper with that one option).

## 3. Built stages whose inputs are not on the pod, and what a scan of their stored sequences would cost

**`SPX1`** is the one built stage none of whose inputs is on the pod: its 304 SPIRE v1 assemblies (75,267,693,104
bytes, `*.fa.gz`, card C-48e part 1) were hard links in `/tmp/if2453/s2/inputs/SPX1/assemblies/` to
`/tmp/if2453/c48e1/dl/`, and both are gone (`/tmp/if2453/c48e1/dl/` holds only the fetch log); no file of their names
is under `/tmp/if2453/s3/`, `/tmp/if2453/c48e1/`, `/stor/scratch/Wilke/if2453/all-the-data/`,
`/stor/scratch/Wilke/if2453/bio-box-staging/` or `/stor/work/Wilke/if2453/all-the-data/` (*measured*, `find` by name,
2026-11-27), and none was archived to `provisional/` (made before `SPX1` was staged). `tower` was not searched (the card
asks for the pod).

What the root stores of it can be scanned instead: the stage's 40 sequence shards
`l0seq/seq_shards/stage=SPX1/*.fa.gz` (BGZF, contig-uid headers) hold exactly the bytes ingest kept, 2,817,316 contigs
and 35,688,628,875 bases in 11,310,273,599 bytes (*measured*, `/tmp/if2453/r23/stored.json`: `contig_index` rows and
`contig_len`, shard sizes). `r19_fasta_scan.py scan` reads such a file as it is (gzip by its first two bytes). At the
rate card R-19 measured on BGZF FASTA wrapped at 80 (the new `LR2`, 21,254,892,153 bases in 598.9 CPU-s = 35.5 Mbp per
CPU-second; the new `NX1`, 34.3; *derived*), **a scan of `SPX1`'s stored sequences would take about 17 CPU-minutes**
(*derived*: 35,688.6 Mbp / 34.3 Mbp/s = 1,040 CPU-s), about 5 minutes of wall at 4 processes over the 40 shards, read
only; not run (the card says so). It would read the root's files and write nothing there; it is the PI's to order.
For scale, all 16 stored stages' shards hold 146,367,224,495 bases in 45,603,443,175 bytes, about 71 CPU-minutes at
the same rate (*derived*); the scans of §2 already cover every other stage through its inputs.

## 4. `GV1` re-staged without the trailing gaps of three contigs

@@SECTION4@@

## 5. What a check of the alphabet at ingest would refuse

The check considered: refuse a record whose sequence holds a byte outside the IUPAC nucleotide letters
`ACGTURYSWKMBDHVN` of either case (the scanner's set; lower case, `N` runs and `U` pass). Ingest today refuses only a
header that is not UTF-8 and a line before the first header (card R-20's reading of `src/biobox/fasta.py`, R-19 §1);
it stores every other byte as a base. Over every file cards R-19, R-20 and R-23 scanned (*measured*, their reports;
nothing under `src/` changed):

| Set | Files | Records | Would the check refuse? | Records refused | Bytes outside IUPAC | Byte values | The file now |
|---|---:|---:|---|---:|---:|---|---|
@@TABLE5@@

So **an alphabet check at ingest would refuse nothing that is an input today**: the only records it would refuse are
in three files that are no longer inputs of any stage (`LR2`'s and `NX1`'s first files, superseded by cards R-19 and
R-20, and `GV1`'s first file, superseded here), and nothing the root stores holds such a byte (§2, as far as the
inputs on the pod show; §3 for `SPX1`). In the files scanned it would have caught exactly the three staging defects
found so far, and in two of the three it would have stopped the build hours earlier than ingest's UTF-8 decode did
(`LR2`) or would (`NX1` and `GV1` would not have stopped at all). What it would cost and where it could go is the PI's
to decide; three places, from the cheapest to the most binding: a check in each staging script before it writes its
report (`r19_fasta_scan.py` is one, standard library only); `s2_check_inputs.py` reading every object (card R-21,
about 10 CPU-minutes per 20 Gbp, R-19 §5); or `biobox.fasta.iter_records` refusing at ingest, which would also guard
inputs that never pass a staging script. A check at ingest would refuse whole records; a `-` at a contig's end (as in
`GV1`) is a trimming question, not a refusal, and is best settled where the record is staged.

## 6. Pod compute used

@@SECTION6@@

## 7. Tests

- `tests/r1/test_r23_gv1_trailing_gaps.py` (15 cases, each run with the normal `PATH` and with `PATH=/usr/bin:/bin`;
  the re-stage needs no pod tool) over a fixture stage in card C-48c's layout: headers and line widths as in an
  archive, `records.tsv` one row per FASTA record in file order, `genomes.json` in `json.dumps` form without a final
  newline, a genome that repeats another's contig set (`repeated_records.tsv`), a contig two genomes share, one shared
  row. Three gapped records: `--` inside the last line, `--` as a line of its own, and the gaps split over two lines.
  Truth by hand: the new FASTA written out byte for byte, the uids computed in the test with `blake3` directly. Cases:
  the gaps go and every other byte stays (FASTA, `records.tsv`, `genomes.json`, the copied tables, objects, the
  `gap_trims.tsv` rows, the counts, the removed byte ranges, which rebuild the new file from the old); the new file
  scans clean with `r19_fasta_scan` and the old one does not; the staged files are only read and a second run is byte
  for byte the same; a new uid an earlier stage holds is reported (DEC-26); a trimmed contig equal to another contig of
  its own genome is refused, one equal to a contig of another genome is reported (DEC-27); refusals (an inner gap,
  another byte that is not a base, a `records.tsv` out of step with the FASTA, a trimmed genome in
  `repeated_records.tsv`, a trimmed contig in `shared_records.tsv`, other trimmed records than expected, existing
  targets); the command line; and the template (§4.4). **Seen failing first** (`wilkcomp02`, before the script
  existed): `FileNotFoundError: [Errno 2] No such file or directory:
  '/tmp/if2453/ci/R-23/scripts/stress/build_stage_pod/r23_gv1_trailing_gaps.py'`; then, with the script and before
  the template, the template case alone failed (`FileNotFoundError: ... stages/templates/GV1.conf`, 2 failed, 28
  passed).
- No existing test was changed: `stages/GV1.conf` (read by `tests/r1/test_stage_driver.py`) is left as it is, and the
  new pre-build template is `stages/templates/GV1.conf` (§4.4).
- @@PYTEST@@
