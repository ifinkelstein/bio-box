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

1. **The catalog holds no byte that is not a base.** The root stores 16 stages. The inputs of 15 of them are still on
   the pod, as the objects themselves (`V2`, `PL2`, `PL3`) or as archived copies (`S1`, `P1` to `P6`, `H1`, `PL1`,
   `C5`, `V1`, `G1`), each with the sha256 the stage's `l0_manifest` recorded at ingest. Scanned (`PL3` by card
   R-19): 21,153 objects, 2,867,052 records, 129.98 Gbp, **all clean** (§2). `SPX1`'s inputs are gone; its 40 stored
   shards (2,817,316 contigs, 35.69 Gbp) were scanned instead, on the orchestrating session's word, and are clean too
   (§3). Reported to the orchestrating session: an interim table of 14 clean stages before part (2) started, and the
   full table (with `P6` and `C5`) while the re-stage ran.
2. **`GV1` re-staged** (§4): the three `--` removed (6 bytes), every other byte of the FASTA the old file's (the
   script's read-back and an independent `cmp`); the new FASTA scans clean. Without its gaps **`Hoko` is the GenBank
   assembly `GCA_002966255-1_ASM296625v1` that `GV1` already holds** (Hokovirus, KY684103 to KY684107, all five
   contigs), so by C-48c's DEC-27 rule it becomes that genome's repeat: **18,498 genomes, 281,801 distinct contigs,
   4,315,082,073 bp** (old 18,499, 281,804, 4,316,061,173). No other stage holds the trimmed contigs; `MIR1`'s GOEV
   Hokovirus reference holds them and already counts them as `GV1`'s. Ownership is unchanged (147,791 contigs held),
   134,010 contigs and 1.640 Gbp new. `s2_check_inputs.py` exit 0, `problems: []`; `stage_est.py conf` exit 0;
   template `stages/templates/GV1.conf` with real `SUM`s; mirrored to `tower` under `sources/GV1/r23_restage/`.
3. **An alphabet check at ingest** (§5) would refuse records in exactly three files, all now superseded (`LR2`'s,
   `NX1`'s and `GV1`'s first files); nothing that is an input today, and nothing in the root.

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
| `C5` | 2026-11-19 17:24 | `/tmp/if2453/s2/inputs/C5/objects.json` | no | 2,032 | 2,781,209,046 | members of `provisional/s2c_phase2/s2c_phase2_genomes.tar` and `s2c_phase2_rest_genomes.tar` |
| `V1` | 2026-11-20 15:28 | `/tmp/if2453/s2/inputs/V1/objects.json` | no | 1 | 5,228,187,536 | `/tmp/if2453/s3/dl/uhgv/votus_full.fna.zst` (the UHGV download), decompressed |
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
| `SPX1` (inputs gone; its 40 stored shards, §3) | 40 | 2,817,316 | 35,688,628,875 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | none | 40 of 40 (= run records) | clean |

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
(*derived*: 35,688.6 Mbp / 34.3 Mbp/s = 1,040 CPU-s). For scale, all 16 stored stages' shards hold 146,367,224,495
bases in 45,603,443,175 bytes, about 71 CPU-minutes at the same rate (*derived*); the scans of §2 cover every other
stage through its inputs.

**The scan was then run**, on the orchestrating session's word after the first report (the card said not to; a dated
line goes into the card at merge): `r19_fasta_scan.py scan-stage` over a list of the 40 shard paths
(`/tmp/if2453/r23/objs_SPX1_shards.json`), unit `r23-scan-spx1`, one process, the pod's Python, the root only read,
05:50:21 to 06:07:49 CST, **17.4 CPU-minutes** (*measured*, `cpu_s` 1,043.9; the estimate above was 17.3, within 1%).

| | Shards | Records | Bases | Not UTF-8 | Not ASCII | Non-IUPAC records | Non-IUPAC bytes | Empty | Before header | Read errors | Status |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `SPX1`, stored | 40 | 2,817,316 | 35,688,628,875 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |

All *measured* (`/tmp/if2453/r23/scan/SPX1__shard-00NN.fa.gz.json`). Every shard's sha256 equals the one its ingest run
record gives (`runs/stage=SPX1/ingest.batch-00NN.json`, 40 of 40), and the records and bases equal `contig_index`'s
2,817,316 rows and the sum of their `contig_len` (`r23_shard_check.py`). **`SPX1`'s stored sequences hold no byte
outside the IUPAC letters, no empty sequence and no header that is not ASCII.**

What this can and cannot show. It scans what was **stored**, not what was read: the bytes of the 2,817,316 contigs
ingest kept (those of 5,000 bp and more of the 304 assemblies, card C-48e part 1), under the store's own headers (the
contig uids), as `SeqStore` serves them. So it settles the question that matters for the catalog (no stored base of
`SPX1` is a byte that is not a base, and none would reach a gene caller), and, since ingest stores a record's bytes as
read, it says the same of every record ingest kept. It cannot show anything about what ingest did not keep: the
contigs under 5,000 bp, the records of a contig another stage owns (stored there, and clean there by §2), and the
inputs' own headers (a header that is not UTF-8 would have stopped ingest, so none was there; one that is UTF-8 but not
ASCII would have passed and is not visible here). Nor can it tell the input files' sha256, since they are gone.

## 4. `GV1` re-staged without the trailing gaps of three contigs

### 4.1 The re-stage

Card R-19 found the last two bytes of `GV1`'s `Hoko_1`, `Hoko_3` and `Hoko_4` to be `--`, and card R-20 (its note §8.3)
that ingest would store them as bases and hash them into the uids. `r23_gv1_trailing_gaps.py restage` reads card C-48c's
staged `GV1` (`/tmp/if2453/s2/inputs/GV1/`, bound read-only) and writes it again under a new FASTA name (module
docstring): each record's bytes are checked against its `records.tsv` row (ordinal, accession, `contig_uid`, length:
all 288,536 agree), a record is changed only when its bytes outside IUPAC are a run of `-` at its end (exactly the
three `--expect`ed records), and the new file is read back beside the old one. Run (unit `r23-gv1`, in
`/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif`, the tree of this branch on `PYTHONPATH`, `CPUQuota=200%`,
16G, the root and `/tmp/if2453/s2/inputs` bound read-only):

```
r23_gv1_trailing_gaps.py restage /tmp/if2453/s2/inputs/GV1 /tmp/if2453/r23/GV1 \
    /stor/scratch/Wilke/if2453/bio-box-staging/GV1/gvmags_v2_all_18kgenomes_r23.fna /tmp/if2453/s2/inputs/GV1 \
    --root /tmp/if2453/r10/root --chain-home /tmp/if2453/s2/inputs --work /tmp/if2453/r23/work_gv1
```

05:54:59 to 05:58:44 CST, exit 0 (*measured*, `/tmp/if2453/r23/logs/`; `stage_report.json` in `/tmp/if2453/r23/GV1/`).
A first run at 05:51 **refused** before writing anything (its temporary FASTA removed): "a genome's contig set became
another genome's (DEC-27)". That is the finding of §4.2; the script now applies C-48c's own rule to it, with a test.

| Record (ordinal) | Old `contig_uid`, length | Removed | New `contig_uid`, length | Offset of the removed bytes in the old file |
|---|---|---|---|---:|
| `Hoko_1` (13,776) | `800cb43e5ffc2b514515ee1703e528f1`, 450,697 | `--` | `346a88ae054c6b00c7d21d8ac9f03fe0`, 450,695 | 493,215,808 |
| `Hoko_3` (13,778) | `8bed7beb341979682ead52146fd853b9`, 407,935 | `--` | `256656ff45e6efd70d853a0ed2ec6f41`, 407,933 | 493,941,039 |
| `Hoko_4` (13,779) | `a88c0ca1c4292085b74d4c34edf09683`, 120,468 | `--` | `6267248e2a5ff1900f1962916799bfe5`, 120,466 | 494,061,516 |

(*measured*, `gap_trims.tsv`). The new uids are those card R-20 computed (its §8.3), and in each record the two bytes
were the end of its last line, so the file loses exactly 6 bytes.

### 4.2 Without its gaps, `Hoko` is a genome `GV1` already holds

The three new uids are not new to `GV1`: they are contigs 0, 2 and 3 of the genome
`GCA_002966255-1_ASM296625v1` (`confirmed/GCA_002966255-1_ASM296625v1.fna`, records 2,659 to 2,663, headers
`…|KY684103__1` to `…|KY684107__1`), and `Hoko`'s two other contigs, `Hoko_2` (`1a5927a7…`, 317,278) and `Hoko_5`
(`1e99521e…`, 29,857), are that genome's contigs 1 and 4 already (*measured*, `grep` of the staged `records.tsv`).
`Hoko` is Hokovirus (GenBank KY684103 to KY684107, the assembly `GCA_002966255.1`); GVMAGs V2 carries it twice, once
from GenBank and once under its name with `--` padding at the ends of three contigs. With the padding gone the two
genomes have the same contig set, and C-48c's DEC-27 rule (a genome whose contig set an earlier genome has is that
genome) makes `Hoko`, the later one in the archive, a repeat of `GCA_002966255-1_ASM296625v1`: it leaves
`genomes.json` and gets the row `13776  Hoko  4b9ec1946d3ff9c1d981fee1b8fec082  GCA_002966255-1_ASM296625v1` in
`repeated_records.tsv`, between the old rows, every old row as it was; its five records stay in the FASTA and in
`records.tsv`, as the 227 repeats C-48c found stay (they are ingested as occurrences). DEC-26 asks whether the
trimmed contigs equal contigs some stage already holds: **no other stage** holds them (`trimmed_uids_held: {}` over the
root's 13 stored stages and the 19 of `r26.11.4`, and no `records.tsv`, `genomes.json` or `shared_records.tsv` of any
staged stage on the pod names them, §4.3); only `GV1` itself does.

### 4.3 The new `GV1` against the old

| | old (C-48c, `gvmags_v2_all_18kgenomes.fna`) | new (R-23, `gvmags_v2_all_18kgenomes_r23.fna`) |
|---|---:|---:|
| Genomes (`genomes.json`) | 18,499 | **18,498** (`Hoko` a repeat, §4.2) |
| Repeated genomes (`repeated_records.tsv`) | 227 | 228 |
| Memberships (`genome_input` rows) | 287,159 | 287,154 |
| FASTA records (`records.tsv` rows) | 288,536 | 288,536 |
| Bases in the FASTA | 4,576,307,473 | 4,576,307,467 |
| Distinct contigs | 281,804 | **281,801** |
| Distinct bases | 4,316,061,173 | **4,315,082,073** |
| Held by the root's 13 stages and the `r26.11.4` chain | 147,791 contigs, 2,674,672,049 bp (H1 6, SPX1 3, V2 17,649, PL3 2, MV1 130,131) | the same |
| New | 134,013 contigs, 1,641,389,124 bp | **134,010 contigs, 1,640,410,024 bp** |
| Bytes outside IUPAC | 6 (`0x2d`) in 3 records | **0** (§4.4) |
| FASTA bytes | 4,639,876,745 | 4,639,876,739 |
| FASTA sha256 | `568f3a312fd78d03544792bfb3a89bde708591865065d72c62fc6c0dbdc0028b` | `a87aad1f9987179d55c6d6134f60a3ce7f0bf63b7b8975c065982dd40b920bfb` |

All *measured* (`stage_report.json`: `old`, `new`, `ownership_old`, `ownership_new`). The differences, *derived*: 6
bytes; 3 distinct contigs (the gapped ones) and 450,697 + 407,935 + 120,468 = 979,100 distinct bases (4,316,061,173 −
979,100 = 4,315,082,073); one genome and its 5 memberships. The ownership counts are the template's `EXPECT` lines
(H1 6, SPX1 3, V2 17649, PL3 2, MV1 130131; *measured*) and equal `stages/GV1.conf`'s (C-53 and A-26). The `held`
column is over the stages built before `GV1` in its build order; the old note's 3.89 Gbp "new after `V2`" (C-48c)
counted only `V2` and the stages built then, and `MV1`'s 130,131 contigs came later (A-26).

Files (*measured*, sha256): `/stor/scratch/Wilke/if2453/bio-box-staging/GV1/gvmags_v2_all_18kgenomes_r23.fna`
4,639,876,739 B `a87aad1f9987179d55c6d6134f60a3ce7f0bf63b7b8975c065982dd40b920bfb`; in `/tmp/if2453/r23/GV1/`:
`genomes.json` `8b2f106b6727a0853d3f0c0b496959d2fd34ecd813626f537830a674011c223d`, `objects.json`
`fbfeb221dc25cc5e6cdcc1a4f4e77bf64838dd8d2b6c341ce03373e350f0ca54` (path
`/tmp/if2453/s2/inputs/GV1/gvmags_v2_all_18kgenomes_r23.fna`; source, release, licence, `redistributable: true` and
`retrieved_at` as C-48c wrote them: the source release is the provider's and did not change), `objects_staging.json`
`e97d7872410303f000e7831bf3ffe7b47e51e98c69fd56dbc029f20f6fcf159f`, `records.tsv`
`816960d618dd66738b5eac6ce714f02a14fbbe6e599ed83022ac6340240f1001`, `repeated_records.tsv`
`98f014a248fd4c9baf2bbd13c850b61cb6ed3fcf25009852a8216bd3cdfc997e`, `gap_trims.tsv`
`ca455f3117e40a4539c2df3e75726b70c3b265a59c5f62d18aca2f817683e748`, and `shared_records.tsv` (`4fb544d5…`) and
`skipped_members.tsv` (`2679a24e…`) byte-identical to C-48c's. The new FASTA is plain, as the old one is.

### 4.4 Checks

All in unit `r23-gv1-check`, 06:00 to 06:02 CST, read only (*measured*, `/tmp/if2453/r23/newgv1/`):

- **The new FASTA scans clean** (`r19_fasta_scan.py scan`, the pod's Python): 288,536 records, 4,576,307,467 bases
  (both the stage's), 0 headers not UTF-8 or not ASCII, 0 bytes outside IUPAC, 0 empty records, no read error, sha256
  `a87aad1f…` as written (`scan/GV1__gvmags_v2_all_18kgenomes_r23.fna.json`).
- **Every other byte is the old file's**, twice over. The re-stage's own read-back: all 288,536 headers equal, the
  288,533 untouched records byte for byte, the three trimmed records the old bytes less their `--`, every record's
  `contig_uid` and length equal to the new `records.tsv` (`other_records_byte_identical: true`, `uid_mismatches: 0`).
  And an oracle outside the script, GNU `cmp` over the four stretches between the removed bytes (`r23_cmp.py`:
  `cmp -s -n LEN -i OLD:NEW`): 493,215,808 + 725,229 + 120,475 + 4,145,815,227 bytes, each exit 0; the new file is the
  old one less exactly the 6 bytes at 493,215,808, 493,941,039 and 494,061,516, each read back as `--` (`cmp.json`).
- **No staged stage names the three new uids other than `GV1` and `MIR1`** (`r23_uid_grep.py`, `grep -F` over every
  `records.tsv`, `genomes.json` and `shared_records.tsv` under `/tmp/if2453`, CI copies and the root's parquet left
  out; `uid_grep.json`): C-48c's `GV1` (its `GCA_002966255-1_ASM296625v1`, and the old uids in `Hoko`), and card C-70's
  `MIR1` (`/tmp/if2453/c70/meta/MIR1/`), whose GOEV reference `Mimiviridae_KY684103_Hokovirus_HKV1_Hokovirus_1` holds the
  same three contigs and already records them as owned by `GV1` (`owner_stage` `GV1` in its `records.tsv` and
  `shared_records.tsv`: they were `GV1`'s through the GenBank genome before this card). So `MIR1`'s staging needs no
  change. The old uids are in no stage but `GV1`.
- **`s2_check_inputs.py GV1 genomes.json objects_staging.json`** (in the image, this branch's tree): exit 0,
  `problems: []`, 18,498 genomes, 287,154 `genome_input` rows, 18,498 distinct genome uids, 0 circular rows, kind
  `virus` (gv-meta), code 11, quality `contig`.
- **`stage_est.py conf stages/templates/GV1.conf --stage GV1 --tree <tree>`** (the pod's Python): exit 0,
  `DEPTH=annotated`, `LIBNAME=pfam LIBVER=38.2 LIBCUT=gathering`, `CRISPR=yes`, `AFTER=V2 PL3 GART2`, the six notes.
  The template's `SUM` lines are the sha256 of the new `genomes.json`, `objects.json` and
  `gvmags_v2_all_18kgenomes_r23.fna` above; its `EXPECT` lines are §4.3's ownership; `JUDGE=genecall 236` is kept, as `stages/GV1.conf`
  has it. The build's cost is not estimated again here: the re-stage removes 979,100 of 1,641,389,124 new bases
  (0.06%, *derived*), well inside AGENTS §6's 20%.
- **Mirror to `tower`** with `scripts/orchestrator/relay_to_tower.sh`, into the new directory
  `tower:/mnt/user/research-bulk/bio-box/sources/GV1/r23_restage/` (nothing else written on `tower`): the nine metadata
  files (53,578,693 bytes, `DONE: 9 verified, 0 in place and skipped, 0 failed`, 06:03:18 CST) and the FASTA
  (4,639,876,739 bytes in 183 s, 25.35 MB/s, sha256 equal to the record, `DONE: 1 verified, 0 in place and skipped,
  0 failed`, 06:06:31 CST) (*measured*, the script's log); `SHA256SUMS` there has 10 lines, every one verified in these
  runs.

**Superseded, not deleted:** `/tmp/if2453/s2/inputs/GV1/` (C-48c's `gvmags_v2_all_18kgenomes.fna`, sha256
`568f3a31…`, and its tables) is superseded by the files above and is not to be built; it was not written to. Card
A-29 (or the card that promotes `GV1`) moves the new files into `/tmp/if2453/s2/inputs/GV1/` (the new FASTA has a
new name; `genomes.json`, `objects.json`, `records.tsv` and `repeated_records.tsv` replace C-48c's, whose fate is the
PI's) and replaces `stages/GV1.conf` with `stages/templates/GV1.conf`. Until then `stages/GV1.conf` still names the old
files, so `stage_driver.sh GV1` would build the old `GV1`: its `SUM` check passes on the old files.

## 5. What a check of the alphabet at ingest would refuse

The check considered: refuse a record whose sequence holds a byte outside the IUPAC nucleotide letters
`ACGTURYSWKMBDHVN` of either case (the scanner's set; lower case, `N` runs and `U` pass). Ingest today refuses only a
header that is not UTF-8 and a line before the first header (card R-20's reading of `src/biobox/fasta.py`, R-19 §1);
it stores every other byte as a base. Over every file cards R-19, R-20 and R-23 scanned (*measured*, their reports;
nothing under `src/` changed):

| Set | Files | Records | Would the check refuse? | Records refused | Bytes outside IUPAC | Byte values | The file now |
|---|---:|---:|---|---:|---:|---|---|
| R-19 §2: the 16 wave stages not yet started (`CG1` … `GART2`) | 44 | 9,414,041 | no | 0 | 0 | none | inputs of the wave |
| R-19 §3: the staged `r26.11.5` members and `PL3`, less the three below | 31 | 9,407,708 | no | 0 | 0 | none | inputs (`PL3` built) |
| `LR2`'s first file (C-66, `lr2_compendia.fna.gz`) | 1 | 575,827 | **yes** | 27,396 | 363,929,092 | compressed data, many values (R-19 §4) | superseded by card R-19 |
| `NX1`'s first file (C-66, `nx1_compendia.fna.gz`) | 1 | 391,103 | **yes** | 69 | 17,393 | `0x2d` only | superseded by card R-20 |
| `GV1`'s first file (C-48c, `gvmags_v2_all_18kgenomes.fna`) | 1 | 288,536 | **yes** | 3 | 6 | `0x2d` only | superseded here (§4) |
| R-19: the new `LR2` (`lr2_compendia_r19.fna.gz`) | 1 | 548,565 | no | 0 | 0 | none | input |
| R-20: the new `NX1` (`nx1_compendia_r20.fna.gz`) | 1 | 391,034 | no | 0 | 0 | none | input |
| R-23 §2: the inputs of 14 built stages (`S1` … `PL2`; `PL3` is in R-19's row) | 21,151 | 2,074,589 | no | 0 | 0 | none | built |
| R-23 §3: `SPX1`'s 40 stored shards | 40 | 2,817,316 | no | 0 | 0 | none | built (what was stored) |
| R-23 §4: the new `GV1` (`gvmags_v2_all_18kgenomes_r23.fna`) | 1 | 288,536 | no | 0 | 0 | none | input |

Records summed from the cards' tables (*derived*); each row's counts *measured* by its card. `LR2`'s first file
also has 20,526 headers that are not UTF-8, which ingest already refuses, and 130 empty records, which an alphabet check
would not catch (an empty record has no byte to refuse; a separate rule would). `MGE1`'s two UTF-8 headers that are not
ASCII are headers, not sequence, and pass either way.

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

Detached `systemd-run --user` units at `nice -n 19`, `MemorySwapMax=0`, `IOSchedulingClass=idle`, each started after
a look at the active units (the wave was in `PL3`'s Pfam scan throughout; no `MV1` or `V3` unit ran; card R-21's
`r21-newlr2` and `r21-gwv1` ran for a while, and then this card kept to two units) and with no `STOP` file in the
root. At most four processes of this card at a time, one CPU each except the re-stage. Times CST (the reports are UTC).

| Unit(s) | What | Cap | Wall | CPU |
|---|---|---|---|---|
| `r23-find` | `find` by name over `/tmp/if2453` and the project's `/stor` areas (§1.1); stopped once the archived copies were found | `CPUQuota=100%` | about 20 min, from 05:11 | at most 20 CPU-minutes (*derived*: wall × 1; mostly waiting on the file system) |
| `r23-scan-1`, `-2`, `-4`, `-6`, `-7` | `V2`, `PL2`, `PL1` `H1` `S1` `P1` `P2`, `V1`, `C5` | `CPUQuota=100%`, 8G each | 05:17:42 to 05:54:47 | 25.9 CPU-minutes in the scanner (*measured*: sum of `cpu_s`, 1,553.5 s) |
| `r23-scan-3` | `G1`, 7,139 files, one process each | same | 05:29:11 to 05:48:09 | 2.0 CPU-minutes in the scanner (*measured*, 118.0 s) and about 4 more in 7,139 Python start-ups (*estimate*) |
| `r23-scan-5` | `P3` to `P6` (its `V1` step stopped at once, §2; stopped by hand at 05:54 before it read `C5` a second time) | same | 05:33:38 to 05:54:07 | 15.4 CPU-minutes (*measured*, 925.4 s) |
| `r23-scan-spx1` | `SPX1`'s 40 stored shards (§3) | same | 05:50:21 to 06:07:49 | 17.4 CPU-minutes (*measured*, 1,043.9 s) |
| `r23-gv1` | the re-stage (§4.1); a first start at 05:51 refused within a minute | `CPUQuota=200%`, 16G | 05:54:59 to 05:58:44 | at most 7.5 CPU-minutes (*derived*: 3.75 min × 2), and under 2 for the refused start (*derived*) |
| `r23-gv1-check` | §4.4 | `CPUQuota=100%`, 8G | 06:00:02 to about 06:02:30 | 1.5 CPU-minutes in the scan (*measured*, 88.0 s), at most 2.5 in all (*derived*) |

Run directly, each a few seconds: the reads of the root's run records, `l0_manifest` and `contig_index` (the last two
in the image), the `r23_match.py` checks, the stream-scan self-tests, and `scripts/pod_pytest.sh` three times (one
process each). In all at most about 1.6 CPU-hours (*derived*: 20 + 25.9 + 2.0 + 4 + 15.4 + 17.4 + 7.5 + 2 + 2.5,
the upper bounds summed); no process came near 300 CPU-minutes. Nothing was written to the root, to
`/tmp/if2453/s2/inputs/`, to `/tmp/if2453/stages/` or over anything on `/stor`; the one new directory on `/stor` is
`/stor/scratch/Wilke/if2453/bio-box-staging/GV1/` with the new FASTA.

## 7. Tests

- `tests/r1/test_r23_gv1_trailing_gaps.py` (16 cases, each run with the normal `PATH` and with `PATH=/usr/bin:/bin`;
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
  passed). The sixteenth case, a trimmed genome equal to an earlier genome becoming its repeat (`Hoko`'s case, §4.2),
  was written after the first run on `GV1` refused with "a genome's contig set became another genome's (DEC-27)" (the
  real data failing the script); with the rule in the script it passes (30 passed, 2 failed: the template).
- No existing test was changed: `stages/GV1.conf` (read by `tests/r1/test_stage_driver.py`) is left as it is, and the
  new pre-build template is `stages/templates/GV1.conf` (§4.4).
- `env BIOBOX_TREE=/tmp/biobox-wt/R-23 sh scripts/pod_pytest.sh R-23 tests/r1/test_r23_gv1_trailing_gaps.py
  tests/r1/test_r19_fasta_scan.py tests/r0`: **290 passed in 10.87 s** (*measured*, `wilkcomp02`, 2026-11-27, one
  process).
