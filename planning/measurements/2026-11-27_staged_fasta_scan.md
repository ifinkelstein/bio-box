# Staged FASTA scanned for bytes that are not sequence, and `LR2` re-staged (card R-19)

- **Date and host:** 2026-11-27, `wilkcomp02`. Every number carries a label (AGENTS §6): *measured* (on `wilkcomp02`,
  that day, by the command or file named), *derived* (formula and inputs), or *estimate*.
- **Scripts:** `scripts/stress/build_stage_pod/r19_fasta_scan.py` (`scan`, `scan-stage`, `table`; standard library
  only), tests `tests/r1/test_r19_fasta_scan.py`; the re-stage is `restage-lr2` of
  `scripts/stress/build_stage_pod/c76_lr2_sag1_nx1_r232.py`, tests in `tests/r1/test_c76_lr2_sag1_nx1_r232.py`
  (the `test_r19_*` cases); template `scripts/stress/build_stage_pod/stages/templates/LR2.conf`.
- **Read only:** `/tmp/if2453/s2/inputs/` (the wave's stages, `GV1`, `PL3`), `/tmp/if2453/stages/`, the root
  `/tmp/if2453/r10/root`, the staged metadata of the `r26.11.5` members and their FASTA on `/stor`, card C-66's
  `/tmp/if2453/c66/` and its four archives in `/stor/scratch/Wilke/if2453/bio-box-staging/LR2/dl/`. Nothing under
  `/tmp/if2453/s2/inputs/` or `/tmp/if2453/stages/` was written; nothing under `/stor` was overwritten or deleted.
- **Written:** the scan reports under `/tmp/if2453/r19/` (`wave/scan/`, `r26_11_5/scan/`, one JSON per file and one
  `<STAGE>.stage.json` per stage); the new `LR2` FASTA `lr2_compendia_r19.fna.gz` on
  `/stor/scratch/Wilke/if2453/bio-box-staging/LR2/` beside the old one; its metadata in `/tmp/if2453/r19/LR2/`; its
  r232 status in `/tmp/if2453/r19/c76/LR2/`.

## 0. Summary

R26_SUMMARY_PLACEHOLDER

## 1. The scan

`r19_fasta_scan.py scan PATH OUT_JSON` reads one FASTA, plain or gzip (taken from the first two bytes, not the name;
concatenated gzip members are read whole), in line pieces of at most 1 MiB, and holds nothing but counters and the
current header. It splits records as `biobox.fasta.iter_records` does (a line starting with `>` opens a record,
trailing CR/LF stripped, blank lines skipped) and reports: records; bases (sequence bytes); headers that are not
UTF-8; headers that are not printable ASCII (bytes `0x20` to `0x7e` and tab); records holding a byte outside the IUPAC
nucleotide letters `ACGTURYSWKMBDHVN` (either case) and the count of such bytes; empty records; lines before the first
header; a read error (a truncated gzip, trailing bytes that are not gzip); the first record of each kind; and the
first offending record with its number, its header, the header before it, its reasons and its first bad byte and
offset. It also gives the file's sha256 (hashed as read). The report is written under a temporary name and renamed.
`scan-stage OBJECTS_JSON OUT_DIR --stage S` runs `scan` on every path of a stage's `objects.json` (or
`objects_staging.json`), each in its own process, one after another, skips a file whose report exists, and writes
`<S>.stage.json`. The runs on the pod used a small worker (`/tmp/if2453/r19/bin/r19_worker.sh`) that takes the stages
in order and claims each by an atomic `mkdir`, so that four workers share one ordered list. Python: the pod's own
`/usr/bin/python3` (3.12.3), since the script needs nothing beyond the standard library.

A "damaged" status means any count above is not zero. Of these, a header that is not UTF-8 stops ingest
(`biobox.fasta` decodes every header), and so would a read error or a line before the first header; a header that is
UTF-8 but not ASCII does not, nor does a byte outside the IUPAC letters (ingest stores the bytes; the gene callers and
`skani` meet them later: card C-81's `skani` stopped on `0xe5` in `LR2`).

## 2. The stages of the running wave that had not started

The wave `r26-11-4-wave-6` was building `PL3` (*measured*, `tail /tmp/if2453/stages/r26-11-4-wave-6.log`, 01:30 CST).
The 16 stages after it were scanned in its order from their `/tmp/if2453/s2/inputs/<STAGE>/objects.json`, read-only,
units `r19-scan-1` to `-4` (01:34:51 to 01:52:45 CST). Reported to the orchestrating session at 01:53 CST.

| Stage | Files | Records | Bases | Not UTF-8 | Not ASCII | Non-IUPAC records | Non-IUPAC bytes | Empty | Before header | Read errors | Status |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `CG1` | 1 | 120,568 | 10,731,821,100 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |
| `MGE1` | 10 | 26,032 | 209,759,220 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | reported: `iceberg3_ICE.fna` |
| `SAT1` | 10 | 46,662 | 649,174,505 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |
| `ECE1` | 8 | 15,262 | 1,717,947,087 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |
| `O1` | 2 | 15,225 | 2,326,329,045 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |
| `O2` | 2 | 2,415 | 338,041,679 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |
| `LPL1` | 1 | 626 | 4,994,414 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |
| `SYM1` | 1 | 2,742 | 229,989,343 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |
| `LRV1` | 2 | 2,559 | 105,497,712 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |
| `LGRT1` | 1 | 378 | 3,026,010 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |
| `MGE2` | 1 | 24,312 | 686,974,567 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |
| `MV1` | 1 | 1,177,430 | 39,984,922,073 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |
| `V3` | 1 | 691,723 | 26,146,357,334 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |
| `LR1` | 1 | 856,719 | 20,629,983,870 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |
| `GART1` | 1 | 5,218,293 | 39,413,154,433 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |
| `GART2` | 1 | 1,213,095 | 11,853,071,408 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |

All *measured* (`/tmp/if2453/r19/wave/scan/<STAGE>.stage.json`, summed per stage; one row per file, with its sha256,
in `r19_fasta_scan.py table /tmp/if2453/r19/wave/scan` and in each file's JSON). **No wave stage holds a header that is not
UTF-8, a byte outside the IUPAC letters, an empty record, a line before the first header or a read error.** The one
report: `MGE1`'s `iceberg3_ICE.fna` has 2 headers that are UTF-8 but not ASCII; the first is record 739,
`ICEberg|1869|ICESpsH35–1|GenBank|CP076721.1|1607033..1676868 Streptococcus parasuis strain H35 chromosome, complete
genome` (an en dash, U+2013). It does not stop ingest. `PL2` is built and `PL3` was past ingest, so neither was
scanned as a wave stage; `PL3` was scanned with the `r26.11.5` list (§3).

## 3. The staged members of `r26.11.5`

Every member that has a template under `stages/templates/` with real `SUM=` lines, scanned from its
`objects_staging.json` (the files where they are now; their `objects.json` names `/tmp/if2453/s2/inputs/<STAGE>/`,
where card A-29 will copy them), and three more: `GV1` (staged in `/tmp/if2453/s2/inputs/GV1/`, not yet built), `WET1`
(staged by card C-81, `/tmp/if2453/c81/meta/WET1/`, its template not yet merged) and `PL3` (building, past ingest).
`ET1` and `LGNRNAP1` have placeholder `SUM`s (all zeros) and no staged files, so there was nothing to scan. Units
`r19-scan-5` to `-8`, 01:46:41 to 02:07 CST.

| Stage | Files | Records | Bases | Not UTF-8 | Not ASCII | Non-IUPAC records | Non-IUPAC bytes | Empty | Before header | Read errors | Status |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `GRARE1` | 1 | 1,659,995 | 11,583,415,844 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |
| `GRARE2` | 1 | 1,253,559 | 10,455,837,641 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |
| `GRARE3` | 1 | 1,341,468 | 10,335,237,414 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |
| `LR2` (old file) | 1 | 575,827 | 21,670,443,395 | 20,526 | 20,526 | 27,396 | 363,929,092 | 130 | 0 | 0 | **damaged** (§4) |
| `SAG1` | 1 | 236,405 | 4,941,715,933 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |
| `NX1` | 1 | 391,103 | 3,398,390,142 | 0 | 0 | 69 | 17,393 | 0 | 0 | 0 | **reported**: 3 alignment files staged as genomes (below) |
| `GV1` | 1 | 288,536 | 4,576,307,473 | 0 | 0 | 3 | 6 | 0 | 0 | 0 | reported: 6 gap bytes `-` in 3 contigs |
| `EGV1` | 1 | 3,306 | 51,192,459 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |
| `ECE2` | 2 | 966 | 17,828,599 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |
| `MIR1` | 6 | 68,647 | 1,185,579,620 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |
| `MV2` | 1 | 101,573 | 3,909,111,248 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |
| `MVEUKV` | 1 | 106,087 | 1,257,239,393 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |
| `KAVHIFI` | 2 | 530,743 | 4,688,660,420 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |
| `INPH2604` | 1 | 386 | 34,889,413 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |
| `SAG2` | 1 | 620,195 | 8,020,962,436 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |
| `GWV1` | 1 | 280,420 | 3,740,796,061 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |
| `SAG3` | 1 | 301,688 | 2,408,106,251 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |
| `WET1` | 1 | 1,749,972 | 12,615,718,580 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |
| `DPN1` | 2 | 135,231 | 1,798,397,715 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |
| `SYM2` | 2 | 89,775 | 868,366,559 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |
| `LRP1` | 1 | 78,737 | 240,689,693 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |
| `MIR2` | 2 | 56,092 | 585,909,240 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |
| `PL3` (building) | 2 | 792,463 | 17,842,081,525 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | clean |

All *measured* (`/tmp/if2453/r19/r26_11_5/scan/<STAGE>.stage.json`). **Other than the old `LR2` file, no staged
input holds a header that is not UTF-8, an empty record, a line before the first header or a read error; nothing else
would stop ingest.** Two stages hold gap bytes, all of them `-` (`0x2d`): a second scan with the byte-value count added
to the script (unit `r19-rescan`, `/tmp/if2453/r19/rescan/scan/`) gives `NX1` `{0x2d: 17,393}` and `GV1`
`{0x2d: 6}`, one distinct value each.

- **`NX1`: three "genomes" are 16S alignments, not genome sequence.** The 69 records are all the records of
  `lost_city:Bipolaricaulota-16S-phylogeny` (18 records of 1,592 columns), `lost_city:Methanosarcinales-16S-phylogeny`
  (34 of 676) and `lost_city:Thermodesulfovibrionales-16S-phylogeny` (17 of 1,665) in `/tmp/if2453/c66/NX1/records.tsv`
  (*measured*, `grep -i -E 'phylogen|16S|align|tree'`: 69 rows, these three genome ids; the first offending record
  274,358, `4a1f9c8c…`, is the first of them). Card C-66 took them from the Lost City zip as MAGs, so `NX1`'s 32
  `lost_city_serpentinizing_mags` genomes are 29 MAGs and 3 alignment files; `NX1` holds 893 genomes, not 896, and
  about 75 kbp of alignment (*derived*: 18 × 1,592 + 34 × 676 + 17 × 1,665 = 79,945 columns, of which 17,393 gaps). By
  the standing rules this is a removal of non-sequence files; it is not this card's to make (`NX1`'s files are not in
  its Files): **the orchestrating session should card it.**
- `GV1`: 6 gap bytes in 3 contigs; the first is `Hoko_1` (record 13,777, 450,697 bytes) at offset 450,695, so the
  record ends in `--`. A source artefact; it does not stop ingest.

No other member holds any of these.

## 4. `LR2`

### 4.1 The old file, scanned (the positive control on real data)

`/stor/scratch/Wilke/if2453/bio-box-staging/LR2/lr2_compendia.fna.gz` (6,329,546,501 bytes, sha256
`89dba3ba4098796fd3042b132089b9ea5393949d48dc7bf2783912e0b1eef377`, the file card C-66 recorded): **damaged**, 575,827
records, 21,670,443,395 bases, **20,526 headers not UTF-8** (the same 20,526 not printable ASCII), **27,396 records
with bytes outside the IUPAC letters** (363,929,092 such bytes), 130 empty records, no read error. The first
offending record is record **504,174** (header `5ddd4555721b18a765049dbb2004e612`, first bad byte `0x46` at offset 0;
the header before it `1e5159ff6725a3d531c0de6247b703db`); the first header that is not UTF-8 is record 504,176 and the
first empty record 504,217 (*measured*, unit `r19-scan-5`, 01:46:41 to 01:56:20 CST,
`/tmp/if2453/r19/r26_11_5/scan/LR2__lr2_compendia.fna.gz.json`). The record count, the 20,526 and the 27,396 equal
card C-81's scan exactly; C-81 put the start at record 504,175, one later (its count of the record before the damage).

### 4.2 The cause: the Bickhart bins went in compressed

C-66's `LR2` was not written by `c66_rt_hotspot_compendia.py stage` but by a session script on the pod,
`/tmp/if2453/c66/run_stage_lr2.py` (20,094 bytes, 2026-11-24 08:05; not in the tree), which imports C-66's module and
streams the four archives in order: EMERGE 2019, Singleton 2021, the Bickhart bins, then the Bickhart metaFlye
assembly for the unbinned contigs. For Singleton it gunzips a member whose name ends in `.gz`
(`if m.name.endswith(".gz"): gzip.open(...)`); for EMERGE and the Bickhart bins it hands the member straight to
`c66._iter_fasta_stream`, which reads it as plain text, **drops every line before the first `>`** and decodes headers
with `errors="replace"`. The members of the bins archive are gzip files:

| Archive (`.../LR2/dl/`) | FASTA members | Their first two bytes | Read by C-66 as |
|---|---:|---|---|
| `EMERGE_MAGs_2019_long-short-hybrid.tar.gz` | 733 `.fna` | 733 `>` | plain (right) |
| `MQ_HQ_MAGs_20200502.tar.gz` | 3,733 `.fa.gz` | 3,733 gzip | gunzipped (right) |
| `hifi_das.bin3c.bins.tar.gz` | 691 `bin3c.*.contigs.fa.gz` + 4 `bin3c.*.contigs.fa` (in `b3c_hifi_dastool/flye4.das_DASTool_bins/`) | 691 gzip, 4 `>` | **plain: 691 wrong**, 4 right |
| `metaflye_hifi_asm.tar.gz` | `assembly.fasta` (3,484,359,955 bytes) | (tar member, plain) | plain (right) |

All *measured* (`/tmp/if2453/r19/bin/r19_inspect_lr2.py`, unit `r19-inspect2`, which streams each archive, reads the
first two bytes of every FASTA member and gunzips the bins; output `/tmp/if2453/r19/lr2_inputs_inspect.json`; the
member list `tar tzvf`, unit `r19-binslist`, `/tmp/if2453/r19/bins_members.txt`). So in step "C. Bickhart 2022 Rumen
HiFi bins" each gzip member became a "MAG" whose "contigs" were the stretches of compressed bytes after a `0x3e` byte
that happened to start a line; the stager wrote them, wrapped at 80, under contig-uid headers, and a wrapped line that
began with `0x3e` became a record of its own with a binary header. That is the whole damage, and it explains every
count (*derived*):

- the damage starts at record 504,174 = 504,173 contigs after EMERGE and Singleton (C-66's log
  `/tmp/if2453/c66/stage_lr2.log`: "after Singleton 2021: genomes=4466 contigs=504173") + 1;
- C-66's bins bases, 426,856,095, are compressed bytes: the 691 gzip members hold 464,705,172 bytes (the rest of each
  member before its first `>` line and the newline bytes were dropped); the real bins hold **1,505,258,010 bases in
  5,489 contigs** over 695 bins;
- C-66 counted 642 bins genomes from 695 members (members with no `>` line gave no record) and 7,031 contigs
  (511,204 − 504,173);
- 575,827 − 555,301 = 20,526 extra records, the headers that are not UTF-8;
- the unbinned set was cut against binned contig uids that were compressed bytes, so only the 4 plain bins were
  seen: "22 already-binned contigs skipped". 5,216 of the 5,489 binned contigs are of 5 kb and more, so about 5,194
  metaFlye contigs that are in bins were staged a second time as unbinned genomes.

The EMERGE, Singleton and metaFlye records are sound: the scan finds nothing before record 504,174, and the
re-stage reproduces C-66's `records.tsv` rows of EMERGE and Singleton row for row (§4.3).

**The inputs are sound.** The four archives on the pod are the files C-66 fetched: their sha256 today
(`acb77d18…`, `0837e380…`, `bb82982f…`, `8d424877…` for EMERGE, Singleton, bins, metaFlye; *measured*, `sha256sum`, unit
`r19-archsha`, `/tmp/if2453/r19/lr2_archives.sha256`) equal those `sources.yaml` recorded at fetch, when the sizes and
the providers' md5 matched. `r19-inspect2` read all four to the end
through `tarfile` and `gzip` with no error, and every bin member decompressed. On `tower` there are no copies of the
four archives (`find /mnt/user/research-bulk/bio-box -maxdepth 4` for their names finds none; the one hit,
`sources_incoming/LR_ARC/metaFlye_assembly.fasta`, is another set). `tower:/mnt/user/research-bulk/bio-box/sources/LR2/`
holds the mirror of C-66's staged files: its `SHA256SUMS` gives `lr2_compendia.fna.gz` the sha256 above, so it is the
same damaged file (not scanned again on `tower`: same bytes). It stays there (nothing on `tower` is deleted); it is
superseded (§4.4).

### 4.3 The re-stage

LR2_RESTAGE_PLACEHOLDER

## 5. Should `s2_check_inputs.py` have caught it?

No, as it is written: `s2_check_inputs.py` parses `genomes.json` and the objects list with the readers of `stage inputs`
and `stage ingest` and checks that **each object path is a file**; it never opens the FASTA (its docstring says so).
C-66's and C-76's checks of `LR2` (`problems: []`) were checks of its tables, which were consistent with each other
(the uids in `genomes.json` were the blake3 of the compressed bytes that were written). The build would have stopped
at `ingest`, at the first header that is not UTF-8 (record 504,176), hours into the stage.

A few lines after `objects = cli._source_objects(opath)` would catch it, at the cost of one read of every object
(about 10 minutes for a stage of 20 Gbp, *estimate* from §2's rates). For a card (not changed here):

```python
from biobox import fasta, ids
want = {u for g in genomes for u in g.contig_uids}
for o in objects:
    for rec in fasta.iter_records(Path(o.path)):        # raises UnicodeDecodeError on a header that is not UTF-8
        if rec.seq.translate(None, b"ACGTURYSWKMBDHVNacgturyswkmbdhvn"):
            problems.append(f"object {o.path}: record {rec.accession}: a byte that is not a nucleotide letter")
        want.discard(ids.contig_uid(rec.seq))
if want:
    problems.append(f"{len(want)} contigs of genomes.json are in no object")
```

with the `UnicodeDecodeError` caught and reported as a problem, and an option to skip it for a re-check of a stage
whose objects were already read. The uid test also catches a FASTA that does not hold the contigs `genomes.json` names.
`r19_fasta_scan.py` is the same check without `biobox` and could be run by a staging script before it writes its
report.

## 6. Pod compute used

COMPUTE_PLACEHOLDER

## 7. Tests

TESTS_PLACEHOLDER
