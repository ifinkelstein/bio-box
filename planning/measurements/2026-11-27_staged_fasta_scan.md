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
  r232 status in `/tmp/if2453/r19/c76/LR2/`; its scan in `/tmp/if2453/r19/newlr2/scan/`; and the mirror of these on
  `tower` in the new directory `sources/LR2/r19_restage/` (nothing else written on `tower`).

## 0. Summary

1. **The wave's 16 unstarted stages (44 files) would all pass ingest** (§2): no header that is not UTF-8, no byte
   outside the IUPAC letters, no empty record, no read error. `MGE1`'s `iceberg3_ICE.fna` has 2 UTF-8 headers that are
   not ASCII (an en dash), which ingest accepts. Reported to the orchestrating session at 01:53 CST.
2. **Of the staged `r26.11.5` members (and `GV1`, `WET1`, `PL3`), only the old `LR2` file is damaged** (§3): 575,827
   records, 20,526 headers not UTF-8, 27,396 records with bytes outside IUPAC, from record 504,174, equal to card
   C-81's counts. Two more are reported: **`NX1` holds three 16S alignment files as genomes** (69 records, 17,393 gap
   bytes `-`; a staging defect for a card) and `GV1` 6 gap bytes in 3 contigs.
3. **Cause** (§4.2): card C-66's stager read the 691 gzip-compressed Bickhart bins (`bin3c.*.contigs.fa.gz`) as plain
   text, so their compressed bytes went into the FASTA as sequence, and the unbinned set took the binned contigs a
   second time. The four archives on the pod are sound (sha256 as at fetch; read to the end); none is on `tower`.
4. **The new `LR2`** (§4.3, §4.4): `lr2_compendia_r19.fna.gz`, 44,064 genomes (733 + 3,733 + 695 bins + 38,903
   unbinned), 548,565 contigs, 21.255 Gbp (old: 49,205, 555,301, 21.672), 21.132 Gbp new; it scans clean;
   `s2_check_inputs.py` exit 0 with `problems: []`; every genome's r232 status recorded (all not visible), none
   removed; template with real `SUM`s; mirrored to `tower` under `sources/LR2/r19_restage/`.
5. **`s2_check_inputs.py` could not have caught it** (§5): it never opens an object; a few lines that read each object
   with `biobox.fasta` and match its contig uids to `genomes.json` would.

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

`c76_lr2_sag1_nx1_r232.py restage-lr2` (card C-76's script, which the card names for the re-staging) follows C-66's
stager step for step (the same four archives in the same order, the same genome ids, quality labels, DTR rule
(DEC-65), DEC-27 checks, 5 kb cut of the unbinned contigs, 80-column wrap, contig-uid headers and `bgzip`), with
these changes:

1. a tar member is gunzipped when its first two bytes are the gzip magic, whatever its name (the fix);
2. every record is refused unless its header is UTF-8, no line precedes the first header and its sequence holds IUPAC
   letters only, so a damaged input stops the re-stage instead of entering the FASTA (none was refused);
3. a genome of EMERGE or Singleton keeps the genetic code of C-66's screen (its `genetic_codes.tsv` row), and the
   bins are screened as C-66 screened them (`screen_genetic_code`, 250 kb); the script compares its EMERGE and
   Singleton `records.tsv` rows with C-66's and reports the first difference;
4. a metaFlye contig whose name is in a Bickhart bin is not taken as unbinned even when its uid differs from the bin's
   (circular in `assembly_info.txt`, its terminal repeat trimmed in the unbinned pass and not in the bin, which C-66's
   MAG rule leaves untrimmed); and a DTR trim is recorded only for a contig that is staged.

Run (unit `r19-lr2-restage`, in `/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif`, the root, `s2/inputs`, the
archives and `/tmp/if2453/c66` bound read-only, `CPUQuota=200%`, one screening process and `bgzip -@ 1`):

```
c76_lr2_sag1_nx1_r232.py restage-lr2 /stor/scratch/Wilke/if2453/bio-box-staging/LR2/dl \
    /tmp/if2453/c66/dl/metadata_MAGs_2019_EMERGE.tsv /tmp/if2453/c66/LR2 /tmp/if2453/r19/LR2 \
    /stor/scratch/Wilke/if2453/bio-box-staging/LR2/lr2_compendia_r19.fna.gz /tmp/if2453/s2/inputs/LR2 \
    --workers 1 --compressor bgzip --bgzip-threads 1 --root /tmp/if2453/r10/root --chain-home /tmp/if2453/s2/inputs
```

02:04 to 03:28:40 CST (*measured*, file times); `bgzip` on one thread was the bound.

**The new `LR2` against the old** (*measured*, `/tmp/if2453/r19/LR2/stage_report.json` and C-66's
`/tmp/if2453/c66/LR2/stage_report.json`):

| | old (`lr2_compendia.fna.gz`, C-66) | new (`lr2_compendia_r19.fna.gz`, R-19) |
|---|---:|---:|
| Genomes | 49,205 | **44,064** |
| EMERGE 2019 / Singleton 2021 MAGs | 733 / 3,733 | 733 / 3,733 (rows identical, 504,801 compared) |
| Bickhart bins | 642 "genomes" of compressed bytes, 426,856,095 "bases" | **695 bins, 1,505,258,010 bp** |
| Bickhart unbinned contigs (`>= 5 kb`) | 44,097, 3,377,465,007 bp | **38,903, 1,881,851,628 bp** |
| Unbinned skipped as binned | 22 | 5,215 by uid + 1 by name = 5,216 (the binned contigs of 5 kb and more) |
| Unbinned under 5 kb dropped | 15,931 (40,088,657 bp) | 15,931 (40,088,657 bp) |
| Distinct contigs | 555,301 | **548,565** |
| Bases | 21,672,103,617 | **21,254,892,153** |
| Circular contigs / DTR trims | 1,863 / 20 | 1,678 / 19 |
| Genetic code 11 / 25 | 49,170 / 35 | 44,029 / 35 |
| Held by the root and the `r26.11.4` stages | 59 contigs, 123,240,285 bp | 59 contigs, 123,240,285 bp (the same `shared_records.tsv`, sha256 `b74015b3…`) |
| New | 555,242 contigs, 21,548,863,332 bp | **548,506 contigs, 21,131,651,868 bp** |
| FASTA bytes | 6,329,546,501 | 5,906,399,819 |
| FASTA sha256 | `89dba3ba4098796fd3042b132089b9ea5393949d48dc7bf2783912e0b1eef377` | `36bb8a4e47213d5bbc44471fa44a7709c483cf1c3a7cff9c946d1fef9ea14ebe` |

The counts agree with what the inspection predicted (*derived*): 4,466 + 695 + (60,050 − 15,931 − 5,216) = 44,064
genomes; 504,173 + 5,489 + 38,903 = 548,565 contigs. The bases fall by 417,211,464 (*derived*): the 426,856,095
compressed bytes go, the bins' 1,505,258,010 bases come in, and the unbinned set loses 1,495,613,379 bp to the bins.
**No genome was removed for r232 membership.** The 5,141 fewer genomes (*derived*) are the 5,194 binned contigs C-66
had also staged as unbinned genomes, less the 53 bins C-66's misreading had dropped.

Files (*measured*, sha256): `/stor/scratch/Wilke/if2453/bio-box-staging/LR2/lr2_compendia_r19.fna.gz` 5,906,399,819 B
`36bb8a4e47213d5bbc44471fa44a7709c483cf1c3a7cff9c946d1fef9ea14ebe`; in `/tmp/if2453/r19/LR2/`: `genomes.json`
`423a8280574ce5fe3dc87e42a3401cd80bbc4977c6eb904c22e1408d6ddce281`, `objects.json`
`14c5565c05659dc6636f834c699dd7f1fff6d45b7de45b0d978d20d59a317985` (path `/tmp/if2453/s2/inputs/LR2/lr2_compendia_r19.fna.gz`),
`objects_staging.json` `71eb6f3900838a1ca5cebadd39de5103f7115e449274098e652cf8c02ce9a2d1`, `records.tsv`
`1a40e438fc0e2e2927d87eb454c8ea1a558169b05a1600970938fdda1d3e51e9`, `genetic_codes.tsv`
`128a1071a718abb318674e5856db185ccf1d34bd07bde0d39e2614275e51af6d`, `trims.tsv`
`f5953424954cf32babd0ada1aa86abdaa8b3a55383d112235289dc8e51396352`, and `repeated_records.tsv`
`06c2a972da2b53aaade81df0f723c984576ce440a64838ac5172e05ef764c432` and `shared_records.tsv`
`b74015b35539dac07cbd7c3688e8d8be702b3b6042cf08c17b538922f68d6e58` (both equal to C-66's). C-66's
`EMERGE_v2_novel_for_LR3.tsv.gz` (the `LR3` list) is not a stage input and stays in `/tmp/if2453/c66/LR2/`.

### 4.4 Checks of the new `LR2`, its r232 status, the template and the mirror

- **The new FASTA scans clean:** 548,565 records (equal to the distinct contigs), 21,254,892,153 bases (equal to the
  stage's), 0 headers not UTF-8 or not ASCII, 0 bytes outside IUPAC, 0 empty records, no read error, sha256
  `36bb8a4e…` as written (*measured*, unit `r19-scan-newlr2`, 03:30:19 to 03:40:18,
  `/tmp/if2453/r19/newlr2/scan/LR2__lr2_compendia_r19.fna.gz.json`).
- **The r232 status of card C-76, kept:** `c76_lr2_sag1_nx1_r232.py run LR2 /tmp/if2453/r19/LR2
  .../lr2_compendia_r19.fna.gz /tmp/if2453/c74/index /tmp/if2453/r19/c76 --root /tmp/if2453/r10/root --chain-home
  /tmp/if2453/s2/inputs` (unit `r19-lr2-r232`, in the image) wrote `/tmp/if2453/r19/c76/LR2/r232_status.tsv`
  (sha256 `d204e98be97b893692ce81ea9bad4d4bffd6953276a0be7cde6fae3eb4e8f34c`) and `r232_report.json`
  (`803aeb252747266231843df2d55cdccfaf491aec699e719b147860902ad8cd4a`): 44,064 genomes, none with an INSDC accession,
  all `not_visible` (695 bins, 38,903 unbinned, 733 EMERGE, 3,733 Singleton), removed 0, as C-76 found for the old
  stage. Ownership over the root's 13 stored stages and the 19 `r26.11.4` stages: 59 contigs held (P1 2, P2 6, P3 44,
  PL3 3, MV1 1, GART1 3; 123,240,285 bp), 548,506 contigs and 21,131,651,868 bp new (the template's `EXPECT` lines).
- **`s2_check_inputs.py LR2 genomes.json objects_staging.json`** (run by the same command): exit 0, `problems: []`,
  44,064 genomes, 549,193 `genome_input` rows.
- **`stage_est.py conf .../stages/templates/LR2.conf --stage LR2 --tree <tree>`** (in the image, unit `r19-stage-est`):
  exit 0, `DEPTH=called`, `GENECALL_MODE=gv-meta-c25`, `CRISPR=yes`, `AFTER=GRARE3`, both notes. The template's
  `SUM` lines are the sha256 of `genomes.json`, `objects.json` and `lr2_compendia_r19.fna.gz` above; the old file's
  `SUM` is gone. Cost by plan §3's rate R2 (93.5 to 135.3 CPU-minutes per new Gbp, as C-76 used): about 1,976 to 2,859
  CPU-minutes (*derived*: × 21.132), 2% under C-76's figure.
- **Mirror to `tower`** with `scripts/orchestrator/relay_to_tower.sh`, into a new directory, because
  `sources/LR2/` holds C-66's files under the same names and nothing on `tower` is overwritten:
  `tower:/mnt/user/research-bulk/bio-box/sources/LR2/r19_restage/`. The nine metadata files (128,507,182 bytes,
  `DONE: 9 verified, 0 in place and skipped, 0 failed`), the two r232 reports (`DONE: 2 verified`) and the FASTA
  (5,906,399,819 bytes in 325 s, 18.17 MB/s, sha256 equal to the record, `DONE: 1 verified, 0 in place and skipped,
  0 failed`), 03:41 to 03:46:53 CST (*measured*, the script's log); `SHA256SUMS` there has 12 lines, every one
  verified in these runs.

**Superseded, not deleted:** `/stor/scratch/Wilke/if2453/bio-box-staging/LR2/lr2_compendia.fna.gz` (sha256
`89dba3ba…`), C-66's metadata in `/tmp/if2453/c66/LR2/`, C-76's status of it in `/tmp/if2453/c76/LR2/`, and their
mirror at the top of `tower:.../sources/LR2/` are superseded by the files above and are **not** inputs of any stage.
Card A-29 takes `LR2` from `/tmp/if2453/r19/LR2/` and `lr2_compendia_r19.fna.gz`. Whether to delete the old file on
`/stor` (6.33 GB) is the PI's.

**Open, not changed here (C-66's rules, kept as the card asks):** every Bickhart bin is labelled `mag_hq`, as C-66
labelled them, though the archive holds all 695 DAS Tool bins and no quality table was read; and a single-contig bin
is circular only if its header says so, so the 218 binned contigs that `assembly_info.txt` marks circular (112 bins
hold one contig) are staged linear in their bins.

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

Detached `systemd-run --user` units at `nice -n 19`, `MemorySwapMax=0`, `IOSchedulingClass=idle`, each started after
a look at the active units (no other card's staging unit was active at any time; no build unit of `MV1` or `V3` ran:
the wave was in `PL3`'s Pfam scan) and with no `STOP` file in the root. At most four processes of this card at a time
(one CPU each for the scans).

| Unit(s) | What | Cap | Wall | CPU |
|---|---|---|---|---|
| `r19-scan-1` to `-4` | the wave's 16 stages, 44 files | `CPUQuota=100%`, 8G each | 01:34:51 to 01:52:45 | 36.0 CPU-minutes (*measured*: sum of `cpu_s`, 2,157.7 s) |
| `r19-scan-5` to `-8` | the 23 `r26.11.5` lists, 32 files | same | 01:46:41 to 02:07 | 47.1 CPU-minutes (*measured*, 2,826.4 s) |
| `r19-rescan` | `NX1` and `GV1` with the byte-value count | same | 02:00:27 to 02:03:30 | 3.0 CPU-minutes (*measured*, 182.2 s) |
| `r19-inspect`, `r19-inspect2` | the four `LR2` archives (§4.2); the first stopped on a `tarfile` stream call | `CPUQuota=100%` | about 2 min, 01:45:54 to 01:47:22 | at most 3.5 CPU-minutes (*derived*: wall × 1) |
| `r19-binslist`, `r19-archsha` | `tar tzvf` of the bins; `sha256sum` of the archives | `CPUQuota=100%` | 10.7 s; 02:00:28 to 02:01:00 | 10.6 CPU-seconds (*measured*); at most 0.5 CPU-minutes (*derived*) |
| `r19-lr2-restage` | `restage-lr2` (§4.3) | `CPUQuota=200%`, 16G | 02:04 to 03:28:40 (84.7 min) | at most 169 CPU-minutes (*derived*: 84.7 × 2); its largest process, `bgzip` on one thread, at most 85 (*derived*) |
| `r19-lr2-r232` | `c76 ... run LR2` (§4.4) | `CPUQuota=200%`, 16G | 03:30:19 to 03:31:12 | at most 2 CPU-minutes (*derived*) |
| `r19-scan-newlr2` | the new FASTA | `CPUQuota=100%` | 03:30:19 to 03:40:18 | 10.0 CPU-minutes (*measured*, 598.9 s) |
| `r19-stage-est` | `stage_est.py conf` | `CPUQuota=100%` | 0.5 s | 0.48 CPU-seconds (*measured*) |

No process came near 300 CPU-minutes. In all at most about 4.5 CPU-hours (*derived*: 36.0 + 47.1 + 3.0 + 3.5 + 0.5 + 169 + 2 + 10.0 = 271
CPU-minutes, the upper bounds summed); the re-stage is most of it.

## 7. Tests

- `tests/r1/test_r19_fasta_scan.py` (17 cases, each run with the normal `PATH` and with `PATH=/usr/bin:/bin`):
  plain and gzip (by magic, also under a plain name), CRLF and blank lines, a header not UTF-8, a UTF-8 header not
  ASCII, tab allowed, bytes outside IUPAC with value, offset and the count by value, empty records, a line before the
  first header, compressed bytes written into a FASTA (the `LR2` case: the first offending record and its header), a
  truncated gzip, trailing bytes that are not gzip, concatenated members, a 16 MiB unwrapped record scanned under
  8 MiB of traced memory, the sha256, `objects.json` paths, the command line (atomic JSON, exit by status) and
  `scan-stage` (one process per file, resume, `table`). Truth by hand. Seen failing first: every case errored on the
  missing script (`FileNotFoundError: .../r19_fasta_scan.py`); the byte-value count was added later, its assertions
  seen failing against the earlier script (`KeyError: 'non_iupac_byte_values'`, 2 failed).
- `tests/r1/test_c76_lr2_sag1_nx1_r232.py`, five new `test_r19_*` cases over fixture archives in the layout of the four
  `LR2` inputs (a gzip member named `.fa.gz` among the bins, a plain one, a gzip Singleton member, a circular binned
  contig whose repeat is trimmed only in the unbinned pass): the magic rule and the strict reader's refusals; the bins
  decompressed and the new FASTA clean by `r19_fasta_scan`; genomes, unbinned counts, kept C-66 codes, objects; byte
  identity when run again; refusal of a byte that is not sequence and of an existing target; a C-66 row that comes out
  different is reported and the command exits 1. Seen failing first against C-76's script: 10 of 10 failed
  (`AttributeError: module 'c76_lr2_sag1_nx1_r232' has no attribute 'restage_lr2'` / `'strict_records'`).
- Migrated (the card's template change): `test_templates_pass_stage_est_and_carry_the_cards_values` now expects `LR2`'s
  FASTA `lr2_compendia_r19.fna.gz`, its `SUM`s in this note and this note among its `NOTE`s; seen failing before the
  template was rewritten (`AssertionError: ... '2026-11-27_staged_fasta_scan.md' in ('...2026-11-26_lr2_sag1_nx1_r232.md',)`).
- `env BIOBOX_TREE=/tmp/biobox-wt/R-19 sh scripts/pod_pytest.sh R-19 tests/r1/test_r19_fasta_scan.py
  tests/r1/test_c76_lr2_sag1_nx1_r232.py tests/r1/test_r232_index.py tests/r0`: **308 passed in 18.34 s** (*measured*, `wilkcomp02`, 2026-11-27, one process).

## 8. `NX1` re-staged without its three 16S alignment files, and what ingest does with a gap byte (card R-20, 2026-11-27)

Added 2026-11-27 by card R-20, from §3's report. Host `wilkcomp02`; labels as at the top of this note. **Read only:**
C-66's `/tmp/if2453/c66/NX1/` and `/stor/scratch/Wilke/if2453/bio-box-staging/NX1/nx1_compendia.fna.gz`, C-76's
`/tmp/if2453/c76/NX1/`, the root, `/tmp/if2453/s2/inputs/` (the chain, and `GV1`'s inputs for §8.5). **Written:**
`/stor/scratch/Wilke/if2453/bio-box-staging/NX1/nx1_compendia_r20.fna.gz` (new name, beside the old file), the
metadata in `/tmp/if2453/r20/NX1/`, the r232 status in `/tmp/if2453/r20/c76/NX1/`, the scan in `/tmp/if2453/r20/scan/`,
`/tmp/if2453/r20/gv1_gaps.json`, logs and three small check scripts in `/tmp/if2453/r20/logs/`, and the mirror on
`tower` in the new directory `sources/NX1/r20_restage/`. Nothing under `/stor` or on `tower` was overwritten or
deleted.

### 8.1 The re-stage

`c76_lr2_sag1_nx1_r232.py restage-nx1` reads C-66's staged FASTA and tables and writes them again under new names
without the genomes it is given (by default the three Lost City alignment files of §3). A record is left out only when
no kept genome holds its contig; every kept record (header and sequence lines, line ends as they are) and every kept
table row is copied byte for byte; `genomes.json` is checked to follow `records.tsv` genome for genome and to be in
C-66's compact form, and is written in the same form; a genome C-66's `skani` step dropped for a removed
representative would refuse the run (none was). Every kept record is refused unless its header is UTF-8 and its bytes
IUPAC letters, and the new file is read back: the sha256 of its decompressed bytes must equal that of the kept records
as read from the old file, and each header must be the `contig_uid` (`biobox.ids`) of its record. Run (unit
`r20-restage`, in `/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif`, `CPUQuota=200%`, `bgzip -@ 1`):

```
c76_lr2_sag1_nx1_r232.py restage-nx1 /tmp/if2453/c66/NX1 \
    /stor/scratch/Wilke/if2453/bio-box-staging/NX1/nx1_compendia.fna.gz /tmp/if2453/r20/NX1 \
    /stor/scratch/Wilke/if2453/bio-box-staging/NX1/nx1_compendia_r20.fna.gz /tmp/if2453/s2/inputs/NX1 \
    --compressor bgzip --bgzip-threads 1
```

04:51:13 to 05:01:05 CST (*measured*, unit start and file times). A first start at 04:49 was refused before it wrote
anything: C-66's `genomes.json` ends in a newline the script did not yet expect (the fixture now has it).

**Removed** (*measured*, `/tmp/if2453/r20/NX1/removed_genomes.tsv`; reason "16S rRNA alignment file of the source
archive, not genome sequence"):

| Genome | Rows (contigs) | Bytes | Gap bytes `-` |
|---|---:|---:|---:|
| `lost_city:Bipolaricaulota-16S-phylogeny` | 18 | 28,656 | 9,109 |
| `lost_city:Methanosarcinales-16S-phylogeny` | 34 | 22,984 | 1,074 |
| `lost_city:Thermodesulfovibrionales-16S-phylogeny` | 17 | 28,305 | 7,210 |
| total | 69 | 79,945 | 17,393 |

The 69 rows are 69 distinct contigs, none of them held by a kept genome or by another stage; 17,393 is §3's count of
gap bytes in the whole file, so no gap byte is left (*measured*, `stage_report.json`). The three were code 11, all
`not_visible` in C-76's table.

**The new `NX1` against the old** (*measured*, `/tmp/if2453/r20/NX1/stage_report.json`, C-66's
`/tmp/if2453/c66/NX1/stage_report.json`):

| | old (`nx1_compendia.fna.gz`, C-66) | new (`nx1_compendia_r20.fna.gz`, R-20) |
|---|---:|---:|
| Genomes | 896 | **893** (Lost City 32 → 29; eHOMD 20, ENA lichen and Samail 678, polar mats 20, rhizosphere 146 unchanged) |
| Distinct contigs (FASTA records) | 391,103 | **391,034** |
| `genome_input` rows (`records.tsv` rows) | 391,117 | **391,048** |
| Bases | 3,398,390,142 | **3,398,310,197** |
| Genetic code 11 / 25 / 4 | 889 / 6 / 1 | 886 / 6 / 1 |
| Held by the root and the `r26.11.4` stages | 5 contigs, 84,773 bp (PL1 2, PL3 2, MV1 1) | the same 5 (`shared_records.tsv` byte-identical, sha256 `14f03620…`) |
| New | 391,098 contigs, 3,398,305,369 bp | **391,029 contigs, 3,398,225,424 bp** |
| FASTA bytes | 1,041,939,827 | 984,388,438 |
| FASTA sha256 | `9da4516e264df7125a845848658915706e7079f04bf01b2fdd7ec1bbb17e3863` | `395a840112ce6e1a113487b0b34a57935a5a6ab033d40b8020e86f4cb4639176` |

The differences are the 69 rows (*derived*: 391,103 − 69 = 391,034; 391,117 − 69 = 391,048; 3,398,390,142 − 79,945 =
3,398,310,197; 3,398,305,369 − 79,945 = 3,398,225,424). The old file is 57,551,389 bytes larger than the new
(*derived*), more than the 79,945 bytes removed: the old file is one plain gzip stream (its header carries a file name and
no `BC` field; *measured*, `xxd -l 18`) and the new one is BGZF from `bgzip`, as R-19's `LR2`; so only decompressed
bytes are compared.

**Byte identity of the kept records:** the decompressed stream of the old file less the 69 records and the
decompressed new file have the same sha256, `02b2e17c1ad4c98afc81a401d0f5db939937a11cb94e48c861978783fc88bbf4`
(`kept_records_byte_identical: true`), and all 391,034 headers are the `contig_uid` of their bytes
(`uid_mismatches: 0`) (*measured*, `stage_report.json`). Every table is C-66's less the removed genomes' rows, its bytes
and CRLF line ends kept: `records.tsv` 391,048 kept / 69 out, `genetic_codes.tsv` 893 / 3, `skani_within_dropped.tsv`
26 / 0, `shared_records.tsv` 5 / 0 (sha256 unchanged), `trims.tsv` and `repeated_records.tsv` header only (sha256
unchanged).

Files (*measured*, sha256): in `/tmp/if2453/r20/NX1/`: `genomes.json`
`e31b9b1970d71939595fe82c0b007a8af30d79dee54789716fa97fb45c8e642b`, `objects.json`
`bfb7726a13d85e1f4d03ff16f4031502fe3187e6582d93e449dd5d48f5d215d7` (path
`/tmp/if2453/s2/inputs/NX1/nx1_compendia_r20.fna.gz`, `source_release` `2026-11-27` as R-19 set it for `LR2`; licence and
`redistributable: true` as C-66 recorded them), `objects_staging.json`
`61af1b92a96d12242069db997a25a66d91691d771e01e2d22fe249c2bcc8059a`, `records.tsv`
`d9dcb5f68dcaa2904c469bbe126cbe5761810a802e8c1081208b8f0d719f1955`, `genetic_codes.tsv`
`3c74c3971117acd58f49dab1cab8a05598aa4fe167f1c252388e646c5e61e9b6`, `skani_within_dropped.tsv`
`68bbce115a8b979b242e7191ada63d6f3dc2775b3685a8c9cfa83a56792fdc32`, `removed_genomes.tsv`
`91b6ef91f7bbb3f006093fe7064407cb2a0232b49a60c176c4f3abb262e0ec4e`, `stage_report.json`
`739a0d686195f5e65d2d04e44a71529edb76938beefa9a7aee94e341ef33a3d6`, and `shared_records.tsv`, `trims.tsv`,
`repeated_records.tsv` equal to C-66's.

### 8.2 Checks, the r232 status, the template and the mirror

- **The new FASTA scans clean:** 391,034 records, 3,398,310,197 bases (both equal to the stage's), 0 headers not UTF-8
  or not ASCII, 0 bytes outside IUPAC, 0 empty records, no read error, sha256 `395a8401…` as written (*measured*, unit
  `r20-scan`, 05:01:39 to 05:03:17 CST, `r19_fasta_scan.py scan`, the pod's `/usr/bin/python3`,
  `/tmp/if2453/r20/scan/NX1__nx1_compendia_r20.fna.gz.json`).
- **The r232 status of C-76, kept for the 893:** `c76_lr2_sag1_nx1_r232.py run NX1 /tmp/if2453/r20/NX1
  .../nx1_compendia_r20.fna.gz /tmp/if2453/c74/index /tmp/if2453/r20/c76 --root /tmp/if2453/r10/root --chain-home
  /tmp/if2453/s2/inputs` (unit `r20-r232`, in the image, 05:03:46 to 05:04:22) wrote
  `/tmp/if2453/r20/c76/NX1/r232_status.tsv` (sha256 `d9c191b3d087b13a7c7c47f696d8206952e47c028f88c1431c35d5e244388cb5`)
  and `r232_report.json` (`608b19d9a548c93f22c048107b1a65536d579f27cf97f2af9f76079350ef9cb6`). The status table is
  **C-76's `/tmp/if2453/c76/NX1/r232_status.tsv` less its three rows of the removed genomes, byte for byte** (896 rows
  − 3 = 893; the sha256 of C-76's table less those rows equals the new one; *measured*, `/tmp/if2453/r20/logs/
  r20_status_cmp.py`). Counts: 473 in GTDB R232 and GlobDB r232, 64 in GTDB R232 only, 161 in neither, **195** not
  visible (was 198); 537 in GTDB R232; 698 with an accession; removed 0 for r232. Ownership over the root's 13 stored
  stages and the 19 `r26.11.4` stages: 5 contigs held (PL1 2, 17,964 bp; PL3 2, 13,775 bp; MV1 1, 53,034 bp), 391,029
  contigs and 3,398,225,424 bp new, equal to `stage_report.json` (the template's `EXPECT` lines).
- **`s2_check_inputs.py NX1 genomes.json objects_staging.json`** (run by the same command): exit 0, `problems: []`,
  893 genomes, 391,048 `genome_input` rows.
- **`stage_est.py conf .../stages/templates/NX1.conf --stage NX1 --tree <tree>`** (`python3 -I`, the pod's host Python):
  exit 0, `DEPTH=called`, default gene-calling mode, `CRISPR=yes`, `AFTER=SAG1`, both notes. The template's `SUM` lines
  are the sha256 of `genomes.json`, `objects.json` and `nx1_compendia_r20.fna.gz` above; the old file's `SUM` is gone.
  Cost by plan §3's rate R1 (24.8 CPU-minutes per new Gbp, as C-76 used): about 84 CPU-minutes (*derived*: 24.8 ×
  3.398), as before.
- **Mirror to `tower`** with `scripts/orchestrator/relay_to_tower.sh`, into the new directory
  `tower:/mnt/user/research-bulk/bio-box/sources/NX1/r20_restage/` (`sources/NX1/` holds C-76's first pass and
  `stage_inputs_full/` holds C-66's files; nothing there is overwritten): the eleven metadata files (70,552,344 bytes,
  `DONE: 11 verified, 0 in place and skipped, 0 failed`), the two r232 reports (`DONE: 2 verified`) and the FASTA
  (984,388,438 bytes in 32 s, 30.76 MB/s, sha256 equal to the record, `DONE: 1 verified`), 05:06:24 to 05:07:18 CST
  (*measured*, the script's log); `SHA256SUMS` there has 14 lines, every one verified in these runs.

**Superseded, not deleted:** `/stor/scratch/Wilke/if2453/bio-box-staging/NX1/nx1_compendia.fna.gz` (sha256 `9da4516e…`),
C-66's metadata in `/tmp/if2453/c66/NX1/`, C-76's status of it in `/tmp/if2453/c76/NX1/`, and their mirror in
`tower:.../sources/NX1/stage_inputs_full/` are superseded by the files above and are **not** inputs of any stage (C-76's
first-pass `nx1_compendia_r232.fna.gz` was superseded already, its note §6). Card A-29 takes `NX1` from
`/tmp/if2453/r20/NX1/` and `nx1_compendia_r20.fna.gz`. Whether to delete the old file on `/stor` (1.04 GB) is the PI's.

### 8.3 What ingest does with a gap byte `-`: it stores it

Read from `src/` at commit `77881e8` (no change to `src/`):

- `src/biobox/fasta.py:47-72`, `iter_records`, the only reader of a source object at ingest: a line that starts with
  `>` opens a record and its header is decoded as UTF-8 (line 65; that decode is what refuses a header that is not
  UTF-8); a line before the first header raises (lines 68-69); **every other non-blank line is appended to the sequence
  as bytes (line 70) with no alphabet check**. Nothing in the module looks at the sequence's letters.
- `src/biobox/l0seq/store.py:554-566` (`ingest`, pass 1) and `:1066-1077` (`_pass1_run`, the batched pass 1 of
  `build-stage`): the length is `len(rec.seq)`, the stored-bytes key `ids.raw_digest(rec.seq)` and the contig key
  `ids.contig_uid(rec.seq)`; the only reason to leave a record out is `min_contig_len` (lines 558-563 and 1069-1074).
- `src/biobox/ids.py:10-12`: `contig_uid(seq) = blake3(seq.upper())`, so the `-` is part of the contig's identity.
- `src/biobox/l0seq/store.py:773-775` (`_write_shard`, pass 2): the record's bytes are stored after a length and
  digest check against pass 1.
- `src/biobox/l0seq/store.py:92-96`, `revcomp`: the complement table maps the IUPAC letters only, so a `-` stays a `-`
  on the minus strand. `src/biobox/stage.py:765-772` (`fasta_counts`, the estimate) counts it as a base.

So **a `-` is neither refused nor dropped: it is stored as a base of the contig**, counted in its length and hashed into
its `contig_uid`; slices and the reverse strand return it. Test (behaviour as it is):
`tests/r1/test_c76_lr2_sag1_nx1_r232.py::test_r20_ingest_stores_a_gap_byte_as_sequence` ingests a fixture with a
record ending in `--` (as `GV1`'s do) and one with an inner `-`: ingest succeeds, `SeqStore.get` returns the bytes with
the gaps, `length_of` counts them, the uid is that of the gapped bytes (the same letters without the gaps are not
stored), the minus strand keeps the `-` in place, and `pyrodigal.GeneFinder(meta=True).find_genes` takes the gapped
sequence without an error (how it scores the `-` is not tested here).

**`GV1`'s three contigs** (*measured*, read only, unit `r20-gv1gaps`, 05:05:35 to 05:06:01 CST, `biobox.fasta` and
`biobox.ids` in the image over `/tmp/if2453/s2/inputs/GV1/gvmags_v2_all_18kgenomes.fna` and its `records.tsv`;
`/tmp/if2453/r20/gv1_gaps.json`): all six gap bytes are the last two bytes of three contigs of one genome, `Hoko`
(`confirmed/Hoko.fna`): `Hoko_1` (450,697 bytes, record 13,777), `Hoko_3` (407,935, record 13,779) and `Hoko_4`
(120,468, record 13,780), each ending `--`. `GV1`'s `records.tsv` gives each the `contig_uid` of its bytes with the gaps
(`800cb43e…`, `8bed7beb…`, `a88c0ca1…`), as ingest computes it, so the stage is consistent and **the build will not stop
on them**. What follows if `GV1` is built as staged: the three contigs enter the store with two trailing `-` each (6
bytes that are not bases), their lengths are 2 too long, and their `contig_uid`s are those of the gapped bytes; a later
correction cannot rewrite them (AGENTS §3.10) and would give three new contigs (the uids without the gaps are
`346a88ae…`, `256656ff…`, `6267248e…`). **Proposed, not done:** before `GV1` is built, re-stage it the way this card
re-staged `NX1`, under a new FASTA name, with the trailing `-` of these three records removed (terminal padding of the
source's FASTA, not sequence), their uids and lengths recomputed in `genomes.json` and `records.tsv`, and the
ownership checked again; every other record byte for byte. If `GV1` must be built first, the six bytes are harmless to
the build and can be recorded instead. A separate question for a card (not here): whether ingest, or
`s2_check_inputs.py` (card R-21), should refuse a byte outside IUPAC; today nothing in the path does.

### 8.4 Pod compute used

Detached `systemd-run --user` units at `nice -n 19`, `MemoryMax=16G`, `MemorySwapMax=0`, `IOSchedulingClass=idle`,
each started after a look at the active units (the wave was in `PL3`'s Pfam scan throughout, no `MV1` or `V3` unit ran;
card R-21's `r21-oldlr2` was the one other staging unit, from about 05:01) and with no `STOP` file in the root. One
unit of this card at a time.

| Unit | What | Cap | Wall | CPU |
|---|---|---|---|---|
| `r20-restage` (first start) | refused on reading `genomes.json` | `CPUQuota=200%` | seconds, 04:49 | under 0.1 CPU-minutes (*estimate*) |
| `r20-restage` | §8.1 | `CPUQuota=200%` | 04:51:13 to 05:01:05 (9.9 min) | at most 19.7 CPU-minutes (*derived*: 9.9 × 2) |
| `r20-scan` | the new FASTA | `CPUQuota=100%` | 05:01:39 to 05:03:17 | 1.65 CPU-minutes (*measured*, `cpu_s` 99.0) |
| `r20-r232` | `c76 ... run NX1` (§8.2) | `CPUQuota=200%` | 05:03:46 to 05:04:22 | at most 1.2 CPU-minutes (*derived*) |
| `r20-gv1gaps` | §8.3 | `CPUQuota=100%` | 05:05:35 to 05:06:01 | at most 0.5 CPU-minutes (*derived*) |

Run directly, each under a second: `r20_gj_diff.py` (§8.1's refusal), `r20_status_cmp.py` (§8.2), `stage_est.py conf`.
In all at most about 23 CPU-minutes (*derived*: 0.1 + 19.7 + 1.65 + 1.2 + 0.5); no process came near 300 CPU-minutes.
The units bound `/tmp/if2453` read-write (the root is under it; the scripts only read it) and the NX1 staging directory
read-write for the new file only.

### 8.5 Tests

- `tests/r1/test_c76_lr2_sag1_nx1_r232.py`, five new `test_r20_*` cases (each with the normal `PATH` and with
  `PATH=/usr/bin:/bin`) over a fixture stage in C-66's layout (CRLF tables, compact `genomes.json` with a final newline,
  contig-uid headers wrapped at 80; three alignment genomes with gap bytes interleaved with MAGs; one contig an alignment
  genome shares with a MAG): the new FASTA equal to the kept records written by hand and clean by `r19_fasta_scan`,
  every table C-66's less the removed rows, the removed table and counts, objects; byte identity when run again and the
  old stage untouched; refusals (an existing target, a genome not in `records.tsv`, a kept record with a gap byte, a
  `skani` representative that is removed, a `genomes.json` out of step); the command line; and the ingest behaviour of
  §8.3. Seen failing first: all 8 runs of the four restage cases failed (`AttributeError: module 'c76_lr2_sag1_nx1_r232' has no
  attribute 'restage_nx1'`; the command line `invalid choice: 'restage-nx1'`). The ingest case passed at once, as it
  must: it states the behaviour of `src/` as it is, which this card does not change.
- Migrated (the card's template change): `test_templates_pass_stage_est_and_carry_the_cards_values` now expects
  `NX1`'s FASTA `nx1_compendia_r20.fna.gz`, its `SUM`s in this note and this note among its `NOTE`s; seen failing before
  the template was rewritten (`AssertionError: assert 'planning/measurements/2026-11-27_staged_fasta_scan.md' in
  ('planning/measurements/2026-11-26_lr2_sag1_nx1_r232.md',)`).
- `env BIOBOX_TREE=/tmp/biobox-wt/R-20 sh scripts/pod_pytest.sh R-20 tests/r1/test_c76_lr2_sag1_nx1_r232.py
  tests/r1/test_r19_fasta_scan.py tests/r1/test_r232_index.py tests/r0`: **318 passed in 17.50 s** (*measured*,
  `wilkcomp02`, 2026-11-27, one process).
