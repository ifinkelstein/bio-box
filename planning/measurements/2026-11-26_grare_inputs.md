# GlobDB r226's rare phyla, unclassified domain and remaining archaea as `GRARE1`, `GRARE2`, `GRARE3`, re-keyed to r232 (card C-75, DEC-79)

Date: 2026-11-26, `wilkcomp02`, 14:22 to 18:40 CST. Every number carries a label (`AGENTS.md` §6): *measured* (host,
date, command or file), *derived* (formula and inputs), *estimate* or *published*. Host and date of every *measured*
number are `wilkcomp02`, 2026-11-26, and the command is the subcommand of `c75_grare_inputs.py` named with it (its JSON
report under `/tmp/if2453/c75/`), unless stated. Nothing was built. Read-only throughout (bound `:ro` into every unit):
the data root `/tmp/if2453/r10/root`, the staged `r26.11.4` inputs `/tmp/if2453/s2/inputs`, the r226 archive directory
`/stor/scratch/Wilke/if2453/fragdef/db/globdb-r226/`, `/tmp/if2453/s3` (card C-41's tables) and `/tmp/if2453/c74` (card
C-74's index and re-key). Written: four new files under `/stor/scratch/Wilke/if2453/bio-box-staging/GRARE/` and
everything else under `/tmp/if2453/c75/`. Nothing entered `/tmp/if2453/s2/inputs/` (card A-29 moves the stages).

Script `scripts/stress/build_stage_pod/c75_grare_inputs.py` (`select`, `stream`, `screen`, `finish`); tests
`tests/r1/test_c75_grare_inputs.py`; templates `scripts/stress/build_stage_pod/stages/templates/GRARE{1,2,3}.conf`.

## Summary

- **Three stages, 12,223 genomes, 4,254,800 distinct contigs (4,255,022 records), 32.374 Gbp, 183 phyla** (*measured*), equal to the
  upper end of the wishlist's Z3 (12,223 genomes, 32.38 Gbp, 180 to 195 phyla;
  `2026-11-25_r26_11_5_wide_and_audit_wishlist.md`). **31.840 Gbp are new** against the root's stored stages, the 19
  `r26.11.4` stages and the earlier GRARE stages (*measured*, §5).
- **r232 status** (card C-74's `rekey`, beside every genome): kept 5,698, replaced by a GTDB R232 representative 3,
  dropped by r232's dereplication 6,522 (*measured*). None is left out for it (the orchestrating session's choice,
  STATUS W8). All 754 genomes of the unclassified domain are `dropped` in r232.
- **Codes** (card C-63's screen of 11, 25 and 4): 12,212 genomes code 11, 11 code 25 (all SPIRE genomes of the
  unclassified domain), none code 4, no shared-contig code conflict (*measured*).
- **`s2_check_inputs.py`**: exit 0 and `problems: []` on all three (with `objects_staging.json`); `stage_est.py conf`
  accepts the three templates (§5).
- **The archaeal references** for cards C-70 (`ASG1`) and C-78 (`DPN1`): all 14,839 archaea of GlobDB r226 (6,967 GTDB,
  7,872 not), 2,166,084 records, 23.381 Gbp, in one bgzip FASTA with its list, from the same single pass (§7).
- **Pod compute** (*measured*, the units' reports): one stream, 3 h 38 min wall, 213.7 CPU-minutes in 5 processes (the
  largest about 77.5, *derived*); the screen 16.3 min wall, 127.6 CPU-minutes at 8 workers; no process near 300
  CPU-minutes (§3, §4).

## 1. Inputs, checked first

| input | bytes | checksum (*measured*, `sha256sum` or the reports) |
|---|---:|---|
| `/stor/scratch/Wilke/if2453/fragdef/db/globdb-r226/globdb_r226_genome_fasta.tar.gz` | 251,371,669,660 (`ls -l`) | md5 `a34a10d72cd688ed518eb6c0e2c30d26` by `md5.result` there ("OK"; fetched 2026-10-23 by aria2c, `README.txt` there), equal to `sources.yaml`; sha256 `c57db3d1…c7b247b` (card C-41) |
| `/tmp/if2453/s3/globdb/tables/non_gtdb_stages.tsv` (162,653 rows, card C-41) | 28,963,968 | `1b70c1f97c3344cfa18cde7f5c36ab24978cbce97f41206e8068d755633bd98f` |
| `/tmp/if2453/s3/globdb/dl/globdb_r226_tax_plus_stats.tsv.gz` (306,260 rows) | 11,811,303 | `3efae383ec08539b56544e128b8447c40a4db1001e791e76da24f53c4ce6cb2c` (equal to `sources.yaml`) |
| `/tmp/if2453/c74/rekey_unbuilt.tsv` (130,960 rows, card C-74 `rekey`) | 7,262,952 | `780225326633577e7e87a24e9f11c029184a8b2a1d84350f1a11a72968f8988f` (equal to the C-74 note) |

The `sources.yaml` entry `globdb_r226` said `stored: local only (/tmp on wilkcomp02); nothing under /stor`; the archive
is on `/stor` as above, and the entry is corrected (the card's Why).

## 2. The selection and the cut (`select`)

**Rule** (parameters in `selection.json`): of the GlobDB r226 non-GTDB species representatives that no stage holds
(not `G1`, not a phylum of `GART1`/`GART2`: `r232_index.unbuilt_ids`, 130,960, equal to C-74's count), take those that are
(a) archaeal, (b) of the unclassified domain, (c) in a phylum with no GTDB genome in GlobDB r226 (from the GTDB rows of
`tax_plus_stats`), or (d) bacterial, in a phylum with fewer than **`minor_max_genomes = 1000`** non-GTDB genomes in GlobDB
r226. The wishlist's Z3 does not state its rule; its two ends are reproduced by (d) at 1,000 (12,223 genomes, 32.38 Gbp:
every phylum up to *Bdellovibrionota*, 919) and at 500 (9,472 genomes, 25.04 Gbp, *derived* from card C-41's
`phyla.tsv`: 12,223 less *Bdellovibrionota* 919, *Nitrospirota* 729, *Armatimonadota* 563, *Marinisomatota* 540). The
wishlist's lower end, 9,481 and 25.05, is 9 genomes and 0.007 Gbp more, which is *Undinarchaeota* (9 genomes, 0.0069
Gbp) counted twice (*derived*; a reading of the difference, not of the wishlist's code). The upper end is taken (DEC-59:
the complete catalog is the goal, and three stages of at most 12 Gbp hold it).

**Order and cut**: archaea, then the unclassified domain, then the phyla with no GTDB genome (by name), then the minor
phyla by their non-GTDB count and name; inside a phylum by (taxonomy, ID). A phylum is kept whole: the phylum whose
genomes start at `b` of the `B = 32,374,490,899` selected bases goes to stage `1 + 3 b // B` (card C-41's genus rule at
phylum grain); a stage over `max_stage_bp = 12,000,000,000`, or empty, is refused.

| stage | genomes | Gbp | records | phyla | groups | phyla, first to last (non-GTDB genomes in GlobDB r226) | r232: kept / replaced / dropped |
|---|---:|---:|---:|---:|---|---|---|
| `GRARE1` | 4,779 | 11.583 | 1,659,995 | 162 | archaea 733 (*Thermoproteota* 724, *Undinarchaeota* 9), unclassified 754, no-GTDB phyla 57, minor 3,235 | *Thermoproteota* … *Latescibacterota* (176) | 1,817 / 3 / 2,959 |
| `GRARE2` | 3,806 | 10.456 | 1,253,559 | 15 | minor 3,806 | *Fusobacteriota* (178) … *Chlamydiota* (391) | 1,983 / 0 / 1,823 |
| `GRARE3` | 3,638 | 10.335 | 1,341,468 | 6 | minor 3,638 | *Desulfobacterota_I* (408), *Desulfobacterota_B*, *Marinisomatota*, *Armatimonadota*, *Nitrospirota*, *Bdellovibrionota* (919) | 1,898 / 0 / 1,740 |
| **all** | **12,223** | **32.374** | **4,255,022** | **183** | | | **5,698 / 3 / 6,522** |

All *measured* (`selection.json`, `select_summary.json`). The three `replaced` genomes are `SPECIV4_31793`,
`SPECIV4_07262` and `SPECIV4_26504` (by `GCA_001516765`, `GCA_046533895`, `GCA_001508935`;
`insdc_accession_in_gtdb_r232`). Dropped, by group: archaea 473 of 733, unclassified 754 of 754, no-GTDB phyla 32 of 57,
minor phyla 5,263 of 10,679 (*measured*). That more than half of these rare genomes are not in GlobDB r232 is r232's own
dereplication; the tables do not say what replaced them (C-74 note §4), so a later r232 cut must deduplicate against
these stages by `contig_uid` and ANI.

## 3. The stream (`stream`): one pass, and the 300 CPU-minute bound

**Plan.** C-63's run record (`/tmp/if2453/c63/stage.log`) gives the wall of its pass, 15,402.7 s, but no CPU time; card
C-41's pass of the same archive (`2026-11-18_globdb_inputs/G1/stream_summary.json`) records **79.121 CPU-minutes** for a
whole pass in one process (decompression in Python's zlib, md5 and sha256 of every byte). This card's pass has the same
shape (the image has no `pigz` or `gzip` binary, so `tarfile` reads the archive with zlib in the one process), with
C-63's genetic-code screen moved out to its own step (C-63 ran it inside its pass, in two worker processes), and the four
outputs compressed by four `bgzip -c` child processes. Expected: every process under 300 CPU-minutes, so no split by
member ranges; the stream also stops itself, leaving no file, if its own CPU time passes 280 CPU-minutes (the build
guard's figure; tested).

**Run** (*measured*, `meta/stream/stream.json`, unit `c75-stream`, `CPUQuota=200%`, `MemoryMax=32G`, `nice 19`):
306,260 genome members seen, every selected genome and every archaeon found (`missing: {}`), every member's records and
bases equal to the tables; **wall 13,076.4 s (3 h 38 min)**; the Python process **28.685 CPU-minutes**, its four `bgzip`
children **185.044 CPU-minutes together**. Per child, *derived* from its share of the 55,947,443,208 FASTA bytes written:
archaea 77.5, `GRARE1` 38.5, `GRARE2` 34.7, `GRARE3` 34.3 CPU-minutes. The pass was bound by `bgzip`, compressing about
4 MB/s of FASTA per process and one process at a time (*measured* 30 minutes in from `ps` and `/proc/<pid>/io`: Python
waiting in `pipe_write`), which is also about C-63's pass rate (51.3 Gbp in 15,402.7 s, 3.3 MB/s, *derived*). A later
user of this code (card C-80) can run `bgzip -@ 2` for a shorter wall; BGZF output does not depend on the thread count.

| file (`/stor/scratch/Wilke/if2453/bio-box-staging/GRARE/`) | genomes | records | bp | bytes | sha256 |
|---|---:|---:|---:|---:|---|
| `grare1.fna.gz` | 4,779 | 1,659,995 | 11,583,415,844 | 3,148,662,820 | `558ffa4e9461fc92e7f450d0c2c47133eeea338157209ea453b462d61d32b9c4` |
| `grare2.fna.gz` | 3,806 | 1,253,559 | 10,455,837,641 | 2,825,119,398 | `7d0d7800017c09a8fb3aa40dfd2eb06c84eb37e364fdcc34253c9921e1fc019b` |
| `grare3.fna.gz` | 3,638 | 1,341,468 | 10,335,237,414 | 2,810,598,717 | `6febb9ef5324c7ede8163c54184fd940e3b1d2084c06ae7e98a11e64b592dda2` |
| `globdb_r226_archaea.fna.gz` | 14,839 | 2,166,084 | 23,380,572,012 | 6,273,649,237 | `28e1ce80d2f1271b37d792ca6fd75941f9e2a536c7ec84940464e9b3cb3c8fca` |

All *measured*. Records are in archive order, each with its original header and line width; genomes are whole.

## 4. The genetic-code screen (`screen`)

C-63's `screen_genetic_code`, unchanged: pyrodigal single mode, code 11 on the longest contigs up to 500 kb, the whole
genome with codes 25 and 4 when `cd_11 < 0.82`, the mean ORF is under 600 bp or the taxonomy is a code-25 marker; code 25
or 4 only when coding density gains at least 0.10 and the mean ORF grows at least 1.25 times (C-63's GART1 setting; no
CPR is in these stages). Unit `c75-screen`, 8 worker processes (`CPUQuota=800%`), each replaced after 200 genomes
(*measured*, `meta/screen/screen.json`): **wall 978.1 s, 7.761 CPU-minutes in the main process, 119.819 in the
workers** (about 62 worker processes, about 2 CPU-minutes each, *derived*).

| stage | code 11 | code 25 | code 4 | tested further | conflicts |
|---|---:|---:|---:|---:|---:|
| `GRARE1` | 4,768 | 11 | 0 | 395 | 0 |
| `GRARE2` | 3,806 | 0 | 0 | 19 | 0 |
| `GRARE3` | 3,638 | 0 | 0 | 8 | 0 |

All *measured*. The 11 code-25 genomes are SPIRE genomes of the unclassified domain (`SPIREOTU_00731420`, `_00919215`,
`_01787455`, `_01872133`, `_01960652`, `_02182012`, `_02194008`, `_02213981`, `_02222546`, `_02300713`, `_02756643`).
For each, codes 4 and 25 give the same coding density, so C-63's rule (4 only if denser than 25) gives 25; their code-25
densities are 0.40 to 0.87, low for a prokaryote, and they may not be bacteria at all (an observation, not a test).

## 5. Ownership, the stage files and the checks (`finish`)

Exact prefix `contig_uid` ownership: the root's stored stages `S1, P1..P6, H1, PL1, C5, V1, G1, SPX1` (their
`contig_index`), then the 19 `r26.11.4` stages in `a26_r26_11_4_chain.py`'s order (`R26_11_4_CHAIN`; their
`genomes.json` in `/tmp/if2453/s2/inputs/`), then the earlier GRARE stages; C-63's `_owners_from_records` does the
lookup (C-79's form). The 19 stages' FASTA are not read again (a26's `compute_chain_ownership` does): the lengths are
this card's own records. *Measured*, `meta/finish.json`:

| stage | distinct contigs | bp | held earlier, contigs by stage | held bp | **new contigs** | **new bp** |
|---|---:|---:|---|---:|---:|---:|
| `GRARE1` | 1,659,970 | 11,582,612,217 | 12,169: SPX1 11,830, V3 320, GART1 6, GART2 4, MV1 4, G1 2, P2 1, H1 1, PL3 1 | 213,250,601 | 1,647,801 | **11,369,361,616** |
| `GRARE2` | 1,253,416 | 10,455,174,844 | 8,589: SPX1 8,529, V3 39, PL3 6, GART1 4, MV1 3, P2 2, V2 2, GRARE1 2, G1 1, GART2 1 | 175,004,986 | 1,244,827 | **10,280,169,858** |
| `GRARE3` | 1,341,414 | 10,334,938,936 | 7,825: SPX1 7,768, V3 34, GRARE2 10, PL3 3, MV1 3, V2 3, GART1 2, GART2 2 | 144,765,201 | 1,333,589 | **10,190,173,735** |
| **all** | **4,254,800** | **32,372,725,997** | **28,583** | **533,020,788** | **4,226,217** | **31,839,705,209** |

Totals *derived* (sums). No genome repeats another's contig set (`repeated_records.tsv` holds only its header); 25, 143
and 54 contigs are listed by two genomes of the same stage (DEC-27). Quality (`mag_hq` / `mag_mq` / `contig`):
`GRARE1` 1,361 / 3,261 / 157; `GRARE2` 1,441 / 2,278 / 87; `GRARE3` 1,281 / 2,259 / 98.

**Checks** (*measured*): `s2_check_inputs.py GRARE{1,2,3} genomes.json objects_staging.json` exit 0, `problems: []`
(4,779 genomes and 1,659,995 rows, codes [11, 25]; 3,806 and 1,253,559, [11]; 3,638 and 1,341,468, [11]).
`stage_est.py conf` on the three templates (`--tree`, so every NOTE and EXPECT note must be a file of the tree): exit 0
for each, `DEPTH=called`, `CRISPR=yes`, `AFTER=GART2`, `GRARE1`, `GRARE2` (*measured*, `/tmp/if2453/c75/est.log`).
`stage_est.py before` is card A-29's, once the stages are in `INPUTS`. The templates carry `DEPTH=called`, `CRISPR=yes`, the default
gene-call mode, `AFTER=GART2`, `GRARE1`, `GRARE2`, and `EXPECT=held_by` lines from the table above.

**Build estimate** (*derived*, `biobox.stage` and `biobox.l1.genecall` constants, C-63's formula): `GRARE1` 280.70
CPU-minutes (ingest 50.03, single-mode gene calling 185.08 in 1 batch, CRISPR 45.58), `GRARE2` 251.73, `GRARE3` 250.36;
**782.79 together**, inside the plan's 620 to 800 (plan §4.1); root 22.93 to 27.06 GB.

| file (`/tmp/if2453/c75/meta/<STAGE>/`) | `GRARE1` sha256 | `GRARE2` sha256 | `GRARE3` sha256 |
|---|---|---|---|
| `genomes.json` | `ca6526002510c2e8bd123844d545ba75d4a69d03dce5734052ac14054f96d078` | `6973abf7e0435f416c0ca7c4d20f3b9156b44ad81d192afa9ef74429a2eed485` | `23e1e586691c90c436518171af2fe697e1794ee8bd0067fb3f6cb5b209d6a603` |
| `objects.json` | `773148ba1920cdda3c28df73942718babad9fc74e51e58632c4d1f686fec2715` | `1709e7bd1cbe921924881366bd19c2ac83a14524dc8c2832c6bb3bfce38b6475` | `38b1de83ccfc54819702d89a2dcd9243a78681d0d46a09c555ab108392bf0f23` |
| `objects_staging.json` | `37e3194153238a2daa83ed4a8aae21ad4fcfa37030a6e2a418700733f3159eb5` | `fba980aceca0b92a923e122c4f8905dd6172aa72468be015925ef3f688f1c215` | `39705d72c6be3320e93a0f0e1fc31e0648e966597fd4a7a8523d1c8003ca08d8` |
| `records.tsv` | `c1884e2788f672581696ccbfeb7b01f7fa2c7b7397a50e02981a7e415ff1f68e` | `eb732f71a02b05c8316ceda5afce58297ad417e2e6b94faec4fc8418be878dad` | `b676d2b482717d1901236ec2aa424d197e6e6ef607cd66585eb2f8d4d966c24c` |
| `genetic_codes.tsv` | `43547a645aa8274c7b536d5c650461d82ec6ad77d6c808bb4d0da4f2af77d682` | `e5572d02bd2e4b9c4733ee4eee1ccf775e1b5a6178749cdaf9b71a31e5085566` | `373445b79e96c444f6dfb8afe37506c42de6b8f1ec5ca8ab3e0413ecf89bc103` |
| `shared_records.tsv` | `daa7cbd764d189aab9da110b1146e9a0101eaec02e60c2ab4d017053c1755063` | `736bcd1837179e2961d298ff90089a4479e549e44892c535a775cdf513dfb7f3` | `9bcef30500c3c0e7dc628a89a2d114c85c1235379ac98f611d4de976d632f89d` |
| `repeated_records.tsv` (header only) | `d109495ee76e4643c181eae183c452c9a84560f97a79bbe22d6630780368c565` | the same | the same |

Bytes (*measured*): `records.tsv` 253,587,041 / 190,113,744 / 204,279,525; `genomes.json` 60,271,173 / 45,535,366 /
48,682,114; every other file under 1 MB. `records.tsv` is C-63's form with four columns added: `group`, `r232_status`
(`kept`, `replaced_by <id>`, `dropped`), `r232_id`, `r232_reason`; `genetic_codes.tsv` adds the same four. The metadata
stays on `/tmp` (the plan's §2 rule for metadata, and AGENTS §7.3 keeps files under 64 MB off `/stor`), the FASTA on
`/stor`; `objects.json` names `/tmp/if2453/s2/inputs/<STAGE>/grare<n>.fna.gz`, where card A-29 puts it.

## 6. The remainder (`CATALOG_GOAL.md` §4.1)

`/tmp/if2453/c75/select/remainder.tsv` (118,737 rows, 26,346,957 bytes, sha256
`803af01eb6f73b952452dac20ecc913178ff6d09d94ac0f7968c20bda0c4ed0f`) lists the unbuilt genomes not selected, each with
its r232 key: **118,737 genomes, 316,639,582,514 bp, 37,085,309 contigs, in the 12 phyla with 1,000 or more non-GTDB
genomes** (kept 64,616, replaced 70, dropped 54,051; *measured*, `selection.json`). 12,223 + 118,737 = 130,960
(*derived*).

| phylum | genomes | Gbp | kept / replaced / dropped |
|---|---:|---:|---|
| *Pseudomonadota* | 31,557 | 89.985 | 18,478 / 20 / 13,059 |
| *Bacillota* | 26,918 | 50.734 | 13,916 / 19 / 12,983 |
| *Bacteroidota* | 20,271 | 57.143 | 11,168 / 16 / 9,087 |
| *Actinomycetota* | 14,829 | 36.625 | 8,566 / 7 / 6,256 |
| *Chloroflexota* | 6,026 | 17.922 | 2,829 / 2 / 3,195 |
| *Acidobacteriota* | 5,397 | 22.731 | 2,653 / 1 / 2,743 |
| *Verrucomicrobiota* | 5,123 | 17.071 | 2,607 / 0 / 2,516 |
| *Desulfobacterota* | 2,523 | 7.921 | 1,111 / 2 / 1,410 |
| *Gemmatimonadota* | 1,695 | 5.715 | 927 / 0 / 768 |
| *Campylobacterota* | 1,685 | 2.859 | 1,039 / 3 / 643 |
| *Bacteroidota_A* | 1,361 | 4.207 | 718 / 0 / 643 |
| *Spirochaetota* | 1,352 | 3.728 | 604 / 0 / 748 |

All *measured* (a read of `remainder.tsv`, `/tmp/if2453/c75/c75_summary.py`, output `select_summary.json`).

## 7. The archaeal references (cards C-70 and C-78)

`globdb_r226_archaea.fna.gz` (§3) holds every genome whose `tax_plus_stats` domain is Archaea: 14,839 genomes (6,967
GTDB R226 representatives and 7,872 others, of which 7,139 are `G1`'s and 733 `GRARE1`'s), 2,166,084 records,
23,380,572,012 bp (*measured*). Its list `/tmp/if2453/c75/meta/archaea/archaea_list.tsv` (14,839 rows, 2,481,692 bytes,
sha256 `2d627c856c45523b97e914303a41e12b0359d4343bdc154314c10fe3bb0836ae`) gives per genome `gtdb`, `first_record` (the
0-based position of its first record in the FASTA, archive order), `records`, `bp`, `headers_start_with_genome` (true
for all 14,839, *measured*: each header is `<genome>_<12 digits>`, SPIRE's with a double underscore) and the lineage;
`/tmp/if2453/c75/select/archaea.tsv` (sha256 `d66cd9f3336abed1f5593248fcad01018bee39f455cc31760284db808dab63cb`) is the
same set as selected, sorted by genome.

## 8. The mirror to `tower`

`scripts/orchestrator/relay_to_tower.sh` (card C-49), one stream through the laptop, after a `--dry-run`, 18:38 to
18:47 CST (*measured*, its log): per stage, the seven metadata files from `/tmp/if2453/c75/meta/<STAGE>/` and the FASTA
from `/stor/.../GRARE/` into `tower:/mnt/user/research-bulk/bio-box/sources/<STAGE>/`; and `selection.tsv`,
`remainder.tsv`, `archaea.tsv`, `selection.json` from `/tmp/if2453/c75/select/` into
`tower:/mnt/user/research-bulk/bio-box/deferred_lists/GRARE/`. Every run ended `DONE: n verified, 0 in place and skipped,
0 failed` (7, 1, 7, 1, 7, 1 and 4 files); every sha256 on `tower` equals the record above, and each directory's
`SHA256SUMS` was written. The FASTA went at 18.1 to 23.8 MB/s (*measured*), about twice card C-49's 9.70. The archaeal
references (`globdb_r226_archaea.fna.gz`, 6.27 GB) are not mirrored: they are not a stage, and the card names only the
three.

## 9. Commands (all on `wilkcomp02`, inside `bio-box_faec1f4.sif`, through `/tmp/if2453/c75/c75_run.sh`)

The tree is a copy of this branch's `src` and `scripts` in `/tmp/if2453/c75/tree`; `c75_run.sh` exits if the data root
holds a `STOP` file, and binds the root, `/tmp/if2453/s2/inputs`, `/tmp/if2453/s3`, `/tmp/if2453/c74` and the archive
directory read-only. Each step is one detached unit, `systemd-run --user --unit=<u> --collect -p MemoryMax=<m> -p
MemorySwapMax=0 -p CPUQuota=<q> nice -n 19 sh /tmp/if2453/c75/c75_run.sh <subcommand> ...`:

| unit | quota | wall | subcommand |
|---|---|---|---|
| `c75-select` | 100% | under 1 min | `select /tmp/if2453/s3/globdb/tables/non_gtdb_stages.tsv /tmp/if2453/s3/globdb/dl/globdb_r226_tax_plus_stats.tsv.gz /tmp/if2453/c74/rekey_unbuilt.tsv /tmp/if2453/c75/select` |
| `c75-stream` | 200% | 3 h 38 min (14:24 to 18:02) | `stream /tmp/if2453/c75/select <archive> /tmp/if2453/c75/meta /stor/scratch/Wilke/if2453/bio-box-staging/GRARE --bgzip --cpu-guard-min 280` |
| `c75-screen` | 800% | 16.3 min (18:04 to 18:21) | `screen /tmp/if2453/c75/select /tmp/if2453/c75/meta <GRARE dir> --workers 8 --screen-max-bp 500000 --max-tasks-per-child 200` |
| `c75-finish` | 200% | about 2.5 min (18:25 to 18:27) | `finish /tmp/if2453/c75/select /tmp/if2453/c75/meta <GRARE dir> --root /tmp/if2453/r10/root --chain-home /tmp/if2453/s2/inputs --work /tmp/if2453/c75/work --threads 2 --memory-limit 8GB` |
| `c75-est` | 100% | seconds | `python scripts/stress/build_stage_pod/stage_est.py conf stages/templates/GRARE<n>.conf --stage GRARE<n> --tree /tmp/if2453/c75/tree`, n = 1, 2, 3 |

At most two staging units ran at once (with `c76-nx1` or `c77-inspect`), except the 1.2-second `c75-summary` (a read
of `selection.tsv` for this note) while `c76-nx1` and `c75-stream` ran. Logs: `/tmp/if2453/c75/{select,stream,screen,
finish}.{log,err}`, `meta/stream/stream.json`, `meta/screen/screen.json`, `meta/finish.json`.
