# GlobDB release 226 species absent from GTDB: metadata, stage proposal, `G1` inputs, and full-set cost (card C-41)

Date: 2026-11-18, `wilkcomp02`, 09:15 to 15:50 CST (15:15 to 21:50 UTC). Every number carries a label: *measured*
(host, date, command or file), *derived* (formula and inputs), *estimate* (with its basis), or *published* (the
provider's figure, with where it is stated). **No build**: nothing was written in the data root `/tmp/if2453/r10/root`
(bound read-only into every command that read it) and nothing under `/stor` (`gtdb_r226_metadata.tar` bound read-only).

Small results are beside this file in [`2026-11-18_globdb_inputs/`](2026-11-18_globdb_inputs/): `fetch_tables_SUMS`
(every fetched metadata/checksum file: bytes, md5, sha256, UTC time), `table_summary.json`, `phyla.tsv`, and
`stage_proposal.tsv` (`c41_tables.py`), and `G1/` (`stream_summary.json` of `c41_stream_stage.py`, `report.json` of
`make_inputs.py prok`, `check.json` of `s2_check_inputs.py G1`, `inputs_SUMS`, `shared_records.tsv`, and `est.json` of
`c41_est.py`). Scripts: `scripts/stress/build_stage_pod/c41_*` and `make_inputs.py prok`.

## Summary

- **Source entry and fetch (`sources.yaml`, `globdb_r226`)**: written from `globdb.org` and the LISC Vienna release
  directory and committed (`82c97d2`) under DEC-50 (1) and DEC-51 **before** any data byte was fetched. `c41_fetch.sh`
  fetched all seven metadata/checksum files (**54,463,226 B**, every size equal to `Content-Length` and every top-level
  md5 equal to `globdb_r226_md5sums.txt`; *measured*, `fetch_tables_SUMS`). In step (4), `c41_stream_stage.py` streamed
  `globdb_r226_genome_fasta.tar.gz` once over a single connection without storing the archive on `/tmp`:
  **251,371,669,660 B**, md5 `a34a10d72cd688ed518eb6c0e2c30d26` (equal to `globdb_r226_md5sums.txt`), sha256
  `c57db3d185ad60c2a7d79cf517bc37450705b499846feda1f6be6b57dc7b247b`, 306,260 member `.fa.gz` files across 31 chunks,
  extracting only `G1`'s 7,139 genomes (**2,903,786,890 B** on disk; `/tmp` free stayed at 641.7–686.2 GiB; *measured*,
  `G1/stream_summary.json`).
- **What the metadata tables say (`c41_tables.py`)**: GlobDB r226 holds **306,260 species representatives**
  (**875,797,776,137 bp**, 72,692,249 contigs; *measured*, matching *published* 306,260): **143,607 from GTDB R226**
  (465.47 Gbp) and **162,653 non-GTDB species** across 13 MAG catalogues (**410,331,702,192 bp = 410.33 Gbp**,
  48,999,006 contigs; led by `mOTU` 63,092 / 154.57 Gbp, `SPIRE` 49,369 / 128.93 Gbp, `RBG` 23,406 / 55.04 Gbp, `TPMC`
  8,644 / 22.05 Gbp, `SMAG` 6,216 / 18.70 Gbp, `GEM` 3,762 / 10.51 Gbp). All **7 GTDB R226 species representatives
  absent from GlobDB r226** (143,614 vs 143,607; one of them, `GB_GCA_022009895.1`, is already built in `P1`–`P6`) have
  GTDB CheckM2 completeness `< 40%` (28.03%–39.05%) and were dropped by GlobDB r226's post-dereplication `< 40%`
  completeness filter (*measured*).
- **Quality, size, and accession overlap of the 162,653 non-GTDB genomes**: **69,905 genomes (42.98%, *derived*) are
  under 2 Mbp** (89.28 Gbp, 21.76% of bases), **0 are under 100 kbp** (minimum 110,704 bp, so 100% run in `single`
  mode), median CheckM2 completeness is **84.82%** (IQR 71.08%–94.03%), median contamination **1.63%** (IQR
  0.52%–3.61%), with **54,145 `mag_hq`** (33.29%), **103,990 `mag_mq`** (63.93%), and **4,518 low-quality** (`contig`,
  2.78%) (*measured*). By accession (`globdb_r226_dictionaries.tar.gz`), **783 non-GTDB genomes** (`GEM` 458, `SPIRE`
  `SPECIV4` 297, `GOMC` 28) carry INSDC `GCA_`/`GCF_` accessions, of which **205** (by 9-digit assembly stem; 171 exact)
  are non-representative member genomes in GTDB R226 (`> 4%` fastANI from GTDB's chosen species representative),
  including **11** (6 exact) already in `P1`–`P6` or `C1`–`C6` (*measured*).
- **Proposed stages (`G1` to `G41`) and stage order**: **41 stages of ~10 Gbp** (`7.98` to `12.08` Gbp; 37 of 41 within
  `9.61`–`10.45` Gbp) in `(taxonomy, ID)` order with every genus (58,412 genera) kept whole (`stage_proposal.tsv`).
  Build **after both `S3` and `S5`** (`S2C` → `S3` → `S5` → `G1`..`G41`), so every GTDB R226 species representative owns
  any shared contig under DEC-26, and **before `S6`** (GTDB non-rep strains).
- **Stage `G1` inputs and overlap (`/tmp/if2453/s2/inputs/G1/`)**: **7,139 genomes** (all 22 non-GTDB archaeal phyla up
  to *Nitrosotalea*), **1,227,222 distinct contigs** (`10,048,202,556` bp; 65 contigs / 3.20 Mbp shared between pairs of
  `G1` MAGs), `problems: []` in `s2_check_inputs.py`, byte-identical across runs (*measured*). **8 contigs (5,699,697
  bp)** are already owned by stored stages (`P2` 1: *Thermococcus aciditolerans* SY113 complete chromosome; `P4` 6:
  *Halarchaeum* sp. CBA1220 chromosome + 5 plasmids; `H1` 1: `huge_phage_2772` binned inside `SMAGOTU_00531`), and **13
  contigs** are shared with `PL2` (*measured*).
- **Build cost — `G1` and the whole 410.33 Gbp set (`G1`..`G41`)**:
  - For **`G1`**: ingest **41.51 CPU-min** (1 process); gene calling **163.79 CPU-min** in 1 process (**no batches**;
    1-in-20 sample measured `0.6916` s/Mbp → **142.38 CPU-min**, `−12.9%` vs code's `0.7941` s/Mbp, within 20%); CRISPR
    **30.46** and **36.25 CPU-min** (sample JVMs measured `143.08` CPU-s/Gbp across 12,006 contigs, confirming 1.23 M
    contigs/stage add no JVM blowup); full Pfam **23,833 CPU-min (397 CPU-h)** in 100 shards; stage-local NR90 (`2c`)
    scans **97.39%** of residues in the 0.49 Gbp sample (`83.1%` at `rN.3` species-rep scale), costing **20,850–24,254
    CPU-min** including `represent`; **~8.1 GB** in root (*measured* and *derived*, `G1/est.json`).
  - For **all 41 stages (`G1`..`G41`, 162,653 genomes, 410.33 Gbp)**:
    - **Depth `annotated` with Pfam by stage-local NR90 representatives (cards B-47a–e)**: **13,791 to 16,703 CPU-h** —
      **11.4 to 13.8 wall-days** on `wilkcomp02` (40/56-thread schedule; **14.4 to 17.4 days** at flat 40 threads) or
      **107.7 to 130.5 node-hours (SUs)** on TACC Lonestar6 (128 cores/node), plus **~354 GB** in the data root
      (*derived* and *estimate*). Because GlobDB is already dereplicated at 96% ANI (one genome per species), stage-local
      NR90 saves **at most ~13% (`~2,100` CPU-h net)** vs full Pfam (**15,904 to 16,417 CPU-h**, **124.2 to 128.3 SUs**),
      and any combined release over `S3 + S5 + G1..G41` (~780–850 M proteins) exceeds `wilkcomp02`'s 200 GB RAM
      clustering cap (`DEC-29`), requiring card C-38c (`linclust`) or a TACC Stampede3 4 TB node.
    - **Depth `called` (DEC-54, depth 1)**: **142.9 to 161.7 CPU-h** (~162–187 CPU-h allowing extra `igr`/`tapes` sort
      time over 49 M contigs) — **2.8 to 4.0 hours** on `wilkcomp02` or **1.12 to 1.26 node-hours (SUs)** on TACC
      Lonestar6 (**~100× cheaper** than depth `annotated`), **~297 GB** in root, with no Pfam scan and no release
      clustering wall (*derived*).

---

## 1. Source entry and fetch (`sources.yaml`)

GlobDB release 226 (Speth et al. 2025, *Bioinformatics Advances*, `doi:10.1093/bioadv/vbaf280`) is served without login
from the University of Vienna LISC fileshare (`https://fileshare.lisc.univie.ac.at/globdb/globdb_r226/`). Its entry
`globdb_r226` was written into [`sources.yaml`](../../sources.yaml) from the public download, methods, and citations
pages and HTTP `HEAD` headers and committed as `82c97d2` under DEC-50 (1) and DEC-51 **before** any file was fetched.

| File (`globdb_r226`) | Bytes | Provider MD5 (`md5sums.txt`) | Our SHA-256 | Mode & wall time |
|---|---:|---|---|---|
| `globdb_r226_md5sums.txt` | 526 | — | `4a5618b165d515727abf1046fbc77120647f41d1c8db6295f397136efa501499` | stored (`dl/`), 1 s |
| `README.txt` | 10,941 | — | `cf57ec97f22bdcf4cbfe0ebf9e9a055a5f358b5ef081a96cedd09a20eaa9843a` | stored (`dl/`), <1 s |
| `globdb_r226_tax_plus_stats.tsv.gz` | 11,811,303 | `3f587611cc014f33ebdc0b3cdeaf5e58` (matched) | `3efae383ec08539b56544e128b8447c40a4db1001e791e76da24f53c4ce6cb2c` | stored (`dl/`), 1 s |
| `globdb_r226_taxonomy.tsv.gz` | 5,832,515 | `232b47dc34f46d1c916888f1c4524514` (matched) | `b6bccdfc2e8784300f9b2462c85b7a5798c0ab88dfa1d20f502617f154a8dae2` | stored (`dl/`), 1 s |
| `globdb_r226_dictionaries.tar.gz` | 3,176,715 | `96d5183fee37de9f122d9468acdeb8bc` (matched) | `9850e0195b66559cdf387d77fe4c7c260c2751a486a28223c1cd5cdde8cb9383` | stored (`dl/`), 2 s |
| `environmental_metadata/high_level_env__pct.tsv` | 20,736,725 | — | `1c9056f7b6cc99d33875f772d1b77fc0987d05fb6df0984cea5b8efabfb80381` | stored (`dl/`), 2 s |
| `environmental_metadata/high_level_env__value.tsv` | 12,894,501 | — | `f3385db4cf5653e7cd2d3ed9c763e286dec7834f3a58f32d749ab28af9c0b569` | stored (`dl/`), 2 s |
| `globdb_r226_genome_fasta.tar.gz` | 251,371,669,660 | `a34a10d72cd688ed518eb6c0e2c30d26` (matched) | `c57db3d185ad60c2a7d79cf517bc37450705b499846feda1f6be6b57dc7b247b` | streamed once (`G1` extracted), 13,569 s |

All *measured* on `wilkcomp02` (`2026-11-18`, [`fetch_tables_SUMS`](2026-11-18_globdb_inputs/fetch_tables_SUMS) and
[`G1/stream_summary.json`](2026-11-18_globdb_inputs/G1/stream_summary.json); `-A curl`, single connection, zero HTTP
refusals). Because GlobDB distributes all 306,260 genome FASTAs as a single 251.4 GB `.tar.gz` (`R226CHUNK000` to
`R226CHUNK030`, each member already a per-genome `<ID>.fa.gz`) while `r26.11.1` was running alongside on `wilkcomp02`,
`c41_stream_stage.py` streamed the response body through an in-process MD5+SHA-256 wrapper into `tarfile.open(...,
mode="r|gz")`, extracting only `G1`'s 7,139 `.fa.gz` files (`2.90 GB` on `/tmp`, 79.12 CPU-min total). Note on a future
multi-stage fetch: once `/stor` has staging space approved or on TACC `/scratch`, a single stream of
`globdb_r226_genome_fasta.tar.gz` (~3.8 h at 17.7 MiB/s) can demultiplex all 41 `G1`..`G41` directories in one pass
rather than streaming per stage.

**Terms (`CC BY-SA 4.0`).** `globdb.org` states: *"the GlobDB propagates the licenses of the underlying data sources, and
is licensed under CC BY-SA 4.0"* (*published*, read 2026-11-18). `redistributable = true`, with share-alike attribution
to GlobDB and the 13 underlying MAG catalogues (`GEM`, Nayfach et al. 2021, derives from published JGI IMG/M metagenomes
under CC BY 4.0; `SPIRE`, Schmidt et al. 2024, states no explicit licence on its download portal).

---

## 2. What the metadata tables say before any sequence

### 2.1 Total species representatives and the 7 GTDB R226 species GlobDB lacks

Both `globdb_r226_tax_plus_stats.tsv.gz` and `globdb_r226_taxonomy.tsv.gz` hold **306,260 rows** (*measured*,
`table_summary.json`, matching *published* 306,260):

- **GTDB R226 representatives in GlobDB**: **143,607 genomes**, **465,466,073,945 bp (465.47 Gbp)**, 23,693,243 contigs
  (6,967 Archaea, 136,640 Bacteria; stored under versionless `GCA_xxxxxxxxx` / `GCF_xxxxxxxxx` identifiers; 143,607 of
  143,607 match a GTDB R226 species representative with 0 base-count mismatches against `gtdb_r226_metadata.tar`;
  *measured*).
- **Non-GTDB MAG species representatives**: **162,653 genomes**, **410,331,702,192 bp (410.33 Gbp)**, 48,999,006 contigs
  (*measured*).

**Why GTDB R226 has 143,614 species representatives (`gtdb_r226_metadata.tar`, *measured*) while GlobDB r226 keeps
143,607.** Cross-referencing all 143,614 GTDB R226 `gtdb_representative == 't'` rows against GlobDB r226 identifies the
exact **7 missing GTDB R226 species representatives** (`table_summary.json` `missing_7_gtdb_reps`, *measured*):

| GTDB R226 accession | GTDB taxonomy (phylum; species) | bp | Contigs | GTDB CheckM2 comp / cont (%) | NCBI level | In `P1`–`P6`? |
|---|---|---:|---:|---|---|---|
| `GB_GCA_022009895.1` | Bacillota; *Qingrenia sp022009895* | 2,648,256 | 1 | 39.05 / 0.47 | Complete Genome | **yes (built)** |
| `GB_GCA_031290575.1` | Elusimicrobiota; *Endomicrobiellum sp031290575* | 409,366 | 7 | 34.48 / 0.02 | Contig | no |
| `GB_GCA_035294065.1` | Thermoproteota; *ECH-B-1 sp035294065* | 325,681 | 5 | 39.02 / 0.00 | Contig | no |
| `GB_GCA_039592855.1` | Bacteroidota; *Chitinophaga sp039592855* | 2,673,383 | 5 | 38.58 / 0.00 | Contig | no |
| `GB_GCA_916438745.1` | Bacillota; *F0422 sp916438745* | 2,155,950 | 1 | 36.38 / 1.91 | Contig | no |
| `GB_GCA_916439735.1` | Bacillota; *AWVT01 sp916439735* | 1,112,174 | 4 | 37.66 / 1.48 | Contig | no |
| `GB_GCA_934201845.1` | Bacillota; *Ligilactobacillus sp934201845* | 548,054 | 3 | 28.03 / 0.00 | Contig | no |

Every one of the 7 has **CheckM2 completeness `< 40%`**, whereas the minimum completeness across all 306,260 genomes in
GlobDB r226 is **40.01%** (`40.57%` among GTDB genomes; *measured*). GlobDB r226's post-dereplication QC filter
(*"genomes with a CheckM2 completeness under 40 % (341)"*, *published*, `website_pages/globdb_r226_methods`) did not
exempt GTDB representatives in r226 (9 GTDB R226 reps have `< 40%` in GTDB's metadata, and 2 of those 9 scored `>=
40.57%` when GlobDB re-ran CheckM2 v1.0.2; an explicit GTDB exemption was added in GlobDB r232). Because our stages `S3`
and `S5` come directly from GTDB R226 (`GB_GCA_022009895.1` is already in `P1`–`P6`), all 7 remain covered.

### 2.2 Non-GTDB species representatives by resource

GlobDB dereplicated 14 resources sequentially at 96% fastANI over ≥50% aligned fraction (*published*, `globdb_r226_methods`).
In `03_SPIREOTU_dictionary.tsv` (107,078 rows), SPIRE v1 comprises both newly assembled MAG clusters (`SPIREOTU_*`,
92,134 in dictionary, **49,072** in GlobDB) and proGenomes3 `specI_v4` clusters (`SPECIV4_*`, 14,944 two-column rows in
dictionary, **297** in GlobDB), summing to **`49,072 + 297 = 49,369`** — equal to the published SPIRE figure in card
C-41 (*measured* and *derived*, `table_summary.json`):

| Resource (dereplication order) | ID prefix(es) | Genomes | Bases (bp) | Gbp | Contigs | Archaea / Bacteria / Uncl. | `< 2 Mbp` genomes (Gbp) | MIMAG HQ / MQ / LQ | Median size (Mbp) | Median comp / cont (%) | Median N50 (kbp) |
|---|---|---:|---:|---:|---:|---|---|---|---:|---|---:|
| **mOTU** (mOTUs v4.0) | `MOTU40_` | 63,092 | 154,570,503,887 | 154.57 | 16,072,757 | 1,752 / 61,226 / 114 | 28,001 (37.98) | 24,499 / 36,927 / 1,666 | 2.14 | 87.21 / 1.22 | 17.45 |
| **RBG** (Rubin et al. 2025) | `BCRBG_` | 23,406 | 55,038,474,780 | 55.04 | 6,807,308 | 2,203 / 21,203 / 0 | 11,363 (12.54) | 9,615 / 13,791 / 0 | 2.06 | 88.16 / 1.33 | 12.93 |
| **SPIRE** (v1 total) | `SPIREOTU_`, `SPECIV4_` | 49,369 | 128,928,894,257 | 128.93 | 13,835,976 | 2,864 / 45,897 / 608 | 20,908 (26.05) | 14,080 / 34,632 / 657 | 2.28 | 82.23 / 2.24 | 12.28 |
| ↳ *SPIRE MAGs* | `SPIREOTU_` | *49,072* | *128,067,504,300* | *128.07* | *13,795,729* | *2,814 / 45,650 / 608* | *20,825 (25.92)* | *13,839 / 34,581 / 652* | *2.28* | *82.04 / 2.25* | *12.18* |
| ↳ *SPIRE proGenomes3* | `SPECIV4_` | *297* | *861,389,957* | *0.86* | *40,247* | *50 / 247 / 0* | *83 (0.13)* | *241 / 51 / 5* | *2.70* | *96.90 / 0.92* | *47.91* |
| **GEM** (Nayfach et al. 2021) | `GEMOTU_` | 3,762 | 10,510,597,009 | 10.51 | 819,910 | 321 / 3,441 / 0 | 1,470 (1.87) | 1,572 / 2,058 / 132 | 2.40 | 86.60 / 0.88 | 22.26 |
| **MGnify** (human oral/gut, marine, etc.) | `MGYG` | 1,732 | 3,360,848,885 | 3.36 | 625,710 | 123 / 1,578 / 31 | 1,040 (1.22) | 277 / 1,426 / 29 | 1.65 | 73.61 / 0.98 | 7.60 |
| **GOMC** (Global Ocean Microbiome) | `GOMCOTU_` | 2,858 | 7,152,877,229 | 7.15 | 1,262,557 | 236 / 2,621 / 1 | 1,232 (1.51) | 705 / 1,975 / 178 | 2.23 | 80.45 / 1.68 | 7.34 |
| **SMAG** (Ma et al. 2023 soil) | `SMAGOTU_` | 6,216 | 18,695,628,770 | 18.70 | 4,242,279 | 277 / 5,939 / 0 | 1,680 (2.53) | 202 / 4,805 / 1,209 | 2.67 | 62.88 / 4.51 | 4.55 |
| **TPMC** (Tibetan Plateau) | `TPMCOTU_` | 8,644 | 22,052,915,258 | 22.05 | 3,922,305 | 82 / 8,562 / 0 | 3,028 (3.99) | 2,054 / 6,044 / 546 | 2.44 | 79.90 / 2.39 | 9.09 |
| **MRGM / HRGM2** (mouse / human gut) | `HRGM2_` (290), `MRGM_` (44) | 334 | 762,520,614 | 0.76 | 44,329 | 0 / 334 / 0 | 67 (0.12) | 265 / 69 / 0 | 2.15 | 97.12 / 0.91 | 42.32 |
| **cFMD** (fermented foods) | `FMDMAG_` | 336 | 983,068,393 | 0.98 | 87,723 | 6 / 330 / 0 | 71 (0.10) | 168 / 155 / 13 | 2.82 | 91.20 / 1.02 | 28.10 |
| **SHGO** (sheep gut) | `SHGOMAG_` | 1,021 | 1,778,071,044 | 1.78 | 260,934 | 3 / 1,018 / 0 | 706 (0.98) | 164 / 789 / 68 | 1.65 | 76.03 / 2.61 | 12.06 |
| **AMXMAG** (anammox) | `AMXMAG_` | 506 | 2,069,169,853 | 2.07 | 220,300 | 2 / 504 / 0 | 43 (0.05) | 270 / 224 / 12 | 3.63 | 93.16 / 1.74 | 19.16 |
| **GFS** (glacier-fed streams) | `GFSMAG_` | 1,377 | 4,428,132,213 | 4.43 | 796,918 | 3 / 1,374 / 0 | 296 (0.32) | 274 / 1,095 / 8 | 3.09 | 81.26 / 2.52 | 7.38 |
| **All 13 non-GTDB resources** | | **162,653** | **410,331,702,192** | **410.33** | **48,999,006** | **7,872 / 154,027 / 754** | **69,905 (89.28)** | **54,145 / 103,990 / 4,518** | **2.21** | **84.82 / 1.63** | **13.37** |
| *GTDB R226 in GlobDB (reference)* | `GCA_`, `GCF_` | *143,607* | *465,466,073,945* | *465.47* | *23,693,243* | *6,967 / 136,640 / 0* | *38,746 (52.44)* | *91,341 / 52,185 / 81* | *2.88* | *94.99 / 0.83* | *46.45* |
| **All GlobDB r226** | | **306,260** | **875,797,776,137** | **875.80** | **72,692,249** | **14,839 / 290,667 / 754** | **108,651 (141.71)** | **145,486 / 156,175 / 4,599** | — | — | — |

All counts, base totals, and quantiles are *measured* (`table_summary.json`); Gbp/Mbp/kbp unit conversions are *derived*.
Three resources (`mOTU`, `SPIRE`, `RBG`) contribute **135,867 (83.53%, *derived*) of the non-GTDB species** and **338.54
Gbp (82.50%, *derived*) of the bases**; adding `TPMC`, `SMAG`, and `GEM` covers **94.99%** of genomes and **94.93%** of
bases (*derived*).

### 2.3 Domains, phyla, and novel phyla

Across all 306,260 GlobDB genomes there are **234 phyla** (`210` Bacteria, `23` Archaea, plus `1` `Unclassified` group;
*measured*, [`phyla.tsv`](2026-11-18_globdb_inputs/phyla.tsv)). The 162,653 non-GTDB species comprise:

- **Bacteria**: **154,027 genomes** (`397.90` Gbp, *derived*: `410.3317 − 11.0168 − 1.4189` Gbp) across 199 bacterial phyla;
- **Archaea**: **7,872 genomes** (**11,016,809,225 bp = 11.02 Gbp**, 1,428,448 contigs, 6,400 `< 2 Mbp`; *measured* and
  *derived* by summing `archaea_all_phyla`) across 22 archaeal phyla — **more archaeal species than GTDB R226 itself**
  (7,872 non-GTDB vs 6,967 GTDB in GlobDB, +113.0%, *derived*);
- **Unclassified domain**: **754 genomes** (**1,418,904,273 bp = 1.42 Gbp**, 290,487 contigs; `SPIREOTU` 608, `mOTU` 114,
  `MGnify` 31, `GOMC` 1; *measured*). Note on parsing: in `globdb_r226_taxonomy.tsv.gz` these 754 rows carry
  `d__Unclassified;p_Unclassified;c__Unclassified;...` (single underscore after `p`), which caused GlobDB's own table
  builder to split `domain` in `globdb_r226_tax_plus_stats.tsv.gz` as `"Unclassified;p_Unclassified"` and shift ranks by
  one column (*measured*).

Top 18 phyla by non-GTDB species count (`table_summary.json`, *measured*; share column *derived*):

| Domain | Phylum | Non-GTDB genomes | Non-GTDB Gbp | Non-GTDB `< 2 Mbp` | GTDB R226 genomes | Non-GTDB share of GlobDB |
|---|---|---:|---:|---:|---:|---:|
| Bacteria | `Pseudomonadota` | 31,557 | 89.99 | 8,590 | 35,839 | 46.8% |
| Bacteria | `Bacillota` | 26,918 | 50.73 | 17,280 | 24,862 | 52.0% |
| Bacteria | `Bacteroidota` | 20,271 | 57.14 | 5,161 | 16,911 | 54.5% |
| Bacteria | `Patescibacteriota` | 15,320 | 11.85 | 15,293 | 6,435 | **70.4%** |
| Bacteria | `Actinomycetota` | 14,829 | 36.63 | 6,256 | 16,922 | 46.7% |
| Bacteria | `Chloroflexota` | 6,026 | 17.92 | 1,621 | 4,182 | 59.0% |
| Bacteria | `Acidobacteriota` | 5,397 | 22.73 | 241 | 3,758 | 59.0% |
| Bacteria | `Verrucomicrobiota` | 5,123 | 17.07 | 956 | 3,351 | 60.5% |
| Bacteria | `Planctomycetota` | 4,205 | 18.75 | 225 | 3,168 | 57.0% |
| Bacteria | `Desulfobacterota` | 2,523 | 7.92 | 487 | 2,192 | 53.5% |
| Bacteria | `Myxococcota` | 2,039 | 12.33 | 30 | 1,283 | 61.4% |
| Archaea | `Thermoproteota` | 1,990 | 3.06 | 1,603 | 1,544 | 56.3% |
| Archaea | `Nanobdellota` | 1,953 | 1.68 | 1,942 | 1,064 | **64.7%** |
| Bacteria | `Cyanobacteriota` | 1,906 | 6.08 | 672 | 2,718 | 41.2% |
| Bacteria | `Gemmatimonadota` | 1,695 | 5.71 | 105 | 1,398 | 54.8% |
| Bacteria | `Campylobacterota` | 1,685 | 2.86 | 1,331 | 1,311 | 56.2% |
| Archaea | `Thermoplasmatota` | 1,384 | 2.23 | 1,179 | 1,268 | 52.2% |
| Bacteria | `Bacteroidota_A` | 1,361 | 4.21 | 135 | 1,314 | 50.9% |

In addition to the 754 `Unclassified` genomes, GlobDB r226 defines **25 provisional phyla with zero GTDB R226
representatives** (23 bacterial, 2 archaeal: `BCRBG_10398` and `SPIREOTU_01873564`), holding 54 genomes and 136.7 Mbp
in total (*measured*, `table_summary.json` `novel_phyla_zero_gtdb`).

### 2.4 Completeness, contamination, fragmentation, and `< 2 Mbp` distributions

Across the 162,653 non-GTDB species representatives (`table_summary.json`, *measured*; percentages *derived*):

| Metric (162,653 non-GTDB genomes) | Min | p05 | p25 | Median (p50) | p75 | p95 | Max | Mean |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Genome size (bp) | 110,704 | 686,985 | 1,432,128 | 2,214,365 | 3,270,431 | 5,390,564 | 18,089,507 | 2,522,743 |
| CheckM2 completeness (%) | 40.01 | 54.03 | 71.08 | 84.82 | 94.03 | 99.91 | 100.00 | 81.61 |
| CheckM2 contamination (%) | 0.00 | 0.04 | 0.52 | 1.63 | 3.61 | 7.94 | 102.55 | 2.50 |
| Quality score (`comp − 5 × cont`) | −413.04 | 29.68 | 55.52 | 71.95 | 86.23 | 97.43 | 100.00 | 69.12 |
| Contig count per genome | 1 | 24 | 93 | 213 | 412 | 870 | 6,769 | 301.2 |
| Contig N50 (bp) | 1,083 | 3,616 | 6,387 | 13,365 | 35,419 | 156,241 | 14,598,830 | 48,519 |
| Largest contig (bp) | 2,497 | 13,696 | 28,226 | 55,709 | 120,642 | 365,796 | 14,598,830 | 113,623 |

- **Under 2 Mbp**: **69,905 genomes (42.98%)**, **89,278,743,644 bp (89.28 Gbp, 21.76% of bases)**: 2,180 in
  `100–<500 kbp`, 19,955 in `500 kbp–<1 Mbp`, and 47,770 in `1–<2 Mbp` (against 43,432 in `2–<3 Mbp`, 38,427 in
  `3–<5 Mbp`, and 10,889 `>= 5 Mbp`; *measured*). Over half of the `< 2 Mbp` genomes sit in three phyla: `Bacillota`
  (17,280), `Patescibacteriota` (15,293 of 15,320 = **99.82%**), and `Pseudomonadota` (8,590), plus 81.3% of all
  non-GTDB Archaea (6,400 of 7,872, including 1,942 of 1,953 `Nanobdellota`).
- **Under 100 kbp (`SINGLE_MODE_MIN_BP`)**: **0 genomes** (minimum `110,704` bp, `BCRBG_` MAG; *measured*). Every one of
  the 162,653 non-GTDB genomes is `>= 100 kbp` and therefore calls genes in **`single` mode** (`genecall.mode_of`).
- **Completeness & contamination bands**:
  - Completeness: `>= 95%`: 35,346 (21.73%); `90–<95%`: 24,849 (15.28%); `70–<90%`: 64,266 (39.51%); `50–<70%`: 35,889
    (22.06%); `40–<50%`: 2,303 (1.42%) (*measured* and *derived*).
  - Contamination: `0%`: 2,919; `(0, 1]%`: 58,086; `(1, 2]%`: 30,223; `(2, 5]%`: 46,279; `(5, 10]%`: 22,884; `> 10%`:
    2,262 (*measured*). Notice that GlobDB r226's post-dereplication QC removed genomes with CheckM2 completeness `<
    40%` or failing GUNC, **not** genomes with CheckM2 contamination `> 10%` or `comp − 5×cont < 50`: **28,323 genomes
    (17.41%)** have `comp − 5×cont < 50` (mostly `SPIREOTU` 12,943, `mOTU` 6,473, `SMAG` 4,432, `TPMC` 2,593), and 3
    `SPECIV4` genomes exceed 100% CheckM2 contamination (`max 102.55%`) (*measured*).
  - **Mapped to `genome_input.assembly_quality`**: **54,145 `mag_hq`** (`comp > 90` and `cont < 5`, MIMAG completeness
    and contamination thresholds without rRNA/tRNA counts, which GlobDB's table does not report; 33.29%), **103,990
    `mag_mq`** (`comp >= 50` and `cont < 10`, 63.93%), and **4,518 `contig`** (`40 <= comp < 50` or `cont >= 10`, 2.78%)
    (*measured*). Zero non-GTDB rows list a `culture_collection` (0 vs 13,239 in GTDB; *measured*).

### 2.5 Genomes already in our stages by accession

1. **GlobDB's 143,607 GTDB representatives**: **8,311** are already built in `P1`–`P6` (8,313 in `P1`–`P6`'s selection)
   and **1,225** are in `C1`–`C6` (**9,536** in S2C overall; 11,972 GTDB species in GlobDB have at least one complete
   genome in `P1`–`P6`'s selection, *measured*). The remaining **134,071** GTDB representatives (*derived*: `143,607 −
   9,536`, plus the 6 unbuilt `< 40%` GTDB reps) belong to `S3` (`< 2 Mbp`) and `S5` (`>= 2 Mbp`).
2. **GlobDB's 162,653 non-GTDB representatives** (via `globdb_r226_dictionaries.tar.gz`, 107,078 `SPIRE` + 45,599 `GEM`
   + 24,195 `GOMC` + 8 other dictionaries, 0 malformed rows; *measured*):
   - **783 non-GTDB genomes** map to INSDC `GCA_`/`GCF_` assembly accessions: **`GEM` 458**, **`SPIRE` (`SPECIV4`) 297**,
     and **`GOMC` 28** (*measured*, `table_summary.json` `non_gtdb_gca_overlap`).
   - **205 of those 783** (by 9-digit numeric assembly stem; **171** by exact versioned accession) are **non-representative
     member genomes in GTDB R226** (*measured*). Why did GlobDB keep them after starting from GTDB R226's species
     representatives? Because GTDB clusters species at ~95% ANI around a type/representative genome, whereas GlobDB ran
     dereplication against GTDB's *species representatives alone* at **96% fastANI**: a GTDB member genome sitting at
     95.0%–95.9% ANI from its GTDB species representative stays unclustered at 96% ANI and becomes a separate GlobDB
     species.
   - **11 of those 783** (by 9-digit stem; **6** exact versioned accession) are **already in our S2C stages**: **6 in
     `P1`–`P6` built** (3 exact, e.g. `GEMOTU_05512` = `GCF_000747315.1` *Corynebacterium*, `SPECIV4_` / `GEMOTU_` stems)
     and **5 in `C1`–`C6`** (3 exact) (*measured*).
   - **0** non-GTDB GlobDB IDs or dictionary original IDs appear in `V1` or `PL2` (`records.tsv`, *measured*). (Exact
     sequence-level `contig_uid` sharing is larger than dictionary accession matching because `mOTU` and `SMAG` do not
     carry `GCA_`/`GCF_` strings in their dictionaries — §4.2 measures 8 stored contigs in `G1` alone.)

---

## 3. Stage proposal (`G1` to `G41`) and placement relative to `S3`, `S5`, and `S6`

### 3.1 Stage cut rule (`G1` to `G41`)

Following card C-25 step 2 (`2026-11-13_s2c_phase2_fetch.md`) and `c36_stages.py`:

- Take all **162,653 non-GTDB genomes** (`410,331,702,192` bp), ordered deterministically by `(taxonomy, ID)` (using the
  full 7-rank string of `globdb_r226_taxonomy.tsv.gz`).
- Target stage count: `n = round(410,331,702,192 / 10,000,000,000) = 41` (*derived*).
- Keep every **genus** (`d__;p__;c__;o__;f__;g__`, **58,412 distinct genera**, largest `g__Prevotella` at 633 genomes /
  1.83 Gbp and `g__Collinsella` at 1,053 genomes / 2.13 Gbp; *measured*) whole in stage
  `k = 1 + floor(41 × bp_before_first_genome_of_genus / total_bp)`, named **`G1` to `G41`**.

Every genus is contiguous along `(taxonomy, ID)` (`contiguous_in_tax_order: true`), giving **41 stages** with min
**7.98 Gbp** (`G14`, right after 1.83 Gbp `g__Prevotella` lands in `G13` at **12.08 Gbp**), mean **10.01 Gbp**, and **37
of 41 stages within `9.61`–`10.45` Gbp** (*measured*, [`stage_proposal.tsv`](2026-11-18_globdb_inputs/stage_proposal.tsv)).
Cutting at **family** coherence instead (`14,320` families, largest `f__Lachnospiraceae` at 8,009 genomes / 19.65 Gbp and
`f__Bacteroidaceae` at 4,380 genomes / 12.66 Gbp) widens stage sizes to **`4.60`–`15.87` Gbp** (`table_summary.json`
`stage_proposal_family_coherent`), so genus coherence is proposed.

Representative stages from the full 41-stage table ([`stage_proposal.tsv`](2026-11-18_globdb_inputs/stage_proposal.tsv),
*measured*):

| Stage | Genomes | Bases (bp) | Gbp | Contigs | Phyla / Families / Genera | `< 2 Mbp` | `mag_hq` | Taxonomic span |
|---|---:|---:|---:|---:|---|---:|---:|---|
| `G1` | 7,139 | 10,051,402,938 | 10.05 | 1,227,287 | 22 / 860 / 3,449 | 5,907 | 1,558 | Archaea: `p__Aenigmatarchaeota` .. `p__Thermoproteota` (`g__Nitrosotalea`) |
| `G2` | 2,954 | 9,981,375,925 | 9.98 | 1,294,457 | 9 / 208 / 989 | 722 | 812 | Archaea (`g__Nitrosotenuis`..end) + Bacteria (`p__Acidobacteriota` start) |
| `G5`–`G8` | 17,540 | 40,009,435,987 | 40.01 | 4,656,510 | Actinomycetota + 15 minor | 8,975 | 5,722 | `p__Actinomycetota` .. `p__Bacillota` (`c__Bacilli`) |
| `G9`–`G13` | 26,166 | 52,132,465,571 | 52.13 | 5,423,970 | Bacillota → Bacteroidota | 15,020 | 7,388 | `p__Bacillota` (`c__Bacilli`, `c__Clostridia`) .. `g__Prevotella` |
| `G14`–`G19` | 20,319 | 57,952,856,647 | 57.95 | 6,322,716 | Bacteroidota + 18 minor | 5,065 | 7,607 | `p__Bacteroidota` (`f__Bacteroidaceae`..end) .. `p__Calditrichota` |
| `G27` | 12,356 | 10,041,914,035 | 10.04 | 1,118,483 | 4 / 1,314 / 6,573 | 12,192 | 5,663 | `p__Omnitrophota` .. `p__Patescibacteriota` (`g__Nanosynbacter`) |
| `G30`–`G39` | 34,853 | 100,080,198,707 | 100.08 | 11,880,256 | Pseudomonadota + 53 minor | 9,340 | 12,225 | `p__Planctomycetota` end, all `p__Pseudomonadota`, .. `p__Verrucomicrobiota` |
| `G41` | 3,290 | 9,996,260,393 | 10.00 | 1,339,740 | 6 / 203 / 823 | 980 | 1,020 | `p__Verrucomicrobiota` end .. `d__Unclassified` (754 genomes) |
| **All `G1`–`G41`** | **162,653** | **410,331,702,192** | **410.33** | **48,999,006** | **200 phyla / 14,320 fam / 58,412 gen** | **69,905** | **54,145** | all non-GTDB GlobDB r226 species |

### 3.2 Where `G1`–`G41` sit relative to `S3`, `S5`, and `S6`

In `STAGES.md`, `S2C` (`P1`–`P6`, `C1`–`C6`) holds complete prokaryotes, `S3` holds GTDB R226 species representatives
under 2 Mbp (~38,750 genomes, ~52.4 Gbp), `S5` holds the remaining GTDB R226 species representatives (`>= 2 Mbp`,
~95,330 unbuilt genomes, ~358 Gbp), and `S6` holds all remaining GTDB R226 non-representative genomes (~563k strain
assemblies).

**Proposed build order: `S2C` → `S3` → `S5` → `G1`..`G41` (and before `S6`).**

1. **After both `S3` and `S5` (DEC-26 ownership)**: ingest stores each distinct sequence once, in the first stage that
   sees it, and that stage owns the contig's gene calls, IGRs, and tape. Building all GTDB R226 species representatives
   (`S3` and `S5`) before `G1`–`G41` guarantees that whenever a GlobDB MAG shares an exact contig (chromosome, chromid,
   plasmid, or binned cobiont) with a GTDB species representative, the GTDB representative — with curated NCBI assembly
   metadata and translation table — owns it. (Even if one only wanted `< 2 Mbp` candidates early alongside `S3`, `G`
   stages mix `< 2 Mbp` and `>= 2 Mbp` species within genera; only `G27`, 98.7% `< 2 Mbp` Patescibacteriota, is almost
   purely small-genome, and `S5` can still share plasmids/integrons with `< 2 Mbp` MAGs.)
2. **Before `S6`**: `G1`–`G41` add **162,653 distinct species** (`< 96%` ANI from every GTDB rep), whereas `S6` adds
   within-species strains (`>= 95%` ANI to `S2C`/`S3`/`S5`), so `G1`–`G41` deliver far more sequence and protein-family
   novelty per Gbp than `S6`. For the **205 non-GTDB GlobDB genomes that are GTDB R226 non-representative INSDC
   assemblies** (§2.5), either exclude those 205 IDs from `G1`–`G41` (leaving them to `S6` under `ncbi_genomes`) or let
   `G1`–`G41` ingest them first as `96%`-ANI divergent species lineages.

---

## 4. Proposed stage `G1`: inputs, shared contigs, and build estimates

### 4.1 Stage `G1` inputs (`/tmp/if2453/s2/inputs/G1/`)

`make_inputs.py prok` was extended (test-first in `tests/r1/test_make_inputs.py`, `@pytest.mark.covers("C-41")`: 6 tests
failed first, then 23 passed on `wilkcomp02`) with `--quality-column`, `--source-name`, `--license`, and `--root`.
Running `c41_run.sh inputs` produced `/tmp/if2453/s2/inputs/G1/`:

| File | Content | SHA-256 |
|---|---|---|
| `/tmp/if2453/s2/inputs/G1/genomes.json` | 7,139 prokaryote genomes, code 11, linear; 1,558 `mag_hq`, 5,414 `mag_mq`, 167 `contig` | `cf0740cbd3bf3896a17b3c38fff8bdc3bc7fa34e15d477318dac884335f25fdd` |
| `/tmp/if2453/s2/inputs/G1/objects.json` | 7,139 `.fa.gz` files (2,903,786,890 B compressed, 10,087,932,724 B FASTA), `source_name = globdb`, `source_release = r226`, `redistributable = true` | `90f2f3d43c9313988abd665b9f471907e6fbf2a1b7cd4c86c3014cadbaee6fe6` |
| `/tmp/if2453/s2/inputs/G1/shared_records.tsv` | the 8 records whose contig an earlier stored stage owns (`P2` 1, `P4` 6, `H1` 1) | `c98328bae923515ca2453157f591d2e7d3816ed4a4aa3b1e5549158a41ab917b` |
| `/tmp/if2453/s3/globdb/genomes/G1/lock.json` | archive stream md5/sha256 + per-genome `.fa.gz` and FASTA sizes, md5, sha256, record and base counts | `84fa96dfa5353f60da413ac26264b9d52af5344e5d7d2a578a6d1a739f579c5c` |
| `/tmp/if2453/s3/globdb/genomes/G1/selection.tsv` | 7,139 stage-ordered rows with CheckM2 stats, mapped `assembly_quality`, resource, and taxonomy | `54c6ab813087881079c9a9e1c200abb62d3114e8d2765d5c738d28271ee343b0` |

All *measured* ([`G1/inputs_SUMS`](2026-11-18_globdb_inputs/G1/inputs_SUMS), [`G1/report.json`](2026-11-18_globdb_inputs/G1/report.json),
[`G1/check.json`](2026-11-18_globdb_inputs/G1/check.json)):

- **Dry check (`s2_check_inputs.py G1`)**: exit 0, `problems: []`, 7,139 genomes = 7,139 distinct `genome_uid`s (0
  repeated genomes, 0 contigs repeated inside one genome), 1,227,287 `genome_input` rows across 1,227,222 distinct
  `contig_uid`s (`10,048,202,556` distinct bp), 0 circular rows, every object file present (*measured*).
- **Determinism**: a second run of `make_inputs.py prok` into scratch gave byte-identical `genomes.json`, `objects.json`,
  `shared_records.tsv`, and JSON report (`cmp -s`, exit 0; *measured*).
- **Ambiguity runs (`N`)**: 3,714 of `G1`'s 7,139 genomes (52.0%, 105,385 records) contain `N` characters (scaffolded
  MAG contigs; *measured*, `G1/stream_summary.json`), which `ids.contig_uid` and `pyrodigal` handle natively.
- **Single `source_name = "globdb"`**: using `"globdb"` on all 7,139 objects writes one `contig_occurrence` partition
  (`source=globdb/part-0.parquet`, ~147 MB) rather than 13 per-catalogue partitions, many under `/stor`'s 64 MB floor,
  while every record's `accession` (`<GlobDB_ID>_000000000001`) preserves the catalogue prefix.

### 4.2 Shared contigs (`contig_uid`) with stored and unbuilt stages

Comparing `G1`'s 1,227,222 distinct contigs against all 9 stages on disk in `/tmp/if2453/r10/root` (`S1`, `P1`–`P6`,
`H1`, `PL1`), against `PL2` and `V1` inputs, and within `G1` itself (`G1/est.json`,
[`G1/shared_records.tsv`](2026-11-18_globdb_inputs/G1/shared_records.tsv), *measured*):

1. **With the 9 stored stages (`S1`..`PL1`)**: **8 distinct contigs (8 records, 5,699,697 bp)** are already stored:
   - **`P2` (1 contig, 2,165,038 bp)**: `MOTU40_033702` (`mOTU`, *Thermococcus*) is a single-contig genome byte-identical
     to `NZ_CP041932.1` (*Thermococcus aciditolerans* strain SY113 chromosome, complete genome).
   - **`P4` (6 contigs, 3,175,875 bp)**: all 6 contigs of `MOTU40_034106` (`mOTU`, *Halarchaeum*) are byte-identical to
     `NZ_CP054707.1`–`NZ_CP054712.1` (*Halarchaeum* sp. CBA1220 chromosome + plasmids `pCBA1220-01` through `-05`).
     (mOTUs v4 incorporates RefSeq/proGenomes isolates under `MOTU40_` IDs; both assemblies were non-representative
     members in GTDB R226 sitting `< 96%` fastANI from GTDB's chosen species rep, so GlobDB kept them. Because `P2` and
     `P4` own every contig of `MOTU40_033702` and `MOTU40_034106`, `genecall` skips calling both genomes:
     `genomes_not_called = 2`.)
   - **`H1` (1 contig, 358,784 bp)**: contig #1466 of `SMAGOTU_00531` (`SMAG`, soil *Methanosarcina* MAG) is
     byte-identical to **`huge_phage_2772`** (Al-Shayeb et al. 2020 huge phage co-binned into the MAG; `H1` keeps its
     `gv-meta` viral calls under DEC-26).
   - `S1`, `P1`, `P3`, `P5`, `P6`, `PL1`: **0** shared contigs.
2. **With unbuilt stage inputs**: **13 contigs** of `G1` appear in `PL2`'s plasmid inputs, and **0** in `V1` (*measured*).
3. **Within `G1` itself**: **65 contigs (3,200,382 bp)** are listed in two MAGs of `G1` (both within one catalogue —
   e.g. `GEMOTU_10750` and `GEMOTU_39339` share 3 contigs totaling 2.65 Mbp; `SPIREOTU_01281252` and `SPIREOTU_01281444`
   share 4 contigs totaling 305.7 kbp — and across catalogues, e.g. `BCRBG_63717` with `SPIREOTU_00094590`, and
   `SMAGOTU_14796` with `SPIREOTU_01823628`; *measured*). Under DEC-27, both MAGs list the contig and the lesser
   `genome_uid` writes its calls.
4. **Net addition of `G1`**: **1,227,214 new contigs** and **10,042,502,859 new bases** (*measured*).

### 4.3 Sample measurements (`c41_est.py`, 1-in-20 sample) and AGENTS §6 rate checks

Calling every 20th callable genome of `G1` (`357` genomes, `58,146` contigs, `490,133,448` bp) with `genecall._call` in
single mode (1 thread), clustering the resulting distinct proteins with `mmseqs easy-linclust` (`--min-seq-id 0.9 -c 0.8
--cov-mode 0 --cluster-mode 2 --kmer-per-seq 80`), and running MinCED (`-minNR 3`, 50 Mbp batches) on a 1-in-100 sample
(`72` genomes, `12,006` contigs, `96,001,784` bp) measured (`G1/est.json`, *measured*):

- **Single-mode gene calling rate**: **`0.6916` CPU-s/Mbp** (`338.992` CPU-s on `490.13` Mbp; `0.6654` s/Mbp on 300
  genomes `< 2 Mbp` and `0.7488` s/Mbp on 57 genomes `>= 2 Mbp`), **12.9% below** `genecall`'s upper-bound constant
  `0.7941` CPU-s/Mbp (`P6`; *derived*: `(0.6916 − 0.7941) / 0.7941 = −12.91%`, within 20%). High contig fragmentation
  (~163 contigs/genome in the sample vs ~2 in `P2`/`P6`) does **not** slow `pyrodigal` single mode.
- **Gene density, partial genes, and protein length**: **`1,115.54` genes/Mbp** (`546,762` genes), **`15.55%` partial
  genes** running off contig edges (`85,033` genes with `dist_to_contig_end_bp == 0`; `16.73%` in `< 2 Mbp` MAGs, up to
  `29.2%`–`40.3%` in fragmented `SMAG`/`MGnify`/`GFS` MAGs), **`546,169` distinct `protein_uid`s (`99.89%` of genes)**,
  and **`267.56` residues per distinct protein** (`146,132,278` residues; *measured*). Notice that while archaeal compact
  genomes and edge-truncated partial ORFs raise `genes/Mbp` (+22.7% vs `P2`'s `908.99`) and lower `aa/protein` (−17.3% vs
  `P2`'s `323.61`), **distinct residues per Mbp of genome are within 3.0% of `P2`** (`298,148` aa/Mbp in `G1`'s sample vs
  `289,364` aa/Mbp in `P2`, *derived*) — and Pfam scan cost scales with residues, not gene count.
- **Proteins already in `S1`..`PL1`**: **366 of 546,169 (`0.07%`)** (`P3` 164, `P4` 128, `P5` 80, `P2` 52, `P6` 46, `P1`
  32; *measured*). Card B-47b exact-protein copying saves `< 0.1%` here (and likely `< 2%` even after `S3`+`S5`, since
  every non-GTDB genome is `< 96%` ANI from every GTDB species rep).
- **Stage-local NR90 clustering (`linclust`)**: **531,355 clusters (`97.29%` of proteins, `97.39%` of residues)**; only
  **4** of 14,814 non-representative members had `|L_m − L_r| > 0.2 L_r` (`--cov-mode 0 -c 0.8` already enforces `min/max
  >= 0.8` on aligned spans; *measured*). Because a 0.49 Gbp sample under-clusters relative to a full 10 Gbp stage, `97.39%`
  is an upper bound; on `P1 + P2` (10.6 Gbp of one-genome-per-species prokaryotes in `rN.3`), NR90 representatives held
  **`83.1%` of residues** (`B-47_inferred_domains.md` §2c, *measured*).
- **MinCED CRISPR rate on fragmented MAGs**: `13.736` child-JVM CPU-s across `96.00` Mbp (`12,006` contigs in two 50 Mbp
  batches; 12 arrays = `125.0` arrays/Gbp) = **`143.08` CPU-s/Gbp** (*measured*). Compared to `stage.CRISPR_CPU_S_PER_GBP =
  182.0` (`P1`'s full `stage crispr` step including Python `SeqStore` reads and DuckDB Parquet sort, vs child JVMs alone
  here: `−21.4%`), batching 50 Mbp per JVM keeps 1.23 M MAG contigs per 10 Gbp comfortably *below* the code's estimate.

### 4.4 Step and storage estimates for stage `G1` (`10.04 Gbp` new)

| Step / metric (`G1`) | Estimate (`G1/est.json`) | Shards / batches | Basis (`P2` 9.60 Gbp measured in `root` for comparison) |
|---|---:|---|---|
| `ingest` (card B-52) | **41.51 CPU-min** | 1 process (`< 240`) | `0.07 + (10.088 GB − 0.0057 GB) × 4.11` (*derived*; `P2` measured `38.34`) |
| `genecall` (`single`, code rate) | **163.79 CPU-min** | 1 process (**no batches**, `< 240`) | `genecall.cut_plan` at `0.7941` s/Mbp × `1.23` (*derived*; `P2` measured `143.29`) |
| `genecall` (at `G1` sample rate) | **142.38 CPU-min** | 1 process | `0.6916` s/Mbp × `1.23` on `10,042.5` Mbp (*derived*) |
| `crispr` (default, `-minNR 3`) | **30.46 CPU-min** (23.95 sample) | 1 process (`< 240`) | `stage.crispr_cpu_min` at `182.0` CPU-s/Gbp (*derived*; `P2` both CRISPR steps `58.28`) |
| `igr` + `tapes` | **~1.7 CPU-min** (~15–30 *est.*) | 1 process each | `P2` per-Gbp rate (`1.63` CPU-min); *estimate* 15–30 CPU-min sorting 1.23 M contigs |
| Distinct proteins / residues | **≤ 11.19 M** / **≤ 2,994 M aa** | — | from sample (`1,115.5` genes/Mbp, `0.9989` distinct/gene, `267.6` aa; `P2`: `8.58 M` / `2,777 M`) |
| `domains-pfam-38.2` (full scan) | **23,833 CPU-min (397.2 CPU-h)** | **100 shards** | `2,994.1 M × 7.96` CPU-min/M (`23,817` after B-47b copy; `P2` measured `20,630` = 343.8 CPU-h) |
| `represent` (B-47c stage NR90) | **1,059 CPU-min (17.6 CPU-h)** | in parts | scaled from `rN.3` NR90 (`992.2` CPU-min at `10.62 M` proteins, exponent `1.244`; *derived*) |
| `domains-pfam-38.2` (B-47 NR90 reps) | **19,792–23,195 CPU-min (330–387 CPU-h)** | **83–97 shards** | `83.1%` (`rN.3` species-rep rate) to `97.39%` (`G1` sample upper bound) of residues (*derived*) |
| `crispr` (`--near-cas`, `-minNR 2`) | **36.25 CPU-min** (28.50 sample) | 1 process | `1.19 ×` default CRISPR (`CRISPR_NEAR_CAS_COST`; *derived*) |
| Root bytes (`annotated` / `called`) | **~8.07 GB** / **~7.27 GB** | — | `P2`'s `0.7763` GB/Gbp (`7.80` GB) + `0.27` GB for 1.23 M contig index/label/occurrence rows (*derived*/*estimate*) |

**Genetic-code note (`Codetta`).** GlobDB r226's table carries no translation-table column, so `selection.tsv` defaults to
code 11 (appropriate for `G1`'s Archaea). Before building bacterial stages — especially `G27`–`G28` (15,320
Patescibacteriota, where Gracilibacteria/Absconditabacteria use code 25) and `G9`–`G18` (mycoplasmas/Clostridia/Bacteroidota
clades with code 4 or 15) — running Codetta (`STAGES.md` S3) or `pyrodigal` code-11 vs code-25/4 coding-density screening
is needed so alternative-code MAGs are not fragmented under code 11.

---

## 5. Full-set cost (`162,653` genomes, `410.33 Gbp`, stages `G1`–`G41`): pod vs TACC

Because `G1` is Archaea-enriched (mean genome `1.41 Mbp`), `c41_est.py` brackets the 41-stage protein corpus between
**`P2`'s bacterial/archaeal composition** (`908.99` genes/Mbp, `0.9837` distinct/gene, `323.61` aa → **366.9 M**
stage-distinct proteins, **118.74 G** residues; *derived*) and **`G1`'s sample composition** (`1,115.54` genes/Mbp,
`0.9989` distinct/gene, `267.56` aa → **457.2 M** stage-distinct proteins, **122.34 G** residues; *derived*). Total
residues differ by only **3.0%**.

Host conversion factors:
- **`wilkcomp02`**: weekday 12 h × 40 + 12 h × 56 = 1,152 thread-h/day and weekend 24 h × 56 = 1,344 thread-h/day →
  **`1,206.9` thread-h/day** weekly average (`profiles/wilkcomp02.yaml`, assuming no other user computes at night), or
  **`960` thread-h/day** at a flat 40-thread cap (*derived*).
- **TACC Lonestar6** (card C-37): 2 × AMD EPYC 7763 Milan = **128 cores/node**, **1 SU per node-hour**, taking 1
  Lonestar6 core-hour ≈ 1 `wilkcomp02` thread-hour (*estimate* until card C-37 step 3 measures the exact ratio on `P1`).

| Option across all 41 stages (`G1`–`G41`, 410.33 Gbp, 49.00 M contigs) | CPU-hours | `wilkcomp02` wall days (40/56 sched · flat 40) | TACC Lonestar6 node-hours (= SUs) | Data-root storage |
|---|---:|---|---:|---:|
| **Depth `called` (DEC-54, depth 1)**: `ingest` (28.3 h) + `genecall` (97.1–111.5 h) + default `crispr` (16.3–20.7 h) + `igr`/`tapes` (1.2 h *derived*, ~20–35 h *est.* for 49 M contigs) | **142.9 – 161.7 h** (~162–187 h w/ sort *est.*) | **0.12 – 0.13 d** (2.8–3.2 h) · **0.15 – 0.17 d** (3.6–4.0 h) | **1.12 – 1.26 SUs** (~1.3–1.5 w/ sort) | **~297 GB** |
| **Depth `annotated`, Pfam by stage-local NR90 reps (cards B-47a–e, 2c+2b)**: depth `called` + 41 `represent` steps (547.9–720.5 h) + Pfam on NR90 reps (`rN.3` 83.1% rate: 13,081–13,478 h; `G1` sample 97.4% bound: 15,330–15,796 h) + near-Cas `crispr` (19.4–24.7 h) | **13,791 – 16,703 h** | **11.4 – 13.8 d** · **14.4 – 17.4 d** | **107.7 – 130.5 SUs** | **~354 GB** (`+24.6 GB` `stage_rep` + `inferred_from`) |
| **Depth `annotated`, full Pfam 38.2 scan (with B-47b copy)**: depth `called` + Pfam on all non-copied stage proteins (15,741–16,230 h) + near-Cas `crispr` (19.4–24.7 h) | **15,904 – 16,417 h** | **13.2 – 13.6 d** · **16.6 – 17.1 d** | **124.2 – 128.3 SUs** | **~330 GB** |

All CPU-hour and storage figures are *derived* from `G1/est.json` (`full_set_410gbp`), with the 49 M-contig `igr`/`tapes`
sort allowance and `stage_rep` row sizes marked *estimate*.

### Three operational conclusions from these numbers

1. **Stage-local NR90 (`2c`) saves little on species-dereplicated MAG stages.** Cards B-47a–e were designed for S2C phase
   2 (`C1`–`C6`), where multiple genomes per species collapse to `h ≈ 45%` of residues at NR90 (`B-47_inferred_domains.md`
   §6). Because GlobDB is already dereplicated to **one genome per species at 96% ANI**, stage-local NR90 keeps
   **83.1%–97.4%** of residues, saving **at most ~2,100 CPU-h (~13%) net** of the 548–721 CPU-h `represent` steps (and can
   even cost slightly *more* than a direct full scan when a stage's NR90 ratio exceeds ~96.5%). Cross-stage search
   against `S3`+`S5` (`2d` in `B-47_inferred_domains.md`) would save modestly more only on conserved core proteins (`>=
   90%` aa identity across `< 96%` nucleotide-ANI species boundaries).
2. **TACC vs `wilkcomp02` for depth `annotated`.** Annotating all 41 stages with Pfam takes **~2–2.5 weeks** of dedicated
   pod time (on top of `S3` and `S5`'s ~410 Gbp, another ~2–2.5 weeks), whereas on **TACC Lonestar6** the entire 410.33
   Gbp Pfam scan is **~108–131 node-hours (SUs)** — under **3 hours of wall time across 50 compute nodes**, trivial inside
   a startup allocation (card C-37).
3. **Depth `called` (`DEC-54`) + the release-clustering memory wall (`DEC-29`, card `C-38c`).** Ingesting `G1`–`G41` at
   depth `called` (`DEC-54`) costs **137–162 CPU-h** (~3–4 hours on `wilkcomp02` or **1.1–1.3 SUs** on Lonestar6) and
   **~297 GB**, making all 162,653 non-GTDB species immediately searchable by sequence/HMM profile (`protein_seq` /
   `dark_orf`), CRISPR spacers, and local gene-neighborhood synteny (`feature`, `igr`, `contig_tapes`) with **zero Pfam
   cost and zero release-clustering cost**. Conversely, putting `G1`–`G41` (~367–457 M proteins) into an `annotated`
   release alongside `S3`+`S5` (~420 M proteins) reaches **~780–850 M proteins**, far past `wilkcomp02`'s 200 GB RAM
   k-mer prefilter ceiling (~74 M proteins at NR90, `DEC-29`), so an `annotated` release requires either **card C-38c
   (`linclust` mode)** or a **TACC Stampede3 4 TB memory node**.
