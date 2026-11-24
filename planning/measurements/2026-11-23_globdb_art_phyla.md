# GlobDB r226 ART-RT-enriched bacterial phyla (`GART1` and `GART2`): genetic-code audit, overlap, and staging (card C-63, DEC-71, DEC-72)

Date: 2026-11-23 to 2026-11-24, `wilkcomp02`. Every number carries a label (`AGENTS.md` §6): *measured* (host, date, command/file), *derived* (formula and inputs), or *published*. Host and date of every *measured* number below are `wilkcomp02`, 2026-11-23/24, unless stated otherwise. **No build and no network fetch**: `/stor/scratch/Wilke/if2453/fragdef/db/globdb-r226/`, `/tmp/if2453/r10/root`, `/tmp/if2453/s2/inputs`, and `/tmp/if2453/c62/MV1` were bound read-only (`:ro`) throughout.

Script: `scripts/stress/build_stage_pod/c63_globdb_art_phyla.py`; tests: `tests/r1/test_c63_globdb_art_phyla.py`.

## Summary

- **Two-stage split (`GART1` and `GART2`)**: `stage_driver.sh` enforces `^[A-Z][A-Za-z0-9]{0,15}$` (no underscores), so the five ART-RT-enriched non-GTDB `GlobDB r226` bacterial phyla (`24,554` species representatives, `51.266 Gbp`, `6,431,388` contigs; *measured* from `non_gtdb_stages.tsv` and `/stor/scratch/Wilke/if2453/fragdef/db/globdb-r226/globdb_r226_genome_fasta.tar.gz`) are staged as two stages:
  - **`GART1`** (*Myxococcota*, *Myxococcota_A*, *Cyanobacteriota*, *Omnitrophota*, *Planctomycetota*): **`9,234`** distinct genomes (`0` repeated), **`5,218,158`** distinct contigs (`5,218,293` `genome_input` rows, `135` contigs shared within `GART1`), **`39,412,657,567` bp (`39.413 Gbp`)**, `100%` `genetic_code = 11` (`167` triggered alt-code checks, `0` UGA-recoded), `3,072` `mag_hq` / `5,893` `mag_mq` / `269` `contig`. Net new after `S1..SPX1, V2, PL2, PL3, CG1, MGE1, MGE2, MV1`: **`5,194,823` contigs, `38,962,258,398` bp (`38.962 Gbp`)** (*measured*).
  - **`GART2`** (*Patescibacteriota* / CPR): **`15,320`** distinct genomes (`0` repeated), **`1,212,969`** distinct contigs (`1,213,095` `genome_input` rows, `125` contigs shared within `GART2`), **`11,852,486,932` bp (`11.852 Gbp`)**, **`14,803` `genetic_code = 11` + `517` `genetic_code = 25`** (`1,175` triggered alt-code checks; `488` code-25 genomes in `c__JAEDAM01` *Absconditabacterales*/*BD1-5*, `5` in `c__Gracilibacteria`, `23` in `c__Minisyncoccia`, `1` conflict-aligned in `c__Microgenomatia`), `7,238` `mag_hq` / `7,875` `mag_mq` / `207` `contig`. Net new after `S1..SPX1, V2..MV1, GART1`: **`1,203,001` contigs, `11,665,537,509` bp (`11.666 Gbp`)** (*measured*).
- **Validation & storage placement**: Both stages pass `s2_check_inputs.py` with exit `0` and `"problems": []`. Consolidated deterministic `bgzip` FASTA archives (`>> 64 MB`, obeying `AGENTS.md` §7.3) are staged under `/stor/scratch/Wilke/if2453/bio-box-staging/GART/` (`gart1_globdb.fna.gz` `10.60 GB`, `gart2_globdb_cpr.fna.gz` `3.20 GB`), and compact metadata (`833 MB` total, `< 1.5 GB`) lives under `/tmp/if2453/c63/{GART1,GART2}/` while `r26.11.3` finishes on `/tmp`.
- **Estimated `--depth called --crispr` build cost** (*derived*, `biobox.stage` / `biobox.l1.genecall`): **`957.86` CPU-min** (`28.05–33.12 GB` root) for `GART1` and **`283.87` CPU-min** (`8.40–9.92 GB` root) for `GART2` — **`1,241.73` CPU-min total** (~31 min wall at 40 threads).

---

## 1. Phylum selection and `C-41` (`G2..G41`) coverage

Selected from `/tmp/if2453/s3/globdb/tables/non_gtdb_stages.tsv` (`162,653` non-GTDB species representatives; `0` genomes `< 100 kb`, so 100% run in `single` mode) and streamed in one pass from `/stor/scratch/Wilke/if2453/fragdef/db/globdb-r226/globdb_r226_genome_fasta.tar.gz` (`251,371,669,660` B, verified `md5.result: OK`; `0` missing genomes, `0` contig/base mismatches; *measured*, `gart_report.json`):

| Stage | Phylum | Non-GTDB genomes | Selected bp | Gbp | Contigs |
|---|---|---:|---:|---:|---:|
| `GART1` | `Planctomycetota` | `4,205` | `18,745,296,093` | `18.745` | `2,393,610` |
| `GART1` | `Myxococcota` | `2,039` | `12,332,625,778` | `12.333` | `1,716,101` |
| `GART1` | `Cyanobacteriota` | `1,906` | `6,084,980,130` | `6.085` | `832,656` |
| `GART1` | `Omnitrophota` | `865` | `1,346,636,116` | `1.347` | `158,165` |
| `GART1` | `Myxococcota_A` | `219` | `903,616,316` | `0.904` | `117,761` |
| **`GART1` total** | **5 phyla** | **`9,234`** | **`39,413,154,433`** | **`39.413`** | **`5,218,293`** |
| **`GART2` total** | **`Patescibacteriota` (CPR)** | **`15,320`** | **`11,853,071,408`** | **`11.853`** | **`1,213,095`** |
| **Combined `GART1 + GART2`** | **6 phyla** | **`24,554`** | **`51,266,225,841`** | **`51.266`** | **`6,431,388`** |

**Overlap with `C-41`'s proposed `G2..G41` stages** (*measured*):
- `GART1` pulls forward `9,234` genomes from `G22` (`1,906`), `G25` (`989`), `G26` (`1,657`), `G27` (`477`), `G28` (`1,977`), `G29` (`1,891`), and `G30` (`337`).
- `GART2` pulls forward `15,320` genomes from `G27` (`11,856`) and `G28` (`3,464`).
- **Remaining unbuilt `G2..G41` non-GTDB species representatives** (*derived*): `155,514 - 24,554 = 130,960` genomes, `400,280,299,254 - 51,266,225,841 = 349,014,073,413` bp (`349.014 Gbp`), `47,771,719 - 6,431,388 = 41,340,331` contigs.

---

## 2. Genetic-code audit (`pyrodigal` code 11 vs code 25 / code 4)

Every `GART2` genome (`15,320` CPR genomes, full genome) and every `GART1` genome (`9,234` genomes, longest-contig `>= 500 kb` pre-screen escalating to full genome whenever triggered) was screened with `pyrodigal.GeneFinder(meta=False)` (`translation_table=11`). Any genome with `cd_11 < 0.82`, `mean_orf_11 < 600 bp`, or taxonomy in `c__Gracilibacteria` / `c__Absconditabacteria` / `o__Absconditabacterales` / `o__BD1-5` was tested on the full genome with `translation_table=25` (and `translation_table=4` + internal `TGA/kb` whenever `cd_25 - cd_11 >= 0.10` and `mean_orf_25 / mean_orf_11 >= 1.25`).

- **`GART1` (`9,234` genomes)**: `167` genomes triggered alt-code testing; **`0`** met UGA-recoding criteria -> **`9,234 / 9,234` (`100%`) `genetic_code = 11`**, `0` within-stage code conflicts (*measured*, `/tmp/if2453/c63/GART1/genetic_codes.tsv`).
- **`GART2` (`15,320` genomes)**: `1,175` genomes triggered alt-code testing; **`517` assigned `genetic_code = 25`** and **`14,803` assigned `genetic_code = 11`** (*measured*, `/tmp/if2453/c63/GART2/genetic_codes.tsv`):

| Class (`Patescibacteriota`) | Order | Code 11 | Code 25 | Median `cd_11` | Median `cd_25` | Median ORF 11 (bp) | Median ORF 25 (bp) | Median internal `TGA`/kb |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| `JAEDAM01` | `Absconditabacterales` | `2` | `256` | `0.556` | `0.893` | `280.8` | `899.5` | `6.83` |
| `JAEDAM01` | `BD1-5` | `2` | `232` | `0.608` | `0.914` | `306.8` | `866.0` | `5.93` |
| `Minisyncoccia` | `UBA6257` | `586` | `12` | `0.742` | `0.910` | `427.8` | `902.1` | `2.33` |
| `Minisyncoccia` | `UBA9973` | `2,906` | `11` | `0.738` | `0.874` | `443.0` | `808.4` | `2.07` |
| `Gracilibacteria` | `GCA-2401425` | `1` | `5` | `0.229` | `0.909` | `212.9` | `952.6` | `19.79` |
| `Microgenomatia` | `Levybacterales` | `265` | `1`* | — | — | — | — | — |
| All other CPR orders | (241 orders) | `11,041` | `0` | — | — | — | — | — |
| **Total `GART2`** | | **`14,803`** | **`517`** | | | | | |

\* **Within-stage shared-contig `genetic_code` conflicts (`DEC-27`)**: `3` contigs in `GART2` were shared between a code-11 MAG and a code-25 MAG (`483491b0…` `BCRBG_66247` [11] vs `BCRBG_66789` [25]; `33f97425…` `BCRBG_00616` [25] vs `BCRBG_56476` [11]; `ba301ae3…` `GEMOTU_03525` [25] vs `GEMOTU_09538` [11]). Because `genome_input` requires genomes sharing a contig to agree on `genetic_code`, each conflict was resolved deterministically to the least `genome_uid`'s code (moving `BCRBG_66789` to `11`, and `BCRBG_56476` / `GEMOTU_09538` — the latter in `Levybacterales` — to `25`).

---

## 3. `contig_uid` overlap against `S1..SPX1, V2, PL2, PL3, CG1, MGE1, MGE2, MV1` (`DEC-26`)

Build-order ownership across all 13 root stages (`S1, P1..P6, H1, PL1, C5, V1, G1, SPX1`) + staged `V2, PL2, PL3, CG1, MGE1, MGE2, MV1` (plus `GART1` ahead of `GART2`; all *measured*, `gart_report.json`):

| Earlier owner stage | `GART1` contigs held | `GART1` bp held | `GART2` contigs held | `GART2` bp held |
|---|---:|---:|---:|---:|
| `SPX1` | `23,277` | `437,395,540` | `9,938` | `186,084,887` |
| `PL3` | `42` | `2,253,480` | `1` | `12,013` |
| `MV1` | `5` | `233,797` | `13` | `551,732` |
| `G1` | `0` | `0` | `10` | `25,974` |
| `V2` | `4` | `234,542` | `1` | `40,117` |
| `GART1` | — | — | `4` | `13,486` |
| `P2` | `3` | `9,659,941` | `0` | `0` |
| `PL2` | `3` | `411,536` | `0` | `0` |
| `H1` | `1` | `210,333` | `1` | `221,214` |
| **Total held earlier** | **`23,335`** | **`450,399,169`** | **`9,968`** | **`186,949,423`** |
| **Net NEW in stage** | **`5,194,823`** | **`38,962,258,398`** | **`1,203,001`** | **`11,665,537,509`** |

---

## 4. Staged files, checksums, and `s2_check_inputs.py`

| Stage | File | Location | Bytes | sha256 |
|---|---|---|---:|---|
| `GART1` | `gart1_globdb.fna.gz` | `/stor/scratch/Wilke/if2453/bio-box-staging/GART/` | `10,596,793,098` | `c679c6316be2f02434b1bcb32c6e188e4a3f140e700f7b9ed9701efd92a07c3b` |
| `GART1` | `genomes.json` | `/tmp/if2453/c63/GART1/` | `188,846,586` | `55fe43800b7d65ada245d3b1b0ab7e890563a38707554b8dc8d2d57ee0bcfd11` |
| `GART1` | `objects.json` | `/tmp/if2453/c63/GART1/` | `336` | `924cca2ee65d08eee079d7e3741edeaa2132a8f6762ddf25ec06f0cfb2383868` |
| `GART1` | `records.tsv` | `/tmp/if2453/c63/GART1/` | `512,199,774` | `702c45e4acee51c5953d439ed498f2ea4279d54d5def9319592708ffc4bf9e15` |
| `GART1` | `repeated_records.tsv` | `/tmp/if2453/c63/GART1/` | `39` | `d109495ee76e4643c181eae183c452c9a84560f97a79bbe22d6630780368c565` |
| `GART1` | `shared_records.tsv` | `/tmp/if2453/c63/GART1/` | `1,701,465` | `ec0f20232ff6fe761c984878a353d6dae9f16a239484b2604f465f5e2612a2f8` |
| `GART1` | `genetic_codes.tsv` | `/tmp/if2453/c63/GART1/` | `1,229,428` | `d618ac3bb903b1acc07f876a7d6b8d0d5abaa0d40907208b20b52a727f12b7c1` |
| `GART2` | `gart2_globdb_cpr.fna.gz` | `/stor/scratch/Wilke/if2453/bio-box-staging/GART/` | `3,204,867,642` | `d6f3675de4e6a97b475d213779f96f1b47c5cce8bf495e2ca575ad4bb31a80c5` |
| `GART2` | `genomes.json` | `/tmp/if2453/c63/GART2/` | `45,310,660` | `ffd935d5b7efb43e41409b4e10fa8f911dab1f1341da47d2d4163f2ff6793898` |
| `GART2` | `objects.json` | `/tmp/if2453/c63/GART2/` | `340` | `b051ce98b58751e0db38391ce220d2db68ee2ef1fd041ac604f3930eedf6ab8a` |
| `GART2` | `records.tsv` | `/tmp/if2453/c63/GART2/` | `120,697,642` | `220b524daacd155be53753095b9d53ae9a8458303788b5147ad348f653cd5ca5` |
| `GART2` | `repeated_records.tsv` | `/tmp/if2453/c63/GART2/` | `39` | `d109495ee76e4643c181eae183c452c9a84560f97a79bbe22d6630780368c565` |
| `GART2` | `shared_records.tsv` | `/tmp/if2453/c63/GART2/` | `711,340` | `b66a036233a54237b2821aad3f528888da877b52efef663a2e34a9e7491d1198` |
| `GART2` | `genetic_codes.tsv` | `/tmp/if2453/c63/GART2/` | `2,115,199` | `7fefc8d81eeb341ccc6efe69817ec714dc2e43844b4db95b6b6913f1345e8599` |

- **Dry checks**: `s2_check_inputs.py GART1` -> exit `0`, `problems: []` (`9,234` genomes, `5,218,293` rows, `codes: [11]`); `s2_check_inputs.py GART2` -> exit `0`, `problems: []` (`15,320` genomes, `1,213,095` rows, `codes: [11, 25]`).
- **`/tmp` usage**: `833 MB` total in `/tmp/if2453/c63/` (`< 1.5 GB`).

---

## 5. Estimated build cost (`--depth called --crispr`) and next step

*Derived* from `biobox.stage` (`INGEST_CPU_MIN = 4.11`, `INGEST_GZ_CPU_MIN = 0.26`, `CRISPR_CPU_S_PER_GBP = 182.0`, `CRISPR_CPU_S_PER_K_CONTIGS = 0.50`) and `biobox.l1.genecall` (`CALL_S_PER_MBP[("single", False)] = 0.7941`, `STEP_OVERHEAD = 1.23`):

| Stage | `ingest` (CPU-min) | `genecall` `single` (CPU-min / batches) | `crispr` default (CPU-min) | Total `--depth called --crispr` | Est. root storage |
|---|---:|---|---:|---:|---:|
| `GART1` (`38.962` new Gbp, `5.195M` new contigs) | `171.14` | `634.27` (`3` batches) | `152.45` | **`957.86` CPU-min** | `28.05–33.12 GB` |
| `GART2` (`11.666` new Gbp, `1.203M` new contigs) | `51.26` | `189.90` (`1` batch) | `42.71` | **`283.87` CPU-min** | `8.40–9.92 GB` |
| **Combined `GART1 + GART2`** | **`222.40`** | **`824.17`** | **`195.16`** | **`1,241.73` CPU-min** | **`36.45–43.04 GB`** |

**Part 2 (once `r26.11.3` clears `/tmp`)**: move `/tmp/if2453/c63/{GART1,GART2}/*` + `/stor/scratch/Wilke/if2453/bio-box-staging/GART/{gart1_globdb,gart2_globdb_cpr}.fna.gz` into `/tmp/if2453/s2/inputs/{GART1,GART2}/` (where `objects.json` already points), mirror to `tower:/mnt/user/research-bulk/bio-box/sources/{GART1,GART2}/`, and add `stages/GART1.conf` + `stages/GART2.conf` (`DEPTH=called`, `CRISPR=yes`).
