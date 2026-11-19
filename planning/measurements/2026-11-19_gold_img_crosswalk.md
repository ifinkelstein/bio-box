# JGI GOLD and IMG/M Environmental Metadata Crosswalk & DuckDB Sidecar (Card C-48b)

Date: 2026-11-19 (`wilkcomp02`, `bio-box-test_faec1f4.sif`, `OMP_NUM_THREADS=1`, DuckDB `threads=1`, `MemoryMax=16G`, `nice 19`).
Zero external network requests were issued under the 2026-11-19 source freeze (`AGENTS.md` §8). Every number carries an `AGENTS.md` §6 label (`*measured*`, `*derived*`, `*estimate*`, or `*published*`).

## 1. On-Disk Inputs and Source-Freeze Note

| Resource | On-disk path (`wilkcomp02`, 2026-11-19) | Status | Rows / Units (*measured*) |
|---|---|---|---|
| `GOLD` v.10 | `/tmp/if2453/s3/tables/gold/` (`goldData.xlsx`, `goldSraBiomeImg.xlsx`, `GOLDs5levelEcosystemClassificationPaths.xlsx`; also bundled in `/stor/work/Wilke/if2453/all-the-data/jgi_20261118/jgi_gold_and_img_m_metadata_20261118.tar`) | On disk | `247,513` Biosamples, `533,108` Organisms, `641,400` Sequencing Projects, `463,257` Analysis Projects, `7,126` SRA-biome APs, `2,426` 5-level ecosystem paths (PII/contact columns stripped at parse time) |
| `IMG/M` | `/tmp/if2453/s3/tables/img_m/img_m_*_{72,70}col.tsv.gz` (5 slices) | On disk | `79,524` Microbiomes (`126,357` total rows across all 5 slices) |
| `SPIRE` | `/tmp/if2453/s3/spire/tables/study_tsv/*.tsv` | On disk | `96,067` samples across `715` studies |
| `GlobDB` r226 | `/tmp/if2453/s3/globdb/dl/environmental_metadata/high_level_env__value.tsv` | On disk | Environmental counts for GlobDB r226 species reps |
| `MGnify` & `ENA`/`INSDC` C-43 tables | `/tmp/if2453/s3/tables/{mgnify,ena,ncbi_metagenomes}`, `/tmp/if2453/s3/c43_work/` | **Absent after 2026-11-17 `/tmp` reboot** | Not re-fetched under the 2026-11-19 source freeze; C-43 recorded totals from `planning/measurements/2026-11-16_run_crosswalk.md` are reported beside the on-disk SPIRE+IMG/M crosswalk below |

## 2. GOLD + IMG/M Identifier Linking (`79,524` IMG/M Microbiomes)

All counts *measured* from `/tmp/if2453/c48b/summary.json` (`scripts/stress/c48b_env_crosswalk.py build`, `wilkcomp02`, 2026-11-19):

- **Total IMG/M microbiome datasets**: **`79,524`** (*measured*).
- **GOLD cross-references**: `79,524` (`100.0%`, *derived*) carry `GOLD Analysis Project ID` (`Ga*`), `GOLD Sequencing Project ID` (`Gp*`), and `GOLD Study ID` (`Gs*`).
- **INSDC / SRA links (IMG/M + GOLD AP/SP combined)**:
  - `NCBI Biosample Accession`: **`56,130`** (`70.58%`, *derived*: `56,130 / 79,524`)
  - `SRA Experiment IDs` (`SRX`/`ERX`/`DRX`): **`44,950`** (`56.52%`, *derived*)
  - `SRA Run IDs` (`SRR`/`ERR`/`DRR`): **`6,954`** (`8.74%`, *derived*)
- **JGI Data Utilization Status (IMG/M enriched by GOLD AP/SP)**:
  - Explicit `Unrestricted`: **`50,699`** (`63.75%`, *derived*; `+329` over IMG/M alone via GOLD AP/SP enrichment)
  - Blank / non-JGI public (not marked `Restricted`): **`19,607`** (*derived*: `70,306 - 50,699`)
  - Total `non-Restricted` (`Unrestricted` + blank): **`70,306`** (`88.41%`, *derived*)
  - Explicit `Restricted`: **`9,218`** (`11.59%`, *derived*; excluded from the unrestricted take-list)
- **Overlap with on-disk `SPIRE` (`96,067` samples)**:
  - IMG/M microbiomes covered by `SPIRE` (matching BioSample or SRA run): **`18,061`** (`22.71%` of all IMG/M microbiomes, *derived*).
  - **Novel explicit-`Unrestricted` IMG/M metagenomes not in `SPIRE`**: **`36,580`** datasets holding **`17.307 Tbp`** (`17,306,870,652,752` bp, *measured*).
  - **Novel `non-Restricted` (`Unrestricted` + blank) IMG/M metagenomes not in `SPIRE`**: **`52,246`** datasets holding **`21.027 Tbp`** (`21,027,108,187,250` bp, *measured*).
  - Because `MGnify` and `ENA` run tables were absent on disk after the 2026-11-17 `/tmp` reboot, `in_mgnify = 0` and `in_insdc = 0` in this offline pass; `36,580` (`17.307 Tbp`) is therefore an **upper bound** on IMG/M metagenomes novel after both `SPIRE` and `MGnify` (to be tightened by passing `--c43-units` once `C-48e` is unheld).

## 3. Novel Unrestricted IMG/M Metagenomes Stratified by GOLD Harmonized Ecosystem

All counts *measured* from `/tmp/if2453/c48b/img_m_crosswalk.parquet` (`wilkcomp02`, 2026-11-19); Tbp figures *derived* (`bp / 1e12`). `redistributable = false` for all IMG/M entries (`DEC-51`).

| Harmonized habitat (GOLD 5-level L1/L2) | Total IMG/M | Explicit `Unrestricted` | `Restricted` | Covered by `SPIRE` | Novel explicit `Unrestricted` units | Novel explicit `Unrestricted` bases (Tbp) | Novel `non-Restricted` units | Novel `non-Restricted` bases (Tbp) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `terrestrial` (`Environmental; Terrestrial`) | `21,722` | `17,162` | `953` | `4,677` | **`13,184`** | **`7.828`** | `16,092` | `8.597` |
| `aquatic` (`Environmental; Aquatic`) | `31,829` | `23,805` | `1,302` | `8,667` | **`16,122`** | **`4.774`** | `21,860` | `6.271` |
| `unlabelled` (`Mixed, ...` / multi-L1) | `960` | `492` | `16` | `175` | **`321`** | **`2.695`** | `770` | `3.260` |
| `host-associated` (`Host-associated`) | `11,279` | `4,536` | `339` | `3,662` | **`3,009`** | **`1.338`** | `7,278` | `2.000` |
| `anthropogenic` (`Engineered`) | `13,617` | `4,646` | `6,608` | `846` | **`3,920`** | **`0.671`** | `6,163` | `0.897` |
| `air` (`Environmental; Air`) | `117` | `58` | `0` | `34` | **`24`** | **`0.001`** | `83` | `0.002` |
| **Total** | **`79,524`** | **`50,699`** | **`9,218`** | **`18,061`** | **`36,580`** | **`17.307`** | **`52,246`** | **`21.027`** |

Key finding: unlike `SPIRE` and `MGnify` (which are dominated by `host-associated` human/animal gut samples; `67,298` of `93,097` SPIRE-only runs and `38,141` of `56,898` MGnify-only runs in `2026-11-16_run_crosswalk.md`), **72.8% of novel unrestricted IMG/M assembled bases (`12.602 / 17.307 Tbp`, *derived*) are environmental `terrestrial` (`7.828 Tbp`) and `aquatic` (`4.774 Tbp`) metagenomes** — directly complementing `SPIRE`/`MGnify` for habitat-stratified tranches in `C-48e`.

## 4. Updated Four-Resource Non-Redundant Run Take-List

Order of preference (`2026-11-16_run_crosswalk.md` §4 + `DEC-51`): `SPIRE` → `MGnify` → `INSDC WGS` → `IMG/M` (4th because `redistributable = false`).

| Order | Resource | Taken (linked) | Taken (unlinked) | Assembled size stated | `redistributable` | Provenance of counts |
|---:|---|---:|---:|---|---|---|
| 1 | `SPIRE` | `85,971` samples (`122,062` runs) | `10,096` samples | `16 Tbp` (*published*, Schmidt 2024) | `true` | *measured*, `2026-11-16_run_crosswalk.md` §4 (96,067 samples re-verified on disk 2026-11-19) |
| 2 | `MGnify` | `37,628` assemblies (`44,177` runs) | `200` assemblies | `1.345 TB` gzipped contigs | `true` | *measured*, `2026-11-16_run_crosswalk.md` §4 (tables absent on disk post-reboot) |
| 3 | `INSDC WGS` | `2,306` assemblies (`2,922` runs) | `2,642` assemblies | `1.289 Tbp` (`550.3 + 739.2 Gbp`, *derived*) | `true` | *measured*, `2026-11-16_run_crosswalk.md` §4 (tables absent on disk post-reboot) |
| 4 | `IMG/M` (`Unrestricted` only) | `36,580` total novel vs SPIRE (`take` + `take_unlinked`; upper bound prior to MGnify/INSDC dedup) | included in `36,580` | `17.307 Tbp` (*measured*) | `false` | *measured*, `/tmp/if2453/c48b/img_m_crosswalk.parquet`, 2026-11-19 |

## 5. Scratch DuckDB/Parquet Environmental Sidecar (`/tmp/if2453/c48b/`)

Built in bounded memory (`OMP_NUM_THREADS=1`, DuckDB `threads=1`, `memory_limit='2GB'`) with deterministic PyArrow Parquet writes (`ROW_GROUP_SIZE = 122,880`, ZSTD, sorted by primary sort key):

| File | Bytes (*measured*) | Rows (*measured*) | SHA-256 (*measured*) |
|---|---:|---:|---|
| `env_sample_dim.parquet` | `10,274,665` | `222,424` (`126,357` `img_m` + `96,067` `spire`) | `7505928849c79e5b975e8e3279a5243be1f8805eeba6304a5b276f923ada106a` |
| `img_m_crosswalk.parquet` | `2,652,788` | `79,524` | `4c4c3835fc99ff1edcda09e1492e38811c5a77bf0330005be5a85f689cdafd84` |
| `genome_env_link.parquet` | `382,733` | `13,113` | `c05360a40cd742798b455d596c5bcbb80bca3e56b636c43fac521c28ab961647` |
| `env_sidecar.duckdb` | `274,432` | views: `env_sample_dim`, `img_m_crosswalk`, `genome_env_link`, `v_genome_env` | (DuckDB catalog over Parquet views) |

Attribute completeness inside `env_sample_dim.parquet` (*measured*, `wilkcomp02`, 2026-11-19):
- `img_m` (`126,357` rows, `redistributable = false`): `90,679` with latitude/longitude (`71.8%`, *derived*), `6,475` with collection temperature (`5.1%`), `5,052` with pH (`4.0%`).
- `spire` (`96,067` rows, `redistributable = true`): `74,053` with latitude/longitude (`77.1%`, *derived*).

## 6. Environmental Metadata Coverage Across Stored Stages and Staged Inputs

All counts *measured* from `/tmp/if2453/c48b/coverage.json` (`scripts/stress/c48b_env_crosswalk.py coverage`, read-only against `/tmp/if2453/r10/root` and `/tmp/if2453/s2/inputs`, `wilkcomp02`, 2026-11-19):

| Stage | Kind | Total units (*measured*) | Linked units (*measured*) | Coverage % (*derived*) | Link method |
|---|---|---:|---:|---:|---|
| `P1` | Stored prokaryote | `250` | `141` | `56.40%` | IMG/M+GOLD `assembly_accession` (`GCA`/`GCF`) |
| `P2` | Stored prokaryote | `2,247` | `1,331` | `59.23%` | IMG/M+GOLD `assembly_accession` (`GCA`/`GCF`) |
| `P3` | Stored prokaryote | `2,671` | `1,181` | `44.22%` | IMG/M+GOLD `assembly_accession` (`GCA`/`GCF`) |
| `P4` | Stored prokaryote | `2,329` | `2,059` | `88.41%` | IMG/M+GOLD `assembly_accession` (`GCA`/`GCF`) |
| `P5` | Stored prokaryote | `2,289` | `1,540` | `67.28%` | IMG/M+GOLD `assembly_accession` (`GCA`/`GCF`) |
| `P6` | Stored prokaryote | `2,180` | `828` | `37.98%` | IMG/M+GOLD `assembly_accession` (`GCA`/`GCF`) |
| **`P1..P6` subtotal** | Stored prokaryote | **`11,966`** | **`7,080`** | **`59.17%`** | IMG/M+GOLD `assembly_accession` |
| `S1` | Stored virus (`INPHARED`) | `27,258` | `3,859` | `14.16%` | IMG/M+GOLD `genbank_accession` |
| `PL1` | Stored plasmid (`PLSDB`) | `68,693` | `0` | `0.00%` | IMG/M+GOLD `genbank_accession` (replicon `NZ_*` IDs; host/env lives in PLSDB metadata) |
| `H1` | Stored virus (`HPGC`) | `7,228` | `0` | `0.00%` | Internal `huge_phage_*` IDs only in `contig_label` |
| `V1` | Staged virus (`UHGV` v1.0) | `168,536` | `168,536` | `100.00%` | Source provenance (`UHGV` human gut → `host-associated`) |
| `PL2` | Staged plasmid (`PlasmidScope` + `Yu 2024`) | `1,063,115` | `636,501` | `59.87%` | `IMG-PR` taxon OID (`568,150 / 699,971 = 81.17%`, *derived*) + `gut_plasmids_yu2024` provenance (`68,350`) + GenBank (`1`) |
| `G1` | Staged prokaryote (`GlobDB` r226) | `7,139` | `6,033` | `84.51%` | `GlobDB` r226 `high_level_env__value.tsv` |

See `planning/design/C-48b_sample_env_dim.md` for the proposed `sample_env_dim` + `genome_env_link` schema addition before `rS0.1`.
