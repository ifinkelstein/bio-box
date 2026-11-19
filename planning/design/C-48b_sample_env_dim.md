# Design Proposal: `sample_env_dim` and `genome_env_link` (Card C-48b)

Status: **Proposed** (before `rS0.1` schema freeze, `AGENTS.md` §3.8; `src/biobox/registry/schemas.py` untouched by card C-48b).

## 1. Motivation (`AGENTS.md` §3.6)

Environmental and sample metadata (GOLD 5-level ecosystem path, `SPIRE` microntology, `MGnify` biome, physicochemical parameters such as pH/temperature/salinity/oxygen/depth/elevation/subsurface/lat/lon/Longhurst code, host/body site, and `redistributable` flag) can be updated or enriched without any nucleotide sequence changing. Per `AGENTS.md` §3.6, mutable sample/genome attributes must **never** be copied into fact tables (`contig_tapes`, `protein_occurrence`, `token_postings`). Instead, they belong in release-partitioned or sidecar dimension tables joined at query time through dimension views.

## 2. Proposed Tables

### 2.1 `sample_env_dim` (`layer="dim"`, `location="dim/sample_env_dim"`, `partition_keys=("release",)`, `small_ok=True`)

- **Primary key**: `(env_uid, release)` where `env_uid = blake3("env:" + source + ":" + source_id)` (16-byte `FIXED_LEN_BYTE_ARRAY(16)`, `AGENTS.md` §3.2).
- **Sort key**: `(source, source_id)` (min/max statistics on sort-key columns, `AGENTS.md` §3.7).
- **Columns**:
  - Identity & cross-references: `env_uid` (`binary(16)`), `source` (`string`: `img_m`, `spire`, `mgnify`, `gold`, `globdb`), `source_id` (`string`), `img_taxon_oid`, `gold_study_id`, `gold_biosample_id`, `gold_sp_id`, `gold_ap_id`, `biosample_acc`, `sra_runs`, `sra_experiments`, `insdc_assembly_acc`, `insdc_genbank_acc`.
  - Harmonized ecology: `eco_l1..eco_l5` (GOLD 5-level path), `ecosystem_path`, `harmonized_habitat` (enum: `host-associated`, `anthropogenic`, `terrestrial`, `aquatic`, `air`, `unlabelled`), `spire_microntology`, `mgnify_biome`, `habitat`.
  - Physicochemical & geographic: `ph` (`float64`), `temperature_c` (`float64`), `temperature_range`, `salinity`, `salinity_conc`, `oxygen_req`, `oxygen_conc`, `depth_m` (`float64`), `elevation_m` (`float64`), `subsurface_m` (`float64`), `latitude` (`float64`), `longitude` (`float64`), `longhurst_code`, `host_name`, `body_site`, `body_subsite`, `geo_location`.
  - Provenance & redistribution (`DEC-51`): `jgi_utilization_status`, `asm_bases` (`int64`), `redistributable` (`bool`, `false` for JGI/IMG-M/IMG-PR/mMGE/GEMs; `true` for GOLD metadata, SPIRE, MGnify, INSDC, GlobDB).
  - **PII exclusion (`AGENTS.md` §8)**: Contact/submitter/PI email columns (`STUDY CONTACTS`, `PROJECT CONTACTS`, `AP CONTACTS`) are stripped at ingest and never stored.

### 2.2 `genome_env_link` (`layer="dim"`, `location="dim/genome_env_link"`, `partition_keys=("release",)`, `small_ok=True`)

- **Primary key**: `(stage, genome_key, release)` (where `genome_key` is `genome_uid` hex for stored stages or source genome/record key for staged inputs).
- **Sort key**: `(stage, genome_key)`.
- **Columns**: `stage`, `genome_key`, `source_accession`, `link_method` (`asm_acc`, `img_oid`, `genbank_acc`, `globdb_env`, `source_provenance`), `env_uid` (`binary(16)`, nullable when provenance-only), `env_source`, `harmonized_habitat`, `ecosystem_path`, `redistributable` (`bool`).

## 3. Query Pattern (Zero Fact-Table Duplication)

Queries join `genome_dim` (or `contig_dim` → `genome_dim`) to `v_genome_env` (`genome_env_link LEFT JOIN sample_env_dim ON env_uid`), filtering by `harmonized_habitat`, `eco_l1..eco_l5`, `ph`, `temperature_c`, or `redistributable = true` before exporting or sharing results.

## 4. Open Decision for PI

Confirm whether `sample_env_dim` and `genome_env_link` should be promoted from the scratch sidecar (`/tmp/if2453/c48b/env_sidecar.duckdb`) into `src/biobox/registry/schemas.py` before `rS0.1` is signed, or remain a release-independent DuckDB sidecar until unbinned metagenome stages (`C-48e`) land.
