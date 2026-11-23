# Planning, Decisions, Design Notes, and Measurements (`planning/`)

This directory holds the **live governance, catalog roadmap, active architecture designs, task queue, and empirical measurement records** of `bio-box`. For the repository-wide entrypoint, see [`../DOCS.md`](../DOCS.md).

## 1. Authority & Live State

| File | Role |
|---|---|
| [`DECISIONS.md`](DECISIONS.md) | **Highest authority (`DEC-01`..`DEC-69`)**: append-only record of PI and tier-`H` architectural, operational, and catalog-scope decisions |
| [`CATALOG_GOAL.md`](CATALOG_GOAL.md) | **Corpus goal, build order, and deferred-data register (`DEC-59`, `DEC-60`, `DEC-69`)**: what is built, what is staged on disk, what selection filters left for later, and which external/niche sources are queued |
| [`tasks/STATUS.md`](tasks/STATUS.md) | **Live operational status**: active release/stage builds on `wilkcomp02`, merged card table, and open PI items |
| [`tasks/TASKS.md`](tasks/TASKS.md) | **Task cards (`A-*`, `B-*`, `C-*`, `D-*`, `E-*`, `F-*`, `W-*`)**: self-contained specifications for code, staging, measurement, and design work |

## 2. Architecture, Search, Sizing & Source Surveys

| Path | Role |
|---|---|
| [`data_layer_and_tests.md`](data_layer_and_tests.md) | L0–L5 table responsibilities, sequence store contracts, and plan test IDs (`C-`, `L-`, `R-`, `Q-`, `A-`, `D-`, `G-`, `K-`) |
| [`search/SEARCH_PLAN.md`](search/SEARCH_PLAN.md) | Custom-HMM search, search packs, locus neighbourhoods, tiers `ST0`–`ST5`, and `S-*` test contracts (supporting round notes in `search/01_*.md`..`07_*.md`) |
| [`tests/TEST_ARCHITECTURE.md`](tests/TEST_ARCHITECTURE.md) | Operation registry (`O-*`), failure modes (`F-*`), metamorphic relations (`MR-*`), synthetic/null worlds, and test rings `R0`–`R5` |
| [`sources_and_sizing.md`](sources_and_sizing.md) | `wilkcomp02` resource envelope, sequence-store byte rates, and stage sizing formulas |
| [`sources_metagenomic_survey.md`](sources_metagenomic_survey.md) | Survey of metagenomic MAG, viral, plasmid, and assembly collections |
| [`sources_from_s41467-026-79223-z.md`](sources_from_s41467-026-79223-z.md) | Cross-check of microbial/metagenomic databases named in recent literature against `bio-box` scope |
| [`footprint_request.md`](footprint_request.md), [`preflight_gaps.md`](preflight_gaps.md) | Shared-server footprint request and initial pre-mortem gap checklist |

## 3. Focused Design Notes (`planning/design/`)

| Design Note | Topic |
|---|---|
| [`design/C-38_import_depths.md`](design/C-38_import_depths.md) | Three stage import depths (`stored`, `called`, `annotated`), `--min-contig-len`, `domain_state`, and on-hit annotation rounds (`DEC-54`) |
| [`design/C-49_tiers.md`](design/C-49_tiers.md) | Storage tiers across `wilkcomp02` `/tmp`, `/stor/work`, and the `tower` NAS (`DEC-66`, `bio-box release retire`) |
| [`design/C-48b_sample_env_dim.md`](design/C-48b_sample_env_dim.md) | Joining GOLD, IMG/M, SPIRE, ENA, and GlobDB environmental/biome metadata into DuckDB dimension tables |
| [`design/B-37_split_clustering.md`](design/B-37_split_clustering.md) | Splitting MMseqs2 NR90/NR50 clustering into bounded `< 300` CPU-minute part processes (`DEC-29`, `DEC-38`) |
| [`design/B-47_inferred_domains.md`](design/B-47_inferred_domains.md) | Representative-based Pfam domain projection vs full scans and cross-stage measured-row copying (`DEC-35`, `DEC-62`) |
| [`design/C-21_shared_contigs.md`](design/C-21_shared_contigs.md) | First-stage contig ownership and multi-genome contig membership (`DEC-26`, `DEC-27`) |
| [`design/0-08_ncrna_models.md`](design/0-08_ncrna_models.md) | Curated Rfam 15.1 `prok-2` (`770` models) intergenic ncRNA scanning with clan competition (`DEC-28`, `DEC-30`) |
| [`design/W-25_flanks.md`](design/W-25_flanks.md), [`design/W-25_flank_tools_review.md`](design/W-25_flank_tools_review.md) | `bio-box nbr flanks` locus context grouping, UPGMA trees, and synteny figures |
| [`design/pLM_remote_system_search.md`](design/pLM_remote_system_search.md) | Design plan for protein language model (pLM) remote homology and multi-gene defense-system search |
| [`design/storage_decisions_draft.md`](design/storage_decisions_draft.md) | Working draft that led to `DEC-65`..`DEC-67` |

*(Queued tier-`H` design notes in [`tasks/TASKS.md`](tasks/TASKS.md): `planning/design/C-50_rna_viruses.md` for RNA viruses, Obelisks, and viroid-like circular RNAs; `planning/design/C-60_microbial_eukaryotes.md` for microbial eukaryotes and giant-virus hosts.)*

## 4. Empirical Measurements (`planning/measurements/`)

- Browse **[`measurements/README.md`](measurements/README.md)** for the complete topical index of stage builds, rolling releases, benchmarks, and scaling studies.
