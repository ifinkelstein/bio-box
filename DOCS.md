# bio-box Documentation Map

Single entrypoint to every document in the repository. All documentation lives directly in the repository tree (never inside ephemeral agent worktrees).

## 1. Start Here (By Role)

| Who you are | Read first |
|---|---|
| **Overview & current capabilities** | [`README.md`](README.md) |
| **Scientist running HMM searches & synteny analyses** | [`GUIDE.md`](GUIDE.md) and [`examples/README.md`](examples/README.md) ([`rnap`](examples/rnap/README.md), [`crispr_cas_neighbourhoods`](examples/crispr_cas_neighbourhoods/README.md)) |
| **LLM agent querying a published release** | [`AGENT_GUIDE.md`](AGENT_GUIDE.md) |
| **Operator building stages, cutting releases, or managing storage on `wilkcomp02`** | [`RUNBOOK.md`](RUNBOOK.md) and [`profiles/wilkcomp02.yaml`](profiles/wilkcomp02.yaml) |
| **PI / planner checking live build progress, catalog scope, or open questions** | [`planning/tasks/STATUS.md`](planning/tasks/STATUS.md), [`planning/CATALOG_GOAL.md`](planning/CATALOG_GOAL.md), and [`planning/DECISIONS.md`](planning/DECISIONS.md) |
| **Developer / agent implementing a task card** | [`AGENTS.md`](AGENTS.md) (`CLAUDE.md`) and [`planning/tasks/TASKS.md`](planning/tasks/TASKS.md) |

## 2. Corpus Scope, Data Provenance & Governance

| Path | What it holds |
|---|---|
| [`planning/DECISIONS.md`](planning/DECISIONS.md) | **Highest authority (`DEC-01`..`DEC-78`)**: append-only log of architectural, operational, storage, and source-scope decisions |
| [`planning/CATALOG_GOAL.md`](planning/CATALOG_GOAL.md) | **Catalog roadmap & deferred-data register**: built stages, staged waves, reversible selection/filter registers, and queued viral, long-read, single-cell, integrative-MGE, niche, organellar, circular-RNA/Obelisk, and microbial-eukaryote sources |
| [`sources.yaml`](sources.yaml) | Approved external data sources, URLs, licences, `redistributable` flags, and verified `sha256`/`md5` checksums |
| [`refdb.yaml`](refdb.yaml), [`golden.yaml`](golden.yaml), [`golden.lock.yaml`](golden.lock.yaml) | Frozen reference libraries (Pfam 38.2, Rfam 15.1) and the 15-genome golden test set |
| [`planning/tasks/STATUS.md`](planning/tasks/STATUS.md) | Running status of live builds, rolling releases, merged cards, and open items |
| [`planning/tasks/TASKS.md`](planning/tasks/TASKS.md) | Full specification of every completed and queued task card |

## 3. Active Architecture, Search & Design Documents

- **[`planning/README.md`](planning/README.md)** — full guide to `planning/`, including:
  - **Data model & layers**: [`planning/data_layer_and_tests.md`](planning/data_layer_and_tests.md), [`src/biobox/registry/schemas.py`](src/biobox/registry/schemas.py)
  - **Search & neighbourhood engine**: [`planning/search/SEARCH_PLAN.md`](planning/search/SEARCH_PLAN.md)
  - **Test rings, contracts & synthetic worlds**: [`planning/tests/TEST_ARCHITECTURE.md`](planning/tests/TEST_ARCHITECTURE.md)
  - **Resource envelope & sizing**: [`planning/sources_and_sizing.md`](planning/sources_and_sizing.md), [`scripts/stress/README.md`](scripts/stress/README.md)
  - **Source surveys**: [`planning/sources_metagenomic_survey.md`](planning/sources_metagenomic_survey.md), [`planning/sources_from_s41467-026-79223-z.md`](planning/sources_from_s41467-026-79223-z.md)
  - **Focused design notes (`planning/design/`)**:
    - Import depths & on-hit rounds: [`planning/design/C-38_import_depths.md`](planning/design/C-38_import_depths.md)
    - Pod `/tmp` + `/stor/work` + `tower` NAS storage tiers: [`planning/design/C-49_tiers.md`](planning/design/C-49_tiers.md)
    - Environmental metadata dimensions: [`planning/design/C-48b_sample_env_dim.md`](planning/design/C-48b_sample_env_dim.md)
    - Split clustering under process caps: [`planning/design/B-37_split_clustering.md`](planning/design/B-37_split_clustering.md)
    - Shared contigs across genomes/stages: [`planning/design/C-21_shared_contigs.md`](planning/design/C-21_shared_contigs.md)
    - Representative vs full Pfam domain scans: [`planning/design/B-47_inferred_domains.md`](planning/design/B-47_inferred_domains.md)
    - Curated intergenic ncRNA (`prok-2`): [`planning/design/0-08_ncrna_models.md`](planning/design/0-08_ncrna_models.md)
    - Locus flank clustering & synteny plots: [`planning/design/W-25_flanks.md`](planning/design/W-25_flanks.md), [`planning/design/W-25_flank_tools_review.md`](planning/design/W-25_flank_tools_review.md)
    - Protein language model (pLM) remote homology & system search: [`planning/design/pLM_remote_system_search.md`](planning/design/pLM_remote_system_search.md)

## 4. Empirical Measurements, Stage Reports & Release Notes

- **[`planning/measurements/README.md`](planning/measurements/README.md)** — topical index of every empirical benchmark, stage build report, and rolling release record, including:
  - **Rolling releases**: [`r26.11.3`](planning/measurements/2026-11-24_release_r26.11.3.md) (`119.90 Gbp`, `13` stages, live `--release latest`), [`r26.11.2`](planning/measurements/2026-11-21_release_r26.11.2.md) (`74.16 Gbp`, `11` stages, live beside it), [`r26.11.1`](planning/measurements/2026-11-16_release_r26.11.1.md) (`59.88 Gbp`, `9` stages, retired under `DEC-66` on 2026-11-25)
  - **Audits**: [`2026-11-25_audit.md`](planning/measurements/2026-11-25_audit.md) (2026-11-22 to 25), [`2026-11-20_audit_of_1117_1120.md`](planning/measurements/2026-11-20_audit_of_1117_1120.md)
  - **Recent built stages**: [`SPX1`](planning/measurements/2026-11-23_stage_SPX1.md) (`35.69 Gbp` SPIRE metagenomes), [`G1`](planning/measurements/2026-11-23_stage_G1.md) (`10.04 Gbp` GlobDB archaea), [`V1`](planning/measurements/2026-11-20_stage_V1.md) (`5.18 Gbp` UHGV phages), [`C5`](planning/measurements/2026-11-20_stage_C5.md) (`9.11 Gbp` complete prokaryotes), [`PL1`](planning/measurements/2026-11-15_stage_PL1.md) (`5.92 Gbp` PLSDB), [`H1`](planning/measurements/2026-11-15_stage_H1.md) (`1.82 Gbp` huge phages)
  - **Storage & harvest inventories**: [`2026-11-20_tower_relay.md`](planning/measurements/2026-11-20_tower_relay.md), [`2026-11-18_jgi_harvest_inventory.md`](planning/measurements/2026-11-18_jgi_harvest_inventory.md)

## 5. Slides & Worked Examples

- **Presentations**: [`presentations/README.md`](presentations/README.md) (`all-the-data-deck.html` and `all-the-data-deck.pptx`)
- **Worked examples (also run as release gates)**: [`examples/README.md`](examples/README.md)

## 6. Archived Initial Specifications & Audits (`docs/archive/`)

Foundational October–November 2026 design documents and pre-build critiques moved out of the repository root to keep live guides front-and-center; superseded sections yield to [`planning/DECISIONS.md`](planning/DECISIONS.md):

- [`docs/archive/README.md`](docs/archive/README.md) — overview of archived files: [`00-spec.md`](docs/archive/00-spec.md), [`PLAN.md`](docs/archive/PLAN.md), [`STAGES.md`](docs/archive/STAGES.md), [`DEPLOYMENT-AND-MALLEABILITY.md`](docs/archive/DEPLOYMENT-AND-MALLEABILITY.md), [`PRIOR-ART.md`](docs/archive/PRIOR-ART.md), and [`docs/archive/audits/`](docs/archive/audits/SYNTHESIS_2.md).
