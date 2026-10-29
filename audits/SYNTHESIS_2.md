# Round-2 audit synthesis (v2.0 → v3.0 decisions)

Date: 2026-10-28. Inputs: five independent audits, each of which first rewrote its own brief (what a top professional would add, rawest material, methods, outside proof) and then executed it.

| Report | Lens | Findings |
|---|---|---|
| [`critique_4_bioinformatics.md`](critique_4_bioinformatics.md) | Annotation practitioner | BX-01..BX-20 |
| [`critique_5_data_science.md`](critique_5_data_science.md) | Statistics / ML | DS-01..DS-17 |
| [`critique_6_architecture.md`](critique_6_architecture.md) | Systems architecture | AR-01..AR-16 + roadmap |
| [`critique_7_database.md`](critique_7_database.md) | Storage engines / DuckDB | DB-01..DB-20 |
| [`critique_8_data_quality.md`](critique_8_data_quality.md) | Provenance / data quality | DQ-01..DQ-19 + gates + golden set |

The product, container and CLI are renamed **`bio-box`** (formerly `koonin-box`).

## 1. What all five agreed on

1. **The v2.0 numbers are not a model of anything.** Contigs (2.5×10⁹ vs 5×10⁸), slots (3.6×10⁹ is 1.44 slots/contig, impossible under the ≥2-feature gate), windows (3.6×10⁹≫4 = 2.25×10⁸, not 2.8×10⁹), `window_postings` (1.8×10¹⁰ × 7.5 B = 135 GB, not 1.35 TB), Logan (hundreds of TB compressed, not 40.8 TB; and 40.8 TB at 400 MB/s is 28 h, not 68 h), GTDB (R226: 143,614 species / 732,475 genomes, not r220), and the 4.7-day budget all fail arithmetic or source checks. Every scale number now comes from one sizing model fed by measured pilot rates.
2. **Stage-0 dedup is wrong in three independent ways.** (a) First sightings are never written to RocksDB, so the second-sighting comparison has nothing to compare against. (b) `POPCNT(64-bit 1-bit MinHash XOR) ≤ 6` is ≈99.3% ANI with SD ≈3.8 bits; 95% ANI is E[popcnt] ≈ 22. (c) An exact hash over syncmers is not "90% identity invariant" (expected Jaccard ≈0.21 at 90% AA ID). The design is also order-dependent, so builds are not reproducible.
3. **Derived statistics break under "malleability".** Appends, remaps and compaction invalidate NR50 membership, subfamily IDs, pair statistics, entropies, q-values, pangenome planes and transition edges; `dense_contig_id` renumbering breaks sidecars. Stable content-derived IDs and an explicit staleness contract are required.
4. **The benchmarks are circular and unmeasurable as written.** DefenseFinder/PADLOC/TADB/AcrDB/Breaker models are inputs, so "rediscovering" CBASS or ykkC subtypes is trivial. "0.0%" and "100%" cannot be measured (0/200 still allows ≈1.5% at 95%).
5. **Build small first.** All five independently asked for a measured pilot before any planetary claim. This is now the organizing principle: see [`../STAGES.md`](../STAGES.md).

## 2. Decision record (applied to all plan files)

| ID | Decision | Source findings |
|---|---|---|
| D-01 | Scope: "public INSDC/SRA-derived assemblies plus named MGE databases at freeze date *D* after filters *F*". No "every locus on Earth". | DQ-10 |
| D-02 | One sizing model (`planning/sizing.yaml`, built from measured S0 rates). Every number in the docs is either *measured (stage, date)*, *derived (formula)*, or *estimate ±40%, pilot-pending*. Release-dependent widths (bitvector bytes = ⌈species/512⌉×64) are computed, not hard-coded. | BX-01, AR-02, AR-09, DB-01, DQ-09, DS-03 |
| D-03 | Hardware facts corrected. Current target: **wilkcomp02** (56 Haswell AVX2 cores, 1 TB RAM shared, CPU only, house rules). Future target: Hetzner AX162 (EPYC 9454P 48 c; base disks 2×1.92 TB, larger NVMe optional; 1 Gbps default, 10 G add-on) + GEX130 (128 GB RAM) for GPU stages. ProstT5 realistic 24–39 GPU-h at the 18M-segment scale, unmeasured. | AR-08, AR-02 |
| D-04 | Stage −1 identity resolution (BioSample/run/GCA↔GCF, then sequence digest, then skani containment ≥0.99 / ANI ≥99.5%) with source precedence. Stage 0b screening (FCS-adaptor + UniVec, phiX/spike-ins, FCS-GX, CheckM2, GUNC, geNomad). Stage-0 dedup becomes deterministic two-pass sort-merge with ≥1,024-bit FracMinHash/containment prefilter and skani verification; `contig_membership` ledger; one count per BioProject per vehicle. | BX-15, BX-16, BX-18, DQ-01..03, DB-04, AR-03, DS-02, DS-03, DS-14 |
| D-05 | Logan is the last stage, compute-to-data on AWS (us-east-1), accessions preselected by SRA metadata and STAT profiles. Never streamed to a lab box. | AR-01, BX-02 |
| D-06 | Taxonomy is explicit per contig: MMseqs2 taxonomy/CAT against the pinned GTDB release for prokaryotes; geNomad ICTV + iPHoP/spacer hosts for viruses and plasmids. Store `tax_rank_resolved`, `tax_confidence`, `tax_method_id`; low-confidence contigs never count toward order dispersion or decoy strata. | BX-03, DQ-04 |
| D-07 | Gene calling: single mode per genome + Codetta code inference (Tier 1); pyrodigal-gv (viruses); meta mode (fragments). Genetic-code host labels corrected (code 15 is phage; code 25 is SR1/BD1-5 lineages, not "all CPR"). smORFs 15–75 aa are *candidates* until promoted (cluster in ≥3 orders with RNAcode or dN/dS < 0.5). Split-gene merges require a slippery motif, Sec/Pyl machinery, or an intron/HEG call; otherwise `CANDIDATE_INDEL`. Every slot carries `evidence_level`. | BX-06..08, DQ-14 |
| D-08 | ncRNA: `cmsearch` with Infernal's own filter pipeline (no extra nhmmer prefilter); orphan CMs `cmcalibrate`d with their own thresholds and deduplicated against Rfam; IGR90 used only to choose what to scan; riboswitch pocket subtype assigned per member after `cmalign`. | BX-04, BX-05 |
| D-09 | Domains: Pfam / NCBIfam / PHROGs / arCOG with each library's own curated cutoffs and provenance; COG via eggNOG; CDD only via RPS-BLAST if used. DefenseFinder, PADLOC, CONJScan, TXSScan run as MacSyFinder *system calls* stored as sidecars, not as tokens. HMMER `Z` fixed and recorded. Nested/CP rules start from Pfam nesting and ECOD boundaries; one repeat-gap threshold (≤35 aa). | BX-09, BX-10, DQ-11 |
| D-10 | Structure: Foldseek-distributed 3Di only (AFDB, ESMAtlas30); transfer only at ≥90% identity to models with pLDDT ≥70 over ≥70% of the envelope; `fold_source` and confidence stored; AlphaFill pockets only for AFDB-model hits. ProstT5 runs on a GPU stage; CPU-only stages have no predicted 3Di. | BX-11, DQ-12, DS-15 |
| D-11 | Mobility regime from geNomad, CheckV, IntegronFinder 2, CONJScan/ICEscreen, ISEScan; `UNKNOWN_SHORT_CONTIG` and `regime_score` instead of defaulting to CORE. | BX-13 |
| D-12 | Catalytic motif stored as packed residues + `site_source` (M-CSA / CDD site / conservation-only); hash only for indexing; Hamming penalty only on curated sites; computed per member, not per representative; γ fitted, not chosen. | BX-12, DB-13, DQ-07, DS-08 |
| D-13 | Statistics: primary ranking is a Poisson-binomial score test `Z` plus a fitted log-odds effect `beta_hat ± se`; geometry, entropy and vehicle terms are separate columns, never subtracted ad hoc. Permutation null (not "target-decoy") with block length ≥25 ORFs, curveball content null, Efron empirical null + local fdr, Yekutieli hierarchical FDR across token levels with one pre-registered primary level. `powered` / `min_detectable_count` columns. Tier 2 uses a replicon-matched geometry test. Discovery unit is a *module*, not a pair. Weights fitted by a documented protocol into `params.yaml`; per-pair sufficient statistics stored so `bio-box rescore` can recompute without the base build. | DS-01, DS-06..08, DS-12, DS-14 |
| D-14 | Gene-tree homoplasy (`H_B|T_A`) is on-demand for top-N pairs, null in the base build. | BX-14 |
| D-15 | Pangenome planes are three-valued with Beta-binomial shrinkage, stored as compressed bitmaps over species IDs (dense only for top ~10⁵ tokens); one quality mask definition (closed genome or ≥3 independent genomes); Pagel-style correlated-evolution test replaces `hole & (hole>>1)`. | DS-13, DB-09, BX-16 |
| D-16 | `gLM-Domain` is a separate later project. When resumed: GraphPart-style balanced partitioning with a giant-component gate, vocabulary chosen on the training split only, gLM2 and frequency-prior baselines, context-shuffled ablation, calibration (ECE). | DS-05, DS-09, DS-17, AR-02 |
| D-17 | Minimum architecture: Parquet + DuckDB + a thin API. Custom indexes (Roaring, offset-indexed tape store) are built only when a measured SLO fails. If Index 1 is built: slot-level bitmaps queried with offsets (or overlapping ORF-index windows), never `orf_mod_8`. Latency targets are per query class, p50/p95, warm and cold, measured on the deployment host. | AR-06, DB-02, DB-03, DB-05..08 |
| D-18 | Schema fixes: badge re-laid out (sorting tags 8 bits; no overlap); `delta_bp` int32; `parent_seg` matches `seg` width; `regulon_locus_id` uint32; token IDs namespaced; `slot_annotations` side table for multi-library hits; provenance tables; enums mapped to Sequence Ontology terms. | DB-14..16, DQ-07, DQ-16 |
| D-19 | Identifiers: content-derived stable UIDs (sequence digest, accession digest, representative-UID clusters) plus release-local dense IDs; `release_id`, `stats_release_id`; `id_history` crosswalk (same / split / merged / retired). Sidecars declare `key_type` + `key_release`. | AR-04, DB-11, DQ-13 |
| D-20 | Malleability contract: deltas are raw-searchable immediately; statistics are "as of `stats_release_id`" with stale flags; `remap-token` supports merges of additive statistics only; parameter changes go through `rescore`. | AR-05, DB-10, DS-12 |
| D-21 | Catalog: read-only DuckDB catalog generated from a manifest and swapped by atomic rename; one writer process. | AR-11, DB-12 |
| D-22 | Operations: object storage is the system of record, RAID-0 is scratch; API binds localhost behind TLS + token; parameterized endpoints; any raw SQL is read-only with external access disabled. | AR-10, AR-12 |
| D-23 | Licensing: `sources.yaml` contract per source with a `redistribution` flag; public packs contain only cleared sources; non-redistributable sources ship as bring-your-own-data recipes. Known blockers: ISfinder (no download without authorization), JGI IMG policy, GSA-Human (do not ingest). | AR-13, DQ-05 |
| D-24 | Benchmarks: B0 engineering oracle harness (brute force equals index); B0-sim Zombi/ALF planted truth; B1 temporal holdout against a dated dictionary snapshot; golden-genome regression set; every criterion as a rate with an exact-binomial 95% CI; pre-registered, hash-committed trap panel curated independently; external baselines (STRING, Coinfinder, gLM2); a hash-timestamped prospective prediction registry. | DS-04, DS-10, DS-11, DS-16, DQ-10, BX-20 |
| D-25 | Orchestration: one workflow DAG (see `planning/data_layer_and_tests.md`) with content-addressed shards, per-shard manifests, deterministic outputs and kill-resume tests. | AR-07 |
| D-26 | Containers: CPU image (Apptainer on wilkcomp02, Docker on Hetzner) and a separate CUDA build image; model weights ship in the data pack, not the image. | AR-14 |
| D-27 | Data accretes in stages (S0 golden genomes → phages → plasmids → small genomes → GTDB reps → …), each with a gate. See [`../STAGES.md`](../STAGES.md). | All five |
| D-28 | Prior art updated: gLM2 / OMG corpus / Gaia (Tatta Bio), FESNov (Rodríguez del Río et al. 2024), Pavlopoulos et al. 2023, Durairaj et al. 2023, Coinfinder, MetaGraph, Logan-search, Bacformer, Evo 2, geNomad. | DS-16, prior-art gaps |

## 3. Disagreements and how they were resolved

| Issue | Positions | Resolution |
|---|---|---|
| Keep Index 1 (Roaring windows)? | AR: drop until needed. DB: fix with slot-level offsets. | Drop from the critical path (D-17); the DB design is kept as the fallback spec. |
| Dedup verification | DS: 1,024-bit sketch + skani. DQ: FracMinHash containment. DB: two-pass sort-merge. | All three combined (D-04): sort-merge structure, containment sketch prefilter, skani verify. |
| Contig count placeholder | DB: use 5×10⁸ until measured. AR: re-derive from a sample. | No placeholder survives: the sizing model's first input is the S0–S4 measured slots/contig. |
| `func_coupling_llr` | DS: not a likelihood ratio; replace. Round 1: primary metric. | Replaced by score-test `Z` + `beta_hat` (D-13); the round-1 components survive as diagnostics. |
