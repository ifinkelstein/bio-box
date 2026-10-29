# `all-the-data`: comparative genomics and synteny engine (`bio-box`) — master specification v3.0

> **Document suite**
> - [`00-spec.md`](00-spec.md): scope, targets and non-negotiable invariants (this file, v3.0).
> - [`STAGES.md`](STAGES.md): the staged accretion plan. Data is added one gated stage at a time, starting with small genomes on wilkcomp02.
> - [`PLAN.md`](PLAN.md): engineering plan (pipeline, schema, statistics, benchmarks).
> - [`DEPLOYMENT-AND-MALLEABILITY.md`](DEPLOYMENT-AND-MALLEABILITY.md): `bio-box` container, release and update model, deployment to wilkcomp02 and later Hetzner.
> - [`planning/sources_and_sizing.md`](planning/sources_and_sizing.md): source ladder, versions, licenses and the sizing model per stage.
> - [`planning/data_layer_and_tests.md`](planning/data_layer_and_tests.md): data layers (L0–L5), physical layout, DuckDB catalog and the test suite that gates each stage.
> - [`PRIOR-ART.md`](PRIOR-ART.md): existing systems, what each already covers, and what this project reuses.
> - [`planning/search/SEARCH_PLAN.md`](planning/search/SEARCH_PLAN.md): ad hoc search (recipes, runner, search pack, Z policy, search tiers ST0–ST5, S-tests). Decision IDs `SD-01`..`SD-21` refer to it.
> - [`planning/tests/TEST_ARCHITECTURE.md`](planning/tests/TEST_ARCHITECTURE.md): test organization (operation registry O-*, failure ranking F-01..F-15, metamorphic relations MR-01..MR-18, workflow manifests, micro-corpus and synthetic worlds, layout tests LY-01..LY-10, rings R0–R5, `test_run` table).
> - Audits: round 1 [`audits/SYNTHESIS.md`](audits/SYNTHESIS.md) (critiques 1–3); round 2 [`audits/SYNTHESIS_2.md`](audits/SYNTHESIS_2.md) with [`critique_4_bioinformatics.md`](audits/critique_4_bioinformatics.md), [`critique_5_data_science.md`](audits/critique_5_data_science.md), [`critique_6_architecture.md`](audits/critique_6_architecture.md), [`critique_7_database.md`](audits/critique_7_database.md), [`critique_8_data_quality.md`](audits/critique_8_data_quality.md). Decision IDs `D-01`..`D-28` below refer to `SYNTHESIS_2.md`.

---

## 1. Mission and scope

Build a self-hostable comparative genomics engine (`bio-box`) that indexes the gene neighbourhoods of prokaryotic, archaeal, viral and plasmid sequences, with protein domains, non-coding RNAs and mobility context on one coordinate tape. The engine answers questions such as "which uncharacterized domains sit next to family X, in how many independent lineages, outside mobile vehicles, and with what evidence that the neighbours are functionally coupled?"

**Scope statement (D-01).** The corpus is the set of public INSDC/SRA-derived assemblies and named MGE databases listed in `sources.yaml`, frozen at a stated date *D*, after the filters *F* defined in `PLAN.md` §2. It is not "every locus on Earth". Every release states *D*, *F*, the source versions and the counts that survived each filter.

**Representations carried on each slot of the tape:**
1. **Sequence-family hierarchy.** Clan → family (Pfam, NCBIfam, PHROGs, arCOG, COG via eggNOG, Rfam) → `subfamily_50_id` (50% identity over the domain envelope) → catalytic motif (packed residues at curated active-site columns, with `site_source`). Each hit keeps its library, model accession, library version and curated-cutoff status.
2. **Structure, where evidence exists (D-10).** Fold assignments (ECOD) come from Foldseek against distributed 3Di databases (AFDB, ESMAtlas30). 3Di is transferred only at ≥90% identity to models with pLDDT ≥70 over ≥70% of the envelope, and the slot records `fold_source` and confidence. ProstT5-predicted 3Di for dark segments is produced only in a GPU stage; CPU-only stages carry no predicted 3Di.
3. **Non-coding and recombination features.** Riboswitches and other structured RNAs (Infernal `cmsearch`, with ligand-pocket subtypes assigned per member), CRISPR arrays, `attC` sites, tRNAs and tmRNAs, terminators, introns and ribozymes.
4. **Biophysical badge.** TM helices, signal peptide class, sorting tags, pI, motif flags, length and repeat count, in a fixed non-overlapping bit layout (`PLAN.md` §3.3).

## 2. Deployment targets

| Target | Role | Constraints that shape the design |
|---|---|---|
| **wilkcomp02** (current) | All stages that fit a shared, CPU-only lab server | 56 physical Haswell cores (AVX2, no AVX-512); 40 threads by day and 56 at night; 1 TB RAM shared, with a 200 GB memory cgroup as the working limit; ~940 GB local `/tmp`; 30 TB NFS scratch shared over 10 GbE (no small-file fan-outs); no scheduler; a watchdog renices any process past ~369 CPU-minutes, so monolithic tools start at `nice 19`; Apptainer available. No GPU on wilkcomp02: ProstT5 and gLM stages are out of scope here. The same pod's wilkcomp03 (4× RTX 6000 Ada) is a GPU option by arrangement with the pod owners. |
| **Rented bare metal, e.g. Hetzner** (ultimate) | Stages that outgrow wilkcomp02; public serving | Hardware facts verified in `DEPLOYMENT-AND-MALLEABILITY.md` (AX162: EPYC 9454P, 48 cores; base disks 2×1.92 TB, larger NVMe optional; 1 Gbps uplink by default, 10 G as a paid add-on; GEX130 GPU node with 128 GB RAM). Deployment must be one container plus a data pack pulled from object storage. |
| **AWS us-east-1** (compute-to-data, last stage only) | Logan contig mining | Logan is hundreds of TB compressed and never streamed to a lab box (D-05). |

All latency, throughput and footprint numbers are measured on the target that will run them. Numbers in these documents are labelled *measured (stage, date)*, *derived (formula)* or *estimate ±40%, pilot-pending* (D-02).

---

## 3. Invariants (v3.0, seventeen)

1. **Data accretes in gated stages (D-27).** The corpus grows one stage at a time (golden reference genomes → complete phages → plasmids → complete prokaryotic genomes (S2C) → small genomes → GTDB species representatives → …), as defined in [`STAGES.md`](STAGES.md). A stage starts only when the previous stage passes its gate: ledger reconciliation, golden-genome regression, brute-force oracle equivalence, accretion equivalence (incremental equals rebuilt-from-scratch), and measured performance on the target host. The gate is a query over recorded test results (invariant 17). Earlier stages are never rewritten to admit a later one.

2. **Raw sequences stay available.** Every source contig's nucleotide sequence, including contigs collapsed into a representative by deduplication, is stored content-addressed by its UID. Protein sequences are stored the same way. Intergenic regions (IGRs) are stored as coordinates and sliced from the stored contig on demand, so their bases are never lost. All of them are retrievable by UID and by coordinate slice (`uid:start-end:strand`). UIDs are BLAKE3-128 content digests. Deduplication, clustering and representative selection add membership rows; they never discard a sequence. Any project may read this store directly without going through the derived tables.

3. **Identity resolution, screening and deterministic deduplication (D-04).**
   - Stage −1 links copies of the same biological sequence across sources: BioSample/run and GCA↔GCF links first, then exact sequence digest, then skani containment ≥0.99 with ANI ≥99.5%. Source precedence is RefSeq complete > GenBank complete > AllTheBacteria > GTDB draft > GlobDB/SPIRE MAG > Logan.
   - Stage 0b screens before gene calling: FCS-adaptor with UniVec, phiX and spike-in removal, FCS-GX, CheckM2 and GUNC on MAGs and SAGs, geNomad.
   - Contig deduplication is a two-pass sort-merge with a ≥1,024-bit FracMinHash containment prefilter and skani verification. Representative choice follows a fixed rule, so the same inputs give the same outputs regardless of input order. Every collapsed copy is a row in `contig_member`; one BioProject contributes at most one count per vehicle.

4. **Gene calling matched to the genome type (D-07).** Complete and HQ genomes use `pyrodigal` single mode per genome with Codetta genetic-code inference. Viruses use `pyrodigal-gv`. Fragments use meta mode and inherit the code from the contig's taxonomy. smORFs (15–75 aa) are *candidates* until promoted: a cluster present in ≥3 orders with RNAcode support or dN/dS < 0.5. Split-gene merges require positive evidence (slippery motif, Sec/Pyl machinery, intron or homing-endonuclease call); otherwise the pair is flagged `CANDIDATE_INDEL`. Every slot carries an `evidence_level`, and low-evidence slots are excluded from default statistics.

5. **Non-coding RNA search uses Infernal as designed (D-08).** `cmsearch` runs with Infernal's own filter pipeline. Orphan covariance models are `cmcalibrate`d, given their own thresholds and deduplicated against Rfam. IGR clustering only chooses what to scan; riboswitch pocket subtypes are assigned per member sequence after `cmalign`, because 90%-identity clusters erase 1–5-nt pocket differences.

6. **Sub-ORF domain segmentation with provenance (D-09, D-12).** Each protein is split into domain segments `(orf_idx, seg_idx)`, including nested hosts, circular permutations and collapsed repeat arrays. Nesting and permutation rules start from Pfam nesting annotations and ECOD boundaries; the gap rules are flagged heuristics. Each library is searched with its own curated cutoffs and one fixed HMMER `Z` per release line (`Z_prot = domZ = 1e8` sequences, SD-04 in [`planning/search/SEARCH_PLAN.md`](planning/search/SEARCH_PLAN.md)). System-level calls (DefenseFinder, CONJScan and TXSScan via MacSyFinder; PADLOC with its own HMM and system-rule engine) are stored as sidecars, not as tape tokens. Catalytic motifs are computed per member sequence, stored as residues with `site_source`, and only curated sites (M-CSA, CDD site annotations) drive pseudoenzyme penalties.

7. **Regulatory topology on the tape.** Junctions between adjacent slots carry `junction_topology_prev` (intra-ORF fusion, co-directional, divergent, convergent, sense–antisense nested, inter-locus break). `regulon_locus_id` (uint32) groups divergons and excludons. Riboswitch and leader attenuators never break `directon_id`. Terminator calls are a score with a clade-specific threshold validated against RegulonDB and term-seq, not a hard break.

8. **Explicit taxonomy, mobility and quality per contig (D-06, D-11).** Prokaryotic contigs get taxonomy from MMseqs2 taxonomy or CAT against the pinned GTDB release; viruses and plasmids get geNomad ICTV taxonomy and an iPHoP or spacer-based host. Each contig stores `tax_rank_resolved`, `tax_confidence` and `tax_method_id`; low-confidence contigs never count toward order-level dispersion or permutation strata. Mobility regime comes from geNomad, CheckV, IntegronFinder 2, CONJScan/ICEscreen and ISEScan, with `UNKNOWN_SHORT_CONTIG` and a `regime_score` instead of defaulting to core. The old `tier` is a derived view over separate `source_class`, `replicon_class`, `assembly_quality` and `molecule_type` columns. Contig-edge distances (`orfs_to_contig_start/end`) are recorded on every slot.

9. **Minimum storage architecture first (D-17, D-21).** The base system is hive-partitioned Parquet, a read-only DuckDB catalog generated from a manifest and swapped by atomic rename, and a thin API. Custom indexes (Roaring bitmaps, offset-indexed tape store, posting tables) are built only when a measured service-level objective fails on the target host. If a slot-proximity index is built, it uses slot-level positions queried with offsets, or overlapping ORF-index windows; never aligned windows that miss boundary pairs, and never modular ORF counters. Latency targets are per query class (p50/p95, warm and cold) and measured, not asserted.

10. **Stable identifiers across releases (D-19).** Every contig, protein, IGR and cluster has a content-derived UID (sequence digest, accession digest, or for clusters the BLAKE3 digest of the sorted member UIDs at creation plus the clustering-parameter hash, with the representative UID stored) next to any release-local dense ID. Every table carries `release_id`, and statistics carry `stats_release_id`. An `id_history` crosswalk records same / split / merged / retired between releases. Sidecars declare `key_type` and `key_release`, and registration refuses keys from another release.

11. **Statistics with a stated null and calibrated error rates (D-13, D-14, D-15).**
    - The primary pair ranking is a Poisson-binomial score test `Z` with a fitted log-odds effect `beta_hat ± se`. Geometry, neighbour entropy and vehicle terms are separate columns, not subtracted into one score.
    - Error control is a permutation null (block length ≥25 ORFs, plus a degree-preserving content null), an Efron empirical null with local fdr, and hierarchical FDR across token levels with one pre-registered primary level. Each pair carries `powered` and `min_detectable_count`.
    - Plasmid- and phage-confined pairs use a replicon-matched geometry test. The reported discovery unit is a module, not a pair.
    - Weights come from a documented fitting protocol stored in `params.yaml`. Per-pair sufficient statistics are stored so results can be rescored without a rebuild. Gene-tree homoplasy is computed on demand for top pairs.
    - Pangenome presence is three-valued with Beta-binomial shrinkage over species, and correlated gain/loss is tested with a Pagel-style likelihood-ratio test.

12. **Envelope- and active-site-gated closure.** An external hit (PDB, Swiss-Prot, PaperBLAST, Fitness Browser) is attributed to a domain segment only if it covers ≥70% of the segment envelope and, for PDB, ≥70% of the envelope maps to resolved `ATOM` residues. Closure is keyed by `(subfamily_50_id, catalytic motif)` and separates structural-genomics-only depositions from biochemically characterized proteins. Pairwise closure (co-crystal, shared publications, co-fitness) is indexed per pair.

13. **Malleability contract (D-20).** New data appended to a release is searchable as raw neighbourhoods immediately. Derived statistics are "as of `stats_release_id`" and flagged stale where an append, remap or new model affects them. `remap-token` handles merges of additive statistics only; splits require a release. Parameter changes go through `bio-box rescore` over stored sufficient statistics.

14. **Provenance, licensing and validation (D-23, D-24).**
    - Every source has a contract in `sources.yaml`: accession and version, retrieval date, checksum, license, and a `redistribution` flag. Public data packs contain only sources cleared for redistribution. Non-redistributable sources (for example ISfinder, restricted JGI data, GSA-Human, which is never ingested) ship as bring-your-own-data recipes. Every slot traces back to a source accession and coordinates, with the tool, version and database release that annotated it.
    - Acceptance tests are stated as rates with exact-binomial 95% confidence intervals. Never "0%" or "100%".
    - The tests include a brute-force oracle, a simulated planted-truth benchmark (Zombi/ALF), a temporal holdout against a dated dictionary snapshot (so curated defense/RNA models cannot leak the answers), a golden-genome regression set, a pre-registered hash-committed trap panel, external baselines (STRING, Coinfinder, gLM2), and a hash-timestamped prospective prediction registry.

15. **Searches are recorded recipe runs (SD-01..SD-05, SD-10, SD-12, SD-17..SD-20).**
    - Every ad hoc search is a `search_run` launched from a hashed, declarative recipe by one runner (`bio-box search run`). `serve` submits and polls; it never runs a search inline.
    - Each run records the recipe hash, argv, tool versions, target pack and shard sha256, query hash, Z, outputs, wall time, peak RSS, CPU-min per process, host and exit status. Seeds carry their source (paper and `file:line`, or database release and entry) and are resolved by sequence; model IDs are resolved against the release at submit time. A positive control is required and must be recovered.
    - E-values use one fixed search space per release line: `Z_prot = domZ = 1e8` sequences for protein search, `Z_nt = 1e5` Mb for Infernal and nhmmer. Every hit stores the bitscore, the E-value at that Z and the real target size. Thresholds are stored as bits.
    - Every hit joins to `protein_occurrence` (every slot, contig, coordinates, strand, genome, taxonomy, species-rep flag, tier), so iteration rounds and curation see neighbourhoods, not just sequences. Counts are reported raw, vehicle-independent, per species representative and per distinct arrangement.
    - A non-zero exit, timeout or OOM is recorded as `failed`. A failed run never reads as zero hits. Every run carries a `max_output_gb` cap and a `partial` policy for fragment genes (SD-20).
    - Outputs are capped at query time and schema-versioned: `--json` NDJSON, where a removed or renamed key needs a version bump (SD-18, S-13). Hit-set networks store a threshold table, so any cut is a lookup (SD-19). Candidates leave for the Domainator workbench as `--format domainator-gb`, with `evalue` at Z_ref and `score` in bits, and exported subsets are filtered on score (SD-17).

16. **Nucleotide search corpora are release artifacts (SD-06, SD-07).** Each release builds, once, a search pack that includes three nucleotide corpora next to the protein targets: IGRs with flank-gene metadata, dark ORFs (stop-to-stop ORFs ≥30 aa in six frames outside called genes, plus all frames in the 3 kb at each contig end) and CRISPR spacers with their array, repeat and contig. They are versioned, checksummed and keyed by stable UIDs like any other release table, so a search can go outside the gene calls without re-extracting the corpus. Whole-contig translated search is a scheduled batch job, enabled only at tiers where it fits. The pack's corpus-wide 3Di target (from ST3) is built only with the encoder chosen by the SD-21 bake-off.

17. **Tests are registered operations, manifests and recorded results ([`planning/tests/TEST_ARCHITECTURE.md`](planning/tests/TEST_ARCHITECTURE.md)).**
    - Every operation (O-INGEST through O-RUN) is registered with a contract (registered input and output schemas, keys, invariants), an oracle (brute force or an independent implementation) and at least one metamorphic relation (MR-01..MR-18). An operation missing any of these is not in the registry and cannot be used.
    - Every workflow (base pipeline step, query class, search recipe pattern) is a manifest over registered operations, with a positive and a negative control and a declared truth source. The harness generates the workflow's tests from the manifest; CI refuses a workflow that needs an unregistered operation.
    - Every test run writes `test_run` rows (test, ring, operation, workflow, params hash, outcome, metrics, code commit, image digest, `release_id`, seed). A stage gate is a query over `test_run` for the stage's release and image; a regression is a diff between runs. Blocking tests use fixed seeds and are never retried into passing.
    - Statistics require constructed truth: planted truth recovered in synthetic worlds, and null worlds with phylogeny and HGT but no coupling giving uniform p-values and FDR at nominal level (MR-18). Golden snapshots catch change, not error, and are not evidence of correctness.
    - Open layout choices (sequence-store format, row-group sizes, sort keys) are decided by row groups read for canonical queries (LY-06), not by latency; SLOs are checked on the chosen layout.

## 4. Deferred, not dropped

- **`gLM-Domain` (D-16)** is a separate project, started only after the tape is stable and a GPU target exists. When resumed it uses balanced homology partitioning with a giant-component gate (not single-linkage connected components). It also requires a vocabulary chosen on the training split only, gLM2 and frequency-prior baselines, a context-shuffled ablation and calibration reporting.
- **Logan (D-05)** is the last stage, run as compute-to-data on AWS with accessions preselected from SRA metadata and STAT taxonomy profiles.
- **Containers (D-26).** `bio-box` ships as a CPU image (Apptainer `.sif` on wilkcomp02, Docker on rented hosts) and a separate CUDA build image. Model weights and reference libraries ship in the data pack, not in the image.
