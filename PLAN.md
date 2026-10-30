# `all-the-data`: Engineering, Ingestion, Data-Layer and Validation Plan (v3.0, after round-2 audit)

> **Provenance.**
> - v2.0 combined the round-1 audits ([`critique_1_biology.md`](audits/critique_1_biology.md), [`critique_2_systems.md`](audits/critique_2_systems.md), [`critique_3_ml_stats.md`](audits/critique_3_ml_stats.md), [`SYNTHESIS.md`](audits/SYNTHESIS.md)).
> - v3.0 applies the round-2 decision record [`audits/SYNTHESIS_2.md`](audits/SYNTHESIS_2.md) (D-01..D-28), which draws on five reports:
>   - [`critique_4_bioinformatics.md`](audits/critique_4_bioinformatics.md) (BX-)
>   - [`critique_5_data_science.md`](audits/critique_5_data_science.md) (DS-)
>   - [`critique_6_architecture.md`](audits/critique_6_architecture.md) (AR-)
>   - [`critique_7_database.md`](audits/critique_7_database.md) (DB-)
>   - [`critique_8_data_quality.md`](audits/critique_8_data_quality.md) (DQ-)
>
> Inline tags such as `(DB-04)` name the finding that changed a section.
>
> The tool stack, container and CLI are named **`bio-box`**.

---

## 0. Scope, staging and acceptance gates

### 0.1 Scope (D-01, DQ-10)
`all-the-data` indexes **public INSDC/SRA-derived prokaryotic, archaeal, viral and plasmid assemblies, plus the named MGE databases in §2.2, as of a freeze date *D*, after the filters *F* in §2.3**. Each release records *D* and *F*. The plan makes no claim to cover "every locus on Earth". It covers what those sources contain after those filters, and records what was rejected and why.

### 0.2 Where numbers come from (D-02)
Every scale number in this document carries one of three labels:

| Label | Meaning |
|---|---|
| **measured** (stage, date) | Taken from a ledger or benchmark run on named hardware. |
| **derived** (formula) | Computed from measured or published inputs; the formula is shown. |
| **estimate ±40%, pilot-pending** | Placeholder until the stage that measures it has run. |

The single source for all scale numbers is `planning/sizing.yaml`. It is regenerated from stage ledgers and checked in CI, and documents cite it rather than restating numbers. Release-dependent widths are computed, not hard-coded. Example: the bitvector row width is ⌈n_species / 512⌉ × 64 B, which is 17,984 B for GTDB R226 (143,614 species) (BX-01).

Round-2 corrections that every later section respects:

| v2.0 claim | Correction | Source |
|---|---|---|
| 2.5×10⁹ contigs, 3.6×10⁹ slots | 1.44 slots per contig contradicts the ≥2-feature gate. Both numbers are withdrawn until measured in S0–S4. | DB-01, AR-09 |
| 2.8×10⁹ windows | 3.6×10⁹ ≫ 4 = 2.25×10⁸. Per-contig padding changes it further. Moot under D-17. | DB-01 |
| `window_postings` 1.35 TB | 1.8×10¹⁰ × 7.5 B = 135 GB. At a realistic 15–25 B/row it is 0.3–0.5 TB. The row count itself is unmeasured. | DB-06 |
| Logan 40.8 TB, 27.3M runs | Hundreds of TB compressed (Logan v1.2, about 38M accessions; verify at freeze). The plan's own arithmetic was also wrong: 40.8 TB at 400 MB/s is 28 h, not 68 h. | BX-01, AR-01 |
| GTDB r220 (113,104 species / 596,859 genomes) | GTDB R226: 143,614 species / 732,475 genomes. Pinned per release in `sources.yaml`. | BX-01, DQ-06 |
| 8.5×10⁹ raw ORFs | Tier 1 alone is about 7–10×10⁹ ORFs before dedup (estimate ±40%). | BX-19, DQ-09 |
| "Ten invariants" | v2.0 actually listed twelve; the v3.0 spec listed fourteen and, with the search invariants 15–16 and the test invariant 17, now lists seventeen. | BX-19, SEARCH_PLAN §8, TEST_ARCHITECTURE |
| 4.7-day end-to-end build | Withdrawn. Replaced by the capacity model in §1.3. | AR-02 |

### 0.3 Staged accretion (D-27)
The build accretes stage by stage, smallest data first. The order (authoritative in [`STAGES.md`](STAGES.md)) is: S−1 preflight; S0 golden genomes; S1 complete phages; S2 plasmids; S2C complete prokaryotic genomes (GTDB R226 `Complete Genome` assemblies, ≤20 per species; SEARCH_PLAN §6.1); S3 small-genome bacteria and archaea; S4 IMG/VR HQ vOTUs (local-only); S5 all GTDB species representatives; S6 all GTDB genomes plus the RefSeq complete flag; S7+ Tier 3 and Logan off-box. Each stage adds data partitions without rewriting earlier ones. Each stage must pass its gate before the next begins.

- Stage ladder, gates and per-stage deliverables: [`STAGES.md`](STAGES.md)
- Sources, counts, licenses and per-stage compute: [`planning/sources_and_sizing.md`](planning/sources_and_sizing.md)
- Data-layer layout, identifiers and the per-layer test suite: [`planning/data_layer_and_tests.md`](planning/data_layer_and_tests.md)

### 0.4 Acceptance gates common to every stage
Tests are organized per [`planning/tests/TEST_ARCHITECTURE.md`](planning/tests/TEST_ARCHITECTURE.md) (authoritative): registered operations (O-*) with a contract, an oracle and metamorphic relations (MR-01..MR-19); workflows as manifests over those operations; test rings R0–R5. Every test run writes rows to the `test_run` table (test, ring, operation, workflow, params hash, outcome, metrics, code commit, image digest, `release_id`, seed). **A gate is a SQL query over `test_run`:** all blocking tests for the stage pass on the stage's `release_id` and image digest. A regression is a diff between runs. Blocking tests use fixed seeds and are never retried into passing. Items 1–8 keep their numbering; TEST_ARCHITECTURE §9 maps their test families onto operations and rings.

1. **Ledger balances.** For every object, `n_in = kept + merged + Σ rejected_by_reason`, reconciled against the source manifest (DQ-08).
2. **Determinism.** The same inputs in shuffled order give identical BLAKE3 hashes for every L0, L1 and L3 table (AR-07, DB-04).
3. **Accretion equivalence.** Building stage *N* on top of stage *N−1* gives the same content hashes as building stages 1..*N* from scratch. L2 cluster tables must match up to the published crosswalk (D-19).
4. **Oracle agreement.** Every query class returns the same answer as a brute-force scan of the unnested slot table (DB-02, DB-03).
5. **Golden-genome checks.** The S0 golden set passes its curated-truth checks (§5, B-G; G-01, G-02, G-04). Golden snapshot hashes detect change, not error; they are not correctness evidence (item 11 is).
6. **Sequence round trip.** Retrieved sequence digests equal the source checksums (§4.0, §5 B-L0).
7. **SLOs.** Every query class meets its p50/p95 targets on the deployment host, warm and cold (§1.2). SLOs are checked on the layout chosen by item 9; they do not choose it.
8. **Search tests.** The stage's search tier (ST0–ST5, §4.7) passes the S-family tests marked as gates (S-01, S-03..S-09, S-02d, S-12, S-13; S-11 from ST3) and its replays, S-02 is reported, the S-10 human dry run is done once per tier from ST1, and the search cost table is measured on the pod (SEARCH_PLAN §7; per-stage list in `planning/data_layer_and_tests.md` §4.10).
9. **Layout tests (LY-01..LY-10).** Every table of the stage's release passes: LY-01 registered schema (16-byte UIDs as `FIXED_LEN_BYTE_ARRAY`, nullability, enum closure); LY-02 declared sort order holds per file; LY-03 min/max statistics on every sort-key column in every row group; LY-04 file and row-group sizes within the `data_layer_and_tests.md` §3.2 bands; LY-05 partition keys in paths equal column values; LY-06 canonical point queries (one contig slice, one neighbourhood, one anchor) read at most N row groups, counted from `EXPLAIN ANALYZE`/profiling output; LY-07 DuckDB, pyarrow and polars read identical rows; LY-08 byte-identical rewrite; LY-09 schema evolution; LY-10 explicit catalog file list with sha256. Open layout choices (K-13 Parquet chunks vs bgzip, row-group sizes, sort keys) are decided at S0 by LY-06 row-group counts on the micro-corpus, not by latency.
10. **Operation and workflow coverage.** Every operation the stage uses is in the registry with contract, oracle and ≥1 MR; its MRs pass on the micro-corpus. Every workflow the stage ships has a manifest whose generated tests pass (schema, MRs, positive and negative controls, truth recovery within bands, S-06 budget, S-13 JSON schema). The coverage matrix (TEST_ARCHITECTURE §5.3) has no gap for the stage.
11. **Constructed truth.** Planted truth in the synthetic worlds is recovered at rates within its bands, and null worlds (phylogeny and HGT, no coupling) give uniform p-values and FDR at nominal level (MR-18; §5 B-NULL) for every statistic the stage reports.
12. **Tests of the tests.** Every S1-severity failure (F-01..F-10) has a test that fails when the failure is injected deliberately. Mutation score on the operation library ≥80% of mutants killed by R1–R2 (*target*; confirmed after the first measurement at S0).

Rings (TEST_ARCHITECTURE §8):

| Ring | Content | Data | Cadence | Gate role |
|---|---|---|---|---|
| R0 static | schema registry, recipe and manifest lint; registry completeness; image manifest | none | every commit | green on the gate commit |
| R1 unit and property | operation functions; Hypothesis properties; stateful accretion machine | generated | every commit | green on the gate commit |
| R2 micro-corpus | full build, oracles, all MRs, manifests, LY tests, fault injection; mutation testing weekly | micro-corpus (~50 synthetic + 3 real genomes) | every merge | green on the gate commit |
| R3 golden | S0 golden set with full model sets; G snapshots; S-02d | golden set | nightly and every image build | passes on the stage's image |
| R4 stage gate | sampled oracles, accretion, determinism, K, S, SLOs | stage data | per stage | items 1–10 |
| R5 calibration and benchmark | null worlds at stage scale, B1–B6, replays, scaling fit | synthetic worlds at scale; real | per tier | items 8, 11; §5 |

---

## 1. Hardware targets, workload SLOs and capacity model

### 1.1 Hardware (D-03, AR-08)

**Current target: `wilkcomp02` (UT Austin BRCF Wilke pod), CPU only.**
Specs and house rules were measured 2026-09-30. Sources: `/Users/ilya/Work/20 manuscripts/Generative AI and IP/research/homology/WILKECOMP.org`, and the watchdog notes recorded 2026-10-15.

| Item | Value | Consequence for this plan |
|---|---|---|
| CPU | 4× Xeon E7-4850 v3 (Haswell, 2.2 GHz), 56 physical / 112 logical cores, AVX2, no AVX-512, 4 NUMA nodes | Use the AVX2 code paths. Measure per-thread rates on this machine; never extrapolate from modern cores (roughly 3× slower per thread). |
| Thread budget | 40 threads by day, 56 at night or weekends; never hyperthreads | Every stage sets `--threads` explicitly. |
| RAM | 1,007 GB physical; house rule <20% by day, <65% overnight | Default cap: a `systemd-run --user --scope -p MemoryMax=200G` cgroup per job. |
| Local disk | `/tmp`, about 940 GB free, ext4, volatile | All tmp space and indexes live here during a run. Results go back with one rsync per chunk. |
| Shared storage | `/stor/scratch/Wilke` (30 TB, not backed up); `/stor/work/Wilke` (75 TB, backed up) over shared NFS on 10 GbE | Few large files only. Never fan out per-genome files on NFS. |
| Watchdog | Renices any single process past about 369 CPU-minutes | Start monolithic tools (MMseqs2, Foldseek, HMMER bulk scans) at `nice 19`. Run per-genome worker fan-outs at `nice 15`. |
| Scheduler | None. Singularity/Apptainer is present. Docker group use is unconfirmed. | The workflow manager runs with the local executor and resource caps. A STOP file pauses the run between chunks. |
| GPU | None on wilkcomp02. wilkcomp03 has 4× RTX 6000 Ada but is out of scope for now. | No ProstT5 stage on this target. 3Di comes only from transfer (D-10). |

**Future target: rented bare metal (Hetzner).** Facts below were corrected in AR-08; re-verify prices and configurations at order time.

| Node | Verified facts | Notes |
|---|---|---|
| AX162-R / AX162-S | AMD EPYC 9454P, 48 cores / 96 threads, Zen 4 with AVX-512. Base storage 2×1.92 TB NVMe; larger NVMe is a priced option. | v2.0 named an EPYC 9654 and 2×7.68 TB as standard. Both were wrong. |
| GEX130 (GPU) | 1× RTX 6000 Ada (48 GB), 128 GB RAM | v2.0 said 256 GB RAM. |
| Network | 1 Gbps default uplink; 10 G is a paid add-on | At 1 Gbps a 3.45 TB pack takes 7.7 h (derived: 3.45×10¹²×8 / 10⁹ s). At 10 Gbps line rate it takes 46 min. |

RTX 6000 Ada: 91 TFLOPS is the FP32 non-tensor peak (AR-02). ProstT5 encoder throughput on 1.8×10⁷ segments is an **estimate of 24–39 GPU-h, pilot-pending**. It assumes about 2.4 GFLOP per residue at 30–50% of the dense FP16 tensor peak. v2.0's 14.3 GPU-h implied more than 80% utilization.

### 1.2 Workload model and SLOs (DB workload model, D-17)
SLOs are stated per query class with explicit output bounds. They are measured on the deployment host under 8 concurrent clients. Tokens are drawn stratified by document frequency: hub (>10⁶), mid (10³–10⁶) and tail (<10³). Targets apply at the stage where the class first exists. A stage fails its gate if a target is missed; the fix is either an index (§4.6) or a revised target with a written reason.

| Class | Example | Served by (minimum architecture) | Output bound | Target warm p50 / p95 | Target cold p95 |
|---|---|---|---|---|---|
| Q1 count | Number of loci with A and B within ±k ORFs, filtered by tier or phylum | DuckDB over `slot_pairs` (§4.2) | scalar | 0.5 / 2 s | 10 s |
| Q2 neighbourhoods | Fetch 1,000 contig tapes | DuckDB over `contig_tapes` sorted by `dense_contig_id` | 1,000 contigs | 0.3 / 1 s | 5 s |
| Q3 anchor profile | Partners of A ranked by `z_score` | `pair_stats` (§4.3) | top 500 | 0.1 / 0.5 s | 3 s |
| Q4 pair drill-down | Loci of (A, B) with geometry | `slot_pairs` | ≤10⁴ rows, sampled | 0.5 / 2 s | 10 s |
| Q5 NOGD hole | `core(A)&core(B)&~pan(C)&masks` | `pangenome_planes` (§4.4) | species set | 0.2 / 1 s | 5 s |
| Q6 co-transmission screen | Top-k X for C over 10⁵ candidates | `transition_edges` | top 100 | 1 / 5 s | 20 s |
| Q7 sequence slice | Retrieve (contig_uid, start, end, strand) | Sequence store (§4.0) | ≤100 kb | 0.05 / 0.2 s | 1 s |
| Q8 ad-hoc SQL | Analytic scans | DuckDB, read-only | none | no SLO; memory-capped and cancellable | none |

**Asynchronous search job classes (SD-01, §4.7).** These are recorded `search_run` jobs, not serve queries. They are **not subject to serve timeouts or the p50/p95 targets above**; `serve` only submits and polls. Their per-query costs are budgets on the search lane (8 threads by day), gated by S-06. Codes, servers and full notes are in SEARCH_PLAN §5; all are seconds at S0.

| Class | Served by | S5 budget (derived) | S6 budget (derived) |
|---|---|---|---|
| Q-FAM | stored library hits via `protein_occurrence`; new HMM by `scan-new-models` on NR90 + projection | seconds (stored); ~0.5 core-h per new HMM | ~1 core-h per HMM |
| Q-ITER | recipe runner + iteration controller on NR50 | jackhmmer ×5 rounds ~5–10 min at 40 threads, ~30–50 min on the lane | ~10–20 min at 40 threads |
| Q-REMOTE | HH-suite3 vs Pfam-HHM, PDB70, hit-local HH DB; HHpred external | seconds–minutes | same |
| Q-STRUCT | external AF3/ColabFold/DALI; local 3Di of NR50 from ST3 | ~11 min per query vs 8e7 reps at 40 threads | same |
| Q-ARCH | `slot_pairs` (`junction_topology = 0`); `segments.faa` re-search | seconds (tokens) | seconds |
| Q-NBR | `bio-box nbr`; `nbr score` positive-set enrichment | seconds per anchor set; full-corpus score minutes–1 h (estimate) | Hetzner for full scans |
| Q-PROF | presence bitvectors, pangenome planes; user tokens on demand | seconds | seconds |
| Q-RNA | cmsearch/nhmmer on `igr.fna`; targeted CMfinder | ~12 min per CM over all IGRs at 40 threads | ~1.5 h |
| Q-NT | `dark_orf.faa`; `bio-box region translate`; spacer → protospacer | minutes (dark ORFs) | DIAMOND blastx over all contigs ~55 h, Hetzner or spot |
| Q-MOTIF, Q-SSN, Q-CLUST, Q-SYS | per-hit profile columns; all-vs-all in a hit set; stored clusters; system sidecars | seconds; minutes for ≤1e5 (SSN); stored | same |

These targets replace v2.0's "<4 ms / <12 ms / <10 ms" claims. Those relied on `rowid IN (...)` pushdown into Parquet, which DuckDB does not document, and on hub anchors of 5–80 MB that are really 6–12 GB (DB-05, DB-06). Full enumeration of hub results is streamed through a cursor and has no latency target (DB-07).

### 1.3 Capacity model (AR-02)
Each row is a per-unit rate multiplied by a number of units. Every rate stays unmeasured until the S0 microbenchmarks on wilkcomp02, which use 2.2 GHz Haswell with AVX2. Unit counts are derived per stage in `planning/sources_and_sizing.md`. Wall time = Σ(rate × units) / threads × 2 (contingency). Rows marked "unbudgeted in v2.0" had no budget line at all.

| Stage | Work unit | Rate | Units (per stage) | Notes |
|---|---|---|---|---|
| Identity resolution + screening (§2.3) | Mbp screened (FCS-adaptor, FCS-GX, CheckM2, GUNC, geNomad) | unmeasured (S0) | Mbp of new input | FCS-GX deferred (needs ~470 GB in RAM); accreted later as a QC sidecar (§2.3). |
| Gene calling | Mbp per thread (pyrodigal single, -gv, meta; × number of codes tested) | **measured S0, 2026-10-30:** 1.2 Mbp/s per Haswell thread (golden set, 43 Mbp in 36.7 s, policy p2, from the image); 3.3 Mbp/s on Apple Silicon (`planning/measurements/2026-10-30_wilkcomp02.md`) | Mbp kept | Codetta per Tier 1 genome is an extra cost. At 1.2 Mbp/s, S5 (400–450 Gbp) is ~95–105 core-h of gene calling (*derived*). |
| smORF promotion | smORF clusters × RNAcode | unmeasured (S0) | candidate smORFs | Unbudgeted in v2.0 (BX-07). |
| Protein clustering | `mmseqs linclust` / `clusterupdate` core-h per 10⁹ sequences | unmeasured (S0) | proteins | Linclust scaled from its paper is about 1,490 core-h per pass at 8.5×10⁹ sequences (derived, AR-02). v2.0 budgeted 105 CPU-h. |
| Domain scan | pyhmmer cells/s/thread (Pfam, NCBIfam, PHROGs, arCOG) | unmeasured (S0) | Exact-unique protein residues × models through S2C; NR90 representative residues × models from S3 on (SD-05) | From S3, hits are projected to members through `cluster_member_aln`; members whose projected bitscore is within 10 bits of the cutoff, or whose alignment covers <90% of the domain envelope, are re-scored exactly (test S-08). Fixed Z recorded: `-Z 1e8 --domZ 1e8` for every protein scan in the release line (SD-04, D-09). This replaces the three earlier statements (per genome, NR50 residues, NR90 reps). |
| System calls | Genomes per thread-hour (MacSyFinder for DefenseFinder/CONJScan/TXSScan; PADLOC separately) | unmeasured (S0) | genomes | Results go to a sidecar (D-09). |
| ncRNA | `cmsearch` Mbp/thread-hour, Rfam prokaryotic subset + orphans | unmeasured (S0) | IGR90 Mbp | No extra nhmmer prefilter (D-08). Fixed Z recorded: `-Z 1e5` (Mb) for every cmsearch and nhmmer run in the release line (SD-04). |
| 3Di transfer | `mmseqs search` vs AFDB/ESMAtlas30 | unmeasured (S0) | NR50 reps | Applied only at ≥90% ID (D-10). |
| ProstT5 | Residues/s on GPU | unmeasured (GPU stage) | dark NR50 residues | Off wilkcomp02. Estimate: 24–39 GPU-h at 1.8×10⁷ segments. |
| Taxonomy | `mmseqs taxonomy` per Mbp | unmeasured (S0) | Tier 2/3 contigs | Unbudgeted in v2.0 (BX-03). |
| Tape + pair build | Rows/s, DuckDB COPY | unmeasured (S0) | slots × partners | |
| Vehicle clustering | MinHash-LSH over flank sets, per pair-locus | unmeasured (S0) | pair-loci | O(n) with banding. v2.0's single linkage was O(n²) (DB-10). |
| Statistics + nulls | Permutation replicates × scan cost | unmeasured (S0) | pairs × replicates | Unbudgeted in v2.0. |
| Gene trees | FAMSA + VeryFastTree per family (≤5,000 leaves) | unmeasured (S0) | top-N pairs only | On demand (D-14). |

---

## 2. Raw data ingestion

### 2.1 Why nucleotide contigs, not pre-called proteins
v2.0 listed five blind spots of protein-only ingestion: missing smORFs, genetic-code mistranslation, programmed frameshifts and introns, shadow ORFs over ncRNAs, and erased intergenic elements. These still justify starting from nucleotides. The raw nucleotide sequence of every ingested contig, including deduplicated members, is retained and retrievable (§4.0).

### 2.2 Sources (D-01, D-23, DQ-05, DQ-06, DQ-15)
Each source has a contract in `sources.yaml` with these fields:
- accession scheme, release/version, retrieval date, checksum;
- record and bp counts filled from the source's own manifest;
- license, plus a `redistribution` flag (`public` | `restricted` | `byod`);
- citation;
- stage of first use.

Counts in this document are illustrative; `sources.yaml` holds the authoritative numbers.

v2.0's single `tier` field is split into four orthogonal fields (DQ-15):
- `source_class`: isolate, MAG, SAG, MGE database, SRA assembly;
- `replicon_class`: chromosome, plasmid, virus, provirus, unknown, assigned by geNomad;
- `assembly_quality`: complete, HQ draft, MQ draft, fragment;
- `molecule_type`: DNA, RNA.

`tier` (1 HQ genomes, 2 MGE/virome/plasmid, 3 unbinned/SRA) remains as a derived view for continuity.

| Tier (derived) | Sources (stage of first use; see STAGES.md) | Current release facts (verify at freeze) | License flag |
|---|---|---|---|
| 1 | Golden set (S0); small-genome GTDB reps (S3); all GTDB R226 species reps (S5); all GTDB R226 genomes plus the RefSeq complete flag (S6); AllTheBacteria, GlobDB, proGenomes4 (S7+) | GTDB R226: 143,614 species / 732,475 genomes. AllTheBacteria: about 2.44M genomes (MIT). GlobDB r226: 306,260 species reps. proGenomes4 supersedes proGenomes3. | GTDB terms to confirm (AR-13) |
| 2 | INPHARED complete phages (S1); PLSDB (S2); IMG/VR v4 HQ vOTUs and optionally IMG/PR (S4, local-only); ICE/IS/integron resources where licensed (S4+) | IMG/VR v4: 15.7M sequences / 8.7M vOTUs | IMG/VR and IMG/PR: JGI policy (restricted). ISfinder: no download without authorization (byod). ICEberg, TnCentral: terms not found. |
| 3 | SPIRE unbinned, JGI/ENA environmental, long-read compendia, SAGs, AllTheBacteria rejects (S7+); Logan (last, off-box) | Logan v1.2: about 38M accessions, hundreds of TB compressed. GORG-Tropics about 12k SAGs; the ">50k SAGs" figure is unsupported. | Logan: no license; NIH GDS policy applies. GSA-Human: controlled access, do not ingest. |

### 2.3 Ingestion pipeline (D-04, D-05, D-06; BX-02, BX-15, BX-16, BX-18, DQ-01..03, DB-04, AR-03, DS-02..04)

The v2.0 single-pass Cuckoo filter, RocksDB and 64-bit sketch sieve is removed. It had four problems:
- First sightings were never written to RocksDB, so the second-sighting comparison had nothing to compare against.
- A 64-bit 1-bit MinHash with `POPCNT ≤ 6` corresponds to about 99.3% ANI, not 95% (expected popcount at 95% ANI is about 22.3, sd 3.8).
- An exact hash over syncmers is not identity-invariant: expected Jaccard is about 0.21 at 90% AA identity.
- Results depended on stream order.

Its replacement:

```text
Stage −1  Identity resolution (per source batch)
  • Join on BioSample / run / GCA↔GCF pair → exact sequence digest (BLAKE3-128 of the
    uppercase sequence) → containment (skani; containment ≥ 0.99 and ANI ≥ 99.5%).
  • Representative precedence: RefSeq complete > GenBank complete > AllTheBacteria >
    GTDB draft > GlobDB/SPIRE MAG > Logan.
  • Every non-representative copy becomes a `member_of` row. Its sequence stays in the
    sequence store (§4.0); only its counting weight changes.
  • Output: source overlap matrix (published per release), contig_uid assignment.

Stage 0a  Run/sample preselection (SRA-derived inputs only)
  • Select runs by SRA metadata (library_strategy, organism, STAT taxonomy profile)
    before any download. Replaces the 64 MB eukaryotic-repeat Bloom filter (BX-02).
Stage 0b  Screening (per contig; rejects recorded with reason codes)
  • FCS-adaptor + UniVec; phiX / lambda / ERCC spike-ins; human/vertebrate host
    screen (minimap2 or sketch containment vs GRCh38/T2T + common hosts; low RAM);
    length/complexity (DUST).
  • FCS-GX (domain of life) deferred (decided 2026-10-29): it needs ~470 GB RAM
    (verify), above the 200 GB scope. It is accreted later as a QC sidecar when a
    ≥512 GB host is available; flagged contigs are masked from counting at the
    next release. Each genome records fcs_gx_status (ncbi_screened |
    screened_here | not_screened); MR-19 tests the mask.
  • CheckM2 + GUNC on every MAG/SAG; GTDB-Tk or skani placement for non-GTDB genomes.
  • geNomad (replicon_class, ICTV taxonomy) on every contig.
  • Chimera flags: coverage breaks where available, mixed lineage along the contig.
  • Obelisk / viroid rescue branch: circular contigs 700–1,500 bp with ribozyme or
    rod-like structure evidence are kept even with 1–2 ORFs (single size range,
    fixing the ≥750 bp vs 700–1,500 bp inconsistency, BX-19).

Stage 0c  Taxonomy (D-06)
  • Prokaryotic contigs: MMseqs2 taxonomy (2bLCA) or CAT against the pinned GTDB proteins.
  • Viral/plasmid contigs: geNomad ICTV + iPHoP / CRISPR-spacer host prediction.
  • Stored: tax_lineage, tax_rank_resolved, tax_confidence, tax_method_id. Contigs below
    the confidence cutoff never count toward order dispersion or null strata.

Stage 0d  Deterministic architecture dedup (two-pass sort-merge, counting only)
  • Pass 1 (stateless): emit (arch_key, biome_bucket, quality_rank, contig_uid,
    ≥1,024-bit FracMinHash/containment sketch) into hash partitions on local disk.
  • Pass 2: sort each partition; within a key, greedy clustering by containment with
    deterministic representatives (quality_rank, then contig_uid); verify candidate
    pairs with skani.
  • Output: contig_member.parquet (member_uid → rep_uid, containment, ani_est).
    Dedup can be redone at a new threshold without re-reading sources. Members keep
    their sequences and coordinates.
  • Counting weight: at most one count per BioProject per vehicle cluster for inference.
    The raw count is kept as a descriptive column (DS-14/DS-09).
```

Logan is processed last, compute-to-data on AWS (us-east-1). Accessions are preselected by Stage 0a; only screened, called and deduplicated outputs leave AWS (D-05, AR-01). The streaming rate, collapse ratio and output volume are unmeasured. They are estimated from a 0.1% accession sample before any commitment.

---

## 3. Annotation stack (`bio-box`)

Every annotation records its tool, version, database release, parameters and evidence level in `slot_provenance` (§4.2) (DQ-07, DQ-11).

### 3.1 Stage 1: genes, smORFs, ncRNA, arbitration and split genes

1. **Gene calling (D-07, BX-06).**
   - **Complete and HQ genomes:** pyrodigal in single mode per genome. The genetic code is chosen by Codetta, or by genome-level comparison of coding density over codes 11, 4, 15 and 25.
   - **Viral contigs:** pyrodigal-gv, with geNomad's code call.
   - **Fragments:** meta mode. The code is inherited from the contig's taxonomy, because meta mode cannot detect codes 15 or 25.
   - **S0 implementation, gene-call policy `p2` (approved 2026-10-29, `planning/DECISIONS.md` DEC-13):**
     - Viral contigs: pyrodigal-gv in meta mode, with the code taken per gene (11, 4 or 15). Until geNomad regimes exist, a source's contract declares its kind (INPHARED, IMG/VR: virus).
     - Standalone plasmids: meta mode; retrained with the host chromosome once hosts arrive at S2C.
     - Prokaryotic genomes ≥100 kb: single mode trained on all contigs together with the genome's code; meta below 100 kb (Prodigal's advice; Prokka uses the same switch).
     - Circular contigs: genes are called on the contig with 50 kb of the other end added on both sides and kept when their start lies on the contig, so an origin-spanning gene is called once with `end > contig_len` (Bakta merges edge genes the same way).
     - 20–100 kb genomes are run in both modes on a sample at S1 and the gene-count and length differences recorded.
     - On the full golden set (`golden.yaml`, 15 genomes, 43 Mbp; ring R3, 2026-10-29): ingest 0.4 s and gene calling 13 s on a laptop, 33,716 genes, 3 of them wrapping an origin; all complete genes match an independent translation under codes 4, 11, 15 and 25. pyrodigal-gv calls code 15 for Lak-A1 and also for **Lak-B8, whose NCBI record says code 11**: under code 11, Lak-B8 has coding density 0.69 and 379 genes ending at TAG, and under code 15 it has 0.895 and none, the same signature as Lak-A1. Source genetic-code metadata is therefore treated as a hypothesis checked against the sequence, never as truth.
     - On the micro-S0 fixture: λ 62 genes and T4 265 (pyrodigal-gv; identical to plain meta mode for these code-11 phages), *M. genitalium* 527 under code 4 (995 much shorter genes under code 11). Every complete gene matches an independent translation (R-04); ≥98% of genes are recovered at mirrored coordinates on the reverse complement (MR-03); rotating *M. genitalium* so the origin cuts a gene recovers ≥99% of genes and calls the cut gene once, whole (MR-04).
   - **Host labels, corrected from v2.0:** code 4 applies to Mycoplasmatota and Spiroplasma. Code 15 applies to some phages (e.g. Crassvirales); Absconditabacteriales use their own code. Code 25 applies to SR1 / BD1-5 lineages, not all CPR, since most Patescibacteria use code 11. The v2.0 claim of code 4 in Hodarchaeales is unverified and removed.
   - **Recorded fields:** `partial_5prime`, `partial_3prime`, `orfs_to_contig_start`, `orfs_to_contig_end`.
2. **smORFs, 15–75 aa (D-07, BX-07).** This range replaces the four inconsistent ranges in v2.0.
   - The five initiation modes from v2.0 generate candidates:
     - pyrodigal initiation model;
     - Shine-Dalgarno or translational coupling;
     - leaderless TATA/BRE or -10;
     - Bacteroidota AT-rich window;
     - profile, TM or conservation rescue.
   - A candidate is promoted to a slot only if its 50%-AA cluster appears in at least 3 independent orders with RNAcode support or dN/dS < 0.5 (Sberro 2019 approach). The promotion step is budgeted (§1.3).
   - Unpromoted candidates go into a sidecar with `evidence_level = candidate`. They are excluded from the default statistics.
   - Validation uses Ribo-seq and Ribo-RET sets (§5 B-G).
3. **ncRNA (D-08, BX-04, BX-05).**
   - **Step 3a.** Extract intergenic regions ≥40 bp plus 25 bp CDS flanks, and cluster them at 90% nt identity (IGR90). The clusters only choose what to scan.
   - **Step 3b.** The CM library is the Rfam prokaryotic and viral subset, derived from Rfam's taxonomy tables. The plan does not assume a count; Rfam 15.1 lists 4,227 families. Add the orphan CMs (Weinberg, Stav, Brewer, Breaker). Each orphan CM is run through `cmcalibrate`, given its own reporting threshold (it has no GA threshold), and deduplicated against Rfam.
   - **Step 3c.** Run `cmsearch` with Infernal's own filter pipeline. There is no extra nhmmer prefilter, which would duplicate Infernal's filters and lose families conserved mainly by structure. `--cut_ga` applies only to Rfam models.
   - **Step 3d.** Assign the pocket subtype per member sequence, not per IGR90 representative. Each hit's cluster members are re-aligned with `cmalign`, then the `cm_pocket_rules.json` columns are read. This matters because a 90% cluster tolerates about 10 nt differences, while ykkC and purine subtypes differ by 1–5 nt.
4. **CRISPR, integrons, tRNAs and terminators.**
   - Minced plus CRISPRDirection, with spacers extracted to `spacer_protospacer_graph.parquet`.
   - IntegronFinder 2.
   - tRNAscan-SE 2.0 and Aragorn.
   - TransTermHP calls are stored as a **score** (`terminator_score`), not as a hard directon break (BX-17). The `3'-operon-stop` and `5'-cis-attenuator` classes are kept, and the attenuator exemption from v2.0 stays. Directon thresholds are validated per clade against RegulonDB and term-seq.
5. **Stage 1.4 overlap arbitration.** The rules from v2.0 stand:
   - purge shadow ORFs that lack an HMM hit when they overlap CRISPR arrays, rRNAs, tRNAs or tmRNAs, or the aptamer cores of confident Infernal hits;
   - keep `uORF_cis_leader`, CreT, and antisense Type I TA pairs.

   Purged ORFs are recorded in the rejection ledger with their coordinates; they are not silently dropped.
6. **Stage 1.5 split-gene reconstruction (BX-08).** Adjacent same-strand ORFs that hit consecutive halves of one profile are merged only with independent evidence:
   - a slippery-site motif (`prfB`, `gpG-gpGT`, IS `OrfAB`);
   - an in-frame TGA with SECIS or selenocysteine machinery, or a TAG with the pyl genes;
   - a Group I/II intron call or HEG.

   Otherwise both pieces are kept and flagged `CANDIDATE_INDEL`. Assembly indels in ONT, short-read or Minia contigs produce the same signature as real frameshifts, and real pseudogenes must not be erased. Merged ORFs keep both exon coordinate sets in `slot_provenance`.

Every slot carries `evidence_level ∈ {curated_model, profile_hit, promoted_smorf, rule_merge, candidate}` (DQ-14).

### 3.2 Stage 2: domains, subfamilies and structure

```text
All called proteins (count: estimate, pilot-pending; Tier 1 alone ≈ 7–10×10⁹ before dedup)
  │
  ├─► Stage 2a: MMseqs2 NR90 → NR50 (incremental across stages with `clusterupdate`, D-19)
  │     • Terminal-window guard for proteins ≥350 aa: N- and C-terminal 120-aa windows must
  │       also meet the identity threshold, so polymorphic-toxin CTD variants stay separate (BIO-01).
  │     • Every member keeps its own sequence and its alignment to the representative in
  │       `cluster_member_aln` (L2; one row per member per level): `member_uid`, `rep_uid`,
  │       `level` (nr90 | nr50 | subfam), qstart, qend, tstart, tend, pident, alnlen and a
  │       compact CIGAR (~40 B/row; ~19 GB at S5, estimate). Downstream residue-level calls
  │       and library hits are projected per member through it (SD-05).
  │     • 3Di transfer (D-10, BX-11): only from Foldseek-distributed 3Di (AFDB, ESMAtlas30
  │       about 36M); only at ≥90% identity to models with pLDDT ≥70 over ≥70% of the envelope.
  │       Per-residue pLDDT is carried; `fold_source = transferred_3Di`.
  │       "ESMAtlas-Full 772M" 3Di is not distributed and is not assumed.
  │
  ├─► Stage 2b: domain dissection (D-09, BX-09, BX-10)
  │     • Libraries: Pfam, NCBIfam, PHROGs, arCOG, each with its own GA/TC cutoffs.
  │       COG via eggNOG. CDD only via RPS-BLAST if used. HMMER Z fixed and recorded:
  │       `-Z 1e8 --domZ 1e8` sequences, one constant for the release line (SD-04). Each
  │       hit stores bitscore, E at that Z and the real target size. Scan unit per SD-05
  │       (exact-unique proteins through S2C; NR90 reps + projection from S3).
  │       Every hit keeps source library, model accession, version and cutoff pass in
  │       `slot_annotations`; there is no merged "115k non-redundant HMM" set.
  │     • DefenseFinder, CONJScan and TXSScan run as MacSyFinder system calls on
  │       finished tapes. PADLOC uses its own engine (HMMER searches against PADLOC-DB,
  │       then its own presence/absence and synteny rules per system model; Payne et al.
  │       2021, doi:10.1093/nar/gkab883), not MacSyFinder. CRISPRCasTyper and
  │       CRISPRCasFinder add CRISPR-Cas typing. All are stored as sidecars, not as tape
  │       tokens. This keeps the B1 temporal holdout clean (§5).
  │     • Polymorphic-toxin CTD slicing (BIO-01), nested/discontinuous stitching (BIO-02),
  │       circular permutation (BIO-03) and repeat collapse (BIO-04) follow v2.0 rules,
  │       but start from Pfam nesting annotations and ECOD domain boundaries. The gap
  │       rules are labeled heuristics. One repeat-gap threshold: ≤35 aa.
  │     • subfamily_50_id: 50% ID / 80% coverage clusters of excised envelopes within each
  │       seq_family_id. Stable key = `cluster_uid` as defined in §4.1 (D-19); the
  │       representative envelope's UID is stored alongside it.
  │     • Catalytic motif (D-12, BX-12, DB-13, DQ-07): residues at curated site columns
  │       stored as `catalytic_motif_packed: uint64` (≤12 residues × 5 bits + 4-bit length)
  │       with `site_source ∈ {M-CSA, CDD_site, conservation_only}`. Computed per member,
  │       not per representative. The hash is a join key only. The Hamming penalty
  │       applies only to curated sites.
  │
  ├─► Stage 2c: predicted structure (GPU stage; not run on wilkcomp02)
  │     • Foldseek's ProstT5 path (`createdb --prostt5-model`) for recurrent dark NR50
  │       segments without a transfer. `fold_source = predicted_3Di`, with head entropy.
  │       Foldseek DBs built from ProstT5 predictions do not support TM-score alignment;
  │       v2.0's `--alignment-type 1` on them is removed (BX-11).
  │     • ECOD assignment by Foldseek with E-value stored. AlphaFill/BioLiP2 pockets only
  │       for hits to AFDB models.
  │     • Measured ECOD_T accuracy on post-cutoff PDB depositions is published with each
  │       release (DS-15).
  │
  └─► Stage 2d: remote homology for multi-member dark clusters
        • `mmseqs result2profile` (no loose files) + hhalign against candidates nominated by
          Foldseek or profile search. The rate is unmeasured.
```

"Recurrent" (≥2 non-identical contigs) needs occurrence counts that exist only after tape construction. Stage 2c therefore runs after a global barrier following Stage 3 tape counts (AR-16).

### 3.3 Biophysical badge (D-18, DB-14)
`biophysical_badge: uint64`, with `badge_version` stored in file metadata and accessor macros in the catalog:

| Bits | Field | Notes |
|---|---|---|
| 0–4 | `tm_count` | 0–31, from TMbed / pyTMHMM |
| 5–8 | `signal_type` | 0 none, 1 Sec/SPI, 2 Tat/SPI, 3 lipoprotein SPII, 4 prepilin SPIII, 5 secretion-system substrate |
| 9–16 | `sorting_tag` (8 flags) | PEP-CTERM, LPXTG, GlyGly-CTERM, PGF-CTERM, archaeosortase PGF, BMC EP, encapsulin CLP, RiPP leader |
| 17–21 | `pi_quantized` | pI 3.0–12.0 in 0.3 steps |
| 22–29 | `motif_flags` | CxxC..CxxC, CxxxCxxC, HEXXH, GxGxxG, DxD, H..H..E, coiled coil ≥28 aa, disorder ≥40% |
| 30–41 | `length_aa` | saturating at 4,095, with an overflow flag in bits 50–53 |
| 42–49 | `repeat_unit_count` | 1 for non-repeat segments |
| 50–53 | edit/overflow flags | frameshift or Sec/Pyl merge, intron splice, CANDIDATE_INDEL, length overflow |
| 54–63 | reserved | |

A round-trip property test covers every field at its maximum (§5 B-0).

---

## 4. Data layer and query architecture

The minimum architecture is **Parquet + DuckDB + a thin API** (D-17, AR-06). Custom index structures are built only when a measured SLO from §1.2 fails at a stage gate, and their designs are kept in §4.6. Layering, file geometry and tests are specified in [`planning/data_layer_and_tests.md`](planning/data_layer_and_tests.md):

| Layer | Contents | Update rule |
|---|---|---|
| L0 | Source snapshots and the sequence store | Append-only |
| L1 | Per-shard annotations | Append-only |
| L2 | Global clusters | Incremental, with crosswalks |
| L3 | Tapes and pair postings | Append by partition; global sort at release |
| L4 | Statistics | Recomputed per `stats_release_id` |
| L5 | Optional serving indexes | Built only if needed |

### 4.0 Sequence store (L0)
Raw nucleotide sequences stay available to any project. That covers every source contig, including deduplication members and contigs rejected at screening (rejects go in a quarantine partition with their reason code), plus called proteins and IGRs. No stage discards a non-representative sequence. Deduplication and clustering change only counting weights and representative pointers.

- **Identity.** `contig_uid` = BLAKE3-128 of the uppercase nucleotide sequence (the single UID hash scheme; `raw_digest` = sha256 of the exact source bytes; L0 objects = sha256). Source accession, version and coordinates are recorded in `contig_occurrence`. `protein_uid` is BLAKE3-128 of the amino-acid sequence. IGRs are stored as coordinates `(contig_uid, start, end, strand)` and sliced from the contig on demand, with `igr_uid` = BLAKE3-128 of the sliced sequence. Identical sequences from different sources share a uid and keep one occurrence row per source copy.
- **Storage.** Primary: Parquet `seq_chunks` (64 KB chunks, byte-exact, sorted by digest) plus `protein_seq`, as specified in `planning/data_layer_and_tests.md` §2A. bgzip FASTA with `.fai`/`.gzi` is generated on demand by `bio-box seq extract` for tools. The S0 bake-off (test K-13) is decided by LY-06 (§0.4 item 9): the format that reads fewer row groups (BGZF blocks for bgzip) for the canonical slice queries on the micro-corpus wins, and R-08 latency is checked afterwards on the winner. zstd-seekable FASTA is not used (no samtools/pysam support).
- **Retrieval.** `bio-box seq get (contig_uid, start, end, strand)` returns the slice. Protein and IGR retrieval works by uid, or by slot through its coordinates.
- **Lineage.** slot → (`contig_uid`, start, end, strand) → `contig_provenance` (source, accession, version, source coordinates, retrieval checksum). For dedup members, `contig_member` maps member to representative without removing the member.
- **Integrity.** Re-reading each shard and recomputing each digest must reproduce the source checksum and the stored uid (§5 B-L0).

### 4.1 Identifiers, tokens and provenance (D-18, D-19, AR-04, DB-11, DB-16, DQ-07, DQ-13)

**Stable keys.**
- `contig_uid`, `protein_uid` and `igr_uid` are content digests.
- Cluster keys (`nr90_uid`, `nr50_uid`, `subfamily_uid`) are `cluster_uid = blake3(sorted member UIDs at creation, cluster_params_hash)`, 128 bits, stable until a split or merge recorded in `id_history`. The representative UID is stored with each cluster, not just a cluster index (as in `planning/data_layer_and_tests.md` §2).
- `slot_uid` is the composite `(contig_uid, slot_ordinal)`, where `slot_ordinal` is the slot's index in the contig's tape lists, stable per `gene_call_version`. No separate slot digest is stored.
- All UIDs are 16 bytes (`FIXED_LEN_BYTE_ARRAY(16)` in Parquet) wherever they appear in L0–L4, including list columns.
- Stable keys are the only legal keys for sidecars and cross-release references.

**Release-local dense IDs.** `dense_contig_id` and `dense_slot_id` are serving IDs (L5 scope). They are used only for serving indexes and for the physical sort order the release writes (the `dense_contig_id` scalar in `contig_tapes` is that sort column, regenerated at every release). They are never a join key in L0–L4, never stored in list columns and never used by sidecars. They are never renumbered within a release, and compaction writes a new release.

**Release fields.**
- `release_id` is on every table.
- `stats_release_id` is on every L4 table.
- `id_history.parquet` maps clusters across releases as same, split, merged or retired.

**Token namespace.** `token_id = level(3 bits) | source(5 bits) | local(24 bits)`, widened to 64-bit if counts require it. The tape column holds the canonical pick under a documented precedence. `slot_annotations(slot_uid, token_id, evalue, bitscore, env_start, env_end, cutoff_pass, library_version)` holds all hits.

**Enums.** All enums (`feature_kind`, `junction_topology`, `regime`, `evidence_level`) map to Sequence Ontology terms in `enum_registry.parquet`. They include rRNA, intron, homing endonuclease, pseudogene and ribozyme (DQ-16).

**Provenance tables.**
- `contig_provenance`
- `slot_provenance`: caller, smORF mode, genetic code, edit rule, evidence level, tool version
- `cluster_member_aln` (member → representative per level, with alignment coordinates; §3.2)
- `contig_member`
- `tool_versions`
- `annotation_provenance`: tool, version, DB release per slot group

### 4.2 Tables

#### `contig_tapes` (one row per contig; sorted by `dense_contig_id`, stored explicitly)
Per-contig scalars:
- `contig_uid`, `dense_contig_id`
- `source_class`, `replicon_class`, `assembly_quality`, `molecule_type`, derived `tier`
- `genome_or_sample_uid`, `bioproject`, `biome_envo` (ENVO / MGnify biome ontology, DQ-18)
- `tax_lineage`, `tax_rank_resolved`, `tax_confidence`, `tax_method_id`
- `is_circular`, `genetic_code`, `n_slots`, `release_id`

Parallel `LIST<T>` columns, one entry per slot:

| Column | Type | Notes |
|---|---|---|
| (no `slot_uids` column) | n/a | The stable slot key is the composite `(contig_uid, slot_ordinal)`; `slot_ordinal` is the list index (§4.1). The earlier `LIST<uint64>` digest of contig_uid + coordinates + kind is withdrawn: it was narrower than the 16-byte UIDs and duplicated the composite key. |
| `orf_idxs` | `LIST<uint16>` | Stitched across validated merges only |
| `seg_idxs_in_orf` | `LIST<uint8>` | |
| `parent_seg_idxs_in_orf` | `LIST<int16>` | −1 for top level; same domain as `seg` (DB-15) |
| `domain_topology_flags` | `LIST<uint8>` | discontinuous host, CP, repeat array, CTD tip, edit type |
| `feature_kinds` | `LIST<uint8>` | SO-mapped enum |
| `strands` | `LIST<int8>` | |
| `start_bps`, `end_bps` | `LIST<uint32>` | Contig coordinates; envelope residues in `env_start_aa`, `env_end_aa` |
| `delta_bp_prev_slots` | `LIST<int32>` | int16 overflowed after large CRISPR arrays and introns (DB-15) |
| `junction_topology_prevs` | `LIST<uint8>` | 0 intra-ORF fusion, 1 co-directional, 2 divergon, 3 convergon/excludon, 4 sense-antisense nested, 5 inter-locus break |
| `directon_ids` | `LIST<uint16>` | Attenuator exemption kept |
| `regulon_locus_ids` | `LIST<uint32>` | Type unified with critique_1 (DQ-16) |
| `orfs_to_contig_starts`, `orfs_to_contig_ends` | `LIST<uint16>` | Censoring denominators |
| `local_regimes` | `LIST<uint8>` | From geNomad / CheckV / IntegronFinder 2 / CONJScan / ICEscreen / ISEScan, including `UNKNOWN_SHORT_CONTIG` (D-11, BX-13) |
| `regime_scores` | `LIST<float16>` | |
| `vehicle_cluster_ids` | `LIST<uint32>` | MinHash-LSH over ±15-ORF flank sets (DB-10) |
| `clan_or_cog_ids`, `seq_family_ids`, `subfamily_ids` | `LIST<uint32>` | Namespaced token IDs |
| `catalytic_motif_packed` | `LIST<uint64>` | Plus `site_sources: LIST<uint8>` |
| `catalytic_intacts` | `LIST<uint8>` | Three-valued: intact, broken, unknown |
| `struct_fold_ids` | `LIST<uint32>` | Null where no fold evidence |
| `fold_sources` | `LIST<uint8>` | transferred_3Di, predicted_3Di, hmm_inherited, none |
| `fold_confidences` | `LIST<float16>` | |
| `evidence_levels` | `LIST<uint8>` | §3.1 |
| `biophysical_badges` | `LIST<uint64>` | §3.3 |
| `protein_uids` | `LIST<FIXED_LEN_BYTE_ARRAY(16)>` | BLAKE3-128 `protein_uid` into the sequence store; null for non-protein slots. The reverse map `protein_uid → slots` is `protein_occurrence` (§4.7). |

Q2 neighbourhood fetches read rows by the explicit sorted `dense_contig_id` column, with row-group statistics pruning (sort order, statistics and row groups read are checked by LY-02, LY-03 and LY-06). v2.0's `rowid IN (...)` claim is withdrawn (DB-05).

#### `slot_pairs` (replaces `window_postings`)
Directed pairs within ±k ORFs, or fused, inside one `regulon_locus_id`. Files are partitioned by anchor-token bucket so a query opens one small footer, and an `(anchor, partner) → (file, first_row_group, n_rows)` offset table is kept (DB-06).
- **Columns:** anchor and partner tokens at family and subfamily level, `orf_offset`, `junction_topology`, flags, `delta_bp` (int32), `tax_order` (with confidence), `vehicle_cluster_id`, `bioproject`, `partner_badge`, `contig_uid` and the `slot_ordinal` of anchor and partner (together the two composite `slot_uid`s).
- **Size** (derived): rows ≈ slots × mean partners per slot (about 6 after edge truncation) × token levels; bytes ≈ rows × 15–25 B. Both factors are measured from S1 onward.
- **Hub routing:** anchors with document frequency above a measured threshold are answered from `pair_stats` for Q3. `slot_pairs` serves drill-down with sampling.

### 4.3 Pair statistics (`pair_stats`) and closure (D-13, DS-01, DS-06..09, DS-12, DS-14, DS-16)

**Unit of discovery.** The primary rows are pairs at the pre-registered primary token level, `subfamily`, with family as a secondary level. A module step thresholds the pair graph into dense modules. **Modules are the reported discovery unit**, with pairs as evidence, and FDR and confidence intervals are computed at the module level by cluster bootstrap. One 5-gene operon yields 10 correlated pairs and must not count as 10 discoveries (DS-09).

**Model.** Let y_m ∈ {0,1} indicate that B lies within ±k ORFs of A at independent locus m. Loci are vehicle-deduplicated, uncensored, capped at one per BioProject per vehicle, and restricted to `tax_confidence` above the cutoff.
- Null: y_m ~ Bernoulli(p_m⁰), with p_m⁰ from the regime- and clade-conditioned background.
- Alternative: logit π_m = logit p_m⁰ + β.

**Columns.**

| Column | Definition |
|---|---|
| `z_score` | Z = Σ(y_m − p_m⁰) / √Σ p_m⁰(1 − p_m⁰): the Poisson-binomial score statistic, streaming-computable. |
| `beta_hat`, `beta_se` | Conditional-logistic MLE of the log odds ratio of adjacency, with cluster-robust SE. This is the effect size, kept separate from evidence. |
| `sum_y`, `sum_p0`, `sum_p0q0` per regime × tier; Δbp 32-bin histogram; topology counts; strand/order counts; phylum and order bitmaps; genus HLL; vehicle count; BioProject count; top-K neighbour histogram | **Sufficient statistics.** Every score is recomputable from these, so `bio-box rescore --params params.yaml` works without the base build and ships in `lite`. |
| `n_veh_indep_events`, `n_distinct_bioprojects`, `n_distinct_biomes` | Dependence-aware counts. |
| `p_b_given_a_uncens`, `p_a_given_b_uncens`, `cond_asymmetry_log2` | As in v2.0. |
| `ctx_entropy_a_bits`, `ctx_entropy_b_bits`, `slot_fold_entropy_ab_bits`, `micro_geometry_score`, `attc_between_frac`, topology fractions | **Features, kept as separate columns.** They are never subtracted from a log Bayes factor. They are combined only by the fitted model below. |
| `p_value`, `local_fdr`, `q_hier` | Two-tier null (DS-07). The far tail uses the parametric score-test distribution, checked against permutation over 10⁻³ to 10⁻⁶ (QQ plot, with the extrapolation error reported). The permutation null keeps block length ≥25 ORFs (at least 4× the window), and the pair must be split across blocks. A curveball/fastball degree-preserving content null serves as a second null, and results are reported under both. An Efron empirical null (central-90% fit, π₀) gives `local_fdr`. Yekutieli hierarchical FDR runs family → subfamily → motif and gives `q_hier`. v2.0's "target-decoy" label is withdrawn: this is a permutation null. |
| `powered`, `min_detectable_count` | Computed from Σp_m⁰ at the operative α. Under the chromosomal regime about 6–9 vehicle-independent co-occurrences are needed, for N_eff of 20–100 (derived, DS-08). Queries distinguish "tested, no evidence" from "untestable". |
| `geometry_z` (Tier 2 / small replicons) | On a small replicon a ±k window spans most of the replicon, so p₀ → 1 and Z → 0. For these, the test is of geometry (adjacency, order, Δbp, strand) against a replicon-matched background pool that carries both genes. They are ranked separately (DS-08). |
| `rank_sd` | Standard deviation of rank over K = 10 annotation realizations, sampled from top-k token and fold calls per slot (DS-12). |
| `coupling_score` | Output of the fitted, monotone-constrained model (§4.3.1). `func_coupling_llr` survives only as a derived diagnostic view. |
| Pair closure | `pdb_cocrystal_status`, `max_pdb_pair_seq_id`, `paperblast_pair_co_pmids`, `feba_pair_cofitness_max`, `feba_pair_shared_phenotype_count`, as in v2.0. FEBA columns are excluded from any feature set used to fit `coupling_score`. |
| `gene_tree_homoplasy_b_on_a` | Null in the base build. Computed on demand for the top-N pairs with FAMSA + VeryFastTree, capped at 5,000 leaves (D-14). |
| `stats_release_id`, `stale_since_epoch` | Staleness contract (D-20). |

**Mergeability.** Counts, histograms, bitmaps and HLLs merge across epochs. Median, entropy, z, q and the fitted score do not; they are recomputed at a stats release (DB-10).

#### 4.3.1 Parameters and fitting protocol (DS-02, DS-14)
- Every constant goes into `params.yaml`, which is hashed into the release manifest. Each constant is tagged `asserted`, `fitted` or `derived`, with a plausible range. The constants include: Jaccard and LSH bands; γ; envelope and ATOM coverage of 0.70; 50%/80% clustering; ±k; entropy thresholds; the IQR scale; containment thresholds; CTD window; repeat gap; leakage cutoffs.
- `coupling_score` is a monotone-constrained GBM, or an L2 logistic regression, over the feature columns. It is trained on half of a gold set and evaluated on the other half, with AUPRC and cluster-bootstrap CIs.
  - Positives: RegulonDB/ODB verified operons, PDB assemblies with real interfaces, MetaCyc/GapMind adjacent steps, MIBiG same-BGC.
  - Negatives: the pre-registered trap panel.
- γ (closure penalty) is fitted by maximum likelihood on M-CSA plus Swiss-Prot `ECO:0000269` mutagenesis labels. The penalty is information-weighted: exp(−γ Σ_j w_j 1[mismatch_j]), with w_j proportional to the column's information content. v2.0's γ = 1.4 was set to clear its own 0.25 pass line and is removed.
- **Sensitivity gate.** Run Morris screening, then Sobol indices, on a stratified sample. Each constant at ±50% must keep Kendall τ ≥ 0.9 on the top-100 modules. Otherwise the constant is published as fragile, and results are reported under both settings.

#### 4.3.2 `subfamily_closure_card` (CLOSE-01, updated)
Keyed by (`subfamily_uid`, `catalytic_motif_packed`), with a rollup view by `seq_family_id`.
- **Gates.** The envelope gate (EnvCov ≥ 0.70) and the resolved-ATOM gate (≥ 0.70) from v2.0 stand. `pdb_is_structural_genomics_only` and `max_swissprot_env_seq_id` with `ECO:0000269` are kept.
- **Isofunctional closure score.** C_iso = max_T [S_seq(seq_id_env) · exp(−γ Σ_j w_j 1[mismatch_j]) · w_evidence(T)], with γ fitted as above. The motif penalty applies only when `site_source` is M-CSA or a CDD site; otherwise the score carries `motif_unverified = true`.
- **Provenance.** Each row keeps its evidence pointers: PDB IDs, UniProt accessions, PMIDs.

### 4.4 Pangenome planes and transitions (D-15, DS-13, DB-09, BX-16)
1. **Three-valued planes per token.** Each species cell is present, absent or unknown. `core_p` and `pan_p` are Beta-binomial posterior means, with a genus-level prior and `n_genomes_in_species` stored as the denominator. Species below the power threshold are `unknown`, not 0. The `mge_accessory` plane is derived from `replicon_class` and the strain fraction.
2. **Storage.** Compressed bitmaps over species leaf IDs, ordered by tree DFS. Dense padded rows are used only for the top ~10⁵ tokens by species count. Row width is derived per release from the species count (17,984 B at R226). A dense 3-plane row costs 3 × row width per token, so dense rows for millions of subfamilies would need tens to hundreds of GB (DB-09).
3. **Quality mask.** One definition: a species counts as high quality if it has a closed chromosome, or ≥3 independent genomes that pass CheckM2/GUNC thresholds. The v2.0 documents disagreed between ≥3 and ≥5. Contig-edge masks are as in v2.0.
4. **NOGD scoring.** Hole sets are computed on the three-valued planes, and unknown cells never count as holes. Candidate replacements are tested with a Pagel-style correlated-evolution likelihood-ratio test over the precomputed gain/loss transitions, which reports a p-value. v2.0's `hole & (hole >> 1)` "coherence check" had no null and is removed.
5. **`transition_edges`.** Fitch/Sankoff transitions are precomputed once per tree release (at <400 B per token) and serve as sufficient statistics for the test.
6. **Order count check.** The number of GTDB orders is checked against any fixed-width order bitmask at build time. If it does not fit, Roaring is used instead.

### 4.5 `gLM-Domain`: separate later project (D-16)
Model training and the masked-slot index are not on the critical path, and no stage depends on them. The requirements recorded for when the project resumes are:
- **Splitting.** GraphPart-style balanced homology partitioning, not single-linkage components. It excludes hub and housekeeping families from the edge definition, and a gate rejects any split whose largest component exceeds 20% of operons. Achieved leakage is reported as a measured bound, not as "0.0%".
- **Temporal split.** An additional temporal split, plus a phylum-plus-homology holdout.
- **Vocabulary.** Chosen on the training split only. The 131,072-token cap with compositional representation of rare tokens is carried over from v2.0.
- **Architecture.** 2D split RoPE, a pairwise geometry bias, and a coarse-to-fine decoder, all carried over from v2.0. The RoPE dimension split is reconciled; v2.0 gave both 256/256 and 48/16.
- **Controls.** A frequency-prior baseline, a context-shuffled ablation, ECE < 0.05 after temperature scaling, risk–coverage curves, metrics stratified by frequency decile, and comparison against gLM2.

### 4.6 Deferred indexes: build only if an SLO fails at a stage gate
These designs are kept so that an SLO failure has a ready answer. None is built by default.

1. **Index 1: proximity bitmaps (corrected design, DB-02, DB-03, DB-08).**
   - **Bitmaps.** `token_to_slots` Roaring bitmaps over local slot IDs, partitioned per (`replicon_class`/tier, epoch) so that each ID space stays well under 2³². Proximity is `A ∧ ⋁_{δ=−k..k}(B ≪ δ)` via `roaring_bitmap_add_offset`. The predicate is defined on ORF index, using an ORF-level bitmap, so modular NRPS/PKS genes do not break it. A `contig_start_slots` bitmap gives the boundary check.
   - **Slot attributes.** A fixed-stride `slot_attrs.u32[slot]` array holds strand, kind, junction, local directon and edge flag.
   - **Removed from v2.0.** Aligned 16-slot windows lost about 22% of true pairs at d = 1..6 (derived: (1+…+6)/(6·16)). `orf_mod_8` aliased d with d ± 8.
   - **Size.** About 2 B per posting × postings × levels. This is an estimate, of about 35 GB at the v2.0 slot count (DB-08); it gets re-derived from measured slots.
   - **Hubs.** Hub bigram bitmaps are defined by a measured document-frequency threshold. v2.0's "2,048 tokens × >2×10⁶ occurrences" needs ≥4.1×10⁹ postings, which exceeds the total slot count (DB-07).
2. **Offset-indexed tape store (DB-05).** `tape_offsets.u64[dense_contig_id]` points into a per-contig, block-compressed, row-oriented heap. It is built only if Q2 misses its target on DuckDB. Lance is the alternative to benchmark.
3. **Vector index.** Belongs to the gLM project (§4.5).

### 4.7 Ad hoc search (SD-01..SD-21)
Authoritative design: [`planning/search/SEARCH_PLAN.md`](planning/search/SEARCH_PLAN.md). This section is a summary; where they differ, SEARCH_PLAN wins. Nothing here is implemented.

**Evidence the design rests on** (SEARCH_PLAN §1). 348/365 synteny hypotheses need a neighbourhood step; fixed-library annotation answers about a third; searching outside gene calls is routine; no 2023–2026 corpus-scale paper runs an iterative profile search over the whole corpus (they reduce first, then search 10^3–10^6 objects); at 1e9 proteins most tools are bound by memory or I/O.

**Recipes and runner (SD-01, invariant 15).** A search is a `search_run` launched from a hashed YAML recipe: seed (accession, HMM, MSA, CM, structure, locus or region, with source), target, tool and version, iterations and stop rules, Z, curation checks, context window, required positive control, outputs. One runner (`bio-box search run`) launches every tool; non-zero exit, timeout and OOM are `failed`, never zero hits. `serve` submits and polls; jobs are async (§1.2).

**Locus context on every hit (SD-02).** `protein_occurrence` (L1, per stage; ~29 GB at S5, estimate) maps each protein to every slot, contig, coordinates, strand, genome, taxonomy, species-rep flag and tier. Hits in `search_hit` join to it.

**Reduce first (SD-03, SD-05).** Iterative rounds run on NR50 representatives; the final round runs on NR90 representatives, projected to members through `cluster_member_aln`, with near-cutoff members re-scored exactly. Through S2C every run also searches all exact-unique proteins as its recall check; from S3 on, exhaustive search is a sampled recall check. Library scans follow the same rule (§1.3).

**Z policy (SD-04).** One search space per release line: `Z_prot = domZ = 1e8` sequences; `Z_nt = 1e5` Mb for Infernal and nhmmer. Every hit stores bitscore, E at Z and real target size. Recipe thresholds are stored as bits; paper E-thresholds are converted at the paper's stated database size. Test S-03 checks that bitscores and E at Z are identical across stages. The constants are an open item for the user (SEARCH_PLAN §9).

**Search pack (SD-06, SD-07, invariant 16).** Built once per release as ≥256 MB shards on NFS with a manifest; runs copy shards to `/tmp`. Contents: NR50 and NR90 FASTA plus an MMseqs2 index, `unique.dmnd`, `segments.faa` (domain envelopes plus uncovered segments ≥50 aa), `igr.fna`/`igr90.fna`, `dark_orf.faa`, `spacers.fna`, pressed models, HH-suite Pfam-HHM and PDB70, and from ST3 a 3Di DB of NR50 reps (`3di_nr50/`) built with the encoder that wins the SD-21 bake-off (ProstT5 vs ESM-2 3B; reseek tested as an alternative aligner). S5 total about 370–530 GB (estimate; row sum in `planning/sources_and_sizing.md`).

**Iteration controller (SD-10).** Per round: search at fixed Z; collapse to NR50 reps and record first-inclusion round; check coverage, required motif, foreign-domain alarm, locus-context retention and reverse search; build the next MSA from hit envelopes only; log hits and profile; stop on the recipe's rule, reverting one round on a drift alarm. Significance is reported at the round of first inclusion.

**Neighbourhood API (SD-08).** `bio-box nbr` takes anchors (run, token or UIDs), a window in genes or bp, strand and operon-gap predicates, a contig-edge policy (exclude, truncate, require N each side) and circular topology, and reports distance in genes and bp with censoring flags. `nbr score` runs positive-set enrichment.

**User families into the statistics (SD-09, SD-11, SD-12).** `bio-box search promote` turns a run's final hit set into a user token (`source_db = user`) as a sparse `user_token_overlay` on the tapes. `bio-box search stats` runs the §4.3 score test against per-release `background_rates`. Registering the model in `user_models.yaml` makes it a library model at the next release. `domain_promiscuity` feeds a seeding blocklist and a fusion down-weight. Counts are reported raw, vehicle-independent, per species rep and per distinct arrangement.

**Local vs external (SD-13, SD-14).** Corpus-linear work that feeds statistics stays local; candidate-linear work (HHpred, DALI, AF3/ColabFold, Foldseek server) is external; crossover about 10^4 queries. Image `bio-box:3.1` adds HH-suite3, BLAST+, DIAMOND, CMfinder, R-scape, PADLOC, CRISPRCasTyper, CRISPRCasFinder, reseek and Domainator (§6). Domainator is the candidate workbench (network viewer, contig plots, structure annotation, GUI); bio-box exports to it with `--format domainator-gb` and does not build its own (SD-17). Search and report commands emit capped, schema-versioned NDJSON with `--json` (SD-18); hit-set networks come with threshold tables (SD-19); runs carry `max_output_gb` and `partial` guardrails (SD-20).

**Interactive lane (SD-15).** Build 32 threads and search 8 by day, 48/8 at night; search shards under 300 CPU-min each at nice 15; `/tmp` quota 150 GB.

**Tiers (SD-16).** Capabilities switch on by data quality, closed replicons first. A tier is certified when its replays pass and its cost table is measured on the pod.

| Tier | Stages | Adds |
|---|---|---|
| ST0 | S0 | runner, provenance, `protein_occurrence`, `nbr`, exhaustive search, IGR/dark-ORF/spacer corpora, Z policy, `--json` and schema rule (SD-18), guardrails (SD-20), `domainator-gb` export (SD-17) |
| ST1 | S1, S2 | NR90/NR50 pack, iteration controller, projection, spacer → protospacer, circular neighbourhoods, cross-tool concordance |
| ST2 | S2C | user-token promotion, on-demand statistics, `nbr score`, `domain_promiscuity`, local HH-suite, batch lane, hit-set networks with threshold tables (SD-19) |
| ST3 | S3, S5 | mandatory edge censoring and absence checks, quality-weighted counts, 3Di encoder bake-off (SD-21) then 3Di of NR50, B1 temporal holdout |
| ST4 | S4, S6 | search host on Hetzner, NR90 projection everywhere, scheduled DIAMOND blastx over contigs |
| ST5 | S7+ | compute-to-data, reduction-only searches, Logan-Search and Serratus as external services |

**Tests (family S, SEARCH_PLAN §7).** S-01 replay recall ≥0.90 (gate); S-02 cross-tool concordance (report); S-03 Z invariance; S-04 controls; S-05 provenance and failure recording; S-06 envelope; S-07 drift alarm; S-08 projection accuracy (≥99% agreement above cutoff + 10 bits); S-09 neighbourhood semantics; S-10 human dry run; S-11 B1 temporal leakage; S-02d Domainator as an independent implementation on the golden genomes (gate: exact agreement above cutoff + 10 bits, S0 and every container rebuild); S-12 recipe expressiveness (three Domainator cases, five synteny hypotheses; gate at S0, on paper); S-13 `--json` schema stability (gate from S0). They join the stage gates (§0.4 item 8).

**New search workflows (TEST_ARCHITECTURE §5).** A new search workflow enters as a recipe plus a workflow manifest: registered operations, outputs, positive and negative controls, truth source and budget class. The harness generates its schema, MR, control, truth-recovery, budget (S-06) and JSON-schema (S-13) tests; no harness code is written per workflow. A workflow that needs an unregistered operation is refused by CI until the operation's registry entry exists. The S-12 recipes are the first manifest set.

---

## 5. Validation and benchmarks (D-24, DS-04, DS-05, DS-10, DS-11, DS-15, DS-17, DQ-10, BX-20)

Every pass criterion is a rate with an exact-binomial (Wilson or Clopper-Pearson) 95% interval, and no criterion is written as 0% or 100%. For example, 48/50 has a Wilson interval of [0.868, 0.985]; 0/200 bounds failure at ≤1.5% (rule of three); 0/100 bounds it at ≤3%. Every benchmark input set is frozen and BLAKE3-hashed before the first scoring run.

**Sources of truth, in order of trust (TEST_ARCHITECTURE §1).** Planted-truth recovery is reported separately from snapshot agreement; only ranks 1–3 count as correctness evidence, and statistics need rank 1.

| Rank | Source | Benchmarks and tests | What a pass shows |
|---|---|---|---|
| 1 | Constructed truth: synthetic worlds with planted truth, and null worlds | B-SIM, B-NULL, manifest truth-recovery tests (TEST_ARCHITECTURE §5.2) | Planted features recovered at measured rates; nulls calibrated |
| 2 | Independent implementations | B-0 brute-force oracles, S-02d (Domainator), LY-07 (readers), pyhmmer vs `hmmsearch` | Two implementations agree |
| 3 | Published and curated results | S-01 replays, B-G curated annotations, B1 holdout | Known answers recovered |
| 4 | Golden snapshots | G-* hashes, G-05 drift | Change, not error |

| ID | Name | What it shows | Pass criterion | First stage |
|---|---|---|---|---|
| **B-L0** | Sequence round trip | The L0 store is lossless | For every shard, recomputed digests equal both the source checksum and the stored `contig_uid`/`protein_uid`/`igr_uid`. Random slices of (contig_uid, start, end, strand) match the source FASTA. | S0 |
| **B-0** | Engineering oracle harness | Indexes and tables return correct answers | Exact set equality with the brute-force unnested scan over 10⁴ stratified queries per class. Badge and descriptor round-trip property tests pass. Shuffled-input builds have identical hashes. Sizing ledger ×(scale) is within 20% of `sizing.yaml`. | S0 |
| **B-G** | Golden-genome truth checks | Annotation layers match curated truth (rank 3; G-* snapshot hashes are rank 4) | Per-layer targets from `critique_8` (golden set) and `critique_4` (truth table): genetic code ≥99% correct; promoted-smORF FDR <10% against Ribo-seq/Ribo-RET; frameshift and intron merges on lambda gpG-T and T4; operon boundaries against RegulonDB/SubtiWiki; domain partitions against ECOD/CATH manual; riboswitch subtypes against Rfam seeds. | S0 |
| **B-SIM** | Simulation with planted truth | Reported q means what it says | `synthworld` (ZOMBI2 backbone; `planning/tests/synthworld_spec.md`), about 5,000 genomes with planted modules, promiscuous passengers and NOGDs, then simulated fragmentation and MAG dropout. Z is N(0,1) on unplanted pairs (QQ slope 1.00 ± 0.02 over the central 99%). Observed FDR is within 2× of nominal on the reliability diagram (a reliability target with signal present; FDR control itself is gated on null worlds in B-NULL). Power curves are reported per regime. Planted passengers are not recovered above the stated FPR. | S1 |
| **B-NULL** | Null-world calibration (MR-18) | The null is not too easy (F-07) | `synthworld` null worlds: the same species tree, HGT and gene-order processes as B-SIM, no modules planted (a null world without phylogeny and HGT is not used). Null p-values uniform (KS test; QQ against permutation); observed FDR at nominal level within binomial-CI tolerance bands at fixed seeds; minimum detectable effect reported per regime. Covers O-STATS, O-ANCHORSTATS and O-NBRSCORE. | S0 (micro-corpus, ring R2, every merge); at stage scale in ring R5 at S2C, S5 and S6, where the statistics change scale (other tiers rerun at micro-corpus scale) |
| **B1** | Temporal holdout | The engine finds systems it was not told about | Build with `dict@2013`: Pfam 27, Rfam 11, COG 2014, TIGRFAM 15; PDB, Swiss-Prot, PaperBLAST and FEBA frozen 2013-12-31; no PADLOC, DefenseFinder, AcrDB, TADB, PHROGs, MIBiG, Breaker orphan set or post-2013 pocket rules. Score recall@15 for systems published after 2013, each with its date, reported with a Wilson CI. `dict@now` is shown as the ceiling and naive hypergeometric as the floor. Systems whose models are in the frozen set are excluded. PADLOC-DB candidate systems (HEC, PDC, including unpublished ones) and any other post-2013 models are excluded from the `dict@2013` build, and the exact DefenseFinder and PADLOC-DB versions are pinned per release with their candidate-system content checked at each update (test S-11, gate at ST3). The positive panel has ≥300 pairs, using MetaCyc/GapMind/MIBiG/PDB-assembly auto-generation. | S5 |
| **B2** | Masked-slot prediction | gLM project only | Deferred with §4.5. | none |
| **B3** | Fragmentation and MAG dropout | NOGD calls survive real assembly quality | Fragment 1,000 complete RefSeq chromosomes into 2.5 kb pieces and simulate 8% core-gene dropout. False-hole rate ≤1% (95% upper bound). Nominal q within 2× of observed. Uncensored P(B\|A) shifts <1%. | S6 |
| **B4** | Subfamily and pocket precision | Subfamilies do not mask paralogs | PF04055, PF00005, COG0720 and RF00442 resolve into the expected subfamilies and pocket subtypes. Rates are reported with CIs. | S5 |
| **B5** | Passenger and boundary traps | MGE cargo does not rank as coupling | 500 pre-registered traps, curated by someone not tuning the score, with an `UNKNOWN` class. At the threshold giving the B1 recall target, FPR on confident negatives is reported with a Wilson CI. Target: upper bound <2.3% at 5/500. Complemented by B-NULL: traps test named passengers; null worlds test the null as a whole. | S2 (plasmid subset), S4 (prophage subset), S5 (full) |
| **B6** | Chimera and pseudoenzyme closure decoys | Closure does not bleed across fused domains | 6A: 200 solved+dark chimeras, false-closure rate ≤1.5% (95% UB). 6B: 100 pseudoenzyme paralogs, ≤3% (95% UB) scored as closed. | S5 |
| **B-EXT** | External baselines | The engine adds something | STRING combined score, Coinfinder and phylogenetic logistic regression on shared genomes. AUPRC with cluster-bootstrap CI on the pre-registered sets. gLM2 belongs to the gLM project. | S5 |
| **B-PRO** | Prospective registry | Predictions hold up later | Commit a hash-timestamped file of the top 50 modules, each with a specific claim. Score after 12–24 months against new PDB assemblies, RB-TnSeq, DefenseFinder/PADLOC releases and third-party AF3/Boltz interface pTM. | First public release |

Ablation tables (naive hypergeometric, raw host-tree homoplasy, each null alone) are produced by running them, with CIs, not by asserting them.

---

## 6. Packaging, malleability and deployment

See [`DEPLOYMENT-AND-MALLEABILITY.md`](DEPLOYMENT-AND-MALLEABILITY.md). Constraints from this plan that it must honour:
- **Containers (D-26).** A CPU image runs as Apptainer on wilkcomp02 and as Docker on Hetzner. A separate CUDA image is used for the GPU stages. Weights ship in the data pack, not the image.
- **Search tools (SD-13).** Image `bio-box:3.1` adds HH-suite3, BLAST+ (including RPS-BLAST for the CDD route), DIAMOND, CMfinder, R-scape, PADLOC, CRISPRCasTyper, CRISPRCasFinder, reseek and Domainator (MIT; candidate workbench, SD-17); ESM-2 3B 3Di encoder in the CUDA image for the SD-21 bake-off. No ad hoc installs; each further tool needs an image release.
- **Test tooling (TEST_ARCHITECTURE §8).** pytest, Hypothesis, pandera, `coverage.py`, mutmut, asv, `zombi2==0.54.0` and the `synthworld` generator ship in a separate `bio-box-test:3.1-cpu` image layered on `bio-box:3.1-cpu`, built from the same lockfile; they are never in the runtime image. Tests therefore exercise the runtime image's exact binaries (`hmmemit`, `cmemit` and `hmmsearch` come from it). `test_run` rows record the runtime digest the tests ran against. Hosts per ring: CI host for R0–R2, wilkcomp02 for R3, stage host for R4–R5 (DEPLOYMENT-AND-MALLEABILITY.md §2).
- **Search pack (SD-06).** Ships per release next to the data pack; edition contents are in `DEPLOYMENT-AND-MALLEABILITY.md`.
- **Malleability contract (D-20).** Deltas are raw-searchable immediately. Statistics are "as of `stats_release_id`" and carry stale flags. `remap-token` merges additive statistics only. Parameter changes go through `bio-box rescore`.
- **Catalog and operations (D-21, D-22).** The catalog is read-only, generated from a manifest and swapped by atomic rename, with one writer. Object storage is the system of record. The API binds to localhost, behind TLS and a token, with parameterized endpoints.
- **Licensing (D-23).** Public packs contain only sources flagged `public` in `sources.yaml`. Everything else ships as bring-your-own-data recipes.
