# Staged accretion plan: from 14 genomes on wilkcomp02 to rented bare metal

Status: plan only. Nothing here has been run. Date: 2026-10-28; S2C, the search tiers and the test-harness gates added 2026-10-29.

This file is authoritative for stage numbering. Detail lives in:
- [`planning/sources_and_sizing.md`](planning/sources_and_sizing.md): sources, versions, licenses, unit-rate formulas, envelope checks.
- [`planning/data_layer_and_tests.md`](planning/data_layer_and_tests.md): layers, IDs, sequence store, test IDs (C-, L-, R-, Q-, A-, D-, G-, K-), orchestration.
- [`PLAN.md`](PLAN.md): pipeline and schema.
- [`DEPLOYMENT-AND-MALLEABILITY.md`](DEPLOYMENT-AND-MALLEABILITY.md): hosts, images, contract.
- [`planning/search/SEARCH_PLAN.md`](planning/search/SEARCH_PLAN.md): ad hoc search design (SD-01..SD-21), search pack, tiers ST0–ST5, S-tests (S-01..S-13, plus S-02d).
- [`planning/tests/TEST_ARCHITECTURE.md`](planning/tests/TEST_ARCHITECTURE.md): operation registry (O-*), failures F-01..F-15, MR-01..MR-19, workflow manifests, micro-corpus, synthetic and null worlds, LY-01..LY-10, rings R0–R5, `test_run`.
- [`audits/SYNTHESIS_2.md`](audits/SYNTHESIS_2.md): decision record D-01..D-28.

---

## 0. How this brief was rewritten

The request was: "define how to scale this up one piece at a time, starting with smaller genomes. Ultimate target is rapid deployment to rented bare metal (Hetzner); for now the target is wilkcomp02, CPU only. Build it in pieces over time; test DuckDB and the other data layers; accrete data in stages."

**1. What a top professional would do that was not asked for**
- Give each stage one scientific question it must answer against outside truth, so the stage tests the science and not only the plumbing.
- Fix the decisions that are expensive to reverse in S0, before any data accretes: stable IDs, the schema, provenance, the sequence store and the license flags.
- Make accretion provable. Building stage *k* incrementally on top of *k−1* must equal a from-scratch build of 0..*k*, hash for hash.
- Use each stage to measure unit rates on the real box and refit the sizing model, so the next stage's cost is predicted rather than guessed.
- Treat the host's house rules as hard constraints.
- Rehearse the move to rented hardware early, while the data is still small.

**2. The rawest material to start from**
- The measured wilkcomp02 envelope: `WILKECOMP.org` and the watchdog notes.
- Each source's own release files: GTDB R226 metadata and trees, the INPHARED `*_data.tsv`, PLSDB release files, IMG/VR v4 files, NCBI `assembly_summary_refseq.txt`.
- The golden genomes' FASTAs, with their curated truth: EcoCyc/RegulonDB, SubtiWiki, Rfam seeds, M-CSA.
- Tool binaries timed on this CPU.
- Real Parquet files built from S0 data, rather than extrapolated byte counts.

**3. Methods and tools**
- Workflow and environment: Snakemake local executor behind a `bio-box run-step` wrapper (systemd-run caps, nice 19, STOP file, staging on local `/tmp`). Apptainer `.sif` built from one OCI image, with `pixi.lock` and `uv.lock`.
- Annotation: pyrodigal / pyrodigal-gv / Codetta, pyhmmer, Infernal `cmsearch`, MMseqs2 ≥18 (`linclust`, `clusterupdate`), geNomad / CheckV, skani.
- Data layer: DuckDB with hive-partitioned Parquet; pandera or dbt-duckdb contracts; brute-force oracles; property tests (hypothesis).
- Measurement: `hyperfine` and `/usr/bin/time -v` microbenchmarks; DuckDB `EXPLAIN ANALYZE`.
- Sizing: power-law scaling fits across stages.

**4. What, other than us, shows the result is right**
- Curated truth sets the pipeline never trains on.
- Source-published counts reconciled by the ledger.
- Brute-force query results equal to indexed results.
- Rebuild-from-scratch hashes equal to accreted hashes.
- The same image and manifest reproducing a stage on a rented Hetzner box with identical hashes.

**Rewritten brief**
> Write a stage ladder S0..S7 for accreting `all-the-data`, starting from 14 golden genomes and growing through phages, plasmids, small genomes, viral collections, GTDB species representatives and all GTDB genomes, on wilkcomp02 (CPU only, within its measured envelope). For each stage give:
> 1. the data increment: named source, version, counts, bytes, license;
> 2. the one new capability and one science question, each with external truth;
> 3. the work broken into independently resumable chunks;
> 4. the data-layer and DuckDB tests that gate it, including sequence-store round-trip and accretion equivalence;
> 5. the unit rates it measures and feeds into the sizing model;
> 6. the go/no-go exit gate.
>
> Fix every hard-to-reverse decision in S0. Add a parallel portability track that rehearses Hetzner deployment on small stages and moves the global steps there when wilkcomp02's envelope runs out. Plan only; execute nothing.

---

## 1. Principles

1. **One stage = one data increment + one capability + one science question + one release.** A stage ends by cutting `release_id = rS<k>.<n>`. Nothing is half-merged.
2. **Chunks inside a stage are shards.** A shard is 500–2,000 genomes (or 1–4 GB of sequence). Every shard is idempotent, content-addressed and resumable, and pausing between shards is always safe (STOP file).
3. **The layers have different update rules** ([data layer §1](planning/data_layer_and_tests.md)):
   - Append-only: L0 sources, L0-seq sequence store, L1 per-genome annotations, L3 coordinates.
   - Incremental with a crosswalk: L2 clusters.
   - Recomputed every release: L4 statistics.
   - Built only if a measured SLO fails: L5 indexes.
4. **Raw sequences are never discarded.** Every source contig is kept in the sequence store, including dedup members and screening rejects (which go to a quarantine partition). The one exception is contigs classified as human or vertebrate host: for those only the digest and reason code are kept (preflight gap #16). Proteins are kept too; IGRs are kept as coordinates and sliced from the stored contig. Each is retrievable by `(uid, start, end, strand)` with lineage back to the source accession. Dedup changes counting weights only.
5. **Accretion must be provable.** Test A-01..A-07: S(k−1) plus the delta equals S(0..k) rebuilt from scratch. L0, L1 and L3 must match by hash; L2 must match up to the crosswalk with ARI ≥ 0.99; L4 must match after recompute.
6. **Measure, then predict.** Every stage writes measured rates into `planning/sizing.yaml`. The next stage's cost estimate must come from that file, and a stage whose measured cost misses its estimate by more than 2× stops for review.
7. **House rules are constraints, not etiquette.**
   - Threads: 40 by day, 56 at night.
   - Memory: 200 GB cgroup per job.
   - Monolithic tools run at nice 19.
   - No small files on NFS; staging happens on local `/tmp`.
   - Durable small outputs go to `/stor/work`; bulk outputs go to `/stor/scratch`.
8. **No public service on wilkcomp02.** `bio-box serve` binds to localhost and is reached over an SSH tunnel.
9. **Search capability scales by data quality, closed high-quality nucleotide data first.** Search tiers ST0–ST5 ([SEARCH_PLAN §6](planning/search/SEARCH_PLAN.md)) switch on with the stages: closed golden genomes, then complete phages and plasmids, then complete prokaryotic chromosomes (S2C), and only then drafts, MAGs and metagenomic viruses. A tier is certified when its replays pass on its own data and its cost table is measured on the pod. Capabilities that matter only for fragmented assemblies (edge censoring, gene-call rescue for absence claims) become mandatory when fragmented data arrives (ST3).
10. **Tests attach to operations; workflows inherit coverage.** Every workflow is a manifest over registered operations (O-*), each with a contract, an oracle and metamorphic relations. A new workflow needs only a manifest and a recipe; a new operation needs a registry entry first. Truth is ranked: constructed truth (synthetic and null worlds), then independent implementations, then published results, then golden snapshots, which catch change, not error. Every gate is a SQL query over `test_run` ([TEST_ARCHITECTURE §1, §5, §8](planning/tests/TEST_ARCHITECTURE.md)).

---

## 2. Stage ladder at a glance

Estimates come from [`sources_and_sizing.md`](planning/sources_and_sizing.md) §4. All rates are unmeasured until S0, and wilkcomp02 delivers about 800 core-h/day.

| Stage | Adds | n | Gbp (est.) | Core-h (est.) | Seq store / other disk | License | Runs on | Search tier |
|---|---|---|---|---|---|---|---|---|
| S0 | Golden set (~14 genomes) + CAMI II subset; harness | ~14 + subset | 0.05 (+1–3) | 100–1,500 | 0.03 (+1) GB / <20 GB | open | wilkcomp02 | ST0 |
| S1 | INPHARED complete phages | 30–40k | 1.7–2.2 | 250–2,500 | 0.9 GB / ~15 GB | open (scripts AGPL) | wilkcomp02 | ST1 |
| S2 | PLSDB plasmids | 72,360 | 5–7 | 800–6,000 | 2.8 GB / ~40 GB | check before redistribution | wilkcomp02 | ST1 (certified) |
| S2C | GTDB R226 complete prokaryotic genomes (`Complete Genome`, ≤20 per species) | est. 15–30k (pin at S−1) | 60–120 | 3,000–22,000 | 28–55 GB / ~80–240 GB | INSDC + GTDB | wilkcomp02 | ST2 |
| S3 | GTDB R226 species reps <2 Mbp | 9–12k | 10–13 | 1,000–10,000 | 5.2 GB / ~60 GB | GTDB CC BY-SA (verify) | wilkcomp02 | ST3 (deferred to S5) |
| S4 | IMG/VR v4 HQ vOTU reps (+ IMG/PR) | 0.23M (+0.7M) | ~12 (+10–30) | 4,500–17,000 | 5.5–19 GB / ~150 GB | JGI policy: local-only | wilkcomp02 | ST4 (deferred to S6) |
| S5 | All GTDB R226 species reps (remaining after S2C and S3) | 143,614 total | 400–450 | 20,000–80,000 | 207 GB / 0.6–0.8 TB | INSDC + GTDB | wilkcomp02, borderline (§4) | ST3 (certified) |
| S6 | All GTDB R226 genomes (remaining after S2C) + RefSeq complete flag | 732,475 total | 2,700–3,100 | 60k–1.6M | +1.34 TB / +2–4 TB | INSDC | split: per-genome work on wilkcomp02, global work on Hetzner | ST4 (certified) |
| S7+ | AllTheBacteria, SPIRE, GlobDB extras, then Logan | millions | Tbp–Pbp | — | — | mixed | Hetzner; Logan on AWS next to its bucket | ST5 |

S2C figures are *derived* in [`sources_and_sizing.md`](planning/sources_and_sizing.md) §4 from the S5 rate and the 0.46 GB/Gbp store model; the count is an *estimate* until the S−1 metadata count. S5 and S6 keep their totals; their per-stage deltas shrink by what S2C (and, for S5, S3) already loaded. Search tiers are from [SEARCH_PLAN §6](planning/search/SEARCH_PLAN.md); "deferred" means the stage loads data whose search acceptance belongs to a later tier gate.

The portability track (§5) runs alongside: P0 at S0, P1 after S2, P2 at S5/S6, P3 at S7.

---

## 3. Stage cards

Each card lists six things:
- **Adds:** the data increment.
- **Capability:** the new pipeline or schema capability.
- **Science check:** the one question this stage answers, judged against outside truth.
- **Chunks:** resumable work units, in order.
- **Gate:** blocking tests ([data layer §4](planning/data_layer_and_tests.md)), measurements written to the sizing model, and the go/no-go call.
- **Search gate:** the search tier, its replays and the S-tests that apply ([SEARCH_PLAN §6–7](planning/search/SEARCH_PLAN.md)). Replay recall (S-01) is scored only on published hits whose genomes are loaded at that stage.

### S−1: Preflight (no corpus data, no pod compute; CI-scale fixtures only)

- **Why:** a pre-mortem ([`planning/preflight_gaps.md`](planning/preflight_gaps.md)) found 22 gaps (29 after the search and test rounds) and 10 cross-document contradictions that would stall S0–S3. Examples: no signed-off first questions, no repo or CI, FCS-GX needing ~470 GB of RAM against the 200 GB cgroup, no pod footprint request, L0 kept only on unbacked scratch, reference DBs not pinned, no kill criteria.
- **Work:**
  1. Sign off the first-10 questions, each with an existing-tool baseline.
  2. Assign roles, time boxes and kill criteria K1–K7.
  3. Send the footprint request to the pod owner, covering wilkcomp03 GPU access and a second account.
  4. Commit the ADRs: IDs (BLAKE3-128), the FCS-GX route (**decided 2026-10-29**: deferred and accreted later as a QC sidecar when ≥512 GB RAM is available; human/vertebrate host screen from S0; `preflight_gaps.md` #4), the human-data rule, the upgrade policy, L0 durability, and the license and name.
  5. Set up the repo and CI with a micro-S0 fixture.
  6. Pull the `.sif` by digest and record `bio-box doctor`.
  7. Port `run-step` from `homology_pod.sh`.
  8. Fill and freeze `golden.yaml`, `refdb.yaml` and `sources.yaml`.
  9. Write `planning/costs.md`.
  10. Count complete genomes from the GTDB R226 metadata (`bac120_metadata_r226.tsv`, `ar53_metadata_r226.tsv`, field `ncbi_assembly_level = Complete Genome`), split into species reps and members, and apply the ≤20-per-species selection. This pins the S2C count, Gbp and cost, and the S5/S6 remainders.
  11. Fetch the replay supplements not yet in `lit/supp/` ([SEARCH_PLAN §9.5](planning/search/SEARCH_PLAN.md)): Doron 2018 S1–S3, Gao 2020 S1–S3, Altae-Tran Data S1–S4, Millman Table S1, Brewer 2021 files 1–2, CRISPR-Cas 2025 Zenodo PSSMs, Domainator Zenodo.
  12. Confirm the Z constants as an ADR (SD-04): `Z_prot = 1e8`, `domZ = 1e8`, `Z_nt = 1e5` Mb, fixed for the release line.
  13. Fix the `bio-box:3.1` tool list (SD-13, SD-17): `3.0` plus HH-suite3, BLAST+, DIAMOND, CMfinder, R-scape, PADLOC, CRISPRCasTyper, CRISPRCasFinder, reseek and Domainator; ProstT5 and the ESM-2 3B 3Di encoder in the CUDA image only. Write the three Domainator case studies and five synteny hypotheses as recipes on paper (S-12).
  14. Record the synteny corpus repair as a dependency ([SEARCH_PLAN §9.4](planning/search/SEARCH_PLAN.md)): about 46 `lit/` files hold a different paper than their name. Replays use supplements and are not blocked; hypotheses drawn from mislabeled texts are re-extracted before they drive experiments.
  15. Specify the test harness before the layout ([`planning/tests/TEST_ARCHITECTURE.md`](planning/tests/TEST_ARCHITECTURE.md)): the schema registry, the operation registry with every O-* entry's contract, oracle and MRs, the workflow manifests for the S-12 recipes, the synthetic-world generator spec (written: `planning/tests/synthworld_spec.md`, ZOMBI2 backbone plus custom code; S−1 fit spike in its §1.4), the `test_run` schema, and the CI host for ring R2. Layout choices still open (K-13, row-group sizes, sort keys) are decided at S0 by LY-06 row-group counts on the micro-corpus.
- **Gate:** the checklist in `preflight_gaps.md` §7 is all green (its test-harness group is item 15), plus items 10–13 and 15 done and item 14 owned. Time box: 3–4 weeks (*estimate*; the earlier 2–3 weeks predates item 15).

### S0: Golden genomes, harness, and the decisions that cannot be revisited cheaply

- **Adds:**
  - The ~14 golden genomes from critique_8 (16 entries with the CAMI II subset and the synthetic spike-in panel): MG1655, *B. subtilis* 168, T4, λ, *M. genitalium*, an SR1 genome, a Lak phage, *M. acetivorans*, *H. volcanii*, *S. acidocaldarius*, H37Rv, *B. thetaiotaomicron*, *S. pyogenes* SF370, an Obelisk. Accessions marked "pin" get chosen here.
  - A CAMI II marine subset.
- **Capability:**
  - The full layer stack at toy scale: L0 contracts (`sources.yaml`), the L0-seq sequence store, L1 Parquet, L2 clusters, L3 tapes and `slot_pairs`, L4 score-test statistics (`pair_stats`), and the read-only DuckDB catalog swapped in by atomic rename.
  - The `bio-box run-step` wrapper, the ledger, and `bio-box doctor`.
- **Decisions frozen here (ADRs):**
  - ID scheme: blake3 content UIDs + dense per-release IDs + `id_history`.
  - Schema v1, including the corrected badge layout, int32 `delta_bp` and namespaced token IDs.
  - Enum registry mapped to Sequence Ontology terms.
  - Sequence store format: Parquet chunks vs bgzip+faidx, decided by bake-off test K-13, scored by LY-06 row-group counts.
  - Partition keys, file sizes, row-group sizes and sort keys (LY-04, LY-06).
  - Release naming.
  - License flags.
- **Science check:** does the pipeline reproduce curated annotations?
  - Gene starts vs EcoCyc / SubtiWiki.
  - Operon boundaries vs RegulonDB transcription units.
  - T4 introns and HEGs stitched; the λ gpG-T frameshift merged.
  - *M. genitalium* code 4 and SR1 code 25 detected by Codetta.
  - Leaderless starts in *H. volcanii*.
  - Riboswitches vs Rfam seeds.
  - Full-Rfam vs production-mode Infernal recall ≥ 0.98.
- **Chunks:**
  1. Build the image, lockfiles and `.sif` (portability step P0).
  2. Build the harness in TEST_ARCHITECTURE §10 order: registries, `synthworld` with truth tables and null worlds, the micro-corpus, oracles, MRs, the S-12 manifests, LY-01..LY-10 on candidate layouts, fault injection, `test_run`.
  3. Fetch the golden set into L0 with checksums.
  4. Run the microbenchmarks: the seven measurements in sizing §6.
  5. Run every annotation step per genome.
  6. Build L2–L4.
  7. Write the oracle implementations of query classes Q1–Q8: brute-force Python over the L1 Parquet.
  8. Run the full test suite (R0–R3).
  9. Cut release `rS0.1`.
- **Gate:**
  - Blocking tests: C-01..C-08, C-10..C-12, L-01..L-03, R-01..R-07 and R-09 (sequence-store round trip: per-record digest, 10k random slices vs `samtools faidx`, minus strand, chunk edges, extract round trip), Q1–Q8 oracle equivalence, D-01..D-04 determinism (thread count, kill-resume, input order), G-01, G-02, G-04, and K-01..K-13 (DuckDB version pin, memory_limit under the cgroup, temp_directory on local disk, thread cap, partition and row-group pruning via EXPLAIN ANALYZE, LIST cost, 1–16 read-only clients, fail-fast on a second writer, NFS vs local).
  - Harness ([TEST_ARCHITECTURE](planning/tests/TEST_ARCHITECTURE.md)): micro-corpus and synthetic worlds (with null worlds) built with truth tables; rings R0–R3 green; LY-01..LY-10 run, and LY-06 decides K-13, row-group sizes and sort keys; first mutation-score measurement (target ≥ 80%, *target*); first null-world calibration of O-STATS (MR-18: uniform p-values, FDR at nominal).
  - Measurements: `R_pyro`, `r_pfam`, `s_cm`, the geNomad rate, MMseqs2 `l` with RSS and tmp bytes, bytes per bp in the sequence store, and DuckDB scan rates.
  - Go/no-go: a `test_run` query returns every blocking test above (and the search gate below) passing on `rS0.1` and the pinned image digest, and the sizing model is refit. If a measured rate pushes S5 beyond 150 days, re-plan before S1.
- **Search gate (ST0):**
  - Switched on: recipe runner and provenance, `protein_occurrence`, `bio-box nbr`, exhaustive search on all proteins (no reduction), IGR, dark-ORF and spacer corpora, six-frame `region translate`, Z policy, cost measurement.
  - Replays: Makarova 2011 defense islands (golden subset); Wolf 2012 arCOG concordance; Weinberg 2010 motifs in *B. subtilis*, *B. theta*, SF370; Shmakov 2019 Set2 (*H. volcanii* seeds); Axen 2014 BMC loci (MG1655 eut/pdu); SYN-143 and SYN-125 end to end.
  - S-tests: S-01, S-02d, S-04, S-05, S-06, S-07, S-09, S-12 (on paper), S-13. S-03 starts at S1 (it compares adjacent stages); S-08 starts at ST1 (no reduction here).

### S1: Complete phages (INPHARED)

- **Adds:** about 30–40k complete phage genomes. The exact count is read from INPHARED's `*_data.tsv` at freeze.
- **Capability:**
  - pyrodigal-gv calling with codes 11 and 15.
  - PHROGs.
  - Intron/HEG and frameshift stitching at scale.
  - Viral regime classification with geNomad.
  - The first real `slot_pairs` and `pair_stats` tables.
  - The first L2 cluster dictionary, with `clusterupdate` exercised.
- **Science check:** recover known phage modules as top-ranked modules, using defense, anti-defense and other annotations held out from the tokens:
  - lysis cassettes (holin, endolysin, spanin);
  - Acr–Aca loci;
  - the tail assembly chaperone frameshift.
- **Chunks:**
  1. Shard by INPHARED host-genus bins.
  2. Run L1 per shard.
  3. Run L2 `linclust`, then `clusterupdate` against S0.
  4. Build L3 and L4.
  5. Run the tests and cut release `rS1.1`.
- **Gate:**
  - Blocking: everything S0 runs, plus **A-01..A-07 accretion equivalence (S0 → S1 incremental vs S0+S1 from scratch)**.
  - Recorded: the first scaling-fit point and R-08 retrieval latency.
  - Harness: gate = R4 on `rS1.1`, as a `test_run` query; new workflows arrive with manifests (TEST_ARCHITECTURE §5.1).
- **Search gate (ST1, first half):**
  - Switched on: NR90/NR50 targets and the search pack, the iteration controller at scale, projection and re-score, spacer → protospacer, circular neighborhoods, first cross-tool concordance.
  - Replays: Yutin 2024 polymerases (23,246 proteins, 18,383 phage genomes); Chaban 2024 catalytic-motif column; Silas 2025 (phage loci); Peters 2017 spacers (phage protospacers).
  - S-tests: everything S0 runs, plus S-02 (report), S-03 (S0 → S1) and S-08.

### S2: Complete plasmids (PLSDB)

- **Adds:** 72,360 plasmids.
- **Capability:**
  - Plasmid backbone vs cargo.
  - MOB/Rep/conjugation calls (MacSyFinder CONJScan, as a sidecar).
  - IntegronFinder 2 `attC` calls.
  - Vehicle clustering on real MGEs.
  - The replicon-matched geometry test for Tier 2 (D-13).
- **Science check:** the score separates real coupling from passengers.
  - Couplings that should rank: known toxin–antitoxin and partition systems (*parAB*, *ccdAB*, *hok/sok*).
  - Passengers that should not: AMR–transposon pairs (`blaTEM`–`tnpA`).
  - This is the first, pre-registered and hash-committed, portion of the Benchmark 5 trap panel. FPR is reported with an exact-binomial CI.
- **Data-layer decisions:**
  - DuckLake spike against the §3.3 criteria.
  - K-07: LIST vs exploded tape.
  - K-12: whether an offset index is needed.
- **Chunks:**
  1. Run L1 by PLSDB host-taxon shard.
  2. Run L2 `clusterupdate`.
  3. Recompute L3 and L4.
  4. Score the trap panel and cut release `rS2.1`.
- **Gate:**
  - Blocking: everything S1 runs, plus the three data-layer decisions recorded as ADRs.
  - Recorded: a 2-point scaling fit.
  - Harness: gate = R4 on `rS2.1` plus R5 for ST1 certification, as a `test_run` query; new workflows arrive with manifests.
- **Search gate (ST1, certified here):**
  - Replays: Faure 2023 Tn7 (plasmid subset); Peters 2017 spacers (plasmid protospacers); the S1 replays re-run on `rS2.1`.
  - S-tests: everything S1 runs (S-03 on S1 → S2), plus S-10, the human dry run, once for ST1. The ST1 cost table is measured on the pod.
- **Portability:** P1 runs after this gate (§5).

### S2C: Complete prokaryotic genomes (GTDB R226, assembly level `Complete Genome`)

Inserted and approved 2026-10-29 ([SEARCH_PLAN SD-16, §6.1](planning/search/SEARCH_PLAN.md)): ≤20 complete genomes per GTDB species, chosen by CheckM2 quality then ANI diversity. Releases are `rS2C.<n>`; no other stage is renumbered. Optional follow-on **S2C-b** (the remaining complete genomes, uncapped; est. 7,500–54,000 core-h for all complete genomes together) is appended after the S2C gate only if the signed-off first-10 questions include strain-level work (defense hotspots, pangenome turnover); it costs an L2/L4 recompute and exercises accretion at scale.

- **Adds:**
  - GTDB R226 genomes, representatives and members, whose `ncbi_assembly_level` is `Complete Genome`, at most 20 per species, chosen by CheckM2 completeness and contamination, then by ANI diversity (skani).
  - Count: *estimate* 15–30k assemblies, pinned by the S−1 metadata count.
  - Complete species reps enter here and are not re-ingested at S5 (nor at S3 if under 2 Mbp). Complete non-rep members enter here and are not re-ingested at S6.
  - Plasmid replicons of these assemblies that PLSDB already supplied at S2 dedup by `seq_uid`; they gain the chromosome context, not a second copy.
- **Capability:**
  - The first closed prokaryotic chromosomes at scale: exact neighborhoods with no contig-edge censoring.
  - Multi-strain species (up to 20 genomes): dependence-aware counts (SD-12) report raw, vehicle-independent, species-rep and arrangement counts.
  - Search (ST2): user-token promotion (`bio-box search promote`) and on-demand anchor statistics against `background_rates`; `nbr score` with editable positive sets; `domain_promiscuity`; Rosetta and architecture queries at scale; local HH-suite (Pfam-HHM, PDB70); the batch query lane.
  - Library scans still run on exact-unique proteins (SD-05); NR90 projection takes over at S3.
- **Science check:** do score-test neighborhoods on closed chromosomes recover published defense islands and CRISPR-linked gene sets? Truth: Makarova 2011 PD2 islands (full) and Shmakov 2019 Sets 1–3 (`Seeds.tsv` → `Vicinity.tsv`), recall ≥ 0.90 on hits in loaded genomes (S-01), every miss classified.
- **Chunks:**
  1. Select genomes from the R226 metadata (the S−1 list, frozen in `sources.yaml`).
  2. Fetch via NCBI datasets in shards.
  3. Run L1.
  4. Run L2 `clusterupdate`.
  5. Build L3 and L4, then the search pack (`bio-box pack build`).
  6. Run the replays and cut release `rS2C.1`.
- **Gate:**
  - Blocking: everything S2 runs, plus A-01..A-07 accretion equivalence (S2 → S2C).
  - Recorded: a 3-point scaling fit (S1, S2, S2C). At 60–120 Gbp, S2C is the largest stage before S5, so its measured rate is the best predictor of S5.
  - Harness: gate = R4 on `rS2C.1` plus R5 for ST2 certification, as a `test_run` query; new workflows arrive with manifests; null-world calibration of O-STATS (MR-18) repeated at S2C scale (kill criterion K7 in `preflight_gaps.md` §6).
- **Search gate (ST2):**
  - Replays: Makarova 2011 defense islands (full); Shmakov 2019 Sets 1–3; Shmakov 2018 CRISPR islands; Wu 2026 Druantia/Zorya (RefSeq loci); Bell 2026 YprA; Jimmy 2020 toxSAS; Atkinson 2011 RSH; the Mutte all-models score matrix.
  - S-tests: S-01 and S-03..S-09 blocking (S-03 on S2 → S2C); S-02 cross-tool concordance reported, including HH-suite on remote cases; S-10 once for ST2. The ST2 cost table is measured on the pod (S-06).
- **Estimates:** *est.* 15–30k assemblies, 60–120 Gbp (at ~4 Mbp per assembly); 3,000–22,000 core-h (*derived*: 60–120 Gbp × 50–180 core-h/Gbp, the S5 rate; SEARCH_PLAN §6.1 quotes 7,500–20,000, narrower than the product of its own inputs); sequence store 28–55 GB (0.46 GB/Gbp); other disk ~80–240 GB (S5 ratio); search pack ~50–105 GB ([sizing §4](planning/sources_and_sizing.md)). Up to ~27 days at 800 core-h/day at the high end.
- **Runs on:** wilkcomp02. MMseqs2 linclust over 6–12×10⁷ proteins needs ~20–40 GB (336 B/seq), inside the 200 GB cgroup.
- **License:** INSDC (open) for genomes; GTDB CC BY-SA (verify) for taxonomy.

### S3: Small-genome bacteria and archaea (GTDB R226 species reps under 2 Mbp)

- **Adds:** about 9–12k species reps from Patescibacteria, DPANN, Mycoplasmatota and endosymbionts, selected by `genome_size` in the R226 metadata. Reps already loaded at S2C (complete assemblies) are not re-ingested.
- **Capability:**
  - Codetta-driven genetic code calls at scale.
  - Taxonomy assignment wired into `tax_order_id` with confidence (D-06).
  - The first tree-ordered presence bitmaps (sparse) and transition edges on the R226 tree.
  - The smORF promotion rule exercised.
- **Science check:** is the missing-step search sane on reduced genomes?
  - Known pathway losses in endosymbionts, such as *Buchnera* amino-acid pathways, are called as losses and not as missing-step holes needing a replacement.
  - Reassigned codes in BD1-5/SR1 are detected.
- **Chunks:**
  1. Select genomes from the metadata.
  2. Fetch via NCBI datasets in shards.
  3. Run L1.
  4. Run L2 `clusterupdate`.
  5. Build the tree-ordered L4 and cut release `rS3.1`.
- **Gate:**
  - Blocking: everything S2C runs, with oracle and accretion checks on a stratified sample.
  - Recorded: a 4-point scaling fit (S1, S2, S2C, S3), which drives the S5 go/no-go.
  - Harness: gate = R4 on `rS3.1`, as a `test_run` query (no R5: ST3 certifies at S5); new workflows arrive with manifests.
- **Search gate (ST3, deferred):** search acceptance for these small-genome reps is part of ST3, certified at the S5 gate ([SEARCH_PLAN §6](planning/search/SEARCH_PLAN.md)); they load here for genetic-code work. At S3 only the regression checks run: S-03 (S2C → S3), S-04, S-05, and the ST0–ST2 replays must not lose recall. Library scans move from exact-unique proteins to NR90 reps with projection (SD-05), so S-08 runs on the S3 members.

### S4: Uncultivated viruses and plasmids (IMG/VR v4 HQ, optionally IMG/PR)

- **Adds:** 231,408 HQ vOTU representatives, optionally plus about 0.7M IMG/PR plasmids.
- **Local-only:** these rows carry `redistributable = false`, and no public pack contains their sequences; packs ship fetch recipes instead.
- **Capability:**
  - Host predictions (from IMG/VR plus iPHoP or spacer links) stored with their confidence.
  - Provirus boundaries.
  - Spacer–protospacer links.
  - An L2 drift test, comparing ARI of incremental vs full recluster on a sample.
- **Science check:** anti-defense hotspots and viral counter-defense modules rank above prophage moron–structural passengers. This is the prophage-boundary portion of the trap panel, with a CI.
- **Chunks:**
  1. Ingest with IMG/VR's shipped geNomad/CheckV calls; recompute them only on a calibration sample.
  2. Run L1.
  3. Run L2.
  4. Run L3 and L4, then cut release `rS4.1`.
- **Gate:**
  - Blocking: everything S3 runs, plus the drift test and the license-flag test (no non-redistributable UID in any exportable pack, search pack included).
  - Harness: gate = R4 on `rS4.1`, as a `test_run` query (no R5: ST4 certifies at S6); new workflows arrive with manifests.
- **Search gate (ST4, deferred):** S4 enters the ladder before S5 for data reasons; its search acceptance waits for ST4 at the S6 gate ([SEARCH_PLAN §6](planning/search/SEARCH_PLAN.md)). At S4 only the regression checks run: S-03 (S3 → S4), S-04, S-05, and no loss of ST0–ST2 replay recall.

### S5: All GTDB R226 species representatives (the last stage that fits wilkcomp02 in full)

- **Adds:** the remaining species reps: 143,614 in total, minus the S3 reps (9–12k) and the complete reps already loaded at S2C (count from S−1). That is ~132k before the S2C subtraction. Complete reps are not re-ingested.
- **Capability:**
  - The full token dictionary (`seq_family`, `subfamily_50`).
  - Tree-wide dispersion.
  - The core pangenome plane over species reps.
  - B1 temporal holdout first run: the dictionary frozen at a 2013 snapshot, scored only on systems published after 2013.
  - Decision per query class on whether any L5 index is needed.
- **Science check:**
  - B1 recall with a CI.
  - The external baselines STRING and Coinfinder on the same pairs.
  - The first hash-timestamped prospective prediction registry entry.
- **Envelope (sizing §5):** borderline on memory (MMseqs2 ≈160 GB), on local `/tmp` (0.5–0.9 TB) and on wall clock (25–100 days). Before starting S5, do these in order:
  1. Exact-sequence dedup before clustering.
  2. Per-phylum `linclust` merged with `clusterupdate`.
  3. HMMER and Infernal in shards that stay under ~300 CPU-minutes per process.
  4. A written overnight memory exception from the pod owners if still needed.
  5. If the S3 scaling fit predicts more than 150 days, run S5's global steps (L2) on a rented box (P2), not on the pod.
- **Chunks:**
  1. Run L1 in about 100 shards over weeks, with STOP checks between shards.
  2. Run L2 per phylum, then merge.
  3. Run L3 and L4, then cut release `rS5.1`.
- **Gate:**
  - Blocking: everything S4 runs.
  - The SLO table measured warm and cold, p50/p95, for every query class.
  - The decision to build L5 or not, recorded as an ADR.
  - Harness: gate = R4 on `rS5.1` plus R5 for ST3 certification, as a `test_run` query; new workflows arrive with manifests; null-world calibration of O-STATS (MR-18) repeated at S5 scale.
- **Search gate (ST3, certified here; covers the S3 reps):**
  - Switched on: edge censoring mandatory; dark-ORF and frameshift checks mandatory for absence claims; quality-weighted counts; ProstT5 3Di of NR50 (GPU, wilkcomp03 or rented); targeted CMfinder; B1 temporal holdout.
  - Replays: Weinberg 2017 motifs; Makarova 2019 archaeal dark matter; Makarova 2024 Asgard ESCRT; Kaur 2022 death-domain conflict systems; Mestre 2020 RTs; Wu 2026 (full); Corbino 2005.
  - S-tests: S-01, S-03..S-09 (S-03 on S4 → S5); S-02 reported; S-10 once for ST3; S-11 temporal leakage (the B1 `dict@2013` build holds no PADLOC-DB candidate systems or other post-2013 models). The ST3 cost table (§5 of SEARCH_PLAN) is measured on the pod.
- **Packaging:** this release becomes the `lite` edition (§6), with the S5 search pack.

### S6: All GTDB R226 genomes + RefSeq complete flag (split across hosts)

- **Adds:**
  - The remaining non-representative members: 588,861 minus the complete members already loaded at S2C (count from S−1), 732,475 genomes in total. Complete members are not re-ingested.
  - The `assembly_level` flag from `assembly_summary_refseq.txt`.
- **Capability:**
  - The three-valued pangenome planes with shrinkage (D-15).
  - The `hq_closed_or_multi_strain` mask.
  - Strain-level vehicle dedup.
  - Benchmark B3, the fragmentation and MAG-dropout robustness test.
- **Split:**
  - On wilkcomp02: per-genome L1, the sequence store, and the per-species pangenome work, all embarrassingly parallel.
  - On Hetzner (P2): global L2 clustering and geNomad on new members.
  - Moving between hosts: L1 and the sequence-store shards move by manifest, and the results come back as L2 Parquet with crosswalks.
- **Science check:** missing-step search with the quality masks reproduces known non-orthologous gene displacements (ThyA/ThyX, MetE/MetH, FolB/PTPS-III) with CIs. B3's false-hole rate is under its pre-registered bound.
- **Gate:**
  - Blocking: everything S5 runs, plus D-05 cross-host hash equality.
  - Recorded: cost per stage on Hetzner.
  - Harness: gate = R4 on the S6 release (on both hosts) plus R5 for ST4 certification, as a `test_run` query; new workflows arrive with manifests; null-world calibration of O-STATS (MR-18) repeated at S6 scale.
- **Search gate (ST4, certified here; covers S4):**
  - Switched on: the search pack exceeds the pod, so search runs on a Hetzner search host (≥512 GB RAM, ≥4 TB NVMe, *estimate*; sized at the S5 gate, §8); NR90 projection everywhere; batch lane with query batching; DIAMOND blastx over contigs as a scheduled job.
  - Replays: Gussow 2020 anti-CRISPR; Silas 2025 (full); Yutin 2018 crAss (S4 data); a FLSHclust-style CRISPR-proximity cluster score on a subset.
  - S-tests: everything S5 runs (S-03 on S5 → S6, and across hosts), plus S-10 once for ST4.
- **Packaging:** this release becomes the `core` edition, with the S6 search pack on Hetzner.

### S7+: Tier 3, then Logan (off wilkcomp02 by design)

- **Order:**
  1. AllTheBacteria (2.44M genomes, MIT).
  2. SPIRE and GlobDB extras.
  3. Logan last.
- **Logan:** compute-to-data on AWS us-east-1. Accessions are preselected by SRA metadata and STAT profiles, and only the reduced outputs (L1 Parquet plus the sequence-store shards for kept contigs) leave AWS.
- **Gate (at every substage):**
  - Identity resolution against S0–S6.
  - Screening: GUNC, CheckM2, human/vertebrate host screen; FCS-GX added when ≥512 GB RAM is available (deferred 2026-10-29, `preflight_gaps.md` #4).
  - A rarefaction check: novel NR50 clusters per Tbp stays above a pre-set threshold, otherwise stop.
  - Harness: R4 at every substage plus R5 for ST5, as a `test_run` query; new workflows arrive with manifests.
- **Search gate (ST5):** compute-to-data; reduction-only searches (cluster, then score by context); Logan-Search and Serratus as external nucleotide services. Replays: Altae-Tran 2023 FLSHclust subset; Domainator 2024 workflows. S-tests: S-01, S-03..S-06 and S-08 at every substage.

---

## 4. Envelope limits and when to move hosts

| Limit | First stage that strains it | Planned response |
|---|---|---|
| 200 GB memory cgroup | S5 (MMseqs2 ≈160 GB) | Exact dedup first; per-phylum `linclust`; written exception; otherwise move L2 to P2 |
| ~940 GB local `/tmp` | S5 (0.5–0.9 TB peak) | `--split-memory-limit`, shard-wise passes, delete tmp between passes |
| Watchdog (~369 CPU-min per PID) | Every MMseqs2 run from S1 on | Launch at nice 19; split HMMER/Infernal into short shard processes |
| Shared NFS (10 GbE) | S6 (~2 TB through NFS) | One throttled stream (`ionice -c3`); bulk sync during off-hours |
| `/stor/scratch` 30 TB, not backed up | S6 (sequence store passes 1 TB) | Second copy of the sequence store on object storage or the lab storage server; `seq_index` and manifests on `/stor/work` |
| CPU-years for geNomad at S6 scale | S6 | Move to Hetzner; use IMG/VR- and GTDB-shipped calls where they exist |

---

## 5. Portability track (towards rapid deployment on rented bare metal)

| Step | When | What | Gate |
|---|---|---|---|
| P0 | S0 | One OCI image → Apptainer `.sif` (pod) and Docker (Hetzner); `pixi.lock` + `uv.lock`; `bio-box doctor` host probe; reference DBs as L0 sources, not in the image | Image digest pinned; S0 runs from the `.sif` only |
| P1 | after S2 | Rent one Hetzner AX-class box for a few days. Run `bootstrap.sh` (OS, RAID, Docker), pull the image by digest, pull the `rS2.1` pack from object storage, then rebuild S0–S2 from manifests | D-05 cross-host hash equality on L0/L1/L3 and L2 via crosswalk; K-tests re-run on the new host; wall time from rental to first query recorded (target: under 2 h for S0–S2) |
| P2 | S5/S6 | Hetzner runs L2 global clustering and geNomad. Object storage becomes the system of record for packs | Round-trip pod → Hetzner → pod with identical manifests |
| P3 | S7 | AWS compute-to-data for Logan; Hetzner serves `planetary` | Cost per novel cluster within budget |

"Rapid deployment" is measured, not asserted: P1 and P2 record the minutes from bare rental to first correct query, broken into OS bootstrap, image pull, pack download (at 1 Gbps and 10 Gbps), catalog open, and first oracle-checked query.

---

## 6. Editions as stage cuts

| Edition | Stages | Notes |
|---|---|---|
| `lite` | S0–S5 (including S2C) | IMG/VR rows (S4) are present as statistics only; their sequences ship as fetch recipes. Ships the S5 search pack ([SEARCH_PLAN §3.2](planning/search/SEARCH_PLAN.md)): NR50, IGR, dark ORF, spacers; S4 sequences are not in it |
| `core` | S0–S6 | Adds multi-strain pangenome planes. The S6 search pack lives on the Hetzner search host |
| `planetary` | S0–S7+ | Built and served on Hetzner |

Every edition carries the sequence store for its stages. Sources that cannot be redistributed are replaced by fetch recipes (D-23).

---

## 7. Stage-number mapping for older drafts

`planning/data_layer_and_tests.md` §4.9 was drafted with an earlier numbering. On 2026-10-29 it was relabeled to this file's numbering using the mapping below, which is kept for reading older drafts and audits:

| data-layer draft | this file |
|---|---|
| S0 golden + harness | S0 |
| S1 small genomes | S1–S3 |
| S2 ~1k-genome pilot decisions (DuckLake, K-07, K-12) | S2 gate |
| S3 species reps (L5 decision) | S5 gate |
| S4 MGE collections at scale (drift test) | S4 gate |
| S5 Hetzner port (D-05) | P1 / P2 |
| (no counterpart) | S2C, inserted 2026-10-29 between S2 and S3 without renumbering; S3–S7+ keep their numbers and the rows above are unchanged |

---

## 8. What is not decided yet (owned by named stages)

| Open item | Decided at |
|---|---|
| Sequence-store primary format (Parquet chunks vs bgzip+faidx) | S0 (K-13, by LY-06 row-group counts) |
| Row-group sizes and sort keys | S0 (LY-06) |
| ~~Zombi fit~~ **Decided 2026-10-29:** ZOMBI2 plus custom `synthworld` ([synthworld_spec.md](planning/tests/synthworld_spec.md) §1); fit spike remaining | S−1 (spike) |
| ~~CI host for ring R2~~ **Decided 2026-10-29 (DEC-06):** GitHub Actions on private `ifinkelstein/bio-box` | done |
| Mutation-score target (≥ 80%, *target*) confirmed after the first measurement (TEST_ARCHITECTURE §11 item 3) | S0 |
| Which lab member runs the flexibility check, and on which workflow (TEST_ARCHITECTURE §11 item 4) | S0 (named; TEST_ARCHITECTURE gives no owner) |
| Golden-set accessions marked "pin" (SR1, Lak, Obelisk) | S0 |
| Infernal production-mode recall and cost | S0 |
| geNomad rate on Haswell; whether to rerun or trust shipped calls | S0, then S4 |
| DuckLake vs generated read-only catalog | S2 |
| LIST vs exploded tapes; offset index | S2 |
| Recluster trigger %, ARI floor, drift bands | S1–S4 measurements |
| Whether any L5 index is built | S5 |
| Whether S5's global steps run on the pod or on Hetzner | S3 scaling fit |
| GTDB taxonomy license, PLSDB redistribution terms, INPHARED current count | Before S1–S3 respectively |
| ~~Approve S2C~~ **Approved 2026-10-29** (≤20 per species; S2C-b optional, decided by the first-10 questions). Complete-genome count still to pin | S−1 (count) |
| ~~Z constants~~ **Decided 2026-10-29 (DEC-04):** `Z_prot = domZ = 1e8`, `Z_nt = 1e5` Mb | done |
| Search host for ST4: size and timing of a Hetzner box that holds the S6 search pack (SEARCH_PLAN §9.3) | S5 gate (from the measured S5 pack), before S6 |
