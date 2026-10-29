# Audit Report 6: Software and Data-Platform Architecture

Auditor: Architecture lens (principal engineer, large scientific data platforms: build orchestration, storage engines, release engineering, operations).
Targets: `00-spec.md`, `PLAN.md` (v2.0), `DEPLOYMENT-AND-MALLEABILITY.md`, `PRIOR-ART.md`, `audits/*`.
Date: 2026-10-28.

Verdict: the scientific decomposition is sound. The platform design is not ready to build. The largest input (Logan) is undersized by roughly 8-9x, the build budget assumes rates that were never measured and that contradict the plan's own hardware numbers, the Stage 0 deduplicator has a logic bug and an order-dependent result, the "malleability" layer updates raw rows but not the global statistics the product actually ranks by, and identifiers are not stable across rebuilds. The fix is not more components. It is fewer components, measured rates, a 1% vertical slice, stable content-derived IDs, and a release process.

---

## Step A. Meta-prompt: a deeper version of "audit the plan as an architect"

1. What a top professional would do that the request does not ask for. Recompute every number from first principles and check the documents against each other before judging design. Then ask which user-facing service level each component exists to meet, and delete the ones that serve no measured need. Name the decisions that are hard to undo (identifier schemes, file formats, licenses, sort orders) and settle those first. Treat the build as a release-engineering problem, with determinism, checkpoints, and provenance, not only as a throughput problem.
2. The rawest material. The S3 bucket listing for `s3://logan-pub/c/` (`aws s3 ls --summarize --recursive`), SRA run metadata and SRA STAT taxonomy tables, Hetzner product pages and configurator, NVIDIA datasheets, tool papers and docs (Linclust, pyrodigal, HMMER, Infernal, ProstT5), DuckDB and CRoaring docs, the license text for each data source, and a 1 to 10 GB sample of real Logan contigs to benchmark against. Not the plan's figures or the prior audits' figures, which are the digested material.
3. Methods and tools. Back-of-envelope and roofline arithmetic (FLOPs/s, bytes/s, IOPS), a capacity spreadsheet driven by measured per-unit rates; microbenchmarks (`hyperfine`, `fio`, `iperf3`, `s5cmd`, `nvidia-smi dmon`); a workflow DAG in Snakemake or Nextflow; DuckDB `EXPLAIN ANALYZE`; architecture decision records (ADRs); failure-mode and effects analysis (FMEA); STRIDE threat modeling; an SPDX-style license matrix.
4. What else would show the result is right. Wall-clock and peak-RAM logs from a 1% slice run; bucket sizes reported by AWS; vendor spec pages; bit-identical outputs (equal manifest hashes) from two rebuilds with different thread counts and input orders; an outside engineer rebuilding the slice from the manifest alone; the benchmark suite (B1 to B6) passing on the slice.

Rewritten prompt:
> Audit the `all-the-data` plan as a principal data-platform architect. Work from primary sources (S3 bucket listings, SRA metadata, vendor spec pages, tool papers and docs, source licenses, and a measured sample of real input), not from figures in the plan or prior audits. (a) Recompute every capacity, throughput, size, and latency figure with the arithmetic shown, and list internal contradictions. (b) Draw the real dataflow DAG, and mark causal dependencies, global-statistic barriers, and sources of nondeterminism. (c) For each bespoke component, name the user-facing SLO it serves and propose the minimum viable architecture. (d) Test each malleability claim against the derived global statistics it would leave stale. (e) Do an FMEA/operations and a STRIDE security review. (f) Build a license and redistribution matrix per source and per edition. (g) Propose a phased roadmap with measurable go/no-go gates, starting with microbenchmarks and a 1% vertical slice. For each finding give severity, location, failure scenario, fix, and the external evidence that would confirm the fix.

---

## Step B. Findings

Severity scale: BLOCKER (plan cannot succeed as written), CRITICAL (wrong results or data loss), MAJOR (large cost or schedule error), MODERATE.

### AR-01 BLOCKER. Logan is about 8-9x larger than the plan assumes
- Location: `00-spec.md` Invariant 2; `PLAN.md` §1.3 Stage 0, §2.2 Tier 3, §2.3; `DEPLOYMENT` §4.1.
- Evidence: the Logan paper reports 0.90 Pbp of contigs at 0.31 PB compressed (4.59 Pbp of unitigs at 2.18 PB) ([Chikhi et al., bioRxiv 2024 / PMC12424806](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12424806/)); a 2025 INRAE presentation lists 385 TB at the contigs prefix ([Peterlongo 2025](https://pepi-ibis.inrae.fr/sites/default/files/AG2025/04_Pierre-Peterlongo_DiffusionOK.pdf)). The plan uses 40.8 TB.
- Arithmetic. Transfer of 310-385 TB from AWS us-east-1 to Hetzner:
  - at the default 1 Gbps port (125 MB/s): 385e12 / 1.25e8 = 3.08e6 s = 35.6 days;
  - at the plan's 400 MB/s: 385e12 / 4e8 = 9.6e5 s = 11.1 days;
  - at 10 Gbps line rate: 3.6 days, which a single transatlantic host rarely sustains.
- The plan is also internally inconsistent: 40.8 TB at 400 MB/s is 40.8e12 / 4e8 = 28.3 h, not 68 h; 68 h at 400 MB/s is 97.9 TB.
- Cardinality knock-on. If only 10% of 0.9 Pbp is prokaryotic or viral sequence in contigs of at least 1.2 kb, that is 9e13 bp, or about 9e10 ORFs at about 1 gene/kb. That is about 10x the plan's "~8.5B raw ORFs". Every downstream count (NR90, NR50, dark segments, slots, windows, postings) is unanchored until measured.
- Fix: (1) Replace "40.8 TB" with the measured listing size and record the date. (2) Move Stage 0 to the data: run the sieve in AWS us-east-1 (AWS Batch on spot instances) and ship only filtered prokaryotic and viral contigs or proteins. The Open Data Sponsorship Program covers egress from sponsored buckets ([AWS Open Data FAQ](https://aws.amazon.com/opendata/faqs/)), but compute in your own account is not covered. (3) Before downloading, preselect accessions with SRA run metadata and SRA STAT taxonomy (library source, organism, k-mer taxonomy) so human, mouse, and plant runs are never decompressed. (4) Re-derive all cardinalities from a measured 0.1% random sample of Logan accessions.
- Validation: `aws s3 ls --no-sign-request --summarize --human-readable --recursive s3://logan-pub/c/`; yield of filtered bp per TB from the sample.

### AR-02 BLOCKER. The 4.7-day build budget is not derived from measured rates and omits whole stages
- Location: `PLAN.md` §1.3; `audits/critique_2_systems.md` §6; `audits/SYNTHESIS.md`.
- The total is a sum of point estimates: 68 + 8.8 + 3.5 + 12.5 + 14.3 + 0.5 + 6.0 = 113.6 h. It assumes perfectly serial stages, zero retries, and zero downloads of Tier 1 and Tier 2 sources or reference databases.
- Stage 0 gene calling. pyrodigal speeds up Prodigal by a third to a half through SIMD connection scoring ([pyrodigal docs](https://pyrodigal.readthedocs.io/en/v3.2.1/)), but it is still on the order of 1 Mbp/s per thread (measure this). On 3e14 bp at 96 threads × 1e6 bp/s: 3e14 / 9.6e7 = 3.1e6 s = 36 days. Trying four genetic codes multiplies gene calling by up to 4x unless codes are chosen per contig from a cheap stop-codon-density test first.
- Stage 2a clustering. Linclust clustered 1.6 billion sequences to 50% identity in 10 h on one 28-core server ([Steinegger & Söding 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6026198/)), about 280 core-hours. Scaled linearly to 8.5e9 sequences, that is about 1,490 core-hours for one linclust pass, before the NR50 pass and before the search of 85M queries against 986M targets for 3Di transfer. The plan budgets 105 CPU-hours for all three. With the AR-01 cardinalities it is 10x more again.
- 3Di transfer source. Foldseek distributes clustered databases (for example `mgnify_esm30`), not 3Di strings for all 772M ESMAtlas entries. Transferring 3Di from the full atlas means downloading the structure tarballs and running `foldseek createdb` yourself. That download and conversion step is absent from the budget.
- Stage 2c ProstT5. The plan lists the RTX 6000 Ada at "~91 TFLOPS FP16/BF16, ~182 TOPS INT8". 91.1 TFLOPS is the FP32 non-tensor peak; the tensor figure is 1,457 TFLOPS FP8 with sparsity, and dense INT8 is 728 TOPS ([Lenovo Press LP1940](https://lenovopress.lenovo.com/lp1940.pdf); [NVIDIA forum](https://forums.developer.nvidia.com/t/l40-vs-rtx-6000-ada-fp16-fp8-throughput/245204)). The ProstT5 encoder (~1.2B parameters) costs about 2 × 1.2e9 = 2.4 GFLOP per residue before attention. The claimed 62,500 aa/s therefore needs 1.5e14 FLOP/s = 150 TFLOPS sustained. By the plan's own 91 TFLOPS figure that is 165% of peak, which is impossible. Against the roughly 364 TFLOPS dense FP16 tensor peak (inferred, not on the datasheet) it is 41% model FLOPs utilization, achievable only with TensorRT or FlashAttention kernels. At a typical 15-25% utilization: 23-38k aa/s, so 3.24e9 aa takes 24-39 GPU-hours, not 14.3.
- Unbudgeted work: gene-tree parsimony `H_{B|T_A}` needs a protein tree per family A (FastTree on up to 1e5-1e6 families); vehicle single-linkage clustering (±15-ORF Jaccard over all loci of every pair, quadratic for hub pairs); circular block-shuffle decoys for `q_emp` (at least one full re-run of the pair-summary pipeline on shuffled tapes); `subfamily_50_id` clustering of excised envelopes; Stage 1.5 split-gene reconstruction; IGR extraction; gLM-Domain training, which is folded into "Stage 3: 6.0 h". Training a 12-layer transformer over billions of slots on one Ada GPU in under 6 h is not credible.
- Fix: replace the §1.3 table with a capacity model: rows are stages, columns are measured per-unit rate (source and date), measured input cardinality, derived CPU-hours, GPU-hours, peak RAM, peak scratch, and a 2x contingency. Stages with no measurement are labeled "unmeasured" rather than given a number. gLM-Domain moves to its own project and budget.
- Validation: microbenchmark logs from M0 (roadmap below); slice wall-clock within 2x of the model.

### AR-03 CRITICAL. The Stage 0 deduplicator has a logic bug, a sizing error, and a mis-stated ANI test
- Location: `PLAN.md` §2.3 Filter 4; `00-spec.md` Invariant 2; `critique_2` §4.2.
- Logic bug. On first sighting the record goes only to `tier3_unique_contigs.bin.zst` ("Zero RocksDB/Bitmap allocation"). On second sighting the code looks up `(stage0_arch_hash, biome_bucket)` in RocksDB, where nothing was ever written, so there is no representative sketch to compare against. Either first sightings are also written to RocksDB (adding about 1e9 × 16 B = 16 GB on NVMe plus write amplification), or the design does not work.
- Sizing error. 1e9 16-bit fingerprints occupy 1e9 × 2 B = 2 GB, not 16 GB. 16 GB holds 8e9 fingerprints. The 26 GB RAM total is built on this.
- Cuckoo false positives. With 4-way buckets and 16-bit fingerprints, the false-positive rate is about 2 × 4 / 65,536 = 1.2e-4. Over 1e9 inserts and more lookups that is about 1.2e5 or more false "seen" verdicts, each sending a genuinely novel contig down the second-sighting path. The plan does not say what happens when that path finds no representative.
- ANI test. A 64-bit one-permutation b-bit MinHash with b = 1 has a collision probability of J + (1 − J)/2. For 95% ANI at k = 15, Mash gives 2J/(1 + J) = e^(−15 × 0.05) = 0.472, so J ≈ 0.31 and the expected mismatch count is 64 × (1 − 0.31)/2 ≈ 22. The rule `POPCNT(xor) <= 6` therefore requires (1 − J)/2 ≤ 0.094, so J ≥ 0.81. That is d ≈ −ln(0.895)/15 ≈ 0.0074, or about 99.3% ANI. The binomial standard deviation of about 3.8 bits makes the call noisy, and Jaccard is the wrong measure for contigs of unequal length, which needs containment. The plan's "≥95% ANI in 1 CPU cycle" is wrong by about 4 ANI points, and it has a high false-negative rate.
- Non-determinism. "First sighting becomes the representative" depends on S3 listing order, thread scheduling, and retries. Two builds of the same inputs produce different representatives, different `dense_contig_id` orderings, and different downstream counts.
- Fix: replace streaming first-seen deduplication with a two-pass, sort-based design. Pass 1 (embarrassingly parallel, idempotent per accession) emits `(arch_hash, contig_key, sketch, length, quality)` records to sharded Parquet. Pass 2 runs an external sort by `arch_hash`, groups, and picks the representative by a deterministic rule (longest, then highest quality, then lexicographic `contig_key`), with ANI or containment verified by a real estimator (a 128-1024-slot sketch, or skani or Mash on candidate groups only). Sorting 1e10 records of 32 B (320 GB) on NVMe is a few hours. No Cuckoo filter or RocksDB is needed.
- Validation: two builds with shuffled input order and different thread counts produce equal manifest hashes; ANI calibration curve of the estimator against skani on 10k contig pairs.

### AR-04 CRITICAL. Identifiers are not stable across rebuilds, so sidecars, remaps, and external references break silently
- Location: `DEPLOYMENT` §3.2, §3.3, §4.2 ("zero-downtime quarterly release"); `PLAN.md` §4.1, §4.2.
- `dense_contig_id` is assigned by sort order `(tier, tax_order_id, stage2_arch_hash)`; `subfamily_50_id`, `DARK_MMSEQS30` IDs, `super_cluster_id`, `vehicle_cluster_id`, and `window_id` are arbitrary integers from clustering runs. A quarterly rebuild or `bio-box compact` renumbers all of them. A sidecar keyed by `subfamily_50_id` or `(dense_contig_id, slot_idx)` then joins against the wrong entities with no error.
- Fix:
  1. Two-layer IDs. Each entity gets a stable, content-derived public ID (`contig_uid` = hash of the source accession and contig name; `protein_uid` = hash of the amino-acid sequence, in the style of UniParc UPI or SEGUID; `cluster_uid` = `protein_uid` of the cluster representative plus the clustering parameter set) and a release-local dense integer for physical layout. Every table carries the dense ID; a crosswalk table maps dense IDs to UIDs per release.
  2. Release IDs. Every file, sidecar, and API response records `release_id` (for example `atd-2027.1`) and a parameter-set hash.
  3. Release-to-release crosswalks for clusters (`old_cluster_uid → new_cluster_uid, overlap_fraction`), as UniRef and GTDB do.
  4. `register-sidecar` requires the sidecar's key namespace and `release_id`, and refuses to attach across releases without an explicit crosswalk join.
- Validation: a sidecar registered on release N and migrated to N+1 yields the same rows for unchanged entities in a regression test.

### AR-05 CRITICAL. Malleability updates rows, not the global statistics the product ranks by
- Location: `00-spec.md` Invariant 12; `DEPLOYMENT` §3.1-3.4; `PLAN.md` §6.2.
- What `append` can update correctly: raw tapes, postings, window bitmaps, and additive counts (`n_total_occurrences`, OR of phylum and order bitmasks).
- What it leaves stale or wrong:
  - NR90 and NR50 membership. New proteins are not clustered against base representatives, and a new protein may change which sequence should represent a cluster. (MMseqs2 `clusterupdate` exists but is not in the plan.)
  - `subfamily_50_id`, `catalytic_motif_hash` reference columns, and the dark-family dictionary.
  - Deduplication. The delta is not deduplicated against the base, so occurrence counts inflate.
  - Non-decomposable pair statistics: `n_veh_indep_events` (single-linkage clusters span epochs), `ctx_entropy_*`, `slot_fold_entropy_ab_bits`, `micro_geometry_score` (IQR is not mergeable from IQRs), `func_coupling_llr` (depends on the regime-conditioned null, which depends on global token frequencies), and `empirical_fdr_qval` (decoys must be regenerated).
  - Index 3: `core_chr_bits` (the ≥50% threshold changes with new strains), `mge_accessory_bits` (<30% threshold), `hq_closed_or_multi_strain_mask`, and `transition_edges` (Fitch reconstruction is global).
  - gLM-Domain vocabulary, embeddings, `slot_substitutions.parquet`, and the USearch index.
- Token remap. Merges can be applied to additive counts. Entropies, LLRs, q-values, and bitvector planes cannot. Splits cannot be applied at all from an overlay keyed by `(old_subfamily_50_id, catalytic_motif_hash)`, because the split criterion (a new alignment) lives at the segment level and needs a segment-level overlay or a rewrite of the tapes.
- Roaring union. Concatenating containers is valid only when epoch `window_id` ranges begin on a multiple of 65,536 windows, so their high 16-bit keys do not overlap. Otherwise containers must be merged. The plan does not require this alignment.
- Fix: state the contract in the spec. Deltas are "raw-searchable immediately; statistics as of release N". Every statistic column carries `stats_release_id`, and the API flags results whose statistics are stale relative to the data. Recomputation is a scheduled release, not an online operation. Private user epochs are search-only overlays that never feed published statistics. Align epoch window ranges to 65,536.
- Validation: append 1% new genomes, then compare query results against a full rebuild; the differences must appear only in columns flagged stale.

### AR-06 MAJOR. Too many bespoke components for the stated needs; latency SLOs are vanity targets
- Location: `00-spec.md` §1 and Invariant 6; `PLAN.md` §1.2, §4.1; `DEPLOYMENT` §2.
- The serving path combines a custom Rust daemon (`atd-engine`), CRoaring, RocksDB, DuckDB, USearch, Arrow Flight, HTTP, a DuckDB UI, an MCP server, a Python CLI, and a custom binary `PackedSlotDescriptor` format. Each is a format and an operational surface to maintain. No user story requires a 4 ms response: the consumers are scientists and agents running analyses where 200 ms to 2 s is acceptable, and the system of record is batch-computed pair statistics.
- Index 1 has defects of its own:
  - Window boundaries. Two tokens in adjacent windows (slots 15 and 16) are 1 slot apart but share no `window_id`, so the "AND proves proximity" claim yields false negatives. Correct queries need `win(A) AND (win(B) OR shift(win(B), ±1))`.
  - `orf_mod_8` aliases distances d and d + 8, and the 4,096-contig block layout in `PLAN.md` differs from the 65,536-window layout in `critique_2`. The two designs were never reconciled.
  - Headroom. Contigs start on window boundaries, so windows ≥ contigs ≈ 2.5e9 against a 2^32 = 4.29e9 limit. Roughly one more Logan-scale append overflows 32-bit Roaring.
- `window_postings` latency. A hub anchor such as PF00005, with about 5e7 occurrences × 12 partners = 6e8 rows × 7.5 B, is 4.5 GB, not "1 to 4 row groups (5-80 MB)". The 3-10 ms guarantee holds only for rare anchors.
- Minimum viable architecture: Parquet files plus DuckDB, with a thin FastAPI (or Arrow Flight) service and MCP tools wrapping a fixed set of parameterized queries. `window_postings` sorted by `(anchor, partner)` answers anchor scans; `contig_tapes` answers neighborhood fetches via an explicit sorted `dense_contig_id` column (DuckDB prunes on row-group min/max statistics; a virtual `file_row_number` is not a stored column and should not be relied on for pruning, per [DuckDB Parquet docs](https://duckdb.org/docs/current/data/parquet/overview)). Add Roaring windows or a custom daemon only when a measured SLO fails on the `core` edition.
- SLOs to adopt: anchor scan p95 < 2 s for non-hub anchors and < 10 s for hubs (or served from precomputed `directon_pair_summary`); neighborhood fetch of 1,000 contigs p95 < 1 s; pair-summary lookup p95 < 100 ms.

### AR-07 MAJOR. No workflow manager, checkpointing, or provenance
- Location: `PLAN.md` §1.3, §2-3 (pipelines described as shell pipes and ASCII diagrams).
- A multi-week build on one RAID-0 node with no DAG, retries, or per-stage manifests will be restarted from scratch after the first failure.
- Fix: a Snakemake or Nextflow DAG with one rule per stage and shard. Inputs are pinned by content hash (source accession lists with dates, database versions such as Pfam 37.x and Rfam 15.0, tool container digests). Outputs are content-addressed Parquet shards with a per-shard manifest of inputs, parameters, tool versions, wall-clock, and peak RSS. Shards are idempotent so retries are safe. Use OpenLineage or plain JSON lineage. The same DAG runs at 0.1%, 1%, and full scale with only a config change.
- Validation: kill -9 at random points during the slice build; resume completes with identical manifests.

### AR-08 MAJOR. Hardware and network facts are wrong or inconsistent
- Location: `00-spec.md` §1; `PLAN.md` §1.1; `DEPLOYMENT` §4.
- CPU. The AX162 ships the EPYC 9454P (48 cores / 96 threads) ([Hetzner AX162 announcement](https://www.hetzner.com/de/news/new-ax162/)). The spec says "EPYC 9654" (96 cores), and `PLAN.md` says "48-96 cores / 96 threads". Pick one; all budgets should assume 48 physical cores.
- Disks. The AX162-R base is 2 × 1.92 TB NVMe with 256 GB RAM; larger drives are configurable add-ons ([Hetzner AX162 product page](https://www.hetzner.com/de/dedicated-rootserver/ax162/)). Treat 2 × 7.68 TB as an option to price and confirm in the configurator.
- GPU node. The GEX130 base has 128 GB RAM and 2 × 1.92 TB NVMe, not 256 GB ([Hetzner GEX130 announcement](https://www.hetzner.com/de/news/gpu-server-gex130/)). The successor GEX131 has an RTX PRO 6000 Blackwell with 96 GB and 256 GB RAM ([Hetzner GEX131](https://www.hetzner.com/de/news/gex131/)), which would also change the ProstT5 estimate.
- Network. Hetzner dedicated servers default to a 1 Gbps port. The 10G uplink is a €43/month add-on with 20 TB of outgoing traffic included, then €1/TB; incoming traffic is free ([Hetzner 10G uplink](https://docs.hetzner.com/robot/dedicated-server/network/10g-uplink/)). The "52 minutes at 10 Gbps" and "1.5 minutes for lite on a MacBook" claims assume 10 Gbps at both ends. At 1 Gbps: planetary 3.45e12 / 1.25e8 = 7.7 h; core 56 min; lite 11 min. At 200 Mbps home broadband, lite takes 57 min.
- Distribution host. Each planetary download is 3.45 TB of egress, so the 20 TB allowance covers about 5 downloads per month. Host release packs on Cloudflare R2 (no egress fees) or a sponsored bucket, not on the query node.
- Fix: one "Hardware reference" table in `PLAN.md` §1.1 with vendor links and verification dates; every budget references it.

### AR-09 MAJOR. Internal size and count contradictions
- Location: `PLAN.md` §1.2, §2.2, §4.1, §4.2; `00-spec.md` §1.
- Contigs against slots. The sieve requires at least 2 features per contig (with a small whitelist), so 2.5e9 contigs imply at least 5e9 slots, not 3.6e9. `critique_2` used 5e8 contigs. One of these is wrong by 5x.
- `window_postings` size. 1.8e10 rows × 7.5 B/row = 135 GB, not 1.35 TB. Alternatively, 3.6e9 slots × 12 partners within ±6 = 4.3e10 directed pairs × 7.5 B = 324 GB. Neither gives 1.35 TB.
- Index 1 postings. `token_to_windows` is said to hold `seq_family_id`, `subfamily_50_id`, and `struct_fold_id`; `critique_2` computed 8.3 GB for two token kinds. With three to five kinds, 3.2 GB is low by 3-6x.
- Tier 1 input size. 1.8e6 assemblies × ~4e6 bp = 7.2e12 bp; gzip FASTA at about 0.28 B/bp gives about 2 TB, not 650 GB. (AllTheBacteria's phylogenetic compression is much smaller, but it needs its own decompression stage.)
- The spec says "Ten Non-Negotiable Invariants" and lists 12.
- Fix: a single generated `sizing.yaml` (cardinalities, bytes/row, compression ratios) from which every size in every document is derived; CI fails if the documents disagree with it.

### AR-10 MAJOR. RAID-0 as the only copy during a multi-week build
- Location: `DEPLOYMENT` §4.2; `PLAN.md` §1.1.
- One NVMe failure loses all intermediate state. For a serving node that is fine (re-pull the pack). For the build it is not.
- Fix: object storage (Hetzner Object Storage, R2, or S3) is the system of record for every completed shard; RAID-0 is scratch and cache only. The Storage Box holds cold checkpoints only (it is not a high-throughput target). Monitoring: Prometheus node_exporter, NVMe SMART and wear metrics, per-rule metrics from the workflow manager, disk-full alerts. A quarterly restore drill from object storage.

### AR-11 MAJOR. The DuckDB catalog cannot be shared the way the CLI implies
- Location: `DEPLOYMENT` §2.3, §3.2, §4.2 ("attaches `catalog.duckdb` in read-write overlay mode").
- DuckDB allows one read-write process or many read-only processes on a database file, not both. With `serve` running, a separate `docker run ... append` or `register-sidecar` container cannot take the write lock, and `serve` cannot see the new views without reopening.
- Fix: the catalog is a declarative manifest (`catalog/release.json` plus `sidecars/*.json`) from which each server process builds an in-memory DuckDB with views at startup or on a reload signal. Mutations go through the daemon's admin API, which writes the manifest atomically. Alternatively adopt a lakehouse table format with snapshots (DuckLake or Apache Iceberg), which provides epochs, time travel, and atomic commits without inventing them.
- Sidecar cost. A hash join against a 2 GB sidecar builds a hash table per query unless it is persisted. "< 1 ms" holds only for small sidecars or pre-joined materializations.

### AR-12 MAJOR. Security: the default deployment exposes arbitrary SQL
- Location: `DEPLOYMENT` §2.3 (`-p 8420:8420`, "HTTP/Arrow-Flight/DuckDB-UI/MCP").
- Unauthenticated DuckDB SQL on a public port allows reading local files (`read_csv('/etc/...')`), writing with `COPY ... TO`, installing extensions, and exfiltrating data through `httpfs`. Combined with the docker socket or a writable `/data`, it is a host compromise.
- Fix: bind to 127.0.0.1 by default and put a TLS reverse proxy with token authentication in front. Expose parameterized query endpoints, not raw SQL. Where SQL is offered, use a per-request read-only connection with `enable_external_access = false`, extension autoload disabled, `lock_configuration = true`, and memory, thread, and timeout limits. Admin operations (`append`, `register-sidecar`, `remap-token`) go on a separate authenticated port. MCP tools get scoped tokens. Run the container as non-root with `/data` mounted read-only for `serve`. Keep a STRIDE table in `DEPLOYMENT.md`.

### AR-13 MAJOR. Licensing and redistribution of derived data packs is unaddressed
- Location: `DEPLOYMENT` §4.1 (editions), `PLAN.md` §2.2.
- Examples that need a decision:
  - GTDB taxonomy is distributed under CC BY-SA 4.0 (confirm on gtdb.ecogenomic.org). Share-alike can propagate to packs that embed it.
  - JGI/IMG data (IMG/VR, IMG/PR, IMG/M) fall under the JGI data policy, where use-restricted legacy data require PI permission to redistribute ([JGI data policy](https://jgi.doe.gov/data-policy-support/legacy-data-policy)).
  - PaperBLAST's `litsearch.db` contains text snippets from papers. Redistributing derived snippets needs review.
  - CNCB/CNGBdb data access terms and Chinese data-export rules need review before redistribution.
  - Each HMM library (Pfam CC0, NCBIfam, CDD, PHROGs, PADLOC, DefenseFinder, TXSScan, MIBiG CC BY 4.0, ECOD) has its own terms.
  - Tools: MMseqs2, Foldseek, and HH-suite are GPLv3, which is compatible with redistribution in a container if source is offered.
- Fix: `LICENSES.md` with a matrix of source, license, attribution text, redistribution allowed (yes, no, or conditional), and which editions include it. The build emits a per-edition license manifest. Sources that cannot be redistributed ship as "bring your own data" recipes that rebuild locally.

### AR-14 MODERATE. Container design couples model weights to the image
- Location: `DEPLOYMENT` §2.1; `00-spec.md` Invariant 11.
- ProstT5 weights determine the 3Di strings in the database, and gLM-Domain weights are tied to the release vocabulary. Both belong in the data pack, versioned with the release, not in the compute image. A CUDA plus PyTorch runtime alone exceeds about 3 GB, so "4.2 GB" covering CUDA, all tools, and 1.7 GB of weights is optimistic.
- Fix: two images. `bio-box-cpu` (query, annotate, append) and `bio-box-build-cuda` (build only). Both are pinned by digest in the release manifest. Weights live in `/data/models/<release_id>/`.

### AR-15 MODERATE. No query language or versioned API is specified
- Location: `PLAN.md` §1.2 ("Gapped Operon Pattern Search"); `DEPLOYMENT` §2.3.
- The core product, pattern queries such as `[DomA, DomB, ncRNA_C] within ±k`, has no grammar, no semantics for strand, gaps, token level (family, subfamily, fold), or tier filters, and no result schema.
- Fix: a small, documented pattern DSL (in the spirit of cblaster or PROSITE-style operon regular expressions) compiled to SQL over `window_postings` and `contig_tapes`; an OpenAPI spec at `/v1`; every response carries `release_id`, `stats_release_id`, the query hash, and pagination. MCP tools are thin wrappers over `/v1`.

### AR-16 MODERATE. Causal ordering is unstated
- Location: `PLAN.md` §3.2 Stage 2c ("occurring in ≥2 non-identical contigs") and §2.3.
- Choosing "recurrent" dark segments needs occurrence counts from the deduplicated tapes (Stage 3 inputs), while Stage 2c runs before Stage 3. Regime labels need hallmark tokens from Stage 2, and vehicle clustering needs regime labels. These are global barriers.
- Fix: draw the DAG explicitly in `PLAN.md` with barrier stages marked (dedup grouping, NR clustering, occurrence counting, regime labeling, vehicle clustering, pair aggregation, decoys, tree reconstruction). Each barrier is a checkpoint and a release-recompute boundary for AR-05.

---

## Proposed phased roadmap with go/no-go gates

| Milestone | Scope | Duration (est.) | Go/no-go gate |
|---|---|---|---|
| M0 Measure | Bucket listings; SRA metadata accession preselection; microbenchmarks of pyrodigal (per code), linclust, hmmsearch, cmsearch with nhmmer prefilter, ProstT5 aa/s, Foldseek, DuckDB scans; `iperf3` and `s5cmd` from us-east-1 to Hetzner; ADRs for IDs, formats, licenses | 2-3 weeks | Capacity model rebuilt from measured rates; total build within budget at 2x contingency, or scope cut (for example, defer Logan) |
| M1 Vertical slice (0.1-1%) | 3-5 orders spanning the benchmark loci (for example Enterobacterales, Bacillales, Bacteroidales, Methanobacteriales) plus their phages and plasmids plus 2 biomes of metagenomes; full DAG in Snakemake or Nextflow; Parquet + DuckDB only; stable UIDs; license manifest | 6-8 weeks | Two rebuilds bit-identical; kill-and-resume identical; B1, B4, B5, B6 on slice systems pass; sizes within 2x of `sizing.yaml` |
| M2 `core` (Tier 1 + Tier 2) | Full HQ genomes and MGEs; release process (`release_id`, crosswalks, stats versioning); `/v1` API with auth | 2-3 months | Query SLOs met with DuckDB alone, or a measured gap justifies Index 1; B1-B3 pass at full Tier 1 |
| M3 Tier 3 without Logan | SPIRE, AllTheBacteria, SAGs, long-read metagenomes | 1-2 months | Rarefaction: novel NR50 clusters per added Tbp still above a preset threshold |
| M4 Logan via compute-to-data | Sieve on AWS Batch spot next to the bucket; ship filtered shards | cost-gated | Cost per novel NR50 cluster within budget; yield curve from the 0.1% sample holds |
| M5 Incremental updates and gLM-Domain | Search-only deltas with stale-stats flags; gLM as a separate project with its own budget | ongoing | Delta-vs-rebuild regression test; gLM benchmarks B2A/B2B |

---

## Sources
- [Logan paper (PMC12424806)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12424806/); [bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2024.07.30.605881.full.pdf); [Peterlongo 2025 slides](https://pepi-ibis.inrae.fr/sites/default/files/AG2025/04_Pierre-Peterlongo_DiffusionOK.pdf); [AWS registry: Logan](https://registry.opendata.aws/pasteur-logan); [IndexThePlanet/Logan](https://github.com/IndexThePlanet/Logan)
- [AWS Open Data FAQ](https://aws.amazon.com/opendata/faqs/)
- [Hetzner AX162 announcement](https://www.hetzner.com/de/news/new-ax162/); [AX162 product page](https://www.hetzner.com/de/dedicated-rootserver/ax162/); [GEX130](https://www.hetzner.com/de/news/gpu-server-gex130/); [GEX131](https://www.hetzner.com/de/news/gex131/); [10G uplink](https://docs.hetzner.com/robot/dedicated-server/network/10g-uplink/); [add-on prices](https://docs.hetzner.com/robot/general/pricing/price-list-for-additional-products)
- [Lenovo Press RTX 6000 Ada](https://lenovopress.lenovo.com/lp1940.pdf); [NVIDIA developer forum](https://forums.developer.nvidia.com/t/l40-vs-rtx-6000-ada-fp16-fp8-throughput/245204)
- [Linclust, Steinegger & Söding 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6026198/)
- [pyrodigal docs](https://pyrodigal.readthedocs.io/en/v3.2.1/); [pyrodigal benchmarks](https://pyrodigal.readthedocs.io/en/stable/_sources/benchmarks.rst.txt)
- [DuckDB Parquet overview](https://duckdb.org/docs/current/data/parquet/overview); [DuckDB Parquet tips](https://duckdb.org/docs/current/data/parquet/tips.html)
- [BFVD / AFDB counts, NAR 2025](https://tagteam.harvard.edu/hub_feeds/4353/feed_items/14028250); [ESM Atlas README](https://huggingface.co/spaces/wenkai/FAPM_demo/blob/ac229dd601851c1e5ad4eaa875c2c8580e94df52/esm_scripts/atlas/README.md)
- [JGI legacy data policy](https://jgi.doe.gov/data-policy-support/legacy-data-policy)

Unverified items to confirm in M0: pyrodigal per-thread Mbp/s; dense FP16 tensor peak of the RTX 6000 Ada; GTDB license terms; CNCB redistribution terms; HMMER MSV cells/s on Zen 4.
