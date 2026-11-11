# Design stress test, 2026-11-09

Purpose: test the four decisions that cannot be changed cheaply after S1 (slot key, release-local
tokens, materialized `slot_pairs`, digest-sorted bgzip shards) from built files, not from the plan
documents. Verdicts and the resulting tasks are in [`planning/tasks/TASKS.md`](../tasks/TASKS.md).
Scripts are in `scripts/stress/`. Code under test: commit `2556ddd`.

Hosts. Laptop: Apple Silicon (arm64), macOS, Python 3.13.12, packages from `uv.lock`, no MMseqs2.
Pod: wilkcomp02, x86-64, test image `bio-box-test_74b3e05.sif` (Python 3.14.10, MMseqs2 18.8cc5c)
with the checkout's HEAD `src/` bound over the image's source. The HEAD image
(`bio-box-test@sha256:48296e3a…c6944c`) could not be pulled: the registry needs a `read:packages`
credential (task 0-02). Every pod command ran as
`systemd-run --user --scope -p MemoryMax=200G nice -n 19 singularity exec --cleanenv …`.

Labels: *measured*, *derived*, *estimate* as in `AGENTS.md` §6.

## 1. Cross-host build (measured)

`scripts/stress/build_s0.py`, stage A (ingest, gene calling, IGRs, tapes) on the golden set:
36,475,789 bp, 21 contigs, 15 genomes, 33,716 genes, 33,442 distinct proteins, 3 origin-wrapping genes.

| | laptop | pod |
|---|---|---|
| manifest hash | `1a3b9f7ffe28ce9cf36aec13da205e3f` | `1a3b9f7ffe28ce9cf36aec13da205e3f` |
| files byte-identical | 15 of 15 | |
| ingest / gene calling, 1 thread | 6.2 s / 13.4 s | 9.1 s / 31.9 s |

The golden set is 36.5 Mbp, not the 43 Mbp recorded on 2026-10-30, so the pod's gene-calling rate
is 1.1 Mbp/s per thread, not 1.2.

Stage B on the pod (clustering at 3 threads, no domain scans): 33,010 NR90 and 29,219 NR50 clusters,
1,347,888 pair rows, clustering 87 s, pairs 13 s, peak RSS 1.27 GB.

Full build on the pod (Pfam 38.2 and NCBIfam 20.0 at 16 threads, clustering at 8 threads, tokens,
pairs), manifest `2548b144db0d7ec9a1cd8074e29a7e8a`, peak RSS 2.6 GB:

| Step | Wall | Output |
|---|---|---|
| Pfam scan | 290 s | 56,786 domains on 27,176 proteins (identical to the 2026-10-30 record) |
| NCBIfam scan | 599 s | 18,361 domains on 12,548 proteins (identical to the 2026-10-30 record) |
| clustering | 48 s | 33,010 NR90 and 29,219 NR50 clusters |
| tokens | 1.5 s | 11,698 family, 29,219 NR50, 33,010 NR90 tokens |
| pairs | 20 s | 1,823,702 rows |

Determinism: `cluster_member` and `cluster_member_aln` are byte-identical between the 3-thread and
the 8-thread build, as are all 15 stage-A files.

With all three token levels `slot_pairs` holds 54.1 rows per gene at 7.80 bytes per row: 422 bytes
per gene, more than the sequence shards (334 bytes per gene). `l1_domains` is 53.8 bytes per row at
2.25 domains per protein.

## 2. Independent tools (measured)

Laptop-built shard copied to the pod and read by tools the project did not write.

| Check | Result |
|---|---|
| `bgzip -t` on the shard | passes |
| `samtools faidx` (1.19.2) on every non-wrapping gene, both strands, using the repo-written indexes | 33,713 of 33,713 slices byte-identical |
| `.fai` and `.gzi` rebuilt by samtools | byte-identical to the repo's |
| Prodigal 2.6.3 binary, single mode, ten prokaryotes | 32,195 of 32,198 bio-box genes exact; the differences are at the ends of three circular contigs, where bio-box calls origin-wrapping genes and Prodigal calls edge fragments |
| RefSeq CDS, *E. coli* K-12 (GCF_000005845.2, non-pseudo) | 4,119 of 4,303 matched by stop (95.7%); 3,889 of those with the same start (94.4%); 199 bio-box genes without a RefSeq stop |

## 3. Storage per unit and projection

Measured on the golden set (bytes on disk, zstd Parquet or bgzip):

| Table | Unit | Bytes per unit |
|---|---|---|
| sequence shards (with indexes) | bp | 0.309 |
| `protein_seq` | distinct protein | 227 (0.73 per aa) |
| `protein_occurrence` | gene | 32.2 |
| `slot_protein` | gene | 19.5 |
| `contig_tapes` | slot | 30.4 |
| `igr` | row (0.835 rows per gene) | 12.7 |
| `cluster_member` | protein and level | 36.5 |
| `cluster_member_aln` | protein and level | 26.2 |
| `tape_tokens`, two levels | gene | 9.0 |
| `slot_pairs`, two levels | row; gene | 7.7; 308 (40 rows per gene) |
| `slot_pairs`, three levels | row; gene | 7.8; 422 (54.1 rows per gene) |
| `l1_domains`, Pfam and NCBIfam | row; gene | 53.8; 120 |

At S0 almost every cluster is a singleton and there are 21 contigs, so `contig_uid` compresses to
nothing in `slot_pairs`. With millions of contigs the 16-byte key no longer dictionary-encodes, so
bytes per row rise toward the plan's own 15–25 (*derived*).

Unit counts from GTDB R226 metadata (`bac120` and `ar53`, *measured*, `scripts/stress/gtdb_counts.py`):

| Stage | Genomes | Gbp | Contigs | Proteins |
|---|---|---|---|---|
| S2C, complete genomes ≤20 per species | 25,458 | 106.4 | 56,387 | 97.7 M |
| S3, species representatives under 2 Mbp | 38,750 | 52.4 | 4.17 M | 53.2 M |
| S5, all species representatives | 143,614 | 465.5 | 24.95 M | 437.2 M |
| S6, all genomes | 732,475 | 2,587.6 | 111.5 M | 2,436 M |

The plan gives S3 as 9–12k genomes and 10–13 Gbp; the metadata gives 38,750 and 52.4 Gbp.
The median genome has 93 contigs (99th percentile 868).

Projection (*derived*: measured bytes per unit × metadata counts; pair rows per gene per level are
between a lower bound that spreads each genome's genes evenly over its contigs and 20):

| Table | S2C | S5 | S6 | Budget in `sources_and_sizing.md` |
|---|---|---|---|---|
| sequence shards | 33 GB | 144 GB | 800 GB | S5 store 207 GB including proteins |
| `protein_seq` (no deduplication assumed) | 22 GB | 99 GB | 553 GB | 0.55 bytes/aa; measured 0.73 |
| `cluster_member`, two levels | 7 GB | 32 GB | 178 GB | S5 "~8 GB membership" |
| tapes | 3 GB | 13 GB | 74 GB | 40 bytes per slot |
| pair rows per level | 1.95 × 10⁹ | 6.2–8.7 × 10⁹ | 37.5–48.7 × 10⁹ | |
| `slot_pairs`, 2.66 levels, at 7.7 bytes per row | 40 GB | 128–178 GB | 0.77–1.0 TB | S5: tapes, pairs and statistics together 25–60 GB |
| `slot_pairs` at 20 bytes per row | 104 GB | 330–465 GB | 2.0–2.6 TB | |
| `token_postings` (proposed), about 2.8 rows per gene at 20 bytes | 5 GB | 24 GB | 136 GB | *estimate*; measure at B-18 |

## 4. Slot key (measured)

`scripts/stress/step4_keys.py` on the *E. coli* K-12 chromosome (4,318 slots), using the repo's own
IGR, tape and pair functions on the changed slot list. "Rewritten" counts rows whose key or content
changed.

| Change | Key | Slot rows rewritten | IGR rows rewritten | Pair rows rewritten (of 86,360, one level) |
|---|---|---|---|---|
| add one ncRNA at slot 3,104 | ordinal | 1,214 | 1,006 | 24,342 |
| | coordinate | 0 | 1 | 110 |
| split the longest gene (slot 1,946) into two segments | ordinal | 2,372 | 1,946 | 47,430 |
| | coordinate | 1 | 2 | 130 |

Under the ordinal key, 86,340 of 86,360 pair primary keys survive the change while naming different
genes. Tables that hold an ordinal as key or reference: `slot_protein`, `protein_occurrence`, `igr`,
`contig_tapes`, `tape_tokens`, `slot_pairs`, `neighborhoods`, `search_region`; 12 source modules.

## 5. Layout at 1,000× the fixture (measured, laptop)

`scripts/stress/synth_scale.py`: 3,000 contigs, 797 Mbp, 854,000 slots, 12 shards of 2²⁶ bases,
44.7 M pair rows in one file of 364 row groups. Sequences went through the repo's `ingest`; gene
tables and tokens are replicated and synthetic (their bytes per row are not representative).

| Query | Read |
|---|---|
| Q1 count, family pair on 1,000 loci | 1 of 364 row groups, 0.38 MB |
| Q1 count, NR90 pair on 10 loci | 1 of 364 |
| Q4 instances of a pair | 1 of 364 |
| Q3 partners of a hub token (26,000 loci) | 5 of 364 |
| 1,000 neighborhoods of ±10 genes, sequence (16.7 Mbp) | 1,146 bgzf blocks, 22.8 MB; `contig_index` scanned whole 3 times and `seq_locator` twice |
| tape and tokens for one neighborhood | the only row group: all 3,000 contigs, 2.9 MB |
| export of 500 whole contigs (135 Mbp) | 2,463 blocks, 48.7 MB (0.36 bytes/bp), all 12 shards |

Limits of the current code:

| Item | Measured | Derived consequence |
|---|---|---|
| `ingest` | 153 s for 797 Mbp (5.2 Mbp/s, one thread); holds every new sequence of the stage in memory | S5 does not fit in memory |
| pair builder | 544 bytes of Python heap per row | S1 (about 3 M genes, 120 M rows at two levels) needs about 65 GB |
| `SeqStore` open | 2.0 s and 414 MB of heap per 10⁶ contigs | S5 about 50 s and 10 GB per connection; S6 about 220 s and 46 GB |
| sequence macros | a literal key reads 1 of 9 row groups of a 10⁶-row index; a joined key reads 9 of 9 | every batched SQL sequence fetch scans the whole index |
| tape row group | 8,192 contigs | at S2C one row group holds about 8,192 chromosomes: 33 M slots, about 1 GB to decode for one neighborhood |

## 6. Mutation testing

Incomplete. mutmut generated 2,937 mutants for `l0seq`, `l1` and `l3` (laptop; tests `tests/r1` and
`tests/r2`; the MMseqs2, samtools and `hmmsearch` tests skip there). It left 594 alive, but its
per-mutant test selection misses real kills (confirmed on one `genecall` mutant), so 594 is an upper
bound. A re-check of the 594 against the full test selection was interrupted on 2026-11-09 and has
no result. The raw survivor lists, one file per module, are in `2026-11-09_mutmut_raw/`; they are
unverified and unclassified. Card D-01 must rerun with the full test selection per mutant before any
D-02 work starts.

## 6a. S1 pilot on INPHARED phages (measured, wilkcomp02, in progress)

Source: `14Apr2025_genomes.fa.gz` (643,340,443 bytes) from `millardlab-inphared.s3.climb.ac.uk`.
Driver: `scripts/stress/pilot_s1/pilot.py`; store kept at `/tmp/if2453/pilot_s1/root` on the pod,
release name `rS1.pilot1`. Results so far (the rest is in `/tmp/if2453/pilot_s1/result.json` when
the run ends):

| Item | Value |
|---|---|
| records; bases | 34,076; 2,084,350,440 bp (median 45,487 bp, maximum 735,411 bp) |
| distinct sequences stored | 27,258; 1,664,109,938 bp (20% of the records are exact duplicates) |
| shards | 2 |
| ingest, 1 thread | 427 s (4.9 Mbp/s of input), peak RSS 5.6 GB |
| phage gene calling (`pyrodigal-gv`, meta mode), 200-genome probe | 0.222 Mbp/s per thread, 5× slower than the 1.1 Mbp/s of single mode on prokaryotes; building a finder per genome costs 0.008 ms and is not the cause |

The run finished on 2026-11-09 (all steps without error; manifest `551b5da7b3f53a5a73a9280df7063e13`):

| Step | Wall | Threads | Peak RSS | Result |
|---|---|---|---|---|
| ingest | 427 s | 1 | 5.6 GB | 2 shards, 0.310 bytes per bp |
| gene calling | 7,032 s (0.237 Mbp/s) | 1 | 6.8 GB | 2,398,087 genes (1.44 per kb), 1,134,550 distinct proteins (47% of genes) |
| IGRs | 17 s | 1 | | 1,542,772 regions, 939,283 of at least 30 nt |
| tapes | 38 s | 1 | 2.7 GB | 27,258 contigs, 23.7 bytes per slot, 4 row groups |
| clustering | 442 s | 16 | 6.0 GB | 534,101 NR90 and 297,082 NR50 clusters |
| tokens | 37 s | 1 | | 831,183 tokens |
| pairs (counted, not written) | 2 s | 1 | | 45,043,224 rows per level, 18.8 per gene |
| pairs, 2% sample written | 24 s | 1 | 2.9 GB for 1.85 M rows | 7.46 bytes per row |
| catalog | 9 s | 1 | | 1 MB DuckDB file over a 1.1 GB store |

Bytes per unit against the golden set: `protein_seq` 185 per distinct protein (227), `protein_occurrence`
30.5 per gene (32.2), `slot_protein` 16.1 (19.5), `cluster_member` 31.0 and `cluster_member_aln` 31.0
per row (36.5, 26.2), IGR 9.7 per row (12.7), tapes 23.7 per slot (30.4), shards 0.310 per bp (0.309).

Gene calling is 94% of the serial time. It is one process on one thread; across 16 threads it
would take about 7.5 minutes (*derived*). Card B-05 now includes parallel gene calling.

The DuckDB catalog step needs numpy for the Arrow sequence function; the `74b3e05` image lacks it
(the HEAD image has it). The pilot used the lock-pinned numpy 2.5.4 from a side directory.

## 6b. Provisional phage database with Pfam domains (measured, wilkcomp02, 2026-11-09)

Store: `/stor/scratch/Wilke/if2453/all-the-data/provisional/rS1.pilot1` (1.1 GB, manifest
`bb554e38067f23841353d94b40521654`, verified; catalog under `catalog/`). Pre-phase-B format.

| Item | Value |
|---|---|
| Pfam 38.2 scan, 17 shard processes of 4 threads, 8 at a time | 271.8 M residues; 1,984 CPU-minutes in all; 121 CPU-minutes and 30 minutes wall per full shard; 0.44 GB peak per shard; none over 250 CPU-minutes |
| rate | 7.3 CPU-minutes per million residues; 1.1 × 10¹⁰ cells per second per thread (golden set: 1.0 × 10¹⁰) |
| result | 605,564 domains on 436,807 of 1,134,550 proteins (38.5%); 6,392 models with a hit |
| family tokens | 5,402; 941,262 of 2,398,087 slots carry one (39%) |
| queries through the catalog | top families 0.09 s; terminase neighbourhoods within 3 genes over 8,148 loci 0.27 s; gene plus 500 bp flanks 0.5 s, equal to `samtools faidx` |

## 6c. Short-contig pilot, CAMI II marine sample 0 (measured, wilkcomp02, 2026-11-09)

Source `marmgCAMI2_sample_0_contigs.tar.gz`, sha256 `e7639d22…1cb55a`; store kept at
`/tmp/if2453/pilot_cami/work/root` (1.6 GB), release `rCAMI.pilot1`. Code at commit `2556ddd`
(before phase B). Every step ran without error.

| Item | CAMI | Phages (§6a) |
|---|---|---|
| contigs; bases; median length | 1,475,972; 909.5 Mbp; 294 bp (91.6% under 1 kb) | 27,258; 1.66 Gbp; 45 kb |
| ingest | 296 s (3.1 Mbp/s), 3.9 GB | 427 s, 5.6 GB |
| gene calling, 1 thread | 2,322 s (0.39 Mbp/s, 636 contigs/s), 7.8 GB peak | 7,032 s (0.24 Mbp/s) |
| genes; per contig; partial | 1,963,951; 1.33; 83.5% run off a contig end | 2,398,087; 88; not counted |
| contigs with no gene | 3.3% | |
| distinct proteins | 98.8% of genes | 47% |
| NR90 singleton clusters | 89% | |
| tapes | 824 s; 52.5 bytes per slot (70 per contig) | 38 s; 23.7 bytes per slot |
| pair rows per gene per level | 2.09; 95% of slots have their window cut by a contig end | 18.8 |
| shard; `.fai` index | 0.403 bytes per bp; 85 bytes per contig, 34% of the shard's size | 0.310; under 1% |
| `SeqStore` open | 7.2 s and 1.15 GB (4.9 s and 782 MB per million contigs; the synthetic estimate was 2 s and 414 MB) | negligible |
| one-contig tape lookup (`tapes_of`) | 181 of 181 row groups, 53 MB read | |
| 1,000 gene fetches, direct | 8.5 s, 1.01 blocks each | |
| 1,000 gene fetches, catalog macro | 15.2 s; `contig_index` scanned whole 3 times and `seq_locator` twice, 326 MB of Parquet read | |

What this shows. With short contigs the cost moves from bases to contigs: the per-contig index
files and the in-memory index dominate, and almost every neighborhood is censored. Cards B-07
(tape lookup that prunes, done), B-13 (index that is looked up, not loaded), B-14 (sequence macros
without index scans) and B-15 (batch ingest) address the measured limits; nothing here calls for a
change to the target model.

## 7. State of the pod and CI at the start (measured)

- CI is red at HEAD: `test_runner_death_is_reported_and_the_worker_does_not_outlive_it` (run 41096564184).
- Pod checkout at HEAD and clean; newest image on the pod is `74b3e05`; scratch holds only the three reference databases; `/stor/work/Wilke/if2453/all-the-data` does not exist; no release, manifest or `test_run` rows exist on the pod.
- This test left its files under `/tmp/if2453/review_s0` on the pod (about 100 MB); nothing was written to `/stor`.
