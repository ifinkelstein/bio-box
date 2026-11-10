# Import ladder with the phase-B code (wilkcomp02, 2026-11-09 and 2026-11-10)

Three stages of growing size, built with `plan/task-cards` at `f9af5b6` (before B-16, B-19 and
C-07), to see how each builder and each query scales. Not a release: a provisional store at
`/stor/scratch/Wilke/if2453/all-the-data/provisional/rL`, work in `/tmp/if2453/ladder`.

Every number is *measured* on wilkcomp02 on 2026-11-09/10 (CST) by `scripts/driver.py`, which
records each step's wall time, process-tree CPU and peak RSS with `os.wait4`, unless it is marked
*derived*. Raw records: `results/<A|B|Bp|C>/`. Genome lists: `selection/`. The driver is not
`bio-box run-step` (it predates C-07) and the sync wrote about 100 files under 64 MB to `/stor`;
both are deviations from AGENTS.md §7 that the orchestrator's brief allowed for this pilot.

## Inputs

| Stage | Source | Genomes | bp stored | Contigs | Genes | Distinct proteins |
|---|---|---|---|---|---|---|
| S1 (ladder A) | INPHARED 14Apr2025 | 27,258 | 1,664,109,938 | 27,258 | 2,398,087 | 1,134,550 |
| P1 (ladder B) | GTDB R226, complete, one per species, seed 20261109 | 250 | 1,016,625,803 | 539 | 927,920 | 918,765 |
| P2 (ladder C) | same selection, next 2,250 | 2,247 (3 failed to download) | 9,596,180,328 | 5,018 | 8,722,859 | 8,580,713 |

FASTA carries no topology, so every contig was stored as linear.

## Releases

| Release | Stages | State |
|---|---|---|
| `rL.1` | S1 | done, verified. Row-for-row equal to the old-code pilot on locators, genes, proteins, IGRs, tapes, clusters and tokens |
| `rL.2` | S1 + P1 | done, verified |
| `rL.2p` | meant to be S1 + P1 with Pfam on P1 | **contaminated, do not use for counts**: see finding 1 |
| `rL.3` | S1 + P1 + P2 | clustering was running at 01:55 CST on 2026-11-10 (unit `ladder-C`) |

## Step timings

| Step | Threads | S1 wall s / peak MB | P1 wall s / peak MB | P2 wall s / peak MB | Rate |
|---|---|---|---|---|---|
| ingest | 1 | 462.5 / 1,904 | 253.5 / 1,238 | 2,597.6 / 1,319 | 3.6 to 4.0 Mbp/s |
| gene calling | 16 | 977.0 / 1,673 (meta mode) | 86.8 / 2,609 | 787.9 / 3,444 | 1.7 Mbp/s phage; 11.7 to 12.2 Mbp/s prokaryote |
| IGRs | 1 | 15.4 / 987 | 10.5 / 662 | 53.6 / 2,178 | 73 to 134 k rows/s |
| tapes | 1 | 22.7 / 1,558 | 12.9 / 930 | 78.7 / 2,443 | 72 to 111 k slots/s |
| clustering (whole release) | 16 | 473.5 / 5,352 (1.13 M proteins) | 1,319.9 / 9,571 (2.05 M) | running (10.6 M) | see finding 3 |
| tokens | 1 | 21.2 / 1,299 | 38.1 / 1,763 | | 87 to 113 k slots/s |
| postings | 1 | 20.0 / 1,353 | 27.9 / 1,749 | | 272 to 287 k rows/s |
| catalog | 1 | 17.6 / 421 | 33.5 / 746 | | about 250 MB/s per hashing pass, two passes |

Pfam on the 918,525 proteins new in P1 (294.5 M residues, 19 shards, 8 at a time): 5,153 s wall,
2,149.9 CPU-min, 7.30 CPU-min per million residues; largest shard 119.1 CPU-min and 617 MB.

## Bytes

| Item | S1 | P1 | P2 |
|---|---|---|---|
| sequence shards, bytes per bp | 0.310 | 0.305 | 0.3055 |
| `feature`, bytes per gene | 24.9 | 26.0 | |
| `feature_protein`, bytes per gene | 40.5 | 43.0 | 42.9 |
| `protein_occurrence`, bytes per gene | 47.5 | 54.4 | 55.5 |
| tapes, bytes per slot | 42.7 | 47.9 | 48.0 |
| `protein_seq`, bytes per protein | 185 | 228.5 | 229.6 |
| `token_postings`, bytes per row | 34.6 | 36.6 | |

## Queries through each catalog

| Query | `rL.1` | `rL.2` | `rL.2p` |
|---|---|---|---|
| open catalog | 0.38 s | 0.40 s | 0.39 s |
| 1,000 random genes ±500 bp | 0.589 s | 0.656 s | 0.746 s |
| index row groups read for that | 1 | 2 | 3 |
| tapes of one contig, row groups read | 1 of 19 | 1 of 26 | 1 of 93 |
| NR50 postings lookup, row groups read | 1 of 47 | 1 of 62 | 1 of 69 |
| neighbours of the terminase family (6,200 anchors) | 0.09 s | 0.19 s | 0.21 s |

Lookups do constant work as the store grows; the sequence fetch adds one index row group per
stage (one `seq_locator` file per stage, DEC-22).

## Findings

1. **A release builder reads every stage on disk.** `rL.2p` was built while stage P2 was being
   written; its tokens, postings, `contig_dim` and manifest took in P2's files. It covers 32,815
   contigs with clusters computed without P2. `rL.1` and `rL.2` were built before P2 existed and
   are clean. A release must name its stages. Card B-25.
2. **`token_postings` is 34.6 to 36.8 bytes per row; the estimate in `postings.py` is 20.** This
   misses the estimate by 73 to 84%, over the 20% threshold of AGENTS.md §6. Card B-26.
3. **Clustering does not scale as written.** Time grew as about n^1.73 between 1.13 M and
   2.05 M proteins and memory is about 4.7 KB per protein (*derived* from the two points). At
   about 95 M proteins (*derived*, 106 Gbp of prokaryotes) that is about 450 GB, over the 200 GB
   scope. Two points are thin; `rL.3` (10.6 M) is the third, and it revises this: see "Stage C
   outcome" below (about n^1.2 in time; memory still over the scope). Card B-24 already covers the
   rewrite and is needed before S2, not after.
4. **Ingest of a large stage breaks the 300 CPU-minute rule.** It is linear at about 3.7 Mbp/s
   in one process: 42.9 CPU-min for 9.6 Gbp, so about 477 CPU-min at 106 Gbp (*derived*). A stage
   must be ingested as several batch processes. Folded into card C-04.
5. `l0/l0_manifest.parquet` is the only file an append-only check finds rewritten between
   releases (27 of 28 and 51 of 52 files identical). Card B-22 already covers it.
6. MMseqs2 at 16 threads reproduced the pilot's clusters and alignments exactly on S1
   (2,269,100 rows). One agreement is not a determinism guarantee across thread counts.
7. Gene calling in the new code is 7.2 times faster than the pilot on the same input (977 s on
   16 workers against 7,032 s on one thread), with 1.7 GB peak against 6.8 GB.

## Stage C outcome (added 2026-11-10 03:00 CST)

Stage C stopped at 02:53 CST: the unit `ladder-C` exited 1 in the postings step. `rL.3` has
ingest, gene calls, IGRs, tapes, clusters and tokens; it has no postings, `contig_dim`, manifest
or catalog, and nothing of it was synced to `/stor`.

- **Clustering at 10.63 M proteins: 9,460 s and 46.1 GB peak** (16 threads; `results/C/cluster.json`
  on the pod, not yet copied here). With the two earlier points this is about n^1.2 in time
  between 2.05 M and 10.63 M proteins, not the n^1.73 the first two points suggested, and about
  4.3 KB per protein in memory (*derived*). At 95 M proteins (*derived*) that is about 36 h and
  about 410 GB: time is workable, memory is still twice the 200 GB scope, so card B-24 stands.
- **Tokens: 216 s, 2.3 GB.**
- **Postings failed after 40 s**: DuckDB `OutOfMemoryException ... (1.8 GiB/1.8 GiB used)` inside
  `postings.run`. This is the first builder to break its bounded-memory rule on real data; the
  earlier stages had no long contigs in the majority. Card B-29.

## Projection to 106 Gbp of prokaryotes (all *derived* from the rates above)

| Step | Estimate |
|---|---|
| ingest | 8 h in one process; must be split |
| gene calling, 16 workers | 2.4 h, 1,600 CPU-min |
| IGRs, tapes, tokens, postings | 10 to 20 min each |
| catalog | 15 min per release |
| Pfam, 48 threads | 3,730 CPU-h, about 78 h |
| clustering | about 36 h and about 410 GB with the current builder (three points); over the 200 GB scope |

## How to read stage C later

`tail -n 20 /tmp/if2453/ladder/logs/C.log`; marker `C.done` or `C.failed` in the same directory;
results in `/tmp/if2453/ladder/results/C/result.json`. To resume after a failure:
`systemd-run --user --unit=ladder-C --collect -p MemoryMax=200G sh /tmp/if2453/ladder/scripts/stage.sh C`
(finished steps are skipped). The Pfam scan of P2 was not started: 2,776.8 M residues at 7.30
CPU-min per million is about 338 CPU-h (*derived*).
