# C-64: Targeted `Logan` SRA Contig Bait-Pull (`LGRT1`) via `logan50` for `ART-RT` and Rare Prokaryotic `RT`s

Date: 2026-11-24  
Host: `wilkcomp02.ccbb.utexas.edu` (read-only against `/tmp/if2453/r10/root`; zero `< 64 MB` files written to `/stor`)  
Task card: `C-64` (`DEC-26`, `DEC-27`, `DEC-50`, `DEC-65`, `DEC-69`, `DEC-71`)  
Script: `scripts/stress/build_stage_pod/c64_logan_rt_bait.py`

## 1. `s3://logan-pub/p/` Layout & Seekable Index Discovery

Direct S3 XML listing (`https://logan-pub.s3.amazonaws.com/?list-type=2&prefix=p/&max-keys=100`, `2026-11-24`, *measured*):

| Object Key | Bytes | ETag | Role |
|---|---:|---|---|
| `p/logan50_c1.0_nonhuman_complete.fa.zst` | `271,342,209,328` (`252.71 GiB`) | `357cde7dd7847a6726500702e2e7125c-8087` | MMseqs2-50%-clustered Prodigal complete non-human proteins (`~2.89B` reps, *derived*) |
| `p/logan50_c1.0_nonhuman_complete_cluster-to-protein.tsv.zst` | `864,901,503,638` (`805.50 GiB`) | `8bb179a8be0cf9ca06ffd107c11cee23-6445` | Seekable zstd frames mapping `<cluster_rep> <member_protein_id>` |
| `p/logan50_c1.0_nonhuman_complete_cluster-to-protein.index.tsv` | `460,706` (`450 KiB`, `13,956` frames) | `b329406bf5ac2e9abab145b3220e2abf` | Sorted `(first_cluster_id, byte_offset)` frame index enabling `O(1)` HTTP `Range:` cluster lookup |
| `c/{acc}/{acc}.contigs.fa.zst` | per-run | per-run | `minia3` (`k = 31`) assembled contigs per SRA run accession |

All `logan50` FASTA headers use `<run_acc>_<contig_num>_<gene_num> # start # end # strand # ...` (and `c64_logan_rt_bait.py` also parses explicit `CL:<cluster_id>` tags where present).

## 2. Rare-RT Bait HMM Panel (`bait_rt12.hmm`)

Built from pinned `refdb.yaml` HMMs (`c64_logan_rt_bait.py build-bait`, `wilkcomp02`, 2026-11-24, *measured*):
- Target panel: `bait_rt12.hmm` (`12` models, sha256 `12257b51eda00d6737ac89d9881d5ed67488ce10b7220daee10d0f3d22cd77e2`)
- Inline ART-RT decoy panel: `bait_art_decoys8.hmm` (`8` `*_REV`/`*_PERM` models from `both.hmm`, sha256 `ba2eda27084d845f7d46b798766289049d63072905725d49424fb0015a20b205`, evaluated only on `both.hmm` hits)

| Model | Length (aa) | Pack (`refdb.yaml`) | Min Bitscore | Notes |
|---|---:|---|---:|---|
| `NTERM`, `CTERM`, `PROXORF`, `DISTORF` | `221, 78, 271, 110` | `lab_art_rt` (`both.hmm`) | `32.0` | `min_dom_frac >= 0.45` + 8-decoy veto |
| `clade11-A`, `clade11-B` | `228, 223` | `spire_retrons_2026` (`toro_types.hmm`) | `40.0` | Proto-ART / clade 11 retrons |
| `RVT-CRISPR-G1`, `G2`, `G3` | `204, 203, 205` | `myrt` (`rvtall.hmm`) | `45.0` | CRISPR-associated RTs |
| `RVT-G2L4` | `197` | `myrt` (`rvtall.hmm`) | `45.0` | Group II-like clade 4 RTs |
| `UG28` (`DRT9__DRT9.hmm`) | `507` | `defensefinder_rt` | `45.0` | Defense-associated RT 9 (`DRT9`/`UG28`) |
| `UG27` (`UG27__UG27.hmm`) | `420` | `defensefinder_rt` | `45.0` | ART-linked `UG27` RT clade |

## 3. Streaming Scan Throughput & Resumable Checkpoint (`c64-scan-1`)

Run detached on `wilkcomp02` (`systemd-run --user --unit=c64-scan-1 --collect`, `6` CPUs, `nice 19`, `MemoryMax=32G`, `MemorySwapMax=0`, `--max-wall-s 600`, streaming HTTP -> `zstd -dc --long=31` -> `pyhmmer` with `0` archive bytes written to `/tmp`, *measured*):

| Metric | Value | Label |
|---|---:|---|
| Compressed bytes streamed (`start_byte = 0`) | `2,626,945,024` B (`2.447 GiB` = `0.9681%` of `271.34 GB`) | *measured* |
| `logan50` protein representatives scanned | `28,375,643` (through `ERR11530903_91159_3`) | *measured* |
| Wall time / CPU time | `600.41 s` (`10.0 min`) / `1,195.82 s` (`19.93 CPU-min`) | *measured* |
| Throughput (`6` CPUs) | `47,260.1` proteins/s (`4.375` compressed MB/s) | *measured* |
| Full `271.34 GB` scan projection at `6` CPUs (`40` overnight CPUs) | `~17.2 h` wall / `~2,059` CPU-min (`~2.6 h` wall at `40` CPUs) | *derived* |
| Total decoy-filtered hits / distinct hit clusters | `4,003` hits across `1,361` distinct clusters | *measured* |
| Hits by model | `RVT-CRISPR-G3: 1,022`, `RVT-CRISPR-G2: 1,021`, `clade11-B: 400`, `RVT-G2L4: 369`, `clade11-A: 364`, `UG28: 363`, `RVT-CRISPR-G1: 359`, `UG27: 105` | *measured* |

## 4. Cluster Resolution, Contig Extraction (`>= 2,000` bp), and `contig_uid` Novelty (`LGRT1`)

1. **Cluster-to-protein resolution (`c64-resolve-1`, *measured*)**: fetched `100` seekable zstd frames of `cluster-to-protein.tsv.zst` via HTTP `Range:` requests (`171 MB` peak RSS), resolving `1,361` hit clusters into **`17,950` distinct SRA runs**. Selected top **`1,500` runs** (`5,130` target contigs) for `LGRT1` and recorded **`696,795` deferred cluster members/contigs** in `/tmp/if2453/c64/resolve/LGRT1_deferred_members.tsv.gz` (`4,785,381` B, sha256 `32eac39b95c4bea26faaefb52403a079c113142ebef7fb4ff6b31cb4c1636bc8`).
2. **Baited contig streaming (`c64-stage-2`, *measured*)**: `4` bounded connections, `max_run_bytes = 100 MB` compressed per run (`792.6 s` wall, `0` failed runs):
   - `736` runs streamed cleanly (`0` failures), `764` `> 100 MB` compressed runs deferred cleanly in `/tmp/if2453/s2/inputs/LGRT1/LGRT1_deferred_large_runs.tsv` (`33,465` B, sha256 `d08a34e1fef746ea03537e5d43963aa031c8f3205307da09b83755129179566e`).
   - Extracted **`378` distinct RT-bearing contigs `>= 2,000` bp** (`3,026,010` bp = `3.03 Mbp`, mean `8,005` bp) across **`310` distinct SRA runs**; dropped `1,737` `< 2,000` bp short contigs (`863,670` bp).
   - `DEC-65` DTR trims: `0` (none circular-flagged); within-stage duplicate `contig_uid`s: `0`.
3. **`contig_uid` novelty vs all 27 built + staged stages (`S1..SPX1, V2, PL2, PL3, CG1, MV1, V3, LR1, MGE1, MGE2, SAT1, ECE1, O1, O2, LPL1, LRV1, SYM1, SAG1, GART1, GART2, LR2, NX1`, *measured*)**:
   - Shared contigs: **`0`** (`0` bp) -> **`378` contigs / `3,026,010` bp (`100.0%`) net new**!
4. **Model breakdown across the `378` staged contigs**: `RVT-CRISPR-G3: 306`, `RVT-CRISPR-G2: 292`, `RVT-G2L4: 124`, `RVT-CRISPR-G1: 120`, `UG28 (DRT9): 110`, `clade11-B: 94`, `clade11-A: 77`, `UG27: 24`.
5. **SRA biome & holobiont distribution (`run_biomes.tsv`, sha256 `ec4dbfa6983141a5b0f70d383f407fc92e33037cace1855e99aee019b549837a`)**:
   - Gut/fecal metagenomes: `human gut metagenome` (`274`), `pig gut metagenome` (`13`), `human feces metagenome` (`10`), `gut metagenome` (`6`), `mouse gut metagenome` (`3`), `feces metagenome` (`1`).
   - Niche holobionts & environmental assemblies absent from MAG/virome compendia: insect/arthropod holobionts (*Ephemera danica* `10`, *Priocnemis perturbator* `7`, *Adalia decempunctata* `2`, *Pediacus dermestoides* `2`, *Exephanes ischioxanthus* `1`, *Amblyjoppa fuscipennis* `1`, *Cheilosia scutellata* `1`, *Cryptus armator* `1`), marine sponge holobiont (*Aplysina aerophoba* `2`), fungal holobiont (*Chlorophyllum brunneum* `1`), `soil metagenome` (`1`), `sediment metagenome` (`1`), `bacterial mixed DNA library` (`2`), `unknown`/`unidentified` (`32`), plus bacterial isolates (*Citrobacter freundii* `2`, *Legionella pneumophila* `1`, *Aeromonas caviae* `1`, *Aeromonas media* `1`).

## 5. Staged `LGRT1` Inputs & Checksums (`s2_check_inputs.py` Verified)

Staged in `/tmp/if2453/s2/inputs/LGRT1/` (`1.3 MB` total; kept on `/tmp` rather than `/stor` per `AGENTS.md` §7.3's `64 MB` `/stor` floor; total C-64 `/tmp` footprint `10.5 MB < 1.0 GB`):

| File | Bytes | SHA-256 |
|---|---:|---|
| `genomes.json` | `50,276` | `17a4bf17b2292d1a671f4a5af19600c3ac0243b8cbc921b1f84347cddbe0f79b` |
| `objects.json` | `275` | `b146639b9a45b7eb44afa0f01628c0397151c9eb787bc968d1eeb2e6c2888e73` |
| `lgrt1_baited.fna.gz` | `944,931` | `e1c6f3ebd3fbcf7a04c0dbc6f76b5abc43f4bd4a4b64b78a255b699e5a839fd9` |
| `records.tsv` | `54,068` | `0cd7b6048a11de5671273293c36fdceec2d957d45acc84332acaede91625570a` |
| `repeated_records.tsv` | `37` | `c578da7b9828e2f77260c3c655f21fe574d4e68b1dcda8f26957623e3c87b127` |
| `shared_records.tsv` | `31` | `e54c390bcc76621b9016fde92e3282f55397b40f612e61b6a4c124fdae111887` |
| `trims.tsv` | `47` | `0283575205551106e61020d371aaf46b624196d132c4793a26dad4f9879c1e3d` |

`s2_check_inputs.py LGRT1 /tmp/if2453/s2/inputs/LGRT1/genomes.json /tmp/if2453/s2/inputs/LGRT1/objects.json` exited `0` with `"problems": []` (`genomes: 378`, `distinct_genome_uids: 378`, `kinds: ["prokaryote"]`, `codes: [11]`, `assembly_quality: ["contig"]`).

## 6. Continuation Plan for `LGRT2`

1. Resume `logan50_c1.0_nonhuman_complete.fa.zst` scanning past byte `2,626,945,024` (`0.9681%` scanned; `99.0319%` remaining, `~2.6 h` wall at `40` overnight threads in `< 200` CPU-min shards).
2. Process the `764` deferred `> 100 MB` compressed SRA runs in `LGRT1_deferred_large_runs.tsv` and the `696,795` deferred cluster members in `LGRT1_deferred_members.tsv.gz`.
