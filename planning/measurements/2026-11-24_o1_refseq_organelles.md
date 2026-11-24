# Stage `O1`: NCBI RefSeq Release 238 Organellar Genomes (Card C-59, `DEC-69`)

- **Date**: 2026-11-24
- **Host**: `wilkcomp02` (container `bio-box-test_faec1f4.sif`, read-only on `/tmp/if2453/r10/root`)
- **Script**: `scripts/stress/build_stage_pod/c59_refseq_organelles.py`
- **Decisions**: `DEC-69` (5), `DEC-65` (2), `DEC-66`, `DEC-50`, `DEC-51`, `DEC-26`, `DEC-27`

## 1. Sources Fetched and Verified (`sources.yaml`)

Fetched sequentially (`1` connection, `<= 2 req/s`, `User-Agent: curl`, `0` refusals) on 2026-11-24 from NCBI RefSeq Release `238` into `/stor/scratch/Wilke/if2453/bio-box-staging/O1/raw/` (every file `>= 185 MiB`, obeying the `64 MB` `/stor` floor; zero `/tmp` scratch used for raw archives). Every file matched NCBI's official MD5 in `release-catalog/release238.files.installed` (*measured*, `wilkcomp02`, 2026-11-24):

| Source | File | Bytes | Provider MD5 (`release238.files.installed`) | Measured SHA-256 |
|---|---|---:|---|---|
| `ncbi_refseq_plastid_r238` | `plastid.1.genomic.gbff.gz` | 658,905,588 | `6087d566829a80292951c0bbc7e9ee32` | `7ff7b6e6ce358d4443fa2e3222014704fd80b7ad4c3001ca44db3b24d5adf586` |
| `ncbi_refseq_plastid_r238` | `plastid.2.genomic.gbff.gz` | 662,721,084 | `2609346fc4e4e4b1b4ac3bdb3985c0d8` | `2f09a4c39ea71865abf85f4899f23653023fb4d929687bef299736cf09918e0f` |
| `ncbi_refseq_plastid_r238` | `plastid.3.genomic.gbff.gz` | 193,833,895 | `eba5517505edf6a23b4c7fc137618f46` | `46d74a06a17f49c1d538145c6bfde3790069171606578f59287a333d8db2cd20` |
| `ncbi_refseq_mitochondrion_r238` | `mitochondrion.1.genomic.gbff.gz` | 403,081,513 | `11660db0524ea3a09468b969bb7cf31f` | `746e340bc7dd1c2a16600b91f017b68a4c1665a15dd56e2ac46b66a421dd4cf8` |

License: INSDC / NCBI RefSeq (public domain), `redistributable: true`.

## 2. RefSeq Release 238 Organellar Audit (`transl_table`, Introns, Homing Endonucleases / RTs)

Full streaming audit across all **33,451 records** (`2,918,161,894` bp = `2.918 Gbp`, *measured*, `c59_refseq_organelles.py audit`, `/tmp/if2453/c59/audit.json`, 2026-11-24). Note: in GenBank flatfiles, CDS features without an explicit `/transl_table=` qualifier use the INSDC standard genetic code (table `1`, predominantly plant mitochondria, shown as `table_1 (INSDC default)` below):

| Compartment & `transl_table` | Records | Bases (bp) | Circular | `>= 100 kb` | CDSs | Spliced CDSs | Introns | Homing/RT ORFs | Records w/ Homing/RT | Routing |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `plastid:table_11` | 15,226 | 2,325,813,312 | 15,142 | 14,966 | 1,315,896 | 201,734 | 29,594 | 14,944 | 14,385 | **`O1`** (15,218 valid DNA; 8 empty/invalid) |
| `mitochondrion:table_11` | 7 | 515,733 | 7 | 1 | 471 | 5 | 23 | 0 | 0 | **`O1`** |
| `mitochondrion:table_1` (INSDC default) | 1,089 | 270,622,450 | 934 | 571 | 38,430 | 5,589 | 3,997 | 778 | 535 | `O1_left_for_c59a.tsv.gz` |
| `mitochondrion:table_2` | 7,808 | 131,036,670 | 7,725 | 0 | 101,515 | 632 | 0 | 0 | 0 | `O1_left_for_c59a.tsv.gz` |
| `mitochondrion:table_5` | 7,437 | 117,176,399 | 7,410 | 0 | 96,457 | 36 | 0 | 0 | 0 | `O1_left_for_c59a.tsv.gz` |
| `mitochondrion:table_4` | 1,403 | 60,769,763 | 1,303 | 118 | 30,655 | 2,960 | 5,271 | 3,371 | 279 | `O1_left_for_c59a.tsv.gz` |
| `mitochondrion:table_3` | 99 | 5,528,475 | 93 | 7 | 1,451 | 309 | 328 | 61 | 17 | `O1_left_for_c59a.tsv.gz` |
| `mitochondrion:table_9` | 331 | 5,035,074 | 331 | 0 | 4,116 | 3 | 0 | 0 | 0 | `O1_left_for_c59a.tsv.gz` |
| `mitochondrion:table_16` | 7 | 577,709 | 6 | 2 | 161 | 24 | 85 | 21 | 3 | `O1_left_for_c59a.tsv.gz` |
| `mitochondrion:table_22` | 11 | 447,612 | 11 | 0 | 189 | 16 | 59 | 32 | 8 | `O1_left_for_c59a.tsv.gz` |
| `mitochondrion:table_13` | 26 | 385,510 | 26 | 0 | 334 | 0 | 0 | 0 | 0 | `O1_left_for_c59a.tsv.gz` |
| `plastid:table_4` | 5 | 224,327 | 5 | 1 | 187 | 0 | 2 | 0 | 0 | `O1_left_for_c59a.tsv.gz` |
| `mitochondrion:table_14` | 2 | 28,860 | 2 | 0 | 24 | 0 | 0 | 0 | 0 | `O1_left_for_c59a.tsv.gz` |
| **Total** | **33,451** | **2,918,161,894** | **32,995** | **15,666** | **1,589,886** | **211,308** | **39,359** | **19,207** | **15,227** | |

## 3. `pyrodigal` Recall Benchmark Against RefSeq CDSs & Homing-Endonuclease / Intron-RT ORFs

Reservoir-sampled (`<= 30` genomes per stratum, fixed seed `59`, circular origin-wrapping via `biobox.l1.genecall._genes`, *measured*, `wilkcomp02`, 2026-11-24):

| Stratum (Sampled Genomes) | Mode | Pred ORFs | All CDS Exact % | All CDS Overlap % | Unspliced Exact % | Unspliced Overlap % | Homing/RT Overlap % |
|---|---|---:|---:|---:|---:|---:|---:|
| **`plastid:table_11`** (`30`, `2,558` CDSs, `29` homing/RT) | **`meta` (`kind="plasmid"`)** | **3,202** | **64.00%** | **86.20%** | **75.82%** | **87.59%** | **100.0%** (29/29) |
| `plastid:table_11` (`30`) | `gv-meta` (`kind="virus"`) | 3,200 | 63.80% | 86.08% | 75.59% | 87.45% | 100.0% (29/29) |
| `plastid:table_11` (`30`) | `prokaryote_policy` / `single` (`kind="prokaryote"`) | 3,798 | 58.37% | 83.35% | 69.15% | 84.34% | 100.0% (29/29) |
| **`mitochondrion:table_11` + `table_1`** (`37`, `1,744` CDSs, `37` homing/RT) | `meta` (`kind="plasmid"`) | 7,364 | 67.55% | 92.95% | 75.55% | 95.64% | 94.59% (35/37) |
| `mitochondrion:table_11` + `table_1` (`37`) | `prokaryote_policy` | 9,902 | 69.15% | 93.46% | 77.34% | 96.21% | 97.30% (36/37) |
| `mitochondrion:table_4` (`30`, `559` CDSs, `13` homing/RT) | `meta` vs `single(t=4)` | 1,004 / 1,192 | 52.06% / 53.85% | 87.30% / 88.73% | 59.27% / 61.30% | 89.82% / 91.85% | 92.31% / 100.0% |
| `mitochondrion:table_3` (`30`, `474` CDSs, `42` homing/RT) | `meta` vs `single(t=3)` | 762 / 1,037 | 3.16% / 58.44% | 33.12% / 88.61% | 4.01% / 74.06% | 33.69% / 94.39% | 38.10% / 88.10% |
| `mitochondrion:table_16` (`7`, `161` CDSs, `21` homing/RT) | `meta` vs `single(t=16)` | 410 / 613 | 10.56% / 57.14% | 50.31% / 85.71% | 12.41% / 67.15% | 52.55% / 90.51% | 57.14% / 90.48% |
| `mitochondrion:table_22` (`11`, `189` CDSs, `32` homing/RT) | `meta` vs `single(t=22)` | 328 / 453 | 12.17% / 53.44% | 73.54% / 86.77% | 13.29% / 58.38% | 77.46% / 90.75% | 93.75% / 84.38% |
| `mitochondrion:table_9` (`30`, `373` CDSs) | `meta` vs `single(t=9)` | 537 / 537 | 0.27% / 0.27% | 33.51% / 33.51% | 0.27% / 0.27% | 33.33% / 33.33% | n/a |
| `mitochondrion:table_13` (`26`, `334` CDSs) | `meta` vs `single(t=13)` | 413 / 413 | 0.30% / 0.30% | 13.77% / 13.77% | 0.30% / 0.30% | 13.77% / 13.77% | n/a |
| `mitochondrion:table_5` (`30`, `388` CDSs) | `meta` vs `single(t=5)` | 394 / 394 | 0.00% / 0.00% | 11.86% / 11.86% | 0.00% / 11.86% | 11.86% / 11.86% | n/a |
| `mitochondrion:table_2` (`30`, `392` CDSs) | `meta` vs `single(t=2)` | 820 / 820 | 0.00% / 0.00% | 7.91% / 7.91% | 0.00% / 0.00% | 7.99% / 7.99% | n/a |

**Choice of `kind="plasmid"` (`meta` mode) for `O1`**: `98.3%` (`14,966 / 15,226`) of `plastid:table_11` genomes are `>= 100 kb`, so `kind="prokaryote"` invokes single-genome training (`finder.train(..., translation_table=11)`). Because plastids hold only `~80–120` CDSs across AT-rich inverted-repeat architectures, `meta` mode (`kind="plasmid"`) outperforms `single` by **+6.67 pp** on unspliced CDS exact recall (`75.82%` vs `69.15%`) and **+2.85 pp** on overall CDS overlap (`86.20%` vs `83.35%`) while achieving **100.0% (29/29)** recall on homing-endonuclease / maturase ORFs and avoiding `~18%` excess intergenic false-positive ORFs (`3,202` vs `3,798`).

**Implications for `C-59a`**: Non-table-11 fungal/protist/algal mitochondria (`table_3`, `table_4`, `table_16`, `table_22`, holding `3,485` homing/RT ORFs) require per-genome `transl_table` training (where `single(t=rec_table)` boosts `table_3` CDS overlap from `33.12% -> 88.61%` and homing/RT from `38.10% -> 88.10%`, and `table_16` from `50.31% -> 85.71%` / `57.14% -> 90.48%`), while compact `< 20 kb` animal mitochondria (`table_2`, `table_5`, `table_9`, `table_13`) cannot be trained by `single` mode and carry `0` homing-endonuclease/RT ORFs.

## 4. Staged `O1` Inputs & Build-Order Deduplication

Staged into `/tmp/if2453/s2/inputs/O1/` (`821,336,457` B = `783.3 MiB` on `/tmp`, well under the `1.5 GB` cap; *measured*, `wilkcomp02`, 2026-11-24, `/tmp/if2453/c59/o1_report.json`):

- **Records**: `15,225` table-11 DNA records (`15,218` plastids, `7` bacterial-like protist mitochondria), `75` repeated identical genomes (`DEC-27`) -> **`15,150` distinct genomes/contigs**, **`2,314,608,899` bp (`2.315 Gbp`)**, **`15,071` circular (`99.48%`)**, `kind="plasmid"`, `genetic_code=11`, byte-identical as deposited (`DEC-65 (2)`).
- **Build-order deduplication** against all 13 stored root stages (`S1, P1..P6, H1, PL1, C5, V1, G1, SPX1`) plus staged `V2, PL2, PL3, CG1, MGE1, SAT1, MGE2, ECE1`:
  - Held by earlier stages: **`1` contig** (`PL2`: `1` contig, `4,343` bp; `0` held by any other stage).
  - Net new: **`15,149` contigs (`99.993%`), `2,314,604,556` bp (`2.315 Gbp`)**.
- **Deferred for `C-59a`** (`O1_left_for_c59a.tsv.gz`, `18,226` rows, `743,688` B, sha256 `fe7a7bf41cef2247b9fc349aae4e7672faa561ce6b7272b7711f1d3dae97f100`):
  - `transl_table_2`: `7,801`; `transl_table_5`: `7,435`; `transl_table_4`: `1,367`; `transl_table_1`: `932`; `transl_table_9`: `331`; `invalid_or_empty_seq`: `214`; `transl_table_3`: `98`; `transl_table_13`: `26`; `transl_table_22`: `11`; `transl_table_16`: `7`; `transl_table_14`: `2`; `no_transl_table`: `2`.

### Staged Files and SHA-256 Checksums (`/tmp/if2453/s2/inputs/O1/`)

| File | Bytes | SHA-256 |
|---|---:|---|
| `genomes.json` | 2,663,714 | `3330a0794ea410f652876be57c58922443732c8a663895e04e929acb3c27ce90` |
| `objects.json` | 560 | `de4cf4829204f9aa495c7f3c109e1935c7d20a61d12f318cc8ec31dce1ff70e1` |
| `ncbi_refseq_plastid_r238_t11.fna.gz` | 813,338,579 | `2296566e813a687469dc9de741830c9b657a00ba62b0446a5e4574f7250836fb` |
| `ncbi_refseq_mitochondrion_r238_t11.fna.gz` | 179,347 | `2ccaf7a0386157dbba4df0337e5849adfe4378885975c9bd7090ced2551fd623` |
| `o1_metadata.tsv` | 2,430,751 | `aa0e537a67e1beef3f8f6c8ed27ebbac8db22627461837459d0f4fd686c50090` |
| `records.tsv` | 1,974,963 | `278ecafc2e74c7f448bf5c8d285bf80343cb00df80d31368d4a10e1bf7444e1a` |
| `repeated_records.tsv` | 4,731 | `e0925c97ce77f0f5b9db6eafe8877a8ff6f77ef0f47cb627af93fd561f6b0fc0` |
| `shared_records.tsv` | 96 | `f41796acf3ffa4990320e84de458c6bbf3b3d6c069fab4d23d5e62b88bcbb959` |
| `O1_left_for_c59a.tsv.gz` | 743,688 | `fe7a7bf41cef2247b9fc349aae4e7672faa561ce6b7272b7711f1d3dae97f100` |

## 5. Verification, Pre-Launch Estimate, and `tower` Mirror

- **`s2_check_inputs.py O1`**: exit `0`, `"problems": []` (`15,150` genomes, `15,071` circular, `kinds: ["plasmid"]`, `codes: [11]`, `2` `.fna.gz` objects).
- **`stage_est.py conf --stage O1`**: exit `0` (`DEPTH=annotated`, `LIBNAME=pfam`, `LIBVER=38.2`, `LIBCUT=gathering`, `CRISPR=yes`, `AFTER=SAT1`).
- **`stage_est.py before`** (read-only on `/tmp/if2453/r10/root`, *measured* 2026-11-24): `"refusals": []`, `0` held by all 13 stored root stages (exit `8` solely because queued `PL2` holding `1` contig is not yet built in `/tmp/if2453/r10/root`); estimates: `ingest` `10.386` CPU-min, `genecall` `346.334` CPU-min (`2` batches: `173.298` + `173.286` CPU-min, `2,497,694.5` genes expected, merge `0.479` CPU-min), `crispr_default` `15.011` CPU-min, `crispr_near` `17.243` CPU-min.
- **`tower` mirror (`DEC-66` `S-2`)**: `relay_to_tower.sh /tmp/if2453/s2/inputs/O1 /mnt/user/research-bulk/bio-box/sources/O1` -> `SHA256SUMS: written, 9 lines (9 verified in this run); DONE: 9 verified, 0 in place and skipped, 0 failed`.
