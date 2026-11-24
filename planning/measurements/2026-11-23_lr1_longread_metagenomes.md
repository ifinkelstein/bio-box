# Public Long-Read (`PacBio HiFi` / `ONT`) Metagenome Assemblies: Crosswalk and Tranche `LR1` (card C-56, DEC-69, DEC-70, DEC-71)

Date: 2026-11-23 to 2026-11-24, `wilkcomp02`. Every number carries a label (`AGENTS.md` §6): *measured* (host, date,
command/file), *derived* (formula and inputs), or *estimate*. Host and date of every *measured* number below are
`wilkcomp02`, 2026-11-23/24, unless stated otherwise. **No stage was built in `/tmp/if2453/r10/root`**:
`/tmp/if2453/r10/root`, `/tmp/if2453/s2/inputs`, and `/tmp/if2453/c62` were bound read-only (`:ro`) throughout.

Script: `scripts/stress/build_stage_pod/c56_longread_metagenomes.py`; tests:
`tests/r1/test_c56_longread_metagenomes.py`.

## Summary

- **Global long-read metagenome inventory** (*measured*, `c56_longread_metagenomes.py crosswalk` over `61,070` ENA
  long-read metagenomic runs, `316,919` ENA `analysis` primary metagenome assemblies, `108,493` ENA `wgs_set` rows,
  `27,219` NCBI GenBank metagenome assemblies, and `79,524` IMG/M microbiome rows):
  - **`1,902` raw long-read primary metagenome assemblies** (`1,074` ENA `ERZ` + `169` NCBI `GCA` + `659` unrestricted
    IMG/M) deduplicated by biosample and run accession to **`1,582` winners** (`320` duplicates removed):
    - **`1,050` open-access ENA/NCBI winners** (`959` ENA `ERZ` + `91` NCBI `GCA`; **`3.042 Tbp`** estimated assembled
      bases; includes Microflora Danica `PRJEB58634` `125` soil ONT assemblies = `254.7 Gbp` and Aalborg University
      `PRJEB122849` `177` assemblies = `2.546 Tbp`).
    - **`532` JGI-only unrestricted IMG/M winners** (**`257.48 Gbp`**; `160` soil/terrestrial `142.9 Gbp`, `210`
      freshwater `75.4 Gbp`, `128` wastewater/sludge `30.9 Gbp`) recorded in `LR_jgi_batch.tsv.gz` for a future PI JGI
      portal fetch batch (`DEC-50 (2)`).
- **Empirical `EMG` vs submitter-assembly finding** (*measured*, probe & initial fetch on `wilkcomp02`, 2026-11-23):
  - `EMG` (MGnify) re-assemblies of long-read-linked biosamples retain only **`11.6%`** (`ERZ2606005`: `277.2 Mbp` of
    `~2.38 Gbp`) to **`21.6%`** (`ERZ841183`: `8.15 Mbp` of `~37.8 Mbp`) of bases at `>= 5,000` bp because MGnify often
    assembles companion short-read runs, whereas submitter-deposited long-read assemblies retain **`95.3%`**
    (`GCA_024708845.1`: `544.8 / 571.9 Mbp`) to **`99.4%`** (`ERZ24888620`: `29.52 / 29.70 Mbp`). Accordingly,
    `c56_longread_metagenomes.py` prioritizes submitter-deposited assemblies (`center_name != "EMG"`) ahead of `EMG`
    re-assemblies.
- **Tranche `LR1` fetched, filtered (`>= 5,000` bp, `DEC-57`), and verified** (`s2_check_inputs.py LR1` exit `0`,
  `problems: []`):
  - **`110` submitter-deposited open assemblies** (`84` NCBI WGS + `26` ENA `ERZ`; `83` ONT, `25` PacBio, `2`
    PacBio+ONT; **`50.3%` `soil_terrestrial`** = `10.384 Gbp`), yielding **`856,719` distinct contigs (`>= 5,000` bp),
    `20,629,983,870` bp (`20.630 Gbp`)**, mean contig length **`24,080` bp** (~`1.9×` `SPX1`'s `12,683` bp), with
    **`24,100` contigs `>= 100 kb` (`7,014,431,268` bp = `7.014 Gbp` = `34.0%` of all bases**, vs `9.3%` in `SPX1`).
  - **Overlap by `contig_uid` against `S1..SPX1, V2, PL2, PL3, CG1, MV1`**: **`94` contigs (`25,736,493` bp = `0.12%`)**
    held earlier (`MV1` `87`, `P3` `3`, `P2` `2`, `P6` `2`; `0` in `SPX1`, `G1`, `V1`, `V2`, `PL1..PL3`, `CG1`, `H1`,
    `S1`, `C5`) -> **`856,625` contigs (`99.99%`) and `20,604,247,377` bp (`20.604 Gbp`, `99.88%`) net NEW**.
  - **Storage placement (`/tmp` disk guard)**: `/tmp` held `248–335 GiB` free (`< 420 GiB`) while `r26.11.3` ran, so
    `lr1_longread_ge5kb.fna` (`20,931,498,312` B = `20.93 GB`, `>= 64 MB` per `AGENTS.md` §7.3) is staged at
    `/stor/scratch/Wilke/if2453/bio-box-staging/LR1/` and compact metadata (`272 MB` total in `/tmp/if2453/c56/`) at
    `/tmp/if2453/c56/LR1/`. Moving to `/tmp/if2453/s2/inputs/LR1/` and adding `stages/LR1.conf` (`GENECALL_MODE=gv-meta`)
    awaits `B-67` merge and `/tmp >= 420 GiB` free.

## 1. Full Long-Read Metagenome Crosswalk (`1,582` deduplicated winners)

*Measured* (`c56_longread_metagenomes.py crosswalk`, `/tmp/if2453/c56/crosswalk/crosswalk_summary.json`):

| Habitat stratum (`DEC-71` priority) | Open ENA/NCBI assemblies | Open est. Gbp | JGI-only IMG/M assemblies | JGI-only Gbp |
|---|---:|---:|---:|---:|
| 1. `soil_terrestrial` | `136` | `273.32` | `160` | `142.88` |
| 2. `wastewater_sludge` | `53` | `11.07` | `128` | `30.92` |
| 3. `freshwater_aquatic` | `33` | `12.18` | `210` | `75.40` |
| 4. `marine_sediment` | `187` | `126.91` | `3` | `0.67` |
| 5. `host_associated` | `385` | `49.72` | `29` | `7.59` |
| 6. `other_unlabelled` | `256` | `2,568.88` | `2` | `0.02` |
| **Total deduplicated winners** | **`1,050`** | **`3,042.07`** | **`532`** | **`257.48`** |

By platform across the `1,050` open winners: `ont` `862` (`2,974.10 Gbp`), `pacbio` `186` (`67.82 Gbp`, of which `136`
submitter-deposited HiFi = `38.84 Gbp`), `pacbio_ont` `2` (`0.15 Gbp`). Across the `532` JGI-only winners: `pacbio` `323`
(`92.72 Gbp`), `hybrid_pacbio` `127` (`139.12 Gbp`), `ont` `62` (`22.58 Gbp`), `hybrid_ont` `20` (`3.07 Gbp`).

## 2. Tranche `LR1` Composition (`110` submitter assemblies)

*Measured* (`/tmp/if2453/c56/LR1/fetch_log.jsonl` + `/tmp/if2453/c56/LR1/input_report.json`):

| Habitat stratum | Staged assemblies | Contigs `>= 5 kb` | Bases `>= 5 kb` (bp) | Share of `LR1` bp |
|---|---:|---:|---:|---:|
| `soil_terrestrial` | `14` | `539,179` | `10,384,078,172` | `50.0%` |
| `host_associated` | `38` | `106,103` | `4,970,574,067` | `23.9%` |
| `wastewater_sludge` | `25` | `111,768` | `1,971,339,640` | `9.5%` |
| `marine_sediment` | `14` | `72,001` | `1,595,059,781` | `7.7%` |
| `freshwater_aquatic` | `16` | `16,583` | `946,339,264` | `4.6%` |
| `other_unlabelled` | `3` | `14,138` | `926,631,017` | `4.5%` |
| **Total (pre-dedup across assemblies)** | **`110`** | **`859,772`** | **`20,794,021,941`** | **`100.0%`** |
| **Distinct staged in `lr1_longread_ge5kb.fna`** | **`110`** | **`856,719`** | **`20,629,983,870`** | — |

By platform: `ont` `83` (`12.994 Gbp`), `pacbio` `25` (`7.656 Gbp`), `pacbio_ont` `2` (`0.144 Gbp`).
By source: `ena_analysis` `26` (`12.359 Gbp`), `ncbi_wgs` `84` (`8.435 Gbp`).
Contig continuity (*measured*): `24,100` distinct contigs `>= 100 kb` (`7,014,431,268` bp = `34.0%` of `LR1` bases);
mean contig length `24,080` bp. `15` transient ENA connection resets / unreachable records during fetch were logged in
`fetch_log.jsonl` and retained in `LR1_left_for_later.tsv.gz`.

## 3. Staged Files, Checksums, and `s2_check_inputs.py`

| File | Location | Bytes | sha256 |
|---|---|---:|---|
| `lr1_longread_ge5kb.fna` | `/stor/scratch/Wilke/if2453/bio-box-staging/LR1/` | `20,931,498,312` | `7164e7ed6ff7179ea62df5c39bd9fd99c10664337536867e319d15e8536d2543` |
| `genomes.json` | `/tmp/if2453/c56/LR1/` | `122,510,820` | `1b238cb8a06f5e76fdc580fbf9a4268885a6592d9e002e59e11b4fadb0474bbe` |
| `objects.json` | `/tmp/if2453/c56/LR1/` | `465` | `47d9cc6d55d392af2e6ff1a338b4dbaa8d7764729cd867638ce894aee9400c23` |
| `fetch_log.jsonl` | `/tmp/if2453/c56/LR1/` | `48,798` | `0898f1fe6da9b13bf007cd0b12ce01863b5d547e4d5a6c9d575c47205df2e76e` |
| `shared_records.tsv` | `/tmp/if2453/c56/LR1/` | `4,081` | `86f28db594f2809d89f1c08a488daffe4448daf13177d1f9bf800cd933bbf32f` |
| `LR1_left_for_later.tsv.gz` | `/tmp/if2453/c56/LR1/` | `60,636` | `de77d44334d186525c6860c569ed31819d9fc3d49a1d7be589528079ddac0320` |
| `LR_jgi_batch.tsv.gz` | `/tmp/if2453/c56/LR1/` | `18,846` | `02df650bdb58be67d0a0a8e8cd1d5a6b8fc206eaa1915619d3f05d442dfce32a` |

- **Dry check**: `s2_check_inputs.py LR1 /tmp/if2453/c56/LR1/genomes.json /tmp/if2453/c56/LR1/objects.json` -> exit `0`,
  `856,719` genomes, `856,719` `genome_input` rows, `kind: ["prokaryote"]`, `codes: [11]`, `assembly_quality: ["contig"]`,
  `redistributable: true`, `problems: []`.

## 4. `contig_uid` Overlap Against Earlier Stages (`S1..SPX1, V2, PL2, PL3, CG1, MV1`)

*Measured* (`/tmp/if2453/c56/LR1/shared_records.tsv` and `input_report.json`):

| Owner stage | Contigs held | Bases held (bp) |
|---|---:|---:|
| `MV1` | `87` | `3,272,920` |
| `P3` | `3` | `8,821,838` |
| `P2` | `2` | `8,033,710` |
| `P6` | `2` | `5,608,025` |
| All other stages (`S1, P1, P4, P5, H1, PL1, C5, V1, G1, SPX1, V2, PL2, PL3, CG1`) | `0` | `0` |
| **Total held earlier** | **`94`** | **`25,736,493`** |
| **Net NEW in `LR1`** | **`856,625`** | **`20,604,247,377`** |

## 5. Deferred Remainder (`DEC-59`, `CATALOG_GOAL.md` §4.1)

- **Open-access remainder (`LR1_left_for_later.tsv.gz`)**: `940` deduplicated open ENA/NCBI long-read primary
  metagenome assemblies (`3,014,815,847,881` bp = `3.015 Tbp` est. raw, including Microflora Danica `PRJEB58634` and
  Aalborg `PRJEB122849` remainder plus `15` transient-reset retries) ordered `lr_rank = 1..940` for `LR2`..
- **JGI-only remainder (`LR_jgi_batch.tsv.gz`)**: `532` unrestricted IMG/M long-read metagenomes (`257.48 Gbp`;
  requires PI JGI login under `DEC-50 (2)`).
- Check (*derived*): `110 + 940 = 1,050` open winners; `1,050 + 532 = 1,582` total deduplicated winners.

## 6. Estimated `--depth called --crispr --genecall-mode gv-meta` Build Cost and Next Step

*Derived* from `MV1`'s `gv-meta` rates (`2026-11-23_mv1_inputs.md` §6) and `DEC-67` (`0.85 GB/Gbp`) for `20.604 Gbp`
net new (`856,625` contigs, `20.93 GB` plain FASTA):
- `ingest`: **`~62`** CPU-min
- `genecall` (`gv-meta`, card `B-67`): **`~1,775`** CPU-min across `~8` batches
- `crispr` (default): **`~65`** CPU-min
- **Total `--depth called --crispr --genecall-mode gv-meta`**: **`~1,900`** CPU-min (~`0.4 h` wall at `80` night
  threads under `DEC-74`)
- **Root storage**: **`~17.5 GB`** (`20.604 Gbp × 0.85 GB/Gbp`, `DEC-67`)

**Next step once `B-67` is merged and `/tmp` free `>= 420 GiB`**: move `/tmp/if2453/c56/LR1/*` +
`/stor/scratch/Wilke/if2453/bio-box-staging/LR1/lr1_longread_ge5kb.fna` into `/tmp/if2453/s2/inputs/LR1/` (updating
`objects.json` path), mirror to `tower:/mnt/user/research-bulk/bio-box/sources/LR1/`, and add `stages/LR1.conf`
(`DEPTH=called`, `MIN_CONTIG_LEN=5000`, `CRISPR=yes`, `GENECALL_MODE=gv-meta`, `AFTER=MV1`).
