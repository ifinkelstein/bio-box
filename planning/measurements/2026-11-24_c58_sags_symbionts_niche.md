# Card C-58: Single-Cell SAGs (`SAG1`), Genome-Reduced Endosymbionts (`SYM1`), Post-`R226` Prokaryotes, and Niche-Habitat MAG Compendia (`NX1`..)

Date: 2026-11-24, `wilkcomp02`. Every number carries a label per `AGENTS.md` §6: *measured* (host, date, command/file), *derived* (formula and inputs), *estimate*, or *published*.
Decision rows: `DEC-69` (scope), `DEC-59` (deferred-data register), `DEC-50` & `DEC-51` (`sources.yaml` before fetch, license & `redistributable` flags), `DEC-26` & `DEC-27` (`contig_uid` ownership and multi-genome membership), `DEC-66` (`tower` mirror and `/tmp` hygiene).
Script and tests: `scripts/stress/build_stage_pod/c58_sags_symbionts_niche.py`, `tests/r1/test_c58_sags_symbionts_niche.py`, `scripts/stress/build_stage_pod/stages/SYM1.conf`.

---

## 1. Executive Summary

1. **NCBI Prokaryotic Assembly Diff vs `GTDB R226` & `GlobDB r226` (`diff-ncbi`, *measured*, `wilkcomp02`, 2026-11-24):**
   Across **`4,177,635`** NCBI RefSeq and GenBank prokaryotic `assembly_summary.txt` rows (`RefSeq Bacteria` `552,701`, `RefSeq Archaea` `3,216`, `GenBank Bacteria` `3,565,774`, `GenBank Archaea` `55,944`; `555,915` paired RefSeq/GenBank deduplicated, `408` partial skipped), **`732,388`** assemblies are already present in `GTDB R226` (`1,109,592` versionless `GCA_`/`GCF_` accessions) or `GlobDB r226` (`159,747` base accessions). The remaining assemblies split cleanly into:
   - **`SYM1` (`< 1 Mbp` CheckM-excluded genome-reduced endosymbionts & obligate symbionts):** `1,885` total assemblies deduplicated to **`323` species representatives (`229.86 Mbp`, `2,722` distinct contigs)**, **`323/323` fetched** (`0` not served), screened with `pyrodigal` (`11` vs `4` vs `25`), staged into `/tmp/if2453/s2/inputs/SYM1/` (`68.6 MB` `.fna.gz`), validated with `s2_check_inputs.py SYM1` (`problems: []`), mirrored `8/8` to `tower:/mnt/user/research-bulk/bio-box/sources/SYM1/`, and configured in `stages/SYM1.conf` (`DEPTH=annotated`, Pfam 38.2 + `--crispr`, `AFTER=ECE1`, `redistributable: true`).
   - **Post-`R226` standard prokaryotes:** `2,887,039` assemblies deduplicated to **`49,640` species representatives (`230.35 Gbp`)**, recorded in full in `/tmp/if2453/c58/diff/SYM1_post_r226_left_for_later.tsv.gz` (`2,469,478` B, sha256 `0e100f7562cc821a3f38e3a304db193928bf5955d356aaafc022adbd01e53af1`; mirrored to `tower:/mnt/user/research-bulk/bio-box/deferred_lists/SYM1/`).
2. **`GORG-Tropics` (`SAG1`) & `GORG-Dark` (`SAG2`) Single-Cell SAGs (*measured*, `wilkcomp02`, 2026-11-24):**
   - `GORG-Tropics` OSF `GORG_v1.fasta.gz` (`1,441,099,680` B, sha256 `ef097ad1dbc4a8da84cc53a52ba828dea230a82445275732085ea2052bd47910`, `CC BY-NC 4.0`, `redistributable: false`) holds **`5,537` classifier-reference marine SAGs** (`236,398` distinct contigs, **`4,941,681,237` bp = `4.942 Gbp`**, `100%` genetic code `11`, `0` conflicts). Consolidated into `/stor/scratch/Wilke/if2453/bio-box-staging/SAG1/sag1_gorg.fna.gz` (`1,503,486,943` B, sha256 `42b500ab64bd734357935120b2c6eed4815d2864cd227cc2b08e7e928770413e`) with metadata in `/tmp/if2453/c58/SAG1/` awaiting `/tmp` headroom after `r26.11.3` (`s2_check_inputs.py SAG1` `problems: []`). Only `954` contigs (`17.92 Mbp`) are held by earlier stages (`S1..ECE1, MV1, V3, LR1, GART1, GART2`) -> **`235,444` contigs / `4.924 Gbp` (`99.64%`) net new**.
   - Across ENA `PRJEB33281` (`GORG-Tropics`, `12,715` assemblies, `8.122 Gbp`), `3,160` are in `GTDB R226`, `1,927` in `GlobDB r226`, `0` in `IUIG` (`3,210` in any -> **`9,505` novel vs all three**). Across ENA `PRJNA1259050` (`GORG-Dark`, `8,088` assemblies / `9,698` read runs, `6.466 Gbp`), **`0` are in `GTDB R226`, `GlobDB r226`, or `IUIG` (`8,088` = `100.0%` novel)**, registered as follow-up stage `SAG2`.
3. **Niche-Habitat & Holobiont MAG Compendia (`NX1`..`NX4`):**
   Ordered in §5 below and `planning/CATALOG_GOAL.md` §4.3 to align cleanly with card `C-66` (which owns `NX1`: `MEER` Hadal Trench, `eHOMD v11.03`, cyanobacterial mats/biocrusts, serpentinizing subsurface, sponges, lichens, rhizosphere), followed by `NX2` (hydrothermal vents & cold seeps), `NX3` (hypersaline/soda/AMD/polar cryosphere beyond `GFS`), and `NX4` (termite/insect/coral holobionts, fermented foods beyond `cFMD`, ancient paleomicrobiomes).

---

## 2. NCBI Prokaryotic Assembly Summary Diff (`wilkcomp02`, 2026-11-24)

### 2.1 Fetched Tables (`sources.yaml`: `ncbi_prokaryote_assembly_summaries`, `ena_gorg_sag_tables`)

| File | Raw Bytes | Disk Bytes | Raw SHA-256 | Label |
|---|---:|---:|---|---|
| `ncbi_refseq_bacteria_assembly_summary.txt.gz` | `245,781,616` | `41,601,362` | `69851b5e09536b778ae9094efd6c43cd52a9aefa3a3132e7fa1c65ffad085e82` | *measured* |
| `ncbi_refseq_archaea_assembly_summary.txt.gz` | `1,497,905` | `305,950` | `52e97f905c9a3a0ef628fab36090ca999c67482d7520f74197c47b5bb4098ef5` | *measured* |
| `ncbi_genbank_bacteria_assembly_summary.txt.gz` | `1,631,746,343` | `206,852,222` | `c3623bd1ad56f6b4e1efcc0118eb30bb9db273a9e0a8026b1cfa0ae36e69a1ba` | *measured* |
| `ncbi_genbank_archaea_assembly_summary.txt.gz` | `25,450,321` | `3,577,102` | `111d7a17bf92a4f920819107117602ae854bb0daf9bb9bd5415d86bd60b3ac76` | *measured* |
| `ena_PRJEB33281_assemblies.tsv` (`12,715` rows) | `1,004,200` | `1,004,200` | `dc9761d9c7d0cadd1e00eccfcdc973347f3912ec8ec89cf23a571d439a3f268a` | *measured* |
| `ena_PRJNA1259050_assemblies.tsv` (`8,088` rows) | `676,035` | `676,035` | `2866ceed6ef7e86fb6f60df2695713b8006faa04f83132c0b02c41c4ec6fefc2` | *measured* |

### 2.2 Diff Counts Against `GTDB R226` & `GlobDB r226`

| Category | Count | Bases (bp) | Label |
|---|---:|---:|---|
| Total NCBI prokaryotic rows (`RefSeq` + `GenBank`) | `4,177,635` | — | *measured* |
| Paired RefSeq/GenBank duplicates removed | `555,915` | — | *measured* |
| Partial genome representations skipped | `408` | — | *measured* |
| Already in `GTDB R226` or `GlobDB r226` | `732,388` | — | *measured* |
| `< 1 Mbp` endosymbionts / obligate symbionts (all assemblies) | `1,885` | — | *measured* |
| ↳ **`SYM1` species representatives (1 best per `species_taxid`)** | **`323`** | `229,989,343` (summary) / `229,862,011` (FASTA) | *measured* |
| Post-`R226` standard prokaryotes (all assemblies) | `2,887,039` | — | *measured* |
| ↳ **Post-`R226` species representatives (`SYM1_post_r226_left_for_later.tsv.gz`)** | **`49,640`** | `230,348,229,209` (`230.35 Gbp`) | *measured* |

---

## 3. Stage `SYM1` (`< 1 Mbp` Genome-Reduced Endosymbionts)

Staged at `/tmp/if2453/s2/inputs/SYM1/` (`66 MiB` on disk) and mirrored `8/8` (`SHA256SUMS`-verified) to `tower:/mnt/user/research-bulk/bio-box/sources/SYM1/`. Verified with `s2_check_inputs.py SYM1` (`problems: []`, exit `0`).

### 3.1 Genetic-Code Audit (`pyrodigal` `11` vs `4` vs `25`)

| Genetic Code | Genomes | Share | Representative Taxa (*measured*, `sym1_metadata.tsv`) |
|---|---:|---:|---|
| Code `11` (standard) | `192` | `59.4%` | *Ca. Tremblaya phenacola* (`171.5 kb`), *Ca. Riesia pediculicola* (`581.9 kb`), *Rickettsiales* TMED endosymbionts |
| Code `4` (`UGA -> Trp`) | `131` | `40.6%` | *Ca. Hodgkinia cicadicola* (`150.3 kb`), *Ca. Nasuia deltocephalincola* (`144.6 kb`), *Ca. Stammera capleta* (`336.8 kb`), *Ca. Pinguicoccus supinus* (`163.2 kb`), *Ca. Organicella extenuata* (`158.2 kb`), *Ca. Shikimatogenerans* spp. (`133–163 kb`), *Mycoplasmataceae* / *Spiroplasmataceae* / *Ureaplasma* |
| Code `25` (`UGA -> Gly`) | `0` | `0.0%` | None among `< 1 Mbp` NCBI endosymbionts (code-25 CPR MAGs are captured in `GART2`, card `C-63`) |
| Within-stage shared-contig conflicts | `0` | — | *measured* |

### 3.2 `contig_uid` Ownership Overlap (`DEC-26` / `DEC-27`)

Checked against `/tmp/if2453/r10/root` (`S1..SPX1`) and staged `V2, PL2, PL3, CG1, MGE1, MGE2, SAT1, ECE1, MV1, V3, LR1, GART1, GART2`:

| Owner Stage | Shared Contigs | Shared Bases (bp) | In `SYM1.conf` (`AFTER=ECE1`) |
|---|---:|---:|---|
| `G1` | 1 | 773,982 | yes |
| `PL1` | 5 | 155,880 | yes |
| `PL2` | 6 | 295,094 | yes |
| **Subtotal (`S1..ECE1`, `SYM1.conf`)** | **12** | **1,224,956** | — |
| `LR1` (not yet built before `ECE1`) | 1 | 993,452 | note only |
| `GART2` (not yet built before `ECE1`) | 2 | 1,834,225 | note only |
| **Total held across all built + staged** | **15** | **4,052,633** | — |
| **Net new (`SYM1`)** | **2,707** (`2,710` at `AFTER=ECE1`) | **225,809,378** (`228,637,055` at `AFTER=ECE1`) | — |

### 3.3 Locked `SYM1` SHA-256 Checksums (`/tmp/if2453/s2/inputs/SYM1/`)

- `genomes.json`: `e649e4a45e38b4b24f0ec5438b7e10b4c7dd70faeba1f6e396e0499c21560527` (`323` distinct genomes, `2,742` `genome_input` rows, `101` circular contigs)
- `objects.json`: `fbcd27988f602867ba22c7f0daf238c7c59b3dd307a0b70fd44d6597af6d5635` (`redistributable: true`)
- `sym1_endosymbionts.fna.gz`: `268dcde2ed88421d24c90916c776e5b5746c12860c79a906613009828eef68d7` (`68,618,115` B)
- `sym1_metadata.tsv`: `dd5cbabc1b3955b2d49514dbd8d51d58f36887b7e110dc82250041a272c5aec0`
- `genetic_codes.tsv`: `f5645fe28e6e15d7e43cb85e84a756b4d55b40d294fc7e7c086556719d640f15`
- `records.tsv`: `e08fdb3c58c6ec11ea6ff5065eddc1249d1ae5d76c34d77da58011e7946f4c87`
- `repeated_records.tsv`: `d109495ee76e4643c181eae183c452c9a84560f97a79bbe22d6630780368c565`
- `shared_records.tsv`: `76b9b8d0a1652ad3eeed26a35894273f6e8955e01c3e57579d3e611aa2576550`

---

## 4. Stage `SAG1` (`GORG-Tropics`) & `GORG-Dark` (`SAG2`) Survey

### 4.1 Accession-Level Overlap Against `GTDB R226`, `GlobDB r226`, and `IUIG` (*measured*, `wilkcomp02`, 2026-11-24)

| Dataset (ENA BioProject) | Total Assemblies | Total Bases (Gbp) | In `GTDB R226` | In `GlobDB r226` | In `IUIG` | In Any | **Novel vs All Three** |
|---|---:|---:|---:|---:|---:|---:|---:|
| **`GORG-Tropics` (`PRJEB33281`)** | `12,715` | `8.122` | `3,160` | `1,927` | `0` | `3,210` | **`9,505` (`74.8%`)** |
| **`GORG-Dark` (`PRJNA1259050`)** | `8,088` (`9,698` read runs) | `6.466` | `0` | `0` | `0` | `0` | **`8,088` (`100.0%`)** |

### 4.2 Stage `SAG1` (`GORG-Tropics` OSF `GORG_v1.fasta.gz`)

- Consolidated FASTA: `/stor/scratch/Wilke/if2453/bio-box-staging/SAG1/sag1_gorg.fna.gz` (`1,503,486,943` B = `1.50 GB`, sha256 `42b500ab64bd734357935120b2c6eed4815d2864cd227cc2b08e7e928770413e`; `5,011,022,093` B uncompressed; `redistributable: false` under OSF `CC BY-NC 4.0`).
- Metadata directory: `/tmp/if2453/c58/SAG1/` (`genomes.json` sha256 `9667599059cddfd012bb437fa4ea87bfeeb81d8cb9f4cdeaccf4cc36fea908ad`, `objects.json` sha256 `e7977369f84893f3419c878226dc20eac3d0d4e48b5d6c6a27993f9de64b2234`, `objects_staging.json`), validated with `s2_check_inputs.py SAG1` (`problems: []`, exit `0`). Awaits `/tmp` headroom after `r26.11.3` to move into `/tmp/if2453/s2/inputs/SAG1/` and add `SAG1.conf` (`DEPTH=called`, `CRISPR=yes`).
- Counts (*measured*): **`5,537` distinct SAGs** (`0` repeated), **`236,398` distinct contigs** (`4,941,681,237` bp = `4.942 Gbp`), **`100%` genetic code `11`** (`0` code 4/25, `0` conflicts).
- `contig_uid` ownership overlap vs `S1..SPX1, V2, PL2, PL3, CG1, MGE1, MGE2, SAT1, ECE1, MV1, V3, LR1, GART1, GART2`: **`954` contigs (`17,920,126` bp)** held earlier (`GART1: 805`, `G1: 99`, `MV1: 40`, `PL3: 9`, `V2: 1`) -> **`235,444` contigs / `4,923,761,111` bp (`4.924 Gbp`, `99.64%`) net new**.
- Remainder deferred in `CATALOG_GOAL.md` §4.1: `7,178` remaining lower-completeness `PRJEB33281` SAGs (`~3.18 Gbp`, *derived*: `8.122 - 4.942`) not included in OSF `GORG_v1.fasta.gz`, plus `8,088` `GORG-Dark` (`PRJNA1259050`, `6.466 Gbp`, `100%` novel) registered as `SAG2`.

---

## 5. Ordered Niche-Habitat & Holobiont MAG Compendia Register (`NX1`..`NX4`)

Aligned with card `C-66` (`NX1`) and excluding `GlobDB r226`'s already-absorbed 13 MAG catalogues (`mOTU`, `SPIRE`, `RBG`, `TPMC`, `SMAG`, `GEM`, `GOMC`, `MGnify`, `GFS`, `SHGO`, `AMXMAG`, `cFMD`, `MRGM/HRGM2`):

| Stage | Owner Card | Compendia | Est. Size | Dedup Rule |
|---|---|---|---|---|
| **`NX1`** | `C-66` | `MEER` Hadal Trench (`OEP004067` / `deepseadb`, `7,564` piezophile MAGs `> 6,000 m`), `eHOMD v11.03` (`955` oral MAGs + `25,388` phages), cyanobacterial mats/biocrusts (`figshare.21932371`), serpentinizing subsurface (`Lost City` `zenodo.5798015`, `Samail` `zenodo.8108239`), marine sponges (`figshare.29695256`), lichens (`figshare.27054937`), rhizosphere (`zenodo.16947911`) | `15–30 Gbp` (*estimate*) | `skani` 96% ANI / 50% AF vs `GTDB R226` + `GlobDB r226`; `contig_uid` dedup |
| **`NX2`** | `C-58` follow-up | Deep-sea hydrothermal vent & cold-seep sediment MAGs (Global Cold Seep catalogs `figshare.22568107` / `figshare.32788992`; Guaymas Basin, Mid-Cayman Rise, Okinawa Trough vent MAGs) | `10–25 Gbp` (*estimate*) | `skani` 96% ANI / 50% AF vs `GTDB R226` + `GlobDB r226` + `NX1`; `contig_uid` dedup |
| **`NX3`** | `C-58` follow-up | Hypersaline/soda-lake, acid mine drainage (`PRJNA666025`/`PRJNA666095`, Iberian Pyrite Belt), and polar/permafrost/cryosphere MAGs beyond `GFS` (Southern Ocean, Antarctic dry valleys, Arctic permafrost outside `LR2`) | `10–25 Gbp` (*estimate*) | `skani` 96% ANI / 50% AF vs `GTDB R226` + `GlobDB r226` + `NX1..NX2`; `contig_uid` dedup |
| **`NX4`** | `C-58` follow-up | Termite-gut, insect, coral, and non-human holobiont MAGs + fermented-food MAGs beyond `cFMD` + ancient paleomicrobiome MAGs (Pleistocene/Neanderthal dental calculus & coprolite MAGs) | `10–20 Gbp` (*estimate*) | `skani` 96% ANI / 50% AF vs `GTDB R226` + `GlobDB r226` + `NX1..NX3`; `contig_uid` dedup |
