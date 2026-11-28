# Staging light `tower`-landed `r26.11.6` cohorts `O4`, `GV2`, and `HROM1` on `wilkcomp02` (card C-86, DEC-80, DEC-84, DEC-85 (4), DEC-88)

- **Date and hosts:** 2026-11-27, `wilkcomp02` and `tower`. Every number carries a label (`AGENTS.md` §6): *measured* (on the named host, 2026-11-27, by the command or file named), *derived* (formula and inputs), or *estimate* (source).
- **Script and tests:** `scripts/stress/build_stage_pod/c86_o4_gv2_hrom1_inputs.py` (`census`, `o4-stage`, `gv2-stage`, `hrom1-stage`, `hrom1-ani`); `tests/r1/test_c86_o4_gv2_hrom1_inputs.py`; templates `scripts/stress/build_stage_pod/stages/templates/{O4,GV2,HROM1}.conf`.
- **Read only:** `tower:/mnt/user/research-bulk/bio-box/sources_incoming/{O3_FUNGAL,batch6b_ncbi/O4_ORG,batch5_gv/GV_*,batch6b_oral/HROM_ORAL,batch5_proks/CODE4_GLMC,batch6b_ncbi/CODE4_ACTINO1}/`, the pod data root `/tmp/if2453/r10/root` (`13` stages `S1..SPX1`), `/tmp/if2453/s2/inputs/` (`19` `r26.11.4` stages `V2..GART2`), the `24` staged `r26.11.5` stages (`GRARE1..GBRMGD`, including `/tmp/if2453/c83/meta/{RMD1,GBRMGD}`), card `C-74`'s index `/tmp/if2453/c74/index/`, and card `C-75`'s archaeal references (`/stor/scratch/Wilke/if2453/bio-box-staging/GRARE/globdb_r226_archaea.fna.gz`, `/tmp/if2453/c75/meta/archaea/archaea_list.tsv`). Nothing was written to `/tmp/if2453/r10/root`, `/tmp/if2453/s2/inputs/`, or `tower`'s `sources_incoming/`; nothing on `/stor` or `tower` was deleted or overwritten.
- **Written (`AGENTS.md` §7.3 64 MiB floor strictly enforced):** raw downloads and staged FASTA objects `>= 67,108,864` B (`64 MiB`) under `/stor/scratch/Wilke/if2453/bio-box-staging/{O4,GV2,HROM1}/` (`10` files, smallest `72,032,344` B, *measured* via `find ... -type f -printf '%s\t%p\n'`); all raw downloads `< 64 MiB`, `< 64 MiB` staged FASTA objects (`o4_fungal_mitogenomes.fna.gz`, `hrom1_code4_actino1.fna.gz`), and stage metadata under `/tmp/if2453/c86/{dl,meta}/{O4,GV2,HROM1}/`; mirrored gently (`--streams 1`, `--dry-run` first, sequential) to `tower:/mnt/user/research-bulk/bio-box/sources/{O4,GV2,HROM1}/`.

---

## 1. Gentle single-stream raw relay from `tower` to `wilkcomp02` (`DEC-88 (1)`)

Checked `pgrep -af relay_to_` locally (`0` active relays, *measured*) and `ssh root@tower uptime` before relaying. All `14` source directories (`20` size-split `SHA256SUMS` lists) were relayed sequentially one directory at a time with `scripts/orchestrator/relay_to_pod.sh --dry-run --streams 1` followed by `--streams 1` (2026-11-27 `15:56` to `16:05 CST`, *measured*): **`272` files (`5,088,493,144` B = `5.088 GB`, *derived* sum of the `20` dry-run byte totals against the card's `5.09 GB` *derived*, `-0.03%` difference), `272` verified, `0` failed**.

| Stage | Subcohort / directory on `tower` | Files (*measured*) | Bytes (*measured*) | Pod destination (`>= 64 MiB` on `/stor`, `< 64 MiB` on `/tmp`) |
|---|---|---:|---:|---|
| `O4` | `O3_FUNGAL` | 2 | 137,723,255 | `Assembled_fungal_mitochondrial_genomes.fasta` (`136,853,604` B) on `/stor/.../O4/dl/O3_FUNGAL/`; `All_fungal_mitochondrial_genomes_accessions.tsv` (`869,651` B) on `/tmp/if2453/c86/dl/O4/O3_FUNGAL/` |
| `O4` | `batch6b_ncbi/O4_ORG/O4_CRYP` | 8 | 43,723,571 | `/tmp/if2453/c86/dl/O4/O4_CRYP/` (`4` paired `.fna`/`.gbff` chunks) |
| `O4` | `batch6b_ncbi/O4_ORG/O4_PROT_MT` | 42 | 265,728,638 | `/tmp/if2453/c86/dl/O4/O4_PROT_MT/` (`21` paired `.fna`/`.gbff` chunks) |
| `O4` | `batch6b_ncbi/O4_ORG/O4_PROT_PT` | 22 | 912,048,408 | `/tmp/if2453/c86/dl/O4/O4_PROT_PT/` (`11` paired `.fna`/`.gbff` chunks; max file `66,634,503` B `< 64 MiB`) |
| `O4` | `batch6b_ncbi/O4_ORG/O4_WEIRD_MT` | 42 | 43,140,331 | `/tmp/if2453/c86/dl/O4/O4_WEIRD_MT/` (`21` paired `.fna`/`.gbff` chunks) |
| `GV2` | `batch5_gv/GV_POLAR_2026` | 1 | 589,306,100 | `GVMAGs_from_polar_regions.zip` on `/stor/.../GV2/dl/GV_POLAR_2026/` |
| `GV2` | `batch5_gv/GV_POLAR_PITOT` | 2 | 427,163,351 | `GVMAGs_fna.tar.gz` (`427,016,646` B) on `/stor/.../GV2/dl/GV_POLAR_PITOT/`; `GVMAGS_Phylogenetic_tree.txt` (`146,705` B) on `/tmp/if2453/c86/dl/GV2/GV_POLAR_PITOT/` |
| `GV2` | `batch5_gv/GV_GKS` | 1 | 72,032,344 | `GKS_genomes.tar.gz` on `/stor/.../GV2/dl/GV_GKS/` |
| `GV2` | `batch5_gv/GV_POLAR_GEVE` | 13 | 154,723,944 | `all_ICEL_viral_contigs.zip` (`122,586,105` B) on `/stor/.../GV2/dl/GV_POLAR_GEVE/`; `12` excised VR `.fna` + metadata files (`32,137,839` B) on `/tmp/if2453/c86/dl/GV2/GV_POLAR_GEVE/` |
| `GV2` | `batch5_gv/GV_POLAR_MARINE` | 98 | 22,253,901 | `/tmp/if2453/c86/dl/GV2/GV_POLAR_MARINE/` (`98` `.fa` GVMAG files) |
| `GV2` | `batch5_gv/GV_WWTP` | 2 | 17,382,066 | `/tmp/if2453/c86/dl/GV2/GV_WWTP/` (`genomes.tar.gz` `17,379,171` B + `README.md` `2,895` B) |
| `HROM1` | `batch6b_oral/HROM_ORAL` | 2 | 2,285,109,613 | `HROM_Genome.All.tar.gz` (`2,283,175,292` B) on `/stor/.../HROM1/dl/HROM_ORAL/`; `HROM-Species-metadata.tsv` (`1,934,321` B) on `/tmp/if2453/c86/dl/HROM1/HROM_ORAL/` |
| `HROM1` | `batch5_proks/CODE4_GLMC` | 3 | 107,727,903 | `merged_MAGs.zip` (`107,629,084` B) on `/stor/.../HROM1/dl/CODE4_GLMC/`; `metadata_mags.csv` + `README.txt` (`98,819` B) on `/tmp/if2453/c86/dl/HROM1/CODE4_GLMC/` |
| `HROM1` | `batch6b_ncbi/CODE4_ACTINO1` | 36 | 13,104,331 | `/tmp/if2453/c86/dl/HROM1/CODE4_ACTINO1/` (`MANIFEST_CODE4_ACTINO1_PRJNA1366607.tsv` `8,151` B + `35` `.fna.gz` assemblies `13,096,180` B) |
| **Total** | **14 directories** | **274** (*derived*; `272` non-sum files) | **5,091,167,756** (*derived*: `O4` `1,402,364,203` B + `GV2` `1,282,861,706` B + `HROM1` `2,405,941,847` B) | |

---

## 2. Raw archive census & cohort verification (`STAGING_WORKER_RULES.md`)

Inspected every relayed archive and flatfile on `wilkcomp02` by magic bytes (`gzip`, `zip`, `tar`, `gbff`, `fasta`), member paths, and IUPAC nucleotide alphabet (`probe_c86_raw.py` -> `/tmp/if2453/c86/work/{o4_probe,gv2_probe}.json`, `check_gks_and_hrom.py`, and `check_gks_pfam.py`; *measured*, `wilkcomp02`, 2026-11-27):

| Cohort / subcohort | Raw contents (*measured*) | What enters staging vs left out (`left_out.tsv` / `skipped_members.tsv`) |
|---|---|---|
| `O4` / `O4_CRYP` | `4` paired `.fna`/`.gbff` chunks: `643` records, `15,273,263` bp (`0` `.fna`/`.gbff` sequence mismatches); `641` plastids + `2` chromatophores (`0` nucleomorphs and `0` nuclear chromosomes in this NCBI query) | **621 staged**; **`22` left out**: `2` RNA (`KM062115.1`, `KM062117.1`), `2` mRNA (`DQ017796.1`, `DQ017797.1`), and `18` *Balanophoraceae* holoparasitic plastid records annotated with NCBI `transl_table=32` (`unsupported_transl_table_32`, unsupported by `pyhmmer.easel.GeneticCode` and `pyrodigal`) |
| `O4` / `O4_PROT_MT` | `21` paired `.fna`/`.gbff` chunks: `4,017` mitochondrial records, `90,476,697` bp (`0` `.fna`/`.gbff` mismatches) | **4,009 staged**; **`8` left out**: `6` RNA kinetoplast transcript records (`MT762281.1..MT762287.1`) and `2` artificial sequences (`CP034511.1`, `CP034536.1`, `not_eukaryota`) |
| `O4` / `O4_PROT_PT` | `11` paired `.fna`/`.gbff` chunks: `2,067` records (`2,066` plastids, `1` chromatophore), `306,525,925` bp (`0` `.fna`/`.gbff` mismatches) | **all 2,067 staged** (`0` left out) |
| `O4` / `O4_WEIRD_MT` | `21` paired `.fna`/`.gbff` chunks: `4,036` mitochondrial/kinetoplast records, `13,370,311` bp (`0` `.fna`/`.gbff` mismatches) | **3,635 staged**; **`401` left out**: `156` RNA + `245` mRNA mitochondrial/kinetoplast transcript/EST records (`moltype_not_dna`) |
| `O4` / `O3_FUNGAL` | `Assembled_fungal_mitochondrial_genomes.fasta`: `2,695` circular fungal mitogenomes (`136,297,649` bp, `0` non-IUPAC records), all `2,695` present in `All_fungal_mitochondrial_genomes_accessions.tsv` (`4,039` rows total: Ascomycota `1,679`, Basidiomycota `634`, Mucoromycota `270`, Zoopagomycota `96`, Chytridiomycota `5`, Blastocladiomycota `3`, `NA` `8`) | **all 2,695 staged** (`0` left out; exact `genome_uid` repeats deduplicated via `DEC-27`) |
| `GV2` / `GV_POLAR_2026` | `GVMAGs_from_polar_regions.zip`: `7,056` entries (`8` directories, `5,481` DNA GVMAG `.fna` files = `466,864` contigs / `1,617,914,493` bp, `1,565` predicted-protein `.proteins.faa` files, `1` `.tsv`, `1` `.xlsx`) | **all 5,481 DNA `.fna` bins staged**; `1,565` `.proteins.faa` + `2` tables skipped in `skipped_members.tsv` |
| `GV2` / `GV_POLAR_PITOT` | `GVMAGs_fna.tar.gz`: `3,304` DNA GVMAG `.fna` files (`264,717` contigs, `1,363,535,396` bp) + `GVMAGS_Phylogenetic_tree.txt` | **all 3,304 DNA `.fna` bins staged**; `1` tree `.txt` skipped in `skipped_members.tsv` |
| `GV2` / `GV_GKS` | `GKS_genomes.tar.gz`: `657` single-cell `tomas_scMAGs/*.fna` (`5,115` contigs, `120,542,242` bp) + `272` alpine-lake `all_MAGs/*.fltr.fa` (`3,563` contigs, `123,008,360` bp). **Pfam 38.2 HMM oracle check (`check_gks_pfam.py`, *measured*):** `25/25` sampled `all_MAGs` (`9,029` proteins) carry NCLDV hallmark domains (`50` `Capsid_NCLDV`, `23` `DNA_pol_B`, `22` `D5_N`, `7` `Pox_A32`) and `0/25` carry bacterial `Ribosomal_S5`/`Ribosomal_S8`/`SecY` (`23/25` of `tomas_scMAGs` carry the same NCLDV hallmarks, `0/25` bacterial) | **all 929 (`657 + 272`) GVMAG bins staged** (`all_MAGs` confirmed genuine giant viruses, not host/bacterial bins) |
| `GV2` / `GV_POLAR_GEVE` | `9` excised confirmed viral region `.fna` files (`493` VRs, `29,288,257` bp across `9` polar algae, incl. `ICEL_confirmed_VR_noflank.fna` `404` VRs / `25,811,477` bp) + `all_ICEL_viral_contigs.zip` (`170` whole nuclear genome scaffolds of *Chlamydomonas* sp. ICE-L `JACBWV010000xxx.1`, `463,777,451` bp, max `162,574,961` bp, whose integrated viral coordinates are listed in `coords_all_viral_contigs`) + `3` documentation/spreadsheet files | **all 493 excised confirmed VRs staged**; **`170` whole *Chlamydomonas* sp. ICE-L nuclear scaffolds left out** in `left_out.tsv` (`eukaryote_host_nuclear_contig_c60`) for card `C-60`; `3` non-FASTA files skipped |
| `GV2` / `GV_POLAR_MARINE` | `98` `.fa` polar marine GVMAGs (`1,811` contigs, `22,167,844` bp) | **all 98 bins read** (`DEC-27` deduplicates exact `genome_uid` overlaps with `GV_POLAR_2026`) |
| `GV2` / `GV_WWTP` | `genomes.tar.gz`: `61` polished NCV genomes in `genomes/ncv/contigs/` (`17,526,882` bp) + `29` VPH/PLV genomes in `genomes/vph_plv/contigs/` (`415,304` bp) + `859` Prokka annotation/protein/graph files (`87` Prokka `.fna`/`.ffn` duplicates, `177` `.faa`, `595` tables/graphs) + `README.md` | **all 90 (`61` NCV + `29` VPH/PLV) contig assemblies staged**; `860` Prokka/protein/graph/`README.md` files skipped in `skipped_members.tsv` |
| `HROM1` / `HROM_ORAL` | `HROM_Genome.All.tar.gz` + `HROM-Species-metadata.tsv`: `5,113` species representatives (`3,426` `NC/HQ` [`486` Isolate + `2,940` MAG, incl. `1,402` `p__Patescibacteriota` CPR of which `1,137` are novel unnamed species] + `1,687` `MQ`) | **all 3,426 `NC/HQ` oral species reps staged** (`--hrom-hq-only`, matching the card's exact `"3,426 HQ oral reps"`); **`1,687` `MQ` species reps left out** in `left_out.tsv` (`hrom_non_hq:MQ`) |
| `HROM1` / `CODE4_GLMC` | `merged_MAGs.zip`: `272` leafhopper symbiont/gut MAGs; UTF-8 BOM `metadata_mags.csv` has `337` rows (`272` matching the zip + `65` `Sample == "NCBI"` comparison rows without FASTA in the zip) | **all 272 zip MAGs staged** (`271` distinct; `GLMC_Bin334` is an exact `genome_uid` repeat of `GLMC_Bin333` in `repeated_records.tsv`) |
| `HROM1` / `CODE4_ACTINO1` | `MANIFEST_CODE4_ACTINO1_PRJNA1366607.tsv` + `35` `.fna.gz` INSDC assemblies (`31` `Scaffold`, `4` `Contig`; *Equivita* / *Tapirivita* / *Gorillivita*, BioProject `PRJNA1366607`) | **all 35 INSDC assemblies staged** (`0` left out) |

---

## 3. Stage `O4`: Protist organelles & fungal mitogenomes (`DEC-65 (2)`, `DEC-80`, `DEC-85 (4)`)

**Rules (`stage_o4`):**
- INSDC organellar records (`O4_ORG`) keep deposited topology (`LOCUS ... circular`) and sequence bytes untouched (`DEC-65 (2)`), `kind="prokaryote"` routed per genome by NCBI `/transl_table` where annotated (`7,745` records) or screened across organellar genetic codes (`1, 3, 4, 11, 16` via `pyrodigal` single-mode coding-density/ORF-length gain on `>= 5 kb` unannotated records, or lineage default on `< 5 kb` kinetoplast minicircles/fragments), verified via `biobox.search.orfs.check_code`.
- `O3_FUNGAL` records (`2,695` circular mitogenomes) are screened across fungal mitochondrial tables `1, 3, 4, 16` (`screen_fungal_mt_code`).
- Split into two staged FASTA objects so INSDC records retain `redistributable: true` while `O3_FUNGAL` carries `redistributable: false` (`DEC-80 (3)`: `sources.yaml` `tower_o3_fungal` notes an unresolved audit licence question).

| Metric | Value | Label & source |
|---|---:|---|
| Qualifying DNA records written | 13,027 (`O4_PROT_MT` 4,009; `O4_WEIRD_MT` 3,635; `O3_FUNGAL` 2,695; `O4_PROT_PT` 2,067; `O4_CRYP` 621) | *measured* (`wilkcomp02`, 2026-11-27, `o4_stage_report.json`) |
| Records left out (`left_out.tsv`) | 431 (`moltype_not_dna:mRNA` 247; `moltype_not_dna:RNA` 164; `unsupported_transl_table_32` 18; `not_eukaryota` 2) | *measured* (`o4_stage_report.json`) |
| Exact `genome_uid` repeats (`repeated_records.tsv`, `DEC-27`) | 2,083 | *measured* (`o4_stage_report.json`) |
| Distinct genomes (`genomes.json`) / distinct contigs | 10,944 / 10,944 (`7,420` circular) | *measured* (`o4_stage_report.json`) |
| Compartments of distinct genomes | `mitochondrion` 8,478; `plastid` 2,463; `chromatophore` 3 | *measured* (`o4_stage_report.json`) |
| Genetic codes of distinct genomes | `4`: 7,157; `11`: 2,424; `1`: 1,181; `16`: 119; `3`: 38; `22`: 20; `2`: 4; `23`: 1 | *measured* (`o4_stage_report.json`) |
| Distinct bases | 544,499,955 bp (0.544 Gbp, *derived*) | *measured* (`o4_stage_report.json`) |
| Held across 56 earlier stages (`S1..GBRMGD`, exact `contig_uid`) | **1,461 contigs, 148,539,583 bp** (*derived* sum): `O1` 815 (`122,750,901` bp), `O2` 643 (`25,716,223` bp), `LPL1` 3 (`72,459` bp); all 53 other stages `0` | *measured* (`o4_stage_report.json`) |
| **Net-new contigs / net-new bases** | **9,483 contigs / 395,960,372 bp** (`0.396 Gbp <= 1.5 Gbp` -> `DEPTH=annotated`, `LIBRARY=pfam 38.2 gathering`, `CRISPR=yes`) | *measured* (`o4_stage_report.json`) |
| `s2_check_inputs.py O4` (full byte pass) | `exit: 0`, `problems: []`, `genomes: 10,944`, `genome_input_rows: 10,944`, `codes: [1, 2, 3, 4, 11, 16, 22, 23]` | *measured* (`o4_stage_report.json`) |

---

## 4. Stage `GV2`: Polar, freshwater, alpine-lake, and WWTP giant viruses (`DEC-65 (1)`, `DEC-80`, `DEC-88`)

**Rules (`stage_gv2`):**
- Multi-contig GVMAG bins (`GV_POLAR_2026`, `GV_POLAR_PITOT`, `GV_GKS`, `GV_POLAR_MARINE`, `GV_WWTP` multi-record bins) enter as one `kind="virus"` genome per bin with contigs whole; single-record viral genomes (`GV_POLAR_GEVE` excised VRs and single-contig `GV_WWTP` NCV/VPH/PLV assemblies) apply `DEC-65 (1)` (`20..1,000` bp exact DTR trimmed once and marked circular: `1` exact `98` bp trim in `trims.tsv`).
- Sampled every `50`th staged genome across all `6` subcohorts (`206` genomes, `63,848,184` bp) with `check_gv_code25` (`sample_gv2_c25.py` -> `gv2_c25_sample.json`, *measured*): `0` code-4-dominant genomes and **`0` code-25-rescued genomes**, fixing `GENECALL_MODE=gv-meta`.
- Net-new bases after deduplicating against all `56` earlier stages (`S1..GBRMGD` + `O4`) measure **`2,267,601,494` bp (`2.268 Gbp > 1.5 Gbp`)**, so `GV2` takes **`DEPTH=called`** (`CRISPR=yes`) under the `1.5 Gbp` depth threshold (matching `DPN1` at `1.792 Gbp` in card `C-78`). `redistributable: false` (`DEC-80 (3)`).

| Metric | Value | Label & source |
|---|---:|---|
| Viral genomes read by subcohort | **10,395**: `GV_POLAR_2026` 5,481; `GV_POLAR_PITOT` 3,304; `GV_GKS` 929; `GV_POLAR_GEVE` 493; `GV_POLAR_MARINE` 98; `GV_WWTP` 90 | *measured* (`wilkcomp02`, 2026-11-27, `gv2_stage_report.json`) |
| Archive members skipped (`skipped_members.tsv`) | 2,512 (`1,565` `.proteins.faa` in `GV_POLAR_2026`, `859` Prokka/protein/graph files + `README.md` in `GV_WWTP`, `1` tree in `GV_POLAR_PITOT`, `3` docs/tables in `GV_POLAR_GEVE`, plus directory/non-FASTA entries) | *measured* (`gv2_stage_report.json`) |
| Records left out (`left_out.tsv`) | 170 (`eukaryote_host_nuclear_contig_c60`: *Chlamydomonas* sp. ICE-L whole nuclear scaffolds in `all_ICEL_viral_contigs.zip`, `463,777,451` bp) | *measured* (`gv2_stage_report.json`, `gv2_probe.json`) |
| Exact `genome_uid` repeats (`repeated_records.tsv`, `DEC-27`) | 111 | *measured* (`gv2_stage_report.json`) |
| Distinct genomes (`genomes.json`) / `genome_input` rows | 10,284 (`1` circular) / 740,762 (`742,653` FASTA records; `33,352` contigs shared across `>= 2` bins, `32` repeated within a bin) | *measured* (`gv2_stage_report.json`) |
| Distinct contigs / distinct bases | 689,523 / 3,155,057,728 bp (3.155 Gbp, *derived*) | *measured* (`gv2_stage_report.json`) |
| Held across 56 earlier stages + `O4` (exact `contig_uid`, build order) | **44,936 contigs, 887,456,234 bp** (*derived* sums): `H1` 1 (`610,033` bp), `SPX1` 34 (`437,114`), `V2` 11,526 (`334,992,906`), `PL3` 6 (`53,056`), `MV1` 3,956 (`146,238,576`), `V3` 11 (`1,550,719`), `GART2` 2 (`8,211`), `GRARE1` 134 (`1,389,909`), `LR2` 1 (`13,145`), `GV1` 23,220 (`299,852,043`), `MIR1` 925 (`7,720,773`), `MIR2` 156 (`1,414,265`), `MV2` 934 (`39,816,219`), `MVEUKV` 4,030 (`53,359,265`); `O4` `0` | *measured* (`gv2_stage_report.json`) |
| **Net-new contigs / net-new bases** | **644,587 contigs / 2,267,601,494 bp** (`2.268 Gbp > 1.5 Gbp` -> `DEPTH=called`, `GENECALL_MODE=gv-meta`, `CRISPR=yes`) | *measured* (`gv2_stage_report.json`) |
| `s2_check_inputs.py GV2` (full byte pass) | `exit: 0`, `problems: []`, `genomes: 10,284`, `genome_input_rows: 740,762`, `kinds: ["virus"]`, `gene_call_mode: {"virus": "gv-meta"}` | *measured* (`gv2_stage_report.json`) |

---

## 5. Stage `HROM1`: Oral/CPR HQ species reps + code-4 sets (`DEC-80`, `DEC-84 (2)`, `DEC-88`)

**Rules (`stage_hrom1` & `hrom1_ani`):**
- Stages all `3,426` `NC/HQ` species representatives of `HROM_ORAL` (`--hrom-hq-only`, recording the `1,687` `MQ` reps in `left_out.tsv`), all `272` leafhopper symbiont/gut MAGs of `CODE4_GLMC`, and all `35` INSDC assemblies of `CODE4_ACTINO1` (`PRJNA1366607`).
- Screens every genome across genetic codes `11`, `4` (`UGA -> Trp`), and `25` (`UGA -> Gly`) via single-mode `pyrodigal` (`screen_prokaryote_code`, `4` workers), resolving `6` within-stage shared-contig code conflicts deterministically (`c58._resolve_within_stage_code_conflicts`).
- Records GTDB R232 and GlobDB r232 membership via `/tmp/if2453/c74/index/` in `r232_status.tsv` (**`0` genomes removed**).
- Runs `skani dist -t 4 --min-af 50` (`ANI >= 96%`, larger `AF >= 50%`) against `NX1` (`893` refs incl. eHOMD v11.03), `SYM1` (`323` endosymbiont refs), `SYM2` (`625` refs), `DPN1` (`1,840` refs), `WET1` (`4,164` refs), `GBRMGD` (`5,283` refs), and card `C-75`'s `r226_archaea` of `HROM1` genera (`32` refs), recording all `41` hits in `hrom1_ani_hits.tsv` and `ani_summary.json` without cutting any genome (`DEC-84 (2)`).
- Split into two staged FASTA objects: `hrom1_oral_glmc.fna.gz` (`redistributable: false`, `DEC-80 (3)`) and `hrom1_code4_actino1.fna.gz` (`INSDC`, `redistributable: true`).

| Metric | Value | Label & source |
|---|---:|---|
| Input genomes staged | **3,733** (`HROM_ORAL` `3,426` HQ; `CODE4_GLMC` `272`; `CODE4_ACTINO1` `35`) | *measured* (`wilkcomp02`, 2026-11-27, `hrom1_stage_report.json`) |
| Genomes left out (`left_out.tsv`) | 1,687 (`hrom_non_hq:MQ` from `HROM_ORAL`) | *measured* (`hrom1_stage_report.json`) |
| Exact `genome_uid` repeats (`repeated_records.tsv`, `DEC-27`) | 1 (`GLMC_Bin334` == `GLMC_Bin333`) | *measured* (`hrom1_stage_report.json`, `repeated_records.tsv`) |
| Distinct genomes (`genomes.json`) / `genome_input` rows | 3,732 (`HROM_ORAL` 3,426; `CODE4_GLMC` 271; `CODE4_ACTINO1` 35) / 267,817 (`267,818` total records; `736` contigs in `>= 2` genomes) | *measured* (`hrom1_stage_report.json`) |
| Screened genetic codes (`3,733` input genomes) | **`11`: 3,638** (`HROM_ORAL` 3,398; `CODE4_GLMC` 237; `CODE4_ACTINO1` 3); **`4`: 76** (`CODE4_GLMC` 35; `CODE4_ACTINO1` 32; `HROM_ORAL` 9); **`25`: 19** (`HROM_ORAL` 19); `code_conflicts` resolved: `6` | *measured* (`hrom1_stage_report.json`) |
| `r232_status.tsv` counts (none removed) | `not_in_r232`: **35** (all `35` post-R232 `CODE4_ACTINO1` `GCA_056388515..9195` assemblies); `not_visible`: **3,698** (`HROM_ORAL` and `CODE4_GLMC` zip MAGs carry study IDs without `GCA_`/`GCF_`); `in_gtdb_r232`: `0`; `in_globdb_r232`: `0` | *measured* (`hrom1_stage_report.json`) |
| `skani >= 96%` ANI / `>= 50%` AF hits (`hrom1_ani_hits.tsv`, recorded not cut) | **41 of 3,731** query genome files (`1.10%`, *derived*): `SYM1` **23** (`13` `CODE4_GLMC` leafhopper symbionts + `10` `HROM_ORAL`), `NX1` **14** (`10` `ehomd:*`, `3` `rhizo_*`, `1` `lichen:*`), `WET1` **2**, `GBRMGD` **1**, `r226_archaea` **1** (`HROM_Genome_0002` vs `GCA_028561875`, `98.47%` ANI); `DPN1` `0`, `SYM2` `0` (`506.6` CPU-s total, *derived* sum) | *measured* (`hrom1_ani_report.json`, `hrom1_ani_hits.tsv`) |
| Distinct contigs / distinct bases | 266,916 / 5,883,730,291 bp (5.884 Gbp, *derived*) | *measured* (`hrom1_stage_report.json`) |
| Held across 56 earlier stages + `O4, GV2` (exact `contig_uid`, build order) | **1,388 contigs, 229,738,373 bp** (*derived* sums): `P1` 4 (`6,575,312` bp), `P2` 15 (`25,722,794`), `P3` 23 (`29,890,477`), `P4` 35 (`73,757,741`), `P5` 26 (`44,213,403`), `P6` 2 (`4,536,307`), `PL1` 12 (`745,750`), `C5` 3 (`8,126,875`), `V2` 1 (`114,551`), `PL2` 4 (`551,531`), `PL3` 45 (`2,915,453`), `MGE1` 2 (`2,547`), `SYM1` 25 (`2,442,941`), `MV1` 10 (`420,882`), `GART1` 138 (`692,870`), `GART2` 1,037 (`28,889,162`), `GRARE2` 6 (`139,777`); `O4` `0`, `GV2` `0` | *measured* (`hrom1_stage_report.json`) |
| **Net-new contigs / net-new bases** | **265,528 contigs / 5,653,991,918 bp** (`5.654 Gbp > 1.5 Gbp` -> `DEPTH=called`, `CRISPR=yes`) | *measured* (`hrom1_stage_report.json`) |
| `s2_check_inputs.py HROM1` (full byte pass) | `exit: 0`, `problems: []`, `genomes: 3,732`, `genome_input_rows: 267,817`, `kinds: ["prokaryote"]`, `codes: [4, 11, 25]` | *measured* (`hrom1_stage_report.json`) |

---

## 6. Staged files (`wilkcomp02` and `tower` mirror)

Every file on `/stor/scratch/Wilke/if2453/bio-box-staging/` is `>= 67,108,864` B (`AGENTS.md` §7.3; verified via `find /stor/scratch/Wilke/if2453/bio-box-staging/{O4,GV2,HROM1} -type f -printf '%s\t%p\n'`, *measured*). All byte counts and SHA-256 digests below are *measured* (`ids.sha256_file` on `wilkcomp02`, 2026-11-27, `o4_stage_report.json`, `gv2_stage_report.json`, `hrom1_stage_report.json`, `hrom1_ani_report.json`).

### 6.1 `O4` (`10` files, `167,680,785` B *derived*)

| File | Pod directory | Bytes (*measured*) | sha256 (*measured*) |
|---|---|---:|---|
| `o4_ncbi_organelles.fna.gz` | `/stor/scratch/Wilke/if2453/bio-box-staging/O4/` | 125,375,557 | `e4e5f6ee758779efb0d668c26b80296761eae14ec521e187c2a406d3f89947bd` |
| `o4_fungal_mitogenomes.fna.gz` | `/tmp/if2453/c86/meta/O4/` (`< 64 MiB`) | 37,958,983 | `771944c742f676ec45214c46c1556d881e49a398af91f3c2d4295fc52cd97be3` |
| `genomes.json` | `/tmp/if2453/c86/meta/O4/` | 1,828,338 | `a7afbc06c26636a83d78f26847c506a3ba340eb7d773cfc11dd8a7cc41da4f22` |
| `objects.json` | `/tmp/if2453/c86/meta/O4/` | 926 | `88553056a24e6d053f2015f11c88a5bb789f897e639ac076f7e8ddcdfea6a394` |
| `objects_staging.json` | `/tmp/if2453/c86/meta/O4/` | 946 | `ea343ee33be718babe5a98291b2fd53dda581df32acddde4811e30a68ecb44c2` |
| `records.tsv` | `/tmp/if2453/c86/meta/O4/` | 2,256,430 | `49592ea55ed59fc10f69d5ce3ef81cb76d79525fe6aa9ba7caf5eeb472ce3caf` |
| `left_out.tsv` | `/tmp/if2453/c86/meta/O4/` | 47,479 | `0027cbd73d348d5e9b06cf231a4560c71d5ae0296e8d1e1108c1d216f352ee07` |
| `repeated_records.tsv` | `/tmp/if2453/c86/meta/O4/` | 124,883 | `9b7826618150080deb07578cd769cd786c7bd53b3af58cfe1d859cd457e567d9` |
| `shared_records.tsv` | `/tmp/if2453/c86/meta/O4/` | 87,203 | `cace4b83aed2d52e5a1cf5db7dd98ce096178784c23c4e48ecf0a9284b5c66b8` |
| `later_records.tsv` | `/tmp/if2453/c86/meta/O4/` | 40 | `35715d5e3e9f7b84789769b5e7130ae87a7e93e2adfa6cae8ac0222b940589fb` |

### 6.2 `GV2` (`11` files, `1,198,145,889` B *derived*)

| File | Pod directory | Bytes (*measured*) | sha256 (*measured*) |
|---|---|---:|---|
| `gv2_giant_viruses.fna.gz` | `/stor/scratch/Wilke/if2453/bio-box-staging/GV2/` | 1,018,521,042 | `4f421ebd83694646e30675d7218336c9e927578a4d368ba57938ef62966f30ec` |
| `genomes.json` | `/tmp/if2453/c86/meta/GV2/` | 27,716,436 | `e99808259568eea62bc697b4b23ef3948d780af256df3a04d508f86de018aeba` |
| `objects.json` | `/tmp/if2453/c86/meta/GV2/` | 633 | `edda8b4003ef52f840229f299fbf9a0e486bd35e654985b1fb143c0bdfaaa98d` |
| `objects_staging.json` | `/tmp/if2453/c86/meta/GV2/` | 654 | `3e372137b48517c93ffb8ec1fc6a14c0b68357dfd2929fcd4ed007d18a81de09` |
| `records.tsv` | `/tmp/if2453/c86/meta/GV2/` | 148,233,023 | `a9fba407d6c4046cbc217d3d134af3dc5b76581a75e9022b09f9878a650ec916` |
| `left_out.tsv` | `/tmp/if2453/c86/meta/GV2/` | 21,618 | `6ab87090c0407c4cc96385ca74c5b759bf35c51a7fa31cb185daea68edb5bc88` |
| `skipped_members.tsv` | `/tmp/if2453/c86/meta/GV2/` | 323,620 | `4cbdc928f05d1a3f88534d1662e515036285c6e8df329989c2e65ee80abe2a16` |
| `repeated_records.tsv` | `/tmp/if2453/c86/meta/GV2/` | 10,313 | `282ff22bceac5f572fe99470992954b50f3c1acaa5c487c78bfed884f7806ca1` |
| `shared_records.tsv` | `/tmp/if2453/c86/meta/GV2/` | 3,318,261 | `2c81df22bee17bf9add57de56e2d0c4c89b171b53a643d2d461d623b51ac535b` |
| `later_records.tsv` | `/tmp/if2453/c86/meta/GV2/` | 40 | `35715d5e3e9f7b84789769b5e7130ae87a7e93e2adfa6cae8ac0222b940589fb` |
| `trims.tsv` | `/tmp/if2453/c86/meta/GV2/` | 249 | `73b6bf883ccadb549adc5aeed6e1add03175ea176c684753cab0101f06f91db0` |

### 6.3 `HROM1` (`14` files, `1,879,907,801` B *derived*)

| File | Pod directory | Bytes (*measured*) | sha256 (*measured*) |
|---|---|---:|---|
| `hrom1_oral_glmc.fna.gz` | `/stor/scratch/Wilke/if2453/bio-box-staging/HROM1/` | 1,829,376,298 | `aee8c620b668650b14be5d070806217c7ef8690bb2cbf1ce88c29c62125df1c7` |
| `hrom1_code4_actino1.fna.gz` | `/tmp/if2453/c86/meta/HROM1/` (`< 64 MiB`) | 13,835,044 | `3ba7708eb4930124dba6a87bad5da095214269ee7a5fd61d13ab84587ce8f495` |
| `genomes.json` | `/tmp/if2453/c86/meta/HROM1/` | 10,040,723 | `7b9752c616afd7298a187e0889e5e21bb2e5bf2780d641f7af1bef34bd849e3e` |
| `objects.json` | `/tmp/if2453/c86/meta/HROM1/` | 969 | `88f2c5e8da978a8f564fbb709a883b66e0dfd77d0dda52ad89cb8fc3aab21c31` |
| `objects_staging.json` | `/tmp/if2453/c86/meta/HROM1/` | 989 | `1ab8c5496fb0bb5ca1238dd6502e9a57b901fc2a88827f48f397c79d472d399e` |
| `records.tsv` | `/tmp/if2453/c86/meta/HROM1/` | 24,812,288 | `c9f0595cc2b6e23a635b956e504631d57555bed70f3816157d093b896dd54bd8` |
| `genetic_codes.tsv` | `/tmp/if2453/c86/meta/HROM1/` | 839,744 | `fad9dc9eb0acbee6ad03bbcae0e8b61f93480a8a59db522d00a04df40667924c` |
| `r232_status.tsv` | `/tmp/if2453/c86/meta/HROM1/` | 733,255 | `0c071bc626a01cb52961a1d359339a7b8b03e77dfbf646d9bc5f1cff10337ef8` |
| `left_out.tsv` | `/tmp/if2453/c86/meta/HROM1/` | 166,445 | `15a6548b24631ffb660c221573264dfec422d26fd72d834cae36777046bca006` |
| `shared_records.tsv` | `/tmp/if2453/c86/meta/HROM1/` | 98,696 | `d11247117a414d700316a0fa11c9ee8c4c08371b61324f563723a466f86cff58` |
| `repeated_records.tsv` | `/tmp/if2453/c86/meta/HROM1/` | 102 | `12c2d4698d8b2cac8abcbf250c1700831615e112503ed428ba6841022c944758` |
| `later_records.tsv` | `/tmp/if2453/c86/meta/HROM1/` | 55 | `ba2f1dbd71ed4260b0ec24d03d968ba0c6251afbeeb9efb134735ed1064d05f4` |
| `hrom1_ani_hits.tsv` | `/tmp/if2453/c86/meta/HROM1/` | 2,086 | `b0a2a43910f62e8629e537baa6ffd2a1e7ca5e2f2e63011677bc4b2aa0a4c896` |
| `ani_summary.json` | `/tmp/if2453/c86/meta/HROM1/` | 1,107 | `95d940b7272ad90b2964609a18d5fa0f37b0723ead0b58b4877e063583c517cc` |

### 6.4 Mirror to `tower` (`DEC-88 (1)`)

Mirrored sequentially (`scripts/orchestrator/relay_to_tower.sh --dry-run --streams 1` then `--streams 1`, one directory at a time, 2026-11-27 `18:50:39` to `18:53:54 CST`, `3,245,734,475` B = `3.246 GB` total *derived* across all `35` staged files) to `tower:/mnt/user/research-bulk/bio-box/sources/{O4,GV2,HROM1}/`:
- `sources/O4/SHA256SUMS`: `10` verified lines (`1` `/stor` FASTA `125,375,557` B at `25.08 MB/s` + `9` `/tmp` FASTA/metadata files `42,305,228` B, `0` failed; *measured*, `tower`, 2026-11-27 `18:50:51 CST`).
- `sources/GV2/SHA256SUMS`: `11` verified lines (`1` `/stor` FASTA `1,018,521,042` B at `20.37 MB/s` + `10` `/tmp` metadata files `179,624,847` B, `0` failed; *measured*, `tower`, 2026-11-27 `18:52:08 CST`).
- `sources/HROM1/SHA256SUMS`: `14` verified lines (`1` `/stor` FASTA `1,829,376,298` B at `21.03 MB/s` + `13` `/tmp` FASTA/metadata/ANI files `50,531,503` B, `0` failed; *measured*, `tower`, 2026-11-27 `18:53:54 CST`).

---

## 7. Summary across all three cohorts

| Stage | Genomes | Distinct contigs | Distinct bp | Net-new bp (after 56 stages) | Depth & mode | r232 status | Left out / skipped |
|---|---:|---:|---:|---:|---|---|---|
| `O4` | 10,944 | 10,944 | 544,499,955 | 395,960,372 | `annotated` (`pfam 38.2 gathering`), codes `1,2,3,4,11,16,22,23` | n/a (eukaryote organelles) | `431` (`411` RNA/mRNA, `18` `transl_table=32`, `2` synthetic) + `2,083` exact repeats |
| `GV2` | 10,284 | 689,523 | 3,155,057,728 | 2,267,601,494 | `called`, `gv-meta` | n/a (giant viruses) | `170` *Chlamydomonas* ICE-L nuclear scaffolds (`C-60`), `2,512` non-DNA/Prokka archive members, `111` exact repeats |
| `HROM1` | 3,732 | 266,916 | 5,883,730,291 | 5,653,991,918 | `called`, `default` codes `11,4,25` | `35` `not_in_r232`, `3,698` `not_visible` (`41` `skani >=96%` hits recorded) | `1,687` `MQ` oral reps, `1` exact repeat |
| **Total** (*derived*) | **24,960** | **967,383** | **9,583,287,974** | **8,317,553,784** | | | |

---

## 8. Pod compute & verification gates

All pod work ran at `nice -n 19` inside `/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif` (`MemoryMax=32G`, `MemorySwapMax=0`, `CPUQuota=400%`), checking `/tmp/if2453/r10/root/STOP` and `systemctl --user list-units` (`<= 1` other staging unit active, `r26-11-4-wave-6` on `MGE1`/`SAT1`, never `MV1`/`V3` `genecall`) before each launch:
- `c86-probe` + `check_gks_and_hrom` + `check_gks_pfam`: `< 5` min wall total (`CPUQuota=200%`–`400%`, *measured*).
- `c86-stage` (`O4` + `GV2`): `45` min wall (`16:38` to `17:23 CST`), peak RSS `2.8 GiB`, `< 50` CPU-minutes (*measured* via `systemctl --user status`).
- `c86-part2` (`GV2` code-25 sample + `HROM1` staging + `HROM1` `skani` `4`-thread ANI check): `70` min wall (`17:36` to `18:46 CST`), peak RSS `3.37 GiB` (`3,617,931,264` B, *measured* via `MemoryPeak`), `~175` CPU-minutes (*derived* from `CPUUsageNSec=9.62e12` ns at `18:24 CST` + `506.6` CPU-s of `skani`), well below the `300` CPU-minute per-process limit.
- `python3 -I scripts/stress/build_stage_pod/stage_est.py conf` on `O4.conf`, `GV2.conf`, and `HROM1.conf`: all three exit `0` (*measured*).
- `pod_pytest.sh C-86 tests/r1/test_c86_o4_gv2_hrom1_inputs.py tests/r1/test_r232_index.py tests/r0`: failed first (`FileNotFoundError` before `c86_o4_gv2_hrom1_inputs.py` existed), then passed 100% under both normal `PATH` and `PATH=/usr/bin:/bin` (*measured*).
