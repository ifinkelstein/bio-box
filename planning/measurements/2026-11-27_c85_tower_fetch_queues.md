# Gentle next-wave `tower` fetch queues (`23..27`), `GlobDB r232` stream slicer, and 10-minute PI portal bundle (Card `C-85`, `DEC-88 (1), (2)`)

- **Date:** `2026-11-27` (`15:28`–`16:15 CST`, `21:28`–`22:15 UTC`)
- **Card / Authority:** `C-85` (`DEC-88 (1), (2)`, `DEC-80`, `DEC-79`, `DEC-85 (4)`, `AGENTS.md` §6 & §8; following cards `C-74`, `C-82`, `R-24`, and `C-84`)
- **Scripts & Tests:** [`scripts/orchestrator/c85_prepare_tower_queues.py`](../../scripts/orchestrator/c85_prepare_tower_queues.py), [`scripts/orchestrator/after_retry20_gentle_chain.sh`](../../scripts/orchestrator/after_retry20_gentle_chain.sh), [`tests/r1/test_c85_prepare_tower_queues.py`](../../tests/r1/test_c85_prepare_tower_queues.py)
- **Host / Safety Compliance:** Zero `fetch_to_tower.sh` downloads were started on `tower` (`chain_r22.sh` PID `1901150` -> `after_chain_retry20.sh` PID `2106634` remains untouched holding `.fetch_to_tower.lock`). All metadata resolutions ran at `1` thread (`OMP_NUM_THREADS=1`, `CPUQuota=150%`, `MemoryMax=16G`, `nice -n 19`) on `wilkcomp02`; temporary HTTP `HEAD` header files on `tower` (`/tmp/c85_*.txt`) were removed immediately after inspection. Every number below carries `*measured*`, `*derived*`, `*published*`, or `*estimate*` (`AGENTS.md` §6).

---

## 1. Summary of Prepared `tower` Job Lists (`23..27`)

All five lists (`41,134` jobs, `8,546,359` bytes = `8.15 MiB` total TSV size in git, *measured* and *derived*, `2026-11-27`) pass `fetch_to_tower.sh --check --sources sources.yaml` individually and jointly alongside all `23` earlier repository lists (`128,492` total jobs across all `28` manifests in `scripts/orchestrator/tower_jobs/*.tsv`: `87,358` earlier + `41,134` new, `0` refusals, `0` duplicate destinations; *measured*, `2026-11-27`):

| List (`scripts/orchestrator/tower_jobs/`) | Target Dir (`sources_incoming/`) | Jobs (*measured*) | TSV Bytes (*measured*) | TSV `sha256` (*measured*) | Sequence Volume | Gentle Flags (`DEC-88 (1)`) |
|---|---|---:|---:|---|---|---|
| `23_r26_11_6_asg1_fna.tsv` | `ASG1_FNA/assemblies` | `403` | `76,611` | `964e8fddc0ec63701409f9895e0c92cbf6da7a0b18aeb8ed53a5ba020f47daf2` | `1,262,845,012` bp (`1.263 Gbp`, *derived*) | `--jobs 2 --pause 1` |
| `24_r26_11_6_dsv70_prjna1244896.tsv` | `DSV70_NCBI/assemblies` | `2,937` | `561,545` | `e8cfcf7c69e7c0974a406358c9e7d0b9a1d39a8ffb23e8b4c0a59e469695b839` | `6,648,188,839` bp (`6.648 Gbp`, *derived*) | `--jobs 2 --pause 1` |
| `25_r26_11_6_mgnify_nonhuman_reps.tsv` | `MGNIFY_REPS/<biome>/<prefix11>` | `22,858` | `4,713,742` | `950f599dfe194df50586274856756f26dd2e33ae12ec372d393a3389e65f5ff3` | `57,485,187,511` bp (`57.485 Gbp`, *derived*) | `--jobs 2 --pause 1` |
| `26_r26_11_6_globdb_r232.tsv` | `GLOBDB_R232` | `1` | `753` | `2ea2ebca86d634bcb1b5afe24a0b26d9441c76a889421ef065c1dc09cbd0b14f` | `292,068,313,974` B (`292.07 GB` *measured* tarball; `110.31 Gbp` in `46,920` new-resource species + `154.61 Gbp` in `47,294` GTDB R232 draft MAG reps, *derived*) | `--jobs 1 --pause 1` |
| `27_r26_11_6_gtdb_r232_delta.tsv` | `GTDB_R232_DELTA/{gtdb_r232_priority,post_r232_complete}/<2d>` | `14,935` | `3,193,708` | `8b30fcc8fb826fc45003a9773dd7e3f6f6f5448ddb55b2ee6b8414da07da1678` | `62,555,870,414` bp (`62.556 Gbp`, *derived*: `36.18 Gbp` + `26.38 Gbp`) | `--jobs 2 --pause 1` |
| **Total (`23..27`)** | — | **`41,134`** (*derived*) | **`8,546,359`** (`8.15 MiB`, *derived*) | — | **`238.26 Gbp` targeted novel species + `154.61 Gbp` GTDB R232 draft MAGs in list `26`** (*derived*) | `ionice -c 3 nice -n 19` |

---

## 2. Per-Queue Resolution & Measured Accounting

### 2.1 List `23` — `ASG1` (*Asgardarchaeota* MAGs of Appler, De Anda & Baker 2026, Card `C-82`)

- **Input inspected:** `Data/Fig1.zip -> 001_Fig1/08_PD/00_ThisStudy.tsv` inside `/stor/scratch/Wilke/if2453/bio-box-staging/ASG1/dl/Data.zip` (`404` study bin IDs, `465` `01_PreviousStudies.tsv`, `309` `02_Outgroup.tsv`; *measured*, `wilkcomp02`, `2026-11-27`).
- **Resolution against NCBI GenBank `assembly_summary.txt` (`/tmp/if2453/c58/tables/ncbi_genbank_archaea_assembly_summary.txt.gz`):**
  - Stripping the clade prefix (`Bald__`, `Gefi__`, `HeGe__`, `HeHe__`, `HeHo__`, `Hela__`, `Loki__`, `Thor__`, `Jord__`, `Rana__`, etc.) and NCBI's `_clean` suffix on `29` decontaminated isolates (`D4993_C5_H4_Bin_142_clean`, `D4994_C39_H2_Bin_6_clean`, ...) maps **`403` of `404` study bins (`99.75%`, *derived*, exceeding the `90%` gate of card `C-82`)** to `latest` public GenBank `GCA_` assemblies across three BioProjects (*measured*, `wilkcomp02`, `2026-11-27`):
    - `PRJNA743900`: `140` MAGs
    - `PRJNA1112871`: `139` MAGs (`110` exact + `29` `*_clean` isolates, released `2026-01-06` by University of Texas at Austin)
    - `PRJNA692327`: `124` unique canonical MAGs (`14` *Helarchaeia*/*Promethearchaeati* MAGs had both `2026-01-06` and `2026-06-04` taxonomy-updated GenBank accessions of identical byte/base length; `map_asg1_study_bins` deterministically selects the newer `2026-06-04` accession).
  - **Unmapped (`1` of `404` = `0.25%`, *derived*):** `HeHo__Meg22_1416_Bin_92` (not present in `ncbi_genbank_{archaea,bacteria}_assembly_summary.txt.gz`, recorded in `/tmp/if2453/c85/aux/asg1_unmapped.txt`).
  - **Span & `r232` overlap:** **`1,262,845,012` bp (`1.263 Gbp`, *derived* sum of `genome_size`)**, matching card `C-70`'s `1.0–1.4 Gbp` *estimate*. **`0` of `403`** are in `GTDB R232` and **`0` of `403`** are in `GlobDB r232` (*measured* against `/tmp/if2453/c74/index/`).

### 2.2 List `24` — `DSV70` (`PRJNA1244896` Hydrothermal Vent MAGs, Card `R-24`)

- **GTDB-Tk tree census (`/tmp/if2453/c85/DSV70/`):** `DSV70_archaea_gtdbtk.unrooted.tree` holds **`1,359` MAG tips** and `DSV70_bacteria_gtdbtk.unrooted.tree` holds **`6,063` MAG tips**, summing to **`7,422` total MAG tips** (`0` `RS_GCF_`/`GB_GCA_` reference tips; *measured* and *derived*, `wilkcomp02`, `2026-11-27`).
- **GenBank `PRJNA1244896` resolution:**
  - `ncbi_genbank_{archaea,bacteria}_assembly_summary.txt.gz` contains **`2,937` `latest` public `GCA_` MAG assemblies** under `PRJNA1244896` (`552` archaea + `2,385` bacteria, released `2025-12-12` to `2026-04-24` by Portland State University, PubMed `41540059`; *measured*, `wilkcomp02`, `2026-11-27`).
  - **`2,937` of `2,937` (`100%`, *derived*)** match the `isolate` tip labels of the two `DSV70` GTDB-Tk trees, spanning **`6,648,188,839` bp (`6.648 Gbp`, *derived*)**, with **`0` in `GTDB R232`** and **`0` in `GlobDB r232`** (*measured*).
  - **Note on the remaining `4,485` tree-only bins (`7,422 - 2,937 = 4,485`, `60.43%` *derived*):** Only `2,937` of the `7,422` study bins (`39.57%`, *derived*) were accepted/released as public `GCA_` assemblies in `PRJNA1244896` (the remainder either failed NCBI's prokaryotic MAG completeness/contamination portal screen or were left solely in EBI BioStudies `S-BSST2227`). Meanwhile, all `70` full deep-sea hydrothermal vent metagenome assemblies (`DSV70_assemblies.tar.gz`, `8,564,570,028` B, `gzip -t` clean, `28.16 GB` uncompressed) are **already intact on `tower`** (`2026-11-27_dsv70_archive_check.md` §2).
  - **Optional salvage command if the `4,485` unsubmitted bins are ever desired outside `PRJNA1244896`:** Download `DSV70_mag.tar.gz` once outside `fetch_to_tower.sh` so the mid-stream `FAIL_GZIP` check does not unlink it, and recover every intact `ustar` member preceding the corrupt deflate block:
    ```bash
    ionice -c 3 nice -n 19 curl -fSL -A curl --proto '=https' --proto-redir '=https' \
      -o /mnt/user/research-bulk/bio-box/sources_incoming/DSV70/DSV70_mag.corrupt_raw.tar.gz \
      https://ftp.ebi.ac.uk/biostudies/fire/S-BSST/227/S-BSST2227/Files/DSV70_mag.tar.gz
    # On pod/tower: stream members with tar -zxf (retaining only bins whose size matches the ustar header)
    ```

### 2.3 List `25` — `MGnify` Non-Human Prokaryotic Species Representatives

- **HTTP endpoint verification on `tower` (`2026-11-27`, `curl -sSIL`):**
  - Across MGnify catalogue versions, `species_catalogue/.../genome/{id}.fna.gz` returns `HTTP 200` on `marine/v3.0` (`504,152` B for `MGYG000444050`) but **`HTTP 404`** on `soil/v1.0` and `chicken-gut/v1.0.1` (which store uncompressed `.fna`), and `all_genomes/.../{id}.fa.gz` returns **`HTTP 404`** (*measured*, `tower`, `2026-11-27`).
  - By contrast, the exact `FTP_download` path recorded in every row of `genomes-all_metadata.tsv` (`https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/{biome}/{version}/all_genomes/{prefix}/{rep}/genomes1/{rep}.gff.gz`) returns **`HTTP 200 OK`** (`Content-Type: application/x-gzip`) across **both `v1.0` (`1,066,684` B vs `2,976,849` B uncompressed `.fna` on `MGYG000490800`) and `v3.0` (`680,727` B on `MGYG000444050`)**, bundling GFF3 CDS features plus the full `##FASTA` nucleotide contigs at `2.8×` compression (*measured* and *derived*).
- **Filtering & deduplication across the `10` non-human catalogues in `MGNIFY2` (`/tmp/if2453/c85/MGNIFY2/`, `106,365,785` B *measured*):**
  - Scanned `263,808` total genome rows -> **`52,342` species representatives** (`Genome == Species_rep`; *measured*).
  - Excluded **`3,277`** already in `GlobDB r226/r232` by `MGYG` ID, **`5,154`** already in `GTDB R226/R232` or `GlobDB` by `GCA_` / ENA WGS master prefix (`CAJVWK01`, etc.), and **`21,053`** already in `GTDB R226/R232` or `GlobDB` by 1-to-1 per-MAG derived ENA/NCBI `BioSample` (`SAMEA...`/`SAMN...`/`ERS...`; verified `50,418` of `50,423` non-empty `Sample_accession` values are unique 1-to-1 MAG BioSamples; the `5` shared metagenome samples are never excluded by BioSample) (*measured*, `wilkcomp02`, `2026-11-27`).
  - **Retained novel non-human species representatives (`22,858` jobs, `57,485,187,511` bp = `57.485 Gbp`, *measured* and *derived*, within `-4.2%` of `C-84`'s `60 Gbp` *estimate*):**
    - `marine_v3.0`: `16,542` reps (`42,292,671,690` bp = `42.29 Gbp`)
    - `cow-rumen_v1.0.1`: `2,470` reps (`5,717,544,253` bp = `5.72 Gbp`)
    - `soil_v1.0`: `1,724` reps (`4,865,403,424` bp = `4.87 Gbp`)
    - `pig-gut_v1.0`: `1,205` reps (`2,523,741,376` bp = `2.52 Gbp`)
    - `chicken-gut_v1.0.1`: `720` reps (`1,552,507,201` bp = `1.55 Gbp`)
    - `zebrafish-fecal_v1.0`: `78` reps (`312,887,456` bp = `0.31 Gbp`)
    - `sheep-rumen_v1.0`: `61` reps (`86,808,845` bp = `0.09 Gbp`)
    - `maize-rhizosphere_v1.0`: `33` reps (`92,723,670` bp = `0.09 Gbp`)
    - `non-model-fish-gut_v2.0`: `18` reps (`32,232,683` bp = `0.03 Gbp`)
    - `honeybee-gut_v1.0.1`: `7` reps (`8,666,913` bp = `0.01 Gbp`)

### 2.4 List `26` — `GlobDB r232` Genome Archive & Bounded-Memory Stream Slicer

- **Measured archive headers (`tower`, `2026-11-27`, `curl -sSIL`):** `HTTP/1.1 200 OK`, exact **`Content-Length: 292,068,313,974` bytes** (`292.07 GB` decimal / `272.01 GiB` binary, *measured* — refining `C-72`'s `299.1 GB` *derived* value where the Apache directory index printed bytes under a `KB` column header), `Last-Modified: Mon, 22 Jun 2026 11:52:25 GMT`, provider MD5 **`c95ef3b3278b9f035e2089bd729b0b48`** (*measured* from `globdb_r232_md5sum.txt`).
- **Target subset IDs on `wilkcomp02` (`/tmp/if2453/c85/probe/globdb_r232_new12_46920_ids.txt`, `sha256 9b9e3c184d8e8cc1a9d6717d44354e87ee13b23853a12f09351646c69860d166`):**
  - **`46,920` non-GTDB species representatives** across the `12` resources new in `GlobDB r232` (`GCMETA` `24,063`, `NGDC` `13,287`, `TPMCS` `3,064`, `QXLSG` `1,635`, `CRBC` `1,572`, `HOGU` `1,398`, `TPLM` `1,117`, `TG2G` `277`, `PREC` `261`, `DAWW` `147`, `CRLG` `82`, `SCSSF` `17`; `1,403` Archaea + `45,517` Bacteria), totaling **`110,310,705,273` bp (`110.311 Gbp`, *derived* from CheckM2 `Genome_Size`, within `-6.5%` of `C-72`'s `118 Gbp` *estimate*)**. (In addition, `12,087` non-GTDB species in the `13` continuing resources carry IDs absent from `r226`, *measured*.)
- **Streaming extraction command (zero full unpack to pod `/tmp`, bounded memory via `tarfile` `r|gz`):**
  ```bash
  python3 scripts/orchestrator/c85_prepare_tower_queues.py extract-globdb-subset \
    --tar-gz /mnt/user/research-bulk/bio-box/sources_incoming/GLOBDB_R232/globdb_r232_genome_fasta.tar.gz \
    --ids /tmp/if2453/c85/probe/globdb_r232_new12_46920_ids.txt \
    --out-dir /stor/scratch/Wilke/if2453/bio-box-staging/GLOBDB232/extracted
  ```

### 2.5 List `27` — `GTDB R232` Priority New Species & Post-`R232` `REFSEQ-NEW` Delta

- **Exact reconciliation of `GTDB R232` vs `GlobDB r232` (`C-74` Parquet index):**
  - `GTDB R232` has `199,923` species clusters vs `143,614` in `GTDB R226` (`+56,309` net cluster count; **`57,109`** R232 species representatives have both a species name and representative accession absent from R226 due to `800` split/renamed clusters; *measured*, `wilkcomp02`, `2026-11-27`).
  - Card `C-84` estimated `~27,854` GTDB-unique new species by subtracting `75,375 - 46,920 = 28,455` from `56,309`. Direct inspection of `/tmp/if2453/c74/index/globdb_r232.parquet` shows two key facts:
    1. **Every single one of the `199,923` GTDB R232 representatives (`100%`) is bundled inside `GlobDB r232` (`is_gtdb == True`)**, so list `26` (`globdb_r232_genome_fasta.tar.gz`) already delivers all `57,109` new GTDB R232 species representatives to `tower` (`CATALOG_GOAL.md` line 202).
    2. Meanwhile, GlobDB's non-GTDB dictionaries only carry `GCA_` accessions for `4` of `25` non-GTDB resources (`C-74` note §3.2), leaving `56,744` new R232 species reps (`56,742` resolvable on NCBI; `2` RefSeq accessions `GCF_048585685`, `GCF_048627035` suppressed/absent in Nov 2026 `assembly_summary.txt`) after subtracting `48` already built (`SYM1`) and `317` in `GTDB R226` non-reps / `GlobDB` INSDC dictionaries (*measured*).
  - **Clean priority split (`14,935` jobs <= `15,000` lines, `3.19 MB` in git, avoiding `47,294` duplicate HTTP fetches of draft MAGs already inside list `26`):**
    1. **`GTDB_R232_DELTA/gtdb_r232_priority/<2d>` (`9,448` jobs, `36,175,339,997` bp = `36.175 Gbp`, *measured* and *derived*):** Every resolvable new GTDB R232 species representative that is **Archaea** (`3,168` total before dedup), a **Complete Genome or Chromosome** (`3,454` total), a **cultured / non-MAG isolate** (`4,139` bacterial + `94` archaeal), or a **single-cell (`SAG`)** (`72`). Fetching these directly from NCBI provides clean INSDC `.fna.gz` assemblies immediately ahead of the `292.07 GB` `GlobDB r232` stream.
    2. **`GTDB_R232_DELTA/post_r232_complete/<2d>` (`5,487` jobs, `26,380,530,417` bp = `26.381 Gbp`, *measured* and *derived*):** Filtering `SYM1_post_r226_left_for_later.tsv.gz` (`49,640` rows from `C-58`) against all `901,341` `GTDB R232` genomes (`19,224` absorbed) and `GlobDB r232` leaves `30,416` genuine **post-R232** NCBI species candidates (`0` in `GTDB R232`, `0` in `GlobDB r232`), of which all **`5,487` Complete Genomes (`5,230`) and Chromosomes (`257`)** (`REFSEQ-NEW`, `DEC-85 (4)`) are queued here.
    3. **Recorded auxiliary lists on `wilkcomp02` (`/tmp/if2453/c85/aux/`):**
       - `gtdb_r232_bacterial_draft_mag_reps_in_globdb_tar.txt` (`47,294` IDs, `154,613,205,012` bp, `sha256 3ab416dc0ddb489a9c30eb821d023c2be7b97299c6bfb5d45869ff39d75cde02`): short-read bacterial Contig/Scaffold MAG species reps of GTDB R232 already inside `26_r26_11_6_globdb_r232.tsv` — extractable via `extract-globdb-subset` with zero extra network traffic.
       - `post_r232_ncbi_draft_species_reps.txt` (`24,929` IDs, `102,386,061,281` bp, `sha256 7ef651f2c34d7a0a13d3c40ab8d67c2e56be4107eefbc5613ff430a4a0483987`): post-R232 draft Scaffold/Contig GenBank candidates deferred behind complete genomes.

---

## 3. Gentle `tower` Chain Launch Command (`DEC-88 (1)`)

**Not launched in this card** (`after_chain_retry20.sh` PID `2106634` is currently queued behind `chain_r22.sh` PID `1901150` on `tower`). After merge and checkout refresh on `tower`, the orchestrating session launches the gentle chain with:

```bash
# 1. Dry-run verification on tower first (0 network requests, 0 writes):
/bin/bash scripts/orchestrator/after_retry20_gentle_chain.sh --dry-run 2106634

# 2. Detached launch behind after_chain_retry20.sh (PID 2106634):
nohup /bin/bash scripts/orchestrator/after_retry20_gentle_chain.sh 2106634 \
  > /mnt/user/research-bulk/bio-box/fetch_logs/2026-11-27_gentle_c85_chain.log 2>&1 &
```

Guarantees enforced by `after_retry20_gentle_chain.sh`:
- Waits cleanly for PID `2106634` (`after_chain_retry20.sh`) to exit before touching `.fetch_to_tower.lock`.
- Runs every `fetch_to_tower.sh` process under **`ionice -c 3`** (idle disk I/O class) and **`nice -n 19`** (lowest CPU scheduling priority).
- Caps multi-file queues (`23, 24, 25, 27`) at **`--jobs 2 --pause 1`** (`<= 2` concurrent `curl` streams, `1 s` pause per worker after each file) and single large tarballs (`26_r26_11_6_globdb_r232.tsv`) at **`--jobs 1 --pause 1`**.
- Checks `/mnt/user/research-bulk/bio-box/sources_incoming/STOP` before waiting and before every list (`touch .../STOP` halts the chain cleanly before the next file/list).

---

## 4. Validation Summary

- `fetch_to_tower.sh --check --sources sources.yaml` on each individual list (`23..27`):
  - `23_r26_11_6_asg1_fna.tsv`: `check: 403 URLs, all covered by sources.yaml; nothing fetched` (exit `0`)
  - `24_r26_11_6_dsv70_prjna1244896.tsv`: `check: 2937 URLs, all covered by sources.yaml; nothing fetched` (exit `0`)
  - `25_r26_11_6_mgnify_nonhuman_reps.tsv`: `check: 22858 URLs, all covered by sources.yaml; nothing fetched` (exit `0`)
  - `26_r26_11_6_globdb_r232.tsv`: `check: 1 URLs, all covered by sources.yaml; nothing fetched` (exit `0`)
  - `27_r26_11_6_gtdb_r232_delta.tsv`: `check: 14935 URLs, all covered by sources.yaml; nothing fetched` (exit `0`)
- `fetch_to_tower.sh --check --sources sources.yaml scripts/orchestrator/tower_jobs/*.tsv` across all `28` repository manifests (`01..27` + `2026-11-26_c72_part1.tsv`): `check: 128492 URLs, all covered by sources.yaml; nothing fetched` (exit `0`, zero destination collisions).
- `bash -n scripts/orchestrator/after_retry20_gentle_chain.sh`: exit `0`.
- `pod_pytest.sh C-85 tests/r1/test_c85_prepare_tower_queues.py tests/r1/test_r08_fetch_to_tower.py tests/r1/test_r15_fetch_ftp.py tests/r0`: `386 passed in 43.00s` (*measured*, `wilkcomp02`, `2026-11-27`, including both normal `PATH` and tools-hidden `PATH`).

---

## 5. Consolidated 10-Minute PI Browser & Portal Bundle (`DEC-88 (2)`, `DEC-50`)

All high-novelty datasets gated on interactive browser SSO, signed portal URLs, or Cloudflare Turnstile (`AGENTS.md` behavioral rules: never bypassed by agents) are consolidated into **three quick browser actions (~10 minutes total)**:

### Step 1 (~3 min) — JGI SSO Session Cookie for `532` Unrestricted Long-Read Metagenomes (`257.48 Gbp` *measured*)

- **Target manifest already prepared on `wilkcomp02`:** `/tmp/if2453/s2/inputs/LR1/LR_jgi_batch.tsv.gz` (`18,846` bytes, `sha256 02df650bdb58be67d0a0a8e8cd1d5a6b8fc206eaa1915619d3f05d442dfce32a`, `532` unrestricted JGI-only `IMG/M` PacBio/ONT metagenome assemblies = **`257.48 Gbp`** *measured* in card `C-56`, `redistributable: false` under JGI Data Utilization Policy).
- **PI browser steps:**
  1. Sign in at `https://data.jgi.doe.gov/` (or `https://genome.jgi.doe.gov/portal/`) and click **Accept** on the JGI Data Usage Policy prompt.
  2. Click Profile icon -> **Copy My Session Token** (or export Netscape `cookies.txt`) and save on `tower` at `/root/.jgi_token` (`chmod 600 /root/.jgi_token`, never committed to git).
  3. Optional while logged into `https://img.jgi.doe.gov/m`: `Find Genomes -> Metagenome Bins -> Public=Yes, Unrestricted, HQ/MQ, Add Date > 2024-02-06` -> export table/cart IDs (`~10k–20k` net-new species, `~25–50 Gbp` *estimate*).

### Step 2 (~5 min) — Biosino `NODE` Signed-URL / FTPS Harvest (`MEER` Hadal Trench, `CPR1`, `GDSV 2.0`)

- **Why:** Anonymous CNGB FTP (`CNP0004890`) exposes raw reads only; the assembled genomes/viromes sit on Biosino `NODE` (`https://www.biosino.org/node/`).
- **PI browser steps:**
  1. Log into `https://www.biosino.org/node/` and open three project pages:
     - **`OEP004067` (`MEER` Mariana Hadal Trench, Xiao *et al.* *Cell* 2025):** `https://www.biosino.org/node/project/detail/OEP004067` — **`7,564` hadal prokaryotic SRGs (`89.4%` novel) + `23,095` trench viral contigs** (*published*, `~20 Gbp` *estimate*).
     - **`OEZ00021651` / `OED842740` (`CPR1` `PBGD` + `PPGD` *Patescibacteria* & patesciphage compendium):** `https://www.biosino.org/node/project/detail/OEZ00021651` — **`24,523` CPR genomes (`13,242` HQ) + `18,876` patesciphages** (*published*).
     - **`OEP003844` (`GDSV 2.0` Global Deep-Sea Sediment Virome):** `https://www.biosino.org/node/project/detail/OEP003844` — **`554,664` benthic DNA vOTUs** (*published*).
  2. Either copy the browser's signed HTTPS download links into `/mnt/user/research-bulk/bio-box/sources_incoming/BIOSINO_PORTAL/urls.txt`, or run `lftp` gently on `tower`:
     ```bash
     mkdir -p /mnt/user/research-bulk/bio-box/sources_incoming/BIOSINO_PORTAL
     ionice -c 3 nice -n 19 lftp -e \
       "set ftp:ssl-force true; open ftps://fms.biosino.org:2122; user <EMAIL> <PASSWORD>; lcd /mnt/user/research-bulk/bio-box/sources_incoming/BIOSINO_PORTAL; mirror --parallel=2 <PROJECT_PATH>; bye"
     ```

### Step 3 (~2 min) — Dryad 1-Click Browser Download (`PEONY 1.3` + `PIGEON 3.0`, `sources.yaml`: `peony_pigeon3_dryad`)

- **Why:** Dryad API v2 `/download` requires an ORCID Bearer Token (`HTTP 401`) and browser stream URLs sit behind Cloudflare/WAF (`HTTP 403` to headless `curl`; never bypassed per `AGENTS.md`). License is **`CC0-1.0` (`redistributable: true`)**.
- **PI browser steps:**
  1. Open `https://datadryad.org/dataset/doi:10.5061/dryad.prr4xgz35` in Chrome (or sign in with ORCID -> Account -> **API Token** and place token in `tower:/root/.dryad_token`).
  2. Download the four files listed under `peony_pigeon3_dryad` in [`sources.yaml`](../../sources.yaml#L1882-L1895) into `tower:/mnt/user/research-bulk/bio-box/sources_incoming/PEONY_PIGEON/` and verify against the published Dryad API v2 SHA-256 checksums (*measured*, `2026-11-24`):
     ```bash
     cd /mnt/user/research-bulk/bio-box/sources_incoming/PEONY_PIGEON
     cat <<'EOF' > SHA256SUMS
     fb474db3a2feb9a0f83fcd760aef2e7bcc3396fc9eab8809e46abb0d687eaf05  PEONY_PIGEON3_metadata_Nov2026.csv
     0107a228f09fafff66ddedcdaa13f23e42d50058e917515115915a9c269b2928  PEONY1.3_Nov26.fa.gz
     f12f09a59bf259384c595ae463450d1a886f11bcaf348f45285b3dc17d9cf5b5  PIGEON3.0.fa.gz
     6e6178f64aa016a0eccf721977bd1ce2f8e4ef0310165249b66aafde648ab53c  README.md
     EOF
     ionice -c 3 nice -n 19 sha256sum -c SHA256SUMS
     ```
