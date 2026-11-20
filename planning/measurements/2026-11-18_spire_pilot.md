# SPIRE 1,000-sample unbinned-contig pilot: habitat composition, depth-1 costs, NR50 novelty and on-the-hit annotation

Card C-42, 2026-11-18 to 2026-11-20 UTC. Executed on `wilkcomp02` (`112` logical CPUs, local NVMe `/tmp`) inside the pinned runtime image `/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif` (`bio-box` branch `task/C-42` from `0948bc8`), strictly inside scratch `/tmp/if2453/s3/spire/` (peak scratch `172 GiB` / `184.7 GB`, cap `200 GB`; `/tmp/if2453/r10/root` and `/stor/scratch/Wilke/if2453/all-the-data/l0/refdb/pfam/38.2` bound read-only; zero bytes written to `/tmp/if2453/r10/root` or `/stor`). Every process ran under `systemd-run --user` with `MemoryMax` (`128G` for Steps 1–4, `180G` with `MemorySwapMax=0` for Steps 5–6) at `nice 19`, within the hour's thread allowance (`1` thread during concurrent builds, `28–56` threads when all other login users used `0.00` CPUs), and every single process stayed strictly under `280` CPU-minutes (`AGENTS.md` §7.2).

**Labels (`AGENTS.md` §6).** Every number below carries one of four labels:
- *measured*: counted directly on `wilkcomp02` (`2026-11-18` to `2026-11-20` UTC) by the commands in §8, with the exact output file under `planning/measurements/2026-11-18_spire_pilot/` named beside it;
- *derived*: explicit arithmetic on named *measured* or *published* inputs;
- *published*: stated by an external paper or web resource, with citation;
- *estimate*: an engineering judgement or power-law extrapolation, with its basis stated.

---

## 0. Up-front `AGENTS.md` §6 `>20%` discrepancies against `sources_metagenomic_survey.md` and `C-38_import_depths.md`

Per `AGENTS.md` §6 and card C-42's acceptance rule, every measured quantity that differs by more than `20%` from `planning/sources_metagenomic_survey.md` (§2 #10, §4.1) or `planning/design/C-38_import_depths.md` (§1 `N1..N24`) is reported here first:

| Quantity | Survey / C-38 value | C-42 value | Ratio (`C-42 / prior`) | Why it differs |
|---|---|---|---|---|
| SPIRE v1 total assembled bases (survey §2 #10) | `16 Tbp` (*published*, Schmidt et al. 2024) | **`22.771 Tbp`** *measured* across `89,888` samples with API `assembly_bp` (`23.348 Tbp` *derived* across all `96,067` samples; `select/select_summary.json`) | **`1.423×` (`+42.3%`)** known; **`1.459×` (`+45.9%`)** estimated | SPIRE's paper rounded / quoted a filtered total; `/api/v1/samples` sums raw per-sample MEGAHIT `assembly_bp` (`>= 200 bp`). |
| Pilot raw assembled bases (`N24`) | `160 Gbp` *estimate* (`16 Tbp / 99,146 × 1,000`) | **`341.501 Gbp`** *measured* across `990` non-empty assemblies (`114.419 Gbp` at `>= 1.5 kb`; `contigs/contigs_summary.json`) | **`2.134×` (`+113.4%`)** unfiltered; **`0.715×` (`-28.5%`)** at `>= 1.5 kb` | Stratified allocation $\propto \sqrt{N_h}$ deliberately upsamples deep environmental strata (`terrestrial` mean `567.5 Mbp`, `aquatic` `419.9 Mbp`) relative to `host-associated:animal host` (`188.8 Mbp`, `64.9%` of SPIRE). |
| Depth-1 (`called`) storage bytes per Gbp (survey §4.1 `0.71 GB/Gbp`; `N1=0.309`, `N2=0.405`, `N3=0.714`) | `0.71 GB/Gbp` (`0.37 GB/Gbp` without nt & IGR) (*derived* from `H1`) | Unfiltered (`SP0`): **`1.3739 GB/Gbp`** (`0.8110` without nt/IGR) *measured* (`depth1_unfiltered/depth1_summary.json`). At `>= 1.5 kb` (`SPA1..SPU1`): **`0.8515 GB/Gbp`** (`N1 = 0.4098`, `N2 = 0.4417`; **`0.4601 GB/Gbp`** without nt/IGR) *measured* (`depth1/depth1_summary.json`) | Unfiltered: **`1.935×` (`+93.5%`)** (`2.192×` no-nt). At `>= 1.5 kb`: **`1.199×` (`+19.9%`)** total; **`1.326×` (`+32.6%`)** `N1`; **`1.244×` (`+24.4%`)** no-nt | Short contigs inflate `.fai`/`contig_index`/`genome_input` (`1.218M` contigs/Gbp unfiltered, `255.1k` contigs/Gbp at `>= 1.5 kb` vs `3,689` in `H1`), and intra-stage protein uniqueness is `93.6%` (`1.034M` distinct proteins/Gbp vs `0.884M` in `H1`). |
| Default CRISPR producer CPU per Gbp (survey §4.1 `3.0`; `N13`) | `3.0 CPU-min/Gbp` (*derived* from `E-12` on isolate/phage contigs) | **`10.245 CPU-min/Gbp`** *measured* on `SPA1`+`SPU1` (`1.485 Gbp`, `302,227` contigs; `depth1/depth1_summary.json`); `>33 CPU-min/Gbp` on `SP0` (`3.93M` contigs) | **`3.415×` (`+241.5%`)** at `>= 1.5 kb` | `src/biobox/l1/crispr.py:566` writes one 1-row Parquet row group per contig into scratch `contigs.parquet` (`302k`–`3.93M` row groups), inflating DuckDB footer scanning and RAM (§4.4). |
| Genes and distinct proteins per Gbp (`N6=2.16M`, `N7=2.19M`, `N9=98.8%` unique, `N16=323.6 aa`) | `2.16M`–`2.19M` genes/Gbp; `323.6 aa` (`P2`) | Unfiltered (`SP0`): `1.812M` genes/Gbp (`77.37%` fragments). At `>= 1.5 kb`: **`1.104M` genes/Gbp**, **`1.006M` distinct proteins/Gbp** (`115.03M` total across `114.37 Gbp`), `33.98%` fragments, mean **`258.90 aa`** *measured* | At `>= 1.5 kb`: **`0.511×` (`-48.9%`)** genes/Gbp vs `N6`; **`0.800×` (`-20.0%`)** aa length vs `N16` | Dropping `< 1.5 kb` contigs removes the `200–1,499 bp` band where Prodigal calls `1.54M–2.85M` truncated ORFs/Gbp (`76%–98%` fragments). |
| Search-pack bytes per Gbp (`N22 = 0.495 GB/Gbp`, `DEC-54 (3)`) | `0.495 GB/Gbp` (*derived* from `rM.2`) | **`1.0100 GB/Gbp`** *measured* on `SPA1` (`unique` `0.1133`, `igr` `0.0879`, `dark_orf` **`0.8087 GB/Gbp`**; `depth1/pack_probe_SPA1.json`) | **`2.040×` (`+104.0%`)** total; **`0.229×` (`-77.1%`)** for `unique` alone | Invariant 16's `dark_orf` target extracts all 6-frame ORFs in the `3 kb` end zones of every linear contig (`80.1%` of the pack on `1.5–5 kb`-dominated assemblies). |
| MMseqs2 `linclust` memory per sequence (`N20`) | `336 B/protein` (*estimate*, MMseqs2 wiki) | **`621.1 B/protein`** *derived* (`96,578.0 MiB` = `101.27 GB` peak RSS on `01_kmermatcher` over `163,060,400` proteins with `--split-memory-limit 110G`; `novelty/linclust_steps.json`) | **`1.848×` (`+84.8%`)** | Actual `kmermatcher` k-mer table + sequence index RSS before split flush is ~`0.62 KB/seq` (bounded cleanly by `--split-memory-limit`). |
| One-HMM search CPU (`N21`) | `8–11 CPU-min` per `10.6M` proteins (*derived* from `W-21` pack runner) | **`2.792 CPU-min` per `10M` proteins** *measured* (`32.113 CPU-min` across `115,030,133` proteins; `tniq/tniq_summary.json`) | **`0.27×–0.37×` (`-63%` to `-73%`)** | Direct streaming of sorted `protein_seq` Parquet into `pyhmmer` skips search-pack zstd shard copy/verification overhead. |

By contrast, the two primary depth-1 CPU rates at `>= 1.5 kb` agree tightly with `sources_metagenomic_survey.md` §4.1 and `C-38_import_depths.md`:
- **Ingest CPU (`>= 1.5 kb`):** **`4.840 CPU-min/Gbp`** *measured* vs `4.5 CPU-min/Gbp` in survey §4.1 (**`1.076×`, `+7.6%`**) and `4.11 CPU-min/Gbp` in `B-52` (`1.178×`, `+17.8%`).
- **Gene calling (`meta` mode, `>= 1.5 kb`):** **`60.662 CPU-min/Gbp`** *measured* (`6,938.11 CPU-min` over `114.374 Gbp`) vs **`60.4 CPU-min/Gbp`** in `C-38` `N12` (**`1.0043×`, `+0.43%`**) and `60–91 CPU-min/Gbp` in survey §4.1.

---

## 1. Source inspection and registration (`sources.yaml` entry `spire_pilot_assemblies`)

Inspection of SPIRE's study pages (`https://spire.embl.de/spire/v1/study/Howe_2014_soil`), API documentation (`https://spire.embl.de/api/`), and FAQ (`https://spire.embl.de/faqs`) on `2026-11-18` established:
1. **Per-sample contig FASTA URL:** `https://spire.embl.de/download_assembly/{sample_id}` (HTTP 302 redirect to `https://swifter.embl.de/~fullam/spire/...` or `https://black.embl.de/~fullam/spire/...`, serving gzipped MEGAHIT contigs `>= 200 bp`).
2. **Sample assembly/binning table URL:** `https://spire.embl.de/api/v1/samples?format=tsv&limit=10000&offset={offset}` (`10` pages of `<= 10,000` rows covering all `96,067` v1 samples, providing `assembly_bp`, `assembly_n50`, `binned_fraction`, `bases`, `reads`, `n_mags`, `n_hq_mags`, `n_mq_mags`, which were absent from C-43's study TSVs).
3. **Per-MAG FASTA URL:** `https://spire.embl.de/download_file/{mag_id}` (gzipped FASTA of individual SPIRE v1 MAGs).
4. **Licence and rate limits:** CC BY-SA 4.0 (`https://spire.embl.de/faqs`, `redistributable: true` under ShareAlike); documented API limit `120 requests/min` (`2 req/s`).

Entry `spire_pilot_assemblies` (`approved: under DEC-50`) was committed to `sources.yaml` (`afbee51`) before any sequence request (`AGENTS.md` §8). Polite fetch summary (`1` connection, `>= 0.5 s` interval, `User-Agent: curl` with zero PII, stopping on 5 consecutive HTTP refusals):
- **Metadata & API tables (`c42_fetch.py tables`):** `730` files (`151,021,997` bytes *measured*, including `27,407,083` bytes across `10` `/api/v1/samples` TSV pages), `2026-11-19T01:07:42Z` to `01:21:16Z`, `0` refusals; log `fetch_tables_log.jsonl` sha256 `155ecae3d4506ec9ef3e1dd97920913dc3f5a745b49b30317158da6202da8526` (*measured*).
- **Pilot sample assemblies (`c42_fetch.py assemblies`):** `1,000` requests (`990` HTTP 200 `.fa.gz` files = **`105,461,158,824` bytes** / `105.461 GB` wire and disk; `10` HTTP 404 where SPIRE hosts no assembly file: `SAMEA103892532`, `SAMEA9531694`, `SAMEA9531781`, `SAMEA9531906`, `SAMEA9531909`, `SAMEA9531927`, `SAMEA9531933`, `SAMN21987496`, `SAMN21987622`, `SAMN21987745`; `0` HTTP refusals), `2026-11-19T01:35:56Z` to `03:35:50Z`; log `fetch_assemblies.jsonl` sha256 `b63a1b40f0b64873eab5e3fcd1a8bc8f175d2ad391d1eef4963805afb9fc5fc6` (*measured*).
- **MAG overlap probes (`c42_fetch.py mag-probe` & `sample-mags`):** `4` MAGs for `SAMN06343983` (`mag_probe_SAMN06343983.json`) plus `7,157` MAG requests across the `339` MAG-bearing `TniQ`-hit samples (`7,029` HTTP 200 `.fa.gz` streamed and unlinked immediately = `5,339,657,208` wire bytes, `128` HTTP 404, `0` HTTP refusals), `2026-11-20T02:24:48Z` to `05:25:38Z`; log `fetch_sample_mags.jsonl` sha256 `b77860809d73f92701c2be6f0d50b7160cf4d62219d38d14e03d95ef3a6106f6` (*measured*).

---

## 2. Full SPIRE v1 habitat inventory and $\propto \sqrt{N_h}$ sample selection (`select/`)

Joining the `715` SPIRE study TSVs (`96,067` samples *measured*), `/api/v1/samples` (`96,067` rows, of which `89,888` carry non-null `assembly_bp` *measured*), and `spire_v1_genome_metadata.tsv.gz` (`1,147,979` v1 MQ/HQ MAGs totaling `2,603,289,669,758 bp` = `2.603 Tbp` *measured*) yields `6` level-1 habitats and `32` level-2 habitats (`select/habitat_level1_all.tsv`, `select/habitat_level2_all.tsv`).

**Selection rule and seed (`c42_select.py`).** Across the $H = 32$ level-2 habitats of size $N_h$ ($\sum_{h=1}^{32} N_h = 96,067$), target $K = 1,000$ samples are allocated by Hamilton's largest-remainder method proportional to $\sqrt{N_h}$ subject to $1 \le k_h \le N_h$. Within each stratum $h$, sample IDs are sorted lexicographically and sampled without replacement via an isolated `random.Random(f"20261118:{habitat2}")` PRNG (`seed = 20261118`). This preserves rare habitats (`terrestrial:tundra` $N_h = 1 \to k_h = 1$, `anthropogenic:waste` $N_h = 8 \to k_h = 3$, `air` $N_h = 17 \to k_h = 4$) while compressing `host-associated:animal host` from `64.9%` of SPIRE (`62,392` samples) down to `23.7%` of the pilot (`237` samples), spanning `320` distinct studies (*measured*, `select/select_summary.json`).

### 2.1 Full SPIRE v1 vs selected pilot by level-1 habitat

*All counts `*measured*` (`select/habitat_level1_all.tsv`, `select/pilot_by_habitat1.tsv`) except `Est. all Tbp` (`*derived*`: stratum mean `assembly_bp` applied to the `6,179` samples lacking API `assembly_bp`; `unlabelled` has `0` rows with API `assembly_bp`).*

| Level-1 habitat | SPIRE samples | With API `assembly_bp` | Known `assembly_bp` (Tbp) | Est. all (Tbp, *derived*) | Mean sample (Mbp) | Median sample (Mbp) | SPIRE v1 MAGs | MAG bp (Tbp) | MAG share of known bp (*derived*) | Weighted MetaBAT2 binned frac | Pilot $k_h$ (`1,000`) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `host-associated` | 63,909 | 62,044 | 12.370 | 12.742 | 199.4 | 117.0 | 933,389 | 2.061 | 16.66% | 31.15% | 279 |
| `aquatic` | 11,141 | 11,028 | 4.631 | 4.678 | 419.9 | 224.9 | 97,341 | 0.224 | 4.83% | 11.85% | 244 |
| `terrestrial` | 7,492 | 7,244 | 4.111 | 4.252 | 567.5 | 116.2 | 48,028 | 0.135 | 3.29% | 11.12% | 215 |
| `anthropogenic` | 9,662 | 9,555 | 1.659 | 1.678 | 173.7 | 33.3 | 54,841 | 0.152 | 9.14% | 18.19% | 199 |
| `unlabelled` | 3,846 | 0 | — | — | — | — | 14,351 | 0.032 | — | — | 59 |
| `air` | 17 | 17 | 0.00049 | 0.00049 | 29.0 | 26.1 | 29 | 0.00013 | 26.33% | 38.10% | 4 |
| **All SPIRE v1** | **96,067** | **89,888** | **22.771** | **23.348** | **253.3** | — | **1,147,979** | **2.603** | **11.43%** | — | **1,000** |

The full 32-row level-2 table is in [`2026-11-18_spire_pilot/select/habitat_level2_all.tsv`](2026-11-18_spire_pilot/select/habitat_level2_all.tsv) and the 1,000-sample manifest is in [`2026-11-18_spire_pilot/select/pilot_selection.tsv`](2026-11-18_spire_pilot/select/pilot_selection.tsv).

---

## 3. Contig length distributions, binned shares and inter-sample repeats (`contigs/`)

Streaming all `990` fetched `.fa.gz` assemblies through `c42_contigs.py` (`453,491,664` FASTA records, **`341,501,413,064 bp` = `341.501 Gbp`** *measured*, `contigs/contigs_summary.json`) established:

### 3.1 Assembled bases by contig-length cutoff and binned fraction per level-1 habitat

*All values `*measured*` (`contigs/contigs_by_habitat1.tsv`, `contigs/contigs_summary.json`; full 32-stratum breakdown in [`2026-11-18_spire_pilot/contigs/contigs_by_habitat2.tsv`](2026-11-18_spire_pilot/contigs/contigs_by_habitat2.tsv)).*

| Habitat | Non-empty / selected | Contigs (`>= 200 bp`) | Total bp (Gbp) | Share `>= 1 kb` | Share **`>= 1.5 kb`** (Gbp) | Share **`>= 5 kb`** (Gbp) | Share **`>= 10 kb`** (Gbp) | Share `>= 100 kb` | Share in SPIRE MQ/HQ MAGs | MetaBAT2 `binned_fraction` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `terrestrial` | 208 / 215 | 175,067,299 | 122.407 | 39.16% | **28.37%** (34.724) | **11.62%** (14.228) | **6.96%** (8.517) | 0.69% | **2.96%** | 11.31% |
| `aquatic` | 244 / 244 | 138,168,521 | 100.308 | 41.86% | **31.46%** (31.557) | **13.90%** (13.941) | **8.50%** (8.526) | 0.82% | **5.73%** | 14.92% |
| `host-associated` | 279 / 279 | 66,359,041 | 62.637 | 55.79% | **46.14%** (28.901) | **27.35%** (17.131) | **19.83%** (12.420) | 3.52% | **12.93%** | 28.64% |
| `anthropogenic` | 199 / 199 | 70,938,720 | 53.169 | 43.24% | **33.38%** (17.749) | **16.67%** (8.861) | **11.29%** (6.002) | 1.75% | **7.92%** | 15.61% |
| `unlabelled` | 56 / 59 | 2,829,099 | 2.864 | 59.67% | **49.88%** (1.428) | **29.49%** (0.844) | **21.63%** (0.620) | 3.83% | **15.96%** | — |
| `air` | 4 / 4 | 128,984 | 0.117 | 57.51% | **51.74%** (0.060) | **38.09%** (0.045) | **32.35%** (0.038) | 14.13% | **32.96%** | 44.56% |
| **Pilot total** | **990 / 1,000** | **453,491,664** | **341.501** | **43.82%** (`149.63 Gbp`) | **33.50%** (`114.419 Gbp`, `29.18M` contigs) | **16.12%** (`55.050 Gbp`, `4.29M` contigs) | **10.58%** (`36.123 Gbp`, `1.50M` contigs) | **1.44%** (`4.93 Gbp`) | **6.49%** (`22.180 Gbp`) | **16.51%** |

**Key findings from §3.1:**
1. **Why MAG catalogues miss `88%–97%` of environmental metagenomic sequence:** Across `terrestrial` and `aquatic` metagenomes, SPIRE's MQ/HQ MAGs capture only **`2.96%`** and **`5.73%`** of assembled bases (`11.31%` and `14.92%` even counting all raw MetaBAT2 bins before quality filtering; *measured*). Even at `>= 10 kb`, terrestrial and aquatic samples hold `8.517 Gbp` and `8.526 Gbp` (`2.35×` and `1.48×` their entire MQ/HQ MAG bases, *derived*).
2. **Verification that MAG contigs are exact `contig_uid` subsets of sample assemblies (`mag_probe_SAMN06343983.json`):** On probe sample `SAMN06343983` (`37,571` contigs, `40,534,684 bp` *measured*), all `4` published SPIRE v1 MAGs (`394` contigs, `15,902,704 bp`, minimum contig length `2,549 bp` *measured*) matched `394 / 394` exact `contig_uid` digests in the sample assembly with `0` contigs shared between MAGs (*measured*).
3. **Comparison with CAMI (`N10`):** On SPIRE's real MEGAHIT assemblies, `< 1 kb` contigs are `87.11%` of records (`395,030,171 / 453,491,664` *derived*) and **`56.18%`** of bases (`191.871 Gbp` *derived*), `1–10 kb` contigs hold **`33.24%`** of bases (`113.507 Gbp` *derived*), and `>= 10 kb` contigs hold **`10.58%`** of bases (`36.123 Gbp` *measured*), versus CAMI's simulated gold-standard assembly (`52.4%`, `27.5%`, `20.0%`; `N10`).

### 3.2 Exact `contig_uid` repeats across samples and overlap with `r26.11.1`

*All values `*measured*` (`contigs/contigs_summary.json`, `contigs/repeats_by_habitat1.tsv`) via 16-bucket spilled DuckDB aggregation over all `453,491,664` records.*

| Metric | All contigs (`>= 200 bp`) | Contigs `>= 1.5 kb` |
|---|---:|---:|
| Total FASTA records | `453,491,664` (`341.501 Gbp`) | `29,180,674` (`114.419 Gbp`) |
| Within-sample duplicate `contig_uid` records | `0` (`0.0000%`) | `0` (`0.0000%`) |
| Distinct forward `contig_uid`s | `453,429,050` (`341.430 Gbp`) | `29,174,924` (`114.374 Gbp`) |
| Contigs appearing in `>= 2` pilot samples | `59,248` (`0.0131%` of contigs) | `5,309` (`0.0182%` of contigs) |
| — within a single level-1 habitat only | `58,441` (`98.64%` of repeats) | `5,205` (`98.04%` of repeats) |
| — shared across `>= 2` level-1 habitats | `807` (`1.36%` of repeats) | **`104`** (`1.96%` of repeats) |
| Redundant bases across samples (exact `contig_uid`) | `71,584,853 bp` (**`0.0210%`** of bp) | `44,636,665 bp` (**`0.0390%`** of bp) |
| Distinct strand-canonical UIDs ($\min(\text{fwd}, \text{rc})$) | `453,428,659` (`341.430 Gbp`) | `29,174,881` (`114.374 Gbp`) |
| Redundant bases across samples (strand-canonical) | `71,858,090 bp` (`0.0210%` of bp) | `44,735,700 bp` (`0.0391%` of bp) |

**Overlap with `r26.11.1` (`/tmp/if2453/r10/root`, stages `S1, P1..P6, H1, PL1`):** Across all `453,429,050` distinct pilot contigs, **only `3` exact `contig_uid`s (`4,054 bp`)** are already stored in `r26.11.1`: `1` contig (`763 bp`) in `P3` and `2` contigs (`3,291 bp`) in `PL1` (`0` in `S1`, `P1`, `P2`, `P4`, `P5`, `P6`, `H1`; *measured*, `contigs/contigs_summary.json`). Even when the same clonal plasmid or phage is present in two metagenomes, independent MEGAHIT assemblies differ by `>= 1` terminal base, so `contig_uid` deduplication eliminates almost zero bases (`0.039%`) across samples; redundancy collapses at the **protein** level (`91.08%` distinct `protein_uid`s across the pilot, and `50%`-identity clustering in §5).

---

## 4. Building at depth 1 (`--depth called`) in scratch (`depth1_unfiltered/`, `depth1/`)

Inputs were constructed by `c42_inputs.py` using one `GenomeInput` per unbinned contig (`genome_uid = blake3("biobox:c42:unbinned_contig:v1:" + contig_uid.hex())`, `kind="prokaryote"`, `genetic_code=11`, `circular=()`, `assembly_quality="contig"`), guaranteeing Prodigal `meta` mode on every contig regardless of length (`DEC-13`, `C-38` §5.3).

### 4.1 Part A: Unfiltered calibration stage `SP0` (`--min-contig-len 0`, `32` samples across all 32 strata)

Built on `1` sample per level-2 stratum (`32` samples, `3,930,592` contigs, `3,226,018,786 bp` = `3.226 Gbp` *measured*, `depth1_unfiltered/depth1_summary.json`):
- **CPU rates (`SP0` unfiltered):** `ingest` **`5.163 CPU-min/Gbp`** (`16.657 CPU-min`), `genecall` **`68.865 CPU-min/Gbp`** (`222.159 CPU-min`), `igr` `0.406 CPU-min/Gbp`, `tapes` `0.542 CPU-min/Gbp` (*measured*). Separately running default `crispr` on `SP0` found `437` arrays (`135.46 / Gbp` *measured*), peaking at `30,244 MiB` RSS and `>105 CPU-min` in DuckDB footer scanning due to the 1-row-group-per-contig bug (§4.4).
- **Storage (`SP0` unfiltered):** depth `stored` = **`0.7568 GB/Gbp`** (`seq_shards` alone = `0.5621 GB/Gbp`, `+39.5%` above CAMI `N5` `0.403`); depth `called` total = **`1.3739 GB/Gbp`** (`+93.5%` above survey §4.1's `0.71 GB/Gbp`); without `seq_shards` and `igr` = **`0.8110 GB/Gbp`** (`+119.2%` above `0.37 GB/Gbp`; *measured*).
- **Genes and fragment share (`dist_to_contig_end_bp == 0`) across all 6 length bands (`SP0`, *measured*):**

| Contig length band | Contigs | Bases (Mbp) | Genes called | Genes / Gbp (*derived*) | Fragment genes (`dist == 0`) | Fragment share (*derived*) |
|---|---:|---:|---:|---:|---:|---:|
| `200–499 bp` | 2,271,170 | 862.94 | 2,458,336 | 2,848,800 | 2,403,413 | **97.77%** |
| `500–999 bp` | 1,191,694 | 793.67 | 1,589,074 | 2,002,191 | 1,441,913 | **90.74%** |
| `1,000–1,499 bp` | 224,263 | 269.52 | 414,571 | 1,538,194 | 316,178 | **76.27%** |
| `1,500–4,999 bp` | 195,753 | 473.63 | 591,192 | 1,248,227 | 294,546 | **49.82%** |
| `5,000–9,999 bp` | 27,059 | 185.48 | 193,726 | 1,044,438 | 39,932 | **20.61%** |
| `>= 10,000 bp` | 20,653 | 640.79 | 599,002 | 934,791 | 27,061 | **4.52%** |
| **All (`>= 200 bp`)** | **3,930,592** | **3,226.02** | **5,845,901** | **1,812,110** | **4,523,043** | **77.37%** |

Contigs `< 1.5 kb` hold `59.71%` of `SP0`'s bases (`1.926 Gbp`) and `76.33%` of its called genes (`4.462M`), of which **`93.27%` (`4.162M`) are truncated fragments** running off a contig edge (*derived*).

### 4.2 Part B: Full 990-sample pilot at `--min-contig-len 1500` (`SPA1..SPU1`, `114.374 Gbp` stored)

Of the `990` non-empty assemblies, `967` hold at least one contig `>= 1,500 bp` (`23` tiny assemblies hold only `< 1,500 bp` contigs; *measured*). Those `967` samples were partitioned by level-1 habitat (`<= 19 Gbp` and `<= 5.5M` contigs per stage) into `9` stages (`SPA1, SPN1, SPQ1, SPQ2, SPH1, SPH2, SPT1, SPT2, SPU1`), enforcing exact cross-stage `contig_uid` ownership (`DEC-26`) via `write_owned_uids_bin` (`104` cross-habitat contigs = `38.03 Mbp` skipped in later stages and tracked in `cross_stage_occurrences.tsv`; *measured*). Default `--crispr` ran on the two stages with `<= 350k` contigs (`SPA1`, `SPU1`; `1.485 Gbp`).

*All values `*measured*` (`depth1/depth1_by_stage.tsv`, `depth1/depth1_summary.json`).*

| Stage | Habitat | Samples | Owned contigs | Stored Gbp | Genes (`M`) | Genes / Gbp | Frag. share | Intra-stage distinct prots (`M`) | Mean aa | Ingest CPU-min/Gbp | Genecall CPU-min/Gbp | CRISPR CPU-min/Gbp (arrays/Gbp) | Stored (`N1`) GB/Gbp | Called (`N3`) GB/Gbp | No-nt/IGR GB/Gbp | Peak step RSS (MiB) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `SPA1` | `air` | 4 | 8,305 | 0.0573 | 0.059 | 1,033,577 | 22.23% | 0.036 | 278.6 | 6.21 | 98.88 | 7.40 (157.0) | 0.3603 | 0.6915 | 0.3269 | 699.2 |
| `SPN1` | `anthropogenic` | 192 | 4,378,245 | 17.7353 | 20.051 | 1,130,557 | 34.09% | 18.772 | 262.9 | 4.49 | 59.27 | — | 0.4055 | 0.8584 | 0.4690 | 3,887.3 |
| `SPQ1` | `aquatic` | 121 | 4,991,014 | 18.9662 | 21.625 | 1,140,211 | 34.75% | 20.742 | 259.5 | 4.80 | 62.22 | — | 0.4138 | 0.8740 | 0.4807 | 4,353.8 |
| `SPQ2` | `aquatic` | 120 | 3,588,517 | 12.5770 | 14.441 | 1,148,196 | 35.94% | 14.086 | 255.9 | 4.92 | 56.19 | — | 0.4222 | 0.8890 | 0.4915 | 3,265.7 |
| `SPH1` | `host-associated` | 178 | 4,058,254 | 18.9262 | 18.987 | 1,003,208 | 27.66% | 16.496 | 251.7 | 4.98 | 57.82 | — | 0.3940 | 0.7802 | 0.3968 | 3,639.6 |
| `SPH2` | `host-associated` | 99 | 1,904,594 | 9.9698 | 10.192 | 1,022,318 | 26.79% | 8.370 | 296.1 | 4.70 | 62.44 | — | 0.3854 | 0.7858 | 0.4061 | 2,471.8 |
| `SPT1` | `terrestrial` | 79 | 5,401,073 | 18.8504 | 21.384 | 1,134,392 | 38.77% | 20.726 | 259.0 | 5.02 | 62.86 | — | 0.4215 | 0.8836 | 0.4867 | 4,867.7 |
| `SPT2` | `terrestrial` | 119 | 4,551,000 | 15.8639 | 17.869 | 1,126,421 | 37.25% | 17.438 | 245.8 | 4.97 | 62.75 | — | 0.4218 | 0.8758 | 0.4782 | 4,017.9 |
| `SPU1` | `unlabelled` | 55 | 293,922 | 1.4279 | 1.687 | 1,181,625 | 25.47% | 1.570 | 254.3 | 4.32 | 68.19 | 10.36 (215.0) | 0.3919 | 0.8551 | 0.4710 | 2,626.2 |
| **Total / weighted** | **9 stages** | **967** | **29,174,924** | **114.3740** | **126.296** | **1,104,237** | **33.98%** | **115.030** *(dedup across all 9)* | **258.9** | **4.840** | **60.662** | **10.245** (212.8) | **0.4098** | **0.8515** | **0.4601** | **4,867.7** |

Across all `29,174,924` owned contigs (`114.374 Gbp`), fragment share drops monotonically by length band (`depth1/depth1_summary.json`, *measured*):
- `1,500–4,999 bp` (`24.88M` contigs, `59.358 Gbp`): `1,230,497` genes/Gbp, **`50.38%` fragments** (`36,794,419 / 73,039,686`).
- `5,000–9,999 bp` (`2.79M` contigs, `18.923 Gbp`): `1,024,769` genes/Gbp, **`20.94%` fragments** (`4,059,646 / 19,391,441`).
- `>= 10,000 bp` (`1.50M` contigs, `36.093 Gbp`): `938,255` genes/Gbp, **`6.08%` fragments** (`2,058,194 / 33,864,809`).

### 4.3 Part C: Search-pack size per Gbp on unbinned contigs (`DEC-54 (3)`, `pack_probe_SPA1.json`)

Building release `rSPA1` and its verified search pack (`30fdb105ef20ed6868415f6c1abe0004`) on `SPA1` (`57,338,749 bp`, `8,305` contigs, `35,880` distinct proteins; *measured*, `depth1/pack_probe_SPA1.json`) gave:
- **`unique` (`part-0.faa.zst`, `35,880` proteins, `9,996,937` aa):** `6,497,671` bytes = **`0.1133 GB/Gbp`** (`11.2%` of pack). Note: `SPA1` has `625,755` distinct proteins/Gbp (`22%` intra-air protein sharing); at the full pilot's `1,005,737` distinct proteins/Gbp (`258.9 aa`), `unique` scales to **`~0.169 GB/Gbp`** (*derived*).
- **`igr` (`part-0.fna.zst` `2.60 MB` + `igr_meta.parquet` `2.44 MB`, `34,205` IGRs):** `5,040,621` bytes = **`0.0879 GB/Gbp`** (`8.7%` of pack).
- **`dark_orf` (`part-0.faa.zst` `33.69 MB` + `dark_orf_meta.parquet` `12.68 MB`, `459,824` dark ORFs, `39,399,560` aa):** `46,368,150` bytes = **`0.8087 GB/Gbp`** (**`80.1%` of pack**!).
- **Total search pack:** `57,909,582` bytes = **`1.0100 GB/Gbp`** (*measured*).

**Direct answer to `DEC-54 (3)` (and `C-38` Q3):** Keeping `dark_orf` in the search pack for `called` unbinned-contig stages multiplies pack size **`5.9×` to `8.9×`** (`1.01 GB/Gbp` vs `0.113–0.169 GB/Gbp` for `unique` proteins alone). For `7.8 Tbp` of `>= 1.5 kb` SPIRE contigs (`23.35 Tbp × 33.5%`, *derived*), `dark_orf` would add **`~6.3 TB`** to the pack on local scratch, whereas `unique` proteins alone require **`~1.3 TB`** before cross-tranche protein deduplication (and `< 1 TB` after deduplication). Exempting `called` stages from `dark_orf` (and optionally `igr`) in `pack.build` is essential before multi-Tbp unbinned stages enter a release.

### 4.4 Three `src/biobox/` defects uncovered at multi-million-contig scale (reported, not edited here)

1. **`src/biobox/l1/crispr.py:566` — 1-row Parquet row groups in scratch `contigs.parquet`:** `contig()` calls `self.writers["contigs"].write_batch(...)` with a 1-row `RecordBatch` on every contig. Because `ParquetWriter.write_batch` without buffering or `write_table(..., row_group_size=...)` flushes small row groups / footers, `SP0` (`3,930,592` contigs) wrote `3.93M` row-group entries (`2.19 GB` scratch file for `3.93M` 24-byte rows) and consumed `30,244 MiB` RSS + `>105 CPU-min` in DuckDB `_publish` (`*measured*`). Fix for a follow-up `L` card: buffer `122,880` contig rows per `RecordBatch` (`AGENTS.md` §4.6) before calling `write_batch`.
2. **`src/biobox/l0seq/lookup.py:51` — `KeyLookup(cache_row_groups=64)` LRU thrash cliff above `7.86M` 1-contig genomes:** `genecall` iterates `genome_input` ordered by `genome_uid = blake3(...)`, which is pseudorandom relative to `seq_locator`'s `raw_digest` sort key. While total stored contigs in a root fit inside `64 × 122,880 = 7,864,320` rows (`<= 64` row groups), `KeyLookup` caches 100% of `seq_locator` in ~`250 MB`. As soon as cumulative contigs cross `7.86M` (observed when `SPN1` `4.38M` + `SPQ1` `4.99M` reached `9.37M` = `78` row groups in a shared root), random `genome_uid` order misses the 64-entry LRU cache on `(78 - 64)/78 = 18%` of lookups, decompressing a `122,880`-row Arrow row group per miss (`>3×` CPU inflation and `>100 GB` aggregate RSS across 35 workers). Building each `<= 5.4M`-contig stage in its own isolated scratch root (`<= 44` row groups) restored 100% cache hit rate (`< 4.9 GB` RSS, `60.66 CPU-min/Gbp`). Fix for a follow-up `L` card: sort 1-contig genomes in `genecall` batch planning by `raw_digest` (or raise `cache_row_groups` / use a compact `mmap` index), which reduces memory to `O(1)` row groups at `S6`'s `111.5M` contigs!
3. **`src/biobox/l3/tapes.py` + `io.write_parts` tail shard `< 64 MiB` (`LY-04`):** On `SPQ2` (`3,588,517` contigs, `14,440,846` slots), `contig_tapes` cut at `max_file_bytes` into `part-0.parquet` (`635,177,986 B`) and `part-1.parquet` (`45,674,018 B < 67,108,864 B`), tripping `small_file_floor_bytes` when a second part barely spills past the 512 MiB target. Lowered only in scratch `/tmp/if2453/s3/spire/profile.yaml` (`small_file_floor_bytes: 1048576`); on `/stor`, `write_parts` should coalesce a trailing `< 64 MiB` part into the preceding part.

---

## 5. NR50 novelty against `r26.11.1` via 50%-identity MMseqs2 `linclust` (`novelty/`)

Streaming all distinct proteins from `r26.11.1` (`S1, P1..P6, H1, PL1`: **`49,102,041`** distinct proteins *measured*) and the `9` pilot depth-1 stages (**`115,030,133`** distinct proteins across `114.419 Gbp` *measured*; exact amino-acid overlap with `r26.11.1` = **`1,071,774`** proteins = `0.932%` of pilot proteins *measured*, `novelty/export_meta.json`) directly into `mmseqs createdb stdin` produced **`163,060,400`** combined distinct proteins (`50.843 GB` uncompressed FASTA equivalent, `0` bytes written as `.faa`).

### 5.1 `linclust` resource scaling at `163.06M` proteins (`novelty/linclust_steps.json`)

Running MMseqs2 50%-identity linear clustering (`--min-seq-id 0.5 -c 0.8 --cov-mode 1 --cluster-mode 2`, `--threads 56`, `--split-memory-limit 110G`, `align` split into `20` balanced index parts via `biobox.l2.cluster.cut`) completed in **`69.05` wall minutes** and **`1,230.91 CPU-min`** (`7.549 CPU-min per 10^6 proteins` *derived*), with every single process strictly `< 156 CPU-min` (`< 280 CPU-min` cap):
- `00_createdb` (stdin): `6.34 CPU-min`, `1,146 MiB` RSS.
- `01_kmermatcher`: **`155.95 CPU-min`** (`403.4 s` wall), **`96,578.0 MiB` (`101.27 GB`) peak RSS** (`621.1 B/protein` *derived*).
- `02_rescore1`: `100.25 CPU-min` (`315.0 s` wall), `74,947.4 MiB` RSS.
- `03_pre_clust` + `04..07` (`subdb`/`filterdb`/`rescore2`): `88.79 CPU-min` combined (`rescore2` = `40.50 CPU-min`, `53,223 MiB` RSS).
- `08_align_part_0000..0019` (`20` shards): **`828.64 CPU-min`** total (`38.2` to `44.9 CPU-min` per shard, max `58.9 s` wall/shard, max `48,112 MiB` RSS). *(Unsharded `align` would have consumed `~829 CPU-min` in one process and been killed by the pod watchdog at `369 CPU-min`.)*
- `09_clust` + `10_mergeclusters` + `11_createtsv`: `50.94 CPU-min` combined.

**Cluster totals (`novelty/novelty_summary.json`, *measured*):**
- Total NR50 `linclust` clusters across `163,060,400` proteins: **`56,757,541`**.
- Clusters containing `>= 1` `r26.11.1` protein: **`15,474,784`** *(note: `r26.11.1`'s full cascade produced `9,291,312` NR50 clusters on the same `49.10M` proteins, showing that 50%-id `linclust` alone produces **`1.666×`** as many clusters as the sensitive `-s 1..3` cascade — previewing `C-38c`)*.
- Shared clusters (`>= 1` `r26.11.1` member and `>= 1` pilot member): **`4,097,867`** (`26.48%` of `r26.11.1`'s `linclust` clusters).
- **Novel pilot-only NR50 clusters (`0` `r26.11.1` members):** **`41,282,757`** (**`360,804.53` novel NR50 families / Gbp** across `114.419 Gbp`).
- **Novel non-singleton NR50 clusters (`>= 2` pilot proteins, `0` `r26.11.1` members):** **`9,466,693`** (**`82,737.35` / Gbp**; `22.93%` of novel clusters).

### 5.2 Novel NR50 families per Gbp and saturation slope by level-1 habitat

*All values `*measured*` (`novelty/novelty_by_habitat1.tsv`, `novelty/accumulation_global.tsv`, `novelty/accumulation_by_habitat1.tsv`). First-half vs second-half marginal rates are averaged over `20` random sample permutations (`seed = 20261118`) within each habitat to measure within-habitat saturation cleanly.*

| Level-1 habitat | Samples (`>= 1.5 kb`) | Pilot bases (`>= 1.5 kb`, Gbp) | Novel NR50 clusters in habitat | Novel non-singleton (`>= 2`) | **Novel NR50 / Gbp** | **Novel non-singleton / Gbp** | 1st-half mean novel / Gbp | 2nd-half mean novel / Gbp | 2nd / 1st half retention (*derived*) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `unlabelled` | 55 | 1.428 | 705,468 | 215,331 | **493,912** | **150,757** | 379,665 | 520,285 | 1.370 |
| `aquatic` | 241 | 31.557 | 14,683,998 | 4,043,875 | **465,323** | **128,147** | 519,389 | 413,208 | **0.796** |
| `anthropogenic` | 192 | 17.749 | 7,543,985 | 2,413,905 | **425,040** | **136,003** | 469,951 | 379,326 | **0.807** |
| `terrestrial` | 198 | 34.724 | 14,147,716 | 4,480,037 | **407,439** | **129,020** | 465,288 | 349,345 | **0.751** |
| `host-associated` | 277 | 28.901 | 7,797,371 | 1,885,730 | **269,797** | **65,248** | 297,281 | 237,266 | **0.798** |
| `air` | 4 | 0.060 | 5,492 | 2,307 | **90,814** | **38,148** | 143,545 | 49,837 | **0.347** |
| **Global union** | **967** | **114.419** | **41,282,757** | **9,466,693** | **360,805** | **82,737** | **395,210** *(steps 1..483)* | **331,148** *(steps 484..967)* | **0.838** |

**Level-2 highlights (`novelty/novelty_by_habitat2.tsv`, *measured*):**
- Within `host-associated`, **`host-associated:animal host`** (`235` samples, `25.432 Gbp`) yields only **`229,724` novel NR50/Gbp** (`49,834` non-singleton/Gbp, dropping to `205,846/Gbp` in the second half) because GTDB (`P1..P6`) and `H1` already sample mammalian gut microbiomes heavily, whereas **`host-associated:plant host`** (`36` samples, `3.333 Gbp`) yields **`568,416` novel NR50/Gbp** (**`182,167` non-singleton/Gbp**, `98.1%` second-half retention: `567.1k -> 556.5k/Gbp`!).
- Highest-yielding major environmental level-2 strata (`>= 1 Gbp` in pilot): **`aquatic:benthic`** (`593,062/Gbp`, `93.5%` retention), **`terrestrial:subsurface`** (`647,903/Gbp`, `93.8%` retention), **`terrestrial:soil`** (`544,475/Gbp`, `211,095` non-singleton/Gbp, `91.9%` retention), **`terrestrial:grassland`** (`517,093/Gbp`, **`241,686` non-singleton/Gbp**), **`terrestrial:forest`** (`514,130/Gbp`), **`aquatic:brackish water`** (`500,702/Gbp`), **`anthropogenic:agriculture`** (`495,850/Gbp`, `229,831` non-singleton/Gbp), **`anthropogenic:wastewater`** (`492,308/Gbp`), and **`aquatic:marine`** (`482,250/Gbp`).
- **Stopping-rule implication (`C-38` §7 item 6):** Stratifying stopping rules at **level 2** rather than level 1 is essential — stopping `host-associated` as one bucket would prematurely cut off `plant host` (`556.5k/Gbp` marginal yield) because `animal host` (`205.8k/Gbp`) dominates sample counts `42:1`. Across environmental level-2 habitats, `114.4 Gbp` reduces marginal yield by only `6%–25%`, meaning diversity is far from saturated at `1,000` samples.

---

## 6. End-to-end `TniQ` (`PF06527`) search + `±10`-gene Pfam 38.2 round (`tniq/`)

Running `c42_tniq.py` over all `9` depth-1 pilot stages (`115,030,133` distinct proteins, `29,781,781,464` aa) demonstrated `DEC-49` / `C-38` §3 ("search first, annotate on the hit") end to end (`tniq/tniq_summary.json`, `tniq/tniq_loci.tsv`):

1. **Bait search (`TniQ`, `PF06527.17`, `M = 142`, GA `27.0` bits, `Z = domZ = 10^8`):**
   - Searched `115,030,133` distinct proteins in **`32.113 CPU-min`** (`1,282.6 s` wall; **`2.792 CPU-min per 10^7 proteins`** *measured*), identifying **`3,172` distinct `TniQ` proteins** at **`3,408` anchor loci** on **`3,161` contigs** across **`378` samples** (*measured*).
2. **Window & CRISPR extraction (`[-10, +10]` strand-oriented slots on `l3/tapes` + MinCED on hit contigs):**
   - Completed in **`5.186 CPU-min`** (`56.15 s` wall *measured*), collecting **`21,189` distinct window proteins** (`6,660,561` aa) — a **`5,428.8×` reduction in proteins** and **`4,471.4×` reduction in residues** relative to scanning the pilot in full (*derived*).
3. **Targeted Pfam 38.2 round scan (`biobox.l1.domains.scan`, gathering cutoffs, `Z = domZ = 10^8`, `56` threads):**
   - Completed in **`51.420 CPU-min`** (**`67.74 s` wall** *measured*; `7.720 CPU-min / 10^6 aa` including HMM load), annotating `14,373` window proteins with **`25,890` Pfam domain rows** (`round_domains.parquet` = `1.31 MB` *measured*).
   - **End-to-end CPU speedup:** `32.113 + 5.186 + 51.420 =` **`88.719 CPU-min`** *measured* vs **`213,833.2 CPU-min`** (`3,563.9 CPU-hours` *derived* at `7.18 CPU-min / 10^6 aa` over `29.782 × 10^9 aa`) for a full Pfam 38.2 scan — a **`2,410.2×` CPU reduction**!

### 6.1 `TniQ` locus classes discovered by habitat and SPIRE MAG membership

*All counts `*measured*` (`tniq/tniq_summary.json`, `tniq/tniq_loci.tsv`). Classification follows Peters et al. 2017 (`examples/crispr_cas_neighbourhoods/README.md`).*

| Locus class | `aquatic` | `terrestrial` | `anthropogenic` | `host-associated` | `unlabelled` | `air` | **Total loci** |
|---|---:|---:|---:|---:|---:|---:|---:|
| `I-F3_CAST` (`Cas5f/6f/7f/8f`, no `Cas1/2/3`) | **5** | 0 | 0 | **1** | 0 | 0 | **6** |
| `V-K_CAST` (`Cas12`-family / `Cas12f1-like_TNB` in window, no `Cas1/2/3`) | 0 | **1** | **2** | **1** | 0 | 0 | **4** |
| `I-B_CAST` (`Cas8b`, no `Cas1/2/3`) | 0 | **2** | 0 | 0 | 0 | 0 | **2** |
| `minimal_Cas_other` (Cas without `Cas1/2/3`) | **13** | **9** | **1** | **2** | 0 | 0 | **25** |
| `Cas_with_Cas1_Cas2_or_Cas3` | 1 | 0 | 0 | 0 | 0 | 0 | **1** |
| `Tn7_non_Cas` (`TnsA/B/C/TniB` without Cas) | 368 | 699 | 821 | 355 | 21 | 12 | **2,276** |
| `TniQ_only_or_truncated` | 212 | 452 | 297 | 124 | 9 | 0 | **1,094** |
| **All `TniQ` loci** | **599** | **1,163** | **1,121** | **483** | **30** | **12** | **3,408** |

Notable architectures in `tniq_loci.tsv` (*measured*):
- All `6` `I-F3_CAST` loci carry the canonical Peters et al. 2017 `+0:TniQ -> Cas5f(Csy2) -> Cas7f(Csy3) -> Cas6f(Csy4)` cascade (`5` in marine/brackish/briny water — e.g. `SAMN08987116` `26,381 bp` with an adjacent MinCED CRISPR array, `SAMEA2622717` `18,312 bp` with `TIR` + `OLD-like_TOPRIM` defense cargo, and three `5,581–6,335 bp` marine/brackish contigs — plus `SAMN19769733` `39,094 bp` in animal host with `-4:TnsB; -3:TniB` and a MinCED array).
- Both `I-B_CAST` loci (`SAMD00324618` `5,776 bp` and `SAMN26502303` `19,836 bp`, terrestrial) carry `TniQ -> Cas6 -> Cas8b`, and the `4` `V-K_CAST` hits carry `TnsA/TnsB/TniB + TniQ` alongside `Cas12f1-like_TNB` (`OrfB_IS605` / TnpB-family) domains (`23.9–58.0 kb`). Every single one of these `12` loci (`12 / 12`, `100%`) resides on a contig `>= 5,000 bp` (`8 / 12` on `>= 10,000 bp`, `4 / 12` in `5,000–9,999 bp`, `0 / 12` below `5,000 bp`; *measured*).

**Binned vs unbinned recovery:** Checking every `TniQ`-hit contig against all `7,029` downloadable SPIRE v1 MQ/HQ MAGs (`1,521,199` distinct MAG `contig_uid`s *measured*) published for the `378` hit-bearing samples (`339` with MAGs; `c42_fetch.py sample-mags`, `tniq/tniq_summary.json`, `tniq/tniq_loci.tsv`) showed that only **`351` of `3,408` loci (`10.30%`)** reside in SPIRE v1 MAGs, while **`3,057` loci (`89.70%`) — including `9 / 12` (`75.0%`) of the `I-F3/I-B/V-K` CAST loci and `20 / 25` (`80.0%`) of the `minimal_Cas_other` loci — exist exclusively in the unbinned contigs** and would have been missed completely by MAG-only catalogues (`GlobDB`, `GEM`, `SPIRE` MAGs; *measured*). Mobile genetic elements carry anomalous GC/tetranucleotide composition and differential read depth relative to host chromosomes, causing MetaBAT2 to exclude `75%–90%` of them from host bins even on `5–58 kb` contigs.

### 6.2 Contig-length truncation of `[-10, +10]` gene windows (`DEC-54 (4)` & `(12)`)

| Contig length band | `TniQ` loci (*measured*) | Share of `>= 1.5 kb` `TniQ` loci (*derived*) | Intact `21`-slot `[-10, +10]` windows (*measured*) | Window truncated on `5'` or `3'` side (*derived*) |
|---|---:|---:|---:|---:|
| `1,500–4,999 bp` | 1,664 | 48.83% | **0** (`0.00%`) | **100.00%** |
| `5,000–9,999 bp` | 740 | 21.71% | **0** (`0.00%`) | **100.00%** |
| `>= 10,000 bp` | 1,004 | 29.46% | **217** (`21.61%`) | **78.39%** |
| **Total (`>= 1,500 bp`)** | **3,408** | **100.00%** | **217** (`6.37%`) | **93.63%** |

---

## 7. Summary of replacements for `C-38_import_depths.md` §1 and inputs for PI decisions (`DEC-54`)

### 7.1 Updated table for `C-38_import_depths.md` §1 (`N1..N24`)

| # | Quantity | `C-38` prior | **C-42 measured / derived (`>= 1.5 kb`, unless noted)** |
|---|---|---|---|
| `N1` | Depth `stored` files / Gbp | `0.309 GB` | **`0.4098 GB/Gbp`** (`0.7568 GB/Gbp` unfiltered `SP0`) *measured* |
| `N2` | Depth `called` added files / Gbp | `0.405 GB` | **`0.4417 GB/Gbp`** (`0.6171 GB/Gbp` unfiltered `SP0`) *derived* |
| `N3` | Depth `called` total files / Gbp | `0.714 GB` | **`0.8515 GB/Gbp`** (`0.4601 GB/Gbp` without nt/IGR; `1.3739 GB/Gbp` unfiltered) *measured* |
| `N5` | `seq_shards` bytes / base | `0.403` (CAMI) | **`0.3510`** at `>= 1.5 kb` (`0.5621` unfiltered `SP0`) *derived* |
| `N6`, `N7` | Genes / Gbp of unbinned contigs | `2.16M` / `2.19M` | **`1.104M` genes/Gbp** at `>= 1.5 kb` (`1.812M` unfiltered `SP0`) *measured* |
| `N9` | Distinct proteins / gene; fragment share | `98.8%`; `83.5%` | **`91.08%`** across 9 stages (`93.62%` intra-stage); **`33.98%` fragments** at `>= 1.5 kb` (`77.37%` unfiltered) *measured* |
| `N10` | Contig length base shares | `<1 kb` `52.4%`, `>=10 kb` `20.0%` | `<1 kb` **`56.18%`**, `1–1.5 kb` **`10.31%`**, `1.5–5 kb` **`17.39%`**, `5–10 kb` **`5.54%`**, `>=10 kb` **`10.58%`** *measured* |
| `N11` | Ingest CPU / Gbp | `4.11` (`B-52`) / `4.5` (survey) | **`4.840 CPU-min/Gbp`** at `>= 1.5 kb` (`5.163` unfiltered) *measured* |
| `N12` | Gene calling (`meta` mode) CPU / Gbp | `60.4 CPU-min/Gbp` | **`60.662 CPU-min/Gbp`** at `>= 1.5 kb` (`68.865` unfiltered) *measured* |
| `N13` | Default CRISPR CPU / Gbp | `3.0 CPU-min/Gbp` | **`10.245 CPU-min/Gbp`** before fixing `crispr.py:566` *measured* |
| `N16` | Mean protein length | `323.6 aa` (`P2`) | **`258.90 aa`** at `>= 1.5 kb` *derived* (`29.782G aa / 115.030M` proteins) |
| `N19`, `N20` | `linclust` NR50 CPU & RAM | `~7.0 CPU-min/10^6`; `336 B/seq` | **`7.549 CPU-min / 10^6 proteins`**; **`621.1 B/protein`** (`101.27 GB` at `163.06M`) *measured* |
| `N21` | One-HMM search CPU | `8–11 CPU-min / 10.6M` | **`2.792 CPU-min / 10M proteins`** (streaming `protein_seq`) *measured* |
| `N22` | Search pack GB / Gbp | `0.495 GB/Gbp` | **`1.0100 GB/Gbp`** full (`0.8087` `dark_orf`, `0.0879` `igr`, **`0.1133–0.169` `unique` only**) *measured* |
| `N24` | 1,000-sample pilot scale | `160 Gbp`, `350M` proteins | **`341.501 Gbp`** raw; **`114.419 Gbp` & `115.030M` distinct proteins** at `>= 1.5 kb` *measured* |

### 7.2 Trade-off table for choosing `--min-contig-len L` (`DEC-54 (12)`, `C-38d`) across all of SPIRE (`23.35 Tbp` *derived*)

*Scaled from measured pilot shares (§3.1, §4.1, §4.2, §6.2) to SPIRE's `23.348 Tbp` estimated total (`*derived*`):*

| Candidate `--min-contig-len L` | Share of SPIRE bases kept | Total SPIRE bases kept (*derived*) | Depth-1 `called` storage at `0.85 / 0.78 / 0.74 GB/Gbp` (*estimate* from band density) | Depth-1 `genecall + ingest` CPU-hours (`65.5 CPU-min/Gbp`, *derived*) | Fragment-gene share (*measured*) | Share of `>= 1.5 kb` `TniQ` hits retained (*measured*) | Multi-gene operon context (`>= 4–5` genes) |
|---|---:|---:|---:|---:|---:|---:|---|
| `0` (unfiltered) | `100.0%` | `23.35 Tbp` | `~32.1 TB` (`1.374 GB/Gbp`) | `~28,800 h` | `77.4%` | — | Dominated by 1-gene fragments (`93.3%` below `1.5 kb`) |
| **`1,500 bp`** (Faure 2023) | **`33.50%`** | **`7.82 Tbp`** | **`~6.66 TB`** (`0.852 GB/Gbp`) | **`~8,540 h`** | **`34.0%`** | **`100.0%`** (`3,408`) | `1.5–5 kb` holds 2–4 genes (`50.4%` edge fragments); captures short cassettes (`Tn7` `tnsD/tniQ` borders) |
| **`5,000 bp`** (`C-38` §4.2 #4) | **`16.12%`** | **`3.76 Tbp`** | **`~2.94 TB`** (`~0.78 GB/Gbp`) | **`~4,110 h`** | **`11.4%`** (`5–10 kb` + `>= 10 kb` weighted) | **`51.17%`** (`1,744`) | `>= 5` genes per contig; halves storage and clustering proteins (`~36M` vs `115M` in pilot) while retaining **`12 / 12` (`100%`)** of `I-F3/I-B/V-K` CAST loci (`4` of which lie in `5.5–6.4 kb` and would be lost at `10 kb`, *measured* in `tniq_loci.tsv`) |
| **`10,000 bp`** | **`10.58%`** | **`2.47 Tbp`** | **`~1.83 TB`** (`~0.74 GB/Gbp`) | **`~2,700 h`** | **`6.08%`** | **`29.46%`** (`1,004`) | Only band with intact `21`-slot `[-10, +10]` windows (`21.6%` intact), but loses `70.5%` of `TniQ` loci |

---

## 8. Reproducible commands and scripts

All scripts live in `scripts/stress/build_stage_pod/c42_*.py` with unit tests in `tests/r1/test_c42_pilot.py`:
- `c42_fetch.py tables | assemblies | mag-probe | sample-mags` — polite rate-limited downloader (`User-Agent: curl`, `>= 0.5 s` interval, atomic `.part` rename, wire/disk MD5 & SHA-256 logs, `STOP`-aware).
- `c42_select.py` — full SPIRE v1 habitat aggregation and deterministic $\propto \sqrt{N_h}$ Hamilton selection (`seed = 20261118`).
- `c42_contigs.py` — streaming contig length/binning profiler and 16-bucket spilled DuckDB inter-sample `contig_uid` + strand-canonical repeat & `r26.11.1` overlap analyzer.
- `c42_inputs.py` & `c42_step4_driver.py` — 1-contig-per-genome meta-mode `genomes.json` / `objects.json` builder with `write_owned_uids_bin` cross-stage deduplication, resumable `SP0` + `SPA1..SPU1` `--depth called` stage builds, and `SPA1` search-pack profiler.
- `c42_novelty.py` & `c42_step56_driver.py` — streaming k-way merge directly into `mmseqs createdb stdin`, sharded 50%-identity `linclust`, spilled DuckDB + `array('I')`/bitset per-habitat novelty and permutation saturation curves.
- `c42_tniq.py` — streaming `pyhmmer` `TniQ` (`PF06527`) search (`Z = domZ = 10^8`), `contig_tapes` `[-10, +10]` window + MinCED CRISPR extraction, sharded Pfam 38.2 round scan (`biobox.l1.domains.scan`), Peters et al. 2017 CAST locus classifier, and SPIRE MAG overlap annotator.
