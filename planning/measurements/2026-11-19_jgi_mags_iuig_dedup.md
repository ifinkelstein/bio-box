# Deduplication and Staging of JGI `GEMs`, `IMG-UMAG` (Prokaryote & Eukaryote MAGs), and `IUIG` after `GTDB` R226 and `GlobDB` r226 (Card C-48d)

Date: 2026-11-19, `wilkcomp02`. Every number carries a label per `AGENTS.md` §6: *measured* (host, date, command/file), *derived* (formula and inputs), *estimate*, or *published*.
**2026-11-19 Source Freeze (`AGENTS.md` §8):** Zero external network downloads were executed. All archives were streamed directly and read-only (`:ro`) from `/stor/work/Wilke/if2453/all-the-data/jgi_20261118/` without unpacking any multi-hundred-GB `.tar` onto `/tmp` (`df -h /tmp`: `533–555 GiB` free throughout, *measured*).
**Scripts & tests:** `scripts/stress/build_stage_pod/c48d_{metadata,stream,skani_stage}.py`, `tests/r1/test_c48d_jgi_dedup.py` (`232 passed` on `wilkcomp02` with `tests/r0` and `tests/r3/test_cross_host.py`, *measured*).

---

## Executive Summary

1. **`GEMs` (`52,515` MAGs, `45,599` OTUs, `jgi_gems_*`) is 100% evaluated by `GlobDB` r226:** Every one of the `52,515` GEM MAGs (`100%` HQ/MQ: `0` below `50%` completeness or above `10%` contamination, *measured*) maps to one of the `45,599` GEM OTUs already present in GlobDB r226's `05_GEMOTU_dictionary.tsv` (`0` un-evaluated OTUs, *measured*). `15,200` GEM MAGs belong directly to OTUs retained in GlobDB r226 (`3,762` `GEMOTU_*` species reps) or GTDB R226, and the remaining `37,315` MAGs (`13,396` OTUs, `34.59 Gbp`) were already evaluated and collapsed at `>= 96%` fastANI into GTDB/mOTUs/SPIRE/RBG by GlobDB r226. **Recommendation:** Do not create separate prokaryote stages for `GEMs` after `S3`/`S5`/`G1..G41`; retain `gems/prophages/virsorter.tsv` (**`24,313` prophages** across `14,257` bins / `23,770` contigs, `686.96 Mbp`, median `23,891` bp; `8,787` directly on GlobDB/GTDB-kept GEM bins, *measured*) as a sidecar coordinate table.
2. **`IUIG` (`143,796` genome rows, `524.83 GB` zip) collapses `95.49%` at metadata level:** Of `126,346` prokaryote rows (`123,665` Bacteria + `2,681` Archaea; plus `15,544` Viruses, `1,185` Plasmids, `721` Eukaryota, *measured*), **`120,654` (`95.49%`, *derived*)** already match GTDB R226 / GlobDB r226 by NCBI assembly accession (`1,738`), NCBI taxid (`112,614`), or species binomial (`6,302`), leaving **`5,692` novel prokaryote candidate `taxon_oid`s** (*measured*). Streaming the first `25,001` archive members skipped `24,455` matched genomes without decompressing their inner `.tar.gz` (`902.36 MB` peak RSS) and extracted `250` novel candidates; `skani` at `95%` ANI / `50%` AF (`13.90 s`, `2` threads) collapsed `64` intra-source duplicates (`25.6%`) into **`186` locked novel species representatives** (`348.23 Mbp`, `redistributable = false`, *measured*).
3. **`IMG-UMAG` prokaryote bins (`408.90 GB` zip, `18,912` unrestricted metagenomes):** `7,247` of `7,304` GEM metagenomes (`99.22%`, *derived*) lie inside `IMG-UMAG`, while **`11,665`** are post-GEM (2020–2024) metagenomes (*derived*: `18,912 − 7,247`). A bounded `500`-bin streaming pilot (`224.62 s`, `906.64 MB` peak RSS, *measured*) showed `100%` HQ/MQ quality (`55` HQ + `441` MQ + `4` unlabelled; median completeness `80.91%`, contamination `1.15%`, `1.31 Gbp`, `136` contigs shared with `PL2`/`G1`/`V1`). Running `skani` (`13.42 s`, `2` threads) on `250` extracted bins against `1,000` `G1` refs dropped `1` GlobDB match (`BCRBG_06715`, `95.34%` ANI) and `3` intra-pilot duplicates, locking **`246` species reps** (`664.23 Mbp`, `redistributable = false`, *measured*).
4. **`IMG-UMAG` eukaryote bins (`38.57 GB` zip) — FULL inventory of all `691` bins:** True O(1) nested tar+gzip streaming processed all **`691` eukaryote bins** (`14,435,447,512` bp = **`14.44 Gbp`**, `1,227,548` contigs) in `2,736.3 s` (`45.6` CPU-min) at **`905.01 MB` peak RSS** — even across a single `5.28 GB` compressed mega-bin (`./3300059473_5929.tar.gz`) (*measured*). Overlap against all stored stages (`S1`, `P1..P6`, `H1`, `PL1`) and staged inputs (`V1`, `PL2`, `G1`) is **`0` shared contigs (`100%` novel)** (*measured*). GFF diagnostics across all `691` bins prove severe prokaryotic-caller intron fragmentation (`0` intron/exon features, `366` bp median CDS, `64.95%` coding density, down to `39.8%` in `85–102 Mbp` bins, *measured*). **Recommendation:** Stage separately as **`E1` (`~10.0 Gbp`) + `E2` (`~4.44 Gbp`)**, `redistributable = false`, at **`--depth stored`** (`C-38a`) after adding `kind = "eukaryote"` to `schemas.py` before `rS0.1`, and defer gene calling until an intron-aware eukaryotic caller (MetaEuk / AUGUSTUS / Helixer) is integrated.

---

## 1. Stage 1: Metadata Deduplication (`c48d_metadata.py`)

Command (`wilkcomp02`, 2026-11-19, `bio-box_faec1f4.sif`, `OMP_NUM_THREADS=1`, `MemoryMax=16G`, output `/tmp/if2453/c48d/meta/metadata_summary.json`):
`c48d_metadata.py` indexed **`732,475` GTDB R226 rows (`143,614` species reps)**, **`306,260` GlobDB r226 reps** (including all `45,599` `05_GEMOTU_dictionary.tsv` rows), and **`126,357` IMG/M rows** (`43,697` IMG OID → `GCA_`/`GCF_` links; `106,842` IMG OID → GOLD ecosystem paths; all *measured*).

### 1.1 `GEMs` (`jgi_gems_meta_prophages_otus_20261118.tar`)

| Metric | Value | Label & Source |
|---|---:|---|
| Total GEM MAGs (`genome_metadata.tsv`) | `52,515` | *measured* (`c48d_metadata.py`, 2026-11-19) |
| HQ/MQ (`completeness >= 50%`, `contamination <= 10%`) | `52,515` (`100.0%`) | *measured* (`0` LQ, `0` `>10%` contam) |
| Total OTUs (`otu_representative_genomes.tsv`) | `45,599` | *measured* (matches `05_GEMOTU_dictionary.tsv`) |
| MAGs in OTUs already kept in GlobDB r226 (`3,762` `GEMOTU_*`) or GTDB R226 | `15,200` | *measured* |
| Remaining MAGs evaluated & collapsed by GlobDB r226 (`05_GEMOTU_dictionary.tsv`) | `37,315` (`13,396` OTUs, `34.59 Gbp`) | *measured* (`23,919` non-rep MAGs + `13,396` top OTU MAGs) |
| GEM OTUs **not** evaluated by GlobDB r226 | **`0`** | *measured* |
| Unique IMG metagenomes in GEMs | `7,304` (`7,247` in IMG-UMAG) | *measured* |
| VirSorter prophages (`prophages/virsorter.tsv`) | **`24,313`** (`14,257` bins, `23,770` contigs) | *measured* (`686,959,122` bp total, median `23,891` bp) |
| ↳ Prophages on GlobDB/GTDB-kept GEM bins | `8,787` | *measured* |
| ↳ Prophages on top MAGs of GlobDB-collapsed OTUs | `4,132` | *measured* |

*Note on file format quirk:* In NERSC's `gems/genomes/genome_metadata.tsv`, the header has 20 column names while each data row has 21 fields (an unnamed GTDB taxonomy column sits between `otu_id` and `ecosystem`). `c48d_metadata.py` detects and corrects this shift automatically.

### 1.2 `IUIG` (`img_genomes_tarballs_metadata_111023.txt`) and `IMG-UMAG` (`unrestricted_metagenomes_for_bin_tarballs.02062024.txt`)

| Dataset / Stratum | Count | Share / Notes | Label |
|---|---:|---|---|
| **IUIG total rows** | `143,796` | `123,665` Bacteria, `2,681` Archaea, `15,544` Viruses, `1,185` Plasmids, `721` Eukaryota | *measured* |
| ↳ IUIG prokaryotes (`Bacteria + Archaea`) | `126,346` | `23,811` Finished, `118,353` Permanent Draft, `1,623` Draft (across all domains) | *measured* |
| ↳ Matched to GTDB R226 / GlobDB r226 | `120,654` (`95.49%`) | `1,738` by `GCA_`/`GCF_`, `112,614` by NCBI taxid, `6,302` by species binomial | *measured* / *derived* |
| ↳ **Novel IUIG prokaryote candidates** | **`5,692`** (`4.51%`) | Top phyla: `Pseudomonadota` `1,507`, `Actinomycetota` `735`, `Bacillota` `560`, `unclassified` `439`, `Bacteroidota` `334`, `Euryarchaeota` `207`, `Cyanobacteriota` `174`, `Chloroflexota` `161`, `Planctomycetota` `126`, `Thermoproteota` `111` | *measured* / *derived* |
| **IMG-UMAG unrestricted metagenomes** | `18,912` | `7,247` shared with GEMs (`99.22%` of GEM metagenomes); **`11,665` post-GEM metagenomes** | *measured* / *derived* |
| ↳ Top IMG-UMAG GOLD ecosystems | — | `Environmental / Aquatic` `7,680`, `Environmental / Terrestrial` `4,344`, `Host-associated / Human` (`1,798 + 784` = `2,582`), `Host-associated / Plants` `1,110`, `Engineered / Built environment` `974`, `Engineered / Bioreactor` `618`, `Engineered / Wastewater + WWTP` `352` | *measured* / *derived* |

---

## 2. Stage 2: Bounded-Memory Streaming & `contig_uid` Overlap (`c48d_stream.py`)

All four archives were streamed directly from `/stor/work/Wilke/if2453/all-the-data/jgi_20261118/` (`:ro`) using `tarfile.open(..., mode="r|")` / nested `zipfile -> tarfile(r|) -> tarfile(r|gz)` with O(1) line-by-line FASTA and GFF parsing (`OMP_NUM_THREADS=1`, `MemoryMax=16G`, `wilkcomp02`, 2026-11-19). Every contig was hashed with `biobox.ids.contig_uid` and checked against all `9` stored stages in `/tmp/if2453/r10/root` (`S1`, `P1..P6`, `H1`, `PL1`) plus staged inputs `V1`, `PL2`, and `G1`.

| Stream Run | Members Seen | Skipped (Un-decompressed) | Processed | Total Bases (Gbp) | Contigs | Shared `contig_uid`s (Breakdown) | Median Size / N50 | Wall Time (s) | Peak RSS (MB) |
|---|---:|---:|---:|---:|---:|---|---|---:|---:|
| **`GEMs` full tar** (`jgi_gems_52515_fna_20261118.tar`) | `52,515` | `39,119` | `13,396` | `34.59` (`34,590,779,522` bp) | `2,241,001` | `4,243` (`PL2`: `2,831`, `G1`: `1,335`, `V1`: `76`, `H1`: `1`) | `2.33 Mbp` / `26.3 kbp` | `1,009.84` | `899.84` |
| **`IUIG` filtered pass** (`jgi_iuig_genomes_20261118.zip`) | `25,001` | `24,455` | `545` | `1.24` (`1,241,437,389` bp) | `108,375` | `186` (`PL2`: `184`, `P3`: `2`) | `1.85 Mbp` / `27.4 kbp` | `1,050.06` | `902.36` |
| **`IMG-UMAG` prok pilot** (`jgi_img_umag_prok_bins_20261118.zip`) | `501` | `0` | `500` (`100%` HQ/MQ) | `1.31` (`1,309,115,785` bp) | `103,353` | `136` (`PL2`: `94`, `G1`: `36`, `V1`: `6`) | `2.36 Mbp` / `15.1 kbp` | `224.62` | `906.64` |
| **`IMG-UMAG` euk FULL** (`jgi_img_umag_euk_bins_20261118.zip`) | **`691`** | `0` | **`691`** | **`14.44`** (`14,435,447,512` bp) | `1,227,548` | **`0`** (`100%` novel across all stages) | `17.60 Mbp` / `10.8 kbp` | `2,736.30` | `905.01` |

All figures above are *measured* (`wilkcomp02`, 2026-11-19, `/tmp/if2453/c48d/stream_*.json` and `stream_*.tsv`). Peak RSS stayed flat at **`899.8–906.6 MB`** across all four runs (of which ~880 MB is the `r10/root + V1 + PL2 + G1` `contig_uid` owner hash map).

### 2.1 `IMG-UMAG` Eukaryote Bin Inventory & Intron/Gene-Calling Diagnostics (`691` Bins, `14.44 Gbp`)

- **Size & assembly distribution (*measured*):** `691` bins spanning `4.70 Mbp` (min) to `102.30 Mbp` (max: `3300030522_3`, `4,138` contigs, N50 `39.6 kbp`), median `17.60 Mbp`, median N50 `10,825` bp, median GC `52.59%`. Largest compressed member is `./3300059473_5929.tar.gz` (`5,278,357,062` B compressed, *measured*), confirming why O(1) inner-tar streaming is mandatory (`MemoryMax=16G`).
- **Prokaryotic gene-caller artifact in JGI GFFs (*measured*):** Every one of the `691` bins was annotated by JGI using prokaryotic ORF callers (`Prodigal v2.6.3` in `587+160` bins, `GeneMark.hmm-2` in `592+80` bins) with **`0` `intron`, `exon`, or `mRNA` features**. Across all `691` bins, median CDS length is only **`366` bp (`122` aa)** (compared to `803–971` bp in `G1`/`P2` prokaryotes) and median apparent coding density is **`64.95%`** (falling to `39.84%–45.48%` in the `85–102 Mbp` bins). Running `pyrodigal` on these bins would fill `protein_seq` with exon-truncated fragments and spurious short ORFs inside introns.
- **Schema & gene-calling decision:** `src/biobox/registry/schemas.py` currently restricts `genome_input.kind` to `("prokaryote", "virus", "plasmid")` (while `_EDIT_TYPES` already reserves `"intron"`). `E1..E2` should be staged separately (`redistributable = false`), ingested at **`--depth stored`** (`C-38a`) once `"eukaryote"` is added to `genome_input.kind` (before `rS0.1` lock), and gene-called only after an intron-aware eukaryotic caller (MetaEuk, AUGUSTUS, or Helixer) is wired into `l1/genecall`.

---

## 3. Stage 3: `skani` 95% ANI Species-Rep Selection & Locked Pilot Manifests (`c48d_skani_stage.py`)

Both pilots ran `skani 0.4.0` (`sketch --separate-sketches`, `dist` vs `1,000` GlobDB `G1` refs, and `triangle -E` all-vs-all at `ANI >= 95.0%`, `AF >= 50.0%`, `2` threads, `wilkcomp02`, 2026-11-19; *measured*):

| Pilot Set | Input Candidates | Wall Time (s) | Dropped vs `G1` Refs (`>=95%` ANI) | Collapsed Intra-Pilot | Locked Novel Species Reps | Locked Bases | Locked Contigs | `redistributable` |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| **`IUIG` pilot** (`/tmp/if2453/c48d/pilot_iuig/`) | `250` | `13.90` | `0` | `64` (`25.6%`) | **`186`** | `348,229,181` bp | `20,087` | **`false`** |
| **`IMG-UMAG` prok pilot** (`/tmp/if2453/c48d/pilot_umag_prok/`) | `250` | `13.42` | `1` (`BCRBG_06715`, `95.34%` ANI) | `3` (`1.2%`) | **`246`** | `664,229,694` bp | `54,142` | **`false`** |

### Locked Pilot SHA-256 Digests (`inputs_SUMS`, *measured*, `wilkcomp02`, 2026-11-19)

**`IUIG` pilot (`/tmp/if2453/c48d/pilot_iuig/inputs_SUMS`):**
- `genomes.json`: `c0c423a2e7dd07972694106a0f78cece0a71c5ddc18e320919e195738e06f82a`
- `lock.json`: `0a6b2e6648a9567d9e81942448f38b06a99ffc4b0ef58dc54eadfb6a942f26f0`
- `objects.json`: `aef715d25277be52c2d04887891cb1507e1e236e917e86b5a058adff78311e98`
- `selection.tsv`: `86774a47c9ef0139d159a37c14fb1e35913ff6302b75917a7109c8a89bda4369`

**`IMG-UMAG` prokaryote pilot (`/tmp/if2453/c48d/pilot_umag_prok/inputs_SUMS`):**
- `genomes.json`: `ac0a00b9b427c5b58d10f4c90c6d9b2cc0f36c38048fbee3f1151e5b6c6b7fce`
- `lock.json`: `42e82febd8e374d6c5df8913231ab0b847401e1834e0b159fd635ff23a26cf25`
- `objects.json`: `295ad1f8dde47cc45ca1202cb921014f4f8cb058952b2a2fb78551932ba2a4a7`
- `selection.tsv`: `db8411dd7887e3e52c7858d321195a974175b945dc3c103d829445fea406b6d9`

---

## 4. Proposed `~10 Gbp` Stages (Order: `S3` → `S5` → `G1..G41` → `J1..Jn` / `E1..E2`)

All JGI-derived stages carry **`redistributable = false`** (`DEC-51`, JGI Data Utilization Policy) and are placed **after** `S3`, `S5` (GTDB R226) and `G1..G41` (GlobDB r226) so public INSDC/GlobDB representatives own every shared contig under `DEC-26`:

| Proposed Stage(s) | Source | Kind | Est. / Measured Genomes | Est. / Measured Size | Recommended Depth | Notes |
|---|---|---|---:|---:|---|---|
| *(none — sidecar only)* | `GEMs` (`52,515` MAGs) | — | `0` new species | `0` Gbp | sidecar (`virsorter.tsv`) | `100%` of `45,599` GEM OTUs already evaluated in GlobDB r226; retain `24,313` VirSorter prophage coordinates |
| **`J1`** | `IUIG` novel prokaryotes | `prokaryote` | `~4,235` reps (*derived*: `5,692 × 0.744`) | `~9.65 Gbp` (*derived*: `4,235 × 2.28 Mbp`) | `called` (or `annotated` w/ B-47) | `1` stage of `~10 Gbp`; `redistributable = false`; `pyrodigal` single mode (screen CPR/DPANN for code 25/4) |
| **`J2..J15`** *(est.)* | `IMG-UMAG` post-GEM prokaryote MAGs | `prokaryote` | `~50k–60k` novel reps (*estimate* after full `skani` vs GTDB+GlobDB) | `~120–140 Gbp` (*estimate*) | `called` (`C-38a`) | `~12–14` stages of `~10 Gbp`; `redistributable = false`; stream full 409 GB zip in 1 pass (~45–60 CPU-min) once `G1..G41` refs are on disk |
| **`E1`, `E2`** | `IMG-UMAG` eukaryote bins (FULL) | `eukaryote` *(new enum)* | **`691` bins** (*measured*) | **`14.44 Gbp`** (`E1`: `~10.0 Gbp`, `E2`: `~4.44 Gbp`, *measured*/*derived*) | **`stored`** (`C-38a`) | `100%` novel (`0` shared contigs); `redistributable = false`; do **not** run `pyrodigal`; await intron-aware eukaryotic caller |
