# C-50: RNA viruses, Obelisks, and viroid-like circular RNAs

A Tier-H design proposal for the PI (card `C-50`, `DEC-60 (2)`, `DEC-69`, card `C-45` question (1)), 2026-11-24.
No code, no `sources.yaml`, no `TASKS.md`, and no `DECISIONS.md` changes are made by this document.
Every number carries an `AGENTS.md` §6 label (`*published*`, `*measured*`, `*derived*`, or `*estimate*`) and its source.
"From the code" means read statically from the repository tree on `task/c50-design`, not executed.

---

## 1. Executive Summary

- **Tiny footprint, outsized biological novelty.** The entire Phase 1 curated RNA & circular-RNA wave (`RNA1` RefSeq `Riboviria` + `VIRE` `Riboviria` + `ELGV` RNA; `OBL1` Obelisks; `VIRD1` ViroidDB + Lee/Neri 2023 cccRNAs + *Ambiviricota*) totals **~0.25–0.35 Gbp** (*estimate*) — less than 0.2% of the DNA corpus — and even Phase 2 (`RNA2` LucaProt species reps + `RNA3` RVMT v3 + Serratus nucleotide contigs) adds only **~2.0–4.0 Gbp** (*estimate*).
- **Zero L0-seq or schema-lock breakage.** `bio-box` already supports multi-contig genomes (`genome_uid` over segments), circular origin-wrapping coordinates (`end > contig_len`), 0-gene contigs in `contig_tapes`, multi-library `l1_domains`, and whole-contig circular `Infernal` scans in `l1/ncrna.py`. Only four targeted extensions are needed:
  1. **Alphabet normalization at staging**: fold `U -> T` / `u -> t` before L0 ingest (recorded in `u_to_t_normalized.tsv`), keeping L0-seq exact-byte round-trips, `revcomp()`, `contig_uid` deduplication, and `pyhmmer`/`pyrodigal` DNA-alphabet assumptions intact.
  2. **Two new `genome_input.kind` / `genome_dim.kind` enum values**: `rna_virus` and `cccrna` (additive under `LY-09`, zero column changes).
  3. **RNA/cccRNA ORF calling**: `pyrodigal-gv` is a prokaryotic/phage RBS caller that truncates eukaryotic RNA polyproteins and drops overlapping/ambisense ORFs (*Obelisks*, *Ambiviricota*, *Ribozyviria*). We recommend a stratified pilot benchmark (`C-50a`) before freezing a deterministic 6-frame ORF caller built on pinned `biobox.search.orfs.six_frame_orfs` (`pyhmmer.easel`, zero new container dependencies) alongside `gv-meta` for prokaryotic RNA phages (*Leviviricetes*).
  4. **A frozen ribozyme CM set (`ribozyme-1`) from already-pinned `Rfam 15.1`**: `prok-2` contains only 3 ribozyme/hairpin models, all intergenic-only (`whole_genome=0`). Scanning ~22 Rfam 15.1 self-cleaving ribozyme & viroid CMs with `whole_genome=1` on RNA/cccRNA stages costs `< 20` CPU-minutes (*estimate*) and captures hammerhead, HDV, twister, twister-sister, hairpin, hatchet, and pistol ribozymes across circular origins and inside/overlapping ORFs.
- **Six PI decisions (§8)** unblock implementation cards `C-50a` through `C-50f`.

---

## 2. Primary-Source Inventory, Overlaps, and Non-Redundant Take

Every source below was verified against primary publications and public repository records on 2026-11-24. Exact file checksums and byte counts marked *verify on fetch* are filled in `sources.yaml` immediately before download per `DEC-50`, `DEC-51`, and `AGENTS.md` §8.

### 2.1 Source Table

| Source | Primary Reference & Deposit URL / DOI | Contents & Unit | Sequence / Species Count | Size (compressed / bp) | License & `redistributable` |
|---|---|---|---|---|---|
| **NCBI RefSeq RNA viruses** (`Riboviria`) | NCBI RefSeq Viral FTP: `https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.{k}.1.genomic.fna.gz` (+ `.gbff.gz` for segment/taxid grouping) | Curated complete genomes & segments across *Riboviria* (`Orthornavirae` + `Pararnavirae`) | `~19,000–19,600` total RefSeq viral records, of which **`~8,600–11,000+`** are *Riboviria* (*published*, NCBI RefSeq 2024–2025) | `~0.10–0.15 Gbp` (*estimate*, `~10k × ~10 kb` mean) | Public domain (INSDC / NCBI); **`redistributable: true`** |
| **LucaProt** | Hou, He et al. 2024 *Cell* (*"Using artificial intelligence to document the hidden RNA virosphere"*). Figshare: `https://doi.org/10.6084/m9.figshare.26298802` (*"The LucaProt-Related Resources"*); code `https://github.com/alibaba/LucaProt` | Metatranscriptomic RNA viral contigs & RdRp-clustered species from `10,487` global metatranscriptomes (`180` supergroups; max length `47,250 nt`) | **`513,134`** viral contigs; **`161,979`** putative species (`>90%` RdRp aa ID), of which **`70,458`** are novel vs RVMT/Serratus/NCBI (*published*) | `~1.0–1.5 Gbp` all contigs; `~0.4–0.6 Gbp` species reps (*estimate*, `~2.5–3 kb` mean) | Data: **CC BY 4.0** (Figshare); code: Apache-2.0; **`redistributable: true`** |
| **RVMT v3** | Neri, Wolf, Roux et al. 2022 *Cell* 185:4023 (*"Expansion of the global RNA virome..."*). Zenodo v3: `https://doi.org/10.5281/zenodo.7368133` (`RVMT_Zenodo_v3.tar.gz`; earlier `10.5281/zenodo.6091357`); portal `riboviria.org`; code `https://github.com/UriNeri/RVMT` | Metatranscriptomic RNA viral contigs, RdRp ORFs, and HMMs from `5,150` metatranscriptomes | **`>2,500,000`** RNA viral contigs; **`~330,000`** RdRps (*published*) | `9.5 GB` (`RVMT_Zenodo_v3.tar.gz`, *published*); curated RdRp subset `~1.1 GB` FASTA (*published*) | Code: MIT; Zenodo data: open-access CC BY (*verify exact badge on fetch*); **`redistributable: true`** (pending badge check) |
| **Serratus / palmDB** | Edgar, Babaian et al. 2022 *Nature* (*"Petabase-scale sequence alignment..."*). `palmDB`: `https://github.com/rcedgar/palmdb` (`sotus.palmprint.faa`, `unique.palmprint.faa`); nucleotide assemblies: `s3://serratus-public/` (e.g. `rdva_v0.2.fa.lz4`) | `palmDB` = **protein-only** `~100 aa` RdRp palmprints; `s3://serratus-public/` = SRA RdRp-hit nucleotide contigs | **`131,957`** novel RNA virus sOTUs (`<90%` palmprint ID, *published*) | `palmDB` `< 50 MB` (protein-only); S3 nucleotide contigs `~0.5–1.5 Gbp` (*estimate*) | Data: **CC0 / CC BY 4.0** (INSDC-derived); code: GPLv3; **`redistributable: true`** |
| **Obelisks** | Zheludev, Edgar, Fire et al. 2024 *Cell* (*"Viroid-like colonists of human microbiomes"*). Stanford Digital Repository: `https://purl.stanford.edu/wb363nt3637` (DOI `10.25740/wb363nt3637`); `VNom` code `https://github.com/Zheludev/VNom`; *S. sanguinis* raw `PRJNA1129866` | `~1 kb` rod-like circular RNAs encoding `Oblin-1` (`~200 aa`) and `Oblin-2` (`~50–60 aa`), some with variant HHR-III ribozymes | **`29,959`** distinct Obelisk clusters (`90%` nt ID, *published*) | **`~30 Mbp`** (*derived*: `29,959 × ~1,000 nt`) | SDR open-access (CC BY / CC BY-NC — *verify exact SDR metadata badge on fetch*); **`redistributable`** set per badge (`DEC-51`) |
| **Viroid-like cccRNAs** (Lee/Neri 2023) | Lee, Neri, Roux, Koonin et al. 2023 *Cell* (*"Mining metatranscriptomes reveals a vast world of viroid-like circular RNAs"*). Zenodo: `https://doi.org/10.5281/zenodo.6859104`; code `https://github.com/Benjamin-Lee/vdsearch` | Ribozyme-bearing circular RNAs (`0.2–5 kb`) from `5,150` metatranscriptomes + `1,344` plant transcriptomes | **`11,378`** viroid-like cccRNAs (*published*) | **`~10–15 Mbp`** (*estimate*, `~1 kb` mean) | Code: MIT; Zenodo data: open-access CC BY (*verify exact badge on fetch*); **`redistributable: true`** (pending badge check) |
| **Ambiviruses** (*Ambiviricota*) | Forgia, Navarro, de la Peña et al. 2023 *Nature Communications* 14:2591 (`https://doi.org/10.1038/s41467-023-38301-2`) Supplementary Data + NCBI GenBank `txid` *Ambiviricota* query | `4–5 kb` circular ssRNA fungal viruses carrying **both** RdRp and self-cleaving ribozymes (ambisense ORFs) | **`378`** RdRp sOTUs (*published*; `~400–600` complete/near-complete GenBank genomes, *estimate*) | **`~2–5 Mbp`** (*estimate*, `~500 × 4.5 kb`) | **CC BY 4.0** (*Nat. Commun.* open access) + INSDC; **`redistributable: true`** |
| **ViroidDB** | Lee et al. 2022 *NAR* / 2023 update (`https://viroids.org/db/latest/all.fasta`, `all.json`); archived Zenodo `https://doi.org/10.5281/zenodo.5202945`; code `https://github.com/Benjamin-Lee/viroiddb` | Curated viroids (`Avsunviroidae`, `Pospiviroidae`), `Ribozyviria` (deltaviruses), circular satRNAs, retrozymes | **`~9,891`** sequences (`9,354` viroids: `5,284` *Avsunviroidae*, `3,980` *Pospiviroidae*, `90` unclass.; `453` deltaviruses; `~35` satRNAs; *published*) | **`~4.1 Mbp`** (*derived*: `9,354 × 350 nt + 453 × 1,700 nt`) | MIT / INSDC; **`redistributable: true`** (pin snapshot date + sha256 under `DEC-11`) |
| **ELGV** (RNA fraction) | Zeng, Almeida et al. 2024 *Nature Communications* (*"A metagenomic catalog of the early-life human gut virome"*). Figshare: `https://doi.org/10.6084/m9.figshare.21901557` (`ELGV_representatives_82141.fna.gz`, `378.3 MB`; `ELGV_all_160478.fna.gz`, `721.5 MB`) | Infant gut DNA + RNA viral catalog (`160,478` sequences, `82,141` species reps); **only `Riboviria` subset enters `C-50`** | `82,141` total reps (`160,478` total, *published*); RNA fraction **`~500–2,000`** reps (*estimate*, mostly eukaryotic/plant/phage RNA viruses) | **`~5–15 Mbp`** RNA fraction (*estimate*) | **CC BY 4.0** (Figshare); **`redistributable: true`** |
| **VIRE `Riboviria`** (already deferred from `C-55`) | Fullam et al. 2026 *NAR* (`https://swifter.embl.de/~fullam/vire/archives/Riboviria_genomes.tar`, already listed in `sources.yaml` under `vire_genomes`) | Metagenomic/metatranscriptomic `Riboviria` genomes explicitly deferred from `V3` (`V3_left_for_later.tsv.gz`, `reason=riboviria_deferred_c50`) | **`8,056`** genomes across **`3,480`** vOTUs (*measured*, card `C-55`, `CATALOG_GOAL.md` §4.1) | **`~20–40 Mbp`** (*estimate*, `3,480` vOTU reps `~10–20 Mbp`) | **CC BY-SA 4.0**; **`redistributable: true`** |

### 2.2 Mutual Overlaps & Deduplication Strategy

1. **`RefSeq Riboviria` vs `LucaProt` / `RVMT` / `Serratus`**: `RefSeq` holds the curated taxonomic backbone (including complete segmented genomes). Staging `RefSeq` first (`RNA1`) guarantees that every RefSeq contig owns its `contig_uid` and `genome_uid` (`DEC-26`).
2. **`LucaProt` vs `RVMT` vs `Serratus` (high SRA overlap)**: All three mined overlapping subsets of public SRA metatranscriptomes (`RVMT`: `5,150` metatranscriptomes; `LucaProt`: `10,487` metatranscriptomes, reporting `70,458` of its `161,979` species as novel over RVMT/Serratus/NCBI, *published*). Because different assemblers (`MEGAHIT` vs `rnaviralSPAdes`) produce slightly different contig ends for the same virus, exact `contig_uid` deduplication catches identical contigs while `NR90`/`NR50` release clustering collapses near-identical RdRps across stages. Taking **`LucaProt` species representatives first (`RNA2`)**, then **`RVMT` v3 RdRp-bearing representatives (`RNA3`)** (excluding `contig_uid`s already in `RNA1`/`RNA2`), and finally **`Serratus` `s3://serratus-public/` nucleotide contigs** (deferring protein-only `palmDB` under `DEC-60 (2)`) maximizes diversity per Gbp while keeping every non-representative contig recorded in `CATALOG_GOAL.md` §4.1 (`DEC-59 (3)`).
3. **`ViroidDB` vs `Lee/Neri 2023` vs `Ambiviricota` vs `Obelisks`**:
   - `Obelisks` (`OBL1`) are phylogenetically distinct from both RdRp viruses and classical viroids (encoding `Oblin-1`/`Oblin-2`, no RdRp); expected overlap with `ViroidDB` and `RNA1..3` is near zero (*estimate*).
   - `ViroidDB`, `Ambiviricota` (Forgia 2023), and `Lee/Neri 2023` overlap on known plant viroids, *Ribozyviria* (deltaviruses), and fungal ambiviruses (~300+ ambivirus-like cccRNAs appear in both Lee/Neri 2023 and Forgia 2023). Precedence order **`ViroidDB -> Ambiviricota -> Lee/Neri 2023`** inside `VIRD1` gives curated INSDC/ICTV records ownership of shared `contig_uid`s (`DEC-27`).
4. **`ELGV` DNA vs RNA**: `ELGV` is predominantly DNA bacteriophages overlapping `V1` (UHGV) and `CG1` (CGVR). Only the `Riboviria` rows of `ELGV_representatives_82141` enter `RNA1`; its DNA fraction and non-representative remainder are recorded in `CATALOG_GOAL.md` §4.1 for a future DNA-virome delta card.
5. **What is deferred (`DEC-59 (3)`)**:
   - `palmDB` (`sotus.palmprint.faa`, `unique.palmprint.faa`): protein-only (`~100 aa` barcodes, no nucleotide contigs) — deferred under `DEC-60 (2)`.
   - `LucaProt` non-representative contigs (`~351,155` of `513,134`, *derived*) and `RVMT` non-representative / non-RdRp contigs (`~2.2M+` of `2.5M`, *derived*): recorded in `left_for_later.tsv.gz` and `CATALOG_GOAL.md` §4.1, recoverable from the same archived tarballs with zero re-fetch.
   - `VIRE` non-representative `Riboviria` members (`4,576` of `8,056`, *derived*: `8,056 - 3,480`): already tracked in `V3_left_for_later.tsv.gz`.

---

## 3. Mapping RNA Viruses, Obelisks, and cccRNAs onto `bio-box`'s Data Model

### 3.1 `U -> T` Alphabet Normalization and `contig_uid`

From the code:
- `src/biobox/fasta.py` preserves exact source bytes; `ids.contig_uid(seq)` hashes `seq.upper()` (`src/biobox/ids.py:10-12`); `store.py:91` defines `_COMP = bytes.maketrans(b"ACGTRYKMBVDHNacgtrykmbvdhn", b"TGCAYRMKVBHDNtgcayrmkvbhdn")` (no `U`/`u`).
- `pyhmmer.easel.Alphabet.dna()` (`src/biobox/search/orfs.py:38`) and `pyrodigal`/`pyrodigal-gv` require DNA alphabet (`T`, not `U`).

**Design rule**: Normalize `U -> T` and `u -> t` **at input staging** (before L0 ingest, in the stage-prep script), and record every modified accession and count of replaced bases in `u_to_t_normalized.tsv` beside `trims.tsv`.
- *Why at staging instead of inside `store.py`*: (a) preserves L0-seq exact-byte round-trip (`raw_digest`), (b) ensures an RNA record deposited with `U` (e.g. a ViroidDB or supplementary FASTA) hashes to the **exact same `contig_uid`** as its cDNA twin in RefSeq/RVMT/LucaProt, and (c) requires zero changes to `store.py`, `ids.py`, `orfs.py`, or `pyrodigal`.

### 3.2 Genome Kinds (`rna_virus`, `cccrna`)

Currently `genome_input.enums["kind"]` and `genome_dim.enums["kind"]` allow `("prokaryote", "virus", "plasmid")` (`src/biobox/registry/schemas.py:820, 909`), and `genecall.KINDS` matches (`src/biobox/l1/genecall.py:169`).

**Design rule**: Add two enum values to `genome_input.kind`, `genome_dim.kind`, and `genecall.KINDS`:
- **`rna_virus`**: Linear or circular RdRp-encoding (or reverse-transcribing *Pararnavirae*) RNA viruses (`RefSeq Riboviria`, `LucaProt`, `RVMT`, `Serratus`, `VIRE Riboviria`, `ELGV` RNA).
- **`cccrna`**: Subviral and hybrid covalently closed circular RNAs (`Obelisks`, `ViroidDB` viroids/satRNAs/deltaviruses, `Lee/Neri 2023` viroid-like cccRNAs, and optionally *Ambiviricota* — or *Ambiviricota* as `rna_virus` with `is_circular=True` + `ribozyme-1`; see PI Question 1 in §8).

*Schema compatibility (`LY-09`, `AGENTS.md` §3.8)*: Widening an `enums` tuple does not add, remove, rename, retype, or reorder any Parquet column. `schema_lock.json` is unchanged; all existing stages and releases remain byte-identical and continue to verify.

### 3.3 Segmented Genomes

`bio-box` already models multi-contig genomes natively (`ids.genome_uid(sorted(contig_uids))`, `genome_input` many-to-many under `DEC-27`):
- **RefSeq segmented RNA viruses** (e.g., *Orthomyxoviridae* 8 segments, *Spinareoviridae*/*Sedoreoviridae* 10–12 segments, *Bunyavirales* 3 segments, *Partitiviridae* 2 segments): Group segments sharing a RefSeq assembly / BioProject / isolate taxid+strain into **one `genome_uid`** with `n_contigs = k` segments, preserving each segment as its own `contig_uid` and recording segment name (`RNA1`, `L`/`M`/`S`, `PB2`, etc.) in the stage metadata sidecar (`rna1_metadata.tsv`).
- **Metatranscriptomic RdRp contigs** (`LucaProt`, `RVMT`, `Serratus`, `VIRE`, `ELGV`): Unlinked contigs are staged as 1-contig `genome_uid`s (identical to `V1`/`V2`/`MV1`), with a sidecar recording supergroup/phylum, RdRp completeness, and circularity.
- **No `5,000 bp` floor**: Explicitly exempt under `DEC-60 (1)` and `DEC-69` (most RNA virus segments and all Obelisks/viroids are `0.2–5 kb`).

### 3.4 Circular RNA Topology, DTR Trimming (`DEC-65`), and Origin-Wrapping ORFs

- **DTR trimming (`DEC-65`)**:
  - Metatranscriptome-assembled circular RNAs (`Obelisks`, `Lee/Neri 2023`, `Ambiviruses`, circular `LucaProt`/`RVMT` contigs) follow `DEC-65 (1)`: exact `20..1,000` bp direct terminal repeats (not low-complexity) have one trailing copy trimmed (logged in `trims.tsv`) and are staged `is_circular=True`; unit-length monomers already trimmed by `VNom`/`vdsearch` (`none` repeat) are staged `is_circular=True` untrimmed.
  - Public curated records (`RefSeq`, `ViroidDB`, GenBank `Ambiviricota`) follow `DEC-65 (2)`: never trimmed, byte-identical as deposited, `is_circular=True` where flagged circular.
- **Origin-wrapping coordinates (`AGENTS.md` §3.1)**:
  - Both `genecall.py:321-328` (`WRAP_FLANK = 50_000`, `flank = min(n, WRAP_FLANK)`) and `orfs.py:114-134` (`text = base * 3`, `start in [0, L)`, `end > L`) already emit `0 <= start < contig_len < end` for origin-wrapping ORFs, and `ncrna.py:552-559` (`_extended(seq, flank)`) already emits `end > contig_len` for origin-wrapping ribozymes.
  - **Stopless circular reading frames (edge case)**: In `orfs.py:123`, `if open_l or open_r or not 3 <= lo < L + 3 or hi - lo > L: continue` silently drops any circular frame that has **zero stop codons around the entire circle** (rolling-circle / super-unit-length translation, observed in a minority of small `~0.3–1.2 kb` cccRNAs and *Avsunviroidae*-like circles). Similarly, `genecall.py:325` drops `partial_begin or partial_end` on circular contigs. **Design rule**: Explicitly flag stopless circular frames in the stage QC report; do **not** emit synthetic infinite/multi-turn proteins into `protein_seq` (which would violate `end <= start + contig_len` in `pack.py:729`), and optionally emit a single unit-length `[0, L)` `dark_orf` flagged `partial=True` if benchmark `C-50a` shows biologically real stopless Oblins/HDAgs.

---

## 4. Gene Calling & Ribozyme Annotation

### 4.1 Why `pyrodigal` / `pyrodigal-gv` Alone Is Insufficient for RNA Genomes

`pyrodigal` and `pyrodigal-gv` (`src/biobox/l1/genecall.py`) are trained on prokaryotic/phage Shine-Dalgarno RBS motifs, prokaryotic start codons (`ATG`/`GTG`/`TTG`), and a **strictly non-overlapping** dynamic programming path. On RNA genomes this causes four systematic failure modes:

1. **Eukaryotic polyprotein truncation & spurious internal starts**: ~85%+ of *Riboviria* infect eukaryotes (no Shine-Dalgarno RBS; translation via Kozak scanning, IRES, or cap-independent elements). `pyrodigal-gv` frequently splits or truncates long polyproteins (*Picornavirales*, *Flaviviridae*, *Nidovirales* up to `47.3 kb`) at internal `GTG`/`TTG` codons with accidental RBS-like upstream purine runs.
2. **Overlapping and ambisense ORFs dropped**:
   - **Obelisks (`~1 kb`)**: `Oblin-1` (`~200 aa`) and `Oblin-2` (`~50–60 aa`) frequently overlap on opposite strands or alternative frames.
   - **Ambiviruses (`4–5 kb`)**: Two large **ambisense** ORFs (`ORF-A` RdRp on `+`, `ORF-B` unknown on `-`), sometimes overlapping near termini.
   - **Ribozyviruses / Deltaviruses (`~1.6–1.7 kb`)**: `HDAg` (`~195–214 aa`) sits on the antigenomic strand and wraps the circular origin in many isolates.
   - **Overlapping viral genes** (*Coronaviridae*, *Orthomyxoviridae*, *Tymovirales*, *Leviviricetes* lysis genes): Prodigal's non-overlapping DP discards the nested/overlapping ORF.
3. **Genetic codes**: Eukaryotic RNA viruses use **standard code 1** (`ATG` start only, `TAA/TAG/TGA` stop); prokaryotic RNA phages (*Leviviricetes*, *Cystoviridae*) use **code 11**; fungal/protist mitochondrial RNA viruses (*Mitoviridae*, *Narnaviridae*) use **code 4** (`TGA = Trp`); ciliate RNA viruses use **code 6** (`TAA/TAG = Gln`). `pyrodigal-gv` tests codes `11, 4, 15` (not code `1`).
4. **Spurious short calls on non-coding viroids (`0.25–0.45 kb`)**: Classical plant viroids (*Pospiviroidae*, *Avsunviroidae*) encode **no proteins**; a bacterial meta-mode caller or unfiltered 6-frame caller with low `min_aa` generates false short CDSs.

### 4.2 Proposed Gene-Calling Architecture (No New Container Binary Required)

`orfipy` is **not** in `container/environment.yml`. However, `src/biobox/search/orfs.py` **already implements** deterministic, circular-aware 6-frame stop-to-stop ORF extraction (`six_frame_orfs`) via pinned `pyhmmer.easel` (supporting codes `1–6, 9–14, 15, 16, 21–25`).

We propose evaluating two clean calling policies in pilot card **`C-50a`** against RefSeq curated CDSs, RVMT *Leviviricetes*/*Mitoviridae*, Zheludev Oblin-1/Oblin-2, Forgia Ambiviruses, and ViroidDB (negative control):

| Mode | Target `kind` | Algorithm | Genetic Code | Polyproteins & Frameshifts | Overlapping / Ambisense ORFs |
|---|---|---|---|---|---|
| **`rna-orf6`** *(recommended default for `rna_virus`)* | `rna_virus` | Deterministic 6-frame stop-to-stop (or `ATG`-anchored) ORF caller on `orfs.six_frame_orfs` with **shadow-ORF suppression**: keep every ORF `>= min_aa` (`60 aa` default, `>= 100 aa` for partial contig-edge ORFs) that is not `>= 80%` antisense-shadowed by a strictly longer ORF (unless its length is `>= 100 aa`, preserving ambisense/overlapping viral genes like `ORF1a`/`ORF1b`, *Ambivirus* `ORF-B`, and nested movement/capsid proteins). | Per-contig best code among **`{1, 4, 11}`** (plus `6` if ciliate flag / `TAA`/`TAG` depletion), chosen by longest-ORF coding fraction (same principle as `RVMT`/`LucaProt` `transeq` + code-4 mitovirus rescue, with `C-50a` calibrating the code-4 vs code-1 LLR threshold just like `B-69`'s code-25 LLR). | **Polyproteins kept intact** as 1 CDS slot (domains split by `l1_domains` `env_start..env_end` & `domain_family_tokens`). **Programmed frameshifts** (`-1`/`+1` PRF) emit both `ORF1a` and `ORF1b` as adjacent overlapping CDS slots (`edit_type="none"`). | **Preserved** (both strands and shifted frames above shadow threshold). |
| **`cccrna-orf6`** *(recommended default for `cccrna`)* | `cccrna` | Circular-aware 6-frame `ATG`-to-stop (or stop-to-stop fallback) caller with `min_aa = 50` (capturing `Oblin-2` at `~50–60 aa`, `Oblin-1` at `~200 aa`, `HDAg` at `~195 aa`, while yielding **0 CDS** on `250–450 nt` non-coding viroids whose longest `ATG`-to-stop ORF is typically `< 45 aa`). | Code **`1`** (standard eukaryotic/viral; code `11` for bacterial-host Obelisks is identical on stop codons and differs only in non-`ATG` starts, which `ATG`-anchoring avoids). | N/A | **Preserved** (`Oblin-1` + `Oblin-2`, `HDAg` antigenome, *Ambivirus* ambisense). |
| **`gv-meta`** *(baseline in `C-50a`)* | `rna_virus` (*Leviviricetes* comparison) | Existing `pyrodigal_gv.ViralGeneFinder(meta=True)`. | `11 / 4 / 15` | Often splits polyproteins; misses code-1 eukaryotic starts. | Drops overlapping/ambisense ORFs. |

*Producer naming (`AGENTS.md` §3.10, `DEC-18`)*: Because `l1/feature`, `l1/feature_protein`, `l1/protein_occurrence`, and `l0seq/protein_seq` are partitioned by `producer` / `gene_call_version`, an RNA caller can either:
- **(Option A — Recommended)**: Run under a new `genecall_mode` (`rna-orf6` / `cccrna-orf6`) recorded in `genecall_plan` while keeping `GENE_CALL_VERSION` unchanged (exactly as `DEC-71` / card `B-67` did for `gv-meta` and `B-69` for `gv-meta-c25`), so **one release-wide `gene_call_version` and one `tape_version` span DNA + RNA stages** with zero changes to `pack.py` or `catalog.py`; OR
- **(Option B)**: Use a distinct `gene_call_version` (`biobox-rna-orf6-v1`), which would require multi-`gene_call_version` support in `pack.py`, `tapes.py`, and `cluster.py`. **Option A is strongly preferred** because `B-67`/`B-69` already established the pattern and `feature_protein.genetic_code` + `genecall_plan` preserve full provenance per stage.

### 4.3 Ribozyme & Viroid ncRNA Annotation (`Infernal` + `Rfam 15.1` `ribozyme-1`)

Checking `src/biobox/l1/rfam_sets/prok-2.tsv` against `planning/design/0-08_rfam15.1_models.tsv` (*measured* statically):
- `prok-2` (`770` models, `DEC-30`) contains **only 3 ribozyme/hairpin models** (`RF02276` `Hammerhead_II`, `RF02678` `Hatchet`, `RF03059` `raiA-hairpin`), and **all 3 have `whole_genome=0`** (intergenic regions only — missing ribozymes inside or overlapping ORFs on compact circular RNAs!).
- Meanwhile, **`Rfam 15.1` (already pinned in `refdb.yaml`, `sha256 e2e636e9...`)** already contains **22+ self-cleaving ribozyme and viroid/RNA-virus structural CMs**:
  - **Hammerhead (6)**: `RF00008` (`Hammerhead_3`), `RF00163` (`Hammerhead_1`), `RF02275` (`Hammerhead_HH9`), `RF02276` (`Hammerhead_II`), `RF02277` (`Hammerhead_HH10`), `RF03152` (`RAGATH-1-hammerhead`)
  - **HDV / Deltavirus-like (5)**: `RF00094` (`HDV_ribozyme`), `RF02682` (`HDV-F-prausnitzii`), `RF00622` (`CPEB3_ribozyme`), `RF01787` (`drz-agam-1`), `RF01788` (`drz-agam-2-2`)
  - **Twister & Twister-sister (4)**: `RF03160` (`twister-P1`), `RF03154` (`twister-P3`), `RF02684` (`Twister-P5`), `RF02681` (`Twister-sister`)
  - **Hairpin (3)**: `RF00173` (`Hairpin`), `RF04190` (`Hairpin-meta1`), `RF04191` (`Hairpin-meta2`)
  - **Hatchet & Pistol (2)**: `RF02678` (`Hatchet`), `RF02679` (`Pistol`)
  - **Viroid & other catalytic/viral (2+)**: `RF00362` (`Pospi_RY`), `RF04188` (`Hovlinc`) (plus optional `RF00173`/`RF00502`/`RF00617`/`RF02359` audited in `C-50a`).

**Design rule**:
1. Freeze a new `rfam_sets/ribozyme-1.tsv` model set (~22 models from `Rfam 15.1`, **all with `whole_genome=1`**) registered in `ncrna.FROZEN_SETS`.
2. `ncrna.py` ALREADY handles circular origin-wrapping (`_extended(seq, flank_whole)`), fixed `Z_nt = 1e5 Mb` (`DEC-04`), `--cut_ga`, and `Rfam-15.1.clanin` clan competition (`DEC-28 (4)`).
3. **Compute cost**: `22` short CMs (`CLEN 45–258`) on `~0.3 Gbp` (Phase 1) or `~3 Gbp` (Phase 1+2) is **`< 20` CPU-minutes for Phase 1** and **`~2–3` CPU-hours for Phase 2** (*estimate*, scaled from `prok-2`'s `300` CPU-hours / `10 Gbp` for `770` models, `DEC-30`).
4. **Known sensitivity gap**: Divergent Obelisk hammerhead variants (Zheludev 2024) and custom `vdsearch` / `Ambivirus` ribozyme CMs can miss Rfam 15.1 GA cutoffs. Adding a second custom CM library (`vdsearch_ribozymes` / `obelisk_hhr`) in `refdb.yaml` is a clean, additive follow-up once `ribozyme-1` is measured in `C-50a`.

---

## 5. Tokens, Domain Libraries, `pack build` (`dark_orf`), and Neighborhoods

### 5.1 Domain Libraries & Family Tokens (`Pfam`, `NeoRdRp2`, `Oblin`)

- **Pfam 38.2** (`30,134` models, already pinned): Covers canonical eukaryotic/prokaryotic viral RdRps, helicases, proteases, capsids, and `HDAg` (`PF01205`), at trivial cost on Phase 1 (`~0.3 Gbp` = `~450` CPU-minutes, *derived* from `1,550` CPU-min/Gbp in `C-38`). However, Pfam misses a large fraction of divergent metatranscriptomic RdRps (`LucaProt` 180 supergroups, `RVMT` novel phyla) and has **0 models for `Oblin-1`/`Oblin-2`** or *Ambivirus* `ORF-B`.
- **Specialized RNA-virus HMM library (additive, zero schema change)**:
  - **`NeoRdRp2`** (`19,394` HMMs, MIT, `redistributable: true`, Zenodo `10.5281/zenodo.10807561` / `10.5281/zenodo.10583247`, Sakaguchi et al. 2024 *Frontiers in Virology*) — comprehensive RdRp HMM library built from `557,197` seed RdRps (*published*).
  - **Custom `oblin-1` / `oblin-2` HMMs**: Built deterministically (`hmmbuild --seed 42`, matching `W-40` `phrogs_rnap` in `refdb.yaml`) from Zheludev et al. 2024 SDR alignments (`wb363nt3637`).
  - Because `l1_domains` is partitioned by `(library, stage)` and `tokens.py:170` explicitly allows a library to be scanned on **only a subset of stages** (`"a library need not be scanned on every stage"`), scanning `neordrp-2.1` or `oblin-1` **only on RNA/cccRNA stages** incurs **zero compute on DNA stages** and automatically populates `token_dict`, `tape_tokens`, and `token_postings`!
- **Ribozyme (`ncrna`) tokens**: Currently `l3/tokens.py:43` assigns `null` tokens to non-CDS slots (`kind != 1`). Extending `tokens.py` (optional follow-up card `C-50e`) so `kind=2` (`ncrna`) slots emit `family` tokens from `l1_ncrna.model_acc` (with `source` = `rfam-15.1-ribozyme1`) would make `"Hammerhead_3 next to Oblin-1"` or `"HDV_ribozyme next to RdRp"` a native O(1) `token_postings` pair query.

### 5.2 `pack build` (`unique`, `igr`, `dark_orf`) Behavior on `0.2–10 kb` RNA Genomes

From `src/biobox/search/pack.py` and `src/biobox/search/orfs.py`:
1. **`unique` proteins**: All RNA/cccRNA proteins enter `unique` and release `NR90`/`NR50` clustering automatically (`MMseqs2` clusters proteins regardless of genome kind).
2. **`igr` target (`IGR_MIN_NT = 30`)**: Non-coding viroids (`0` CDS, `250–450 nt`) have their entire circular genome emitted as one `igr` record (`igrs_of_contig` returns `(0, n)` for a 0-gene contig), making them directly searchable by nucleotide/CM queries against `--target igr`.
3. **`dark_orf` target (`END_BP = 3,000`, `MIN_AA = 30`, `MAX_GENE_OVERLAP = 0.5`) — CRITICAL INTERACTION**:
   - In `orfs.py:180-182`, any ORF within `3,000 bp` of a linear contig end (`at_end = not circular and (o.start < 3000 or o.end > L - 3000)`) **bypasses `MAX_GENE_OVERLAP`**.
   - On linear RNA virus contigs `<= 6,000 bp` (the majority of `LucaProt`/`RVMT`/`Serratus`), **every 6-frame ORF `>= 30 aa` in all 6 frames** (even 100% overlapping a called polyprotein) is kept in `dark_orf`!
   - *Impact*: On `~0.3 Gbp` (Phase 1) this is negligible (`~1–2M` dark ORFs, `< 0.3 GB` in the pack, *estimate*); on `~3 Gbp` (Phase 2) it adds `~15–25M` dark ORFs (`~1.5–2.5 GB` in the pack, *estimate*) — acceptable and actually **useful** for recovering missed overlapping/frameshift ORFs via `bio-box search --target dark_orf`, with no code change required.
   - *0-gene contig code check*: Verify in `C-50b` that `pack.py`'s `_contig_codes` defaults cleanly (code `1` or `11`) for 0-CDS viroids so `pack build` never fails on a contig with zero `feature_protein` rows.

### 5.3 What Neighborhood Queries (`nbr`) Mean on `1–10 kb` Genomes

In `src/biobox/search/nbr.py`, `Policy(window="10g", edge="truncate", circular="auto")` behaves as follows on `1–10 kb` RNA genomes:
- **Linear `1–10 kb` contigs (1–4 genes)**: `window=10g` returns the **entire contig**, with `FLAG_CENSORED` (bit 0) set on every row. `edge=exclude:3kb` or `edge=require:3` would drop **>90% of RNA virus contigs** (*estimate*) and should **not** be used on RNA queries.
- **Circular `0.3–5 kb` contigs (`Obelisks`, `Ambiviruses`, `Ribozyviria`)**: `circular="auto"` wraps around the circle with zero censoring, returning the full circular gene + ribozyme layout.
- **Three biological neighborhood modes for RNA genomes**:
  1. **Intra-polyprotein domain architecture**: Answered by `l1_domains` + `tape_tokens.domain_family_tokens` + `fusions_of` (e.g., `Protease + Helicase + RdRp` order on one CDS slot).
  2. **Whole-contig gene + ribozyme synteny**: Answered by `nbr --window 10g --edge truncate` on a tape that includes both `cds` and `ribozyme-1` `ncrna` slots (e.g., `Oblin-1 <-> Oblin-2 <-> Hammerhead_3`, or *Ambivirus* `RdRp <-> ORF-B <-> HHR/HDV`).
  3. **Cross-segment co-occurrence**: For segmented RefSeq genomes (`RNA1`), segments are separate contigs sharing one `genome_uid`. `nbr` stays per-contig; cross-segment partner queries join on `genome_input.genome_uid` in DuckDB.

---

## 6. Recommended Phased Staging Order

| Phase / Stage | Sources Included | `kind` & Topology | Depth & Annotations | Est. Contigs / Size | Est. Build Cost |
|---|---|---|---|---|---|
| **Phase 0** (`C-50a`, `C-50b`) | Benchmark sample (`~500` RefSeq + `500` RVMT + `500` Obelisks + `500` ViroidDB/Ambiviruses) | Test only (scratch) | Compare `gv-meta` vs `rna-orf6` vs `cccrna-orf6` + `ribozyme-1` | `~2,000` contigs (`~10 Mbp`, *estimate*) | `< 15` CPU-min (*estimate*) |
| **Phase 1a: `RNA1`** | NCBI RefSeq `Riboviria` + `VIRE` `Riboviria` vOTU reps (`3,480`) + `ELGV` RNA reps | `rna_virus` (segmented `genome_uid` for RefSeq; linear/circular per source) | **`annotated`**: Pfam 38.2 + `ribozyme-1` (+ `NeoRdRp2` if approved) | `~15k–20k` contigs; **`~0.15–0.25 Gbp`** (*estimate*) | `~400–600` CPU-min (*estimate*) |
| **Phase 1b: `OBL1`** | Zheludev et al. 2024 Obelisks (`29,959` clusters) | `cccrna`, `is_circular=True` (`DEC-65` DTR trim) | **`annotated`**: Pfam 38.2 + `ribozyme-1` + `oblin` HMMs | `29,959` contigs (*published*); **`~30 Mbp`** (*derived*) | `~60–100` CPU-min (*estimate*) |
| **Phase 1c: `VIRD1`** | `ViroidDB` (`~9,891`) + *Ambiviricota* (`~500`) + Lee/Neri 2023 cccRNAs (`11,378`), deduped in that precedence order | `cccrna`, `is_circular=True` (`DEC-65` rules) | **`annotated`**: Pfam 38.2 + `ribozyme-1` | `~18k–21k` net contigs; **`~15–20 Mbp`** (*estimate*) | `~40–80` CPU-min (*estimate*) |
| **Phase 2a: `RNA2`** | `LucaProt` `161,979` species representatives (`~351k` non-reps deferred to `CATALOG_GOAL.md` §4.1) | `rna_virus` (1-contig `genome_uid`s) | **`called`** (or **`annotated`** with Pfam 38.2 + `NeoRdRp2` — only `~0.5 Gbp`!) | `~162k` contigs; **`~0.4–0.6 Gbp`** (*estimate*) | `~100` CPU-min (`called`) or `~1,000` (`annotated`, *estimate*) |
| **Phase 2b: `RNA3`** | `RVMT` v3 RdRp-bearing reps (`~125k–330k`) + `Serratus` `s3://serratus-public/` nucleotide contigs (deduped vs `RNA1`+`RNA2`) | `rna_virus` (1-contig `genome_uid`s) | **`called`** (hit-window Pfam via `bio-box annotate round`, `DEC-54`) | `~300k–600k` net contigs; **`~1.0–2.5 Gbp`** (*estimate*) | `~250–600` CPU-min (*estimate*) |

*Note on CRISPR*: `minced` (`--crispr`) is **skipped** on `rna_virus` and `cccrna` stages (RNA viruses and viroids do not carry prokaryotic CRISPR spacer arrays).

---

## 7. Draft `sources.yaml` Entries (For Review — Not Yet Added to `sources.yaml`)

These entries follow the exact schema of `sources.yaml` (`DEC-50`, `DEC-51`, `AGENTS.md` §8). Exact file names, byte sizes, and checksums marked `<verify on fetch>` are populated by the staging cards (`C-50c`..`C-50f`) before downloading.

```yaml
  - name: ncbi_refseq_riboviria
    what: NCBI RefSeq viral genomic FASTA and GenBank flatfiles, filtered to realm Riboviria; stage RNA1
    url: https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.{k}.1.genomic.fna.gz
    url_gbff: https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.{k}.genomic.gbff.gz
    checksum: RefSeq release catalog / MD5SUMS in https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/
    records: "~8,600-11,000+ Riboviria records across curated complete/segmented genomes (published)"
    license: public domain (NCBI / INSDC)
    redistributable: true
    limits: one connection; at most 3 requests per second; stop after 5 refusals in a row
    approved: pending PI decision on C-50

  - name: lucaprot_rna_viruses
    what: LucaProt global RNA virosphere (Hou, He et al. 2024 Cell), 513,134 contigs / 161,979 species reps; stage RNA2
    url: https://ndownloader.figshare.com/files/{file_id}
    record: https://doi.org/10.6084/m9.figshare.26298802
    files:
      - {file: "<verify exact species/contig FASTA filenames on Figshare API>", bytes: 0, md5: "<verify>"}
    records: "513,134 RNA viral contigs; 161,979 species reps at >90% RdRp ID (published)"
    license: CC BY 4.0
    redistributable: true
    limits: one connection; stop after 5 refusals in a row
    approved: pending PI decision on C-50

  - name: rvmt_v3
    what: RNA Viruses in MetaTranscriptomes v3 (Neri et al. 2022 Cell), >2.5M contigs / ~330k RdRps; stage RNA3
    url: https://zenodo.org/records/7368133/files/{file}?download=1
    record: https://doi.org/10.5281/zenodo.7368133
    files:
      - {file: RVMT_Zenodo_v3.tar.gz, size: "9.5 GB (published)", md5: "<verify on Zenodo API>"}
    records: ">2,500,000 RNA virus contigs; ~330,000 RdRps (published)"
    license: "<verify exact Zenodo license badge on fetch; code MIT>"
    redistributable: true   # update if Zenodo badge restricts redistribution (DEC-51)
    limits: one connection; stop after 5 refusals in a row
    approved: pending PI decision on C-50

  - name: serratus_public_contigs
    what: Serratus petabase-scale SRA RdRp nucleotide contig assemblies (Edgar et al. 2022 Nature); stage RNA3 (palmDB protein-only deferred under DEC-60)
    url: https://serratus-public.s3.amazonaws.com/{key}
    record: https://serratus.io and https://github.com/ababaian/serratus
    records: "131,957 novel RNA virus sOTUs (published); nucleotide contig subset deduplicated against RNA1+RNA2+RVMT"
    license: CC0 / INSDC (open-access SRA derivatives)
    redistributable: true
    limits: one connection; stop after 5 refusals in a row
    approved: pending PI decision on C-50

  - name: obelisks_zheludev2024
    what: Global metatranscriptomic Obelisks (Zheludev et al. 2024 Cell), 29,959 circular RNA clusters and Oblin alignments; stage OBL1
    url: https://stacks.stanford.edu/file/druid:wb363nt3637/{file}
    record: https://purl.stanford.edu/wb363nt3637  # DOI 10.25740/wb363nt3637
    files:
      - {file: "<verify exact sequence/table filenames on SDR record>", bytes: 0, md5: "<verify>"}
    records: "29,959 Obelisk clusters at 90% nt ID, ~1 kb circular RNAs (published)"
    license: "<verify exact SDR license metadata on fetch>"
    redistributable: true   # set false if SDR specifies CC BY-NC (DEC-51)
    limits: one connection; stop after 5 refusals in a row
    approved: pending PI decision on C-50

  - name: viroid_cccrna_lee2023
    what: Metatranscriptomic viroid-like cccRNAs (Lee, Neri et al. 2023 Cell), 11,378 ribozyme-bearing circular RNAs; stage VIRD1
    url: https://zenodo.org/records/6859104/files/{file}?download=1
    record: https://doi.org/10.5281/zenodo.6859104
    files:
      - {file: "<verify exact FASTA/archive filenames on Zenodo API>", bytes: 0, md5: "<verify>"}
    records: "11,378 viroid-like cccRNAs (published)"
    license: "<verify exact Zenodo license badge on fetch; vdsearch code MIT>"
    redistributable: true   # update if Zenodo badge restricts redistribution (DEC-51)
    limits: one connection; stop after 5 refusals in a row
    approved: pending PI decision on C-50

  - name: ambiviricota_forgia2023
    what: Fungal ambiviruses / Ambiviricota (Forgia et al. 2023 Nat Commun + NCBI GenBank Ambiviricota), 4-5 kb circular ssRNA viruses with RdRp + ribozymes; stage VIRD1
    url: https://doi.org/10.1038/s41467-023-38301-2  # Supplementary Data + NCBI E-utilities txid Ambiviricota
    records: "378 RdRp sOTUs (published)"
    license: CC BY 4.0 (Nat Commun) / public domain (INSDC)
    redistributable: true
    limits: one connection; at most 3 requests per second on NCBI E-utilities
    approved: pending PI decision on C-50

  - name: viroiddb
    what: ViroidDB curated viroids, deltaviruses (Ribozyviria), and circular satRNAs (Lee et al. viroids.org); stage VIRD1
    url: https://viroids.org/db/latest/{file}
    record: https://viroids.org and https://doi.org/10.5281/zenodo.5202945
    files:
      - {file: all.fasta, bytes: 0, sha256: "<record snapshot date + sha256 on fetch; mirror to tower per DEC-11>"}
      - {file: all.json, bytes: 0, sha256: "<record snapshot date + sha256 on fetch>"}
    records: "~9,891 sequences (9,354 viroids, 453 deltaviruses, ~35 satRNAs; published)"
    license: MIT / INSDC
    redistributable: true
    limits: one connection; stop after 5 refusals in a row
    approved: pending PI decision on C-50

  - name: elgv_rna_fraction
    what: Early-Life Gut Virome (Zeng, Almeida et al. 2024 Nat Commun), Riboviria subset of 82,141 species reps; stage RNA1 (DNA fraction deferred)
    url: https://ndownloader.figshare.com/files/{file_id}
    record: https://doi.org/10.6084/m9.figshare.21901557
    files:
      - {file: ELGV_representatives_82141.fna.gz, size: "378.3 MB (published)", md5: "<verify on Figshare API>"}
    records: "82,141 species reps total (160,478 all, published); Riboviria fraction extracted by metadata"
    license: CC BY 4.0
    redistributable: true
    limits: one connection; stop after 5 refusals in a row
    approved: pending PI decision on C-50
```

*(Note: `vire_genomes` `Riboviria_genomes.tar` is already listed in `sources.yaml` under `DEC-50`; `NeoRdRp2` Zenodo `10.5281/zenodo.10807561` belongs in `refdb.yaml` if Question 4 is approved.)*

---

## 8. Questions for the PI & Follow-Up Cards

| # | Decision Question | Options | Recommendation |
|---|---|---|---|
| **Q1** | **Genome kinds (`schemas.py`)**: How should RNA viruses, Obelisks, Ambiviruses, and viroids be typed in `genome_input.kind` / `genome_dim.kind`? | (a) Add two kinds: `rna_virus` and `cccrna` (*Ambiviricota* in `cccrna` with `is_circular=True` or `rna_virus` with `is_circular=True`); (b) Reuse existing `virus` kind for all RNA viruses and add only `cccrna`. | **(a) Add `rna_virus` and `cccrna`** (with *Ambiviricota* as `cccrna`, `is_circular=True` so all ribozyme-bearing circular RNAs share one kind and pipeline; `rna_virus` stays cleanly separable from DNA phages/giant viruses in `genome_dim`). |
| **Q2** | **RNA & cccRNA gene calling**: Keep `pyrodigal-gv` (`gv-meta`) everywhere, or add deterministic 6-frame ORF modes (`rna-orf6`, `cccrna-orf6`) built on `biobox.search.orfs`? | (a) Benchmark `gv-meta` vs `rna-orf6` / `cccrna-orf6` on a `~2,000`-contig stratified sample (`C-50a`), then adopt `rna-orf6` + `cccrna-orf6` under the existing `GENE_CALL_VERSION` (following `B-67`/`B-69`); (b) Use `gv-meta` as-is (accepting polyprotein splits and missed Oblin-2 / ambisense ORFs); (c) Add a new container binary (`orfipy` / `VIGOR4`). | **(a) Run `C-50a` benchmark and adopt `rna-orf6` + `cccrna-orf6` on pinned `pyhmmer.easel`** — zero new container dependencies, preserves intact polyproteins, overlapping/ambisense ORFs, and 0-CDS viroids. |
| **Q3** | **Ribozyme annotation (`Infernal`)**: How should self-cleaving ribozymes be annotated on RNA/cccRNA stages? | (a) Freeze `rfam_sets/ribozyme-1.tsv` (~22 `Rfam 15.1` ribozyme/viroid CMs, all `whole_genome=1`) and run `l1/ncrna.py` on `RNA1`, `OBL1`, `VIRD1` (and optionally `RNA2`/`RNA3`); (b) Skip ribozyme CMs for now; (c) Also import custom `vdsearch`/Obelisk ribozyme CMs into `refdb.yaml` immediately. | **(a) Freeze `ribozyme-1` from pinned `Rfam 15.1` now** (`< 20` CPU-min for Phase 1, zero external CM downloads), and measure Obelisk/Ambivirus recall in `C-50a` to decide whether custom `vdsearch`/Obelisk CMs are needed later. |
| **Q4** | **Specialized HMM libraries (`refdb.yaml`)**: Scan RNA/cccRNA stages with `Pfam 38.2` only, or add `NeoRdRp2` (`19,394` HMMs, MIT) + `Oblin-1`/`Oblin-2` HMMs on RNA/cccRNA stages only? | (a) `Pfam 38.2` + `oblin` HMMs (tiny: 2 models) on Phase 1 (`RNA1`, `OBL1`, `VIRD1`), plus `NeoRdRp2` on `RNA1`–`RNA3` (scanned on RNA stages only — zero cost on DNA stages); (b) `Pfam 38.2` only at L1, keeping `NeoRdRp2`/`Oblin` as user search panels (`bio-box search`). | **(a) Add `oblin` (2 HMMs) + `NeoRdRp2` (`19,394` HMMs, MIT) for RNA/cccRNA stages only** — gives `token_dict` / `token_postings` real RdRp & Oblin family tokens where Pfam is blind, with zero compute on DNA stages. |
| **Q5** | **Ribozyme (`ncrna`) tokens on `token_postings`**: Should `l3/tokens.py` emit `family` tokens for `kind=2` (`ncrna`) slots? | (a) Keep `tokens.py` unchanged for Phase 1 (ribozymes live in `l1_ncrna` and on `contig_tapes` for `nbr` / SQL joins); (b) Extend `tokens.py` (`C-50e`) so `ncrna` slots carry `family` tokens from `l1_ncrna.model_acc`. | **(a) Defer `ncrna` token emission to follow-up card `C-50e`** so Phase 1 (`RNA1`, `OBL1`, `VIRD1`) ships without touching `tokens.py` / `postings.py` on the main DNA release train; revisit once `prok-2` (`DEC-42`) and `ribozyme-1` are both live. |
| **Q6** | **Phased staging & `Serratus` / `palmDB` scope**: Approve the `RNA1 -> OBL1 -> VIRD1` (Phase 1, `~0.3 Gbp`, `annotated`) followed by `RNA2` (`LucaProt` species reps) + `RNA3` (`RVMT` v3 + `Serratus` S3 nucleotide contigs, deferring protein-only `palmDB` per `DEC-60 (2)`) order? | (a) Approve Phase 1 (`RNA1`, `OBL1`, `VIRD1`) + Phase 2 (`RNA2`, `RNA3`) in that order, with `DEC-50`/`DEC-51` fetch rules; (b) Stage `LucaProt` (`RNA2`) ahead of `OBL1`/`VIRD1`. | **(a) Phase 1 (`RNA1`, `OBL1`, `VIRD1`) first, then Phase 2 (`RNA2`, `RNA3`)** — Phase 1 is `< 0.3 Gbp` total, validates both `rna_virus` and `cccrna` pipelines end-to-end in under 2 hours on the pod, and gives curated ground-truth ownership (`DEC-26`) before metatranscriptomic contigs arrive. |

### Follow-Up Cards (To Be Added to `TASKS.md` After PI Decision)

- **`C-50a` (M, Pod: yes)**: Stratified RNA/cccRNA gene-calling & ribozyme benchmark (`~2,000` contigs from RefSeq `Riboviria`, RVMT *Leviviricetes*/*Mitoviridae*, Zheludev Obelisks, Forgia *Ambiviricota*, and ViroidDB) comparing `gv-meta` vs `rna-orf6` / `cccrna-orf6` and measuring `Rfam 15.1` `ribozyme-1` recall.
- **`C-50b` (M, Pod: yes)**: Code enablement — widen `genome_input.kind` / `genome_dim.kind` (`rna_virus`, `cccrna`), add `rna-orf6` / `cccrna-orf6` modes to `genecall.py`, freeze `rfam_sets/ribozyme-1.tsv` in `ncrna.py`, verify 0-CDS viroid handling in `pack.py` (`dark_orf`/`igr`), and add `U -> T` staging helper with unit tests.
- **`C-50c` (M, Pod: yes)**: Fetch, verify (`sources.yaml` under `DEC-50`/`DEC-51`), stage, and build Phase 1 **`RNA1`**, **`OBL1`**, and **`VIRD1`** (`DEPTH=annotated`, Pfam 38.2 + `ribozyme-1` + optional `oblin`/`NeoRdRp2`), mirror to `tower`, and record `ELGV` DNA/non-rep remainders in `CATALOG_GOAL.md` §4.1.
- **`C-50d` (M, Pod: yes)**: Fetch, deduplicate against `RNA1`/`VIRD1`, stage, and build Phase 2 **`RNA2`** (`LucaProt` `161,979` species reps) and **`RNA3`** (`RVMT` v3 + `Serratus` S3 nucleotide contigs), recording all non-representative remainders in `CATALOG_GOAL.md` §4.1.
- **`C-50e` (M, optional)**: Emit `family` tokens in `l3/tokens.py` / `token_postings` for `kind=2` (`ncrna`) slots (`ribozyme-1` and `prok-2`), enabling direct `token_postings` pair queries between CDS domains and ribozymes/ncRNAs.
