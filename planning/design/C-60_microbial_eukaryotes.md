# C-60 design: microbial eukaryotes and giant-virus hosts (`DEC-69`)

Card `C-60` (tier **H**, design note; the PI decides from §8). `DEC-69` (7) brings unicellular/microbial eukaryotes and giant-virus (`GV1`, `EGV1`, `MIR1`) hosts into scope alongside `IMG-UMAG` eukaryote bins (`E1`/`E2`, `691` bins / `14.44 Gbp` already held on `/stor/work`). This note surveys every named source with primary-source-verified counts, URLs, licences, and overlaps (§2), designs the spliced-CDS data model under `AGENTS.md` §3–§4 (§3), evaluates eukaryotic gene-calling architectures within `wilkcomp02`'s CPU-only container bounds (§4), estimates storage and Pfam CPU budgets per `10 Gbp` (§5), proposes a staged ingest order with draft `sources.yaml` entries (§6–§7), and ends with six PI decision questions and follow-up implementation cards (§8–§9). It changes no code, schema, or `sources.yaml`.

This note reads the codebase on `task/c60-design` (`src/biobox/registry/schemas.py`, `src/biobox/ids.py`, `src/biobox/l1/genecall.py`, `src/biobox/l1/igr.py`, `src/biobox/l3/tapes.py`, `container/environment.yml`), `planning/DECISIONS.md` (`DEC-04` through `DEC-75`), `planning/CATALOG_GOAL.md`, `planning/tasks/TASKS.md` (`C-60`, `C-68`), and `planning/measurements/2026-11-19_jgi_mags_iuig_dedup.md`. No pod compute or local Python execution was run for this note.

Labels (`AGENTS.md` §6): *measured* (host, date, source note), *derived* (formula and inputs), *estimate* (explicit basis), *published* (primary citation or repository record).

---

## 0. The problem in ten lines

1. Today `genome_input.kind` and `genome_dim.kind` allow only `("prokaryote", "virus", "plasmid")` (`src/biobox/registry/schemas.py:820,909`), blocking even `--depth stored` ingest of the `691` `IMG-UMAG` eukaryote bins already on `/stor/work`.
2. Running `pyrodigal` on eukaryotic genomes shatters multi-exon CDSs at spliceosomal introns and fills `protein_seq` with intronic junk: JGI's prokaryotic GFFs on `IMG-UMAG`'s `691` eukaryote bins have `0` introns/exons, `366 bp` median CDS length (`122 aa`), and apparent coding density down to `39.8%` (*measured*, `2026-11-19_jgi_mags_iuig_dedup.md` §2.1).
3. The same failure hits intron-rich eukaryotic dsDNA viruses (`MIR1` *Mirusviricota* / *Okeanovirales*, whose nuclear-replicating genes carry dense spliceosomal introns and even trans-spliced terminases; *published*, Medvedeva et al. 2025).
4. `bio-box`'s core invariant (`DEC-18`, `AGENTS.md` §3.3) keys every feature by `feature_uid = blake3(contig_uid, start, end, strand, kind)` with a 1:1 `feature_protein` mapping (`PK: (feature_uid, gene_call_version)`), one slot per feature on `contig_tapes`, and `igr.py` defining IGRs as gaps between `feature` spans.
5. Representing each exon as a separate `feature` would blow up `contig_tapes` (`±10`-gene windows reduced to 1–2 genes), break 1:1 `feature_protein`, inflate `token_postings` 5–10×, and dump millions of spliceosomal introns into `igr` and `dark_orf` search packs.
6. Representing each spliced CDS as **one gene-span `feature`** (`start = min(exon_start)`, `end = max(exon_end)`, `edit_type = "intron"` when `n_exons > 1`) plus an optional **`l1_exons`** table preserves `igr.py`, `tapes.py`, `protein_seq`, `l1_domains`, `l2/cluster`, and `token_postings` with **zero** downstream schema breaks.
7. `bio-box.sif` (`container/environment.yml`) is CPU-only (`wilkcomp02`: `112` logical threads, `40`/`80` cap under `DEC-74`, `750G` RAM scope under `DEC-73`, `<= 280` CPU-min/process); `wilkcomp03` GPU is not approved (`DEC-16`), making GPU-dependent `Helixer` impractical for bulk CPU-only builds.
8. High-value reference hosts (`Parabasalia`, RefSeq `protozoa`/`fungi`, JGI `PhycoCosm`/`MycoCosm`, *Mirusviricota*) already ship curated spliced GFF3/GBFF gene models that can enter via a deterministic **`gff-import`** producer (re-translated and verified against `SeqStore`), while unannotated eMAGs (`E1`/`E2`, `SMAGs`, `TOPAZ`, `OMD`, `Tagirdzhanova`) need a CPU-only spliced caller (`miniprot` / `MetaEuk`).
9. Meanwhile, card `C-68` (`ET1`: `MMETSP` + `1KP` spliced cDNA transcriptomes) bypasses introns entirely for immediate eukaryotic protein/domain search, and its translated proteins can feed `miniprot`/`MetaEuk` as a pinned reference database alongside `EukProt v3`.
10. Because microbial eukaryotes have lower gene density than prokaryotes (`~3.9M` vs `~9.3M` genes per `10 Gbp`, *derived*), their Pfam 38.2 cost (`~190–215` CPU-h / `10 Gbp`) and `called`-depth storage (`~5.5–6.5 GB` / `10 Gbp`) are **smaller per Gbp** than `SPX1`, making compact host and reference stages affordable at full `--depth annotated`.

---

## 1. Baseline numbers used in this note

| # | Quantity | Value | Label and source |
|---|---|---|---|
| N1 | `depth="stored"` storage per `10 Gbp` | `3.1–4.0 GB` (`~3.3 GB` nominal for eMAGs) | *derived* from `H1` `0.309 GB/Gbp` (`2026-11-15_stage_H1/root_added.tsv`, `C-38` N1) and `CAMI` `0.403 GB/Gbp` (`C-38` N5); eMAGs have `10–40 kb` N50 |
| N2 | Ingest CPU rate | `4.11` CPU-min / Gbp (+ gzip term) | *measured* fit (`2026-11-16_ingest_estimate.md`, `C-38` N11, `DEC-56`) |
| N3 | Pfam 38.2 scan rate | `7.18` CPU-min / `10^6` residues | *measured* on `H1` (`C-24`, `C-38` N14); consistent with `C5` (`8,899.7` CPU-min / `9.11 Gbp`, `DEC-62`) |
| N4 | `IMG-UMAG` eukaryote bins on disk | `691` bins, `1,227,548` contigs, `14.435 Gbp`, `0` shared contigs | *measured*, `wilkcomp02`, 2026-11-19 (`2026-11-19_jgi_mags_iuig_dedup.md` §2.1); median `17.60 Mbp`, N50 `10.8 kbp` |
| N5 | `Tagirdzhanova 2024` lichen eMAGs | `325` eukaryotic MAGs, `10.16 Gbp` | *measured*, `wilkcomp02`, 2026-11-24 (`2026-11-24_lr2_nx1_inputs.md`, `CATALOG_GOAL.md` §4.1) |
| N6 | `Tara` `SMAGs` gene density | `~394` genes / Mbp (`~3.94M` genes / `10 Gbp`) | *derived* from *published* `>10M` predicted protein-coding genes over `~25.4 Gbp` (Delmont et al. 2022) |
| N7 | Mean microbial eukaryotic protein length | `~420–480` aa (`450` aa nominal) | *estimate* (~`1.4×` prokaryotic `323.6` aa *measured* in `P2`, `C-38` N16; consistent with eukaryotic multi-domain proteins) |
| N8 | `wilkcomp02` build envelope | `40` day / `80` night threads, `750G` RAM scope, `280` CPU-min/process hard guard | *measured* / policy (`DEC-73`, `DEC-74`, `AGENTS.md` §7.2) |

---

## 2. Source survey: counts, sizes, licences, and overlap

All external URLs, counts, and licences below were primary-source-checked via public documentation on 2026-11-24. Prokaryotic MAGs bundled in mixed studies (`TOPAZ`'s `4,022`, `OMD`'s `~34,667`, `Tagirdzhanova`'s `675`) are explicitly excluded from `C-60` because they are already covered by `GlobDB r226`, `SAG1`, `SPIRE v1`, and `NX1`.

| Source | Primary URL / Archive | Genomes / Bins | Size (Gbp) | Contigs | Licence & `redistributable` | Overlap & Deduplication Notes |
|---|---|---|---|---|---|---|
| **`IMG-UMAG` eukaryote bins** (`E1`/`E2`) | Local `/stor/work/Wilke/if2453/all-the-data/jgi_20261118/jgi_img_umag_euk_bins_20261118.zip` (`38.57 GB`) | `691` bins (*measured*) | `14.44 Gbp` (*measured*) | `1,227,548` (*measured*) | JGI Data Policy (`redistributable: false`) | `0` shared contigs vs `S1..PL1, V1, PL2, G1` (*measured*, `C-48d`); check `contig_uid` + `skani` vs `SMAGs`/`TOPAZ` (some IMG aquatic metagenomes overlap Tara/public SRA) |
| **`Tara` Oceans `SMAGs`** (Delmont et al. 2022, *Cell Genomics*, `10.1016/j.xgen.2022.100123`) | `https://www.genoscope.cns.fr/tara/` (SMAGs portal + Zenodo/SEANOE supplements) | `713` (`683` eMAGs + `30` SAGs, *published*) | `~25.2–25.4 Gbp` (*published*; `10 Mbp–1.3 Gbp` range) | `~2–4M` (*estimate*) | Open-access Genoscope / CC BY (`redistributable: true`; verify archive header) | Canonical *Tara Oceans* eukaryotic co-assembly set; **baseline** against which `TOPAZ` and `OMD` eMAGs are deduplicated |
| **`TOPAZ`** (Alexander, Hu et al. 2023, *mBio*, `10.1128/mbio.01676-23`) | `https://osf.io/gm564/` (`10.17605/OSF.IO/GM564`) | `988` eMAGs (*published*; `4,022` prok excluded) | `~20–35 Gbp` (*estimate*, `~20–35 Mbp`/eMAG) | `~2–5M` (*estimate*) | Open OSF / CC BY 4.0 (`redistributable: true`; verify OSF file licence) | Built from *Tara Oceans* `0.8–2,000 µm` particle-associated metagenomes via `EukHeist`; substantial species overlap with `Tara` `SMAGs` (same underlying Tara runs, different binning) — `skani` `>= 95%` ANI dedup against `SMAGs` required |
| **`OMD` eukaryotic ocean MAGs** (Paoli et al. 2022, *Nature* 607:111–118) | ENA `PRJEB45951`, Zenodo `10.5281/zenodo.4474310` | `148` eMAGs (*published*; prok MAGs excluded) | `~2.5–4.5 Gbp` (*estimate*) | `~0.2–0.5M` (*estimate*) | CC BY 4.0 / INSDC (`redistributable: true`) | Subset of *Tara*/Malaspina ocean bins; `skani` `>= 95%` ANI + `contig_uid` check vs `SMAGs` and `TOPAZ` |
| **`Tagirdzhanova 2024` lichen eMAGs** (*PLoS Biology*, Nov 2024) | ENA `PRJEB77567` (325 rows in `/tmp/if2453/c66/dl/lichen_euk_deferred.json`) | `325` eMAGs (*measured*, `C-66`; `294` fungal + `32` algal = `326` *published*) | `10.16 Gbp` (*measured*, `C-66`) | `~0.5–1.0M` (*estimate*) | CC BY 4.0 / INSDC (`redistributable: true`) | `675` prokaryotic lichen MAGs already staged in `NX1`; `0%` expected overlap with marine eMAGs; fungal mycobionts (*Cladonia*, *Lecanorales*) + *Trebouxia* algae |
| **Gaïa et al. 2026 *Parabasalia*/*Trichomonas* hosts** (*Nature*, `10.1038/s41586-026-11283-z`) | NCBI Datasets (`GCA_026262505.1` / `GCF_026262505.1` + `8` GenBank parabasalian accessions) | `9` genomes (*published*) | `~0.8–1.2 Gbp` (*estimate*; *T. vaginalis* G3 `175–180 Mbp`, *published*) | `< 50k` (*estimate*; G3 has `6` chr scaffolds) | INSDC public (`redistributable: true`) | Natural gut protist hosts of `EGV1` (*Egovirales*) and endogenous polintons; **intron-poor** (*T. vaginalis* G3 has `~63` active introns, *published*), making them ideal Wave-1 validation genomes |
| **`REMAG` & 2025–2026 polar/marine eMAGs** (+ `Kavagutti 2023` freshwater eMAGs) | `https://github.com/danielzmbp/remag`, Zenodo `10.5281/zenodo.23324293`; Figshare `20431332` (`euk-MAGs-fna.tar.gz`) | Tool + benchmark eMAGs (`~50–200` bins, *estimate*) | `~1–5 Gbp` (*estimate*) | `~0.1–0.5M` (*estimate*) | MIT (`REMAG`) / CC BY 4.0 (`Kavagutti`) (`redistributable: true`) | **Important distinction:** `REMAG` (HyenaDNA + Barlow Twins binner) is primarily an open-source binning tool/benchmark rather than a global catalog; ingest its published benchmark eMAGs + `Kavagutti` `euk-MAGs-fna.tar.gz` (already deferred in `CATALOG_GOAL.md` §4.1) and any 2025–2026 polar/marine eMAG releases |
| ***Mirusviricota* (*Okeanovirales*) intron-bearing viruses** (Medvedeva et al. 2025/2026) | Figshare `10.6084/m9.figshare.28955240` (genomes + intron-informed gene models) | `1,257` genomes (*published*) | `~0.4–0.55 Gbp` (*estimate*, `CATALOG_GOAL.md` §4.3) | `~5–15k` (*estimate*) | CC BY 4.0 (`redistributable: true`) | Contigs enter in `MIR1` (`C-70`) as `kind="virus"`; *Okeanovirales*/*Styxvirales*/*Demutovirales* carry dense spliceosomal introns + trans-spliced terminases that fragment under `pyrodigal-gv`. Layer Medvedeva's intron-informed GFF models onto `MIR1`'s contigs via `C-60`'s `gff-import` producer (zero duplicate nucleotides under `DEC-26`) |
| **`EukProt v3`** (Richter et al. 2022, *Peer Community J.* 2:e56) | Figshare `10.6084/m9.figshare.12417881.v3` | `993` spp total: `375` genomes + `56` SAGs (`431` genomic) vs `498` bulk + `47` SC transcriptomes + `17` ESTs (`562` RNA, *published*) | `~1–2 GB` protein FASTAs; `126` assembled transcriptomes (*published*) | N/A (proteins) | CC BY 4.0 (`redistributable: true`) | **Critical finding:** `EukProt v3`'s flagship files are **amino-acid protein FASTAs**, which violate `DEC-60 (2)`'s nucleotide-contig rule. Use `EukProt v3` in two compliant ways: (a) route its `126` nucleotide `EukProt_assembled_transcriptomes.v03` into `C-68` (`ET1`); (b) pin `EukProt v3` proteins in `refdb.yaml` as the homology reference DB for `miniprot`/`MetaEuk`, while fetching underlying nucleotide assemblies of the `431` genomic species from INSDC/JGI |
| **NCBI RefSeq `protozoa` & `fungi`** | `https://ftp.ncbi.nlm.nih.gov/genomes/refseq/{protozoa,fungi}/assembly_summary.txt` + NCBI Datasets v2 | `99` `protozoa` spp (*published*); `~500–800` `fungi` RefSeq spp (*estimate*, 7 phyla / 33 classes *published*) | `protozoa` `~2–4 Gbp`, `fungi` `~20–30 Gbp` (*estimate*) | `~0.2–1.0M` (*estimate*) | INSDC public (`redistributable: true`) | Curated **1 reference genome per species** with RNA-seq/evidence-backed RefSeq `.gff.gz` spliced gene models; zero overlap with `IMG-UMAG`/`SMAGs` eMAGs; overlaps `MycoCosm`/`EukProt` genomic species (RefSeq takes precedence under `DEC-26` because it is `redistributable: true`) |
| **JGI `PhycoCosm` & `MycoCosm`** | `https://phycocosm.jgi.doe.gov/`, `https://mycocosm.jgi.doe.gov/`, `https://genome.jgi.doe.gov/` (PI JGI login, `DEC-50 (2)`) | `PhycoCosm` `~300` datasets; `MycoCosm` `~3,867` fungal genomes (*published*, late 2026) | `PhycoCosm` `~15–30 Gbp`, `MycoCosm` `~120–150 Gbp` (*estimate*) | `~2–8M` (*estimate*) | JGI Data Policy; embargo/use-restricted entries `redistributable: false` | Requires PI JGI portal session (`DEC-50 (2)`); ships curated JGI GFF3 gene models; deduplicate against RefSeq `fungi`/`protozoa` and `IUIG`'s `721` Eukaryota rows (`C-48d`) before staging |

**Total estimated unique microbial eukaryote corpus (post-dedup):**
- **Open-access Waves 1–3 (`EPAR1` + `EPRO1` + `EFUN1` + `ELIC1` + `EOMD1` + `SMAGs` + `TOPAZ` net-new + `REMAG`/`Kavagutti`)**: **`~70–105 Gbp`** (*estimate*).
- **With JGI `E1`/`E2` (`14.44 Gbp`, *measured*) + `PhycoCosm` + `MycoCosm` (net-new)**: **`~180–260 Gbp`** (*estimate*, matching `CATALOG_GOAL.md` §4.3's `50–300 Gbp` range).

---

## 3. Data model for spliced eukaryotic CDSs (`AGENTS.md` §3 & §4)

### 3.1 Genome kind and enum widening

`src/biobox/registry/schemas.py` currently restricts `genome_input.kind` (line 820) and `genome_dim.kind` (line 909) to `("prokaryote", "virus", "plasmid")`, and `src/biobox/l1/genecall.py:169` mirrors `KINDS = ("prokaryote", "virus", "plasmid")`.

- **Change**: add `"eukaryote"` to `genome_input.kind`, `genome_dim.kind`, and `genecall.KINDS`.
- **Compatibility**: strictly additive enum widening; allowed before `rS0.1` lock (`AGENTS.md` §3.8), leaves all existing Parquet files and `schema_lock.json` column lists unchanged (only `json_schema_lock.json` enum domain updates via `scripts/update_schema_lock.py`).
- **Already present**: `feature_protein.edit_type` (`schemas.py:35`) **already defines `"intron"`** in `_EDIT_TYPES = ("none", "frameshift", "selenocysteine", "pyrrolysine", "intron")`. No change needed there.

### 3.2 Span-level `feature` + optional `l1_exons` table (Option A, recommended) vs exon-level `feature` (Option B)

| Axis | **Option A (Recommended): 1 gene-span `feature` per spliced CDS + optional `l1_exons` table** | Option B (Rejected): 1 `feature` row per exon |
|---|---|---|
| `feature_uid` (`ids.py:34`) | `blake3(contig_uid, min(exon_start), max(exon_end), strand, KIND["cds"])` — unchanged formula | Different `feature_uid` per exon; no stable gene-level `feature_uid` |
| `feature_protein` (`PK: (feature_uid, gene_call_version)`) | Clean **1:1** mapping `(feature_uid -> protein_uid)` with `edit_type = "intron"` (`n_exons > 1`) or `"none"` (`n_exons == 1`) | Broken: 10 exons map to 1 `protein_uid` (10:1), duplicating protein occurrences 10× |
| `contig_tapes` (`l3/tapes.py`) | **1 slot = 1 gene**; `±10`-gene neighborhood windows (`nbr`), `orfs_to_contig_starts/ends`, and `token_postings` work **unchanged** | 10-exon gene consumes 10 slots; destroys `±10` neighborhood windows and inflates `token_postings` 5–10× |
| `igr` extraction (`l1/igr.py`) | `igrs_of_contig()` merges overlapping `feature` spans and emits gaps **between gene spans only** — **introns are NEVER extracted as IGRs**, with **zero code changes** to `igr.py` | Every spliceosomal intron becomes an `igr` row and enters `dark_orf` search packs! |
| Exact exon/intron recovery | Stored in new table `l1_exons` (partitioned by `gene_call_version, stage`), joined by `feature_uid` when nucleotide exon slicing or intron analysis is needed | Directly in `feature`, at the cost of breaking every downstream L1/L3 table |

**Coordinate limits check (`schemas.py`):**
- `feature.start`, `feature.end`, `igr.start`, `igr.end` are `int64` (`COORD`).
- `contig_tapes.start_bps` and `end_bps` are `list<uint32>` (max `4,294,967,295` bp = `4.29 Gbp` per single contig). No microbial eukaryote chromosome approaches `4.29 Gbp` (even the `1.3 Gbp` `SMAGs` assemblies are fragmented across scaffolds), so `uint32` does not overflow.
- `contig_tapes.orfs_to_contig_starts/ends` are `list<uint16>` saturating at `65,535` with `slot_flags` bit 0 set (`tapes.py:106`); safe for any eukaryotic chromosome.

### 3.3 Proposed `l1_exons` table specification (`AGENTS.md` §3.3, §3.7)

Following `l1_crispr`'s proven pattern (`schemas.py:380-410`), where `ordinal` sits beside `feature_uid` for display/assembly order and is **never** part of the primary key (`AGENTS.md` §3.3):

- **Table name**: `l1_exons`, `layer = "L1"`, `location = "l1/exons"`
- **Partition keys**: `("gene_call_version", "stage")` (written only by spliced/eukaryotic callers; prokaryotic stages write no `l1_exons` files and require zero backfills, `AGENTS.md` §3.10)
- **Primary key**: `("feature_uid", "region_uid", "gene_call_version")` — where `region_uid = ids.region_uid(contig_uid, start, end, strand)` keys each exon by coordinates, never by ordinal
- **Sort key**: `("contig_uid", "feature_uid", "exon_ordinal")` — enables O(1) row-group pruning by `contig_uid` (`AGENTS.md` §4.4)
- **Columns**:
  - `feature_uid` (`UID`, non-null): parent gene-span `feature_uid`
  - `contig_uid` (`UID`, non-null)
  - `region_uid` (`UID`, non-null): `ids.region_uid(contig_uid, start, end, strand)`
  - `exon_ordinal` (`uint32`, non-null): 0-based exon index in **transcription order** (`5' -> 3'` on gene's `strand`), valid only beside `feature_uid`
  - `start` (`COORD`, non-null), `end` (`COORD`, non-null): 0-based half-open contig coordinates
  - `strand` (`STRAND`, non-null): `+1` or `-1`
  - `phase` (`uint8`, non-null): GFF3 reading-frame phase (`0`, `1`, or `2` bases to skip at exon 5' end)
  - `is_partial_5p` (`bool`), `is_partial_3p` (`bool`): whether the exon is truncated at a contig boundary or missing start/stop codon

### 3.4 Canonical isoform rule, nested genes, and trans-splicing

1. **Canonical single isoform per locus**: Because `feature_protein` has `PK = (feature_uid, gene_call_version)`, two alternative-splicing isoforms sharing the exact same `(start, end, strand)` would collide. Both `gff-import` and *ab initio*/homology callers MUST select **one canonical isoform per gene locus** deterministically: longest translated protein (`len(aa) DESC`), then fewest internal stops/ambiguous residues (`X_count ASC`), then fewest exons (`n_exons ASC`), then lexicographically smallest `protein_uid` (`protein_uid ASC`). Dropped isoforms are counted in the `genecall` ledger (`reason_code = "alt_isoform_dropped"`).
2. **Nested / intronic genes**: Eukaryotes occasionally encode a second gene inside a long intron (sense or antisense). Because `feature` rows are keyed by `(contig_uid, start, end, strand, kind)` and `contig_tapes` orders slots by `(start, end, strand, kind)` (`tapes.py:220`), nested genes produce distinct `feature_uid`s and receive `junction_topology = 4` (nested antisense) or negative `delta_bp` in `pair_geometry`, which `tapes.py` and `igr.py` already handle cleanly.
3. **Trans-spliced genes** (e.g., *Mirus_G_0318* terminase in *Okeanovirales*, where exons lie on opposite strands or separate loci): a single `[start, end, strand)` span cannot represent opposite-strand trans-splicing. In v1, `gff-import` rejects cross-strand/cross-contig trans-spliced gene records into the ledger (`reason_code = "trans_spliced_unsupported"`) or records their largest contiguous cis-span if annotated as separate CDS parts; this affects `< 0.01%` of genes (*estimate*).
4. **Spliced translation verification**: Every eukaryotic CDS stored in `protein_seq` is the **spliced** amino-acid sequence (exons concatenated in `exon_ordinal` order, reverse-complemented if `strand == -1`, translated with the genome's `genetic_code` — e.g. Table 1 default, or ciliate Tables `6`/`10`/`29`/`30`). Thus `protein_uid = blake3(spliced_aa)`, `l1_domains` (Pfam 38.2), `l2/cluster` (MMseqs2 NR90/NR50), `token_postings`, and `search_hit` work **without a single line of downstream change**.

---

## 4. Eukaryotic gene-calling architecture (`bio-box.sif`, CPU-only `wilkcomp02`)

`wilkcomp02` is a CPU-only host (`112` logical Xeon E7-4850 v3 threads; `40` day / `80` night thread cap under `DEC-74`; `750G` RAM scope under `DEC-73`; `280` CPU-min per-process cap under `AGENTS.md` §7.2; no approved GPU under `DEC-16`). `bio-box.sif` (`container/environment.yml`) already ships `mmseqs2=18.8cc5c`, `hmmer=3.4`, `infernal=1.1.5`, `samtools=1.25`, `skani=0.4.0`, and `minced=0.4.2`.

### 4.1 Caller comparison

| Approach / Tool | Paradigm | Bioconda / CPU-only? | Speed (est. per `10 Gbp`) | Peak RAM | Strengths for `C-60` | Weaknesses / Risks |
|---|---|---|---|---|---|---|
| **`--depth stored`** (already in `C-38a`) | Ingest sequences only, no calling | Yes (built-in) | `~0.7` CPU-h (*derived*, N2) | `< 4 GB` | Zero image change; locks `E1`/`E2` (`14.44 Gbp`) into `l0seq` immediately | No proteins/tapes until called later |
| **`gff-import`** (new `biobox.l1.gff_import` producer) | Parse curated GFF3/GBFF exons, re-translate from `SeqStore` | Yes (pure Python) | `~1–2` CPU-h (*estimate*) | `< 4 GB` | **Gold-standard accuracy** for RefSeq `protozoa`/`fungi`, `Parabasalia`, `PhycoCosm`/`MycoCosm`, and Medvedeva *Mirusviricota*; 100% deterministic; verifies every exon against `SeqStore` | Only applies to sources that ship curated spliced GFF3s (not unannotated eMAGs) |
| **`miniprot`** (Heng Li 2023) | Ultrafast C/SIMD protein-to-genome spliced aligner | Yes (`bioconda::miniprot`, tiny C binary) | `~30–120` CPU-h (*estimate*, DB-size dependent) | `8–32 GB` per batch | **10–30× faster than `MetaEuk`**, exact canonical splice-site scoring, tiny image footprint, handles fragmented eMAGs cleanly, deterministic when sorted | Homology-only: needs a pinned protein DB (`refdb.yaml`, e.g., `EukProt v3` + `ET1` + `UniRef50` eukaryote subset); misses lineage-specific orphan genes with no homolog |
| **`MetaEuk`** (Levy Karin et al. 2020, *Microbiome*) | MMseqs2 6-frame search + DP exon chaining | Yes (`bioconda::metaeuk`) | `~300–600` CPU-h (*derived* from *published* `~0.5 Mbp/min` multi-core rate) | `32–128 GB` (scales with ref DB) | Purpose-built for fragmented eukaryotic MAGs; shares MMseqs2 ecosystem; no species training | Slower and less splice-boundary accurate than `miniprot` (*published* BUSCO/Compleasm benchmarks); requires large reference DB; misses orphans |
| **`Augustus`** / **`GALBA`** (`miniprot` + `Augustus`) | gHMM *ab initio* (or `miniprot`-hinted) | Yes (`bioconda::augustus`) | `~1,500–5,000+` CPU-h (*estimate*) | `4–16 GB` / worker | Recovers *de novo* orphan genes when trained | Raw `Augustus` fails on diverse uncultivated eMAGs lacking species parameters; `GALBA` is too slow for `100+ Gbp` of fragmented eMAGs |
| **`GlimmerHMM`** | Older GHMM *ab initio* | Yes (`bioconda::glimmerhmm`) | `~200–500` CPU-h (*estimate*) | `< 4 GB` | Fast, small binary | Requires organism-specific trained matrices (only ~7 model species ship); poor accuracy on uncultivated marine protists |
| **`Helixer`** | Deep learning (BiLSTM/CNN + HMM) | **GPU-bound** (PyTorch/TF); CPU inference `20–50×` slower | `~2,000–8,000` CPU-h on CPU (*estimate*) | `32–128 GB` | Strong *ab initio* on fungi/plants/animals without homology DB | Pulls multi-GB ML stack into `bio-box.sif`; **no protist/SAR/excavate model**; infeasible on CPU-only `wilkcomp02` without `wilkcomp03` GPU approval (`DEC-16`) |
| **`C-68` `ET1` (`MMETSP` + `1KP`)** *(complement)* | Spliced cDNA transcriptomes -> `pyrodigal`/ORF caller (`table 1, 6, 10, 29, 30`) | Yes (built-in) | `~10` CPU-h (*estimate*) | `< 8 GB` | **Immediate** eukaryotic protein & domain diversity without spliceosomal introns; feeds `miniprot` ref DB | No genomic synteny (`1` transcript = `1` contig); complements rather than replaces `C-60` |

### 4.2 Recommended three-tier eukaryotic calling architecture

1. **Tier 1 (`--depth stored` immediately)**: Add `kind="eukaryote"` so `E1`/`E2` (`14.44 Gbp`) can be ingested at `--depth stored` right away without waiting for any caller binary.
2. **Tier 2 (`gff-import` producer for curated genomes)**: Build a pure-Python deterministic producer (`gene_call_version = "gff-import-v1"`) that reads a genome's GFF3/GBFF beside its FASTA, selects one canonical isoform per locus (§3.4), validates every exon's coordinates against `SeqStore`, re-translates the spliced CDS using `genome.genetic_code` (cross-checking against the provider's `.faa` where present and logging mismatches in the ledger), and writes `feature`, `feature_protein`, `protein_occurrence`, `protein_seq`, and `l1_exons`. This immediately unlocks **Wave 1–2 (`EPAR1`, `EMIR1`, `EPRO1`, `EFUN1`)** and **Wave 4 (`PhycoCosm`, `MycoCosm`)** with curated accuracy and zero new container binaries.
3. **Tier 3 (`miniprot` CPU-only spliced caller for unannotated eMAGs, verified by a small spike against `MetaEuk`)**: Add `miniprot` (pinned in `container/environment.yml`) and a pinned eukaryotic reference protein database in `refdb.yaml` (built from `EukProt v3` + `ET1` `MMETSP`/`1KP` + RefSeq `protozoa`/`fungi` proteins, `~15–25M` NR90 eukaryotic proteins, *estimate*). Run a bounded spike (`C-60c`) on 10 RefSeq protist/fungal genomes (using RefSeq GFF3 as independent oracle, `AGENTS.md` §5.2) + 20 `IMG-UMAG` `E1` bins comparing `miniprot` vs `MetaEuk` on recall, exact splice-boundary precision, CPU-min/Gbp, and peak RSS.
   - **Batching under `AGENTS.md` §7.2**: Follow `genecall.py`'s existing `plan -> batch -> merge` pattern (`B-48`/`B-50`), cutting `genome_input` in `genome_uid` order into whole-genome batches sized to `<= 220` CPU-min (leaving margin under the `280` CPU-min hard guard). For mega-genomes (`> 200 Mbp` like large `SMAGs`), split by contig chunks and merge deterministically by `(contig_uid, start, end, strand)`.
   - **Why `miniprot` over `Helixer` for v1**: Zero GPU requirement, tiny C binary in `bio-box.sif`, covers deep-branching protists (`SAR`, `Excavata`, `Haptista`, `Cryptista`) that `Helixer` has no trained model for, and stays comfortably under `wilkcomp02`'s `750G` RAM scope (`DEC-73`) and `280` CPU-min/process cap.

---

## 5. Storage and Pfam CPU budget per `10 Gbp` of microbial eukaryotes

All figures below are *derived* from `bio-box`'s measured rates (§1) or explicitly labelled *estimate*.

| Metric (per `10 Gbp` stored) | Prokaryote isolate (`P1..P6` / `C5`) | Unbinned metagenome (`SPX1`, `>= 5 kb`) | **Microbial eukaryote (`C-60`)** | Basis & Label |
|---|---:|---:|---:|---|
| Genes / `10 Gbp` | `~9.3M` (*derived*, `DEC-29`) | `~10.0M` (*measured*, `SPX1`: `35.76M / 35.69 Gbp`) | **`~3.5–4.0M`** (`3.94M` nominal) | *derived* from `Tara` `SMAGs` (`>10M / 25.4 Gbp`, N6) |
| Mean protein length | `323.6` aa (*measured*, `P2`) | `~260–300` aa (*estimate*) | **`~450` aa** (`420–480` aa) | *estimate* (N7) |
| Total residues / `10 Gbp` | `~3.01B` aa (*derived*) | `~2.8B` aa (*estimate*) | **`~1.58–1.80B` aa** (`1.77B` nominal) | *derived*: `3.94M × 450` aa |
| Storage: `depth="stored"` | `~3.1 GB` (*derived*, N1) | `~3.5–4.0 GB` (*derived*, N5) | **`~3.3 GB`** (`3.1–3.6 GB`) | *estimate* from N1/N5 (`10–40 kb` N50) |
| Storage: `depth="called"` (incl. `l1_exons`) | `~7.1 GB` (*derived*, `C-38` N3) | `~8.5 GB` (*measured*, `DEC-67`) | **`~5.8 GB`** (`5.5–6.5 GB`) | *derived*/*estimate*: `3.3 GB` stored + `~1.7 GB` `protein_seq`/L1/tapes/IGRs + `~0.8 GB` `l1_exons` (`~4.5` exons/gene × `~45 B`) |
| Storage: `depth="annotated"` (full Pfam 38.2) | `~7.4 GB` (*derived*) | N/A (`called` only) | **`~6.1 GB`** (`5.8–6.8 GB`) | *estimate*: `called` + `~0.25–0.30 GB` `l1_domains` |
| Pfam 38.2 CPU-hours / `10 Gbp` | `~163` CPU-h (*derived*, `C5` `DEC-62`) | `~330` CPU-h (*estimate*) | **`~190–215` CPU-h** (`212` CPU-h nominal) | *derived*: `1,773M` aa × `7.18` CPU-min/M aa / `60` (N3) |
| Pfam 38.2 wall-hours / `10 Gbp` on `wilkcomp02` | `~2.0 h` night / `~4.1 h` day | — | **`~2.7 h` night (`80`t) / `~5.3 h` day (`40`t)** | *derived* from `212` CPU-h and `DEC-74` thread caps |
| Release-pack / clustering delta per `10 Gbp` | `9.3M` proteins | `9.4M` distinct proteins | **`~3.5–3.9M` distinct proteins** | *derived* (`~40%` of prokaryotic clustering/token load per Gbp!) |

**Budget takeaway:** A `10 Gbp` microbial eukaryote stage is **cheaper** to cluster (`~40%` as many proteins) and **smaller on disk** (`~5.8–6.1 GB` vs `7.1–8.5 GB`) than a `10 Gbp` prokaryotic or metagenomic stage, and its full Pfam 38.2 scan finishes in a single `~2.7`-hour overnight run on `wilkcomp02`! Therefore:
- Compact high-value host/reference stages (`<= 15 Gbp` each: `EPAR1`, `EPRO1`, `ELIC1`, `EOMD1`) should build at **`--depth annotated`** (full Pfam 38.2; CRISPR off for eukaryotes since CRISPR-Cas is prokaryotic, saving `~30` CPU-min / `10 Gbp`).
- Large multi-stage compendia (`SMAGs` `~25.4 Gbp`, `TOPAZ`, `E1`/`E2` `14.44 Gbp`, `EFUN1..`, `MycoCosm` `~120+ Gbp`) should build at **`--depth called`** with hit-neighborhood Pfam supplied on demand via `bio-box annotate round` (`DEC-54`, `DEC-75`).

---

## 6. Recommended staging order

| Wave | Stage(s) | Source(s) | Est. Size | Depth & Caller | Why this order |
|---|---|---|---:|---|---|
| **Wave 0** *(unblocked now / minimal schema)* | `ET1` (`C-68`), `E1`, `E2` | `MMETSP` + `1KP` + `EukProt v3` `126` transcriptomes (`C-68`); `IMG-UMAG` `691` euk bins (`E1`/`E2` on `/stor/work`) | `~14.44 Gbp` (`E1`/`E2`, *measured*) + `ET1` | `ET1`: `called`/`annotated` (cDNA); `E1`/`E2`: **`stored`** | Zero caller wait; `ET1` immediately enables eukaryotic `ssRNAP` & immune-homolog search while `E1`/`E2` lock into `l0seq` |
| **Wave 1** *(intron-poor & curated hosts)* | `EPAR1`, `EMIR1`, `EPRO1` | Gaïa 2026 `9` *Parabasalia* hosts (`~1 Gbp`), Medvedeva `1,257` *Mirusviricota* intron models (`0` new bp over `MIR1`), RefSeq `protozoa` `99` spp (`~3 Gbp`) | `~4 Gbp` (*estimate*) | **`annotated`** via **`gff-import`** | Highest biological yield per CPU-h (`EGV1`/`MIR1`/`GV1` natural hosts + core protist immune homologs); exercises `l1_exons` + `gff-import` on curated, intron-poor & intron-rich truth sets |
| **Wave 2** *(curated RefSeq fungi & non-marine eMAGs)* | `EFUN1..3`, `ELIC1`, `EOMD1` | RefSeq `fungi` (`~20–30 Gbp`, `gff-import`), `Tagirdzhanova 2024` `325` lichen eMAGs (`10.16 Gbp`, *measured*), `OMD` `148` ocean eMAGs (`~3.5 Gbp`) | `~34–44 Gbp` (*estimate*) | `EFUN`: `called` (`gff-import`); `ELIC1`/`EOMD1`: `annotated` (`miniprot`) | First real-world test of `miniprot` eMAG calling on `ELIC1`/`EOMD1` (`redistributable: true`); fungal immune effectors (NLRs, gasdermins, cGAS/STING, viperins) |
| **Wave 3** *(large ocean eMAGs & `E1`/`E2` upgrade)* | `ESMAG1..3`, `ETOP1..2`, `E1`/`E2` upgrade | `Tara` `SMAGs` (`713` genomes, `~25.4 Gbp`), `TOPAZ` (`988` eMAGs, `skani`-deduped vs `SMAGs`/`OMD`), `IMG-UMAG` `E1`/`E2` (`14.44 Gbp`) | `~55–70 Gbp` (*estimate*) | **`called`** (`miniprot`) | Major marine protist / giant-virus host diversity; `SMAGs` enters first so it owns shared *Tara* species before `TOPAZ` and `E1`/`E2` (`redistributable: false`) |
| **Wave 4** *(JGI login batch & 2025–2026 eMAGs)* | `EPHY1..2`, `EMYC1..`, `EREM1` | JGI `PhycoCosm` (`~300` spp) & `MycoCosm` (`~3,867` spp, PI login, `DEC-50 (2)`), `REMAG` benchmark + `Kavagutti` `euk-MAGs` + polar eMAGs | `~100–150 Gbp` net new (*estimate*) | **`called`** (`gff-import` for JGI, `miniprot` for eMAGs) | Deferred behind public waves because JGI portals require PI login (`DEC-50 (2)`) and carry `redistributable: false` |

---

## 7. Draft `sources.yaml` entries (for PI approval; not yet added to `sources.yaml`)

```yaml
  # --- DRAFT C-60 entries (do not fetch until PI approves C-60 decisions) ---
  - name: ncbi_datasets_parabasalia_hosts
    what: 9 Parabasalia/Trichomonas anaerobic gut protist host genomes cited in Gaia et al. 2026 (Egovirales hosts, incl. GCA_026262505.1); stage EPAR1
    url: https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accessions}/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,PROT_FASTA
    record: https://doi.org/10.1038/s41586-026-11283-z
    license: INSDC (public)
    redistributable: true
    note: draft C-60; verify all 9 GenBank accessions from Gaia et al. 2026 supplementary table before fetch

  - name: ncbi_refseq_protozoa_fungi
    what: NCBI RefSeq reference genomes (FASTA + GFF3 + protein FAA) for unicellular protozoa (~99 spp) and fungi; stages EPRO1, EFUN1..
    url: https://ftp.ncbi.nlm.nih.gov/genomes/refseq/{division}/assembly_summary.txt
    url_genome: https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accessions}/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,PROT_FASTA
    license: INSDC / NCBI RefSeq (public)
    redistributable: true
    note: draft C-60; division in {protozoa, fungi}; filter to reference/representative genomes (1 per species)

  - name: tagirdzhanova_2024_lichen_emags
    what: 325 eukaryotic fungal/algal lichen MAGs from Tagirdzhanova et al. 2024 PLoS Biology (PRJEB77567, 10.16 Gbp); stage ELIC1
    url: https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accessions}/download?include_annotation_type=GENOME_FASTA
    record: https://www.ebi.ac.uk/ena/browser/view/PRJEB77567
    license: CC BY 4.0 / INSDC
    redistributable: true
    note: draft C-60; 325 accessions already enumerated in /tmp/if2453/c66/dl/lichen_euk_deferred.json

  - name: omd_paoli2022_euk_mags
    what: 148 eukaryotic ocean MAGs from Ocean Microbiomics Database (Paoli et al. 2022 Nature, PRJEB45951); stage EOMD1
    url: https://zenodo.org/records/4474310/files/{filename}
    record: https://doi.org/10.5281/zenodo.4474310
    license: CC BY 4.0
    redistributable: true
    note: draft C-60; prokaryotic OMD MAGs excluded (already in GlobDB r226 / SAG1 / SPIRE v1)

  - name: tara_oceans_smags_delmont2022
    what: Tara Oceans eukaryotic SMAGs (683 eMAGs + 30 SAGs = 713 genomes, ~25.4 Gbp; Delmont et al. 2022 Cell Genomics); stages ESMAG1..3
    url: https://www.genoscope.cns.fr/tara/localdata/data/SMAGs/{filename}
    record: https://doi.org/10.1016/j.xgen.2022.100123
    license: open-access (verify exact Genoscope/Zenodo archive licence header on fetch)
    redistributable: true
    note: draft C-60; verify exact Genoscope/SEANOE/Zenodo tarball filenames before fetch

  - name: topaz_alexander2023_emags
    what: TOPAZ 988 particle-associated Tara Oceans eukaryotic MAGs (Alexander, Hu et al. 2023 mBio); stages ETOP1..2
    url: https://osf.io/gm564/
    record: https://doi.org/10.17605/OSF.IO/GM564
    license: CC BY 4.0 (verify OSF project licence on fetch)
    redistributable: true
    note: draft C-60; exclude 4,022 prokaryotic MAGs; skani 95% ANI dedup against SMAGs and OMD

  - name: mirusviricota_medvedeva2025_introns
    what: 1,257 Mirusviricota genomes and intron-informed gene models (Medvedeva et al. 2025/2026); stage EMIR1 producer over MIR1 contigs
    url: https://figshare.com/articles/dataset/28955240
    record: https://doi.org/10.6084/m9.figshare.28955240
    license: CC BY 4.0
    redistributable: true
    note: draft C-60; reuse MIR1 nucleotides (DEC-26), import intron-informed GFF models via gff-import

  - name: eukprot_v3_refdb
    what: EukProt v3 predicted proteomes (993 eukaryotic species) + 126 assembled transcriptomes (Richter et al. 2022)
    url: https://figshare.com/articles/dataset/12417881
    record: https://doi.org/10.6084/m9.figshare.12417881.v3
    license: CC BY 4.0
    redistributable: true
    note: draft C-60; proteins go to refdb.yaml (miniprot/MetaEuk reference only, DEC-60 (2)); 126 nucleotide transcriptomes route to C-68 (ET1)

  - name: jgi_phycocosm_mycocosm
    what: JGI PhycoCosm (~300 algal/protist genomes) and MycoCosm (~3,867 fungal genomes) with JGI GFF3 gene models; Wave 4
    url: https://genome.jgi.doe.gov/portal/
    license: JGI Data Utilization Policy (embargo/use-restricted entries non-redistributable)
    redistributable: false
    note: draft C-60; requires PI JGI login (DEC-50 (2)); deduplicate against RefSeq protozoa/fungi and IUIG
```

---

## 8. PI decision questions and recommendations

| Q# | Decision needed | Options | Recommendation | Impact |
|---|---|---|---|---|
| **Q1** | **Spliced-CDS data model & `"eukaryote"` kind** | **(a)** Add `kind="eukaryote"`, store 1 gene-span `feature` per canonical spliced CDS (`edit_type="intron"` when `n_exons > 1`), add optional `l1_exons` table (`PK: (feature_uid, region_uid, gene_call_version)`), leave `igr`/`tapes`/`l2`/`l3` unchanged; **(b)** 1 `feature` per exon | **Option (a)** | Preserves `±10`-gene tape windows, prevents introns from polluting `igr`/`dark_orf`, zero backfill on built stages (`AGENTS.md` §3.10) |
| **Q2** | **Eukaryotic calling stack on CPU-only `wilkcomp02`** | **(a)** 3-tier: `--depth stored` for `E1`/`E2` now + pure-Python `gff-import` (re-translated & verified against `SeqStore`) for curated genomes + CPU-only `miniprot` (benchmarked vs `MetaEuk` in a small spike) for unannotated eMAGs; **(b)** `GALBA` (`miniprot`+`Augustus`); **(c)** `Helixer` (requires `wilkcomp03` GPU approval) | **Option (a)** | Zero GPU dependency; gold-standard accuracy on RefSeq/Parabasalia/JGI/Mirusviricota; fits `<= 280` CPU-min/process and `750G` RAM |
| **Q3** | **Intron-bearing eukaryotic viruses (`MIR1` *Mirusviricota*)** | **(a)** Layer Medvedeva's intron-informed GFF models onto `MIR1`'s contigs as a second `gff-import` producer (zero duplicate nucleotides under `DEC-26`); **(b)** keep `pyrodigal-gv` calls only | **Option (a)** | Recovers intact spliceosomal-intron genes (terminase, MCP, RNAP) in *Okeanovirales* without re-storing sequences |
| **Q4** | **`EukProt v3` disposition** | **(a)** Keep `EukProt v3` protein-only FASTAs out of `l0seq` (`DEC-60 (2)`), pin them in `refdb.yaml` as the `miniprot`/`MetaEuk` homology DB, and route its `126` nucleotide transcriptomes to `C-68` (`ET1`); **(b)** wait to fetch all `431` underlying genomic assemblies before calling eMAGs | **Option (a)** | Strictly honors `DEC-60 (2)` while immediately boosting `miniprot` protist coverage across all `993` eukaryotic taxa |
| **Q5** | **Import depth & Pfam 38.2 policy for eukaryotes** | **(a)** `--depth annotated` (full Pfam 38.2, `--no-crispr`) for compact Wave 1–2 host/reference stages (`<= 15 Gbp`, `~190–215` CPU-h / `10 Gbp`), `--depth called` + on-demand `annotate round` for large Wave 3–4 compendia (`SMAGs`, `TOPAZ`, `E1`/`E2`, `MycoCosm`); **(b)** `--depth called` for all eukaryotes | **Option (a)** | Full Pfam on `EPAR1`/`EPRO1`/`ELIC1`/`EOMD1` costs `< 1` overnight run total while keeping `100+ Gbp` waves fast |
| **Q6** | **Staging order & JGI login batch** | **(a)** Wave 0 (`ET1` + `E1`/`E2` `stored`) -> Wave 1 (`EPAR1`, `EMIR1`, `EPRO1`) -> Wave 2 (`EFUN1..`, `ELIC1`, `EOMD1`) -> Wave 3 (`ESMAG1..3`, `ETOP1..2`, `E1`/`E2` `called`) -> Wave 4 (`PhycoCosm`/`MycoCosm` PI JGI batch + `REMAG`/polar eMAGs); **(b)** call `E1`/`E2` before RefSeq/`EPAR1` | **Option (a)** | Validates `l1_exons` and `gff-import` on curated, intron-poor hosts (`Parabasalia`, RefSeq `protozoa`) before scaling to fragmented eMAGs |

---

## 9. Follow-up implementation cards (blocked on PI decision row)

1. **`C-60a` (L, schema & `stored`-depth `E1`/`E2`)**: Add `"eukaryote"` to `genome_input.kind`, `genome_dim.kind`, `genecall.KINDS`, and add `l1_exons` `TableSpec` to `src/biobox/registry/schemas.py`; regenerate locks via `scripts/update_schema_lock.py`; stage and ingest `IMG-UMAG` `E1`/`E2` (`14.44 Gbp`) at `--depth stored`.
2. **`C-60b` (M, `gff-import` producer + Wave 1 `EPAR1`, `EMIR1`, `EPRO1`)**: Implement `src/biobox/l1/gff_import.py` (canonical single-isoform selection, `SeqStore` exon re-translation check, `l1_exons` writer, `edit_type="intron"`, ledger accounting) with R1/R2 tests; stage `EPAR1` (`9` *Parabasalia* hosts), `EMIR1` (*Mirusviricota* intron models), and `EPRO1` (RefSeq `protozoa`, `99` spp) at `--depth annotated`.
3. **`C-60c` (M, `miniprot` vs `MetaEuk` spike & container pin)**: Benchmark `miniprot` vs `MetaEuk` (with `EukProt v3` + `ET1` ref DB) on 10 RefSeq protist/fungal genomes (RefSeq GFF3 oracle) + 20 `E1` bins on `wilkcomp02`; record recall, splice precision, CPU-min/Gbp, and peak RSS; pin the winning caller in `container/environment.yml` and ref DB in `refdb.yaml`.
4. **`C-60d` (M, Wave 2 `EFUN1..`, `ELIC1`, `EOMD1`)**: Wire the `C-60c` caller into `build-stage` (`plan -> batch -> merge`, `<= 220` CPU-min/batch); stage RefSeq `fungi` (`gff-import`), `Tagirdzhanova 2024` `325` lichen eMAGs (`10.16 Gbp`), and `OMD` `148` ocean eMAGs.
5. **`C-60e` (M, Wave 3 `ESMAG1..3`, `ETOP1..2`, and `E1`/`E2` upgrade)**: Fetch and `skani`-deduplicate `Tara` `SMAGs` (`713` genomes) and `TOPAZ` (`988` eMAGs); build at `--depth called`; upgrade `E1`/`E2` from `stored` to `called`.
6. **`C-60f` (M, Wave 4 JGI `PhycoCosm`/`MycoCosm` + `REMAG`/polar eMAGs)**: After PI JGI portal download (`DEC-50 (2)`), deduplicate `PhycoCosm`/`MycoCosm` against RefSeq and `IUIG`, stage with `gff-import` (`redistributable: false`), and add `REMAG` benchmark / `Kavagutti` / 2025–2026 polar eMAGs.
