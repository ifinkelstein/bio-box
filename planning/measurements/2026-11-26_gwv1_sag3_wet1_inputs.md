# Groundwater viruses (`GWV1`), OMZ single cells (`SAG3`) and wetland MAGs (`WET1`): relay, contents, staging (card C-81, DEC-79)

- **Date and host:** 2026-11-26 (`WET1`: 2026-11-26 22:55 to 2026-11-27, after card R-18), `wilkcomp02`. Every
  number carries a label (AGENTS §6): *measured* (on `wilkcomp02`, those days, by the command or file named),
  *derived* (formula and inputs), *estimate*, or *published*.
- **Script:** `scripts/stress/build_stage_pod/c81_gwv1_sag3_wet1.py` (`gwv1-stage`, `gwv1-ani`, `sag3-census`,
  `sag3-stage`, `sag3-ani`, `wet1-stage`, `wet1-ani`); tests `tests/r1/test_c81_gwv1_sag3_wet1.py`; templates
  `scripts/stress/build_stage_pod/stages/templates/GWV1.conf`, `SAG3.conf` and `WET1.conf`.
- **Read only:** the relayed files (§1), card C-74's index `/tmp/if2453/c74/index/` and GTDB R232's metadata
  `/tmp/if2453/c74/GTDB_R232_META/`, the data root `/tmp/if2453/r10/root`, `/tmp/if2453/s2/inputs/` (the `r26.11.4`
  stages and `GV1`), and the staged metadata of the `r26.11.5` stages on the pod (§2.3). Nothing was written to the
  root or to `/tmp/if2453/s2/inputs/`; nothing under `/stor` was overwritten; the GlobDB r226 archive was not read.
- **Written:** `gwv1_gwvc.fna.gz`, `sag3_omz.fna.gz` and `wet1_gwm.fna.gz` (each over 64 MB) to
  `/stor/scratch/Wilke/if2453/bio-box-staging/{GWV1,SAG3,WET1}/`, with `GWV1`'s 78 MB and `WET1`'s 143 MB `records.tsv`
  (§4); the other metadata to `/tmp/if2453/c81/meta/{GWV1,SAG3,WET1}/`; reports,
  ANI tables and inspection files to `/tmp/if2453/c81/work/`.

## 0. Summary

| | `GWV1` | `SAG3` | `WET1` |
|---|---:|---:|---:|
| Set | Groundwater Virome Catalogue, 280,420 vOTUs of 5 kb and more | OMZ single-cell compendium, 3,570 assembled SAGs | GWM wetlands, 4,164 dereplicated MAGs |
| Genomes staged | 280,420 | 3,570 cells, 3,568 genomes (2 repeated contig sets, DEC-27) | 4,164 |
| Contigs (distinct) | 280,420 | 299,181 | 1,749,962 |
| bp staged | 3,740,796,061 | 2,405,536,267 distinct (2,408,106,251 in all records) | 12,615,693,814 distinct (12,615,718,580 in all records) |
| Held by the root or `r26.11.4` (exact `contig_uid`) | 70 contigs, 18,912,716 bp (all `H1`) | 1,280 contigs, 28,838,169 bp (`GART1`, `G1`, `GART2`, `PL3`, `MV1`, `V2`) | 64 contigs, 193,761 bp (`GART1` 58, `GART2` 6) |
| Held by the `r26.11.5` stages on the pod before it | 0 | 52,307 contigs, 852,786,641 bp (`SAG1` 50,613: 1,100 GORG-Tropics cells; `GRARE1` 690; `MV2` 2; `MVEUKV` 2) | 127 contigs, 356,472 bp (`GRARE1` 66, `GRARE2` 37, `GRARE3` 24) |
| New bp (after both) | 3,721,883,345 | 1,523,911,457 | 12,615,143,581 |
| r232 status | (viruses: not asked) | in GTDB R232 574 (of which in GlobDB r232 439), not in GTDB R232 2,855, not visible 141 | not in GTDB R232 4,164 (R232 holds no genome of `PRJNA1430164`) |
| ANI, recorded (nothing cut) | 320 of 11,216 sampled vOTUs (2.9 %) meet the vOTU threshold with `MV1`, `V2`, `V3` or `MV2` | 1,739 of 3,568 meet ANI 96 / AF 50 with `SAG1`, `SAG2` or r226 archaea of the same genus (1,104 are the same cell) | 79 of 4,164 (1.9 %) meet ANI 96 / AF 50 with `NX1`, `KAVHIFI`, `DPN1`, same-genus `GART1/2`, `GRARE1..3` or r226 archaea (`LR2`'s FASTA is damaged, §5.5) |
| Left out | nothing (no repeat, no refusal) | nothing but 2 repeated contig sets listed once (DEC-27) | nothing (1,540 MAGs of the table were dereplicated away by the providers) |
| `s2_check_inputs.py` | exit 0, `problems: []`, `gene_call_mode` `{"virus": "gv-meta"}` | exit 0, `problems: []`, codes `[11, 25]` | exit 0, `problems: []`, codes `[11, 25]` |
| Template | `DEPTH=called`, `CRISPR=yes`, `AFTER=INPH2604`, `EXPECT=held_by H1 70` | `DEPTH=called`, `CRISPR=yes`, `AFTER=GWV1`, six `EXPECT` lines (§3.3) | `DEPTH=called`, `CRISPR=yes`, `AFTER=SAG3`, `EXPECT` for `GART1` 58 and `GART2` 6 |

All *measured* (the scripts' reports `/tmp/if2453/c81/work/{gwv1,sag3,wet1}_report.json`). Against the C-72 note's
estimates (*estimate*: `GWV1` 3.6 to 4.0 Gbp, `SAG3` 2.3 to 2.5 Gbp, `WET1` 12 to 13 Gbp): `GWV1` 3.741 Gbp and
`WET1` 12.616 Gbp are inside; `SAG3`'s 2.408 Gbp
is inside, but its new bases, 1.524 Gbp, are 34 % under the lower bound (*derived*: 1 − 1.524 / 2.3), because 1,100
of its cells are `SAG1`'s (§3.3). The estimate was of the set's size, which is met; the shortfall is reported here
(AGENTS §6) and was not taken as a stop, since nothing in the card's plan depends on `SAG3`'s new bases.

## 1. The bytes on the pod

`scripts/orchestrator/relay_to_pod.sh` (card C-73) from the laptop, `--dry-run` first, each list cut from the
directory's `SHA256SUMS` on `tower` (taken to the laptop with `scp`), files of 64 MiB and more to `/stor`, the rest
to `/tmp` (*measured*, the script's output, 2026-11-26 CST; every sha256 equal to `tower`'s record, which `tower`
took after the provider's md5 matched, card C-72 part 1):

| File | Bytes | Where on the pod | Time | sha256 |
|---|---:|---|---|---|
| `GWVC.vOTUs.fna.gz` | 1,166,850,604 | `/stor/scratch/Wilke/if2453/bio-box-staging/GWV1/dl/` | 130 s (8.98 MB/s), 18:55 | `377d5f59b3d28a8a9f2218009042f673b54fdb86a5cf20dcc41c87484e26bef0` |
| `GWVC.vOTUs.information.txt` | 62,086,246 | `/tmp/if2453/c81/dl/GWV1/` | 18:55 | `82e33893eaaf0b4f03e686fd1fd848d764e7243da818ee5e5cf5000eac563264` |
| `GWVC.vOTUs.virus-host.linkage.txt` | 59,543,590 | `/tmp/if2453/c81/dl/GWV1/` | 3 s, 18:55 | `24126ef80ef19463f7232a30549e31902f3572c80a1574ab9228aafa9e381bb3` |
| `OMZ_Compendium_SAGs_Full.tar.gz` | 747,416,178 | `/stor/scratch/Wilke/if2453/bio-box-staging/SAG3/dl/` | 79 s (9.46 MB/s), 18:57 | `a2a89c5de07a9045e283fc214b1322384873d3e030c810f748eded9a67a4351b` |
| `Table_S1_Metadata-template-Bio-Med-SAGdescriptor-OMZ_May_23_2023.xlsx` | 4,433,905 | `/tmp/if2453/c81/dl/SAG3/` | 1 s, 18:57 | `9f6d87f2528a1863f23461d7453fbebc8b829d5648ae3ff17caa62d574a64e64` |

## 2. `GWV1`

### 2.1 What the files hold

Read before staging (*measured*, a read-only pass `/tmp/if2453/c81/bin/inspect.py`, `work/inspect/inspect.json`):

- `GWVC.vOTUs.fna.gz`: **280,420 records, 3,741,506,338 bp**, no repeated name; one record per vOTU representative;
  names are MEGAHIT contig names with the well (`k141_2088200_length_42584_cov_21.0005_B10`). By length: 172,810
  records of 5 to 10 kb, 45,709 of 10 to 15 kb, 208 of 200 kb and more. The publication's 280,420 vOTUs of 5 kb and
  more (*published*, C-72 note §2.4): equal.
- `GWVC.vOTUs.information.txt`: one row per record (280,420 rows; every record has exactly one row, checked by the
  stager): CheckV quality (Complete 5,366, High 6,092, Medium 15,669, Low 210,488, Not-determined 42,805),
  MIUViG quality, completeness and its method (`DTR (high-confidence)` 4,615, `DTR (medium-confidence)` 676, ITR 33,
  provirus 42, AAI-based 139,477, HMM-based 92,772, NA 42,805), geNomad taxonomy (*Caudoviricetes* 239,998 of
  280,420; *Nucleocytoviricota* about 6,500), lifestyle, well, zone, BioProject `PRJNA858913`, BioSample and run.
- `GWVC.vOTUs.virus-host.linkage.txt`: virus-host links (vOTU, MAG, method such as `[CRISPR]`, GTDB r207 lineage of
  the host). It describes hosts, not sequence: relayed and recorded, not staged.
- `CPR.DPANN.genomes.tar.gz` of the same record was not fetched by card C-72 part 1 (its `sources.yaml` note), so it is
  not on `tower` and not read here; it stays in `CATALOG_GOAL.md` §4.

### 2.2 Staging (`c81_gwv1_sag3_wet1.py gwv1-stage`, report `work/gwv1_report.json`)

Rule: every record, in the file's order, `kind="virus"`, `genetic_code` 11 (pyrodigal-gv chooses its own),
`assembly_quality` `vOTU`; **DEC-65 (1)**, as card C-79 applied it to MetaVR (`c62_metavr_inputs.apply_dec65`): the
5,291 records whose CheckV completeness method is a direct terminal repeat were classified by their ends.

| | Records | bp | Label |
|---|---:|---:|---|
| DTR records, `exact` 20 to 1,000 bp repeat: second copy trimmed, staged circular | 5,210 | 710,277 trimmed | *measured* |
| DTR records, `partial` (a near repeat too): linear, untrimmed | 79 | | *measured* |
| DTR records, `low_complexity`: linear, untrimmed | 2 | | *measured* |
| Every other record (ITR, provirus, AAI, HMM, NA): linear, untrimmed | 275,129 | | *derived*: 280,420 − 5,291 |
| Genomes listed (DEC-27 repeated contig sets) | 280,420 (0) | 3,740,796,061 | *measured* |

`s2_check_inputs.py GWV1` on `objects_staging.json`: exit 0, `problems: []`, 280,420 genomes and `genome_input` rows,
codes `[11]`, kinds `["virus"]`, `gene_call_mode` `{"virus": "gv-meta"}` (*measured*). The licence is the Zenodo
record's `cc-by-4.0` (read 2026-11-26 by card C-72 part 1), so the objects carry `redistributable: true` (DEC-80 (3)).

### 2.3 Overlap by exact `contig_uid`

Owners in build order (the method of cards C-77 and C-79): the root's 13 stored stages `S1, P1..P6, H1, PL1, C5, V1,
G1, SPX1` (`contig_index`, default variants), then the 19 `r26.11.4` stages from their `genomes.json` under
`/tmp/if2453/s2/inputs/` in `a26_r26_11_4_chain.R26_11_4_CHAIN` order (so `V2`, `MV1` and `V3` among them); then, as
later holders, the `r26.11.5` stages on the pod before `GWV1` in plan §5's order: `GRARE1..3`
(`/tmp/if2453/c75/meta/`), `LR2`, `NX1` (`/tmp/if2453/c66/`), `SAG1` (`/tmp/if2453/c58/SAG1/`), `GV1`, `EGV1`
(`/tmp/if2453/c69/EGV1/`), `SAG2` (`/tmp/if2453/c77/meta/SAG2/`), `MV2`, `MVEUKV` (`/tmp/if2453/c79/`), `KAVHIFI`,
`INPH2604` (`/tmp/if2453/c80/`). `ECE2`, `MIR1`, `ASG1`, `DPN1`, `SYM2`, `LRP1`, `MIR2`, `MGE4` were not staged on the
pod when this ran: their overlap is card A-29's.

| Holder | Contigs | bp |
|---|---:|---:|
| `H1` (the huge phages, root) | **70** | 18,912,716 |
| `MV1`, `V2`, `V3` (`r26.11.4`) | **0** | 0 |
| Every other root and `r26.11.4` stage | 0 | 0 |
| `MV2` and every other `r26.11.5` stage listed above | **0** | 0 |
| **New** | **280,350** | **3,721,883,345** |

All *measured*. The 70 contigs `H1` holds (`shared_records.tsv`) are the template's only `EXPECT` line; since `MV2`
holds none of `GWV1`, no `EXPECT=held_by MV2` line is written (the line would read 0, which `stage_est.py` does not
take as a holder). The catch of C-72 §5.5 ("`GWV1` against `MV1`, `V2`, `V3`") finds no byte-identical vOTU: the reads
of `PRJNA858913` were released with the paper (2024), after IMG/VR v4.1 and before MetaVR; the vOTU-level test is §2.4.

### 2.4 The vOTU threshold on a sample (card C-48c's method)

`c81_gwv1_sag3_wet1.py gwv1-ani` (*measured*, `/tmp/if2453/c81/work/gwv1_ani/summary.json` and one skani table per
chunk, `<STAGE>.<n>.tsv`): the sample is every 25th non-repeated `GWV1` genome in FASTA order, **11,216 vOTUs,
149,851,979 bp**; the references are the whole of `V2` (`imgvr_v4.1_hc_hq_votu_reps.fna` and `inovirus_roux2019.fna`),
`MV2` (`mv2.fna`), `V3` (`vire_v3_dna_votus.fna`) and `MV1` (`metavr_mv1.fna`), about 80 GB of FASTA. Card C-48c's
skani settings through `c39_ani.skani` (`skani dist -c 30 -m 200 --min-af 0 --qi --ri`, from the image, 8 threads),
with the roles turned round so that skani holds the small side's sketches: the sample is skani's reference set and
each held stage's records are the queries, read in chunks of about 4 Gbp (21 skani processes in all). The vOTU
threshold is ANI >= 95 % over >= 85 % of the shorter sequence (rule `shorter`; read as the larger of the two aligned
fractions, since both count the same aligned bases, so no length of the held record is needed); `both` and
`near_identical` as `c39_ani`.

| Held stage | `shorter` | `both` | `near_identical` | skani runs | CPU-seconds (*measured*) |
|---|---:|---:|---:|---:|---:|
| `MV1` | 270 | 9 | 1 | 10 | 3,418.4 (largest run 1,638.3) |
| `V3` | 118 | 4 | 1 | 7 | 557.9 |
| `V2` | 102 | 6 | 1 | 3 | 839.1 |
| `MV2` | 16 | 0 | 0 | 1 | 29.8 |
| **any of them** | **320 (2.85 %)** | **12** | **1** | 21 | 4,845.2 |

| `shorter`, any held stage, by CheckV quality | Queries | Meet |
|---|---:|---:|
| Complete | 210 | 8 |
| High-quality | 247 | 7 |
| Medium-quality | 636 | 18 |
| Low-quality | 8,366 | 276 |
| Not-determined | 1,757 | 11 |

By length: under 10 kb 205 of 6,917; 10 to 30 kb 87 of 3,272; 30 to 100 kb 27 of 990; 100 kb and more 1 of 37. All
*measured*; percentages *derived* (meet / queries). Reading: about 2.9 % ± 0.3 % of `GWV1`'s vOTUs (*derived*: two
binomial standard errors at n = 11,216) are, at the vOTU threshold, a vOTU or a fragment of a virus that MetaVR, IMG/VR
or VIRE already holds, almost all as the shorter of the pair (only 12 cover each other both ways), against 0 exact
`contig_uid` matches; the C-72 note's estimate was under 10 % (*estimate*). They are staged: DEC-26 deduplicates by
exact `contig_uid` only, and the card records the sample, it does not cut by it.

## 3. `SAG3`

### 3.1 What the archive and Table S1 hold

- `OMZ_Compendium_SAGs_Full.tar.gz` (*measured*, `sag3-census` and the read-only pass of §2.1,
  `work/inspect/sag3_tar_members.tsv`): **3,570 entries, all `Final_Catalog/<name>.fna`, 2,469,972,902 bytes**
  uncompressed; no directory entry, no other member. The publication's 3,570 sequenced, assembled and decontaminated
  SAGs (of 5,129 sorted, *published*, C-72 note §2.6): equal. Member names carry the cell's name and the steps that
  made the file: 3,497 plain `<name>.fna`, 43 `_output_clean`, 18 `_contigs`, 10 `_contigs_output_clean`, one
  `_contigs_output_clean_output_clean`, one `_contigs_output_clean_2000`, and `<cell>_AB-905_clean`/`_AB-908_clean`
  pairs. Contig names come in two forms: `<cell>_NODE_...` or `<cell>_<...>contig...`, and IMG's
  `<taxon oid>_Ga<...>_<n>_..._SCGC_<cell>...` with no cell name in front.
- Table S1 (12 sheets, read with `zipfile` and `xml.etree`): `Outputs-Genome_Sequencing` names 3,571 cells and gives
  BioSample (3,525 rows) and GenBank WGS (2,440 rows) accessions and IMG ids; `Outputs_GTDB-TK` holds GTDB-Tk v2.1.0's
  classification (R207) of 3,570 genomes; `Outputs-Genome_CheckM_Statistic` and `Outputs-Genome_Quality_Classifi` the
  CheckM 1.2.1 completeness and the quality class. **The GTDB-Tk sheet cuts six pairs of names**: `AB-754_C22`,
  `AB-754_D02`, `AB-754_D23`, `AB-754_E20`, `AB-754_G03` and `AB-754_G05` each appear twice, for the cells
  `<name>_AB-905` and `<name>_AB-908` that the Sequencing sheet and the archive name in full. The stager therefore takes
  names from both sheets and gives a GTDB-Tk classification only to a name that sheet lists once; the 12 cells of those
  pairs carry none (their code screen runs without a taxonomic hint).

### 3.2 Matching members to Table S1, and the accession test (`sag3-stage`, report `work/sag3_report.json`)

| | Cells | Label |
|---|---:|---|
| Member stem equal to a Table S1 name (`exact`) | 3,036 | *measured* |
| Longest Table S1 name followed in the stem by `_` (`prefix`: `_output_clean`, `_contigs`, `_AB-905_clean`, ...) | 505 | *measured* |
| A name two members would take (`tie`) | 0 | *measured* |
| No name fits (`none`; staged under the member's stem, not visible to the test) | 29 | *measured* |
| Table S1 names no member took | 101 of 3,642 | *measured* (the six short GTDB-Tk names, and Sequencing rows of cells not in the archive) |
| Cells with a BioSample in Table S1 / with a WGS accession | 3,429 / 2,384 | *measured* |
| Cells without a GTDB-Tk classification | 41 | *measured* |

The accession test cannot use C-74's `lookup` directly: Table S1 gives WGS project accessions (`AFHZ00000000`,
`CABWYC000000000`) and BioSamples, not `GCA_` assemblies. GTDB R232's own metadata (`ar53`, `bac120`, columns
`ncbi_wgs_master` and `ncbi_biosample`, read from `/tmp/if2453/c74/GTDB_R232_META/`) maps the WGS project's letters,
then the BioSample, to a GTDB R232 genome, and `lookup` places that accession. GlobDB r232 is then seen only through
GTDB R232 (its GTDB resource); for a cell GTDB R232 does not hold there is no `GCA_` to ask GlobDB with. **No cell is
removed for its status** (the standing rule of 2026-11-26); `r232_status.tsv` holds it for every cell, with GTDB-Tk's
R207 classification and fastANI reference from Table S1 beside it.

| r232 status | Cells | bp | Label |
|---|---:|---:|---|
| in GTDB R232 and GlobDB r232 (`in_gtdb_and_globdb_r232`) | 439 | 503,430,519 | *measured* |
| in GTDB R232 only (`in_gtdb_r232`) | 135 | 163,583,677 | *measured* |
| an INSDC accession that no GTDB R232 genome carries (`not_in_gtdb_r232`) | 2,855 | 1,638,852,540 | *measured* |
| no accession in Table S1, or no Table S1 row (`not_visible`) | 141 | 102,239,515 | *measured* |
| **All** | **3,570** | **2,408,106,251** | |

All 574 cells found in GTDB R232 were found by their WGS project (0 by BioSample alone); 439 of them are GTDB R232
species representatives (*measured*). 574 of 3,429 cells visible to the test (16.7 %, *derived*) are in GTDB R232,
under the C-72 note's estimate of 20 to 50 % represented at species level (*estimate*; that estimate was of species,
this count is of the cells themselves).

### 3.3 Staging and overlap

Rule (card C-77's): the cells in name order for the metadata, the archive's order for the FASTA; whole contigs; the
code screened per cell (`c58.screen_genetic_code`, 11, 25, 4, the GTDB-Tk classification as the hint); a contig name
that does not begin with the cell's name and `_` is given the cell's name and `|` (132,751 of 301,688 records, *measured*;
every name in the file unique, refused otherwise); `kind="prokaryote"`, `assembly_quality` `contig`; one BGZF FASTA.

| | Value (*measured*) |
|---|---:|
| Cells / FASTA records | 3,570 / 301,688 |
| Genomes listed (DEC-27: 2 repeated contig sets listed once) | 3,568 |
| Distinct contigs (bp) | 299,181 (2,405,536,267) |
| Contigs in more than one genome | 351 (`s2_check_inputs.py` accepts them, DEC-27) |
| Genetic codes | code 11: 3,569; code 25: 1; conflicts 0 |
| Held by the root or `r26.11.4` (exact `contig_uid`, build order) | 1,280 contigs, 28,838,169 bp: `GART1` 833 (15,795,258 bp), `G1` 321 (9,515,941), `GART2` 61 (1,779,565), `PL3` 38 (647,026), `MV1` 26 (1,058,931), `V2` 1 (41,448) |
| New after the root and `r26.11.4` | 297,901 contigs, 2,376,698,098 bp |
| Held by `r26.11.5` stages before it on the pod | `SAG1` 50,613 contigs (840,766,748 bp), `GRARE1` 690 (11,886,838), `MV2` 2 (100,544), `MVEUKV` 2 (32,511); `GWV1` and the rest 0 |
| **New after both** | **246,594 contigs, 1,523,911,457 bp** |
| `s2_check_inputs.py SAG3` (with `objects_staging.json`) | exit 0, `problems: []`, 3,568 genomes, 300,422 `genome_input` rows, codes `[11, 25]` |

**`SAG3` repeats 1,100 of `SAG1`'s GORG-Tropics cells.** The contigs `SAG1` holds belong to 1,100 `SAG3` cells, all
with GORG-Tropics names (`AG-3xx-...`, BioProject `PRJEB33281`), and 1,089 of them are held whole (*measured*, a
read-only count over `later_records.tsv` and `records.tsv`). The compendium gathered the GORG-Tropics cells sorted from
oxygen-deficient waters with its own. Nothing is removed for it: the cells are listed in `SAG3`'s `genomes.json`
(DEC-27) and `SAG1`, earlier in the chain, owns their contigs (DEC-26), so the build stores them once; card A-29 writes
the `EXPECT=held_by SAG1` line. `GRARE1` holds 7 `SAG3` cells whole (archaeal or unclassified-domain SAGs that GlobDB
r226 carries). So of `SAG3`'s 2.408 Gbp, 1.524 Gbp is new (*measured* above), 34 % under the lower bound of the C-72
note's 2.3 to 2.5 Gbp (*derived*: 1 − 1.524 / 2.3; an *estimate* of the set's size that did not know of the shared
cells). The set's own size, 2.408 Gbp, is
inside that estimate.

`redistributable: true`: both figshare articles state CC BY 4.0 (read 2026-11-26 by card C-72 part 1; DEC-80 (3)).

### 3.4 ANI against `SAG1`, `SAG2` and the r226 archaea of the same genera (recorded, not cut)

`c81_gwv1_sag3_wet1.py sag3-ani` (*measured*, `/tmp/if2453/c81/work/sag3_ani/summary.json`, `sag3_ani_hits.tsv`,
`sag3_vs_refs.tsv`): one FASTA per genome written for `SAG3` (3,568 genomes, a contig taken by its first genome in
`records.tsv`), `SAG1` (5,537) and `SAG2` (9,698), and for the GlobDB r226 archaea of card C-75's reference set
(`/stor/scratch/Wilke/if2453/bio-box-staging/GRARE/globdb_r226_archaea.fna.gz` with
`/tmp/if2453/c75/meta/archaea/archaea_list.tsv`) whose GTDB genus is one that GTDB-Tk gives a `SAG3` archaeal cell:
15 genera (`Nitrosopelagicus`, `Nitrosopumilus`, `Poseidonia`, `Thalassarchaeum`, `MGIIa-I`, `MGIIa-K1`, `MGIIa-L1`,
`MGIIa-L2`, `MGIIb-O1`, `MGIIb-O2`, `MGIIb-O3`, `MGIIb-O5`, `CG-Epi1`, `ARS21`, `CABZYC01`), 1,009 genomes. One
`skani dist -t 6 --min-af 50` (skani's default sketch), queries against all references; a genome meets when
ANI >= 96 over an aligned fraction >= 50 % (the larger of the two, as card C-66's dedup); the best reference is kept.

| `SAG3` genomes meeting a reference | Genomes | Label |
|---|---:|---|
| the same cell in `SAG1` (same name), sharing exact contigs (§3.3) | 1,099 | *measured* |
| the same cell name in `SAG1`, no exact contig shared (a re-assembly) | 5 | *measured* |
| another `SAG1` cell (same species, a different cell) | 384 | *measured* |
| a `SAG2` cell | 196 | *measured* |
| a GlobDB r226 archaeon of the same genus | 55 | *measured* |
| **any** | **1,739 of 3,568** (48.7 %, *derived*) | |
| of them at ANI >= 99.9 | 1,119 | *measured* |

skani used 351.9 CPU-seconds (*measured*, children's rusage). **Nothing is cut by it** (the standing rule: an ANI
match to a held stage is recorded; the card does not say `SAG3` is deduplicated against `SAG1` or `SAG2`, and a
different cell of the same species is a different genome). Reading: beyond the 1,104 GORG-Tropics cells it repeats,
635 `SAG3` cells (*derived*: 384 + 196 + 55) fall in a species `SAG1`, `SAG2` or GlobDB r226's archaea already has.
The card's bacterial same-genus references from GlobDB r226 would need the r226 archive stream, which was not run
for the reason of §5.3 (CATALOG_GOAL §4.1).

## 4. Staged files

| File | Where | Bytes | sha256 |
|---|---|---:|---|
| `gwv1_gwvc.fna.gz` | `/stor/scratch/Wilke/if2453/bio-box-staging/GWV1/` | 1,183,215,792 | `e387439dbb6b4838d271325837e62ec5d783591304e022f80df56052eb9e16eb` |
| `genomes.json` | `/tmp/if2453/c81/meta/GWV1/` | 38,314,260 | `b920b2c5ed68082ffca013328166a014ee16f4228181fdde16b1be9074c6b84d` |
| `objects.json` | same | 426 | `3dc1500e68c418e51a4495355603bdf8d1e31cb9dd0caa2898920fa04f2a3cd7` |
| `objects_staging.json` | same | 447 | `5d39975ea97f5ec5b2d86fcaedc0b6825403adca418b4181ff3e12eae4a0e868` |
| `records.tsv` | `/stor/scratch/Wilke/if2453/bio-box-staging/GWV1/` (over 64 MB) | 78,157,462 | `d7c65e3c4f09c6d71be86dd77ed42b732855425d16246217231b42221b53730b` |
| `repeated_records.tsv` | same | 39 | `d109495ee76e4643c181eae183c452c9a84560f97a79bbe22d6630780368c565` |
| `shared_records.tsv` | same | 6,526 | `a5edd262fb904e4594e7e439923faf357a514775e17446d60aa6d0ff83715a52` |
| `later_records.tsv` | same | 40 | `35715d5e3e9f7b84789769b5e7130ae87a7e93e2adfa6cae8ac0222b940589fb` |
| `trims.tsv` | same | 942,800 | `b65f35f08a9a1e9fa5263f3337b3de3a2929fc186acaff3df2a79951525a18b2` |
| `sag3_omz.fna.gz` | `/stor/scratch/Wilke/if2453/bio-box-staging/SAG3/` | 759,962,632 | `2c62483fe32de8915f6795df9cca28b4e7eff5d7e68413d4af474c1852491b8d` |
| `genomes.json` | `/tmp/if2453/c81/meta/SAG3/` | 11,196,968 | `9e7c895453ace069ebe0f9205e582f4a8edb93a21fbd12e8e8faf729251330d8` |
| `objects.json` | same | 419 | `a0534151a99c3191bba77ad0e9b6ac6471030454bf1e5e4f5f1310f3efa8ca40` |
| `objects_staging.json` | same | 440 | `7ae69a31782d1826ae44258d6c762063b9d533ef56b9408ba90ce6de9a2655ac` |
| `records.tsv` | same | 24,490,279 | `75483af3acf17f4c08a3778ad6f47e4328a8aa32942a4e88d1cfe69e5612310e` |
| `repeated_records.tsv` | same | 162 | `2c3b31d860d41f6f779b4d7afad6d6e541a181c917aa02e6ae36de0c1515aaa5` |
| `shared_records.tsv` | same | 82,757 | `32cae5c995f8d05408790b266bf0cf594cb08f632d205d6498def8e58c134919` |
| `later_records.tsv` | same | 3,297,388 | `fd811537bdc18247b6917b45b23f4c9f52c82811d97c1a4a94372e4f4d8514aa` |
| `genetic_codes.tsv` | same | 624,764 | `00a3e33892b79d7d9156f95aa76874b4f6458e9b28627a98c82b10a9cd04cb76` |
| `sag_match.tsv` | same | 373,216 | `4936a65edfa707c546fcc2eab4985c13c98dabb8182805571411d5b46bad54f3` |
| `r232_status.tsv` | same | 662,415 | `3fc100d20d16f1a989fd042453a090e1ce6fe96ffefe971085e66a9105ce02fc` |
| `wet1_gwm.fna.gz` | `/stor/scratch/Wilke/if2453/bio-box-staging/WET1/` | 3,923,325,226 | `e77e9f6afa156b2706167c0d3149178ba35c1fde7aa10b44571a65b4c4b7bba9` |
| `records.tsv` | same (over 64 MB) | 142,712,258 | `1de4ab81bfaee24366ce4ae0bc92c3e5ba369d93ae023e039c5adb0865e85c14` |
| `genomes.json` | `/tmp/if2453/c81/meta/WET1/` | 63,444,540 | `b931b90abd9c98065add177a9dde5257ab5875ebd53adad52bc5e8ed94056b2b` |
| `objects.json` | same | 462 | `2d4c6a24ab834954997f9d8fc88b91a9397a89e1fbc19281cc935f8e6cd55a00` |
| `objects_staging.json` | same | 483 | `764fbb0b8fcd73e124b5c1f57ef3fba05896d006143274d966afb880e1f58684` |
| `repeated_records.tsv` | same | 39 | `d109495ee76e4643c181eae183c452c9a84560f97a79bbe22d6630780368c565` |
| `shared_records.tsv` | same | 4,208 | `99c0413f8e9178806c14f768cd1c4458d348fee95578be360538dca756468a47` |
| `later_records.tsv` | same | 8,420 | `8121ed5613103a5389c0b384267fb53c4ae2b94cdf24d83b365ad1d3c422e2c5` |
| `genetic_codes.tsv` | same | 745,885 | `01dcb9d25a7a13c72a5fceb28672a0b0494ca6bc34a0a55777b8af84ed9870d7` |
| `r232_status.tsv` | same | 822,613 | `868d1657b217eac601e928c84263e70207b933d6bf1eed2f868ae1b3f380ce4e` |

All *measured* (`ids.sha256_file`, the scripts' reports). `objects.json` names `/tmp/if2453/s2/inputs/<STAGE>/<FASTA>`
(card A-29 copies there); `objects_staging.json` the file on `/stor`. The templates' `SUM` lines are the rows for
`genomes.json`, `objects.json` and the FASTA. Every metadata file is under 64 MB except `GWV1/records.tsv`
(78,157,462 B), which was therefore moved, after the ANI step that reads it, to
`/stor/scratch/Wilke/if2453/bio-box-staging/GWV1/records.tsv` (`cp -n`, `sha256sum` equal to the row above, then the
`/tmp` copy removed; *measured*, 20:29 CST); card A-29 takes it from there. `stage_est.py conf` on both templates (`python -I`,
`--tree`, in the image): exit 0 (*measured*).

**Mirror to `tower`** (`scripts/orchestrator/relay_to_tower.sh`, card C-49, from the laptop, one stream, *measured*,
the script's log, 2026-11-26): `GWV1`: the FASTA, 1,183,215,792 B in 96 s (12.33 MB/s), then the eight metadata files,
`DONE: 1 verified` and `DONE: 8 verified, 0 failed`, `SHA256SUMS` on `tower` 9 lines;
`SAG3`: the FASTA, 759,962,632 B in 42 s (18.09 MB/s), then the ten metadata files, `DONE: 1 verified` and
`DONE: 10 verified, 0 failed`, `SHA256SUMS` 11 lines. `WET1` (2026-11-27): the FASTA and `records.tsv` from `/stor`,
`DONE: 2 verified, 0 in place and skipped, 0 failed` (`records.tsv` 142,712,258 B in 6 s), then the eight metadata
files, `DONE: 8 verified, 0 in place and skipped, 0 failed`, `SHA256SUMS` 10 lines; the mirror ran while the ANI step
of §5.5 ran, which writes only to `work/` and changes no staged file. All to
`tower:/mnt/user/research-bulk/bio-box/sources/{GWV1,SAG3,WET1}/`; no name there has a space.

## 5. `WET1`

### 5.1 The relay, and the two names of the archive

`WET1`'s genome archive is, on `tower` and at the provider, **`Dereplicated genomes.zip`** (3,885,785,718 bytes,
sha256 `3b63926fb66fcff2522d47e34f2f12128a25bdbb39e97c6917a390d8d9d70f8b` in `tower`'s `SHA256SUMS`). The first
`relay_to_pod.sh --dry-run` (2026-11-26 18:54 CST, *measured*) refused the space, exit 2, before any remote call:
`relay_to_pod: list line 2: refused path (only letters, digits, . _ = + , - and /; no . or .. or leading -)`. Card
R-18 then gave the relay `--space-as-underscore`, the change this note had proposed. With it (*measured*, the
script's output, 2026-11-26): the dry run printed the map, `copy Dereplicated genomes.zip 3885785718 (pod name
Dereplicated_genomes.zip)`; the run copied 3,885,785,718 bytes in 237 s (16.40 MB/s), sha256 equal to the record, to
`/stor/scratch/Wilke/if2453/bio-box-staging/WET1/dl/`. **On the pod the archive is `Dereplicated_genomes.zip`, the
relay's mapping of the provider's name**; nothing on `tower` was renamed or linked, and `tower`'s mirror of the staged
files (§4) has no space in any name. `MAG_meta.xlsx` (806,764 bytes, sha256
`b3ef27facf8e36eb6ef4bcbaf1304d34e6a28c8545c6fa0feb4ed635ce6e00f0`) went to `/tmp/if2453/c81/dl/WET1/` (0 s).

### 5.2 What the archive and the table hold

- `Dereplicated_genomes.zip` (*measured*, a read-only `zipfile` listing in the image, and `wet1_members` in the stage
  report): **4,165 entries: one directory `dereplicated_genomes/` and 4,164 files `dereplicated_genomes/MAGnnnnn.fa`**,
  12,847,878,267 bytes uncompressed, nothing else (no `__MACOSX`). The publication's 4,164 dereplicated genomes
  (*published*, C-72 note §2.4): equal. Contig names are each sample's assembly names (`BYD1_k127_11492162`), without
  the MAG id.
- `MAG_meta.xlsx`: one sheet, a title row ("Table S2. The detailed information of all MAGs used in the study") and a
  header row, then **5,704 MAGs** (the publication's 5,704 before dereplication: equal), each with its id
  (`MAG00001`), bin name (`BYD1_bin.10`), sample (`SMP0001`), GTDB-Tk classification, CheckM completeness and
  contamination, quality score and MIMAG tier. 1,540 of them are not in the archive (dereplicated away, *derived*:
  5,704 − 4,164). **The table gives no accession**: no BioSample, no GenBank assembly.

### 5.3 The accession test, and why the GlobDB r226 stream was not run

With no accession in the table, the test asks GTDB R232's own metadata (`ar53`, `bac120`) for every genome of
BioProject `PRJNA1430164` (the card's) and matches its `ncbi_isolate` or `ncbi_strain_identifiers` to the MAG's bin
name; a match would go through C-74's `lookup` as for `SAG3`. **GTDB R232 holds no genome of `PRJNA1430164`**
(0 rows, *measured*), which fits a project published in June 2026, after R232's April freeze. So every MAG is
`not_in_gtdb_r232` (4,164, 12,615,718,580 bp); none is removed. GlobDB r232 cannot be asked without an INSDC accession:
the r232 status rests on C-74's index and R232's metadata, as for the other stages.

The card's one stream of the GlobDB r226 archive was there to fetch same-genus references for `skani` 96 / AF 50 as a
**dedup**, that is, a cut. The standing rules have turned every such cut into a record (no genome in hand is left out
for r226 or r232 membership; an ANI match is recorded, not cut), and the card names no other use of those references.
So, as the orchestrating session asked, **the r226 archive was not streamed** (about 3.6 h, *estimate*, for a table
nobody would act on); the comparison with GlobDB r226 rests on C-74's index and on the r226 genomes already on the pod:
the `GART1`, `GART2` and `GRARE1..3` genomes of the same genera, and C-75's r226 archaea (§5.5).

### 5.4 Staging and overlap (`wet1-stage`, report `work/wet1_report.json`)

Rule (card C-77's, as for `SAG3`): the MAGs in id order for the metadata, the archive's order for the FASTA; whole
contigs; the code screened per MAG (`c58.screen_genetic_code`, 11, 25, 4, the GTDB-Tk classification as the hint);
every contig name given the MAG id and `|` (all 1,749,972, *measured*: none begins with its MAG id; every name in the
file unique, refused otherwise); `kind="prokaryote"`, `assembly_quality` from the MIMAG tier (`mag_hq` or `mag_mq`;
another tier is refused); one BGZF FASTA.

| | Value (*measured*) |
|---|---:|
| MAGs / FASTA records | 4,164 / 1,749,972 |
| Genomes listed (DEC-27 repeated sets) | 4,164 (0) |
| Distinct contigs (bp) | 1,749,962 (12,615,693,814); 10 contigs in more than one MAG |
| Genetic codes | code 11: 4,160; code 25: 4; conflicts 0 |
| Held by the root or `r26.11.4` (exact `contig_uid`, build order) | 64 contigs, 193,761 bp: `GART1` 58 (169,130 bp), `GART2` 6 (24,631) |
| Held by `r26.11.5` stages before it on the pod | 127 contigs, 356,472 bp: `GRARE1` 66 (182,953), `GRARE2` 37 (91,903), `GRARE3` 24 (81,616) |
| **New after both** | **1,749,771 contigs, 12,615,143,581 bp** |
| `s2_check_inputs.py WET1` (with `objects_staging.json`) | exit 0, `problems: []`, 4,164 genomes, 1,749,972 `genome_input` rows, codes `[11, 25]` |

The later holders read, in plan §5's order: `GRARE1..3`, `LR2`, `SAG1`, `NX1`, `GV1`, `EGV1`, `ECE2`, `MIR1`
(`/tmp/if2453/c70/meta/`), `SAG2`, `DPN1`, `SYM2`, `LRP1`, `MIR2` (`/tmp/if2453/c78/meta/`, `DPN1` after its second
run of 22:59), `MV2`, `MVEUKV`, `KAVHIFI`, `INPH2604`, `GWV1`, `SAG3`; `ASG1` and `MGE4` had no metadata on the pod.
`redistributable: true` (the Zenodo record states `cc-by-4.0`, read by card C-72 part 1; DEC-80 (3)). `WET1`'s
`records.tsv` (142,712,258 B, over 64 MB) was moved, as `GWV1`'s, to `/stor/scratch/Wilke/if2453/bio-box-staging/WET1/`
(`cp -n`, sha256 equal, `/tmp` copy removed).

### 5.5 ANI against MAG stages on the pod and r226 genomes of the same genera (recorded, not cut)

`c81_gwv1_sag3_wet1.py wet1-ani` (*measured*, `/tmp/if2453/c81/work/wet1_ani/summary.json`, `wet1_ani_hits.tsv` and
one skani table per group, 2026-11-27 00:55 to 01:24): one FASTA per `WET1` MAG as queries; each reference group in
its own `skani dist -t 8 --min-af 50` process (default sketch), so that no process holds every sketch; a MAG meets when
ANI >= 96 over an aligned fraction >= 50 % (the larger of the two, card C-66's rule); the best hit over groups is kept.
The 4,164 MAGs carry 1,636 GTDB genera (MAG_meta's GTDB-Tk names).

| Reference group | How chosen | Genomes | MAGs meeting (in the group) | Best hit there | skani CPU-s |
|---|---|---:|---:|---:|---:|
| `NX1` | whole | 896 | 3 | 1 | 90.5 |
| `KAVHIFI` (bacterial MAGs) | whole | 2,214 | 0 | 0 | 93.1 |
| `DPN1` (DPANN genomes) | whole | 1,840 | 0 | 0 | 78.0 |
| `GART1` | GlobDB r226 genus (`non_gtdb_stages.tsv`) one of `WET1`'s | 889 | 6 | 6 | 105.6 |
| `GART2` | same | 1,002 | 2 | 2 | 71.3 |
| `GRARE1` | same | 369 | 17 | 13 | 73.1 |
| `GRARE2` | same | 303 | 8 | 8 | 76.4 |
| `GRARE3` | same | 645 | 20 | 20 | 81.7 |
| GlobDB r226 archaea (card C-75's set) | GTDB genus one of `WET1`'s | 2,007 | 36 | 29 | 92.9 |
| **any** | | 10,165 | | **79 of 4,164 (1.9 %)** | 762.6 |

All *measured*; the percentage *derived*. **Nothing is cut by it.** Reading: few wetland MAGs fall in a species of the
r226 genomes the pod holds (the r226 genomes of the same genera that are not on the pod were not compared, §5.3), which
agrees with the publication's 3,639 of 4,164 putative novel species against GTDB (*published*).

**`LR2` is not among the references, because its staged FASTA is damaged** (found here; not this card's to mend). The
first two `wet1-ani` runs stopped in `LR2`'s file, first on a sequence byte `0xe5`, then on a header that is not UTF-8.
`/stor/scratch/Wilke/if2453/bio-box-staging/LR2/lr2_compendia.fna.gz` has the sha256 card C-66 recorded
(`89dba3ba4098796fd3042b132089b9ea5393949d48dc7bf2783912e0b1eef377`, *measured* by `sha256sum`), so it is the file as
staged. A read-only scan (`/tmp/if2453/c81/bin/scan_lr2.py`, unit `c81-lr2-scan`, 1 CPU, 00:55 to 01:04; output
`/tmp/if2453/c81/logs/lr2_scan.out`) found, of 575,827 records: **20,526 headers that are not UTF-8, 27,396 records with
bytes that are not letters, and 26,484 records whose `contig_uid` is in no genome of `LR2`'s `records.tsv`**. All of
them start at record 504,175 (the last good header before it is a hex `contig_uid`, `d497c8fd...`). The headers and
sequences from there on look like compressed binary data written into the FASTA. `LR2`'s `records.tsv` lists 555,301
contigs; the FASTA holds 575,827 records. `WET1`'s ownership test (§5.4) read `LR2`'s `genomes.json`, which is intact,
not this FASTA. `LR2`'s build would stop at ingest on these records (*inferred*: `biobox.fasta` decodes every header).

## 6. Pod compute used

Detached `systemd-run --user` units at `nice -n 19`, `MemoryMax=32G` (8G for the two small ones), `MemorySwapMax=0`,
inside `/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif`, each started after a look at the active units (at most
one other staging unit, `c70-mir1-stage2`, so never three in all; none while `MV1` or `V3` called genes: the build of wave-6 was in
`PL2`'s rehearsal) and with no `STOP` file in the root:

| Unit | What | Cap | Wall (*measured*, file times) | CPU |
|---|---|---|---|---|
| `c81-inspect` | read-only look at both sets (§2.1, §3.1) | `CPUQuota=100%` | 18:58:14 to 18:59:33 | at most 1.3 CPU-minutes (*derived*: wall × 1) |
| `c81-xlsx` | Table S1's sheets as TSV | `CPUQuota=100%` | 4.9 s | 4.9 CPU-seconds (*measured*, `systemd-run --wait`) |
| `c81-sag3-census` | `sag3-census` | `CPUQuota=100%` | 17.2 s | 16.8 CPU-seconds (*measured*) |
| `c81-gwv1` | `gwv1-stage` | `CPUQuota=200%` | 19:14:45 to 19:30:25 (15.7 min) | about 16 CPU-minutes (*measured*: 14.7 at 14 min, `systemctl status`) |
| `c81-gwv1-ani`, `c81-sag3` (first tries) | refused at once (a stale tree copy; Table S1's repeated GTDB-Tk names, §7) | | seconds | negligible |
| `c81-sag3` | `sag3-stage` (3 code-screen workers, DuckDB 2 threads) | `CPUQuota=400%` | 19:34:57 to 19:46:39 (11.7 min) | at most 47 CPU-minutes (*derived*: 11.7 × 4); 8.9 at 3.4 min (*measured*) |
| `c81-gwv1-ani` | `gwv1-ani`, 21 skani processes of 8 threads | `CPUQuota=800%` | 19:37:33 to 20:26:53 (49.3 min) | skani 80.8 CPU-minutes (*measured*, children's rusage), largest process 27.3 |
| `c81-sag3-ani` | `sag3-ani`, one skani process of 6 threads | `CPUQuota=600%` | about 19:48 to 19:59:34 | skani 5.9 CPU-minutes (*measured*) |
| `c81-wet1-xlsx` | `MAG_meta.xlsx` as TSV | `CPUQuota=100%` | 1.8 s | 1.8 CPU-seconds (*measured*) |
| `c81-wet1` | `wet1-stage` (3 code-screen workers, DuckDB 2 threads) | `CPUQuota=400%` | 23:07 to about 00:01 (54 min) | at most 216 CPU-minutes (*derived*: 54 × 4; several processes) |
| `c81-wet1-ani`, `-ani2`, `-ani3` | stopped in `LR2`'s damaged FASTA (§5.5): a byte `0xe5`, a record in no genome, a header not UTF-8 | `CPUQuota=800%` | 00:01 to 00:16, 00:19 to 00:34, 00:36 to 00:52 | 12.9 CPU-minutes at 00:11 for the first (*measured*); the others alike |
| `c81-lr2-scan` | read-only scan of `LR2`'s FASTA | `CPUQuota=100%` | 00:55 to 01:04 | at most 9 CPU-minutes (*derived*) |
| `c81-wet1-ani4` | `wet1-ani` without `LR2`, 9 skani processes of 8 threads | `CPUQuota=800%` | 00:55 to 01:24 | skani 12.7 CPU-minutes (*measured*), largest process 1.8 |

No process came near 300 CPU-minutes. The card's estimate (8 threads for skani, 2 to 4 h, *estimate*, with an r226
stream that was not run, §5.3): about 1.7 h of wall for `GWV1` and `SAG3` (*derived*: 18:58 to 20:27, units
overlapping) and about 2.5 h for `WET1` (*derived*: 22:55 to 01:24, of which 50 min were the three stopped ANI runs). The wave's `PL3` Pfam build ran beside the `WET1` units; at most one unit of
this card ran at a time then, and no other staging unit was active.

## 7. Tests

`tests/r1/test_c81_gwv1_sag3_wet1.py`, constructed truth throughout: contig uids by blake3 in the test; GWV1 fixture
records built from a known core and repeat (an exact 30 bp DTR trimmed and circular; a DTR-flagged record with no
repeat circular and untrimmed; the same kind of repeat on a record CheckV calls AAI-based left linear and whole; ITR
linear; a repeated contig set listed once; contigs held by an owner and by a later stage); refusals that leave nothing
behind (a record the table lacks, a table row with no record, a length that differs, an existing output, a later stage
without `genomes.json`); the vOTU rules and the sample hits on a hand-written skani table; the chunked FASTA reader;
skani itself on 20 kb fixtures (skipped without skani); a Table S1 written here as raw XML (three sheets in a workbook
order unlike their file numbers, shared and inline strings, a GTDB-Tk name given twice as in the real table, a
Sequencing sheet whose column moved, refused); the member-to-name match (longest name, ties match nothing); the tar
census; the SAG3 stage end to end twice, byte-identical, with the r232 status through a hand-written GTDB metadata
table and C-74 index (via WGS, via BioSample, absent, not visible), the code screen (a *Mycoplasma* hint gives code 4),
the FASTA in archive order with renamed contigs, owner and later-holder tables and `s2_check_inputs.py` exit 0; the
r226 archaea picked by genus; per-genome FASTA files; the best ANI hit at the thresholds; and both templates through
`stage_est.read_conf` with their `SUM`s found in this note and their `EXPECT` holders among the root, `r26.11.4` and
`MV2`. Every test runs with the normal `PATH` and with `PATH=/usr/bin:/bin`.

Seen failing first on the pod (*measured*, `wilkcomp02`, 2026-11-26, `scripts/pod_pytest.sh C-81
tests/r1/test_c81_gwv1_sag3_wet1.py`, with the script not yet in the tree): **37 failed, 1 skipped** (no
`c81_gwv1_sag3_wet1.py`, no templates). The six-pair cut of the GTDB-Tk sheet was found by the first run on the real
table (`refused: ... Outputs_GTDB-TK names AB-754_C22 twice`); the fixture was given that case first and the SAG3
tests seen failing (16 failed: 12 because the stager refused the fixture's table, 4 the templates not yet filled)
before the reader was changed. The `WET1` tests (MAG_meta under its title row, the zip census, the BioProject
crosswalk, the stage end to end twice with `s2_check_inputs.py`, a member the table lacks refused, genus selection,
`merge_best`, the template) were written after card R-18 and seen failing first (*measured*, 2026-11-26: **14
failed**, `module 'c81_gwv1_sag3_wet1' has no attribute 'mag_meta'` and the like); then the MIMAG tier
(`mag_hq`/`mag_mq`, which `s2_check_inputs.py` requires; 6 failed first), a `genome_id` column (cards C-58, C-66; 2
failed first) and a reference record holding the byte `0xe5`, which the first `wet1-ani` run met in `LR2`'s FASTA
(written now as bytes, never decoded; 2 failed first), and a reference record in no genome, counted instead of
refused when the caller asks (2 failed first). Final
(*measured*, same host and day, `env BIOBOX_TREE=/tmp/biobox-wt/C-81 sh scripts/pod_pytest.sh C-81
tests/r1/test_c81_gwv1_sag3_wet1.py tests/r1/test_r232_index.py tests/r0`, one process): **285 passed, 1 skipped in 26.36 s** for `GWV1` and `SAG3`; with `WET1`, after the merge of card R-18
(2026-11-27): **307 passed, 1 skipped in 34.71 s** (the skip: the skani test under the tools-hidden `PATH`).
