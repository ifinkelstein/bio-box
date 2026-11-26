# Sources of diversity the earlier surveys missed, and a review of the `r26.11.5` roadmap (card C-72)

- **Date:** 2026-11-26. Every record below was read on this date unless another date is given.
- **Asked by the PI (2026-11-26):** "separately run an agent to explore additional sources of diversity that proir explorations and
  current ingests may have missed. be deep and thorough. where else can we get sources of new and super diverse metagenomes? higher
  quality, no login preferred bu not essential. i can log in as needed." Added the same day: "review roadmap as well. r26.11.5 should be
  modified if higher diversity/priority datasets can be staged for ingest".
- **Method:** read-only. Web search; providers' own record pages and API listings (GTDB and GlobDB directory indexes and release
  notes, the EBI MGnify FTP tree and its READMEs, the Zenodo records API, the figshare articles API, Europe PMC full text of open
  articles); two `ls` calls on `tower` (`sources_incoming/` and `sources/`). Nothing was downloaded beyond directory indexes, README
  text and API metadata; no request carried a name, an address or a detail of this project; no login, token or cookie was used (the
  `nature.com` article pages answer with a redirect to a cookie-setting `idp.nature.com` flow, which was not followed; Europe PMC served
  the same articles' text).
- **What was read first** (the job was what these missed): `planning/CATALOG_GOAL.md`; the four notes of 2026-11-25
  (`r26_11_5_diversity_survey`, `r26_11_5_deep_collections`, `r26_11_5_wide_and_audit_wishlist`, `pi_login_pulls`) and
  `2026-11-25_tower_sources_incoming.md`; `2026-11-18_globdb_inputs.md`; `sources.yaml` (161 entries); cards C-55 to C-71; the top
  two sections of `planning/tasks/STATUS.md`; `planning/DECISIONS.md` rows DEC-67 to DEC-77; the stage confs and templates.
- **Labels (AGENTS §6):** *measured* (where and how read), *published* (the provider's own statement on its record page),
  *derived* (formula given), *estimate*. Byte counts copied from an API or directory index are *measured* (read from that listing on
  2026-11-26).

## 0. Summary

1. **The largest miss is a release, not a paper: GlobDB r232 (June 2026) and GTDB R232 (April 2026).** Every plan in the repository
   dereplicates against, and cuts tranches from, GlobDB r226 and GTDB R226. GlobDB r232 holds 346,233 species representatives from 26
   resources (r226: 306,260 from 14); 46,920 of its genomes come from 12 resources that r226 did not include, led by gcMeta 2025
   (24,063) and the CNCB-NGDC Genome Warehouse direct submissions (13,287) (*derived* from its `globdb_r232_dataset_list.tsv`, §2.1).
   GTDB R232 has 199,923 species clusters, 56,309 more than R226 (*derived*, §2.2). Both are anonymous downloads of a few large
   files. The note of 2026-11-18 mentions that r232 exists (one clause, on a QC exemption) but no card, survey or remainder row sizes
   it or plans for it.
2. **Six new open sources, each one file or one FTP tree:** the Groundwater Virome Catalogue (280,420 vOTUs, one 1.17 GB file,
   CC BY 4.0); the OMZ single-cell compendium (3,570 SAGs, one 747 MB file, CC BY 4.0); the GWM global wetland catalogue (4,164
   dereplicated MAGs, 3,639 putative novel species, one 3.89 GB zip, CC BY 4.0); MGnify `marine` v3.0 (released 2026-11-09: 11,313
   species new since v2.0); MGnify `marine-sediment` v1.0 (6,158 species); MGnify `marine-eukaryotes` v1.0 (74 species, for card C-60).
3. **Five duplication catches against GlobDB:** the glacier-fed-stream assemblies `GFS1` are the project that GlobDB r226 already
   carries as `GFSMAG` (1,377 non-GTDB species on `/stor`); and `CFMD1`, `CRBC_ISO`, `QXLSG1` (K5) and several smaller catalogues
   proposed for separate fetches are now inside GlobDB r232 at species level (§2.7).
4. **Roadmap (§5):** the order in the handoff of 2026-11-26 07:15 leaves out eleven sets that are already on `tower` or on the pod and
   that the cards or the survey put in `r26.11.5` (`GV1`, `EGV1`, `ECE2`, `MIR1`, `ASG1`, `SAG2`, `LRP1`, `MIR2`, `SYM2`, `DPN1`,
   `MGE4`), and it puts in sets that cannot be ready without thousands of fetches (`SPX2`, `O4`, `RMD1`, `LGRT2`'s run fetches). The
   proposed `r26.11.5` keeps everything that is on disk and high in novelty, adds those eleven and four of the new finds, moves seven
   members to `r26.11.6`, and comes to about 108 to 120 Gbp (*derived* from estimates, §5.4). **It is a proposal for the PI.**

## 1. The scales used for every candidate

Novelty per Gbp, against what is held (built, staged, landed on `tower`, or queued there). The bands are anchored on the measured
yields in `CATALOG_GOAL.md` §2 (new NR50 clusters per Gbp: huge phages about 489,000, SPIRE aquatic 465,323, terrestrial 407,439,
host-associated 269,797, isolate plasmids about 34,000; GTDB complete strain genomes about 27,000 by DEC-70). Placing a candidate in a
band is an *estimate* unless the provider states its overlap.

| Band | New NR50 per Gbp (*estimate*) | Typical members |
|---|---|---|
| N5 | above 400,000 | divergent viruses outside gut, rare-phylum and unclassified MAGs, dark-ocean SAGs, unbinned contigs of thin habitats |
| N4 | 250,000 to 400,000 | environmental MAG catalogues of species absent from GlobDB, long-read environmental MAGs |
| N3 | 100,000 to 250,000 | catalogues that partly overlap GlobDB, host-associated MAGs, compact curated MGEs of known families |
| N2 | 30,000 to 100,000 | complete isolates of new species in well-sampled clades, isolate plasmids and phages |
| N1 | under 30,000 | strain copies, human-gut species already held |

Ease of staging:

| Class | Meaning |
|---|---|
| E1 | on the pod now (`/stor` or `/tmp`), staged or in an archive already verified there |
| E2 | landed on `tower` |
| E3 | queued in a `tower` batch, or an anonymous fetch of at most ten files under about 30 GB |
| E4 | anonymous, but one stream over 100 GB or hundreds to thousands of files, or a selection that needs a new script |
| E5 | login, access request, embargo, web-application firewall, licence not stated, or data location not verified |

## 2. Findings by axis

Each record gives: stage name proposed; what; quality; novelty against held; access; licence as the record states it; effort; grab
priority and release. "Held" means the stages and fetches listed in the card's brief and in `CATALOG_GOAL.md` §3.

### 2.1 Repositories: GlobDB r232 (`G232-A`, `G232-B`, and the `G2..G41` remainder re-cut)

- **What:** GlobDB release 232, "released in June 2026", "based on the GTDB release 232 and contains 346,233 genomes, from 26 data
  sources" (*published*, `https://globdb.org/` and `https://globdb.org/news` post of 2026-06-25). Sequential species dereplication as
  in r226; the news post of 2026-03-07 says genomes from additional sources are added "if they cover species not already in INSDC
  databases or earlier GlobDB releases" (*published*).
- **Composition** (*measured*, `https://fileshare.lisc.univie.ac.at/globdb/globdb_r232/globdb_r232_dataset_list.tsv`, columns
  "species reps in dataset" and "number of genomes in GlobDB r232"): GTDB 199,923; mOTU 42,566; SPIRE 29,737; **GCMETA 24,063** (of
  143,000 in gcMeta 2025, NAR `doi:10.1093/nar/gkaf1115`); **NGDC GWH 13,287** (of 30,433); RBG 6,780; TPMC 5,393; GEM 4,680; MGnify
  4,187; **TPMCS 3,064** (Tibetan Plateau sediments, Asgardarchaeota-rich, `doi:10.34133/research.1213`); GOMC 2,500; SMAG 2,273;
  **QXLSG 1,635**; **CRBC 1,572**; **HOGU 1,398** (horse gut); **TPLM 1,117** (Tibetan lakes); HRGM2 472; AMXMAG 465; **TG2G 277**
  (glaciers); **PREC 261** (Pearl River estuary); cFMD 234; **DAWW 147**; **CRLG 82** (rivers); MRGM 59; SHGO 44; **SCSSF 17**. The
  twelve in bold are not among r226's thirteen MAG resources (`2026-11-18_globdb_inputs.md` §2.2 table; `sources.yaml`
  `globdb_r226.records`). Their genomes in r232: **46,920** (*derived*: 24,063 + 13,287 + 3,064 + 1,635 + 1,572 + 1,398 + 1,117 + 277
  + 261 + 147 + 82 + 17). r226's `GFS` (1,377) is not in r232's list.
- **Size:** `globdb_r232_genome_fasta.tar.gz` 292,068,313 K in the index (*measured*, directory index, 2026-06-22), that is about
  299.1 GB (*derived*: × 1,024, the index's K). Whole release about 1.04 Tbp (*estimate*: r226's 875.80 Gbp per 251,371,669,660 bytes =
  3.484 bp per byte, × 299.1 GB). The 46,920 new-resource genomes: about 118 Gbp (*estimate*: × 2.52 Mbp, r226's non-GTDB mean, 410.33
  Gbp / 162,653). Small tables: `globdb_r232_taxonomy.tsv.gz` 8,071 K, `globdb_r232_checkm2.tsv.gz` 12,896 K,
  `globdb_r232_dictionaries.tar.gz` 3,266 K, `globdb_r232_dataset_list.tsv` 4 K (*measured*, index).
- **Quality:** CheckM2 table for every genome (`globdb_r232_checkm2.tsv.gz`); r226's floor was completeness 40% and r232 adds an
  explicit GTDB exemption (`2026-11-18_globdb_inputs.md` line 140). MIMAG tiers not tabulated for r232 here (not read).
- **Novelty:** the 46,920 genomes of the new resources are, by the inclusion rule above, species not in GTDB R232, not in the
  resources dereplicated before them, and not in earlier GlobDB releases (*published* rule; the count is *derived*). Against what is
  held (GlobDB r226's 162,653 non-GTDB genomes on `/stor`, of which `G1`, `GART1`, `GART2` are staged or built) they are new at species
  level. N4 overall; N5 for the archaea, the unclassified-domain genomes and the phyla absent from held stages (*estimate*).
  GCMETA and NGDC GWH are mostly Chinese soil, water, sediment and animal-gut projects that no earlier survey sized.
- **Access:** `https://fileshare.lisc.univie.ac.at/globdb/globdb_r232/{file}`, anonymous HTTPS (the r226 archive came from the same
  host in one stream of 3 h 46 min, `sources.yaml` `globdb_r226.fetched`). One genome tarball; ten top-level files in all.
- **Licence:** "the GlobDB propagates the licenses of the underlying data sources, and is licensed under CC BY-SA 4.0" (*published*,
  `globdb.org`, read 2026-11-26). The underlying resources' own terms are not checked here; GWH direct submissions are INSDC-like
  public records (not verified per record).
- **Effort:** tables E3 (four small files); genomes E4 (one 299 GB stream with member selection, as `c41_stream_stage.py` did for r226).
- **Grab:** tables now (they are needed to dedupe every other pick in this note and in the 2026-11-25 notes); `G232-A` (the
  new-resource genomes in archaea, unclassified domain and phyla or classes absent from held stages, capped near 12 Gbp) in
  `r26.11.5` if the PI approves the stream; `G232-B` (the rest of the 46,920, about 105 Gbp *estimate*, in two or three tranches) from
  `r26.11.6`. The 130,960-genome `G2..G41` remainder from r226 should be re-keyed to r232's IDs before any further `G` tranche, so that
  genomes r232 replaced with a better GTDB R232 representative are recorded once.

### 2.2 Repositories: GTDB R232 (metadata first; base of `S3`, `S5`, `S6` and of the `SYM1` and `REFSEQ-NEW` lists)

- **What:** GTDB R11-RS232, announced 2026-04-15: 901,341 genomes, 199,923 species clusters ("a 39% increase") (*published*, GTDB
  forum `https://forum.gtdb.ecogenomic.org/t/announcing-gtdb-r11-rs232/826`; release notes
  `https://data.gtdb.aau.ecogenomic.org/releases/release232/232.0/RELEASE_NOTES.txt`). Bacteria 878,998 genomes and 189,801 species;
  Archaea 22,343 and 10,122; 186 phyla (*published*, `https://gtdb.ecogenomic.org/stats/r232`). Against R226 (143,614 species; 136,646
  bacterial, 6,968 archaeal): **56,309 more species clusters, 3,154 of them archaeal** (*derived*: differences).
- **Size:** `gtdb_genomes_reps_r232.tar.gz` 178.97 G; `bac120_metadata_r232.tsv.gz` 275.7 M; `ar53_metadata_r232.tsv.gz` 6.7 M
  (*measured*, directory indexes of `release232/232.0/` and `genomic_files_reps/`).
- **Novelty:** most of the new species are MAGs that GlobDB r226 already carried as non-GTDB genomes; GlobDB r232 contains every GTDB
  R232 representative, so the genomes are best taken from the GlobDB r232 stream, not from GTDB's tarball. N2 to N3 for the isolate
  part (*estimate*). The forum reports a fault in R232's species boundaries (WitChi filtering pulled genomes out of *E. coli*; a fix is
  promised for the next release) (*published*, forum thread 827).
- **Access / licence:** anonymous HTTPS; GTDB data are CC BY-SA 4.0 (`sources.yaml` `gtdb_metadata`, not re-read for R232).
- **Effort:** E3 for the two metadata tables.
- **Grab:** the metadata in `r26.11.5` preparation: it re-keys `SYM1_post_r226_left_for_later.tsv.gz` (49,640 rows) and the
  `REFSEQ-NEW` list (5,036 complete representatives absent from R226), and every "absent from GTDB R226" claim in the 2026-11-25 notes.

### 2.3 Repositories: MGnify genome catalogues released after the earlier surveys' list (`MGM3`, `MGS1`, `MGE-EUK1`)

The FTP tree `https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/` lists twenty catalogues (*measured*, index). The 2025
survey's `MGNIFY2` named marine v2.0, soil v1.0, rumen, pig, chicken, HoloFood, rhizosphere, fish and honeybee, and its batch kept only
ten metadata tables. Not named anywhere in the repository:

| Proposed stage | Catalogue (directory, date in the index) | Contents (*published*, its README) | Novelty (*estimate*) | Grab |
|---|---|---|---|---|
| `MGM3` | `marine/v3.0/` (2026-11-09) | 215,647 genomes in 24,524 species; **164,813 MAGs added, 11,313 new species** since v2.0; genomes from 1,971 studies; `added_genomes.tsv` 4.4 M, `genomes-all_metadata.tsv` 84 M | N3 to N4: v3.0 is two weeks old; its new species overlap GOMC, OMD v2 and SPIRE in unknown part (30 to 60% new vs GlobDB r232, *estimate*); about 28 Gbp for the 11,313 (*estimate*: × 2.5 Mbp) | `r26.11.6`; E4 (one `.fna.gz` per species directory under `species_catalogue/`) |
| `MGS1` | `marine-sediment/v1.0/` (2025-11-18) | 10,126 genomes in 6,158 species | N4: sediment is thin in held stages; about 15 Gbp (*estimate*) | `r26.11.6`; E4 (6,158 files) |
| `MGE-EUK1` | `marine-eukaryotes/v1.0/` (2026-10-27) | 385 eukaryotic genomes in 74 species from 68 studies | N5 but blocked on card C-60 (eukaryote kind, intron-aware calling) | with C-60; E3 (74 files) |
| (none) | `human-gut/v3.0/` (2026-10-21) | 286,621 genomes, 6,461 species, 1,864 new species | N1 to N2 | later, strain pass |
| (none) | `barley-rhizosphere/v2.0/` (2026-01-09), `maize-rhizosphere`, `tomato-rhizosphere` (2025-08-29), `human-skin` (2025-07-25), `human-skin-eukaryotes` (2026-10-28) | small (barley `genomes-all_metadata.tsv` 43 K) | N3 (rhizosphere) | fold into `PLNT1`'s dedup; later |

Licence: none of the READMEs read states one; EMBL-EBI's terms of use apply (not read). Record as "not stated in README".

### 2.4 Environments

**Groundwater and aquifers: the Groundwater Virome Catalogue (`GWV1`).**
- Wu, Liu et al. 2024, *Nat. Commun.* 15, "Unveiling the unknown viral world in groundwater" (`doi:10.1038/s41467-024-51230-y`,
  `PMC11310336`): **280,420 species-level vOTUs of 5 kb and more** from 607 monitored wells in seven geo-environmental zones; 312,741
  viral contigs before 95% ANI clustering; CheckV 5,366 complete, 6,092 high-quality and 15,669 medium-quality; 107,610 vOTUs of 10 kb
  and more; lengths 5 kb to 543.1 kb; "over 99%" of the viruses and about 95% of viral clusters novel (*published*, Europe PMC full text).
- Zenodo `10.5281/zenodo.11230969` "The Groundwater Virome Catalogue (GWVC)", published 2024-05-21, `cc-by-4.0`, open; eleven files,
  `GWVC.vOTUs.fna.gz` 1,166,850,604 bytes, `GWVC.vOTUs.information.txt` 62,086,246, `GWVC.vOTUs.virus-host.linkage.txt` 59,543,590,
  `CPR.DPANN.genomes.tar.gz` 43,339,586 (*measured*, `https://zenodo.org/api/records/11230969`).
- About 3.6 to 4.0 Gbp (*estimate*: gzip FASTA at 3.1 to 3.4 bp per byte). N5: groundwater viruses are absent from every held viral
  stage except as fragments of `V2`/`MV1` (IMG/VR v4.1 and MetaVR draw on public metagenomes; the reads are `PRJNA858913`, released with
  this paper, so overlap with MetaVR is possible and is the first thing to measure, *estimate* under 10%). Hosts span 119 phyla
  including CPR and DPANN.
- The 2025 survey's "Wave B" names `GPDC`/`PRJNA858913` as an aquifer *defensome*; the virome deposit itself is not named in any note.
- Effort E3 (one sequence file). **Grab: `r26.11.5`.** `DEPTH=called`, `gv-meta`, DEC-65 terminal-repeat handling.

**Wetlands, paddies and peat: the GWM global wetland resource (`WET1`).**
- "An integrated global resource of wetland microbiomes linking environmental metadata, community profiles, and genome-resolved
  metabolic traits", *Sci. Data*, published 2026-06-04 (`doi:10.1038/s41597-026-07581-w`, `PMC13547364`): 1,962 samples (129 newly
  sequenced; 567 soil, 194 sediment, 1,201 water); **5,704 MAGs from 251 samples (187 high-quality, 5,517 medium-quality; QS =
  completeness − 5 × contamination > 50), dereplicated to 4,164 genomes; 3,639 putative novel species by GTDB-Tk; 3,671 of 4,164
  lacking a species match in GTDB** (*published*, Europe PMC full text). MAGs also in NCBI `PRJNA1430164`.
- Zenodo `10.5281/zenodo.20263536`, published 2026-05, `cc-by-4.0`, open: `Dereplicated genomes.zip` 3,885,785,718 bytes, `Genome
  sequences.zip` 5,289,781,056, `MAG_meta.xlsx` 806,764, and annotation archives (*measured*, Zenodo API). The article itself is
  CC BY-NC-ND 4.0; the data record says CC BY 4.0.
- About 12 to 13 Gbp for the dereplicated set (*estimate*: 4,164 × 3 Mbp). N4: novelty is stated against GTDB, not GlobDB; part of its
  251 samples are public paddy and wetland runs that SMAG or SPIRE may have binned (*estimate*: 50 to 70% new species vs GlobDB r232).
- Effort E3 (one zip). **Grab: `r26.11.5`**, `skani` 96% against GlobDB r226 on disk and r232's table.

**Oxygen-minimum zones (water column): see §2.6 (`SAG3`).**

**Marine sediment, estuaries: MGnify `marine-sediment` v1.0 (§2.3, `MGS1`).** Also found and not taken: Maeke et al. 2025 coastal
lagoon, salt-marsh and estuary dataset (1,372 MAGs, *published* per search abstract, record not opened); Pearl River estuary `PREC`,
South China Sea cold seep `SCSSF` and the river catalogue `CRLG` are inside GlobDB r232.

**Wastewater: global activated sludge (`WWS2`, unverified).** "Metagenome-resolved global microbial diversity and function in
activated-sludge wastewater treatment systems", *Nature Water*, 2026-01-23, `doi:10.1038/s44221-025-00576-8`: more than 24,500 genomes
from almost 8,000 species from 125 new and 703 public sludge samples (*published*, NTU news page). Europe PMC has no record and the
news page names no repository. **Unverified; overlaps MiDAS (`MIDAS_LR`) and Singleton 2021 (`LR2`).** Later.

**Lakes.** A 2026 conference abstract reports "over 72,000 freshwater bacterial and archaeal MAGs" from a 10-year single-lake series
(*published*, `meetingorganizer.copernicus.org/WBF2026/WBF2026-597.html`); no release is stated: watch, not grab. Buck et al. 2021
(*Sci. Data* 8:131; 12,665 MAGs, 3,640 mOTUs from 41 stratified boreal and subarctic lakes) was not in any note; its MAGs are in ENA
and probably in mOTUs v4 and GTDB R232 (*estimate*, not checked); low priority. Tibetan lakes (`TPLM`, `QXLSG`) are in GlobDB r232.

**The other environments of the brief** are covered in §3 (what was looked for and not found) or by the 2026-11-25 notes.

### 2.5 Lineages

No new lineage-specific deposit beyond the 2026-11-25 notes was found that is open and sized. Two indirect gains: GlobDB r232's
`TPMCS` (3,064 genomes from Tibetan Plateau sediments, described by its paper as rich in *Asgardarchaeota*) and GTDB R232's 3,154
new archaeal species; both bear on `ASG1`, `ASG2` and `ARC1`, which should be deduplicated against r232, not r226. The GWVC deposit
carries `CPR.DPANN.genomes.tar.gz` (43 MB) with CPR- and DPANN-linked viruses.

### 2.6 Data types: single cells (`SAG3`)

- "A compendium of bacterial and archaeal single-cell amplified genomes from oxygen deficient marine waters", *Sci. Data* 2023
  (`doi:10.1038/s41597-023-02222-y`, `PMC10224968`): **5,129 SAGs, of which 3,570 sequenced, assembled and decontaminated**, from
  Saanich Inlet, the Northeastern Subarctic Pacific, the North Pacific and South Atlantic Subtropical Gyres, the Eastern Tropical North
  and South Pacific, the Eastern South Pacific coastal upwelling, Benguela and the Baltic Sea (*published*, Europe PMC full text).
- figshare `10.6084/m9.figshare.20459526` "File 5. Single Cell Assemblies from the OMZ Compendium", published 2023-04-06, licence
  "CC BY 4.0", one file `OMZ_Compendium_SAGs_Full.tar.gz` 747,416,178 bytes (*measured*, figshare API). Metadata with IMG, BioSample and
  GenBank accessions: figshare `10.6084/m9.figshare.20481603` (Table S1, not opened).
- About 2.3 to 2.5 Gbp (*estimate*: × 3.1 to 3.4). N4 to N5: SAGs of oxygen-deficient waters are not in `SAG1` (GORG-Tropics, sunlit
  tropics) or `SAG2` (GORG-Dark, aphotic oxic ocean); part of the compendium is in GenBank and so in GTDB R232 (*estimate*: 20 to 50%
  of SAGs already represented at species level). Dedup by Table S1's GenBank accessions and `skani`.
- Effort E3. **Grab: `r26.11.5`.**
- Microbe-seq (Zheng et al. 2022, *Science*: 21,914 SAGs co-assembled into 76 species from one person's gut) was looked at: human gut,
  N1; not proposed.

### 2.7 Duplication catches against GlobDB (affecting sets already landed, queued or planned)

| Set | Where | Overlap found | What to do |
|---|---|---|---|
| `GFS1` (Vanishing Glaciers, `PRJNA781406`, 3,387 NCBI assemblies, 12.01 Gbp, batch 4) | `tower` manifest | GlobDB r226 already holds this project's MAGs as `GFSMAG`: 1,377 non-GTDB species, 4.43 Gbp, on `/stor` (*measured*, `2026-11-18_globdb_inputs.md` §2.2) | take only genomes whose species `GFSMAG` and GTDB lack; the 14 eukaryotic assemblies over 20 Mb (handoff) are a second filter |
| `CFMD1` (cFMD, 11.23 GB landed) | `tower` | GlobDB r226 holds 336 and r232 234 cFMD species (*measured*, both tables) of cFMD's 1,144 SGBs | species-level dedup; most of the 10,899 MAGs are strain copies (N1 to N2) |
| `CRBC_ISO` (crop root, batch 6b) | queued | GlobDB r232 holds 1,572 of CRBC's 2,330 MAG species | the isolates may be new; the MAGs mostly are not |
| `QXLSG1` (K5, Qinghai-Xizang lakes) | proposed 2026-11-25 | GlobDB r232 holds 1,635 of its 2,742 species | take from the r232 stream, not separately |
| `TG2G`, `TPLM`, `TPMCS`, `PREC`, `HOGU`, `SCSSF`, `DAWW`, `CRLG` | not yet proposed elsewhere | all inside GlobDB r232 | take from the r232 stream |

### 2.8 Programmes

- **Earth Microbiome Project / EMP500:** raw metagenomes under ENA `PRJEB42019`; no MAG release found. Nothing new.
- **Australian Microbiome:** the project page lists an Australian MAG database as a goal; no release found. Nothing new.
- **Bio-GO-SHIP / bioGEOTRACES:** 971 surface-ocean metagenomes (Larkin et al. 2021, *Sci. Data*), no MAG catalogue found; they reach
  the corpus through SPIRE, OMD and Logan. Nothing new.
- **NEON via NMDC:** NMDC's portal reports 95 studies, 16,361 metagenomes, 411 metatranscriptomes and 153,411,782,343,495 bytes of data
  (*measured*, `https://data.microbiomedata.org/api/stats`); bins per study were not counted. The 2025 note's `NMDC-NEON` (K10) stands.
- **gcMeta and CNCB-NGDC GWH:** reach the corpus at species level through GlobDB r232 (§2.1); a separate fetch would add strain
  copies only.
- **Tara successors, Malaspina, Mission Microbiomes, Earth BioGenome side products, Microflora Danica, national soil and water
  programmes:** nothing beyond the 2026-11-25 notes (`RMD1`, `OMD2`, `DTOL1`, `MFD_LR`, `GFS1`, `SODA1`) was found.

## 3. What was looked for and not found

| Axis | Searched for | Result (2026-11-26) |
|---|---|---|
| Deep subsurface, aquifers | global groundwater MAG catalogue 2025 | none; a 2025 *ISME Commun.* synthesis of about 4,000 existing deep-subsurface MAGs and SAGs (not a new deposit); the groundwater **virome** was found (`GWV1`) |
| Polar, sea ice | Southern Ocean or Antarctic sea-ice MAG catalogue | none beyond 7 MAGs (*MRA* 2025); MOSAiC (`MOS1`) already held |
| Marine sediment, subseafloor | global sediment catalogue 2025 to 2026 | none outside MGnify `marine-sediment` (found) and cold-seep sets already held |
| Air, built environment | global indoor or air MAG catalogue | none; Hong Kong built environment (860 MAGs, 2024) and a five-city air virome (2025), both small |
| Insects, fish, wild vertebrates | catalogues of thousands of species | none beyond `iSymBase`, `EHI1`, MGnify fish (held or queued) |
| Ancient | ancient MAG catalogue of thousands | none; permafrost studies of 33 to 60 MAGs |
| Sponges | global sponge MAG catalogue | none beyond 1,188 MAGs (2021) |
| Lakes | global lake catalogue | none released; one single-lake 72,000-MAG series announced, not released |
| Long-read | 2026 HiFi or ONT ocean virome or circular-MAG catalogue | none beyond the 2025 notes' sets |
| Eukaryotes | 2025 to 2026 eMAG catalogue of thousands | none; MGnify `marine-eukaryotes` v1.0 (74 species) found |
| SAGs | 2026 environmental SAG release | none beyond the OMZ compendium (2023, found) |
| Soil viruses | catalogue of a million vOTUs | none beyond GSV Atlas and PIGEON (held or login) |
| Mega-resources | IMG/VR v5; SPIRE v2; PLSDB after the held version | IMG/VR's successor is MetaVR (held); no SPIRE update; PLSDB's 2025 paper (72,360 records) is the version held |
| Aggregators | `UniMAG v1.0` (*Sci. Data* 2026, 581,345 MAGs per a search abstract) | record not reachable without the `nature.com` cookie flow; Europe PMC has no entry: unverified; an aggregator, low novelty expected |

## 4. Rankings of the new finds, and a merged grab order

### 4.1 By novelty per Gbp (*estimate*, §1 bands)

1. `GWV1` Groundwater Virome Catalogue: N5.
2. `G232-A` GlobDB r232 new-resource archaea, unclassified and absent phyla: N5.
3. `SAG3` OMZ SAG compendium: N4 to N5.
4. `MGE-EUK1` MGnify marine eukaryotes: N5, blocked on C-60.
5. `MGS1` MGnify marine-sediment v1.0: N4.
6. `WET1` GWM wetlands: N4.
7. `G232-B` remaining GlobDB r232 new-resource genomes: N4.
8. `MGM3` MGnify marine v3.0 new species: N3 to N4.
9. `WWS2` global activated sludge (unverified): N3.
10. GTDB R232 isolates beyond GlobDB: N2 to N3.

### 4.2 By ease

1. GlobDB r232 tables (four files, about 24 MB): E3.
2. GTDB R232 metadata (two files, 282 MB): E3.
3. `GWV1` (one 1.17 GB file): E3.
4. `SAG3` (one 747 MB file): E3.
5. `WET1` (one 3.89 GB zip): E3.
6. `MGE-EUK1` (74 files): E3, but C-60.
7. `G232-A` / `G232-B` (one 299 GB stream with member selection): E4.
8. `MGS1` (6,158 files): E4.
9. `MGM3` (about 11,313 files): E4.
10. `WWS2`: E5 (data location not found).

### 4.3 Merged grab order

| # | Item | Release | Why here |
|---|---|---|---|
| 1 | GlobDB r232 tables + GTDB R232 metadata | `r26.11.5` preparation | tiny; needed to dedupe every pick in this note and the 2026-11-25 notes |
| 2 | `GWV1` | `r26.11.5` | highest novelty, one file |
| 3 | `SAG3` | `r26.11.5` | single cells of a habitat no stage has, one file |
| 4 | `WET1` | `r26.11.5` | 3,639 putative novel species, one file |
| 5 | `G232-A` (stream, capped about 12 Gbp) | `r26.11.5` if the PI approves the stream, else `r26.11.6` | the largest open pool of new species found |
| 6 | `MGS1` | `r26.11.6` | sediment; 6,158 files |
| 7 | `MGM3` | `r26.11.6` | after dedup against r232 and `OMD2` |
| 8 | `G232-B` | `r26.11.6`, `r26.11.7` | about 105 Gbp, two or three tranches |
| 9 | `MGE-EUK1` | with C-60 | eukaryote kind first |
| 10 | `WWS2` | later | data location to be found |

## 5. Roadmap review and proposed `r26.11.5`

### 5.1 Where the `r26.11.5` plan is written, and where those places disagree

| Place | What it says for `r26.11.5` |
|---|---|
| `STATUS.md`, handoff of 2026-11-26 07:15, H2 step 3 | "stage `r26.11.5` in the PI-approved novelty-first order (`G-RARE1..3` → Coral `RMD1`+`GBR_MGD` → `MV-EUKV`+`MV2` → `KAV-HIFI`+`INPHARED-2604`+`O4` → `SPX2`+`LGRT2` → `LR2+SAG1+NX1`+`REFSEQ-NEW`+`tower` long-read `cMAGs`)" |
| `STATUS.md`, section of 07:25 (current) | nothing on `r26.11.5`; `r26.11.4` is building `V2` alone (W0) |
| `DECISIONS.md` | no row decides `r26.11.5`'s members. DEC-67 (1) put `C1..C4 + C6` there (superseded by DEC-70); DEC-71 (2) put `MV1 + LR1 + MGE1 + SAT1/V3 + G_ART1` there (moved into `r26.11.4` by DEC-72) |
| Card C-69 | `EGV1` "for `r26.11.5`" right alongside deferred `GV1` |
| Card C-70 | `ECE2`, `MIR1`, `ASG1` "for `r26.11.5–r26.11.6`"; templates give the chain `EGV1 -> ECE2 -> MIR1 -> ASG1` |
| Card C-71 | `GHM1`, `VHM1`, `MGE3`, `LR3`, `ARC1`, `NX3` "for `r26.11.6+`" |
| `CATALOG_GOAL.md` §2 | item 3 puts `GV1`, `EGV1` and the C-70 sets in "`r26.11.5–r26.11.6`"; item 4 puts C-71's sets later |
| `2026-11-25_r26_11_5_diversity_survey.md` §3 | Tier 1: `LR2 + SAG1 + NX1 + GV1` (1a), `EGV1 + ECE2 + MIR1 + ASG1` (1b), `DPN1`, `SYM2`, `LRP1`, `MIR2 + MGE4`, `O3`; Tier 2: `SAG2`, `LR-COMP`, `MFD-LR`, `MIDAS-LR`, `GHM1`, `NX2`, `OMD2`, the virome wave |
| `2026-11-25_r26_11_5_wide_and_audit_wishlist.md` | zero-fetch first: `LR2 + SAG1 + NX1`, `REFSEQ-NEW`, `G-RARE1..3`, `MV2 + MV-EUKV`, `KAV-HIFI`, `INPHARED-2604`, `J1`; then `SPX2`, `LGRT2`, `O4` |
| Stage confs | only templates exist for `EGV1`, `ECE2`, `MIR1`, `ASG1`, `ET1`, `LGNRNAP1`; `GV1.conf` is in `stages/`; no conf for `LR2`, `SAG1`, `NX1`, `SAG2`, any `G-RARE` or anything on `tower` |

Disagreements: (1) the handoff calls its order "PI-approved", but no DECISIONS row or recorded PI word in the files read approves it;
(2) the handoff leaves out `GV1`, `EGV1`, `ECE2`, `MIR1`, `ASG1`, which cards C-69 and C-70 and `CATALOG_GOAL.md` §2 put in
`r26.11.5`, and the landed compact sets of the survey's Tier 1 (`SAG2`, `DPN1`, `SYM2`, `LRP1`, `MIR2`, `MGE4`); (3) the survey puts
`LR2 + SAG1 + NX1 + GV1` first, the handoff last; (4) the Great Barrier Reef set has three descriptions: C-71 "`GBR-MGD`, 2,698
MAGs"; `pi_login_pulls.md` "`GBR-SED` (`zenodo:24187765`, 2,698 MAGs) embargo until 2027-02-01"; the `tower` audit "`GBR_MGD` (one
5,283-MAG tarball)" with licence field `{}`; (5) the handoff's `REFSEQ-NEW` (complete isolates) sits in the release DEC-70 meant for
diverse environmental waves, where DEC-70 argued isolates yield about 27,000 NR50 per Gbp.

### 5.2 Every current member, scored and placed

State: "pod" = on `wilkcomp02` now; "tower: landed" = directory present and listed as landed in `2026-11-25_tower_sources_incoming.md`
or in the handoff; "tower: queued" = in a job list whose batch has not started (batches 5, 6 and 7 have no directory in
`sources_incoming/`; `batch6b_euks/` and `batch6b_ncbi/` exist, *measured*, `ls` 2026-11-26); "not fetched".

| Member | Planned by | Size | Novelty | Ease | State |
|---|---|---|---|---|---|
| `G-RARE1..3` (GlobDB r226 rare phyla, unclassified, remaining archaea) | handoff, wishlist Z3 | 25.05 to 32.38 Gbp (*derived*, wishlist) | N5 | E1 | pod (`/stor` tarball, verified) |
| `RMD1` Tara Pacific coral | handoff | 4,224 genomes, 11.6 GB (*published*, handoff) | N4 | E3 | tower: queued (batch 6b) |
| `GBR_MGD` | handoff, C-71 | 2,698 or 5,283 MAGs (records disagree) | N4 | E5 (licence `{}`; embargo per login note) | tower: directory present (batch 3, complete 01:55 CST); licence unresolved |
| `MV2` + `MV-EUKV` (MetaVR) | handoff, wishlist Z4 | 3.909 Gbp (*measured*) + 0.8 to 1.8 Gbp (*estimate*) | N4 to N5 | E1 | pod |
| `KAV-HIFI` | handoff, wishlist Z5 | 4.607 Gbp MAGs + 85.39 Mbp NCLDV + 0.274 Mbp mtDNA (*measured*, wishlist) | N4 | E1 | pod (`/tmp`) |
| `INPHARED-2604` | handoff, wishlist Z6 | 458 phages, 39.21 Mbp (*measured*, wishlist) | N2 | E1 | pod and tower |
| `O4` non-RefSeq organelles | handoff, wishlist U3 | 10,763 records, 0.7 to 1.3 Gbp (*estimate*) | N4 | E3/E4 (10,763 NCBI accessions) | tower: queued (batch 6b) |
| `SPX2` SPIRE plant host + briny water + anthropogenic | handoff, wishlist U1 | about 38 Gbp at 5 kb (*estimate*, wishlist) | N5 (plant host 568,000 NR50 per Gbp, *measured* in the pilot) | E4 (about 1,500+ per-sample fetches) | not fetched |
| `LGRT2` Logan multi-bait sweep | handoff, wishlist U2 | 5 to 15 Gbp (*estimate*) | N5 (targeted) | E4 (overnight scan + hundreds of run fetches) | not fetched |
| `LR2` | handoff, survey 1a | 21.549 Gbp new (*measured*) | N4 | E1 | pod (`/stor`, verified) |
| `SAG1` | handoff, survey 1a | 4.924 Gbp new (*measured*) | N5 | E1 | pod (`/stor`, verified) |
| `NX1` | handoff, survey 1a | 3.398 Gbp new (*measured*) | N4 | E1 | pod (`/stor`, verified) |
| `REFSEQ-NEW` | handoff, wishlist Z2 | 5,036 complete reps, 25.812 Gbp (*measured*, wishlist) | N2 | E1 | pod (`/stor`) |
| tower long-read cMAGs (`MFD_LR`, `MIDAS_LR`, `LR_COMPOST`, `LR_ARC`, `LR_REVIO`, batch 6 long-read set) | handoff | `MFD_LR` 20.79 GB landed (55 to 65 Gbp *estimate*); `MIDAS_LR` 3 of 5 files; batch 6 141.82 GB queued | N4 | E2 to E3 | tower: landed in part |
| `GV1` | C-69, survey 1a | 3.89 Gbp new (*measured*) | N4 | E1 | pod and tower (staged) |
| `EGV1` | C-69 | about 51.3 Mbp (*derived*, card) | N5 | E2 | tower: landed (24,753,281 bytes) |
| `ECE2` | C-70 | under 0.2 Gbp (*estimate*) | N4 | E2 | tower: landed (30 files) |
| `MIR1` | C-70 | 0.4 to 0.55 Gbp (*estimate*) | N5 | E2 | tower: landed (10 files) |
| `ASG1` | C-70 | 1.0 to 1.4 Gbp (*estimate*) | N5 | E2 | tower: landed (`Data.zip` 4,961,607,078 bytes; contents unknown) |
| `SAG2` GORG-Dark | survey 2a | 6.466 Gbp, 100.0% novel (*measured*, card C-58) | N5 | E2 | tower: landed (2,471,331,739 bytes) |
| `DPN1`, `SYM2`, `LRP1`, `MIR2`, `MGE4` | survey 1c to 1f | about 2 to 3.5 Gbp together (*estimate*) | N4 to N5 | E2 | tower: landed (`DPN1` 9.99 GB, larger than its stated contents: audit first) |

### 5.3 Before and after (proposed)

| Member | Before (handoff order) | After (proposed) | Reason |
|---|---|---|---|
| `G-RARE1..3` | `r26.11.5`, first | `r26.11.5`, first | on the pod, N5; cut after reading GlobDB r232's taxonomy so that genomes r232 moved into GTDB R232 are noted |
| `MV2` + `MV-EUKV` | `r26.11.5` | `r26.11.5` | on the pod, N4 to N5 |
| `KAV-HIFI` + `INPHARED-2604` | `r26.11.5` | `r26.11.5` | on the pod; small |
| `LR2` + `SAG1` + `NX1` | `r26.11.5`, last | `r26.11.5`, second | staged and verified; the survey's first tier |
| `GV1` + `EGV1` + `ECE2` + `MIR1` + `ASG1` | not listed | **enter** `r26.11.5` | staged or landed; cards C-69 and C-70 put them here; N4 to N5 |
| `SAG2` | not listed | **enter** `r26.11.5` | landed; 100% novel dark-ocean SAGs, measured |
| `DPN1`, `SYM2`, `LRP1`, `MIR2`, `MGE4` | not listed | **enter** `r26.11.5` | landed, compact, N4 to N5 (after the audit of R-11 reads their contents) |
| `GWV1` | not known | **enter** `r26.11.5` | new find; one open file; N5 |
| `SAG3` | not known | **enter** `r26.11.5` | new find; one open file; N4 to N5 |
| `WET1` | not known | **enter** `r26.11.5` | new find; one open file; N4 |
| `G232-A` | not known | **enter** `r26.11.5`, conditional | one anonymous 299 GB stream; exceptional pool of new species; needs a `sources.yaml` entry and the PI's word on the stream |
| `RMD1` | `r26.11.5` | `r26.11.6` (or `.5` if batch 6b lands before the cut) | queued behind batches 5 and 6 (about 313 GB, *derived*: 93.18 + 220.12) |
| `GBR_MGD` | `r26.11.5` | `r26.11.6` or later | licence not stated; the three records disagree on what it is and whether it is embargoed |
| `O4` | `r26.11.5` | `r26.11.6` | 10,763 accessions, queued in batch 6b |
| `SPX2` | `r26.11.5` | `r26.11.6`, first | about 1,500 per-sample fetches; highest measured yield, so it leads the next release |
| `LGRT2` | `r26.11.5` | `r26.11.6` | an overnight scan and hundreds of run fetches; the ART-RT emphasis of DEC-71 was served in `r26.11.4` |
| `REFSEQ-NEW` | `r26.11.5` | `r26.11.6` | complete isolates, N2; re-key against GTDB R232 first |
| tower long-read cMAGs | `r26.11.5` | `r26.11.6` as `LR3` | `MFD_LR` alone is 55 to 65 Gbp (*estimate*); dedupe against GlobDB r232 (mOTUs v4 and RBG draw on the same runs) |

### 5.4 Size of the proposed `r26.11.5`

| Part | Gbp new |
|---|---|
| `G-RARE1..3` | 25.05 to 32.38 (*derived*, wishlist) |
| `MV2` + `MV-EUKV` | 4.7 to 5.7 (*measured* + *estimate*) |
| `KAV-HIFI` + `INPHARED-2604` | 4.73 (*measured*, wishlist; sum of the parts) |
| `LR2` + `SAG1` + `NX1` | 29.87 (*measured*) |
| `GV1` + `EGV1` + `ECE2` + `MIR1` + `ASG1` | 5.4 to 6.0 (*measured* 3.89 + *estimate* 1.5 to 2.1) |
| `SAG2` | 6.47 (*measured*) |
| `DPN1`, `SYM2`, `LRP1`, `MIR2`, `MGE4` | 2 to 3.5 (*estimate*) |
| `GWV1` + `SAG3` + `WET1` | 18 to 19.5 (*estimate*) |
| `G232-A` | up to 12 (cap) |
| **Total** | **about 108 to 120 Gbp** (*derived*: sums of the rows' low and high ends, `G232-A` at its cap) |

Without `G232-A` it is about 96 to 108 Gbp. For comparison, `r26.11.4`'s nineteen stages bring the corpus to 287.03 Gbp (DEC-77); card
A-28 expects `/tmp` scratch limits beyond about 400 Gbp, which this total reaches (*derived*: 287 + 108 to 120 = 395 to 407).

Every entrant is on the pod, landed on `tower`, or one anonymous file (`G232-A`: one stream). Nothing in it needs a login or thousands
of accessions.

### 5.5 Decontamination and duplication catches, against the members they touch

| Catch (source) | Members touched | Action before staging |
|---|---|---|
| `MIR1/gaia2023_GOEV_Genomic_DB.zip`: 660 of 1,817 (36.3%) overlap `GV1`; 54 *Poxviridae* and 13 *Asfarviridae* references (handoff, wishlist Part IV) | `MIR1`, `GV1` (both in the proposed `r26.11.5`) | exact `contig_uid` dedup against `GV1`; drop the *Poxviridae* and *Asfarviridae* isolate references |
| `GFS1`: 14 eukaryotic assemblies over 20 Mb (handoff); and this note's catch: GlobDB r226 `GFSMAG` holds the same project's 1,377 non-GTDB species | `GFS1` (not in `r26.11.5`) | eukaryote filter and species dedup against `GFSMAG` before any stage |
| `WWV1`: RNA viruses, prokaryotic MAGs, proteins and DNA vOTUs mixed (handoff) | `WWV1` (not in `r26.11.5`) | split by realm; proteins are out of scope |
| `uhgv_extended.fna.gz`: 195,719 of 212,778 (92.0%) already in `V3` (handoff) | none in `r26.11.5` | only the 17,059 non-`V3` genomes (0.606 Gbp) if ever staged |
| `jgi_img_pr`: 0 plasmids beyond `PL3` (handoff) | none in `r26.11.5` | sidecar metadata only |
| `CFMD1`, `CRBC_ISO`, `QXLSG1` and eight smaller catalogues are inside GlobDB r232 (§2.7) | none in `r26.11.5` | take from the r232 stream or dedupe against its table |
| `ASG1`: GTDB R232's 3,154 new archaeal species and GlobDB r232's `TPMCS` may hold some of Appler 2026's MAGs | `ASG1` | `skani` 96% against GTDB R232 and GlobDB r232 archaea, not only `G1` |
| `DPN1`: 9,991,035,723 bytes landed for "256 DPANN + 3,662 DPANN/MGEs" (*measured*, tower note) | `DPN1` | stage only genome FASTA after the R-11 audit reads the archive list |
| `GWV1`: reads `PRJNA858913` are public, so MetaVR may contain some vOTUs | `GWV1` | exact `contig_uid` and 95% ANI check against `MV1`, `V2`, `V3` |

### 5.6 What the proposal would amend (for the PI to ratify)

- The handoff's H2 step 3 order (`STATUS.md`, 2026-11-26 07:15), which is recorded there as PI-approved; no DECISIONS row holds it.
- Card C-69 (no move: `EGV1` stays in `r26.11.5`; its fetch is done), card C-70 (`ECE2`, `MIR1`, `ASG1` fixed to `r26.11.5` rather
  than "`r26.11.5–r26.11.6`"; `ASG1` dedup widened to r232), card C-71 (`GBR-MGD` held on licence; `ARC1`, `NX3` dedup against r232).
- `CATALOG_GOAL.md` §2 items 3 and 4 (order), when the PI decides; this card adds only the §4.3 rows (§6 below).
- A new DECISIONS row is needed to record the PI's choice (this card does not write it).

## 6. Changes made by this card outside this note

- `planning/tasks/TASKS.md`: card C-72 after C-71; dated lines in C-69, C-70 and C-71.
- `planning/CATALOG_GOAL.md` §4.3: rows for GlobDB r232 and GTDB R232 (the GTDB and GlobDB remainder rows of §4.1 speak of R226 and
  r226 only, and "every other GTDB genome" is no longer the whole of GTDB), and for the six new open sources.

## 7. Sources that need the PI's login, an access request, or a browser

None of this note's new finds needs a login. The list of `2026-11-25_pi_login_pulls.md` stands unchanged (JGI SSO; Biosino NODE;
Dryad ORCID for `PEONY`/`PIGEON`; CNCB/NGDC for `GCSV`; the Zenodo access request for `20336193`; CyVerse and firewall-blocked
sites). Two additions:
- `nature.com` article pages (data-availability sections of `UniMAG v1.0`, the *Nature Water* sludge catalogue) answer with a
  cookie-setting redirect; a browser read by the PI would settle both records.
- `GBR_MGD`: the PI's reading of the record's licence (field `{}` in the job list).

## 8. What could not be verified

- GlobDB r232's per-resource base counts, MIMAG tiers and per-genome source table (none published; the dictionaries archive maps IDs
  for some resources only); the 118 Gbp and 1.04 Tbp figures are estimates.
- Whether GTDB R232 includes Appler et al. 2026's Asgard MAGs (its METHODS state no NCBI cutoff date).
- MGnify licences (not stated in the READMEs read); the byte sizes of MGnify `marine` v3.0's and `marine-sediment` v1.0's species
  directories (not listed without opening thousands of directories).
- `UniMAG v1.0` (only a search abstract); the *Nature Water* activated-sludge catalogue's repository; the Lake Zurich 72,000-MAG
  series (not released); Maeke et al. 2025 coastal dataset (record not opened); Buck et al. 2021 lakes' presence in GTDB R232 or
  mOTUs v4.
- The OMZ compendium's overlap with GTDB R232 (Table S1 not opened) and the GWVC's overlap with MetaVR.
- The state of `tower` beyond the two `ls` calls: whether `GBR_MGD`, `LR_ARC`, `LR_REVIO` and batch 4's sets landed whole is for
  card R-11.

## 9. Records read (all 2026-11-26)

- `https://globdb.org/`, `https://globdb.org/downloads`, `https://globdb.org/news`;
  `https://fileshare.lisc.univie.ac.at/globdb/globdb_r232/` (index) and `.../globdb_r232_dataset_list.tsv`.
- `https://forum.gtdb.ecogenomic.org/t/announcing-gtdb-r11-rs232/826`; `https://gtdb.ecogenomic.org/stats/r232`;
  `https://data.gtdb.aau.ecogenomic.org/releases/release232/232.0/` (index, `RELEASE_NOTES.txt`, `METHODS.txt`,
  `genomic_files_reps/` index).
- `https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/` and the indexes and READMEs of `marine/v3.0/`,
  `marine-sediment/v1.0/`, `marine-eukaryotes/v1.0/`, `human-gut/v3.0/`, `barley-rhizosphere/v2.0/`.
- `https://zenodo.org/api/records/11230969` (GWVC); `https://zenodo.org/api/records/20263536` (GWM);
  `https://api.figshare.com/v2/articles/20459526` (OMZ SAGs).
- Europe PMC REST: `PMC11310336`, `PMC13547364`, `PMC10224968` full text; search for the *Nature Water* DOI (no hit).
- `https://data.microbiomedata.org/api/stats`; `https://meetingorganizer.copernicus.org/WBF2026/WBF2026-597.html`; the NTU news
  page on the *Nature Water* paper.
- `tower`: `ls /mnt/user/research-bulk/bio-box/sources_incoming` and `ls /mnt/user/research-bulk/bio-box/sources`.
