# One assembly per sequencing run: SPIRE, MGnify, INSDC (NCBI WGS) and what IMG/M would add

Card C-43, 2026-11-16. Tables only: no sequence was fetched. Every table was fetched on wilkcomp02 into
`/tmp/if2453/s3/tables/` in the runtime image `bio-box_faec1f4.sif` (one connection at a time, nice 19, 2 GB
scopes), and the crosswalk was built there by `scripts/stress/c43_crosswalk.py build` (16 GB scope, one thread).
Each URL, file, size, sha256 and fetch time (the pod's clock, UTC) is in `sources.yaml`, under `spire_tables`,
`ncbi_metagenome_assembly_summary`, `spire_study_api`, `mgnify_api` and `ena_portal_tables`. The last three
entries were written under DEC-50 before any request to them.

Labels. *measured*: counted by us from the fetched tables (wilkcomp02, 2026-11-16, the `build` command in §7,
and its output file named next to the number). *published*: the provider's own figure, with where it is stated.
*derived*: arithmetic on the numbers named. *estimate*: a judgement. "Read bases" means the `base_count` of the
run in ENA's run table: what was sequenced, not what was assembled. Assembled bases are given only where a
table states them.

## 1. Resources and how a unit was linked to its runs

A unit is what a resource assembled: a SPIRE sample, an MGnify assembly (ERZ), or an INSDC assembly (GCA or
WGS set). Units were linked to runs only through their tables' own cross-reference columns, in this order:
(1) the runs the unit's own table names (MGnify's `runs`, ENA's WGS `run_accession`, a SPIRE sample id that is
itself a run accession, 580 of them); (2) its sample accession (primary BioSample or secondary ERS/SRS/DRS)
looked up in ENA's run table (`read_run`, `library_source="METAGENOMIC"` without amplicons, 2,140,720 runs),
then, for samples still unmatched, in ENA's runs of those samples by any library source (75 lookups of 100
samples; amplicon runs dropped); (3) otherwise the unit is unlinked. Nothing was matched by name.

| Resource | Units | By its own run list | By sample | Unlinked | Distinct runs | Read bases of those runs | Assembled bases stated |
|---|---|---|---|---|---|---|---|
| SPIRE (study tables of `spire.embl.de`) | 96,067 samples | 580 | 85,391 | 10,096 | 122,062 | 613.9 Tbp | none in any table |
| MGnify (API listing plus ENA's analyses by MGnify) | 83,313 assemblies | 80,865 | 2,242 | 206 | 85,962 | 336.3 Tbp | 2.67 TB of compressed contig files for 79,102 |
| INSDC primary metagenome assemblies (ENA's WGS table; NCBI's assembly table marks its subset) | 6,463 | 1 | 3,702 | 2,760 | 4,990 | 45.7 Tbp | 1.457 Tbp for 6,459 |

All *measured* (`resources.tsv`). Read bases are unknown for 9,976 of MGnify's runs: runs that MGnify names but
that ENA's metagenomic run table does not hold.

Notes on the units:

- **SPIRE.** Its studies page lists 715 studies and 96,067 samples (*measured*, `studies.html`), against
  99,146 samples of 739 studies in SPIRE's paper (*published*, Schmidt et al. 2024, doi:10.1093/nar/gkad943;
  3.1% fewer samples). No SPIRE table gives run accessions or assembled bases: the study tables hold
  `sample_id, mags, lat, lon, microntology`, and the README covers only the MAG, cluster and microntology
  tables. The 16 Tbp of contigs is SPIRE's own total (*published*, survey §2 #10). The 10,096 unlinked
  samples carry a study's own names, not INSDC accessions: MetaSUB 4,116, HMP 1,146, MetaHIT 816, Vatanen
  2016 812, Kenny 2020 624, and others (*measured*, `units.tsv`).
- **MGnify.** The API lists 51,739 assemblies (*measured*, 207 pages; the API's own count matched). ENA lists
  79,102 primary metagenome assemblies submitted by MGnify (`center_name="EMG"`; *measured*,
  `analysis_emg_assemblies.tsv`), and 31,574 of them are not in the API listing. Both are counted here as
  MGnify's (*derived*: 83,313 − 51,739). ENA indexes MGnify's assemblies under no WGS set and no assembly record:
  261 matched a WGS set, and 2 are in ENA's assembly table. So their bases are not stated anywhere we read. What
  is stated is the size of their gzipped contig files, 2,672,559,908,917 bytes for the 79,102 (*measured*).
  Habitat comes from MGnify's sample listing (biome lineage). For 18,230 assemblies whose sample is not in that
  listing, it comes from the INSDC organism name of the run (§3).
- **INSDC.** ENA's table of WGS sets under the NCBI taxon "metagenomes" has 108,493 sets (*measured*). It is
  mostly bins: 102,261 are typed MAG (262.9 Gbp), 79 are clone or isolate, and 6,153 carry no type (1.387 Tbp).
  The untyped sets are taken as primary assemblies, and 261 of them are MGnify's own (counted under MGnify).
  NCBI's `genbank/metagenomes/assembly_summary.txt` has 27,219 assemblies (*measured*), and 20,696 of them are
  ENA-typed MAGs. So NCBI's "metagenomes" group is not a list of primary assemblies. Of the INSDC units, 1,324
  are under 10 Mbp but hold only 4.5 Gbp. They are likely untyped bins (*estimate*) and are kept and flagged
  (§4).

## 2. Runs by the resources that hold them, and habitat

Each run is counted once, under the combination of resources that assembled it and under the habitat of the
first of them in the order SPIRE, MGnify, INSDC. *Measured*, `runs_by_combination_and_habitat.tsv`; read bases
in Tbp.

| Held by | host-associated | aquatic | terrestrial | anthropogenic | air | unlabelled | All runs | Read bases |
|---|---|---|---|---|---|---|---|---|
| **SPIRE only** | 67,298 | 8,525 | 8,213 | 5,791 | 17 | 3,253 | **93,097** | 478.0 |
| SPIRE + MGnify | 21,820 | 3,833 | 329 | 449 | 0 | 802 | 27,233 | 126.2 |
| SPIRE + INSDC | 12 | 87 | 125 | 0 | 0 | 0 | 224 | 3.7 |
| SPIRE + MGnify + INSDC | 1,253 | 255 | 0 | 0 | 0 | 0 | 1,508 | 6.0 |
| **MGnify only** | 38,141 | 3,485 | 790 | 1,042 | 11 | 13,429 | **56,898** | 200.9 (9,976 runs without bases) |
| MGnify + INSDC | 76 | 221 | 9 | 17 | 0 | 0 | 323 | 3.2 |
| **INSDC only** | 782 | 902 | 663 | 436 | 0 | 152 | **2,935** | 32.8 |
| All | | | | | | | 182,218 | 850.8 |

The unlinked units are beside this table, not in it, because they have no run: SPIRE 10,096 samples, MGnify
206 assemblies, INSDC 2,760 assemblies (§1). Read as (*derived* from the table): MGnify alone holds 56,898 runs
that SPIRE lacks (201 Tbp of reads), 61% as many runs as SPIRE alone (93,097). INSDC alone adds 2,935 runs.
Only 28,741 runs (15.8%) were assembled by both SPIRE and MGnify.

## 3. Habitat labels

- SPIRE: its own microntology tags per sample, from the study tables. SPIRE tags a sample with several terms
  and pulls `aquatic` into sewage and mangroves, so one label is chosen by a stated order. The first top level
  present wins: host-associated, anthropogenic, terrestrial, aquatic, air. Other categories (age group,
  salinity, temperature and so on) describe a sample and are not habitats. SPIRE's second level is in
  `spire_habitat_level2.tsv` for card C-42; for example host-associated:animal host has 62,392 samples and
  aquatic:marine 5,340 (*measured*).
- MGnify and INSDC: mapped to those five labels only through the table `2026-11-16_run_crosswalk/habitat_map.tsv`
  (5 MGnify lineage prefixes, 128 INSDC organism names), which is part of this note. Anything else is
  `unlabelled`. MGnify's biome lineage is used where its sample listing gives one; otherwise the INSDC organism
  name of the assembly's run (e.g. "human gut metagenome") goes through the same organism rows, and the unit's
  detail says `insdc:<name>`. Labels left unlabelled, with counts, are in `unlabelled_sources.tsv`: the largest
  are MGnify units with no biome and no mapped organism, `root:Mixed` (6,652 MGnify units), and the generic
  INSDC names "metagenome" (685) and "viral metagenome" (208).
- The map is my reading of names (*estimate*): for example "mine" and "tailings" go to anthropogenic, as in
  SPIRE's own tree, and the hosts of organismal metagenomes (algae, moss, lichen, fungus) go to
  host-associated.

## 4. Order of preference, and what to take from each resource after SPIRE

The order: **SPIRE first** (one pipeline, MEGAHIT throughout, with curated habitat labels), **then MGnify** (one
institution's pipelines; biome labels on most samples; its own assemblies of public runs), **then the other INSDC
primary assemblies** (many submitters and pipelines, habitat only from an organism name). Within a resource,
units are taken in accession order. A unit is **taken** when none of its runs is taken yet and none of its
samples was assembled by a taken unit. It is **covered** when all its runs are taken, or its sample was. It is
**partial** when some of its runs are taken and some are not: it is not taken, because it would repeat reads.
It is **take_unlinked** when it has no run to check and its sample is new. The tests in
`tests/r1/test_c43_crosswalk.py` state these rules on constructed tables.

| Resource | Taken | Runs | Read bases | Assembled size stated | Take, unlinked | Covered | Partial |
|---|---|---|---|---|---|---|---|
| SPIRE | 85,971 samples | 122,062 | 613.9 Tbp | (16 Tbp in all, *published*) | 10,096 samples | 0 | 0 |
| MGnify | 37,628 assemblies | 44,177 | 179.5 Tbp | 1.345 TB gzipped contigs | 200 | 45,299 | 186 (32,625 runs) |
| INSDC | 2,306 assemblies | 2,922 | 32.5 Tbp | 550.3 Gbp | 2,642 (739.2 Gbp) | 1,480 | 35 |

*Measured*, `decisions_by_habitat.tsv` (also by habitat). Of MGnify's 45,299 covered assemblies, 30,199 share a
run with SPIRE, 452 share a sample with SPIRE, and 14,648 are a second MGnify assembly of a run already taken.
MGnify reassembled the same runs (*measured*, `units.tsv`). The 186 partial MGnify assemblies are co-assemblies
over 32,625 runs that overlap a taken unit. Of the INSDC assemblies taken or taken unlinked, 1,324 are under
10 Mbp (4.5 Gbp in all). The other 3,624 hold 1.285 Tbp (*measured*). 165 INSDC units (208 Gbp) were
submitted by the DOE Joint Genome Institute (*measured*, NCBI `asm_submitter`).

The run lists for each resource stay on the pod (large):

| File (in `/tmp/if2453/s3/c43_work/final2/`) | Bytes | sha256 |
|---|---|---|
| `take_spire.tsv` (uid, decision, runs, habitat) | 7,207,435 | `f0b90a90e9555546e50d08ccf8c6f1912b517d5784f94464c3a3a1380a1e179e` |
| `take_mgnify.tsv` | 3,353,304 | `0189667f131430493dfb8c24c153d462e29da4fdb4f43439071b72e1773f7898` |
| `take_insdc.tsv` | 360,908 | `f2e77bf172ec35ad388ad281210e73091351b11c43f6f0201e7e9bc4b90010b0` |
| `units.tsv` (every unit, its runs, samples, decision, bases, habitat) | 27,823,444 | `1e39e4fb4bf40e39f1c655402a5373e467fb34fc4fdba3197360b86ae900f4fb` |

The build was run twice into separate directories, and the outputs are identical (`diff -r`, *measured*).

**What this means for card C-42 and later (a recommendation).** Choose C-42's pilot from `take_spire.tsv`. All
SPIRE samples are kept, including those not linked to a run, since SPIRE is first. After SPIRE, MGnify is the
large addition: 37,628 assemblies, 1.35 TB compressed. The largest share is host-associated: 28,398 assemblies
and 130 Tbp of reads (*measured*). INSDC adds a little on top: 2,306 linked assemblies (550 Gbp) and 2,642
unlinked ones (739 Gbp). The unlinked ones are taken only by sample and could repeat reads that another
resource assembled under a sample accession ENA does not connect.

## 5. IMG/M and GOLD: what the crosswalk lacks, and what the PI would need to fetch

No JGI host was tried (DEC-50 (2)). Without IMG/M and GOLD:

- **IMG/M's own assemblies are missing.** Only 165 JGI-submitted metagenome assemblies are in INSDC's primary
  set (§4). IMG/M holds tens of thousands of metagenomes. The survey quotes 75.1 billion metagenome genes
  (*published*, IMG 2022), and Faure 2023 searched 6.93 Tbp of NCBI plus JGI contigs (*published*, survey §1),
  against 1.46 Tbp of INSDC primary assemblies here. So most of JGI's assembled bases are probably outside this
  crosswalk (*estimate*). Its runs and habitats are unknown.
- **Run accessions for JGI-sequenced samples.** Many JGI metagenomes have SRA runs (JGI deposits reads), but
  the link from an IMG taxon OID to an SRA run is in GOLD, not in any table read here. Without it, an IMG/M
  assembly cannot be matched to SPIRE's or MGnify's assembly of the same reads. SPIRE's studies named `JGI_*`
  (for example `JGI_soil_forest_LTER`, 359 unlinked samples; `JGI_darkmatter`, 114) cannot be linked either.

Tables the PI would fetch with the JGI login (no sequence):

1. **GOLD's public data export** (`goldData.xlsx` from GOLD's Downloads page): the sheets *Study*,
   *Biosample*, *Sequencing Project* and *Analysis Project*, for the analysis project → sequencing project →
   SRA/BioSample accessions and IMG taxon OID links, and GOLD's five-level ecosystem classification.
2. **IMG/M's list of metagenomes** (IMG Genome Browser export of all public metagenome assemblies): taxon OID,
   GOLD analysis project ID, ecosystem path, assembled bases (or "Genome Size"), scaffold count, sequencing
   center, and the JGI data-use status (for the `redistributable` flag).

With those two, `c43_crosswalk.py` needs one more reader (IMG/M units linked through GOLD's SRA accessions,
habitat from GOLD's ecosystem through a stated map). IMG/M would go fourth in the order, after INSDC, because
its terms limit redistribution (DEC-51).

## 6. Checks and differences from the plan

- The card estimated the tables at "a few GB". What was fetched is 0.58 GB (*derived*: the sum of the sizes
  in `sources.yaml`). It is smaller, and nothing depends on the figure.
- SPIRE's listed samples are 3.1% fewer than its paper's (§1).
- MGnify's proteins and the survey's other MGnify figures were not checked: not part of this card.
- Two requests in the `includeAccessions` form were answered HTTP 400. The entry was corrected to ENA's query
  form before it was used (`sources.yaml`). No other refusal (403, 429, 5xx) occurred.
- ENA's assembly table (`result=assembly`, 50,963 rows) was fetched and not used. It holds 2 of MGnify's
  assemblies and fewer primary assemblies than the WGS table.

## 7. Commands

On wilkcomp02, in `/tmp/if2453/s3/` (scripts in `c43_scripts/`, tables in `tables/<source>/`, outputs in
`c43_work/`). Every command ran in a `systemd-run --user --scope -p MemoryMax=... nice -n 19 singularity exec
--cleanenv --bind /tmp/if2453 /tmp/if2453/images/bio-box_faec1f4.sif python -I`, and the paged fetches ran in
transient `--collect` units named `c43-*`:

```
c43_fetch.py URL OUT [--md5 HEX]                       # each single table, md5 checked where SPIRE gives one
c43_crosswalk.py spire-urls tables/spire/studies.html c43_work/spire_study_urls.tsv
c43_pages.py list c43_work/spire_study_urls.tsv tables/spire/study_tsv --interval 1.0
c43_pages.py next 'https://www.ebi.ac.uk/metagenomics/api/v1/assemblies?format=json&page_size=250&page=1' tables/mgnify/assemblies --interval 0.5
c43_pages.py next 'https://www.ebi.ac.uk/metagenomics/api/v1/samples?format=json&page_size=250&experiment_type=assembly&page=1' tables/mgnify/samples_assembly --interval 0.5
c43_crosswalk.py lookup-urls c43_work/pass1/units.tsv c43_work/lookup_urls2.tsv
c43_pages.py list c43_work/lookup_urls2.tsv tables/ena/read_run_by_sample --interval 1.0
c43_crosswalk.py build --spire-dir tables/spire/study_tsv \
  --read-run tables/ena/read_run_metagenomic.tsv tables/ena/read_run_by_sample/batch_*.tsv \
  --insdc-assemblies tables/ena/wgs_set_metagenomes_typed.tsv \
  --ncbi-summary tables/ncbi_metagenomes/assembly_summary.txt \
  --mgnify-assemblies tables/mgnify/assemblies --mgnify-samples tables/mgnify/samples_assembly \
  --mgnify-ena-analyses tables/ena/analysis_emg_assemblies.tsv \
  --habitat-map c43_scripts/habitat_map.tsv --insdc-include-types 'primary metagenome,' --out-dir c43_work/final2
```

The small outputs are copied beside this note in `2026-11-16_run_crosswalk/`: `resources.tsv`,
`runs_by_combination_and_habitat.tsv`, `decisions_by_habitat.tsv`, `insdc_not_units.tsv`,
`spire_habitat_level2.tsv`, `unlabelled_sources.tsv` and `habitat_map.tsv`.
