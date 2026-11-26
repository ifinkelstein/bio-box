# Card C-72, part 1: the five sources of DEC-79 (4), listed for the fetch to `tower`

Written 2026-11-26 (records read about 08:12 to 08:15 CST, 14:12 to 14:15 GMT). Nothing was fetched. Every size and md5 below is *measured*: read from the provider's own record on 2026-11-26 from the listing named for the source. No request carried a name, an email address or a project detail (curl `-A curl`, WebFetch). Fifteen files, 6,247,225,972 bytes in all (*derived*: the sum of the fifteen sizes below; 5.82 GiB). The job list is `scripts/orchestrator/tower_jobs/2026-11-26_c72_part1.tsv`; the entries are `tower_gwv1`, `tower_sag3`, `tower_wet1`, `tower_globdb_r232_tables`, `tower_gtdb_r232_meta` in `sources.yaml`. Every URL in them is one exact file (no `{file}` form), so no other file or record of the same host is covered. The entries hold the provider's sizes and md5 under `files_listed`; `files` (with our sha256) is written after the fetch, as the not-yet-fetched entries before them wrote nothing for it.

## 1. GWV1: Groundwater Virome Catalogue (Zenodo 11230969)

Record `https://zenodo.org/api/records/11230969` (published 2024-05-21), read 2026-11-26. Licence as the record states it: `cc-by-4.0` (an id only, no title).

| file | bytes | md5 |
|---|---|---|
| GWVC.vOTUs.virus-host.linkage.txt | 59,543,590 | 3e346e34ffedde104e16cde4b398eae6 |
| GWVC.vOTUs.information.txt | 62,086,246 | 4878fa7168ad1078752a7dc01dc1487c |
| GWVC.vOTUs.fna.gz | 1,166,850,604 | d249dd5821d411e195178bf1749c525a |

Not taken: `CPR.DPANN.genomes.tar.gz` (43,339,586 bytes, md5 3a1e3cb0ed70897fada95d5c499ccbc1). The record's description does not say what it holds: it states only that the data come from ultra-deep metagenomic sequencing of groundwater from 607 monitoring wells in China and that "Descriptions of files are provided in the file 'README.txt'". The brief took the file only if the description showed nucleotide genomes; it does not, and reading `0README.txt` (954 bytes) is a read of a data file the brief does not allow. If the orchestrator wants it, one line in the job list and one URL line in `tower_gwv1` add it (43.3 MB). Also not taken: `GWVC.vOTUs.spacer.hit.table.out`, `Caudoviricetes.77markers.tar.gz`, `Groundwater.CPR.viruses.spacer.hit.table.out`, `Spacer.squences.fna`, `Groundwater.CPR.viruses.linked.to.non-CPR.host.spacer.hit.table.out`, `Groundwater.DPANN.viruses.spacer.hit.table.out`, `0README.txt` (not in the brief's list).

## 2. SAG3: OMZ single-cell compendium (figshare 20459526 and 20481603)

Articles `https://api.figshare.com/v2/articles/20459526` (published 2023-04-06, "File 5. Single Cell Assemblies from the OMZ Compendium") and `.../20481603` (published 2023-05-23, "File 1. OMZ Compendium Metadata"), read 2026-11-26. Licence as both articles state it: CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/). Each article holds one file; none is left out.

| file | bytes | md5 | URL |
|---|---|---|---|
| Table_S1_Metadata-template-Bio-Med-SAGdescriptor-OMZ_May_23_2023.xlsx | 4,433,905 | 004c2b180c199f4f39c463d888b33d6e | https://ndownloader.figshare.com/files/40838783 |
| OMZ_Compendium_SAGs_Full.tar.gz | 747,416,178 | 49e76ad707b718cd4e804e43b6d0cb94 | https://ndownloader.figshare.com/files/37394083 |

## 3. WET1: GWM wetland catalogue (Zenodo 20263536)

Record `https://zenodo.org/api/records/20263536` (published 2026-05), read 2026-11-26. Licence as the record states it: `cc-by-4.0` (an id only). The record's md5 values carry the prefix `md5:`; the list gives the bare digest.

| file | bytes | md5 |
|---|---|---|
| MAG_meta.xlsx | 806,764 | 54b567ff549a0e926f54c798609d1657 |
| Dereplicated genomes.zip | 3,885,785,718 | c2c471cb64ee73a223bfcef33277fe6d |

The second file's name has a space; the job list's file field holds it, and the URL is `.../files/Dereplicated%20genomes.zip/content`, exactly as the record gives it. Not taken (as the brief says): `Genome sequences.zip` (5,289,781,056), `Gene sequences.zip` (4,043,066,869), `Functional annotations.zip` (1,549,644,406), `Taxonomic profiles.zip` (37,564,189), `MAG_abundance.tsv` (17,011,380), `Source_Info.xlsx` (66,878), `Sample_meta.xlsx` (280,621).

## 4. GLOBDB_R232_TABLES: GlobDB r232 (https://fileshare.lisc.univie.ac.at/globdb/globdb_r232/)

Directory index and `globdb_r232_md5sum.txt` read 2026-11-26. The index prints sizes in a column labelled KB that is not in KB (it prints `8071 KB` for a file of 8,071,276 bytes), so the bytes come from `Content-Length` of a HEAD request on each file (*measured*, same day). Licence: CC BY-SA 4.0 as `planning/measurements/2026-11-26_missed_diversity_sources.md` section 2.1 quotes globdb.org (read 2026-11-26 by the survey of card C-72); this card did not read globdb.org again.

| file | bytes | md5 (from the md5 list) |
|---|---|---|
| globdb_r232_md5sum.txt | 519 | none (the list does not list itself) |
| globdb_r232_dataset_list.tsv | 4,358 | none (not in the list) |
| globdb_r232_dictionaries.tar.gz | 3,266,383 | 9a493f06c4b41a7b6b63eff6a650a77d |
| globdb_r232_taxonomy.tsv.gz | 8,071,276 | 7dd4a6b2c3bc3c7727118b3f772f7a4e |
| globdb_r232_checkm2.tsv.gz | 12,896,811 | 61da09cb728e7fe4c5116de54e2f173f |

Not taken: `globdb_r232_genome_fasta.tar.gz`, `globdb_r232_protein_fasta.tar.gz`, `globdb_r232_protein_annotations.tar.gz`, `globdb_r232_gff_cog.tar.gz`, `globdb_r232_anvio_dbs.tar.gz` (the brief's exclusions) and the folder `taxonomic_profiling/` (not named by the brief).

## 5. GTDB_R232_META: GTDB release 232.0 (https://data.gtdb.aau.ecogenomic.org/releases/release232/232.0/)

Directory index and the release's md5 list read 2026-11-26. The md5 list is `MD5SUM.txt` (the brief wrote `MD5SUM`; no file of that name exists, a request for it returned HTTP 404). Bytes from `Content-Length` of a HEAD request on each file (*measured*). Licence as entry `gtdb_metadata` records it: CC BY-SA 4.0 ("GTDB's stated licence; not re-verified"); the release directory states none, and this card did not re-verify it. The host is `data.gtdb.aau.ecogenomic.org` as the brief gives it (the r226 entry names `data.gtdb.ecogenomic.org`).

| file | bytes | md5 (from MD5SUM.txt) |
|---|---|---|
| MD5SUM.txt | 3,950 | none: it lists itself (b289d22dd671b2d69117ae3104fbf40c), a value that cannot be checked against a file that holds it |
| ar53_metadata_r232.tsv.gz | 7,015,869 | 0f51fe6f01b4ccbf0d9bab6c07df87af |
| bac120_metadata_r232.tsv.gz | 289,043,801 | 0ba4237077b65cfc5556e1d8be797485 |

## Checks

`sh scripts/orchestrator/fetch_to_tower.sh --check --sources sources.yaml scripts/orchestrator/tower_jobs/2026-11-26_c72_part1.tsv`: each of the fifteen lines printed `match <its own entry> <url>` (for example `match tower_wet1 https://zenodo.org/api/records/20263536/files/Dereplicated%20genomes.zip/content`); the last line:

```
check: 15 URLs, all covered by sources.yaml; nothing fetched
```

Refusal, with five invented URLs of the same hosts (another Zenodo record id, a Zenodo file of the right record that the list does not take, a figshare id one above a listed one, a GlobDB archive the brief excludes, a GTDB release 233 path), exit 2:

```
fetch_to_tower: refuse.tsv line 1: URL not covered by any sources.yaml entry: https://zenodo.org/api/records/11230970/files/GWVC.vOTUs.fna.gz/content
fetch_to_tower: refuse.tsv line 2: URL not covered by any sources.yaml entry: https://zenodo.org/api/records/20263536/files/Genome%20sequences.zip/content
fetch_to_tower: refuse.tsv line 3: URL not covered by any sources.yaml entry: https://ndownloader.figshare.com/files/37394084
fetch_to_tower: refuse.tsv line 4: URL not covered by any sources.yaml entry: https://fileshare.lisc.univie.ac.at/globdb/globdb_r232/globdb_r232_genome_fasta.tar.gz
fetch_to_tower: refuse.tsv line 5: URL not covered by any sources.yaml entry: https://data.gtdb.aau.ecogenomic.org/releases/release233/233.0/bac120_metadata_r233.tsv.gz
fetch_to_tower: 5 refusal(s); nothing was fetched
```

## Conflict with an existing test

`tests/r1/test_r08_fetch_to_tower.py::test_the_tree_sources_yaml_parses_and_every_tower_entry_has_its_record` asserts that every `tower_` entry's `approved` starts with `by the PI, 2026-11-25 ` and that every item of `files` has a sha256. The five entries were approved on 2026-11-26, so the first assertion fails for `tower_gwv1`; the second is avoided by keeping the sizes and md5 under `files_listed` until the fetch gives a sha256. The orchestrating session then named the test in this card's files for that one assertion: it now requires `by the PI, ` and a date `YYYY-MM-DD` and a space. Nothing else in the test changed.
