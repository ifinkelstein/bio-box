# Census and Staging of Pod-Resident JGI Uncultivated Prokaryote Archives (`IUIG` + `IMG-UMAG`) as `JGIUG1` (Card C-87, DEC-80, DEC-88)

Date: 2026-11-27, `wilkcomp02`, 16:09 to 22:43 CST (`tower` mirror completed 2026-11-28 00:01:20 CST: `12` metadata files `471,197,303` B + `jgiug1_genomes.fna.gz` `11,290,730,624` B in `582` s at `19.40` MB/s, `--streams 1`, `13` lines verified in `SHA256SUMS`, *measured*). Every number carries a label (`AGENTS.md` §6): *measured* (host, date, command or report), *derived* (formula and inputs), *estimate*, or *published*. Unless stated otherwise, every *measured* number is on `wilkcomp02`, 2026-11-27, from `scripts/stress/build_stage_pod/c87_jgi_ug_inputs.py` (`census`, `reselect`, `stage`, `reconcile-codes`) inside `/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif` (`census_report.json` and `stage_report.json` under `/tmp/if2453/c87/logs/`). Nothing was built. All pod-resident archives under `/stor/work/Wilke/if2453/all-the-data/jgi_20261118/`, the data root `/tmp/if2453/r10/root`, `/tmp/if2453/s2/inputs`, `/tmp/if2453/c74`, `/tmp/if2453/c75`, `/tmp/if2453/s3`, and `/tmp/if2453/c83` were bound `:ro` throughout (`AGENTS.md` §7.5); zero `tower` I/O occurred during census and staging (`DEC-88 (1)`).

Script: `scripts/stress/build_stage_pod/c87_jgi_ug_inputs.py`; tests: `tests/r1/test_c87_jgi_ug_inputs.py`; template: `scripts/stress/build_stage_pod/stages/templates/JGIUG1.conf`.

## Summary

- **Full stream-census of both archives (`933.72 GB` compressed):** **`318,606` prokaryote assemblies/bins (`992,807,095,242` bp = `992.81 Gbp`)** (*derived*: `126,346` `IUIG` prokaryote assemblies [`481.61 Gbp`] + `192,260` `IMG-UMAG` prokaryote bins [`511.19 Gbp`], both *measured*). `IUIG` additionally holds `17,450` non-prokaryote bundles (`15,544` `Viruses`, `1,185` `Plasmid:*`, `721` `Eukaryota`, *measured*; logged in `iuig_non_prok_skipped.tsv` and left for viral/plasmid/eukaryote cards).
- **`r232` cross-join (`C-74` index + GEM OTU crosswalk):** **`136,000` prokaryote assemblies/bins (`365.49 Gbp`)** have `r232_status = not_in_r232` (`8,053` in `IUIG` [`6.37%`], `127,947` in `IMG-UMAG` [`66.55%`], *measured* & *derived*).
- **`JGIUG1` (`<= 35.00 Gbp` highest-novelty tranche):** **`21,524` distinct `genome_uid`s** (`21,646` streamed non-`r232` `HQ`/`MQ` genomes [`1,397` `IUIG` + `20,249` `uMAG`; `6,685` `Archaea` + `14,961` `Bacteria`] less `122` exact `DEC-27` contig-set repeats), **`2,890,974` distinct contigs (`2,938,090` FASTA records), `34,657,341,931` distinct bp (`34.657 Gbp`)** (*measured*). **`34,447,066,468` bp (`34.447 Gbp`, `2,875,624` contigs, `99.39%`) are net new** after all `56` held and staged stages `S1..GBRMGD` (*measured* & *derived*, §5).
- **Remaining reserve (`JGIUG_left_for_later.tsv.gz`):** **`295,882` `HQ`/`MQ` candidates (`957.77 Gbp`)**, including **`114,005` strictly novel `not_in_r232` spillover genomes (`330,478,507,298` bp = `330.48 Gbp`)** in deterministic rank order so `JGIUG2+` (`~9`–`10` tranches at `35 Gbp`) is 1-command ready (*measured* & *derived*).
- **Genetic codes & validation:** `21,205` genomes code `11`, `319` code `25` (CPR/`Patescibacteria` UGA->Gly readthrough), `0` code `4`, `1` shared-contig component conflict reconciled via Union-Find (*measured*). `s2_check_inputs.py` full byte pass exited `0` with `problems: []` (`0` non-IUPAC or gap bytes across all `35,529,418,406` record bases); `stage_est.py conf` exited `0` (*measured*).
- **Pod compute** (*measured*): `c87-census` (`CPUQuota=250%`, `MemoryMax=32G`, `nice 19`) took `5,635.6` s (`93.9` min) wall, `80.973` CPU-min (`IUIG` worker, `402 MB` RSS) + `42.439` CPU-min (`uMAG` worker, `610 MB` RSS), peak RSS `2.23 GB`; `c87-stage` (`CPUQuota=400%`, `MemoryMax=32G`, `nice 19`, `3` recycled workers with `max_tasks_per_child=250`) took `14,119.4` s (`3 h 55 min`) wall, `220.901` CPU-min self (`2.70 GB` RSS); `c87-reconcile` (`s2_check_inputs.py` byte pass) took `19.91` CPU-min (`1.21 GB` RSS). Every process stayed well under the `300` CPU-minute limit (`AGENTS.md` §7.2).

## 1. Pod-Resident Input Archives & On-Disk Wrapper Structure

All input archives under `/stor/work/Wilke/if2453/all-the-data/jgi_20261118/` were harvested on 2026-11-18/19 (card `C-48a`, `planning/measurements/2026-11-18_jgi_harvest_inventory.md`) and accessed strictly read-only (`:ro`):

| Input archive / index | Bytes (*measured*) | SHA-256 (*measured*, card C-48a / C-74 / C-75) | Structure (*measured*, `ZipFile.infolist()` / `tarfile`) |
|---|---:|---|---|
| `jgi_iuig_genomes_20261118.zip` | `524,828,206,618` | `397fe93cfc07b7c1187240d796763ec93a4817f20c747d97eb78fdb41e954f7f` | `2` ZIP entries: `Download_2543422_File_Manifest.csv` (`302` B) + Deflate-compressed `Custom_MPIIUIG/2023-11-14_1.0/aggregate_tarball.tar` (`525,205,166,080` B uncompressed), holding `143,796` `<taxon_oid>.tar.gz` bundles |
| `jgi_img_umag_prok_bins_20261118.zip` | `408,896,139,104` | `bbbe2fa9b1c7887d661383ce93229390c2f8fe7f913729c9cadea3426a240bfb` | `2` ZIP entries: `Download_2543407_File_Manifest.csv` (`321` B) + Deflate-compressed `CusIMGUMAG/IMG_Data/IMG_prokaryote_binning_files.tar` (`409,281,208,320` B uncompressed), holding `192,260` `./<mg_id>_<bin>.tar.gz` bundles |
| `jgi_inovirus_umag_iuig_meta_20261118.zip` | `88,697,809` | `f8b546e2a50ed9074f3e14bf75852d7b95dd59cdf0ee9b39a259a723e286e393` | Holds `img_genomes_tarballs_metadata_111023.txt` (`143,796` data rows) and `unrestricted_metagenomes_for_bin_tarballs.02062024.txt` (`18,912` metagenome rows) |
| `jgi_gems_meta_prophages_otus_20261118.tar` | `44,091,832,320` | `ffc147dbdf956710225e6f9bb0de8b5cd88c3d5290097ce3597da5f4bb0ad06c` | `gems/otus/genome_to_otu.tsv` (`576,561` rows) and `gems/otus/genome_to_checkm.tsv` |
| `jgi_gold_and_img_m_metadata_20261118.tar` | `405,145,600` | `d63d7bba593d96d32dddafbbecf2ccf20ca38aeba2fbe07597ada1dfce95a7bc` | `img_m/img_m_{archaea,bacteria_finished,microbiomes}_70col.tsv.gz` (`126,357` rows across 5 slices) |
| `/tmp/if2453/c74/index/{gtdb_r232,globdb_r232,globdb_r232_dictionary,globdb_r226_keys}.parquet` | `51,833,962` (*derived* sum) | Card `C-74` (`2026-11-26_r232_index.md`) | `715,230` GTDB R232 + `317,533` GlobDB r232 reps + `1,729,886` GlobDB r232 dictionary rows |

**Note on archive member counts vs prior `C-84` / `CATALOG_GOAL.md` estimates (`STAGING_WORKER_RULES.md` §Bytes):**
1. Both JGI Data Portal `.zip` files wrap a single giant Deflate-compressed `.tar` (`aggregate_tarball.tar` and `IMG_prokaryote_binning_files.tar`) rather than storing individual genomes as top-level ZIP directory entries. `c87_jgi_ug_inputs.py` supports both direct per-genome ZIP layouts (`O(1)` seek via `zf.open(info)`) and nested JGI `zip -> outer .tar -> inner <id>.tar.gz` layouts (streaming the outer `.tar` sequentially while skipping inner `.tar.gz` decompression for unselected entries).
2. **`IUIG` (`143,796` total members, `126,346` prokaryotes, *measured*):** `img_genomes_tarballs_metadata_111023.txt` lists `143,796` rows, matching all `143,796` `<taxon_oid>.tar.gz` members (`0` missing). Of these, `126,346` are prokaryotes (`123,665` `Bacteria` + `2,681` `Archaea`: `19,147` `Finished`, `105,713` `Permanent Draft`, `1,483` `Draft`, `3` blank status, *measured*, `2026-11-19_jgi_mags_iuig_dedup.md`), and `17,450` are non-prokaryotes (`15,544` `Viruses`, `1,104` `Plasmid:Bacteria`, `42` `Plasmid:Archaea`, `39` `Plasmid:Eukaryota`, `721` `Eukaryota`, *measured*). Prior notes (`C-84`) had quoted `~100,455` uncultivated prokaryote rows as an estimate before full archive enumeration; all `126,346` prokaryote assemblies are censused here, while the `17,450` non-prokaryote entries are recorded in `iuig_non_prok_skipped.tsv` (`Viruses`/`Plasmid:*` covered by `IMG/VR` and `IMG/PR`; `721` `Eukaryota` deferred to card `C-60`).
3. **`IMG-UMAG` (`192,260` prokaryote bins across `18,878` metagenomes, *measured*):** Prior notes (`C-84` / `2026-11-19_jgi_mags_iuig_dedup.md`) explicitly estimated `~131,807` bins by extrapolating from a `500`-bin prefix probe of `c48d_stream.py`. Full streaming of `IMG_prokaryote_binning_files.tar` found exactly `192,260` unrestricted `HQ`/`MQ` prokaryote bins across `18,878` of the `18,912` metagenomes in `unrestricted_metagenomes_for_bin_tarballs.02062024.txt`. Because each `<mg_id>_<bin>.tar.gz` bundle carries `mbin_datafile_<mg_id>.txt` describing every bin of that metagenome (`Bin Quality`, `Bin Lineage`, `GTDB-TK lineage`, CheckM completeness/contamination, base count), caching by `mg_id` opened only `19,924` of `192,260` inner `.tar.gz` bundles (`10.36%`, *derived*) during census (`42.439` CPU-min, *measured*).

## 2. Full Census Cross-Joined Against `GTDB R232`, `GlobDB r232`, and Held Stages

All `318,606` prokaryote entries (`census_report.json`, sha256 `c15a78c1178ff32a8d5f18d338e5213d2e140cd666be1bf7e68740f24c671f88`, *measured*) were cross-joined against `/tmp/if2453/c74/index/` (`GTDB R232` assembly/BioSample/WGS/taxid/species + `GlobDB r232` `IMG`/`GEM` original IDs and `OTU-*` cluster IDs via `gems/otus/genome_to_otu.tsv`) and held GlobDB r226 stages (`G1`, `GART1`, `GART2`, `GRARE1..3` via `globdb_r226_keys.parquet`):

| Metric / Classification | `IUIG` prokaryotes (*measured*) | `IMG-UMAG` bins (*measured*) | Combined total (*derived*) |
|---|---:|---:|---:|
| **Total prokaryote archive entries** | **`126,346`** | **`192,260`** | **`318,606`** |
| **Total uncompressed bases (`bp`)** | `481,612,250,566` (`481.61 Gbp`) | `511,194,844,676` (`511.19 Gbp`) | `992,807,095,242` (`992.81 Gbp`) |
| Domain: `Archaea` | `2,681` | `9,745` | `12,426` |
| Domain: `Bacteria` | `123,665` | `182,515` | `306,180` |
| Quality: `HQ` (`>=90%` comp, `<=5%` cont, or `Finished`) | `69,000` | `56,814` | `125,814` |
| Quality: `MQ` (`>=50%` comp, `<=10%` cont, or unscored `>=100 kbp` draft) | `56,268` | `135,446` | `191,714` |
| Quality: `LQ` (`<100 kbp` unscored `IUIG` fragments) | `1,078` | `0` | `1,078` |
| `r232`: **`not_in_r232` (net-new to GTDB R232 & GlobDB r232)** | **`8,053`** (`6.37%`) | **`127,947`** (`66.55%`) | **`136,000`** (`42.69%`) |
| `r232`: `globdb_r232_collapsed_gem_otu` (in GlobDB r232 dictionary via GEM `OTU-*`) | `53,940` | `39,144` | `93,084` |
| `r232`: `gtdb_r232_taxid_match` (strain/isolate match to GTDB R232 NCBI taxid) | `42,530` | `0` | `42,530` |
| `r232`: `in_gtdb_r232` (direct GCA/GCF/BioSample match in GTDB R232) | `14,097` | `4,908` | `19,005` |
| `r232`: `in_gtdb_and_globdb_r232` (representative in both GTDB R232 & GlobDB r232) | `5,990` | `7,006` | `12,996` |
| `r232`: `gtdb_r232_species_match` (named species binomial in GTDB R232) | `490` | `7,662` | `8,152` |
| `r232`: `in_globdb_r232` (direct `IMG`/`GEMOTU` rep in GlobDB r232) | `1,246` | `5,593` | `6,839` |
| Already in GlobDB r226 stages (`G1` / `GART1` / `GART2` / `GRARE1` / `GRARE2` / `GRARE3`) | `112 / 40 / 25 / 25 / 9 / 19` (`230`) | `198 / 206 / 439 / 177 / 67 / 86` (`1,173`) | `310 / 246 / 464 / 202 / 76 / 105` (`1,403`) |

All `1,078` `IUIG` entries classified as `LQ` had no CheckM entry in `genome_to_checkm.tsv`, were not `Finished`, and measured `< 100,000` bp directly from their FASTA (`min_unscored_bp = 100,000`; across all `255,862` CheckM-scored genomes in `IUIG + uMAG`, `0` scored `HQ`/`MQ` genomes had `< 100,000` bp, *measured*).

## 3. `JGIUG1` Prioritization (`<= 35.00 Gbp`) & `JGIUG_left_for_later.tsv.gz`

Candidates were ranked deterministically (`_rank_and_partition_candidates`):
1. Strict novelty: `r232_status == "not_in_r232"` and not in `G1`/`GART1..2`/`GRARE1..3` ahead of all `r232`-covered or held entries;
2. Lineage novelty group: (0) DPANN & Asgard archaea (`archaea_dpann_asgard`), (1) all other `Archaea` (`archaea_other`), (2) CPR / `Patescibacteria` (`cpr_patescibacteria`), (3) unclassified-phylum `Bacteria` (`bacteria_unclassified_phylum`), (4) rare & candidate bacterial phyla (`bacteria_rare_or_candidate_phylum`), (5) ART-enriched phyla, (6) major environmental phyla, (7) major host/human phyla;
3. Quality within each lineage group: `HQ` (`q_ord = 0`), then CheckM-scored `MQ` (`q_ord = 1`), then unscored `>= 100 kbp` `IUIG` draft `MQ` (`q_ord = 2`), ordered by `-quality_score`, `-bp`, `genome_id`. (Sorting `HQ` before `MQ` *within* each lineage group rather than globally prevents `86.70 Gbp` of `HQ` common-phylum `Bacteria` from displacing `6,630` `MQ` CPR/`Patescibacteria` and `1,175` `MQ` DPANN/Asgard archaea whose CheckM1 completeness is capped below `90%` by lineage-specific marker loss.)

| Partition | Genomes (*measured*) | Uncompressed `bp` (*measured*) | Breakdown (*measured*) |
|---|---:|---:|---|
| **`JGIUG1` selected (`selection.tsv`)** | **`21,646`** | **`34,999,996,495`** (`35.000 Gbp`) | `100%` `not_in_r232`; `1,397` `IUIG` + `20,249` `uMAG`; `6,685` `Archaea` + `14,961` `Bacteria`; `3,336` `HQ` + `18,310` `MQ`; `archaea_dpann_asgard` `1,212`, `archaea_other` `5,473`, `cpr_patescibacteria` `6,630`, `bacteria_unclassified_phylum` `1,600`, `bacteria_rare_or_candidate_phylum` `6,731` |
| **Deferred HQ/MQ (`JGIUG_left_for_later.tsv.gz`)** | **`295,882`** | **`957,769,488,232`** (`957.770 Gbp`) | **`114,005` novel `not_in_r232` spillover (`330,478,507,298` bp = `330.48 Gbp`, `spillover_stage_cap`)** + `181,877` `r232`/held reserves (`globdb_r232_collapsed_gem_otu` `92,785`, `gtdb_r232_taxid_match` `41,850`, `in_gtdb_r232` `18,997`, `in_gtdb_and_globdb_r232` `12,977`, `gtdb_r232_species_match` `8,103`, `in_globdb_r232` `5,762`, `held_G1/GART/GRARE` `1,403`) |

Top 15 phyla in `JGIUG1` (`selection.tsv`, *measured*): `Patescibacteria` (`6,147`), `Crenarchaeota` (`1,484`), `Unclassified` (`1,285`), `Thermoplasmatota` (`1,206`), `Thermoproteota` (`1,026`), `Halobacterota` (`808`), `Nanoarchaeota` (`717`), `Bdellovibrionota` (`676`), `unclassified` (`583`), `Armatimonadota` (`582`), `Nitrospirota` (`462`), `Eremiobacterota` (`346`), `Euryarchaeota` (`261`), `Halobacteriota` (`245`), `Micrarchaeota` (`245`).

## 4. Stream-Staging, Genetic-Code Screening & Component Reconciliation

Unit `c87-stage` streamed only the `21,646` selected genomes directly from `jgi_iuig_genomes_20261118.zip` (`1,397` entries) and `jgi_img_umag_prok_bins_20261118.zip` (`20,249` entries) into BGZF `jgiug1_genomes.fna.gz` (`11,290,730,624` B = `10.52 GiB`, *measured*), prefixing every FASTA header `<genome_id>|<orig_header>` for global uniqueness (`2,938,090` records, `0` missing genomes, `0` non-IUPAC or gap-containing genomes dropped; `left_out.tsv` header-only, `35` B, *measured*).

- **Genetic code screen (`11, 25, 4`, `screen_max_bp = 500,000`):** `6,677` genomes triggered alternative-table testing (`tested_alt = true`). Across the `21,524` distinct `genome_uid`s, **`21,205` received code `11`**, **`319` received code `25`** (CPR/`Patescibacteria` UGA->Gly recoding), and **`0` received code `4`** (*measured*).
- **Transitive shared-contig component reconciliation (`DEC-27`):** `28,748` contigs appear in `>= 2` genomes within `JGIUG1` (`listing_contigs_in_several_genomes`, *measured*), including one multi-hop Patescibacteria bin component (`3300025839_42`, `3300025837_50`, `3300004212_28`, `3300004185_21`, `3300025835_40`). Union-Find component resolution (`_resolve_component_codes`, `reconcile-codes`) assigned every component the initial code of its least-`genome_uid` member (`1` genome conflict resolved: `3300025837_50` initial `11` -> component leader `3300025839_42` code `25` via `1eb26b2f8c4acf49ed0f5fac4e5a0694`, *measured*).
- **`DEC-27` repeated contig sets:** `122` genomes duplicated an earlier streamed genome's exact `contig_uid` set (`repeated_records.tsv`) and `4,973` within-genome duplicate records were collapsed by `make_inputs.Listing`, yielding **`21,524` distinct `genome_uid`s** and **`2,923,722` genome–contig memberships** (*measured*). Assembly quality distribution: `mag_mq` `17,388`, `mag_hq` `3,315`, `contig` `810`, `scaffold` `9`, `complete` `2` (*measured*).

## 5. Exact `contig_uid` Ownership (`S1..GBRMGD`) and Verification Checks

Exact 16-byte BLAKE3 `contig_uid` ownership was checked via DuckDB (`--threads 2`, `--memory-limit 8GB`) against all `13` stored root stages (`S1..SPX1`), all `19` `r26.11.4` stages (`V2..GART2`), and all `24` staged `r26.11.5` stages (`GRARE1..WET1`, `RMD1`, `GBRMGD`):

| Ownership scope | Distinct contigs (*measured*) | Distinct bases `bp` (*measured*) | Breakdown by stage (`contigs` [`bp`], *measured*) |
|---|---:|---:|---|
| **Total `JGIUG1` distinct contigs** | **`2,890,974`** | **`34,657,341,931`** | `21,524` distinct genomes (`2,938,090` FASTA records) |
| Held by root (`13` stages) + `r26.11.4` (`19` stages) | `2,271` | `59,716,939` | `G1` `1,006` (`19,061,839`), `GART2` `597` (`6,039,804`), `MV1` `338` (`13,115,933`), `PL3` `217` (`4,904,188`), `V2` `72` (`3,428,762`), `LR1` `18` (`259,413`), `SPX1` `16` (`189,751`), `P6` `3` (`5,421,348`), `GART1` `1` (`3,060`), `P3` `1` (`3,300,546`), `P5` `1` (`3,947,083`), `V1` `1` (`45,212`) |
| **New after root + `r26.11.4`** | **`2,888,703`** | **`34,597,624,992`** | `99.83%` of distinct bases (*derived*) |
| Held by `r26.11.5` (`24` stages `GRARE1..GBRMGD`) | `13,079` | `150,558,524` | `DPN1` `5,903` (`45,200,576`), `SAG3` `2,806` (`43,997,242`), `GRARE1` `2,105` (`22,710,071`), `GRARE3` `831` (`13,735,747`), `GRARE2` `610` (`12,012,371`), `SAG1` `458` (`6,862,294`), `RMD1` `147` (`1,856,639`), `MV2` `77` (`2,882,108`), `NX1` `71` (`756,930`), `SAG2` `36` (`216,462`), `MIR2` `25` (`245,113`), `MVEUKV` `6` (`56,678`), `GV1` `4` (`26,293`) |
| **Net new after all `56` stages (`S1..GBRMGD`)** | **`2,875,624`** | **`34,447,066,468`** | **`99.39%` of distinct bases (`34.447 Gbp` net new, *derived*)** |

- **`s2_check_inputs.py JGIUG1 genomes.json objects_staging.json`:** **`exit: 0`, `problems: []`** (`21,524` genomes, `2,923,722` `genome_input_rows`, codes `[11, 25]`, full byte pass over `jgiug1_genomes.fna.gz`: `2,938,090` records, `35,529,418,406` bases, `0` `contigs_in_no_object`, `0` `records_not_in_genomes`, `0` `empty_records`, `0` `gap_bytes`, `0` `other_non_iupac_bytes`, `1,194.839` s wall, `1,194.717` CPU-s = `19.912` CPU-min, `1,209,655,296` B peak RSS, *measured*).
- **`stage_est.py conf`:** `python -I scripts/stress/build_stage_pod/stage_est.py conf scripts/stress/build_stage_pod/stages/templates/JGIUG1.conf --stage JGIUG1 --tree <worktree>` exited `0` under both normal `PATH` and `PATH=/usr/bin:/bin` (*measured*).

## 6. Staged File Inventory, SHA-256 Digests & Gentle `tower` Mirror (`DEC-88 (1)`)

Every file written to `/stor/scratch/Wilke/if2453/bio-box-staging/JGIUG1/` exceeds the `64 MiB` (`67,108,864` B) floor (`AGENTS.md` §7.3). During `c87-stage`, `records.tsv` (`342,655,432` B, sha256 `2bd4a8affc5a40db4cfd1ee5c42ba45eef748a34a86c88cc78e650c5f1975b72`) was placed on `/stor/scratch/Wilke/if2453/bio-box-staging/JGIUG1/` prior to Union-Find component code reconciliation; per `AGENTS.md` §7.5 ("Never delete or overwrite anything under `/stor`"), that file was left untouched in place and the reconciled `records.tsv` (`5e34bc757431d6a0c289e5dbfeab7d5dd724fab1326d0cf42d6fbb8a7182ce07`) was written under `/tmp/if2453/c87/meta/JGIUG1/records.tsv` and mirrored to `tower`.

| File | Pod location | Bytes (*measured*) | SHA-256 (*measured*) | Mirrored to `tower:/mnt/user/research-bulk/bio-box/sources/JGIUG1/` |
|---|---|---:|---|---|
| `jgiug1_genomes.fna.gz` | `/stor/scratch/Wilke/if2453/bio-box-staging/JGIUG1/` | `11,290,730,624` | `56119e7280b1478cdf9dae7d1a741ee1e9abb711e057d14ba88dc91296bbe424` | Yes (`--streams 1`) |
| `genomes.json` | `/tmp/if2453/c87/meta/JGIUG1/` | `107,557,082` | `1852fe01a01546caca8d89781803c2e83fd08b338cde5fb741004bacd884614c` | Yes (`--streams 1`) |
| `objects.json` | `/tmp/if2453/c87/meta/JGIUG1/` | `475` | `d9f43d2b6f57caa04a160de6040d0d52c6e05fbc34d6100e1cf9a49d512d8929` | Yes (`--streams 1`) |
| `objects_staging.json` | `/tmp/if2453/c87/meta/JGIUG1/` | `496` | `18e0bc299a2d8db6f198d39483eaf17dc6436e090e6b9f56f0a9b7b4b2f4bcdd` | Yes (`--streams 1`) |
| `records.tsv` (reconciled) | `/tmp/if2453/c87/meta/JGIUG1/` | `342,655,432` | `5e34bc757431d6a0c289e5dbfeab7d5dd724fab1326d0cf42d6fbb8a7182ce07` | Yes (`--streams 1`) |
| `genetic_codes.tsv` | `/tmp/if2453/c87/meta/JGIUG1/` | `4,169,196` | `f9535ab49f786e9f2778a6d2a4aa0369803c063b85c6c1537ec194eeac4e6423` | Yes (`--streams 1`) |
| `r232_status.tsv` | `/tmp/if2453/c87/meta/JGIUG1/` | `3,131,863` | `2d35a9f74b38b6e312044391b06d624c8be51a0d85875c4510f28ac597043e20` | Yes (`--streams 1`) |
| `JGIUG_left_for_later.tsv.gz` | `/tmp/if2453/c87/meta/JGIUG1/` | `12,408,407` | `dcfd638a2d6acf10435ed00dbde608c9592bfc4018b9071dd554d99132931b8c` | Yes (`--streams 1`) |
| `later_records.tsv` | `/tmp/if2453/c87/meta/JGIUG1/` | `946,054` | `a0d4a236d2522dda6a9584e624591a379f24d427005795acbd70e67c0664fefe` | Yes (`--streams 1`) |
| `shared_records.tsv` | `/tmp/if2453/c87/meta/JGIUG1/` | `315,442` | `cf3ceec60f84c872a4a3e0a3e23a066407ae07108e4c49b8a0caddbede387d91` | Yes (`--streams 1`) |
| `repeated_records.tsv` | `/tmp/if2453/c87/meta/JGIUG1/` | `8,138` | `98af7070a40584a8e0b71da027855487ee1748982c018595169bb20435a991d6` | Yes (`--streams 1`) |
| `left_out.tsv` | `/tmp/if2453/c87/meta/JGIUG1/` | `35` | `4d6eafbf6b5ac7af077aa26511b328b958af24b446b0b93ae1b0d71705eb4a5d` | Yes (`--streams 1`) |
| `stage_report.json` | `/tmp/if2453/c87/meta/JGIUG1/` | `4,683` | `2f3d2e52c5ab3ecfa4df791c9426346f08028b385167c31abb21781d4a6ba766` | Yes (`--streams 1`) |
| `selection.tsv` | `/tmp/if2453/c87/census/` | `5,514,615` | `fe462e362e950f3c6ce6390ea9389a003ad113eccea1f491ea7e40dbc1e4e2ff` | Census table on pod |
| `iuig_non_prok_skipped.tsv` | `/tmp/if2453/c87/census/` | `871,765` | `0d3efd740debb4ef3a8c36b6a00e184630725a12b43ce7694f0dc62f48c55b68` | Census table on pod |
| `census_summary.json` | `/tmp/if2453/c87/census/` | `4,590` | `c15a78c1178ff32a8d5f18d338e5213d2e140cd666be1bf7e68740f24c671f88` | Census report on pod |
