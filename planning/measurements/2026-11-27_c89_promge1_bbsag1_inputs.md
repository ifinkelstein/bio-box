# Staging `tower`-landed `r26.11.6` cohorts `PROMGE1` and `BBSAG1` on `wilkcomp02` (card C-89, DEC-65, DEC-80, DEC-84 (2), DEC-88)

- **Date and hosts:** 2026-11-27 and 2026-11-28, `wilkcomp02` and `tower`. Every number carries a label (`AGENTS.md` §6): *measured* (on the named host, 2026-11-27 or 2026-11-28, by the command or file named), *derived* (formula and inputs), or *estimate* (source).
- **Script and tests:** `scripts/stress/build_stage_pod/c89_promge1_bbsag1_inputs.py` (`census`, `promge1-stage`, `bbsag1-stage`, `bbsag1-ani`); `tests/r1/test_c89_promge1_bbsag1_inputs.py`; templates `scripts/stress/build_stage_pod/stages/templates/{PROMGE1,BBSAG1}.conf`.
- **Read only:** `tower:/mnt/user/research-bulk/bio-box/sources_incoming/{batch6_mges/{PROMGE2,GSPR,MFD_PLAS,ECE_HODOR,SAT_P4P2,ECE_METH},batch6b_ncbi/{PICJI_SAT,BORG_ONT},batch6_longread_sags/{BBSAG20,HAIMA_LR_SAG,DEATHVALLEY_SAG}}/`, the pod data root `/tmp/if2453/r10/root` (`13` stages `S1..SPX1`), `/tmp/if2453/s2/inputs/` (`19` `r26.11.4` stages `V2..GART2`), the `24` staged `r26.11.5` stages (`GRARE1..GBRMGD`), the `3` merged card `C-86` `r26.11.6` stages (`O4, GV2, HROM1` under `/tmp/if2453/c86/meta/`), card `C-74`'s index `/tmp/if2453/c74/index/`, and card `C-75`'s archaeal references (`/stor/scratch/Wilke/if2453/bio-box-staging/GRARE/globdb_r226_archaea.fna.gz`, `/tmp/if2453/c75/meta/archaea/archaea_list.tsv`). Nothing was written to `/tmp/if2453/r10/root`, `/tmp/if2453/s2/inputs/`, or `tower`'s `sources_incoming/`; nothing on `/stor` or `tower` was deleted or overwritten.
- **Written (`AGENTS.md` §7.3 `64 MiB` floor strictly enforced):** raw downloads and staged FASTA / `records.tsv` objects `>= 67,108,864` B (`64 MiB`) under `/stor/scratch/Wilke/if2453/bio-box-staging/{PROMGE1,BBSAG1}/` (`5` staged files on `/stor`, smallest `304,154,499` B, *measured* via `promge1_stage_report.json` and `bbsag1_stage_report.json`); `< 64 MiB` staged FASTA (`bbsag1_deathvalley.fna.gz` `31,760,498` B) and stage metadata under `/tmp/if2453/c89/meta/{PROMGE1,BBSAG1}/`; mirrored gently (`--streams 1`, `--dry-run` first, sequential when `pgrep -fl relay_to_` on `tower` was empty) to `tower:/mnt/user/research-bulk/bio-box/sources/{PROMGE1,BBSAG1}/`. Temporary `/tmp/if2453/c89/{dl,work}` scratch was removed immediately after each stage finished (`381 GiB` free on `/tmp/if2453`, *measured*).

---

## 1. Gentle single-stream raw relay from `tower` to `wilkcomp02` (`DEC-88 (1)`)

Before every `relay_to_pod.sh` invocation, verified `ssh root@tower pgrep -fl relay_to_` (`0` active relays, *measured*) and `ssh root@tower uptime`, ran `scripts/orchestrator/relay_to_pod.sh --dry-run --streams 1` first, and relayed strictly `--streams 1` one directory list at a time (routing files `>= 67,108,864` B to `/stor/scratch/Wilke/if2453/bio-box-staging/<STAGE>/dl/` and files `< 67,108,864` B to `/tmp/if2453/c89/dl/<STAGE>/`):
- **`PROMGE1` (`8` source directories, `28` non-checksum files):** `36,804,768,535` B (`36.805 GB`, *derived* sum across `PROMGE2` `5` files, `GSPR` `3` files, `MFD_PLAS` `3` files, `ECE_HODOR` `6` files, `SAT_P4P2` `1` file, `ECE_METH` `1` file, `PICJI_SAT` `7` files, `BORG_ONT` `2` files; *measured*, `wilkcomp02`, 2026-11-27, `/tmp/if2453/c89/meta/PROMGE1_census/census.json`), `28` verified, `0` failed.
- **`BBSAG1` (`3` source directories, `90` non-checksum files):** `77,521,133,831` B (`77.521 GB`, *derived* sum of `BBSAG20` `86` files `75,877,683,870` B + `HAIMA_LR_SAG` `3` zip archives `1,612,093,240` B + `DEATHVALLEY_SAG` `1` zip archive `31,356,721` B; *measured*, `wilkcomp02`, 2026-11-28, `/tmp/if2453/c89/meta/BBSAG1_census/summary.json`), `90` verified, `0` failed.

---

## 2. Raw archive census & cohort verification (`STAGING_WORKER_RULES.md`)

Inspected every relayed archive and flatfile on `wilkcomp02` by magic bytes (`gzip`, `zip`, `tar`, `fasta`, `text`, `binary`), member paths, Prodigal `-d` header signatures, and IUPAC nucleotide alphabet (`c89_promge1_bbsag1_inputs.py census` -> `/tmp/if2453/c89/meta/{PROMGE1_census,BBSAG1_census}/{census.json,summary.json}`, *measured*, `wilkcomp02`, 2026-11-27 and 2026-11-28):

| Stage / subcohort | Raw contents (*measured*) | What enters staging vs left out / skipped (`left_out.tsv` / `skipped_members.tsv`) |
|---|---|---|
| `PROMGE1` / `PROMGE2` | `mge100.fna.gz` (`33,912,610,622` B): `5,665,904` DNA records (`107,714,624,028` bp = `107.715 Gbp` *derived*, `0` non-IUPAC) across `855,323` host genomes in `genomes_taxonomy.tsv.gz` (`673,876` `proGenomes3` isolates [`5,143` `specI` species reps + `668,733` non-rep clinical resequencings] + `181,447` `SPIRE` metagenomic MAGs), plus `mge100_clusters.tsv.gz`, `nested_mges.tsv.gz`, `genomes_quality.tsv.gz`. Raw MGE breakdown: `246,897` `cellular_only` (`mge=cellular:1` host `xer` recombinase loci, `3,549,115,070` bp) + `5,419,007` non-cellular MGEs (`1,030,411` `mge_type=nested` composites `33,150,511,754` bp + `4,388,596` `non-nested` elements `71,014,997,204` bp) | **1,915,224 non-cellular, non-nested MGEs staged** across all `181,447` `SPIRE` MAGs and all `5,143` `proGenomes3` `specI` species representatives, plus `20,923` `mge100` cluster representatives rescued via `mge100_clusters.tsv.gz` (orchestrator-approved Option 3b keeping single `PROMGE1 <= 40 Gbp`; full `mge100.fna.gz` stays archived on `tower`); **`3,750,680` records summarized in `left_out.tsv`** (`246,897` `promge2_cellular_recombinase_non_mge`, `1,030,411` `promge2_nested_composite`, `2,473,372` `promge2_pg3_nonrep_host`); `4` `.tsv.gz` tables recorded in `skipped_members.tsv` |
| `PROMGE1` / `GSPR` | `all_soil_plasmids_derep.fna.gz` (`98,281` soil plasmids, `2,528,128,530` bp) + `soil_plasmid_db_parquet.tar.gz` (`19` Parquet/table members) + `accessions.txt` | **all 98,281 DNA records staged** (`98,277` distinct after `DEC-27`); `2` non-FASTA files recorded in `skipped_members.tsv` |
| `PROMGE1` / `MFD_PLAS` | `mfd_plasmids_drep_and_members.fasta` (`9,265` ONT plasmids, `243,392,607` bp) + `pfd_metadata.tsv` (`44` `Inverted terminal repeat` linear plasmids, rest `Circular` / `Direct terminal repeat`) + `MFD01716.bin.1.47.fasta` (`2,844,778` B: `3,124` Prodigal `-d` single-gene CDS nucleotide sequences `# start # end # strand # ID=...`) | **all 9,265 plasmid contigs staged** (`9,133` distinct; `44` TIR linear plasmids preserved untrimmed & linear per `DEC-65`); **Prodigal `.ffn` `MFD01716.bin.1.47.fasta` caught by census (`prodigal_cds_ffn`) and skipped** (`prodigal_cds_ffn_not_plasmid_contig` in `skipped_members.tsv`) |
| `PROMGE1` / `ECE_HODOR` | `Hodor.38.genomes.fasta` (`38` curated >300 kb megaplasmids), `Hodor.prophages.fasta` (`43` curated prophages), `prophage01..03.all.pure.fa` (`10,394` prophages), `draft.Hodor.genomes.tar.gz` (`1,625` draft Hodor MAG bins, `28,520` contigs): `38,995` DNA records (`612,519,981` bp) | **all 38,995 DNA records staged** (`11,658` distinct `kind="plasmid"` genomes after `DEC-27` deduplication across `Hodor.prophages.fasta` and `prophage01..03`) |
| `PROMGE1` / `SAT_P4P2` | `P4P2_DereplicatedGenomes.zip`: `4,596` DNA replicons in `Replicons/*.fna` (`3,085` P4 satellites + `1,511` P2 prophages, `85,919,008` bp), `4,596` protein `.prt` files in `Proteins/*.prt`, `9,203` `__MACOSX`/`.DS_Store` entries | **all 4,596 DNA replicons staged**; `4,596` `.prt` protein files and `9,203` `macosx` entries recorded in `skipped_members.tsv` |
| `PROMGE1` / `ECE_METH` | `Mp_and_ECE_genomes.fa.gz`: `26` PacBio HiFi assemblies (`13` `HMp-Plasmid-like` ECEs + `1` `Mp-ICE` + `12` host *Methanoperedens* chromosomes) | **all 14 ECEs (`13` plasmid-like + `1` ICE) staged**; **`12` host *Methanoperedens* chromosomes recorded in `left_out.tsv`** (`host_archaeal_chromosome_not_ece`) |
| `PROMGE1` / `PICJI_SAT` | `MANIFEST_PICJI_SAT_PRJNA1424553.tsv` + `6` INSDC complete assemblies (`17` replicons: `6` chromosomes + `11` plasmids, `29,450,714` bp) | **all 17 replicons staged** byte-identically (`DEC-65 (2)`, `redistributable: true`) |
| `PROMGE1` / `BORG_ONT` | `MANIFEST_BORG_ONT_PRJNA1119519.tsv` + `1` INSDC Nanopore Borg assembly (`GCA_040223915.1`, `4` linear contigs, `2,426,186` bp) | **1 Borg genome (`4` linear contigs) staged** byte-identically (`DEC-65 (2)`, `redistributable: true`) |
| `BBSAG1` / `BBSAG20` | `86` files (`75,877,683,870` B): `summary.tsv` (`32,395` rows = **`18,071` MIMAG MQ/HQ** [`completeness >= 50.0`, `contamination < 10.0`: `17,202` SAGs (`2,915` HQ + `14,287` MQ) + `869` MAGs (`433` HQ + `436` MQ)] + **`14,324` LQ** [`13,649` SAG:LQ + `675` MAG:LQ]), `83` `*sag.tar.gz` + `mag.tar.gz` (`32,395` `.fna.gz` + `.faa.gz`/`.ffn.gz`/`.gff.gz` annotations; disambiguated `(sample, sagmag, number)` so `QLF001sag-00001` and `QLF001mag-00001` never collide), plus `bbsag20_plasmid.fna.gz` (`9,045` extracted plasmid records: `6,818` from MQ/HQ + `2,227` from LQ) | **all 18,071 MIMAG MQ/HQ SAGs and MAGs staged** (`100%` 1-to-1 match against `summary.tsv`); **`14,324` LQ SAGs/MAGs recorded in `left_out.tsv`** (`bbsag20_low_quality`); `.faa.gz`/`.ffn.gz`/`.gff.gz` annotations and `bbsag20_plasmid.fna.gz` (`plasmid_subset_of_bbsag20_genomes`) recorded in `skipped_members.tsv` |
| `BBSAG1` / `HAIMA_LR_SAG` | `3` zip archives (`1,612,093,240` B): `SAGs.zip` (`4,069` cold-seep microfluidic SAGs), `long-read-MAGs.zip` (`568` long-read MAGs incl. `200` single-contig `cbin` `cMAG`s), `short-read-MAGs.zip` (`313` short-read MAGs) = `4,950` DNA assemblies (`4,700,204` contigs, `5,174,108,436` bp) | **all 4,950 SAGs and MAGs staged** (`195` single-contig `cbin` `cMAG`s confirmed circular under `DEC-65 (1)`/`(3)`, `4` `partial` and `1` `long` terminal repeats kept linear per `DEC-65 (3)`) |
| `BBSAG1` / `DEATHVALLEY_SAG` | `BLM1-MAGs_SAGs_2011_2021.zip` (`31,356,721` B): `57` DNA assemblies (`7,336` contigs, `101,600,100` bp) + `2` duplicate entries of `information regarding the MAGs and SAGs of BLM1.xlsx` | **all 57 groundwater SAGs and MAGs staged** in separate `bbsag1_deathvalley.fna.gz` (`redistributable: false` under `DEC-80 (2)-(3)`); `2` `.xlsx` entries recorded in `skipped_members.tsv` |

---

## 3. Stage `PROMGE1`: Prokaryotic MGEs, plasmids, satellites, and Borgs (`DEC-65`, `DEC-80 (3)`, `DEC-88`)

**Rules (`stage_promge1`):**
- All `8` subcohorts stage as `kind="plasmid"`, `genetic_code=11` (called in metagenomic `meta` mode by `biobox.l1.genecall`).
- Non-INSDC single-contig circular/DTR candidates (`GSPR`, `MFD_PLAS` circular/DTR, `ECE_HODOR` single-record elements, `ECE_METH` plasmid-like ECEs) apply `DEC-65 (1)` (`20..1,000` bp exact DTR trimmed once and marked circular: `2,725` exact DTR trims removing `326,328` bp; `599` `partial` and `45` `long` left linear per `DEC-65 (3)`). `MFD_PLAS` `Inverted terminal repeat` / `Linear` replicons (`44` plasmids) and INSDC `PICJI_SAT` / `BORG_ONT` records (`DEC-65 (2)`) keep exact deposited sequence bytes.
- Split into two staged FASTA objects by redistributability (`DEC-80 (3)`): `promge1_promge2.fna.gz` (`redistributable: false`) and `promge1_plasmids_eces.fna.gz` (`redistributable: true`, `7` open CC-BY / INSDC cohorts).

| Metric | Value | Label & source |
|---|---:|---|
| Qualifying DNA records written | **2,066,396** (`PROMGE2` 1,915,224; `GSPR` 98,281; `ECE_HODOR` 38,995; `MFD_PLAS` 9,265; `SAT_P4P2` 4,596; `PICJI_SAT` 17; `ECE_METH` 14; `BORG_ONT` 4) | *measured* (`wilkcomp02`, 2026-11-28, `promge1_stage_report.json`) |
| MGE class breakdown (`records.tsv`) | `promge2:is_tn` 883,237; `promge2:mi` 578,532; `promge2:phage_like` 290,175; `promge2:ce` 107,184; `soil_plasmid` 98,281; `promge2:phage` 50,736; `hodor_draft` 28,520; `hodor_prophage` 10,437; `mfd_plasmid` 9,265; `promge2:integron` 5,360; `p4_satellite` 3,085; `p2_prophage` 1,511; `hodor_ece` 38; `mp_ece` 13; `picji_plasmid` 11; `picji_chromosome` 6; `borg_ont` 4; `mp_ice` 1 | *measured* (`promge1_stage_report.json`) |
| Excluded records (`left_out.tsv`) | `12` `ECE_METH` host *Methanoperedens* chromosomes + `3` compact `PROMGE2` summary rows covering `3,750,680` excluded `mge100.fna.gz` records (`246,897` cellular `3,549,115,070` bp; `1,030,411` nested `33,150,511,754` bp; `2,473,372` non-rep `pg3` `40,816,988,566` bp) | *measured* (`promge1_stage_report.json`) |
| Skipped archive/dir entries (`skipped_members.tsv`) | 13,807 (`4,596` `SAT_P4P2` `.prt` protein files, `9,203` `macosx` entries, `1` Prodigal `.ffn` `MFD01716.bin.1.47.fasta`, `7` metadata tables/archives) | *measured* (`promge1_stage_report.json`) |
| Exact `genome_uid` repeats (`repeated_records.tsv`, `DEC-27`) | 578 | *measured* (`promge1_stage_report.json`) |
| Distinct genomes (`genomes.json`) / `genome_input` rows | **2,038,920** (`PROMGE2` 1,915,224; `GSPR` 98,277; `ECE_HODOR` 11,658; `MFD_PLAS` 9,133; `SAT_P4P2` 4,596; `PICJI_SAT` 17; `ECE_METH` 14; `BORG_ONT` 1; `7,270` circular) / **2,065,799** (`2,756` shared within-stage contig occurrences) | *measured* (`promge1_stage_report.json`) |
| INSDC host/assembly `r232_status.tsv` (`10,823` accessions, none cut) | `in_gtdb_r232`: 7,189; `in_gtdb_and_globdb_r232`: 2,691; `not_in_r232`: 943 (incl. all `6` `PICJI_SAT` `GCF_058350605..` assemblies; `BORG_ONT` `GCA_040223915.1` is in `GlobDB r232` / `GTDB R232`) | *measured* (`promge1_stage_report.json`) |
| Distinct contigs / distinct bases | 2,063,640 / 33,663,421,208 bp (33.663 Gbp, *derived*) | *measured* (`promge1_stage_report.json`) |
| Held across 59 earlier stages (`S1..HROM1`, exact `contig_uid`, build order) | **33,620 contigs, 1,495,125,808 bp** (*derived* sums across `27` stages): `S1` 1 (`54,286` bp), `P1` 17 (`3,720,999`), `P2` 81 (`20,608,914`), `P3` 198 (`35,887,900`), `P4` 90 (`20,917,553`), `P5` 73 (`10,371,468`), `P6` 15 (`5,132,241`), `H1` 2 (`742,912`), `PL1` 1,059 (`175,498,053`), `C5` 2 (`123,675`), `V1` 19 (`646,538`), `SPX1` 19 (`346,821`), `V2` 88 (`3,345,048`), `PL2` 118 (`22,912,796`), `PL3` 31,096 (`1,174,382,866`), `MGE1` 17 (`59,279`), `SAT1` 3 (`38,101`), `ECE1` 13 (`2,515,184`), `MGE2` 29 (`884,297`), `MV1` 104 (`4,471,525`), `V3` 94 (`3,187,946`), `LR1` 470 (`9,009,698`), `GART1` 1 (`6,120`), `LR2` 5 (`114,934`), `LRP1` 1 (`4,148`), `MV2` 4 (`128,509`), `WET1` 1 (`13,997`); all `32` other stages (incl. `O4, GV2, HROM1`) `0` | *measured* (`promge1_stage_report.json`) |
| **Net-new contigs / net-new bases** | **2,030,020 contigs / 32,168,295,400 bp** (`32.168 Gbp <= 40 Gbp` single-stage cap, `> 1.5 Gbp` -> `DEPTH=called`, `CRISPR=yes`) | *measured* (`promge1_stage_report.json`) |
| `s2_check_inputs.py PROMGE1` (full byte pass) | `exit: 0`, `problems: []`, `genomes: 2,038,920`, `genome_input_rows: 2,065,799`, `kinds: ["plasmid"]`, `codes: [11]` | *measured* (`promge1_stage_report.json`) |

---

## 4. Stage `BBSAG1`: Oral/gut SAG-gel, Haima cold-seep, and Death Valley single-cell & MAG cohorts (`DEC-65`, `DEC-80`, `DEC-84 (2)`, `DEC-88`)

**Rules (`stage_bbsag1` & `bbsag1_ani`):**
- Filters `BBSAG20` via `summary.tsv` (`32,395` rows) to the **`18,071` MIMAG MQ/HQ genomes** (`completeness >= 50.0` and `contamination < 10.0`: `17,202` SAGs + `869` MAGs), recording all `14,324` LQ genomes in `left_out.tsv` (`bbsag20_low_quality`) and skipping `bbsag20_plasmid.fna.gz` (`plasmid_subset_of_bbsag20_genomes`); stages all `4,950` `HAIMA_LR_SAG` cold-seep SAGs and long/short-read MAGs (`DEC-65` DTR screen on the `200` single-contig `cbin` `cMAG`s: `195` circular, `4` `partial` and `1` `long` kept linear) and all `57` `DEATHVALLEY_SAG` groundwater SAGs and MAGs.
- Screens every genome across prokaryotic genetic codes `11`, `4` (`UGA -> Trp`), and `25` (`UGA -> Gly`) on up to `500` kb of the longest contigs (`screen_prokaryote_code`, `4` recycled workers): human oral *Saccharimonadia* (`g__Nanosynbacter`, `g__Saccharimonas`) natively pass standard code `11` (`cd_11 ~ 0.91–0.94`, `mean_orf_11 ~ 850–940` bp); `QLS010sag-00146` (*Absconditabacterales* `g__HOT-345`, `2,476,402` bp, `534` contigs) screens as `uga_recoded_code25` (`cd_11 = 0.66533 -> cd_25 = 0.88988`, `mean_orf_11 = 371.05 -> 895.87` bp, `3.626` internal `TGA`/kb) and shares `5` minor co-binned contigs (`code_conflicts = 5`, out of `133,851` contigs shared across `>= 2` genomes) with standard `code 11` `QLS010` SAGs of smaller `genome_uid`, resolving to `code 11` under `DEC-27` (`l0seq.inputs._rows` shared-contig agreement invariant).
- Net-new bases across all `23,078` genomes after deduplicating against the `60`-stage prefix (`59` prior stages + `PROMGE1`) measure **`39,968,407,467` bp (`39.968 Gbp <= 40.000 Gbp`)**, fitting cleanly in a single stage `BBSAG1` (`stages/templates/BBSAG1.conf`) without requiring a two-stage split.
- Split into two staged FASTA objects by redistributability (`DEC-80 (2)-(3)`): `bbsag1_sags_mags.fna.gz` (`BBSAG20` CC BY 4.0 + `HAIMA_LR_SAG` MIT, `redistributable: true`) on `/stor`, and `bbsag1_deathvalley.fna.gz` (`DEATHVALLEY_SAG` unstated Zenodo licence, `redistributable: false`, `31,760,498` B `< 64 MiB`) on `/tmp`.

| Metric | Value | Label & source |
|---|---:|---|
| Input genomes staged | **23,078** (`BBSAG20` 18,071 MQ/HQ; `HAIMA_LR_SAG` 4,950; `DEATHVALLEY_SAG` 57) | *measured* (`wilkcomp02`, 2026-11-28, `bbsag1_stage_report.json`) |
| Genomes left out (`left_out.tsv`) | **14,324** (`bbsag20_low_quality`: `13,649` LQ SAGs + `675` LQ MAGs from `BBSAG20`) | *measured* (`bbsag1_stage_report.json`) |
| Archive/member groups skipped (`skipped_members.tsv`) | 171 (`84` tarball non-`.fna.gz` annotation summaries, `84` tarball LQ `.fna.gz` summaries, `1` `bbsag20_plasmid.fna.gz`, `2` duplicate `.xlsx` entries in `DEATHVALLEY_SAG`) | *measured* (`bbsag1_stage_report.json`) |
| Exact `genome_uid` repeats (`repeated_records.tsv`, `DEC-27`) | 0 | *measured* (`bbsag1_stage_report.json`) |
| Distinct genomes (`genomes.json`) / `genome_input` rows | **23,078** (`195` circular `cbin` `cMAG`s) / **10,091,780** (`133,851` contigs appearing in `>= 2` genomes; `0` repeated within a genome) | *measured* (`bbsag1_stage_report.json`) |
| Genetic codes (`genetic_codes.tsv`) | **`11`: 23,078** (`1` *Absconditabacterales* `QLS010sag-00146` screened `uga_recoded_code25`, resolved to `11` across `5` shared-contig conflicts under `DEC-27`) | *measured* (`bbsag1_stage_report.json`, `genetic_codes.tsv`) |
| `r232_status.tsv` counts (none removed) | `not_visible`: **23,078** (study-local SAG/MAG identifiers carry no `GCA_`/`GCF_` assembly accessions) | *measured* (`bbsag1_stage_report.json`) |
| `skani >= 96%` ANI / `>= 50%` AF hits (`ani_hits.tsv`, recorded not cut per `DEC-84 (2)`) | **5,487 of 23,078** queries (`23.78%`, *derived*; `23,076` sketched): `HROM1` **5,451** (`5,459` group hits: human oral `BBSAG20` single cells matching `HROM1` oral species reps, kept for strain-level MGE/defense islands per `DEC-84 (2)`), `GBRMGD` **13** (`14` group hits), `NX1` **9** (`11`), `SAG2` **5** (`8`), `SAG3` **4** (`8`), `DPN1` **3** (`3`), `WET1` **2** (`2`), `SAG1` `0`, `r226_archaea` `0` (`1,956.5` CPU-s = `32.6` CPU-min *derived* sum across sketch + `8` `skani dist` calls) | *measured* (`bbsag1_ani_report.json`, `ani_hits.tsv`) |
| Distinct contigs / distinct bases | 9,738,658 / 39,968,467,080 bp (39.968 Gbp, *derived*) | *measured* (`bbsag1_stage_report.json`) |
| Held across 60 earlier stages (`S1..PROMGE1`, exact `contig_uid`, build order) | **36 contigs, 59,613 bp** (*derived* sums): `P6` 1 (`1,665` bp), `PL2` 1 (`1,141`), `MGE1` 21 (`31,903`), `LRP1` 2 (`4,594`), `SAG3` 2 (`117`), `HROM1` 9 (`20,193`); all `54` other stages (incl. `PROMGE1`) `0` | *measured* (`bbsag1_stage_report.json`) |
| **Net-new contigs / net-new bases** | **9,738,622 contigs / 39,968,407,467 bp** (`39.968 Gbp <= 40 Gbp` single-stage cap, `> 1.5 Gbp` -> `DEPTH=called`, `CRISPR=yes`) | *measured* (`bbsag1_stage_report.json`) |
| `s2_check_inputs.py BBSAG1` (full byte pass) | `exit: 0`, `problems: []`, `genomes: 23,078`, `genome_input_rows: 10,091,780`, `kinds: ["prokaryote"]`, `codes: [11]` | *measured* (`bbsag1_stage_report.json`) |

---

## 5. Staged files (`wilkcomp02` and `tower` mirror)

Every staged file on `/stor/scratch/Wilke/if2453/bio-box-staging/{PROMGE1,BBSAG1}/` (excluding raw `dl/` inputs) is `>= 67,108,864` B (`AGENTS.md` §7.3; smallest is `PROMGE1/records.tsv` at `304,154,499` B, *measured*). All byte counts and SHA-256 digests below are *measured* (`ids.sha256_file` on `wilkcomp02`, 2026-11-28, `promge1_stage_report.json`, `bbsag1_stage_report.json`, `bbsag1_ani_report.json`).

### 5.1 `PROMGE1` (`13` files, `11,398,066,549` B *derived*)

| File | Pod directory | Bytes (*measured*) | sha256 (*measured*) |
|---|---|---:|---|
| `promge1_promge2.fna.gz` | `/stor/scratch/Wilke/if2453/bio-box-staging/PROMGE1/` | 9,719,767,929 | `5975a3bc283722dd0ff6c231b99da74750410bf8983f67c47ca6409e579ea81b` |
| `promge1_plasmids_eces.fna.gz` | `/stor/scratch/Wilke/if2453/bio-box-staging/PROMGE1/` | 1,083,491,406 | `c973710c6ea537cc6772aa633f404ad6fdc71f8240ef02d2e18a2f4ef5d51442` |
| `records.tsv` | `/stor/scratch/Wilke/if2453/bio-box-staging/PROMGE1/` | 304,154,499 | `ee786f9f5e44bced9f0f84b2180d021a5c578ca9b3ac7b35b7009d801e1f0aa1` |
| `genomes.json` | `/tmp/if2453/c89/meta/PROMGE1/` | 284,639,245 | `e3671c170904f72f9e6988008ffc9857abe257c5d5e23efe21dc7830bf3fdfef` |
| `objects.json` | `/tmp/if2453/c89/meta/PROMGE1/` | 1,344 | `f7d0366bacb325f4ef2a721ecd3ee40f453794ac24d37809d93e9941a8bcb343` |
| `objects_staging.json` | `/tmp/if2453/c89/meta/PROMGE1/` | 1,386 | `1077eecde010e5edf1a8d144a7b86ff76946ac4690055f4e5f80b9d298825f20` |
| `shared_records.tsv` | `/tmp/if2453/c89/meta/PROMGE1/` | 2,847,277 | `7e8d3744963f1769e532f88f57816c0caf39f734e47adf0049a591667be17cfb` |
| `skipped_members.tsv` | `/tmp/if2453/c89/meta/PROMGE1/` | 1,555,362 | `066ddd3d45307ad836478be916354af062e87386bb41a2ee9197479cfabd65ad` |
| `r232_status.tsv` | `/tmp/if2453/c89/meta/PROMGE1/` | 1,143,250 | `4dc65995bf71abf61220c003c56c0931bc1328116d67e3dae0e573900ee53417` |
| `trims.tsv` | `/tmp/if2453/c89/meta/PROMGE1/` | 415,463 | `292343c99f839c85473f35bb135c6a397520a98cf693dd49446aaf04bf86c59b` |
| `repeated_records.tsv` | `/tmp/if2453/c89/meta/PROMGE1/` | 47,066 | `95b18652c48d237c4ce491871754dc49e937f7d45293bf4ab10517a0e2429f06` |
| `left_out.tsv` | `/tmp/if2453/c89/meta/PROMGE1/` | 2,282 | `eeee2eb950150ab508650ba88bd2e45a97723acb5e6fb283630f913da25a379c` |
| `later_records.tsv` | `/tmp/if2453/c89/meta/PROMGE1/` | 40 | `35715d5e3e9f7b84789769b5e7130ae87a7e93e2adfa6cae8ac0222b940589fb` |

### 5.2 `BBSAG1` (`16` files, `15,244,892,455` B *derived*)

| File | Pod directory | Bytes (*measured*) | sha256 (*measured*) |
|---|---|---:|---|
| `bbsag1_sags_mags.fna.gz` | `/stor/scratch/Wilke/if2453/bio-box-staging/BBSAG1/` | 13,686,495,129 | `8cfa9002d5b18a34556c617ca14fd121702cd71724e4080688c7459ed8ed4c0e` |
| `records.tsv` | `/stor/scratch/Wilke/if2453/bio-box-staging/BBSAG1/` | 1,150,873,874 | `d2df7009bb9c54eb989d0353aa2b5e14e375e3604db0c7668989a80771b21306` |
| `bbsag1_deathvalley.fna.gz` | `/tmp/if2453/c89/meta/BBSAG1/` (`< 64 MiB`) | 31,760,498 | `98c79f14039d8af56ab851ae292534e8e127d96da65ecaa42d1b3c67deb37b83` |
| `genomes.json` | `/tmp/if2453/c89/meta/BBSAG1/` | 365,780,446 | `a3db25d5fe27f888b48afc0090b838d8d4505ff824810bbf44487a85ea1d3379` |
| `objects.json` | `/tmp/if2453/c89/meta/BBSAG1/` | 1,129 | `74c46745948787d938e2e2ac2d10c356524c23f302573faca8e3f4c71c0bc298` |
| `objects_staging.json` | `/tmp/if2453/c89/meta/BBSAG1/` | 1,149 | `9f10d4ac6a3995ca221aa8787d73219c788be0af619a4f5413ccad149876d71f` |
| `genetic_codes.tsv` | `/tmp/if2453/c89/meta/BBSAG1/` | 4,552,781 | `b6ce4525011e041900246581fa7976e38789bfd041687b6a9c55fe3c666e53bf` |
| `r232_status.tsv` | `/tmp/if2453/c89/meta/BBSAG1/` | 4,171,811 | `433466ee2f888188cb7dbe1df8bfa384c884a91e24da52044cae61871f661366` |
| `left_out.tsv` | `/tmp/if2453/c89/meta/BBSAG1/` | 946,800 | `d810400ba763d56ae17fe8bf45fae26e66fdd6f9d4896491e8b32b346d7b7d55` |
| `skipped_members.tsv` | `/tmp/if2453/c89/meta/BBSAG1/` | 15,985 | `e398ea2e1cea4007130d133a06ad659007b44b045e873afa453deb49095ff308` |
| `shared_records.tsv` | `/tmp/if2453/c89/meta/BBSAG1/` | 12,024 | `d022c7072d9660955db723ce7f8ebd32b0d522b4eb9c193427e9721aaa1fefcb` |
| `trims.tsv` | `/tmp/if2453/c89/meta/BBSAG1/` | 102 | `3b6006a823ad5d30eb7116c77f3fd53d5430acfb2e9ca0df4a7dd97a8baab1dd` |
| `later_records.tsv` | `/tmp/if2453/c89/meta/BBSAG1/` | 55 | `ba2f1dbd71ed4260b0ec24d03d968ba0c6251afbeeb9efb134735ed1064d05f4` |
| `repeated_records.tsv` | `/tmp/if2453/c89/meta/BBSAG1/` | 39 | `d109495ee76e4643c181eae183c452c9a84560f97a79bbe22d6630780368c565` |
| `ani_hits.tsv` | `/tmp/if2453/c89/meta/BBSAG1/` | 279,215 | `1dea48083839ca469fd4f2c0a264fb0ac810c47f533c922fb3589675816cf5cb` |
| `ani_summary.json` | `/tmp/if2453/c89/meta/BBSAG1/` | 1,418 | `d2a54f6c86254347eee6caea5f62161a69dba0f6c1f11fead7d61059fe93f520` |

### 5.3 Mirror to `tower` (`DEC-88 (1)`)

Mirrored sequentially (`scripts/orchestrator/relay_to_tower.sh --dry-run --streams 1` then `--streams 1`, strictly one stream at a time after checking `ssh root@tower pgrep -fl relay_to_` and `ssh root@tower uptime`) to `tower:/mnt/user/research-bulk/bio-box/sources/{PROMGE1,BBSAG1}/`:
- `sources/PROMGE1/SHA256SUMS`: `13` verified lines (`11,398,066,549` B *derived*, `0` failed; *measured*, `tower`, 2026-11-28).
- `sources/BBSAG1/SHA256SUMS`: `16` verified lines (`15,244,892,455` B *derived*, `0` failed; *measured*, `tower`, 2026-11-28).

---

## 6. Summary across both cohorts

| Stage | Genomes | Distinct contigs | Distinct bp | Net-new bp (after 59 / 60 stages) | Depth & mode | r232 / ANI status | Left out / skipped |
|---|---:|---:|---:|---:|---|---|---|
| `PROMGE1` | 2,038,920 | 2,063,640 | 33,663,421,208 | 32,168,295,400 | `called`, `default` (`plasmid` -> `meta`, code `11`) | `9,880` INSDC host accs in R232, `943` `not_in_r232` (none cut) | `3,750,680` `PROMGE2` cellular/nested/non-rep + `12` host archaeal chrs; `1` Prodigal `.ffn` + `4,596` `.prt` skipped; `578` exact repeats |
| `BBSAG1` | 23,078 | 9,738,658 | 39,968,467,080 | 39,968,407,467 | `called`, `default` (`prokaryote` code `11`) | `23,078` `not_visible`; `5,487` `skani >= 96%` hits recorded (none cut) | `14,324` `BBSAG20` LQ genomes left out; `bbsag20_plasmid.fna.gz` + annotations skipped; `0` exact repeats |
| **Total** (*derived*) | **2,061,998** | **11,802,298** | **73,631,888,288** | **72,136,702,867** | | | |

---

## 7. Pod compute & verification gates

All pod work ran at `nice -n 19`, `IOSchedulingClass=idle` inside `/stor/home/if2453/bio-box/images/bio-box-test_faec1f4.sif` (`MemoryMax=28G`, `MemorySwapMax=0`, `CPUQuota=400%`–`800%`), verifying `/tmp/if2453/s2/STOP` absence, `<= 1` other active `c8*` staging unit, and `MV1`/`V3` not in `genecall` before every launch:
- `c89-promge1-census` & `c89-bbsag1-census`: magic-byte and IUPAC nucleotide census across all `28` `PROMGE1` and `90` `BBSAG1` files (`CPUQuota=400%`, *measured*).
- `c89-promge1-stage`: staged `2,038,920` `PROMGE1` genomes (`1.46 GiB` Python RSS, `~164` CPU-min in staging process + `s2_check_inputs.py`, *measured*).
- `c89-bbsag1-stage`: staged `23,078` `BBSAG1` genomes (`10,091,780` `genome_input` rows) with recycled `pyrodigal` / BGZF workers (`max_tasks_per_child=500`, `<= 16.5` CPU-min and `212 MiB` RSS per worker; main Python process `< 40` CPU-min and `2.55 GiB` RSS, `posix_fadvise` cache eviction keeping anonymous RAM `< 1 GiB` during `42 GiB` spooling, *measured* via `ps` and `memory.stat`), finishing right as `MGE2` ended.
- `c89-bbsag1-ani`: pre-sketched `BBSAG1` query genomes once (`skani sketch --separate-sketches -t 8`, `417.1` CPU-s), immediately removed temporary query `.fna` files, and ran `skani dist -t 8 --min-af 50` against `8` reference stages + `r226_archaea` (`1,956.5` CPU-s = `32.6` CPU-min total `skani` CPU, `54` min wall `10:11`–`11:05 CST`, *measured*).
- `stage_est.py conf` on `PROMGE1.conf` and `BBSAG1.conf`: both exit `0` (*measured*).
- `tests/r1/test_c89_promge1_bbsag1_inputs.py`, `tests/r1/test_r232_index.py`, and `tests/r0` on `wilkcomp02`: failed first (`FileNotFoundError` before `c89_promge1_bbsag1_inputs.py` existed), then passed 100% under both normal `PATH` and `PATH=/usr/bin:/bin` (*measured*).
