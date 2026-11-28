# Card R-26 — Synchronize `sources.yaml`, `planning/CATALOG_GOAL.md`, and `STAGING_WORKER_RULES.md`

- **Date**: `2026-11-28`
- **Card**: `R-26` (tier `M`, `DEC-50`, `DEC-59`, `DEC-80`, `DEC-88`, `DEC-89`, `AGENTS.md` §6, §8)
- **Hosts**: `tower` (`root@tower`, read-only inspection of `SOURCE_URLS.tsv` and `SHA256SUMS` under `/mnt/user/research-bulk/bio-box/`) and `wilkcomp02` (`if2453@wilkcomp02.ccbb.utexas.edu`, inside `/stor/home/if2453/bio-box/images/bio-box-test_faec1f4.sif` under `MemoryMax=8G`, `CPUQuota=200%`, `nice 19`)

## 1. Summary

Card `R-26` reconciles the repository's external-source ledger (`sources.yaml`), living catalog goal (`planning/CATALOG_GOAL.md`), and worker rulebook (`planning/tasks/STAGING_WORKER_RULES.md`) with:

1. All files landed or actively streaming on `tower` across `/mnt/user/research-bulk/bio-box/sources_incoming/` (Stream 1 + Portal 1), `/mnt/user/research-bulk/bio-box/sources_incoming_s2/` (Stream 2, `DEC-88` / `DEC-89`), and `/mnt/user/research-bulk/bio-box/sources/` (mirrored staged inputs), as of `2026-11-28T18:06:34Z` (*measured*, `tower`, `2026-11-28`, `r26_collect_tower.sh`).
2. All `r26.11.6` cohorts staged on `wilkcomp02` and mirrored to `tower:/mnt/user/research-bulk/bio-box/sources/`: `O4`, `GV2`, `HROM1` (`C-86`), `JGIUG1` (`C-87`), and `PROMGE1`, `BBSAG1` (`C-89`).
3. Pre-fetch URL-form entries for Portal 2 (`Biosino NODE`: `MEER` `OEP004067`, `CPR1` `OEZ00021651` / `OED842740`, `GDSV 2.0` `OEP003844`, `DSR-RNA` `OEP002537`) under `DEC-89 (2)` and `DEC-80 (3)`.
4. A deterministic verifier and reconciler CLI (`scripts/orchestrator/r26_sync_sources_ledger.py`) covered by `tests/r1/test_r26_sync_sources_ledger.py`.

## 2. Live `tower` Inventory and `sources.yaml` Reconciliation (`2026-11-28T18:06:34Z`)

Read-only inventory collection on `tower` (`2026-11-28T18:06:34Z`) scanned `1,034` directory records (*measured*, `tower`, `2026-11-28`, `r26_collect_tower.sh` over `sources_incoming/`, `sources_incoming_s2/`, and `sources/`).

- **Pre-R-26 `sources.yaml` baseline**: `232` entries (*measured*, `wilkcomp02`, `2026-11-28`, `r26_sync_sources_ledger.py` on `HEAD:sources.yaml`), `1,061` top-level `http(s)/ftp` URL forms (*measured*), and `109` entries carrying `"nothing yet"` or `"will be tower:..."` (*measured*, `r26_reconcile_helper.py`).
- **Reconciled `tower` entries**: `97` previously `"nothing yet"` / `"will be tower:..."` entries with landed files on `tower` plus `4` previously partial Batch 3 entries whose retries completed (`tower_soda1`, `tower_euk_y1000`, `tower_rna_viroid`, `tower_segc1`) = `101` entries updated with exact landed file counts, byte totals, date ranges, `SHA256SUMS` digests, and actual `sources_incoming/` vs `sources_incoming_s2/` storage roots (*derived*: `97 + 4 = 101`; *measured* via `r26_reconcile_helper.py` + `tower_r26_11_6_ncbi_fna` brace-expansion reconciliation).
- **Post-R-26 `sources.yaml` totals**: `236` entries (*derived*: `232` baseline + `4` Biosino `NODE` entries = `236`; *measured* via `r26_sync_sources_ledger.py verify`) and `1,070` top-level `http(s)/ftp` URL forms (*measured* via `r26_sync_sources_ledger.py verify`; `1,061` baseline forms preserved verbatim + `9` `https://` forms across the `4` Biosino `NODE` entries, *derived*: `1,061 + 9 = 1,070`, plus `4` `ftps://` forms validated by `verify_sources_yaml`).

### Key Stream 1, Stream 2, and Portal 1 Rollups on `tower`

All counts, byte totals, date ranges, and `SHA256SUMS` digests below are *measured* on `tower` as of `2026-11-28T18:06:34Z` from `SOURCE_URLS.tsv` and `SHA256SUMS` (`r26_collect_tower.sh`); expected queue counts are *measured* line counts of `scripts/orchestrator/tower_jobs/*.tsv` (`2026-11-27_c85_tower_fetch_queues.md`, `DEC-89`).

| `sources.yaml` entry | `tower` directory | Landed files | Queued files | Landed bytes (`B`) | Date range (UTC) | Status / `SHA256SUMS` sha256 |
|---|---|---:|---:|---:|---|---|
| `tower_eemc1` | `sources_incoming/batch7_cngb/EEMC1/` | `27,587` | `75,085` | `20,624,910,975` | `2026-11-27..2026-11-28` | Active Stream 1 (`--jobs 4`; top `SHA256SUMS` `f146be5ab678fb9c963ac73afda171d518c41d97c4ed762485a25defa4f42d0a`, `97` subdirs) |
| `tower_gpb_phage` | `sources_incoming/batch7_cngb/GPB_PHAGE/` | `104` | `104` | `1,931,266` | `2026-11-27` | Complete (top `SHA256SUMS` `9a7a7f3c8d8b1aad6aff2dcea05f88db6cc1d3ca46063f9abf6c4c08d2b54374`; `assemblies` `7925e286f162854089590f7f2a97602b153c5195c944e244d4b717baf69b0228`) |
| `tower_r26_11_6_ncbi_fna` (`ASG1_FNA`) | `sources_incoming_s2/ASG1_FNA/assemblies/` | `403` | `403` | `376,522,191` | `2026-11-28` | Complete (`6330ad1173ee0d19aa2027db532dab9e50a31a325d27935ce39447b02db7c3d6`) |
| `tower_r26_11_6_ncbi_fna` (`DSV70_NCBI`) | `sources_incoming_s2/DSV70_NCBI/assemblies/` | `2,937` | `3,354` | `1,966,027,238` | `2026-11-28` | Queue finished (`417` obsolete `.1` versions superseded on NCBI FTP, *derived*: `3,354 - 2,937 = 417`; `70ef5f806ac5481ea64fe24e7a34c71f170b8f018b9434d8f7d920fa7f8caf02`) |
| `tower_r26_11_6_ncbi_fna` (`GTDB_R232_DELTA`) | `sources_incoming_s2/GTDB_R232_DELTA/` | `14,935` | `14,965` | `18,109,778,780` | `2026-11-28` | Queue finished (`9,448` / `9,478` in `gtdb_r232_priority` [`30` suppressed/replaced, *derived*: `9,478 - 9,448 = 30`] across `10` subdirs; `5,487` / `5,487` in `post_r232_complete` across `100` subdirs) |
| `tower_mgnify_nonhuman_reps` | `sources_incoming_s2/MGNIFY_REPS/` | `7,011` | `22,858` | `6,128,979,288` | `2026-11-28` | Active Stream 2 (`--jobs 2`, `611` subdirs populated as of `18:06:34Z`) |
| `jgi_portal_authenticated_r26_12` | `sources_incoming/JGI_LR532/` | `4` | `532` | `761,883,993` | `2026-11-28` | Active Portal 1 (`<= 2` streams, `65e0655359e3e0ae95920ca26d0ae520849415726be1bdd0579e06e93529fd73`) |
| `tower_soda1` | `sources_incoming/SODA1/` | `960` | `960` | `699,551,455` | `2026-11-26` | Complete (`959` assemblies + manifest; top `d47fdc293c03949bdf5927b156a0e5cf8978ed776bee976bdb7f9fbf5742e8bd`) |
| `tower_euk_y1000` | `sources_incoming/EUK_Y1000/` | `2` | `2` | `5,107,548,243` | `2026-11-26..2026-11-27` | Complete (`8a0c663d74d28f268314df599c62f78185529aa5a7e080cfe9c87d001b540616`) |
| `tower_rna_viroid` | `sources_incoming/RNA_VIROID/` | `6` | `6` | `4,648,462,047` | `2026-11-26` | Complete (`339e625eb10f72de87757a1d998a0a0be6235d6a5b819776b389672b8ba7d127`) |
| `tower_segc1` | `sources_incoming/SEGC1/` | `2` | `2` | `37,329,559,604` | `2026-11-26` | Complete (`20fb3d353ec1896cab961b95f78a5f6b988861aaea54528888f7c125458e87d5`) |
| `tower_globdb_r232_tables` | `sources_incoming/GLOBDB_R232_TABLES/` | `5` | `5` | `24,239,347` | `2026-11-26` | Complete (`64effd1fdf0cf55980c8c57a9690ad881d85bd3d2cacbb4a820e85d9fdd9a1e1`) |
| `tower_gtdb_r232_meta` | `sources_incoming/GTDB_R232_META/` | `3` | `3` | `296,063,620` | `2026-11-26` | Complete (`9298fc3aa3d547a0a8b4da13e53710de46ca670f22d87cb7a344fd3f69b58bf9`) |

### Entries Remaining Pre-Fetch (`0` Landed Files on `tower` as of `2026-11-28T18:06:34Z`)

`12` pre-existing entries (*measured* `0` landed files in `r26_tower_inventory.tsv`) plus the `4` newly added Portal 2 Biosino `NODE` entries (`16` total, *derived*: `12 + 4 = 16`) legitimately retain `"nothing yet"` / `"will be tower:..."`:

- Queued behind `MGNIFY_REPS` on Stream 2 (`1`): `tower_globdb_r232_genomes` (`26_r26_11_6_globdb_r232.tsv`, `292,068,313,974` B provider size, *measured* from `globdb_r232_md5sum.txt` / HTTP `Content-Length` in `C-85`).
- Unfetched Batch 6 / Batch 6b manifest items (`11`): `tower_caulimo_eve`, `tower_starship_plus`, `tower_obelisk_ext`, `tower_rdsdb30`, `tower_rider_rna`, `tower_c60_rcg`, `tower_ccomp_hifi`, `tower_c60_topaz`, `tower_batch6b_rescue_candi1`, `tower_midas_full`, `tower_batch6b_euks_c60_topaz_mags`.
- Pre-fetch Portal 2 Biosino `NODE` entries (`4`, `DEC-89 (2)`): `biosino_node_meer_oep004067`, `biosino_node_cpr1_oez00021651`, `biosino_node_gdsv2_oep003844`, `biosino_node_dsrrna_oep002537`.

## 3. `r26.11.6` Staged Cohorts Recorded in `sources.yaml` and `CATALOG_GOAL.md`

`fetched_pod:` and `staged:` provenance blocks were added to all `24` `sources.yaml` entries feeding the `6` `r26.11.6` staged cohorts (`2 + 6 + 3 + 2 + 8 + 3 = 24` entries, *derived*), plus a spreadsheet-only exclusion note on `tower_gv_geve_euglena`. Stage metrics are *measured* on `wilkcomp02` in `2026-11-27_c86_o4_gv2_hrom1_inputs.md`, `2026-11-27_c87_jgi_ug_inputs.md`, and `2026-11-27_c89_promge1_bbsag1_inputs.md`; mirror `SHA256SUMS` digests are *measured* on `tower` (`2026-11-28T18:06:34Z`).

| Stage | Card | `sources.yaml` entries annotated | Distinct genomes | Distinct contigs | Distinct bases (`bp`) | Net-new contigs | Net-new bases (`bp`) | `tower:.../sources/<STAGE>/SHA256SUMS` sha256 |
|---|---|---:|---:|---:|---:|---:|---:|---|
| `O4` | `C-86` | `2` (`tower_o3_fungal`, `tower_batch6b_ncbi_o4_org`) | `10,944` | `10,944` | `544,499,955` | `9,483` | `395,960,372` | `60938e4eadd1f572c3626f369129f6f5e19b86b70afa527972f9d37a8692ab1e` |
| `GV2` | `C-86` | `6` (`tower_gv_polar_pitot`, `tower_gv_polar_2026`, `tower_gv_polar_marine`, `tower_gv_polar_geve`, `tower_gv_gks`, `tower_gv_wwtp`) | `10,284` | `689,523` | `3,155,057,728` | `644,587` | `2,267,601,494` | `25d990a7fcaf6212c48db64ca18113371ccbd9cec358c9eb74b0c68e18aa6e8c` |
| `HROM1` | `C-86` | `3` (`tower_hrom_oral`, `tower_code4_glmc`, `tower_code4_actino1`) | `3,732` | `266,916` | `5,883,730,291` | `265,528` | `5,653,991,918` | `1db893c1efdade3b6f6f22a7ec6ea12e242f08437929b5fc53cc7b6715722d4e` |
| `JGIUG1` | `C-87` | `2` (`jgi_img_umag`, `jgi_iuig`) | `21,524` | `2,890,974` | `34,657,341,931` | `2,875,624` | `34,447,066,468` | `8e4e115560ba372541848feccb5f6cb74b43822a14277836457e8845d0ed6888` |
| `PROMGE1` | `C-89` | `8` (`tower_promge2`, `tower_gspr`, `tower_mfd_plas`, `tower_ece_hodor`, `tower_sat_p4p2`, `tower_ece_meth`, `tower_picji_sat`, `tower_borg_ont`) | `2,038,920` | `2,063,640` | `33,663,421,208` | `2,030,020` | `32,168,295,400` | `3277dc72c16f54a520a5b0daadb598be267bbdb40b670627c26c43852617c761` |
| `BBSAG1` | `C-89` | `3` (`tower_bbsag20`, `tower_haima_lr_sag`, `tower_deathvalley_sag`) | `23,078` | `9,738,658` | `39,968,467,080` | `9,738,622` | `39,968,407,467` | `118307574b1af0bd6af32a48304396cbeddea8eef4dd1436752792c503659977` |
| **Combined `r26.11.6` staged (`6` cohorts)** | `C-86`, `C-87`, `C-89` | **`24`** | **`2,108,482`** (*derived* sum) | **`15,660,655`** (*derived* sum) | **`117,872,518,193`** (*derived* sum, `117.87 Gbp`) | **`15,563,864`** (*derived* sum) | **`114,901,323,119`** (*derived* sum, `114.90 Gbp`) | `6` mirrored stage directories on `tower` |

## 4. `planning/CATALOG_GOAL.md` & `planning/tasks/STAGING_WORKER_RULES.md` Coverage

- **`planning/CATALOG_GOAL.md` §2**: Added `DEC-88` (`2026-11-27`) and `DEC-89` (`2026-11-28`) status block covering `r26.11.5` build queue, `r26.11.6` staged/in-flight cohorts (`O4`, `GV2`, `HROM1`, `JGIUG1`, `PROMGE1`, `BBSAG1`, `C-90`..`C-92`), Stream 2 (`sources_incoming_s2/`), Portal 1 (`JGI_LR532`), and Portal 2 (`Biosino NODE`) + Dryad.
- **`planning/CATALOG_GOAL.md` §3**: Added explicit stage rows for all `20` previously unrepresented `.conf` stage stems (`GRARE1`, `GRARE2`, `GRARE3`, `EGV1`, `ECE2`, `MIR1`, `DPN1`, `SYM2`, `LRP1`, `MIR2`, `MV2`, `MVEUKV`, `KAVHIFI`, `INPH2604`, `ET1`, `LGNRNAP1`, `O4`, `GV2`, `HROM1`, `JGIUG1`) plus `PROMGE1` and `BBSAG1` (`C-89`), bringing §3 coverage to all `56` distinct `.conf` stage stems (`80` `.conf` files: `50` in `stages/` and `30` in `stages/templates/`, *measured*, `wilkcomp02`, `2026-11-28`, `r26_sync_sources_ledger.py verify`).
- **`planning/CATALOG_GOAL.md` §4.1**: Registered every deferred or excluded remainder from `ET1`, `O4`, `GV2`, `HROM1`, `JGIUG1` (`IUIG` + `IMG-UMAG`), `PROMGE1`, `BBSAG1`, and `C-85` (`ASG1`, `DSV70`, `MGnify`, `GlobDB r232`, `GTDB R232`) with exact file paths, counts, sha256s, and pickup instructions.
- **`planning/CATALOG_GOAL.md` §4.3**: Updated status rows for `GlobDB r232`, `GTDB R232`, `MGnify`, `Global Soil Virus Atlas`, `JGI sets on tape`, `RNA viruses`, `VIRE`, `Public long-read ... JGI LR532`, `SAGs`, `Post-R226/R232`, `Biosino NODE`, `DSV70`, `Niche-habitat`, `Organelles`, `EGV1`, `ASG1`, and `Tier-2`.
- **`planning/tasks/STAGING_WORKER_RULES.md`**: Added closing rule requiring `sources.yaml` and `planning/CATALOG_GOAL.md` synchronization checked with `python scripts/orchestrator/r26_sync_sources_ledger.py verify` at every fetch completion and staging-card merge.

## 5. Verification on `wilkcomp02`

All checks were executed on `wilkcomp02` (`2026-11-28`) inside `/stor/home/if2453/bio-box/images/bio-box-test_faec1f4.sif` (`MemoryMax=8G`, `CPUQuota=200%`, `nice 19`, `OMP_NUM_THREADS=2`):

1. `python -I scripts/orchestrator/r26_sync_sources_ledger.py verify` -> `OK: sources=236 url_forms=1070 conf_stages=56 remainders_checked=12` (exit `0`, *measured*).
2. `python -I scripts/orchestrator/r26_sync_sources_ledger.py verify --baseline-sources /tmp/if2453/r26_baseline_sources.yaml --inventory-tsv /tmp/if2453/r26_tower_inventory.tsv` -> `OK: sources=236 url_forms=1070 conf_stages=56 remainders_checked=12` (`0` missing/altered baseline URL forms and `0` stale `"nothing yet"` / `"will be tower:"` entries on landed `tower` directories, exit `0`, *measured*).
3. `/bin/bash scripts/orchestrator/fetch_to_tower.sh --check --sources sources.yaml` across all `28` `.tsv` job manifests in `scripts/orchestrator/tower_jobs/` -> `check: 128492 URLs, all covered by sources.yaml; nothing fetched` (`0` refusals, exit `0`, *measured*).
4. `pod_pytest.sh R-26 tests/r1/test_r26_sync_sources_ledger.py tests/r1/test_r08_fetch_to_tower.py tests/r0` -> `340 passed` (`0` failures, `0` warnings, *measured*).
