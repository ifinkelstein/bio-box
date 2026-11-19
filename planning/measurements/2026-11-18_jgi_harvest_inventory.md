# 2026-11-18 JGI, NERSC and Zenodo Harvest Inventory (Card C-48a)

Date: 2026-11-19 (harvest executed 2026-11-18T22:24Z to 2026-11-19T11:38Z; metadata bundle and `/stor/work` verification completed 2026-11-19 on `wilkcomp02`).
Closes the deferred JGI acquisition reminder (`DEC-50 (2)`). All numbers follow `AGENTS.md` §6 labels (`*measured*`, `*derived*`, `*published*`).

## 1. Verified `/stor/work` Archive Inventory (11 Archives)

Directory: `/stor/work/Wilke/if2453/all-the-data/jgi_20261118/` on `wilkcomp02`.
Verification command: `scripts/stress/c48a_verify_stor.py --stor-dir /stor/work/Wilke/if2453/all-the-data/jgi_20261118 --sums-tsv /tmp/if2453/s3/dl/jgi/STOR_SHA256SUMS.tsv --sources-yaml sources.yaml --hash-max-bytes 500000000` inside `bio-box-test_faec1f4.sif` (`OMP_NUM_THREADS=1`, 2026-11-19).

- Every archive is `>= 67,108,864` bytes (`64 MiB`, `AGENTS.md` §7.3); smallest archive is `88,697,809` bytes (`84.59 MiB`, *measured*).
- Zero `.cookie.cfg` or temp/partial files exist in `/stor/work/Wilke/if2453/all-the-data/jgi_20261118/`, inside any archive member list, or under `/tmp/if2453/s3/` (*measured*, `find /tmp/if2453/s3 -name '*cookie*'`, 2026-11-19).
- All 10 `jgi-bulk-dl` archives passed `unzip -tq` (for `.zip`) or `tar` creation + readback `sha256sum` after `dd bs=16M conv=fsync` (`bulk.log` / `STOR_SHA256SUMS.tsv`, *measured*, wilkcomp02, 2026-11-18 to 2026-11-19); `jgi_gold_and_img_m_metadata_20261118.tar` was copied via `dd bs=16M conv=fsync` + readback `sha256sum` + `mv -n` (*measured*, wilkcomp02, 2026-11-19).
- Total archived bytes: **`1,137,921,831,192` bytes (`1.138 TB`, `1,059.77 GiB`)** (*derived*: sum of the 11 measured archive byte sizes below).

| Archive | Source entry (`sources.yaml`) | Bytes (*measured*) | Members (*measured*) | Uncompressed member bytes (*measured*) | Archived UTC (*measured*) | SHA-256 (*measured*) |
|---|---|---:|---:|---:|---|---|
| `gvmags_v2_20261118.tar` | `gvmags_v2` | `2,190,807,040` | `3` | `2,190,793,407` | `2026-11-19T10:42:14Z` | `45a3c34ca91e43e3651cd1d2f960bd913b3a5fc3e926efec32f298a29ee0d383` |
| `jgi_gems_52515_fna_20261118.tar` | `jgi_gems_nersc` | `39,502,684,160` | `52,515` | `39,408,621,598` | `2026-11-19T11:38:20Z` | `8230a5bc7a756e9fb88d58cfd98921236634f77025afb1b4ecc0aa708476eb0c` |
| `jgi_gems_meta_prophages_otus_20261118.tar` | `jgi_gems_nersc` | `44,091,832,320` | `10` | `44,091,811,260` | `2026-11-19T11:17:57Z` | `ffc147dbdf956710225e6f9bb0de8b5cd88c3d5290097ce3597da5f4bb0ad06c` |
| `jgi_gold_and_img_m_metadata_20261118.tar` | `jgi_gold_downloads`, `jgi_img_m_tables` | `405,145,600` | `16` | `405,132,113` | `2026-11-19T13:11:22Z` | `d63d7bba593d96d32dddafbbecf2ccf20ca38aeba2fbe07597ada1dfce95a7bc` |
| `jgi_img_pr_20261118.zip` | `jgi_img_pr` | `5,551,950,068` | `4` | `5,693,548,235` | `2026-11-18T23:05:00Z` | `7a51cb6f1d73a070270fa75f0ae8571028c9428dd3ecba2a592e94f82fabd949` |
| `jgi_img_umag_euk_bins_20261118.zip` | `jgi_img_umag` | `38,572,617,388` | `2` | `38,595,932,479` | `2026-11-19T00:21:48Z` | `45d75d286988eab07b83bfbfd377f50ef0cd5b1e3f5e8443954905688491e2ba` |
| `jgi_img_umag_prok_bins_20261118.zip` | `jgi_img_umag` | `408,896,139,104` | `2` | `409,281,208,641` | `2026-11-19T04:50:10Z` | `bbbe2fa9b1c7887d661383ce93229390c2f8fe7f913729c9cadea3426a240bfb` |
| `jgi_img_vr_v4.1_hc_20261118.zip` | `jgi_img_vr` | `24,978,867,388` | `5` | `26,970,047,951` | `2026-11-18T23:15:37Z` | `b887ac945abaf4f8e9422af9414edaf28c983ed444173f4aa595e05cf8b352d9` |
| `jgi_img_vr_v4_full_20261118.zip` | `jgi_img_vr` | `48,814,883,689` | `5` | `54,062,614,689` | `2026-11-18T23:52:46Z` | `a4075d9c28ef31a4796321725e753d87641404cd9ed28f6a57ded37de45bd8a8` |
| `jgi_inovirus_umag_iuig_meta_20261118.zip` | `jgi_inovirus`, `jgi_img_umag`, `jgi_iuig` | `88,697,809` | `7` | `126,354,384` | `2026-11-18T23:16:08Z` | `f8b546e2a50ed9074f3e14bf75852d7b95dd59cdf0ee9b39a259a723e286e393` |
| `jgi_iuig_genomes_20261118.zip` | `jgi_iuig` | `524,828,206,618` | `2` | `525,205,166,382` | `2026-11-19T10:29:56Z` | `397fe93cfc07b7c1187240d796763ec93a4817f20c747d97eb78fdb41e954f7f` |

## 2. Environmental & Genome Metadata Bundle Summary

Bundle: `jgi_gold_and_img_m_metadata_20261118.tar` (`405,145,600` bytes, `16` members, deterministic GNU tar with `uid=0`, `gid=0`, `mtime=2026-11-18T00:00:00Z`, *measured*, wilkcomp02, 2026-11-19).

### 2.1 IMG/M Converted Tables (`scripts/stress/c48a_img_m_tsv.py convert`)

All 5 on-disk `img_m_*_raw.json` exports were converted to deterministic UTF-8 `.tsv.gz` (`mtime=0`, rows sorted by numeric `IMG Genome ID`, leading checkbox column dropped, `#==#` markup and `<font>` HTML tags stripped, PII/email column guard passed with 0 PII columns). Total IMG/M rows converted: **`126,357`** (*derived*: sum of the 5 measured slices below).

| File | Rows (*measured*) | Unique OIDs (*measured*) | Columns (*measured*) | `.tsv.gz` bytes (*measured*) | SHA-256 (*measured*) |
|---|---:|---:|---:|---:|---|
| `img_m_microbiomes_72col.tsv.gz` | `79,524` | `79,524` | `72` | `5,770,667` | `de93af9458889052dcdc4f5a57f55e6239950f336d697d8c738be7325fad9b02` |
| `img_m_bacteria_finished_70col.tsv.gz` | `23,079` | `23,079` | `70` | `2,678,649` | `241a6f9a04280a09828ff6412bc0301c5f2d27112d6859bff891bfefa1fe7084` |
| `img_m_viruses_70col.tsv.gz` | `20,023` | `20,023` | `70` | `1,169,055` | `0e127e63d4cbcb241716dea09a0e7e2c80c00474da037d133ed9303ed4c1c015` |
| `img_m_archaea_70col.tsv.gz` | `3,140` | `3,140` | `70` | `257,042` | `741043f7a5b42e243c74a938eb2fdb4fa3b7260733b80a77c91adc392ae2ff03` |
| `img_m_eukaryota_70col.tsv.gz` | `591` | `591` | `70` | `66,013` | `276e3ed3b2cde4bf998864a2547f35c741185377074cf31c4133102fbb0a6b9e` |

### 2.2 Cross-Check of On-Disk Metadata Row Counts

All counts below are *measured* directly from the on-disk files on `wilkcomp02` (2026-11-19) by streaming XML `<row>` counts (GOLD `.xlsx`) or line counts minus header (TSV/TXT):

| Source | File / Sheet | Data rows (*measured*) | Notes |
|---|---|---:|---|
| `GOLD` v.10 | `goldData.xlsx` — `Study` | `63,817` | `240,297,340` bytes (`SHA256SUMS` verified) |
| `GOLD` v.10 | `goldData.xlsx` — `Biosample` | `247,513` | |
| `GOLD` v.10 | `goldData.xlsx` — `Organism` | `533,108` | |
| `GOLD` v.10 | `goldData.xlsx` — `Sequencing Project` | `641,400` | |
| `GOLD` v.10 | `goldData.xlsx` — `Analysis Project` | `463,257` | |
| `GOLD` v.10 | `goldSraBiomeImg.xlsx` — `SRA biomes with IMG annotations` | `7,126` | `658,230` bytes |
| `GOLD` v.10 | `GOLDs5levelEcosystemClassificationPaths.xlsx` — `site data` | `2,426` | `84,357` bytes |
| `GOLD` v.10 | `goldCVs.xlsx` | `67` sheets | `188,050` bytes (controlled vocabularies) |
| `IMG-UMAG` | `unrestricted_metagenomes_for_bin_tarballs.02062024.txt` | `18,912` | inside `jgi_inovirus_umag_iuig_meta_20261118.zip` (`6,668,197` bytes) |
| `IUIG` | `img_genomes_tarballs_metadata_111023.txt` | `143,797` | inside `jgi_inovirus_umag_iuig_meta_20261118.zip` (`33,340,727` bytes) |
| `GEMs` | `genomes/genome_metadata.tsv` | `52,515` | matches `52,515` `.fna.gz` members in `jgi_gems_52515_fna_20261118.tar` |
| `GEMs` | `otus/otu_representative_genomes.tsv` | `45,599` | species-level OTUs (`576,561` rows in `genome_to_otu.tsv`) |
| `GEMs` | `prophages/virsorter.tsv` | `24,313` | VirSorter prophage coordinates |

## 3. JGI Custom Dataset Portal Inventory (19 Datasets Enumerated 2026-11-18)

On 2026-11-18 the JGI Data Portal (`files.jgi.doe.gov/search/?api_version=2`) listed 19 `Custom_MPI-*` datasets while JGI's tape-restore service was offline (unavailable since 2026-11-09, *measured* 2026-11-18). Under the 2026-11-19 source freeze (`AGENTS.md` §8), no new portal queries were issued by `C-48a`; the table below records the status of the disk-resident sets harvested on 2026-11-18/19 and the tape-purged sets deferred to `C-48e` (with NERSC/Zenodo mirrors used where available):

| Status (2026-11-18/19) | Dataset / Mirror | Action taken | Disposition |
|---|---|---|---|
| **Online (JDP disk)** | `Custom_MPI-IMG_VR` (v4.1 HC + v4.0 full) | Harvested & SHA-256 locked on `/stor/work` (`73.79 GB`, *derived*) | Ready for `C-48c` (`V2`) |
| **Online (JDP disk)** | `Custom_MPI-IMG_PR` | Harvested & SHA-256 locked on `/stor/work` (`5.55 GB`, *measured*) | Ready for `C-48c` (`PL2`) |
| **Online (JDP disk)** | `Custom_MPI-IMG-UMAG` (meta + euk + prok bins) | Harvested & SHA-256 locked on `/stor/work` (`447.56 GB`, *derived*) | Ready for `C-48d` |
| **Online (JDP disk)** | `Custom_MPI-IUIG` (meta + aggregate tarball) | Harvested & SHA-256 locked on `/stor/work` (`524.83 GB`, *measured*) | Ready for `C-48d` |
| **Online (JDP disk)** | `Custom_MPI-Inovirus` | Harvested & SHA-256 locked in `jgi_inovirus_umag_iuig_meta_20261118.zip` | Ready for `C-48c` |
| **Tape-purged on JDP → NERSC mirror** | `GEMs` (`52,515` MAGs, `45,599` OTUs, `24,313` prophages) | Harvested from `portal.nersc.gov/GEM/` & SHA-256 locked (`83.59 GB`, *derived*) | Ready for `C-48d` |
| **Tape-purged on JDP → Zenodo mirror** | `GVMAGs V2` (`18k` giant viruses, doi:10.5281/zenodo.20126600) | Harvested from Zenodo (3 MD5s + tar SHA-256 verified, `2.19 GB`, *measured*) | Ready for `C-48c` (`GV1`) |
| **Tape-purged on JDP (deferred)** | `Custom_MPI-SMC_2025` (`source_data_202604150843.csv` `750 MB`, `smc_gffs.tar.gz` `6.11 TB`, *measured* 2026-11-18) | Held under 2026-11-19 source freeze | Deferred to `C-48e` post-outage |
| **Tape-purged on JDP (deferred)** | `Custom_MPI-IMG_NR` (`img-nr_20240916.tar.gz`, `58.14 GB`, *measured* 2026-11-18) | Held under 2026-11-19 source freeze & `C-45 (3)` protein-only decision | Deferred to `C-48e` |
| **Tape-purged / out-of-scope (10 remaining JDP `Custom_MPI-*` sets)** | Older/specialized or tape-purged JDP custom datasets (enumerated 2026-11-18; portal raw JSON not retained on `/tmp/if2453/s3/dl/jgi/`) | Not fetched; zero external queries issued under 2026-11-19 source freeze | Re-check only if `C-48e` is unheld |
