# What the Antigravity sessions staged on `tower` and for `r26.11.4`: an audit (card R-11)

Read on 2026-11-26 from 07:30 to about 08:00 CST. Everything on `tower` was read over `ssh -o BatchMode=yes root@tower`, and
everything on the pod over `ssh if2453@wilkcomp02.ccbb.utexas.edu`. Both hosts were read-only: `ls`, `find -printf`, `cat`/`scp` of job lists,
logs, scripts and manifests, `file`, `od -N8`, `grep`, and on `tower` alone `xz -t` of one file, `tar -tf` of two, all at `nice -n 19`.
Nothing was hashed on either host. On `tower` the `tar -tf` of `NX2_CAVE/MAGs_MQ.tar.zst` wrote its listing to `/tmp/nx2cave_mq.list`, and the same
command removed it. That was the only write on either host, and it touched no data directory.
Local work: `fetch_to_tower.sh --check` and small Python/awk text processing of the copies, in the audit worker's session scratch
directory. There, `tower/` holds the copies and `scripts/` the processing: `compare.py`, `provider.py`, `classify.py`, `confs.py`, `podsizes.py`,
`pseudo_lists.py`, and `gen_sources.py` and `fix_entries.py`, which wrote this card's change to `sources.yaml`. The directory is not kept in
the tree; the copies can be taken again from `tower` with the same `scp`. No outside URL was fetched.

## Verdict

| Batch | Proper data? | In one line |
|---|---|---|
| 1 (45 files) | **Mostly yes** | Every file landed. 43 of 45 have the provider's md5. Licences are only in the stager's report. MIR1's GOEV zip holds Poxviridae and Asfarviridae isolates, and about 36% of it is already in GV1. |
| 2 (41) | **Mostly yes** | Every file landed, and every size we could compare equals the provider's. DSV1's DSVv1 part and parts of TSV1 and WWV1 are tables. GHM1_VHM1 has 0 VHM1 viruses (per the Antigravity scout). |
| 3 (94) | **Partly** | 92 of 94 landed. The Y1000 genome zip was deleted by a test that failed, although its md5 matched. ViroidDB failed. About 30 GB that landed is not genome sequence: LucaProt model files, Sandpiper tables, MiDAS KO/KEGG tables and proteins, ArchaeaHQ proteins, HUMAnN profiles, code. |
| 4 part A (321) | **Partly; still running** | 316 landed. DSV70's MAGs failed. SEGC1's 36 GB MAG archive is partial, and OMD2 is queued. 130 of the files are VEuPathDB GFFs, with no md5 and no licence record. NX2_CAVE's `.tar.zst` files are plain tar. |
| 4 part B (11,821 NCBI assemblies) | **Not fetched yet** | GFS1 contains 14 eukaryotic assemblies, four of them midges. The ATAC1 BioProject is in dispute. |
| 5 (189) | **Not fetched** | Listed by R-08. Licences match the provider records the logs kept for 175 of 189 rows. |
| 6 (766) | **Not fetched; no PI approval by name** | Before R-11 none of its URLs was in `sources.yaml`. Includes unpinned GitHub files, PDB structures, CDS-only sets and records with no licence. |
| 6b (15 + 62,125 assemblies + 134 nuccore chunks + 988 OSF folders) | **Not fetched; approved only in part** | The PI's "1. y" covers the items that were named to him. CRBC_ISO, HROM_ORAL, C60_TOPAZ, MIDAS_FULL and six NCBI stages were added after he answered. |
| 7 (75,189 ftp:// rows) | **Not fetched; held** | EEMC1 follows the PI's pointer; GPB_PHAGE does not. It runs over cleartext FTP with 4 parallel connections, in a list format that `fetch_to_tower.sh` refuses. |
| The 19 stages of `r26.11.4` on the pod | **Yes** | All inputs are present. Each conf's `SUM` sha256 is in its note and each size matches. Every source is listed. None of the contamination or duplication catches applies to these 19. Nothing from `tower` is on the pod. |

## 1. Per batch: planned, landed, failed

The landed list comes from every `SOURCE_URLS.tsv`, `SHA256SUMS` and `find -type f -printf '%P\t%s'` under `sources_incoming/` (measured, `tower`,
2026-11-26 07:33 CST). It is compared with the job lists line by line (`scripts/compare.py`). In every `SOURCE_URLS.tsv`, each file's recorded
size equals its size on disk and each file appears in `SHA256SUMS`. No landed file is missing from the job lists, and none is zero-length.
The only `.partial` file is `SEGC1/.01.MAG.tar.xz.partial`.

| Batch, list (written CST) | Jobs | Planned bytes | Landed | Landed bytes (*measured*) | Failed / partial |
|---|---|---|---|---|---|
| 1 `_batch1_jobs.tsv` (11-25 15:54) | 45 | not recorded | 45 | 14,709,028,668 | none |
| 2 `_batch2_part1/2.tsv` (16:04) | 41 | 39,399,733,817 for the 37 with a provider size (*derived*, `api_metadata.json`, `zenodo_nersc_metadata.json`) | 41 | 42,718,358,051 | none |
| 3 `_batch3.tsv` (16:54) | 94 | 185,314,074,532 (*derived*, `batch3_resolved.json`, all 94) | 92 | 180,566,597,584 | `EUK_Y1000/y1000p_genome_files.zip` FAIL_ZIP; `RNA_VIROID/Benjamin-Lee/viroiddb-2021-06-06.zip` FAIL_CURL (`_batch3_detached.log`) |
| 4A `_batch4.tsv` (20:04) | 321 | about 230.5 GB (*estimate*, the batch-4 worker) | 316 | 47,303,247,503 | `DSV70_mag.tar.gz` FAIL_CURL (curl 28 twice, then 56); `SEGC1/01.MAG.tar.xz` partial at 14,159,605,760 bytes (07:32); `OMD2` 3 files queued |
| 4B six `MANIFEST_*.tsv` (20:05) | 11,821 assemblies | 42,412,110,240 bases (R-08) | 0 | 0 | not started |
| 5 `_batch5.tsv` (20:10) | 189 | 93,033,600,431 for 175 rows (*derived*, `api_cache/`) | 0 | 0 | waiting |
| 6 `_batch6.tsv` (20:54) | 766 | 175,104,030,819 for 282 rows (*derived*, `batch6_probe.json`); the worker said 220.12 GB | 0 | 0 | waiting |
| 6b `_batch6b.tsv` + manifests (21:26) | 15 + 62,125 assemblies (11,288 NCBI + 50,837 ENA WGS) + 134 efetch chunks + 988 OSF folders | the worker said 124 to 129 GB; total_bp 33,820 Mbp NCBI + 205,551 Mbp MiDAS (`_summary.json`) | 0 | 0 | waiting |
| 7 `_batch7_cngb_urls.tsv` (21:37) | 75,189 (GPB_PHAGE 104, EEMC1 75,085) | the `bp` column sums to 181,337,158,637 bases (not bytes); the worker said 52 GB | 0 | 0 | waiting; `_batch7_cyverse.tsv` has 0 rows (Turnstile) |

Total landed: **494 files, 285,297,231,806 bytes** (*derived*, the sum of the four rows above). That is 265.7 GiB. Adding the
14.16 GB partial gives 278.9 GiB, which is the `du -sh` figure of "278 GB" in STATUS.

The queue was running as follows (measured, `ps -eo pid,ppid,etime,args --forest`, 07:33). Batch 4's runner (2887304) runs
`_fetch_batch2.sh _batch4.tsv`, and its curl is writing `01.MAG.tar.xz`. Batches 5 (2893762), 6 (2957811), 6b (3008415) and 7 (3023133) each
wait in a `kill -0` loop, in the order 4, 5, 6, 6b, 7.

**Integrity.**
- 193 of the 494 landed files were compared with the provider's md5 before they were kept, and no `FAIL_MD5` line appears in any log.
  The other 301 have no provider md5: the ECE2 GitHub archive, GSV1, 4 SYM2, 4 DSV70, 2 LR_COMPOST and 289 VEuPathDB files (*derived*).
- Size against the provider's record: all 129 landed files of batches 2 and 3 that have a provider size match it (`scripts/provider.py`).
  The 4 SYM2 files match the NERSC `Content-Length` the batch-2 worker recorded.
- Sample tests on `tower`:
  - `xz -t SEGC1/02.BGC_final.tar.xz`: OK.
  - `tar -tf GSV1/35846_...tar`: 7 entries (the README, 3 sequence/metadata `.gz` files, a CSV).
  - `tar -tf NX2_CAVE/MAGs_MQ.tar.zst`: 838 entries.
  - `od -c -N8` on 9 large zips: each starts `PK\003\004` (`file` reports "data" because they are zip64).
  - `grep -l '<html'` on every file under 2 MB: none.
  - `grep -L '^##gff-version'` on the 130 VEuPathDB GFFs: every one has the header.
- **Defects:**
  - `NX2_CAVE/MAGs_HQ.tar.zst` and `MAGs_MQ.tar.zst` are plain POSIX tar, not zstd (`file`).
  - `NX2_CSMC/non-redundant_genes_annotations.tsv.gtz` starts with `GTZ!`. It is GTZ-compressed, and no landing check covered it.
  - `EUK_Y1000` failed `unzip -tq` after its md5 had matched. In `_fetch_batch2.sh` the md5 test comes before the archive test (lines 45 to 62),
    so a file identical to the provider's was deleted. The likely cause is the `unzip` on `tower` (*inferred*).

## 2. Is it the proper data?

Classified by file name (`scripts/classify.py`; archives were not opened):

| Kind | Files | Bytes (*measured*) |
|---|---|---|
| nucleotide FASTA | 213 | 62,505,213,038 |
| archives (mostly genome/MAG archives; some below) | 68 | 175,782,587,032 |
| gene catalogue | 1 | 16,388,225,385 (`NX2_CSMC/non-redundant_genes.fa.gz`) |
| protein | 5 | 15,373,766,530 (`MIDAS_LR/..._proteins.faa` 14.8 GB, `WWV1/AS_vPCs.fa`, two RdRp `.faa`, a LucaProt structure tarball) |
| tables/text | 64 | 12,049,687,418 (Sandpiper 6.98 GB, MiDAS KO/KEGG 2.47 GB, abundance and supplementary tables) |
| annotations | 132 | 2,720,603,577 (130 VEuPathDB GFFs, the GTZ file, Y1000 GFF3) |
| documents, code, model | 11 | 477,148,826 (`RNA_LUCA/pytorch_model.bin`, tokenizer/config, `.DS_Store`, `starship-scan.py`, `starbase.sqlite`, two `.docx`, a PDF, a `.sh`) |

**Not genome sequence, landed:**
- **Models and code:** RNA_LUCA (6.1 GB of model artefacts and embedding CSVs); CANDI1 (code only); SODA1 (a .docx).
- **Tables:** SANDPIPER (7.0 GB of tables); CFMD1 `cFMD_humann.tar.gz` (1.4 GB of profiles); FUTM1 (2 xlsx, no MAGs); DSV1's DSVv1 part (tables only).
- **Proteins:** ARC1 `faa.zip` (7.6 GB); MIDAS_LR proteins and KO/KEGG tables (17.3 GB).
- **Gene catalogue:** NX2_CSMC.
- **Mixed or other:** RNA_VIROID (analysis and software tarballs and a 1.0 GB bowtie2 table; ViroidDB itself failed); WWV1 (mixed realms and proteins); EUK_STAR (a script and a SQLite database).

**Queued, not genome sequence:**
- Batch 4: OMD2 `SC_G_NR100` (145 GB). The worker calls it genome scaffolds, but the Antigravity audit calls it an ORF gene catalogue.
- Batch 6: `RDSDB30` (RdRp PDB structures); `C60_RCG` (CDS + GTF only); `CCOMP_HIFI` (contigs inside anvi'o `.db.zip`); `WWB_HIC` (one CSV);
  `RMD1` (16 tables); `C60_SMAG` EggNOG; supplementary xlsx in many sets.
- Batch 6b: CANDI1 `gandi-classifier-db.tar.gz` (16.2 GB per the worker, a classifier database).
- Batch 7: GPB_PHAGE metadata.

**Contamination and duplication (the handoff's catches, checked where cheap):**
- **GFS1:** 14 of its 3,387 manifest assemblies are over 20 Mb and eukaryotic, among them `Chironomidae sp.` (25.5 to 50.8 Mb), diatoms and
  `Ochromonas` (*measured* from `batch4_resolved/GFS1_PRJNA781406.tsv`, which is the same size as the copy on `tower`). Filter them before ingestion.
- **WWV1** mixes RNA viruses, prokaryotic MAGs, DNA phages and proteins (file names).
- **MIR1 GOEV** holds 54 Poxviridae and 13 Asfarviridae isolates, and about 36.3% of it overlaps GV1 (the Antigravity scout; not re-measured).
- `uhgv_extended` (92% already in V3) and `jgi_img_pr` (0 new against PL3) are about future sources. V3 and PL3 are the holders, so nothing in `r26.11.4` is affected.
- **Duplications within tower or the tree:**
  - MIDAS_FULL (50,837 ENA WGS sets, 205.6 Gbp) against MIDAS_LR's core-genera HQ MAGs: probably overlapping (*inferred* from the names).
  - EMERGE (batch 5) is the same archive (same md5, `caf267d4...`) as tree entry `emerge_stordalen_v2_0_1`, whose fetch was deferred to LR3.
  - MFD_LR may overlap Microflora Danica/LR entries already in the tree (*inferred*).
  - RMD1 is split across batch 6 (tables) and 6b (genomes). WWB_HIC is split across 6 (CSV) and 6b (NCBI assemblies).
- **ATAC1:** the Antigravity audit says PRJEB60556 is the wrong BioProject, and the batch-6b worker says the alternative, PRJEB57745, has no assemblies. Unresolved.

## 3. `sources.yaml`

**Before R-11** (`fetch_to_tower.sh --check --sources sources.yaml`, the tree at `5d4f57d`, on the copied lists):
- Batch 6: 766 of 766 URLs refused ("not covered").
- Batch 6b: 14 of 15 refused. Only the ViroidDB retry was covered, by `tower_rna_viroid`.
- 6b manifests, one sample per stage and host:
  - The NCBI assemblies are covered, by the generic form of `ncbi_ftp_genome`.
  - MIDAS_FULL's `ftp.ebi.ac.uk/pub/databases/ena/wgs/public/...` URLs are refused.
  - The TOPAZ OSF folder URLs are refused.
  - The efetch chunks are sent as POST. Their GET equivalent is covered by the open `{ids}` efetch forms of other entries; the first match is `beaulaurier_aloha_afvgs`.
- Batch 7: every URL refused as "refused URL" (`ftp://`). With `https://` put in place of `ftp://` they are still "not covered", and the list's six-column layout and file names containing `/` are refused as well.
- Batches 1 to 5: covered, except the known `/` in the name of `_batch3.tsv` line 103.

Two facts govern how much this check means:
- The queued runners call `_fetch_batch2.sh`, `_fetch_batch6b_*.sh` and `_fetch_batch7_cngb.sh`. **None of them reads `sources.yaml`**, so R-08's
  rule does not govern any queued fetch.
- The open NCBI genome and efetch forms admit any assembly and any nuccore id.

**After R-11** (branch `task/R-11`): all 1,627 URLs are covered. The run checked batches 1 to 6b, `_batch3.tsv` with the one file name made
flat, one sample per stage of the 6b manifests, the 988-folder TOPAZ form sample and all 134 efetch chunks
(`check: 1627 URLs, all covered by sources.yaml; nothing fetched`). Batch 7 is not listed and is held for the PI.

**Errors in R-08's entries, corrected on the branch:**
- **Licences of batch 1:** written as "not recorded", although the batch-1 stager's report states one for each source (CC BY 4.0, ECE2 also
  MIT / CC BY-SA 4.0, GSV1 CC0). The entries now cite that report and say that the provider's response was not kept.
- **Record IDs:**
  - `tower_cfmd1` names zenodo:10891046, but its files are from 13285428.
  - `tower_futm1`'s `what` contained CFMD1's comment.
  - `tower_sandpiper` (10547493 against 21287531) and `tower_rna_viroid` (5202945 against 5202946) name one record while their URLs are of a later version
    (the Zenodo API answers with the latest version, *inferred*).
- **Status of 21 entries** (batch-3 rest and batch 4) said "nothing yet". `files`, `fetched` and `stored` are now updated: 378 landed files are listed with their bytes, sha256 and md5. Three of them, all MIDAS_LR, were already listed.
- **Notes added on:**
  - Content and failures: Y1000, DSV70, NX2_CAVE (plain tar), NX2_CSMC (gene catalogue, GTZ), WWV1, GFS1, MIR1 GOEV, ATAC1, OMD2, EUK_STAR,
    RNA_LUCA, CANDI1, SODA1, MGE4, DSV1, FUTM1, GHM1_VHM1.
  - Access terms written as licences: LR_COMPOST ("ONT-Open-Data"), SYM2 ("open-public-NERSC"), VEuPathDB (a worker's claim).
  - O3_FUNGAL: the figshare record says CC BY 4.0; the scout's NC-ND flag is unsourced.
- Every licence of batches 2, 3 and 5 that could be compared equals what the provider record states. Those records were the ones the Antigravity workers saved:
  `api_metadata.json`, `zenodo_nersc_metadata.json`, `batch3_resolved.json`, `api_cache/`, `batch6_probe.json`.

**What cannot be confirmed from the logs:**
- **Licences not confirmed:**
  - Batch 1: the stager's statement only.
  - GBR_MGD: the Zenodo licence field is empty (`{}`).
  - DEATHVALLEY_SAG and PALEO_HELIUS: the Zenodo licence field is empty.
  - VEuPathDB, PROMGE2, C60_SMAG, C60_EUKCEN and HROM_ORAL: a worker's label, no record.
  - RMD1 tarball and CANDI1 db: "open-academic".
  - CNGBdb sets: "CNGBdb-open", "CNGBdb-CC-BY-4.0".
- **No sha of record before fetch:** GitHub raw files on `main` or `master` branches (AMBIRNA, MOPUP_MICRO).
- **Download URLs unknown in advance:** TOPAZ's 988 folders give their download URLs only when listed at run time.

**Approval words (from `~/.gemini/antigravity-cli/history.jsonl` and the orchestrator transcript):**

| Time (CST, 11-25) | The PI's words | Covers |
|---|---|---|
| 15:47 | "agents can also start to stage the nucleotides for next releases, esp on tower..." | batches 1, 2 |
| 16:51 | "queubatch 3" | batch 3 |
| 19:45 / 19:47 (one combined input at 19:47:39 in the transcript) | "make sure to download everything per the plan above" / "launch more workers ... (even those that may require logins)..." | batches 4, 5 (in general words) |
| 20:05, 20:26 | "run even more diversity scouts..."; "launch scouts to also inspect and audit sources ... add those to the wish list" | scouts only. **No PI word names batch 6**, which the orchestrator chained at 20:54 |
| 21:01 | "1. y, 2. decide based on novelty, maybe coral, 3. open each in separate chrome tabs i will log in and tell you" (answering "Queue a small detached Batch 6b...?" with its named items) | 6b items named in the question: RNA_VIROID, CANDI1, O4 organelles, MARSUPIAL1, TERMITE1, NHP-GUT1, AMD-CPRV1, CAVE-AERO1, GBAT-RNA, RMD1. **Not** CRBC_ISO, HROM_ORAL, C60_TOPAZ, BORG_ONT, PICJI_SAT, CODE4_ACTINO1, AZOAMIC_SYM, TICK_RNA, WWB_HIC, MIDAS_FULL |
| 21:13 | "this can be direct, no? https://db.cngb.org/data_resources/project/CNP0007106/" | EEMC1 (batch 7). **Not** GPB_PHAGE |
| 21:24 | "dont bother, for everything that needs a login, make a note for me to do that later" | the logins were dropped |

## 4. The 19 stages of `r26.11.4` on the pod

The wave order comes from the confs' `AFTER=` lines: V2, PL2, PL3, CG1, MGE1, SAT1, ECE1, O1, O2, LPL1, SYM1, LRV1, LGRT1, MGE2, MV1, V3, LR1,
GART1, GART2 (GV1 is deferred). Measured on `wilkcomp02` 2026-11-26 at about 07:50, with `find /tmp/if2453/s2/inputs/<19 stages> -type f -printf '%p\t%s'`,
compared with the confs and notes in the tree (`scripts/confs.py`, `scripts/podsizes.py`):
- All 19 input directories exist. `/tmp/if2453/s2/inputs/` holds 20 directories: these 19 and GV1.
- Every conf has `SUM=` lines: 3 to 12 per stage, 95 in all. Every file they name is on the pod.
- Every `NOTE=` file exists in the tree, and every `SUM` sha256 appears in its note.
- Every size is in the notes, except `genomes.json` and `objects.json` of ECE1 and SYM1, whose sha256 are in their notes. No checksum was recomputed.
- Each source is listed in `sources.yaml`:
  - V2: `jgi_img_vr`, `jgi_inovirus`
  - PL2/PL3: `plasmidscope`, `gut_plasmids_yu2024`
  - CG1: `cgvr`
  - MGE1: `iceberg3`, `isfinder_sequences`, `mgedb_mobileelementfinder`, `plasann_tncentral`, `parnanen_mgedb`
  - SAT1: `inphared_9aug2025`, `ncbi_refseq_viral_r238`, `ncbi_eutils_c55_delta`
  - ECE1: `ncbi_eutils_c65_ece1`, `zenodo_7816421_bgc_plasmids`, `pfeifer_*`
  - O1/O2: `ncbi_refseq_plastid_r238`, `ncbi_refseq_mitochondrion_r238`
  - LPL1: `ncbi_eutils_c59a_lpl1`
  - SYM1: `ncbi_datasets_c58_symbionts`
  - LRV1: `beaulaurier_aloha_afvgs`, `kavagutti_freshwater_hifi_phages`
  - LGRT1: `logan_contigs_baited`
  - MGE2: `jgi_gems_nersc`
  - MV1: `metavr`
  - V3: `vire_genomes`
  - LR1: `ena_longread_metagenome_assemblies`
  - GART1/2: `globdb_r226`
- None of the catches applies to these 19. GOEV/GV1 concerns the deferred GV1. `uhgv_extended`/V3 and `jgi_img_pr`/PL3 show the
  future sources to be redundant, not the stages. GFS1, WWV1 and MIR1 are tower sets that no stage uses.
- Nothing from `tower` reached the pod.
- Minor: V2.conf says "inputs_SUMS beside the note", but `V2/` has no `inputs_SUMS`. Only ECE1, MGE1, MGE2, PL2 and PL3 have one.

## 5. Written by hand on `tower` or the pod; logins

**On `tower`, not in the tree** (copied to `scratchpad/tower/`):
- Fetchers: `_fetch_one.sh`, `_fetch_batch2.sh`, `_fetch_batch4_ncbi.sh`, `_fetch_batch6b_ncbi.sh` (efetch POST, ENA WGS),
  `_fetch_batch6b_topaz.sh` (OSF WaterButler listing, then downloads found at run time), `_fetch_batch7_cngb.sh` (`xargs -P 4`, `ftp://`, md5 where given).
- Runners: `_run_batch{2,3,4,5,6,6b,7_cngb_cyverse}_detached.sh`.
- Job lists, skipped lists, pid files, and the manifests under `batch6b_ncbi/` and `batch6b_euks/C60_TOPAZ/`.
- None of these checks `sources.yaml`, and none honours a `STOP` file. `_fetch_batch2.sh` still uses a substring match to skip files.

**On the pod:** only card working directories. `c62/run_extract.sh`, `w40/run_panel.sh` and others date from 11-21 to 11-24 (`find -newermt`).
No `tower` data is there.

**Logins, tokens and addresses:** I searched every Antigravity `*.jsonl` transcript under `brain/` (`rg -i` for `mailto=`, `email=`,
`api_key=`, `Authorization: Bearer`, `--cookie`, `token=` and e-mail domains) and found none of these in any request. The only cookies are
`Set-Cookie` lines in Dryad's 401 responses. The record of logins:
- At 21:04 the PI reported logging into Zenodo in his own Chrome.
- At 21:11 the orchestrator launched a CDP-controlled Chrome with login tabs.
- No login is recorded in that browser. The PI dropped logins at 21:24.
- No JGI session token was pasted.

The Antigravity resolvers sent many metadata requests to hosts that `sources.yaml` does not list, before or without listing them: the Zenodo, figshare, OSF and Mendeley APIs,
`ftp.cngb.org`, `ngdc.cncb.ac.cn`, `data.cyverse.org`, Dryad, `unblockme.cyverse.org`. That is outside AGENTS §8 as written.

## Decisions that are the PI's

1. **The queue:** whether batches 5, 6, 6b and 7 run as queued. They fetch with scripts that ignore `sources.yaml` and the `STOP` file, so the
   only way to stop them is to end their processes. Two items within batch 4 also need a word:
   - OMD2: 145 GB. Is it a gene catalogue or scaffolds?
   - SEGC1: its 36 GB `01.MAG.tar.xz` is downloading now.
2. **Batch 6:** approve it or drop it. No PI word names it.
3. **Batch 6b:** approve or drop the items that were not in the question of 20:58: CRBC_ISO, HROM_ORAL, C60_TOPAZ (download URLs found at run time),
   BORG_ONT, PICJI_SAT, CODE4_ACTINO1, AZOAMIC_SYM, TICK_RNA, WWB_HIC, and MIDAS_FULL (50,837 sets, 205.6 Gbp, overlapping MIDAS_LR).
4. **Batch 7:** approve it or hold it. EEMC1 needs a decision on cleartext `ftp://` with 4 parallel connections. GPB_PHAGE was never named.
   Not listed in `sources.yaml`.
5. **Licences to confirm:** GBR_MGD, DEATHVALLEY_SAG, PALEO_HELIUS, the VEuPathDB sets, PROMGE2 (BY-SA), SOIL_TBP_ONT (BY-SA), C60_SMAG,
   C60_EUKCEN, HROM_ORAL, the RMD1 tarball, the CANDI1 db, LR_COMPOST, SYM2, the CNGBdb sets, the six batch-1 sources, O3_FUNGAL's NC-ND flag,
   and the MIT-licensed data of RNA_RVMT and HAIMA_LR_SAG.
6. **Non-sequence files:** keep them as sidecars or drop them. RNA_LUCA model, SANDPIPER, CFMD1 humann, the MIDAS_LR tables and proteins,
   ARC1 faa, NX2_CSMC, the CANDI1 and EUK_STAR code, RDSDB30, C60_RCG.
7. **Rescues:**
   - Y1000 genomes: refetch and test with another unzip.
   - DSV70 MAGs: refetch.
   - ViroidDB: already queued in 6b.
8. **Before ingestion:** decontaminate GFS1, WWV1 and MIR1 GOEV; treat NX2_CAVE as tar; settle ATAC1's BioProject.
9. **The fetcher:** whether the queued fetches should be moved onto `scripts/orchestrator/fetch_to_tower.sh`. That script checks `sources.yaml`,
   honours `STOP` and allows https only.
