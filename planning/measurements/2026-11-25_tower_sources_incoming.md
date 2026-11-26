# Sources staged on `tower` on 2026-11-25: what is there, where it came from, and how it is listed

Card R-08 (audit of 2026-11-25, section 2.1 (3)). Everything below was read from `tower` (the PI's NAS) over `ssh -o BatchMode=yes
root@tower` with read commands only (`ls`, `cat`, `wc`, `find`, `awk` over a manifest, `date`, `ps`), 2026-11-25 20:09 to 20:23 CST.
Nothing was fetched from any provider, hashed, written, deleted, stopped or started on `tower` by this card. The files copied out by
`cat` are the job lists, the logs, the two fetch scripts and every `SOURCE_URLS.tsv` and manifest.

## 1. What is on `tower`

*Measured*, host `tower`, 2026-11-25 20:11 to 20:12 CST, `ls`, `cat SOURCE_URLS.tsv`, `find . -mindepth 2 -type f -printf '%P\t%s\n'`:

- `/mnt/user/research-bulk/bio-box/sources_incoming/` holds 36 source directories (31 with a `SOURCE_URLS.tsv` and a `SHA256SUMS`; five
  with only a `MANIFEST_<bioproject>.tsv`: `ATAC1`, `PLNT1`, `NGOM1`, `GFS1`, `DTOL1`; `SODA1` has both). The audit counted 32 at 19:50.
- **119 files have landed, 138,734,291,679 bytes** (sum of the bytes column of the 31 `SOURCE_URLS.tsv`; 129.2 GiB, *derived*:
  bytes / 2^30). The audit's `du -sh` of 132 GB at 19:50 was taken on disk, not on these sizes.
- Job lists: `_batch1_jobs.tsv` 45 jobs, `_batch2_part1.tsv` 23, `_batch2_part2.tsv` 18, `_batch3.tsv` 94, `_batch4.tsv` 321, `_batch5.tsv`
  189: **690 jobs** (*measured*, `grep -vc '^#'` on the copies). `_batch4.tsv` and `_batch5.tsv` were written at 20:04 and 20:10 CST, after the audit's
  19:50 reading. Part B of batch 4 is `_fetch_batch4_ncbi.sh`, which reads six `MANIFEST_*.tsv` (11,821 assemblies, base_count sum
  42,412,110,240 in all, *derived*: sum of the base_count columns, see section 4) instead of a job list.
- **Landed:** batch 1 (45 files: `EGV1`, `ECE2`, `MIR1`, `GSV1`, `ASG1`, `ET1`), batch 2 (41 files; `_batch2_detached.log` ends `BATCH2 COMPLETE`
  17:41:03 CST, 22 `SKIP` lines from the restart of 16:49 and no `FAIL`), and batch 3 as far as `midas_core_genera_kegg_annotation.tsv`
  (33 files; log last line 2026-11-26T01:41:55Z).
- **Not landed:** the rest of batch 3 (two `MIDAS_LR` files and then `LR_ARC` to `ARC1`, 59 files, so 61 of the 94 jobs of batch 3), all of batch 4 and all of batch 5. At 20:22 CST
  process 2685305 (`_fetch_batch2.sh ... _batch3.tsv`, 2 h 41 min old) was running curl on `midas_core_genera_hqmag.tar.gz`, whose
  `.midas_core_genera_hqmag.tar.gz.partial` held 7,353,286,656 bytes at 20:11 and 9,838,280,704 at 20:22 (*measured*, `find -printf`).
  A `.partial` is not a landed file and is not counted. `_run_batch4_detached.sh` and `_run_batch5_detached.sh` wait for the batch-3 process
  to end (`_batch4_detached.log`: `BATCH4 WAITING FOR BATCH3`).
- No `FAIL_*` line is in either log. Batch 1 has no log (it ran in an interactive session).

## 2. Who asked, and what was approved

The PI's words, with their times, from the orchestrator's prompt history as quoted in the audit's evidence file
(`agy_audit_evidence.md`, section A); each is a request to stage or download, none lists a source or a URL:

| Time (CST, 2026-11-25) | The PI's words, as recorded | What followed |
|---|---|---|
| 15:47 | "stage the nucleotides ... on tower" | batch 1 (45 files) and batch 2 (41 files) |
| 16:51 | "queubatch 3" | batch 3 (94 jobs, chosen by `resolve_batch3.py` against the Zenodo and figshare APIs) |
| 19:47 | "draft a card re renting hardware ... make sure to download everything per the plan above ... launch more workers" | batches 4 and 5 (written 20:04 and 20:10) |

The sources and the URLs were chosen by the orchestrator's workers (`r5-tower-stager`, `r5-tower-batch2`, `r5-tower-batch4`, `r5-giant-viruses`,
`r5-phages-deep`, `r5-prok-portals`, `r5-prok-novel`); the PI approved the staging, not the list (the audit's finding). DEC-69 (2026-11-23) lists
the categories the catalog may grow into and says every new external URL still enters `sources.yaml` before fetch (DEC-50, DEC-51, AGENTS.md
section 8). Every `tower_*` entry in `sources.yaml` carries one of the three quotations under `approved`. Nothing in the record approves
`legacy.tritrypdb.org`, `legacy.microsporidiadb.org`, `legacy.giardiadb.org`, `stacks.stanford.edu`, `ont-open-data.s3.amazonaws.com`,
`sunagawalab.ethz.ch` or `data.mendeley.com` by name; they came with "download everything per the plan".

## 3. `sources.yaml`

83 entries `tower_<directory>` were appended (78 job-list sources and the 5 manifest-only NCBI directories) and `gsv_atlas` gained
`fetched_tower` and `stored_tower` (its archive landed on `tower` as `GSV1`: 960,593,920 bytes, equal to the Content-Length of 2026-11-23, sha256
`8ff20350e9986f2a99efccf4ae7fb5ccd7f18d8899635a4e816fb6b2ad987df0`). `version: 1` is unchanged. Each entry gives the URL form
(`{file}`, `{strain}` one path segment, `{path}` more), the files that landed with bytes, our sha256 and the provider's md5, the licence as the job
list states it, `fetched` with what had and had not landed at 20:12 CST, `stored` and `approved`.

The URL forms are as narrow as the lists allow: a Zenodo form names its record, and every figshare file is listed by its exact URL (209 `url`/`url_N` lines of the
form `https://ndownloader.figshare.com/files/<id>`, one per figshare file of the job lists), as the entries of this file did before this card. No form leaves a figshare file id
open; a figshare id no job list names is refused (tests `test_a_figshare_file_no_job_list_names_is_refused` and `test_no_url_form_of_the_tree_leaves_a_figshare_file_id_open`).
The first version of this card used one `{file_id}` form for all of figshare; the orchestrator's review removed it. The Zenodo forms still allow any file of a listed record
(`{file}`), the VEuPathDB forms any strain and file of a listed host, and `{path}` forms any path under a listed directory.

## 4. Per source

`files` is landed of listed; `bytes` is the sum of the landed files (*measured* on `tower`, `SOURCE_URLS.tsv`, 2026-11-25 20:12 CST); `md5` is the
number of listed files that carry a provider md5 in the job list (the script deletes a file whose md5 differs, and no `FAIL_MD5` is in a log, so every
landed file with an md5 matched it); the licence is the job list's text. Batch 1's list and `SOURCE_URLS.tsv` have no licence column.

| Directory | List | Files | Bytes | md5 | Licence as recorded | State |
|---|---|---|---|---|---|---|
| EGV1 | batch1_jobs | 1 of 1 | 24,753,281 | 1 of 1 | none recorded | landed |
| ECE2 | batch1_jobs | 30 of 30 | 144,363,488 | 29 of 30 | none recorded | landed |
| MIR1 | batch1_jobs | 10 of 10 | 592,796,107 | 10 of 10 | none recorded | landed |
| GSV1 | batch1_jobs | 1 of 1 | 960,593,920 | 0 of 1 | none recorded | landed |
| ASG1 | batch1_jobs | 1 of 1 | 4,961,607,078 | 1 of 1 | none recorded | landed |
| ET1 | batch1_jobs | 2 of 2 | 8,024,914,794 | 2 of 2 | none recorded | landed |
| MGE4 | batch2_part1 | 4 of 4 | 14,325,074 | 4 of 4 | CC BY 4.0 | landed |
| O3_FUNGAL | batch2_part1 | 2 of 2 | 137,723,255 | 2 of 2 | CC BY 4.0 | landed |
| LRP1 | batch2_part1 | 1 of 1 | 249,389,715 | 1 of 1 | CC BY 4.0 | landed |
| DSV1 | batch2_part1, batch3 | 6 of 6 | 480,083,229 | 6 of 6 | CC BY 4.0 | landed |
| MIR2 | batch2_part1 | 5 of 5 | 1,131,772,787 | 5 of 5 | CC BY 4.0 | landed |
| SYM2 | batch2_part1 | 4 of 4 | 3,318,624,234 | 0 of 4 | open-public-NERSC | landed |
| DPN1 | batch2_part1 | 5 of 5 | 9,991,035,723 | 5 of 5 | CC BY 4.0 | landed |
| NX2_EXTREME | batch2_part2 | 5 of 5 | 931,542,531 | 5 of 5 | CC BY 4.0 | landed |
| TSV1 | batch2_part2 | 5 of 5 | 3,464,699,731 | 5 of 5 | CC BY 4.0 | landed |
| SAG2 | batch2_part2 | 1 of 1 | 2,471,331,739 | 1 of 1 | CC BY 4.0 | landed |
| RV1 | batch2_part2 | 4 of 4 | 6,923,237,585 | 4 of 4 | CC BY 4.0 | landed |
| GHM1_VHM1 | batch2_part2 | 3 of 3 | 13,608,323,353 | 3 of 3 | CC BY 4.0 | landed |
| CANDI1 | batch3 | 2 of 2 | 3,329,163 | 2 of 2 | cc-by-4.0 | landed |
| OIL1 | batch3 | 1 of 1 | 2,793,687,427 | 1 of 1 | CC BY 4.0 | landed |
| CFMD1 | batch3 | 2 of 2 | 11,226,893,620 | 2 of 2 | cc-by-4.0 | landed |
| FUTM1 | batch3 | 2 of 2 | 940,829 | 2 of 2 | CC BY 4.0 | landed |
| GXM1 | batch3 | 1 of 1 | 2,224,655,817 | 1 of 1 | CC BY 4.0 | landed |
| FLHC1 | batch3 | 1 of 1 | 32,555,050 | 1 of 1 | cc-by-4.0 | landed |
| MOS1 | batch3 | 1 of 1 | 5,087,895,233 | 1 of 1 | CC BY 4.0 | landed |
| MACA1 | batch3 | 1 of 1 | 4,808,741,144 | 1 of 1 | cc-by-4.0 | landed |
| NX2_CSMC | batch3 | 4 of 4 | 18,121,506,967 | 4 of 4 | CC BY 4.0 | landed |
| WWV1 | batch3 | 9 of 9 | 13,741,179,918 | 9 of 9 | cc-by-4.0 | landed |
| SODA1 | batch3 | 1 of 1 | 38,806 | 1 of 1 | CC BY 4.0 | landed |
| MFD_LR | batch3 | 1 of 1 | 20,792,044,056 | 1 of 1 | cc-by-4.0 | landed |
| MIDAS_LR | batch3 | 3 of 5 | 2,469,706,025 | 5 of 5 | cc-by-4.0 | partly landed |
| LR_ARC | batch3 | 0 of 6 | 0 | 6 of 6 | CC BY 4.0 | queued, nothing landed |
| LR_REVIO | batch3 | 0 of 2 | 0 | 2 of 2 | cc-by-4.0 | queued, nothing landed |
| EUK_Y1000 | batch3 | 0 of 2 | 0 | 2 of 2 | CC0 | queued, nothing landed |
| EUK_STAR | batch3 | 0 of 2 | 0 | 2 of 2 | cc-by-4.0 | queued, nothing landed |
| RNA_LUCA | batch3 | 0 of 30 | 0 | 30 of 30 | CC BY 4.0 | queued, nothing landed |
| RNA_RVMT | batch3 | 0 of 1 | 0 | 1 of 1 | mit-license | queued, nothing landed |
| RNA_VIROID | batch3 | 0 of 6 | 0 | 6 of 6 | other-open, cc-by-4.0 | queued, nothing landed |
| SANDPIPER | batch3 | 0 of 5 | 0 | 5 of 5 | cc-by-4.0 | queued, nothing landed |
| GBR_MGD | batch3 | 0 of 1 | 0 | 1 of 1 | {} | queued, nothing landed |
| ARC1 | batch3 | 0 of 4 | 0 | 4 of 4 | CC BY 4.0 | queued, nothing landed |
| AVGD | batch4 | 0 of 3 | 0 | 3 of 3 | cc-by-4.0 | queued, nothing landed |
| RNA_OBL | batch4 | 0 of 6 | 0 | 6 of 6 | CC BY 4.0 | queued, nothing landed |
| RNA_UPVATLAS | batch4 | 0 of 1 | 0 | 1 of 1 | cc-by-4.0 | queued, nothing landed |
| NX2_CAVE | batch4 | 0 of 8 | 0 | 8 of 8 | CC BY 4.0 | queued, nothing landed |
| DSV70 | batch4 | 0 of 5 | 0 | 0 of 5 | CC0 | queued, nothing landed |
| LR_COMPOST | batch4 | 0 of 4 | 0 | 2 of 4 | ONT-Open-Data | queued, nothing landed |
| EUK_VEUPATH_TRITRYP | batch4 | 0 of 157 | 0 | 0 of 157 | CC BY 4.0 | queued, nothing landed |
| EUK_VEUPATH_MICROSPORIDIA | batch4 | 0 of 106 | 0 | 0 of 106 | CC BY 4.0 | queued, nothing landed |
| EUK_VEUPATH_GIARDIA | batch4 | 0 of 26 | 0 | 0 of 26 | CC BY 4.0 | queued, nothing landed |
| SEGC1 | batch4 | 0 of 2 | 0 | 2 of 2 | cc-by-4.0 | queued, nothing landed |
| OMD2 | batch4 | 0 of 3 | 0 | 3 of 3 | CC BY 4.0 | queued, nothing landed |
| batch5_gv/GV_POLAR_PITOT | batch5 | 0 of 2 | 0 | 2 of 2 | CC BY 4.0 | queued, nothing landed |
| batch5_gv/GV_POLAR_2026 | batch5 | 0 of 1 | 0 | 1 of 1 | CC BY 4.0 | queued, nothing landed |
| batch5_gv/GV_POLAR_MARINE | batch5 | 0 of 98 | 0 | 98 of 98 | CC BY 4.0 | queued, nothing landed |
| batch5_gv/GV_POLAR_GEVE | batch5 | 0 of 13 | 0 | 13 of 13 | CC BY 4.0 | queued, nothing landed |
| batch5_gv/GV_GKS | batch5 | 0 of 1 | 0 | 1 of 1 | cc-by-4.0 | queued, nothing landed |
| batch5_gv/GV_WWTP | batch5 | 0 of 2 | 0 | 2 of 2 | cc-by-4.0 | queued, nothing landed |
| batch5_gv/GV_GEVE_EUGLENA | batch5 | 0 of 1 | 0 | 1 of 1 | cc-by-4.0 | queued, nothing landed |
| batch5_phages/GSVC | batch5 | 0 of 2 | 0 | 2 of 2 | cc-by-4.0 | queued, nothing landed |
| batch5_phages/GHV | batch5 | 0 of 1 | 0 | 1 of 1 | cc-by-4.0 | queued, nothing landed |
| batch5_phages/CGVD | batch5 | 0 of 2 | 0 | 2 of 2 | cc-by-4.0 | queued, nothing landed |
| batch5_phages/TYMEFLIES_V | batch5 | 0 of 1 | 0 | 1 of 1 | CC BY 4.0 | queued, nothing landed |
| batch5_phages/BPHAGE | batch5 | 0 of 2 | 0 | 2 of 2 | cc-by-4.0 | queued, nothing landed |
| batch5_phages/HOVD | batch5 | 0 of 4 | 0 | 4 of 4 | cc-by-4.0 | queued, nothing landed |
| batch5_phages/METHANOGEN_V | batch5 | 0 of 4 | 0 | 0 of 4 | CC BY 4.0 | queued, nothing landed |
| batch5_proks/SYM3 | batch5 | 0 of 1 | 0 | 1 of 1 | cc-by-4.0 | queued, nothing landed |
| batch5_proks/EHI1 | batch5 | 0 of 3 | 0 | 3 of 3 | cc-zero | queued, nothing landed |
| batch5_proks/EMERGE | batch5 | 0 of 2 | 0 | 2 of 2 | cc-by-4.0 | queued, nothing landed |
| batch5_proks/ASG2 | batch5 | 0 of 23 | 0 | 23 of 23 | CC BY 4.0, cc-by-4.0 | queued, nothing landed |
| batch5_proks/DESERT_QAIDAM | batch5 | 0 of 2 | 0 | 2 of 2 | cc-by-4.0 | queued, nothing landed |
| batch5_proks/DESERT_MCMURDO | batch5 | 0 of 3 | 0 | 3 of 3 | cc-by-4.0 | queued, nothing landed |
| batch5_proks/DESERT_MONASH | batch5 | 0 of 1 | 0 | 1 of 1 | CC BY 4.0 | queued, nothing landed |
| batch5_proks/CODE4_GLMC | batch5 | 0 of 3 | 0 | 3 of 3 | cc-by-4.0 | queued, nothing landed |
| batch5_proks/MBGG_TIT | batch5 | 0 of 2 | 0 | 2 of 2 | cc-by-4.0 | queued, nothing landed |
| batch5_proks/ACIDO_SPG | batch5 | 0 of 4 | 0 | 4 of 4 | CC BY 4.0, cc-by-4.0 | queued, nothing landed |
| batch5_proks/SUBSURF_WT29A | batch5 | 0 of 1 | 0 | 1 of 1 | cc-by-4.0 | queued, nothing landed |
| batch5_proks/MGNIFY2 | batch5 | 0 of 10 | 0 | 0 of 10 | CC BY 4.0 | queued, nothing landed |

Manifest-only NCBI directories (batch 4 part B; *measured*, `awk -F'\t' 'NR>1{n++;s+=$5}END{print n, s}'` on each manifest): `ATAC1` 47 assemblies,
37,894,439 bases; `PLNT1` 423, 1,364,663,067; `SODA1` 959, 2,366,049,887; `NGOM1` 1,313, 2,558,619,378; `GFS1` 3,387, 12,011,030,043; `DTOL1` 5,692,
24,073,853,426 (11,821 assemblies, 42,412,110,240 bases in all, *derived*: the sum of the six; compressed sizes are not known). Their licence column is the literal
`INSDC` that `_fetch_batch4_ncbi.sh` writes.

Whether the provider's md5 was checked: for 113 of the 119 landed files, yes (*derived*: 119 less the six below; the md5 is in the list and the script compares it before the
file is kept; no `FAIL_MD5` is in a log). The six landed files without one (`-` in the list) are `ECE2`'s GitHub archive, `GSV1`'s Globus archive and the four `SYM2` NERSC archives.
Queued files with no md5 are `DSV70`, `METHANOGEN_V`, two `LR_COMPOST` files, `MGNIFY2` and the VEuPathDB files; for those the record will be our sha256, the size and the script's archive test.

## 5. Every URL of the job lists is covered by an entry

Command (`scripts/orchestrator/fetch_to_tower.sh --check` reads `sources.yaml` and the lists and fetches nothing), run from the card's worktree on the six
job lists copied from `tower` (*measured*, local, 2026-11-25):

```
bash scripts/orchestrator/fetch_to_tower.sh --check --sources sources.yaml batch1.tsv batch2p1.tsv batch2p2.tsv batch3.tsv batch4.tsv batch5.tsv
check: 690 URLs, all covered by sources.yaml; nothing fetched          # batch3.tsv with the one change below
bash scripts/orchestrator/fetch_to_tower.sh --check --sources sources.yaml ncbi_jobs.tsv    # the 11,821 manifest rows, as <acc>_genomic.fna.gz
check: 11821 URLs, all covered by sources.yaml; nothing fetched
```

The first run on the unchanged `_batch3.tsv` was **refused** (exit 2) at line 103: the job `RNA_VIROID  Benjamin-Lee/viroiddb-2021-06-06.zip` has a `/` in its
file name. `_fetch_batch2.sh` as run would have failed this job (it runs `curl -o ".Benjamin-Lee/viroiddb-2021-06-06.zip.partial"` in a directory with no
`Benjamin-Lee/`, so `FAIL_CURL`; the loop continues). The check above ran on a copy in which the file name is `viroiddb-2021-06-06.zip`; the URL of that
job is covered. A list given to `fetch_to_tower.sh` must not contain such a name (the whole run is refused before any request); the job needs another
file name in the list on `tower`. The test `test_every_source_of_the_2026_11_25_job_lists_is_covered_by_the_trees_sources_yaml` repeats the match
for one real URL per source and per URL form (92 URLs), against the tree's `sources.yaml`; the counts above were re-run on 2026-11-25 after the figshare change and are unchanged (690 and 11,821).

## 6. Differences between the scripts as run and `fetch_to_tower.sh`

- `_fetch_one.sh` (batch 1) sent `-A curl`; `_fetch_batch2.sh` (batches 2 to 5) sent curl's default User-Agent. `fetch_to_tower.sh` sends `-A curl`.
- Both followed redirects to any host; `fetch_to_tower.sh` adds `--proto =https --proto-redir =https` and matches only the URL given.
- The as-run `SKIP` test was `grep -Fq "  <file>" SHA256SUMS` (a substring); it is now an exact name match. `mv -n` that finds the name taken
  now removes the `.partial` (the as-run script left it).
- New: the `sources.yaml` rule; checks of every list field before the first request; a lock directory so two runs cannot write one destination; the `STOP` file
  of AGENTS.md section 7.4; the archive test is skipped, not failed, when `unzip`, `gzip` or `tar` is missing.
- Not archive-tested by either script (no case for the extension): `.gtz`, `.rar`, `.xz`, `.bz2`, `.zst`, `.xlsx`, `.docx`. `NX2_CSMC/non-redundant_genes_annotations.tsv.gtz` (1,733,266,167 bytes,
  landed) has an extension (`.gtz`) that no case matches, so it had only the size and md5 checks, as the list's md5 matched.
- `_run_batch5_detached.sh` is a third wrapper (it waits for the batch-3 and batch-4 processes). Batch 5's stage names are `batch5_gv/<NAME>`, `batch5_phages/<NAME>`
  and `batch5_proks/<NAME>`: directories two levels down; `fetch_to_tower.sh` accepts that form.

## 7. What the batch-5 worker left out (`_batch5_skipped.tsv`, read as is)

Twelve records were skipped or cut down by the worker: a 27.18 GB raw nanopore BAM (`GV_GEVE_EUGLENA`); a Drosophila-virus record that is not the adintovirus one
(`VIR_ADINTO`); an unrelated physics PDF (`VIR_PREPLASMI`); a UMich Deep Blue record that answered HTTP 403 behind Cloudflare (`GV_MYCODNA`, "not bypassed per AGENTS.md
section 8"); ten host MAG tarballs of `GSVC` (8.84 GB); `BPHAGE` (two xlsx files staged, no FASTA in the record); `DESERT_ATACAMA` (SBML models only); `GBR_SED` (embargoed until 2027-02-01);
`MBGG_MIGRATORY` (figshare private link, AWS WAF challenge); `ACIDO_SPG_Z2` (a CSV only); `PARAMO_CPR` (proteins and annotations only); `MGNIFY2` (per-genome FTP tree, only the ten
`genomes-all_metadata.tsv` kept, 108.3 MB). None of these URLs is in `sources.yaml` or a job list.

## 8. Fetched or queued files that are not sequence data or metadata of sequences

For the PI; nothing is deleted. "Landed" is on `tower` now; "queued" is in a job list and not yet fetched. Contents of archives were not opened (nothing was hashed or
unpacked); the judgement is from names and the list's comments.

Landed:

- `MGE4/Samples_name.pdf` (31,942 bytes): a PDF.
- `SODA1/Scripts_for_sequence_analysis.docx` (38,806): a Word document of scripts.
- `CANDI1/mobilome-recovery-workflow.zip` (92,244) and `CANDI1/CAnDi-classifier.zip` (3,236,919): a workflow and a classifier (code, possibly model weights). `CANDI1` has no sequence file.
- `CFMD1/cFMD_humann.tar.gz` (1,418,741,688): HUMAnN functional profiles of food metagenomes (by name), not sequences.
- `TSV1/vOTU-abundance-profile.txt` (64,531,775) and `WWV1/vOTU_table.xlsx` (56,217,852), `WWV1/mOTU_table.xlsx` (2,683,534): abundance tables.
- `MIDAS_LR/midas_core_genera_ko_frequency_table.tsv` (527,672,584) and `midas_core_genera_kegg_annotation.tsv` (1,941,822,954): KEGG-orthology frequency and annotation tables
  (function, not sequence; the annotation table may count as metadata of proteins).
- `NX2_CSMC/tree_files.zip` (9,907), `NX2_CSMC/codes.txt` (5,508), `NX2_CSMC/non-redundant_genes_annotations.tsv.gtz` (1,733,266,167): trees, a code list, gene annotations.
- `MIR1/medvedeva2025_Supplemental_Tables.zip` (171,531,061), `MIR2/PLV_framework_supplementary_data.tar.gz` (265,729,817), `MIR2/EVE_annotations_table.zip` (26,563,179),
  `SYM2/source_data.tar.gz` (254,624,379): supplementary tables and source data of papers.
- `DSV1/Supplementary_Table_1.xlsx`, `Supplementary_Table_2.xlsx`, `Supplementary_table_3.xlsx`, `DSVv1_publicDB_novelty_10422_matched.csv`: supplementary tables (probably vOTU metadata).
- Readmes: `GHM1_VHM1/README.txt`, `MIR1/bellas2023_Readme.txt`, `MIR2/README.txt`.
- `ASG1/Data.zip` (4,961,607,078) and `OIL1/rename_bins.zip` (2,793,687,427): contents unknown (the second may be renamed bins).
- `ECE2/pro-tycheposons-<commit>.tar.gz` (134,806,177): a GitHub source archive (code and data of the repository).

Queued in the job lists (by name): `RNA_LUCA` (30 files, nearly all artefacts of the LucaProt model: **`pytorch_model.bin`** (a model checkpoint), tokenizer and config files, `vocab.txt`, `logs.txt`, a `.DS_Store`, twelve
`*.csv` embedding or label tables, benchmark and result tarballs; the RdRP sequence sets inside `Known_RdRPs.tar.gz` and `00self_sequecing_300aa.fasta.tar.gz` are the likely sequence data); `EUK_STAR/starship-scan.py` (a script) and
`starbase_v0.1.0-pre.sqlite`; `LR_ARC/Workflow-...docx`, `MAG500_METABOLIC_V4.0.result.xlsx`; `AVGD/detecting_archaeal-viruses.sh`; `RNA_VIROID` (`viroiddb` GitHub archive, `analysis`,
`vdsearch` software, `supp`, `data` tarballs, a bowtie2 results table); `RNA_RVMT/RVMT_Zenodo_V4.zip` (a mixed archive); `RNA_OBL` (`Serratus_dir_020924.tar.gz`, `tables.tar.gz`, a supplement table);
`SANDPIPER` (five SRA metadata and summary tables, no sequence); `NX2_CAVE` (six xlsx tables: ARG, BGC, DRAM, CheckM2, classification, metadata); `DSV70` (two tree files);
`LR_COMPOST` (two quality tables); `SEGC1/02.BGC_final.tar.xz`; `OMD2` (a metadata table and a cluster table); `GV_POLAR_PITOT/GVMAGS_Phylogenetic_tree.txt`; `GV_POLAR_GEVE`
(an `.ods` table and two description files); `ASG2` (`Annotations.txt`, `phylogenome_files.rar`, `geNomad-output.zip`); `EHI1` (a readme and a software release zip); `ACIDO_SPG`
(a GitHub archive and tables); `HOVD` (xlsx annotations), `GV_WWTP/README.md`, `CODE4_GLMC` readme and metadata, `MBGG_TIT` xlsx, `MGNIFY2` (ten metadata tables), and 130 VEuPathDB
`.gff` annotation files (`EUK_VEUPATH_*`; annotations of sequences, not sequences). `CANDI1`, `FUTM1` and `SODA1` (the document above) hold no sequence file.

## 9. Open points for the PI

- **No licence recorded** for the six batch-1 sources (`EGV1`, `ECE2`, `MIR1`, `GSV1`, `ASG1`, `ET1`); the `GSV1` archive is the Global Soil Virus Atlas, which `gsv_atlas` records as CC0.
  `GBR_MGD` (one 5,283-MAG tarball, queued) has the licence field `{}` in `_batch3.tsv`: the record states none. `RNA_VIROID` carries `other-open` for one file.
  `LR_COMPOST` carries `ONT-Open-Data` and `SYM2` `open-public-NERSC`, which are access terms, not licences.
- Any licence is the list's text, not checked against the record by this card (AGENTS.md section 8: nothing was fetched).
- The `RNA_VIROID` job `Benjamin-Lee/viroiddb-2021-06-06.zip` in `_batch3.tsv` (line 103) has a `/` in its file name: it will fail in the running batch (`FAIL_CURL`) and needs a new file name in the list on `tower` before `fetch_to_tower.sh` can read that list.
- Nothing on `tower` was deleted, moved or re-hashed; whether the non-sequence files of section 8 should stay is the PI's decision.
- The check in section 5 passes for every queued URL, but a pass says the URL form is listed, not that the PI approved the host (section 2).
