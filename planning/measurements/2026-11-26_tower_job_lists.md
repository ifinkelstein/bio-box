# Batches 6 and 6b as job lists for `fetch_to_tower.sh` (card R-12)

Written 2026-11-26, 08:10 to 08:40 CST. The 19 lists are in
[`scripts/orchestrator/tower_jobs/`](../../scripts/orchestrator/tower_jobs/README.md). The README there gives each list's
run order and state.

## What was read, and how

- **Copies of `tower`'s lists and manifests.** Card R-11's worker took them at 07:33 to 07:44 CST: `_batch6.tsv`,
  `_batch6b.tsv`, `_batch7_cngb_urls.tsv`, `_batch7_cyverse.tsv`, `_batch3.tsv`, `_batch4.tsv`, the skipped lists,
  `batch6b_ncbi/*/MANIFEST_*.tsv`, `batch6b_ncbi/_summary.json` and `batch6b_euks/C60_TOPAZ/MANIFEST_OSF_FOLDERS.tsv`.
  The runners of batches 6, 6b and 7 were ended at 08:09 CST, before any of them had fetched, so these copies are what
  they would have run.
- **Provider records the Antigravity workers saved on the laptop.** All were read, none was requested again:
  - `~/.gemini/antigravity-cli/brain/9169394d-…/scratch/batch6_probe.json`, the batch-6 worker's record. It holds the
    Zenodo and figshare record of every batch-6 set from those hosts (licence, file sizes, md5), the Edmond record of
    `PALEO_CALC1`, and headers, listings or git trees for `PROMGE2`, `MOPUP_MICRO`, `AMBIRNA`, `C60_SMAG` and
    `C60_EUKCEN`.
  - `~/.gemini/antigravity-cli/brain/2d6b228e-…/scratch/batch3_resolved.json`, the batch-3 record (Y1000, ViroidDB).
- **The batch-6b worker kept no record** of the Zenodo, decodebiome or ETH responses it read (its scratch holds scripts and
  manifests only).
- **`tower`, read-only over `ssh -o BatchMode=yes root@tower`, at about 08:20 CST:**
  `find /mnt/user/research-bulk/bio-box/sources_incoming -name SHA256SUMS -printf '%P\t%s\t%TY-%Tm-%Td %TH:%TM\n'`
  (52 files), and `cat` of the `SHA256SUMS` of `RNA_VIROID/`, `EUK_Y1000/` and `DSV70/`.
- **Processing.** Each job was classified by its file name, the job list's comment line and the saved record. No
  archive was opened and no outside URL was requested. The scripts are in this card's scratch directory, not in the tree:
  `r12/build.py` wrote the lists, `r12/check_all.sh` ran the checks, and `r12/summ6.py`, `probe_lic.py`, `show.py` and
  `entries.py` produced the tallies below.

## The lists, in run order

Every list passes `fetch_to_tower.sh --check --sources sources.yaml`. The command, run from the worktree root on the
laptop for each list:

```
sh scripts/orchestrator/fetch_to_tower.sh --check --sources sources.yaml scripts/orchestrator/tower_jobs/<list>.tsv
```

Labels for the table:
- Jobs: *measured*, the count `--check` prints.
- Bytes: *derived*, summed from the provider records named.
- Bases: *derived*, summed from the manifests' `base_count`, or the worker's `total_bp` for the efetch stages.
- GB gzip: *estimate*, the batch-6b worker's `est_gz_gb` in `_summary.json`.

| # | List | Release | Jobs | Size | Last line of `--check` (exit 0) |
|---|---|---|---|---|---|
| 1 | `01_r26.11.6_b6b_efetch_O4_WEIRD_MT.tsv` | r26.11.6 | 42 | 13,370,311 bases | `check: 42 URLs, all covered by sources.yaml; nothing fetched` |
| 2 | `02_r26.11.6_batch6_rmd1.tsv` | r26.11.6 | 14 | 15,133,291 bytes (`batch6_probe.json`) | `check: 14 URLs, all covered by sources.yaml; nothing fetched` |
| 3 | `03_r26.11.6_b6b_efetch_O4_CRYP.tsv` | r26.11.6 | 8 | 15,273,263 bases | `check: 8 URLs, all covered by sources.yaml; nothing fetched` |
| 4 | `04_r26.11.6_b6b_efetch_O4_PROT_MT.tsv` | r26.11.6 | 42 | 90,476,697 bases | `check: 42 URLs, all covered by sources.yaml; nothing fetched` |
| 5 | `05_r26.11.6_b6b_efetch_O4_PROT_PT.tsv` | r26.11.6 | 22 | 306,525,925 bases | `check: 22 URLs, all covered by sources.yaml; nothing fetched` |
| 6 | `06_later_b6b_ncbi_TICK_RNA.tsv` | later | 68 | 460,261 bases; 0.0001 GB gzip | `check: 68 URLs, all covered by sources.yaml; nothing fetched` |
| 7 | `07_later_b6b_ncbi_AZOAMIC_SYM.tsv` | later | 3 | 1,011,055 bases; 0.0003 GB gzip | `check: 3 URLs, all covered by sources.yaml; nothing fetched` |
| 8 | `08_later_b6b_ncbi_BORG_ONT.tsv` | later | 1 | 2,426,186 bases; 0.0007 GB gzip | `check: 1 URLs, all covered by sources.yaml; nothing fetched` |
| 9 | `09_later_b6b_ncbi_PICJI_SAT.tsv` | later | 6 | 29,450,714 bases; 0.0085 GB gzip | `check: 6 URLs, all covered by sources.yaml; nothing fetched` |
| 10 | `10_later_b6b_efetch_TICK_RNA_NUC.tsv` | later | 20 | 11,605,883 bases | `check: 20 URLs, all covered by sources.yaml; nothing fetched` |
| 11 | `11_later_b6b_ncbi_CODE4_ACTINO1.tsv` | later | 35 | 44,091,220 bases; 0.0128 GB gzip | `check: 35 URLs, all covered by sources.yaml; nothing fetched` |
| 12 | `12_later_batch6b_direct.tsv` | later | 1 | 16,460,717 bytes (`batch3_resolved.json`) | `check: 1 URLs, all covered by sources.yaml; nothing fetched` |
| 13 | `13_later_b6b_ncbi_AMD_CPRV1.tsv` | later | 1 | 292,431,986 bases; 0.0848 GB gzip | `check: 1 URLs, all covered by sources.yaml; nothing fetched` |
| 14 | `14_later_b6b_ncbi_CAVE_AERO1.tsv` | later | 1,457 | 4,739,904,181 bases; 1.3746 GB gzip | `check: 1457 URLs, all covered by sources.yaml; nothing fetched` |
| 15 | `15_later_b6b_ncbi_TERMITE1.tsv` | later | 2,312 | 4,758,152,229 bases; 1.3799 GB gzip | `check: 2312 URLs, all covered by sources.yaml; nothing fetched` |
| 16 | `16_later_b6b_ncbi_WWB_HIC.tsv` | later | 552 | 10,007,995,396 bases; 2.9023 GB gzip | `check: 552 URLs, all covered by sources.yaml; nothing fetched` |
| 17 | `17_later_b6b_ncbi_MARSUPIAL1.tsv` | later | 6,853 | 14,365,306,472 bases; 4.1659 GB gzip | `check: 6853 URLs, all covered by sources.yaml; nothing fetched` |
| 18 | `18_later_rescue.tsv` | later | 2 | 4,731,016,231 bytes for Y1000 (`batch3_resolved.json`) + about 4.44 GB for DSV70 (*estimate*, the batch-4 worker's statement) | `check: 2 URLs, all covered by sources.yaml; nothing fetched` |
| 19 | `19_later_batch6.tsv` | later | 690 | 119,469,480,803 bytes (`batch6_probe.json`; every kept job has a size) | `check: 690 URLs, all covered by sources.yaml; nothing fetched` |

Totals:
- 12,129 jobs (*derived*, the sum of the column).
- 124,232,091,042 bytes where a provider record gives the size (*derived*, lists 2, 12, 18 for Y1000, and 19).
- 34,241,229,700 bases in NCBI assemblies (*derived*). That is 9.93 GB gzip (*estimate*, the worker's figures summed).
- 437,252,079 bases by efetch (*derived*).
- 4.44 GB for DSV70 (*estimate*).

R-11's note gives the NCBI assemblies of 6b as "33,820 Mbp". The manifests and `_summary.json` both sum to 34,241 Mbp, so
R-11's figure is 1.2% lower (*derived*). That is within AGENTS §6's 20%.

**Release.** DEC-79 places `RMD1` and `O4` in `r26.11.6`, and lists 1 to 5 hold them. No decision row places any other
set of batches 6 or 6b, or the two rescues, in a release. Those are marked `later` and come after `r26.11.6`. DEC-79 also
moves "the long-read sets on `tower`" to `r26.11.6` as `LR3`. It does not say whether batch 6's long-read sets belong
to `LR3`, so they stay in list 19. These are `BBMO_HIFI`, `RHIZO_HIFI`, `SOIL_TBP_ONT`, `SCHONBUCH_ONT`, `MLW_CMAG`
and `HAIMA_LR_SAG`.

**Order within a release** is by size, smallest first. The size is the provider's bytes, the worker's gzip estimate for
NCBI assemblies, and bases for efetch stages. Inside lists 2 and 19 the sets are ordered by their bytes, smallest
first. Lists 1 to 5 are small enough that the order among them hardly matters.

**Which entry each URL matched** (the `match` lines of `--check`, tallied by `r12/entries.py`):
- Batch 6 matched its own `tower_*` entries, written by card R-11.
- The ViroidDB retry matched `tower_rna_viroid`, the rescues `tower_dsv70` and `tower_euk_y1000`.
- The NCBI assemblies matched `ncbi_ftp_genome`, the first of several entries with the same NCBI genome form.
- The efetch URLs matched `beaulaurier_aloha_afvgs` (FASTA) and `ncbi_eutils_c59a_lpl1` (GenBank). They are built exactly
  as `url` and `url_2` of R-11's `tower_batch6b_ncbi_o4_org` give them. The script names the first entry in file order
  whose form matches, and those two entries come first with the same forms.
- `sources.yaml` was not changed: no kept URL needed a form that R-11 had not written.

## How the manifests became five-column lines

**NCBI assemblies** (ten 6b stages; MIDAS_FULL left out):
- One line per manifest row: name `batch6b_ncbi/<STAGE>/assemblies`, file `<accession>_genomic.fna.gz`, the manifest's
  https URL, md5 `-`, licence `INSDC`.
- The directory, file name and licence are what `_fetch_batch6b_ncbi.sh` would have used. R-11's entries state them in
  the same form.
- The two manifests of `TERMITE1` and of `WWB_HIC` are joined into one list per stage. No accession repeats (the builder
  checks), and each list's job count equals the stage's `queued` in `_summary.json`.

**Efetch chunks** (four `O4` stages and `TICK_RNA_NUC`):
- `_fetch_batch6b_ncbi.sh` sent each chunk as an HTTP POST with `db`, `rettype`, `retmode` and `id`.
- A GET URL can be written, and these lists use it. The ids are decimal UIDs joined by commas. Every character is in
  the script's URL rule (`[A-Za-z0-9._~:/?#@%&=+,;!*()-]`), and the script sets no length limit. The longest URL is
  2,304 characters (*measured*, `r12/build.py`). Each chunk holds at most 200 ids. NCBI's E-utilities documentation
  asks for POST above about 200 ids; that is from the writer's knowledge, and the page was not fetched for this card.
- Each chunk gives two jobs, as in the old script:
  - `<STAGE>_chunk_<nnnn>.fna`, with `rettype=fasta`;
  - `<STAGE>_chunk_<nnnn>.gbff`, with `rettype=gbwithparts`.
  The GenBank records are the INSDC records of the same nucleotide sequences, with their features. The old script
  fetched both.
- Directories as the old script had them: `batch6b_ncbi/O4_ORG/<STAGE>` and `batch6b_ncbi/TICK_RNA/nuccore`.
- **Differences from the old run:**
  - `fetch_to_tower.sh` does not gzip, so the files land uncompressed, and their names drop the old `.gz`.
  - The `.fna` files get the script's first-byte `>` test. The `.gbff` files get no content test, whereas the old
    script checked for a first `LOCUS` line. An NCBI error page served with HTTP 200 would land as a `.gbff`.
  - No key, tool or email field is sent.
  - The script's default pause of 1 s keeps the requests under E-utilities' three per second without a key.

## Kept, and on whose word

**Named in the question the PI answered at 21:01 CST on 2026-11-25** ("1. y, ..."; R-11 note, table of approval words):
- the ViroidDB retry (list 12);
- `O4_PROT_MT`, `O4_PROT_PT`, `O4_WEIRD_MT`, `O4_CRYP` (lists 1, 3 to 5);
- `MARSUPIAL1`, `TERMITE1`, `AMD_CPRV1`, `CAVE_AERO1` (lists 13 to 15, 17);
- `RMD1`: its 14 metadata tables (list 2). Its genome tarball is held for the PI on licence (below).

The question also named three more:
- `CANDI1`: left out, as a classifier database.
- `NHP-GUT1` and `GBAT-RNA`: the batch-6b worker wrote no manifest for either (not in `_summary.json`), so there is
  nothing to list.

**Kept under the PI's "go ahead with your recommendations for tower" of today (08:08 CST), never asked about by
name.** All of these are plain nucleotide sequence or its metadata:
- Batch 6b stages never in the question:
  - `BORG_ONT` (1 assembly), `AZOAMIC_SYM` (3), `PICJI_SAT` (6), `CODE4_ACTINO1` (35);
  - `TICK_RNA` (68 assemblies, and 1,809 nuccore records in 10 chunks);
  - `WWB_HIC` (552 assemblies; with its Hi-C table from batch 6).
- Batch 6, which no PI word names. The recommendation he approved today was to rebuild it holding only sequence data.
  Its 35 kept sets, beside `RMD1`, are in list 19:
  - `HIBC`, `C60_NUC`, `BTCN_VIR1`, `ECE_METH`, `C60_MICRO`, `WRVIRDB`, `SAT_P4P2`, `BBMO_HIFI`, `WWB_HIC` (the table);
  - `OILV1`, `RHIZO_SSRNA`, `EMCG_MAR`, `PROT_FUNG_SYM1`, `KERM1`, `RPGUTDB`, `PF_VLP1`, `WRGV1`, `ECE_HODOR`;
  - `MFD_PLAS`, `FELMGDB_CULT`, `CCZ1`, `WILDBOAR1`, `WUS_HS1`, `PALEO_CALC1`, `ET2_PARABAS`, `SOIL_TBP_ONT`, `GSPR`;
  - `HAIMA_LR_SAG`, `SCHONBUCH_ONT`, `MLW_CMAG`, `AWIGEN2_AFR`, `C60_BBMO`, `RHIZO_HIFI`, `WRGM1`, `BBSAG20`.

**The rescues** (list 18), refetched on today's recommendation:
- `DSV70_mag.tar.gz` (batch 4, `FAIL_CURL`).
- `y1000p_genome_files.zip` (batch 3, `FAIL_ZIP` after its md5 matched). A comment line above it says it must not
  run before card R-13 is merged.

**Kept, with something the PI may want to see:**
- `SOIL_TBP_ONT`: the record states `cc-by-sa-4.0`. R-11 listed this share-alike term among the licences to confirm.
- `HAIMA_LR_SAG`: the record states `MIT` for data.
- `ET2_PARABAS`: the comment says assemblies *and* transcriptomes, in one `FigShare.tar.gz`.
- `C60_BBMO`: the SAG tarball holds BRAKER3 annotations beside the scaffolds.
- `C60_MICRO`: `rsbl20230398_si_002.zip` is called "supplementary genome data".
- `RMD1`: its 12 `supplementary-table-<n>.xlsx` are kept as the record's metadata tables, though their names do not
  say what they hold. Its two BGC tables are left out.
- `WWB_HIC` (6b): `PRJNA1295499`'s 27 assemblies are whole metagenome assemblies (up to 664,789,256 bases each, *derived*
  from the manifest's `base_count`), not MAGs.
- `AMD_CPRV1`: one assembly of 292,431,986 bases. The worker skipped two older versions.
- Licences recorded per line: as the job list states them; `INSDC` for NCBI lines.

## Left out, by kind

The counts are jobs (*measured*, counted by `r12/show.py` from what the builder set aside). Bytes are *derived* from
`batch6_probe.json` where it gives them. No left-out job is in any list (the builder writes only kept jobs). A search of the job
lines for the names of every left-out and held set and file kind finds none (*measured*): the command below prints nothing
and exits 1.

```
rg -c "^[^#].*(MIDAS_FULL|RDSDB30|C60_TOPAZ|CCOMP_HIFI|C60_RCG|\.db\.zip|README|\.cds\.|gandi|cngb|mge_profile|\.gff|\.gtf|EggNog|antismash|biosynthetic|mobileOG|MacSyFinder|Source_data|DEATHVALLEY|PALEO_HELIUS|AMBIRNA|MOPUP|PROMGE2|C60_SMAG|C60_EUKCEN|HROM|CRBC|RMDv1-genomes|ftp://)" scripts/orchestrator/tower_jobs/
```

Every left-out job is accounted for:
- Batch 6: 766 jobs (*measured*) = 704 kept (lists 2 and 19) + 18 held on licence + 44 left out by kind.
- Batch 6b's direct jobs: 15 = 1 kept + 7 held + 7 left out.

| Kind | Jobs | Bytes | What and why |
|---|---|---|---|
| `MIDAS_FULL` | 50,837 ENA WGS sets | none recorded; 205,550,965,600 bases and about 59.6 GB gzip by the worker (*estimate*) | The whole MiDAS catalogue, not dereplicated. It probably overlaps `MIDAS_LR` (R-11). The PI's word: leave it out. |
| PDB structures | 1 | 587,257,550 | `RDSDB30/Rider_RDSDB30.tar.gz` |
| CDS-only and protein-only sets | 2 | 4,254,722,762 | `C60_RCG`: `RCG_450.cds.zip` and the `.gtf.zip` of those gene models. The scaffolds are on CNCB, not in the list. |
| anvi'o databases | 11 | 48,265,760,662 | `CCOMP_HIFI`: six `.db.zip` and the five summary tables of those databases. |
| models, code, documents | 13 | 295,109 for the 10 with a size | 10 READMEs and file listings (`AMBIRNA`, `MOPUP_MICRO`, `OBELISK_EXT`, `SCHONBUCH_ONT`, `RIDER_RNA` ×2, `RHIZO_HIFI`, `PROMGE2`, `BBSAG20`, `HROM_ORAL`); `MLW_CMAG/Table_S33_readme_table_rev.xlsx`; the MacSyFinder model of `SAT_P4P2`; `CANDI1/gandi-classifier-db.tar.gz` (16.209 GB by the worker, *estimate*; a classifier database; also no licence record) |
| abundance, profile and figure source-data tables | 2 | 627,650,043 | `PROMGE2/mge_profile.tsv.gz`, `KERM1/Source_data.zip` |
| annotations | 12 | 238,878,923 for the 10 with a size | GFFs of `C60_NUC` (2) and `C60_SMAG`; EggNOG of `C60_SMAG`; `C60_MICRO/File_Collection_S14.zip` (GeneMark-ES); `STARSHIP_PLUS/Genome_annotations.zip` and `Starfish_output.zip`; `RMD1`'s antiSMASH summary and biosynthetic-gene table; `MFD_PLAS/pre-processed-mobileOGs.tsv`; `BTCN_VIR1/RdRp_motif_collection.xlsx`; `AMBIRNA/cmsearch.txt` |
| `C60_TOPAZ` | 5 tables + 988 folders | none recorded | The PI's word. The folders' download URLs are found only at run time. The five tables describe only those MAGs, and no provider record of their licence was kept. |
| unclear by name | 5 | 3,073,274,711 | `CAULIMO_EVE/Supplementary_data.tar.gz`, `STARSHIP_PLUS/Starfish_input.zip`, `OBELISK_EXT/Supplement_hsObl_v2.1.zip`, and `RIDER_RNA`'s two `.tar.zst` datasets. Neither the name nor the comment says they hold or describe nucleotide sequence. For the PI. |
| batch 7 | 75,189 rows of `_batch7_cngb_urls.tsv` (`GPB_PHAGE` 104, `EEMC1` 75,085, from R-11); `_batch7_cyverse.tsv` has 0 job rows | the `bp` column sums to 181,337,158,637 bases (R-11) | The PI's word. `EEMC1` is held until it has an https route. |

"Left out" means not in a list. Nothing on `tower` is deleted or changed by this card. The non-sequence files already
landed in batches 1 to 4 stay, as the recommendation says.

## Held for the PI: no record states a licence

These jobs hold sequence or its metadata and would be kept, but no saved provider record states a licence for them. The
empty or missing licence is the reason; their licence column is shown as the job list gives it.

| Set (batch) | Jobs | Bytes (*derived*, record) | Job list's licence | Why held |
|---|---|---|---|---|
| `DEATHVALLEY_SAG` (6) | 1 | 31,356,721 | `open-zenodo-verify` | the Zenodo record 21729218 states an empty licence |
| `PALEO_HELIUS` (6) | 1 | 58,519,131 | `open-zenodo-verify` | the Zenodo record 8328078 states an empty licence |
| `AMBIRNA` (6) | 3 | 6,427,276 | `CC BY 4.0` | GitHub raw on branch `main`; the record is a git tree with no licence |
| `MOPUP_MICRO` (6) | 4 | 68,223,371 | `open-github-verify` | GitHub raw on branch `master`; the record is a git tree with no licence |
| `PROMGE2` (6) | 5 | 34,439,287,617 for 3 of 5 | `CC BY-SA 4.0` | the record holds HTTP headers only; the licence is the worker's label |
| `C60_SMAG` (6) | 2 | none recorded | `CC BY 4.0` | the record is a directory listing |
| `C60_EUKCEN` (6) | 2 | none recorded | `CC BY 4.0` | the record is a directory listing |
| `HROM_ORAL` (6b) | 2 | none recorded (the worker said 2.285 GB, *estimate*) | `cc-by-4.0` | no provider record was kept |
| `CRBC_ISO` (6b) | 4 | none recorded (the worker said 9.161 GB, *estimate*) | `cc-by-4.0` | no provider record was kept. Its md5s show the worker read Zenodo record 14091751. |
| `RMD1` genome tarball (6b) | 1 | none recorded (the worker said 11.572 GB, *estimate*) | `open-academic` | `sunagawalab.ethz.ch` has no record, and `open-academic` is a label. Its tables (list 2) have a cc-by-4.0 record. The skipped list says ENA `PRJEB47249` holds no usable alternative. |

25 jobs in all (*measured*).

## Lists that write into directories already landed on `tower`

A job already in its directory's `SHA256SUMS` is skipped at run time. These lists write into directories that already
hold landed files (*measured*, the `find` and `cat` above at about 08:20 CST):

| List | Directory | Landed there | This list's jobs landed? |
|---|---|---|---|
| 12 | `RNA_VIROID/` | 5 files of batch 3 | no: `viroiddb-2021-06-06.zip` is not in its `SHA256SUMS` |
| 18 | `DSV70/` | 4 files of batch 4 | no: `DSV70_mag.tar.gz` is not in its `SHA256SUMS` |
| 18 | `EUK_Y1000/` | `y1000p_gff3_files.tar.gz` | no: `y1000p_genome_files.zip` is not in its `SHA256SUMS` |

No other list's directory has a `SHA256SUMS` on `tower`. Batch 5's list (`_batch5.tsv`, still waiting behind batch 4)
names none of these directories.

## For whoever runs them

- Run in the order of the prefix, on `tower`, as the README says.
- Do not run list 18 before card R-13 is merged.
- `fetch_to_tower.sh` takes its lock in the destination, and the old `_fetch_batch2.sh` runs of batches 4 and 5 do not
  take it. A list started while those run would share the line with them; it would not collide with them.

## 2026-11-26, card R-14: the licence hold lifted (DEC-80)

DEC-80 (the PI, about 08:41 CST: "all personal use, no worries re license i approve") lifts the hold of the section
"Held for the PI" above. The 25 jobs held there are now
[`scripts/orchestrator/tower_jobs/20_later_licence_released.tsv`](../../scripts/orchestrator/tower_jobs/20_later_licence_released.tsv).
`r12/build14.py` in the card's scratch directory wrote the list from R-12's record of the held jobs. Nothing was fetched, and `tower` was not
read again.

**What the list holds.**
- Each line is copied unchanged from `_batch6.tsv` or `_batch6b.tsv`, the licence field included (`open-zenodo-verify`,
  `open-github-verify`, `open-academic`, or the worker's `CC BY 4.0`, `cc-by-4.0`, `CC BY-SA 4.0`).
- DEC-80 makes no licence. For each set the licence stays as the earlier table says: an empty field in the record, or no
  record kept, so the job list's text is a label only. The sets are `redistributable: false` (DEC-80 (3), DEC-51).
- The `tower_*` entries of `sources.yaml` carry no `redistributable` field. This card's Files do not allow adding one,
  so the `false` is stated here and in the list's header.
- Order, as R-12 ordered its lists: first the `RMD1` genome tarball, which is `r26.11.6` by DEC-79 (a comment line says
  so), then the other nine sets smallest first. The size is the provider's bytes, else the worker's size statement, and
  sets with neither come last (`C60_EUKCEN`, `C60_SMAG`).
- `AMBIRNA` (`main`) and `MOPUP_MICRO` (`master`) each have a comment line saying their URLs name a GitHub branch, not
  a pinned commit.

**Jobs and bytes.**
- 25 jobs (*measured*, `--check`).
- 34,603,814,116 bytes for the 12 jobs with a size in `batch6_probe.json` (*derived*). Of that, `PROMGE2`'s
  `mge100.fna.gz` is 33,912,610,622 bytes (*derived*, the record's `Content-Length`).
- No record of size for 13 jobs: two of `PROMGE2`'s tables, and `C60_SMAG`, `C60_EUKCEN`, `HROM_ORAL`, `CRBC_ISO` and
  `RMD1`. For the last three the workers stated about 2.285, 9.161 and 11.572 GB, 23.018 GB together (*estimate*).

**Sequence rule.** R-12's rule was applied again. The 25 were already the jobs that rule kept: genome, MAG, SAG, contig
and viral sequence, and tables that describe them (lineage, species and genome metadata, MGE clusters and nesting, the
SMAG supplemental tables). None is left out. What R-12 left out of these sets by kind stays out, and so do the other
holds: `MIDAS_FULL`, `C60_TOPAZ`, PDB, CDS-only, anvi'o, models and code, batch 7. That covers `PROMGE2`'s profile
table and README, `C60_SMAG`'s GFF and EggNOG, the READMEs of `AMBIRNA`, `MOPUP_MICRO` and `HROM_ORAL`, and
`AMBIRNA/cmsearch.txt`.

**`sources.yaml`: unchanged.** No `note` or `license` text of these ten entries says the set is held for the PI on licence.
Their notes say what the record states, for example "states the licence as none (an empty field)" or "no provider record
of the licence was kept", and that stays true. Two kinds of text say "Held for the PI by card R-11", and neither is a
licence hold, so the card does not let this change them:
- the `approved` lines of the batch-6 entries: no PI word named batch 6;
- the `MIDAS_FULL` and TOPAZ notes.

**Checks** (run on the laptop from `/tmp/biobox-wt/R-14`; full output in the card's scratch directory):

```
sh scripts/orchestrator/fetch_to_tower.sh --check --sources sources.yaml scripts/orchestrator/tower_jobs/20_later_licence_released.tsv
-> check: 25 URLs, all covered by sources.yaml; nothing fetched
sh scripts/orchestrator/fetch_to_tower.sh --check --sources sources.yaml scripts/orchestrator/tower_jobs/*.tsv
-> check: 12169 URLs, all covered by sources.yaml; nothing fetched
```

All lists together are 21 files and 12,169 jobs (*measured*). That is R-12's 19 lists (12,129), this list (25), and
`2026-11-26_c72_part1.tsv` (15) of card C-72, merged on the integration branch. The script refuses any name and file
named twice across the lists it is given, so no job repeats across the 21.
