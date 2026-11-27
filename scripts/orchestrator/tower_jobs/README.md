# Job lists for `fetch_to_tower.sh` (card R-12)

Each `.tsv` here is a job list in the five-column form of
[`../fetch_to_tower.sh`](../fetch_to_tower.sh): `<name> <file> <url> <md5 or -> <licence>`, tab-separated. They
replace the held batches 6 and 6b on `tower`, whose runners were ended at 08:09 CST on 2026-11-26 before they had fetched
anything, and hold the two rescues of batches 3 and 4. What was kept, what was left out and why, and the `--check` result
of every list are in
[`planning/measurements/2026-11-26_tower_job_lists.md`](../../../planning/measurements/2026-11-26_tower_job_lists.md).

The lists were built from read-only copies of `tower`'s lists and manifests. Nothing was fetched to build them. No list
has been run.

## Order and state

The prefix is the run order: first the release (DEC-79 puts `RMD1` and `O4` in `r26.11.6`; no decision row assigns
the other sets to a release, so they come after it, marked `later`), then the smallest list first. Bytes are the
provider's where its record gives them. For NCBI assemblies only bases are known: the manifests give `base_count`, and
gzip bytes are the worker's estimate. Efetch files land uncompressed.

Labels for the table: jobs are *measured*, counted in each list on 2026-11-26. Bytes are *derived*, summed from the
provider records the workers saved (`batch6_probe.json`, `batch3_resolved.json`). Bases are *derived*, summed from the
manifests' `base_count` or the worker's `total_bp`. The 4.44 GB is an *estimate*, from the batch-4 worker's statement.

| List | Source batch | Jobs | Size | State |
|---|---|---|---|---|
| `01_r26.11.6_b6b_efetch_O4_WEIRD_MT.tsv` | 6b, chunk manifest | 42 | 13,370,311 bases | not run |
| `02_r26.11.6_batch6_rmd1.tsv` | 6, `RMD1` tables | 14 | 15,133,291 bytes | not run |
| `03_r26.11.6_b6b_efetch_O4_CRYP.tsv` | 6b, chunk manifest | 8 | 15,273,263 bases | not run |
| `04_r26.11.6_b6b_efetch_O4_PROT_MT.tsv` | 6b, chunk manifest | 42 | 90,476,697 bases | not run |
| `05_r26.11.6_b6b_efetch_O4_PROT_PT.tsv` | 6b, chunk manifest | 22 | 306,525,925 bases | not run |
| `06_later_b6b_ncbi_TICK_RNA.tsv` | 6b, NCBI manifest | 68 | 460,261 bases | not run |
| `07_later_b6b_ncbi_AZOAMIC_SYM.tsv` | 6b, NCBI manifest | 3 | 1,011,055 bases | not run |
| `08_later_b6b_ncbi_BORG_ONT.tsv` | 6b, NCBI manifest | 1 | 2,426,186 bases | not run |
| `09_later_b6b_ncbi_PICJI_SAT.tsv` | 6b, NCBI manifest | 6 | 29,450,714 bases | not run |
| `10_later_b6b_efetch_TICK_RNA_NUC.tsv` | 6b, chunk manifest | 20 | 11,605,883 bases | not run |
| `11_later_b6b_ncbi_CODE4_ACTINO1.tsv` | 6b, NCBI manifest | 35 | 44,091,220 bases | not run |
| `12_later_batch6b_direct.tsv` | 6b, `_batch6b.tsv` (the ViroidDB retry) | 1 | 16,460,717 bytes | not run |
| `13_later_b6b_ncbi_AMD_CPRV1.tsv` | 6b, NCBI manifest | 1 | 292,431,986 bases | not run |
| `14_later_b6b_ncbi_CAVE_AERO1.tsv` | 6b, NCBI manifest | 1,457 | 4,739,904,181 bases | not run |
| `15_later_b6b_ncbi_TERMITE1.tsv` | 6b, NCBI manifests (2) | 2,312 | 4,758,152,229 bases | not run |
| `16_later_b6b_ncbi_WWB_HIC.tsv` | 6b, NCBI manifests (2) | 552 | 10,007,995,396 bases | not run |
| `17_later_b6b_ncbi_MARSUPIAL1.tsv` | 6b, NCBI manifest | 6,853 | 14,365,306,472 bases | not run |
| `18_later_rescue.tsv` | 4 (`DSV70_mag.tar.gz`), 3 (`y1000p_genome_files.zip`) | 2 | 4,731,016,231 bytes + about 4.44 GB | not run; **not before card R-13 is merged** |
| `19_later_batch6.tsv` | 6, `_batch6.tsv` | 690 | 119,469,480,803 bytes | not run |
| `20_later_licence_released.tsv` | 6 and 6b: the 25 jobs R-12 held on licence, released by DEC-80 (card R-14) | 25 | 34,603,814,116 bytes for 12 of 25 jobs; no record of size for 13 (about 23.0 GB by the workers' statements, *estimate*, for 7 of them) | not run; its `RMD1` line is `r26.11.6`; licence fields unchanged; `redistributable: false` (DEC-80 (3)) |
| `21_later_batch7_GPB_PHAGE.tsv` | 7, `_batch7_cngb_urls.tsv`, `GPB_PHAGE` (CNGBdb CNP0004979), released by DEC-81 (5) (card R-15) | 104 | 6,518,982 bases in 102 assemblies; no byte size recorded | not run; `ftp://` (no https route); licence field `CNGBdb-open` as listed; `redistributable: false` (DEC-80 (3)) |
| `22_later_batch7_EEMC1.tsv` | 7, `_batch7_cngb_urls.tsv`, `EEMC1` (CNGBdb CNP0007106), released by DEC-81 (5) (card R-15) | 75,085 | 181,337,158,637 bases in 75,082 assemblies; no byte size recorded (about 52 GB by the batch-7 worker, *estimate*) | not run; `ftp://`; **75,085 files: run with `--pause 0`, and after lists 01 to 21**; licence field `CNGBdb-CC-BY-4.0` as listed; `redistributable: false` |

Lists 21 and 22 (card R-15) are batch 7, the two CNGBdb sets. `ftp.cngb.org` refuses https (port 443, measured
2026-11-26), so their URLs are `ftp://`, which `fetch_to_tower.sh` takes only where the `sources.yaml` form a URL matches
is itself `ftp://` (DEC-81 (5)); curl is then held to `ftp` (`--proto =ftp --proto-redir =ftp`). `GPB_PHAGE` (104
files) runs first. `EEMC1` is 75,085 files, one FTP transfer each: run it with `--pause 0` (the default pause of 1 s
alone would add about 21 h, *derived*: 75,085 s), and only after the lists of cards R-12 and R-14 (01 to 20) and list 21.
Its 74,999 bins land in 256 directories `batch7_cngb/EEMC1/assemblies/<2 hex>/`, as the batch-7 list laid them out, each
with its own `SHA256SUMS` and `SOURCE_URLS.tsv`. Its `--check` took 2.83 s on the laptop (*measured*) because the two
batch-7 entries stand first in `sources.yaml`; with them last it had not ended after 33 min.

The `EUK_Y1000` line of `18_later_rescue.tsv` must not run before card R-13 is merged. Until then the script removes a
file whose md5 matched when its archive test fails, and that is how batch 3 lost this file. Run the rescue list only after
the merge.

## Running

On `tower`, from a checkout of the tree that has these lists:

```
sh scripts/orchestrator/fetch_to_tower.sh --check --sources sources.yaml scripts/orchestrator/tower_jobs/<list>.tsv
sh scripts/orchestrator/fetch_to_tower.sh --sources sources.yaml scripts/orchestrator/tower_jobs/<list>.tsv
```

### Several connections (`--jobs N`, card R-22)

`--jobs N` (1 to 8, default 1) fetches N files of the lists at a time, each by its own worker process over its own
connection, inside the one run that holds the destination's lock. Use it for a list of many small files from one host,
where the time goes into opening a connection for each file rather than into moving bytes: `EEMC1` (list 22) landed 325
files in its first 124 minutes at one connection, about 23 s a file (*measured*, `tower`, 2026-11-27, card R-22's
"Why"). For a few large files it gains little. A provider's stated connection limit comes first: for CNGB and NCBI ftp
stay at 4 unless the provider says more. `--pause` is taken by each worker after each file it lands, so keep `--pause 0`
for list 22.

```
sh scripts/orchestrator/fetch_to_tower.sh --sources sources.yaml --pause 0 --jobs 4 scripts/orchestrator/tower_jobs/22_later_batch7_EEMC1.tsv
```

Every rule of one connection holds per file. The log keeps one whole line per file, in the order files end, not in list
order, and the run's last line counts `ok`, `skip` and `fail` over all workers. A `STOP` file ends each worker before its
next file. `TERM` or `INT` ends every worker and its transfer, removes their `.partial` files and releases the lock; a
rerun fetches what is missing. `--jobs` above 1 needs bash 4.3 or later and is refused on an older bash.

The lock (`.fetch_to_tower.lock`) names the run, its workers and each open transfer (`pid`, `workers`, `curl.<K>`), and
a new run refuses it while any of them lives. Pids can be reused, so after a run was killed with `SIGKILL` look for
leftover workers and `curl` processes (`ps`) before you clear a lock by hand. curl's own error messages (`-S`) go to
stderr outside the line lock: each is a whole line, but their order against the run's lines is not kept. A run killed
with `SIGKILL` also leaves its temporary directory (`fetch_to_tower.XXXXXX` under `$TMPDIR` or `/tmp`) behind.

A job already in its directory's `SHA256SUMS` is skipped. Two lists write into directories that already hold landed
files: `12_later_batch6b_direct.tsv` (`RNA_VIROID/`) and `18_later_rescue.tsv` (`DSV70/`, `EUK_Y1000/`). None of
their jobs has landed. The lock of `fetch_to_tower.sh` does not see the old `_fetch_batch2.sh` runs of batches 4 and 5.
