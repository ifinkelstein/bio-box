# S2C phase 2: fetch of the rest of the remainder (card C-25 step 2; 2026-11-13)

Step (2) of card C-25, for the part of the S2C phase-2 remainder that card C-27 did not fetch: the genomes in
species with at most 20 complete genomes. Card C-27 had already fetched and archived the 4,617 picked in the 243
species over 20 (`2026-11-13_s2c_skani_selection.md` section 6). Every number is labelled *measured* (host,
date, command), *derived* (formula and inputs) or *estimate*. Times are CST. Host: wilkcomp02. Written only under
`/tmp/if2453/c25` (new subdirectories `rest`, `batches`, `bin_fetch`, `genomes`, `archive`, and new files in
`logs`), `/tmp/if2453/ci/C-25-2`, `/tmp/if2453/ci-tmp/C-25-2-*`, and one new file on `/stor`.

## Result in one paragraph

All 8,868 genomes were asked of NCBI Datasets in 4 batches of about 9.7 Gbp, in one chain of detached units from
11:06 to 11:43. **8,860 were served and verified** against NCBI's md5: 38,636,894,022 bp, and every genome's bp
equals GTDB's `genome_size`. 11,217,930,326 bytes are kept as `.fna.gz`. **8 were not served** (0.09%,
35,671,166 bp); the most in one batch was 5 of 2,052 (0.24%), so no batch came near the 5% stop. There was no
refusal, no retry after a network error and no md5 mismatch. They are archived on `/stor` as one deterministic tar:
`provisional/s2c_phase2/s2c_phase2_rest_genomes.tar`, 11,234,908,160 bytes, sha256
`94cc747e86c619f80404941508ca5b8b3cffcf5d7302426f32830bd3db7c19df`. Together with C-27's `s2c_phase2_genomes.tar`,
`/stor` now holds the whole phase-2 remainder less 8 genomes: **13,477 genomes, 56,177,990,736 bp** (*derived*).
For the build, 6 stages of about 9.4 Gbp are proposed below, each holding whole species in GTDB taxonomy order.
No stage inputs were made.

## 1. The list

The input is the ANI selection's remainder, `/tmp/if2453/c27/out_ani/remainder.tsv`. Its sha256 read on the pod
at 11:03 was `48e6072846f6b40e396eea2639db61a2c02528be18e704feafaacc7462b93098`, equal to the C-27 note's value.

The derivation (*measured*, 11:05) ran in unit `c25-rest`, inside the image `/tmp/if2453/images/bio-box_64f1bf0.sif`
with `python -I`. The script `/tmp/if2453/c25/rest/derive_rest.py` (sha256 `62d730eb…`, 25 lines, stdlib only)
keeps the rows whose `ani_pick` is empty, sorted by accession, with the same columns. It asserts that every kept
row has `species_complete_genomes` ≤ 20 and every dropped row has more than 20.

| | Genomes | bp |
|---|---|---|
| remainder.tsv | 13,485 | 56,213,661,902 |
| with an `ani_pick` (C-27, already archived) | 4,617 | 17,541,096,714 |
| **without: this step's list, `rest.tsv`** | **8,868** | **38,672,565,188** |

- Both counts equal the brief's numbers exactly. The bp sums close: 17,541,096,714 + 38,672,565,188 =
  56,213,661,902 (*derived*).
- The 8,868 are in 2,945 species. `rest.tsv` has sha256 `f1b0ca84c6b38f1cd37119fe00afa5fddd5365a55724b2bbf388779967996d2f`.
- An independent check with the host's `awk` gave the same count and the same bp: `8868 38672565188`. It
  counted rows whose `ani_pick` field is empty and summed their `genome_size`.

## 2. How the fetch ran

- **Scripts.** The tree was synced to `/tmp/if2453/ci/C-25-2` by `pod_pytest.sh`. `c27_fetch.py` and `r50_fetch.py`
  were copied unchanged from there to `/tmp/if2453/c25/bin_fetch/`. `c27_fetch.py` was not changed, so source, URL,
  client, pace and verification are exactly C-27's (and phase 1's): the `ncbi_datasets_genome` URL from
  `sources.yaml`, 20 accessions per request, one connection, at most 3 requests per second, the same generic user
  agent, no key and no email. An accession missing from its batch is asked alone once. Files are checked against
  `md5sum.txt` and kept as `.fna.gz`. A per-batch `lock.json` records the FASTA's sha256 and the `.fna.gz`'s
  `gz_sha256`. Stop rules: five refusals in a row, or STOP, or more than 5% of a batch not served.
- **Split.** `c27_fetch.py split rest.tsv /tmp/if2453/c25/batches --batch-gbp 10` (unit `c25-split`) gave n =
  round(38.67) = 4 batches. There is no B00, because the 7 phase-1 genomes are not retried again (DEC-31: once,
  already done in C-27). The result is committed as `2026-11-13_s2c_phase2_fetch/split.json` (sha256 `dc9590d9…`).
- **The launcher: a pod-side copy of `c27_fetch.sh`, changed by `sed`.** The committed launcher has its working
  directory (`/tmp/if2453/c27`) and unit names (`c27-fetch-*`) written into it, and runs the host's `python3`.
  The brief asks for `/tmp/if2453/c25`, unit names starting `c25-`, Python from the image and `CPUQuota=100%`.
  `c27_fetch.sh` is not among the files this step may change, so the copy used was produced on the pod as
  `sed -f c25_fetch.sed c27_fetch.sh > c25_fetch.sh` (in `/tmp/if2453/c25/bin_fetch/`; result sha256
  `0c52c91d…`). The sed script, verbatim:

  ```
  s#/tmp/if2453/c27#/tmp/if2453/c25#g
  s#^BIN=\$W/bin$#BIN=$W/bin_fetch\nPY="singularity exec --cleanenv --bind /tmp/if2453 --bind /stor/scratch/Wilke/if2453/all-the-data:/stor/scratch/Wilke/if2453/all-the-data:ro /tmp/if2453/images/bio-box_64f1bf0.sif python -I"#
  s#python3 -I#$PY#g
  s#c27-fetch-#c25-fetch-#g
  s#c27_fetch\.sh#c25_fetch.sh#g
  s#CPUQuota=200%#CPUQuota=100%#
  ```

  The read-only bind of the data root matters. Without it the image cannot see `/stor`. Checked at 11:05:
  `os.path.isdir('/stor/scratch/Wilke/if2453/all-the-data')` is `False` in the image with only
  `--bind /tmp/if2453`, and `True` with the bind. Without it, the fetch script's STOP check (`r50_fetch.STOP`,
  read before every request batch) would never have seen a STOP file. The launcher also checks STOP on the host
  before each batch. The image's Python has its CA bundle (`/opt/conda/ssl/cert.pem`, 118 certificates).
- **Units.** `c25-fetch-01` to `-04`, one at a time, each launching the next when it ends DONE. Each ran under
  `systemd-run --user`, `nice 19`, `MemoryMax=8G`, `MemorySwapMax=0`, `CPUQuota=100%`. Before every launch the
  launcher checked STOP, that no other `c25-fetch-*` unit was active, and that `/tmp` would keep at least 250 GB
  free after the batch's expected size. Started with `sh /tmp/if2453/c25/bin_fetch/c25_fetch.sh launch B01 chain`.

## 3. Per batch

*Measured*: each batch's `genomes/<B>/DONE` and `fetch_result.json`, and `logs/launch.log`.

| Batch | Asked | Served and verified | Not served | bp verified | Downloaded bytes | Kept `.fna.gz` bytes | FASTA bytes | Requests | Wall s |
|---|---|---|---|---|---|---|---|---|---|
| B01 | 2,365 | 2,363 | 2 | 9,665,591,597 | 2,814,700,801 | 2,812,050,920 | 9,786,872,564 | 121 | 548.2 |
| B02 | 2,268 | 2,267 | 1 | 9,666,709,723 | 2,812,356,800 | 2,809,882,112 | 9,787,945,114 | 115 | 569.9 |
| B03 | 2,183 | 2,183 | 0 | 9,663,315,646 | 2,812,377,580 | 2,810,112,228 | 9,784,504,910 | 110 | 561.2 |
| B04 | 2,052 | 2,047 | 5 | 9,641,277,056 | 2,788,019,362 | 2,785,885,066 | 9,762,144,681 | 108 | 509.2 |
| all | 8,868 | 8,860 | 8 | 38,636,894,022 | 11,227,454,543 | 11,217,930,326 | 39,121,467,269 | 454 | 2,188.5 |

- The "all" row is *derived* (sums).
- Each batch had 0 refusals (`refusal_codes` empty), 0 retries after a network error, 0 md5 refetches, and 0
  genomes asked alone that then came back (`fallback_single` empty). 0 genomes had bp different from GTDB's
  `genome_size` (`bp_vs_gtdb_mismatch` empty). 19,103 FASTA records in all.
- Batch lists and locks (sha256): `batch_B01.tsv` `78d52f4f…`, `B02` `3b5a3f9a…`, `B03` `63c9a4ed…`, `B04`
  `eb3819a3…`; `lock.json` B01 `aa9f9c2b…`, B02 `e53590b3…`, B03 `045fb0e1…`, B04 `6ae75e02…`.
- Rates (*derived*): 243 genomes per minute (8,860 / 2,188.5 s × 60) and 5.13 MB/s downloaded, about C-27's
  240 per minute and 4.83 MB/s. Kept bytes per bp: 0.2903 (11,217,930,326 / 38,636,894,022), against the
  launcher's 0.30 used for the disk check and C-27's 0.2930.
- Kept bytes against what C-27's rate predicts: 38,672,565,188 bp × 0.2930 = 11.33 GB (*estimate*); measured
  11.22 GB, 1% under.

### Not served (8; *measured*, `fetch_result.json` `failed`)

Every one is listed as "not returned by Datasets v2 (alone or in its batch); FTP fallback off". Five of them
were asked alone and drew a response that was not a zip (`single_error File is not a zip file` in
`logs/fetch.log`); the other three drew a zip without their FASTA. None is a phase-1 genome (no `species_rank` 1). Per DEC-31
they were not retried, and they stay out. The list is committed as `2026-11-13_s2c_phase2_fetch/not_served.txt`
(sha256 `18956d4c…`).

| Batch | Accession | Species | bp (GTDB) | Rank in species / complete genomes |
|---|---|---|---|---|
| B01 | GCA_001714745.1 | *Lactobacillus* sp945980025 | 2,159,316 | 8 / 8 |
| B01 | GCA_902502875.2 | *Neorhizobium banfieldiae* | 4,737,957 | 3 / 3 |
| B02 | GCF_003315755.1 | *Clostridium butyricum* | 4,588,510 | 6 / 20 |
| B04 | GCF_039519595.1 | *Klebsiella planticola* | 6,421,810 | 12 / 14 |
| B04 | GCF_041228055.1 | *Xanthomonas hydrangeae* | 5,579,236 | 5 / 6 |
| B04 | GCF_041484315.1 | *Xanthomonas codiaei* | 5,389,326 | 2 / 2 |
| B04 | GCF_902502825.2 | *Neorhizobium flavum* | 4,749,058 | 2 / 2 |
| B04 | GCF_964212645.1 | *Thiodubiliella endoseptemdiera* | 2,045,953 | 2 / 2 |

- GCA_902502875.2 and GCF_902502825.2 sit next to two of phase 1's 7 not served (GCA_902502775.2,
  GCA_902502895.2), and C-27's 26 included runs of neighbouring accessions. Both point to withdrawn or
  suppressed submissions. This is likely, not checked.
- Those species keep their other genomes. Where a species is only 2 genomes (three of the eight), it now keeps
  only its phase-1 genome.

### Disk

*Measured*, `logs/disk.log` (`df -Pk /tmp`, `du -sb genomes`):

- 478,460,804 KiB free on `/tmp` at the start (11:06) and 462,203,004 KiB at the end (11:42). Other cards share
  the disk, and the floor in force was 250 GB.
- `genomes/` held 11,226,671,399 bytes at the end.
- The local tar adds 11,234,908,160 bytes; `archive/` holds 11,235,714,805 bytes with its listing. At 11:51
  `/tmp` had 449,786,336 KiB free (460.6 GB; *measured*, `df -Pk /tmp`).

## 4. The archive

The command, in unit `c25-archive`, in the image, from 11:45 to 11:48:

`python -I c27_fetch.py archive rest/served.tsv archive/s2c_phase2_rest_genomes.tar genomes/B01..B04/lock.json --where-nonempty accession --prefix s2c_phase2_rest`

- *Measured* under `/usr/bin/time -v`: 3:34.8 wall, 94.1 s user and 43.9 s system, peak RSS 67,300 kB.
- `served.tsv` holds the rows of `rest.tsv` that are in a lock: 8,860 rows, sha256 `097890b5…`. It was written
  by `/tmp/if2453/c25/rest/c25_post.py keep` (sha256 `6ff92c50…`, stdlib only, in the image; unit `c25-keep`).
  The archive refuses a kept genome that is in no lock, so the 8 not served had to be left out first.
  `--where-nonempty accession` then keeps every row.
- Each file was checked against its lock's `gz_sha256` before it was written.
- **8,860 genomes, 38,636,894,022 bp, 11,217,930,326 bytes of `.fna.gz`.** The tar has 8,865 members:
  `s2c_phase2_rest/<batch>/<acc>.fna.gz`, the 4 batches' `lock.json`, and `s2c_phase2_rest/kept.tsv`.
- Deterministic, as in C-27: GNU tar, members sorted by name, mode 0644, mtime 0, owner 0 and no owner names.
  The host's `tar --numeric-owner -tvf` lists 8,865 members, all `-rw-r--r-- 0/0` with mtime 0.
- Tar: **11,234,908,160 bytes, sha256 `94cc747e86c619f80404941508ca5b8b3cffcf5d7302426f32830bd3db7c19df`**.
- **Copy to `/stor`** (unit `c25-archive-copy`, host `rsync`, script `bin_fetch/c25_copy.sh`). The script refuses
  to start if STOP, the temporary name or the target exists. It copies to
  `/stor/scratch/Wilke/if2453/all-the-data/provisional/s2c_phase2/.s2c_phase2_rest_genomes.tar.tmp`, reads the
  sha256 back there, compares it with the value above, then renames with `mv -n` to
  `s2c_phase2_rest_genomes.tar`. *Measured*: rsync ran 11:48:54 to 11:50:03 (69 s, 163 MB/s, *derived*), the
  sha256 read back was `94cc747e…`, equal, and the file was renamed at 11:51:00. The directory now holds
  `s2c_phase2_genomes.tar` (5,134,489,600 bytes, unchanged) and `s2c_phase2_rest_genomes.tar` (11,234,908,160
  bytes).
- Nothing else on `/stor` was written, overwritten or deleted. `s2c_phase2_genomes.tar` (C-27) is unchanged.
- Recorded in `sources.yaml` (`selection_s2c_phase2_rest`, `fetched_s2c_phase2_rest`, `stored_s2c_phase2_rest`).

## 5. Proposal: stages of about 10 Gbp for the build (not made)

This is a proposal only. No stage inputs were created. What gets split is the whole phase-2 remainder less the 8
not served: **13,477 genomes, 56,177,990,736 bp** (*derived*: 13,485 − 8; 56,213,661,902 − 35,671,166). They sit
in the two tars: 4,617 in `s2c_phase2_genomes.tar` and 8,860 in `s2c_phase2_rest_genomes.tar`.

**Ordering: by GTDB taxonomy string, then `species_rank`, then accession, with each species kept whole in one
stage.** The number of stages is n = round(total_bp / 10 Gbp) = 6. A species goes to stage
1 + floor(n × bp_before_its_first_genome / total_bp), the rule `c27_fetch.py split` uses for genomes, applied to
species.

Why this order, not accession order (the order of the fetch batches):

- **The representatives-only Pfam of DEC-32 (card B-47) saves more.** A stage's proteins are clustered at NR90
  before its scan, and new members join older clusters only at the next release (B-47's design question). With
  all of a species' strains in one stage, its near-identical proteins fall in that stage's clusters, so only
  the representatives are scanned. In accession order every species' strains spread over all six stages, and
  each stage scans its own representative of the same cluster.
- **Step (3) gets its strain-rich case in one stage:** within-stage contig sharing and the fraction of proteins
  new at NR90.
- **Taxonomic order keeps related species together,** so a stage is closer to a coherent block (several stages
  are one phylum: Bacillota, Pseudomonadota). It changes nothing in what a release holds.
- **The cost is uneven stage sizes,** because a species cannot be split: the largest species adds at most about
  80 Mbp (*E. coli*'s 19 picks; *estimate* from 19 × 4.2 Mbp). Measured spread: 9.23 to 9.50 Gbp.

*Measured* with `c25_post.py stages remainder.tsv not_served.txt stage_proposal.tsv` (unit `c25-stages`, in the
image). The table is committed as `2026-11-13_s2c_phase2_fetch/stage_proposal.tsv` (sha256 `45a1fd00…`). The
column "in species over 20" counts the genomes already archived by C-27.

| Stage | Genomes | bp | Species | In species over 20 | Phyla | First species | Last species |
|---|---|---|---|---|---|---|---|
| S2C-1 | 1,783 | 9,370,307,157 | 579 | 361 | 11 | *Archaeoglobus fulgidus* | *Anoxybacillus_A rupiensis* |
| S2C-2 | 3,102 | 9,371,800,221 | 526 | 1,273 | 1 | *Geobacillus icigianus* | *Hungateiclostridium thermocellum* |
| S2C-3 | 2,470 | 9,502,585,395 | 775 | 646 | 28 | *Aristaeella* sp900322155 | *Bradyrhizobium diazoefficiens* |
| S2C-4 | 1,972 | 9,231,662,797 | 458 | 627 | 1 | *Bradyrhizobium elkanii* | *Edwardsiella piscicida* |
| S2C-5 | 2,032 | 9,392,003,904 | 306 | 988 | 1 | *Edwardsiella tarda* | *Vibrio vulnificus_B* |
| S2C-6 | 2,118 | 9,309,631,262 | 541 | 722 | 5 | *Vibrio zhugei* | *Akkermansia* sp004167605 |
| all | 13,477 | 56,177,990,736 | 3,185 | 4,617 | | | |

- The "all" row is *derived* (sums). Species count distinct GTDB taxonomy strings: 2,945 of at most 20 plus 243
  over 20 is 3,188, and the 3 species whose only remainder genomes were not served (*Xanthomonas codiaei*,
  *Neorhizobium flavum*, *Thiodubiliella endoseptemdiera*) drop out.
- Phylum order follows GTDB's strings: Archaea first, then bacterial phyla alphabetically, so S2C-2, S2C-4 and
  S2C-5 are single phyla.
- One real alternative, for the PI: put the 243 strain-rich species (4,617 genomes, 17.5 Gbp) in two stages of
  their own and the 8,860 others in four. This would isolate the case DEC-32 measures in step (3), at the price
  of two very redundant stages. Not computed.

## Pod use

- **Units**, all `systemd-run --user` at `nice 19` with `CPUQuota=100%` and a memory cap:
  - inside the image: scopes `c25-rest`, `c25-split`, `c25-stages-dry` (the proposal before the fetch, without
    the not-served list), `c25-keep` and `c25-stages`, each 1 CPU for seconds;
  - the fetch chain `c25-fetch-01..04`, 11:06 to 11:43: one process at a time, network-bound, about 0.5 CPU
    (*measured* once: 37.9 CPU-s in the first 73 s of B01, `CPUUsageNSec`);
  - `c25-archive` (3.6 min, 0.64 CPU, in the image) and `c25-archive-copy` (host `rsync` and `sha256sum`).
- **Threads.** At most one CPU of ours at any time. Before the launch the account's CPU, summed from
  `ps -u if2453 -o pcpu=` over processes at least 5 s old, was 34.3 (11:06). It was 32.9 before the archive
  (11:44). The total stayed under 40.
- **Tests:** `pod_pytest.sh C-25-2` (one pytest process, in the test image).
- **Network:** requests from the pod only from `c25-fetch-*`, to the listed Datasets URL. 454 requests in all.
- **Read only:** `/tmp/if2453/c27/out_ani/remainder.tsv`.
- **Units not touched:** `r50-*`, `ci-*`, `ex2*`, `c27-*`, `w27-*`.

## What is left on `/tmp` that the PI could clear

- `/tmp/if2453/c25/genomes/` (11,226,671,399 bytes): the fetched `.fna.gz` files and locks, now in
  `s2c_phase2_rest_genomes.tar` on `/stor`. Keep it while the stages are built from local disk.
- `/tmp/if2453/c25/archive/` (the local copy of the tar and its listing, 11,235,714,805 bytes): can go once the
  `/stor` copy is trusted.
- C-27's (not this step's): `/tmp/if2453/c27/genomes/` (23,840 genomes not kept, about 29.7 GB) and its sketches
  (13.6 GB), as that note says.
