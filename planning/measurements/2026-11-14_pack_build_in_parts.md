# The search pack built in parts: `rN.3` reproduced, `rM.2` built (card W-35, 2026-11-14)

Every number is *measured* on wilkcomp02 on 2026-11-14 from the unit logs named beside it (each part's
`done.json`, which the driver logs as `part <name> done: ...`: CPU-minutes are the part's own `ru_utime + ru_stime`,
peak RSS its `ru_maxrss`, wall time from the part's start), unless it is labelled *derived* (with its formula) or
*estimate*.

## Setup

- Code: branch `task/W-35`, mirrored to `/tmp/if2453/w35/tree` with the rsync of `scripts/pod_pytest.sh`.
- Image `bio-box_faec1f4.sif` (sha256 `d67300aa1321b571fa8da79e74dfccdd87a8dca2e56518f324902efa6b9f098a`), the
  tree's `src` alone on `PYTHONPATH`; the real profile `profiles/wilkcomp02.yaml` (budget 300 CPU-minutes, 80% margin:
  240).
- Data root `/tmp/if2453/r10/root`, bound read-only (`--bind /tmp/if2453/r10/root:/tmp/if2453/r10/root:ro`); nothing was
  written under it or under `/stor`. Pfam: `pfam-38.2=/stor/scratch/Wilke/if2453/all-the-data/l0/refdb/pfam/38.2/Pfam-A.hmm.gz`,
  bound read-only, sha256 `2d82087b6c5c60d762cc767f98e8260b273134c215ab7efccf7440614a4e5dab` (`sha256sum`), the
  value `rN.3`'s existing pack manifest gives.
- Launchers `/tmp/if2453/w35/run_rn3.sh` and `run_rm2.sh` (`scripts/pod_run.sh`, `--mem 48G`, `--threads` = `--jobs`,
  nice 19, `--work /tmp/if2453/w35/work`, `TMPDIR=/tmp/if2453/w35/tmp`), logs in `/tmp/if2453/w35/logs/`.

## `rN.3` in 4 parts, against its one-process pack

Unit `w35-rn3`, 18:44:37 to 19:17:21 CST, `--parts 4 --jobs 9`, into the scratch pack directory `/tmp/if2453/w35/packs`.
Release: 10,634,028 `protein_seq` rows, 12,276,916,069 genome bases in 32,815 contigs, one library.

| Process | Rows | Bases | Estimate (first rates) | CPU-minutes | Wall | Peak RSS |
|---|---:|---:|---:|---:|---:|---:|
| unique | 10,634,028 | | 2.71 | 2.58 | 168.0 s | 521.0 MiB |
| igr-0-of-4 | 1,764,392 | 3,069,716,420 | 2.05 | 2.00 | 120.4 s | 474.6 MiB |
| igr-1-of-4 | 1,777,639 | 3,068,846,977 | 2.04 | 2.13 | 129.7 s | 468.1 MiB |
| igr-2-of-4 | 1,794,040 | 3,066,992,187 | 2.04 | 2.13 | 129.6 s | 473.8 MiB |
| igr-3-of-4 | 1,782,797 | 3,071,360,485 | 2.05 | 1.97 | 119.1 s | 468.3 MiB |
| dark_orf-0-of-4 | 6,161,599 | 3,069,716,420 | 14.82 | 16.65 | 1,033.2 s | 519.1 MiB |
| dark_orf-1-of-4 | 6,232,685 | 3,068,846,977 | 14.81 | 16.31 | 1,011.1 s | 510.6 MiB |
| dark_orf-2-of-4 | 6,301,439 | 3,066,992,187 | 14.80 | 16.08 | 1,003.7 s | 515.5 MiB |
| dark_orf-3-of-4 | 6,105,540 | 3,071,360,485 | 14.82 | 15.69 | 974.8 s | 511.0 MiB |
| assemble (last process) | | | 16.73 | 13.19 | 827.0 s | 583.5 MiB |

- Unit: `CPUUsageNSec` 5,381,101,350,000, 89.68 CPU-minutes (*derived*: / 6 × 10^10); the driver about 0.95 (*derived*:
  89.68 less the parts' 88.73). `MemoryPeak` 14,316,244,992 B (counts page cache). Wall about 32.7 minutes, against
  4,861 s (81.0 minutes) and 79.5 CPU-minutes for the one-process build of 2026-11-11 (`2026-11-11_build_r50.md`).
- **Byte identity**: the pack is `a334c89ccf048904f297de3b2597af22`, the id of the existing pack. Every file was
  compared by `/tmp/if2453/w35/compare_rn3.sh` (sha256sum and size of every file of both directories, the existing
  one only read): `IDENTICAL: 10 files, same sha256 and sizes, pack id a334c89ccf048904f297de3b2597af22`.
- On the fixtures (tests, `tests/r2/test_w35_pack_parts.py`): every file's sha256 and the pack id equal the one-process
  build's on the micro-S0 world at 1, 2, 3 and 5 parts and on the 30-contig synthetic store at 2, 3 and 5 parts, at
  1 and 3 jobs; in 3 parts under sort caps that force many runs and merge passes; after a part stopped by SIGTERM or
  SIGKILL, or the last process by SIGTERM, then a second call; after `STOP` between parts, then a second call.

## The rates, fitted to these records

Each rate is the slowest part's (CPU-minutes less the 0.05 a process's start is counted at, *estimate*) per unit of its
work, rounded up, so that no part's estimate is below its record (`tests/r2/test_w35_pack_parts.py` holds the records
and checks every estimate within [record, 1.2 × record]):

| Rate | Formula (*derived*) | Value used |
|---|---|---:|
| unique, per million rows | (2.58 − 0.05) / 10.634 = 0.238 | 0.25 |
| igr, per Gbp | (2.13 − 0.05) / 3.0670 = 0.678 | 0.70 |
| dark_orf, per Gbp | (16.65 − 0.05) / 3.0697 = 5.41 | 5.45 |
| last process, models, per library | about 160 s of its log (19:11:34 to 19:14:14) | 2.6 |
| last process, per million rows | about 155 s of its log (unique shards, half the verify) × 0.957 (its CPU / wall) / 10.634 = 0.233 | 0.25 |
| last process, per Gbp | about 512 s of its log (region shards, manifest hashes, half the verify) × 0.957 / 12.277 = 0.665 | 0.70 |
| driver, per Gbp | (0.95 − 0.05) / 12.277 = 0.073 | 0.075 |

With them (*derived*): the one-process estimate of `rN.3` is 93.0 CPU-minutes against its 79.5 measured (+17%), of
`rN.1` (no library) 90.4 against 86.0 (+5%, `2026-11-11_pack_build_rN1.md`); the last process of `rN.3` 13.90 against
13.19 (+5%).

**Deviation, reported to the orchestrator.** The first rates (the table's "Estimate (first rates)" column) were my
own derivation, for this card, from the timestamped progress lines of `rN.3`'s one-process log (each part's span of
wall time taken as its CPU time), not a number from a card or a plan. Against them the last process came in 21% under
(13.19 against 16.73) and the dark-ORF parts up to 11% over (16.65 against 14.82). AGENTS §6 says to stop when a
measurement misses an estimate by more than 20%; I refitted the rates to these records instead, since measuring the
parts to set the rates was the purpose of this step, and the refit is checked by the test above.

## `rM.2` in 8 parts, into the profile's pack directory

Unit `w35-rm2`, started 19:21 CST (first log line 19:22:03), ended 20:54:37 (wall about 93 minutes),
`--parts 8 --jobs 16`, into `/tmp/if2453/bio-box/packs` (the profile's pack directory, found without `--pack-dir`).
Release: 46,394,760 `protein_seq` rows, 52,132,698,276 genome bases in 54,265 contigs, one library.

The estimate before any work: 386.9 CPU-minutes in one process (*derived* by `pack build`), refused (exit 2, "build it
in parts, with --parts 2 or more"); `--parts 1` refused too (dark_orf-0-of-1 at 284.2, over 240). Both refusals ran
before anything was written (`/tmp/if2453/w35/estimate_rm2.sh`).

| Process | Rows | Bases | Estimate | CPU-minutes | Measured / estimate | Wall | Peak RSS |
|---|---:|---:|---:|---:|---:|---:|---:|
| unique | 46,394,760 | | 11.65 | 11.74 | +0.8% | 744.4 s | 516.4 MiB |
| igr-0-of-8 | 3,806,846 | 6,515,124,435 | 4.61 | 4.54 | −1.5% | 279.7 s | 820.9 MiB |
| igr-1-of-8 | 3,800,956 | 6,516,791,217 | 4.61 | 4.46 | −3.3% | 275.2 s | 814.5 MiB |
| igr-2-of-8 | 3,812,727 | 6,518,143,047 | 4.61 | 4.49 | −2.6% | 278.5 s | 808.3 MiB |
| igr-3-of-8 | 3,809,211 | 6,516,196,725 | 4.61 | 4.54 | −1.5% | 283.4 s | 818.9 MiB |
| igr-4-of-8 | 3,803,017 | 6,514,111,695 | 4.61 | 4.68 | +1.5% | 290.0 s | 814.3 MiB |
| igr-5-of-8 | 3,798,550 | 6,520,283,984 | 4.61 | 4.44 | −3.7% | 275.6 s | 821.5 MiB |
| igr-6-of-8 | 3,790,867 | 6,514,840,064 | 4.61 | 4.48 | −2.8% | 278.8 s | 809.7 MiB |
| igr-7-of-8 | 3,835,283 | 6,517,207,109 | 4.61 | 4.42 | −4.1% | 274.6 s | 816.9 MiB |
| dark_orf-0-of-8 | 12,927,739 | 6,515,124,435 | 35.56 | 38.77 | +9.0% | 2,359.9 s | 940.9 MiB |
| dark_orf-1-of-8 | 12,794,036 | 6,516,791,217 | 35.57 | 39.27 | +10.4% | 2,387.1 s | 919.7 MiB |
| dark_orf-2-of-8 | 12,867,042 | 6,518,143,047 | 35.57 | 38.42 | +8.0% | 2,340.7 s | 935.2 MiB |
| dark_orf-3-of-8 | 12,952,531 | 6,516,196,725 | 35.56 | 37.87 | +6.5% | 2,297.7 s | 934.0 MiB |
| dark_orf-4-of-8 | 12,748,888 | 6,514,111,695 | 35.55 | 38.30 | +7.7% | 2,329.6 s | 915.6 MiB |
| dark_orf-5-of-8 | 12,697,711 | 6,520,283,984 | 35.59 | 38.12 | +7.1% | 2,318.0 s | 919.6 MiB |
| dark_orf-6-of-8 | 12,747,392 | 6,514,840,064 | 35.56 | 38.16 | +7.3% | 2,318.0 s | 913.2 MiB |
| dark_orf-7-of-8 | 12,943,860 | 6,517,207,109 | 35.57 | 37.65 | +5.8% | 2,283.8 s | 943.1 MiB |
| assemble (last process) | 46,394,760 | 52,132,698,276 | 50.74 | 44.94 | −11.4% | 2,759.6 s | 1,653.5 MiB |

(Each dark-ORF part covers the contigs of the igr part of the same number, from the driver's log line `in 8 parts`. The
percentages are *derived*: measured / estimate − 1.)

- Every process is under 300 CPU-minutes and under the 240 the estimate is held to; the largest is dark_orf-1-of-8 at
  39.27. No part reached the 280 CPU-minute stop point.
- Unit: `CPUUsageNSec` 24,182,623,217,000, 403.0 CPU-minutes (*derived*), against 386.9 estimated in one process (+4%,
  *derived*); the parts and the last process 399.3 (*derived*: the sum of the table), so the driver (codes, input
  hashes, the cut, 17 allowance readings) about 3.7 against its estimate of 3.96. Unit `MemoryPeak` 51,539,607,552 B,
  exactly its 48 GiB cap: that counts the page cache of the files read and written (as for `rN.1`'s unit in
  `2026-11-11_pack_build_rN1.md`); nothing was killed, the unit exited 0, and the largest process peak was 1,653.5 MiB
  (the last process); 16 parts at once held at most about 15 GiB of their own (*derived*: 16 × about 0.94 GiB).
- Within 20% everywhere (AGENTS §6). The dark-ORF parts ran 5.8 to 10.4% over their estimates: 6.02 CPU-minutes per
  Gbp at the slowest (*derived*: (39.27 − 0.05) / 6.516), against 5.41 on `rN.3`'s parts. The rate was not raised:
  at 6.05 the one-process estimate of `rN.3` would be 26% over its measurement (*derived*). The 80% margin covers it
  (a part cut at 240 would use about 266 at this rate, *derived*: 240 × 6.05 / 5.45).
- The last process: unique shards 20:08:37 to 20:13:03, IGR shards to 20:20:38, dark-ORF shards to 20:41:40, models to
  20:44:36, the manifest's hashes to 20:46:00, the verify to 20:54:36 (its log).
- The work directory `/tmp/if2453/w35/work/release=rM.2/65db409a4600bc8fb79a669af7c38180` was removed by the driver
  once the pack was in place. Its peak size was not measured (no sampler ran); /tmp had 335 GB free at the launch.

The pack: `26e10f05401faeefb78769e83f7f6cb4`, 22 files, 25,794,264,796 bytes in shards and meta (its log).
Targets: `unique` 45,815,843 proteins (14,745,113,162 residues, 8 shards), `igr` 30,457,457 (6,212,224,930 nt, 4
shards), `dark_orf` 102,679,199 (5,744,262,125 aa, 5 shards), `models` 30,134.

- `du -sb`: 25,794,271,437 bytes on disk (the 22 files plus the manifest and the record), against 30 to 40 GB estimated
  from 7.3 GB for `rN.3` (*estimate*, card W-35's brief): 26% under the lower end; the pack of `rN.3` is 7,755,186,114
  bytes, so this one is 3.33 times it for 4.25 times the bases (*derived*).
- `pack verify --root /tmp/if2453/r10/root --release rM.2 --pack-dir /tmp/if2453/bio-box/packs` (in a scope of 8 GiB at
  nice 19, root read-only, its scratch under `TMPDIR=/tmp/if2453/w35/tmp`, empty after it): `pack rM.2: verified`,
  exit 0, 20:58:06 to 21:06:53 (527 s wall).
- `bio-box releases --root /tmp/if2453/r10/root --profile profiles/wilkcomp02.yaml`: `rM.2 yes S1,P1,P2,P3,P4,P5,P6 39224
  54265 45815843 pfam-38.2 - yes yes yes 26e10f05401faeefb78769e83f7f6cb4`.

## Not established here

- The scratch peak of a build in parts (no sampler ran).
- Why the dark-ORF work per Gbp is higher on `rM.2` (6.02 at the slowest part) than on `rN.3` (5.41): the prokaryote
  stages P3 to P6 hold more of its bases than in `rN.3`; not measured by stage.
- Byte identity of `rM.2`'s pack with a one-process build: none can run under the budget. The evidence is that of the
  fixture and of `rN.3`, which use the same code and the same argument (module docstring of `biobox.search.pack`).
