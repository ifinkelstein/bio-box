# The NAS `tower` and the relay through the laptop: capacities and rates (card C-49)

Card C-49, 2026-11-20. What `tower` can hold, how fast files reach it from the pod, and what the first mirror of the
root's durable tars took. No transfer was started for this note. Its numbers come from the card's own first
measurements (the orchestrator, from the laptop, 2026-11-20, `planning/tasks/TASKS.md` card C-49, last paragraph),
from `STATUS.md` (open item 19, handoff section 6e), and from the log of the first mirror, `relay_backups.log`, written
by the orchestrator's script `relay_backups.sh` in its session's scratch directory (not in the repository; its general
form is now `scripts/orchestrator/relay_to_tower.sh`). Every number is labelled as AGENTS §6 asks. The design that
uses these numbers is `planning/design/C-49_tiers.md`.

## 1. Capacities

| What | Value | Label and source |
|---|---|---|
| `tower` array | 32 TB, 18 TB free | *measured*, laptop, 2026-11-20, `df` (card C-49) |
| `tower` data disks, free | 1.8, 7.3 and 8.8 TB | *measured*, the same (card C-49 (d)) |
| `tower` cache disk | 954 GB, 242 GB free | *measured*, the same |
| `tower` cores, memory | 8 cores, 31 GB | *measured*, the same, `nproc`, `free` |
| Share `research-bulk` | cache off: writes go to the array; already holds other projects' bulk data | *measured*, read-only (card C-49 (d)) |
| Share `research` | prefers the cache disk (242 GB free), which a large copy would fill | *measured*, the same |
| Layout made on `tower` | `/mnt/user/research-bulk/bio-box/{root_backups,sources,tranches}/` | STATUS handoff 6e |
| Laptop disk free | 26 GiB | *measured*, laptop, 2026-11-20, `df` (card C-49 (b)); so nothing may land on it |
| Pod local disk (`/tmp`, the data root's place) | 1.1 TB in all, 613 GB free; cleared at the reboot of 2026-11-17 | *measured*, wilkcomp02, 2026-11-20 10:22, `df` (STATUS open item 19) |
| Pod `/stor/work` (the group's) | 68 TB, 19 TB free; ours 1.2 TB (JGI archives 1.14 TB, the root's tars 107 GB) | *measured*, wilkcomp02, 2026-11-20, `df`, `du` (STATUS open item 19) |
| Pod `/stor/scratch`, ours | 164 GB (163 GB of it `provisional/`) | *measured*, the same |

Who reaches whom (card C-49): the laptop reaches the pod (through `wilkcomp03` off campus) and `tower` (on the PI's
home network, as `root@tower` with the laptop's own key). The pod cannot reach `tower`. `tower` holds no key for the
pod (the PI: "no tower pull key").

## 2. Rates

| Measurement | Value | Label and source |
|---|---|---|
| One file, one stream, pod to laptop to `tower` (`scp -3`), `bio-box_faec1f4.sif` | 821,764,096 bytes in 84.75 s: **9.70 MB/s** | *measured*, laptop, 2026-11-20 (card C-49 (c)); rate *derived*: bytes / seconds |
| Its sha256 read back on `tower` | `d67300aa…f098a`, equal to the pod's; file and directory then removed | *measured*, the same |
| Two streams at once, copy phase: `root_rM2_20261114/l0seq.tar` | 28,241,715,200 bytes in 6,543 s: 4.32 MB/s | *measured*, `relay_backups.log` line 4; rate *derived* |
| Two streams at once, copy phase: `root_rM2_20261114/l2.tar` | 27,176,785,920 bytes in 6,843 s: 3.97 MB/s | *measured*, `relay_backups.log` line 6; rate *derived* |
| Two streams together | about 8.3 MB/s (4.32 + 3.97); the orchestrator's reading during the run was about 8.4 | *derived* from the two rows above, which overlapped for 6,543 s; so a second stream does not add to one stream's 9.70 |
| Read-back of `l0seq.tar` on `tower` (`sha256sum`, then `mv -n`) | 3,084 s for 28,241,715,200 bytes: 9.16 MB/s, while `l2.tar` was still being copied | *derived*: the log's OK time 13:21:26 less the copy's end (10:40:59 + 6,543 s = 12:30:02); bytes / seconds |
| Read-back of `l2.tar` on `tower` | 3,609 s for 27,176,785,920 bytes: 7.53 MB/s, while `l0seq.tar`'s read-back and then `l3.tar`'s copy ran | *derived*: 13:35:11 less 12:35:02 |
| Each of the two tars from start to verified in place | 9,627 s and 10,452 s | *derived* from the log's start (10:40:59) and OK lines |
| The first two tars together, start to verified | 55,418,501,120 bytes in 10,452 s: 5.30 MB/s | *derived* |

Not measured: which hop limits the 9.70 MB/s (the pod's side, the internet path, or the laptop's home link); the rate
in the other direction, `tower` to laptop to pod, which would use the home connection's upload and may be slower;
the hash taken as the stream arrives on `tower` (the repository script's default, section 4).

## 3. The first mirror: the root's seven durable tars

Started by the orchestrator at 10:40:59 on 2026-11-20 from the laptop (`caffeinate -i bash relay_backups.sh`, two
streams, `scp -3`, then `sha256sum` read back on `tower` and `mv -n`), from `/stor/work/Wilke/if2453/all-the-data/root_backups`
to `tower:/mnt/user/research-bulk/bio-box/root_backups`; 106,785,105,920 bytes in seven tars (*measured*, card
C-49 and the audit note `2026-11-20_audit_of_1117_1120.md` §2).

| Tar | Bytes | Copy started | Verified in place | Copy (s) | Copy rate | Label |
|---|---|---|---|---|---|---|
| `root_rM2_20261114/l0seq.tar` | 28,241,715,200 | 10:40:59 | 13:21:26 | 6,543 | 4.32 MB/s | *measured*, log; rate *derived* |
| `root_rM2_20261114/l2.tar` | 27,176,785,920 | 10:40:59 | 13:35:11 | 6,843 | 3.97 MB/s | *measured*, log; rate *derived* |
| `root_rM2_20261114/l3.tar` | 15,692,646,400 | 13:21:26 | 15:15:19 | 3,768 | 4.16 MB/s | *measured*, log; rate *derived* |
| `added_r26_11_1.tar` | 20,908,011,520 | 13:35:12 | 15:59:56 | 4,738 | 4.41 MB/s | *measured*, log; rate *derived* |
| `root_rM2_20261114/l1_and_records.tar` | 9,860,751,360 | 15:15:19 | 16:08:42 | 1,649 | 5.98 MB/s | *measured*, log; rate *derived* |
| `root_added_20261116/added_pl1.tar` | 3,545,118,720 | 15:59:56 | 16:19:33 | 641 | 5.53 MB/s | *measured*, log; rate *derived* |
| `root_added_20261115/added_h1_crispr_labels.tar` | 1,360,076,800 | 16:19:33 | 16:25:45 | 209 | 6.51 MB/s | *measured*, log; rate *derived* |

**The mirror ended at 16:56:54 with `DONE: all 7 in place; SHA256SUMS written on tower`** (*measured*, the log; added
by the orchestrator when it ended). Every tar's sha256, read back on `tower` before it was put in place, equals the
recorded value, and the seven were read a second time there for `SHA256SUMS` (16:25 to 16:56) with the same seven
values. `du -sb` on `tower` gives 106,785,106,597 bytes for the directory (the seven tars and the 677-byte
`SHA256SUMS`). From the first byte to the last tar verified in place: 5 h 44 min 46 s for 106,785,105,920 bytes,
**5.16 MB/s end to end** (*derived*); to `DONE`, 6 h 16 min.

**The repository script's first real use** (the orchestrator, 16:57:31 to 16:57:54, `relay_to_tower.sh` with
`BIOBOX_POD_JUMP` set, one stream, default hashing on arrival): three deferred-data lists that were on the pod's
`/tmp` only, from `/tmp/if2453/s2/inputs` to `tower:/mnt/user/research-bulk/bio-box/deferred_lists`:
`V2/V2_left_for_later_v4.0_putative.tsv.gz` (106,448,343 bytes in 12 s, 8.87 MB/s),
`V2/V2_left_for_later_hc.tsv.gz` (58,947,247 bytes in 6 s, 9.82 MB/s) and `PL2/skipped_records.tsv` (17,445,698
bytes in 3 s), each with its sha256 equal to the one taken on the pod, then `SHA256SUMS` written (3 lines)
(*measured*, the script's log; the directory listed on `tower` afterwards). A dry run first printed the same three
files and 182,841,288 bytes. So the script works against the real hosts, and on files of 60 to 100 MB it reaches
verified-in-place at about 9 MB/s, where the first mirror's read-back halved it. Its rate on a file of tens of GB is
still to be measured at the next backup tar.

**The script on a multi-GB file, the same evening** (*measured*, the script's log, 19:38:33 to 19:42:06): `added_v1.tar`
(stage `V1`'s durable tar, 4,231,157,760 bytes) from `/stor/work/.../root_backups` to `tower`'s `root_backups` in
211 s, **20.05 MB/s to verified in place**, sha256 `0a285dc3…a4ec` equal to the record; `SHA256SUMS` there extended
to 8 lines. Also `uhgv_metadata.tsv.zst` (38,978,745 bytes, 2 s) to `deferred_lists/V1/`. That is twice the 9.70 MB/s
of the morning's single-stream test and four times the first mirror's end-to-end 5.16 MB/s. One file at one hour of
the evening: whether the morning's rate was the link at that hour, or the two streams and read-backs sharing it, is
not separated. For planning keep 9.70 MB/s; a tar of 20 to 30 GB at the next release will say more.

**Release `r26.11.2`'s 25.2 GB tar (2026-11-22, 06:09:59 to 06:31:30 CST)** (*measured*, `relay_r26_11_2.log`):
`added_r26_11_2.tar` (25,200,107,520 bytes, 278 members) streamed from `/stor/work/Wilke/if2453/all-the-data/root_backups`
to `tower:/mnt/user/research-bulk/bio-box/root_backups` in **1,287 s**, **19.58 MB/s to verified in place** (*derived*:
25,200,107,520 / 1,287 / 1e6), sha256 `d44999e2384042ff2c665a86d579483d95d163f3c045aea0c7127b133e3bfa21` equal to the
`/stor/work` record, `SHA256SUMS` extended to 10 lines (`added_c5.tar` was line 9). That confirms `added_v1.tar`'s
20.05 MB/s within 2.3% across a 6× larger archive: `ssh pod cat | ssh tower 'tee .partial | sha256sum'` sustains ~19.6 to
20.1 MB/s end to end when idle, while 9.70 MB/s remains a conservative daytime planning floor.

**Stage `SPX1`'s 75.67 GB staged inputs (2026-11-22, 06:32:54 to 07:34:27 CST, DEC-66 `S-2`)** (*measured*, `relay_spx1.log`):
all 309 staged files of `/tmp/if2453/s2/inputs/SPX1` (`genomes.json`, `objects.json`, `samples.tsv`, `input_report.json`,
`assemblies.sha256`, and 304 `.fa.gz` assemblies; 75,672,465,511 bytes in all) streamed to
`tower:/mnt/user/research-bulk/bio-box/sources/SPX1` in **3,693 s** (61.55 min), **20.49 MB/s end to end** (*derived*:
75,672,465,511 / 3,693 / 1e6), every file verified equal to its pod sha256 on arrival and `SHA256SUMS` written (309 lines,
`0 failed`).

**Stage `GV1`'s 4.69 GB staged inputs (2026-11-22, 07:39:12 to 07:42:24 CST, DEC-66 `S-2`)** (*measured*, `relay_gv1.log`):
all 7 staged files of `/tmp/if2453/s2/inputs/GV1` (`gvmags_v2_all_18kgenomes.fna` at 4,639,876,745 B in 182 s = 25.49 MB/s,
plus `records.tsv`, `genomes.json`, `shared_records.tsv`, `repeated_records.tsv`, `objects.json`, `skipped_members.tsv`;
4,693,449,816 bytes in all) streamed to `tower:/mnt/user/research-bulk/bio-box/sources/GV1` in **192 s**, **24.45 MB/s end to
end** (*derived*: 4,693,449,816 / 192 / 1e6), `SHA256SUMS` written (7 lines, `0 failed`).

**Stage `G1`'s 2.96 GB staged inputs (2026-11-22, 07:46:10 to 07:48:04 CST, DEC-66 `S-2`)** (*measured*, `relay_g1.log`):
`g1_globdb_files.tar` (deterministic tar of `/tmp/if2453/s3/globdb/genomes/G1/files/`'s 7,139 `.fa.gz` archives:
2,909,265,920 B, sha256 `5b4a2a91fcec430d1999b93f817aad84972f5b2c6ea69d99202671baf8bd2406`, 105 s = 27.71 MB/s) plus
`genomes.json`, `objects.json`, `shared_records.tsv`, `lock.json`, `selection.tsv`, and `stream_summary.json` (2,961,546,429 B
in all) streamed to `tower:/mnt/user/research-bulk/bio-box/sources/G1` in **114 s**, **25.98 MB/s end to end** (*derived*:
2,961,546,429 / 114 / 1e6), `SHA256SUMS` written (7 lines, `0 failed`). Also mirrored `C1`..`C6`'s 12 staged `genomes.json` and
`objects.json` lists (6,173,661 B in 12 s, `SHA256SUMS` 12 lines; `relay_c_inputs.log`).

**Stage `V2`'s 9.63 GB DEC-65 re-staged inputs (2026-11-22, 08:22:48 to 08:30:30 CST, DEC-66 `S-2`)** (*measured*, `relay_v2.log`):
all 8 DEC-65 staged files of `/tmp/if2453/s2/inputs/V2` (`imgvr_v4.1_hc_hq_votu_reps.fna` at 9,458,228,347 B in 442 s =
21.40 MB/s, `inovirus_roux2019.fna` 50,500,520 B, `records.tsv`, `genomes.json`, `trims.tsv`, `shared_records.tsv`,
`repeated_records.tsv`, `objects.json`; 9,634,042,117 B in all) streamed to `tower:/mnt/user/research-bulk/bio-box/sources/V2`
in **462 s**, **20.85 MB/s end to end** (*derived*: 9,634,042,117 / 462 / 1e6), `SHA256SUMS` written (8 lines, `0 failed`).

**Stage `PL2`'s 16.68 GB DEC-65 staged inputs (2026-11-22, 09:29:53 to 10:03:58 CST, DEC-66 `S-2`)** (*measured*, `relay_pl2.log`):
all 15 staged files of `/tmp/if2453/s2/inputs/PL2` (`gut_plasmids_yu2024.fasta`, the 7 public `plasmidscope_*.fasta` sets,
`genomes.json`, `objects.json`, `records.tsv`, `trims.tsv`, `skipped_records.tsv`, `stage_report.json`, and
`gut_plasmids_yu2024_circular_non_representatives.tsv`; 16,675,968,181 B in all) streamed to
`tower:/mnt/user/research-bulk/bio-box/sources/PL2` in **2,045 s**, **8.15 MB/s end to end** (*derived*:
16,675,968,181 / 2,045 / 1e6), `SHA256SUMS` written (15 lines, `0 failed`).

**Stage `PL3`'s 18.43 GB DEC-65 staged inputs (2026-11-22, 10:05:08 to 10:17:52 CST, DEC-66 `S-2`)** (*measured*, `relay_pl3.log`):
all 8 staged files of `/tmp/if2453/s2/inputs/PL3` (`plasmidscope_IMG-PR.fasta` at 17,157,795,650 B in 705 s = 24.34 MB/s,
`plasmidscope_mMGE.fasta` at 1,008,445,970 B in 39 s = 25.86 MB/s, `records.tsv`, `genomes.json`, `trims.tsv`,
`stage_report.json`, `objects.json`, `skipped_records.tsv`; 18,429,944,704 B in all) streamed to
`tower:/mnt/user/research-bulk/bio-box/sources/PL3` in **764 s**, **24.12 MB/s end to end** (*derived*:
18,429,944,704 / 764 / 1e6), `SHA256SUMS` written (8 lines, `0 failed`). Every queued stage's `/tmp` inputs (`G1`, `SPX1`,
`V2`, `GV1`, `PL2`, `PL3`, and `C1..C6` lists: 128,073,589,808 B across 366 files) are now verified on `tower`.

## 4. What the first mirror showed, and the change in the repository script

With `scp -3` the bytes cross the link once and are then read again from the array to be hashed: at about 9 MB/s
while a copy runs (section 2), so each file costs its copy time plus most of that again before it is in place, and
the two streams' reads and writes compete on the array's disks. Two streams moved 55.4 GB in 10,452 s, 5.30 MB/s from
start to verified (*derived*, section 2), against 9.70 MB/s for one stream's copy alone.

`scripts/orchestrator/relay_to_tower.sh` therefore streams `ssh pod cat` into `ssh tower 'tee .partial | sha256sum'`:
`tower` hashes the bytes as they arrive, the script compares that hash with the record and the `.partial`'s size
with the pod's, and only then moves it into place. No byte is read back by default; `--reverify` keeps the read-back.
At one stream the expected rate to verified is then the copy rate, about 9.70 MB/s (*estimate*: it assumes that
hashing on `tower` keeps pace with the stream, which one core does many times over, and that one stream alone copies
as fast as the measured single file; neither is measured for this script). Measure it at the first real use and
record it here; if it misses 9.70 MB/s by more than 20%, stop and report (AGENTS §6).

From 9.70 MB/s (*derived*: bytes / 9.70e6 / 3,600): 100 GB in 2.9 hours, 1 TB in 28.6 hours, 3 TB in 85.9 hours
(3.6 days), with the laptop awake and on the home network throughout. With `--reverify`, add the read-back at about
10 MB/s per file (*measured* for one file read while two ran, card C-49): 100 GB in about 2.8 hours more (*derived*).
