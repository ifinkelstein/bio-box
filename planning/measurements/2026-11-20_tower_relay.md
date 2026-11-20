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
| `root_rM2_20261114/l3.tar` | 15,692,646,400 | 13:21:26 | not yet at 15:15 | | | *measured*, log |
| `added_r26_11_1.tar` | 20,908,011,520 | 13:35:12 | not yet at 15:15 | | | *measured*, log |
| `root_rM2_20261114/l1_and_records.tar` | 9,860,751,360 | queued | | | | sizes: STATUS (2026-11-14 durable copy) |
| `root_added_20261116/added_pl1.tar` | 3,545,118,720 | queued | | | | size *derived*: the seven tars' total less the other six |
| `root_added_20261115/added_h1_crispr_labels.tar` | 1,360,076,800 | queued | | | | STATUS (2026-11-15) |

**The mirror was still running when this note was written** (the log read at 15:15 CST, 2026-11-20; its last line
`13:35:12 added_r26_11_1.tar: try 1: copy started`): two of the seven tars, 55,418,501,120 bytes, were verified in
place on `tower`; the other five, 51,366,604,800 bytes (*derived*), were being copied or queued. The orchestrator
expected it to end at about 16:00. It ends with `SHA256SUMS` written on `tower` and a `DONE` line in the log; the
remaining rows of the table, and the rate of the run as a whole, are to be added from the log then.

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
