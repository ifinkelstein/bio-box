# The search pack of `rN.1`, built by the streaming builder (card W-16, 2026-11-10 to 11)

Every number is *measured* on wilkcomp02 by the command or sampler named beside it, unless it is labelled
*derived* or *estimate*.

## Setup

- Code: branch `task/W-16` at `4f3237f`, mirrored to `/tmp/if2453/w16/ci/tree` with the rsync of
  `scripts/pod_pytest.sh` (`BIOBOX_POD_CI=/tmp/if2453/w16/ci sh scripts/pod_pytest.sh tree tests/r0/test_runbook.py`).
- Data: root `/tmp/if2453/r10/root`, release `rN.1` (stages `S1`, `P1`, `P2`; 12,276,916,069 bp, 32,815
  contigs; no domain library), read from its manifest. The root was bound read-only into the image
  (`--bind /tmp/if2453 --bind /tmp/if2453/r10:/tmp/if2453/r10:ro`); nothing was written under it, no STOP or
  HALT file was made, no `r10-*` unit was touched. Nothing under `/stor` was written.
- Image `bio-box-test_74b3e05.sif`, sha256 `c30be58557d756095346308c58cc3098af4bbe66404f67ab87be3f64f147c898`.
- Profile `/tmp/if2453/w16/profile.yaml`: `profiles/wilkcomp02.yaml` with `tmp: /tmp/if2453/w16/tmp`.
- Pod at launch (22:57 and 23:16 CST): load average 6 to 8; our other units were the 10 Gbp rebuild's.

```
bash /tmp/if2453/w16/ci/tree/scripts/pod_run.sh --profile /tmp/if2453/w16/profile.yaml \
  --image /stor/home/if2453/bio-box/images/bio-box-test_74b3e05.sif --digest sha256:c30be585...c898 \
  --mem 32G --threads 4 --bind /tmp/if2453 --bind /tmp/if2453/r10:/tmp/if2453/r10:ro \
  --env PYTHONPATH=/tmp/if2453/w16/ci/tree/src:/tmp/if2453/pilot_s1/pydeps --env OMP_NUM_THREADS=1 \
  --log /tmp/if2453/w16/logs/w16-pack.log w16-pack pack build --root /tmp/if2453/r10/root --release rN.1 \
  --pack-dir /tmp/if2453/w16/packs --profile /tmp/if2453/w16/profile.yaml
```

Samplers, each a `systemd-run --user` unit at nice 19 beside the build (scripts not committed, in
`/tmp/if2453/w16/`): `acct.sh` wrote the unit's `MemoryCurrent`, `MemoryPeak`, `CPUUsageNSec` and the cgroup's
`anon` and `file` bytes every 30 s to `logs/w16-pack.acct` (176 samples); `hwm.sh` wrote `VmHWM` and `VmRSS` of
each process of the unit every 15 s to `logs/w16-pack.hwm`.

## Result

The build ended `pod_run: w16-pack exited 0`, with pack `71917a4167c7aebe1b7ff44d18a99306` in
`/tmp/if2453/w16/packs/release=rN.1/`, and nothing else in that directory (no `.building-*`, no
`.scratch-*`). `bio-box pack verify` (unit `w16-verify`, same options) printed `pack rN.1: verified` and exited 0.
The pack is left in place for the next step, a real search.

| Build | Value |
|---|---|
| Wall time | 5,282 s, 1 h 28 min 2 s (`ActiveEnterTimestamp` 23:16:20 CST to the last log line, 00:44:22) |
| CPU | 86.0 CPU-minutes (`CPUUsageNSec` 5,161,180,595,000 at the last sample, within 30 s of the end); under the 300 CPU-minute rule |
| Peak memory of the build process | 767,980 kB = 750.0 MiB (`VmHWM` of the `python` process at its last sample, 13 s or less before it ended; `VmHWM` is the process's own high-water mark) |
| Largest anonymous memory of the unit | 717,324,288 B (cgroup `anon`, largest of 176 samples) |
| Unit `MemoryPeak` | 8,631,373,824 B: page cache of the files read and written (cgroup `file` 7,680,327,680 B at that sample), not the process's memory |
| Size | 5,510,724,073 B in 8 files (`du -sb` and the sum of the file sizes below) |
| `pack verify` | about 2 min wall (sampler, 00:52:47 to 00:54:47); `VmHWM` 331,408 kB = 323.6 MiB at its last sample |

Before W-16 the same build filled its 32 GiB cap in the dark-ORF part and was stopped after 9,528 s and
155.1 CPU-minutes without a pack (`planning/measurements/2026-11-10_trial_search_rN1.md`).

### Parts (from the timestamped progress lines in `/tmp/if2453/w16/logs/w16-pack.log`)

| Part | From | To | Rows | Note |
|---|---|---|---|---|
| genetic codes, input hashes | 23:16:23 | 23:17:13 | | `feature_protein` merged by contig; sha256 of every input |
| `unique` | 23:17:13 | 23:20:34 | 10,634,028 rows, 10,620,684 proteins | 163 sorted runs, one 64-way merge pass |
| `igr` | 23:20:34 | 23:32:39 | 9,491,097 IGR rows, 7,118,868 of at least 30 nt | collected by 23:27:25; 109 runs |
| `dark_orf` | 23:32:39 | 00:44:03 | 12,048,866 genes; 24,801,263 dark ORFs | ORFs of 32,815 contigs 23:33:42 to 00:33:25; 379 runs, sorted and written by 00:44:03 |
| models, manifest | 00:44:03 | 00:44:22 | 0 libraries | |

### Files

| File | Bytes | Records |
|---|---:|---:|
| `unique/part-0.faa.zst` | 1,242,095,090 | 6,143,374 (with part-1: 10,620,684, 3,340,668,639 aa) |
| `unique/part-1.faa.zst` | 905,151,551 | |
| `igr/part-0.fna.zst` | 606,018,363 | 7,118,868 (1,416,471,110 nt) |
| `igr/igr_meta.parquet` | 530,288,062 | 7,118,868 rows, 58 row groups |
| `dark_orf/part-0.faa.zst` | 1,313,881,169 | 23,601,162 |
| `dark_orf/part-1.faa.zst` | 66,819,238 | 1,200,101 (24,801,263 in all, 1,388,657,518 aa) |
| `dark_orf/dark_orf_meta.parquet` | 846,467,302 | 24,801,263 rows, 202 row groups |
| `pack_manifest.json` | 3,298 | |

The five files the pre-W-16 build finished in the trial have exactly these sizes there (unique 1,242,095,090
and 905,151,551; igr 606,018,363 and 530,288,062; dark ORF shards 1,313,881,169 and 66,819,238 B). On the
fixtures the two builders give byte-identical packs with the same `pack_id` (tests
`test_w16_*_equals_the_in_memory_oracle`; and a one-off probe that built the fixture, its 20,000-byte-shard
variant and the two-producer release with the old and the new module in one session: same `pack_id`, every
file equal).

## A first attempt, stopped (the reader defect it found)

The first launch (22:57:47, tree at `a6aa9dd`) reached `VmHWM` 1,809,888 kB = 1,767.5 MiB during the `unique`
part's 64-way merge (`logs/w16-pack-try1.hwm`). A probe in the image
(`/tmp/if2453/w16/readahead_probe.py`, not committed) showed why: 64 readers of 64-row-group files of 26 MB,
each asked for one batch, held 1,031,922,496 B of Arrow memory with `ParquetFile.iter_batches` and
27,760,512 B with `read_row_group`. pyarrow's `iter_batches` reads files ahead. The build was stopped with
`systemctl --user stop w16-pack` at 23:05:22 in the `igr` part: its log ends `stopped by a signal (exit 143);
removed .../.building-171683-03e9d04b and .../.scratch-171683-03e9d04b`, `pod_run: w16-pack received SIGTERM`,
`pod_run: w16-pack exited 143`, and the release directory was empty, as card W-16's signal rule says. Commit
`4f3237f` reads every file one row group at a time; the second launch above is that code.

## Against the estimates

The pack is larger than `planning/measurements/2026-11-10_search_pack_size.md` derived from the fixture: about
3.4 GB at 10 Gbp there, so about 4.17 GB at this release's 12.28 Gbp (*derived*, scaled by bases), against
5.51 GB measured, 32% more. The excess is in the nucleotide targets: 2,020 dark ORFs per Mbp here against 972
on the fixture (*derived*: 24,801,263 / 12,276.9 Mbp and 611 / 0.6286 Mbp), and IGR text twice the derived
0.31 GB zstd. The fixture (lambda and *M. genitalium*) is not representative of 5,557 prokaryote contigs. The
`unique` target is within its derived upper bound (2.15 GB against about 2.0 GB at 10 Gbp, 2.46 GB at 12.28).
