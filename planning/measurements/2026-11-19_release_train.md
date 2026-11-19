# Card C-46: Release train rehearsal and real-root dry-run (wilkcomp02, 2026-11-19)

Card C-46 implements the rolling release-train driver (`DEC-52`, `DEC-53`, `DEC-54`, `B-49`, `W-35`):

- `scripts/stress/build_stage_pod/release_train.sh` (modes `--dry-run`, `run`, `rehearse`)
- `scripts/stress/build_stage_pod/train_plan.py` (naming rule, stage-selection and completeness check via `biobox.stage.stage_incomplete`, clustering estimate / `B-49` reuse detection, thread plan by the hour)
- `scripts/stress/build_stage_pod/train_est.py` (pre-launch estimate and `DEC-43` / `DEC-53` rate guard)
- `scripts/stress/build_stage_pod/train_check.py` (before/after listing diff, `root_added.tsv`, and `announcement.md` draft)
- `tests/r1/test_c46_release_train.py` (unit tests for naming, stage selection, completeness, `--dry-run`, rate guard, and post-release check)

Runtime image: `/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif` by digest
`sha256:d67300aa1321b571fa8da79e74dfccdd87a8dca2e56518f324902efa6b9f098a`. All times are CST.

## 1. End-to-end rehearsal on the fixture (2026-11-19, `wilkcomp02`)

Command (*measured*, `wilkcomp02`, 2026-11-19, `<= 4` threads total, `/tmp/if2453/r10/root` bound read-only):

```bash
WORK=/tmp/if2453/c46 TREE=/tmp/if2453/ci/C-46 \
  bash /tmp/if2453/ci/C-46/scripts/stress/build_stage_pod/release_train.sh rehearse
```

Fixture setup in `/tmp/if2453/c46/rh/root`:

1. Built stage `S0` (lambda, `called`) and stage `S1` (*M. genitalium*, `called`) with initial release `r26.11.1`
   over `S0,S1` (`--cluster-parts 2`, 2 threads): `release r26.11.1: verified (stages S0,S1, 41 files, 7 artifacts,
   manifest hash 1d98a189639af04f395aa53450723623)` (*measured*).
2. Built stage `S2` (T4, `called`) to completion, and stage `S3` at `stored` depth (`inputs` + `ingest` only) so the
   driver had both a newly complete stage and an ineligible `stored`-only stage.
3. Ran `release_train.sh --dry-run --date 2026-11-19 --parts 2 --no-live-cpu`:
   - `Next release: r26.11.2 (date 2026-11-19)`
   - `Latest defined release: r26.11.1`
   - `Stages (3): S0,S1,S2`, `New stages (1): S2`, `Action: cut`
   - `Excluded stages (1; left out, never waited for): S3 (stored): stored only (depth=stored; DEC-54 requires depth >= called)`
   - `Clustering: 0.02 CPU-min one-process (*estimate*) from 854 protein_seq rows; mode=in parts; cut={"one_process_cpu_min_max": 20000.0, "parts": 2}`
4. Ran `release_train.sh run --date 2026-11-19 --parts 2 --no-live-cpu` (`T=2 J=1 CPUS=2 Q=4 PACK_PARTS=2 PACK_JOBS=2`):
   - `build-stage` skipped all stage steps of `S2` and ran all release steps (`r26.11.2`) under `r50_watch.py` and
     `train_est.py rates`.
   - Post-build verification (*measured*):
     - `release r26.11.2: verified (stages S0,S1,S2, 56 files, 10 artifacts, manifest hash e1124ff5e1bef40e4f4a5dffb0b8e3d9)`
     - `release r26.11.1: verified (stages S0,S1, 41 files, 7 artifacts, manifest hash 1d98a189639af04f395aa53450723623)` (unchanged)
   - Search pack in 2 parts (`c46-rh-pack`): `pack r26.11.2: verified` (`pack_id=96fd14a3daaae665b8713449e35ac439`, *measured*).
   - Worked examples hook (`scripts/release_examples.sh r26.11.2`): executed end-to-end and wrote `summary.tsv`
     (`92 missing` as expected on the minimal fixture, which holds none of the examples' controls).
   - `train_check.py`: `278` files before, `478` after; `200` new (`1,578,982` bytes), `0` gone, `0` changed, `0` stray;
     wrote `root_added.tsv` (200 entries with size and sha256) and `announcement.md`; `checks: all pass` (*measured*).
5. Post-cut dry-runs (*measured*):
   - `--dry-run`: `Next release: r26.11.3`, `Latest defined release: r26.11.2`, `New stages (0)`, `Action: noop`.
   - `--dry-run --same-stages --parts 2`: `Action: cut`, `Clustering: REUSE from r26.11.2 (card B-49; no MMseqs2 run; protein_seq rows: 854)`.

Large rehearsal scratch (`/tmp/if2453/c46/rh`) was removed after verification.

## 2. Read-only dry-run against the real data root (`/tmp/if2453/r10/root`, 2026-11-19)

The first production release cut into `/tmp/if2453/r10/root` by `release_train.sh` will run once stage `C5` completes
(card C-47). Today, `r26.11.1` (card C-35) already covers all 9 complete stages (`S1,P1,P2,P3,P4,P5,P6,H1,PL1`), so
read-only dry-runs confirm both the noop path and the `--same-stages` `B-49` reuse path (*measured*, `wilkcomp02`,
2026-11-19, 05:17 CST):

### Default `--dry-run`

```text
Next release: r26.11.2 (date 2026-11-19)
Latest defined release: r26.11.1
Stages (9): S1,P1,P2,P3,P4,P5,P6,H1,PL1
New stages (0): (none)
Action: noop (0 new complete stages since latest release; pass --same-stages to cut anyway)
Excluded stages: none
Clustering: 15058.29 CPU-min one-process (*estimate*) from 49,902,840 protein_seq rows; mode=in parts; cut={"one_process_cpu_min_max": 20000.0, "part_cpu_min": 240.0}
Thread plan by the hour: day (08:00-20:00 weekdays) <= 40; night/weekend (20:00-08:00 weekdays & weekends (when other users < 0.5 cores)) <= 56; mem_scope=200G, cpu_min_per_process=300
Now (2026-11-19 05:17 CST): cap=56, allowance=52 (ours=4.13 cores, others=0.02 cores)
```

### `--dry-run --same-stages --parts 24`

```text
Next release: r26.11.2 (date 2026-11-19)
Latest defined release: r26.11.1
Stages (9): S1,P1,P2,P3,P4,P5,P6,H1,PL1
New stages (0): (none)
Action: cut
Excluded stages: none
Clustering: REUSE from r26.11.1 (card B-49; no MMseqs2 run; protein_seq rows: 49902840)
Thread plan by the hour: day (08:00-20:00 weekdays) <= 40; night/weekend (20:00-08:00 weekdays & weekends (when other users < 0.5 cores)) <= 56; mem_scope=200G, cpu_min_per_process=300
```

| Metric | Value | Label |
|---|---|---|
| Rehearsal `r26.11.2` manifest hash | `e1124ff5e1bef40e4f4a5dffb0b8e3d9` | *measured* (`release verify`, `wilkcomp02`, 2026-11-19) |
| Rehearsal `r26.11.1` manifest hash (before & after) | `1d98a189639af04f395aa53450723623` | *measured* (`release verify`, `wilkcomp02`, 2026-11-19) |
| Rehearsal `r26.11.2` pack id | `96fd14a3daaae665b8713449e35ac439` | *measured* (`pack build` + `pack verify`, `wilkcomp02`, 2026-11-19) |
| Rehearsal files added (`root_added.tsv`) | `200` files, `1,578,982` bytes (`0` gone/changed/stray) | *measured* (`train_check.py`, `wilkcomp02`, 2026-11-19) |
| Real-root `protein_seq` rows across 9 stages | `49,902,840` | *measured* (Parquet footers via `release_proteins`, `wilkcomp02`, 2026-11-19) |
| Real-root one-process clustering estimate | `15,058.29` CPU-min (`0.0` with `B-49` reuse over same stages) | *estimate* (`biobox.stage.cluster_cpu_min`) |
