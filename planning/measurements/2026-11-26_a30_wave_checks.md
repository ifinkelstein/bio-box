# What a wave's rehearsal and release verification cost, and what sharing them saves (card A-30, DEC-81 (3))

Card A-30 lets the stages of one wave share the rehearsal and the verification of every release before a build
(DEC-81 (3)). Its step (1) asks for the time of each part of `rehearse`, `before` and `after` first, so that the saving
of each share is known before the design. This note gives those parts, then the saving for a wave of 18 and of 25
stages.

## 1. Source and method

*Measured*, wilkcomp02: the driver logs of the three stages built with `stage_driver.sh` so far, read on 2026-11-26 and
copied with `scp` (read-only; nothing on the pod was run or touched): `/tmp/if2453/stages/{g1,spx1,v2}/driver.log`,
the rehearsals' inner logs `/tmp/if2453/stages/{g1,spx1,v2}/rh/drv/driver.log`, `/tmp/if2453/stages/spx1/queue.out`
and `/tmp/if2453/stages/r26-11-4-wave.out`. The lines were selected with `rg -o` on bounded patterns, for example

```
rg -o "^\S+ \S+ CST \[\w+ (before|after|rehearse|check|backup)\] stage driver: mode \w+ (\(pid|done)" driver.log
rg -o "^\S+ \S+ CST \[G1 [a-z]+\] (verify \S+: exit \d+|snapshot \S+: \d+ files|rehearse: [a-z ]{0,40})" driver.log
```

Each part is the difference of two log timestamps (millisecond resolution): a mode from its `stage driver: mode M
(pid` line to its `mode M done` line; a release's verification from the line before it to its `verify R: exit 0`
line; the snapshot from the last `verify` line to the `snapshot` line (in `before`) or from the mode's start to the
`snapshot` line (in `after`, which takes it first). The dates: G1 2026-11-22 (two attempts; the first stopped in `run`
on the rate guard), SPX1 2026-11-23, V2 2026-11-25 (its `run` has not passed yet, so V2 has no `after`).

## 2. The parts (*measured*, minutes unless marked s)

### `rehearse`

| Part | G1 (2nd) | SPX1 | V2 |
|---|---|---|---|
| set-up (caps, the hour's cap, earlier logs moved) | 3.4 s | 3.4 s | 3.4 s |
| the fixture, stage `E` and release `rRH.1` in scratch (unit `<s>-rh-e`) | 28.19 | 24.93 | 28.79 |
| inner `estimate` | 4.7 s | 5.1 s | 5.0 s |
| inner `before` (verify of `rRH.1`, snapshot) | 4.3 s | 4.0 s | 4.5 s |
| inner `run` (unit `<s>-rh-1`) | 4.00 | 1.93 | 4.02 |
| inner `after`, `check`, `backup` | 13.4 s | 14.1 s | 14.1 s |
| relaunch that must publish nothing (`<s>-rh-again`) | 18.5 s | 13.5 s | 12.8 s |
| **the mode** | **32.94** | **27.53** | **33.48** |

G1's first attempt took 32.94 as well (06:28:56.324 to 07:01:52.459). The card gives 36 for V2; its log gives 33.48
(7% lower; the queue's lines of V2 in `r26-11-4-wave.out`, 15:34:42 to 16:08:11, agree with the log). G1 and V2 are annotated (Pfam 38.2,
CRISPR); SPX1 is called with `MIN_CONTIG_LEN=5000`. Stage `E` is 86 to 91% of a rehearsal (*derived*: 28.19 / 32.94,
24.93 / 27.53, 28.79 / 33.48 are 0.86, 0.91 and 0.86).

### `before`

| Part | G1 (1st) | G1 (2nd) | SPX1 | V2 |
|---|---|---|---|---|
| `release verify` of every release (8 for G1 and SPX1, 9 for V2) | 13.73 | 12.09 | 12.28 | 19.71 |
| snapshot (listing, sha256 of every other stage's file) | 2.98 | 2.97 | 3.44 | 4.89 |
| **the mode** | **16.72** | **15.07** | **15.73** | **24.60** |

Per release, V2 (s): `r26.11.1` 149.7, `r26.11.2` 280.4, `r26.11.3` 500.0, `rM.2` 133.6, `rN.1` 31.3, `rN.1a` 8.3,
`rN.1b` 10.1, `rN.2` 34.4, `rN.3` 35.0. G1 (2nd, s): `r26.11.1` 212.8, `r26.11.2` 259.3, `rM.2` 132.8, `rN.1` 32.2,
`rN.1a` 8.1, `rN.1b` 10.2, `rN.2` 34.9, `rN.3` 35.1. The verification grows with each release published: `r26.11.3`
alone added 8.33 to V2's.

### `after`, `check`, `backup`

| Part | G1 | SPX1 |
|---|---|---|
| `after`: snapshot | 2.96 | 3.56 |
| `after`: `release verify` of every release | 12.10 | 12.19 |
| `after`: `bio-box releases` | 4.4 s | 4.1 s |
| **`after`** | **15.13** | **15.82** |
| `check` (the passing run) | 1.95 | 1.88 |
| `backup` | 3.37 | 11.61 |

## 3. What each share saves

- A shared rehearsal saves the whole `rehearse` mode but the share's own work: the code's sha256 over about 6 MB
  (`du -sh src tests/fixtures/micro_s0 scripts/stress/build_stage_pod`: 2.7 M, 428 K, 2.8 M, *measured* on the laptop
  2026-11-26) and a byte comparison per earlier marker, a few seconds (*estimate*; not timed on the pod).
- A shared `before` saves the verification part only; its snapshot is still taken (the baseline of `check`), and so
  is the sha256 of every other stage's file. *Derived*: 12.09 to 19.71 per stage; 19.71 with the 9 releases the
  stages of `r26.11.4` see.
- `after` is not shared (card A-30 (3)): 15.13 to 15.82 per stage stay. `check` and `backup` stay.

## 4. The saving of a wave (*derived*)

R_a = 33.21, the mean of the two annotated rehearsals (G1 2nd 32.94, V2 33.48); R_c = 27.53, the one called rehearsal
(SPX1); R = 31.32, the mean of the three; V = 19.71, V2's verification of 9 releases.

**18 stages: the rest of `r26.11.4`** (`PL2` to `GART2` in the order of card A-26). By their confs they are of four
kinds: annotated with Pfam 38.2 gathering and CRISPR (11: `PL2`, `PL3`, `MGE1`, `SAT1`, `ECE1`, `O1`, `O2`, `LPL1`,
`SYM1`, `LRV1`, `MGE2`), called with CRISPR (5: `CG1`, `MV1`, `V3`, `GART1`, `GART2`), called with `gv-meta`
(`LGRT1`), called with `MIN_CONTIG_LEN=5000` and `gv-meta-c25` (`LR1`). One rehearsal of each kind runs; 14 are shared.
Every `before` but the first is shared: 17.

- rehearsals: 10 x R_a + 4 x R_c = 332.1 + 110.1 = 442.2 min
- verifications: 17 x V = 335.1 min
- **total: 777.3 min = 12.96 h**, if the wave runs without a resume. Each resume is a new wave and costs one more
  rehearsal of each kind still to come and one more full verification.

**25 stages** (the wave of `r26.11.5`; its kinds are not fixed yet): saving = (25 - k) x R + 24 x V for k kinds. With
four kinds, as in `r26.11.4`: 21 x 31.32 + 24 x 19.71 = 657.7 + 473.0 = **1,130.7 min = 18.85 h**. V will be larger by
then: `r26.11.4` will be one more release to verify in each `before` (*estimate*: about 8 to 10 min more, as
`r26.11.3`'s 8.33), which adds to the saving of each shared `before` by the same amount.

DEC-81 (3) put "up to about 20 h saved" to the PI for the 18-stage wave (*estimate*). The 12.96 h here are less
because card A-30 keeps every `after` in full (18 x about 12.1 min of verification, 3.6 h, *derived* from G1 and SPX1)
and runs the first rehearsal of each kind (4 x about 31 min, 2.1 h); with those two the figure would be about 18.7 h.

## 5. What a shared check still proves and what it no longer does

A shared rehearsal rests on an earlier full rehearsal, in the same wave, of the same kind, under the same tree, code
sha256, image and digest, profile and limits. It no longer exercises the driver on this stage's own conf; nothing in
a rehearsal reads more of the conf than the kind (the rehearsal's conf is written from `DEPTH`, the library, `CRISPR`,
`MIN_CONTIG_LEN`, `GENECALL_MODE` and the `NOTE` lines, which `read_conf` has already checked against the tree for
this stage).

A shared `before` rests on an earlier full verification of the same wave and root: every file of the root has the
same path, size and mtime as when it was verified, every other stage's file the same sha256 (hashed again now), each
manifest the same sha256, and no file is newer than the record. What it no longer re-reads is the content of the files
outside the stages' partitions (a release's own `release=` files, archived records): a change that kept a file's size
and mtime would be found only by the next `after`, which verifies every release in full.
