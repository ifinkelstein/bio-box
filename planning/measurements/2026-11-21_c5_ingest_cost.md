# What stage C5's ingest spent its CPU-minutes on (card B-58, amended 2026-11-21)

Question: `C5`'s ingest cost 45.803 CPU-minutes against the code's 38.057 (ratio 1.2035, DEC-56). Card B-58 first said
the 7.75 CPU-minutes over were zlib decompression at 0.84 CPU-min per GB. A reading test (`getrusage` of `gzip.open`
on `SPX1` files, 2026-11-21) found 0.09 per GB. This note measures what the ingest step itself spends, by arm.

All numbers are *measured* unless labelled *derived* or *estimate*. Host wilkcomp02, 2026-11-21 02:46 to 05:39 CST,
script `scripts/stress/b58_ingest_arms.py` (modes `inputs`, `scale-inputs`, `run`, `report`), logs in
`/tmp/if2453/b58/logs/` (`arms.log`, `arms2.log`, `b58-<arm>-<rep>.log`, one `pod_run.sh` unit each) and
`/tmp/if2453/b58/results.jsonl`. The figure of an arm is `cpu_min` of the run record of its ingest step
(`runs/stage=BA/ingest.json`), the quantity the C-47 rate guard judged. Each arm was a `build-stage --depth stored` in a
fresh scratch root, a detached `pod_run.sh` unit of 2 threads at `nice 19` in `bio-box_faec1f4.sif` (digest
`sha256:d67300aa...`, the image C5 was built in), one arm at a time, `STOP` checked before each. The data root was never
written; nothing under `/stor`. The roots were removed after each arm.

## 1. The samples

| Sample | Objects | FASTA | Rule |
|---|---|---|---|
| `c5` (arms gz, plain, few, view, held) | 204 | 944,466,060 bytes (276,224,309 on disk as `.gz`), 510 contigs | `C5/objects.json[::10]` (index 0, 10, 20, ...), one genome per object |
| `c5x3` (arm scale) | 678 | 3,120,188,970 bases, about 3.16e9 bytes (*derived*: bytes per base of `c5`); 923,764,326 on disk | `objects.json[::3]` |
| `spx` | 4 `.fa.gz` of `SPX1` | 1,237,247,281 bases (379,769,300 on disk) | the middle 4 files in size order (ranks n/2-2 to n/2+1): `SAMN05421999`, `SAMN08323206`, `SAMN05416669`, `SAMN19433754` |

The code's estimate of the `c5` sample is 0.07 + 0.944466 x 4.11 = 3.952 CPU-min (*derived*, `ingest_cpu_min`, nothing held).

## 2. Arms (ingest CPU-min; the load average at the start of the arm in brackets)

| Arm | What differs | rep 1 | rep 2 | rep 3 | mean | mean / estimate |
|---|---|---|---|---|---|---|
| gz | `.fna.gz` objects (as C5) | 4.4646 (38) | 4.2027 (27) | 4.2291 (38) | 4.2988 | 1.088 |
| plain | the same genomes as plain copies | 3.9943 (37) | 4.2102 (21) | 4.2284 (52) | 4.1443 | 1.049 |
| few | the same bases as one plain object | 3.8344 (21) | 4.1851 (57) | 4.2043 (56) | 4.0746 | 1.031 |
| view | `gz` in a view of the 9 earlier stages | 4.0466 (23) | 3.9525 (35) | 4.3381 (56) | 4.1124 | 1.041 |
| held, stage A | `gz` as stage A of a root | 4.0956 (20) | 4.1417 (22) | 4.4781 (55) | 4.2385 | 1.073 |
| held, stage B | the same objects again: every contig held | 0.3333 | 0.3725 | 0.3803 | 0.3620 | |
| gz, batch 2 | as gz, interleaved with scale (reps 4 to 6) | 4.4539 (55) | 4.5305 (56) | 4.4613 (56) | 4.4819 | 1.134 |
| scale | `objects[::3]`, gz, 3.3 times the sample | 14.6850 (34) | 14.8125 (39) | 14.5195 (41) | 14.6723 | 1.124 (of 13.054, *derived*) |

(`results.jsonl` has the active units and the load at the start and end of every arm. The pod's other units were the
release train `c46-train-rel-1` throughout and, from rep 2 on, the `v2-rh-*` and `zx1-rh-*` rehearsals of other cards:
load 35 to 57.)

What the arms say, by the card's candidates:

1. **gzip against plain.** Pairs from the same rep: gz minus plain is +0.470, -0.008 and +0.001 CPU-min, mean +0.155 per
   0.944 GB, 0.164 per GB of FASTA (*derived*); against all six gz runs (mean 4.390) it is 0.26 per GB. The spread of one
   arm across reps (3.83 to 4.46) is as large as the difference. Two passes read each object (`store`: pass 1 every
   record, `_write_shard` in pass 2 the objects with new sequences again), so the reading tests predict 2 x 0.083 = 0.17
   per GB for `SPX1`'s files and 2 x (0.266 - 0.045) = 0.44 per GB for C5's own (*derived*, reading tests of
   2026-11-21). Decompression costs about 0.16 to 0.26 CPU-min per GB, which on C5 is 1.6 to 2.5 of the 7.75
   CPU-minutes (*derived*), not all of it.
2. **Cost per object.** `few` (one object) is 1.7% under `plain` (204 objects): no per-object cost (4.07 against 4.14,
   within the spread).
3. **Held contigs.** Stage B, all 510 contigs held: 0.362 CPU-min per 0.944 GB, 0.383 per GB, 8.5% of what a stored
   byte costs (4.24 / 0.944 = 4.49 per GB). The estimate charges them 0. C5 holds 267,273,180 bases: 0.10 CPU-min
   (*derived*). Not the cause. The `view` arm (the earlier stages' `contig_index` read by the anti-join) is 1.041 of the
   estimate, like `plain`: not the cause.
4. **Load.** `plain`, `few`, `view` and `held A` average 4.04 in rep 1 (load 20 to 37) and 4.31 in rep 3 (load 52 to
   56): 6.7% more (*derived*). C5's ingest ran at load 30 (section 3).
5. **Scale.** `scale` costs 4.64 CPU-min per GB of FASTA (14.672 / 3.159, *derived*), the `gz` arms 4.55 (batch 1) and
   4.75 (batch 2): the cost per byte does not grow with the size of the stage.

Does a cause explain C5's record to within 15%? The scale arm (a third of C5's own `.fna.gz` files) at its cost per
byte puts C5 at 4.644 x (9.510 - 0.267 x 0.915) + 0.07 = 43.1 CPU-min (*derived*, held bases at 8.5% of a stored
byte), 5.9% under the record of 45.803. That is within 15% but names no single cause: C5-type genomes in plain copies
already cost 4.39 per GB (`plain`: 4.144 / 0.944), at the high end of the 3.94 to 4.41 the nine fitted stages ranged
over (card B-52), and decompression adds 0.16 to 0.26 per GB on top. A decompression term alone brings C5 to
38.057 + 0.164 x 9.51 = 39.6 (ratio 1.156, outside 15%, *derived*) with the paired difference, or 38.057 + 0.26 x 9.51 =
40.5 (1.131, inside) with the unpaired one: the arms do not decide between them, because the spread from run to run
(about 0.2 CPU-min, 5%) is as large as the effect. **No term is written** (card: "a term fitted to one record with no
mechanism is not written"). The arms rule out a per-object cost, held contigs, the earlier stages' index and the size
of the stage. They leave gz decompression of 0.16 to 0.26 per GB (a mechanism, measured twice), the genome type of C5
(4.39 per GB in plain) and load.

## 3. What else ran beside C5's ingest, and the threads (read-only: `driver.log`, run records)

The ingest run records (`runs/stage=<S>/ingest.json`, `threads` and `wall_s`, read 2026-11-21): threads asked of ingest
were `S1` 16, `P1` 16, `P2` 16, `P3` 16, `P4` 16, `P5` 16, `P6` 32, `H1` 36, `PL1` 36, `V1` 6, `C5` 6. Each ingest ran
as one process (CPU time over wall: `S1` 0.997, `P2` 0.998, `P3` 0.964, `PL1` 1.031, `H1` 1.047, `V1` 1.017, `C5` 1.045,
*derived*), so the thread count is not what moved C5: `V1` also had 6 and ran at 1.021 of its estimate.

Beside C5's ingest (`/tmp/if2453/c47/driver.log`, 2026-11-19 11:22 to 12:08): the run-step line at 11:24:10 says "ours
elsewhere 17.18" CPUs (unit `c42-step4`, with `c47-driver-run`), other users 0.00 CPUs, load average 30.01 (the driver's
snapshot at 11:22:45: 26.03). For `V1`: ours elsewhere 3.84, load 21.55. For the stages the rate was fitted on (`S1` to
`PL1`, 2026-11-10 to 16) the drivers' logs are not under `/tmp/if2453/{c47,c48c,c24,c26,r10,r50}/logs` (no
`run-step ... ingest: allowance` line there, *measured* 2026-11-21) and their run records do not say what ran beside
them. The comparison with the fitted stages cannot be made from the records.

## 4. What the arms say for `SPX1` (304 gzip objects, 256.9 GB of FASTA, 85.23% of bases under 5 kb)

Arms `spx` (every contig kept) and `spx5000` (`--min-contig-len 5000`; 1,237,247,281 bases read, 273,549,222 kept,
963,698,059 read and dropped): `spx` 6.0604, 6.7426, 6.6227 (mean 6.4752); `spx5000` 1.6454, 1.7576, 1.7557 (mean 1.7196).
The code estimates 0.07 + 1.2372 x 4.11 = 5.155 for both (*derived*): `spx` is 1.256 of it, `spx5000` 0.334.

Two-point model, cost = a x (Gbp read) + b x (Gbp kept) (*derived*): b = (6.4752 - 1.7196) / 0.96370 = 4.94 CPU-min per
Gbp stored; a = (1.7196 - 4.94 x 0.27355) / 1.2372 = 0.299 CPU-min per Gbp read and not stored (it carries gzip
decompression, the hashing of every record and the 0.07 start-up).

`SPX1` at `--min-contig-len 5000` (*derived*): 14.77% of 256.9 Gbp stored = 37.94 Gbp x 4.94 = 187 CPU-min, plus 256.9
x 0.299 = 77 read and dropped: **about 264 CPU-minutes, against the code's 1,056** (the code counts every byte as
stored). With C5-type genomes' stored rate (4.39 per Gbp) in place of `SPX1`'s many small contigs (4.94): 167 + 77 = 244.
So the code is about four times too high for `SPX1`, and a guard would judge the step at 0.25 (DEC-53: under 20% is
recorded, not stopped); the ingest fits in 1 or 2 processes of 240 CPU-minutes, not the 5 shards the code's rate gives.
How the first-pass shards should be cut is a card for the PI. This rests on 4 files and 3 reps with a spread of about
5% per arm, and the fraction kept in the sample (22.1%) is not the stage's (14.77%).

## 5. Open points

- A decompression term (0.16 to 0.26 per GB) is real but explains about a quarter to a third of C5's excess. Writing it
  and refitting the rate with `C5` and `V1` among the fitted stages (eleven) changes every estimate: not byte for byte,
  so not this card's.
- Whether the estimate should price bases read and dropped separately (a = 0.299 per Gbp) from bases stored
  (b = 4.94) is a design choice: `--min-contig-len` is a case the estimate does not know.
- The arms ran 2 threads; C5 and V1 had 6 and the fitted stages 16 to 36. Ingest used about one CPU in every record,
  so no effect is expected, but it was not measured (the cap allowed 2).
