# Can a term bring both CRISPR estimates within 15% of the eleven built stages? (card B-60)

Answer: no term with a stated mechanism and an input known before the step starts does. Stage `C5` is the miss, and no
fit that includes `C5` moves it inside 15% for the near-Cas producer. Card B-60 ends here: no constant in `src/biobox/stage.py`
changes, no test is added or migrated, `RUNBOOK.md` keeps its rates.

All numbers are *measured* unless labelled *derived*. Read-only reads of the data root `/tmp/if2453/r10/root` on
wilkcomp02, 2026-11-21, with DuckDB inside `bio-box-test_faec1f4.sif` at `nice 19` (one process, 1 thread; nothing was
run beside the release train, no MinCED). Sources: `runs/stage=<S>/crispr-minced-0.4.2-w8-nr3-r23x47-s26x50.json` and
`...-nr2-r23x47-s26x50-cas1kb-caspfam1-pfam-38.2.json` (`cpu_min`, `started_at`, `wall_s`, `threads`), `l0seq/contig_index`
and `l0seq/genome_input` through each stage's `ingest.json` and `inputs.json` (the rule of `stage_bases` and
`crispr_flank_bases`), `l1/crispr/producer=*/stage=<S>` (arrays: distinct `region_uid`), and the step ledgers. The scratch
scripts that did this are in `/tmp/if2453/ci-tmp/B-60-meas/` (`table.py`, `calls.py`, `fit.py`, `fit2.py`, `*.json`); they
are not part of the repository.

## 1. The table

Scanned bases = stored bases + circular wrap-flank bases. Ratio = run record / code's estimate (`crispr_cpu_min`
at 168.1 CPU-s per Gbp, times 1.15 for near-Cas; *derived*). "Calls" is the sum of `n_in` over the near-Cas ledger's
contig rows (two-repeat calls before the filter). Arrays are default-producer arrays.

| Stage | Scanned Gbp | Contigs | Circular (flank Gbp) | Default arrays | Near calls | Default CPU-min | Ratio | Near CPU-min | Ratio | Default CPU-s/Gbp |
|---|---|---|---|---|---|---|---|---|---|---|
| S1 | 1.664 | 27,258 | 0 | 2,517 | 7,430 | 5.085 | 1.091 | 5.487 | 1.023 | 183.4 |
| P1 | 1.017 | 539 | 0 | 18,288 | 8,923 | 3.085 | 1.083 | 3.383 | 1.033 | 182.1 |
| P2 | 9.596 | 5,018 | 0 | 124,321 | 78,999 | 26.489 | 0.985 | 31.794 | 1.028 | 165.6 |
| P3 | 9.949 | 9,107 | 0 | 237,776 | 106,771 | 29.195 | 1.047 | 33.386 | 1.042 | 176.1 |
| P4 | 9.976 | 4,005 | 0 | 126,179 | 79,182 | 27.585 | 0.987 | 32.606 | 1.014 | 165.9 |
| P5 | 9.968 | 4,111 | 0 | 118,008 | 75,767 | 27.505 | 0.985 | 32.361 | 1.008 | 165.6 |
| P6 | 9.963 | 4,227 | 0 | 109,088 | 76,972 | 27.938 | 1.001 | 32.796 | 1.022 | 168.3 |
| H1 | 1.820 | 6,714 | 0 | 15,304 | 27,108 | 4.925 | 0.966 | 5.876 | 1.002 | 162.4 |
| PL1 | 12.130 | 61,647 | 56,548 (6.205) | 49,752 | 62,507 | 30.938 | 0.910 | 35.513 | 0.909 | 153.0 |
| V1 | 6.997 | 168,358 | 22,515 (1.819) | 15,439 | 36,548 | 20.560 | 1.049 | 21.893 | 0.971 | 176.3 |
| C5 | 9.109 | 2,653 | 0 | 62,928 | 30,939 | 21.361 | **0.837** | 22.257 | **0.758** | **140.7** |

Before (the code as it stands): 20 of 22 records are within 10% of their estimate (|record - estimate| / estimate), two
are not, both `C5`: default 16.3% under, near-Cas 24.2% under. After: unchanged, because no constant changes.

Beside the step. Cards C-33's records for `S1` and `P1` to `P6` (2026-11-15) ran as up to six units of 4 threads at once,
so the time-weighted threads of our other records overlapping a step are 5.7 to 15.1 (*derived* from `started_at` and
`wall_s` of all run records of the root). `H1`, `PL1`, `V1` and `C5` have no overlapping record of ours. The pod's other
users and the release train are in no run record, so the load at each step is not recorded (the load average read at the
start of this work was 54.7, 44.7, 44.0). Card B-58 measured 6.7% between a load of about 30 and about 55. The
data do not show a load effect larger than the spread: the stages that ran beside five others of ours (`P2` to `P6`:
165.6 to 176.1) are not above the ones that ran alone (`H1` 162.4, `PL1` 153.0, `C5` 140.7, `V1` 176.3).

## 2. What the data support

Fits are least squares on the relative error (rows divided by the record), on the eleven stages. "LOO" fits ten and predicts
the eleventh. Degrees of freedom = 11 minus the constants.

Default producer (CPU-seconds, Gbp scanned, thousands of contigs `k`, thousands of default arrays `a`):

| Model | Constants | dof | Worst in-sample | Worst LOO | `C5` in-sample / LOO (record / prediction) |
|---|---|---|---|---|---|
| 165.3 per Gbp | 1 | 10 | 14.9% (`C5` 0.851) | 16.7% | 0.851 / 0.833 |
| 162.3 per Gbp + 0.70 per k contigs | 2 | 9 | 13.4% (`C5` 0.866) | 15.7% (`C5` 0.843) | 0.866 / 0.843 |
| 159.4 per Gbp + 0.59 per k arrays | 2 | 9 | 14.4% | 19.2% | 0.861 / 0.838 |
| 139.5 per Gbp + 1.69 per k contigs + 1.84 per k arrays | 3 | 8 | 7.9% | 13.4% | 0.921 / 0.888 |
| 166.3 per Gbp - 1.27 per k circular contigs | 2 | 9 | 15.4% | 17.7% | 0.846 / 0.823 |

Near-Cas producer, absolute fit (CPU-seconds):

| Model | Constants | dof | Worst in-sample | Worst LOO | `C5` in-sample / LOO |
|---|---|---|---|---|---|
| 186.7 per Gbp | 1 | 10 | 21.5% | 24.3% | 0.785 / 0.757 |
| 185.5 per Gbp + 0.26 per k contigs | 2 | 9 | 21.0% | 24.9% | 0.790 / 0.751 |
| 161.0 per Gbp + 3.52 per k calls | 2 | 9 | 15.2% | 22.4% (`H1`, `C5`) | 0.848 / 0.776 |
| 153.7 per Gbp + 0.75 per k contigs + 4.04 per k calls | 3 | 8 | 12.6% | 23.2% | 0.874 / 0.770 |
| 156.4 per Gbp + 1.50 per k contigs + 2.38 per k default arrays | 3 | 8 | 15.4% | 22.5% | 0.846 / 0.775 |
| 144.5 per Gbp + 1.39 per k contigs + 1.50 per k arrays + 2.88 per k calls | 4 | 7 | 11.2% | 36.7% | 0.888 / 0.783 |

Findings.

1. **The default producer's miss is not a per-contig cost.** A per-contig constant (a JVM start per batch of contigs, a
   parse per contig) is the card's first suggestion and has a mechanism, but it is 0.70 ms per contig, and it moves `C5` by
   1.5 points (0.851 to 0.866). `C5` has the fewest contigs of the large stages (2,653, 3.4 Mbp each); it is low per base, which
   a per-contig term cannot do. The term helps `S1` and `V1` (27 and 168 thousand contigs) and hurts `P1`, `P3`: it is a
   different shape of error, not a repair, and it fails LOO (15.7%).
2. **The only default model that holds every stage inside 15% in LOO uses the number of arrays found, which is not
   known before the step runs** (it is the step's output). An estimate that needs its own result cannot guard a launch, so it
   is not usable. It also needs three constants on eleven points (8 degrees of freedom), and its leave-one-out worst is
   13.4% (`P3` 0.885, `V1` 0.866, `C5` 0.888), a margin of 1.6 points on the 15% bar.
3. **Near-Cas: no model comes inside 15% in LOO.** The cost over the default step (18 to 318 CPU-s) runs 5.9 to 33.2 CPU-s per
   Gbp and 1.7 to 4.4 ms per two-repeat call; neither is constant, and `C5` is at the low end of both (5.9 and 1.7). The
   candidate mechanisms (the nr2 search finds more calls, the filter reads the gene tape and domain rows, a DuckDB
   connection per run) are plausible, and none of them has a measurement here that separates it: the records hold no
   per-phase times. `C5`'s ratio to the best near-Cas model stays at 0.77 to 0.79 in LOO whatever the other terms.
4. **`C5` is a 15% under-spender per Gbp on both producers**, 140.7 (default) and 146.6 (near-Cas) CPU-s per Gbp against
   medians of 165.9 and 196.1. Its stage is the only one of complete prokaryote genomes in the table (9.11 Gbp, 3.4 Mbp
   mean contig); nothing in the records separates "complete genomes", "nothing of ours
   beside it" or the hour from it (`H1`, `PL1` also ran alone and read 162.4 and 153.0). One stage cannot support a constant
   (AGENTS §6 asks for a stated mechanism; none is stated).

## 3. Decision

- No term is written. A term built on inputs known before the step (bases, contigs, circular flank) moves `C5` by at most 3 points of 16
  (default) or 3 of 24 (near-Cas, 0.758 to 0.790 in-sample), and pulls stages that are inside 10% today toward the
  15% bar (default with contigs: `P1` 1.083 to 1.119, `P3` 1.047 to 1.081).
- Lowering the rate to follow `C5` is refused. Under DEC-53 a step more than 20% *under* its estimate is recorded and costs a
  refit card; a step *over* its estimate stops a build at night. Lowering `CRISPR_CPU_S_PER_GBP` by 10% would put `S1` at 1.21 and
  `P1` at 1.20 of their estimates, over the guard. The estimate is therefore left high for `C5`-like stages: it over-estimates
  them by 19% (default, 1 / 0.837) and 32% (near-Cas, 1 / 0.758), *derived*.
- Stages that cannot be brought in: `C5`, default (16.3% under) and near-Cas (24.2% under; the estimate is high by 32%), for
  the reasons in section 2. Not changed: `S1`, `P1` to `P6`, `H1`, `PL1`, `V1`; all within 10% of their estimates for both producers.
- What would change this (not done here): further stages of complete genomes (`SPX1`, later `C`-stages) that read below 150
  CPU-s per Gbp would make `C5` a class, and a per-phase timing of one near-Cas run (MinCED, the tape read, the filter),
  which needs a pod arm and a card, would give the near-Cas extra a mechanism. Either is a card for a later stage.
