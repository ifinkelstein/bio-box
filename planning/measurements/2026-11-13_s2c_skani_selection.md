# S2C phase 2: skani run and the ANI-diversity selection (card C-27 step 4; 2026-11-13)

Step (4) of card C-27 (DEC-31). Steps (1) to (3) are in `2026-11-13_s2c_skani_fetch.md`. Every number is
labelled *measured* (host, date, command), *derived* (formula and inputs) or *estimate*. Times are CST. Host:
wilkcomp02. skani 0.4.0 ran only from the image `/tmp/if2453/images/bio-box_64f1bf0.sif`
(`ghcr.io/ifinkelstein/bio-box@sha256:2faf4b8de020430d3c5c9dace5569445fdb0eb0057088972793874a799e25bea`), the
tree's scripts (stdlib only) run as `python -I` from the image under `--cleanenv`, synced to
`/tmp/if2453/ci/C-27-4`. Written
only under `/tmp/if2453/c27`, `/tmp/if2453/ci/C-27-4`, `/tmp/if2453/ci-tmp/C-27-4-*` (pytest) and, for the
archive, the new directory `provisional/s2c_phase2/` on `/stor`.

## Result in one paragraph

The driver as written from skani's README would have failed on the first species: skani 0.4.0's `sketch` writes
one consolidated database unless `--separate-sketches` is given (fixed, with a test on real skani output and a
stand-in binary that now behaves as 0.4.0 does). On one species (*Xanthomonas campestris*, 104 genomes) skani
cost 0.065 CPU-s per genome sketched and 11.4 CPU-ms per pair; the whole run, all 28,700 candidates that have a
sequence and the 499,130 farthest-first pairs plus the quality top 20's 46,170, took **100.8 CPU-minutes** in four
single-thread processes of 19 to 29 CPU-minutes (*measured*), against the 4 to 17 CPU-hour *estimate*. skani
reported every pair asked for, and no pair's alignment fraction came near 15% (smallest larger fraction 56.47%), so
**the alignment-fraction rule removed nothing**. The ANI rule keeps phase 1 as a subset and changes **2,347 of the
4,617 genomes** picked after rank 1 in the 243 species (241 species change); the quality-only picks held **near-clones (ANI >
99.9) in 230 of 243 species, 6,930 of their 46,170 pairs**; the ANI picks hold them in 79 species, 2,895 pairs,
nearly all in species that are clonal as a whole. The median over species of the minimum pairwise ANI among the 20
falls from 97.59 to 97.12. The remainder stays 13,485 genomes and becomes 56.21 Gbp (55.91 before). The 4,617
genomes phase 2 keeps from the fetched candidates are archived on `/stor` as one 5.13 GB tar.

## 1. The driver against the real binary

`skani --help`, `skani sketch --help`, `skani dist --help` and `skani triangle --help` in the image, and a run on
four real genomes (`/tmp/if2453/c27/skani_probe/`, 09:09; three *Legionella pneumophila*: GCF_001752745.1 and
GCF_001752785.1, 3,447,085 bp each, GCA_000277065.1, 2,682,626 bp; and GCA_000007385.1 of another species), all
*measured* 2026-11-13:

| Question | skani 0.4.0 | The driver as written from the README |
|---|---|---|
| `sketch -l LIST -o DIR` | writes **one consolidated database** (`sketches.db2`, `index.db2`, `markers.bin*`, `README.txt`); per-genome `<FASTA file name>.sketch` files only with `--separate-sketches` | expected per-genome `.sketch` files: **every species would have failed** ("sketches missing"). Fixed: `--separate-sketches` |
| output folder | refused if it exists ("Output directory exists ... Exiting", exit 1) | the driver already removes it first |
| `.fna.gz` | read directly (sketch and dist gave the same values from `.fna.gz` as from the sketches) | as assumed |
| columns of `dist` and `triangle -E` | `Ref_file Query_file ANI Align_fraction_ref Align_fraction_query Ref_name Query_name`; the file columns hold the FASTA path a sketch was made from, not the `.sketch` path | as assumed (README) |
| which fraction is whose | `Align_fraction_ref` is the reference's own aligned fraction: GCA_000277065.1 (2.68 Mbp) against GCF_001752745.1 (3.45 Mbp) prints 99.91 for itself and 77.75 for the larger genome; 100 × 2,682,626 / 3,447,085 = 77.82 (*derived*) | as assumed |
| `--min-af` | "Only output ANI values where **one** genome has aligned fraction **>** than this value [default: 15]", that is the larger of the two fractions, strictly; `--both-min-af` (off) would ask both | the driver used skani's default, so a pair under 15% would have been absent, indistinguishable from a pair under the screen. Changed: dist and triangle run with `--min-af 0`, and `c25_select.py --min-af 15` applies the rule (a pair is unreported when its larger fraction is under 15; a value of exactly 15.00 counts as reported there and was hidden by skani: see section 5) |
| an unreported pair | **absent** from the output (no row), also with `-s 0 --min-af 0` for the other species' genome; the driver writes it as a row with an empty `ani` (`--expect-*`, `fill_expected`) | as assumed |
| the screen | `-s 80`: pairs of approximately < 80% identity are not computed | not used by the driver |
| `-t` | threads; for one `dist` query it splits the references: on the species of section 3, 2 threads halved the wall time (29.6 s to 15.3 s) at the same CPU (28.5 and 28.8 CPU-s) and gave a byte-identical pair table (sha256 `a04ddbda…`) | as assumed |

Tests, first failing then passing (`sh scripts/pod_pytest.sh C-27-4 tests/r1/test_c25_select.py` with
`BIOBOX_POD_IMAGE=/tmp/if2453/images/bio-box-test_64f1bf0.sif`):

- `test_skani_normalise_real_output` (new): skani's own output from the probe, embedded in the test (two `dist`
  rows, three `triangle -E` rows; 1.8 kB), normalised; the fractions on the right genome (checked against the
  genome sizes, an independent check), the absent pair filled as empty when asked for. It passed at once: the
  normaliser's reading of the columns was right.
- `test_skani_driver_run_with_stand_in_binary` (changed): the stand-in binary now behaves as 0.4.0 does
  (consolidated database without `--separate-sketches`, refuses an existing folder, prints only pairs with one
  fraction > `--min-af`, prints FASTA paths); the test asserts `--separate-sketches` and `--min-af 0` on every call
  and that low-fraction pairs reach the table with their values. Failed first with the driver unchanged:
  `RuntimeError: sp000_W: sketches missing for ['GCA_126287286.1', ...]`; later also asserts the CPU record.
- `test_ani_unavailable_genomes_are_left_out`, `test_skani_driver_plan_leaves_out_unavailable` (new): see section
  2. Failed first: `c25_select.py: error: unrecognized arguments: --unavailable`, and `c27_skani.py plan` exit 2.
- `test_ani_species_median_near_clones_and_table_counts` (new): the columns of section 5 against a constructed
  world. Failed first: `KeyError: 'ani_median_pairwise_ani'`.
- `test_fetch_archive_keeps_the_picked_genomes` (new): section 6. Failed first: `c27_fetch.py archive` exit 2.

Before: `18 passed`; with the new tests and the code unchanged: `4 failed, 18 passed`; after: `23 passed` (with the archive test of section 6).
Whole suite on the pod in the test image (`pod_pytest.sh C-27-4 --shards 4`, 10:02 to 10:54): 1,737 passed, 20
skipped, 0 failed (shards: 427 passed 2 skipped; 582 passed 4 skipped; 435 passed 11 skipped; 293 passed 3
skipped).

## 2. `paths.tsv`

The 243 rank-1 genomes are phase 1's, on local disk as `.fna` under `/tmp/if2453/r50/genomes/P3..P6/` (193) and
`/tmp/if2453/ladder/genomes/P1, P2/` (50), read in place (read only); no tar was opened. The other candidates are
this card's `.fna.gz` under `/tmp/if2453/c27/genomes/B01..B12/`.

`c27_skani.py paths /tmp/if2453/c27/ani/paths.tsv` over the six phase-1 locks and the twelve batch locks
(*measured*, 09:12): 40,423 accessions (all of phase 1's stages P1 to P6 and the batches), sha256
`53165671c52825a541690752aa9eb7f78747c4722bb27b944251712df1f7fc52`. Against the 28,726 candidates
(`sel/ani_candidates.tsv`, sha256 `d49ee8b9…`): **28,700 have a path and the file exists** (28,457 batch files and
243 phase-1 files); **26 have none, exactly the 26 not served** (`ani/not_served.txt`, from the batches'
`fetch_result.json`, sha256 `a5c730d0c7572f4caa1936eb2c82eb4b26e03f48e3743b487d86c45e6be90d64`); no phase-1
genome is missing.

The 26 have no sequence, so no ANI. Left as they were, the pick would count every pair with them as unreported
(the largest distance) and so pick them first. They are listed with the new `--unavailable` (in `plan` and in
`c25_select.py --ani`): in a species over the cap they are left out of the farthest-first candidates; a rank-1
genome cannot be (refused); the quality rule they are compared with is then the first 20 of the genomes that have a
sequence, and `selection.json` counts the difference to the full quality-only list separately.

## 3. One species first: *Xanthomonas campestris*

104 candidates, 526,376,077 bp, 5.06 Mbp per genome (above the candidates' mean of 4.17 Mbp, *derived*:
119,700,310,300 / 28,726), no unavailable genome. `c27_skani.py plan` on that species alone, then
`run-shard` under `/usr/bin/time -v` in unit `c27-measure-t1` (`CPUQuota=100%`, `nice 19`, 09:19:04 to 09:19:34)
and again with 2 threads (`c27-measure-t2`). *Measured*:

| | 1 thread | 2 threads |
|---|---|---|
| wall (process) | 29.64 s | 15.34 s |
| CPU (user + system; unit `CPUUsageNSec`) | 28.50 s; 28.52 s | 28.79 s; 28.80 s |
| `sketch` (104 genomes) | 6.741 s | 3.242 s |
| 19 `dist` rounds (1,786 pairs) | 20.330 s | 10.581 s |
| `triangle -E` of the quality top 20 (190 pairs; 153 not already in the rounds) | 1.825 s | 0.884 s |
| peak RSS (largest process) | 67,200 kB | 70,700 kB |
| pairs in the table; unreported | 1,939; 0 | the same bytes |

- Sketch: **0.065 CPU-s per genome** (6.741 / 104), 12.8 ns per bp (*derived*), against the *estimate* of 0.5 to
  2 CPU-s: 8 to 31 times less.
- Pairs: **11.4 CPU-ms per pair** in the rounds (20.330 / 1,786; each round loads the remaining genomes' sketches
  from disk, so this includes the loading), 9.6 in the triangle (1.825 / 190), against the *estimate* of 1 to 10
  ms: 14% above its upper end.
- Sketch size on disk: 55,016,880 bytes for 104 genomes, 0.1045 bytes per bp, 529 kB per genome (*measured*,
  `du` of the `.sketch` files); for all 119.6 Gbp about 12.5 GB (*derived*), on `/tmp`.

**Projection** (*derived*): sketching 119,571,741,801 bp × 12.8 ns = 1,530 CPU-s = 25.5 CPU-min; pairs (499,130
farthest-first + 46,170 quality top) × 11.4 ms = 6,220 CPU-s = 103.7 CPU-min; **about 129 CPU-min, 2.2
CPU-hours**, against the 4 to 17 CPU-hour *estimate* and far under the 60 CPU-hour stop. Per process, the largest
species (*E. coli*, 4,003 genomes, 76,057 pairs) is about 19 CPU-min (*derived*).

The measured costs miss the card's *estimate* by more than 20% (AGENTS §6, §9), downward for sketching and
slightly upward per pair; the brief's stop rule for this step was a projection over 60 CPU-hours, so the run went
on and the miss is reported here. **All pairs** within the species would be 15,817,516 pairs × 11.4 ms ≈ 50
CPU-hours (*derived*; *E. coli* alone 8,010,003 pairs, 25 CPU-hours, to be split into some ten processes), 24
times the farthest-first work: **only the farthest-first pairs (and the quality top 20's triangle) were
computed.**

Shards were replanned from the measured costs plus 25% (`--sketch-cpu-s 0.08 --pair-ms 14`, now the defaults)
under a budget of 45 CPU-min per process (a sixth of the 300 CPU-minute rule, to keep shards short and room for a
fivefold misestimate): 4 shards, estimated 44.9, 44.9, 44.9 and 30.8 CPU-min (shard 0 holds *S. aureus*, *E.
coli* and *K. pneumoniae*), 165.5 in all (*estimate*).

## 4. The run

`c27_skani.py plan sel/ani_candidates.tsv ani/paths.tsv /tmp/if2453/c27/skani --threads 1 --budget-cpu-min 45
--image <the image> --unavailable ani/not_served.txt` (09:23; `plan.json` sha256 `7c5d1a48…`, `members.tsv`
`6a7a7992…`): 243 species, 28,700 genomes, 119,571,741,801 bp, 499,130 farthest-first pairs (C-27's 499,624 less
19 rounds × 26 genomes not served), 46,170 quality-top pairs; all pairs would be 15,817,516.

**Which pairs.** Only the pairs the greedy farthest-first pick consults: in each species, 19 sequential rounds,
`skani dist` from the newest pick to every genome not yet picked (so every pair among the 20 picks is computed),
plus one `skani triangle -E` of the 20 best-quality genomes for the comparison with the quality rule. All pairs
would be about 24 times dearer (section 3) and say nothing more about the pick.

**Launch.** By day, so at most 4 of our threads; the account was then at about 33 to 34 CPUs (the r10 Pfam scan's
eight 4-thread shards). A small gate unit (`c27-skani-gate`, a shell loop, `CPUQuota=50%`) started one
single-thread shard a minute while the account's CPU, summed from `ps -u if2453 -o pcpu=` over processes at least
5 s old, plus one stayed at or under 40 (the measuring `ps` itself shows several hundred percent at age 0 s, so
the unfiltered sum read 38 to 41 at 09:18 to 09:24 while `top` read 32 to 35). Shards `c27-skani-0..3` were
started at 09:25:08, 09:26:29, 09:27:50 and 09:29:10 (the gate's sums: 34, 36, 35, 36), each `systemd-run --user`,
`nice 19`, `MemoryMax=32G`, `MemorySwapMax=0`, `CPUQuota=100%`, `singularity exec` of the image. From 09:29 to
09:51 we ran 4 threads; none at night.

*Measured* (each shard's `DONE`: `cpu_s` is user + system of the process and its skani children):

| Shard | Species | Genomes | Wall s | CPU s | CPU-min | Peak RSS of a child |
|---|---|---|---|---|---|---|
| 0 | 3 (*S. aureus*, *E. coli*, *K. pneumoniae*) | 7,790 | 1,920.7 | 1,725.4 | 28.8 | 2,247,168 kB |
| 1 | 13 | 7,793 | 1,839.7 | 1,662.0 | 27.7 | 991,872 kB |
| 2 | 67 | 7,784 | 1,634.8 | 1,500.3 | 25.0 | 155,688 kB |
| 3 | 160 | 5,333 | 1,283.6 | 1,162.6 | 19.4 | 50,316 kB |
| all | 243 | 28,700 | | 6,050.3 | **100.8** | |

- Every process stayed under a tenth of the 300 CPU-minute rule. *E. coli* alone took 904.4 CPU-s (15.1 CPU-min)
  against 19 *derived* in section 3. The whole run was 22% under the projection of section 3 (129 CPU-min) and 39%
  under the plan's padded *estimate* (165.5).
- The largest skani process held 2.2 GB (a `dist` round over *E. coli*'s 4,003 sketches).
- Work directory `/tmp/if2453/c27/skani`: 13,614,384,499 bytes (sketches and round files; kept). `/tmp` had 470 GB
  free at 09:58 (`df`).

`c27_skani.py collect` (unit `c27-select-ani`, 09:58): **`/tmp/if2453/c27/ani/ani_pairs.tsv`, 513,636 pairs,
sha256 `2c95c71038b0de851435bbcc7317d48cf99a1324e8abadd860a5467d792a68c1`** (columns genome_a, genome_b, ani, af_a,
af_b; the 499,130 consulted pairs and 14,506 quality-top pairs the rounds had not covered). *Measured* on it:

- **Unreported pairs: 0.** skani printed every pair asked for: no pair fell under its screen.
- **Alignment fraction:** the larger fraction of a pair is 56.47% at the least (1st percentile 72.77, median
  89.06); the smaller is 40.11% at the least (median 84.95); 326 pairs have a smaller fraction under 50%. No pair is
  near the 15% rule, so **it removed nothing**, in any species; the exactly-15.00 case of section 1 does not occur.
- ANI from 93.26 to 100.00 (median 98.37); 331 pairs under 95; 3,733 printed as 100.00.

## 5. The selection

`python -I c25_select.py /tmp/if2453/c25/gtdb /tmp/if2453/c25/phase1/all.tsv /tmp/if2453/c27/out_ani
--phase1-not-served <the 7> --ani ani/ani_pairs.tsv --require-pairs --unavailable ani/not_served.txt` (unit
`c27-select-ani`, in the image; 36.1 s, peak RSS 260,328 kB, *measured*). Inputs: the GTDB R226 metadata as in
`sources.yaml`, phase 1's `all.tsv` (sha256 `e68d2576…`). Every pair the pick consulted was in the table
(499,130 consulted, 0 absent, 0 unreported). The comparison is with C-25's quality-only outputs
(`/tmp/if2453/c27/sel/`, byte-identical to C-25's).

Outputs under `/tmp/if2453/c27/out_ani/` (*measured*, `sha256sum`):

| File | sha256 |
|---|---|
| `s2c_selected.tsv` | `928b6e71b045233a1ad084cff45186246a4c84895223a4904e044f73886dcda1` |
| `remainder.tsv` | `48e6072846f6b40e396eea2639db61a2c02528be18e704feafaacc7462b93098` |
| `by_phylum.tsv` | `b3a99bef530603218ffca0843b2498b4feb93de760ba1d4319d515b348681902` |
| `ani_species.tsv` (committed beside this note) | `1f6e061ba8e10300651d8b81563cc07f4153e2ebe3f000ef707499d993ec7fd5` |
| `selection.json` (committed beside this note) | `44951fc8f37b41702f80df3aafa5f213c2318a8a104787bf21ffd8df9b7a43cf` |

**Phase 1 stays a subset** (asserted by the script; `phase1.not_in_selection` is empty, as it was for the quality
rule). Every species over the cap still has 20 picks.

**Genomes that differ from the quality-only list:** 2,347 (`selected_differ_from_quality_all_ranks`: in the ANI
selection and not in C-25's `s2c_selected.tsv`; the same number leaves). Measured per species against the first
20 genomes that have a sequence: 2,343, in 241 of the 243 species; per species median 9, 90th percentile 17, most
19 (nearest rank over `ani_species.tsv`, *derived*). The four more against C-25's list are ANI picks among the
genomes that move into the first 20 when the 8 not-served genomes of C-25's top 20 are left out (last part of this
section).

The ten largest species (*measured*, `ani_species.tsv`; ANI values as skani printed them; "near-clone" = pairs of
the 190 with ANI > 99.9):

| Species | Complete | Not served | Differ | Min ANI, quality | Min ANI, ANI rule | Median, quality | Median, ANI rule | Near-clone pairs, quality | Near-clone pairs, ANI rule |
|---|---|---|---|---|---|---|---|---|---|
| *Escherichia coli* | 4,007 | 4 | 19 | 97.72 | 96.20 | 98.34 | 97.22 | 61 | 0 |
| *Klebsiella pneumoniae* | 2,312 | 0 | 19 | 98.81 | 95.18 | 99.02 | 98.65 | 12 | 0 |
| *Salmonella enterica* | 1,833 | 0 | 19 | 98.18 | 95.04 | 98.58 | 98.09 | 13 | 0 |
| *Staphylococcus aureus* | 1,475 | 0 | 19 | 97.96 | 96.52 | 98.97 | 97.80 | 13 | 0 |
| *Pseudomonas aeruginosa* | 920 | 0 | 19 | 98.97 | 97.14 | 99.28 | 98.61 | 3 | 0 |
| *Bordetella pertussis* | 792 | 0 | 18 | 99.76 | 98.22 | 99.97 | 98.92 | 139 | 2 |
| *Acinetobacter baumannii* | 750 | 0 | 18 | 97.41 | 97.12 | 97.88 | 97.68 | 3 | 0 |
| *Campylobacter_D jejuni* | 522 | 0 | 19 | 97.77 | 94.76 | 98.51 | 97.19 | 15 | 0 |
| *Mycobacterium tuberculosis* | 505 | 0 | 18 | 99.09 | 98.80 | 99.91 | 99.74 | 103 | 0 |
| *Bacillus velezensis* | 496 | 0 | 19 | 98.31 | 97.14 | 99.28 | 98.06 | 8 | 0 |

**Does the rule matter?** Distribution over the 243 species (nearest-rank quantiles over `ani_species.tsv`,
*derived*; the median is the lower median of the 190 pairs; every one of the 46,170 pairs of either rule was
reported):

| Per species, among the 20 picked | Rule | p0 | p10 | p25 | p50 | p75 | p90 | p100 |
|---|---|---|---|---|---|---|---|---|
| minimum pairwise ANI | quality only | 93.77 | 95.69 | 96.77 | 97.59 | 98.44 | 99.09 | 99.98 |
| minimum pairwise ANI | ANI rule | 93.26 | 95.18 | 96.03 | 97.12 | 98.09 | 98.80 | 99.97 |
| median pairwise ANI | quality only | 95.69 | 97.42 | 98.01 | 98.61 | 99.03 | 99.80 | 100.00 |
| median pairwise ANI | ANI rule | 95.55 | 97.04 | 97.68 | 98.36 | 98.86 | 99.50 | 100.00 |
| near-clone pairs (of 190) | quality only | 0 | 1 | 3 | 10 | 30 | 90 | 190 |
| near-clone pairs (of 190) | ANI rule | 0 | 0 | 0 | 0 | 3 | 26 | 190 |

- **The quality-only picks were not diverse in the sense that matters most:** 230 of the 243 species had at least
  one pair of near-clones (ANI > 99.9) among their 20, 6,930 of 46,170 pairs in all (15%); in half the species 10
  or more of the 190 pairs. Under the ANI rule: 79 species, 2,895 pairs (6.3%).
- The near-clones left under the ANI rule are in species that are clonal as a whole: in 20 species even the most
  distant pair of the 20 ANI picks is at 99% or more (3 over 99.9: *Chlamydia muridarum*, *Bordetella holmesii*,
  *Francisella orientalis*, all 190 pairs near-clones under both rules), where no diverse 20 exist.
- In spread the change is modest: per species the minimum pairwise ANI moved by a median of −0.12 points (10th
  percentile −1.52, extreme −4.58); it fell in 158 species, stayed in 71 and rose in 14 (greedy farthest-first
  maximises each pick's distance to the earlier picks, not the minimum over all pairs). In the middle half of the
  species the ANI picks span down to 96.0 to 98.1, so 20 quality-ranked genomes already reach much of a species'
  range; what they add is repeats of one strain.

**Pairs unreported or under the 15% rule:** none, in any species (`table_unreported_empty` and
`table_unreported_af` are 0 in every row; `table_counts`: 0 empty, 0 below `--min-af`).

**The 26 genomes not served** (in 9 species: *Xanthomonas euvesicatoria* 8 of 46, *X. cissicola* 7 of 122, *E.
coli* 4 of 4,007, *Streptococcus dysgalactiae* 2 of 45, one each in *Klebsiella variicola*, *Vibrio cholerae*,
*V. harveyi*, *Acinetobacter nosocomialis* (22, so 21 left to choose from), *X. arboricola*). They were left out
of the pick; 8 of them were in their species' quality top 20, so the quality-only list could not have been
fetched as it stood. Every species still had more than 20 with a sequence, so each kept 20; what is lost is
whatever diversity the 26 carried (not measurable without their sequence).

**The new remainder** (selected minus phase 1; *measured*, `selection.json`): **13,485 genomes, 56,213,661,902
bp** (52,154,809 proteins by NCBI's counts), against C-25's 13,485 and 55,910,526,156 bp: the same count, 303 Mbp
(0.54%) more. The whole S2C selection: 25,458 genomes, 106,712,605,589 bp (quality only: 106,409,469,843). Of the
remainder, 4,617 genomes (17,541,096,714 bp) are in the 243 species and already fetched (section 6); the other
8,868 genomes, 38,672,565,188 bp (*derived*: the differences), are in species with at most 20 and are still to be
fetched.

## 6. The archive

`python -I c27_fetch.py archive out_ani/remainder.tsv archive/s2c_phase2_genomes.tar genomes/B01..B12/lock.json
--where-nonempty ani_pick --prefix s2c_phase2` (unit `c27-archive`, in the image, 09:59:22 to 10:01:07, 51.4 CPU-s,
*measured*): the remainder's rows with an `ani_pick` (picked in a species over the cap), each file checked against
its lock's `gz_sha256` first. **4,617 genomes** (243 × 19), 17,541,096,714 bp, 5,125,639,925 bytes of `.fna.gz`;
4,630 members: `s2c_phase2/<batch>/<acc>.fna.gz`, each batch's `lock.json` cut to the kept genomes, and
`s2c_phase2/kept.tsv`; GNU tar, members sorted, mode 0644, mtime 0, owner 0 and no names.

- Tar: **5,134,489,600 bytes, sha256 `af3ae7c348328e78e7091e7790ba20e3103a1f0bef21eef1be18d36c2e36fcb9`**.
- Copied (unit `c27-archive-copy`, host `rsync`, 10:01:19 to 10:02:00) to
  `/stor/scratch/Wilke/if2453/all-the-data/provisional/s2c_phase2/.s2c_phase2_genomes.tar.tmp`; its sha256 read
  back there equalled the above; renamed with `mv -n` to `s2c_phase2_genomes.tar`. The directory is new; nothing
  else on `/stor` was written, overwritten or deleted. Recorded in `sources.yaml` (`stored_s2c_phase2`).
- The 23,840 fetched genomes not kept (28,457 − 4,617) stay in `/tmp/if2453/c27/genomes/` until the PI says they
  may be removed; the local copy of the tar stays in `/tmp/if2453/c27/archive/`.

## Pod use

- Units, all `systemd-run --user`, `nice 19`, `MemorySwapMax=0`, inside the image: `c27-probe` (scope, 2 CPUs,
  seconds), `c27-paths`, `c27-mplan`, `c27-mplan2`, `c27-skani-plan` (scopes, 1 CPU, seconds), `c27-measure-t1`
  (1 CPU, 30 s), `c27-measure-t2` (2 CPUs, 15 s), `c27-skani-gate` (a shell loop), `c27-skani-0..3` (1 CPU each,
  09:25 to 09:58), `c27-select-ani` (1 CPU, 40.1 CPU-s), `c27-afsum`, `c27-archive`, `c27-archive-copy` (1 CPU).
  At most 4 of our threads at once, all by day; the account's total stayed at or under 40 by the gate's sum.
- Test runs: `pod_pytest.sh C-27-4` under `/tmp/if2453/ci/C-27-4`, with `BIOBOX_POD_PYDEPS` set to the tree's own
  `src` (an empty value falls back to `/tmp/if2453/pilot_s1/pydeps`: `${BIOBOX_POD_PYDEPS:-...}` substitutes for an
  empty value too), so only the image's locked environment and the tree were on the path. The whole suite ran as 4
  pytest processes from 10:02, after the skani run, when the account was at 34 to 35 CPUs.
- Read only: phase 1's genomes and locks under `/tmp/if2453/r50/genomes` and `/tmp/if2453/ladder/genomes`, the GTDB
  metadata and phase 1's list under `/tmp/if2453/c25`. No network request was made.

## Open questions for the PI

- Whether to keep the ANI selection as it is: it removes most near-clone pairs (230 species to 79; 6,930 pairs to
  2,895) at the cost of a modest change in spread (median minimum pairwise ANI 97.59 to 97.12) and 303 Mbp more
  (0.54%); phase 1 stays a subset.
- The alignment-fraction rule never acted (smallest larger fraction 56.47%); it can stay as stated or be dropped.
- When the 23,840 genomes not kept (about 29.7 GB of `.fna.gz`, *derived*: 34,788,621,074 − 5,125,639,925 bytes)
  and the 13.6 GB of sketches may be removed from `/tmp`.
- Whether the 26 not served should be retried later (nine species; 8 were in their quality top 20).
