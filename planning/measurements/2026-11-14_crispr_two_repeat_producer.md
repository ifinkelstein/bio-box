# The producer of two-repeat CRISPR arrays beside Cas genes (card E-15, DEC-40)

Date: 2026-11-14, wilkcomp02. Every number carries a label: *measured* (host, date, command or record), *derived*
(formula and inputs) or *estimate*.

## What it is

A second CRISPR array producer, `minced-0.4.2-w8-nr2-r23x47-s26x50-cas1kb-caspfam1-pfam-38.2`
(`biobox.l1.crispr.run_near_cas`; `bio-box stage crispr --near-cas --library pfam --library-version 38.2`; in
`build-stage --crispr` after the stage's Pfam scan when a `pfam` library is given). MinCED runs exactly as for the
default producer but at `-minNR 2`; a call is kept only when

- (a) it shares no base with an array of the default producer (`minced-0.4.2-w8-nr3-r23x47-s26x50`) on its contig,
- (b) the gap max(0, gs - e, s - ge) between its span [s, e) and some Cas gene's [gs, ge) is at most 1,000 bases, a
  Cas gene being a gene of the stage's gene calls (`feature_protein`) whose protein has a domain of one of the 31 Pfam
  families of `src/biobox/l1/cas_families/cas-pfam-1.tsv` in the stage's `pfam-38.2` scan (`model_acc` without its
  version),
- (c) no gene of the stage has gs <= s and e <= ge,

each read on a circular contig at the shifts -n, 0 and +n. The rejected calls are counted in the ledger by the first
rule they fail, in the order (a), (b), (c). The family list is the CRISPR worked example's `cas` roles, the list card
E-13 measured with; frozen by its sha256, so another list is another producer. No schema changed.

## Setup

- Code: branch `task/E-15`, mirrored to `/tmp/if2453/e15/tree` (the producer's source byte-identical to the branch's,
  `sha256sum` of `crispr.py`, `cli.py`, `stage.py`, `provenance.py` on both sides).
- Image: the runtime image of `faec1f4`, `/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif`, digest
  `sha256:d67300aa1321b571fa8da79e74dfccdd87a8dca2e56518f324902efa6b9f098a` (checked by `pod_run.sh --digest` and by
  the post-run script's `sha256sum`).
- Data: the data root `/tmp/if2453/r10/root`, bound read-only. A scratch view `/tmp/if2453/e15/view` holds a symlink
  per stage directory of what the producers read (`contig_index`, `seq_locator`, `seq_shards`, `genome_input`,
  `store_params.json`, `feature_protein` of `pyrodigal-3.7.1-gv-0.3.2-p2`, `l1_domains` of `pfam-38.2`) for S1, P1
  and P2 (`/tmp/if2453/e15/bin/mkview.py`); the producers' files, ledgers and run records went into it.
- Launches: one driver (`/tmp/if2453/e15/bin/run_all.sh`) ran, one at a time, the units `e15-nr3-{p1,p2,s1}` (the
  default producer, as cards E-12 and E-13 ran it) and then `e15-near-{p1,p2,s1}`: `pod_run.sh --profile ... --image
  ... --digest ... --bind /tmp/if2453 --bind /tmp/if2453/r10/root:/tmp/if2453/r10/root:ro --env
  PYTHONPATH=/tmp/if2453/e15/tree/src --mem 16G --threads 3 e15-near-p1 run-step --profile ... --root
  /tmp/if2453/e15/view --stage P1 --step crispr-minced-0.4.2-w8-nr2-r23x47-s26x50-cas1kb-caspfam1-pfam-38.2 --threads
  3 -- python -m biobox.cli stage crispr --near-cas --library pfam --library-version 38.2 --root /tmp/if2453/e15/view
  --stage P1 --cpus 3 --tmp ... --seq-stages P1,S1,P2`. Three threads each, `nice 19`, beside the 50 Gbp build.
- After them, one unit `e15-post` (MemoryMax 16G, CPUQuota 300%, nice 19) ran `/tmp/if2453/e15/bin/post.sh`: the
  default producer's bytes against card E-12's, the brute-force oracle, the scratch view of `rN.3` with both
  producers and `nbr` and the example's array checks on it. Nothing was written under the root or under
  `/tmp/if2453/e12`, `/tmp/if2453/e13`, `/tmp/if2453/ex2`.

## Calls and verdicts (*measured*, wilkcomp02, 2026-11-14, the run-step log lines `{"result": ...}` of `e15-near-*`)

| | P1 | P2 | S1 |
|---|---|---|---|
| bases | 1,016,625,803 | 9,596,180,328 | 1,664,109,938 |
| MinCED calls at 2 repeats | 3,379 | 30,682 | 3,613 |
| rejected: overlaps a default array (a) | 294 | 2,354 | 219 |
| rejected: no Cas gene within 1 kb (b) | 3,076 | 28,243 | 3,360 |
| rejected: wholly inside a gene (c) | 3 | 13 | 0 |
| rejected: longer than the circle | 0 | 0 | 0 |
| **kept** | **6** | **72** | **34** |
| kept per Gbp (*derived*, kept / bases) | 5.90 | 7.50 | 20.43 |

The calls at 2 and the overlaps are card E-13's (3,379 and 30,682 calls; 3,379 - 3,085 = 294 and 30,682 - 28,328 =
2,354 overlapping an array of 3, *derived* from its table). Rule (c) rejects few because rule (b) comes first: most of
E-13's 13,416 extra calls inside genes are far from Cas genes. On P1 + P2 the producer keeps 78 on 10.61 Gbp, 7.35 per
Gbp (*derived*: 78 / 10.613), card E-13's rule A5. S1 (phages) was not part of E-13's rule measurement: its 34 kept
calls, 20.4 per Gbp, have no measured false-positive bound (an open question below).

## Oracle: the kept calls equal card E-13's rule A5 set (*measured*, `/tmp/if2453/e15/bin/e15_oracle.py`, unit `e15-post`)

A brute-force filter written apart from the producer: plain Python loops, row group by row group, over card E-13's
stored calls at 2 repeats (`/tmp/if2453/e13/view`, `minced-0.4.2-w8-nr2-r23x47-s26x50`, read only), its arrays at 3,
the root's gene calls and Pfam scan of P1 and P2, contig topology from `genome_input` and lengths from the calls'
`l1_crispr` rows, with the 31 families read from the example's `q01_loci_by_domains.sql` (not from the producer's
list). Three sets of (stage, contig_uid, start, end) compared:

| | P1 | P2 | together |
|---|---|---|---|
| this oracle | 6 | 72 | 78 |
| the producer's `feature` rows | 6 | 72 | 78 |
| card E-13's rule A5 (`2026-11-14_crispr_min_repeats/extra_near_cas.tsv`, `near_1_kb = 1`, `inside_gene = 0`) | 6 | 72 | 78 |

All three sets are equal; no call on one side only. The oracle also re-derived E-13's intermediate counts: 3,085 and
28,328 extra calls, 1,383 and 12,033 of them inside a gene, 927,920 and 8,722,859 genes, 1,552,653 and 14,608,733 domain
rows on P1 and P2. Cost 97.4 s wall, 97.2 CPU-s, 159 MB (the script's own record).

## The default producer is unchanged (*measured*, `/tmp/if2453/e15/bin/cmp_default.sh`)

The nine files of the default producer written in this view (feature, `l1_crispr` and ledger of S1, P1 and P2) have
the sha256 of card E-12's (`/tmp/if2453/e12/view`): `default producer: the 9 files are byte-identical to card E-12's`;
for example P2's `l1_crispr` part `d54c1530...31ea` and feature part `57985603...07a4`. The filtered step reads them
and writes only its own partitions (also tested on constructed contigs and in `build-stage`).

## Cost (*measured*, run records `runs/stage=<S>/crispr-<producer>.json` of the view)

| | P1 default | P1 filtered | P2 default | P2 filtered | S1 default | S1 filtered |
|---|---|---|---|---|---|---|
| CPU-minutes (`cpu_min`) | 3.109 | 3.355 | 26.723 | 31.790 | 5.670 | 6.479 |
| wall seconds (`wall_s`) | 66.2 | 71.7 | 539.6 | 647.8 | 122.9 | 138.6 |
| peak RSS of one process (MB, `peak_rss_bytes`) | 339.8 | 431.1 | 429.8 | 1,018.7 | 536.6 | 594.3 |

The filtered step costs 1.08, 1.19 and 1.14 times the default step (*derived*: 3.355 / 3.109, 31.790 / 26.723, 6.479 /
5.670), against 1.12 to 1.13 for MinCED at 2 alone (card E-13): the filter's reads and joins add about 5% on P2
(*derived*: 1.190 / 1.131, E-13's P2 ratio 30.800 / 27.235).
The estimate's factor `CRISPR_NEAR_CAS_COST` in `src/biobox/stage.py` is the largest ratio, 1.19 (the first draft's
1.13, E-13's ratio, was 5% under P2's: within 20%, so no stop); with it the step is refused above 66.5 Gbp in one
stage (*derived*: 0.8 × 300 × 60 / (182.0 × 1.19)), against 9.6 Gbp for P2. The 1.02 GB peak on P2 is the filter's
DuckDB connection (one thread, `DUCK_MEMORY_LIMIT` 2 GB, spilling to scratch) holding the 30,682 calls and the genes of
their contigs; MinCED's JVMs have the fixed 512 MiB heap. systemd's account of the P1 unit: MemoryPeak 1,521,586,176
bytes (the unit's cgroup, page cache included), CPUUsageNSec 206,630,909,000 (`pod_run: e15-near-p1 usage`).

## A scratch view of rN.3 with both producers (*measured*, `/tmp/if2453/e15/bin/mkview_rel.py`, unit `e15-post`)

`/tmp/if2453/e15/relview`, made as card E-14's view: symlinks to the root's directories, the root's tables of `l1`
and `l1/feature`'s producers, both `minced-` producers' `feature` and `l1_crispr` files and crispr ledgers of S1, P1
and P2 from `/tmp/if2453/e15/view`; its own `manifests/release=rN.3.json` (the root's with 18 entries added) and
catalog. C-08 (`provenance.check_provenance`) on that manifest: no problem over its 78 files, the six real crispr
ledgers included. 11.2 s wall, 379 MB (`/usr/bin/time -v`).

**The control Tn7-like element** (`bio-box nbr <uid> --release rN.3 --window 10g --format tsv`, uids from `bio-box
find`, `/tmp/if2453/e15/bin/e15_real.py`; all on `NZ_AP026781.1`):

| Locus | Window 10g (bases) | Arrays within it |
|---|---|---|
| TniQ 426340-427525 | 413722-444027 | 422557-422647, 2 repeats, 1 spacer, repeat `TTATTATCAGCTACCTATTCGGCAGTTCAC`, the filtered producer, `delta_bp` 3,693, `side_rel` 1, not partial |
| Cas7f 423377-424433 | 411868-440339 | the same array, `delta_bp` 730 |
| complete I-F Cas7f 687857-688841 | 676009-697441 | as card E-14: 679543-681134 (27 repeats) and 689534-690762 (21 repeats), the default producer |

The control's array is the one card E-12 found by hand at a minimum of 2 and card E-13 lists (0-based
[422557, 422647), 128 bases from Cas6f).

**The example's two runs** (`/tmp/if2453/ex2/rN.3-b/runs`, read-only; `nbr <run> --window 10g --format tsv` and `nbr
summarize ... --format tsv` on the view):

| Run | loci | with an array | of the default producer | of the filtered producer | arrays (default, filtered) |
|---|---|---|---|---|---|
| TniQ | 635 | 22 | 13 | 9 | 13, 8 |
| Cas7f | 157 | 146 | 138 | 43 | 173, 45 |

`nbr summarize` writes one array row per producer with these counts (for example `crispr_array
minced-0.4.2-w8-nr2-r23x47-s26x50-cas1kb-caspfam1-pfam-38.2 635 9 0.0142 0 8` for the TniQ run).

**The example's array checks** (`bio-box sql --release rN.3 --runs ... sql/checks_arrays.sql` on the view):

| check | value | with the default producer only (card E-14) |
|---|---|---|
| `arrays_tniq_cas7f_loci_with_array` | 14 of the 20 I-F3 TniQ loci | 7 |
| `arrays_tniq_cas7f_loci_with_near_cas_array` | 7 | (new) |
| `arrays_tniq_control_window`, `arrays_cas7f_control_window` | 1, 1 | 0, 0 |
| `arrays_control_cas6f_nearest_bp` | 128 | 256,168 |
| `arrays_cas7f_tniq_loci_with_array` | 13 | 6 |
| `arrays_cas7f_complete_not_virus_loci_with_array` | 78 of 78 | 78 of 78 |
| `arrays_cas7f_virus_loci_with_array` | 46 | 45 |
| `arrays_pdam_complete_if_window`, `arrays_pdam_complete_if_window_has_array` | 2, true | 2, true |
| `arrays_icp1_cas7f_window` | 1 | 1 |
| `arrays_producers` | both names | the default's |

14 of the 20 Tn7-like TniQ loci have an array within the window: 7 at three repeats and 7 more from the filtered
producer, card E-13's count; the two producers' arrays are disjoint at these loci (7 + 7 = 14). Cost of the real-data
pass (`e15_real.py`, all commands in one process): 71.0 s wall, 157.5 CPU-s, 4.19 GB peak (`/usr/bin/time -v`), the
largest step `nbr summarize` of the TniQ run as in card E-14.

## Fixture tests (*measured*, wilkcomp02, 2026-11-14, `pytest -rA` in the test image of `faec1f4`, CPUQuota 300%)

- Constructed contigs (`tests/r2/test_e15_near_cas_producer.py`): 14 planted cases, 6 kept, 5 rejected far, 2
  inside a gene, 1 overlapping a three-repeat array; the kept set equals the test's own reading of the rules.
- Pruning (`test_the_filters_lookups_by_contig_read_the_row_groups_of_those_contigs`, 48 contigs, 8 rows per row
  group): `feature_protein` 16 of 60 row groups scanned, `l1_domains` 8 of 60, the default producer's `feature` 4 of
  4, each equal to the row groups whose footer range holds the looked-up keys (`e15 prune` line).
- Memory (`tests/r2/test_e15_near_cas_memory.py`, `scripts/stress/synth_scale.py` at 10 and 100 copies, 7,974,810
  and 79,748,100 bases): peak 339,083,264 and 339,070,976 bytes (`runstep.step_peak`, `child_maxrss`: the JVM), under
  the bound of 768 MiB with no growth. Limitation: the synthetic store's calls (26 and 262) all fail rule (b) or (a),
  so it keeps none and the write of kept rows is not loaded by it; the real P2 run (72 kept, 30,682 calls) is the
  larger load, 1.02 GB, under the DuckDB cap.

## Open

- The rate of kept calls on S1 (phages, 20.4 per Gbp, *derived*) is about three times that of P1 + P2; card E-13's
  false-positive bound was computed for P1 + P2 only. Whether phage arrays beside Cas genes are as often real (phages
  such as ICP1 carry CRISPR-Cas systems) is not measured.
