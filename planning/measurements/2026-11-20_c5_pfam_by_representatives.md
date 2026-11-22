# Pfam by NR90 representatives on stage C5: DEC-35 (6) stops on missed pairs (card C-47 part 2)

Date: 2026-11-20 23:05 CST to 2026-11-21 00:24 CST, wilkcomp02 (the pod's clock). Every number carries a label:
*measured* (host, date, command or file), *derived* (formula and inputs) or *estimate*.

**DEC-35 (6) tripped on missed pairs.** Inferring each member's Pfam 38.2 domains from its stage-local NR90
representative, under the rule as decided (DEC-32, DEC-35: a member is scanned when its length differs from its
representative's by more than 20%, or when a representative domain lies less than 80% inside the stored alignment),
misses **1.0721% of the true (protein, model) pairs** (72,286 of 6,742,372), above the 1% limit. Wrong pairs are
**0.9565% of reported pairs** (64,417 of 6,734,503), just under theirs (*measured*, wilkcomp02, 2026-11-21 00:13,
`c47_driver.sh evaluate`, `eval.json`). **Inference by stage-local NR90 representatives is not adopted. Cards B-47a and
B-47c to B-47e, and stages C1 to C4 and C6, are held for the PI.** This note gives the PI the result, where the
disagreements come from, and what alternative rules would cost, measured on the same data. It recommends nothing and
nothing in `src/` was changed.

Files beside this note in [`2026-11-20_c5_pfam_by_representatives/`](2026-11-20_c5_pfam_by_representatives/):
`eval.json` (the evaluation), `breakdown.json` (the breakdown and the sweep), `cm.files` and `aln.files` (the cluster
tables the two read, in the scratch output root on the pod), `units` (the clustering's step units).

## 1. How it was measured

- **Truth**: C5's full Pfam 38.2 scan, as built (card C-47 part 1; `planning/measurements/2026-11-20_stage_C5.md`):
  for each of C5's 3,972,512 distinct proteins its own rows, measured in C5 or copied from the stage that measured it
  (card B-47b; every copied row equals its source's, *measured* there). A pair is a distinct (protein, model).
- **Clusters**: C5's 3,972,512 distinct proteins clustered at NR90 by the release clustering's own code and parameters
  (`biobox.l2.cluster` in parts: `release cluster plan`, `step`, `finish`; MMseqs2 `--min-seq-id 0.9 -c 0.8 --cov-mode 1
  --cluster-mode 2`, the members aligned to their representatives), reading the data root bound read-only and writing
  under the scratch output root `/tmp/if2453/c47/eval/out` as release `c47eval5`, a name that exists only there. Nothing
  was defined or written in the data root. Result: 870,646 NR90 clusters (and 237,759 NR50, made because `finish`
  writes both levels; not used) (*measured*, `driver.log`, 2026-11-21 00:07). The plan's estimate was 645.9
  CPU-minutes for both levels (*estimate*, the code's `release_cpu_min`); the 60 step units (`c47-ev-001`, the plan, to
  `c47-ev-060`, `finish`; 8 threads each) used **248.29 CPU-minutes** in all, the largest 99.10 (`c47-ev-019`), peak
  memory 5,283,405,824 bytes (*measured*, the units' `CPUUsageNSec` and `MemoryPeak` in `driver.log`): **0.384 of the
  estimate** (*derived*). That is more than 20% under; it is reported here as DEC-53 asks of an under-run, not a stop.
  The clustering estimate is fitted on releases of one genome per species; a strain-rich stage's redundancy is removed
  cheaply by linclust before the expensive prefilter, which it does not model.
- **The rule evaluated** (`c47_nr90_pfam_eval.py eval`): a representative is scanned; a member whose length differs
  from its representative's by more than 20% is scanned ("length"); a member without an alignment, or with any
  representative domain whose envelope lies less than 0.8 inside [tstart, tend), is scanned ("unmapped"); every other
  member takes its representative's models ("inferred"). Right, wrong and missed are counted on inferred members; a
  scanned protein's pairs are its own, all right.
- **The breakdown** (`c47_nr90_pfam_eval.py breakdown`, 2026-11-21 00:21 to 00:24, unit `c47-breakdown`): the same
  tables, plus each model's gathering thresholds from `Pfam-A.hmm.gz` 38.2 (30,134 models with a `GA` line,
  *measured*) and each pair's best sequence and domain bit score. Its rule-as-built reproduces `eval.json`'s totals
  exactly (`CONSISTENT_WITH_EVAL: True`: scanned proteins and residues, true, right, wrong and missed pairs).
  Test-first on a hand-worked fixture: `tests/r1/test_c47_breakdown.py`.

## 2. The result (*measured*, `eval.json`)

| | proteins | residues |
|---|---|---|
| all of C5 | 3,972,512 | 1,350,839,880 |
| representatives | 870,646 | |
| scanned for length | 367,127 | |
| scanned, unmapped | 22,472 | |
| inferred | 2,712,267 | |
| **scanned in all** | **1,260,245 (31.72%)** | **334,284,889 (24.75%)** |
| scanned after B-47b's copy (proteins earlier stages measured are copied, not scanned) | 1,036,895 | 278,447,246 (20.61%) |

| pairs | count | rate |
|---|---|---|
| true (C5's full scan) | 6,742,372 | |
| reported (right + wrong) | 6,734,503 | |
| right | 6,670,086 | |
| **wrong** (inferred, not in the member's scan) | **64,417** | **0.9565% of reported** |
| **missed** (in the member's scan, not inferred) | **72,286** | **1.0721% of true: above 1%, DEC-35 (6) stops** |

Percentages *derived* from the counts. DEC-32's basis, for comparison (card 0-09b, on `rN.1b`, release-level NR90): 71%
of residues scanned, 0.4% of reported pairs wrong (*measured* there). C5 scans far fewer residues (24.75%) and has
more than twice the wrong rate.

- **By phylum**: every C5 genome is *Pseudomonadota* (Enterobacterales to *Vibrio*), so the by-phylum table has one
  row, equal to the totals.
- **B-47b's copy**: 703,241 of C5's own distinct proteins (204,247,443 residues) are in an earlier stage's scan and
  copied (*measured*). After both savings, 278,447,246 residues (20.61%) would be scanned.
- **f, the share of a representative domain inside the alignment** (rows on members that pass the length rule; rows,
  and in how many the member's own scan has the model): f = 1: 5,365,490, 5,301,944 (98.82%); [0.8, 1): 129,623,
  126,358 (97.48%); [0.5, 0.8): 14,076, 11,341 (80.57%); (0, 0.5): 5,165, 1,394 (26.99%); 0: 7,193, 938 (13.04%)
  (counts *measured*, percentages *derived*); the same shape as on `rN.3` (B-47 design note §3).

## 3. Where the disagreements come from (*measured*, `breakdown.json`)

### 3.1 Cause

- **Missed, 72,286**: in every one the representative has no row for the model at all (72,286). None comes from a
  representative row that was not mapped (0: a member with any representative domain under f 0.8 is scanned, DEC-35
  (4)). In 4,547 of them the representative has no domain at all.
- **Wrong, 64,417**: in every one the member's own scan lacks a model its representative has.

So all disagreements are models that one of two near-identical proteins carries above its gathering threshold and the
other does not; none is a failure of the coordinate mapping.

### 3.2 How far the hit sits from the gathering threshold

For each disagreeing pair, the protein that has the hit (the member for a missed pair, the representative for a
wrong pair), its best bit score minus the model's GA, by bin; and, for comparison, the same for the inferred pairs that
are right (the member's own hit). Counts *measured*; shares *derived* from them.

| bit score − GA | missed, domain | missed, sequence | wrong, domain | wrong, sequence | right inferred, domain | right inferred, sequence |
|---|---|---|---|---|---|---|
| [0, 1) | 30,910 | 8,167 | 29,291 | 7,460 | 70,383 | 14,479 |
| [1, 2) | 15,679 | 15,891 | 13,981 | 15,017 | 80,593 | 43,232 |
| [2, 5) | 15,070 | 25,121 | 14,017 | 22,784 | 231,087 | 186,206 |
| [5, 10) | 5,432 | 10,046 | 4,492 | 9,731 | 318,089 | 299,207 |
| [10, 25) | 3,411 | 9,433 | 2,057 | 7,353 | 711,411 | 714,307 |
| ≥ 25 | 1,784 | 3,628 | 579 | 2,072 | 3,652,040 | 3,806,172 |
| total | 72,286 | 72,286 | 64,417 | 64,417 | 5,063,603 | 5,063,603 |

- **Domain score within 2 bits of GA: 64.4% of missed and 67.2% of wrong pairs, against 3.0% of right inferred
  pairs**; within 5 bits: 85.3%, 88.9% and 7.5% (*derived*: (30,910 + 15,679) / 72,286; (29,291 + 13,981) / 64,417;
  (70,383 + 80,593) / 5,063,603; and with the [2, 5) bin). The reading holds: most disagreements are hits near the
  threshold, which two near-identical sequences fall either side of. The full scan calls each protein either way by
  the same threshold; this is the threshold's own sensitivity, not an error of the inference. About 1 in 7 disagreeing
  pairs (14.7% of missed, 11.1% of wrong) is 5 bits or more above GA (*derived*).

### 3.3 Identity to the representative, length difference, model

By the member's identity to its representative (MMseqs2's `pident` of the stored alignment; inferred members only):

| identity | members | true pairs | reported | wrong | missed | wrong rate | missed rate |
|---|---|---|---|---|---|---|---|
| < 0.90 | 59,163 | 110,239 | 109,759 | 2,078 | 2,558 | 1.89% | 2.32% |
| [0.90, 0.95) | 982,663 | 1,910,249 | 1,905,690 | 34,770 | 39,329 | 1.82% | 2.06% |
| [0.95, 0.97) | 431,660 | 784,790 | 783,187 | 10,464 | 12,067 | 1.34% | 1.54% |
| [0.97, 0.99) | 734,276 | 1,323,812 | 1,322,708 | 12,619 | 13,723 | 0.95% | 1.04% |
| [0.99, 1) | 490,294 | 987,110 | 987,070 | 4,401 | 4,441 | 0.45% | 0.45% |
| 1.0 | 14,211 | 19,689 | 19,606 | 85 | 168 | 0.43% | 0.85% |

(counts *measured*, rates *derived*). The 59,163 members under 0.90 identity are NR90 members whose clustering
identity (MMseqs2's own, over its alignment) passed 0.9 while the stored alignment's `pident` is lower; they are a
small share. **The error falls steeply with identity: members under 0.95 identity are 38% of inferred members and
carry 57% of wrong and 58% of missed pairs** (*derived*: (59,163 + 982,663) / 2,712,267; (2,078 + 34,770) / 64,417;
(2,558 + 39,329) / 72,286).

By length difference (|L_member − L_rep| / L_rep; inferred members only):

| length difference | members | wrong | missed | wrong rate | missed rate |
|---|---|---|---|---|---|
| 0 | 1,747,422 | 36,268 | 35,747 | 1.15% | 1.14% |
| (0, 0.05] | 704,968 | 20,936 | 25,078 | 1.38% | 1.64% |
| (0.05, 0.10] | 148,401 | 4,303 | 7,474 | 1.55% | 2.67% |
| (0.10, 0.15] | 71,909 | 1,968 | 2,547 | 1.63% | 2.10% |
| (0.15, 0.20] | 39,567 | 942 | 1,440 | 1.62% | 2.46% |

(rates *derived*: wrong / reported, missed / true in each bin; reported and true per bin in `breakdown.json`). Length
matters less than identity: equal-length members already have both rates above 1%.

By model: the disagreements are spread thinly. The 20 models with most missed pairs carry 15.8% of them, led by
PF13304 (1,117, 1.55%), PF00083 (938, 1.30%) and PF00155 (729, 1.01%); the 20 with most wrong pairs carry 14.9% of
them, led by PF00083 (926, 1.44%), PF13304 (717, 1.11%) and PF00890 (560, 0.87%) (*measured* counts and shares; the
sums *derived*; full lists in `breakdown.json`). No single family explains the rate.

### 3.4 The coordinate block, corrected

The evaluation's coordinate block (`eval.json`: 69.4% of a sample within 5 residues, p95 1,373) joined a
representative's domain row to every member row of the same model, so a model occurring twice in a protein paired
mismatched instances. The breakdown matches each mapped envelope to the member's own row of the model that overlaps it
most. On the first 200,000 representative rows of inferred members (of 5,436,875): 197,009 matched; 2,956 with no
overlapping instance of the model in the member (the wrong pairs among them); 35 unmappable (no aligned column in the
envelope). Of the matched, **94.3% have both ends within 5 residues of the member's own envelope and 96.4% within 10;
medians 0, p95 2 (start) and 3 (end)** (*measured*). The mapping itself is sound; `eval.json`'s figure was an artefact
of the join.

## 4. The options for the PI, with their measured cost (*measured*, `breakdown.json` `sweep`)

Every rule scans a member when its length differs by more than the length cut, when a representative domain lies less
than the f cut inside the alignment, or, if given, when its identity to the representative is under the identity cut.
Rates: wrong of reported, missed of true. Residues scanned of C5's 1,350,839,880; "after B-47b" leaves out proteins an
earlier stage measured. 48 rules were measured (4 length × 4 identity × 3 f); these are the ones the card names and
the cheapest that bring both rates to at most 1% and at most 0.5%.

| length cut | identity cut | f cut | wrong | missed | proteins scanned | residues scanned | after B-47b |
|---|---|---|---|---|---|---|---|
| **20%** | none | **0.8** (as decided) | **0.957%** | **1.072%** | 31.72% | **24.75%** | **20.61%** |
| 15% | none | 0.8 | 0.943% | 1.051% | 32.72% | 25.47% | 21.20% |
| 10% | none | 0.8 | 0.913% | 1.013% | 34.53% | 26.98% | 22.40% |
| 5% | none | 0.8 | 0.849% | 0.902% | 38.27% | 30.74% | 25.48% |
| 20% | none | 0.9 | 0.934% | 1.065% | 32.18% | 25.18% | 21.01% |
| 20% | none | 1.0 | 0.874% | 1.015% | 34.64% | 27.69% | 23.09% |
| 20% | 95% | 0.8 | 0.409% | 0.451% | 57.95% | 53.73% | 44.44% |
| 20% | 97% | 0.8 | 0.254% | 0.272% | 68.82% | 65.07% | 54.12% |
| 20% | 99% | 0.8 | 0.067% | 0.068% | 87.30% | 84.44% | 71.33% |
| **15%** | none | **1.0** | **0.866%** | **0.998%** | 35.35% | **28.20%** | **23.50%** | ← cheapest with both ≤ 1% |
| **20%** | **95%** | 0.8 | **0.409%** | **0.451%** | 57.95% | **53.73%** | **44.44%** | ← cheapest with both ≤ 0.5% |

(the rates *derived* from `breakdown.json`'s counts; "cheapest" is the rule with fewest residues scanned among the 48
that meet both levels, ties broken by the larger length cut). Other rules near the 1% line: 10% / none / 1.0 (0.848%,
0.967%; 29.33% of residues) and 5% / none / 0.8 (0.849%, 0.902%; 30.74%).

What the table says, without choosing:

- **Within the rule's own form (length and f), the 1% line is reached only just**: 15% length with f 1.0 brings missed
  to 0.998% at 28.2% of residues (+3.5 points over the rule as decided, *derived*). The margin is thin; on the next
  strain stage it may not hold.
- **An identity cut is what moves both rates**: scanning members under 95% identity halves both (0.41%, 0.45%) but
  scans 53.7% of residues (44.4% after B-47b), 29 points more than the rule as decided. 99% brings both to 0.07% but
  scans 84.4%: little saving over a full scan.
- **None of these removes the threshold effect** (§3.2): two thirds of the disagreements are hits within 2 bits of GA.
  A rule built on that (scan a member whenever its representative has a hit within a few bits of GA, or has no hit
  within a few bits below it) is a different design; it is not measurable from the full scan's rows, which keep only
  hits above GA, and was not measured here.
- **Release-level representatives (DEC-35 (7))** are the third option the decision rows name; their comparison is owed
  (§5).

## 5. DEC-35 (7): comparison with release-level representatives (`r26.11.2`, measured 2026-11-21)

Run once `r26.11.2` (`S1,P1,P2,P3,P4,P5,P6,H1,PL1,C5,V1`) published and verified its manifest (`46454d51b321845b9c7bffe15efaa165`), read-only on `/tmp/if2453/r10/root`, inside `bio-box_faec1f4.sif` (unit `c47-eval-rel`, 21:46 to 21:55 CST on 2026-11-21, 2 threads, `MemoryMax=16G`, `nice -n 19`, 502.5 s; output [`2026-11-20_c5_pfam_by_representatives/eval_release.json`](2026-11-20_c5_pfam_by_representatives/eval_release.json)), with `--release-cm` and `--release-aln` taken from `r26.11.2`'s manifest file list (never a glob):

```
python c47_nr90_pfam_eval.py eval --root /tmp/if2453/r10/root --stage C5 \
  --earlier S1,P1,P2,P3,P4,P5,P6,H1,PL1,V1 \
  --stage-cm $(cat /tmp/if2453/c47/eval/cm.files) --stage-aln $(cat /tmp/if2453/c47/eval/aln.files) \
  --release-cm $(cat /tmp/if2453/c47/eval/rel_cm.files) --release-aln $(cat /tmp/if2453/c47/eval/rel_aln.files) \
  --selection /tmp/if2453/c47/stages/C5/selection.tsv --objects /tmp/if2453/s2/inputs/C5/objects.json \
  --out /tmp/if2453/c47/eval/eval_release.json --tmp /tmp/if2453/c47/eval/eval_release_tmp --threads 2 --memory 12GB
```

| Clustering scope | Reps (*measured*) | Inferred members (*measured*) | Scanned proteins (after B-47b) (*measured*) | Scanned residues % (after B-47b) (*measured*) | Wrong pairs / reported (*measured*, % *derived*) | Missed pairs / true (*measured*, % *derived*) | DEC-35 (6) stop? |
|---|---:|---:|---:|---:|---:|---:|---|
| Stage-local NR90 (`C5` alone) | 870,646 | 2,712,267 | 1,260,245 (1,036,895) / 31.72% | **24.7464%** (**20.6129%**) | 64,417 / 6,734,503 (0.9565%) | 72,286 / 6,742,372 (**1.0721%**) | **yes** (missed > 1%) |
| Release-level NR90 (`r26.11.2`) | 714,266 | 2,843,695 | 1,128,817 (935,767) / 28.42% | **21.5519%** (**18.0404%**) | 69,431 / 6,731,947 (**1.0314%**) | 79,856 / 6,742,372 (**1.1844%**) | **yes** (both > 1%) |

- **`DEC-35 (7)` residue gap**: stage-local scans **`3.1945` percentage points** more residues than release-level (`24.7464%` vs `21.5519%`; **`2.5725` percentage points** after card B-47b's copy, `20.6129%` vs `18.0404%`, *measured* in `dec35_7_comparison`), **well below `DEC-35 (7)`'s 10-percentage-point revisit threshold** (`revisit_triggered_gt_10pp: false`).
- **Error rates under release-level NR90**: because `156,380` stage-local `C5` representatives become non-representative members of earlier-stage NR90 clusters in `r26.11.2` (`2,843,695` vs `2,712,267` inferred members), release-level representatives infer 3.19 percentage points more residues and **err more on both sides**: **`1.1844%` missed** (`79,856 / 6,742,372`) and **`1.0314%` wrong** (`69,431 / 6,731,947`), both exceeding the 1% `DEC-35 (6)` limit. Thus release-level representatives neither save >10 percentage points of scan work nor clear the 1% error gate (confirming DEC-62's choice of full scans for `C1` to `C4` and `C6`).

## 6. Pod use

| Unit | CST | What | Threads |
|---|---|---|---|
| `c47-evalplan`, `c47-ev-001` | 11-20 23:05 | `release cluster plan` (3,972,512 proteins into the MMseqs2 database) | 1 to 8 |
| `c47-evalcluster`, `c47-ev-002` to `c47-ev-060` | 11-20 23:06 to 11-21 00:07 | the clustering steps and `finish`, one unit each (none near the 280 CPU-minute stop; largest 99.10) | 8 |
| `c47-evaluate` | 11-21 00:12 to 00:14 | `c47_nr90_pfam_eval.py eval` (101.8 s) | 2 |
| `c47-breakdown` | 11-21 00:21 to 00:24 | `c47_nr90_pfam_eval.py breakdown` (191.7 s) | 2 |
| `c47-eval-rel` | 11-21 21:46 to 21:55 | `c47_nr90_pfam_eval.py eval` with `r26.11.2` release-level NR90 (`DEC-35 (7)`, 502.5 s) | 2 |

All beside the release train `c46-train-rel-1`, inside the account's cap of 56 threads. Nothing under the data root or
`/tmp/if2453/c46` was touched; nothing was written under `/stor`. Left on the pod: `/tmp/if2453/c47/eval/` (the
scratch output root with the cluster tables, the work directory, `eval.json`, `breakdown.json`, the logs).
