# Context-free gene calling (policy p3) on stage P1, 2026-11-12 (card C-21f)

Card C-21f asks how gene calls change if they stop depending on the genome (policy p3 of
`planning/design/C-21_shared_contigs.md` §2b, A3). The answer feeds the PI's decision on DEC-28.

**No truth set was used.** The card allows reference annotations only through a `sources.yaml` entry, and there is
none, so nothing was fetched. Every agreement figure below compares p3 with **today's stored calls (p2)**. It is not
accuracy: where p2 and p3 differ, this note cannot say which one is right. As a weak hint only, §4 uses the stage's
stored Pfam hits.

## Summary

- **Chromosomes: 99.78% of genes identical** (857,267 of 859,139). The design note expected more than 99% (*estimate*),
  so the criterion "chromosome calls agree to within about 1%" holds in aggregate. But 23 of the 250 chromosomes are
  below 99%, and the lowest is at 94.9%.
- **Other contigs of 100 kb or more** (secondary chromosomes and megaplasmids, trained on themselves under p3):
  **96.39% identical**. Agreement grows with length: the median contig is at 91.0% for 100 to 200 kb and 98.3% for
  1 to 3 Mb. 113 of 113 contigs change.
- **Contigs under 100 kb** (meta under p3, genome-trained single under p2): **86.14% identical**, 10.49% start-only,
  3.38% lost, 6.14% gained.
- **All of P1**: 99.47% identical, 0.40% start-only, 0.13% lost, 0.19% gained. Genes go from 927,920 to 928,470.
- **Lost genes are mostly short and have no domain.** 19.3% of them carry a stored Pfam hit, against 81.6% of all
  stored genes. 237 lost genes with a Pfam hit are 0.026% of the stage's genes.
- **Cost: p3 is +1.1% CPU on P1**, because only 0.545% of P1's bases are in contigs under 100 kb. Meta costs 12.2
  times genome-trained single per base when finding genes alone, and 4.1 times single with its training spread over
  the bases. Both figures miss the note's *derived* 5.8 by more than 20% (§5).
- **The whole measurement cost 28.0 CPU-minutes**, 53% under the card's "about 1 CPU-hour" (*estimate*), and within
  9% of the note's own derivation, 2 × 15.39 = 30.8 (§5).

## 1. Setup

Every number in §2 to §5 is *measured* on wilkcomp02 on 2026-11-12, from 07:57 to 08:15 CST, unless it is labelled
*derived*. The pod's load average was 32 to 37 from other jobs. The account's CPU use was 3,318% to 3,793% at each
launch, under the 4,000% day cap (`ps -u if2453 -o pcpu=`, summed by `launch.sh`).

- **Image**: `/stor/home/if2453/bio-box/images/bio-box-test_74b3e05.sif` with pyrodigal 3.7.1 and pyrodigal-gv 0.3.2
  (`envcheck.py`, unit `c21f-check`). The tree was `plan/task-cards` at `c3dc901` plus the scripts below, rsynced to
  `/tmp/if2453/ci/c21f`.
- **Units**: each job ran as `systemd-run --user --unit c21f-<name> -p MemoryMax=16G -p MemorySwapMax=0
  -p CPUQuota=200%`, at `nice -n 19`, via `singularity exec --cleanenv`. The data root `/tmp/if2453/r10/root` was
  bound read-only, and `envcheck.py` confirmed from inside the image that it is not writable. All output went to
  `/tmp/if2453/c21f`. At most two units ran at once.
- **Input**: stage P1 of `/tmp/if2453/r10/root`, read-only. It holds 250 genomes (GTDB R226, complete, one per
  species), 539 contigs and 1,016,625,803 bp. Every genome has kind `prokaryote`, every genome is at least 100 kb,
  and no contig is marked circular in `genome_input`. P1 therefore has no viral, plasmid-kind or under-100-kb genome
  and no origin-wrapping gene, so this measurement covers none of those cases. Under p3, viral and small genomes
  would be called as today anyway.
- **Stored calls**: `l1/feature/producer=pyrodigal-3.7.1-gv-0.3.2-p2/stage=P1` (927,920 genes), written on
  2026-11-10 by run record `runs/stage=P1/genecall.json` (15.39 CPU-min, 16 workers).

Scripts, in `scripts/stress/c21f_context_free/`:

- `recall.py ROOT P1 OUT {p2|p3}` re-calls every genome, one at a time, read-only. p2 goes through the product's own
  `genecall._finder`. p3 uses gv-meta for viral genomes; for any other contig of 100 kb or more, `GeneFinder(meta=False)`
  trained on that contig alone with the genome's code; otherwise meta. Genes are found by the product's
  `genecall._genes`. It records per-contig thread CPU (`time.thread_time`) for training and for finding.
- `compare.py` compares genes per contig class: identical (start, end, strand); start-only (same strand and stop);
  lost; gained. It checks that each side has one gene per (contig, strand, stop), so the four categories add up.
- `detail.py` gives lengths, the stored Pfam-38.2 hits of lost genes, start shifts, and per-contig agreement by
  length bin.
- `launch.sh` starts each unit, checking the STOP files, the two-unit limit and the account's CPU; `envcheck.py` is
  the image check above.

Raw records are in `planning/measurements/2026-11-12_context_free_genecall/`: `compare.json`, `detail.json`, the
`timed.py` records `c21f-*.json`, and the comparison log.

## 2. Control: the p2 re-call equals the stored calls

| class | contigs | stored genes | re-called | identical | start only | lost | gained |
|---|---|---|---|---|---|---|---|
| all | 539 | 927,920 | 927,920 | 927,920 (100%) | 0 | 0 | 0 |

The environment reproduces the stored calls gene for gene, so every difference in §3 comes from the policy alone.

## 3. p3 against the stored calls

Contig classes, as defined in `compare.py`:

- `chromosome`: the longest contig of a genome, itself at least 100 kb;
- `other_ge100k`: every other contig of at least 100 kb;
- `lt100k_linear`: contigs under 100 kb. All of them are linear in P1.

Percentages are of the stored genes in the class.

| class | contigs | changed contigs | bases | stored genes | p3 genes | identical | start only | lost | gained |
|---|---|---|---|---|---|---|---|---|---|
| chromosome | 250 | 79 | 939,743,485 | 859,139 | 858,861 | 857,267 (99.78%) | 1,375 (0.16%) | 497 (0.06%) | 219 (0.03%) |
| other_ge100k | 113 | 113 | 71,341,752 | 63,097 | 63,768 | 60,821 (96.39%) | 1,738 (2.75%) | 538 (0.85%) | 1,209 (1.92%) |
| lt100k_linear | 176 | 155 | 5,540,566 | 5,684 | 5,841 | 4,896 (86.14%) | 596 (10.49%) | 192 (3.38%) | 349 (6.14%) |
| all | 539 | 347 | 1,016,625,803 | 927,920 | 928,470 | 922,984 (99.47%) | 3,709 (0.40%) | 1,227 (0.13%) | 1,777 (0.19%) |

Per contig, the share of a contig's stored genes that are identical under p3, by length:

| class | length | contigs | minimum | median | contigs under 99% |
|---|---|---|---|---|---|
| chromosome | 200 to 500 kb | 4 | 94.86% | 98.91% | 2 |
| chromosome | 500 kb to 1 Mb | 11 | 99.15% | 100% | 0 |
| chromosome | 1 to 3 Mb | 78 | 98.30% | 100% | 5 |
| chromosome | 3 Mb or more | 157 | 97.34% | 100% | 16 |
| other_ge100k | 100 to 200 kb | 50 | 70.16% | 90.96% | 49 |
| other_ge100k | 200 to 500 kb | 27 | 84.55% | 94.89% | 26 |
| other_ge100k | 500 kb to 1 Mb | 15 | 92.41% | 96.73% | 15 |
| other_ge100k | 1 to 3 Mb | 18 | 93.56% | 98.34% | 14 |
| other_ge100k | 3 Mb or more | 3 | 98.24% | 98.95% | 2 |
| lt100k_linear | under 10 kb | 45 | 0% | 88.89% | 28 |
| lt100k_linear | 10 to 100 kb | 131 | 0% | 88.89% | 117 |

Under p2, a chromosome's training set also includes its genome's other replicons. Dropping them leaves 171 of 250
chromosomes unchanged. A contig of 100 kb to 1 Mb trained alone differs from the genome-trained call on a median of
3% to 9% of its genes, depending on length, and on up to 30%. Whether the self-trained or the genome-trained call is closer to the truth cannot be decided here.

## 4. What changes: lengths and stored Pfam hits

Pfam-38.2 hits are read from `l1/domains/library=pfam-38.2/stage=P1`, at the gathering cutoff. That partition has no
rows for P1 proteins that S1 already held, so the domain shares below are lower bounds.

| class | stored genes: median nt, share with a Pfam hit | lost: n, median nt, with a Pfam hit | gained: n, median nt |
|---|---|---|---|
| chromosome | 825, 82.0% | 497, 237, 73 (14.7%) | 219, 198 |
| other_ge100k | 834, 77.9% | 538, 225, 103 (19.1%) | 1,209, 243 |
| lt100k_linear | 676.5, 68.9% | 192, 211.5, 61 (31.8%) | 349, 219 |
| all | 825, 81.6% | 1,227, 228, 237 (19.3%) | 1,777, 237 |

Start-only changes. "p3 shorter" means p3's start lies downstream of the stored start. "Cuts a stored domain" means
the shift in codons passes the first residue of the protein's first stored Pfam envelope.

| class | start-only | p3 shorter | p3 longer | median shift (nt) | p3 shorter and cuts a stored domain |
|---|---|---|---|---|---|
| chromosome | 1,375 | 939 | 436 | 45 | 175 |
| other_ge100k | 1,738 | 710 | 1,028 | 45 | 143 |
| lt100k_linear | 596 | 225 | 371 | 36 | 57 |
| all | 3,709 | 1,874 | 1,835 | 45 | 375 |

Both the lost and the gained genes are short: median 228 and 237 nt, against 825 nt for all genes. Lost genes carry
a Pfam hit a quarter as often as genes overall (*derived*: 19.3 / 81.6). These are mostly marginal ORFs, where the
two models disagree, not well-supported genes. Without a truth set, it is open which side is right about the 237 lost
genes with a hit and the 375 start changes that cut a stored domain.

## 5. CPU cost

Thread CPU seconds of training and finding, from `recall.py`, summed by class:

| class | bases | p2: train + find (s) | p3: train + find (s) | p3 / p2 |
|---|---|---|---|---|
| chromosome | 939,743,485 | 458.6 + 332.0 | 421.5 + 331.1 | 0.95 |
| other_ge100k | 71,341,752 | 0 + 25.3 | 29.9 + 25.4 | 2.18 |
| lt100k_linear | 5,540,566 | 0 + 1.5 | 0 + 18.2 (meta) | 12.18 |
| all | 1,016,625,803 | 458.6 + 358.8 = 817.4 | 451.4 + 374.7 = 826.1 | 1.011 |

Under p2, a genome's training time is booked on its chromosome's row.

Whole processes (`timed.py` rusage, one thread each):

| unit | wall | CPU-min | peak RSS |
|---|---|---|---|
| c21f-p2 | 831.2 s | 13.85 | 615 MiB |
| c21f-p3 | 839.8 s | 13.98 | 586 MiB |
| c21f-check, c21f-compare, c21f-detail | 2.9, 4.5, 4.1 s | 0.05, 0.08, 0.07 | 248, 339, 456 MiB |

**Derived from the tables above:**

- p3 against p2 on P1: +1.06% (826.1 / 817.4 thread-s) and +0.99% (13.98 / 13.85 process CPU-min).
- Share of P1's bases under 100 kb: 0.545% (5,540,566 / 1,016,625,803). This share is what decides the meta-mode cost
  of p3. On a stage with many plasmids or draft contigs it would be larger. Open item 4 of the design note (the share
  on the 50 Gbp root) is not measured here.
- Meta against single per base, on the same 176 contigs: 3.285 s/Mbp (18.20 / 5.541) against 0.270 s/Mbp
  (1.494 / 5.541) for finding with a genome-trained model, 12.2 times.
- Meta against single with training included: P1's whole p2 call costs 0.804 s/Mbp (817.4 / 1,016.6), so meta is
  4.09 times that.
- **Both miss the note's 5.8 (*derived*) by more than 20% (AGENTS §6): +110% and −30%.** The note compared different
  stages and modes (gv-meta on S1 against single on P2, whole-step CPU). The two measured figures bracket it, and
  neither matches it. At P1's 0.545% of bases under 100 kb, the difference changes p3's total cost by under 1%.
- Total CPU of this card: 28.04 CPU-min (0.05 + 13.85 + 13.98 + 0.08 + 0.07).
  - **Against the card's "about 1 CPU-hour" (*estimate*): −53%, a miss of more than 20%.** The estimate overstates its
    own basis: the note derives it as two re-calls of P1 at 15.39 CPU-min, which is 30.8 CPU-min, not an hour.
  - Against that 30.8: −9%.
  - Each re-call ran on one thread in 13.9 CPU-min, against 15.39 CPU-min for the stored 16-worker call, which also
    wrote the tables.

## 6. What this does not settle

- **Which calls are right.** p3's agreement with p2 is not accuracy. Judging p3 needs a truth set, as the design
  note's §2b suggests: genes matched by stop codon against RefSeq CDS for non-chromosome contigs. That needs a PI entry
  in `sources.yaml`.
- **Contigs P1 lacks.** Circular contigs, origin-wrapping genes, plasmid-kind genomes and genomes under 100 kb do not
  occur in P1, so these numbers say nothing about them.
- **The 100 kb threshold for self-training.** Agreement on 100 to 200 kb contigs trained alone has a median of 91%.
  If the PI prefers fewer changed calls on mid-size replicons, a higher threshold for self-training, with meta below
  it, is one option. This note does not measure it.
