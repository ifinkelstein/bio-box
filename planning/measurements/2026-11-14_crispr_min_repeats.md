# CRISPR arrays at a minimum of two repeats against three (card E-13)

Date: 2026-11-14, wilkcomp02. Every number carries a label: *measured* (host, date, command or record), *derived*
(formula and inputs) or *estimate*. Small tables are in `2026-11-14_crispr_min_repeats/` beside this note. No default
was changed: the producer of every release stays `minced-0.4.2-w8-nr3-r23x47-s26x50`; which minimum a release uses is
the PI's.

## Summary

Every number of this summary is *measured* (wilkcomp02, 2026-11-14, the commands of the setup) or *derived* from
those, with its source and formula in the section it comes from.

- At `-minNR 2` MinCED does not return "the arrays of 3 plus some arrays of 2". It returns a different and much larger
  set: on P1 + P2 (10.61 Gbp) 34,061 arrays against 2,577 at 3; 31,413 of them overlap no array found at 3 (the
  *extra calls*), and only 13,990 of those have two repeats: 17,423 have three or more (up to 816), which the default
  run did not report. 210 of the 2,577 arrays found at 3 come back at 2 with other coordinates, 13 with no overlapping
  array at all (*measured*, below). So "a minimum of 2 everywhere" would replace the current calls, not extend them.
- The extra calls are mostly not CRISPR arrays: 42.7% lie wholly inside a called gene (5.8% of the arrays found at 3),
  40.2% are in genomes with no Cas domain at all (4.4%), 1.10% lie within 10 genes of a Cas gene (68.1%), and only
  208 of 31,413 (0.66%) carry a repeat like one of their own genome's arrays of 3. Near Cas genes they are enriched
  only 1.4 to 1.6 times over random placement in their own genome; the arrays found at 3 are enriched 33 times.
- The shuffled null (dinucleotides kept in 1 kb windows, repeats destroyed) gives 6 calls at 2 and none at 3 on
  12.28 Gbp of S1, P1 and P2, all six on P2 and none of two repeats: 0.57 per Gbp of P1 + P2 against 2,960 extra
  calls per Gbp in the real sequence. Chance repeats of the genomes' composition are not where the extra calls come
  from; real repeats are, above all repeated and low-complexity genes.
- The worked example gains what the card is about: of its 20 Tn7-like loci (TniQ with Cas7f within 10 genes), 7 have
  an array within 10 genes at 3 and 14 at 2; the 7 that gain are six 2-repeat arrays and one 3-repeat array, each
  117 to 133 bases from Cas6f, with type I-F-like repeats, among them the control's (`NZ_AP026781.1`, 422557-422647).
- Proposed rule (for the PI): keep the default producer at 3 everywhere and add a second producer at 2 whose calls
  are kept only when they overlap no array of 3, lie within 1 kb of a gene with a Cas domain, and lie not wholly
  inside a called gene (rule A5 below): it adds 7.35 arrays per Gbp of P1 + P2 (78), of which at most 2.96 per Gbp
  (31.4) are false by the random-placement bound, and it recovers all 7 Tn7-like gains. Two repeats everywhere
  (rule B) would add 2,960 per Gbp, of which about 2,926 per Gbp lie far from any Cas gene.

## Setup

- Code: branch `task/E-13` (this card), mirrored to `/tmp/if2453/ci/E-13` by `scripts/pod_pytest.sh E-13 ...`; the
  producer takes the minimum as `run(..., min_repeats=k)` (step 1 of the card), the default unchanged.
- Image: the runtime image of `faec1f4`, `/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif`, digest
  `sha256:d67300aa1321b571fa8da79e74dfccdd87a8dca2e56518f324902efa6b9f098a` (RUNBOOK §0; checked by `pod_run.sh
  --digest` at every launch). MinCED 0.4.2.
- Data: read-only. The data root `/tmp/if2453/r10/root` was bound read-only into every unit. The scratch root
  `/tmp/if2453/e13/view` was made as card E-12 made its own (`python scripts/stress/e13_crispr_min_repeats.py view
  /tmp/if2453/r10/root /tmp/if2453/e13/view S1,P1,P2`: a symlink per stage directory of `contig_index`,
  `seq_locator`, `seq_shards`, `genome_input`, and `store_params.json`). The producers' files, ledgers and run records
  went into the view. Nothing was written under the root.
- Profile: `profiles/wilkcomp02.yaml` with `tmp: /tmp/if2453/e13/tmp`.
- Launches: one driver unit (`e13-driver`, then `e13-driver2` to resume after a fix of the analysis) ran, one at a
  time, the units below; each at 4 threads, `nice 19`, a memory cap (`top` before each launch: under 2 CPUs of the
  account in use, cap 56 on a Saturday). Two trial analyses on the P1-only null ran by hand in a scope beside the
  null (4 threads, 20 s, before it began; and 2 threads, 52 s, while it ran: 6 threads in all for that minute).
  - `e13-nr3-{p1,p2,s1}`: the default producer exactly as card E-12 ran it: `pod_run.sh --profile ... --image ...
    --digest ... --bind /tmp/if2453 --bind /tmp/if2453/r10/root:/tmp/if2453/r10/root:ro --env
    PYTHONPATH=/tmp/if2453/ci/E-13/src --mem 16G --threads 4 e13-nr3-p1 run-step --profile ... --root
    /tmp/if2453/e13/view --stage P1 --step crispr-minced-0.4.2-w8-nr3-r23x47-s26x50 --threads 4 -- python -m
    biobox.cli stage crispr --root /tmp/if2453/e13/view --stage P1 --cpus 4 --tmp ... --seq-stages P1,S1,P2`.
  - `e13-nr2-{p1,p2,s1}`: the same with `--step crispr-minced-0.4.2-w8-nr2-r23x47-s26x50` and the command `python
    scripts/stress/e13_crispr_min_repeats.py produce --root /tmp/if2453/e13/view --stage P1 --cpus 4 --min-repeats 2
    --tmp ... --seq-stages P1,S1,P2` (it calls the producer with `min_repeats=2` and writes into `$BIOBOX_STAGING`;
    `stage.py` and `cli.py` are not this card's).
  - `e13-null`: `systemd-run --user --wait --pipe --collect --unit=e13-null -p MemoryMax=16G -p MemorySwapMax=0 -p
    CPUQuota=400% nice -n 19 singularity exec ... python scripts/stress/e13_crispr_min_repeats.py null --root
    /tmp/if2453/r10/root --stages S1,P1,P2 --cpus 4 --tmp /tmp/if2453/e13/scratch-null --out /tmp/if2453/e13/null`
    (a P1-only trial `e13-null-p1` first).
  - `e13-analyse`: the same prefix (24G) with `analyse --root /tmp/if2453/r10/root --release rN.3 --view
    /tmp/if2453/e13/view --null /tmp/if2453/e13/null --out /tmp/if2453/e13/out --tmp ... --threads 4`. It reads
    `rN.3`'s catalog (opened read-only from a copy in scratch) and the view's producer files.

## The default producer is unchanged (byte identity with card E-12)

The nine files of the default producer written here (feature, `l1_crispr` and ledger of S1, P1 and P2) have the
sha256 of the files card E-12 wrote on 2026-11-13 with the code before this card (`sha256sum` of
`/tmp/if2453/e12/view/...` and `/tmp/if2453/e13/view/...`, *measured*, wilkcomp02, 2026-11-14); for example P2's
`l1_crispr` part `d54c1530...31ea` and feature part `57985603...07a4` in both.

## Cost (*measured*, run records `runs/stage=<S>/crispr-<producer>.json` of the view)

| | P1 nr3 | P1 nr2 | P2 nr3 | P2 nr2 | S1 nr3 | S1 nr2 |
|---|---|---|---|---|---|---|
| CPU-minutes (`cpu_min`) | 3.028 | 3.415 | 27.235 | 30.800 | 5.166 | 5.805 |
| wall seconds (`wall_s`) | 49.6 | 54.8 | 414.4 | 467.5 | 90.7 | 96.9 |
| peak RSS of one process (MB, `peak_rss_bytes`) | 339.8 | 344.5 | 439.8 | 510.1 | 540.6 | 540.5 |

At 2 the step costs 1.12 to 1.13 times what it costs at 3 (*derived*: 3.415 / 3.028, 30.800 / 27.235, 5.805 /
5.166), within 20% of the card's expectation (the same order as 3.1, 26.6 and 5.1, card E-12). The default's rerun
took 0.98, 1.02 and 1.01 times E-12's CPU-minutes.

The null (unit `e13-null`, 4 worker processes) took 160.7 CPU-minutes and 40.3 minutes of wall time, memory peak 2.1
GB (systemd's "CPU time consumed: 2h 40min 42.468s", "Service runtime: 40min 21.391s", *measured*); of it 102.0
CPU-minutes were the shuffle in Python (`shuffle_cpu_s` 6,118.7) and the rest MinCED, twice. The analysis
(`e13-analyse`) took 1.1 CPU-minutes and 47 s (systemd), reading `rN.3`'s 12,048,866 genes.

## Arrays at 2 and at 3 (*measured*, `analyse`, `summary.json` → `per_stage`)

| | P1 | P2 | S1 |
|---|---|---|---|
| bases | 1,016,625,803 | 9,596,180,328 | 1,664,109,938 |
| arrays at 3 | 286 | 2,291 | 213 |
| arrays at 2 | 3,379 | 30,682 | 3,613 |
| of them with 2 repeats | 1,293 | 12,757 | 2,116 |
| arrays of 3 found at 2 with the same coordinates and repeats | 267 | 2,100 | 200 |
| arrays of 3 found at 2 with other coordinates | 19 | 191 | 13 |
| of these, with no overlapping array at 2 | 1 | 12 | 1 |
| arrays at 2 that overlap an array of 3 but differ from it | 27 | 254 | 19 |
| **extra calls** (arrays at 2 overlapping no array of 3) | 3,085 | 28,328 | 3,394 |
| genomes with an extra call | 231 of 250 | 2,045 of 2,247 | 2,592 of 27,258 |

Repeats of the extra calls (P1 + P2, *measured*): 2 in 13,990, 3 in 6,684, 4 in 3,652, 5 to 9 in 5,482, 10 or more
in 1,605; the largest has 816 repeats of `GGGGGTGTGGGGGGGGTGTGGTG` over 59.8 kb inside a called gene of a P2 genome.
The largest ones are mostly of this kind, low-complexity or tandem-repeat protein genes: all twelve largest (80 to
816 repeats) overlap a called gene and seven lie wholly inside one; one of them (87 identical 37-base repeats, P2,
[670838, 677181)) lies within 10 genes of a Cas gene and looks like a real array the default run missed. Why MinCED at 2 reports arrays of 3 or more that it does not report at
3 was not investigated (it would need MinCED's source); the counts above are what it does. The spacer-length limit
`-maxSL 50` is not a hard limit in MinCED's output: the first spacer of an extra call reaches 224 bases and of an
array at 3, 109 (*measured*).

## The extra calls on P1 and P2, on the release's genes and domains

Definitions are those of the script's docstring. "Cas gene": a called gene (`protein_occurrence`) whose protein has a
domain of one of the 31 Pfam families of the CRISPR example's role list
(`examples/crispr_cas_neighbourhoods/sql/q01_loci_by_domains.sql`, every role whose name starts with `cas`: cas1
PF01867; cas2 PF09827, PF09707; cas3 PF21384, PF18019, PF21802, PF18395, PF22590, PF30455; cas4 PF01930; cas5
PF09704; cas5f PF09614; cas6 PF09559, PF01881, PF10040; cas6f PF09618; cas7 PF05107; cas7f PF09615; cas8b PF24121,
PF24122; cas8c PF09709; cas8e PF09481; cas8f PF09611; cas9 PF16593, PF16592, PF16595, PF22702; cas10 PF22335,
PF12469; cas11e PF09485; cas12a PF21918) in `l1_domains`, library `pfam-38.2` (the only one of `rN.3`). `rN.3` has
12,048,866 genes on S1, P1 and P2, 7,628 of them Cas genes; 139 of 250 P1 genomes and 1,227 of 2,247 P2 genomes
have one (*measured*). "Within 10 genes": at most 9 called genes lie wholly between the array and the Cas gene;
"within 1 gene": none does (the next gene, or one the array overlaps); "within 10 kb", "within 1 kb": the gap in
bases. Circular contigs are read across the origin. The expected counts are the second null (section below).

*Measured* (wilkcomp02, 2026-11-14, `analyse`, `summary.json` → `extra_p1_p2`, `extra_of_2_repeats_p1_p2`,
`extra_of_3_or_more_p1_p2`, `nr3_p1_p2`); "expected" rows are the second null (*derived*, Σ f).

| | extra calls, P1 + P2 | of 2 repeats | of 3 or more | arrays at 3 (reference) |
|---|---|---|---|---|
| calls | 31,413 | 13,990 | 17,423 | 2,577 |
| within 10 genes of a Cas gene | 346 (1.10%) | 162 | 184 | 1,756 (68.1%) |
| expected by random placement | 248.1 | 121.7 | 126.4 | 53.3 |
| within 10 kb of a Cas gene | 308 | 155 | 153 | 1,724 |
| expected by random placement | 215.3 | 105.8 | 109.5 | 47.8 |
| within 1 gene of a Cas gene | 112 | 56 | 56 | 1,279 |
| expected by random placement | 70.1 | 33.9 | 36.2 | 19.2 |
| within 1 kb of a Cas gene | 94 | 45 | 49 | 1,287 |
| expected by random placement | 57.4 | 27.6 | 29.8 | 16.5 |
| wholly inside a called gene | 13,416 (42.7%) | 4,727 | 8,689 | 150 (5.8%) |
| overlapping a called gene | 16,783 | 6,244 | 10,539 | 1,217 |
| in a genome with no Cas gene at all | 12,624 (40.2%) | 5,135 | 7,489 | 113 (4.4%) |
| repeat like one of the genome's arrays of 3 (identity ≥ 0.8) | 208 | 82 | 126 | (all) |
| within 10 genes and not inside a gene | 242 | 132 | 110 | 1,740 |
| expected by random placement | 136.4 | 80.5 | 55.9 | 51.5 |
| within 1 kb and not inside a gene | 78 | 38 | 40 | 1,276 |
| expected by random placement | 31.4 | 18.2 | 13.1 | 16.1 |

By stage (*measured*, `extra_by_stage`): P1 3,085 extra calls, 45 within 10 genes of a Cas gene (32.4 expected), 33
within 10 kb (27.2), 1,383 inside a gene, 1,170 in genomes without a Cas gene; P2 28,328, 301 (215.7), 275 (188.1),
12,033 and 11,454.

Reading. Observed over expected near Cas genes is 1.39 at 10 genes, 1.43 at 10 kb, 1.60 at 1 gene and 1.64 at 1 kb
for the extra calls (*derived*: 346 / 248.1, 308 / 215.3, 112 / 70.1, 94 / 57.4), against 32.9 at 10 genes for the
arrays found at 3 (1,756 / 53.3). Most extra calls behave like repeats placed with no regard to Cas genes; a small
excess near them (about 98 at 10 genes, 37 at 1 kb, *derived*: observed − expected) is where the real minimal
arrays are.

## Repeat length, spacer length and repeat identity (*measured*, `summary.json` → `distributions`)

Per array: the length of its first repeat, of its first spacer, and the identity of its first two repeats (matching
bases over the length when the two have one length, else 1 - Levenshtein distance / longer length); for arrays of 3
or more also the mean identity of every repeat to the column-majority consensus. Quantiles min / p10 / p25 / median
/ p75 / p90 / max.

| | extra calls, P1 + P2 (31,413) | extra calls within 10 genes of a Cas gene (346) | arrays at 3, P1 + P2 (2,577) |
|---|---|---|---|
| repeat length | 23 / 23 / 25 / 28 / 35 / 41 / 47 | 23 / 23 / 25 / 28 / 34 / 38 / 47 | 23 / 26 / 28 / 30 / 36 / 37 / 47 |
| spacer length | 26 / 26 / 26 / 34 / 47 / 55 / 224 | 26 / 26 / 30 / 36 / 43 / 51 / 104 | 26 / 30 / 32 / 34 / 37 / 40 / 109 |
| identity of the first two repeats | 0.38 / 0.77 / 0.92 / 1.00 / 1.00 / 1.00 / 1.00 | 0.67 / 0.78 / 0.96 / 1.00 / 1.00 / 1.00 / 1.00 | 0.39 / 0.75 / 0.88 / 1.00 / 1.00 / 1.00 / 1.00 |
| first two identical | 21,269 (67.7%) | 250 (72.3%) | 1,459 (56.6%) |
| identity in [0.9, 1) / [0.8, 0.9) / < 0.8 | 2,820 / 3,308 / 4,016 | 29 / 30 / 37 | 402 / 324 / 392 |
| consensus identity (arrays of 3 or more) | n 17,423: 0.51 / 0.87 / 0.90 / 0.93 / 0.97 / 0.99 / 1.00 | n 184: 0.74 / 0.89 / 0.91 / 0.94 / 0.97 / 0.99 / 1.00 | n 2,577: 0.63 / 0.94 / 0.97 / 0.99 / 1.00 / 1.00 / 1.00 |

Reading. The identity of the first two repeats does not separate the extra calls from real arrays: identical pairs
are more common among the extra calls (67.7%) than among the arrays at 3 (56.6%), whose first repeat is often a
degenerate end repeat. A filter on it (rule B1) removes a third of the extra calls without making the rest look like
arrays. What does separate them is the spacer length (a quarter of the extra calls sit at MinCED's minimum of 26, and
a tenth above 55; the arrays at 3 lie between 30 and 40 in 80% of cases), the consensus identity of the longer calls
(median 0.93 against 0.99), the position inside genes and the absence of Cas genes. Spacer length is a lead for a
later filter, not measured as a rule here. (The quantiles are of the "lower" kind: a p25 of 26 means at least a
quarter of the extra calls have a first spacer of 26 bases, MinCED's minimum; at least a tenth have one of 55 or
more.)

## The worked example's Tn7-like loci (*measured*, `tn7_loci.tsv`)

A locus: a gene with a TniQ domain (PF06527) with a Cas7f gene (PF09615) within 10 genes of it by gene rank on its
contig. 20 loci, as the example reports on `rN.3` (all in P1 and P2). An array is "at" a locus when it lies within 10
genes of the TniQ gene.

| | loci |
|---|---|
| TniQ loci with a Cas7f within 10 genes | 20 |
| with an array at 3 | 7 |
| with an array at 2 | 14 |
| that gain an array at 2 and have none at 3 | 7 |

The seven loci that gain (coordinates 0-based half-open on the contig; distance to the nearest Cas6f gene in bases;
all seven are within 1 gene and within 1 kb of a Cas gene and not inside a gene, so rules A4 and A5 keep them; none
has a repeat like one of its genome's arrays of 3):

| stage | contig | TniQ | gained array | repeats | repeat | to Cas6f |
|---|---|---|---|---|---|---|
| P1 | NZ_CP058972.1 | 1203829-1205005 | 1208806-1208894 | 2 | GTGTACTGCCGAGTAGGCAGCTGAGAAA | 122 |
| P2 | CP063663.1 | 325295-326525 | 330436-330520 | 2 | GTGACCTGCCGTATAGGCAGCTGA | 130 |
| P2 | NZ_AP025476.1 | 1066620-1067832 | 1071707-1071789 | 2 | CCTGCCGCATAGGCAGCTAAGAA | 133 |
| P2 | NZ_AP025507.1 | 871322-872516 | 867537-867625 | 2 | ATTATCAGCTACCTATTCGGCAGTTTAC | 127 |
| P2 | NZ_AP026781.1 (control) | 426340-427525 | 422557-422647 | 2 | TTATTATCAGCTACCTATTCGGCAGTTCAC | 128 |
| P2 | NZ_CP011034.1 | 2303561-2304737 | 2299618-2299762 | 3 | CCTCTGCCTTTTCGGCAGTTAAA | 117 |
| P2 | NZ_CP066149.1 | 2902819-2904034 | 2907924-2908011 | 2 | GTGTCCTGCCGAATAGGCAGCTAAGAA | 129 |

The six 2-repeat arrays carry repeats of the type I-F kind: four of the form `...CTGCCG..TAGGCAGCT...`, two
(`NZ_AP025507.1` and the control) `...ATCAGCTACCTATTCGGCAGTT.AC`, the family of the control's 5-repeat array on its
chromosome 1 (card E-12); the gained 3-repeat array of `NZ_CP011034.1` (`CCTCTGCCTTTTCGGCAGTTAAA`, of the second
family) is one the default run did not report although it has three repeats. The control (`NZ_AP026781.1`) gains the array card
E-12 found by hand with `minced -minNR 2` (1-based 422558-422647, here 0-based [422557, 422647)). Six loci have no
array at either minimum.

## The two nulls

### Null 1: shuffled contigs (*measured*, `null_summary.json`)

Every contig of S1, P1 and P2 (the same contigs and bases as above), its sequence shuffled within non-overlapping
windows of 1,000 bases (the last shorter) by the Altschul-Erickson doublet shuffle, seeded per contig from seed 13
and its contig_uid: each window keeps its multiset of dinucleotides and its first and last base exactly, so its base
composition and the contig's dinucleotide counts are kept exactly (checked for every window as it was made, and by
`tests/r1/test_crispr_plant.py`). MinCED with the producer's settings at 2 and at 3 on the shuffled records (linear
records; no flanks for circular contigs).

| | S1 | P1 | P2 | P1 + P2 |
|---|---|---|---|---|
| contigs | 27,258 | 539 | 5,018 | 5,557 |
| bases | 1,664,109,938 | 1,016,625,803 | 9,596,180,328 | 10,612,806,131 |
| arrays at 3 | 0 | 0 | 0 | 0 |
| arrays at 2 | 0 | 0 | 6 | 6 |
| arrays at 2 per Gbp (*derived*) | 0 | 0 | 0.63 | 0.57 |
| real extra calls per Gbp, for comparison (*derived*) | 2,039.5 | 3,034.5 | 2,952.0 | 2,959.9 |

The six chance calls (`null_extra_arrays.tsv` on the pod) all have 3 or more repeats of 23 to 47 bases with first
spacers of 37 to 659 bases, identity of the first two repeats 0.65 to 1.00; all six overlap a called gene's position
and none lies within 10 genes or 10 kb of a Cas gene (0.06 expected by the second null). So the chance rate of the
genomes' composition is at most 0.0002 of the extra calls (*derived*: 0.57 / 2,959.9). The shuffle keeps composition
but not the long-range structure of genes (repeated domains, tandem repeats, paralogues), which is where the extra
calls are; it is a null for "two repeats by chance", not for "repeats that are not CRISPR".

### Null 2: real repeats, placed at random in their own genome (*derived*)

For an array of length L in genome g, f = Σ_c |N_c(L)| / Σ_c (n_c − L + 1) over g's contigs c with n_c ≥ L, where
n_c is the contig's length and N_c(L) the set of start positions s in [0, n_c − L] at which [s, s + L) lies in the
window beside a Cas gene of c (an interval per Cas gene: (b_K − L + 1, e_K − 1) for "within K genes", where b_K is
the K-th largest start among the genes ending at or before the Cas gene's start and e_K the K-th smallest end among
those starting at or after its end; [C.start − B − L, C.end + B] for "within B bases"; their union). The expected
number of calls in the window, were they placed at random in their genomes, is Σ f over the calls
(`placement_share` in the script). It keeps each call's real length and genome and so the genome's Cas density;
it does not keep a preference for intergenic or coding DNA (the calls wholly inside genes are counted apart).

## Proposed rule, with its measured error rates

Each rule's false calls are bounded from the second null: were every call that passes the rule's other filters false
and placed at random in its genome, Σ f of them would fall in the window, so the false calls a windowed rule keeps
are at most that sum (*derived*, an upper bound; the true calls kept are at least adds − bound). For the rules
without a window the estimate of the false calls is the kept calls far from any Cas gene (neither within 10 genes nor
10 kb): it assumes that a true array of 2 repeats sits beside Cas genes (orphan arrays exist, so it overstates).
"Matches an array of 3" counts kept calls whose repeat has identity ≥ 0.8 (either strand) with the repeat of an array
found at 3 in the same genome: an independent sign of a real CRISPR repeat, which misses true minimal arrays whose
repeat family has no array of 3 in the genome. Per Gbp: of P1 + P2, 10.613 Gbp (*measured*, sum of contig lengths).

*Measured* counts (`rules.tsv`), *derived* bounds and rates.

| rule | what it keeps (beside every array of the default producer) | adds | adds per Gbp | wrongly adds (estimate) | per Gbp | basis | inside a gene | repeat like an array of 3 | Tn7-like gains kept (of 7) | shuffled-null calls kept |
|---|---|---|---|---|---|---|---|---|---|---|
| A | extra calls within 10 genes of a Cas gene | 346 | 32.60 | ≤ 248.1 | ≤ 23.38 | random placement | 104 | 41 | 7 | 0 |
| A1 | as A, not wholly inside a called gene | 242 | 22.80 | ≤ 136.4 | ≤ 12.85 | random placement | 0 | 38 | 7 | 0 |
| A2 | extra calls within 10 kb of a Cas gene | 308 | 29.02 | ≤ 215.3 | ≤ 20.29 | random placement | 85 | 40 | 7 | 0 |
| A3 | as A, the first two repeats identical | 250 | 23.56 | ≤ 176.9 | ≤ 16.67 | random placement | 62 | 30 | 6 | 0 |
| A4 | extra calls next to a Cas gene (1 gene), not inside a gene | 90 | 8.48 | ≤ 38.3 | ≤ 3.61 | random placement | 0 | 18 | 7 | 0 |
| **A5** | **extra calls within 1 kb of a Cas gene, not inside a gene** | **78** | **7.35** | **≤ 31.4** | **≤ 2.96** | random placement | 0 | 17 | 7 | 0 |
| B | every extra call (a minimum of 2 everywhere) | 31,413 | 2,959.91 | 31,053 | 2,925.99 | far from Cas | 13,416 | 208 | 7 | 6 (0.57 per Gbp) |
| B1 | as B, the first two repeats identical | 21,269 | 2,004.09 | 21,005 | 1,979.21 | far from Cas | 8,084 | 149 | 6 | 1 |
| B2 | as B, not wholly inside a called gene | 17,997 | 1,695.78 | 17,747 | 1,672.23 | far from Cas | 0 | 172 | 7 | 2 |

The Tn7-like column is from `tn7_loci.tsv` and the gained arrays' rows of `extra_arrays.tsv` (all seven pass every
window and are not inside a gene; the six 2-repeat gains have identical repeats and first spacers of 30 to 36 bases,
the 3-repeat gain of `NZ_CP011034.1` a first-pair identity of 0.74, so A3 and B1 lose it; *measured*). Rule B (and B1, B2) also changes the default's own calls if it replaces
the producer: 210 of its 2,577 arrays on P1 + P2 come back with other coordinates and 13 vanish (section "Arrays at 2
and at 3"); as a second producer beside the default it does not.

Recommendation (the PI decides; no default was changed):

1. Do not lower the minimum to 2 everywhere (rules B, B1, B2): it multiplies the calls by 13 (34,061 against 2,577),
   changes the coordinates of 8% of the current arrays, and nearly everything it adds lies far from Cas genes,
   inside genes or in genomes without Cas genes.
2. If minimal arrays are wanted (the CRISPR example's question), add a second producer at 2 and keep from it only the
   calls that overlap no array of the default producer, lie within 1 kb of a gene with a Cas domain and not wholly
   inside a called gene (rule A5): +7.35 arrays per Gbp, at most 2.96 per Gbp false, all 7 Tn7-like gains kept. The
   next-gene window (A4) is about the same (+8.48, at most 3.61 false per Gbp). The 10-gene window (A1) keeps more
   true arrays (at least 106, against 47 for A5) but also more false ones (at most 12.85 per Gbp, 56% of what it
   adds).
3. The filter needs the release's gene calls and Pfam scan, so it is a step after both (a later card: in
   `build-stage` after the domains, or in the release), and its output is a `feature` producer of its own (for
   example `minced-0.4.2-w8-nr2-r23x47-s26x50-cas1kb`), never mixed with the default's partitions.

## Files

`2026-11-14_crispr_min_repeats/` (copied from `/tmp/if2453/e13/out` and `/tmp/if2453/e13/null` on the pod):
`summary.json` (every count of this note), `rules.tsv`, `tn7_loci.tsv`, `null_summary.json`, `extra_near_cas.tsv`
(the extra calls within 10 genes or 10 kb of a Cas gene, with their annotations). Kept on the pod only (too large
for the repository): `extra_arrays.tsv` (all 31,413 extra calls) and `null_extra_arrays.tsv`, sha256 in
`SHA256SUMS`.
