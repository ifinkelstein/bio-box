# The filtered CRISPR producer on phages: its chance rate on S1 (card E-16)

Date: 2026-11-15, wilkcomp02 (times from the pod's clock). Every number carries a label: *measured* (host, date,
command or record), *derived* (formula and inputs) or *estimate*. Small tables are in
`2026-11-15_crispr_near_cas_phages/` beside this note. A measurement only: no code under `src/` changed, no default
changed, the producer (`minced-0.4.2-w8-nr2-r23x47-s26x50-cas1kb-caspfam1-pfam-38.2`, DEC-40, card E-15) is as it was.

## Outcome

**The rule holds on phages as it is.** Under card E-13's null (each call at 2 repeats placed at random in its own
genome, the producer's three rules applied where it lands) the filter keeps 2.43 calls on S1 by chance, 1.46 per
Gbp (95% interval of the expectation 1.44 to 1.48 per Gbp), against the 34 it keeps, 20.4 per Gbp: 14.0 times
chance; and all 34 are one real CRISPR-Cas locus, the type I-F system of the *Vibrio* phage ICP1 and its kin, in 33
genomes: none lies beside a lone Cas4 or Cas3-domain gene. The phage rate is three times the prokaryote rate because
INPHARED holds 30 near-identical ICP1 isolates, not because phages attract chance calls: per Gbp the chance rate on
S1 is lower than on P1 + P2 (1.46 against 2.30 by the same null; 2.09 against 2.96 by E-13's bound).

## Summary

All *measured* on wilkcomp02, 2026-11-15, by the run of the setup (unit `e16-measure`, `summary.json`), or *derived*
from it where the formula is given.

- **Chance rate on S1** (1.664 Gbp; 3,613 calls at 2 repeats, 213 default arrays, 34 kept):

  | null | what is placed | expected kept by chance | per Gbp | Monte Carlo (10,000 draws), 95% interval of the mean, per Gbp | kept / expected |
  |---|---|---|---|---|---|
  | **N2** (the card's) | every call at 2, the three rules at the new place | **2.433** (sd 1.53) | **1.462** | 1.442 to 1.478 | **14.0** |
  | N3 | the 3,394 calls that pass rule (a) where they are | 1.006 (sd 0.99) | 0.605 | 0.586 to 0.609 | 33.8 |
  | E13 (card E-13's bound) | the 1,348 calls of N3 not inside a gene; only rule (b) at the new place | 3.474 (sd 1.80) | 2.087 | 2.079 to 2.121 | 9.8 (*derived*, 34 / 3.474) |

  Observed 34 = 20.43 per Gbp (*derived*, 34 / 1.664). P(K >= 34) under N2 is 4.8e-30 (exact Poisson-binomial); no
  draw of 10,000 reached 34 (the largest: 10). Excess over chance 31.6 calls (*derived*, 34 - 2.433).
- **The 34** (`kept_s1.tsv`): all in group (a), a phage with a Cas system (five Cas roles, cas1, cas3, cas5f, cas6f,
  cas7f, on five genes within 10 genes of the call); 0 in (b) beside a lone Cas4-like or Cas3-domain gene; 0 in (c).
  They are 31 calls in 30 ICP1 isolates (header "ICP1") and 3 in *Vibrio* phages JSF13, JSF14 and JSF17, which carry
  the same system (the worked example's "ICP1 isolates and the JSF phages"). 33 are the same array: 2 or 3 repeats
  (32 with 3, 1 with 2), 59 to 61 bases from the Cas3 gene (`Cas3-like_C_2`, PF22590); on 32 of their 33 contigs the
  default producer's array of the same repeat family lies 2,798 to 2,801 bases away. The 34th is an 8-repeat array
  87 bases from Cas6f on the 33rd contig (MW794172), where the default producer found no array at all.
- **Families**: the kept calls are beside `Cas3-like_C_2` (PF22590, 33) and `Cas_Csy2` (Cas5f, PF09614, 33) within
  1 kb, and `Cas_Csy4` (Cas6f, 1) and `Cas_Csy3` (Cas7f, 1). The commonest Cas family in phages, `Cas_Cas4` (PF01930:
  285 genes in 285 of the 27,258 phage genomes, 52% of S1's 544 Cas genes), attracts 0.43 calls by chance (0.26 per
  Gbp) and has 0 kept beside it.
- **Oracle checks** (both pass): (1) the three rules applied by this script to the 3,613 calls at 2 keep exactly the
  producer's 34 on S1 (and its 6 on P1 and 72 on P2), call for call, with the producer's rejection counts per rule;
  (2) card E-13's quantity on P1 + P2 is reproduced exactly: 31.3614 calls, 2.955 per Gbp (E-13: 31.36, 2.96), every
  one of its 17,997 per-call shares equal at six decimals, Monte Carlo 31.32 (31.21 to 31.42).
- **For scale**, the default producer's arrays (3 repeats) placed in the same way: on S1 49 of 213 lie within 1 kb
  of a Cas gene against 3.71 expected (13.2 times) and 49 within 10 genes against 8.92 (5.5 times); on P1 + P2 1,756
  within 10 genes against 53.34 (32.9 times, card E-13's number, reproduced).

## Setup

- Code: branch `measure/E-16` (from `plan/task-cards` at `933eb81`), mirrored to `/tmp/if2453/ci/E-16` by
  `scripts/pod_pytest.sh E-16 ...`. The script is `scripts/stress/e16_crispr_near_cas_phages.py` (new); card E-13's
  script is unchanged (its `Windows` and `placement_share` are imported for the 10-gene windows and as a second
  reading of the 1 kb share).
- Image: the test image of `faec1f4`, `/stor/home/if2453/bio-box/images/bio-box-test_faec1f4.sif` (the default of
  `scripts/pod_pytest.sh`), with the tree's `src` on `PYTHONPATH`.
- Data, all read-only, each table by the explicit path of its partition (no glob over `l1/feature` or `l1/crispr`,
  while card C-33's agent wrote new CRISPR producer files into the root):
  - from card E-13's view `/tmp/if2453/e13/view`: the calls at 2 (`l1/feature/producer=minced-0.4.2-w8-nr2-r23x47-s26x50/stage=S`
    and its `l1/crispr`), the default arrays (`...-nr3-...`);
  - from card E-15's view `/tmp/if2453/e15/view`: the kept calls (`l1/feature` and `l1/crispr` of
    `producer=minced-0.4.2-w8-nr2-r23x47-s26x50-cas1kb-caspfam1-pfam-38.2/stage=S`);
  - from the data root `/tmp/if2453/r10/root`: `l1/feature_protein/gene_call_version=pyrodigal-3.7.1-gv-0.3.2-p2/stage=S`,
    `l1/domains/library=pfam-38.2/stage=S`, `l0seq/genome_input/stage=S`, `l0seq/contig_index/stage=S`, and
    `l0seq/contig_occurrence/stage=S1/source=inphared` (record names), for S = S1, P1, P2;
  - card E-13's per-call table `/tmp/if2453/e13/out/extra_arrays.tsv`.
  The 31 Cas families and their roles are read from the worked example's `sql/q01_loci_by_domains.sql` (the script
  checks them equal to card E-13's `CAS_FAMILIES`); the producer's frozen list `cas-pfam-1.tsv` holds the same 31.
- Launch (`top` first: the other agent's MinCED processes at about 16 to 24 CPUs; this card 1 thread):
  `systemd-run --user --unit=e16-measure --collect -p MemoryMax=24G -p MemorySwapMax=0 -p CPUQuota=100% ... nice -n 19
  singularity exec --cleanenv --env OMP_NUM_THREADS=1 --env PYTHONPATH=/tmp/if2453/ci/E-16/src --bind /tmp/if2453
  bio-box-test_faec1f4.sif python /tmp/if2453/ci/E-16/scripts/stress/e16_crispr_near_cas_phages.py measure --root
  /tmp/if2453/r10/root --calls /tmp/if2453/e13/view --kept /tmp/if2453/e15/view --tree /tmp/if2453/ci/E-16 --out
  /tmp/if2453/e16/out --stages S1,P1,P2 --draws 10000 --e13-extra /tmp/if2453/e13/out/extra_arrays.tsv`, 02:06 to
  02:08 CST. Cost (*measured*, the script's record): 111.0 s wall, 110.9 CPU-s (S1 1.2, P1 12.6, P2 80.2 CPU-s for
  the null, the rest reading), peak RSS 1,197 MB. A rerun of S1 alone (`--out /tmp/if2453/e16/rerun_s1`) gave the
  same numbers to the last digit.
- Written only under `/tmp/if2453/e16/` and `/tmp/if2453/ci/E-16`, `/tmp/if2453/ci-tmp/E-16*` (tests). Nothing under
  the root, `/tmp/if2453/e12`, `e13`, `e15` or `/stor`.

## The null: what is randomized and what is held fixed

Randomized: the start of each call. Held fixed: each call's length L and its genome, so the number of calls per
genome and their lengths; the genome's contigs with their lengths and topology; its called genes; which genes are Cas
genes and of which families; the default producer's arrays. The start is uniform over every place where the call
fits in its genome: on a linear contig of n bases the n - L + 1 starts, on a circular one the n starts (the call
running across the origin), a contig weighted by its number of starts. No contig of S1, P1 or P2 is marked circular
in `genome_input` (*measured*: 0 of 27,258, 539 and 5,018), so the placement is linear everywhere and the circle
rule is exercised only by the test. Every S1 genome is one contig (27,258 contigs, 27,258 genomes, *measured*).

At the new place the producer's three rules are applied (its docstring; card E-15's note): (a) no base shared with a
default array of the contig, (b) a gap of at most 1,000 bases to a Cas gene (a called gene whose protein has a
domain of one of the 31 families in the stage's `pfam-38.2` scan), (c) not wholly inside a called gene.

The probability p_i that call i is kept is computed two ways: exactly from the geometry (in the space of starts, the
union of the Cas genes' windows [g_s - L - 1000, g_e + 1000] less the starts that share a base with a default array,
[d_s - L + 1, d_e - 1], and those wholly inside a gene, [g_s, g_e - L], as integer intervals), and by Monte Carlo
(10,000 places per call, the rules applied from the coordinates by the same function that reproduces the producer's
34, below). Each call has its own random stream (PCG64, seeded from 16 and its `feature_uid`), so nothing depends on
the order of the calls. A draw is one placement of every call; the number kept in a draw, K, is a sum of independent
Bernoulli(p_i): its expectation is the sum of p_i, its law a Poisson-binomial.

Three readings (the card asks for the first; the others bracket it):

- **N2**: every call at 2 repeats is placed, the three rules applied at the new place: the calls the producer would
  keep if no call had anything to do with Cas genes.
- **N3**: only the calls that pass rule (a) where they are (card E-13's "extra calls"; the others are default arrays
  called again) are placed, the three rules applied.
- **E13**: card E-13's quantity behind "at most 2.96 per Gbp": over the extra calls that are not wholly inside a gene
  where they are, the sum of f, the share of their genome's starts (every contig read as linear) within 1 kb of a
  Cas gene. Rules (a) and (c) are not applied at the new place, which is why it is an upper bound on what the rule
  keeps by chance (a call placed inside a gene would be rejected by (c); phage genomes are mostly coding).

Not kept by the null: any preference of real repeats for intergenic or coding DNA beyond what rule (c) does, and the
redundancy of the genomes (30 ICP1 isolates count as 30 genomes).

## Oracle checks

1. **The rules against the producer** (*measured*, `summary.json` → `stages.<S>.oracle_kept_equal`): applied to the
   real calls at 2 (card E-13's stored calls), the script's rules keep the producer's calls exactly, as sets of
   (contig_uid, start, end):

   | | S1 | P1 | P2 |
   |---|---|---|---|
   | calls at 2 | 3,613 | 3,379 | 30,682 |
   | rejected by (a) / (b) / (c), this script | 219 / 3,360 / 0 | 294 / 3,076 / 3 | 2,354 / 28,243 / 13 |
   | the producer's ledger (card E-15's table) | 219 / 3,360 / 0 | 294 / 3,076 / 3 | 2,354 / 28,243 / 13 |
   | kept, this script | 34 | 6 | 72 |
   | kept, the producer's `feature` rows | 34 | 6 | 72 |
   | on one side only | 0 | 0 | 0 |

   The Cas genes counted are 544, 822 and 6,262 (*measured*), 7,628 in all: card E-13's count from `rN.3`.
2. **Card E-13's bound reproduced** (*measured*, `summary.json` → `p1_p2`): "at most 2.96 per Gbp" is E-13's sum of f
   (null 2 of its note: the share of a call's genome's starts within 1 kb of a Cas gene, every contig linear) over
   the extra calls of P1 + P2 not wholly inside a gene. Here: the same 17,997 calls (none on one side only), each f
   equal to E-13's stored `f_1_kb` at six decimals (largest difference 0.0), sum 31.3614, 2.9551 per Gbp (*derived*,
   / 10.6128 Gbp), against E-13's 31.36 and 2.96. By Monte Carlo (the call placed, only rule (b) read): 31.32, 95%
   interval of the mean 31.21 to 31.42 (2.941 to 2.961 per Gbp), which holds the exact value.
3. **Exact against Monte Carlo** (*measured*): in all nine (stage, null) pairs and for P1 + P2 the exact expectation
   lies inside the Monte Carlo 95% interval of the mean (for example S1 N2: 2.4334 in 2.3999 to 2.4599).
4. **Constructed genomes** (`tests/r1/test_e16_near_cas_null.py`, 14 tests): a planted Cas4 gene, another gene and
   a default array at known places, the kept share counted by hand (1,500 of 149,802 starts on two linear contigs,
   2,200 of 10,000 on a circular one), equal to the exact function; the Monte Carlo within four standard errors of it;
   E-13's window share by hand (3,101 of 149,802; 2,502 of 9,901) equal to this script's and to E-13's own
   `placement_share`; the rules at seven hand-placed spans; the Poisson-binomial tail against enumeration. Seen
   failing first (14 errors: the script did not exist).

## The chance rate, by stage (*measured*, `summary.json` → `stages.<S>.nulls`; per Gbp *derived*)

| | S1 | P1 | P2 | P1 + P2 |
|---|---|---|---|---|
| Gbp | 1.664 | 1.017 | 9.596 | 10.613 |
| kept by the producer | 34 | 6 | 72 | 78 |
| kept per Gbp | 20.43 | 5.90 | 7.50 | 7.35 |
| N2: expected by chance (exact) | 2.433 | 3.458 | 20.915 | 24.373 |
| N2: per Gbp | 1.462 | 3.401 | 2.180 | 2.297 |
| N2: Monte Carlo per Gbp, 95% interval of the mean | 1.442 to 1.478 | 3.351 to 3.422 | 2.167 to 2.185 | 2.283 to 2.301 |
| N2: draws, 2.5% / median / 97.5% | 0 / 2 / 6 | 0 / 3 / 7 | 13 / 21 / 30 | 15 / 24 / 34 |
| N2: kept / expected | 14.0 | 1.74 | 3.44 | 3.20 |
| N2: P(K >= kept), exact | 4.8e-30 | 0.136 | 1.7e-18 | |
| N3: expected (per Gbp) | 1.006 (0.605) | 2.266 (2.228) | 12.173 (1.269) | 14.438 (1.360) |
| E13: expected (per Gbp) | 3.474 (2.087) | 4.054 (3.988) | 27.307 (2.846) | 31.361 (2.955) |
| calls at 2 in a genome with a Cas gene | 128 of 3,613 | 2,203 of 3,379 | 19,113 of 30,682 | |
| genomes with a Cas gene | 355 of 27,258 | 139 of 250 | 1,227 of 2,247 | |

Reading. Phages are gene-dense, but the chance rate per Gbp is lower on S1 than on P1 + P2 by every reading, because
only 128 of S1's 3,613 calls (3.5%) sit in one of the 355 phage genomes (1.3%) that carry a Cas gene, against 63%
of P1 and P2's calls (*derived*: (2,203 + 19,113) / (3,379 + 30,682)), and most of a phage genome is inside genes,
where rule (c) rejects a call (2,046 of S1's 3,394 extra calls lie wholly inside a gene where they are, *derived*:
3,394 - 1,348). On P1 alone the producer's 6 kept calls are not distinguishable from chance (3.46 expected, P = 0.14);
on P2 and on P1 + P2 they are 3.4 and 3.2 times chance. N2 on P1 + P2 (2.30 per Gbp) is below E-13's bound (2.96),
as it should be.

Where S1's chance mass lies (*measured*, `/tmp/if2453/e16/bin/mass.py` over `calls_S1.tsv`): only 60 of the 27,258
contigs can keep a call by chance (p_i > 0 for some call); 1.42 of the 2.43 expected under N2 (and 2.69 of the 3.47
under E13) are the 71 calls of the 33 ICP1-like genomes themselves, placed elsewhere in their own genome, and 1.01
lies in the other 27 genomes.

## The 34 kept calls (*measured*, `kept_s1.tsv`)

Group (a): two or more distinct Cas roles on two or more genes within 10 genes of the call (card E-13's "within 10
genes": at most 9 called genes lie wholly between); (b): otherwise, every Cas gene within 1 kb is of the role `cas4`
or `cas3` (a lone Cas4-like or Cas3-domain gene); (c): the rest. **34 in (a), 0 in (b), 0 in (c).**

Every one of the 34 is on a linear contig of a *Vibrio* phage (the header's host word), with Cas roles cas1, cas3,
cas5f, cas6f and cas7f on 5 genes within 10 genes, cas3 and cas5f within 1 kb (cas6f and cas7f for MW794172's
8-repeat call). The nearest Cas gene is the `Cas3-like_C_2` (PF22590, role cas3) gene, 2,793 bases, protein
`b9832d5c49f4ee66af60d0322f83e27e` in 32 rows (`429bc337a871674a27c3b317d3a84777` in MH310934), and for MW794172's
first row the `Cas_Csy4` (PF09618, cas6f) gene, protein `c46212eafda2e188710f0f2ed39f2b9d`. "Default array" is the
default producer's array on the same contig (`/tmp/if2453/e16/bin/defaults_near.py`), with its repeats and its gap
to the kept call. Coordinates 0-based half-open.

| phage (accession, header name) | contig bp | kept call | rep / sp | consensus repeat | nearest Cas gene, strand | gap bp | default array (repeats, gap) |
|---|---|---|---|---|---|---|---|
| KY883638 Vibrio phage JSF13 | 128,814 | 116740-116889 | 3 / 2 | TGTTTAAGCAGCCTATGCGGCTGCTAACT | cas3 116949-119742, -1 | 60 | 113615-113942 (6, 2,798) |
| KY883639 Vibrio phage JSF14 | 125,096 | 29394-29543 | 3 / 2 | AGTTAGCAGCCGCATAGGCTGCTTAAACA | cas3 26541-29334, +1 | 60 | 32341-32668 (6, 2,798) |
| KY883640 Vibrio phage JSF17 | 125,174 | 26474-26622 | 3 / 2 | GTTAGCAGCCGCATAGGCTGCTTAAACA | cas3 23620-26413, +1 | 61 | 29420-29988 (10, 2,798) |
| MH310934 Vibrio phage ICP1_2006_E | 128,298 | 50178-50262 | 2 / 1 | TAAGCAGCCTATGCGGCTGCTAAC | cas3 50321-53114, -1 | 59 | 47108-47377 (5, 2,801) |
| MH310935 Vibrio phage ICP1_2011_B | 125,128 | 48449-48598 | 3 / 2 | TGTTTAAGCAGCCTATGCGGCTGCTAACT | cas3 48658-51451, -1 | 60 | 45264-45651 (7, 2,798) |
| MN419153 Vibrio phage ICP1_2017_F_Mathbaria | 129,621 | 50247-50396 | 3 / 2 | TGTTTAAGCAGCCTATGCGGCTGCTAACT | cas3 50456-53249, -1 | 60 | 47062-47449 (7, 2,798) |
| MW794150 ICP1 isolate 2015_Dha_A | 125,188 | 48508-48656 | 3 / 2 | TGTTTAAGCAGCCTATGCGGCTGCTAAC | cas3 48717-51510, -1 | 61 | 45261-45710 (8, 2,798) |
| MW794151 ICP1 isolate 2016_Dha_A | 125,047 | 48467-48616 | 3 / 2 | TGTTTAAGCAGCCTATGCGGCTGCTAACT | cas3 48676-51469, -1 | 60 | 45342-45669 (6, 2,798) |
| MW794153 ICP1 isolate 2017_Dha_AA | 127,275 | 48705-48853 | 3 / 2 | TGTTTAAGCAGCCTATGCGGCTGCTAAC | cas3 48914-51707, -1 | 61 | 45263-45907 (11, 2,798) |
| MW794154 ICP1 isolate 2017_Dha_AB | 124,986 | 48329-48478 | 3 / 2 | TGTTTAAGCAGCCTATGCGGCTGCTAACT | cas3 48538-51331, -1 | 60 | 45264-45531 (5, 2,798) |
| MW794155 ICP1 isolate 2017_Dha_AC | 124,986 | 48329-48478 | 3 / 2 | as above | cas3 48538-51331, -1 | 60 | 45264-45531 (5, 2,798) |
| MW794156 ICP1 isolate 2017_Dha_AD | 124,968 | 48329-48478 | 3 / 2 | as above | cas3 48538-51331, -1 | 60 | 45264-45531 (5, 2,798) |
| MW794157 ICP1 isolate 2017_Dha_AE | 124,986 | 48406-48555 | 3 / 2 | as above | cas3 48615-51408, -1 | 60 | 45341-45608 (5, 2,798) |
| MW794158 ICP1 isolate 2017_Dha_B | 124,986 | 48406-48555 | 3 / 2 | as above | cas3 48615-51408, -1 | 60 | 45341-45608 (5, 2,798) |
| MW794159 ICP1 isolate 2017_Dha_C | 125,106 | 48526-48675 | 3 / 2 | as above | cas3 48735-51528, -1 | 60 | 45341-45728 (7, 2,798) |
| MW794160 ICP1 isolate 2017_Dha_D | 125,101 | 48521-48670 | 3 / 2 | as above | cas3 48730-51523, -1 | 60 | 45264-45723 (8, 2,798) |
| MW794161 ICP1 isolate 2017_Dha_E | 125,101 | 48521-48670 | 3 / 2 | as above | cas3 48730-51523, -1 | 60 | 45262-45723 (8, 2,798) |
| MW794162 ICP1 isolate 2017_Dha_N | 124,968 | 48329-48478 | 3 / 2 | as above | cas3 48538-51331, -1 | 60 | 45264-45531 (5, 2,798) |
| MW794163 ICP1 isolate 2017_Dha_O | 124,985 | 48329-48478 | 3 / 2 | as above | cas3 48538-51331, -1 | 60 | 45264-45531 (5, 2,798) |
| MW794164 ICP1 isolate 2017_Dha_P | 124,986 | 48329-48478 | 3 / 2 | as above | cas3 48538-51331, -1 | 60 | 45264-45531 (5, 2,798) |
| MW794166 ICP1 isolate 2017_Dha_S | 124,987 | 48329-48478 | 3 / 2 | as above | cas3 48538-51331, -1 | 60 | 45264-45531 (5, 2,798) |
| MW794170 ICP1 isolate 2017_Dha_Y | 124,984 | 48327-48476 | 3 / 2 | as above | cas3 48536-51329, -1 | 60 | 45262-45529 (5, 2,798) |
| MW794171 ICP1 isolate 2017_Dha_Z | 124,984 | 48327-48476 | 3 / 2 | as above | cas3 48536-51329, -1 | 60 | 45262-45529 (5, 2,798) |
| MW794172 ICP1 isolate 2017_Mat_I | 129,620 | 46983-47448 | 8 / 7 | TCTTTAAGCAGCCTATGCGGCTGCTAAC | cas6f 47535-48042, -1 | 87 | none on the contig |
| MW794172 ICP1 isolate 2017_Mat_I | 129,620 | 50246-50395 | 3 / 2 | TGTTTAAGCAGCCTATGCGGCTGCTAACT | cas3 50455-53248, -1 | 60 | none on the contig |
| MW794173 ICP1 isolate 2017_Mat_K | 130,453 | 51079-51228 | 3 / 2 | as above | cas3 51288-54081, -1 | 60 | 47894-48281 (7, 2,798) |
| MW794177 ICP1 isolate 2018_Mat_159 | 129,607 | 50232-50381 | 3 / 2 | as above | cas3 50441-53234, -1 | 60 | 47047-47434 (7, 2,798) |
| MW794179 ICP1 isolate 2018_Mat_164 | 128,568 | 50109-50258 | 3 / 2 | as above | cas3 50318-53111, -1 | 60 | 46984-47311 (6, 2,798) |
| MW794180 ICP1 isolate 2018_Mat_166 | 129,265 | 50109-50258 | 3 / 2 | as above | cas3 50318-53111, -1 | 60 | 46984-47311 (6, 2,798) |
| MW794182 ICP1 isolate 2018_Mat_170 | 129,607 | 50177-50326 | 3 / 2 | as above | cas3 50386-53179, -1 | 60 | 46992-47379 (7, 2,798) |
| MW794183 ICP1 isolate 2018_Mat_B | 129,629 | 50177-50326 | 3 / 2 | as above | cas3 50386-53179, -1 | 60 | 46992-47379 (7, 2,798) |
| MW794188 ICP1 isolate 2019_Mat_005 | 129,336 | 50235-50384 | 3 / 2 | as above | cas3 50444-53237, -1 | 60 | 47050-47437 (7, 2,798) |
| MW794190 ICP1 isolate 2019_Mat_C | 129,632 | 50180-50329 | 3 / 2 | as above | cas3 50389-53182, -1 | 60 | 46995-47382 (7, 2,798) |
| MZ292654 Vibrio phage ICP1 | 124,986 | 48329-48478 | 3 / 2 | as above | cas3 48538-51331, -1 | 60 | 45264-45531 (5, 2,798) |

"as above": `TGTTTAAGCAGCCTATGCGGCTGCTAACT`. The consensus is the column majority of the repeats (a tie goes to the
earlier repeat, so for two repeats that differ it is the first). The repeats of the default arrays on these contigs
are of the same family (`AATAAAGCAGCCTGTACGGCTGCTAAC` and its kin: the core `AAGCAGCC...GGCTGCTAAC` in both). In
JSF14 and JSF17 both arrays read as the reverse complement, because those two genomes are stored in the other
orientation (their Cas3 gene is on +1).

What is recognised. These are the worked example's phage I-F loci: "All 46 phage Cas7f loci (*Vibrio* phage ICP1
isolates and the JSF phages) hold Cas1 and Cas3 and no TniQ" (`examples/crispr_cas_neighbourhoods/README.md`,
section 3); all 34 kept calls are at those loci (31 calls in 30 genomes named ICP1, 3 in JSF13, JSF14, JSF17). The
layout, one array found at 3 repeats on one side of the cas genes and a second, shorter array of the same repeat
family 60 bases from *cas3* that MinCED at 3 does not report, fits the published ICP1 system with two CRISPR loci
(Seed et al. 2013, *Nature* 494:489; not re-read for this card, so no claim here rests on it). Why MinCED at a
minimum of 2 reports these 3-repeat arrays (and MW794172's 8-repeat array) that it does not report at 3 is card
E-13's open finding (its note, "Arrays at 2 and at 3"); it is not investigated here.

## Per Cas family (*measured*, `families.tsv`; per Gbp *derived*)

Expected under N2: the sum over calls of the share of kept starts within 1 kb of a gene of that family (a call beside
two families counts for both, so the column does not add up to the total); its Monte Carlo reading agrees to within
0.01 calls for every family on S1 and within 0.03 on P1 and P2 (`monte_carlo_N2` in `families.tsv`). The 21
families with the most expected chance calls; the other ten are named below the table.

| family | name | role | S1 Cas genes (genomes) | S1 expected by chance | per Gbp | S1 kept beside | P1 + P2 Cas genes | P1 + P2 expected | per Gbp | P1 + P2 kept beside |
|---|---|---|---|---|---|---|---|---|---|---|
| PF22590 | Cas3-like_C_2 | cas3 | 52 (52) | 0.970 | 0.583 | 33 | 1,277 | 8.185 | 0.771 | 10 |
| PF09614 | Cas_Csy2 | cas5f | 46 (46) | 0.804 | 0.483 | 33 | 105 | 0.346 | 0.033 | 0 |
| PF01867 | Cas_Cas1 | cas1 | 47 (47) | 0.775 | 0.466 | 0 | 1,198 | 6.475 | 0.610 | 30 |
| PF09615 | Cas_Csy3 | cas7f | 46 (46) | 0.626 | 0.376 | 1 | 111 | 0.387 | 0.037 | 9 |
| PF09618 | Cas_Csy4 | cas6f | 46 (46) | 0.595 | 0.357 | 1 | 118 | 0.346 | 0.033 | 9 |
| PF01930 | Cas_Cas4 | cas4 | 285 (285) | 0.433 | 0.260 | 0 | 635 | 4.219 | 0.398 | 7 |
| PF01881 | Cas_Cas6_C | cas6 | 2 (2) | 0.128 | 0.077 | 0 | 210 | 1.440 | 0.136 | 3 |
| PF09704 | Cas_Cas5d | cas5 | 1 (1) | 0.105 | 0.063 | 0 | 586 | 2.929 | 0.276 | 2 |
| PF09827 | CRISPR_Cas2 | cas2 | 18 (18) | 0.033 | 0.020 | 0 | 892 | 3.587 | 0.338 | 19 |
| PF21384 | Cas3_I-F_Cas2 | cas3 | 1 (1) | 0.006 | 0.003 | 0 | 90 | 0.370 | 0.035 | 1 |
| PF09611 | Cas_Csy1 | cas8f | 1 (1) | 0.004 | 0.002 | 0 | 77 | 0.273 | 0.026 | 0 |
| PF18019 | Cas3_HD | cas3 | 1 (1) | 0.000 | 0.000 | 0 | 527 | 3.869 | 0.365 | 5 |
| PF22335 | Cas10-Cmr2_palm2 | cas10 | 0 | 0 | 0 | 0 | 386 | 3.039 | 0.286 | 7 |
| PF10040 | CRISPR_Cas6 | cas6 | 0 | 0 | 0 | 0 | 200 | 2.154 | 0.203 | 9 |
| PF05107 | Cas_Cas7 | cas7 | 0 | 0 | 0 | 0 | 297 | 1.770 | 0.167 | 3 |
| PF09709 | Cas_Csd1 | cas8c | 0 | 0 | 0 | 0 | 205 | 1.398 | 0.132 | 1 |
| PF09481 | CRISPR_Cse1 | cas8e | 0 | 0 | 0 | 0 | 240 | 1.164 | 0.110 | 0 |
| PF12469 | Cmr2_N | cas10 | 0 | 0 | 0 | 0 | 98 | 1.094 | 0.103 | 3 |
| PF09485 | CRISPR_Cse2 | cas11e | 0 | 0 | 0 | 0 | 234 | 0.880 | 0.083 | 0 |
| PF18395 | Cas3_C | cas3 | 0 | 0 | 0 | 0 | 115 | 0.697 | 0.066 | 1 |
| PF09707 | Cas_Cas2CT1978 | cas2 | 0 | 0 | 0 | 0 | 200 | 0.590 | 0.056 | 10 |

The ten families with fewer Cas genes in P1 + P2 (the four Cas9 families PF22702, PF16593, PF16592, PF16595; PF09559
Cas6; PF21802 Cas3-like_C; PF21918 cas12a; PF24121, PF24122 Cas8b; PF30455 Cas3_4th) have no gene in S1 and each
attracts under 0.5 chance calls on P1 + P2. S1 holds 544 Cas genes (*measured*), 285 of them Cas4 (*measured*); S1's
genomes with a Cas3-like_C_2, Cas5f, Cas1, Cas7f or Cas6f gene are the ICP1-like phages and a few more (52, 46,
47, 46, 46 genomes).

Reading. Families with many genes in phages attract more chance calls only where those phages also have calls at 2:
Cas4, in 285 phage genomes, attracts 0.43 (0.26 per Gbp), less than each family of the ICP1 system, whose genomes
carry 71 calls between them. Nothing in S1 asks for a family to be excluded: no kept call is beside a lone Cas4 or
Cas3-domain gene, and the most a single family could add by chance on S1 is under one call (0.97, Cas3-like_C_2).

## The default producer's arrays, for scale (*measured*, `summary.json` → `stages.<S>.default_arrays`)

The same placement, the arrays of 3 repeats instead of the calls at 2, observed against expected near a Cas gene
(within 1 kb: this script's share, equal to E-13's `placement_share` to three decimals; within 10 genes: E-13's
`placement_share`):

| | S1 | P1 | P2 | P1 + P2 |
|---|---|---|---|---|
| arrays | 213 | 286 | 2,291 | 2,577 |
| within 1 kb of a Cas gene: observed / expected | 49 / 3.711 | 145 / 1.972 | 1,142 / 14.54 | 1,287 / 16.51 |
| enrichment | 13.2 | 73.5 | 78.5 | 77.9 (*derived*) |
| within 10 genes: observed / expected | 49 / 8.922 | 205 / 6.379 | 1,551 / 46.96 | 1,756 / 53.34 |
| enrichment | 5.5 | 32.1 | 33.0 | 32.9 (*derived*; E-13: 1,756 / 53.3 = 32.9) |

On phages most default arrays (164 of 213, *derived*) are far from any Cas gene, so even the default producer is less
enriched there; the filtered producer's 14-fold enrichment on S1 is the same order as the default's 13-fold within
1 kb.

## Reading, and the three outcomes

- *The rule holds on phages as it is*: yes. The chance expectation on S1 is 2.43 calls by the card's null (at most
  3.47 by E-13's bound), 7% of the 34 kept (*derived*, 2.433 / 34; 10% by the bound); the 34 are, by inspection, one
  real locus with five Cas roles in 33 genomes of one phage group; 0 are beside a lone Cas4 or Cas3-domain gene.
- *Exclude named families on virus contigs*: not supported. No family attracts as much as one chance call on S1, and
  the family the card suspected (Cas4, the commonest in phages) has none kept and 0.43 expected. So no new card and
  no new filter tag.
- *The rate cannot be bounded at this size*: no. The expected chance count is small and known to within 0.03 calls
  (Monte Carlo 95% interval 2.40 to 2.46), and its exact law gives P(K >= 34) = 4.8e-30.

The rate per Gbp (20.4) is a fact of the stage's redundancy, not of the filter: the 30 ICP1 isolates are 3.80 Mbp
of one phage (*derived*: the sum of their 30 contig lengths in the table above, 3,801,374), 0.23% of S1, and hold 31
of the 34 kept calls. A stage of other viruses (card C-24's huge phages, or metagenomic viruses) has other genomes; the
script takes `--stages` and reruns in seconds per Gbp of phages there (S1: 1.2 CPU-s for the null, *measured*).

## Open

- MinCED at 2 reports 3-repeat and 8-repeat arrays of ICP1 that MinCED at 3 does not (card E-13's finding, now seen
  in 33 phage genomes); the filtered producer recovers them, which is why its S1 calls are mostly not two-repeat
  arrays (2 repeats in 1 of 34). Not this card's question.
- On P1 alone the producer's kept calls (6) do not exceed the N2 chance expectation significantly (3.46, P = 0.14,
  *measured*); P1 + P2 together are 3.2 times chance. Noted for the PI; the PI's rule is about P1 + P2 together.

## Files

`2026-11-15_crispr_near_cas_phages/` (copied from `/tmp/if2453/e16/out` on the pod): `summary.json` (every number
of this note but the two probes'), `kept_s1.tsv` (the 34 with every column of the table above and the Cas roles
within 1 kb and 10 genes), `families.tsv` and `roles.tsv` (per stage and family or role: Cas genes, genomes, expected
under N2 and N3, Monte Carlo, kept beside), `SHA256SUMS`. Kept on the pod only: `calls_S1.tsv`, `calls_P1.tsv`,
`calls_P2.tsv` (per call: the rules where it is, p_i, E-13's f) and `draws_<S>.npz` (K per draw for N2, N3, E13),
sha256 in `SHA256SUMS`. The two probes (`defaults_near.py`, `mass.py`) are in `/tmp/if2453/e16/bin/`.
