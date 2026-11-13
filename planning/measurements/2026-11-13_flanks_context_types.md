# `nbr flanks` on mobile elements: context types and figure colours (card W-31)

*Measured*, wilkcomp02, 2026-11-13, 12:41 to 13:05 CST, tree `task/W-31` mirrored to `/tmp/if2453/w31/tree`, image
`bio-box-test_64f1bf0.sif` (sha256 `3d77df6b...fbda4`), by `2026-11-13_flanks_context_types/run.sh` (units
`w31-driver` to `w31-driver4`; each `bio-box` command through `scripts/pod_run.sh` in its own unit `w31-<name>` at 2
threads, 16 GB, nice 19; each Python script in a `systemd-run --scope` of 2 CPUs and 16 GB, as the CRISPR example's
`run.sh` does). Read only, bound read-only into the image: the root `/tmp/if2453/r10/root` (release `rN.3`), the
CRISPR example's outputs `/tmp/if2453/ex2/rN.3-b` (runs TniQ `ef8d86b8...`, Cas7f `47db5988...`) and the RNA polymerase
run `599bea5d...` of `2026-11-12_flanks_rN3.md` (`/tmp/if2453/w25/runs`, still there). Outputs in
`/tmp/if2453/w31/out`. The account was at 34 to 43 CPUs of other work (a Pfam scan and other agents' tests); every
step waited for room under the 40-thread day cap.

Every number below is *measured* by that run unless labelled otherwise. The tables here are copied from
`measure.py`'s output (`types.tsv`, `colours.tsv`, `heights.tsv`, `summaries.json` in the directory beside this note).

## What was compared

Six `nbr flanks` views, with the settings the example and the flanks note used:

| view | anchors | settings | loci | groups | classes |
|---|---|---|---|---|---|
| cas7f_all_nr50 | Cas7f run | `--window 6g --level nr50 --group-by none` (example, figure 3) | 157 | 1 | 19 cast, 46 phage, 92 prok |
| cas7f_all_pfam | Cas7f run | `--window 6g --level pfam --group-by none` (example, extra) | 157 | 1 | same |
| tniq_cas_pfam | the 28 TniQ loci with a Cas domain | `--window 10g --level pfam --group-by none` (example, figure 1) | 28 | 1 | 27 cast, 1 prok |
| tniq_nr50 | TniQ run | `--window 10g --level nr50` (example, step 5a) | 635 | 543 | 33 cast, 602 prok |
| tniq_all_nr50 | TniQ run | as tniq_nr50 with `--group-by none` | 635 | 1 | same |
| rnap_nr50 | RNA polymerase run | `--window 5g --level nr50` (flanks note) | 2,581 | 432 | 2,555 phage, 26 prok |

Classes (`measure.py`): *phage*, genome kind virus; *cast*, a Tn7-like CRISPR-associated locus: the anchor or a
flanking gene has a TniQ domain and a flanking gene a domain whose Pfam name starts with `Cas`; *prok*, the rest. This
is broader than the example's `q03` (27 TniQ loci with a Cas domain within 10 genes against 33 here, since any `Cas*`
family counts), so "cast" is a superset of the example's I-F3 class.

`measure.py` runs the unflagged `nbr flanks` through `biobox.cli.main`, captures the loci, and recomputes in-process,
with the functions `nbr flanks` itself uses, the context types under each cut rule and the colour of every flanking
gene under each colouring; its fixed-cut types equal the CLI's `groups.tsv` in all six views (asserted). The flagged
CLI runs of step 4 gave the same type counts (for example `--auto-cut gap` on cas7f_all_nr50: 115 types at cut 0.1).

## The defaults are untouched

- Fixture: `fixture_identity.py` ran seven unflagged `nbr flanks` calls (default; `--figure --export gff`; `--level
  pfam --figure`; `--group-by none --distance combined --figure`; `--level nr90 --distance a --cut 0.3`; and two on the
  card's new constructed release with `--figure`) once with the code before the card (`origin/plan/task-cards`,
  `d9d4970`, src only) and once with the card's code: `diff -r` of the two output trees (123 files, figures included)
  was empty.
- `rN.3`: the example's three `nbr flanks` calls, rerun unflagged with the card's code, are byte-identical to the
  example's outputs written by the code before the card: `tniq_nr50` 2,720 files, `cas7f_all_nr50` 12 files (figure
  included), `tniq_cas_pfam` 11 files (`diff -r -q`, empty).
- In the tests, the flags' default values (`--colour-weight locus --colour-fallback none`) give the unflagged bytes,
  and `--auto-cut` leaves the matrices and the tree byte-identical to the unflagged call's.

## (a) Context types: a cut chosen from the tree

Rules: fixed cut 0.5 (today), 0.7 and 0.8 for reference; `--auto-cut gap` (the middle of the largest gap between the
distinct merge distances below 1 and 1 itself); `--auto-cut N` (the least cut leaving at most N types, never at
distance 1). Columns: types; loci in types of one locus; the five largest types (loci); cast loci: types they fall in
(largest type); phage loci: types (largest); types mixing a class with another.

| view | rule | cut | types | loci alone | largest | cast | phage | mixed |
|---|---|---|---|---|---|---|---|---|
| cas7f_all_nr50 | fixed 0.5 | 0.5 | 108 | 102 | 45,2,2,2,2 | 18 (2) | 2 (45) | 0 |
| | fixed 0.8 | 0.8 | 66 | 42 | 45,7,7,6,4 | 14 (3) | 2 (45) | 0 |
| | gap | 0.1 | 115 | 114 | 43,1,1,1,1 | 19 (1) | 4 (43) | 0 |
| | N=10 | 0.9986 | 18 | 12 | 72,46,13,10,2 | 6 (13) | 1 (46) | 0 |
| | N=20 | 0.9738 | 20 | 12 | 72,46,10,7,3 | 8 (7) | 1 (46) | 0 |
| cas7f_all_pfam | fixed 0.5 | 0.5 | 152 | 147 | 2,2,2,2,2 | 18 (2) | 46 (1) | 0 |
| | fixed 0.8 | 0.8 | 32 | 18 | 51,45,8,7,6 | 13 (3) | 2 (45) | 1 |
| | gap | 0.45 | 156 | 155 | 2,1,1,1,1 | 19 (1) | 46 (1) | 0 |
| | N=10 | 0.8767 | 10 | 7 | 121,18,11,1,1 | 2 (18) | 1 (46) | 1 |
| tniq_cas_pfam | fixed 0.5 | 0.5 | 28 | 28 | 1,1,1,1,1 | 27 (1) | - | 0 |
| | fixed 0.8 | 0.8 | 24 | 20 | 2,2,2,2,1 | 23 (2) | - | 0 |
| | gap | 0.7298 | 27 | 26 | 2,1,1,1,1 | 26 (2) | - | 0 |
| | N=10 | 0.9185 | 10 | 4 | 7,7,4,2,2 | 10 (7) | - | 1 |
| tniq_all_nr50 | fixed 0.5 | 0.5 | 567 | 502 | 3,3,3,2,2 | 33 (1) | - | 1 |
| | fixed 0.8 | 0.8 | 487 | 380 | 7,6,6,6,6 | 27 (2) | - | 2 |
| | gap | 0.1405 | 630 | 625 | 2,2,2,2,2 | 33 (1) | - | 0 |
| | N=10 | 1.0 (below 1) | 124 | 80 | 409,22,16,8,7 | 8 (14) | - | 7 |
| tniq_nr50 (543 groups) | fixed 0.5 | 0.5 | 614 | 596 | 3,3,3,2,2 | 33 (1) | - | 1 |
| | gap | 0 to 0.987, median 0.872 | 558 | 503 | 5,4,4,4,3 | 29 (3) | - | 3 |
| | N=10 | 0 (median) | 633 | 631 | 2,2,1,1,1 | 33 (1) | - | 0 |
| rnap_nr50 (432 groups) | fixed 0.5 | 0.5 | 598 | 346 | 227,141,132,121,65 | - | 572 (227) | 0 |
| | fixed 0.7 | 0.7 | 463 | 266 | 306,214,170,152,121 | - | 437 (306) | 0 |
| | fixed 0.8 | 0.8 | 441 | 251 | 306,214,175,152,122 | - | 415 (306) | 0 |
| | gap | 0.286 to 0.921, median 0.694 | 447 | 256 | 306,214,175,152,121 | - | 422 (306) | 0 |
| | N=10 | 0 to 0.603, median 0 | 874 | 632 | 227,132,125,108,74 | - | 848 (227) | 0 |

(The fixed 0.7 rows and N=20 of the other views are in `types.tsv`.) "1.0 (below 1)": the highest merge below 1, which
rounds to 1.0 at 4 decimals. tniq_all_nr50 has about 125 merges at distance 1 (`heights.tsv`), which the rule never
makes, so N=10 and N=20 both stop at the same 124 types.

Why the largest gap fails on these trees (`heights.py` on the four one-group trees, `heights.tsv`):

| tree | merges | at 0 | below 0.5 | 0.5 to 1 | at 1 | three largest gaps [from, to] |
|---|---|---|---|---|---|---|
| cas7f_all_nr50 | 115 | 1 | 7 | 91 | 17 | 0.200 [0, 0.200], 0.100 [0.400, 0.500], 0.086 [0.200, 0.286] |
| cas7f_all_pfam | 156 | 0 | 1 | 155 | 0 | 0.100 [0.400, 0.500], 0.088 [0.500, 0.588], 0.051 [0.615, 0.667] |
| tniq_cas_pfam | 27 | 0 | 0 | 27 | 0 | 0.040 [0.710, 0.750], 0.038 [0.750, 0.788], 0.033 [0.855, 0.888] |
| tniq_all_nr50 | 632 | 0 | 66 | 441 | 125 | 0.059 [0.111, 0.170], 0.045 [0.050, 0.095], 0.039 [0.261, 0.300] |

The merge distances of these mobile-element trees are spread evenly from 0.5 to 1 with no gap wider than 0.04 to 0.1
there; the largest gap lies among the few low merges (near-identical loci), so the gap rule cuts low and makes more
types than the fixed cut. The trees have no gap to find. On the phage groups of the RNA polymerase run, which do have
tight clades, the gap rule gives 447 types (256 loci alone) against 598 (346), about what a fixed 0.8 gives (441, 251),
and still never mixes phage and prokaryote loci.

A target number of types makes the one-group views readable: on cas7f_all_nr50 `--auto-cut 10` gives three large types,
the complete type I-F systems (72 loci), the ICP1 phages (all 46 in one type) and 13 of the 19 Tn7-like I-F3 loci,
with no type mixing the three classes (`cas7f_all_nr50_n10_fb_pfam_w_nr90.svg`, types 1, 2 and 3). It does so by
cutting at 0.97 to 0.999, where a type means loci that share something on average, not a conserved context; on the
TniQ loci in one group it lumps 409 of 635 loci into one type; and in grouped views (hundreds of groups of a few loci)
N is larger than most groups, so the rule leaves every leaf apart (633 types against 614).

## (b) Colours: weighting loci before the 10% threshold, and (c) the Pfam fallback

Columns: flanking genes coloured of all (by the Pfam fallback); the Cas genes of cast loci coloured; the genes of phage
loci coloured; Cas genes of prok loci coloured; colours in the run's legend palette.

| view | weight | fallback | coloured (by Pfam) | cast Cas genes | phage genes | prok Cas genes | colours |
|---|---|---|---|---|---|---|---|
| cas7f_all_nr50 (1,878 genes) | locus (today) | none | 596 (0) | 0 of 39 | 532 of 546 | 61 of 430 | 13 |
| | leaf | none | 106 (0) | 0 of 39 | 0 of 546 | 103 of 430 | 4 |
| | genome | none | 611 (0) | 0 of 39 | 532 of 546 | 76 of 430 | 14 |
| | anchor-nr90 | none | 106 (0) | 0 of 39 | 0 of 546 | 103 of 430 | 4 |
| | locus | pfam | 993 (397) | 37 of 39 | 535 of 546 | 397 of 430 | 22 |
| | anchor-nr90 | pfam | 599 (493) | 37 of 39 | 139 of 546 | 399 of 430 | 12 |
| cas7f_all_pfam (1,878) | locus | (level pfam) | 682 | 37 of 39 | 227 of 546 | 397 of 430 | 10 |
| | anchor-nr90 | (level pfam) | 591 | 37 of 39 | 139 of 546 | 394 of 430 | 8 |
| tniq_cas_pfam (560) | locus | (level pfam) | 157 | 72 of 85 | - | - | 24 |
| | genome | (level pfam) | 139 | 68 of 85 | - | - | 21 |
| tniq_all_nr50 (12,553) | any | none | 0 | 0 of 91 | - | - | 0 |
| | any | pfam | 1,090 (1,090) | 0 of 91 | - | - | 6 |
| tniq_nr50, 543 groups (12,553) | any | any | 12,553 | 91 of 91 | - | - | 8,516 |
| rnap_nr50, 432 groups (24,473) | locus | none | 23,266 | - | 23,006 of 24,213 | - | 3,157 |
| | leaf | none | 23,383 | - | 23,123 of 24,213 | - | 3,194 |
| | locus | pfam | 23,637 (371) | - | 23,377 of 24,213 | - | 3,171 |

Readings:

- Today, at NR50 with every Cas7f locus in one group, the 42 near-identical ICP1 loci (one leaf, one anchor NR90
  cluster) hold 532 of the 596 coloured genes, and no Cas gene of a Tn7-like locus is coloured (0 of 39). One
  prokaryote Cas1 cluster (`d10f704c`, in 40% of the loci) is coloured, so "every prokaryote Cas gene is grey" in the
  example's README is nearly but not exactly right: 61 of 430 are coloured.
- Weighting alone (leaf or anchor NR90) removes the phage loci's monopoly (0 phage genes coloured) but colours no
  more Cas genes of the Tn7-like loci (0 of 39): each species' Cas7f, Cas5f and Cas6f clusters are rare however the
  loci are weighted. Genome weighting changes almost nothing here (611 against 596): the 42 ICP1 loci are on 42
  genomes, and it collapses only loci that share a genome. Under the default grouping (one group per anchor NR90
  cluster) anchor-nr90 weighting is a no-op by construction, and leaf weighting moves the RNA polymerase run by 0.5%.
- The Pfam fallback is what colours the Tn7-like loci's Cas genes: 37 of 39 (the other two carry Cas families held by
  fewer than 10% of the loci), and 397 of 430 prokaryote Cas genes, as the Pfam level does (37 of 39). With anchor-nr90
  weighting as well, the ICP1 clusters lose their own colours (139 of the 546 phage genes stay coloured, all by their
  Pfam family) and the figure needs 12 colours instead of 22
  (`cas7f_all_nr50_fb_pfam_w_nr90.svg`: Cas6f `Cas_Csy4`, Cas5f `Cas_Csy2`, Cas8f `Cas_Csy1`, Cas1, the Cas2-Cas3
  fusions and TniQ each in their colour, in 99%, 97%, 73%, 71%, 26% to 14% and 14% of the loci, weighted).
- The fallback has a cost: a family can carry two colours, its common cluster's own and the Pfam colour of its rare
  clusters (in that figure, Cas1 as cluster `d10f704c` and as "Pfam Cas_Cas1"; Cas5f as two clusters and as "Pfam
  Cas_Csy2"). It cannot help where the family itself is rare: across all 635 TniQ loci in one group no NR50 cluster
  reaches 10%, and the Cas families, in 33 of 635 loci (5%), stay grey with the fallback too.
- Cost: the flags add no measurable time. The flagged Cas7f figure runs took 43.9 to 46.0 CPU-s and at most 4.08 GB
  (unit `MemoryPeak`, which counts page cache) against 46.1 CPU-s and 4.10 GB unflagged (`pod_run` usage lines).

## Proposed defaults (the PI's choice; nothing was changed)

1. **Keep `--cut 0.5`** as the default. No cut chosen from the tree is better on the mobile-element trees, which have
   no gap above 0.5; the largest-gap rule made more types than 0.5 on three of the four one-group views (115 against
   108, 156 against 152, 630 against 567) and one fewer on the fourth (27 against 28, for 28 loci), and on phage groups
   it is close to a fixed 0.8. Recommend
   **dropping `--auto-cut gap`** before it is documented (or keeping it undocumented), and keeping **`--auto-cut N`**
   as an option for one-group views (`--group-by none`) of a few hundred loci, documented as coarse (on
   cas7f_all_nr50, N=10 separates the complete I-F systems, the phages and the Tn7-like loci without mixing them, at a
   cut of 0.9986). For mobile elements the guide should say to read the tree and the figure, not the types.
2. **`--colour-fallback pfam` on by default** at `--level nr50` and `nr90` when the release has a domain library: it is
   the change that colours the Cas genes of the Tn7-like loci (37 of 39 against 0) and costs nothing measurable; the
   price is that one family may show in two colours, which the legend names.
3. **`--colour-weight anchor-nr90` on by default** (leaf gives the same numbers here): near-identical loci count once,
   so 42 copies of one phage no longer take the colours (0 against 532 phage genes coloured without the fallback);
   under the default grouping it changes nothing. Not `genome`: it collapses only loci that share a genome and changed
   15 genes of 1,878 here.

With 2 and 3 the figure of the example's step 5b would be `cas7f_all_nr50_fb_pfam_w_nr90.svg`. Taking them would
change the bytes of every `--figure` output at nr50 or nr90 (the tables are untouched); `GUIDE.md` §6, `AGENT_GUIDE.md`
and the example's README would then need the new defaults (not done: those files are not this card's).

## Files here

`2026-11-13_flanks_context_types/`: `run.sh` (the driver), `measure.py` (the comparison), `fixture_identity.py` (the
unflagged identity check), `heights.py` (merge distances of a tree), `types.tsv` and `colours.tsv` (every rule and
colouring of every view, more columns than above), `heights.tsv`, `summaries.json` (each view's arguments and class
counts), and two figures of the Cas7f loci at NR50 in one group (`--window 6g`): `cas7f_all_nr50_fb_pfam_w_nr90.svg`
(`--colour-fallback pfam --colour-weight anchor-nr90`, cut 0.5) and `cas7f_all_nr50_n10_fb_pfam_w_nr90.svg` (the same
with `--auto-cut 10`). Today's figure of the same loci is `examples/crispr_cas_neighbourhoods/figures/cas7f_all_nr50.svg`.
Both were converted to PNG and checked by eye.
