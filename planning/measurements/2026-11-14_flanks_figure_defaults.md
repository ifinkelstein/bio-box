# `nbr flanks` figure defaults of DEC-41 (card W-32)

*Measured*, wilkcomp02, 2026-11-14 (a Saturday), 15:52 to 17:20 CST, unless labelled otherwise. Code: branch
`task/W-32` mirrored to `/tmp/if2453/ci/W-32` by `scripts/pod_pytest.sh`, frozen per driver as `/tmp/if2453/w32/tree`
(run.sh, the card's second commit) and `/tmp/if2453/w32/tree2` (run2.sh, the third); the code before the card is the
`src` of `4469fab`, copied to `/tmp/if2453/w32/base/src` before any edit of `src/`. Runtime image
`bio-box_faec1f4.sif` (sha256 `d67300aa...b9f098a`, checked by `pod_run.sh --digest`); tests and BASE_TABLES in the
test image `bio-box-test_faec1f4.sif` (sha256 `70263973...cfc2ddf4`). Release `rN.3` of `/tmp/if2453/r10/root`, bound
read-only; everything written under `/tmp/if2453/w32`. Every `bio-box` step at 2 threads, one at a time, in its own
capped unit (`pod_run.sh`, 16 GB, nice 19) under the detached drivers `w32-driver` and `w32-driver2`.

## What the card changed

- At `--level nr50` and `nr90`, when the release has a domain library, `--colour-weight anchor-nr90` and
  `--colour-fallback pfam` are the defaults of `--figure`; `locus` and `none` stay as explicit values. Otherwise, and
  at `--level pfam`, the defaults are `locus` and `none`.
- "Has a domain library": the release's manifest (`manifests/release=<R>.json`) lists a file of table `l1_domains`
  (with `--library L`, one whose `library` partition is L): `flanks.manifest_libraries`, `flanks.colour_defaults`. A
  release with only a `release_def` counts as having none, so its library is never found by listing directories.
  Without one, `--figure` at nr50 or nr90 prints one line on stderr (`note: release R has no manifest, ...` or
  `note: the manifest of release R lists no domain library (no l1_domains file): the figure keeps --colour-weight
  locus --colour-fallback none ...`). `rN.3`'s manifest lists `l1_domains` files of `pfam-38.2` for S1, P1 and P2.
- The one-colour rule (module docstring of `flanks.py`). A cluster's family is the commonest Pfam architecture of
  its genes in the run, the name the legend gives it. An architecture is a stand-in when, in some group, it reaches
  10% (weighted) and colours a gene of a cluster that is grey there. A cluster coloured in some group whose family
  is a stand-in is drawn in the architecture's Pfam colour in every group. It gets no colour or legend row of its
  own; the legend's "Pfam A" row gives A's weighted share, or the cluster's share where no gene of the group carries
  A. Two common clusters of one family stay two colours.
- `--auto-cut gap` is refused (exit 2) with a message naming `--auto-cut N`; `--auto-cut N` and the fixed cut 0.5
  (the default) stay.
- At `--level pfam` nothing changed: the defaults there were and stay `locus` and `none`. The fallback is refused
  at pfam, since card W-31, because the families are architectures already. An explicit `--colour-weight` is still
  honoured (test `test_w32_at_level_pfam_the_defaults_do_not_change`). The two Pfam-level figures of the CRISPR
  example are byte-identical to the code's before the card (below).

## The tables keep their bytes

Only `figure.svg` and `figure.pdf` change, and `README.txt` when a figure is drawn with weighted or Pfam colours:
it gains the paragraph on the options and one line naming the colours used (52 lines before, 72 now). `groups.tsv`
names a colour option only when it is given on the command line with a value other than `locus` or `none`, as
before.

On constructed releases (`tests/r2/test_w27_flanks_figure.py`, `BASE_TABLES`): six calls (nr50 with `--export gff`;
nr90 in one group with the combined distance; pfam; nr50 without a figure; the 42 mobile loci in one group; the
one-colour release), each with `--figure` but one. For each, a sha256 over the sorted `path<TAB>sha256` lines of every
file but the figures (README.txt separately) was written by `base_digests.py` with the code of `4469fab`
(`/tmp/if2453/w32/base_digests`, 2026-11-14 15:53). The card's code gives the same digests for every call, and the
same README.txt where the colours do not change (pfam, no figure).

On `rN.3`, every `nbr flanks` call of the two examples, run with the card's code and with the code before the card
on the same runs (`compare.sh`, `compare2.sh`; sha256 of every file):

| call | files compared (besides README.txt and figures) | tables | README.txt | figures |
|---|---|---|---|---|
| CRISPR `tniq_nr50` (nr50, 543 groups, no figure) | 2,719 | identical | identical | none |
| CRISPR `cas7f_all_nr50` (nr50, one group) | 9 | identical | differs (gains the colour paragraph) | 2 differ |
| CRISPR `tniq_cas_pfam` (pfam) | 8 | identical | identical | 2 identical |
| CRISPR `cas7f_all_pfam` (pfam) | 9 | identical | identical | 2 identical |
| RNA polymerase (nr50, 432 groups) | 2,164 | identical | differs (as above) | 864 differ |

With `--colour-weight locus --colour-fallback none` the card's code wrote every file of `cas7f_all_nr50` (12) and of
the RNA polymerase call (3,029, figures included) byte-identical to the code before the card's unflagged call. The
base code's figures equal the committed figures of both examples (`cas7f_all_nr50`, `tniq_cas_pfam`,
`cas7f_all_pfam`, `t7_group_c1bf639f`, SVG and PDF where checked), so the committed figures were what that code
draws. The final code's `cas7f_all_nr50` call equals run.sh's (12 files).

## A defect found on the pod, and fixed

run.sh's RNA polymerase example failed in `nbr flanks` with the card's second commit: `KeyError:
'pfam:S-AdoMet_lyase'` in `figure_plan` (`/tmp/if2453/w32/rnap/logs/flanks.log`). A cluster absorbed by the
one-colour rule was coloured in a group where none of its genes carries its family's architecture, so the legend had
no share for the row. The third commit gives the row the cluster's own share there. The test
`test_w32_the_rule_holds_where_the_clusters_genes_carry_no_domain` failed first without the guard (`KeyError:
'pfam:Cas1'`), on a constructed release whose common Cas1 cluster has no domain in one group. run2.sh then ran the
RNA polymerase example again, with that code, as a new run (`c93ee2943d3f18b595be3a03841798ef`) and compared on it.

## The examples' checks with the card's code

- CRISPR (`run.sh`, OUT `/tmp/if2453/w32/crispr`, runs TniQ `91479292...`, Cas7f `4763088c...`): `check: 50 ok, 14
  n/a, 0 FAIL, 0 differs, 0 new, 0 missing, 0 unexpected; invariants hold`. The 14 n/a are the CRISPR array checks
  (`rN.3` has no arrays). `recorded/rN.3.tsv` is unchanged and matched.
- RNA polymerase (`run.sh`, OUT `/tmp/if2453/w32/rnap2`): `check: 42 ok, 0 FAIL, 0 differs, 0 new, 0 missing, 0
  unexpected; invariants hold`.

## What the regenerated figures draw

From the plans `render_figure` receives (`figure_stats.py`, run2.sh), defaults against the old colours, on the drawn
rows only (60 per figure):

| figure | setting | flanking genes drawn | coloured by a cluster | coloured by a Pfam colour | grey | colours |
|---|---|---|---|---|---|---|
| `cas7f_all_nr50` | before (locus, none) | 720 | 50 | 0 | 670 | 13 |
| | DEC-41 (anchor-nr90, pfam) | 720 | 0 | 251 | 469 | 8 |
| `t7_group_c1bf639f` | before | 545 | 419 | 0 | 126 | 13 |
| | DEC-41 | 545 | 127 | 338 | 80 | 14 |

- `cas7f_all_nr50`: the legend was 13 NR50 clusters, at 27% to 40% of the 157 loci. One was Cas1 (`d10f704c`); the
  others were mostly clusters of the 42 near-identical ICP1 loci, several without a Pfam name. Now it is 8 Pfam
  families with weighted shares: `Cas_Csy4` 99%, `Cas_Csy2` 97%, `Cas_Csy1` 73%, `Cas_Cas1` 71%, three Cas2-Cas3
  fusions (26%, 18%, 14%) and `TniQ` 14%. The Cas genes of the Tn7-like loci are in colour. Card W-31 measured 37 of
  their 39 coloured with these two settings. The one-colour rule changes only which colour a coloured gene takes,
  never whether it is grey (*derived* from the rule), so the count stands.
- `t7_group_c1bf639f`: in one group per anchor NR90 cluster the weighting changes nothing. The fallback with the
  rule draws five conserved neighbours as their Pfam architecture instead of their cluster, because rarer clusters
  of the same architecture take its colour in other groups: the ligase 96%, `T7-like_gp12` 95%, `Phage_T7_1_1` 95%,
  `BC10` 89% and the kinase 88%. It adds `Phage_T7_GP04`, `DUF3310` and `HNH_3`. 46 more genes are coloured (465
  against 419).

Reading: on the CRISPR figure the defaults do what DEC-41 wanted. On grouped runs such as the RNA polymerase one, the
rule's "throughout" replaces cluster ids by architecture names in groups where the cluster alone would have been
unambiguous. That is the rule as decided; the cluster ids stay in `conservation.tsv`.

## Cost

*Measured*, units' `CPUUsageNSec`. On one group, unchanged: the CRISPR `cas7f_all_nr50` call took 49.5 CPU-s with the
card's code and 47.5 with the code before it; `tniq_nr50` (no figure) 50.8 and 53.9; `tniq_cas_pfam` 30.3 and 30.7.
On 432 groups, dearer: the RNA polymerase call with its 432 figures took 305.1 CPU-s with the defaults, against 259.1
before the card and 259.4 with the old colours (+18%, *derived*: 305.1 / 259.1). The extra work is the weighted and
per-architecture shares of every group (`run_palette`, `arch_shares`). Peak unit memory 2.8 to 4.4 GB in every
call. Both drivers took 73 minutes of wall time together (16:05:53 to 17:19:15), with no wait for room.

## Files here

`2026-11-14_flanks_figure_defaults/`: `base_digests.py` (BASE_TABLES with the code before the card), `run.sh` and
`run2.sh` (the drivers), `compare.sh` and `compare2.sh` (the comparisons), `figure_stats.py` (the table above). The
outputs stay on the pod under `/tmp/if2453/w32` (`compare.txt`, `sums/`, `stats_*.json`).
