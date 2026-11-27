# The release-diversity figure: an audit of its numbers, and the planned releases (card A-31)

- **Date:** 2026-11-27. Branch `task/A-31`. Host for every *measured* number of this card: `wilkcomp02`.
- **Files:** `scripts/plot_release_diversity.py` (collector, checks, renderer), `planning/measurements/release_diversity.json`
  (the data; schema 3), `tests/r1/test_a31_release_diversity.py`, the figure
  `planning/measurements/2026-11-27_release_diversity/release_diversity.png` and `.svg`.
- **Labels (AGENTS §6):** *measured* (host, date, command), *derived* (formula and inputs), *estimate*. Outside
  benchmarks are labelled *measured* with the publisher's citation and `re_verified: false`: nothing was fetched
  (AGENTS §8) and none was checked again.
- **The old figure and its note are unchanged** (`2026-11-24_release_r26.11.3/source_diversity_r26.11.3.png`, section 9 of
  `2026-11-24_release_r26.11.3.md`): they are the record of `r26.11.3`.

## 0. Regenerate

On the pod, in the image, one thread (the tree mirrored to `$TREE` with `rsync`, as `scripts/pod_pytest.sh` does):

```
singularity exec --cleanenv --env OMP_NUM_THREADS=1 --bind /tmp/if2453 /stor/home/if2453/bio-box/images/bio-box-test_faec1f4.sif python $TREE/scripts/plot_release_diversity.py render --out $TREE/planning/measurements/2026-11-27_release_diversity/release_diversity
```

`check` prints the table the figure draws and every problem of the data file (exit 1 on any); `fill` recomputes every
formula after an input changes; `tower` prints the `tower` remainder; `root` re-reads the data root (section 1).

## 1. The audit's read of the root

*Measured*, `wilkcomp02`, 2026-11-27 07:45 CST, unit `a31-root` (`systemd-run --user --collect -p MemoryMax=8G
-p MemorySwapMax=0`, `nice -n 19`, test image `bio-box-test_faec1f4.sif`, the root bound read-only, DuckDB at 2 threads
and 4 GB): `plot_release_diversity.py root --root /tmp/if2453/r10/root --built V2,PL2`. It ran 16 s (*measured*, the
unit's own timings: 8.4 s and 7.7 s for the two cluster reads, under 0.1 s per stage). It reads each release's
manifest (`manifests/release=<id>.json`), the `contig_index` files the manifest names for each stage (default-variant
rows: count and sum of `contig_len`), the manifest's `cluster_member` row counts, and, where the release's
`cluster_member` files are still in the root (`r26.11.2`, `r26.11.3`), the rows per level and the rows whose member is
its cluster's representative (one per cluster). `rN.3`, `rM.2` and `r26.11.1` are archived (`archived/release=<id>.json`):
their cluster files are gone, their manifests and stage tables stay. `V2` and `PL2` are built and in no manifest yet: their
`contig_index` partition directory was listed (the one place this card reads a directory, read-only).

| Stage | Default-variant contigs | Rows | Default-variant bases | Earlier record of the same number |
|---|---:|---:|---:|---|
| `S1` | 27,258 | 27,258 | 1,664,109,938 | `2026-11-12_ncrna_set.md` |
| `P1` | 539 | 539 | 1,016,625,803 | `2026-11-14_crispr_min_repeats.md` |
| `P2` | 5,018 | 5,018 | 9,596,180,328 | the same |
| `P3` | 9,107 | 9,107 | 9,949,077,522 | none found (sum P1..P6 equals `rM.2`'s 52,132,698,276 of `2026-11-14_pack_build_in_parts.md`) |
| `P4` | 4,005 | 4,005 | 9,976,456,950 | as `P3` |
| `P5` | 4,111 | 4,111 | 9,967,711,632 | as `P3` |
| `P6` | 4,227 | 4,227 | 9,962,536,103 | as `P3` |
| `H1` | 6,714 | 6,714 | 1,819,831,338 | `2026-11-15_stage_H1.md` |
| `PL1` | 61,647 | 61,647 | 5,924,552,875 | `CATALOG_GOAL.md` §2 |
| `C5` | 2,653 | 2,653 | 9,109,234,757 | `2026-11-21_release_r26.11.2.md` |
| `V1` | 168,358 | 168,358 | 5,178,493,077 | the same |
| `G1` | 1,227,214 | 1,227,214 | 10,042,502,859 | `2026-11-24_release_r26.11.3.md` |
| `SPX1` | 2,817,316 | 2,817,316 | 35,688,628,875 | the same |
| `V2` (built 2026-11-26) | 243,466 | 243,915 | 8,199,089,073 | card A-26's new bases; card R-17's note |
| `PL2` (built 2026-11-26) | 112,562 | 112,562 | 1,974,983,002 | card A-26's new bases |

| Release | Stages | Stored bases (sum of its stages) | Contigs | `cluster_member` rows | Distinct proteins | NR90 clusters | NR50 clusters |
|---|---|---:|---:|---:|---:|---:|---:|
| `rN.3` | `S1 P1 P2` | 12,276,916,069 | 32,815 | 21,241,368 | 10,620,684 (rows / 2) | files gone | files gone |
| `rM.2` | `S1 P1..P6` | 52,132,698,276 | 54,265 | 91,631,686 | 45,815,843 (rows / 2) | files gone | files gone |
| `r26.11.1` | `+ H1 PL1` | 59,877,082,489 | 122,626 | 98,204,082 | 49,102,041 (rows / 2) | files gone | files gone |
| `r26.11.2` | `+ C5 V1` | 74,164,810,323 | 293,637 | 115,560,352 | 57,780,176 | 35,475,280 | 10,091,519 |
| `r26.11.3` | `+ G1 SPX1` | 119,895,942,057 | 4,338,167 | 204,675,750 | 102,337,875 | 71,208,987 | 22,155,373 |

Every number in these two tables equals the release note or `c35_check.py` output it was compared with: the stored bases
of all five releases (`rN.3` and `rM.2` as the sums of their stages in the notes above), the proteins of all five, and the
NR90 and NR50 counts of `r26.11.2` and `r26.11.3` (`2026-11-21_release_r26.11.2/summary.json`, `2026-11-24_release_r26.11.3.md`
§5.1). For the three archived releases the cluster counts are those `c35_check.py` and card B-42 recorded when the files
existed: `rN.3` 2,966,200 NR50 and 8,278,537 NR90 (`2026-11-11_build_r50.md`, A2); `rM.2` 8,192,065 and 31,183,876, and
`r26.11.1` 9,291,312 and 33,117,773 (`2026-11-16_release_r26.11.1/summary.json`). Each stage's index rows equal its
default-variant contigs except `V2`, whose 449 other rows are non-default variants of contigs other stages hold (card
R-17); the stored bases count default variants only, as every release note does.

## 2. Differences the audit found (old value, new value)

The old data file was `release_diversity.json` schema 2 (commit `bbc6de1`); "old figure" is
`source_diversity_r26.11.3.png`.

| # | What | Old | New | Why |
|---|---|---|---|---|
| 1 | `rM.2` NR90 clusters (data only; the old figure does not draw NR90) | 30.89 M, no label, source "`2026-11-16_release_r26.11.1.md`; README.md" | **31,183,876**, *measured* (`c35_check.py`, `2026-11-16_release_r26.11.1/summary.json`, `clusters.rM.2.nr90`) | the cited note says 31,183,876 (its §6.2 table); 30.89 M is in no file of the tree and is close to the 30,867,052 clusters that kept their members |
| 2 | `rN.3` NR90 clusters | null | **8,278,537**, *measured* (card B-42's comparison, `2026-11-11_build_r50.md`) | recorded, left out |
| 3 | Date labels of `rM.2`, `r26.11.1`, `r26.11.2`, `r26.11.3` | Nov 15, Nov 19, Nov 22, Nov 25 (no rule stated, no source) | **Nov 14, Nov 18, Nov 21, Nov 24** | one rule: the day the release's manifest was written (*measured*, `ls -l /tmp/if2453/r10/root/manifests`, 2026-11-27). The old labels were closeout days (worked examples, announcement) without a source. `rN.3`'s Nov 11 is the same by both rules |
| 4 | `r26.11.4` column, states | the whole bar hatched "queued", 287.0 Gbp, including the 119.9 Gbp already released | **119.896 released, 10.174 built (`V2`, `PL2`, *measured* on the root), 156.964 staged**; total 287.034 Gbp (*derived*, sum) | the old bar drew measured sequence as queued. The class sums of the column are unchanged (isolate 59.578, uncultured 60.891, virus 83.200, plasmid 24.417, metagenome 56.295, organelle 2.653 Gbp) and equal card A-26's new bases |
| 5 | `r26.11.4` NR50 families | about 64 M, open square, "estimate (DEC-77)" | **not drawn** | DEC-77 does not state 64 M. The number is the hybrid design's table (`planning/design/hybrid_linclust_nr90_cascade_nr50.md` §2: cascade NR50 representatives "~24%" of 265 M proteins, an *estimate*); the measured ratio of `r26.11.3` is 21.6% (*derived*: 22,155,373 / 102,337,875), and it rose from 17.5% in `r26.11.2`, so no stated formula with measured inputs gives the NR50 count of an unbuilt release |
| 6 | `r26.11.4` distinct proteins | about 265 M, one open circle, "(DEC-77)" | **265.1 M with its range 252.8 to 286.2 M**, *estimate* | DEC-77 states "~265M"; the formula is the hybrid design's §2: 102,337,875 + (0.90, 0.974, 1.10) M per new Gbp x 167.138 Gbp (*derived*; 0.974 measured on `G1` and `SPX1`). The old point hid the range the design states |
| 7 | Old footnote: "rN.3 to r26.11.3 measured on wilkcomp02 (c35_check.py / release notes)" | | `rN.3` and `rM.2` predate `c35_check.py` (card C-35); their numbers are the B-42 comparison's and `c35_check.py`'s on `r26.11.1` | wording |
| 8 | Class of `GART1`, `GART2` | "Uncultured archaea & MAGs" | the same stages, the class renamed "Uncultured prokaryotes (MAGs, single cells, symbionts)" | `GART1` and `GART2` are bacteria; no stage changes class, so no number changes |
| 9 | Legend stage lists | the virus class named `V2`, `V3`, `CG1`, `MV1`, `SAT1` (then queued) and not `LRV1`, `MGE2`, which it held | rules in the data file (section 4), not lists | the legend did not say which stages it covered |
| 10 | The index `D` | "D = exp(H) effective classes (single-domain DB = 1.0, max = 6.0)" | the same formula and values; caption "D = exp(H) over 6 classes: 1.0 for one, 6.0 for 6 equal" | the formula is right (Hill number of order 1 over the six classes' shares of stored Gbp in each column, every state counted); "single-domain DB" was wrong (one class, not one domain). Values: `rN.3` 1.49, `rM.2` 1.15, `r26.11.1` 1.71, `r26.11.2` 1.88, `r26.11.3` 3.51, `r26.11.4` 4.89 (*derived*, `check`), equal to the old figure's to its one decimal |
| 11 | Panel d, new NR50 per stored Gbp | H1 489.4k, SPX1 267.4k, G1 266.0k, V1 135.0k, PL1 33.9k, C5 10.5k | the same: 489,360; 267,394; 266,041; 135,037; 33,876; 10,463 (*derived*: `c35_check.py`'s new NR50 clusters with a member of the stage, over its stored Gbp of section 1) | values agree. Not said before: each stage is measured against the release before its own (`H1`, `PL1` against `rM.2`; `C5`, `V1` against `r26.11.1`; `G1`, `SPX1` against `r26.11.2`), so the yields are against different bases; the panel now says so |
| 12 | `RUNBOOK.md`'s command | `uv run --project /Users/ilya/projects/dataviz-python ...` on the laptop | on the pod, in the image (section 0) | AGENTS §5.5: no compute on the laptop. The script no longer reads the laptop's `dataviz-python`; its theme is the same Wilke style written out |
| 13 | `r26.11.3` callout "+120% NR50" | | the same, now *derived* in the data: (22,155,373 / 10,091,519 - 1) x 100 = 119.5, drawn as 120 | none |

Everything else in the old file agreed with its sources: every stored Gbp by class of the five released releases (to the
three decimals it gave), every protein and NR50 count, the counterfactual's rate (10,463 per Gbp) and the `r26.11.4`
class sums.

## 3. Totals by release and state

Gbp (*derived* by `plot_release_diversity.py check`, sums of the stages' bases; cumulative, each column holds every
earlier release). States: released (in a released release, *measured*), built (*measured* on the root, in no release),
staged (*measured* new bases of the staging notes, not built), planned (*estimate*, drawn at the middle of its range).

| Column | Released | Built | Staged | Planned | Total | Range of the total | Isolate | Uncultured | Virus | Plasmid | Metagenome | Organelle | D |
|---|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|
| `rN.3` | 12.277 | 0 | 0 | 0 | 12.277 | | 10.613 | 0 | 1.664 | 0 | 0 | 0 | 1.49 |
| `rM.2` | 52.133 | 0 | 0 | 0 | 52.133 | | 50.469 | 0 | 1.664 | 0 | 0 | 0 | 1.15 |
| `r26.11.1` | 59.877 | 0 | 0 | 0 | 59.877 | | 50.469 | 0 | 3.484 | 5.925 | 0 | 0 | 1.71 |
| `r26.11.2` | 74.165 | 0 | 0 | 0 | 74.165 | | 59.578 | 0 | 8.662 | 5.925 | 0 | 0 | 1.88 |
| `r26.11.3` | 119.896 | 0 | 0 | 0 | 119.896 | | 59.578 | 10.043 | 8.662 | 5.925 | 35.689 | 0 | 3.51 |
| `r26.11.4` | 119.896 | 10.174 | 156.964 | 0 | 287.034 | | 59.578 | 60.891 | 83.200 | 24.417 | 56.295 | 2.653 | 4.89 |
| `r26.11.5` | 119.896 | 10.174 | 259.687 | 44.365 | 434.122 | 427.450 to 440.794 | 59.578 | 196.006 | 94.916 | 24.675 | 56.295 | 2.653 | 4.15 |
| `r26.11.6` | 119.896 | 10.174 | 259.687 | 192.377 | 582.134 | 564.962 to 599.306 | 85.390 | 269.206 | 94.916 | 24.675 | 104.295 | 3.653 | 4.09 |

What each release adds (*derived*, differences of the columns): `r26.11.4` 167.138 Gbp (10.174 built, 156.964 staged);
`r26.11.5` 147.088 Gbp, of which 102.723 staged and 44.365 planned (`RMD1` 35.0, `GBRMGD` 9.365); `r26.11.6` 148.012 Gbp
planned (137.512 to 158.512). The 102.723 Gbp staged for `r26.11.5` is inside DEC-79's 96 to 108 Gbp (*derived* there
from estimates), so no 20% stop applies; the coral sets came after that estimate (DEC-85).

**`r26.11.5` is an upper bound.** Each staging card computed its new bases against the root, the 19 `r26.11.4` stages and
those `r26.11.5` stages that were on the pod when it ran (the data file says which, stage by stage); none computed
against every stage before it in the chain. Card A-29's exact ownership over the chain replaces them (section 7). `RMD1`
and `GBRMGD` are estimates: card C-83 is staging them today.

## 4. Classes, and where the new sources go

The six classes of the old figure are kept, so `D` stays comparable; each has a written rule in the data file (`classes[].rule`):
isolate genomes (complete prokaryotes); uncultured prokaryotes; viruses; plasmids and mobile elements; unbinned
metagenome contigs; organelles. The new sources:

- **Single cells** (`SAG1`, `SAG2`, `SAG3`) go with the MAGs as uncultured prokaryotes: their genomes are of organisms
  not in culture, as a MAG's are; a class of their own would split one kind of source by method.
- **Giant viruses** (`GV1`, `EGV1`, `MVEUKV`) and **mirusviruses and polinton-like viruses** (`MIR1`, `MIR2`) are viruses.
- **Reef and wetland MAGs** (`RMD1`, `GBRMGD`, `WET1`) and the holobiont and niche sets (`NX1`, `KAVHIFI`) are uncultured
  prokaryotes, as `G1` and the `GART` sets are.
- **A stage that mixes kinds goes to the class that holds most of its bases.** `LR2`: MAGs 19.4 of its 21.3 Gbp (*derived*,
  card R-19's note: 21.255 less 1.882 Gbp of unbinned contigs), uncultured. `SYM2`: prokaryote bins 0.49 of 0.87 Gbp
  (*measured*, card C-78), uncultured. `KAVHIFI`: MAGs 4.607 of 4.659 Gbp (*measured*, wishlist and card C-80),
  uncultured.
- `ECE2` and `LRP1` are plasmids and mobile elements; `DPN1` uncultured; `INPH2604`, `GWV1`, `MV2` viruses; for
  `r26.11.6`, `SPX2` and `LGRT2` unbinned metagenome contigs, `O4` organelles, `REFSEQ-NEW` isolates, `LR3`, `G232-A`
  and `ASG1` uncultured.

## 5. Projections: drawn and left out

- **Drawn:** the distinct proteins of `r26.11.4`, by the hybrid design's stated formula (section 2, row 6): 265.1 M,
  252.8 to 286.2 M (*estimate*), on the figure as an open circle with its range.
- **Drawn:** the isolate-only counterfactual (as before): `r26.11.2`'s NR50 families plus 10,463 per Gbp (`C5`'s rate,
  *derived*) times the Gbp added after `r26.11.2`, through `r26.11.4`: 12.3 M (*derived*: 10.092 + 0.010463 x
  (287.034 - 74.165)).
- **Left out:** NR50 families of every unbuilt release (no stated formula, row 5); proteins of `r26.11.5` and `r26.11.6`
  (the design's 350 M and 500 M are bare *estimates*, and its yield per Gbp was stated for `r26.11.4`'s sources); the
  figure says "not projected" in those columns.

## 6. Benchmarks (not re-verified)

The values and citations are the old data file's; nothing was fetched (AGENTS §8) and none was checked again.

| Line | Value | Citation and date as the file states | Against the documents of this tree |
|---|---|---|---|
| OpenGenome 1 | 300 Gbp | Nguyen et al., *Science* 2024 (Evo), prokaryote, phage and plasmid corpus | `docs/archive/PRIOR-ART.md` cites Evo and gives OpenGenome2's 8.8 Tbp, not OpenGenome 1's size; no contradiction |
| NCBI ClusteredNR | about 286 M (90% identity, 90% coverage; full nr over 595 M) | NCBI BLAST default, Aug 2025 | `planning/search/05_external_modern.md` gives 281 M nr queries in DIAMOND's 2021 run: an older nr, not a contradiction |
| UniRef100 | 227.09 M | UniProt 2026_04 | **flag:** `05_external_modern.md` (ColabFold) gives 277.5 M UniRef100 for an earlier release; a fall of 50 M between releases is not explained by any document here and should be re-verified before the figure is shown |
| UniRef50 | 41.24 M (41,243,524) | UniProt 2026_04 | **flag:** the tree cites 53.6 M (PLMSearch 2024, `05_external_modern.md` and `planning/design/pLM_remote_system_search.md`), 51 M (2023) and 37.5 to 39 M (2021) for earlier releases; 41.24 M is below the 2024 figure: re-verify |

## 7. What changes when cards A-29 and C-83 land

- **Card A-29** writes its ownership JSON with `a29_r26_11_5_chain.py ownership --out` (each stage's `new_bases`). Copy it
  into the tree beside its note and set `"r26_11_5_ownership"` in `release_diversity.json` to its path (or pass
  `--ownership PATH` for a look first). Its `new_bases` then replace the staging notes' for every `r26.11.5` stage it
  holds, those stages become *measured*, and the footnote drops "upper bound". One command applies it:

  ```
  python scripts/plot_release_diversity.py render --ownership <path of A-29's ownership JSON> --out planning/measurements/2026-11-27_release_diversity/release_diversity
  ```

- **Card C-83** stages `RMD1` and `GBRMGD`. A-29 takes them into the chain's end (STATUS H1), so its ownership JSON
  carries them and the same command draws them as staged with measured bases. Before A-29, a C-83 note can replace the
  two estimates by hand: set each stage's `state` to `staged` and its `bases` to the note's new bases with `path` and
  `quote`, then `check`.
- A release that is built: move its stages to `released`, add the release's record with `stored_bases`, `proteins`,
  `nr50`, `nr90`, run `root` for the audit section, then `check` and `render`.

## 8. The `tower` remainder

*Measured* as `sources.yaml` records it (`plot_release_diversity.py tower`, wilkcomp02, 2026-11-27): of 158 `tower_*`
entries, 23 are taken by a stage of the figure or are reference tables (`release_diversity.json`, `tower_remainder.assigned`);
of the rest, **36 sets with recorded files hold 152,669,082,210 bytes of FASTA and archives and 12,473,780,413 bytes of
other files**, and **99 sets are listed with no file recorded**. Bases are not known: the archives are not opened, and
the audit of card R-11 found about 30 GB of what landed is not genome sequence (`2026-11-26_tower_ingest_audit.md` §2).
This is a lower bound of what is on `tower`: `sources.yaml` records the files card R-11 read on 2026-11-26, and lists
fetched since (the 20-list run of 2026-11-26, `ok=12152`, STATUS H0) are not in it. `EEMC1` (75,085 files) is being
fetched; its list's `bp` column sums to 181,337,158,637 bases (*derived*, `sources.yaml`), and it is in no release. The
figure gives the remainder as a box in panel b, in bytes, not on the Gbp axis.

## 9. The figure

Four panels and the Wilke / Okabe-Ito style are kept. Changes, and why:

- **Eight release columns and three states.** States are drawn by fill, colours stay the source classes: solid
  (released), cross-hatched (built, measured), light hatched (staged), dashed outline with no fill (planned). A second
  legend names the states. The figure is wider (16.5 in) so that eight columns keep their labels.
- **Panel b stacks by state, then by class** (released at the bottom): the measured part of each bar is one block, so
  measured, staged and planned read at a glance; panel a keeps each class in one block (its question is the share), with
  the state as the fill inside it. A whisker gives the range of a total that holds estimates; the label gives it too.
- **Panel c** draws measured releases only, the `r26.11.4` protein estimate with its range, and "not projected" where no
  formula applies; each benchmark line is labelled with its value and date.
- **Panel d** says which release each stage was measured against.
- No em dashes or en dashes in the figure's text (the test checks it); sentence case; `D` and the release-date rule in the
  footnote.

Rendered on the pod in the image at one thread (section 0); the `.svg` has no date (`metadata={"Date": None}`) and a
fixed hash salt, so it is the same on every run (the test renders twice and compares). The test does not compare the
committed image.

## 10. 2026-11-27 update (card A-32, DEC-85, DEC-87): exact 24-stage `r26.11.5`, verified UniRef benchmarks, and built wave-6 stages

- **Date and hosts:** 2026-11-27, `wilkcomp02` (root audit, 24-stage ownership, checks, 1-thread figure render) and
  `tower` (UniProt plain-text release notes). Branch `task/A-32`.

### 10.1 Wave-6 stages newly built in the root (`a32-root`)

*Measured*, `wilkcomp02`, 2026-11-27 16:39 CST, unit `a32-root` (`systemd-run --user --unit=a32-root --wait --collect
-p MemoryMax=8G -p MemorySwapMax=0 -p CPUQuota=200% -p IOSchedulingClass=idle`, `nice -n 19`, test image
`bio-box-test_faec1f4.sif`, root bound read-only, `19.9 s` wall, `38.2 s` CPU):
`plot_release_diversity.py root --root /tmp/if2453/r10/root --built V2,PL2,PL3,CG1,MGE1 --threads 2 --memory 4GB`.
Since §1 (`a31-root`), wave 6 completed `PL3`, `CG1` and `MGE1` (`after`, `check`, and `backup` complete; `SAT1` was
actively building in `sat1-1.service` and stays `staged` in the figure):

| Stage | Default-variant contigs | Rows | Default-variant bases | Against card A-26's exact ownership |
|---|---:|---:|---:|---|
| `PL3` | 745,873 | 746,483 | 16,282,680,966 | equal (`16,282,680,966` bp; `610` non-default variant rows) |
| `CG1` | 120,428 | 120,428 | 10,727,099,394 | equal (`10,727,099,394` bp) |
| `MGE1` | 13,847 | 13,994 | 148,953,649 | equal (`148,953,649` bp; `147` non-default variant rows) |

In `r26.11.4` (`287.034 Gbp` cumulative), the built portion rises from `10.174 Gbp` (`V2, PL2`) to **`37.333 Gbp`**
(`37,332,806,084` bp, *derived* sum of `V2, PL2, PL3, CG1, MGE1`) and the staged portion moves from `156.964 Gbp` to
**`129.805 Gbp`** (`129,804,807,455` bp, *derived* difference from `167,137,613,539` bp); class totals and `D = 4.89`
are unchanged.

### 10.2 UniProt Release `2026_04` verification of `UniRef50` and `UniRef100` (`DEC-87`)

Following `DEC-87` (`"i approve adding uniprot url to sources.yaml"`), entry `uniprot_relnotes` was appended to
`sources.yaml` first (`AGENTS.md` §8) and the three plain-text release notes were fetched gently one at a time on
`tower` (`curl -fsSL -A curl --proto =https --proto-redir =https`, UTC `2026-11-27T22:38:21Z` to `2026-11-27T22:42:10Z`)
into `/mnt/user/research-bulk/bio-box/sources_incoming/UNIPROT_RELNOTES/` alongside `SHA256SUMS` and `SOURCE_URLS.tsv`:

| File | Bytes (*measured*, `tower`, `wc -c`) | sha256 (*measured*, `tower`, `sha256sum`) | Published contents (UniProt Release `2026_04`, `04-Nov-2026`, `CC BY 4.0`) |
|---|---:|---|---|
| `relnotes.txt` | 1,125 | `764d77997cab2d053555fe859aabeb946bef4a93359c30681f560b84ed436e86` | `UniProtKB` `154,222,205` entries (`Swiss-Prot` `576,073`, `TrEMBL` `153,646,132`); `UniRef100` `227,086,028` entries; `UniRef90` `125,724,408` entries; `UniRef50` `41,243,524` entries |
| `uniref50.release_note` | 311 | `9ff782280c206c85067aeeb3f005eca659f887b153c1f282b41e725d9199c97c` | `Release: 2026_04, 04-Nov-2026`; `Number of clusters: 41,243,524` (`41.243524 M`, *derived* `÷ 10^6`) |
| `uniref100.release_note` | 312 | `f3ddf0e7fdc51c5d526cfb7fc412164344883792bb61300628eb95fce4dca12b` | `Release: 2026_04, 04-Nov-2026`; `Number of clusters: 227,086,028` (`227.086028 M`, *derived* `÷ 10^6`) |

This resolves both flags of §6: `227,086,028` (`227.1 M`) for `UniRef100` and `41,243,524` (`41.2 M`) for `UniRef50`
are the exact UniProt `2026_04` (`04-Nov-2026`) counts (smaller than 2023-2024 releases because UniProt pruned redundant
proteomes, leaving `154,222,205` `UniProtKB` entries). Both benchmarks in `release_diversity.json` now carry
`re_verified: true`, exact million-scaled values (`227.086028` and `41.243524`), and their `sha256` digests.

### 10.3 Exact 24-stage `r26.11.5` ownership (`a32-own`) and redrawn totals

With card C-83's staged `RMD1` (`38,301,644,358` new bp, *measured*) and `GBRMGD` (`10,238,417,378` new bp, *measured*)
appended to `R26_11_5_CHAIN` and exact 24-stage prefix ownership committed at
[`planning/measurements/2026-11-27_release_diversity/r26_11_5_ownership.json`](2026-11-27_release_diversity/r26_11_5_ownership.json)
(`a32-own`, `0` problems, `0` discrepancies, `151,204,587,185` new bp across all 24 stages; see
`2026-11-27_a29_r26_11_5_confs.md` §8), `r26_11_5_ownership` in `release_diversity.json` points to the committed
ownership JSON, `r26_11_5_upper_bound` is `false`, and all 24 stages of `r26.11.5` are drawn as `staged` (*measured*
exact new bases, `0` planned Gbp in `r26.11.5`). Updated cumulative columns (*derived* by `plot_release_diversity.py
check`):

| Column | Released | Built | Staged | Planned | Total | Range of the total | Isolate | Uncultured | Virus | Plasmid | Metagenome | Organelle | D |
|---|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|
| `r26.11.4` | 119.896 | 37.333 | 129.805 | 0 | 287.034 | | 59.578 | 60.891 | 83.200 | 24.417 | 56.295 | 2.653 | 4.89 |
| `r26.11.5` | 119.896 | 37.333 | 281.009 | 0 | 438.238 | | 59.578 | 200.167 | 94.871 | 24.675 | 56.295 | 2.653 | 4.12 |
| `r26.11.6` | 119.896 | 37.333 | 281.009 | 148.012 | 586.250 | 575.750 to 596.750 | 85.390 | 273.367 | 94.871 | 24.675 | 104.295 | 3.653 | 4.07 |

