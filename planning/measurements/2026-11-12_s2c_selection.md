# S2C selection: the ≤20-per-species complete genomes (card C-25, step 1; wilkcomp02, 2026-11-12)

Step (1) of card C-25 only. Nothing was fetched or built. Every number is *measured* on wilkcomp02 on 2026-11-12 (CST)
by the run below unless it is labelled *derived* (formula given) or *estimate*. Nothing was written under `/stor`;
the GTDB archive there was only read.

## Result in one paragraph

By quality order alone, the S2C set is **25,458 genomes and 106,409,469,843 bp**, exactly the plan's 25,458 and
106.4 Gbp (0% difference). **All 11,973 phase-1 genomes are in it** (each is its species' rank 1). The remainder
(the set minus phase 1's 11,973) is **13,485 genomes and 55,910,526,156 bp**, exactly the card's 13,485 and 55.9 Gbp
(0%). No number differs from the card's by more than 20%. The ANI-diversity step that is missing can only act in
the 243 species with more than 20 complete genomes. There it could change at most 3,368 of the 25,458 selected
genomes (13.2%), and so at most 3,368 of the 13,485 remainder genomes (25.0%), if each species' best genome is kept.
Whatever genomes it picked, S2C would stay between 106.01 and 107.28 Gbp.

## Inputs and commands

- Metadata: the durable copy `/stor/scratch/Wilke/if2453/all-the-data/provisional/s1_inputs/gtdb_r226_metadata.tar`.
  Its sha256 is `d686a3bf3541eeb4610798da39a0174c0661d15a3485ef19d850e10aab77a817`, equal to `sources.yaml` and the
  archive's `MANIFEST.sha256`. It was unpacked to `/tmp/if2453/c25/gtdb/`. `bac120_metadata_r226.tsv.gz` has sha256
  `45f1614a…bcbf` and `ar53_metadata_r226.tsv.gz` has `f12310f2…b6e5`, both equal to `sources.yaml` (full digests
  in `2026-11-12_s2c_selection/selection.json`).
- Phase 1's list: `planning/measurements/2026-11-11_r50_inputs/selection/all.tsv` (11,973 genomes, sha256
  `e68d2576067d4d83f77a50384865de5470e6dff1a42716a06148cdf26c0efeaf`), copied to `/tmp/if2453/c25/phase1/all.tsv`.
- Scripts: `scripts/stress/build_stage_pod/c25_select.py`, with its wrapper `c25_select.sh`. The script uses phase 1's
  metadata reader. That reader was moved, unchanged in behaviour, into the function `complete_genomes` of
  `r50_select.py`. The wrapper reruns `r50_select.py` and checks that its `all.tsv` still equals phase 1's
  committed file (`cmp`: equal).
- Run (host `python3` 3.12.3, standard library only), 16:29:05 to 16:30:07 CST:

  ```
  systemd-run --user --unit c25-select --collect -p MemoryMax=8G -p MemorySwapMax=0 -p CPUQuota=100% \
    -p StandardOutput=append:/tmp/if2453/c25/logs/select.log -p StandardError=append:/tmp/if2453/c25/logs/select.log \
    nice -n 19 sh /tmp/if2453/c25/bin/c25_select.sh
  ```

  Inside the wrapper: `python3 -I c25_select.py /tmp/if2453/c25/gtdb /tmp/if2453/c25/phase1/all.tsv
  /tmp/if2453/c25/out --cap 20 --phase1-not-served <the 7 accessions of the r50 inputs record>`. It took 29.8 s of
  wall time at one CPU, with a peak RSS of 252,672 kB (`/usr/bin/time -v`).

## The rule

- **Complete**: `ncbi_assembly_level = 'Complete Genome'`, and no other field is used. This is phase 1's definition
  (`r50_select.py`) and the plan's (`scripts/stress/gtdb_counts.py`). 49,324 genomes qualify, in 11,973 species.
- **Species**: the full `gtdb_taxonomy` string.
- **Order within a species**: `checkm2_completeness − 5 × checkm2_contamination`, descending. A value that does not
  parse would sort last, but none of the 49,324 lacks a value or has a non-finite one (`complete_genomes_checkm2_unparsed` 0,
  `complete_genomes_checkm2_nonfinite` 0). Every tie is broken by GTDB accession, ascending, with its `RS_`/`GB_`
  prefix. This orders the same way as the bare accession, because `GB_` sorts before `RS_` and `GCA` before `GCF`.
- **Cut**: ranks 1 to 20 are selected. No ANI step is applied (skani is not in the image and was not installed).
- Rank 1 is phase 1's pick. The script rebuilt the rank-1 rows in phase 1's columns, and they equal phase 1's
  `all.tsv` byte for byte (`rank1_equals_file: true`). The whole rule is the plan's SQL in `gtdb_counts.py`
  (`row_number() OVER (PARTITION BY gtdb_taxonomy ORDER BY comp - 5 * cont DESC NULLS LAST, accession) <= 20`),
  and it gives the plan's counts exactly (next table).

## The set

| Count | This run | Card / plan | Difference |
|---|---|---|---|
| S2C genomes | **25,458** | 25,458 | 0% |
| S2C bp (GTDB `genome_size`) | **106,409,469,843** | 106.4 Gbp | 0% |
| S2C contigs, proteins (GTDB `contig_count`, `protein_count`) | 56,387; 97,651,068 | 56,387; 97.7 M | 0% |
| phase 1 (rank 1) genomes, bp | 11,973; 50,498,943,687 | 11,973; 50.5 Gbp | 0% |
| remainder genomes | **13,485** | about 13,485 | 0% |
| remainder bp | **55,910,526,156** | about 55.9 Gbp | 0% |
| remainder contigs, proteins (GTDB) | 29,332; 51,597,124 | | |

Other counts for the S2C set: 706 archaea. Translation table 11 / 4 / 25 / missing: 23,419 / 545 / 5 / 1,489.
9,557 of the selected genomes are GTDB species representatives (`gtdb_representative = t`).

The remainder has 164 archaea. Its translation tables (11 / 4 / missing) are 12,853 / 354 / 278. The missing ones
are GenBank-only accessions; phase 1 gave such genomes code 11.

Files, on the pod only (`/tmp` is cleared on reboot; the script rebuilds them byte for byte from the archive):

| File | Rows | Bytes | sha256 |
|---|---|---|---|
| `/tmp/if2453/c25/out/s2c_selected.tsv` | 25,458 | 5,231,154 | `9a67b230bc898874a8aed6f566a253a3207bd71f63ce412a1dfcf2edd0a7beb5` |
| `/tmp/if2453/c25/out/remainder.tsv` | 13,485 | 2,767,119 | `7ae727a8fd5fa6be223584901201871d2bbea982cd45cca95d5b842091968c51` |
| `/tmp/if2453/c25/out/by_phylum.tsv` (committed) | 80 | 3,666 | `3b3e831c1f29df3acc94f8715bca55c08277f213e7fd1ce498f0060daae931b9` |
| `/tmp/if2453/c25/out/selection.json` (committed) | | 4,581 | `04dcd5bc005e282a23cd398e9e1978db3aa05656839e9c1db497c6a1637df194` |

Columns of the two lists: phase 1's (`accession`, `genome_size`, `contig_count`, `protein_count`,
`translation_table`, `gtdb_taxonomy`, `gtdb_accession`), then `species_rank`, `species_complete_genomes`,
`quality_score`, `checkm2_completeness`, `checkm2_contamination`, `gtdb_representative`. Both are sorted by accession.

## Phase 1 against the set

- **Every one of the 11,973 phase-1 genomes is in the set** (`not_in_selection: []`). Each is rank 1 of its species,
  and the cap keeps rank 1. Of the 11,966 that were actually built, none is missing either.
- The 7 that NCBI did not serve (`GCA_900637225.1`, `GCA_902153285.2`, `GCA_902502775.2`, `GCA_902502895.2`,
  `GCA_963667335.1`, `GCA_963921365.2`, `GCF_040719245.2`) are in the set as rank 1 of their species. They are not
  in the remainder.
  - Five of them are the only complete genome of their species. Unless a later fetch gets them, those 5 species
    have no genome in S2C.
  - `GCA_902502775.2` has 3 complete genomes in its species and `GCA_902502895.2` has 2. The others of those two
    species (2 and 1 genomes) are in the remainder.
  - Whether the 7 are tried again in step (2) is a fetch question.
- **Under an ANI rule.** If such a rule starts from the best-quality genome, phase 1 stays a subset. If it could drop
  the best-quality genome of an over-cap species, at most 243 phase-1 genomes could fall out (one per species with
  more than 20 complete genomes; *derived*).

## The remainder

By how many complete genomes the species holds (`by_species_size` in `selection.json`):

| Complete genomes in the species | Species | S2C genomes | S2C bp | Remainder genomes | Remainder bp |
|---|---|---|---|---|---|
| 1 | 8,785 | 8,785 | 37,006,319,407 | 0 | 0 |
| 2 | 1,424 | 2,848 | 12,227,101,087 | 1,424 | 6,150,018,913 |
| 3–5 | 970 | 3,532 | 15,480,667,068 | 2,562 | 11,308,171,909 |
| 6–10 | 370 | 2,777 | 12,084,497,529 | 2,407 | 10,495,929,606 |
| 11–20 | 181 | 2,656 | 11,476,412,465 | 2,475 | 10,718,444,760 |
| 21–50 (capped at 20) | 147 | 2,940 | 11,046,339,793 | 2,793 | 10,501,503,956 |
| 51–100 (capped) | 47 | 940 | 3,453,265,883 | 893 | 3,282,395,369 |
| 101 or more (capped) | 49 | 980 | 3,634,866,611 | 931 | 3,454,061,643 |
| all | 11,973 | 25,458 | 106,409,469,843 | 13,485 | 55,910,526,156 |

The species with the most complete genomes are *E. coli* (4,007), *K. pneumoniae* (2,312), *S. enterica* (1,833),
*S. aureus* (1,475), *P. aeruginosa* (920), *B. pertussis* (792), *A. baumannii* (750), *C. jejuni* (522),
*M. tuberculosis* (505) and *B. velezensis* (496). Each contributes 20.

By GTDB phylum: 42 of the 80 phyla with a complete genome have a genome in the remainder. Full table:
`2026-11-12_s2c_selection/by_phylum.tsv`. The largest:

| Phylum | Species | Remainder genomes | Remainder bp | Remainder proteins (GTDB) |
|---|---|---|---|---|
| Pseudomonadota | 4,508 | 6,395 | 30,227,246,671 | 27,791,083 |
| Bacillota | 1,993 | 3,547 | 10,886,358,675 | 10,488,652 |
| Actinomycetota | 1,992 | 1,563 | 8,761,524,307 | 7,877,056 |
| Bacteroidota | 1,067 | 631 | 2,514,698,021 | 2,160,288 |
| Campylobacterota | 217 | 285 | 519,852,611 | 516,015 |
| Spirochaetota | 145 | 262 | 631,556,455 | 563,383 |
| Fusobacteriota | 49 | 137 | 321,448,099 | 301,580 |
| Cyanobacteriota | 225 | 114 | 473,781,902 | 450,778 |
| Chlamydiota | 29 | 105 | 117,579,461 | 100,479 |
| Halobacteriota | 285 | 67 | 244,410,417 | 234,955 |
| 32 other phyla | | 379 | 1,212,069,537 | 1,112,855 |

The last row is *derived* as the total minus the rows above.

**Expected proteins.** This is the number of gene calls; distinct proteins will be fewer, because strains of one
species share many identical proteins (step 3 measures that).

- GTDB's own `protein_count` over the 13,485 remainder genomes: **51,597,124** (*measured*, metadata), 0.9229 per kbp.
- At phase 1's built density: **about 51.1 M genes** (*derived*: 55,910,526,156 bp × 10,029,713 genes /
  10,965,710,288 bp = 51,138,186). The density inputs are stages P1 (927,920 genes, 1,016,625,803 bp;
  `2026-11-10_build_stage_pod.md` and STATUS C-10) and P3 (9,101,793 genes, 9,949,084,485 bp; `2026-11-11_build_r50.md`).
- The same stages hold 9,884,583 distinct proteins (counted per stage), which gives about 50.4 M at that density
  (*derived*, same formula). That ratio of distinct proteins to genes (98.6%) comes from stages with one genome per
  species. The remainder adds up to 19 more strains per species, so its distinct count will be lower. Step (3)
  measures by how much.
- At phase 1's GTDB density (46,053,944 proteins / 50,498,943,687 bp): about 51.0 M (*derived*).

## What the missing ANI step could change

The plan orders by CheckM2 quality and then picks for ANI diversity. A species with 20 or fewer complete genomes
keeps all of them under either rule. Only species with more than 20 can differ:

| Count | Value |
|---|---|
| species with more than 20 complete genomes | **243** |
| complete genomes they hold | **28,726** |
| genomes selected from them (20 each) | 4,860 (of which 243 rank 1, in phase 1) |
| remainder genomes from them | 4,617 (17,237,960,968 bp; 34.2% of the remainder's genomes, 30.8% of its bp, *derived*) |
| most selected genomes that could differ, best genome kept: Σ min(19, n − 20) | **3,368** (13.2% of S2C, 25.0% of the remainder, *derived*) |
| most that could differ, best genome not kept: Σ min(20, n − 20) | 3,497 |
| S2C bp over every choice that keeps the best genome (min, quality order, max) | 106,013,943,010; 106,409,469,843; 107,275,959,885 |
| remainder bp over the same choices (min, max) | 55,514,999,323; 56,777,016,198 (−0.7% and +1.5%, *derived*) |

The genome count of S2C and of the remainder does not depend on the rule. Their bp change by at most the range
above. The only thing an ANI step changes is *which* strains of the 243 large species are kept.

**How much the quality order decides.**

- In 102 of the 243 species, the genome at rank 20 and the one at rank 21 have the same score, so the cut falls
  inside a tie. 485 selected genomes were chosen inside such ties by accession alone.
- In 5 species, the top 21 all share one score.
- In these large species, quality mostly separates genomes by small contamination differences. Beyond that the
  pick is by accession number, which favours older submissions. It is not a diversity choice.

**Diversity proxies in the metadata.**

- The R226 metadata has no ANI column (`metadata_columns_named_ani: []`). `gtdb_genome_representative` names the
  species representative but carries no distance to it.
- Proxies that would need no new tool:
  - prefer distinct `ncbi_bioproject` / `ncbi_biosample`, so that one submitter's near-identical isolates are not
    all taken;
  - spread over `gc_percentage` and `genome_size`;
  - spread over `ncbi_country`, `ncbi_isolation_source` or `ncbi_date`.
- None of these measures ANI. They could reasonably replace the accession tie-break at the cut (485 genomes), but
  that is a change to the plan's rule, so it is **not adopted** here.
- GTDB's separate per-species ANI files are not in `sources.yaml` and were not fetched.

## Test

`tests/r1/test_c25_select.py` (covers D-01, D-04; ring R1) runs the script on hand-made metadata (ties at equal
scores, missing CheckM2 values, non-complete genomes, a species at exactly the cap, `GB_`/`RS_` prefixes, an
archaeal file) plus 400 seeded random rows, with `--cap 3`. It checks:

- the result against a brute-force oracle: a genome's rank is one plus the number of genomes of its species that
  beat it, by explicit pairwise comparison, with no sort;
- that the output is byte-identical when the input rows are shuffled;
- the phase-1 reasons and the not-served list;
- the ANI counts and the bp bounds against an enumeration of every allowed subset (`itertools.combinations`).

It failed first (the script did not exist; exit status 2) and passes on the pod in the image:
`sh scripts/pod_pytest.sh C-25 tests/r0/test_runbook.py tests/r1/test_c25_select.py tests/r1/test_make_inputs.py`
gave `45 passed`.

## Pod use

- One unit, `c25-select`: about 62 s at one CPU, nice 19, 8 G cap.
- Two short runs of the test in the image through `pod_pytest.sh`.
- One `awk | sort` over the output to list the largest species.
- Nothing ran on the laptop except one line of arithmetic for the *derived* numbers above.
- Written only under `/tmp/if2453/c25` and `/tmp/if2453/ci/C-25`. `/tmp/if2453/ladder/genomes/P{1,2}/lock.json`
  was read (by `r50_select.py`'s recheck). Nothing under `/tmp/if2453/r10`, `/tmp/if2453/r50` or `/stor` was
  changed, and no network request was made.
