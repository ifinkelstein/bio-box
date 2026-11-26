# MetaVR's `>= 30 kb` spillover (`MV2`) and eukaryotic-virus slice (`MVEUKV`): selection, one stream, staging (card C-79, DEC-79)

- **Date and host:** 2026-11-26, `wilkcomp02`. Every number carries a label (AGENTS §6): *measured* (on `wilkcomp02`,
  2026-11-26, by the command named in the section, unless stated otherwise), *derived* (formula and inputs) or
  *estimate*. Counts quoted from card C-62's note say so.
- **Script:** `scripts/stress/build_stage_pod/c79_metavr_mv2_eukv.py` (`profile`, `slice`, `stage`, `ani`); tests
  `tests/r1/test_c79_metavr_mv2_eukv.py`; templates `scripts/stress/build_stage_pod/stages/templates/MV2.conf` and
  `MVEUKV.conf`.
- **Read only:** `/stor/scratch/Wilke/if2453/fragdef/db/metavr/` (`METAVR_main_table.parquet` 1,087,432,195 B,
  `METAVR.fna.bgz` 71,552,362,108 B, the sizes card C-62 recorded, *measured* by `ls -la`), C-62's
  `MV1_left_for_later.tsv.gz` (298,672,948 B) on `/stor/scratch/Wilke/if2453/bio-box-staging/MV1/`, the data root
  `/tmp/if2453/r10/root` and `/tmp/if2453/s2/inputs/` (bound `:ro` in the image). Nothing was written to the root,
  to `/tmp/if2453/s2/inputs/` or to `/tmp/if2453/stages/`; nothing under `/stor` was overwritten.
- **Written:** `mv2.fna` and `mveukv.fna` (each over 64 MB) to `/stor/scratch/Wilke/if2453/bio-box-staging/{MV2,MVEUKV}/`;
  metadata to `/tmp/if2453/c79/{MV2,MVEUKV}/`; profiles and the ANI sample to `/tmp/if2453/c79/work/`.

## 1. Summary

| | `MV2` | `MVEUKV` |
|---|---:|---:|
| Rule | the 98,920 `mv2_spillover` rows of `MV1_left_for_later.tsv.gz`, `mv2_rank` 1..98,920 | one representative per untouched vOTU of *Nucleocytoviricota*, *Maveriviricetes*, *Polintoviricetes*, `viral_confidence` High, length >= 5,000 bp (§3) |
| UViGs selected (table `length` sum) | 98,920 (3,909,189,032 bp) | 105,902 (1,257,254,288 bp) |
| FASTA records written | 101,573 | 106,087 |
| Genomes listed (DEC-27 repeated sets) | 98,830 (90) | 105,877 (25) |
| Distinct contigs (bp) | 101,483 (3,905,501,982) | 106,061 (1,256,895,787) |
| DEC-65 exact DTRs trimmed (bp) | 669 (77,784) | 133 (14,895) |
| Circular genomes | 669 | 133 |
| Held by the root or `r26.11.4` stages | 1,164 contigs, 44,501,988 bp | 63 contigs, 886,270 bp |
| New after the root and `r26.11.4` | 100,319 contigs, 3.861 Gbp | 105,998 contigs, 1.256 Gbp |
| Also held by `GV1` (exact `contig_uid`) | 4,376 contigs, 98,597,781 bp | 3,812 contigs, 55,488,128 bp |
| New after the root, `r26.11.4`, `GV1` (and `MV2` for `MVEUKV`) | 95,943 contigs, 3.762 Gbp | 102,186 contigs, 1.201 Gbp |
| `s2_check_inputs.py` (with `objects_staging.json`) | exit 0, `problems: []`, `gene_call_mode` `{"virus": "gv-meta"}` | exit 0, `problems: []`, `gene_call_mode` `{"virus": "gv-meta"}` |
| Conf template | `DEPTH=called`, `CRISPR=yes`, `AFTER=MGE4` | `DEPTH=annotated`, `LIBRARY=pfam 38.2 gathering`, `CRISPR=yes`, `AFTER=MV2` |

All *measured* (`c79_metavr_mv2_eukv.py stage`, report `/tmp/if2453/c79/work/stage_report.json`). The `MV2` count
and bases equal card C-62's `mv2_spillover` remainder (98,920 UViGs, 3,909,189,032 bp, *measured* there): 0 % apart.
`MVEUKV`'s 1.257 Gbp is inside the plan's 0.8 to 1.8 Gbp (*estimate*, plan §4.7) once the rule of §3 is applied;
the card's rule without it was not (§3). The GV1 ANI sample (§6): **31.7 %** of the sampled `MVEUKV` representatives
meet the vOTU threshold with a GV1 contig (rule `shorter`), against 3.6 % exact `contig_uid` overlap.

## 2. The table's schema and taxonomy column (read before the rule; `profile`)

`METAVR_main_table.parquet` (*measured*, DuckDB `DESCRIBE`, `/tmp/if2453/c79/work/profile.json`): 24,435,662 rows,
34 columns, as `README.upstream.txt` states (`uvig`, `taxon_oid`, `discovery`, `discovery_doi`, `uvig_topology`,
`is_concatemer`, `quality`, `checkv_contamination`, `completeness_method`, `completeness`, `genomad_score`, `cds_count`,
`trna_count`, `total_gene_count`, `length`, `votu`, `n_virus_hallmarks`, `viral_confidence`, `ictv_taxonomy`,
`ictv_taxonomy_method`, `genome_type`, `source`, `genetic_code`, `host_taxonomy`, `host_taxonomy_method`, `Latitude`,
`Longitude`, `Isolation_Country`, `Geographic_Location`, `Ecosystem`, `Ecosystem_Category`, `Ecosystem_Type`,
`Ecosystem_Subtype`, `Specific_Ecosystem`).

- The taxonomy column is `ictv_taxonomy`, a lineage of eight ranked tokens
  `r__<realm>;k__<kingdom>;p__<phylum>;c__<class>;o__;f__;g__;s__` (empty names where unassigned).
- `ictv_taxonomy_method`: geNomad 22,491,460; none 767,477 (`ictv_taxonomy` NULL or empty); MMSeqs Taxonomy 748,381;
  RVMT 346,643; vOTU 49,442; Imported from NCBI 18,355; GVClass 12,994; Viroids 910.
- **No row of `ictv_taxonomy` or `ictv_taxonomy_method` names Mirusviricota** (0 rows matching `%mirus%`, case
  ignored; *measured*, `slice`). This is a fact about the table's taxonomy column, which follows the ICTV taxonomy
  geNomad and the other methods assign; it is not a finding that MetaVR holds no mirusviruses, which would sit
  unclassified or under another lineage here. `MVEUKV` cannot take them by this column.
- *Tectiliviricetes* (in *Preplasmiviricota* with the virophages and PLVs) holds tectiviruses, corticoviruses,
  turriviruses and adenoviruses (lineages listed in `profile.json`): bacterial, archaeal and animal viruses. It is not
  taken; the eukaryotic-virus classes of *Preplasmiviricota* named by the card are *Maveriviricetes* (virophages) and
  *Polintoviricetes* (polintons, PLVs).

## 3. The slice's size, and the rule

**The card's rule as first written** (one representative per vOTU of the named taxa, by C-62's order, among rows
that are not concatemers and have contamination at most 10 %, in vOTUs none of whose UViGs MV1 or MV2 took)
counted (*measured*, `slice`, `/tmp/if2453/c79/work/slice.json`):

| Taxon | vOTUs | bp |
|---|---:|---:|
| *Nucleocytoviricota* | 799,347 | 6,517,736,078 |
| *Maveriviricetes* | 39,498 | 274,805,303 |
| *Polintoviricetes* | 6,677 | 47,566,493 |
| *Mirusviricota* | 0 | 0 |
| **All** | **845,522** | **6,840,107,874** |

By `viral_confidence` (same file): High 128,345 vOTUs and 1,336,822,492 bp (*Nucleocytoviricota* 116,489 and
1,246,923,385; *Maveriviricetes* 10,164 and 77,231,880; *Polintoviricetes* 1,692 and 12,667,227), Low 717,177 and
5,503,285,382. By length: 591,921 vOTUs and 5,747,876,790 bp at 5 kb and more; 191,254 and 3,003,810,511 at 10 kb and
more; 37,927 and 918,214,337 at 20 kb and more. By quality: Low-quality 706,469, Not-determined 139,052 and one
Medium-quality. The 6.84 Gbp is 3.8 times the top of the plan's 0.8 to 1.8 Gbp (*derived*: 6.840 / 1.8), so the card
stopped and reported it (AGENTS §6, §9).

**The rule decided** (the orchestrating session, 2026-11-26, from that table; the module docstring states it): of
those rows, only `viral_confidence = 'High'` and `length >= 5,000` (DEC-60's floor for metagenome contigs) qualify,
and each untouched vOTU with a qualifying row gives one, by C-62's order (quality tier, capped completeness,
contamination, length, `uvig`); ranked by taxon (*Nucleocytoviricota*, *Maveriviricetes*, *Polintoviricetes*), then
the same key. The filter is applied before the representative is chosen, so a vOTU whose best row is of low
confidence but which holds a qualifying row is taken by that row. `DEPTH=annotated`.

Result (*measured*, `stage`):

| Taxon | vOTUs | bp (table `length`) |
|---|---:|---:|
| *Nucleocytoviricota* | 96,780 | 1,178,212,188 |
| *Maveriviricetes* | 7,862 | 68,075,114 |
| *Polintoviricetes* | 1,260 | 10,966,986 |
| **`MVEUKV`** | **105,902** | **1,257,254,288** |

By quality: Low-quality 105,230 (1,250,991,466 bp), Not-determined 672 (6,262,822 bp). Untouched vOTUs with a row of
the three taxa: 887,189 (*measured*); selected + remainder (§7) = 105,902 + 486,502 + 253,118 + 41,667 = 887,189
(*derived*, the script refuses otherwise). The `slice` table above counts per taxon (a vOTU with rows of two taxa
counts in each: 888,879 untouched by taxon against 887,189 distinct, *measured*).

## 4. One stream, DEC-65 and DEC-27

`c79_metavr_mv2_eukv.py stage` read `METAVR.fna.bgz` once (`bgzip -dc` in the image) for both stages: 24,635,130
records scanned, 207,660 written (*measured*), the record count C-62 measured for the whole file. 0 selected UViGs
missing from the FASTA, 0 length mismatches between the table and the FASTA for either stage (*measured*).

- **DEC-65 (1)** on single-record `Direct terminal repeat` UViGs, the code of `c62_metavr_inputs.apply_dec65`:
  `MV2` exact 669 (77,784 bp trimmed, staged circular), partial 13 (linear, untrimmed); `MVEUKV` exact 133 (14,895
  bp), partial 9, low complexity 2.
- **DEC-27:** `MV2` 90 repeated contig sets (`repeated_records.tsv`), 386 multi-record (Virus MAG) genomes;
  `MVEUKV` 25 repeated sets, 49 multi-record genomes, 1 contig in two genomes that agree on kind, code and topology
  (`s2_check_inputs.py` accepts it).
- **Order:** `genomes.json` and `records.tsv` in rank order (`mv2_rank` for `MV2`, the rank of §3 for `MVEUKV`); the
  FASTA in the archive's order. `kind="virus"`, `genetic_code` 11 (the gv-meta caller chooses its own code),
  `assembly_quality` `vOTU`.
- **Objects:** `objects.json` names `/tmp/if2453/s2/inputs/<STAGE>/<file>` (card A-29 copies there, and A-26's
  rewrite keeps such a file byte for byte); `objects_staging.json` names the file on `/stor`. Both
  `redistributable: false` (JGI IMG-derived, as MV1).

## 5. Staged files

| File | Where | Bytes | sha256 |
|---|---|---:|---|
| `mv2.fna` | `/stor/scratch/Wilke/if2453/bio-box-staging/MV2/` | 3,915,276,060 | `5a706ef42466f73488d5b9914a8a9a5ded42bcaf715b9902d63d6802e3767bb7` |
| `genomes.json` | `/tmp/if2453/c79/MV2/` | 13,559,134 | `136fc16f335d5345a8b12db191d4476d4e9d20e18b4407e0dcbf1a6cc1054fcb` |
| `objects.json` | `/tmp/if2453/c79/MV2/` | 568 | `68858a4f4fc3b7556e1f4506ba8c0a807bc9d5d150e094cbf80ee6541d998ee1` |
| `objects_staging.json` | `/tmp/if2453/c79/MV2/` | 589 | `b2ed288a622f4120ae6aecf7125d0c5ae57a254975da5c36f12e260432eddb29` |
| `records.tsv` | `/tmp/if2453/c79/MV2/` | 31,056,804 | `6296614defe56a903c467c79eae555d7533fd18588138221dd3aaef3ab7f3199` |
| `repeated_records.tsv` | `/tmp/if2453/c79/MV2/` | 8,771 | `3628f7d9df58051e7f62176924a90b5f6206895bc1174c85d26780c1a0b326a0` |
| `shared_records.tsv` | `/tmp/if2453/c79/MV2/` | 83,098 | `a52686dcdb5f06c6dee79b39655e11c71b3b304fa5d713f75bf95a33e28be052` |
| `later_records.tsv` | `/tmp/if2453/c79/MV2/` | 329,306 | `a4d40805dcc1608b5a551a0b1d2d8ab5b0366411441e07115dedb351dfc74a4c` |
| `trims.tsv` | `/tmp/if2453/c79/MV2/` | 120,091 | `766c1dadb60a67ebd5260fe7ded6ccb4b17446d5912882320d0a56586e377a4f` |
| `mveukv.fna` | `/stor/scratch/Wilke/if2453/bio-box-staging/MVEUKV/` | 1,263,638,675 | `8ec7b906090ccbf82258d3e323e0bd8d3c360831c9c442bf0cf6b46c8112151f` |
| `genomes.json` | `/tmp/if2453/c79/MVEUKV/` | 14,410,454 | `a9b7ffb8760e323c5cba91038d383deef0ea4df8d90027d891a1773b1246d3d2` |
| `objects.json` | `/tmp/if2453/c79/MVEUKV/` | 577 | `bcf00cad7e723caa75d88b65def4b5a60f13f90f09344bf9a8c60e8840afedd0` |
| `objects_staging.json` | `/tmp/if2453/c79/MVEUKV/` | 598 | `c6e11c43dc87305673fcf137457e463bcd053a9069906994593b82fb29a632d2` |
| `records.tsv` | `/tmp/if2453/c79/MVEUKV/` | 31,366,099 | `95c82f820fe40313776b5214aa25fb94d0b0de55387932fec82b5d85566fd45f` |
| `repeated_records.tsv` | `/tmp/if2453/c79/MVEUKV/` | 2,456 | `4fa7a6add3dd847f98d6f4f5e68de35d087beadd3876ee62aa7582a49b6b36cf` |
| `shared_records.tsv` | `/tmp/if2453/c79/MVEUKV/` | 4,540 | `65ce9e99fdbb8114b7ba0ff204acee0d0a1d92205fdfb0b64eceda19a70912c1` |
| `later_records.tsv` | `/tmp/if2453/c79/MVEUKV/` | 275,104 | `2c4ec4c1da5e56f3cce8ef44cf8519ebef90cd884bb30a3bdfbb4c53ccdac518` |
| `trims.tsv` | `/tmp/if2453/c79/MVEUKV/` | 23,901 | `9ee9e639e4cf083197d3ff6e91c1e632f2f990e4f66c582ed5e505b7a43fd8ca` |
| `MVEUKV_left_for_later.tsv.gz` | `/tmp/if2453/c79/MVEUKV/` | 10,138,171 | `6210883be377598b095513a707e850dd4ecba8f7d19531c1d95fe98a735b6efe` |

All *measured* (sha256 by the script, `ids.sha256_file`). The conf templates' `SUM` lines are the rows for
`genomes.json`, `objects.json` and the FASTA. The metadata files are under 64 MB and stay on `/tmp` (AGENTS §7.3).

## 6. Overlap

**Exact `contig_uid`, by build order** (*measured*, `stage`: the root's 13 stored stages `S1, P1..P6, H1, PL1, C5,
V1, G1, SPX1` from `contig_index` with `is_default_variant`, not the `V2` partition being built; then the 19
`r26.11.4` stages from their `genomes.json` under `/tmp/if2453/s2/inputs/` in `a26_r26_11_4_chain.R26_11_4_CHAIN`
order; the first holder owns). These are the templates' `EXPECT=held_by` lines:

| Holder | `MV2` contigs | `MV2` bp | `MVEUKV` contigs | `MVEUKV` bp |
|---|---:|---:|---:|---:|
| `S1` | 1 | 35,082 | 0 | 0 |
| `PL1` | 3 | 94,564 | 0 | 0 |
| `V1` | 1,092 | 41,594,192 | 7 | 141,142 |
| `G1` | 15 | 678,508 | 0 | 0 |
| `V2` | 17 | 669,923 | 17 | 232,583 |
| `PL2` | 1 | 31,073 | 0 | 0 |
| `PL3` | 1 | 57,461 | 0 | 0 |
| `MGE2` | 2 | 64,145 | 0 | 0 |
| `MV1` | 13 | 654,018 | 39 | 512,545 |
| `LR1` | 15 | 477,865 | 0 | 0 |
| `GART1` | 3 | 108,876 | 0 | 0 |
| `GART2` | 1 | 36,281 | 0 | 0 |
| **All** | **1,164** | **44,501,988** | **63** | **886,270** |

Every other stage of the root and of `r26.11.4` holds 0 of either. **Later in the `r26.11.5` chain** (plan §5:
`GV1` before `MV2`, `MV2` before `MVEUKV`), for card A-29: `GV1` holds 4,376 contigs (98,597,781 bp) of `MV2` and
3,812 (55,488,128 bp) of `MVEUKV` that no earlier stage holds; `MV2` holds 0 of `MVEUKV`. They are listed in each
stage's `later_records.tsv`. `MIR1` is not on the pod (plan §1); its overlap is card A-29's, by chain order.

**`MVEUKV` against `GV1` at the vOTU threshold, on a sample** (*measured*, `c79_metavr_mv2_eukv.py ani`, card
C-48c's method through `c39_ani`: `skani dist -c 30 -m 200 --min-af 0 --qi --ri` from the image, 3 threads; the
threshold ANI >= 95 % over >= 85 % of the shorter sequence, rule `shorter`; `both` and `near_identical` as `c39_ani`).
The sample is every 25th of the 105,828 single-record, non-repeated `MVEUKV` genomes in FASTA order: 4,233 queries,
49,292,924 bp; the reference is all of `GV1`'s FASTA (288,536 records, 4,576,307,473 bp; skani skips records under
500 bp). `/tmp/if2453/c79/work/ani/summary.json`, `a.tsv`.

| Queries | `shorter` | `both` | `near_identical` |
|---|---:|---:|---:|
| all 4,233 | **1,342 (31.7 %)** | 300 (7.1 %) | 225 (5.3 %) |
| *Nucleocytoviricota* 3,853 | 1,335 (34.6 %) | 299 | 224 |
| *Maveriviricetes* 321 | 7 (2.2 %) | 1 | 1 |
| *Polintoviricetes* 59 | 0 | 0 | 0 |
| under 10 kb 2,276 | 603 (26.5 %) | 112 | 79 |
| 10 to 30 kb 1,950 | 734 (37.6 %) | 188 | 146 |
| 30 to 100 kb 7 | 5 | 0 | 0 |

Percentages *derived* (meet / queries). Reading: about a third of the giant-virus fragments of `MVEUKV` lie, at the
vOTU threshold, inside a `GV1` contig (GVMAGs and MetaVR draw on the same IMG/M assemblies), though only 3.6 % are
byte-identical contigs (*derived*: 3,812 / 106,061). By the sample, about 31.7 % ± 1.4 % of `MVEUKV`'s single-record
genomes are such fragments (*derived*: two binomial standard errors at n = 4,233). They are staged: DEC-26 dedups by
exact `contig_uid` only, and the card records the sample, it does not cut by it.

## 7. The remainder (DEC-59)

MetaVR's remainder after `MV1`, `MV2` and `MVEUKV` is C-62's `MV1_left_for_later.tsv.gz` (23,419,058 rows) less
`MV2`'s 98,920 rows and `MVEUKV`'s 105,902 rows: 23,214,236 UViGs (*derived*). Its eukaryotic-virus part, one row
per vOTU (the vOTU's best eligible row, or its best row when none is eligible), is
`/tmp/if2453/c79/MVEUKV/MVEUKV_left_for_later.tsv.gz` (781,287 rows; *measured*):

| Reason | *Nucleocytoviricota* | *Maveriviricetes* | *Polintoviricetes* | All vOTUs | bp |
|---|---:|---:|---:|---:|---:|
| `low_confidence_ge5kb` (best eligible row at 5 kb and more, `viral_confidence` Low) | 464,147 | 19,043 | 3,312 | **486,502** | **4,496,427,184** |
| `lt5kb` (best eligible row under 5 kb) | 238,429 | 12,581 | 2,108 | **253,118** | **1,090,078,533** |
| `no_eligible_member` (only concatemers or contamination over 10 %) | 41,078 | 524 | 65 | 41,667 | 345,844,855 |

All *measured* (`stage`). A later stage `MVEUKV2` can take the `low_confidence_ge5kb` part at `DEPTH=called` (about
4.50 Gbp, about 420 to 610 CPU-minutes at class R2 of the plan's §3, *derived*: 4.496 × 93.5 and × 135.3). No
Mirusviricota can be taken by `ictv_taxonomy` (§2). The row is in `planning/CATALOG_GOAL.md` §4.1.

## 8. Build cost (*derived*, the plan's §3 rates; lower bounds, unpriced steps excluded)

- `MV2` (`called`, gv-meta): class R2, 93.5 to 135.3 CPU-minutes per new Gbp, × 3.762 to 3.861 new Gbp = about
  **350 to 520 CPU-minutes** (plan §4.7: 365 to 529).
- `MVEUKV` (`annotated`): R2 + R3, (93.5 + 2,090) to (135.3 + 2,740) per new Gbp, × 1.201 to 1.256 new Gbp = about
  **2,600 to 3,600 CPU-minutes** (plan §4.7: 1,700 to 4,900).

## 9. Pod compute used

All in detached `systemd-run --user` units at `nice -n 19`, `MemoryMax=32G`, `MemorySwapMax=0`, inside
`/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif` (the image card C-48c used), checked beforehand that no other
staging unit ran (only the `r26-11-4-wave-4` and `v2-4` build units were active) and that no `STOP` file existed:

| Unit | What | Cap | Wall (*measured*, file times) | CPU |
|---|---|---|---|---|
| `c79-profile` | schema, taxonomy values | `CPUQuota=200%`, DuckDB 2 threads | under 3 min (output 08:46:22) | at most 6 CPU-minutes (*derived*: wall × 2) |
| `c79-slice` | the slice count of §3 | same | under 2 min (output 08:48:42) | at most 4 |
| `c79-stage` | selection, one stream of `METAVR.fna.bgz`, both stages, ownership, `s2_check_inputs.py` | `CPUQuota=400%` (python, `bgzip -dc`, DuckDB 2 threads) | 20.7 min (12:30:30 to 12:51:13) | at most 83 CPU-minutes (*derived*: 20.7 × 4); the unit was collected before its `CPUUsageNSec` was read (150 s at 12:32, *measured*) |
| `c79-ani` | the GV1 sample | `CPUQuota=400%`, skani 3 threads | 1.8 min (12:52:52 to 12:54:41) | skani 147 CPU-seconds (*measured*, children's rusage) |

No process came near 300 CPU-minutes. The card's estimate of 2 to 3 h of wall at 4 threads (*estimate*) was not
needed: one pass took 20.7 minutes.

## 10. Tests

`tests/r1/test_c79_metavr_mv2_eukv.py` (an independent brute-force oracle for the selection over shuffled rows;
contig_uid by blake3; the stage end to end on a fixture table and FASTA, twice, byte-identical; refusals; the ANI
summary on a hand-written skani table; a skani run that skips without skani; the templates through
`stage_est.read_conf` with their `SUM`s found in this note). Seen failing first on the pod (7 failed in 3.13 s: no
`select_lists`, `Refused`, `stage`, `ani_summary`, `ani_sample`, and no templates), then passing (*measured*,
wilkcomp02, 2026-11-26, `scripts/pod_pytest.sh C-79 tests/r1/test_c79_metavr_mv2_eukv.py tests/r0`): **231 passed in
15.02 s**; with the pod's tools hidden (`PATH=/usr/bin:/bin`, the same image and selection): **230 passed, 1 skipped**
(the skani test). `stage_est.py conf` on `MV2.conf` and `MVEUKV.conf` (`python -I`, `--tree`): exit 0 for both.

## 11. Mirror to `tower`

`scripts/orchestrator/relay_to_tower.sh` (card C-49), from the laptop, one stream, four runs, to
`tower:/mnt/user/research-bulk/bio-box/sources/MV2/` and `.../MVEUKV/` (*measured*, 2026-11-26, the script's log):
`mv2.fna` 3,915,276,060 B in 158 s (24.78 MB/s), `mveukv.fna` 1,263,638,675 B in 45 s (28.08 MB/s), then the eight
metadata files of `MV2` and the nine of `MVEUKV` (with `MVEUKV_left_for_later.tsv.gz`); every file's sha256 on `tower`
equal to §5's, `DONE: 1, 1, 8, 9 verified, 0 failed`; `SHA256SUMS` written in each directory (9 and 10 lines). Both
stay `redistributable: false`: the mirror is the PI's own NAS, not a redistribution.
