# Which long-tail Rfam models hit on a sample of the corpus, 2026-11-12 (card 0-08a)

Card 0-08a, M1 of `planning/design/0-08_ncrna_models.md` (DEC-28's follow-up): full Rfam 15.1 on whole contigs of
a 100-genome prokaryote sample of release `rN.3` (stages P1 and P2, 2,497 genomes, stratified by GTDB R226
phylum) plus 300 phages of stage S1, with the producer's `cmsearch` settings, to decide which long-tail models join
the frozen set. Every number is *measured* on wilkcomp02 in the image `bio-box_80d4c55.sif` (Infernal 1.1.5,
Rfam 15.1 `Rfam.cm.gz`, sha256 `e2e636e9…8168`, equal to `refdb.yaml`; checked on the pod, `lib/Rfam.cm.gz.sha256`)
by the commands below, unless labelled *derived* (formula given) or *estimate*. The scan started 2026-11-11
21:58 CST (2026-11-12 UTC, the run records' dates).

**Status: complete.** All 27 processes finished on 2026-11-12 at 19:09 CST; the results section below was written by
`rfam_sample_report.py` from the scan's outputs, and the reading at the end by the orchestrating session.

**Measured against the estimate (AGENTS §6).** The note's M1 cost is 81 CPU-hours (*derived* there from 1,903
CPU-hours per 10 Gbp × an *estimated* 425 Mbp). The sample is 412,032,974 bases (*measured*), and the first 4
processes ran at 1,702 to 2,108 bases per CPU-second (*measured*, `runs/s00.json` to `s03.json`), against 1,454 in
0-09c: projected over the whole sample about 60 CPU-hours, about 25% under the estimate (*derived*: 9.21 CPU-hours
/ 62,675,593 bases × 412,032,974). The deviation is in the cheap direction and per-genome rates differ by 24% among
these four processes; the final figure is in the summary table.

<!-- results:begin -->

### Summary

| quantity | value | label |
|---|---|---|
| genomes scanned | 100 prokaryotes, 300 phages | sample rule |
| bases scanned | 412,032,974 | *measured*, wilkcomp02.ccbb.utexas.edu, 2026-11-12 |
| CPU (cmsearch, all processes) | 70.17 CPU-hours | *measured*, wilkcomp02.ccbb.utexas.edu, 2026-11-12, timed.py |
| bases per CPU-second | 1,631 | *derived*: bases / CPU-s |
| full library, CPU-hours per 10 Gbp | 1,703 | *derived*: CPU-s / bases x 1e10 / 3600 (0-09c: 1910) |
| whole sample's CPU against the note's 81 CPU-hour *estimate* (M1) | 70.2 CPU-hours (-13%) | *derived*: CPU-s / bases scanned x 412,032,974 sample bases / 3600 |
| largest process | 222.6 CPU-minutes | *measured*, wilkcomp02.ccbb.utexas.edu, 2026-11-12 |
| hits at or above GA | 11,912 | *measured*, wilkcomp02.ccbb.utexas.edu, 2026-11-12 |
| families with a hit (all / prokaryotes / phages) | 358 / 341 / 29 | *measured*, wilkcomp02.ccbb.utexas.edu, 2026-11-12 |
| ... of them in prok-1 / (a-) / excluded (c) / long tail (d) | 167 / 4 / 5 / 182 | *measured*, wilkcomp02.ccbb.utexas.edu, 2026-11-12 |
| families at 8 and at 100 prokaryotes (mean of 200 orders) | 82.8 and 341.0 | *derived* from measured hits |
| new families from the last prokaryote added (mean) | 1.92 | *derived* |
| long-tail (d) models that hit: candidates to add | 182 | *measured*, wilkcomp02.ccbb.utexas.edu, 2026-11-12 |
| excluded (c) models that hit | 5 models, 8 hits | *measured*, wilkcomp02.ccbb.utexas.edu, 2026-11-12 |
| rule-(b) models of prok-1 with no hit | 256 | *measured*, wilkcomp02.ccbb.utexas.edu, 2026-11-12 |
| hits mostly inside genes / outside every margin-200 region | 664 / 412 of 11912 | *measured*, wilkcomp02.ccbb.utexas.edu, 2026-11-12 |
| proposed prok-2: models (whole genome) | 770 (34) | rule in item 6 |
| prok-2 cost, CPU-hours per 10 Gbp | 299.7 (prok-1: 254.5) | *derived*: 254.5 x 291.4 / 247.4 |

### The sample

Seed 20261112; 2,497 prokaryote genomes in P1 and P2, GTDB R226 match exact 2,497; 55 phyla present (8 archaeal); phages every 90th of S1 by genome_uid. 676 contigs (0 circular); bases: prokaryotes 393,171,996, phages 18,860,978, all 412,032,974 (*measured*: contig_index lengths).

| stratum (GTDB domain:phylum) | genomes available | sampled | contigs | bases |
|---|---|---|---|---|
| Archaea:Asgardarchaeota | 2 | 1 | 1 | 1,418,269 |
| Archaea:Halobacteriota | 65 | 2 | 5 | 6,273,922 |
| Archaea:Iainarchaeota | 1 | 1 | 1 | 1,241,428 |
| Archaea:Methanobacteriota | 7 | 1 | 1 | 1,746,697 |
| Archaea:Methanobacteriota_B | 7 | 1 | 1 | 1,785,673 |
| Archaea:Nanobdellota | 1 | 1 | 1 | 490,885 |
| Archaea:Thermoplasmatota | 5 | 1 | 1 | 1,878,916 |
| Archaea:Thermoproteota | 25 | 2 | 5 | 4,723,118 |
| Bacteria:4484-113 | 1 | 1 | 1 | 4,112,261 |
| Bacteria:Acidobacteriota | 21 | 2 | 11 | 10,505,769 |
| Bacteria:Actinomycetota | 432 | 9 | 18 | 58,608,140 |
| Bacteria:Aquificota | 3 | 1 | 1 | 1,552,775 |
| Bacteria:Armatimonadota | 2 | 1 | 2 | 3,988,786 |
| Bacteria:Babelota | 1 | 1 | 1 | 1,123,447 |
| Bacteria:Bacillota | 372 | 8 | 10 | 28,314,618 |
| Bacteria:Bacteroidota | 233 | 5 | 7 | 22,295,072 |
| Bacteria:Bacteroidota_A | 13 | 1 | 1 | 2,763,181 |
| Bacteria:Bdellovibrionota | 4 | 1 | 1 | 3,601,631 |
| Bacteria:Bdellovibrionota_B | 3 | 1 | 1 | 2,806,830 |
| Bacteria:CAIWAD01 | 1 | 1 | 1 | 3,721,866 |
| Bacteria:Caldisericota | 1 | 1 | 1 | 1,558,103 |
| Bacteria:Campylobacterota | 43 | 2 | 3 | 4,152,562 |
| Bacteria:Chlamydiota | 6 | 1 | 1 | 1,101,914 |
| Bacteria:Chloroflexota | 29 | 2 | 91 | 5,363,840 |
| Bacteria:Chrysiogenota | 3 | 1 | 1 | 2,526,590 |
| Bacteria:Cloacimonadota | 7 | 1 | 2 | 3,814,504 |
| Bacteria:Coprothermobacterota | 1 | 1 | 1 | 1,424,912 |
| Bacteria:Cyanobacteriota | 33 | 2 | 2 | 9,492,689 |
| Bacteria:Deinococcota | 11 | 1 | 6 | 4,785,116 |
| Bacteria:Desulfobacterota | 43 | 2 | 3 | 8,717,329 |
| Bacteria:Desulfobacterota_I | 20 | 1 | 1 | 2,172,708 |
| Bacteria:Electryoneota | 1 | 1 | 35 | 5,165,167 |
| Bacteria:Elusimicrobiota | 3 | 1 | 1 | 1,588,979 |
| Bacteria:Fibrobacterota | 1 | 1 | 1 | 5,836,944 |
| Bacteria:Fusobacteriota | 15 | 1 | 2 | 3,346,458 |
| Bacteria:Gemmatimonadota | 4 | 1 | 1 | 3,977,883 |
| Bacteria:Hydrogenedentota | 1 | 1 | 24 | 6,145,620 |
| Bacteria:JAAXHH01 | 1 | 1 | 1 | 3,813,436 |
| Bacteria:Latescibacterota | 2 | 1 | 1 | 5,944,037 |
| Bacteria:Lernaellota | 1 | 1 | 9 | 6,306,852 |
| Bacteria:Marinisomatota | 1 | 1 | 1 | 2,335,625 |
| Bacteria:Myxococcota | 9 | 1 | 1 | 12,320,079 |
| Bacteria:Nitrospinota | 1 | 1 | 1 | 3,921,641 |
| Bacteria:Nitrospirota | 6 | 1 | 1 | 4,156,062 |
| Bacteria:Omnitrophota | 2 | 1 | 1 | 1,974,201 |
| Bacteria:Patescibacteriota | 26 | 2 | 2 | 2,263,726 |
| Bacteria:Planctomycetota | 21 | 1 | 1 | 6,690,754 |
| Bacteria:Poribacteria | 1 | 1 | 1 | 5,062,467 |
| Bacteria:Pseudomonadota | 954 | 18 | 36 | 85,102,639 |
| Bacteria:Spirochaetota | 31 | 2 | 67 | 5,583,178 |
| Bacteria:Sumerlaeota | 1 | 1 | 1 | 3,289,105 |
| Bacteria:Synergistota | 1 | 1 | 1 | 2,203,405 |
| Bacteria:Thermotogota | 6 | 1 | 1 | 1,915,344 |
| Bacteria:UBA8248 | 1 | 1 | 2 | 2,078,579 |
| Bacteria:Verrucomicrobiota | 11 | 1 | 1 | 4,086,264 |
| Viruses:phage (S1) | 27,258 | 300 | 300 | 18,860,978 |

### The runs

| process | genomes | records | bases | CPU-min | wall s | bases per CPU-s | peak RSS MiB | hits | families |
|---|---|---|---|---|---|---|---|---|---|
| s00 | 7 | 10 | 16,087,072 | 157.58 | 4,787 | 1,702 | 4,110 | 380 | 29 |
| s01 | 5 | 5 | 17,138,714 | 144.30 | 4,366 | 1,980 | 4,331 | 339 | 29 |
| s02 | 4 | 5 | 14,918,440 | 136.06 | 4,134 | 1,827 | 3,433 | 357 | 45 |
| s03 | 4 | 18 | 14,531,367 | 114.89 | 3,493 | 2,108 | 2,756 | 298 | 30 |
| s04 | 5 | 91 | 16,830,868 | 150.77 | 4,574 | 1,861 | 3,712 | 328 | 28 |
| s05 | 4 | 30 | 15,162,008 | 140.81 | 4,273 | 1,795 | 3,884 | 335 | 33 |
| s06 | 5 | 112 | 17,656,806 | 178.66 | 5,435 | 1,647 | 4,036 | 389 | 30 |
| s07 | 5 | 5 | 16,815,989 | 167.83 | 5,076 | 1,670 | 4,922 | 537 | 51 |
| s08 | 8 | 10 | 17,141,954 | 222.63 | 6,727 | 1,283 | 5,097 | 736 | 47 |
| s09 | 5 | 13 | 15,659,580 | 159.66 | 4,836 | 1,635 | 5,386 | 428 | 73 |
| s10 | 4 | 12 | 14,348,258 | 187.17 | 5,669 | 1,278 | 4,504 | 500 | 55 |
| s11 | 2 | 2 | 16,125,291 | 155.59 | 4,724 | 1,727 | 4,882 | 346 | 37 |
| s12 | 3 | 3 | 14,829,489 | 180.87 | 5,490 | 1,366 | 5,103 | 630 | 102 |
| s13 | 3 | 4 | 17,559,685 | 166.42 | 5,032 | 1,759 | 3,432 | 569 | 90 |
| s14 | 4 | 4 | 14,269,656 | 141.02 | 4,298 | 1,687 | 4,932 | 585 | 62 |
| s15 | 4 | 4 | 17,043,699 | 170.75 | 5,162 | 1,664 | 3,564 | 484 | 73 |
| s16 | 4 | 6 | 15,455,577 | 182.67 | 5,562 | 1,410 | 4,273 | 702 | 80 |
| s17 | 3 | 4 | 14,875,531 | 155.87 | 4,742 | 1,591 | 3,921 | 547 | 83 |
| s18 | 3 | 5 | 17,265,429 | 200.33 | 6,074 | 1,436 | 4,768 | 569 | 89 |
| s19 | 4 | 5 | 13,562,706 | 119.95 | 3,645 | 1,885 | 2,932 | 269 | 32 |
| s20 | 1 | 2 | 8,981,933 | 83.98 | 2,558 | 1,783 | 2,606 | 206 | 45 |
| s21 | 1 | 1 | 10,812,875 | 100.71 | 3,079 | 1,789 | 2,598 | 169 | 47 |
| s22 | 3 | 8 | 15,908,517 | 165.09 | 5,026 | 1,606 | 4,575 | 336 | 53 |
| s23 | 3 | 9 | 13,578,089 | 163.98 | 4,962 | 1,380 | 4,452 | 496 | 57 |
| s24 | 3 | 3 | 11,569,552 | 120.76 | 3,680 | 1,597 | 3,771 | 363 | 42 |
| s25 | 56 | 58 | 17,990,498 | 177.61 | 5,402 | 1,688 | 4,546 | 346 | 36 |
| s26 | 247 | 247 | 15,913,391 | 164.50 | 4,985 | 1,612 | 886 | 668 | 25 |

All *measured*, wilkcomp02.ccbb.utexas.edu, 2026-11-12, `timed.py` records `runs/s<k>.json`; bases per CPU-s *derived* (bases / CPU-s).

### 1. Families with a hit as genomes are added

Prokaryotes (*measured*, wilkcomp02.ccbb.utexas.edu, 2026-11-12; *derived*: the number of distinct models with a hit in the first n genomes, over 200 orders seeded `0-08a:curve:<i>`):

| genomes | families (mean) | min | max |
|---|---|---|---|
| 8 | 82.8 | 31 | 166 |
| 10 | 95.4 | 43 | 175 |
| 20 | 148.9 | 70 | 226 |
| 30 | 188.8 | 94 | 251 |
| 40 | 220.7 | 116 | 278 |
| 50 | 248.2 | 189 | 306 |
| 60 | 270.3 | 229 | 313 |
| 70 | 290.6 | 251 | 323 |
| 80 | 309.3 | 271 | 336 |
| 90 | 325.5 | 295 | 340 |
| 100 | 341.0 | 341 | 341 |

The 8 golden genomes of 0-09c had 258 families; 93 of them hit here. By class, prokaryotes: prok-1 162, a- 4, c 5, d 170.

Phages (same method):

| genomes | families (mean) | min | max |
|---|---|---|---|
| 25 | 6.1 | 2 | 17 |
| 50 | 10.2 | 2 | 21 |
| 100 | 16.4 | 5 | 26 |
| 150 | 20.4 | 10 | 28 |
| 200 | 24.0 | 17 | 29 |
| 250 | 26.8 | 21 | 29 |
| 300 | 29.0 | 29 | 29 |

By class, phages: prok-1 9, a- 0, c 0, d 20. Families hit in phages and in no prokaryote of the sample: AsxR, GOLLD, Hatchet, IS009, Lambda_thermo, MicC, PhotoRC-II, RAGATH-20, RAGATH-21, RAGATH-7, RAGATH-9, STnc100, Spy779816, int-alpA, sar, skipping-rope, sok.

By GTDB phylum (prokaryotes):

| phylum | genomes | families | only in this phylum | long-tail (d) families |
|---|---|---|---|---|
| Archaea:Asgardarchaeota | 1 | 12 | - | - |
| Archaea:Halobacteriota | 2 | 23 | CRISPR-DR24, L4-Archaeoglobi, S15-Halobacteria, S15-Methanomicrobia, leuA-Halobacteria, sR-tMet, sR3, sot0042 | L4-Archaeoglobi, S15-Halobacteria, S15-Methanomicrobia, eL15-Euryarchaeota, leuA-Halobacteria, sot0042 |
| Archaea:Iainarchaeota | 1 | 12 | - | - |
| Archaea:Methanobacteriota | 1 | 18 | CRISPR-DR62, flpD, sR2 | flpD |
| Archaea:Methanobacteriota_B | 1 | 37 | CRISPR-DR6, HgcE, HgcF, HgcG, Pab160, Pab19, Pab91, SscA, sR22, sR23, sR28, sR34, sR36, sR43, sR46, sR51, sR55, sR58, sR60, snoPyro_CD, snoR9 | HgcE, HgcF, HgcG, Pab160, Pab19, Pab91, SscA, eL15-Euryarchaeota |
| Archaea:Nanobdellota | 1 | 11 | - | - |
| Archaea:Thermoplasmatota | 1 | 14 | - | - |
| Archaea:Thermoproteota | 2 | 16 | sR14, sR20, sR5 | - |
| Bacteria:4484-113 | 1 | 16 | - | - |
| Bacteria:Acidobacteriota | 2 | 25 | - | - |
| Bacteria:Actinomycetota | 9 | 59 | 6C, ASdes, ASpks, Actino-pnp, F6, L17-Actino-Proteobacteria, L31-Actinobacteria, L31-Corynebacteriaceae, Ms_AS-1, Ms_IGR-5, SAM-IV, Scr1601, Scr2736, Scr3202, Scr3920, Scr4115, Scr5676, Scr6106, Streptomyces-metH, Streptomyces-metK, atpB, che1, dauA, gntR-DTE, mraW, msiK, nrdJ, raiA-hairpin, scr4677, scr5239, sodF_sRNA | ALIL, L17-Actino-Proteobacteria, L31-Actinobacteria, L31-Corynebacteriaceae, Ms_IGR-5, Scr1601, Scr2736, Scr3202, Scr3920, Scr4115, Scr5676, Scr6106, Sernc350, Streptomyces-metH, Streptomyces-metK, atpB, che1, dauA, gntR-DTE, msiK, nrdJ, raiA-hairpin, scr4677, scr5239, sodF_sRNA |
| Bacteria:Aquificota | 1 | 14 | - | - |
| Bacteria:Armatimonadota | 1 | 18 | - | - |
| Bacteria:Babelota | 1 | 13 | - | - |
| Bacteria:Bacillota | 8 | 79 | Bacillaceae-1, BsrC, BsrF, BsrG, L13_leader, L19_leader, L20_leader, L21_leader, L31-Firmicutes, OLE, PreQ1, RAGATH-18, RNaseP_bact_b, Ref64, RsaE, S10-Clostridia, S4-Clostridia, SpF59_sRNA, THF, antisense_SR1, epsC, glycine-GGAnGA, group-II-D1D4-5, gut-2, ncr1575, pan, pemK, raiA, rli30, sRNA_0030, sRNA_1300, uup, ylbH | OLE, RAGATH-18, Ref64, SECIS_4, STnc490, gut-2, raiA, sRNA_0030, sRNA_1300, uup |
| Bacteria:Bacteroidota | 5 | 32 | 6S-Flavo, Bacteroid-trp, Bacteroidales_small_SRP, Flavo-1, L19-Flavobacteria, RT-17, S15-Flavobacteria, S16-Flavobacteria, S4-Flavobacteria, SAM-II, whalefall-1 | Bacteroid-trp, Flavo-1, L19-Flavobacteria, RT-17, S15-Flavobacteria, S16-Flavobacteria, S4-Flavobacteria, whalefall-1 |
| Bacteria:Bacteroidota_A | 1 | 26 | CRISPR-DR38, Chlorobi-1, Chlorobi-RRM, Lnt, S6-S18-Chlorobi, SAM-Chlorobi | Chlorobi-1, Chlorobi-RRM, Lnt, S6-S18-Chlorobi, SAM-Chlorobi |
| Bacteria:Bdellovibrionota | 1 | 17 | - | - |
| Bacteria:Bdellovibrionota_B | 1 | 15 | - | HEARO |
| Bacteria:CAIWAD01 | 1 | 19 | - | - |
| Bacteria:Caldisericota | 1 | 19 | - | - |
| Bacteria:Campylobacterota | 2 | 20 | HPnc0260, purD | purD |
| Bacteria:Chlamydiota | 1 | 14 | IhtA, Plasmid_R1162 | IhtA, Plasmid_R1162 |
| Bacteria:Chloroflexota | 2 | 22 | - | - |
| Bacteria:Chrysiogenota | 1 | 20 | - | - |
| Bacteria:Cloacimonadota | 1 | 19 | - | - |
| Bacteria:Coprothermobacterota | 1 | 23 | 2dG-II | SECIS_4 |
| Bacteria:Cyanobacteriota | 2 | 34 | CRISPR-DR28, CRISPR-DR5, CyVA-1, Cyano-1, Glutamine, IsrR, NsiR1, NsiR4, Plastid-ndhA-ndhB-2, Yfr1 | CyVA-1, Cyano-1, HEARO, IsrR, NsiR1, NsiR4, Yfr1 |
| Bacteria:Deinococcota | 1 | 21 | c-di-GMP-II, c-di-GMP-II-GAG, c-di-GMP-II-GCG | - |
| Bacteria:Desulfobacterota | 2 | 27 | c4-2 | PYLIS_3, c4-2 |
| Bacteria:Desulfobacterota_I | 1 | 16 | - | - |
| Bacteria:Electryoneota | 1 | 24 | DUF3800-VIII | - |
| Bacteria:Elusimicrobiota | 1 | 13 | - | - |
| Bacteria:Fibrobacterota | 1 | 22 | - | ALIL |
| Bacteria:Fusobacteriota | 1 | 28 | NiCo, S4-Fusobacteriales | S4-Fusobacteriales |
| Bacteria:Gemmatimonadota | 1 | 21 | - | - |
| Bacteria:Hydrogenedentota | 1 | 21 | - | nhaA-I |
| Bacteria:JAAXHH01 | 1 | 18 | - | - |
| Bacteria:Latescibacterota | 1 | 22 | - | - |
| Bacteria:Lernaellota | 1 | 23 | Acido-Lenti-1 | Acido-Lenti-1 |
| Bacteria:Marinisomatota | 1 | 17 | - | - |
| Bacteria:Myxococcota | 1 | 17 | - | Sernc350, nhaA-I |
| Bacteria:Nitrospinota | 1 | 21 | - | - |
| Bacteria:Nitrospirota | 1 | 17 | - | - |
| Bacteria:Omnitrophota | 1 | 14 | - | - |
| Bacteria:Patescibacteriota | 2 | 13 | - | - |
| Bacteria:Planctomycetota | 1 | 18 | - | - |
| Bacteria:Poribacteria | 1 | 18 | - | - |
| Bacteria:Pseudomonadota | 18 | 177 | AAC_AAD_leader, Aar, AbsR11, Alpha_RBS, AsdA, Atu_At1, Atu_C6, Atu_C8, BASRCI27, Betaproteobacteria_toxic_sRNA, BjrC1505, C4, CRISPR-DR4, CcsR1, Cis52_sRNA, CrcZ, DABA-DC-AT, GcvB, Hammerhead_II, His_leader, L2-Alphaproteobacteria, L31-Gammaproteobacteria, MH_s36, MH_s7, NrsZ, P1, P10, P11, P14, P15, P17, P18, P24, P26, P27, P31, P34, P36, P4, P9, PhrS, Polynucleobacter-1, Proteo-phage-1, PrrB_RsmZ, PrrF, Pseudomon-1, Pseudomon-Rho, Pseudomon-groES, Pseudomonadales-1, PtaRNA1, Qrr, ROSE, ROSE_3, RT-12, RgsA, Rhodo-rpoB, RpsF_leader, RsmY, RyhB, S15, SAH, SAM-SAH, SdsR_RyeB, SorX, Spot_42, SrbA, THF-II, Thr_leader, TwoAYGGAY, V_AS7, V_AS9, Xanthine, Xoo1, Xoo2, Xoo5, Xoo8, YrlA, abiF, alpha_tmRNA, ar14, ar15, ar45, ar7, ar9, asX1, asX4, asponA, c4-a1b1, carA, chrB-a, chrB-b, ffh, gabT, gltS, gyrA, icd-II, isrK, ivy-DE, narK, proV, radC, rimP, rmf, rne-II, rnk_pseudo, rpfG, rpsL_pseudo, rsmX, sRNA-Xcc1, sX11, sX12, sX13, sX14, sX15, sX2, sX4, sX5, sX6, sX7, sX9, serC, speF, sucA, sucA-II, sucC, suhB, sul1, t44, terC, tmRNA_permuted, traJ-II | ALIL, Aar, AbsR11, Atu_At1, Atu_C6, Atu_C8, BASRCI27, Betaproteobacteria_toxic_sRNA, BjrC1505, C4, CcsR1, Cis52_sRNA, DABA-DC-AT, Hammerhead_II, L2-Alphaproteobacteria, MH_s36, MH_s7, NrsZ, PYLIS_3, PhrS, Polynucleobacter-1, Proteo-phage-1, PrrB_RsmZ, PrrF, Pseudomon-1, Pseudomon-Rho, Pseudomon-groES, Pseudomonadales-1, PtaRNA1, Qrr, ROSE, ROSE_3, RT-12, RgsA, Rhodo-rpoB, RsmY, SECIS_4, STnc490, SorX, SrbA, TwoAYGGAY, Xoo1, Xoo2, Xoo5, Xoo8, YrlA, abiF, ar14, ar15, ar45, ar7, ar9, asX1, asX4, asponA, c4-a1b1, carA, chrB-a, chrB-b, ffh, gabT, gltS, gyrA, icd-II, isrK, ivy-DE, narK, nhaA-I, proV, radC, rmf, rne-II, rpfG, rsmX, sRNA-Xcc1, sX11, sX12, sX13, sX14, sX15, sX2, sX4, sX5, sX6, sX7, sX9, serC, speF, sucA, sucA-II, sucC, suhB, sul1, terC, traJ-II |
| Bacteria:Spirochaetota | 2 | 25 | RT-14, SNORA17, TD-1, TD-2 | RT-14, SECIS_4, TD-1, TD-2 |
| Bacteria:Sumerlaeota | 1 | 19 | - | - |
| Bacteria:Synergistota | 1 | 26 | - | ALIL, SECIS_4 |
| Bacteria:Thermotogota | 1 | 21 | - | - |
| Bacteria:UBA8248 | 1 | 18 | - | - |
| Bacteria:Verrucomicrobiota | 1 | 15 | - | - |

### 2. Long-tail (d) models that hit: candidates to add

Hits and genomes *measured*; CPU-s is the model's own `# CPU time` line summed over this scan's processes (*measured*), its share of all models' lines; CPU-h/10 Gbp *derived* from the model's 8-genome lines (design note formula) for the intergenic (margin 200, `--notrunc`) and whole-genome scans. New loci: hits that no prok-1 model's hit overlaps on the same strand.

| accession | name | description | hits | genomes (prok/phage) | new loci | max score / GA | CPU-s | share | h/10 Gbp igr / wg | phyla |
|---|---|---|---|---|---|---|---|---|---|---|
| RF00043 | Plasmid_R1162 | R1162-like plasmid antisense RNA | 1 | 1 (1/0) | 1 | 73.3 / 55.0 | 41.3 | 0.016% | 0.10 / 0.27 | Chlamydiota |
| RF00058 | HgcF | HgcF RNA (Pab35) | 1 | 1 (1/0) | 1 | 156.0 / 70.0 | 104.7 | 0.042% | 0.21 / 0.65 | Methanobacteriota_B |
| RF00060 | HgcE | HgcE RNA (Pab105) | 1 | 1 (1/0) | 1 | 87.4 / 50.0 | 140.1 | 0.056% | 0.21 / 0.54 | Methanobacteriota_B |
| RF00063 | SscA | SscA RNA | 1 | 1 (1/0) | 1 | 85.8 / 60.0 | 67.4 | 0.027% | 0.19 / 0.37 | Methanobacteriota_B |
| RF00064 | HgcG | HgcG RNA (Pab40) | 1 | 1 (1/0) | 1 | 188.5 / 77.0 | 162.6 | 0.065% | 0.49 / 1.08 | Methanobacteriota_B |
| RF00166 | PrrB_RsmZ | PrrB/RsmZ RNA family | 4 | 4 (4/0) | 4 | 133.1 / 44.0 | 75.4 | 0.030% | 0.22 / 0.50 | Pseudomonadota |
| RF00195 | RsmY | RsmY RNA family | 4 | 4 (4/0) | 4 | 125.1 / 57.0 | 56.4 | 0.022% | 0.16 / 0.38 | Pseudomonadota |
| RF00262 | sar | sar RNA | 1 | 1 (0/1) | 1 | 97.8 / 43.0 | 29.9 | 0.012% | 0.08 / 0.19 | phage |
| RF00378 | Qrr | Qrr RNA | 8 | 1 (1/0) | 8 | 100.9 / 50.0 | 41.9 | 0.017% | 0.15 / 0.32 | Pseudomonadota |
| RF00435 | ROSE | Repression of heat shock gene expression (ROSE) element | 2 | 2 (2/0) | 2 | 73.2 / 45.0 | 38.5 | 0.015% | 0.08 / 0.21 | Pseudomonadota |
| RF00444 | PrrF | PrrF RNA | 8 | 4 (4/0) | 8 | 162.7 / 90.0 | 50.9 | 0.020% | 0.26 / 0.46 | Pseudomonadota |
| RF00517 | serC | serC leader | 2 | 2 (2/0) | 2 | 69.3 / 35.0 | 21.7 | 0.009% | 0.06 / 0.14 | Pseudomonadota |
| RF00518 | speF | speF leader | 1 | 1 (1/0) | 1 | 94.1 / 37.0 | 57.3 | 0.023% | 0.17 / 0.49 | Pseudomonadota |
| RF00519 | suhB | Makes More Granules Regulator RNA (mmgR) | 5 | 2 (2/0) | 5 | 64.6 / 41.0 | 41.6 | 0.017% | 0.13 / 0.23 | Pseudomonadota |
| RF00628 | RgsA | RgsA sRNA | 4 | 4 (4/0) | 4 | 204.3 / 88.4 | 89.6 | 0.036% | 0.31 / 0.70 | Pseudomonadota |
| RF01069 | purD | purD RNA motif | 1 | 1 (1/0) | 1 | 115.2 / 41.0 | 36.2 | 0.014% | 0.10 / 0.25 | Campylobacterota |
| RF01070 | sucA | SucA RNA motif | 1 | 1 (1/0) | 1 | 64.4 / 40.0 | 39.8 | 0.016% | 0.10 / 0.30 | Pseudomonadota |
| RF01071 | OLE | Ornate Large Extremophilic RNA | 1 | 1 (1/0) | 1 | 544.5 / 50.0 | 143.1 | 0.057% | 0.24 / 1.15 | Bacillota |
| RF01116 | Yfr1 | Cyanobacterial functional RNA 1 | 2 | 2 (2/0) | 2 | 61.8 / 38.0 | 27.1 | 0.011% | 0.07 / 0.15 | Cyanobacteriota |
| RF01394 | isrK | isrK Hfq binding RNA | 9 | 6 (3/3) | 9 | 74.7 / 50.0 | 53.7 | 0.021% | 0.20 / 0.40 | Pseudomonadota, phage |
| RF01405 | STnc490 | STnc490 Hfq binding RNA | 3 | 2 (2/0) | 3 | 113.2 / 56.0 | 48.1 | 0.019% | 0.15 / 0.42 | Bacillota, Pseudomonadota |
| RF01419 | IsrR | Antisense RNA which regulates isiA expression | 2 | 1 (1/0) | 2 | 160.7 / 52.0 | 102.8 | 0.041% | 0.35 / 0.89 | Cyanobacteriota |
| RF01497 | ALIL | ALIL pseudoknot | 16 | 9 (6/3) | 16 | 121.5 / 50.0 | 100.9 | 0.040% | 0.26 / 0.57 | Actinomycetota, Fibrobacterota, Pseudomonadota, Synergistota, phage |
| RF01673 | PhrS | PhrS | 5 | 5 (5/0) | 5 | 219.3 / 37.0 | 245.9 | 0.098% | 0.95 / 2.34 | Pseudomonadota |
| RF01687 | Acido-Lenti-1 | Acido-Lenti-1 RNA | 5 | 1 (1/0) | 5 | 43.2 / 37.0 | 18.9 | 0.008% | 0.06 / 0.12 | Lernaellota |
| RF01692 | Bacteroid-trp | Bacteroidete tryptophan peptide leader RNA | 1 | 1 (1/0) | 1 | 56.0 / 49.0 | 65.2 | 0.026% | 0.30 / 0.57 | Bacteroidota |
| RF01695 | C4 | C4 antisense RNA | 104 | 12 (8/4) | 104 | 86.0 / 45.0 | 68.4 | 0.027% | 0.13 / 0.42 | Pseudomonadota, phage |
| RF01696 | Chlorobi-1 | Chlorobi-1 RNA | 3 | 1 (1/0) | 3 | 88.3 / 50.0 | 31.9 | 0.013% | 0.08 / 0.18 | Bacteroidota_A |
| RF01697 | Chlorobi-RRM | Chlorobi-RRM RNA | 2 | 1 (1/0) | 2 | 79.1 / 45.0 | 33.4 | 0.013% | 0.10 / 0.27 | Bacteroidota_A |
| RF01701 | Cyano-1 | Cyano-1 RNA | 2 | 1 (1/0) | 2 | 52.9 / 41.0 | 23.9 | 0.009% | 0.08 / 0.17 | Cyanobacteriota |
| RF01705 | Flavo-1 | Flavo-1 RNA | 12 | 1 (1/0) | 12 | 65.6 / 42.0 | 62.9 | 0.025% | 0.29 / 0.53 | Bacteroidota |
| RF01711 | Lnt | Lnt RNA | 1 | 1 (1/0) | 1 | 76.2 / 55.0 | 39.4 | 0.016% | 0.14 / 0.30 | Bacteroidota_A |
| RF01717 | PhotoRC-II | PhotoRC-II RNA | 3 | 3 (0/3) | 3 | 75.4 / 39.0 | 29.7 | 0.012% | 0.08 / 0.24 | phage |
| RF01718 | Polynucleobacter-1 | Polynucleobacter-1 RNA | 1 | 1 (1/0) | 1 | 95.5 / 37.0 | 31.6 | 0.013% | 0.08 / 0.21 | Pseudomonadota |
| RF01719 | Pseudomon-1 | Pseudomon-1/ErsA RNA | 5 | 5 (5/0) | 5 | 144.3 / 50.0 | 35.5 | 0.014% | 0.11 / 0.21 | Pseudomonadota |
| RF01720 | Pseudomon-Rho | Pseudomon-Rho RNA | 4 | 4 (4/0) | 4 | 112.8 / 45.0 | 152.0 | 0.060% | 0.21 / 0.58 | Pseudomonadota |
| RF01721 | Pseudomon-groES | Pseudomon-groES RNA | 4 | 4 (4/0) | 4 | 98.8 / 50.0 | 208.5 | 0.083% | 0.47 / 1.92 | Pseudomonadota |
| RF01724 | SAM-Chlorobi | SAM-Chlorobi RNA | 1 | 1 (1/0) | 1 | 96.2 / 35.0 | 32.3 | 0.013% | 0.08 / 0.20 | Bacteroidota_A |
| RF01731 | TwoAYGGAY | TwoAYGGAY RNA | 53 | 4 (4/0) | 53 | 62.3 / 43.0 | 86.3 | 0.034% | 0.18 / 0.35 | Pseudomonadota |
| RF01737 | flpD | flpD RNA | 1 | 1 (1/0) | 1 | 68.3 / 40.0 | 25.5 | 0.010% | 0.06 / 0.15 | Methanobacteriota |
| RF01738 | gabT | gabT RNA | 3 | 3 (3/0) | 3 | 71.2 / 45.0 | 38.5 | 0.015% | 0.07 / 0.20 | Pseudomonadota |
| RF01740 | gyrA | gyrA RNA | 6 | 4 (4/0) | 6 | 76.4 / 45.0 | 80.1 | 0.032% | 0.12 / 0.35 | Pseudomonadota |
| RF01747 | msiK | msiK RNA | 8 | 8 (8/0) | 8 | 64.6 / 33.0 | 31.2 | 0.012% | 0.07 / 0.18 | Actinomycetota |
| RF01754 | radC | radC RNA | 5 | 2 (2/0) | 5 | 89.5 / 40.0 | 338.6 | 0.134% | 0.40 / 1.40 | Pseudomonadota |
| RF01755 | rmf | rmf RNA | 4 | 4 (4/0) | 4 | 124.1 / 45.0 | 117.8 | 0.047% | 0.23 / 0.84 | Pseudomonadota |
| RF01756 | rne-II | rne-II RNA | 4 | 4 (4/0) | 4 | 267.2 / 50.0 | 160.7 | 0.064% | 0.36 / 0.86 | Pseudomonadota |
| RF01758 | sucA-II | sucA-II RNA | 4 | 4 (4/0) | 4 | 189.1 / 70.0 | 81.0 | 0.032% | 0.14 / 0.50 | Pseudomonadota |
| RF01759 | sucC | sucC RNA | 4 | 4 (4/0) | 4 | 85.9 / 50.0 | 41.2 | 0.016% | 0.10 / 0.26 | Pseudomonadota |
| RF01760 | traJ-II | traJ-II RNA | 2 | 2 (2/0) | 2 | 112.6 / 40.0 | 32.4 | 0.013% | 0.09 / 0.19 | Pseudomonadota |
| RF01762 | whalefall-1 | Whalefall-1 RNA | 1 | 1 (1/0) | 1 | 37.8 / 37.0 | 22.6 | 0.009% | 0.07 / 0.15 | Bacteroidota |
| RF01793 | ffh | ffh sRNA | 3 | 3 (3/0) | 3 | 62.0 / 35.0 | 38.0 | 0.015% | 0.08 / 0.21 | Pseudomonadota |
| RF01811 | PtaRNA1 | Plasmid transferred anti-sense RNA | 1 | 1 (1/0) | 1 | 78.2 / 47.0 | 58.1 | 0.023% | 0.11 / 0.28 | Pseudomonadota |
| RF01990 | SECIS_4 | Selenocysteine insertion sequence 4 | 10 | 5 (5/0) | 10 | 64.6 / 38.0 | 20.4 | 0.008% | 0.06 / 0.14 | Bacillota, Coprothermobacterota, Pseudomonadota, Spirochaetota, Synergistota |
| RF02032 | GOLLD | Giant, ornate, lake- and Lactobacillales-derived (GOLLD) RNA | 5 | 5 (0/5) | 5 | 176.4 / 150.0 | 197.3 | 0.078% | 0.44 / 1.92 | phage |
| RF02033 | HEARO | HNH endonuclease-associated RNA and ORF (HEARO) RNA | 6 | 2 (2/0) | 6 | 113.4 / 83.0 | 37.5 | 0.015% | 0.08 / 0.23 | Bdellovibrionota_B, Cyanobacteriota |
| RF02076 | STnc100 | Gammaproteobacterial sRNA STnc100 | 3 | 3 (0/3) | 3 | 206.9 / 47.0 | 77.4 | 0.031% | 0.21 / 0.57 | phage |
| RF02144 | rsmX | rsmX | 4 | 2 (2/0) | 4 | 101.8 / 56.0 | 104.2 | 0.041% | 0.36 / 0.65 | Pseudomonadota |
| RF02221 | sRNA-Xcc1 | sRNA-Xcc1 | 15 | 4 (2/2) | 15 | 81.5 / 53.0 | 30.9 | 0.012% | 0.08 / 0.19 | Pseudomonadota, phage |
| RF02222 | sX2 | Xanthomonas sRNA sX2 | 2 | 2 (2/0) | 2 | 116.4 / 40.0 | 43.4 | 0.017% | 0.11 / 0.29 | Pseudomonadota |
| RF02223 | sX4 | Proteobacterial sRNA sX4 | 2 | 2 (2/0) | 0 | 195.7 / 70.0 | 110.5 | 0.044% | 0.36 / 0.66 | Pseudomonadota |
| RF02224 | sX5 | Xanthomonadaceae sRNA sX5 | 2 | 2 (2/0) | 2 | 82.5 / 50.0 | 25.0 | 0.010% | 0.06 / 0.14 | Pseudomonadota |
| RF02225 | sX6 | Proteobacterial sRNA sX6 | 2 | 2 (2/0) | 2 | 307.3 / 40.0 | 1979.2 | 0.785% | 0.90 / 5.63 | Pseudomonadota |
| RF02226 | sX7 | Xanthomonas sRNA sX7 | 2 | 2 (2/0) | 2 | 114.5 / 55.0 | 92.8 | 0.037% | 0.12 / 0.36 | Pseudomonadota |
| RF02228 | sX9 | Xanthomonas sRNA sX9 | 81 | 2 (2/0) | 81 | 95.1 / 44.0 | 89.4 | 0.035% | 0.11 / 0.39 | Pseudomonadota |
| RF02230 | sX11 | Proteobacterial sRNA sX11 | 2 | 2 (2/0) | 2 | 153.4 / 70.0 | 56.5 | 0.022% | 0.08 / 0.22 | Pseudomonadota |
| RF02231 | sX12 | Xanthomonas sRNA sX12 | 1 | 1 (1/0) | 1 | 107.5 / 106.0 | 78.5 | 0.031% | 0.09 / 0.27 | Pseudomonadota |
| RF02232 | sX13 | Xanthomonadaceae sRNA sX13 | 2 | 2 (2/0) | 2 | 142.2 / 50.0 | 113.6 | 0.045% | 0.10 / 0.29 | Pseudomonadota |
| RF02233 | sX14 | Xanthomonas sRNA sX14/Xoo3 | 2 | 2 (2/0) | 2 | 119.2 / 80.0 | 46.4 | 0.018% | 0.09 / 0.23 | Pseudomonadota |
| RF02234 | sX15 | Xanthomonadaceae sRNA sX15 | 15 | 2 (2/0) | 15 | 170.0 / 55.0 | 127.7 | 0.051% | 0.12 / 0.70 | Pseudomonadota |
| RF02235 | asX1 | Xanthomonadaceae sRNA asX1 | 2 | 2 (2/0) | 2 | 129.1 / 70.0 | 90.8 | 0.036% | 0.13 / 0.39 | Pseudomonadota |
| RF02238 | asX4 | Xanthomonas sRNA asX4/Xoo4 | 2 | 2 (2/0) | 2 | 314.9 / 100.0 | 293.6 | 0.117% | 0.13 / 1.18 | Pseudomonadota |
| RF02240 | Xoo1 | Xanthomonadaceae sRNA Xoo1 | 2 | 2 (2/0) | 2 | 108.7 / 70.0 | 62.3 | 0.025% | 0.12 / 0.38 | Pseudomonadota |
| RF02241 | Xoo2 | Xanthomonadaceae sRNA Xoo2 | 2 | 2 (2/0) | 2 | 106.7 / 45.0 | 102.4 | 0.041% | 0.15 / 0.53 | Pseudomonadota |
| RF02242 | Xoo5 | Xanthomonadaceae sRNA Xoo5 | 2 | 2 (2/0) | 2 | 143.9 / 50.0 | 428.0 | 0.170% | 0.48 / 2.11 | Pseudomonadota |
| RF02243 | Xoo8 | Proteobacterial sRNA Xoo8 | 2 | 2 (2/0) | 2 | 349.7 / 90.0 | 428.7 | 0.170% | 0.71 / 2.35 | Pseudomonadota |
| RF02276 | Hammerhead_II | Hammerhead ribozyme (type II) | 9 | 3 (1/2) | 9 | 78.0 / 37.0 | 61.8 | 0.025% | 0.14 / 0.34 | Pseudomonadota, phage |
| RF02278 | Betaproteobacteria_toxic_sRNA | Betaproteobacteria toxic small RNA | 2 | 2 (2/0) | 2 | 63.9 / 50.0 | 27.4 | 0.011% | 0.08 / 0.18 | Pseudomonadota |
| RF02342 | ar7 | Alphaproteobacterial sRNA ar7 | 1 | 1 (1/0) | 1 | 118.6 / 60.0 | 140.1 | 0.056% | 0.19 / 0.52 | Pseudomonadota |
| RF02343 | ar9 | Alphaproteobacterial sRNA ar9 | 1 | 1 (1/0) | 1 | 120.7 / 45.0 | 54.1 | 0.021% | 0.13 / 0.37 | Pseudomonadota |
| RF02344 | ar14 | Alphaproteobacterial ar14 | 4 | 1 (1/0) | 4 | 93.1 / 58.8 | 135.4 | 0.054% | 0.09 / 0.31 | Pseudomonadota |
| RF02345 | ar15 | Alphaproteobacterial ar15 | 2 | 1 (1/0) | 2 | 50.8 / 40.0 | 36.8 | 0.015% | 0.10 / 0.24 | Pseudomonadota |
| RF02347 | ar45 | Alphaproteobacterial sRNA ar45 | 1 | 1 (1/0) | 1 | 153.3 / 40.0 | 55.2 | 0.022% | 0.09 / 0.31 | Pseudomonadota |
| RF02356 | BjrC1505 | Alphaproteobacterial sRNA BjrC1505 | 1 | 1 (1/0) | 1 | 74.9 / 50.0 | 42.9 | 0.017% | 0.14 / 0.33 | Pseudomonadota |
| RF02375 | Aar | Aar sRNA | 1 | 1 (1/0) | 1 | 113.2 / 60.0 | 68.0 | 0.027% | 0.25 / 0.54 | Pseudomonadota |
| RF02399 | NsiR1 | Nitrogen stress-induced RNA 1 | 13 | 1 (1/0) | 13 | 58.2 / 39.0 | 27.3 | 0.011% | 0.11 / 0.20 | Cyanobacteriota |
| RF02471 | Ms_IGR-5 | Actinobacteria sRNA Ms_IGR-5 | 6 | 6 (6/0) | 6 | 122.2 / 45.0 | 84.3 | 0.033% | 0.19 / 0.48 | Actinomycetota |
| RF02496 | Atu_At1 | Rhizobiales sRNA Atu_At1 | 1 | 1 (1/0) | 1 | 100.1 / 60.0 | 83.3 | 0.033% | 0.18 / 0.53 | Pseudomonadota |
| RF02500 | Atu_C6 | EcpR1 | 1 | 1 (1/0) | 1 | 105.9 / 60.0 | 94.2 | 0.037% | 0.11 / 0.37 | Pseudomonadota |
| RF02502 | Atu_C8 | Rhizobiales sRNA Atu_C8 | 1 | 1 (1/0) | 1 | 94.5 / 50.0 | 77.0 | 0.031% | 0.19 / 0.44 | Pseudomonadota |
| RF02510 | PYLIS_3 | Pyrrolysine insertion sequence mttB | 2 | 2 (2/0) | 2 | 60.9 / 38.0 | 53.5 | 0.021% | 0.16 / 0.50 | Desulfobacterota, Pseudomonadota |
| RF02523 | ROSE_3 | Repression of heat shock gene expression (ROSE) element | 3 | 3 (3/0) | 3 | 78.7 / 43.0 | 46.2 | 0.018% | 0.10 / 0.27 | Pseudomonadota |
| RF02553 | YrlA | Y RNA-like | 3 | 2 (2/0) | 2 | 82.7 / 69.0 | 148.6 | 0.059% | 0.86 / 1.08 | Pseudomonadota |
| RF02569 | IhtA | IhtA sRNA | 1 | 1 (1/0) | 1 | 120.8 / 76.0 | 32.7 | 0.013% | 0.08 / 0.22 | Chlamydiota |
| RF02600 | BASRCI27 | Brucella sRNA CI27 | 1 | 1 (1/0) | 1 | 63.9 / 59.8 | 61.4 | 0.024% | 0.11 / 0.33 | Pseudomonadota |
| RF02605 | scr5239 | Streptomyces sRNA scr5239 | 3 | 3 (3/0) | 3 | 172.1 / 64.0 | 98.2 | 0.039% | 0.15 / 0.45 | Actinomycetota |
| RF02608 | AbsR11 | Acinetobacter sRNA 11 | 1 | 1 (1/0) | 1 | 185.0 / 158.0 | 38.5 | 0.015% | 0.10 / 0.26 | Pseudomonadota |
| RF02647 | Cis52_sRNA | Ruegeria cis52 sRNA | 2 | 2 (2/0) | 2 | 129.1 / 56.0 | 1233.0 | 0.489% | 0.70 / 3.25 | Pseudomonadota |
| RF02656 | sot0042 | Sense overlapping transcript RNA 0042 (sot) | 1 | 1 (1/0) | 1 | 91.4 / 88.0 | 69.3 | 0.027% | 0.23 / 0.65 | Halobacteriota |
| RF02673 | scr4677 | Streptomyces sRNA 4677 | 6 | 3 (3/0) | 6 | 197.6 / 41.0 | 111.0 | 0.044% | 0.16 / 0.57 | Actinomycetota |
| RF02677 | NsiR4 | Nitrogen stress-induced RNA 4 | 1 | 1 (1/0) | 1 | 42.9 / 40.0 | 27.9 | 0.011% | 0.08 / 0.20 | Cyanobacteriota |
| RF02678 | Hatchet | Hatchet ribozyme | 1 | 1 (0/1) | 1 | 92.5 / 51.0 | 70.7 | 0.028% | 0.26 / 0.48 | phage |
| RF02695 | NrsZ | Nitrogen regulated small RNA | 5 | 4 (4/0) | 5 | 152.0 / 76.9 | 98.2 | 0.039% | 0.14 / 1.75 | Pseudomonadota |
| RF02703 | AsxR | Anti stx2 sRNA | 2 | 2 (0/2) | 2 | 88.1 / 63.0 | 24.0 | 0.010% | 0.08 / 0.19 | phage |
| RF02735 | Sernc350 | Saccharopolyspora sRNA 350 | 12 | 6 (6/0) | 12 | 168.1 / 50.0 | 1056.3 | 0.419% | 0.29 / 5.09 | Actinomycetota, Myxococcota |
| RF02784 | SorX | Singlet oxygen resistance RNA X | 1 | 1 (1/0) | 1 | 77.4 / 46.6 | 53.4 | 0.021% | 0.09 / 0.22 | Pseudomonadota |
| RF02790 | sodF_sRNA | sodF sRNA | 2 | 2 (2/0) | 2 | 74.2 / 46.0 | 70.5 | 0.028% | 0.09 / 0.23 | Actinomycetota |
| RF02791 | CcsR1 | Conserved CCUCCUCCC motif stress-induced RNA 1 | 1 | 1 (1/0) | 1 | 103.8 / 65.3 | 100.3 | 0.040% | 0.17 / 0.46 | Pseudomonadota |
| RF02794 | Pab19 | Pab19 RNA | 1 | 1 (1/0) | 1 | 93.8 / 55.0 | 95.2 | 0.038% | 0.15 / 0.41 | Methanobacteriota_B |
| RF02795 | Pab91 | Pab91 RNA | 1 | 1 (1/0) | 1 | 82.7 / 71.0 | 28.8 | 0.011% | 0.06 / 0.13 | Methanobacteriota_B |
| RF02796 | Pab160 | Pab160 RNA | 1 | 1 (1/0) | 1 | 93.0 / 47.4 | 39.5 | 0.016% | 0.09 / 0.22 | Methanobacteriota_B |
| RF02812 | asponA | anti ponA RNA | 4 | 4 (4/0) | 4 | 170.5 / 113.0 | 500.2 | 0.199% | 0.30 / 2.08 | Pseudomonadota |
| RF02827 | Scr6106 | Streptomyces RNA 6106 | 1 | 1 (1/0) | 1 | 61.2 / 52.5 | 181.0 | 0.072% | 0.19 / 0.90 | Actinomycetota |
| RF02828 | Scr3920 | Streptomyces RNA 3920 | 3 | 3 (3/0) | 3 | 151.8 / 88.0 | 102.7 | 0.041% | 0.12 / 0.46 | Actinomycetota |
| RF02829 | Scr4115 | Streptomyces RNA 4115 | 4 | 4 (4/0) | 4 | 102.1 / 58.0 | 144.3 | 0.057% | 0.17 / 0.68 | Actinomycetota |
| RF02830 | Scr1601 | Streptomyces RNA 1601 | 3 | 3 (3/0) | 3 | 170.4 / 100.0 | 140.5 | 0.056% | 0.17 / 0.69 | Actinomycetota |
| RF02831 | Scr2736 | Streptomyces RNA 2736 | 2 | 2 (2/0) | 2 | 95.2 / 43.0 | 56.9 | 0.023% | 0.09 / 0.24 | Actinomycetota |
| RF02832 | Scr5676 | Streptomyces RNA 5676 | 3 | 3 (3/0) | 3 | 143.9 / 56.0 | 126.7 | 0.050% | 0.20 / 0.54 | Actinomycetota |
| RF02833 | Scr3202 | Streptomyces RNA 3202 | 1 | 1 (1/0) | 1 | 77.7 / 52.0 | 88.7 | 0.035% | 0.10 / 0.39 | Actinomycetota |
| RF02847 | Ref64 | Enterococcus sRNA 64 | 1 | 1 (1/0) | 1 | 53.4 / 44.3 | 22.6 | 0.009% | 0.07 / 0.17 | Bacillota |
| RF02862 | sRNA_1300 | Enterococcus sRNA 1300 | 2 | 2 (2/0) | 2 | 266.3 / 48.8 | 78.7 | 0.031% | 0.20 / 0.59 | Bacillota |
| RF02864 | sRNA_0030 | Enterococcus sRNA 30 | 1 | 1 (1/0) | 1 | 194.9 / 123.0 | 46.1 | 0.018% | 0.13 / 0.31 | Bacillota |
| RF02872 | SrbA | sRNA regulator of biofilms A | 4 | 4 (4/0) | 4 | 132.9 / 69.0 | 170.5 | 0.068% | 0.14 / 0.59 | Pseudomonadota |
| RF02878 | MH_s7 | Mesorhizobail RNA 7 | 1 | 1 (1/0) | 1 | 244.9 / 221.0 | 308.4 | 0.122% | 0.24 / 1.53 | Pseudomonadota |
| RF02882 | MH_s36 | Mesorhizobail RNA 36 | 1 | 1 (1/0) | 1 | 66.3 / 46.0 | 39.7 | 0.016% | 0.09 / 0.22 | Pseudomonadota |
| RF02916 | atpB | atpB RNA | 2 | 2 (2/0) | 2 | 86.2 / 37.4 | 69.9 | 0.028% | 0.10 / 0.35 | Actinomycetota |
| RF02921 | RT-14 | RT-14 RNA | 1 | 1 (1/0) | 1 | 72.3 / 48.1 | 75.0 | 0.030% | 0.31 / 0.73 | Spirochaetota |
| RF02924 | skipping-rope | skipping-rope RNA | 2 | 2 (0/2) | 2 | 50.9 / 36.7 | 30.3 | 0.012% | 0.07 / 0.18 | phage |
| RF02935 | che1 | che1 RNA | 4 | 4 (4/0) | 4 | 97.7 / 54.1 | 83.7 | 0.033% | 0.18 / 0.41 | Actinomycetota |
| RF02944 | c4-2 | c4-2 RNA | 1 | 1 (1/0) | 1 | 47.4 / 41.0 | 26.7 | 0.011% | 0.08 / 0.18 | Desulfobacterota |
| RF02951 | DABA-DC-AT | DABA-DC-AT RNA | 1 | 1 (1/0) | 1 | 143.0 / 56.1 | 29.5 | 0.012% | 0.09 / 0.24 | Pseudomonadota |
| RF02965 | CyVA-1 | CyVA-1 RNA | 1 | 1 (1/0) | 1 | 56.4 / 41.4 | 44.8 | 0.018% | 0.13 / 0.36 | Cyanobacteriota |
| RF02982 | gltS | gltS RNA | 1 | 1 (1/0) | 1 | 74.5 / 47.5 | 54.5 | 0.022% | 0.15 / 0.37 | Pseudomonadota |
| RF02989 | gntR-DTE | gntR-DTE RNA | 1 | 1 (1/0) | 1 | 119.7 / 63.6 | 324.1 | 0.129% | 0.40 / 1.56 | Actinomycetota |
| RF02990 | gut-2 | gut-2 RNA | 1 | 1 (1/0) | 1 | 46.1 / 39.8 | 22.5 | 0.009% | 0.07 / 0.13 | Bacillota |
| RF02996 | int-alpA | int-alpA RNA | 1 | 1 (0/1) | 1 | 36.6 / 35.8 | 34.6 | 0.014% | 0.07 / 0.19 | phage |
| RF02999 | ivy-DE | ivy-DE RNA | 3 | 3 (3/0) | 3 | 93.5 / 41.4 | 43.8 | 0.017% | 0.06 / 0.19 | Pseudomonadota |
| RF03001 | leuA-Halobacteria | leuA-Halobacteria RNA | 1 | 1 (1/0) | 1 | 73.5 / 47.9 | 43.2 | 0.017% | 0.09 / 0.25 | Halobacteriota |
| RF03016 | RT-12 | RT-12 RNA | 1 | 1 (1/0) | 1 | 77.1 / 44.2 | 43.4 | 0.017% | 0.12 / 0.31 | Pseudomonadota |
| RF03020 | RT-17 | RT-17 RNA | 1 | 1 (1/0) | 1 | 44.0 / 40.4 | 45.2 | 0.018% | 0.15 / 0.33 | Bacteroidota |
| RF03023 | rpfG | rpfG RNA | 2 | 2 (2/0) | 2 | 101.3 / 75.9 | 68.3 | 0.027% | 0.12 / 0.42 | Pseudomonadota |
| RF03032 | narK | narK RNA | 4 | 3 (3/0) | 4 | 83.2 / 43.2 | 64.2 | 0.025% | 0.17 / 0.69 | Pseudomonadota |
| RF03034 | nrdJ | nrdJ RNA | 1 | 1 (1/0) | 1 | 81.9 / 40.3 | 29.7 | 0.012% | 0.08 / 0.20 | Actinomycetota |
| RF03044 | Proteo-phage-1 | Proteo-phage-1 RNA | 1 | 1 (1/0) | 1 | 42.2 / 42.0 | 48.4 | 0.019% | 0.11 / 0.31 | Pseudomonadota |
| RF03045 | proV | proV RNA | 1 | 1 (1/0) | 1 | 151.7 / 88.5 | 50.2 | 0.020% | 0.11 / 0.31 | Pseudomonadota |
| RF03046 | Pseudomonadales-1 | Pseudomonadales-1 RNA | 1 | 1 (1/0) | 1 | 118.6 / 68.0 | 42.6 | 0.017% | 0.13 / 0.33 | Pseudomonadota |
| RF03048 | RAGATH-20 | RAGATH-20 RNA | 1 | 1 (0/1) | 1 | 122.4 / 101.4 | 29.2 | 0.012% | 0.09 / 0.23 | phage |
| RF03057 | nhaA-I | nhaA-I RNA | 6 | 6 (6/0) | 6 | 53.8 / 37.0 | 37.6 | 0.015% | 0.08 / 0.21 | Hydrogenedentota, Myxococcota, Pseudomonadota |
| RF03058 | sul1 | sul1 RNA | 3 | 2 (2/0) | 3 | 80.2 / 36.2 | 55.7 | 0.022% | 0.12 / 0.39 | Pseudomonadota |
| RF03059 | raiA-hairpin | raiA-hairpin RNA | 3 | 3 (3/0) | 3 | 74.0 / 43.7 | 45.0 | 0.018% | 0.08 / 0.20 | Actinomycetota |
| RF03060 | uup | uup RNA | 3 | 2 (2/0) | 3 | 64.6 / 36.8 | 21.7 | 0.009% | 0.07 / 0.15 | Bacillota |
| RF03063 | Streptomyces-metK | Streptomyces-metK RNA | 3 | 3 (3/0) | 3 | 122.0 / 44.2 | 110.7 | 0.044% | 0.19 / 0.64 | Actinomycetota |
| RF03064 | RAGATH-18 | RAGATH-18 RNA | 4 | 2 (1/1) | 4 | 72.2 / 44.9 | 26.4 | 0.010% | 0.09 / 0.19 | Bacillota, phage |
| RF03067 | terC | terC RNA | 8 | 5 (5/0) | 8 | 58.7 / 39.9 | 29.3 | 0.012% | 0.06 / 0.17 | Pseudomonadota |
| RF03072 | raiA | raiA RNA | 1 | 1 (1/0) | 1 | 200.6 / 101.3 | 43.6 | 0.017% | 0.08 / 0.29 | Bacillota |
| RF03074 | Rhodo-rpoB | Rhodo-rpoB RNA | 1 | 1 (1/0) | 1 | 63.7 / 33.6 | 45.5 | 0.018% | 0.10 / 0.25 | Pseudomonadota |
| RF03076 | Streptomyces-metH | Streptomyces-metH RNA | 4 | 4 (4/0) | 4 | 108.2 / 50.0 | 85.7 | 0.034% | 0.24 / 0.57 | Actinomycetota |
| RF03078 | chrB-a | chrB-a RNA | 2 | 2 (2/0) | 2 | 68.5 / 43.9 | 38.4 | 0.015% | 0.08 / 0.17 | Pseudomonadota |
| RF03085 | abiF | abiF RNA | 3 | 3 (1/2) | 3 | 47.5 / 37.0 | 19.3 | 0.008% | 0.06 / 0.13 | Pseudomonadota, phage |
| RF03086 | chrB-b | chrB-b RNA | 1 | 1 (1/0) | 1 | 39.0 / 36.5 | 23.3 | 0.009% | 0.06 / 0.13 | Pseudomonadota |
| RF03097 | RAGATH-21 | RAGATH-21 RNA | 1 | 1 (0/1) | 1 | 83.9 / 70.3 | 39.5 | 0.016% | 0.10 / 0.31 | phage |
| RF03129 | L31-Corynebacteriaceae | L31-Corynebacteriaceae ribosomal protein leader | 1 | 1 (1/0) | 1 | 60.5 / 36.0 | 28.7 | 0.011% | 0.06 / 0.14 | Actinomycetota |
| RF03130 | L19-Flavobacteria | L19-Flavobacteria ribosomal protein leader | 1 | 1 (1/0) | 1 | 44.5 / 31.5 | 17.1 | 0.007% | 0.05 / 0.12 | Bacteroidota |
| RF03132 | L17-Actino-Proteobacteria | L17-Actino-Proteobacteria ribosomal protein leader | 4 | 4 (4/0) | 4 | 70.4 / 34.9 | 29.5 | 0.012% | 0.07 / 0.19 | Actinomycetota |
| RF03133 | S4-Fusobacteriales | S4-Fusobacteriales ribosomal protein leader | 1 | 1 (1/0) | 1 | 71.8 / 48.8 | 21.3 | 0.008% | 0.07 / 0.17 | Fusobacteriota |
| RF03134 | S6-S18-Chlorobi | S6-S18-Chlorobi ribosomal protein leader | 1 | 1 (1/0) | 1 | 69.3 / 47.4 | 40.5 | 0.016% | 0.08 / 0.23 | Bacteroidota_A |
| RF03135 | L4-Archaeoglobi | L4-Archaeoglobi ribosomal protein leades | 1 | 1 (1/0) | 1 | 67.1 / 34.5 | 44.2 | 0.018% | 0.07 / 0.23 | Halobacteriota |
| RF03137 | S15-Methanomicrobia | S15-Methanomicrobia ribosomal protein leader | 1 | 1 (1/0) | 1 | 58.6 / 42.2 | 23.6 | 0.009% | 0.06 / 0.16 | Halobacteriota |
| RF03138 | S16-Flavobacteria | S16-Flavobacteria ribosomal protein leader | 1 | 1 (1/0) | 1 | 59.2 / 35.7 | 19.4 | 0.008% | 0.05 / 0.13 | Bacteroidota |
| RF03141 | S4-Flavobacteria | S4-Flavobacteria ribosomal protein leader | 1 | 1 (1/0) | 1 | 51.5 / 35.2 | 18.6 | 0.007% | 0.06 / 0.15 | Bacteroidota |
| RF03142 | L2-Alphaproteobacteria | L2-Alphaproteobacteria ribosomal protein leader | 1 | 1 (1/0) | 1 | 52.6 / 37.8 | 24.2 | 0.010% | 0.06 / 0.16 | Pseudomonadota |
| RF03143 | S15-Halobacteria | S15-Halobacteria ribosomal protein leader | 1 | 1 (1/0) | 1 | 74.4 / 49.0 | 31.7 | 0.013% | 0.07 / 0.19 | Halobacteriota |
| RF03144 | eL15-Euryarchaeota | eL15-Euryarchaeota ribosomal protein leader | 2 | 2 (2/0) | 2 | 40.0 / 33.7 | 21.5 | 0.009% | 0.06 / 0.14 | Halobacteriota, Methanobacteriota_B |
| RF03146 | RAGATH-9 | RAGATH-9 RNA | 1 | 1 (0/1) | 1 | 142.5 / 75.2 | 55.5 | 0.022% | 0.19 / 0.45 | phage |
| RF03155 | TD-1 | TD-1 RNA | 35 | 1 (1/0) | 35 | 136.3 / 64.5 | 70.8 | 0.028% | 0.24 / 0.46 | Spirochaetota |
| RF03158 | L31-Actinobacteria | L31-Actinobacteria ribosomal protein leader | 6 | 6 (6/0) | 6 | 55.5 / 35.0 | 45.3 | 0.018% | 0.08 / 0.25 | Actinomycetota |
| RF03159 | RAGATH-7 | RAGATH-7 RNA | 1 | 1 (0/1) | 1 | 169.4 / 111.8 | 53.7 | 0.021% | 0.16 / 0.43 | phage |
| RF03161 | S15-Flavobacteria | S15-Flavobacteria ribosomal protein leader | 1 | 1 (1/0) | 1 | 71.8 / 46.5 | 41.9 | 0.017% | 0.17 / 0.33 | Bacteroidota |
| RF03162 | TD-2 | TD-2 RNA | 47 | 1 (1/0) | 47 | 92.2 / 48.7 | 47.8 | 0.019% | 0.17 / 0.31 | Spirochaetota |
| RF03163 | c4-a1b1 | a1, b1 targets of C4 antisense RNA | 42 | 7 (5/2) | 42 | 66.8 / 46.7 | 33.6 | 0.013% | 0.10 / 0.24 | Pseudomonadota, phage |
| RF04189 | icd-II | icd-II ncRNA motif | 1 | 1 (1/0) | 1 | 125.9 / 99.0 | 41.5 | 0.016% | 0.13 / 0.35 | Pseudomonadota |
| RF04192 | carA | carA ncRNA motif | 1 | 1 (1/0) | 1 | 212.0 / 100.0 | 102.1 | 0.041% | 0.46 / 0.73 | Pseudomonadota |
| RF04311 | dauA | dauA ncRNA motif | 1 | 1 (1/0) | 1 | 80.7 / 45.0 | 84.3 | 0.033% | 0.17 / 0.41 | Actinomycetota |

182 models; together 7.07% of this scan's per-model CPU (*derived*).

### 3. Excluded (c) models that hit

| accession | name | rule | hits | genomes (prok/phage) | new loci | max score / GA | share |
|---|---|---|---|---|---|---|---|
| RF00065 | snoR9 | c1 | 1 | 1 (1/0) | 1 | 108.3 / 90.0 | 0.028% |
| RF00560 | SNORA17 | c1 | 1 | 1 (1/0) | 1 | 32.1 / 29.0 | 0.016% |
| RF01856 | Protozoa_SRP | c7 | 2 | 2 (2/0) | 0 | 72.2 / 65.0 | 0.040% |
| RF02163 | sR-tMet | c1 | 1 | 1 (1/0) | 1 | 66.4 / 50.0 | 0.059% |
| RF04357 | Plastid-ndhA-ndhB-2 | c9 | 3 | 1 (1/0) | 0 | 60.9 / 53.0 | 0.009% |

Every hit (*measured*; gene: the gene calls the hit overlaps, either strand, with their Pfam domains; in genes: fraction of the hit's bases covered by gene calls):

| model | rule | genome | phylum | contig | start-end | strand | score - GA | E-value | in genes | genes overlapped (start-end strand: Pfam) | prok-1 hits overlapping |
|---|---|---|---|---|---|---|---|---|---|---|---|
| snoR9 | c1 | GCF_002214485.1 | Methanobacteriota_B | 2435a8abb338 | 1,773,433-1,773,560 | + | +18.3 | 4.7e-22 | 0.08 | 1,773,550-1,774,171 -: HAD+Hydrolase+Hydrolase_3 | - |
| SNORA17 | c1 | GCA_963676605.1 | Spirochaetota | 6c703a4dc3b5 | 6,037-6,138 | - | +3.1 | 4.9 | 1.00 | 5,353-6,541 -: Sulfatase | - |
| sR-tMet | c1 | GCF_001971705.1 | Halobacteriota | 43caf3ce2867 | 296,017-296,090 | - | +16.4 | 2e-08 | 0.00 | - | - |
| Protozoa_SRP | c7 | GCF_000220645.1 | Methanobacteriota | 3e60e709f0b6 | 1,545,870-1,546,113 | + | +7.2 | 1.1e-10 | 0.00 | - | Archaea_SRP, Bacteria_small_SRP |
| Protozoa_SRP | c7 | GCF_029633915.1 | Halobacteriota | f28b7ee3f008 | 2,081,706-2,081,944 | - | +1.4 | 3.1e-09 | 0.00 | - | Archaea_SRP |
| Plastid-ndhA-ndhB-2 | c9 | GCA_022701275.1 | Cyanobacteriota | 920e6dfb6344 | 170,280-170,321 | - | +7.9 | 3.6e-09 | 0.05 | 170,319-172,041 -: GIIM+HNH+RVT_1+RVT_N | Intron_gpII |
| Plastid-ndhA-ndhB-2 | c9 | GCA_022701275.1 | Cyanobacteriota | 920e6dfb6344 | 3,989,641-3,989,682 | + | +7.9 | 3.6e-09 | 0.05 | 3,987,921-3,989,643 +: GIIM+HNH+RVT_1+RVT_N | Intron_gpII |
| Plastid-ndhA-ndhB-2 | c9 | GCA_022701275.1 | Cyanobacteriota | 920e6dfb6344 | 5,037,711-5,037,752 | + | +3.0 | 8.7e-08 | 0.68 | 5,037,724-5,038,279 +: HNH+HNH_5 | Intron_gpII |

The eukaryotic rRNA models dropped by (a-), for the record:

| accession | name | hits | new loci (no prok-1 hit on the same strand) |
|---|---|---|---|
| RF00002 | 5_8S_rRNA | 109 | 0 |
| RF01960 | SSU_rRNA_eukarya | 313 | 0 |
| RF02542 | SSU_rRNA_microsporidia | 313 | 0 |
| RF02543 | LSU_rRNA_eukarya | 313 | 0 |

### 4. Rule-(b) models of prok-1 that never hit

CPU-s: this scan's per-model lines of the models with no hit (*measured*); h/10 Gbp: their cost in prok-1 as the producer scans them (*derived*, design-note formula).

| rule | models | with a hit | never hit | CPU-s of those | h/10 Gbp of those |
|---|---|---|---|---|---|
| b3 | 1 | 0 | 1 | 130.4 | 0.20 |
| b4 | 4 | 2 | 2 | 566.4 | 5.98 |
| b5 | 1 | 1 | 0 | 0.0 | 0.00 |
| b6 | 2 | 1 | 1 | 58.4 | 0.20 |
| b7 | 61 | 8 | 53 | 1200.7 | 3.72 |
| b8 | 6 | 6 | 0 | 0.0 | 0.00 |
| b9 | 33 | 23 | 10 | 318.7 | 0.87 |
| b10 | 188 | 23 | 165 | 14539.5 | 68.46 |
| b11 | 38 | 14 | 24 | 448.3 | 1.44 |

- (b3) never hit, 1: RF02357 RNaseP-T

- (b4) never hit, 2: RF01850 beta_tmRNA, RF01851 cyano_tmRNA

- (b6) never hit, 1: RF01693 Bacteroidales-1

- (b7) never hit, 53: RF01316 CRISPR-DR3, RF01320 CRISPR-DR7, RF01321 CRISPR-DR8, RF01322 CRISPR-DR9, RF01323 CRISPR-DR10, RF01324 CRISPR-DR11, RF01325 CRISPR-DR12, RF01326 CRISPR-DR13, RF01327 CRISPR-DR14, RF01328 CRISPR-DR17, RF01329 CRISPR-DR15, RF01331 CRISPR-DR18, RF01332 CRISPR-DR19, RF01333 CRISPR-DR20, RF01334 CRISPR-DR21, RF01336 CRISPR-DR23, RF01339 CRISPR-DR27, RF01340 CRISPR-DR29, RF01341 CRISPR-DR30, RF01342 CRISPR-DR32, RF01343 CRISPR-DR33, RF01344 CRISPR-DR34, RF01345 CRISPR-DR35, RF01346 CRISPR-DR36, RF01347 CRISPR-DR37, RF01349 CRISPR-DR40, RF01350 CRISPR-DR41, RF01351 CRISPR-DR42, RF01352 CRISPR-DR43, RF01353 CRISPR-DR44, RF01354 CRISPR-DR45, RF01355 CRISPR-DR26, RF01357 CRISPR-DR31, RF01358 CRISPR-DR39, RF01359 CRISPR-DR46, RF01360 CRISPR-DR47, RF01361 CRISPR-DR48, RF01362 CRISPR-DR49, RF01363 CRISPR-DR50, RF01364 CRISPR-DR51, RF01365 CRISPR-DR52, RF01366 CRISPR-DR53, RF01367 CRISPR-DR54, RF01368 CRISPR-DR55, RF01369 CRISPR-DR56, RF01370 CRISPR-DR57, RF01371 CRISPR-DR58, RF01373 CRISPR-DR60, RF01374 CRISPR-DR61, RF01376 CRISPR-DR63, RF01377 CRISPR-DR64, RF01378 CRISPR-DR65, RF01379 CRISPR-DR66

- (b9) never hit, 10: RF01510 2dG-I, RF01704 Glutamine-II, RF01767 SAM-III, RF01826 SAM_V, RF02680 PreQ1-III, RF02885 SAM_VI, RF03013 nadA, RF03071 Sodium, RF03164 NiCo-AnGGAG, RF03168 c-di-GMP-I-UAU

- (b10) never hit, 165: RF00615 LhrA, RF00616 LhrC, RF01456 VrrA, RF01457 rli22, RF01458 rli23, RF01459 rliE, RF01460 rliH, RF01461 rli24, RF01462 rli26, RF01463 rli27, RF01464 rliA, RF01465 rli31, RF01466 rli34, RF01467 rli36, RF01468 rli32, RF01469 rli33, RF01471 rliB, RF01472 rli40, RF01473 rli41, RF01475 rli45, RF01476 rliF, RF01477 rli43, RF01478 rli47, RF01479 rli48, RF01480 rli52, RF01481 rli53, RF01483 rli56, RF01484 rli59, RF01485 rli61, RF01486 rli62, RF01487 rliI, RF01488 rli49, RF01489 sbrA, RF01490 rli51, RF01491 rli54, RF01493 rli37, RF01494 rliD, RF01519 CC0196, RF01520 CC0734, RF01521 CC1840, RF01527 CrfA, RF01528 CC3513, RF01530 CC3664, RF01665 P13, RF01677 P35, RF01678 P37, RF01680 P5, RF01682 P8, RF01683 P6, RF01699 Clostridiales-1, RF01808 MicX, RF01843 neisseria_FSE, RF01867 CC2171, RF02048 STnc30, RF02049 STnc460, RF02050 STnc470, RF02054 STnc420, RF02055 STnc380, RF02056 STnc390, RF02059 STnc50, RF02062 STnc361, RF02063 STnc350, RF02065 STnc340, RF02066 STnc320, RF02067 STnc310, RF02070 STnc300, RF02071 STnc280, RF02072 STnc590, RF02073 STnc260, RF02075 STnc230, RF02077 STnc220, RF02078 STnc210, RF02080 STnc170, RF02088 STnc510, RF02270 nse_sRNA, RF02385 sau-13, RF02386 sau-19, RF02387 sau-27, RF02388 sau-30, RF02389 sau-31, RF02390 sau-41, RF02392 sau-53, RF02393 sau-59, RF02394 sau-63, RF02395 sau-66, RF02396 sau-5949, RF02397 sau-5971, RF02398 sau-6072, RF02404 P33, RF02414 rli60, RF02415 rliG, RF02418 Spd-sr07, RF02425 SpF01_sRNA, RF02426 SpF03_sRNA, RF02427 SpF10_sRNA, RF02428 SpF11_sRNA, RF02429 SpF14_sRNA, RF02430 SpF19_sRNA, RF02431 SpF22_sRNA, RF02432 SpF25_sRNA, RF02433 SpF36_sRNA, RF02434 SpF39_sRNA, RF02435 SpF41_sRNA, RF02436 SpF43_sRNA, RF02437 SpF44_sRNA, RF02438 SpF51_sRNA, RF02439 SpF56_sRNA, RF02441 SpF61_sRNA, RF02442 SpF66_sRNA, RF02443 SpR08_sRNA, RF02444 SpR10_sRNA, RF02445 SpR14_sRNA, RF02446 SpR18_sRNA, RF02447 SpR19_sRNA, RF02448 SpR20_sRNA, RF02467 Ms_IGR-2, RF02468 Ms_IGR-4, RF02470 Ms_IGR-8, RF02537 TarA, RF02538 TarB, RF02567 Va-907, RF02568 UptR, RF02577 tsr24, RF02578 tsr25, RF02579 tsr26, RF02580 tsr31, RF02581 tsr32, RF02582 tsr33, RF02583 Teg23, RF02627 Ssr1, RF02652 StyR-3, RF02653 StyR-143, RF02662 BLi_r0872, RF02676 Esr41, RF02736 Rev11, RF02737 Rev13, RF02738 Rev24, RF02739 Rev41, RF02740 Fwd6_3p_UTR, RF02741 Rev44, RF02742 Rev72, RF02745 Rev42_5p_UTR, RF02799 CssrB, RF02813 PA5194_thermometer, RF02818 V_AS5, RF02821 V_IGR5, RF02834 VqmR, RF02861 IR_84, RF02877 NmsRb, RF02887 mgtC_leader, RF02888 BtsR1, RF02889 Sr006, RF02894 S35, RF02895 S414, RF02896 S774, RF02897 S808, RF02909 rli117, RF02937 Clostridiales-2, RF02942 Clostridiales-3, RF03091 Clostridium-PBP, RF03112 Staphylococcus-1, RF03531 n00280, RF03533 SQ1002_RNA, RF03534 sRNA71, RF03535 TEG147

- (b11) never hit, 24: RF01119 sR32, RF01120 sR33, RF01121 sR38, RF01122 sR39, RF01123 sR35, RF01126 sR41, RF01127 sR42, RF01129 sR44, RF01131 sR47, RF01132 sR48, RF01134 sR30, RF01135 sR24, RF01142 sR19, RF01143 sR16, RF01146 sR15, RF01147 sR12, RF01149 sR10, RF01150 sR11, RF01152 sR1, RF01276 sR53, RF01297 sR40, RF01303 sR49, RF01306 sR52, RF01310 sR7

Rule-(a) models of prok-1 (a hit in the 8 golden genomes) with no hit in this sample: 165: RF00014 DsrA, RF00018 CsrB, RF00033 MicF, RF00034 RprA, RF00035 OxyS, RF00039 DicF, RF00040 rne5, RF00077 SraB, RF00078 MicA, RF00079 OmrA-B, RF00081 ArcZ, RF00082 SraG, RF00083 GlmZ_SraJ, RF00084 CsrC, RF00101 SraC_RyeA, RF00110 RybB, RF00112 CyaR_RyeE, RF00113 SIB_RNA, RF00115 McaS, RF00116 C0465, RF00117 C0719, RF00118 rydB, RF00119 C0299, RF00122 GadY, RF00124 IS102, RF00125 IS128, RF00126 ryfA, RF00128 GlmY_tke1, RF00368 sroB, RF00369 sroC, RF00370 sroD, RF00371 sroE, RF00372 sroH, RF00382 DnaX, RF00391 RtT, RF00505 RydC, RF00512 Leu_leader, RF00513 Trp_leader, RF00534 SgrS, RF00552 rncO, RF01054 preQ1-II, RF01056 Magnesium-II, RF01065 23S-methyl, RF01125 sR4, RF01137 sR21, RF01141 sR18, RF01144 sR17, RF01312 sR9, RF01335 CRISPR-DR22, RF01338 CRISPR-DR25, RF01400 istR, RF01401 rseX, RF01407 STnc560, RF01408 sraL, RF01470 rli38, RF01492 rli28, RF01517 iscRS, RF01691 Bacillus-plasmid, RF01707 JUMPstart, RF01708 L17DE, RF01709 Lacto-rpoB, RF01732 asd, RF01748 nuoG, RF01764 azaaromatic, RF01769 greA, RF01771 rnk_leader, RF01776 RatA, RF01779 AS1726, RF01780 AS1890, RF01783 b55, RF01796 frnS, RF01798 g2, RF01809 symR, RF01813 rdlD, RF01829 sR6, RF01830 StyR-44, RF01832 ROSE_2, RF01859 Phe_leader, RF01988 SECIS_2, RF02029 sraA, RF02030 tp2, RF02031 tpke11, RF02051 STnc450, RF02052 STnc630, RF02053 STnc430, RF02057 STnc40, RF02060 STnc410, RF02064 STnc370, RF02068 STnc480, RF02069 STnc70, RF02074 STnc240, RF02079 STnc180, RF02081 STnc550, RF02082 STnc540, RF02083 OrzO-P, RF02084 STnc130, RF02099 rivX, RF02273 FsrA, RF02341 ncrMT1302, RF02348 tracrRNA, RF02370 Trp_leader_2, RF02371 PyrG_leader, RF02372 PyrC_leader, RF02373 PyrD_leader, RF02377 SurA, RF02378 SurC, RF02384 FasX, RF02401 ClpQY_promoter, RF02419 Spd-sr37, RF02449 ncr1015, RF02450 ncr1175, RF02451 ncr1241, RF02453 ncr952, RF02454 ncr982, RF02465 Ms_AS-5, RF02495 ohsC_RNA, RF02524 sagA, RF02525 SSRC30, RF02526 SSRC34_1, RF02527 SSRC38, RF02528 SSRC41, RF02529 SSRC8_1, RF02564 naRNA4, RF02566 Ms1, RF02571 Mcr7, RF02590 Spy1186876, RF02591 Spy1786666, RF02592 Spy392987, RF02616 SSRC8_2, RF02617 SSRC10, RF02618 SSRC34_2, RF02640 MOSES4, RF02641 Spy490483c, RF02642 Spy491311c, RF02643 Spy491738, RF02644 Spy490380c, RF02654 MicL, RF02659 ncRv12659, RF02694 RalA, RF02700 HtrA_thermometer, RF02744 Rev39_5p_UTR, RF02746 AsrC, RF02750 ES003, RF02751 ES036, RF02752 ES056, RF02753 ES173, RF02754 ES205, RF02755 ES222, RF02756 ES239, RF02757 Erse, RF02768 Ysr155_RyfD, RF02773 TrxA_thermometer, RF02782 CpoB_ybgF_thermometer, RF02797 PssrA, RF02798 CssrA, RF02822 Srn266, RF02850 Ysr276, RF02855 Ysr251, RF02886 npcTB_6715, RF02891 AgrB, RF02892 SR6_antitoxin, RF02929 algC, RF03039 Peptidase-S11, RF03131 L25-Gammaproteobacteria, RF04312 trpB

### 5. Hits inside genes, by family

Of 11,912 hits, 442 lie wholly and 664 at least half inside gene calls; 412 do not lie wholly inside an intergenic region with the 200-base margin (*measured* hits against the stage's gene calls; regions by `biobox.l1.ncrna.regions`). Families with any such hit:

| accession | name | class | whole genome in prok-1 | hits | wholly in genes | mostly in genes | outside regions | h/10 Gbp wg / igr |
|---|---|---|---|---|---|---|---|---|
| RF00001 | 5S_rRNA | prok-1 (a) | no | 326 | 10 | 14 | 1 | 0.25 / 0.13 |
| RF00002 | 5_8S_rRNA | (a-) | - | 109 | 10 | 12 | 8 | 0.92 / 0.49 |
| RF00005 | tRNA | prok-1 (a) | no | 5893 | 170 | 253 | 19 | 0.69 / 0.60 |
| RF00010 | RNaseP_bact_a | prok-1 (a) | no | 81 | 1 | 6 | 4 | 1.96 / 0.70 |
| RF00013 | 6S | prok-1 (a) | no | 46 | 0 | 1 | 0 | 0.35 / 0.13 |
| RF00023 | tmRNA | prok-1 (a) | no | 92 | 0 | 3 | 3 | 6.99 / 5.59 |
| RF00028 | Intron_gpI | prok-1 (a) | no | 30 | 0 | 2 | 2 | 3.73 / 1.79 |
| RF00029 | Intron_gpII | prok-1 (a) | no | 128 | 6 | 15 | 4 | 0.14 / 0.07 |
| RF00043 | Plasmid_R1162 | d | - | 1 | 0 | 1 | 0 | 0.27 / 0.10 |
| RF00050 | FMN | prok-1 (a) | no | 42 | 0 | 1 | 0 | 0.33 / 0.13 |
| RF00059 | TPP | prok-1 (a) | no | 113 | 5 | 7 | 0 | 0.25 / 0.12 |
| RF00080 | Manganese | prok-1 (a) | no | 11 | 1 | 2 | 0 | 0.41 / 0.19 |
| RF00114 | S15 | prok-1 (a) | no | 8 | 0 | 4 | 0 | 1.65 / 0.78 |
| RF00140 | Alpha_RBS | prok-1 (a) | no | 9 | 0 | 1 | 0 | 0.39 / 0.16 |
| RF00162 | SAM | prok-1 (a) | no | 86 | 1 | 3 | 0 | 0.33 / 0.16 |
| RF00169 | Bacteria_small_SRP | prok-1 (a) | no | 86 | 4 | 6 | 0 | 0.25 / 0.11 |
| RF00174 | Cobalamin | prok-1 (a) | no | 188 | 3 | 4 | 3 | 0.38 / 0.21 |
| RF00177 | SSU_rRNA_bacteria | prok-1 (a) | no | 313 | 0 | 0 | 8 | 6.63 / 5.27 |
| RF00230 | T-box | prok-1 (a) | no | 119 | 0 | 4 | 0 | 0.94 / 0.75 |
| RF00373 | RNaseP_arch | prok-1 (a) | no | 12 | 0 | 2 | 1 | 0.96 / 0.53 |
| RF00379 | c-di-AMP | prok-1 (a) | no | 36 | 0 | 2 | 0 | 0.24 / 0.11 |
| RF00504 | Glycine | prok-1 (a) | no | 41 | 1 | 1 | 0 | 0.27 / 0.11 |
| RF00516 | ylbH | prok-1 (a) | no | 3 | 1 | 1 | 0 | 0.57 / 0.20 |
| RF00560 | SNORA17 | c1 | - | 1 | 1 | 1 | 1 | 0.27 / 0.10 |
| RF01051 | c-di-GMP-I | prok-1 (b9) | no | 41 | 0 | 1 | 0 | 0.23 / 0.09 |
| RF01055 | Molybdenum | prok-1 (a) | no | 14 | 0 | 1 | 0 | 0.27 / 0.11 |
| RF01130 | sR46 | prok-1 (b11) | no | 1 | 0 | 1 | 0 | 0.26 / 0.10 |
| RF01274 | sR45 | prok-1 (b11) | no | 2 | 0 | 1 | 0 | 0.24 / 0.08 |
| RF01317 | CRISPR-DR4 | prok-1 (a) | no | 110 | 2 | 2 | 0 | 0.12 / 0.06 |
| RF01411 | BsrF | prok-1 (a) | no | 3 | 0 | 2 | 0 | 0.33 / 0.09 |
| RF01419 | IsrR | d | - | 2 | 2 | 2 | 2 | 0.89 / 0.35 |
| RF01497 | ALIL | d | - | 16 | 13 | 14 | 11 | 0.57 / 0.26 |
| RF01668 | P10 | prok-1 (b10) | no | 2 | 1 | 1 | 1 | 0.46 / 0.12 |
| RF01670 | P17 | prok-1 (b10) | no | 3 | 0 | 3 | 3 | 1.84 / 0.47 |
| RF01671 | P18 | prok-1 (b10) | no | 3 | 3 | 3 | 0 | 3.25 / 0.48 |
| RF01673 | PhrS | d | - | 5 | 0 | 5 | 0 | 2.34 / 0.95 |
| RF01690 | Bacillaceae-1 | prok-1 (a) | no | 97 | 4 | 12 | 1 | 0.23 / 0.10 |
| RF01695 | C4 | d | - | 104 | 1 | 3 | 0 | 0.42 / 0.13 |
| RF01696 | Chlorobi-1 | d | - | 3 | 2 | 2 | 0 | 0.18 / 0.08 |
| RF01719 | Pseudomon-1 | d | - | 5 | 1 | 1 | 0 | 0.21 / 0.11 |
| RF01731 | TwoAYGGAY | d | - | 53 | 1 | 8 | 1 | 0.35 / 0.18 |
| RF01762 | whalefall-1 | d | - | 1 | 1 | 1 | 0 | 0.15 / 0.07 |
| RF01766 | cspA | prok-1 (a) | no | 13 | 0 | 11 | 0 | 3.53 / 1.32 |
| RF01781 | ASdes | prok-1 (a) | yes | 7 | 7 | 7 | 7 | 0.60 / 0.20 |
| RF01782 | ASpks | prok-1 (a) | yes | 92 | 90 | 90 | 90 | 0.68 / 0.17 |
| RF01804 | Lambda_thermo | prok-1 (a) | no | 3 | 0 | 3 | 0 | 0.25 / 0.10 |
| RF01831 | THF | prok-1 (b9) | no | 3 | 1 | 1 | 0 | 0.17 / 0.08 |
| RF01852 | tRNA-Sec | prok-1 (a) | no | 55 | 6 | 11 | 2 | 0.69 / 0.43 |
| RF01854 | Bacteria_large_SRP | prok-1 (a) | no | 29 | 0 | 4 | 0 | 0.56 / 0.21 |
| RF01857 | Archaea_SRP | prok-1 (a) | no | 9 | 0 | 1 | 0 | 0.72 / 0.27 |
| RF01959 | SSU_rRNA_archaea | prok-1 (a) | no | 313 | 0 | 0 | 8 | 7.50 / 6.15 |
| RF01960 | SSU_rRNA_eukarya | (a-) | - | 313 | 0 | 0 | 8 | 15.65 / 10.83 |
| RF01998 | group-II-D1D4-1 | prok-1 (b8) | no | 76 | 21 | 21 | 3 | 0.23 / 0.09 |
| RF01999 | group-II-D1D4-2 | prok-1 (b8) | no | 11 | 2 | 2 | 0 | 0.44 / 0.14 |
| RF02001 | group-II-D1D4-3 | prok-1 (b8) | no | 99 | 7 | 13 | 4 | 0.60 / 0.24 |
| RF02003 | group-II-D1D4-4 | prok-1 (a) | no | 19 | 0 | 1 | 0 | 0.25 / 0.11 |
| RF02005 | group-II-D1D4-6 | prok-1 (b8) | no | 35 | 2 | 2 | 0 | 0.56 / 0.27 |
| RF02111 | IS009 | prok-1 (a) | no | 1 | 1 | 1 | 1 | 2.16 / 0.49 |
| RF02194 | HPnc0260 | prok-1 (a) | yes | 1 | 0 | 1 | 0 | 3.20 / 2.46 |
| RF02225 | sX6 | d | - | 2 | 0 | 2 | 0 | 5.63 / 0.90 |
| RF02228 | sX9 | d | - | 81 | 2 | 4 | 2 | 0.39 / 0.11 |
| RF02345 | ar15 | d | - | 2 | 1 | 1 | 1 | 0.24 / 0.10 |
| RF02463 | Ms_AS-1 | prok-1 (b10) | no | 1 | 1 | 1 | 0 | 1.53 / 0.23 |
| RF02510 | PYLIS_3 | d | - | 2 | 2 | 2 | 2 | 0.50 / 0.16 |
| RF02514 | 5_ureB_sRNA | prok-1 (a) | yes | 24 | 24 | 24 | 24 | 6.35 / 0.40 |
| RF02540 | LSU_rRNA_archaea | prok-1 (a) | no | 313 | 0 | 1 | 53 | 23.48 / 19.17 |
| RF02541 | LSU_rRNA_bacteria | prok-1 (a) | no | 313 | 0 | 1 | 53 | 21.86 / 17.31 |
| RF02542 | SSU_rRNA_microsporidia | (a-) | - | 313 | 0 | 0 | 8 | 6.89 / 5.19 |
| RF02543 | LSU_rRNA_eukarya | (a-) | - | 313 | 0 | 1 | 46 | 48.32 / 32.22 |
| RF02589 | Spy779816 | prok-1 (a) | no | 2 | 0 | 1 | 0 | 0.25 / 0.11 |
| RF02647 | Cis52_sRNA | d | - | 2 | 2 | 2 | 2 | 3.25 / 0.70 |
| RF02656 | sot0042 | d | - | 1 | 0 | 1 | 0 | 0.65 / 0.23 |
| RF02674 | AsdA | prok-1 (a) | yes | 2 | 2 | 2 | 2 | 0.43 / 0.13 |
| RF02735 | Sernc350 | d | - | 12 | 12 | 12 | 12 | 5.09 / 0.29 |
| RF02795 | Pab91 | d | - | 1 | 0 | 1 | 0 | 0.13 / 0.06 |
| RF02812 | asponA | d | - | 4 | 0 | 4 | 4 | 2.08 / 0.30 |
| RF02819 | V_AS7 | prok-1 (b10) | no | 1 | 1 | 1 | 0 | 0.93 / 0.20 |
| RF02820 | V_AS9 | prok-1 (b10) | no | 2 | 2 | 2 | 2 | 0.36 / 0.14 |
| RF02832 | Scr5676 | d | - | 3 | 0 | 2 | 0 | 0.54 / 0.20 |
| RF02833 | Scr3202 | d | - | 1 | 1 | 1 | 0 | 0.39 / 0.10 |
| RF02882 | MH_s36 | d | - | 1 | 0 | 1 | 0 | 0.22 / 0.09 |
| RF02908 | rli30 | prok-1 (b10) | no | 1 | 1 | 1 | 1 | 0.42 / 0.16 |
| RF02916 | atpB | d | - | 2 | 0 | 1 | 0 | 0.35 / 0.10 |
| RF02951 | DABA-DC-AT | d | - | 1 | 0 | 1 | 0 | 0.24 / 0.09 |
| RF03057 | nhaA-I | d | - | 6 | 0 | 1 | 0 | 0.21 / 0.08 |
| RF03058 | sul1 | d | - | 3 | 0 | 1 | 0 | 0.39 / 0.12 |
| RF03059 | raiA-hairpin | d | - | 3 | 2 | 2 | 0 | 0.20 / 0.08 |
| RF03135 | L4-Archaeoglobi | d | - | 1 | 1 | 1 | 1 | 0.23 / 0.07 |
| RF03146 | RAGATH-9 | d | - | 1 | 0 | 1 | 0 | 0.45 / 0.19 |
| RF03155 | TD-1 | d | - | 35 | 1 | 1 | 1 | 0.46 / 0.24 |
| RF03161 | S15-Flavobacteria | d | - | 1 | 1 | 1 | 0 | 0.33 / 0.17 |
| RF03162 | TD-2 | d | - | 47 | 3 | 5 | 2 | 0.31 / 0.17 |
| RF03167 | c-di-GMP-I-GGC | prok-1 (b9) | no | 21 | 0 | 1 | 0 | 0.24 / 0.07 |
| RF04183 | Bacteroidales_small_SRP | prok-1 (b5) | no | 5 | 1 | 1 | 0 | 0.22 / 0.08 |
| RF04357 | Plastid-ndhA-ndhB-2 | c9 | - | 3 | 0 | 1 | 0 | 0.16 / 0.07 |

The 23 whole-genome models of prok-1 in this sample:

| accession | name | hits | outside regions |
|---|---|---|---|
| RF00382 | DnaX | 0 | 0 |
| RF01779 | AS1726 | 0 | 0 |
| RF01781 | ASdes | 7 | 7 |
| RF01782 | ASpks | 92 | 90 |
| RF01843 | neisseria_FSE | 0 | 0 |
| RF01988 | SECIS_2 | 0 | 0 |
| RF01989 | SECIS_3 | 4 | 0 |
| RF02194 | HPnc0260 | 1 | 0 |
| RF02376 | antisense_SR1 | 2 | 0 |
| RF02465 | Ms_AS-5 | 0 | 0 |
| RF02514 | 5_ureB_sRNA | 24 | 24 |
| RF02592 | Spy392987 | 0 | 0 |
| RF02654 | MicL | 0 | 0 |
| RF02662 | BLi_r0872 | 0 | 0 |
| RF02674 | AsdA | 2 | 2 |
| RF02736 | Rev11 | 0 | 0 |
| RF02737 | Rev13 | 0 | 0 |
| RF02738 | Rev24 | 0 | 0 |
| RF02741 | Rev44 | 0 | 0 |
| RF02742 | Rev72 | 0 | 0 |
| RF02746 | AsrC | 0 | 0 |
| RF02757 | Erse | 0 | 0 |
| RF02850 | Ysr276 | 0 | 0 |

### 6. Proposed prok-2

`planning/design/0-08a_prok-2_proposed.tsv`: prok-1's 588 models with their flags, plus the 182 long-tail models of item 2; a member is scanned whole genome if prok-1 says so or if more than half of its hits here lie outside every margin-200 region (11: RF01419 IsrR, RF01497 ALIL, RF01670 P17, RF02111 IS009, RF02510 PYLIS_3, RF02647 Cis52_sRNA, RF02735 Sernc350, RF02812 asponA, RF02820 V_AS9, RF02908 rli30, RF03135 L4-Archaeoglobi). 770 models, 34 whole genome.

Cost (*derived*): CPU-hours per 10 Gbp = sum over the set of the model's 8-genome CPU line (margin 200 `--notrunc` if intergenic, whole genome if whole genome) / 18,106,231 x 1e10 / 3,600 gives 247.4 for prok-1 (the note's 247.4) and 291.4 for prok-2; scaled to prok-1's measured 254.5: 254.5 x 291.4 / 247.4 = **299.7**.


<!-- results:end -->

## Method

**The sample** (`rfam_sample.py sample`; reproducible from the seed 20261112 and the sorted accession lists):

- Prokaryotes: the 2,497 genomes of stages P1 (250) and P2 (2,247) of `/tmp/if2453/r10/root`, read through the
  current code (`genome_input`, `contig_index`). A genome's accession is the file name of the source object at
  the same index of its stage's `objects.json` as the genome's entry in `genomes.json` (the two lists were made
  together by `make_inputs.py` and have equal lengths); for every sampled genome this pairing was checked by
  hashing its FASTA records to the genome's contig uids. Each genome's GTDB R226 phylum is read from
  `bac120_metadata_r226.tsv.gz` and `ar53_metadata_r226.tsv.gz` (`gtdb_taxonomy`), matched by accession (all
  2,497 matched exactly). Allocation of the 100 over the 55 phyla present (8 archaeal): one each, the remaining
  45 in proportion to the phylum's genome count (largest remainder, ties by phylum name, never more than a phylum
  has); within a phylum `random.Random("0-08a:20261112:<domain>:<phylum>").sample(<sorted accessions>, k)`.
- Phages: the 27,258 genomes of S1 sorted by `genome_uid` (hex), every 90th from the first (`27,258 // 300`).
- The sample's 676 contigs (none circular) were exported to `/tmp/if2453/m08a/fa/s<k>.fna` (`rfam_sample.py
  fasta`, through `SeqStore`; every record's BLAKE3 contig uid rechecked against its name: 676 records,
  412,032,974 bases, 0 mismatches).
- Gene calls on those contigs (`l1/feature`, producer `pyrodigal-3.7.1-gv-0.3.2-p2`, stages S1, P1, P2) and the
  Pfam 38.2 domains of their proteins (`l1/feature_protein`, `l1/domains`) were exported to
  `sample/genes.jsonl` (`rfam_sample.py genes`): 379,579 genes on all 676 contigs, 289,697 with a Pfam domain
  (*measured*).

**The plan.** Genomes in the order prokaryotes by accession, then phages by genome_uid, packed greedily into
processes of at most 18,000,000 bases (a genome never split): 27 processes of 8,981,933 to 17,990,498 bases
(`sample/plan.json`). At 1,450 bases per CPU-second (0-09c, *measured*) the largest is 207 CPU-minutes, and 227
at E-01's slower 1,319 (*derived*), under the 300 CPU-minute limit (AGENTS §7.2).

**The scan** is the producer's whole-contig argv (`biobox.l1.ncrna._cmsearch_argv`: `cmsearch --cpu 2 -Z 100000
--cut_ga --noali --tformat fasta --tblout <tbl>`; not `--notrunc`, not intergenic) with `-o <file>` instead of
`/dev/null`, so that each model's own `# CPU time` line is kept, over all 4,227 models of `Rfam.cm`, one process
per planned part. Hits are mapped to contigs by `infernal_igr.py collect` (the producer's rules) and per-model CPU
read by `infernal_igr.py permodel` (both from card 0-09c). CPU is user + system of the `cmsearch` process
(`timed.py`).

**Pod limits kept.** Each process is a unit `m08a-s<k>` with `MemoryMax=16G`, `MemorySwapMax=0`,
`CPUQuota=200%`, `nice 19`, in the image; at most 2 at a time (4 threads), beside the 50 Gbp build. Before every
launch the driver sums `ps -u if2453 -o pcpu=` (the ps line left out) and launches only when that plus 400% is at
most 4000% (the account's day cap, and its cap whenever the load average is above 8); otherwise it polls every 4
minutes. Nothing is written under `/tmp/if2453/r10/root` or `/stor`; the driver only reads for a `STOP` file there.

## Commands

On the laptop: `git -C <worktree> checkout -b measure/0-08a origin/plan/task-cards` (at `4648286`); the tree was
mirrored to `/tmp/if2453/ci/m08a` by `sh scripts/pod_pytest.sh m08a tests/r0/test_runbook.py` (34 passed), and on
the pod copied to the frozen `/tmp/if2453/m08a/tree` (`rsync -a --delete`), from which every step runs. `IMG` is
`/tmp/if2453/images/bio-box_80d4c55.sif`, `BIN` is `/tmp/if2453/m08a/tree/scripts/stress/annotation_cost`, and
each python step runs as
`systemd-run --user --scope -p MemoryMax=8G -p MemorySwapMax=0 -p CPUQuota=100% nice -n 19 singularity exec
--cleanenv --bind /tmp/if2453 --env PYTHONPATH=/tmp/if2453/m08a/tree/src $IMG python $BIN/timed.py <record> -- python ...`:

```bash
sha256sum /tmp/if2453/golden/rfam/Rfam.cm.gz > /tmp/if2453/m08a/lib/Rfam.cm.gz.sha256
gunzip -c /tmp/if2453/golden/rfam/Rfam.cm.gz > /tmp/if2453/m08a/lib/Rfam.cm
python $BIN/rfam_sample.py sample /tmp/if2453/r10/root /tmp/if2453/r10/inputs /tmp/if2453/b18/gtdb /tmp/if2453/m08a   # 0.66 CPU-min
python $BIN/rfam_sample.py fasta  /tmp/if2453/r10/root /tmp/if2453/m08a /tmp/if2453/m08a/lib/Rfam.cm                   # 0.16 CPU-min
python $BIN/rfam_sample.py genes  /tmp/if2453/r10/root /tmp/if2453/m08a                                                 # 0.23 CPU-min
systemd-run --user --unit=m08a-driver --collect --quiet -p MemoryMax=4G -p MemorySwapMax=0 -p CPUQuota=100% \
  -p StandardOutput=append:/tmp/if2453/m08a/driver.log -p StandardError=append:/tmp/if2453/m08a/driver.log \
  nice -n 19 bash $BIN/rfam_sample_driver.sh
```

The driver (`rfam_sample_driver.sh`) launches each process as `rfam_sample_shard.sh <k>` in unit `m08a-s<k>`:
`timed.py runs/s<k>.json -- cmsearch --cpu 2 -Z 100000 --cut_ga --noali --tformat fasta --tblout out/s<k>.tbl -o
out/s<k>.out lib/Rfam.cm fa/s<k>.fna`, then `infernal_igr.py collect out/s<k>.tbl fa/s<k>.map.json
out/s<k>.hits.jsonl` and `infernal_igr.py permodel out/s<k>.out out/s<k>.permodel.json`, then `out/s<k>.ok`. When
all 27 are done it runs `rfam_sample_report.py /tmp/if2453/m08a /tmp/if2453/m08a/tree /tmp/if2453/m08a/report`
and writes `DONE`.

## How to finish (if the scan was still running when this was committed)

1. Is it done? On the pod: `ls /tmp/if2453/m08a/DONE /tmp/if2453/m08a/FAILED`; progress is
   `tail /tmp/if2453/m08a/driver.log` and `ls /tmp/if2453/m08a/out/*.ok | wc -l` (27 when complete);
   `systemctl --user list-units 'm08a-*'` shows the driver and the running processes.
2. `FAILED` names the failed processes (their logs: `logs/s<k>.log`; `out/s<k>.fail` holds the exit code). Run the
   driver again (the `systemd-run` line above): done processes are skipped, failed ones retried. To pause: `touch
   /tmp/if2453/m08a/HALT` (running processes finish; remove the file and relaunch the driver to resume). To use 4
   processes once the build has finished and the pod is otherwise idle: `echo 4 > /tmp/if2453/m08a/MAXP`.
3. When `DONE` exists, `/tmp/if2453/m08a/report/` holds `report.json`, `results.md`, `0-08a_prok-2_proposed.tsv`
   and this file with its results filled in (`2026-11-12_rfam_sample.md`). Copy the last two into the repository
   (`planning/measurements/` and `planning/design/`), write the reading paragraph below, and commit. To rebuild
   the report by hand (also mid-scan with `--partial`), run the scope line above with `python
   $BIN/rfam_sample_report.py /tmp/if2453/m08a /tmp/if2453/m08a/tree /tmp/if2453/m08a/report [--partial]`.

## Reading

(The orchestrating session's interpretation, 2026-11-12; the numbers are the report's, *measured* unless labelled.)

1. **`prok-1` finds about 92% of the loci on a phylum-stratified sample.** Of 11,912 hits at or above GA, 10,003 are
   `prok-1` models' (*derived*: 11,912 − 853 long-tail − 8 excluded − 1,048 of the four eukaryotic rRNA models, which
   sit on bacterial rRNA loci that `prok-1`'s own rRNA models find). The long tail adds 850 loci no `prok-1` model
   covers: 10,003 of 10,853 distinct loci is 92.2% (*derived*). The design note's 99.4% was measured on the 8 golden
   genomes, from which the set was drawn; this is the first figure on genomes the set has not seen.
2. **The long tail does not saturate at 100 genomes.** 182 long-tail models hit, 91 of them in one genome only, 14 in
   five or more; the 100th prokaryote still adds 1.92 new families on average, and the curve from 80 to 100 genomes
   (309 to 341) has not bent. Only 93 of the 8 golden genomes' 258 families hit here at all. A set frozen from any
   sample of this size will keep missing family-specific sRNAs in lineages the sample did not hold; at 12,000 species
   the number of further models with a hit is unknown, and could be several hundred (*estimate*, from the slope).
3. **What the tail is.** Mostly lineage-specific sRNAs and leaders (Pseudomonadota: RsmY, RsmZ, PrrF, Qrr, RgsA;
   archaeal Hgc RNAs; cyanobacterial Yfr), a few large structured RNAs (OLE), and 20 models with phage hits. None
   shares a clan with a `prok-1` hit at the same locus, so they are additions, not relabels.
4. **`prok-2` costs 18% more.** 770 models (34 scanned on whole genomes) at 299.7 CPU-hours per 10 Gbp against 254.5
   (*derived*). Eleven of the added models have most of their hits outside every intergenic region and so join the
   whole-genome list.
5. **Rule (b) models that never hit are not dead weight.** 165 of the 188 "sRNAs of model prokaryotes" models had
   no hit here and cost 68.5 CPU-hours per 10 Gbp (*derived*), but they are species-specific (Listeria `rli`,
   Caulobacter `CC`, Bacillus and Staphylococcus sRNAs) and the sample holds one genome per phylum stratum, not those
   species. Dropping them would save 27% of the cost and lose them where they exist. Not recommended.
6. **Excluded models behaved as expected.** The five excluded models with a hit (8 hits) are eukaryotic snoRNAs, a
   protozoan SRP and a plastid motif at or just above GA; nothing argues for readmitting them.

**For the PI.** DEC-28 says the long-tail models that hit on this sample are added once and the set's version
changes: that is `prok-2` (`planning/design/0-08a_prok-2_proposed.tsv`). The sample adds one fact DEC-28 did not have:
the tail keeps growing with the corpus, so "once" freezes a floor, not a ceiling. Three ways to go: (1) adopt `prok-2`
as decided (recommended: +18% cost, about 8% more loci on this sample, a stable set); (2) stay on `prok-1` and accept
92%; (3) replace the sample rule by a class rule (every Rfam family whose seed is bacterial, archaeal or phage), which
needs a costing first: the 3,969 models with no hit on the golden genomes carry 84% of full Rfam's 1,703 CPU-hours per
10 Gbp. Nothing in the built stages has ncRNA yet, so the choice costs no rework.

## What was not measured, and why

- Intergenic cost of the added models: the scan is whole genome (the truth scan); the cost of prok-2 is *derived*
  from the 8-genome per-model lines of 0-09c, as prok-1's was.
- Whether each hit of an excluded or long-tail model is a real RNA: the report gives coordinates, score against
  GA, overlapping gene calls with their Pfam domains and overlapping prok-1 hits; the judgement is the PI's.
- Clan competition is not applied: every hit at or above GA is counted (card E-08 applies it). The report says
  which candidate hits overlap a prok-1 hit of the same clan.
- Genomes outside P1 and P2, metagenome contigs, and stage S1's phages beyond the 300 sampled.

## Files

Scripts in `scripts/stress/annotation_cost/`: `rfam_sample.py` (sample, export, gene calls),
`rfam_sample_driver.sh` and `rfam_sample_shard.sh` (the scan), `rfam_sample_report.py` (every table), reusing
`infernal_igr.py` and `timed.py`. Outputs stay on the pod under `/tmp/if2453/m08a` (`sample/`, `fa/` the exported
FASTA, `runs/`, `out/` tables, gzipped `-o` outputs, hits and per-model CPU, `logs/`, `report/`).
