# The compact `tower` sets `DPN1`, `SYM2`, `LRP1`, `MIR2`, `MGE4`: what they hold, and their staging (card C-78, DEC-79, DEC-80)

- **Date and host:** 2026-11-26, `wilkcomp02` and `tower`. Every number carries a label (AGENTS §6): *measured* (on
  the host named, 2026-11-26, by the command or file named), *derived* (formula and inputs), *estimate*, or
  *published*.
- **Script:** `scripts/stress/build_stage_pod/c78_compact_tower_sets.py` (`census`, `extract`, `lrp1`, `mir2`, `sym2`,
  `dpn1-*`); tests `tests/r1/test_c78_compact_tower_sets.py`; templates
  `scripts/stress/build_stage_pod/stages/templates/{LRP1,MIR2,SYM2,DPN1}.conf`.
- **Read only:** `tower:/mnt/user/research-bulk/bio-box/sources_incoming/{DPN1,SYM2,LRP1,MIR2,MGE4}/` (`ls -la`,
  `SOURCE_URLS.tsv`, `SHA256SUMS`, `unzip -l`, `tar -tvzf`, `grep` of FASTA headers), the data root
  `/tmp/if2453/r10/root`, `/tmp/if2453/s2/inputs/` (the `r26.11.4` stages and `GV1`), the staged `r26.11.5` stages
  (`/tmp/if2453/c75/meta/GRARE{1,2,3}`, `/tmp/if2453/c66/{LR2,NX1}`, `/tmp/if2453/c58/SAG1`, `/tmp/if2453/c69/EGV1`,
  `/tmp/if2453/c77/meta/SAG2`), card C-74's index `/tmp/if2453/c74/index/`, card C-75's archaeal references. Nothing
  was written to the root, to `/tmp/if2453/s2/inputs/` or under `sources_incoming/` on `tower`; nothing under `/stor`
  was overwritten.
- **Written:** relayed downloads (64 MiB and more) to `/stor/scratch/Wilke/if2453/bio-box-staging/<STAGE>/dl/`, smaller
  ones to `/tmp/if2453/c78/dl/<STAGE>/`; staged FASTA to `/stor/scratch/Wilke/if2453/bio-box-staging/<STAGE>/`;
  metadata to `/tmp/if2453/c78/meta/<STAGE>/`; reports to `/tmp/if2453/c78/work/`.

## 1. What each set holds (the card's first deliverable)

Read on `tower` (`ls -la`, `SOURCE_URLS.tsv`, `unzip -l`, `tar -tvzf`, `grep -c '>'`; *measured*, 2026-11-26 19:00 to
19:10 CST) and, for the files relayed, by `c78_compact_tower_sets.py census` on the pod (every member named; a FASTA
member is `dna` when every record is IUPAC nucleotide letters, else `protein`, `alignment` or `not_fasta`; report
`/tmp/if2453/c78/work/census*.json`).

| Set | Files on `tower` (bytes, *measured*) | What they hold (*measured*) | The plan's warning | Verdict |
|---|---|---|---|---|
| `MGE4` | `Donor_strain-pKJK5_AtraI.fasta` 4,696,601; `RHB23_traJ_reconverted.fasta` 5,052,879; `Transconjugant.fasta` 4,543,652; `Samples_name.pdf` 31,942 | three isolate genome assemblies, one per file (headers `EPI005213_contig_1..3`, `EPI005211_contig_1..2`, `EPI005212_contig_1..117`): the donor 4,550,895 + 54,066 + 14,584 bp (a chromosome, a plasmid the size of pKJK5, about 54 kb, and a small contig), the recipient RHB23 4,860,572 + 109,433 bp, the transconjugant 117 contigs, 4,466,681 bp, the largest 418,022 bp | "lab strains of one experiment rather than mobile elements" | **true; not staged** (§2) |
| `LRP1` | `metagenome_plasmidome_plasmid2.fasta` 249,389,715 | **78,738 contigs**, 241,335,107 bp, 1,000 to 224,924 bp each: plasmidome (circular-DNA-enriched metagenome) assemblies of many samples (`DP-Sample<n>-...`, `18097D-03-<n>-...`), SPAdes `NODE_` names; 2,503 headers say `_type_circular`, 2,896 are `RNODE_` records | none | **staged** (§3); the record says "10,942 plasmids" (plan §4.6, survey 1e), the file holds 78,738 contigs, 69,004 of them under 5 kb |
| `MIR2` | `genomes_polinton_virophages.zip` 8,293,636; `PLV_framework_supplementary_data.tar.gz` 265,729,817; `protist_genomes_assembled.tar.gz` 831,184,917; `EVE_annotations_table.zip` 26,563,179; `README.txt` 1,238 | the zip: one FASTA, 1,193 PLV and virophage genomes, 27,836,611 bp (6,044 to 63,893 bp; names are protist-assembly coordinates, `N-HFG070301-C1323_2898-26282`, as `README.txt` says); the framework tar: 966 members, 19 of them genome DNA FASTA, 59,078 records, 633,726,614 bp (`fasta_unsorted_PLV/IMGVR_PLV_all.fasta`, 52,658 IMG/VR UViGs, 548,922,055 bp; `fasta_unsorted_PLV/Austrian_lake_metavirome_PLVs.fasta`, 3,328 contigs, 25,263,516 bp; 17 `annotated_PLVs/*.fasta`, 3,092 records, 59.5 Mbp, many of them copies of the unsorted ones), the rest protein FASTA (9), alignments (294), HMMs, trees, GFF, structures, tables; the protist tar: eleven `<species>_1kb.fasta.gz` assemblies of cryptophytes and chrysophytes (12,495,681 to 194,758,266 bytes each); `EVE_annotations_table.zip` a table | "the large archive of `MIR2` looks like protist host genomes" | **true for `protist_genomes_assembled.tar.gz`** (eukaryote host assemblies, left for card C-60); the zip and the framework's DNA are viruses: **staged** (§4) |
| `SYM2` | `assemblies_fna.tar.gz` 2,789,033,552; `bins_fna.tar.gz` 161,076,450; `viruses.fna.tar.gz` 113,889,853; `source_data.tar.gz` 254,624,379 | `assemblies_fna`: 453 whole assemblies of single-protist libraries (`LKH<n>_P<plate>_<Am, Ci, Co, Rh>.fna`: 114 `Am`, 262 `Ci`, 19 `Co`, 58 `Rh`; 9,640,806,649 bytes uncompressed, one file empty, the largest 161,530,688), the hosts' and their symbionts' DNA together; `bins_fna`: **625 bins**, 27,816 contigs, 515,578,694 bp; `viruses.fna`: one FASTA, **64,460 viral contigs**, 389,958,735 bp, 500 to 459,844 bp, grouped by library (`LKH041_P01_Ci_virus|NODE_...`); `source_data`: the paper's tables, notebooks, trees and figures (one DNA FASTA of 813 18S rRNA genes) | "the large archive of `SYM2` looks like protist host genomes" | **true for `assemblies_fna.tar.gz`** (left for card C-60); bins and viral contigs **staged** (§5) |
| `DPN1` | `Nanobdellota_Zen.zip` 9,013,602,602; `3662_DPANN_archaeal_genomes.tar.gz` 976,420,873; `DPANN_viruses_and_unclassified_MGEs_assembled_from_metagenomes.fna` 517,244; two metadata tables (493,344 and 1,660) | the tar: **3,662 genome FASTA** (`3662.DPANN.genomes/<GCA_..., IMG id>_genomic.fna`), 3,265,594,791 bytes; the MGE FASTA: 26 contigs; the zip: 98,490 entries, of which **genome FASTA are 259 files** (`genomes/complete_208/`, 208 closed genomes, and `genomes/noncircular_51/`, 51; 318,637,440 bytes) and the 243 NCBI genomes the study compared them with (`ncbi_all/hq_genomes/NCBI_GCA_*.fna`, 243,293,749 bytes); everything else is analysis: an ESMFold structure archive `deposit_artifacts/esmfold_pdbs.tar.gz` (6,051,186,772 bytes) and its Foldseek database (0.31 GB), OrthoFinder output (47,123 protein `.fa`, 42,859 `.tsv`, 2.44 GB), KofamScan tables (2.35 GB), proteomes (`.faa`, about 1.0 GB), genes (`.ffn`, `.gff`), trees, alignments, HMMs, figures | "9.99 GB is far more than 256 DPANN genomes and their elements" | **true**: 0.56 GB of the zip is genome FASTA and 8.4 GB is structures, orthologue tables and annotations (§6) |

## 2. `MGE4` leaves `r26.11.5`

The three FASTA files are the whole-genome assemblies of one conjugation experiment (donor `EPI005213` with the
plasmid pKJK5, recipient `RHB23` with `traJ`, transconjugant `EPI005212`; record `zenodo.22794128`; `Samples_name.pdf`
is the sample sheet). They are isolate genomes of lab strains, not a catalogue of mobile elements; the survey's
"PIMEs + capsid satellites + AVGD" (plan §10) is not what landed. By the card's rule `MGE4` is not staged and leaves
`r26.11.5` (DEC-79 (1): it moves to `r26.11.6` without a row; here there is nothing of it to stage there either). The
plasmid contigs (`EPI005213_contig_2`, 54,066 bp; `EPI005211_contig_2`, 109,433 bp) are not split out: the card names
no such rule. Nothing of `MGE4` was relayed to the pod. In the chain, `MV2` takes `AFTER=MIR2` (plan §5).

## 3. `LRP1`

**Bytes.** `metagenome_plasmidome_plasmid2.fasta` to `/stor/scratch/Wilke/if2453/bio-box-staging/LRP1/dl/` by
`scripts/orchestrator/relay_to_pod.sh` (card C-73), `--dry-run` first, the list `tower`'s `SHA256SUMS` (copied with
`scp`): 249,389,715 bytes in 9 s (27.71 MB/s), sha256 `26fc65b932b9b9276ce7966877732f4163a8ba5facba52d8de3deda3db18bf65`
equal to the record (*measured*, the relay's log, 19:08 CST).

**Rule** (`lrp1`): every DNA record is one genome, `kind="plasmid"`, code 11, header as deposited; DEC-65 (1)
(metagenome-assembled): an exact 20..1,000 bp terminal repeat is trimmed once and the record staged circular
(`assembly_quality` `complete`); a record whose header says `circular` and has no repeat would be circular as it
stands (none was); every other record is linear and whole (`contig`). A record whose sequence repeats an earlier one is
listed once (DEC-27).

| | Value (*measured*, `lrp1`, report `/tmp/if2453/c78/work/lrp1_report.json`) |
|---|---:|
| Records read | 78,738 |
| Genomes (one record each) | 78,737 (1 exact repeat skipped: `skipped_records.tsv`) |
| Bases | 240,689,693 bp |
| Terminal repeats (`c48c_virus_inputs.ends`, classes) | `exact` 4,025 trimmed (639,721 bp removed), `partial` 1,505, `long` 448, `low_complexity` 1, `none` 72,758 (one more `exact` in the repeat skipped) |
| Circular genomes | 4,025 (of 2,503 headers that say `_type_circular`; every circular one had an exact repeat) |
| Records under 5,000 bp | 69,004 (kept: plan §2 keeps plasmids whole) |
| Held by the root's 13 stored stages | 0 |
| Held by the `r26.11.4` stages (exact `contig_uid`) | **13 contigs, 40,741 bp**: `PL3` 9 (33,788 bp), `PL2` 2 (2,536 bp), `MGE1` 2 (4,417 bp); `ECE1`, `LPL1`, `PL1` 0 |
| Held only by an `r26.11.5` stage before it | 2 contigs, 3,104 bp: `GRARE2` 1, `NX1` 1 (`later_records.tsv`) |
| New | 78,724 contigs, 240,648,952 bp (0.241 Gbp; 78,722 and 240,645,848 bp after the later holders) |
| `s2_check_inputs.py` (with `objects_staging.json`) | exit 0, `problems: []`, 78,737 genomes and `genome_input` rows, kind `plasmid`, code 11 |

`ECE2` (card C-70) was being restaged while this ran, so `LRP1`'s exact overlap with it is card A-29's (chain order),
as are `DPN1` and `SYM2`, which come before `LRP1` in the chain and were staged after it here (smallest first). Against
the plan: new bases 0.241 Gbp against "under 0.3" (*estimate*); the record count differs (78,738 contigs against
"10,942 plasmids", survey 1e), which is recorded, not resolved: the file is the assembly, not a curated plasmid list.
`redistributable: false` (DEC-80 (3): "CC BY 4.0 as the job list states it", no provider record read).

| File | Where | Bytes | sha256 |
|---|---|---:|---|
| `lrp1_urban_water_plasmidome.fna` | `/stor/scratch/Wilke/if2453/bio-box-staging/LRP1/` | 248,738,773 | `e28a44bb6d26789208705986ae052f03444fe2585f6169da539365483598797a` |
| `genomes.json` | `/tmp/if2453/c78/meta/LRP1/` | 11,168,080 | `146a72757c2a32d18766a25b8e40abfc7154439f036920e998f1bf4de9754437` |
| `objects.json` | same | 488 | `aae43066322dd3b1b648786f11665f373639217ab2cace1a4fa8dab6f439cafe` |
| `objects_staging.json` | same | 509 | `98e7e4cf836f29f2df764c50fa6860512f48ea07612a192f565701b54c4d54a7` |
| `records.tsv` | same | 18,687,479 | `0067605abb6f36c6fe8c6ef41cd197afd27bccdacc03cfbdd9a65fd3c5e43501` |
| `trims.tsv` | same | 940,022 | `731533bc7dfcd91f0a752f94fd76c94cebea1cad5083b513e0ef69a506ff3c3a` |
| `shared_records.tsv` | same | 1,418 | `9c6f28b954fa5900a2907a37520d091a4e131782281bdc2c90ee7b2d57f9cbee` |
| `later_records.tsv` | same | 235 | `97598d022c0e071b0fdc7b9387f053df253a9b505538bb6751a4a4ab1cb86311` |
| `skipped_records.tsv` | same | 216 | `6e2e816b44ce436f1ab98ed4a14fb4d9aa8dde26efbe1e2414a2d57e6e1b1a4e` |
| `repeated_records.tsv` (header only) | same | 39 | `d109495ee76e4643c181eae183c452c9a84560f97a79bbe22d6630780368c565` |

All *measured* (`ids.sha256_file`, the script's report). `objects.json` names
`/tmp/if2453/s2/inputs/LRP1/lrp1_urban_water_plasmidome.fna` (card A-29 copies there). Template: `DEPTH=annotated`,
`LIBRARY=pfam 38.2 gathering`, `CRISPR=yes`, `AFTER=SYM2`, `EXPECT=held_by` for `MGE1` 2, `PL2` 2, `PL3` 9.
Build cost (*derived*, plan §3: R2 for meta-mode calling, about 93.5 to 135.3, plus R3 2,090 to 2,740 CPU-minutes per
new Gbp, × 0.2406): about **525 to 692 CPU-minutes**, inside the plan's "at most about 820".
`stage_est.py conf` on the template (`--tree`, in the image): exit 0 (*measured*, unit `c78-est-lrp1`).

**Mirror** (*measured*, the relay's log, 19:26 CST): `scripts/orchestrator/relay_to_tower.sh`, `--dry-run` first, two
runs to `tower:/mnt/user/research-bulk/bio-box/sources/LRP1/`: the FASTA (248,738,773 B in 9 s, 27.64 MB/s) and the nine
metadata files of the table; every sha256 on `tower` equal to the table's, `DONE: 1 verified` and `DONE: 9 verified, 0
failed`; `SHA256SUMS` there has 10 lines.

## 4. `MIR2`

**Bytes** (*measured*, the relay's log, 19:08 to 19:09 CST): `genomes_polinton_virophages.zip` (8,293,636 B) and
`README.txt` to `/tmp/if2453/c78/dl/MIR2/`, `PLV_framework_supplementary_data.tar.gz` (265,729,817 B, 9.16 MB/s) to
`/stor/scratch/Wilke/if2453/bio-box-staging/MIR2/dl/`; each sha256 equal to `tower`'s `SHA256SUMS`.
`protist_genomes_assembled.tar.gz` (the eleven protist assemblies) and `EVE_annotations_table.zip` (a table) were
listed on `tower` and not relayed.

**Rule** (`mir2`): every member of the two archives is read in archive order; a FASTA member whose first record is DNA
is staged, one `kind="virus"` genome per record (header as deposited; DEC-65 (1) as for `LRP1`, no source flag); every
other member is listed in `skipped_members.tsv` with its kind (alignment 294, protein 9, not FASTA 626, directories
18). A record whose sequence an earlier one already has is listed in `skipped_records.tsv` (DEC-27); most such are the
`annotated_PLVs/` copies of UViGs the unsorted file already holds.

| | Value (*measured*, `mir2`, report `/tmp/if2453/c78/work/mir2_report.json`) |
|---|---:|
| Bulzu 2025 (`genomes_polinton_virophages.fasta`) | 1,193 genomes, 27,836,543 bp |
| Bellas 2026 framework | 54,899 genomes, 558,072,697 bp, from 59,078 DNA records (4,179 exact repeats skipped) |
| of which `IMGVR_PLV_all.fasta` | 51,461 genomes, 530,519,848 bp (IMG/VR UViGs that carry a PLV-type major capsid protein, the authors' unsorted search) |
| of which `Austrian_lake_metavirome_PLVs.fasta` | 3,321 genomes, 25,178,268 bp |
| of which `annotated_PLVs/*.fasta` (17 files) | 117 genomes, 2,374,581 bp not in the unsorted files |
| All | 56,092 genomes, 585,909,240 bp; 1,335 exact terminal repeats trimmed (152,007 bp) and staged circular |
| Held by the root (exact `contig_uid`) | `G1` 1 (11,479 bp) |
| Held by the `r26.11.4` stages | `MV1` 3,792 (61,462,988 bp), `V2` 1,522 (35,565,939 bp) |
| Held only by an `r26.11.5` stage before it | `GV1` 95 (1,156,898 bp); `MIR1` (card C-70) was not on the pod: card A-29's |
| New | 50,777 contigs, 488,868,834 bp (0.489 Gbp; 487,711,936 bp after `GV1`) |
| `s2_check_inputs.py` | exit 0, `problems: []`, 56,092 genomes, kind `virus` (`gv-meta`), code 11 |

Against the plan: 0.489 new Gbp against "under 0.5 with `MGE4`" (*estimate*; `MGE4` adds nothing). Most of the bases
(530.5 of 585.9 Mbp, 90.5 %, *derived*) are IMG/VR UViGs; 5,314 of MIR2's contigs are already held by `MV1` or `V2`
by exact `contig_uid` (which file each came from is in `shared_records.tsv` and `record_sources.tsv`), and the other
UViGs are ones `V2` (one high-quality representative per vOTU) did not take; the UViGs' provider's
terms are JGI's (as entry `jgi_img_vr`), the framework's record states CC BY 4.0; both stay `redistributable: false`.

**A breach of AGENTS §7.3, reported:** the Bulzu FASTA (`mir2_bulzu_2025_plv_virophage.fna`, 28,333,081 B) was written
to `/stor/scratch/Wilke/if2453/bio-box-staging/MIR2/` with the framework FASTA, although it is under 64 MiB. It is left
where it is (nothing under `/stor` is deleted) and named there by `objects_staging.json`; the script was then changed
(§9) so that a FASTA under 64 MiB is written beside the metadata on `/tmp`, as `SYM2` and `DPN1` were.

| File | Where | Bytes | sha256 |
|---|---|---:|---|
| `mir2_bellas_2026_plv_framework.fna` | `/stor/scratch/Wilke/if2453/bio-box-staging/MIR2/` | 561,748,052 | `88ea664a7de971be13b927811a371585467634dd1353383b3fd9a1ec23d566db` |
| `mir2_bulzu_2025_plv_virophage.fna` | same | 28,333,081 | `382c77a48fa82d8cdffd2ff73204b6ce5895425bde570c15df075325d523f0b0` |
| `genomes.json` | `/tmp/if2453/c78/meta/MIR2/` | 7,788,756 | `10639d1bc059d555f84670b7be37e74240527836d08c588b24a1f697c64bdcff` |
| `objects.json` | same | 1,251 | `a89224942974e7181a7328166eeb1ff79c73864b72aac28b1e7af1a1cb9b1101` |
| `objects_staging.json` | same | 1,293 | `4f2ac996696a9d2031f0c3bd5151200c9eb5f32a35e506aa9db90a72e3ab1e01` |
| `records.tsv` | same | 13,025,188 | `a8809546681f43c7f3f0a35f675f0240bfa7ae6b9d2021c075b03c82129fb951` |
| `record_sources.tsv` | same | 7,138,391 | `0015e61b0ea437436bfd5154dc17e4562ee872e951a14077b9b35505f5acf3da` |
| `trims.tsv` | same | 274,079 | `ea8b6bd81d176ed0136d3f30d8cd1ebcee5859fc29bb7679c186e0bd9ad8bb23` |
| `shared_records.tsv` | same | 535,250 | `32d961115b8843c9e6f92788ca9c47cb3837e1fce6f072e1df9e6f76de7bad7f` |
| `later_records.tsv` | same | 9,652 | `45aab991ac0bae07620e853c6cbc902a9dc0a4ec9a2b60192f0681e86cfad85d` |
| `skipped_records.tsv` | same | 781,529 | `5dd336f7f0552525bafe62374b895bab15d1ffb3c5ab1b6ad71b3093e4c487af` |
| `skipped_members.tsv` | same | 81,510 | `05b5530c70e4536612b1bec356eb56b51c2dc037fed67b548cb42771966cbe20` |
| `repeated_records.tsv` (header only) | same | 39 | `d109495ee76e4643c181eae183c452c9a84560f97a79bbe22d6630780368c565` |

All *measured*. Template: `DEPTH=annotated`, `GENECALL_MODE=gv-meta`, `LIBRARY=pfam 38.2 gathering`, `CRISPR=yes`,
`AFTER=LRP1`, `EXPECT=held_by` `G1` 1, `V2` 1,522, `MV1` 3,792. Build cost (*derived*, plan §3 R2 93.5 to 135.3 plus R3
2,090 to 2,740 CPU-minutes per new Gbp, × 0.4889): about **1,067 to 1,406 CPU-minutes**, the plan's "at most about
1,400" at its upper end. `stage_est.py conf`: exit 0 (*measured*, unit `c78-est-mir2`).

**Mirror** (*measured*, the relay's log, 19:33 CST): two runs of `relay_to_tower.sh` to
`tower:/mnt/user/research-bulk/bio-box/sources/MIR2/`: the two FASTA and the eleven metadata files of the table;
`DONE: 2 verified` and `DONE: 11 verified, 0 failed`; `SHA256SUMS` there has 13 lines.

**Left on `tower`** (`CATALOG_GOAL.md` §4.1): `protist_genomes_assembled.tar.gz` (eleven cryptophyte and chrysophyte
assemblies of at least 1 kb, including their bacterial and viral contigs, as `README.txt` says) for card C-60;
`EVE_annotations_table.zip` (a table). Inside the framework archive, the non-DNA members (`skipped_members.tsv`).

## 5. `SYM2`

**Licence.** DEC-80 lifted the hold: no record states a licence ("open-public-NERSC" is the job list's access label,
R-11); the objects say so and carry `redistributable: false`.

**Bytes** (*measured*, the relay's log, 19:09 to 19:10 CST): `bins_fna.tar.gz`, `viruses.fna.tar.gz` and
`source_data.tar.gz` (the paper's tables, for the bins' metadata and the viruses' final list) to
`/stor/scratch/Wilke/if2453/bio-box-staging/SYM2/dl/`, each sha256 equal to `tower`'s `SHA256SUMS`;
`assemblies_fna.tar.gz` (the host assemblies) was listed on `tower` and not relayed. Twelve tables were copied out of
`source_data.tar.gz` byte for byte to `/tmp/if2453/c78/dl/SYM2/tables/` (`extract`); two are read: `metadata_mags.tab`
(999 MAGs) and `viruses/final_virus_df.tsv` (the authors' 7,392 complete and high-quality viruses, 44,382,967 bp).

**Rule** (`sym2`): each bin is one genome, its records whole, quality from the table's MIMAG (`HQ` `mag_hq`, `MQ`
`mag_mq`, else `contig`); a bin whose `symbclass` is *Nucleocytoviricota* is a giant-virus genome, `kind="virus"`;
every other bin is `kind="prokaryote"` with its code screened as `SYM1`'s (card C-58's `screen_genetic_code`, the
whole genome, the authors' GTDB classification as the hint); each record of `viruses.fna` is one `kind="virus"`
genome under DEC-65 (1) (the authors' `topology` `DTR` as the source's circular flag), except a contig that a bin
already holds (a contig must have one kind, `s2_check_inputs.py`) and an exact repeat (DEC-27).

**A correction in the run.** The first run (unit `c78-sym2`) staged all 625 bins as prokaryotes; reading its
`sym2_bins.tsv` showed 74 bins with `symbclass` *Nucleocytoviricota* (22,908,692 bp, classified "Unclassified" by
GTDB-Tk). The script was changed (a test first, seen failing), the first run's metadata moved aside on `/tmp`
(`/tmp/if2453/c78/meta/SYM2.run1/`), and the stage run again (`c78-sym2b`) with `--reuse-fasta`: the two FASTA
already on `/stor` had the same bytes as the new ones and were kept as they were, never rewritten (`fasta_reused` in
the report). Only `genomes.json`, `genetic_codes.tsv`, `r232_status.tsv` and `sym2_bins.tsv` changed.

| | Value (*measured*, `sym2`, report `/tmp/if2453/c78/work/sym2b_report.json`) |
|---|---:|
| Bins | 625, 27,816 records, 515,578,694 bp (576 in `metadata_mags.tab`; 49 not, all small "Unclassified" bins, 19,687,698 bp) |
| of which prokaryote | 551, 492,670,002 bp (*derived*): Patescibacteria 215, free-living 201, Chlamydiae 38, alpha-symbionts 22, gamma-symbionts 15, TM6 10, *Amoebophilus* 1, no class 49 (the table's `symbclass`); MIMAG HQ 68, MQ 155, LQ 353 (all 625) |
| of which giant viruses (`kind="virus"`) | 74, 22,908,692 bp |
| Codes of the 551 prokaryote bins | 11: 546; 25: 5; 4: 0 |
| Viral contigs read | 64,460 |
| Viral contigs staged | 61,959, 352,787,865 bp: 7,183 of the authors' final set (5,430 complete, 27,800,417 bp; 1,753 high-quality, 14,337,452 bp) and 54,776 others (310,649,996 bp) |
| Viral contigs not staged | 2,296 held by a bin (`contig_in_bin`), 205 exact repeats (`skipped_records.tsv`) |
| DEC-65 (1) | 8,303 exact terminal repeats trimmed (972,295 bp), staged circular |
| Genomes | 62,584 (625 bins and 61,959 viral contigs); 89,775 `genome_input` rows; 31 contigs listed by two bins (DEC-27) |
| Held by the root, the `r26.11.4` stages (`SYM1` included) | **0** |
| Held only by an `r26.11.5` stage before it | `GV1` 1,050 viral contigs, 20,968,682 bp |
| New | 89,741 contigs, 867,530,342 bp (0.868 Gbp; 846,561,660 after `GV1`) |
| r232 status | all 551 prokaryote bins `not_visible` (no INSDC accession) |
| `s2_check_inputs.py` | exit 0, `problems: []`, kinds `prokaryote`, `virus`, codes 11, 25 |

Against the plan: 0.847 to 0.868 new Gbp against 0.5 to 0.8 (*estimate*, survey 1d), 6 to 9 % over its upper end
(*derived*), inside the 20 % of AGENTS §6. The card names `viruses.fna`; it holds every viral contig the authors
called, of which their final set is 7,392: all are staged, and `sym2_viruses.tsv` marks which are in the final set
(`in_final_set`), so a later card can narrow it. Build cost (*derived*, plan §3: R1 24.6 to 24.8 × 0.4927 plus R2
93.5 to 135.3 × 0.3757 plus R3 2,090 to 2,740 × 0.8675): about **1,860 to 2,440 CPU-minutes**, the plan's 1,000 to
2,200 at and over its upper end.

| File | Where | Bytes | sha256 |
|---|---|---:|---|
| `sym2_protist_symbiont_bins.fna` | `/stor/scratch/Wilke/if2453/bio-box-staging/SYM2/` | 525,681,148 | `9b26ace425ff5f2fa9b9748490bbb44aa377c51b8e8885d20db3c5bd0a3ff9ba` |
| `sym2_protist_viruses.fna` | same | 356,338,487 | `103f6ca82ea23ef8c7243cdb07c213e48a4d0f8262570a88fe70b6ead2e3622a` |
| `genomes.json` | `/tmp/if2453/c78/meta/SYM2/` | 9,917,131 | `0a0ce08bd55413796ec631c532875a227777df1b06d9392d12fc8332565aaca4` |
| `objects.json` | same | 1,126 | `2fdfcd03d2a93459f35ffb00cafdb7693ac6b8cdb6d92785fe62fc99127f44f0` |
| `objects_staging.json` | same | 1,168 | `c57ea6124aa241f4bc33060a671f4d6d19a2bb90c10f675c99d85e9bd0b2149a` |
| `records.tsv` | same | 19,182,535 | `0808c8336cbe1c84876accff4684d2e72f0872965c2286a87130eb32f54d39c3` |
| `sym2_bins.tsv` | same | 123,897 | `68010d5d8bb284543831a921f9798dcaaeaf0b76ef4cad3bf3ffa9cdf268bfea` |
| `sym2_viruses.tsv` | same | 5,610,101 | `75deb3fc0914f09d3b0dc42717bf9b2aa796d35e8965247c7ce6158d3b39a335` |
| `genetic_codes.tsv` | same | 31,850 | `3aae03aaefd75af4e711dd46a760e960676603fb83dbbb76bfc3b78f52074539` |
| `r232_status.tsv` | same | 25,279 | `3545a597ae8f62546807fe20790d63b1c64fb5198206d9f23e5cb7d00fe4f6a5` |
| `trims.tsv` | same | 1,665,861 | `91b4992c6b80ce34f66404f5336dfdc06b186dda5b4f50d50bb17ada5bae03a6` |
| `later_records.tsv` | same | 66,017 | `ac05f1aeb8a59048a8bd6d7007dfef81d272811f6a9eff003c817cd15e060457` |
| `skipped_records.tsv` | same | 287,993 | `5884e9eb5e4176daefdba884571d431ffa886bed1b871fc8c14bc848500545ed` |
| `skipped_members.tsv` | same | 54 | `a95f28ca51520a9655e223a73649ab327c2f4ad8045ccfdf04edd56be0623de7` |
| `shared_records.tsv` (header only) | same | 38 | `2f9af65a096679c5376b81c7ecf6f62bb2ef3b99b972693d68c9dbb7586534df` |
| `repeated_records.tsv` (header only) | same | 39 | `d109495ee76e4643c181eae183c452c9a84560f97a79bbe22d6630780368c565` |

All *measured*. Template: `DEPTH=annotated`, `LIBRARY=pfam 38.2 gathering`, `CRISPR=yes`, `AFTER=DPN1`, no `EXPECT`
line (nothing held by the root or `r26.11.4`).

**Mirror** (*measured*, the relay's log, 20:21 CST): two runs of `relay_to_tower.sh` to
`tower:/mnt/user/research-bulk/bio-box/sources/SYM2/`: the two FASTA and the fourteen metadata files of the table;
`DONE: 2 verified` and `DONE: 14 verified, 0 failed`; `SHA256SUMS` there has 16 lines. The first run's metadata
(`SYM2.run1/`) was not mirrored.

**Left on `tower`** (`CATALOG_GOAL.md` §4.1): `assemblies_fna.tar.gz`, the 453 single-protist library assemblies
(9.64 GB uncompressed), host and symbiont DNA together, for card C-60 (its non-eukaryote part is what the bins and
viruses already took); `source_data.tar.gz`'s notebooks, trees and figures (not sequence).
