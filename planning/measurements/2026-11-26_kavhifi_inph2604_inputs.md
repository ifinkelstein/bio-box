# Kavagutti's HiFi MAGs and NCLDVs (`KAVHIFI`) and INPHARED's 7Apr2026 delta (`INPH2604`): selection, staging, overlap (card C-80, DEC-79, DEC-80)

- **Date and host:** 2026-11-26, `wilkcomp02`, 14:03 to 14:25 CST. Every number carries a label (AGENTS §6): *measured*
  (on `wilkcomp02`, 2026-11-26, by the command named in the section), *derived* (formula and inputs) or *published*.
- **Script:** `scripts/stress/build_stage_pod/c80_kavagutti_inphared.py` (`kav`, `inph`); tests
  `tests/r1/test_c80_kavagutti_inphared.py`; templates `scripts/stress/build_stage_pod/stages/templates/KAVHIFI.conf` and
  `INPH2604.conf`.
- **Read only:** `/tmp/if2453/c67/dl/kavagutti/Support-data.gz` (1,675,816,300 B, the size and sha256
  `3cd7e07da707d65cb454df78632227e34c07f420a1ee006751248f92e4884bb5` that `sources.yaml` records, *measured* by `ls -l`),
  `/stor/scratch/Wilke/if2453/fragdef/db/inphared/7Apr2026/` (`7Apr2026_genomes.fa.gz` 714,521,461 B, sha256
  `331cb0487b5ed60863faab38c6be1496e5ec015a5c0f787dd9ce865212f110c9` per its `MANIFEST.json`/`SHA256SUMS`), the data root
  `/tmp/if2453/r10/root`, the staged inputs of the `r26.11.4` chain under `/tmp/if2453/s2/inputs/`, card C-74's index
  `/tmp/if2453/c74/index/`, and the staged sets of `r26.11.5` stages before these two (`LR2`, `SAG1`, `NX1`, `GV1`,
  `MV2`, `MVEUKV`). Nothing was written to the root, to `/tmp/if2453/s2/inputs/` or to `/tmp/if2453/stages/`; nothing
  under `/stor` was overwritten or deleted.
- **Written:** `/stor/scratch/Wilke/if2453/bio-box-staging/KAVHIFI/` (the two FASTA files and `deferred/euk-MAGs-fna.tar.gz`,
  each over 64 MB); `/tmp/if2453/c80/KAVHIFI/` and `/tmp/if2453/c80/INPH2604/` (metadata, and the 35.5 MB INPH2604 FASTA);
  `/tmp/if2453/c80/{outer_list.txt,bac_list.txt,x/,work/,kav_report.json,inph_report.json}` (listings and scratch). One file
  over 64 MB stayed on `/tmp`: `KAVHIFI/records.tsv` (89,586,620 B, a sidecar no build reads); a copy to `/stor` was refused
  by the session's permission check and not retried.

## 1. Summary

| | `KAVHIFI` | `INPH2604` |
|---|---:|---:|
| Source | Kavagutti et al. 2023, Figshare 20431332: `bac-MAGs-fnas.tar.gz`, `NCLDV-data.tar.gz` | INPHARED 7Apr2026 `7Apr2026_genomes.fa.gz` |
| Records read | 2,214 MAG files (527,434 records) and 3,309 NCLDV files (3,309 records) | 36,025 records |
| Records staged | 530,743 | 386 |
| Genomes listed (DEC-27 repeated sets) | 5,522 (1) | 353 (33) |
| Distinct contigs (bp) | 529,287 (4,659,064,333) | 353 (31,534,280) |
| Kinds | 2,214 `prokaryote` (all code 11), 3,308 `virus` | 353 `virus` (code 11; gv-meta) |
| DEC-65 | NCLDV: 11 exact DTRs trimmed (1,637 bp), 11 circular genomes | none: database records byte-identical (DEC-65 (2)) |
| Held by the root or `r26.11.4` stages | 16 contigs, 138,278 bp (`SPX1`) | 8 contigs, 533,414 bp (`MV1`) |
| Held by `LR2`, `SAG1`, `NX1`, `GV1`, `MV2`, `MVEUKV` (and `KAVHIFI` for `INPH2604`) | 0 | 0 |
| New after all of those | 529,271 contigs, 4,658,926,055 bp | 345 contigs, 31,000,866 bp |
| `s2_check_inputs.py` | exit 0, `problems: []`, `gene_call_mode` `{"prokaryote": null, "virus": "gv-meta"}` | exit 0, `problems: []`, `gene_call_mode` `{"virus": "gv-meta"}` |
| Conf template | `DEPTH=called`, `CRISPR=yes`, `AFTER=MVEUKV` | `DEPTH=annotated`, `LIBRARY=pfam 38.2 gathering`, `CRISPR=yes`, `AFTER=KAVHIFI` |
| Licence (as `sources.yaml` records it) | CC BY 4.0, `redistributable: true` | INSDC records (public), `redistributable: true` |

All *measured* (`c80_kavagutti_inphared.py kav` and `inph`; reports `/tmp/if2453/c80/kav_report.json`,
`inph_report.json`). Against the card's figures (wishlist Z5, Z6, all *measured* there): 2,214 MAGs and 4.607 Gbp
(4,607,298,152 bp) equal; 3,309 NCLDV files (the card's 3,310 counts the directory entry `NCLDV/NCLDV-fna/` of the
listing, 3,310 lines, *measured* by `grep -c`) and 81.36 Mbp of records against 85.39 Mbp, 4.7 % lower
(*derived*: 81,362,268 / 85,390,000); the 10 mitochondrial genomes and 0.274 Mbp (274,163 bp) equal; the 28
eukaryotic MAGs and 560.96 Mbp (560,955,601 bp) equal. INPHARED: 386 records and 34.89 Mbp (34,889,413 bp) against
458 and 39.21 Mbp, 15.7 % and 11.0 % lower (*derived*), inside the 20 % of AGENTS §9; 39 records of 200 kb and more
(9,870,200 bp) against 38, and 11 *Wolbachia* records against 11. The card's 458 is
36,025 - 32,628 - 2,939 (*derived*: the wishlist's arithmetic on three counts), not a set difference; the diff here is
by accession and by sequence (§4), and uses neither of the two disagreeing counts of SAT1's part (plan §4.8).

## 2. The archive

`Support-data.gz` is one tar.gz of `figshare/`: `phages-fna.tar.gz` (taken by `LRV1`, card C-67), `euk-MAGs-fna.tar.gz`
(173,569,286 B), `bac-MAGs-fnas.tar.gz` (1,428,293,652 B), `mitochondrias-fna.tar.gz` (84,999 B), `NCLDV-data.tar.gz`
(26,443,303 B), `bacterial-trees-data.tar.gz` and a `README` (listing, `tar -tzv`, *measured*; its README says "2214
MAGs", "3309 NCLDV contigs", "679 complete phage genomes"). The script reads it once as a stream (`tarfile` in
`r|gz` mode) and takes parts by member name:

| Member | Taken as | Count (*measured*) |
|---|---|---|
| `bac-MAGs-fnas.tar.gz` (`MAGs-fnas/*.fasta`) | one `prokaryote` genome per file, quality `contig` (the archive has no completeness table, as `GV1`), linear (no topology is given) | 2,214 files, 527,434 records, 4,607,298,152 bp |
| `NCLDV-data.tar.gz` (`NCLDV/NCLDV-fna/*.fna`) | one `virus` genome per file (called by gv-meta in any case), DEC-65 (1) exact 20 to 1,000 bp terminal repeat trimmed once and staged circular, as card C-67 did for the phages | 3,309 files and records, 81,362,268 bp; DTR classes exact 11, partial 7, none 3,291 |
| `euk-MAGs-fna.tar.gz` | not staged: card C-60 (eukaryotic calling). Kept byte for byte at `/stor/.../KAVHIFI/deferred/euk-MAGs-fna.tar.gz`, listed in `euk_deferred.tsv` | 28 files, 147,950 records, 560,955,601 bp |
| `mitochondrias-fna.tar.gz` | not staged (§3). Kept at `/tmp/if2453/c80/KAVHIFI/deferred/`, listed in `mito_deferred.tsv` | 10 files (one record each), 274,163 bp |
| `phages-fna.tar.gz`, `bacterial-trees-data.tar.gz`, `NCLDV/NCLDV-tree/` | not taken (the phages are in `LRV1`; the trees and alignments are not genomes) | |

The NCLDV files are named `HRS-<sample>-C<n>`; the mitochondria `HRS-EL3-C80_Synura_synuroidea`, `..._Pleodorina_starrii`,
`..._Teleaulax_amphioxeia`, `..._Copepoda`, `..._Stramenopiles` (names only).

## 3. Genetic codes, the mitochondria and the accession test

- **Codes of the MAGs.** Each MAG is screened with card C-66's `screen_genetic_code` (pyrodigal on the longest contigs
  up to 500 kb, codes 11, 25, 4; no taxonomy hint). All 2,214 stay code 11, reason `code11_pass` (`genetic_codes.tsv`,
  *measured*). The screen would have taken 25 or 4 only where coding density and open-reading-frame length gained.
- **The 10 mitochondrial genomes are deferred, not staged with per-record codes.** The archive's own table gives no
  genetic code, only a file name with the host taxon (the five taxa above: a stramenopile, a chlorophyte, a
  cryptophyte, a copepod, an unclassified stramenopile). Card C-59a's `O2` routes by the `transl_table` of the RefSeq
  record, which these files do not have, and a code chosen here from the host name would be a guess. They are 0.274 Mbp
  (*measured*) and wait for a later organelle stage (`O3`), recorded in `planning/CATALOG_GOAL.md` §4.1.
- **The accession test sees none of the MAGs.** By card C-74 only some resources of GlobDB r232 carry INSDC accessions
  (AMXMAG, GEM, GOMC, SPIRE), and GTDB R232 is indexed in full by accession. The test is `r232_index.lookup`, applied to
  every `GCA_`/`GCF_` accession in a MAG's file name or record headers. **0 of 2,214 MAGs carry one** (file names are
  `HRS-<sample>-bin-<n>`, headers `contig` names; `mag_accessions.tsv`, *measured*): 0 could be looked up in GTDB R232 or
  GlobDB r232, 2,214 could not be seen by the index, and none was dropped by it. No ANI test was made (the card's
  rule). The ENA accessions of the project `PRJEB52406` are not in the archive, and fetching them is outside this card
  (AGENTS §8: no entry lists the ENA assembly records), so whether GTDB R232 holds these MAGs under an INSDC accession
  stays unknown here.

## 4. The INPHARED diff

`7Apr2026_genomes.fa.gz` has 36,025 records (the `.fai` and the data table agree, *measured*). Known accessions: the
versionless first header token of every record the root's `contig_label` holds for stage `S1`, and of every row of
`SAT1`'s `records.tsv` (69,066 accessions, *measured*); known sequences: the `contig_uid` of the root's `S1` default
variants and of `SAT1`'s `genomes.json`. A record is staged when its accession is not known **and** its `contig_uid` is
held by neither. The four cells over the 36,025 records (*measured*):

| | sequence held by `S1` or `SAT1` | sequence not held |
|---|---:|---:|
| accession known | 29,247 | 131 |
| accession not known | 6,261 | **386** |

The 6,261 "accession not known, sequence held" are the same sequences under another accession (RefSeq `NC_` copies of
GenBank records, renamed records); the 131 "accession known, sequence not held" are records whose sequence changed.
Neither group is staged; both crossed cells are listed in `skipped_records.tsv` (6,392 rows). The 386 staged records are
byte-identical to the archive (header, line width and sequence). Of them, 33 repeat another staged record's sequence and
are listed once as a genome (DEC-27; `repeated_records.tsv`), so 353 genomes and 353 distinct contigs. 39 are of 200 kb and
more (9,870,200 bp), 11 name *Wolbachia* in their description, 0 name a synthetic construct (*measured*).

## 5. Staged files

| File | Where | Bytes | sha256 (*measured*, by the script) |
|---|---|---:|---|
| `kavagutti_hifi_bacterial_mags.fna` | `/stor/scratch/Wilke/if2453/bio-box-staging/KAVHIFI/` | 4,692,735,019 | `5c310ee70e183ecf4a5132fd85a2d738824a93d521e3d75bc38471bc3942dcd5` |
| `kavagutti_hifi_ncldv.fna` | same | 82,769,235 | `02275a2997b6dee6989cd0df16f7a7bbc3979fba8fcadc2798fe9232eaa4c6c1` |
| `deferred/euk-MAGs-fna.tar.gz` | same (the 28 eukaryotic MAGs, for card C-60) | 173,569,286 | `f9ea899e78b1751cd4e49f693a8bf3f3d27e527fe194bcd630676767d0ecd8fe` |
| `genomes.json` | `/tmp/if2453/c80/KAVHIFI/` | 19,681,400 | `1efb3b86fc62f5cdf18b397cf561b0ddea0c029c61beba37739650afd7e5ef6c` |
| `objects.json` | same | 1,005 | `b3461aaad23eb1d380f7135045b75ae5555343939cc38a10617a502edaf50dcf` |
| `objects_staging.json` | same | 1,047 | `a88411f2ece88b7b716ac8f0a2ea3528a934165c2fdd2a242df3680adc7dfdd3` |
| `records.tsv` | same (stays on `/tmp`, §header) | 89,586,620 | `f14164c4c782935ae4015352f579df675a9b351e65a71310670c8fed4a61c344` |
| `repeated_records.tsv` | same | 107 | `9297a722525c617a47ac5b678bdcdc86e4baa7cc446092248719375b8eca25be` |
| `shared_records.tsv` | same | 1,799 | `2e7a876e226b7507e562c595339e8ad928f41077d0bd84d905816e75ac6dbe89` |
| `later_records.tsv` | same | 33 | `887d06e335810fda96b778fac18eb04f5d4b4e39cb47c0a3b9573d32f79fc426` |
| `trims.tsv` | same | 1,526 | `af5ab7defce1a00fa07245d3d6f0a50bdc2632d73c8cda6c5533362dd99b84e8` |
| `genetic_codes.tsv` | same | 105,984 | `8459efda11d94d03312d6cef44170a76121f6d936ac5701b28e103ebab6fea94` |
| `mag_accessions.tsv` | same | 41,467 | `54c869e1fdc7d624264f1fd6f96d0da7691c1760ce3f96ef5deb50971f892d2b` |
| `skipped_members.tsv` | same (header only: no member skipped) | 30 | `eda3234c939198ea38d06d137bd43906a0a0bd23d7684629519d3f727ac34bca` |
| `euk_deferred.tsv` | same | 2,668 | `139fffe83248049dc17e844958d8b1353c082bdab95b966354afbd42162f34dc` |
| `mito_deferred.tsv` | same | 490 | `ec6009a39928aca544c57885913adb2f8e0baa240d9f2e8775b9938fc4f41b7c` |
| `deferred/mitochondrias-fna.tar.gz` | same | 84,999 | `4bc99a12ca14f84e39412853178ad510b348ff28f8c019ec1d19198a5db6b76a` |
| `inphared_7apr2026_delta.fna` | `/tmp/if2453/c80/INPH2604/` (35.5 MB, under 64 MB) | 35,496,233 | `de2c3531ac6930a3da2c4f77f0fdb333e8fb16ec970f6b0089efd910d4aa1c15` |
| `genomes.json` | same | 49,420 | `cacb54f32459b70880d3bdbfe80018ca59cd5891b063f228d25b874265ff06b0` |
| `objects.json` | same | 476 | `ce4e5007a13d32a79048b305076380dbd3beb653334d302a551ae5c982be1497` |
| `objects_staging.json` | same | 470 | `376ed369ab84cf3dc15bff2840cc8a0a8a6ba3c8656dddbe2331c46fc410faf7` |
| `records.tsv` | same | 57,088 | `0086cc0351bc4706a0d0b7436bc9bb4271eabcac867d77ca87390a20fcdf96da` |
| `repeated_records.tsv` | same | 1,836 | `821ba96de0c8020aa221a51b1bcba4a3d68afb5cd85b484efe14973e89851fdf` |
| `shared_records.tsv` | same | 439 | `6e7737689509af0b377d1fad2a76d6d5f920576150de61a40c71090de652b920` |
| `later_records.tsv` | same | 33 | `887d06e335810fda96b778fac18eb04f5d4b4e39cb47c0a3b9573d32f79fc426` |
| `skipped_records.tsv` | same | 576,104 | `25b28759c24986d82cdf99531060e0e432250f33abd54d1155722b7fd029cc00` |
| `trims.tsv` | same | 102 | `3b6006a823ad5d30eb7116c77f3fd53d5430acfb2e9ca0df4a7dd97a8baab1dd` |

`objects.json` of each stage names `/tmp/if2453/s2/inputs/<STAGE>/<file>` (card A-29 copies there);
`objects_staging.json` names the file where it is now. The conf templates' `SUM` lines are the rows for `genomes.json`,
`objects.json` and the FASTA files. The FASTA of `KAVHIFI` is in the archive's order (MAG files, then NCLDV files);
`genomes.json` and `records.tsv` are sorted by category and genome name.

## 6. Overlap (exact `contig_uid`)

By build order (*measured*, `kav` and `inph`): the root's 13 stored stages `S1, P1..P6, H1, PL1, C5, V1, G1, SPX1` from
`contig_index` with `is_default_variant`; then the 19 `r26.11.4` stages from their `genomes.json` under
`/tmp/if2453/s2/inputs/` in `a26_r26_11_4_chain.R26_11_4_CHAIN` order; the first holder owns. These are the templates'
`EXPECT=held_by` lines.

| Holder | `KAVHIFI` contigs | bp | `INPH2604` contigs | bp |
|---|---:|---:|---:|---:|
| `SPX1` | 16 (29 records in 12 MAGs) | 138,278 | 0 | 0 |
| `MV1` | 0 | 0 | 8 | 533,414 |
| every other stage of the root and `r26.11.4` | 0 | 0 | 0 | 0 |

Later in the `r26.11.5` chain (plan §5), for card A-29: `GV1`, `LR2`, `SAG1`, `NX1`, `MV2`, `MVEUKV` hold 0 contigs of
`KAVHIFI` (so the NCLDV files and MAGs do not repeat `GV1` or `MVEUKV`'s giant viruses byte for byte) and 0 of `INPH2604`;
`KAVHIFI` holds 0 of `INPH2604`. The sets of the staging cards C-75 to C-78 that are not yet on the pod
(`GRARE1..3`, `EGV1`, `ECE2`, `MIR1`, `ASG1`, `SAG2`, `DPN1`, `SYM2`, `LRP1`, `MIR2`, `MGE4`) are card A-29's to test by
chain order. Of `KAVHIFI`'s genomes, 774 contigs are listed in more than one genome (DEC-27; MAGs sharing a contig
sequence), and one NCLDV file repeats another's sequence (`HRS-ES2-C2780` is `HRS-ES1-C3014`, `repeated_records.tsv`).

## 7. The remainder (DEC-59)

Recorded in `planning/CATALOG_GOAL.md` §4.1: the 28 eukaryotic MAGs (147,950 records, 560,955,601 bp, `euk-MAGs-fna.tar.gz`,
sha256 `f9ea899e...ecd8fe`) for card C-60; the 10 mitochondrial genomes (274,163 bp) for a later organelle stage; the
INPHARED 7Apr2026 records not staged (35,508 whose sequence `S1` or `SAT1` holds, 131 whose accession is known and sequence
changed; `skipped_records.tsv` lists the 6,392 of the two crossed cells; the 29,247 held under a known accession are in
`S1` or `SAT1`); the Kavagutti phages are in `LRV1` (card C-67).

## 8. Build cost (*derived*, the plan's §3 rates; lower bounds, unpriced steps excluded)

- `KAVHIFI` (`called`): MAGs, class R1, 24.8 CPU-minutes per new Gbp, × 4.439 new Gbp (*derived*: 4,577,725,195 distinct
  MAG bp less 138,278 held by `SPX1`) = 110.1; NCLDV, class R2, 93.5 to 135.3 per new Gbp, × 0.0813 Gbp = 7.6 to 11.0;
  together about **118 to 121 CPU-minutes** (plan §4.8: about 125).
- `INPH2604` (`annotated`): R2 + R3, (93.5 + 2,090) to (135.3 + 2,740) per new Gbp, × 0.0310 new Gbp = about **68 to 89
  CPU-minutes** (plan §4.8: 85 to 110, for 39.2 Mbp).

## 9. Pod compute used

All in detached `systemd-run --user` units at `nice -n 19`, `MemoryMax=16G`, `MemorySwapMax=0` (the staging runs inside
`/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif`, `PYTHONPATH` the tree's `src`), after looking at the active units
(the build `r26-11-4-wave-4`/`v2-4`, and one staging unit of card C-76 at a time) and checking that no `STOP` file existed:

| Unit | What | Cap | Wall (*measured*, file times) | CPU |
|---|---|---|---|---|
| `c80-list`, `c80-x1`, `c80-x2` | listings and small extractions of the archive | `CPUQuota=100%` | about 1 to 3 min each | at most 8 CPU-minutes together (*derived*: wall x 1) |
| `c80-kav` | one stream of `Support-data.gz`, screen of 2,214 MAGs (2 workers), ownership, `s2_check_inputs.py` | `CPUQuota=400%` | 14:12 to about 14:21 (9 min) | at most 36 CPU-minutes (*derived*: wall x 4); 453 CPU-seconds used after 8 minutes (*measured*, `systemctl show`) |
| `c80-inph` | one stream of `7Apr2026_genomes.fa.gz`, diff, ownership, `s2_check_inputs.py` | `CPUQuota=200%` | 14:22 to about 14:24 (2 min) | at most 4 CPU-minutes (*derived*) |

No process came near 300 CPU-minutes; no more than three processes of the unit ran at once. The card's estimate of under 1 hour at
4 threads (*estimate*) held.

## 10. Tests

`tests/r1/test_c80_kavagutti_inphared.py`: constructed truth written out in the file (a fake archive with the six members;
the diff over eight records and the four cells; contig_uid by blake3 directly): the parts taken by member name; kinds,
codes (the screen stubbed for code 25) and DEC-65; the deferral of the eukaryotic MAGs and the mitochondria, byte for
byte; the accession test through the stubbed lookup and through `r232_index.lookup` over a hand-written index; ownership
lists; byte identity of database records; twice-run determinism; `s2_check_inputs.py` on both stages; the refusals
(archive without the MAG member, nothing new); the default screen; the templates through `stage_est.read_conf`, their
`SUM`s found in this note. Each test runs with the normal `PATH` and with `PATH=/usr/bin:/bin`. Seen failing first on the
pod (26 failed: no script), then passing (*measured*, wilkcomp02, 2026-11-26, `scripts/pod_pytest.sh C-80
tests/r1/test_c80_kavagutti_inphared.py tests/r1/test_r232_index.py tests/r0`): **270 passed in 34.45 s** (each C-80 test
twice, normal and tools-hidden `PATH`). `stage_est.py conf` on both templates (`python -I`, `--tree` the mirrored tree,
in the image): exit 0 for both.

## 11. Mirror to `tower`

`scripts/orchestrator/relay_to_tower.sh` (card C-49), from the laptop, one stream, three runs, to
`tower:/mnt/user/research-bulk/bio-box/sources/KAVHIFI/` and `.../INPH2604/` (*measured*, 2026-11-26 14:34 to 14:36, the
script's log): `KAVHIFI`'s two FASTA files and `deferred/euk-MAGs-fna.tar.gz` (4,949,073,540 B, `--dry-run` first), then its 14
metadata files (with `records.tsv` and `deferred/mitochondrias-fna.tar.gz`), then `INPH2604`'s 10 files; every file's sha256
on `tower` equal to §5's; `DONE: 3, 14, 10 verified, 0 failed`; `SHA256SUMS` written in each directory (17 lines in
`KAVHIFI`, 10 in `INPH2604`). Both stay `redistributable: true` as their entries record; the mirror is the PI's own NAS.
