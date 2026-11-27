# The coral-reef sets: Great Barrier Reef genomes (`GBR_MGD`, stage `GBRMGD`) and Tara Pacific's Reef Microbiomics Database (`RMD1`): relay, contents, staging (card C-83, DEC-85)

- **Date and host:** 2026-11-27, `wilkcomp02` (the pod) and `tower` (read only, then the mirror). Every number carries
  a label (AGENTS §6): *measured* (on `wilkcomp02` that day, by the command or file named), *derived* (formula and
  inputs), *estimate*, or *published*.
- **Script:** `scripts/stress/build_stage_pod/c83_coral_inputs.py` (`census`, `gbrmgd-stage`, `rmd1-plan`,
  `rmd1-stage`, `ani`, `overlap`); tests `tests/r1/test_c83_coral_inputs.py`; templates
  `scripts/stress/build_stage_pod/stages/templates/GBRMGD.conf` and `RMD1.conf`.
- **Read only:** the relayed files (§1), card C-74's index `/tmp/if2453/c74/index/` and GTDB R232's metadata
  `/tmp/if2453/c74/GTDB_R232_META/`, the data root `/tmp/if2453/r10/root`, `/tmp/if2453/s2/inputs/` (the `r26.11.4`
  stages), and the staged metadata of the `r26.11.5` stages on the pod (§2.2). Nothing was written to the root or to
  `/tmp/if2453/s2/inputs/`; nothing under `/stor` was overwritten or deleted.
- **Written:** the archives (`dl/`) and the staged FASTA to `/stor/scratch/Wilke/if2453/bio-box-staging/{GBRMGD,RMD1}/`;
  metadata to `/tmp/if2453/c83/meta/{GBRMGD,RMD1}/`; census tables, reports and ANI tables to `/tmp/if2453/c83/work/`.

## 0. What each archive holds (read before anything was staged)

`c83_coral_inputs.py census` reads an archive once, streaming, and writes one row per member: its compression told by
its first bytes (never by its name), and its content told by its decompressed bytes (records, bases, nucleotide
letters, gap bytes, other bytes); a member is `dna` only when every sequence byte is an IUPAC nucleotide letter
(*measured*, `/tmp/if2453/c83/work/{GBRMGD,RMD1}/members.tsv` and `census.json`).

| Archive | Members | By kind | Genomes by kind | Staged | Left out |
|---|---:|---|---|---|---|
| `prokaryote_bins_5283.tar.bz2` (`GBR_MGD`) | 5,283 files and one directory `prokaryote_bins/` | 5,283 `dna` (plain FASTA, none compressed inside), 467,187 records, 10,244,764,294 bp; no protein, alignment, table or other file | 5,283 prokaryote metagenome bins (MAGs: the file names are binners' names, `<site>__metabat2_refined_<n>`, `__rosella_refine_<n>`, `__vamb_<n>`; the headers say `[group=P]`); no eukaryote, virus, isolate or single cell | all 5,283 | nothing |
| `RMDv1-genomes.tar.gz` (`RMD1`) | 13,549 files and one directory `RMDv1-genomes/` | 13,549 `dna` (plain FASTA, none compressed inside), 4,850,322 records, 40,388,129,994 bp; no protein, alignment, table or other file | 13,549 prokaryote genomes, each a row of `genomes-aggregated-summary.tsv` (GTDB-Tk R207 and CheckM): 1,524 Tara Pacific coral MAGs (`TPAC-CORAL-MAGS`), 11,922 MAGs of other reef studies (`EXTERNAL-MAGS`), 103 isolate genomes (`EXTERNAL-ISOG`); no eukaryote (every member is a row of the prokaryote table), no virus, no single cell | all 13,549 planned; **staging not finished at this commit** (§3) | nothing |

The census table of `GBR_MGD` has sha256 `857de8909e7d46725c2e15b40505ac71a9a02f258f8cc17d4161203451dd0dc6`
(1,320,530 B). Its headers are the assembly's contig names with the bin's name in brackets
(`IMOS__P__hoskyn_2__contig_3485_pilon:p [group=P] [bin:hoskyn_2__rosella_refine_806_0] [length=622590bp]`, some with
`[proviral_region:virsorter=...]`): Pilon-polished assemblies of reef-water samples. The archive holds **no table**:
no completeness, contamination, taxonomy, sample or accession. The Zenodo record may hold other files (eukaryote bins,
viruses, tables); only this one was listed and fetched (`sources.yaml` `tower_gbr_mgd`), and nothing else of the
record was asked for here (AGENTS §8).

## 1. The bytes on the pod

`scripts/orchestrator/relay_to_pod.sh` (card C-73) from the laptop, `--dry-run` first, each list `tower`'s
`SHA256SUMS` taken to the laptop with `scp` (*measured*, the script's output, 2026-11-27 CST; every sha256 equal to
`tower`'s record):

| File | Bytes | Where on the pod | Time | sha256 |
|---|---:|---|---|---|
| `prokaryote_bins_5283.tar.bz2` | 2,967,595,725 | `/stor/scratch/Wilke/if2453/bio-box-staging/GBRMGD/dl/` | 103 s (28.81 MB/s), ending 07:34 | `000d87436b9a408ac7a09b5e70c9d92e42d30312ead648b1d2e55dddca56bea0` |
| `RMDv1-genomes.tar.gz` | 11,571,956,924 | `/stor/scratch/Wilke/if2453/bio-box-staging/RMD1/dl/` | 466 s (24.83 MB/s), ending 07:43 | `aa4b6631e5c6de2bf8e0cd8828170bcb3d0138d545f749700396ebad8a6fa21a` |
| the 14 tables of `tower_rmd1` (`genomes-aggregated-summary.tsv` 4,864,682 B, `RMD-motus-membership.tsv`, `supplementary-table-1..12.xlsx`) | 15,133,291 | `/tmp/if2453/c83/dl/RMD1/` (with `SHA256SUMS`, 14 lines) | 07:35 | as `tower`'s `batch6_symbionts/RMD1/SHA256SUMS` (`genomes-aggregated-summary.tsv` `8af9e1806be793d2d83a4277ef4357425f8e7ebb87b8feaab2ceabfcc8c44dd6`) |

The relays ran one after another; `tower`'s fetch chain was running and nothing else was started there.

## 2. `GBRMGD`

### 2.1 Staging (`c83_coral_inputs.py gbrmgd-stage`, report `/tmp/if2453/c83/work/gbrmgd_report.json`)

Rule: every `dna` member of `prokaryote_bins/` is a genome, its id the file's stem; the bins in id order for the
metadata, the archive's order for the FASTA; whole contigs; every contig name given the bin's id and `|` (all 467,187:
none begins with its bin's id); a header that is not UTF-8, an empty record or a byte that is not a nucleotide letter
refuses the stage (none did); `kind="prokaryote"`; the genetic code screened per bin (`c58.screen_genetic_code`, 11,
25, 4) with no taxonomic hint, since the archive gives none; `assembly_quality` `mag_mq` for every bin, because no
table of completeness was fetched (a label, not a measurement of the bins' quality); one BGZF FASTA.

| | Value (*measured*) |
|---|---:|
| Bins / FASTA records | 5,283 / 467,187 |
| Genomes listed (DEC-27 repeated contig sets) | 5,283 (0) |
| Distinct contigs (bp) | 466,872 (10,238,417,378); 125 contigs in more than one bin (bins of different binners on one assembly; `s2_check_inputs.py` accepts them, DEC-27) |
| Genetic codes | code 11: 5,283; conflicts 0 |
| r232 status | `not_visible`: 5,283 (10,244,764,294 bp): no bin carries an accession |
| Held by the root or `r26.11.4` (exact `contig_uid`, build order) | 0 |
| Held by the `r26.11.5` stages staged before it (§2.2) | 0 |
| **New after both** | **466,872 contigs, 10,238,417,378 bp** |
| `s2_check_inputs.py GBRMGD` (with `objects_staging.json`, the byte pass included) | exit 0, `problems: []`, 5,283 genomes, 467,187 `genome_input` rows, codes `[11]`, kinds `["prokaryote"]` |

Depth: 10.238 Gbp new is well over the card's 1.5 Gbp, so `DEPTH=called`, `CRISPR=yes`, default gene-call mode;
`AFTER=RMD1` (DEC-85 (1)). No contig is held by the root or `r26.11.4`, so the template has no `EXPECT` line.
`redistributable: false`: the record states no licence (the job list copied its licence field as `{}`; DEC-80 (2),
(3)). `stage_est.py conf` on the template, in the image: exit 0 (*measured*). Against the card's description
(5,283 prokaryote MAGs, *published* in the `sources.yaml` entry): equal.

### 2.2 The owners and later holders read

Owners in build order (the method of cards C-77, C-79 and C-81): the root's 13 stored stages `S1, P1..P6, H1, PL1, C5,
V1, G1, SPX1` (`contig_index`, default variants), then the 19 `r26.11.4` stages from their `genomes.json` under
`/tmp/if2453/s2/inputs/` in `a26_r26_11_4_chain.R26_11_4_CHAIN` order. Later holders, the `r26.11.5` stages staged on
the pod, in plan §5's order: `GRARE1..3` (`/tmp/if2453/c75/meta/`), `LR2` (`/tmp/if2453/r19/LR2/`, card R-19's
re-stage), `SAG1` (`/tmp/if2453/c58/SAG1/`), `NX1` (`/tmp/if2453/r20/NX1/`, card R-20's), `GV1`
(`/tmp/if2453/r23/GV1/`, card R-23's), `EGV1`, `ECE2`, `MIR1`, `SAG2`, `DPN1`, `SYM2`, `LRP1`, `MIR2`, `MV2`,
`MVEUKV`, `KAVHIFI`, `INPH2604`, `GWV1`, `SAG3`, `WET1` (their cards' metadata directories). `ASG1` and `MGE4` are out
of `r26.11.5` (DEC-85 (3)).

### 2.3 ANI against the staged prokaryote sets of `r26.11.5` and the r226 archaea (recorded, not cut)

`c83_coral_inputs.py ani` (*measured*, `/tmp/if2453/c83/work/gbrmgd_ani/summary.json`, `ani_hits.tsv` and one skani
table per group, unit `c83-gbr-ani`): one FASTA per bin as queries; each reference group in its own `skani dist -t 8
--min-af 50` process (default sketch); a bin meets when ANI >= 96 over an aligned fraction >= 50 % (the larger of the
two, card C-66's rule); the best hit over groups is kept. The bins carry no taxonomy, so every group is whole: `LR2`
by its MAGs (`assembly_quality` `mag_hq` or `mag_mq`, not its unbinned contigs), `SYM2` and `KAVHIFI` by their
bacterial-bin FASTA, and all 14,839 genomes of card C-75's r226 archaea.

| Reference group | Genomes | Bins meeting (in the group) | Best hit there | skani CPU-s |
|---|---:|---:|---:|---:|
| `GRARE1` | 4,779 | 29 | 8 | 128.3 |
| `GRARE2` | 3,806 | 2 | 2 | 120.3 |
| `GRARE3` | 3,638 | 58 | 57 | 122.9 |
| `LR2` (MAGs) | 5,159 | 0 | 0 | 168.4 |
| `SAG1` | 5,537 | 1,604 | 1,483 | 182.5 |
| `NX1` | 893 | 0 | 0 | 73.0 |
| `SAG2` | 9,698 | 105 | 18 | 106.6 |
| `DPN1` | 1,840 | 3 | 0 | 66.8 |
| `SYM2` (bins) | 625 | 0 | 0 | 58.3 |
| `KAVHIFI` (bacterial MAGs) | 2,214 | 0 | 0 | 84.3 |
| `SAG3` | 3,568 | 548 | 92 | 94.4 |
| `WET1` | 4,164 | 0 | 0 | 128.2 |
| r226 archaea (card C-75) | 14,839 | 541 | 530 | 295.2 |
| **any** | | | **2,190 of 5,283 (41.5 %)** | 1,629.2 |

All *measured*; the percentage *derived*. **Nothing is cut by it.** Reading: the reef-water bins are largely of
species that the surface-ocean single cells of `SAG1` (GORG-Tropics) and the marine archaea of GlobDB r226 already
hold, which is what open-water bins of the Great Barrier Reef would be; none is byte-identical to a held contig
(§2.1). `RMD1` as a reference group for `GBRMGD` is not yet run (§3).

## 3. `RMD1` (not finished at this commit)

**Done:** the relay (§1); the census (§0: 13,549 DNA members, 4,850,322 records, **40,388,129,994 bp**, *measured*;
inside the card's 10 to 60 Gbp, so its stop is not in play; the session's estimate was 30 to 40 Gbp: 1 % over its upper
end); the plan (`rmd1-plan`, *measured*, `/tmp/if2453/c83/work/rmd1_plan.json`), below.

| Plan of `RMD1` | Genomes | bp |
|---|---:|---:|
| MAGs (`TPAC-CORAL-MAGS` 1,524, `EXTERNAL-MAGS` 11,922), quality from CheckM: `mag_hq` (over 90 % complete, under 5 % contaminated) 6,319, `mag_mq` 7,127 | 13,446 | 39,907,128,492 |
| isolates (`EXTERNAL-ISOG`), `assembly_quality` `scaffold` (the table gives no assembly level) | 103 | 481,001,502 |
| r232 status: isolates by their BioSample through GTDB R232's `ncbi_biosample`, then C-74's `lookup`: `in_gtdb_and_globdb_r232` 52, `in_gtdb_r232` 51 (all 103 found) | 103 | |
| r232 status: MAGs `not_visible` (a MAG's BioSample in the table is its metagenome's, shared by every MAG of that sample; no GTDB R232 genome carries any of these BioSamples) | 13,446 | |
| table rows not in the archive / members not in the table / left out | 0 / 0 / 0 | |

Rule of the stage (as `GBRMGD`'s, §2.1): `kind="prokaryote"`, the code screened per genome with GTDB-Tk's lineage as
the hint, contig names (`<genome>-scaffold_<n>`) given the genome's id and `|`, one BGZF FASTA
`rmd1_reef_genomes.fna.gz`; `redistributable: false` (the tarball's job list says `open-academic`, the worker's label,
not a licence; the tables' Zenodo record states `cc-by-4.0`; DEC-80 (2), (3)). The screen workers are replaced after
500 genomes each (`max_tasks_per_child`), so that no process passes 300 CPU-minutes.

**Not done at this commit:** the staging itself (unit `c83-rmd1-stage`, started 2026-11-27 09:17:19 CST, still running at this commit, writing
`/stor/scratch/Wilke/if2453/bio-box-staging/RMD1/rmd1_reef_genomes.fna.gz` and `/tmp/if2453/c83/meta/RMD1/`, report
and `s2_check_inputs.py` result in `/tmp/if2453/c83/work/rmd1_report.json`); its overlap counts; its ANI step and
`GBRMGD` against `RMD1` (`/tmp/if2453/c83/bin/rmd1_ani.sh`, not started); `RMD1.conf`; `stage_est.py conf`; the move
of its `records.tsv` (over 64 MB) to `/stor`; the mirror to `tower`; the `sources.yaml` and `CATALOG_GOAL.md` lines
of `RMD1`.

## 4. Staged files

| File | Where | Bytes | sha256 |
|---|---|---:|---|
| `gbrmgd_bins.fna.gz` | `/stor/scratch/Wilke/if2453/bio-box-staging/GBRMGD/` | 3,189,105,489 | `46e74ae0c87b6e529afe684b1560bed61364de5d200cde232ea1eba1f7e75a7d` |
| `genomes.json` | `/tmp/if2453/c83/meta/GBRMGD/` | 17,384,013 | `42485e2ba213ff648c860d410f58bde007b75a24c94eeaf5d4bf2b6f696a2c4f` |
| `objects.json` | same | 444 | `b58627809f7e6f6dd55b783710efd5ef67baf8f50b2592e8bf2f28c43185d7b9` |
| `objects_staging.json` | same | 465 | `e140a51f91ccb5e1e2fbfabd3bb746b84c5cd4fe164d0067ad2b1eba7184f343` |
| `records.tsv` | same | 45,518,147 | `384b9160c220caa10e892ac41331630fb878e8aacd973ba627816876ed233e4e` |
| `genetic_codes.tsv` | same | 559,452 | `d1998710f73156e2c8c07a89b89fd09cfd8dd72d465464b9bdf563adf3bd7904` |
| `r232_status.tsv` | same | 430,883 | `eb08faf27085e9e250ceb81213f6cdcde2b18a39e03ea36fcb775862cd2b76e2` |
| `repeated_records.tsv` | same | 39 | `d109495ee76e4643c181eae183c452c9a84560f97a79bbe22d6630780368c565` |
| `shared_records.tsv` | same | 60 | `45d2b7bc910a65fe91a4e3623653345589dae4ffa7a7a4499641e13d5cea1b06` |
| `later_records.tsv` | same | 55 | `ba2f1dbd71ed4260b0ec24d03d968ba0c6251afbeeb9efb134735ed1064d05f4` |
| `left_out.tsv` | same | 33 | `4ad9c0ef6fa328dd54a1a79b7720ec0e031daaa22c48567e1490df6d84a07efc` |

All *measured* (`ids.sha256_file`, the stage reports). `objects.json` names `/tmp/if2453/s2/inputs/<STAGE>/<FASTA>`
(card A-29 copies there); `objects_staging.json` the file on `/stor`. The templates' `SUM` lines are the rows for
`genomes.json`, `objects.json` and the FASTA.

**Mirror to `tower`** (`scripts/orchestrator/relay_to_tower.sh`, card C-49, from the laptop, one stream, *measured*,
the script's log, 2026-11-27): `GBRMGD`: the FASTA, 3,189,105,489 B in 131 s (24.34 MB/s), `DONE: 1 verified`, then
the ten metadata files, `DONE: 10 verified, 0 in place and skipped, 0 failed`, `SHA256SUMS` on `tower` 11 lines; to
`tower:/mnt/user/research-bulk/bio-box/sources/GBRMGD/`.

## 5. Pod compute used

Detached `systemd-run --user` units at `nice -n 19`, `MemorySwapMax=0`, inside
`/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif`, each started after a look at the active units (one unit of
this card at a time while two of cards A-29 and B-75 were active; never more than two of this card's), with no
`STOP` file in the root and no build unit of `MV1` or `V3` calling genes (the wave was in `PL3`'s Pfam scan):

| Unit | What | Cap | Wall (*measured*, unit and file times) | CPU |
|---|---|---|---|---|
| `c83-gbr-census` | the first census, with card R-19's line reader: stopped by hand after 4 minutes at about 1.35 MB/s of the archive (*measured*, `/proc/<pid>/io`), for the block reader of §0 | `CPUQuota=100%`, 8G | 07:40:16 to about 07:44 | at most 4 CPU-minutes (*derived*: wall × 1; 1.18 at 1 min 10 s, *measured*, `systemctl status`) |
| `c83-gbr-census2` | `census` of `GBR_MGD` | `CPUQuota=100%`, 8G | 07:44:35 to 07:56:09 | at most 11.6 CPU-minutes (*derived*: wall × 1) |
| `c83-gbr-stage` | `gbrmgd-stage` (3 code-screen workers, DuckDB 2 threads) | `CPUQuota=400%`, 32G | 07:56:28 to 09:01:40 | at most 261 CPU-minutes over several processes (*derived*: 65.2 × 4) |
| `c83-est-gbr` | `stage_est.py conf GBRMGD` | `CPUQuota=100%` | 0.6 s | 0.6 CPU-seconds (*measured*) |
| `c83-rmd1-census` | `census` of `RMD1` | `CPUQuota=100%`, 8G | 09:03:19, 7 min 56 s | 7.8 CPU-minutes, 337 MB peak RSS (*measured*, `/usr/bin/time -v`) |
| `c83-gbr-ani` | `ani` of `GBRMGD`, 13 skani processes of 8 threads | `CPUQuota=800%`, 32G | 49 min 38 s | 62.7 CPU-minutes (*measured*, `/usr/bin/time -v`), largest skani process 4.9 |
| `c83-rmd1-plan` | `rmd1-plan` | `CPUQuota=100%`, 16G | under 10 minutes | at most 10 CPU-minutes (*derived*) |
| `c83-rmd1-stage` | `rmd1-stage` (3 recycled code-screen workers, DuckDB 2 threads); **running at this commit** | `CPUQuota=400%`, 32G | from 09:17:19 | 241 CPU-minutes in all processes at the handover (*measured*, `CPUUsageNSec`) |

No process came near 300 CPU-minutes.

## 6. Tests

`tests/r1/test_c83_coral_inputs.py`, constructed truth throughout: contig uids by blake3 in the test; fixture archives
written here (`tar.bz2` and `tar.gz`) whose members' compression is hidden behind names that say otherwise (a gzip
member named `.fa`, a plain one named `.fa.gz`), with a bzip2 and an xz member, a protein FASTA, an alignment with gap
bytes, a table and an empty file; the block reader against `biobox.fasta.iter_records` at block sizes of 5, 7, 64
bytes and 4 MiB, with CRLF line ends and a last line without its newline; the record reader's refusals (a gap byte,
another byte, an empty record, sequence before a header, a header that is not UTF-8); the GBRMGD plan and stage end to
end twice, byte-identical, with a bin that repeats another's contig set (DEC-27), an owner and a later holder, names
prefixed, the FASTA in archive order, `s2_check_inputs.py` exit 0; a stage refused when a member's bytes are not what
the census read, leaving nothing; `genomes_where`, the r226 archaea split, the overlap count, and skani itself on 60 kb
fixtures (skipped without skani); the RMD1 table and stage (an isolate found in GTDB R232 by its BioSample, one not
found, a MAG not visible with the count of R232 genomes of its BioSample, a *Mycoplasma* MAG screened to code 4);
and the templates through `stage_est.read_conf` with their `SUM`s found in this note. Every test runs with the
normal `PATH` and with `PATH=/usr/bin:/bin`.

Seen failing first on the pod (*measured*, `wilkcomp02`, 2026-11-27, `env BIOBOX_TREE=/tmp/biobox-wt/C-83 sh
scripts/pod_pytest.sh C-83 tests/r1/test_c83_coral_inputs.py tests/r1/test_r232_index.py tests/r0`, with the script
not yet in the tree): **8 failed, 244 passed** (`FileNotFoundError: ... c83_coral_inputs.py`), the census tests; the
template test, written before the template held its values: **4 failed** (`not KEY=VALUE: '@SUMS@'`, `no template
RMD1.conf`). The staging, block-reader and ANI tests were written after the code they test (a deviation from AGENTS
§5.1, reported). With `GBRMGD` staged (same command, one process): **283 passed, 1 skipped in 27.33 s** (the skip: the
skani test under the tools-hidden `PATH`).
