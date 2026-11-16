# S2C phase 2: the inputs of stages C1 to C6, and what they share (card C-36)

Date: 2026-11-15, wilkcomp02, 20:23 to 20:54 CST (the pod's clock). Every number carries a label: *measured* (host,
date, command or file), *derived* (formula and inputs) or *estimate*. **No build**: nothing was written in the data
root `/tmp/if2453/r10/root` (bound read-only into every command that read it) and nothing under `/stor` (the two
archives were read in place). Card C-26 part 1 (`2026-11-15_plsdb_inputs.md`) is the model for this note.

Small results are committed beside this file in [`2026-11-15_s2c_phase2_inputs/`](2026-11-15_s2c_phase2_inputs/):
`stages/` (`assign.json`, `stage_table.tsv`, `stages.tsv` with every genome's stage, each stage's `verify.json`),
`inputs/` (each stage's `make_inputs.py` report `C<k>.json`, its dry check `C<k>.check.json` with the object list
summarised, `SUMS`), `analyze/` (`summary.json`, and per stage `C<k>_contigs_elsewhere.tsv`: every contig another
stage stores or lists, or that several genomes of the stage list), `rehearsal/` (the sample's list and report, the
three units' logs, `check.json`), and the pod's `driver.log`.

## Summary

- **Stages reproduced exactly.** The proposal's rule, written again in `c36_stages.py` from the fetch note's text
  (not copied from `c25_post.py`), gives a per-stage table **byte-identical** to `stage_proposal.tsv` (sha256
  `45a1fd00…e756` both; *measured*, `stages/assign.json` `equal_to_proposal: true`): 13,477 genomes, 56,177,990,736
  bp, n = 6. The fetched set (the two archives' `kept.tsv`) is exactly the remainder less the 8 not served.
- **Inputs** in `/tmp/if2453/s2/inputs/C1/` to `C6/` on the pod, made by `make_inputs.py prok` unchanged, one stage
  lock and selection each. Every genome's `.fna.gz` matched its lock's `gz_sha256` and size, and its FASTA the lock's
  `sha256`, `md5`, bytes, record count and bp, which equal GTDB's `genome_size` (13,477 of 13,477; *measured*,
  `c36_verify.py`). 13,444 genomes listed (33 accessions repeat another accession's exact contig set and are that
  genome, DEC-27), 30,024 records, 29,790 distinct contigs.
- **Shared, by `ids.contig_uid`** (*measured*, `c36_analyze.py`): **84 contigs are already stored by S1, P1 to P6**
  (none by H1), among them whole chromosomes of phase-1 assemblies; **13,289 are in PL1's inputs** (13,237 of them
  PL1 would store), 1.83 Gbp, plasmids and some chromids; **one** (820 bp) is in two phase-2 stages (C4, C5); 120 are
  listed by several genomes of one stage. Built in the order stored stages, PL1, C1 to C6, the six stages add
  **16,470 contigs and 54,132,065,968 bases**.
- **Estimates** (the code's functions; *derived* unless said): ingest 41.9 to 43.1 CPU-minutes per stage in one
  process; gene calling 150.2 to 153.8 CPU-minutes per stage in one process, single mode, **no batches**
  (under 240 = 0.8 × 300); distinct proteins 3.6 M (C5) to 4.5 M (C3) (*estimate*, from the rehearsal's ratios);
  the full Pfam scan 9,961 (C5) to 12,428 (C3) CPU-minutes per stage in 42 to 52 shards (*estimate*).
- **Build first: C5** (Enterobacterales to *Vibrio*): the most genomes per species (6.64), the largest share of bases
  in genomes after the first of their species (85.4%), the most genomes in species over 20 (988 of 2,032), and the
  fewest estimated distinct proteins per gene (0.418).
- **Rehearsal** (the first 300 genomes of C5 in stage order, as stage `C5R` in a scratch view of S1, P1 to P6, H1;
  run-step in the image by digest): ingest stored the 968 contigs no stored stage holds; gene calling skipped **3
  contigs, exactly the 3 counted** (listed in two genomes of the sample; none of the sample's contigs is stored by
  S1, P or H1); 1,361,815 genes, 548,153 distinct proteins; ingest 6.48 CPU-min (estimate 6.61), gene calling 19.98
  (estimate 23.78). `problems: []`.

## 1. The stages (card (1))

`c36_stages.py assign REMAINDER NOT_SERVED PROPOSAL EXTRACT OUT` (unit `c36-assign`, 20:30:02 to 20:30:05):

- REMAINDER `/tmp/if2453/c27/out_ani/remainder.tsv`, sha256 `48e60728…3098` read on the pod before use, equal to the
  C-27 note's. NOT_SERVED the committed `2026-11-13_s2c_phase2_fetch/not_served.txt` (8 accessions).
- The rule: rows less the not served, ordered by `gtdb_taxonomy`, `species_rank` (as a number), `accession`;
  n = round(56,177,990,736 / 10^10) = 6; a species goes whole to stage 1 + floor(n × bp before its first genome /
  total bp). The script refuses when the fetched accessions (both archives' `kept.tsv`, each row matched to its batch
  lock) differ from the remainder less the not served, when an accession is in two locks, or when a lock's bp is not
  GTDB's `genome_size`. None happened: 4,617 genomes from `s2c_phase2_genomes.tar`, 8,860 from
  `s2c_phase2_rest_genomes.tar` (*measured*).
- Its table is the proposal's byte for byte (`differing_lines: []`); Ck is S2C-k.

| Stage | Genomes | bp | Species | In species over 20 | Genomes per species | bp after the first of its species | Records | FASTA bytes | `.fna.gz` bytes |
|---|---|---|---|---|---|---|---|---|---|
| C1 | 1,783 | 9,370,307,157 | 579 | 361 | 3.08 | 66.5% | 2,929 | 9,487,676,919 | 2,665,666,784 |
| C2 | 3,102 | 9,371,800,221 | 526 | 1,273 | 5.90 | 82.6% | 6,694 | 9,489,508,963 | 2,752,990,698 |
| C3 | 2,470 | 9,502,585,395 | 775 | 646 | 3.19 | 67.3% | 5,689 | 9,621,864,080 | 2,766,477,471 |
| C4 | 1,972 | 9,231,662,797 | 458 | 627 | 4.31 | 77.0% | 4,605 | 9,347,446,902 | 2,671,558,355 |
| C5 | 2,032 | 9,392,003,904 | 306 | 988 | 6.64 | 85.4% | 5,228 | 9,509,853,070 | 2,781,209,046 |
| C6 | 2,118 | 9,309,631,262 | 541 | 722 | 3.92 | 75.7% | 4,879 | 9,426,417,771 | 2,705,667,897 |
| all | 13,477 | 56,177,990,736 | 3,185 | 4,617 | | | 30,024 | 56,882,767,705 | 16,343,570,251 |

All *measured* (`stages/assign.json`, `stages/C<k>/verify.json`); the "all" row and the two ratio columns are
*derived* (sums; genomes / species; bp of genomes that are not the first of their species in stage order / bp).

## 2. The inputs (card (2))

**Extraction** (unit `c36-extract`, 20:23:22 to 20:29:50, host `sha256sum` and `tar`): the sha256 of each archive on
`/stor` equalled `sources.yaml` (`af3ae7c3…fcb9`, `94cc747e…df`); both were extracted into
`/tmp/if2453/s2/phase2.tmp`, renamed to **`/tmp/if2453/s2/phase2/`: 16,359,091,501 bytes, 13,495 files** (13,477
`.fna.gz`, 16 batch locks, 2 `kept.tsv`; *measured*, `du -sb`, `find`). `/tmp` free: 651,377,684 KiB before,
635,125,244 after (*measured*, `df -Pk`); the floor of 400 GB was checked before.

**Per-stage lock and selection** (`c36_stages.py assign`): `stages/C<k>/selection.tsv` holds the remainder's rows of
the stage in stage order; `stages/C<k>/lock.json` each genome's entry of its batch lock with `file` the extracted
`.fna.gz` (the fetch path kept as `fetched_file`, the member as `archived_as`) and the sha256 of every batch lock it
took entries from. These two files are what `make_inputs.py prok` reads in phase 1's form (one lock per stage there
because phase 1's fetch batches were its stages).

**Checks against the fetch** (`c36_verify.py`, six scopes `c36-verify-C<k>` at once, 20:30:27 to 20:35:58): each
`.fna.gz` against its lock's `gz_bytes` and `gz_sha256`, its decompressed bytes against `bytes`, `sha256` and
`md5`, the records `biobox.fasta.iter_records` reads against `records`, their bases against `bp` and GTDB's
`genome_size`. **0 problems in 13,477 genomes** (*measured*). It also wrote each record's `ids.contig_uid` and length
(`contigs.tsv`, the oracle of sections 3 and 6).

**The inputs** (`make_inputs.py prok LOCK SELECTION OUT --source-release datasets-v2-2026-11-13`, six scopes
`c36-inputs-C<k>`, 20:36:15 to 20:38:59), with phase 1's rule for each field:

- kind `prokaryote` and assembly quality `complete` for every genome (fixed by the mode, as phase 1; every genome
  here is GTDB-complete by the selection);
- genetic code from the selection's `translation_table` (GTDB R226 `ncbi_translation_table`), 11 when empty: code 4
  for 353 genomes of C2 (one more is a repeated genome), empty for 313 genomes (24, 36, 117, 36, 40, 60), which take 11;
- topology linear for every contig (the FASTA carries none; phase 1 the same); gene calling then uses no wrap flanks;
- one source object per genome, the extracted `.fna.gz` (source `ncbi_genomes`, licence "INSDC public data",
  redistributable, `retrieved_at` the lock's `retrieved_utc`, 2026-11-13T03:23:05Z to 17:42:51Z). The release label
  is the UTC date of every retrieval in both fetches, as phase 1's `datasets-v2-2026-11-11`.

Dry check (`s2_check_inputs.py C<k>`, scopes `c36-check-C<k>`, 20:39:2x to 20:39:48): exit 0 for all six, genomes
equal to distinct genome uids, 0 circular rows, kind `[prokaryote]`, quality `[complete]`, every object a file
(*measured*, `inputs/C<k>.check.json`).

| `/tmp/if2453/s2/inputs/` | Genomes | Contigs | Memberships | Contigs in several genomes | Genomes repeated | Objects | sha256 `genomes.json` | sha256 `objects.json` |
|---|---|---|---|---|---|---|---|---|
| C1 | 1,776 | 2,918 | 2,922 | 4 | 7 | 1,783 | `32ec62f5a50646824a2b6fbd2f3d373902fe73bfef4d5973c688b8b5e89909b9` | `f8a56a5ce4375a6f9bb0e9f7ac9c73e371607941809d2dac1e999b201c4688de` |
| C2 | 3,096 | 6,675 | 6,687 | 11 | 6 | 3,102 | `7e3f4cfdd6aad5643a9bca7d74b28f12f208c46c6d447514c7ff0c42ef755194` | `deedad548ad8ea08f20fc00ce55aa13b137aec9eb3ae4c1725c6cb6454582248` |
| C3 | 2,456 | 5,579 | 5,671 | 47 | 14 | 2,470 | `4564b44ce1a65dbe943b1b3fe0865f7049ee234e9c88eb9dba1e30edf37b6db3` | `74106398f7e94d459aa9a767053978bfca3639d9000ecfcde9de07269160b34e` |
| C4 | 1,970 | 4,571 | 4,602 | 23 | 2 | 1,972 | `d8fb4df529fa81b55a6df805af038bba13cd2d8d7486104328fc7b2d0738c1e8` | `46891a790939c115f5578f9219661d4f671269a4ffa5a2ddba0625b2368e2988` |
| C5 | 2,029 | 5,203 | 5,224 | 13 | 3 | 2,032 | `1c65060c80ef120d4550d702a76c4a9c5c3eff17ee25f2694d5a172873468573` | `09d742d70a05528e0bd91b735605b63e0f67ad3758e9f1d02107be30bf5483ad` |
| C6 | 2,117 | 4,845 | 4,878 | 22 | 1 | 2,118 | `605ab8c841f2e92f4dfa31a7160da60eff8c7c0551f963a0ba061435a39789b3` | `cb04620345090ec5eae819baf3e72fce987af0f0942ebf0664333d22c1dbaeaa` |

All *measured* (`inputs/C<k>.json`, `inputs/SUMS`). "Genomes repeated" are accessions whose contig set is an earlier
accession's: listed once, every object still ingested (DEC-27; `make_inputs.py`'s rule since card C-21a). The 33 pairs
are in `analyze/summary.json` (`identical_genomes`); most are a `GCF_9023…`/`GCF_9028…` accession beside an older
assembly of the same bytes, or GCA/GCF neighbours (likely resubmissions; not checked).

**The objects are `.fna.gz`, not plain FASTA as phase 1's were** (phase 1's fetch kept plain `.fna`). Kept so because
they are the archived bytes (an object's sha256 in `l0_manifest` will be the archive member's, the lock's
`gz_sha256`), no 57 GB decompressed copy is made, ingest reads `.gz` (`fasta._open`), and `build-stage` sizes ingest
by decompressed bytes (`stage.object_sizes`, `fasta_bytes`). The rehearsal measured ingest of 300 such objects at
6.48 CPU-minutes against 6.61 estimated (section 6).

## 3. What the stages share (card (3))

`c36_analyze.py` (unit `c36-analyze`, 20:40:1x, 6.6 s): ownership from the `contig_index` of S1, P1 to P6 and H1,
read by name (default variants, DEC-26; PL1, being built beside this, is never read from the root); PL1's contigs from
its inputs' `genomes.json` (card C-26 part 1). Check: the 7,046 PL1 contigs owned by a stored stage equal the set of
PL1's `shared_records.tsv` (`owned_by_stored_equals_its_shared_records: true`), so PL1 would store 61,647, as C-26
reported. Build order assumed: the stored stages, then PL1, then C1 to C6.

| Stage | Distinct contigs | Bases | Stored by S1, P, H1 (contigs; bases; by stage) | In PL1's inputs (contigs; bases; of them stored by S1/P) | In another phase-2 stage | Listed by several genomes of the stage (contigs; bases; most genomes) | **Adds after PL1, in order** (contigs; bases) |
|---|---|---|---|---|---|---|---|
| C1 | 2,918 | 9,354,298,875 | 13; 41,107,061; P2 3, P3 6, P4 1, P6 3 | 1,033; 199,069,722; 1 | 0 | 4; 162,942; 2 | **1,873; 9,114,195,377** |
| C2 | 6,675 | 9,354,111,682 | 12; 4,622,210; P4 3, P5 5, P6 4 | 3,230; 222,738,775; 10 | 0 | 11; 123,359; 3 | **3,443; 9,126,873,544** |
| C3 | 5,579 | 9,444,719,398 | 25; 26,818,896; P1 1, P2 4, P3 4, P4 6, P5 4, P6 6 | 2,248; 610,051,759; 12 | 0 | 47; 10,515,825; 8 | **3,318; 8,808,991,551** |
| C4 | 4,571 | 9,217,916,766 | 14; 10,130,560; P2 5, P3 5, P4 4 | 1,930; 403,397,602; 12 | 1 (C5; 820 bp) | 23; 2,982,568; 3 | **2,639; 8,806,059,807** |
| C5 | 5,203 | 9,376,507,937 | 7; 8,542,539; P1 1, P2 3, P4 2, S1 1 | 2,547; 258,792,564; 4 | 1 (C4; 820 bp) | 13; 659,867; 4 | **2,653; 9,109,234,757** |
| C6 | 4,845 | 9,304,426,824 | 13; 425,185; P2 5, P4 1, P6 7 | 2,301; 137,715,892; 13 | 0 | 22; 1,087,756; 5 | **2,544; 9,166,710,932** |
| all six | 29,790 | 56,051,980,662 | 84 | 13,288 distinct (13,289 by stage) | 1 | 120 | **16,470; 54,132,065,968** |

All *measured* (`analyze/summary.json`); the "all six" row is the script's own union, and the sums of the per-stage
"adds" equal it (*derived*). Notes:

- **H1 owns none** of these contigs (H1 holds phages over 200 kbp). S1 owns one (a phage in a C5 genome).
- **The 84 stored contigs include whole chromosomes** (for example an 8,983,317 bp contig of a C1 genome stored by
  P6, and 3.9, 3.6 and 3.3 Mbp ones by P3): a phase-2 assembly that repeats the bytes of a phase-1 genome's
  chromosome. 52 of the 84 are also in PL1's inputs (plasmids of related isolates, as card C-20 found).
- **PL1's inputs hold 13,237 contigs that PL1 would store first** (*derived*: 13,289 − 52), 1.83 Gbp; some are over
  1 Mbp (PLSDB's records include chromids and megaplasmids up to 3.5 Mbp here). After PL1 is built, the C stages'
  gene calling skips them (DEC-26): those calls are PL1's, made in meta mode, while the C genome around them is called
  in single mode.
- **Within a stage**, 120 contigs are in more than one genome (C3 most: 47, one in 8 genomes); the genome of least
  `genome_uid` owns each (DEC-27).
- **Genomes not called at all** (every contig stored before them; `genecall.stage_items`' rule): C1 12, C2 2, C3 9,
  C4 2, C5 2, C6 0 (*derived*, section 4); they stay genomes of their stage in `genome_input`.

**Against PL1 as stored** (added 2026-11-16, 06:52 CST): the table above takes PL1's contigs from its inputs because
PL1 was not stored when it was made. PL1's build has since published its `inputs` and `ingest` steps in the root
(run records `runs/stage=PL1/inputs.json`, `ingest.json`, 2026-11-15 21:06 and 21:32) and was then stopped. The same
script with `--pl1-stored` (unit `c36-analyze3`) read that `contig_index`: **61,647 default variants, exactly PL1's
input contigs less the 7,046 a stored stage owns** (0 only in the root, 0 only expected), and in each phase-2 stage
the contigs PL1 stores are exactly the ones counted from the inputs: 1,032, 3,220, 2,236, 1,918, 2,543, 2,288
(*measured*, `analyze/summary.json` `pl1_stored`). Every other number of the summary was unchanged by the rerun.

## 4. Estimates per stage (card (4))

From the code's own functions (`c36_analyze.py`): ingest `stage.ingest_cpu_min(FASTA bytes)` (42.9 CPU-min per
9,596,180,328 bytes plus 0.07, *measured* 2026-11-10) against `SHARD_MARGIN` × 300 = 240; gene calling
`genecall.cut_plan` over the stage's genomes in `genome_uid` order with `stage_items`' rule (a genome is called when a
contig of it is stored by the stage and its least genome is this one), each called genome costing
`genome_cost_us` (single mode 0.7941 CPU-s per Mbp × `STEP_OVERHEAD` 1.23) and expecting `genome_genes` (917.4 per
Mbp), `batch_cpu_min` 240 as `build-stage` passes it; modes from `mode_of`. All *derived*.

| Stage | Ingest CPU-min (batches?) | Gene calling: genomes called / not | Mode | Bases called | One process CPU-min (batches?) | Genes expected | Distinct proteins (*estimate*) | Residues (*estimate*) | Pfam full: CPU-min; shards (*estimate*) |
|---|---|---|---|---|---|---|---|---|---|
| C1 | 42.5 (no) | 1,764 / 12 | single | 9,313,428,041 | 151.9 (no) | 8,544,139 | 4,462,581 | 1,554,049,244 | 12,370; 52 |
| C2 | 42.5 (no) | 3,094 / 2 | single | 9,349,738,726 | 152.5 (no) | 8,577,450 | 3,716,956 | 1,294,392,668 | 10,303; 43 |
| C3 | 43.1 (no) | 2,447 / 9 | single | 9,432,772,486 | 153.8 (no) | 8,653,626 | 4,483,495 | 1,561,332,250 | 12,428; 52 |
| C4 | 41.9 (no) | 1,968 / 2 | single | 9,212,611,695 | 150.2 (no) | 8,451,650 | 3,923,462 | 1,366,306,459 | 10,876; 46 |
| C5 | 42.6 (no) | 2,027 / 2 | single | 9,369,040,755 | 152.8 (no) | 8,595,158 | 3,593,325 | 1,251,339,444 | 9,961; 42 |
| C6 | 42.2 (no) | 2,117 / 0 | single | 9,305,791,557 | 151.7 (no) | 8,537,133 | 4,027,648 | 1,402,588,235 | 11,165; 47 |

- **Ingest and gene calling each fit one process** in every stage; no genome is under 100 kbp, so no meta mode. The
  largest genome costs at most 0.24 CPU-minutes. The bases called count every contig of a called genome (it is
  called whole; DEC-34 (4)); the merge of a batched step would add about 1 CPU-minute, not needed.
- **Genes** are `GENES_PER_MBP["single"]` × bases called. Cross-check: GTDB's `protein_count` sums to 8,493,692,
  9,156,142, 8,744,320, 8,485,527, 8,774,401 and 8,468,134 (*measured*, `assign.json`), within 7% of each.
- **Distinct proteins** (*estimate*): genes on the first genome of each species (genes × the stage's share of bases
  in those genomes) times 0.8891, plus the other genes times 0.3376, the two ratios *measured* on the rehearsal's 300
  genomes (36 species): distinct proteins per gene on the species' first genomes, and new distinct proteins per gene
  on the later genomes (section 6). The rehearsal's species are Enterobacterales; stages of other taxa may differ,
  and C5's own ratios fit C5 best. At the first-genome ratio alone (one genome per species) the stages would have
  7.5 to 7.7 M distinct proteins each.
- **Residues** = distinct × 348.24, the rehearsal's residues per distinct protein (*measured*). **Pfam full** =
  `stage.scan_cpu_min(residues)` (7.96 CPU-min per M residues) in shards of `stage.shard_residues` = 30,150,753
  residues (*derived*). Under DEC-32 and DEC-35 the phase-2 stages scan representatives only and copy rows of
  proteins already scanned: 16.4% of the rehearsal's distinct proteins are already in a stored stage's
  `protein_seq` (*measured*; each species' phase-1 genome is in P1 to P6).

## 5. Which stage to build first (card (5))

**C5.** Of the six it is richest in strains on every measure taken (*measured* `assign.json`; *derived* ratios):

| Measure | C1 | C2 | C3 | C4 | **C5** | C6 |
|---|---|---|---|---|---|---|
| genomes per species | 3.08 | 5.90 | 3.19 | 4.31 | **6.64** | 3.92 |
| share of bases after the first genome of its species | 66.5% | 82.6% | 67.3% | 77.0% | **85.4%** | 75.7% |
| genomes in species over 20 (of the stage) | 20.2% | 41.0% | 26.2% | 31.8% | **48.6%** | 34.1% |
| genomes in species with 10 or more genomes in the stage | 36.6% | 63.2% | 41.4% | 52.9% | **72.1%** | 44.7% |
| estimated distinct proteins per gene (section 4) | 0.522 | 0.433 | 0.518 | 0.464 | **0.418** | 0.472 |

C5 holds 2,032 genomes of 306 species, 87 of them with 10 or more genomes in the stage (*Escherichia*, *Klebsiella*,
*Salmonella*, *Vibrio* among them). C2 is second (Bacillota, 3,102 genomes in 526 species). C5 also has the most
contigs in PL1's inputs after C2 (2,547), so it exercises the skip of PL1's contigs once PL1 is stored.

## 6. Rehearsal (card (6))

`c36_rehearse.sh` (`STAGE=C5`, units `c36-rh-driver-prep`, `-run`, `-check3`), 20:40:4x to 20:54.

- **The view** `/tmp/if2453/c36/rh/root` (`c26_view.py ROOT VIEW S1,P1,P2,P3,P4,P5,P6,H1`, host `python3 -I`): 56
  links to the stored stages' directories (`contig_index`, `seq_locator`, `seq_shards`, `genome_input`, gene calls,
  runs) and `store_params.json`; the sample stage's files land in the view's own directories.
- **The sample**: the first 300 rows of `stages/C5/selection.tsv` (stage order: *Edwardsiella tarda* to *Erwinia*
  sp036865545, whose first 5 genomes only), 36 species, 1,445,377,659 bp; its inputs by `make_inputs.py prok` with the
  stage's rule: 300 genomes, 968 contigs, 971 memberships, 3 contigs in two genomes (*measured*,
  `rehearsal/sample.json`, `report.json`).
- **The launches**: three `bio-box run-step`s through `pod_run.sh` (image `bio-box_faec1f4.sif`, digest
  `sha256:d67300aa…f098a` checked), in `build-stage`'s argv, profile with the real 300 CPU-minute budget, 8 threads:

  | Step (unit) | wall s | CPU-min | peak RSS (run record) | result | estimate |
  |---|---|---|---|---|---|
  | inputs (`c36-rh-inputs`) | 3.4 | 0.07 | 311 MiB | 300 genomes, 968 contigs, 971 rows, 3 shared | |
  | ingest (`c36-rh-ingest`) | 385.5 | 6.48 | 1,258 MiB | 300 objects, 971 records, 968 new sequences, 2 shards | 6.61 (`ingest_cpu_min` of 1,463,528,965 bytes) |
  | genecall (`c36-rh-genecall`, 7 workers, `--seq-stages C5R,S1,P1,...,H1`) | 190.5 | 19.98 | 1,309 MiB | 1,361,815 genes, 548,153 distinct proteins, mode single, 3 skipped | 23.78 (cost of the 300 genomes, `genome_cost_us` plus a process) |

  All *measured* (run records, unit logs in `rehearsal/`); the estimates *derived*. Measured against estimate:
  ingest −2%, gene calling −16%, genes +2.7% (1,361,815 against 1,325,990 = 917.4 × 1,445.38 Mbp), all within
  AGENTS §6's 20%.
- **The checks** (`c36_rh_check.py`, oracle the sample's `contigs.tsv` and `genomes.json` and the stored stages'
  `contig_index`; `problems: []`, *measured*, `rehearsal/check.json`):
  - ingest's `contig_index` holds exactly the 968 contigs no stored stage holds (all of the sample's: **none of the
    first 300 genomes of C5 has a contig that S1, P1 to P6 or H1 stores**; C5's 7 such contigs are later in the
    stage);
  - gene calling's skipped contigs are **exactly the 3 expected**, each the contig two genomes of the sample list,
    skipped in the genome of larger `genome_uid` with the other as owner (DEC-27); 0 owned by a stored stage, as
    counted;
  - no `feature_protein` row on a stored stage's contig; one row per `feature_uid` (1,361,815); `protein_seq` rows
    equal the distinct proteins.
  - The sample's 603 contigs in PL1's inputs (47,196,369 bases) were stored and called by `C5R`, since PL1 is not in
    the view; in the real build after PL1 they are skipped as PL1's.
- **Ratios for section 4** (*measured*): 160,267 genes and 142,497 distinct proteins on the 36 species' first
  genomes (0.8891); 1,201,548 genes and 405,656 new distinct proteins on the other 264 (0.3376); 348.24 residues per
  distinct protein (190,887,760 / 548,153); 942.3 genes per Mbp; 89,973 of the 548,153 distinct proteins (16.4%) are
  in a stored stage's `protein_seq`.
- Scratch: the view's own files of `C5R`, 21 files, 729,722,772 bytes (*measured*, `find … -printf %s`).

## 7. Pod use and scratch

| Unit | CST | What | Threads |
|---|---|---|---|
| `c36-extract` | 20:23:22 to 20:29:50 | sha256 of both archives, `tar -xf` into `/tmp/if2453/s2/phase2` | 1 |
| `c36-assign` (scope `c36-assign-py`) | 20:30:02 to 20:30:05 | `c36_stages.py assign` | 1 |
| `c36-verify` (scopes `c36-verify-C1` to `C6`) | 20:30:27 to 20:35:58 | `c36_verify.py`, six at once | 6 |
| `c36-inputs` (scopes `c36-inputs-C1` to `C6`) | 20:36:15 to 20:38:59 | `make_inputs.py prok`, six at once | 6 |
| `c36-check` (scopes `c36-check-C<k>`) | 20:39:2x to 20:39:48 | `s2_check_inputs.py`, `sha256sum` | 1 |
| `c36-analyze`, `c36-analyze2` | 20:40:1x; 20:53 | `c36_analyze.py` (the second with the rehearsal's ratios) | 1 |
| `c36-rh-driver-prep` | 20:40:4x to 20:41:17 | the view, the sample, its inputs | 1 |
| `c36-rh-driver-run`: `c36-rh-inputs`, `-ingest`, `-genecall` | 20:41:25 to 20:52:04 | the three run-steps | 8 |
| `c36-rh-driver-check`, `-check2`, `-check3` | 20:52 to 20:54 | `c36_rh_check.py` (the first two stopped on defects of the script, fixed) | 2 |
| `c36-analyze3` | 2026-11-16 06:52 | `c36_analyze.py --pl1-stored` (PL1's stored `contig_index`, read only) | 1 |

- Every process in a `systemd-run --user` unit or scope with a memory cap, `nice 19` (`ionice` 7 for the input
  modes); Python in the runtime image `bio-box_faec1f4.sif` with the tree's `src` alone on the path and the data root
  bound read-only. Host tools: `sha256sum`, `tar` (extraction), `python3 -I` for `c26_view.py` (symlinks only), as card
  C-26 did. At most 8 threads of ours at once; the account's load was 10 to 17 (`uptime`, 20:30 and 20:40, the PL1
  agent's rehearsal beside ours), under the night cap of 56. No `STOP` file appeared in the data root.
- Tests: `pod_pytest.sh C-36` (the test image of `faec1f4`); the whole suite at `4e89725` as two detached shards
  (`detached_suite.sh C-36 2`, units `ci-C-36-0`, `-1`, 20:58 to about 23:02): 1,217 passed and 7 skipped, 1,173 passed
  and 16 skipped (*measured*, `ci-tmp/C-36-detached/log.*`). The session was interrupted by the laptop sleeping at
  about 21:50 and resumed at 06:51 on 2026-11-16; the pod's state was read first and no step was repeated.
- Disk: `/tmp` free 651,377,684 KiB at 20:23, 635,125,244 after extraction, 630,460,860 at 20:54 (*measured*,
  `df -Pk`); never near the 400 GB floor.
- **Left on the pod**: `/tmp/if2453/s2/phase2/` (16.4 GB, the extracted genomes: the inputs' objects, needed by the
  build), `/tmp/if2453/s2/inputs/C1/` to `C6/` (the inputs, about 1 MB each), `/tmp/if2453/c36/` (stage locks and
  selections, contigs tables, the rehearsal view, about 0.76 GB). Nothing in the data root, nothing under `/stor`.
  The fetch directories `/tmp/if2453/c25/genomes` and `/tmp/if2453/c27/genomes` (16 GB of the same files) are no
  longer needed by these inputs.
