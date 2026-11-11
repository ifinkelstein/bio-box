# Infernal on intergenic regions against whole-genome scanning, 2026-11-11 (card 0-09c)

For PI cards 0-09 (option e) and 0-08 (question 2, cost). Every number is *measured* on wilkcomp02 on
2026-11-11 in the image `bio-box-test_74b3e05.sif` (Infernal 1.1.5, Rfam 15.1 `Rfam.cm.gz`, sha256
`e2e636e9…8168`, equal to `refdb.yaml`), by the commands in "Commands", unless labelled *derived* (formula
given) or *estimate*. CPU is user + system of the `cmsearch` process (`timed.py`, rusage of the waited
child, all threads). Every `cmsearch` ran with `--cpu 2` in a unit capped at `MemoryMax=16G`,
`CPUQuota=200%`, `nice 19`, at most two at a time.

## The options (18.1 Mbp of 8 genomes; truth = the producer's whole-genome scan)

CPU-hours per 10 Gbp is *derived*: `CPU-s of the runs / genome bases they cover × 1e10 / 3600`, with
genome bases = 18,106,231 (all eight genomes) or 9,031,568 (batch B only). The fraction of bases is the
bases put into `cmsearch` over the genome bases. Recall and precision are against the 1,120 truth hits in
258 families; a truth hit is recovered by a hit of the same model on the same strand that overlaps at
least half of it (the E-02 rule).

| option | bases scanned | CPU-min | CPU-h per 10 Gbp | speed-up | hit recall | family recall | hit precision | spike-ins (E-02) | what it loses |
|---|---|---|---|---|---|---|---|---|---|
| whole genome, producer settings (truth) | 1.004 | 207.5 | **1,910** | 1 | 1 | 1 | 1 | 97/100 | (3 spike-ins, Infernal's filters at Z = 1e5) |
| whole genome, `--rfam` | 1.004 | 202.6 | 1,865 | 1.02 (noise) | 1.000 | 1.000 | 1.000 | 97/100 | nothing: the same pipeline, hit for hit |
| intergenic, margin 0 | 0.118 | 99.2 | **913** | 2.09 | 0.906 | 0.860 | 0.929 | 97/100 | 27 of 364 tRNAs, 7 of 78 riboswitches/T-box, 7 rRNA model hits, 62 of 377 sRNAs; 36 families entirely |
| intergenic, margin 50 | 0.196 | 134.9 | 1,241 | 1.54 | 0.935 | 0.922 | 0.903 | 97/100 | 17 tRNAs, 5 riboswitches, 7 rRNA model hits, 44 sRNAs; 20 families |
| intergenic, margin 200 | 0.411 | 133.8 | **1,231** | 1.55 | 0.983 | 0.965 | 0.997 | 97/100 | 1 tRNA, 18 sRNAs (antisense RNAs and elements inside coding sequence); 9 families |
| intergenic, margin 0, `--notrunc` | 0.118 | PENDING | | | | | | not run | |
| intergenic, margin 200, `--notrunc` | 0.411 | PENDING | | | | | | not run | |
| whole genome, `--hmmonly` (batch B, 9.0 Mbp) | 1.005 | PENDING | | | | | | PENDING | |
| reduced model set (258 models, scaling only) | 1.004 | PENDING | | | | | | (panel subset) | |
| the E-02 panel's 20 models (scaling only) | 1.004 | 2.2 | 20.4 | | 486/486 of the panel families' truth hits | | 1.000 | | |

The headline: **scanning only intergenic regions does not cut Infernal's cost in proportion to the bases.**
At margin 0 it searches 11.8% of the bases for 48% of the CPU; at margins 50 and 200 it searches 19.6%
and 41.1% for 65% and 64%. The reason is measured, not guessed: for every one of the 14,144 short records
`cmsearch` re-searches both ends for truncated hits, and on records shorter than about twice a model's
window that re-search covers the whole record again: summed over the 4,227 models it re-searched
21.1 G residues at margin 0, more than the 18.0 G residues of the main search, against 9.6 M for the
whole genomes (table "Runs"). And `--rfam` buys nothing at our Z: with `-Z 1e5` the default pipeline
already uses the Rfam-level filter thresholds (every filter stage of every model passed the same number
of windows with and without `--rfam`).

## Genomes and store

Eight genomes of the golden set (`golden.yaml`), chosen to span GC content, an archaeon and phages, in two
batches of about 9 Mbp so that each whole-genome `cmsearch` stays under 300 CPU-minutes (E-01's rate gave
about 115 CPU-minutes per batch, *derived*; measured 103 to 104). *B. subtilis* is in because it is the
riboswitch-rich reference; *M. tuberculosis* for high GC and its many antisense RNAs; *S. acidocaldarius*
is the archaeon (dense tRNAs); lambda and T4 are phages, LAK-A1 a jumbo phage with genetic code 15.

| genome | batch | kind | bases | GC |
|---|---|---|---|---|
| ecoli_mg1655 | A | bacterium | 4,641,652 | 0.508 |
| bsub_168 | A | bacterium | 4,215,606 | 0.435 |
| lambda | A | phage (linear) | 48,502 | 0.499 |
| t4 | A | phage (linear) | 168,903 | 0.353 |
| mtb_h37rv | B | bacterium | 4,411,532 | 0.656 |
| spyogenes_sf370 | B | bacterium | 1,852,433 | 0.385 |
| sacidocaldarius_dsm639 | B | archaeon | 2,225,959 | 0.367 |
| lak_a1 | B | phage (code 15) | 541,644 | 0.259 |
| all | | | 18,106,231 | |

The store is new, built for this card in `/tmp/if2453/m09c/root` with the current tree
(`origin/plan/task-cards` at 2556ddd plus this card's scripts), stage `M9C`: `ingest`, `genome_input`,
`biobox.l1.genecall.run` (producer `pyrodigal-3.7.1-gv-0.3.2-p2`: 17,671 genes, one crossing an origin) and
`biobox.l1.igr.run` (IGR version `t-ca596dd85040-igr1`: 14,144 runs, 10,703 of at least 30 nt,
2,127,711 bp), 22.5 s wall and 0.50 CPU-minutes. I did not use `build-stage`: it would also run the Pfam
scan and the release steps, which this card does not need; the two calls are the same builder functions
`build-stage` calls. No earlier golden store with gene calls was reused (the E-01 root has no gene calls);
the E-01 root's `l1_ncrna` rows were read, read only, as a check (below).

## The producer's settings, read from the code

`src/biobox/l1/ncrna.py`, `_cmsearch_argv`: `cmsearch --cpu <n> -Z 100000 --cut_ga --noali --tformat fasta
--tblout <tbl> -o /dev/null <cm> <fna>`. So: Z fixed at 1e5 Mb (DEC-04), each model's GA as the reporting
threshold, Infernal's default ("search space size-dependent") filter pipeline, not `--rfam`, not `--max`,
not `--nohmmonly` (so the 347 Rfam 15.1 models with no base pairs run in HMM-only mode, as the default
does), truncated hits allowed at sequence ends (the default). A circular contig is searched with 5,746
bases of flank on each side (the library's largest W or MAXL) and a hit is kept where its leftmost base
lies on the contig. This card's runner (`infernal_igr.py fasta wg` and `collect`) does exactly that, and
adds only `-o <file>` instead of `/dev/null` so that each model's own CPU line is kept.

Check of the runner against the producer: the E-01 producer run (2026-11-10, `stage ncrna` at `--cpu 8`,
read from `/tmp/if2453/e01/measure/root`) has 1,120 `l1_ncrna` rows on these eight genomes; this card's
whole-genome default run has 1,120 hits, all 1,120 with the same model, strand, start and end (*measured*,
`infernal_report.py`, `e01_check`). The thread count does not change hits.

`cmsearch -h` in the image (Infernal 1.1.5) lists the acceleration options `--max`, `--nohmm`, `--mid`,
`--default` ("run search space size-dependent pipeline [default]"), `--rfam` ("set heuristic filters at
Rfam-level (fast)"), `--hmmonly`, `--FZ <x>` ("set filters to defaults used for a search space of size
<x> Mb"), and `--notrunc`, `--anytrunc`, `--toponly`; `--devhelp` adds `--nohmmonly` ("never run HMM-only
mode, not even for models with 0 basepairs") and the per-stage thresholds `--F1` to `--F6`. Rfam's own
genome annotation recipe, as I remember the Rfam documentation (not fetched; no download is allowed here),
is `cmscan --rfam --cut_ga --nohmmonly` with `-Z` set to the genome size and clan competition; `--nohmmonly`
makes a scan slower, not faster, so it is not a cheaper variant and was not run.

## 1. Truth: whole genome, producer settings

| run | genomes | records | bases searched | CPU-min | wall s | CPU/wall | bases per CPU-s | peak RSS MiB | hits | families |
|---|---|---|---|---|---|---|---|---|---|---|
| wgA_default | A | 4 | 9,097,647 | 104.21 | 3,188.5 | 1.96 | 1,451 | 4,774 | 740 | 193 |
| wgB_default | B | 4 | 9,077,536 | 103.29 | 3,168.9 | 1.96 | 1,457 | 4,027 | 380 | 96 |
| both | 8 genomes | 8 | 18,175,183 | 207.50 | | | **1,454** | | 1,120 | 258 |

- Bases per CPU-second (*derived*): genome bases / CPU-s: 9,074,663 / 6,252.4 = 1,451 (A),
  9,031,568 / 6,197.7 = 1,457 (B), 18,106,231 / 12,450.1 = 1,454 overall. E-01 measured 1,319 over the whole
  golden set at `--cpu 8` (*measured* 2026-11-10); 10% apart, within the 20% rule (AGENTS §6).
- CPU-hours per 10 Gbp (*derived*): 1e10 / 1,454.3 / 3,600 = **1,910** (E-01: 2,105).
- Peak RSS at `--cpu 2`: 4.0 to 4.8 GiB for the whole genomes, 5.5 to 6.4 GiB for the intergenic runs
  (E-01: 12.3 GiB at `--cpu 8`).
- Fixed cost per `cmsearch` call (*measured*, `tiny_default`: 1,000 seeded random bases against all 4,227
  models): 0.477 CPU-minutes = 28.6 CPU-s. E-01's two-point fit had put it at 615 CPU-s (labelled an
  *estimate* there); the direct measurement is 21 times smaller, so the fixed cost is negligible beside a
  9 Mbp batch (0.5% of 6,250 CPU-s).
- 1,120 hits at or above GA in 258 families; by RNA class (the script's heuristic from model names and the
  Rfam description line, `infernal_report.py rna_class`): tRNA 364, sRNA and other 377, rRNA 214 (nine
  rRNA models over every operon: each operon is hit by several SSU and LSU models, card 0-08 question 1),
  riboswitch or T-box 78, cis leader or motif 38, CRISPR repeat 23, housekeeping (tmRNA, RNase P, SRP, 6S)
  20, intron 6.

## 2. Intergenic regions only

Records: the stage's IGR rows (every run of at least 1 nt not covered by a gene call, both strands'
genes counted), each widened by the margin into the flanking genes on both sides, clipped at linear
ends, taken modulo the length on circular contigs, overlapping or touching runs merged (also across the
origin), cut from the + strand (cmsearch searches both strands), each record named by contig and start so
that a hit maps back to `start + offset` (reduced by the contig length past the origin, so a hit crossing
it has `end > contig_len`). Same argv and the same `-Z 100000` as the truth, so E-values and GA decisions
are comparable hit for hit; a hit's bit score depends only on the sequence it is found in.

| run | margin | records | bases searched | fraction | CPU-min | wall s | CPU/wall | residues searched (all models) | residues re-searched for truncated hits | hits | families |
|---|---|---|---|---|---|---|---|---|---|---|---|
| igr_m0_default | 0 | 14,144 | 2,127,711 | 0.1175 | 99.22 | 3,194.9 | 1.86 | 17,987,668,794 | 21,136,384,064 | 1,093 | 229 |
| igr_m50_default | 50 | 14,120 | 3,541,835 | 0.1956 | 134.87 | 4,265.1 | 1.90 | 29,942,673,090 | 30,903,404,424 | 1,160 | 241 |
| igr_m200_default | 200 | 11,582 | 7,439,788 | 0.4109 | 133.78 | 4,207.0 | 1.91 | 62,895,967,752 | 27,069,821,034 | 1,104 | 249 |
| whole genome (A + B) | | 8 | 18,175,183 | 1.0038 | 207.50 | 6,357.4 | 1.96 | 153,652,997,082 | 19,139,712 | 1,120 | 258 |

Residue counts are the sums of each model's "Target sequences" and "re-searched for truncated hits" lines
in the `-o` output (*measured*). Speed-up against the truth (*derived*, 12,450.1 CPU-s / run CPU-s): 2.09
(margin 0), 1.54 (50), 1.55 (200). Margin 200 costs no more than margin 50 although it searches twice the
bases, because merging leaves fewer, longer records (11,582 against 14,120) with less of their length in
the re-searched ends.

Against the truth:

| margin | truth hits recovered | at the exact coordinates | lost | lost: not scanned (inside genes) | lost: cut at a region boundary | extra hits (no truth match) | of them truncated (5' or 3') | families recovered | families lost entirely |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 1,015 (0.906) | 867 | 105 | 51 | 54 | 78 | 77 | 222 of 258 | 36 |
| 50 | 1,047 (0.935) | 928 | 73 | 26 | 47 | 113 | 112 | 238 | 20 |
| 200 | 1,101 (0.983) | 1,090 | 19 | 16 | 3 | 3 | 3 | 249 | 9 |

No truth hit that lay wholly inside a scanned record was lost, at any margin ("scanned whole, not
reported": 0). Every loss is either a hit with no base scanned or a hit cut by a record boundary to less
than half its length (or below GA). The extra hits are fragments of the same loci: nearly all are flagged
truncated by cmsearch, and most are SSU and LSU rRNA model hits on pieces of rRNA operons split by
spurious gene calls (margin 0: 69 of 78; margin 50: 107 of 113).

Lost hits by RNA class (lost of truth):

| class | truth | margin 0 | margin 50 | margin 200 |
|---|---|---|---|---|
| tRNA | 364 | 27 | 17 | 1 |
| rRNA (model hits) | 214 | 7 | 7 | 0 |
| riboswitch or T-box | 78 | 7 | 5 | 0 |
| cis leader or motif | 38 | 2 | 0 | 0 |
| housekeeping (tmRNA, RNase P, SRP, 6S) | 20 | 0 | 0 | 0 |
| intron | 6 | 0 | 0 | 0 |
| CRISPR repeat | 23 | 0 | 0 | 0 |
| sRNA and other | 377 | 62 | 44 | 18 |

Why they are lost (from the gene calls each lost hit overlaps, `compare` output `lost[].genes`):

- **tRNAs under short spurious ORFs.** At margin 0, 14 tRNAs had no base scanned and 13 were cut; every one
  overlaps a gene call, all but one on the opposite strand, 114 to 804 nt long (two longer, 1,560 and
  2,151): pyrodigal calls short ORFs across tRNA genes, most often in the archaeon *S. acidocaldarius* (14
  of the 27). Margin 50 still loses 17; margin 200 loses one (*E. coli* tRNA at 1,287,246, under a 567 nt
  opposite-strand call, cut to 18 of 82 bases).
- **rRNA split by a spurious ORF.** In *S. pyogenes* a 153 nt call inside each of six 23S genes splits
  the LSU_rRNA_eukarya model hit (2,881 nt, about 1,250 bits) into two pieces of about 1,360 nt, each
  scoring above GA but covering less than half: counted lost, and the pieces counted extra. One *E. coli*
  5S rRNA lies under a 324 nt opposite-strand call. All recovered at margin 200.
- **Riboswitches and T-box leaders overlapped by the start of the next gene.** At margin 0: *B. subtilis*
  SAM (under a 162 nt call), T-box (156 nt), Manganese (gene start called into the leader); *M. tuberculosis*
  TPP and SAM-IV and two Glycine aptamers. All recovered at margin 200.
- **RNAs that lie inside coding sequence by nature** are lost at every margin: cis-antisense RNAs of
  *M. tuberculosis* (ASpks 5, ASdes 2, AS1726, Ms_AS-5), *E. coli* AsdA, Ysr276, 5_ureB_sRNA (2), DicF, the
  dnaX frameshift element (DnaX), Erse. These are the 16 "not scanned" of margin 200 and its 9 families
  lost entirely: DnaX, AS1726, ASdes, ASpks, Ms_AS-5, 5_ureB_sRNA, AsdA, Erse, Ysr276.
- Families lost entirely at margin 0 (36): Manganese, sroD, DnaX, SgrS, SAM-IV, sR28, sR5, Actino-pnp, cspA,
  AS1726, AS1890, ASdes, ASpks, Lambda_thermo, RsaE, sR6, SECIS_2, SECIS_3, STnc630, STnc370, STnc180,
  Ms_AS-5, 5_ureB_sRNA, Spy490380c, MicL, AsdA, RalA, AsrC, ES205, Erse, TrxA_thermometer,
  CpoB_ybgF_thermometer, PssrA, Ysr276, Ysr251, npcTB_6715. At margin 50 (20): Manganese, DnaX, SAM-IV,
  Actino-pnp, cspA, AS1726, ASdes, ASpks, RsaE, Ms_AS-5, 5_ureB_sRNA, MicL, AsdA, RalA, AsrC, ES205, Erse,
  Ysr276, Ysr251, npcTB_6715.

PENDING: `--notrunc` on margins 0 and 200 (cost and loss without the truncated-hit re-search).

## 3. Cheaper variants of the whole-genome scan

### (a) `--rfam`

`cmsearch` with the producer's argv plus `--rfam`, both batches:

| run | CPU-min | wall s | hits | same as default |
|---|---|---|---|---|
| wgA_rfam | 102.31 | 3,125.2 | 740 | all 740, same coordinates |
| wgB_rfam | 100.25 | 3,068.8 | 380 | all 380 |

CPU is 1.8% and 2.9% below the default runs, which is the size of run-to-run variation on a shared host,
not a speed-up: for batch B every "Windows passing ... filter" and "Total ... hits reported" line of
all 4,227 models (40,535 lines) is byte-identical between the default and `--rfam` outputs, and so is
every filter threshold the output prints ("expected" P-values, for example SSU_rRNA_bacteria: SSV 0.06,
Viterbi 0.02, Forward 2e-4, CYK 1e-4, in both). So at `-Z 1e5` the default size-dependent pipeline is
already the Rfam-level one: the producer already runs the fast heuristics, and the 3 spike-ins E-02 loses
are lost to them (item 4). A cheaper pipeline would need thresholds stricter than Rfam's (`--F1` and
following, or `--hmmonly`), which costs sensitivity that Rfam's curators did not accept.

### (a') `--hmmonly` (batch B)

PENDING.

### (b) Reduced model set: how cost scales with the number of models

Not a recommendation; the scaling only. Models: every family with a truth hit in these genomes (258) plus
the E-02 panel's 20 (all 20 already among the 258), so 258 models, fetched with `cmfetch -f` from
`Rfam.cm`. And the panel's 20 alone.

| models | run | CPU-s measured | CPU-s predicted from the full run's per-model lines | share of the full library's CPU |
|---|---|---|---|---|
| 4,227 | wgA_default + wgB_default | 12,450.1 | 12,403.2 (sum of all) | 1 |
| 258 | wgA_subset + wgB_subset | PENDING | 2,029.4 | 0.163 (predicted) |
| 20 | wgAll_panel20 | 133.2 | 135.1 | 0.0107 |

Cost is additive over models: each model's own "# CPU time" line in the full run predicts the 20-model run
within 1.4%. It does not scale with the count alone: the 258 models are 6.1% of the library and about
16% of its CPU, because they include the long rRNA models (the library's models named rRNA cost 866 CPU-s
of 12,403 in the full run, 7%; LSU_rRNA_eukarya alone 315 CPU-s, the most expensive model) and tmRNA and
RNase P. The 3,969 models with no hit in these genomes cost 84% of the CPU (*derived*, 1 - 2,029.4 /
12,403.2). By class (sum of per-model CPU in the full run, *measured*): sRNA and other 10,572 CPU-s
(85%), rRNA 866, cis leader or motif 601, housekeeping 111, riboswitch or T-box 102, intron 79, CRISPR 63,
tRNA 9.

## 4. The E-02 spike-in panel under each variant

The planted contig is built by `infernal_spikein.py emit`, which imports `PANEL`, `_emit` and `_plant` from
`tests/r3/test_g04_spikeins.py` (not copied) and fetches the 20 models as the test does, so it is the
test's contig (68,794 bp) and truth (100 planted sequences). The whole linear contig is one
intergenic-style record (`region.fna`), so a scan of it isolates the cmsearch settings; it is scanned with
the panel's 20 models at `-Z 100000` (the test's own justification: GA is a bit score and Z is fixed, so a
hit's decision does not depend on the other models).

| variant | recovered | exact 95% interval | hits on background | missed |
|---|---|---|---|---|
| producer settings (`g04_default`) | 97/100 | 0.915 to 0.994 | 0 | RF00005 18758-18837 (+), RF00005 29073-29145 (-), RF00029 27517-27600 (-) |
| `--rfam` (`g04_rfam`) | 97/100 | 0.915 to 0.994 | 0 | the same three |
| `--hmmonly` (`g04_hmmonly`) | PENDING | | | |

Both reproduce the test's 97 of 100 and its three misses exactly. Since `--rfam` is the default here, the
3% that E-02 attributes to "Infernal's filter pipeline at Z = 1e5" is the Rfam-level pipeline's own loss.

## 5. Reading (my interpretation)

PENDING the last runs; draft: Intergenic scanning is not the cheap cut card 0-09 option (e) assumed. Its
cost is set by the number of records as much as by the bases, so the margin that keeps recall (200 bases:
98.3% of hits, 96.5% of families, all tRNA, rRNA, riboswitch, housekeeping and leader hits but one tRNA)
saves only a third of the CPU (1,231 against 1,910 CPU-hours per 10 Gbp); margin 0 halves it but loses one
tRNA in 13 and one riboswitch in 11, mostly to spurious short ORFs, so a gene caller's mistakes would
become ncRNA losses. What intergenic scanning can never find, at any margin, is the cis-antisense RNAs
and elements inside coding sequence (9 families here, 16 hits in *M. tuberculosis* and *E. coli*). The
`--rfam` switch is already on in effect. The cost is in the 3,969 models that hit nothing here (84% of
CPU), so the larger lever is the model set (card 0-08 option 2b), which cuts cost in proportion to the
models' own CPU, not to their number; and the two levers multiply only partly, because the per-record
overhead of intergenic records is per model too.

## What was not measured, and why

- Rfam's own `--nohmmonly` recipe and clan competition: slower, not a cheaper variant; clan competition is
  card 0-08 question 1, not cost.
- Intergenic plus `--rfam`: not run; `--rfam` is identical to the default at this Z (item 3a), so it would
  repeat the intergenic rows. Lane 1's `igr_m200_rfam` was stopped at its start for this reason and lane
  2's `--rfam` intergenic jobs were dropped (logs `lane1.log`, `lane2.log`).
- Metagenome contigs, other gene callers, and eukaryotic genomes: not in this subset.
- The spike-in panel under the intergenic margins in its real sense (planted RNAs next to genes): the card
  asked for the whole contig as one region; the panel contig has no genes.

## Commands

On the laptop: `git -C <worktree> checkout -b measure/0-09c origin/plan/task-cards`; the tree (without
`.git` and `planning/`) was copied to `/tmp/if2453/m09c/tree` with `rsync -a --delete`.

On the pod, each lane in its own unit (scripts in `scripts/stress/annotation_cost/`):

```bash
systemd-run --user --unit=m09c-lane<k> --collect --quiet -p MemoryMax=16G -p MemorySwapMax=0 -p CPUQuota=200% \
  -p StandardOutput=append:/tmp/if2453/m09c/logs/lane<k>.log -p StandardError=append:/tmp/if2453/m09c/logs/lane<k>.log \
  nice -n 19 bash /tmp/if2453/m09c/tree/scripts/stress/annotation_cost/infernal_lane.sh <k>
```

Lanes 1 and 2 started at 07:03 CST; lane 2 was stopped at 10:58 while it waited and replaced by lane 3;
lane 1 was replaced by lane 4 at the start of its last job (see "not measured"). Before each job the lane
summed `ps -u if2453 -o pcpu=` (without the ps line) and started only when that plus 200% was at most
4000% on weekdays 08:00 to 20:00 CST, else 5600%; the account was at 997% to 4,419% at the starts (every
start in the logs), and no job had to wait. Every job runs in the image:

```bash
singularity exec --cleanenv --bind /tmp/if2453 --env PYTHONPATH=/tmp/if2453/m09c/tree/src:/tmp/if2453/pilot_s1/pydeps \
  --env OMP_NUM_THREADS=2 --pwd /tmp/if2453/m09c/work bio-box-test_74b3e05.sif python .../timed.py runs/<name>.json -- <command>
```

- `prepare`: `python infernal_igr.py prepare <tree> /tmp/if2453/golden /tmp/if2453/m09c/root /tmp/if2453/m09c/scratch`
- library: `gunzip -c /tmp/if2453/golden/rfam/Rfam.cm.gz > lib/Rfam.cm`
- records: `python infernal_igr.py fasta root lib/Rfam.cm wg A|B|all fa/wg_<b>.fna fa/wg_<b>.map.json` and
  `python infernal_igr.py fasta root lib/Rfam.cm igr all fa/igr_m<m>.fna fa/igr_m<m>.map.json <m>` for m = 0, 50, 200
- every scan: `cmsearch --cpu 2 -Z 100000 --cut_ga --noali --tformat fasta --tblout out/<name>.tbl -o out/<name>.out [flags] <cm> <records>.fna`,
  then `python infernal_igr.py collect out/<name>.tbl <records>.map.json out/<name>.hits.jsonl` and
  `python infernal_igr.py permodel out/<name>.out out/<name>.permodel.json`

| name | cm | records | flags |
|---|---|---|---|
| tiny_default | Rfam.cm | fa/tiny (1,000 seeded random bases) | |
| wgA_default, wgB_default | Rfam.cm | fa/wg_A, fa/wg_B | |
| igr_m0_default, igr_m50_default, igr_m200_default | Rfam.cm | fa/igr_m0, m50, m200 | |
| wgA_rfam, wgB_rfam | Rfam.cm | fa/wg_A, fa/wg_B | `--rfam` |
| wgB_hmmonly | Rfam.cm | fa/wg_B | `--hmmonly` |
| igr_m0_notrunc, igr_m200_notrunc | Rfam.cm | fa/igr_m0, fa/igr_m200 | `--notrunc` |
| wgA_subset, wgB_subset | lib/subset.cm (`infernal_igr.py subset`, then `cmfetch -f -o lib/subset.cm lib/Rfam.cm lib/subset.keys`) | fa/wg_A, fa/wg_B | |
| wgAll_panel20 | work/g04/panel.cm | fa/wg_all | |
| g04_default, g04_rfam, g04_hmmonly | work/g04/panel.cm | work/g04/region (`infernal_spikein.py emit <tree> lib/Rfam.cm work/g04`) | none, `--rfam`, `--hmmonly` |

- spike-in scores: `python infernal_spikein.py score <tree> work/g04 out/g04_<v>.hits.jsonl`
- E-01 check: `python infernal_igr.py e01 root /tmp/if2453/e01/measure/root out/e01_truth.hits.jsonl`
- comparison and tables: `python infernal_report.py /tmp/if2453/m09c /tmp/if2453/m09c/report.json` (in a
  `systemd-run --user --scope -p MemoryMax=4G -p CPUQuota=100%` scope at nice 19)

Outputs stay on the pod under `/tmp/if2453/m09c` (`runs/` run records, `out/` tables, `-o` outputs, hits,
`out/merged/*.compare.json`, `report.json`); nothing was written under `/stor`.
