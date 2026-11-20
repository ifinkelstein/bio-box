# SPIRE tranche 1 (stage `SPX1`): 304 aquatic and terrestrial samples chosen, fetched and staged

Card C-48e part 1, 2026-11-20 (the pod's clock, UTC). Host `wilkcomp02`, image
`/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif` with the tree of branch `task/C-48e-1` on `PYTHONPATH`
(`/tmp/if2453/c48e1/tree`). Every pod job was a detached unit (`c48e1-fetch`, `c48e1-inputs`, `c48e1-rerun`,
`c48e1-est`) at `nice 19`, `MemoryMax=32G`, `MemorySwapMax=0`, one or two CPUs, the data root bound read-only;
nothing was written to `/tmp/if2453/r10/root` or to `/stor`, and no stage was built. Scripts:
`scripts/stress/build_stage_pod/c48e_select.py`, `c48e_fetch.py`, `c48e_inputs.py`, `c48e_est.py`, `c48e_run.sh`
(the units' modes); tests `tests/r1/test_c48e_tranche.py`. Files beside this note in
[`2026-11-20_spire_tranche_1/`](2026-11-20_spire_tranche_1/): `select/` (the selection) and `out/` (the fetch
log, the inputs' report, check and sums, the measurements, costs and remainder).

Labels (AGENTS §6): *measured* (wilkcomp02, 2026-11-20, the command or file named), *derived* (formula and inputs),
*estimate*, *published*.

## 0. Summary

- **Selected** (*estimate*, before the fetch): 304 samples of 19 aquatic and terrestrial level-2 strata and 189
  studies, 243.05 Gbp assembled, 2,985,971 contigs and 35.467 Gbp at >= 5 kb, against caps of 3.0 M contigs and
  40 Gbp (`select/select_summary.json`).
- **Fetched** (*measured*, `out/fetch_log.jsonl`): 304 of 304 at the first try, 0 HTTP 404, 0 retries, 0 refusals;
  75,267,693,104 bytes gzip, 256,948,810,775 bytes of FASTA; 17:00:58 to 18:50:54 UTC.
- **Staged** in `/tmp/if2453/s2/inputs/SPX1/` (*measured*, `out/inputs_report.json`): the 304 assemblies whole
  (325,805,489 records, 243,047,507,220 bp), 2,829,368 one-contig genomes (the distinct contigs of 5,000 bp or more,
  35,885,221,307 bp). Dry check exit 0, no problem; a second run gave every file byte for byte.
- **Measured against the estimate** at >= 5 kb: contigs 0.948, bases 1.012 (both within 20%); both caps hold
  (2,829,368 <= 3,000,000 contigs; 35.885 <= 40 Gbp). No sample was dropped.
- **Overlap** by `contig_uid` of the 2,829,368 long contigs: 0 with every stage in the root; 12,052 contigs
  (196.6 Mbp) with the staged inputs `G1`, 2 with `PL2`, 7 with `V2`, 0 with `V1`, `PL1`, `C1` to `C6`.
- **Cost of the build** at `--depth called --min-contig-len 5000` (*derived*): 3,196 CPU-minutes by the code's
  estimates (ingest 1,056, gene calling 2,022 in 9 batches plus 3.75 for the merge, CRISPR 100.5, IGR 5.2, tapes
  7.6, inputs 0.8), 3,813 at the pilot's measured rates; about 26.4 GB added to the root. Search pack 33.2 GB with
  `dark_orf`, 9.7 GB without it.
- **Left for later**: 17,968 aquatic and terrestrial samples in tranche order (8.50 Tbp assembled; about 83.3 M
  contigs and 958.6 Gbp at >= 5 kb, *estimate*) and 361 without a published assembly size; 77,434 samples of other
  habitats, as counts. Under 5 kb in the staged assemblies: 99.13% of records, 85.23% of bases (outside the corpus,
  DEC-60 (1)).

## 1. Source, terms and the fetch

Entry `spire_tranche_1` of `sources.yaml` was written before the first request (DEC-50, DEC-58) and completed after.
URL `https://spire.embl.de/download_assembly/{sample_id}`, the pilot's (entry `spire_pilot_assemblies`); terms CC
BY-SA 4.0 (`https://spire.embl.de/faqs`, read 2026-11-18 for the pilot; `redistributable: true`). One connection,
at least 1.0 s between request starts, User-Agent `curl`; a sample tried at most three times on a transient error;
HTTP 403 or 429 stops the whole fetch; `/tmp` kept at 400 GB free or more and the fetched files at 150 GB or less.

*Measured* (`out/fetch_log.jsonl`, sha256 `e4b75b45d05ce335dfe4138b3b34563c90cfe12424da55850ba6cfe48aec4cf3`,
`out/fetch_summary.json`): 304 HTTP 200 at try 1; every answer gzip; served from `spire.embl.de` itself (the pilot
saw redirects to `black.embl.de` and `swifter.embl.de`; none was seen today); largest file 4,938,946,579 bytes
(`SAMN08779656`); 6,580 s of transfer, 11.4 MB/s (*derived*). `/tmp` free went from 678.5 GB at the start to
600.4 GB at the end of the fetch, never under 400 GB (*measured*, `df`).

## 2. The selection (`c48e_select.py run`, seed 20261120)

The rule (the module's docstring; `select/select_summary.json` has it word for word):

1. Universe: the 96,067 SPIRE v1 study samples (`c42_select.merge_sample_tables`, as the pilot); habitats `aquatic`
   and `terrestrial`: 18,633 samples in 19 level-2 strata. 361 of them have no published assembly size
   (`assembly_bp` absent in the API table): they are not ordered but listed apart, because all seven such aquatic or
   terrestrial samples the pilot requested answered HTTP 404 (*measured*, the pilot's
   `dl/assemblies/fetch_assemblies.jsonl`). 18,272 are ordered.
2. Within a stratum, studies in turn: studies in the order of sha256("20261120:study:<study>"), each study's samples
   in the order of sha256("20261120:sample:<sample_id>"), one sample from each study per round. A pure function of
   the table: no order depends on the reading order.
3. Strata interleaved by Webster's divisor method on weights sqrt(N_h) (N_h every sample of the stratum, as the pilot
   counted): the next sample comes from the stratum with the largest sqrt(N_h) / (2 k_h + 1). The pilot used
   Hamilton's largest remainder for a fixed 1,000; a divisor method never takes a sample back as the tranche grows,
   so every prefix is a sqrt(N_h) allocation and tranche 2 continues the same sequence.
4. Contigs and bases at >= 5 kb per sample (*estimate*): the pilot's measured values for a sample the pilot fetched
   (12 of the 304), else the ratio of sums over the 15 pilot aquatic or terrestrial samples nearest in ln(N50) (the
   API's `assembly_n50`), times the sample's published assembly size. N50 is the predictor because the share at
   >= 5 kb runs from 0.6% to 68% across the pilot's N50 octiles (*measured* by a one-off exploration script on the
   pod, not committed, over the pilot's `contigs/sample_contig_stats.tsv`; octiles of the 444 aquatic and terrestrial samples with a published
   size), and only from 5% to 36% between the strata of 5 or more pilot samples but the level-1-only `aquatic` one
   (66%, six samples of 1 Mbp).
5. The tranche is the longest prefix of the order whose estimated contigs and bases at >= 5 kb both stay within the
   caps (3,000,000 and 40,000,000,000); `--exclude <earlier selection.tsv>` skips samples an earlier tranche took.

**The estimator, checked before use** (`select/estimate_validation.json`, `c48e_select.py validate`, *measured*):
each of the pilot's 452 aquatic and terrestrial samples estimated from the samples of the other 128 studies. With
k = 15 the estimated totals over all are 0.996 (contigs) and 0.997 (bases) of the measured; over 2,000 resamples of
300 samples the 5th to 95th percentiles are 0.968 to 1.026 (contigs) and 0.958 to 1.034 (bases). (k = 5, 25, 40 in
the same file; 15 has the smallest bias.)

**Totals against the caps** (`select/select_summary.json`, *estimate*): 304 samples, 162 aquatic and 142
terrestrial; 189 studies; 243,047,507,220 bp assembled; 2,985,971 contigs and 35,467,228,896 bp at >= 5 kb. The
next sample in the order, `SAMN07631787` (`terrestrial:wetland`, rank 305), would have added 22,961 contigs and
put the total over 3.0 M. The contig cap binds; the bases end at 88.7% of theirs. The selection was run twice
(`select/` and `/tmp/if2453/c48e1/select_rerun/`): every output has the same sha256 (*measured*).

Per stratum (`select/by_stratum.tsv`; *estimate* at 5 kb, *measured* assembled bases):

| Stratum | N_h | Taken | Assembled Gbp | Est. contigs >= 5 kb | Est. Gbp >= 5 kb |
|---|---:|---:|---:|---:|---:|
| aquatic (level 1 only) | 46 | 4 | 0.005 | 168 | 0.003 |
| aquatic:benthic | 1,652 | 26 | 14.577 | 201,435 | 2.439 |
| aquatic:brackish water | 515 | 14 | 16.448 | 201,995 | 2.287 |
| aquatic:briny water | 55 | 5 | 1.861 | 38,880 | 0.529 |
| aquatic:fresh water | 2,392 | 31 | 17.449 | 294,859 | 3.688 |
| aquatic:ice | 48 | 4 | 1.196 | 20,777 | 0.270 |
| aquatic:lentic | 718 | 17 | 3.301 | 43,063 | 0.501 |
| aquatic:littoral | 30 | 3 | 1.180 | 23,275 | 0.329 |
| aquatic:lotic | 345 | 12 | 2.761 | 99,249 | 1.212 |
| aquatic:marine | 5,340 | 46 | 29.381 | 372,638 | 4.446 |
| terrestrial (level 1 only) | 247 | 10 | 15.099 | 314,334 | 4.079 |
| terrestrial:arid | 272 | 10 | 2.589 | 26,979 | 0.327 |
| terrestrial:forest | 1,846 | 27 | 48.573 | 370,378 | 3.885 |
| terrestrial:grassland | 337 | 12 | 14.015 | 103,525 | 1.068 |
| terrestrial:permafrost | 837 | 18 | 20.554 | 230,504 | 2.682 |
| terrestrial:soil | 1,874 | 27 | 19.405 | 222,070 | 2.621 |
| terrestrial:subsurface | 1,748 | 26 | 18.469 | 274,422 | 3.447 |
| terrestrial:tundra | 1 | 1 | 0.015 | 20 | 0.0002 |
| terrestrial:wetland | 330 | 11 | 16.167 | 147,400 | 1.653 |
| **All** | **18,633** | **304** | **243.048** | **2,985,971** | **35.467** |

## 3. The staged inputs: `/tmp/if2453/s2/inputs/SPX1/`

Name `SPX1`, as the orchestrator proposed: the pilot's scratch stages were `SP0`, `SPA1`, `SPN1`, `SPQ1`, `SPQ2`,
`SPH1`, `SPH2`, `SPT1`, `SPT2`, `SPU1` (`SP` + a habitat letter + a number; `X` was only the scheme's fallback for an
unknown habitat and never used), so `SPX1` cannot collide with them, and one stage spans two habitats. Kind and mode
as the pilot's depth-1 stages; to be built at `--depth called --min-contig-len 5000`.

The form (`c48e_inputs.py build`, the pilot's `c42_inputs.py` form with the length rule of ingest, card C-38d):

- `assemblies/<sample_id>.fa.gz`: the 304 fetched files **whole**, hard links of the fetched files (no copy). The
  length filter is ingest's: `build-stage --min-contig-len 5000` drops the 322,975,207 shorter records and records
  them and their bases in the ingest ledger (`below:5000`, `n_bases_rejected`).
- `objects.json`: 304 source objects, paths in this directory, `spire` `v1`, CC BY-SA 4.0, redistributable,
  `retrieved_at` the fetch time; in sample order.
- `genomes.json`: **one genome per distinct contig of 5,000 bp or more** (2,829,368), kind `prokaryote`, code 11, no
  circular contig, quality `contig`, in `contig_uid` order. Per contig, as the pilot, and not per sample: a sample is
  a mixture of organisms, and one genome per sample would be called in Prodigal's single mode trained on the mixture.
  Contigs under 100 kb are called in meta mode (2,809,610 contigs, 32.557 Gbp); the 19,758 contigs of 100 kb or more
  (3.328 Gbp) in single mode, each trained on itself (`genecall.mode_of`). Ingest requires `genome_input` to list
  exactly the contigs at or over the length (`store._check_genome_input`), which is why only those are genomes.
- `samples.tsv` (per sample: file bytes, sha256, records and bases in all, under and over 5 kb, over 10 and 100 kb),
  `assemblies.sha256`, `input_report.json`. No time, host or path of the run in them.

`inputs_SUMS` (*measured*, `out/inputs_SUMS`):

| File | sha256 |
|---|---|
| `genomes.json` | `407e8b6eb8bb82d1b753fd71f21446b923dd793d61f5d5e4af903a22d4ac25c7` |
| `objects.json` | `de4889d0a48cb4e685b6f3b4d1783865cdd681cc7849eb1760de3e800e8a5341` |
| `samples.tsv` | `428bab8640994e667bd411221cb4811cdba00d7b1cebd5d9c2f8b57811262a77` |
| `input_report.json` | `6cd4258a9d084c370cb66494fe2274db930f18dbfab48ba59518886cd80d2184` |
| `assemblies.sha256` (the 304 assemblies' sha256) | `c12c797b858337f030f3e8993f95f32d0f99cd32a4e590fa63cfb1dde42bc222` |

The directory holds 75,672,465,511 bytes (*measured*, `du -sb`). The distinct long contigs with their lengths are
also in `/tmp/if2453/c48e1/work/contigs_ge5000.parquet` (2,829,368 rows, sha256
`74d710ea8232eb9a5bb6835a353f6a5b8f05c514f2dd985c74eb60c7bb34cec0`; scratch, not committed).

**Dry check** (`s2_check_inputs.py SPX1`, the `inputs` and `ingest` steps' own readers; `out/check.json`,
`out/check_full.json`): exit 0, 2,829,368 genomes, 2,829,368 `genome_input` rows and distinct genome uids, 0
circular rows, kinds `[prokaryote]`, codes `[11]`, quality `[contig]`, 304 objects all files, `problems: []`.

**Determinism** (unit `c48e1-rerun`, `/tmp/if2453/c48e1/logs/rerun.log`): the build again into scratch with the same
`--final-dir` gave `genomes.json`, `objects.json`, `samples.tsv`, `input_report.json`, `assemblies.sha256`, the
contigs Parquet file and the printed report byte for byte (`cmp`), and the same 304 assemblies; the scratch copy was
then removed.

**Bounded memory** (AGENTS §4.1): one record at a time; the long contigs' uids went through a scratch Parquet file and
one DuckDB connection (1 thread, 1 GB, spilling). The unit took 54 minutes on one CPU (18:51 to 19:45).

## 4. Measured against the estimate (`out/measured_vs_estimate.tsv`, `out/per_sample_measured.tsv`)

*Measured*: 325,805,489 records and 243,047,507,220 bp (equal, sample by sample, to the API's `assembly_bp`);
2,830,282 records of 5,000 bp or more (35,905,121,922 bp), 2,829,368 distinct (35,885,221,307 bp); 907 contigs in
two or more samples (914 repeated records), none twice in a sample; 969,912 records of 10 kb or more (23.307 Gbp).

| Level | Est. contigs | Meas. contigs | Ratio | Est. Gbp | Meas. Gbp | Ratio |
|---|---:|---:|---:|---:|---:|---:|
| aquatic (162 samples) | 1,296,339 | 1,273,629 | 0.983 | 15.704 | 16.545 | 1.054 |
| terrestrial (142) | 1,689,632 | 1,556,653 | 0.921 | 19.763 | 19.360 | 0.980 |
| **All (304)** | **2,985,971** | **2,830,282** | **0.948** | **35.467** | **35.905** | **1.012** |

Both totals are within 20% of the estimate. The contigs are 5.2% under it, outside the validation's 5th percentile
(0.968): the deficit is in a few large samples whose nearest-N50 neighbours had more long contigs per base than they
do. Three strata are more than 20% off (AGENTS §6), each from one or two samples: `terrestrial:forest` (contigs 0.784,
bases 0.752: `SAMN08779656`, 15.6 Gbp assembled, estimated 70,104 contigs and 658 Mbp, measured 33,573 and 265 Mbp),
`aquatic:lotic` (contigs 0.751, bases 1.004: `SAMN11532453`, estimated 57,037 contigs, measured 33,377 of the same
total length, so longer contigs), and the level-1-only `aquatic` stratum (contigs 1.470 on 168 estimated: four
samples of 4.9 Mbp in all). A stratum of 3 to 46 samples is not what the estimator was checked on (300-sample sets).

**Caps** on the measured values: 2,829,368 distinct contigs (2,830,282 summed by sample) <= 3,000,000, and
35.885 Gbp <= 40 Gbp. No sample is dropped (`out/est.json`, `drop_from_end: []`).

## 5. Overlap by `contig_uid` (`out/est.json`, `overlap`; read 2026-11-20 about 19:46 UTC)

Of the 2,829,368 distinct contigs of 5 kb or more (*measured*):

| Against | Contigs shared | bp shared |
|---|---:|---:|
| the root's stages on disk (`C5`, `H1`, `P1` to `P6`, `PL1`, `S1`, `V1`; 293,637 default-variant contigs in all) | 0 | 0 |
| staged inputs `G1` (GlobDB r226 non-GTDB archaea) | 12,052 | 196,592,432 |
| staged inputs `PL2` | 2 | 33,456 |
| staged inputs `V2` | 7 | 198,779 |
| staged inputs `V1`, `PL1`, `C1` to `C6` | 0 | 0 |

(`GV1` had no inputs directory when this was read.) The 12,052 contigs shared with `G1` are contigs of SPIRE MAGs
that GlobDB took as species representatives: a MAG's contigs are exact contigs of its sample's assembly (the pilot's
`mag_probe`). Whichever stage is built first owns them (DEC-26); see open questions.

The root's contigs after `SPX1` (*derived*): 293,637 on disk + 2,829,368 = 3.12 M, and about 5.7 M with the stages
the orchestrator has queued (about 2.9 M by its count), under the 7.86 M of card B-56.

## 6. Cost of building `SPX1` at `--depth called --min-contig-len 5000` (`out/est.json`, `costs`)

| Step | Code's estimate (CPU-min, *derived*) | Pilot's measured rate (CPU-min, *derived*) | How |
|---|---:|---:|---|
| inputs | — | 0.8 | pilot 0.27 CPU-min per million genomes |
| ingest | 1,056.1 | 1,254.9 | code: `stage.ingest_cpu_min(256,948,810,775 FASTA bytes, 0 held)`; pilot: `SP0`'s 5.163 CPU-min per Gbp read from `.fa.gz`; 1,271.1 with card C-47's gzip factor 1.2035 |
| gene calling | 2,022.0 in 9 batches (largest 224.7) + merge 3.75 | 2,176.9 | code: `genecall.cut_plan` over the 2,829,368 one-contig genomes in genome_uid order, batches at 240 CPU-min; pilot: 60.662 per Gbp at >= 1.5 kb |
| CRISPR (default) | 100.5 | 367.7 | code: 168.1 CPU-s per Gbp (after card B-55); pilot: 10.245 per Gbp before card B-55's fix |
| IGR | — | 5.2 | pilot 0.144 per Gbp |
| tapes | — | 7.6 | pilot 0.213 per Gbp |
| **Total** | **3,196** (code, with the pilot's rates where the code has none) | **3,813** | |

A measured sample (below) called genes at 3.133 CPU-s per Mbp in meta mode, 6.4% over the code's 2.945 (*measured*
vs. constant): the gene-call estimate holds within 20%.

**Ingest in parts.** `stage.pass1_bytes` cuts pass 1 at 17,100,339,219 bytes of `.gz` per shard (240 CPU-min at
4.11 per GB of FASTA, at this stage's 3.414 bytes of FASTA per byte of gzip): 5 shards (*derived*). The largest holds
55.29 Gbp; at the pilot's `.gz` rate (5.163 per Gbp) it would take 285.5 CPU-min: under the 300 CPU-minute budget
of a process, but above the 240 the cut aims at. Card B-58 (the gzip term in `ingest_cpu_min`) is open; see open
questions.

**GB added to the root** (*derived*): a least-squares fit of the pilot's measured stage bytes (nine `>= 1.5 kb`
stages and `SP0`) on stored bases and contigs gives 0.684 B per base + 653 B per contig for depth `called`
(largest error over the ten 12.6%): **26.4 GB** for 35.885 Gbp and 2.829 M contigs; 12.4 GB of it at depth `stored`
(0.317 B per base + 364 B per contig, largest error 2.6%). At the pilot's 0.8515 GB per Gbp at >= 1.5 kb it would be
30.6 GB.

**Search pack** (DEC-54 (3), the pilot's §4.3). Measured sample (*measured*, `out/est.json`, `sample`): the 11,093
long contigs whose `contig_uid` starts with byte 0 (1/256; 141,555,873 bp), called one contig at a time as the step
calls them: 141,589 genes (1.000 M per Gbp), 11.8% fragments, 141,562 distinct proteins (41,857,836 aa), 73,774
IGRs of 30 nt or more, 1,203,463 dark ORFs (78,746,995 aa; 1,080,194 of them touch a contig's 3 kb end zone).
Scaled to 35.885 Gbp with the pilot's measured pack bytes per residue and per IGR (stage `SPA1`: 0.650 B per unique
aa, 1.177 B per dark-ORF aa, 147.4 B per IGR) (*derived*): unique proteins 6.9 GB, IGRs 2.76 GB, dark ORFs 23.49 GB:
**33.15 GB with `dark_orf`, 9.66 GB without** (0.92 and 0.27 GB per Gbp). The pilot's 1.01 GB per Gbp would give
36.2 GB. Dark ORFs are 71% of the pack even at >= 5 kb, mostly in the end zones (90% of them).

## 7. What the tranche leaves for later (DEC-59 (3))

(a) **Aquatic and terrestrial samples not taken**: `out/remainder_aquatic_terrestrial.tsv.gz` (18,329 rows; sha256
`2266555bd93b1b098b47f78aa7d5f316c8e6de850cd94fe87f5a80711bd8c63b`; committed, 469 KB), in the order a second
tranche takes them (`order_rank` 305 on), with sample id, study, habitat, stratum, N_h, rank in stratum, assembled
bases, N50, the estimate at >= 5 kb and its source, and a status: 17,968 `not_taken` (8,498,642,140,659 bp
assembled; 83,343,010 contigs and 958.6 Gbp at >= 5 kb, *estimate*) and 361 `no_assembly_size` at the end. Tranche 2
is `c48e_select.py run ... --exclude planning/measurements/2026-11-20_spire_tranche_1/select/selection.tsv`, which
continues the same order. By stratum: `select/by_stratum.tsv` (columns `left_*`).

(b) **Other habitats** (`out/other_habitats.tsv`; *measured* counts and assembled bases, *estimate* at >= 5 kb by the
pilot's measured share in the habitat):

| Habitat | Samples | With assembly size | Assembled Tbp | Est. Gbp >= 5 kb | Est. contigs >= 5 kb |
|---|---:|---:|---:|---:|---:|
| host-associated | 63,909 | 62,044 | 12.370 | 3,383.1 | 227.3 M |
| anthropogenic | 9,662 | 9,555 | 1.659 | 276.6 | 20.6 M |
| unlabelled | 3,846 | 0 | — | — | — |
| air | 17 | 17 | 0.0005 | 0.19 | 7,804 |

(c) **Under 5 kb in the staged assemblies** (*measured*, `out/below_length_by_stratum.tsv`): 322,975,207 of
325,805,489 records (99.13%) and 207,142,385,298 of 243,047,507,220 bp (85.23%); aquatic 81.2% of bases, terrestrial
87.5%; by stratum from 23.7% (the four tiny level-1-only `aquatic` samples) and 56.0% (`aquatic:lotic`) to 94.0%
(`terrestrial:forest`) and 99.0% (`terrestrial:tundra`, one sample). These are
outside the corpus (DEC-60 (1)) and are not a remainder to fetch; ingest records them in its ledger.

## 8. Clean-up

`objects.json` names the hard links under `/tmp/if2453/s2/inputs/SPX1/assemblies/`, not `/tmp/if2453/c48e1/dl/`,
so the 304 names under `dl/` were removed (`find ... -links 2 -delete`, only files with a second link); no byte was
freed: the 75.27 GB of assemblies are held by the staged inputs, which are the only copy. Left in
`/tmp/if2453/c48e1/` (100.5 MB, *measured*): `dl/fetch_log.jsonl` and its sum, the selection and its rerun, `out/`,
`work/contigs_ge5000.parquet`, `logs/`, the units' `DONE` markers. `/tmp` had 597.2 GB free at the end (*measured*).

## 9. Proposed register rows (`planning/CATALOG_GOAL.md` §4; not edited by this card)

§4.1 (left out by a selection inside a source we hold):

| Source | Taken | Left for later, in full | Where the remainder is recorded | How to pick it up |
|---|---|---|---|---|
| SPIRE v1 unbinned contigs, tranche 1 (card C-48e part 1) | 304 aquatic and terrestrial samples, staged as `SPX1` (2,829,368 contigs of 5 kb or more, 35.885 Gbp; inputs in `/tmp/if2453/s2/inputs/SPX1/`, not built) | the other 17,968 aquatic and terrestrial samples with a published size (8.50 Tbp assembled; about 83.3 M contigs, 958.6 Gbp at >= 5 kb, *estimate*) and 361 without one; every sample of the other habitats (host-associated 63,909, anthropogenic 9,662, unlabelled 3,846, air 17) | `planning/measurements/2026-11-20_spire_tranche_1/out/remainder_aquatic_terrestrial.tsv.gz` (in tranche order), `out/other_habitats.tsv`, `select/by_stratum.tsv` | `c48e_select.py run --exclude <earlier selection.tsv>` for the next tranche in the same order; fetch by sample id from `https://spire.embl.de/download_assembly/{sample_id}` (`c48e_fetch.py`); stage with `c48e_inputs.py` |

Change to the existing §4.1 row "SPIRE pilot": 12 of the pilot's samples were fetched again and staged in `SPX1`
(`est_source` `measured_pilot` in `select/selection.tsv`).

§4.2 (left out by a filter):

| Filter | Decision | What it leaves out | How to pick it up |
|---|---|---|---|
| Least contig length 5,000 bp, stage `SPX1` | DEC-57, DEC-60 (1) | 322,975,207 records, 207,142,385,298 bp (85.23% of the staged bases) | Not deferred: outside the corpus; the ingest ledger of `SPX1` will record them (`below:5000`) |

§4.3 (sources not yet fetched), row "SPIRE v1 unbinned contigs": "first tranche fetched and staged (`SPX1`, 304
samples, card C-48e part 1); 17,968 aquatic and terrestrial samples and all other habitats left (§4.1)".

## 10. Discrepancies and open questions

1. **DEC-57 and `CATALOG_GOAL.md` §4.2 say 4.81 M contigs of 5 kb or more in the pilot; the pilot's own per-sample
   table gives 4,292,779** (*measured*: the sum of `n_ge_5000` over `contigs/sample_contig_stats.tsv`, 990 samples;
   the pilot note §3.1 says 4.29 M). 4.81 M looks like a transcription error; the share of bases (16.12%, 55.050
   Gbp) is right. This card's estimates used the per-sample table, not the ratio of the brief.
2. **`G1` and `SPX1` share 12,052 contigs.** Built first, `G1` owns them and calls them as parts of MAGs (single mode
   for a genome of 100 kb or more); built first, `SPX1` owns them, one-contig meta-mode calls. The orchestrator's
   order decides.
3. **Ingest's pass-1 shards on `.gz` inputs**: the largest of the 5 is 285.5 CPU-min at the pilot's measured `.gz`
   rate, under 300 but over the 240 the cut aims at (card B-58 open). A lower `cpu_min_per_process` in the build's
   profile for the ingest, or card B-58 first, keeps the margin.
