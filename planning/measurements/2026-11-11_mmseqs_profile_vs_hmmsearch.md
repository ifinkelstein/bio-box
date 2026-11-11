# MMseqs2 profile search against the producer's `hmmsearch` at Pfam GA (card 0-09a, wilkcomp02, 2026-11-10/11)

Card 0-09 (**H**, PI decision) option (d): a faster first pass for Pfam. This file measures how much of the
producer's Pfam 38.2 domain table MMseqs2 recovers on stage P1, and at what cost. Every number is
*measured* on wilkcomp02 between 2026-11-10 23:10 and 2026-11-11 01:00 CST by the commands below, unless
labelled *derived* (formula given) or *estimate*. No code under `src/` or `tests/` changed; nothing was
written under `/tmp/if2453/r10` or `/stor`; nothing was downloaded.

## Summary

Truth: the producer's `l1_domains` rows of stage P1 for `pfam-38.2` (1,552,653 domain rows,
1,447,450 (protein, model) pairs, 918,765 proteins, 294,516,869 residues). MMseqs2 18.8cc5c, proteins as
queries against 30,134 profiles built from Pfam-A.hmm alone (route below), `-e 100`, 4 threads.

| | `hmmsearch` (producer) | MMseqs2 `-s 5.7` | MMseqs2 `-s 7.5` |
|---|---|---|---|
| CPU-min, stage P1 | 2,191.21 (9 shards) | 23.13 | 79.82 |
| speed-up (*derived*, 2,191.21 / CPU-min) | 1 | 94.7x | 27.5x |
| speed-up with the one-off profile build (23.86 CPU-min) | 1 | 46.6x | 21.1x |
| pair recall, any hit (E <= 100) | 100% (truth) | 55.4% | 75.0% |
| best F1 over the bits threshold (bits >=, precision, recall) | - | 0.689 (32, 95.2%, 54.0%) | 0.814 (32, 93.8%, 71.8%) |
| recall at precision >= 99% (bits >=) | - | 49.7% (39) | 62.6% (40) |
| recall / precision at E <= 1e-3 | - | 49.0% / 99.3% | 62.4% / 99.1% |
| recall of pairs >= 50 bits above GA (701,508 pairs) | - | 78.4% | 92.7% |
| recall of pairs 0 to 5 bits above GA (101,206 pairs) | - | 17.5% | 35.0% |
| truth domains covered by a hit of any model (stand-in for clans) | - | 65.5% | 83.2% |
| median abs. start / end difference to the envelope (residues) | - | 1 / 1 | 1 / 1 |
| prefilter by protein (E <= 100): `hmmsearch` residues left, recall ceiling | - | 71.6%, 80.8% | 83.3%, 93.1% |

CPU-min of MMseqs2 = search + `convertalis` (22.80 and 79.49) + `createdb` of the proteins (0.33).
Precision is agreement with the truth, not correctness: 85% of the calls that are not truth pairs at the
best-F1 threshold sit on a truth domain of another model (below).

**Reading (my interpretation, not a measurement).** With profiles derived from the Pfam HMMs (the only route
open without a download), MMseqs2 is 27 to 95 times cheaper than the producer's `hmmsearch` but recovers
only 55% to 75% of its (protein, model) pairs, and it loses most of what sits near the gathering cutoff
(18% to 35% of pairs within 5 bits of GA; even 7% of pairs 50 or more bits above GA are missed at
`-s 7.5`). Raising sensitivity past `-s 7.5` does not help: MMseqs2's own exhaustive mode on a 2,000-protein
sample recovers 72.8% against 76.2% for `-s 7.5`, so the ceiling is the converted profiles and MMseqs2's
single local alignment score, not the k-mer prefilter. Used as a prefilter it saves little: proteins with
any MMseqs2 hit hold 83% of the residues, so `hmmsearch` work falls by 17% at a cost of 7% of truth pairs.
As measured here, option (d) does not replace `hmmsearch` for the corpus-wide family index; it could serve
a clearly labelled, lower-recall bulk tier ("Pfam-like, MMseqs2") where 62% recall at 99% agreement is
acceptable. One untested reason for caution in the other direction: MMseqs2's published Pfam database is
built from the Pfam-A.full alignments (`convertmsa` + `msa2profile`), whose profiles are sharper than
profiles sampled from entropy-weighted HMMs; that route needs Pfam-A.full or Pfam-A.seed, which are not in
`refdb.yaml` and were not fetched, so the recall of the standard MMseqs2 Pfam database is not measured
here. If the PI wants it, adding Pfam-A.seed to `refdb.yaml` and repeating these commands is about two
hours of pod time.

## Setup

- Host wilkcomp02 (112 CPUs, 1007 GiB RAM). Image `/stor/home/if2453/bio-box/images/bio-box-test_74b3e05.sif`
  (`sha256:c30be58557d756095346308c58cc3098af4bbe66404f67ab87be3f64f147c898`, the digest in the producer's run
  records): MMseqs2 `18.8cc5c` (`mmseqs version`), HMMER 3.4 (`hmmemit -h`), Python 3.14.10, pyarrow 26.0.0,
  DuckDB 2.0.0. The image has no `/usr/bin/time`; costs come from `os.wait4` (`timed.py`), as `run-step` does.
- Data root `/tmp/if2453/r10/root` (built by the tree `/tmp/if2453/ci/r10`, commit `740d049`), read only.
  Stage P1's Pfam step was published at 05:14 UTC (`runs/stage=P1/domains-pfam-38.2.json`, exit 0).
- Pfam 38.2: `/stor/scratch/Wilke/if2453/all-the-data/l0/refdb/pfam/38.2/Pfam-A.hmm.gz` (bound read-only),
  30,134 models (`pfam_models.py models`). The directory holds only `Pfam-A.hmm.gz`: no `Pfam-A.clans.tsv`,
  `Pfam-A.hmm.dat`, `Pfam-A.seed` or `Pfam-A.full`.
- Scripts: `scripts/stress/annotation_cost/` (this branch), mirrored to `/tmp/if2453/m09a/bin/`. Scratch
  `/tmp/if2453/m09a/{data,work,runs,logs,test}`.
- Every job ran through `launch.sh`: unit `m09a-<name>` (`systemd-run --user --collect -p MemoryMax=32G
  -p MemorySwapMax=0 -p CPUQuota=<threads x 100>%`), `nice -n 19 singularity exec --cleanenv --bind /tmp/if2453
  --bind <pfam dir>:ro ... python timed.py runs/m09a-<name>.json -- <command>`, at most 4 threads, one job of
  4 threads at a time (two jobs together only when their threads summed to at most 4). A transient service
  with these properties instead of `--scope`, so the job is detached from the ssh session; the memory, swap
  and CPU caps are the same. Before every launch `launch.sh` summed `ps -u if2453 -o pcpu=` (the `ps`
  process itself left out: its lifetime pcpu reads several hundred %) and started only if that plus the
  job's threads x 100 was at most 5600% (all launches were between 23:26 and 00:55 CST). Every launch and
  its reading is in `/tmp/if2453/m09a/logs/launch.log`: readings 3,992% to 4,870%; one launch waited (23:23,
  5,525% with the `ps` line counted) before the `ps` line was excluded. The smoke tests of the profile routes
  (`/tmp/if2453/m09a/test`, 3 models, 60 queries) ran with `systemd-run --user --scope -p MemoryMax=32G
  -p MemorySwapMax=0 nice -n 19 singularity exec ...`, 1 thread. No process passed 300 CPU-minutes: the
  largest was one P1 `-s 7.5` shard, 40.47 CPU-min for the whole shard.
- At the same time the account ran the r10 Pfam scan of stage P2 (12 shard processes of 4 threads) and other
  agents' tests (account total about 40 to 49 CPUs, above). The producer's P1 shards ran under similar
  sharing (8 shards of 4 threads at once). CPU-minutes are comparable; wall times are not rates.

## Truth

`export_stage.py /tmp/if2453/r10/root P1 /tmp/if2453/m09a/data/P1` (unit `m09a-export-P1`, 0.23 CPU-min).
Files are the producer's own lists, not globs: the `protein_seq` file named in the stage's shard plan
(`l1/domain_scan/library=pfam-38.2/stage=P1/plan.json`), and the `l1_domains` file named in the step's run
record, checked against the record's sha256:

- `l0seq/protein_seq/gene_call_version=pyrodigal-3.7.1-gv-0.3.2-p2/stage=P1/protein_prefix=all/part-0.parquet`,
  209,946,839 bytes, sha256 `c9268d41...683d`: 918,765 proteins, 294,516,869 residues (both equal the plan's).
- `l1/domains/library=pfam-38.2/stage=P1/part-0.parquet`, 80,777,664 bytes, sha256 `fb4a545b...e539`:
  1,552,653 rows (the same count as the C-10 build of P1). The producer (`biobox.l1.domains.scan` at
  `740d049`) runs `pyhmmer.hmmsearch` with `bit_cutoffs="gathering"`, `Z = domZ = 1e8`, and keeps only
  included hits and included domains; a row is one domain (envelope `env_start`, `env_end`, 0-based
  half-open), with `seq_bitscore` and `dom_bitscore`.

Exported: `proteins.fa` (uid hex, sequence with `*` removed, as `scan` does), `proteins.tsv`, `truth.tsv`
(uid, model accession, env and ali coordinates, sequence and domain bits). 1,447,450 truth pairs; 750,383
proteins with a domain, holding 88.8% of the residues; 67,534 pairs with more than one domain; no pair
below GA (margin < 0: 0).

`hmmsearch` cost: the producer's nine shard run records `runs/stage=P1/domains-pfam-38.2.shard-0000{0..8}.json`
(`biobox.cli stage domains shard ... --cpus 4`, unit `r10-p2-P1-3`, started 03:18 and 04:19 UTC): 245.47,
247.78, 244.26, 246.79, 247.20, 246.47, 245.90, 246.61, 220.75 CPU-min, sum **2,191.21** CPU-min (printed by
`python3 costs.py /tmp/if2453/m09a /tmp/if2453/r10/root P1`), shard walls 3,291 to 3,727 s at 4 threads, peak
RSS 1.06 to 1.28 GiB. Plan and merge add 0.13 and 0.13 CPU-min. Rate 7.44 CPU-min per million residues
(*derived*, 2,191.21 / 294.52), 6.5% under the 7.96 of `2026-11-10_build_stage_pod.md` (within AGENTS §6's
20%). Stage S1 for reference: 2,042.37 CPU-min for 271,831,237 residues (7.51 per million, *derived*).

## Profiles: what MMseqs2 18.8cc5c can read

`mmseqs -h` lists three ways to make a profile database: `convertmsa` ("Convert Stockholm/PFAM MSA file to a
MSA DB") then `msa2profile` ("Convert a MSA DB to a profile DB"); `convertprofiledb` ("Convert a HH-suite HHM
DB to a profile DB"); and `databases` ("List and download databases"), which downloads and is forbidden here.
Nothing reads HMMER3 files. The two local routes tried:

1. **HMMER3 to HHM, then `convertprofiledb`: does not work.** `pfam_models.py hhm` writes each model as HHM
   text (HMMER's -ln p converted to HH-suite's -1000 log2 p; consensus from the CONS column) into an
   MMseqs2-format DB. `mmseqs convertprofiledb three_hhm three_prof` accepted it and wrote a profile of the
   right length (2,784 bytes = 121 positions x 23 + 1 for the 121-state PF26733.1), but the PSSM
   (`mmseqs profile2pssm`) is not the model's (first position of PF26733.1, where V is the most probable
   residue: V scored -22, I +30; the consensus column is blank from the second position on), and `mmseqs search qdb three_prof ...` died in the prefilter:
   `Index table: counting k-mers ... Segmentation fault`, `Error: Prefilter died`. The HHM layout
   MMseqs2 expects is not documented in its help and its source is not on the pod; I did not pursue it.
2. **Alignments sampled from the HMMs, then `convertmsa` + `msa2profile`: used.** `build_emit.sh`:
   `gzip -dc Pfam-A.hmm.gz | hmmemit -a -N 200 --seed 1 - | python pfam_models.py sto Pfam-A.hmm.gz >
   emit.sto` (200 sequences emitted by each model's core model, aligned on its match states; the filter adds
   `#=GF ID`/`#=GF AC` lines and rewrites insert columns' `.` as `-` and lower case as upper case, because
   MMseqs2 reads the rows as FASTA, where `.` is not a gap: without it the 121-state PF26733.1 became a
   profile of about 266 columns in the smoke test); `mmseqs convertmsa emit.sto msa --identifier-field 1`;
   `mmseqs msa2profile msa prof --match-mode 1 --match-ratio 0.5 --threads 3` (the settings of MMseqs2's own
   Pfam-A.full download). Result: 30,134 profiles (`prof.index`), sha256 of `prof` `908c33d6...4ffb`.
   Profile length against the model's LENG (from `profile2consensus`): 1,536 equal, 28,576 longer,
   22 shorter, mean absolute difference 4.76 positions: insert columns that more than half of the sampled
   sequences use became profile columns. On 60 sequences emitted from 3 models (`hmmemit -p -N 20 --seed 2`),
   16 were found at `-s 7.5 -e 100`.

Cost of route 2 (unit `m09a-build-emit200`, 3 threads; `data/prof_emit200/build-*.json`): hmmemit and filter
6.72 CPU-min (360 s), `convertmsa` 0.49 (29 s), `msa2profile` 16.63 (335 s, peak RSS 3.5 GiB); total
**23.86** CPU-min, 726 s. Once per Pfam release.

## Searches

All with `mm_search.sh` (`mmseqs search <q> <prof> aln tmp -s <s> -e 100 --max-seqs 300 --threads 4`, then
`mmseqs convertalis ... --format-output query,target,qstart,qend,tstart,tend,evalue,bits`), proteins as
queries. P1's runs used `mm_shards.sh` (two consecutive halves of the query DB's keys, `mmseqs createsubdb`
for sequences and headers, one search each, hits concatenated) so that no process could pass 300 CPU-min.
`mmseqs createdb proteins.fa qdb` (unit `m09a-createdb-P1`): 0.33 CPU-min, 22 s.

| unit | input | setting | shards: CPU-min | CPU-min | wall s | peak RSS MiB | hits |
|---|---|---|---|---|---|---|---|
| `m09a-P1-s5.7` | P1, 918,765 proteins | `-s 5.7` | 11.57, 11.16 | 22.80 | 354 | 1,720 | 985,636 |
| `m09a-P1-s7.5` | P1 | `-s 7.5` | 40.47, 38.95 | 79.49 | 1,207 | 3,973 | 1,433,040 |
| `m09a-sub-seq-s5.7` | P1 shard-0 first 20,000 (6,396,249 res.) | `-s 5.7` | - | 2.06 | 33 | 1,536 | 26,360 |
| `m09a-sub-seq-s7.5` | same 20,000 | `-s 7.5` | - | 7.27 | 112 | 3,755 | 44,475 |
| `m09a-sub-seq-s7.5-m1000` | same 20,000 | `-s 7.5 --max-seqs 1000` | - | 8.09 | 125 | 3,751 | 44,885 |
| `m09a-sub-prof-s5.7` | same 20,000 | profiles as queries, `-s 5.7 --max-seqs 20000` | - | 0.85 | 14 | 659 | 34,439 |
| `m09a-sub2k-s7.5` | P1 shard-0 first 2,000 (656,149 res.) | `-s 7.5` | - | 5.62 | 87 | 3,760 | 4,560 |
| `m09a-sub2k-exh` | same 2,000 | `-s 7.5 --exhaustive-search` | - | 2.07 | 32 | 623 | 5,721 |

MMseqs2 reports one alignment per (protein, profile) (no `--alt-ali`), so hits = pairs. The search cost on
the full stage is far below the subsets' per-residue rate (P1 `-s 5.7`: 0.077 CPU-min per million residues;
20,000-protein subset: 0.32, *derived*): building the profile k-mer index (about 25 s at 4 threads at
`-s 5.7`, unit log) is paid once per search, so small inputs are dominated by it.

Subset diagnostics (compare.py on the subsets, `summarize.py`):

| run | truth pairs | recall, any hit | best F1 (bits >=, P, R) | R at P >= 0.99 |
|---|---|---|---|---|
| sub-seq-s5.7 | 31,505 | 56.1% | 0.694 (32, 95.3%, 54.5%) | 51.1% |
| sub-seq-s7.5 | 31,505 | 75.8% | 0.815 (32, 93.7%, 72.0%) | 63.8% |
| sub-seq-s7.5-m1000 | 31,505 | 75.9% | 0.815 (32, 93.6%, 72.1%) | 63.8% |
| sub-prof-s5.7 (profiles as queries) | 31,505 | 51.9% | 0.659 (32, 95.3%, 50.4%) | 47.1% |
| sub2k-s7.5 | 3,042 | 76.2% | 0.823 (32, 94.4%, 73.0%) | 66.9% |
| sub2k-exh (MMseqs2 exhaustive mode) | 3,042 | 72.8% | 0.802 (32, 94.5%, 69.6%) | 64.0% |

So: `--max-seqs 300` is not limiting; proteins as queries beat profiles as queries; and MMseqs2's
`--exhaustive-search` (which for a profile target ran the profiles as queries against all 2,000 proteins,
`prefilter ... -k 0 --max-seqs 2000`, then `align`) is not better than `-s 7.5`. The 20,000-subset numbers
match the full stage's within 1 point, so the full-stage table is not a sampling accident.

Where `-s 7.5` misses strong hits (20,000-protein subset, pairs >= 50 bits above GA, `missed.py`): 1,137 of
15,346 missed (recall 92.6%), in 435 of 3,895 models; 102 models lost every such pair. The most-missed
models are large membrane and regulator families: LysR_substrate (49 of 121 missed), MFS_1 (48 of 160),
BPD_transp_1 (46 of 75), EamA (37 of 43), Aminotran_1_2 (30 of 77), AA_permease_2 (23 of 30), ADH_zinc_N
(21 of 25). Recall does not depend on model length (91% to 96% in every length bin).

## Comparison on stage P1

Definitions (`compare.py`). A truth pair is a (protein, model) with at least one producer row; its margin is
`seq_bitscore - GA_seq` (the sequence cutoff that admits it). An MMseqs2 pair's score is its bits. At threshold
t, pairs with bits >= t are called; recall = called truth pairs / 1,447,450; precision = called truth pairs /
called pairs. The curve runs over every distinct bits value; with `-e 100` MMseqs2 reports nothing below
about 27 bits, so the curve's loose end is the `-e 100` cut.

**Curves** (bits >= t: recall, precision):

| t | 25 | 30 | 32 | 35 | 40 | 50 | 60 | 80 | 100 |
|---|---|---|---|---|---|---|---|---|---|
| `-s 5.7` | 55.4, 82.2 | 54.7, 91.8 | 54.0, 95.2 | 52.3, 97.9 | 49.2, 99.3 | 42.2, 99.6 | 36.0, 99.7 | 25.1, 99.8 | 18.0, 99.8 |
| `-s 7.5` | 74.9, 77.0 | 73.4, 89.5 | 71.8, 93.8 | 68.4, 97.4 | 62.6, 99.1 | 50.9, 99.6 | 41.5, 99.7 | 27.2, 99.8 | 18.7, 99.8 |

E-value <= e (recall, precision), MMseqs2's E-value against the 30,134-profile DB:

| e | 1e-10 | 1e-5 | 1e-3 | 1e-2 | 0.1 | 1 | 10 | 100 |
|---|---|---|---|---|---|---|---|---|
| `-s 5.7` | 34.8, 99.7 | 44.8, 99.6 | 49.0, 99.3 | 51.1, 98.7 | 53.1, 96.9 | 54.6, 92.8 | 55.2, 86.0 | 55.4, 81.3 |
| `-s 7.5` | 39.7, 99.7 | 55.1, 99.5 | 62.4, 99.1 | 66.2, 98.4 | 70.1, 96.3 | 73.1, 90.9 | 74.6, 82.0 | 75.0, 75.7 |

Best F1: 0.689 at bits >= 32 (`-s 5.7`), 0.814 at bits >= 32 (`-s 7.5`). Recall at precision >= 99%: 49.7%
(bits >= 39) and 62.6% (bits >= 40). At precision >= 95%: 54.0% and 71.1%.

**Recall by margin above GA** (truth pairs per bin; recall with any hit / bits >= 32 / E <= 1e-3):

| margin (bits) | truth pairs | `-s 5.7` | `-s 7.5` |
|---|---|---|---|
| [0, 1) | 8,585 | 17.5 / 8.7 / 1.0 | 32.4 / 14.4 / 1.2 |
| [1, 2) | 21,276 | 16.5 / 9.0 / 0.8 | 32.5 / 15.9 / 1.1 |
| [2, 5) | 71,345 | 17.8 / 11.2 / 1.3 | 36.0 / 20.3 / 1.8 |
| [5, 10) | 104,297 | 21.9 / 16.6 / 3.0 | 43.3 / 31.1 / 4.7 |
| [10, 20) | 168,755 | 28.5 / 25.9 / 10.7 | 53.0 / 47.2 / 17.8 |
| [20, 50) | 371,684 | 43.9 / 43.3 / 38.1 | 71.1 / 69.9 / 60.3 |
| >= 50 | 701,508 | 78.4 / 78.3 / 77.7 | 92.7 / 92.5 / 91.6 |

Hits near GA are the ones lost: at E <= 1e-3 almost none of the 101,206 pairs within 5 bits of GA are found.

**Boundaries** (called truth pairs; MMseqs2's alignment `[qstart-1, qend)` against the pair's truth domain
that overlaps it most): `-s 7.5`, 1,084,987 pairs: median absolute difference 1 residue at the start and 1 at
the end against the envelope (90th percentile 13 and 21), 0 and 1 against the HMMER alignment coordinates;
signed medians 0 and 0; median length ratio MMseqs2 / envelope 0.985; 99.9% overlap the domain. `-s 5.7`:
the same medians, 90th percentiles 11 and 17. Boundaries are not the problem. Because MMseqs2 gives one
alignment per pair, of the 1,552,653 truth domain rows only 70.7% (`-s 7.5`; 52.2% at `-s 5.7`) are
overlapped by their pair's alignment: repeats (67,534 pairs with more than one domain) lose their extra copies
unless `--alt-ali` is set (not measured).

**Calls that are not truth pairs.** At bits >= 32, `-s 7.5`: 68,231 such pairs; 58,110 (85%) overlap a truth
domain of another model on the same protein (most likely a family of the same clan, or the same domain
scored by a related model); 3,834 are on proteins with no truth domain at all. `-s 5.7`: 39,727, 34,252
(86%), 2,006.

**Clans: not measured.** Clan membership needs `Pfam-A.clans.tsv` (or the CL lines of `Pfam-A.hmm.dat`),
not on the pod and not in `refdb.yaml`. As a family-agnostic stand-in for clan collapsing: a truth domain is
covered when an MMseqs2 alignment of any model on its protein overlaps it by at least half of the shorter of
the two; an alignment is on a truth region when it so overlaps any truth domain of its protein.

| | covered truth domains: any hit / bits >= 32 / E <= 1e-3 | alignments on a truth region, same thresholds |
|---|---|---|
| `-s 5.7` | 65.5% / 63.4% / 58.8% | 96.2% / 99.2% / 99.9% |
| `-s 7.5` | 83.2% / 80.3% / 73.2% | 94.8% / 98.9% / 99.9% |

This is an upper bound on clan-level recall (a hit of an unrelated model on the region counts too).

## Prefilter option (*derived* from the `-s 7.5` and `-s 5.7` tables, not run)

(i) **By protein**: proteins with any MMseqs2 hit at E <= e go to `hmmsearch` against all models; the rest are
not scanned. Work left = their residues / 294,516,869 (`hmmsearch` time is close to linear in residues at a
fixed model set: the producer's shards of 32.9 M residues took 244 to 248 CPU-min each, 31.5 M took 221);
recall ceiling = truth pairs on those proteins / 1,447,450 (`hmmsearch` on a kept protein reproduces its
truth rows exactly, since at Z = 1e8 a protein's hits do not depend on the other proteins, MR-12).

| E <= | `-s 7.5`: proteins kept, work left, ceiling | `-s 5.7`: proteins kept, work left, ceiling |
|---|---|---|
| 100 | 689,139, 83.3%, 93.1% | 572,858, 71.6%, 80.8% |
| 10 | 679,599, 82.2%, 92.7% | 565,648, 70.7%, 80.3% |
| 1 | 666,272, 80.7%, 92.0% | 556,350, 69.7%, 79.4% |
| 0.1 | 652,140, 79.2%, 90.7% | 546,731, 68.6%, 78.3% |
| 1e-3 | 619,154, 75.8%, 87.0% | 524,867, 66.1%, 75.3% |

Total at `-s 7.5`, E <= 100: 0.833 x 2,191.21 + 79.82 = 1,905 CPU-min, 13% less than `hmmsearch` alone, for
6.9% fewer truth pairs. A perfect protein filter could save at most 11.2% of residues here (proteins with a
truth domain hold 88.8%), so on complete prokaryote genomes this option cannot save much.

(ii) **By (protein, model) pair**: only the pairs MMseqs2 calls go to `hmmsearch`. The recall ceiling is
MMseqs2's own pair recall (75.0% at `-s 7.5`, E <= 100; 1,433,040 pairs kept). The dynamic-programming work
is tiny: the sum of protein length x model length over the kept pairs is 0.0079% of 294,516,869 residues x
4,754,065 match states (`pair_work_fraction` in `compare.json`). But `hmmsearch` as the producer runs it (one
model against a block of sequences) would need a different target set per model, so its per-call overhead,
not the DP cells, would set the cost; that is not derivable from these tables. The recall ceiling alone rules
it out as a lossless prefilter.

## What could not be measured, and why

- The standard MMseqs2 Pfam database (profiles from Pfam-A.full alignments): needs `mmseqs databases` or
  Pfam-A.full / Pfam-A.seed, a download; not in `refdb.yaml` (AGENTS §8, the card's limits). Its recall
  could differ from the HMM-sampled profiles measured here, in either direction.
- HMMER3 to MMseqs2 via HH-suite HHM: `convertprofiledb` accepted the converted file but produced an unusable
  profile and the prefilter crashed (above). No HH-suite tool (`hhmake`, `hhconvert`) is in the image.
- Clan-collapsed comparison: no clan file on the pod (the any-model region stand-in is above).
- Repeats: one alignment per pair; `--alt-ali` not run.
- The prefilter options were derived, not run (as the card asks).
- Wall time is not compared as a rate: both scans shared the pod with 40 to 49 other CPUs of the account.

## Commands, in order (all from the pod host; `L=/tmp/if2453/m09a/bin/launch.sh`, `B=/tmp/if2453/m09a/bin`,
`D=/tmp/if2453/m09a/data`, `W=/tmp/if2453/m09a/work`, `PF=/stor/scratch/Wilke/if2453/all-the-data/l0/refdb/pfam/38.2/Pfam-A.hmm.gz`)

```
M09A_THREADS=1 bash $L models python $B/pfam_models.py models $PF $D/models.tsv
M09A_THREADS=3 bash $L build-emit200 sh $B/build_emit.sh $PF $D/prof_emit200 200 1
M09A_THREADS=1 bash $L export-P1 python $B/export_stage.py /tmp/if2453/r10/root P1 $D/P1
M09A_THREADS=1 bash $L createdb-P1 mmseqs createdb $D/P1/proteins.fa $D/P1/qdb
M09A_THREADS=1 bash $L export-P1sub python $B/export_stage.py /tmp/if2453/r10/root P1 $D/P1sub --shard 0 --head 20000
M09A_THREADS=1 bash $L createdb-P1sub mmseqs createdb $D/P1sub/proteins.fa $D/P1sub/qdb
bash $L sub-seq-s5.7 sh $B/mm_search.sh $D/P1sub/qdb $D/prof_emit200/prof $W/sub-seq-s5.7 5.7 100 seq 300
bash $L sub-seq-s7.5 sh $B/mm_search.sh $D/P1sub/qdb $D/prof_emit200/prof $W/sub-seq-s7.5 7.5 100 seq 300
bash $L sub-prof-s5.7 sh $B/mm_search.sh $D/P1sub/qdb $D/prof_emit200/prof $W/sub-prof-s5.7 5.7 100 prof 20000
bash $L sub-seq-s7.5-m1000 sh $B/mm_search.sh $D/P1sub/qdb $D/prof_emit200/prof $W/sub-seq-s7.5-m1000 7.5 100 seq 1000
bash $L missed-sub-seq-s7.5 python $B/missed.py $D/P1sub $D/models.tsv $W/sub-seq-s7.5/hits.tsv
bash $L P1-s5.7 sh $B/mm_shards.sh $D/P1/qdb $D/prof_emit200/prof $W/P1-s5.7 5.7 100 2
bash $L P1-s7.5 sh $B/mm_shards.sh $D/P1/qdb $D/prof_emit200/prof $W/P1-s7.5 7.5 100 2
M09A_THREADS=1 bash $L export-P1sub2k sh -c "python $B/export_stage.py /tmp/if2453/r10/root P1 $D/P1sub2k --shard 0 --head 2000 && mmseqs createdb $D/P1sub2k/proteins.fa $D/P1sub2k/qdb -v 1"
bash $L sub2k-s7.5 sh $B/mm_search.sh $D/P1sub2k/qdb $D/prof_emit200/prof $W/sub2k-s7.5 7.5 100 seq 300
M09A_EXTRA=--exhaustive-search bash $L sub2k-exh sh $B/mm_search.sh $D/P1sub2k/qdb $D/prof_emit200/prof $W/sub2k-exh 7.5 100 seq 300
# for each run R and its export E:
bash $L cmp2-R python $B/compare.py $D/E $D/models.tsv $W/R/hits.tsv $W/R/compare.json --label R
python3 $B/summarize.py $W/P1-s5.7/compare.json $W/P1-s7.5/compare.json
python3 $B/costs.py /tmp/if2453/m09a /tmp/if2453/r10/root P1
```

Every run record is `/tmp/if2453/m09a/runs/m09a-<name>.json` (and `data/prof_emit200/build-*.json`,
`work/P1-s*/shard-*.json`); every unit log `/tmp/if2453/m09a/logs/m09a-<name>.log`; every comparison
`/tmp/if2453/m09a/work/<run>/compare.json`.
