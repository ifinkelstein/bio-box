# Scanning cluster representatives only: saving and inference error (card 0-09b, wilkcomp02, 2026-11-11)

Card 0-09 (**H**, PI decision) option (b): scan only the NR90 or NR50 representatives with Pfam and infer
the domains of the other members. This file measures what that saves and how often the inferred call is
wrong. Every number is *measured* on wilkcomp02 on 2026-11-11 between 07:03 and 07:20 CST by the commands at
the end, unless labelled *derived* (formula given) or *estimate*. "Measured" here means counted by
`representatives.py` from the producer's tables; the Pfam calls themselves are the producer's `hmmsearch`
rows (Pfam 38.2, gathering cutoff), not rerun. No code under `src/` or `tests/` changed; nothing was written
under `/tmp/if2453/r10` or `/stor`; nothing was downloaded.

## Summary

Release `rN.1b` (stages S1, 1,134,550 INPHARED phage proteins, and P1, 918,765 proteins of 250 complete
prokaryotes, one per GTDB species; 2,053,075 distinct proteins, 566,303,726 residues, 2,031,144 true
(protein, model) pairs). "Right" and "missed" are fractions of the true pairs (the producer's own calls on
every protein); "wrong" is the fraction of reported pairs (scanned plus inferred) that the protein's own scan
does not have. Hybrid: representatives are scanned, and so is every member whose length differs from its
representative's by more than 20%. CPU-minutes are *derived*: residues not scanned / 1e6 x 7.44 (0-09a,
`2026-11-11_mmseqs_profile_vs_hmmsearch.md`, P1) or x 7.96 (`2026-11-10_build_stage_pod.md`).

| level | scope | reps only: residues scanned, right, wrong, missed | CPU-min saved (7.44 / 7.96) | hybrid: residues scanned, right, wrong, missed | CPU-min saved (7.44 / 7.96) |
|---|---|---|---|---|---|
| NR90 | S1 | 41.7%, 98.90%, 4.27%, 1.10% | 1,179 / 1,262 | 45.1%, 99.09%, 0.98%, 0.91% | 1,111 / 1,188 |
| NR90 | P1 | 94.0%, 99.88%, 0.41%, 0.12% | 132 / 141 | 94.3%, 99.89%, 0.12%, 0.11% | 126 / 134 |
| NR90 | S1+P1 | 68.9%, 99.60%, 1.54%, 0.40% | 1,311 / 1,403 | 70.7%, 99.66%, 0.37%, 0.34% | 1,236 / 1,323 |
| NR50 | S1 | 21.7%, 94.08%, 11.68%, 5.92% | 1,583 / 1,693 | 32.5%, 95.44%, 3.80%, 4.56% | 1,365 / 1,461 |
| NR50 | P1 | 49.9%, 96.23%, 6.62%, 3.77% | 1,097 / 1,174 | 55.8%, 96.81%, 2.88%, 3.19% | 969 / 1,037 |
| NR50 | S1+P1 | 36.4%, 95.61%, 8.11%, 4.39% | 2,680 / 2,867 | 44.6%, 96.41%, 3.14%, 3.59% | 2,335 / 2,498 |

At 10 Gbp (`rN.1`, clusters only; P2's Pfam is not published, so no error is measured there): NR90
representatives hold 78.1% of the residues and NR50 representatives 26.6%; scanning only them would save
5,446 / 5,827 (NR90) and 18,240 / 19,515 (NR50) CPU-minutes of the 24,855 / 26,592 that a full scan costs
(*derived*, 3,340.67 M residues x 7.44 or 7.96). The hybrid scans 78.9% (NR90) and 41.2% (NR50).

Where inference fails: a member more than 20% shorter than its representative has its representative's set
in 68.8% (NR90) and 64.3% (NR50) of cases, and 45.3% (NR90) and 37.5% (NR50) of the pairs inferred for it are
wrong; members within 20% of the representative's length get 97.8% (NR90) and 89.9% (NR50) equal sets, with
1.69% and 5.71% of inferred pairs wrong. Where the sets agree, the representative's envelopes mapped through
the stored member alignment land within 2 residues of the member's own at the 95th percentile (NR90; 8 and
12 at NR50); copied without mapping, 23 and 33 (NR90), 80 and 90 (NR50).

**Reading (my interpretation, not a measurement).** NR90 inference is nearly free of error but saves
little where genomes are one per species: 6% of P1's residues in `rN.1b`, about 18% of P1's and P2's at
10 Gbp. It saves 58% of S1's, because INPHARED keeps many isolates of one phage species: S1 is the redundant
case here. NR50 saves two to three times more (64% of S1+P1's residues, 73% at 10 Gbp) but errs in both
directions: 8.1% of reported pairs are false and 4.4% of true pairs are missed. The hybrid rule is worth its
cost: rescanning length-mismatched members adds 1.8 points of residues at NR90 and 8.2 at NR50,
and cuts false pairs by a factor of 4 (1.54% to 0.37%) and 2.6 (8.1% to 3.1%). What it does not fix is the
misses among full-length members (1.6% of their true pairs at NR90, 6.5% at NR50), where the member has a
domain the representative lacks or the reverse; I expect these are mostly calls close to the gathering
cutoff that flip between homologs, but that was not measured (the score margin of the missed and false
pairs is a short follow-up on the same inputs). If option (b) is taken, NR90 plus the hybrid rule loses about 0.3% of
pairs in each direction for a saving that grows with redundancy; NR50 plus the hybrid is a lower tier
(3.1% false, 3.6% missed) that should be labelled as inferred, as card 0-09 proposes. Two costs outside this
measurement: scanning by representative means the release's clustering must precede the scan (today the
scan is per stage and the clustering per release), and representatives change when a stage is added (in
`rN.1`, 576,327 P1 proteins have an NR50 representative from another stage), so an incremental build must
scan new representatives that were members before.

Strain collections (*estimate*, not measured): exact duplicates are already scanned once, because a protein
is keyed by its sequence; a strain collection's distinct proteins are then mostly allelic variants at 95% to
100% identity and full length, so the non-representative fraction at NR90 would be far above P1's 6% (S1, at
53% of proteins, is closer to that case), and the error would follow the high-identity and
within-20%-length rows here (at NR90, 95.8% to 96.6% equal sets at identity 0.95 to 0.999, 97.8% within 20%
of length). Draft assemblies add genes truncated at contig ends, which fall in the shorter-than-20% class
that the hybrid rule rescans.

## Setup

- Host wilkcomp02 (112 CPUs). Image `/stor/home/if2453/bio-box/images/bio-box-test_74b3e05.sif`
  (`sha256:c30be58557d756095346308c58cc3098af4bbe66404f67ab87be3f64f147c898`, the digest in the producer's run
  records), Python 3.14, pyarrow, DuckDB 2.0.0; `PYTHONPATH=/tmp/if2453/ci/r10/src:/tmp/if2453/pilot_s1/pydeps`
  (the tree that built the root, commit `740d049`).
- Data root `/tmp/if2453/r10/root`, read only. While this ran, the root's phase-2 driver (`r10-p2-driver3`,
  `r10-p2-P2-7`) kept running; nothing of it was touched.
- Scripts: `scripts/stress/annotation_cost/representatives.py` (the measurement), `rep_tables.py` (prints the
  tables below from its JSON), and `launch.sh` (0-09a's launcher, now taking `M09_CARD=m09b` for this card's
  scratch directory, unit prefix and CPU gate; card 0-09a's behaviour is unchanged). Mirrored to
  `/tmp/if2453/m09b/bin/`. Scratch `/tmp/if2453/m09b/{bin,work,runs,logs}`.
- Every job ran through `launch.sh`: unit `m09b-<name>` (`systemd-run --user --collect -p MemoryMax=32G
  -p MemorySwapMax=0 -p CPUQuota=200%`), `nice -n 19 singularity exec --cleanenv ... python timed.py
  runs/m09b-<name>.json -- <command>`, 2 threads each, at most two jobs (4 threads) at once. Before every
  launch it summed `ps -u if2453 -o pcpu=` (the `ps` line left out) and started only if that plus the job's
  threads x 100 was at most 4000% between 08:00 and 20:00 CST on weekdays and 5600% otherwise; all launches
  were between 07:03 and 07:15 CST (cap 5600%). Readings in `/tmp/if2453/m09b/logs/launch.log`: 3,936%,
  4,283%, 4,412%, 4,261%, 4,246%; one launch waited (07:11, 6,142%) and started at the next poll. The table
  renders and the schema inspection ran in `systemd-run --user --scope -p MemoryMax=32G -p MemorySwapMax=0
  nice -n 19 singularity exec ...`, 1 thread.
- Cost: `m09b-infer-rN.1b-3` 2.07 CPU-min, 93.8 s wall, peak RSS 4,440 MiB; `m09b-saving-rN.1-2` 3.19 CPU-min,
  118.1 s, 7,891 MiB (`runs/*.json`). No process came near 300 CPU-minutes.

## Inputs (the producer's file lists, each checked by sha256)

- Proteins: the `protein_seq` files of the release manifest (`manifests/release=<R>.json`). `rN.1b`:
  stage S1 `.../stage=S1/protein_prefix=all/part-0.parquet` (1,134,550 rows, sha256 `2e617231...`) and P1
  `part-0.parquet` (918,765 rows, `c9268d41...683d`, the file 0-09a used); 240 proteins are in both stages, so
  2,053,075 distinct. `rN.1` adds P2 `part-0` and `part-1` (4,585,679 and 3,995,034 rows); 10,620,684
  distinct, 13,302 in two stages. Length is the sequence without `*` (what `domains.scan` scans); it equals
  the table's `length` column for every row, and no sequence has `*` before its last residue, so the
  clustering's alignment coordinates (made on the unstripped sequence) and the domain coordinates agree.
- Clusters: the manifest's `cluster_member` and `cluster_member_aln` files. `rN.1b`
  `l2/cluster_member/release=rN.1b/part-0.parquet` (`8bc4461d...c709`) and `l2/cluster_member_aln/...`
  (`709d9ba1...4fd1`), 4,106,150 rows each; `rN.1` `6764c023...94f0` and `c91ab1fc...ef4dd`, 21,241,368 rows
  each. Every protein is a member exactly once per level and every representative is its own member (checked).
  Clusters: `rN.1b` 1,395,712 (NR90) and 761,833 (NR50); `rN.1` 8,278,537 and 2,966,200, the counts card 0-09
  quotes.
- Domains: each stage's `l1_domains` file named in `runs/stage=<S>/domains-pfam-38.2.json` (exit 0): S1
  `l1/domains/library=pfam-38.2/stage=S1/part-0.parquet` (605,564 rows, `1973c6f2...d9ea`) and P1
  (1,552,653 rows, `fb4a545b...e539`). The 240 shared proteins were scanned in both stages and their rows agree
  (0 differ); distinct rows 2,158,015 on 1,187,061 proteins, 2,031,144 (protein, model) pairs.
- What NR90 and NR50 mean (`src/biobox/l2/cluster.py` at `740d049`): MMseqs2 `cluster --min-seq-id 0.9 -c 0.8
  --cov-mode 1 --cluster-mode 2` on all distinct proteins of the release's stages; NR50 is `--min-seq-id 0.5`
  with the same coverage on the NR90 representatives, merged so every protein maps to an NR50 representative.
  Coverage mode 1 requires the alignment to cover 80% of the member (not of the representative), and cluster
  mode 2 is greedy by length, so the representative is the longest sequence and a member may be a fragment of
  it. Clusters are made per release, over all its stages: a stage's member may have its representative in
  another stage (the "rep in another stage" column).
- Identity: `cluster_member` holds no identity; `cluster_member_aln` does. It is a re-alignment of every
  member to its representative (`mmseqs align -a -e inf --min-seq-id 0 -c 0 --alignment-mode 3`), with
  `pident` (MMseqs2's fident), the member's span `qstart, qend`, the representative's span `tstart, tend`
  (0-based half-open) and the CIGAR (query = representative). Every member has one (0 missing). Because it is
  a re-alignment with thresholds off, 13,854 NR90 members read below 0.9 on it.
- CIGAR convention, checked independently: on the first 5,000 non-representative members by uid at each
  level, walking the CIGAR over the two sequences with `I` consuming the representative and `D` the member
  ends exactly at both stored ends and reproduces the stored `pident` (identical residues in M columns /
  `alnlen`, within 0.0015) for 5,000 of 5,000 at NR90 and at NR50; with `I` and `D` swapped, 3,652 and 1,746
  of 5,000. All 1,948,605 non-representative CIGARs (657,363 at NR90, 1,291,242 at NR50) span exactly their
stored coordinates under that convention.

## Definitions (`representatives.py`)

- D(p): the distinct model accessions with an `l1_domains` row on protein p. Architecture: the accessions
  ordered by `env_start`, `env_end`, accession, repeats kept.
- For a member m (not its own representative) with representative r, inference reports D(r). Pairs right
  |D(m) & D(r)|, false |D(r) - D(m)|, missed |D(m) - D(r)|. Classes (exclusive): equal (both empty included);
  subset (D(m) a proper subset of D(r), D(m) empty included); superset (D(r) a proper subset, D(r) empty
  included); disjoint (both non-empty, nothing shared); overlap (the rest).
- Member tables: "false of inferred" = sum of false / sum of |D(r)| over members; "missed of true" = sum of
  missed / sum of |D(m)|. Per cluster: the same fractions within each cluster (pairs pooled in it), then the
  mean over clusters that have at least one member in the scope, so each cluster counts once.
- Decision table: scanned = the scope's representatives (plus, for the hybrid, members with
  |L_m - L_r| > 0.2 L_r); right = true pairs of scanned proteins + right pairs of inferred members; right +
  missed = the scope's true pairs. A scope is the proteins of one stage (a protein in two stages counts in
  both) or all; a member whose representative is in the other stage still counts as inferred.
- Length classes: shorter (L_m < 0.8 L_r), within 20%, longer (L_m > 1.2 L_r); the hybrid rescans the first
  and the last.

## 1. Saving

`rN.1b` (clusters over S1 and P1 together):

| level | scope | proteins | residues | representatives | not representatives: proteins | residues | members whose rep is in the other stage | hybrid: residues scanned |
|---|---|---|---|---|---|---|---|---|
| NR90 | S1 | 1,134,550 | 271,831,237 | 533,751 | 53.0% | 58.3% | 543 | 45.1% |
| NR90 | P1 | 918,765 | 294,516,869 | 862,131 | 6.2% | 6.0% | 966 | 94.3% |
| NR90 | all | 2,053,075 | 566,303,726 | 1,395,712 | 32.0% | 31.1% | 0 | 70.7% |
| NR50 | S1 | 1,134,550 | 271,831,237 | 294,149 | 74.1% | 78.3% | 8,903 | 32.5% |
| NR50 | P1 | 918,765 | 294,516,869 | 467,789 | 49.1% | 50.1% | 7,605 | 55.8% |
| NR50 | all | 2,053,075 | 566,303,726 | 761,833 | 62.9% | 63.6% | 0 | 44.6% |

`rN.1` (clusters over S1, P1 and P2; 10 Gbp):

| level | scope | proteins | residues | representatives | not representatives: proteins | residues | members whose rep is in another stage | hybrid: residues scanned |
|---|---|---|---|---|---|---|---|---|
| NR90 | S1 | 1,134,550 | 271,831,237 | 530,978 | 53.2% | 58.6% | 4,971 | 44.8% |
| NR90 | P1 | 918,765 | 294,516,869 | 758,431 | 17.5% | 17.5% | 138,516 | 83.1% |
| NR90 | P2 | 8,580,713 | 2,776,784,358 | 6,993,637 | 18.5% | 18.8% | 135,187 | 81.7% |
| NR90 | all | 10,620,684 | 3,340,668,639 | 8,278,537 | 22.1% | 21.9% | 0 | 78.9% |
| NR50 | S1 | 1,134,550 | 271,831,237 | 286,280 | 74.8% | 79.0% | 36,252 | 32.0% |
| NR50 | P1 | 918,765 | 294,516,869 | 267,465 | 70.9% | 71.9% | 576,327 | 42.7% |
| NR50 | P2 | 8,580,713 | 2,776,784,358 | 2,413,675 | 71.9% | 73.0% | 607,717 | 41.9% |
| NR50 | all | 10,620,684 | 3,340,668,639 | 2,966,200 | 72.1% | 73.4% | 0 | 41.2% |

CPU-minutes not spent at 10 Gbp (*derived*, residues not scanned / 1e6 x 7.44 or 7.96): NR90 reps only
732.0 M residues, 5,446 / 5,827; hybrid 705.2 M, 5,247 / 5,613. NR50 reps only 2,451.6 M, 18,240 / 19,515;
hybrid 1,965.4 M, 14,623 / 15,645. Per stage, `rN.1`: P2 alone NR90 3,890 / 4,162 and NR50 15,084 / 16,138.
For scale (*measured*, run records): the clustering itself cost 284.10 CPU-min for `rN.1b` and 2,195.91 for
`rN.1`; it is built for every release already.

How redundant these sets are: P1 and P2 are one complete genome per GTDB species (the ladder's selection,
`2026-11-11_r50_inputs.md`), so redundancy is low: P1 alone has 6.2% non-representative proteins at NR90, and
it rises to 17.5% once P2 (2,250 more species) joins the clustering, because P1 proteins then find
relatives at 90% identity among other species. S1 is not one
per species: INPHARED keeps every complete phage genome, including many isolates of one species, and 53% of
its proteins are NR90 members. The card's estimate for `rN.1` (22% and 72% of proteins, from the cluster
counts) is reproduced: 22.1% and 72.1%.

## 2. Inference error (`rN.1b`)

| level | scope | members | clusters | equal | subset | (member has none) | superset | (rep has none) | disjoint | overlap | architecture equal | false of inferred | missed of true | per cluster: equal, false, missed |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| NR90 | S1 | 600,799 | 153,532 | 95.5% | 3.44% | 1.71% | 0.93% | 0.46% | 0.04% | 0.07% | 95.0% | 6.67% | 1.76% | 97.3%, 4.45%, 2.22% |
| NR90 | P1 | 56,634 | 38,263 | 90.5% | 6.74% | 1.81% | 2.44% | 0.19% | 0.05% | 0.28% | 89.1% | 5.56% | 1.70% | 91.5%, 3.70%, 1.23% |
| NR90 | all | 657,363 | 191,405 | 95.1% | 3.72% | 1.72% | 1.06% | 0.43% | 0.04% | 0.09% | 94.5% | 6.43% | 1.75% | 96.1%, 4.15%, 1.84% |
| NR50 | S1 | 840,401 | 124,525 | 89.8% | 6.49% | 2.87% | 3.08% | 1.41% | 0.27% | 0.39% | 88.1% | 13.74% | 7.04% | 94.5%, 8.87%, 8.13% |
| NR50 | P1 | 450,976 | 118,051 | 77.4% | 12.77% | 1.48% | 7.71% | 0.54% | 0.10% | 1.98% | 74.1% | 11.33% | 6.59% | 79.9%, 7.76%, 5.33% |
| NR50 | all | 1,291,242 | 238,501 | 85.5% | 8.69% | 2.38% | 4.70% | 1.10% | 0.21% | 0.94% | 83.2% | 12.24% | 6.76% | 87.6%, 8.02%, 6.12% |

"Equal" includes members where neither has a domain: 49.3% of NR90 members overall. Among members where
either has a domain, the set is equal for 90.4% (S1), 89.7% (P1), 90.3% (all) at NR90 and 77.9%, 76.1%, 77.0%
at NR50. Where the sets are equal, the ordered architecture is also equal for 99.4% of NR90 members (all).

**By identity to the representative** (`cluster_member_aln.pident`; S1+P1; per stage in the JSON and in
`/tmp/if2453/m09b/work/tables.md`):

| level | pident | members | equal | subset | superset | disjoint or overlap | false of inferred | missed of true |
|---|---|---|---|---|---|---|---|---|
| NR90 | < 0.80 (re-alignment)* | 145 | 49.7% | 44.1% | 1.4% | 4.8% | 55.8% | 9.4% |
| NR90 | [0.80, 0.90) | 13,709 | 91.6% | 6.62% | 1.52% | 0.26% | 10.28% | 2.63% |
| NR90 | [0.90, 0.95) | 232,137 | 94.6% | 3.74% | 1.51% | 0.17% | 6.37% | 2.43% |
| NR90 | [0.95, 0.99) | 313,466 | 95.8% | 3.26% | 0.86% | 0.09% | 5.97% | 1.48% |
| NR90 | [0.99, 1.00) | 84,980 | 96.6% | 2.74% | 0.58% | 0.07% | 3.84% | 0.77% |
| NR90 | 1.00 | 12,926 | 81.5% | 17.52% | 0.62% | 0.39% | 37.10% | 2.64% |
| NR50 | < 0.50 (re-alignment) | 56,802 | 74.3% | 15.90% | 6.98% | 2.81% | 18.23% | 8.87% |
| NR50 | [0.50, 0.60) | 409,972 | 79.7% | 11.17% | 7.30% | 1.86% | 12.26% | 7.87% |
| NR50 | [0.60, 0.70) | 245,158 | 83.8% | 9.12% | 5.77% | 1.27% | 11.44% | 7.09% |
| NR50 | [0.70, 0.80) | 155,565 | 87.5% | 7.61% | 3.99% | 0.94% | 11.09% | 5.97% |
| NR50 | [0.80, 0.90) | 141,921 | 91.5% | 6.02% | 2.17% | 0.34% | 11.34% | 4.15% |
| NR50 | [0.90, 0.95) | 99,481 | 92.8% | 5.28% | 1.64% | 0.31% | 12.89% | 4.13% |
| NR50 | [0.95, 0.99) | 139,935 | 94.1% | 4.76% | 1.02% | 0.14% | 11.35% | 2.54% |
| NR50 | [0.99, 1.00) | 33,973 | 95.5% | 3.69% | 0.71% | 0.12% | 7.04% | 1.39% |
| NR50 | 1.00 | 8,435 | 81.6% | 16.83% | 0.98% | 0.58% | 41.77% | 5.11% |

\* *derived* by pooling the JSON's four NR90 bins below 0.80 ([0, 0.5) 12 members, [0.5, 0.6) 6, [0.6, 0.7) 11,
[0.7, 0.8) 116): classes 72, 64, 2 and 7 of 145 members; pairs false 110 of 197 inferred, missed 9 of 96 true.
Identity 1.00 means identical over the aligned part: these are mostly fragments of the representative, hence
the high false rate.

**By cluster size** (members, representative included; S1+P1):

| level | size | members | equal | subset | superset | disjoint or overlap | false of inferred | missed of true |
|---|---|---|---|---|---|---|---|---|
| NR90 | 2 | 99,904 | 96.2% | 2.53% | 1.19% | 0.09% | 4.20% | 1.70% |
| NR90 | 3-4 | 118,222 | 96.2% | 2.81% | 0.90% | 0.09% | 6.06% | 1.72% |
| NR90 | 5-9 | 142,768 | 96.2% | 2.76% | 0.90% | 0.10% | 5.70% | 1.73% |
| NR90 | 10-99 | 262,540 | 94.5% | 4.40% | 1.00% | 0.13% | 7.35% | 1.65% |
| NR90 | >= 100 | 33,929 | 87.9% | 9.25% | 2.38% | 0.50% | 8.67% | 2.37% |
| NR50 | 2 | 97,689 | 87.6% | 6.98% | 4.45% | 0.99% | 11.08% | 6.84% |
| NR50 | 3-4 | 156,919 | 87.6% | 6.98% | 4.50% | 0.95% | 10.70% | 6.87% |
| NR50 | 5-9 | 232,133 | 88.2% | 6.60% | 4.32% | 0.90% | 10.23% | 6.67% |
| NR50 | 10-99 | 643,295 | 86.1% | 8.17% | 4.60% | 1.18% | 12.11% | 6.65% |
| NR50 | >= 100 | 161,206 | 75.8% | 16.47% | 6.02% | 1.69% | 16.41% | 7.10% |

**By length against the representative** (fragments and fusions):

| level | scope | class | members | equal | subset | superset | disjoint or overlap | false of inferred | missed of true |
|---|---|---|---|---|---|---|---|---|---|
| NR90 | S1 | shorter by > 20% | 55,482 | 70.5% | 27.57% | 1.10% | 0.87% | 44.69% | 4.41% |
| NR90 | S1 | within 20% | 545,317 | 98.1% | 0.99% | 0.91% | 0.03% | 1.68% | 1.56% |
| NR90 | P1 | shorter by > 20% | 5,295 | 51.7% | 45.25% | 1.34% | 1.72% | 48.48% | 4.02% |
| NR90 | P1 | within 20% | 51,339 | 94.5% | 2.77% | 2.55% | 0.19% | 1.72% | 1.59% |
| NR90 | all | shorter by > 20% | 60,773 | 68.8% | 29.11% | 1.12% | 0.94% | 45.31% | 4.35% |
| NR90 | all | within 20% | 596,590 | 97.8% | 1.14% | 1.05% | 0.04% | 1.69% | 1.56% |
| NR50 | S1 | shorter by > 20% | 151,785 | 71.0% | 24.26% | 2.72% | 1.98% | 39.18% | 9.18% |
| NR50 | S1 | within 20% | 688,592 | 93.9% | 2.58% | 3.16% | 0.36% | 5.51% | 6.58% |
| NR50 | P1 | shorter by > 20% | 71,731 | 50.1% | 40.68% | 4.43% | 4.83% | 36.08% | 7.58% |
| NR50 | P1 | within 20% | 379,236 | 82.6% | 7.50% | 8.33% | 1.56% | 5.83% | 6.44% |
| NR50 | all | shorter by > 20% | 223,502 | 64.3% | 29.53% | 3.27% | 2.89% | 37.48% | 8.29% |
| NR50 | all | within 20% | 1,067,707 | 89.9% | 4.32% | 5.00% | 0.79% | 5.71% | 6.49% |
| NR50 | all | longer by > 20% | 33 | 87.9% | 3.03% | 9.09% | 0.00% | 8.33% | 12.00% |

No NR90 member is more than 20% longer than its representative (37 are longer at all; at NR50 7,078 are
longer, 33 by more than 20%), as cluster mode 2 makes the longest sequence the representative. Fusions
therefore show up as representatives with extra domains, which is the "subset" class of their shorter
members.

## 3. Coordinates where the sets agree (`rN.1b`)

Members whose set equals the representative's and is not empty; each member domain is paired with the
representative's domain of the same model whose envelope, mapped through the stored alignment, lies
closest. Residues: absolute start and end differences against the member's own envelope, mapped and copied
without mapping. Fractions: |s_m / L_m - s_r / L_r| (and the same for ends), on the same pairs.

| level | scope | domains | mapped start / end: median, p95 | both ends within 5 (mapped) | copied start / end: median, p95 | length fraction start / end: median, p95 |
|---|---|---|---|---|---|---|
| NR90 | S1 | 351,888 | 0, 2 / 0, 2 | 95.5% | 0, 24 / 0, 34 | 0, 0.080 / 0, 0.064 |
| NR90 | P1 | 93,905 | 0, 1 / 0, 2 | 95.6% | 0, 21 / 0, 27 | 0, 0.050 / 0, 0.036 |
| NR90 | all | 445,736 | 0, 2 / 0, 2 | 95.5% | 0, 23 / 0, 33 | 0, 0.073 / 0, 0.057 |
| NR50 | S1 | 397,809 | 0, 9 / 0, 14 | 86.1% | 1, 89 / 4, 103 | 0.006, 0.228 / 0.007, 0.216 |
| NR50 | P1 | 583,138 | 0, 7 / 0, 11 | 87.8% | 6, 74 / 8, 83 | 0.014, 0.169 / 0.012, 0.141 |
| NR50 | all | 980,854 | 0, 8 / 0, 12 | 87.1% | 4, 80 / 6, 90 | 0.011, 0.188 / 0.010, 0.168 |
| NR90 | all, member shorter by > 20% | 18,152 | 0, 3 / 0, 8 | 90.7% | 11, 326 / 55, 397 | 0.113, 0.538 / 0.072, 0.472 |
| NR50 | all, member shorter by > 20% | 95,347 | 0, 11 / 0, 19 | 84.5% | 32, 349 / 48, 363 | 0.138, 0.479 / 0.087, 0.456 |

The stored alignment is needed: copying the representative's coordinates, or scaling them by length, is off
by tens of residues at the 95th percentile, and by hundreds for fragments. With architecture also equal, the
mapped 95th percentiles are 1 / 2 (NR90) and 5 / 8 (NR50) residues.

## 4. The decision, in full

The counts behind the summary table (`rN.1b`; `infer-rN.1b.json`, `inference.<level>.scopes.<scope>.decision`).
The summary's fractions are *derived* from them: residues scanned / the scope's residues (section 1), right /
true, false / (right + false), missed / true.

| level | scope | mode | true pairs | members scanned too | residues scanned | pairs right | pairs false | pairs missed |
|---|---|---|---|---|---|---|---|---|
| NR90 | S1 | reps only | 583,892 | 0 | 113,323,368 | 577,442 | 25,773 | 6,450 |
| NR90 | S1 | hybrid | 583,892 | 55,482 | 122,547,487 | 578,586 | 5,752 | 5,306 |
| NR90 | P1 | reps only | 1,447,450 | 0 | 276,763,388 | 1,445,714 | 5,897 | 1,736 |
| NR90 | P1 | hybrid | 1,447,450 | 5,295 | 277,637,960 | 1,445,902 | 1,677 | 1,548 |
| NR90 | all | reps only | 2,031,144 | 0 | 390,054,194 | 2,022,959 | 31,669 | 8,185 |
| NR90 | all | hybrid | 2,031,144 | 60,773 | 400,152,460 | 2,024,291 | 7,429 | 6,853 |
| NR50 | S1 | reps only | 583,892 | 0 | 59,082,752 | 549,322 | 72,663 | 34,570 |
| NR50 | S1 | hybrid | 583,892 | 151,809 | 88,306,656 | 557,267 | 22,014 | 26,625 |
| NR50 | P1 | reps only | 1,447,450 | 0 | 147,083,697 | 1,392,907 | 98,741 | 54,543 |
| NR50 | P1 | hybrid | 1,447,450 | 71,740 | 164,242,013 | 1,401,222 | 41,541 | 46,228 |
| NR50 | all | reps only | 2,031,144 | 0 | 206,145,204 | 1,942,042 | 171,397 | 89,102 |
| NR50 | all | hybrid | 2,031,144 | 223,535 | 252,525,625 | 1,958,301 | 63,550 | 72,843 |

What the hybrid adds, in points of the scope's residues (*derived*, difference of the two modes' residues
scanned / the scope's residues): NR90 S1 3.39, P1 0.30, all 1.78; NR50 S1 10.75, P1 5.83, all 8.19. At NR90
over S1+P1 it leaves 596,590 inferred members with 97.8% equal sets, 1.69% of their inferred pairs false and
1.56% of their true pairs missed (`decision.hybrid.remaining_members`). For scale, the producer's full Pfam
scan of these two stages cost 2,042.37 (S1) + 2,191.21 (P1) = 4,233.58 CPU-min (*measured*, its shard run
records, summed in 0-09a); the NR90 reps-only saving of 1,311 to 1,403 CPU-min is 31% to 33% of it.

## What could not be measured, and why

- Inference error at 10 Gbp (`rN.1`): stage P2's Pfam is not published (its scan is running in the root),
  so only the saving is measured there. `representatives.py --domains` on `rN.1` will measure it once
  `runs/stage=P2/domains-pfam-38.2.json` exists.
- A strain collection: none is built; the expectation above is an *estimate*.
- Identity from the clustering's own alignment: not stored; `cluster_member_aln.pident` comes from a
  re-alignment with thresholds off, which is why a few NR90 members read below 0.9.
- Why full-length members miss or gain domains (score margin to the gathering cutoff of those pairs): not
  measured here; the producer's rows carry `seq_bitscore`, so it is a short addition to the same script.
- Clan-level agreement: no clan file on the pod (as in 0-09a); a member and representative that hit two
  families of one clan count as a false and a missed pair here.

## Commands, in order

Laptop: `scp scripts/stress/annotation_cost/{representatives.py,rep_tables.py,launch.sh,timed.py}
if2453@wilkcomp02.ccbb.utexas.edu:/tmp/if2453/m09b/bin/`. On the pod (`L=/tmp/if2453/m09b/bin/launch.sh`,
`B=/tmp/if2453/m09b/bin`, `W=/tmp/if2453/m09b/work`, `R=/tmp/if2453/r10/root`,
`S="systemd-run --user --scope --quiet -p MemoryMax=32G -p MemorySwapMax=0 nice -n 19 singularity exec
--cleanenv --bind /tmp/if2453 --env PYTHONPATH=/tmp/if2453/ci/r10/src:/tmp/if2453/pilot_s1/pydeps
/stor/home/if2453/bio-box/images/bio-box-test_74b3e05.sif"`):

```
$S python $B/peek_tables.py      # unit m09b-inspect-2: schemas and 4 rows of the four tables (throwaway, not kept)
M09_CARD=m09b M09A_THREADS=2 bash $L infer-rN.1b python $B/representatives.py $R rN.1b $W/infer-rN.1b.json --domains --threads 2
#   exit 1: DuckDB 2.0 has no Connection.fetch_record_batch; fixed (to_arrow_reader) and rerun:
M09_CARD=m09b M09A_THREADS=2 bash $L infer-rN.1b-2 python $B/representatives.py $R rN.1b $W/infer-rN.1b.json --domains --threads 2
M09_CARD=m09b M09A_THREADS=2 bash $L saving-rN.1 python $B/representatives.py $R rN.1 $W/saving-rN.1.json --threads 2
#   both exit 0; then the hybrid's residues were added to the saving step and both rerun (final):
M09_CARD=m09b M09A_THREADS=2 bash $L infer-rN.1b-3 python $B/representatives.py $R rN.1b $W/infer-rN.1b.json --domains --threads 2
M09_CARD=m09b M09A_THREADS=2 bash $L saving-rN.1-2 python $B/representatives.py $R rN.1 $W/saving-rN.1.json --threads 2
$S python $B/rep_tables.py $W/infer-rN.1b.json $W/saving-rN.1.json > $W/tables.md     # unit m09b-tables-2
```

The decision table of `infer-rN.1b-2` equals that of `infer-rN.1b-3` digit for digit (the script is
deterministic: inputs sorted, no clock or random input). Final outputs: `infer-rN.1b.json` sha256
`0db085f6...52f1`, `saving-rN.1.json` sha256 `3c46a0b8...13ce`; run records `/tmp/if2453/m09b/runs/m09b-*.json`,
unit logs `/tmp/if2453/m09b/logs/m09b-*.log`, gate readings `/tmp/if2453/m09b/logs/launch.log`.
