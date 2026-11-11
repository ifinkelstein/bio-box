# Second trial: a custom-HMM search with neighbourhoods on `rN.1`, end to end (2026-11-11)

Question: now that `rN.1` has a search pack (card W-16), does one `bio-box search hmm` command run a real
custom-HMM search with gene neighbourhoods against the 10 Gbp release, does it give the answer the first
trial made by hand, and what does it cost? This is a trial, not a task card: no code under `src/` or
`tests/` changed. Every number is *measured* on wilkcomp02, 2026-11-11, between 01:37 and 01:58 CST, by the
command named beside it, unless it is labelled *derived*.

## Result in one paragraph

Yes. `bio-box search hmm --hmm RNA_pol.hmm --release rN.1 --bits 30 --neighbors 5 --positive
NC_001604:3170-5822` ran in a 32G, 4-CPU `pod_run.sh` unit and ended `ok` in 534 s wall and about 12.6
CPU-minutes, with every output written. Its answer equals the hand-made one exactly: the same 1,957 proteins
at 30 bits or more (none in one set and not the other, every score equal to 0.1 bit), 1,931 in phages and 26
in prokaryotes, T7 gp1 recovered at 367.7 bits, the top ten NR50 neighbour clusters the same clusters with
the same locus counts, and all 2,581 locus tables equal to an independent SQL window. The cost is memory: the
search worker peaked at 26.3 GiB resident (28,203,290,624 B) because it holds the whole 3.34 G-residue target
in Python, and the runner itself peaked at 5.65 GiB building a set of every protein uid of the target to
check the positive control. That bears out the suspicion behind card W-19 (there is no W-19 card in the
repository yet). The search also re-verifies the whole 5.5 GB pack on every submit, and writes a scratch
directory into the shared pack directory while it does so.

## Setup

- Code: branch `trial/search-rN1-2` = `origin/plan/task-cards` at `950ffc2` with `origin/task/W-16` merged
  (`3a2bd7b`), mirrored to `/tmp/if2453/ci/trial2` by `sh scripts/pod_pytest.sh trial2 tests/r0/test_runbook.py`
  (`32 passed in 4.15s`).
- Data, read-only: root `/tmp/if2453/r10/root`, release `rN.1` (no domain library, no `member_cluster` or
  `contig_label` files); pack `/tmp/if2453/w16/packs/release=rN.1/71917a4167c7aebe1b7ff44d18a99306/` (W-16,
  5,510,724,073 B). Root bound read-only into the image in every command
  (`--bind /tmp/if2453 --bind /tmp/if2453/r10:/tmp/if2453/r10:ro`).
- Image `bio-box-test_74b3e05.sif`, sha256 `c30be58557d756095346308c58cc3098af4bbe66404f67ab87be3f64f147c898`,
  source bound in by `PYTHONPATH=/tmp/if2453/ci/trial2/src:/tmp/if2453/pilot_s1/pydeps` (RUNBOOK section 0).
- Profile `/tmp/if2453/trial2/profile.yaml`: `profiles/wilkcomp02.yaml` with `tmp: /tmp/if2453/trial2/tmp`.
- Query: the first trial's `RNA_pol.hmm` (Pfam PF00940.25, sha256 `32e6f0f0...0e87e`), copied to
  `/tmp/if2453/trial2/work/`. Reference answer: the first trial's `hmmsearch_all.tbl` (HMMER 3.4, Z = 1e8, -T 30).
- Pod: load average 42 to 48 throughout (the Pfam scan `r10-p2-*` of another build). Before each launch,
  `load.sh` summed `ps -u if2453 -o pcpu=` without the `ps` line: 4,494.6%, 4,969.6%, 4,414.0%, 4,450.8%,
  4,435.6% and 5,092.5% (01:38 to 01:58), each plus 400% under the night cap of 5,600%. With the load above 8
  the runner's thread cap was the day cap (40); the search asked for 4.
- Short commands ran through `/tmp/if2453/trial2/bb.sh` (or `py.sh` for a script): `/usr/bin/time -v` around
  `systemd-run --user --scope -p MemoryMax=32G -p MemorySwapMax=0 -p CPUQuota=400% nice -n 19 singularity exec ...`.
  The search ran through `scripts/pod_run.sh --mem 32G --threads 4` (nice 19 from the profile). A sampler unit
  (`trial2-sample`, nice 19, `/tmp/if2453/trial2/sample.sh`) wrote the unit's `MemoryCurrent`, `MemoryPeak`,
  `CPUUsageNSec`, cgroup `anon` and `file`, and `VmHWM`/`VmRSS` of every process in the unit, every 15 s, to
  `/tmp/if2453/trial2/logs/trial2-hmm.hwm` (1,670 lines, 01:39:50 to 01:48:53).

## 1. Finding the positive control (card W-18)

| Command | Last line | Cost |
|---|---|---|
| `bio-box find --root $ROOT --release rN.1 --accession NC_001604` | `NC_001604 719747a8fef97b2380955d92845861db 37369 39130 +1 f32d6c446a0450da9ad00862268a1baa eaa7ca5520ea5599998aab5d98a42b53` (49 genes, as the first trial's SQL) | 4.37 s wall, 4.42 s CPU, 333,100 kB max RSS |
| `bio-box find ... --accession NC_001604:3170-5822` | `NC_001604 719747a8fef97b2380955d92845861db 3170 5822 +1 39f9ef263e234b012d2671634468546e 2eef240afefddcd2c8ac7660eb50e384` | 4.49 s, 4.58 s, 331,088 kB |
| `bio-box find ... --protein 2eef240afefddcd2c8ac7660eb50e384` | `PP384410 72db861b0d1d12614d21bcdb80220dc1 3171 5823 +1 56cdfa046249671dec838cfa47976412 2eef240a...` (33 genes on 33 contigs carry T7 gp1's protein) | 4.14 s, 4.21 s, 330,064 kB |

Each printed on stderr `note: release rN.1 has no contig_label files (built before card B-27); reading
contig_occurrence instead (the slow path: every row group is scanned)`. On this release the slow path costs
about 4 s. The protein uid is the one the first trial found by hand.

## 2. The one command (card W-07)

```
bash /tmp/if2453/ci/trial2/scripts/pod_run.sh --profile /tmp/if2453/trial2/profile.yaml \
  --image /stor/home/if2453/bio-box/images/bio-box-test_74b3e05.sif --digest sha256:c30be585...c898 \
  --mem 32G --threads 4 --bind /tmp/if2453 --bind /tmp/if2453/r10:/tmp/if2453/r10:ro \
  --env PYTHONPATH=/tmp/if2453/ci/trial2/src:/tmp/if2453/pilot_s1/pydeps --env OMP_NUM_THREADS=1 \
  --log /tmp/if2453/trial2/logs/trial2-hmm.log trial2-hmm search hmm --root /tmp/if2453/r10/root --release rN.1 \
  --hmm /tmp/if2453/trial2/work/RNA_pol.hmm --bits 30 --neighbors 5 --positive NC_001604:3170-5822 \
  --out /tmp/if2453/trial2/out/rnap --profile /tmp/if2453/trial2/profile.yaml --pack-dir /tmp/if2453/w16/packs \
  --cpus 4 --mem-gb 30
```

Launched 01:39:48. The log (`/tmp/if2453/trial2/logs/trial2-hmm.log`), in full apart from repeated notes:

```
note: release rN.1 has no contig_label files (built before card B-27); reading contig_occurrence instead (the slow path: every row group is scanned)
search run: 87021d87e0f41e453c14d35763d51d17 running; pending record /tmp/if2453/trial2/tmp/search/if2453/runs/pending/run=87021d87e0f41e453c14d35763d51d17.json
note: release rN.1 has no member_cluster files (built before card B-27); reading cluster_member instead (the slow path: every row group is scanned)
... (9 slow-path notes in all, one or two per export)
87021d87e0f41e453c14d35763d51d17	ok
recipe	/tmp/if2453/trial2/out/rnap/recipe.yaml
hits	/tmp/if2453/trial2/out/rnap/hits.tsv	2581 genes of 1957 proteins
fasta	/tmp/if2453/trial2/out/rnap/hits.faa
neighbors	/tmp/if2453/trial2/out/rnap/neighbors.tsv	levels nr50,nr90, window 5g
loci	/tmp/if2453/trial2/out/rnap/loci	2581 tables
note	release rN.1 has no domain library: the neighbour summary is at the cluster levels nr50,nr90 only, and pfam_architecture is -
pod_run: trial2-hmm usage 
pod_run: trial2-hmm exited 0
```

`search status` from another shell: at 01:43 `{"run_id": "87021d87e0f41e453c14d35763d51d17", "status":
"running"}` (3.95 s wall, 328,016 kB); after the close `{"run_id": "87021d87e0f41e453c14d35763d51d17",
"status": "ok", "error": null, "rounds": 1, "final_round": 1, "stop_reason": "max_rounds"}`. Run record
(`run.parquet`): `controls_json` positive `recovered: true` for `2eef240a...`, negative control
(`shuffled_seed`) 0 hits, pass; `warnings_json` `[]`; `target_size` 10,620,684.

### Cost

| Search | Value | Source |
|---|---|---|
| Wall time, whole command | 534 s (8 min 54 s) | *derived*: sampler's first sample 01:39:50 to the log's last write 01:48:44 (`stat`) |
| Wall time, the run | 341.7 s | `run.parquet` `wall_s` (pending record to close) |
| CPU, whole unit | 12.64 CPU-minutes at least (758,315,688,000 ns) | `CPUUsageNSec` at the last sample, 01:48:38, 6 s before the end. The unit's own totals line is empty (rough edge 6) |
| CPU, the worker | 9.00 CPU-minutes | `run.parquet` `cpu_min_max_pid` |
| Peak memory, worker | 27,542,276 kB = 28,203,290,624 B = 26.27 GiB | `VmHWM` of the worker (`python -c import json, os, sys`), equal to `run.parquet` `peak_rss_bytes` |
| Peak memory, runner (`search hmm` process) | 5,920,924 kB = 5.65 GiB | `VmHWM`, reached between 01:41:51 and 01:42:21 |
| Unit memory | `MemoryPeak` 31,688,663,040 B of the 34,359,738,368 B cap; cgroup `anon` at most 29,148,450,816 B | sampler; `memory.events` `max 0`, `oom_kill 0` at 01:47 |

Phases, from the samples (*derived* from the times and `CPUUsageNSec` deltas):

| From | To | What | CPU |
|---|---|---|---|
| 01:39:50 | about 01:41:40 | positive lookup (about 4 s) and `pack.verify_dir` at submit: sha256 of the 5.5 GB pack, record counts of every shard, and the distinct proteins of the release through a sorter; RSS 0.3 GB | about 1 CPU |
| 01:41:40 | 01:42:30 | `resolve` builds `target_ids`, a Python set of all 10.6 M protein uids, by decompressing both `unique` shards; runner RSS to 5.65 GiB | about 1 CPU |
| 01:42:30 | | pending record written: from here `search status` sees the run (160 s after launch) | |
| 01:42:30 | about 01:42:45 | `_stage` copies the two `unique` shards (2,147,246,641 B) into the work directory | |
| about 01:42:45 | about 01:46:07 | worker `engine.load_target`: a dict of 10.6 M `str`, then a `TextSequenceBlock`, then its digitized copy; RSS 4.3 GB to 27.5 GB | 1 CPU |
| 01:46:07 | 01:48:08 | `hmmsearch` (pyhmmer) of the model and of the shuffled negative control over the target | 330 s in 121 s, 2.7 CPUs |
| 01:48:12 | 01:48:44 | run closed (`run.parquet`), then the four exports in the same process: `hits.tsv` +4.0 s, `hits.faa` +6.3 s, `neighbors.tsv` +11.5 s, `loci/` +10.0 s (file times) | |

The worker held 8.4 bytes of memory per target residue (*derived*: 28,203,290,624 B / 3,340,668,639 aa in
`unique`). The first trial planned this search in a 16G scope, which would have killed it. With `--mem-gb 30`
under an inherited 32 GiB cap the runner also sets RLIMIT_AS at 30 GB (runner module docstring, step 3), so
the margin was about 1.8 GB (*derived*: 30 × 10^9 - 28,203,290,624 B of resident memory). The address space,
which is larger than resident memory, is not measured here.

### The exports run separately (RUNBOOK section 7)

Each in its own capped scope (`bb.sh`), on run `87021d87...`, writing into `/tmp/if2453/trial2/out/sep/`:

| Command | Last line | Wall | CPU (user + system) | Max RSS |
|---|---|---|---|---|
| `search hits $RUN` | `protein 55bf72b132489fa724396089834ae0a1 375.85125732421875 364.69476318359375 3.6289827575699195e-105 449 861 2 keep` (row 50) | 7.78 s | 11.63 s | 1,294,212 kB |
| `search hits $RUN --format tsv --out hits.tsv` | (stdout empty; stderr the `contig_label` note) | 6.40 s | 12.30 s | 681,268 kB |
| `search hits $RUN --format faa --out hits.faa` | (same) | 8.07 s | 19.75 s | 1,272,644 kB |
| `nbr $RUN --window 5g --format tsv --out nbr` | (same); 2,581 files | 12.55 s | 19.36 s | 1,279,104 kB |
| `nbr summarize $RUN --window 5g --by nr50` | `nr50  true      0 0.0 0 0 genome_uid` (the empty "unannotated" row); stderr the `member_cluster` note | 10.76 s | 18.40 s | 1,322,636 kB |
| `nbr summarize $RUN --window 5g --by pfam,nr50` (RUNBOOK's levels) | line 20 `pfam  true  ... 2581 1.0 24473 2526 genome_uid ... 1 0.9919`, exit 0 | 10.77 s | 18.58 s | 1,279,972 kB |

The separate `hits.tsv`, `hits.faa` and the 2,581 locus files are byte-identical to those of `search hmm`
(`diff`, `cmp`, `diff -r`: no output), and the `nr50` rows of `nbr summarize --by nr50` equal those of
`neighbors.tsv`. All exports together are under one CPU-minute (*derived*: sum of the CPU column, 100.0 s).

## 3. Against the hand-made answer

Checked by `/tmp/if2453/trial2/check.py` (in the image, 8.23 s, 1,526,576 kB): `hits.tsv` against the
first trial's `hmmsearch_all.tbl`, and against the catalog (`catalog.open_catalog` with a local copy and temp
directory under `/tmp/if2453/trial2`).

| Check | `search hmm` | Hand-made (2026-11-10) | Equal |
|---|---|---|---|
| Proteins at 30 bits or more | 1,957 | 1,957 | yes: 0 only in `hits.tsv`, 0 only in `hmmsearch` |
| Full-sequence score | | | max difference 0.000 over 1,957 (both to 0.1 bit); lowest 30.3 (`8d198546...`), then 30.5, 31.3 |
| Gene rows | 2,581 | 2,581 (`protein_occurrence` of the hit proteins) | yes, same (protein, contig, start, end, strand) set |
| Phages (`S1`): proteins, genes, contigs | 1,931, 2,555, 2,504 | 1,931, 2,555, 2,504 | yes |
| Prokaryotes (`P1`, `P2`) | 26, 26, 22 | 26, 26, 22 | yes |
| T7 gp1 `2eef240a...` | in `hits.tsv` at NC_001604 3170-5822 +1, 367.7 bits, E 4.3e-106; `recovered: true` | 367.7 bits, E 4.3e-106 | yes |
| FASTA | 1,957 records, the hit set; 0 sequences differ from `protein_seq` | | yes |
| Hits in the run's `hits.parquet` | 1,984 rows: 1,930 proteins with 1 row, 27 with 2; all `keep` | | |

The engine is pyhmmer in the runner and the HMMER 3.4 binary by hand, both at Z = 1e8; the scores agree to
the printed precision.

### Neighbours (card W-17)

`rN.1` has no `member_cluster`. `nbr summarize` took the slow path and said so on stderr: `note: release
rN.1 has no member_cluster files (built before card B-27); reading cluster_member instead (the slow path:
every row group is scanned)`. In the first trial the same command reported every neighbour `unannotated`
with exit 0. Now `n_keys` is `{"nr50":3405,"nr90":6482}` and the clusters are named by uid:

| Rank | NR50 cluster | Rep | Rep aa | Loci (`neighbors.tsv`) | Loci (first trial, q5.sql) | Neighbours (both) | Modal offset |
|---|---|---|---|---|---|---|---|
| 1 | `fa4cb451607a4dfdd797df9a5e462694` | `a074c68c...` | 377 | 663 | 663 | 667 | 4 |
| 2 | `b03b93e4638525d3babbccfea1d6cd7b` | `4e512cc9...` | 76 | 493 | 493 | 493 | 2 |
| 3 | `5fe1410ce04d10c0e759661268fcbd6a` | `0b150c83...` | 379 | 481 | 481 | 485 | -1 |
| 4 | `11696bd113cafd2c989b6902980b05bb` | `03b45b3e...` | 116 | 310 | 310 | 310 | 5 |
| 5 | `eb9c25c3c916325df8c65576c9194736` | `444ce148...` | 90 | 299 | 299 | 299 | 3 |
| 6 | `1b667be750e69bc8b4a80e7ce846427f` | `f82eab81...` | 347 | 285 | 285 | 285 | 5 |
| 7 | `eb4f50cb168eba454cd3a4bca985bb7d` | `fcbc3584...` | 94 | 263 | 263 | 263 | 2 |
| 8 | `d6b99b4b297c63773deb00f53ab3b790` | `717b5f15...` | 88 | 253 | 253 | 253 | -2 |
| 9 | `dd4299690e8d61227f321e586ec5a145` | `00757f77...` | 187 | 251 | 251 | 254 | 1 |
| 10 | `1ca31585f9533ee2ef90dde3a9760a26` | `4ca65325...` | 216 | 233 | 233 | 298 | -4 |

Same clusters, representatives, lengths, loci and neighbour counts, in the same order. Header: 1,957 anchor
proteins, 2,581 loci on 2,526 contigs and genomes, 463 censored loci, 24,473 neighbours. The first trial's
SQL counted 24,472 because it counted each neighbour protein once per locus; in one locus
(`locus_ON864052_825a5134...`) two genes carry the same protein (`6bc41626...`) at two offsets. Counted per
gene, the independent window also gives 24,473, and the neighbour multiset of each of the 2,581 locus tables
equals that window's (`check.py`: `2581 of 2581`).

In T7's own locus (`loci/locus_NC_001604_39f9ef26....tsv`) the neighbours at -1 (2020-3100, 359 aa, cluster
`5fe1410c`, rank 3) and +3 (6474-7554, 359 aa, cluster `fa4cb451`, rank 1) sit where T7's protein kinase
(gene 0.7) and DNA ligase (gene 1.3) are. So the two commonest neighbours are probably those families. No
file says so (rough edge 14).

## 4. The outputs as a scientist reads them

- `hits.tsv` (22 columns): `accession header contig_len start end strand length seq_bitscore seq_evalue_zref
  pfam_architecture protein_uid nr90 nr50 nr50_rep genome_uid tier contig_uid feature_uid ordinal tape_version
  dist_to_contig_end_bp n_records`. Coordinates are 0-based, end exclusive (T7 gp1 3170-5822 is GenBank
  3171..5822), but nothing in the file says so. `tier` is empty in all 2,581 rows (`contig_dim.tier` is NULL
  for all 32,815 contigs of `rN.1`). `pfam_architecture` is `-` throughout (no library; the log says so).
  `n_records` (1 to 14) is the number of source records of the contig (`export.py`), which is easy to read as
  hit records. `dist_to_contig_end_bp` is 0 for 17 hit genes that end 1 or 2 bp from a contig end (for
  example AP013538 42996-44703 of 44,705). In the release it is 0 for 32,312 genes, every one of them within
  3 bp of an end. It seems to mean "runs off the edge" (protein `869087d0...` there begins with N, not M), not
  a distance.
- `hits.faa`: one line per sequence. The header is `>2eef240a... length=883 seq_bitscore=367.7
  seq_evalue_zref=4.3e-106 n_genes=33 example=AY264774:3171-5823(+1) AY264774 Enterobacteria phage T7,
  complete genome.`: the example is the first of 33 genomes, not the positive control's NC_001604, and
  `3171-5823` is 0-based but looks like a 1-based GenBank range.
- `search hits` (table): 50 rows of 1,984 with no line saying it stopped at 50; the protein is `target_uid`;
  `seq_bitscore` and `bitscore` (the best domain's) side by side at 13 digits (`397.3708190917969`);
  `evalue_zref` says nothing of the Z it refers to (1e8).
- Locus tables: a comment line with the header, a key line (`contig_uid`, `anchor_feature_uid`,
  `anchor_ordinal`, `anchor_protein_uid`, `tape_version`, `policy_hash`, `censored`), then `offset offset_rel
  ordinal feature_uid start end strand protein_uid length_aa pfam_architecture nr90 nr50 hit_seq_bitscore
  delta_bp crosses_origin`. `offset` follows coordinates and `offset_rel` follows the anchor's strand (anchor
  on -1: offset -5 is `offset_rel` 5). `neighbors.tsv`'s `offsets` and `modal_offset` are the strand-relative
  ones (cluster `5fe1410c` is at -1 in 476 of 481 loci), but their name does not say so. `delta_bp` is
  unsigned.
- `neighbors.tsv`: `key` and `cluster_uid` are the same value; `model_name` is empty; a cluster is a 32-hex
  uid with a representative uid and length, and no name, accession or product. Strand is `+1` in `find`, `1`
  in `hits.tsv` and the locus tables, `(+1)` in `hits.faa`.
- Ordinals: every ordinal has its `feature_uid` beside it. `hits.tsv` has `ordinal` and `feature_uid` (0 rows
  with one and not the other); every locus table has `ordinal` and `feature_uid` in each row (0 rows without)
  and `anchor_ordinal` with `anchor_feature_uid` in its key line; the run's `neighborhoods.parquet` has
  `anchor_ordinal` with `anchor_feature_uid` and `ordinal` with `feature_uid`; `neighbors.tsv`, `hits.faa` and
  the `search hits` table hold no ordinal (AGENTS 3.3).

## 5. Clean-up

- Root: a listing of its top two levels (`find -mindepth 1 -maxdepth 2 -printf '%y %T %p'`, 56 entries) at
  01:37:32 before anything ran and at 01:56:02 after the last run are identical (`diff` prints nothing), and
  `find /tmp/if2453/r10/root -newer root_before.txt -maxdepth 4` printed nothing. So in this window nothing
  at depth 4 or less changed, neither ours nor the scan's: the scan's `l1/` and `runs/stage=P2` entries kept
  their 01:22 and earlier times. Had they changed, the scan's entries would have been told apart by path
  (`l1/`, `runs/stage=P2`, `.incoming/_build-stage.P2.*`, written by the `r10-p2-*` units) against ours (none
  expected: the root was bound read-only into every command). No STOP or HALT file was found (`ls -d`, 01:38
  and 01:50), and none was made.
- Pack: sizes and modification times of every file under `/tmp/if2453/w16/packs` are unchanged. The
  directory `release=rN.1` changed its time from 00:54:33 to 01:41:47: `pack.verify_dir`, run by every search
  submit, makes a scratch sorter directory `.verify-<pid>-<hex>` beside the pack (`_count_source`,
  `src/biobox/search/pack.py`) and removes it. None is left.
- Left on the pod, 30 MB in all (`du -sh /tmp/if2453/trial2`): `/tmp/if2453/trial2/out/rnap` (the
  `search hmm` output, 8,599,617 B), `/tmp/if2453/trial2/out/sep` (the separate exports),
  `/tmp/if2453/trial2/tmp/search/if2453/runs/run=87021d87e0f41e453c14d35763d51d17` (the run, 1,341,398 B with
  its pending record `runs/pending/run=87021d87....json`, which stays after the close), `logs/` (unit log,
  samples, every command's stdout, stderr and `time -v`), `work/` (copies of the first trial's HMM, table and
  SQL, the root and pack listings), `cat/` (the catalog copy), the helper scripts. The work directory of the
  run was removed by the run. The tree `/tmp/if2453/ci/trial2` is 16 MB. No `trial2-*` unit is loaded.
  Nothing under `/stor` was written.

## 6. Commands and last lines

All on wilkcomp02 through `ssh -o BatchMode=yes`; `$B` is `bash /tmp/if2453/trial2/bb.sh <name>` and `$W` is
`--root /tmp/if2453/r10/root --profile /tmp/if2453/trial2/profile.yaml --pack-dir /tmp/if2453/w16/packs`.

| Time | Command | Last line |
|---|---|---|
| 01:37 | (laptop) `sh scripts/pod_pytest.sh trial2 tests/r0/test_runbook.py` | `32 passed in 4.15s` |
| 01:38 | `bash load.sh` | `... load 44.03 43.11 43.81 ours 4494.6% cap 5600% launch yes` |
| 01:39 | `$B find1 find --root ... --release rN.1 --accession NC_001604` | `exit 0` (49 genes) |
| 01:39 | `$B find2 find ... --accession NC_001604:3170-5822` | `NC_001604 719747a8... 3170 5822 +1 39f9ef26... 2eef240a...` |
| 01:39 | `$B find3 find ... --protein 2eef240afefddcd2c8ac7660eb50e384` | `PP384410 ... 3171 5823 +1 ... 2eef240a...` (33 genes) |
| 01:39 | `bash .../pod_run.sh ... trial2-hmm search hmm ...` (section 2) | `pod_run: log /tmp/if2453/trial2/logs/trial2-hmm.log; watch: systemctl --user status trial2-hmm` |
| 01:39 | `systemd-run --user --unit=trial2-sample --collect --quiet nice -n 19 bash sample.sh trial2-hmm logs/trial2-hmm.hwm` | sampler's last line `01:48:53 unit trial2-hmm state 'inactive': stop` |
| 01:43 | `$B status1 search status 87021d87... $W` | `{"run_id": "87021d87e0f41e453c14d35763d51d17", "status": "running"}` |
| 01:48 | (unit log) | `pod_run: trial2-hmm exited 0` |
| 01:51 | `$B status2 search status 87021d87... $W` | `{... "status": "ok", "error": null, "rounds": 1, "final_round": 1, "stop_reason": "max_rounds"}` |
| 01:51 | `$B hits_table search hits $R $W`, `--format tsv`, `--format faa`, `nbr $R $W --window 5g --format tsv`, `nbr summarize $R $W --window 5g --by nr50`, `--by pfam,nr50` | see section 2, separate exports |
| 01:54 | `bash py.sh check check.py out/rnap tmp/search/if2453/runs/run=87021d87...` | `hits.tsv rows with ordinal but no feature_uid: 0` |
| 01:55 | `bash py.sh q1 q.py` (tier, `dist_to_contig_end_bp`) | `('869087d0abf96a74643c911105e9c224', 'C', 'N', 568)` |
| 01:56 | `systemd-run --user --unit=trial2-echo --collect -p StandardOutput=append:... /bin/sh -c 'use="A=1 B=2 "; echo "[${use% }] [$use]"'` | `[] [A=1 B=2 ]` |
| 01:56 | root and pack listings, `diff` | root: no output; pack: only the `release=rN.1` directory time |
| 01:58 | `bash py.sh nrec nrec.py` (what `n_records` counts) | `T7 gp1 parquet: [(1, 449, 861, 359.88592529296875, 367.7415466308594, 'exact')]` |

## What a user would trip over

1. The search worker loads the whole target into Python (`engine.load_target`: a dict of `str`, a `TextSequenceBlock` and its digitized copy): 26.27 GiB resident for 10.6 M proteins (8.4 B per residue, derived), filled at one CPU for about 3.4 of the run's 5.7 minutes. A 16G cap fails, and a release several times larger would not fit 200G. Search shard by shard or in bounded blocks (the W-19 suspicion, confirmed).
2. `runner.resolve` decompresses every `unique` shard into a set of all 10.6 M uids only to check that the positive control is in the target: the runner reached 5.65 GiB before the worker started. Look up the few positive uids while streaming, or from an index.
3. Every search submit runs the full `pack.verify_dir` (sha256 of 5.5 GB, record counts, and the release's distinct proteins through a sorter), about 110 s at one CPU before the run exists. Verify once per pack and record it, or check only the shards the target needs.
4. That verify writes a scratch `.verify-<pid>-<hex>` directory into the shared pack's `release=<R>/` directory. A user who can read but not write the shared pack cannot search, and a SIGKILL leaves the directory behind. Put the scratch in the user's work directory.
5. `search hmm` has no run id, pending record or `search status` for its first 160 s, and its log says nothing between the first slow-path note and `search run: ... running`. Write the pending record first, and print a line per phase (lookup, verify, stage, load, search, export).
6. `pod_run.sh`'s totals line (card W-18) is always empty, `pod_run: trial2-hmm usage `: systemd expands `${use% }` in the unit's command line as its own variable. A test unit running `echo "[${use% }] [$use]"` logged `[] [A=1 B=2 ]`. Escape it as `$${use% }` or avoid `${...}`, and add a test that reads a real unit's line.
7. `nbr summarize --by pfam,nr50`, the command RUNBOOK section 7 prints, exits 0 on a release with no domain library and reports every neighbour as `pfam` unannotated (2,581 loci, share 1.0) with no note. `search hmm` switches to `nr50,nr90` by itself; the separate command should refuse `pfam` there, or say why.
8. `tier` is empty in every row of `hits.tsv` (the release's `contig_dim.tier` is NULL everywhere), with no note.
9. `dist_to_contig_end_bp` is 0 for genes 1 or 2 bp from a contig end (17 hit rows; 32,312 genes of the release). It apparently flags genes that run off the edge, which the name does not say. Document it, or add an edge or partial flag and keep the distance.
10. `ACC:START-END` (in `--positive`, `find` and the `hits.faa` `example=`) is 0-based and end-exclusive, but biologists read it as a GenBank range: T7 gp1, GenBank 3171..5822, is `3170-5822`, and on AY264774 `3171-5823`. No output file says which convention it uses. Put it in the TSV headers and the FASTA, or print 1-based ranges in the `ACC:START..END` style.
11. The `hits.faa` `example=` is the first of a protein's genomes in sort order (AY264774 of 33 for T7 gp1), not the positive control's NC_001604, so the control cannot be found by its accession in the FASTA.
12. Strand is written three ways: `+1` (`find`), `1`/`-1` (`hits.tsv`, locus tables) and `(+1)` (`hits.faa`).
13. Offsets: the locus tables have `offset` (coordinate order) and `offset_rel` (anchor's strand), while `neighbors.tsv` has `offsets` and `modal_offset`, which are strand-relative, without `_rel` in their names. A user reading a minus-strand locus gets the sign backwards.
14. Clusters have no names. `neighbors.tsv` gives the uid twice (`key`, `cluster_uid`), `rep_uid`, `rep_length` and an empty `model_name`. Without Pfam, a user cannot tell that the top two clusters are probably DNA ligase and protein kinase except by `find --protein <rep>` and an outside lookup. Add an example `ACC:range` and the source annotation (product) of each representative. `hits.tsv` calls the representative `nr50_rep`, `neighbors.tsv` calls it `rep_uid`.
15. The `search hits` table prints 50 of 1,984 rows with no note that it stopped, calls the protein `target_uid` (`protein_uid` elsewhere), prints scores to 13 digits, and shows `seq_bitscore` and `bitscore` (best domain) without saying which is which.
16. `n_records` in `hits.tsv` (1 to 14) is the number of source records of the contig. It is not documented in the file, and reads as hit records.
17. The pending record `runs/pending/run=<id>.json` stays after the run closes. RUNBOOK section 9 lists it, but does not say whether a closed run's record may be removed.
18. Each export repeats the slow-path note (9 in one `search hmm` log), and the exports print nothing on success (no row or file count). One note per command and a one-line summary per file would read better.
