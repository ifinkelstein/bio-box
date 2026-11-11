# Trial: a custom-HMM search with neighbourhoods on `rN.1` (2026-11-10)

Question: can a scientist, today, run a real custom-HMM search with gene neighbourhoods against the
10 Gbp release `rN.1` on the pod, with the current code? This is a trial, not a task card: no code under
`src/` or `tests/` changed. Every number is *measured* on wilkcomp02, 2026-11-10, between 19:28 and
22:16 CST, by the command named beside it, unless it is labelled *derived*.

Setup:

- Code: `origin/plan/task-cards` at `f7c511e` (after B-27, W-03, W-13), mirrored to `/tmp/if2453/ci/trial`
  with the rsync of `scripts/pod_pytest.sh` (`sh scripts/pod_pytest.sh trial tests/r0/test_runbook.py`:
  `31 passed in 4.16s`).
- Data: root `/tmp/if2453/r10/root`, release `rN.1`, stages `S1` (phages, inphared: 27,258 contigs,
  1,664,109,938 bp), `P1` (539 contigs, 1,016,625,803 bp) and `P2` (5,018 contigs, 9,596,180,328 bp);
  12,276,916,069 bp in all (*derived*, sum of the three). Built with the tree at `740d049`, before B-27.
  No domain library in the release (`l1_domains` has no file in its manifest).
- Image `bio-box-test_74b3e05.sif`, sha256 `c30be585...c898`. Every command ran in a
  `systemd-run --user` scope or a `scripts/pod_run.sh` unit at nice 19, memory cap 32G (16G for the short
  ones), CPU quota at most 8 CPUs, with the source bound in by `PYTHONPATH` (RUNBOOK section 0).
- The root was treated as read-only and, to make sure, bound read-only into the image
  (`--bind /tmp/if2453/r10:/tmp/if2453/r10:ro` after `--bind /tmp/if2453`). Every output went to
  `/tmp/if2453/trial`. Nothing under `/stor` was written; Pfam was read from
  `/stor/scratch/.../pfam/38.2/Pfam-A.hmm.gz` through a read-only bind.
- Profile `/tmp/if2453/trial/profile.yaml`: `profiles/wilkcomp02.yaml` with `tmp: /tmp/if2453/trial/tmp`.
- Query: Pfam `RNA_pol` (PF00940.25, 449 match states), extracted with `hmmfetch` in the image.
  Cutoff 30 bits sequence score at Z = 1e8. Positive control: T7 RNA polymerase gp1 (NP_041960.1).

## Result in one paragraph

The current code reads `rN.1`: `release verify` passes, the catalog opens, `contig_label`, the
postings and the sequence macros answer, and `export proteins` writes all 10,620,684 proteins. The
search path stops at its first step: `pack build` does not fit the 32G memory cap. It filled the cap
2 h 3 min in, while collecting dark ORFs, and then thrashed for 36 min without writing a byte,
and the trial stopped it after 2 h 39 min and 155 CPU-minutes. So `search run`, `search status`,
`search hits` and `nbr <run>` could not be tried. `nbr summarize` given the hit proteins as uids runs,
but on this release it is silently wrong: it reports every one of 24,473 neighbours as `unannotated` at
`nr50` (exit 0, no warning), because the release has no `member_cluster` table (B-27 came after it was
built). By hand (hmmsearch from the image on the exported FASTA, and SQL on the catalog) the answer is
1,957 proteins at 30 bits or more: 1,931 in phages, 26 in prokaryotes; T7 gp1 scores 367.7 bits.

## 1. Reading `rN.1` with the current code

| Command | Last line | Time |
|---|---|---|
| `bio-box release verify --root /tmp/if2453/r10/root --release rN.1` | `release rN.1: verified (stages S1,P1,P2, 52 files, 34 artifacts, manifest hash 84c2c057b09643bb4e9925938703095d)` | 60.4 s wall, 41.8 s user, 5.8 s system, 327 MB peak RSS (`/usr/bin/time -v`) |
| `python trial.py sql <catalog> q1.sql` (`catalog.open_catalog`, a local copy of the catalog in `/tmp/if2453/trial/cat`) | `open_catalog ok in 0.44 s` | |

`trial.py` (in `/tmp/if2453/trial/work`) opens the root's `catalog/release=rN.1.duckdb` with
`catalog.open_catalog(path, local_dir=/tmp/if2453/trial/cat, temp_directory=/tmp/if2453/trial/tmp/duck)`
and runs SQL statements from a file. Without `temp_directory`, `open_catalog` puts DuckDB's spill
directory beside the catalog, under the data root (rough edge 3).

What answered (q1.sql to q5.sql, each in under 8 s):

- `catalog_info`: release `rN.1`, root `/tmp/if2453/r10/root`, tape version `t-ca596dd85040`, DuckDB 2.0.0.
- `contig_label` (in this catalog a view over `contig_occurrence`, the pre-B-27 form): 32,815 rows;
  `WHERE accession LIKE 'NC_001604%'` gives contig `719747a8fef97b2380955d92845861db`,
  `NC_001604 Enterobacteria phage T7, complete genome.` (0.08 s).
- `protein_occurrence` on that contig: 49 genes; the 883 aa gene at 3170-5822 (+), N-terminus
  `MNTINIAKNDFSDIELAAIPFNTLADHYGE`, is T7 gp1 (GenBank gene 1 is 3171..5822, 1-based): protein uid
  `2eef240afefddcd2c8ac7660eb50e384`. `contig_len_of` gives 39,937 and `seq_get(..., 3170, 3200, 1)`
  gives `ATGAACACGATTAACATCGCTAAGAACGAC`.
- `cluster_member` by member uid (0.18 s): T7 gp1 is in NR50 cluster `5ea5af7357bbdd6d2ad9e3ddfe56c605`
  and NR90 cluster `c1bf639f92687b483b659b4dadf7c52c`.
- Postings: `token_of('nr50', <cluster>)` = 4611686018428485233 (1.61 s); `token_rows` of it: 1,064
  postings on 1,008 contigs (0.24 s); `partners_of(<token>, 5)`: top partner on 664 loci (7.73 s).
- `token_dict` holds `nr50` (2,966,200) and `nr90` (8,278,537) tokens only: no `family` level.
- `member_cluster`: `Catalog Error: Table with name member_cluster does not exist!` The release has no
  `member_cluster` or `contig_label` files (B-27 lookup tables); its manifest lists 20 tables.

Would a release built by the current code have what is missing? Yes: ingest writes `contig_label`
(`src/biobox/l0seq/store.py`, hook `contig_label`) and clustering writes `member_cluster`
(`src/biobox/l2/cluster.py`, `table("member_cluster")`). Rebuilding `rN.1` was out of scope (it writes
the root).

## 2. The search path of RUNBOOK section 7

### Pack build

```
bash /tmp/if2453/ci/trial/scripts/pod_run.sh --profile /tmp/if2453/trial/profile.yaml --image $IMAGE \
  --digest sha256:c30be585...c898 --mem 32G --threads 8 --bind /tmp/if2453 \
  --bind /tmp/if2453/r10:/tmp/if2453/r10:ro --env PYTHONPATH=/tmp/if2453/ci/trial/src:/tmp/if2453/pilot_s1/pydeps \
  trial-pack pack build --root /tmp/if2453/r10/root --release rN.1 --profile /tmp/if2453/trial/profile.yaml
```

Started 19:35:53. No `--model`: the release has no domain library, and `pack build` accepted that
(the `--model` of RUNBOOK section 7 is needed only for a library the release has; an `hmm_file` seed
does not need the pack's models). The pack directory went outside the root, to
`/tmp/if2453/trial/tmp/packs/release=rN.1/.building-2540568-a86914a7/` (W-03 works).

Progress, from the build directory's file times and the unit's accounting (`systemctl --user show
trial-pack -p CPUUsageNSec -p MemoryCurrent -p MemoryPeak`, sampled every 10 s into
`/tmp/if2453/trial/work/pack_acct.txt`):

| Phase | Ended | Files | Unit memory at the end |
|---|---|---|---|
| `unique` | 19:39 | `part-0.faa.zst` 1,242,095,090 B, `part-1.faa.zst` 905,151,551 B | |
| `igr` | 20:30 | `part-0.fna.zst` 606,018,363 B, `igr_meta.parquet` 530,288,062 B | 22.1 GB peak (anonymous memory fell to 16.2 GB after it) |
| `dark_orf` collection | 21:31 (directory made) | | 33,351,327,744 B in the sample at 21:28 |
| `dark_orf` shards | 21:33:27 | `part-0.faa.zst` 1,313,881,169 B, `part-1.faa.zst` 66,819,238 B | 34,359,476,224 B at 21:38:44, the 32 GiB cap (34,359,738,368 B) |
| `dark_orf/dark_orf_meta.parquet` | never | | anonymous memory 34,204,147,712 B from 21:48 on, page cache 0.3 MB |

After the last shard the process wrote nothing for 41 minutes (`write_bytes` in `/proc/<pid>/io` stayed
at 4,664,262,656). It ran at about one CPU, thrashing: the cgroup's `memory.events` counted 934,924
`max` events, `pgmajfault` rose from 524,267 (22:04) to 790,903 (22:14), and it read about 1.2 GB from
disk every 30 s (`read_bytes`) with no write; `oom_kill` stayed 0, so the kernel never ended it. The
step it was in when stopped is, from the order in `_build_dark_orf`, `_write_meta`, which builds one
Arrow table from the whole list of per-ORF dicts. The trial stopped it at 22:14:43 with
`systemctl --user stop trial-pack`; the log's last lines are `pod_run: trial-pack received SIGTERM` and
`pod_run: trial-pack exited 143`.

| Pack build | Value |
|---|---|
| Wall time | 9,528 s (2 h 38 min 48 s; `ActiveEnterTimestamp` 19:35:55 to the stop at 22:14:43) |
| CPU | 9,305.4 CPU-seconds, 155.1 CPU-minutes (`CPUUsageNSec` 9,305,374,709,000 just before the stop) |
| Peak memory | the cap, 34,359,738,368 B (`MemoryPeak`) |
| On disk when stopped | 4,664,253,473 B in `.building-2540568-a86914a7/` (`du -sb`), no `pack_manifest.json` |
| Left behind | the `.building-...` directory: SIGTERM ends Python without running the `except BaseException` cleanup. Removed by exact path afterwards |

The pack is not finished, so its full size is unknown. The three finished target files are already
4.66 GB (*derived*: unique 2,147,246,641 + igr 606,018,363 + 530,288,062 meta + dark ORF shards
1,380,700,407 = 4,664,253,473 B).

Why: `_build_igr` and `_build_dark_orf` (`src/biobox/search/pack.py`) collect every region's metadata
row and sequence in Python lists for the whole release and sort them before writing the first shard.
That is the unbounded-memory pattern AGENTS section 4.1 forbids. The phage-only store of 2026-11-09
(1.66 Gbp) built in 495 s and 632 MB on disk; this release is 7.4 times larger (*derived*,
12,276.9 / 1,664.1 Mbp).

### What could not be tried

Every later step needs the pack. After the stop,

```
bio-box search run /tmp/if2453/trial/work/recipe.yaml --root /tmp/if2453/r10/root --release rN.1 \
  --profile /tmp/if2453/trial/profile.yaml --cpus 4 --mem-gb 16 --wait
```

(in a 16G, 4-CPU scope) exits 2 with `refused: no search pack of release rN.1 in the pack directory
/tmp/if2453/trial/tmp/packs; build it with: bio-box pack build --root /tmp/if2453/r10/root --release rN.1
--pack-dir /tmp/if2453/trial/tmp/packs --model <library>-<version>=<HMM file> ...`. So none of these ran: `search run --wait`, `search status`,
`search hits` (table, `--format tsv`, `--format faa`) and `nbr <run> --window 5g`. The recipe that
would have run is `/tmp/if2453/trial/work/recipe.yaml`:

```yaml
recipe: rna_pol_pf00940_30bits
written_against: rN.1
seed: {hmm_file: /tmp/if2453/trial/work/RNA_pol.hmm, label: Pfam RNA_pol PF00940.25}
target: unique
tool: hmmsearch
iterations: {max: 1, inclusion_bits: 30}
context: {window: 5g}
controls:
  positive:
    - {protein_uid: 2eef240afefddcd2c8ac7660eb50e384, label: "T7 RNA polymerase gp1 (NP_041960.1)"}
```

The positive control had to be given as a protein uid found by hand (section 1): a
`protein_accession` control is valid on paper but not executable (`RUN_POSITIVE_KINDS` in
`src/biobox/search/recipe.py`), and there is no command that maps a contig accession and coordinates
to a protein uid.

### `nbr summarize` without a run

`nbr summarize` also takes protein uids and `--release`, and needs no pack. Given the 1,957 hit uids of
section 3:

```
bio-box nbr summarize <1,957 uids, comma-separated> --root /tmp/if2453/r10/root --release rN.1 \
  --profile /tmp/if2453/trial/profile.yaml --window 5g --by nr50 --top 10
```

Exit 0, 14.7 s wall, 10.5 s user, 1.23 GB peak RSS. Header: 1,957 anchor proteins, 2,581 anchor loci on
2,526 contigs and genomes, 463 censored loci, 24,473 neighbours, `n_keys {"nr50": 0}`. The only row:
`nr50 unannotated 2581 1.0 24473` with `same_strand_share 0.9919`. Every neighbour is "unannotated"
because `nbr.summarize` reads cluster keys from `member_cluster`, and a missing table becomes an
empty view (`_store_view`), not an error. A SQL check (q5.sql) finds 0 of 24,472 neighbours without an
NR50 cluster in `cluster_member`. The answer is wrong and looks plausible: rough edge 1.

## 3. The answer by hand (independent of the runner)

| Step | Command | Result |
|---|---|---|
| Query HMM | `hmmfetch -o RNA_pol.hmm /refdb/Pfam-A.hmm.gz RNA_pol` (image, Pfam bound read-only) | `Retrieved HMM RNA_pol.`; sha256 `32e6f0f0b1ef15d23d8ead23095f0b2a5a3323551b899bba99218d75f830e87e` |
| Export | `bio-box export proteins --root ... --release rN.1 --threads 2 --out /tmp/if2453/trial/work/proteins.faa` | exit 0; 10,620,684 records, 3,712,392,579 B; 88.2 s wall, 69.6 s user, 35.9 s system, 2.28 GB peak RSS |
| Search | `hmmsearch --cpu 6 -Z 1e8 --domZ 1e8 -T 30 --noali --tblout hmmsearch_all.tbl RNA_pol.hmm proteins.faa` (HMMER 3.4, image) | 1,957 proteins; 26.3 s wall, 161.1 s user (2.7 CPU-minutes, *derived*), 88 MB peak RSS |
| Control | `grep 2eef240a... hmmsearch_all.tbl` | T7 gp1: E = 4.3e-106, 367.7 bits (best domain 359.9): a hit well above 30 bits |
| Split | q4.sql: hit uids joined to `protein_occurrence` (its `stage`) and `contig_dim` | below |

| | Proteins | Genes | Contigs (= genomes) |
|---|---|---|---|
| Phages (`S1`) | 1,931 | 2,555 | 2,504 |
| Prokaryotes (`P1` 8, `P2` 18) | 26 | 26 | 22 |
| All | 1,957 | 2,581 | 2,526 |

No protein is in both a phage and a prokaryote. The hits fall in 121 NR50 and 432 NR90 clusters. The
phage count equals the 1,931 Pfam hits of the 2026-11-09 walk-through on the phage store `rS1.pilot1`
(same model, same cutoff, same Z), so the phage part of `rN.1` reproduces it. The prokaryote hits score
254 to 392 bits, in chromosomes of *Neorhizobium*, *Brevundimonas* (two), *Xanthomonas*,
*Pseudoblastomonas*, *Variovorax*, *Pseudoalteromonas*, *Pseudomonas*, *Pyruvatibacter*, *Pandoraea*,
*Vibrio* and several MAGs (top 15 by score in q4.log): full-length polymerases, likely prophages, not
checked further.

### Top ten NR50 clusters within 5 genes (q5.sql)

An independent window, not `nbr`: for each of the 2,581 hit genes, the genes of its contig ordered by
(start, end, strand), offsets -5..-1 and 1..5, wrapping on circular contigs (`contig_dim.is_circular`),
cut at the ends of linear ones; each neighbour protein counted once per locus; clusters from
`cluster_member` at `nr50`. 24,472 neighbour rows (`nbr summarize` counted 24,473).

| Rank | NR50 cluster | Representative | Rep length (aa) | Loci | Share of 2,581 loci |
|---|---|---|---|---|---|
| 1 | `fa4cb451607a4dfdd797df9a5e462694` | `a074c68ca618ce89e6055c94765a53fb` | 377 | 663 | 0.2569 |
| 2 | `b03b93e4638525d3babbccfea1d6cd7b` | `4e512cc9a9fbde6aad298c9fdec6a1b2` | 76 | 493 | 0.1910 |
| 3 | `5fe1410ce04d10c0e759661268fcbd6a` | `0b150c8369b70d24f45a5f37825ed33b` | 379 | 481 | 0.1864 |
| 4 | `11696bd113cafd2c989b6902980b05bb` | `03b45b3ee8dbd84c3e64ed751c691c31` | 116 | 310 | 0.1201 |
| 5 | `eb9c25c3c916325df8c65576c9194736` | `444ce148420dfb916d5398561e626122` | 90 | 299 | 0.1158 |
| 6 | `1b667be750e69bc8b4a80e7ce846427f` | `f82eab8195089f4a1b30bdab5e07c50c` | 347 | 285 | 0.1104 |
| 7 | `eb4f50cb168eba454cd3a4bca985bb7d` | `fcbc3584eb0cd3744deca413645aa75b` | 94 | 263 | 0.1019 |
| 8 | `d6b99b4b297c63773deb00f53ab3b790` | `717b5f15f7fa62ac518ceb6700ad77f3` | 88 | 253 | 0.0980 |
| 9 | `dd4299690e8d61227f321e586ec5a145` | `00757f77b9ff3b5b400512a14c9614b3` | 187 | 251 | 0.0972 |
| 10 | `1ca31585f9533ee2ef90dde3a9760a26` | `4ca6532568f642c0f5404b7b8cb2e1d1` | 216 | 233 | 0.0903 |

Shares are *derived* (loci / 2,581). Without a domain library the clusters have no names; the
2026-11-09 walk-through named the commonest neighbours of the same family (`Phage_T7_1_1`, DNA ligase,
`Head-tail_con`, T7 kinase) with Pfam domains, which `rN.1` does not have yet. `partners_of` on T7 gp1's
own NR50 token gives a similar head (664, 501, 484 loci), from the postings path.

## 4. Commands and last lines

All run on wilkcomp02 through `ssh -o BatchMode=yes`, in `/tmp/if2453/trial/work`, with `$SCOPE` =
`systemd-run --user --scope --quiet -p MemoryMax=<16G or 32G> -p MemorySwapMax=0 -p CPUQuota=<N>00% nice -n 19`
and `$SING` = `singularity exec --cleanenv --env PYTHONPATH=/tmp/if2453/ci/trial/src:/tmp/if2453/pilot_s1/pydeps
--env OMP_NUM_THREADS=2 --bind /tmp/if2453 --bind /tmp/if2453/r10:/tmp/if2453/r10:ro --bind /sys/fs/cgroup:/sys/fs/cgroup:ro $IMAGE`.

| Time | Command | Last line |
|---|---|---|
| 19:28 | `date; uptime; ...; top` (load and our CPU) | load 48.22; our processes 24.2 CPUs |
| 19:29 | `sh scripts/pod_pytest.sh trial tests/r0/test_runbook.py` (laptop) | `31 passed in 4.16s` |
| 19:31 | `$SCOPE $SING python -m biobox.cli release verify --root /tmp/if2453/r10/root --release rN.1` | `release rN.1: verified (stages S1,P1,P2, 52 files, 34 artifacts, manifest hash 84c2c057b09643bb4e9925938703095d)` |
| 19:32 | `$SCOPE $SING python trial.py sql /tmp/if2453/r10/root/catalog/release=rN.1.duckdb q1.sql` | `Catalog Error: Table with name member_cluster does not exist!` (the last statement; the others answered) |
| 19:33 | same with q2.sql (T7 genes) | `# 2 rows in 0.04 s` |
| 19:34 | same with q3.sql (clusters, postings, sequence) | `ATGAACACGATTAACATCGCTAAGAACGAC` |
| 19:35 | `$SCOPE singularity exec ... hmmfetch -o RNA_pol.hmm /refdb/Pfam-A.hmm.gz RNA_pol` | `Retrieved HMM RNA_pol.` |
| 19:35 | `bash .../pod_run.sh ... trial-pack pack build --root ... --release rN.1 --profile ...` | `pod_run: log /tmp/if2453/trial/tmp/logs/trial-pack.log; watch: systemctl --user status trial-pack` |
| 21:24 | `$SCOPE $SING python -m biobox.cli export proteins --root ... --release rN.1 --threads 2 --out proteins.faa` | exit 0 (no output; 10,620,684 records) |
| 21:26 | `$SCOPE singularity exec ... hmmsearch --cpu 6 -Z 1e8 --domZ 1e8 -T 30 --noali --tblout hmmsearch_all.tbl -o hmmsearch_all.out RNA_pol.hmm proteins.faa` | exit 0; 1,957 rows in the table |
| 21:27 | same `trial.py sql` with q4.sql (phage / prokaryote split) | `# 15 rows in 0.04 s` |
| 21:28 | `$SCOPE $SING python -m biobox.cli nbr summarize $(cat hmmsearch_uids.txt) --root ... --release rN.1 --profile ... --window 5g --by nr50 --top 10` | `nr50  true ... 2581 1.0 24473 2526 genome_uid ... 1 0.9919` |
| 21:30 | same `trial.py sql` with q5.sql (independent window) | `neighbours_without_nr50 0` |
| 21:38 | (unit sample) | `MemoryCurrent` 34,359,476,224 B: at the cap |
| 22:14 | `systemctl --user stop trial-pack` | log: `pod_run: trial-pack exited 143` |
| 22:14 | `rm -r /tmp/if2453/trial/tmp/packs/release=rN.1/.building-2540568-a86914a7` | exit 0 |
| 22:15 | `$SCOPE $SING python -m biobox.cli search run recipe.yaml --root ... --release rN.1 --profile ... --cpus 4 --mem-gb 16 --wait` | `refused: no search pack of release rN.1 in the pack directory /tmp/if2453/trial/tmp/packs; build it with: ...` (exit 2) |

The SQL files, `trial.py`, the logs and the hit table are in `/tmp/if2453/trial/work` on the pod.
Left there: `proteins.faa` (3,712,392,579 B, the export), `hmmsearch_all.tbl` and `.out`, the uid
lists, `pack_acct.txt` (the unit samples), `RNA_pol.hmm`, `recipe.yaml`; the catalog copy in
`/tmp/if2453/trial/cat`; the unit log in `/tmp/if2453/trial/tmp/logs/trial-pack.log`; 3.5 GB in all
(`du -sh /tmp/if2453/trial`). No `trial-*` unit is running. Nothing was written under
`/tmp/if2453/r10` or `/stor`.

## What a user would trip over

1. `nbr summarize --by nr50` on a release without `member_cluster` reports every neighbour as `unannotated` with exit 0: a missing lookup table must be an error that names the table and the card (B-27), not an empty view.
2. `pack build` holds every IGR and dark ORF of the release in Python lists before writing (AGENTS 4.1): on `rN.1` it filled a 32G cap and thrashed in `_write_meta` after 155 CPU-minutes, never finishing; stream regions to shards per batch of contigs, with a memory test at 12 Gbp.
3. `catalog.open_catalog` defaults DuckDB's `temp_directory` to `<catalog path>.tmp`, under the data root: a reader who cannot write the root fails on the first spill, and one who can writes into it; default it to the profile's `tmp`.
4. There is no command that finds a protein by contig accession and coordinates (or by an NCBI protein accession): the positive control had to be found with hand-written SQL over `contig_label` and `protein_occurrence`; add `bio-box find` (accession, header text, coordinates) printing protein uids.
5. A `protein_accession` positive control is valid on paper but refused at submit ("not executable"): either resolve it through `contig_label` and the release's protein table, or say in the recipe docstring which kinds run.
6. `pack build` prints nothing until it ends (no phase, no count, no progress) in a run that takes hours; log each target as it starts and ends, with counts.
7. The unit's CPU and memory totals are lost when a `pod_run.sh` unit ends (`--collect`): write `CPUUsageNSec` and `MemoryPeak` to the log's last line, so a user can report what a command cost without sampling it.
8. Old releases (`rN.1`) are readable but not searchable with all features: `release verify` passes and says nothing about missing lookup tables; have it (or `pack build`) list the registry tables the release lacks that current readers need.
9. `nbr summarize` takes anchors only as a comma-separated argument (1,957 uids made a 64,581-byte argument); accept a file of uids or a hit table.
10. RUNBOOK section 7 says `--model` "must match a domain library of the release" but not that it may be left out when the release has none, nor that an `hmm_file` seed needs no model in the pack; the refusal of `search run` without a pack also prints `--model <library>-<version>=<HMM file>` as if it were required.
12. A thrashing process at its cgroup cap is never OOM-killed (`oom_kill` 0 after 36 minutes at the cap, about one CPU busy): `pod_run.sh` units could set `MemoryHigh` below `MemoryMax` or a watchdog on `memory.pressure`, so a stuck build fails fast instead of burning CPU-minutes.
11. A pack build ended from outside (SIGTERM from `systemctl --user stop`, as here; SIGKILL; OOM) leaves its `.building-<pid>-<hex>` directory behind (its cleanup is an `except BaseException`); RUNBOOK section 9 should list it among the leftovers to remove by exact path.
