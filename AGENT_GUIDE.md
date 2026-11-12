# Agent guide: answering a scientist's question with bio-box

For an LLM agent with a shell that **uses** bio-box releases. Agents that change the code follow `AGENTS.md`
instead. Every command, flag, table, column and macro named here is in the code; the source is cited per section.
State of the data: 2026-11-11. It changes; §2 says how to find the current state. Measured numbers are from
wilkcomp02 and cite their file: [T1] `planning/measurements/2026-11-10_trial_search_rN1.md`, [T2]
`.../2026-11-11_trial_search_rN1_second.md`, [T3] `.../2026-11-11_search_memory_rN1.md` (and its W-21 addendum).

## 1. Scope and authority

You may: read releases, run searches, export sequences and tables, query a release's catalog read-only.
You may not: build, define or change a release; write anything under a data root or under `/stor`; fetch from
the network (AGENTS §8); start a service; use Docker; edit `~/bio-box` on the pod; run a binary outside the image.

Rules of `AGENTS.md` that bind you (read the sections named):

- Shell (AGENTS §1.8): one command per call; no `&&`, `;`, `||`, pipes; quote paths. A new call keeps no
  variables: in the commands below write every `$NAME` out in full (table in §3). Anything with quotes inside
  quotes goes in a script file, never `python -c` (RUNBOOK §5).
- No compute on the laptop. Everything runs on the pod `wilkcomp02` (reach it with `ssh -o BatchMode=yes` and
  the host in `scripts/pod_pytest.sh`); `bio-box` and the tools run inside the image, in a memory-capped scope
  at nice 19. Read AGENTS §7 whole: it all applies.
- Threads (AGENTS §7.2): 40 by day (08:00 to 20:00 Mon to Fri), 56 at night and weekends only when nothing
  else computes. **The cap is a total over every process of this account**: builds, tests, other agents, you.
  Before a launch, count what runs (§3 step 0); yours plus theirs must stay under the cap of the hour.
- Memory: every command in a scope with `MemoryMax` (200G at most) and `MemorySwapMax=0`; one process under
  300 CPU-minutes. Units: a new name per launch with your own prefix (`qrnap-hmm-1`); never stop a unit you did
  not start. Scratch: `$SCR` = `/tmp/if2453/bio-box/q/<task>` (AGENTS §7.3); never delete outside it.
- Data model (AGENTS §3): coordinates 0-based half-open; strand `+1`/`-1`; an ordinal is never a key and never
  reported without its `feature_uid`; token ids are release-local: never report or store one (DEC-23).
- Numbers (AGENTS §6): every number is *measured* (host, date, command), *derived* (formula, inputs) or *estimate*.
- **Stop and report** (AGENTS §9) when a command behaves unexpectedly, a `STOP` file exists, a check
  disagrees (§6), or a step would break a rule above. Report what you ran and its last lines; do not work around it.

## 2. Orient: which releases exist and what they hold

As of 2026-11-11 (source: the coordinator and `planning/tasks/STATUS.md`): root `/tmp/if2453/r10/root` on
wilkcomp02 (local disk: gone after a reboot) holds

| Release | Content | Built by | Search pack |
|---|---|---|---|
| `rN.1` | 27,258 phage + 2,497 complete prokaryote genomes, stages `S1,P1,P2`; 10,620,684 proteins; NR90, NR50; no Pfam | tree `740d049` (before B-27, B-21, B-39) | `/tmp/if2453/w19/packs` (with `pack_verified.json`; an older copy without it at `/tmp/if2453/w16/packs`), until the operator copies it into the profile's pack directory `/tmp/if2453/bio-box/packs` (card W-22) |
| `rN.2` | the same with Pfam 38.2 | tree `740d049` | none |
| `rN.3` | the same stages, release steps rebuilt with current code (card C-17); being built | current | none |

Find the current state; do not trust the table above. `releases` (card W-23; `src/biobox/cli.py`,
`release_summary`) prints one line per release of the root, with a manifest or only a `release_def`: stages,
genomes, contigs, distinct proteins (from the manifest's file list and Parquet footers; genomes read one column of
`contig_dim`), domain libraries (`pfam-38.2(not P3)` when a library covers some of its stages only), ncRNA
producers, whether it has `member_cluster`, `contig_label` and the `domain_family_tokens` column (`no`: an older
release, slow paths), and the pack it has in the profile's pack directory (`-` for none). `--json` gives one object
per release.

The host profile names the host's one pack directory: its `tmp` plus `/packs`, `/tmp/if2453/bio-box/packs` on
wilkcomp02 (card W-22; `_packs` in `src/biobox/cli.py`, `default_pack_base` in `src/biobox/search/runner.py`).
`releases`, `search`, `nbr` and `pack verify` find a release's pack there with only `--profile`; `--pack-dir` names
another directory. If `releases` shows `-` for a release that should have a pack, look for one with the `find`
below; as of 2026-11-12 `rN.1`'s is still in `/tmp/if2453/w19/packs`, so add `--pack-dir /tmp/if2453/w19/packs` to
every `releases`, `search` and `nbr` command until the operator has copied it.

```bash
bash "$SCR/img.sh" python -m biobox.cli releases --root /tmp/if2453/r10/root --profile "$TREE/profiles/wilkcomp02.yaml"
find /tmp/if2453 -maxdepth 5 -name pack_manifest.json
ls -l /tmp/if2453/images
```

Then read the top of `planning/tasks/STATUS.md` in the repo. A release exists when `manifests/release=<R>.json`
does; a pack of release R is `<pack dir>/release=<R>/<pack id>/`.

Verify the release you will use (60.4 s, 327 MB peak on `rN.1`; measured, 2026-11-10, [T1]; `img.sh` is in §3):

```bash
bash "$SCR/img.sh" python -m biobox.cli release verify --root /tmp/if2453/r10/root --release rN.1
```

The last line is `release rN.1: verified (stages S1,P1,P2, 52 files, ... manifest hash <hex>)` (or `N
problems`); keep the hash. `warning: ... no member_cluster files (built before card B-27)` means slow paths.
`not in manifest: l1/domains/library=pfam-38.2/...` lines on `rN.1` are card C-19 (its stages were scanned
later), not damage: say so in the report. Any other problem: stop.

What a release holds (`catalog_files`, §5): `l1_domains` files mean a domain library (`library=pfam-38.2`);
`member_cluster` and `contig_label` are card B-27's lookup tables; no release has `l1_ncrna`. Which code built
it matters: an older release prints `note: release R has no ... files (built before card B-27); reading ...
instead (the slow path ...)` and is right but slower (about 4 s for `find` on `rN.1`, measured, [T2]). Family
postings and `partners_of` families are complete only after cards B-21 and B-39 (`rN.3`).

## 3. The standard workflow (worked example: Pfam RNA_pol, T7 gp1 control)

Placeholders (write them out in every call): `$SCR` your scratch; `$TREE` a checkout of current code on the pod
(mirror one from the laptop with `sh scripts/pod_pytest.sh <task> tests/r0/test_runbook.py`, which copies the
tree to `/tmp/if2453/ci/<task>`); `$IMAGE` the `.sif` (the trials used `bio-box-test_74b3e05.sif`, sha256
`c30be58557d756095346308c58cc3098af4bbe66404f67ab87be3f64f147c898`, with the source bound in, RUNBOOK §11);
`$DIGEST` `sha256:` plus `sha256sum "$IMAGE"`; `$RUN` a run id.

Once, write `$SCR/img.sh` (runs one short command in the image, capped; the binds of [T2]):

```sh
#!/bin/bash
IMAGE=/path/to/the.sif                      # $IMAGE
TREE=/tmp/if2453/ci/<task>                  # $TREE; the image runs this source
exec systemd-run --user --scope --quiet -p MemoryMax=16G -p MemorySwapMax=0 -p CPUQuota=400% nice -n 19 \
  singularity exec --cleanenv --env PYTHONPATH="$TREE/src:/tmp/if2453/pilot_s1/pydeps" --env OMP_NUM_THREADS=2 \
  --bind /tmp/if2453 --bind /tmp/if2453/r10:/tmp/if2453/r10:ro \
  --bind /stor/scratch/Wilke/if2453/all-the-data/l0/refdb:/refdb:ro --bind /sys/fs/cgroup:/sys/fs/cgroup:ro \
  "$IMAGE" "$@"
```

**0. Preflight.** No `STOP` file (`No such file or directory` is the good answer); load and this account's
CPU (sum the `%CPU` column of `top.txt`, over 100 = threads in use); the units running; free memory and space.

```bash
ls -d /tmp/if2453/r10/root/STOP
uptime
top -b -n 1 -u "$USER" -o %CPU -w 200 > "$SCR/top.txt"
systemctl --user list-units --type=service --type=scope --no-pager
free -g
df -h /tmp
```

**1. Positive control** (required: a protein the search must find). Name it as `ACCESSION:START..END` (GenBank:
1-based, both ends included, as the scientist's record gives it) or `ACCESSION:START-END` (bio-box: 0-based,
half-open; GenBank `3171..5822` is `3170-5822`); do not convert by hand. Any other spelling exits 2 and shows both
forms. Code: `parse_span` and `find` in `src/biobox/cli.py`.

```bash
bash "$SCR/img.sh" python -m biobox.cli find --root /tmp/if2453/r10/root --release rN.1 --accession NC_001604
bash "$SCR/img.sh" python -m biobox.cli find --root /tmp/if2453/r10/root --release rN.1 --accession NC_001604:3171..5822
bash "$SCR/img.sh" python -m biobox.cli find --root /tmp/if2453/r10/root --release rN.1 --protein 2eef240afefddcd2c8ac7660eb50e384
```

Output: a header line, then `accession contig_uid start end strand feature_uid protein_uid span genbank_span` per
gene. With the 0-based form, the second printed `NC_001604 719747a8fef97b2380955d92845861db 3170 5822 +1
39f9ef263e234b012d2671634468546e 2eef240afefddcd2c8ac7660eb50e384` in 4.49 s (measured, [T2], before the two span
columns existed); it now ends `NC_001604:3170-5822 NC_001604:3171..5822`. Exit 2 when it names no gene.

**2. The query.** The scientist's HMM file, or one model from Pfam on the pod (read-only):

```bash
bash "$SCR/img.sh" hmmfetch -o "$SCR/RNA_pol.hmm" /refdb/pfam/38.2/Pfam-A.hmm.gz RNA_pol
```

**3. The pack.** A search needs a pack of the same release in the profile's pack directory,
`/tmp/if2453/bio-box/packs` (§2; `releases` names it). None: `search hmm` exits 2, names the pack directory it
looked in and prints the `pack build` command. Do not build one without asking: `rN.1`'s took 5,282 s, 86.0
CPU-minutes, 5.51 GB (measured, wilkcomp02, unit `w16-pack`, row W-16 of STATUS.md). Prefer a pack with
`pack_verified.json`; without it every search first verifies the whole pack, about 130 s more (derived from
samples, [T3] §4). The pack directory may be read-only for you: a search, `pack verify` and a `pack build` that
finds the pack built write nothing into it (card W-22).

**4. Run** in a `pod_run.sh` unit (4 CPUs, 32G; one long line; flags in `scripts/pod_run.sh` and `_hmm_parser`):

```bash
bash "$TREE/scripts/pod_run.sh" --profile "$TREE/profiles/wilkcomp02.yaml" --image "$IMAGE" --digest "$DIGEST" --mem 32G --threads 4 --bind /tmp/if2453 --bind /tmp/if2453/r10:/tmp/if2453/r10:ro --env PYTHONPATH="$TREE/src:/tmp/if2453/pilot_s1/pydeps" --log "$SCR/logs/qrnap-hmm-1.log" qrnap-hmm-1 search hmm --root /tmp/if2453/r10/root --release rN.1 --hmm "$SCR/RNA_pol.hmm" --bits 30 --neighbors 5 --positive NC_001604:3170-5822 --out "$SCR/out/rnap" --profile "$TREE/profiles/wilkcomp02.yaml" --cpus 4 --mem-gb 30
```

`--out` must be new or empty and outside the root; `--model NAME` picks one model of a multi-model file.

**5. Watch.** The log names the run id (`search run: <id> running; ...`) within about 10 s (measured, [T3]).

```bash
tail -n 5 "$SCR/logs/qrnap-hmm-1.log"
bash "$SCR/img.sh" python -m biobox.cli search status "$RUN" --root /tmp/if2453/r10/root --profile "$TREE/profiles/wilkcomp02.yaml"
bash "$SCR/img.sh" python -m biobox.cli search cancel "$RUN" --root /tmp/if2453/r10/root --profile "$TREE/profiles/wilkcomp02.yaml"
```

`status` prints JSON: `running`, then `ok`, `failed`, `stopped`, `oom`, `timeout` or `output_cap` with `error`,
`rounds`, `final_round`, `stop_reason`. Done when the log ends `pod_run: qrnap-hmm-1 exited 0` after the line
`pod_run: ... usage MemoryPeak=<bytes> CPUUsageNSec=<ns>` (the cost to report). Exit 1: the run ended other
than `ok` and `--out` holds only `recipe.yaml`; exit 2: refused before the run (read the `refused:` line).

**6. Read.** The example ended `2581 genes of 1957 proteins` and `2581 tables`; wall 244.6 s, worker peak
1.78 GiB, unit `MemoryPeak` 4,570,214,400 B, at least 9.53 CPU-minutes (measured, wilkcomp02, 2026-11-11, [T3]
W-21 addendum). 1,931 proteins in phages, 26 in prokaryotes; T7 gp1 367.7 bits, E 4.3e-106 ([T2] §3). The top
NR50 neighbour cluster `fa4cb451607a4dfdd797df9a5e462694` is in 663 of 2,581 loci ([T2]); by position it is
probably T7 DNA ligase, which no file says (rN.1 has no Pfam).

## 4. Output contract (`src/biobox/search/export.py`, `src/biobox/cli.py`)

`--out DIR` holds `README.txt` (`export.README_TEXT`: read it), `recipe.yaml` (`search run DIR/recipe.yaml` reruns
it), `hits.tsv`, `hits.faa`, `neighbors.tsv`, `loci/`. The run is `<profile tmp>/search/<user>/runs/run=<id>/`:
`run.parquet`, `hits.parquet`, `profiles.parquet`, `neighborhoods.parquet`.

- Conventions: `start`,`end` 0-based half-open; a gene across the origin of a circular contig has `end` above
  the contig length. Strand `+1`/`-1`. Uids are 32 lowercase hex digits (16 bytes). Scores in bits to one
  decimal; `seq_evalue_zref` at the fixed Z = 1e8 (DEC-04), two significant digits. Empty cell = no value.
- `hits.tsv`, one row per gene of a hit protein: `accession header contig_len start end strand length
  seq_bitscore seq_evalue_zref pfam_architecture protein_uid nr90 nr50 nr50_rep genome_uid [tier] contig_uid
  feature_uid ordinal tape_version dist_to_contig_end_bp n_records stage kind n_genomes genome_uids`. `tier` appears only when some gene
  has one (none does yet). `stage` is the build stage of the gene (`S1` phages, `P1`, `P2` prokaryotes on `rN.*`),
  `kind` its genome's kind from the stage's inputs (`virus`, `prokaryote`, `plasmid`; empty if none recorded). `pfam_architecture` is `-` without a library. `n_records` counts source records of the
  contig, not hits. `dist_to_contig_end_bp` 0 = partial gene at a linear end; contig length on a circular contig.
  `n_genomes` and `genome_uids` (card C-21c): every genome of the release that carries the contig; `genome_uid` is
  its one representative. The JSON occurrences and the locus tables' facts lines carry the same two.
- `hits.faa`: `>protein_uid length= seq_bitscore= seq_evalue_zref= n_genes= example=ACC:start-end strand=+1
  header`, sequence on one line. `example` is the control's own contig for the control, else the first gene.
- `neighbors.tsv`: `# key<TAB>json` header lines (`n_anchor_proteins`, `n_anchor_loci`, `n_censored_loci`,
  `n_neighbors`, `n_keys`, `policy`, ...), then one table per level. pfam: `level key unannotated library
  model_acc model_name`; nr50/nr90: `level key unannotated rep_uid rep_length model_name` (`key` is the cluster uid; card W-24); then
  `n_loci share_of_loci n_neighbors n_genomes genome_basis offsets modal_offset same_strand_share`. `offsets`
  (`offset:count,...`) and `modal_offset` are relative to the anchor's strand.
- `loci/locus_<accession>_<anchor_feature_uid>.tsv`: `# <header>`, then `# contig_uid=... anchor_feature_uid=...
  anchor_ordinal=... anchor_protein_uid=... tape_version=... policy_hash=... censored=...`, then `offset
  offset_rel ordinal feature_uid start end strand protein_uid length_aa pfam_architecture nr90 nr50
  hit_seq_bitscore delta_bp crosses_origin`. `offset` is coordinate order; `offset_rel` follows the anchor's
  strand (positive downstream). `censored=true`: a linear contig end cut the window.
- `nbr flanks RUN --out DIR` (card W-25, `src/biobox/search/flanks.py`, design `planning/design/W-25_flanks.md`):
  `README.txt`, `groups.tsv` (`# key<TAB>json` settings; `group n_loci n_leaves n_leaves_in_tree sampled
  n_context_types`), `context_types.tsv` (one row per locus: `feature_uid group context_type leaf in_tree contig_uid
  accession start end strand anchor_protein_uid anchor_nr90 anchor_nr50 stage kind genome_uid n_genomes genome_uids
  tax_order cut_up cut_down n_up n_down`), `conservation.tsv` (`group context_type key unannotated n_loci_type n_loci
  share_of_loci n_neighbors n_genomes offsets modal_offset same_strand_share`, offsets along the anchor's strand),
  `groups/<group>/{leaves.tsv,distance_j.tsv,distance_a.tsv,positions.tsv,tree.nwk}` (leaves named by the least
  `feature_uid` of their loci), `gff/` with `--export gff`. J = 1 − Jaccard of flank tokens (id, side, strand
  relative to the anchor), A = 1 − Jaccard of adjacent token pairs; a gene with no id at `--level` is a token of its
  own. Defaults: `--window 5g --level nr50 --group-by anchor-nr90 --max-loci 2000 --cut 0.5 --distance j`; only gene
  windows. `sampled=true`: only the 2,000 leaves with the least names are in the tree; the rest have `in_tree=false`
  and no context type: say so. On phages prefer `nr50`: at `pfam` unannotated genes split nearly every locus into
  its own type (measured, `planning/measurements/2026-11-12_flanks_rN3.md`). Without `--profile`, give `--runs` and
  `--pack-dir`.
- Name a gene by `feature_uid` (or accession and range), a protein by `protein_uid`, a family by `model_acc`, a
  cluster by its uid **with the release** (cluster uids change between releases, 0-07); never a bare ordinal or a token.
- Every command given `--out` (and `search hmm` for each file in its DIR) prints to stderr `wrote PATH: N rows`
  (`sequences` for FASTA, `records` for JSON lines; `N tables, M rows` for `loci/`): report those counts.
- One sequence without exporting the release: `seq --protein UID [UID ...]` (or `--proteins-file`), `seq --contig
  ACC:START..END` or `ACC:START-END` or a contig uid, `--strand -1`; read by key (`export.protein_seqs`,
  `cli.contig_slice`). Headers: `>uid release=R length=N`, or `>contig_uid release=R accession=A start=S end=E
  strand=+1 contig_len=N` (0-based, half-open).

```bash
bash "$SCR/img.sh" python -m biobox.cli seq --root /tmp/if2453/r10/root --release rN.1 --protein 2eef240afefddcd2c8ac7660eb50e384
bash "$SCR/img.sh" python -m biobox.cli seq --root /tmp/if2453/r10/root --release rN.1 --contig NC_001604:3171..5822 --strand -1 --out "$SCR/out/t7_gene1_rc.fna"
```

## 5. Questions the one command does not answer: the catalog

The catalog `catalog/release=<R>.duckdb` holds views over the release's files and macros
(`src/biobox/catalog.py`: `macro_sql`, `pair_macro_sql`, `partners_of_sql`; tables in
`src/biobox/registry/schemas.py`). Query it with `bio-box sql` (card W-23; `_sql` and `sql_connection` in
`src/biobox/cli.py`): it opens the catalog with `catalog.open_catalog` (a plain `duckdb` connection fails on the
sequence macros and `partners_of`) from a copy in a scratch directory under the profile's `tmp` (or the system's
temporary directory), where DuckDB also spills, never under the root. It runs only statements DuckDB's parser
calls SELECT (also `SHOW`, `DESCRIBE`, `SUMMARIZE`, `FROM ...`), on a connection that cannot write: the catalog is
attached read-only and external access is off, with only the release's and the runs' files readable, so `COPY ...
TO`, `CREATE`, `ATTACH`, `INSTALL`, `LOAD`, `SET` and `PRAGMA` are refused (exit 2). Each statement's rows go to
stdout (a TSV table with its column line, a blank line between statements; `--format json` for one statement, one
object per row) or to `--out FILE` outside the root; stderr has, per statement, its row count and the row groups it
read per table (`catalog.profile_query`). `--runs` attaches closed runs (a run directory, or a runs directory: its
runs of this release); with `--profile` and no `--runs`, your runs directory. Write the statements in a file
(`;`-separated) or give one with `-e`:

```bash
bash "$SCR/img.sh" python -m biobox.cli sql --root /tmp/if2453/r10/root --release rN.1 --profile "$TREE/profiles/wilkcomp02.yaml" "$SCR/q1.sql"
bash "$SCR/img.sh" python -m biobox.cli sql --root /tmp/if2453/r10/root --release rN.1 -e "SELECT key, value FROM catalog_info ORDER BY key"
```

Uids are BLOBs: write a literal as `unhex('2eef...')`; `sql` prints them as hex. One example per object
(numbers measured on `rN.1`, [T1] §1); put a `LIMIT` on anything that may return many rows (the rows of a statement
are held in memory before they are printed):

```sql
SELECT key, value FROM catalog_info ORDER BY key;              -- release_id, manifest_hash, root, tape_version, stages
SELECT "table", count(*) AS files, sum("rows") AS n FROM catalog_files GROUP BY 1 ORDER BY 1;  -- what it holds
SELECT * FROM slots_of(unhex('2eef240afefddcd2c8ac7660eb50e384'));   -- protein_occurrence: genes of a protein
SELECT protein_of(unhex('2eef240afefddcd2c8ac7660eb50e384'));        -- its sequence
SELECT * FROM domains_of(unhex('2eef240afefddcd2c8ac7660eb50e384')); -- l1_domains rows (rN.2, rN.3)
SELECT level, cluster_uid, rep_uid FROM member_cluster WHERE member_uid = unhex('2eef240afefddcd2c8ac7660eb50e384');
  -- rN.3; on rN.1 and rN.2 the same columns from cluster_member (every row group)
SELECT contig_uid, accession, original_header, n_records FROM contig_label WHERE accession = 'NC_001604';
SELECT contig_len_of(unhex('719747a8fef97b2380955d92845861db')),
       seq_get(unhex('719747a8fef97b2380955d92845861db'), 3170, 3200, 1);  -- 39937, ATGAACACGATTAACATCGCTAAGAACGAC
SELECT genome_uid, is_circular, tier, tax_order FROM contig_dim_now() WHERE contig_uid = unhex('719747a8fef97b2380955d92845861db');
SELECT count(*) AS postings, count(DISTINCT contig_uid) AS contigs
  FROM token_rows(token_of('nr50', '5ea5af7357bbdd6d2ad9e3ddfe56c605'));  -- 1,064 on 1,008
SELECT d.level, d.key, p.n_loci, p.n_pairs FROM partners_of(token_of('nr50', '5ea5af7357bbdd6d2ad9e3ddfe56c605'), 5) p
  JOIN token_dict d ON d.token_id = p.partner_token AND d.release = catalog_release() ORDER BY p.n_loci DESC LIMIT 20;
SELECT q1_count(token_of('nr50', '5ea5af7357bbdd6d2ad9e3ddfe56c605'), token_of('nr50', 'fa4cb451607a4dfdd797df9a5e462694'), 5, NULL);
SELECT * FROM fusions_of(token_of('family', 'PF00940.25'), token_of('family', 'PF00000.0')) LIMIT 20;
SELECT o.stage, count(DISTINCT o.protein_uid) AS proteins, count(*) AS genes, count(DISTINCT o.contig_uid) AS contigs,
       count(DISTINCT d.genome_uid) AS genomes
  FROM protein_occurrence o LEFT JOIN contig_dim_now() d ON d.contig_uid = o.contig_uid
  WHERE o.protein_uid IN (SELECT target_uid FROM search_hit WHERE target_kind = 'protein') GROUP BY 1 ORDER BY 1;
SELECT run, status, recipe_hash, pack_id, ran_against, z_prot, wall_s, cpu_min_max_pid, peak_rss_bytes FROM search_run;
SELECT DISTINCT lower(hex(target_uid)) AS "#protein_uid" FROM search_hit WHERE target_kind = 'protein' ORDER BY 1;
```

The last three need the run attached (`--runs` with the run directory). The stage split is also in `hits.tsv`
(columns `stage` and `kind`); the query adds genomes and contigs. The last, alone in its own file `$SCR/uids.sql`
and run with `--out "$SCR/uids.txt"`, writes a file `nbr --anchors-file` reads (its column line starts with `#`, so
it is skipped as a comment; `sql` refuses `COPY ... TO`). Keys for `token_of`: a cluster uid as lowercase hex at
`nr50`/`nr90`, a `model_acc` as `l1_domains` stores it at `family` (take it from `domains_of` of a known
protein; `PF00000.0` above is a placeholder). `token_of` turns a stable key into a token inside the query; the
token id itself never leaves the query.

**Pruning rule (AGENTS §4.4).** A lookup by a key a table is sorted by reads one or two row groups per key; `sql`
prints on stderr what each statement read: report it with any number from SQL. Cheap: by `protein_uid` (`slots_of`,
`protein_of`, `domains_of`, `member_cluster`), by `contig_uid` (`contig_label`, `contig_len_of`, `seq_get`,
`contig_dim_now()`, `tapes_of`), by stable key (`token_of`, `token_rows`, `partners_of` with few loci).
Expensive (every row group; seconds to minutes at 10 Gbp, an estimate; say so): by accession or header
(`contig_label` is sorted by contig, W-15), by model in `l1_domains` (sorted by protein), any lookup before B-27
(slow path), joins with thousands of keys, whole-table aggregates. `partners_of` refuses a token with more than
100,000 loci. Not supported: taxonomy and tiers (`contig_dim.tax_order` and `tier` are null: no such step yet,
`src/biobox/dim.py`), enrichment or p-values (cards F-01 to F-03), ncRNA (E-07), clusters across releases (0-07).

## 6. Verify before you report

Both trials checked the answer with a tool independent of the runner: `hmmsearch` (HMMER 3.4, in the image) on
every protein of the release, at the same Z and cutoff. Export (88.2 s, 2.28 GB peak, 3,712,392,579 B for
`rN.1`; measured, [T1]), then search (26.3 s wall, 2.7 CPU-minutes at `--cpu 6`; measured, [T1]):

```bash
bash "$SCR/img.sh" python -m biobox.cli export proteins --root /tmp/if2453/r10/root --release rN.1 --threads 2 --out "$SCR/proteins.faa"
bash "$SCR/img.sh" hmmsearch --cpu 4 -Z 1e8 --domZ 1e8 -T 30 --noali --tblout "$SCR/check.tbl" -o "$SCR/check.out" "$SCR/RNA_pol.hmm" "$SCR/proteins.faa"
bash "$SCR/img.sh" python "$SCR/check.py" "$SCR/out/rnap/hits.tsv" "$SCR/check.tbl" 2eef240afefddcd2c8ac7660eb50e384
```

```python
"""$SCR/check.py HITS_TSV TBLOUT CONTROL_UID: proteins and scores of a run against hmmsearch --tblout."""
import csv
import sys

ours, theirs = {}, {}
with open(sys.argv[1], newline="") as fh:
    for r in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
        ours[r["protein_uid"]] = float(r["seq_bitscore"])
for line in open(sys.argv[2]):
    if not line.startswith("#"):
        f = line.split()
        theirs[f[0]] = float(f[5])  # target name, ..., full-sequence score
both = ours.keys() & theirs.keys()
print("run", len(ours), "hmmsearch", len(theirs), "only run", len(ours.keys() - theirs.keys()),
      "only hmmsearch", len(theirs.keys() - ours.keys()))
print("max score difference", max((abs(ours[u] - theirs[u]) for u in both), default=None))
print("control", sys.argv[3], ours.get(sys.argv[3]), theirs.get(sys.argv[3]))
```

Expect what both trials found ([T2] §3): the same protein set (0 only in either), scores equal to 0.1 bit, the
control in both. Then split the hits by source with the `stage` and `kind` columns of `hits.tsv` (phages are stage
`S1`, prokaryotes `P1`, `P2`), and genomes per stage with the stage query of §5. Neighbours can be checked the same way: an independent SQL window (each contig's genes
by start, end, strand) gave the same 2,581 locus tables ([T2] §3). **If a check disagrees, report both answers
and the commands; do not choose one.** Remove the 3.7 GB FASTA when done: `rm "$SCR/proteins.faa"`.

## 7. Report to the scientist

Lead with the answer in one or two sentences, then the evidence:

- Every number labelled *measured* (wilkcomp02, date, command) or *derived* (formula, inputs). Counts with their
  unit, split by source: proteins (distinct sequences), genes, contigs, genomes (`contig_dim.genome_uid`, or every
  genome that carries a hit contig: `genome_uids` of `hits.tsv`, the `genome_member` view; they differ only where a
  contig is shared, card C-21c).
- To reproduce: the release and its manifest hash, the run id, recipe hash and pack id (`search_run`), the image
  sha256, the commit of `$TREE` (the laptop's `git rev-parse HEAD` when mirrored), the query file's `sha256sum`.
- The independent check (§6) and its result; what `sql` says on stderr each SQL statement read.
- What the release cannot tell: no Pfam on `rN.1` (clusters have uids, no names); no ncRNA; no taxonomy; no
  tiers; neighbour counts are counts, not significance (no enrichment statistic exists).
- Where the outputs are (`$SCR` on the pod's `/tmp`: gone after a reboot; offer to copy them to the laptop).

## 8. Traps

Each is stated where it applies above; this is the checklist. The card or decision that tracks it is in brackets.

- Releases built before B-27 (`rN.1`, `rN.2`) answer by slow paths with a `note:` line: right, slower (W-17, C-17).
- Before B-21/B-39 (`rN.1`, `rN.2`) a slot posts only its canonical family: an empty `fusions_of` there is not
  "never fused", and `partners_of` families are canonical only (B-21, B-39).
- `nbr summarize --by pfam` without a library drops the level with a note; a pfam row `unannotated` at share 1.0
  means code older than W-20 (W-20).
- Old releases fail `release verify` once a library is scanned on their stages (C-19).
- A pack is per release; only `rN.1` has one. A pack built before W-19 is verified in full by every search (W-19).
- Packs live in the profile's pack directory, `/tmp/if2453/bio-box/packs`; `rN.1`'s is still in
  `/tmp/if2453/w19/packs` until the operator copies it (§2; W-22). `pack verify` writes nothing into a pack
  directory (its scratch is in `<profile tmp>/search/<user>/work`, W-22); never run `pack build` without asking.
- `/tmp` is local disk: roots, packs, runs and your outputs vanish on a reboot (RUNBOOK §4).
- A contig shared by genomes of two stages belongs to the first; the later genome is one contig short in
  `contig_dim`. Not in `rN.*`; it is in the 50 Gbp set (DEC-26, C-21).
- `search run` in a `pod_run.sh` unit needs `--wait` (W-13). A search or `nbr` command without `--profile`, or
  `--pack-dir` and `--runs`, uses `/tmp/atd`, which is not ours (W-03).
- `search hits` without `--format` shows 50 rows (`--top_n 0` for all) (W-20).
- `--positive 'TEXT'` (a header substring) nearly always matches many genes and is refused with a list (W-15).
- Cluster uids, `ACC:START-END` against `ACC:START..END` (two dots: GenBank, 1-based), offsets and the thread
  total: §4, §3 step 1, §1 (0-07, W-20, W-23, AGENTS §7.2).
- `sql` refuses anything but SELECT, `COPY ... TO` included: write results with `--out` (W-23).

## 9. Decision table

| The scientist asks | Do |
|---|---|
| Homologs of a family (an HMM) | §3 on the release with a pack; check with §6 |
| Homologs of a sequence they give | Write `$SCR/recipe.yaml` (below); `pod_run.sh ... search run "$SCR/recipe.yaml" --root ... --release rN.1 --profile ... --cpus 4 --mem-gb 30 --wait`; then the exports below; `tool: jackhmmer` with `max: 3` to iterate |
| What is next to them | `neighbors.tsv` and `loci/`. With Pfam names: summarize the same proteins on `rN.2`/`rN.3` (no pack needed, below) |
| Which hits share a genomic context | `nbr flanks "$RUN" ... --out "$SCR/out/flanks"` (§4); report context types and `conservation.tsv` rows, and the level and cut used |
| Are families A and B ever fused | `rN.3`: `fusions_of` (§5). `rN.2`: `pfam_architecture` of a family-A search's hits, or a self-join of `l1_domains` on `protein_uid` (a full scan: say so) |
| Export sequences | Hits: `hits.faa`. Some proteins or a contig region: `seq` (§4). Every protein: `export proteins` (§6) |
| How many genomes have it | The stage split query of §5 (`count(DISTINCT genome_uid)`); define genome as `contig_dim.genome_uid` (one per contig) or count every carrier through `genome_member` or `hits.tsv`'s `genome_uids`, and say which. By source only: `stage`, `kind` of `hits.tsv` |
| What releases exist | `releases` (§2) |
| Taxonomy of the hits | Not available (no taxonomy in any release yet). Give `accession` and `header` of `hits.tsv`, labelled as source headers |
| Compare two releases | `release verify --root ... --release rN.3 --previous rN.2`; run the same search on both (a pack each) and compare `protein_uid` sets |
| Is it significant / enriched | No statistic exists (cards F-01 to F-03): report counts and shares only, and say so |
| Run on the larger build | Check STATUS.md for the root and a pack; no pack: ask before building one; the search memory does not grow with the release (W-19) |

```yaml
recipe: user_seq_phmmer
written_against: rN.1
seed: {sequence: MNTINIAKNDFSDIELAAIPFNTLADHYGE}      # the scientist's protein sequence
target: unique
tool: phmmer
iterations: {max: 1, inclusion_bits: 30}
context: {window: 5g}
controls:
  positive: [{protein_uid: 2eef240afefddcd2c8ac7660eb50e384, label: "NC_001604:3170-5822"}]
```

A positive control is `protein_uid` or `sequence` (the exact sequence of a protein in the release); recipe
format in `src/biobox/search/recipe.py`. Exports of a run, and a summary of given proteins on another release:

```bash
bash "$SCR/img.sh" python -m biobox.cli search hits "$RUN" --root /tmp/if2453/r10/root --profile "$TREE/profiles/wilkcomp02.yaml" --format tsv --out "$SCR/out/seq/hits.tsv"
bash "$SCR/img.sh" python -m biobox.cli search hits "$RUN" --root /tmp/if2453/r10/root --profile "$TREE/profiles/wilkcomp02.yaml" --format faa --out "$SCR/out/seq/hits.faa"
bash "$SCR/img.sh" python -m biobox.cli nbr "$RUN" --root /tmp/if2453/r10/root --profile "$TREE/profiles/wilkcomp02.yaml" --window 5g --format tsv --out "$SCR/out/seq/loci"
bash "$SCR/img.sh" python -m biobox.cli nbr summarize "$RUN" --root /tmp/if2453/r10/root --profile "$TREE/profiles/wilkcomp02.yaml" --window 5g --by nr50,nr90 --format tsv --out "$SCR/out/seq/neighbors.tsv"
bash "$SCR/img.sh" python -m biobox.cli nbr summarize --anchors-file "$SCR/uids.txt" --root /tmp/if2453/r10/root --release rN.2 --profile "$TREE/profiles/wilkcomp02.yaml" --window 5g --by pfam,nr50 --format tsv --out "$SCR/out/neighbors_pfam.tsv"
bash "$SCR/img.sh" python -m biobox.cli nbr flanks "$RUN" --root /tmp/if2453/r10/root --profile "$TREE/profiles/wilkcomp02.yaml" --window 5g --level nr50 --out "$SCR/out/flanks"
```

`uids.txt` is the last query of §5, with `AND run = '<run id>'` added to keep one run, written with `sql --out`;
`--profile` attaches your runs of the release (protein uids are content hashes, the same in every release over
these stages):

```bash
bash "$SCR/img.sh" python -m biobox.cli sql --root /tmp/if2453/r10/root --release rN.1 --profile "$TREE/profiles/wilkcomp02.yaml" "$SCR/uids.sql" --out "$SCR/uids.txt"
```

If the summary reports 0 anchor loci, pass `--tape-version` with the release's `tape_version` from `catalog_info`.
