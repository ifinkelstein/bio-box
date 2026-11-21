# C-49: what the NAS `tower` should hold, and how data reaches it

A proposal for the PI (card C-49, DEC-60 (3), STATUS open item 19), 2026-11-20. No code here except the relay script
the card asks for (`scripts/orchestrator/relay_to_tower.sh`). Rates and capacities are in
`planning/measurements/2026-11-20_tower_relay.md`. Numbers carry AGENTS §6's labels; a number copied from another note
names it. "From the code" means read in the tree of this branch, not run.

## 1. In short

- `tower` can hold a great deal (18 TB free, *measured*) and can compute almost nothing for us (8 cores, 31 GB, on a
  home network the pod cannot reach). Every file reaches it through the laptop at about 9.70 MB/s (*measured*), a
  terabyte in about 29 hours (*derived*), and only while the laptop is awake at home. The other direction is not
  measured.
- **No file that a live release lists, and no file of a built stage, can live only on `tower`.** Every reader opens
  the root's files by path on the pod; a build of the next stage reads every earlier stage's index files and silently
  leaves out a stage whose files are absent (section 3.3). So `tower` cannot be a working tier for the corpus; it can be
  a copy, and a cold store for what no command reads.
- **Recommended: tier A now (the second durable copy, running today) and the copying half of tier B** (source
  archives, staged inputs and the deferred-data lists to `tower` as copies; the PI decides later whether `/stor`'s
  copies go). Tier C (old releases' large files off the pod) frees 61 to 76 GB today (*derived*) for a set of code
  changes and a change of the rule that a root is append-only: not now; decide when the pod's free space calls for
  it.

## 2. How a file reaches `tower`

The route the PI chose: pod to laptop to `tower`, with nothing landing on the laptop (26 GiB free, *measured*). No key
on `tower` reaches the pod (the PI, 2026-11-20). `scripts/orchestrator/relay_to_tower.sh` (this card) does it:

- input: a source directory on the pod, a destination directory on `tower`, and a list of `<path><TAB><sha256>` lines
  or a `SHA256SUMS` file; every path and checksum is checked before the first remote call;
- per file: skip if in place (never overwrite); stream `ssh pod cat` into `ssh tower 'tee .partial | sha256sum'`;
  compare the hash with the record and the size with the pod's; `mv -n` into place; on a mismatch remove the
  `.partial` and try again (three tries by default); at the end write or extend `SHA256SUMS` on `tower` with the files
  verified in this run; append to a log; exit 1 if any file failed;
- `--streams N` (default 1), `--dry-run` (what and how many bytes; writes nothing), `--reverify` (read every file back
  on `tower` as well).

**Why the stream is hashed on arrival.** The first mirror (`scp -3`, then `sha256sum` of the file on `tower`) read
every byte back from the array at 7.5 to 9.2 MB/s while copies ran, so two streams together delivered 5.30 MB/s from
start to verified against 9.70 MB/s for one stream's copy (*derived* from its log; the measurement note §2). Hashing
the stream as `tee` writes it costs one process on `tower` and no extra read. What it checks: the bytes `tower`
received equal the record, and the file on disk has the pod's size. What it does not check: that the disk later
returns the bytes it was given; a silent write error would pass. `--reverify` covers that at the cost of reading every
byte again (100 GB: about 2.8 hours more, *derived*); the array's parity check covers it for the disks. I chose the
stream for the default because one copy of every file is already on `/stor` (tier A), so a fault found later costs a
copy, not data; for a file whose `/stor` copy the PI means to remove (tier B), run `--reverify` once before.

**Limits of the script.** A file interrupted mid-stream starts again from zero: a 524.8 GB archive
(`jgi_iuig_genomes_20261118.zip`, *(card C-48d)*) needs about 15 hours uninterrupted at 9.70 MB/s (*derived*), more than
the laptop's "about 7 hours" on 2026-11-20. Resuming by byte offset is possible (`tail -c +K` on the pod, `cat >>` on
`tower`) but the hash would then need a read-back of the first part; not done here (open question 4). Many small files
cost about six `ssh` calls each through the jump host: for a directory of many files, make one tar on the pod first
(in `/tmp`, by a card that may run on the pod) and relay the tar.

## 3. Each kind of file

For each kind: (a) which commands read it and how often, from the code; (b) whether a search, `nbr`, `find`, `export`,
`release verify`, a stage build or a release build needs it on the machine that computes; (c) what breaks, exactly, if
it is only on `tower`; (d) its size.

### 3.1 Source archives (`/stor/work/.../jgi_20261118/`, `/stor/scratch/.../provisional/`)

- (a) Read by the scripts that make a stage's inputs (`scripts/stress/build_stage_pod/make_inputs.py`, the `c36_*`,
  `c48c_*` and `c48d_*` drivers), once per stage they feed, and again when a deferred remainder is picked up
  (`planning/CATALOG_GOAL.md` §4: "no fetch needed: the same member of `jgi_img_vr_v4.1_hc_20261118.zip`, by id").
  No `bio-box` command reads them; no manifest lists them.
- (b) None of search, `nbr`, `find`, `export`, `release verify`, a stage build or a release build. Only input making.
- (c) Input making for a stage not yet staged, or a remainder not yet taken, cannot start until the archive is back on
  the pod: relayed back (rate not measured) or fetched again from the provider (JGI's need the PI's session, card
  C-48). Nothing else breaks.
- (d) JGI archives: 1,137,921,831,184 bytes (*measured*, `2026-11-20_audit_of_1117_1120.md` §2). `provisional/` on
  `/stor/scratch`: 163 GB, stage inputs and older roots together (*measured*, STATUS open item 19; not split here).
  Per Gbp: not meaningful (compressed archives).

### 3.2 Staged stage inputs (`/tmp/if2453/s2/inputs/<stage>/`)

- (a) `build-stage`'s `inputs` and `ingest` steps (`stage.py`: `--genomes`, `--objects`; `input_bytes` and
  `object_sizes` read the objects for the estimate at each launch); a resumed ingest reads them again. After ingest the
  root holds the sequences and their source sha256.
- (b) A stage build until its ingest is published; nothing after.
- (c) Before ingest: the build cannot start. After: nothing breaks; they are needed only to rebuild the stage from
  scratch, and can be made again from the archives "at a cost of hours" (audit §5).
- (d) 38 GB for `V1`, `PL2`, `G1` and `C1` to `C6` (*measured*, audit §5), on `/tmp` only.

### 3.3 A stage's files in the root, by layer

Partitioned by `stage=` (registry: `src/biobox/registry/schemas.py`).

- **`l0seq`** (`seq_shards`, `seq_locator`, `contig_index`, `contig_occurrence`, `contig_label`, `genome_input`,
  `igr`, `protein_seq`).
  - (a) Every later stage's ingest reads **every** stage's `seq_locator` and `contig_index` found on disk
    (`store.locator_files`, `store.index_files`: a glob over `stage=*`), so a contig is stored once (DEC-26). Gene
    calling of every later stage is given every stage whose `seq_locator` exists (`stage.with_seq_stages`,
    `store.stages_on_disk`) and reads the sequences of those that own a contig of its genomes. Every release's
    clustering reads every listed stage's `protein_seq`; the pack build reads `protein_seq`, `igr` and the sequence
    shards; `seq`, `export`, `find` and `search --positive` read shards and labels through the catalog.
  - (b) Stage build: yes, for every earlier stage. Release build, pack build, `find`, `seq`, `export`, `nbr`,
    `release verify`: yes. A search over a built pack reads the pack, then the catalog for its hits.
  - (c) **Silent.** A stage whose `l0seq/seq_locator/stage=<S>/` is absent is not in `store.stages_on_disk`: the next
    stage's ingest does not compare against it and stores its contigs again, and gene calling does not see it. No
    step fails. The release over both then holds one contig in two stages, which `dim.owner_problem` was written to
    refuse only when a listed stage lacks its owner. Every other reader fails loudly (a missing file).
  - (d) Stage `H1`: 0.309 GB per Gbp at depth `stored` (*derived*, `planning/design/C-38_import_depths.md` N1), almost
    all of it the sequence shards: the other depth-0 files are about 0.9 MB per Gbp for this stage of few contigs
    (*derived*: N1's 561,924,885 bytes less the shards' 560,348,485, C-38 §6, over 1.82 Gbp). The `l0seq` layer of
    `S1` and `P1` to `P6`: 28.24 GB over 52.1 Gbp, 0.54 GB per Gbp (*derived*: `l0seq.tar`, 28,241,715,200 bytes,
    STATUS 2026-11-14; 52.1 Gbp *measured*, card C-33, as cited in C-38 N22).
- **`l1`** (`feature`, `feature_protein`, `protein_occurrence`, `domains`, `crispr`, `ncrna`; plans `genecall_plan`,
  `domain_scan`).
  - (a) The release's token and postings steps (tapes and domains); `nbr` and `domain_state` (`catalog`, card
    C-38b); a later stage's Pfam scan copies the rows an earlier stage measured (`domains.sources`, card B-47b) from
    every stage on disk with a scan plan.
  - (b) Release build, `nbr`, exports with `pfam_architecture`, `release verify`: yes. Stage build: only as the
    source of copied Pfam rows.
  - (c) Release build and readers fail on the missing file. A later stage's scan finds no source and scans the
    proteins again (Pfam about 1,550 CPU-minutes per Gbp, *derived*, C-38 N14): a cost, not a wrong answer.
  - (d) Stage `H1` at depth `called`: 0.714 GB per Gbp in all its layers, and 0.018 GB per Gbp more for Pfam on a
    phage stage (*derived*, C-38 N3, N4). Unbinned contigs of 1.5 kb and more: 0.85 GB per Gbp at `called`
    (*measured*, `2026-11-18_spire_pilot.md`).
- **`l3/tapes`** (`contig_tapes`, by stage and tape version). (a) The release's token step; `nbr` and `find` through
  the catalog. (b) Release build, `nbr`, `find`, `release verify`. (c) Fails on the missing file. (d) Inside the
  figures above.
- **Records** (`runs/stage=<S>/`, `ledger/stage=<S>/`, `l0/l0_manifest`). (a) `build-stage` on every launch (resume),
  the release train's choice of complete stages (`stage.stage_incomplete`), `built_by`. (b) Stage and release
  builds. (c) A stage without its run records looks unbuilt: the train leaves it out. (d) Small: all records,
  catalogs, dims and `l1` of the first seven stages together were 9.86 GB (`l1_and_records.tar`, *measured*, STATUS
  2026-11-14), the records 26 MB of it (STATUS, the same).

### 3.4 A release's `release=` files

`cluster_member`, `cluster_member_aln`, `member_cluster`, `token_dict`, `tape_tokens`, `token_postings`,
`contig_dim`, `genome_dim`, `ledger/release=<R>`, `release_def`, the manifest and the catalog.

- (a) Every query on the release: `open_catalog` (from `generate_catalog`, whose views name each listed file by an
  absolute path under the build root, `catalog_info.root`; `open_catalog(root=...)` moves the sequence shards only).
  `release verify` reads every listed file (`catalog.verify_manifest`, sha256 of each). The release train
  (`scripts/stress/build_stage_pod/release_train.sh`, `verify_rel`) verifies **every** release with a manifest before
  and after each cut, so today a cut reads the seven releases before it and the eight after. The next release takes this one's clustering when
  its proteins are the same (`l2/cluster.reuse_sources`, `reuse`, card B-49), by copying its cluster files. Novelty
  checks of a new release read the previous one's clusters (`c35_check.py`).
- (b) Search: the hits' annotation and `--positive` go through the catalog; `nbr`, `find`, `export`, `release
  verify`: yes. Release build: the previous release's cluster files when its clustering is reused.
- (c) `release verify` prints `missing: <path>` and exits 1; the release train then refuses to cut ("an earlier release
  did not verify before launch"). A query that scans a table with an absent file fails in DuckDB with an I/O error
  (from the code: a view reads its whole file list; not run in this card). `--release latest` (`catalog.latest_release`)
  still picks the release if it has a manifest and a verified pack. A reuse of its clustering fails, and the next cut
  clusters again: 15,312.59 CPU-minutes for `r26.11.1` (*measured*, `2026-11-16_release_r26.11.1.md`).
- (d) `r26.11.1`: 930 files, 20,907,311,935 bytes added to the root by the release (*measured*, the same note §5);
  over about 59.84 Gbp (*derived*: 52.1 + 1.82 + 5.92, the bases of `rM.2`, `H1` and `PL1`'s records) that is 0.35 GB per
  Gbp per release. `rM.2`: at least 17 GB (*measured*, DEC-52). The clustering files of the seven releases before
  `r26.11.1`: 27,176,785,920 bytes (`l2.tar`, *measured*, STATUS 2026-11-14); their tokens and postings are inside
  `l3.tar`'s 15,692,646,400 bytes with the stages' tapes, not split.

### 3.5 Search packs (`/tmp/if2453/bio-box/packs/release=<R>/<pack id>/`)

- (a) Every `search` (`runner.prepare`: `pack.find` in the profile's pack directory; `runner.check`:
  `pack.search_check`, sizes and the build's record; `runner._stage`: copies the target's shards into the run's work
  directory and checks their sha256); `nbr` on a run (card W-24 refuses without the run's pack); `bio-box releases`;
  `--release latest` for a search.
- (b) Search and `nbr`: yes. `release verify` and builds: no. A pack is a cache of its release: `pack build` makes it
  again, byte for byte (card W-35).
- (c) `search` refuses with the command that builds the pack; `latest` passes over the release for a search.
- (d) `rN.3` 7,755,186,114 bytes, `rM.2` 25,794,264,796, `r26.11.1` 29,293,903,636 (*measured*,
  `2026-11-14_pack_build_in_parts.md`, `2026-11-16_release_r26.11.1.md`); 62.84 GB in all (*derived*). Per Gbp: 0.49
  for genome stages (*derived*: 29.29 / 59.84), 1.01 for unbinned contigs (*measured*, SPIRE pilot). A build:
  498.86 CPU-minutes for `r26.11.1` (*measured*), 112 minutes of wall time at 12 jobs (*derived* from its unit's times).

### 3.6 The durable tars (`/stor/work/Wilke/if2453/all-the-data/root_backups/`)

- (a) Read by no command. They are the root's restore: `tar -C <new root> -xf` each in order (STATUS 2026-11-14), after
  a reboot clears `/tmp` (it did on 2026-11-17).
- (b) No. (c) Nothing breaks; they are copies. (d) Seven tars, 106,785,105,920 bytes (*measured*, audit §2): the whole
  root at `r26.11.1` (*derived*: 85,874,224,242 bytes before the release plus 20,907,311,935 added, within 4 MB of the
  tars, which add headers). Each later stage or release adds one tar of its new files.

### 3.7 The "left for later" lists (`planning/CATALOG_GOAL.md` §4)

- (a) Read by no command; a later card picks a remainder up by these identifiers (DEC-59 (3)).
- (b) No. (c) Nothing breaks now; losing one loses the record of what to fetch later.
- (d) Not measured here. Two are on `/tmp` only and so are lost at the next reboot:
  `/tmp/if2453/s2/inputs/V2/V2_left_for_later_hc.tsv.gz` and `V2_left_for_later_v4.0_putative.tsv.gz` (sha256 in
  `CATALOG_GOAL.md`), and `PL2/skipped_records.tsv`. Files under 64 MB cannot go to `/stor` (AGENTS §7.3); `tower`
  has no such rule.

### 3.8 User results (`<bulk>/results/` on `/stor/scratch`, DEC-42)

- (a) `bio-box results list`, `results get` (card W-28, W-33), by the user who archived the run.
- (b) No build reads them. (c) `results list` and `get` read the profile's `bulk` only; a bundle only on `tower` is
  not found. (d) One bundle so far, 164,628,480 bytes (*measured*, RUNBOOK Appendix A, W-28).

### 3.9 Which command needs which kind on the pod

| Command | Archives | Staged inputs | Stage files | `release=` files | Pack | Tars, lists, results |
|---|---|---|---|---|---|---|
| `search` (and `search hits`) | no | no | through the catalog for hits and `--positive` | yes | yes | no |
| `nbr` | no | no | yes (tapes, domains, labels) | yes | yes (the run's) | no |
| `find`, `seq`, `export` | no | no | yes | yes | no | no |
| `release verify` | no | no | yes (every listed file) | yes | no | no |
| stage build | no | until ingest | **every earlier stage's index**, and shared contigs' sequences | no | no | no |
| release build | no | no | every listed stage | the reused release's clusters | no | no |
| pack build | no | no | yes | yes | writes it | no |
| input making | yes | writes them | `contig_index` of stored stages, for the overlap reports | no | no | writes the lists |

## 4. Three tierings

At 9.70 MB/s (*measured*, one stream), relay hours are bytes / 9.70e6 / 3,600 (*derived*). The other direction is not
measured; where a way back crosses it, its time is given at the same rate and marked as a lower bound only if the
home upload is the slower link (open question 1).

### Tier A: `tower` as the second durable copy only (running today)

- **What moves**: nothing moves; copies. The seven tars (106.79 GB) and every later `added_<stage>.tar` and
  `added_<release>.tar`, into `root_backups/` with its `SHA256SUMS`.
- **GB freed**: none.
- **Relay hours**: 3.1 for today's tars (*derived*; the first mirror with two streams and a read-back, measured note
  §3, takes longer); per release cut of `r26.11.1`'s size about 0.6 (20.9 GB, *derived*).
- **What a user loses**: nothing.
- **Way back**: from `/stor/work` as today (local). `tower`'s copy is used only if `/stor/work`'s is lost: relayed
  back, 3.1 hours for the whole root at 9.70 MB/s if the reverse is as fast (not measured).
- **What changes**: no code, no manifest, no profile. RUNBOOK "Mirror to tower" (this card). Each backup step (the
  release train's `root_added.tsv` and tar, card C-46; the stage drivers' backups) is followed by one
  `relay_to_tower.sh` run with the tar's line; a later card could have the train write that list line.

### Tier B: A, plus source archives, staged inputs and the deferred-data lists

- **What moves**: copies first, of the JGI archives, `provisional/`'s stage-input tars, the staged inputs on `/tmp`
  (one tar per stage, made on the pod), and the left-for-later lists. Then, **for each source whose stages are built
  and verified and whose remainder is recorded in `CATALOG_GOAL.md` §4**, the PI may remove its copy on `/stor`
  (agents never delete under `/stor`). Read each such file back once with `--reverify` before.
- **GB freed**: up to 1,137.9 GB on `/stor/work` (the JGI archives), where the footprint draft says under 50 GB; **none
  of their stages is built today** (`V2`, `GV1` staged; `J1` and the IUIG set not staged; STATUS, CATALOG_GOAL §3), so
  nothing can go yet. Up to 163 GB of `/stor/scratch/.../provisional/` (not split here: which of it is still read is a
  question for the card that frees it). On `/tmp`, 38 GB of staged inputs once their stages are ingested.
- **Relay hours**: JGI archives 32.6, `provisional/` 4.7, staged inputs 1.1; 38.4 hours in all (*derived*), one
  stream, laptop awake at home. The 524.8 GB IUIG archive alone is 15 hours in one piece (section 2).
- **What a user loses**: nothing a search or a release reads. A builder loses quick access to an archive whose `/stor`
  copy was removed: picking up a remainder (for example IUIG's 120,654 strain genomes, CATALOG_GOAL §4.1) first needs
  the archive back on the pod.
- **Way back**: relay `tower` to pod, at least 15 hours for the IUIG archive at 9.70 MB/s (rate of that direction not
  measured), or fetch again from the provider (the JGI sets need the PI's JGI session, card C-48).
- **What changes**: no code, no manifest (archives are in none). `sources.yaml` would record the second copy's place
  for each archive (a field beside its sha256; who edits `sources.yaml` and in which card is the PI's call), and
  `CATALOG_GOAL.md` §4's "how to pick it up" names `tower` where `/stor`'s copy has gone. A runbook step for the way
  back (relay to `/tmp`, check sha256 against `sources.yaml`).
- **Not `tower`'s only use for this**: `/stor/scratch` was agreed for 1.5 to 2 TB and holds 164 GB (STATUS open item
  19); moving the JGI archives from `/stor/work` to `/stor/scratch` would end the overage against the footprint with
  no relay and no loss of access. That is the PI's to weigh against the pod owner's wishes.

### Tier C: A and B, plus the large files of releases older than the newest N

- **What moves**: for every release but the newest N, its `release=` files (clusters, tokens, postings, dims, its
  release ledger) and its search pack leave the pod. Its manifest, catalog, `release_def` and the stage files stay
  (every stage is in the newest release, DEC-52 (3)). The `release=` files are **already on `tower` and on
  `/stor/work` inside the durable tars** once tier A is done (`l2.tar`, `l3.tar`, `added_r26_11_1.tar`), so "moving"
  them is removing them from `/tmp`; packs are not in the tars (they live outside the root) and are either relayed
  (33.55 GB for `rN.3` and `rM.2`, 1.0 hour) or simply built again when wanted.
- **GB freed** (on `/tmp`, with N = 1, keeping `r26.11.1`): the old packs 33.55 GB, the old clustering files 27.18
  GB, and the old tokens and postings inside `l3.tar` (at most 15.69 GB, which also holds the stages' tapes): 60.7 to
  76.4 GB of `/tmp`'s 613 GB free (*derived*). Later, for each release that ages out: about 0.35 GB per Gbp of release
  files plus its pack, 0.49 GB per Gbp for genome stages and 1.01 for unbinned contigs (*derived*, section 3.4, 3.5).
  At 100 Gbp a release that is about 84 GB (*derived*: 100 × (0.35 + 0.49)).
- **Relay hours**: none for the root's files (in the tars); 1.0 for the two old packs if they are kept rather than
  rebuilt.
- **What a user loses**: `--release rM.2` (and every older one): `search`, `nbr`, `find`, `export`, `releases`
  listing it as usable, `release verify` of it, and the worked examples on it. `latest` is unaffected.
- **Way back**: extract the release's files from the tars on `/stor/work` (local; a tar holds a whole layer, so
  `l2.tar`'s 27 GB is read for one release's clusters; a `tar -x` of named members), check them with `release verify`,
  build its pack again (498.86 CPU-minutes for a release of `r26.11.1`'s size, *measured*; about 2 hours at 12 jobs).
  From `tower` instead: the tars relayed back first (`l2.tar` and `l3.tar`, 42.87 GB: 1.2 hours, *derived*, if the
  reverse is as fast as 9.70 MB/s).
- **What changes**:
  1. A rule: RUNBOOK §9 says "never delete any file of a data root: a root is append-only", and DEC-52 (4) says
     "earlier releases stay and still verify". Removing an old release's files needs the PI's decision row. AGENTS
     §3.10 (a stage never rewrites an earlier stage's file) is not touched: no file is rewritten.
  2. A record that a release is archived, append-only in the root: a new `ArtifactSpec` in
     `src/biobox/registry/schemas.py` (AGENTS §3.7), for example `archived/release=<R>.json`, listing each removed file
     by path and sha256 and the tar that holds it (by the tar's sha256, with both copies' places).
  3. `catalog.verify_manifest` and `cli._verify`: an archived release checks the files it kept and reports the
     archived ones as archived, with exit 0, never as `missing:`; a file that is neither present nor in the record
     still fails.
  4. `release_train.sh` (`verify_rel`, the before and after loops) and `train_check.py`: verify archived releases by
     3, so a cut is not refused.
  5. `catalog.latest_release` and `manifest_releases`, `cli` `releases`: mark archived releases; `latest` never
     resolves to one.
  6. `l2/cluster.reuse_sources`: an archived release is not a source (else a cut fails after its plan).
  7. `catalog.open_catalog` and `runner.prepare` (through `pack.find`): refuse an archived release with the way back
     in the message, instead of a DuckDB I/O error.
  8. Manifests do not change. Paths stay relative to the root, and the archived copy is the same relative path inside
     a tar, so the manifest stays the one list of what the release is (AGENTS §3.9: readers use it, never a glob);
     what is new is a record of where a file is now. A release that lists files in two places is avoided. No profile
     key is needed, since the pod never reads `tower`; if one were wanted, `profile.Profile` must gain the field first
     (`profile.load` refuses unknown keys).
  9. Card C-46's release train would also need the number N, and DEC-52's open point ("how many earlier releases keep
     their large files and packs on the pod's local disk") answered.

### What is not proposed: tranches of built stages on `tower`

The card's first idea ("tranches can move there") does not fit the code or the decisions as they stand. A built stage
is in every later release (DEC-52 (3)), and every later stage's ingest must see its index (DEC-26, section 3.3); a stage
whose files left the pod would leave the release line and, worse, be stored again by the next stage without an error.
For it to work: a "parked" state of a stage, recorded in the root; DEC-52 changed so a release may leave a parked stage
out; ingest and gene calling reading a parked stage's `seq_locator` and `contig_index`, kept on the pod (small for
genome stages, section 3.3 (d); larger for unbinned contigs, which have 255 thousand contigs per Gbp at 1.5 kb,
*measured*, SPIRE pilot), while its shards and the rest go; readers refusing it. That is item 19's question (where a
corpus of several TB lives), and SPIRE at 5 kb would be about 3.2 TB at depth `called` (*derived*: 3.76 Tbp, card
C-42, × 0.8515 GB per Gbp, measured at 1.5 kb), 3.8 days of relay each way (*derived*). `tower` can keep such a
tranche cold; it cannot serve it.

## 5. Recommendation

1. **Tier A, kept up.** It is running; once its log says `DONE`, mirror each new backup tar the same day with
   `relay_to_tower.sh` (RUNBOOK "Mirror to tower").
2. **Tier B's copies, in this order**: the left-for-later lists and the staged inputs that are on `/tmp` only (small;
   a reboot loses them; `/stor` refuses files under 64 MB); then the JGI archives, the largest first only if the
   laptop can stay awake for 15 hours (else the others first), one stream. Removing any `/stor` copy waits for the PI,
   for the source's stages to be built and verified, and for a `--reverify` of the file on `tower`.
3. **Tier C not now.** It frees 61 to 76 GB of 613 GB free today (*derived*) for nine changes and a change of rule.
   Decide DEC-52's open point first; when `/tmp`'s free space falls under about two releases' worth, write the card.
   Until then the cheaper step is to drop old packs (a cache, rebuilt byte for byte), which needs only the PI's word
   that they may go (they are not in the root).

What I would not do:

- Put on `tower` alone any file a live release lists, or any file of a built stage (section 3.3: the next build goes
  silently wrong).
- Have any build or reader on the pod read from `tower` (no route; AGENTS §7.1 has outputs made only on the pod).
- Run two streams (they share the link: 8.3 MB/s together against 9.70 for one, *derived* and *measured*).
- Use the share `research` (its cache disk has 242 GB free, *measured*).
- Let `tower` hold the only copy of anything before the PI has decided on a copy of `tower` itself.
- Make a key for `tower` to pull from the pod while the laptop route suffices (the PI declined it; section 7 has the
  line in case).

## 6. Open questions

1. The rate from `tower` back to the pod (home upload) is not measured; every "way back" above assumes 9.70 MB/s. One
   small file back, with the PI's word, would settle it.
2. Which hop limits 9.70 MB/s is not measured (the card's own open point).
3. `tower`'s parity (single or dual) and its parity-check schedule are not read; they decide how much `--reverify` is
   worth.
4. Files larger than one laptop session (the 524.8 GB IUIG archive: 15 hours) need either a session that long or a
   resumable stream (byte offset plus a read-back of the first part to hash it). A card if tier B's copies are wanted.
5. Whether the rate of the repository script (hash on arrival, one stream) is the 9.70 MB/s expected (*estimate*):
   measure at its first real use and record it in the measurement note.

## 7. For the PI to set up or decide (not done by the agent)

1. **Retention on `tower`.** Proposed: keep every tar and archive, never delete or overwrite, one `SHA256SUMS` per
   directory (the script keeps it); look at free space at each mirror (18 TB free today; a release cut of `r26.11.1`'s
   size adds about 21 GB, *derived*).
2. **`/stor`'s copy once `tower` holds a verified one.** Agents never delete under `/stor`. Proposed rule: a source
   archive's `/stor/work` copy may go when (a) its stages are built and verified, (b) its remainder is in
   `CATALOG_GOAL.md` §4, (c) `relay_to_tower.sh --reverify` has passed on it, and (d) a second copy of `tower` exists
   or the PI accepts one copy. The root's tars stay on `/stor/work` in every case: they are the restore.
3. **A restricted pull key, if wanted later** (the PI said "no tower pull key" for now). The PI would make two keys on
   `tower` (one to pass through `wilkcomp03`, one to read on `wilkcomp02`) and add two lines to the account's
   `~/.ssh/authorized_keys` (one file serves both hosts if the home directory is shared on `/stor`, not checked). The
   lines, with the keys' public halves and two addresses filled in:

   ```text
   restrict,port-forwarding,permitopen="wilkcomp02.ccbb.utexas.edu:22",command="/bin/false",from="<tower's public address>" ssh-ed25519 <jump key> tower-jump-bio-box
   restrict,command="/usr/bin/rrsync -ro /stor/work/Wilke/if2453/all-the-data/root_backups/",from="<wilkcomp03's address as wilkcomp02 sees it>" ssh-ed25519 <pull key> tower-pull-bio-box
   ```

   `restrict` turns off forwarding, the agent, X11 and a terminal; the first line then allows only a tunnel to the
   pod's port 22 and no command; the second allows only `rsync` reads under the backup directory (`rrsync -ro`, which
   ships with `rsync` 3.2 and later; whether the pod has `/usr/bin/rrsync` is not checked). For the source archives a
   third line of the same form names `/stor/work/Wilke/if2453/all-the-data/jgi_20261118/`. I would give no key to the
   data root on `/tmp`: the tars hold it. **What risk remains**: anyone who takes the pull key from `tower` (a home
   machine whose `root` runs other services) can read everything under the named directory, including the JGI data
   whose terms say "not redistributed"; if the home address is not stable, `from=` on the first line must go, and the
   jump key then works from anywhere (the second line still limits it to reads); every pull is a process of our
   account on the pod outside `run-step` (I/O, not CPU, but AGENTS §7 names the pod's guests' rules, and the pod
   owner should know); UT's rules on unattended keys are the PI's to check.
4. **A second copy of `tower`.** Its parity protects against a failed disk, not against fire, theft, a power event or
   a mistaken command at home. While `/stor` keeps the first copy of everything (tier A), none is needed; before any
   `/stor` copy is removed (tier B), one is: an offsite or cloud copy of `bio-box/`, or TACC's Corral storage if card
   C-37a's request is granted (5,000 GB asked).
5. **DEC-52's open point** (how many earlier releases keep their large files and packs on the pod), which decides
   whether and when tier C is written as a card.
6. **Item 19 itself**: whether built tranches should ever leave the pod (section 4, "What is not proposed").
