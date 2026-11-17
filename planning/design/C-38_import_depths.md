# C-38 design: importing at a stated depth, and annotating on the hit (DEC-49)

Card C-38 (design half, tier H). DEC-49 keeps ideas 1 to 4 of `planning/sources_metagenomic_survey.md` §4.2 and
keeps the nucleotides. This note designs those four items, records in one paragraph what keeping nucleotides costs,
and lays out ideas 6 to 9 without recommending anything that depends on them. It decides nothing; the PI decides
from §9. It reads the code on `plan/task-cards` at `7f28896` (`registry/schemas.py`, `stage.py`, `catalog.py`,
`l1/domains.py`, `l2/cluster.py`, `l3/tokens.py`, `l3/postings.py`, `ledger.py`, `l0seq/inputs.py`,
`search/pack.py`, `search/runner.py`, `search/nbr.py`, `search/export.py`, `cli.py`) and the measurement notes
named beside each number. No measurement was made for this note.

Labels (AGENTS §6): *measured* (the note that measured it), *derived* (formula and inputs), *estimate* (its basis),
*published* (the source). Every number of the survey's §4.1 is quoted with the survey's own label. §1 lists every
number this note uses and the card that will replace it.

Words used throughout, each in one sense:
- **stage**: a set of source objects built into the data root by `build-stage` (`STAGES.md`); its files sit under
  `stage=<S>` partitions.
- **depth**: how far a stage was built: *stored* (sequences ingested), *called* (genes, proteins, intergenic regions,
  tapes, CRISPR arrays of the default producer), *annotated for a library* (every distinct protein of the stage has
  measured or inferred rows of that library, cards B-47, B-47b). §2 says why *clustered* is a property of a release
  and not of a stage.
- **unbinned contigs**: a metagenome assembly's contigs, whether or not a bin holds them (SPIRE, MGnify, INSDC primary
  metagenome assemblies; the crosswalk of card C-43).
- **round**: one scan of a library over a chosen set of proteins of stages that are not annotated for it (§3).
- **domain state** of a protein for a library, in a release: *measured* (it has measured rows), *none* (it was
  scanned and has no row), *inferred* (card B-47c's rows), *not scanned*.

## 0. The problem in ten lines

1. Today a stage is built all or nothing: `build-stage` runs every step it is given, and a stage with a
   `--library` is scanned in full (`stage.steps`).
2. A depth already exists in the code, unnamed: a stage built without `--library` has no `l1_domains` files, and
   `catalog.partial_libraries`, `release verify` and `bio-box releases` report it as a library missing on some
   stages (card C-19).
3. For unbinned contigs a full Pfam scan costs about 20 times gene calling (survey §4.1, *derived*): 413,000
   CPU-hours for SPIRE's 16 Tbp, 1.2 years of the pod at 40 threads.
4. Such contigs also carry about 2.2 million genes per Gbp, against about 0.93 million distinct proteins per Gbp in
   our genome stages (§1), so every per-protein cost (clustering, the pack, tokens) grows faster than bases.
5. A reader cannot tell a protein scanned with no domain from a protein never scanned: both print `-` in
   `pfam_architecture` (`search/export.py`) and both have no `l1_domains` row.
6. Annotating only around search hits needs a home for rows that belong to no single stage and arrive between
   releases; AGENTS §3.10 forbids rewriting a stage's file, and `catalog._scope` (card B-33) knows only stage,
   release and version partitions.
7. Rows made around hits are concentrated where someone already looked, so they would bias any family-level count
   that treats them as a sample of the corpus.
8. DEC-24 makes every release cluster all its proteins from scratch with our cascade; the cascade's measured cost
   grows as n^1.244 and its memory reached 178 GB at 45.8 million proteins (§1). A release with even the 1,000-sample
   pilot of card C-42 would not fit the pod.
9. A filter at ingest changes what a stage stores, which touches invariant 2 of `00-spec.md` §3 and the ledger,
   whose rows count records and not bases.
10. The design must keep a release a function of its stages (AGENTS §4.3, card B-25) and every rule of AGENTS §3
    and §4, or name the rule it would change as a question for the PI.

## 1. The numbers this note uses

| # | Quantity | Value | Label and source | Replaced by |
|---|---|---|---|---|
| N1 | Files at depth stored, per Gbp stored | 0.309 GB | *derived*: H1's depth-0 files (sequence shards with `.fai` and `.gzi`, `contig_index`, `contig_label`, `contig_occurrence`, `genome_input`, `seq_locator`, `l0_manifest`, ingest ledger), 561,924,885 bytes, over 1,819,831,338 bases; sizes *measured* in `2026-11-15_stage_H1/root_added.tsv` | C-42 (4) |
| N2 | Files added by depth called, per Gbp | 0.405 GB | *derived*: H1's `protein_seq`, gene-caller `feature`, `feature_protein`, `protein_occurrence`, tapes, IGRs, default CRISPR producer and their ledgers, 736,912,325 bytes, same source | C-42 (4) |
| N3 | Files at depth called in all, per Gbp | 0.714 GB | *derived*: N1 + N2; the survey's 0.71 | C-42 (4) |
| N4 | Pfam files added at depth annotated, per Gbp | 0.018 GB (phages) | *derived*: H1's `l1_domains` and its ledger, 32,271,623 bytes, same source; 24.1% of H1's proteins have a domain (*measured*, C-24) | C-42 does not measure it (no Pfam) |
| N5 | Sequence shard bytes per base | 0.403 (CAMI) against 0.310 (phages); `.fai` 85 bytes per contig, 34% of the shard | *measured*, `2026-11-09_design_stress_test.md` §6c, CAMI II marine sample 0 | C-42 (4) |
| N6 | Genes per Gbp of unbinned contigs | 2.16 million | *derived*: 1,963,951 genes / 0.9095 Gbp, CAMI II marine sample 0 (*measured*, same §6c) | C-42 (4) |
| N7 | Proteins per Gbp, SPIRE | 2.19 million | *derived*: 35 billion proteins / 16 Tbp (*published*, survey §2 #10) | C-42 (4) |
| N8 | Distinct proteins per Gbp, genome stages | 0.93 million | *derived*, DEC-29 | |
| N9 | Distinct proteins as a share of genes, CAMI | 98.8%; genes running off a contig end 83.5%; 1.33 genes per contig; 95% of slots have their window cut by a contig end | *measured*, `2026-11-09_design_stress_test.md` §6c | C-42 (4) |
| N10 | CAMI contig lengths | 91.6% of contigs under 1 kb; in a 100,000-contig sample, contigs under 1 kb hold 52.4% of bases, 1 to 10 kb 27.5%, 10 kb or more 20.0% | *measured*, §6c and `2026-11-10_genecall_short_contigs.md` (inputs). CAMI's file is its gold-standard assembly of simulated reads (`anonymous_gsa.fasta.gz`), not a MEGAHIT assembly like SPIRE's | C-42 (3) |
| N11 | Ingest | 4.11 CPU-minutes per 10^9 bytes of FASTA not held, plus 0.07 | *measured* fit, `2026-11-16_ingest_estimate.md` (the survey's 4.5 is the rate before card B-52) | C-42 (4) |
| N12 | Gene calling, meta mode | 60.4 CPU-minutes per Gbp of linear contigs | *derived*: 2.9450 s per Mbp (*measured*, `2026-11-15_genecall_batches.md`, PL1's linear sample) × 1.23 step overhead × 1,000 / 60 | C-42 (4) |
| N13 | CRISPR arrays, default producer | 3.0 CPU-minutes per Gbp | *derived*: 182.0 CPU-s per Gbp (*measured*, card E-12, `stage.CRISPR_CPU_S_PER_GBP`) | C-42 (4) |
| N14 | Pfam scan | 7.18 CPU-minutes per million residues; about 1,550 per Gbp | *measured* on H1 (C-24); per Gbp *derived* in survey §4.1 (2,817 / 1.82) | C-37 (3) for another host |
| N15 | Fixed cost of one Pfam scan process | 1.533 CPU-minutes | *measured*, H1 shard 13 (DEC-43) | |
| N16 | Mean protein length, P2 | 323.6 residues | *measured*, B-47 note §6 | C-42 (4) for unbinned contigs |
| N17 | Our clustering cascade, release | `284.10 × (n / 2,053,075)^1.244` CPU-minutes; 14,087 measured at 45.8 million proteins | *derived* power law (`stage.cluster_cpu_min`); *measured* `rM.2` (card B-49) | C-38c (§10) |
| N18 | Cascade memory | 178 GB peak at 45.8 million proteins (`rM.2`) | *measured*, quoted in card C-37 | C-38c |
| N19 | Share of the cascade's CPU in linclust | 7.4% of the NR90 clustering at 10.6 million proteins; about 6% of both levels on a 100,000-protein sample | *derived*: 73.0 of 992 CPU-minutes, `rN.3`'s run records (B-47 note §11); B-37 note §2 (kmermatcher 4.7, linclust `align` 9.5 and 14.6, `rescorediagonal` 1.2 CPU-s of 489.8) | C-38c |
| N20 | Linclust memory | 336 bytes per sequence | *estimate*: `planning/sources_and_sizing.md` §4, the MMseqs2 wiki formula "recalled, verify" | C-38c |
| N21 | One HMM search | 8 to 11 CPU-minutes over 10.6 million proteins | *derived*, survey §4.1 from card W-21 | C-42 (6) |
| N22 | Release files, search pack | at least 17 GB of release files and a 25.8 GB pack for `rM.2`'s 52.1 Gbp: at least 0.32 and 0.495 GB per Gbp | *measured* sizes (DEC-52; card W-35), bases *measured* (card C-33); per Gbp *derived* | C-42 for unbinned contigs |
| N23 | Pod disk | local disk 1.1 TB, data root 80 GB; `/stor/scratch` 22 TB shared | *measured*, survey §4.1 (`df`, 2026-11-16) | |
| N24 | C-42's pilot | about 160 Gbp, about 350 million proteins | *estimate* (survey: 16 Tbp / 99,146 samples × 1,000); proteins *derived* at N7 | C-42 (3), (4) |

C-37 replaces every TACC figure (all assume a Lonestar6 core does what a pod thread does). C-42 replaces every figure
labelled "C-42" above; until then the metagenomic numbers rest on CAMI's simulated sample and on SPIRE's published
totals. C-38c (§10) replaces N17 to N20 above 45.8 million proteins, where nothing has been measured.

## 2. Depth as a recorded property (item 1)

### 2.1 What the depths are, and why "clustered" is not a stage's

| Depth | Steps of `build-stage` (`stage.steps`) | What a release gets from the stage |
|---|---|---|
| stored | `inputs`, `ingest` | sequences, `contig_index`, labels; no tape |
| called | stored, then `genecall`, `igr`, `tapes`, the default CRISPR producer | proteins, tapes, IGRs, arrays; NR90 and NR50 tokens from the release's clustering |
| annotated for library L | called, then `domains-<L>` (full, copied, or by representatives under DEC-35); for Pfam also the filtered CRISPR producer (card E-15, which needs the stage's Pfam rows); ncRNA (`features` slot) is a producer of its own and counts the same way | family tokens and postings for L on the stage's slots |

DEC-49 names four depths: stored, called, clustered, annotated. Under DEC-24 clustering belongs to a release: every
release clusters all of its proteins from scratch, and card B-51 refuses tokens for a release with protein slots and
no clustering. So every called stage in a release is clustered by that release, and no stage file records a
clustering. Stage-local clustering exists only as card B-47a's `represent` step, which serves annotation by
representatives and is part of "annotated" here. **The note treats "clustered" as a property of a release**, and
asks the PI to confirm it (Q1).

ncRNA is annotation too, and dearer than Pfam per base: 358 CPU-hours per 10 Gbp with `prok-2` (*derived*, DEC-30,
`stage.NCRNA_SET_COST`), about 2,150 CPU-minutes per Gbp against Pfam's 1,550. At depth called it is off.

### 2.2 Where depth is recorded: three options

**(a) Declared, in a stage table.** A new `TableSpec` `stage_def` (partition `stage`, `small_ok`), written once by
the `inputs` step, with the depth the stage is built to.
- *Strains:* AGENTS §3.10. A stage later annotated in full (the PI may annotate a pilot after all) cannot change its
  declared depth: the file is written once. A second file would need a version key and a rule for which one wins.
- *Old stages:* S1 to `PL1` have no such file; adding one later is the case card C-34 had to special-case
  (`STAGE_ADDED_TABLES`).

**(b) Derived from the files a release lists.** `catalog.stage_depths(manifest)` (new, beside `partial_libraries`)
reads the manifest's file list: a stage with `protein_seq` and tapes is called; with `l1_domains` of library L it is
annotated for L (and by representatives when card B-47c's `inferred_from` is set on any of its rows, which the
function reads from the manifest's listed files only when asked); with an ncRNA producer's files it has ncRNA.
- *Rules:* none strained. The manifest is already the record (AGENTS §3.9); a stage annotated later gains files
  under its `stage=` partition, which card C-19 already lets a later release take and an earlier one ignore.
- *Old manifests:* the function works on every manifest written so far; `MANIFEST_VERSION` does not change.
- *What it cannot say:* whether a stage without Pfam was built at depth called on purpose or awaits its scan. A
  release does not need to know; the build does (below).

**(c) Derived, and also written into the manifest** as a field `stage_depths`. The same facts as (b), easier to read
by eye, at the cost of a manifest version and a field that must agree with the file list.

**Recommendation: (b).** The build's intent is an argument of `build-stage` (`--depth stored|called|annotated`,
§2.3), and what a release holds is derived from its manifest. Depth can only grow, by adding files, never by
rewriting one.

### 2.3 What each part of the system does with depths

- **`build-stage`.** `--depth called` builds the list of `stage.steps` without the domain steps, the filtered CRISPR
  producer and the ncRNA producer, and refuses `--library`, `--ncrna-cm` with it; `--depth annotated` requires a
  `--library`; `--depth stored` stops after ingest. Default: `annotated` when a library is given, else `called`,
  which is what the code does today. Card C-46's driver judges a stage complete against the step list of the depth
  it was built to, so the function that makes the list takes the depth (`stage.steps(prof, b)` already does, through
  `Build`).
- **A release.** It takes stages at depth called or more. A stage at depth stored has no tape, and every release
  builder (tokens, postings, `dim`, the pack) reads tapes per stage; it waits until called (Q2). A release that mixes
  called and annotated stages is expected (DEC-49), and is built as today: clustering over all proteins (DEC-24),
  NR90 and NR50 tokens on every slot, family tokens only on slots of stages annotated for the library (`tokens.
  _domain_files` already reads a library's files stage by stage, "a library need not be scanned on every stage").
- **Family-level statistics in a mixed release.** A slot of a called stage carries no family token. Counted as "no
  family", it would read as evidence of absence and lower every family's frequency. Rule proposed: a family-level
  count states its universe, and the default universe is the slots of the stages annotated for that library;
  cluster-level counts use every slot. Invariant 11's null models are not built yet; the rule is for whoever builds
  them, and for `nbr summarize` now (it prints "not scanned", §3.4).
- **The search pack.** It takes the `unique` proteins of every stage of the release, so a search reaches called
  proteins with no change (`search/pack.py`: every `protein_seq` row). The pack's `igr` and `dark_orf` targets are,
  by invariant 16 of `00-spec.md` §3, built for the whole release. For unbinned contigs they are large: dark ORFs are
  every stop-to-stop ORF of 30 codons or more outside called genes plus all frames in the 3 kb at each contig end,
  and CAMI's median contig is 294 bp (*measured*, the source of N9), so every frame of most short contigs is in
that end zone. Whether a called
  stage's IGRs and dark ORFs enter the pack is a change to invariant 16, so it is the PI's (Q3); the recommendation
  keeps the invariant as written until C-42 measures the pack's bytes per Gbp of unbinned contigs.
- **`release verify`.** Unchanged in what it checks. Its warning for a library on some stages only (card C-19)
  becomes one information line, `depths: S1,P1..PL1 annotated pfam-38.2; SP1 called`, since in a mixed release it is
  no longer a sign of a missed scan.
- **The manifest.** Unchanged (option b).
- **`bio-box releases`.** A column `depths`, in the same form, from `stage_depths`; with rounds (§3) a column
  `rounds` with their number.

### 2.4 What it does to the rules, tests and cost

- *Rules:* AGENTS §3.7 (no new table under option b), §3.9 (readers keep the manifest's lists), §3.10 (depth grows
  by new files only) hold. Invariant 1 of `00-spec.md` §3 (gated stages) is untouched: a stage at depth called passes
  its own gate.
- *Tests* (AGENTS §5; each fails first): (1) on the fixture, a release over one called and one annotated stage builds,
  verifies, and its family tokens are those of the annotated stage alone, checked against a brute-force count over
  the domain files written separately; (2) `stage_depths` of every manifest in `tests/` gives what the stage list and
  libraries say, and of an old manifest (a copy of `rM.2`'s file list) "annotated pfam-38.2" for all seven stages;
  (3) `build-stage --depth called --library ...` and `--depth stored` in a release step are refused, exit 2, nothing
  run; (4) the step list at each depth is a prefix-and-subset of today's full list, so a stage built called and later
  given `--library` reaches byte-identical files to a stage built annotated at once (manifest hash compared).
- *Cost:* none at build time; depth called saves the scan, at 1,550 CPU-minutes per Gbp (N14).

## 3. Annotation on the hit (item 2)

### 3.1 The path

1. `bio-box search hmm` (or `search run`) over release R. The pack's `unique` target holds the proteins of every
   stage, called ones included.
2. Its hits, and the slots within N genes of each on the release's tapes (`nbr`, the window of the run).
3. Of those proteins, the ones whose domain state in R for library L is *not scanned*: a **round** scans them with L.
4. The rows are kept in the data root, so the next search near the same proteins finds them scanned, and the next
   release lists them.
5. The run's own neighbour tables (`neighbors.tsv`, `loci/`) are written with R's domains plus the round's.

Steps 3 to 5 are new. The rest is today's `search hmm` (card W-07).

### 3.2 Where a round's rows live: four options

**(a) In the stages' own partitions.** A new table `round_domains` partitioned (`library`, `stage`, `round`): each
round writes, for each stage that holds a scanned protein, a file under that stage.
- *Scope:* no new rule. Rule (a) of card B-33 gives the files to every release that lists the stage, and card C-19's
  owner key (table, stage, library, round) makes an earlier release ignore them.
- *Strains:* the file count. A round's proteins come from every stage the hits fall in; SPIRE in stages of about
  10 Gbp would be 1,600 stages (*derived*: 16 Tbp / 10 Gbp), so one round could write 1,600 small files. A protein
  held by two stages must be assigned to one by a rule that depends on the release.
- *Lookups:* a key lookup reads one row group per round and stage file (AGENTS §4.4 fails as rounds accumulate).

**(b) One partition per round, and a scope rule for it.** Two tables, each partitioned (`round`, `library`):
`round_domains` (the columns of `l1_domains`), one row per included domain, and `round_coverage` (`protein_uid`,
`n_rows`), one row per protein the round covered, `n_rows` 0 for a protein with no domain. A small table
`release_round` (partition `release`, `small_ok`; columns `position`, `round`) lists the rounds a release takes,
written once by `define_release` beside `release_def`. New scope rule (e) in `catalog._scope`: a file with a
`round=<id>` partition belongs to every release whose `release_round` lists that round.
- *Lookups:* the release step `rounds` (after `cluster`, before `tokens`) writes the union of the release's rounds,
  one file per table sorted by `protein_uid`, under the alternative partition keys (`release`, `library`), so a lookup
  by protein reads one row group (AGENTS §4.4). This is card B-27's form: a lookup table re-sorted beside its source
  (`member_cluster` beside `cluster_member`). The two key sets differ in their first key (`round`, `release`), as
  `TableSpec.alt_partition_keys` requires.
- *Strains:* `_scope`, `_in_release`, `ReleaseFiles._layout`, `build_manifest` and `verify_manifest` learn the new
  kind; `define_release` writes a second small file. Nothing written earlier changes.

**(c) In the search run's directory only.** The round's rows go into the run's outputs, outside the data root (card
W-03), and are archived with it (DEC-33).
- *Strains:* none of AGENTS §3. But rows are not reused: the next search scans the same proteins again, and no
  release sees them. DEC-49 says "those rows kept"; this keeps them only as a record.

**(d) In the searched release's `release=` partition.** Refused by rule (b) of card B-33: R's manifest exists, and
a new file there is "not in manifest" for R.

**Recommendation: (b).** It is the one option that keeps one file per round, prunes, and adds no file to a published
release. The round's id is content-derived: `r-` and the first 16 hex digits of the sha256 over the library's name and
version, the HMM file's sha256, cutoff, Z, and the sorted protein uids requested. The same request is the same round,
and asking again is a no-op.

### 3.3 What a round scans, and how rounds accumulate without rescanning

- **Requested:** the protein uids of the slots in the window whose domain state in the searched release R is *not
  scanned* (§3.4), sorted. R is fixed, so the list is a function of the run and the window. The request is a plain
  list, so the same step serves a search's neighbours and any other list (§4.4's representatives).
- **Left out of the scan, but copied:** proteins of the request already measured on disk with the same HMM sha256,
  cutoff, Z and tool versions (card B-47b's `COMPARED`): in a stage annotated for L after R was cut, or in an earlier
  round, read from that stage's or round's plan record.
- **Copied, not rescanned:** as card B-47b does for a stage, the round's file holds the rows of every requested
  protein: scanned now, or copied from the stage or round that measured it. At a fixed Z a protein's rows do not
  depend on the other proteins (MR-12), so the round's files are a function of its request alone, byte-identical
  whichever stages and rounds were made before (AGENTS §4.3). The split between scanned and copied
  goes in the round's plan and shard records only (a `domain_scan`-like artifact under `round=<id>`, `in_manifest`
  False). The price is duplicate rows across rounds that overlap: about 54 bytes per row (`domains.EST_ROW_BYTES`),
  small beside N2.
- **One source per protein:** a reader of a release takes a protein's rows from an annotated stage when one holds it,
  else from any round (their rows are equal), as the B-47 note §1.1 proposes for inferred rows.
- **Shards and budget:** the round plans, scans in shards and merges exactly as `domains.plan`, `scan_shard` and
  `merge` do, each a `run-step` under 300 CPU-minutes (AGENTS §7.2). It runs on the pod only as a `run-step` of the
  build account (AGENTS §7.1), never as part of a user's search process: a search user may be unable to write the data
  root (card W-03). `search hmm --annotate` runs it when the caller can write the root, else prints the command.
- **Ledger:** `ledger/round=<id>/step=domains-<L>/`, a third key set of the ledger (beside (`stage`, `step`) and
  (`release`, `step`)), with the rows of a stage scan (`n_in` requested, kept with a domain, `no_hit`) plus
  `settled_elsewhere` for the proteins left out. It sits outside `ledger.CHAIN`'s `domains-*` check, whose `n_in` is a
  stage's gene count.

### 3.4 How a reader tells "scanned, no domain" from "never scanned"

Today it cannot (§0 line 5). Proposed: a catalog view `domain_state(protein_uid, library)` with four values.

| State | Rule |
|---|---|
| measured | a measured row in a listed `l1_domains` file or in the release's `round_domains` |
| inferred | only inferred rows (card B-47c's `inferred_from` set) |
| none | the protein is in a stage the release lists as annotated for L (§2.2: every protein of that stage was scanned or inferred; card B-47b decides "measured" from records, never from rows, and the same holds here), or in the release's `round_coverage` with `n_rows` 0 |
| not scanned | none of the above |

The stage clause needs no per-protein row: a stage annotated for L covers every protein in its `protein_seq`.
Exports print the state where they now print `-`: `pfam_architecture` becomes `-` for *none* and `?` for *not
scanned* (`search hits`, `nbr`, `nbr summarize`, `nbr flanks`), and the README of `search hmm`'s output says so.

### 3.5 Rounds, family tokens and statistics

Round rows are concentrated around past hits. Counted in family tokens and postings they would make every family near
an earlier search look enriched near it. Inferred rows (DEC-35 (5)) cover whole stages and carry no such bias; round
rows do. **Proposed: round rows are kept out of family tokens and postings by default**, and shown in neighbour
tables, flanks, exports and `domains_of` with a column `source` (`stage`, `round`). Family statistics then stay over
the annotated stages (§2.3). This is the PI's (Q5).

### 3.6 With card B-47b's copies and card B-47's inferred rows

- A round never requests a protein of a stage of R annotated for L, whether its rows there are measured, copied
  (B-47b) or inferred (B-47c): its state is settled. A protein measured by a stage built after R is copied, not
  scanned; one only inferred there is scanned (inferred rows are never copied).
- Rounds are sources for a later stage's scan (card B-47b): a stage built later and annotated in full copies the
  round-measured proteins' rows instead of scanning them. `domains.sources` learns to read a round's plan record and
  `round_coverage` as it reads a stage's plan; the copied rows are equal (MR-12).
- Inferred rows are never copied into a round and never infer from round rows: inference stays a stage property
  (DEC-35 (1)).

### 3.7 A contig an annotated stage shares with a called stage

Under DEC-26 the stage that stored a contig first owns its genes. If a called stage stores contig c first (an
unbinned contig of a sample) and a later annotated stage's genome holds c (a MAG of GlobDB, a plasmid), c's genes are
the called stage's and get no Pfam: the annotated stage writes no feature rows for c, so its `protein_seq` and its
scan miss them. With rolling releases (DEC-52 (2)) stages are built in the order they come, so this will happen.
Proposed: an annotated stage's scan plan also takes the proteins of contigs its genomes hold and a stage at a lower
depth owns (read from that stage's `feature_protein`, as gene calling already reads other stages with
`--seq-stages`), copying those already measured. Its domain file then holds rows for proteins not in its own
`protein_seq`, which the B-47 note §2a called odd but which no reader forbids: tokens join by `protein_uid` (Q8).

### 3.8 Rules, tests and cost of the recommendation

- *Rules:* §3.7 (four new `TableSpec`s: `round_domains`, `round_coverage`, `release_round`, and the round's plan
  artifact; the ledger's third key set), §3.9 (readers use the release's merged files, listed in its manifest),
  §3.10 (rounds are new files; the release's merged copy is under its own `release=`), §4.1 (the merge streams sorted
  runs, as `domains.merge`), §4.3 (§3.3), §4.4 (the merged lookup). The new scope rule (e) changes no earlier
  release's scope: none lists a round.
- *Tests:* (1) oracle: a planted fixture (a hit, neighbours with and without domains, a neighbour in an annotated stage,
  a protein in two stages), the round's rows against an `hmmsearch` run of the image's binary over the requested
  proteins, pair by pair; (2) two rounds that overlap, made in either order, give byte-identical files, and the second
  scans only the proteins the first did not request (counts from the plan records); (3) `domain_state` on the fixture
  for all four states, against a brute-force reading of the files written separately; (4) a release taking the
  rounds verifies, an earlier release still verifies with its hash unchanged, a stray file under `round=` that no
  release lists is reported only where rule (e) says; (5) `profile_query` of a lookup by protein in a release with
  ten rounds reads one row group of the merged files; (6) family tokens and postings of a release are byte-identical
  with and without rounds (Q5 as proposed); (7) memory test on the 100× synthetic store for the round's plan and merge.
- *Cost (all derived; window N genes on each side):* a hit's window holds at most 2N + 1 proteins. At N = 10 and
  323.6 residues (N16) that is at most 6,796 residues, 0.049 CPU-minutes at 7.18 per million (N14), plus 1.533 per
  scan process (N15). 10,000 hits cost at most about 490 CPU-minutes, against 260 to 360 CPU-minutes for the search
  itself over C-42's pilot (N21 scaled to N24's 350 million proteins) and 450 to 580 CPU-hours over all of SPIRE
  (survey §4.1). Short contigs make it smaller: CAMI's contigs carry 1.33 genes each (N9). The rows are about 54 bytes
  each. C-42 (6) measures one round end to end.

## 4. Protein families without Pfam (item 3)

### 4.1 What exists, and what is new

Every release already has cluster tokens at NR90 and NR50 on every slot, and `nbr summarize` uses `nr50,nr90` when a
release has no library (`cli._search_hmm`). Neighbour statistics on cluster tokens therefore need no new table. What
idea 3 adds is a clustering that scales to unbinned contigs, and Pfam on the representatives of the families that
matter.

### 4.2 Why the cascade does not scale, under DEC-24

DEC-24: every release clusters all its proteins from scratch, no `clusterupdate`. DEC-52: a release is cut after every
completed stage. Together, every release pays a full clustering of every protein.

| Proteins | Cascade CPU-minutes (N17) | Cascade memory, linear in proteins from N18 | Linclust both levels, CPU-minutes | Linclust memory (N20) |
|---|---|---|---|---|
| 45.8 M (`rM.2`, measured) | 14,087 *measured* | 178 GB *measured* | about 640 *estimate*: 14 CPU-minutes per million, twice N19's NR90 rate of 6.9 (*derived*: 73.0 / 10.62), assuming the NR50 linclust costs what the NR90 one does | 15 GB *estimate* |
| 350 M (C-42's pilot added, N24) | about 170,000 *derived*, extrapolated 7.6 times past the largest measurement | about 1.36 TB *estimate* | about 4,900 *estimate* | about 118 GB *estimate* |
| 35 billion (SPIRE, *published*) | about 5.2e7 *derived*, extrapolated | about 136 TB *estimate* | about 490,000 *estimate* | about 11.8 TB *estimate* |

So the first release that takes the pilot cannot run our cascade on the pod (200 GB scopes, AGENTS §7.2) or on one
Lonestar6 node (256 GB, card C-37); a Stampede3 node of 4 TB might hold it (C-37 lists them). Linclust alone would
hold the pilot within a 200 GB scope by N20's estimate, and nothing we have holds all of SPIRE in one clustering. DEC-24
then stops the design at a question (Q10): a release of several billion proteins needs a clustering that is split by
proteins and merged, which is the incremental clustering DEC-24 rules out, or no such release.

### 4.3 Options for a release's clustering

**(a) The cascade, as now.** Exact to what we have; out of reach above a few hundred million proteins (§4.2).

**(b) Linclust alone, both levels, for every protein of a release.** The cascade without its prefilter steps: the
`kmermatcher`, `rescorediagonal`, `align` and `clust` calls of `l2/cluster.workflow`'s redundancy phase, then members
aligned to representatives as now. A new `cluster_params_hash`, so every `cluster_uid` changes once (cluster ids hash
the member set and the parameters, invariant 10; the lineage table of DEC-42 (5) maps old to new).
- *Lost against the cascade:* linclust joins a protein only to the longest sequence of a group sharing one of its 21
  selected k-mers in a reduced alphabet (B-37 note §1), and has no sensitive prefilter; at 50% identity it will leave
  members the cascade's three prefilter steps find (`-s 1`, `2`, `3`). Expected: more NR50 clusters, more singletons
  (the cascade already gives 89% NR90 singletons on CAMI, N9), lower counts per token and so less power in neighbour
  statistics. **None of this is measured** (B-47 note §11 item 4); card C-38c measures it.
- *Rules:* §4.3 holds (linclust is deterministic as the cascade is, `--shuffle 0`, same threads rule). DEC-24 holds.

**(c) The cascade over annotated stages, linclust over called stages.** Two methods in one release; a protein's
token level means two things. Rejected for that, and the annotated part also outgrows the cascade at GlobDB's size
(about 400 Gbp, survey §3, *estimate*).

**(d) Stage-local clustering joined at the release.** Each stage clustered once, the release merging stage
representatives. This is incremental: the release's clusters depend on how proteins were split into stages. DEC-24
rules it out.

**Recommendation:** measure (b) against (a) first (card C-38c, on `r26.11.1`'s proteins, read-only, in scratch), then
choose one method for the release line (Q9). If the loss is acceptable, (b) for every release from the first that
takes an unbinned stage; one method throughout, so that NR50 means one thing in every release.

### 4.4 Pfam on representatives of the families that matter

Under (b), Pfam can go to the representatives of NR50 clusters above a size, or near hits, as a round (§3): the request
is the list of representatives, the rows are measured rows of those proteins. A member is then shown "in a cluster
whose representative has Pfam X" by a query-time join (`member_cluster` → the representative's state), labelled so,
never written as rows and never in family tokens. Writing inferred rows from release clusters would make domains a
property of a release, which DEC-35 (1) chose against. Cost: one representative per cluster above the size; the size
is chosen from C-38c's and C-42's cluster-size distribution (Q14).

### 4.5 NR90, NR50 and rolling releases

Idea 3 adds no third level: the release's NR50 is the family. With rolling releases each release reclusters (DEC-24);
card B-49's reuse applies only when no protein changed, which a release adding an unbinned stage never meets. At the
pilot's size, (b) costs about 4,900 CPU-minutes a release (*estimate*, §4.2), at least 21 processes at 80% of the
300 CPU-minute budget when split as card B-41 splits the cascade's calls (*derived*: 4,900 / 240); (a) costs about 170,000 (*derived*, extrapolated). Cluster
uids churn between releases as members join (card 0-07's lineage table follows them).

- *Tests for (b) when chosen:* on the fixture, linclust-only clusters are byte-identical at 1 and 4 threads and with
  input order shuffled; every NR90 member's stored alignment satisfies 90% identity and 80% member coverage against
  its representative, checked by a separate brute-force aligner on a sample (Biopython's `PairwiseAligner`, or the
  image's `mmseqs align` used as an oracle only); a release clustered both ways on the fixture differs only in
  `cluster_params_hash`-bearing files.

## 5. The filter at ingest (item 4)

### 5.1 Options

**(a) A least contig length L, applied by ingest.** Ingest reads every record and hashes it already; a record under
L bases is not stored. The stage's inputs (`genomes.json` from the input script) are made by the same rule, and
ingest checks that the two agree: a contig of `genome_input` shorter than L, or a stored contig not in it, is
refused (one streamed comparison on `contig_uid`, as `held_bases` reads `genome_input` now).
- *Ledger:* ingest's rows are per source object, in records (`ledger.py`: "per L0 object: records in; kept =
  byte sequences stored new; merged = records whose bytes were already stored"). A row with reason `below:<L>` and
  `n_rejected` the records dropped. Bases are not in the ledger today; two nullable columns are appended,
  `n_bases_in` and `n_bases_rejected` (allowed by LY-09 before and after the first signed release, AGENTS §3.8), set
  on ingest's rows. L-01's sum per object stays in records.
- *Determinism:* the rule reads a record's own length only, so the result is the same in any input order or build
  order (AGENTS §4.3).

**(b) A least gene count, after gene calling.** The contig is stored and called, then left off the tape. It saves
tapes, tokens and postings but not storage or calling, and its proteins would still be in `protein_seq`, the pack and
the release's clustering unless removed again. Rejected: it saves least where cost is.

**(c) Store everything, call only above L.** The short contigs at depth stored. Cheaper than (a) to reverse, dearer in
storage: short contigs cost most per base (N5: `.fai` 85 bytes per contig, 34% of CAMI's shard).

**Recommendation: (a), by length.** A gene count needs genes; length is its proxy (CAMI: 1.33 genes per contig at a
median of 294 bp, N9). L is chosen from C-42's measured shares of bases over 1.5, 5 and 10 kb per habitat (card C-42
says the numbers come first; Q12). Every release states its filters and the counts that survived them (D-01's scope
statement, `00-spec.md` §1); the ledger carries them.

### 5.2 A contig another stage already stores

The length rule is applied before the check for contigs held by other stages, so the decision does not depend on
which stage was built first.
- Shorter than L, held by an earlier stage: the filtering stage drops it (counted under `below:<L>`, not as held);
  the earlier stage keeps owning it.
- Shorter than L, stored later by a stage without the filter (a plasmid stage): stored there as usual.
- L or longer, held by an earlier stage: DEC-26 as today; the filtering stage records the skip.
- Consequence: a dropped record leaves no `contig_occurrence` row, so "in how many samples does this contig occur"
  is undercounted below L.

### 5.3 What it changes, tests and cost

- *Invariant 2 of `00-spec.md` §3* says every source contig's sequence is stored. A filtered contig is not. If the
  PI reads D-01's "filters F" as placing it outside the corpus, the invariant holds as written; else the invariant
  changes. Q11.
- *Unit of a genome.* An unbinned sample is not a genome: one `genome_uid` per sample would let gene calling use
  single mode on a sample above 100 kb (DEC-13), and one per contig gives `genome_dim` a row per contig (as the CAMI
  pilot did). This is card C-42's input script to settle (meta mode for every unbinned contig, DEC-13 "fragments"),
  and counts "by genome" (card C-21c) must say which.
- *Tests:* (1) a FASTA with records on both sides of L, in two input orders: equal files, and the ledger's dropped
  records and bases equal a brute-force count written separately in the test; (2) a contig below L held by another
  stage is counted under `below:<L>`, not as held, built in both orders; (3) inputs that disagree with the rule are
  refused before anything is written; (4) a stage without the filter gives today's files and ledger byte for byte
  (fixture manifest hash unchanged, the new columns null); (5) memory test on the 100× synthetic store.
- *Cost:* none; it removes work. On CAMI's sample, L = 1 kb would drop 91.6% of contigs and 52.4% of bases (N10);
  SPIRE's MEGAHIT contigs will differ, which C-42 (3) measures.

## 6. Nucleotides are kept (item 5, dropped as a proposal)

Storage is planned at the full figure: 0.714 GB per Gbp at depth called (N3), likely more for short contigs (N5). For
the record only: the reference form (proteins, coordinates and arrays kept; sequences refetched from the source by
sample, contig and sha256) would have dropped the sequence shards and IGRs, 0.308 and 0.038 GB per Gbp (*derived*, H1:
560,348,485 and 70,023,748 bytes over 1.8198 Gbp), leaving 0.367 GB per Gbp: 5.9 TB instead of 11.4 TB for SPIRE's
16 Tbp (*derived*).

## 7. What this note lays out and does not recommend (ideas 6 to 9)

Nothing above depends on these. Each needs a measurement first.

**Where a corpus of several TB lives (idea 9, card item 6).** At N3, the pilot is about 114 GB at depth called, plus at
least 52 GB of release files and 79 GB of pack per release (*derived*, N22 × 160 Gbp; more per base for unbinned
contigs); all of SPIRE is about 11.4 TB at depth called. The pod's local disk is 1.1 TB (N23).
- Options: (i) the depth-called part of the root moved to `/stor/scratch` (shared 22 TB; DEC-16 approved about 1.5 to
  2 TB through S5, so this needs a new footprint request); (ii) TACC storage (`$WORK` 1 TB per user, `$SCRATCH`
  purged after ten days, card C-37), which cannot hold the root; (iii) the depth-called stages in a second data root
  with its own releases, searched beside the main one.
- Rules touched: AGENTS §7.3 (files under 64 MB on `/stor`), §7.5 (nothing under `/stor` overwritten), §3.9 (a
  release's files listed by one manifest, under one root).
- To measure: bytes per Gbp at depth called for unbinned contigs, and the pack's (C-42 (4)); NFS read rates for tape
  and sequence lookups (not yet carded).

**Stopping by habitat (idea 6, card item 7).** Take a habitat until new NR50 families per Gbp fall under a fraction of
its first value.
- Options: the fraction; per habitat level 1 or level 2 of SPIRE's ontology (C-43's `spire_habitat_level2.tsv`); the
  clustering that defines "new" (cascade or linclust, §4).
- To measure: the curve of new families per Gbp as samples are added, per habitat (C-42 (5)). It also depends on Q9.

**Importing from Logan on the hit (idea 7).** Search where Logan lives, pull the contigs of the runs that hit, and
build them as a small stage that records its query.
- Fits this design as a called stage (§2) whose neighbours get a round (§3). Open: what a stage records about the
  query that made it (`sources.yaml` entry per pull, the run that selected the runs).
- To measure: what Logan publishes per run, the cost of a protein search over its metagenome contigs, a trial of at
  most 1,000 runs (C-44).

**The machine for each depth (idea 8).** Ingest, calling and scans are many independent pieces; a clustering of
hundreds of millions of proteins needs one large memory.
- Options: the pod for all of depth called; TACC's ordinary nodes for scans and rounds; Stampede3's 4 TB nodes for a
  cascade above the pod's 200 GB (§4.2).
- Rules touched: AGENTS §7.1 says outputs are produced only through `run-step` in the image on the pod; results made
  elsewhere entering the data root is a decision C-37 (4) raises.
- To measure: a Lonestar6 core against a pod thread, and byte identity across hosts (C-37 (3)); clustering memory
  above 45.8 million proteins (C-38c).

## 8. Recommendation

1. **Depth is derived from what a release lists, and chosen at build time.** `build-stage --depth stored|called|
   annotated`; `catalog.stage_depths(manifest)` reports it; `release verify` and `bio-box releases` print it; no new
   stage table and no manifest version. "Clustered" is a property of the release (DEC-24). A release takes stages at
   depth called or more. A mixed release is normal; family-level counts state their universe, by default the stages
   annotated for that library.
2. **Annotation on the hit is a round.** One partition per round, `round=<content id>/library=<L>`, holding
   `round_domains` and `round_coverage`; a release lists the rounds it takes in `release_round` (a new scope rule
   beside card B-33's) and writes one merged, sorted copy under its own `release=` for lookups. A round scans only
   proteins whose state is not settled, copies rows measured by earlier rounds, and its files depend on its request
   alone. Rounds run as `run-step` steps of the build account. Round rows are shown in neighbour tables and exports and
   kept out of family tokens and postings.
3. **Every reader can tell "none" from "not scanned"** through `domain_state`, and exports print `-` and `?`.
4. **An annotated stage also scans the proteins of contigs a called stage owns in its genomes**, copying measured
   ones.
5. **Clustering: measure linclust alone against the cascade before choosing** (C-38c). The cascade cannot cluster a
   release that includes C-42's pilot on the pod; nothing we have clusters all of SPIRE at once, and DEC-24 rules out
   the split clustering that would. That limit is the PI's (Q10).
6. **The filter is a least contig length at ingest**, applied with the inputs by one rule, before the held-contig
   check; its records and bases go in ingest's ledger through two appended columns. L is chosen from C-42.
7. **Nucleotides kept**; storage planned at 0.714 GB per Gbp or more.
8. **Order of work:** C-38a (depth), C-38b (domain state), C-38d (filter) can start now; C-38c (the clustering
   measurement) needs `r26.11.1`; C-38e and C-38f (rounds) after C-38b; C-38g after C-38e; C-42 measures what the
   thresholds need. None of these blocks C-39 to C-41, which build genome-like stages at depth annotated.

Risks: the round's scope rule is the first new release scope since card B-33, and every reader of `l1_domains` must
learn round rows (tokens, postings, `nbr`, the exports, the catalog views, `domains.sources`); one missed reader shows a
protein as not scanned or, if it reads rounds into family tokens, biases them (test (6) of §3.8). Linclust's loss is
unmeasured; if it is large, the pilot has no clustering that fits the pod.

## 9. Questions only the PI can answer

1. Is "clustered" a property of a release (DEC-24), so a stage's depth is stored, called or annotated? *yes / no*
2. Does a release take only stages at depth called or more? *yes / no*
3. Does a release's search pack hold the IGRs and dark ORFs of its called stages too (invariant 16 as written)?
   *yes / proteins only*
4. The window of a round: genes on each side of a hit? *10 / 5 / other*
5. Are round rows kept out of family tokens and postings, and shown only in neighbour tables and exports? *yes / no*
6. Does a release take every round complete when it is cut, with no release cut for rounds alone? *yes / no*
7. Do rounds run only as build-account steps, never inside a user's search? *yes / no*
8. Does an annotated stage scan the proteins of contigs that a called stage owns in its genomes? *yes / no*
9. After card C-38c, one clustering method for the whole release line, chosen then? *yes / no*
10. When a release's proteins cannot be clustered at once on any machine we have, does DEC-24 stand (no such
    release), or is it revisited? *stand / revisit*
11. Is a contig dropped by the length filter outside the corpus (D-01's filters), so invariant 2 does not cover it?
    *yes / no*
12. The least contig length, chosen from C-42's measured shares: *1.5 kb / 5 kb / other*
13. Filter by length only, not by gene count? *yes / no*
14. Pfam on representatives of large NR50 clusters as rounds, members shown "by cluster" at query time and never
    written as inferred rows? *yes / no*

## 10. Cards that follow from the recommendation (text for `TASKS.md`)

```
### C-38a Stage depth: chosen at build, derived for a release
**M** · Depends: C-38 (PI answers Q1, Q2) · Design: `planning/design/C-38_import_depths.md` §2
Files: `src/biobox/stage.py` (`--depth stored|called|annotated` on `Build` and the step list; refusals), `src/biobox/catalog.py` (`stage_depths(manifest)`; a release refuses a stage without tapes, naming its depth), `src/biobox/cli.py` (`build-stage --depth`; `release verify` prints a `depths:` line in place of card C-19's warning; `bio-box releases` column `depths`), `RUNBOOK.md` (build-stage usage), tests in `tests/r1` and `tests/r2`
Do, test first: (1) the step list at each depth; `--depth called` with `--library` or `--ncrna-cm`, and `--depth stored` given a release step, refused (exit 2, nothing run). (2) A stage built called, then scanned with the library, has files byte-identical to the same stage built annotated at once; fixture manifest hash compared. (3) `stage_depths` on every test manifest and on a copy of `rM.2`'s file list (all seven annotated pfam-38.2); on a fixture release over one called and one annotated stage, family tokens equal a brute-force count over the annotated stage's domain files written separately in the test. (4) A release over a stage at depth stored is refused. (5) The function the release train (card C-46) uses to judge a stage complete takes the depth.
Accept: the tests; `release verify` and `bio-box releases` on the pod, read-only, give the same depths for every published release.

### C-38b Domain state: "none" is not "not scanned"
**M** · Depends: C-38a · Design: §3.4
Files: `src/biobox/catalog.py` (view `domain_state(protein_uid, library)`: measured, inferred, none, not scanned; stage coverage from the manifest's annotated stages, round coverage once C-38f exists), `src/biobox/search/export.py`, `src/biobox/search/nbr.py` (`pfam_architecture` `-` for none, `?` for not scanned; a `domain_state` column in hits and locus tables; README text), `GUIDE.md` where it explains `pfam_architecture`, tests
Do, test first: a fixture release with a called stage and an annotated stage, and a protein in both: each state checked against a brute-force reading of the files written separately in the test; the exports of a release whose stages are all annotated are byte-identical but for the new column (W-07's byte-identity test migrated for that column only); a lookup of one protein's state reads O(1) row groups (`profile_query`).
Accept: the tests.

### C-38c Linclust alone against our cascade, on `r26.11.1`'s proteins
**M** · Pod: yes (read-only on the data root; scratch only; at most 16 threads, under the hour's cap) · Depends: C-35 (`r26.11.1` built) · Design: §4.2, §4.3
Files: `scripts/stress/c38c_*` (new), one measurement note with its directory
Do: (1) In scratch, cluster `r26.11.1`'s distinct proteins with the linclust calls of `l2/cluster.workflow`'s redundancy phase alone, at NR90 and at NR50 over the NR90 representatives, with the release's parameters, in parts under the 300 CPU-minute budget (card B-41's cut). (2) Report against the release's clusters: clusters and singletons per level; for each cascade NR50 cluster, how many linclust clusters its members fall in; members the cascade joined and linclust did not, by identity band; CPU-minutes and peak RSS per call, against N19 and N20 of the design note; the memory per protein. (3) Run the two worked examples' neighbour summaries (`examples/`, the TniQ and RNA polymerase searches) at NR50 with each clustering and report what changes in their tables. (4) From (2), the CPU and memory at 350 million and at 1 billion proteins, labelled as extrapolations.
Accept: the note; the 20% rule of AGENTS §6 applied to N19 and N20.

### C-38d The least contig length at ingest
**M** · Depends: C-38 (PI answers Q11, Q13); the value of L from C-42 (Q12), not needed to build the mechanism · Design: §5
Files: `src/biobox/l0seq/store.py` (ingest drops records under `--min-contig-len`, before the held check; refuses inputs that disagree), `src/biobox/stage.py` (`build-stage --min-contig-len` passed to ingest; the ingest estimate unchanged), `src/biobox/registry/schemas.py` (`ledger.n_bases_in`, `ledger.n_bases_rejected`, appended, nullable), the locks by `scripts/update_schema_lock.py`, `src/biobox/ledger.py` (the reason `below:<L>`; L-01 unchanged in records), `scripts/stress/build_stage_pod/make_inputs.py` (the same rule for the inputs it writes), `RUNBOOK.md`, tests
Do, test first: §5.3's tests (1) to (5).
Accept: the tests; the fixture manifest hash unchanged without the option.

### C-38e Rounds: tables, scope and release listing
**M** · Depends: C-38 (PI answers Q5, Q6), C-38b · Design: §3.2, §3.8
Files: `src/biobox/registry/schemas.py` (`round_domains`, `round_coverage`, partitioned (`round`, `library`) with alternative keys (`release`, `library`); `release_round`, partition `release`, `small_ok`; the ledger's key set (`round`, `step`); an artifact for a round's plan and shard records, `in_manifest` False), the locks by `scripts/update_schema_lock.py`, `src/biobox/catalog.py` (scope rule (e) in `_scope`, `_in_release`, `ReleaseFiles._layout`, `build_manifest`, `verify_manifest`; `define_release` writes `release_round` with every complete round), tests
Do, test first: §3.8 tests (4) and (5) on hand-written round files; every published-release manifest copy in the tests verifies unchanged; a release built before any round ignores round files; a round file no release lists is not reported by any release's verify.
Accept: the tests.

### C-38f The round step: scan a list of proteins, copy what was measured
**M** · Pod: yes (tests; one round on a scratch release) · Depends: C-38e, B-47b · Design: §3.3, §3.6
Files: `src/biobox/l1/domains.py` (a round's plan, shards and merge over a requested list; the round id; rounds as sources in `sources`), `src/biobox/cli.py` (`bio-box annotate round --release R --library NAME VERSION HMM CUTOFF --proteins FILE`; `search hmm --annotate`), `src/biobox/stage.py` only to launch the shards as scan shards are launched, `RUNBOOK.md`, tests
Do, test first: §3.8 tests (1), (2) and (7); the request holds only proteins not scanned in the searched release; the plan copies rows measured by a stage or an earlier round on disk and scans the rest; ledger rows as §3.3; each process under the budget.
Accept: the tests; on the pod, in a scratch view, one round over the neighbours of one search of the TniQ example, with its CPU-minutes against §3.8's estimate.

### C-38g Readers of rounds: the release's merged copy, tokens kept clean
**M** · Depends: C-38f · Design: §3.2 (lookups), §3.5
Files: `src/biobox/l1/domains.py` or a new module for the release step `rounds` (the merged, sorted copy under `release=`), `src/biobox/stage.py` (the step after `cluster`), `src/biobox/catalog.py` (`domain_state` and `domains_of` read the merged copy; a `source` column), `src/biobox/search/export.py`, `src/biobox/search/nbr.py`, `src/biobox/l3/tokens.py` and `src/biobox/l3/postings.py` (assert they do not read round rows), tests
Do, test first: §3.8 tests (3), (5) and (6); `nbr summarize` shows round rows with `source round`.
Accept: the tests; the fixture manifest hash of a release without rounds unchanged.

### C-38h An annotated stage scans the proteins of contigs a called stage owns
**M** · Depends: C-38 (Q8), C-38a, B-47b · Design: §3.7
Files: `src/biobox/l1/domains.py` (the plan adds those proteins, read from the owner's `feature_protein`; copies measured ones), `src/biobox/stage.py` (passes the stages on disk, as `--seq-stages`), tests
Do, test first: a called stage storing a contig that a later annotated stage's genome holds: the annotated stage's domain file holds rows for that contig's proteins, equal to an `hmmsearch` run of them; built in either order the release's family tokens are the same.
Accept: the tests.
```

A card for Pfam on representatives of large NR50 clusters (§4.4) waits for Q9 and Q14 and for C-38c's cluster sizes;
it would be a round whose request is a list of representatives (C-38f's step) and a query-time view in the catalog.

## 11. Not determined here, and what settles each

1. Genes, proteins and bytes per Gbp of unbinned contigs, gene-calling cost in meta mode, the share of genes that are
   fragments, the shares of bases over 1.5, 5 and 10 kb: card C-42 (3), (4). Until then N5 to N10 rest on CAMI's
   simulated gold-standard assembly and on SPIRE's published totals.
2. The value of L: Q12, from C-42 (3).
3. What linclust alone loses against the cascade, and clustering CPU and memory above 45.8 million proteins: C-38c.
4. New NR50 families per Gbp per habitat, the input of a stopping rule: C-42 (5), after Q9.
5. The cost of a round end to end and of a search over the pilot: C-42 (6) and C-38f.
6. Whether a Lonestar6 core does what a pod thread does, and results from another host entering the root: C-37.
7. What Logan holds per run and the cost of searching it: C-44.
8. The pack's bytes per Gbp of unbinned contigs, which decides Q3 in practice: C-42 (4) with a pack built in scratch.
