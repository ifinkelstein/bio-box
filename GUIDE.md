# Searching a release with your own HMM

For a scientist who wants hits and gene neighbourhoods, not operations. State as of 2026-11-11; releases,
paths and costs will change. The operator's version of every command is `RUNBOOK.md` (sections 0, 7, 10, 11).
Numbers carry a label: *measured* (with the file it comes from) or *derived* (with the formula).

## 1. What this is

A release is a fixed snapshot of a set of genomes: genes called, identical proteins stored once, proteins
clustered at 90% and 50% identity (NR90, NR50) and, in some releases, Pfam domains; once built it never changes.
A search runs your profile HMM (HMMER3 `hmmsearch`, through pyhmmer) against every distinct protein of one
release and keeps the proteins at or above your bit-score cutoff. It will not run without a positive control, a
protein you know it must find, and it also searches a shuffled copy of your model, which must find nothing. For
every gene that carries a hit protein it reads the genes on either side, in contig order. You get a hit table, the
hit sequences, one table per hit locus, and counts of which clusters (and Pfam families) sit next to your hits.

What exists, on the pod `wilkcomp02`, in the data root `/tmp/if2453/r10/root` (a search pack is a read-only copy
of a release's proteins laid out for searching, one per release, shared; the pod's packs belong in one place, the
host profile's pack directory `/tmp/if2453/bio-box/packs`, card W-22):

| Release | Contents | Search pack |
|---|---|---|
| `rN.1` | 27,258 phage genomes and 2,497 complete prokaryote genomes; 12,048,866 genes; 10,620,684 distinct proteins; NR50 2,966,200 and NR90 8,278,537 clusters; no Pfam (*measured*, `planning/tasks/STATUS.md`, "10 Gbp rebuild"; `planning/measurements/2026-11-10_trial_search_rN1.md`) | yes; until the operator copies it into the profile's pack directory, in `/tmp/if2453/w19/packs` (section 2) |
| `rN.2` | the same genomes with Pfam 38.2 domains | none yet: ask the operator for one (RUNBOOK section 7, `pack build`) |
| `rN.3` | the same genomes, built by the current code (card C-17) | being built |

## 2. Before you start

- Everything runs on the pod, in `bash`, inside the project's container image. Paste the block of RUNBOOK
  section 0 once per shell, with every value marked `real:`. It defines `bio-box` (one short command in its own
  memory-capped scope) and `$TREE`, `$IMAGE`, `$DIGEST`, `$PYTHONPATH_POD`, `$PROFILE`, `$BBTMP`.
- `$TREE` must contain card W-21 (merged into `plan/task-cards` on 2026-11-11); `main` and `~/bio-box` on the pod
  are older. Ask the operator which tree is current.
- The pod is shared: a search goes through `scripts/pod_run.sh` with at most a few threads (4 below). The 40
  threads by day and 56 at night are a total over everything we run (RUNBOOK section 10).

Then your own directories; keep every file you pass under `/tmp/if2453`, which is bound into the image.

```bash
export ROOT=/tmp/if2453/r10/root
export WORK=/tmp/if2453/pi-search
mkdir -p "$WORK"
cp /tmp/if2453/trial2/work/RNA_pol.hmm "$WORK/"   # Pfam PF00940, T7-like RNA polymerase: the worked example
```

Where the search pack is: the commands below pass `--profile "$PROFILE"` and so find a release's pack in the
profile's pack directory, `/tmp/if2453/bio-box/packs` (`$BBTMP/packs`; card W-22); you never name it. Until the
operator has copied `rN.1`'s pack there (as of 2026-11-12 it is in `/tmp/if2453/w19/packs`), add `--pack-dir
/tmp/if2453/w19/packs` to each `releases`, `search` and `nbr` command; without it a search says `no search pack
of release rN.1 in the pack directory /tmp/if2453/bio-box/packs` and exits 2, and `releases` shows `-` for the
pack. The pack directory may be read-only for you: a search writes nothing into it.

The table in section 1 will age; `releases` lists what the root holds now, one line per release: its stages,
genomes, contigs and distinct proteins, its domain libraries (Pfam or none), whether it has the lookup tables of
the current code (`member_cluster`, `contig_label`, `domain_family_tokens`; `no` means slower lookups, section 8),
and the search pack it has in the profile's pack directory (card W-23):

```bash
bio-box releases --root "$ROOT" --profile "$PROFILE"
```

## 3. Choose a positive control

A positive control is one protein of the release that your model must find; if the run misses it, the run ends
`failed` (`positive control not recovered`) with no tables, so a wrong cutoff or model shows up as a failure, not a
short list. Pick a protein whose family you trust and that scores well above the cutoff.

You name it by one gene, in either of two forms. GenBank's, copied from a record: `ACCESSION:START..END`, 1-based
with both ends included, so T7 gene 1 (GenBank `NC_001604`, 3171..5822) is `NC_001604:3171..5822`. Or bio-box's
own, which every table it writes uses: `ACCESSION:START-END`, 0-based and half-open (`START` is the GenBank start
minus 1, `END` the GenBank end), so the same gene is `NC_001604:3170-5822`. Two dots or one dash decide which;
anything else (`3171.-5822`, `3171...5822`) is refused with both spellings of what you may have meant. Find it with
`find`:

```bash
bio-box find --root "$ROOT" --release rN.1 --accession NC_001604
bio-box find --root "$ROOT" --release rN.1 --accession NC_001604:3171..5822
bio-box find --root "$ROOT" --release rN.1 --protein 2eef240afefddcd2c8ac7660eb50e384
```

Columns: accession, contig_uid, start, end, strand, feature_uid, protein_uid, then the gene in both forms, `span`
(`NC_001604:3170-5822`) and `genbank_span` (`NC_001604:3171..5822`). The first lists the 49 genes of the T7 genome;
the second the one gene, with its protein uid; the third every gene that carries that protein (33, in 33 T7
isolates). Each took about 4 s on `rN.1` and printed an expected `note: ... slow path` line (*measured*,
`planning/measurements/2026-11-11_trial_search_rN1_second.md` section 1, with the 0-based form; section 8 below).

`--positive` also accepts text from a contig's header (case-sensitive), such as `'Enterobacteria phage T7'`. On a
real release that text names many contigs and every gene on them, so the command stops with exit 2 and lists up to
20 matching genes, each as `ACCESSION:START-END(strand)`, the same as `ACCESSION:START..END`, protein uid and
header. Copy one gene from that list.

## 4. Run a search

One command, in a capped unit (`pi-rnap` is the unit's name; use a new name for every launch):

```bash
bash "$TREE/scripts/pod_run.sh" --profile "$PROFILE" --image "$IMAGE" --digest "$DIGEST" --mem 32G --threads 4 --bind /tmp/if2453 --bind /tmp/if2453/r10:/tmp/if2453/r10:ro --env PYTHONPATH="$PYTHONPATH_POD" pi-rnap search hmm --root "$ROOT" --release rN.1 --profile "$PROFILE" --hmm "$WORK/RNA_pol.hmm" --bits 30 --neighbors 5 --positive NC_001604:3170-5822 --out "$WORK/rnap" --cpus 4 --mem-gb 30
```

Before `pi-rnap` are the unit's limits: 32 GB of memory, 4 CPUs, the image, and the release bound read-only so
nothing can write into it. After it, the search:

| Flag | Meaning |
|---|---|
| `--root`, `--release` | the data root, and the release in it to search |
| `--profile` | the host profile: thread and memory caps, where runs are kept, and the pack directory where the release's search pack is (`/tmp/if2453/bio-box/packs`) |
| `--hmm` | your HMM file (HMMER3 format) |
| `--model` | the model to use, by name or accession, when the file holds several |
| `--bits` | inclusion cutoff, in bits, on the full-sequence score (E-value cutoffs are refused) |
| `--neighbors` | genes on each side of each hit in the neighbour tables (default 5) |
| `--positive` | the positive control, `ACCESSION:START..END` (GenBank), `ACCESSION:START-END` or header text (section 3) |
| `--out` | a new or empty directory, outside the data root, for the results |
| `--cpus`, `--mem-gb` | search threads, and its memory limit in GB; never above the profile's caps |
| `--library`, `--threads` | the domain library for `pfam_architecture` when there are several; export threads (default 4) |
| `--runs`, `--work`, `--pack-dir` | where runs and their scratch go, and where the pack is, if not where the profile says (section 2) |
| `--cpu-min`, `--wall-s` | limits of the search itself (defaults 300 CPU-minutes, 6 hours) |

The worked example is this exact search. It found 1,957 proteins at 30 bits or more, 1,931 in phages and 26 in
prokaryotes, on 2,581 genes; T7 gp1 scored 367.7 bits; the same 1,957 proteins as `hmmsearch` (HMMER 3.4) run by
hand on all exported proteins, scores equal to 0.1 bit (*measured*, second trial file, section 3).

Cost with the current code, 4 CPUs on a loaded pod: 244.6 s for the whole command, 9.53 CPU-minutes, the search
process at most 1.78 GiB and the unit 4.57 GB including file cache (*measured*,
`planning/measurements/2026-11-11_search_memory_rN1.md`, addendum W-21). Memory no longer grows with the
release (card W-19). A pack without a verification record is checked in full first, about 130 s more (same
file, section 4).

The log is `$BBTMP/logs/pi-rnap.log`. About 6 s after launch it has a line `search run: <run id> running`
(*measured*, same file); that id names the run. Watch and stop:

```bash
tail -n 5 "$BBTMP/logs/pi-rnap.log"
grep "search run:" "$BBTMP/logs/pi-rnap.log"
RUN=87021d87e0f41e453c14d35763d51d17   # the id from your own log line; this one is the trial's
bio-box search status "$RUN" --root "$ROOT" --profile "$PROFILE"
bio-box search cancel "$RUN" --root "$ROOT" --profile "$PROFILE"
```

`status` prints `"status": "running"`, then `ok` or another status (`failed`, `stopped`, `oom`, `timeout`) with
the reason. `cancel` (or `systemctl --user stop pi-rnap`) ends the run `stopped`. A finished log has one line per
output file with what it holds (`wrote .../hits.tsv: 2581 rows`, `.../hits.faa: 1957 sequences`, `.../loci: 2581
tables, ... rows`; every command given `--out` prints such a line), then the run id and `ok`, one line per file, and
`pod_run: pi-rnap exited 0`. `exited 2`: refused before the search (positive
control, missing pack, `--out` not empty; the log says which). `exited 1`: the run ended other than `ok`, and
`--out` holds only `recipe.yaml`.

## 5. What comes out

In `--out` (`$WORK/rnap`); the run itself (its record and raw hits) stays in
`$BBTMP/search/$USER/runs/run=<run id>/`:

| File | What it is for |
|---|---|
| `README.txt` | every file and column below, written by the command itself |
| `recipe.yaml` | the exact search, to run again or change (section 7) |
| `hits.tsv` | one row per gene that carries a hit protein (2,581 rows for 1,957 proteins in the example) |
| `hits.faa` | one sequence per hit protein, named by protein uid |
| `neighbors.tsv` | which clusters (and Pfam families) sit near the hits, counted over all hit loci |
| `loci/` | one table per hit gene and its neighbours: `locus_<accession>_<feature_uid>.tsv` |

Conventions in every file. Coordinates (`start`, `end`) are 0-based and half-open: the first base is 0 and `end` is
the first base after the gene, so a gene is `end - start` bases long. Strand is `+1` or `-1`. On a circular contig
a gene across the origin has `end` greater than the contig length. Every id is 32 hex digits. A `feature_uid`
names one gene; an `ordinal` (a gene's position in its contig's gene list) is valid only inside its `tape_version`,
so use `feature_uid` to name a gene.

`hits.tsv`, the columns you will use:

| Column | Meaning |
|---|---|
| `accession`, `header` | the contig's source record (the first, when identical sequences were merged) |
| `start`, `end`, `strand`, `contig_len`, `length` | the gene; the contig's length in bases; the protein's in amino acids |
| `seq_bitscore` | full-sequence score in bits, one decimal |
| `seq_evalue_zref` | E-value at a fixed reference database size of 10^8 sequences, so it compares across releases; at this release's 10,620,684 proteins the E-value would be about 0.106 times it (*derived*: 10,620,684 / 10^8) |
| `pfam_architecture` | Pfam domains along the protein, joined by `+`, a repeated domain written once; `-` for none, and `-` everywhere on a release without Pfam (`rN.1`) |
| `protein_uid` | the protein; identical sequences share it |
| `nr90`, `nr50`, `nr50_rep` | the protein's clusters at 90% and 50% identity (ids of this release only, not comparable across releases), and the NR50 representative's protein uid (its sequence: section 7) |
| `dist_to_contig_end_bp` | 0: the gene caller marked the gene partial at a contig end (it runs off it); otherwise the distance to the nearer end; the contig length for a circular contig |
| `n_records` | how many input records carried this identical contig sequence; not a hit count |
| `genome_uid`, `contig_uid`, `feature_uid`, `ordinal`, `tape_version` | the genome, the contig and the gene |
| `stage`, `kind` | the build stage the gene comes from (on `rN.*`: `S1` phages, `P1` and `P2` prokaryotes) and its genome's kind (`virus`, `prokaryote`, `plasmid`), so hits split by source with a filter on one column |

A locus table (`loci/*.tsv`): line 1 the contig header, line 2 facts (`anchor_feature_uid`, `anchor_protein_uid`,
`censored=true` when a contig end cut the window), then one row per gene. `offset` counts genes from the hit in
coordinate order (lower coordinates negative); `offset_rel` counts along the hit's own strand (positive is
downstream of the hit), so on a `-1` hit the two have opposite signs. `hit_seq_bitscore` is filled when that
neighbour is itself a hit; `delta_bp` is the gap to the hit in bases (no sign); `crosses_origin` marks a gene across
the origin.

## 6. Reading neighbourhoods

A locus is one hit gene with its window; the hit is not its own neighbour. `neighbors.tsv` starts with
`# key<TAB>value` totals (example: 2,581 loci on 2,526 contigs, 463 cut by a contig end, 24,473 neighbour genes;
*measured*, second trial file, section 3), then one table per level (`pfam`, `nr50`, `nr90`), 50 rows at most.

| Column | Meaning |
|---|---|
| `key`, `cluster_uid` | the cluster (the same value twice), or the Pfam family at the `pfam` level |
| `rep_uid`, `rep_length` | the cluster's representative protein and its length |
| `model_name` | on a release with Pfam, the domain name most members carry; empty otherwise |
| `n_loci`, `share_of_loci` | hit loci with at least one such neighbour, and the share of all loci |
| `n_neighbors` | (locus, neighbour gene) pairs; a gene in the windows of two hits counts twice |
| `n_genomes` | genomes with such a neighbour |
| `offsets`, `modal_offset` | `offset:count` pairs and the commonest offset, along the hit's strand (as `offset_rel`, though the name lacks `_rel`) |
| `same_strand_share` | the share of those neighbours on the hit's strand |
| `unannotated` | `true` on the row of neighbours with no cluster or family at that level |

Read rows by `n_loci`. In the example the top NR50 cluster (`fa4cb451...`) is in 663 of 2,581 loci (25.7%,
*derived*), commonest at +4; the third (`5fe1410c...`) sits at -1 in 476 of its 481 loci (*measured*, second trial
file). On `rN.1` (no Pfam) a cluster is only an id. To name one, look up its representative: `bio-box find
--protein <rep_uid>` gives an accession and coordinates, and `bio-box seq --protein <rep_uid>` the sequence
(section 7). In T7's own
locus the -1 and +3 neighbours fall where T7's protein kinase (gene 0.7) and DNA ligase (gene 1.3) are, so those two
clusters are probably those families; no file says so. On `rN.2` the neighbours carry Pfam names.

These are counts, not significance: no enrichment statistic exists yet, and identical genomes count separately
(T7 gp1 alone is in 33 genomes). Look at one locus:

```bash
ls "$WORK"/rnap/loci/locus_NC_001604_*
column -t -s $'\t' "$WORK/rnap/loci/locus_NC_001604_39f9ef263e234b012d2671634468546e.tsv"
```

## 7. Going further with the separate commands

`search hmm` is these commands in a row; each can be run alone on a finished run, with its own window or format.
A higher cutoff needs no new run: filter `seq_bitscore` in `hits.tsv`. A lower one does. Copy the recipe, change
`inclusion_bits`, and run it (a new unit name; `--wait` is required inside a unit):

```bash
cp "$WORK/rnap/recipe.yaml" "$WORK/rnap25.yaml"
sed -i 's/inclusion_bits: 30.0/inclusion_bits: 25.0/' "$WORK/rnap25.yaml"
bash "$TREE/scripts/pod_run.sh" --profile "$PROFILE" --image "$IMAGE" --digest "$DIGEST" --mem 32G --threads 4 --bind /tmp/if2453 --bind /tmp/if2453/r10:/tmp/if2453/r10:ro --env PYTHONPATH="$PYTHONPATH_POD" pi-rnap25 search run "$WORK/rnap25.yaml" --root "$ROOT" --release rN.1 --profile "$PROFILE" --cpus 4 --mem-gb 30 --wait
```

Then, with `RUN` set to the new run id from `$BBTMP/logs/pi-rnap25.log`, the tables, the hit sequences, and the
neighbours with a wider window (`10g` is 10 genes on each side; `5kb` is 5,000 bases):

```bash
bio-box search hits "$RUN" --root "$ROOT" --profile "$PROFILE" --format tsv --out "$WORK/rnap25_hits.tsv"
bio-box search hits "$RUN" --root "$ROOT" --profile "$PROFILE" --format faa --out "$WORK/rnap25_hits.faa"
bio-box nbr summarize "$RUN" --root "$ROOT" --profile "$PROFILE" --window 10g --by nr50,nr90 --out "$WORK/rnap25_neighbors_10g.tsv"
bio-box nbr "$RUN" --root "$ROOT" --profile "$PROFILE" --window 10g --format tsv --out "$WORK/rnap25_loci_10g"
```

Each took 6 to 13 s and at most 1.3 GB of memory on `rN.1` (*measured*, second trial file, section 2). On a
release with Pfam use `--by pfam,nr50`. `nbr summarize --top 0` lists every cluster, not 50.

Below: one protein's sequence and its genes (`seq` reads the protein by its uid: it takes several uids, or
`--proteins-file` with one per line, and `--out FILE`); a stretch of a contig, by accession in either coordinate
form, or by contig uid, `--strand -1` for the reverse complement (an end past the contig's length crosses the
origin of a circular contig); every protein of the release as FASTA named by protein uid, which you need only to
search all of them with another tool (88.2 s, 3.7 GB, 2.28 GB of memory on `rN.1`, *measured*,
`planning/measurements/2026-11-10_trial_search_rN1.md`); and a check that the release is intact (it re-reads every
file; about 60 s on `rN.1`, *measured*, the same file).

```bash
REP=2eef240afefddcd2c8ac7660eb50e384   # T7 gp1 here; put a rep_uid from neighbors.tsv or nr50_rep from hits.tsv
bio-box seq --root "$ROOT" --release rN.1 --protein "$REP"
bio-box find --root "$ROOT" --release rN.1 --protein "$REP"
bio-box seq --root "$ROOT" --release rN.1 --contig NC_001604:3171..5822
bio-box seq --root "$ROOT" --release rN.1 --contig NC_001604:3171..5822 --strand -1 --out "$WORK/t7_gene1_rc.fna"
bio-box export proteins --root "$ROOT" --release rN.1 --out "$WORK/proteins.faa"
bio-box release verify --root "$ROOT" --release rN.1
```

A FASTA header from `seq` names the uid and the release, and for a stretch of contig its accession, `start` and
`end` (0-based, half-open), `strand` and the contig's length. Any other question about the release (how many hits
per genome, which genes sit on a contig, what a cluster holds) is a SELECT on its catalog with `bio-box sql`; the
agent guide (`AGENT_GUIDE.md` section 5) lists the tables and gives examples. It runs nothing that writes, and
says on stderr how much of each table each statement read:

```bash
bio-box sql --root "$ROOT" --release rN.1 -e "SELECT key, value FROM catalog_info ORDER BY key"
```

## 8. Limits and rough edges today

- The data is on the pod's local disk (`/tmp`), not yet copied to permanent storage; a reboot may clear it. No card
  yet (`STATUS.md`, "Plan for larger imports", step 1).
- Only `rN.1` has a search pack. Until the operator copies it into the profile's pack directory
  (`/tmp/if2453/bio-box/packs`, card W-22), it is in a work directory (`/tmp/if2453/w19/packs`) and every
  `releases`, `search` and `nbr` command needs `--pack-dir /tmp/if2453/w19/packs` (section 2); `rN.2` needs one
  built (RUNBOOK section 7).
- `rN.1` and `rN.2` were built before card B-27, so cluster and contig-name lookups read whole tables and print a
  `note: ... slow path` line; answers are the same, a few seconds slower. `rN.3` will not need it (card C-17).
- Looking up `--positive`, `find --accession` or `seq --contig ACCESSION` reads every contig name of the release,
  about 4 s on `rN.1` (card W-15); `seq --contig` with a contig uid does not.
- Clusters have no names or products; on a release without Pfam a neighbour is only an id (no card yet).
- NR50 and NR90 ids change between releases when a cluster gains a member (PI card 0-07).
- A contig whose identical sequence appears in two genomes of different build stages belongs to the first genome
  that brought it, with that genome's gene calls (DEC-26; card C-21).
- No non-coding RNA genes are in any release yet: neighbours are protein-coding genes only (card E-07, PI card 0-08).
- No statistics for neighbour enrichment: counts are counts, not significance, and redundant genomes are not
  down-weighted (cards F-01 to F-03).
- `release verify` on `rN.1` exited 1, naming Pfam files added later for `rN.2`; `rN.1` itself is unchanged. Seen
  with the code that built it; whether current code still does it is card C-19.

## 9. Where to ask, and the operator documentation

Ask the operator who runs the builds for a pack, a release, or a current `$TREE`; a failure is worth a task card
in `planning/tasks/TASKS.md`. The operator's documents: `RUNBOOK.md` (every command, the limits, recovery and
clean-up), `AGENTS.md` section 7 (rules for the shared pod), `planning/tasks/STATUS.md` (what exists now), and the
`README.txt` that each search writes beside its results.
