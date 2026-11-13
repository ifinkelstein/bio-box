# bio-box (all-the-data)

**Status, 2026-11-13:** three releases over the 10 Gbp set exist (search `rN.3`), the 50 Gbp set is clustered and
its Pfam scans are running, and no release has been signed yet. The running record is
[`planning/tasks/STATUS.md`](planning/tasks/STATUS.md); where this file and that one differ, that one is newer.

Numbers here carry the labels [`AGENTS.md`](AGENTS.md) §6 requires: *measured* (on the pod unless stated, with the
date and the record cited beside it), *derived* (with its formula) or *estimate*.

## What bio-box is

bio-box is a staged, content-addressed and reproducible database of microbial and phage genomes: their sequences,
gene calls, proteins, Pfam domains, non-coding RNA, protein clusters and gene order. It is built to answer discovery
questions by genomic context: search the whole collection for a protein family with your own profile HMM, then see
what sits next to every hit, which neighbours recur, and which loci share a context.

It is for scientists who want hits and gene neighbourhoods ([`GUIDE.md`](GUIDE.md)) and for LLM agents that answer
their questions with it ([`AGENT_GUIDE.md`](AGENT_GUIDE.md)). It is not a web service: it is a command-line tool that
runs on a shared university server, where no network service is started ([`AGENTS.md`](AGENTS.md) §7.5). It is not a
general annotation pipeline: it annotates only the genomes in its own stages, with tools and reference libraries
frozen for the release line (DEC-12, [`refdb.yaml`](refdb.yaml)), and has no command to annotate a genome you bring.

No licence has been granted (DEC-07). Data follow each source's licence ([`sources.yaml`](sources.yaml)).

## What you can do with it today

Every command runs on the shared server inside the project's container image, through the shell function and the
capped-unit launcher that [`RUNBOOK.md`](RUNBOOK.md) §0 and §7 define. [`GUIDE.md`](GUIDE.md) walks through each one.

| Command | What it does |
|---|---|
| `bio-box search hmm` | runs your HMM (HMMER3, through pyhmmer) against every distinct protein of a release, with a required positive control and a shuffled-model negative control; writes hits, their FASTA, a neighbour summary, one table per locus and a recipe to rerun (GUIDE §4, §5) |
| `bio-box nbr` | the genes around each hit of a run, by a window in genes or bases (GUIDE §7) |
| `bio-box nbr summarize` | which NR50 and NR90 clusters and Pfam families sit near the hits, counted over loci (GUIDE §6) |
| `bio-box nbr flanks` | groups a run's loci by their flanking genes into context types, with a UPGMA tree per group; `--figure` draws a synteny figure per group (GUIDE §6) |
| `bio-box find` | the genes an accession (GenBank or 0-based coordinates) or a protein uid names (GUIDE §3) |
| `bio-box seq` | a protein by uid, or a contig slice by accession and coordinates, as FASTA (GUIDE §7) |
| `bio-box sql` | read-only SELECT statements on a release's catalog (GUIDE §7; AGENT_GUIDE §5) |
| `bio-box releases` | the releases of a data root: stages, counts, libraries, lookup tables, search pack (GUIDE §2) |
| `bio-box results` | `archive`, `list` and `get` keep a finished run on durable storage (DEC-33; GUIDE §5). The durable directory does not exist yet, so these wait on an open item for the PI |

Also there: `search run|status|cancel|hits`, `export proteins` and `release verify` (GUIDE §7).

### Worked example: RNA polymerase on `rN.3`

Pfam's T7-like RNA polymerase model (`RNA_pol`, PF00940) at 30 bits, with T7 gene 1 as the positive control. The
search ran in a unit capped at 4 CPUs and 32 GB, launched as GUIDE §4 shows:

```bash
bio-box search hmm --root "$ROOT" --release rN.3 --profile "$PROFILE" --hmm RNA_pol.hmm --bits 30 --neighbors 5 --positive NC_001604:3170-5822 --out rnap --cpus 4 --mem-gb 30
bio-box nbr flanks "$RUN" --root "$ROOT" --profile "$PROFILE" --window 5g --level nr50 --figure --out rnap_flanks
```

| Result | Value | Label and record |
|---|---|---|
| proteins at 30 bits or more | 1,957, on 2,581 genes | *measured*, 2026-11-11, [`2026-11-11_build_r50.md`](planning/measurements/2026-11-11_build_r50.md) §A3 |
| the same set by `hmmsearch` (HMMER 3.4) run by hand on every exported protein | 1,957, none different; scores equal to 0.1 bit; 1,931 in phages, 26 in prokaryotes | *measured* on `rN.1`, 2026-11-11, [`2026-11-11_trial_search_rN1_second.md`](planning/measurements/2026-11-11_trial_search_rN1_second.md) |
| T7 gene 1 | 367.7 bits, architecture `RPOL_N+RNA_pol`; shuffled model 0 hits | *measured*, 2026-11-11, `2026-11-11_build_r50.md` §A3 |
| neighbourhoods, 5 genes each side | 2,581 loci on 2,526 contigs, 463 cut by a contig end, 24,473 neighbour genes | *measured*, same record |
| commonest Pfam neighbours, by loci | Phage_T7_1_1 1,032; DNA_ligase_A_M 1,014; DNA_ligase_C 989; Head-tail_con 753; Phage_T7_kinase 680 | *measured*, same record |
| cost of `search hmm` | 319 s for the whole command; search worker peak 2.51 GiB | *measured*, same record |
| `nbr flanks` at NR50 | 432 groups, 598 context types, none holding both phage and prokaryote loci; 15.0 s and 1.41 GB without the figures, 273.0 s with them | *measured*, 2026-11-12, [`2026-11-12_flanks_rN3.md`](planning/measurements/2026-11-12_flanks_rN3.md) |

The figure of the T7-like group is
[`nr50_4a0062a1_t7_like.svg`](planning/measurements/2026-11-12_flanks_rN3/figures/nr50_4a0062a1_t7_like.svg). Counts
are counts, not significance: no enrichment statistic exists yet, and identical genomes count separately (GUIDE §6).

## What is in it

The data are not in this repository. Built stages and releases live in one data root on the local disk of a shared
university server (the pod, `wilkcomp02` in [`AGENTS.md`](AGENTS.md) §7); that root has not yet been copied to the
server's durable storage, and a reboot can clear local disk (STATUS handoff). Fetched inputs are archived on the
server's network storage, each file with a sha256 in [`sources.yaml`](sources.yaml); reference libraries are pinned in
[`refdb.yaml`](refdb.yaml). The repository holds code, plans, decisions, measurement records and tests.

**Stages** (stage numbering and plan: [`STAGES.md`](STAGES.md), "Status and mapping")

| Stage | Source | Size | Annotations | State |
|---|---|---|---|---|
| S0 golden set | 15 reference genomes ([`golden.yaml`](golden.yaml)) | | test harness only | used by every test ring; release `rS0.1` not cut (card E-06) |
| `S1` | INPHARED phages, 14Apr2025 | 34,076 records, 27,258 distinct sequences (*measured*, `sources.yaml`) | genes, Pfam 38.2, clusters | built |
| `P1`, `P2` | GTDB R226, best complete genome per species: 250 and 2,247 genomes | 10.6 Gbp with `S1` (*measured*, `STAGES.md`) | genes, Pfam 38.2, clusters | built |
| `P3` to `P6` | the same selection, the other 9,469 genomes | 39.9 Gbp (*measured*, `STAGES.md`); the selection is 11,973 genomes, 11,966 fetched, 50.5 Gbp (*measured*, [`2026-11-11_r50_inputs.md`](planning/measurements/2026-11-11_r50_inputs.md)) | genes, clusters; Pfam: `P3` and `P4` done, `P5` and `P6` running | building |
| `H1` (planned) | Huge Phage Genome Collection | 7,230 genomes, 1,961,987,211 bp, 516 already in `S1` (*measured*, `sources.yaml`) | none | fetched; card C-24 |
| plasmids (planned) | PLSDB 2024_05_31_v2 | 72,556 records, 7,297,669,596 bp, 7,139 already contigs of `P1` to `P6` (*measured*, `sources.yaml`) | none | fetched; card C-26 |
| S2C phase 2 (planned) | GTDB R226, up to 20 complete genomes per species, chosen for strain diversity with skani in the 243 species that have more | 13,485 genomes selected, 56.2 Gbp; 13,477 fetched, 56.18 Gbp (*measured* and *derived*, [`2026-11-13_s2c_skani_selection.md`](planning/measurements/2026-11-13_s2c_skani_selection.md), [`2026-11-13_s2c_phase2_fetch.md`](planning/measurements/2026-11-13_s2c_phase2_fetch.md)) | none | fetched and archived; stage split and build to come; card C-25 |

**Releases** (names are provisional: none is signed, so none is `rS<k>.<n>`, `STAGES.md`)

| Release | Stages | Contents | State |
|---|---|---|---|
| `rN.1` | `S1`, `P1`, `P2` | 12,048,866 genes, 10,620,684 distinct proteins, 2,966,200 NR50 and 8,278,537 NR90 clusters; no Pfam (*measured*, 2026-11-10, STATUS "10 Gbp rebuild") | verified; built before the lookup tables |
| `rN.2` | the same | the same with Pfam 38.2 | verified |
| `rN.3` | the same | the same with Pfam 38.2 and the lookup tables, built by current code; clusters equal to `rN.2` row for row (*measured*, `2026-11-11_build_r50.md` §A3) | verified; has a search pack: **the one to use** |
| `rN.1a`, `rN.1b` | intermediate releases of the same rebuild | | verified |
| `rM.1` | `S1`, `P1` to `P6` | 54,265 contigs, 39,224 genomes, 45,815,843 proteins, 31,183,876 NR90 and 8,192,065 NR50 clusters (*measured*, 2026-11-12, STATUS handoff and card B-46) | clustered; its later release steps ran clean in scratch (card B-46) and are to be rerun in the root |
| `rM.2` | `S1`, `P1` to `P6` | Pfam on every stage | not built: follows the scans |

Non-coding RNA: the producer and its frozen model set exist (DEC-28, DEC-30), but no built stage has been scanned, so no
release has ncRNA features.

## How it is built

Each stage is built by one command, `bio-box build-stage` (RUNBOOK §2), which runs every step below as a
`bio-box run-step` under a memory cap, a thread cap and a 300 CPU-minute budget per process, and writes a run record
and ledger rows for each step (AGENTS §4.5, §7).

1. **Sequence store** (`src/biobox/l0seq/`): every contig once, as bgzip FASTA shards named by content digest, with
   `.fai` and `.gzi`; one shard per ingest batch; a locator per stage (DEC-17, DEC-21).
2. **Gene calls and intergenic regions** (`src/biobox/l1/genecall.py`, `igr.py`): pyrodigal for prokaryotes,
   pyrodigal-gv for viruses, genes across a circular origin called once (DEC-13); a contig shared by several genomes
   is called once (DEC-26, DEC-27).
3. **Domains** (`src/biobox/l1/domains.py`): Pfam 38.2 at its gathering cutoffs with a fixed Z (DEC-04), in shards.
4. **Non-coding RNA** (`src/biobox/l1/ncrna.py`): Infernal on intergenic regions with a curated Rfam 15.1 set,
   `prok-2`, and clan competition (DEC-28, DEC-30).
5. **Clusters** (`src/biobox/l2/cluster.py`): MMseqs2 at 90% and 50% identity, every release from scratch (DEC-24),
   in parts under the CPU budget; a call no number of parts fits runs as one process under its own cap (DEC-29).
6. **Tapes** (`src/biobox/l3/tapes.py`): per contig, its features in order, keyed by `feature_uid` (DEC-18), in row
   groups bounded by slots (DEC-22).
7. **Tokens and postings** (`src/biobox/l3/tokens.py`, `postings.py`): release-local token ids (DEC-23) and one
   posting per slot and token; neighbour questions are answered at query time, with no table of pairs (DEC-19).
8. **Dimension tables** (`src/biobox/dim.py`): `contig_dim` and `genome_dim` hold what can change without the
   sequence changing (DEC-20); a membership view records every genome that holds a contig (DEC-27).
9. **Manifest** (`src/biobox/catalog.py`): the explicit file list of a release with hashes; readers use it, never
   globs (AGENTS §3.9).
10. **Catalog** (`src/biobox/catalog.py`): a read-only DuckDB database generated from the manifest, with lookup macros.
11. **Search pack** (`src/biobox/search/pack.py`): a release's proteins, intergenic regions and dark ORFs laid out for
    searching, kept outside the data root.

Every table is declared once in `src/biobox/registry/schemas.py` (AGENTS §3.7).

| Property | Rule | Defined in |
|---|---|---|
| Content-addressed identifiers | contigs, proteins and features are named by BLAKE3-128 digests; a gene's position in a list is never a key | AGENTS §3.2, §3.3; DEC-18; `src/biobox/ids.py` |
| A release is a function of its stages | a release names its stages; it reclusters from scratch; membership depends on the set of stages, not their order | DEC-24, DEC-27; card B-25 |
| Byte-identical rebuilds | same inputs, same bytes, on any host, thread count or input order; the pod and GitHub Actions give the same fixture manifest hash (*measured* at commit `132230f`, card D-04). Exception: the generated catalog file differs in about 65 bytes of DuckDB internals between builds (*measured*, card B-34); its wording is an open item for the PI | AGENTS §4.3; `tests/r3/test_cross_host.py` |
| Bounded memory | a builder holds one genome or one batch at a time and writes through one streaming writer | AGENTS §4.1, §4.2 |
| Every lookup prunes | a lookup by key reads O(1) row groups, also for keys from a joined table | AGENTS §4.4 |
| No rewrites | a stage never rewrites an earlier stage's files; index tables are partitioned by stage; files are renamed into place | AGENTS §3.9, §3.10; DEC-21 |

The specification and plan behind these are [`00-spec.md`](00-spec.md), [`PLAN.md`](PLAN.md) and
[`DEPLOYMENT-AND-MALLEABILITY.md`](DEPLOYMENT-AND-MALLEABILITY.md); they hold superseded text, and
[`planning/DECISIONS.md`](planning/DECISIONS.md) outranks them.

## Milestones

**Done**

| Date | Milestone | Record |
|---|---|---|
| 2026-10-28 | Specification, plan and stage ladder written after two rounds of independent audits; preflight gaps listed | [`audits/SYNTHESIS_2.md`](audits/SYNTHESIS_2.md), [`planning/preflight_gaps.md`](planning/preflight_gaps.md), `STAGES.md` |
| 2026-10-29 | Preflight decisions DEC-01 to DEC-16 (S2C scope, Z constants, gene-calling policy, repository and CI, server footprint); repository started | `planning/DECISIONS.md` |
| 2026-10-30 | First S0 gate report on the golden set; sequence-store format chosen by bake-off (DEC-17) | [`2026-10-30_s0_gate.md`](planning/measurements/2026-10-30_s0_gate.md), [`2026-10-30_k13_bakeoff.md`](planning/measurements/2026-10-30_k13_bakeoff.md) |
| 2026-11-09 | Design stress test; data model moved to coordinate keys and postings (DEC-18 to DEC-23); task cards and agent rules written; first card merged (A-01) | [`2026-11-09_design_stress_test.md`](planning/measurements/2026-11-09_design_stress_test.md), `planning/tasks/TASKS.md`, `AGENTS.md`, STATUS merged table |
| 2026-11-10 | Import ladder with the code of the time, `rL.3` (10,620,684 proteins, *measured*); current code cannot open it | [`2026-11-10_import_ladder/rL3.md`](planning/measurements/2026-11-10_import_ladder/rL3.md) |
| 2026-11-10 | 10 Gbp rebuild with current code: `rN.1`, equal to `rL.3` in every gene, protein, IGR and cluster row (*measured*) | STATUS, "10 Gbp rebuild" |
| 2026-11-11 | First end-to-end custom-HMM search with neighbourhoods, on `rN.1`: the hand-made answer exactly | `2026-11-11_trial_search_rN1_second.md` |
| 2026-11-11 | `rN.3` built by current code, clustered in parts, searched with Pfam names on the neighbours | `2026-11-11_build_r50.md` |
| 2026-11-11 | 50 Gbp inputs fetched and its build started; `main` takes 87 cards (pull requests 1 and 2) | `2026-11-11_r50_inputs.md`; STATUS, "State at 20:45 CST, 2026-11-11" |
| 2026-11-11 to 12 | Shared contigs: a contig belongs to the stage that stored it first and may belong to several genomes (DEC-26, DEC-27; cards C-20, C-21a to C-21c) | `planning/DECISIONS.md`; STATUS rows |
| 2026-11-11 to 12 | Curated ncRNA set: `prok-1` (588 models, DEC-28), then `prok-2` (770 models, DEC-30); producer with clan competition (cards E-07 to E-09) | [`2026-11-12_rfam_sample.md`](planning/measurements/2026-11-12_rfam_sample.md), [`2026-11-12_ncrna_prok2.md`](planning/measurements/2026-11-12_ncrna_prok2.md) |
| 2026-11-12 | `rM.1` clustered (45,815,843 proteins, *measured*) | STATUS handoff; card B-46 |
| 2026-11-12 | Split clustering equals the published clusters of `rN.1b` (card B-42) | [`2026-11-12_split_clustering_rN1b.md`](planning/measurements/2026-11-12_split_clustering_rN1b.md) |
| 2026-11-12 | Locus-context comparison: `nbr flanks` and its synteny figure (cards W-25, W-27) | `2026-11-12_flanks_rN3.md` |
| 2026-11-12 to 13 | Synthetic-world generator, sessions 1 to 7: backbone, taxonomy and plasmids, families and coding sequences, intergenic regions and RNA motifs, phages and CRISPR arrays, fusions, call errors and genetic codes 4 and 25, assembly classes, contamination and duplicates (card F-02) | STATUS rows F-02 (1) to (7) |
| 2026-11-12 | S2C phase-2 candidates fetched: 28,457 of 28,483 genomes (*measured*) | `2026-11-13_s2c_skani_fetch.md` |
| 2026-11-13 | Image built from `main` (pull request 4): 1,720 tests passed, 22 skipped inside it (*measured* on GitHub Actions, run 43380434543) | STATUS handoff |
| 2026-11-13 | S2C phase 2 selected for strain diversity with skani: near-clone pairs among a species' 20 fall from 6,930 to 2,895 of 46,170 (*measured*) | [`2026-11-13_s2c_skani_selection.md`](planning/measurements/2026-11-13_s2c_skani_selection.md) |
| 2026-11-13 | Second worked example and release check: CRISPR-Cas systems carried by Tn7-like transposons (Peters et al. 2017) on `rN.3`; all 19 invariants hold (*measured*) | [`examples/crispr_cas_neighbourhoods/README.md`](examples/crispr_cas_neighbourhoods/README.md) |
| 2026-11-13 | Worked examples are release checks: one command runs both on a release; on `rN.3` 92 checks ok, none failed (*measured*) | [`examples/README.md`](examples/README.md), [`2026-11-13_release_examples_rN3.md`](planning/measurements/2026-11-13_release_examples_rN3.md) |

**Now** (2026-11-13; STATUS handoff)

- The Pfam scans of `P5` and `P6`; `P3`'s finished on 2026-11-12 and `P4`'s on 2026-11-13. Then `rM.2`. The scans are
  expected to end about 2026-11-14 (*estimate*, STATUS).
- Tests on the server run in the image built from `main` since 2026-11-13 (card C-29); the 50 Gbp build keeps the
  older image until it is done.
- S2C phase 2 is fetched (13,477 of 13,485 genomes, 2026-11-13); its split into stages is an open item for the PI,
  and its build follows the 50 Gbp release (card C-25).
- Worked examples are run on a release before it is announced (AGENTS §5.7; `scripts/release_examples.sh`): the next
  release to check is `rM.2`.
- The synthetic-world generator's last step (emit), and the first build of a synthetic world by the normal builders,
  compared with its truth (card F-02, session 8).
- CRISPR arrays as features: MinCED approved by the PI on 2026-11-13 (card E-12); it needs the next image.
- Merged on 2026-11-13, not yet used by a build: the night thread allowance read from other users' CPU (card C-28).

**Next, in order**

1. `rM.2`, the 50 Gbp release with Pfam on every stage; verify it, copy the data root to durable storage, and finish
   `rM.1`'s release steps (STATUS handoff).
2. Huge phages (stage `H1`, card C-24), then PLSDB plasmids (card C-26).
3. S2C phase 2, up to 20 genomes per species: the selection, the fetch of the genomes kept (card C-25), and Pfam by
   cluster representatives with the rest inferred and labelled (DEC-32; design B-47). Whether gene calls become a
   function of the contig alone is decided before this stage is built (card C-21f).
4. ncRNA annotation of the built stages with `prok-2`; a dedicated tRNA caller is a later card (E-10).
5. The statistics phase: the `pair_stats` design, the synthetic-world sessions that remain, then calibration on null
   worlds (cards F-01 to F-03). No enrichment statistic exists today.
6. The first signed release, `rS0.1` (card E-06), after which locked schema columns can only be appended to (AGENTS
   §3.8).
7. Work beyond the shared server: scan shards on a larger system (card C-18, deferred), contamination screening with
   FCS-GX on a host with 512 GB of memory or more (DEC-02, deferred), and the rented hosts the plan names for its later
   stages (`STAGES.md` §4, §5; not started). Clustering a release far beyond the 50 Gbp set also needs an MMseqs2
   without the index abort (card B-44, deferred).

Several of these wait on decisions; they are listed under "Open items for the PI" in
[`planning/tasks/STATUS.md`](planning/tasks/STATUS.md).

## How the work is done

Work is split into task cards ([`planning/tasks/TASKS.md`](planning/tasks/TASKS.md)), each naming the only files it
may change, the commands that must pass and when to stop. Cards tagged **H** are design questions or decisions for
the PI, recorded in [`planning/DECISIONS.md`](planning/DECISIONS.md), which outranks every plan document. LLM agents
take one card each, on a branch of its own, under the rules in [`AGENTS.md`](AGENTS.md): write the test first and see
it fail, check against an independent oracle (an external tool, a separately written brute-force loop or a published
table, never the function's own earlier output), run tests on the server in the pinned image, and label every number.
An orchestrating session audits each report, merges it into `plan/task-cards`, reruns the affected tests on the merged
tree, checks GitHub Actions and records the card in [`planning/tasks/STATUS.md`](planning/tasks/STATUS.md). `main`
takes `plan/task-cards` by pull request, and each push to `main` builds the container images.

## Repository map

| Path | What |
|---|---|
| `README.md` | this file |
| `AGENTS.md`, `CLAUDE.md` | rules for agents that change the code (`CLAUDE.md` points to `AGENTS.md`) |
| `GUIDE.md`, `AGENT_GUIDE.md` | how a scientist, and an agent, search a release |
| `RUNBOOK.md` | operator commands: build, pause, resume, verify, search, clean up, limits |
| `STAGES.md` | the stage ladder and how what exists maps to it |
| `00-spec.md`, `PLAN.md`, `DEPLOYMENT-AND-MALLEABILITY.md` | specification, pipeline and schema plan, hosts and images (partly superseded) |
| `PRIOR-ART.md` | prior art |
| `sources.yaml`, `refdb.yaml` | approved data sources with what was fetched; pinned reference libraries |
| `golden.yaml`, `golden.lock.yaml` | the golden genome set and its lock |
| `pyproject.toml`, `uv.lock` | Python package and locked dependencies |
| `container/` | `Dockerfile` (runtime and test images) and `environment.yml` (native tools) |
| `profiles/` | host profile of the shared server: thread, memory and file-size limits |
| `src/biobox/` | the package: layers `l0seq`, `l1`, `l2`, `l3`, `search`, `synthworld`, the registry, catalog, CLI |
| `tests/` | rings `r0` (static), `r1` (unit), `r2` (fixture builds), `r3` (golden set); `registry`, `oracles`, `fixtures` |
| `scripts/` | test and run launchers for the server, schema-lock update, stress and build drivers, orchestration helpers |
| `planning/` | decisions, task cards and status, measurements, design notes, search and test plans, sizing |
| `audits/` | the audit rounds behind the plan |
| `presentations/` | the group-meeting deck |
| `.github/workflows/` | `ci`, `cross-host`, `image`, `mutation` |

## Requirements and tests

The package needs Python 3.11 or later and [uv](https://docs.astral.sh/uv/); dependencies are locked in `uv.lock`. The
native tools (MMseqs2, HMMER, Infernal, samtools, skani, and Prodigal as a test oracle only) come from
`container/environment.yml` and are used only inside the image. GitHub Actions builds two images on each push to
`main`, `ghcr.io/ifinkelstein/bio-box` and `ghcr.io/ifinkelstein/bio-box-test`, tagged by commit; the server pulls them
with `singularity` by digest (`container/Dockerfile`).

What CI runs (`.github/workflows/ci.yml`); tests that need MMseqs2, HMMER, Infernal or samtools skip outside the image:

```bash
uv sync --python 3.14
uv run pytest -q
```

Inside the test image (`.github/workflows/image.yml`):

```bash
docker run --rm --init ghcr.io/ifinkelstein/bio-box-test:<commit> pytest -q
```

On the shared server, with an account there (RUNBOOK §8): `sh scripts/pod_pytest.sh <name> [test paths]` mirrors the
working tree and runs pytest in the image; `--shards 8` splits the suite, and `--golden` sets `BIOBOX_GOLDEN_DIR` so
ring R3 runs. The project's agents run tests only there (AGENTS §5.5). The cross-host check builds the fixture on
GitHub Actions (`.github/workflows/cross_host.yml`) and compares its manifest hash with the server's
(`tests/r3/test_cross_host.py`).
