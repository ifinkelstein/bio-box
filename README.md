# bio-box (all-the-data)

**Documentation map:** **[`DOCS.md`](DOCS.md)** indexes every user guide, operator runbook, catalog register, design note, measurement report, and archived spec in one place.

**Status, 2026-11-25:** the 119.90 Gbp rolling release `r26.11.3` (thirteen stages `S1`, `P1` to `P6`, `H1`, `PL1`, `C5`,
`V1`, `G1` `annotated` Pfam 38.2, and `SPX1` `called`; *measured*,
[`2026-11-24_release_r26.11.3.md`](planning/measurements/2026-11-24_release_r26.11.3.md)) is **the live release to
search** (`--release latest` resolves to `r26.11.3`), with `r26.11.2` (74.16 Gbp) live beside it (`DEC-66` `S-4`;
pre-rolling releases `rN.1a`..`rM.2` are durably archived on `/stor/work` and `tower` via `bio-box release retire`,
card C-54). `r26.11.3` holds `3,121,518` genomes, `4,338,167` contigs, `119,876,311` genes, `102,337,875` distinct
proteins, `71,208,987` NR90 and `22,155,373` NR50 clusters, and verified search pack `2c94c11d…` (`89.01 GB`, *measured*).
No release has been signed yet. The running record is [`planning/tasks/STATUS.md`](planning/tasks/STATUS.md); where this
file and that one differ, that one is newer.

Numbers here carry the labels [`AGENTS.md`](AGENTS.md) §6 requires: *measured* (on the pod unless stated, with the
date and the record cited beside it), *derived* (with its formula) or *estimate*.

## What bio-box is

bio-box is a staged, content-addressed and reproducible database of microbial and phage genomes: their sequences,
gene calls, proteins, Pfam domains, non-coding RNA, protein clusters and gene order. It is built to answer discovery
questions by genomic context: search the whole collection for a protein family with your own profile HMM, then see
what sits next to every hit, which neighbours recur, and which loci share a context.

**The goal is a complete catalog of all useful sequenced data** (the PI, 2026-11-20, DEC-59). Sources are taken in
order of the new protein families they add per Gbp, because storage and CPU are short; that is the order of work,
not the scope, and whatever a stage leaves out (non-representative genomes, short contigs, a habitat stopped early)
is recorded for a later pass in [`planning/CATALOG_GOAL.md`](planning/CATALOG_GOAL.md).

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
| `bio-box results` | `archive`, `flush`, `list` and `get` keep a finished run with its outputs under the profile's `bulk` area, `<bulk>/results` (DEC-33, DEC-42; GUIDE §5). The first bundle, the RNA polymerase run of `rN.3`, is on `/stor` since 2026-11-14 |

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
university server (the pod, `wilkcomp02` in [`AGENTS.md`](AGENTS.md) §7), with verified `.tar` backups of the stages
and releases through `r26.11.3` on the server's durable storage (`/stor/work/Wilke/if2453/all-the-data/root_backups/`,
mirrored to the PI's NAS `tower`, STATUS handoff). Fetched inputs are archived on the server's network storage, each file with a sha256 in
[`sources.yaml`](sources.yaml); reference libraries are pinned in [`refdb.yaml`](refdb.yaml). The repository holds
code, plans, decisions, measurement records and tests.

**Stages** (stage numbering and plan: [`STAGES.md`](docs/archive/STAGES.md), "Status and mapping")

| Stage | Source | Size | Annotations | State |
|---|---|---|---|---|
| S0 golden set | 15 reference genomes ([`golden.yaml`](golden.yaml)) | | test harness only | used by every test ring; release `rS0.1` not cut (card E-06) |
| `S1` | INPHARED phages, 14Apr2025 | 34,076 records, 27,258 distinct sequences (*measured*, `sources.yaml`) | genes, Pfam 38.2, both CRISPR producers, labels, clusters | built |
| `P1`, `P2` | GTDB R226, best complete genome per species: 250 and 2,247 genomes | 10.6 Gbp with `S1` (*measured*, `STAGES.md`) | genes, Pfam 38.2, both CRISPR producers, labels, clusters | built |
| `P3` to `P6` | the same selection, the other 9,469 genomes | 39.9 Gbp (*measured*, `STAGES.md`); the selection is 11,973 genomes, 11,966 fetched, 50.5 Gbp (*measured*, [`2026-11-11_r50_inputs.md`](planning/measurements/2026-11-11_r50_inputs.md)) | genes, Pfam 38.2, both CRISPR producers, labels, clusters | built |
| `H1` | Huge Phage Genome Collection | 7,228 genomes, 6,714 contigs and 1,819,831,338 bp stored (514 shared with `S1`; *measured*, [`2026-11-15_stage_H1.md`](planning/measurements/2026-11-15_stage_H1.md)) | genes, Pfam 38.2, both CRISPR producers, labels, clusters | built (card C-24) |
| `PL1` | PLSDB 2024_05_31_v2 | 68,693 genomes, 61,647 contigs and 5,924,552,875 bp stored (7,046 shared with `S1` and `P1` to `P6`; *measured*, [`2026-11-15_stage_PL1.md`](planning/measurements/2026-11-15_stage_PL1.md)) | genes, Pfam 38.2, both CRISPR producers, labels, clusters | built (card C-26) |
| `C5`, `C1` to `C4`, `C6` (S2C phase 2) | GTDB R226, up to 20 complete genomes per species, chosen for strain diversity with skani in the 243 species that have more | 13,477 genomes fetched, 56,177,990,736 bp (*measured*, [`2026-11-15_s2c_phase2_inputs.md`](planning/measurements/2026-11-15_s2c_phase2_inputs.md)); `C5` stores 2,653 contigs and 9,109,234,757 bp across 2,029 genomes (2,550 shared with earlier stages; *measured*, [`2026-11-20_stage_C5.md`](planning/measurements/2026-11-20_stage_C5.md)) | `C5`: genes, Pfam 38.2, both CRISPR producers, labels, clusters | `C5` built (card C-47); `C1` to `C4` and `C6` inputs archived on `/stor` and `tower`, deferred after `r26.11.4` in favour of long-read metagenomes and high-diversity environmental waves (DEC-70) |
| `V1` | UHGV v1.0 gut viral vOTU representatives | 168,536 genomes, 168,358 contigs and 5,178,493,077 bp stored (178 shared with `H1` and `S1`; *measured*, [`2026-11-20_stage_V1.md`](planning/measurements/2026-11-20_stage_V1.md)) | genes, Pfam 38.2, both CRISPR producers, labels, clusters | built (card C-48c part 1) |
| `G1` | GlobDB r226 non-GTDB archaeal species representatives (`G2`..`G41` deferred in [`planning/CATALOG_GOAL.md`](planning/CATALOG_GOAL.md)) | 7,139 genomes, 1,227,214 contigs and 10,042,502,859 bp (10.04 Gbp) stored, 11,189,212 genes, 11,053,665 distinct proteins (*measured*, [`2026-11-23_stage_G1.md`](planning/measurements/2026-11-23_stage_G1.md)) | genes, Pfam 38.2 (15,056,794 rows), both CRISPR producers (1,773 + 127 arrays), labels | built & backed up (cutting in `r26.11.3`) |
| `SPX1` | SPIRE v1 unbinned aquatic & terrestrial metagenome contigs `>= 5,000` bp (tranche 1: 304 samples across 189 studies) | 2,829,368 genomes, 2,817,316 contigs and 35,688,628,875 bp (35.69 Gbp) stored (12,052 held by `G1`), 35,759,143 genes, 34,039,158 distinct proteins (33,512,765 new; *measured*, [`2026-11-23_stage_SPX1.md`](planning/measurements/2026-11-23_stage_SPX1.md)) | genes (`meta`), default CRISPR producer (5,333 arrays), labels (`--depth called`, DEC-54) | built & backed up (cutting in `r26.11.3`) |
| `V2`, `GV1`, `PL2`, `PL3` | IMG/VR v4.1 HC + Inovirus (`V2`), GVMAGs V2 (`GV1`), Yu2024 + 7 public PlasmidScope sets (`PL2`), IMG-PR + mMGE (`PL3`) | `V2` 8.20 Gbp new, `GV1` 3.89 Gbp new, `PL2` 1.98 Gbp new, `PL3` 16.28 Gbp new (*measured*, cards C-52, C-53, DEC-65) | `V2` `called`; `GV1`, `PL2` `annotated` (Pfam 38.2); `PL3` `called` | staged under DEC-65 & mirrored to `tower`; queued for `r26.11.4` (DEC-67) |

**Releases** (names are provisional: none is signed, so none is `rS<k>.<n>`, `STAGES.md`)

Release names (DEC-52, 2026-11-16). The releases below were named by hand: `rN.*` is the 10 Gbp set and `rM.*` the
50 Gbp set, the number counting builds of the same set. Every release from now on is named `r<YY>.<MM>.<n>`, the
n-th release cut in that month (`r26.11.1` is the first, over nine stages), and releases are rolling: one is cut
whenever a stage has been completed since the last one, over every complete stage. Older releases keep their
names and stay usable. `rS<k>.<n>` is kept for a signed release at a gate of the stage plan. Every command that reads
or searches a release takes `--release latest` for the newest one (card W-37): the greatest name by this scheme
(the older names sort before it, in the order they were built, `rM.2` last) that has a manifest and, for a search, a
verified search pack; the
command prints the name it took and records that name, and `bio-box releases` marks it in its column `latest`.

| Release | Stages | Contents | State |
|---|---|---|---|
| `rN.1` | `S1`, `P1`, `P2` | 12,048,866 genes, 10,620,684 distinct proteins, 2,966,200 NR50 and 8,278,537 NR90 clusters; no Pfam (*measured*, 2026-11-10, STATUS "10 Gbp rebuild") | verified; built before the lookup tables |
| `rN.2` | the same | the same with Pfam 38.2 | verified |
| `rN.3` | the same | the same with Pfam 38.2 and the lookup tables, built by current code; clusters equal to `rN.2` row for row (*measured*, `2026-11-11_build_r50.md` §A3) | verified; has a search pack; both worked examples hold (2026-11-13); the release to use for the 10 Gbp set |
| `rN.1a`, `rN.1b` | intermediate releases of the same rebuild | | verified |
| `rM.1` | `S1`, `P1` to `P6` | the same stages and proteins as `rM.2`; its cluster files are byte-identical to `rM.2`'s (*measured*, 2026-11-14, STATUS handoff) | clustered only, left so on purpose (it would be a second copy of `rM.2`; STATUS open items 1 and 6) |
| `rM.2` | `S1`, `P1` to `P6` | 54,265 contigs, 39,224 genomes, 45,815,843 distinct proteins, 31,183,876 NR90 and 8,192,065 NR50 clusters, Pfam 38.2 on every stage (*measured*, 2026-11-14, STATUS handoff; the cluster counts are `rM.1`'s, whose cluster files are the same bytes) | verified (2026-11-14); search pack `26e10f05…` (25.8 GB); both worked examples hold, 37 invariants (*measured*, 2026-11-15, [`2026-11-14_release_examples_rM2.md`](planning/measurements/2026-11-14_release_examples_rM2.md)); superseded as `--release latest` by `r26.11.1`, still available |
| `r26.11.1` | `S1`, `P1` to `P6`, `H1`, `PL1` | 59.88 Gbp stored across 122,626 contigs (7,046 shared across genomes) and 114,629 genomes; 57,086,886 genes, 49,102,041 distinct proteins, 33,117,773 NR90 and 9,291,312 NR50 clusters, Pfam 38.2, both CRISPR producers (16,067 default 3-repeat and 493 filtered 2-repeat near-Cas arrays) and `contig_label` on every stage (*measured*, 2026-11-18/19, [`2026-11-16_release_r26.11.1.md`](planning/measurements/2026-11-16_release_r26.11.1.md)) | verified (manifest hash `202b4cad…`); search pack `f8ab9bd5…` (29.29 GB); both worked examples hold, all 39 invariants (*measured*, 2026-11-19); superseded as `--release latest` by `r26.11.2`, still available |
| `r26.11.2` | `S1`, `P1` to `P6`, `H1`, `PL1`, `C5`, `V1` | 74.16 Gbp stored across 293,637 contigs (9,593 shared across genomes) and 285,014 genomes; 72,927,956 genes, 57,780,176 distinct proteins, 35,475,280 NR90 and 10,091,519 NR50 clusters, Pfam 38.2, both CRISPR producers (18,553 default 3-repeat and 589 filtered 2-repeat near-Cas arrays) and `contig_label` on every stage (*measured*, 2026-11-21/22, [`2026-11-21_release_r26.11.2.md`](planning/measurements/2026-11-21_release_r26.11.2.md)) | verified (manifest hash `46454d51…`); search pack `02dec831…` (36.01 GB); both worked examples hold, all 39 invariants (*measured*, 2026-11-22); superseded as `--release latest` by `r26.11.3`, still available |
| `r26.11.3` | `S1`, `P1` to `P6`, `H1`, `PL1`, `C5`, `V1`, `G1`, `SPX1` | 119.90 Gbp (*derived*) stored across 4,338,167 contigs (21,710 shared across genomes) and 3,121,518 genomes; 119,876,311 genes, 102,337,875 distinct proteins (44,557,699 new vs `r26.11.2`), 71,208,987 NR90 and 22,155,373 NR50 clusters, Pfam 38.2 on 12 `annotated` stages (`SPX1` `called`), 25,659 default 3-repeat and 716 filtered 2-repeat near-Cas CRISPR arrays, and `contig_label` on all 13 stages (*measured*, 2026-11-23/25, [`2026-11-24_release_r26.11.3.md`](planning/measurements/2026-11-24_release_r26.11.3.md)) | verified (manifest hash `956c237f…`); search pack `2c94c11d…` (89.01 GB); both worked examples hold, all 39 invariants including the CRISPR array checks (*measured*, 2026-11-25): **the one to use** (`--release latest`) |

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

The specification and plan behind these are [`00-spec.md`](docs/archive/00-spec.md), [`PLAN.md`](docs/archive/PLAN.md) and
[`DEPLOYMENT-AND-MALLEABILITY.md`](docs/archive/DEPLOYMENT-AND-MALLEABILITY.md); they hold superseded text, and
[`planning/DECISIONS.md`](planning/DECISIONS.md) outranks them.

## Milestones

**Done**

| Date | Milestone | Record |
|---|---|---|
| 2026-10-28 | Specification, plan and stage ladder written after two rounds of independent audits; preflight gaps listed | [`audits/SYNTHESIS_2.md`](docs/archive/audits/SYNTHESIS_2.md), [`planning/preflight_gaps.md`](planning/preflight_gaps.md), `STAGES.md` |
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
| 2026-11-13 | The synthetic-world generator runs end to end (card F-02, session 8), and a synthetic world was built by the normal builders for the first time: 99.0% of 93,883 true genes found by their stop, duplicates and shared contigs handled as designed (*measured*) | [`2026-11-13_synthworld_emit.md`](planning/measurements/2026-11-13_synthworld_emit.md) |
| 2026-11-13 | Worked examples are release checks: one command runs both on a release; on `rN.3` 92 checks ok, none failed (*measured*) | [`examples/README.md`](examples/README.md), [`2026-11-13_release_examples_rN3.md`](planning/measurements/2026-11-13_release_examples_rN3.md) |
| 2026-11-13 to 14 | CRISPR arrays: MinCED in the images, arrays as features (2,790 arrays on the 10 Gbp stages, *measured*), shown beside a locus by `bio-box nbr` (cards E-12, E-13, E-14) | [`2026-11-13_crispr_arrays.md`](planning/measurements/2026-11-13_crispr_arrays.md), [`2026-11-14_arrays_in_neighbourhoods.md`](planning/measurements/2026-11-14_arrays_in_neighbourhoods.md) |
| 2026-11-14 | Pfam 38.2 scanned on every stage of the 50 Gbp set (the last scan, `P6`, published 04:52 CST) | STATUS handoff |
| 2026-11-14 | A second CRISPR producer for the short arrays beside Cas genes (DEC-40, card E-15): MinCED at two repeats, kept only within 1 kb of a Cas gene, outside genes and away from the default arrays; 78 calls on 10.6 Gbp, and 14 of the 20 Tn7-like TniQ loci of the worked example now show an array (*measured*) | [`2026-11-14_crispr_two_repeat_producer.md`](planning/measurements/2026-11-14_crispr_two_repeat_producer.md) |
| 2026-11-14 | Search results are kept under the `bulk` area (DEC-42, card W-33); the first bundle, the RNA polymerase run of `rN.3` with its neighbourhood tables (164.6 MB, *measured*), is on durable storage | STATUS row W-33 |
| 2026-11-14 to 15 | 50 Gbp release `rM.2` (`S1`, `P1` to `P6`: 39,224 genomes, 45,815,843 proteins, *measured*), its sharded search pack (25.8 GB, card W-35), and both worked examples passing all 37 invariants (card W-36); CRISPR arrays and `contig_label` files added for all seven stages (cards C-33 part 1, C-34) | [`2026-11-14_pack_build_in_parts.md`](planning/measurements/2026-11-14_pack_build_in_parts.md), [`2026-11-14_release_examples_rM2.md`](planning/measurements/2026-11-14_release_examples_rM2.md), [`2026-11-15_label_files_old_stages.md`](planning/measurements/2026-11-15_label_files_old_stages.md) |
| 2026-11-15 to 16 | Stage `H1` (huge phages: 6,714 contigs, 1.82 Gbp stored, *measured*, card C-24) and stage `PL1` (PLSDB plasmids: 61,647 contigs, 5.92 Gbp stored, *measured*, card C-26) built with batched gene calling (cards B-48, B-50) and cross-stage Pfam domain copying (card B-47b); clustering reuse when proteins match (card B-49); S2C phase-2 inputs `C1` to `C6` prepared (card C-36) | [`2026-11-15_stage_H1.md`](planning/measurements/2026-11-15_stage_H1.md), [`2026-11-15_stage_PL1.md`](planning/measurements/2026-11-15_stage_PL1.md), [`2026-11-15_s2c_phase2_inputs.md`](planning/measurements/2026-11-15_s2c_phase2_inputs.md) |
| 2026-11-16 to 18 | Import depths (`stored`, `called`, `annotated`), `--min-contig-len` at ingest, shared-contig domain coverage, `domain_state` (`measured`, `inferred`, `none`, `not scanned`), and `--release latest` (DEC-52, DEC-54; cards C-38, W-37, C-38a, C-38b, C-38d, C-38h); inputs and overlap measured for `V1` (`UHGV`), `PL2` (`Yu2024` + `PlasmidScope` + `IMG/PR`), `GlobDB` `G1`..`G41` and the `SPIRE`/`MGnify`/`INSDC` run crosswalk (cards C-39 to C-41, C-43) | [`planning/design/C-38_import_depths.md`](planning/design/C-38_import_depths.md), [`2026-11-16_run_crosswalk.md`](planning/measurements/2026-11-16_run_crosswalk.md), [`2026-11-16_virus_votu_inputs.md`](planning/measurements/2026-11-16_virus_votu_inputs.md), [`2026-11-16_pl2_inputs.md`](planning/measurements/2026-11-16_pl2_inputs.md), [`2026-11-18_globdb_inputs.md`](planning/measurements/2026-11-18_globdb_inputs.md) |
| 2026-11-18 to 19 | First rolling release `r26.11.1` over `S1`, `P1` to `P6`, `H1` and `PL1` (59.88 Gbp, 114,629 genomes, 49,102,041 distinct proteins, 33,117,773 NR90 and 9,291,312 NR50 clusters, 16,067 + 493 CRISPR arrays, *measured*): verified, backed up onto `/stor/work`, search-packed (29.29 GB) and passed all 39 worked-example invariants (card C-35) | [`2026-11-16_release_r26.11.1.md`](planning/measurements/2026-11-16_release_r26.11.1.md) |
| 2026-11-19 to 22 | Stages `C5` (9.11 Gbp, 2,029 complete prokaryote genomes, card C-47) and `V1` (5.18 Gbp, 168,536 UHGV gut viral genomes, card C-48c part 1) built; representative-based Pfam domain inference evaluated against `C5`'s full scan and `r26.11.2` NR90 clusters (1.07% and 1.18% missed pairs, stopping at DEC-35 (6); DEC-62 chooses full scans for `C1` to `C4` and `C6`); second rolling release `r26.11.2` cut by `release_train.sh` over all 11 stages (74.16 Gbp, 285,014 genomes, 57,780,176 distinct proteins, 35,475,280 NR90 and 10,091,519 NR50 clusters, 18,553 + 589 CRISPR arrays, *measured*), verified, backed up onto `/stor/work`, mirrored to `tower`, search-packed (36.01 GB) and passed all 39 worked-example invariants (cards C-46, A-15 to A-20, B-54, B-61 to B-65) | [`2026-11-20_stage_C5.md`](planning/measurements/2026-11-20_stage_C5.md), [`2026-11-20_stage_V1.md`](planning/measurements/2026-11-20_stage_V1.md), [`2026-11-20_c5_pfam_by_representatives.md`](planning/measurements/2026-11-20_c5_pfam_by_representatives.md), [`2026-11-21_release_r26.11.2.md`](planning/measurements/2026-11-21_release_r26.11.2.md) |

| 2026-11-23 to 25 | Stages `G1` (10.04 Gbp GlobDB archaea, `annotated`) and `SPX1` (35.69 Gbp SPIRE unbinned `>= 5 kb`, `called`) cut into rolling release `r26.11.3` over 13 stages (119.90 Gbp *derived*, 3,121,518 genomes, 4,338,167 contigs, 119,876,311 genes, 102,337,875 distinct proteins, 71,208,987 NR90 and 22,155,373 NR50 clusters, 25,659 + 716 CRISPR arrays, *measured*), verified, backed up onto `/stor/work`, search-packed (89.01 GB) and passed all 39 worked-example invariants (cards B-68, B-70, W-41, A-25b, W-42); `DEC-77` hybrid clustering (`B-71`, `B-72`, `A-27`) merged for `r26.11.4` | [`2026-11-23_stage_G1.md`](planning/measurements/2026-11-23_stage_G1.md), [`2026-11-23_stage_SPX1.md`](planning/measurements/2026-11-23_stage_SPX1.md), [`2026-11-24_release_r26.11.3.md`](planning/measurements/2026-11-24_release_r26.11.3.md) |

**Now** (2026-11-25; STATUS handoff)

- `r26.11.3`, the 119.90 Gbp (*derived*) thirteen-stage rolling release (`S1`, `P1` to `P6`, `H1`, `PL1`, `C5`, `V1`,
  `G1`, and `SPX1`), is the release to use since 2026-11-25 (`--release latest` resolves to it). It holds 3,121,518
  genomes, 4,338,167 contigs (21,710 shared across genomes), 102,337,875 distinct proteins (44,557,699 not in
  `r26.11.2`), 71,208,987 NR90 and 22,155,373 NR50 clusters, Pfam 38.2 on 12 `annotated` stages (`SPX1` `called`),
  `contig_label` on all 13 stages, and 25,659 default 3-repeat plus 716 filtered 2-repeat near-Cas CRISPR arrays
  (*measured*, [`2026-11-24_release_r26.11.3.md`](planning/measurements/2026-11-24_release_r26.11.3.md)). Every unchanged
  cluster of `r26.11.2` kept its exact `cluster_uid` (34,062,365 of 34,062,365 at NR90 and 7,802,761 of 7,802,761 at
  NR50, *measured*). Its search pack (`2c94c11d5ddacf0e1f0dc8a2724ae1a6`, 89.01 GB) is verified, both worked examples
  pass all 39 invariants (91 of 155 TniQ+Cas7f loci carry a CRISPR array; 3,728 RNA polymerase hit proteins across
  4,392 genes, *measured*), and its additions (`43,744,563,200` B) are backed up onto
  `/stor/work/Wilke/if2453/all-the-data/root_backups/added_r26_11_3.tar` (`tower` mirror pending).
- External sequence fetching is paused (`STATUS.md` handoff): work focuses on ingesting the genomes, plasmids, viruses
  and metagenomes already staged on disk into the DuckDB/Parquet corpus through `stage_queue.sh` and `release_train.sh`
  (cards A-12, C-46).
- Tests, builds and searches on the server run in the images built from `main` at pull request 5 (cards C-29, C-32).
- The synthetic world's seed-stability and null-module positive-control adjustments wait on open item 12 for the PI
  (STATUS handoff), holding the statistics track (`F-01` to `F-03`).

**Next, in order**

1. Finish rolling release **`r26.11.3`** over the 13 built stages (`S1..V1`, `G1`, `SPX1`: `119.89 Gbp`, `119.88M` genes,
   `102.34M` distinct proteins, *measured*), then build Batch 2 (`V2` `IMG/VR` v4.1 HC + `Inovirus`, `GV1` `GVMAGs V2`,
   `PL2` redistributable plasmids, and `PL3` `IMG-PR` + `mMGE` -> **`r26.11.4`**, DEC-65, DEC-67).
2. Jump directly after `r26.11.4` to **public long-read (`PACBIO_SMRT` / `OXFORD_NANOPORE`) metagenome assemblies**
   (`LR1`.., card C-56) and high-diversity viral, satellite, integrative-MGE, single-cell `GORG` SAG, organellar, and
   `SPIRE` `>= 5 kb` metagenome waves (cards C-55, C-57, C-58, C-59, `SPX2`..), deferring strain-rich `C1` to `C4` and
   `C6` to the later strain-completion pass (`DEC-69`, `DEC-70`, [`planning/CATALOG_GOAL.md`](planning/CATALOG_GOAL.md)).
3. Fast multi-tool query indexes alongside the search pack (MMseqs2/Foldseek/pLM representative search) and on-hit
   annotation rounds (`domain_rounds`, cards C-38e to C-38g) on `--depth called` metagenomic stages.
4. ncRNA annotation of the built stages with `prok-2`, plus RNA viruses/Obelisks/viroid-like cccRNAs (card C-50) and
   microbial eukaryotes (card C-60).
5. The statistics phase: the `pair_stats` design, the synthetic-world sessions that remain, then calibration on null
   worlds (cards F-01 to F-03). No enrichment statistic exists today.
6. The first signed release, `rS0.1` (card E-06), after which locked schema columns can only be appended to (AGENTS
   §3.8).
7. CPU-heavy scans and larger-scale clustering on TACC Lonestar6 once allocation and SSH access are active (cards
   C-37a, C-37, C-38c).

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
| [`README.md`](README.md), [`DOCS.md`](DOCS.md) | overview (this file) and the repository-wide documentation map |
| [`GUIDE.md`](GUIDE.md), [`AGENT_GUIDE.md`](AGENT_GUIDE.md) | how a scientist, and an agent, search and query a release |
| [`RUNBOOK.md`](RUNBOOK.md) | operator commands on `wilkcomp02`: build, pause, resume, verify, retire, search, clean up, limits |
| [`AGENTS.md`](AGENTS.md), `CLAUDE.md` | rules for agents that change the code (`CLAUDE.md` points to `AGENTS.md`) |
| [`planning/`](planning/README.md) | [`DECISIONS.md`](planning/DECISIONS.md), [`CATALOG_GOAL.md`](planning/CATALOG_GOAL.md), [`tasks/`](planning/tasks/STATUS.md), [`design/`](planning/README.md#3-focused-design-notes-planningdesign), [`measurements/`](planning/measurements/README.md), [`search/`](planning/search/SEARCH_PLAN.md), [`tests/`](planning/tests/TEST_ARCHITECTURE.md) |
| [`docs/archive/`](docs/archive/README.md) | initial Oct–Nov 2026 specs (`00-spec.md`, `PLAN.md`, `STAGES.md`, `DEPLOYMENT-AND-MALLEABILITY.md`), `PRIOR-ART.md`, and `audits/` |
| [`examples/`](examples/README.md) | worked scientific searches (`rnap`, `crispr_cas_neighbourhoods`) run as release checks |
| [`presentations/`](presentations/README.md) | group-meeting HTML and PowerPoint decks |
| [`sources.yaml`](sources.yaml), [`refdb.yaml`](refdb.yaml) | approved data sources with checksums; pinned reference libraries |
| `golden.yaml`, `golden.lock.yaml` | the 15-genome golden test set and its lock |
| `pyproject.toml`, `uv.lock` | Python package and locked dependencies |
| `container/`, `profiles/` | `Dockerfile` + `environment.yml`; `wilkcomp02` host profile |
| `src/biobox/`, `tests/`, `scripts/` | package (`l0seq`, `l1`, `l2`, `l3`, `search`, `synthworld`, registry, catalog, CLI), test rings `r0`–`r3`, and pod/build drivers |
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
