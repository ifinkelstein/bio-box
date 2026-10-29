# Search plan: how ad hoc searches run on the bio-box layers

Status: **planning only** (2026-10-29). Nothing here is implemented, installed or run. This file is authoritative for search; the plan files (`00-spec.md`, `PLAN.md`, `STAGES.md`, `DEPLOYMENT-AND-MALLEABILITY.md`, `planning/data_layer_and_tests.md`, `planning/sources_and_sizing.md`) point here.

Inputs:

| file | what it holds |
|---|---|
| `00_brief.md` | meta answers, rewritten prompt, recipe schema, question-type codes |
| `01_synteny_workload.md`, `synteny_hypothesis_query_types.tsv` | all 365 synteny hypotheses classified by the search each needs; replay candidates; provenance contract |
| `02_recipes_2022_2026.md`, `recipes_2022_2026.tsv` | 177 recipes from 42 corpus papers |
| `03_recipes_2017_2021.md`, `recipes_2017_2021.tsv` | 197 recipes from 75 papers |
| `04_recipes_classic.md`, `recipes_1993_2016.tsv` | 206 recipes from 123 stems (1993–2016) |
| `05_external_modern.md`, `recipes_external.tsv` | 45 recipes from modern papers outside the corpus; tool throughput table; cost per query at 1e7–1e9 proteins |
| `06_infra_inventory.md` | what the plan already serves, holes and inconsistencies, sizes, envelope |

Numbers below are labeled: *measured* (by tool authors), *derived* (arithmetic shown in `05` or `06`), *estimate*. None is measured on wilkcomp02 yet.

---

## 1. What the evidence says (the eight findings the design rests on)

1. **Neighborhoods, not iterative search, are the commonest step.** 348/365 synteny hypotheses (95%) need a neighborhood step (`01` §6). Q-NBR is the most frequent code in every literature slice: 38 rows in 2022–2026, 34 in 2017–2021, 38 in 1993–2016 (`02` §2, `03` §2, `04` §2).
2. **Fixed-library annotation answers about a third of the workload.** 133/365 hypotheses are served by precomputed library hits plus neighborhoods. The rest need custom profiles (89), nucleotide-level work (120), IGR sequence (104), catalytic-residue checks (131, often for absence), or handling of promiscuous domains (143) (`01` §4).
3. **Searching outside the gene calls is routine.** About 85% of RODEO lasso precursors were not GenBank CDS; Kannan re-called 10 billion ORFs over 3.17 Tbp; tBLASTn-without-BLASTp hits are used to flag missed ORFs and pseudogenes; 15–20% of gene starts are wrong (`03` §1 summary, `02` §7.8, `04` §6.7). The raw contig store must be searchable directly, not only its protein table.
4. **Modern corpus-scale discovery reduces before it searches.** No 2023–2026 paper runs an iterative profile search over a whole corpus. They cluster, bracket a locus, subtract an expected partner or predict structure, then run the expensive tool on 10^3–10^6 objects. FLSHclust (8.8 Tbp, 8e9 proteins → 499.9M clusters at 30% → ~130,000 context-enriched clusters → 188 systems) is the canonical case, and its scoring step is the plan's tape + pair-statistics design (`05` §a Pattern A).
5. **Iteration stops by fiat, and curation between rounds uses locus context.** Fixed rounds (3–6, sometimes 10), inclusion 1e-3 to 1e-5 (0.01 for PSI-BLAST), manual vetting every round; only a few papers state a semantic stop rule ("stop when the last hit loses its locus context", "the last iteration before the first known false positive", ">98% retained, ≤10 new false positives") (`02` §4, `03` §4, `04` §5).
6. **E-values are rarely made comparable.** Only the Koonin lab pins database size (`-dbsize 20000000`, `2x10^6`), Price and Spang use 1e8, Fanzor pins `-Z`; thresholds are bitscores about as often as E-values (`02` §4, `03` §4, `04` §6.3). At 1e9 targets an inclusion E of 0.01 means a score about 1,000× stricter than against NR-2003 (`04` §6.3).
7. **At 1e9 proteins most search tools are bound by memory or I/O, not arithmetic** (`05` §c.9, *derived*): hmmsearch with one profile 5–10 min; jackhmmer ×5 ~50 min (12× less against 50%-identity clusters); default MMseqs2 index 2.1 TB (300 GB with `--prefilter-mode 1`); DIAMOND under 16 GB; Foldseek holds ~7e7 structures in 200 GB; translated search over all contigs needs a 42 TB index, over IGRs plus contig ends ~40 min. Batching queries makes MMseqs2 ~1,200× cheaper per query.
8. **Keep corpus-linear work local and candidate-linear work external.** Clustering, neighborhood extraction, enrichment statistics, first-pass homology and rule-based system calls are always local. HHpred, DALI, AF3/ColabFold, webFlaGs and EFI-EST/GNT are used per candidate on servers (17 papers in 2017–2021 used HHpred) (`05` §d, `03` §6, `02` §6).

---

## 2. Design decisions (SD-01..SD-21)

| ID | decision | why (evidence) |
|---|---|---|
| SD-01 | **Search is a recorded job, not a query.** Every search is a `search_run` launched from a declarative recipe (§4), run by the `bio-box search` runner, and written as packed Parquet with full provenance. `serve` never runs a search inline; it can submit and poll. | Minute-to-hour runtimes exceed any serve timeout; the synteny contract shows how to make results auditable (`01` §1.4; `06` §5) |
| SD-02 | **Every hit carries locus context.** Hit rows join to `protein_occurrence` (protein → every slot, contig, coordinates, strand, genome, taxonomy, species-rep flag, tier). Iteration rounds and curation see neighborhoods, not just sequences. | Curation between rounds is by tree plus neighborhood in every careful paper (`02` §10) |
| SD-03 | **Reduce first.** The default target is NR50 representatives; the final round runs on NR90 representatives and is projected to members through stored member alignments; members near the cutoff are re-scored exactly. The NR90/NR50 targets are built from ST1 on. Through S2C (≤5e7 exact-unique proteins, *estimate*) every run also searches all exact-unique proteins as its recall check; from S3 on, reduction is the default and exhaustive search is a sampled recall check. | `05` §c: NR50 cuts jackhmmer ~12×; projection needs alignments that the plan claims but does not store (`06` §3.B) |
| SD-04 | **One fixed search-space size per release line.** Protein searches use `Z_prot = 1e8` sequences and `domZ = 1e8`; nucleotide searches (Infernal, nhmmer) use `Z_nt = 1e5` Mb. Every hit stores the bitscore, the E-value at the fixed Z, and the real target size. Recipe thresholds are stored as bits; E-thresholds from papers are converted to bits at the paper's database size when it is stated. | Resolves the A-02 conflict (stored E-values change with a per-stage Z, `06` §3.C); matches the literature's 1e8 convention (`04` §6.3) |
| SD-05 | **Library scans: exact-unique proteins through S2C, NR90 representatives from S3 on, projected to members.** This replaces the three inconsistent statements (per genome, NR50, NR90). Members whose projected bitscore is within 10 bits of the cutoff, or whose alignment covers < 90% of the domain envelope, are re-scored exactly. Projected hits are release-scoped (`l1/domains_projected/release=<R>`); A-02 byte identity applies to exact rows, and projected rows are checked for call agreement under A-03. | `06` §3.C item 15; test S-08 checks projection |
| SD-06 | **Per-release search pack.** Each release builds, once, the search targets in §3.2 as ≥256 MB files on NFS with a manifest. Runs copy the needed shards to `/tmp`. | The plan generates FASTA on demand only; extracting all S5 proteins at the R-08 rate takes ~16 h (`06` §3.A, *derived*) |
| SD-07 | **Nucleotide search is first-class and tiered by data quality.** Three nucleotide corpora ship in the search pack: IGRs (with flank metadata), dark ORFs (stop-to-stop ORFs outside called genes, plus contig-end regions), and CRISPR spacers. Whole-contig translated search is a batch job reserved for tiers where it fits (§6). | Finding 3; `05` §c.8 |
| SD-08 | **Neighborhood fetch is an API with explicit semantics.** Window in genes or bp, strand and operon predicates, contig-edge policy (exclude, truncate, require N each side), circular topology, and censoring flags. Results report distance in genes and bp. | Papers differ on every one of these and state them as methods (`02` §7.3, `03` §7.1) |
| SD-09 | **User families can enter the statistics.** `bio-box search promote` turns a run's final hit set into a user token (namespace `source_db = user`) as a sparse overlay on the tapes. On-demand anchor statistics use the same score test against the release's background rates. Registering the model in `user_models.yaml` makes it a library model at the next release. | Today custom hits never reach `slot_pairs` or `pair_stats` (`06` §3.A item 4) |
| SD-10 | **Iteration is controlled, logged and bounded.** Per round: collapse hits to NR50 representatives; automated checks (coverage, required motif, foreign-domain alarm from stored library hits, locus-context retention, reverse search for new seeds); stop rule from the recipe; significance reported at the round of first inclusion. | `04` §6.1–6.2, `02` §4 |
| SD-11 | **Promiscuous domains are handled by data, not by hand.** A per-release `domain_promiscuity` table (number of distinct fusion partners and neighbor families per token, vehicle-deduplicated) feeds a blocklist for transitive seeding and a down-weight for fusion inference. | `04` §4.3 and §6.2; Prolinks hypergeometric correction |
| SD-12 | **Counts are dependence-aware.** Every search summary reports raw loci, vehicle-independent loci, species-representative count and distinct arrangements (Rousset), using the plan's existing vehicle clusters and BioProject caps. | `05` §e; `04` §6.4 |
| SD-13 | **Remote homology and structure stay external per candidate; local only where corpus-wide.** Image adds HH-suite3 with local Pfam-HHM and PDB70; hit-local HH databases are built only for hit clusters (10^2–10^4 HMMs). A corpus-wide 3Di database of NR50 representatives (ProstT5, GPU) is built at tier ST3. AF3/ColabFold, DALI and the Foldseek server are per-candidate external. | `05` §d; HH build cost per cluster is unpublished (`05` §c.5) |
| SD-14 | **Rule of thumb for local vs external:** local if cost is linear in the corpus and the output feeds statistics; external if cost is linear in candidates a human reviews. Crossover about 10^4 queries. | `05` §d |
| SD-15 | **An interactive lane on the pod.** By day the build uses 32 threads and the search lane 8; at night 48/8. Search jobs are split into per-process shards under 300 CPU-min each (watchdog limit ~369), at nice 15, with a `/tmp` quota of 150 GB. During L2 barriers the lane reads from NFS instead of `/tmp`. | One 40-thread process hits the watchdog in ~9 min (`06` §5) |
| SD-16 | **Search capability scales by tier, starting with closed, high-quality nucleotide data** (§6), with a new stage S2C (complete prokaryotic genomes) inserted before the draft and MAG stages. | User direction (2026-10-29); Koonin-lab practice uses 24,757 and 47,545 complete genomes (`02` §5) |
| SD-17 | **Domainator is the candidate workbench; bio-box does not build one.** Domainator (MIT, conda) ships in `bio-box:3.1`. Its similarity-network viewer (`.dsnv`), `plot_contigs`, structure annotation, profile comparison and browser GUI (`domainator_server`) serve candidate-scale review. `bio-box nbr --format domainator-gb` and `bio-box search hits --format domainator-gb` write GenBank with `Domainator` features and its qualifiers (`program`, `database`, `description`, `evalue`, `score`, `name`, `identity`, `cds_id`, `rstart`, `rend`, `rlen`), with `evalue` at Z_ref and `score` in bits. Domainator's own filters mix E-value and score and default to `Z = 1000`, so exported subsets are filtered on score. | Johnson et al. 2025 NAR (doi:10.1093/nar/gkae1175); repo `docs/file_formats.md`, `docs/limitations_and_FAQ.md` |
| SD-18 | **Agent-friendly, capped, schema-stable outputs.** Every `bio-box search`, `nbr` and report command has `--json` (NDJSON, one object per record, streamable). Lists are capped at query time (`--top_n`, default 50; member lists only on request), so one stored result serves any cap. Schema rule: added keys are compatible, removed or renamed keys are breaking and need a version bump (test S-13). `serve` endpoints and MCP tools wrap the same JSON. | Domainator `docs/agent_workflows.md` (`ssn_navigator`, `--json`, `--top_n`) |
| SD-19 | **Hit-set similarity networks at every threshold from one pass.** For a hit set, build the exact maximum spanning forest plus k nearest neighbours in a streaming pass, and store a threshold table (cluster count per distinct cut, merge events). Any cut is a lookup, so subfamily thresholds can be chosen and compared without recomputing. Edge weight is an EFI-style alignment score from DIAMOND or MMseqs2 bitscores, which does not depend on database size. Implementation calls Domainator's `build_ssn` / `matrix_report` code rather than reimplementing it. Network output is a `.dsnv` bundle plus a Parquet threshold table. | Domainator `docs/algorithms.md` (MSF ∪ kNN, threshold tables); fills the Q-SSN gap (`06` §2) |
| SD-20 | **Guardrails on every run.** `max_output_gb` (default 25 GB) per run, enforced by the runner, protects NFS and the small-file waiver. A `partial` field (`exclude`, `include`, `flag`) handles fragment genes and ties into the contig-edge flags in `protein_occurrence`. For nucleotide profile comparisons, default score cutoffs are scaled by the per-column maximum (about 1.39 bits for DNA/RNA against about 3.0 for protein), and the scaled value is recorded. | Domainator `--max_output_gb`, `--partial exclude`; `docs/limitations_and_FAQ.md` (nucleotide profile scores) |
| SD-21 | **3Di encoder bake-off before the corpus-wide structure build.** At ST3, ProstT5 and the ESM-2 3B 3Di encoder (Domainator's route) are compared on a fixed set of NR50 representatives for throughput (CPU and GPU), Foldseek recall against AFDB-derived 3Di, and SCOPe/ECOD benchmark sensitivity. reseek is tested as an alternative aligner. The winner builds `3di_nr50/`. | Domainator `docs/esm_3b_foldseek.md`; `05` §c.6 |

---

## 3. What the infrastructure gains

### 3.1 New tables

| table | layer | one row per | columns | size (S5, *estimate*) |
|---|---|---|---|---|
| `protein_occurrence` | L1 (derived per stage) | slot that carries a protein | `protein_uid`, `slot_uid`, `genome_uid`, `contig_uid`, start, end, strand, `gene_call_version`, `tax_order`, `is_species_rep`, `tier`, `assembly_quality`, `dist_to_contig_end_bp` | ~6e8 rows × ~48 B ≈ 29 GB |
| `cluster_member_aln` | L2 | member protein per clustering level | `member_uid`, `rep_uid`, `level` (nr90/nr50/subfam), qstart, qend, tstart, tend, pident, alnlen, compact CIGAR | ~19 GB per level (4.8e8 × ~40 B, upper bound); ~40–60 GB for all three levels |
| `domain_promiscuity` | L4 | token | distinct fusion partners, distinct neighbor families, vehicle-deduplicated counts, blocklist flag | small |
| `background_rates` | L4 | token × regime × clade × tier | per-window occurrence probability used as p0 for on-demand anchor statistics | small–medium |
| `search_run` | search | run | recipe hash, argv, tool versions, target pack and shard sha256, Z, query hash, outputs, wall, peak RSS, CPU-min per PID, host, exit code, status (`ok`, `failed`, `timeout`, `oom`, `stopped`; a `stopped` run keeps finished shards and resumes when rerun) | tiny |
| `search_hit` | search | hit per round | `run_id`, round, target kind (protein, segment, dark ORF, IGR, contig region, spacer), target UID, bitscore, E at Z_ref, target size, envelope, alignment coordinates on the target, alignment coordinates on the profile (`hmm_from`, `hmm_to`, `model_len`, as Domainator's `rstart`/`rend`/`rlen`), first-inclusion round, curation label and reason | per run |
| `search_profile` | search | profile per round | HMM/CM/PSSM blob (zstd), MSA (Stockholm, zstd), number of sequences, effective number | per run |
| `user_token_overlay` | L3 | slot hit by a promoted family | `slot_uid`, `user_token_id`, segment envelope, `run_id` | sparse |

Sidecar rule clarified: coverage means the fraction of the sidecar's keys that resolve in the release, not the fraction of the release's keys present in the sidecar. Sparse hit lists are therefore accepted (`06` §3.A item 5).

### 3.2 Search pack (built once per release; ≥256 MB shards on NFS)

| target | content | tools | S5 size (*estimate*) | S6 size (*estimate*) |
|---|---|---|---|---|
| `nr50.faa` shards + MMseqs2 DB (`--prefilter-mode 1` index) | NR50 representatives | pyhmmer, jackhmmer, MMseqs2 | 33–56 Gaa; ~20–35 GB zstd FASTA; ~35–60 GB index | 40–80 Gaa |
| `nr90.faa` shards | NR90 representatives | final round, projection | 80–110 Gaa | 115–200 Gaa |
| `unique.dmnd` | all exact-unique proteins | DIAMOND (first pass, recall check) | ~100–150 GB | ~330 GB; Hetzner only |
| `segments.faa` | domain envelopes, plus uncovered residue segments ≥50 aa (Pattern H) | envelope re-search | ~1.1× NR90 residues | — |
| `igr.fna` + `igr90.fna` | IGRs ≥30 nt (Weinberg convention) with flank genes, orientation, GC, distance to contig end. The base Rfam scan (PLAN §3.1) covers IGRs ≥40 bp; 30–39 nt IGRs are search-only | cmsearch, nhmmer, BLASTn, CMfinder | 48–54 Gbp; ~14 GB zstd | 324–372 Gbp |
| `dark_orf.faa` | stop-to-stop ORFs ≥30 aa in six frames outside called genes, plus all frames in the 3 kb at each contig end | MMseqs2, DIAMOND, pyhmmer | ~30 Gaa (*estimate*; measure at S0) | ~200 Gaa |
| `spacers.fna` + index | CRISPR spacers (union of MinCED and CRISPRCasFinder calls) with array, repeat and contig | protospacer search | small | small |
| `models/` | pressed library HMMs, Rfam CMs, user models | hmmscan, cmscan | small | small |
| `hh/` | Pfam-A HHM, PDB70 | HH-suite3 | ~30–60 GB (*estimate*) | same |
| `3di_nr50/` (from ST3) | 3Di of NR50 representatives (encoder chosen by the SD-21 bake-off: ProstT5 or ESM-2 3B) | Foldseek, reseek | ~50 GB without Cα | ~60–120 GB |

S5 search pack total on NFS: about 370–530 GB (*estimate*; the row sum in `planning/sources_and_sizing.md` "Search pack sizing"), against 30 TB of scratch. The NR50/NR90 rows use the 4.8e8-protein upper bound and run 15–25% high against the reconciled S5 range of 120–135 Gaa. Per-run working set in `/tmp`: 60–120 GB (NR50 plus IGRs).

### 3.3 New commands (planned; none exists)

```
bio-box search run <recipe.yaml>           # submit; returns run_id; async
bio-box search status|log|cancel <run_id>
bio-box search hits <run_id> [--round N] [--members] [--with-context ±N genes|bp]
                    [--json] [--top_n 50] [--format parquet|tsv|domainator-gb]   # SD-17, SD-18
bio-box search ssn <run_id> [--score efi] [--knn 10]  # MSF ∪ kNN + threshold table, .dsnv bundle (SD-19)
bio-box search promote <run_id> --name <family>      # user token overlay (SD-09)
bio-box search stats <run_id|user_token>             # on-demand anchor statistics against background_rates
bio-box nbr --anchors <run_id|token|uids> --window 10g|20kb --strand any|same --operon-gap 50
            --edge exclude:3kb|truncate|require:5 --circular auto     # neighborhood fetch (SD-08)
            [--json] [--format parquet|domainator-gb]
bio-box nbr score --positive-set <file|token list> [--exclude-anchor-family]   # Doron/Millman-style
bio-box region translate <contig_uid:start-end> [--six-frame] [--min-orf 30]
bio-box pack build|verify <release>                  # search pack (SD-06)
```

---

## 4. The recipe

A recipe is YAML, hashed, stored with the run, and re-runnable on a later release.

```yaml
recipe: rosetta_uri_dnaq_v1
written_against: rS2C.1            # release the recipe was written for
seed:
  - kind: protein_accession        # or hmm | msa | cm | structure | locus | region
    value: Rv2191                  # verified against the release's identity tables at submit
    source: {paper_id: 1999_Aravind_Walker_Koonin_DNA_Repair_Domains, evidence: "lit/...txt:NNN"}
split_domains: true                # search each envelope separately (fusions)
target: nr50                       # nr50 | nr90 | unique | segments | dark_orf | igr | contigs | spacers
tool: {name: jackhmmer, engine: pyhmmer, version: from_image}
iterations: {max: 5, inclusion_bits: 27, stop: [no_new_reps, context_lost, drift_alarm]}
z: release_default                 # Z_prot = 1e8
curation:
  min_query_coverage: 0.6
  required_motif: {columns: [D, E, D], source: M-CSA}     # optional
  foreign_domain_alarm: true       # stop if a new library family enters >X% of round-k hits
  promiscuity: down_weight         # SD-11
context: {window: 10g, strand: any, edge: require:5}
partial: flag                      # exclude | include | flag fragment genes (SD-20)
max_output_gb: 25                  # runner-enforced cap (SD-20)
controls:
  positive: [protein_accession: MJ0613]    # required; the run fails if the positive is not recovered
  negative: shuffled_seed
outputs: [hits, msa_per_round, hmm_per_round, neighborhoods, summary]
```

Contract, adopted from the synteny project (`01` §1.4) with its gaps closed:
- One runner launches every tool; non-zero exit raises; timeouts and OOM are `failed`, never "0 hits".
- Every hit row carries `run_id`; curated keep/drop labels are labels with a reason, never numbers.
- Seeds from papers carry `paper_id` and `file:line`; seeds from databases carry the database release and entry; accessions are resolved by sequence, not identifier (GI numbers are retired).
- Model IDs are resolved against the release's name table at submit time; unknown IDs are refused.
- A positive control is required; the negative defaults to `shuffled_seed` (MR-09); degenerate result columns (all zero, all identical) raise a warning.
- The recipe's controls are the controls its workflow manifest reuses (below), so the run and its tests check the same things. Where a synthetic world covers the recipe's operations, synthetic controls are preferred: a planted-truth positive and a null-world negative (TEST_ARCHITECTURE §6.2); real-data controls stay as a second check.
- No fallback provenance; inputs are identified by content hash.
- The release the recipe was written against and the one it ran against are both recorded.

**Sibling manifest.** Every recipe has a workflow manifest next to it (TEST_ARCHITECTURE §5.1; its example, `rosetta_fusion_census`, is the manifest for this recipe). The manifest lists the operations used (`O-SEARCH`, `O-ITER`, `O-PROJECT`, `O-NBR`, `O-ARCH`, `O-ANCHORSTATS`, `O-EXPORT` here), the registered output schemas, the controls (the recipe's own plus synthetic ones), planted or replay truth, and the §5 budget class. The harness generates the tests from it (TEST_ARCHITECTURE §5.2).

### 4.1 Iteration controller (SD-10)

Per round k:
1. Search the target with the current profile at the fixed Z; keep hits above `inclusion_bits`.
2. Collapse hits to NR50 representatives; record which are new in round k (first-inclusion round).
3. Checks: query coverage; required motif; foreign-domain alarm (a library family not in the seed appears in more than a set fraction of new hits); locus-context retention (fraction of new hits whose neighborhood shares at least one family with round-1 neighborhoods); reverse search of new representatives against the round-1 profile set.
4. Build the next MSA from hit envelopes only (`hmmalign` to the profile, trimmed to match states), so a fused partner domain never enters the alignment.
5. Log the round (hits, alarms, profile) as rows in `search_hit` and `search_profile`.
6. Stop on any configured rule; on a drift alarm, revert to round k−1.

Final round: search NR90 representatives, project to members through `cluster_member_aln`, re-score near-cutoff members exactly (SD-05).

### 4.2 Worked example: a Rosetta-stone check (SYN-143, Uri–DnaQ and Uri–EndoIII)

1. Recipe seeds Rv2191 (H37Rv, an S0 golden genome) with `split_domains: true`, giving GIY-YIG (PF01541), DEDD (PF00929) and HhH-GPD (PF00730) envelopes.
2. For families with library models, no search is needed: `slot_pairs` with `junction_topology = 0` gives fused counts, and topology counts give the split forms (Q-ARCH from stored data, seconds).
3. For divergent members, jackhmmer per envelope on NR50 (S5: ~1–2 min per profile per round, *derived* from `05` §c.1 at 1.7e8 sequences); final round on NR90 with projection.
4. `bio-box nbr` on all hits with `--window 10g --edge require:5`; `bio-box search stats` against background rates; counts reported raw, vehicle-independent, species-rep and by arrangement.
5. Classify each locus: fused (non-overlapping envelopes in one protein), split-adjacent, split-distant, single. Down-weight promiscuous partners; check fused calls against frameshift and read-through flags.
6. `bio-box seq get --slot ... --flank 10000 --nt` for loci to plot; candidates go to HHpred and AF3 externally.

---

## 5. Question types: what serves each, and cost

Costs are per query on the search lane (8 threads by day; 40-thread figures from `05` §c scaled ×5 where stated). S0 is seconds for everything.

| code | served by | S5 (*derived*) | S6 (*derived*) | notes |
|---|---|---|---|---|
| Q-FAM | stored library hits (L1) via `protein_occurrence`; new HMMs by `scan-new-models` on NR90 + projection | seconds (stored); ~0.5 core-h per new HMM | ~1 core-h per HMM | `scan-new-models` gains the projection and re-score rule |
| Q-ITER | recipe runner, iteration controller, NR50 target | jackhmmer 5 rounds ~5–10 min at 40 threads; ~30–50 min on the lane | ~10–20 min at 40 threads | MMseqs2 3-iteration profile for many seeds at once (batch lane) |
| Q-REMOTE | HH-suite3 vs Pfam-HHM and PDB70 locally; hit-local HH DB; HHpred external | seconds–minutes | same | no corpus-wide HH DB until a workload needs >10^4 profile queries |
| Q-STRUCT | external AF3/ColabFold/DALI per candidate; Foldseek vs AFDB on server; Domainator `structure_domainate`/`structure_search` for candidate sets; local 3Di of NR50 from ST3 | ~11 min per query vs 8e7 reps at 40 threads | same | 3Di build ~0.7–1.2k GPU-h at S5 NR50 with ProstT5 (*derived* from 41 proteins/s); encoder chosen by the SD-21 bake-off |
| Q-ARCH | `slot_pairs` (`junction_topology = 0`) for tokens; `segments.faa` re-search and user tokens otherwise | seconds (tokens) | seconds | promiscuity table applied |
| Q-NBR | `bio-box nbr` over tapes; `nbr score` for positive-set enrichment | seconds per anchor set (Q2 fetch); full-corpus positive-set score is a batch scan of the tapes (*estimate* minutes–1 h) | Hetzner for full scans | edge and circular semantics explicit |
| Q-PROF | presence bitvectors and pangenome planes; on-demand for user tokens | seconds | seconds | tree-weighted counts |
| Q-RNA | cmsearch/nhmmer on `igr.fna`; targeted CMfinder on context-anchored IGR subsets | ~12 min per CM over all IGRs at 40 threads | ~1.5 h | calibrate with Rfam-calibrated CMs or Infernal 1.1.5 `--split`; `cmcalibrate` needs ~14 GB per thread |
| Q-NT | `dark_orf.faa` protein search; `bio-box region translate`; spacer → protospacer streamed scan; contig BLASTn for small tiers | minutes (dark ORFs); ~40 min translated over IGRs + edges at S6 scale | DIAMOND blastx over all contigs ~55 h, Hetzner or spot | whole-contig tblastn not planned at S5+ |
| Q-MOTIF | per-hit alignment to profile columns from `search_hit`; curated sites from L1 | seconds | seconds | absence calls require dark-ORF and frameshift checks |
| Q-SSN | `bio-box search ssn`: all-vs-all inside a hit set (DIAMOND/MMseqs2), EFI-style score, MSF ∪ kNN with a threshold table (SD-19); `.dsnv` bundle for the Domainator viewer; Cytoscape export | minutes for ≤1e5 | same | EFI-EST stays external for single families |
| Q-CLUST | NR90/NR50/subfamily/IGR90 clusters; uncovered-segment clusters | stored | stored | `PRIOR-ART.md` novelty claim narrowed (Domainator 2024) |
| Q-SYS | DefenseFinder, PADLOC, CRISPRCasTyper, MacSyFinder sidecars | stored | stored | PADLOC-DB candidate systems excluded from the B1 temporal holdout |

---

## 6. Tiered implementation: closed, high-quality nucleotide data first

Search capabilities switch on tier by tier. A tier is certified when its replays (§7) pass on that tier's data and its cost table is measured on the pod. Data quality rises from closed replicons to drafts and MAGs; capabilities that only matter for fragmented assemblies (edge censoring, gene-call rescue for absence claims) become mandatory when fragmented data arrives.

| tier | data (quality) | stages | capabilities switched on | replays (acceptance) | scale |
|---|---|---|---|---|---|
| **ST0** | 14 closed golden genomes + spike-ins; first test bed is the micro-corpus and synthetic worlds with planted truth and null worlds (TEST_ARCHITECTURE §6.1–6.2, ring R2), then the golden set (R3) | S0 | recipe runner and provenance; `protein_occurrence`; `bio-box nbr`; exhaustive search on all proteins (no reduction); IGR, dark-ORF and spacer corpora; six-frame region translate; Z policy; `--json` outputs and the schema rule (SD-18); run guardrails (SD-20); `domainator-gb` export; cost measurement (cells/s on Haswell) | Domainator as an independent implementation of domain calls and neighborhoods on the golden genomes (S-02); the three Domainator case studies and five synteny hypotheses written as recipes (S-12); Makarova 2011 defense islands (golden subset); Wolf 2012 arCOG concordance; Weinberg 2010 motifs in *B. subtilis*, *B. theta*, SF370; Shmakov 2019 Set2 (*H. volcanii* seeds); Axen 2014 BMC loci (MG1655 eut/pdu); SYN-143 and SYN-125 end to end | ~4.5e4 proteins |
| **ST1** | complete phages and plasmids (closed, mostly circular) | S1, S2 | NR90/NR50 targets and search pack; iteration controller at scale; projection and re-score; spacer → protospacer; circular neighborhoods; first cross-tool concordance (pyhmmer, MMseqs2, DIAMOND) | Yutin 2024 (23,246 polymerases, 18,383 phage genomes); Peters 2017 spacers; Faure 2023 Tn7 (plasmid subset); Chaban 2024 catalytic-motif column; Silas 2025 (phage loci) | 2–8e6 proteins |
| **ST2** | complete prokaryotic chromosomes (assembly level Complete Genome, ≤20 per GTDB species, chosen by quality then diversity) | **S2C (new)** | user-token promotion and on-demand anchor statistics; `nbr score` with editable positive sets; `domain_promiscuity`; Rosetta and architecture queries at scale; local HH-suite (Pfam-HHM, PDB70); batch query lane; hit-set similarity networks with threshold tables (SD-19) | Domainator case 1, CARF uncovered segments (13,116-genome complete set, Zenodo 10.5281/zenodo.10989173); Makarova 2011 defense islands (full); Shmakov 2019 Sets 1–3; Shmakov 2018 CRISPR islands; Wu 2026 Druantia/Zorya (RefSeq loci); Bell 2026 YprA; Jimmy 2020 toxSAS; Atkinson 2011 RSH; Mutte all-models score matrix | est. 15–30k assemblies, 60–120 Gbp, 5.5–11e7 proteins (~2–4e7 exact-unique) |
| **ST3** | drafts and MAGs: small-genome reps and the remaining species reps | S3, S5 | edge censoring mandatory; dark-ORF and frameshift checks mandatory for absence claims; quality-weighted counts; 3Di encoder bake-off (SD-21), then 3Di of NR50 (GPU); targeted CMfinder; B1 temporal holdout | Domainator case 3, orphan methyltransferases to cryptic RM systems (~366 GB good-quality GenBank; S5 part, full at S6); Weinberg 2017 motifs; Makarova 2019 archaeal dark matter; Makarova 2024 Asgard ESCRT; Kaur 2022 death-domain conflict systems; Mestre 2020 RTs; Wu 2026 full; Corbino 2005 | ~4.8e8 proteins at S5 |
| **ST4** | metagenomic viruses (local-only) and all GTDB genomes | S4, S6 | search pack exceeds the pod: search host on Hetzner (≥512 GB RAM, ≥4 TB NVMe, *estimate*); NR90 projection everywhere; batch lane with query batching; DIAMOND blastx over contigs as a scheduled job | Domainator case 2, P4-like satellite hotspots in 2,718 *E. coli* genomes (needs S6; S2C's 20-per-species cap gives only a partial check); Domainator case 3 full; Gussow 2020 anti-CRISPR; Silas 2025 full; Yutin 2018 crAss (S4); a FLSHclust-style CRISPR-proximity cluster score on a subset | ~3e9 proteins, ~1e9 unique |
| **ST5** | AllTheBacteria, SPIRE, Logan | S7+ | compute-to-data; reduction-only searches (cluster, then score by context); Logan-Search and Serratus as external nucleotide services | Altae-Tran 2023 FLSHclust subset | 10^9–10^10 proteins |

S4 enters the ladder before S5 for data reasons (STAGES.md); its search acceptance waits for ST4. S3's small-genome reps are loaded at S3 for genetic-code work; their search acceptance is part of ST3.

### 6.1 New stage S2C: complete prokaryotic genomes

- **Adds:** GTDB R226 genomes (representatives and members) whose NCBI assembly level is `Complete Genome`, at most 20 per species (selected by CheckM2 completeness and contamination, then by ANI diversity). Complete species representatives enter here and are not re-ingested at S5; complete members enter here and are not re-ingested at S6.
- **Count:** from `bac120_metadata_r226.tsv` and `ar53_metadata_r226.tsv` (`ncbi_assembly_level`), an S−1 task. *Estimate* 15–30k assemblies, 60–120 Gbp, 3,000–22,000 core-h at the S5 rate (50–180 core-h/Gbp, *derived* from `planning/sources_and_sizing.md` S5 row).
- **Why here:** closed chromosomes give exact neighborhoods with no contig-edge censoring, so neighborhood statistics and replays are clean before drafts and MAGs arrive; this is the corpus Koonin-lab searches use (24,757 and 47,545 complete genomes).
- **Release names:** `rS2C.<n>`. Existing stage numbers are unchanged.
- **Gate:** the ST2 replays pass; cross-tool concordance reported; cost table measured.

---

## 7. Tests (family S, added to the stage gates)

| ID | test | op (TEST_ARCHITECTURE §2) | ring (§8) | gate |
|---|---|---|---|---|
| S-01 | Replay recall: for each tier's replays, the fraction of published hits (mapped by sequence) that lie in loaded genomes and are recovered. Target ≥ 0.90; every miss is classified (threshold, gene call, absent from the build, error in the published list). Precision is reported, not gated. | O-SEARCH, O-ITER, O-PROJECT, O-NBR (whole workflow) | R5, per tier | B |
| S-02 | Cross-tool concordance on 20 fixed seeds: pyhmmer vs MMseqs2 profile vs DIAMOND (and HH-suite for remote cases). Jaccard of representative hit sets at matched bitscore reported; disagreements sampled and adjudicated. From ST0, Domainator (`domainate.py`, `select_by_cds.py`) is run as an independent implementation on the golden genomes: per-CDS domain calls and neighborhood membership are compared with bio-box L1 and `bio-box nbr` at the same models and score cutoffs. | O-SEARCH (differential: pyhmmer, MMseqs2, DIAMOND); S-02d: O-SCAN, O-NBR | R4 from S1; S-02d R3 | R (report); the Domainator comparison is a separate blocking test, S-02d (`data_layer_and_tests.md` §4.10): exact agreement on calls above cutoff + 10 bits at S0 and at every container rebuild |
| S-03 | Z invariance: the same recipe at stage k and k+1 gives identical bitscores and identical E at Z_ref for targets present in both. | O-SEARCH (MR-17, MR-12) | R4, each stage boundary | B |
| S-04 | Controls: recipes without a positive control are refused; the positive must be recovered; the shuffled-seed negative returns no more hits than its expectation at Z_ref. | O-SEARCH (MR-09) | R1 (refusal), R2 (manifest controls), R4 | B |
| S-05 | Provenance: every `search_hit` row resolves to a `search_run`; an injected OOM and an injected timeout are recorded as `failed`, never as zero hits. | O-RUN (F-10) | R2 (fault injection), R4 | B |
| S-06 | Envelope: each question type's cost on the pod is inside its budget in §5; no search process exceeds 300 CPU-min; `/tmp` quota respected. | O-RUN (F-15) | R4 on the search host | B |
| S-07 | Drift alarm: a seed known to corrupt (a promiscuous-domain fusion) triggers the foreign-domain alarm and the revert. | O-ITER (planted traps, `truth_traps`) | R2, R4 | B |
| S-08 | Projection accuracy: rep-projected hits vs exact re-score on a sample of 10^5 members agree on ≥99% of calls above cutoff + 10 bits; all near-cutoff calls use the exact score. | O-PROJECT (F-04, MR-08) | R2 (all members exact), R4 (10^5 sample) | B |
| S-09 | Neighborhood semantics: synthetic contigs test each edge policy, circular wrap, strand and operon-gap predicates, and censoring flags. | O-NBR (MR-03, MR-04, MR-05; F-01, F-09) | R2, R4 | B |
| S-10 | Human dry run: a lab member runs one synteny hypothesis end to end from a written recipe without help. | one workflow manifest, all its ops | R5, per tier (human) | human, once per tier from ST1 |
| S-11 | Temporal leakage: the B1 `dict@2013` build contains no PADLOC-DB candidate systems or other post-2013 models. | O-SCAN, O-SYS | R5 (B1 holdout) | B (at ST3) |
| S-12 | Recipe expressiveness: the three Domainator case studies (CARF uncovered segments, Psu–integrase hotspots, orphan-MTase subtraction) and five synteny hypotheses of different chains (e.g. SYN-143 Q-ARCH, SYN-114 Q-RNA, SYN-162 Q-NBR, SYN-254 Q-MOTIF, SYN-343 Q-SYS→Q-ITER) are written in the recipe language with no escape to ad hoc scripts. Any step that cannot be expressed is a recipe-language gap to close before ST1. | workflow manifests (§4 sibling manifest) | R0 lint at S0; R5 at the replay's tier | B (at ST0, on paper; executed at the replay's tier) |
| S-13 | Schema stability: `--json` outputs of `search`, `nbr` and reports are checked against a versioned schema; a removed or renamed key without a version bump fails. | O-EXPORT | R2, every merge and image build | B (from ST0) |

Rings (TEST_ARCHITECTURE §8): R0 static lint; R1 unit and property; R2 micro-corpus and synthetic worlds; R3 golden set, nightly on the pod; R4 stage gate on stage data; R5 calibration, benchmarks and replays, per tier. The gate column is unchanged; rings say where and how often a test runs before its gate.

**Workflow manifests.** The S-12 recipes (three Domainator case studies, five synteny hypotheses) become the first workflow manifests (TEST_ARCHITECTURE §5, §10 step 4). Every new search workflow enters as a recipe plus its sibling manifest (§4); the harness then generates its schema, MR, control, truth-recovery, budget (S-06) and JSON-schema (S-13) tests with no harness code. A workflow that needs an operation missing from the registry is refused until the operation is registered (TEST_ARCHITECTURE §5.2).

---

## 8. Changes to the plan files

| file | change |
|---|---|
| `00-spec.md` | Invariant 15: searches are recorded recipe runs with provenance, fixed Z and locus context on every hit. Invariant 16: nucleotide search corpora (IGRs, dark ORFs, spacers) are first-class release artifacts. |
| `PLAN.md` | New §4.7 "Ad hoc search" pointing here; base-scan unit per SD-05 (fix PLAN:118); Z per SD-04; `protein_uids`/`slot_uids` widths fixed to 16-byte UIDs in tapes (dense IDs only in L5); `contig_membership` → `contig_member`; `protein_member` replaced by `cluster_member_aln`; PADLOC described as its own engine; invariant count 16 (17 after the test-architecture round adds invariant 17). |
| `STAGES.md` | Insert S2C card and ladder row; add search tier column and per-stage search gate; S−1 adds the complete-genome count, replay supplement fetch, Z constants and image tool list; §8 open items. |
| `DEPLOYMENT-AND-MALLEABILITY.md` | CLI from §3.3; image `bio-box:3.1` adds HH-suite3, BLAST+, DIAMOND, CMfinder, R-scape, PADLOC, CRISPRCasTyper, CRISPRCasFinder; search pack in editions (`lite` ships the S5 pack: NR50, IGR, dark ORF, spacers; `core` ships S6 pack on Hetzner); interactive lane; sidecar coverage definition. |
| `planning/data_layer_and_tests.md` | Tables in §3.1; search pack layout; S-family tests; fix §2A.4 stage labels; Z policy; scan unit. |
| `planning/sources_and_sizing.md` | S2C row; search-pack sizing; NR90/NR50 estimates labeled; S5 residue figure reconciled. |
| `PRIOR-ART.md` | FLSHclust, Domainator 2024 (uncovered-segment clustering is published; the novelty is the combination with RNA features and dependence-aware statistics), Rousset 2022 framework bracketing. Domainator replay stages corrected (case 1 at ST2, case 3 at ST3/ST4, case 2 at ST4) and its role as the candidate workbench (SD-17). |
| Domainator round (2026-10-29) | SD-17..SD-21 and S-12, S-13 added here; `data_layer_and_tests.md` gets the `search_hit` profile coordinates, SSN threshold table, guardrail fields and S-12/S-13; `DEPLOYMENT-AND-MALLEABILITY.md` gets Domainator in the image, the `--json`/`--top_n`/`--format domainator-gb` conventions and `bio-box search ssn`. |
| `planning/preflight_gaps.md` | Add search gaps G-23..G-28. |
| Test-architecture round (2026-10-29) | `planning/tests/TEST_ARCHITECTURE.md` added (operations O-*, failures F-*, MR-01..MR-19, manifests, micro-corpus and synthetic worlds, LY-*, rings R0–R5, `test_run`). Here: §7 gains op and ring columns; S-12 recipes become the first manifests; §4 recipes get a sibling manifest and controls it reuses (synthetic preferred); ST0 tests first on the micro-corpus. `DEPLOYMENT-AND-MALLEABILITY.md` gets the `bio-box-test` image, test hosts per ring and `test_run` storage. |

---

## 9. Open items for the user

1. ~~Approve S2C~~ **Approved 2026-10-29:** ≤20 complete genomes per species. Optional S2C-b (the remaining complete genomes) only if the first-10 questions include strain-level work.
2. ~~Z constants~~ **Decided 2026-10-29 (DEC-04):** `Z_prot = domZ = 1e8`, `Z_nt = 1e5` Mb.
3. **Search host for ST4:** size and timing of a Hetzner box that holds the S6 pack.
4. **Synteny corpus repair:** about 46 `lit/` files hold a different paper than their name (9 in 2022–2026, ~25 in 2017–2021, 12 in 1993–2016), and the synteny audit counts are zero for all 365 rows. Replays use supplements, which are unaffected, but the hypotheses drawn from mislabeled texts need re-extraction before they drive experiments.
5. **Fetch list for replays:** Doron 2018 S1–S3, Gao 2020 S1–S3, Altae-Tran Data S1–S4, Millman Table S1, Brewer 2021 files 1–2, CRISPR-Cas 2025 Zenodo PSSMs, Domainator Zenodo data (10.5281/zenodo.10989173) and analysis scripts (Zenodo 14056380, GitHub `nebiolabs/domainator_examples`) (`03` §9, `02` §9, `05` §f).
6. **Test-architecture items** (generator: ZOMBI2 plus custom code, decided; R2 CI host: GitHub Actions, decided; mutation-score target; flexibility-check owner): see `planning/tests/TEST_ARCHITECTURE.md` §11.
