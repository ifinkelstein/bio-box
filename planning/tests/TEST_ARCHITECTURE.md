# Test architecture: how the suite is organized as the database layout starts

Status: **planning only** (2026-10-29). Nothing here is implemented. This file organizes the tests already defined in `planning/data_layer_and_tests.md` §4 (families C, L, R, Q, A, D, G, K, S and the SLOs) and `PLAN.md` §5 (B1–B6), and adds what they lack. Test IDs in those files keep their meaning; this file adds an index over them.

---

## 0. How this brief was rewritten

**1. What a top professional would do that was not asked.**
- Organize tests around **operations**, not layers or stages. Every planned workflow is a composition of a small set of operations (slice, translate, call, scan, cluster, project, tape, pair, search, iterate, neighborhood, architecture, statistics, network, export, catalog, accrete). If each operation has a contract, an oracle and metamorphic relations, a new workflow built from tested operations inherits coverage. Only a new operation needs new test code.
- Rank failures by **silent wrong answers**, not crashes. A crash is found the first day; an off-by-one strand, a dense ID that leaks across releases, or a null that is too easy produces confident, wrong biology for months.
- Build **synthetic worlds with planted truth** (gene families, operons, fusions, RNA motifs, spacers, duplicates, contamination, fragmentation), including **null worlds** with phylogeny and HGT but no coupling, so the statistics can be calibrated, not just exercised.
- Use **metamorphic relations** where no ground truth exists: transformations of the input with a known effect on the output (reverse-complement, rotate, cut, duplicate, rename, add unrelated genomes, shard differently).
- **Test the tests.** Mutation testing shows the suite catches injected bugs; a lab member adding a workflow without editing the harness shows the suite is flexible.
- Write the harness **before** the layout. Layout choices still open (Parquet chunks vs bgzip, row-group sizes, sort keys) are decided by tests that count row groups read, not by argument.
- Store **test results as data** (a `test_run` table), so a stage gate is a query and a regression is a diff.

**2. Rawest material.**
- The 625 recipe rows in `planning/search/recipes_*.tsv` and the 365 classified hypotheses (`synteny_hypothesis_query_types.tsv`): the operations real workflows use, with their parameters.
- Raw tool output on the golden genomes: pyrodigal GFF, `--domtblout` from hmmsearch, Infernal `--tblout`, MinCED and CRISPRCasFinder output, MMseqs2 cluster TSVs. These are the outputs whose parsing breaks, not summaries of them.
- Published hit lists in supplements (`planning/search/01` §5, `02`–`05` replay sections): external truth.
- Format and tool specifications: the Parquet format spec (row groups, statistics, FIXED_LEN_BYTE_ARRAY), DuckDB docs on pruning and `EXPLAIN ANALYZE`, HMMER and Infernal user guides on `-Z`/`--domZ`, the genetic-code tables.
- Issue trackers of the tools in the image (MMseqs2 `clusterupdate` issue #961, DuckDB Parquet and concurrency issues, pyhmmer, Infernal): real failure modes found by others.

**3. Methods and tools.**
- `pytest` with markers per ring and operation; `Hypothesis` for property-based tests and a `RuleBasedStateMachine` for the accretion ledger (append, recluster, rescore, compact, remap).
- A schema registry as code (pyarrow schemas plus semantic constraints) checked with Pandera and plain DuckDB SQL assertions (dbt-duckdb is not used; `data_layer_and_tests.md` §4).
- Brute-force oracles in pure Python/SQL on the micro-corpus (already the Q family's approach).
- Differential testing against independent implementations: pyhmmer vs the `hmmsearch` binary, DuckDB vs pyarrow vs polars readers, `bio-box nbr` vs Domainator (S-02d), MMseqs2 vs DIAMOND vs pyhmmer (S-02).
- Simulation: Zombi (species trees, genome evolution with duplication, transfer, loss, inversion and transposition; sequence mode via Pyvolve) for gene order and families; `hmmemit` from real Pfam HMMs for family members at controlled divergence; `cmemit` from Rfam CMs for planted RNA motifs; planted CRISPR arrays with spacers matching planted phage contigs; a fragmenter for assembly quality; a mutator for gene-call errors (frameshifts, wrong starts). Fit of Zombi to be checked at S−1 (ALF is the fallback).
- Statistical calibration: simulation-based checks that null p-values are uniform (KS test, QQ against permutation), FDR at nominal level, power curves under planted signal; tolerance bands from binomial CIs with fixed seeds.
- DuckDB `EXPLAIN ANALYZE` and profiling JSON to count files and row groups read; `parquet-tools`/pyarrow metadata to check sort order, statistics and sizes.
- Fault injection: low `MemoryMax` scopes, `timeout`, `kill -9` during write and during catalog rename, `/tmp` filled to quota, a read-only NFS mount, the STOP file, a simulated watchdog renice.
- Mutation testing with `mutmut` on the operation library; coverage with `coverage.py`; performance history with `asv` (airspeed velocity) or `pytest-benchmark`.

**4. What, other than us, shows the result is right.**
- Independent implementations agree (brute-force oracles, pyarrow vs DuckDB, the `hmmsearch` binary vs pyhmmer, Domainator vs `bio-box nbr`).
- Planted truth in synthetic worlds is recovered at measured rates, and null worlds give calibrated p-values.
- Published hit lists are recovered (S-01).
- Injected mutants are killed (mutation score target below).
- A lab member adds a new workflow by writing only a manifest and a recipe, and it gets full coverage without changes to the harness.

**Rewritten prompt.**
> Design the test suite for the bio-box data layout and planned workflows, before any layout is built. Derive an operation registry from the recipe tables and hypotheses; give every operation a contract, an oracle and metamorphic relations, so workflows composed from tested operations inherit coverage and a new workflow needs only a manifest. Rank failure modes by silent wrong answers and map each to tests. Specify a seeded synthetic-world generator with planted truth and null worlds, a micro-corpus that runs end to end in CI, the test rings from unit to stage gate, layout tests that decide the open layout questions by counting row groups read, fault injection, statistical calibration, mutation testing and a `test_run` results table. Map every existing test ID into this structure and name the gaps. Success is shown by independent implementations agreeing, planted truth recovered, published hit lists recovered, mutants killed, and a lab member adding a workflow without touching the harness.

---

## 1. Principles

1. **Operations are the unit of coverage.** A workflow is a DAG of registered operations (§2). Tests attach to operations; workflows add controls and truth (§5).
2. **Every operation has three kinds of evidence:** a contract (types, keys, invariants), an oracle (brute force or independent implementation), and at least one metamorphic relation (§4).
3. **Silent wrong answers first** (§3). Crash tests are cheap and come along; they are not where effort goes.
4. **Truth comes from four places, in this order of trust:** constructed truth (synthetic worlds), independent implementations, published results, golden snapshots. Golden snapshots catch change, not error.
5. **Fixed seeds, no retries.** Blocking tests are deterministic. Statistical tests use fixed seeds and tolerance bands from binomial CIs. A flaky test is quarantined with an expiry date, never retried into passing.
6. **Results are data.** Every test run writes `test_run` rows (§8); gates are SQL over it.
7. **Planning-only constraint respected.** This document specifies; nothing is implemented until S−1 signs off.

---

## 2. Operation registry

Each operation has an ID, inputs and outputs (registered schemas), its oracle, its metamorphic relations (§4) and the existing tests that already cover it.

| op | what it does | oracle | MRs | existing tests |
|---|---|---|---|---|
| O-INGEST | source object → `l0`, `contig_occurrence`, `seq_chunks` | re-parse the source independently; sha256 of every record | MR-01, MR-02, MR-14 | R-01, R-02, C-01, L-01 |
| O-IDENT | identity resolution and dedup (L0.5, `contig_member`) | brute-force pairwise containment on micro-corpus | MR-06, MR-07, MR-01 | C-*, A-* (dedup), D-* |
| O-SLICE | `seq get` contig, region, strand; `seq_slice` macro | Python slicing of the source FASTA | MR-03, MR-04, MR-14 | R-03, R-05, R-08 |
| O-TRANSLATE | region and six-frame translation, genetic codes | Biopython translation tables | MR-03, MR-13 | R-04 |
| O-EXTRACT | bulk FASTA export, `--members`, `--segment` | re-ingest and compare digests | MR-01, MR-16 | R-07 |
| O-CALL | gene calling, smORFs, split-gene arbitration | planted genes in synthetic worlds; golden GFF | MR-03, MR-04, MR-13 | G-*, R-04 |
| O-SCAN | library HMM/CM scans, Z policy | `hmmsearch`/`cmsearch` binaries on the same inputs | MR-12, MR-15, MR-16 | G-*, A-02, S-03 |
| O-ARRAY | CRISPR arrays and spacers | planted arrays; MinCED vs CRISPRCasFinder union | MR-03, MR-04 | G-* |
| O-SYS | DefenseFinder, PADLOC, CRISPRCasTyper sidecars | tool's own test data; planted systems | MR-01 | G-*, S-11 |
| O-CLUSTER | NR90/NR50/subfamily, `cluster_member_aln`, `clusterupdate` | full recluster on micro-corpus; pairwise identity check of members | MR-01, MR-08, MR-16 | A-*, C-11, D-* |
| O-PROJECT | rep → member projection and near-cutoff re-score | exact re-score on all members (micro) or a sample (stage) | MR-08, MR-12 | S-08, A-03 |
| O-TAPE | contig tapes, segments, junction topology, directons | brute-force rebuild from L1 in Python | MR-03, MR-04, MR-05, MR-11 | Q2 oracle, G-* |
| O-PAIRS | `slot_pairs`, window enumeration | brute-force window enumeration | MR-03, MR-05, MR-11 | Q1, Q4 oracles |
| O-SEARCH | recipe search against a pack target | `hmmsearch`/`mmseqs`/`diamond` binaries; exhaustive search on `unique` | MR-09, MR-10, MR-12, MR-16 | S-01, S-02, S-03, S-04 |
| O-ITER | iteration controller | hand-built rounds on micro-corpus; planted divergent members | MR-09, MR-10 | S-07 |
| O-NBR | neighborhood fetch with window, strand, operon, edge, circular rules | brute-force Python over coordinates | MR-03, MR-04, MR-05 | S-09, S-02d |
| O-NBRSCORE | positive-set enrichment | brute-force counting on micro-corpus | MR-06, MR-18 | (new) |
| O-ARCH | fusion and architecture classification | planted fusions and splits | MR-11, MR-03 | (new) |
| O-STATS | pair statistics, nulls, FDR, modules | brute-force score test; simulation calibration | MR-06, MR-08, MR-18 | Q3 oracle, B1–B6 |
| O-ANCHORSTATS | on-demand statistics for one anchor or user token | equals O-STATS when the anchor is a library token | MR-06, MR-18 | (new) |
| O-PROMOTE | user token overlay | overlay rows equal the run's final hit slots | MR-01, MR-08 | C-03 |
| O-PROF | presence profiles and pangenome planes | brute-force presence matrix | MR-06, MR-07 | Q5/Q6 oracles |
| O-SSN | hit-set network, MSF ∪ kNN, threshold table | full graph and union-find on micro-corpus; Domainator `build_ssn` | MR-01, MR-16 | (new) |
| O-EXPORT | `--json`, `--format domainator-gb`, TSV | parse back and compare; Domainator reads it | MR-01 | S-13 |
| O-CATALOG | manifest → catalog, atomic swap, read-only open | list files independently; compare views | MR-01 | C-*, K-05, K-06, K-08, K-09 |
| O-SERVE | query classes Q1–Q8 | brute-force oracle per class | per class | Q1–Q8, SLOs |
| O-ACCRETE | append stage, recluster, rescore, compact, remap | from-scratch rebuild | MR-08, MR-01 | A-01..A-07, L-* |
| O-RUN | runner: provenance, guardrails, failure states | fault injection | none | S-05, S-06 |

A new operation enters the registry only with all four columns filled. The registry lives as data (`tests/registry/operations.yaml`) and the coverage matrix (§5.3) is computed from it.

---

## 3. Failure modes ranked by harm

Severity: **S1** silent wrong biology that ships; **S2** wrong numbers that are flagged or visible; **S3** crash or refusal. Effort follows severity.

| # | failure | severity | where it hides | tests that catch it |
|---|---|---|---|---|
| F-01 | Coordinate or strand error (off by one, half-open mixups, minus-strand slices, origin-spanning genes) | S1 | O-SLICE, O-TAPE, O-NBR | MR-03, MR-04, R-03, R-05, Hypothesis properties on `seq_slice`/`revcomp` |
| F-02 | ID mismatch across releases (dense IDs used as join keys, UID width truncation, collision) | S1 | O-CATALOG, sidecars, overlays | C-02 (uniqueness and collision), C-03, schema registry type checks (16-byte UIDs, LY-01), MR-08 |
| F-03 | Dedup over- or under-collapse | S1 | O-IDENT | MR-06, MR-07, planted duplicates, brute-force containment |
| F-04 | Projection error from representatives to members | S1 | O-PROJECT | S-08, exact re-score on micro-corpus |
| F-05 | E-values drift with database size; thresholds not comparable | S1 | O-SCAN, O-SEARCH | S-03, MR-12 |
| F-06 | Stale statistics served as current | S1 | O-ACCRETE, O-SERVE | A-*, `stale` flag tests, MR-08 |
| F-07 | Null too easy (phylogeny or clonality ignored), so FDR is not controlled | S1 | O-STATS | null worlds (§6.2), MR-18, B5 |
| F-08 | Gene-call errors create false fusions or splits, or false absences | S1 | O-CALL, O-ARCH | planted frameshifts and wrong starts, dark-ORF checks, MR-11 |
| F-09 | Contig-edge censoring wrong (absence claimed where the window was cut) | S1 | O-NBR, O-STATS | MR-05, S-09, B3 |
| F-10 | Tool failure read as zero hits | S1 | O-RUN | S-05, fault injection |
| F-11 | Partial write visible to readers (crash mid-write or mid-swap) | S2 | O-CATALOG | `kill -9` during write and rename; K-09 (writer exclusivity); LY-10 |
| F-12 | Non-redistributable data in an edition | S2 (legal) | editions | license tests on every pack and edition manifest |
| F-13 | Nondeterminism (thread count, DuckDB version, shard count) | S2 | all writers | D-01..D-05, MR-01, MR-16 |
| F-14 | Layout regression (sort lost, statistics missing, small files) | S3 | all writers | layout tests (§7) |
| F-15 | Envelope breach (memory, CPU-min, `/tmp`) | S3 | O-RUN | S-06, SLO tests |

---

## 4. Metamorphic relations

A metamorphic relation (MR) states how the output must change when the input is transformed. MRs apply wherever an operation declares them (§2) and are run automatically on the micro-corpus and on the golden set.

| MR | transformation | required effect |
|---|---|---|
| MR-01 | permute input order (records, files, genomes) | byte-identical outputs at fixed DuckDB version and threads |
| MR-02 | rename accessions or headers | identical UIDs and annotations; only occurrence rows change |
| MR-03 | reverse-complement a contig | same hits and neighborhoods; coordinates mirrored; strands flipped; previous/next junctions swapped |
| MR-04 | rotate a circular contig | same genes and neighborhoods; no edge flags; the origin-spanning gene handled |
| MR-05 | cut a contig at position p | only windows crossing p change; those gain edge flags; under `edge: exclude` the count drops by exactly those loci |
| MR-06 | add an exact duplicate genome | raw counts change; vehicle- and BioProject-deduplicated counts, Z and q-values do not |
| MR-07 | add a 99.9% ANI near-duplicate | collapsed for counting, kept in storage; R-01 round trip holds for both |
| MR-08 | add unrelated genomes (accretion) | exact hits on existing targets unchanged (bitscore and E at Z_ref); statistics change only where flagged stale |
| MR-09 | shuffle the seed (composition preserved) | hits no more than expected false positives at Z_ref |
| MR-10 | plant one more family member in a new genome | recall includes it; counts rise by exactly one |
| MR-11 | split a fused ORF into two adjacent ORFs | architecture flips fused → split-adjacent; topology counts move by exactly one |
| MR-12 | search a subset of the target | a subset of hits with identical bitscores |
| MR-13 | recode a genome to genetic code 4 or 25 | translation follows the code; code detection picks it |
| MR-14 | lowercase soft-mask or IUPAC variants | same `contig_uid`, different `raw_digest`; exact round trip |
| MR-15 | search the reverse complement with a nucleotide model | same hits on the other strand |
| MR-16 | change shard count or thread count | identical results |
| MR-17 | rerun a recipe on a later release | same results on unchanged targets (generalizes S-03) |
| MR-18 | null world: relabel genes within phylogenetic and HGT constraints | p-values uniform; FDR at nominal level |
| MR-19 | add a contamination mask (an FCS-GX or host-screen sidecar flags contigs) | sequences, tapes and UIDs unchanged; only statistics involving masked contigs change, and they change exactly as if those contigs had been excluded at ingest; `fcs_gx_status` updated |

---

## 5. Workflows: how new ones get coverage

### 5.1 Workflow manifest

Every workflow (base pipeline step, query class, or search recipe pattern) is registered with a manifest next to its recipe:

```yaml
workflow: rosetta_fusion_census
ops: [O-SEARCH, O-ITER, O-PROJECT, O-NBR, O-ARCH, O-ANCHORSTATS, O-EXPORT]
recipe: recipes/rosetta_uri_dnaq_v1.yaml
outputs: [search_hit, neighborhoods, arch_calls]     # registered schemas
controls:
  positive: {synthetic: planted_fusion_set_A, real: MJ0613}
  negative: {synthetic: null_world_1, transform: MR-09}
truth: {synthetic: planted_fusions, replay: null}    # replay optional; justify if none
budget_class: Q-ARCH                                 # SEARCH_PLAN §5
```

### 5.2 What the harness generates from a manifest

- Schema contract tests for each declared output.
- The MRs of every listed operation, run on the micro-corpus with this workflow's recipe.
- Positive and negative control tests.
- A truth-recovery test (recall and precision against planted truth, with bands).
- A budget test (S-06) at the declared class.
- A JSON schema test (S-13) for its exports.

A new workflow made only of registered operations therefore needs a manifest and a recipe, and no harness code. A workflow that needs a new operation needs that operation's registry entry (§2) first; CI refuses the workflow until it exists.

### 5.3 Coverage matrix

Computed from the registry and manifests: rows are workflows, columns are operations, cells say which tests ran and passed. Gaps are reported per stage. The initial workflow set is the 13 question types (Q-FAM … Q-SYS), the base pipeline steps (ingest through L4), the eight query classes, and the S-12 recipes (three Domainator cases, five synteny hypotheses).

---

## 6. Test data

### 6.1 Micro-corpus (CI, minutes)

- About 50 synthetic genomes from one Zombi species tree, plus three small real genomes (λ, T4, *M. genitalium*), built end to end in under 15 minutes on a laptop with a reduced model set (about 200 Pfam HMMs, 20 Rfam CMs).
- Every layer, the search pack, and every query class are built; brute-force oracles run on all of it.
- Seeded and deterministic; outputs hashed and recorded as the micro-corpus manifest.

### 6.2 Synthetic worlds (planted truth)

The generator (`synthworld`, versioned, seeded) writes genomes plus truth tables:

| planted feature | how | truth table |
|---|---|---|
| gene families | Zombi gene trees; sequences by `hmmemit` from real Pfam HMMs at set divergence | `truth_family` |
| modules (true coupling) | co-transferred gene blocks with set strength | `truth_module` |
| fusions and splits | ORF merge and split events on chosen lineages | `truth_fusion` |
| RNA motifs in IGRs | `cmemit` from Rfam CMs upstream of chosen families | `truth_rna` |
| CRISPR arrays and protospacers | arrays with spacers drawn from planted phage contigs | `truth_spacer_target` |
| duplicates and near-duplicates | exact copies and 99.9% ANI variants across "sources" | `truth_dup` |
| contamination | foreign contigs in some genomes | `truth_contam` |
| fragmentation | cuts to a target N50 | `truth_edges` |
| gene-call errors | frameshifts, wrong starts, missed small ORFs | `truth_callerr` |
| hard cases | promiscuous domains (HTH, P-loop) fused widely; MGE cargo hitchhiking | `truth_traps` |

**Null worlds** use the same tree, HGT and gene-order processes with no modules planted. They calibrate O-STATS (MR-18): p-value uniformity, FDR control, and the minimum detectable effect. A null world without phylogeny and HGT would be too easy and is not used.

### 6.3 Golden set and replays

- The S0 golden genomes, CAMI II subset and spike-ins (`STAGES.md` S0), with snapshot hashes for G-* tests.
- Published replays by tier (SEARCH_PLAN §6, S-01).

---

## 7. Layout tests (decide the open layout questions)

These run on the micro-corpus at S−1/S0, and on each stage's real data at its gate.

| ID | test | decides or protects |
|---|---|---|
| LY-01 | Every table matches its registered schema: types (16-byte UIDs as FIXED_LEN_BYTE_ARRAY), nullability, enum closure | F-02 |
| LY-02 | Declared sort order holds within each file (checked, not assumed) | pruning, F-14 |
| LY-03 | Min/max statistics present for every sort-key column in every row group | pruning |
| LY-04 | File size and row-group size within the bands in `data_layer_and_tests.md` §3.2; no file under 64 MB on NFS outside the waivers. Waivers: search runs and user overlays (`data_layer_and_tests.md` §2B.3), `test_run` session files (§8), and metadata files declared `small_ok: true` in the schema registry (`l0_manifest.parquet`, `seq_locator.parquet`, `id_history.parquet`, ledger files, release manifests, run and test indexes, small L4 tables), each bounded to a declared count per release | NFS rule |
| LY-05 | Partition keys in paths match column values | correctness of pruning |
| LY-06 | (generalizes K-05/K-06 to every canonical query) Canonical point queries (one contig slice, one neighborhood, one anchor) read at most N units, counted from `EXPLAIN ANALYZE` / profiling output, not wall time. The unit is the Parquet row group; for a bgzip candidate it is BGZF blocks decompressed plus index seeks, so the two layouts are compared on bytes decoded per query | row-group sizing, sort keys, K-13 (Parquet chunks vs bgzip) |
| LY-07 | DuckDB, pyarrow and polars read identical rows and checksums from every table | reader interoperability |
| LY-08 | The same inputs written twice give identical bytes (fixed DuckDB version and threads) | F-13 |
| LY-09 | Schema evolution: adding a nullable column leaves older catalogs readable; removing or retyping one fails CI | stability |
| LY-10 | The catalog lists files explicitly (no globs) and every manifest file exists with its sha256 | F-11 |

LY-06 replaces latency arguments in the open layout decisions (K-13 sequence-store format, row-group sizes, sort keys, and K-07 `LIST` tapes vs exploded slot tables): the layout that decodes fewer units for the canonical queries wins, and latency SLOs are checked afterwards. Whether to build an L5 index (K-12) stays SLO-driven, because L5 exists only when an SLO fails.

---

## 8. Rings, cadence and the results table

| ring | content | data | where | cadence | budget |
|---|---|---|---|---|---|
| R0 static | schema registry lint, recipe and manifest lint, registry completeness, image manifest | none | laptop, CI | every commit | < 1 min |
| R1 unit and property | operation functions; Hypothesis properties; stateful accretion machine on toy data | generated | laptop, CI | every commit | < 5 min |
| R2 micro-corpus | full build; oracles; all MRs; workflow manifests; LY tests; fault injection; mutation testing (weekly) | micro-corpus | CI runner or laptop | every merge | < 30 min |
| R3 golden | S0 golden set with real tools and full model sets; G snapshots; S-02d | golden set | pod | nightly and every image build | < 6 h at nice 19 |
| R4 stage gate | sampled oracles, accretion, determinism, K tests, S tests, SLOs; one scaling-fit point per stage (S1, S2, S3 included) | stage data | pod or Hetzner | per stage | as in STAGES |
| R5 calibration and benchmark | null-world calibration at stage scale where the statistics change scale (S2C, S5, S6; other tiers rerun it at micro-corpus scale), B1–B6, replays, the scaling fit over the R4 points | synthetic worlds at stage scale; real | pod or Hetzner | per tier | as in SEARCH_PLAN |

**Where the tooling lives.** Test tooling (pytest, Hypothesis, pandera, mutmut, asv, Zombi with Pyvolve, ALF as fallback) ships in a separate `bio-box-test` image layered on `bio-box:3.1`, so the runtime image stays unchanged; `hmmemit` and `cmemit` come from the HMMER and Infernal already in the runtime image (`DEPLOYMENT-AND-MALLEABILITY.md` §2.1). R0–R2 run on the CI host (open item §11.2); R3 runs nightly on the pod at nice 19 inside the run-step envelope; R4–R5 run on the stage host; multi-client tests run at night (`planning/preflight_gaps.md` #20).

**Results table.** Every test run writes `test_run` rows: `test_id`, ring, operation, workflow, params hash, outcome, metrics (JSON), duration, code commit, `runtime_image_digest` (the `bio-box:3.1` binaries under test), `test_image_digest` (the `bio-box-test` layer), `release_id`, host, seed. Stored as packed Parquet under `$ATD_ROOT/tests/`. Rows from an off-pod CI host are written as Parquet artifacts and pulled by the pod-side single writer (rsync over SSH from the pod), then registered; the transport is fixed with the CI-host choice (§11 item 2). A stage gate is a SQL query over `test_run` (all blocking tests for the stage pass on the stage's release and image). Regressions are diffs between runs.

Further columns and rules (first written in `data_layer_and_tests.md` §4.12, adopted here):
- `started_at` (timestamp) orders runs for regression diffs; `outcome` is one of `pass`, `fail`, `error`, `skip`.
- The gate specification is data: `tests/registry/gates.yaml`, beside `operations.yaml`, generated from the stage tables in `data_layer_and_tests.md` §4.9–§4.10, with ring requirements expanded to every registered test of those rings. It loads as `gate_spec(stage, test_id, gate)`; a stage passes when no blocking `gate_spec` row lacks a `pass` row for the stage's release and both image digests.
- **Quarantine.** A quarantined test is listed in `tests/registry/quarantine.yaml` with an owner, a reason and an expiry date. Before expiry the gate query treats it as recorded (W), not blocking; at expiry it is blocking again and fails until fixed. S1-severity tests (§3) cannot be quarantined.

**Tests of the tests.**
- Mutation score on the operation library: ≥ 80% of mutants killed by R1–R2 (*target*; measured at S0).
- Every S1-severity failure (§3) has at least one test that fails when the failure is injected deliberately (seeded fault injection list maintained with the registry).
- Flexibility check: a lab member adds one new workflow from a manifest and recipe only; the coverage matrix fills without harness edits.

---

## 9. Mapping of existing tests and gaps

| existing family | where it fits | gaps it leaves |
|---|---|---|
| C-01..C-12 contracts | R0/R2 schema registry and O-* contracts | schema registry not yet code; LY-01 added |
| L-01..L-04 ledger | O-ACCRETE, O-RUN | stateful testing of operation sequences (Hypothesis machine) |
| R-01..R-09 sequence round trip | O-INGEST, O-SLICE, O-TRANSLATE, O-EXTRACT oracles | property tests on coordinates (F-01) |
| Q1–Q8 oracles | O-SERVE oracles | run on synthetic worlds as well as golden data |
| A-01..A-07 accretion | O-ACCRETE, MR-08 | random operation orders (stateful) |
| D-01..D-05 determinism | MR-01, MR-16, LY-08 | shard-count invariance for search |
| G-01..G-05 golden | R3 snapshots | snapshots catch change, not error; paired with planted truth |
| K-01..K-13 DuckDB | LY tests, O-CATALOG | row-group counting (LY-06) as the decision rule |
| S-01..S-13, S-02d search | O-SEARCH, O-ITER, O-NBR, O-PROJECT, O-RUN, workflow manifests | none new; S-12 becomes the first manifest set |
| SLOs | R4 | performance history across commits (`asv`) |
| B1–B6 benchmarks | R5 | null-world calibration complements B5 |

**New in this file:** the operation registry; the failure ranking; MR-01..MR-19; workflow manifests and the coverage matrix; the micro-corpus and the synthetic-world generator with null worlds; LY-01..LY-10; rings; the `test_run` table; mutation testing; deliberate fault injection for every S1 failure; the flexibility check.

---

## 10. Build order (when implementation is approved)

1. Schema registry (R0) and the operation registry with empty cells.
2. Synthetic-world generator and the micro-corpus, with truth tables.
3. Brute-force oracles for O-SLICE, O-TAPE, O-PAIRS, O-NBR, O-STATS.
4. MR harness and the first manifests (S-12 recipes).
5. Layout tests LY-01..LY-10, run against candidate layouts to settle K-13, row-group sizes and sort keys.
6. Fault injection and the `test_run` table.
7. R3 on the pod with the golden set; then S0 proceeds.

---

## 11. Open items

1. Zombi fit for gene-order and sequence simulation at micro-corpus size (fallback ALF); owner S−1.
2. CI host for R2: laptop, GitHub Actions or a pod-side runner (licenses and image size decide); owner S−1.
3. Mutation-score target confirmed after the first measurement at S0.
4. Which lab member runs the flexibility check, and on which workflow.
