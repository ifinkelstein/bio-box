# `all-the-data`: Packaging (`bio-box`), Data-Pack Malleability, and Deployment (v3.0)

> **Status:** v3.0, revised after the round-2 audits ([`audits/SYNTHESIS_2.md`](audits/SYNTHESIS_2.md), decisions D-02, D-03, D-17, D-19..D-26). Supersedes v2.0 of this file.
> **Related:** [`STAGES.md`](STAGES.md) (staged accretion plan), [`planning/data_layer_and_tests.md`](planning/data_layer_and_tests.md) (layers, layout, test suite), [`planning/sources_and_sizing.md`](planning/sources_and_sizing.md) (sources and sizing per stage), [`planning/search/SEARCH_PLAN.md`](planning/search/SEARCH_PLAN.md) (ad hoc search; authoritative for SD-01..SD-21, planned), [`planning/tests/TEST_ARCHITECTURE.md`](planning/tests/TEST_ARCHITECTURE.md) (test organization, rings, `test_run`; planned).
> **Naming:** the product, container and CLI are `bio-box` (formerly `koonin-box`).
> **Numbers:** every size and time in this file is either *measured (stage, date)*, *derived (formula shown)*, or *estimate, pilot-pending (±40%)*. None of the v2.0 figures have been measured yet.

## Design goals

1. One pinned, versioned software environment that runs the same way on the current lab host (wilkcomp02, CPU only) and on rented bare metal (Hetzner).
2. Data packs that are built in stages, verified by manifest, and attachable without a rebuild.
3. A precise contract for what changes when data is added, and what becomes stale.

The v2.0 targets of "< 4 ms queries", "< 3 s cold start" and "52-minute planetary bootstrap" are withdrawn. Latency targets are now per query class, measured on the deployment host (see [`planning/data_layer_and_tests.md`](planning/data_layer_and_tests.md)).

---

## 0. Deployment target 1 (current): wilkcomp02, CPU only

Source of facts: `WILKECOMP.org` (probed 2026-09-30) and the watchdog notes of 2026-10-15. The host is a shared UT-Austin BRCF/RCTF pod machine. We are guests on it, and the house rules are hard constraints.

### 0.1 Envelope

| Item | Value | Consequence for `bio-box` |
|---|---|---|
| CPU | 4× Xeon E7-4850 v3 (Haswell, 2.2 GHz), 56 physical / 112 logical, AVX2, no AVX-512 | Build every binary with AVX2 dispatch. Per-thread speed is about one third of a modern core, so every rate must be measured here, not extrapolated. |
| RAM | 1,007 GB physical, shared | Stay under 20% by day and under 65% only overnight for a few hours. Default cap: `systemd-run --user --scope -p MemoryMax=200G`. |
| GPU | None | No ProstT5, no gLM training, no predicted 3Di. GPU stages run elsewhere (wilkcomp03 by arrangement, or a GEX130). |
| Local disk | `/tmp` ~940 GB free, local ext4, volatile | All tool temp directories, MMseqs2 indices and DuckDB `temp_directory` go here. Never hold the only copy of anything. |
| Shared storage | `/stor/work/Wilke` (75 TB, backed up) and `/stor/scratch/Wilke` (30 TB, not backed up), one NFS server over 10 GbE shared by three hosts | Durable small outputs go to `/stor/work/Wilke/if2453/`: manifests, ledgers, release tables, `sources.yaml`, and the sequence-store manifest (§3.10). Regenerable bulk data goes to `/stor/scratch/Wilke/if2453/`: source snapshots, Parquet shards, MMseqs2 DBs, and the sequence store itself, which can be rebuilt from L0 snapshots and the manifest. |
| Scheduler | None. No SLURM, no modules, no cgroup quotas | Orchestration supplies its own resource caps (§0.3). |
| Limits | `ulimit -v` 500 GB soft; open files 4,096 hard (8,192 allowed) | Chunk MMseqs2 targets; raise `ulimit -n` in the launcher. |
| Containers | Apptainer/Singularity present; Docker group use not confirmed | Run `bio-box` as an Apptainer `.sif` only. |
| Network | Outbound internet works; storage server reachable from campus/VPN only | Fetch sources directly on the host; large transfers go to the storage server from campus. |

### 0.2 House rules, as the launcher enforces them

- **Threads:** physical cores only. Cap at 40 threads by day; 56 at night and on weekends, and only when `ps`/`top` shows no other compute. Every tool receives an explicit thread flag; MMseqs2 defaults to 112 otherwise. Numerical Python workers also get `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS` and `MKL_NUM_THREADS` set, because a tool's own worker flag is not a thread cap. DefenseFinder `-w 1` drove load to 67 on 2026-10-15.
- **Priority:** monolithic multi-threaded tools (MMseqs2, Foldseek, HMMER/pyhmmer bulk scans, Infernal bulk scans) start at `nice 19`. The hourly watchdog renices any single process that passes about 369 CPU-minutes and mails rctf-support. Starting at 19 leaves it nothing to change. Per-genome fan-out workers run at `nice 15`. All jobs run under `ionice -c 2 -n 7`.
- **Memory:** every job runs inside a `systemd-run --user --scope -p MemoryMax=<cap>` scope. DuckDB `memory_limit` is set below the scope cap.
- **I/O:** stage inputs to local `/tmp`, compute there, and write results back to NFS as one archive or one Parquet shard per chunk with one `rsync`. No per-genome or per-query files on NFS. Watch `nfsstat -c` retransmits.
- **Pausability:** every orchestrated stage checks a `STOP` file between chunks.
- **Hygiene:** every directory over 300 GB carries a README stating what it is, who owns it, until when, and what is safe to delete.

### 0.3 How `bio-box` runs here

```bash
# Build once (on a workstation or in CI), copy the image to scratch
apptainer build bio-box_cpu_v3.0.sif docker://ghcr.io/<org>/bio-box:3.0-cpu

# Every stage is launched through the workflow profile for this host,
# which wraps each rule in systemd-run, nice/ionice and the thread caps above
bio-box run --profile wilkcomp02 --stage S1 --release rS1.1
```

- The workflow profile (`profiles/wilkcomp02/`) encodes the caps: `max_threads_day=40`, `max_threads_night=56`, `mem_scope=200G`, `tmp=/tmp/if2453/bio-box`, `durable=/stor/work/Wilke/if2453/all-the-data`, `bulk=/stor/scratch/Wilke/if2453/all-the-data`, `nice_monolithic=19`, `nice_fanout=15`.
- **Project checkout on the pod: `~/bio-box`** (decided 2026-10-29). It holds the git repository (plan docs, code, registries, recipes, `profiles/`) and the `.sif` images. It holds no data: bulk data stays under `bulk=` and durable small outputs under `durable=`, so home-directory quota and the NFS small-file rule are not at risk.
- **No public network service.** `bio-box serve` binds to `127.0.0.1` only. Users reach it through an SSH tunnel (`ssh -L 8420:127.0.0.1:8420 if2453@wilkcomp02.ccbb.utexas.edu`). There is no Docker daemon, no exposed port, and no MCP server reachable from outside the host. The Domainator GUI (`domainator_server`, image 3.1) follows the same rule: started per user on a `127.0.0.1` port with `--data-dir` pointing at an exported subset in the user's project directory, reached through the same SSH tunnel, never at the store itself.
- Interactive queries use DuckDB directly against the read-only catalog (§3.8), with `threads` and `memory_limit` set to the daytime caps.

### 0.4 What wilkcomp02 is for, and where it stops

wilkcomp02 builds and serves the early stages of [`STAGES.md`](STAGES.md): golden genomes, phages, plasmids, complete prokaryotic genomes (S2C), small genomes and GTDB species representatives. It also measures the per-unit rates that feed `planning/sizing.yaml`. The stage at which its envelope stops being sufficient (local `/tmp`, 30 TB scratch, the 200 GB scope, the watchdog, shared NFS) is identified in [`planning/sources_and_sizing.md`](planning/sources_and_sizing.md). Stages past that point move to rented hardware (§4).

---

## 1. The rule that does not change: compute image and data pack are separate

Reference libraries, Parquet tables and indexes are never baked into image layers. `bio-box` ships as two kinds of artifact.

1. **Compute images** (software only; §2).
2. **Data packs**: content-addressed, manifest-verified directories of Parquet files, reference model libraries (pressed HMMs, calibrated CMs, Foldseek DBs) and model weights. Mounted at `/data`.

Model weights (ProstT5, any future gLM) belong to the data pack, versioned with the release that used them, and are not part of any image.

---

## 2. Compute images (D-26)

### 2.1 Images, one lockfile

| Image | Contents | Used on |
|---|---|---|
| `ghcr.io/<org>/bio-box:3.0-cpu` (+ `bio-box_cpu_v3.0.sif`) | pyrodigal / pyrodigal-gv, Codetta, pyhmmer / HMMER 3.4, Infernal 1.1.5, MMseqs2 (≥ 18, static AVX2 build), Foldseek (CPU), geNomad, CheckV, IntegronFinder 2, MacSyFinder (DefenseFinder, CONJScan, TXSScan models in the data pack), ISEScan, CheckM2, GUNC, skani, MinCED, tRNAscan-SE 2, Aragorn, TransTermHP, cath-resolve-hits, FAMSA, VeryFastTree, DuckDB, Python environment, `bio-box` CLI and API | wilkcomp02 (Apptainer), Hetzner AX-class (Docker or Apptainer) |
| `ghcr.io/<org>/bio-box:3.0-cuda` | The CPU image plus CUDA runtime, PyTorch/ONNX runtime for ProstT5 encoder inference, Foldseek GPU build | GPU build stages only (GEX130, wilkcomp03 by arrangement) |
| `ghcr.io/<org>/bio-box:3.1-cpu` / `3.1-cuda` (planned, SEARCH_PLAN SD-13) | `3.0` plus the search tools: HH-suite3, BLAST+ (RPS-BLAST for CDD), DIAMOND, CMfinder, R-scape, PADLOC, CRISPRCasTyper, CRISPRCasFinder, and Domainator (MIT; conda package plus its Rust parser extension) as the candidate workbench: `.dsnv` network viewer, `plot_contigs`, structure annotation, profile comparison, `domainator_server` GUI (SD-17). reseek for the SD-21 aligner test. ProstT5 and the ESM-2 3B 3Di encoder stay in the CUDA image only (SD-21 bake-off). | Same hosts as `3.0` |
| `ghcr.io/<org>/bio-box-test:3.1-cpu` (planned, [`planning/tests/TEST_ARCHITECTURE.md`](planning/tests/TEST_ARCHITECTURE.md)) | A layer on `bio-box:3.1-cpu` with test tooling only: pytest, Hypothesis, pandera, coverage.py, mutmut, asv, Zombi with Pyvolve (sequence mode), ALF as the fallback simulator (TEST_ARCHITECTURE §11 item 1). HMMER `hmmemit` and Infernal `cmemit` come from the base image and are not added again. | CI host (R0–R2), wilkcomp02 (R3), stage host (R4–R5); §2.3 |

- Nothing is installed ad hoc on either host ([`planning/data_layer_and_tests.md`](planning/data_layer_and_tests.md) §6). A tool missing from the image needs an image release before anyone can use it.

- Both images are built from one lockfile (pixi or uv lock plus pinned binary digests), so tool versions are identical across hosts. The image digest is recorded in every release manifest.
- **Test tooling is never in the runtime image.** It lives in `bio-box-test`, built from the same lockfile (a `test` feature) on top of the runtime image's digest, so tests exercise the exact binaries a release uses. Releases and runs record the runtime digest; `test_run` rows record the runtime digest the tests ran against (§2.3).
- The v2.0 "4.2 GB image with automatic CPU/GPU fallback" is withdrawn. CPU inference of ProstT5 is not used for bulk work. Its throughput on Haswell has not been measured, and the bulk structural stage is a GPU stage by design.
- Exact image sizes are unmeasured.

### 2.2 CLI subcommands

```bash
bio-box run --profile <host> --stage <Sn> --release <id>   # build or extend a stage through the workflow DAG
bio-box verify --release <id>                                # manifest, ledger and contract checks
bio-box pull-db --release <id> --edition <lite|core|planetary>
bio-box serve --bind 127.0.0.1:8420                          # read-only API (§5)
bio-box append <fasta> --tier <n> --tag <label>              # search-only delta (§3.3)
bio-box register-sidecar <file.parquet> --key-type <type> --key-release <id>   # (§3.4)
bio-box remap-token --merge <old_uid,...> --into <uid>       # merges only (§3.5)
bio-box scan-new-models <models> [--z release_default]        # incremental annotation; Z fixed per SD-04 (§3.6)
bio-box rescore --params params.yaml                          # recompute statistics from stored sufficient statistics (§3.7)
bio-box compact --into-release <id>                           # promote deltas into a new release
bio-box seq get --contig|--slot|--protein|--igr <uid> [--start S --end E --strand +|-] [--flank N] [--nt|--aa]
                                                              # point lookup by content UID (§3.10)
bio-box seq extract --uids uids.txt | --sql "SELECT ..." --out out.fa.gz [--bgzip --faidx] [--members|--reps-only]
                                                              # bulk export, incl. dedup members; requests sorted by (file, row_group)
bio-box seq lineage <uid|accession>                           # slot → contig → source → L0 object → contract

# Planned, search (SEARCH_PLAN §3.3; none exists yet)
bio-box search run <recipe.yaml>                              # submit; returns run_id; async (§4.6)
bio-box search status|log|cancel <run_id>
bio-box search hits <run_id> [--round N] [--members] [--with-context ±N genes|bp]
                    [--json] [--top_n 50] [--format parquet|tsv|domainator-gb]   # SD-17, SD-18
bio-box search ssn <run_id> [--score efi] [--knn 10]         # MSF ∪ kNN + threshold table, .dsnv bundle (SD-19)
bio-box search promote <run_id> --name <family>               # user token overlay (SD-09)
bio-box search stats <run_id|user_token>                      # on-demand anchor statistics against background_rates
bio-box nbr --anchors <run_id|token|uids> --window 10g|20kb --strand any|same --operon-gap 50
            --edge exclude:3kb|truncate|require:5 --circular auto   # neighborhood fetch (SD-08)
            [--json] [--format parquet|domainator-gb]
bio-box nbr score --positive-set <file|token list> [--exclude-anchor-family]
bio-box region translate <contig_uid:start-end> [--six-frame] [--min-orf 30]
bio-box pack build|verify <release>                           # per-release search pack (SD-06)
```

**Output conventions for search and report commands (SEARCH_PLAN SD-18, adapted from Domainator's agent interface).**
- `--json` writes NDJSON to stdout (one object per record) or to a path; it sits alongside the Parquet outputs, never replaces them.
- Lists are capped at query time (`--top_n`, default 50); member lists are returned only when asked for. One stored result serves any cap.
- Schemas are versioned: added keys are compatible, removed or renamed keys are breaking and need a version bump (test S-13). `serve` endpoints and MCP tools return the same JSON.
- `--format domainator-gb` writes GenBank with `Domainator` features and the qualifiers `program`, `database`, `description`, `evalue`, `score`, `name`, `identity`, `cds_id`, `rstart`, `rend`, `rlen`, so results open directly in Domainator's tools. `evalue` is at Z_ref and `score` is in bits; filter exported subsets on `score`, because Domainator's own tools default to `Z = 1000` (SD-17).
- Recipes carry `max_output_gb` (default 25) and `partial: exclude|include|flag` (SD-20); the runner enforces both.

There is no `--gpus` auto-detection. GPU stages run with the CUDA image on hosts that have one.

### 2.3 Where tests run, and where their results go (planned; TEST_ARCHITECTURE §8)

All rings run the `bio-box-test` image (§2.1).

| Ring | Host | Cadence | Envelope |
|---|---|---|---|
| R0 static, R1 unit and property, R2 micro-corpus | CI host: laptop, GitHub Actions or a pod-side runner (open, TEST_ARCHITECTURE §11 item 2) | R0–R1 every commit; R2 every merge; mutation testing weekly | R2 < 30 min (*target*); micro-corpus and synthetic worlds only, no store data. `/tmp`-full and read-only-NFS faults use a size-limited or read-only mount on the CI host, never the shared pod `/tmp` |
| R3 golden | wilkcomp02 | nightly and every image build | inside the `bio-box run-step` envelope: own `MemoryMax` scope, `nice 19`, `ionice -c 2 -n 7`, shards < 300 CPU-min each, STOP file honored; threads come out of the night build share (48), not added to it; < 6 h (*target*) |
| R4 stage gate | stage host (wilkcomp02 or Hetzner, per `STAGES.md`); search tests on the search host: pod lane through S5, Hetzner from S6 (`planning/data_layer_and_tests.md` §4.10) | per stage | as the stage's build steps (`STAGES.md`); multi-client tests (K-* with 1–16 DuckDB clients, SLO load) at night only (`planning/preflight_gaps.md` #20) |
| R5 calibration, benchmarks, replays | stage host (pod or Hetzner) | per tier | as SEARCH_PLAN §5–6; null-world calibration at stage scale is budgeted in `planning/sizing.yaml` like any stage step |

**Results.** Every test session writes `test_run` rows (`test_id`, ring, operation, workflow, params hash, outcome, metrics JSON, duration, code commit, `runtime_image_digest`, `test_image_digest`, `release_id`, host, seed) as packed Parquet under `$ATD_ROOT/tests/` (`$ATD_ROOT` is the data root, §3.9), one file per session through `run-step` local staging (one `rsync`, atomic rename) and compacted per stage, so the NFS small-file rule holds. The single writer registers the files in `manifest.json`, and the generated catalog exposes `test_run` as a view (§3.8). A stage gate is a SQL query over that view (every blocking test for the stage passes on the stage's `release_id` and image digest); a regression is a diff between two sessions. R0–R2 rows from an off-pod CI host reach `$ATD_ROOT/tests/` through the same writer; the transport depends on the CI-host choice.

---

## 3. Malleability: what can change without a rebuild, and what becomes stale (D-19, D-20, D-21)

v2.0 claimed full malleability. That claim was wrong. Appends, family remaps and compaction invalidate derived global statistics, and they renumbered the IDs that sidecars were keyed by. v3.0 replaces the claim with a contract.

### 3.1 The contract

1. **Raw content accretes.** New contigs, proteins and annotations are searchable as soon as their delta is registered.
2. **Statistics are versioned.** Every statistic carries `stats_release_id`. It is correct "as of" that release and is flagged `stale=true` for any key whose inputs changed after it.
3. **Identifiers are stable.** User-facing keys never change meaning within or across releases. Changes are expressed through a crosswalk.
4. **Nothing silently re-keys.** Any operation that would change the meaning of an existing key is refused or routed through a new release with a crosswalk.

### 3.2 Identifiers

| ID | Definition | Stability |
|---|---|---|
| `contig_uid` | Digest of the normalized nucleotide sequence (BLAKE3, truncated to 128 bits) | Stable forever |
| `protein_uid` | Digest of the amino-acid sequence | Stable forever |
| `source_record_uid` | Digest of `(source, accession, version)` | Stable forever |
| `cluster_uid` (NR90, NR50, subfamily) | `blake3(sorted member UIDs at creation, cluster_params_hash)`, 128 bits ([`planning/data_layer_and_tests.md`](planning/data_layer_and_tests.md) §2). The representative's UID is stored as an attribute, not used as the key | Stable until a split or merge; changes recorded in `id_history` |
| `token_uid` | Namespaced `(level, library, accession, library_version)` | Stable; library upgrades produce new UIDs with crosswalk rows |
| `dense_contig_id`, `dense_slot_id` | Release-local integers for compact storage and joins | Valid only inside one `release_id`. Never exposed as a sidecar key. |

- `id_history.parquet` records, per release transition, each cluster or token as `same`, `split`, `merged` or `retired`, with its successor(s).
- Every table row carries `release_id`. Every statistic carries `stats_release_id`.

### 3.3 Mechanism 1: search-only deltas (`bio-box append`)

1. The base release is read-only.
2. `append` runs the same per-genome pipeline as the base build (identity resolution, screening, gene calling, annotation; [`STAGES.md`](STAGES.md) layer L1). It deduplicates against the base by `contig_uid` and containment before allocating anything.
3. New proteins are assigned to existing clusters with `mmseqs clusterupdate` semantics where possible. Unassigned proteins form provisional clusters flagged `provisional=true`.
4. New sequences, including collapsed duplicates, are appended to the sequence store (§3.10). The delta itself is written as new Parquet partitions under `epochs/<release_id>/delta_<n>/`, plus a manifest entry. Dense IDs in a delta are allocated after the base's maximum. If a bitmap index exists, epoch window ranges start on multiples of 65,536 so containers never straddle epochs.
5. **Stale after `append`:** NR50 and subfamily membership of affected clusters, dedup collapse counts, vehicle clusters, entropies, Δbp distributions, coupling statistics and q-values, pangenome planes, transition edges, and any model trained on the base. Each affected key is flagged `stale=true` until the next `compact` or `rescore`.
6. Pair counts are additive and are merged on the fly. Medians, entropies, test statistics and q-values are not additive and are never "added" across epochs.
7. **Ad hoc searches (§4.6) include appended and private epochs only when the recipe asks for them.** The epochs searched are recorded in `search_run`. Z does not change when they are included: it is the release-line constant (SEARCH_PLAN SD-04), and each hit stores the real target size beside the E-value at that Z.

### 3.4 Mechanism 2: sidecar tables (`bio-box register-sidecar`)

1. A sidecar is a Parquet file keyed by a declared `key_type` (`contig_uid`, `protein_uid`, `cluster_uid`, `token_uid`, `(token_uid_a, token_uid_b)`, or `(contig_uid, slot_ordinal)`) and a declared `key_release`.
2. Registration checks the schema, checks key coverage against the target release (warns below 95% coverage, refuses below 50%), and refuses keys from a release other than the target unless a crosswalk is supplied. **Coverage** is the fraction of the sidecar's keys that resolve in the release, not the fraction of the release's keys present in the sidecar, so sparse hit lists are accepted (SEARCH_PLAN §3.1).
3. Sidecars are joined after the main filter, never before. A join builds a hash table on the sidecar side. For a 2 GB sidecar that takes seconds, not "< 1 ms". Large slot-level sidecars (up to the number of slots in the release) are stored sorted by key and probed rather than hash-joined. Join cost is measured per query class in the test suite.
4. Defense, secretion and conjugation system calls (DefenseFinder, PADLOC, CONJScan, TXSScan) are stored as sidecars of this kind, not as tokens (D-09).

### 3.5 Mechanism 3: token remapping, merges only (`bio-box remap-token`)

- A remap that **merges** tokens (for example, a DUF recognized as a member of an existing clan) writes `token_remap_overlay.parquet` and an `id_history` row. Additive statistics (counts) for the merged token are recomputed by summation. Non-additive statistics are flagged `stale`.
- A **split** cannot be done by overlay. It requires re-assigning every member slot. Splits go through `compact` into a new release, with `id_history` recording `split`.
- The v2.0 claim of "< 200 ms" remaps that keep all statistics valid is withdrawn.

### 3.6 Mechanism 4: incremental model scans (`bio-box scan-new-models`)

- Scan unit per SEARCH_PLAN SD-05: exact-unique proteins through S2C; from S3 on, new HMMs are scanned against the stored NR90 representatives and hits are projected to members through `cluster_member_aln` (member-to-representative alignments, L2). Members whose projected bitscore is within 10 bits of the cutoff, or whose alignment covers < 90% of the domain envelope, are re-scored exactly. Test S-08 checks projection.
- `Z` is one constant across the release line (SEARCH_PLAN SD-04), not a per-release value: HMMER `Z = domZ = 1e8` sequences; Infernal and nhmmer `Z = 1e5` Mb. It is passed explicitly. Each hit stores the bitscore, the E-value at that Z, and the real target size, so E-values from incremental scans stay comparable to the base build and across stages (S-03). The constants await approval (SEARCH_PLAN §9 item 2).
- The v2.0 "45 seconds" figure is unmeasured and likely optimistic by 4–20× (AR audit, cell-count arithmetic). The runtime is measured per stage.
- Output is a sidecar or delta partition with full provenance (library, version, cutoff type, `Z`).

### 3.7 Mechanism 5: statistical rescoring (`bio-box rescore`)

- Coupling statistics depend on parameters (`params.yaml`: thresholds, priors, fitted weights; D-13). v2.0 baked them into the pair summary, which made `lite` and `core` statistics permanently frozen.
- Each release stores per-pair **sufficient statistics**: Σp₀, Σp₀q₀ and Σy per regime × tier stratum; Δbp histograms; vehicle and BioProject counts; and neighbor sketches. With these stored, `bio-box rescore --params new.yaml` recomputes `Z`, effect sizes and FDR without the base contig tapes.
- A rescore writes a new `stats_release_id`. The params hash, seed and code commit are recorded in the manifest.

### 3.8 Catalog: read-only and generated (D-21)

- DuckDB allows one read-write process or many read-only processes per database file. The v2.0 design, where `serve`, `append` and `register-sidecar` all opened `catalog.duckdb` read-write, was not possible.
- v3.0: the source of truth is `manifest.json` (partitions, schemas, checksums, `release_id`, sidecar registrations). A single writer process (the orchestrator, or `bio-box` admin commands holding a lock file) updates the manifest. It then generates a fresh `catalog.duckdb` containing only views and macros over the Parquet partitions, and swaps it in with an atomic `rename`.
- `serve` and interactive users open the catalog read-only. A swap takes effect on their next connection.
- Whether DuckLake replaces this generated catalog is evaluated in [`planning/data_layer_and_tests.md`](planning/data_layer_and_tests.md).

### 3.9 Directory layout

```
<data-root>/                              # /stor/scratch/Wilke/if2453/all-the-data on wilkcomp02; /mnt/nvme/all-the-data on Hetzner
├── manifest.json                         # system-of-record index of everything below (atomic rename)
├── catalog.duckdb                        # generated, read-only views/macros
├── sources.yaml                          # per-source contracts, licenses, redistribution flags (§6)
├── LICENSES/  CITATION.cff
├── l0/refdb/<name>/<version>/            # reference DBs: pressed HMMs, calibrated CMs, Foldseek DBs, weights; sha256, treated as L0 sources
├── l0seq/                                # raw sequence store (§3.10); append-only, shared by all releases
│   ├── manifest.json                     # durable copy also kept on /stor/work on wilkcomp02
│   ├── seq_chunks/stage=<S>/digest_prefix=<00..ff>/part-<k>.parquet       # (raw_digest, chunk_ix, seq), 64 KB chunks, incl. all dedup members
│   ├── contig_occurrence/stage=<S>/source=<name>/part-<k>.parquet         # one row per source record
│   ├── protein_seq/gene_call_version=<v>/stage=<S>/protein_prefix=<00..ff>/part-<k>.parquet
│   ├── igr/igr_version=<v>/stage=<S>/part-<k>.parquet                     # coordinates only; bases sliced from seq_chunks
│   ├── seq_locator.parquet               # raw_digest → file, row_group, first_chunk_row; rebuilt per stage
│   └── recipes/<source>.yaml             # fetch recipes for sequences not redistributable (§6)
├── releases/<release_id>/
│   ├── L0_sources/  L1_annotations/  L2_clusters/  L3_tapes/  L4_stats/
│   ├── id_history.parquet  ledger.parquet  release_manifest.json
│   └── epochs/delta_<n>/                 # search-only deltas after the release
├── search_pack/release=<R>/              # planned: per-release search targets, ≥256 MB shards + pack_manifest.json (§4.6; SEARCH_PLAN §3.2)
├── search/run=<id>/                      # planned: packed outputs of one search run (run, hits, profiles, neighborhoods; optional ssn); search/run_index.json
├── tests/                                # planned: test_run results, packed Parquet, registered in the catalog (§2.3)
└── sidecars/<name>/<key_release>/*.parquet
```

This tree is the edition (packaging) view. The physical layout on the pod and on Hetzner is `planning/data_layer_and_tests.md` §3.1 (`l0/`, `l1/annot_version=…`, `l2/release=…`, `l3/`, `l4/release=…`, `search_pack/`, `search/`), which is canonical; `releases/<release_id>/` is the per-release manifest view over those directories, not a second copy.

Tool-interop FASTA (bgzip + `.fai` + `.gzi`) is not part of the store. `bio-box seq extract` generates it on demand into `/tmp` or a project directory. Layout follows [`planning/data_layer_and_tests.md`](planning/data_layer_and_tests.md) §2A.1.

### 3.10 Sequence store: raw sequences stay available to every project

Requirement: the nucleotide sequence of every ingested contig, including members collapsed by deduplication, stays retrievable, along with every protein and IGR sequence. Any downstream project can then re-annotate, re-cluster or re-assemble without going back to the original sources.

- **Keying.** Contig bytes are stored once per `raw_digest` (sha256 of the exact bytes) and mapped to `contig_uid` (§3.2); proteins once per `protein_uid`. IGRs are coordinates `(contig_uid, start, end, strand, left_slot, right_slot, igr_version)`, not a second copy of the bases. The same sequence arriving from several sources is stored once. `contig_occurrence` records every source occurrence, and `contig_member` records dedup membership. Dedup decides what is *counted*, never what is *kept*.
- **Format.** The primary store is Parquet: `seq_chunks` holds 65,536-byte chunks sorted by `(raw_digest, chunk_ix)`, zstd, about 1,024 chunks per row group; `protein_seq` is sorted by `protein_uid`. `seq_locator.parquet` maps each `raw_digest` to its file and row group, so a slice reads one row group. bgzip FASTA with `.fai`/`.gzi` is generated on demand for tools that need FASTA; zstd-seekable is not used (no samtools/pysam support). The S0 bake-off (test K-13) decides whether Parquet chunks or bgzip shards are primary.
- **Immutability.** The store is append-only and shared across releases. A release references sequences by UID and never copies them. Sequences are removed only for licensing or human-read removal, and such removals are logged in the manifest.
- **Access.** `bio-box seq get` handles point lookups. `bio-box seq extract` handles bulk export by UID list or `--sql` query, including dedup members and region slices. DuckDB views (`v_seq_contig`, `v_seq_protein`, `v_seq_igr`) join sequences to any table by UID.
- **Integrity.** On write and in `bio-box verify`, the digest of each stored sequence must equal its UID. The test suite includes a round-trip check: extract, re-digest, compare.
- **Placement on wilkcomp02.** Shards live on `/stor/scratch/Wilke/if2453/all-the-data/l0seq/`, where they are regenerable from L0 snapshots. `l0seq/manifest.json` and per-shard checksums are copied to `/stor/work/Wilke/if2453/all-the-data/` so the store can be rebuilt and verified. On Hetzner, object storage is the system of record (§4.2).
- **Size.** Estimate, pilot-pending. It scales with total ingested bp after exact-sequence dedup, not after architecture dedup.

---

## 4. Deployment target 2 (ultimate): rented bare metal (Hetzner)

### 4.1 Hardware facts (verified by the round-2 architecture audit; prices and options change, so re-check at order time)

| Role | Product | Facts | Notes |
|---|---|---|---|
| CPU build and query node | AX162-R / AX162-S | AMD EPYC 9454P, 48 cores / 96 threads, Zen 4 with AVX-512; 256 GB DDR5 ECC (more optional); base storage 2× 1.92 TB NVMe, larger NVMe drives as paid options | v2.0 cited an EPYC 9654 and 2× 7.68 TB as standard. Both were wrong. |
| GPU build node | GEX130 | Xeon Gold 5412U (24 cores); 128 GB RAM; 1× RTX 6000 Ada 48 GB | v2.0 said 256 GB RAM. RTX 6000 Ada's 91 TFLOPS is the FP32 non-tensor figure. ProstT5 time at the 18M-segment scale is estimated at 24–39 GPU-hours, unmeasured. |
| Network | Dedicated server uplink | 1 Gbps by default; 10 Gbps uplink is a paid add-on with a traffic allowance | All v2.0 download times assumed 10 Gbps. |

### 4.2 Storage roles

- **Object storage is the system of record:** Hetzner Object Storage, Cloudflare R2 or a sponsored bucket. Releases are published there as content-addressed chunks plus `manifest.json`.
- **Local NVMe RAID-0 is scratch and cache.** A single-drive failure loses the array. Recovery means re-pulling from object storage, never restoring from a peer. User-private epochs are synced to object storage on every change, because RAID-0 is not a backup.
- The query node does not host the public download.

### 4.3 Provisioning

1. `deploy/hetzner-bootstrap.sh`:
   - Install Ubuntu 24.04 via `installimage`.
   - `mdadm` RAID-0 across the data NVMe drives, XFS, mounted `noatime`.
   - `vm.max_map_count` raised; transparent hugepages set to `madvise`.
   - Install Docker or Apptainer; the same image digest as wilkcomp02.
   - node-exporter, SMART alerting and Prometheus scraping.
2. `bio-box pull-db --release <id> --edition <e>`:
   - Parallel chunk download (`s5cmd`); BLAKE3 verification per chunk; unpack into a staging directory.
   - Atomic swap of the `current` symlink.
   - Changed chunks only on upgrade. After any release that renumbers dense IDs, most chunks change, so delta downloads help less than v2.0 claimed.
3. `bio-box verify --release <id>`: manifest, ledger and contract checks before `serve` starts.

### 4.4 Download time (derived: `size × 8 / rate`, at full line rate; sizes are estimates, pilot-pending)

| Edition | Size without search pack (estimate) | Search pack (estimate, derived) | Total | 1 Gbps | 10 Gbps |
|---|---|---|---|---|---|
| `lite` (S0–S5) | ~0.3–1.0 TB (sequence store ~0.2 TB + derived tables; v2.0's 85 GB is withdrawn) | ~0.12–0.19 TB (S5 pack subset, see note) | ~0.4–1.2 TB | ~55 min – 2.7 h | ~5–16 min |
| `core` (S0–S6) | ~4–6 TB | ~0.8–1.0 TB (S6 pack, see note) | ~4.8–7 TB | ~11–16 h | ~64–93 min |
| `planetary` (S7+) | unknown until the S7 rarefaction check | ST5 is reduction-only; no pack sized | — | — | — |

Search-pack note (*derived* from SEARCH_PLAN §3.2 rows; the zstd ratios ~0.6 B/aa and ~0.27 B/nt are the S5 NR50 and IGR rows' own ratios, the MMseqs2 index ~1.07 B/aa likewise):
- `lite` subset: NR50 FASTA 20–35 GB + NR50 index 35–60 GB + IGR 14 GB + dark ORF ~18 GB (30 Gaa × 0.6) + `hh/` 30–60 GB + spacers and models (small) ≈ 117–187 GB. The full S5 pack, including `nr90`, `unique.dmnd` and `segments`, is ~370–530 GB (*estimate*); the rest is not shipped in `lite`.
- `core` S6 pack: NR50 24–48 GB + index 43–86 GB + NR90 69–120 GB + `unique.dmnd` ~330 GB + IGR 87–100 GB + dark ORF ~120 GB + `3di_nr50` 60–120 GB + `hh/` 30–60 GB ≈ 0.76–0.98 TB.

Real throughput depends on the source bucket's egress and on the server's disk write rate. Measure it during the first rehearsal (see [`STAGES.md`](STAGES.md)).

### 4.5 Start-up

- Parquet is read through DuckDB's buffer manager and page cache. It is not "zero-copy". Bitmap indexes, if built, are zero-copy only in CRoaring's frozen serialization format.
- `madvise(MADV_WILLNEED)` is a hint. Warming a few GB of index files takes seconds. Everything else starts cold, and first-query latency reflects that.
- `serve` is ready when `verify` passes and the catalog is attached. Cold and warm latencies are reported separately in the test suite.

### 4.6 Search on the pod and on Hetzner (planned; SEARCH_PLAN SD-01, SD-06, SD-15)

Searches are recorded jobs, not queries. `bio-box search run` submits a recipe and returns a `run_id`; the runner executes it asynchronously, outside `serve` and its query timeouts (§5). `serve` can submit a run and poll its status and hits; it never runs a search inline. Runtimes are minutes to hours (SEARCH_PLAN §5).

**Interactive lane on wilkcomp02 (SD-15).** The §0.2 thread caps are split between the build and one search lane, not raised:

| Period | Build threads | Search lane threads | Total |
|---|---|---|---|
| Day | 32 | 8 | 40 |
| Night and weekends | 48 | 8 | 56 |

- Search jobs run as per-process shards under 300 CPU-min each (watchdog limit ~369 CPU-min), at `nice 15`, inside their own `systemd-run` scope.
- `/tmp` quota for the lane: 150 GB. A run copies the search-pack shards it needs from NFS to `/tmp`; the S5 per-run working set is 60–120 GB (NR50 plus IGRs, *estimate*).
- During L2 barriers the lane reads shards from NFS instead of `/tmp`.
- **Per-run outputs are packed.** Each run writes `search_run`, `search_hit`, `search_profile` and its neighborhoods (plus, optionally, the SSN table and `.dsnv` bundle) under `search/run=<id>/` on NFS with one `rsync` (§0.2 I/O rule; file list in `planning/data_layer_and_tests.md` §2B.3). Exports (`--format domainator-gb`, `tsv`, `--json`) go to a project directory or `/tmp`, never into the run directory. No per-hit, per-round or per-query loose files. Small-file waiver: run files are often < 64 MB, so the NFS rule is waived for at most 6 files per run (run, hits, profiles, neighborhoods, and the optional `ssn.parquet` and `ssn.dsnv`; ~10^3 runs per year, *estimate*, keeping `search/` under 10^4 files), and for one overlay file per promoted family (`planning/data_layer_and_tests.md` §2B.3).
- **Search pack (SD-06).** `bio-box pack build` writes each release's targets once, as ≥256 MB shards on NFS with a manifest. S5 pack: ~370–530 GB (*estimate*, row sum in `planning/sources_and_sizing.md`) against 30 TB of scratch.

**Hetzner (tier ST4 and later).** From S4/S6 the search pack exceeds the pod. ST4 needs a dedicated search host on Hetzner: ≥512 GB RAM and ≥4 TB NVMe (*estimate*; size and timing are an open item, SEARCH_PLAN §9 item 3). The base AX162 (256 GB RAM, 2× 1.92 TB NVMe, §4.1) therefore needs its RAM and NVMe options. `unique.dmnd` at S6 (~330 GB) and DIAMOND blastx over all contigs (~55 h, *derived*) run there as scheduled jobs.

**SQL console and `enable_external_access=false` (§5).** The optional SQL console cannot open files outside the release directory, so it cannot read `search/run=<id>/` by path. When a run closes, the runner adds it to `search/run_index.json` by atomic rename (`planning/data_layer_and_tests.md` §2B.3), and the regenerated catalog exposes the runs listed there as views over explicit file lists (§3.8). The console reads registered run outputs through the catalog only; a new run becomes visible on the next connection after the swap. Outputs of runs that include private epochs follow the private-epoch rule in §5.

---

## 5. Security and API (D-22)

| Threat (STRIDE) | Control |
|---|---|
| Spoofing / unauthorized access | Bind to `127.0.0.1`. Remote access via SSH tunnel (wilkcomp02) or a TLS reverse proxy with token auth (Hetzner). Admin commands are local CLI only, never exposed over HTTP. |
| Tampering with data packs | BLAKE3 per chunk, signed `release_manifest.json`, read-only mounts for releases. |
| Information disclosure (private epochs) | Private epochs live in a separate directory with their own manifest, are never included in published packs, and are excluded from `serve` unless explicitly enabled. |
| Denial of service | Query timeouts, row limits and per-client concurrency limits. DuckDB `memory_limit` and `threads` are set below host caps. Searches are async jobs on the search lane, outside these timeouts; `serve` only submits and polls them (§4.6). |
| Elevation via SQL | No raw SQL endpoint by default. The API exposes parameterized query classes (`/v1/...`). An optional SQL console opens DuckDB read-only with `enable_external_access=false`, `lock_configuration=true`, and file access limited to the release directory. |
| Running as root | Containers run as an unprivileged user. Apptainer is rootless by default. |

Every API response carries `release_id`, `stats_release_id`, the query-class version and a query hash.

---

## 6. Licensing and redistribution (D-23)

- `sources.yaml` holds one contract per source: name, URL, release/version, retrieval date, checksum, license text reference, a `redistribution` flag (`yes`, `no`, `attribution`, `share-alike`, `unknown`) and citation.
- **Public data packs contain only sources whose `redistribution` flag permits it.** Sources that cannot be redistributed ship as bring-your-own-data recipes: a script that fetches the data under the user's own agreement and builds the matching partitions locally.
- Every pack includes `LICENSES/` and `CITATION.cff`.
- Known constraints, to be confirmed against current source terms before any public release:

| Source | Constraint |
|---|---|
| ISfinder | No bulk download without written authorization. Use ISEScan models, or obtain permission. |
| JGI IMG/M, IMG/VR, IMG/PR | JGI data policy applies. Legacy-data restrictions and PI contact requirements mean these are bring-your-own by default. |
| GSA-Human (CNCB) | Controlled access. Do not ingest. |
| GTDB | Check the license terms, possibly share-alike, before redistributing derived taxonomy. |
| Logan | No license stated; NCBI/NIH data-sharing policy applies to the underlying runs. Human-read removal is mandatory before anything derived is published. |
| ICEberg, TnCentral, PaperBLAST, Fitness Browser | No explicit license found. Treat as `unknown` until clarified. |
| AllTheBacteria | MIT. |
| HMM / CM libraries | Per-library terms recorded in `l0/refdb/<name>/<version>/LICENSE`. |

---

## 7. Editions are outputs of stages

Editions are no longer designed top-down. Each is a named cut of the releases produced by the stages in [`STAGES.md`](STAGES.md). The mapping below follows the stage ladder in that file. If the two disagree, [`STAGES.md`](STAGES.md) is authoritative.

| Edition | Built from stages | Contents | Size | Runs on |
|---|---|---|---|---|
| `lite` | S0–S5 (golden genomes, complete phages, complete plasmids, complete prokaryotic genomes (S2C), small genomes, IMG/VR HQ as statistics only, all GTDB R226 species representatives) | Contig tapes for those sources, cluster dictionaries, pair statistics with sufficient statistics, core plane over species reps; S5 search pack (NR50, IGR, dark ORF, spacers, models, `hh/`) | Estimate, pilot-pending; search pack ~0.12–0.19 TB (§4.4) | Laptop, wilkcomp02, small VM |
| `core` | S0–S6 (adds all GTDB R226 genomes and the RefSeq complete flag) | Adds multi-strain pangenome planes and quality masks; S6 search pack | Estimate, pilot-pending; search pack ~0.8–1.0 TB (§4.4) | Hetzner AX-class; search on the ST4 search host (§4.6); wilkcomp02 for per-genome work only (§0.4) |
| `planetary` | S7+ (Tier 3 metagenomic contigs; Logan built compute-to-data on AWS) | Adds unbinned and SRA-derived contigs | Estimate, pilot-pending | Hetzner AX-class with larger NVMe |

Every edition includes the sequence store (§3.10) for its stages: contigs (including dedup members), proteins and IGRs. Where a source is not redistributable (§6), the edition carries that source's fetch recipe (`l0seq/recipes/<source>.yaml`, listing UIDs, accessions and the build command) in place of the sequences. `bio-box seq` then reports those UIDs as `fetch-required` until the user runs the recipe.

Search in `lite`: IMG/VR (S4) sequences are not redistributable, so they are absent from the shipped search pack. A search over `lite` reports S4 targets as `fetch-required` coverage in its summary rather than silently missing them. `core` ships the S6 pack, which runs on Hetzner (§4.6).

Each edition carries the release manifest of the stages it came from. An edition can be upgraded in place only when the target release's `id_history` maps every key in the installed release.
