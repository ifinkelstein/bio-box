# Data layer and test plan for staged accretion

Status: plan only. Nothing here has been run. Every number marked *(measure)* is a placeholder until the stage that measures it.

Scope: how data accretes stage by stage without rewriting earlier stages, and which tests gate each stage. The deployment target is `wilkcomp02` (CPU-only, shared, no scheduler). The ultimate target is a rented bare-metal box (Hetzner). The container and CLI are named `bio-box`.

Inputs: audits `critique_4` through `critique_8` (stable IDs, deterministic dedup, stale-statistics contract, single-writer DuckDB, sizing ledger, golden genomes) and `WILKECOMP.org` (pod envelope).

---

## 0. Operating envelope this design must fit

| constraint | value on wilkcomp02 | design consequence |
|---|---|---|
| CPU | 56 physical Haswell cores, AVX2, 2.2 GHz; 40 threads by day, 56 at night/weekend | thread cap is a profile variable, never hard-coded; AVX2 builds only |
| memory | 1,007 GB physical, but a 200 GB `systemd-run --user --scope -p MemoryMax=` cgroup | DuckDB `memory_limit` ≤ 150 GB inside a 200 GB scope; MMseqs2 query chunking |
| address space | `ulimit -v` 500 GB soft | never mmap a target DB larger than ~400 GB in one process |
| local disk | `/tmp` ~940 GB ext4, volatile, not shared | all hot work (MMseqs tmp, DuckDB `temp_directory`, shard assembly) on `/tmp`; never the only copy |
| shared disk | `/stor/scratch/Wilke` 30 TB NFS over 10 GbE, shared by 3 hosts, not backed up; `/stor/work` backed up | NFS holds only large immutable files (≥64 MB); one rsync per finished shard |
| small files | banned on NFS | no per-genome files on NFS; per-genome outputs are packed into shard Parquet on `/tmp` first |
| watchdog | renices any single PID after ~369 CPU-minutes | monolithic multi-threaded tools launch at `nice 19`; per-genome workers are short-lived at `nice 15` |
| open files | 4,096 hard (raise to 8,192 with `ulimit -n`) | cap concurrent Parquet writers; MMseqs split counts |
| scheduler | none | the orchestrator itself enforces caps, STOP file and load checks |

---

## 1. Layered architecture

Each layer reads only from layers below it and writes immutable, content-hashed files. Higher layers can be recomputed from lower ones, so a bug in L4 never forces re-running L1.

| layer | content | update mode | recompute trigger | owner of truth |
|---|---|---|---|---|
| **L0 Source snapshots + contracts** | Raw downloaded FASTA/GFF/metadata exactly as fetched; `sources.yaml` contract per source (accession list, release, retrieval date, sha256, license, redistribution flag) | **Append-only.** A new source release is a new snapshot directory, never an overwrite | Never recomputed. Superseded snapshots are retired by manifest, not deleted, until the next release is signed off | sha256 of every object, recorded in `l0_manifest.parquet` |
| **L0-seq Raw sequence store** (§2A) | Every source contig's nucleotide sequence, including every dedup member, stored once per distinct byte sequence, with an occurrence row per source record. Protein sequences for every called protein. IGRs and slots as coordinates, sliced on demand | **Append-only.** A new stage adds new shards and occurrence rows; existing shards are never rewritten | Never recomputed. A new gene-calling version adds a new protein table version; nucleotide shards are untouched | `raw_digest` (sha256 of exact bytes) per sequence; `l0_manifest` sha256 per source object |
| **L0.5 Identity resolution** | `genome_identity.parquet`: assembly UID, BioSample, GCA↔GCF links, sequence digest, containment-cluster representative, precedence rank, `member_of` | **Incremental**: a new genome is matched against existing reps; a new rep never displaces an old one within a release | Change in the precedence policy, or a new source that contains existing genomes (then rerun cross-source matching for the affected clusters only) | Deterministic rule set (`identity_policy.yaml` hash) |
| **L1 Per-genome annotation** | Gene calls, smORF candidates, ncRNA/CRISPR/tRNA hits, HMM domain hits with per-library cutoffs (scan unit per §2B.5), MGE/regime calls (geNomad etc.), QC (CheckM2/GUNC), taxonomy. One row-set per genome, packed into shards | **Append-only per genome.** A genome's L1 output depends only on its own sequence plus pinned tool and DB versions (projected domain hits from S3 on are release-scoped, §2B.5) | A tool or reference DB version bump creates a new `annot_version`, with L1 recomputed for all genomes under that version. The old version stays readable until the release is cut | `(genome_uid, annot_version)` |
| **L1-derived `protein_occurrence`** (§2B.1) | Protein-first index: one row per slot that carries a protein, with genome, contig, coordinates, strand, taxonomy, species-rep flag, tier, assembly quality, distance to contig end. Serves `protein_uid → slots` and the locus context on every search hit | **Appended per stage**: one partition per `(gene_call_version, stage)`; earlier partitions untouched | Built at the end of each stage's L1 from `slot_protein` + slot coordinates + L0.5; a new `gene_call_version` builds a new partition set | `(slot_uid, gene_call_version)` |
| **L2 Global clustering** | NR90/NR50 protein clusters, domain-envelope subfamilies, IGR clusters, plus `cluster_member.parquet`, `cluster_member_aln` (member → representative alignment per level, §2B.1) and `id_history.parquet` crosswalks | **Incremental** within a release line (`mmseqs clusterupdate`, MMseqs2 ≥ 18-8cc5c for the issue #961 fix); **full recluster** at major releases. New members get alignment rows at `clusterupdate` | (a) cumulative new sequences > 20% of the existing set *(tunable)*; (b) clustering-parameter change; (c) a drift test fails (§4.6). Every recluster emits a crosswalk old→new with `same/split/merged/retired` and regenerates `cluster_member_aln` for the affected clusters | Cluster UID = `blake3(sorted member UIDs at creation, cluster_params_hash)` (§2); the representative UID is stored as an attribute, not just the cluster index |
| **L3 Contig tapes + slot pairs** | `contig_tapes` (one row per contig; slot lists) and `slot_pairs` (anchor, partner, offset, topology; PLAN §4.2), keyed by stable UIDs and joined to L2 tokens | **Append-only for tapes** (a new genome adds contigs). **Token columns are recomputed** when L2 changes: tapes store `protein_uid`/`segment_uid` (stable), and token IDs are a release-time join | New genomes append tape shards. An L2 recluster regenerates only the token-projection columns (`tape_tokens_<release>`), never the coordinate tape | `(contig_uid)` for coordinates; `(contig_uid, release_id)` for tokens |
| **L3 `user_token_overlay`** (§2B.1) | Sparse rows: slots hit by a family promoted from a search run (`bio-box search promote`, token namespace `source_db = user`), with segment envelope and `run_id` | **Append-only per promotion**; never edits the tapes | Re-projected only when a new `gene_call_version` changes slots (rerun the recipe on the new release). Read by on-demand statistics only; it enters release L3/L4 when the model is registered in `user_models.yaml` and scanned at the next release, so the layer order holds within a release | `(slot_uid, user_token_id, run_id)` |
| **L4 Statistics** | Pair summaries, conditional probabilities, entropies, vehicle clusters, LLR, decoy-calibrated q-values, pangenome planes, tree transitions, closure cards, `domain_promiscuity` and `background_rates` (§2B.1) | **Recomputed per release**; mergeable counts may be updated incrementally between releases, but the published values are a full recompute | Any change in L1–L3 inside the release scope. Between releases, deltas are marked `stale=true` and carry `stats_release_id` of the last full computation | `stats_release_id` |
| **L5 Serving indexes (optional)** | Roaring postings, offset-indexed tape store, USearch, frozen catalogs | **Rebuilt per release** from L3/L4 | Only built if a measured SLO in §4.8 fails with DuckDB-over-Parquet alone. Rebuilt whenever L3/L4 for that release changes | Derived; disposable |
| **Search pack (per release, derived)** (§2B.2) | Search targets as ≥256 MB shards plus `pack_manifest.json`: NR50, NR90, exact-unique, segments, IGR, IGR90, dark ORFs, spacers, models, HH databases, 3Di of NR50 (from ST3) | **Built once per release**, after `barrier_L2_<release>` and before the release cut; immutable | A new release (any L0-seq, L1 or L2 change in the release scope). Never patched between releases | `pack_id` (§2) |
| **Search runs** (§2B.3) | `search_run`, `search_hit`, `search_profile`, plus the run's neighborhoods; one packed directory per run | **Append-only per run**; a run directory is never edited after it closes | Never recomputed. A rerun on a later release is a new `run_id` linked by `recipe_hash` | `run_id` (§2) |

Rules that make accretion safe:
1. **A layer never mutates a file it has published.** Changes are new files plus a new manifest.
2. **Coordinates are stable; annotations and tokens are versioned.** Things that change when knowledge changes (family assignment, cluster membership, regime calls, statistics) live in columns keyed by `release_id` or `annot_version`, separate from things that do not change (sequence, coordinates).
3. **Global steps are explicit barriers.** L2 and L4 never run on a partial L1, and the DAG names these barriers (`barrier_L2_<release>`, `barrier_L4_<release>`).

---

## 2. Identifier scheme

| ID | definition | stability | width |
|---|---|---|---|
| `source_uid` | `sha256(source_name, source_release)[:8]` | permanent | uint64 |
| `genome_uid` | `blake3(sorted contig sequence digests)`, 128 bits | permanent; identical assemblies from two sources collapse | 16 bytes |
| `contig_uid` | `blake3(canonical uppercase sequence)`, 128 bits, with the source accession kept as an attribute | permanent | 16 bytes |
| `protein_uid` | `blake3(translated AA sequence, including U/O and frameshift-stitched sequence)`, 128 bits | permanent | 16 bytes |
| `segment_uid` | `blake3(protein_uid, env_start, env_end, segmentation_version)`, 128 bits | stable per segmentation version | 16 bytes |
| `slot_uid` | `(contig_uid, slot_ordinal_in_contig)` | stable per gene-calling version | composite |
| `cluster_uid` | `blake3(sorted member UIDs at creation, cluster_params_hash)`, 128 bits. The representative's `protein_uid` (or `segment_uid`) is stored as an attribute, never as the key | stable until a split or merge, recorded in `id_history` | 16 bytes |
| `token_id` | namespaced: `level(3 bits) | source_db(5 bits) | local(24 bits)`; for Pfam/Rfam the accession is the key and `local` comes from a frozen dictionary | stable for external accessions; release-scoped for dark clusters (resolved through `cluster_uid`) | uint32 |
| `user_token_id` | a `token_id` in namespace `source_db = user`; `local` is assigned at `bio-box search promote` and registered in `l3/user_overlay/user_tokens.parquet` (`user_token_id`, name, `run_id`, hit-set hash, owner) | stable across releases; becomes a library token (new namespace) when the model is registered in `user_models.yaml` | uint32 |
| `region_uid` | `blake3(contig_uid, start, end, strand)`, 128 bits; names a nucleotide search target (IGR, dark ORF, contig region, spacer) in `search_hit` | permanent (coordinates are stable) | 16 bytes |
| `run_id` | `blake3(recipe_hash, release_id manifest hash, image digest, submit time in ns, host)`, 128 bits; printed as hex | permanent; a rerun of the same recipe is a new `run_id` linked by `recipe_hash` | 16 bytes |
| `pack_id` | `blake3(canonical bytes of pack_manifest.json)`, 128 bits, one per release | permanent | 16 bytes |
| `dense_contig_id`, `dense_slot_id`, `window_id` | release-local contiguous integers for serving indexes only | **never** exposed as a join key outside a release; regenerated at every release | uint32 (partitioned by tier/epoch if > 2^31) |
| `release_id` | `rS<stage>.<n>` (e.g. `rS1.1`, `rS2C.1`, per STAGES.md) plus the manifest hash | permanent | text + hash |
| `stats_release_id` | `release_id` of the last full L4 recompute | permanent | text |
| `annot_version` | hash of (container digest, tool versions, DB versions, parameters) for L1 | permanent | text |

Collision policy: UIDs are 128-bit BLAKE3, so the birthday bound is negligible even at 10^12 items (64-bit UIDs would reach ~2.7% collision probability at 10^9.5 contigs, which S7 could exceed). Dense release-local uint64 IDs are used for compact joins inside a release. The uniqueness test (§4.1) still checks for collisions; any collision is salted and recorded in `uid_collisions.parquet`. `raw_digest` (exact bytes) and L0 object checksums stay sha256.

`cluster_uid` has one definition, the row above. `PLAN.md` (PLAN:359) and `DEPLOYMENT-AND-MALLEABILITY.md` (DEPLOY:130) are being aligned to it (`planning/search/06_infra_inventory.md` §3 item 18).

`id_history.parquet`: `(entity_type, old_uid, new_uid, relation ∈ {same, split, merged, retired, renamed}, from_release, to_release, evidence)`. Sidecars must declare `key_type` and `key_release`. `bio-box register-sidecar` refuses a key from a different release unless it is routed through `id_history`, and it warns if key coverage is below 95%. Coverage is the fraction of the sidecar's keys that resolve in the release (§2B.6), not the fraction of the release's keys present in the sidecar.

---

## 2A. L0-seq: raw sequence store (first-class layer)

Requirement: any project can retrieve the raw nucleotide sequence of every source contig and every dedup member, not just representatives, plus the protein and IGR sequences, sliced by `(contig_uid, start, end, strand)`. Dedup and clustering decide what is *counted*, never what is *kept*.

### 2A.1 What is stored and how

| object | stored as | dedup in storage | notes |
|---|---|---|---|
| contig nucleotide sequence | **`seq_chunks`**: Parquet rows `(raw_digest, chunk_ix, seq)`. Each row holds a 65,536-byte chunk of the sequence exactly as in the source (case and IUPAC preserved), with a `contig_len` column. Sorted by `(raw_digest, chunk_ix)`, zstd, ~1,024 chunks per row group | Byte-identical sequences are stored once (`raw_digest`). Near-identical members (99% ANI dedup, same strain from two sources) are **all stored** | `contig_uid` = digest of the uppercase sequence (§2). `raw_digest` = sha256 of the exact bytes, so case variants are kept distinct and exact round-trip stays possible |
| source record occurrence | **`contig_occurrence`**: `(source_uid, source_release, l0_object_sha256, record_ordinal, accession, original_header, raw_digest, contig_uid, line_width)` | one row per source record, never collapsed | This is what reproduces any source file and answers "which sources contain this contig" |
| dedup/cluster membership | `contig_member` (`member contig_uid → representative contig_uid, method, ANI/containment, release_id`) and `genome_identity` (L0.5) | n/a | Representatives are a view; members stay first-class |
| protein sequence | **`protein_seq`**: Parquet `(protein_uid, protein_digest, length, seq, gene_call_version)`, sorted by `protein_uid`. Selenocysteine `U`, pyrrolysine `O` and frameshift-stitched sequences are stored as reconstructed | stored once per distinct sequence | `slot_protein` maps `(slot_uid, gene_call_version) → protein_uid` plus genetic code and edit type, so the same protein in 10^4 genomes costs one row of sequence |
| IGR and slot nucleotide | **coordinates only**: `igr(contig_uid, start, end, strand, left_slot, right_slot, igr_version)` and slot coordinates in L1/L3 | n/a | sliced from `seq_chunks` on demand; no second copy of the bases |
| tool-interop FASTA | bgzip FASTA + `.fai` + `.gzi` per shard, **generated on demand** by `bio-box seq extract` into `/tmp` or a project directory | derived, disposable | MMseqs2/Infernal/HMMER read streamed FASTA. Nothing downstream depends on a persistent FASTA copy |

Why Parquet chunks are primary and bgzip+faidx is secondary:
- **SQL-native.** The store is queryable and sliceable from DuckDB with no extension or UDF (§2A.3), which is the "retrievable by any project" path.
- **Same layout, partitioning and manifest machinery as every other layer**, so manifests, the ledger, determinism and accretion tests apply unchanged.
- **Content-addressed storage dedup is natural.** The same plasmid from PLSDB, RefSeq and IMG/PR is stored once and has three occurrence rows.
- **bgzip+faidx stays available.** It is the de facto interop format, widely supported (samtools, pysam, pyfaidx), and is produced by `extract` from the Parquet store.

Test K-13 compares both on the pod at S0 and is decided by LY-06 (§4.11): for the canonical slice queries, the format that reads fewer storage units (Parquet row groups vs BGZF blocks, with bytes decoded) wins; extract throughput and bytes/bp are recorded; R-08 latency is checked afterwards, not used to decide. If bgzf wins, the primary flips to bgzf shards (keeping the same `raw_digest` and occurrence tables), recorded as an ADR. zstd-seekable FASTA is rejected because samtools, pysam and pyfaidx don't support it.

Physical layout:
```
l0seq/seq_chunks/stage=<S>/digest_prefix=<00..ff>/part-<k>.parquet   # 256 MB–1 GB files; prefix buckets keep files balanced
l0seq/contig_occurrence/stage=<S>/source=<name>/part-<k>.parquet
l0seq/protein_seq/gene_call_version=<v>/stage=<S>/protein_prefix=<00..ff>/part-<k>.parquet
l0seq/igr/igr_version=<v>/stage=<S>/part-<k>.parquet
l0seq/seq_locator.parquet   # (raw_digest → file, row_group, first_chunk_row): rebuilt per stage, tiny, keeps slice reads to one row group
```
Once a stage has written its shards, it never touches an earlier stage's files. A sequence first seen in S3 that already exists from S1 adds only occurrence rows (`raw_digest` already present).

### 2A.2 Lineage: slot → contig → source

```
slot_uid (contig_uid, slot_ordinal) ──► slot coordinates (start, end, strand, gene_call_version)
   ├─► slot_protein ──► protein_seq (protein_uid)
   └─► contig_uid ──► contig_occurrence (1..n rows) ──► source_uid, source_release, accession,
                         original_header, record_ordinal ──► l0_object_sha256 ──► sources.yaml contract
                     └─► contig_member ──► representative / other members (release-scoped)
```
Coordinates are 0-based half-open on the stored contig. Because stored bytes equal source bytes, source coordinates equal stored coordinates (offset 0); no liftover is needed. `bio-box seq lineage <uid|accession>` prints this chain as JSON.

### 2A.3 Retrieval APIs

**CLI (in the `bio-box` container):**
```
bio-box seq get --contig <contig_uid|accession> [--start S --end E --strand +|-] [--fasta]
bio-box seq get --slot <contig_uid:ordinal> [--flank N] [--nt|--aa]
bio-box seq get --protein <protein_uid>
bio-box seq get --segment <segment_uid> [--aa | --nt --slot <contig_uid:ordinal>]
bio-box seq get --igr <contig_uid:left_slot>
bio-box seq slots --protein <protein_uid> [--species-reps] [--tier T] [--limit N]   # reverse lookup via protein_occurrence
bio-box seq extract --uids uids.txt | --sql "SELECT ..." --out out.fa.gz [--segment] [--bgzip --faidx] [--members|--reps-only]
bio-box seq lineage <uid|accession>
bio-box region translate <contig_uid:start-end> [--six-frame] [--min-orf 30]
```
- `extract` sorts requests by `(file, row_group)` before reading, so a batch of 10^5 slices is one sequential pass per touched file, not 10^5 random reads.
- `--sql` is the canonical extract flag. The `--query <view> --kind` form in `DEPLOYMENT-AND-MALLEABILITY.md` (DEPLOY:104) is not adopted (`06_infra_inventory.md` §3 item 22).
- **Segment mode.** `--segment` slices the domain envelope `[env_start, env_end)` from `protein_seq`. In `extract`, `--segment` means the UIDs (or the `--sql` result column `uid`) are `segment_uid`s. `--nt` needs `--slot`, because one protein has different codons in different slots.
- **Reverse lookup `protein_uid → slots`.** `seq slots` reads `protein_occurrence`, which is sorted by `protein_uid`, so one lookup touches one row group per file (zone maps, K-06). A search hit list is protein-first; this is the join that gives every hit its locus context (SEARCH_PLAN SD-02).
- **`region translate`** translates a contig interval from `seq_chunks` under the genetic code recorded for the contig; `--six-frame` emits stop-to-stop ORFs of at least `--min-orf` aa in all six frames. It is the on-demand form of the `dark_orf` pack target (§2B.2).

**DuckDB (pure SQL, in the generated catalog):**
```sql
-- forward-strand slice, 0-based half-open
CREATE MACRO seq_slice(cuid, s, e) AS (
  SELECT substr(string_agg(seq, '' ORDER BY chunk_ix),
                s - min(chunk_ix) * 65536 + 1, e - s)
  FROM seq_chunks sc JOIN contig_index ci USING (raw_digest)
  WHERE ci.contig_uid = cuid
    AND chunk_ix BETWEEN s // 65536 AND (e - 1) // 65536);
CREATE MACRO revcomp(x) AS reverse(translate(x, 'ACGTRYKMBVDHNacgtrykmbvdhn', 'TGCAYRMKVBHDNtgcayrmkvbhdn'));
CREATE MACRO seq_get(cuid, s, e, strand) AS
  CASE WHEN strand = -1 THEN revcomp(seq_slice(cuid, s, e)) ELSE seq_slice(cuid, s, e) END;
-- domain envelope of a segment, 0-based half-open on the protein
CREATE MACRO segment_seq(suid) AS (
  SELECT substr(p.seq, s.env_start + 1, s.env_end - s.env_start)
  FROM segments s JOIN protein_seq p USING (protein_uid)
  WHERE s.segment_uid = suid);
-- reverse lookup: every slot that carries a protein
CREATE MACRO slots_of(puid) AS TABLE
  SELECT * FROM protein_occurrence WHERE protein_uid = puid;
```
`contig_index` maps `contig_uid → raw_digest`. Where several raw variants share a `contig_uid`, the default variant is the earliest occurrence by precedence, and the other variants are reachable by `raw_digest`. Python clients may also call `atd.seq.get(...)`, which uses `seq_locator` to read exactly one row group. Both paths must return identical bytes (R-03).

**Export to FASTA:** `extract`, as above. The `--bgzip --faidx` output is directly usable by samtools, pysam, pyfaidx, MMseqs2 `createdb` and Infernal.

### 2A.4 Sizing per stage *(estimates; measure at each stage)*

Assumptions: zstd-compressed ASCII DNA ≈ 0.27–0.30 bytes/bp *(measure, K-13)*; proteins ≈ 0.6 bytes/aa after zstd *(measure)*; coding fraction ~88%; protein storage after exact-sequence dedup is ≤ 1/3 of nucleotide bytes.

Stage labels follow `STAGES.md` §2 (S2C from `planning/search/SEARCH_PLAN.md` §6.1). Rows S0–S4 are per stage increment. S5 and S6 are sized over their whole GTDB set; the increment is smaller by the genomes already ingested at S2C and S3, which are not re-ingested. `seq_chunks` ≈ 0.27–0.30 × bp, `protein_seq` ≤ 1/3 of that, occurrence + IGR + locator ≈ 0.02 × bp (from the earlier S3 row: 10 GB at 500 Gbp). All *estimates*.

| stage | content (STAGES.md) | input bp | `seq_chunks` | `protein_seq` | occurrence + IGR + locator |
|---|---|---|---|---|---|
| S0 | ~14 golden genomes + CAMI II subset + spike-ins | ~0.05 Gbp (+1–3 CAMI) | ~15–20 MB (+0.3–0.9 GB) | ~6 MB (+0.1–0.3 GB) | < 10 MB (+< 60 MB) |
| S1 | INPHARED complete phages (30–40k) | 1.7–2.2 Gbp | ~0.5–0.7 GB | ≤ 0.2 GB | ~0.04 GB |
| S2 | PLSDB plasmids (72,360) | 5–7 Gbp | ~1.4–2.1 GB | ≤ 0.7 GB | ~0.1 GB |
| S2C | GTDB R226 `Complete Genome` assemblies, ≤ 20 per species (est. 15–30k) | 60–120 Gbp | ~16–36 GB | ≤ 12 GB | ~1.2–2.4 GB |
| S3 | GTDB R226 species reps < 2 Mbp (9–12k), minus those already in S2C | 10–13 Gbp | ~2.7–3.9 GB | ≤ 1.3 GB | ~0.2–0.3 GB |
| S4 | IMG/VR v4 HQ vOTU reps (231,408) (+ IMG/PR 699,973), local-only | ~12 (+10–30) Gbp | ~3.2–3.6 GB (+2.7–9 GB) | ≤ 1.2 GB (+≤ 3 GB) | ~0.2 GB (+0.2–0.6 GB) |
| S5 | all GTDB R226 species reps (143,614) | 400–450 Gbp | ~110–135 GB | ≤ 36–45 GB (≤ 1/3 rule); 66–74 GB at `sources_and_sizing.md` 0.165 GB/Gbp, being reconciled (`06_infra_inventory.md` §3 item 26) | ~8–9 GB |
| S6 | all GTDB R226 genomes (732,475) | ~2.7–3.1 Tbp | ~0.8–0.9 TB (near-identical members stored in full) | ~0.1–0.2 TB (heavy exact dedup) | ~30 GB |
| S7+ | AllTheBacteria, SPIRE, Logan | Tbp–Pbp | off-box; not sized here | | |

Everything through S6 fits within the 30 TB pod scratch with room for L1–L4. On the pod, the store lives on NFS as ≥256 MB files; hot shards for a job are copied to `/tmp`.

### 2A.5 Tests (added to §4)

| ID | op (MR) | test | gate |
|---|---|---|---|
| R-01 | O-INGEST (MR-14) | Per-record digest: for every source record, sha256 of the sequence bytes parsed directly from the source file (whose object sha256 matches the contract) equals the stored `raw_digest` and the digest of the sequence reassembled from `seq_chunks` | B |
| R-02 | O-INGEST (MR-14) | Whole-object round trip: reconstruct each L0 object's normalized FASTA from `contig_occurrence` (source order, original headers, stored line width) and `seq_chunks`. Its sha256 must equal the sha256 of the decompressed, line-normalized source computed at ingest. Exact byte round trip of the compressed file isn't required, since gzip bytes aren't reproducible | B |
| R-03 | O-SLICE (MR-14) | Random-slice equality: per stage, 10,000 random `(contig, start, end, strand)` slices compared to `samtools faidx` (with `-i` for minus strand) on the decompressed source FASTA. The sample must include start = 0, end = len, length-1 slices, empty slices, slices spanning one and several chunk boundaries, minus-strand slices, IUPAC/lowercase regions, and both members and representatives. The CLI, the DuckDB macro and the Python path must all return identical bytes | B |
| R-04 | O-TRANSLATE, O-CALL | Protein consistency: for unedited genes, `protein_seq` equals the translation of the slot's nucleotide slice under its recorded genetic code. For edited genes (frameshift, Sec/Pyl, intron), it equals the recorded reconstruction and the edit coordinates resolve. Every `slot_protein` row resolves | B |
| R-05 | O-SLICE | IGR consistency: IGR coordinates lie strictly between their flanking slots, and the slice equals the source slice | B |
| R-06 | O-INGEST, O-IDENT | Lineage completeness: every slot, every dedup member and every protein resolves to at least one `contig_occurrence`, then to an L0 object, then to a contract. No orphan chunks (every `raw_digest` has at least one occurrence) | B |
| R-07 | O-EXTRACT | Extract round trip: `bio-box seq extract` output re-ingested gives identical per-record digests, and `samtools faidx` accepts the `.fai`/`.gzi` | B |
| R-08 | O-SLICE, O-EXTRACT | Retrieval latency on wilkcomp02, warm and cold, p50/p95, scaling-fit across stages. Targets: single slice ≤ 50 ms warm / ≤ 500 ms cold; full 5 Mbp contig ≤ 1 s warm; batch extract of 10^4 slices ≤ 30 s; 10^6 proteins ≤ 2 min | W; B if > 10× target |
| R-09 | O-INGEST, O-ACCRETE (MR-08) | Append-only: entering a stage leaves the content hashes of every earlier `l0seq` file unchanged. A sequence already stored adds occurrence rows only | B |
| K-13 | O-SLICE (LY-06) | Format bake-off at S0: Parquet chunks vs bgzf+faidx; decided by LY-06 storage units read per canonical slice (§2A.1, §4.11); extract throughput, bytes/bp and slice latency recorded; decides the primary format | decision |

---

## 2B. Search layer

Requirement: an ad hoc sequence, profile, CM or region search runs as a recorded job over per-release targets, with a fixed search-space size and locus context on every hit. `planning/search/SEARCH_PLAN.md` is authoritative (SD-01..SD-21); this section is its data-layer view. `serve` never runs a search inline; it can submit and poll (SD-01).

### 2B.1 Tables

Types: UIDs are 16-byte BLAKE3 (§2); coordinates are 0-based half-open.

| table | layer | key (sort order) | columns | S5 (*estimate*) | S6 (*estimate*) |
|---|---|---|---|---|---|
| `protein_occurrence` | L1-derived, per stage | `(slot_uid, gene_call_version)`; sorted by `protein_uid` | `protein_uid`, `slot_uid`, `genome_uid`, `contig_uid`, `start`, `end` (int64), `strand` (int8), `gene_call_version`, `tax_order` (uint32), `is_species_rep` (bool), `tier` (uint8), `assembly_quality` (enum), `dist_to_contig_end_bp` (uint32) | ~6e8 rows × ~48 B ≈ 29 GB | ~3e9 × ~48 B ≈ 144 GB (*derived*) |
| `cluster_member_aln` | L2 | `(member_uid, level, release_id)`; sorted by `(level, rep_uid, member_uid)` | `member_uid`, `rep_uid`, `level` (`nr90`, `nr50`, `subfam`), `qstart`, `qend`, `tstart`, `tend` (uint32), `pident` (float32), `alnlen` (uint32), compact CIGAR (binary run-length). Representatives carry a self row | ~4.8e8 × ~40 B ≈ 19 GB (nr90 level) | ~1e9 exact-unique × ~40 B ≈ 40 GB per level (*derived*) |
| `domain_promiscuity` | L4 | `(token_id, release_id)` | distinct fusion partners, distinct neighbor families, both raw and vehicle-deduplicated, `blocklist` (bool) | small | small |
| `background_rates` | L4 | `(token_id, regime, clade, tier, release_id)` | per-window occurrence probability `p0` (float64), window definition, `n_windows`; p0 for on-demand anchor statistics | small–medium | small–medium |
| `search_run` | search | `run_id` | `recipe_hash`, recipe YAML, argv, tool versions, image digest, `written_against` and `ran_against` `release_id`, `pack_id`, sha256 of every shard read, Z constants used, query hash, output paths, wall s, peak RSS, CPU-min per PID, host, exit code, `status` (`ok`, `failed`, `timeout`, `oom`, `stopped`, `output_cap`), `partial_policy` (`exclude`, `include`, `flag`), `max_output_gb`, scaled nucleotide cutoffs if any (SEARCH_PLAN SD-20) | tiny | tiny |
| `search_hit` | search | `(run_id, round, target_kind, target_uid, env_start)` | `run_id`, `round`, `target_kind` (`protein`, `segment`, `dark_orf`, `igr`, `contig_region`, `spacer`), `target_uid` (`protein_uid`, `segment_uid` or `region_uid`), `bitscore` (float32), `evalue_zref` (float64), `target_size` (real, uint64), envelope `env_start`/`env_end`, alignment coordinates on the target `ali_start`/`ali_end`, alignment coordinates on the profile `hmm_from`/`hmm_to` and `model_len` (uint32; Domainator's `rstart`/`rend`/`rlen`), `is_partial` (bool, from `protein_occurrence` edge flags), `first_inclusion_round`, `score_kind` (`exact`, `projected`), curation label and reason (text, never a number) | per run | per run |
| `search_ssn` | search | `(run_id, edge_i, edge_j)` | maximum spanning forest ∪ kNN edges over the run's representative hit set, `weight` (EFI-style score from DIAMOND or MMseqs2 bitscores), `is_mst` (bool) | per run | per run |
| `search_ssn_thresholds` | search | `(run_id, threshold)` | one row per distinct edge weight: `cluster_count`, `largest_component`, merge events; any cut is a lookup (SEARCH_PLAN SD-19) | per run | per run |
| `search_profile` | search | `(run_id, round)` | `kind` (`hmm`, `cm`, `pssm`), profile blob (zstd), MSA (Stockholm, zstd), `n_seq`, `n_eff` | per run | per run |
| `user_token_overlay` | L3 | `(slot_uid, user_token_id, run_id)` | `slot_uid`, `user_token_id`, segment envelope `env_start`/`env_end`, `run_id` | sparse | sparse |

`search_hit` joins to `protein_occurrence` on `target_uid = protein_uid` (protein kinds) or to `seq_chunks` coordinates through `region_uid` (nucleotide kinds), so every hit carries contig, coordinates, strand, genome, taxonomy, species-rep flag and tier (SD-02).

### 2B.2 Search pack (per release; SD-06)

```
$ATD_ROOT/search_pack/release=<R>/
  pack_manifest.json        # per shard: target, path, sha256, bytes, n_seq, n_residues,
                            # source tables + their manifest hashes, builder tool + version; pack_id = blake3(this file)
  nr50/      part-<k>.faa.zst + MMseqs2 DB (--prefilter-mode 1 index)   # default target (SD-03)
  nr90/      part-<k>.faa.zst                                           # final round, projection
  unique/    part-<k>.dmnd                                              # all exact-unique proteins; DIAMOND
  segments/  part-<k>.faa.zst                                           # domain envelopes + uncovered segments ≥ 50 aa
  igr/       part-<k>.fna.zst + igr_meta.parquet                        # IGRs ≥ 30 nt; flank genes, orientation, GC, distance to contig end
  igr90/     part-<k>.fna.zst                                           # IGR90 representatives
  dark_orf/  part-<k>.faa.zst + dark_orf_meta.parquet                   # stop-to-stop ORFs ≥ 30 aa outside called genes, all frames in 3 kb at each contig end
  spacers/   part-<k>.fna.zst + index                                   # MinCED ∪ CRISPRCasFinder spacers with array, repeat, contig
  models/    pressed library HMMs, Rfam CMs, user models
  hh/        Pfam-A HHM, PDB70                                          # from S2C (ST2)
  3di_nr50/  3Di of NR50 representatives (ProstT5 or ESM-2 3B, by SEARCH_PLAN SD-21 bake-off)   # from ST3
```
- Every shard is ≥ 256 MB except the last of each target, so a target smaller than 256 MB (all of S0) is one file.
- Built by `bio-box pack build <release>` after `barrier_L2_<release>`; verified by `bio-box pack verify <release>` (C-12) before the release manifest is signed.
- Targets per stage: S0 has `unique`, `igr`, `dark_orf`, `spacers` and `models` (no reduction, ST0); from S1, `nr50`, `nr90`, `segments` and `igr90`; `hh` from S2C; `3di_nr50` from ST3. At S6, `unique` lives on the Hetzner search host only.
- A run copies the shards it needs to `/tmp/atd/search/<run_id>/`, checks each sha256 against the manifest, and records them in `search_run`. During L2 barriers the lane reads shards from NFS instead (§5).

### 2B.3 Run outputs

```
$ATD_ROOT/search/run=<run_id>/
  run.parquet               # the one search_run row
  hits.parquet              # search_hit, all rounds
  profiles.parquet          # search_profile, all rounds
  neighborhoods.parquet     # bio-box nbr output for the final hit set
  ssn.parquet               # optional: search_ssn + search_ssn_thresholds (bio-box search ssn)
  ssn.dsnv                  # optional: Domainator viewer bundle, generated from ssn.parquet
```
- Exports are derived and disposable, written to a project directory or `/tmp`, never into the run directory: `--format domainator-gb` (GenBank with `Domainator` features; `evalue` at Z_ref, `score` in bits), `tsv`, and `--json` NDJSON (SEARCH_PLAN SD-17, SD-18).
- **Output cap.** The runner stops a run whose outputs pass `max_output_gb` (default 25 GB) and records `status = output_cap`; partial outputs are kept and flagged.
- Written in `/tmp/atd/<job_uid>/`, packed, then one `rsync` and an atomic rename (§5 local staging). A failed, OOM or timed-out run still writes `run.parquet` with its status; it never writes an empty `hits.parquet` as if it had zero hits.
- **Small-file waiver.** Run files are often < 64 MB. The waiver is bounded: at most 6 files per run (run, hits, profiles, neighborhoods, optional `ssn.parquet` and `ssn.dsnv`), and runs are ~10^3 per year (*estimate*), so `search/` stays < 10^4 files. The same waiver covers `l3/user_overlay/user_token=<id>.parquet`, one file per promoted family.
- The catalog exposes `search_run`, `search_hit` and `search_profile` as views over the run directories listed in `search/run_index.json`, which the runner rewrites by atomic rename when a run closes.

### 2B.4 Z policy (SD-04)

| constant | value | applies to |
|---|---|---|
| `Z_prot` | 1e8 sequences | HMMER `hmmsearch`/`jackhmmer`/`phmmer` `-Z` (pyhmmer `Z=`) |
| `domZ` | 1e8 | HMMER `--domZ` |
| `Z_nt` | 1e5 Mb | Infernal `cmsearch -Z`, `nhmmer -Z` |

- The constants are fixed for the whole release line (every stage, every release). They are recorded in the release manifest and in every `search_run`. Changing them starts a new release line and a new `annot_version`.
- Every hit, library or ad hoc, stores the bitscore, `evalue_zref` (E at the fixed Z) and the real target size. E at the real size is `evalue_zref × target_size / Z_ref`. For a tool with no sequence-count Z option (MMseqs2, DIAMOND, BLAST+, whose size options are in residues), `evalue_zref = E_native × Z_ref / target_size`, with `target_size` in sequences.
- Library hits keep their GA/TC bitscore thresholds, which do not depend on Z. Recipe thresholds are stored as bits. E-thresholds from papers are converted to bits at the paper's database size when it is stated; otherwise at Z_ref, and flagged.
- `target_size` is stored in `search_hit` and `search_run` only, never in L1, because it changes with the stage.

**A-02 vs Z, resolved** (`06_infra_inventory.md` §3.C item 16). A-02 requires byte-identical L1 row-sets between the incremental and from-scratch builds. A per-stage or per-release Z would change stored E-values at every stage, so A-02 would fail on any stage boundary. With a constant Z, a library hit's bitscore and `evalue_zref` depend only on the protein sequence, the model and the pinned tool version, so L1 E-values are stage-invariant and A-02 byte identity holds with E-values stored in L1. S-03 checks the same invariance for search runs.

### 2B.5 Base-scan unit and projection (SD-05)

| stages | unit scanned by each library | hits on other proteins |
|---|---|---|
| S0, S1, S2, S2C | every exact-unique protein (distinct `protein_uid`), once | none; all hits are exact. Per-genome L1 rows are the join through `slot_protein` |
| S3 on | NR90 representatives (new representatives only after `clusterupdate`) | projected to members through `cluster_member_aln` (level `nr90`) |

This replaces the three earlier statements: per genome (L1 row in §1), NR50 residues (PLAN:118) and NR90 representatives (`sources_and_sizing.md`:85-86) (`06_infra_inventory.md` §3 item 15).

**Projection rule.** The member's envelope is the representative's envelope mapped through the member's alignment, and its projected bitscore is the representative's (`score_kind = projected`). A member is re-scored exactly (`score_kind = exact`) when either:
1. its projected bitscore is within 10 bits of the cutoff (GA/TC for a library model, `inclusion_bits` for a recipe), or
2. its alignment covers < 90% of the representative's domain envelope.

A member with no `cluster_member_aln` row is a C-11 failure, never a silent skip. The same rule applies to the final round of an ad hoc search (NR90 representatives → members) from S1. S-08 checks projection accuracy.

Projected rows depend on which protein is the representative, which is an L2 choice. They are therefore keyed by `release_id` (like L3 tokens), not by `(genome_uid, annot_version)`, and A-02 byte identity applies to exact rows only. Projected rows are compared under A-03: calls (hit or no hit, family, envelope within 5 residues *(tunable)*) must agree after the crosswalk.

### 2B.6 Sidecar coverage

Coverage of a sidecar = (sidecar keys that resolve in the release) / (sidecar keys). `register-sidecar` warns below 95% and refuses below 50% (DEPLOY:149). The fraction of the release's keys present in the sidecar is reported, never gated, so sparse hit lists from `scan-new-models` or a search run are accepted (`06_infra_inventory.md` §3 item 5).

### 2B.7 Sizes *(estimates; measure from S0)*

| item | S5 | S6 |
|---|---|---|
| `nr50` (FASTA zstd + MMseqs2 index) | 33–56 Gaa; ~20–35 GB + ~35–60 GB | 40–80 Gaa |
| `nr90` | 80–110 Gaa | 115–200 Gaa |
| `unique.dmnd` | ~100–150 GB | ~330 GB; Hetzner only |
| `segments` | ~1.1× NR90 residues | not sized |
| `igr` + `igr90` | 48–54 Gbp; ~14 GB zstd | 324–372 Gbp |
| `dark_orf` | ~30 Gaa (measure at S0) | ~200 Gaa |
| `spacers`, `models` | small | small |
| `hh` | ~30–60 GB | same |
| `3di_nr50` | ~50 GB without Cα | ~60–120 GB |
| **pack total on NFS** | **~370–530 GB** (row sum, `planning/sources_and_sizing.md`) | **~0.8–1.1 TB** (*derived*: S6 column at ~0.6 B/aa and ~0.28 B/bp zstd, plus `unique.dmnd`, `hh`, `3di_nr50` and a 2× aa MMseqs2 index). Fits on NFS, but the S6 working set exceeds the 200 GB scope and `/tmp`, so the S6 search host is on Hetzner (≥ 512 GB RAM, ≥ 4 TB NVMe, *estimate*; SEARCH_PLAN ST4) |
| per-run working set in `/tmp` | 60–120 GB (NR50 + IGRs), inside the 150 GB lane quota (§5) | same on the Hetzner search host |
| `protein_occurrence` | ~29 GB | ~144 GB |
| `cluster_member_aln` | ~19 GB (nr90 level) | ~40 GB per level |

---

## 3. Physical layout

### 3.1 Directory tree (identical on the pod and on Hetzner; only the root differs)

```
$ATD_ROOT/                                   # /stor/scratch/Wilke/if2453/all-the-data on the pod; /mnt/nvme/all-the-data on Hetzner (DEPLOYMENT §3.9)
  l0/source=<name>/release=<rel>/part-<k>.tar.zst      # ≥ 256 MB packs; never one file per genome
  l0/contracts/sources.yaml, l0_manifest.parquet
  l0seq/{seq_chunks,contig_occurrence,protein_seq,igr}/...   # raw sequence store, layout in §2A.1
  l1/annot_version=<v>/kind=<genes|domains|ncrna|mge|qc|tax>/stage=<S>/shard=<nnnnn>.parquet
  l1/domains_projected/release=<R>/shard=<nnnnn>.parquet  # projected hits from S3 on (§2B.5)
  l1/protein_occurrence/gene_call_version=<v>/stage=<S>/part-<k>.parquet
  l2/release=<R>/{nr90,nr50,subfam,igr}/...parquet; cluster_member/; cluster_member_aln/; id_history.parquet
  l3/tapes/stage=<S>/shard=<nnnnn>.parquet               # coordinates (stable)
  l3/tokens/release=<R>/shard=<nnnnn>.parquet            # token projection (versioned)
  l3/slot_pairs/release=<R>/anchor_bucket=<bb>/part-<k>.parquet
  l3/user_overlay/user_tokens.parquet; user_token=<id>.parquet   # one file per promoted family; same waiver as search/
  l4/release=<R>/{pair_stats,closure,pangenome,transitions,domain_promiscuity,background_rates}/...
  l5/release=<R>/...                                      # only if built
  search_pack/release=<R>/{nr50,nr90,unique,segments,igr,igr90,dark_orf,spacers,models,hh,3di_nr50}/...   # ≥ 256 MB shards (§2B.2)
  search_pack/release=<R>/pack_manifest.json
  search/run=<run_id>/{run,hits,profiles,neighborhoods}.parquet   # packed per run (§2B.3)
  search/run_index.json
  ledger/stage=<S>/step=<name>.parquet
  tests/test_run/<yyyy-mm-dd>-<k>.parquet   # packed test results (§4.12); small-file waiver there
  tests/test_run_index.json
  manifests/release=<R>.json   ->  current.json (symlink swapped by atomic rename)
  catalog/release=<R>.duckdb   (read-only, generated)
```

### 3.2 File and row-group sizing

| table | sort key inside file | target file size | row-group size | rationale |
|---|---|---|---|---|
| L1 shards | `(genome_uid, contig_uid, start_bp)` | 256 MB – 1 GB | ~122,880 rows (DuckDB default) | about 500–2,000 genomes per shard keeps the file count under 10^4 per stage and suits NFS |
| L3 tapes | `(tax_order_bucket, contig_uid)` | 512 MB – 1 GB | 8,192–16,384 contigs | smaller row groups bound the bytes decoded per neighborhood fetch; tuned by LY-06 row-group counts (§4.11) |
| L3 pairs | `(anchor_token, partner_token, tax_order_bucket)`; partitioned by `anchor_bucket` | ≤ 1 GB | 122,880 | zone maps on the leading sort key prune anchor scans; hub anchors get their own bucket |
| L4 tables | primary key | ≤ 1 GB | 122,880 | small; mostly fully scanned |
| `protein_occurrence` | `(protein_uid, genome_uid, slot_uid)` | 256 MB – 1 GB | 122,880 | protein-first reverse lookup prunes on zone maps |
| `cluster_member_aln` | `(level, rep_uid, member_uid)` | 256 MB – 1 GB | 122,880 | projection reads all members of a representative in one range |
| `search_hit` | `(round, target_kind, target_uid)` | one file per run | 122,880 | per-run; read whole |

Sort keys and row-group sizes above are starting values. LY-06 row-group counts on the micro-corpus and S0 data decide them at S0 (§4.11); LY-02..LY-04 then check them on every stage.

All Parquet is written by DuckDB `COPY ... (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE n)` with explicit `ORDER BY`, so files are reproducible given the same DuckDB version and thread count. Determinism is tested in §4.5. Partitions use hive layout so DuckDB prunes on partition keys ("filters on the partition keys are automatically pushed down into the files"). Partition keys (`stage`, `gene_call_version`, `igr_version`, `digest_prefix`, ...) live only in the path, never as file columns, so readers never see two copies that could disagree; a primary key may include a partition key. The schema registry (`src/biobox/registry/schemas.py`) enforces this in ring R0. Row-group pruning relies on Parquet min/max statistics (zonemaps), which DuckDB uses to skip row groups. Page-index skipping is not assumed.

### 3.3 Catalog: generated read-only DuckDB file (chosen), DuckLake as a gated alternative

**Choice:** each release has a JSON manifest that lists every file with its sha256, row count, schema hash and layer. A generator script turns the manifest into `catalog/release=<R>.duckdb`, which contains only views and macros over explicit file lists (no globs, so a stray file can't leak in). The file is copied to local `/tmp` and opened with `access_mode='READ_ONLY'`. `current.json` is a symlink swapped by `rename(2)`, which is atomic on a single filesystem.

Why:
- DuckDB allows one read-write process or many read-only processes on a database file. The docs warn to "exercise extra caution when accessing a DuckDB database file in a shared directory", naming network-attached storage. A read-only catalog copied to local disk sidesteps both NFS locking and the single-writer limit. `serve` and N analysis clients can then all read.
- Building a release is a batch job with one writer (the orchestrator), so the build needs no multi-writer catalog.
- Rebuilding the catalog from the manifest makes it disposable. The manifest is the system of record and is small and diffable in git.

**DuckLake (evaluated, not adopted for S0–S2):** DuckLake v1.0 was released in April 2026, with backward-compatibility guarantees. It provides snapshots, time travel and schema evolution without rewriting Parquet, plus a "frozen DuckLake" read-only mode. With DuckDB as the catalog it allows a single client only, and multi-writer use needs PostgreSQL or another ACID catalog. It also has data inlining for small appends, which addresses the small-file problem. It is attractive for L4 time travel across releases. Adopt it only if a spike at S2 passes all of these:
1. Frozen-catalog reads from 8 concurrent clients on the pod.
2. Snapshot time travel reproduces release R−1 query results byte-identically.
3. No degradation over the manifest approach in the §4.7 tests.
4. The catalog (SQLite on local disk during build) survives a kill-and-resume.

Until then, release-to-release comparison uses two catalogs side by side.

---

## 4. Test suite

Tests run under `pytest`, with DuckDB SQL assertions and `pandera` schemas. dbt-duckdb is not used because its model/test DAG duplicates the orchestrator. Each test has an ID and a gate level: **B** = blocks the stage, **W** = warns and is recorded in the release notes.

**Organization.** [`planning/tests/TEST_ARCHITECTURE.md`](tests/TEST_ARCHITECTURE.md) is authoritative for how these tests are structured: an operation registry (O-*) that every workflow is composed from, failures ranked by silent wrong answers (F-*), metamorphic relations (MR-01..MR-19), workflow manifests that generate tests for new workflows, a micro-corpus and synthetic worlds with planted truth and null worlds, layout tests LY-01..LY-10, test rings R0–R5, and the `test_run` results table. The IDs below keep their meaning; that file maps each family into the structure (§9) and names the gaps.

**Conventions in the tables below.** The "op (MR)" column names the operation each test covers (TEST_ARCHITECTURE §2, §9) and, in parentheses, the metamorphic relation (§4) or layout test (§7) the test instantiates, where it does. "all writers" means every operation that writes Parquet. Property tests use Hypothesis; mutation testing uses `mutmut` (§4.12). Every test result is a `test_run` row (§4.12), and each gate in §4.9 and §4.10 is a query over that table.

### 4.1 Schema and data contracts (every layer, every stage)

| ID | op (MR) | check | gate |
|---|---|---|---|
| C-01 | all writers (LY-01) | Each Parquet file's schema hash equals the registered schema for `(table, schema_version)`; types match exactly (no silent int64→float) | B |
| C-02 | O-IDENT, O-CATALOG | Non-null, unique primary keys; full-digest uniqueness for truncated UIDs (collision check) | B |
| C-03 | O-CATALOG, O-PROMOTE | Every foreign key resolves (slot→contig→genome→source; segment→protein; member→cluster; `search_hit`→`search_run`; `user_token_overlay`→slot, run) | B |
| C-04 | O-TAPE | Parallel `LIST` columns in a tape row have equal lengths | B |
| C-05 | O-CALL, O-TAPE, O-SCAN | Coordinate bounds: `0 ≤ start < end ≤ contig_len`; `env_end ≤ protein_len`; strands in {+1,−1} | B |
| C-06 | all writers (LY-01) | Enum columns hold only registered values (feature kinds mapped to Sequence Ontology terms) | B |
| C-07 | O-PAIRS | Saturating fields (e.g., `delta_bp` if narrowed) carry an overflow flag whenever clipped; share of clipped values < 0.1% | W |
| C-08 | all writers, O-RUN | Every row carries `annot_version`/`release_id` provenance, and the tool versions in it equal the container lock | B |
| C-09 | O-CALL, O-SCAN | Distribution checks against the previous stage, per tier and per order: gene density, mean protein length, GC, fraction of slots with an HMM hit | W, or B if > 5σ |
| C-10 | O-CALL, O-NBR | `protein_occurrence` completeness: every `slot_protein` row of the stage's `gene_call_version` appears exactly once (by `slot_uid`) with the same `protein_uid`, and every `protein_occurrence` row resolves to a `slot_protein` row; coordinates, strand and `dist_to_contig_end_bp` equal the slot's | B |
| C-11 | O-CLUSTER, O-PROJECT | `cluster_member_aln` coverage: every member at each level (`nr90`, `nr50`, `subfam`) of `cluster_member` has exactly one alignment row to its representative, representatives included as self rows; alignment coordinates within both sequence lengths; CIGAR length equals `alnlen` | B |
| C-12 | O-SEARCH (pack), O-CATALOG | Search pack manifest verification (`bio-box pack verify`): every shard's sha256 and byte count equal `pack_manifest.json`; `pack_id` equals the manifest hash; per-target sequence counts equal their source (`nr50` = NR50 clusters, `nr90` = NR90 clusters, `unique` = distinct `protein_uid`, `igr` = IGRs ≥ 30 nt, `spacers` = spacer calls); every shard ≥ 256 MB except the last of each target; no file on disk outside the manifest | B |

### 4.2 Ledger reconciliation

Every step writes `ledger/stage=S/step=X.parquet` with rows `(object_uid, n_in, n_kept, n_merged, n_rejected, reason_code)`.

| ID | op (MR) | check | gate |
|---|---|---|---|
| L-01 | O-INGEST, O-RUN | For each object, `n_in = n_kept + n_merged + Σ n_rejected` | B |
| L-02 | O-INGEST | Every L0 object in the stage's input list appears in the first ledger, with none missing and none extra (checked against the L0 manifest) | B |
| L-03 | O-RUN, O-ACCRETE | Downstream `n_in` equals upstream `n_kept` along the DAG | B |
| L-04 | O-INGEST, O-IDENT | The breakdown by reason code is within the expected band (e.g., FCS-GX contamination ≤ X% per source) | W |

### 4.3 Brute-force oracle equivalence (every query class, every stage at S0–S2; a sampled 1,000 queries per class from S3 on)

For each query class, a naive reference implementation in Python or plain DuckDB, with no indexes, runs over the raw L1/L3 rows. Its results must equal the production query path exactly: set equality, plus equal counts and aggregates. Every query class is an O-SERVE oracle (TEST_ARCHITECTURE §2); the column names the operation underneath.

Oracles run on three data sets, not on golden data alone: the micro-corpus (R2, every merge; TEST_ARCHITECTURE §6.1), synthetic worlds with planted truth and null worlds (R2 at micro scale, R5 at stage scale; §6.2), and the golden set and stage data (R3, R4).

Q1–Q8 follow `PLAN.md` §1.2. An earlier draft of this file numbered them differently (`06_infra_inventory.md` §3 item 20); its gapped-pattern and closure-lookup classes have no PLAN number and are kept as QG and QC. The lab questions numbered Q1–Q10 in `planning/preflight_gaps.md` are a separate list.

| query class | op (MR) | oracle |
|---|---|---|
| Q1 count (loci with A and B within ±k ORFs, filtered by tier or phylum) | O-PAIRS | nested loop over tapes, count |
| Q2 neighbourhoods (fetch k contig tapes) | O-TAPE | read tape shards fully, filter in Python |
| Q3 anchor profile (all partners within ±k ORFs of a token) | O-STATS, O-PAIRS | nested loop over tapes |
| Q4 pair drill-down (all instances of `(A,B)`) | O-PAIRS | as Q3 |
| Q5 NOGD hole over pangenome planes | O-PROF | per-species set logic from L1 directly |
| Q6 co-transmission over tree transitions | O-PROF | Fitch on the tree, recomputed per query |
| Q7 sequence slice | O-SLICE | `samtools faidx` on the decompressed source (R-03) |
| Q8 ad-hoc SQL (sample set) | O-SERVE | run once against an unindexed copy |
| QG gapped pattern `[A, B, C]` within window w, with strand/topology constraints | O-PAIRS, O-TAPE | sliding window over every tape; must include pairs that straddle any window boundary |
| QC closure lookup | O-STATS | envelope-overlap join done row by row |

**Search job classes are asynchronous.** They are submitted (`bio-box search run`, or `serve` submit and poll), run on the interactive lane (§5), and are never under `serve` timeouts or the §4.8 SLOs. Their correctness and cost are gated by family S (§4.10), not by oracles.

| job class | command | typical cost (SEARCH_PLAN §5, *derived*) |
|---|---|---|
| J1 recipe search (Q-FAM new model, Q-ITER, Q-RNA, Q-NT, Q-SSN) | `bio-box search run <recipe.yaml>` | minutes to ~1 h on the lane at S5 |
| J2 positive-set neighborhood score over the whole corpus | `bio-box nbr score` | minutes to 1 h at S5 |
| J3 on-demand anchor statistics for a run or user token | `bio-box search stats` | *(measure)* |
| J4 search pack build | `bio-box pack build <release>` | batch; part of the release DAG |

`bio-box nbr` on a stored hit set is a Q2 fetch and stays synchronous.

Tokens are sampled stratified by frequency decile, so hub tokens are always in the sample. Gate: **B** on any mismatch.

### 4.4 Accretion-equivalence

Applies at every stage boundary from S1 on. Path A builds `S_{n-1}` and then applies `S_n` incrementally. Path B builds `S_{n-1} ∪ S_n` from scratch. Both use the same inputs and the same container. Every A test is O-ACCRETE (TEST_ARCHITECTURE §2) checked at one layer, and is MR-08 (add unrelated genomes) with a from-scratch oracle.

| ID | op (MR) | layer | required equality | gate |
|---|---|---|---|---|
| A-01 | O-IDENT (MR-08) | L0 / L0.5 | Identical manifests; identical identity clusters and representatives (precedence is deterministic) | B |
| A-02 | O-CALL, O-SCAN (MR-08) | L1 | Byte-identical content hash per `(genome_uid, annot_version)` row-set, ignoring shard boundaries (hash over sorted rows, not files). Holds with stored E-values because Z is a release-line constant (§2B.4); covers exact domain-hit rows and `protein_occurrence` | B |
| A-03 | O-CLUSTER, O-PROJECT (MR-08) | L2 | Equal up to `id_history`: for each Path B cluster, the crosswalked Path A cluster has adjusted Rand index ≥ 0.99 *(tunable)* over shared members; no member unassigned. Projected domain hits (§2B.5) agree on calls after the crosswalk | B |
| A-04 | O-TAPE (MR-08) | L3 coordinates | Byte-identical sorted-row hashes | B |
| A-05 | O-TAPE (MR-08) | L3 tokens | Identical after mapping through the crosswalk | B |
| A-06 | O-STATS (MR-08) | L4 mergeable counts | Equal exactly | B |
| A-07 | O-STATS (MR-08) | L4 non-mergeable (LLR, q, entropies) | Path A is marked `stale` until recompute; after recompute, equal to Path B within float tolerance 1e-6 relative | B |

At S0–S2, Path B is cheap and runs every time. From S3 on, it runs on a 1–5% stratified sample of genomes, re-sampled each stage.

**Stateful accretion machine (complements A-01..A-07).** A-01..A-07 test one operation order per stage boundary. A Hypothesis `RuleBasedStateMachine` (TEST_ARCHITECTURE §0 item 3, §8 ring R1) draws random sequences of five rules: append (a genome batch or stage), recluster (`clusterupdate` or full), rescore (L4 recompute), compact (repack shards), remap (`id_history` crosswalk). It runs on toy data in R1 (every commit) and on the micro-corpus in R2 (every merge). Invariants after every step:
- C-02, C-03 and L-01 hold.
- Current state equals a from-scratch build of the same inputs under the A-02..A-06 equalities.
- Non-mergeable L4 values are `stale` after any append or recluster until the next rescore (A-07, F-06).
- Compact and remap change no sorted-row hash outside the crosswalk.

Failing sequences are shrunk by Hypothesis and kept as fixed regression cases.

### 4.5 Determinism

| ID | op (MR) | check | gate |
|---|---|---|---|
| D-01 | all writers (LY-08) | Same inputs + same container digest + same thread count → identical sorted-row hashes for every table | B |
| D-02 | all writers (MR-16) | Same inputs at a different thread count (e.g., 8 vs 40) → identical sorted-row hashes. File bytes may differ, so content hashes are compared, not file hashes | B |
| D-03 | O-RUN, O-CATALOG (F-11) | Kill-and-resume: kill the orchestrator at a random step and resume; manifests identical to an uninterrupted run | B |
| D-04 | O-IDENT, O-CLUSTER (MR-01) | Input-order permutation: shuffle the genome order; representatives and IDs unchanged | B |
| D-05 | all writers | Cross-host (pod vs Hetzner, S5): identical sorted-row hashes for L1–L4 on the regression sample | B at S5 |

D-02 and D-04 are MR-16 and MR-01 run on stage data; R2 runs the same MRs on the micro-corpus for every operation that declares them. Shard-count invariance for search (MR-16 on O-SEARCH) is not yet covered by a D test (TEST_ARCHITECTURE §9).

### 4.6 Golden-genome regression (ring R3: nightly, every container rebuild, every stage)

This is the fixed set from `critique_8` §golden: MG1655, *B. subtilis* 168, T4, lambda, *M. genitalium*, an SR1 genome, a Lak phage, *M. acetivorans*, *H. volcanii*, *S. acidocaldarius*, H37Rv, *B. thetaiotaomicron*, *S. pyogenes* SF370, an Obelisk, a CAMI II subset and a synthetic spike-in panel. Accessions are pinned in `golden.yaml`.

| ID | op (MR) | check | gate |
|---|---|---|---|
| G-01 | O-CALL | Gene-call agreement with curated annotation (EcoCyc/SubtiWiki/PGAP): start/stop recall and precision reported with a 95% CI; must not drop more than 1 point from the previous release | B |
| G-02 | O-CALL, O-TRANSLATE | Known special cases recovered: lambda gpG-T frameshift, T4 introns (td, nrdB, nrdD), *M. genitalium* code 4, SR1 code 25, `fdhF` Sec, *M. acetivorans* Pyl | B |
| G-03 | O-TAPE | RegulonDB transcription units vs `directon`/`regulon_locus` boundaries: agreement reported per topology class | W |
| G-04 | O-SCAN | Known riboswitch subtypes in the spike-in panel assigned correctly | B |
| G-05 | O-STATS | Drift: Jaccard of the top-1,000 pairs vs the previous release ≥ 0.9 on the golden set, else explain | W |

Golden snapshots detect change, not error: comparisons to the previous release (G-01's 1-point rule, G-05) say the output moved, not that either release is right. Each G test is therefore paired with planted truth in synthetic worlds (TEST_ARCHITECTURE §6.2) for the same operation:

| G test | planted truth |
|---|---|
| G-01 | planted genes; `truth_callerr` (frameshifts, wrong starts, missed small ORFs) |
| G-02 | `truth_callerr` frameshifts; MR-13 recoding to code 4 and 25 |
| G-03 | O-TAPE brute-force rebuild (no planted transcription units) |
| G-04 | `truth_rna` (`cmemit` motifs in IGRs) |
| G-05 | `truth_module` recovered in planted worlds; null worlds give no modules at nominal FDR (MR-18) |

### 4.7 DuckDB behavior tests (run at S0, then at every DuckDB version bump and on each new host)

Part of ring R4 when triggered. K-05 and K-06 are the partition and row-group cases of LY-06, which extends them to every canonical query (§4.11).

| ID | op (MR) | test | pass condition |
|---|---|---|---|
| K-01 | O-CATALOG | Version pin: `SELECT version()` equals the lock; Parquet written by version N is read identically by N and N+1 (forward check before any bump) | exact |
| K-02 | O-RUN | `SET memory_limit='150GB'` inside a 200 GB scope; a deliberately oversized `GROUP BY` spills to `temp_directory` instead of being OOM-killed | completes; peak RSS < 200 GB (measured with `systemd-cgtop` / `memory.peak`) |
| K-03 | O-RUN | `SET temp_directory='/tmp/atd_duck'` (local), never NFS; the test asserts no spill files appear under `/stor` | zero spill bytes on NFS |
| K-04 | O-RUN | `SET threads=N` honored: CPU usage tracks N for 8/20/40 | ±10% |
| K-05 | O-CATALOG (LY-06) | Hive partition pruning: a query filtered on `stage=`/`anchor_bucket=` reads only matching files (`EXPLAIN ANALYZE` file count; `strace -e openat` count) | files read = expected |
| K-06 | O-CATALOG (LY-06) | Row-group pruning: on the sorted pairs table, an anchor filter reads ≤ the expected row groups (from `parquet_metadata()` min/max) | measured ≤ expected + 1 |
| K-07 | O-TAPE, O-SERVE | `LIST` column cost: bytes decoded and row groups read per Q2 fetch (LY-06) vs list length and number of list columns projected, with Q2 latency recorded alongside; LY-06 decides whether tapes stay `LIST<T>` or become exploded slot tables, and the Q2 SLO is checked on the chosen layout | recorded curves; decision at S2 |
| K-08 | O-CATALOG, O-SERVE | Read-only concurrency: 1, 4, 8 and 16 clients on the same read-only catalog copy on local `/tmp`; no lock errors; p95 degradation recorded | no errors |
| K-09 | O-CATALOG | Writer exclusivity: a second read-write open of a build catalog fails fast with a clear error, never corrupting it | fails fast |
| K-10 | O-SERVE | NFS vs local: the same Q1–Q4 on Parquet stored on NFS vs copied to `/tmp`; cold (after `echo 3 > drop_caches` is not allowed, so use a fresh never-read file copy) and warm | recorded; decides what is staged locally |
| K-11 | all writers (MR-16, LY-08) | Sorted write determinism: `COPY ... ORDER BY` under different thread counts gives identical sorted-row hashes (feeds D-02) | equal |
| K-12 | O-SERVE | Large `IN`-list and semi-join behavior for k = 10/1,000/100,000 contig UIDs (decides whether Q2 needs an offset index, i.e. L5) | recorded |

### 4.8 Performance SLOs and scaling fit

Measured on wilkcomp02 under the day profile (40 threads, 150 GB DuckDB limit), with frequency-stratified tokens.

Classes and targets follow `PLAN.md` §1.2, which is authoritative where this file's earlier draft differed (Q6 was 1 s / 5 s here). Rows with no PLAN target keep this file's.

| class | p95 warm target | p95 cold target |
|---|---|---|
| Q1 count | ≤ 2 s | ≤ 10 s |
| Q2 neighbourhoods, 1,000 contigs | ≤ 1 s | ≤ 5 s |
| Q3 anchor profile from `pair_stats` | ≤ 0.5 s | ≤ 3 s |
| Q3 raw anchor scan over postings, median-frequency token | ≤ 2 s | ≤ 10 s |
| Q3 raw anchor scan, hub token (top 0.1%) | ≤ 20 s or served from L4 summary | ≤ 60 s |
| Q4 pair drill-down | ≤ 2 s | ≤ 10 s |
| Q5 NOGD hole | ≤ 1 s | ≤ 5 s |
| Q6 co-transmission screen | ≤ 5 s | ≤ 20 s |
| Q7 sequence slice (≤ 100 kb) | ≤ 0.2 s | ≤ 1 s (R-08 also applies) |
| Q8 ad-hoc SQL | no SLO; memory-capped and cancellable | none |
| QG gapped pattern, 3 tokens | ≤ 5 s | ≤ 30 s |
| J1–J4 search jobs | asynchronous; no SLO; budgets in SEARCH_PLAN §5, gated by S-06 | none |

- "Cold" means the files have not been read since being copied to the host. "Warm" means the second consecutive run.
- Each stage records p50 and p95 over at least 200 queries per class, plus bytes read and CPU seconds.
- **Scaling fit:** for each class, fit `latency = a · size^b` across stages S1…S_n, using total rows or bytes in the touched tables. Before entering a stage, extrapolate to the next stage. If the extrapolated p95 breaks the SLO, the stage gate requires either an L5 index design for that class or a scope cut. L5 is built only on this evidence.
- These are pod targets for an interactive research tool. Hetzner targets get set after the S5 port, from the same harness.

### 4.9 Stage gates (summary)

Stage numbers follow [`STAGES.md`](../STAGES.md) (S0 golden, S1 phages, S2 plasmids, S2C complete genomes, S3 small genomes, S4 IMG/VR, S5 GTDB species reps, S6 all GTDB genomes, S7+ off-box). The table was relabeled from the earlier draft numbering on 2026-10-29 using the `STAGES.md` §7 mapping: the pilot decisions (DuckLake, K-07, K-12) sit at the S2 gate, the L5 decision at the S5 gate, the drift test at the S4 gate, and the Hetzner port checks at portability steps P1/P2.

Each row is evaluated as the §4.12 gate query: the stage passes when every blocking test listed for it has `test_run` rows on the stage's release and image, all with outcome `pass`. Ring requirements (§4.12): S0 needs R0–R3 green; every later stage needs R0–R3 and R4 on its own data; R5 runs at the tier-certifying stages listed in §4.10.

| stage | must pass (B) | recorded (W) |
|---|---|---|
| S0 | C-01..C-08, C-10..C-12, L-01..L-03, R-01..R-07, R-09, Q1–Q8 oracle (full), D-01..D-04, G-01, G-02, G-04, K-01..K-13; R0–R3 green; LY-01..LY-10 on the micro-corpus and S0 data (LY-06 decides K-13, row-group sizes and sort keys); micro-corpus built with its manifest hashes recorded; mutation score measured; null-world calibration first pass run | G-03, C-09 baseline, SLO and R-08 baselines; mutation score vs the ≥ 80% target (TEST_ARCHITECTURE §11 item 3); null-world p-value uniformity and FDR |
| S1 | all S0 + A-01..A-07 (full Path B); R4; LY-01..LY-10 on stage data | scaling fit and R-08, first point |
| S2 | all S1 + DuckLake spike decision + K-07 decision (by LY-06 bytes decoded per Q2 neighborhood fetch, Q2 latency checked afterwards) + K-12 decision (L5 index, SLO-driven) | 2-point scaling fit |
| S2C | all S2 + null-world calibration at stage scale (R5) | 3-point fit |
| S3 | all S2C, oracle and accretion on a stratified sample | 4-point fit |
| S4 | all S3 + an L2 drift test (ARI vs full recluster on a sample) | drift report |
| S5 | all S4 + L5 decision per class from the scaling fit; null-world calibration at stage scale (R5) | cost per stage |
| S6 | all S5 on the split hosts (per-genome work on wilkcomp02, global work on Hetzner); null-world calibration at stage scale (R5) | cost per stage |
| P1/P2 (Hetzner port, `STAGES.md` §5) | D-05 cross-host, K-* rerun on the new host, SLOs re-baselined | |
| S7+ | all S6 on the off-box host | |

Family S gates (§4.10) are listed by `STAGES.md` numbering in the second table of §4.10, not in this table.

### 4.10 Search tests (family S, from SEARCH_PLAN §7)

Gate types: **B** blocks the stage; **R** is reported in the release notes (equivalent to W); **human** is a sign-off by a lab member. Tiers ST0–ST5 are SEARCH_PLAN §6 (ST0 = S0; ST1 = S1, S2; ST2 = S2C; ST3 = S3, S5; ST4 = S4, S6; ST5 = S7+). "Where" is the search host: the pod lane (§5) through S5, the Hetzner search host from ST4 (S6).

| ID | op (MR) | test | where and when | gate |
|---|---|---|---|---|
| S-01 | O-SEARCH | Replay recall: for each tier's replays, the fraction of published hits (mapped by sequence) that lie in loaded genomes and are recovered. Target ≥ 0.90; every miss is classified (threshold, gene call, absent from the build, error in the published list). Precision is reported, not gated. | search host; at each tier's certifying stage, on that tier's replays (SEARCH_PLAN §6) | B |
| S-02 | O-SEARCH | Cross-tool concordance on 20 fixed seeds: pyhmmer vs MMseqs2 profile vs DIAMOND (and HH-suite for remote cases). Jaccard of representative hit sets at matched bitscore reported; disagreements sampled and adjudicated. | search host; every stage from S1 (first at ST1) | R (report) |
| S-02d | O-SCAN, O-NBR | Independent implementation: Domainator (`domainate.py`, `select_by_cds.py`) on the golden genomes with the same models and score cutoffs; per-CDS domain calls and neighborhood membership compared with L1 and `bio-box nbr`. | pod, ring R3; S0, then nightly and at every container rebuild | B: exact agreement on calls above cutoff + 10 bits; differences near the cutoff reported |
| S-03 | O-SEARCH, O-SCAN (MR-17) | Z invariance: the same recipe at stage k and k+1 gives identical bitscores and identical E at Z_ref for targets present in both. | pod; every stage boundary from S1 (fixed recipe set rerun on the new release) | B |
| S-04 | O-SEARCH (MR-09) | Controls: recipes without a positive control are refused; the positive must be recovered; the shuffled-seed negative returns no more hits than its expectation at Z_ref. | runner unit tests plus the stage harness; every stage from S0 | B |
| S-05 | O-RUN | Provenance: every `search_hit` row resolves to a `search_run`; an injected OOM and an injected timeout are recorded as `failed`, never as zero hits. | stage harness with fault injection; every stage from S0 | B |
| S-06 | O-RUN | Envelope: each question type's cost on the pod is inside its budget in SEARCH_PLAN §5; no search process exceeds 300 CPU-min; `/tmp` quota respected. | the host that serves search (pod through S5, Hetzner from S6); every stage from S0 | B |
| S-07 | O-ITER | Drift alarm: a seed known to corrupt (a promiscuous-domain fusion) triggers the foreign-domain alarm and the revert. | search host; every stage from S0 | B |
| S-08 | O-PROJECT | Projection accuracy: rep-projected hits vs exact re-score on a sample of 10^5 members agree on ≥99% of calls above cutoff + 10 bits; all near-cutoff calls use the exact score. | search host; every stage from S1 (ad hoc final round); library base scan from S3 (§2B.5) | B |
| S-09 | O-NBR (MR-03, MR-04, MR-05) | Neighborhood semantics: synthetic contigs test each edge policy, circular wrap, strand and operon-gap predicates, and censoring flags. | stage harness on synthetic contigs; every stage from S0 | B |
| S-10 | workflow (any manifest) | Human dry run: a lab member runs one synteny hypothesis end to end from a written recipe without help. | search host; once per tier from ST1 (recorded at S1, S2C, S5, S6) | human, once per tier from ST1 |
| S-11 | O-SCAN, O-STATS | Temporal leakage: the B1 `dict@2013` build contains no PADLOC-DB candidate systems or other post-2013 models. | B1 holdout build; at ST3 (S3 and S5 gates) | B (at ST3) |
| S-12 | workflow manifests (TEST_ARCHITECTURE §5.1) | Recipe expressiveness: the three Domainator case studies and five synteny hypotheses of different chains (SYN-143, SYN-114, SYN-162, SYN-254, SYN-343) are written in the recipe language with no escape to ad hoc scripts; each gap found is closed before ST1. | on paper at S0; each recipe executed at its replay tier | B (at S0) |
| S-13 | O-EXPORT | Schema stability: `--json` outputs of `search`, `nbr` and reports match a versioned schema; a removed or renamed key without a version bump fails; added keys pass. | CI on every container rebuild, from S0 | B |

Search gates by stage (`STAGES.md` numbering; C-10..C-12 are in §4.1 and run at every stage):

| stage | tier | must pass (B) | recorded (R / human) |
|---|---|---|---|
| S0 | ST0 | C-10..C-12 (partial pack), S-02d, S-04, S-05, S-06, S-07, S-09, S-12, S-13, S-01 on ST0 replays; R0–R3 green; LY-01..LY-10; micro-corpus built; mutation score measured; null-world calibration first pass (§4.9) | cost table measured on the pod; mutation score and calibration values |
| S1 | ST1 | all S0 + S-03, S-08; R4 | S-02 (first), S-10 (once for ST1) |
| S2 | ST1 | all S1 + S-01 on ST1 replays; R5 (ST1 certifying) | S-02 |
| S2C | ST2 | all S2 + S-01 on ST2 replays; R5 (ST2 certifying) | S-02 (cross-tool concordance required in the S2C gate, SEARCH_PLAN §6.1), S-10 |
| S3 | ST3 | all S2C + S-08 on the library base scan, S-11 | S-02 |
| S4 | ST4 | all S3; search acceptance for S4 waits for ST4 (S6) | S-02 |
| S5 | ST3 | all S4 + S-01 on ST3 replays, S-11; R5 (ST3 certifying) | S-02, S-10 |
| S6 | ST4 | all S5 on the Hetzner search host + S-01 on ST4 replays; R5 (ST4 certifying) | S-02, S-10 |
| S7+ | ST5 | all S6 on the off-box host | S-02 |

Ring rule for this table: "all Sk" inherits R0–R3, R4 and LY-01..LY-10 (rerun on the new stage's data) but not R5. R5 (null-world calibration at stage scale, B1–B6, replays, scaling fit; TEST_ARCHITECTURE §8) is required only where marked, which are the stages where S-01 certifies a tier's replays from ST1 on; at S0 the R5 content is limited to the null-world first pass and S-01 on ST0 replays. Each row is the §4.12 gate query.

### 4.11 Layout tests (family LY, TEST_ARCHITECTURE §7)

Defined in TEST_ARCHITECTURE §7; one line each here. Run on the micro-corpus at S−1/S0 (ring R2) and on each stage's data at its gate (R4). Gate: B at S0 and at every later stage (§4.9).

| ID | op (MR) | test | decides or protects |
|---|---|---|---|
| LY-01 | all writers | every table matches its registered schema: types (16-byte UIDs as FIXED_LEN_BYTE_ARRAY), nullability, enum closure | F-02; extends C-01, C-06 |
| LY-02 | all writers | declared sort order holds within each file | pruning, F-14 |
| LY-03 | all writers | min/max statistics present for every sort-key column in every row group | pruning |
| LY-04 | all writers | file and row-group sizes within the §3.2 bands; no NFS file < 64 MB outside the waivers (§2B.3, §4.12) | NFS rule (§0) |
| LY-05 | O-CATALOG | partition keys in paths equal column values | pruning correctness |
| LY-06 | O-SLICE, O-NBR, O-SERVE | canonical point queries (one contig slice, one neighborhood, one anchor) read ≤ N row groups, counted from `EXPLAIN ANALYZE` / profiling output, not wall time | K-13, row-group sizes, sort keys |
| LY-07 | O-CATALOG | DuckDB, pyarrow and polars read identical rows and checksums from every table | reader interoperability |
| LY-08 | all writers | same inputs written twice give identical bytes (fixed DuckDB version and threads) | F-13; D-01, K-11 |
| LY-09 | O-CATALOG | adding a nullable column leaves older catalogs readable; removing or retyping one fails CI | schema stability |
| LY-10 | O-CATALOG | catalog lists files explicitly (no globs); every manifest file exists with its sha256 | F-11; §3.3 |

LY-06 generalizes K-05 (partition pruning) and K-06 (row-group pruning) from one table each to every canonical query. At S0, LY-06 row-group counts decide K-13 (§2A.1), the row-group sizes and the sort keys (§3.2). Latency (R-08, §4.8) is checked afterwards and does not decide.

### 4.12 Rings and results (TEST_ARCHITECTURE §8)

| ring | content | data | where | cadence | budget |
|---|---|---|---|---|---|
| R0 static | schema registry, recipe and manifest lint; registry completeness; image manifest | none | laptop, CI | every commit | < 1 min |
| R1 unit and property | operation functions; Hypothesis properties; stateful accretion machine on toy data (§4.4) | generated | laptop, CI | every commit | < 5 min |
| R2 micro-corpus | full build; oracles; all MRs; workflow manifests; LY tests; fault injection; mutation testing (weekly) | micro-corpus | CI runner or laptop | every merge | < 30 min |
| R3 golden | golden set with real tools and full model sets; G tests; S-02d | golden set | pod | nightly and every image build | < 6 h at nice 19 |
| R4 stage gate | sampled oracles, A, D, K, LY and S tests, SLOs | stage data | pod or Hetzner | per stage | as in STAGES |
| R5 calibration and benchmark | null-world calibration at scale, B1–B6, replays (S-01), scaling fit | synthetic worlds at stage scale; real | pod or Hetzner | per tier (certifying stages, §4.10) | as in SEARCH_PLAN |

**`test_run` table.** One row per test execution.

| column | type | note |
|---|---|---|
| `test_id` | text | e.g. `C-01`, `LY-06`, `MR-03`, `S-10` |
| `ring` | enum `R0`..`R5` | |
| `operation` | text | O-* from the registry, or `all` |
| `workflow` | text, nullable | workflow manifest name (TEST_ARCHITECTURE §5.1) |
| `params_hash` | 16 bytes | BLAKE3 of the test's parameters |
| `outcome` | enum | `pass`, `fail`, `error`, `skip` (*values proposed here*) |
| `metrics` | JSON | e.g. row groups read (LY-06), recall (S-01), KS statistic (null worlds); signer for human gates (S-10) |
| `duration_s` | float64 | |
| `code_commit` | text | |
| `image_digest` | text | |
| `release_id` | text, nullable | null for R0–R2, which run without a stage release |
| `host` | text | |
| `seed` | uint64 | fixed per test; no retries |
| `started_at` | timestamp | *added here*; orders runs for regression diffs |

**Storage.** Packed Parquet under `$ATD_ROOT/tests/test_run/` (§3.1): one file per R3–R5 invocation, one file per day for all R0–R2 rows pushed from the laptop or CI host. `tests/test_run_index.json` lists the files and is rewritten by atomic rename when a file closes, as `search/run_index.json` is (§2B.3); the catalog exposes `test_run` as a view over that list (LY-10).

**Small-file waiver.** Test result files are often < 64 MB. The waiver is bounded: one file per R3–R5 invocation and one per day for R0–R2, ~10^3 files per year (*estimate*), so `tests/` stays < 10^4 files.

**Gate query.** `gate_spec(stage, test_id, gate)` (`tests/registry/gates.yaml`) is generated from the §4.9 and §4.10 tables, with ring requirements expanded to every registered test of those rings. A stage passes when this returns no rows:
```sql
SELECT g.test_id
FROM gate_spec g
LEFT JOIN test_run t
  ON t.test_id = g.test_id
 AND t.image_digest = $image
 AND (t.release_id = $release OR t.ring IN ('R0', 'R1', 'R2'))
WHERE g.stage = $stage AND g.gate = 'B'
GROUP BY g.test_id
HAVING count(t.test_id) = 0 OR NOT bool_and(t.outcome = 'pass');
```
Regressions are diffs between two runs' rows for the same `test_id` and `params_hash`.

**Tests of the tests** (TEST_ARCHITECTURE §8): mutation score ≥ 80% of mutants killed by R1–R2 (*target*; measured at S0); every S1-severity failure (F-01..F-10) has a test that fails when it is injected; one flexibility check (a lab member adds a workflow from a manifest and recipe only).

---

## 5. Orchestration

**Recommendation: Snakemake (≥ 8)** with a custom executor wrapper and two profiles, `pod` and `hetzner`.

Snakemake over Nextflow:
- **No scheduler is needed.** Snakemake's local executor with `--cores` and `--resources mem_gb=… cpu_min=…` is a natural fit for one shared box. Nextflow's local executor works, but its strength is cluster and cloud executors, which are absent here.
- **Rules are Python, and the test suite and ledger code are Python.** The ledger writer, the STOP check and the manifest generator are plain functions shared between rules and pytest.
- **Content-based rerun triggers** (`--rerun-triggers mtime params input code software-env`) line up with the `annot_version` hash.
- **Its weaknesses are manageable.** The DAG of a 10^6-genome run is too large if each genome is a job, so jobs are shards of 500–2,000 genomes, keeping the DAG at ≤ 10^4 jobs per stage. Snakemake's metadata directory (`.snakemake/`) lives on local `/tmp` or `/stor/work`, never on scratch, because it is many small files.

Required behaviors, all implemented in a single wrapper `bio-box run-step` that every rule calls:

| behavior | implementation |
|---|---|
| resource caps | `systemd-run --user --scope -p MemoryMax=${mem}G -p CPUQuota=${cpu*100}% --` wraps every rule. The profile sets `mem` and `cpu`; the scope's `memory.peak` is written to the ledger |
| nice policy | rule tag `monolithic: true` (MMseqs2, Foldseek, HH-suite, DuckDB global builds) → `nice -n 19 ionice -c2 -n7`; per-genome worker shards → `nice -n 15`. The watchdog then has nothing to change |
| thread hygiene | the wrapper exports `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=NUMEXPR_NUM_THREADS=$threads` and passes explicit `--threads`, so Python tools like MacSyFinder cannot oversubscribe (the problem in `WILKECOMP` §9.3b) |
| day/night profile | the wrapper reads `profile.yaml` (`day: 40 threads`, `night: 56`) by local time and checks `who`/load. If other users' load exceeds 6 cores, it waits before starting a monolithic step. With the search lane on, the build gets 32 threads by day and 48 at night; the lane holds the other 8 |
| interactive search lane (SEARCH_PLAN SD-15) | `bio-box search run` jobs go through the same wrapper with `lane: search`: 8 threads by day, 8 at night, reserved out of the 40/56 caps. Each job is split into per-process shards of < 300 CPU-min (watchdog ~369) at `nice -n 15`; a shard that would exceed 300 CPU-min is refused at plan time, and a step that cannot be sharded that small (e.g. MMseqs2 iterative over S6) runs on the batch side as `monolithic` at nice 19. `/tmp` quota 150 GB per lane under `/tmp/atd/search/`. During L2 barriers (`barrier_L2_<release>`, when the build needs `/tmp`) the lane reads pack shards from NFS instead of copying them. The lane does not wait on other users' load the way batch monolithic steps do: its 8 threads are already reserved, so it starts at once. It honors the STOP file between shards (exit 75; finished shards are kept and the run resumes on rerun). `serve` submits and polls lane jobs; it never waits on them |
| STOP file | before each job and between chunks of chunked tools, the wrapper checks `$ATD_ROOT/STOP`; if present it exits with code 75 (EX_TEMPFAIL). Snakemake marks the job incomplete and rerun resumes |
| local staging | each job works in `/tmp/atd/<job_uid>/`, writes outputs there, verifies the checksum, then does **one** `rsync --checksum` to NFS and an atomic rename into the final path. Partial outputs never appear on NFS |
| small-file guard | the wrapper refuses to rsync a directory with > 1,000 files; the job must pack into Parquet or tar.zst first |
| ledger | the wrapper appends the ledger row (in/kept/rejected, RSS peak, CPU-seconds, wall time, versions) for every job |
| Hetzner switch | `--profile profiles/hetzner` changes: root path, thread caps (all cores), memory cap, NVMe local = final storage (no rsync hop, staging is a same-filesystem rename), and an optional object-store sync for the system-of-record copy. Rules are unchanged |

---

## 6. Portability

| item | pod | Hetzner |
|---|---|---|
| runtime | Apptainer/Singularity `.sif` built from the same OCI image (singularity is installed on the pod; Docker group membership is not assumed to be for general use) | Docker (or Apptainer, for parity tests) |
| image | `ghcr.io/<org>/bio-box:<version>-cpu` (DEPLOYMENT §2.1), AVX2 baseline; no CUDA in the CPU image | same digest; an optional AVX-512 build is a separate tag, and D-05 must pass across tags |
| environment lock | `pixi.lock` for bioconda/conda-forge binaries (MMseqs2, HMMER, Infernal, pyrodigal, geNomad, CheckM2, …) plus `uv.lock` for pure-Python packages; both locks hashed into `annot_version` | same locks |
| reference DBs | not in the image; fetched into `l0/refdb/<name>/<version>/` with sha256 and treated as L0 sources | same |
| DuckDB | pinned in `pixi.lock`; the Python package and CLI must report the same version (K-01) | same |
| host probe | `bio-box doctor` records CPU flags, cores, RAM, cgroup limits, `ulimit -v/-n`, filesystem types of `$ATD_ROOT` and `/tmp`, and DuckDB version into the ledger at the start of every run | same |

The image is built once per version from a recipe, and its digest is the identity. The pod pulls the `.sif` by digest. Nothing is installed ad hoc on either host.

---

## 7. Open items for the staged roadmap to own

- Thresholds marked *(tunable)* need values from S1/S2 measurement: the recluster trigger %, ARI floor, and C-09 bands.
- Which stage first introduces geNomad/CheckM2 (needs their DBs in L0) depends on runtime measured on Haswell.
- Whether `contig_tapes` stays `LIST<T>` or becomes an exploded slot table is decided at S2 by K-07 (LY-06 bytes decoded per Q2 fetch; latency checked afterwards). Whether Q2 needs an L5 offset index is decided by K-12 against the SLO.
- DuckLake adoption is decided at S2 from the §3.3 spike criteria.
- The Z constants (§2B.4) and stage S2C await user approval (SEARCH_PLAN §9 items 1–2). NR90/NR50 counts and all search-pack sizes in §2B.7 are estimates until S0–S3 measure them.
- Test-architecture items owned in TEST_ARCHITECTURE §11: (1) Zombi fit at micro-corpus size, fallback ALF (S−1); (2) CI host for R2 (S−1); (3) mutation-score target confirmed after the first S0 measurement; (4) who runs the flexibility check, on which workflow.
- LY-04 waivers: resolved in TEST_ARCHITECTURE §7 (metadata files declared `small_ok: true` in the schema registry, each bounded to a declared count per release). The registry entries are written at S−1.
- `gate_spec` (`tests/registry/gates.yaml`), quarantine (`tests/registry/quarantine.yaml`) and the `test_run` columns from §4.12 are adopted in TEST_ARCHITECTURE §8.

Sources: [DuckDB concurrency](https://duckdb.org/docs/2.0/connect/concurrency.html), [DuckDB Parquet overview](https://duckdb.org/docs/2.0/data/parquet/overview.html), [DuckDB hive partitioning](https://duckdb.org/docs/2.0/data/partitioning/hive_partitioning.html), [DuckLake FAQ](https://ducklake.select/faq.html), [DuckLake docs](https://www.ducklake.select/docs/stable), pod envelope from `WILKECOMP.org` and `fragment-defense/old-plan.2026-10-18.bak` §9.
