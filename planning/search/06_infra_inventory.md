# 06. What the infrastructure plan already provides for ad hoc sequence search

Status: audit only, 2026-10-29. No plan file was edited. Every claim cites `file:line` from the plan as it stands. Paths are relative to `all-the-data/`. Abbreviations: `spec` = `00-spec.md`, `PLAN` = `PLAN.md`, `STAGES` = `STAGES.md`, `DEPLOY` = `DEPLOYMENT-AND-MALLEABILITY.md`, `PRIOR` = `PRIOR-ART.md`, `DL` = `planning/data_layer_and_tests.md`, `SZ` = `planning/sources_and_sizing.md`, `PF` = `planning/preflight_gaps.md`, `SYN2` = `audits/SYNTHESIS_2.md`, `deck` = `presentations/all-the-data-deck.html`. Question-type codes (Q-FAM … Q-SYS) come from `planning/search/00_brief.md:62-79`.

Labels: **plan** = stated in a plan file; **derived** = arithmetic on plan numbers, formula shown; **estimate** = this audit's own guess, with the reasoning given; **recalled** = a figure from outside the corpus that was not re-checked.

**Headline.** The plan's word "search" means querying the precomputed tapes, not running a sequence search. Examples: `spec:73` ("searchable as raw neighbourhoods"), `DEPLOY:97` ("search-only delta"), `DEPLOY:118`. The plan stores everything a search would need: every protein, contig and IGR (`spec:43`, `DL:72-85`), and the binaries for HMMER, Infernal, MMseqs2 and Foldseek (`DEPLOY:83`). What it lacks is a command, persistent per-release target databases, a defined E-value regime for user searches, and any path for user hits back into tokens, `slot_pairs` or `pair_stats`.

---

## 1. Capability inventory

### 1.1 Sequence store (L0-seq) and retrieval

| Capability | What the plan says | Evidence |
|---|---|---|
| Invariant | Every contig (dedup members and screening rejects included), every protein, and every IGR as coordinates. Retrievable by UID and by `uid:start-end:strand`. UIDs are BLAKE3-128. Any project may read the store directly. | `spec:43`; `PLAN:347`; `STAGES:69` |
| Exception | Human or vertebrate-host contigs keep only their digest and reason code | `STAGES:69`; `PF:32` |
| Contig storage | bgzip FASTA shards with `.fai`/`.gzi`, records named by `raw_digest` (DEC-17, 2026-10-30; was `seq_chunks` Parquet) | `DL` §2A.1 |
| Occurrence | `contig_occurrence` (source, accession, header, `raw_digest`, `contig_uid`) | `DL:81` |
| Proteins | `protein_seq(protein_uid, protein_digest, length, seq, gene_call_version)`, sorted by `protein_uid`; U/O and frameshift-stitched sequences stored as reconstructed | `DL:83` |
| Slot → protein | `slot_protein (slot_uid, gene_call_version) → protein_uid`, plus code and edit type | `DL:83` |
| IGRs | Coordinates only: `igr(contig_uid, start, end, strand, left_slot, right_slot, igr_version)` | `DL:84`; `PF:50` (C5 resolved) |
| Tool-interop FASTA | bgzip, `.fai` and `.gzi`, **generated on demand** into `/tmp` or a project dir. "Nothing downstream depends on a persistent FASTA copy." | `DL:85`; `DEPLOY:193`; `PLAN:350` |
| Physical layout | `l0seq/seq_chunks/stage=/digest_prefix=/…`, `protein_seq/gene_call_version=/stage=/protein_prefix=/…`, `igr/…`, `seq_locator.parquet` | `DL:96-102` |
| Alternative layout (conflicts with the row above) | `sequences/contigs/part-<shard>.parquet (contig_uid, seq zstd)`, `proteins/…`, `igrs/…`, `fasta/<kind>/…` | `DEPLOY:189-194` |
| CLI: point lookup | `bio-box seq get --contig/--slot [--flank N] [--nt/--aa]/--protein/--igr` | `DL:120-123`; `DEPLOY:103` |
| CLI: bulk export | `bio-box seq extract --uids …/--sql "SELECT …" --out out.fa.gz [--bgzip --faidx] [--members/--reps-only]`. Requests are sorted by `(file,row_group)`. | `DL:124,127`; `DEPLOY:104` (spells the flag `--query <view> --kind <contig/protein/igr>`) |
| CLI: lineage | `bio-box seq lineage <uid/accession>` prints the chain slot → contig → source → L0 object → contract | `DL:114,125`; `PLAN:352` |
| DuckDB | Macros `seq_slice`, `revcomp`, `seq_get`; `contig_index`; views `v_seq_contig/v_seq_protein/v_seq_igr`; Python `atd.seq.get` | `DL:131-142`; `DEPLOY:209` |
| Latency SLO (R-08) | One slice ≤50 ms warm / ≤500 ms cold; 5 Mbp contig ≤1 s; 10⁴ slices ≤30 s; **10⁶ proteins ≤2 min** | `DL:172` |
| Integrity tests | R-01..R-07 and R-09 block S0; R-03 compares 10k random slices with `samtools faidx` | `DL:165-173`; `STAGES:159` |
| Format decision | K-13 bake-off at S0: Parquet chunks vs bgzf | `DL:93,174`; `STAGES:137,373` |

### 1.2 What is precomputed per protein, segment and IGR

| Item | Plan | Evidence |
|---|---|---|
| HMM libraries | Pfam, NCBIfam, PHROGs, arCOG, each with its own GA/TC cutoffs. COG comes via eggNOG. CDD only via RPS-BLAST, "if used". There is no merged HMM set. | `PLAN:278-281`; `spec:54`; `SYN2:35` |
| PHROGs scope | Run only on viral, plasmid and prophage proteins | `SZ:31` |
| Fixed Z | "a fixed, recorded HMMER `Z`". Elsewhere "fixed per release and passed explicitly". The value and definition (sequences in which DB?) are not stated. | `spec:54`; `PLAN:118,279`; `DEPLOY:162`; `deck:256` |
| Hit table | `slot_annotations(slot_uid, token_id, evalue, bitscore, env_start, env_end, cutoff_pass, library_version)` holds all hits. The tape column holds one canonical pick. | `PLAN:369` |
| Reporting threshold for stored hits | Not stated. `cutoff_pass` implies that some sub-GA hits are kept, but no inclusion E-value is given. | `PLAN:369` |
| Unit actually scanned | **Inconsistent:** per genome in L1 (`DL:36`); "NR50 residues × models" (`PLAN:118`); "Pfam on NR90" and "Pfam on new NR90 reps only" (`SZ:85,86`) | see §3 |
| Segments | `(orf_idx, seg_idx)` covering nested hosts, circular permutations, collapsed repeats (gap ≤35 aa) and CTD tips. Pfam `NE` and ECOD boundaries; heuristics flagged. | `spec:54`; `PLAN:285-288`; tape columns `PLAN:396-400` |
| `segment_uid` | `blake3(protein_uid, env_start, env_end, segmentation_version)` | `DL:57` |
| Subfamily | `subfamily_50_id`: 50% ID / 80% coverage clusters of excised envelopes within each `seq_family_id`. Key = digest of the representative envelope. | `PLAN:289-290`; `spec:22` |
| Catalytic motif | `catalytic_motif_packed` (≤12 residues), `site_source` ∈ {M-CSA, CDD_site, conservation_only}, computed per member; `catalytic_intacts` | `PLAN:291-295,412-413` |
| Badge | TM count, signal type, 8 sorting tags, pI, 8 fixed motif flags, length, repeat count, edit flags | `PLAN:315-327` |
| Structure | ECOD via Foldseek against distributed 3Di (AFDB, ESMAtlas30). Transfer only at ≥90% ID, pLDDT ≥70 over ≥70% of the envelope. `struct_fold_ids`, `fold_sources`, `fold_confidences`. | `spec:23`; `PLAN:272-275,302,414-416`; `PLAN:121` (`mmseqs search` vs AFDB/ESMAtlas30) |
| ProstT5 | GPU stage only, for recurrent dark NR50 segments. No predicted 3Di on CPU-only stages. Runs after a global barrier. | `PLAN:297-301,312`; `spec:23,31`; `PLAN:81` |
| Remote homology | Stage 2d: `mmseqs result2profile` + `hhalign` against candidates nominated by Foldseek or profile search, for multi-member dark clusters. Rate unmeasured. No output table named. | `PLAN:307-309` |
| ncRNA | IGRs ≥40 bp + 25 bp CDS flanks → IGR90. CMs: Rfam prokaryotic/viral subset + `cmcalibrate`d orphan CMs. `cmsearch` with native filters; `--cut_ga` for Rfam only. Pocket subtype per member after `cmalign`. | `PLAN:239-242`; `spec:52`; `PLAN:120` |
| CRISPR, tRNA, integrons | MinCED + CRISPRDirection, spacers → `spacer_protospacer_graph.parquet`; IntegronFinder 2; tRNAscan-SE, Aragorn; TransTermHP as a score | `PLAN:243-247` |
| smORFs | 15–75 aa candidates; promoted only if present in ≥3 orders with RNAcode or dN/dS < 0.5. Unpromoted candidates go to a sidecar. | `PLAN:228-236` |
| Purged and split ORFs | Purged shadow ORFs stay in the rejection ledger with coordinates. Unsupported split-gene merges are flagged `CANDIDATE_INDEL`. | `PLAN:248-258` |
| Evidence level | `evidence_level ∈ {curated_model, profile_hit, promoted_smorf, rule_merge, candidate}` | `PLAN:260,417` |
| External closure | PDB, Swiss-Prot, PaperBLAST and FEBA hits, envelope- and ATOM-gated. **The search tool and target DBs are not named.** | `spec:71`; `PLAN:452,466-470` |

### 1.3 Clusters (L2)

| Item | Plan | Evidence |
|---|---|---|
| NR90 → NR50 | MMseqs2, incremental with `clusterupdate` (MMseqs2 ≥18-8cc5c). Terminal-window guard for proteins ≥350 aa. | `PLAN:267-269`; `DL:37` |
| Member alignments | "Every member keeps its own sequence and its alignment to the representative (`protein_member` table)". `scan-new-models` projects hits "through the stored cluster alignments". | `PLAN:270-271,376`; `DEPLOY:161` |
| L2 layout | `l2/release=<R>/{nr90,nr50,subfam,igr}/…; cluster_member/; id_history.parquet`. **No alignment table, no `protein_member`.** | `DL:188`; `DL:37` |
| Membership size | NR90/NR50 membership Parquet at 16 B/protein, about 8 GB at S5; MMseqs2 DBs about 300 GB, **transient** on `/tmp` | `SZ:96` |
| Cluster UID | **Three definitions:** digest of the representative sequence (`PLAN:359`); representative `protein_uid` + level (`DEPLOY:130`); `blake3(sorted member UIDs at creation, params)` (`DL:37,59`) | §3 |
| `id_history` | same / split / merged / retired (/ renamed) | `DL:68`; `DEPLOY:134`; `spec:62` |
| Recluster triggers | >20% new sequences, a parameter change, or a failed drift test; crosswalk emitted | `DL:37` |
| IGR90 | 90% nt identity; "only chooses what to scan" | `PLAN:239`; `spec:52` |
| Contig dedup | `contig_membership` (`PLAN:205`; `spec:48`), also called `contig_member` (`DL:82`) | |
| Vehicle clusters | MinHash-LSH over ±15-ORF flanks | `PLAN:410,125` |
| Dark tokens | "release-scoped for dark clusters (resolved through `cluster_uid`)" | `DL:60` |

### 1.4 Tapes, `slot_pairs` and `pair_stats` (L3/L4)

| Item | Plan | Evidence |
|---|---|---|
| `contig_tapes` | One row per contig with parallel `LIST` columns: `slot_uids`, `orf_idxs`, `seg_idxs_in_orf`, `parent_seg_idxs_in_orf`, `domain_topology_flags`, `junction_topology_prevs`, `directon_ids`, `regulon_locus_ids`, edge distances, regimes, `vehicle_cluster_ids`, `clan_or_cog_ids/seq_family_ids/subfamily_ids`, motifs, folds, badges, `protein_uids`. Contig scalars include `tax_lineage` and confidence. | `PLAN:383-421` |
| Junction topology | 0 intra-ORF fusion, 1 co-directional, 2 divergon, 3 convergon/excludon, 4 sense–antisense nested, 5 inter-locus break | `PLAN:404`; `spec:56` |
| Tokens | `token_id = level(3) / source(5) / local(24)`. For Pfam/Rfam, `local` comes from a **frozen dictionary**. Token columns are a release-time join, never stored on the coordinate tape. | `PLAN:369`; `DL:60,38` |
| `slot_pairs` | Directed pairs within ±k ORFs, **or fused**, inside one `regulon_locus_id`. Columns: family and subfamily tokens, `orf_offset`, `junction_topology`, `delta_bp`, `tax_order`, vehicle, bioproject, `partner_badge`, slot UIDs. Partitioned by anchor bucket, with an offset table. | `PLAN:423-427` |
| Naming drift | `window_postings` (`STAGES:171`; `SZ:81`), `pair_postings`/`l3/pairs` (`DL:38,191`), `pair_summary` (`DL:192`) vs `slot_pairs`/`pair_stats` (`PLAN:423,429`) | |
| `pair_stats` | Score test `Z`, `beta_hat`, sufficient statistics per regime × tier, Δbp histogram, **topology counts**, phylum/order bitmaps, genus HLL, entropies, `p_value/local_fdr/q_hier`, `powered`, `geometry_z`, `rank_sd`, `coupling_score`, closure, staleness | `PLAN:437-454` |
| Discovery unit | Module, not pair; primary level `subfamily` | `PLAN:431` |
| Pangenome | Three-valued planes per token, Q5 NOGD, `transition_edges`, Pagel LRT | `PLAN:473-478` |
| Gene trees | FAMSA + VeryFastTree on demand for the top-N pairs, ≤5,000 leaves | `PLAN:127,453`; `SYN2:40` |

### 1.5 Malleability commands

| Command | What it does | Search relevance | Evidence |
|---|---|---|---|
| `bio-box append <fasta>` | Runs the per-genome pipeline, dedups against the base, assigns clusters with `clusterupdate` semantics, appends to the store. Writes deltas under `epochs/<release_id>/delta_<n>/`. | Adds *genomes*, not hits | `DEPLOY:97,137-144` |
| `register-sidecar` | Parquet keyed by `contig_uid`, `protein_uid`, `cluster_uid`, `token_uid`, a token pair, or `(contig_uid, slot_ordinal)`, plus `key_release`. Coverage: warns below 95%, refuses below 50%. Joined **after** the main filter. | The only declared route for a user hit list; see §3 on coverage | `DEPLOY:98,146-151`; `DL:68` |
| `scan-new-models <models> --z` | New HMMs or CMs vs stored **NR50 reps**. Hits projected to members through stored alignments; near-cutoff hits re-checked on members. Z / Infernal DB size fixed per release. Output: a sidecar or delta with library, version, cutoff type and Z. | Closest thing to an ad hoc family search, but hits do not become tokens | `DEPLOY:100,159-164`; `SYN2:35` |
| `remap-token --merge` | Merges only, via an overlay; additive counts re-summed, others marked stale. Splits require a release. | Can pool existing tokens; cannot create a token from a hit list | `DEPLOY:99,153-157` |
| `rescore --params` | Recomputes Z/FDR from stored sufficient statistics. Ships in `lite`. | Statistics only | `DEPLOY:101,166-170`; `PLAN:443` |
| `compact` | Promotes deltas into a new release | | `DEPLOY:102` |

### 1.6 Catalog, serve, query classes

| Item | Plan | Evidence |
|---|---|---|
| Catalog | Generated from `manifest.json`: views and macros over explicit file lists, no globs. Read-only, copied to local `/tmp`, swapped by `rename(2)`. One writer. | `DL:212-217`; `DEPLOY:172-177` |
| `bio-box serve` | `127.0.0.1:8420`, reached over an SSH tunnel. Parameterized `/v1/...` query classes; timeouts, row limits, per-client concurrency. No raw SQL by default; the optional console runs with `enable_external_access=false` and `lock_configuration=true`. | `DEPLOY:57,96,268-275`; `STAGES:78` |
| Interactive DuckDB | Read-only catalog with daytime `threads`/`memory_limit` | `DEPLOY:58`; `DL:16` |
| Query classes (PLAN) | Q1 count, Q2 neighbourhoods, Q3 anchor profile, Q4 drill-down, Q5 NOGD, Q6 co-transmission, Q7 sequence slice, Q8 ad hoc SQL (no SLO) | `PLAN:96-105` |
| SLOs (PLAN, warm p50/p95, cold p95) | Q1 0.5/2/10 s; Q2 0.3/1/5 s; Q3 0.1/0.5/3 s; Q4 0.5/2/10 s; Q5 0.2/1/5 s; Q6 1/5/20 s; Q7 0.05/0.2/1 s; Q8 none | `PLAN:98-105` |
| Query classes (DL, **different**) | Q1 neighbourhood fetch, Q2 anchor scan, **Q3 gapped pattern `[A,B,C]`**, Q4 drill-down, Q5 pathway hole, Q6 co-transmission, **Q7 closure lookup**, Q8 ad hoc | `DL:262-271` |
| SLOs (DL, p95 warm/cold) | Q1 1/5 s; Q2 median token 2/10 s; Q2 hub 20/60 s; Q3 gapped 5/30 s; Q4 2/10 s; Q5/Q6 1/5 s | `DL:334-341` |
| Lab questions | Also numbered Q1–Q10. Q1 = bait-family ±10-gene neighbourhoods with FASTA; Q5 = subfamily export for SSN/GNN. | `PF:63-80`; `deck:640-647` |
| **No sequence-search query class** | None of Q1–Q8 takes a sequence, HMM, MSA, structure or CM as input | `PLAN:96-105`; `DL:262-271` |

### 1.7 Run-step wrapper, image and envelope hooks

| Item | Plan | Evidence |
|---|---|---|
| `bio-box run-step` | `systemd-run --scope MemoryMax/CPUQuota`; `monolithic: true` → nice 19 + ionice; fan-out → nice 15; thread env vars; day/night caps; **waits if other users' load > 6 cores**; STOP file → exit 75; `/tmp/atd/<job_uid>/` staging; one `rsync`; **refuses to rsync >1,000 files**; ledger row | `DL:374-386`; `DEPLOY:36-43` |
| HH-suite marked monolithic | Listed with MMseqs2, Foldseek and DuckDB | `DL:379` |
| CPU image | pyhmmer/HMMER 3.4 (so `phmmer`, `jackhmmer`, `hmmsearch`, `nhmmer` binaries), Infernal 1.1.5, MMseqs2 ≥18, Foldseek CPU, FAMSA, VeryFastTree, MacSyFinder, MinCED, … **No HH-suite, BLAST+/RPS-BLAST, DIAMOND, eggNOG-mapper, PADLOC, CMfinder, CRISPRCasFinder/Typer.** | `DEPLOY:83` |
| CUDA image | ProstT5 runtime, Foldseek GPU | `DEPLOY:84` |
| No ad hoc installs | "Nothing is installed ad hoc on either host" | `DL:401` |
| Data pack | "pressed HMMs, calibrated CMs, Foldseek DBs and model weights" under `ref_models/<library>/<version>/` | `DEPLOY:71,187` |
| Alternative refdb path (conflicts with the row above) | `l0/refdb/<name>/<version>/` | `DL:397` |
| refdb inventory | Pfam, NCBIfam, PHROGs, arCOG, eggNOG, Rfam, orphan CMs, MacSyFinder/PADLOC, geNomad, CheckV, CheckM2, GUNC, FCS-GX, GTDB, iPHoP, Foldseek AFDB/ESMAtlas30 ("hundreds of GB"), truth sets. **No PDB seqres, Swiss-Prot, UniRef, pdb70/HHM DBs, CDD.** | `PF:84-108` |

### 1.8 Editions

| Edition | Stages | Ships | Evidence |
|---|---|---|---|
| `lite` | S0–S5 | Contig tapes, cluster dictionaries, pair stats with sufficient statistics, core plane over species reps, and the sequence store for its stages. IMG/VR (S4) is statistics only, with sequences as fetch recipes (`fetch-required`). Size about 0.3–1.0 TB. Runs on "Laptop, wilkcomp02, small VM". | `STAGES:280,346,350`; `DEPLOY:250,305,309` |
| `core` | S0–S6 | Adds multi-strain planes and quality masks; about 4–6 TB; Hetzner | `STAGES:300,347`; `DEPLOY:251,306` |
| `planetary` | S7+ | Hetzner | `STAGES:348`; `DEPLOY:307` |
| Search DBs in any edition | **Not listed.** Only "Foldseek DBs" in the generic pack definition (reference DBs, not the corpus). | `DEPLOY:71` |

---

## 2. Coverage by question type

Served = a stated table, command or query class answers it. Partial = possible through Q8 SQL or `seq extract` plus a tool run by hand, or only for library tokens. Missing = no stated path.

### Q-FAM: census of a known family with a fixed model
- **Served, for the four libraries.** `slot_annotations` holds every hit with E-value, bitscore, envelope and `cutoff_pass` (`PLAN:369`). GA/TC cutoffs are per library (`PLAN:278-281`). Taxonomy is on `contig_tapes` (`PLAN:388`). A count with geography or taxonomy filters is Q8 SQL, or Q1 count for pairs (`PLAN:98`).
- **Partial.**
  - PHROGs hits exist only on viral, plasmid and prophage proteins (`SZ:31`).
  - COG comes via eggNOG (`PLAN:279`), but eggNOG-mapper is not in the image (`DEPLOY:83`).
  - The stored threshold below GA is unstated, so "GA plus relaxed E-value" censuses may be impossible (`PLAN:369`).
  - If scans run on NR50 or NR90 representatives (`PLAN:118`; `SZ:85-86`) rather than per genome (`DL:36`), member-level census counts depend on projection through alignments that the L2 layout does not store (`DL:188`).
- **Missing.**
  - TIGRFAM appears only as part of NCBIfam (`PLAN:278`).
  - CDD (`PLAN:280`): RPS-BLAST is not in the image.
  - Any other library (KOfam, SMART, PANTHER) needs `scan-new-models` (`DEPLOY:159-164`), whose hits are sidecars, not tokens (`DEPLOY:164`).

### Q-ITER: iterative profile expansion (jackhmmer, PSI-BLAST, MMseqs2 iterative)
- **Served:** nothing.
- **Partial.**
  - `jackhmmer` and `mmseqs search --num-iterations` are in the image (`DEPLOY:83`).
  - The target has to be produced with `bio-box seq extract --kind protein` (`DEPLOY:104`; `DL:124`). At the R-08 rate (10⁶ proteins ≤ 2 min, `DL:172`), extracting the S5 set of about 4.8×10⁸ proteins takes about 16 h (derived: 480 × 2 min).
  - `--reps-only` (`DL:124`) is ambiguous: it does not say whether it means contig-dedup reps or NR90/NR50 reps.
- **Missing.**
  - A persistent NR90/NR50 FASTA or pyhmmer digital DB.
  - A per-release MMseqs2 search DB. The only MMseqs2 DBs are transient (`SZ:96`).
  - PSI-BLAST/DELTA-BLAST (no BLAST+ in the image, `DEPLOY:83`).
  - A rule for which iteration's hits count, and how Z is set (§3).

### Q-REMOTE: profile–profile and remote homology (HHpred, HHblits, hhsearch)
- **Served:** nothing.
- **Partial.**
  - Stage 2d (`PLAN:307-309`) runs `result2profile` + `hhalign` internally for dark clusters, but its output table is unnamed and its rate unmeasured.
  - HH-suite is tagged monolithic in the wrapper (`DL:379`) but **is absent from the image** (`DEPLOY:83`).
- **Missing.**
  - HHM databases: no pdb70, Pfam-HHM, ECOD-HHM or UniRef30 in `refdb.yaml` (`PF:84-108`).
  - An MSA or HMM store for corpus clusters or user profiles.
  - A command to run a user profile against cluster profiles.

### Q-STRUCT: structure search
- **Served:** lookup by fold. `struct_fold_ids`, `fold_sources` and `fold_confidences` on tapes (`PLAN:414-416`) answer "which slots carry ECOD fold X", for transferred 3Di only.
- **Partial.**
  - Foldseek CPU is in the image (`DEPLOY:83`), and AFDB/ESMAtlas30 3Di DBs are refdb entries (`PF:107`). A user can therefore search a query structure against AFDB.
  - Mapping AFDB hits back to corpus proteins needs the AFDB accession behind each transfer. Tapes store only fold ID, source and confidence (`PLAN:414-416`), so the accession is not stored.
- **Missing.**
  - A Foldseek DB of the corpus (NR50 3Di).
  - CPU stages carry no predicted 3Di (`spec:23`; `PLAN:81`). ProstT5 is GPU-only, after a barrier (`PLAN:297-312`).
  - No GPU on wilkcomp02 (`DEPLOY:28`).
  - Single-query ProstT5 on CPU is not addressed. Only bulk use is excluded (`DEPLOY:87`).

### Q-ARCH: domain architecture and fusion (Rosetta stone)
- **Served, for library-annotated families.**
  - `slot_pairs` includes fused pairs with `junction_topology = 0` (`PLAN:404,424-425`).
  - `pair_stats` carries topology counts and fractions (`PLAN:443,446`), so "fused in some genomes, adjacent in others" is countable per token pair.
  - Segment order comes from `orf_idxs/seg_idxs_in_orf/parent_seg_idxs_in_orf/domain_topology_flags` (`PLAN:396-400`).
  - Gene-call error guards: `CANDIDATE_INDEL` (`PLAN:258`) and `evidence_level` (`PLAN:260`).
- **Partial.**
  - The tape carries only the canonical pick under a precedence (`PLAN:369`). Architectures built from overlapping libraries need `slot_annotations` through Q8.
  - The gapped `[A,B,C]` pattern class exists only in DL (`DL:266,339`), not in PLAN's classes (`PLAN:96-105`).
- **Missing, for custom-HMM hits.** There is no path into `slot_pairs`:
  - tokens come from a frozen dictionary (`DL:60`);
  - `scan-new-models` writes a sidecar or delta (`DEPLOY:164`);
  - sidecars are joined after the main filter (`DEPLOY:150`).

  The only partial route is `remap-token --merge`, when the hit set is exactly a union of existing tokens (`DEPLOY:155-156`).

### Q-NBR: neighbourhood, guilt by association
- **Served, for token anchors.** Q2 tape fetch (`PLAN:99,421`), Q3 anchor profile (`PLAN:100`), Q4 drill-down (`PLAN:101`), lab Q1 (`PF:69`). Neighbours carry tokens, vehicle, regime, badge and edge distances (`PLAN:404-419`). Enrichment comes as `Z`/`q`/`powered` (`PLAN:441-448`).
- **Partial, for user hits.**
  - A user can map a hit list `protein_uid → slot_uid` via `slot_protein` (`DL:83`) and fetch tapes (Q2), so a neighbour census is possible.
  - Nothing states a sort order or index on `slot_protein` for reverse lookup by `protein_uid`.
  - Q2 is keyed by `dense_contig_id` (`PLAN:421`), which is release-local (`DEPLOY:132`).
  - `slot_pairs` stops at `regulon_locus_id` (`PLAN:424`), so a ±10-gene census (`PF:69`) needs Q2, not `slot_pairs`. ±k is a parameter (`PLAN:459`) whose value is not stated.
- **Missing.** No enrichment statistics with a user hit set as anchor. `pair_stats` exists only for token pairs (`PLAN:431`).

### Q-PROF: phylogenetic profile, co-occurrence, co-evolution
- **Served, for tokens.** Phylum/order bitmaps and genus HLL in `pair_stats` (`PLAN:443`); pangenome planes and Q5 (`PLAN:473-476,103`); `transition_edges` and Q6 (`PLAN:477,104`); Pagel LRT (`PLAN:476`). Stage availability: sparse planes at S3 (`STAGES:220`), core plane at S5 (`STAGES:259`), three-valued planes at S6 (`STAGES:288`).
- **Partial, for user hits.** Presence by taxon is computable with Q8 from `contig_tapes.tax_lineage` (`PLAN:388`), but without the Beta-binomial shrinkage, masks or power columns (`PLAN:473-475`).
- **Missing.** Mirror-tree and residue co-evolution (DCA). Gene trees serve only homoplasy for top-N pairs (`PLAN:453`).

### Q-RNA: structured RNA motif discovery and search
- **Served.** Rfam prokaryotic/viral subset plus orphan CMs on IGR90 with native filters; per-member pocket subtypes (`PLAN:239-242`); production-mode recall checked against full Rfam at S0, target ≥0.98 (`STAGES:148`; `SZ:130`). Lab Q8, riboswitch upstream of a transporter (`PF:76`).
- **Partial.**
  - A user CM can go through `scan-new-models` (`DEPLOY:161`). Its target is stated as "NR50 representatives", which is wrong for CMs, and an IGR90 target is not named.
  - IGR90 sequences must be sliced each time, since IGRs are coordinates only (`DL:84`).
- **Missing.**
  - CMfinder or any motif-discovery tool (not in the image, `DEPLOY:83`).
  - Search of ncRNAs in IGRs under 40 bp, inside CDS or antisense. The scan unit is IGR only (`PLAN:239`).
  - IGR90 FASTA or a persistent CM target DB.
  - A conflict on IGR bounds: "+25 bp CDS flanks" (`PLAN:239`) vs "coordinates lie strictly between their flanking slots" (`DL:169`).

### Q-NT: nucleotide-level search (tblastn, unannotated ORFs, nhmmer, repeats, CRISPR)
- **Served.**
  - Raw contigs, members included, are retrievable (`spec:43`; `DL:80`).
  - CRISPR arrays and spacers (`PLAN:244`); spacer–protospacer links at S4 (`STAGES:242`).
  - smORF candidates in a sidecar (`PLAN:236`).
  - Purged ORFs with coordinates in the ledger (`PLAN:252`).
- **Partial.** `nhmmer` is in HMMER (`DEPLOY:83`), but a contig target means `seq extract --kind contig` of 400–450 Gbp at S5 (`SZ:85`).
- **Missing.**
  - Six-frame or unannotated-ORF protein search over contigs: no tblastn, no `mmseqs search --search-type 2/4` DB.
  - Persistent nucleotide search DBs.
  - Tools for spacer-to-protospacer matching (not in the image).
  - PRIOR-ART excludes SRA-scale nucleotide search (`PRIOR:65`).

### Q-MOTIF: catalytic residue or motif filtering on hits
- **Served, for library models with curated sites.** `catalytic_motif_packed`, `site_source`, `catalytic_intacts` (`PLAN:291-295,412-413`); eight fixed badge motif flags (`PLAN:323`).
- **Partial.** A regex over `protein_seq` via Q8 (`DL:83`; `PLAN:105`) is a full scan with no SLO. At S5, `protein_seq` is about 45 GB (`DL:155`, old stage label) or about 66–74 GB (derived from `SZ:58`).
- **Missing.**
  - Residue-at-column extraction for user HMM hits: `scan-new-models` output lists provenance only, with no alignments (`DEPLOY:164`).
  - `CDD_site` needs RPS-BLAST, which is not in the image (`PLAN:280,293`; `DEPLOY:83`).

### Q-SSN: sequence similarity and genome neighbourhood networks
- **Partial.** Lab Q5 exports a subfamily's members with neighbourhoods, sequences and taxonomy (`PF:73`; `deck:644`). EFI-GNT is the named baseline (`PRIOR:32`).
- **Missing.**
  - All-vs-all edge computation or storage.
  - Envelope-level export: `seq get/extract` have no `--segment` option (`DL:120-124`), although `segment_uid` exists (`DL:57`).
  - Persistent alignments to representatives (§3).

### Q-CLUST: corpus-scale clustering for new families
- **Served.** NR90/NR50 (`PLAN:267`), `subfamily_50` (`PLAN:289`), IGR90 (`PLAN:239`), vehicle clusters (`PLAN:410`), crosswalks (`DL:68`), the novel-NR50-per-Tbp rarefaction gate (`STAGES:312`), and cross-checks against Pavlopoulos/FESNov (`PRIOR:40-41`).
- **Partial.** Remote-homology grouping of dark clusters (Stage 2d, `PLAN:307-309`), with no output table.
- **Missing.**
  - Structural clustering of the corpus (no corpus 3Di on CPU, `spec:23`).
  - Profile-based deep clustering (e.g. 30% or HMM–HMM).
  - Sensitive re-clustering of a user's hit set as a first-class operation.

### Q-SYS: rule-based system detection
- **Served.** DefenseFinder, PADLOC, CONJScan and TXSScan run on finished tapes and are stored as sidecars, not tokens (`PLAN:282-284`; `spec:54`; `DEPLOY:151`). They are excluded from the B1 dictionary (`PLAN:512`).
- **Partial.** Sidecar keys may be `(contig_uid, slot_ordinal)` (`DEPLOY:148`). The key chosen for system calls is not stated in PLAN. critique_4 proposed `(dense_contig_id, slot range)` (`audits/critique_4_bioinformatics.md:132`), which `DEPLOY:132` forbids.
- **Missing.**
  - A command to run a user-supplied MacSyFinder model.
  - CRISPRCasTyper.
  - PADLOC is described as a MacSyFinder call (`PLAN:282`), but PADLOC has its own engine *(recalled, verify)* and is not in the image (`DEPLOY:83`).

---

## 3. Inconsistencies and holes

**A. Missing core capability**
1. **No ad hoc search command or query class.** The CLI list (`DEPLOY:93-105`) and both query-class lists (`PLAN:96-105`; `DL:262-271`) take tokens, UIDs or SQL, never a sequence, HMM, MSA, CM or structure.
2. **No per-release search target databases.**
   - FASTA exists only on demand (`DL:85`; `DEPLOY:193`).
   - MMseqs2 DBs are transient (`SZ:96`).
   - The pack's "Foldseek DBs" are reference DBs (`DEPLOY:71`).
   - Absent everywhere: NR90/NR50 FASTA, pyhmmer digital DB, MMseqs2 search/profile DB, HH-suite DB, corpus Foldseek DB, IGR90 FASTA, contig nucleotide DB.
3. **No MSA or HMM store** for corpus clusters (Stage 2d profiles are internal, `PLAN:308`) or for user profiles.
4. **No way to register user hits as tokens.**
   - Tokens come from frozen dictionaries (`DL:60`).
   - `scan-new-models` writes sidecars (`DEPLOY:164`).
   - Sidecars are joined after filtering (`DEPLOY:150`), so a user family never reaches `slot_pairs`, `pair_stats`, planes or `rescore`.
   - `remap-token` merges existing tokens only (`DEPLOY:155-156`).
5. **Sidecar coverage rule versus sparse hit lists.** "warns below 95% coverage, refuses below 50%" (`DEPLOY:149`; `DL:68`). If coverage means the fraction of the release's keys present, every hit list would be refused. If it means the fraction of sidecar keys that resolve, it passes. The text does not say which.
6. **No six-frame or unannotated-ORF search over contigs**, and no persistent nucleotide DB (§2 Q-NT).
7. **Reverse lookup `protein_uid → slots` is unspecified.** `slot_protein` is defined slot-first (`DL:83`). A hit list is protein-first, and at S6 one protein may map to many slots, given the 35% unique estimate (`SZ:113`).
8. **No envelope-level extraction.** `segment_uid` exists (`DL:57`), but `seq get/extract` have no segment mode (`DL:120-124`).
9. **External closure search is unspecified.** The tool is not named, and PDB seqres and Swiss-Prot are absent from refdb (`spec:71`; `PLAN:452`; `PF:84-108`).

**B. Member-to-representative alignments**

10. PLAN says members keep "their alignment to the representative (`protein_member` table)" (`PLAN:270-271`), and `scan-new-models` projects "through the stored cluster alignments" (`DEPLOY:161`).
11. The L2 layout lists only `nr90/nr50/subfam/igr`, `cluster_member/` and `id_history` (`DL:188`); there is no alignment table.
12. `protein_member` appears only in PLAN (`PLAN:271,376`).
13. The sizing gives membership at 16 B/protein (`SZ:96`), which is too small for alignments (CIGAR or start/end).
14. critique_2's Stage 2d cost assumes a persistent `NR90_to_NR50_alnDB` (`audits/critique_2_systems.md:117`). PLAN says "no loose files" and never names it (`PLAN:308`).

**C. Where the base scan runs, and E-values**

15. **Scan unit, three versions:** per genome in L1 (`DL:36`), NR50 residues (`PLAN:118`), NR90 reps (`SZ:85-86`). This decides whether hits are exact per member or projected.
16. **Z is undefined.** "Fixed, recorded" (`spec:54`) and "fixed per release" (`DEPLOY:162`; `deck:256`), but no value or definition: sequences in NR50? all proteins? the stage?
    - **Conflict with accretion equivalence.** E-values in L1 change with Z. A-02 requires L1 to be byte-identical between incremental and from-scratch builds (`DL:282`). A per-release Z changes stored E-values across stages, so Z must be one constant across the release line, or E-values must be kept out of L1.
    - **Infernal.** The search space ("database size", `DEPLOY:162`) for orphan CMs with their own thresholds (`PLAN:240`) needs the same rule.
17. **E-value comparability for ad hoc searches is undefined.** Library hits are thresholded by GA/TC bitscores, which do not depend on Z (`PLAN:278`). User searches have no GA, so their E-values depend on the target (NR50 vs all proteins vs a stage) and the Z passed. Nothing states the Z for a user search, how iteration inclusion is set, or how a hit count at S5 compares with one at S6 (`SZ:126` only says to "record Z").

**D. Naming and schema drift that a search path would touch**

18. Cluster UID has three definitions (`PLAN:359`; `DEPLOY:130`; `DL:37,59`).
19. `slot_uid` is a `uint64` digest of contig, coordinates and kind (`PLAN:395`) in one place and a composite `(contig_uid, slot_ordinal)` in another (`DL:58`; `DEPLOY:148`). `protein_uids` is `LIST<uint64>` (`PLAN:419`), but `protein_uid` is BLAKE3-128 (`PLAN:349`; `DL:56`).
20. Query-class numbering differs between PLAN and DL: Q3 is anchor profile vs gapped pattern, Q7 is sequence slice vs closure lookup (`PLAN:96-105`; `DL:262-271`). The lab questions are also numbered Q1–Q10 (`PF:67-78`).
21. Store layout and name differ: `l0seq/seq_chunks` + `seq_locator` (`DL:96-101`) vs `sequences/contigs/part-*.parquet` (`DEPLOY:189-194`) vs `seq_index.parquet` (`SZ:48`).
22. The extract flags differ: `--sql` (`DL:124`) vs `--query <view> --kind` (`DEPLOY:104`).
23. The refdb path differs: `ref_models/` (`DEPLOY:187`) vs `l0/refdb/` (`DL:397`).
24. Pair tables are named `window_postings` (`STAGES:171`; `SZ:81`), `pair_postings`/`pair_summary` (`DL:38,192`), and `slot_pairs`/`pair_stats` (`PLAN:423,429`).
25. **Old stage labels in the sequence-store sizing.** `DL:150-157` uses the pre-STAGES numbering: S1 = phages + plasmids + reduced genomes, S3 = GTDB reps. `STAGES:356-365` remaps only DL §4.9, not §2A.4.
26. **Protein residues at S5 disagree.** About 4.8×10⁸ proteins (`SZ:95,104`) × 330 aa (`SZ:41`) = 158 Gaa (derived). `SZ:106` says 145 Gaa. 0.30 × 400–450 Gbp = 120–135 Gaa (derived, `SZ:41`). Protein store size is 45 GB (`DL:155`) vs 0.165 GB/Gbp × 400–450 = 66–74 GB (derived, `SZ:58`).

**E. Tools and data that a search would need but the image or refdb lacks**

27. HH-suite: used in Stage 2d (`PLAN:308`) and tagged in the wrapper (`DL:379`), but not in the image (`DEPLOY:83`).
28. BLAST+/RPS-BLAST: CDD route (`PLAN:280`), `CDD_site` (`PLAN:293`). Not in the image.
29. eggNOG-mapper: COG (`PLAN:279`). Not in the image.
30. PADLOC (`PLAN:282`; `PF:98`), CMfinder, DIAMOND, CRISPRCasTyper: not in the image.
31. No ad hoc installs are allowed (`DL:401`), so each missing tool needs an image release.

**F. Editions**

32. **How `lite` users search is unstated.** `lite` ships the sequence store (`DEPLOY:309`) and runs on a laptop (`DEPLOY:305`), but no search DB. A laptop user must extract about 4.8×10⁸ proteins (about 16 h at the R-08 rate, derived) and build their own DB.
    - IMG/VR proteins are `fetch-required` (`STAGES:346`; `DEPLOY:309`), so searches in `lite` silently miss S4 sequences.
    - It is unclear whether IMG/VR tapes ship. "Contig tapes for those sources" (`DEPLOY:305`) conflicts with "statistics only" (`STAGES:346`).
33. **`core` sequences only on Hetzner.** `core` runs on Hetzner (`DEPLOY:306`), so S6-scale ad hoc search is off-pod by design. No search sizing exists for that host.

**G. User epochs**

34. `append` deltas are "search-only" (`DEPLOY:97`), and private epochs are excluded from `serve` unless enabled (`DEPLOY:270`). It is not stated whether ad hoc searches should include appended or private epochs, or how Z changes when they do.

---

## 4. Sizes per stage (search-target inputs)

Stage counts and Gbp come from `SZ:80-87` and `STAGES:88-95`. Plan constants: aa ≈ 0.30 × bp; 330 aa/protein (`SZ:41`); IGR ≥40 bp ≈ 12% of bp for prokaryotes and 8% for phages (`SZ:42`). S5 includes the S3 genomes (`STAGES:255`); S6 is cumulative (`SZ:86`).

| Stage | Genomes / records | Contig bp | Residues | Proteins | Exact-unique proteins | IGR bp (≥40 bp) |
|---|---|---|---|---|---|---|
| S0 | ~14 + CAMI subset (`SZ:80`) | 0.05 Gbp (+1–3 CAMI) (`SZ:80`) | 15 Maa (+0.3–0.9 Gaa) *derived* | ~4.5×10⁴ (+0.9–2.7×10⁶) *derived* | not given | ~6 Mbp (+0.12–0.36 Gbp) *derived, 12%* |
| S1 | 30–40k phages (`SZ:81`) | 1.7–2.2 Gbp (`SZ:81`) | 0.51–0.66 Gaa *derived* | 1.5–2.0×10⁶ *derived at 330 aa*; 2.2–3.3×10⁶ *estimate at 200–230 aa phage mean* | not given | 0.14–0.18 Gbp *derived, 8%* |
| S2 | 72,360 plasmids (`SZ:82`) | 5–7 Gbp (`SZ:82`) | 1.5–2.1 Gaa *derived* | 4.5–6.4×10⁶ *derived* | not given | 0.6–0.84 Gbp *derived, 12%* |
| S3 | 9–12k reps <2 Mbp (`SZ:83`) | 10–13 Gbp (`SZ:83`) | 3.0–3.9 Gaa *derived* | 0.9–1.2×10⁷ *derived* | not given | 1.2–1.6 Gbp *derived* |
| S4 | 231,408 vOTUs (+699,973 IMG/PR) (`SZ:84`) | ~12 Gbp (+10–30) (`SZ:84`) | 3.6 Gaa (+3–9) *derived* | ~1.1×10⁷ (+0.9–2.7×10⁷) *derived* | not given | ~1.0 Gbp *derived, 8%* (+1.2–3.6 Gbp at 12%) |
| S5 | 143,614 reps (`SZ:85`) | 400–450 Gbp (`SZ:85`) | 145 Gaa (`SZ:106`); 120–135 Gaa *derived*; 158 Gaa *derived from 4.8×10⁸* | **~4.8×10⁸** (`SZ:95,104`); 6×10⁸ slots (`SZ:97`) | not given | 48–54 Gbp *derived* |
| S6 | 732,475 genomes (`SZ:86`) | 2,700–3,100 Gbp cum. (`SZ:86`) | 810–930 Gaa *derived* | **~3×10⁹** (`SZ:113`); Tier 1 7–10×10⁹ pre-dedup (`PLAN:43`) | **~1×10⁹** (35%, `SZ:113`); 0.1–0.2 TB `protein_seq` (`DL:157`) | 324–372 Gbp *derived* |

**NR90 and NR50 counts are not given in any plan file.** critique_2's 300–800M NR90 and 30M dark NR50 (`audits/critique_2_systems.md:19,116`) are at the withdrawn v2 planetary scale (`PLAN:38`). The labeled estimates below are ±2× and must be replaced by S0–S3 measurements.

| Stage | NR90 (estimate) | NR50 (estimate) | Reasoning |
|---|---|---|---|
| S0 | ~0.9× proteins | ~0.8× | 14 phylogenetically distant genomes; little redundancy |
| S1 | 0.3–0.5× → 0.5–1.5×10⁶ | 0.15–0.3× → 0.2–1.0×10⁶ | INPHARED has many near-identical phages per host genus |
| S2 | 0.2–0.4× → 0.9–2.6×10⁶ | 0.1–0.25× | PLSDB is dominated by near-identical Enterobacteriaceae plasmids |
| S3 | 0.7–0.9× → 0.6–1.1×10⁷ | 0.4–0.6× | Distinct small genomes in divergent lineages (DPANN, CPR) |
| S5 | 0.5–0.7× → 2.4–3.4×10⁸ | 0.2–0.35× → 1.0–1.7×10⁸ | Species reps are ≥5% ANI apart, but conserved orthologs across congeners exceed 90% AA identity. UniRef50/UniRef90 ≈ 0.3 and AFDB 214M → ~52M at 50% are the outside anchors *(recalled, verify)* |
| S6 | 3.5–6×10⁸ (S5 × 1.3–2) | 1.2–2.5×10⁸ (S5 × 1.2–1.5) | The 588,861 added genomes are strains of S5 species. They add mostly <90% variants of accessory genes and few new 50% families. |

**Derived target-DB footprints** (sizes are estimates; the formulas are shown):
- pyhmmer digital DB ≈ 1 B/residue plus per-sequence overhead *(estimate)*:
  - S5 all proteins about 135–160 GB;
  - S5 NR50 about 33–56 GB (1.0–1.7×10⁸ × 330 aa);
  - S6 exact-unique about 330 GB.
- MMseqs2 DB ≈ 2× aa bytes (`SZ:96`):
  - S5 about 300 GB (`SZ:96`);
  - S6 exact-unique about 0.66 TB;
  - S6 all about 1.6–1.9 TB.
- Contig nucleotide target: 400–450 Gbp at S5 (`SZ:85`), 2.7–3.1 Tbp at S6 (`SZ:86`).
- Per-model `hmmsearch` cost from the plan's Pfam rate (≈73 core-h/Gaa across about 22k Pfam models, `SZ:29`): about 3.3×10⁻³ core-h/Gaa per model. That gives about 0.45–0.5 core-h per HMM per pass at S5 and about 1.1 core-h at S6 exact-unique (*derived*). The anchor is a Pfam-vs-proteins benchmark, so single-model cost is unmeasured and scales with model length (`SZ:29,126`).

---

## 5. wilkcomp02 envelope constraints for search

| Constraint | Value | Evidence | Effect on an interactive search |
|---|---|---|---|
| Memory cap | `systemd-run --scope MemoryMax=200G` per job; <20% of RAM by day; DuckDB `memory_limit` ≤150 GB | `PLAN:76`; `DEPLOY:27,40`; `DL:16,318`; `SZ:13` | An in-RAM pyhmmer DB of all S5 proteins (about 135–160 GB, §4) nearly fills a scope. S6 cannot be held in RAM and must stream. NR50 at S5 (about 33–56 GB) fits. |
| Address space | `ulimit -v` 500 GB soft; "never mmap a target DB larger than ~400 GB in one process" | `DEPLOY:32`; `DL:17` | The S6 MMseqs2 DB (0.66–1.9 TB) must be split across processes |
| `/tmp` | ~940 GB, volatile, about 800 GB usable per job; S5 build peak 500–900 GB | `PLAN:77`; `DEPLOY:29`; `SZ:14,106`; `PF:30` | A persistent search DB on `/tmp` competes with the S5 build and can vanish on reboot. NFS-resident DBs are slower to mmap (K-10, `DL:326`). |
| NFS small files | No per-genome **or per-query** files on NFS; only ≥64 MB files; run-step refuses >1,000 files | `DEPLOY:41`; `DL:19-20,384`; `SZ:17` | Per-query outputs (HMMs, MSAs, hit tables, per-hit alignments) must be packed into one Parquet or tar per query, or kept on `/tmp` |
| Watchdog | Renices any PID past ~369 CPU-min and mails rctf-support. Monolithic tools start at nice 19; fan-out at nice 15 and <300 CPU-min per PID. | `PLAN:79`; `DEPLOY:39`; `SZ:18`; `STAGES:322` | One 40-thread process reaches 369 CPU-min in about 9 min wall (derived: 369/40). A 5-iteration jackhmmer at S5 is about 2.3–2.5 core-h ≈ 140–150 CPU-min (derived, §4) and stays under the limit. HHblits, MMseqs2 iterative or Foldseek over S6 would not, so they must start at nice 19. |
| Threads | 40 by day, 56 at night, physical cores only, explicit `--threads` and BLAS env vars | `PLAN:75`; `DEPLOY:38`; `DL:15,380` | The build plans for 800 core-h/day, about 70% duty (`SZ:11`), over 25–100 days at S5 (`SZ:108`). Interactive searches compete for the same 40 threads, and no reservation for them is modeled. |
| No scheduler | Local executor; run-step applies the caps itself | `PLAN:80`; `DEPLOY:31`; `DL:23` | There is no queue for interactive jobs. Run-step's "wait if other users' load > 6 cores" (`DL:381`) and the STOP file (`DL:382`) can block or kill an interactive search. Run-step is designed for batch work. |
| Apptainer only | `.sif` only; no ad hoc installs | `DEPLOY:33`; `DL:394,401` | Missing tools (HH-suite, BLAST+, CMfinder; §3 E) need an image release before anyone can use them |
| No GPU | ProstT5 and Foldseek-GPU only on wilkcomp03 by arrangement, or on a GEX130 | `DEPLOY:28`; `spec:31`; `PF:25` | Q-STRUCT queries without a structure need ProstT5. CPU use for a single query is unaddressed (`DEPLOY:87`). |
| `serve` | localhost:8420 via an SSH tunnel; parameterized endpoints; timeouts, row limits, concurrency caps; SQL console without external access | `DEPLOY:57,268-273`; `STAGES:78` | No search endpoint exists. Minute-scale searches would exceed any query timeout, so search has to be an async job and not a `/v1` call. With `enable_external_access=false`, the SQL console cannot read a user's hit file. |
| I/O | One NFS server, 10 GbE, shared by three hosts; `ionice -c 2 -n 7` | `DEPLOY:30,39`; `SZ:17` | Streaming about 66–80 GB of compressed S5 proteins takes at least 1 min per pass at line rate (derived). Extracting them at the R-08 target (`DL:172`) takes about 16 h. |
| Open files | 4,096 hard (8,192 via `ulimit -n`) | `DEPLOY:32`; `DL:22` | Caps sharded target DBs and parallel readers |
