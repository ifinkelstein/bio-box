# Audit Report 8: Data Quality, Provenance, Lineage and Licensing

**Auditor:** Data-quality and provenance lead (large public biological data integration)
**Target documents:** `00-spec.md` (v2.0), `PLAN.md` (v2.0), `DEPLOYMENT-AND-MALLEABILITY.md`, `PRIOR-ART.md`, `audits/*`
**Verdict:** The plan has no data-quality layer. It specifies how bytes move and how statistics are computed, but not what enters the system, how much of it is contaminated or duplicated, how any row traces back to a source record, how counts reconcile across stages, or under what terms the result can be redistributed. Every downstream statistic in `STAT-01/02` (`n_independent_orders`, `N_eff^veh`, `P(B|A)`, `q_emp`) is a count over inputs that are, as specified, double-counted across sources, unscreened for contamination and chimeras, and taxonomically unlabelled for most of Tier 3. The fixes below are mostly cheap relative to the build; they are a precondition for the benchmarks meaning anything.

---

## Step A: Deepening the task

**(1) What a top professional would do that the request does not ask for.**
Treat the plan as a data product and audit it against data contracts, not prose. That means writing down, per source, the exact release, object inventory, checksum, license and redistribution terms; computing the expected overlap matrix between sources before any code is written; defining a conservation ledger so that every input record ends in exactly one of {kept, merged-into-X, rejected-for-reason-R}; and building a golden regression set whose correct answers come from curated external resources. They would also check the plan's numbers against the sources' own published statistics, because a plan that misstates its inputs will mis-size everything.

**(2) The rawest material to start from.**
Not the plan's source table, but: the sources' own release manifests and stats pages (GTDB `RELEASE_NOTES.txt` and `stats/r226`, Logan GitHub release notes and the `logan-pub` bucket listing, AllTheBacteria release metadata, IMG/VR v4 paper Table 1, GlobDB paper Table 1); NCBI Assembly and SRA run metadata (BioSample/BioProject/run accessions, `library_strategy`, `library_source`, organism); license and policy pages (JGI data policy, ISfinder About page, AWS Open Data registry entries, GSA-Human policy); and published contamination/chimerism base rates (FCS-GX, GUNC).

**(3) Methods and tools.**
Accession-level identity resolution (GCA/GCF pairing, BioSample → run mapping via SRA metadata in BigQuery/Athena), sequence digests (SHA-512/24 or GA4GH refget/seqcol digests) for exact duplicates, containment sketches (sourmash/skani/sylph) for near-duplicates; NCBI FCS-GX and FCS-adaptor, UniVec, phiX and spike-in screens; CheckM2 + GUNC for genome quality and chimerism; geNomad for replicon class; contig domain classifiers (Tiara/Whokaryote, or FCS-GX taxonomy) for eukaryotic removal; GTDB-Tk/skani at a single pinned GTDB release; MMseqs2 `taxonomy` with LCA confidence for Tier 3 contigs; pandera or Great Expectations or dbt-duckdb tests for schema and invariants; OpenLineage-style run metadata; stratified random audit samples with binomial confidence intervals; CAMI II gold-standard metagenomes for end-to-end Tier 3 truth.

**(4) What, other than us, shows the result is right.**
Curated truth sets that the pipeline never saw: EcoCyc/RegulonDB transcription units for *E. coli* K-12, SubtiWiki/DBTBS for *B. subtilis* 168, RefSeq curated annotations of model phages, Rfam seed alignments, Swiss-Prot experimental entries, Ribo-seq-supported smORFs, CAMI II assemblies with known genome-of-origin, and a prospective set of systems published after the data freeze. Also: the sources' own published counts (reconciliation must match them), and an independent rebuild of a 1% sample by a second operator that produces byte-identical outputs.

**Rewritten prompt.**
> Audit `all-the-data` as a data product. For every source in `PLAN.md` §2.2, write a data contract (release, retrieval method, object inventory and checksum, record count from the source's own stats, license and redistribution terms, citation), and verify the plan's numbers against the source. Build the expected cross-source overlap matrix and specify identity resolution (accession, sequence digest, containment sketch) and source precedence. Specify contamination, chimera and artifact screens with tools, thresholds and known base rates. Specify record-level lineage so any slot resolves to source accession, coordinates, caller, mode and tool versions. Specify a per-stage rejection ledger with a conservation invariant and per-source reconciliation. Specify schema and invariant tests. Define a golden-genome regression set with external truth sources and pass metrics, plus identifier stability, release versioning and drift monitoring. Replace absolute claims with measured rates and confidence intervals. Rank findings by how much they corrupt the coupling statistics, and give each a concrete fix and a validation that does not depend on the pipeline itself.

---

## Step B: Findings

Severity scale: **CRITICAL** (corrupts primary statistics or blocks lawful release), **MAJOR** (biases results or breaks traceability), **MODERATE** (correctness or operability debt).

### DQ-01 Cross-source duplication is unhandled; the same genome enters 3 to 6 times (CRITICAL)
- **Where:** `PLAN.md` §2.2 table, §2.3; `00-spec.md` §1.
- **Failure:** Tier 1 lists GTDB (596,859 genomes), RefSeq Complete, GlobDB, proGenomes3 and curated assemblies as "~1.8M assemblies". These are not disjoint: GlobDB is built on GTDB species representatives plus SPIRE, mOTUs and others ([GlobDB paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12629238/)); proGenomes and RefSeq Complete are drawn from NCBI Assembly, which is also GTDB's source ([GTDB R226 release notes](https://data.gtdb.aau.ecogenomic.org/releases/latest/RELEASE_NOTES.txt)). Tier 3 then adds AllTheBacteria, a uniform reassembly of SRA isolate reads, many of whose GenBank assemblies are already in GTDB ([AllTheBacteria preprint](https://www.biorxiv.org/content/10.1101/2024.03.08.584059.full.pdf)), and Logan, which assembles nearly all public SRA runs, including those same isolates and the metagenomes from which SPIRE, IMG/M and IMG/VR were derived ([Logan](https://github.com/IndexThePlanet/Logan)). PLSDB plasmids reappear inside RefSeq complete genomes; IMG/VR viruses are mined from IMG/M metagenomes that are also in Logan.
  A strain sequenced once can therefore appear as (GTDB GCA) + (RefSeq GCF) + (proGenomes) + (AllTheBacteria reassembly) + (Logan contigs of the same run) + (PLSDB for its plasmid). The Stage 0 sieve deduplicates only Tier 3 by architecture hash and only within a biome bucket; it does not know that two records are the same biological sample. `n_total_occurrences`, `arch_occurrence_count`, `pan_any_bits` support counts and every `P(B|A)` denominator inherit this multiplicity non-uniformly (well-studied pathogens are inflated most).
- **Fix:** Add **Stage −1 Identity Resolution** before the sieve:
  1. Map every record to `(biosample, run_set, assembly_accession_pair)` using NCBI Assembly (GCA↔GCF pairs) and SRA run metadata. One BioSample = one biological sample regardless of how many assemblies exist.
  2. Exact duplicates by canonical sequence digest (e.g., refget-style SHA-512/24 of uppercase, strand-canonical sequence).
  3. Near-duplicates by containment, not Jaccard (`skani`/`sylph`/sourmash containment ≥ 0.99 at ≥ 99.5% ANI), because Logan fragments are contained in finished assemblies.
  4. Source precedence for the retained representative: RefSeq Complete > GenBank complete > AllTheBacteria > GTDB draft > GlobDB/SPIRE MAG > Logan contig. Every other copy becomes a `member_of` row, never a new occurrence.
  5. Publish the overlap matrix (sources × sources, by BioSample and by bases) in each release.
- **Validation:** Count of unique BioSamples per tier must be reproducible from NCBI metadata alone. Spot check 200 random Tier 3 Logan contigs from isolate runs: each must resolve to its Tier 1 assembly when one exists.

### DQ-02 No contamination, artifact or chimera screening anywhere in the pipeline (CRITICAL)
- **Where:** `PLAN.md` §2.3 (Filters 1–4), §3.1; `00-spec.md` invariants 1–2.
- **Failure:** The only non-prokaryote screen is a 64 MB eukaryotic repeat Bloom filter. Known base rates make this inadequate:
  - NCBI found contaminants in about one in ten prokaryote genome submissions in 2022 ([FCS-GX, Genome Biology 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC10246020)).
  - GUNC estimated 15–30% of CheckM-prefiltered "high-quality" MAGs and about 5% of GenBank/RefSeq genomes are undetected chimeras ([GUNC, Genome Biology 2021](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02393-0)).
  - Logan assembles all public SRA, which includes Illumina phiX spike-ins, vectors and adapters, host DNA (human, mouse, plant), and reagent contaminants. A phiX174 contig will be assembled in a very large number of runs and will dominate Tier 2-like phage statistics unless removed.
  - Non-repeat eukaryotic contigs ≥ 1,200 bp pass the repeat filter, and `pyrodigal --meta` will call ORFs in them, producing fake "prokaryotic" families and adjacencies.
  - Chimeric metagenomic contigs create adjacencies that never existed in any genome. Adjacency is the primary signal of the engine, so chimeras are not noise: they are false positives in the exact quantity being ranked (`junction_topology_prev`, Rosetta-stone fusions at assembly junctions).
- **Fix:** Insert a **Stage 0b Screen** before gene calling:
  1. FCS-adaptor (adapter/vector, UniVec) and an exact phiX/Illumina spike-in/lambda-control k-mer screen; trim or drop.
  2. Domain-of-life classification per contig (FCS-GX with its curated database, or Tiara/Whokaryote plus sketch hits to GTDB vs. eukaryotic references). Keep eukaryotic contigs only if a prokaryotic/viral region is explicitly excised and flagged (`host_associated_symbiont` rescue must be evidence-based, not SD-motif based).
  3. Assembly-level CheckM2 completeness/contamination and GUNC CSS for every MAG and SAG; store per assembly, propagate to contigs.
  4. Contig-level chimera evidence for Tier 3: coverage discontinuity (Logan provides k-mer abundance; `count_logan_tig_coverage` exists in Bioconda), GUNC-style per-contig lineage homogeneity, and a rule that a novel adjacency supported by a single contig in a single BioSample cannot contribute to `N_eff^veh`.
  5. Library-strategy filter from SRA metadata: exclude `AMPLICON`, `ChIP-Seq`, `Bisulfite-Seq`, targeted capture; flag `RNA-Seq`/metatranscriptome as `molecule_type = RNA`.
- **Validation:** Spike-in experiment: inject known phiX, UniVec, human and *Saccharomyces* fragments and 1,000 synthetic chimeras (two genomes joined at random intergenic sites) into a 10,000-contig test set; report sensitivity and specificity with 95% CIs. Re-measure on CAMI II assemblies, where genome-of-origin per contig is known.

### DQ-03 Stage 0 dedup logic loses representatives and members, and its ANI threshold is miscalibrated (CRITICAL)
- **Where:** `PLAN.md` §2.3 Filter 4; `critique_2_systems.md` §4.2.
- **Failure:** Three defects:
  1. On first sighting the contig is written to `tier3_unique_contigs.bin.zst` and its fingerprint to the Cuckoo filter, but nothing is written to RocksDB. On second sighting the flow looks up `(stage0_arch_hash, biome_bucket)` in RocksDB, which has no entry. The plan does not say what happens then, so either every second sighting becomes a new representative (no dedup) or it is dropped (data loss). The Cuckoo filter is keyed on hash only, while RocksDB is keyed on hash plus biome, so the two disagree by construction.
  2. `POPCNT(sketch_a XOR sketch_b) ≤ 6` on a 64-bit one-permutation b-bit sketch does not mean 95% ANI. With 1-bit values, the expected mismatch fraction is (1 − J)/2; 6/64 implies J ≈ 0.81, which for k = 15 corresponds to about 99.3% ANI (Mash formula). 95% ANI at k = 15 gives J ≈ 0.30 and about 22 mismatched bits. With 64 bits the binomial variance is large, so the threshold also produces many random false merges and misses. And Jaccard on contigs of different lengths understates identity for fragments contained in longer contigs.
  3. Collapsed contigs keep only `(rep_contig_id, sample_id, biome_id)`. Their own coordinates, sequence and gene calls are discarded, so a user can never see the member that actually carried a polymorphic CTD or a CRISPR spacer, and lineage from a reported count back to source records is broken (see DQ-05).
- **Fix:** Write every first sighting to RocksDB (the Cuckoo filter is only a negative-lookup accelerator). Key both structures identically. Replace the 64-bit sketch test with a containment estimate from a FracMinHash sketch (scaled 1/1000) or a direct banded alignment check on candidate pairs, and calibrate the threshold empirically against skani ANI on 10,000 random candidate pairs. Keep a member ledger: `contig_member(member_source_accession, member_digest, rep_dense_contig_id, ani_est, containment, rule_id)`.
- **Validation:** Calibration curve of the merge test against skani ANI with reported false-merge and missed-merge rates. Conservation check: members + representatives + rejected = input contigs for every Logan object.

### DQ-04 Tier assignment conflates source, replicon type and quality (MAJOR)
- **Where:** `PLAN.md` §2.2; `00-spec.md` invariant 5; `contig_tapes.tier`.
- **Failure:** Tier is a single `uint8` per contig but means three different things. A plasmid contig inside a RefSeq complete genome is "Tier 1" by source and "Tier 2" by replicon. AllTheBacteria isolates (uniform, QC'd) sit in Tier 3 next to unbinned Logan contigs and single-cell genomes (high chimera and contamination rates). SAGs are listed as Tier 3 with no quality gate. Statistics stratified by tier then mix high-quality isolates with the noisiest data.
- **Fix:** Replace `tier` with orthogonal columns: `source_class` (isolate_complete, isolate_draft, MAG, SAG, unbinned_contig, virus_db, plasmid_db), `replicon_class` (chromosome, plasmid, virus, provirus, unknown; from geNomad plus source annotation), `assembly_quality` (CheckM2 completeness/contamination, GUNC CSS, N50, contig count) and `molecule_type` (DNA/RNA). Keep `tier` as a derived view for compatibility.
- **Validation:** Cross-tab of `source_class` × `replicon_class` published per release; every complete-genome plasmid must have `replicon_class = plasmid`.

### DQ-05 No record-level lineage from a slot to its source (MAJOR)
- **Where:** `PLAN.md` §4.2 (`contig_tapes` schema), §3.2; `00-spec.md` invariants 2–3.
- **Failure:** A slot carries `nr90_protein_ids`, so its annotation comes from the NR90 representative, not the slot's own protein. Envelope boundaries are projected from NR50 via CIGAR, and `catalytic_motif_hash` and `catalytic_intact` are computed on the representative. At 90% (or 50%) identity a member can differ at precisely the active-site columns, so the Level 4 token, the closure score `C_iso` and Benchmark 6B inherit the representative's residues. Nothing records which source accession, contig name and coordinates a slot came from, which gene-caller mode produced it, which genetic code was chosen, or which Stage 1.5 rule merged it.
- **Fix:** Add provenance tables, not wider rows:
  - `contig_provenance(dense_contig_id, source_db, source_release, source_accession, source_contig_name, source_coords, seq_digest, biosample, run_accessions[], ingest_run_id)`.
  - `slot_provenance(dense_contig_id, slot_idx, caller ∈ {pyrodigal, smorf_mode1..5, infernal, minced, ...}, genetic_code, edit_rule_id, evidence_code, tool_version_id)`, with evidence codes mapped to the Evidence and Conclusion Ontology where possible.
  - `protein_member(member_digest, nr90_id, nr50_id, subfamily_50_id, aln_identity, aln_cigar)` and compute `catalytic_motif_hash` per unique member sequence (cheap: it is a column lookup on an existing alignment).
  - A `tool_versions` table (container digest, tool version, database release, parameters hash) referenced by `tool_version_id`.
- **Validation:** For 1,000 random slots, an automated resolver must return the source accession and coordinates, refetch the source sequence, re-translate, and match the stored protein digest. Target 100% resolvable; report the rate.

### DQ-06 No rejection ledger or stage-by-stage reconciliation (MAJOR)
- **Where:** `PLAN.md` §1.3, §2.3, §3.
- **Failure:** Filters drop "78% of short unitig fragments" and reject contigs with "< 2 features", but no counts are persisted by source and reason. S3 transfers at 400 MB/s for 68 h will hit transient errors; without an object inventory and per-object completion records, a partially ingested Logan release is indistinguishable from a complete one. Note also that 40.8 TB at 400 MB/s is about 28 h, not 68 h, so the plan's own timing does not reconcile.
- **Fix:** A `ledger` table written by every stage: `(stage, source_db, source_object, n_in, n_kept, n_merged, n_rejected, reason_code, bytes_in, bytes_out, checksum_in, run_id)`. Enforce the invariant `n_in = n_kept + n_merged + Σ n_rejected` per object and stage; a release is blocked if any object is missing or unbalanced. Use the S3 inventory (or a pinned listing) as the expected-object manifest.
- **Validation:** Per-source totals reconcile to the source's own published counts (GTDB stats page, AllTheBacteria release metadata, IMG/VR release tables) within a documented tolerance.

### DQ-07 Tier 3 taxonomy is unspecified, yet dispersion statistics depend on it (CRITICAL)
- **Where:** `PLAN.md` §4.2 (`tax_order_id` in sort key and schema), §4.3 (`n_independent_orders`, decoys "within order"); `PRIOR-ART.md` breakthrough 4.
- **Failure:** `tax_order_id` is used for sorting, dispersion counts, decoy stratification and benchmarks, but nothing says how an unbinned 3 kb Logan contig gets an order. Contig-level classifiers resolve a minority of short contigs to order, misassign HGT islands to the donor lineage, and are biased toward well-represented clades. Misassignment inflates `n_independent_orders` for exactly the mobile loci the stats are trying to discount. Unassigned contigs either get a sentinel (and become one giant fake "order") or are dropped (biasing toward known taxa).
- **Fix:** Assign taxonomy with one pinned method and release (e.g., MMseqs2 `taxonomy` LCA against GTDB at the same release as the tree, or GTDB-Tk/skani for bins). Store `tax_rank_resolved`, `tax_confidence`, `tax_method_id`. Compute dispersion only over contigs resolved at the required rank with confidence ≥ threshold, and report how many were excluded. Never use a sentinel order as a category in counts. Use a single GTDB release across all sources: AllTheBacteria, GlobDB and the tree currently use different releases.
- **Validation:** On CAMI II and on Tier 1 genomes fragmented to Tier 3 lengths, measure order-level precision and recall of contig taxonomy; run the coupling statistics with and without low-confidence contigs and report rank stability of the top 1,000 pairs.

### DQ-08 Licensing and redistribution are unexamined; the planetary edition may not be legally publishable (CRITICAL)
- **Where:** `DEPLOYMENT-AND-MALLEABILITY.md` §4 (public S3/R2 editions), `PLAN.md` §6.3.
- **Failure:** The plan streams derived data packs to anyone. Verified terms:
  - **ISfinder:** its About page states the database may not be downloaded without written authorization ([ISfinder About](https://tncentral.ncc.unesp.br/ISfinder/about.php)). Bulk ingestion and redistribution of IS sequences or IS-derived HMMs needs permission.
  - **IMG/VR, IMG/PR, IMG/M:** governed by JGI data policy. Legacy policy for proposals accepted after Nov 2018 lifts use restrictions on publication or two years after public release, asks users to contact the PI before publishing and to acknowledge JGI ([JGI legacy data policy](https://jgi.doe.gov/data-policy-support/legacy-data-policy), [JGI data policy](https://jgi.doe.gov/data-policy-support/data-policy)). Some IMG datasets may still be within restriction windows at freeze time.
  - **Logan:** the AWS registry lists no license; it points to NCBI policy and the NIH Genomic Data Sharing Policy ([AWS registry](https://registry.opendata.aws/pasteur-logan)). Logan contigs from human-derived samples contain human sequence; redistributing it in derived packs is a privacy and policy issue regardless of license.
  - **AllTheBacteria:** listed as MIT on the AWS registry ([AWS registry](https://registry.opendata.aws/allthebacteria)).
  - **CNCB GSA-Human:** controlled access by default, with DAC approval and Chinese human genetic resource regulations ([GSA-Human policy](https://ngdc.cncb.ac.cn/gsa-human/policy)). Do not ingest from GSA-Human; for non-human GSA/CNGBdb data, record per-dataset terms.
  - **ICEberg 3.0, TnCentral, PaperBLAST, Fitness Browser:** no license statement found ([ICEberg Database Commons](https://ngdc.cncb.ac.cn/databasecommons/database/id/513), [TnCentral About](https://tncentral.ncc.unesp.br/about), [PaperBLAST GitHub](https://github.com/morgannprice/PaperBLAST)). Unknown terms mean redistribution of copies is not cleared; derived counts are lower risk than sequences or text.
  - **Benefit sharing:** sequence data derived from biodiversity may fall under the CBD digital sequence information (DSI) mechanism agreed at COP16 (2024). Check before any commercial use.
- **Fix:** A `sources.yaml` data contract per source with `license`, `redistribution ∈ {allowed, attribution, derived-only, permission-required, prohibited, unknown}`, `citation`, `contact`. The edition builder includes a source's sequences only if redistribution is `allowed`/`attribution`; `derived-only` sources contribute counts and identifiers but no sequence; `permission-required`/`unknown` are excluded from public packs until cleared. Ship `LICENSES/` and `CITATION.cff` inside every pack. Human-read removal (DQ-02) is mandatory before any public pack.
- **Validation:** A release check that fails if any file in a public pack derives from a source whose redistribution flag is not cleared (lineage from DQ-05 makes this computable).

### DQ-09 The plan's source facts are stale or wrong (MAJOR)
- **Where:** `PLAN.md` §2.2; `00-spec.md` §1; `PRIOR-ART.md`.
- **Failure and corrections (verified):**
  | Plan says | Source says | Consequence |
  |---|---|---|
  | GTDB r220, 113,104 species, 596,859 genomes | Current is R10-RS226 (16 Apr 2025): 143,614 species clusters, 732,475 genomes ([GTDB stats r226](https://gtdb.ecogenomic.org/stats/r226)) | Bitvector width `14,144 B` is hardcoded to r220; must be a function of the pinned release (r226 needs 143,614 bits → 17,984 B at 64-B padding) |
  | GlobDB r220 "~145,000" | GlobDB r226: 306,260 species representatives ([GlobDB](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12629238/)) | Tier 1 size and overlap with GTDB misestimated |
  | IMG/VR v4: 15.7M contigs / 5.5M vOTUs | >15M genomes/fragments, 8.7M vOTUs ([IMG/VR v4, NAR 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC9825611)) | Vehicle dedup by vOTU sized wrong |
  | Logan 27.3M runs, 40.8 TB | Logan v1 covered a Dec 2023 freeze; v1.1 recomputed contigs to fix spurious repetition; v1.2 covers a Dec 2025 freeze ([Logan](https://github.com/IndexThePlanet/Logan)) | Must pin a version ≥ v1.1; v1 contigs carry a known artifact |
  | proGenomes3 (~1M) | proGenomes4 reports ~2M genomes ([PubMed 41261732](https://pubmed.ncbi.nlm.nih.gov/41261732/)) | Version not pinned |
  | AllTheBacteria ~1.9M | v0.1 had 1,932,812 assemblies; AWS registry describes an update covering SRA isolates to Aug 2024 ([preprint](https://www.biorxiv.org/content/10.1101/2024.03.08.584059.full.pdf), [registry](https://registry.opendata.aws/allthebacteria)) | Version not pinned |
  | GORG-Tropics + GORG-Dark, ">50,000 SAGs" | GORG-Tropics ~12,000 SAGs ([Bigelow](https://www.bigelow.org/news/articles/2019-12-12.html)); GORG-Dark ~9,000 in a 2025 abstract only ([ASLO 2025](https://aslo.secure-platform.com/2025/solicitations/16/sessiongallery/1928/application/8828)) | Count unsupported; GORG-Dark availability unverified |
- **Fix:** Every number in §2.2 becomes a field in `sources.yaml` populated from the source's manifest at freeze, with the URL and retrieval date. Prose references the contract, not a literal. Parameterize all release-dependent sizes (bitvector width, species count, tree edges).
- **Validation:** CI recomputes each count from the downloaded manifest and fails on drift beyond tolerance.

### DQ-10 Internal count inconsistencies imply most "contigs" have one slot (MAJOR)
- **Where:** `PLAN.md` §1.2, §2.2, §4.1; `critique_2_systems.md` §2.
- **Failure:** The plan states 2.5 × 10⁹ post-sieve contigs and 3.6 × 10⁹ slots (1.44 slots per contig), while the sieve rejects contigs with < 2 features unless the longest ORF ≥ 250 aa, and the mission claims "multi-feature" loci. critique_2 uses 5 × 10⁸ contigs. The Tier rows sum slots (480M + 320M + 2.8B = 3.6B), but Tier 1 alone at ~1.8M assemblies × ~4,000 genes is ~7 × 10⁹ ORFs before dedup. These cannot all be true, and the index sizes, window counts and `uint32` headroom depend on which is.
- **Fix:** Derive every scale number from one sizing model (inputs per source × measured survival rate per gate from a 0.1% pilot), and publish it as a table with units. Re-derive index sizes from it.
- **Validation:** Pilot run on a stratified 0.1% sample of each source; extrapolated counts with CIs replace hand estimates.

### DQ-11 Annotation comparability: E-values, thresholds and model provenance are lost (MAJOR)
- **Where:** `PLAN.md` §3.2 Stage 2b ("115k non-redundant HMMs"), `DEPLOYMENT-AND-MALLEABILITY.md` §3.4 (`scan-new-models`).
- **Failure:** Collapsing Pfam, COG, CDD, NCBIfam and others into "non-redundant" profiles discards each library's own gathering/trusted cutoffs and model accession. HMMER E-values scale with the search-space size Z; an incremental 50-model scan reports E-values on a different Z than the base build, so thresholds like `E < 1e-5` in Stage 1.4 are not comparable across runs. Reference library versions (Pfam, Rfam 15.0, COG, PHROGs, DefenseFinder, PADLOC, ECOD, AFDB, ESMAtlas) are not pinned anywhere in the data model.
- **Fix:** Keep every hit with `(model_db, model_accession, model_version, bitscore, ga_cutoff_passed)`; choose representatives for dedup in a separate mapping table rather than deleting models. Fix and record Z (`--Z` / `Z=` in pyhmmer) for all scans, including incremental ones. Rank by bitscore relative to model cutoffs, not raw E-values. Pin each reference DB in `tool_versions`.
- **Validation:** Re-scan a 1M-sequence sample with the incremental path and the full path; hit sets must match exactly.

### DQ-12 Transferred 3Di strings carry no structure-confidence filter (MAJOR)
- **Where:** `PLAN.md` §3.2 Stage 2a; `00-spec.md` §1 item 2.
- **Failure:** 3Di states are copied from AFDB/ESMAtlas/BFVD predictions via alignment at ≥ 50% identity. Many ESMAtlas and AFDB models have low-confidence regions; 3Di over low-pLDDT residues is noise, and fold assignments built on it propagate to `struct_fold_id`, `G_leak`, and gLM training labels.
- **Fix:** Transfer per-residue pLDDT with the 3Di string; mask residues with pLDDT < 70; require ≥ 70% of the envelope above threshold to assign `struct_fold_id`; store `source_structure_id`, `source_db_version`, `mean_plddt`.
- **Validation:** On 1,000 domains with experimental PDB structures, compare ECOD assignment from transferred 3Di vs. from the PDB structure, stratified by pLDDT.

### DQ-13 Identifiers are not stable across releases, and sidecars/remaps silently invalidate stats (MAJOR)
- **Where:** `PLAN.md` §4.1 (`dense_contig_id` ordered by `(tier, tax_order_id, stage2_arch_hash)`); `DEPLOYMENT-AND-MALLEABILITY.md` §3.1–3.3.
- **Failure:** `dense_contig_id` depends on taxonomy and architecture hash, so any taxonomy update, new source or re-clustering renumbers every contig. `subfamily_50_id` and `seq_family_id` for dark clusters are cluster ordinals that change with every reclustering. Sidecars keyed on these IDs then join to the wrong entities after an upgrade, with no error. `remap-token` changes token identity in 200 ms but leaves precomputed `func_coupling_llr`, entropies, bitvectors, `transition_edges` and q-values computed under the old mapping. `append` epochs add counts additively, but entropy, single-linkage vehicle clusters and q-values are not additive, and `append` does not check whether the new genome is already in the base (double counting; see DQ-01).
- **Fix:**
  - Stable external IDs: content digests for sequences and contigs; cluster IDs with a release-to-release mapping table (`id_history(old_id, new_id, release, relation ∈ {same, split, merged, retired})`, as UniRef and GTDB do for clusters).
  - Every sidecar declares `key_type`, `key_release`, `producer`, `created_at`; `register-sidecar` validates key uniqueness, coverage (fraction of keys resolved, warn below 95%) and release compatibility.
  - Derived-stat tables carry `computed_at_release`, `token_map_version`; `remap-token` and `append` mark affected rows `stale = true` and queries surface staleness. Non-additive stats are recomputed for affected tokens or flagged.
  - `append` runs Stage −1 identity resolution against the base before allocating IDs.
- **Validation:** Upgrade test: register a sidecar on release N, upgrade to N+1, and verify joins through `id_history` reproduce the same biological entities for a sampled 1,000 keys.

### DQ-14 Schema has overflow, enum and parallel-array risks with no data tests (MODERATE)
- **Where:** `PLAN.md` §3.3 (badge bits), §4.2 schema.
- **Failure:** `delta_bp_prev_slots` is `int16` (max 32,767 bp); gaps larger than that occur on long contigs with uncalled regions and will wrap silently. `length_aa_bin` is 12 bits (4,095 aa); NRPS/PKS megasynthases and collapsed repeat super-segments exceed it. `regulon_locus_ids` is `uint16` in PLAN but `uint32` in critique_1. `feature_kinds` has no codes for rRNA, intron, homing endonuclease, pseudogene, ribozyme or Obelisk, and critique_2 uses a different 4-bit encoding. Parallel `LIST<T>` columns can drift in length with no check.
- **Fix:** Saturating encodings with an overflow flag bit; an enum registry table (`enum_name, code, label, so_term, introduced_release`) mapped to Sequence Ontology terms; one canonical schema file (Arrow schema + JSON Schema) from which docs and code are generated. Data tests (pandera on Polars or dbt-duckdb tests) run on every partition:
  - all `LIST` columns per row have length `n_slots`;
  - `0 ≤ start < end ≤ contig_length`; strands ∈ {−1, +1};
  - `orfs_to_contig_start + orfs_to_contig_end + 1 = n_orfs` per ORF;
  - `directon_ids`, `regulon_locus_ids` non-decreasing; no window spans two contigs;
  - every foreign key resolves (including across editions; see DQ-15);
  - enum values ∈ registry; counts non-negative; probabilities ∈ [0, 1].
- **Validation:** Tests are blocking in the build; the release report lists test counts and failures.

### DQ-15 The lite edition ships summaries whose evidence is absent (MODERATE)
- **Where:** `DEPLOYMENT-AND-MALLEABILITY.md` §4.1.
- **Failure:** `lite` includes the 3-tier `directon_pair_summary` and `rep_dense_contig_id` pointers, but only Tier 1 species-rep tapes. Most evidence pointers dangle, so users cannot inspect why a pair scored high. A summary that cannot be audited is not reviewable science.
- **Fix:** Every summary row carries `evidence_edition_min` and up to N evidence pointers that resolve in `lite` (preferring Tier 1 representatives); the UI states when evidence requires `core`/`planetary`.
- **Validation:** FK test per edition: the fraction of resolvable evidence pointers is reported and ≥ 1 resolvable pointer exists for every row shown by default.

### DQ-16 Benchmarks leak ground truth through the annotation libraries, and absolute pass criteria are unmeasurable (MAJOR)
- **Where:** `PLAN.md` §5; `audits/SYNTHESIS.md` item 25; `00-spec.md` invariant 9.
- **Failure:** Benchmark 1 asks the engine to "rediscover" CBASS, Thoeris, Gabija, etc., while Stage 2b annotates with DefenseFinder, PADLOC and Pfam families that were defined from those discoveries. Recall then measures whether the HMM library contains the answer. Criteria like "0.0% false hole calls", "100% (200/200)", "zero HGT leakage guaranteed" and "every multi-feature genetic locus on Earth" are not measurable claims: 0/200 failures still leaves a 95% upper bound of about 1.5% (rule of three), and coverage of "Earth" is bounded by what was sequenced and deposited.
- **Fix:** Leave-system-out annotation: for each benchmark system, rebuild the relevant slots with the system's defining models removed (and models dated after the discovery), then test recall. Add a prospective set: systems published after the data freeze date. State every criterion as an estimated rate with a 95% CI and a minimum sample size. Replace "every locus on Earth" with "public INSDC/SRA-derived assemblies as of freeze date D, after filters F, with coverage reported per source".
- **Validation:** Benchmarks are scored by a script against a frozen, versioned truth file maintained separately from the pipeline.

### DQ-17 Metadata quality (biome, host, sample) is assumed (MODERATE)
- **Where:** `PLAN.md` §2.3 (`biome_bucket`), §4.2 (`genome_or_sample_id`, biome, host).
- **Failure:** SRA BioSample attributes are free text and inconsistently populated; biome and host labels drive dedup buckets and stratification. Tier 2 phage/plasmid host assignments (needed for `pan_any_bits` per species) are predictions (CRISPR spacer matches, iPHoP-like tools) with substantial error.
- **Fix:** Harmonize biome to a controlled vocabulary (ENVO or the MGnify biome ontology) with `biome_source ∈ {curated, inferred, missing}`. Store host predictions with method and confidence; include MGE contigs in `pan_any_bits` only above a host-confidence threshold, and report sensitivity.
- **Validation:** Manual audit of 300 random BioSamples for biome label accuracy; host-assignment precision measured on phages and plasmids with isolation hosts.

### DQ-18 No release versioning, drift monitoring or reproducibility check (MODERATE)
- **Where:** `DEPLOYMENT-AND-MALLEABILITY.md` §4.2 ("quarterly release"), `PLAN.md` §6.
- **Failure:** `manifest.blake3.json` gives integrity, not meaning. Nothing records which source releases, tool versions and parameters produced a pack, or how results changed between releases.
- **Fix:** A `release_manifest.json` (semver; source contract versions; container digest; tool_versions; parameter hashes; build date; data freeze date). A drift report per release: per-source and per-tier counts, family frequency shifts, top-1,000 LLR pair Jaccard and Spearman rank correlation versus the previous release, benchmark deltas, schema diff. Reproducibility: rebuild a 1% stratified sample twice from scratch with the pinned container and diff outputs (expect byte-identical, or document nondeterministic steps).
- **Validation:** Release blocked if drift exceeds documented thresholds without an explanatory changelog entry.

### DQ-19 smORF and Stage 1.5 calls have no evidence tiers and will dominate counts with false positives (MAJOR)
- **Where:** `PLAN.md` §3.1 items 2, 5, 6.
- **Failure:** A 5-mode smORF caller across 52 TB will produce a very large number of short ORFs; most short ORFs in random sequence are not genes. Frameshift and intron stitching merges are rule-based inferences. All enter `contig_tapes` as equal slots and shift `orf_idx`, `directon_id` and `P(B|A)`.
- **Fix:** Per-slot `evidence_level` (e.g., 1 = profile hit above GA, 2 = conserved across ≥ 3 independent BioSamples with dN/dS < 0.5, 3 = motif-only). Default statistics use levels 1–2; level 3 is queryable but excluded from coupling stats unless requested. Same for Stage 1.5 merges, with the rule ID recorded.
- **Validation:** Precision of smORF calls against Ribo-seq-supported smORFs in *E. coli* K-12 and *Haloferax volcanii*; precision of merges against the golden set below.

---

## Data quality gates per stage

| Stage | Gate | Check | Tool / method | Threshold / action | Ledger reason code | Blocking |
|---|---|---|---|---|---|---|
| −2 Acquisition | Contract | Source release pinned, license and redistribution recorded | `sources.yaml` | Missing field → stop | `NO_CONTRACT` | Yes |
| −2 Acquisition | Completeness | Downloaded objects = manifest/S3 inventory; checksums match | `s5cmd` + inventory diff, BLAKE3/MD5/ETag | Any missing or mismatched → retry, then stop | `MISSING_OBJECT`, `CHECKSUM_FAIL` | Yes |
| −2 Acquisition | Count reconciliation | Record counts vs. source stats page | Script | > 0.5% deviation → investigate | `COUNT_DRIFT` | Yes |
| −1 Identity | Sample identity | BioSample/run/assembly mapping | NCBI Assembly + SRA metadata | Same BioSample → one entity | `SAME_BIOSAMPLE` | No |
| −1 Identity | Exact duplicate | Canonical sequence digest | SHA-512/24 | Duplicate → member row | `EXACT_DUP` | No |
| −1 Identity | Near duplicate | Containment and ANI | skani / sylph / sourmash | Containment ≥ 0.99 and ANI ≥ 99.5% → member | `CONTAINED_DUP` | No |
| −1 Identity | Library type | SRA `library_strategy`/`library_source` | Metadata filter | Amplicon, ChIP, bisulfite, capture → drop; RNA → flag | `EXCLUDED_LIBRARY` | No |
| 0b Screen | Adapter/vector | UniVec, adapters | FCS-adaptor | Trim or drop | `ADAPTOR_VECTOR` | No |
| 0b Screen | Spike-ins | phiX174, lambda control, ERCC | Exact k-mer screen | Drop | `SPIKE_IN` | No |
| 0b Screen | Foreign / host | Domain of life, human | FCS-GX (or Tiara + sketch hits) | Eukaryotic contig → drop unless excised prokaryotic region; human → drop always | `EUK_HOST`, `HUMAN` | No |
| 0b Screen | Genome quality | Completeness, contamination, chimerism per assembly | CheckM2, GUNC | Record; MAG/SAG with contamination > 5% or GUNC CSS > 0.45 → excluded from pangenome planes | `LOW_QUALITY_ASM`, `CHIMERIC_ASM` | No |
| 0b Screen | Contig chimera | Coverage break, lineage heterogeneity | Logan k-mer abundance, per-contig LCA | Flag `chimera_suspect`; single-contig adjacencies excluded from `N_eff^veh` | `CHIMERA_SUSPECT` | No |
| 0 Sieve | Length/complexity | ≥ 1,200 bp, DUST | Rust filter | As specified; record counts | `SHORT`, `LOW_COMPLEXITY` | No |
| 0 Sieve | Dedup | Containment-based merge, member ledger | FracMinHash + calibrated threshold | Calibration report required | `ARCH_DUP` | Yes (calibration) |
| 1 Calling | Gene-call sanity | Coding density, stop-codon usage, code choice | pyrodigal stats | Density < 50% or > 98% on prokaryotic contig → flag | `ODD_CODING_DENSITY` | No |
| 1 Calling | Evidence level | smORF and merge evidence tiers | Rule engine | Level 3 excluded from default stats | — | No |
| 1 Calling | Golden set | Regression genomes (table below) | Test harness | Any metric below target → stop | — | Yes |
| 2a Cluster | Member lineage | Every protein maps to NR90, NR50, subfamily | Join test | 100% mapped | `UNMAPPED_PROTEIN` | Yes |
| 2b Annotate | Model provenance | Hit carries db, accession, version, GA pass, fixed Z | Schema test | Missing → fail | — | Yes |
| 2b Annotate | Active site per member | `catalytic_motif_hash` from member, not rep | Test on sample | Mismatch rate reported | — | No |
| 2c Structure | Confidence | pLDDT carried; low-confidence masked | Schema + filter | `struct_fold_id` only if ≥ 70% envelope pLDDT ≥ 70 | `LOW_PLDDT` | No |
| 2 Taxonomy | Tier 3 taxonomy | Rank, confidence, method recorded; single GTDB release | MMseqs2 taxonomy / GTDB-Tk / skani | Below confidence → excluded from dispersion counts | `TAX_UNRESOLVED` | No |
| 3 Build | Schema/invariants | Data tests in DQ-14 | pandera / dbt-duckdb | Any failure → stop | — | Yes |
| 3 Build | Conservation | `n_in = kept + merged + rejected` per stage and source | Ledger check | Any imbalance → stop | — | Yes |
| 3 Build | FK integrity | All pointers resolve, per edition | SQL tests | Report resolvable fraction | — | Yes (core/planetary) |
| Release | Licensing | Every packed file's sources cleared for redistribution | Lineage + `sources.yaml` | Uncleared → excluded | `LICENSE_BLOCK` | Yes |
| Release | Drift | Counts, top-pair stability, benchmarks vs. previous release | Drift report | Exceeds thresholds without changelog → stop | — | Yes |
| Release | Reproducibility | 1% rebuild diff | Pinned container | Non-identical outputs documented | — | Yes |
| Append / sidecar | Identity and key checks | Dedup against base; key coverage and release | `append`, `register-sidecar` | Coverage < 95% → warn; stale stats flagged | `APPEND_DUP`, `SIDECAR_KEY_MISS` | No |

---

## Golden-genome regression set

Every build runs these end to end through Stage 3. Truth comes from curated resources outside the pipeline. Accessions should be pinned in `golden.yaml` with sequence digests; the ones marked "pin" need an accession chosen at freeze.

| # | Genome | Accession | What it tests | External truth | Pass metric (target) |
|---|---|---|---|---|---|
| 1 | *E. coli* K-12 MG1655 | NC_000913.3 | Operons, divergons, `directon_id`, prfB +1 frameshift, fdhF/fdoG/fdnG selenoproteins, rhs loci, leader peptides (trpL, hisL, pheL) | EcoCyc/RegulonDB transcription units; RefSeq annotation | Transcription-unit adjacency recall ≥ 0.95; prfB and 3 selenoproteins merged; leader peptides present as `uORF_cis_leader` |
| 2 | *B. subtilis* 168 | NC_000964.3 | Riboswitches and T-boxes with attenuator terminators (no `directon_id` break), ykkC guanidine-I, SAM-I, glmS | SubtiWiki, DBTBS, Rfam annotations | ≥ 95% of known riboswitches detected and linked to downstream gene in same `directon_id` |
| 3 | Phage T4 | NC_000866.4 | Group I introns in td, nrdB, nrdD with homing endonucleases | RefSeq annotation, literature | All three host genes stitched (`INTRON_SPLICED`), HEGs parented |
| 4 | Phage lambda | NC_001416.1 | gpG/gpGT −1 programmed frameshift | RefSeq annotation | Merged as one logical ORF with `PROGRAMMED_FRAMESHIFT` |
| 5 | *Mycoplasma genitalium* G37 | NC_000908.2 | Genetic code 4 selection | RefSeq (code 4) | Code 4 chosen; no UGA-truncated genes vs. reference |
| 6 | SR1/Gracilibacteria closed genome (e.g., SR1 RAAC1, Campbell et al. 2013) | pin | Genetic code 25 (UGA → Gly) | Published annotation | Code 25 chosen; gene count within 5% of reference |
| 7 | Lak megaphage or other amber-recoded phage (Devoto et al. 2019) | pin | Code 15 (UAG reassignment) | Published annotation | Code 15 chosen on affected contigs |
| 8 | *Methanosarcina acetivorans* C2A | NC_003552.1 | Pyrrolysine (mtmB/mtbB/mttB) readthrough | RefSeq, literature | Pyl genes merged with `SEC_OR_PYL_READTHROUGH` |
| 9 | *Haloferax volcanii* DS2 | NC_013967.1 | Leaderless transcripts, archaeal smORFs | Published TSS and Ribo-seq maps | smORF precision ≥ 0.8 at Ribo-seq-supported set; leaderless starts recovered |
| 10 | *Sulfolobus acidocaldarius* DSM 639 | NC_007181.1 | Archaeal leaderless initiation, TATA/BRE mode | Published TSS maps | Start-site agreement ≥ 0.9 |
| 11 | *Mycobacterium tuberculosis* H37Rv | NC_000962.3 | Actinobacterial leaderless genes | Published TSS maps (e.g., Cortes 2013; Shell 2015) | Leaderless starts recovered ≥ 0.85 |
| 12 | *Bacteroides thetaiotaomicron* VPI-5482 | NC_004663.1 | SD-independent AT-rich initiation (mode 4) | RefSeq; Ribo-seq where available | Start-site agreement ≥ 0.9; no SD-mode bias |
| 13 | *Streptococcus pyogenes* SF370 | NC_002737.2 | Cas9 nested RuvC/HNH (discontinuous domain), CRISPR array arbitration | UniProt Q99ZW2 domain map; CRISPRCasdb | One stitched RuvC with HNH parented; `catalytic_intact = true`; no shadow ORFs in array |
| 14 | Obelisk reference (Zheludev et al. 2024) | pin | Circular RNA rescue branch | Published sequence | Retained with Oblin-1 call |
| 15 | CAMI II marine and strain-madness assemblies | CAMI II | Tier 3 chimera, taxonomy, dedup behavior | CAMI gold standard genome-of-origin | Chimera detection sensitivity/specificity, order-level taxonomy precision, dedup false-merge rate reported with CIs |
| 16 | Spike-in contamination panel (synthetic) | generated | phiX, UniVec, human, yeast fragments; synthetic chimeras | Known by construction | Removal sensitivity ≥ 0.99 for spike-ins; chimera flag sensitivity reported |

---

## Hyperbolic or unmeasurable claims to rewrite

| Location | Claim | Replace with |
|---|---|---|
| `00-spec.md` §1 | "indexes every multi-feature genetic locus on Earth" | "indexes multi-feature loci from public INSDC/SRA-derived assemblies frozen on date D, after filters F; per-source coverage in the release report" |
| `00-spec.md` §1 | "entire sequenced ... biosphere" | "the public sequenced fraction available under the source contracts in `sources.yaml`" |
| `PLAN.md` §4.5 | "guaranteeing 0.0% HGT operon leakage and 0.0% remote-homology token leakage" | "leakage measured by cross-split homology search; report rate and 95% CI" |
| `PLAN.md` §5 B3 | "0.0% false hole calls" | "false hole rate ≤ x% (95% upper bound) on n simulated holes" |
| `PLAN.md` §5 B6 | "100% (200/200)" | "failure rate with 95% upper bound (0/200 → ≤ 1.5%)" |
| `PRIOR-ART.md` exec summary | "No." (nothing like this exists) | A scoped claim naming the specific combination not offered elsewhere, with dates of the survey |

---

## Out of scope, noted for other auditors
- `PackedSlotDescriptor` stores `orf_mod_8`, which cannot verify ORF distances up to ±6 without ambiguity across windows (architecture).
- The 16 GB Cuckoo filter for 10⁹ 16-bit fingerprints is about 8× the needed ~2 GB (systems).
- Decoy stratification "within order" inherits DQ-07 taxonomy error (statistics).
