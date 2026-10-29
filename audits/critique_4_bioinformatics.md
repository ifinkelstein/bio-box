# Audit Report 4: Bioinformatics and Annotation Practice Review of `all-the-data` v2.0

Auditor: annotation practitioner lens (someone who has run pyrodigal, Infernal, MMseqs2, Foldseek, HMMER, geNomad, CheckM2, GTDB-Tk on millions of genomes).
Targets: `00-spec.md`, `PLAN.md` v2.0, `DEPLOYMENT-AND-MALLEABILITY.md`, `PRIOR-ART.md`, `audits/*`.
Verdict: the v2.0 plan fixes real v1.0 problems, but it fixed them by adding new rules on paper. Several of those rules contradict each other, rely on tool behavior that doesn't exist, or quote input sizes that are out of date by a factor of 5 to 10. Nothing in the plan says how any annotation layer will be measured against ground truth before the 3.45 TB build. Read the numeric guarantees as targets, not facts.

---

## Step A: deepening the task

1. What a top professional would do that wasn't asked for. Before writing a pipeline, they would run a pilot of about 1,000 genomes plus about 10,000 Logan runs end to end with off-the-shelf tools. They would measure each annotation layer's precision and recall against curated truth sets, and only then decide which bespoke rules (5-mode smORFs, CP stitching, attenuator exemption) are worth building. They would also reconcile every input size against the current release before budgeting compute.
2. Rawest material. Use the actual inputs, not the plan's tables or the earlier critiques' digest of them: the Logan S3 listing plus SRA run metadata (BioSample, organism, library strategy, NCBI STAT k-mer taxonomy profile per run), GTDB release files (`bac120_metadata`, trees), and raw Pfam/Rfam/NCBIfam model files with their GA/TC cutoffs and nesting annotations. Plus a few hundred genomes with experimental annotation (RegulonDB, term-seq, Ribo-seq, PDB-mapped domains).
3. Methods and tools. Gene calling: pyrodigal single mode per genome, meta mode or pyrodigal-gv for short or viral contigs, Codetta for genetic code. RNA: Infernal `cmsearch` with its native filter pipeline and `cmcalibrate`. Domains: pyhmmer with GA/TC cutoffs, cath-resolve-hits, Pfam nesting data. Structure: Foldseek with ProstT5 `createdb`. MGEs and viruses: geNomad, CheckV, IntegronFinder 2, MacSyFinder (CONJScan, DefenseFinder, TXSScan), ISEScan. Genome quality and taxonomy: CheckM2, GUNC, GTDB-Tk, skani or sylph, MMseqs2 `taxonomy`, CAT/BAT. Hosts: iPHoP and CRISPR spacer matching. Truth sets: CAMI II.
4. What shows it is right, other than us. External truth sets and benchmark data:
   - RegulonDB operons and transcription units
   - term-seq terminator maps (Dar et al. 2016)
   - Ribo-seq and Ribo-RET smORF sets (Weaver 2019; Sberro 2019)
   - ECOD/CATH manual domain partitions for discontinuous domains
   - M-CSA catalytic residues
   - Rfam seed alignments
   - geNomad/CheckV benchmark sets
   - CAMI II simulated metagenomes, for Tier 3 taxonomy
   - PGAP annotations of RefSeq reference genomes
   - Fitness Browser cofitness, held out
   - concordance of results across two independent GTDB releases

Rewritten prompt:
> Audit the `all-the-data` plan as a microbial genome annotation practitioner. (a) Reconcile every input dataset size, version, and format against the current primary release, and recompute the downstream counts. (b) For each annotation layer (genes and smORFs, genetic code, ncRNA, terminators and directons, domains including nested, CP, and repeat cases, active sites, 3Di and fold, MGE regime, taxonomy and host, genome quality), name the standard tool a practitioner would use. State whether the plan's bespoke rule beats it, contradicts it, or relies on behavior the tool doesn't have. (c) Find internal contradictions across files. (d) For every layer, name the external truth set and metric that must pass in a pilot (about 1,000 genomes plus about 10,000 Logan runs) before the planetary build. Cite primary sources; mark unverified numbers.

---

## Step B: findings

Severity scale: CRITICAL (wrong science or an unbuildable stage), MAJOR (large measurable error), MODERATE.

### BX-01 [CRITICAL] Input volumes and versions are stale. The Logan figure is roughly 10× too small.
- Files: `00-spec.md` §2 Inv. 2; `PLAN.md` §1.3, §2.2, §2.3; `DEPLOYMENT` §4.1.
- Facts:
  - Logan contigs live at `s3://logan-pub/c/<acc>/<acc>.contigs.fa.zst`. The project describes their total as "hundreds of terabytes compressed, tens of millions of files" ([Logan Contigs.md](https://github.com/IndexThePlanet/Logan/blob/main/Contigs.md)). The v1.2 release covers the December 2025 SRA freeze, about 38M accessions and 87 Pbp ([Logan README](https://github.com/IndexThePlanet/Logan)). The plan uses 40.8 TB and 27.3M runs, which describes the v1 2023 snapshot.
  - GTDB R226 (April 2025) has 143,614 species representatives among 732,475 genomes ([GTDB r226 stats](https://gtdb.ecogenomic.org/stats/r226)). The plan hard-codes r220: 113,104 species and 596,859 genomes.
  - GlobDB r226 has 306,260 species representatives ([GlobDB paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12629238/)). The plan says about 145,000.
  - AllTheBacteria's current release has 2,440,377 genomes ([preprint](https://www.biorxiv.org/content/10.1101/2024.03.08.584059.full.pdf)). The plan says 1.9M.
  - proGenomes4 (about 2M genomes) has superseded proGenomes3 ([PubMed 41261732](https://pubmed.ncbi.nlm.nih.gov/41261732/)).
  - IMG/VR v4 has about 15.7M UViGs in 8.7M vOTUs, not 5.5M vOTUs ([Camargo 2023 NAR](https://pmc.ncbi.nlm.nih.gov/articles/PMC9825611)).
- Failure: Stage 0's "68 h at 400 MB/s" assumes 40.8 TB. At several hundred TB, the same rate gives about 3 to 4 weeks of streaming alone. Every downstream count (2.5B contigs, 8.5B ORFs, 85M NR50) inherits the error. The fixed 14,144-byte bitvector rows are wrong for R226: 143,614 bits padded to 64 B gives 17,984 B.
- Fix: add `data_sources.lock.yaml` pinning release, URL, byte size, and record count per source, generated by a script from S3 listings and release files. Derive bitvector width from `n_species` at build time; never hard-code it. Re-derive the §1.3 budget from measured per-run statistics (see BX-02).
- Validate: the script's totals must match the provider's published totals within 1%.

### BX-02 [CRITICAL] Most of Logan is not microbial. Filtering should start from SRA metadata, not a 64 MB Bloom filter.
- Files: `PLAN.md` §2.3 Filter 2.
- Failure: SRA is dominated by human, mouse, plant, and model-organism data. Filter 2 (a "64-MB Euk-Repeat Syncmer Bloom Filter" plus a Shine-Dalgarno whitelist) is applied after downloading and decompressing every byte. An SD motif occurs by chance roughly every few hundred bp, so the whitelist will pass much eukaryotic sequence. Meanwhile real symbionts in host runs are a tiny fraction of the bytes.
- Fix:
  - Pre-select runs using SRA metadata: library strategy, source, organism, BioSample `env_*` fields. Use the NCBI SRA STAT k-mer taxonomy profile per run, available in the SRA BigQuery tables (verify field names), and keep runs with a prokaryotic or viral fraction above a threshold.
  - For host-dominated runs that pass, classify contigs with a real classifier: Kraken2 against a GTDB-plus-host index, Tiara or Whokaryote for eukaryotic contigs, and geNomad for viral and plasmid contigs.
  - Record per-run `stat_prok_frac` so tiers and biomes are traceable.
- Validate: on 1,000 random Logan runs, report the bytes kept and the recall of known symbiont contigs (Buchnera, Wolbachia) against a full scan.

### BX-03 [CRITICAL] Taxonomy for Tier 2 and Tier 3 is undefined, yet `tax_order_id` drives sorting, decoys, bitmasks, and partitions.
- Files: `PLAN.md` §4.1 (sort key), §4.2 (`window_postings` sort), §4.3 (order bitmask), STAT-02 decoys; `critique_3` §1.2.
- Failure: unbinned Minia contigs (median 2 to 6 ORFs per `critique_3`) have no taxonomy. Phage and plasmid contigs have no GTDB order. Without a defined assignment, `tax_order_id` defaults silently, and `n_independent_orders` and within-order decoys become meaningless for about 80% of slots.
- Fix:
  - Bacterial and archaeal contigs: assign `tax_lineage` per contig with MMseqs2 `taxonomy` (2bLCA) against GTDB R226 proteins, or CAT. Store `tax_rank_resolved` and `tax_confidence`.
  - Viral contigs: geNomad taxonomy (ICTV) plus a host prediction (iPHoP or spacer matches) stored as `host_tax_order_id` with `host_confidence`.
  - Plasmids: host from PLSDB or IMG/PR metadata, or CRISPR and k-mer matching.
  - Add explicit `UNASSIGNED` buckets. Count orders only from confidently assigned slots, and report the unassigned fraction next to every dispersion count.
- Validate: CAMI II marine and strain-madness datasets for order-level precision and recall at contig lengths of 1 to 10 kb.

### BX-04 [CRITICAL] IGR90 clustering erases the riboswitch subtypes that BIO-11 promises to resolve.
- Files: `PLAN.md` §3.1 item 3a vs 3d; `00-spec.md` §1.3; Benchmark 4.
- Failure: ykkC subtypes (Guanidine-I, ppGpp, PRPP) and purine G/A/2'-dG differ at 1 to 5 pocket nucleotides. 90% identity on a 100 nt aptamer tolerates about 10 differences. Members of one IGR90 cluster can therefore carry different ligand specificities. The pocket subtyper sees only the representative, so the subtype label is projected wrongly onto every member.
- Fix: use IGR90 only to choose what to scan. After a CM hit on a representative, re-align every member IGR in that cluster to the hit CM with `cmalign` (cheap, since only hit clusters are realigned) and apply `cm_pocket_rules.json` per member. Store the subtype per slot, not per cluster.
- Validate: ykkC, purine, and GEMM subtype seed sets from the Breaker lab papers (Sherlock 2017–2019) must have per-member subtype accuracy of at least 98%.

### BX-05 [MAJOR] An nhmmer prefilter in front of Infernal duplicates Infernal's own filter and drops structure-dominant families. Orphan CMs have no GA thresholds.
- Files: `PLAN.md` §3.1 item 3b/3c; `critique_2` §1.1.
- Failures:
  1. `cmsearch` already runs a profile-HMM filter cascade (MSV, Viterbi, Forward, then HMM-banded CYK). An external nhmmer pass at p ≤ 1e-3 is a second, uncalibrated HMM filter. It loses families whose signal is mainly covariation with weak primary conservation (some T-box regions, long ncRNAs like OLE/GOLLD/HEARO). The quoted speed-up is unverified.
  2. `--cut_ga` only works for models with GA lines. The roughly 550 "orphan" Weinberg/Breaker models compiled from `.sto` have none. Many Weinberg motifs are already in Rfam 14 and 15, so counting both double-counts them.
  3. "1,450 non-eukaryotic Rfam 15.0 models" is unverified. Rfam 15.1 lists 4,227 families ([rfam.org](https://rfam.org/)), and the prokaryotic subset must be derived from the Rfam taxonomy tables, not assumed.
- Fix:
  - Drop the external nhmmer stage. Run `cmsearch --cpu N --tblout` with default or `--rfam` filters on IGR90.
  - Build models with `cmbuild` and calibrate with `cmcalibrate`. Use E ≤ 1e-5 thresholds for non-Rfam models, set per model from scans of shuffled IGRs.
  - Deduplicate orphan models against Rfam by seed overlap.
  - Derive the prokaryotic or viral model list from `Rfam.full_region` taxonomy.
- Validate: on 500 complete genomes, compare recall of full `cmsearch` against the prefiltered pipeline, per family. Accept the shortcut only if per-family loss is under 1%.

### BX-06 [MAJOR] Gene calling uses one mode for everything. Genetic-code assignments are partly wrong.
- Files: `PLAN.md` §3.1 item 1; §2.3 Filter 3; `00-spec.md` Inv. 1.
- Facts:
  - Table 25 (UGA→Gly) was described in SR1 (Campbell 2013 PNAS). A Codetta survey on GTDB r220 places it in the Patescibacteria order BD1-5. Absconditibacterales use a distinct code ([summary with refs](https://en.wikipedia.org/wiki/Candidate_division_SR1_and_gracilibacteria_code); Shulgina & Eddy 2021 eLife, Codetta). Most Patescibacteria use code 11.
  - Code 15 (TAG→Gln) is a phage code, common in Crassvirales (Borges 2022 Nat Microbiol; Peters 2022 Nat Commun).
  - The plan writes "15 (UAG→Gln in phages / Absconditabacteria)" and "25 (UGA→Gly in Patescibacteria / CPR)". Both host attributions are wrong or too broad. The claim of code 4 in Hodarchaeales is unverified.
  - Prodigal meta mode can't pick code 15 or 25. pyrodigal-gv adds code-15 phage models and giant-virus models ([pyrodigal-gv](https://github.com/althonos/pyrodigal-gv)), and geNomad reports code 11, 4, or 15 per sequence.
- Failure: running `--meta` on Tier 1 genomes discards per-genome training, which is the most accurate Prodigal mode. Running coding-density comparisons on short Tier 3 contigs is noise.
- Fix:
  - Tier 1: pyrodigal single mode per genome. Pick the code with Codetta, or by coding-density comparison over codes 11, 4, 15, and 25 at genome level.
  - Tier 2 and viral contigs: pyrodigal-gv meta mode, taking geNomad's genetic code call.
  - Tier 3: meta mode; inherit the code from the contig's taxonomic assignment (BX-03) when it is confident.
  - Store `genetic_code` and `gene_caller_mode` per contig.
- Validate: PGAP CDS concordance on 500 RefSeq genomes (start and stop exact-match rates). For code 15, use the Peters 2022 metaproteomics-supported crAss-like phage genomes.

### BX-07 [MAJOR] The 5-mode smORF caller has no false-discovery control. Its dN/dS and RNAcode rescue isn't budgeted.
- Files: `PLAN.md` §3.1 item 2; Stage 1.4 Rule B.
- Failure: random intergenic ORFs of 15 to 75 aa occur at several per kb. Modes 3 and 4 (any ATG within 0 to 6 bp of a TATA/-10-like hexamer; any AT-rich window) will call large numbers of spurious peptides, especially in AT-rich genomes, where those motifs are common by chance. Mode 5 needs multi-contig alignments for RNAcode or dN/dS, which no stage builds or budgets. Size ranges also disagree across files: 15–70 aa (spec), 15–75 aa (plan), 12–60 aa (Rule B), 20–75 aa (critique_1).
- Fix:
  - Keep Mode 1 (pyrodigal with a lowered minimum length) and Mode 5 (profile hits) as primary.
  - Treat Modes 2 to 4 as candidates. Promote a candidate to a slot only if its cluster (MMseqs2 at 50% AA over smORFs) appears in at least 3 independent orders with RNAcode or dN/dS < 0.5 support. This is the Sberro 2019 approach. Budget the clustering and RNAcode stage explicitly.
  - Store `smorf_evidence` as a bitmask.
  - Pick one size range.
- Validate: Ribo-seq and Ribo-RET smORF sets for E. coli, B. subtilis, H. volcanii, and Bacteroides (Weaver 2019 mBio; Sberro 2019 Cell). Report precision against shuffled-intergenic decoys, targeting FDR < 10% for promoted smORFs.

### BX-08 [MAJOR] Split-gene merging can't tell programmed frameshifts from assembly indels or pseudogenes.
- Files: `PLAN.md` §3.1 item 6 (Stage 1.5).
- Failure: ONT and short-read metagenome assemblies, and Minia contigs, carry indel errors that produce exactly the "adjacent same-strand ORFs hitting consecutive halves of one HMM" signature. Merging them inflates `catalytic_intact` and erases real pseudogenes. That corrupts NOGD holes and the Benchmark 6B pseudoenzyme logic.
- Fix:
  - Merge as a programmed frameshift only when a known slippery motif is present (prfB `CUU UGA C`; IS `A AAA AAC`). Requires `ribosomal_slippage`-style motif scans.
  - Or merge when the same frameshift position is conserved in at least 3 independent genomes.
  - Otherwise set `biological_edit_type = CANDIDATE_INDEL` and keep both fragments. Use Selenoprofiles or bSECISearch-style SECIS detection for Sec.
- Validate: prfB, IS3 family, and selenoprotein sets in RefSeq complete genomes; simulated indels at ONT error rates on 100 genomes, targeting false merges under 1%.

### BX-09 [MAJOR] The domain database mix assumes things the tools don't provide.
- Files: `PLAN.md` §3.2 Stage 2b; `critique_2` §1.3.
- Facts:
  - CDD is distributed as RPS-BLAST PSSMs, not HMMs, so pyhmmer can't run it.
  - COG 2020 is distributed via CDD PSSMs. HMMs must be built from COG alignments or taken from eggNOG.
  - DefenseFinder, PADLOC, TXSScan, and CONJScan are system-level MacSyFinder models with gene-order rules. Running their HMMs alone doesn't call systems.
  - "E < 1e-5" cutoffs ignore curated GA (Pfam) and TC (NCBIfam) thresholds.
- Failure: Stage 2b as written can't execute for CDD and COG. Defense annotations would be per-gene hits, not systems. Thresholds would be inconsistent across databases.
- Fix:
  - Run Pfam-A, NCBIfam, PHROGs, and arCOG via pyhmmer with `--cut_ga` or `--cut_tc` where the database provides them.
  - Get COG via eggNOG-mapper or COG-built HMMs. Use CDD via MMseqs2 profile search or RPS-BLAST only for the site annotations needed by `catalytic_motif_hash`.
  - Run DefenseFinder, PADLOC, CONJScan, and TXSScan as MacSyFinder on the tapes after Stage 3. Store system calls as a sidecar keyed by `(dense_contig_id, slot range)`.
- Validate: on DefenseFinder's published RefSeq set, reproduce system calls (Tesson 2022) with concordance of at least 98%.

### BX-10 [MAJOR] The nested, CP, and repeat logic overlaps existing Pfam and ECOD curation. Thresholds disagree.
- Files: `PLAN.md` §3.2 Stage 2b (BIO-02/03/04); `critique_1` BIO-01/03/04.
- Failures:
  - Pfam already encodes allowed nesting (the `NE` field in `Pfam-A.seed`/`.hmm.dat`) and has separate families for many permuted MTases. cath-resolve-hits already accepts discontinuous hits, and ECOD and CATH publish manual discontinuous partitions.
  - The plan's own thresholds disagree: nesting gap of 35 aa or more, repeat-merge gap of 15 aa or less (plan) versus 35 aa (critique_1), and CP HMM overlap under 15%.
  - The CP stitcher will also fire on genuine tandem duplications of half-domains. Without a check, it would also merge those.
- Fix:
  - Use Pfam `NE` nesting rules plus ECOD domain boundaries as the primary source. Apply the gap heuristics only for pairs with no curated rule, and flag those `heuristic=true`.
  - CP merge only when the combined alignment covers at least 70% of the HMM, both halves score above GA, and no curated permuted family scores higher.
  - Set one repeat gap threshold.
- Validate: ECOD manual partitions of PDB chains (discontinuous and nested), CPDB circular-permutation pairs, and Cas9/AspRS exemplars. Report boundary F1 within ±20 residues.

### BX-11 [MAJOR] The 3Di, ProstT5, Foldseek, and AlphaFill claims don't match tool capabilities.
- Files: `PLAN.md` §3.2 Stage 2a/2c; `00-spec.md` §1.2; `critique_2` §1.2.
- Facts:
  - `foldseek databases` distributes Alphafold/UniProt, UniProt50, Proteome, Swiss-Prot, ESMAtlas30 (the about 36M high-quality clustered subset), and PDB ([Foldseek README](https://github.com/steineggerlab/foldseek)). There's no ready-made 3Di for "ESMAtlas-Full 772M". Getting it means downloading roughly 15 TB of structures for the v0 set and running Foldseek `createdb` on them, filtered by pLDDT. Low-confidence models give unreliable 3Di.
  - ProstT5-derived databases contain predicted 3Di only and "do not support TM-score/LDDT" output (Foldseek README). So the plan's `--alignment-type 1` (TM-align) on ProstT5-only segments isn't possible. The claim that this mode "supports circular permutations" is unverified.
  - AlphaFill transplants ligands into AFDB models only. Dark NR50 segments with only ProstT5 3Di have no coordinates to transplant into.
  - Transferring 3Di states through a sequence alignment CIGAR at ≥50% identity is an approximation the plan doesn't validate. The "~72% coverage" and "62,500 aa/s" figures are unverified.
- Fix:
  - Use Foldseek's own ProstT5 path (`createdb --prostt5-model … --gpu 1`) for all NR50 segments without an AFDB or ESMAtlas30 hit. Measure throughput in the pilot.
  - Restrict alignment transfer to hits at ≥90% identity and ≥90% coverage to models with pLDDT ≥ 70.
  - Assign ECOD by 3Di+AA E-value only. Mark `fold_evidence ∈ {STRUCTURE, TRANSFER, PROSTT5}`.
  - Compute AlphaFill pockets only where an AFDB model exists, and state that limit.
- Validate: on 10,000 AFDB proteins with known ECOD, compare ECOD T-group agreement for ProstT5-only, transfer, and real 3Di.

### BX-12 [MAJOR] `catalytic_motif_hash` is brittle and undefined for most families.
- Files: `PLAN.md` §3.2, §4.2; `critique_3` §3.2.
- Failure: exact `xxhash32` over active-site residues splits D/E or S/T conservative variants, and a single sequencing error flips the hash. M-CSA and CDD site annotations cover a small fraction of families. For dark families, "top-8 information-content columns" are conservation columns, not active sites, yet the plan treats a mismatch there as pseudoenzyme evidence.
- Fix: store the `catalytic_motif_str` aligned to a family reference, with `site_source ∈ {MCSA, CDD_SITE, ALPHAFILL, CONSERVATION}`. Keep the hash only as an index key over a reduced alphabet (for example, D=E, S=T). Apply the pseudoenzyme penalty only when `site_source ∈ {MCSA, CDD_SITE}`.
- Validate: Benchmark 6B pseudoenzymes, plus active paralogs with conservative substitutions, targeting false "inactive" calls under 2%.

### BX-13 [MAJOR] MGE regime classification uses ad hoc ±15-ORF rules instead of established classifiers, and has no UNKNOWN state.
- Files: `PLAN.md` §4.2 `local_regimes`; `critique_3` §1.2 Component A.
- Failure: hallmark-window rules mislabel defense islands and satellites (PICIs), miss degraded prophages, and can't work on short Tier 3 contigs. A median of 2 to 6 ORFs means no ±15 window exists, so these contigs default to `CORE`. That is the exact error STAT-01 tries to prevent.
- Fix:
  - geNomad for plasmid, virus, and provirus scores and boundaries, plus CheckV for viral completeness and host-region trimming.
  - IntegronFinder 2 for integrons, CONJScan (MacSyFinder) and ICEfinder/ICEscreen for ICEs and IMEs, ISEScan for IS elements.
  - Add `UNKNOWN_SHORT_CONTIG` and record `regime_score`.
  - Restrict Tier 3 regime assignment to contigs at least 10 kb long, or to those where geNomad gives a confident call.
- Validate: geNomad's published benchmark, IMG/PR and PLSDB plasmids, prophage coordinates from PHASTER and Prophage-DB in RefSeq complete genomes. Report per-class F1.

### BX-14 [MAJOR] Gene trees for `H_{B|T_A}` are required but not budgeted.
- Files: `PLAN.md` §4.3 `gene_tree_homoplasy_b_on_a`; `critique_3` §1.2.
- Failure: this needs a protein phylogeny for every family A, including families with more than 700k members (Radical SAM). FastTreeMP is in the container but no stage builds trees or reconciles them. The §1.3 budget omits the step, which is likely the most expensive in the whole build.
- Fix: compute it on demand, only for the top-N pairs after `func_coupling_llr` ranking. Use FAMSA + VeryFastTree on subfamily NR50 representatives, capped at 5,000 leaves. Leave the column null in the base build, and add a budget line from pilot timing.
- Validate: timing and stability, measured as Robinson-Foulds variance across subsampling replicates, on 100 families.

### BX-15 [MAJOR] Overlap between sources double-counts genomes across tiers.
- Files: `PLAN.md` §2.2.
- Failure:
  - GTDB, RefSeq, GlobDB (built on GTDB plus 13 catalogues), proGenomes, and AllTheBacteria overlap heavily.
  - SPIRE MAGs overlap GTDB MAGs.
  - IMG/VR is mined from IMG/M, which the plan ingests separately.
  - Logan re-assembles the raw reads behind many Tier 1 assemblies, so one organism appears in Tier 1 and again in Tier 3.
  - The Stage 0 dedup works within a biome bucket and only on Tier 3, so it can't catch these duplicates. They inflate occurrence, phylum, and order counts.
- Fix:
  - First, dedup at accession level via BioSample/SRA run links from each source's metadata.
  - Second, at genome level with skani or sylph, at ≥99.5% ANI with high aligned fraction, applying a precedence order: RefSeq complete, then GTDB rep, then other Tier 1, then Tier 3.
  - Prophages already inside Tier 1 chromosomes should link to their IMG/VR record (`same_as`), not be counted twice.
  - Record `source_set` per contig.
- Validate: report the duplicate fraction per source pair, and check that it matches published overlaps such as GlobDB's composition table.

### BX-16 [MAJOR] Tier 1 and pangenome gating needs genome-quality tools the plan never runs.
- Files: `PLAN.md` §2.2, §4.4; BIO-08/09.
- Failure: the `hq_closed_or_multi_strain_mask` and core/pan planes depend on completeness, contamination, and chimerism, but no stage runs CheckM2 or GUNC. MAGs systematically lose plasmids, rRNA operons, and MGEs, so `pan_any_bits` built from MAGs underestimates accessory genes. Non-GTDB genomes (AllTheBacteria, proGenomes, Asgard and CPR assemblies) need GTDB-Tk or skani placement on the R226 tree before they can set bits. The mask is defined as "≥ 3 genomes" in critique_1 and "≥ 5" in PLAN §4.4.
- Fix:
  - Add a Stage 0.5: CheckM2 and GUNC on every Tier 1 assembly, then GTDB-Tk `classify_wf` (skani ANI screen) for genomes not already in GTDB.
  - Store `completeness`, `contamination`, `gunc_pass`, and `gtdb_species`.
  - Use one threshold for the mask.
  - Weight `pan_any` evidence from MAGs separately from isolate evidence.
- Validate: Benchmark 3 using real closed genomes subsampled into synthetic MAGs with CheckM2-estimated completeness.

### BX-17 [MAJOR] Terminator-based directon rules rely on TransTermHP outside its training domain.
- Files: `PLAN.md` §3.1 item 4; BIO-10.
- Failure: TransTermHP predicts intrinsic terminators tuned on E. coli and B. subtilis-like hairpins. It has known weak sensitivity in AT-rich clades and archaea, and it doesn't see Rho-dependent termination. The `feature_kind = 9` attenuator exemption depends on Infernal leader calls, which miss unannotated leaders.
- Fix: treat terminators as soft evidence (`terminator_score`) feeding `junction_topology_prev`, not as hard breaks. Report directon agreement with an operon predictor as a check (for example, operon-mapper-style intergenic distance plus conservation).
- Validate: RegulonDB transcription units (E. coli), term-seq maps (Dar 2016: B. subtilis, E. faecalis, L. monocytogenes), and an archaeal TSS/term-seq dataset (H. volcanii, S. acidocaldarius). Report directon-boundary precision and recall per clade.

### BX-18 [MODERATE] The Stage 0 sketch thresholds don't mean what the plan says.
- Files: `PLAN.md` §2.3 Filter 4.
- Failures:
  - A single 32-bit per-ORF syncmer hash is not "90% AA ID invariant". At 90% identity, the expected k=10 k-mer Jaccard is about 0.21, so a min-hash matches about 21% of the time. Murphy-10 raises that somewhat; the exact value is unverified.
  - For the 64-bit one-permutation b-bit sketch with b=1, ≤6 mismatches requires P(bit equal) ≥ 0.906, which gives J ≈ 0.81. At k=15, that is about 99% ANI, not 95%. The standard deviation is about 4 bits, so the test is noisy.
  - Jaccard on contigs of different lengths (fragments of one genome) is low even at 100% identity. Containment is needed.
  - On a second sighting, RocksDB has no representative yet, because the first sighting was written only to the Cuckoo filter. The insert step is missing.
- Fix: replace with a proper sketch (FracMinHash or containment via sourmash/skani logic, at least 128 to 256 hashes). State the measured ANI threshold. Write the representative on first sighting.
- Validate: on simulated contig pairs at 90 to 100% ANI with variable lengths, plot the ROC.

### BX-19 [MODERATE] Internal contradictions across files.
- `00-spec.md` §2 says "Ten Non-Negotiable Invariants" but lists 12.
- Obelisk rescue: ≥750 bp (spec, BIO-07 checklist) vs 700–1,500 bp (PLAN §2.3). Obelisks average about 1.1 kb (Zheludev 2024).
- Counts:
  - Contigs: 2.5B contigs vs 3.6B slots implies 1.44 slots per contig, which contradicts the "≥2 features" gate. critique_2 uses 5×10^8 contigs.
  - Windows: `window_id` range 2.8×10^9 (PLAN) vs 2.25×10^8 (critique_2). Because contigs start on window boundaries, window count is at least contig count.
  - NR90: about 220M (PLAN §3.2) vs 350–800M or 500M (critique_2). 8.5B raw ORFs is less than Tier 1 alone (about 2.4M genomes × about 4,000 genes ≈ 10B).
- IGR90 coverage: `-c 0.90` (PLAN) vs `-c 0.80` (critique_2).
- RoPE split: 256/256 dims (PLAN) vs 48/16 per head (critique_3).
- `directon_ids` and `regulon_locus_ids` are `uint16` (PLAN) but `uint32` (BIO-10). A 10 Mb chromosome can exceed 65,535 directons only in edge cases, but the types should match.
- `PRIOR-ART.md` still describes v1.0 features as current: Z-order Parquet, `gLM-lite`, `dom_idx_in_orf`, `ProstT5` on every dark protein.
- Hard guarantees such as "0.0% false hole calls" and "100% (200/200)" can't be verified as absolutes. Restate them as measured rates with confidence intervals.

### BX-20 [MODERATE] There's no pilot, truth-set gate, or annotation provenance.
- Failure: no layer has an acceptance test before the 4.7-day build, and slots carry no record of which tool version produced them. That breaks the malleability story: a `scan-new-models` result can't be compared with the base without provenance.
- Fix: add PLAN §0, "Pilot and acceptance gates". Run about 1,000 GTDB genomes (stratified by phylum, including CPR, DPANN, Asgard), 500 RefSeq complete genomes, 10,000 Logan runs, and the CAMI II datasets through the full stack, with each BX validation listed as a gate. Add `annotation_provenance.parquet` mapping `(layer, tool, version, db_release, params_hash)` to an id, and store that id per slot group.

---

## Ground-truth matrix (pilot gates)

| Layer | Truth set | Metric and gate |
|---|---|---|
| CDS calls | PGAP on 500 RefSeq complete genomes | stop-exact ≥ 98%, start-exact ≥ 85% |
| Genetic code | Codetta survey; Peters 2022 code-15 phages | per-genome code accuracy ≥ 99% |
| smORFs | Ribo-seq/Ribo-RET (Weaver 2019; Sberro 2019) | promoted-smORF FDR < 10% |
| ncRNA / subtype | Rfam seed; Breaker subtype sets | family recall ≥ 99% vs full cmsearch; subtype accuracy ≥ 98% |
| Directons / terminators | RegulonDB; term-seq (Dar 2016) | boundary F1 reported per clade |
| Domain boundaries | ECOD/CATH manual partitions; CPDB | boundary F1 within ±20 aa |
| Active sites | M-CSA | site-residue recall ≥ 95% |
| Fold | AFDB proteins with ECOD | T-group agreement per `fold_evidence` |
| MGE regime | geNomad benchmark; PLSDB; Prophage-DB | per-class F1 |
| Taxonomy (Tier 3) | CAMI II | order-level precision ≥ 90% at stated recall |
| Genome QC | CheckM2/GUNC on known-quality sets | agreement with published values |
| Dedup | simulated ANI pairs | ROC; stated threshold |
