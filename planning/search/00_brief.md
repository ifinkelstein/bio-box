# Search workloads: brief

Status: working brief for the search-infrastructure round (2026-10-29). Outputs of this round live in `planning/search/`.

## Why this round exists

The plan stores every contig, protein and IGR, and precomputes Pfam/NCBIfam/PHROGs/arCOG hits, clusters, tapes and pair statistics. It does not say how a person runs an ad hoc search: a jackhmmer from a seed, a profile-profile check, a structure search, an RNA motif search, a neighborhood census around new hits, and how those hits flow back into the tapes and statistics. The question for this round is whether the infrastructure serves the searches people actually run, and what has to change if not.

## Meta answers

**1. What a top professional would do that was not asked.**
- Treat search as a measured workload, not a feature. Derive the query taxonomy from what discovery papers actually did (methods sections and supplements), not from what tools exist.
- Extract a structured recipe per paper: seed, tool and parameters, iterations and inclusion threshold, target database and its size, search unit (protein, domain envelope, six-frame nucleotide, IGR, structure), curation, context step, downstream analysis, and what the paper published that would allow a replay.
- Pick a replay set of published searches whose hit lists, HMMs or alignments are in the supplements, and make "replay recall against the published hit list" an acceptance test for each stage.
- Decide what not to build: which steps stay on public services (HHpred, Foldseek server, CDD, EFI-EST/GNT, webFlaGs) and which must be local because they need the whole corpus.
- Plan the failure modes that dominate real searches: profile drift, promiscuous domains, gene-call errors creating false fusions or splits, clonal inflation, and E-values that are not comparable across database sizes.

**2. Rawest material.**
- The methods sections and supplementary files (parameter tables, HMMs, MSAs, accession lists) of the papers in `../bioinformatics-synteny/lit/` and `lit/supp/`, newest first.
- The 365 hypotheses in `../bioinformatics-synteny/data/synteny_hypotheses.sqlite`, especially `computational_testability_spec` and `first_bioinformatic_experiment`: these are the concrete queries the lab wants to run.
- Tool manuals and source: HMMER 3.4 user guide, pyhmmer docs, MMseqs2 wiki, HH-suite3 wiki, Foldseek README, Infernal user guide.
- The plan files themselves (`PLAN.md`, `STAGES.md`, `DEPLOYMENT-AND-MALLEABILITY.md`, `planning/*.md`), not summaries of them.

**3. Methods and tools.**
- Structured extraction into a recipe table with file:line evidence; a query taxonomy with counts; a replay-set selection.
- Search tools as practiced: HMMER3 (phmmer, jackhmmer, hmmsearch, hmmscan, nhmmer), pyhmmer, PSI-BLAST and DELTA-BLAST, MMseqs2 (iterative profile search, linclust, taxonomy, GPU search), DIAMOND, HH-suite3 (hhblits, hhsearch, HHpred), Foldseek and ProstT5, Infernal and CMfinder for RNA, CRISPRCasFinder/MinCED for arrays, MacSyFinder/DefenseFinder/PADLOC, EFI-EST/GNT sequence similarity networks, FlaGs/GCsnap/clinker for neighborhoods, AlphaFold2/3 and multimer screens downstream.
- Cost modeling from published throughput (cells/s, residues/s) against the stage sizes in `planning/sources_and_sizing.md`.

**4. What, other than us, shows the result is right.**
- Replay recall and precision against hit lists published in the papers' supplements, at the stage where the target taxa first enter.
- Concordance across independent tools on the same query (HMMER vs MMseqs2 profile vs HHblits; Foldseek for the dark cases).
- Measured wall time and memory on wilkcomp02 at S0, inside the envelope (200 GB cgroup, local `/tmp`, watchdog).
- A lab member running one real hypothesis from the synteny database end to end from the written recipe, without help.

## Rewritten prompt

> Treat `../bioinformatics-synteny` as the workload source: read its spec, report, scripts, tests and SQLite database, and classify all 365 hypotheses by the kind of search each needs. From the full text and supplements of its corpus (2022–2026 first, then 2017–2021, then older), and from modern papers outside it, extract a structured recipe for every search a paper actually ran, with file:line evidence. Build a query taxonomy with counts, and choose a replay set of published searches whose hit lists, HMMs or alignments are public. Map each recipe step onto the bio-box layers and say what is served, what is missing and what it costs on wilkcomp02 at S0, S5 and S6. Decide which steps must be local and which stay on public services. Then rewrite the infrastructure plan files to add what is missing, with acceptance tests: replay recall against published hit lists, cross-tool concordance, and measured cost inside the pod envelope.

## Additions from the user (2026-10-29)

- Planning only. Nothing in this round is implemented, installed or run on wilkcomp02.
- The search capability scales incrementally as tiers, starting from high-quality nucleotide datasets (complete, closed, curated assemblies) and widening only after each tier's replay and cost tests pass. The tiers go into the planning docs and map onto the stage ladder in `STAGES.md`.

## Recipe schema (shared by all literature agents)

One row per search a paper actually ran (a paper may have several). Columns:

| column | meaning |
|---|---|
| `paper_id` | file stem in `lit/`, or DOI for papers outside the corpus |
| `year` | |
| `question_type` | taxonomy code below |
| `seed` | what the search started from: sequence, HMM, MSA, structure, RNA motif, locus, gene family; and where it came from |
| `tool` | tool, version and the key parameters as written (E-value, inclusion threshold, iterations, coverage, bitscore cutoffs, `-Z`) |
| `iterations` | number of rounds and how convergence or stopping was decided |
| `target_db` | name, version, size (sequences, residues, genomes) as stated |
| `search_unit` | protein, domain envelope, six-frame nucleotide, contig, IGR, structure, CRISPR spacer |
| `curation` | how hits were filtered: length, coverage, catalytic residues, HHpred confirmation, tree placement, manual |
| `context_step` | neighborhood window (genes or bp), tool, how neighbors were annotated |
| `clustering` | how hits were grouped (MMseqs2, CD-HIT, SSN, tree) |
| `downstream` | trees, structures, AlphaFold, experiments |
| `scale` | number of hits, genomes or loci; compute or wall time if reported |
| `published_artifacts` | what the supplement contains that would let us replay the search (hit list with accessions, HMM, MSA, Newick, none) |
| `replayable` | yes / partial / no |
| `evidence` | `file:line` (or page) for the methods statement |

### Question-type codes

| code | question |
|---|---|
| Q-FAM | census of a known family with a fixed model (Pfam/TIGRFAM/COG HMM) |
| Q-ITER | iterative profile expansion from a seed (jackhmmer, PSI-BLAST, MMseqs2 iterative) |
| Q-REMOTE | profile-profile or remote homology (HHpred, HHblits, hhsearch) |
| Q-STRUCT | structure search (Foldseek, DALI, AF/ESMFold models) |
| Q-ARCH | domain architecture or fusion (Rosetta stone, domain lego) |
| Q-NBR | neighborhood or guilt-by-association (anchor then neighbors, enrichment, islands) |
| Q-PROF | phylogenetic profile or co-occurrence, mirror tree, co-evolution |
| Q-RNA | structured RNA motif discovery or search (CMfinder, Infernal, IGR comparison) |
| Q-NT | nucleotide-level search (tblastn on contigs, unannotated ORFs, nhmmer, repeats, CRISPR arrays and spacers) |
| Q-MOTIF | catalytic residue or sequence motif filtering on hits |
| Q-SSN | sequence similarity network or genome neighborhood network (EFI-EST/GNT) |
| Q-CLUST | corpus-scale clustering to find new families (FLSHclust, MMseqs2 cascaded, Foldseek clusters) |
| Q-SYS | system-level detection with rules (MacSyFinder, DefenseFinder, PADLOC, CRISPRCasTyper) |
| Q-OTHER | anything else; describe |
