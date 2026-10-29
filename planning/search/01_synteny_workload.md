# 01. The synteny project as a search workload

Status: written 2026-10-29 for the search-infrastructure round (`00_brief.md`). Source project: `../bioinformatics-synteny` (read-only). Per-hypothesis classification: [`synteny_hypothesis_query_types.tsv`](synteny_hypothesis_query_types.tsv).

## How this was done, and how far to trust it

- **Facts from the DB and files are exact.** Counts that come from the SQLite DB, supplement files or scripts were produced by a query or a script, and the query is given next to each number. "Spec" means `hypotheses.computational_testability_spec`. "First experiment" means `opportunity_ranking.first_bioinformatic_experiment`.
- **The question-type codes are judgments.** Five reader agents applied one shared rubric (the codes in `00_brief.md`, plus the flag definitions in the TSV header). Each read 73 hypotheses in full: spec, first experiment, syntenic architecture, families and accessions. Read the per-code counts as estimates with roughly ±10% inter-reader drift. The y/n flags record what the text literally asks for. `needs_custom_model` and `promiscuous_domain` are the readers' judgments.
- **The TSV has extra columns.** The 12 requested columns come first. They are followed by `seed_count`, `needs_custom_model`, `needs_iteration`, `needs_remote`, `needs_rna`, `needs_catalytic`, `promiscuous_domain`, `needs_neighborhood`, `window`, `taxon_scope_class`, `tools_named`, `dbs_named` and `paper_id`.
- **`named_models` came from a regex.** It is a regex pull of model accessions (`PF#####`, `TIGR#####`, `COG####`, `arCOG#####`, `RF#####`, `cd#####`, `IPR######`, `K#####`) from every text field of the hypothesis. It is not a check that the accession is correct; see §1.3.

---

## 1. What the project is, and what its contract teaches

### 1.1 The project

The project builds a corpus of 235 comparative-genomics papers from 1993 to 2026 across six lineages: Koonin/Makarova, Aravind, Weinberg/Breaker, Overbeek/Osterman/Haft, Sorek and Bondy-Denomy (`00-spec.md:18-25`). From that corpus it extracts **365 "atomic" synteny hypotheses** into `data/synteny_hypotheses.sqlite`. Each hypothesis is audited for experimental follow-up and ranked as an opportunity (`REPORT.md:5-10`).

**Tables.** `papers` (235), `hypotheses` (365), `followup_audit` (365), `opportunity_ranking` (365) and `paper_supplements` (235). These counts come from `select count(*)` on each table.

**Status split** (`select experimental_status,count(*) from followup_audit group by 1`):

| status | n |
|---|---:|
| OPEN | 99 |
| OPEN_ISH | 111 |
| HEAVY_EXPERIMENTAL_FOLLOWUP | 74 |
| HISTORICAL_BENCHMARK_SOLVED | 81 |

**Corpus files.** Full text exists as 118 PDF and 117 XML/BioC-derived `.txt` files (`corpus_manifest.json`, file_path suffix counts). Supplements were captured for 133 papers: 2,200 files, 1,562.5 MB (`lit/SUPPLEMENT_JOURNAL.md:5-9`). They are missing for 83 papers (54 publisher-only, 29 paywalled/JS-gated) and do not exist for 19 papers (`select supplement_capture_status,count(*) from paper_supplements group by 1`).

**Deductive types are not a closed vocabulary.** The DB holds 22 distinct `deductive_type` strings. The six enumerated in `00-spec.md:32` do not appear verbatim; the data uses "Guilt-by-Association" (134), "NOGD" (51), "Rosetta-Stone-Fusion" (49), "Promoter-Riboswitch-Cis" (36), "Island-Enrichment" (24) and others (`select deductive_type,count(*) from hypotheses group by 1`). Any consumer should treat the field as free text.

**Each hypothesis is a ready-made query.** It carries an anchor family, a partner family, accessions, a syntenic architecture, a testability spec and a first experiment, for example SYN-143 at `REPORT.md:120-134`. The specs are 115,567 characters in total and the first experiments 86,698 characters (`select sum(length(...))`), about 320 and 240 characters per hypothesis.

### 1.2 The provenance / zero-fabrication contract as built

| element | where | what it guarantees |
|---|---|---|
| `run_tool` is the only way to call an external tool | `00-spec.md:10`; `src/synteny/provenance.py:110-203` | It records argv, cwd, return code, wall time, the tool's own `--version`, and sha256 of declared inputs and outputs to `provenance/<run_id>.json` (`provenance.py:131-200`). It raises on a non-zero exit (`:191-192`) and when a declared output is missing (`:194-199`). |
| No literals that could become data | `provenance.py:17-18` | The module has no RNG and no numeric literals. |
| Every results column declares its source | `provenance.py:213-279`, `write_table` | Each column is `tool`, `input`, `derived` (needs a formula) or `label`. A numeric `label` column is refused (`:260-264`). Every TSV gets a `.provenance.json` sidecar. |
| Static tests | `tests/test_no_fabrication.py` | They ban `random` and `default_rng` (`:47-70`), direct `subprocess`/`os.system` outside provenance.py (`:74-100`), `check=False` (`:103-111`), and literal containers holding more than 12 floats (`:115-134`). They require a sidecar for every `results/*.tsv` (`:137-146`). |
| Corpus integrity | `tests/test_corpus_integrity.py` | Every paper file exists and its hash matches the DB (`:36-60`). Every `verbatim_quote` is found in the paper's text after whitespace and case normalisation (`:63-100`). Every audit row has an existing provenance record (`:96-97`). |
| Hypothesis ingestion | `scripts/03_build_database_and_audit.py:209-248` | A hypothesis whose quote cannot be resolved to exact text is dropped (`:239-242`). |

### 1.3 Where the contract held and where it did not

These findings bear directly on how user searches should be logged.

1. **Quotes are verified; every other field is unverified analyst text.**
   - The hypothesis rows come from hand-written or agent-written `data/hypotheses_batch_*.json` files. No script in `scripts/` or `src/` writes them; grep for `hypotheses_batch` finds only the reader at `03_build_database_and_audit.py:209-215`.
   - The sidecar declares `anchor_family`, `partner_family` and `accessions_mentioned` as kind `input`, "from paper text" (`03_build_database_and_audit.py:531`; `results/synteny_hypotheses_all.tsv.provenance.json`). Nothing checks this.
   - I checked every model accession against the source paper's own text (scratch script, exact substring):

     | accession type | hypothesis-model pairs | found in source text |
     |---|---:|---:|
     | `PF` | 449 | 156 |
     | `COG` | 212 | 128 |
     | `RF` | 50 | 7 |
     | `TIGR` | 101 | 79 |
     | `arCOG` | 57 | 53 |

   - 170 of the 241 hypotheses that name a model name at least one that is not in their source paper.
   - Most of these IDs look right from domain knowledge (SYN-143: PF01541 GIY-YIG, PF00929 DEDD, PF00730 HhH-GPD). They were nevertheless supplied from memory, not from the text. At least one is malformed: SYN-265 "PF0466" (`REPORT.md:31`).
   - Only 1,539 of 2,139 `accessions_mentioned` tokens appear verbatim in the source text.
   - **Adopt:** resolve every model ID and accession against the release's own name table (Pfam-A.hmm.dat, Rfam family.txt, COG/arCOG definitions, NCBI) when a search is submitted. Reject the search, or flag it, when the name disagrees.

2. **The tool-derived audit columns are degenerate, and the tests cannot see it.**
   - `specific_exp_citing_count = 0` and `pdb_hit_count = 0` for all 365 hypotheses, and `pdb_ids = 'None'` for all 365 (`select experimental_status,count(*),sum(specific_exp_citing_count=0),sum(pdb_hit_count=0) from followup_audit group by 1`).
   - That holds even for the 74 HEAVY_EXPERIMENTAL_FOLLOWUP rows. The curated `key_followup_papers` text names PDB entries in 5 rows.
   - The cause is over-specific queries. For example, SYN-001's OpenAlex query is `("MG252 cysS methyltransferase" OR "RlmB YacO 23S rRNA Gm2251 cysS") AND (<experimental terms>)` restricted to works citing the paper (`provenance/1793242743140-968-curl.json`, field `command`).
   - Each call returned 0 with exit 0, was logged correctly, and passed every test.
   - As a result `opportunity_score` (`03_build_database_and_audit.py:360-370`) depends only on the analyst-supplied `experimental_status` and `in_padloc_or_defensefinder` labels. It takes 11 distinct values, and 87 hypotheses tie at the 99.5 cap (`select count(distinct opportunity_score), sum(opportunity_score=99.5) from opportunity_ranking`).
   - **Adopt:** provenance proves that a call happened, not that it measured anything. Every user search recipe must carry at least one positive control, a known member of the family or locus. The run is marked failed if the control is not recovered. Any result column that is constant across all rows gets a warning.

3. **Fallback provenance.**
   - Papers with no OpenAlex record inherit another paper's run ID: `fallback_run_id` at `03_build_database_and_audit.py:165`, `:190` and `:200-201`.
   - 23 papers share a run ID with another paper, against 219 distinct run IDs (`select count(distinct openalex_run_id) from papers`).
   - 12 papers have `doi='NA'`.
   - Year and first author were parsed from the filename, with year 2015 as the default (`:172-174`).
   - **Adopt:** never fill a provenance field with a placeholder. A missing record is an error, not a default.

4. **File identity is not checked against content.**
   - Integrity is a sha256 of whatever was downloaded. Nothing checks that the file is the paper its stem names.
   - Examples where the text is a different paper:
     - `lit/2002_Anantharaman_Koonin_Aravind_RNA_Metabolism_Proteins_Evolution_NAR.txt` begins "Expanded Genomic Sampling Refines … Sulfur Metabolisms in the Desulfobulbales". SYN-148, the aprAB/qmoABC hypothesis, was extracted from it.
     - `2026_Pittman_…Type_IV_C_CRISPR…` contains a Monkeypox LAMP paper (SYN-260).
     - `2001_Wolf_…Directons…` contains a Xanthomonas TALome paper.
     - `2019_Patz_…Tailocins…` has a Cancer Research DOI (`10.1158/0008-5472.can-19-0049`).
   - A keyword-overlap screen (scratch script; stem keywords vs. the first 600 characters of text) flags 56 papers. About half of those are false alarms from PMC or "HHS Public Access" banners. By hand I count at least 12 that are clearly a different paper.
   - In 25 papers the year in the stem disagrees with `papers.year`, and in 14 the first author disagrees (scratch script over `papers`).
   - **Adopt:** downloaded inputs, including reference databases and published hit lists used for replay, are identified by content as well as hash. Examples: the DB's own release string, the first FASTA header, a title check for a paper. The identity assertion is stored next to the hash.

5. **Releases drift.**
   - The specs pin GTDB r220: 90 hypotheses say "GTDB r220" and 172 say "GTDB" (`… where spec||first_exp like '%GTDB r220%'`).
   - The infrastructure plan targets R226 (`STAGES.md:91-94`). No hypothesis mentions R226.
   - **Adopt:** a recipe stores the release it was written against. Execution records the release it ran against, and the difference is shown.

### 1.4 Concrete contract for user searches

The parts to take over directly:

| rule | origin in synteny | form in the search layer |
|---|---|---|
| One wrapper launches every search tool (hmmsearch, jackhmmer, mmseqs, hhblits, foldseek, cmsearch) | `run_tool` | A `search_run` record holds argv, tool version, DB name, release and sha256 of the DB index, query file hash, output hash, wall time, peak RSS, host and cgroup, and the exit code. A non-zero exit raises. |
| Every hit row carries its run ID | `ColumnSource(kind="tool", run_ids=…)` | Hit tables carry `run_id` per row. Neighborhood and tape tables derived from hits carry the parent `run_id` list. |
| Column sources are declared, and numeric labels are refused | `write_table` | Hit tables declare which columns are tool output (E-value, bitscore), derived (coverage, `-Z`-normalised E-value) or curated label (keep/drop with a reason). A human curation is a label and can never become a number. |
| No silent failure | ban on `check=False`; `ToolFailure` | Searches that time out or hit the memory limit are recorded as failed, never as "0 hits". |
| Quote-level grounding for seeds | `verify_quote_in_text` | A seed taken from a paper stores the `paper_id`, `file:line` and the accession string, verified against the text. A seed taken from a DB stores the DB release and entry. |

The additions synteny lacked, from §1.3:

- resolve model IDs against the release's name table when the search is submitted;
- require a positive control for every search;
- warn when a result column is degenerate;
- never use fallback provenance;
- verify input identity by content;
- record the release a recipe was written against and the one it ran against.

---

## 2. Classification of all 365 hypotheses

### 2.1 Primary code

`primary_code` is the single search question the test hinges on (TSV column `primary_code`).

| primary | n | % | OPEN | OPEN_ISH | HEAVY | HIST |
|---|---:|---:|---:|---:|---:|---:|
| Q-NBR | 143 | 39 | 38 | 41 | 40 | 24 |
| Q-PROF | 59 | 16 | 8 | 22 | 12 | 17 |
| Q-ARCH | 50 | 14 | 18 | 9 | 6 | 17 |
| Q-STRUCT | 29 | 8 | 12 | 10 | 3 | 4 |
| Q-NT | 23 | 6 | 2 | 11 | 5 | 5 |
| Q-RNA | 22 | 6 | 10 | 4 | 1 | 7 |
| Q-MOTIF | 18 | 5 | 4 | 9 | 2 | 3 |
| Q-SSN | 7 | 2 | 5 | 1 | 1 | 0 |
| Q-OTHER | 5 | 1 | 0 | 2 | 2 | 1 |
| Q-SYS | 4 | 1 | 0 | 1 | 2 | 1 |
| Q-FAM | 2 | <1 | 0 | 0 | 0 | 2 |
| Q-ITER | 2 | <1 | 1 | 1 | 0 | 0 |
| Q-CLUST | 1 | <1 | 1 | 0 | 0 | 0 |
| Q-REMOTE | 0 | 0 | | | | |

The five Q-OTHER rows:

- SYN-004 is a whole-genome gene-order (DCJ) reconstruction.
- SYN-075 is an intragenic gene-conversion scan on codon alignments.
- SYN-196 needs replication-strand attributes, not homology.
- SYN-200 is a eukaryotic promoter hypothesis.
- SYN-280 is downstream metabolic choke-point inference (TSV `notes`).

How primary code lines up with deductive type (TSV cross-tab):

| deductive type | n | Q-NBR | Q-PROF | Q-ARCH | Q-RNA |
|---|---:|---:|---:|---:|---:|
| Guilt-by-Association | 134 | 81 | | | |
| NOGD | 51 | 15 | 22 | | |
| Rosetta-Stone-Fusion | 49 | | | 38 | |
| Promoter-Riboswitch-Cis | 36 | | | | 13 |

Promoter-Riboswitch-Cis also has 8 Q-NT, 8 Q-NBR and 6 others. The deductive type predicts the search, but not tightly.

### 2.2 How often each code appears anywhere in the chain

| code | hypotheses with the code in the chain | % |
|---|---:|---:|
| Q-NBR | 324 | 89 |
| Q-STRUCT | 267 | 73 |
| Q-FAM | 262 | 72 |
| Q-PROF | 138 | 38 |
| Q-MOTIF | 126 | 35 |
| Q-ARCH | 98 | 27 |
| Q-NT | 89 | 24 |
| Q-RNA | 57 | 16 |
| Q-CLUST | 29 | 8 |
| Q-SYS | 19 | 5 |
| Q-OTHER | 19 | 5 |
| Q-SSN | 13 | 4 |
| Q-REMOTE | 11 | 3 |
| Q-ITER | 9 | 2 |

**Q-STRUCT is mostly the last step, not a search.** It is the last step in 230 of the 267 chains that contain it. In those it means AlphaFold3 or AF-Multimer modelling and docking of the found complex. Structure is used as a *search* (Foldseek, or "structure as query") in only 23 hypotheses (TSV `tools_named`/`notes` match "foldseek"). Infrastructure consequence: GPU structure prediction of a handful of complexes per hypothesis is a downstream job; structure search against a corpus-scale structure DB is a minority need (§4).

**Chain shape.**

- Chains are 2–7 steps long: 4 steps in 193, 3 in 78, 5 in 70, 6 in 13, 2 in 10, 7 in 1.
- The entry step is Q-FAM in 254 (70%), Q-RNA in 27, Q-NBR in 25, Q-ARCH in 23, Q-NT in 12, Q-SYS in 9 and Q-ITER in 8.
- There are 212 distinct chains.

### 2.3 Chains

Most common full chains:

| chain | n |
|---|---:|
| Q-FAM→Q-NBR→Q-MOTIF→Q-STRUCT | 19 |
| Q-FAM→Q-NBR→Q-PROF→Q-STRUCT | 18 |
| Q-FAM→Q-ARCH→Q-NBR→Q-STRUCT | 14 |
| Q-FAM→Q-NBR→Q-STRUCT | 12 |
| Q-FAM→Q-NBR→Q-NT→Q-STRUCT | 11 |
| Q-FAM→Q-PROF→Q-NBR→Q-STRUCT | 8 |
| Q-ARCH→Q-NBR→Q-STRUCT | 7 |
| Q-FAM→Q-NBR→Q-STRUCT→Q-MOTIF | 6 |
| Q-FAM→Q-PROF→Q-NBR | 5 |
| Q-FAM→Q-ARCH→Q-NBR→Q-MOTIF→Q-STRUCT | 5 |
| Q-FAM→Q-PROF→Q-NBR→Q-MOTIF→Q-STRUCT | 5 |
| Q-ARCH→Q-NBR→Q-NT→Q-STRUCT | 5 |
| (200 other chains, ≤4 each) | 250 |

**Search skeletons.** Dropping the terminal modelling step (Q-STRUCT) and Q-OTHER collapses the 212 chains to 161 skeletons:

| skeleton | n | what the infrastructure must do |
|---|---:|---|
| Q-FAM→Q-NBR→Q-MOTIF | 27 | fixed-model hits → neighbors → residue check on hits |
| Q-FAM→Q-NBR→Q-PROF | 26 | hits → neighbors → co-occurrence statistics over genomes |
| Q-FAM→Q-ARCH→Q-NBR | 17 | domain-level hits on one ORF (fusion) → neighbors |
| Q-FAM→Q-NBR | 15 | plain anchor-and-neighbors census |
| Q-FAM→Q-NBR→Q-NT | 14 | neighbors, then IGR / ORF / codon-overlap work on nucleotides |
| Q-FAM→Q-PROF→Q-NBR | 13 | phyletic pattern first, then locus check |
| Q-FAM→Q-ARCH→Q-NBR→Q-MOTIF | 8 | |
| Q-ARCH→Q-NBR | 7 | |
| Q-RNA→Q-NBR→Q-PROF | 5 | CM search → flanking genes → profile |
| Q-ARCH→Q-NBR→Q-NT | 5 | |
| Q-FAM→Q-NBR→Q-CLUST | 4 | |
| Q-RNA→Q-NBR→Q-MOTIF | 4 | |
| Q-RNA→Q-NBR | 4 | |
| Q-ITER→Q-NBR→Q-MOTIF | 3 | |

**Reading.** The dominant workload is *anchor by a fixed model, pull the neighborhood, then do one of three things*: filter on residues, compute co-occurrence across genomes, or go down to nucleotides. That shape is exactly "Q-FAM→Q-NBR→X". It is what the precomputed annotation, neighborhood tapes and pair statistics already target, so the core of the plan fits the workload. Section 4 is about the X.

### 2.4 Seed kind

`seed_kind` is the kind of thing the search starts from.

| seed_kind | n | of which `needs_custom_model=y` |
|---|---:|---:|
| fixed_model (a named Pfam/TIGR/COG/arCOG/Rfam) | 190 | 13 |
| named_family_no_model (gene names only) | 60 | 38 |
| locus_operon (a whole cluster as the seed) | 33 | 12 |
| rna_motif | 29 | 7 |
| specific_protein (locus tag / accession) | 26 | 11 |
| domain_combination | 15 | 4 |
| other | 10 | 4 |
| structure | 2 | 0 |

---

## 3. Tools, databases, models and taxa named in the specs

### 3.1 By SQL over the spec and first experiment

Counts are of hypotheses, n = 365. The query is:

`with t as (select h.computational_testability_spec||' '||o.first_bioinformatic_experiment x from hypotheses h join opportunity_ranking o using(hypothesis_id)) select sum(x like '%TERM%') from t`

**Databases and releases**

| term | n | note |
|---|---:|---|
| GTDB (any) | 172 | |
| GTDB r220 | 90 | No R226, r214 or r207 mentions. |
| RefSeq | 28 | |
| metagenom* | 14 | |
| IMG/VR | 12 | 6 say "IMG/VR v4" (TSV `dbs_named`). |
| PDB | 11 | 44 including family text (regex `\bPDB\b`). |
| AFDB / AlphaFold DB | 7 | regex |
| PLSDB / plasmid DB | 4 | regex |
| SPIRE | 4 | TSV `dbs_named` |
| UHGG, Logan, MGnify, INPHARED | 0–2 each | MGnify and INPHARED: 0 (regex). |

**Models**

| term | n | note |
|---|---:|---|
| Pfam / PF0… / PF1… | 96 | |
| COG | 80 | |
| Rfam / RF0… / RF1… | 32 | |
| TIGR | 21 | |
| arCOG | 12 | |
| KEGG / KO | 4 | regex |
| InterPro | 4 | regex |
| CDD | 2 | |
| NCBIfam, PHROGs, eggNOG | 0 | regex |

**Search tools**

| term | n | note |
|---|---:|---|
| "HMM" (any) | 34 | |
| EFI / SSN | 19 | |
| HMMER family | 28 | regex: HMMER, hmmsearch, hmmscan, phmmer, jackhmmer, nhmmer, "profile HMM" |
| Infernal / cmsearch | 13 | |
| R-scape | 5 | regex |
| antiSMASH | 6 | regex |
| FlaGs / GCsnap / clinker / CORASON | 6 | regex |
| MMseqs | 4 | |
| DefenseFinder / PADLOC | 4 | |
| HHpred / HHblits | 3 | |
| tblastn | 2 | |
| jackhmmer | 1 | |
| PSI-BLAST, DIAMOND, nhmmer, CMfinder, DALI, CRISPRCasTyper | 0 | |

**Structure**

| term | n | note |
|---|---:|---|
| AlphaFold (any) | 226 | AF3 alone: 220 (regex) |
| docking / dock | 68 | |
| AF-Multimer | 63 | |
| Foldseek | 22 | |

**Analysis**

| term | n | note |
|---|---:|---|
| catalytic / active site | 102 | 158 with "motif / inactive / degenerate" (regex) |
| co-occurrence / phylogenetic profil | 61 | 109 with presence/absence, anticorrelation, complementarity (regex) |
| intergenic / 5'-UTR / promoter | 51 | 77 with "leader / upstream region / IGR" (regex) |
| mirror tree | 10 | 27 with co-evolution / DCA / tree congruence (regex) |
| six-frame / unannotated | 6 | 11 with "small ORF / missed ORF" (regex) |
| CRISPR spacer / array | 14 | 22 with "protospacer / repeat array" (regex) |

### 3.2 What the reader agents extracted (TSV `tools_named`, `dbs_named`, distinct per hypothesis)

**Tools**

| tool | n |
|---|---:|
| AlphaFold3 | 225 |
| AlphaFold-Multimer | 27 |
| Foldseek | 22 |
| Infernal / cmsearch / cmbuild | ~22 |
| MirrorTree | 5 |
| R-scape | 5 |
| PGLS | 4 |
| EFI-GNT | 4 |
| HHpred | 4 |
| FlaGs | 4 |
| MMseqs2 | 4 |
| DefenseFinder | 3 |
| GREMLIN | 3 |
| cblaster | 3 |
| antiSMASH | 3 |
| DeepTMHMM | 3 |
| one or two each: MacSyFinder, ESMFold, SIMBAL, DCA, RODEO, ARTS, CORASON, RNAfold, MEME, jackhmmer, tBLASTn, PAML/HyPhy | |

**Databases**

| database | n |
|---|---:|
| Pfam | 91 |
| GTDB r220 | 89 |
| GTDB (no release) | 84 |
| COG | 30 |
| RefSeq | 29 |
| PDB | 20 |
| TIGRFAM | 12–15 |
| UniProt | 8 |
| InterPro | 7 |
| arCOG | 7 |
| IMG/VR | 12 |
| AFDB | 5 |
| KEGG | 5 |
| SPIRE | 4 |
| ECOD | 3 |
| PLSDB | 3 |
| IMG/M | 3 |
| one or two each: SCOP, SFLD, UniRef50/90, MIBiG, BacDive, Logan, SRA, UHGG v2 | |

**Hypotheses that name a published supplement as an input.** A few hypotheses name a supplement file as the seed input, for example SYN-343 "published CorA.tar.gz tree and FASTA" and SYN-364 "published media-1.xlsx". These are natural replay hooks (§5).

### 3.3 Models: named, missing and custom

- **Named models** (regex `named_models` column): 242/365 hypotheses name at least one model accession.

  | model type | hypotheses naming one |
  |---|---:|
  | PF | 175 |
  | COG | 76 |
  | RF | 28 |
  | TIGR | 25 |
  | arCOG | 13 |
  | IPR | 6 |
  | K (KEGG) | 3 |

- **No model at all:** 123/365 name no model accession.
- **Custom profile needed** (`needs_custom_model=y`): 89/365 (24%), 60 of them OPEN or OPEN_ISH. These are families named only by gene name, or newly defined in the source paper, so a profile must be built from seeds before any census. Examples:
  - SYN-148: aprAB/qmoABC by gene name.
  - SYN-221: an orphan gltD RNA motif. The CM must be built from 95 representatives in a .docx supplement.
  - SYN-268: hypermodification enzymes; the spec says "build HMMs".
  - SYN-326: the EH_Signature family.
  - SYN-329: MazG-C, defined only by a DxYRxHDxxH motif.
- **Model ID present but unverified:** see §1.3. Only 35% of PF IDs (156/449) and 14% of RF IDs (7/50) occur in the hypothesis's own source text. The search layer must resolve these IDs against the release before running.

### 3.4 Taxa

**Scope class** (TSV `taxon_scope_class`):

| scope | n |
|---|---:|
| all_bacteria | 100 |
| clade (phylum to genus) | 92 |
| all_prokaryotes | 78 |
| mixed (e.g. Actinomycetota + Methanococcales; host + phage) | 51 |
| all_archaea | 19 |
| phage_virus | 9 |
| single_genus | 6 |
| plasmid_MGE | 6 |
| metagenome | 4 |

The 92 clade-scoped hypotheses name specific clades: Bacillota, Actinomycetota, Desulfobulbales, Rhizobiales, Streptococcus, Sulfolobales, Halobacteria, Thermoplasmatota, Asgard and others.

**Stated genome counts.** These are rare. About 30 hypotheses state a genome or locus count. That is an estimate, counted after removing regex matches that were really release numbers ("r220 archaeal genomes") or accession fragments (scratch `scale_counts.py`). Examples:

| count stated | hypothesis |
|---|---|
| 5,000 representative bacterial genomes | SYN-001 |
| 5,000 GTDB MAGs | SYN-021 |
| 2,000 Bacillota | SYN-002 |
| 2,000 Lactobacillales | SYN-289 |
| 1,599 Streptococcus | SYN-024..028 |
| 10,000 GTDB Actinomycetota | SYN-085 |
| 20,000 complete bacterial genomes | SYN-103 |
| 50,000 GTDB genomes | SYN-061, 068, 070, 278 |
| 500,000 proteomes | SYN-062 |
| 25 genomes | SYN-263 |
| 80 Buchnera | SYN-004 |

---

## 4. What needs more than precomputed fixed-library annotation

**Baseline served by the precomputed layers.** The precomputed layers are Pfam/NCBIfam/PHROGs/arCOG hits, clusters, neighborhood tapes and pair statistics. They can answer a hypothesis if all of these hold:

- the primary code is Q-FAM, Q-NBR, Q-PROF or Q-ARCH;
- no custom model, iteration, remote homology, RNA, nucleotide-level or IGR step is needed.

**133/365 (36%)** pass that test. Only **21** remain once structure modelling and catalytic-residue checks are excluded too. So 94% of hypotheses need *something* beyond the fixed library (353/365 with any flag set). Most of that "something" is cheap and per-hit, not a new corpus-scale index.

| need | n (open + open-ish) | what it requires of the infrastructure | example IDs |
|---|---|---|---|
| **Custom seed / profile build** | 89 (60) | Build an HMM or CM from user seeds, a published MSA, or the first-round hits. Run it against the stage's protein (or IGR) set, and store the model, its build inputs and its hash as a first-class object next to Pfam. | SYN-010, 013, 014, 148, 166, 221, 225, 241, 268, 326, 329, 346, 355–358, 364 |
| **Iteration** (explicit jackhmmer, or a profile from one protein expanded) | 7 explicit; 9 with Q-ITER in the chain; **91** counting custom-model rows, which imply at least one build-search-rebuild cycle | An iterative loop with a per-round inclusion threshold, a frozen per-round hit set, and a drift check (does the round-n profile still hit the seed and the positive control?). The `-Z` must be fixed to the stage size so E-values are comparable across rounds and stages. | SYN-013, 014 (jackhmmer from BH3254), 091, 107, 151, 179, 278 (CM rebuild) |
| **Remote homology** (HHpred/HHblits, fold assignment) | 12 (5) | Profile-profile search on a few hundred to a few thousand dark ORFs per hypothesis. Small enough to leave on public HHpred, or to run hhblits locally against UniRef30/PDB70 (decision for a later file). | SYN-019, 020, 089, 090, 093, 138, 139, 179, 192, 235, 284, 291 |
| **Structure search** (Foldseek or structure-as-query) | 23 (12) | Foldseek over AFDB/ESM-Atlas or a predicted structure set of dark ORFs in found loci. Query sets are small (tens to thousands per hypothesis). | SYN-014, 048, 090, 125, 139, 160, 186, 209, 235, 284, 291, 326, 329 |
| **Structure modelling** (AF3 / AF-Multimer / docking) | 272 flagged (180) | A downstream GPU queue, not search: complexes of the found anchor + partner, often with RNA or DNA ligands (AF3). Out of scope for the search layer but must accept hit-set IDs as input. | most rows; e.g. SYN-001 (AF3 protein–RNA), SYN-143 |
| **Structured RNA** (CM build, cmsearch, R-scape, RNAfold) | 60 (31); Q-RNA primary 22 | IGR extraction with strand and distance to the downstream start codon. Infernal over IGRs, not whole contigs, with stage-sized `-Z`. CM build from Stockholm. Covariation check. Rfam IDs (28 hypotheses) must be resolved to the current Rfam release. | SYN-112–120, 221–226, 278–288, 340 |
| **IGR without an RNA structure** (promoters, operators, terminators, TF motifs) | 47 (30) | Store the IGR as a first-class sequence with flanking-gene IDs. Scan for DNA motifs (MEME/FIMO-style PWM, bipartite sigma-54 motifs, LexA boxes), intergenic distance, and bidirectional-promoter geometry. | SYN-002, 012–015, 063, 126, 143, 154, 191, 196, 197, 242, 351 |
| **Nucleotide-level** (six-frame, unannotated or small ORFs, tblastn, codon overlaps, CRISPR arrays/spacers, repeats, TIR/TSD, DGR) | 120 (64); Q-NT primary 23 | Contig sequence must be searchable directly, not only through called proteins. Needs tblastn or six-frame hmmsearch around loci, short-ORF calling, CRISPR array detection, spacer-to-protospacer matching against S1/S2/S4, and repeat finding. SYN-067 (missed mftA ORF) is an explicit positive control proposed in the DB. | SYN-027 (frameshift alleles), 053, 067, 088, 155, 159, 171, 176, 177, 242, 245, 260, 333, 358 |
| **Catalytic-residue / motif checks on hits** | 131 (77); Q-MOTIF primary 18 | Align hits to a reference (hmmalign to the model or a structure-anchored MSA) and report residue identity at named positions per hit. Many tests are *absence* tests ("catalytically inactive", "eroded metal site"), so the per-position report must distinguish absent from not aligned. | SYN-143 (catalytic triads), 262 (ZorE PD-(D/E)xK), 329, 351 (absent catalytic residues), 364 (eroded YprA nuclease) |
| **Promiscuous domains** (HTH, ABC, AAA+, P-loop, Rossmann/MTase/SDR, GNAT, HAD, radical SAM, RR/HK, TPR) | 143 (81) | The anchor or partner returns floods of hits. The workload needs architecture-aware filters (co-domain required, length window, fused-vs-split), neighbor-conditioned counting, and caps on neighborhood fan-out. Without them a Q-FAM→Q-NBR census on SYN-001's COG0566 (SpoU) or SYN-306's DEDD exonucleases returns the whole family. | SYN-001, 063, 099, 100, 106, 123, 125, 143, 149, 299, 364 |
| **Domain-level / fusion calls** (Q-ARCH in chain) | 98 | Domain envelopes per ORF from the fixed library (served). Also split-vs-fused equivalence (SYN-299, 343) and circular permutation (SYN-052), which the per-protein best-hit view loses. | SYN-052, 062, 102, 143, 149, 299, 306, 343 |
| **Absence queries** | — (readers flagged SYN-018, 048, 049, 343 and others in notes) | "Loci with X and Y but no Z within the window" requires the neighborhood tapes to support negation with a known denominator. | SYN-018, 048, 049, 122 (no proteasome), 343 |
| **Systems rules** (Q-SYS in chain) | 19 (5) | DefenseFinder / PADLOC / MacSyFinder / CRISPRCasTyper models run as precomputed system calls. | SYN-073, 093, 155, 157, 161, 220, 231, 339, 360 |
| **SSN / GNN** (Q-SSN in chain) | 13 (7) | All-vs-all on a hit set of up to about 10^4 proteins, plus a neighborhood join. Can stay on EFI-EST/GNT for small sets, or run locally with MMseqs2 all-vs-all. | SYN-022, 030, 031, 034, 199, 336, 337, 345 |
| **Corpus clustering** (Q-CLUST in chain) | 29 (10) | Uses the precomputed clusters (served). New-family discovery from dark neighbors needs re-clustering of a hit-neighborhood set. | SYN-024, 071, 072, 282 |

**Summary of extra capabilities.** The infrastructure needs four extras to cover the workload:

- (a) user-built HMM and CM objects, with an iterative loop and a fixed `-Z`;
- (b) first-class IGR and contig-level nucleotide search around loci;
- (c) per-hit residue reports from an alignment to the model;
- (d) promiscuity-aware architecture filters and negation over neighborhoods.

Structure search and remote homology are small-volume needs: 23 and 12 hypotheses, with small query sets. AF3 modelling is the largest single "tool", named in 220–225 hypotheses, but it is downstream of search.

---

## 5. Replay candidates from `lit/supp/`

### 5.1 What the supplements contain

**Method.** I walked all 133 supplement directories (2,200 files) with scratch scripts that sniff file heads, open `.xlsx` sheets, list zip and tar members without extracting them into the source repo, and regex-count accessions. The counts below come from those scripts. The `.xls` files (BIFF) were scanned as raw bytes only.

**By file type.**

| type | files | MB |
|---|---:|---:|
| images (gif 818, jpg 797, eps 14, tif 13, webp 6, jpeg 4) | 1,652 | |
| pdf | 142 | 319 |
| xlsx | 77 | 94 |
| xls | 40 | 43 |
| txt | 34 | 123 |
| html/htm | 50 | |
| zip | 28 | 301 |
| `.ali` alignments | 25 | |
| doc/docx | 39 | |
| tar.gz/tgz/gz/7z | 26 | |
| csv | 13 | |
| sto | 8 | |
| fa/fasta/fna | 8 | |
| tre | 3 | |
| hmm | 1 | |

Images are three-quarters of all files.

**Sequence artefacts inside archives.** Most of them live inside archives:

- **HMMs:** 2 papers.
  - `2021_Sutter…/41467_2021_24126_MOESM7_ESM.zip` holds `combined_hmms.hmm` (1,220 HMMs, 147 MB) and `ep_hmms.hmm` (104 HMMs).
  - `2014_Axen…/pcbi.1003898.s004.hmm` holds 1 BMC HMM.
- **RNA Stockholm alignments:** 4 papers.
  - Weinberg 2010: `gb-2010-11-3-r31-S4.GZ`, a tar of 121 `.sto`, plus S5 with 121 cleaned `.sto`, plus S10 with the IGR pipeline source code (`bliss-scripts/add_igr.perl`, `add_genomecontext.perl`, `add_terminators.perl`).
  - Weinberg 2017: `gkx699_supp.zip` holds 224 motif-only `.sto`, 224 annotated `.sto` and 52 bSCRAP `.sto`.
  - Corbino 2005: 5 `.sto` (serC, speF, suhB, ybhL, metA).
  - Sherlock 2018 PRPP: 2 `.sto`.
- **Protein MSAs:** about 15 papers. Examples:
  - Wolf 1999: 22 aaRS `.ali`.
  - Makarova 2009 pAgo and TA: 3 + 5 `.ali`.
  - Makarova 2011 CRISPR: S1/S5/S6 multi-family alignments.
  - Yutin 2012 tubulin: 89-sequence MSAs.
  - Axen 2014: 133- and 352-sequence MSAs.
  - Makarova 2014 CARF: SF3/SF6/SF7.
  - Shmakov 2018: 128 `icity*.sr` profile alignments.
  - Jimmy 2020 toxSAS: `sd05.txt`.
  - Kaur/Iyer 2022: 14 `.aln`.
  - Faure 2023: TnsC 6,384 and TniQ 4,916 aligned FASTA.
- **FASTA:** 7 papers.
  - Wolf 2012: `ar14.fa.gz`, 394,110 archaeal proteins.
  - Shmakov 2018: Cas10 4,968; CysH 11,786; STAND 2,389; RT 1,196; CorA 547; CARF 50; 16S 3,308.
  - Shmakov 2019: Vicinity.faa ×3.
  - Yutin 2018 crAss: 2,684 proteins.
  - Peters 2017: 804 spacers (nt).
  - Gussow 2020: 3,718 proteins.
  - Makarova 2024 ESCRT: 12,399 proteins.
- **Trees:** 7 papers. Shmakov 2018 (7 FastTree), Mestre 2020 (6 Newick, 544–2,253 leaves), Faure 2023 (TnsC 6,383, TniQ 4,915, TnsF 1,092 leaves), Makarova 2024 (2 IQ-TREE), Sutter 2021 (6 phyloXML), Wolf 2012 and Koonin 2001 (genome trees).
- **Machine-readable hit or locus lists with accessions:** about 28 papers (estimate from the per-file accession counts). One more paper uses IMG OIDs only: Ofir 2018 DISARM, with IMG gene and taxon OIDs that need mapping.

**Not replayable from the supplement alone.** This covers 83 papers with gated or publisher-only supplements, and papers whose supplement is PDF-only. Examples: Makarova 2025 CRISPR classification (3 PDFs), Lee 2026 PRRs (18 MB PDF), Yee 2026 END nucleases (34 KB PDF), Kretsch 2025. These need PDF table extraction first.

**Two caveats for any replay.**

1. **Strain-level lists against species-level stages.** Hit lists are strain-level GenBank/RefSeq assemblies. At S5 (species representatives), recall must be scored only on hits whose assembly is a GTDB R226 species representative, or by species-level mapping. At S6 it is scored on hits whose assembly is in GTDB.
2. **Retired GI numbers.** Lists from 2009–2018 identify proteins by GI number (for example Shmakov 2018 `*.faa` headers, Makarova 2011 islands, Atkinson 2011 `gi|…|ref|YP_…|`). NCBI has retired GIs, so map these by sequence (exact or ≥95% MMseqs2 match against the stage proteome), not by identifier.

### 5.2 Golden-genome (S0) footholds

S0 holds 14 genomes (`STAGES.md:128`): MG1655, *B. subtilis* 168, T4, λ, *M. genitalium*, SR1, Lak, *M. acetivorans*, *H. volcanii*, *S. acidocaldarius*, H37Rv, *B. thetaiotaomicron*, *S. pyogenes* SF370 and an Obelisk.

I scanned every supplement for these genomes' accessions and names (scratch `golden_scan.py`; counts are string mentions, not de-duplicated hits). Supplements with substantial S0 content:

| supplement | golden-genome mentions |
|---|---|
| Makarova 2011 `all_PD2_islands.txt` (defense islands) | MG1655 328, H37Rv 974, *M. acetivorans* 932, *S. acidocaldarius* 741, *B. theta* 578, SF370 161, *B. subtilis* 110, *M. genitalium* 10 |
| Shmakov 2018 `Islands.tsv` (CRISPR-linked islands) | *S. acidocaldarius* 5,087, MG1655 657, H37Rv 399, *M. acetivorans* 267, *H. volcanii* 103 |
| Wolf 2012 `ar14.fa.gz` + `ar14.arCOG.csv` (arCOG assignments) | about 4,000–4,500 proteins each for *M. acetivorans* and *H. volcanii*; *S. acidocaldarius* also present |
| Weinberg 2010 S4 `.sto` (RNA motif members) | *B. subtilis* 65, *B. theta* 32, SF370 10, MG1655 6, H37Rv 5 |
| Shakya 2017 `pogs.txt` (phage orthologous groups) | T4 2,629, λ 2,343 |
| Shmakov 2019 Set2 `Vicinity.tsv` | *H. volcanii* 96, and *H. volcanii* seeds in `Seeds.tsv` |
| Atkinson 2011 RSH table | one or two entries per golden genome (RelA/SpoT/Rel) |
| Makarova 2014 CARF DataSheet1 | *S. acidocaldarius* 17 |
| Axen 2014 BMC Dataset S1 | MG1655 3 |
| Bell 2026 YprA | *S. acidocaldarius* 7, *H. volcanii* 4, *M. acetivorans* 4 |

Several OPEN hypotheses are anchored in S0 genomes:

| hypothesis | S0 genome | what it concerns |
|---|---|---|
| SYN-001 | *M. genitalium* | MG252 |
| SYN-143 | H37Rv | the Rv2191 fusion |
| SYN-125 | *M. acetivorans* | MA2710–2713 |
| SYN-024..028 | Streptococcus (S0 has SF370) | |

These are natural S0 smoke tests.

### 5.3 Ranked replay list

Ranked by usefulness as an acceptance test: a clear machine-readable hit list, public NCBI accessions, diversity of question types, and an early stage.

**Stage column.** It gives the first stage that holds the relevant taxa in bulk. S0 means the list has entries in golden genomes, usable as a smoke test.

**Size column.** File size on disk. Row counts come from the openpyxl or line counts in the scratch scripts.

**Paths.** Paths are relative to `bioinformatics-synteny/lit/supp/`. Hypotheses are those from the same paper (`select paper_id, group_concat(hypothesis_id) from hypotheses group by paper_id`).

| # | paper (hypotheses) | supplement file(s) | contents | size | codes exercised | IDs | stages |
|---|---|---|---|---|---|---|---|
| 1 | Wu 2026 Druantia III / Zorya II (SYN-262) | `2026_Wu…Druantia…/media-2.xlsx` | Table S1: 13,772 rows of Druantia proteins with RefSeq assembly, contig, WP_ ID, PADLOC HMM name and coordinates. Table S2: 16,358 rows of neighborhood genes with COG/Pfam/KEGG IDs, the defense-enrichment set. | 3.4 MB | Q-SYS (PADLOC) → Q-NBR enrichment → Q-ARCH (ZorE extension) | GCF/WP (public) | S5; S6 for full recall |
| 2 | Shmakov 2019 functionally linked genes (no hypothesis rows) | `2019_Shmakov_Faure…/Set{1,2,3}.ProtocolFiles.tar.gz`, `Set*.Database.tar.gz` | For each of 3 sets: `Seeds.tsv` (36–124 seed loci with contig and coordinates), `Vicinity.tsv/.faa` (523–2,197 neighbor proteins), `Relevance.tsv` (cluster relevance scores), and the BLAST DB used (93 MB). This is the full seed → neighborhood → enrichment protocol. | 0.1–0.4 MB protocol + 25–93 MB DB per set | Q-NBR → Q-CLUST → Q-PROF (relevance) | GenBank contig + coordinates; protein GIs | S0 (*H. volcanii* in Set2) → S5 |
| 3 | Makarova 2011 defense islands (SYN-162, 163) | `2011_Makarova_Koonin_Defense_Islands/Makarova_DI_add.tar.gz` | `all_PD2_islands.txt` (224,026 lines of defense islands as gene runs with GI, coordinates and COG/cluster IDs), `all_Mob_islands.txt` (323,017 lines), `mk.xls.csv` (176,046 lines), Readme | 12 MB packed | Q-FAM (COG) → Q-NBR islands → Q-PROF | GI + NC_ (map by sequence) | **S0** (all golden bacteria and archaea present) → S5 |
| 4 | Weinberg 2010 104 motifs (SYN-279 crcB, 280 pfl, 281 msiK) | `2010_Weinberg…/gb-2010-11-3-r31-S4.GZ` (121 `.sto`), `S5.GZ` (cleaned), `S10.GZ` (pipeline code) | One Stockholm alignment per motif; members carry nt accessions and coordinates (3,311 distinct). Also the authors' IGR/genome-context pipeline. | 2.2 MB packed / 23 MB unpacked | Q-RNA (cmbuild → cmsearch on IGRs) → Q-NBR (downstream gene) | RefSeq v25 NC_ | **S0** (*B. subtilis*, *B. theta*, SF370) → S5 |
| 5 | Weinberg 2017 224 motifs (SYN-284, 285, 286, 340) | `2017_Weinberg…/gkx699_supp.zip` → `File009.zip` (224 motif `.sto`), `File010.zip` (224 annotated `.sto`), `File013.zip` (52 bSCRAP `.sto`), `SupplFile-6-Excel-updated.xlsx` | Includes `drum`, `DUF1646`, `DUF3800-I…XI`, `ldcC`, `nhaA-I/II`, `skipping-rope` for the four hypotheses. Some motifs are from metagenomes (cow-rumen-1…4). | 44 MB zip | Q-RNA → Q-NBR → Q-PROF (DUF3800 flanking) | RefSeq + metagenome IDs | S3/S5; metagenome motifs out of ladder until S7 |
| 6 | Yutin 2024 DNA polymerase swapping in Caudoviricetes (SYN-258) | `2024_Yutin…/12985_2024_2482_MOESM1_ESM.xlsx` | `genome_tree` (18,383 phage genomes with GenBank IDs and genus), `polymerases` (23,246 proteins with PolA/PolB/PolC type and cluster), `4_clades` | 3.6 MB | Q-FAM (PolA/B/C census) → Q-PROF (mutual exclusivity / swapping) | GenBank nt + protein | **S1** (INPHARED); T4/λ at S0 |
| 7 | Faure 2023 Tn7 target selectors (SYN-247, 251) | `2023_Faure…/mmc2.zip` | `LOCIannot.lst` (32,060 lines of annotated loci with assembly and contig), `TNSC.csv` (about 6,384 TnsC with sequence), `TNIQ.csv`, `TNSF_vin.csv`, aligned FASTA (TnsC 6,384; TniQ 4,916), 3 trees | 6.8 MB zip | Q-ITER (TnsC/TniQ expansion) → Q-NBR (target-selector typing) → Q-ARCH | GCA + contig | S2 (plasmids) / S5; SF370 mentions at S0 |
| 8 | Shmakov 2018 CRISPR gene neighborhoods (SYN-343) | `2018_Shmakov…/Islands.tar.gz` (`Islands.tsv`, 848,314 lines, 115 MB), `{Cas10,CysH,STAND,RT,CorA,CARF}.tar.gz` (FASTA + FastTree), `NewProfiles.tar.gz` (128 `icity*.sr` profile alignments), `icity-profiles.xlsx` | CRISPR-linked islands with per-gene COG/cd/cluster IDs and island subtype; protein families with trees; new "CRISPRicity" profiles | 25 MB packed | Q-SYS (CRISPR subtype) → Q-NBR → Q-PROF (CRISPRicity) → Q-ITER (new profiles) | GI + GenBank contig | **S0** (*S. acidocaldarius*, MG1655, H37Rv) → S5 |
| 9 | Bell 2026 YprA helicases (SYN-364) | `2026_Bell…YprA…/media-1.xlsx`, `media-2.xlsx` | 51,793 rows of YprA-like helicases (type ARMADA/BRIGADE/TALON…, class, species, genome accession, locus tag, protein ID); annotation statistics | 3.1 MB | Q-FAM/Q-ITER census → Q-NBR (system typing) → Q-MOTIF (nuclease-site erosion) | GenBank nt + protein | S5; small S0 subset (archaea) |
| 10 | Kaur/Iyer 2022 death-domain / cyclic-nucleotide conflict, eLife 70394 (SYN-273, 355–358) | `2022_Iyer…/elife-70394-fig2-data2.zip` (14 `.aln`: nSTAND1 540 seqs, nSTAND3 638, GT4 375, …), `10…21.*.txt` (neighborhoods and architectures per family: AP-GTPase 1,416 proteins, GREAB-PIN 1,394, nSTAND1 552, SPRY 368, …), profile-profile p-value CSVs | Seed MSAs plus hit neighborhoods; architectures as domain strings | 17.6 MB | Q-REMOTE (profile-profile) → Q-ARCH → Q-NBR | WP/GCA (public) | S5/S6 |
| 11 | Makarova 2014 CARF/WYL (SYN-172, 173) | `2014_Makarova…CARF…/DataSheet1.XLSX` (1,441 CARF proteins), `DataSheet8.XLSX` (5,354 neighbor rows with CDD and Cas labels), `DataSheet2.XLSX` (HHpred), `DataSheet3/6/7.ZIP` (alignments) | Iterative (PSI-BLAST-era) family census with neighborhoods and HHpred | 9.8 MB | Q-ITER → Q-REMOTE → Q-NBR | GenBank protein IDs | S0 (*S. acidocaldarius* 17) → S5 |
| 12 | Jimmy 2020 toxSAS (SYN-106, 107, 346) | `2020_Jimmy…/pnas.1916617117.sd01.xls` (28,681 protein and 15,764 assembly accessions), `sd04.xlsx` (per-subfamily neighbor clusters: RelV 1,874 rows, RelP 1,561, CapRel 1,535, …), `sd05.txt` (alignments), `sd02.xls` (24,049 assemblies) | RSH/SAS census + TA-pair neighborhoods | 26 MB + 7 MB | Q-FAM → Q-NBR (antitoxin partner) → Q-PROF | WP/GCF | S5/S6; S0 small |
| 13 | Axen 2014 BMC loci (SYN-040..048, 331, 332) + Sutter 2021 BMC survey (SYN-049..051, 351, 352) | `2014_Axen…/pcbi.1003898.s001.xlsx` (16,984 rows: locus type, scaffold, GI, Pfam), `s004.hmm` (BMC HMM), `s003/s006.fasta` (MSAs); `2021_Sutter…/MOESM7_ESM.zip` (1,220 + 104 HMMs), `MOESM6` (6 phyloXML trees), `MOESM8` (HMM name table) | A custom-HMM library plus a locus-type truth table | 2.8 MB + 84 MB | Q-FAM with **custom HMMs** → Q-NBR locus typing → Q-CLUST | GI + NC_ (Axen) | S0 (MG1655 eut/pdu loci) → S5 |
| 14 | Bhoobalan-Chitty 2024 archaeal-virus anti-defense (SYN-253, 361) | `2024_Bhoobalan…/MOESM8_ESM.xlsx` (Loci: 48,689 rows Sulfolobales neighborhoods with arCOG; Summary 116 clusters), `MOESM4_ESM.xlsx` (256 ADG predictions in archaeal viruses with protein accessions and sequences; motif-scoring sheet) | Virus–host guilt by association plus a promoter motif (MOESM5) | 4.6 + 0.07 MB | Q-NBR → Q-MOTIF (promoter/HTH motif) → Q-NT | GCF/WP + YP_ | **S0** (*S. acidocaldarius*); viruses at S4 (IMG/VR) or GenBank |
| 15 | Makarova 2024 ESCRT (SYN-256, 362, 363) | `2024_Makarova…ESCRT/Asgard_ESCRT_loci.fasta` (12,399 proteins), `mbio.00335-24-s0009.xlsx` (6,024 + 13,635 loci rows with arCOG), 2 IQ-TREE trees, `s0010.xlsx` (HHpred evidence) | Archaeal (Asgard) locus census | 4.7 + 1.5 MB | Q-FAM (arCOG) → Q-NBR → Q-REMOTE | GCA + IMG (mixed) | S5 (Asgard MAGs) |
| 16 | Atkinson 2011 RSH superfamily (SYN-102..104, 324) | `2011_Atkinson…/pone.0023479.s005.xlsx` | 2,078 RSH proteins (domain, group, organism, subfamily, `gi|…|ref|YP_…|`) | 135 KB | Q-FAM census (classic) → Q-ARCH (subfamily domain composition) | gi + RefSeq protein | **S0** (one or two per golden genome) → S5 |
| 17 | Silas 2025 injected anti-restriction proteins (no hypothesis rows) | `2025_Silas…/41467_2025_63056_MOESM8_ESM.xlsx` | Loci 35,318 rows with GCA, contig, profile hits; Clusters 924; Selected 16 with HHpred | 3.0 MB | Q-ITER → Q-NBR (prophage) → Q-REMOTE | GCA/WP | S1 (phages) / S5 (prophages); λ at S0 |
| 18 | Peters 2017 Tn7-CRISPR (no hypothesis rows) | `2017_Peters…/Supplemental_Table_S2.xlsx` (278,810 rows of TniQ/TnsA loci), `Automatically_identified_spacers.fna` (804 spacers) | Locus table plus spacers: the only clean **spacer → protospacer** replay | 18.7 MB + 69 KB | Q-NT (spacer matching) → Q-NBR | GenBank nt | S1/S2 for protospacer targets; S5/S6 for loci |
| 19 | Gussow 2020 anti-CRISPR ML (SYN-239, 350) | `2020_Gussow…/MOESM4_ESM.txt` (Acr candidates with directon, assembly, CRISPR-Cas context), `MOESM5_ESM.txt` (66 MB cluster table), `MOESM8_ESM.txt` (3,718 proteins) | Guilt by association with HTH (Aca) neighbors | 69 MB | Q-NBR → Q-CLUST → Q-PROF (self-targeting) | WP + GenBank nt | S5/S6 (prophages) |
| 20 | Mestre 2020 defense-associated RTs (SYN-296, 297, 348) | `2020_Mestre…/gkaa1149_supplemental_files.zip` | 6 Newick trees (RT 1,927 leaves with WP IDs; ATPase, Toprim, DUF3800, NDT, PRTase) plus Tables S1–S3 (8 MB) | 56 MB | Q-ITER → Q-NBR → Q-ARCH | WP | S5 |
| 21 | Chaban 2024 thermophilic-virus RNAPs (SYN-254) | `2024_Chaban…/41467_2023_44630_MOESM7_ESM.xlsx` | 3,230 rows of viral proteins with in-house profile ID, E-value and **"potential catalytic motif"** | 63 KB | Q-FAM (custom profile) → Q-MOTIF | GenBank protein/genome | S1/S4 |
| 22 | Corbino 2005 (SYN-112 metA, 113 serC/speF, 114 ybhL) | `2005_Corbino…/gb-2005-6-8-r70-S2…S6.sto` | 5 small RNA alignments (6–70 nt accessions each) | 41 KB | Q-RNA → Q-NBR | RefSeq NC_ | S5 (alpha-proteobacteria) |
| 23 | Makarova 2019 archaeal dark matter (SYN-349) | `2019_Makarova…Dark_Matter/dark-suppl.cls.csv` | 56,366 genes in "dark" clusters with locus ID, GCF/GCA, protein, coordinates, cluster, arCOG; `gen.tab` (524 genomes) | 5 MB | Q-CLUST → Q-NBR (dark islands) | GCF/GCA/WP | S5 (archaea); *H. volcanii* 154 rows |
| 24 | Wolf 2012 arCOGs (no open hypotheses) | `2012_Wolf…/ar14.fa.gz` (394,110 proteins), `ar14.arCOG.csv` (arCOG assignment per protein) | A fixed-model assignment truth set | 75 MB + 36 MB | Q-FAM concordance (our arCOG hmmsearch vs theirs) | GI (map by sequence) | **S0** (*M. acetivorans*, *H. volcanii*, *S. acidocaldarius*) |
| 25 | Ofir 2018 DISARM (SYN-292, 341) | `2018_Ofir…/Supplementary_Table_2.xlsx` (Class 1: 280, Class 2: 44, unclear: 31 systems with per-gene IMG OIDs, coordinates, gene distance), `Table_1.xlsx` (1,371 DUF1998 genes) | A clean system-level hit list, but IMG-only IDs | 0.15 MB | Q-NBR → Q-SYS | IMG OIDs (needs IMG mapping) | S5/S6 after ID mapping |

**Also usable, with lower priority:**

- Vassallo 2022: 23 systems with contig and CDS accessions, plus an HHpred table of 515 rows. It is a Q-REMOTE concordance check (SYN-302).
- Pinilla-Redondo 2020: Aca/Acr homolog tables with WP accessions (SYN-241).
- Makarova 2009 TA and pAgo: `.ali` files and `.xls` (SYN-156..158, 320).
- Abby 2016 TXSScan: `s2.xls` with 4,811 nt accessions, a Q-SYS replay (SYN-095..097).
- Rogozin 2002: `neighbor.txt`, COG arrays (SYN-150, 311, 312).
- Shakya 2017: phage orthologous groups, 52k proteins, S1 (SYN-191).
- Wolf 1999: aaRS seed alignments (SYN-144, 307, 308).
- Yutin 2018 crAss: 2,684 proteins and metagenome contig sources, S4.

### 5.4 Suggested acceptance set (one per capability)

| capability under test | replay | stage | pass metric |
|---|---|---|---|
| Fixed-model census and assignment concordance | #24 Wolf 2012 arCOG; #16 Atkinson RSH | S0 then S5 | arCOG agreement per protein; RSH recall on rep genomes |
| Rule-based systems + neighborhood enrichment | #1 Wu 2026 Druantia/Zorya | S5 | recall of Table S1 loci in reps; overlap of enriched neighbor families with Table S2 |
| Seed → neighborhood protocol | #2 Shmakov 2019 Set2 | S0 (*H. volcanii*) then S5 | Vicinity set recall given `Seeds.tsv` |
| Island detection | #3 Makarova 2011 PD2 islands | S0 | island boundary and gene-membership overlap on the golden genomes |
| Structured RNA from published Stockholm | #4 Weinberg 2010 (msiK/crcB/pfl), #22 Corbino 2005 | S0 then S5 | cmsearch recall of alignment members (by coordinates); downstream-gene concordance |
| Custom HMM library | #13 Sutter 2021 + Axen 2014 | S0 (MG1655 eut/pdu) then S5 | locus-type recall against Dataset S1 |
| Iterative expansion | #7 Faure 2023 TnsC/TniQ; #11 CARF | S2/S5 | recall of `TNSC.csv` members by round; drift check |
| Phage-scale census | #6 Yutin 2024 polymerases | S1 | recall of 23,246 polymerases; PolA/B/C typing concordance |
| Nucleotide (spacer) search | #18 Peters 2017 spacers | S1/S2 | protospacer hits vs published targets |
| Residue/motif filter | #21 Chaban 2024 catalytic-motif column; #9 YprA | S1 / S5 | per-hit motif call agreement |
| Remote homology / profile-profile | #10 Kaur/Iyer 2022 p-value tables; Vassallo HHpred table | local subset | rank concordance of hits |

---

## 6. Query patterns

All numbers are from the TSV unless marked as an estimate.

- **Seed count.** This is the number of distinct protein or RNA families to search as seeds.
  - The median is **3** and the mean 5.3.
  - Distribution: 1 seed in 41 hypotheses, 2 in 90, 3 in 85, 4 in 59, 5 in 30, 6–10 in 43, and more than 10 in 14.
  - 84% (305/365) need 1–5 seeds.
  - The large outliers are whole-catalogue scans: SYN-282 (67 SCRAP families) and three rows of 100 or more.
  - Implication: a typical hypothesis is a handful of profile searches followed by joins, not one search.
- **Iteration.**
  - Explicit iteration (jackhmmer, CM rebuild, or a profile grown from one protein): **7/365 (2%)**. Q-ITER appears in 9 chains.
  - Counting hypotheses that need a custom profile built from seeds or first-round hits, which implies at least one build-search-rebuild cycle: **91/365 (25%)**.
  - The specs almost never name an iterative tool: jackhmmer 1, PSI-BLAST 0 (SQL). They describe outcomes ("all homologs", "diverged paralogs"), so iteration is an unstated requirement the infrastructure must supply.
- **Neighborhoods.**
  - **348/365 (95%)** need a neighborhood step, and Q-NBR appears in 324 chains (89%).
  - Only 17 hypotheses are pure architecture, profile or RNA tests with no locus context.
- **Window size.** For the 348 that need a neighborhood (classified from TSV `window`):

  | window type | n |
  |---|---:|
  | operon / locus / island stated qualitatively | 124 |
  | adjacent or immediately up/downstream | 101 |
  | explicit gene count | 41 |
  | explicit kb | 30 |
  | explicit bp (IGR or adjacency) | 18 |
  | unstated | 16 |
  | genome-wide co-occurrence | 2 |
  | specific cassettes (two- to five-gene) | about 16 |

  Explicit gene windows have a median of **±5 genes** (range 2–20). Explicit kb windows have a median of **±10 kb** (range 2–50; 15–50 kb are defense or transposon islands).

  Default the neighborhood API to ±10 genes / ±10 kb, with operon (directon) boundaries and the IGR sequence returned. About 95% of stated windows fit inside ±10 genes or ±20 kb (estimate).
- **Taxon scope.** Breadth dominates: all_bacteria 100, all_prokaryotes 78 and all_archaea 19, so 197/365 (54%) want a domain-wide or wider scan.
  - Clade-restricted: 92; mixed: 51; single genus: 6.
  - Phage/virus, plasmid/MGE and metagenome together: 19 (5%).
  - Stated genome counts, when present (about 30 hypotheses, estimate), sit mostly at 10^3–10^4 genomes, with a few at 5×10^4–5×10^5 (§3.4).
  - So most hypotheses are answerable at S5 (143,614 species representatives). The clade-restricted ones are first answerable at S3 for small-genome clades, or at S5.
  - The 19 phage, plasmid and metagenome hypotheses need S1, S2, S4 or later.
  - Many specs say "GTDB r220" (90), so the stage index must map r220 taxon names to R226.
- **Nucleotide and IGR are not marginal.** 120 hypotheses (33%) need nucleotide-level work and 104 (28%) need IGR sequence, including 47 that need IGRs for promoter or operator motifs rather than structured RNA. A protein-only search layer would leave about a third of hypotheses unservable.
- **Residue checks and promiscuity are the main curation burdens.** 131 hypotheses (36%) hinge on catalytic residues or motifs, often their *absence*. 143 (39%) involve a promiscuous domain family. Both are per-hit filters on top of a family census, so they belong in the hit-table schema: per-position residue calls and architecture strings. They should not be separate tools.
- **Structure is downstream.**
  - AlphaFold3 or AF-Multimer modelling is requested by 220–225 hypotheses, almost always as the last step.
  - Structure *search* (Foldseek) is used by 23 and profile-profile search (HHpred) by 12.
  - Budget GPU modelling as a per-hypothesis queue of a few complexes (estimate: 1–10 per hypothesis). It is not a corpus-scale structure index.
- **The workload is a join engine more than a search engine.** The typical chain is Q-FAM→Q-NBR followed by a profile, motif or nucleotide step. So the binding constraints are:
  - how fast hits are joined to neighborhoods and IGRs;
  - how fast co-occurrence is computed over 10^5 genomes;
  - how well the search layer brings in user-built profiles and new hits, with provenance, so that the tapes and pair statistics can be recomputed for them.

  Raw search throughput is not the binding constraint.
