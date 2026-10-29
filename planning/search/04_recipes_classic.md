# Classic search recipes, 1993–2016

Status: literature extraction for the search-infrastructure round (2026-10-29). Companion table: `recipes_1993_2016.tsv` (206 rows, schema from `00_brief.md`). Corpus: every `lit/` stem in `../bioinformatics-synteny/lit/` whose year is 1993–2016 (123 stems, `.new.xml` variants folded into their stem).

Evidence is cited as `lit/<file>:<line>`, relative to `bioinformatics-synteny/`. For PDF-derived `.txt` files with short lines, a cited line can be the first line of a window of up to four lines that holds the statement. Every quoted fragment in the TSV was machine-checked against the cited file. Where a paper did not state a value, the table says `not stated`.

## 1. Corpus condition (read this before using the slice)

About a quarter of the slice does not contain the paper its stem names. This affects which papers can be replay targets.

| Problem | Stems |
|---|---|
| Full text is an unrelated paper (no usable methods) | 1996_Lawrence (B. subtilis SpoIIE paper), 2001_Wolf_Rogozin directons (Xanthomonas TALomes; the supplement `Operon.doc` and `genomes.html` is the real directon paper's), 2002_Anantharaman RNA metabolism (Desulfobulbales sulfur; the supplement is the small-molecule-binding-domain tables), 2002_Leipe (2016 EF-G/EF4/BipA review), 2009_Hanson (minimal-genome and curli texts spliced together), 2009_deCrecyLagard tRNA (OBO ontology survey), 2011_Anantharaman macrolide (cardiac surgery), 2011_Breaker (herpesvirus miRNAs), 2012_Haft exosortase (leukemia kinases), 2016_Hanson (2021 bioRxiv, outside the slice), 2016_Pawluk Aca (lake ecosystem) |
| Duplicates another stem | 2006_Iyer_Lake = 2006_Iyer_Burroughs; 2006_Haft PEP-CTERM BMCMicrobiol = 2006_Haft PPP; 2009_Meyer = 2009_Weinberg metagenomes; 2010_Makarova_Yutin Cdv = 2010_Makarova_Yutin vesicle; 2013_Leyn = 2016_Leyn |
| Abstract only | 2004_Barrick, 2006_Pukatzki, 2007_Hebbeln, 2016_Hallberg (and 2010_Aoki `.txt`; its `.xml` is full) |
| A related paper by the same group, not the named one (used as found) | 2005_Selengut holds Haft et al. 2013 (TIGRFAMs and Genome Properties in 2013); 2005_Price NAR holds Price et al. 2006 (operons and the leading strand); 2011_Makarova CRISPR `.txt` is the 2011 Nat Rev Microbiol opinion, and its `.xml` is the 2011 Biol Direct Cas-unification paper (methods taken from the `.xml`); 2015_Makarova CRISPR supplement and `.xml` are Chylinski 2014 (type II); 2012_Yu holds a phosalacine/bialaphos cluster paper; 2013_Jorda holds a 2016 nanocage design paper; 2013_Russell holds Russell 2014 (Bacteroidetes T6SS); 2014_Salomon holds Calder 2014 (VPA1380); 2014_Khoroshkin holds the 2016 paper; 2015_Grell holds the 2018 SkfB structure; 2015_Wetmore `.txt` is a hypochlorite paper (the `.xml` is correct); 2016_Thiaville holds Zhang et al. 2016 PNAS (DUF1537 kinases, built on SSN/GNN) |

So 90 of the 123 stems contain usable search methods, and 33 got a single `Q-OTHER` row that records the defect. Recommendation for the corpus owner: re-fetch the 12 wrong texts. Directons (Wolf 2001), defense-island precursors, anti-CRISPR (Pawluk 2016), the selfish operon (Lawrence and Roth 1996) and P-loop GTPases (Leipe 2002) are all influential papers that are missing.

## 2. Counts per question type

The counts exclude the 33 defect rows. A paper can contribute several rows.

| Code | Rows | Papers | What it was in this era |
|---|---:|---:|---|
| Q-NBR | 38 | 36 | gene neighborhoods, conserved pairs and directons, defense islands, operon prediction, anchor-then-neighbors |
| Q-ITER | 37 | 37 | PSI-BLAST to convergence (dominant); jackhmmer from about 2012; iterative HMM rebuilds (TIGRFAM, RSH) |
| Q-PROF | 19 | 19 | phyletic patterns (COG/arCOG), hypergeometric and mutual-information profiles, PPP, Jaccard profiles, PCA of patterns |
| Q-NT | 17 | 15 | tblastn on unfinished or unannotated genomes, six-frame IGR scans, CRISPR array finding, spacer BLAST, TFBS PWM scans |
| Q-RNA | 11 | 9 | IGR clustering, then CMfinder, then CM search with refinement; Infernal/RAVENNA homolog expansion |
| Q-ARCH | 10 | 10 | domain architecture by case-by-case PSI-BLAST plus PSSM/HMM libraries; Rosetta Stone fusion; fused vs unfused models |
| Q-REMOTE | 10 | 10 | HHpred against PDB, Pfam, CDD and COG libraries for fold or family links |
| Q-FAM | 9 | 9 | fixed-model censuses (COG PSSMs by RPS-BLAST, Pfam, TIGRFAM, modified Pfam HMMs) |
| Q-SSN | 6 | 6 | EFI-EST/GNT, SFLD/Pythoscape networks |
| Q-CLUST | 5 | 5 | COG/KOG/arCOG construction, single-linkage within EC numbers, ClusterFinder BGCs |
| Q-SYS | 5 | 5 | SEED subsystems, Genome Properties, Cas PSSM typing, T6SS COG rules, MacSyFinder/TXSScan |
| Q-MOTIF | 5 | 5 | helicase motif patterns, HEPN motif vetting, SIMBAL, PTPS motif classes |
| Q-STRUCT | 1 | 1 | DALI (Iyer 2003); otherwise structure entered only as HHpred-to-PDB or DaliLite confirmation |
| Q-OTHER (real searches) | 3 | 3 | HGT by best-hit taxon (Omelchenko 2003), island conservation (Makarova 2011), RB-TnSeq operon polarity (Wetmore 2015) |

Tool mentions in the tool and iteration columns: PSI-BLAST 48, HHpred 16, BLASTN 10, BLASTCLUST 9, jackhmmer 8, tblastn 7, Infernal 5, CMfinder 4, RAVENNA 4, PWM scans 4, hmmsearch 4, RPS-BLAST 6, PPP 3, MacSyFinder 1. Replayability over all 206 rows: 62 yes, 89 partial, 55 no.

## 3. Database scale by era (as stated in the papers)

| Year | Stated target size | Evidence |
|---|---|---|
| 1999 | aaRS complement "from 12 complete bacterial genomes"; COG system ">20 genomes" | lit/1999_Wolf_Aravind_Grishin_Koonin_Aminoacyl_tRNA_Synthetases_Evolution_GenomeRes.txt:515, :725 |
| 2002 | 31 genomes, 2890 COGs | lit/2002_Rogozin_Makarova_Wolf_Koonin_Connected_Gene_Neighborhoods_NAR.txt:34 |
| 2003 | 66 genomes, 138,458 proteins, 4873 COGs; STRING on 43 genomes | lit/2003_Tatusov_Koonin_COG_Database_KOGs.txt:44-45; lit/2003_vonMering_Huynen_Jaeggi_Schmidt_Bork_Snel_STRING_Database_NAR.txt:233 |
| 2004 | Prolinks: 83 organisms, 10 million links | lit/2004_Bowers_Pellegrini_Thompson_Fierro_Yeates_Eisenberg_ProLinks_Database_GenomeBiol.txt:7 |
| 2005 | SEED: 383 organisms, 173 subsystems; BLISS IGR all-vs-all on 116 genomes | lit/2005_Overbeek_Begley_Choudhuri_SEED_Subsystems_Genomic_Neighborhoods.txt:3; lit/2005_Corbino_Barrick_Lim_Breaker_Riboswitch_Candidates_Alpha_Proteobacteria.txt:65 |
| 2006 | PPP over 280 genomes | lit/2006_Haft_Paulsen_Selengut_PPP_PEP_CTERM_Exopolysaccharide_Sorting.txt:33 |
| 2007 | 41 archaeal genomes (arCOGs); IGR plus environmental data of "roughly 3.2 billion nucleotides" | lit/2007_Makarova_Sorokin_Novichkov_Wolf_Koonin_arCOGs_Archaeal_Comparative_Genomics_BiolDirect.txt:78; lit/2007_Weinberg_Barrick_Ruzzo_Breaker_22_Candidate_Structured_RNAs.txt:43 |
| 2008 | 750 complete genomes (Sept 2008); 503 RefSeq proteomes (Oct 2008), 450 after de-duplication | lit/2009_Makarova_Wolf_Koonin_Type2_Toxin_Antitoxin_Mobile_Stress_Islands_BiolDirect.txt:259; lit/2010_Wuichet_Zhulin_Origins_Diversification_Chemotaxis_Signal_Transduction_Operons.txt:62 |
| 2009–10 | 976 replicons (Genome Reviews, Sept 2007); 774 bacterial genomes; 718 KEGG genomes; 1,154 IMG genomes (2010) | lit/2009_Boyer_Fichant_Berthod_Vandenbrouck_Attree_T6SS_Conserved_Gene_Clusters_BMCGenomics.txt:98; lit/2010_Fozo_Makarova_Shabalina_Koonin_Storz_Type_I_sRNA_Toxin_Antitoxin_NAR.txt:1; lit/2010_Omelchenko_Galperin_Wolf_Koonin_Non_Homologous_Isofunctional_Enzymes.txt:1781; lit/2014_Cimermancic_Medema_Claesen_Fischbach_Insights_Secondary_Metabolism_BGC_Atlas_Cell.txt:45 |
| 2010–11 | 1,055 genomes and 3,468,563 proteins (March 2010), 383 representatives; PPP over 1,450 genomes; 1,072 proteomes | lit/2011_Makarova_Koonin_Defense_Islands.txt:177-203; lit/2011_Haft_Mycofactocin_Radical_SAM_Peptide_Cofactor_Operon.txt:93; lit/2011_Atkinson_Tenson_Hauryliuk_RSH_Superfamily_ppGpp_Synthetases_Hydrolases.txt:20 |
| 2012–13 | 1,348 genomes (Ter phyletic PCA); RegPrecise regulons propagated to 640 genomes; RefSeq Feb 2013, 2,262 genomes | lit/2012_Anantharaman_Iyer_Aravind_Ter_Dependent_Stress_Membrane_Synteny.txt:191; lit/2013_Novichkov_Kazakov_Ravcheev_Leyn_Rodionov_RegPrecise_3_Comparative_Regulomics_BMCGenomics.txt:30; lit/2014_Makarova_Koonin_Aravind_CARF_WYL_Domains.txt:579-582 |
| 2014–15 | 2,751 complete genomes (1 Feb 2014); 2,643 RefSeq genomes (Nov 2013); NCBI WGS searched by tblastn | lit/2015_Makarova_Koonin_Updated_Evolutionary_Classification_CRISPR.txt:175-178; lit/2016_Abby_Cury_Guglielmini_Rocha_Protein_Secretion_Systems_TXSScan.txt:704; lit/2015_Shmakov_Koonin_Zhang_Discovery_Class2_CRISPR.txt:751-752 |
| 2015–16 | PubSEED 11,411 genomes (Nov 2015); UniProt 2016_02 ">60 million sequences" | lit/2016_Zallot_Harrison_deCrecyLagard_Functional_Annotations_Missing_Genes_Synteny.txt:86; lit/2016_Thiaville_Kellner_Dedon_IwataReuyl_deCrecyLagard_7_Deazaguanine_DNA_Modification_dpd_Operon.txt:120 |

The classic recipes were tuned on 10^1 to 10^3 genomes. NR was the default protein target and its release was almost never stated. Zhang 2012 is the exception ("NR database from May 23rd 2012"), and only because a referee demanded it (lit/2012_Zhang_Iyer_Aravind_Polymorphic_Toxin_Systems.txt:401-405).

## 4. Canonical recipes as step sequences

### 4.1 Iterative transitive profile search (Aravind/Iyer/Burroughs, Koonin/Makarova)

1. **Seed selection.** Take multiple diverged representatives, not one sequence ("multiple representatives from each of the above mentioned structurally characterized versions", lit/2006_Iyer_Burroughs_Aravind_Prokaryotic_Ubiquitin_Antecedents.txt:35; "several diverged representatives within the family", lit/2011_Makarova_Koonin_CRISPR_Cas_Evolution_Classification.xml:9). Seeds can be single sequences or curated alignments.
2. **Split multidomain queries first.** Define domain boundaries, then search each domain separately: "case-by-case dissection ... with the sequences of individual domains as queries" (lit/1999_Aravind_Walker_Koonin_DNA_Repair_Domains.txt:148-150); aaRS non-core domains "cut out of the alignments" (lit/1999_Wolf_Aravind_Grishin_Koonin_Aminoacyl_tRNA_Synthetases_Evolution_GenomeRes.txt:727); toxin C-terminal domain boundaries defined before seeding (lit/2012_Zhang_Iyer_Aravind_Polymorphic_Toxin_Systems.txt:45); Iyer 2003 split the SGC extension into N- and C-terminal parts after first-round hits (quoted by the referee at lit/2012_Zhang_Iyer_Aravind_Polymorphic_Toxin_Systems.txt:403).
3. **Composition check.** Mask with SEG, or turn composition statistics off when no bias is found: "If no such bias was detected, searches were run with the composition-based statistics turned off, in order to maximize sensitivity" (lit/2003_Iyer_Koonin_Aravind_RNA_Polymerase_DPBB_Evolution.txt:2484-2488). The Koonin lab defaults to no composition statistics and no low-complexity filter (Ago 2009, DI 2011, CARF 2014, Shmakov 2015).
4. **Iterate.** Run PSI-BLAST with inclusion h = 0.01 to convergence (sections 5 and 6 list the variants).
5. **Vet every round.** Newly included sequences are "examined for being false positives and the search was continued with the same e-value threshold only if the profile was uncorrupted" (lit/2012_Zhang_Iyer_Aravind_Polymorphic_Toxin_Systems.txt:353). Sub-threshold hits carrying the conserved motif are added by hand, and false positives are removed (lit/2014_Makarova_Koonin_Aravind_CARF_WYL_Domains.txt:589-593; lit/2006_Makarova_Koonin_Prokaryotic_RNAi_CRISPR_Hypothesis.txt:2104-2106; lit/2013_Anantharaman_Makarova_Koonin_Aravind_HEPN_Superfamily.txt:2072).
6. **Stop.** Stop at convergence, or "until the last iteration before the first known false positives appear" (lit/2009_Makarova_Wolf_Koonin_Type2_Toxin_Antitoxin_Mobile_Stress_Islands_BiolDirect.txt:255), or after a fixed number of rounds (3 by default in 1999, lit/1999_Aravind_Walker_Koonin_DNA_Repair_Domains.txt:147).
7. **Report significance from the round of first inclusion.** Later rounds inflate e-values for sequences already in the profile (lit/1999_Wolf_Aravind_Grishin_Koonin_Aminoacyl_tRNA_Synthetases_Evolution_GenomeRes.txt:727; lit/1999_Aravind_Walker_Koonin_DNA_Repair_Domains.txt:141). Papers report links as "e = 10-3, iteration 2", and some links appear only late ("e = 0.01, iteration 11", lit/2006_Iyer_Burroughs_Aravind_Prokaryotic_Ubiquitin_Antecedents.txt:37; "iteration 13", lit/2012_Anantharaman_Iyer_Aravind_Ter_Dependent_Stress_Membrane_Synteny.txt:125).
8. **Go transitive.** Re-seed with new hits: "Each retrieved sequence was used as the query for additional searches until no new sequences could be detected" (lit/2006_Makarova_Koonin_Prokaryotic_RNAi_CRISPR_Hypothesis.txt:2102-2104). Re-seed with borderline motif-bearing hits ("two-pronged search strategy", lit/2013_Anantharaman_Makarova_Koonin_Aravind_HEPN_Superfamily.txt:204-205). Re-seed at library level ("integrated into the initial library ... iterated this procedure for several rounds", lit/2016_Zhang_Burroughs_Vidal_Aravind_Transposons_To_Toxins_Polymorphic_Conflict_Systems_NAR.txt:20). Two quality controls ride on this step:
   - Reverse searches confirm links ("Reverse searches with the C-terminal region of cytochrome F", lit/2003_Iyer_Koonin_Aravind_RNA_Polymerase_DPBB_Evolution.txt:2116).
   - Recovery order counts as evidence ("prior to other members of the GHKL superfamily", lit/2008_Iyer_Abhiman_Aravind_MutatRes_Prokaryotic_DNA_Repair_Modification_Neighborhoods.txt:9).
9. **Confirm remote links.** Use HHpred against PDB, Pfam and CDD libraries, secondary-structure prediction (PHD, then JPred), and DALI/DaliLite. Report HHpred probability or p-value (lit/2012_Anantharaman_Iyer_Aravind_Ter_Dependent_Stress_Membrane_Synteny.txt:55; lit/2015_Krupovic_Koonin_Viral_Nucleocapsid_From_Cas4_Nuclease.txt:178-211).
10. **Cluster and align.** Run BLASTCLUST "for both classification and culling of nearly identical sequences" (lit/2013_Iyer_Zhang_Burroughs_Aravind_Phage_DNA_Hypermodification_Synteny.txt:15), align with Kalign, MUSCLE or PCMA plus manual correction, and build an HMM or PSSM for the next campaign. Pfam is "used as a guide" and augmented with newly found divergent members (lit/2015_Burroughs_Zhang_Schaeffer_Aravind_Nucleotide_NAD_Signaling_Networks.txt:21).

### 4.2 Domain-lego architecture mapping

1. For every protein recovered in 4.1, determine the complete architecture "using a combination of case-by-case PSI-BLAST searches and searches against libraries of position specific score matrices (PSSMs) or HMMs of previously characterized protein domains" (lit/2006_Iyer_Burroughs_Aravind_Prokaryotic_Ubiquitin_Antecedents.txt:49).
2. Find regions not covered by Pfam, cluster them, align, and search again to define new domains, building "a comprehensive library of domain architectures covering the entire protein space of the identified systems" (lit/2015_Burroughs_Zhang_Schaeffer_Aravind_Nucleotide_NAD_Signaling_Networks.txt:21). The 2016 Crn work isolated such regions and ran PSI-BLAST and JACKHMMER on them (lit/2016_Zhang_Burroughs_Vidal_Aravind_Transposons_To_Toxins_Polymorphic_Conflict_Systems_NAR.txt:33).
3. Add TM and signal-peptide calls (TMHMM, Phobius).
4. Merge architecture and neighborhood co-occurrence into one directed domain network, with edges N to C for fusions and 3' to 5' for neighbors, and edge weight equal to species-level counts. Render with Cytoscape (Kamada–Kawai) (lit/2015_Burroughs_Zhang_Schaeffer_Aravind_Nucleotide_NAD_Signaling_Networks.txt:21).
5. Use the resulting templates (for example, trafficking domain plus pre-toxin plus C-terminal toxin, with an immunity gene next door) to recognize new members (lit/2012_Zhang_Iyer_Aravind_Polymorphic_Toxin_Systems.txt:353-355).
6. Handle promiscuous domains in one of three ways:
   - Mask them before clustering. KOGs masked WD40, ANK, TPR, LRR, IG and others with RPS-BLAST + CDD (lit/2003_Tatusov_Koonin_COG_Database_KOGs.txt:1271-1289).
   - Choose queries that avoid them. Kapitonov used the RuvC-only family "to avoid multiple false-positive hits to numerous HNH nucleases" (lit/2015_Kapitonov_Makarova_Koonin_ISC_Transposons_Cas9_Homologs.txt:138).
   - Drop them from statistics. The DI positive set lost a family "prone to produce false positives due to the presence of a ubiquitous helix-turn-helix (HTH) domain" (lit/2011_Makarova_Koonin_Defense_Islands.txt:240-243).

### 4.3 Rosetta Stone / fusion inference

1. **Precomputed variant (Prolinks).** BLAST every protein of a genome against nr. Call a fusion link when "two nonhomologous proteins both align over at least 70% of their sequence to different portions of a third protein" (lit/2004_Bowers_Pellegrini_Thompson_Fierro_Yeates_Eisenberg_ProLinks_Database_GenomeBiol.txt:47).
2. **Promiscuity correction.** Compute a hypergeometric P from k' Rosetta proteins, homolog counts n and m, and nr size N, so that proteins with many homologs ("kinase or zinc finger domains") are not linked by chance (lit/2004_Bowers_Pellegrini_Thompson_Fierro_Yeates_Eisenberg_ProLinks_Database_GenomeBiol.txt:49-53).
3. **Orthology variant (STRING 2003).** "fused proteins are not recognized by homology, but rather by orthology of the fused parts to other, non-fused proteins". The score is non-redundant observations divided by species with either group, calibrated against KEGG co-membership (lit/2003_vonMering_Huynen_Jaeggi_Schmidt_Bork_Snel_STRING_Database_NAR.txt:151-154, :169-175). STRING 2000 also counts fused genes as part of a run (lit/2000_Snel_Lehmann_Bork_Huynen_STRING_Gene_Context_Server.txt:53-54).
4. **Manual variant (Koonin, Aravind).** Fusions are read from architectures in COG context and checked against phyletic patterns and adjacency (lit/2000_Galperin_Koonin_Whos_Your_Neighbor.txt:537, :574). Separate fused-form models are built where needed: TIGR03882 for the cyclodehydratase fused to the docking scaffold (lit/2010_Haft_Basu_Mitchell_TOMM_Thiazole_Oxazole_RiPP_Operons.txt:80).

### 4.4 Neighborhood census (anchor, then neighbors, then statistics)

1. **Anchor set.** The anchor can come from a profile search (4.1), a fixed model (COG PSSM, TIGRFAM, Pfam), or a tblastn hit for unannotated genes.
2. **Window.** Windows varied by lab:
   - ±5 genes (lit/2012_Zhang_Iyer_Aravind_Polymorphic_Toxin_Systems.txt:355; lit/2011_Makarova_Koonin_Defense_Islands.txt:246-247)
   - ±5–10 genes for initial screening (lit/2014_Makarova_Koonin_Aravind_CARF_WYL_Domains.txt:616-617)
   - ±10 genes for enrichment (lit/2009_Makarova_Koonin_Prokaryotic_Argonautes_Defense.txt:1275)
   - within 3 genes (lit/2010_Wuichet_Zhulin_Origins_Diversification_Chemotaxis_Signal_Transduction_Operons.txt:70)
   - ±10 kb, merged (lit/2014_Axen_Erbilgin_Kerfeld_Taxonomy_Bacterial_Microcompartment_Loci.txt:308-318)
   - ±20 kb with re-calling by GeneMarkS (lit/2015_Shmakov_Koonin_Zhang_Discovery_Class2_CRISPR.txt:753-755)
   - 6 kb around RNAs (lit/2015_Weinberg_Kim_Chen_Breaker_New_Classes_Riboswitches_Synteny.txt:75)
3. **Operon and directon rule.** Rules also varied:
   - Same strand with intergenic distance <100 nt (lit/2009_Makarova_Koonin_Prokaryotic_Argonautes_Defense.txt:325-328; lit/2009_Makarova_Wolf_Koonin_Type2_Toxin_Antitoxin_Mobile_Stress_Islands_BiolDirect.txt:255).
   - <200 nt in RegPredict CRONs (lit/2016_Leyn_Kazakov_Novichkov_Rodionov_Transcriptional_Regulons_Proteobacteria_MicrobGenom.txt:3412-3413).
   - Runs on one strand not interrupted by >300 bp (lit/2000_Snel_Lehmann_Bork_Huynen_STRING_Gene_Context_Server.txt:52-53).
   - A per-genome directon spacer limit at "the 95% quantile of the distribution of the actual spacer lengths" (lit/2002_Rogozin_Makarova_Wolf_Koonin_Connected_Gene_Neighborhoods_NAR.txt:74).
   - Head-to-head pairs checked for a shared bidirectional promoter (lit/2012_Zhang_Iyer_Aravind_Polymorphic_Toxin_Systems.txt:355).
4. **Annotate neighbors.**
   - Aravind lab: BLASTCLUST all neighbor proteins, then label clusters by domain architecture (TASS Perl package).
   - Koonin lab: COG PSSMs plus Pfam/CDD for every gene (lit/2011_Makarova_Koonin_Defense_Islands.txt:181-186).
   - Others: Pfam hmmsearch "with a loose e-value cutoff of 0.01" (lit/2014_Axen_Erbilgin_Kerfeld_Taxonomy_Bacterial_Microcompartment_Loci.txt:323-326).
   - Unassigned neighbors go through two-round blastclust (lit/2011_Makarova_Koonin_Defense_Islands.txt:223-232).
5. **Conservation and enrichment statistics.**
   - Conserved pairs in ≥3 genomes with a Monte Carlo null (lit/2002_Rogozin_Makarova_Wolf_Koonin_Connected_Gene_Neighborhoods_NAR.txt:78), or ≥10 genomes as a conservative operon proxy (lit/2003_Omelchenko_Makarova_Wolf_Koonin_Mosaic_Operons_In_Situ_Displacement.txt:3916-3919).
   - Shuffled-genome null for co-occurrence in runs (lit/2000_Snel_Lehmann_Bork_Huynen_STRING_Gene_Context_Server.txt:140-145).
   - 10,000 pseudo-genome shuffles (lit/2008_Iyer_Abhiman_Aravind_MutatRes_Prokaryotic_DNA_Repair_Modification_Neighborhoods.txt:32).
   - Binomial per window combined by a Fisher omnibus test (lit/2009_Makarova_Koonin_Prokaryotic_Argonautes_Defense.txt:1270-1294).
   - DI family Z-score with Bonferroni P < 0.05, islands defined as runs of non-negative genes bounded by negative genes, and positive/negative sets refined over two passes (lit/2011_Makarova_Koonin_Defense_Islands.txt:211-258).
6. **Iterate from context.** Neighbors become new seeds: "Each of the families belonging to the conserved neighborhoods were used as starting points for further PSI-BLAST searches" (lit/2006_Iyer_Burroughs_Aravind_Prokaryotic_Ubiquitin_Antecedents.txt:49; lit/2014_Makarova_Koonin_Aravind_CARF_WYL_Domains.txt:341; lit/2006_Makarova_Koonin_Prokaryotic_RNAi_CRISPR_Hypothesis.txt:202-204).
7. **De-bias sampling.** Pick one representative genome per genus (383 of 1,055, lit/2011_Makarova_Koonin_Defense_Islands.txt:200-208; 659 of 2,262, lit/2014_Makarova_Koonin_Aravind_CARF_WYL_Domains.txt:582-586). Drop near-identical species by ribosomal distance below 0.023 (lit/2010_Wuichet_Zhulin_Origins_Diversification_Chemotaxis_Signal_Transduction_Operons.txt:62). Exclude close strains from PWM training (lit/2013_Ravcheev_Best_Sernova_Kazanov_Novichkov_Rodionov_Polysaccharide_Utilization_Regulons_BMCGenomics.txt:1448-1453).
8. **Rule-based systems.** Systems follow from the census:
   - T6SS: COG PSSM coverage ≥80% and E ≤ 1e-6, with fewer than 5 intervening genes and at least 5 distinct COGs (lit/2009_Boyer_Fichant_Berthod_Vandenbrouck_Attree_T6SS_Conserved_Gene_Clusters_BMCGenomics.txt:98).
   - MacSyFinder: mandatory, accessory and forbidden components, a quorum, and inter_gene_max_space (lit/2016_Abby_Cury_Guglielmini_Rocha_Protein_Secretion_Systems_TXSScan.txt:715-736).
   - Cas typing: strict PSSM hits, then relaxed hits allowed only in the neighborhood (lit/2015_Makarova_Koonin_Updated_Evolutionary_Classification_CRISPR.txt:179-186).

### 4.5 IGR motif discovery (Breaker/Weinberg/Ruzzo)

1. **Choose the IGR set.** Options, in historical order:
   - Gene-anchored: "potential 5' UTRs of genes classified in the Conserved Domain Database", grouped per domain and per phylum group (lit/2007_Weinberg_Barrick_Ruzzo_Breaker_22_Candidate_Structured_RNAs.txt:31-35).
   - Gene-independent: all IGRs "of at least 30 nucleotides" (lit/2010_Weinberg_Wang_Ruzzo_Breaker_104_Candidate_Structured_RNAs.txt:99).
   - Context-targeted: IGRs near ribozyme-associated elements, "~7 million base pairs" (lit/2015_Weinberg_Kim_Chen_Breaker_New_Classes_Riboswitches_Synteny.txt:17).
   In every option, IGRs are extended 50 nt into flanking ORFs to absorb annotation errors (lit/2007_Weinberg_Barrick_Ruzzo_Breaker_22_Candidate_Structured_RNAs.txt:33). Runs are made per taxon or per metagenome source "to reduce false positives" (lit/2010_Weinberg_Wang_Ruzzo_Breaker_104_Candidate_Structured_RNAs.txt:99).
2. **Mask known RNAs.** Use Rfam, tRNAscan-SE, rRNA WU-BLAST, and RAVENNA with published riboswitch alignments (lit/2010_Weinberg_Wang_Ruzzo_Breaker_104_Candidate_Structured_RNAs.txt:97).
3. **Cluster IGRs.** Run BLASTN "-W 7 -G 2 -E 2 -q -2 -m 8". Keep scores between S and 100, merge nodes that overlap by at least 5 nt, and keep transitive clusters with at least 3 members. Scan S from 35 to 85 against known-RNA recovery and take the union (lit/2010_Weinberg_Wang_Ruzzo_Breaker_104_Candidate_Structured_RNAs.txt:101-105).
4. **Discover motifs.** Run CMfinder (v0.2 in 2007) on each cluster.
5. **Expand by homology and refine.** Search with RaveNnA/Infernal CMs: in 2007, "E-value cutoff of 10", global mode, ML-heuristic filters. CMfinder then realigns with the new homologs (lit/2007_Weinberg_Barrick_Ruzzo_Breaker_22_Candidate_Structured_RNAs.txt:33). Later work used the null3 model (lit/2010_Weinberg_Wang_Ruzzo_Breaker_104_Candidate_Structured_RNAs.txt:105) and iterative Infernal 1.0/1.1 rounds (lit/2013_Nelson_Sudarsan_Furukawa_Weinberg_Breaker_ydaO_Cyclic_di_AMP_Riboswitch.txt:77).
6. **Filter automatically.** Reject motifs with "no covarying base-pair positions", G+C below 24%, reverse-complement overlap of 30% or more, or "fewer than six positions that were >=97% conserved". Use GSC weighting, and do not call a pair covarying when more than 10% non-canonical pairs are present (lit/2010_Weinberg_Wang_Ruzzo_Breaker_104_Candidate_Structured_RNAs.txt:105-107; lit/2009_Weinberg_Perreault_Meyer_Breaker_Exceptional_Structured_ncRNAs_Metagenomes.txt:41).
7. **Curate by hand.** Edit in RALEE, check terminators (Rnall), and look up the downstream genes. Use targeted searches to reduce false positives in a large database: subset by clade, search IGRs only, or run tblastn for the associated gene and then search 2 kb upstream (lit/2007_Weinberg_Barrick_Ruzzo_Breaker_22_Candidate_Structured_RNAs.txt:41-43).
8. **Validate.** In-line probing and reporter fusions.

### 4.6 Missing-gene and NOGD discovery (Overbeek/Osterman/Rodionov/de Crécy-Lagard, Haft PPP, Koonin NISE)

1. **Frame the problem.** Build a populated subsystem: functional roles by genomes, with gaps where a role is expected but "can not be located by similarity" (lit/2005_Overbeek_Begley_Choudhuri_SEED_Subsystems_Genomic_Neighborhoods.txt:223).
2. **Physical clustering.** Find candidate genes that cluster on the chromosome with known pathway genes. Weak homologs in the cluster can then be annotated, "Because these were weak homologs they could not be accurately characterized without considering the chromosomal neighborhood" (lit/2005_Overbeek_Begley_Choudhuri_SEED_Subsystems_Genomic_Neighborhoods.txt:85). Later work ran SEED clustering on 58 archaeal genomes (lit/2012_Phillips_Swairjo_IwataReuyl_deCrecyLagard_Archaeosine_Synthase_ArcS_RaSEA_NOGD.txt:48) and 11,411 genomes (lit/2016_Zallot_Harrison_deCrecyLagard_Functional_Annotations_Missing_Genes_Synteny.txt:86).
3. **Phyletic complementarity.** Look for candidates present where the known enzyme is absent: identical or complementary patterns (lit/2010_Makarova_Yutin_Bell_Koonin_Cell_Division_Vesicle_Formation_Archaea_BiolDirect.txt:30; lit/2014_Makarova_Krupovic_Koonin_Evolution_Replicative_DNA_Polymerases_Archaea.txt:799-800) and the COG phyletic-pattern search tool (lit/2003_Tatusov_Koonin_COG_Database_KOGs.txt:899).
4. **Partial phylogenetic profiling (Haft).**
   - The profile is the set of genomes carrying a marker (a TIGRFAM above its trusted cutoff, or more than two PEP-CTERM copies).
   - For each protein in a YES genome, walk its BLAST hit list in order of significance, keeping one hit per species.
   - At each depth, score the binomial probability of the YES count, and keep the best depth as the protein's score and its candidate family cutoff (lit/2006_Haft_Paulsen_Selengut_PPP_PEP_CTERM_Exopolysaccharide_Sorting.txt:117-121; lit/2011_Haft_Mycofactocin_Radical_SAM_Peptide_Cofactor_Operon.txt:91-95).
   - SIMBAL applies the same idea to every subsequence window to localize function-determining residues (lit/2011_Haft_Mycofactocin_Radical_SAM_Peptide_Cofactor_Operon.txt:99).
5. **Regulon support (Rodionov).**
   - Build the training set from upstream regions of genes clustered with the TF, and find the motif with SignalX or RegPredict Discover.
   - Scan genomes with the PWM at a threshold equal to the lowest training-site score (lit/2008_Rodionov_Li_Osterman_NAD_Metabolism_Transcriptional_Regulons_Synteny.txt:20), or 10% below it (lit/2014_Khoroshkin_Leyn_VanSinderen_Rodionov_Transcriptional_Regulation_Carbohydrate_Utilization_Bifidobacteria.txt:459-466).
   - Keep sites that pass the consistency check across orthologous operons (CRONs, iterated to convergence), then refine the PWM and rescan (lit/2016_Leyn_Kazakov_Novichkov_Rodionov_Transcriptional_Regulons_Proteobacteria_MicrobGenom.txt:3412-3423).
6. **Analogous-enzyme census (NISE).** Single-link proteins that share an EC number at E 0.01 and discard clusters that share a SCOP fold (lit/2010_Omelchenko_Galperin_Wolf_Koonin_Non_Homologous_Isofunctional_Enzymes.txt:1798-1815).
7. **SSN/GNN step (EFI, from 2011).** Build an SSN at a chosen alignment score, then a GNN to find conserved neighbors of each isofunctional cluster (lit/2016_Zallot_Harrison_deCrecyLagard_Functional_Annotations_Missing_Genes_Synteny.txt:90; lit/2016_Thiaville_Kellner_Dedon_IwataReuyl_deCrecyLagard_7_Deazaguanine_DNA_Modification_dpd_Operon.txt:303-306).
8. **Validate.** Complementation, enzyme assays, metabolomics.

## 5. Recurring parameter conventions

| Parameter | Convention(s) seen | Example evidence |
|---|---|---|
| PSI-BLAST inclusion | 0.01 is the default in both NCBI labs; 0.005 (Crn 2016); 0.1 for proteins under 100 aa (TA 2009); 1.0 with expect 100, PAM70, word size 2 for very short peptides (type I toxins); 0.5 for a 49-aa leader peptide; 0.001 default and 7e-4 tuned (Zhulin) | lit/2006_Iyer_Burroughs_Aravind_Prokaryotic_Ubiquitin_Antecedents.txt:133; lit/2016_Zhang_Burroughs_Vidal_Aravind_Transposons_To_Toxins_Polymorphic_Conflict_Systems_NAR.txt:20; lit/2009_Makarova_Wolf_Koonin_Type2_Toxin_Antitoxin_Mobile_Stress_Islands_BiolDirect.txt:255; lit/2010_Fozo_Makarova_Shabalina_Koonin_Storz_Type_I_sRNA_Toxin_Antitoxin_NAR.txt:29; lit/2010_Haft_Basu_Mitchell_TOMM_Thiazole_Oxazole_RiPP_Operons.txt:78; lit/2010_Wuichet_Zhulin_Origins_Diversification_Chemotaxis_Signal_Transduction_Operons.txt:66 |
| Rounds | "to convergence" is the norm; reported convergence after 3 to 12 iterations; 6 iterations gave 1,227 Cas9 homologs; fixed 3 rounds (1999 survey); jackhmmer with 5 rounds and a bits-per-residue threshold sweep (Meeske) | lit/2010_Wuichet_Zhulin_Origins_Diversification_Chemotaxis_Signal_Transduction_Operons.txt:66-84; lit/2015_Kapitonov_Makarova_Koonin_ISC_Transposons_Cas9_Homologs.txt:141; lit/2016_Meeske_Sham_Bernhardt_Rudner_SEDS_Peptidoglycan_Polymerases_Synteny.txt:323-327 |
| Statistics flags | composition-based statistics off and low-complexity filter off (Koonin lab); composition statistics on (Zhang 2012, Haft TOMM) | lit/2014_Makarova_Koonin_Aravind_CARF_WYL_Domains.txt:587-589; lit/2012_Zhang_Iyer_Aravind_Polymorphic_Toxin_Systems.txt:403 |
| Fixed database size | "effective database size of 10^8" for BLAST and RPS-BLAST; PSI-BLAST "fixing the database size to 10^8" | lit/2005_Price_Huang_Arkin_Alm_Operon_Formation_Driven_by_Coregulation_GenomeRes.txt:412-415; lit/2015_Spang_Saw_Jorgensen_ZarembaNiedzwiedzka_Ettema_Lokiarchaeota_ESCRT_GTPase_Nature.txt:74 |
| Fixed-model cutoffs | TIGRFAM trusted cutoff; COG PSSM ≥80% coverage with E ≤ 1e-6; HMMER i-evalue < 1e-3 with >50% profile coverage (MacSyFinder); E ≤ 1e-5 "chosen based on manual inspection" (BMC); CDD E 1e-4; RSH HMMs at E-2 for gathering, then E-4 and E-5 per domain | lit/2009_Boyer_Fichant_Berthod_Vandenbrouck_Attree_T6SS_Conserved_Gene_Clusters_BMCGenomics.txt:98; lit/2016_Abby_Cury_Guglielmini_Rocha_Protein_Secretion_Systems_TXSScan.txt:733-734; lit/2014_Axen_Erbilgin_Kerfeld_Taxonomy_Bacterial_Microcompartment_Loci.txt:265-269; lit/2011_Atkinson_Tenson_Hauryliuk_RSH_Superfamily_ppGpp_Synthetases_Hydrolases.txt:20-22 |
| Orthology | bidirectional best hit at E < 0.01 by Smith–Waterman (STRING); BBH with 75% coverage and E < 1e-5 (Tamames, Price); symmetric best-hit triangles (arCOG); COGNITOR | lit/2000_Snel_Lehmann_Bork_Huynen_STRING_Gene_Context_Server.txt:79; lit/2001_Tamames_Evolution_Gene_Order_Conservation_Prokaryotes.txt:78; lit/2007_Makarova_Sorokin_Novichkov_Wolf_Koonin_arCOGs_Archaeal_Comparative_Genomics_BiolDirect.txt:80 |
| Clustering | BLASTCLUST "-L 0.5 -S 1.75" for neighbors and "-L 0.75 -S 1.0" for TA representatives; 90% coverage and 90% identity; two rounds (1.5 bits/position with 0.8 one-way coverage, then 0.7 with 0.7 two-way); 95% identity representatives used as the next seeds; TribeMCL at E 1e-6 with I = 2 | lit/2009_Makarova_Wolf_Koonin_Type2_Toxin_Antitoxin_Mobile_Stress_Islands_BiolDirect.txt:255-259; lit/2015_Kapitonov_Makarova_Koonin_ISC_Transposons_Cas9_Homologs.txt:94-95; lit/2011_Makarova_Koonin_Defense_Islands.txt:227-232; lit/2010_Fozo_Makarova_Shabalina_Koonin_Storz_Type_I_sRNA_Toxin_Antitoxin_NAR.txt:29; lit/2009_Boyer_Fichant_Berthod_Vandenbrouck_Attree_T6SS_Conserved_Gene_Clusters_BMCGenomics.txt:106 |
| Phylogenetic profile | presence when BLAST E < 1e-10, scored by hypergeometric P; mutual information with collapse of concordant taxa; PPP best-depth binomial; Jaccard distance with mean linkage | lit/2004_Bowers_Pellegrini_Thompson_Fierro_Yeates_Eisenberg_ProLinks_Database_GenomeBiol.txt:21-25; lit/2003_vonMering_Huynen_Jaeggi_Schmidt_Bork_Snel_STRING_Database_NAR.txt:157-167; lit/2009_Boyer_Fichant_Berthod_Vandenbrouck_Attree_T6SS_Conserved_Gene_Clusters_BMCGenomics.txt:108 |
| Remote homology | HHpred probability or p-value reported per link, for example "P = 94.5" against CDD, with the hit depending on the library | lit/2015_Krupovic_Koonin_Viral_Nucleocapsid_From_Cas4_Nuclease.txt:178-211 |
| Short and missed genes | tblastn over contigs "because these regions are often not annotated"; six-frame IGR translation (16–70 aa); re-calling with GeneMarkS around anchors; 15–20% wrong start calls observed | lit/2011_Poole_Diner_Aoki_Low_Hayes_Rhs_Proteins_Polymorphic_Toxin_Immunity_Synteny.txt:737-740; lit/2010_Fozo_Makarova_Shabalina_Koonin_Storz_Type_I_sRNA_Toxin_Antitoxin_NAR.txt:33; lit/2016_Zallot_Harrison_deCrecyLagard_Functional_Annotations_Missing_Genes_Synteny.txt:86 |
| CRISPR arrays | 20-mer anchors extended while information content ≥1.5 bits, then BLASTN (word 7, E 0.001), keeping repeats of 20–50 nt with spacers of 15–60 nt; later PILER-CR or CRISPRfinder; nearest-array cut-off of 530 bp; spacers searched by MEGABLAST (word 11, E 0.01) | lit/2006_Makarova_Koonin_Prokaryotic_RNAi_CRISPR_Hypothesis.txt:2117-2129, :2231; lit/2015_Makarova_Koonin_Updated_Evolutionary_Classification_CRISPR.txt:1363-1368 |
| SSN | "alignment score" = -log10 of 2^-bitscore × (query length × subject length), independent of database size; Option A collects up to 5,000 hits at E < 1e-5; the threshold is swept by eye | lit/2015_Gerlt_Bouvier_Davidson_Imker_Sadkhin_Zallot_EFI_EST_SSN_GNT_BBA.txt:34, :86-91 |

## 6. What breaks or changes at 10^5–10^6 genomes

1. **Manual vetting every round.** The defining quality control of 4.1 is a person reading each round (Zhang 2012, CARF 2014, Fozo 2010, the mycofactocin work "with manual review at each iterative step"). Hit lists grow roughly in proportion to genome count, so a round that held 10^2 new sequences in 2010 holds 10^4 to 10^5 now. Three changes follow:
   - Collapse each round to cluster representatives before review. Fozo already re-seeded from 95%-identity representatives.
   - Automate the checks the papers did by eye: required catalytic motif present (as HEPN and CARF vetting did), domain coverage, no new foreign Pfam domain entering the alignment (a corruption alarm), and an HHpred-style confirmation for borderline hits.
   - Keep a per-round log so a person can audit a sample instead of every hit.
2. **Transitive and context-seeded search explodes.** Two rules saturated at 10^2–10^3 genomes: "each retrieved sequence used as the query ... until no new sequences" (CRISPR 2006), and "neighbors become seeds" (Iyer 2006, CARF 2014). At 10^6 genomes the frontier grows combinatorially, and promiscuous domains (HTH, P-loop, HNH, TPR) bridge unrelated families. Needed:
   - de-duplication of seeds by cluster;
   - a bounded frontier with a provenance graph of seed, hit and round;
   - a promiscuous-domain blocklist (the papers did this ad hoc: the HTH family dropped from the DI positive set, the RuvC-only query, KOG masking);
   - a reverse-search check before any new seed is promoted.
3. **E-values no longer mean the same thing.** Inclusion at 0.01 against NR-2003 (about 10^6 sequences) is roughly 1,000 times more permissive in raw score than the same E against 10^9 sequences. Some authors already fixed the effective database size (Price: 10^8; Spang: 10^8). EFI moved to a database-independent alignment score, and TIGRFAM uses bit-score cutoffs. The infrastructure should store bit scores, run HMMER, MMseqs2 and Infernal with a fixed `-Z`/`--domZ` or Infernal `-Z`, and treat recipe E-thresholds as bit thresholds at a reference size.
4. **Clonal inflation breaks the statistics.** DI Z-scores, binomial windows, hypergeometric profiles and PPP binomials all assume independent genomes. The papers controlled this by hand with representative subsets: 383 of 1,055, 659 of 2,262, 450 of 503, one strain per species (Boyer), close strains excluded (Rodionov). At 10^6 genomes (MAG-heavy, about 10^4–10^5 species) this must be systematic: species-representative weighting (GTDB clusters), STRING-style collapse of concordant clades, and per-species de-duplication of neighborhoods before counting.
5. **All-vs-all steps are quadratic.** These include COG/arCOG construction, Tamames' pairwise gene-order comparison, Rogozin's conserved-pair enumeration, PPP over "pre-computed BLAST results on 1450 ... genomes", NISE within-EC all-vs-all, and BLISS IGR all-vs-all. At scale they must be replaced by precomputed clusters (MMseqs2 linclust/cascaded), orthology from cluster membership plus synteny, and profile vectors over clusters. PPP has to be re-expressed over cluster hierarchies: walk the cluster tree instead of a BLAST hit list.
6. **Curated resources do not grow with the corpus.** SEED subsystems, TIGRFAM trusted cutoffs and RegPrecise CRON curation are expert-per-family. They stay useful as fixed rules (Q-FAM/Q-SYS), but the recipes that created them ("manual review of ... cutoff scores") cannot be rerun at 10^6 genomes. Treat their published outputs as gold sets, not as procedures to repeat.
7. **Gene-call errors grow with MAG share.** Fozo's tblastn and six-frame IGR scans, Poole's tblastn for unannotated CT modules, Shmakov's re-calling within ±20 kb, and Zallot's "15% to 20%" wrong starts show that classic discovery already depended on searching outside the annotation. At scale this needs a contig-level protein search service (six-frame or ORF-level index) and uniform re-calling. Truncated loci at contig edges must be flagged, as Axen excluded loci "truncated by the beginning or end of the scaffold".
8. **RNA discovery depends on targeting.** All-IGR BLAST clustering was feasible at about 3.2 Gnt (2007). At 10^6 genomes the IGR space is about 10^11–10^12 nt, so the 2015 strategy is the scalable one: gene- or context-anchored IGR subsets (6 kb around anchors, about 7 Mbp in total). CM searches need Infernal's filter pipeline and a fixed `-Z`.
9. **Fixed windows meet fragmented assemblies.** ±5 genes and ±10 kb windows still compute cheaply, but MAG contigs cut windows short. Neighborhood statistics should count only windows fully contained in a contig, or model truncation.

## 7. Ten best replay candidates (published lists, alignments or models in `lit/supp/`)

| # | Paper | Search to replay | Artifact (path under `bioinformatics-synteny/lit/supp/`) | Why |
|---|---|---|---|---|
| 1 | 2011_Makarova_Koonin_Defense_Islands | DI enrichment: 383 representative genomes, ±5 genes, island = run bounded by negative genes, Bonferroni Z | `2011_Makarova_Koonin_Defense_Islands/Makarova_DI_add.tar.gz` (holds `all_PD2_islands.txt`, `all_Mob_islands.txt`, `mk.xls.csv`, `mk.yls.csv`, `Readme.txt`) | Fully specified statistics and published island lists; the canonical Q-NBR acceptance test |
| 2 | 2009_Makarova_Wolf_Koonin_Type2_Toxin_Antitoxin_Mobile_Stress_Islands_BiolDirect | PSI-BLAST (h 0.01 or 0.1 by length, stop before first false positive), then co-directed neighbors with IGR <100 bp, then BLASTCLUST and a census over 750 genomes | `.../1745-6150-4-19-S1.xls` … `S14.xls` (S14 is 4.7 MB), alignments `S2/S3/S5/S8/S9.zip`, `S7.txt` | Large published census with stated parameters for every step |
| 3 | 2010_Weinberg_Wang_Ruzzo_Breaker_104_Candidate_Structured_RNAs | IGR BLAST clustering, CMfinder, CM expansion and filters | `.../gb-2010-11-3-r31-S4.GZ` (raw `sto/` alignments), `S5.GZ` (`clean-sto/`), `S10.GZ` (pipeline source code), `S11.TAB` | Stockholm alignments for 104 motifs plus code; recall test for Q-RNA |
| 4 | 2006_Iyer_Burroughs_Aravind_Prokaryotic_Ubiquitin_Antecedents | transitive PSI-BLAST/HMMER (h 0.01, to convergence) from beta-grasp, E2 and JAB seeds; neighborhood gi lists | `.../prok_ub_supplement_file_3.html` (3 MB: "hit lists for PSI-BLAST and HMM searches", with starting points, gis and e-values), `_file_1.html` (neighborhoods, architectures, alignments), `_file_2.html` (gis with genomic positions) | The only classic Aravind paper with seed-by-seed hit lists; direct Q-ITER recall test |
| 5 | 2015_Kapitonov_Makarova_Koonin_ISC_Transposons_Cas9_Homologs | single seed WP_044523571.1, PSI-BLAST default, 6 iterations, then 1,227 Cas9 homologs | `.../JB.00783-15_zjb999093943so1.pdf` (Table S1) | One seed, one iteration count, one published hit count; tests E-value drift between NR 2014 and now |
| 6 | 2014_Makarova_Koonin_Aravind_CARF_WYL_Domains | PSI-BLAST (0.01, composition and filters off, motif vetting per round), ±5–10-gene neighborhoods, co-occurrence over 659 genomes | `.../DataSheet1.XLSX` … `DataSheet10.ZIP` | Tables plus alignments; tests automated motif vetting against the published curated set |
| 7 | 2011_Atkinson_Tenson_Hauryliuk_RSH_Superfamily_ppGpp_Synthetases_Hydrolases | Pfam seed, subgroup HMMs, 2 HMM iterations over 1,072 proteomes, per-domain E filters, 30 subgroups | `.../pone.0023479.s005.xlsx`, `s006.xls`, `s001-s004.pdf` | 2,077-sequence classified set; tests iterative HMM expansion and subgroup assignment |
| 8 | 2014_Axen_Erbilgin_Kerfeld_Taxonomy_Bacterial_Microcompartment_Loci | modified PF00936 HMM by hmmsearch (E ≤ 1e-5), ±10 kb loci, Pfam 27 at 0.01, LoClass | `.../pcbi.1003898.s004.hmm` (model), `s003.fasta`, `s006.fasta`, `s001.xlsx` and `s002.xlsx` (loci), `s005.cys` (network) | Model plus locus tables; clean Q-FAM-to-Q-NBR replay |
| 9 | 2009_Boyer_Fichant_Berthod_Vandenbrouck_Attree_T6SS_Conserved_Gene_Clusters_BMCGenomics, with 2016_Abby_Cury_Guglielmini_Rocha_Protein_Secretion_Systems_TXSScan | T6SS COG rule (≥80% coverage, E ≤ 1e-6, <5 intervening genes, ≥5 COGs) vs MacSyFinder models | Boyer `.../1471-2164-10-104-S1.pdf` … `S8.doc`, `S7.tgz`; Abby `.../srep23080-s1.pdf`, `s2.xls`, `s3.xls` | Two generations of rule-based detection on overlapping genomes; Q-SYS concordance test |
| 10 | 2002_Rogozin_Makarova_Wolf_Koonin_Connected_Gene_Neighborhoods_NAR | conserved COG pairs (≤2 intervening genes, ≥3 genomes), EM_TRAILS arrays, single-linkage neighborhoods | `.../edges1.dat`, `mtriplet.dat`, `purgedmaxwalks.dat`, `neighbor.txt`, `end_max_walks.pl`, `purge1.pl`, `Em_walks.doc`, `readme` | Code plus intermediate data; tests the neighborhood-graph stage at historical scale |

Alternates:
- `2005_Corbino_Barrick_Lim_Breaker_Riboswitch_Candidates_Alpha_Proteobacteria/gb-2005-6-8-r70-S2.sto` … `S6.sto` (5 small Stockholm alignments, a quick Infernal unit test).
- `2011_Haft_Varghese_GlyGly_CTERM_Rhombosortase_Synteny/pone.0028886.s004.doc` (complete GlyGly-CTERM list from 108 genomes; PPP and HMM replay with public TIGR03501).
- `2012_Wolf_Makarova_Yutin_Koonin_Updated_arCOGs_Complex_Archaeal_Ancestor_BiolDirect/ar14.arCOG.csv` + `ar14.fa.gz` (full arCOG membership; Q-CLUST concordance).
- `2006_Makarova_Koonin_Prokaryotic_RNAi_CRISPR_Hypothesis/1745-6150-1-7-S1.xls`, `S2.txt`, `S3.txt` (Cas tables and motif alignments).
- `1999_Wolf_Aravind_Grishin_Koonin_Aminoacyl_tRNA_Synthetases_Evolution_GenomeRes/aars_*.pure.ali` (per-aaRS alignments).

Replay caution: the classic hit lists use GI numbers, which NCBI retired. Each replay first needs a GI-to-accession/UniProt mapping. Recall should be measured at the stage whose genome count matches the paper's era (for example, the 383 DI genomes, the 750 TA genomes, the 2,751 CRISPR genomes), and only then at full scale.

## 8. The 30 most influential papers in the slice (covered in depth)

Iterative search and domain lego: 1999_Aravind_Walker; 1999_Wolf_Aravind aaRS; 2000_Aravind UDG; 2003_Iyer DPBB; 2006_Iyer_Burroughs Ub; 2009_Burroughs E1; 2012_Zhang polymorphic toxins; 2013_Anantharaman HEPN; 2013_Iyer phage hypermodification; 2015_Burroughs nucleotide networks; 2016_Zhang Crn.

Koonin/Makarova: 2002_Rogozin neighborhoods; 2003_Tatusov COGs; 2006_Makarova CRISPR; 2007_Makarova arCOGs; 2009_Makarova Ago; 2009_Makarova TA; 2010_Fozo type I TA; 2011_Makarova Defense Islands; 2014_Makarova CARF/WYL; 2015_Makarova CRISPR classification; 2015_Shmakov Class 2.

Context methods: 2000_Snel STRING; 2003_vonMering STRING; 2001_Tamames; 2004_Bowers Prolinks.

RNA: 2007_Weinberg 22 RNAs; 2010_Weinberg 104 RNAs.

Missing genes: 2005_Overbeek SEED; 2006_Haft PPP; 2016_Leyn RegPredict regulons.

Every other usable paper has at least one row. The 33 defect stems are listed in section 1.
