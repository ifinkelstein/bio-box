# Prior art: what already exists, what `bio-box` reuses, and what is actually new (v3.0)

## Summary

Several systems already cover large slices of this project, so the honest claim is narrow:

- **Genomic context search at scale exists.** Gaia (gLM2 embeddings over 85M protein clusters), STRING, GeCoViz, AnnoTree, EFI-GNT, FlaGs and GCsnap all show neighbourhoods. Gaia and gLM2 learn context from an open metagenomic corpus of 3.1 Tbp.
- **Dark protein families from metagenomes have been catalogued with genomic context.** Pavlopoulos et al. 2023 report 106,198 novel families; FESNov reports about 404,000 families with context-based function prediction.
- **Structure-level dark matter has been mapped.** AFDB clusters and Durairaj et al. 2023 did this, and Foldseek/ProstT5 make structure search routine.
- **Planetary sequence search exists.** Logan-Search and MetaGraph cover it at the nucleotide/k-mer level.
- **Mobile-element and defense classification are solved well enough to reuse.** geNomad, CheckV, IntegronFinder, MacSyFinder with DefenseFinder/CONJScan, and PADLOC cover them.
- **Co-occurrence statistics corrected for phylogeny exist.** Coinfinder does this for pangenomes.
- **Corpus-scale context-scored clustering exists.** FLSHclust (Altae-Tran et al. 2023) clustered 8 billion proteins and ranked clusters by CRISPR proximity with an effective-sample-size correction.
- **Sub-ORF neighbourhood analysis is published.** Domainator (Johnson et al., NAR 2025; preprint 2024) clusters the residue not covered by any annotated domain, extracts framework-bracketed hotspots and subtracts expected partners. Rousset et al. 2022 bracket variable cargo between two framework genes and count distinct arrangements.

What `bio-box` adds is the combination, not any single layer:

1. one coordinate tape where protein domains (sub-ORF, including uncovered segments, as Domainator already does per query), structured RNAs and recombination sites sit side by side with junction topology;
2. co-occurrence statistics that account for mobile vehicles, contig edges and taxonomy uncertainty, with stated nulls and calibrated FDR;
3. evidence-gated closure against PDB, Swiss-Prot, PaperBLAST and Fitness Browser at the domain-envelope level;
4. stable IDs and a release model that lets the corpus accrete in stages on modest hardware.

The plan reuses existing tools for every layer where one exists (`PLAN.md` §3) and does not re-implement them. Whether the combination finds anything the existing systems miss is an empirical question, answered by the benchmarks in `PLAN.md` §5 against the external baselines listed below.

---

## Comparison matrix

| Platform / project | Group, year | Scale | Unbinned / MGE content | Unit of context | Structure | ncRNA in context | Experimental closure | Engine | What `bio-box` reuses or must beat |
|---|---|---|---|---|---|---|---|---|---|
| **STRING v12** | Szklarczyk, von Mering, Bork; 2023 | 12,535 organisms | No | Whole-protein orthologs (eggNOG) | AFDB links | No | Experiments, text mining | Postgres + flat files | Baseline for B1: STRING combined score on the same pairs. |
| **GeCoViz** | Botas, Huerta-Cepas et al.; 2022 | ~12k reference genomes | No | eggNOG/Pfam/KEGG neighbourhoods | No | No | No | Web API | Visual reference; small genome set. |
| **AnnoTree / AnnoView** | Mendler, Doxey et al.; 2019–2024 | ~85k GTDB genomes | No | Pfam/TIGRFAM/KEGG per genome | No | No | No | MySQL | Distribution of family presence across GTDB; reusable as a sanity check for pangenome planes. |
| **EFI-GNT** | Gerlt, Zallot et al.; 2015–2024 | UniProt/ENA | No | ±10 neighbours of an SSN | AlphaFold links | No | Swiss-Prot flag | Per-job SQLite | Per-query neighbourhood networks; minutes to hours per job. |
| **FlaGs / webFlaGs** | Saha, Atkinson et al.; Bioinformatics 2021 (online Dec 2020) | NCBI proteins on demand | No | Flanking-gene clusters + tree | No | No | No | Per-job | Per-query neighbourhood conservation with a tree; reference for visual output. |
| **GCsnap / GCsnap2** | Pereira et al., J Mol Biol 2021; GCsnap2 preprint 2025 | Any NCBI/UniProt genomes on demand | No | Protein-coding neighbourhoods | Links | No | No | Python; HPC version in 2025 | Per-query context comparison; a reference for output formats. |
| **Gaia (on gLM2)** | Tatta Bio; bioRxiv 2024, Science Advances 2025 | 85M protein clusters from 131,744 genomes | Partly (MGnify/IMG-derived) | gLM2 contextual embeddings, HNSW search | Structure predictions shown | No | No | Vector DB (HNSW) | Closest existing system. Must be compared directly on context-aware retrieval. |
| **gLM2 + OMG corpus** | Cornman, Hwang et al.; bioRxiv 2024, ICLR 2025 | 3.1 Tbp, 3.3B CDS from IMG + MGnify | Yes | Mixed-modality: protein + intergenic DNA tokens | Implicit | Intergenic DNA tokens | No | Transformer; open weights | Baseline for masked-slot prediction (if `gLM-Domain` cannot beat gLM2, it is not built). OMG is a candidate source corpus; check its license. |
| **gLM** | Hwang, Girguis et al.; Nat Commun 2024 | ~7M MGnify contigs, 15–30 genes | Yes | Whole-protein ESM2 embeddings | Implicit | No | No | Masked transformer | Superseded by gLM2 as a baseline. |
| **Bacformer** | Wiatrak et al.; bioRxiv 2025 (preprint) | >1.3M bacterial genomes | No | Ordered proteome as a sequence | Implicit | No | No | Transformer | Operon and protein-interaction prediction baseline. Preprint-level claims. |
| **Evo / Evo 2 (OpenGenome, OpenGenome2)** | Nguyen, Hie et al., Science 2024; Brixi et al., bioRxiv 2025 / Nature 2026 | Evo 2: 8.8 Tbp across all domains, 1 Mbp context | Includes phage and plasmid sequence | Single nucleotides | No | Implicit | No | StripedHyena / StripedHyena 2 | Generative nucleotide models. Not a relational query engine; usable as a scoring baseline for intergenic features. |
| **Pavlopoulos et al. (Novel Metagenome Protein Families)** | Kyrpides lab and consortium; Nature 2023 | 26,931 metagenomes, 1.17B novel proteins → 106,198 families | Yes | Families annotated with taxonomy, habitat, neighbourhood | Predicted models for some families | No | No | Static release | Source of dark-family definitions and a cross-check for `subfamily_50_id` clusters. |
| **FESNov** | Rodríguez del Río, Huerta-Cepas et al.; Nature 2024 | ~150k environmental genomes → ~404k novel families | Yes (uncultivated taxa) | Families with genomic context and structure | Structural alignments | No | 2 families tested experimentally | Static catalogue | Must be reproduced: context-based function prediction for 32% of families. A B1-style recall test against FESNov predictions. |
| **AFDB clusters** | Barrio-Hernandez, Steinegger, Beltrao et al.; Nature 2023 | 214M AFDB models → 2.3M structural clusters | Via UniProt | Monomer structure clusters | Yes | No | No | Foldseek clusters | Fold assignment source for AFDB-covered proteins. |
| **Durairaj et al.** | Schwede lab; Nature 2023 | AFDB-wide sequence similarity network | Via UniProt | Monomer families and folds | Yes | No | No | Static release | Dark-fold reference set. |
| **Foldseek / ProstT5 / AFDB / ESMAtlas / BFVD** | Steinegger, Rost et al.; 2022–2025 | 214M AFDB; ESMAtlas30 distributed as 3Di; 351k BFVD | Via ESMAtlas / BFVD | Monomer 3Di | Yes | No | PDB100 matches | Memory-mapped C++ | Reused directly for fold assignment (D-10). Has no genomic context of its own. |
| **SocialGene** | Chase, Medema et al.; 2023 | Custom RefSeq/MIBiG subsets | User supplied | Pfam domains in a property graph | No | No | MIBiG links | Neo4j | Graph model that does not scale to 10⁸+ proteins on one node. |
| **GlobDB** | Ettema, Spang et al.; 2024–2025 | r226: ~306k species representatives | MAGs only | Per-genome annotations | No | No | No | Per-genome files | Candidate Tier 1 source. Overlaps GTDB, so identity resolution is required. |
| **SPIRE** | Bork lab; 2024 | ~100k metagenomes, 1.16M MAGs, 35B proteins | Assemblies held; MAGs indexed | eggNOG/Pfam | No | No | No | REST + Parquet | Candidate Tier 3 source; no cross-contig synteny statistics. |
| **Logan + Logan-Search** | Chikhi et al.; bioRxiv 2024, 2025 version | 27M SRA accessions (Dec 2023 freeze), 44.1 Pbp raw → 0.9 Pbp contigs | Yes | Unitigs/contigs; k-mer index over unitigs (Logan-Search, minutes per query) | No | No | No | k-mer index on AWS | Last-stage source, compute-to-data only (D-05). Logan-Search is a k-mer index; LexicMap is an unrelated aligner and was misattributed in v1. |
| **MetaGraph** | Karasikov, Kahles, Rätsch et al.; Nature 2025 | 18.8M sequence sets, 210B amino acids | Yes | Annotated de Bruijn graph | No | No | No | Compressed k-mer index | Reference for "which samples contain this sequence" queries; reuse instead of building a nucleotide index. |
| **Coinfinder** | Whelan, McInerney et al.; Microbial Genomics 2020 | Single-species or clade pangenomes | Via input | Gene presence/absence, phylogeny-aware | No | No | No | C++ / R | Baseline for co-occurrence statistics corrected for phylogeny (B1 ablation). |
| **geNomad** | Camargo, Kyrpides et al.; Nat Biotechnol 2024 (online 2023) | Classifier | Yes | Plasmid / virus / chromosome per contig, ICTV taxonomy | No | No | No | Python + MMseqs2 | Reused directly for the mobility regime and viral taxonomy (D-11). |
| **DefenseFinder / PADLOC / MacSyFinder** | Tesson, Bernheim et al.; Payne, Jackson et al.; 2021–2025 | Classifiers | Yes | Rule-based system calls | No | PADLOC-DB includes CMs (e.g. retron ncRNA) | Literature-backed models | DefenseFinder: Python on MacSyFinder. PADLOC: its own engine (HMMER searches plus its own presence/absence and synteny rules per system model), not MacSyFinder | System calls stored as sidecars, never as tape tokens. Excluded from the temporal-holdout dictionary (D-24). PADLOC-DB v2.0.0 adds hundreds of candidate systems (HEC, PDC), some unpublished, so each release's candidate content is checked against the B1 holdout (test S-11). |
| **FLSHclust** | Altae-Tran, Zhang et al.; Science 2023 | 8.8 Tbp, 8×10⁹ proteins, 10.2M CRISPR arrays | Yes (metagenomic contigs ≥2 kbp) | LSH clusters at 90/50/30% identity (1.3×10⁹ / 646.4M / 499.9M), scored by array proximity (3 kb operon, 10 kb repeat window) divided by N_eff | AlphaFold2 per candidate | CRISPR arrays only | Wet-lab validation of selected systems | Custom LSH clustering on a 2-node grid | The canonical "cluster, then score by context" recipe: its scoring step is the tape + pair-statistics design (§4.3). Replay: re-derive the naive and enhanced scores from the Zenodo loci (~130,000 clusters at enhanced ≥0.35, N_eff ≥3; 188 systems). No HMMs or MSAs released, so the replay is partial. |
| **Domainator** | Johnson et al., New England Biolabs; NAR 2025 (gkae1175), repo v0.12.0, MIT | Per-workflow on downloaded GenBank: 13,116 complete genomes (CARF), 2,718 *E. coli* genomes (~30 GB, P4 hotspots), ~366 GB of good-quality genomes (RM) | Via input | Domain-annotated neighbourhoods; uncovered-segment extraction; framework-bracketed hotspots; partner subtraction; similarity networks with a threshold table | ESM-2 3B 3Di + Foldseek (or reseek) per neighbourhood | Infernal CMs over contigs, no typed RNA features on a tape | 1 of 2 predicted REs active (Gba686I); Lamassu variant active in vivo | Composable CLI, Python library, browser GUI; GenBank as the carrier format; reference DBs held in memory | Publishes uncovered-region clustering (1,844 CARF-only proteins → 1,524 fragments ≥50 aa → 94 clusters), hotspot extraction (910 regions, Jaccard > 0.7) and orphan-partner subtraction (205,039 → 124,853 → 2,628 contigs). No null model or significance test; counts are per contig; E-values use `Z = 1000` by default. Replays (data on Zenodo 10.5281/zenodo.10989173, scripts on Zenodo 14056380): case 1 at ST2 (S2C), case 3 at ST3–ST4 (S5–S6), case 2 at ST4 (S6; S2C's 20-per-species cap gives only a partial check). bio-box borrows its profile-coordinate fields, agent JSON conventions, network threshold tables and run guardrails, and uses it as the candidate workbench (SEARCH_PLAN SD-17..SD-21) and as an independent implementation in tests (S-02d domain calls and neighbourhoods; `build_ssn` as the O-SSN oracle, TEST_ARCHITECTURE §2). |
| **Rousset et al. (P2/P4 hotspots)** | Rousset, Bikard et al.; Cell Host Microbe 2022 | 20,125 *E. coli* GenBank genomes; 21,738 RefSeq genomes for the census | Yes (prophages, satellites) | Cargo between two framework genes (P4: `psu`–`int`; P2: `gpA`–`gpQ`), clustered (MMseqs2 `-c 0.6`) into distinct arrangements | No | No | 7 of 18 P4 systems defend; 4 of 9 P2 candidates (8 systems plus *old-tin*) defend | BLASTP anchor + window; MacSyFinder census | Framework-gene bracketing and arrangement counting (5,251 P4 loci → 541 clusters / 318 arrangements). Arrangement counts are one of the dependence-aware counts every search summary reports (SD-12). Bit-score thresholds (32/27/27), not E-values. |

---

## What is different, stated narrowly

1. **Sub-ORF units on a shared tape with ncRNA.** Most context tools treat a protein as one token and ignore structured RNAs. Sub-ORF analysis itself is not new: Domainator (2024) already extracts the residue not covered by any annotated domain and clusters it, and brackets hotspots between framework genes. What `bio-box` adds is the combination: domain segments, uncovered segments, riboswitches (with pocket subtypes), CRISPR arrays and `attC` sites on the same ordered tape with junction topology, corpus-wide rather than per query, so domain fusion (same ORF) and operon adjacency (neighbouring ORFs) are queried the same way and feed the dependence-aware statistics in item 2. gLM2 includes intergenic DNA but not typed RNA features; whether typed features add recall is tested, not assumed.
2. **Statistics that account for vehicles, edges and taxonomy uncertainty.** STRING and Coinfinder correct for phylogeny; FLSHclust divides its context score by N_eff and Rousset counts arrangements rather than loci. `bio-box` adds collapsing of shared mobile vehicles, contig-edge censoring, explicit per-contig taxonomy confidence and a replicon-matched test for MGE-confined pairs, with calibrated FDR (`PLAN.md` §4.3). The claim is checked against simulated planted truth, null worlds with phylogeny and HGT but no coupling (Zombi, ALF as fallback; `planning/tests/TEST_ARCHITECTURE.md` §6.2, MR-18), a temporal holdout and published replays (S-01).
3. **Envelope-level closure.** Closure is attributed per domain envelope and per active-site constellation, so a dark domain fused to a solved one is not marked solved.
4. **Staged, reproducible accretion on modest hardware.** The corpus grows stage by stage on a CPU-only lab server, with stable content-derived IDs, ledgers and release crosswalks ([`STAGES.md`](STAGES.md)). Rented bare metal is used only once a stage outgrows the lab server.

## Modern pattern: reduce before search

No 2023–2026 corpus-scale discovery paper read for the search plan runs an iterative profile search over its whole corpus (`planning/search/05_external_modern.md` §a). Each reduces first, then runs the expensive tool on 10³–10⁶ objects:

| Pattern | Reduction | Example | What `bio-box` provides |
|---|---|---|---|
| A. Cluster, then score by context | 8×10⁹ proteins → 499.9M clusters at 30% → ~130,000 context-enriched clusters → 188 systems | FLSHclust (Altae-Tran 2023) | NR90/NR50/subfamily clusters (L2); tapes and `pair_stats` as the scoring step |
| B. Bracket a locus with two framework genes | BLASTP anchor at E ≤ 1e-10, next 10 genes → 5,251 P4 loci → 318 arrangements | Rousset 2022; Domainator case 2 | `bio-box nbr` windows; arrangement counts (SD-12) |
| C. Subtract the expected partner | 205,039 MTase neighbourhoods → 124,853 orphans → 2,628 contigs by 3Di search | Domainator case 3 | Neighbourhood predicates over stored library hits; 3Di of NR50 from tier ST3, encoder chosen by the SD-21 bake-off |
| H. Subtract annotated domains, cluster the residue | 1,844 proteins → 1,524 uncovered fragments ≥50 aa → 94 clusters | Domainator case 1 | `segments.faa` (domain envelopes plus uncovered segments ≥50 aa) in the search pack |

`bio-box` adopts the same order (SD-03 in [`planning/search/SEARCH_PLAN.md`](planning/search/SEARCH_PLAN.md)): iterative rounds on NR50 representatives, a final round on NR90 projected to members, and exhaustive search of all proteins only as a recall check. Candidate-linear work (HHpred, DALI, AF3/ColabFold) stays external.

## What `bio-box` does not attempt

- Nucleotide-level search of the SRA: use Logan-Search or MetaGraph.
- Learned context embeddings as the primary product: Gaia/gLM2 already provide them. A `gLM-Domain` model is built only if it beats gLM2 on held-out masked-slot prediction (D-16).
- Re-deriving mobile-element, defense or secretion classifiers: geNomad, MacSyFinder models and PADLOC are reused.

## Sources

- [OMG / gLM2 (ICLR 2025)](https://proceedings.iclr.cc/paper_files/paper/2025/hash/dbb40375daa3c39f1c098b12608aed12-Abstract-Conference.html)
- [Gaia (Science Advances 2025, PMC)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12180486/)
- [FESNov (Nature 2024, PMC)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10849945/)
- [Pavlopoulos et al. 2023](https://pmiweb.ornl.gov/publication/pavlopoulos-et-al-2023-nature/)
- [Bacformer preprint](https://www.biorxiv.org/content/10.1101/2025.07.20.665723.full.pdf)
- [Logan preprint](https://www.biorxiv.org/content/10.1101/2024.07.30.605881.full.pdf)
- [Logan GitHub](https://github.com/IndexThePlanet/Logan)
- [MetaGraph (Nature 2025, PMC)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12657231/)
- [Evo 2 manuscript](https://arcinstitute.org/manuscripts/Evo2.pdf)
- [GCsnap](https://publikationen.uni-tuebingen.de/xmlui/handle/10900/121509)
- [GCsnap2](https://arxiv.org/abs/2505.02195)
- [FlaGs (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8189683)
- [FLSHclust: Altae-Tran et al. 2023, Science (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10910872); DOI 10.1126/science.adi1910 (verify); data and scripts Zenodo `10.5281/zenodo.8371343`
- [Domainator, Johnson et al. 2025, NAR 53:gkae1175 (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11754643), doi:10.1093/nar/gkae1175, PMID 39657740; code [nebiolabs/domainator](https://github.com/nebiolabs/domainator) (v0.12.0, MIT); workflows [nebiolabs/domainator_examples](https://github.com/nebiolabs/domainator_examples), data Zenodo `10.5281/zenodo.10989173`, scripts Zenodo `14056380`
- [Rousset et al. 2022, Cell Host Microbe (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9122126/)
- [PADLOC: Payne et al. 2021, Nucleic Acids Res](https://doi.org/10.1093/nar/gkab883); web server [Payne et al. 2022 (PMC)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9252829/); [PADLOC-DB](https://github.com/padlocbio/padloc-db)

Durairaj 2023, Barrio-Hernandez 2023, Coinfinder 2020 and geNomad 2024 are cited from the round-2 audits and were not re-checked here.
