# Modern large-scale sequence search outside the local corpus (2020–2026)

Status: external-literature round for the search-infrastructure brief (`00_brief.md`), 2026-10-29.
Scope: how discovery groups and tool authors actually run searches at corpus scale, what the published
throughput and memory numbers are, and what those numbers imply for ~732k GTDB genomes
([GTDB R10-RS226: 732,475 genomes, 143,614 species clusters](https://gtdb.ecogenomic.org/stats/r226)),
~3e9 proteins / ~1e9 unique, on 56 cores with a 200 GB cap.

Companion file: `recipes_external.tsv` (one row per search, schema per `00_brief.md`).
Papers already covered by `../../../bioinformatics-synteny/lit/INDEX.md` (Gao 2020, Doron 2018,
Millman 2020/2022, Vassallo 2022, Mestre 2020, Pinilla-Redondo 2020, Shmakov/Makarova/Koonin CRISPR
series, Tesson 2022, Breaker/Weinberg ncRNA series, EFI 2015/2019) are **not** re-extracted here.
Overlaps with `PRIOR-ART.md` (Gaia/gLM2, FESNov, Pavlopoulos, MetaGraph, SPIRE, geNomad, Coinfinder)
are not repeated either; this file adds the *search* layer those systems sit on.

Every number below carries a URL. Numbers I computed are marked **[derived]** with the arithmetic shown.
Numbers the authors measured are marked **[measured]**.

---

## (a) Modern recipe patterns, as step sequences

Ten patterns cover essentially every external search in `recipes_external.tsv`. The scale-relevant
observation is that **none of the 2023–2026 discovery papers runs an iterative profile search against the
whole corpus.** They all reduce first — by clustering, by bracketing a locus, by subtracting an expected
partner, or by predicting structure — and only then run the expensive tool on 10^3–10^6 objects.

### Pattern A — Cluster the whole corpus, then score clusters by context
*Altae-Tran 2023 (FLSHclust), the canonical corpus-scale discovery recipe.*

1. Assemble the corpus: 8.8 Tbp of prokaryotic genomic + metagenomic contigs (NCBI, WGS, JGI), drop
   metagenomic contigs < 2 kbp → **8 billion proteins, 10.2 million CRISPR arrays**.
2. Call genes (GeneMark); call arrays with four finders over six runs (PILERCR, CRT, CRISPRFinder, CRONUS),
   deduplicate by interval keeping the best array per interval.
3. Cluster *all* proteins at 90%, then 50%, then 30% identity with an LSH algorithm, not an alignment-first
   one: L = 24 hash functions, k = 12, 3 dropped positions, 2 tolerated mismatches, MCMC-optimised hash
   family, two representatives per bucket, graph simplification, community detection, greedy clustering
   → **1.3e9 / 646.4M / 499.9M clusters**.
4. Define operons per contig; call a protein operon "array-associated" if it lies within **3000 bp** of an array.
5. Score each cluster: naive score = weighted fraction of non-redundant members in an array-associated
   operon ÷ N_eff; enhanced score repeats this using divergent direct repeats found within a **10 kbp**
   window of each member. Known *cas* genes score 0.44 naive / 0.72 enhanced.
6. Cut at enhanced ≥ 0.35 with N_eff ≥ 3 → **~130,000 of 16 million eligible clusters**.
7. Curate manually → **188 previously unreported CRISPR-linked modules**.
8. Only now run expensive tools: HHpred pairwise scores → UPGMA dendrogram of Cas7 relatives; AlphaFold2;
   wet-lab validation.

Source: <https://pmc.ncbi.nlm.nih.gov/articles/PMC10910872> (Methods summary).

**What this means for the plan:** the enrichment statistic, not the search, is the discovery engine. The
plan's tapes and pair statistics *are* steps 4–6. Step 3 is the prerequisite the plan already has
(`subfamily_50_id`-style clusters). The acceptance test is whether a replay reproduces the ~130,000
surviving clusters from the published cluster table.

### Pattern B — Bracket a hypervariable locus with two conserved framework genes
*Rousset 2022; Domainator 2024 case 2.*

1. Pick two genes that reliably flank the variable cargo. P4-like satellites: `psu` and `int`.
   P2-like phages: `gpA` and `gpQ`.
2. BLASTP the anchor at E ≤ 1e-10; scan the next 10 downstream genes by annotation for the second
   framework gene; keep everything between them.
3. → 5,251 P4 loci in 26% of 20,125 *E. coli* GenBank genomes; 18,150 P2 loci.
4. Cluster the cargo proteins (MMseqs2 `-c 0.6`): 10,860 P4 proteins → 541 clusters / 318 arrangements;
   66,498 P2 proteins → 846 clusters / 1,650 arrangements.
5. Curate arrangements by occurrence count (≥ 5 for P4 covering 94.4% of loci; ≥ 10 for P2 covering 82%);
   drop P2 arrangements containing only proteins < 120 aa.
6. Annotate domains with HHpred / HH-suite3, reporting Pfam matches above **80% probability**; Phobius for
   TM and signal peptides.
7. Clone and challenge: 18 P4 systems tested → 7 defended; 8 P2 systems + *old-tin* → 4 defended.
8. Separately, turn the best hit into a rule-based census: MacSyFinder 1.0.2 with PF13175 / PF13304 /
   PF14491 at **bit-score thresholds 32 / 27 / 27** (not E-values) across 21,738 RefSeq genomes.

Sources: <https://pmc.ncbi.nlm.nih.gov/articles/PMC9122126/> (STAR Methods);
<https://pmc.ncbi.nlm.nih.gov/articles/PMC11754643> (Domainator case 2, 910 regions, Jaccard > 0.7
arrangement clustering, 2,718 genomes ≈ 30 GB).

**What this means for the plan:** this is cheap — a BLASTP anchor plus a window — and it is the single
most productive discovery recipe in the anti-phage literature. It needs locus coordinates and
framework-gene calls, which the tape already provides. The arrangement/Jaccard clustering step is the
piece most likely missing.

### Pattern C — Subtract the expected partner to find the orphan, then search by structure
*Domainator 2024 case 3 — the cleanest published pLM-3Di + Foldseek neighbourhood recipe.*

1. Seed with REBASE Gold Standard DNA methyltransferases; `phmmer` prokaryotic genomes
   (~366 GB of high-quality GenBank genomes) for MTases.
2. Extract the MTase **± 4 ORFs**.
3. **Subtract**: discard every neighbourhood that contains a `phmmer`-detectable restriction endonuclease
   → 205,039 MTase neighbourhoods become **124,853 orphan MTases**.
4. Recode the remaining neighbour proteins to Foldseek 3Di with **ESM-2 3B** (no structure prediction).
5. Foldseek-search the 3Di for Type II restriction endonucleases → **2,628 contigs** with putative
   cryptic RM systems.
6. SSN the predicted REs, inspect representatives, test two → one active (Gba686I).

Source: <https://pmc.ncbi.nlm.nih.gov/articles/PMC11754643> (case 3).

**What this means for the plan:** "expected partner is absent" is a first-class query the plan should
support directly, and it is the strongest argument for storing a pLM-derived 3Di string alongside every
protein (see the Foldseek/ProstT5 costs in §c.6).

### Pattern D — Structure-first dark-matter triage of the whole structure universe
*Barrio-Hernandez 2023; Lau 2024 (TED).*

1. MMseqs2 pre-cluster 214,684,311 AFDB sequences at **50% identity / 90% overlap**, highest-pLDDT member
   as representative → **52,327,413 clusters**.
2. Foldseek cluster the representatives: **E ≤ 0.01**, structural alignment overlap **90% of both**
   sequences, Linclust-style m = 300 min-hash k-mers, k = 10, ungapped diagonal pass then
   Gotoh–Smith–Waterman, cascaded prefilter → **18,661,407 clusters**; **2,302,908** non-singleton after
   fragment removal (13.05 proteins/cluster); 13,012,338 singletons. **129 h on 64 cores** [measured].
3. Define "dark" by successive exclusion: Foldseek vs PDB at **E ≤ 0.1** → 1,135,118 clusters with no PDB
   match; MMseqs2 vs Pfam at `-s 7.5 --max-seqs 100000 -e 0.001` → 883,788; then UniProt Pfam/TIGRFAM
   annotation → **711,705 (30.9%) dark clusters**.
4. Mine the dark set for remote matches to known folds. This is how bacterial gasdermin (PDB 7N51),
   archaeal/clostridial AIM2-like, BPI-like and CD4/TNFRSF4-like proteins surfaced.

Source: <https://pmc.ncbi.nlm.nih.gov/articles/PMC10584675/>.
TED adds domain-level resolution: Chainsaw + Merizo + UniDoc boundaries, CATH assignment by Foldseek
modes plus Foldclass (explicitly *not* the HMM-plus-curation CATH protocol), yielding **~365 million
domains over 188 million unique sequences, 77% similar to known superfamilies, 7,427 candidate novel-fold
PDB files** (<https://zenodo.org/record/13908086>;
<https://www.biorxiv.org/content/10.1101/2024.03.18.585509.full.pdf>).

### Pattern E — Fan out from seed structures, filter out everything sequence search already finds
*Duan 2024.*

1. Seeds: all known Acrs (100 experimentally validated) plus 46 anti-CBASS, 1,905 anti-RecBCD,
   911 anti-restriction, 10 anti-Thoeris proteins.
2. Target: IMG/VR v3 — **66,585,678 viral proteins** — plus 331,144 MAGs.
3. Reduce hard before predicting structure: length ≤ 200 aa (43,089,031 proteins) → USEARCH 95% identity
   → **15,028,625 clusters** → **BLASTP novelty filter** (drop > 30% identity or E < 0.1 to any known Acr)
   → 14,098,963 → pI > 7 discarded → CDD at E 0.01 → 8,662,777 → final catalogue **7,178,464**.
4. Predict **285,000** structures with AlphaFold2 v2.1.2 monomer, 5 models, best by mean pLDDT.
5. Structure-search with **US-align v20220511**, default parameters, **TM ≥ 0.6** (anti-defense arm also
   requires the shorter sequence ≥ 70% of the longer) → **920 candidates** with TM ≥ 0.6 to two
   different Acrs; 42 phage proteins structurally but not sequence-similar to 41 counter-defense proteins.
6. Validate: AcrVA5Bsp (Cas12a inhibition) and, via the parallel defense-island arm, BxaP.

Source: <https://pmc.ncbi.nlm.nih.gov/articles/PMC10799925/>.

**Note:** this paper used US-align, not Foldseek. Foldseek would be ~4,000× faster at comparable
sensitivity (<https://pmc.ncbi.nlm.nih.gov/articles/PMC10869269/>), which is the single clearest
cost-saving substitution available to a replay.

### Pattern F — Bait-complex screen: predict who binds the immune protein
*Yirmiya 2025.*

1. Baits: bacterial Thoeris TIR and CBASS cGAS-like proteins.
2. Screen **> 30 million phage proteins** for predicted binders whose interface covers the bait active site.
3. Keep high-confidence predictions; test in vivo and in vitro.
4. → **7 families** of Thoeris/CBASS inhibitors, thousands of genes; crystal structure of a Thoeris
   TIR–inhibitor complex (PDB 9B7D) matched the prediction; one phage TIR inhibitor blocks human and
   plant TIRs, one cGAS-like inhibitor blocks human cGAS.

Source: <https://www.weizmann.ac.il/molgen/sorek/sites/molgen.sorek/files/2025-06/Yirmiya_Cell_2025.pdf>
(DOI 10.1016/j.cell.2024.12.035). *Methods detail (pair counts, score thresholds, GPU hours) could not be
retrieved in this round — flagged as a gap.*

**What this means for the plan:** pair prediction at 3e7 scale is a GPU budget, not a CPU one, and it is
outside the 56-core envelope. The plan's contribution is the *candidate reduction* that makes the
GPU screen affordable (neighbourhood, co-occurrence, phage-encoded filter).

### Pattern G — SSN component triage with taxon-split iterative search
*Durairaj 2023.*

1. Nodes: **6,136,321** UniRef50 representatives with an AFDB v4 model at mean pLDDT > 90 (AFDB90).
2. All-vs-all MMseqs2 (release 13-45111); edge if the match covers **≥ 50% of both full-length sequences**
   at **E < 1e-4**; edge weight ∝ E-value; **max 4 outbound edges per node** (the key trick that keeps the
   graph tractable); direction ignored → 4,270,404 nodes / 10,339,158 edges → 688,852 communities →
   **242,876 components with ≥ 2 nodes**, largest holding ~50% of AFDB90; 1,865,917 unconnected.
3. Score "functional brightness" from InterPro + UniProtKB/UniParc families + disorder + coiled-coil, with
   names containing *putative / hypothetical / uncharacterized / DUF* scored zero. Dark ≤ 5%, bright ≥ 95%.
4. For a dark component: take ~94 random members, MUSCLE v5.1, then **three separate PSI-BLAST searches of
   eight rounds each over nr, split eukaryotic / archaeal / bacterial** — explicitly to stop one clade
   dominating the profile. Filter to 95% identity, cluster with CLANS at E 1e-10.
5. Context: **GCsnap v1.0.17 with a window of four flanking genes**; neighbour families by MMseqs2 at
   E < 1e-4; contexts grouped with PaCMAP + DBSCAN; keep only families present in **≥ 30% of contexts**.
6. Function hypothesis: DeepFRI v1.0.0 on up to 100 AFDB structures per cluster and context family.
7. Complex: AlphaFold-Multimer; then wet-lab (TumE–TumA toxin–antitoxin).
8. For a different component, build a Pfam family: iterative HMMER v3.3, **inclusion threshold 27 bits by
   default, adjusted manually down to recruit and up to drop matches contradicted by AlphaFold2 models**;
   confirm with HHblits, **3 iterations, inclusion E 1e-3, against UniRef30 and PDB70**; Foldseek v7
   TM-align mode and DaliLite v5 → PF22187 / DUF6946, clan CL0236.

Source: <https://pmc.ncbi.nlm.nih.gov/articles/PMC10584680>.

**Three transferable details:** the 4-outbound-edge cap; the taxon-split iterative search as a drift guard;
and the bit-score threshold adjusted against structure models rather than an E-value cut.

### Pattern H — Subtract annotated domains, cluster the residue
*Domainator 2024 case 1.*

1. Take 1,844 CARF-only proteins (from a catalogue of 6,665 CARF proteins across 13,116 genomes).
2. `hmmscan` against Pfam v36 plus curated CARF effector profiles.
3. Extract the **largest region not covered by any domain**, keep if ≥ 50 aa → **1,524 fragments**.
4. SSN the fragments → **94 clusters**; top 12 hold 1,391 of 1,524, largest holds 916.
5. AlphaFold2 (local ColabFold) on two random members per cluster.

Source: <https://pmc.ncbi.nlm.nih.gov/articles/PMC11754643> (case 1).

**What this means for the plan:** this is exactly the sub-ORF unit the plan claims as novel
(`PRIOR-ART.md` §1) — and it is already published practice. The plan's distinctive claim has to be the
*combination* with RNA features and co-occurrence statistics, not the uncovered-region clustering itself.

### Pattern I — Stream the planet, invert the index
*Serratus 2022; Logan 2024/2025.*

1. Build a small protein reference, not a big one. Serratus: **14,653 RdRP centroids at 90% aa identity,
   7.1 MB total**. Logan: 213 PAZy seeds expanded through nr and HMM filtering to 2,735,790 sequences,
   reduced to **~1.05 million non-redundant at 90% identity**.
2. **Invert query and target** so the index stays resident and the petabase-scale data streams past it.
   Logan: "The query sequence(s) were provided as indexed references to DIAMOND2, and Logan contigs were
   streamed accession by accession as queries."
3. Run DIAMOND with cheap, explicitly tuned flags. Serratus: `--mmap-target-index --target-indexed
   --masking 0 --mid-sensitive -s 1 -c1 -p1 -k1 -b 0.75` (DIAMOND v2.0.6).
   Logan: `-b 0.4 --masking 0 -s 1 --sensitive` (DIAMOND2 v2.1.9).
4. Post-filter with a domain-aware discriminator rather than an E-value: Serratus ran **palmscan `--hicon`**
   (reported FDR 0.001) to separate RdRPs from reverse transcriptases.
5. Cluster hits into taxonomic tiers (Serratus: 90% / 75% / 40% for species / genus / family-like OTUs).

Sources: <https://www.biorxiv.org/content/10.1101/2020.08.07.241729v2.full>;
<https://pmc.ncbi.nlm.nih.gov/articles/PMC12424806/>.

### Pattern J — Rule-based system census as the annotation substrate
*DefenseFinder 2024; PADLOC-DB; MacSyFinder v2.*

1. Curated profile HMMs (and covariance models, in PADLOC) per component.
2. HMM search, then completeness and co-localisation/synteny rules per system model.
3. Bit-score thresholds per profile, not global E-values (see Rousset's 32/27/27).
4. DefenseFinder grew from **60 systems at v1.0.0 to 152**
   (<https://peercommunityjournal.org/articles/10.24072/pcjournal.470/>); it now accepts nucleic-acid
   FASTA directly, calling genes with **pyrodigal v3.0.1**. PADLOC-DB v2.0.0 folds in Millman 2022
   (HMMs from DefenseFinder), Vassallo 2022, Juk, HEC and hundreds of PDC candidate systems — including
   unpublished ones (<https://github.com/padlocbio/padloc-db>).

**What this means for the plan:** `PRIOR-ART.md` already decides to store these as sidecars, never as tape
tokens. The PADLOC-DB "unpublished candidate systems" caveat is a live leakage hazard for the temporal
holdout (D-24) and should be checked per release, not once.

---

## (b) Tool performance table

Shared assumptions for all **[derived]** rows: average protein length 300 aa; query profile length
M = 300; 40 usable cores of a modern AVX2 x86 server; sequential NVMe read 1 GB/s. Where a measured
number used a different query length or core count, the normalisation is shown.

| Tool / workload | Number | Status | Hardware & database | Source URL + section |
|---|---|---|---|---|
| hmmsearch (HMMER 3.4) | `# Mc/sec: 19136.43`, `CPU time: 3.67u 0.14s 00:00:03.81 Elapsed: 00:00:01.59` for fn3 (M=149) vs 567,483 seqs / 204,940,973 residues | measured | default `--cpu 2`; UniProt/SwissProt | [HMMER 3.4 Userguide, hmmsearch tutorial statistics summary](http://eddylab.org/software/hmmer/Userguide.pdf) |
| hmmsearch per-core throughput | **8,014 Mc/s/core**; **53.8 M residues/s/core at M=149** | derived | 149 × 204,940,973 = 3.054e10 cells ÷ 3.81 s CPU = 8.01e9 cells/s; 204,940,973 ÷ 3.81 = 5.38e7 res/s | same |
| hmmsearch per-core cross-check | **8.14e9 cells/s/core** | derived | HMMER3 measured median 2 min 57 s for one iteration over UniProt 15M seqs on a single 2.9 GHz Xeon core; 300 × (15e6 × 320) = 1.44e12 cells ÷ 177 s. Agrees with the line above to 2% | [HHblits 2012, Results run times](https://zhanggroup.org/teaching/literature/Soding_2012_HHblits.pdf) |
| hmmsearch filter pass rates | MSV 0.0368 (expected 0.02), bias 0.0302, Viterbi 0.00411 (expected 0.001), Forward 0.00199 (expected 1e-5) | measured | same tutorial run — note all three stages pass **2–200× more** than their nominal rate on real data | [HMMER 3.4 Userguide](http://eddylab.org/software/hmmer/Userguide.pdf) |
| nhmmer (DNA) | `# Mc/sec: 3004.03`, CPU 0.02 s for M=80 vs 660,000 residues → **2,640 Mc/s/core** | measured + derived | 80 × 660,000 = 5.28e7 cells ÷ 0.02 s | [HMMER 3.4 Userguide, nhmmer section](http://eddylab.org/software/hmmer/Userguide.pdf) |
| HMMER3 pipeline design | MSV ~500× faster than serial full Forward; Forward/Backward ~30× slower than MSV; HMMER3 2–3× slower than NCBI BLAST, ~3× faster than WU-BLAST, up to 140× faster than HMMER2, ~600× faster than SAM | measured | one core of 2.66 GHz Intel X5550 | [Eddy 2011, Results](https://pmc.ncbi.nlm.nih.gov/articles/PMC3197634/) |
| jackhmmer defaults | max **5** iterations (`-N`); inclusion `--incE 0.001` and `--incdomE 0.001`; stops early on convergence; `-Z`/`--domZ` override the search-space size used for E-values | documented | — | [HMMER 3.4 Userguide, jackhmmer + option reference](http://eddylab.org/software/hmmer/Userguide.pdf) |
| jackhmmer per-iteration cost | "Additional iterations took roughly the same amount of time as the first iteration" | measured | HHblits paper, Supplementary Fig. 6 | [HHblits 2012](https://zhanggroup.org/teaching/literature/Soding_2012_HHblits.pdf) |
| pyhmmer vs hmmsearch | runtime reduced **72%**: "totaling only 27 h where HMMER took 97 h"; hmmscan of 1M proteins in **~2.5 h vs >7 days** estimated, 32 threads; parallel code fraction 35% → 96% | measured | 6-core i7-10710U (hmmsearch run); proGenomes2 + Pfam 33.1/35.0 | [Larralde & Zeller 2023, Results](https://pmc.ncbi.nlm.nih.gov/articles/PMC10159651/) |
| pyhmmer scaling | 2,100 proteins vs Pfam 33.1 (18,259 HMMs): **40.2 s at 11+1 threads** (v0.4.5) | measured | i7-10710U, 6c/12t, warm SSD | [pyhmmer benchmarks page](https://pyhmmer.readthedocs.io/en/stable/guide/benchmarks.html) |
| MMseqs2 prefilter memory | **M = (7 · N · L + 8 · a^k) bytes**; UniProtKB, 55M seqs × 350 aa, "needs about 71 GB of main memory"; split: `((7 · N · L)/#chunks + 21^k · 8)` | documented | — | [MMseqs2 wiki, memory requirements](https://github.com/soedinglab/MMseqs2/wiki) |
| MMseqs2 clustering memory | **M = (6 · N · r) bytes**, r = average alignment results per sequence | documented | — | same |
| MMseqs2 batch throughput | "Searching with 637000 protein sequences against 30 Mio Uniprot sequences took around **12 minutes on 16 cores**" | measured | 16 cores | same |
| MMseqs2 amortised cost | **1.1e-6 core-s per target residue** for a 637k-query batch, 1 iteration | derived | 12 min × 16 = 11,520 core-s ÷ (30e6 × 350 = 1.05e10 residues) | same |
| MMseqs2 alignment stage | "roughly **2 giga cell updates per second (GCUPS)**"; "computing 23 Mio. alignments on 8 cores takes 2 minutes"; linclust of all UniProtKB in under an hour | documented | 8 cores | same |
| MMseqs2 CPU per-query (large batch) | per-query time plateaus at **0.17 s** at 100,000 queries vs 30,430,281 sequences → **21.8 core-s/query**, i.e. **2.06e-9 core-s per target residue** | measured + derived | 2 × AMD EPYC 7742, 128 physical cores / 256 threads, 1 TB RAM; 0.17 × 128 = 21.8 core-s; ÷ (3.04e7 × 350) | [Kallenborn 2025, Results](https://pmc.ncbi.nlm.nih.gov/articles/PMC12510879/) |
| MMseqs2-GPU throughput | **11.3 TCUPS** peak on 1 L40S (18.4× the 2×64-core CPU); **67.5 TCUPS** on 8 L40S (110×); gapless synthetic 13.5 / 102 TCUPS; prior FPGA 1.7, prior K40 0.4 TCUPS | measured | L40S / A100 / H100 / L4 | same |
| MMseqs2-GPU memory | index drops from **~7 bytes/residue to 1 byte/residue**; the k-mer approach "needed up to 2 TB of RAM"; database need not fit GPU memory (streams from host RAM/disk, 64.9% of in-memory throughput) | measured | — | same |
| MMseqs2 vs jackhmmer, single batch | MMseqs2-GPU on 1 L40S **177×** faster than JackHMMER and 6.4× faster than BLAST (figure caption says ~16× — discrepancy in the paper) | measured | 6,370 queries vs 3.04e7 sequences | same |
| Iterative sensitivity (ROC1, SCOPe-derived) | 1 / 2 / 3 iterations — MMseqs2-CPU `-s 8.5` **0.391 / 0.606 / 0.665**; MMseqs2 gapless **0.400 / 0.612 / 0.669**; JackHMMER **0.400 / 0.614 / 0.685**; BLAST 0.401; PSI-BLAST 0.591 (2 it); DIAMOND ultra-sensitive 0.409 | measured | 6,370 queries vs 30,430,281 refs | same |
| MMseqs2-GPU energy / cost | 4 L40S **80.7×** more energy-efficient than JackHMMER, 2.1× vs MMseqs2-CPU; AWS c7a.32xlarge (128 cores) $6.57/h, g6e.xlarge (1 L40S) $1.86/h, g6e.48xlarge (8 L40S) $30.131/h (24 Sept 2024) | measured | — | same |
| DIAMOND v2, nr vs UniRef50 | **17.77 h** ultra-sensitive, **5.42 h** very-sensitive, on **520 nodes × 40 cores = 20,800 cores** (41,600 threads, 192 GB/node); BLASTP estimated **~2 months** on the same machine | measured | 281M nr queries vs 37.5M UniRef50 (39M in the HPC run) | [Buchfink 2021, Results](https://pmc.ncbi.nlm.nih.gov/articles/PMC8026399/) |
| DIAMOND cost per query-target pair | ultra-sensitive **2.85e10 pairs/core-h**; very-sensitive **9.35e10 pairs/core-h** (3.3× cheaper) | derived | 281e6 × 37.5e6 = 1.054e16 pairs ÷ (17.77 × 20,800 = 369,616 core-h); and ÷ (5.42 × 20,800 = 112,736) | same |
| DIAMOND memory | "**less than 16 GB of RAM**" in ultra-sensitive mode | measured | — | same |
| DIAMOND speedup / sensitivity | 80× (ultra-sensitive, BLASTP-matching) to ~8,000× (fast) vs BLASTP; 12–15× vs MMseqs2 release 11; 6–8× vs DIAMOND v0.7.12; BLASTP reference AUC1 = 0.622; seed shapes 2×w10 / 16×w8 / 14×w7 / 64×w7 | measured | 1.71M SCOPe-annotated queries | same |
| DIAMOND output volume | 23.1 billion alignments, **~1,100 GB uncompressed**, keeping top 100 hits per sequence | measured | — | same |
| DIAMOND DeepClust | **19 billion sequences in 18 days on 27 HPC systems** → 1.7 billion clusters, 32% non-singleton; 544 million clusters cover 94% of all known proteins | measured (preprint) | — | [Buchfink 2023 preprint abstract](https://www.biorxiv.org/content/10.1101/2023.01.24.525373.full.pdf) |
| HHblits per-query | 2 iterations through UniProt20 (2.6M HMMs / 15M seqs): **median 31 s, average 1 min 13 s** on a single Xeon 2.9 GHz core. One PSI-BLAST iteration over UniProt 15M: median 1 min 7 s. HMMER3: median 2 min 57 s | measured | single core | [HHblits 2012, Results](https://zhanggroup.org/teaching/literature/Soding_2012_HHblits.pdf) |
| HHblits prefilter cell rate | step 1 processes **4–5 DP cells per CPU clock cycle**, step 2 ~1.3; worked example `300 × 320 × 2.6e6 / (4.5 × 2.9 GHz) = 18 s`, "about 25% of the average time needed for the entire HHblits search" | measured | 2.9 GHz core | same, Online Methods |
| hhblits per target HMM | **2.85e-6 core-s per target HMM per iteration** | derived | HH-suite3 AVX2 is 2.1× faster than 2.0.16 at 2 iterations, so median 31 s → 14.8 core-s; ÷ (2.6e6 HMMs × 2 iterations) | measured inputs from the two HH papers |
| HH-suite3 speedups | hhblits3 **~10× faster than PSI-BLAST, ~20× faster than HMMER3**; AVX2 vs HMMER3 **14× / 20× / 29×** at 1 / 2 / 3 iterations; vs PSI-BLAST 9× / 10× / 11×; hhsearch3 3.2× (SSE2) / 4.2× (AVX2) over 2.0.16, rising to 5.0× / 7.4× for queries > 1,000 residues | measured | single core of 2 × Xeon E5-2640v3, 128 GB; uniprot20 2015_06 = 7,313,957 clusters; PDB70 = 35,000 HMMs, mean length 234 | [Steinegger 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6744700) |
| HH-suite memory | backtrace + cell-off storage **0.9 GB (SSE2) or 1.8 GB (AVX2) per thread**; HHsearch fixed overhead ~20 s per query | measured | — | same |
| HH-suite requirements | "HH-suite typically requires at least **4GB of main memory**"; estimate `query_len × max_db_seq_len × threads × (4 or 8) bytes` + the same at 8 bytes + 1 GB; keep databases on RAM disk or local SSD; recommended **1–4 iterations** (default 2); **"One iteration!"** against a custom genome database | documented | — | [HH-suite wiki](https://github.com/soedinglab/hh-suite/wiki) |
| HH custom database build | `ffindex_build` → `hhconsensus` via `ffindex_apply_mpi` → (skip `addss.pl`) → `hhmake` via `ffindex_apply_mpi` → `cstranslate_mpi -f -x 0.3 -c 4 -I a3m -i <db>_a3m -o <db>_cs219` → `ffindex_order` | documented | per-cluster build cost **not published** — must be measured locally | same |
| Foldseek search | SARS-CoV-2 RdRp (**942 residues**) vs AFDB v1 (365,198 models): **6 s** Foldseek, **33 h** TM-align, **10 d** Dali | measured | AMD EPYC 7702P 64-core, 1,024 GB RAM, 64 threads | [van Kempen 2024, Results](https://pmc.ncbi.nlm.nih.gov/articles/PMC10869269/) |
| Foldseek unit cost | **1.05e-3 core-s per target structure** for a 942-residue query; **1.12e-6 core-s per (target structure × query residue)** | derived | 6 s × 64 threads = 384 core-s ÷ 365,198 structures; ÷ 942 residues | same |
| Foldseek speed/sensitivity | >4,000× faster than TM-align and Dali and >21,000× faster than CE on SCOPe40; ~184,600× / 23,000× faster than Dali / TM-align on AFDB v1; sensitivity **86% / 88% / 133%** of Dali / TM-align / CE | measured | same | same |
| Foldseek parameters | k = 6 or 7 by database size; prefilter keeps ungapped alignments scoring ≥ **15 bits**; local alignment weights 1.4 (aa) and 2.1 (3Di); default `-s 9.5 -e 10 --max-seqs 2000`; Foldseek-TM `--alignment-mode 1 --max-seqs 4000`; E-values from a 492-parameter network, corrected by raising to the power **0.32**; score 51 ⇒ 50% true-match probability | measured | — | same |
| Foldseek memory | **(6 bytes Cα + 1 3Di byte + 1 AA byte) × database residues**; "**The 54M AFDB50 entries require 151GB**"; `--sort-by-structure-bits 0` cuts this to **35 GB** (rankings/scores change, E-values do not) | documented | — | [Foldseek README](https://github.com/steineggerlab/foldseek) |
| Foldseek 200 GB capacity | **~7.1e7 structures with Cα**, **~3.1e8 without** | derived | 151 GB ÷ 54e6 = 2.80 kB/entry → 200/0.00280 GB; 35 GB ÷ 54e6 = 0.65 kB/entry → 200/0.00065 | same |
| Foldseek server RAM (published service) | 191 GB (AFDB/UniProt50), 110 GB (ESM Atlas HQ Clu30), 3.8 GB (Proteome), 3.4 GB (Swiss-Prot), 2.0 GB (PDB100), 1.4 GB (CATH50) | measured | — | [van Kempen 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC10869269/) |
| Foldseek GPU prefilter | prefilter on an RTX 4090 is **4× faster than a 64-core CPU** | measured | needs Ampere or newer for full speed | [Foldseek README](https://github.com/steineggerlab/foldseek) |
| Foldseek-GPU (in MMseqs2-GPU paper) | **4×** faster than CPU k-mer on 1 L40S, **27.3×** on 8 GPUs; 8 GPUs 7× faster than 1 | measured | 6,370 structures from AFDB50 v4 (53.6M entries, mean 264 3Di letters) | [Kallenborn 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12510879/) |
| Foldseek-Multimer | **billions of complex pairs in 11 hours**; 3–4 orders of magnitude faster than US-align; detection 97.4% (MM-TM) / 95.8% (MM) vs 97.6% US-align at TM ≥ 0.65; chain pairing matches in > 99% of cases. For context, QSalign all-vs-all of ~100,000 complexes "took several months using 100 threads" | measured | — | [Kim 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC11903335/) |
| ProstT5 3Di prediction | *M. jannaschii* (1,787 proteins): **40 min 6 s** on a MacBook Pro M1 CPU, **43.5 s** on an RTX A5000; ColabFold **48 h on GPU** for the same proteome (extrapolated **112 days** on the same CPU); ~10 proteins/s for the human proteome on one GPU | measured | — | [Heinzinger 2024, Results](https://pmc.ncbi.nlm.nih.gov/articles/PMC11616678/) |
| ProstT5 throughput | **41 proteins/s on one RTX A5000**; **0.74 proteins/s on an M1 CPU** | derived | 1,787 ÷ 43.5 s; 1,787 ÷ 2,406 s | same |
| ProstT5 sensitivity | SCOPe40 **superfamily ROC-AUC**: predicted 3Di 0.45, encoder+CNN 3Di **0.47**, experimental-structure 3Di 0.49, MMseqs2 sequence search **0.06** | measured | Cα refinement disabled for predictions | same |
| ProstT5 vs structure prediction | "about **400-4000x**" faster than predicted structures via ColabFold; one to two further orders of magnitude from GPUs | documented | — | [Foldseek README](https://github.com/steineggerlab/foldseek) |
| Foldseek structural clustering | 52M AFDB50 structures clustered in **129 h on 64 cores** ≈ **8,256 core-hours**, i.e. **0.57 core-s per structure** | measured + derived | 129 × 64; 8,256 × 3600 ÷ 52.3e6 | [Barrio-Hernandez 2023, Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC10584675/) |
| Infernal 1.1 cmsearch | **1.5 s/Mb/query** on the benchmark, **0.6 s/Mb/query** on a several-gigabase 15-genome database; benchmark completed in **0.44 h** vs 49.31 h (v1.0.2) and 4,359 h unfiltered (≈ 100× and 10,000×); HMMER3 0.02 h and BLASTN 0.01 h but "significantly less sensitive"; HMM filter stages take ~1/3 of runtime | measured | — | [Nawrocki & Eddy 2013, Results](https://pmc.ncbi.nlm.nih.gov/articles/PMC3810854/) |
| cmsearch example | 1 CM of 72 consensus positions, 5,874,406 residues searched (both strands): `CPU time: 1.44u 0.03s 00:00:01.47 Elapsed: 00:00:02.91`, "about 3 seconds (wall clock) running on 4 cores" | measured | — | [Infernal 1.1.5 Userguide, cmsearch tutorial](http://eddylab.org/infernal/Userguide.pdf) |
| cmsearch filter design | filter thresholds **scale with search-space size** ("This differs from HMMER3 which always uses the same filter thresholds") | documented | — | same |
| cmcalibrate memory | eukaryotic LSU rRNA model RF02543 "requires about **14Gb of RAM per thread** and roughly **30 hours on 12 threads (so 168Gb RAM total)**"; `--split`/`--merge` added in 1.1.5 to partition without MPI | documented | — | same, "Parallelizing calibration of large models" |
| ColabFold MSA | 3-iteration MMseqs2 profile search, **40–60× faster** than the AlphaFold2 jackhmmer+HHblits pipeline, which "can take up to hours for a single protein" and needs >2 TB of storage; 1,762-protein proteome in **113 min**; MMseqs2 max resident set **308 GB** | measured | AMD EPYC 7402P 24-core, 512 GB; server 2×14-core E5-2680v4, 768 GB | [Mirdita 2022, Methods 2.2/2.7.4](https://pmc.ncbi.nlm.nih.gov/articles/PMC9184281/) |
| ColabFold MSA per protein | **3.85 s wall, ~92 core-s per protein** for a 3-iteration search over UniRef30 + ColabFoldDB | derived | 113 min ÷ 1,762; × 24 cores | same |
| ColabFold database sizes | UniRef30 = 29.3M cluster consensus sequences standing in for 277.5M UniRef100; BFD+MGnify reduced 2.5 billion → 513 million, "requiring only 84 GB RAM for headers and sequences"; ColabFoldDB = 209,335,865 representatives / 738,695,580 members | measured | E-value cutoff 0.1; diversity filters `--max-seq-id 0.95`, `--filter-min-enable 1000 --diff 3000` | same |
| PLMSearch (pLM retrieval) | SCOPe40-test all-vs-all (2,207 proteins, 4,870,849 pairs): **PLMSearch 4 s, SS-predictor 10 s, TM-align 11,303 s**; one web-server query **~0.15 min vs Swiss-Prot (568K)** and **~1.1 min vs UniRef50 (53.6M)**; remote-homology recall over 1,105 remote pairs: MMseqs2 183, BLASTp 203, Foldseek 934, Foldseek-TM 940, **PLMSearch 1,087**; AUROC MMseqs2 → PLMSearch: family 0.318 → 0.928, superfamily 0.050 → 0.826, fold 0.002 → 0.438 | measured | 56-core Xeon E5-2680 v4, 256 GB (benchmarks); 64 × Xeon E5-2682 v4, 512 GB, **CPU only** (server) | [Liu 2024, Results](https://pmc.ncbi.nlm.nih.gov/articles/PMC10981738/) |
| DHR (dense homolog retrieval) | ">10% increase in sensitivity" overall, ">56% ... at the superfamily level"; "up to 22 times faster than ... PSI-BLAST and DIAMOND and up to 28,700 times faster than HMMER"; used as a **prefilter feeding JackHMMER** for MSA construction | institutional summary, not verified against Methods | — | [CUHK summary of Nat Biotechnol 2024 (10.1038/s41587-024-02353-6)](https://www.cse.cuhk.edu.hk/news/achievements/prof-yu-lis-group-develops-a-ground-breaking-method-dense-homolog-retrieval-dhr-for-the-rapid-and-sensitive-detection-of-remote-protein-homologs-from-extensive-databases/) |
| Logan translated search | DIAMOND2 v2.1.9, `-b 0.4 --masking 0 -s 1 --sensitive`, "each take approximately **11 hours on 60,000 vCPUs**"; contigs 0.90 Pbp (0.31 PB compressed) from 26.8M accessions; unitigs 4.59 Pbp (2.18 PB compressed) | measured | AWS/Azure cloud | [Logan, Methods 3.2/3.5](https://pmc.ncbi.nlm.nih.gov/articles/PMC12424806/) |
| Logan unit cost | **~660,000 vCPU-hours**, i.e. **0.733 vCPU-h per Gbp of contigs** for a ~1.05M-sequence protein reference | derived | 11 h × 60,000; ÷ 9.0e5 Gbp | same |
| Logan-Search | k = 31; **~1 PB on-disk index** of ~2e15 k-mers from 23.4M accessions; **~6 min per 1 kb query on 12 small 4-vCPU VMs**; built in ~10 days on 625 Azure Standard_D32d_v5 VMs (~45,000 tasks) | measured | — | same, Methods 3.3 |
| Logan assembly | ~30 h wall clock, peak **2.18 million CPU cores**, ~**30 million CPU hours** | measured | — | same |
| Serratus | 5,686,715 SRA runs (**10.2 petabases**) in **11 days**, peak **22,250 vCPUs**; **$0.0062/run** nucleotide (bowtie2) and **$0.0042/run** translated (DIAMOND); total **$23,980**; 883,502 RdRP contigs → 140,208 sOTUs | measured | ~4,300 c5.xlarge spot (align), ~1,300 r5.xlarge (download), ~60 c5.large (merge) | [Serratus preprint, Methods](https://www.biorxiv.org/content/10.1101/2020.08.07.241729v2.full) |
| Serratus cost per petabase | **$2,351 per Pbp** = **$2.35 per Tbp** for a 7.1 MB (14,653-sequence) protein reference | derived | $23,980 ÷ 10.2 Pbp | same |
| FLSHclust corpus & benchmark | 8.8 Tbp, **8 billion proteins**, 10.2M arrays; 1.3e9 / 646.4M / 499.9M clusters at 90 / 50 / 30% identity; mean 2.0 non-redundant proteins per 30% cluster; benchmarks on a 2-node grid, 64 CPUs, 416 GB RAM, 2 TB SSD per node; MMseqs2 on full UniRef50 (51M proteins) "required substantially more compute resources to complete within a week" | measured | — | [Altae-Tran 2023, Methods summary](https://pmc.ncbi.nlm.nih.gov/articles/PMC10910872) |
| GCsnap2 Cluster | **22× improvement in execution time** over GCsnap1 Desktop; "hundreds of thousands of input sequences in HPC clusters" via `mpi4py.futures` | measured | — | [arXiv 2505.02195](https://arxiv.org/abs/2505.02195) |
| BLAST E-value statistics | `E = K · m · n · exp(−λS)`, with *n* the total length of the database: **E scales linearly with database size at fixed score**; composition-based statistics apply a per-database-sequence adjustment so "identical alignments can receive different scores"; conditional compositional score matrix adjustment (Yu & Altschul 2005) is the current NCBI protein BLAST default | documented | — | [NCBI BLAST documentation](https://https.ncbi.nlm.nih.gov/books/NBK62051/) |
| Profile drift (PSI-BLAST) | "**5–9% of alignments at iteration 5 are non-homologous**" in full-protein searches; homologous over-extension (HOE) occurs **1.5–3× as often** as non-homologous error at iteration 5; boundary masking improves specificity **4–8×** with little sensitivity loss (hard families: 74.6% coverage at 1.4% FP vs 77.3% at 5% FP) | measured | Pfam-derived libraries with/without partial domains | [Gonzalez & Pearson 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2853128/) |
| Profile drift (query seeding) | query seeding cuts false-positive alignment frequency **>20-fold vs psiblast and jackhmmer** and FDR **5–20×**; on far50 after 10 iterations sensitivity fell 93.9% (jackhmmer) → 86.9% while **FDR fell 13.9% → 0.2%**; median identity of bottom-quartile alignments rose 13.3% → 19.0% | measured | RPD3 library (~598,000 proteins) from Pfam 28 + RefProtDom; 100 far50 + 100 far66 queries | [Pearson, Li & Lopez 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC5605230/) |
| PSI-Search vs jackhmmer | PSI-Search produces **2–4× fewer false positives than JackHMMER** but is **~5% less sensitive** | measured | — | [Li 2012 PSI-Search](https://pmc.ncbi.nlm.nih.gov/articles/PMC3371869) |

---

## (c) Cost per query on 40 cores at 1e7, 1e8, 1e9 proteins and 3e12 nt

**Assumptions.** Average protein 300 aa, so 1e7 / 1e8 / 1e9 proteins = 3e9 / 3e10 / 3e11 residues.
Query profile M = 300. Forty cores comparable to the 2.25–2.9 GHz AVX2 server cores in the measured
sources. NVMe sequential read 1 GB/s. Six-frame translation of 3e12 nt yields **6e12 aa residues**
(6 frames × 3e12/3 codons), i.e. 20× the residue count of the 1e9-protein set.

Every cell is **[derived]** from a **[measured]** constant in §b. The "binding constraint" column is the
point of the exercise: at 1e9, five of the eight workloads are limited by memory or I/O, not arithmetic.

### c.1 hmmsearch, one profile

| targets | cells | core-s | 40-core wall (compute) | I/O floor | binding constraint |
|---|---|---|---|---|---|
| 1e7 (3e9 res) | 300 × 3e9 = 9.0e11 | 9.0e11 ÷ 8.01e9 = **112** | **2.8 s** | 3 GB → 3 s | compute, cache-resident |
| 1e8 (3e10 res) | 9.0e12 | **1,124** | **28 s** | 30 GB → 30 s | balanced |
| 1e9 (3e11 res) | 9.0e13 | **11,236** (3.12 core-h) | **281 s** | 300 GB → 300 s | **I/O**; realistic wall 5–10 min |

Constant: **8.01e9 cells/s/core** from `# Mc/sec: 19136.43` ÷ 2.39 effective cores
(3.81 s CPU / 1.59 s elapsed), independently reproduced by the HHblits-paper HMMER3 timing to 2%.

Scale-out note: a full Pfam-class scan is ~20,000 profiles. 20,000 × 3.12 core-h = **62,400 core-h =
65 days on 40 cores** for 1e9 proteins. That is the whole argument for precomputing family hits once per
release and for using pyhmmer's threading model (measured 72% runtime reduction, 97 h → 27 h).

### c.2 jackhmmer, 5 iterations, one query

Per-iteration cost ≈ one hmmsearch at M = query length, because "additional iterations took roughly the
same amount of time as the first iteration" (HHblits 2012).

| targets | core-s | 40-core wall (compute) | + 5× I/O | binding constraint |
|---|---|---|---|---|
| 1e7 | 5 × 112 = **560** | **14 s** | +15 s | compute |
| 1e8 | 5 × 1,124 = **5,620** (1.56 core-h) | **141 s** | +150 s | balanced |
| 1e9 | 5 × 11,236 = **56,180** (15.6 core-h) | **23 min** | +25 min → **~50 min** | **I/O + profile memory** |

Cross-check: the HHblits paper measured HMMER3 at median 2 min 57 s per iteration over UniProt 15M seqs
(4.8e9 residues) on one core. Five iterations = 885 core-s; scaled to 3e11 residues (× 62.5) =
**55,300 core-s** — within 2% of 56,180.

Unmodelled cost at 1e9: a permissive query can include >1e6 sequences in the round-2 alignment, and
profile construction plus the stored alignment then dominate memory. Mitigations with published support:
tighten `--incE` below the 0.001 default; run against a 50%-identity clustered subset (FLSHclust measured
8e9 proteins → 646.4M clusters at 50%, a **12.4× reduction**, cutting the 1e9 case to ~4.5 core-h);
or use the DHR pattern — retrieve a candidate subset, then jackhmmer only that.

### c.3 MMseqs2 3-iteration profile search

Constant: **2.06e-9 core-s per target residue** per iteration, from the measured 0.17 s/query at
100,000 queries on 128 cores against 30,430,281 sequences (0.17 × 128 = 21.8 core-s ÷ 1.06e10 residues).

**Single query** (the interactive case):

| targets | compute core-s (3 it) | 40-core compute | default index (7 B/res) | gapless index (1 B/res) | I/O-bound wall, 3 it |
|---|---|---|---|---|---|
| 1e7 | 3 × 6.2 = **18.5** | **0.5 s** | 21 GB — fits | 3 GB | ~2 s |
| 1e8 | 3 × 62 = **185** | **4.6 s** | 210 GB — **over cap**, 2 splits | 30 GB | ~30 s (gapless) |
| 1e9 | 3 × 618 = **1,854** (0.52 core-h) | **46 s** | **2.1 TB** — 11+ splits | **300 GB** — 2 splits | **~1.75 h** default; **~15 min** gapless |

The 1e9 default-index figure is `2.1 TB ÷ 1 GB/s × 3 iterations = 6,300 s`. With `--prefilter-mode 1`
(measured 1 byte/residue in the MMseqs2-GPU paper) it is `300 GB × 3 = 900 s`. **Compute is 46 s. Storage
is 15 minutes to 1.75 hours.** At this scale MMseqs2 is a storage-layout problem.

**Batched** (the right way to run it):
amortised **1.1e-6 core-s per target residue** for a 637k-query batch, from "637000 ... against 30 Mio ...
around 12 minutes on 16 cores". Scaling to 3e11 residues: `11,520 core-s × (3e11 ÷ 1.05e10) = 329,000
core-s = 91.4 core-h` per iteration for 637,000 queries → **2.3 h on 40 cores per iteration, 6.9 h for
three**. Per query that is **0.52 core-s per iteration — about 1,200× cheaper than the single-query
figure.** Batching queries is the single largest lever in this whole document.

### c.4 DIAMOND, one query

| targets | ultra-sensitive core-s | very-sensitive core-s | conservative (measured-benchmark) core-s | 40-core wall (conservative) |
|---|---|---|---|---|
| 1e7 | 1e7 ÷ 2.85e10 core-h = **1.26** | 0.39 | 18 | 0.5 s |
| 1e8 | **12.6** | 3.9 | 177 | 4.4 s |
| 1e9 | **126** | 38.5 | **1,770** | **44 s** |

Two routes disagree by 14×, and both are reported honestly here. The DIAMOND v2 HPC run gives
**2.85e10 pairs/core-h** (ultra-sensitive) from 1.054e16 pairs ÷ 369,616 core-h. The MMseqs2-GPU
benchmark measured DIAMOND at 0.42 s/query on 128 cores against 3.04e7 sequences = **53.8 core-s/query**,
i.e. **1.77e-6 core-s per target sequence**. The modes, block sizes (`-b`) and per-node memory differ.
**Plan with the conservative column.**

DIAMOND is the only tool in this table with no memory problem at 1e9: "**less than 16 GB of RAM**" in
ultra-sensitive mode. Under a 200 GB cap that makes it the natural first-pass tool, with MMseqs2 or
jackhmmer applied to the reduced candidate set. Its sensitivity cost is small: DIAMOND ultra-sensitive
ROC1 **0.409** vs BLAST 0.401 and JackHMMER 0.400 at one iteration.

### c.5 hhblits against a custom DB of cluster HMMs

Cluster count from the FLSHclust ratio (8e9 proteins → 499.9M clusters at 30% identity = **1 per 16**):
1e7 → 6.3e5 HMMs; 1e8 → 6.3e6; 1e9 → 6.3e7.
Search constant: **2.85e-6 core-s per target HMM per iteration** (median 31 s for 2 iterations over 2.6e6
HMMs on one 2.9 GHz core, ÷ 2.1 for HH-suite3 AVX2 = 14.8 core-s, ÷ 5.2e6 HMM-iterations).

| targets | cluster HMMs | core-s (3 it) | 40-core wall | storage / memory | binding constraint |
|---|---|---|---|---|---|
| 1e7 | 6.3e5 | **5.4** | **0.13 s** | cs219 ~0.2 GB; a3m ~GB | trivial |
| 1e8 | 6.3e6 | **54** | **1.3 s** | cs219 ~1.9 GB; a3m ~tens of GB | build cost |
| 1e9 | 6.3e7 | **539** | **13 s** | cs219 ~19 GB; a3m proportional to 3e11 alignment residues → **hundreds of GB** | **build cost + a3m storage** |

Per-thread memory is also non-trivial: 1.8 GB (AVX2) × 40 threads = **72 GB**, over a third of the cap,
before any database is loaded. HH-suite's own guidance is `-cpu 1 -maxmem 1` to reduce this, and to keep
the database on RAM disk or local SSD.

**Honest gap:** the search arithmetic says corpus-wide hhblits is seconds. The *build* (one
`hhconsensus` + `hhmake` + `cstranslate` per cluster, via `ffindex_apply_mpi`) has **no published
per-cluster timing** in any source checked, so I will not invent one — it must be measured locally on
1e4 clusters and extrapolated. The published practice is instructive: Durairaj built PF22187 by hand and
verified it with HHblits against the *stock* UniRef30 and PDB70 databases, not against a custom corpus
database. Build the HH database for the 1e5–1e6 clusters that matter, not for 6.3e7.

### c.6 Foldseek search against predicted structures

Constant: **1.12e-6 core-s per (target structure × query residue)** → **3.36e-4 core-s per target
structure** for a 300-residue query.

| targets | core-s | 40-core wall | memory with Cα (8 B/res) | memory without Cα (~1.85 B/res) | binding constraint |
|---|---|---|---|---|---|
| 1e7 | **3,360** (0.93 core-h) | **84 s** | 24 GB | 5.6 GB | compute |
| 1e8 | **33,600** (9.3 core-h) | **14 min** | 240 GB — **over cap** | 56 GB | memory unless `--sort-by-structure-bits 0` |
| 1e9 | **336,000** (93 core-h) | **2.3 h** | 2.4 TB | 560 GB | **memory — infeasible as monomers** |

Derived 200 GB capacity: **~7.1e7 structures with Cα** (151 GB ÷ 54e6 = 2.80 kB/entry) or **~3.1e8
without** (35 GB ÷ 54e6 = 0.65 kB/entry). So the corpus-scale configuration is forced: run Foldseek
against **50%-identity cluster representatives**, which for 1e9 unique proteins is ~8e7 representatives
(1 per 12.4) = 224 GB with Cα or **52 GB without** — comfortably inside the cap with
`--sort-by-structure-bits 0`. AFDB50 at 54M entries / 151 GB is the published precedent.

**The real cost is generating the structural alphabet, not searching it.** ProstT5 encoder+CNN measured at
41 proteins/s on one RTX A5000: 8e7 representatives = **545 GPU-hours**; all 1e9 proteins =
**6,800 GPU-hours**. ColabFold measured at 48 h for 1,787 proteins = 96.7 GPU-s/protein: 8e7
representatives = **2.1e6 GPU-hours** — a ~4,000× difference, matching Foldseek's own "400-4000x" claim.
Sensitivity cost of using ProstT5 instead of structures is small at the level that matters here:
SCOPe40 superfamily ROC-AUC **0.47 (ProstT5 CNN) vs 0.49 (experimental 3Di) vs 0.06 (MMseqs2)**.
**Decision: ProstT5 3Di for the corpus; real or predicted 3D only for curated candidates.**

### c.7 cmsearch on IGRs

Constant: **0.6 s/Mb/query** on a multi-gigabase database (the right regime), **1.5 s/Mb/query** on the
benchmark. IGR fraction is a *parameter* — 12% is used below as a placeholder the lab must measure on its
own corpus, not a cited figure.

| target | Mb | core-s per CM @0.6 | 40-core wall per CM | @1.5 s/Mb |
|---|---|---|---|---|
| full 3e12 nt | 3.0e6 | **1.8e6** (500 core-h) | **12.5 h** | 31 h |
| 12% IGR slice (3.6e11 nt) | 3.6e5 | **2.16e5** (60 core-h) | **1.5 h** | 3.75 h |

Full-Rfam scale (~4,100 families): 4,100 × 60 core-h = **246,000 core-h = 256 days on 40 cores** over the
IGR slice; **5.8 years** over the full 3e12 nt. One CM over the IGR slice (1.5 h) is an interactive query;
one CM over everything (12.5 h) is an overnight job; all of Rfam over everything is a rented-bare-metal
job run once per release.

**The binding constraint is calibration, not search.** `cmcalibrate` for RF02543 needs **14 GB per thread
and 168 GB total on 12 threads**; 40 threads would need ~560 GB, far over the cap. Either take calibrated
CMs from Rfam, or use Infernal 1.1.5's `--split`/`--merge`. Note also that cmsearch's filter thresholds
*scale with search-space size*, unlike HMMER3's fixed thresholds — so cmsearch results on a 3.6e11 nt IGR
set are not directly comparable to the same CM run on one genome.

### c.8 tblastn / MMseqs2 translated search on contigs (3e12 nt)

Six-frame residue count: **6e12 aa**.

**MMseqs2 TBLASTN-like** (nucleotide target database; `--search-type 2` for TBLASTX-like):

| configuration | compute core-s (3 it) | index size | splits under 200 GB | I/O-bound 40-core wall |
|---|---|---|---|---|
| all 6e12 aa, default index | 3 × 12,360 = **37,080** (10.3 core-h) | **42 TB** | ~210 | infeasible |
| all 6e12 aa, gapless index | same | **6 TB** | 30 | 6,000 s/it × 3 = **5 h** |
| IGRs + contig ends only (~7.5e11 aa) | 3 × 1,545 = **4,635** | **750 GB** | 4 | 750 s/it × 3 = **~40 min** |

Compute is 10 core-hours. Storage is 42 TB or 6 TB. The plan already stores called proteins, so
six-frame search exists only to catch gene-call misses — restrict it to IGRs and contig ends and it
becomes a 40-minute job instead of an impossible one.

**DIAMOND blastx over contigs** (the measured route): **0.733 vCPU-h per Gbp** → 3,000 Gbp =
**2,200 vCPU-h = 55 h on 40 cores** for a ~1.05M-sequence protein reference, with the query/target
inversion Logan used. Cloud cost anchors: Serratus's **$2.35 per Tbp** (derived: $23,980 ÷ 10.2 Pbp) for a
7.1 MB reference → **~$7 for 3 Tbp**; Logan's 1.05M-sequence reference is far larger and reportedly cost
~$10k for 0.90 Pbp (talk-slide figure, not in the paper) → ~$33 scaled to 3 Tbp. For a *single* protein
query the cost is dominated by streaming 3 TB once (~1 h at 1 GB/s), so single-query and
hundred-query runs cost almost the same. **Batch.**

**tblastn proper**: no modern published throughput at this scale exists in anything checked. The right
order of magnitude comes from DIAMOND's measured 80× speedup over BLASTP at matched sensitivity and the
paper's estimate that BLASTP would need ~2 months on 20,800 cores for nr-vs-UniRef50. tblastn over
3e12 nt on 40 cores is a multi-month job. **Do not plan for it.**

### c.9 Summary: what fits in the envelope

| workload, 1e9 proteins or 3e12 nt | 40-core wall | peak RAM | verdict under 200 GB |
|---|---|---|---|
| hmmsearch, 1 profile | 5–10 min | < 5 GB | interactive |
| jackhmmer, 5 iterations | ~50 min | GB-scale, grows with hit count | overnight; prefer clustered target |
| MMseqs2 3-it profile, single query | 15 min (gapless) – 1.75 h (default) | 200 GB in splits | fits with `--prefilter-mode 1` |
| MMseqs2 3-it profile, 637k-query batch | 6.9 h | as above | the right way to run it |
| DIAMOND, 1 query | 44 s | **< 16 GB** | interactive; best first pass |
| hhblits vs 6.3e7 cluster HMMs | 13 s search | 72 GB threads + DB | search trivial, **build unmeasured** |
| Foldseek vs 8e7 cluster reps (ProstT5 3Di) | ~11 min | 52 GB (no Cα) | fits; 545 GPU-h to build |
| cmsearch, 1 CM over 12% IGRs | 1.5 h | < 10 GB | interactive; **calibration needs 168 GB** |
| MMseqs2 translated, IGRs + edges | ~40 min | 200 GB in 4 splits | fits |
| DIAMOND blastx, all contigs | 55 h | < 16 GB | weekend job, or $7–33 on spot |

---

## (d) What modern groups keep local, and what they leave on web services

The pattern across every 2023–2026 paper read here is the same: **keep what the whole corpus is needed
for; rent what a human looks at once.**

**Kept local, always.**
- *The clustering.* Every discovery paper builds its own clusters (FLSHclust LSH at 30%; MMseqs2 50/90 in
  Barrio-Hernandez; USEARCH 95% in Duan; MMseqs2 `-c 0.6` in Rousset). No service offers clustering of
  your corpus, and the cluster is the unit the statistics run on.
- *The neighbourhood extraction.* 3000 bp windows, ±4 ORFs, framework-gene bracketing — all local, all
  cheap, all dependent on having coordinates.
- *The enrichment statistic.* FLSHclust's N_eff-normalised score, Rousset's arrangement counts,
  Durairaj's "families present in ≥30% of contexts". This is where the discovery happens, and it cannot
  be outsourced because it depends on the corpus's own null.
- *First-pass homology at corpus scale.* DIAMOND (< 16 GB RAM) or MMseqs2 with a gapless prefilter.
- *Rule-based system calls.* DefenseFinder, PADLOC, MacSyFinder, run locally over everything, stored as
  sidecars.

**Left to web services.**
- *HHpred / HHblits against stock databases.* Durairaj verified PF22187 with HHblits against the stock
  UniRef30 and PDB70 — not a custom corpus HH database. Rousset annotated hotspot cargo with HHpred at a
  **80% probability** threshold. Both are per-candidate, dozens-of-queries workloads. The reason is in
  §c.5: HH databases are expensive to build and cheap to search, so you only build them for something you
  will query thousands of times.
- *Foldseek server.* Published RAM for the service instance is **191 GB for AFDB/UniProt50 and 110 GB for
  ESM Atlas HQ Clu30** — i.e. the server exists precisely because that does not fit on a lab box
  alongside anything else. Use it for one-off candidate triage; keep a local AFDB50-sized or
  ProstT5-3Di-of-cluster-representatives database for systematic work.
- *AlphaFold2 / AlphaFold-Multimer.* Everyone uses GPUs they do not own, or Colab. The CPU-only MSA stage
  is the part worth localising, and ColabFold's measured 40–60× speedup is exactly that substitution.
- *EFI-EST / EFI-GNT.* Per-query SSN and genome-neighbourhood networks, with published sequence caps
  (release notes describe limits in the 10^4–3×10^5 range depending on input type,
  <https://efi.igb.illinois.edu/efi-est/notes.php>). Good for one enzyme family, structurally incapable of
  corpus scale — which is the gap `PRIOR-ART.md` identifies.
- *Nucleotide search of the SRA.* Logan-Search (**~6 min per 1 kb query on 12 small VMs**, 1 PB index) and
  MetaGraph. `PRIOR-ART.md` already rules this out of scope; the cost numbers confirm it — the index
  alone is 1 PB and took ~10 days on 625 VMs to build.
- *Translated search of the SRA.* Logan's measured **11 h on 60,000 vCPUs** per sweep is a
  compute-to-data, rent-it job. Serratus's **$0.0042 per run** is the number to quote when someone asks.

**The decision rule this suggests:** localise anything whose cost is linear in the corpus and whose output
feeds the statistics; rent anything whose cost is linear in the number of *candidates* a human will look
at. The crossover sits at roughly 10^4 queries: below that, services win on total cost of ownership;
above it, the local corpus wins because you are paying per query to a service and per release to yourself.

---

## (e) Failure modes and how modern papers guard against them

**1. Profile drift and homologous over-extension.** The best-quantified failure in the whole field.
Gonzalez & Pearson measured **5–9% of alignments non-homologous at PSI-BLAST iteration 5** in full-protein
searches, with HOE occurring **1.5–3× as often** as outright non-homologous error, and showed HOE drops
substantially when partial domains are removed from the library. Pearson, Li & Lopez then measured
jackhmmer's FDR at **13.9% after 10 iterations on far50**, reduced to **0.2%** by query seeding, at a cost
of 93.9% → 86.9% sensitivity.
*Guards actually used:* boundary masking / query seeding (4–8× and >20× specificity gains measured);
taxon-split iteration (Durairaj ran **three separate PSI-BLAST searches over eukaryotic, archaeal and
bacterial nr** rather than one over all of nr); bit-score inclusion thresholds adjusted against structure
models rather than E-value cuts (Durairaj's hmmsearch 27-bit default "lowered to find more homologues and
raised to remove false matches flagged by AlphaFold2 models"); and capping iterations at 2–4 (HH-suite's
own recommendation, versus jackhmmer's default 5).
*For the plan:* the Pearson RPD3 library (~598,000 proteins, far50/far66 query sets) is a ready-made,
published harness for measuring drift on our own stack. Use it as an acceptance test, not an assumption.

**2. E-values that are not comparable across database sizes.** `E = K·m·n·exp(−λS)` makes E scale
linearly with *n*, the total database length. A hit at E = 1e-10 against 1e7 proteins and one at 1e-10
against 1e9 proteins differ by two orders of magnitude in evidential weight. Infernal makes this worse in
a different direction: its **filter thresholds themselves scale with search-space size**, "unlike HMMER3
which always uses the same filter thresholds", so the *recall* of a cmsearch changes with database size,
not just the reported E-value.
*Guards:* HMMER's `-Z` and `--domZ` exist for exactly this — assert a fixed search space so E-values stay
comparable across releases and across subsets. Rousset sidestepped E-values entirely for the system
census, using **per-profile bit-score thresholds (32 / 27 / 27)**. Bit scores are the portable quantity;
E-values are not.
*For the plan:* fix `-Z` to a declared constant per release, record it in the ledger, and store bit
scores alongside E-values for every hit. Otherwise a stage-5 to stage-6 corpus growth silently moves every
threshold in the system.

**3. Clonal and taxonomic inflation.** A family present in 10,000 near-identical *E. coli* genomes looks
10,000× more "conserved" than one in 10 diverse archaea.
*Guards:* FLSHclust divides its score by **N_eff, an effective sample size**, and requires **N_eff ≥ 3** —
the single most important line in that paper's methods. Rousset counted **arrangements** (distinct cargo
configurations), not loci, and required ≥ 5 or ≥ 10 occurrences. Durairaj capped **outbound edges per node
at 4** so a hyper-represented family cannot dominate the SSN.
*For the plan:* the co-occurrence statistics already plan for phylogeny correction; the missing piece is
an N_eff-equivalent applied at the *cluster* level, before any enrichment test.

**4. Composition bias.** Composition-based statistics apply a per-database-sequence scoring adjustment, so
"identical alignments can receive different scores, based upon the compositions of the sequences they
involve". HMMER3 runs a dedicated **bias filter** between MSV and Viterbi, and a **null2** correction in
scoring.
*Guards:* leave the bias filter on; Serratus's DIAMOND run used `--masking 0` for speed but then applied
**palmscan** as a domain-aware discriminator with a reported **FDR of 0.001** rather than trusting the
alignment score. Duan discarded candidates with **pI > 7** — a crude but effective composition filter for
the Acr search.
*For the plan:* never accept a hit on score alone for low-complexity-prone families (coiled-coil,
TM-rich, repeat proteins); require a domain-model or motif confirmation step.

**5. Filter pass rates that are nothing like their nominal values.** The HMMER 3.4 tutorial's own example
shows MSV passing **0.0368 of targets where 0.02 is expected**, Viterbi **0.00411 vs 0.001**, and Forward
**0.00199 vs 1e-5** — the Forward filter passes **200× more than nominal** on real data. Runtime and memory
estimates built on nominal pass rates will be wrong by orders of magnitude for real, correlated databases.
*Guard:* benchmark on the actual corpus, not on shuffled sequences. (Eddy's own speed figures used
"randomized (shuffled) target sequences" for one panel and real TrEMBL for another, for this reason.)

**6. Promiscuous domains and false fusions.** A protein carrying AAA+, TIR, HEPN or a DUF that appears in
hundreds of unrelated architectures will link everything to everything in a fusion or neighbourhood graph.
*Guards:* Domainator's CARF recipe inverts the problem — it **subtracts all annotated domains** and
clusters only the uncovered residue, so promiscuous domains cannot drive the clustering. Rousset required
custom fusion HMMs with high thresholds (**80** for AAA_15+DUF4435) distinct from the component thresholds
(32, 27). Durairaj's brightness score assigns **zero coverage** to names containing "putative",
"hypothetical", "uncharacterized" or "DUF", so a DUF annotation cannot make a dark family look bright.
*For the plan:* maintain an explicit promiscuity statistic per domain model (number of distinct
architectures, number of distinct neighbour families) and down-weight or exclude the top tail in any
fusion or guilt-by-association inference. This is cheap to compute from the tape and is the guard most
likely to be missing.

**7. Gene-call errors creating phantom fusions and splits.** Six-frame and nucleotide-level evidence is
the only check. DefenseFinder 2024 now calls genes itself with **pyrodigal v3.0.1** rather than trusting
submitted annotation. Serratus assembled micro-contigs and ran palmscan on those, not on reads.
*For the plan:* this is why the §c.8 IGR-and-contig-edge translated search matters even though the corpus
stores called proteins — it is the audit, not the primary search.

**8. Benchmarks that do not measure what you need.** The SCOPe-derived ROC1 numbers in §b are tight —
JackHMMER 0.685 vs MMseqs2 gapless 0.669 at three iterations — but they are measured on 6,370 queries
against 30M sequences with reversed-sequence decoys. That says little about a 1e9-protein search where the
decoys are real proteins with real compositional structure, and nothing at all about neighbourhood
recall. PLMSearch's remote-homology recall framing is closer to useful (**1,087 of 1,105** remote pairs vs
MMseqs2's 183), but it is a 2,207-protein set.
*Guard used by the better papers:* benchmark against a published hit list from a prior study, at the
stage where the target taxa first enter. That is exactly the replay test `00_brief.md` specifies, and §f
lists the candidates.

**9. Candidate-set leakage into temporal holdouts.** PADLOC-DB v2.0.0 includes "several Hma–Embedded
Candidate (HEC) systems and hundreds of additional Phage Defence Candidate (PDC) systems" from
**unpublished work**. Any temporal-holdout evaluation that uses a PADLOC release newer than the holdout
date is contaminated.
*Guard:* pin the exact DefenseFinder and PADLOC-DB versions per release in the ledger, and check the
candidate-system content of each new version rather than assuming it only adds published systems.

**10. Reported numbers that disagree with themselves.** Worth stating because it recurs: the Logan paper
reports "1.12 billion matching sequences" in the main text and "an additional 6.5 billion hits" at
E < 1e-8 in the Methods; the MMseqs2-GPU paper says "6.4× faster than BLAST" in the text and "~16× faster
than BLAST" in the figure caption; Barrio-Hernandez's Methods cite Pfam v35.0 while Data availability says
v34.0; DIAMOND's UniRef50 is 37.5M in Methods and 39M in the HPC run.
*Guard:* for any replay, derive the number from the deposited artifact, not the prose.

---

## (f) Best external replay candidates

Ranked by how much of the search is actually reconstructable from public artifacts.

**Tier 1 — complete, end-to-end replays.**

1. **Domainator 2024** — `https://github.com/nebiolabs/domainator_examples` plus Zenodo
   `10.5281/zenodo.10989173` (also record 14056380) carry the **full command-line workflows for all three
   cases**, plus a Colab notebook. This is the single best replay target in the external literature: it
   exercises sub-ORF domain subtraction (case 1), framework-gene hotspot extraction (case 2), and
   orphan-partner subtraction with ESM-2 3Di + Foldseek (case 3) — three of the plan's distinctive claims
   — on datasets of 2,718 genomes (~30 GB) and ~366 GB, which fit the lab server. Acceptance targets:
   1,524 fragments → 94 clusters; 910 hotspot regions with cluster 1 = 259; 205,039 → 124,853 → 2,628.
   Source: <https://pmc.ncbi.nlm.nih.gov/articles/PMC11754643>.

2. **Barrio-Hernandez 2023 AFDB clusters** — full cluster membership and metadata at
   <https://cluster.foldseek.com/> and <https://afdb-cluster.steineggerlab.workers.dev>, with every
   parameter published (E 0.01, 90% overlap both sides, m = 300, k = 10). Acceptance targets: 52,327,413
   MMseqs2 clusters; 18,661,407 structural clusters; 2,302,908 non-singleton; **711,705 dark**. Also
   the only published wall-clock anchor for structural clustering at scale (**129 h on 64 cores**), which
   is directly comparable to the lab server. Source:
   <https://pmc.ncbi.nlm.nih.gov/articles/PMC10584675/>.

3. **Durairaj 2023 AFDB90 SSN** — Zenodo <https://zenodo.org/record/8121336> (SSN + metadata, CC-BY 4.0),
   interactive atlas <https://uniprot3d.org/atlas/AFDB90v4>, code
   <https://github.com/ProteinUniverseAtlas/AFDB90v4>, and — crucially — the resulting **Pfam accessions
   PF21784 / PF21785 / PF21786 (clan CL0395) and PF22187 (clan CL0236)**. A Pfam accession is the ideal
   replay artifact: the seed alignment and HMM are public and versioned. Acceptance targets: 4,270,404
   nodes / 10,339,158 edges; 688,852 communities; 242,876 components. Source:
   <https://pmc.ncbi.nlm.nih.gov/articles/PMC10584680>.

**Tier 2 — hit lists and structures published, pipeline partially reconstructable.**

4. **Duan 2024 structure-guided Acr/anti-defense screen** — ~**300,000 predicted structures** on figshare
   `10.6084/m9.figshare.24948714.v1`, Supplementary Data 1–7, and code at
   <https://github.com/EmiolaLab/ExtractGenes> (Zenodo `10.5281/zenodo.10403005`). The replay with real
   value: re-run the structure search with **Foldseek instead of US-align** on the deposited structures
   and measure recall against the published 920 candidates and 42 anti-defense matches. If Foldseek
   recovers them at ~4,000× lower cost, that result alone changes the plan's structure-search budget.
   Source: <https://pmc.ncbi.nlm.nih.gov/articles/PMC10799925/>.

5. **Rousset 2022 prophage hotspots** — Tables S1 and S5 give **P4 and P2 locus accessions with system
   assignments**, S2 the tested systems, S3 the HHpred Pfam predictions. Replay: re-derive 5,251 P4 loci
   from 20,125 GenBank *E. coli* genomes with the published BLASTP E ≤ 1e-10 anchor, then reproduce the
   541 clusters / 318 arrangements with MMseqs2 `-c 0.6`. This is the cheapest high-value replay in the
   list and tests the plan's locus-extraction layer directly. Source:
   <https://pmc.ncbi.nlm.nih.gov/articles/PMC9122126/>.

6. **Altae-Tran 2023 FLSHclust** — Zenodo `10.5281/zenodo.8371343` (scripts and **all redundant genomic
   loci for all identified systems**), Table S1 (TSV, 38.4 MB), Data S1–S6 including FASTA. The corpus
   itself (8.8 Tbp) is not redistributed and no HMMs or MSAs were released, so the replay is partial:
   take the deposited loci and re-derive the **naive and enhanced scores** (3000 bp array proximity,
   10 kbp repeat window, N_eff ≥ 3, enhanced ≥ 0.35) to check that the plan's enrichment machinery
   reproduces the published ranking. Source: <https://pmc.ncbi.nlm.nih.gov/articles/PMC10910872>.

7. **TED (Lau 2024)** — per-sequence domain assignments on Zenodo (<https://zenodo.org/record/13908086>),
   covering **188 million unique sequences / ~365 million domains**, plus 7,427 candidate novel-fold PDB
   files. The right replay for the plan's sub-ORF tape: compare TED domain boundaries against
   Pfam/InterPro envelopes on the same proteins and quantify how many extra segments structure-based
   segmentation adds. Source: <https://www.biorxiv.org/content/10.1101/2024.03.18.585509.full.pdf>.

**Tier 3 — benchmark harnesses rather than discovery replays (use as acceptance tests).**

8. **Pearson RPD3 + far50/far66** — ~598,000 proteins from Pfam 28 + RefProtDom, 100 far50 and 100 far66
   queries, with published FDR and sensitivity for psiblast, jackhmmer and psisearch2 at up to 10
   iterations. The only published harness that measures *profile drift* directly. Target to beat:
   jackhmmer FDR 13.9% at 10 iterations on far50.
   <https://pmc.ncbi.nlm.nih.gov/articles/PMC5605230/>.

9. **MMseqs2-GPU sensitivity set** — 6,370 queries against 30,430,281 reference sequences (3.4M UniProt +
   27M reversed), with ROC1 published for MMseqs2 (1/2/3 iterations), JackHMMER, BLAST, PSI-BLAST and
   DIAMOND. Reusable as the cross-tool concordance test `00_brief.md` asks for.
   <https://pmc.ncbi.nlm.nih.gov/articles/PMC12510879/>.

10. **Foldseek benchmark data** — <https://wwwuser.gwdg.de/~compbiol/foldseek>, with SCOPe40
    (11,211 structures) and AFDB v1 (365,198 models) and published per-query timings (6 s / 33 h / 10 d
    for Foldseek / TM-align / Dali on a 942-residue query, 64 threads). The direct way to calibrate our
    own Foldseek cost model. <https://pmc.ncbi.nlm.nih.gov/articles/PMC10869269/>.

11. **PLMSearch SCOPe40-test** — 2,207 proteins, 4,870,849 pairs, with published AUROC at family /
    superfamily / fold for MMseqs2 and PLMSearch and remote-pair recall for MMseqs2, BLASTp, Foldseek,
    Foldseek-TM and PLMSearch. The right test for deciding whether an embedding index earns its place
    alongside the plan's existing layers. <https://pmc.ncbi.nlm.nih.gov/articles/PMC10981738/>.

12. **PADLOC-DB and DefenseFinder models** — <https://github.com/padlocbio/padloc-db> (HMMs + covariance
    models, versioned releases, precomputed RefSeq results) and the DefenseFinder model set
    (<https://peercommunityjournal.org/articles/10.24072/pcjournal.470/>). Replay for the Q-SYS layer, and
    the place to audit candidate-system leakage into temporal holdouts (see §e.9).

13. **Serratus and Logan** — <https://serratus.io>, `s3://lovelywater/`, the Logan AWS Open Data registry
    and <https://logan-search.org>. Not replays of a search we would run, but the only published
    cost models for petabase translated and k-mer search, and the numbers to quote when someone proposes
    doing it locally: **$0.0042 per SRA run**, **11 h on 60,000 vCPUs**, **1 PB index**.

**What is missing from the external record, and worth noting as a gap.** Three things I could not find a
public artifact for, despite looking: (i) per-cluster `hhmake`/`cstranslate` build timings for custom HH
databases, which is the dominant unknown in §c.5; (ii) the Methods-level parameters of the
Yirmiya 2025 bait-complex screen (pair counts, score thresholds, GPU hours) — the paper's scale is
published but the pipeline is not retrievable from the sources reached in this round; and (iii) any
modern published throughput for `tblastn` at 10^12-nt scale, which is itself the finding: nobody runs it,
so the plan should not either.
