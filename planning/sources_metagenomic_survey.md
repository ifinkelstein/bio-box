# Metagenomic sources beyond the current ladder: a survey

Date: 2026-11-16. Asked for by the PI: the best, most comprehensive collection of metagenomic data without too much
redundancy; incomplete contigs and short reads too; viral and bacterial, from diverse niches. This is a survey, not a
decision: nothing here is fetched, and no source is added to `sources.yaml`.

Labels. Every count below is *published* by the resource named (paper or site, with its year), read through web
searches on 2026-11-16; none is measured by us. Where a figure is my own arithmetic it is marked *estimate*. Where
two sources disagree, both are given.

## 1. What the discovery literature searched

From the methods of the papers in the sibling project's `bioinformatics-synteny/lit/` (236 with text; this project
has no `lit/` folder). The recent discovery papers, the kind of result bio-box is built to reach, used:

| Paper | What it searched |
|---|---|
| Millman 2020 and 2022, Bernheim 2021 (Sorek lab) | 38,167 bacterial and archaeal genomes from IMG (October 2017) |
| Vassallo 2022 | 844,603 genomes from GenBank and RefSeq |
| Zhang 2025 (Kiwa); Wu 2026 | PADLOC's database of 223,143 prokaryotic genomes; RefSeq |
| Benler 2021 (Tn7) | NCBI whole-genome shotgun (WGS) nucleotide database |
| Özcan 2021 (Cas7-11) | bacterial and metagenomic data from NCBI, JGI and ENA |
| Faure 2023 (Tn7 target selectors) | contigs of NCBI GenBank and WGS plus JGI: **521,828,662 contigs over 1.5 kb, 6.93 Tbp** |
| Altae-Tran 2021, Faure 2025, Tordoff 2025 (Zhang and Koonin labs) | "genomic and metagenomic databases" (the same NCBI plus JGI contig collection) |
| Breaker lab 2022 to 2023 | RefSeq (44,270 genomes) plus environmental sequences |
| Makarova 2024 | UniRef30 and MGnify proteins, for profile searches |

So the working collection of the field's largest searches is unbinned contigs of NCBI WGS and JGI, about 7 Tbp. The
data root holds 0.06 Tbp (*measured*, 2026-11-16).

## 2. Sources, in increasing size and diversity

| # | Source | What it is | Size (published) | Niches | Overlap with what we have or plan | Access |
|---|---|---|---|---|---|---|
| 1 | UHGV (2025 preprint) | human gut viruses, built from the earlier gut catalogues | 873,994 genomes; 168,536 vOTUs | human gut | replaces MGV, GPD and the other gut catalogues it merges | JGI site |
| 2 | Global Soil Virus Atlas (2024) | uncultivated viral genomes from 2,953 soil metagenomes | 616,935 genomes; 38,508 vOTUs | soil | part may be in IMG/VR | paper's deposit |
| 3 | IMG/VR v4 (2023) | viral genomes and fragments, DNA and RNA, proviruses included | over 15 million; 8.7 million vOTUs; 231,408 vOTUs with a high-quality member | every IMG biome | our stage S4 takes only the high-quality vOTU representatives | JGI policy; local only |
| 4 | IMG/PR (2023) | plasmids from genomes, metagenomes and metatranscriptomes | 699,973 | every IMG biome | PLSDB (stage `PL1`) holds isolate plasmids | JGI policy |
| 5 | VIRE (EMBL, 2025) | viral genomes from over 100,000 public metagenomes | over 1.7 million genomes | SPIRE's habitats | unknown against IMG/VR; not checked | EMBL |
| 6 | RNA viruses: LucaProt (2024), Serratus (2022), RVMT (2022) | RNA-virus contigs and RdRps from metatranscriptomes and all of SRA | LucaProt: 161,979 species from 10,487 metatranscriptomes; Serratus: about 132,000 novel RdRps; RVMT not checked | global | none: we hold no RNA viruses | papers' deposits |
| 7 | **GlobDB r226** (2025) | species-dereplicated prokaryote genomes from 14 resources, GTDB first | **306,260 species representatives**, against GTDB R226's 143,614 | every biome the MAG catalogues cover | by construction none within itself (96% ANI); the 143,607 GTDB ones are our S3 and S5 | open; one download per release |
| 8 | Biome MAG catalogues: GEM (2021), OMD (2022), mOTUs database (2025), UHGG, soil, glacier and others | all genomes behind the species of #7, strains included | GEM 52,515 MAGs; OMD 34,799 genomes; mOTUs 2.83 million MAGs and 919,090 isolate and single-cell genomes in 124,295 species | per catalogue | heavy with GTDB and with each other; #7 is their dereplicated union | open |
| 9 | AllTheBacteria (2024, increment 2025-05) | every bacterial isolate assembly, uniformly assembled | 2,440,377 at the 2024-08 increment | isolates; a few pathogen species dominate | high: strain-level copies of species we hold | open |
| 10 | **SPIRE** (2024) | 99,146 metagenomes of 739 studies, assembled by one pipeline, with curated habitat labels | **16 Tbp of contigs; about 35 billion proteins; 1.16 million MAGs** | 92,134 species-level clusters new at the time | its MAG species are in #7; its unbinned contigs are in nothing we plan | EMBL; per sample |
| 11 | MGnify assemblies and proteins | EBI's assemblies of public metagenomes | proteins: 1.66 billion in release 2026_07 (the 2023 paper says over 2.4 billion; not reconciled) | 28,000 or more sources | the same SRA samples as #10 and #12, assembled by another pipeline | EBI |
| 12 | NCBI WGS and IMG/M contigs | what §1's papers searched | 6.93 Tbp over 1.5 kb (Faure 2023); IMG holds 75.1 billion metagenome genes (2022) | every biome | the same samples as #10 and #11 in part | NCBI open; JGI partly restricted (legacy use terms) |
| 13 | **Logan v1** (2024 preprint) | every public SRA run of December 2023 assembled: unitigs and contigs | 27 million runs, 50 Pbp of reads; **0.90 Pbp of contigs** (385 TB by one talk, 0.3 PB by another); 4.59 Pbp of unitigs (2.18 PB) | everything sequenced | contains the reads behind #10 to #12; one assembly per run, so no duplicate assemblies | AWS open data; compute next to it |
| 14 | Raw SRA short reads | the reads themselves | 50 Pbp (December 2023); MetaGraph indexes about 67 Pbp | everything | fully behind #13 | NCBI, AWS |

## 2b. Plasmids and other mobile elements (added at the PI's word, "and plasmids too")

| Source | What it is | Size (published) | Overlap | Access |
|---|---|---|---|---|
| PLSDB (2025 update) | complete plasmids of isolates, curated from NCBI | 72,360 entries; our stage `PL1` is its 2024_05_31_v2 release, 72,556 records (*measured*, card C-26) | this is what we hold | open |
| Human gut plasmids (Yu, Fogarty and Eren, 2024) | plasmids predicted by PlasX in gut metagenomes | 68,350 non-redundant plasmids in 1,169 plasmid systems | not in PLSDB (metagenomic); not checked against IMG/PR | paper's deposit |
| IMG/PR (2023) | plasmids found in genomes, metagenomes and metatranscriptomes, with ecosystem and host labels | 699,973 | its isolate part overlaps PLSDB; the metagenomic part is new | JGI policy |
| PlasmidScope (2024) | a union of ten plasmid repositories, annotated | 852,600 plasmids; 25.2 million genes | contains PLSDB and IMG/PR among its ten; duplicates across repositories are not removed by sequence as far as its abstract says | open site |
| ICEberg 3.0 (2023) | integrative and conjugative elements, curated | 2,065 ICEs, 607 IMEs, 275 CIMEs | these lie inside genomes we hold or will hold: an annotation, not new sequence | open |

Plasmids that are fragments (most of those in metagenomes) arrive with the contigs of §2 (#10 to #13) and are
found there by prediction, which the plan already has (geNomad), so the catalogues above matter most for complete
or near-complete plasmids and for their habitat and host labels. Not checked today: mMGE, the phage-plasmid sets,
and plasmidome studies of wastewater, soil and the ocean.

## 3. What I would do (a recommendation, not a decision)

The redundancy rule that follows from the table: **genomes by species representative; contigs by one assembly per
sequencing run; viruses by vOTU representative; and our own keys do the rest** (`contig_uid` is the hash of the
sequence, so the same contig from two sources is stored once; proteins collapse at NR90).

1. **Viruses and plasmids from diverse niches first** (#1, #2, #3's high-quality set as vOTU representatives; for
   plasmids IMG/PR, then the gut plasmids of §2b, then whatever PlasmidScope holds that neither has, matched by
   sequence hash): small, and the widest diversity of mobile elements per base. This widens stage S4 from IMG/VR
   alone and makes IMG/PR a firm part of it instead of an option. A few tens of Gbp (*estimate*).
2. **GlobDB's non-GTDB species** (#7): 162,653 representatives that GTDB does not have (*derived*: 306,260 less
   143,607), the uncultivated majority, with no species-level duplicate. Taken with S3 and S5 it makes the stage
   "every known prokaryote species". About 400 Gbp (*estimate*: 162,653 at 2.5 Mbp a MAG).
3. **Unbinned contigs, one source: SPIRE** (#10), sampled by habitat first. A pilot of about 1,000 samples spread
   over its habitat labels (about 160 Gbp, *estimate*: 16 Tbp over 99,146 samples) measures what the plan already
   asks of every later stage: new NR50 clusters per Tbp, per habitat. That number decides how much more to take.
   MGnify and WGS (#11, #12) are then taken only for samples SPIRE lacks, matched by run accession.
4. **Short reads: never as reads.** Logan (#13) is the non-redundant form of all of SRA. Use its search index (or
   MetaGraph) to find the runs that hold a family of interest, and pull the contigs of metagenome runs that no
   assembly above covers. This is the step that needs compute next to the data, or TACC.
5. **Left out unless asked**: AllTheBacteria and the full strain sets of #8 (redundant for diversity); RNA viruses
   (#6) until gene calling and the token model for RNA genomes are decided, since everything built so far assumes
   DNA.

Open points this survey did not settle: the overlap of VIRE with IMG/VR; SPIRE's download form and size per sample;
which GlobDB members are also in SPIRE's own species set (GlobDB took 49,369 of SPIRE's 107,078); the terms of
each JGI dataset; and current totals for MGnify assemblies and IMG/M contigs.

## 4. How to stage the imports when annotating everything is out of reach

Asked by the PI on 2026-11-16 ("we probably can't annotate them in any reasonable way but think thru and give me
ideas"). What follows is reasoning from our own measured rates; the conclusions are suggestions, and card C-38
turns them into a design to decide.

### 4.1 What each depth costs

Rates are *measured* on wilkcomp02 (their cards named); everything scaled from them is *derived*, and assumes a
metagenomic contig costs what our genomes cost per base, which the pilot of card C-42 measures.

| Depth | What it gives | Cost per Gbp | 16 Tbp (all of SPIRE's contigs) |
|---|---|---|---|
| 0. Stored | sequences by content hash, with source and habitat labels | ingest about 4.5 CPU-minutes (42.9 per 9.6 Gbp, the code's constant) | 1,200 CPU-hours |
| 1. Called | genes, proteins, gene order (tapes), CRISPR arrays | gene calling in meta mode 60 to 91 CPU-minutes (2.95 to 4.42 s per Mbp × 1.23, card B-48); arrays 3 (182 s per Gbp, card C-24) | 16,000 to 24,000 CPU-hours |
| 2. Clustered | protein families by sequence alone | not measured at this scale; our cascade does not scale (14,087 CPU-minutes for 45.8 M proteins, exponent 1.244, card B-36); a linear clustering (linclust) is the candidate | unknown; needs memory in the TB range |
| 3. Annotated | Pfam domains on every distinct protein | about 1,550 CPU-minutes (H1: 2,817 for 1.82 Gbp, card C-24) | 413,000 CPU-hours |

- Annotation costs about 20 times what calling genes costs. On the pod at 40 threads 413,000 CPU-hours is 1.2
  years; at TACC it is about 3,200 node-hours if a Lonestar6 core does what a pod thread does (card C-37 measures
  that). For Logan's 0.9 Pbp it is about 180,000 node-hours. So: depth 3 for everything is out on the pod, a
  budget question at TACC for a SPIRE-sized set, and out anywhere for Logan.
- **Storage binds before CPU.** Stage H1's files at depth 1 are 0.71 GB per Gbp, 0.37 without the nucleotide
  sequences and IGRs (*derived* from the sizes in card C-24's file list). So 1 Tbp is about 710 GB, or 370 GB
  without sequences; 16 Tbp is 11 TB, or 6 TB. The pod's local disk is 1.1 TB in all and the data root lives on
  it (80 GB today); `/stor/scratch` is a shared 22 TB (*measured*, `df`, 2026-11-16).
- One search costs little: an HMM over `rN.1`'s 10.6 M proteins took about 8 to 11 CPU-minutes (*derived* from
  card W-21's 190 to 245 s at 2.6 CPUs). At that rate one model over SPIRE's 35 billion proteins is about 450 to
  580 CPU-hours, four node-hours at TACC: not interactive, but a routine job.

### 4.2 Ideas, in the order I would try them

1. **Separate importing from annotating, as named depths.** A stage records its depth; a release and a search say
   which depths they cover. Genome-like sets (viruses, plasmids, species representatives) go to depth 3 as now.
   Unbinned contigs stop at depth 1.
2. **Search first, annotate on the hit.** The discovery papers of §1 work this way: find the bait in everything,
   pull each hit's neighbourhood, annotate only the neighbours. For us: `search hmm` runs over the depth-1
   proteins; Pfam then runs on the proteins within N genes of a hit, and those rows are kept, so the annotated
   part of the corpus grows where questions were asked. Cost follows the number of hits, not the size of the
   corpus.
3. **Cluster before annotating, and annotate representatives of families that matter.** With a linear clustering
   the family of a protein is known without Pfam, so neighbour statistics can run on cluster tokens alone. Pfam
   then goes only to representatives of clusters above a size, or near hits. This is card B-47's idea taken one
   step further.
4. **Filter at the door.** A contig with fewer than about five genes carries no context. Dropping contigs under
   about 5 kb (the papers of §1 used 1.5 kb) with the dropped bases counted in the ledger could remove a large
   share of an assembly; the pilot measures how much.
5. **Keep proteins and coordinates; refer to the nucleotides.** For unbinned contigs, call genes and CRISPR arrays
   as the sequence streams through, keep the proteins, the gene order and the arrays, and record where the
   sequence can be fetched again (source, sample, contig, sha256) instead of storing it. Halves the storage. It
   changes an invariant of the spec (the corpus holds its sequences), so it is the PI's decision.
6. **Take a habitat until it stops paying.** Import samples habitat by habitat and track new NR50 families per
   Gbp; stop a habitat when the rate falls under a set fraction of its first value. The stage plan already asks
   this of every stage after S6.
7. **Import on the hit, from Logan.** Do not import Logan. Search it where it lives (its own index for close
   matches; a protein search over its contigs is a cloud-sized job to be costed), then pull only the contigs of
   the runs that hit, as a small stage that records the query that brought it in.
8. **Use the right machine for each depth.** Gene calling and scans are many independent pieces: TACC's ordinary
   nodes. A clustering of billions of proteins needs one very large memory: Stampede3's 4 TB nodes.
9. **Decide where a multi-TB corpus lives before importing one.** The pod's local disk cannot hold it. The
   candidates are a larger footprint on `/stor/scratch` (ours to ask for, `planning/footprint_request.md`) or
   TACC storage; either way the data root, or the depth-1 part of it, moves off `/tmp`.

### 4.3 The order of the cards

Sources that can still be annotated in full come first, because they are small and widen diversity most per base:
viruses (C-39), plasmids (C-40), then GlobDB's species (C-41; about 400 Gbp, which is 57 node-hours of scan at
TACC by the rates above, or about eight days of the pod at 40 threads). Then the unbinned contigs, as a measured
pilot before any commitment (C-42, with the design C-38 and the run-accession crosswalk C-43). Logan is probed,
not imported (C-44). What stays out, and why, is card C-45.
