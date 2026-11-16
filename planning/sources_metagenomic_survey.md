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
