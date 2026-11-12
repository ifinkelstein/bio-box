# W-25: tools for comparing the flanks of a set of hits

Review for card W-25, 2026-11-12. The question: which published method should `bio-box` use to say which loci of a
search share a genomic context, and how should it come in. The review starts from Flankophile (the card's source),
follows its citations in both directions in OpenAlex, and adds the neighbourhood tools that do not cite it.

Sources. Every paper is named by its OpenAlex work ID (`W…`). Citation counts are OpenAlex `cited_by_count`,
*measured* 2026-11-12 (api.openalex.org, no key, plain User-Agent). Methods and scale numbers come from each paper's
full text, read through Europe PMC (open-access copies) or the publisher's page. The section is named where the paper
gives one. Licences and code locations are given only when the paper states them. No repository was cloned or read,
except GCsnap's README, which was read because its paper was not reachable.

## 1. What was searched

| Step | Query | Works examined |
|---|---|---|
| Anchor | `/works/doi:10.1128/spectrum.02413-23` | W4389625868 (Flankophile, 17 citations) |
| Backward | its `referenced_works` | 32 of 32. The only flank tool among them is Flanker (W3176511630). The rest are dependencies (Prokka, KMA, CD-HIT, ggtree, treeio, ape, Snakemake, SeqKit, BEDTools) and AMR epidemiology. |
| Forward | `filter=cites:W4389625868` | 17 of 17. All are applications (AMR epidemiology, for example blaTEM in *S.* Typhi, W4403448080, 6,079 genomes) plus one review (W4414107466). None is a new method. |
| Related | its `related_works` | 10. They add ARGprofiler (W4391952941), ARGContextProfiler (W4410571464) and Sarand (W7168296042), all AMR read and assembly-graph context tools. |
| Second hop, forward | citers of FlaGs (165; the 120 matching tool, pipeline, software or server), GCsnap (19), Flanker (23), clinker (132 matching synteny, neighbourhood or genomic context; top 40 by citations) | The second hop found GeCoViz, AnnoView, LoVis4u, NetSyn, GCsnap2, GNAT, zol/fai, the BGC-similarity benchmark, Gaia and DiGAlign. |
| Searches | 25 `search=` queries: FlaGs/webFlaGs, GCsnap, Flanker, clinker, CAGECAT, CORASON, EFI-GNT, pyGenomeViz, LoVis4u, gggenomes, Domainator, MGEfinder, zol, cblaster, BiG-SLiCE, BiG-SCAPE 2, iLund4u, ProFaNA, MultiGeneBlast, and generic terms for neighbourhood clustering by protein family | pyGenomeViz has no paper in OpenAlex. Only software deposits and papers that use it came back, so its row below rests on its name only. MGEfinder concerns insertions in reads, not neighbourhoods, and is left out. |

## 2. Candidates

"Unit" is what a gene or locus is reduced to before two loci are compared. "Cites" is the OpenAlex count.

| Tool (work) | Unit | Distance or score between loci | Grouping, tree | Scale stated | Outputs | Code, licence | Cites |
|---|---|---|---|---|---|---|---|
| **Flankophile** 2023 (W4389625868) | nucleotides of the flanks | KMA `dist` on k-mers (Jaccard, k-mer Hamming, cosine or chi-square; Jaccard with k = 16 in the paper). Three matrices per group: flank, hit, and flank plus hit (Methods, "calculation of distance matrices"). | Hits are grouped by CD-HIT-EST clusters of their closest references (0.98 in the paper). UPGMA (`hclust` "average") in R. | 2,057 isolates and 273 metagenomes. 7,563 ARG hits, of which those with ≥1,500 bp flank on both sides went into 165 clusters (Results). Run on 40 cores and 180 GB. | gggenes/ggtree synteny plots with variant and metadata bars; tables | Bitbucket `genomicepidemiology/flankophile`; Snakemake and conda; licence not stated | 17 |
| **Flanker** 2021 (W3176511630) | nucleotides, one window or windows stepped away from the gene | Mash distance (MinHash of k-mer Jaccard), default threshold 0.001 (Methods, "Flanker") | Connected components of the thresholded graph (single linkage) give "flank patterns". Stepping the window (for example 100 bp steps to 5,000 bp) gives a "Flankergram" that shows where contexts diverge. | 50 KPC plasmids; 425 and 442 short-read collections (Results) | flank FASTA, cluster tables | GitHub `wtmatlock/flanker`, MIT (Introduction) | 24 |
| **FlaGs / webFlaGs** 2020 (W3087148269) | proteins of N flanking genes | none between loci. Neighbour proteins are grouped by jackhmmer (user E-value and iterations), and each group is shown with its conservation. | Rows ordered on a phylogenetic tree of the query proteins (ETE 3), with flanks drawn as pennants | input: RefSeq accessions (about 170 M proteins in RefSeq in 2020, §2) | to-scale vector figure, tree, description table | GitHub `GCA-VH-lab/FlaGs`; webflags.se; licence not stated | 165 |
| **GCsnap** 2021 (W3138507280); **GCsnap2** 2025 (W4411490480) | proteins, 4 flanking genes each side by default (README) | families from all-against-all BLASTp or MMseqs2 (README; thresholds not available). In advanced mode, context types come from PaCMAP clustering of family composition. | Context types; ordering by taxonomy or the target's tree | GCsnap2: "hundreds of thousands of input sequences" on HPC, "22×" faster (abstract) | interactive Bokeh HTML, tables, publication figures | GitHub `JoanaMPereira/GCsnap`; licence not stated | 20; 0 |
| **NetSyn** 2026 (W7214700820) | MMseqs2 families (30% identity over 80% length), window of 11 genes (5 each side) | Synteny score = (GS/2)·(GS/GT), where GS is the number of genes in the shared synteny block and GT all genes of both contexts including gap genes. Gaps of up to 3 genes are allowed; strand is ignored. | Network of anchors with edges at score ≥ 3; MCL, Walktrap, Louvain and Infomap, each run separately | 15,452 candidates gave 13,341 nodes and 310 Louvain clusters. Graph files above 10,000 sequences may exceed Gephi, Cytoscape and the browser. | GraphML, interactive HTML (network, contexts, family table) | GitHub `labgem/netsyn`; licence not stated | 2 |
| **BiG-SCAPE 1 / CORASON** 2019 (W2990471452) | strings of Pfam domains | d = 1 − (w_J·JI + w_DSS·DSS + w_AI·AI). JI is the Jaccard of domain types. AI, the adjacency index, is the Jaccard of adjacent domain pairs. DSS is domain sequence similarity. The weights are calibrated per BGC class against compound similarity; the default is J 0.2, DSS 0.75, AI 0.05 (Online Methods, "Weight optimization"). "Glocal" mode compares the longest common subcluster, extended, and is used when a cluster is fragmented at a contig edge. Strand is not considered. | Network at cutoff 0.3, then affinity propagation (damping 0.9) into families. CORASON builds a FastTree tree of core domains. | 363 strains validated against metabolomics | SVG per cluster, interactive HTML, networks | git.wageningenur.nl (1.x) | 1,067 |
| **BiG-SCAPE 2 / BiG-SLiCE 2** 2026 (W7125601577) | Pfam domains; BiG-SLiCE uses vectors of pHMM scores | BiG-SCAPE 2: as 1, with the alignment unit changed from CDS to domain, a "local" mode, and a preference for an LCS that holds core domains. BiG-SLiCE 2: cosine-like distance on L2-normalised vectors. | Affinity propagation, with fewer families in dense components (density ≥ 0.85). BiG-SLiCE uses BIRCH, which is near-linear. | 292,050 protoclusters: BiG-SCAPE 2 in 89.2 h on 128 cores and about 1 TB, BiG-SLiCE 2 in 6.1 h on 64 cores and about 175 GB (Results, Fig. 3d) | SQLite database, interactive UI | GitHub `medema-group/BiG-SCAPE` and `bigslice`, **AGPL-3.0** (Code availability) | 58 (BiG-SLiCE 1: 235) |
| **Domainator** 2024 (W4405221537) | Pfam or other domain annotations, or orthogroups, per neighbourhood | `compare_contigs`: Jaccard index (unique annotations shared / unique in either) and the adjacency index of BiG-SCAPE (Methods) | Similarity network (edges at Jaccard > 0.7 in Fig. 4); "Jaccard index trees" by UPGMA on 1 − Jaccard (Results, P4 hotspots) | 205,039 neighbourhoods extracted (MT hits ± 4 CDS) from about 366 GB of GenBank (Results, third example) | GenBank, tables, SSN for Cytoscape, trees | GitHub `nebiolabs/domainator`; code licence not stated (paper CC BY-NC). **Already on the planned `bio-box:3.1` tool list** (STAGES.md, SD-13). | 7 |
| **clinker** 2021 (W3121652846) | proteins, global alignment (BLOSUM62, gap −10/−0.5) | S = h + i·s: h is the cumulative identity of linked genes, s the number of shared contiguous pairs, i = 0.5 (adapted from MultiGeneBlast). Gene groups by single linkage of links above an identity threshold. | Ward clustering of the S matrix orders the clusters | not stated; all-against-all alignment, meant for tens of clusters | interactive HTML (clustermap.js, D3), SVG | GitHub `gamcil/clinker`, MIT | 1,669 (CAGECAT front end: 175) |
| **LoVis4u** 2025 (W4407987156) | MMseqs2 protein groups, **or a user table of predefined groups** | proteome composition similarity = fraction of shared homologous groups | average-linkage clustering orders loci; core and accessory genes per cluster | 78 phages and 13,630 proteins in 50 s on an M1 laptop (Results) | vector PDF (ReportLab), similarity matrix | GitHub `art-egorov/lovis4u`; licence not stated | 56 |
| **zol / fai** 2025 (W4407150452) | ortholog groups (InParanoid-like on DIAMOND, MCL) across cluster instances | none between loci. Per ortholog group: conservation (share of instances), consensus order and direction, Tajima's D, entropy, dN/dS. Dereplication with skani. | — | 5,052 *epa* loci: 25 h and 101.6 GB without dereplication, 2.9 h and 2.3 GB with it (30 threads) | XLSX report, `cgc` PDF schematic in consensus order, `cgcg` graph | GitHub `Kalan-Lab/zol`; licence not stated | 17 |
| **GeCoViz** 2022 (W4282934256) | eggNOG, KEGG or Pfam terms | none between loci. Vertical conservation score (VCS) = share of shown genomes that hold the term near the anchor (section "Genomic context panel"). | Rows on a species tree | 12,221 genomes, 42,542,377 genes (Results) | web; SVG and PDF; table | web server | 29 |
| **EFI-GNT** 2023 (W4321207605) | Pfam families within ±20 genes (Results, "EFI-GNT") | none between loci. Co-occurrence of Pfam families with each SSN cluster gives a genome neighbourhood network. | SSN clusters (sequence) | UniProt 2022_04: 230,496,503 sequences (Introduction) | GNN, GND diagrams | GitHub `EnzymeFunctionInitiative`; licence not stated | 435 (2019 paper: 808) |
| **AnnoView** 2024 (W4396940192) | KEGG, Pfam or TIGRFAM labels, ±10 kb | clustering by "similarity of neighboring gene content" (formula in its supplement, not read) | — | 30,238 bacterial and 1,672 archaeal genomes; up to 500 neighbourhoods per view (Results) | web browser | GitHub `satellite830/AnnoView` (workflow) | 21 |
| **ProFaNA** 2023 (W4385066400) | Pfam domains near a query domain | Poisson enrichment per neighbour domain, with optional sub-taxon averaging (resampling to the largest subtaxon) against uneven sampling | — | 45,479 SelO neighbourhoods (Results) | web tables | web server | 5 |
| gggenomes 2024 (W4404649747), pyGenomeViz (no paper) | — | figures only | — | — | R (ggplot2); Python | CRAN; PyPI | 74; — |

Also seen and set aside. iLund4u (W4403499985) annotates hotspots from proteome communities over 873K phage genomes
and 696K plasmids; it suits a later "hotspot" card. Gaia/gLM2 (W4411472741) uses embeddings from a genomic language
model, which is neither deterministic across versions nor buildable here. Pandoomain (W4414765600) runs an HMM search,
takes ±12 genes, and builds a UPGMA tree on the Jaccard of domain architectures. CAGEcleaner (W4411711523) and
SyntenyQC (W4416226586) dereplicate hit sets before plotting. pangraph-based structural diversity (W4390032317),
ARGprofiler, ARGContextProfiler and Sarand all work on nucleotides or assembly graphs.

The one independent benchmark of locus-similarity methods (W4416503985, BGCs against product structure) found that
every method correlates only moderately (Spearman's ρ 0.19 to 0.47, Table 1) and that BiG-SCAPE clusters best after
the data errors were corrected (§3.4). No benchmark of gene-neighbourhood (non-BGC) context distances was found.

## 3. Judgement against what bio-box needs

The needs: (1) work from stored tables (protein uid, Pfam, NR50 and NR90, gene order, strand) without re-annotating;
(2) handle thousands of anchors per search over tens of millions of proteins; (3) deterministic; (4) a licence that
allows adopting the method; (5) publishable figures.

| Tool | (1) our tables | (2) scale | (3) deterministic | (4) method adoptable | Verdict |
|---|---|---|---|---|---|
| Flankophile | no: nucleotides, Prokka, ABRicate | per gene group; KMA is fast | yes | yes (algorithm), licence unstated | Source of the shape (group by anchor, UPGMA, plot with metadata), not of the distance. Its three trees (flank, hit, both) are worth copying. |
| Flanker | no: nucleotides | yes | yes | MIT | Copy one idea: distance as a function of how far from the anchor (the Flankergram), which shows where contexts diverge. |
| FlaGs | partly: jackhmmer groups ≈ our clusters | slow (jackhmmer, NCBI fetch) | mostly | unstated | Display idea only: anchor tree beside the flanks, conservation shown per neighbour group. |
| GCsnap / GCsnap2 | partly (accepts "protein classification maps") | GCsnap2 yes, on MPI | depends on an embedding (PaCMAP); not checked | unstated | Context types and interactive HTML are good ideas; the embedding is not needed. |
| NetSyn | yes, given families | network outputs fail above about 10,000 | community detection varies by method | unstated | Its gap-tolerant score is the only published distance built for neighbourhoods (not BGCs) that tolerates insertions. Worth adopting as the order-aware variant if the adjacency index proves too weak. |
| **BiG-SCAPE** | **yes: it already works on Pfam strings** | 292,050 records in 89 h on 128 cores. Too slow as all-pairs at corpus scale, but W-25 runs per group with at most K loci. | AP is deterministic given input order | AGPL-3.0: reimplement, do not vendor | **Source of the distance**: Jaccard plus the adjacency index on domain or cluster tokens, and the glocal rule for loci cut by a contig end. |
| **Domainator** | **yes: domain or orthogroup content per neighbourhood** | 205,039 neighbourhoods extracted | yes | unstated, but in our planned image | **The closest precedent**: Jaccard and adjacency index on neighbourhoods, then UPGMA. Once in `bio-box:3.1`, its `compare_contigs` is an independent oracle for W-25's distances. |
| clinker | no: aligns sequences | tens of loci | yes | MIT | Use for a user's small subset on exported GenBank; not inside W-25. |
| LoVis4u | **yes, with a predefined-groups table** (our NR50 ids) | 78 genomes in 50 s; figure size is the limit | yes | unstated | Best candidate to run as-is on our exports for a figure of tens to a few hundred loci. |
| zol | no: recomputes ortholog groups | thousands, with dereplication | yes | unstated | Copy its per-group report: conservation and consensus order and direction of each neighbour family. |
| GeCoViz, EFI-GNT, ProFaNA | yes in principle (Pfam terms) | precomputed web resources | — | — | Copy the per-neighbour conservation score (VCS = share of loci). ProFaNA's sub-taxon averaging belongs to the statistics cards (F-01 to F-03), not W-25. |
| AnnoView, webFlaGs, CAGECAT | web services on public genomes | — | — | — | Not usable on our corpus. |

Scale check for W-25 as written. K = 2,000 loci per group gives 1,999,000 pairs, a condensed float64 matrix of about
16 MB and a full square matrix of 32 MB (*derived*: 2,000·1,999/2 pairs × 8 bytes; 2,000² × 8 bytes). UPGMA on that
is quick and bounded. All-pairs comparison does not scale to the corpus, and does not need to: each group holds one
search's hits. The card's figure of "1,957 anchors" for the RNA polymerase run counts proteins. The same run has 2,581
hit loci on 2,526 contigs (GUIDE.md §6, *measured* there), so with `--group-by none` it is above K = 2,000. The test
should say which count it means.

## 4. Recommendation

### (a) What W-25 should implement

1. **Tokens.** Per flanking gene, a token at the chosen level: the NR50 or NR90 cluster id, or at the `pfam` level the
   protein's Pfam architecture string (`pfam_architecture` in `hits.tsv`, domains in order, a repeat written once).
   This follows BiG-SCAPE's Pfam strings (W2990471452) and Domainator's domain content (W4405221537). A gene with no
   annotation at that level gets a token unique to that gene, so two unannotated genes never count as shared. Without
   this rule, "unannotated" would be the commonest shared token on a release without Pfam. Each token also carries
   its strand relative to the anchor (same or opposite) and its side (up or down). That makes the comparison aware of
   operons and of which side a gene is on, which BiG-SCAPE and NetSyn do not do (both ignore strand).
2. **Distance, default: 1 − Jaccard** of the two loci's token sets (both sides together, side-tagged). This is the
   measure of Domainator, LoVis4u (shared share of groups) and Pandoomain, and the protein-level counterpart of
   Flankophile's k-mer Jaccard. "Weighted Jaccard" in the card is ambiguous. Use the multiset (min/max of counts)
   form only if a test shows that repeated families matter, and say so in the note.
3. **Order-aware variant: the adjacency index (AI)** of BiG-SCAPE, the Jaccard of the sets of adjacent token pairs
   along each oriented flank, reported as a second matrix and offered as a combined
   d = 1 − (w·J + (1 − w)·AI) with w stated in the note. Prefer AI over "the share of positions that match", which
   the card suggests. One inserted gene (an IS element, a single cargo gene) shifts every position beyond it, so
   positional matching drops to near 0 for the rest of that flank. AI loses only the pairs that touch the insertion
   (*derived* from the definitions). If AI is still too strict on real data, NetSyn's gap-tolerant block score
   (gaps ≤ 3 genes) is the published fallback. No paper calibrates these weights for gene neighbourhoods. BiG-SCAPE's
   weights were fitted to compound similarity for BGCs (J 0.2, DSS 0.75, AI 0.05 by default), so do not copy them;
   test w against planted truth.
4. **Censoring.** Where a contig end cuts a flank, compare the two loci only over the gene positions both of them
   observe on that side. This is BiG-SCAPE's rule of switching to glocal comparison when a cluster sits at a contig
   edge. Report per pair how many positions were compared, and mark loci with fewer than a stated minimum. Flankophile
   instead drops hits with less than 1,500 bp of flank, which removed about nine in ten of its pig hits, most of them from
   metagenomes (Results, Fig. 2B).
5. **Grouping and tree.** Group by the anchor's own cluster, as in the card (Flankophile's grouping by its closest
   reference). Before any pairwise step, collapse loci with identical token lists into one leaf with a count and the
   list of their `feature_uid`s; strain-redundant genomes make this common (zol, W4407150452: 101.6 GB to 2.3 GB with
   dereplication). Then build a UPGMA tree on the distance matrix (Flankophile; Domainator's Jaccard trees), with
   inputs sorted by `feature_uid` so that ties and output are deterministic.
6. **Flat context types** (missing from the card). Cut the tree at a stated distance and number the clades
   deterministically. This gives one `context_type` per locus in a table. It is what Flanker's flank patterns,
   GCsnap's context types and BiG-SCAPE's families give, and what a user filters on.
7. **Per-type conservation table** (missing from the card). For each context type and each neighbour token: the share
   of loci that hold it (GeCoViz's VCS), its modal offset along the anchor's strand, the share on the same strand,
   and the count of distinct genomes. Together these give a consensus order and direction as zol reports them. Most of
   this already exists for a whole run in `nbr summarize`; this is the same aggregate per context type.
8. **Figure.** As the card says: arrows coloured by family, rows in tree order, metadata columns. Two additions: a
   context-type colour bar, and grey for neighbour families present in fewer than a set share of the group (FlaGs and
   LoVis4u colour only conserved or variable groups, so a few hundred rows stay readable).
9. **Exports for outside viewers** (optional, small). One GFF3 or GenBank file per locus with `protein_uid` and the
   token as qualifiers, plus a TSV of protein to group. LoVis4u (predefined groups), clinker and gggenomes read these
   directly, and Domainator `compare_contigs` can use them as an oracle.

### (b) Running a tool as-is

None is needed for W-25, and none should be built into the release path. Two are worth using on W-25's exports, by
hand, for a figure of a subset:

- **LoVis4u**: it takes predefined protein groups, so our NR50 or Pfam tokens are used as they are. Its output is a
  vector PDF, and it ran 78 genomes in 50 s on a laptop. It needs the per-locus GFF or GenBank and the group table from
  (a) 9, and `pip install lovis4u` (it embeds MMseqs2). This runs off the pod on the small export, or in the image if
  the PI adds it.
- **clinker**: interactive HTML for a few dozen loci chosen from one context type (MIT). It realigns every protein, so
  only for small sets. Same export, `pip install clinker`.

**Domainator** is already planned for `bio-box:3.1` (STAGES.md, SD-13). It is the right external oracle for W-25's
Jaccard and adjacency index (AGENTS §5.2) once the image has it. Until then the oracle is a separate brute-force loop.
Any image change is the PI's decision. No tool here needs to go into `container/environment.yml` for W-25 itself.

### (c) Proposed text for card W-25

Replace the **Do** and **Accept** paragraphs with the following. Keep the **Why** and **Not proposed** paragraphs, and
add to **Why** one sentence: "Method: distance on flank tokens as in Domainator (Jaccard and adjacency index, from
BiG-SCAPE) and UPGMA as in Flankophile; review in `planning/design/W-25_flank_tools_review.md`."

> Do (design note first, half a page, then test first): `bio-box nbr flanks RUN --window Ng --level nr50|nr90|pfam
> [--group-by anchor-nr90|anchor-nr50|none] [--max-loci K] [--cut D] [--export gff]`. Each locus is the anchor's
> flanks as two ordered token lists (upstream, downstream) oriented to the anchor's strand. Each token is (level id,
> side, strand relative to the anchor). At `pfam` the id is the protein's Pfam architecture. A gene without an id at
> the level gets a token unique to it (never shared). Censoring by a linear contig end is stated per locus.
> Distances, both reported: J = 1 − Jaccard of the two loci's token sets (default), and A = 1 − adjacency index (the
> Jaccard of adjacent oriented token pairs, BiG-SCAPE). The combined distance 1 − (w·J_sim + (1 − w)·A_sim) is offered
> with w stated in the note and chosen on the planted test, not copied from BiG-SCAPE. Where a contig end cuts a flank,
> a pair is compared only over the positions both loci observe on that side, and the number of positions compared is
> written per pair. Loci are grouped by the anchor's own cluster. Within a group, loci with identical token lists
> collapse to one leaf with a count before any pairwise step. Per group: the distance matrices (TSV, leaves named by
> `feature_uid`), a UPGMA tree (Newick, inputs sorted by `feature_uid`), context types from cutting the tree at D (one
> row per locus: `feature_uid`, group, `context_type`), and a conservation table per context type and neighbour token
> (share of loci, modal offset along the anchor's strand, share on the same strand, distinct genomes through
> `genome_member`). Also one figure (SVG and PDF, matplotlib; no R): genes as arrows coloured by family, families
> below a set share in grey, rows in tree order, columns for context type, genome kind, stage, taxonomy and the
> anchor's NR90. Optional `--export gff`: one GFF3 per locus with `protein_uid` and the token as attributes, plus a
> protein-to-group TSV, readable by LoVis4u, clinker and Domainator. Bounded: at most K leaves per group after
> collapsing (default 2,000; above it, a deterministic sample by `feature_uid`, said in the output). Deterministic
> output (AGENTS §4.3). Every row names `feature_uid` beside any ordinal (AGENTS §3.3).
>
> Accept: tests on a constructed release with planted contexts. (i) Two flank types around one anchor family, one
> shared by phages and prokaryotes, are recovered as two context types. (ii) A pair with the same families in a
> rearranged order is separated by A and not by J. (iii) One inserted gene changes A by at most the pairs it touches.
> (iv) A censored locus is compared only over the positions it observes. (v) Unannotated neighbours do not make two
> loci similar. The oracle for J and A is a brute-force loop written separately, and Domainator `compare_contigs` on
> the export once the image has it. On `rN.3`, the RNA polymerase run (state whether the count is 1,957 proteins or
> 2,581 loci) is grouped with time and peak memory recorded (Pod: yes, read-only).

## 5. Limits of this review

- Code licences are given only where the paper states them. FlaGs, GCsnap, NetSyn, LoVis4u, zol and Domainator do not
  state one in the text read. Check the repository before vendoring anything; the recommendation does not vendor.
- Last release or commit dates are not visible from the papers or OpenAlex. The newest archived versions the papers
  name are BiG-SCAPE 2.0.0-beta.8 (Zenodo 17865391) and NetSyn (Zenodo 22878461). GCsnap2 (2025) shows that GCsnap is
  being developed further.
- AnnoView's clustering formula and GCsnap's family thresholds are in supplements or documentation that were not
  read.
- Full texts were read through Europe PMC and publisher pages, which reach beyond OpenAlex; no other service was used.
