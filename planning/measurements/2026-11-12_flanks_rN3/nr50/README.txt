Files written by `bio-box nbr flanks` (README.txt; card W-25, planning/design/W-25_flanks.md)

A locus is one anchor gene (a gene of a hit protein) and its window of N genes on each side. Each flanking gene is a
token: its id at the level (NR50 or NR90 cluster uid of the release, or the protein's Pfam architecture), its side
(u: upstream, d: downstream, along the anchor's strand) and its strand relative to the anchor (+ same, - opposite).
A gene with no id at the level is a token of its own (~ and its feature_uid): it is never shared, so unannotated
neighbours never make two loci similar. A side is cut when a linear contig ends inside the window.

groups.tsv          # key<TAB>json facts of the run, then one row per group: the anchor cluster (or all), its loci,
                    its leaves (loci with equal token lists and cuts collapse into one leaf), the leaves in the tree,
                    sampled = more leaves than --max-loci (the leaves with the least names are in the tree), and the
                    number of context types.
context_types.tsv   one row per locus: feature_uid (the anchor gene), group, context_type (empty when the locus's
                    leaf is not in the tree), leaf (the least feature_uid of the leaf's loci), in_tree, the anchor
                    gene (contig_uid, accession, start, end, strand, 0-based half-open), anchor_protein_uid, the
                    anchor's nr90 and nr50, stage, kind, genome_uid, n_genomes and genome_uids (every genome that
                    holds the contig), tax_order (empty: no release has taxonomy yet), cut_up, cut_down (a linear
                    contig end cut that side), n_up, n_down (flanking genes seen on each side).
conservation.tsv    per group, context type and neighbour id (key): n_loci_type (the type's loci), n_loci (loci with
                    the id anywhere in their window), share_of_loci, n_neighbors (genes), n_genomes (distinct genomes
                    of those loci's contigs), offsets (offset:count along the anchor's strand, negative upstream),
                    modal_offset (ties: nearer the anchor, then downstream), same_strand_share. One row per type
                    with unannotated=true gathers the genes with no id.
groups/<group>/     leaves.tsv: one row per leaf (n_loci, in_tree, tree_order, context_type, cuts, tokens outward
                    from the anchor, its loci's feature_uids); distance_j.tsv, distance_a.tsv: J = 1 - Jaccard of the
                    token sets, A = 1 - Jaccard of the adjacent token pairs (anchor, t1), (t1, t2), ... per side, each
                    pair of leaves compared over the positions both observe on each side (4 decimals; 1 when nothing
                    is compared); positions.tsv: those positions per pair; tree.nwk: UPGMA (average linkage, leaves
                    weighted by their loci) on the distance --distance names, node heights = merge distance / 2.
                    Context types are the clades merged at a distance <= --cut, numbered by loci, then leaf name.
gff/                (--export gff) one GFF3 per locus, locus_<accession>_<feature_uid>.gff3: a CDS row per gene of
                    the window, 1-based inclusive, ID=feature_uid, protein_uid, token, offset_rel, anchor; and
                    protein_groups.tsv, protein_uid to its id at the level (its own uid when it has none).

No ordinal is written: genes are named by feature_uid. Ids of clusters are of this release only.
