# W-25: grouping a search's loci by their flanks (design note)

Card W-25; method from `W-25_flank_tools_review.md` §4 (Domainator's Jaccard and BiG-SCAPE's adjacency index on
tokens, Flankophile's grouping and UPGMA), reimplemented from the papers. No code is copied.

**Command.** `bio-box nbr flanks RUN|UID,... --root R [--release R] --out DIR [--window 5g] [--level nr50|nr90|pfam]
[--group-by anchor-nr90|anchor-nr50|none] [--max-loci 2000] [--cut 0.5] [--distance j|a|combined] [--export gff]`.
Loci are `nbr.fetch`'s with edge `truncate`, circular `auto`; only gene windows (`Ng`) are accepted, since positions
are compared. DIR is new or empty, outside the root, filled under a temporary name and renamed into place.

**Tokens.** Per flanking gene of a locus: (id, side, strand relative to the anchor). Side `u` (offset_rel < 0) or `d`
(> 0); strand `+` when the gene's strand equals the anchor's, else `-`. The id is the gene's protein's NR50 or NR90
cluster uid of the release, or at `pfam` its architecture as `hits.tsv` writes it (domain names in order, a run of one
name once, `+`-joined; library as `--library` or the only one). A gene with no id at the level (no cluster row, no
domain, a non-CDS slot) gets a token unique to that gene in that locus (`~` + anchor feature_uid + feature_uid), so it
is never shared, not even by two loci whose windows hold the same gene. Each side is a list ordered outward from the
anchor (|offset_rel| ascending); this orients every locus to its anchor's strand.

**Censoring.** A side is cut when the contig is linear (tape `is_circular` false) and the side holds fewer than N
genes; `obs(side)` is its number of genes when cut, else N (a small circular contig observes all it has). A pair (X, Y)
is compared over `m = min(obs_X, obs_Y)` positions on each side: each side's list is truncated to its first `m`
entries before both distances. `positions(X, Y) = m_u + m_d`, written per pair. Each locus row says `cut_up`,
`cut_down`, `n_up`, `n_down`.

**Distances.** On the truncated lists. J = 1 − |S_X ∩ S_Y| / |S_X ∪ S_Y|, S the set of tokens of both sides (sets: a
token repeated on a side counts once; the multiset form is not used, review §4a.2). A = the same Jaccard on the sets of
adjacent pairs: on each side, with `^` for the anchor, the ordered pairs (^, t1), (t1, t2), …, (t_{m−1}, t_m). Pairs are
ordered (inner, outer), possible because every list is oriented to its anchor; a side of m genes has m pairs. An empty
union (m_u = m_d = 0) gives J = A = 1 with 0 positions: nothing compared is never similar. Combined:
1 − (w·(1 − J) + (1 − w)·(1 − A)). w is chosen on the planted test (below): the largest w on the grid 0, 0.05, …, 1 at
which every planted between-type distance exceeds every within-type distance, i.e. as much weight on the default J as
still separates a rearranged context. The test computes it from the oracle and asserts the constant equals it: w =
0.2 (*derived* on the planted world of `tests/r2/test_w25_flanks.py`: the tightest pairs are T1 against its rearranged
copy, 0.667·(1 − w), and the insertion copy against the unannotated copy, 0.571 − 0.238·w, so the margin is
0.095 − 0.429·w, positive below w = 0.222). On that world J alone (w = 1) cannot separate the rearranged copy (J = 0).

**Grouping and collapsing.** Groups by the anchor protein's own cluster (`anchor-nr90`, default; `anchor-nr50`), or
one group `all`. In a group, loci with equal (up list, down list, cut_up, cut_down) collapse to one leaf with a count;
its name is the least feature_uid of its loci. Leaves are sorted by that name; above K leaves (default 2,000) the K
with the least names are kept (a deterministic sample by feature_uid that keeps a leaf with probability rising with
its loci), the rest are written with `in_tree=false` and no context type, and the output says so.

**UPGMA.** Average linkage on the chosen distance, each leaf weighted by its count, so the tree equals UPGMA on the
uncollapsed loci. Each step merges the pair of active clusters with the least distance; ties go to the pair whose
lesser members' ranks (in the sorted leaf order) are least, lexicographically. A merged cluster keeps the lower rank.
Node height = merge distance / 2; branch length = parent height − child height; the left child is the lower rank.

**Context types.** Cut at merge distance D (default 0.5, `--cut`): clusters formed by merges ≤ D are types, numbered
1, 2, … per group by loci descending, then least leaf name.

**Conservation** (per group, context type and level id, over the type's loci in the tree, full windows): `n_loci`,
`share_of_loci`, `n_neighbors`, `offsets` and `modal_offset` (offset_rel; ties to the offset nearer the anchor, then
downstream, as `nbr summarize`), `same_strand_share`, `n_genomes` (distinct genomes of the loci's contigs: `contig_dim`
genome and every `genome_input` genome of the release, i.e. `genome_member`). All unique tokens form one
`unannotated=true` row.

**Outputs** (DIR): `README.txt`; `groups.tsv` (`# key<TAB>json` facts, then group, n_loci, n_leaves,
n_leaves_in_tree, sampled, n_context_types); `context_types.tsv` (one row per locus: feature_uid, group,
context_type, leaf, in_tree, contig_uid, accession, start, end, strand, anchor_protein_uid, anchor_nr90, stage, kind,
genome_uid, n_genomes, tax_order, cut_up, cut_down, n_up, n_down); `conservation.tsv`; per group
`groups/<group>/leaves.tsv` (leaf, n_loci, tree_order, in_tree, context_type, cut and observed counts, tokens,
feature_uids), `distance_j.tsv`, `distance_a.tsv`, `positions.tsv` (square, leaves by name, 4 decimals),
`tree.nwk`. `--export gff`: `gff/locus_<accession>_<feature_uid>.gff3` (one CDS row per gene of the window, 1-based,
`ID=feature_uid`, `protein_uid`, `token`, `offset_rel`) and `gff/protein_groups.tsv` (protein_uid, group). No ordinal
appears anywhere; genes are named by feature_uid.

**Bounds and determinism.** One group's matrices at a time (three K×K arrays and the UPGMA working copy: 112 MB at
K = 2,000, *derived*: 2·2,000²·8 + 2,000²·4 + 2,000²·8 bytes; 124 MB *measured* on rN.3,
`planning/measurements/2026-11-12_flanks_rN3.md`). Pure Python and NumPy over explicitly sorted inputs; DuckDB
reads are ordered; no clock, path or hash order in any output.

**The figure** was added by card W-27 (2026-11-12, `--figure`; matplotlib approved by the PI).

**Not done.** No nucleotide distance, no bp windows, no
significance, no NetSyn gap-tolerant score, no Domainator oracle (not in the image), no GenBank export.
