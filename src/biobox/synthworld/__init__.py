"""synthworld: the synthetic-world generator (planning/tests/synthworld_spec.md; card F-02).

Implemented so far: the backbone (`backbone`, spec §2.1 S1-S2 and §2.2; F-02 sessions 1 and 2), the seed
derivation (`seeds`, spec §6.1), family design and protein-coding sequences (`profiles`, `codons`,
`sequences`, spec §2.1 S3-S4 and §2.3; session 3), and IGRs, regulatory elements, RNA motifs and replicon
sequences (`intergenic`, `rna`, spec §2.1 S5, §2.4 and §2.5; session 4), and phages, prophages and CRISPR-Cas
loci with their spacer targets (`mobile`, spec §2.1 S6; session 5), and gene fusions, splits and
non-orthologous displacement (`fusions`, spec §2.1 S7; session 6), and gene-call errors, genetic codes 4
and 25, assembly classes and contigs (`assembly`, spec §2.1 S8-S9; session 7), and contamination,
duplicates and sources (`sources`, spec §2.1 S10; session 7), and the emitted world: per-source FASTA, the
ingest manifest and build-stage inputs, truth in emitted-contig coordinates and the world manifest (`emit`,
spec §2.1 S11; session 8), and the null worlds N0-full, N0-strict and N1 with the N2 seeds (`nulls`, spec §4.1;
session 9). `backbone` is the only module that imports ZOMBI2 or reads its files.
"""
