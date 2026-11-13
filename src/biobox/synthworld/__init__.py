"""synthworld: the synthetic-world generator (planning/tests/synthworld_spec.md; card F-02).

Implemented so far: the backbone (`backbone`, spec §2.1 S1-S2 and §2.2; F-02 sessions 1 and 2), the seed
derivation (`seeds`, spec §6.1), family design and protein-coding sequences (`profiles`, `codons`,
`sequences`, spec §2.1 S3-S4 and §2.3; session 3), and IGRs, regulatory elements, RNA motifs and replicon
sequences (`intergenic`, `rna`, spec §2.1 S5, §2.4 and §2.5; session 4), and phages, prophages and CRISPR-Cas
loci with their spacer targets (`mobile`, spec §2.1 S6; session 5). `backbone` is the only module that
imports ZOMBI2 or reads its files.
"""
