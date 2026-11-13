"""synthworld: the synthetic-world generator (planning/tests/synthworld_spec.md; card F-02).

Implemented so far: the backbone (`backbone`, spec §2.1 S1-S2 and §2.2; F-02 sessions 1 and 2), the seed
derivation (`seeds`, spec §6.1), and family design and protein-coding sequences (`profiles`, `codons`,
`sequences`, spec §2.1 S3-S4 and §2.3; session 3). `backbone` is the only module that imports ZOMBI2 or reads
its files.
"""
