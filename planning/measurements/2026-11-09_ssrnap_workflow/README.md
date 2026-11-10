# Workflow walk-through: single-subunit RNA polymerases from a custom HMM (2026-11-09)

An agent ran the question "find T7-like single-subunit RNA polymerases with my own HMM and show me
their neighborhoods" end to end on the provisional phage store `rS1.pilot1` (27,258 genomes,
1,134,550 distinct proteins, Pfam 38.2 domains), with the code at commit `2556ddd`. All numbers are
*measured* on wilkcomp02, 2026-11-09, at 8 threads. The gaps it found are cards W-01 to W-08 in
[`planning/tasks/TASKS.md`](../../tasks/TASKS.md). `scripts/` holds every script it used; `export/`
holds the hit table, three example loci and the neighbor summary.

## What works today

HMM in hand → export proteins to FASTA (4.7 s) → `hmmsearch -Z 1e8 --domZ 1e8` (5–8 s) → one DuckDB
join to genes, genomes and clusters → `biobox.search.nbr.fetch` for ±5 genes → TSV and FASTA
(4.3 s for mapping, neighborhoods and export). About 30 s of compute and about 150 lines of ad hoc
code, which needed knowledge of internal tables (uids as hex FASTA names, deduplicating
`contig_occurrence`, list columns in tapes, the `level` filter on `cluster_member`).

## Results

| Item | Value |
|---|---|
| models | Pfam `RNA_pol` PF00940 (stand-in for a user's HMM); a profile built from 8 diverse seeds found in the store (`hmmalign` + `hmmbuild`; the image has no multiple aligner) |
| cutoff | sequence score ≥ 30 bits (E about 1e-3 at Z = 1e8); Pfam's GA of 22 bits is E about 0.3 at this Z |
| hits, union of both models | 2,097 distinct proteins, 2,790 genes, 2,533 genomes, 482 NR90 and 125 NR50 clusters |
| seed-built model vs Pfam at 30 bits | 2,089 vs 1,931 proteins; 1,923 shared |
| full-length class (700–1,300 aa) | 1,571 proteins in 2,072 genomes; median length 822 aa |
| partial hits under 500 aa | 418 proteins; 145 contigs carry adjacent partial pairs (split polymerases) |
| known phages recovered | T7, T3, SP6, K11, phiKMV, gh-1, LKA1, phiYeO3-12, K1F, phiA1122, Berlin, P60, Syn5, P-SSP7: all found |
| not found | N4's 3,500 aa virion polymerase and 252 relatives score under 10 bits with both profiles (a different family; found by its own Pfam domains) |
| false positives at 30 bits | none seen: no DNA polymerase A or multi-subunit RNAP protein in the hit set |
| genome position of full-length genes | two modes: 49% in the first fifth of the genome (T7-like early region), a second mode at 0.40–0.46 (phiKMV-like) |
| commonest neighbors within 5 genes | `Phage_T7_1_1` (48% of loci, +1/+2), DNA ligase (48%, +3/+4), `Head-tail_con` (35%), T7 kinase (32%, −1), `Endonuclease_7` (30%); 60% of distinct neighbor proteins have no Pfam domain |

Observations to follow up, not conclusions: two 1,463 aa `Phage_Mu_F`+`RNA_pol` proteins in
Roseobacter phages that also carry a normal polymerase; 24 full-length proteins of 704–844 aa that
score only 32–70 bits; two identical 888 aa polymerases in one Enterobacter phage genome.

## The intended path (search pack + recipe + runner)

Once worked around, the runner returned exactly the proteins `hmmsearch` did (1,931 for the Pfam
model), with positive and negative controls passing, and a three-round jackhmmer run from T7 gp1
found 2,217 with the drift alarm quiet. It stopped in these places:

1. The search pack and the runs can only live under the data root, which a user cannot write.
2. The pack build crashes on genetic code 15 (`pyhmmer.easel.GeneticCode(15)` does not exist); 61
   phage contigs use it. With code 15 mapped to 11 the pack builds in 495 s (632 MB; dark ORFs are
   421 s of that).
3. A user's own `.hmm` file cannot be a seed: seeds must be models already in the pack's libraries.
4. Under a memory limit the worker failed to start threads (194 threads on the 112-core host; the
   Arrow pool is sized to the host, not to the run's CPU limit), then hung instead of exiting, and
   an outside SIGKILL was recorded as `oom`.
5. `--format` on `bio-box search hits` is documented but not implemented; exports carry uids only.
6. With no domain table present the drift alarm is evaluated against nothing and says nothing.
