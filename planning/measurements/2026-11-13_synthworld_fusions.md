# synthworld S7 (gene fusions, splits, NOGD) on `micro`: time, memory and realized values

Card F-02, session 6 (branch `task/F-02-6`). Code: `src/biobox/synthworld/fusions.py` (its docstring states the
model), the `fusions` section of `presets/micro.yaml`, and `mobile.PHAGE_SECTIONS` (the phage world reads only
the S1-S5 sections, so a later section leaves S6 unchanged).

Host and command for every *measured* number below: wilkcomp02, 2026-11-13, one process at nice 19 inside
`/tmp/if2453/images/bio-box-test_64f1bf0.sif`, with the account at about 34 to 37 threads of its 40-thread day
cap (a 32-thread Pfam scan and other agents' tests);
`env BIOBOX_POD_IMAGE=/tmp/if2453/images/bio-box-test_64f1bf0.sif BIOBOX_POD_PYDEPS=/tmp/if2453/ci/F-02-6/src
sh scripts/pod_pytest.sh F-02-6 -s tests/r2/test_synthworld_fusions.py` (the second run with every synthworld test and `tests/r0` in the same call), `micro` default seed 20261112. In that second run S1-S6 measured 36.55 s (backbone), 78.71 s (sequences), 45.91 s (replicons) and 10.57 s (mobile) in their own fixtures.

## Time and memory

| quantity | value | label |
|---|---|---|
| S7 wall time, `fusions.build` in the test fixture | 4.96 s, 5.09 s (two runs) | measured |
| S7 wall time in a child process (`PYTHONHASHSEED=4242`, one thread) | 6.06 s, 6.33 s | measured |
| S7 peak RSS in that child (VmHWM) | 210 MB, 213 MB | measured |
| time left for S7-S11 before this session | about 12 s | derived (session 5: 180 s target - about 157 s for S1-S5 - about 11 s for S6) |
| time left for S8-S11 after S7 | about 6 to 7 s | derived (12 s - 5 to 6 s) |

No document gives an estimate for S7 alone, so there is no 20% comparison. Most of S7 is reading and writing the
38 edited host genomes again (about 80 Mb, *estimate* from 38 of 46 genomes of about 2.2 Mb each); the 8 genomes
without an edit are copied byte for byte. S9 (assembly) and S11 (emit) will each rewrite every genome, and at
S6's and S7's cost per rewrite (5 to 11 s, measured above and in the S6 note) the 3-minute target cannot hold
for S8 to S11 as written: about 3 min 15 s to 3 min 30 s is the *estimate* for S1-S11 if S8, S9 and S11 each
cost about what S7 does and S10 about 1 s. Either the target moves to about 3.5 min, or S8 to S11 share one
rewrite (S8 edits and S9 cuts planned per genome and applied in one pass).

## Realized values (measured)

| quantity | value |
|---|---|
| events | 15 fusions, 5 splits, 2 NOGD (§5.1: 15, 5, 2) |
| instances (event x clade genome) | 81: 52 fused ORFs, 14 split pairs, 15 NOGD genes |
| clade sizes | 2 (8 events), 3 (6), 4 (1), 5 (4), 6 (2), 11 (1) |
| genomes with an edit | 38 of 46 |
| genes after S7 | 90,242: host 88,969, prophage 1,058, Cas 120, fused 52, split halves 28, NOGD 15 |
| S6 genes replaced by event genes | 133 (2 per fused ORF, 1 per split pair, 1 per NOGD gene) |
| S6 genes, insertions, prophages, arrays, spacer slots that read the same after the edits | 90,147 genes, 46, 16, 30, 3,462 |
| S5 features (through S6's insertions and S7's edits) that read the same | 171,797 elements, 1,078 motifs, 82,371 IGRs |
| S5 features an edit touched (all inside a fused pair's junction) | 75 elements, 0 motifs, 34 IGRs |
| every gene's bases translate to its protein_uid | 90,242 of 90,242 |
| S7 digest | `db09abf3fe79ae34a73e3ab0664a9b13`, byte-identical in the child |
| earlier stage digests, before and after S7 | backbone `025fdd0d27fdc65e678782216ecd21fb`, sequences `34ebce1324413a261acd874d7ee4fe3e`, replicons `55bda97155d75148e5cdca7b90a1d7be`, mobile `0bd7f0db5691de58ea5e1a0df1e6e570` (session 5's measured S6 digest) |
