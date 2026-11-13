# synthworld S8-S10 (call errors, codes 4 and 25, assembly, contamination, duplicates) on `micro`: time, memory and realized values

Card F-02, session 7 (branch `task/F-02-7`). Code: `src/biobox/synthworld/assembly.py` (S8 and S9 in one rewrite
of the genomes) and `src/biobox/synthworld/sources.py` (S10); their docstrings state the model. Preset sections
`callerr`, `assembly`, `contamination`, `duplicates` of `presets/micro.yaml`. Tests:
`tests/r1/test_synthworld_assembly_units.py`, `tests/r1/test_synthworld_sources_units.py`,
`tests/r2/test_synthworld_assembly.py` (one fixture builds S1-S10, so CI builds the world once).

Host and command for every *measured* number below: wilkcomp02, 2026-11-13, one process at nice 19 inside
`/tmp/if2453/images/bio-box-test_64f1bf0.sif`, with the account at about 27 to 33 threads of its 40-thread day
cap (a 32-thread Pfam scan and other agents' tests);
`env BIOBOX_POD_IMAGE=/tmp/if2453/images/bio-box-test_64f1bf0.sif BIOBOX_POD_PYDEPS=/tmp/if2453/ci/F-02-7/src
sh scripts/pod_pytest.sh F-02-7 -s <test files>`, `micro` default seed 20261112. Runs: **A** (13:25, S1-S9 tests
only), **B** (13:41, S1-S10), **C** (13:47, every synthworld test and `tests/r0` in one call, after the S10 fix
below; 453 passed, 6 skipped, 2 failed, both test-side and fixed as described at the end), **D** (14:13, the
changed test files again).

## Time and memory

| quantity | A | B | C | label |
|---|---|---|---|---|
| S1-S2 backbone (fixture) | 34.95 s | 36.48 s | 36.78 s | measured |
| S3-S4 sequences | 76.01 s | 78.35 s | 74.25 s | measured |
| S5 replicons | 42.79 s | 43.13 s | 42.16 s | measured |
| S6 mobile | 10.57 s | 10.66 s | 10.62 s | measured |
| S7 fusions | 5.85 s | 4.83 s | 5.01 s | measured |
| S8-S9 assembly (fixture) | 17.61 s | 17.71 s | 18.83 s | measured |
| S8-S9 in a child (`PYTHONHASHSEED=4242`, one thread) | 17.72 s | 19.20 s | 18.53 s | measured |
| S8-S9 peak RSS in that child (VmHWM) | 178 MB | 184 MB | 183 MB | measured |
| S10 sources (fixture) | - | 6.70 s | 7.22 s | measured |
| S10 in a child (one thread) | - | 6.85 s | 6.96 s | measured |
| S10 peak RSS in that child | - | 129 MB | 128 MB | measured |
| S1-S9 in one run | 187.78 s | 191.16 s | 187.64 s | derived (sum of the measured stage times above, one process) |
| S1-S10 in one run | - | 197.86 s | 194.86 s | derived (same) |

S8 and S9 are not separable in time: S8's edits and S9's cuts are applied in one pass per genome (the
orchestrator's decision for this session), so the table gives them together. Against the session-6 note's
*estimate* (S8 and S9 at about 5 to 6 s each, S10 about 1 s, from S7's cost per rewrite) the measured S8-S9 is
about 50% to 75% above 10 to 12 s and S10 about 6 to 7 times its 1 s (more than the 20% of AGENTS §6, reported).
Both estimates priced one rewrite each; the excess is probably the per-gene Python work (S8's candidate search
and the code-4/25 recoding of 28,820 edits, S9's placement of 90,282 genes on contigs, S10's synonymous-site
table), which this session did not profile (*estimate*). S11 (not built) is not measured; with it,
S1-S11 will exceed 3.5 minutes (*estimate*: 198 s measured for S1-S10 plus S11's one write of every source, about
5 s at S7's cost per rewrite).

## Realized values, S8 (measured, runs A, B and C identical)

| quantity | value |
|---|---|
| genetic codes | code 4: g031, g038, g039 (one clade, size in [2, 6]); code 25: g015; the other 42 genomes code 11 |
| call errors (§5.1 about 200) | 30 frameshift insertions, 30 deletions, 30 upstream and 30 downstream alternative starts, 40 smORFs, 40 internal stops |
| alt_code rows (genes holding an in-frame TGA in the code-4/25 genomes) | 6,551 |
| edits (the lift from S7) | 28,820: 30 + 30 frameshifts, 40 internal stops, 30 start recodes, 88 bases for the 30 upstream starts, 40 smORF insertions, 26,039 recoded codons' bases, 2,189 recoded stops, 334 -4 overlaps made -1 |
| internal Trp encoded TGA in the code-4 genomes | 0.898 (preset 0.9; codons sharing a base with another gene are not recoded) |
| genes in the S8 frame | 90,282 (S7's 90,242 and 40 smORFs); every one translates to its protein_uid under its genome's code (frameshifts: the intended protein, with the indel undone; pseudogenes: an in-frame stop) |
| call-error status after S9 | 6,080 whole on an emitted contig, 268 cut by a break, 403 on dropped contigs (alt_code rows included) |

**pyrodigal 3.7.1** (the external caller; trained per genome on the S8 frame): on g024 (code 11, the genome with the most call errors, 10) none of
the 5 frameshifted or internal-stop genes is called whole, and all 5 have a call that ends at the true stop (the
3' piece); 2 of 4 downstream and 1 of 1 upstream alternative starts are called with the true start. On g031
(code 4), of the 1,772 genes holding an internal TGA, a call ending at the true stop exists for 0.432 with
translation table 11 and 0.995 with table 4.

## Realized values, S9 (measured)

| quantity | value |
|---|---|
| assembly classes (§5.1 40/40/20%) | 18 complete, 18 draft, 10 MAG |
| contigs | 2,431: 2,110 emitted, 92 under 500 bp, 229 MAG dropout; 36 complete circular, each with a gene across its origin |
| repeat breaks | 764 (copies of a vehicle chassis family, or equal-length copies >= 1 kb at >= 98% identity) |
| N50 / target over draft and MAG genomes | geometric mean 1.206, range 0.884 to 1.708 |
| reversed contigs | 1,269 of 2,431 |
| S7 genes on the contigs | 81,306 whole and reading their S7 bases, 1,842 cut by a break, 1,924 on dropped contigs, 7,094 changed by S8 (their call error or recoding) |
| S5 features through four lifts | elements 171,776 whole, 21 cut, 4,128 on dropped contigs, 75 gone with S7's fusions; IGRs 82,125 whole, 176 cut, 70 changed by an S8 upstream start or smORF, 34 gone with S7; motifs 1,078 whole |
| S6 features | spacer slots 3,462 whole; arrays 24 whole, 6 recoded (code 4/25); prophages 9 whole, 7 cut |
| S8-S9 digest | `2f425158dc38f2302c1d02e9ebed30da`, byte-identical in the child (another hash seed, one thread) |
| earlier stage digests, before and after S8-S10 | backbone `025fdd0d27fdc65e678782216ecd21fb`, sequences `34ebce1324413a261acd874d7ee4fe3e`, replicons `55bda97155d75148e5cdca7b90a1d7be`, mobile `0bd7f0db5691de58ea5e1a0df1e6e570`, fusions `db09abf3fe79ae34a73e3ab0664a9b13` |

## Realized values, S10 (measured, run C)

| quantity | value |
|---|---|
| contaminated hosts (§5.1: 3) and donor | g013, g016, g030, each from g000 (the farthest genome from each on the species tree) |
| contaminant contigs | g013: 3 contigs, 47,303 bp, GC 0.466 (host contigs median 0.531); g016: 4, 25,462 bp, GC 0.461 (0.517); g030: 3, 47,101 bp, GC 0.476 (0.441) |
| mock FCS-GX sidecar | flags exactly the 10 contaminant contigs, EXCLUDE, whole length |
| duplicates (§5.1: 2 exact, 2 near) | 2 exact (equal contig_uid per contig); 2 near, g008 at ANI 0.999001 and g031 (a code-4 genome) at 0.999000 over their emitted contigs, from 1,574 intergenic and 1,561 synonymous changes; 2 contigs lie wholly inside a gene cut at both edges and have no site to change (run D) |
| S10 digest | `770007bde8615f41d9dd5f2eb912e753` (runs C and D), byte-identical in the child |

Run B found a bug before run C: a near duplicate could change a base of a gene cut by a contig edge as if it were
intergenic (the coverage counted whole gene pieces only). Fixed in `sources.synonymous_sites` (every gene piece
counts for coverage) and covered by the R1 oracle; run B's S10 digest `0d773cbe0ff67c352034d60c6ba95bf1` is
therefore not the final one.

Run C's two failures were in the tests. (1) `tests/r1/test_synthworld_fusions_units.py::
test_phage_world_preset_leaves_the_later_sections_out` lists the preset sections later than S6. The four new
sections are left out of the phage world as they should be, but the list did not name them; it now does and
asserts each is absent (migrated, not loosened). (2) The near-duplicate test asked every contig for
round(0.001 L) changes. A contig lying wholly inside a gene cut at both edges has no site (2 such contigs). The
test now requires exactly that condition when a contig has fewer changes, and the docstring says so. Run D (all
40 tests of the changed files) passed: S8-S9 18.31 s in the fixture and 17.38 s / 183 MB in a child; S10 6.67 s
and 6.59 s / 128 MB; S1-S10 209.31 s (derived, sum).
