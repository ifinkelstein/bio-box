# synthworld S6 (phages, prophages, CRISPR arrays) on `micro`: time, memory and realized values

Card F-02, session 5 (branch `task/F-02-5`). Code: `src/biobox/synthworld/mobile.py` (its docstring states the
model), `backbone.build_phages`, the `phages`, `prophages` and `crispr` sections of `presets/micro.yaml`.

Host and command for every *measured* number below: wilkcomp02, 2026-11-13, one process at nice 19 inside
`/tmp/if2453/images/bio-box-test_64f1bf0.sif` (image of `64f1bf0`, which carries zombi2 0.54.0), with the account
at about 30 to 35 threads of its 40-thread day cap (Pfam scans and other agents' tests);
`env BIOBOX_POD_IMAGE=/tmp/if2453/images/bio-box-test_64f1bf0.sif BIOBOX_POD_PYDEPS=/tmp/if2453/ci/F-02-5/src
sh scripts/pod_pytest.sh F-02-5 -s tests/r2/test_synthworld_mobile.py`, `micro` default seed 20261112.

## Time and memory

| quantity | value | label |
|---|---|---|
| S6 wall time, `mobile.build` in the test fixture | 10.23 s, 10.95 s, 10.84 s (three runs) | measured |
| S6 wall time in a child process (`PYTHONHASHSEED=4242`, one thread) | 11.85 s, 11.46 s, 11.75 s | measured |
| S6 peak RSS in that child (VmHWM) | 141 MB (first phage setting), 166 MB, 166 MB | measured |
| time left for S6-S11 before this session | about 23 s | derived (session 4: 180 s target - about 157 s for S1-S5) |
| time left for S7-S11 after S6 | about 12 s | derived (23 s - about 11 s) |

No earlier document gives an estimate for S6 alone, so there is no 20% comparison to make; S6 takes about half
of what was left for S6 to S11 (derived), and the four stages still to come have about 12 s (derived). The
largest parts of S6 are rewriting the host replicons (about 100 Mb read and written again) and the phage set's
own S3 to S5.

## Phage set (the nested world `mobile/phages/`)

Backbone rates were scanned on the default seed (8 settings, backbone only, through a temporary test that was
not committed; measured): duplication, transfer and loss at 0.01, 0.05, 0.1 per copy gave 12 phages of 31 to 63
genes; halved (0.005, 0.025, 0.05; origination 0.015 x 60 per lineage, derived to balance) they give 49 to 71
genes. Other settings scanned: (loss, duplication, transfer, root families) (0.05, 0.005, 0.025, 64) 53-89 genes;
(0.03, 0.003, 0.015, 60) 57-75; (0.03, 0.003, 0.015, 64) 57-83; (0.02, 0.002, 0.01, 62) 52-65 with every family
shared; (0.05, 0.005, 0.025, 66) 52-86; (0.03, 0.003, 0.015, 66) 55-78.

With the chosen setting (measured): genes per phage 49, 49, 50, 51, 58, 59, 66, 68, 68, 70, 71, 71 (§5.1: 40-80);
phage lengths 54,249 to 74,726 bp; 70 families, 0.87 of them in two or more phages; strand switches about 30
per 100 genes; every phage gene translates to its protein; the phage set's HMM panel file is byte-identical to
the host's.

## Prophages and CRISPR-Cas loci

| quantity | value | label |
|---|---|---|
| prophages, instances | 4 prophages in 16 host genomes | measured |
| prophage identity to its phage | 0.9682 to 0.9987 | measured |
| arrays, instances | 10 arrays in 30 host genomes | measured |
| spacers per array instance | 8 to 22, median 16 (§5.1: 5-40) | measured |
| spacers acquired (all, lost ones included) | 316: 16 decoys; PAM present in 0.927 of the phage-derived; mismatches 0/1/2/3/4/5 in 142/62/23/39/21/13 | measured |
| spacer x target rows | 3,462: phage 1,373, prophage 2,070, none 19 | measured |
| brute-force check of the target set | agrees on 73 spacers (680 targets) | measured |
| host genes, S5 elements, motifs and IGRs that read the same after the insertions | 89,102 genes and 255,355 features after 46 insertions | measured |
| stage digest | `0bd7f0db5691de58ea5e1a0df1e6e570`, the same in every run and in the child | measured |
| chance that a random 30-mer lies within 5 mismatches of a given one | 3.0e-11 | derived: sum over k <= 5 of C(30, k) 3^k / 4^30 |
| expected chance hits of one spacer in the host outside prophages | about 6e-3 | derived: 3.0e-11 x 2e8 positions (about 100 Mb, both strands) |
