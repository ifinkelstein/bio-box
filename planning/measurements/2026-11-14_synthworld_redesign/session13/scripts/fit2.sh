#!/bin/sh
# fit2 (card F-02 session 13, DEC-48 (1)): after fit1 (partners at the top level: the partner's share of the other
# lineages 0.98-0.99 at the default seed, above the decoy test's 0.9, and decoys still dying before they transfer at
# 6-12 seeds up to rate 1.2), partners at the middle level (0.03) and mixed, and session 12's partners between its
# rates; 4 processes in one unit, one scan at a time.
F=/tmp/if2453/f0213/scripts/fit1.sh
O=/tmp/if2453/f0213/scan2
sh $F $O l1a [4,12,13,21,22] 4 2.0 4.0 6.0
sh $F $O l1b [5,12,14,22,23] 4 2.0 4.0 6.0
sh $F $O mix [6,24,12,21,22] 4 4.0
sh $F $O s12 [1,9,11,19,22] 4 6.0
