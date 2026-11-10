#!/bin/sh
# A shadow data root under /tmp/if2453/ssrnap/shadow: symlinks to the read-only durable copy of
# the pilot store (provisional/rS1.pilot1, with Pfam 38.2 l1/domains), so pack.build and
# runner.run (which insist on writing under <root>/search_pack and <root>/search) write only here.
set -e
R=/stor/scratch/Wilke/if2453/all-the-data/provisional/rS1.pilot1
S=/tmp/if2453/ssrnap/shadow
mkdir -p "$S/l1"
for d in l0 l0seq l2 l3; do ln -sfn "$R/$d" "$S/$d"; done
for d in protein_occurrence slot_protein domains; do ln -sfn "$R/l1/$d" "$S/l1/$d"; done
cat /tmp/if2453/ssrnap/hmm/RNA_pol.hmm /tmp/if2453/ssrnap/hmm/RPOL_N.hmm /tmp/if2453/ssrnap/hmm/ssRNAP_seed8.hmm > /tmp/if2453/ssrnap/hmm/query_lib.hmm
ls -la "$S" "$S/l1"
