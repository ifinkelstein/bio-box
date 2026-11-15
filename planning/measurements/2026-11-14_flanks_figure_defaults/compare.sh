#!/bin/bash
# Card W-32: compare the outputs run.sh wrote (beside this file). ON THE POD HOST, after run.sh. Prints one line per
# comparison; writes sha256 lists under /tmp/if2453/w32/sums and the figures to copy under /tmp/if2453/w32/figs.
set -uo pipefail
W=/tmp/if2453/w32
B=$W/base_out
O=$W/oldflags
EXC=$W/tree/examples/crispr_cas_neighbourhoods/figures
EXR=$W/tree/examples/rnap/figures
T7=c1bf639f92687b483b659b4dadf7c52c
mkdir -p "$W/sums" "$W/figs/crispr" "$W/figs/rnap"

sums() {  # sums DIR NAME: sha256 of every file but the figures, sorted by path
  (cd "$1" && find . -type f ! -name figure.svg ! -name figure.pdf | LC_ALL=C sort | xargs -d '\n' sha256sum) \
    > "$W/sums/$2.sha256"
  echo "$(wc -l < "$W/sums/$2.sha256")"
}

tables() {  # tables NAME NEW BASE: every file but the figures (README.txt by itself), new against base
  local n_new n_base
  n_new=$(sums "$2" "$1.new")
  n_base=$(sums "$3" "$1.base")
  if diff -q <(grep -v ' ./README.txt$' "$W/sums/$1.new.sha256") <(grep -v ' ./README.txt$' "$W/sums/$1.base.sha256") \
    > /dev/null; then
    echo "$1: tables identical ($((n_new - 1)) files besides README.txt and the figures, sha256)"
  else
    echo "$1: TABLES DIFFER"
    diff <(grep -v ' ./README.txt$' "$W/sums/$1.new.sha256") <(grep -v ' ./README.txt$' "$W/sums/$1.base.sha256") \
      | head -n 20
  fi
  if cmp -s "$2/README.txt" "$3/README.txt"; then
    echo "$1: README.txt identical"
  else
    echo "$1: README.txt differs ($(wc -l < "$3/README.txt") lines before, $(wc -l < "$2/README.txt") now)"
  fi
  echo "$1: $n_base and $n_new files (with README.txt) in base and new; figures: $(find "$2" -name 'figure.*' | wc -l)"
}

figures() {  # figures NAME A B: every figure of A against B, byte for byte
  local same=0 differ=0 f
  for f in $(cd "$2" && find . -name 'figure.*' | LC_ALL=C sort); do
    if cmp -s "$2/$f" "$3/$f"; then same=$((same + 1)); else differ=$((differ + 1)); fi
  done
  echo "$1: figures identical $same, different $differ"
}

one() {  # one NAME FILE COMMITTED
  if cmp -s "$2" "$3"; then echo "$1: identical to the committed file"; else echo "$1: differs from the committed file"; fi
}

tables tniq_nr50 "$W/crispr/flanks/tniq_nr50" "$B/tniq_nr50"
tables cas7f_all_nr50 "$W/crispr/flanks/cas7f_all_nr50" "$B/cas7f_all_nr50"
tables tniq_cas_pfam "$W/crispr/flanks/tniq_cas_pfam" "$B/tniq_cas_pfam"
tables cas7f_all_pfam "$W/crispr/flanks/cas7f_all_pfam" "$B/cas7f_all_pfam"
tables rnap "$W/rnap/flanks" "$B/rnap"
figures "cas7f_all_nr50 new/base" "$W/crispr/flanks/cas7f_all_nr50" "$B/cas7f_all_nr50"
figures "tniq_cas_pfam new/base" "$W/crispr/flanks/tniq_cas_pfam" "$B/tniq_cas_pfam"
figures "cas7f_all_pfam new/base" "$W/crispr/flanks/cas7f_all_pfam" "$B/cas7f_all_pfam"
figures "rnap new/base" "$W/rnap/flanks" "$B/rnap"
# the card's code under the old colours against the code before the card: every file, figures and README.txt
for x in cas7f_all_nr50 rnap; do
  if diff -r -q "$O/$x" "$B/$x" > /dev/null; then
    echo "$x old flags/base: every file identical ($(find "$O/$x" -type f | wc -l) files, figures included)"
  else
    echo "$x old flags/base: DIFFER"
    diff -r -q "$O/$x" "$B/$x" | head -n 20
  fi
done
# the code before the card against the committed figures (made on rN.3 by the code of their cards)
one "base cas7f_all_nr50.svg" "$B/cas7f_all_nr50/groups/all/figure.svg" "$EXC/cas7f_all_nr50.svg"
one "base cas7f_all_nr50.pdf" "$B/cas7f_all_nr50/groups/all/figure.pdf" "$EXC/cas7f_all_nr50.pdf"
one "base tniq_cas_pfam.svg" "$B/tniq_cas_pfam/groups/all/figure.svg" "$EXC/tniq_cas_pfam.svg"
one "base cas7f_all_pfam.svg" "$B/cas7f_all_pfam/groups/all/figure.svg" "$EXC/cas7f_all_pfam.svg"
one "base t7_group.svg" "$B/rnap/groups/$T7/figure.svg" "$EXR/t7_group_c1bf639f.svg"
one "base t7_group.pdf" "$B/rnap/groups/$T7/figure.pdf" "$EXR/t7_group_c1bf639f.pdf"
# the new figures, named as the examples name them
for x in cas7f_all_nr50 cas7f_all_pfam tniq_cas_pfam; do
  cp "$W/crispr/flanks/$x/groups/all/figure.svg" "$W/figs/crispr/$x.svg"
  cp "$W/crispr/flanks/$x/groups/all/figure.pdf" "$W/figs/crispr/$x.pdf"
done
cp "$W/rnap/flanks/groups/$T7/figure.svg" "$W/figs/rnap/t7_group_c1bf639f.svg"
cp "$W/rnap/flanks/groups/$T7/figure.pdf" "$W/figs/rnap/t7_group_c1bf639f.pdf"
for f in "$W"/figs/crispr/* "$W"/figs/rnap/*; do
  sha256sum "$f"
done
# the examples' checks
tail -n 3 "$W/crispr/check.txt"
tail -n 3 "$W/rnap/check.txt"
cat "$W/logs/rc.txt"
