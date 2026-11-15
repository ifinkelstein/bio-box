#!/bin/bash
# Card W-32: the comparisons of run2.sh's outputs (as compare.sh does for run.sh's). ON THE POD HOST, after run2.sh.
set -uo pipefail
W=/tmp/if2453/w32
T7=c1bf639f92687b483b659b4dadf7c52c
EXR=$W/tree2/examples/rnap/figures
mkdir -p "$W/sums" "$W/figs/rnap"

sums() {
  (cd "$1" && find . -type f ! -name figure.svg ! -name figure.pdf | LC_ALL=C sort | xargs -d '\n' sha256sum) \
    > "$W/sums/$2.sha256"
  wc -l < "$W/sums/$2.sha256"
}

n_new=$(sums "$W/rnap2/flanks" rnap2.new)
n_base=$(sums "$W/base_out/rnap2" rnap2.base)
if diff -q <(grep -v ' ./README.txt$' "$W/sums/rnap2.new.sha256") <(grep -v ' ./README.txt$' "$W/sums/rnap2.base.sha256") \
  > /dev/null; then
  echo "rnap2: tables identical ($((n_new - 1)) files besides README.txt and the figures, sha256)"
else
  echo "rnap2: TABLES DIFFER"
fi
cmp -s "$W/rnap2/flanks/README.txt" "$W/base_out/rnap2/README.txt" && echo "rnap2: README.txt identical" \
  || echo "rnap2: README.txt differs ($(wc -l < "$W/base_out/rnap2/README.txt") lines before, $(wc -l < "$W/rnap2/flanks/README.txt") now)"
echo "rnap2: $n_base and $n_new files with README.txt; figures $(find "$W/rnap2/flanks" -name 'figure.*' | wc -l)"
same=0
differ=0
for f in $(cd "$W/rnap2/flanks" && find . -name 'figure.*' | LC_ALL=C sort); do
  if cmp -s "$W/rnap2/flanks/$f" "$W/base_out/rnap2/$f"; then same=$((same + 1)); else differ=$((differ + 1)); fi
done
echo "rnap2 new/base: figures identical $same, different $differ"
if diff -r -q "$W/oldflags/rnap2" "$W/base_out/rnap2" > /dev/null; then
  echo "rnap2 old flags/base: every file identical ($(find "$W/oldflags/rnap2" -type f | wc -l) files, figures included)"
else
  echo "rnap2 old flags/base: DIFFER"
  diff -r -q "$W/oldflags/rnap2" "$W/base_out/rnap2" | head -n 20
fi
for x in svg pdf; do
  cmp -s "$W/base_out/rnap2/groups/$T7/figure.$x" "$EXR/t7_group_c1bf639f.$x" \
    && echo "base t7_group.$x: identical to the committed file" || echo "base t7_group.$x: differs from the committed file"
done
if diff -r -q "$W/final/cas7f_all_nr50" "$W/crispr/flanks/cas7f_all_nr50" > /dev/null; then
  echo "cas7f_all_nr50 final/run.sh: every file identical ($(find "$W/final/cas7f_all_nr50" -type f | wc -l) files)"
else
  echo "cas7f_all_nr50 final/run.sh: DIFFER"
fi
cp "$W/rnap2/flanks/groups/$T7/figure.svg" "$W/figs/rnap/t7_group_c1bf639f.svg"
cp "$W/rnap2/flanks/groups/$T7/figure.pdf" "$W/figs/rnap/t7_group_c1bf639f.pdf"
cp "$W/final/cas7f_all_nr50/groups/all/figure.svg" "$W/figs/crispr/cas7f_all_nr50.svg"
cp "$W/final/cas7f_all_nr50/groups/all/figure.pdf" "$W/figs/crispr/cas7f_all_nr50.pdf"
sha256sum "$W"/figs/crispr/* "$W"/figs/rnap/*
grep -v '^ok' "$W/rnap2/check.txt"
cat "$W/logs/rc.txt"
