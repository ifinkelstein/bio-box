#!/bin/sh
for t in mafft muscle clustalo famsa kalign hmmbuild hmmsearch hmmalign jackhmmer hmmscan hmmpress mmseqs esl-reformat esl-sfetch seqkit; do
  printf "%s " "$t"
  command -v "$t" || echo MISSING
done
