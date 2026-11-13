#!/bin/bash
# Card 0-08b: put first on PATH as `cmsearch` (inside the image) so that each cmsearch call the producer makes is
# timed by timed.py (CPU-minutes, wall, peak RSS) into $M08B_REC/<label>-<igr|whole>-<ns>.json, and its `-o
# /dev/null` becomes `-o <same stem>.out` so each model's own "# CPU time" line is kept (as card 0-09c did). Every
# other argument is the producer's (`biobox.l1.ncrna._cmsearch_argv`). A call without --tblout (`cmsearch -h`, read
# for the version) goes to the binary untouched.
REAL=/opt/conda/bin/cmsearch
case " $* " in *" --tblout "*) ;; *) exec "$REAL" "$@" ;; esac
D=${M08B_REC:?M08B_REC}
mkdir -p "$D"
kind=whole
for a in "$@"; do [ "$a" = --notrunc ] && kind=igr; done
base="$D/${M08B_LABEL:-x}-$kind-$(date +%s%N)"
args=()
prev=""
for a in "$@"; do
  if [ "$prev" = "-o" ] && [ "$a" = /dev/null ]; then args+=("$base.out"); else args+=("$a"); fi
  prev=$a
done
exec python "${M08B_TIMED:?M08B_TIMED}" "$base.json" -- "$REAL" "${args[@]}"
