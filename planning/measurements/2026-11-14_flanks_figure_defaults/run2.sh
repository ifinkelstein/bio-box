#!/bin/bash
# Card W-32, second driver: run.sh's RNA polymerase example failed in nbr flanks with the code it froze
# (KeyError 'pfam:S-AdoMet_lyase', fixed in the card's third commit). With the final code: the RNA polymerase example
# again (OUT rnap2), its nbr flanks call with the code before the card and with the old colours on that run, the
# CRISPR example's NR50 figure call again (against run.sh's output), the comparisons, and figure_stats.py on both
# NR50 figures. ON THE POD HOST, detached (unit w32-driver2), read-only on the root; 2 threads, one step at a time.
set -euo pipefail

W=/tmp/if2453/w32
export TREE=$W/tree2
BASESRC=$W/base/src
export PROFILE=$TREE/profiles/wilkcomp02.yaml
export IMAGE=/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif
export DIGEST=sha256:d67300aa1321b571fa8da79e74dfccdd87a8dca2e56518f324902efa6b9f098a
export ROOT=/tmp/if2453/r10/root
export RELEASE=rN.3
export THREADS=2
export BIND=/tmp/if2453
export MEM=16G
ME=w32-driver2
T7=c1bf639f92687b483b659b4dadf7c52c

rm -rf "$TREE"
cp -a /tmp/if2453/ci/W-32 "$TREE"
. "$TREE/examples/_shared/launch.sh"
[ ! -e "$ROOT/STOP" ] || die "a STOP file exists in $ROOT"

run_id() {
  sed -n '/search run: [0-9a-f]\{32\}/{s/.*search run: \([0-9a-f]\{32\}\).*/\1/p;q;}' "$1"
}

bbase() {  # as bb, with the code before the card
  local name=$1
  shift
  local unit="$UNIT-$name-$$"
  local log="$OUT/logs/$name.log"
  room
  bash "$TREE/scripts/pod_run.sh" --profile "$PROFILE" --image "$IMAGE" --digest "$DIGEST" --mem "$MEM" \
    --threads "$THREADS" --bind "$BIND" --bind "$ROOT:$ROOT:ro" --env PYTHONPATH="$BASESRC" \
    --env TMPDIR="$OUT/tmp" --log "$log" "$unit" "$@" > "$OUT/logs/$name.launch"
  until grep -q "^pod_run: $unit exited" "$log"; do
    sleep 5
  done
  grep -q "^pod_run: $unit exited 0$" "$log" || die "$name failed; see $log"
  note "$name: $(grep "^pod_run: $unit usage" "$log" || true)"
}

note "rnap example (final code)"
set +e
OUT=$W/rnap2 UNIT=w32-rnap2 bash "$TREE/examples/rnap/run.sh" > "$W/logs/rnap2.out" 2> "$W/logs/rnap2.err"
echo "rnap2 run.sh exited $?" >> "$W/logs/rc.txt"
set -e
RUN_RNAP=$(run_id "$W/rnap2/logs/search.log")
RUN_CAS7F=$(run_id "$W/crispr/logs/search-cas7f.log")
[ -n "$RUN_RNAP" ] && [ -n "$RUN_CAS7F" ] || die "no run ids"
note "runs: rnap $RUN_RNAP cas7f $RUN_CAS7F"

export OUT=$W/extra2 UNIT=w32-y
mkdir -p "$OUT/logs" "$OUT/tmp"
FC="--root $ROOT --profile $PROFILE --runs $W/crispr/runs --threads $THREADS"
FR="--root $ROOT --profile $PROFILE --runs $W/rnap2/runs --threads $THREADS"
bbase base-rnap2 nbr flanks "$RUN_RNAP" $FR --window 5g --level nr50 --figure --out "$W/base_out/rnap2"
bb old-rnap2 nbr flanks "$RUN_RNAP" $FR --window 5g --level nr50 --figure --colour-weight locus \
  --colour-fallback none --out "$W/oldflags/rnap2"
bb final-cas7f-nr50 nbr flanks "$RUN_CAS7F" $FC --window 6g --level nr50 --group-by none --figure \
  --out "$W/final/cas7f_all_nr50"

# figure_stats.py: the plans of the two NR50 figures, defaults against the old colours
STATS="$TREE/planning/measurements/2026-11-14_flanks_figure_defaults/figure_stats.py"
tool "$OUT/logs/stats-cas7f.time" env PYTHONPATH="$TREE/src" TMPDIR="$OUT/tmp" python "$STATS" "$W/stats/cas7f" all -- \
  "$RUN_CAS7F" $FC --window 6g --level nr50 --group-by none > "$W/stats_cas7f.json"
tool "$OUT/logs/stats-rnap.time" env PYTHONPATH="$TREE/src" TMPDIR="$OUT/tmp" python "$STATS" "$W/stats/rnap" "$T7" -- \
  "$RUN_RNAP" $FR --window 5g --level nr50 > "$W/stats_rnap.json"
note "done"
