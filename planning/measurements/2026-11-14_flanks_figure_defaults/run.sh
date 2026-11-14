#!/bin/bash
# Card W-32: the two worked examples on rN.3 with the card's code, the same nbr flanks calls with the code before the
# card (src of 4469fab), and the comparisons (planning/measurements/2026-11-14_flanks_figure_defaults.md). Run ON THE
# POD HOST as a detached unit (w32-driver), read-only on the root:
#
#   systemd-run --user --unit w32-driver --collect -p MemoryMax=2G -p MemorySwapMax=0 \
#     -p StandardOutput=append:/tmp/if2453/w32/driver.log -p StandardError=append:/tmp/if2453/w32/driver.log \
#     --setenv=USER=if2453 /usr/bin/nice -n 19 bash /tmp/if2453/w32/scripts/run.sh
#
# Every bio-box step runs through scripts/pod_run.sh at 2 threads (THREADS), one at a time, in the runtime image
# bio-box_faec1f4.sif by digest. Writes only under /tmp/if2453/w32.
set -euo pipefail

W=/tmp/if2453/w32
export TREE=$W/tree                       # a frozen copy of the card's tree (mirrored by scripts/pod_pytest.sh)
BASESRC=$W/base/src                       # src of 4469fab
export PROFILE=$TREE/profiles/wilkcomp02.yaml
export IMAGE=/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif
export DIGEST=sha256:d67300aa1321b571fa8da79e74dfccdd87a8dca2e56518f324902efa6b9f098a
export ROOT=/tmp/if2453/r10/root
export RELEASE=rN.3
export THREADS=2
export BIND=/tmp/if2453
export MEM=16G
ME=w32-driver

rm -rf "$TREE"
cp -a /tmp/if2453/ci/W-32 "$TREE"
mkdir -p "$W/logs"
. "$TREE/examples/_shared/launch.sh"
[ ! -e "$ROOT/STOP" ] || die "a STOP file exists in $ROOT"

run_id() {
  sed -n '/search run: [0-9a-f]\{32\}/{s/.*search run: \([0-9a-f]\{32\}\).*/\1/p;q;}' "$1"
}

bbase() {  # bbase NAME ARGS...: as bb, with the code before the card
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

# 1. The examples, with the card's code (each run.sh's own steps and checks).
note "crispr example"
set +e
OUT=$W/crispr UNIT=w32-crispr bash "$TREE/examples/crispr_cas_neighbourhoods/run.sh" > "$W/logs/crispr.out" \
  2> "$W/logs/crispr.err"
echo "crispr run.sh exited $?" >> "$W/logs/rc.txt"
note "rnap example"
OUT=$W/rnap UNIT=w32-rnap bash "$TREE/examples/rnap/run.sh" > "$W/logs/rnap.out" 2> "$W/logs/rnap.err"
echo "rnap run.sh exited $?" >> "$W/logs/rc.txt"
set -e

RUN_TNIQ=$(run_id "$W/crispr/logs/search-tniq.log")
RUN_CAS7F=$(run_id "$W/crispr/logs/search-cas7f.log")
RUN_RNAP=$(run_id "$W/rnap/logs/search.log")
TV=$(sed -n 's/^# tape_version\t"\(.*\)"$/\1/p' "$W/crispr/tniq/neighbors.tsv")
[ -n "$RUN_TNIQ" ] && [ -n "$RUN_CAS7F" ] && [ -n "$RUN_RNAP" ] && [ -n "$TV" ] || die "no run ids"
note "runs: tniq $RUN_TNIQ cas7f $RUN_CAS7F rnap $RUN_RNAP"

# 2. The extra figure of the CRISPR example (its README, section 2), with the card's code.
export OUT=$W/extra UNIT=w32-x
mkdir -p "$OUT/logs" "$OUT/tmp"
FC="--root $ROOT --profile $PROFILE --runs $W/crispr/runs --threads $THREADS"
FR="--root $ROOT --profile $PROFILE --runs $W/rnap/runs --threads $THREADS"
bb new-cas7f-pfam nbr flanks "$RUN_CAS7F" $FC --window 6g --level pfam --group-by none --figure \
  --out "$W/crispr/flanks/cas7f_all_pfam"

# 3. Every nbr flanks call of the two examples with the code before the card, and the two nr50 calls with the
#    card's code under the old colours (--colour-weight locus --colour-fallback none).
B=$W/base_out
bbase base-tniq-nr50 nbr flanks "$RUN_TNIQ" $FC --window 10g --level nr50 --out "$B/tniq_nr50"
bbase base-cas7f-nr50 nbr flanks "$RUN_CAS7F" $FC --window 6g --level nr50 --group-by none --figure \
  --out "$B/cas7f_all_nr50"
bbase base-tniq-cas-pfam nbr flanks --anchors-file "$W/crispr/uids_tniq_cas.txt" --release "$RELEASE" \
  --tape-version "$TV" $FC --window 10g --level pfam --group-by none --figure --out "$B/tniq_cas_pfam"
bbase base-cas7f-pfam nbr flanks "$RUN_CAS7F" $FC --window 6g --level pfam --group-by none --figure \
  --out "$B/cas7f_all_pfam"
bbase base-rnap nbr flanks "$RUN_RNAP" $FR --window 5g --level nr50 --figure --out "$B/rnap"
O=$W/oldflags
bb old-cas7f-nr50 nbr flanks "$RUN_CAS7F" $FC --window 6g --level nr50 --group-by none --figure \
  --colour-weight locus --colour-fallback none --out "$O/cas7f_all_nr50"
bb old-rnap nbr flanks "$RUN_RNAP" $FR --window 5g --level nr50 --figure --colour-weight locus \
  --colour-fallback none --out "$O/rnap"

# 4. The comparisons (sha256 lists, then compare.sh's lines) and the figures to copy.
bash "$W/scripts/compare.sh" > "$W/compare.txt" 2>&1 || true
note "done"
