#!/bin/bash
# Card W-31 measurement driver. Run ON THE POD HOST (launched as unit w31-driver), as the CRISPR example's run.sh does:
# every bio-box command through scripts/pod_run.sh (a capped unit, nice 19, 2 threads), every python script in a
# capped scope; it waits for room under the account's thread cap and stops on a STOP file in the root.
#
#   TREE=/tmp/if2453/w31/tree BASE=/tmp/if2453/w31/base OUT=/tmp/if2453/w31/out bash run.sh
#
# TREE: the card's code; BASE: the code before the card (src) with the card's tests; OUT: new scratch directory.
# Read only: the root, /tmp/if2453/ex2 (the CRISPR example's runs and outputs) and /tmp/if2453/w25 (the RNA
# polymerase run of planning/measurements/2026-11-12_flanks_rN3.md); both are bound read-only into the image.
# STEPS (default "1 2 3 4 5") runs only the steps named; OUT must be new when step 1 runs. The note's run: all five
# steps, as "1 2 3 4" (unit w31-driver, which stopped at step 3 on a bug of measure.py since fixed, then unit
# w31-driver2 with steps 3 and 4, then called SKIP_PLAIN=1), then STEPS=5 (unit w31-driver3, and again as unit
# w31-driver4 after the figure's legend and title texts of the options were shortened; step 5's outputs removed first;
# that run predates the `rm -f "$log"` in bb, so its three units found the old logs and ran at the same time, 6
# threads for about 40 s).
set -euo pipefail

: "${TREE:?}" "${BASE:?}" "${OUT:?}"
ROOT=/tmp/if2453/r10/root
RELEASE=rN.3
IMAGE=/tmp/if2453/images/bio-box-test_64f1bf0.sif
DIGEST=sha256:3d77df6b7f0b894b5371e42eaf077dfda09ebc97b1f78f2ed1d2cfc5fa9fbda4
PROFILE="$TREE/profiles/wilkcomp02.yaml"
EX=/tmp/if2453/ex2/rN.3-b
THREADS=2
MEAS="$TREE/planning/measurements/2026-11-13_flanks_context_types"
RO=(--bind "$ROOT:$ROOT:ro" --bind /tmp/if2453/ex2:/tmp/if2453/ex2:ro --bind /tmp/if2453/w25:/tmp/if2453/w25:ro)

die() { echo "run.sh: $*" >&2; exit 2; }
note() { echo "run.sh: $(date '+%F %T') $*" >&2; }

[ ! -e "$ROOT/STOP" ] || die "a STOP file exists in $ROOT"
step() { case " ${STEPS:-1 2 3 4 5} " in *" $1 "*) return 0 ;; esac; return 1; }
! step 1 || [ ! -e "$OUT" ] || die "OUT $OUT exists"
mkdir -p "$OUT/logs" "$OUT/tmp" "$OUT/flanks" "$OUT/measure" "$OUT/fixture"

cap() {
  local h d
  h=$(date +%-H)
  d=$(date +%u)
  if [ "$d" -le 5 ] && [ "$h" -ge 8 ] && [ "$h" -lt 20 ]; then echo 40; else echo 56; fi
}

room() {
  local total v c
  while true; do
    [ ! -e "$ROOT/STOP" ] || die "a STOP file appeared in $ROOT"
    total=0
    for v in $(ps -u "$USER" -o pcpu=); do
      total=$((total + ${v%.*}))
    done
    c=$(cap)
    if [ $(((total + 99) / 100 + THREADS)) -le "$c" ]; then
      return
    fi
    note "waiting for room: account at ${total}% CPU, cap $c threads"
    sleep 60
  done
}

bb() {  # bb NAME ARGS...: one bio-box command (the card's tree) in unit w31-NAME; waits for it
  local name=$1
  shift
  local unit="w31-$name"
  local log="$OUT/logs/$name.log"
  room
  rm -f "$log"  # a log left by an earlier run of the step would end the wait below at once
  bash "$TREE/scripts/pod_run.sh" --profile "$PROFILE" --image "$IMAGE" --digest "$DIGEST" --mem 16G \
    --threads "$THREADS" --bind /tmp/if2453 "${RO[@]}" --env PYTHONPATH="$TREE/src" \
    --env TMPDIR="$OUT/tmp" --log "$log" "$unit" "$@" > "$OUT/logs/$name.launch"
  until grep -q "^pod_run: $unit exited" "$log"; do
    sleep 5
  done
  grep -q "^pod_run: $unit exited 0$" "$log" || die "$name failed; see $log"
  note "$name: $(grep "^pod_run: $unit usage" "$log" || true)"
}

tool() {  # tool LOG SRC ARGS...: python of the image with SRC first on the path, in a capped scope
  local log=$1 src=$2
  shift 2
  room
  systemd-run --user --scope --quiet -p MemoryMax=16G -p MemorySwapMax=0 -p CPUQuota=$((THREADS * 100))% \
    nice -n 19 /usr/bin/time -v -a -o "$log.time" singularity exec --cleanenv --env PYTHONPATH="$src" \
    --env OMP_NUM_THREADS="$THREADS" --env TMPDIR="$OUT/tmp" --bind /tmp/if2453 "${RO[@]}" --pwd "$OUT" "$IMAGE" \
    "$@" > "$log" 2>&1
}

note "tree $TREE, base $BASE, out $OUT"

# 1. Unflagged outputs on the fixture: the code before the card against the card's code.
if step 1; then
  tool "$OUT/logs/fixture-base" "$BASE/src" python "$MEAS/fixture_identity.py" "$BASE" "$OUT/fixture/base"
  tool "$OUT/logs/fixture-new" "$TREE/src" python "$MEAS/fixture_identity.py" "$TREE" "$OUT/fixture/new"
  if diff -r -q "$OUT/fixture/base/runs" "$OUT/fixture/new/runs" > "$OUT/fixture/diff.txt"; then
    note "fixture: unflagged outputs identical ($(find "$OUT/fixture/new/runs" -type f | wc -l) files)"
  else
    note "fixture: DIFFERENCES, see $OUT/fixture/diff.txt"
  fi
fi

# 2. Unflagged outputs on rN.3: the example's three nbr flanks calls with the card's code against the example's
#    outputs (written by the code before the card, image bio-box-test_64f1bf0).
RUN_TNIQ=ef8d86b8e0e4b756a8b1d63d9eb92aab
RUN_CAS7F=47db59881f7779cb86d7e35984b06705
RUN_RNAP=599bea5d471c84759a0eec02e2ca9225
TV=$(sed -n 's/^# tape_version\t"\(.*\)"$/\1/p' "$EX/tniq/neighbors.tsv")
[ -n "$TV" ] || die "no tape_version"
F="--root $ROOT --profile $PROFILE --runs $EX/runs --threads $THREADS"
FR="--root $ROOT --profile $PROFILE --runs /tmp/if2453/w25/runs --threads $THREADS"
TNIQ_NR50=(nbr flanks "$RUN_TNIQ" $F --window 10g --level nr50)
TNIQ_ALL_NR50=(nbr flanks "$RUN_TNIQ" $F --window 10g --level nr50 --group-by none)
CAS7F_NR50=(nbr flanks "$RUN_CAS7F" $F --window 6g --level nr50 --group-by none)
CAS7F_PFAM=(nbr flanks "$RUN_CAS7F" $F --window 6g --level pfam --group-by none)
TNIQ_CAS_PFAM=(nbr flanks --anchors-file "$EX/uids_tniq_cas.txt" --release "$RELEASE" --tape-version "$TV" $F
               --window 10g --level pfam --group-by none)
RNAP_NR50=(nbr flanks "$RUN_RNAP" $FR --window 5g --level nr50)
if step 2; then
  bb plain-tniq "${TNIQ_NR50[@]}" --out "$OUT/flanks/tniq_nr50"
  bb plain-cas7f "${CAS7F_NR50[@]}" --figure --out "$OUT/flanks/cas7f_all_nr50"
  bb plain-tniq-cas "${TNIQ_CAS_PFAM[@]}" --figure --out "$OUT/flanks/tniq_cas_pfam"
  for d in tniq_nr50 cas7f_all_nr50 tniq_cas_pfam; do
    if diff -r -q "$EX/flanks/$d" "$OUT/flanks/$d" > "$OUT/flanks/$d.diff"; then
      note "rN.3 $d: identical to the example's output ($(find "$OUT/flanks/$d" -type f | wc -l) files)"
    else
      note "rN.3 $d: DIFFERENCES, see $OUT/flanks/$d.diff"
    fi
  done
fi

# 3. The comparison, in-process, per run.
M() { tool "$OUT/logs/measure-$1" "$TREE/src" python "$MEAS/measure.py" "$OUT/measure" "$1" "$2" -- "${@:3}"; }
if step 3; then
  M cas7f_all_nr50 Cas7f "${CAS7F_NR50[@]}"
  M cas7f_all_pfam Cas7f "${CAS7F_PFAM[@]}"
  M tniq_cas_pfam TniQ "${TNIQ_CAS_PFAM[@]}"
  M tniq_nr50 TniQ "${TNIQ_NR50[@]}"
  M tniq_all_nr50 TniQ "${TNIQ_ALL_NR50[@]}"
  M rnap_nr50 other "${RNAP_NR50[@]}"
  note "measure: done"
fi

# 4. Figures and tables with the flags, through the CLI.
if step 4; then
  bb fl-cas7f-gap "${CAS7F_NR50[@]}" --figure --auto-cut gap --out "$OUT/flanks/cas7f_all_nr50_gap"
  bb fl-cas7f-nr90 "${CAS7F_NR50[@]}" --figure --colour-weight anchor-nr90 --out "$OUT/flanks/cas7f_all_nr50_w_nr90"
  bb fl-cas7f-genome "${CAS7F_NR50[@]}" --figure --colour-weight genome --out "$OUT/flanks/cas7f_all_nr50_w_genome"
  bb fl-cas7f-pfam "${CAS7F_NR50[@]}" --figure --colour-fallback pfam --out "$OUT/flanks/cas7f_all_nr50_fb_pfam"
  bb fl-cas7f-both "${CAS7F_NR50[@]}" --figure --colour-weight anchor-nr90 --colour-fallback pfam --auto-cut gap \
    --out "$OUT/flanks/cas7f_all_nr50_all3"
  bb fl-cas7f-pfam-gap "${CAS7F_PFAM[@]}" --figure --auto-cut gap --out "$OUT/flanks/cas7f_all_pfam_gap"
  bb fl-tniq-cas-gap "${TNIQ_CAS_PFAM[@]}" --figure --auto-cut gap --out "$OUT/flanks/tniq_cas_pfam_gap"
  bb fl-tniq-all-gap "${TNIQ_ALL_NR50[@]}" --auto-cut gap --out "$OUT/flanks/tniq_all_nr50_gap"
  bb fl-rnap-gap "${RNAP_NR50[@]}" --auto-cut gap --out "$OUT/flanks/rnap_nr50_gap"
fi

# 5. The merge distances of the four one-group trees, and the figures of the note: both colour options at the
#    fixed cut, and with --auto-cut 10.
if step 5; then
  tool "$OUT/logs/heights" "$TREE/src" python -I "$MEAS/heights.py" \
    "$OUT"/measure/{cas7f_all_nr50,cas7f_all_pfam,tniq_cas_pfam,tniq_all_nr50}/cli_plain/groups/all/tree.nwk
  cp "$OUT/logs/heights" "$OUT/heights.tsv"
  bb fl-cas7f-fb-nr90 "${CAS7F_NR50[@]}" --figure --colour-weight anchor-nr90 --colour-fallback pfam \
    --out "$OUT/flanks/cas7f_all_nr50_fb_pfam_w_nr90"
  bb fl-cas7f-n10 "${CAS7F_NR50[@]}" --figure --colour-weight anchor-nr90 --colour-fallback pfam --auto-cut 10 \
    --out "$OUT/flanks/cas7f_all_nr50_n10_fb_pfam_w_nr90"
  bb fl-tniq-cas-n10 "${TNIQ_CAS_PFAM[@]}" --figure --auto-cut 10 --out "$OUT/flanks/tniq_cas_pfam_n10"
fi
note "done"
