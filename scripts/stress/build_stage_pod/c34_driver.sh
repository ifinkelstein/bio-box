#!/bin/bash
# Card C-34: the label files of stages S1, P1 and P2 of the data root (ingested before card B-27), each stage one
# `run-step` process of `bio-box stage labels` through scripts/pod_run.sh, with the runtime image of `faec1f4` and the
# tree's src (the form of card C-33's c33_driver.sh). Host side, bash. Copies of this script, c34_view.py, c34_check.py
# and c34_oracle.py run from /tmp/if2453/c34/bin, and the tree is a copy in /tmp/if2453/c34/tree (the card's code) and
# /tmp/if2453/c34/oldtree (the code before the card, src only), so nothing synced later changes them. Each mode is
# launched in its own small unit, for example
#
#   systemd-run --user --unit=c34-driver-rehearse --collect -p MemoryMax=4G -p MemorySwapMax=0 \
#     -p StandardOutput=append:/tmp/if2453/c34/driver.out -p StandardError=append:/tmp/if2453/c34/driver.out \
#     nice -n 19 bash /tmp/if2453/c34/bin/c34_driver.sh rehearse
#
# Modes, in the order the card runs them:
#   rehearse  A scratch root C/rehearsal (c34_view.py rehearsal: links to the contig_occurrence and ingest ledger of
#             S1, P1, P2; its own runs/, ledger/stage=<S>/, l0seq/contig_label), a scratch profile with its own tmp,
#             and the exact step command of the real run for each stage with --root the scratch root (units
#             c34-rh-<s>). Then a view of the whole root with those label files in place (c34_view.py verify,
#             C/verifyview) and `release verify` of every release with a manifest there, by the code before the card
#             (C/oldtree) and by the card's (C/tree): C/rehearsal-verify/{old,new}-<R>.log. Writes nothing in ROOT.
#   oracle    A re-ingest of S1's inputs (the objects.json its ingest run record names, read-only) into C/oracle/root
#             with the card's code, `stage ingest` in unit c34-oracle-s1; then c34_oracle.py: the label files against
#             the rehearsal's S1 files, and the re-ingest's contig_occurrence against the root's, by sha256. Writes
#             C/oracle/OK when they are the same.
#   before    Refuses when ROOT/STOP exists, a build unit is active, or C/rehearsal/OK or C/oracle/OK is missing;
#             `release verify` of every release with a manifest (card's code, root bound read-only), and the snapshot
#             into C/before: listing.tsv (path, bytes, mtime of every file of ROOT) and stages.sha256 (every file under
#             a stage=S1, stage=P1 or stage=P2 directory). Writes C/before/DONE when every release verifies.
#   run       Needs C/before/DONE. `stage labels` on S1, P1, P2, one unit at a time (c34-<s>), T threads each (DuckDB
#             `--cpus T`, run-step `--threads T`, pod_run `--threads T --mem 8G`). Before each launch: ROOT/STOP (stop),
#             the step's run record (skip) or intent file (stop), and the account's allowance for a T-thread launch
#             (`python -m biobox.runstep`, card C-28; wait while it is under T). A unit that does not end `exited 0`
#             with its run record published stops the run; nothing in the root is cleaned. Writes C/run/DONE.
#   after     The same snapshot into C/after, `release verify` of every release again, `bio-box releases`, and the
#             lookups `find --accession NC_001604` and `NZ_AP026781.1` on rM.2 and rN.3 (also run in before, to
#             compare byte for byte).
#   check     c34_check.py in the image: before against after, the new files against the rehearsal's, and
#             C/check/root_added.tsv (copied to C/root_added.tsv).
#   seven     A view of the root with its own release_def (c34_view.py seven, C/view7), the release rC34.v defined
#             there over the seven stages, and `bio-box releases` of the view. No `release define` in ROOT.
set -u

C=/tmp/if2453/c34
ROOT=/tmp/if2453/r10/root
TREE=$C/tree
OLD=$C/oldtree
PROFILE=$TREE/profiles/wilkcomp02.yaml
IMAGE=/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif
DIGEST=sha256:d67300aa1321b571fa8da79e74dfccdd87a8dca2e56518f324902efa6b9f098a
BIN=$(cd "$(dirname "$0")" && pwd)
LOGS=$C/logs
LOG=$C/driver.log
STAGES=(S1 P1 P2)
SEVEN=S1,P1,P2,P3,P4,P5,P6
RELEASES=(rN.1a rN.1b rN.1 rN.2 rN.3 rM.2)
T=${T:-2}
CPU_MAX=${CPU_MAX:-60}
mkdir -p "$LOGS" "$C/scratch"

log() { echo "$(date '+%F %T %Z') $*" >> "$LOG"; }
fail() { log "FAILED: $*"; mkdir -p "$C/$MODE"; echo "$*" > "$C/$MODE/FAILED"; exit 1; }

bb() {  # a short command in the image, in its own capped scope at nice 19, the root bound read-only; $1 tree
  local tree=$1 what=$2
  shift 2
  case "$what" in
    cli) set -- python -m biobox.cli "$@" ;;
    py) set -- python "$@" ;;
  esac
  systemd-run --user --scope --quiet -p MemoryMax=16G -p MemorySwapMax=0 -p CPUQuota=200% nice -n 19 \
    singularity exec --cleanenv --env PYTHONPATH="$tree/src" --env OMP_NUM_THREADS=2 --env PYTHONUNBUFFERED=1 \
    --bind /tmp/if2453 --bind "$ROOT:$ROOT:ro" --bind /sys/fs/cgroup:/sys/fs/cgroup:ro "$IMAGE" "$@"
}

lower() { echo "$1" | tr 'A-Z' 'a-z'; }

snapshot() {  # $1 directory
  local d=$1
  mkdir -p "$d"
  find "$ROOT" ! -type d -printf '%P\t%s\t%T@\n' | LC_ALL=C sort > "$d/listing.tsv"
  find "$ROOT" ! -type d \( -path '*/stage=S1/*' -o -path '*/stage=P1/*' -o -path '*/stage=P2/*' \) -printf '%P\n' \
    | LC_ALL=C sort > "$d/stages.files"
  (cd "$ROOT" && systemd-run --user --scope --quiet -p MemoryMax=2G -p MemorySwapMax=0 nice -n 19 \
    xargs -d '\n' sha256sum < "$d/stages.files") > "$d/stages.sha256"
  log "snapshot $d: $(wc -l < "$d/listing.tsv") files, $(wc -l < "$d/stages.sha256") of S1, P1, P2 hashed, $(awk -F '\t' '{s += $2} END {printf "%.0f", s}' "$d/listing.tsv") bytes"
}

verify_all() {  # $1 tree, $2 root, $3 directory, $4 label: every release; returns 1 if any fails
  local tree=$1 root=$2 d=$3 tag=$4 r rc bad=0
  mkdir -p "$d"
  for r in "${RELEASES[@]}"; do
    bb "$tree" cli release verify --root "$root" --release "$r" > "$d/$tag-$r.log" 2>&1
    rc=$?
    log "verify $tag $r on $root: exit $rc: $(grep -E '^release ' "$d/$tag-$r.log" | tail -n 1)"
    [ "$rc" = 0 ] || bad=1
  done
  return $bad
}

allowance() { bb "$TREE" py -m biobox.runstep --profile "$PROFILE" --threads "$T" 2>> "$LOG"; }

reads() {  # $1 directory: two lookups by accession on the published releases over S1, P1, P2 (card W-36's controls)
  local d=$1 r acc
  for r in rM.2 rN.3; do
    for acc in NC_001604 NZ_AP026781.1; do
      bb "$TREE" cli find --root "$ROOT" --release "$r" --accession "$acc" > "$d/find-$r-$acc.out" 2> "$d/find-$r-$acc.err"
      log "find $r $acc: exit $?: $(wc -l < "$d/find-$r-$acc.out") lines"
    done
  done
}

step_argv() {  # $1 root, $2 stage, $3 scratch base
  echo stage labels --root "$1" --stage "$2" --cpus "$T" --tmp "$3/$2.labels"
}

launch() {  # $1 unit, $2 root, $3 stage, $4 profile, $5 scratch base
  local unit=$1 root=$2 s=$3 prof=$4 sb=$5
  mkdir -p "$sb/$s.labels"
  log "launch $unit: load $(cut -d' ' -f1-3 /proc/loadavg); account CPU $(top -b -n 1 -u if2453 | awk 'NR > 7 {s += $9} END {print s}')%; df /tmp $(df -h /tmp | awk 'NR == 2 {print $4}') free"
  log "launch $unit: run-step --root $root --stage $s --step labels --threads $T -- python -m biobox.cli $(step_argv "$root" "$s" "$sb")"
  # shellcheck disable=SC2046
  bash "$TREE/scripts/pod_run.sh" --profile "$prof" --image "$IMAGE" --digest "$DIGEST" --mem 8G --threads "$T" \
    --bind /tmp/if2453 --env PYTHONPATH="$TREE/src" --log "$LOGS/$unit.log" "$unit" run-step --profile "$prof" \
    --root "$root" --stage "$s" --step labels --threads "$T" -- python -m biobox.cli \
    $(step_argv "$root" "$s" "$sb") >> "$LOG" 2>&1
}

exit_of() { sed -n "s/^pod_run: $1 exited \([0-9]*\)$/\1/p" "$LOGS/$1.log" | tail -n 1; }

unit_cpu_min() {
  local ns
  ns=$(systemctl --user show -p CPUUsageNSec --value "$1" 2> /dev/null)
  [[ "$ns" =~ ^[0-9]+$ ]] || ns=0
  echo $((ns / 60000000000))
}

wait_one() {  # $1 unit: wait until it ends, guarding its CPU
  local u=$1 m
  sleep 5
  while systemctl --user is-active --quiet "$u"; do
    m=$(unit_cpu_min "$u")
    if [ "$m" -gt "$CPU_MAX" ]; then
      log "guard: $u at $m CPU-minutes, over $CPU_MAX; stopping it"
      systemctl --user stop "$u"
      fail "$u passed $CPU_MAX CPU-minutes and was stopped"
    fi
    sleep 15
  done
  log "$u ended: $(tail -n 2 "$LOGS/$u.log" | tr '\n' ' ')"
}

rehearse() {
  local RH=$C/rehearsal RPROF=$C/rehearsal-profile.yaml u rc s
  [ -e "$RH/OK" ] && { log "rehearsal: OK present; nothing to do"; return; }
  if [ ! -d "$RH" ]; then
    python3 -I "$BIN/c34_view.py" rehearsal "$ROOT" "$RH" "$(IFS=,; echo "${STAGES[*]}")" >> "$LOG" 2>&1 || fail "c34_view.py rehearsal"
  fi
  sed "s#^tmp:.*#tmp: $C/rehearsal-tmp#" "$PROFILE" > "$RPROF"
  for s in "${STAGES[@]}"; do
    [ -e "$RH/runs/stage=$s/labels.json" ] && continue
    { [ -e "$RH/STOP" ] || [ -e "$ROOT/STOP" ]; } && fail "STOP present before the rehearsal of $s"
    u=c34-rh-$(lower "$s")
    launch "$u" "$RH" "$s" "$RPROF" "$C/scratch-rh"
    wait_one "$u"
    rc=$(exit_of "$u")
    [ "$rc" = 0 ] || fail "rehearsal: $u exit '$rc' (read $LOGS/$u.log)"
  done
  if [ ! -d "$C/verifyview" ]; then
    python3 -I "$BIN/c34_view.py" verify "$ROOT" "$RH" "$C/verifyview" >> "$LOG" 2>&1 || fail "c34_view.py verify"
  fi
  verify_all "$OLD" "$C/verifyview" "$C/rehearsal-verify" old
  log "rehearsal: the code before the card on the verify view: $(grep -l 'not in manifest' "$C"/rehearsal-verify/old-*.log | wc -l) of ${#RELEASES[@]} releases report files not in their manifest"
  verify_all "$TREE" "$C/verifyview" "$C/rehearsal-verify" new || fail "rehearsal: a release does not verify on the view with the card's code"
  touch "$RH/OK"
  log "rehearsal: OK"
}

oracle() {
  local O=$C/oracle u=c34-oracle-s1 rc objs
  [ -e "$O/OK" ] && { log "oracle: OK present; nothing to do"; return; }
  [ -e "$C/rehearsal/OK" ] || fail "no $C/rehearsal/OK"
  objs=$(python3 -I -c 'import json,sys; a=json.load(open(sys.argv[1]))["argv"]; print(a[a.index("--objects") + 1])' "$ROOT/runs/stage=S1/ingest.json")
  log "oracle: S1's objects: $objs"
  mkdir -p "$O/scratch"
  if [ ! -e "$O/root/l0/l0_manifest/stage=S1" ]; then
    log "launch $u: stage ingest --root $O/root --stage S1 --objects $objs --tmp $O/scratch"
    bash "$TREE/scripts/pod_run.sh" --profile "$PROFILE" --image "$IMAGE" --digest "$DIGEST" --mem 16G --threads "$T" \
      --bind /tmp/if2453 --bind "$ROOT:$ROOT:ro" --env PYTHONPATH="$TREE/src" --log "$LOGS/$u.log" "$u" stage ingest \
      --root "$O/root" --stage S1 --objects "$objs" --tmp "$O/scratch" >> "$LOG" 2>&1
    wait_one "$u"
    rc=$(exit_of "$u")
    [ "$rc" = 0 ] || fail "oracle: $u exit '$rc' (read $LOGS/$u.log)"
  fi
  bb "$TREE" py "$BIN/c34_oracle.py" --stage S1 --step "$C/rehearsal" --ingest "$O/root" --root "$ROOT" > "$O/oracle.out" 2>&1
  rc=$?
  log "oracle: exit $rc: $(tail -n 1 "$O/oracle.out")"
  [ "$rc" = 0 ] || fail "oracle: the label files differ (read $O/oracle.out)"
  touch "$O/OK"
}

before() {
  local d=$C/before act
  [ -e "$d/DONE" ] && { log "before: DONE present; nothing to do"; return; }
  [ -e "$ROOT/STOP" ] && fail "$ROOT/STOP exists"
  [ -e "$C/rehearsal/OK" ] || fail "no $C/rehearsal/OK"
  [ -e "$C/oracle/OK" ] || fail "no $C/oracle/OK"
  act=$(systemctl --user list-units --type=service --type=scope --state=active --no-legend --plain | awk '{print $1}' | paste -sd' ' -)
  log "before: active units of the account: $act"
  if systemctl --user list-units --state=active --no-legend --plain 'r50-*' 'c24-*' 'c25-*' 'c26-*' 'b47-*' 'c33-*' | grep -q .; then
    fail "a build unit is active: $act"
  fi
  verify_all "$TREE" "$ROOT" "$d" new || fail "a release does not verify before the first write"
  bb "$TREE" cli releases --root "$ROOT" > "$d/releases.tsv" 2> "$d/releases.err"
  reads "$d"
  snapshot "$d"
  touch "$d/DONE"
}

run() {
  local s u rc a
  [ -e "$C/before/DONE" ] || fail "no $C/before/DONE: run the before mode first"
  [ -e "$C/run/DONE" ] && { log "run: DONE present; nothing to do"; return; }
  mkdir -p "$C/run"
  sed "s#^tmp:.*#tmp: $C/tmp-real#" "$PROFILE" > "$C/real-profile.yaml"
  for s in "${STAGES[@]}"; do
    if [ -e "$ROOT/runs/stage=$s/labels.json" ]; then
      log "$s: run record present; skipped"
      continue
    fi
    [ -e "$ROOT/runs/stage=$s/labels.publishing.json" ] && fail "$s: an intent file is present"
    [ -e "$ROOT/STOP" ] && fail "STOP present; no further launches: $(head -c 300 "$ROOT/STOP")"
    while :; do
      a=$(allowance)
      [[ "$a" =~ ^[0-9]+$ ]] && [ "$a" -ge "$T" ] && break
      log "allowance '$a' under $T threads; waiting a minute before $s"
      sleep 60
    done
    u=c34-$(lower "$s")
    log "allowance $a; launching $u"
    launch "$u" "$ROOT" "$s" "$C/real-profile.yaml" "$C/scratch"
    wait_one "$u"
    rc=$(exit_of "$u")
    log "$u: exit '$rc': $(grep -E 'published|failed|refused|STOP' "$LOGS/$u.log" | tail -n 1 | cut -c1-300); $(grep -E 'pod_run: .* usage' "$LOGS/$u.log" | tail -n 1)"
    if [ "$rc" != 0 ] || [ ! -e "$ROOT/runs/stage=$s/labels.json" ]; then
      fail "$u did not publish; no further launches (the root is left as the step left it)"
    fi
  done
  touch "$C/run/DONE"
}

after() {
  local d=$C/after
  mkdir -p "$d"
  snapshot "$d"
  verify_all "$TREE" "$ROOT" "$d" new || log "after: a release does NOT verify"
  bb "$TREE" cli releases --root "$ROOT" > "$d/releases.tsv" 2> "$d/releases.err"
  log "after: releases: $(wc -l < "$d/releases.tsv") lines"
  reads "$d"
  touch "$d/DONE"
}

check() {
  mkdir -p "$C/check"
  bb "$TREE" py "$BIN/c34_check.py" --root "$ROOT" --before "$C/before" --after "$C/after" --rehearsal "$C/rehearsal" \
    --out "$C/check" --stage S1 --stage P1 --stage P2 > "$C/check/check.out" 2>&1
  local rc=$?
  cp "$C/check/root_added.tsv" "$C/root_added.tsv"
  log "check: exit $rc: $(tail -n 1 "$C/check/check.out")"
  return $rc
}

seven() {
  local V=$C/view7
  if [ ! -d "$V" ]; then
    python3 -I "$BIN/c34_view.py" seven "$ROOT" "$V" >> "$LOG" 2>&1 || fail "c34_view.py seven"
    bb "$TREE" cli release define --root "$V" --release rC34.v --stages "$SEVEN" > "$C/seven-define.log" 2>&1 \
      || fail "release define in the view (read $C/seven-define.log)"
  fi
  local tag=${1:-now}
  bb "$TREE" cli releases --root "$V" > "$C/seven-releases-$tag.tsv" 2> "$C/seven-releases-$tag.err"
  log "seven ($tag): $(tail -n 1 "$C/seven-releases-$tag.tsv")"
}

MODE=${1:-}
log "c34 driver: mode '$MODE' (pid $$), tree $TREE, image $IMAGE"
case "$MODE" in
  rehearse) rehearse ;;
  oracle) oracle ;;
  before) before ;;
  run) run ;;
  after) after ;;
  check) check ;;
  seven) seven "${2:-now}" ;;
  *) echo "usage: c34_driver.sh rehearse|oracle|before|run|after|check|seven [tag]" >&2; exit 2 ;;
esac
log "c34 driver: mode '$MODE' done"
