#!/bin/bash
# Card C-33 part 1 (DEC-40): the two CRISPR array producers on the seven built stages of the data root, each stage
# and producer one `run-step` process through scripts/pod_run.sh, with the runtime image of `faec1f4` and the tree's
# src. Host side, bash. Copies of this script, c33_view.py and c33_check.py run from /tmp/if2453/c33/bin, so the
# synced tree is not touched while they run. Each mode is launched in its own small unit, for example
#
#   systemd-run --user --unit=c33-driver-run --collect -p MemoryMax=4G -p MemorySwapMax=0 \
#     -p StandardOutput=append:/tmp/if2453/c33/driver.out -p StandardError=append:/tmp/if2453/c33/driver.out \
#     nice -n 19 bash /tmp/if2453/c33/bin/c33_driver.sh run
#
# Modes, in the order the card runs them:
#   rehearse  A scratch root /tmp/if2453/c33/rehearsal (c33_view.py: symlinks to the sequences of every stage and to
#             S1's gene calls and Pfam scan; its own runs/, ledger/, l1/feature, l1/crispr), a scratch copy of the
#             profile with its own tmp, and the exact commands of the real run for S1 (default producer, then the
#             filtered one) with --root the scratch root; then the sha256 of the six new files against
#             /tmp/if2453/e15/view and of the default's three against /tmp/if2453/e12/view. Writes nothing in ROOT.
#   before    Refuses when ROOT/STOP exists or a build unit (r50-*, c24-*, c25-*, c26-*, b47-*) is active; then
#             `release verify` of rM.2 and rN.3 (root bound read-only) and the snapshot into C/before: listing.tsv
#             (path, bytes, mtime of every file of ROOT) and p1s1.sha256 (every file under a stage=P1 or stage=S1
#             directory: the tables' files, run records, ledgers). Writes C/before/DONE when both verify.
#   run       Needs C/before/DONE and C/rehearsal/OK. The default producer on S1, P1, ..., P6, then the filtered one on
#             each, JOBS units at once of T threads (T MinCED processes, `--cpus T`; run-step `--threads T`;
#             pod_run `--threads T --mem 16G`). Before each launch: ROOT/STOP (stop launching), the step's run
#             record (skip) or intent file (stop), and the account's allowance for a T-thread launch (`python -m
#             biobox.runstep`, card C-28; wait while it is under T). While units run, once a minute: a unit of
#             this driver over CPU_MAX CPU-minutes (systemd's CPUUsageNSec) is stopped and the driver fails. A unit
#             that does not end `exited 0` with its run record published stops new launches; the running ones go on
#             to their end; nothing in the root is cleaned. Writes C/run/DONE or C/run/FAILED.
#   after     The same snapshot into C/after, then `release verify` of rM.2 and rN.3 again.
#   check     c33_check.py in the image: before against after, the oracle against the two views on S1, P1, P2, the
#             per-stage counts and costs, and C/check/root_added.tsv (copied to C/root_added.tsv).
#
# The step's argv is the scratch views' (cards E-12, E-15: `stage crispr [--near-cas --library pfam
# --library-version 38.2] --root R --stage S --cpus N --tmp D --seq-stages ...`), with --seq-stages as build-stage
# gives it (`with_seq_stages` in src/biobox/stage.py): the stage itself, then every other stage whose seq_locator is in
# the root, sorted by name. A contig another stage owns is skipped and named in the ledger (DEC-26, DEC-27).
set -u

C=/tmp/if2453/c33
ROOT=/tmp/if2453/r10/root
TREE=${TREE:-/tmp/if2453/ci/C-33}
PROFILE=$TREE/profiles/wilkcomp02.yaml
IMAGE=/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif
DIGEST=sha256:d67300aa1321b571fa8da79e74dfccdd87a8dca2e56518f324902efa6b9f098a
BIN=$(cd "$(dirname "$0")" && pwd)
LOGS=$C/logs
LOG=$C/driver.log
STAGES=(S1 P1 P2 P3 P4 P5 P6)
NR3=minced-0.4.2-w8-nr3-r23x47-s26x50
NEAR=minced-0.4.2-w8-nr2-r23x47-s26x50-cas1kb-caspfam1-pfam-38.2
T=${T:-4}
JOBS=${JOBS:-5}
CPU_MAX=${CPU_MAX:-280}
E12=/tmp/if2453/e12/view
E15=/tmp/if2453/e15/view
mkdir -p "$LOGS" "$C/scratch"

log() { echo "$(date '+%F %T %Z') $*" >> "$LOG"; }
fail() { log "FAILED: $*"; mkdir -p "$C/$MODE"; echo "$*" > "$C/$MODE/FAILED"; exit 1; }

bb() {  # a short command in the image, in its own capped scope at nice 19, the root bound read-only
  local what=$1
  shift
  case "$what" in
    cli) set -- python -m biobox.cli "$@" ;;
    py) set -- python "$@" ;;
  esac
  systemd-run --user --scope --quiet -p MemoryMax=16G -p MemorySwapMax=0 -p CPUQuota=200% nice -n 19 \
    singularity exec --cleanenv --env PYTHONPATH="$TREE/src" --env OMP_NUM_THREADS=2 --env PYTHONUNBUFFERED=1 \
    --bind /tmp/if2453 --bind "$ROOT:$ROOT:ro" --bind /sys/fs/cgroup:/sys/fs/cgroup:ro "$IMAGE" "$@"
}

seq_stages() {  # $1 root, $2 stage: as build-stage's with_seq_stages
  local s out=$2
  for s in $(ls "$1/l0seq/seq_locator" | sed -n 's/^stage=//p' | LC_ALL=C sort); do
    [ "$s" = "$2" ] || out=$out,$s
  done
  echo "$out"
}

lower() { echo "$1" | tr 'A-Z' 'a-z'; }

snapshot() {  # $1 directory
  local d=$1
  mkdir -p "$d"
  find "$ROOT" ! -type d -printf '%P\t%s\t%T@\n' | LC_ALL=C sort > "$d/listing.tsv"
  find "$ROOT" ! -type d \( -path '*/stage=P1/*' -o -path '*/stage=S1/*' \) -printf '%P\n' | LC_ALL=C sort > "$d/p1s1.files"
  (cd "$ROOT" && systemd-run --user --scope --quiet -p MemoryMax=2G -p MemorySwapMax=0 nice -n 19 \
    xargs -d '\n' sha256sum < "$d/p1s1.files") > "$d/p1s1.sha256" || fail "snapshot $d: sha256sum failed"
  [ "$(wc -l < "$d/p1s1.sha256")" = "$(wc -l < "$d/p1s1.files")" ] || fail "snapshot $d: not every file hashed"
  log "snapshot $d: $(wc -l < "$d/listing.tsv") files, $(wc -l < "$d/p1s1.sha256") of P1 and S1 hashed, $(awk -F '\t' '{s += $2} END {printf "%.0f", s}' "$d/listing.tsv") bytes"
}

verify() {  # $1 release, $2 directory
  bb cli release verify --root "$ROOT" --release "$1" > "$2/verify-$1.log" 2>&1
  local rc=$?
  log "verify $1: exit $rc: $(grep -E '^release ' "$2/verify-$1.log" | tail -n 1)"
  return $rc
}

allowance() { bb py -m biobox.runstep --profile "$PROFILE" --threads "$T" 2>> "$LOG"; }

step_argv() {  # $1 kind, $2 root, $3 stage, $4 scratch base: the bio-box arguments of the step command
  local extra=()
  [ "$1" = near ] && extra=(--near-cas --library pfam --library-version 38.2)
  echo stage crispr "${extra[@]}" --root "$2" --stage "$3" --cpus "$T" --tmp "$4/$3.$(step_of "$1")" \
    --seq-stages "$(seq_stages "$2" "$3")"
}

step_of() { if [ "$1" = near ]; then echo "crispr-$NEAR"; else echo "crispr-$NR3"; fi; }

launch() {  # $1 unit, $2 kind, $3 root, $4 stage, $5 profile, $6 scratch base
  local unit=$1 kind=$2 root=$3 s=$4 prof=$5 sb=$6 step
  step=$(step_of "$kind")
  mkdir -p "$sb/$s.$step"
  log "launch $unit: load $(cut -d' ' -f1-3 /proc/loadavg); account CPU $(top -b -n 1 -u if2453 | awk 'NR > 7 {s += $9} END {print s}')%; df /tmp $(df -h /tmp | awk 'NR == 2 {print $4}') free"
  log "launch $unit: run-step --root $root --stage $s --step $step --threads $T -- python -m biobox.cli $(step_argv "$kind" "$root" "$s" "$sb")"
  # shellcheck disable=SC2046
  bash "$TREE/scripts/pod_run.sh" --profile "$prof" --image "$IMAGE" --digest "$DIGEST" --mem 16G --threads "$T" \
    --bind /tmp/if2453 --env PYTHONPATH="$TREE/src" --log "$LOGS/$unit.log" "$unit" run-step --profile "$prof" \
    --root "$root" --stage "$s" --step "$step" --threads "$T" -- python -m biobox.cli \
    $(step_argv "$kind" "$root" "$s" "$sb") >> "$LOG" 2>&1
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
    sleep 30
  done
  log "$u ended: $(tail -n 2 "$LOGS/$u.log" | tr '\n' ' ')"
}

shas() {  # $1 dir, $2 stage, $3 producer: sha256 of the producer's files of the stage, relative paths
  (cd "$1" && find "l1/feature/producer=$3/stage=$2" "l1/crispr/producer=$3/stage=$2" \
    "ledger/stage=$2/step=crispr-$3" -name '*.parquet' | LC_ALL=C sort | xargs -d '\n' sha256sum)
}

rehearse() {
  local RH=$C/rehearsal RPROF=$C/rehearsal-profile.yaml kind u rc ok=1 view p
  [ -e "$RH/OK" ] && { log "rehearsal: OK present; nothing to do"; return; }
  if [ ! -d "$RH" ]; then
    python3 -I "$BIN/c33_view.py" "$ROOT" "$RH" S1 "$(IFS=,; echo "${STAGES[*]}")" >> "$LOG" 2>&1 || fail "c33_view.py"
  fi
  sed "s#^tmp:.*#tmp: $C/rehearsal-tmp#" "$PROFILE" > "$RPROF"
  for kind in nr3 near; do
    [ -e "$RH/runs/stage=S1/$(step_of "$kind").json" ] && continue
    [ -e "$RH/STOP" ] || [ -e "$ROOT/STOP" ] && fail "STOP present before the rehearsal's $kind step"
    u=c33-rh-$kind-s1
    printf '%s\t%s\t%s\t%s\n' "$u" S1 "$kind" "$(step_of "$kind")" >> "$C/units-rehearsal.tsv"
    launch "$u" "$kind" "$RH" S1 "$RPROF" "$C/scratch-rh"
    wait_one "$u"
    rc=$(exit_of "$u")
    [ "$rc" = 0 ] || fail "rehearsal: $u exit '$rc' (read $LOGS/$u.log)"
  done
  mkdir -p "$C/rehearsal-cmp"
  for kind in nr3 near; do
    if [ "$kind" = near ]; then p=$NEAR; else p=$NR3; fi
    shas "$RH" S1 "$p" > "$C/rehearsal-cmp/rehearsal-$kind.sha256"
    for view in "$E15" "$E12"; do
      [ "$kind" = near ] && [ "$view" = "$E12" ] && continue
      shas "$view" S1 "$p" > "$C/rehearsal-cmp/$(basename "$(dirname "$view")")-$kind.sha256"
      if cmp -s "$C/rehearsal-cmp/rehearsal-$kind.sha256" "$C/rehearsal-cmp/$(basename "$(dirname "$view")")-$kind.sha256"; then
        log "rehearsal: S1 $kind: $(wc -l < "$C/rehearsal-cmp/rehearsal-$kind.sha256") files byte-identical to $view"
      else
        log "rehearsal: S1 $kind: DIFFERENT from $view"
        ok=0
      fi
    done
  done
  [ "$ok" = 1 ] || fail "rehearsal: not byte-identical to the views (read $C/rehearsal-cmp)"
  touch "$RH/OK"
  log "rehearsal: OK"
}

before() {
  local d=$C/before act
  [ -e "$d/DONE" ] && { log "before: DONE present; nothing to do"; return; }
  [ -e "$ROOT/STOP" ] && fail "$ROOT/STOP exists"
  act=$(systemctl --user list-units --type=service --type=scope --state=active --no-legend --plain | awk '{print $1}' | paste -sd' ' -)
  log "before: active units of the account: $act"
  if systemctl --user list-units --state=active --no-legend --plain 'r50-*' 'c24-*' 'c25-*' 'c26-*' 'b47-*' | grep -q .; then
    fail "a build unit is active: $act"
  fi
  mkdir -p "$d"
  verify rM.2 "$d" || fail "rM.2 does not verify before the first write"
  verify rN.3 "$d" || fail "rN.3 does not verify before the first write"
  snapshot "$d"
  touch "$d/DONE"
}

after() {
  local d=$C/after ok=1
  mkdir -p "$d"
  rm -f "$d/DONE"  # withdrawn until this `after` has passed
  snapshot "$d"
  verify rM.2 "$d" || { log "after: rM.2 does NOT verify"; ok=0; }
  verify rN.3 "$d" || { log "after: rN.3 does NOT verify"; ok=0; }
  [ "$ok" = 1 ] || fail "after: a published release failed verify (read $d/verify-*.log)"
  touch "$d/DONE"
}

run_phase() {  # $1 kind: every stage, JOBS at once
  local kind=$1 queue=("${STAGES[@]}") running=() still u s rc a n bad=0 step
  step=$(step_of "$kind")
  while [ ${#queue[@]} -gt 0 ] || [ ${#running[@]} -gt 0 ]; do
    still=()
    for u in "${running[@]}"; do
      if systemctl --user is-active --quiet "$u"; then
        if [ "$(unit_cpu_min "$u")" -gt "$CPU_MAX" ]; then
          log "guard: $u over $CPU_MAX CPU-minutes; stopping it"
          systemctl --user stop "$u"
          bad=1
        fi
        still+=("$u")
        continue
      fi
      rc=$(exit_of "$u")
      s=$(awk -F '\t' -v u="$u" '$1 == u {print $2}' "$C/units.tsv" | tail -n 1)
      log "$u ended, exit '$rc': $(grep -E 'published|failed|refused|STOP' "$LOGS/$u.log" | tail -n 1 | cut -c1-300); $(grep -E 'pod_run: .* usage' "$LOGS/$u.log" | tail -n 1)"
      if [ "$rc" != 0 ] || [ ! -e "$ROOT/runs/stage=$s/$step.json" ]; then
        log "$u did not publish; no further launches (the root is left as the step left it)"
        bad=1
      fi
    done
    running=("${still[@]+"${still[@]}"}")
    if [ "$bad" = 0 ] && [ ${#queue[@]} -gt 0 ] && [ ${#running[@]} -lt "$JOBS" ]; then
      s=${queue[0]}
      if [ -e "$ROOT/runs/stage=$s/$step.json" ]; then
        log "$s $kind: run record present; skipped"
        queue=("${queue[@]:1}")
        continue
      fi
      [ -e "$ROOT/runs/stage=$s/$step.publishing.json" ] && { log "$s $kind: an intent file is present"; bad=1; continue; }
      if [ -e "$ROOT/STOP" ]; then
        log "STOP present; no further launches: $(head -c 300 "$ROOT/STOP")"
        bad=2
        continue
      fi
      a=$(allowance)
      if ! [[ "$a" =~ ^[0-9]+$ ]] || [ "$a" -lt "$T" ]; then
        log "allowance '$a' under $T threads; waiting a minute before $s $kind"
        sleep 60
        continue
      fi
      n=$(($( [ -e "$C/units.tsv" ] && wc -l < "$C/units.tsv" || echo 0) + 1))
      u=c33-$kind-$(lower "$s")-$n
      printf '%s\t%s\t%s\t%s\n' "$u" "$s" "$kind" "$step" >> "$C/units.tsv"
      log "allowance $a; launching $u"
      launch "$u" "$kind" "$ROOT" "$s" "$PROFILE" "$C/scratch"
      running+=("$u")
      queue=("${queue[@]:1}")
      sleep 10
      continue
    fi
    [ "$bad" != 0 ] && [ ${#running[@]} -eq 0 ] && break
    sleep 30
  done
  return $bad
}

run() {
  local rc
  [ -e "$C/before/DONE" ] || fail "no $C/before/DONE: run the before mode first"
  [ -e "$C/rehearsal/OK" ] || fail "no $C/rehearsal/OK: the rehearsal has not passed"
  [ -e "$C/run/DONE" ] && { log "run: DONE present; nothing to do"; return; }
  mkdir -p "$C/run"
  run_phase nr3
  rc=$?
  [ "$rc" = 0 ] || fail "the default producer's phase ended with $rc (1 a failed or stopped unit, 2 STOP)"
  run_phase near
  rc=$?
  [ "$rc" = 0 ] || fail "the filtered producer's phase ended with $rc (1 a failed or stopped unit, 2 STOP)"
  touch "$C/run/DONE"
}

check() {
  [ -e "$C/after/DONE" ] || fail "check: no $C/after/DONE"
  mkdir -p "$C/check"
  rm -f "$C/check/DONE"  # withdrawn until this check has passed
  bb py "$BIN/c33_check.py" --root "$ROOT" --before "$C/before" --after "$C/after" --units "$C/units.tsv" \
    --logs "$LOGS" --out "$C/check" --view "e15=$E15" --view "e12=$E12" > "$C/check/check.out" 2>&1
  local rc=$?
  log "check: exit $rc: $(tail -n 1 "$C/check/check.out")"
  [ "$rc" = 0 ] || fail "check: c33_check.py exited $rc (read $C/check/check.out)"
  [ -s "$C/check/root_added.tsv" ] || fail "check: no $C/check/root_added.tsv"
  cp "$C/check/root_added.tsv" "$C/root_added.tsv" || fail "check: cp root_added.tsv"
  touch "$C/check/DONE"
}

MODE=${1:-}
log "c33 driver: mode '$MODE' (pid $$), tree $TREE, image $IMAGE"
case "$MODE" in
  rehearse) rehearse ;;
  before) before ;;
  run) run ;;
  after) after ;;
  check) check ;;
  *) echo "usage: c33_driver.sh rehearse|before|run|after|check" >&2; exit 2 ;;
esac
# The mode's own status is the driver's exit status (card A-11; on 2026-11-20 c47_driver.sh `check` returned 1 and the
# driver exited 0, the status of the log line that followed). A mode that fails through `fail` exits 1 itself.
rc=$?
if [ "$rc" = 0 ]; then
  log "c33 driver: mode '$MODE' done"
else
  log "c33 driver: mode '$MODE' FAILED with status $rc"
fi
exit "$rc"
