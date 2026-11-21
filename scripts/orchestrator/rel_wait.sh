#!/bin/bash
# Orchestrator, 2026-11-20: cut the next rolling release (card C-46's release_train.sh) by itself once a stage's
# queued build (unit c47-queue: rebefore, run, after, check, backup) has ended with every mode passed, so that the
# night is used without the laptop.
#
#   rel_wait.sh [DRIVER_DIR [SUITE]]
#     DRIVER_DIR  the stage driver's directory (default /tmp/if2453/c47, stage C5); its queue unit is
#                 <basename>-queue (override: QUEUE_UNIT)
#     SUITE       the frozen tree's name (default rel1120): the tree /tmp/if2453/ci/SUITE, the shards' logs in
#                 /tmp/if2453/ci-tmp/SUITE-detached, the units ci-SUITE-*
#   REL_ROOT, REL_W, REL_TREE, REL_SUITE_DIR override the data root, the waiter's directory, the tree and the logs.
#
# It starts nothing unless ALL of these hold:
#   - no STOP file in the data root;
#   - the queue unit has ended and EVERY mode's success marker is there: before/DONE, run/DONE, after/DONE,
#     check/DONE and backup/DONE in DRIVER_DIR, and DRIVER_DIR/backup.tsv. The queue's log line is never the
#     evidence (card A-11: on 2026-11-20 a driver exited 0 after its check failed, the queue wrote "all passed", and
#     this waiter started a release on that line); the line is only printed when the markers are not all there;
#   - no build unit of ours is active (c47-*, c48c-*, c46-train-*, <basename of DRIVER_DIR>-*, and <stage>-<n> of every
#     stages/<STAGE>.conf of the tree: card A-12's stage driver), looked for after the queue, again just before the
#     dry run and again just before the run;
#   - the whole test suite of the frozen tree (units ci-SUITE-0..7) has ended and every shard's last line says
#     "passed" and not "failed" or "error".
# Then: release_train.sh --dry-run (logged), then release_train.sh run. The train has its own guards (STOP, free
# space, 280 CPU-minutes per process, DEC-43/DEC-53 rate guard, release verify of the new and every earlier release,
# pack verify, the worked examples). It never deletes or moves anything of an earlier release. Nothing is announced:
# the orchestrator audits first.
set -u
DRIVER_DIR=${1:-/tmp/if2453/c47}
SUITE_NAME=${2:-rel1120}
W=${REL_W:-/tmp/if2453/c46/wait}
ROOT=${REL_ROOT:-/tmp/if2453/r10/root}
TREE=${REL_TREE:-/tmp/if2453/ci/$SUITE_NAME}
SUITE=${REL_SUITE_DIR:-/tmp/if2453/ci-tmp/$SUITE_NAME-detached}
Q4=$DRIVER_DIR/queue.out
QUEUE_UNIT=${QUEUE_UNIT:-$(basename "$DRIVER_DIR")-queue}
MARKERS="before/DONE run/DONE after/DONE check/DONE backup/DONE backup.tsv"
say() { echo "$(date '+%F %T %Z') rel_wait: $*"; }

say "start (pid $$); waiting for $QUEUE_UNIT to end with every marker of $DRIVER_DIR"
t0=$(date +%s)
k=0
while :; do
  if [ -e "$ROOT/STOP" ]; then say "STOPPED: $ROOT/STOP exists; nothing started"; exit 3; fi
  if ! systemctl --user is-active --quiet "$QUEUE_UNIT"; then
    missing=""
    for m in $MARKERS; do
      case "$m" in
        *.tsv) [ -s "$DRIVER_DIR/$m" ] || missing="$missing $m" ;;
        *) [ -e "$DRIVER_DIR/$m" ] || missing="$missing $m" ;;
      esac
    done
    if [ -z "$missing" ]; then
      say "$QUEUE_UNIT ended and every marker is there ($MARKERS)"
      break
    fi
    say "STOPPED: $QUEUE_UNIT ended without the marker(s)$missing in $DRIVER_DIR; nothing started. Its last line: $(tail -n 1 "$Q4" 2> /dev/null | cut -c1-300)"
    exit 4
  fi
  if [ $(( $(date +%s) - t0 )) -gt 172800 ]; then say "STOPPED: waited over 48 h; nothing started"; exit 5; fi
  [ $((k % 30)) = 0 ] && say "waiting: $QUEUE_UNIT active; last: $(tail -n 1 "$Q4" 2> /dev/null | cut -c1-160)"
  k=$((k + 1))
  sleep 60
done

# the build units of ours: the older names, the stage directory's own (<basename>-*), and every stage of the tree's
# stage driver (stages/<STAGE>.conf: units <stage in lower case>-<n>, card A-12); looked for again just before the
# dry run and just before the run, since the suite's wait can last hours
build_units() {
  local pats=('c47-*' 'c48c-*' 'c46-train-*' "$(basename "$DRIVER_DIR")-*") f s
  for f in "$TREE"/scripts/stress/build_stage_pod/stages/*.conf; do
    [ -e "$f" ] || continue
    s=$(basename "$f" .conf)
    [[ "$s" =~ ^[A-Z][A-Za-z0-9]{0,15}$ ]] && pats+=("${s,,}-[0-9]*")
  done
  systemctl --user list-units --type=service --state=active,activating,deactivating --no-legend --plain "${pats[@]}" \
    | awk '{print $1}' | paste -sd' ' -
}
no_build_units() {
  act=$(build_units)
  if [ -n "$act" ]; then say "STOPPED: a build unit is active ($act); nothing started"; exit 6; fi
}
no_build_units

# the suite of the frozen tree must have ended, every shard passed
n=0
while systemctl --user list-units --type=service --state=active,activating --no-legend --plain "ci-$SUITE_NAME-*" | grep -q .; do
  n=$((n + 1))
  if [ "$n" -gt 360 ]; then say "STOPPED: the suite was still running 6 h after C5 ended; nothing started"; exit 7; fi
  [ $((n % 15)) = 1 ] && say "waiting for the suite (ci-$SUITE_NAME-*) to end"
  sleep 60
done
for s in 0 1 2 3 4 5 6 7; do
  last=$(tail -n 1 "$SUITE/log.$s" 2>/dev/null)
  case "$last" in
    *failed*|*error*) say "STOPPED: suite shard $s: $last; nothing started"; exit 7 ;;
    *passed*) say "suite shard $s: $last" ;;
    *) say "STOPPED: suite shard $s has no summary line ('$(echo "$last" | cut -c1-120)'); nothing started"; exit 7 ;;
  esac
done
if [ -e "$ROOT/STOP" ]; then say "STOPPED: $ROOT/STOP exists; nothing started"; exit 3; fi
no_build_units

say "dry run:"
ROOT=$ROOT TREE=$TREE bash "$TREE/scripts/stress/build_stage_pod/release_train.sh" --dry-run > "$W/dry_run.txt" 2>&1
rc=$?
cat "$W/dry_run.txt"
if [ "$rc" != 0 ]; then say "STOPPED: the dry run exited $rc; nothing started"; exit 8; fi
if ! grep -q '^Action: cut' "$W/dry_run.txt"; then say "STOPPED: the dry run does not say 'Action: cut'; nothing started"; exit 8; fi
if [ -e "$ROOT/STOP" ]; then say "STOPPED: $ROOT/STOP exists; nothing started"; exit 3; fi
no_build_units

say "release_train.sh run: start (Q=56: the unit's CPU quota; build-stage's own allowance holds the cap of the hour)"
Q=56 ROOT=$ROOT TREE=$TREE bash "$TREE/scripts/stress/build_stage_pod/release_train.sh" run > "$W/train.out" 2>&1
rc=$?
say "release_train.sh run: exit $rc; read /tmp/if2453/c46/train/driver.log and $W/train.out"
tail -n 5 "$W/train.out"
exit "$rc"
