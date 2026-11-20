#!/bin/bash
# Orchestrator, 2026-11-20: cut the next rolling release (card C-46's release_train.sh) by itself once stage C5's
# queued build (unit c47-queue: rebefore, run, after, check, backup) has ended with every mode passed, so that the
# night is used without the laptop. It starts nothing unless ALL of these hold:
#   - no STOP file in the data root;
#   - c47-queue has ended and its log has the line "done: rebefore, run, after, check and backup all passed", and
#     /tmp/if2453/c47/backup.tsv exists;
#   - no build unit of ours is active (c47-*, c48c-*, c46-train-*);
#   - the whole test suite of the frozen tree (units ci-rel1120-0..7) has ended and every shard's last line says
#     "passed" and not "failed" or "error".
# Then: release_train.sh --dry-run (logged), then release_train.sh run. The train has its own guards (STOP, free
# space, 280 CPU-minutes per process, DEC-43/DEC-53 rate guard, release verify of the new and every earlier release,
# pack verify, the worked examples). It never deletes or moves anything of an earlier release. Nothing is announced:
# the orchestrator audits first.
set -u
W=/tmp/if2453/c46/wait
ROOT=/tmp/if2453/r10/root
TREE=/tmp/if2453/ci/rel1120
SUITE=/tmp/if2453/ci-tmp/rel1120-detached
Q4=/tmp/if2453/c47/queue.out
say() { echo "$(date '+%F %T %Z') rel_wait: $*"; }

say "start (pid $$); waiting for c47-queue to end with its done line"
t0=$(date +%s)
k=0
while :; do
  if [ -e "$ROOT/STOP" ]; then say "STOPPED: $ROOT/STOP exists; nothing started"; exit 3; fi
  if ! systemctl --user is-active --quiet c47-queue; then
    if grep -q 'c47_queue: done: rebefore, run, after, check and backup all passed' "$Q4" && [ -s /tmp/if2453/c47/backup.tsv ]; then
      say "c47-queue ended with every mode passed"
      break
    fi
    say "STOPPED: c47-queue ended without its done line; nothing started. Its last line: $(tail -n 1 "$Q4" | cut -c1-300)"
    exit 4
  fi
  if [ $(( $(date +%s) - t0 )) -gt 172800 ]; then say "STOPPED: waited over 48 h; nothing started"; exit 5; fi
  [ $((k % 30)) = 0 ] && say "waiting: c47-queue active; last: $(tail -n 1 "$Q4" | cut -c1-160)"
  k=$((k + 1))
  sleep 60
done

act=$(systemctl --user list-units --type=service --state=active,activating,deactivating --no-legend --plain \
  'c47-*' 'c48c-*' 'c46-train-*' | awk '{print $1}' | paste -sd' ' -)
if [ -n "$act" ]; then say "STOPPED: a build unit is active ($act); nothing started"; exit 6; fi

# the suite of the frozen tree must have ended, every shard passed
n=0
while systemctl --user list-units --type=service --state=active,activating --no-legend --plain 'ci-rel1120-*' | grep -q .; do
  n=$((n + 1))
  if [ "$n" -gt 360 ]; then say "STOPPED: the suite was still running 6 h after C5 ended; nothing started"; exit 7; fi
  [ $((n % 15)) = 1 ] && say "waiting for the suite (ci-rel1120-*) to end"
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

say "dry run:"
ROOT=$ROOT TREE=$TREE bash "$TREE/scripts/stress/build_stage_pod/release_train.sh" --dry-run > "$W/dry_run.txt" 2>&1
rc=$?
cat "$W/dry_run.txt"
if [ "$rc" != 0 ]; then say "STOPPED: the dry run exited $rc; nothing started"; exit 8; fi
if ! grep -q '^Action: cut' "$W/dry_run.txt"; then say "STOPPED: the dry run does not say 'Action: cut'; nothing started"; exit 8; fi

say "release_train.sh run: start (Q=56: the unit's CPU quota; build-stage's own allowance holds the cap of the hour)"
Q=56 ROOT=$ROOT TREE=$TREE bash "$TREE/scripts/stress/build_stage_pod/release_train.sh" run > "$W/train.out" 2>&1
rc=$?
say "release_train.sh run: exit $rc; read /tmp/if2453/c46/train/driver.log and $W/train.out"
tail -n 5 "$W/train.out"
exit "$rc"
