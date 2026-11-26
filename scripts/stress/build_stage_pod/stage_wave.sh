#!/bin/bash
# Card R-07: a wave of stage builds, in the tree. The one place that says which stages a wave builds and in what order,
# so that a stopped wave can be read afterwards and resumed by a command that is in the tree.
#
#   bash stage_wave.sh [--resume-at STAGE --from-mode MODE] LOG STAGE...
#
# 1. Before anything runs, the arguments are checked: every STAGE is a stage name, none is listed twice, the list is not
#    empty; with --resume-at, the stage is the first of the list and --from-mode names one of estimate, rehearse,
#    before, run, after, check, backup. A refusal exits 2 and runs and logs nothing.
# 2. For each stage in the order given: the data root's STOP file is looked for first (exit 3); then
#    `stage_queue.sh STAGE` runs, and its exit status and the backup marker of the stage are read. With --resume-at,
#    the first stage is not run through the queue: `stage_driver.sh STAGE <mode>` runs for that mode and each mode
#    after it, each needing its own <mode>/DONE, with the STOP file looked for before each.
# 3. The wave ends at once, with a non-zero status and one log line that names the stage, the mode and the status, when
#    a stage or mode exits non-zero (that status is the wave's), when a mode that exited 0 left no DONE marker, or when
#    the STOP file exists (exit 3, a STOP is never lifted here). A stage that ended without its backup marker (exit 6)
#    is never counted as built. A failure never looks like progress: no later stage starts.
# Each start and end is one dated line in LOG, also printed; nothing else here carries a timestamp.
# 4. Card A-30 (DEC-81 (3)): the wave makes itself a new name, w<16 hex digits> from /dev/urandom (never one from the
#    environment), and exports it as STAGE_WAVE to every queue and driver it runs, so that the stages of this wave, and
#    only they, may share their rehearsal and the verification before a build (stage_driver.sh says on what). The name
#    and whether sharing is on are in the start line. STAGE_SHARE=0 in the wave's environment turns every share off;
#    without a name (urandom unreadable) nothing is shared.
#
# It launches nothing through systemd; the operator's line is in RUNBOOK.md, "A wave of stages". Variables, as in
# stage_queue.sh and stage_driver.sh: STAGE_BIN (the directory of the two scripts; default: this script's), STAGES_HOME
# (the stages' work directory; default /tmp/if2453/stages), STAGE_ROOT (the data root, whose STOP file is read;
# default /tmp/if2453/r10/root). The stage scripts keep every guard of theirs; the wave adds none.
set -u

BIN=${STAGE_BIN:-$(cd "$(dirname "$0")" && pwd)}
STAGES_HOME=${STAGES_HOME:-/tmp/if2453/stages}
DATAROOT=${STAGE_ROOT:-/tmp/if2453/r10/root}
MODES="estimate rehearse before run after check backup"
LOG=""

usage() {
  echo "usage: stage_wave.sh [--resume-at STAGE --from-mode MODE] LOG STAGE..." >&2
  exit 2
}
refuse() {
  echo "stage_wave: REFUSED: $*; nothing was run" >&2
  exit 2
}
say() {
  line="$(date '+%F %T %Z') stage_wave: $*"
  echo "$line" >> "$LOG"
  echo "$line"
}
lower() { echo "$1" | tr 'A-Z' 'a-z'; }

RESUME=""
FROM=""
while [ $# -gt 0 ]; do
  case "$1" in
    --resume-at) [ $# -ge 2 ] || usage; RESUME=$2; shift 2 ;;
    --from-mode) [ $# -ge 2 ] || usage; FROM=$2; shift 2 ;;
    --) shift; break ;;
    -*) usage ;;
    *) break ;;
  esac
done
[ $# -ge 1 ] || usage
LOG=$1
shift
[ $# -ge 1 ] || refuse "no stage listed"
[ -n "$RESUME" ] && [ -z "$FROM" ] && refuse "--resume-at needs --from-mode"
[ -z "$RESUME" ] && [ -n "$FROM" ] && refuse "--from-mode needs --resume-at"

seen=" "
for s in "$@"; do
  echo "$s" | grep -Eq '^[A-Z][A-Za-z0-9]{0,15}$' || refuse "'$s' is not a stage name"
  case "$seen" in *" $s "*) refuse "stage $s is listed twice" ;; esac
  seen="$seen$s "
done
if [ -n "$RESUME" ]; then
  [ "$RESUME" = "$1" ] || refuse "--resume-at $RESUME must be the first stage of the list (the list starts with $1)"
  case " $MODES " in *" $FROM "*) ;; *) refuse "--from-mode $FROM is not one of: $MODES" ;; esac
fi
[ -f "$BIN/stage_queue.sh" ] && [ -f "$BIN/stage_driver.sh" ] || refuse "no stage_queue.sh and stage_driver.sh in $BIN"
mkdir -p "$(dirname "$LOG")" || refuse "cannot make the directory of $LOG"

stopped_at() {  # STAGE MODE: the data root's STOP ends the wave
  if [ -e "$DATAROOT/STOP" ]; then
    say "STOPPED at stage $1 mode $2: $DATAROOT/STOP exists; nothing more is started (exit 3)"
    exit 3
  fi
}

WAVE_ID=$(od -An -N8 -tx1 /dev/urandom 2> /dev/null | tr -d ' \n')
if echo "$WAVE_ID" | grep -Eq '^[0-9a-f]{16}$'; then WAVE_ID=w$WAVE_ID; else WAVE_ID=""; fi
STAGE_WAVE=$WAVE_ID
export STAGE_WAVE
case "${STAGE_SHARE-1}" in
  1) SHARING="sharing on" ;;
  *) SHARING="sharing off (STAGE_SHARE=${STAGE_SHARE-})" ;;
esac
[ -n "$WAVE_ID" ] || SHARING="sharing off (no wave name could be made from /dev/urandom)"

say "start (pid $$): stages $*; resume ${RESUME:-none}${FROM:+ from mode $FROM}; scripts $BIN; stages in $STAGES_HOME; STOP at $DATAROOT/STOP; wave ${WAVE_ID:-none}; $SHARING"

first=1
for s in "$@"; do
  c=$STAGES_HOME/$(lower "$s")
  if [ "$first" = 1 ] && [ -n "$RESUME" ]; then
    first=0
    started=0
    for m in $MODES; do
      [ "$m" = "$FROM" ] && started=1
      [ "$started" = 1 ] || continue
      stopped_at "$s" "$m"
      say "stage $s mode $m: start (resume)"
      bash "$BIN/stage_driver.sh" "$s" "$m"
      rc=$?
      if [ "$rc" != 0 ]; then
        say "STOPPED at stage $s mode $m: exit $rc; read $c/driver.log and $c/$m/"
        exit "$rc"
      fi
      if [ ! -e "$c/$m/DONE" ]; then
        say "STOPPED at stage $s mode $m: exit 0 but no $c/$m/DONE; read $c/driver.log (exit 6)"
        exit 6
      fi
      say "stage $s mode $m: end, exit 0, marker present"
    done
    say "stage $s: end (resumed from $FROM), passed"
    continue
  fi
  first=0
  stopped_at "$s" queue
  say "stage $s: start"
  bash "$BIN/stage_queue.sh" "$s"
  rc=$?
  if [ "$rc" != 0 ]; then
    say "STOPPED at stage $s mode queue: exit $rc; read $c/queue.out and $c/driver.log"
    exit "$rc"
  fi
  if [ ! -e "$c/backup/DONE" ]; then
    say "STOPPED at stage $s mode backup: the queue exited 0 but left no $c/backup/DONE; read $c/queue.out (exit 6)"
    exit 6
  fi
  say "stage $s: end, exit 0, backup marker present"
done
say "done: every stage of the wave passed ($*)"
