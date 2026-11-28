#!/bin/bash
# Card A-12: one stage driver for every stage, in place of a copy per stage (c47_driver.sh, c48c_v1_driver.sh, ...).
#
#   bash stage_driver.sh STAGE MODE        (on the pod host; each mode in a small detached unit, or stage_queue.sh)
#
# What differs between stages is in stages/<STAGE>.conf (read and checked by `stage_est.py conf`): the inputs
# directory and the sha256 of its files as the stage's note recorded them, the depth, a least contig length, the
# library or none, CRISPR or not, the stages it comes after, the figures its note says a step is judged against, the
# contigs earlier stages hold, and records a decision row accepts. Nothing else in this file names a stage, a release,
# an input sum or a thread count.
#
# Modes (each exits with its own status; each writes <mode>/DONE, holding the conf's sha256, only on success, after
# withdrawing its own marker and those of the modes after it; each refuses without the markers of the modes before it,
# and refuses a marker made for another conf):
#   estimate  The inputs against the conf's sums (sha256sum -c); `stage_est.py before`: the code's estimates from the
#             inputs (exit 5: a contig is owned by a stage parked half built; 7: build-stage would refuse; 8: an
#             EXPECT line differs from the contigs that stage holds).
#   rehearse  The fixture (tests/fixtures/micro_s0), in scratch under the stage's directory: stage E and a release
#             over it built in a scratch root, then a stage of the conf's depth, length, library and CRISPR choice
#             driven through this driver's own modes estimate, before, run, after, check and backup (into a scratch
#             backup directory), then a relaunch that must publish nothing. In a wave (below), the marker of another
#             stage's full rehearsal of the same wave is taken instead when every fact it rests on is equal.
#   before    No active build unit in the root; `release verify` of every release the root has a manifest for
#             (manifests/release=<R>.json, read from the root); the listing of the root without the stage's own
#             files, and the sha256 of every other stage's files. In a wave (below), the snapshot is taken first, and
#             the verification is taken from the newest full one of the wave when the root is provably unchanged.
#   run       The inputs' sums again; `build-stage` through pod_run.sh in the image by digest. The unit's CPUQuota is
#             the profile's cap of the hour (runstep.launch_cap) less the driver's helpers' 2, set at launch and again
#             with `systemctl set-property` whenever the hour's cap changes. STOP before the launch and while it runs;
#             the guards: r50_watch.py (a process at 280 CPU-minutes, /tmp under its floor, scratch over its ceiling)
#             and `stage_est.py rates` (DEC-43, DEC-53; any exit but 0 stops the unit, and a stop writes run/FAILED,
#             which only a person moves aside). The unit counts as ended only on its own exit line or a real
#             inactive answer of systemd; after WATCH_LIMIT_S (900) without one readable answer the driver fails,
#             saying the unit's state is unknown and that it may still be running. A trap stops the unit if the
#             driver exits while it runs; the helpers run in the background so that a TERM reaches the trap at once
#             (a SIGKILL of the driver leaves the unit running unguarded: RUNBOOK "Building a stage"). At the end,
#             the unit's exit 0 and the stage complete at its depth (`stage.stage_incomplete`).
#   after     The listing and sha256 again, the same releases verified with the same manifest hashes, `bio-box
#             releases`.
#   check     `stage_check.py` on the two snapshots, its report in check/<n>/: root_added.tsv (every added file with
#             size and sha256).
#   backup    One tar of exactly the added files, written to $BACKUP_DIR (or, when under $TAR_FLOOR_BYTES and
#             $SMALL_BACKUP_DIR names a non-/stor directory, to $SMALL_BACKUP_DIR if at least $SMALL_TAR_FLOOR_BYTES)
#             under a dot-temporary name with dd conv=fsync, read back, sha256 equal, mv -n; never over an existing
#             file; at least 64 MiB under /stor.
#   info      (for stage_queue.sh) the stages this one comes after and the active build units; writes no marker.
#
# Environment (defaults are the pod's): STAGES_HOME (/tmp/if2453/stages; the stage's directory is <home>/<stage in
# lower case>, override STAGE_DIR), STAGE_ROOT (the data root), STAGE_CONF (not in the data root), TREE, PROFILE,
# IMAGE, DIGEST, BACKUP_DIR (/stor/work/Wilke/if2453/all-the-data/root_backups; the rehearsal sets its own scratch
# one), SMALL_BACKUP_DIR (unset by default; optional non-/stor directory for tars under TAR_FLOOR_BYTES and at least
# SMALL_TAR_FLOOR_BYTES, default 10240, card R-25), T (6) and CPUS (4) (build-stage --threads and --scan-cpus, at most
# the profile's smaller cap), Q (lowers the unit's quota), MEM (the profile's mem_scope, at most 750G in the data root),
# CPU_MAX (280, never above), START_MIN_GB (205), TMP_MIN_GB (200) and SCRATCH_MAX_GB (100), which the data root never
# loosens.
# Unit names: <stage in lower case>-<n> for the builds, <stage>-rh-* for the rehearsal, <stage>-queue for the queue.
#
# A wave (card A-30, DEC-81 (3)): STAGE_WAVE, a name stage_wave.sh makes anew for each wave and passes to every stage it
# runs (never set it by hand). Without it, nothing below happens and every mode runs as it did before. With it, and
# unless STAGE_SHARE is anything but 1 (STAGE_SHARE=0 turns every share off), two checks may be taken from an earlier
# stage of the same wave, each only when every fact compared is equal, each logged with what was taken from where on
# which equality, and each run in full on any difference, missing marker or unreadable value:
#   rehearse  another stage's rehearse/DONE whose second line is "rehearsal=full" and whose facts (line 4 on) are
#             byte-equal to this stage's: the wave, the tree (realpath), the code's sha256 (code_sha), the image and
#             its digest, the profile's sha256, the conf's kind (DEPTH, MIN_CONTIG_LEN, GENECALL_MODE, LIBNAME, LIBVER,
#             LIBCUT, CRISPR) and the rehearsal's limits (RH_T, RH_Q, CPU_MAX, TMP_MIN_GB, START_MIN_GB, SCRATCH_MAX_GB,
#             MEM). The taking stage's DONE holds its own conf's sha256, "rehearsal=shared from stage <S> (<path>)" and
#             the same facts; a shared marker is never taken again.
#   before    the newest <stage>/{after,before}/verified.share under STAGES_HOME of this wave (written only after a full
#             verification): its mode's DONE made for the conf it names; the same root, tree, code sha256, image and
#             digest; its listing.all.tsv and earlier.sha256 unchanged since (their sha256 recorded); this mode's
#             listing.all.tsv (name, size, mtime of every file of the root) byte-equal to it; every line of its
#             earlier.sha256 (the sha256 of every other stage's file) in this mode's; no file of the root newer than
#             the record; the same releases, each manifest with the sha256 recorded. Each verify-<R>.log then says the
#             verification was taken, and from where, and carries the manifest hash that `after` compares.
#   after     always verifies every release in full, and records it for the next stage.
# Nothing else is shared: the inputs' sums, estimate, every guard of run, check, backup, and FAILED, REFUSED, STOPPED
# and STOP are as without a wave. A rehearsal's own modes never share.
set -u
set -o pipefail

STAGE=${1:-}
MODE=${2:-}
usage() { echo "usage: stage_driver.sh STAGE estimate|rehearse|before|run|after|check|backup|info" >&2; exit 2; }
[[ "$STAGE" =~ ^[A-Z][A-Za-z0-9]{0,15}$ ]] || usage
case "$MODE" in estimate|rehearse|before|run|after|check|backup|info) ;; *) usage ;; esac
LC=${STAGE,,}

BIN=$(cd "$(dirname "$0")" && pwd)
SELF=$BIN/$(basename "$0")
TREE=${TREE:-$(cd "$BIN/../../.." && pwd)}
REAL_ROOT=$(realpath -m /tmp/if2453/r10/root)
ROOT=$(realpath -m "${STAGE_ROOT:-$REAL_ROOT}")
DATAROOT=0
[ "$ROOT" = "$REAL_ROOT" ] && DATAROOT=1
STAGES_HOME=${STAGES_HOME:-/tmp/if2453/stages}
C=${STAGE_DIR:-$STAGES_HOME/$LC}
OWN_CONF=$BIN/stages/$STAGE.conf
CONF=${STAGE_CONF:-$OWN_CONF}
REALPROF=${PROFILE:-$TREE/profiles/wilkcomp02.yaml}
PROF=$C/profile.yaml
IMAGE=${IMAGE:-/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif}
DIGEST=${DIGEST:-sha256:d67300aa1321b571fa8da79e74dfccdd87a8dca2e56518f324902efa6b9f098a}
PYPATH=$TREE/src
PFAMDIR=/stor/scratch/Wilke/if2453/all-the-data/l0/refdb/pfam/38.2
PFAM=$PFAMDIR/Pfam-A.hmm.gz
STOR_RO=/stor/scratch/Wilke/if2453/all-the-data
BACKUP_DIR=${BACKUP_DIR:-/stor/work/Wilke/if2453/all-the-data/root_backups}
TAR_FLOOR_BYTES=${TAR_FLOOR_BYTES:-67108864}
SMALL_BACKUP_DIR=${SMALL_BACKUP_DIR:-}
SMALL_TAR_FLOOR_BYTES=${SMALL_TAR_FLOOR_BYTES:-10240}
STOR_FLOOR=67108864  # AGENTS §7.3: no file under 64 MB onto /stor
T=${T:-6}
CPUS=${CPUS:-4}
SCAN_JOBS=${SCAN_JOBS:-}
Q=${Q:-}
MEM=${MEM:-}
CPU_MAX=${CPU_MAX:-280}
TMP_MIN_GB=${TMP_MIN_GB:-200}
START_MIN_GB=${START_MIN_GB:-205}
SCRATCH_MAX_GB=${SCRATCH_MAX_GB:-100}
RH_T=${RH_T:-2}
RH_Q=${RH_Q:-4}
REHEARSAL=${STAGE_REHEARSAL:-0}
WAVE=${STAGE_WAVE:-}
SHARE=${STAGE_SHARE-1}
STEM=$LC
if [ "$REHEARSAL" = 1 ] && [ -n "${STAGE_UNIT_STEM:-}" ]; then STEM=$STAGE_UNIT_STEM; fi
[[ "$STEM" =~ ^[a-z0-9][a-z0-9-]{0,30}$ ]] || usage
LOGS=$C/logs
LOG=$C/driver.log
SCRATCH=$C/tmp/scratch
CHAIN="estimate before run after check backup"
RUNNING_UNIT=""
RUNNING_LOGS=""
CONF_SHA=""

mkdir -p "$C" "$LOGS" "$C/tmp" || exit 1

log() { echo "$(date '+%F %T.%3N %Z') [$STAGE $MODE] $*" >> "$LOG"; }
say() { echo "stage_driver $STAGE $MODE: $*" >&2; }
fail() {  # a failure after work began: <mode>/FAILED (run/FAILED blocks a new launch until a person moves it aside)
  log "FAILED: $*"
  say "FAILED: $*"
  [ "$MODE" = info ] && exit 1
  mkdir -p "$C/$MODE"
  echo "$*" > "$C/$MODE/FAILED"
  exit 1
}
refuse() {  # a refusal before any work: <mode>/REFUSED, which blocks nothing
  log "REFUSED: $*"
  say "REFUSED: $*"
  [ "$MODE" = info ] && exit 1
  mkdir -p "$C/$MODE"
  echo "$*" > "$C/$MODE/REFUSED"
  exit 1
}
stop_file() {  # a STOP in the root this mode works on, or in the data root: exit 3, no marker
  local r
  for r in "$ROOT" "$REAL_ROOT"; do
    if [ -e "$r/STOP" ]; then
      log "STOP: $r/STOP exists ($(head -c 200 "$r/STOP" 2> /dev/null)); nothing started"
      say "STOP: $r/STOP exists; nothing started"
      mkdir -p "$C/$MODE"
      echo "$r/STOP" > "$C/$MODE/STOPPED"
      exit 3
    fi
  done
}
withdraw() {  # this mode's marker and those of the modes after it, before any work
  local m seen=0
  for m in $CHAIN; do
    [ "$m" = "$1" ] && seen=1
    [ "$seen" = 1 ] && rm -f "$C/$m/DONE"
  done
  mkdir -p "$C/$1"
  rm -f "$C/$1/REFUSED" "$C/$1/STOPPED"
  [ "$1" = run ] || rm -f "$C/$1/FAILED"
}
mark_done() { echo "$CONF_SHA" > "$C/$1/DONE" || fail "could not write $C/$1/DONE"; }
need() {  # the marker of mode $1, made for this conf
  [ -e "$C/$1/DONE" ] || refuse "no $C/$1/DONE"
  [ "$(head -c 64 "$C/$1/DONE")" = "$CONF_SHA" ] \
    || refuse "$C/$1/DONE was made for another conf than $CONF (sha256 $CONF_SHA); run $1 again"
}

environment() {  # what the data root never allows
  local v n u
  [ "$DATAROOT" = 1 ] || return 0
  [ "$REHEARSAL" = 1 ] && refuse "a rehearsal never works in the data root $REAL_ROOT"
  # rehearse only looks for STOP in the data root and builds in its own scratch root: it may take another conf; no
  # mode that reads or writes the data root (estimate, before, run, after, check, backup) ever does
  if [ "$MODE" != rehearse ]; then
    [ "$(realpath -m "$CONF")" = "$(realpath -m "$OWN_CONF")" ] \
      || refuse "in the data root the conf is the tree's own $OWN_CONF, not $CONF"
    [ "$(realpath -m "$REALPROF")" = "$(realpath -m "$TREE/profiles/wilkcomp02.yaml")" ] \
      || refuse "in the data root the profile is the tree's own $TREE/profiles/wilkcomp02.yaml, not $REALPROF"
  fi
  [[ "$TMP_MIN_GB" =~ ^[0-9]+$ ]] && [ "$TMP_MIN_GB" -ge 200 ] || refuse "TMP_MIN_GB=$TMP_MIN_GB is under 200 GB"
  [[ "$START_MIN_GB" =~ ^[0-9]+$ ]] && [ "$START_MIN_GB" -ge 205 ] || refuse "START_MIN_GB=$START_MIN_GB is under 205 GB"
  [[ "$SCRATCH_MAX_GB" =~ ^[0-9]+$ ]] && [ "$SCRATCH_MAX_GB" -le 100 ] \
    || refuse "SCRATCH_MAX_GB=$SCRATCH_MAX_GB is over 100 GB"
  if [ -n "$MEM" ]; then
    [[ "$MEM" =~ ^([0-9]+)([KMGT]?)$ ]] || refuse "MEM=$MEM: expected a number of bytes with K, M, G or T"
    n=${BASH_REMATCH[1]}
    u=${BASH_REMATCH[2]}
    case "$u" in K) v=$((n * 1024)) ;; M) v=$((n * 1024 ** 2)) ;; G) v=$((n * 1024 ** 3)) ;; T) v=$((n * 1024 ** 4)) ;; *) v=$n ;; esac
    [ "$v" -le $((750 * 1024 ** 3)) ] || refuse "MEM=$MEM is above 750G"
  fi
  return 0
}

read_conf() {  # the conf, checked by stage_est.py against the tree; its values are assigned here, never evaluated
  local k v
  [ -f "$CONF" ] || refuse "no conf $CONF (a stage whose note leaves a value open has none)"
  python3 -I "$BIN/stage_est.py" conf "$CONF" --stage "$STAGE" --tree "$TREE" --sums "$C/conf.sums" > "$C/conf.env" \
    2> "$C/conf.err" || refuse "conf $CONF: $(tr '\n' ' ' < "$C/conf.err")"
  CONF_SHA=$(sha256sum "$CONF" | cut -c1-64)
  [[ "$CONF_SHA" =~ ^[0-9a-f]{64}$ ]] || refuse "could not hash $CONF"
  IN="" DEPTH="" MINLEN="" GENECALL_MODE="" LIBNAME="" LIBVER="" LIBCUT="" CRISPR="" AFTER="" NOTES=""
  while IFS='=' read -r k v; do
    case "$k" in
      INPUTS) IN=$v ;; DEPTH) DEPTH=$v ;; MIN_CONTIG_LEN) MINLEN=$v ;; GENECALL_MODE) GENECALL_MODE=$v ;;
      LIBNAME) LIBNAME=$v ;; LIBVER) LIBVER=$v ;;
      LIBCUT) LIBCUT=$v ;; CRISPR) CRISPR=$v ;; AFTER) AFTER=$v ;; NOTES) NOTES=$v ;;
    esac
  done < "$C/conf.env"
  [ -n "$IN" ] && [ -n "$DEPTH" ] && [ -n "$CRISPR" ] || refuse "conf $CONF: a value is missing"
  LIBARGS=()
  if [ -n "$LIBNAME" ]; then
    [ "$LIBNAME $LIBVER" = "pfam 38.2" ] || refuse "no HMM file is known for library $LIBNAME $LIBVER"
    LIBARGS=(--library "$LIBNAME" "$LIBVER" "$PFAM" "$LIBCUT")
  fi
  BUILD_ARGS=(--stage "$STAGE" --depth "$DEPTH")
  [ -n "$MINLEN" ] && BUILD_ARGS+=(--min-contig-len "$MINLEN")
  [ -n "$GENECALL_MODE" ] && BUILD_ARGS+=(--genecall-mode "$GENECALL_MODE")
  BUILD_ARGS+=(--genomes "$IN/genomes.json" --objects "$IN/objects.json" "${LIBARGS[@]}")
  [ "$CRISPR" = yes ] && BUILD_ARGS+=(--crispr)
  return 0
}

check_sums() {  # FILE fail|refuse: the inputs against the conf's sums, the result into FILE
  (cd "$IN" && sha256sum -c "$C/conf.sums") > "$1" 2>&1 \
    || "$2" "the inputs in $IN differ from the conf's sums: $(tr '\n' ' ' < "$1" | cut -c1-600)"
  log "the $(wc -l < "$C/conf.sums") input files match the conf's sums"
}

prof_key() { sed -n "s/^$1:[[:space:]]*\([^#[:space:]]*\).*/\1/p" "$REALPROF" | head -n 1; }

caps() {  # the thread counts against the profile: nothing here is above the profile's caps
  local v day night
  [ -f "$REALPROF" ] || refuse "no profile $REALPROF"
  day=$(prof_key max_threads_day)
  night=$(prof_key max_threads_night)
  [[ "$day" =~ ^[1-9][0-9]*$ ]] && [[ "$night" =~ ^[1-9][0-9]*$ ]] || refuse "profile $REALPROF: no thread caps"
  LOW=$(( day < night ? day : night ))
  HIGH=$(( day > night ? day : night ))
  MEM=${MEM:-$(prof_key mem_scope)}
  for v in T CPUS RH_T RH_Q CPU_MAX; do
    [[ "${!v}" =~ ^[1-9][0-9]*$ ]] || refuse "$v=${!v}: not a positive integer"
  done
  [ -z "$Q" ] || [[ "$Q" =~ ^[1-9][0-9]*$ ]] || refuse "Q=$Q: not a positive integer"
  [ -z "$SCAN_JOBS" ] || [[ "$SCAN_JOBS" =~ ^[1-9][0-9]*$ ]] || refuse "SCAN_JOBS=$SCAN_JOBS: not a positive integer"
  [ -z "$Q" ] || [ "$Q" -le "$HIGH" ] || refuse "Q=$Q is above the profile's night cap $HIGH"
  [ "$RH_Q" -le "$HIGH" ] || refuse "RH_Q=$RH_Q is above the profile's night cap $HIGH"
  for v in T CPUS RH_T; do
    [ "${!v}" -le "$LOW" ] || refuse "$v=${!v} is above the profile's smaller cap $LOW (allowed at any hour)"
  done
  [ -z "$SCAN_JOBS" ] || [ $((SCAN_JOBS * CPUS)) -le "$LOW" ] || refuse "SCAN_JOBS x CPUS is above $LOW"
  [ "$CPU_MAX" -le 280 ] || refuse "CPU_MAX=$CPU_MAX: the process guard stops at 280 CPU-minutes at most"
  [ -n "$MEM" ] || refuse "profile $REALPROF: no mem_scope"
}

profile_copy() {  # the profile with this driver's own tmp
  [ -f "$REALPROF" ] || refuse "no profile $REALPROF"
  sed "s#^tmp:.*#tmp: $C/tmp#" "$REALPROF" > "$PROF" || fail "could not write $PROF"
}

bb() {  # a short command in the image, in a capped scope (2 CPUs: HELPER_THREADS), the root bound read-only
  local what=$1 binds=(--bind /tmp/if2453 --bind "$ROOT:$ROOT:ro")
  shift
  case "$what" in
    cli) set -- python -m biobox.cli "$@" ;;
    py) set -- python "$@" ;;
  esac
  [ -d "$PFAMDIR" ] && binds+=(--bind "$PFAMDIR:$PFAMDIR:ro")
  [ -d "$STOR_RO" ] && binds+=(--bind "$STOR_RO:$STOR_RO:ro")
  systemd-run --user --scope --quiet -p MemoryMax=16G -p MemorySwapMax=0 -p CPUQuota=200% nice -n 19 \
    singularity exec --cleanenv --env PYTHONPATH="$PYPATH" --env OMP_NUM_THREADS=2 --env PYTHONUNBUFFERED=1 \
    "${binds[@]}" --bind /sys/fs/cgroup:/sys/fs/cgroup:ro "$IMAGE" "$@"
}

quota_now() {  # QUOTA_NOW: the profile's cap of this hour less the helpers' 2, at most Q when given (a background call)
  local q
  QUOTA_NOW=""
  bgcall "$C/tmp/quota.out" bb py "$BIN/stage_est.py" cap --profile "$REALPROF" --quota || return 1
  q=$(head -n 1 "$C/tmp/quota.out")
  [[ "$q" =~ ^[1-9][0-9]*$ ]] || return 1
  [ -n "$Q" ] && [ "$Q" -lt "$q" ] && q=$Q
  QUOTA_NOW=$q
}

lock_patterns() {  # the units that write into this root: every build of ours in the data root, else this stage's own
  local f s
  echo "$STEM-[0-9]*"
  [ "$DATAROOT" = 1 ] || return 0
  printf '%s\n' 'r50-*' 'c24-*' 'c26-*' 'c33-*' 'c34-*' 'c35-*' 'c46-*' 'c47-c5-*' 'c48c-*' 'b47-*'
  for f in "$BIN"/stages/*.conf; do
    [ -e "$f" ] || continue
    s=$(basename "$f" .conf)
    [[ "$s" =~ ^[A-Z][A-Za-z0-9]{0,15}$ ]] && echo "${s,,}-[0-9]*"
  done
  return 0
}

build_units() {
  local pats
  mapfile -t pats < <(lock_patterns | LC_ALL=C sort -u)
  systemctl --user list-units --type=service --state=active,activating,deactivating,reloading --no-legend --plain \
    "${pats[@]}" | awk '{print $1}' | paste -sd' ' -
}

snap_top() {
  log "uptime: $(uptime)"
  log "df /tmp: $(df -h /tmp | tail -n 1)"
  log "active units: $(systemctl --user list-units --type=service --type=scope --state=active --no-legend --plain | awk '{print $1}' | paste -sd' ' -)"
}

free_gb() { df -B1G --output=avail /tmp | tail -n 1 | tr -d ' '; }
exit_of() { sed -n "s/^pod_run: $1 exited \([0-9]*\)$/\1/p" "$2/$1.log" 2> /dev/null | tail -n 1; }
hash_of() { sed -n 's/.*manifest hash \([0-9a-f]*\).*/\1/p' "$1" 2> /dev/null | tail -n 1; }

unit_ended() {  # UNIT LOGS: 0 on the unit's own exit line or a real answer of systemd that it is inactive or failed;
  # 1 on a real answer that it runs; 2 on no answer (the user bus or systemctl failed): it counts as running
  local s
  [ -n "$(exit_of "$1" "$2")" ] && return 0
  if ! s=$(systemctl --user show -p ActiveState --value "$1" 2> /dev/null) || [ -z "$s" ]; then
    log "no answer from systemd about $1; it counts as running"
    return 2
  fi
  case "$s" in inactive|failed) return 0 ;; *) return 1 ;; esac
}

# How long the driver goes on without one readable answer about its unit before it gives up (the user bus may be
# gone): it then says the unit's state is unknown and that it may still be running, and exits non-zero. Never claims
# it stopped. 15 minutes: longer than a busy bus's stalls, short enough that a queue does not hang for hours.
WATCH_LIMIT_S=${WATCH_LIMIT_S:-900}
STOP_LIMIT_S=${STOP_LIMIT_S:-1800}  # how long stop_unit asks before it says the unit would not stop
NOANSWER_SINCE=""

watch_answer() {  # RC of unit_ended: 0 ended, 1 running; on 2 (no answer) give up once WATCH_LIMIT_S passed without one
  case "$1" in
    0) return 0 ;;
    1) NOANSWER_SINCE="" ; return 1 ;;
  esac
  [ -n "$NOANSWER_SINCE" ] || NOANSWER_SINCE=$(date +%s)
  if [ $(( $(date +%s) - NOANSWER_SINCE )) -ge "$WATCH_LIMIT_S" ]; then
    RUNNING_UNIT=""
    fail "no readable answer from systemd about unit $2 for ${WATCH_LIMIT_S} s: its state is unknown and it may still be running (systemctl --user list-units '$STEM-[0-9]*'; a STOP file in the root ends it at a step boundary)"
  fi
  return 1
}

stop_unit() {  # UNIT LOGS WHY: stop it and see it gone before saying so; fail loudly at its bound, never claiming it
  local k=0 t0 rc
  t0=$(date +%s)
  log "stopping $1: $3"
  while :; do
    unit_ended "$1" "$2"
    rc=$?
    [ "$rc" = 0 ] && break
    watch_answer "$rc" "$1"
    if [ $(( $(date +%s) - t0 )) -ge "$STOP_LIMIT_S" ]; then
      RUNNING_UNIT=""
      fail "$3; then unit $1 was asked to stop for ${STOP_LIMIT_S} s and is still not seen gone: it may still be running"
    fi
    if [ $((k % 15)) = 0 ]; then
      systemctl --user stop "$1" || log "systemctl --user stop $1 failed; asking again"
    fi
    k=$((k + 1))
    sleep 2
  done
  RUNNING_UNIT=""
  fail "$3 (unit $1 stopped)"
}

bgcall() {  # OUT CMD...: CMD in the background, its output into OUT, waited for here, so that a TERM to the driver runs
  # its trap at once (bash runs a trap only after a foreground command returns, and systemd kills the driver 90 s after
  # a stop by default)
  local out=$1
  shift
  "$@" > "$out" 2>> "$LOG" &
  BGPID=$!
  wait "$BGPID"
}

on_exit() {  # the driver is leaving while its unit runs (a signal, a crash): stop the unit, as far as it can
  local rc=$? k=0
  [ -n "${BGPID:-}" ] && kill "$BGPID" 2> /dev/null
  if [ -n "$RUNNING_UNIT" ]; then
    log "the driver exits (status $rc) while $RUNNING_UNIT runs: stopping it"
    while ! unit_ended "$RUNNING_UNIT" "$RUNNING_LOGS" && [ "$k" -lt 60 ]; do
      [ $((k % 15)) = 0 ] && { systemctl --user stop "$RUNNING_UNIT" || log "systemctl --user stop failed"; }
      k=$((k + 1))
      sleep 2
    done
    unit_ended "$RUNNING_UNIT" "$RUNNING_LOGS" || log "$RUNNING_UNIT may still run: stop it by hand"
  fi
}
trap on_exit EXIT
trap 'exit 143' TERM
trap 'exit 130' INT
trap 'exit 129' HUP

releases() {  # the root's published releases, as catalog.manifest_releases reads them: manifests/release=<R>.json
  [ -d "$ROOT/manifests" ] || return 0
  find "$ROOT/manifests" -maxdepth 1 -type f -name 'release=*.json' -printf '%f\n' \
    | sed -e 's/^release=//' -e 's/\.json$//' | LC_ALL=C sort
}

verify() {
  bb cli release verify --root "$ROOT" --release "$1" > "$2/verify-$1.log" 2>&1
  local rc=$?
  log "verify $1: exit $rc: $(grep -E '^release ' "$2/verify-$1.log" | tail -n 1)"
  return $rc
}

snapshot() {  # DIR before|after: the listing (before: without the stage's own files) and every other stage's sha256
  local d=$1 n_all n_left
  mkdir -p "$d"
  find "$ROOT" ! -type d -printf '%P\t%s\t%T@\n' > "$d/listing.raw" || fail "snapshot $d: find failed"
  LC_ALL=C sort "$d/listing.raw" > "$d/listing.all.tsv" || fail "snapshot $d: sort failed"
  rm -f "$d/listing.raw"
  if [ "$2" = before ]; then
    awk -F '\t' -v s="/stage=$STAGE/" 'index("/" $1 "/", s) == 0' "$d/listing.all.tsv" > "$d/listing.tsv" \
      || fail "snapshot $d: filter failed"
  else
    cp "$d/listing.all.tsv" "$d/listing.tsv" || fail "snapshot $d: copy failed"
  fi
  awk -F '\t' -v s="/stage=$STAGE/" 'index("/" $1 "/", s) == 0 && ("/" $1) ~ /\/stage=[^\/]+\// {print $1}' \
    "$d/listing.all.tsv" > "$d/earlier.files" || fail "snapshot $d: list of other stages' files failed"
  (cd "$ROOT" && systemd-run --user --scope --quiet -p MemoryMax=2G -p MemorySwapMax=0 nice -n 19 ionice -c 3 \
    xargs -r -d '\n' sha256sum < "$d/earlier.files") > "$d/earlier.sha256" || fail "snapshot $d: sha256sum failed"
  [ "$(wc -l < "$d/earlier.sha256")" = "$(wc -l < "$d/earlier.files")" ] || fail "snapshot $d: not every file hashed"
  n_all=$(wc -l < "$d/listing.all.tsv")
  n_left=$(wc -l < "$d/listing.tsv")
  log "snapshot $d: $n_all files, $n_left listed ($((n_all - n_left)) of stage=$STAGE left out), $(wc -l < "$d/earlier.sha256") of other stages hashed"
}

# --------------------------------------------------------------------------------------------------------------
# A wave's shared checks (card A-30, DEC-81 (3)); see the header. Every function here answers "no" on any doubt.

share_on() {  # 0 when this mode may take a check from an earlier stage of its wave; never in a rehearsal's own modes
  [ -n "$WAVE" ] || return 1
  [ "$REHEARSAL" = 1 ] && return 1
  if [ "$SHARE" != 1 ]; then
    log "share: off (STAGE_SHARE=$SHARE): every check of this mode runs in full"
    return 1
  fi
  if ! [[ "$WAVE" =~ ^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$ ]]; then
    log "share: off (STAGE_WAVE '$WAVE' is not a wave's name): every check of this mode runs in full"
    return 1
  fi
  return 0
}

sha_of() { sha256sum "$1" 2> /dev/null | cut -c1-64; }
man_sha() { sha_of "$ROOT/manifests/release=$1.json"; }

code_sha() {  # CODE_SHA: the sha256 of a list "<sha256>  <tree|bin>/<path>", one line per file, sorted in byte order, of
  # every file under $TREE/src and $TREE/tests/fixtures/micro_s0 but __pycache__/ and *.pyc, every file beside this
  # script ($BIN, not stages/), and $TREE/scripts/pod_run.sh: the code a rehearsal and a verification run. The list is
  # kept in tmp/code.sha256. 1 when a part is missing or cannot be read.
  local out=$C/tmp/code.sha256
  CODE_SHA=""
  [ -d "$TREE/src" ] && [ -d "$TREE/tests/fixtures/micro_s0" ] && [ -f "$TREE/scripts/pod_run.sh" ] || return 1
  (cd "$TREE" && find src tests/fixtures/micro_s0 -type f ! -name '*.pyc' ! -path '*/__pycache__/*' -print0 \
    | xargs -0 -r sha256sum | sed 's#  #  tree/#') > "$out.part" || return 1
  (cd "$BIN" && find -L . -maxdepth 1 -type f -print0 | xargs -0 -r sha256sum | sed 's#  \./#  bin/#') >> "$out.part" \
    || return 1
  (cd "$TREE" && sha256sum scripts/pod_run.sh | sed 's#  #  tree/#') >> "$out.part" || return 1
  LC_ALL=C sort "$out.part" > "$out" || return 1
  rm -f "$out.part"
  CODE_SHA=$(sha_of "$out")
  [[ "$CODE_SHA" =~ ^[0-9a-f]{64}$ ]] || { CODE_SHA=""; return 1; }
}

rh_facts() {  # FILE: what a rehearsal of this conf rests on, one fact a line; 1 when one cannot be read
  local prof
  code_sha || return 1
  prof=$(sha_of "$REALPROF")
  [[ "$prof" =~ ^[0-9a-f]{64}$ ]] || return 1
  printf '%s\n' "wave=$WAVE" "tree=$(realpath -m "$TREE")" "code=$CODE_SHA" "image=$IMAGE" "digest=$DIGEST" \
    "profile=$prof" \
    "kind=DEPTH=$DEPTH MIN_CONTIG_LEN=$MINLEN GENECALL_MODE=$GENECALL_MODE LIBNAME=$LIBNAME LIBVER=$LIBVER LIBCUT=$LIBCUT CRISPR=$CRISPR" \
    "limits=RH_T=$RH_T RH_Q=$RH_Q CPU_MAX=$CPU_MAX TMP_MIN_GB=$TMP_MIN_GB START_MIN_GB=$START_MIN_GB SCRATCH_MAX_GB=$SCRATCH_MAX_GB MEM=$MEM" \
    > "$1"
}

rh_share() {  # 0 when another stage's full rehearsal of this wave rests on the same facts: its marker is taken, with them
  local want=$C/rehearse/facts f src="" s="" n=0
  if ! rh_facts "$want"; then
    rm -f "$want"
    log "rehearse: not shared: a fact of this rehearsal (the code's or the profile's sha256) could not be read; rehearsing in full"
    return 1
  fi
  for f in "$STAGES_HOME"/*/rehearse/DONE; do
    [ -f "$f" ] || continue
    [ "$(realpath -m "$f")" = "$(realpath -m "$C/rehearse/DONE")" ] && continue
    n=$((n + 1))
    [ "$(sed -n 2p "$f")" = "rehearsal=full" ] || continue
    [[ "$(head -n 1 "$f")" =~ ^[0-9a-f]{64}$ ]] || continue
    s=$(sed -n '3s/^stage=//p' "$f")
    [[ "$s" =~ ^[A-Z][A-Za-z0-9]{0,15}$ ]] || continue
    tail -n +4 "$f" | cmp -s - "$want" || continue
    src=$f
    break
  done
  if [ -z "$src" ]; then
    log "rehearse: not shared: no full rehearsal of wave $WAVE by another stage rests on the same facts ($n markers read under $STAGES_HOME); rehearsing in full on: $(paste -sd ';' "$want")"
    return 1
  fi
  { echo "$CONF_SHA"; echo "rehearsal=shared from stage $s ($src)"; echo "stage=$STAGE"; cat "$want"; } \
    > "$C/rehearse/DONE.tmp" && mv -f "$C/rehearse/DONE.tmp" "$C/rehearse/DONE" || fail "could not write $C/rehearse/DONE"
  log "rehearse: SHARED from stage $s ($src): its full rehearsal of wave $WAVE rests on the same facts, each equal: $(paste -sd ';' "$want"); nothing is rehearsed here"
  return 0
}

nos() { log "before: not shared${REC:+ ($REC)}: $*; every release is verified in full"; }
fld() { sed -n "s/^$1=//p" "$REC" | head -n 1; }

verify_record() {  # D MODE: after a full verification in a wave, the facts a later stage's before compares; never fails
  local d=$1 r s h l e tmp=$1/verified.share.tmp
  rm -f "$d/verified.share" "$tmp"
  if ! code_sha; then
    log "$2: no record for wave $WAVE: the code's sha256 could not be taken"
    return 0
  fi
  l=$(sha_of "$d/listing.all.tsv")
  e=$(sha_of "$d/earlier.sha256")
  if ! [[ "$l" =~ ^[0-9a-f]{64}$ ]] || ! [[ "$e" =~ ^[0-9a-f]{64}$ ]]; then
    log "$2: no record for wave $WAVE: its listing or earlier.sha256 could not be hashed"
    return 0
  fi
  printf '%s\n' "wave=$WAVE" "stage=$STAGE" "mode=$2" "conf=$CONF_SHA" "root=$ROOT" "tree=$(realpath -m "$TREE")" \
    "code=$CODE_SHA" "image=$IMAGE" "digest=$DIGEST" "listing=$l" "earlier=$e" > "$tmp" || return 0
  while IFS= read -r r; do
    s=$(awk -v r="$r" '$2 == r {print $1}' "$d/manifests.sha256" 2> /dev/null)
    h=$(hash_of "$d/verify-$r.log")
    if ! [[ "$s" =~ ^[0-9a-f]{64}$ ]] || ! [[ "$h" =~ ^[0-9a-f]+$ ]]; then
      rm -f "$tmp"
      log "$2: no record for wave $WAVE: release $r has no manifest sha256 ('$s') or manifest hash ('$h')"
      return 0
    fi
    echo "release=$r $s $h" >> "$tmp" || { rm -f "$tmp"; return 0; }
  done < "$d/releases.list"
  mv -f "$tmp" "$d/verified.share" || { rm -f "$tmp"; log "$2: no record for wave $WAVE: could not write it"; return 0; }
  log "$2: recorded for wave $WAVE in $d/verified.share: $(wc -l < "$d/releases.list") releases verified in full, listing sha256 $l, earlier.sha256 sha256 $e, code sha256 $CODE_SHA"
}

verify_share() {  # D: 0 when the newest record of this wave proves every release of D/releases.list verified on the root as
  # it is now (D's snapshot taken); each verify-<R>.log then says so. Else 1, its reason logged.
  local d=$1 t f src st md r s h line k n
  REC=""
  while IFS=$'\t' read -r t f; do
    [ "$(sed -n 's/^wave=//p' "$f" | head -n 1)" = "$WAVE" ] && { REC=$f; break; }
  done < <(find "$STAGES_HOME" -mindepth 3 -maxdepth 3 -type f -name verified.share \( -path '*/after/*' -o -path '*/before/*' \) \
    -printf '%T@\t%p\n' 2> /dev/null | LC_ALL=C sort -t "$(printf '\t')" -k1,1gr)
  [ -n "$REC" ] || { nos "no record of a full verification in wave $WAVE under $STAGES_HOME"; return 1; }
  src=$(dirname "$REC")
  for k in wave stage mode conf root tree code image digest listing earlier; do
    [ -n "$(fld "$k")" ] || { nos "the record is incomplete (no $k)"; return 1; }
  done
  st=$(fld stage)
  md=$(fld mode)
  [[ "$st" =~ ^[A-Z][A-Za-z0-9]{0,15}$ ]] && [[ "$md" =~ ^(after|before)$ ]] && [ "$(basename "$src")" = "$md" ] \
    || { nos "the record is malformed (stage '$st', mode '$md')"; return 1; }
  for k in conf code listing earlier; do
    [[ "$(fld "$k")" =~ ^[0-9a-f]{64}$ ]] || { nos "the record is malformed ($k)"; return 1; }
  done
  [ -f "$src/DONE" ] && [ "$(head -c 64 "$src/DONE")" = "$(fld conf)" ] \
    || { nos "stage $st's $md/DONE is missing or made for another conf than the record names"; return 1; }
  [ "$(fld root)" = "$ROOT" ] || { nos "it was made for root $(fld root), not $ROOT"; return 1; }
  [ "$(fld tree)" = "$(realpath -m "$TREE")" ] || { nos "it was made from tree $(fld tree), not $(realpath -m "$TREE")"; return 1; }
  [ "$(fld image)" = "$IMAGE" ] && [ "$(fld digest)" = "$DIGEST" ] \
    || { nos "the image or its digest differs (recorded $(fld image) $(fld digest), now $IMAGE $DIGEST)"; return 1; }
  code_sha || { nos "the code's sha256 could not be taken"; return 1; }
  [ "$(fld code)" = "$CODE_SHA" ] || { nos "the code's sha256 differs (recorded $(fld code), now $CODE_SHA)"; return 1; }
  [ "$(sha_of "$src/listing.all.tsv")" = "$(fld listing)" ] \
    || { nos "its listing.all.tsv no longer has the sha256 it recorded"; return 1; }
  [ "$(sha_of "$d/listing.all.tsv")" = "$(fld listing)" ] \
    || { nos "the root's listing (name, size and mtime of every file) differs from the one it recorded"; return 1; }
  [ "$(sha_of "$src/earlier.sha256")" = "$(fld earlier)" ] \
    || { nos "its earlier.sha256 no longer has the sha256 it recorded"; return 1; }
  LC_ALL=C sort "$src/earlier.sha256" > "$d/share.recorded.sha256" && LC_ALL=C sort "$d/earlier.sha256" > "$d/share.now.sha256" \
    && LC_ALL=C comm -23 "$d/share.recorded.sha256" "$d/share.now.sha256" > "$d/share.differ" \
    || { nos "the sha256 of the other stages' files could not be compared"; return 1; }
  [ ! -s "$d/share.differ" ] \
    || { nos "the sha256 of $(wc -l < "$d/share.differ") files of other stages that stage $st recorded differ or are gone (first: $(head -n 1 "$d/share.differ" | cut -c1-200))"; return 1; }
  f=$(find "$ROOT" ! -type d -newer "$REC" -print -quit) || { nos "the root's files could not be compared with the record's time"; return 1; }
  [ -z "$f" ] || { nos "a file of the root is newer than the record: ${f#"$ROOT"/}"; return 1; }
  sed -n 's/^release=\([^ ]*\) .*/\1/p' "$REC" | LC_ALL=C sort > "$d/share.releases" \
    || { nos "the record's releases could not be read"; return 1; }
  cmp -s "$d/share.releases" "$d/releases.list" \
    || { nos "it records the releases $(paste -sd ' ' "$d/share.releases"), the root has $(paste -sd ' ' "$d/releases.list")"; return 1; }
  while IFS= read -r r; do
    line=$(awk -v r="release=$r" '$1 == r' "$REC")
    [ "$(printf '%s\n' "$line" | wc -l)" = 1 ] || { nos "the record names release $r more than once"; return 1; }
    s=$(echo "$line" | awk '{print $2}')
    h=$(echo "$line" | awk '{print $3}')
    [[ "$s" =~ ^[0-9a-f]{64}$ ]] && [[ "$h" =~ ^[0-9a-f]+$ ]] || { nos "the record's line of release $r is malformed"; return 1; }
    [ "$(man_sha "$r")" = "$s" ] || { nos "the manifest of release $r has sha256 '$(man_sha "$r")', the record $s"; return 1; }
  done < "$d/releases.list"
  n=$(wc -l < "$d/releases.list")
  while IFS= read -r r; do
    line=$(awk -v r="release=$r" '$1 == r' "$REC")
    s=$(echo "$line" | awk '{print $2}')
    h=$(echo "$line" | awk '{print $3}')
    echo "release $r: not verified here: shared from stage $st mode $md of wave $WAVE ($REC), which verified it with manifest sha256 $s, manifest hash $h" \
      > "$d/verify-$r.log" || fail "could not write $d/verify-$r.log"
    log "verify $r: SHARED from stage $st mode $md (wave $WAVE): manifest sha256 $s equal, manifest hash $h"
  done < "$d/releases.list"
  log "before: SHARED the verification of $n releases from stage $st mode $md of wave $WAVE ($REC); equal: root $ROOT, tree $(fld tree), code sha256 $CODE_SHA, image $IMAGE $DIGEST, the root's listing (sha256 $(fld listing): name, size and mtime of each of its $(wc -l < "$d/listing.all.tsv") files), the sha256 of the $(wc -l < "$src/earlier.sha256") files of other stages that stage $st recorded, each release's manifest sha256; no file of the root newer than the record"
  return 0
}

# --------------------------------------------------------------------------------------------------------------

estimate() {
  local d=$C/estimate rc
  withdraw estimate
  read_conf
  profile_copy
  check_sums "$d/SUMS.check" fail
  bb py "$BIN/stage_est.py" before --conf "$CONF" --stage "$STAGE" --root "$ROOT" --profile "$PROF" --inputs "$IN" \
    --tmp "$C/tmp/est" --out-plan "$d/genecall_plan.json" > "$d/before.json" 2> "$d/before.err"
  rc=$?
  log "estimate: exit $rc: $(tr -d '\n ' < "$d/before.json" | cut -c1-1500) $(tail -n 2 "$d/before.err")"
  [ "$rc" = 5 ] && fail "a contig of $STAGE is owned by a stage that is not gene-called, or by two: $(tail -n 1 "$d/before.err")"
  [ "$rc" = 7 ] && fail "build-stage would refuse the stage: $(tail -n 1 "$d/before.err")"
  [ "$rc" = 8 ] && fail "the contigs earlier stages hold differ from the note's: $(tail -n 1 "$d/before.err")"
  [ "$rc" = 0 ] || fail "stage_est.py before exited $rc (read $d/before.err)"
  [ -s "$d/before.json" ] || fail "stage_est.py before wrote nothing"
  mark_done estimate
}

before() {
  local d=$C/before r ok=1 act
  withdraw before
  read_conf
  rm -f "$d/verified.share"
  act=$(build_units) || fail "could not list the active units"
  [ -z "$act" ] || refuse "a build unit is active in the root: $act"
  snap_top
  releases > "$d/releases.list" || fail "could not list $ROOT/manifests"
  [ -s "$d/releases.list" ] || refuse "no published release in $ROOT/manifests: is $ROOT the data root?"
  if [ -d "$ROOT/release_def" ]; then
    find "$ROOT/release_def" -maxdepth 1 -type d -name 'release=*' -printf '%f\n' | sed 's/^release=//' \
      | LC_ALL=C sort | LC_ALL=C comm -23 - "$d/releases.list" > "$d/defined_without_manifest.list"
    [ -s "$d/defined_without_manifest.list" ] \
      && log "before: defined in release_def without a manifest, not verified: $(paste -sd' ' - < "$d/defined_without_manifest.list")"
  fi
  if share_on; then  # in a wave: the snapshot first, then the verification, taken (verify_share) or in full and recorded
    while IFS= read -r r; do
      [[ "$r" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || fail "a release name in $ROOT/manifests is not one: '$r'"
    done < "$d/releases.list"
    snapshot "$d" before
    if verify_share "$d"; then
      mark_done before
      return 0
    fi
    : > "$d/manifests.sha256" || fail "could not write $d/manifests.sha256"
    while IFS= read -r r; do
      echo "$(man_sha "$r") $r" >> "$d/manifests.sha256"
      verify "$r" "$d" || ok=0
    done < "$d/releases.list"
    [ "$ok" = 1 ] || fail "a published release failed verify before launch (read $d/verify-*.log)"
    verify_record "$d" before
    mark_done before
    return 0
  fi
  while IFS= read -r r; do
    [[ "$r" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || fail "a release name in $ROOT/manifests is not one: '$r'"
    verify "$r" "$d" || ok=0
  done < "$d/releases.list"
  [ "$ok" = 1 ] || fail "a published release failed verify before launch (read $d/verify-*.log)"
  snapshot "$d" before
  mark_done before
}

next_unit() {
  local n
  n=$(( $( [ -e "$C/units" ] && wc -l < "$C/units" || echo 0) + 1 ))
  echo "$STEM-$n"
}

stage_records() { find "$ROOT/runs/stage=$STAGE" -maxdepth 1 -name '*.json' ! -name '*.publishing.json' 2> /dev/null | wc -l; }

rates_now() {  # RATES_LINE, RATES_RC: the rate guard once (a background call), its time in seconds logged
  local t0
  t0=$(date +%s.%N)
  bgcall "$C/tmp/rates.out" bb py "$BIN/stage_est.py" rates --conf "$CONF" --stage "$STAGE" --root "$ROOT" \
    --profile "$PROF" --objects "$IN/objects.json" --before "$C/estimate/before.json" --tmp "$C/tmp/rates"
  RATES_RC=$?
  RATES_LINE=$(cat "$C/tmp/rates.out" 2> /dev/null)
  log "rates call: exit $RATES_RC in $(awk -v a="$t0" -v b="$(date +%s.%N)" 'BEGIN {printf "%.1f", b - a}') s"
  echo "$RATES_LINE" >> "$C/run/rates.jsonl"
}

stopped_root() { [ -e "$ROOT/STOP" ] || [ -e "$REAL_ROOT/STOP" ]; }

run() {
  local unit free act rc line nrec=-1 n k=0 foreign=0 a sj=() binds=(--bind /tmp/if2453) quota hour h
  withdraw run
  read_conf
  [ -e "$C/run/FAILED" ] && refuse "$C/run/FAILED exists ($(head -c 300 "$C/run/FAILED")): a stop stays a stop; a person reads it and moves it aside (run/FAILED.<unit>) before a new launch"
  need estimate
  need before
  if [ "$REHEARSAL" = 1 ]; then
    [ "$DATAROOT" = 0 ] || refuse "a rehearsal never builds in the data root"
  else
    need rehearse
  fi
  for a in $AFTER; do
    [ -e "$STAGES_HOME/${a,,}/backup/DONE" ] || refuse "stage $STAGE comes after $a, and $STAGES_HOME/${a,,}/backup/DONE is not there"
  done
  caps
  profile_copy
  check_sums "$C/run/SUMS.check" refuse  # again: nothing is launched on inputs that changed since the estimate
  act=$(build_units) || fail "could not list the active units"
  [ -z "$act" ] || refuse "a build unit is active in the root: $act"
  free=$(free_gb)
  [[ "$free" =~ ^[0-9]+$ ]] && [ "$free" -ge "$START_MIN_GB" ] || refuse "/tmp has '$free' GB free, under $START_MIN_GB"
  quota_now || refuse "could not read the profile's cap of this hour (stage_est.py cap)"
  quota=$QUOTA_NOW
  hour=$(date +%H)
  snap_top
  unit=$(next_unit)
  echo "$unit" >> "$C/units"
  touch "$C/run/since"
  [ -n "$SCAN_JOBS" ] && sj=(--scan-jobs "$SCAN_JOBS")
  [ -n "$LIBNAME" ] && binds+=(--bind "$PFAMDIR:$PFAMDIR:ro")
  log "launch $unit (pod_run --threads $quota --mem $MEM): build-stage --profile $PROF --root $ROOT ${BUILD_ARGS[*]} --threads $T ${sj[*]} --scan-cpus $CPUS"
  RUNNING_UNIT=$unit
  RUNNING_LOGS=$LOGS
  bash "$TREE/scripts/pod_run.sh" --profile "$PROF" --image "$IMAGE" --digest "$DIGEST" --mem "$MEM" --threads "$quota" \
    "${binds[@]}" --env "PYTHONPATH=$PYPATH" --log "$LOGS/$unit.log" "$unit" \
    build-stage --profile "$PROF" --root "$ROOT" "${BUILD_ARGS[@]}" --threads "$T" "${sj[@]}" --scan-cpus "$CPUS" \
    >> "$LOG" 2>&1 || { RUNNING_UNIT=""; fail "pod_run.sh refused $unit"; }
  sleep 10
  while :; do
    unit_ended "$unit" "$LOGS"
    rc=$?
    [ "$rc" = 0 ] && break
    watch_answer "$rc" "$unit"
    if [ "$foreign" = 0 ] && stopped_root; then
      foreign=1
      log "STOP appeared while $unit runs; build-stage ends at the next step boundary"
    fi
    if [ $((k % 6)) = 0 ]; then
      bgcall "$C/tmp/guard.out" python3 -I "$BIN/r50_watch.py" --root "$ROOT" --scratch "$SCRATCH" \
        --since "$C/run/since" --state "$C/run" --unit-prefix "$STEM-" --tmp-min-gb "$TMP_MIN_GB" \
        --scratch-max-gb "$SCRATCH_MAX_GB" --cpu-max "$CPU_MAX" --stages "$STAGE"
      rc=$?
      line=$(cat "$C/tmp/guard.out" 2> /dev/null)
      echo "$line" >> "$C/run/guard.jsonl"
      [ "$rc" = 0 ] || stop_unit "$unit" "$LOGS" "guard (r50_watch.py exit $rc): $(echo "$line" | cut -c1-600)"
      [ $((k % 60)) = 0 ] && log "guard: $(echo "$line" | cut -c1-400); log: $(grep -E 'run-step|build-stage' "$LOGS/$unit.log" | tail -n 1 | cut -c1-240)"
    fi
    h=$(date +%H)
    if [ $((k % 6)) = 0 ] || [ "$h" != "$hour" ]; then  # the cap changes on the hour: looked at within one poll
      hour=$h
      if quota_now; then
        if [ "$QUOTA_NOW" != "$quota" ]; then
          if systemctl --user set-property --runtime "$unit" "CPUQuota=$((QUOTA_NOW * 100))%"; then
            log "the cap of the hour changed: $unit's CPUQuota from $((quota * 100))% to $((QUOTA_NOW * 100))%"
            quota=$QUOTA_NOW
          else
            log "systemctl set-property $unit CPUQuota=$((QUOTA_NOW * 100))% failed; trying again at the next poll"
          fi
        fi
      else
        log "could not read the cap of this hour; the quota stays $((quota * 100))%"
      fi
    fi
    n=$(stage_records)
    if [ "$n" != "$nrec" ]; then
      nrec=$n
      rates_now
      log "rates ($n records of $STAGE): exit $RATES_RC: $(echo "$RATES_LINE" | cut -c1-500)"
      [ "$RATES_RC" = 0 ] || stop_unit "$unit" "$LOGS" "rate guard (stage_est.py rates exit $RATES_RC; 3: a judged step more than 50% over its estimate, DEC-64): $(echo "$RATES_LINE" | cut -c1-900)"
    fi
    k=$((k + 1))
    sleep 10
  done
  RUNNING_UNIT=""
  rc=$(exit_of "$unit" "$LOGS")
  log "$unit ended, exit '$rc'; usage: $(grep -E 'pod_run: .* usage' "$LOGS/$unit.log" | tail -n 1)"
  rates_now
  log "rates at the end: exit $RATES_RC: $(echo "$RATES_LINE" | cut -c1-600)"
  [ "$RATES_RC" = 0 ] || fail "rate guard at the end (stage_est.py rates exit $RATES_RC): $(echo "$RATES_LINE" | cut -c1-900)"
  if [ "$foreign" = 1 ] || stopped_root; then
    log "STOPPED: a STOP file ended $unit at a step boundary (exit '$rc'), its published records judged; remove STOP and run again to resume"
    echo "STOP during $unit" > "$C/run/STOPPED"
    exit 3
  fi
  [ "$rc" = 0 ] || fail "$unit exited '$rc', expected 0 (read $LOGS/$unit.log)"
  bb py "$BIN/stage_est.py" complete --conf "$CONF" --stage "$STAGE" --root "$ROOT" > "$C/run/complete.json" 2>&1 \
    || fail "stage $STAGE is not complete at depth $DEPTH: $(cut -c1-600 "$C/run/complete.json")"
  mark_done run
}

after() {  # always verifies in full; in a wave it records that verification for the next stage's before (card A-30)
  local d=$C/after r ok=1 rc h0 h1 act sh=0
  withdraw after
  read_conf
  rm -f "$d/verified.share"
  need run
  need before
  [ -s "$C/before/releases.list" ] || refuse "no $C/before/releases.list"
  act=$(build_units) || fail "could not list the active units"
  [ -z "$act" ] || refuse "a build unit is active in the root: $act"
  snapshot "$d" after
  releases > "$d/releases.list" || fail "could not list $ROOT/manifests"
  cmp -s "$C/before/releases.list" "$d/releases.list" \
    || fail "the root's published releases differ from before: before $(paste -sd' ' - < "$C/before/releases.list"); now $(paste -sd' ' - < "$d/releases.list")"
  if share_on; then
    sh=1
    : > "$d/manifests.sha256" || fail "could not write $d/manifests.sha256"
  fi
  while IFS= read -r r; do
    [ "$sh" = 1 ] && echo "$(man_sha "$r") $r" >> "$d/manifests.sha256"
    verify "$r" "$d" || ok=0
    h0=$(hash_of "$C/before/verify-$r.log")
    h1=$(hash_of "$d/verify-$r.log")
    if [ -z "$h1" ] || [ "$h0" != "$h1" ]; then
      log "after: release $r manifest hash '$h1', before '$h0'"
      ok=0
    fi
  done < "$d/releases.list"
  [ "$ok" = 1 ] || fail "a published release failed verify after the build or its manifest hash changed"
  bb cli releases --root "$ROOT" > "$d/releases.txt" 2>&1
  rc=$?
  log "after: bio-box releases exit $rc"
  [ "$rc" = 0 ] || fail "bio-box releases exited $rc (read $d/releases.txt)"
  [ "$sh" = 1 ] && verify_record "$d" after
  mark_done after
}

check() {  # its report in check/<n>/, one numbered directory per run, kept
  local n d rc ex=()
  withdraw check
  read_conf
  need after
  profile_copy
  n=$(find "$C/check" -mindepth 1 -maxdepth 1 -type d -regex '.*/[0-9]+' | wc -l)
  d=$C/check/$((n + 1))
  mkdir -p "$d" || fail "could not make $d"
  [ "$REHEARSAL" = 1 ] && [ -n "$LIBNAME" ] && ex=(--expect-copy)
  bb py "$BIN/stage_check.py" --conf "$CONF" --stage "$STAGE" --root "$ROOT" --before "$C/before" --after "$C/after" \
    --objects "$IN/objects.json" --out "$d" --tmp "$C/tmp/check" --profile "$PROF" \
    --before-est "$C/estimate/before.json" "${ex[@]}" > "$d/check.out" 2>&1
  rc=$?
  log "check ($d): exit $rc: $(tail -n 1 "$d/check.out" | cut -c1-600)"
  [ "$rc" = 0 ] || fail "stage_check.py exited $rc (read $d/check.out)"
  [ -s "$d/root_added.tsv" ] && [ "$(wc -l < "$d/root_added.tsv")" -gt 1 ] || fail "no file listed in $d/root_added.tsv"
  if [ -e "$C/backup.tsv" ] && [ -e "$C/root_added.tsv" ]; then
    cmp -s "$C/root_added.tsv" "$d/root_added.tsv" \
      || fail "$d/root_added.tsv differs from $C/root_added.tsv, the list the backup was made from"
    log "check: root_added.tsv is byte-identical to the list the backup was made from"
  else
    cp "$d/root_added.tsv" "$C/root_added.tsv" || fail "could not copy root_added.tsv"
  fi
  mark_done check
}

backup() {
  local name=added_$LC.tar
  local dst_dir tar tmp local_tar=$C/$name rec_tar want_tar backup_real small_real sz sha1 sha2 n
  withdraw backup
  read_conf
  need check
  [ -s "$C/root_added.tsv" ] || refuse "no $C/root_added.tsv"
  while [[ "$SMALL_BACKUP_DIR" == ?*/ ]]; do SMALL_BACKUP_DIR=${SMALL_BACKUP_DIR%/}; done
  [[ "$TAR_FLOOR_BYTES" =~ ^[0-9]+$ ]] || refuse "TAR_FLOOR_BYTES=$TAR_FLOOR_BYTES: not a number"
  [[ "$SMALL_TAR_FLOOR_BYTES" =~ ^[0-9]+$ ]] || refuse "SMALL_TAR_FLOOR_BYTES=$SMALL_TAR_FLOOR_BYTES: not a number"
  if [ -n "$SMALL_BACKUP_DIR" ]; then
    small_real=$(realpath -m -- "$SMALL_BACKUP_DIR") || refuse "could not resolve SMALL_BACKUP_DIR=$SMALL_BACKUP_DIR"
    case "$SMALL_BACKUP_DIR" in
      /stor|/stor/*) refuse "SMALL_BACKUP_DIR=$SMALL_BACKUP_DIR: a file under 64 MiB is never written onto /stor (AGENTS §7.3)" ;;
    esac
    case "$small_real" in
      /stor|/stor/*) refuse "SMALL_BACKUP_DIR=$SMALL_BACKUP_DIR ($small_real): a file under 64 MiB is never written onto /stor (AGENTS §7.3)" ;;
    esac
  fi
  if [ -s "$C/backup.tsv" ]; then  # written before: checked, never written again
    rec_tar=$(cut -f1 "$C/backup.tsv")
    sz=$(cut -f2 "$C/backup.tsv")
    want_tar=$BACKUP_DIR/$name
    if [ -n "$SMALL_BACKUP_DIR" ] && [ "$rec_tar" = "$SMALL_BACKUP_DIR/$name" ]; then
      want_tar=$SMALL_BACKUP_DIR/$name
    fi
    [[ "$sz" =~ ^[0-9]+$ ]] && [ "$rec_tar" = "$want_tar" ] && [ -f "$want_tar" ] \
      && [ "$(stat -c %s "$want_tar")" = "$sz" ] \
      || fail "$C/backup.tsv names a copy that is not at $want_tar with its size"
    log "backup: $want_tar was written before ($(tr '\t' ' ' < "$C/backup.tsv")); nothing written"
    mark_done backup
    return 0
  fi
  backup_real=$(realpath -m -- "$BACKUP_DIR") || refuse "could not resolve BACKUP_DIR=$BACKUP_DIR"
  case "$BACKUP_DIR" in
    /stor|/stor/*)
      [ "$TAR_FLOOR_BYTES" -ge "$STOR_FLOOR" ] || refuse "a file under 64 MiB is never written onto /stor (AGENTS §7.3)"
      [ -d "$BACKUP_DIR" ] || refuse "$BACKUP_DIR does not exist (nothing is made under /stor)" ;;
  esac
  case "$backup_real" in
    /stor|/stor/*)
      [ "$TAR_FLOOR_BYTES" -ge "$STOR_FLOOR" ] || refuse "a file under 64 MiB is never written onto /stor (AGENTS §7.3)"
      [ -d "$BACKUP_DIR" ] || refuse "$BACKUP_DIR does not exist (nothing is made under /stor)" ;;
  esac
  [ ! -e "$BACKUP_DIR/$name" ] && [ ! -L "$BACKUP_DIR/$name" ] \
    || refuse "$BACKUP_DIR/$name already exists; it is never overwritten"
  [ ! -e "$BACKUP_DIR/.$name.tmp" ] && [ ! -L "$BACKUP_DIR/.$name.tmp" ] \
    || refuse "$BACKUP_DIR/.$name.tmp exists (an earlier copy that did not finish; a person reads it)"
  if [ -n "$SMALL_BACKUP_DIR" ]; then
    [ ! -e "$SMALL_BACKUP_DIR/$name" ] && [ ! -L "$SMALL_BACKUP_DIR/$name" ] \
      || refuse "$SMALL_BACKUP_DIR/$name already exists; it is never overwritten"
    [ ! -e "$SMALL_BACKUP_DIR/.$name.tmp" ] && [ ! -L "$SMALL_BACKUP_DIR/.$name.tmp" ] \
      || refuse "$SMALL_BACKUP_DIR/.$name.tmp exists (an earlier copy that did not finish; a person reads it)"
  fi
  awk -F '\t' 'NR > 1 {print $1}' "$C/root_added.tsv" > "$C/added_paths.txt" || fail "could not read root_added.tsv"
  n=$(wc -l < "$C/added_paths.txt")
  [ "$n" -gt 0 ] || fail "root_added.tsv lists no file"
  rm -f "$local_tar"
  tar -C "$ROOT" -cf "$local_tar" -T "$C/added_paths.txt" || fail "tar of the added files failed"
  [ "$(tar -tf "$local_tar" | wc -l)" = "$n" ] || fail "the tar does not hold the $n added files"
  sz=$(stat -c %s "$local_tar")
  if [ "$sz" -ge "$TAR_FLOOR_BYTES" ]; then
    dst_dir=$BACKUP_DIR
    case "$backup_real" in
      /stor|/stor/*) ;;
      *)
        case "$BACKUP_DIR" in
          /stor|/stor/*) ;;
          *) mkdir -p "$dst_dir" || fail "could not make $dst_dir" ;;
        esac ;;
    esac
  elif [ -n "$SMALL_BACKUP_DIR" ]; then
    [ "$sz" -ge "$SMALL_TAR_FLOOR_BYTES" ] \
      || fail "$local_tar is $sz bytes, under the small-backup floor of $SMALL_TAR_FLOOR_BYTES"
    dst_dir=$SMALL_BACKUP_DIR
    mkdir -p "$dst_dir" || fail "could not make $dst_dir"
    case "$(realpath -e -- "$dst_dir")" in
      /stor|/stor/*) fail "SMALL_BACKUP_DIR=$dst_dir resolved under /stor (AGENTS §7.3)" ;;
    esac
  else
    fail "$local_tar is $sz bytes, under the floor of $TAR_FLOOR_BYTES"
  fi
  tar=$dst_dir/$name
  tmp=$dst_dir/.$name.tmp
  [ ! -e "$tar" ] && [ ! -L "$tar" ] || refuse "$tar already exists; it is never overwritten"
  [ ! -e "$tmp" ] && [ ! -L "$tmp" ] || refuse "$tmp exists (an earlier copy that did not finish; a person reads it)"
  sha1=$(sha256sum "$local_tar" | cut -d' ' -f1)
  [[ "$sha1" =~ ^[0-9a-f]{64}$ ]] || fail "the sha256 of $local_tar is not 64 hex digits: '$sha1'"
  log "backup: local tar $local_tar: $n files, $sz bytes, sha256 $sha1 -> $tar"
  dd if="$local_tar" of="$tmp" bs=16M conv=fsync 2>> "$LOG" || fail "the copy to $tmp failed (left for a person)"
  sha2=$(sha256sum "$tmp" | cut -d' ' -f1)
  [[ "$sha2" =~ ^[0-9a-f]{64}$ ]] || fail "the sha256 read back from $tmp is not 64 hex digits: '$sha2' (left for a person)"
  [ "$sha1" = "$sha2" ] || fail "read back sha256 $sha2, not $sha1 (left $tmp for a person)"
  mv -n "$tmp" "$tar" || fail "mv -n $tmp $tar failed"
  [ -f "$tar" ] && [ ! -e "$tmp" ] || fail "$tmp is still there after mv -n (did $tar appear meanwhile?)"
  printf "%s\t%s\t%s\t%s\n" "$tar" "$sz" "$sha2" "$n" > "$C/backup.tsv"
  log "backup: $tar written ($n files, $sz bytes, sha256 $sha2)"
  rm -f "$local_tar"
  mark_done backup
}

rehearse() {  # the fixture, in scratch: E and a release over it, then the stage through this driver's own modes; in a
  # wave, another stage's full rehearsal of the wave under the same facts is taken instead (rh_share, card A-30 (2))
  local rh=$C/rh rs="${STAGE}R" ml=() pbind=(--bind /tmp/if2453) m rc eargs args old rq cap n sh=0
  read_conf
  if [ -e "$C/rehearse/DONE" ] && [ "$(head -c 64 "$C/rehearse/DONE")" = "$CONF_SHA" ]; then
    log "rehearse: $C/rehearse/DONE present for this conf (remove it to rehearse again)"
    return 0
  fi
  # a rehearsal of this conf is new: the markers of the modes that rest on it are withdrawn
  for m in rehearse run after check backup; do rm -f "$C/$m/DONE"; done
  mkdir -p "$C/rehearse"
  rm -f "$C/rehearse/FAILED" "$C/rehearse/REFUSED" "$C/rehearse/STOPPED"
  caps
  [[ "$rs" =~ ^[A-Z][A-Za-z0-9]{0,15}$ ]] || refuse "the rehearsal's stage name $rs is too long"
  rm -f "$C/rehearse/facts" "$C/rehearse/facts.end" "$C/rehearse/DONE.tmp"
  if share_on; then
    sh=1
    rh_share && return 0
  fi
  quota_now || refuse "could not read the profile's cap of this hour (stage_est.py cap)"
  cap=$QUOTA_NOW
  rq=$(( RH_Q < cap ? RH_Q : cap ))
  # the earlier rehearsal's logs are kept, moved aside; its scratch (roots, inputs) is removed
  if [ -e "$rh" ] || compgen -G "$LOGS/$LC-rh-*.log" > /dev/null; then
    n=$(find "$C" -mindepth 1 -maxdepth 1 -type d -name 'rh-logs.*' | wc -l)
    old=$C/rh-logs.$((n + 1))
    mkdir -p "$old" || fail "could not make $old"
    [ -e "$rh/drv/driver.log" ] && mv "$rh/drv/driver.log" "$old/inner-driver.log"
    [ -d "$rh/drv/logs" ] && mv "$rh/drv/logs" "$old/inner-logs"
    [ -d "$rh/drv/check" ] && mv "$rh/drv/check" "$old/inner-check"
    [ -e "$rh/fixture.json" ] && mv "$rh/fixture.json" "$old/"
    for m in "$LOGS/$LC-rh-"*.log; do [ -e "$m" ] && mv "$m" "$old/"; done
    log "rehearse: the earlier rehearsal's logs moved to $old"
  fi
  rm -rf "$rh"  # this driver's own scratch of an earlier rehearsal (its logs kept above)
  mkdir -p "$rh/tmp" "$rh/durable" "$rh/bulk" "$rh/root" "$rh/backups" || fail "could not make $rh"
  sed -e "s#^tmp:.*#tmp: $rh/tmp#" -e "s#^durable:.*#durable: $rh/durable#" -e "s#^bulk:.*#bulk: $rh/bulk#" \
    -e "s#^small_file_floor_bytes:.*#small_file_floor_bytes: 1000#" "$REALPROF" > "$rh/profile.yaml" \
    || fail "could not write $rh/profile.yaml"
  # stage E is built with a budget of 1 CPU-minute, so that its steps run in pieces and its scan in shards with a plan:
  # only a scan with a shard plan is a source of copies (domains._source, card B-47b), as card C-26's rehearsal did
  sed -e "s#^cpu_min_per_process:.*#cpu_min_per_process: 1#" "$rh/profile.yaml" > "$rh/profile-e.yaml" \
    || fail "could not write $rh/profile-e.yaml"
  [ -n "$MINLEN" ] && ml=(--min-contig-len "$MINLEN")
  bb py "$BIN/stage_est.py" fixture --tree "$TREE" --out "$rh/inputs" --stage "$rs" --earlier E "${ml[@]}" \
    > "$rh/fixture.json" 2>> "$LOG" || fail "rehearse: the fixture's inputs (stage_est.py fixture)"
  log "rehearse: fixture: $(cat "$rh/fixture.json")"
  [ -n "$LIBNAME" ] && pbind+=(--bind "$PFAMDIR:$PFAMDIR:ro")
  eargs=(--stage E --release rRH.1 --stages E --depth "$DEPTH" --genomes "$rh/inputs/E/genomes.json"
         --objects "$rh/inputs/E/objects.json" "${LIBARGS[@]}")
  [ "$CRISPR" = yes ] && eargs+=(--crispr)
  rh_launch "$LC-rh-e" build-stage --profile "$rh/profile-e.yaml" --root "$rh/root" "${eargs[@]}" \
    --threads "$RH_T" --scan-jobs 1 --scan-cpus "$RH_T" --cluster-parts 2
  [ -s "$rh/root/manifests/release=rRH.1.json" ] || fail "rehearse: the scratch root has no release rRH.1"
  {
    echo "# the rehearsal's stage, made by stage_driver.sh rehearse from $CONF"
    echo "STAGE=$rs"
    for m in $NOTES; do echo "NOTE=$m"; done
    echo "INPUTS=$rh/inputs/R"
    for m in genomes.json objects.json stage.fna; do echo "SUM=$(sha256sum "$rh/inputs/R/$m" | cut -d' ' -f1) $m"; done
    echo "DEPTH=$DEPTH"
    if [ -n "$LIBNAME" ]; then echo "LIBRARY=$LIBNAME $LIBVER $LIBCUT"; fi
    echo "CRISPR=$CRISPR"
    if [ -n "$MINLEN" ]; then echo "MIN_CONTIG_LEN=$MINLEN"; fi
    if [ -n "$GENECALL_MODE" ]; then echo "GENECALL_MODE=$GENECALL_MODE"; fi
  } > "$rh/$rs.conf" || fail "could not write $rh/$rs.conf"
  for m in $CHAIN; do
    STAGE_ROOT="$rh/root" STAGE_DIR="$rh/drv" STAGE_CONF="$rh/$rs.conf" PROFILE="$rh/profile.yaml" \
      BACKUP_DIR="$rh/backups" TAR_FLOOR_BYTES=0 SMALL_BACKUP_DIR= STAGE_REHEARSAL=1 STAGE_UNIT_STEM="$LC-rh" \
      T="$RH_T" CPUS="$RH_T" SCAN_JOBS=1 Q="$rq" STAGE_WAVE= bash "$SELF" "$rs" "$m" >> "$LOG" 2>&1
    rc=$?
    [ "$rc" = 0 ] || fail "rehearse: mode $m of the rehearsal's stage exited $rc (read $rh/drv/driver.log)"
    [ -e "$rh/drv/$m/DONE" ] || fail "rehearse: mode $m exited 0 and left no $rh/drv/$m/DONE"
    log "rehearse: mode $m passed"
  done
  args=(--stage "$rs" --depth "$DEPTH")
  [ -n "$MINLEN" ] && args+=(--min-contig-len "$MINLEN")
  [ -n "$GENECALL_MODE" ] && args+=(--genecall-mode "$GENECALL_MODE")
  args+=(--genomes "$rh/inputs/R/genomes.json" --objects "$rh/inputs/R/objects.json" "${LIBARGS[@]}")
  [ "$CRISPR" = yes ] && args+=(--crispr)
  find "$rh/root" ! -type d ! -name STOP -printf '%P\t%s\t%T@\n' | LC_ALL=C sort > "$rh/again-before.tsv" \
    || fail "rehearse: listing"
  rh_launch "$LC-rh-again" build-stage --profile "$rh/drv/profile.yaml" --root "$rh/root" "${args[@]}" \
    --threads "$RH_T" --scan-jobs 1 --scan-cpus "$RH_T"
  grep -q ' published ' "$LOGS/$LC-rh-again.log" && fail "rehearse: the relaunch published a step"
  find "$rh/root" ! -type d ! -name STOP -printf '%P\t%s\t%T@\n' | LC_ALL=C sort > "$rh/again-after.tsv" \
    || fail "rehearse: listing"
  cmp -s "$rh/again-before.tsv" "$rh/again-after.tsv" || fail "rehearse: the relaunch changed the scratch root"
  log "rehearse: the relaunch published nothing and left the scratch root unchanged"
  if [ "$sh" = 1 ] && [ -s "$C/rehearse/facts" ] && rh_facts "$C/rehearse/facts.end" \
    && cmp -s "$C/rehearse/facts" "$C/rehearse/facts.end"; then
    { echo "$CONF_SHA"; echo "rehearsal=full"; echo "stage=$STAGE"; cat "$C/rehearse/facts"; } > "$C/rehearse/DONE.tmp" \
      && mv -f "$C/rehearse/DONE.tmp" "$C/rehearse/DONE" || fail "could not write $C/rehearse/DONE"
    log "rehearse: a later stage of wave $WAVE may take this rehearsal under the same facts: $(paste -sd ';' "$C/rehearse/facts")"
    return 0
  fi
  [ "$sh" = 1 ] && log "rehearse: the facts of this rehearsal were not read the same at its start and end: its marker serves this stage only"
  mark_done rehearse
}

rh_launch() {  # one build in the rehearsal's scratch root; a STOP in the data root stops it
  local unit=$1 rc
  shift
  RUNNING_UNIT=$unit
  RUNNING_LOGS=$LOGS
  bash "$TREE/scripts/pod_run.sh" --profile "$rh/profile.yaml" --image "$IMAGE" --digest "$DIGEST" --mem 16G \
    --threads "$rq" "${pbind[@]}" --env "PYTHONPATH=$PYPATH" --log "$LOGS/$unit.log" "$unit" "$@" >> "$LOG" 2>&1 \
    || { RUNNING_UNIT=""; fail "rehearse: pod_run.sh refused $unit"; }
  sleep 5
  while :; do
    unit_ended "$unit" "$LOGS"
    rc=$?
    [ "$rc" = 0 ] && break
    watch_answer "$rc" "$unit"
    [ -e "$REAL_ROOT/STOP" ] && stop_unit "$unit" "$LOGS" "rehearse: $REAL_ROOT/STOP appeared"
    sleep 5
  done
  RUNNING_UNIT=""
  rc=$(exit_of "$unit" "$LOGS")
  [ "$rc" = 0 ] || fail "rehearse: $unit exited '$rc' (read $LOGS/$unit.log)"
}

info() {
  read_conf
  echo "STAGE=$STAGE"
  echo "DIR=$C"
  echo "AFTER=$AFTER"
  echo "BUSY=$(build_units)"
}

if [ "$MODE" = info ]; then
  environment
  info
  exit $?
fi
log "stage driver: mode $MODE (pid $$), conf $CONF, root $ROOT, tree $TREE, image $IMAGE"
stop_file
environment
# The mode's own status is the driver's exit status (card A-11; STATUS "6m": on 2026-11-20 a driver's `check` returned
# 1 and the driver exited 0, the status of a log line after it). A mode that fails through `fail` exits 1 itself.
"$MODE"
rc=$?
if [ "$rc" = 0 ]; then
  log "stage driver: mode $MODE done"
else
  log "stage driver: mode $MODE FAILED with status $rc"
fi
exit "$rc"
