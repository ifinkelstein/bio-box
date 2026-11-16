#!/bin/bash
# Card C-26 part 2 (amended 2026-11-15): stage PL1, the PLSDB plasmids, built into the data root, the stage only, with
# Pfam 38.2 and both CRISPR producers, under the release name rM.3 (S1, P1 to P6, H1, PL1; it has no release_def and
# none is written). Host side, bash. Copies of this script and c26_est.py, c26_check.py, c26_view.py, c26_sample.py,
# make_inputs.py, r50_watch.py run from /tmp/if2453/c26/bin, so the synced tree is not touched while they run. Each mode
# is launched in its own small unit, for example
#
#   systemd-run --user --unit=c26-driver-run --collect -p MemoryMax=4G -p MemorySwapMax=0 \
#     -p StandardOutput=append:/tmp/if2453/c26/driver.out -p StandardError=append:/tmp/if2453/c26/driver.out \
#     nice -n 19 bash /tmp/if2453/c26/bin/c26_driver.sh run
#
# The build is card C-24's (c24_driver.sh) for stage PL1, with three differences:
# - Gene calling is over the margin, so build-stage runs it in batches, each a run-step with its record (cards B-48,
#   B-50): genecall.plan, genecall.batch-<j>, then the merge `genecall`. The rate check judges each batch against the
#   estimate in its published plan.
# - The Pfam scan's plan copies the rows of the proteins that S1, P1 to P6 and H1 measured (card B-47b) and scans the
#   rest in shards.
# - Once PL1 is gene-called every stage of rM.3 is, so the cluster gate no longer stops build-stage; the launch ends by
#   a STOP this driver writes as soon as the stage's last step (the filtered CRISPR producer) has its intent file or
#   run record (part 1, planning/measurements/2026-11-15_plsdb_inputs.md §5): run-step writes the intent file only
#   after its own post-command STOP check, so the step's work is kept, and build-stage, or run-step before the release's
#   first step, then exits 75 with no file of the release written. A relaunch whose stage steps are all recorded
#   writes the STOP when the log says the last step was skipped. The STOP holds a token line this driver records in
#   C/stop_tokens; the driver removes a STOP only when its content is one of those lines. Any other STOP in the root
#   at any time: the driver writes nothing, waits for the unit to end (build-stage stops itself at the next step
#   boundary), and fails, so that it is reported. Second guard: if a path of the release (runs/stage=rM.3,
#   release_def/release=rM.3, l2/cluster_plan/release=rM.3) appears, or a `biobox.cli release` process runs in the
#   build's unit, the unit is stopped at once (run-step forwards the SIGTERM and publishes nothing).
#
# Modes, in the order the card runs them:
#   rehearse  In a scratch root C/rh/root (a view of the data root: P1 and P2 linked with their Pfam 38.2 scans; the
#             data root bound read-only in every command): part 1's sample of 315 PLSDB records (c26_sample.py, the
#             same records), its inputs (make_inputs.py plasmid), and the real command's form for the sample stage PLR
#             under release rC26.3 over P1,P2,PLR, with --library pfam 38.2 <the real Pfam file> gathering --crispr, a
#             scratch profile whose cpu_min_per_process is RH_BUDGET (default 1) so that gene calling runs in batches
#             (checked before the launch with c26_est.py genecall) and the scan in shards; with the STOP watcher. It
#             must end with exit 75 at rC26.3's cluster step with every PLR stage step recorded, the gene-call plan and
#             batch records, a scan plan that copies rows from P1 or P2, both CRISPR producers, nothing of rC26.3, and
#             c26_check.py passing on it (the copied rows equal the sources'); the same launch again must skip every
#             step, stop the same way and change no file. Writes C/rh/OK.
#   estimate  The inputs' sha256 against part 1's record, then c26_est.py before on the real inputs: ingest, both CRISPR
#             steps, the gene-call plan (batches each under 80% of the budget), part 1's expectations, into
#             C/estimate/before.json.
#   before    Refuses when ROOT/STOP exists or a build unit is active; then `release verify` of every published
#             release (root bound read-only) into C/before/verify-<R>.log, and the snapshot into C/before: listing.tsv
#             (path, bytes, mtime of every file of ROOT), earlier.sha256 (every file under a stage=<S> directory of S1,
#             P1 to P6 and H1, the stages PL1 shares contigs with or copies rows from). Writes C/before/DONE.
#   run       Needs C/before/DONE, C/rh/OK and C/estimate/before.json. Before the launch: no STOP (one of ours from an
#             earlier launch is removed), no build unit, at least START_MIN_GB free on /tmp, the allowance of a
#             T-thread launch at least T. One unit c26-pl1-<n> (pod_run.sh --threads Q, the profile's 200G),
#             build-stage --threads T --scan-jobs J --scan-cpus CPUS, and the STOP watcher. While it runs, once a
#             minute, r50_watch.py (units c26-pl1-*, run records of PL1 since the launch): a step over CPU_MAX
#             CPU-minutes, /tmp under TMP_MIN_GB free or scratch over 100 GB stops the unit; and each time PL1 has a
#             new run record, c26_est.py rates --expect: a step estimated at 24 CPU-minutes or more off its estimate
#             by more than 20%, or genes, distinct proteins or the scan's residues more than 20% outside part 1's
#             range, stops the unit the same way (AGENTS §6, DEC-43); when the only such steps are gene-call
#             batches, which run at once and end together, the other batches are let publish and this driver's STOP
#             then ends the launch before the merge (no batch's work is thrown away). A stop removes this driver's
#             STOP if it wrote one and fails; the same command resumes. Writes C/run/DONE when the unit ends as
#             expected.
#   after     The snapshot into C/after, `release verify` of every published release again, `bio-box releases`.
#   check     c26_check.py: before against after, root_added.tsv (copied to C/root_added.tsv), the steps against their
#             estimates, the data the card asks for, the copied rows against their sources, into C/check.
set -u

C=/tmp/if2453/c26
ROOT=/tmp/if2453/r10/root
TREE=${TREE:-/tmp/if2453/ci/C-26-2}
REALPROF=$TREE/profiles/wilkcomp02.yaml
PROFILE=$C/profile.yaml
IMAGE=/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif
DIGEST=sha256:d67300aa1321b571fa8da79e74dfccdd87a8dca2e56518f324902efa6b9f098a
PYPATH=$TREE/src
PFAMDIR=/stor/scratch/Wilke/if2453/all-the-data/l0/refdb/pfam/38.2
PFAM=$PFAMDIR/Pfam-A.hmm.gz
IN=/tmp/if2453/s2/inputs/PL1
META=/tmp/if2453/s2/plsdb/meta
STAGE=PL1
REL=rM.3
ALL=S1,P1,P2,P3,P4,P5,P6,H1,PL1
EARLIER=S1,P1,P2,P3,P4,P5,P6,H1
RELEASES="rN.1a rN.1b rN.1 rN.2 rN.3 rM.2"
T=${T:-36}
J=${J:-9}
CPUS=${CPUS:-4}
Q=${Q:-40}
CPU_MAX=${CPU_MAX:-280}
TMP_MIN_GB=${TMP_MIN_GB:-200}
START_MIN_GB=${START_MIN_GB:-205}
POLL=${POLL:-0.5}
RH=$C/rh
RHSTAGE=PLR
RHREL=rC26.3
RHALL=P1,P2,PLR
RH_BUDGET=${RH_BUDGET:-1}
RH_T=${RH_T:-8}
RH_J=${RH_J:-2}
BIN=$(cd "$(dirname "$0")" && pwd)
LOGS=$C/logs
LOG=$C/driver.log
SCRATCH=$C/tmp/scratch
NR3=minced-0.4.2-w8-nr3-r23x47-s26x50
NEAR=minced-0.4.2-w8-nr2-r23x47-s26x50-cas1kb-caspfam1-pfam-38.2
STEPS_LAST="crispr-$NEAR"
# sha256 of the inputs, from part 1's record (planning/measurements/2026-11-15_plsdb_inputs.md §2)
SUMS_IN="69ec5621edcef675398d1c3cbd07989ccdd80b96a3dd94110ab387f4494229e9  genomes.json
c70d2ec6470394707f976c5f5aef86f05c8c04f8e0e6d86934ff84707b51150a  objects.json
3d04f91367f2d582dff9837c5ba98c2a457844e293e03512d8fe5880d20396bb  repeated_records.tsv
4567837327e4709e3189e99bbbd97ad2a6ed034ca0b9d28bf6217d37589ef297  shared_records.tsv
e654c106d8a15a2d5d55a8b680087377cb84e654e63d41bca181fc69a23e587a  plsdb_2024_05_31_v2_sequences.fasta"
mkdir -p "$LOGS" "$C/tmp"
touch "$C/stop_tokens"

log() { echo "$(date '+%F %T.%3N %Z') $*" >> "$LOG"; }
fail() { log "FAILED: $*"; mkdir -p "$C/$MODE"; echo "$*" > "$C/$MODE/FAILED"; exit 1; }

profile() {  # the real profile with tmp under C (one line differs)
  sed "s#^tmp:.*#tmp: $C/tmp#" "$REALPROF" > "$PROFILE"
  log "profile $PROFILE: $(diff "$REALPROF" "$PROFILE" | grep -c '^>') line differs from $REALPROF"
}

bb() {  # a short command in the image, in its own capped scope at nice 19, the root, Pfam and PLSDB's metadata read-only
  local what=$1
  shift
  case "$what" in
    cli) set -- python -m biobox.cli "$@" ;;
    py) set -- python "$@" ;;
  esac
  systemd-run --user --scope --quiet -p MemoryMax=16G -p MemorySwapMax=0 -p CPUQuota=200% nice -n 19 \
    singularity exec --cleanenv --env PYTHONPATH="$PYPATH" --env OMP_NUM_THREADS=2 --env PYTHONUNBUFFERED=1 \
    --bind /tmp/if2453 --bind "$ROOT:$ROOT:ro" --bind "$PFAMDIR:$PFAMDIR:ro" --bind "$META:$META:ro" \
    --bind /sys/fs/cgroup:/sys/fs/cgroup:ro "$IMAGE" "$@"
}

build_units() {  # active units of builds that write into ROOT (this card's driver units are c26-driver-*)
  systemctl --user list-units --type=service --state=active --no-legend --plain \
    'r50-*' 'c24-h1-*' 'c25-*' 'c26-pl1-*' 'c26-rh-*' 'c33-*' 'c34-*' 'c35-*' 'b47-*' | awk '{print $1}' | paste -sd' ' -
}

snap_top() {
  log "uptime: $(uptime)"
  log "account CPU now (top, %CPU of if2453): $(top -b -n 1 -u if2453 | awk 'NR > 7 {s += $9} END {print s}')"
  top -b -n 1 -o %CPU -w 150 | sed -n '1,5p;7,22p' >> "$LOG"
  log "df /tmp: $(df -h /tmp | tail -n 1)"
  log "active units: $(systemctl --user list-units --type=service --type=scope --state=active --no-legend --plain | awk '{print $1}' | paste -sd' ' -)"
}

free_gb() { df -B1G --output=avail /tmp | tail -n 1 | tr -d ' '; }

exit_of() { sed -n "s/^pod_run: $1 exited \([0-9]*\)$/\1/p" "$LOGS/$1.log" | tail -n 1; }

ours_stop() {  # $1 root: its STOP holds one of this driver's token lines
  [ -e "$1/STOP" ] && [ -s "$C/stop_tokens" ] && grep -qxF -f "$C/stop_tokens" "$1/STOP"
}

remove_our_stop() {  # $1 root
  if ours_stop "$1"; then
    rm "$1/STOP"
    log "removed $1/STOP, written by this driver ($(cat "$C/stop_tokens" | tail -n 1))"
  elif [ -e "$1/STOP" ]; then
    log "$1/STOP is NOT this driver's; left as it is: $(head -c 300 "$1/STOP")"
    return 1
  fi
  return 0
}

verify() {  # $1 root, $2 release, $3 directory
  bb cli release verify --root "$1" --release "$2" > "$3/verify-$2.log" 2>&1
  local rc=$?
  log "verify $2 ($1): exit $rc: $(grep -E '^release ' "$3/verify-$2.log" | tail -n 1)"
  return $rc
}

snapshot() {  # $1 directory
  local d=$1
  mkdir -p "$d"
  find "$ROOT" ! -type d -printf '%P\t%s\t%T@\n' | LC_ALL=C sort > "$d/listing.tsv"
  find "$ROOT" ! -type d -regextype posix-extended -regex ".*/stage=(${EARLIER//,/|})/.*" -printf '%P\n' \
    | LC_ALL=C sort > "$d/earlier.files"
  (cd "$ROOT" && systemd-run --user --scope --quiet -p MemoryMax=2G -p MemorySwapMax=0 nice -n 19 ionice -c 3 \
    xargs -d '\n' sha256sum < "$d/earlier.files") > "$d/earlier.sha256"
  log "snapshot $d: $(wc -l < "$d/listing.tsv") files, $(wc -l < "$d/earlier.sha256") of $EARLIER hashed, $(awk -F '\t' '{s += $2} END {printf "%.0f", s}' "$d/listing.tsv") bytes"
}

# --- the launch and the STOP watcher (rehearsal and real build) -----------------------------------------------------

release_guard() {  # $1 unit, $2 state dir, $3 why
  [ -e "$2/release_guard" ] && return
  echo "$3" > "$2/release_guard"
  log "SECOND GUARD: $3; stopping $1 at once"
  systemctl --user stop "$1"
}

watch_stop() {  # $1 root, $2 stage, $3 release, $4 unit, $5 unit log, $6 state dir (holds `since`)
  local root=$1 stage=$2 rel=$3 unit=$4 ulog=$5 st=$6 tok i=0 pid d
  local rec="$root/runs/stage=$stage/$STEPS_LAST.json" intent="$root/runs/stage=$stage/$STEPS_LAST.publishing.json"
  while systemctl --user is-active --quiet "$unit"; do
    if [ ! -e "$st/stop_written" ] && [ ! -e "$st/foreign_stop" ]; then
      if [ -e "$intent" ] || { [ -e "$rec" ] && [ "$rec" -nt "$st/since" ]; } \
         || grep -q "$stage\.$STEPS_LAST skipped: run record and outputs present" "$ulog" 2> /dev/null; then
        tok="c26 driver STOP $unit $(date '+%F %T.%N %Z') pid $$"
        echo "$tok" >> "$C/stop_tokens"
        if ( set -C; echo "$tok" > "$root/STOP" ) 2> /dev/null; then
          echo "$tok" > "$st/stop_written"
          log "watcher: the last stage step's intent file, run record or skip line appeared; wrote $root/STOP ($tok)"
        else
          touch "$st/foreign_stop"
          log "watcher: $root/STOP exists and this watcher did not write it: $(head -c 300 "$root/STOP")"
        fi
      fi
    fi
    if [ -e "$root/STOP" ] && [ ! -e "$st/foreign_stop" ] && ! ours_stop "$root"; then
      touch "$st/foreign_stop"
      log "watcher: a STOP not written by this driver is in $root: $(head -c 300 "$root/STOP"); writing nothing, waiting for $unit to end"
    fi
    i=$((i + 1))
    if [ $((i % 4)) = 0 ]; then
      for d in "runs/stage=$rel" "release_def/release=$rel" "l2/cluster_plan/release=$rel"; do
        [ -e "$root/$d" ] && release_guard "$unit" "$st" "path $root/$d exists"
      done
      for pid in $(pgrep -u "$(id -u)" -f 'biobox[.]cli release '); do
        grep -q "/$unit.service" "/proc/$pid/cgroup" 2> /dev/null \
          && release_guard "$unit" "$st" "a release command runs in $unit: pid $pid: $(tr '\0' ' ' < "/proc/$pid/cmdline" | cut -c1-300)"
      done
    fi
    sleep "$POLL"
  done
}

launch() {  # $1 root, $2 profile, $3 unit, $4 mem, $5 quota, $6 state dir, rest: build-stage arguments after --root
  local root=$1 prof=$2 unit=$3 mem=$4 quota=$5 st=$6 ro=()
  shift 6
  [ "$root" = "$ROOT" ] || ro=(--bind "$ROOT:$ROOT:ro")
  mkdir -p "$st"
  rm -f "$st/stop_written" "$st/foreign_stop" "$st/release_guard"
  touch "$st/since"
  log "launch $unit (pod_run --threads $quota --mem $mem): build-stage --profile $prof --root $root $*"
  bash "$TREE/scripts/pod_run.sh" --profile "$prof" --image "$IMAGE" --digest "$DIGEST" --mem "$mem" \
    --threads "$quota" --bind /tmp/if2453 "${ro[@]}" --bind "$PFAMDIR:$PFAMDIR:ro" --env "PYTHONPATH=$PYPATH" \
    --log "$LOGS/$unit.log" "$unit" build-stage --profile "$prof" --root "$root" "$@" >> "$LOG" 2>&1
  sleep 1
  watch_stop "$root" "$(arg --stage "$@")" "$(arg --release "$@")" "$unit" "$LOGS/$unit.log" "$st" &
  WATCHER=$!
}

arg() {  # $1 option, rest: arguments: the option's value
  local o=$1
  shift
  while [ $# -gt 1 ]; do
    [ "$1" = "$o" ] && { echo "$2"; return; }
    shift
  done
}

ended() {  # $1 root, $2 unit, $3 stage, $4 release, $5 state dir: how the launch ended; 0 when as expected
  local root=$1 unit=$2 stage=$3 rel=$4 st=$5 rc n step bad=0
  wait "$WATCHER" 2> /dev/null
  rc=$(exit_of "$unit")
  log "$unit ended, exit '$rc'; $(grep -E 'STOP|refused' "$LOGS/$unit.log" | tail -n 3 | tr '\n' ' ' | cut -c1-700)"
  log "$unit usage: $(grep -E 'pod_run: .* usage' "$LOGS/$unit.log" | tail -n 1)"
  [ -e "$st/release_guard" ] && { log "the second guard fired: $(cat "$st/release_guard")"; bad=1; }
  if [ -e "$st/foreign_stop" ]; then
    log "a STOP not written by this driver appeared in $root during $unit: it is left in place"
    return 3
  fi
  [ -e "$st/stop_written" ] || { log "the watcher wrote no STOP"; bad=1; }
  [ "$rc" = 75 ] || { log "exit '$rc', not 75"; bad=1; }
  grep -qE "(STOP exists; (stopping before|not starting)|STOP appeared; not publishing) $rel\.cluster" "$LOGS/$unit.log" \
    || { log "no line saying the launch stopped at $rel's cluster step"; bad=1; }
  grep -qE "STOP appeared; not publishing $rel\.cluster" "$LOGS/$unit.log" \
    && log "NOTE: $rel's cluster plan step ran its command before the STOP; it published nothing"
  grep -qE "run-step: $rel\.[^ ]* published" "$LOGS/$unit.log" && { log "a release step published"; bad=1; }
  for step in inputs ingest genecall genecall.plan "crispr-$NR3" igr tapes domains-pfam-38.2 domains-pfam-38.2.plan \
              "crispr-$NEAR"; do
    [ -e "$root/runs/stage=$stage/$step.json" ] || { log "no run record $stage/$step"; bad=1; }
  done
  n=$(find "$root/runs/stage=$stage" -maxdepth 1 -name 'genecall.batch-*.json' ! -name '*.publishing.json' | wc -l)
  log "$stage: $n gene-call batch records; $(find "$root/runs/stage=$stage" -maxdepth 1 -name '*.shard-*.json' | wc -l) scan shard records"
  [ "$n" -ge 2 ] || { log "fewer than 2 gene-call batch records"; bad=1; }
  [ -e "$root/l1/genecall_plan/stage=$stage/plan.json" ] || { log "no l1/genecall_plan/stage=$stage/plan.json"; bad=1; }
  find "$root/runs/stage=$stage" -maxdepth 1 -name '*.publishing.json' | grep -q . && { log "an intent file is left"; bad=1; }
  n=$(find "$root" -path "*release=$rel*" -o -path "*stage=$rel*" -o -name "*release=$rel*" | wc -l)
  log "paths naming $rel in $root: $n; in $C/tmp and $RH/tmp: $(find "$C/tmp" "$RH/tmp" -name "*$rel*" 2> /dev/null | wc -l)"
  [ "$n" = 0 ] || bad=1
  remove_our_stop "$root" || bad=1
  return $bad
}

# --- rehearse -------------------------------------------------------------------------------------------------------

rh_listing() {  # $1 file: path, size and mtime of every entry of the view (links listed, not followed), STOP left out
  find "$RH/root" ! -type d ! -name STOP -printf '%P\t%s\t%T@\n' | LC_ALL=C sort > "$1"
}

rehearse() {
  local rc st=$RH/run args n unit
  [ -e "$RH/OK" ] && { log "rehearsal: OK present; nothing to do"; return; }
  mkdir -p "$RH/inputs" "$RH/tmp" "$RH/durable" "$RH/bulk" "$RH/check" "$st"
  if [ ! -e "$RH/root" ]; then
    python3 -I "$BIN/c26_view.py" "$ROOT" "$RH/root" P1,P2 pfam-38.2 >> "$LOG" 2>&1 || fail "rehearsal: view"
  fi
  if [ ! -e "$RH/inputs/genomes.json" ]; then
    python3 -I "$BIN/c26_sample.py" "$IN" "$META/nuccore.csv" "$IN/plsdb_2024_05_31_v2_sequences.fasta" \
      "$RH/inputs/sample.fasta" "$RH/inputs/sample.tsv" P1,P2 >> "$LOG" 2>&1 || fail "rehearsal: sample"
    cmp -s "$RH/inputs/sample.tsv" "$BIN/part1_sample.tsv" && log "rehearsal: the sample is part 1's, record for record" \
      || fail "rehearsal: the sample differs from part 1's sample.tsv"
    bb py "$BIN/make_inputs.py" plasmid "$RH/inputs/sample.fasta" "$META/nuccore.csv" "$RH/inputs" --root "$RH/root" \
      > "$RH/inputs/report.json" 2>> "$LOG" || fail "rehearsal: make_inputs.py"
    log "rehearsal: sample inputs: $(tr -d '\n' < "$RH/inputs/report.json" | cut -c1-600)"
  fi
  # the real caps (40 by day, 56 at night): the launch takes --threads RH_T, and the allowance subtracts what else of
  # ours runs (a cap of RH_T itself held the launch while another agent's tests used 6 CPUs)
  cat > "$RH/profile.yaml" << EOF
max_threads_day: 40
max_threads_night: 56
mem_scope: 32G
tmp: $RH/tmp
durable: $RH/durable
bulk: $RH/bulk
nice_monolithic: 19
nice_fanout: 15
day_start_hour: 8
day_end_hour: 20
cpu_min_per_process: $RH_BUDGET
cluster_cpu_min_max: 20000
small_file_floor_bytes: 1000
require_mem_cap: true
EOF
  bb py "$BIN/c26_est.py" genecall --root "$RH/root" --genomes "$RH/inputs/genomes.json" \
    --fasta "$RH/inputs/sample.fasta" --stage "$RHSTAGE" --budget "$RH_BUDGET" --out "$RH/check/genecall_plan.json" \
    > "$RH/check/genecall_estimate.json" 2>> "$LOG" || fail "rehearsal: c26_est.py genecall"
  n=$(python3 -I -c 'import json,sys; print(len(json.load(open(sys.argv[1]))["batches"]))' "$RH/check/genecall_estimate.json")
  log "rehearsal: gene-call plan of $RHSTAGE at a budget of $RH_BUDGET: $n batches; $(tr -d '\n ' < "$RH/check/genecall_estimate.json" | cut -c1-500)"
  [ "$n" -ge 2 ] || fail "rehearsal: the gene-call plan at a budget of $RH_BUDGET has $n batch(es): lower RH_BUDGET or enlarge the sample"
  args=(--stage "$RHSTAGE" --release "$RHREL" --stages "$RHALL" --genomes "$RH/inputs/genomes.json"
        --objects "$RH/inputs/objects.json" --library pfam 38.2 "$PFAM" gathering --crispr --threads "$RH_T"
        --scan-jobs "$RH_J" --scan-cpus 4)
  # first launch
  unit=c26-rh-1
  if [ ! -e "$st/OK-$unit" ]; then
    remove_our_stop "$RH/root" || fail "rehearsal: a STOP not ours in $RH/root"
    rh_listing "$RH/check/before-listing.tsv"
    launch "$RH/root" "$RH/profile.yaml" "$unit" 32G "$RH_T" "$st/$unit" "${args[@]}"
    while systemctl --user is-active --quiet "$unit"; do sleep 5; done
    ended "$RH/root" "$unit" "$RHSTAGE" "$RHREL" "$st/$unit"
    rc=$?
    [ "$rc" = 0 ] || fail "rehearsal: $unit did not end as expected (code $rc; read $LOGS/$unit.log)"
    touch "$st/OK-$unit"
  fi
  # the scan's plan copies rows from P1 or P2
  python3 -I -c 'import json,sys; c=json.load(open(sys.argv[1]))["copy"]; print(json.dumps({k: c[k] for k in ("proteins","residues","rows","proteins_with_hit")} | {"sources": [s["stage"] for s in c["sources"]], "passed_over": c["passed_over"]})); sys.exit(0 if c["proteins"] > 0 and c["sources"] else 1)' \
    "$RH/root/l1/domain_scan/library=pfam-38.2/stage=$RHSTAGE/plan.json" > "$RH/check/copy.json" \
    || fail "rehearsal: the scan's plan copies nothing: $(cat "$RH/check/copy.json")"
  log "rehearsal: the scan's plan copies: $(cat "$RH/check/copy.json")"
  # the checks of the real build on the rehearsal (c26_check.py itself is rehearsed)
  mkdir -p "$RH/check/before" "$RH/check/after"
  cp "$RH/check/before-listing.tsv" "$RH/check/before/listing.tsv"
  rh_listing "$RH/check/after/listing.tsv"
  bb py "$BIN/c26_check.py" --root "$RH/root" --before "$RH/check/before" --after "$RH/check/after" \
    --log "$LOGS/c26-rh-1.log" --objects "$RH/inputs/objects.json" --out "$RH/check/out" --tmp "$RH/tmp/check" \
    --stage "$RHSTAGE" --release "$RHREL" --earlier P1,P2 --expect-owned "$RH/inputs/shared_records.tsv" \
    > "$RH/check/check.out" 2>&1
  rc=$?
  log "rehearsal: c26_check.py exit $rc: $(tail -n 1 "$RH/check/check.out" | cut -c1-600)"
  [ "$rc" = 0 ] || fail "rehearsal: c26_check.py failed (read $RH/check/check.out)"
  # the same launch again: every step skipped, the same end, no file changed
  unit=c26-rh-2
  if [ ! -e "$st/OK-$unit" ]; then
    remove_our_stop "$RH/root" || fail "rehearsal: a STOP not ours in $RH/root"
    rh_listing "$RH/check/listing-again-before.tsv"
    launch "$RH/root" "$RH/profile.yaml" "$unit" 32G "$RH_T" "$st/$unit" "${args[@]}"
    while systemctl --user is-active --quiet "$unit"; do sleep 5; done
    ended "$RH/root" "$unit" "$RHSTAGE" "$RHREL" "$st/$unit"
    rc=$?
    [ "$rc" = 0 ] || fail "rehearsal: $unit did not end as expected (code $rc; read $LOGS/$unit.log)"
    n=$(grep -c 'skipped: run record and outputs present' "$LOGS/$unit.log")
    log "rehearsal: $unit: $n steps skipped; published lines: $(grep -c ' published ' "$LOGS/$unit.log")"
    grep -q ' published ' "$LOGS/$unit.log" && fail "rehearsal: the second launch published a step"
    grep -q "run-step: $RHSTAGE\.genecall skipped" "$LOGS/$unit.log" || fail "rehearsal: the second launch did not skip gene calling"
    rh_listing "$RH/check/listing-again-after.tsv"
    cmp -s "$RH/check/listing-again-before.tsv" "$RH/check/listing-again-after.tsv" \
      || fail "rehearsal: the second launch changed the view"
    log "rehearsal: $unit left the view's listing unchanged"
    touch "$st/OK-$unit"
  fi
  touch "$RH/OK"
  log "rehearsal: OK"
}

# --- estimate, before, after ----------------------------------------------------------------------------------------

estimate() {
  local d=$C/estimate rc
  mkdir -p "$d"
  profile
  echo "$SUMS_IN" > "$d/SUMS.expected"
  (cd "$IN" && systemd-run --user --scope --quiet -p MemoryMax=2G -p MemorySwapMax=0 nice -n 19 ionice -c 3 \
    sha256sum -c "$d/SUMS.expected") > "$d/SUMS.check" 2>&1
  rc=$?
  log "estimate: the inputs' sha256 against part 1's: exit $rc: $(tr '\n' ' ' < "$d/SUMS.check")"
  [ "$rc" = 0 ] || fail "estimate: the inputs differ from part 1's record"
  bb py "$BIN/c26_est.py" before --root "$ROOT" --profile "$PROFILE" --inputs "$IN" --stage "$STAGE" \
    --out-plan "$d/genecall_plan.json" > "$d/before.json" 2> "$d/before.err"
  rc=$?
  log "estimate: exit $rc: $(tr -d '\n ' < "$d/before.json" | cut -c1-1500)"
  [ "$rc" = 0 ] || fail "estimate: c26_est.py before (read $d/before.err)"
  python3 -I -c 'import json,sys; d=json.load(open(sys.argv[1])); g=d["genecall"]; m=d["margin_cpu_min"]; bad=[b for b in g["batches"] if b["estimate_cpu_min"] > m]; print(len(g["batches"]), "batches, largest", max(b["estimate_cpu_min"] for b in g["batches"]), "margin", m, "refusal", g["refusal"]); sys.exit(1 if bad or g["refusal"] or len(g["batches"]) < 2 else 0)' \
    "$d/before.json" >> "$LOG" 2>&1 || fail "estimate: a gene-call batch over the margin, a refusal, or no batches"
}

before() {
  local d=$C/before r ok=1 act
  [ -e "$d/DONE" ] && { log "before: DONE present; nothing to do"; return; }
  [ -e "$ROOT/STOP" ] && fail "$ROOT/STOP exists: $(head -c 300 "$ROOT/STOP")"
  act=$(build_units)
  [ -z "$act" ] || fail "a build unit is active: $act"
  mkdir -p "$d"
  snap_top
  for r in $RELEASES; do verify "$ROOT" "$r" "$d" || ok=0; done
  [ "$ok" = 1 ] || fail "a published release does not verify before the first write (read $d/verify-*.log)"
  snapshot "$d"
  touch "$d/DONE"
}

after() {
  local d=$C/after r
  mkdir -p "$d"
  [ -e "$ROOT/STOP" ] && log "after: $ROOT/STOP exists: $(head -c 300 "$ROOT/STOP")"
  snapshot "$d"
  for r in $RELEASES; do verify "$ROOT" "$r" "$d" || log "after: $r does NOT verify"; done
  bb cli releases --root "$ROOT" > "$d/releases.txt" 2>&1
  log "after: bio-box releases exit $?: $(wc -l < "$d/releases.txt") lines"
  touch "$d/DONE"
}

# --- run ------------------------------------------------------------------------------------------------------------

next_unit() {
  local n
  n=$(($( [ -e "$C/units" ] && wc -l < "$C/units" || echo 0) + 1))
  echo "c26-pl1-$n"
}

stage_records() { find "$ROOT/runs/stage=$STAGE" -maxdepth 1 -name '*.json' ! -name '*.publishing.json' 2> /dev/null | wc -l; }

stop_unit() {  # $1 unit, $2 reason
  log "stopping $1: $2"
  systemctl --user stop "$1"
  while systemctl --user is-active --quiet "$1"; do sleep 2; done
  wait "$WATCHER" 2> /dev/null
  remove_our_stop "$ROOT"
  fail "$2 (unit $1 stopped; the root is left as the steps left it; the same command resumes)"
}

drain_stop() {  # $1 unit, $2 reason: a gene-call batch off its estimate. The batches run at once and end together, so
  # killing the unit would throw away the others' work: let every batch publish (or end), then write this driver's STOP,
  # which ends the launch before the merge (a step boundary, RUNBOOK §3), and fail
  local unit=$1 why=$2 nb have run pid tok
  log "rate miss on a gene-call batch: $why; letting the other batches publish, then a STOP before the merge"
  nb=$(python3 -I -c 'import json,sys; print(len(json.load(open(sys.argv[1]))["batches"]))' \
    "$ROOT/l1/genecall_plan/stage=$STAGE/plan.json")
  while systemctl --user is-active --quiet "$unit"; do
    have=$(find "$ROOT/runs/stage=$STAGE" -maxdepth 1 -name 'genecall.batch-*.json' ! -name '*.publishing.json' | wc -l)
    run=""
    for pid in $(pgrep -u "$(id -u)" -f 'biobox[.]cli stage genecall batch'); do
      grep -q "/$unit.service" "/proc/$pid/cgroup" 2> /dev/null && run="$run $pid"
    done
    if [ "$have" -ge "$nb" ] || [ -z "$run" ]; then
      tok="c26 driver STOP $unit $(date '+%F %T.%N %Z') pid $$ (rate miss)"
      echo "$tok" >> "$C/stop_tokens"
      ( set -C; echo "$tok" > "$ROOT/STOP" ) 2> /dev/null && log "wrote $ROOT/STOP ($have of $nb batches recorded)" \
        || log "$ROOT/STOP exists already: $(head -c 300 "$ROOT/STOP")"
      break
    fi
    sleep 2
  done
  while systemctl --user is-active --quiet "$unit"; do sleep 2; done
  wait "$WATCHER" 2> /dev/null
  log "$unit ended, exit '$(exit_of "$unit")': $(grep -E 'STOP' "$LOGS/$unit.log" | tail -n 2 | tr '\n' ' ' | cut -c1-500)"
  remove_our_stop "$ROOT"
  fail "$why (the launch ended at a step boundary after the batches; the same command resumes)"
}

run() {
  local unit a free act rc line nrec=-1 n k=0
  [ -e "$C/before/DONE" ] || fail "no $C/before/DONE: run the before mode first"
  [ -e "$RH/OK" ] || fail "no $RH/OK: the rehearsal has not passed"
  [ -s "$C/estimate/before.json" ] || fail "no $C/estimate/before.json: run the estimate mode first"
  [ -e "$C/run/DONE" ] && { log "run: DONE present; nothing to do"; return; }
  mkdir -p "$C/run"
  profile
  remove_our_stop "$ROOT" || fail "$ROOT/STOP exists and is not this driver's"
  act=$(build_units)
  [ -z "$act" ] || fail "a build unit is active: $act"
  free=$(free_gb)
  [ "$free" -ge "$START_MIN_GB" ] || fail "/tmp has $free GB free, under $START_MIN_GB: not starting"
  a=$(bb py -m biobox.runstep --profile "$PROFILE" --threads "$T" 2>> "$LOG")
  [[ "$a" =~ ^[0-9]+$ ]] && [ "$a" -ge "$T" ] || fail "the allowance of a $T-thread launch is '$a'"
  snap_top
  log "free ${free} GB; allowance $a for $T threads"
  unit=$(next_unit)
  echo "$unit" >> "$C/units"
  touch "$C/run/since"
  launch "$ROOT" "$PROFILE" "$unit" 200G "$Q" "$C/run/$unit" --stage "$STAGE" --release "$REL" --stages "$ALL" \
    --genomes "$IN/genomes.json" --objects "$IN/objects.json" --library pfam 38.2 "$PFAM" gathering --crispr \
    --threads "$T" --scan-jobs "$J" --scan-cpus "$CPUS"
  sleep 9
  while systemctl --user is-active --quiet "$unit"; do
    if [ -e "$C/run/$unit/foreign_stop" ] && [ ! -e "$C/run/$unit/foreign_logged" ]; then
      touch "$C/run/$unit/foreign_logged"
      log "run: a STOP not ours is in the root; waiting for $unit to end (build-stage stops at its next step boundary)"
    fi
    if [ $((k % 6)) = 0 ]; then
      line=$(python3 -I "$BIN/r50_watch.py" --root "$ROOT" --scratch "$SCRATCH" --since "$C/run/since" --state "$C/run" \
        --unit-prefix c26-pl1- --tmp-min-gb "$TMP_MIN_GB" --scratch-max-gb 100 --cpu-max "$CPU_MAX" --stages "$STAGE")
      rc=$?
      echo "$line" >> "$C/run/guard.jsonl"
      [ "$rc" = 3 ] && stop_unit "$unit" "guard: $line"
      [ $((k % 60)) = 0 ] && log "guard: $(echo "$line" | cut -c1-400); unit log: $(grep -E 'run-step|build-stage' "$LOGS/$unit.log" | tail -n 1 | cut -c1-240)"
    fi
    n=$(stage_records)
    if [ "$n" != "$nrec" ]; then
      nrec=$n
      line=$(bb py "$BIN/c26_est.py" rates --root "$ROOT" --objects "$IN/objects.json" --stage "$STAGE" --expect 2>> "$LOG")
      rc=$?
      echo "$line" >> "$C/run/rates.jsonl"
      log "rates ($n records of $STAGE): $(echo "$line" | cut -c1-400)"
      [ "$rc" = 3 ] && stop_unit "$unit" "a measurement more than 20% from its estimate (AGENTS §6, DEC-43): $(echo "$line" | cut -c1-900)"
      [ "$rc" = 4 ] && drain_stop "$unit" "a gene-call batch more than 20% from its estimate (AGENTS §6, DEC-43): $(echo "$line" | cut -c1-900)"
    fi
    k=$((k + 1))
    sleep 10
  done
  line=$(bb py "$BIN/c26_est.py" rates --root "$ROOT" --objects "$IN/objects.json" --stage "$STAGE" --expect 2>> "$LOG")
  echo "$line" >> "$C/run/rates.jsonl"
  log "rates at the end: $(echo "$line" | cut -c1-600)"
  ended "$ROOT" "$unit" "$STAGE" "$REL" "$C/run/$unit"
  rc=$?
  [ "$rc" = 3 ] && fail "a STOP not written by this driver appeared in the root during $unit; it is left in place (read $LOGS/$unit.log)"
  [ "$rc" = 0 ] || fail "$unit did not end as expected (read $LOGS/$unit.log and $LOG)"
  log "run: PL1's stage steps done; the launch stopped at $REL's cluster step by this driver's STOP, which is removed"
  touch "$C/run/DONE"
}

# --- check ----------------------------------------------------------------------------------------------------------

check() {
  local unit rc
  mkdir -p "$C/check"
  unit=$(tail -n 1 "$C/units")
  bb py "$BIN/c26_check.py" --root "$ROOT" --before "$C/before" --after "$C/after" --log "$LOGS/$unit.log" \
    --objects "$IN/objects.json" --out "$C/check" --tmp "$C/tmp/check" --stage "$STAGE" --release "$REL" \
    --earlier "$EARLIER" --expect-owned "$IN/shared_records.tsv" --expect-counts > "$C/check/check.out" 2>&1
  rc=$?
  cp "$C/check/root_added.tsv" "$C/root_added.tsv"
  log "check: exit $rc: $(tail -n 1 "$C/check/check.out" | cut -c1-600)"
  return $rc
}

MODE=${1:-}
log "c26 driver: mode '$MODE' (pid $$), tree $TREE, image $IMAGE"
case "$MODE" in
  rehearse) rehearse ;;
  estimate) estimate ;;
  before) before ;;
  run) run ;;
  after) after ;;
  check) check ;;
  *) echo "usage: c26_driver.sh rehearse|estimate|before|run|after|check" >&2; exit 2 ;;
esac
log "c26 driver: mode '$MODE' done"
