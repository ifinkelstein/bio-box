#!/bin/bash
# Card C-24 (amended 2026-11-15): stage H1, the huge phages, built into the data root, the stage only, with Pfam 38.2,
# both CRISPR producers and its label files (ingest writes them); no release step. Host side, bash. Copies of this
# script and c24_prep.py, c24_est.py, c24_check.py, c24_probe.py, r50_watch.py run from /tmp/if2453/c24/bin, so the
# synced tree is not touched while they run. Each mode is launched in its own small unit, for example
#
#   systemd-run --user --unit=c24-driver-run --collect -p MemoryMax=4G -p MemorySwapMax=0 \
#     -p StandardOutput=append:/tmp/if2453/c24/driver.out -p StandardError=append:/tmp/if2453/c24/driver.out \
#     nice -n 19 bash /tmp/if2453/c24/bin/c24_driver.sh run
#
# The build is phase B1 of r50_driver.sh for one stage: `build-stage --stage H1 --release rM.3 --stages <the nine>`.
# build-stage has no stage-only mode: it runs the stage's steps (inputs, ingest with its contig_label, gene calling,
# the default CRISPR producer, IGRs, tapes, the Pfam scan in shards, the filtered CRISPR producer) and then refuses at
# rM.3's cluster step, exit 2, because stage PL1 (card C-26, not built) has no gene-call run record: the gate of card
# B-36 (`_cluster_gate` in src/biobox/stage.py). That refusal, with H1's last stage step recorded and no file of a
# release step of rM.3, is the expected end; anything else is a failure. build-stage writes no release_def before the
# cluster step (the cluster step defines the release), and `release define` refuses a stage that is not ingested, so
# rM.3 is named by the command's --stages only (the rehearsal shows both).
#
# Modes, in the order the card runs them:
#   rehearse  The fixture of RUNBOOK Appendix B in a scratch root C/rh/root (its own profile, tmp and units c24-rh-*):
#             stage S0 and release rC24.1 built whole; then stage S1 (T4) with the real command's form, release
#             rC24.3 over S0,S1 and the missing stage PLX, --library pfam 38.2 <the real Pfam file> gathering --crispr;
#             it must end with the cluster gate's refusal, exit 2, every S1 step recorded, nothing of rC24.3 written;
#             the same launch again must skip every step and refuse the same way; `release define` over S0,S1,PLX
#             must refuse; `release verify` of rC24.1 must pass before and after with the same manifest hash.
#             Writes C/rh/OK. Nothing is written in ROOT (it is bound read-only in every short command).
#   estimate  c24_est.py before: the code's estimates of ingest and the CRISPR steps, and the scan and gene calling
#             scaled from stage S1 (derived), into C/estimate/before.json.
#   before    Refuses when ROOT/STOP exists or a build unit is active; then `release verify` of every published
#             release (root bound read-only) into C/before/verify-<R>.log, and the snapshot into C/before: listing.tsv
#             (path, bytes, mtime of every file of ROOT), s1.sha256 (every file under a stage=S1 directory). Writes
#             C/before/DONE when every release verifies.
#   run       Needs C/before/DONE and C/rh/OK. Before the launch: no STOP, no build unit, at least START_MIN_GB free
#             on /tmp, the allowance of a T-thread launch (`python -m biobox.runstep`) at least T. One unit
#             c24-h1-<n> (pod_run.sh --threads Q, the profile's 200G), build-stage --threads T --scan-jobs J
#             --scan-cpus CPUS. While it runs, once a minute, r50_watch.py (units c24-h1-*, run records of H1 since
#             the launch): a step over CPU_MAX CPU-minutes, /tmp under TMP_MIN_GB free or scratch over 100 GB stops
#             the unit (systemctl --user stop: the running steps publish nothing; nothing in the root is cleaned) and
#             the driver fails; and each time H1 has a new run record, c24_est.py rates: a step more than 20% from its
#             estimate stops the unit the same way (AGENTS §6). Writes C/run/DONE when the unit ends as expected.
#   after     The snapshot into C/after, `release verify` of every published release again, `bio-box releases`.
#   e16       The chance check of card E-16 on H1, read-only on ROOT: a scratch view C/e13view of H1's and S1's
#             sequence tables (e13_crispr_min_repeats.py view), the default producer there (units c24-e16-nr3: its
#             files must be byte-identical to the root's) and MinCED at 2 repeats (c24-e16-nr2,
#             e13_crispr_min_repeats.py produce), then e16_crispr_near_cas_phages.py measure --stages H1 into C/e16.
#   check     c24_check.py: before against after, root_added.tsv (copied to C/root_added.tsv), the steps against their
#             estimates, the data the card asks for, into C/check.
set -u

C=/tmp/if2453/c24
ROOT=/tmp/if2453/r10/root
TREE=${TREE:-/tmp/if2453/ci/C-24}
REALPROF=$TREE/profiles/wilkcomp02.yaml
PROFILE=$C/profile.yaml
IMAGE=/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif
DIGEST=sha256:d67300aa1321b571fa8da79e74dfccdd87a8dca2e56518f324902efa6b9f098a
# The test image of the same commit, only for the rehearsal's fixture inputs (c24_prep.py imports the tests' helpers,
# which import pytest; the runtime image has no pytest). Nothing that writes a stage runs in it.
TESTIMAGE=/stor/home/if2453/bio-box/images/bio-box-test_faec1f4.sif
PYPATH=$TREE/src
PFAMDIR=/stor/scratch/Wilke/if2453/all-the-data/l0/refdb/pfam/38.2
PFAM=$PFAMDIR/Pfam-A.hmm.gz
IN=/tmp/if2453/s2/inputs/H1
STAGE=H1
REL=rM.3
ALL=S1,P1,P2,P3,P4,P5,P6,H1,PL1
RELEASES="rN.1a rN.1b rN.1 rN.2 rN.3 rM.2"
# Derived (planning/measurements/2026-11-11_s2_inputs.md): the 7,228 distinct sequences less the 514 stage S1 stores.
H1_BASES=1819831338
T=${T:-36}
J=${J:-9}
CPUS=${CPUS:-4}
Q=${Q:-40}
CPU_MAX=${CPU_MAX:-280}
TMP_MIN_GB=${TMP_MIN_GB:-200}
START_MIN_GB=${START_MIN_GB:-220}
BIN=$(cd "$(dirname "$0")" && pwd)
LOGS=$C/logs
LOG=$C/driver.log
SCRATCH=$C/tmp/scratch
NR3=minced-0.4.2-w8-nr3-r23x47-s26x50
NR2=minced-0.4.2-w8-nr2-r23x47-s26x50
NEAR=minced-0.4.2-w8-nr2-r23x47-s26x50-cas1kb-caspfam1-pfam-38.2
mkdir -p "$LOGS" "$C/tmp"

log() { echo "$(date '+%F %T %Z') $*" >> "$LOG"; }
fail() { log "FAILED: $*"; mkdir -p "$C/$MODE"; echo "$*" > "$C/$MODE/FAILED"; exit 1; }

profile() {  # the real profile with tmp under C (one line differs)
  sed "s#^tmp:.*#tmp: $C/tmp#" "$REALPROF" > "$PROFILE"
  log "profile $PROFILE: $(diff "$REALPROF" "$PROFILE" | grep -c '^>') line differs from $REALPROF"
}

bb() {  # a short command in the image, in its own capped scope at nice 19, the root and Pfam bound read-only
  local what=$1
  shift
  case "$what" in
    cli) set -- python -m biobox.cli "$@" ;;
    py) set -- python "$@" ;;
  esac
  systemd-run --user --scope --quiet -p MemoryMax=16G -p MemorySwapMax=0 -p CPUQuota=200% nice -n 19 \
    singularity exec --cleanenv --env PYTHONPATH="$PYPATH" --env OMP_NUM_THREADS=2 --env PYTHONUNBUFFERED=1 \
    --bind /tmp/if2453 --bind "$ROOT:$ROOT:ro" --bind "$PFAMDIR:$PFAMDIR:ro" --bind /sys/fs/cgroup:/sys/fs/cgroup:ro \
    "${IMG:-$IMAGE}" "$@"
}

build_units() {  # active units of builds that write into ROOT, other than this card's driver
  systemctl --user list-units --type=service --state=active --no-legend --plain \
    'r50-*' 'c24-h1-*' 'c25-*' 'c26-*' 'c33-*' 'c34-*' 'c35-*' 'b47-*' | awk '{print $1}' | paste -sd' ' -
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

wait_plain() {  # $1 unit: wait until it ends (short rehearsal and view units)
  sleep 5
  while systemctl --user is-active --quiet "$1"; do sleep 15; done
  log "$1 ended: $(tail -n 2 "$LOGS/$1.log" | tr '\n' ' ')"
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
  find "$ROOT" ! -type d -path '*/stage=S1/*' -printf '%P\n' | LC_ALL=C sort > "$d/s1.files"
  (cd "$ROOT" && systemd-run --user --scope --quiet -p MemoryMax=2G -p MemorySwapMax=0 nice -n 19 \
    xargs -d '\n' sha256sum < "$d/s1.files") > "$d/s1.sha256"
  log "snapshot $d: $(wc -l < "$d/listing.tsv") files, $(wc -l < "$d/s1.sha256") of S1 hashed, $(awk -F '\t' '{s += $2} END {printf "%.0f", s}' "$d/listing.tsv") bytes"
}

# --- rehearse -------------------------------------------------------------------------------------------------------

rh_launch() {  # $1 unit, $2 mem, rest: build-stage arguments after --root
  local unit=$1 mem=$2
  shift 2
  log "launch $unit: build-stage --profile $RH/profile.yaml --root $RH/root $*"
  bash "$TREE/scripts/pod_run.sh" --profile "$RH/profile.yaml" --image "$IMAGE" --digest "$DIGEST" --mem "$mem" \
    --threads 4 --bind /tmp/if2453 --bind "$ROOT:$ROOT:ro" --bind "$PFAMDIR:$PFAMDIR:ro" --env "PYTHONPATH=$PYPATH" \
    --log "$LOGS/$unit.log" "$unit" build-stage --profile "$RH/profile.yaml" --root "$RH/root" "$@" >> "$LOG" 2>&1
  wait_plain "$unit"
}

rehearse() {
  RH=$C/rh
  local rc s1=(--stage S1 --release rC24.3 --stages S0,S1,PLX --genomes "$C/rh/inputs/genomes-S1.json"
    --objects "$C/rh/inputs/objects-S1.json" --library pfam 38.2 "$PFAM" gathering --crispr --threads 4
    --scan-jobs 1 --scan-cpus 4) step n
  [ -e "$RH/OK" ] && { log "rehearsal: OK present; nothing to do"; return; }
  mkdir -p "$RH/tmp" "$RH/durable" "$RH/bulk" "$RH/check"
  cat > "$RH/profile.yaml" << EOF
max_threads_day: 4
max_threads_night: 4
mem_scope: 16G
tmp: $RH/tmp
durable: $RH/durable
bulk: $RH/bulk
nice_monolithic: 19
nice_fanout: 15
day_start_hour: 8
day_end_hour: 20
cpu_min_per_process: 300
small_file_floor_bytes: 1000
require_mem_cap: true
EOF
  if [ ! -s "$RH/inputs/objects-S1.json" ]; then
    IMG=$TESTIMAGE bb py "$BIN/c24_prep.py" "$TREE" "$RH" >> "$LOG" 2>&1 || fail "rehearsal: c24_prep.py"
  fi
  if [ ! -e "$RH/root/runs/stage=rC24.1/catalog.json" ]; then
    rh_launch c24-rh-s0 8G --stage S0 --release rC24.1 --stages S0 --genomes "$RH/inputs/genomes-S0.json" \
      --objects "$RH/inputs/objects-S0.json" --library test 1 "$RH/inputs/lib.hmm" gathering --threads 4 \
      --scan-jobs 1 --scan-cpus 4
    [ "$(exit_of c24-rh-s0)" = 0 ] || fail "rehearsal: c24-rh-s0 exit '$(exit_of c24-rh-s0)'"
  fi
  verify "$RH/root" rC24.1 "$RH/check" || fail "rehearsal: rC24.1 does not verify before the S1 launch"
  mv "$RH/check/verify-rC24.1.log" "$RH/check/verify-rC24.1-before.log"
  find "$RH/root" ! -type d -printf '%P\t%s\t%T@\n' | LC_ALL=C sort > "$RH/check/listing-before.tsv"
  rh_launch c24-rh-s1 16G "${s1[@]}"
  rc=$(exit_of c24-rh-s1)
  log "rehearsal: c24-rh-s1 exit '$rc'; $(grep 'build-stage: refused' "$LOGS/c24-rh-s1.log" | tail -n 1)"
  [ "$rc" = 2 ] || fail "rehearsal: c24-rh-s1 exit '$rc', not 2"
  grep -q "build-stage: refused before rC24.3.cluster: release rC24.3: stage PLX has no gene-call run record" \
    "$LOGS/c24-rh-s1.log" || fail "rehearsal: no cluster-gate refusal naming PLX"
  for step in inputs ingest genecall "crispr-$NR3" igr tapes domains-pfam-38.2 "crispr-$NEAR"; do
    [ -e "$RH/root/runs/stage=S1/$step.json" ] || fail "rehearsal: no run record S1/$step"
  done
  [ -e "$RH/root/runs/stage=rC24.3" ] && fail "rehearsal: runs/stage=rC24.3 exists: a release step ran"
  [ -e "$RH/root/release_def/release=rC24.3" ] && fail "rehearsal: rC24.3 has a release_def"
  n=$(find "$RH/root" -path '*release=rC24.3*' -o -path '*stage=rC24.3*' | wc -l)
  [ "$n" = 0 ] || fail "rehearsal: $n paths of rC24.3 in the root"
  find "$RH/root" ! -type d -printf '%P\t%s\t%T@\n' | LC_ALL=C sort > "$RH/check/listing-after.tsv"
  rh_launch c24-rh-s1-again 16G "${s1[@]}"
  rc=$(exit_of c24-rh-s1-again)
  log "rehearsal: c24-rh-s1-again exit '$rc'; $(grep -c 'skipped: run record and outputs present' "$LOGS/c24-rh-s1-again.log") steps skipped; $(grep 'build-stage: refused' "$LOGS/c24-rh-s1-again.log" | tail -n 1)"
  [ "$rc" = 2 ] || fail "rehearsal: the second launch exit '$rc', not 2"
  grep -q "run-step: S1\.genecall skipped" "$LOGS/c24-rh-s1-again.log" || fail "rehearsal: the second launch ran gene calling again"
  find "$RH/root" ! -type d -printf '%P\t%s\t%T@\n' | LC_ALL=C sort > "$RH/check/listing-again.tsv"
  cmp -s "$RH/check/listing-after.tsv" "$RH/check/listing-again.tsv" || fail "rehearsal: the second launch changed the root"
  bb cli release define --root "$RH/root" --release rC24.3x --stages S0,S1,PLX > "$RH/check/define.log" 2>&1
  rc=$?
  log "rehearsal: release define over S0,S1,PLX: exit $rc: $(tail -n 1 "$RH/check/define.log")"
  [ "$rc" = 2 ] || fail "rehearsal: release define over a stage with no files exit '$rc', not 2"
  verify "$RH/root" rC24.1 "$RH/check" || fail "rehearsal: rC24.1 does not verify after the S1 launches"
  [ "$(grep -E '^release ' "$RH/check/verify-rC24.1.log")" = "$(grep -E '^release ' "$RH/check/verify-rC24.1-before.log")" ] \
    || fail "rehearsal: rC24.1's verify line changed"
  LC_ALL=C comm -23 "$RH/check/listing-before.tsv" "$RH/check/listing-after.tsv" > "$RH/check/changed.tsv"
  [ -s "$RH/check/changed.tsv" ] && fail "rehearsal: files of the root before the S1 launch changed or went"
  touch "$RH/OK"
  log "rehearsal: OK"
}

# --- estimate, before, after ----------------------------------------------------------------------------------------

estimate() {
  mkdir -p "$C/estimate"
  profile
  bb py "$BIN/c24_est.py" before --root "$ROOT" --profile "$PROFILE" --objects "$IN/objects.json" --stage "$STAGE" \
    --bases "$H1_BASES" --like S1 > "$C/estimate/before.json" 2> "$C/estimate/before.err"
  log "estimate: exit $?: $(tr -d '\n ' < "$C/estimate/before.json" | cut -c1-600)"
}

before() {
  local d=$C/before r ok=1 act
  [ -e "$d/DONE" ] && { log "before: DONE present; nothing to do"; return; }
  [ -e "$ROOT/STOP" ] && fail "$ROOT/STOP exists"
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
  echo "c24-h1-$n"
}

h1_records() { find "$ROOT/runs/stage=$STAGE" -maxdepth 1 -name '*.json' ! -name '*.publishing.json' 2> /dev/null | wc -l; }

stop_unit() {  # $1 unit, $2 reason
  log "stopping $1: $2"
  systemctl --user stop "$1"
  sleep 10
  fail "$2 (unit $1 stopped; the root is left as the steps left it; the same command resumes)"
}

run() {
  local unit a free act rc line nrec=-1 n k=0
  [ -e "$C/before/DONE" ] || fail "no $C/before/DONE: run the before mode first"
  [ -e "$C/rh/OK" ] || fail "no $C/rh/OK: the rehearsal has not passed"
  [ -e "$C/run/DONE" ] && { log "run: DONE present; nothing to do"; return; }
  mkdir -p "$C/run"
  profile
  [ -e "$ROOT/STOP" ] && fail "$ROOT/STOP exists"
  act=$(build_units)
  [ -z "$act" ] || fail "a build unit is active: $act"
  free=$(free_gb)
  [ "$free" -ge "$START_MIN_GB" ] || fail "/tmp has $free GB free, under $START_MIN_GB: not starting"
  a=$(bb py -m biobox.runstep --profile "$PROFILE" --threads "$T" 2>> "$LOG")
  [[ "$a" =~ ^[0-9]+$ ]] && [ "$a" -ge "$T" ] || fail "the allowance of a $T-thread launch is '$a'"
  snap_top
  log "free ${free} GB; allowance $a for $T threads"
  touch "$C/run/since"
  unit=$(next_unit)
  echo "$unit" >> "$C/units"
  log "launch $unit (pod_run --threads $Q): build-stage --profile $PROFILE --root $ROOT --stage $STAGE --release $REL --stages $ALL --genomes $IN/genomes.json --objects $IN/objects.json --library pfam 38.2 $PFAM gathering --crispr --threads $T --scan-jobs $J --scan-cpus $CPUS"
  bash "$TREE/scripts/pod_run.sh" --profile "$PROFILE" --image "$IMAGE" --digest "$DIGEST" --threads "$Q" \
    --bind /tmp/if2453 --bind "$PFAMDIR:$PFAMDIR:ro" --env "PYTHONPATH=$PYPATH" --log "$LOGS/$unit.log" "$unit" \
    build-stage --profile "$PROFILE" --root "$ROOT" --stage "$STAGE" --release "$REL" --stages "$ALL" \
    --genomes "$IN/genomes.json" --objects "$IN/objects.json" --library pfam 38.2 "$PFAM" gathering --crispr \
    --threads "$T" --scan-jobs "$J" --scan-cpus "$CPUS" >> "$LOG" 2>&1
  sleep 10
  while systemctl --user is-active --quiet "$unit"; do
    line=$(python3 -I "$BIN/r50_watch.py" --root "$ROOT" --scratch "$SCRATCH" --since "$C/run/since" --state "$C/run" \
      --unit-prefix c24-h1- --tmp-min-gb "$TMP_MIN_GB" --scratch-max-gb 100 --cpu-max "$CPU_MAX" --stages "$STAGE")
    rc=$?
    echo "$line" >> "$C/run/guard.jsonl"
    [ "$rc" = 3 ] && stop_unit "$unit" "guard: $line"
    k=$((k + 1))
    [ $((k % 10)) != 1 ] || log "guard: $(echo "$line" | cut -c1-400); unit log: $(grep -E 'run-step|build-stage' "$LOGS/$unit.log" | tail -n 1 | cut -c1-240)"
    n=$(h1_records)
    if [ "$n" != "$nrec" ]; then
      nrec=$n
      line=$(bb py "$BIN/c24_est.py" rates --root "$ROOT" --objects "$IN/objects.json" --stage "$STAGE" 2>> "$LOG")
      rc=$?
      echo "$line" >> "$C/run/rates.jsonl"
      log "rates ($n records of $STAGE): $(echo "$line" | cut -c1-300)"
      [ "$rc" = 3 ] && stop_unit "$unit" "a step more than 20% from its estimate (AGENTS §6): $(echo "$line" | cut -c1-600)"
    fi
    sleep 50
  done
  rc=$(exit_of "$unit")
  log "$unit ended, exit '$rc'; $(grep -E 'pod_run: .* usage' "$LOGS/$unit.log" | tail -n 1)"
  log "gene-call result: $(grep '"step": "genecall"' "$LOGS/$unit.log" | tail -n 1 | cut -c1-800)"
  line=$(bb py "$BIN/c24_est.py" rates --root "$ROOT" --objects "$IN/objects.json" --stage "$STAGE" 2>> "$LOG")
  echo "$line" >> "$C/run/rates.jsonl"
  log "rates at the end: $(echo "$line" | cut -c1-300)"
  [ -e "$ROOT/runs/stage=$REL" ] && fail "runs/stage=$REL exists: a release step of $REL ran"
  if [ "$rc" = 2 ] && [ -e "$ROOT/runs/stage=$STAGE/crispr-$NEAR.json" ] \
     && grep -q "build-stage: refused before $REL\.cluster: release $REL: stage PL1 has no gene-call run record" "$LOGS/$unit.log"; then
    log "run: H1's stage steps done; build-stage stopped at $REL's cluster gate as expected: $(grep "refused before $REL.cluster" "$LOGS/$unit.log" | tail -n 1)"
    touch "$C/run/DONE"
    return
  fi
  fail "$unit exit '$rc' is not the expected cluster-gate refusal (read $LOGS/$unit.log)"
}

# --- e16, check -----------------------------------------------------------------------------------------------------

view_step() {  # $1 unit, $2 step, rest: the step command after `--`
  local unit=$1 step=$2
  shift 2
  [ -e "$C/e13view/runs/stage=$STAGE/$step.json" ] && { log "e16: $step in the view: record present"; return 0; }
  mkdir -p "$C/tmp/view-scratch/$STAGE.$step"
  log "launch $unit: run-step --root $C/e13view --stage $STAGE --step $step --threads 4 -- $*"
  bash "$TREE/scripts/pod_run.sh" --profile "$PROFILE" --image "$IMAGE" --digest "$DIGEST" --mem 16G --threads 4 \
    --bind /tmp/if2453 --bind "$ROOT:$ROOT:ro" --env "PYTHONPATH=$PYPATH" --log "$LOGS/$unit.log" "$unit" \
    run-step --profile "$PROFILE" --root "$C/e13view" --stage "$STAGE" --step "$step" --threads 4 -- "$@" >> "$LOG" 2>&1
  wait_plain "$unit"
  [ "$(exit_of "$unit")" = 0 ] || fail "e16: $unit exit '$(exit_of "$unit")'"
}

e16() {
  local V=$C/e13view rc
  mkdir -p "$C/e16"
  profile
  if [ ! -d "$V" ]; then
    bb py "$TREE/scripts/stress/e13_crispr_min_repeats.py" view "$ROOT" "$V" "$STAGE,S1" >> "$LOG" 2>&1 || fail "e16: view"
  fi
  view_step c24-e16-nr3 "crispr-$NR3" python -m biobox.cli stage crispr --root "$V" --stage "$STAGE" --cpus 4 \
    --tmp "$C/tmp/view-scratch/$STAGE.crispr-$NR3" --seq-stages "$STAGE,S1"
  view_step c24-e16-nr2 "crispr-$NR2" python "$TREE/scripts/stress/e13_crispr_min_repeats.py" produce --root "$V" \
    --stage "$STAGE" --cpus 4 --min-repeats 2 --tmp "$C/tmp/view-scratch/$STAGE.crispr-$NR2" --seq-stages "$STAGE,S1"
  rm -rf "$C/e16/out"
  bb py "$TREE/scripts/stress/e16_crispr_near_cas_phages.py" measure --root "$ROOT" --calls "$V" --kept "$ROOT" \
    --tree "$TREE" --out "$C/e16/out" --stages "$STAGE" --seq-stages "$STAGE,S1" --draws 10000 > "$C/e16/measure.log" 2>&1
  rc=$?
  log "e16: measure --stages $STAGE exit $rc: $(tail -n 2 "$C/e16/measure.log" | tr '\n' ' ' | cut -c1-400)"
}

check() {
  local unit
  mkdir -p "$C/check"
  unit=$(tail -n 1 "$C/units")
  bb py "$BIN/c24_check.py" --root "$ROOT" --before "$C/before" --after "$C/after" --log "$LOGS/$unit.log" \
    --objects "$IN/objects.json" --out "$C/check" --tmp "$C/tmp/check" --stage "$STAGE" --release "$REL" \
    --view "$C/e13view" > "$C/check/check.out" 2>&1
  local rc=$?
  cp "$C/check/root_added.tsv" "$C/root_added.tsv"
  log "check: exit $rc: $(tail -n 1 "$C/check/check.out")"
  return $rc
}

MODE=${1:-}
log "c24 driver: mode '$MODE' (pid $$), tree $TREE, image $IMAGE"
case "$MODE" in
  rehearse) rehearse ;;
  estimate) estimate ;;
  before) before ;;
  run) run ;;
  after) after ;;
  e16) e16 ;;
  check) check ;;
  *) echo "usage: c24_driver.sh rehearse|estimate|before|run|after|e16|check" >&2; exit 2 ;;
esac
log "c24 driver: mode '$MODE' done"
