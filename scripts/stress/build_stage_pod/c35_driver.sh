#!/bin/bash
# Card C-35 (amended 2026-11-16, DEC-52): release r26.11.1, the first rolling release, over the nine stages of the data
# root in the order they were built, S1,P1,P2,P3,P4,P5,P6,H1,PL1. (Cards C-24 and C-26 and their records call it rM.3:
# that name was only an argument of their launches, and nothing in the root carries it.) Host side, bash. Copies of
# this script and c35_est.py, c35_calib.py, c35_check.py, c24_prep.py, r50_watch.py run from /tmp/if2453/c35/bin, and
# the build reads the tree copied to TREE (/tmp/if2453/c35/tree), so later syncs of the tree do not touch them while
# they run. Each mode is launched in its own small unit, for example
#
#   systemd-run --user --unit=c35-driver-run --collect -p MemoryMax=4G -p MemorySwapMax=0 \
#     -p StandardOutput=append:/tmp/if2453/c35/driver.out -p StandardError=append:/tmp/if2453/c35/driver.out \
#     nice -n 19 bash /tmp/if2453/c35/bin/c35_driver.sh run
#
# The build is the one `build-stage` command of the card, in the runtime image of faec1f4 by digest (RUNBOOK §0):
#
#   build-stage --stage PL1 --release r26.11.1 --stages S1,P1,P2,P3,P4,P5,P6,H1,PL1 --genomes/--objects <PL1's inputs>
#     --library pfam 38.2 <Pfam 38.2> gathering --crispr --threads T --scan-jobs J --scan-cpus CPUS --cluster-parts 24
#
# Every stage step of the list (PL1's inputs to its filtered CRISPR producer) has its run record, so each is skipped;
# the release steps run: clustering in parts (card B-41; no earlier release qualifies for card B-49's reuse, since none
# is over these stages: the `estimate` mode prints the code's reason for each), tokens, postings, dim, manifest,
# catalog. Why --cluster-parts 24 and not the default cut (planning/measurements/2026-11-16_release_r26.11.1.md):
# rM.2 was cut so; its calls cost 0.78 to 1.34 times the cut's estimates (card C-22's model at 45.8 M proteins), so a
# default cut, whose parts are estimated at up to 240 CPU-minutes, would give parts of up to about 320 for the calls the
# model underestimates; at 24 parts the largest part is projected at 169 CPU-minutes. J x CPUS = 2 x 18: the NR90
# step's prefilter parts each hold the k-mer table of every target (56.4 GB peak RSS on rM.2, about 60 GB here,
# derived), so 2 at once (about 121 GB) stay well under the unit's 200G where 3 (about 181 GB plus page cache) would not;
# 36 threads fit the night cap while the account's other work uses up to 20 (17 was seen by day during the rehearsal);
# by day the allowance (40 less ours elsewhere) runs one part, and T = 18 keeps a one-process step from being held.
#
# Modes, in the order the card runs them:
#   rehearse  The fixture of RUNBOOK Appendix B in a scratch root C/rh/root (own profile, tmp and units c35-rh-*;
#             the data root is bound read-only in every command): stage S0 and release rC35.0 built whole (Pfam 38.2,
#             --crispr); stage S1's steps alone (release rC35.gate over S0,S1,PLX ends at the cluster gate, exit 2);
#             then the real command's form, release rC35.rh over S0,S1 with every stage step skipped and the release
#             in parts (--cluster-parts 2), with the guard and the rate check; `release verify` of rC35.rh and rC35.0
#             (its hash unchanged); its pack in parts (--parts 2) and `pack verify`; the release check as far as the
#             fixture allows (it holds none of the examples' controls, so their invariants are missing: exit 1); the
#             same launch again (every step skipped, the view's listing unchanged); c35_check.py on it. Writes C/rh/OK.
#   estimate  c35_calib.py on rM.2 (its clustering call by call against its cut's estimates), then c35_est.py before:
#             proteins, genes, the code's decisions (in parts, the cut, the reuse reasons), every step's projection,
#             the memory and scratch plan, into C/estimate/before.json; fails when a part is projected at 240 or more.
#   before    Refuses when ROOT/STOP exists or a build unit is active; `release verify` of the six published releases
#             (root bound read-only) into C/before/verify-<R>.log, each exit 0; listing.tsv (path, bytes, mtime of
#             every file of ROOT). Writes C/before/DONE.
#   run       Needs C/before/DONE, C/rh/OK, C/estimate/before.json. Before the launch: no STOP, no HALT, no build unit,
#             at least START_MIN_GB free on /tmp, `ps`/`top` into the log, the allowance of a T-thread launch at least
#             T. One unit c35-rel-<n> (pod_run.sh --threads Q, the profile's 200G, nice 19). While it runs: every 10 s
#             the unit and STOP; once a minute r50_watch.py (units c35-rel-*: a process or record over CPU_MAX
#             CPU-minutes, /tmp under TMP_MIN_GB free, scratch over SCRATCH_MAX_GB) and at each new run record of the
#             release c35_est.py rates (DEC-43, DEC-53: a step estimated at 24 CPU-minutes or more that exceeds its
#             estimate by more than 20% stops the unit at once; steps more than 20% under are recorded in rates.jsonl
#             and reported without stopping): either stops the unit at once (systemctl --user stop: the running steps
#             publish nothing) and the driver fails; the same command resumes. A STOP not written by this driver:
#             nothing is written, the driver waits for the unit (build-stage stops at its next step boundary) and fails.
#             Writes C/run/DONE when the unit ends with exit 0, "all 15 steps done", every stage step skipped.
#   after     The listing into C/after; `release verify` of r26.11.1 and of the six (their hashes against C/before);
#             `bio-box releases`.
#   pack      The search pack in parts (card W-35) into the profile's pack directory /tmp/if2453/bio-box/packs, unit
#             c35-pack (pack build --parts PACK_PARTS --jobs PACK_JOBS), then `pack verify`.
#   examples  scripts/release_examples.sh r26.11.1 (RUNBOOK §6) into C/examples, unit c35-ex.
#   check     c35_check.py: before against after, root_added.tsv (also C/root_added.tsv), the steps against their
#             estimates, clusters, new proteins at NR90 and NR50 per stage, cluster_uid stability against rM.2, arrays
#             per producer, label files, the DEC-27 checks, into C/check.
set -u

C=/tmp/if2453/c35
ROOT=/tmp/if2453/r10/root
TREE=${TREE:-$C/tree}
REALPROF=$TREE/profiles/wilkcomp02.yaml
PROFILE=$C/profile.yaml
IMAGE=/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif
DIGEST=sha256:d67300aa1321b571fa8da79e74dfccdd87a8dca2e56518f324902efa6b9f098a
# The test image of the same commit, only for the rehearsal's fixture inputs (c24_prep.py imports the tests' helpers,
# which import pytest). Nothing that writes a stage or a release runs in it.
TESTIMAGE=/stor/home/if2453/bio-box/images/bio-box-test_faec1f4.sif
PYPATH=$TREE/src
PFAMDIR=/stor/scratch/Wilke/if2453/all-the-data/l0/refdb/pfam/38.2
PFAM=$PFAMDIR/Pfam-A.hmm.gz
IN=/tmp/if2453/s2/inputs/PL1
STAGE=PL1
REL=r26.11.1
ALL=S1,P1,P2,P3,P4,P5,P6,H1,PL1
EARLIER=rM.2
EARLIER_STAGES=S1,P1,P2,P3,P4,P5,P6
RELEASES="rN.1a rN.1b rN.1 rN.2 rN.3 rM.2"
NSTEPS=15
PARTS=${PARTS:-24}
T=${T:-18}
J=${J:-2}
CPUS=${CPUS:-18}
Q=${Q:-40}
CPU_MAX=${CPU_MAX:-280}
TMP_MIN_GB=${TMP_MIN_GB:-150}
SCRATCH_MAX_GB=${SCRATCH_MAX_GB:-250}
START_MIN_GB=${START_MIN_GB:-250}
PACK_PARTS=${PACK_PARTS:-10}
PACK_JOBS=${PACK_JOBS:-12}
PACKDIR=/tmp/if2453/bio-box/packs
RH=$C/rh
BIN=$(cd "$(dirname "$0")" && pwd)
LOGS=$C/logs
LOG=$C/driver.log
SCRATCH=$C/tmp/scratch
NR3=minced-0.4.2-w8-nr3-r23x47-s26x50
NEAR=minced-0.4.2-w8-nr2-r23x47-s26x50-cas1kb-caspfam1-pfam-38.2
STAGE_STEPS="inputs ingest genecall crispr-$NR3 igr tapes domains-pfam-38.2 crispr-$NEAR"
mkdir -p "$LOGS" "$C/tmp"
touch "$C/stop_tokens"

log() { echo "$(date '+%F %T %Z') $*" >> "$LOG"; }
fail() { log "FAILED: $*"; mkdir -p "$C/$MODE"; echo "$*" > "$C/$MODE/FAILED"; exit 1; }

profile() {  # the real profile with tmp under C (one line differs): the build's logs and scratch go there
  sed "s#^tmp:.*#tmp: $C/tmp#" "$REALPROF" > "$PROFILE"
  log "profile $PROFILE: $(diff "$REALPROF" "$PROFILE" | grep -c '^>') line differs from $REALPROF"
}

bb() {  # a short command in the image, in its own capped scope at nice 19, the root and Pfam read-only
  local what=$1 img=${IMG:-$IMAGE}
  shift
  case "$what" in
    cli) set -- python -m biobox.cli "$@" ;;
    py) set -- python "$@" ;;
  esac
  systemd-run --user --scope --quiet -p MemoryMax=16G -p MemorySwapMax=0 -p CPUQuota=200% nice -n 19 \
    singularity exec --cleanenv --env PYTHONPATH="$PYPATH" --env OMP_NUM_THREADS=2 --env PYTHONUNBUFFERED=1 \
    --bind /tmp/if2453 --bind "$ROOT:$ROOT:ro" --bind "$PFAMDIR:$PFAMDIR:ro" \
    --bind /sys/fs/cgroup:/sys/fs/cgroup:ro "$img" "$@"
}

build_units() {  # active units of builds that write into ROOT (this card's driver units are c35-driver-*)
  systemctl --user list-units --type=service --state=active --no-legend --plain \
    'r50-*' 'c24-h1-*' 'c25-*' 'c26-pl1-*' 'c33-*' 'c34-*' 'c35-rel-*' 'c35-pack*' 'b47-*' | awk '{print $1}' | paste -sd' ' -
}

snap_top() {
  log "uptime: $(uptime)"
  log "account CPU now (top, %CPU of if2453): $(top -b -n 1 -u if2453 | awk 'NR > 7 {s += $9} END {print s}')"
  log "other users' CPU now (top, %CPU of users other than if2453 and root): $(top -b -n 1 | awk 'NR > 7 && $2 != "if2453" && $2 != "root" {s += $9} END {print s + 0}')"
  top -b -n 1 -o %CPU -w 150 | sed -n '1,5p;7,22p' >> "$LOG"
  log "df /tmp: $(df -h /tmp | tail -n 1)"
  log "active units: $(systemctl --user list-units --type=service --type=scope --state=active --no-legend --plain | awk '{print $1}' | paste -sd' ' - | cut -c1-900)"
}

free_gb() { df -B1G --output=avail /tmp | tail -n 1 | tr -d ' '; }

exit_of() { sed -n "s/^pod_run: $1 exited \([0-9]*\)$/\1/p" "$LOGS/$1.log" | tail -n 1; }

ours_stop() {  # $1 root: its STOP holds one of this driver's token lines
  [ -e "$1/STOP" ] && [ -s "$C/stop_tokens" ] && grep -qxF -f "$C/stop_tokens" "$1/STOP"
}

remove_our_stop() {  # $1 root
  if ours_stop "$1"; then
    rm "$1/STOP"
    log "removed $1/STOP, written by this driver"
  elif [ -e "$1/STOP" ]; then
    log "$1/STOP is NOT this driver's; left as it is: $(head -c 300 "$1/STOP")"
    return 1
  fi
  return 0
}

verify() {  # $1 root, $2 release, $3 directory, [$4 profile file for the image's bind of a scratch root]
  bb cli release verify --root "$1" --release "$2" > "$3/verify-$2.log" 2>&1
  local rc=$?
  log "verify $2 ($1): exit $rc: $(grep -E '^release ' "$3/verify-$2.log" | tail -n 1)"
  return $rc
}

hash_of() { sed -n 's/.*manifest hash \([0-9a-f]*\).*/\1/p' "$1" | tail -n 1; }

listing() {  # $1 root, $2 file: path, bytes and mtime of every file (links listed, not followed), STOP left out
  find "$1" ! -type d ! -name STOP -printf '%P\t%s\t%T@\n' | LC_ALL=C sort > "$2"
}

# --- one launch of build-stage, its guard and its end ---------------------------------------------------------------

launch() {  # $1 root, $2 profile, $3 unit, $4 mem, $5 quota, rest: build-stage arguments after --root
  local root=$1 prof=$2 unit=$3 mem=$4 quota=$5 ro=()
  shift 5
  [ "$root" = "$ROOT" ] || ro=(--bind "$ROOT:$ROOT:ro")
  log "launch $unit (pod_run --threads $quota --mem $mem): build-stage --profile $prof --root $root $*"
  bash "$TREE/scripts/pod_run.sh" --profile "$prof" --image "$IMAGE" --digest "$DIGEST" --mem "$mem" \
    --threads "$quota" --bind /tmp/if2453 "${ro[@]}" --bind "$PFAMDIR:$PFAMDIR:ro" --env "PYTHONPATH=$PYPATH" \
    --log "$LOGS/$unit.log" "$unit" build-stage --profile "$prof" --root "$root" "$@" >> "$LOG" 2>&1
  sleep 2
}

stop_unit() {  # $1 unit, $2 reason
  log "stopping $1: $2"
  systemctl --user stop "$1"
  while systemctl --user is-active --quiet "$1"; do sleep 2; done
  fail "$2 (unit $1 stopped; the root is left as the steps left it; the same command resumes)"
}

watch() {  # $1 root, $2 release, $3 unit, $4 state dir, $5 estimate JSON, $6 unit prefix, $7 scratch: the guard and the
  # rate check while the unit runs; returns 3 when a STOP not ours appeared (the unit then ends at a step boundary)
  local root=$1 rel=$2 unit=$3 st=$4 est=$5 prefix=$6 scr=$7 k=0 n nrec=-1 line rc foreign=0
  while systemctl --user is-active --quiet "$unit"; do
    if [ "$foreign" = 0 ] && [ -e "$root/STOP" ] && ! ours_stop "$root"; then
      foreign=1
      log "a STOP not written by this driver is in $root: $(head -c 300 "$root/STOP"); writing nothing, waiting for $unit to end"
    fi
    if [ $((k % 6)) = 0 ]; then
      line=$(python3 -I "$BIN/r50_watch.py" --root "$root" --scratch "$scr" --since "$st/since" --state "$st" \
        --unit-prefix "$prefix" --tmp-min-gb "$TMP_MIN_GB" --scratch-max-gb "$SCRATCH_MAX_GB" --cpu-max "$CPU_MAX" \
        --stages "$rel")
      rc=$?
      echo "$line" >> "$st/guard.jsonl"
      [ "$rc" = 3 ] && stop_unit "$unit" "guard: $line"
      [ $((k % 60)) = 0 ] && log "guard: $(echo "$line" | cut -c1-500); unit log: $(grep -E 'run-step|build-stage' "$LOGS/$unit.log" | tail -n 1 | cut -c1-300)"
    fi
    n=$(find "$root/runs/stage=$rel" -maxdepth 1 -name '*.json' ! -name '*.publishing.json' 2> /dev/null | wc -l)
    if [ "$n" != "$nrec" ]; then
      nrec=$n
      line=$(python3 -I "$BIN/c35_est.py" rates --root "$root" --release "$rel" --est "$est" 2>> "$LOG")
      rc=$?
      echo "$line" >> "$st/rates.jsonl"
      [ $((n % 20)) = 0 ] && log "rates ($n records of $rel): $(echo "$line" | cut -c1-500)"
      [ "$rc" = 3 ] && stop_unit "$unit" "a judged step more than 20% over its estimate (AGENTS §6, DEC-43, DEC-53): $(echo "$line" | cut -c1-1200)"
      [ "$rc" != 0 ] && [ "$rc" != 3 ] && log "c35_est.py rates exited $rc: $(echo "$line" | cut -c1-300)"
    fi
    k=$((k + 1))
    sleep 10
  done
  line=$(python3 -I "$BIN/c35_est.py" rates --root "$root" --release "$rel" --est "$est" --out "$st/rates-end.json" 2>> "$LOG")
  rc=$?
  log "rates at the end (exit $rc): $(echo "$line" | cut -c1-900)"
  [ "$rc" = 3 ] && fail "a judged step more than 20% over its estimate at the end of $unit (AGENTS §6, DEC-43, DEC-53): $(echo "$line" | cut -c1-1200)"
  [ "$foreign" = 1 ] && return 3
  return 0
}

ended_release() {  # $1 root, $2 unit, $3 stage, $4 release: 0 when the launch ended as expected (every stage step
  # skipped, every release step done, exit 0)
  local root=$1 unit=$2 stage=$3 rel=$4 rc bad=0 s
  rc=$(exit_of "$unit")
  log "$unit ended, exit '$rc'; $(grep -E 'steps done|refused|STOP' "$LOGS/$unit.log" | tail -n 2 | tr '\n' ' ' | cut -c1-600)"
  log "$unit usage: $(grep -E 'pod_run: .* usage' "$LOGS/$unit.log" | tail -n 1)"
  log "$unit allowance lines: $(grep -c 'allowance\|allows' "$LOGS/$unit.log"); first: $(grep -m 1 -E 'allowance' "$LOGS/$unit.log" | cut -c1-300)"
  [ "$rc" = 0 ] || { log "exit '$rc', not 0"; bad=1; }
  grep -q "build-stage: stage $stage and release $rel: all $NSTEPS steps done" "$LOGS/$unit.log" \
    || { log "no line 'all $NSTEPS steps done'"; bad=1; }
  for s in $STAGE_STEPS; do
    grep -q "run-step: $stage\.$s skipped" "$LOGS/$unit.log" || { log "stage step $stage.$s not skipped"; bad=1; }
  done
  grep -qE "run-step: $stage\.[^ ]* published" "$LOGS/$unit.log" && { log "a stage step of $stage published"; bad=1; }
  for s in cluster.plan cluster tokens postings dim manifest catalog; do
    [ -e "$root/runs/stage=$rel/$s.json" ] || { log "no run record $rel/$s"; bad=1; }
  done
  find "$root/runs" -name '*.publishing.json' | grep -q . && { log "an intent file is left"; bad=1; }
  return $bad
}

# --- rehearse -------------------------------------------------------------------------------------------------------

rh_wait() { while systemctl --user is-active --quiet "$1"; do sleep 5; done; }

rehearse() {
  local rc unit s0 s1g s1 n est=$RH/check/est.json
  [ -e "$RH/OK" ] && { log "rehearsal: OK present; nothing to do"; return; }
  mkdir -p "$RH/tmp" "$RH/durable" "$RH/bulk" "$RH/check" "$RH/state"
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
cpu_min_per_process: 300
cluster_cpu_min_max: 20000
small_file_floor_bytes: 1000
require_mem_cap: true
EOF
  if [ ! -s "$RH/inputs/objects-S1.json" ]; then
    IMG=$TESTIMAGE bb py "$BIN/c24_prep.py" "$TREE" "$RH" >> "$LOG" 2>&1 || fail "rehearsal: c24_prep.py"
  fi
  s0=(--stage S0 --release rC35.0 --stages S0 --genomes "$RH/inputs/genomes-S0.json" --objects "$RH/inputs/objects-S0.json"
      --library pfam 38.2 "$PFAM" gathering --crispr --threads 4 --scan-jobs 2 --scan-cpus 2)
  s1g=(--stage S1 --release rC35.gate --stages S0,S1,PLX --genomes "$RH/inputs/genomes-S1.json"
       --objects "$RH/inputs/objects-S1.json" --library pfam 38.2 "$PFAM" gathering --crispr --threads 4 --scan-jobs 2
       --scan-cpus 2)
  s1=(--stage S1 --release rC35.rh --stages S0,S1 --genomes "$RH/inputs/genomes-S1.json"
      --objects "$RH/inputs/objects-S1.json" --library pfam 38.2 "$PFAM" gathering --crispr --threads 4 --scan-jobs 2
      --scan-cpus 2 --cluster-parts 2)
  # 1. stage S0 and the earlier release rC35.0, whole
  if [ ! -e "$RH/root/runs/stage=rC35.0/catalog.json" ]; then
    launch "$RH/root" "$RH/profile.yaml" c35-rh-s0 16G 8 "${s0[@]}"
    rh_wait c35-rh-s0
    [ "$(exit_of c35-rh-s0)" = 0 ] || fail "rehearsal: c35-rh-s0 exit '$(exit_of c35-rh-s0)'"
  fi
  verify "$RH/root" rC35.0 "$RH/check" || fail "rehearsal: rC35.0 does not verify"
  cp "$RH/check/verify-rC35.0.log" "$RH/check/verify-rC35.0-before.log"
  # 2. stage S1's steps alone: the cluster gate of a release naming a stage that is not built ends the launch
  if [ ! -e "$RH/root/runs/stage=S1/crispr-$NEAR.json" ]; then
    launch "$RH/root" "$RH/profile.yaml" c35-rh-s1 16G 8 "${s1g[@]}"
    rh_wait c35-rh-s1
    rc=$(exit_of c35-rh-s1)
    [ "$rc" = 2 ] && grep -q "refused before rC35.gate.cluster: .*PLX" "$LOGS/c35-rh-s1.log" \
      || fail "rehearsal: c35-rh-s1 exit '$rc', not the cluster gate's refusal naming PLX"
  fi
  # 3. the real command's form: every stage step skipped, the release in parts, the guard and the rate check
  IMG=$IMAGE bb py "$BIN/c35_calib.py" --root "$RH/root" --release rC35.0 --out "$RH/check/calib.json" > /dev/null 2>> "$LOG" \
    || fail "rehearsal: c35_calib.py"
  bb py "$BIN/c35_est.py" before --root "$RH/root" --release rC35.rh --stages S0,S1 --earlier rC35.0 \
    --earlier-stages S0 --calib "$RH/check/calib.json" --parts 2 --profile "$RH/profile.yaml" --jobs 2 --cpus 2 \
    --tmp "$RH/tmp/est" > "$est" 2>> "$LOG" || fail "rehearsal: c35_est.py before"
  log "rehearsal: estimate: $(python3 -I -c 'import json,sys; d=json.load(open(sys.argv[1])); print(json.dumps({"measured": {k: v for k, v in d["measured"].items() if k != "per_stage"}, "code": d["code"]}))' "$est" | cut -c1-900)"
  unit=c35-rh-1
  if [ ! -e "$RH/state/OK-$unit" ]; then
    remove_our_stop "$RH/root" || fail "rehearsal: a STOP not ours in $RH/root"
    listing "$RH/root" "$RH/check/before-listing.tsv"
    mkdir -p "$RH/state/$unit"
    touch "$RH/state/$unit/since"
    launch "$RH/root" "$RH/profile.yaml" "$unit" 16G 8 "${s1[@]}"
    watch "$RH/root" rC35.rh "$unit" "$RH/state/$unit" "$est" c35-rh- "$RH/tmp/scratch"
    [ $? = 3 ] && fail "rehearsal: a STOP not ours appeared during $unit"
    ended_release "$RH/root" "$unit" S1 rC35.rh || fail "rehearsal: $unit did not end as expected (read $LOGS/$unit.log)"
    n=$(find "$RH/root/runs/stage=rC35.rh" -name 'cluster.*.part-*.json' | wc -l)
    log "rehearsal: $unit: $n clustering part records"
    [ "$n" -ge 2 ] || fail "rehearsal: the release was not clustered in parts"
    listing "$RH/root" "$RH/check/after-listing.tsv"
    touch "$RH/state/OK-$unit"
  fi
  verify "$RH/root" rC35.rh "$RH/check" || fail "rehearsal: rC35.rh does not verify"
  verify "$RH/root" rC35.0 "$RH/check" || fail "rehearsal: rC35.0 does not verify after rC35.rh"
  [ "$(hash_of "$RH/check/verify-rC35.0.log")" = "$(hash_of "$RH/check/verify-rC35.0-before.log")" ] \
    || fail "rehearsal: rC35.0's manifest hash changed"
  # 4. the pack in parts, then verified
  rh_pack || fail "rehearsal: the pack"
  # 5. the release check, as far as the fixture allows
  rh_examples
  # 6. the same launch again: every step skipped, nothing changed
  unit=c35-rh-2
  if [ ! -e "$RH/state/OK-$unit" ]; then
    listing "$RH/root" "$RH/check/again-before.tsv"
    mkdir -p "$RH/state/$unit"
    touch "$RH/state/$unit/since"
    launch "$RH/root" "$RH/profile.yaml" "$unit" 16G 8 "${s1[@]}"
    watch "$RH/root" rC35.rh "$unit" "$RH/state/$unit" "$est" c35-rh- "$RH/tmp/scratch"
    rc=$(exit_of "$unit")
    log "rehearsal: $unit exit '$rc': $(grep -c 'skipped: run record and outputs present' "$LOGS/$unit.log") steps skipped, $(grep -c ' published ' "$LOGS/$unit.log") published"
    [ "$rc" = 0 ] || fail "rehearsal: the second launch exit '$rc'"
    grep -q ' published ' "$LOGS/$unit.log" && fail "rehearsal: the second launch published a step"
    grep -q "all $NSTEPS steps done" "$LOGS/$unit.log" || fail "rehearsal: the second launch did not say all $NSTEPS steps done"
    listing "$RH/root" "$RH/check/again-after.tsv"
    cmp -s "$RH/check/again-before.tsv" "$RH/check/again-after.tsv" || fail "rehearsal: the second launch changed the root"
    touch "$RH/state/OK-$unit"
  fi
  # 7. the checks of the real build, on the rehearsal
  bb py "$BIN/c35_check.py" --root "$RH/root" --before "$RH/check/before-listing.tsv" \
    --after "$RH/check/after-listing.tsv" --release rC35.rh --stages S0,S1 --earlier rC35.0 --earlier-stages S0 \
    --est "$est" --out "$RH/check/out" --tmp "$RH/tmp/check" > "$RH/check/check.out" 2>&1
  rc=$?
  log "rehearsal: c35_check.py exit $rc: $(tail -n 1 "$RH/check/check.out" | cut -c1-600)"
  [ "$rc" = 0 ] || fail "rehearsal: c35_check.py failed (read $RH/check/check.out)"
  touch "$RH/OK"
  log "rehearsal: OK"
}

rh_pack() {
  local unit=c35-rh-pack rc
  [ -e "$RH/state/OK-pack" ] && return 0
  bash "$TREE/scripts/pod_run.sh" --profile "$RH/profile.yaml" --image "$IMAGE" --digest "$DIGEST" --mem 16G --threads 4 \
    --bind /tmp/if2453 --bind "$ROOT:$ROOT:ro" --bind "$PFAMDIR:$PFAMDIR:ro" --env "PYTHONPATH=$PYPATH" \
    --log "$LOGS/$unit.log" "$unit" pack build --root "$RH/root" --release rC35.rh --model "pfam-38.2=$PFAM" \
    --profile "$RH/profile.yaml" --parts 2 --jobs 2 >> "$LOG" 2>&1
  sleep 2
  rh_wait "$unit"
  rc=$(exit_of "$unit")
  log "rehearsal: $unit exit '$rc': $(grep -E 'estimate|written|pack ' "$LOGS/$unit.log" | tail -n 2 | tr '\n' ' ' | cut -c1-500)"
  [ "$rc" = 0 ] || return 1
  bb cli pack verify --root "$RH/root" --release rC35.rh --profile "$RH/profile.yaml" > "$RH/check/pack-verify.log" 2>&1
  rc=$?
  log "rehearsal: pack verify exit $rc: $(tail -n 1 "$RH/check/pack-verify.log")"
  [ "$rc" = 0 ] || return 1
  touch "$RH/state/OK-pack"
}

rh_examples() {
  local rc
  [ -e "$RH/state/OK-examples" ] && return 0
  rm -rf "$RH/examples"
  TREE=$TREE PROFILE=$RH/profile.yaml IMAGE=$IMAGE DIGEST=$DIGEST ROOT=$RH/root OUT=$RH/examples UNIT=c35-rh-ex \
    bash "$TREE/scripts/release_examples.sh" rC35.rh > "$RH/check/examples.log" 2>&1
  rc=$?
  log "rehearsal: release_examples.sh rC35.rh exit $rc (the fixture holds none of the controls: 1 expected): $(tail -n 1 "$RH/check/examples.log" | cut -c1-400)"
  [ -s "$RH/examples/summary.tsv" ] || fail "rehearsal: the release check wrote no summary.tsv"
  log "rehearsal: summary.tsv: $(wc -l < "$RH/examples/summary.tsv") rows; verdicts: $(cut -f6 "$RH/examples/summary.tsv" | sort | uniq -c | tr '\n' ' ')"
  touch "$RH/state/OK-examples"
}

# --- estimate, before, after ----------------------------------------------------------------------------------------

estimate() {
  local d=$C/estimate rc
  mkdir -p "$d"
  profile
  python3 -I "$BIN/c35_calib.py" --root "$ROOT" --release "$EARLIER" --out "$d/calib_$EARLIER.json" > /dev/null \
    || fail "estimate: c35_calib.py"
  bb py "$BIN/c35_est.py" before --root "$ROOT" --release "$REL" --stages "$ALL" --earlier "$EARLIER" \
    --earlier-stages "$EARLIER_STAGES" --calib "$d/calib_$EARLIER.json" --parts "$PARTS" --profile "$PROFILE" \
    --jobs "$J" --cpus "$CPUS" --scratch-peak-gb 148.8 --tmp "$C/tmp/est" > "$d/before.json" 2> "$d/before.err"
  rc=$?
  log "estimate: exit $rc: $(python3 -I -c 'import json,sys; d=json.load(open(sys.argv[1])); p=d["projection"]; print(json.dumps({"measured": {k: v for k, v in d["measured"].items() if k != "per_stage"}, "code": d["code"], "projection": {k: v for k, v in p.items() if k not in ("calls", "steps", "factors")}}))' "$d/before.json" | cut -c1-2500)"
  [ "$rc" = 0 ] || fail "estimate: c35_est.py before (read $d/before.err)"
  python3 -I -c 'import json,sys; p=json.load(open(sys.argv[1]))["projection"]; sys.exit(0 if p["every_part_under_240"] else 1)' "$d/before.json" \
    || fail "estimate: a clustering part is projected at 240 CPU-minutes or more"
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
  [ "$ok" = 1 ] || fail "a published release does not verify before the launch (read $d/verify-*.log)"
  listing "$ROOT" "$d/listing.tsv"
  log "before: $(wc -l < "$d/listing.tsv") files, $(awk -F '\t' '{s += $2} END {printf "%.0f", s}' "$d/listing.tsv") bytes"
  touch "$d/DONE"
}

after() {
  local d=$C/after r ok=1 h0 h1
  mkdir -p "$d"
  [ -e "$ROOT/STOP" ] && log "after: $ROOT/STOP exists: $(head -c 300 "$ROOT/STOP")"
  listing "$ROOT" "$d/listing.tsv"
  log "after: $(wc -l < "$d/listing.tsv") files, $(awk -F '\t' '{s += $2} END {printf "%.0f", s}' "$d/listing.tsv") bytes"
  verify "$ROOT" "$REL" "$d" || ok=0
  for r in $RELEASES; do
    verify "$ROOT" "$r" "$d" || ok=0
    h0=$(hash_of "$C/before/verify-$r.log")
    h1=$(hash_of "$d/verify-$r.log")
    [ -n "$h1" ] && [ "$h0" = "$h1" ] || { log "after: $r's manifest hash '$h1' is not '$h0'"; ok=0; }
  done
  bb cli releases --root "$ROOT" --profile "$REALPROF" > "$d/releases.txt" 2>&1
  log "after: bio-box releases exit $?: $(grep -F "$REL" "$d/releases.txt" | head -n 3 | tr '\n' ' ' | cut -c1-600)"
  [ "$ok" = 1 ] || fail "a release does not verify, or an earlier release's hash changed (read $d/verify-*.log)"
  touch "$d/DONE"
}

# --- run ------------------------------------------------------------------------------------------------------------

next_unit() {
  local n
  n=$(($( [ -e "$C/units" ] && wc -l < "$C/units" || echo 0) + 1))
  echo "c35-rel-$n"
}

run() {
  local unit a free act rc
  [ -e "$C/before/DONE" ] || fail "no $C/before/DONE: run the before mode first"
  [ -e "$RH/OK" ] || fail "no $RH/OK: the rehearsal has not passed"
  [ -s "$C/estimate/before.json" ] || fail "no $C/estimate/before.json: run the estimate mode first"
  [ -e "$C/run/DONE" ] && { log "run: DONE present; nothing to do"; return; }
  [ -e "$C/HALT" ] && { log "run: $C/HALT present; not launching"; exit 0; }
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
  log "free ${free} GB; allowance $a for $T threads; parts $J x $CPUS threads, --cluster-parts $PARTS"
  unit=$(next_unit)
  echo "$unit" >> "$C/units"
  mkdir -p "$C/run/$unit"
  touch "$C/run/$unit/since"
  launch "$ROOT" "$PROFILE" "$unit" 200G "$Q" --stage "$STAGE" --release "$REL" --stages "$ALL" \
    --genomes "$IN/genomes.json" --objects "$IN/objects.json" --library pfam 38.2 "$PFAM" gathering --crispr \
    --threads "$T" --scan-jobs "$J" --scan-cpus "$CPUS" --cluster-parts "$PARTS"
  watch "$ROOT" "$REL" "$unit" "$C/run/$unit" "$C/estimate/before.json" c35-rel- "$SCRATCH"
  rc=$?
  [ "$rc" = 3 ] && fail "a STOP not written by this driver appeared in the root during $unit; it is left in place (read $LOGS/$unit.log)"
  ended_release "$ROOT" "$unit" "$STAGE" "$REL" || fail "$unit did not end as expected (read $LOGS/$unit.log and $LOG)"
  log "run: $REL built; scratch peak $(cat "$C/run/$unit/scratch_peak" 2> /dev/null)"
  touch "$C/run/DONE"
}

# --- pack, examples, check ------------------------------------------------------------------------------------------

pack() {
  local unit rc
  [ -e "$C/after/DONE" ] || fail "no $C/after/DONE: run the after mode first"
  [ -e "$ROOT/STOP" ] && fail "$ROOT/STOP exists"
  mkdir -p "$C/pack"
  unit=c35-pack
  if ! ls "$PACKDIR/release=$REL"/*/meta 2> /dev/null | grep -q .; then
    snap_top
    bash "$TREE/scripts/pod_run.sh" --profile "$REALPROF" --image "$IMAGE" --digest "$DIGEST" --mem 48G \
      --threads "$PACK_JOBS" --bind /tmp/if2453 --bind "$ROOT:$ROOT:ro" --bind "$PFAMDIR:$PFAMDIR:ro" \
      --env "PYTHONPATH=$PYPATH" --log "$LOGS/$unit.log" "$unit" pack build --root "$ROOT" --release "$REL" \
      --model "pfam-38.2=$PFAM" --profile "$REALPROF" --parts "$PACK_PARTS" --jobs "$PACK_JOBS" >> "$LOG" 2>&1
    sleep 2
    while systemctl --user is-active --quiet "$unit"; do
      python3 -I "$BIN/r50_watch.py" --root "$ROOT" --scratch /tmp/if2453/bio-box/scratch/pack --since "$C/after/DONE" \
        --state "$C/pack" --unit-prefix c35-pack --tmp-min-gb "$TMP_MIN_GB" --scratch-max-gb "$SCRATCH_MAX_GB" \
        --cpu-max "$CPU_MAX" --stages "$REL" >> "$C/pack/guard.jsonl"
      [ $? = 3 ] && stop_unit "$unit" "guard: $(tail -n 1 "$C/pack/guard.jsonl")"
      sleep 60
    done
    rc=$(exit_of "$unit")
    log "pack: $unit exit '$rc': $(grep -E 'estimate|written|pack id' "$LOGS/$unit.log" | tail -n 3 | tr '\n' ' ' | cut -c1-700)"
    log "pack: $unit usage: $(grep -E 'pod_run: .* usage' "$LOGS/$unit.log" | tail -n 1)"
    [ "$rc" = 0 ] || fail "pack build exit '$rc' (read $LOGS/$unit.log)"
  fi
  bb cli pack verify --root "$ROOT" --release "$REL" --profile "$REALPROF" > "$C/pack/verify.log" 2>&1
  rc=$?
  log "pack: pack verify exit $rc: $(tail -n 1 "$C/pack/verify.log" | cut -c1-400)"
  [ "$rc" = 0 ] || fail "pack verify exit $rc"
  touch "$C/pack/DONE"
}

examples() {
  local rc
  [ -e "$C/pack/DONE" ] || fail "no $C/pack/DONE: run the pack mode first"
  [ -e "$C/examples/summary.tsv" ] && { log "examples: summary.tsv present; nothing to do"; return; }
  snap_top
  TREE=$TREE PROFILE=$REALPROF IMAGE=$IMAGE DIGEST=$DIGEST ROOT=$ROOT OUT=$C/examples UNIT=c35-ex PFAM=$PFAM \
    bash "$TREE/scripts/release_examples.sh" "$REL" > "$C/examples.log" 2>&1
  rc=$?
  log "examples: release_examples.sh $REL exit $rc: $(tail -n 1 "$C/examples.log" | cut -c1-500)"
  log "examples: verdicts: $(cut -f6 "$C/examples/summary.tsv" | sort | uniq -c | tr '\n' ' ')"
  [ "$rc" = 0 ] || fail "the release check: an invariant failed or is missing (read $C/examples/summary.tsv)"
}

check() {
  local rc
  mkdir -p "$C/check"
  python3 -I "$BIN/c35_est.py" rates --root "$ROOT" --release "$REL" --est "$C/estimate/before.json" \
    --out "$C/check/rates.json" > /dev/null
  bb py "$BIN/c35_check.py" --root "$ROOT" --before "$C/before/listing.tsv" --after "$C/after/listing.tsv" \
    --release "$REL" --stages "$ALL" --earlier "$EARLIER" --earlier-stages "$EARLIER_STAGES" \
    --est "$C/estimate/before.json" --out "$C/check" --tmp "$C/tmp/check" --dec27 > "$C/check/check.out" 2>&1
  rc=$?
  cp "$C/check/root_added.tsv" "$C/root_added.tsv"
  log "check: exit $rc: $(tail -n 1 "$C/check/check.out" | cut -c1-600)"
  return $rc
}

MODE=${1:-}
log "c35 driver: mode '$MODE' (pid $$), tree $TREE, image $IMAGE"
case "$MODE" in
  rehearse) rehearse ;;
  estimate) estimate ;;
  before) before ;;
  run) run ;;
  after) after ;;
  pack) pack ;;
  examples) examples ;;
  check) check ;;
  *) echo "usage: c35_driver.sh rehearse|estimate|before|run|after|pack|examples|check" >&2; exit 2 ;;
esac
log "c35 driver: mode '$MODE' done"
