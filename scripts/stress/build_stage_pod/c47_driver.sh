#!/bin/bash
# Card C-47: build stage C5 (Enterobacterales..Vibrio, 2,032 genomes, 9.39 Gbp) with a full Pfam 38.2 scan
# and both CRISPR producers into /tmp/if2453/r10/root, stage only (`build-stage --stage C5 --depth annotated`,
# omitting --release per C-38a).
#
# Modes:
#   restore   Step 0: verify NO STOP and /tmp > 400 GB free; extract s2c_phase2 tars if absent;
#             recreate /tmp/if2453/s2/inputs/C1..C6 via c47_restore.py + make_inputs.py prok;
#             verify all 12 sha256s against 2026-11-15_s2c_phase2_inputs/inputs/SUMS.
#   rehearse  Stage-only rehearsal (`--depth annotated`, no `--release`) on 20 C5 genomes in a scratch view
#             of P1,P2,PL1 with pfam-38.2 linked; checks exit 0, B-47b copy, c47_check.py, and idempotent relaunch.
#   estimate  Verify C5 inputs sha256 and compute pre-launch CPU estimates (`c47_est.py before`).
#   before    Verify all 7 published releases (`rN.1a..r26.11.1`) and snapshot root listing + earlier sha256.
#   run       Launch `build-stage --stage C5 --depth annotated ... --threads 6 --scan-cpus 4` (omitting --scan-jobs
#             so live runstep.allowance scales automatically), with r50_watch.py (280 CPU-min) + c47_est.py rates
#             (DEC-43 + DEC-53: stop ONLY on judged ratio > 1.20).
#   after     Re-snapshot root listing + earlier sha256, re-verify all 7 published releases, `bio-box releases`.
#   check     Run `c47_check.py` (before vs after, root_added.tsv, steps, contigs, genes, B-47b copy, CRISPR).
#   backup    Write /stor/work/Wilke/if2453/all-the-data/root_backups/added_c5.tar (>= 64 MiB) via dot-temp +
#             dd bs=16M conv=fsync + readback sha256sum + mv -n.
set -u

C=/tmp/if2453/c47
ROOT=/tmp/if2453/r10/root
TREE=${TREE:-/tmp/if2453/ci/C-47}
REALPROF=$TREE/profiles/wilkcomp02.yaml
PROFILE=$C/profile.yaml
IMAGE=/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif
DIGEST=sha256:d67300aa1321b571fa8da79e74dfccdd87a8dca2e56518f324902efa6b9f098a
PYPATH=$TREE/src
PFAMDIR=/stor/scratch/Wilke/if2453/all-the-data/l0/refdb/pfam/38.2
PFAM=$PFAMDIR/Pfam-A.hmm.gz
STOR_S2C=/stor/scratch/Wilke/if2453/all-the-data/provisional/s2c_phase2
TAR1=$STOR_S2C/s2c_phase2_genomes.tar
TAR1_SHA=af3ae7c348328e78e7091e7790ba20e3103a1f0bef21eef1be18d36c2e36fcb9
TAR2=$STOR_S2C/s2c_phase2_rest_genomes.tar
TAR2_SHA=94cc747e86c619f80404941508ca5b8b3cffcf5d7302426f32830bd3db7c19df
GTDB_TAR=/stor/scratch/Wilke/if2453/all-the-data/provisional/s1_inputs/gtdb_r226_metadata.tar
X=/tmp/if2453/s2/phase2
IN_ALL=/tmp/if2453/s2/inputs
IN=$IN_ALL/C5
STAGE=C5
EARLIER=S1,P1,P2,P3,P4,P5,P6,H1,PL1
RELEASES="rN.1a rN.1b rN.1 rN.2 rN.3 rM.2 r26.11.1"
T=${T:-6}
CPUS=${CPUS:-4}
Q=${Q:-40}
CPU_MAX=${CPU_MAX:-280}
TMP_MIN_GB=${TMP_MIN_GB:-200}
START_MIN_GB=${START_MIN_GB:-205}
FLOOR_GB=${FLOOR_GB:-400}
RH=$C/rh
RHSTAGE=C5R
RH_N=${RH_N:-2}
RH_T=${RH_T:-4}
RH_CPUS=${RH_CPUS:-4}
BIN=$(cd "$(dirname "$0")" && pwd)
LOGS=$C/logs
LOG=$C/driver.log
SCRATCH=$C/tmp/scratch
BACKUP_DIR=/stor/work/Wilke/if2453/all-the-data/root_backups
BACKUP_TAR=$BACKUP_DIR/added_c5.tar
NR3=minced-0.4.2-w8-nr3-r23x47-s26x50
NEAR=minced-0.4.2-w8-nr2-r23x47-s26x50-cas1kb-caspfam1-pfam-38.2
C5_GENOMES_SHA=1c65060c80ef120d4550d702a76c4a9c5c3eff17ee25f2694d5a172873468573
C5_OBJECTS_SHA=09d742d70a05528e0bd91b735605b63e0f67ad3758e9f1d02107be30bf5483ad

mkdir -p "$LOGS" "$C/tmp"

log() { echo "$(date '+%F %T.%3N %Z') $*" >> "$LOG"; }
fail() { log "FAILED: $*"; mkdir -p "$C/$MODE"; echo "$*" > "$C/$MODE/FAILED"; exit 1; }

profile() {
  sed "s#^tmp:.*#tmp: $C/tmp#" "$REALPROF" > "$PROFILE"
  log "profile $PROFILE written from $REALPROF"
}

bb() {
  local what=$1
  shift
  case "$what" in
    cli) set -- python -m biobox.cli "$@" ;;
    py) set -- python "$@" ;;
  esac
  systemd-run --user --scope --quiet -p MemoryMax=16G -p MemorySwapMax=0 -p CPUQuota=200% nice -n 19 \
    singularity exec --cleanenv --env PYTHONPATH="$PYPATH" --env OMP_NUM_THREADS=2 --env PYTHONUNBUFFERED=1 \
    --bind /tmp/if2453 --bind "$ROOT:$ROOT:ro" --bind "$PFAMDIR:$PFAMDIR:ro" \
    --bind /stor/scratch/Wilke/if2453/all-the-data:/stor/scratch/Wilke/if2453/all-the-data:ro \
    --bind /sys/fs/cgroup:/sys/fs/cgroup:ro "$IMAGE" "$@"
}

build_units() {
  systemctl --user list-units --type=service --state=active --no-legend --plain \
    'r50-*' 'c24-*' 'c26-*' 'c33-*' 'c34-*' 'c35-*' 'c46-*' 'c47-c5-*' 'c48c-*' | awk '{print $1}' | paste -sd' ' -
}

snap_top() {
  log "uptime: $(uptime)"
  log "df /tmp: $(df -h /tmp | tail -n 1)"
  log "active units: $(systemctl --user list-units --type=service --type=scope --state=active --no-legend --plain | awk '{print $1}' | paste -sd' ' -)"
}

free_gb() { df -B1G --output=avail /tmp | tail -n 1 | tr -d ' '; }
exit_of() { sed -n "s/^pod_run: $1 exited \([0-9]*\)$/\1/p" "$LOGS/$1.log" | tail -n 1; }

verify() {
  bb cli release verify --root "$1" --release "$2" > "$3/verify-$2.log" 2>&1
  local rc=$?
  log "verify $2 ($1): exit $rc: $(grep -E '^release ' "$3/verify-$2.log" | tail -n 1)"
  return $rc
}

snapshot() {
  local d=$1
  mkdir -p "$d"
  find "$ROOT" ! -type d -printf '%P\t%s\t%T@\n' | LC_ALL=C sort > "$d/listing.tsv"
  find "$ROOT" ! -type d -regextype posix-extended -regex ".*/stage=(${EARLIER//,/|})/.*" -printf '%P\n' \
    | LC_ALL=C sort > "$d/earlier.files"
  (cd "$ROOT" && systemd-run --user --scope --quiet -p MemoryMax=2G -p MemorySwapMax=0 nice -n 19 ionice -c 3 \
    xargs -d '\n' sha256sum < "$d/earlier.files") > "$d/earlier.sha256"
  log "snapshot $d: $(wc -l < "$d/listing.tsv") files, $(wc -l < "$d/earlier.sha256") of $EARLIER hashed, $(awk -F '\t' '{s += $2} END {printf "%.0f", s}' "$d/listing.tsv") bytes"
}

restore() {
  local free need sha k s
  [ -e "$ROOT/STOP" ] && fail "$ROOT/STOP exists"
  free=$(df -Pk /tmp | awk 'NR==2 {print $4}')
  log "restore: /tmp free ${free} KiB"
  if [ ! -d "$X" ]; then
    need=$(( ($(stat -c %s "$TAR1") + $(stat -c %s "$TAR2")) / 1024 + FLOOR_GB * 1000 * 1000 * 1000 / 1024 ))
    [ "$free" -gt "$need" ] || fail "/tmp has $free KiB free; need $need KiB"
    for t in "$TAR1:$TAR1_SHA" "$TAR2:$TAR2_SHA"; do
      sha=$(sha256sum "${t%%:*}" | cut -d' ' -f1)
      log "restore: sha256 ${t%%:*}: $sha"
      [ "$sha" = "${t##*:}" ] || fail "sha256 mismatch on ${t%%:*}"
    done
    rm -rf "$X.tmp"
    mkdir -p "$X.tmp"
    for t in "$TAR1" "$TAR2"; do
      [ -e "$ROOT/STOP" ] && fail "$ROOT/STOP appeared"
      tar -C "$X.tmp" --no-same-owner --no-same-permissions -xf "$t" || fail "tar -xf $t"
      log "restore: extracted $t ($(find "$X.tmp" -type f | wc -l) files)"
    done
    mv "$X.tmp" "$X"
    log "restore: $X ready ($(du -sb "$X" | cut -f1) bytes)"
  else
    log "restore: $X already exists"
  fi
  mkdir -p "$C/stages" "$C/inputs_log"
  if [ ! -f "$C/stages/C5/lock.json" ]; then
    bb py "$BIN/c47_restore.py" restore \
      "$TREE/planning/measurements/2026-11-15_s2c_phase2_inputs/stages/stages.tsv" \
      "$GTDB_TAR" "$X" "$C/stages" > "$C/stages/restore.json" 2> "$C/stages/restore.err" \
      || fail "c47_restore.py restore failed"
    log "restore: c47_restore.py ok"
  fi
  for k in 1 2 3 4 5 6; do
    s="C$k"
    if [ ! -f "$IN_ALL/$s/genomes.json" ]; then
      rm -rf "$IN_ALL/$s.tmp"
      bb py "$BIN/make_inputs.py" prok "$C/stages/$s/lock.json" "$C/stages/$s/selection.tsv" "$IN_ALL/$s.tmp" \
        --source-release datasets-v2-2026-11-13 > "$C/inputs_log/$s.json" 2> "$C/inputs_log/$s.err" \
        || fail "make_inputs.py prok failed for $s"
      mv "$IN_ALL/$s.tmp" "$IN_ALL/$s"
      log "restore: built $IN_ALL/$s"
    fi
  done
  (cd "$IN_ALL" && sha256sum -c "$TREE/planning/measurements/2026-11-15_s2c_phase2_inputs/inputs/SUMS") \
    > "$C/inputs_log/SUMS.check" 2>&1 || fail "SUMS mismatch: $(cat "$C/inputs_log/SUMS.check")"
  log "restore: all 12 C1..C6 input sha256s match SUMS byte-for-byte"
  touch "$C/restore.DONE"
}

rh_listing() {
  find "$RH/root" ! -type d ! -name STOP -printf '%P\t%s\t%T@\n' | LC_ALL=C sort > "$1"
}

rehearse() {
  local unit rc n
  [ -e "$RH/OK" ] && { log "rehearsal: OK present"; return; }
  [ -f "$C/stages/C5/lock.json" ] || fail "rehearsal: run restore first"
  mkdir -p "$RH/sample" "$RH/inputs" "$RH/tmp" "$RH/durable" "$RH/bulk" "$RH/check" "$RH/run"
  if [ ! -e "$RH/root" ]; then
    python3 -I "$BIN/c26_view.py" "$ROOT" "$RH/root" "$EARLIER" pfam-38.2 >> "$LOG" 2>&1 || fail "rehearsal: view"
  fi
  if [ ! -e "$RH/inputs/genomes.json" ]; then
    bb py "$BIN/c47_restore.py" sample "$C/stages/C5" "$RH_N" "$RH/sample" > "$RH/sample/report.json" 2>> "$LOG" \
      || fail "rehearsal: sample"
    bb py "$BIN/make_inputs.py" prok "$RH/sample/lock.json" "$RH/sample/selection.tsv" "$RH/inputs" \
      --source-release datasets-v2-2026-11-13 > "$RH/inputs/report.json" 2>> "$LOG" \
      || fail "rehearsal: make_inputs.py"
    log "rehearsal: sample inputs: $(tr -d '\n' < "$RH/inputs/report.json" | cut -c1-400)"
  fi
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
  bb py "$BIN/c47_est.py" before --root "$RH/root" --profile "$RH/profile.yaml" --inputs "$RH/inputs" \
    --stage "$RHSTAGE" > "$RH/check/before.json" 2>> "$LOG" || fail "rehearsal: before est"
  local args=(
    --stage "$RHSTAGE" --depth annotated
    --genomes "$RH/inputs/genomes.json" --objects "$RH/inputs/objects.json"
    --library pfam 38.2 "$PFAM" gathering --crispr
    --threads "$RH_T" --scan-cpus "$RH_CPUS"
  )
  unit=c47-rh-1
  if [ ! -e "$RH/run/OK-$unit" ]; then
    rh_listing "$RH/check/before-listing.tsv"
    bash "$TREE/scripts/pod_run.sh" --profile "$RH/profile.yaml" --image "$IMAGE" --digest "$DIGEST" --mem 32G \
      --threads "$RH_T" --bind /tmp/if2453 --bind "$ROOT:$ROOT:ro" --bind "$PFAMDIR:$PFAMDIR:ro" \
      --env "PYTHONPATH=$PYPATH" --log "$LOGS/$unit.log" "$unit" \
      build-stage --profile "$RH/profile.yaml" --root "$RH/root" "${args[@]}" >> "$LOG" 2>&1
    while systemctl --user is-active --quiet "$unit"; do sleep 5; done
    rc=$(exit_of "$unit")
    log "rehearsal: $unit exit '$rc'"
    [ "$rc" = 0 ] || fail "rehearsal: $unit exited '$rc', expected 0 (stage-only --depth annotated)"
    for step in inputs ingest genecall "crispr-$NR3" igr tapes domains-pfam-38.2 domains-pfam-38.2.plan "crispr-$NEAR"; do
      [ -e "$RH/root/runs/stage=$RHSTAGE/$step.json" ] || fail "rehearsal: missing run record $step"
    done
    touch "$RH/run/OK-$unit"
  fi
  python3 -I -c 'import json,sys; c=json.load(open(sys.argv[1]))["copy"]; print(json.dumps({k: c[k] for k in ("proteins","residues","rows","proteins_with_hit")} | {"sources": [s["stage"] for s in c["sources"]]})); sys.exit(0 if c["proteins"] > 0 and c["sources"] else 1)' \
    "$RH/root/l1/domain_scan/library=pfam-38.2/stage=$RHSTAGE/plan.json" > "$RH/check/copy.json" \
    || fail "rehearsal: scan plan copied nothing"
  log "rehearsal: copy: $(cat "$RH/check/copy.json")"
  mkdir -p "$RH/check/before" "$RH/check/after"
  cp "$RH/check/before-listing.tsv" "$RH/check/before/listing.tsv"
  rh_listing "$RH/check/after/listing.tsv"
  bb py "$BIN/c47_check.py" --root "$RH/root" --before "$RH/check/before" --after "$RH/check/after" \
    --log "$LOGS/c47-rh-1.log" --objects "$RH/inputs/objects.json" --out "$RH/check/out" --tmp "$RH/tmp/check" \
    --stage "$RHSTAGE" --earlier "$EARLIER" --before-est "$RH/check/before.json" > "$RH/check/check.out" 2>&1 \
    || fail "rehearsal: c47_check.py failed (read $RH/check/check.out)"
  log "rehearsal: c47_check.py passed"
  unit=c47-rh-2
  if [ ! -e "$RH/run/OK-$unit" ]; then
    rh_listing "$RH/check/again-before.tsv"
    bash "$TREE/scripts/pod_run.sh" --profile "$RH/profile.yaml" --image "$IMAGE" --digest "$DIGEST" --mem 32G \
      --threads "$RH_T" --bind /tmp/if2453 --bind "$ROOT:$ROOT:ro" --bind "$PFAMDIR:$PFAMDIR:ro" \
      --env "PYTHONPATH=$PYPATH" --log "$LOGS/$unit.log" "$unit" \
      build-stage --profile "$RH/profile.yaml" --root "$RH/root" "${args[@]}" >> "$LOG" 2>&1
    while systemctl --user is-active --quiet "$unit"; do sleep 5; done
    rc=$(exit_of "$unit")
    [ "$rc" = 0 ] || fail "rehearsal: $unit exited '$rc'"
    grep -q ' published ' "$LOGS/$unit.log" && fail "rehearsal: relaunch published a step"
    rh_listing "$RH/check/again-after.tsv"
    cmp -s "$RH/check/again-before.tsv" "$RH/check/again-after.tsv" || fail "rehearsal: relaunch changed view"
    log "rehearsal: $unit skipped all steps and left view unchanged"
    touch "$RH/run/OK-$unit"
  fi
  touch "$RH/OK"
  log "rehearsal: OK"
}

estimate() {
  local d=$C/estimate rc
  mkdir -p "$d"
  profile
  printf "%s  genomes.json\n%s  objects.json\n" "$C5_GENOMES_SHA" "$C5_OBJECTS_SHA" > "$d/SUMS.expected"
  (cd "$IN" && sha256sum -c "$d/SUMS.expected") > "$d/SUMS.check" 2>&1 \
    || fail "estimate: C5 inputs sha256 mismatch"
  bb py "$BIN/c47_est.py" before --root "$ROOT" --profile "$PROFILE" --inputs "$IN" --stage "$STAGE" \
    --out-plan "$d/genecall_plan.json" > "$d/before.json" 2> "$d/before.err"
  rc=$?
  log "estimate: exit $rc: $(tr -d '\n ' < "$d/before.json" | cut -c1-1200)"
  [ "$rc" = 0 ] || fail "estimate: c47_est.py before failed"
}

before() {
  local d=$C/before r ok=1 act
  [ -e "$d/DONE" ] && { log "before: DONE present"; return; }
  [ -e "$ROOT/STOP" ] && fail "$ROOT/STOP exists"
  act=$(build_units)
  [ -z "$act" ] || fail "a build unit is active: $act"
  mkdir -p "$d"
  snap_top
  for r in $RELEASES; do verify "$ROOT" "$r" "$d" || ok=0; done
  [ "$ok" = 1 ] || fail "a published release failed verify before launch"
  snapshot "$d"
  touch "$d/DONE"
}

next_unit() {
  local n
  n=$(($( [ -e "$C/units" ] && wc -l < "$C/units" || echo 0) + 1))
  echo "c47-c5-$n"
}

stage_records() { find "$ROOT/runs/stage=$STAGE" -maxdepth 1 -name '*.json' ! -name '*.publishing.json' 2> /dev/null | wc -l; }

stop_unit() {
  log "stopping $1: $2"
  systemctl --user stop "$1"
  while systemctl --user is-active --quiet "$1"; do sleep 2; done
  fail "$2 (unit $1 stopped)"
}

run() {
  local unit free act rc line nrec=-1 n k=0
  [ -e "$C/before/DONE" ] || fail "run: no $C/before/DONE"
  [ -e "$RH/OK" ] || fail "run: no $RH/OK"
  [ -s "$C/estimate/before.json" ] || fail "run: no $C/estimate/before.json"
  [ -e "$C/run/DONE" ] && { log "run: DONE present"; return; }
  mkdir -p "$C/run"
  profile
  [ -e "$ROOT/STOP" ] && fail "$ROOT/STOP exists"
  act=$(build_units)
  [ -z "$act" ] || fail "a build unit is active: $act"
  free=$(free_gb)
  [ "$free" -ge "$START_MIN_GB" ] || fail "/tmp has $free GB free, under $START_MIN_GB"
  snap_top
  unit=$(next_unit)
  echo "$unit" >> "$C/units"
  touch "$C/run/since"
  log "launch $unit (pod_run --threads $Q --mem 200G): build-stage --stage $STAGE --depth annotated --threads $T --scan-cpus $CPUS (no --release, no --scan-jobs)"
  bash "$TREE/scripts/pod_run.sh" --profile "$PROFILE" --image "$IMAGE" --digest "$DIGEST" --mem 200G \
    --threads "$Q" --bind /tmp/if2453 --bind "$PFAMDIR:$PFAMDIR:ro" --env "PYTHONPATH=$PYPATH" \
    --log "$LOGS/$unit.log" "$unit" \
    build-stage --profile "$PROFILE" --root "$ROOT" --stage "$STAGE" --depth annotated \
    --genomes "$IN/genomes.json" --objects "$IN/objects.json" \
    --library pfam 38.2 "$PFAM" gathering --crispr \
    --threads "$T" --scan-cpus "$CPUS" >> "$LOG" 2>&1
  sleep 10
  while systemctl --user is-active --quiet "$unit"; do
    if [ $((k % 6)) = 0 ]; then
      line=$(python3 -I "$BIN/r50_watch.py" --root "$ROOT" --scratch "$SCRATCH" --since "$C/run/since" \
        --state "$C/run" --unit-prefix c47-c5- --tmp-min-gb "$TMP_MIN_GB" --scratch-max-gb 100 \
        --cpu-max "$CPU_MAX" --stages "$STAGE")
      rc=$?
      echo "$line" >> "$C/run/guard.jsonl"
      [ "$rc" = 3 ] && stop_unit "$unit" "guard: $line"
      [ $((k % 60)) = 0 ] && log "guard: $(echo "$line" | cut -c1-400); log: $(grep -E 'run-step|build-stage' "$LOGS/$unit.log" | tail -n 1 | cut -c1-240)"
    fi
    n=$(stage_records)
    if [ "$n" != "$nrec" ]; then
      nrec=$n
      line=$(bb py "$BIN/c47_est.py" rates --root "$ROOT" --objects "$IN/objects.json" --stage "$STAGE" \
        --before "$C/estimate/before.json" 2>> "$LOG")
      rc=$?
      echo "$line" >> "$C/run/rates.jsonl"
      log "rates ($n records of $STAGE): $(echo "$line" | cut -c1-500)"
      [ "$rc" = 3 ] && stop_unit "$unit" "DEC-53 stop (judged step > +20%): $(echo "$line" | cut -c1-900)"
      [ "$rc" = 4 ] && stop_unit "$unit" "DEC-53 stop (genecall batch > +20%): $(echo "$line" | cut -c1-900)"
    fi
    k=$((k + 1))
    sleep 10
  done
  line=$(bb py "$BIN/c47_est.py" rates --root "$ROOT" --objects "$IN/objects.json" --stage "$STAGE" \
    --before "$C/estimate/before.json" 2>> "$LOG")
  echo "$line" >> "$C/run/rates.jsonl"
  log "rates at end: $(echo "$line" | cut -c1-600)"
  rc=$(exit_of "$unit")
  log "$unit ended, exit '$rc'; usage: $(grep -E 'pod_run: .* usage' "$LOGS/$unit.log" | tail -n 1)"
  [ "$rc" = 0 ] || fail "$unit exited '$rc', expected 0"
  for step in inputs ingest genecall "crispr-$NR3" igr tapes domains-pfam-38.2 domains-pfam-38.2.plan "crispr-$NEAR"; do
    [ -e "$ROOT/runs/stage=$STAGE/$step.json" ] || fail "missing run record $STAGE/$step"
  done
  touch "$C/run/DONE"
}

after() {
  local d=$C/after r ok=1
  mkdir -p "$d"
  snapshot "$d"
  for r in $RELEASES; do verify "$ROOT" "$r" "$d" || ok=0; done
  bb cli releases --root "$ROOT" > "$d/releases.txt" 2>&1
  log "after: bio-box releases exit $?; verify_ok=$ok"
  [ "$ok" = 1 ] || fail "after: a published release failed verify"
  touch "$d/DONE"
}

check() {
  local unit rc
  mkdir -p "$C/check"
  unit=$(tail -n 1 "$C/units")
  bb py "$BIN/c47_check.py" --root "$ROOT" --before "$C/before" --after "$C/after" --log "$LOGS/$unit.log" \
    --objects "$IN/objects.json" --out "$C/check" --tmp "$C/tmp/check" --stage "$STAGE" \
    --earlier "$EARLIER" --before-est "$C/estimate/before.json" > "$C/check/check.out" 2>&1
  rc=$?
  cp "$C/check/root_added.tsv" "$C/root_added.tsv"
  log "check: exit $rc: $(tail -n 1 "$C/check/check.out" | cut -c1-600)"
  return $rc
}

backup() {
  local local_tar="$C/added_c5.tar" tmp_tar="$BACKUP_DIR/.added_c5.tar.tmp" sz sha1 sha2
  [ -s "$C/root_added.tsv" ] || fail "backup: missing $C/root_added.tsv"
  [ ! -e "$BACKUP_TAR" ] || fail "backup: $BACKUP_TAR already exists"
  awk -F'\t' 'NR > 1 {print $1}' "$C/root_added.tsv" > "$C/added_paths.txt"
  tar -C "$ROOT" -cf "$local_tar" -T "$C/added_paths.txt" || fail "backup: tar failed"
  sz=$(stat -c %s "$local_tar")
  [ "$sz" -ge 67108864 ] || fail "backup: $local_tar is $sz bytes (< 64 MiB)"
  sha1=$(sha256sum "$local_tar" | cut -d' ' -f1)
  dd if="$local_tar" of="$tmp_tar" bs=16M conv=fsync 2>> "$LOG" || fail "backup: dd failed"
  sha2=$(sha256sum "$tmp_tar" | cut -d' ' -f1)
  [ "$sha1" = "$sha2" ] || { rm -f "$tmp_tar"; fail "backup: readback sha256 $sha2 != $sha1"; }
  mv -n "$tmp_tar" "$BACKUP_TAR" || fail "backup: mv -n failed"
  [ -f "$BACKUP_TAR" ] || fail "backup: $BACKUP_TAR missing after mv"
  printf "%s\t%s\t%s\n" "$BACKUP_TAR" "$sz" "$sha2" > "$C/backup.tsv"
  log "backup: $BACKUP_TAR written ($sz bytes, sha256 $sha2)"
  rm -f "$local_tar"
}

MODE=${1:-}
log "c47 driver: mode '$MODE' (pid $$), tree $TREE, image $IMAGE"
case "$MODE" in
  restore) restore ;;
  rehearse) rehearse ;;
  estimate) estimate ;;
  before) before ;;
  run) run ;;
  after) after ;;
  check) check ;;
  backup) backup ;;
  *) echo "usage: c47_driver.sh restore|rehearse|estimate|before|run|after|check|backup" >&2; exit 2 ;;
esac
log "c47 driver: mode '$MODE' done"
