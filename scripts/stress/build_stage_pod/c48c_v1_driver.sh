#!/bin/bash
# Card C-48c part 1: build stage V1 (UHGV v1.0 vOTU representatives: 168,536 genomes, one contig each, kind virus,
# gene calling gv-meta; card C-39's inputs) into /tmp/if2453/r10/root as a stage only: `build-stage --stage V1
# --depth annotated` with Pfam 38.2 in full (card B-47b's copies of rows earlier stages measured stay on) and
# --crispr, no --release (card C-38a). In the form of card C-47's c47_driver.sh.
#
# Stage C5 is parked half built in the root (its inputs and ingest published; card C-47). Nothing of C5 is written;
# its files are hashed before and after with the other stages', and the launch is refused when a contig of V1 is
# one C5 owns (`c48c_v1_est.py before`, exit 5).
#
# Modes (each in a small detached unit of its own, `systemd-run --user --unit=c48c-drv-<mode> --collect`):
#   estimate  Verify V1's six input files against card C-39's sha256; `c48c_v1_est.py before`: the code's
#             estimates, the contigs each stage owns (exit 5 if C5 owns one), the comparison with C-39's §5 (exit 6
#             on a difference over 20% the code does not explain).
#   rehearse  Stage-only rehearsal in a scratch view (S1,P1..P6,H1,PL1 with their Pfam 38.2 scans, and C5's
#             l0seq tables and runs, linked; the root bound read-only) on a sample of V1 (`c48c_v1_est.py sample`:
#             every 1000th genome, 12 circular, the first records S1 and H1 own): exit 0, every step recorded, the
#             B-47b copy, `c48c_v1_check.py --no-expect`, and a relaunch that publishes nothing and leaves the view
#             unchanged.
#   before    Verify the seven published releases; snapshot the root's listing and the sha256 of every file of
#             S1,P1..P6,H1,PL1,C5.
#   run       Launch `build-stage --stage V1 --depth annotated ... --threads 6 --scan-cpus 4` (no --release, no
#             --scan-jobs: the allowance of the hour, less ours elsewhere, sizes the shards) through pod_run.sh at a
#             quota of 40; watch r50_watch.py (a process at 280 CPU-minutes, /tmp under 200 GB free, scratch over
#             100 GB) and `c48c_v1_est.py rates` (DEC-43, DEC-53: stop only when a judged step is more than 20% over;
#             under is recorded). A stop stays a stop: there is no accepted-record list.
#   after     Snapshot again, verify the seven releases again, `bio-box releases`.
#   check     `c48c_v1_check.py`: before against after (only stage=V1 files added, nothing earlier changed), steps,
#             contigs, genes, B-47b copies, CRISPR arrays of both producers, C-39's counts; root_added.tsv.
#   backup    One uncompressed tar of exactly the added files to
#             /stor/work/Wilke/if2453/all-the-data/root_backups/added_v1.tar (>= 64 MiB): dot-temporary name,
#             dd bs=16M conv=fsync, read back, sha256 equal, mv -n. Refused when added_v1.tar exists.
set -u

C=/tmp/if2453/c48c
ROOT=/tmp/if2453/r10/root
TREE=${TREE:-/tmp/if2453/ci/C-48c-1}
REALPROF=$TREE/profiles/wilkcomp02.yaml
PROFILE=$C/profile.yaml
IMAGE=/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif
DIGEST=sha256:d67300aa1321b571fa8da79e74dfccdd87a8dca2e56518f324902efa6b9f098a
PYPATH=$TREE/src
PFAMDIR=/stor/scratch/Wilke/if2453/all-the-data/l0/refdb/pfam/38.2
PFAM=$PFAMDIR/Pfam-A.hmm.gz
IN=/tmp/if2453/s2/inputs/V1
STAGE=V1
EARLIER=S1,P1,P2,P3,P4,P5,P6,H1,PL1
PARKED=C5
HASHED=$EARLIER,$PARKED
RELEASES="rN.1a rN.1b rN.1 rN.2 rN.3 rM.2 r26.11.1"
T=${T:-6}
CPUS=${CPUS:-4}
Q=${Q:-40}
CPU_MAX=${CPU_MAX:-280}
TMP_MIN_GB=${TMP_MIN_GB:-200}
START_MIN_GB=${START_MIN_GB:-205}
RH=$C/rh
RHSTAGE=V1R
RH_T=${RH_T:-4}
RH_CPUS=${RH_CPUS:-4}
BIN=$(cd "$(dirname "$0")" && pwd)
LOGS=$C/logs
LOG=$C/driver.log
SCRATCH=$C/tmp/scratch
BACKUP_DIR=/stor/work/Wilke/if2453/all-the-data/root_backups
BACKUP_TAR=$BACKUP_DIR/added_v1.tar
NR3=minced-0.4.2-w8-nr3-r23x47-s26x50
NEAR=minced-0.4.2-w8-nr2-r23x47-s26x50-cas1kb-caspfam1-pfam-38.2
STEPS="inputs ingest genecall crispr-$NR3 igr tapes domains-pfam-38.2 domains-pfam-38.2.plan crispr-$NEAR"
# card C-39's note, table of §3 (measured, V1/inputs_SUMS)
V1_SUMS="752ac47ff81cca4432d88574a0183427fe76cb2d1b7148e90a0e4ce85c398240  genomes.json
dc7bed72ef15ccb93f41da05a8d600f8b05cad9312993fe2723c2e5761864d6c  objects.json
31ef5a640048d4bae72e947c380fb00ac0268f057114250c613199155e87a7cf  records.tsv
4142a62fae1678e3272cf972ff05eaa52c5ddaa4ea49e196d68529365e109cf6  repeated_records.tsv
96216d41981855ecc6d1aad0094628d6421844180a5c4e96c757374bd9c8634d  shared_records.tsv
9c595d262d89100f25ebe40964ee31912b26d9b64c6b5f453daf73b931071356  uhgv_v1.0_votus_full.fna"

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
    'r50-*' 'c24-*' 'c26-*' 'c33-*' 'c34-*' 'c35-*' 'c46-*' 'c47-c5-*' 'c48c-v1-*' | awk '{print $1}' | paste -sd' ' -
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
  find "$ROOT" ! -type d -regextype posix-extended -regex ".*/stage=(${HASHED//,/|})/.*" -printf '%P\n' \
    | LC_ALL=C sort > "$d/earlier.files"
  (cd "$ROOT" && systemd-run --user --scope --quiet -p MemoryMax=2G -p MemorySwapMax=0 nice -n 19 ionice -c 3 \
    xargs -d '\n' sha256sum < "$d/earlier.files") > "$d/earlier.sha256"
  log "snapshot $d: $(wc -l < "$d/listing.tsv") files, $(wc -l < "$d/earlier.sha256") of $HASHED hashed, $(awk -F '\t' '{s += $2} END {printf "%.0f", s}' "$d/listing.tsv") bytes"
}

estimate() {
  local d=$C/estimate rc
  mkdir -p "$d"
  [ -e "$ROOT/STOP" ] && fail "$ROOT/STOP exists"
  profile
  printf '%s\n' "$V1_SUMS" > "$d/SUMS.expected"
  (cd "$IN" && sha256sum -c "$d/SUMS.expected") > "$d/SUMS.check" 2>&1 \
    || fail "estimate: V1 inputs differ from card C-39's sha256: $(tr '\n' ' ' < "$d/SUMS.check")"
  log "estimate: the six V1 input files match card C-39's sha256"
  bb py "$BIN/c48c_v1_est.py" before --root "$ROOT" --profile "$PROFILE" --inputs "$IN" --stage "$STAGE" \
    --forbid "$PARKED" --out-plan "$d/genecall_plan.json" > "$d/before.json" 2> "$d/before.err"
  rc=$?
  log "estimate: exit $rc: $(tr -d '\n ' < "$d/before.json" | cut -c1-1500) $(tail -n 2 "$d/before.err")"
  [ "$rc" = 5 ] && fail "estimate: a contig of V1 is owned by $PARKED (do not launch): $(tail -n 1 "$d/before.err")"
  [ "$rc" = 6 ] && fail "estimate: an estimate is more than 20% from C-39's without an explanation"
  [ "$rc" = 0 ] || fail "estimate: c48c_v1_est.py before failed (exit $rc)"
  touch "$d/DONE"
}

rh_listing() {
  find "$RH/root" ! -type d ! -name STOP -printf '%P\t%s\t%T@\n' | LC_ALL=C sort > "$1"
}

rh_launch() {
  local unit=$1
  shift
  bash "$TREE/scripts/pod_run.sh" --profile "$RH/profile.yaml" --image "$IMAGE" --digest "$DIGEST" --mem 32G \
    --threads "$RH_T" --bind /tmp/if2453 --bind "$ROOT:$ROOT:ro" --bind "$PFAMDIR:$PFAMDIR:ro" \
    --env "PYTHONPATH=$PYPATH" --log "$LOGS/$unit.log" "$unit" \
    build-stage --profile "$RH/profile.yaml" --root "$RH/root" "$@" >> "$LOG" 2>&1
  sleep 5
  while systemctl --user is-active --quiet "$unit"; do sleep 5; done
}

rehearse() {
  local unit rc t
  [ -e "$RH/OK" ] && { log "rehearsal: OK present"; return; }
  [ -e "$ROOT/STOP" ] && fail "$ROOT/STOP exists"
  mkdir -p "$RH/tmp" "$RH/durable" "$RH/bulk" "$RH/check" "$RH/run"
  if [ ! -e "$RH/root" ]; then
    python3 -I "$BIN/c26_view.py" "$ROOT" "$RH/root" "$EARLIER" pfam-38.2 >> "$LOG" 2>&1 || fail "rehearsal: view"
    # the parked stage as the real root holds it (its sequences and index, no gene calls): read-only links
    for t in l0seq/contig_index l0seq/seq_locator l0seq/seq_shards l0seq/genome_input runs; do
      [ -d "$ROOT/$t/stage=$PARKED" ] || fail "rehearsal: no $ROOT/$t/stage=$PARKED"
      ln -s "$ROOT/$t/stage=$PARKED" "$RH/root/$t/stage=$PARKED" || fail "rehearsal: link $t/stage=$PARKED"
    done
    log "rehearsal: view $RH/root made ($EARLIER with pfam-38.2, and $PARKED's l0seq tables and runs)"
  fi
  if [ ! -e "$RH/inputs/genomes.json" ]; then
    rm -rf "$RH/inputs"
    bb py "$BIN/c48c_v1_est.py" sample --inputs "$IN" --out "$RH/inputs" > "$RH/sample.json" 2>> "$LOG" \
      || fail "rehearsal: sample"
    log "rehearsal: sample: $(cat "$RH/sample.json")"
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
  bb py "$BIN/c48c_v1_est.py" before --root "$RH/root" --profile "$RH/profile.yaml" --inputs "$RH/inputs" \
    --stage "$RHSTAGE" --forbid "$PARKED" > "$RH/check/before.json" 2>> "$LOG"
  rc=$?
  # exit 6 (C-39's figures are for the whole stage, not a sample) is expected here; 5 (a contig C5 owns) is not
  [ "$rc" = 0 ] || [ "$rc" = 6 ] || fail "rehearsal: before est exit $rc"
  local args=(
    --stage "$RHSTAGE" --depth annotated
    --genomes "$RH/inputs/genomes.json" --objects "$RH/inputs/objects.json"
    --library pfam 38.2 "$PFAM" gathering --crispr
    --threads "$RH_T" --scan-cpus "$RH_CPUS"
  )
  unit=c48c-rh-1
  if [ ! -e "$RH/run/OK-$unit" ]; then
    rh_listing "$RH/check/before-listing.tsv"
    rh_launch "$unit" "${args[@]}"
    rc=$(exit_of "$unit")
    log "rehearsal: $unit exit '$rc'"
    [ "$rc" = 0 ] || fail "rehearsal: $unit exited '$rc', expected 0 (stage-only --depth annotated)"
    for t in $STEPS; do
      [ -e "$RH/root/runs/stage=$RHSTAGE/$t.json" ] || fail "rehearsal: missing run record $t"
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
  bb py "$BIN/c48c_v1_check.py" --root "$RH/root" --before "$RH/check/before" --after "$RH/check/after" \
    --log "$LOGS/c48c-rh-1.log" --objects "$RH/inputs/objects.json" --out "$RH/check/out" --tmp "$RH/tmp/check" \
    --stage "$RHSTAGE" --earlier "$HASHED" --before-est "$RH/check/before.json" --no-expect \
    > "$RH/check/check.out" 2>&1 || fail "rehearsal: c48c_v1_check.py failed (read $RH/check/check.out)"
  log "rehearsal: c48c_v1_check.py passed"
  unit=c48c-rh-2
  if [ ! -e "$RH/run/OK-$unit" ]; then
    rh_listing "$RH/check/again-before.tsv"
    rh_launch "$unit" "${args[@]}"
    rc=$(exit_of "$unit")
    [ "$rc" = 0 ] || fail "rehearsal: $unit exited '$rc'"
    grep -q ' published ' "$LOGS/$unit.log" && fail "rehearsal: relaunch published a step"
    rh_listing "$RH/check/again-after.tsv"
    cmp -s "$RH/check/again-before.tsv" "$RH/check/again-after.tsv" || fail "rehearsal: relaunch changed view"
    log "rehearsal: $unit skipped all steps and left the view unchanged"
    touch "$RH/run/OK-$unit"
  fi
  touch "$RH/OK"
  log "rehearsal: OK"
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
  echo "c48c-v1-$n"
}

stage_records() { find "$ROOT/runs/stage=$STAGE" -maxdepth 1 -name '*.json' ! -name '*.publishing.json' 2> /dev/null | wc -l; }

stop_unit() {
  log "stopping $1: $2"
  systemctl --user stop "$1"
  while systemctl --user is-active --quiet "$1"; do sleep 2; done
  fail "$2 (unit $1 stopped)"
}

rates_now() {
  bb py "$BIN/c48c_v1_est.py" rates --root "$ROOT" --objects "$IN/objects.json" --stage "$STAGE" \
    --before "$C/estimate/before.json" 2>> "$LOG"
}

run() {
  local unit free act rc line nrec=-1 n k=0
  [ -e "$C/before/DONE" ] || fail "run: no $C/before/DONE"
  [ -e "$RH/OK" ] || fail "run: no $RH/OK"
  [ -e "$C/estimate/DONE" ] || fail "run: no $C/estimate/DONE (the C5 gate and the estimates)"
  [ -e "$C/run/DONE" ] && { log "run: DONE present"; return; }
  [ -e "$C/run/FAILED" ] && fail "run: $C/run/FAILED exists (a stop stays a stop; read it)"
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
        --state "$C/run" --unit-prefix c48c-v1- --tmp-min-gb "$TMP_MIN_GB" --scratch-max-gb 100 \
        --cpu-max "$CPU_MAX" --stages "$STAGE")
      rc=$?
      echo "$line" >> "$C/run/guard.jsonl"
      [ "$rc" = 3 ] && stop_unit "$unit" "guard: $line"
      [ $((k % 60)) = 0 ] && log "guard: $(echo "$line" | cut -c1-400); log: $(grep -E 'run-step|build-stage' "$LOGS/$unit.log" | tail -n 1 | cut -c1-240)"
    fi
    n=$(stage_records)
    if [ "$n" != "$nrec" ]; then
      nrec=$n
      line=$(rates_now)
      rc=$?
      echo "$line" >> "$C/run/rates.jsonl"
      log "rates ($n records of $STAGE): $(echo "$line" | cut -c1-500)"
      [ "$rc" = 3 ] && stop_unit "$unit" "DEC-53 stop (a judged step more than 20% over its estimate): $(echo "$line" | cut -c1-900)"
    fi
    k=$((k + 1))
    sleep 10
  done
  line=$(rates_now)
  rc=$?
  echo "$line" >> "$C/run/rates.jsonl"
  log "rates at end: $(echo "$line" | cut -c1-600)"
  [ "$rc" = 3 ] && fail "DEC-53: a judged step more than 20% over its estimate at the end: $(echo "$line" | cut -c1-900)"
  rc=$(exit_of "$unit")
  log "$unit ended, exit '$rc'; usage: $(grep -E 'pod_run: .* usage' "$LOGS/$unit.log" | tail -n 1)"
  [ "$rc" = 0 ] || fail "$unit exited '$rc', expected 0"
  for t in $STEPS; do
    [ -e "$ROOT/runs/stage=$STAGE/$t.json" ] || fail "missing run record $STAGE/$t"
  done
  touch "$C/run/DONE"
}

after() {
  local d=$C/after r ok=1
  [ -e "$C/run/DONE" ] || fail "after: no $C/run/DONE"
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
  [ -e "$C/after/DONE" ] || fail "check: no $C/after/DONE"
  mkdir -p "$C/check"
  unit=$(tail -n 1 "$C/units")
  bb py "$BIN/c48c_v1_check.py" --root "$ROOT" --before "$C/before" --after "$C/after" --log "$LOGS/$unit.log" \
    --objects "$IN/objects.json" --out "$C/check" --tmp "$C/tmp/check" --stage "$STAGE" \
    --earlier "$HASHED" --before-est "$C/estimate/before.json" > "$C/check/check.out" 2>&1
  rc=$?
  cp "$C/check/root_added.tsv" "$C/root_added.tsv"
  log "check: exit $rc: $(tail -n 1 "$C/check/check.out" | cut -c1-600)"
  [ "$rc" = 0 ] || fail "check: c48c_v1_check.py exit $rc (read $C/check/check.out)"
  touch "$C/check/DONE"
}

backup() {
  local local_tar="$C/added_v1.tar" tmp_tar="$BACKUP_DIR/.added_v1.tar.tmp" sz sha1 sha2 n
  [ -e "$C/check/DONE" ] || fail "backup: no $C/check/DONE"
  [ -s "$C/root_added.tsv" ] || fail "backup: missing $C/root_added.tsv"
  [ ! -e "$BACKUP_TAR" ] || fail "backup: $BACKUP_TAR already exists"
  [ ! -e "$tmp_tar" ] || fail "backup: $tmp_tar already exists"
  awk -F'\t' 'NR > 1 {print $1}' "$C/root_added.tsv" > "$C/added_paths.txt"
  n=$(wc -l < "$C/added_paths.txt")
  tar -C "$ROOT" -cf "$local_tar" -T "$C/added_paths.txt" || fail "backup: tar failed"
  [ "$(tar -tf "$local_tar" | wc -l)" = "$n" ] || fail "backup: the tar does not hold the $n added files"
  sz=$(stat -c %s "$local_tar")
  [ "$sz" -ge 67108864 ] || fail "backup: $local_tar is $sz bytes (< 64 MiB)"
  sha1=$(sha256sum "$local_tar" | cut -d' ' -f1)
  log "backup: local tar $local_tar: $n files, $sz bytes, sha256 $sha1"
  dd if="$local_tar" of="$tmp_tar" bs=16M conv=fsync 2>> "$LOG" || fail "backup: dd failed"
  sha2=$(sha256sum "$tmp_tar" | cut -d' ' -f1)
  [ "$sha1" = "$sha2" ] || fail "backup: readback sha256 $sha2 != $sha1 (left $tmp_tar for the PI)"
  mv -n "$tmp_tar" "$BACKUP_TAR" || fail "backup: mv -n failed"
  [ -f "$BACKUP_TAR" ] || fail "backup: $BACKUP_TAR missing after mv"
  [ ! -e "$tmp_tar" ] || fail "backup: $tmp_tar still there after mv -n (did $BACKUP_TAR appear meanwhile?)"
  printf "%s\t%s\t%s\t%s\n" "$BACKUP_TAR" "$sz" "$sha2" "$n" > "$C/backup.tsv"
  log "backup: $BACKUP_TAR written ($n files, $sz bytes, sha256 $sha2)"
  rm -f "$local_tar"
}

MODE=${1:-}
log "c48c v1 driver: mode '$MODE' (pid $$), tree $TREE, image $IMAGE"
case "$MODE" in
  estimate) estimate ;;
  rehearse) rehearse ;;
  before) before ;;
  run) run ;;
  after) after ;;
  check) check ;;
  backup) backup ;;
  *) echo "usage: c48c_v1_driver.sh estimate|rehearse|before|run|after|check|backup" >&2; exit 2 ;;
esac
log "c48c v1 driver: mode '$MODE' done"
