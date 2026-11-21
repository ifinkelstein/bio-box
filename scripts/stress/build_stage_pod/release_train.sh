#!/bin/bash
# Card C-46: the release train (DEC-52, DEC-53, DEC-54, B-49, W-35).
#
# Cuts the next rolling release from every complete stage (`called` or `annotated`) in the data root,
# preserving the latest defined release's stage order and appending newly complete stages in `inputs.json`
# timestamp order. Incomplete and `stored`-only stages are named and left out, never waited for.
# Never deletes or moves anything of an earlier release.
#
# Usage (on the pod host):
#   bash scripts/stress/build_stage_pod/release_train.sh --dry-run [--root R] [--profile P] [--date YYYY-MM-DD] \
#        [--release R] [--parts K] [--same-stages]
#   bash scripts/stress/build_stage_pod/release_train.sh run [--root R] [--profile P] [--date YYYY-MM-DD] \
#        [--release R] [--parts K] [--same-stages]
#   bash scripts/stress/build_stage_pod/release_train.sh rehearse
#
# Environment overrides:
#   ROOT        data root (default: /tmp/if2453/r10/root)
#   WORK        driver working directory outside ROOT (default: /tmp/if2453/c46/train)
#   TREE        code tree checkout (default: repo root of this script)
#   PROFILE     host profile YAML (default: $TREE/profiles/wilkcomp02.yaml)
#   IMAGE       runtime .sif (default: /stor/home/if2453/bio-box/images/bio-box_faec1f4.sif)
#   DIGEST      expected sha256: digest of IMAGE
#   PFAMDIR     Pfam directory to bind read-only
#   PFAM        Pfam-A.hmm.gz path for pack build
#   T, J, CPUS, Q, PARTS, PACK_PARTS, PACK_JOBS, MEM, PACK_MEM
#   CPU_MAX (280), TMP_MIN_GB (150), SCRATCH_MAX_GB (250), START_MIN_GB (250)
set -u

BIN=$(cd "$(dirname "$0")" && pwd)
TREE=${TREE:-$(cd "$BIN/../../.." && pwd)}
ROOT=${ROOT:-/tmp/if2453/r10/root}
WORK=${WORK:-/tmp/if2453/c46/train}
REALPROF=${PROFILE:-$TREE/profiles/wilkcomp02.yaml}
IMAGE=${IMAGE:-/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif}
DIGEST=${DIGEST:-sha256:d67300aa1321b571fa8da79e74dfccdd87a8dca2e56518f324902efa6b9f098a}
TESTIMAGE=${TESTIMAGE:-/stor/home/if2453/bio-box/images/bio-box-test_faec1f4.sif}
PYPATH=$TREE/src
PFAMDIR=${PFAMDIR:-/stor/scratch/Wilke/if2453/all-the-data/l0/refdb/pfam/38.2}
PFAM=${PFAM:-$PFAMDIR/Pfam-A.hmm.gz}
REAL_ROOT=/tmp/if2453/r10/root

T=${T:-18}
J=${J:-2}
CPUS=${CPUS:-18}
Q=${Q:-40}
MEM=${MEM:-200G}
PARTS=${PARTS:-}
PACK_PARTS=${PACK_PARTS:-10}
PACK_JOBS=${PACK_JOBS:-12}
PACK_MEM=${PACK_MEM:-48G}
EX_THREADS=${EX_THREADS:-2}
CPU_MAX=${CPU_MAX:-280}
TMP_MIN_GB=${TMP_MIN_GB:-150}
SCRATCH_MAX_GB=${SCRATCH_MAX_GB:-250}
START_MIN_GB=${START_MIN_GB:-250}
ALLOW_EXAMPLES_FAIL=${ALLOW_EXAMPLES_FAIL:-0}
UNIT_PREFIX=${UNIT_PREFIX:-c46-train}

MODE=""
DRY_RUN=0
DATE_ARG=""
EXPLICIT_REL=""
SAME_STAGES=0
NO_LIVE_CPU=0

while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run|dry-run) DRY_RUN=1; MODE=dry-run; shift ;;
    run|rehearse) MODE=$1; shift ;;
    --root) ROOT=$2; shift 2 ;;
    --work) WORK=$2; shift 2 ;;
    --profile) REALPROF=$2; shift 2 ;;
    --date) DATE_ARG=$2; shift 2 ;;
    --release) EXPLICIT_REL=$2; shift 2 ;;
    --parts) PARTS=$2; shift 2 ;;
    --same-stages) SAME_STAGES=1; shift ;;
    --no-live-cpu) NO_LIVE_CPU=1; shift ;;
    -h|--help)
      sed -n '2,25p' "$0"
      exit 0
      ;;
    *) echo "release_train.sh: unknown argument '$1'" >&2; exit 2 ;;
  esac
done
[ -n "$MODE" ] || MODE=dry-run

LOGS=$WORK/logs
LOG=$WORK/driver.log
PROF_COPY=$WORK/profile.yaml
SCRATCH=$WORK/tmp/scratch

log() {
  mkdir -p "$WORK"
  echo "$(date '+%F %T %Z') $*" >> "$LOG"
}

fail() {
  log "FAILED: $*"
  echo "release_train: FAILED: $*" >&2
  exit 1
}

bb_ro() {
  # Short read-only command in the runtime image (capped scope, nice 19, OMP_NUM_THREADS=2).
  local what=$1 img=${IMG:-$IMAGE} ro_binds=()
  shift
  case "$what" in
    cli) set -- python -m biobox.cli "$@" ;;
    py) set -- python "$@" ;;
  esac
  [ -d "$ROOT" ] && ro_binds+=(--bind "$ROOT:$ROOT:ro")
  [ -d "$REAL_ROOT" ] && [ "$ROOT" != "$REAL_ROOT" ] && ro_binds+=(--bind "$REAL_ROOT:$REAL_ROOT:ro")
  [ -d "$PFAMDIR" ] && ro_binds+=(--bind "$PFAMDIR:$PFAMDIR:ro")
  systemd-run --user --scope --quiet -p MemoryMax=16G -p MemorySwapMax=0 -p CPUQuota=200% nice -n 19 \
    singularity exec --cleanenv --env PYTHONPATH="$PYPATH" --env OMP_NUM_THREADS=2 --env PYTHONUNBUFFERED=1 \
    --bind /tmp/if2453 "${ro_binds[@]}" --bind /sys/fs/cgroup:/sys/fs/cgroup:ro "$img" "$@"
}

make_profile() {
  mkdir -p "$WORK/tmp" "$LOGS"
  sed "s#^tmp:.*#tmp: $WORK/tmp#" "$REALPROF" > "$PROF_COPY"
}

build_units() {
  systemctl --user list-units --type=service --state=active --no-legend --plain \
    'r50-*' 'c24-h1-*' 'c25-*' 'c26-pl1-*' 'c33-*' 'c34-*' 'c35-rel-*' 'c35-pack*' 'b47-*' "${UNIT_PREFIX}-rel-*" \
    | awk '{print $1}' | paste -sd' ' -
}

free_gb() { df -B1G --output=avail /tmp | tail -n 1 | tr -d ' '; }
exit_of() { sed -n "s/^pod_run: $1 exited \([0-9]*\)$/\1/p" "$LOGS/$1.log" | tail -n 1; }
hash_of() { sed -n 's/.*manifest hash \([0-9a-f]*\).*/\1/p' "$1" | tail -n 1; }

listing() {
  find "$1" ! -type d ! -name STOP -printf '%P\t%s\t%T@\n' | LC_ALL=C sort > "$2"
}

verify_rel() {
  local root=$1 rel=$2 outdir=$3
  ROOT=$root bb_ro cli release verify --root "$root" --release "$rel" > "$outdir/verify-$rel.log" 2>&1
  local rc=$?
  log "verify $rel ($root): exit $rc: $(grep -E '^release ' "$outdir/verify-$rel.log" | tail -n 1)"
  return $rc
}

plan_args() {
  PLAN_ARGS=(--root "$ROOT" --profile "$REALPROF")
  [ -n "$DATE_ARG" ] && PLAN_ARGS+=(--date "$DATE_ARG")
  [ -n "$EXPLICIT_REL" ] && PLAN_ARGS+=(--release "$EXPLICIT_REL")
  [ -n "$PARTS" ] && PLAN_ARGS+=(--parts "$PARTS")
  [ "$SAME_STAGES" = 1 ] && PLAN_ARGS+=(--same-stages)
  [ "$NO_LIVE_CPU" = 1 ] && PLAN_ARGS+=(--no-live-cpu)
}

do_dry_run() {
  plan_args
  bb_ro py "$BIN/train_plan.py" "${PLAN_ARGS[@]}"
  return $?
}

stop_unit() {
  local unit=$1 reason=$2
  log "stopping $unit: $reason"
  systemctl --user stop "$unit"
  while systemctl --user is-active --quiet "$unit"; do sleep 2; done
  fail "$reason (unit $unit stopped; root left at last published step; rerun to resume)"
}

watch_unit() {
  local root=$1 rel=$2 unit=$3 st=$4 est=$5 prefix=$6 scr=$7 k=0 n nrec=-1 line rc foreign=0
  while systemctl --user is-active --quiet "$unit"; do
    if [ "$foreign" = 0 ] && [ -e "$root/STOP" ]; then
      foreign=1
      log "STOP appeared in $root: $(head -c 300 "$root/STOP"); waiting for $unit to end at step boundary"
    fi
    if [ $((k % 6)) = 0 ]; then
      line=$(python3 -I "$BIN/r50_watch.py" --root "$root" --scratch "$scr" --since "$st/since" --state "$st" \
        --unit-prefix "$prefix" --tmp-min-gb "$TMP_MIN_GB" --scratch-max-gb "$SCRATCH_MAX_GB" --cpu-max "$CPU_MAX" \
        --stages "$rel")
      rc=$?
      echo "$line" >> "$st/guard.jsonl"
      [ "$rc" = 3 ] && stop_unit "$unit" "guard: $line"
    fi
    n=$(find "$root/runs/stage=$rel" -maxdepth 1 -name '*.json' ! -name '*.publishing.json' 2> /dev/null | wc -l)
    if [ "$n" != "$nrec" ]; then
      nrec=$n
      line=$(python3 -I "$BIN/train_est.py" rates --root "$root" --release "$rel" --est "$est" 2>> "$LOG")
      rc=$?
      echo "$line" >> "$st/rates.jsonl"
      [ "$rc" = 3 ] && stop_unit "$unit" "rate guard >20% over estimate (DEC-43, DEC-53): $(echo "$line" | cut -c1-1000)"
    fi
    k=$((k + 1))
    sleep 5
  done
  line=$(python3 -I "$BIN/train_est.py" rates --root "$root" --release "$rel" --est "$est" --out "$st/rates-end.json" 2>> "$LOG")
  rc=$?
  log "rates at end of $unit (exit $rc): $(echo "$line" | cut -c1-900)"
  [ "$rc" = 3 ] && fail "rate guard >20% over estimate at end of $unit (DEC-43, DEC-53): $(echo "$line" | cut -c1-1000)"
  [ "$foreign" = 1 ] && return 3
  return 0
}

do_run() {
  mkdir -p "$WORK" "$LOGS" "$WORK/before" "$WORK/estimate" "$WORK/run" "$WORK/after" "$WORK/pack" "$WORK/check"
  # $WORK/DONE is the success marker of a run (card A-11): written last, only when every step passed. A waiter tests
  # it, never a log line or an exit status alone; a run withdraws the marker of an earlier one before it starts.
  rm -f "$WORK/DONE"
  make_profile
  log "release_train run: root=$ROOT work=$WORK profile=$REALPROF"

  [ -e "$ROOT/STOP" ] && fail "$ROOT/STOP exists: $(head -c 300 "$ROOT/STOP")"
  local act free
  act=$(build_units)
  [ -z "$act" ] || fail "active build unit(s): $act"
  free=$(free_gb)
  [ "$free" -ge "$START_MIN_GB" ] || fail "/tmp has ${free} GB free, under ${START_MIN_GB} GB"

  # 1. Plan (rules 1 & 2)
  plan_args
  local plan_json=$WORK/plan.json
  bb_ro py "$BIN/train_plan.py" "${PLAN_ARGS[@]}" --json > "$plan_json" || exit $?
  bb_ro py "$BIN/train_plan.py" "${PLAN_ARGS[@]}" > "$WORK/plan.txt" || { cat "$WORK/plan.txt"; fail "train_plan.py failed"; }
  cat "$WORK/plan.txt"

  local rel stages latest nothing_new earlier_rels last_stage last_depth last_genomes last_objects last_crispr
  rel=$(python3 -I -c 'import json,sys; print(json.load(open(sys.argv[1]))["release"])' "$plan_json")
  stages=$(python3 -I -c 'import json,sys; print(",".join(json.load(open(sys.argv[1]))["stages"]))' "$plan_json")
  latest=$(python3 -I -c 'import json,sys; print(json.load(open(sys.argv[1]))["latest_release"] or "")' "$plan_json")
  nothing_new=$(python3 -I -c 'import json,sys; print("1" if json.load(open(sys.argv[1]))["nothing_new"] else "0")' "$plan_json")
  earlier_rels=$(python3 -I -c 'import json,sys; print(" ".join(json.load(open(sys.argv[1]))["earlier_manifest_releases"]))' "$plan_json")

  if [ "$nothing_new" = 1 ] && [ "$SAME_STAGES" != 1 ]; then
    echo "release_train: 0 new complete stages since ${latest:-none}; nothing to cut (pass --same-stages to cut anyway)"
    log "noop: 0 new complete stages since ${latest:-none}"
    return 0
  fi
  [ -n "$stages" ] || fail "no complete stages found in $ROOT"

  # 2. Before: verify every earlier manifest release & snapshot listing
  local r ok=1
  for r in $earlier_rels; do
    verify_rel "$ROOT" "$r" "$WORK/before" || ok=0
  done
  [ "$ok" = 1 ] || fail "an earlier release did not verify before launch (see $WORK/before/verify-*.log)"
  listing "$ROOT" "$WORK/before/listing.tsv"
  log "before: $(wc -l < "$WORK/before/listing.tsv") files in $ROOT"

  # 3. Estimate before launch
  local est_args=(before --root "$ROOT" --release "$rel" --stages "$stages" --profile "$PROF_COPY" --tmp "$WORK/tmp/est")
  [ -n "$latest" ] && est_args+=(--earlier "$latest")
  [ -n "$PARTS" ] && est_args+=(--parts "$PARTS")
  bb_ro py "$BIN/train_est.py" "${est_args[@]}" > "$WORK/estimate/before.json" 2> "$WORK/estimate/before.err" \
    || fail "train_est.py before failed (see $WORK/estimate/before.err)"

  # 4. Build release steps via build-stage (stage steps for last_stage all skip)
  local bs_args=()
  while IFS= read -r line; do
    bs_args+=("$line")
  done < <(python3 -I - "$plan_json" <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))
ls = d["last_stage"]
args = [
    "--stage", ls["stage"],
    "--release", d["release"],
    "--stages", ",".join(d["stages"]),
    "--depth", ls["depth"],
    "--genomes", ls["genomes"],
    "--objects", ls["objects"],
]
if ls["depth"] == "annotated":
    for lib in ls.get("library_specs", []):
        args += ["--library", lib["name"], lib["version"], lib["hmm"], lib.get("cutoff") or "gathering"]
    if ls.get("ncrna_cm"):
        args += ["--ncrna-cm", ls["ncrna_cm"], "--ncrna-set", ls.get("ncrna_set") or "prok-2"]
if ls.get("crispr"):
    args.append("--crispr")
for a in args:
    print(a)
PY
)
  bs_args+=(--threads "$T" --scan-jobs "$J" --scan-cpus "$CPUS")
  [ -n "$PARTS" ] && bs_args+=(--cluster-parts "$PARTS")

  local unit="${UNIT_PREFIX}-rel-1"
  mkdir -p "$WORK/run/$unit"
  touch "$WORK/run/$unit/since"
  local ro_binds=()
  [ -d "$REAL_ROOT" ] && [ "$ROOT" != "$REAL_ROOT" ] && ro_binds+=(--bind "$REAL_ROOT:$REAL_ROOT:ro")
  [ -d "$PFAMDIR" ] && ro_binds+=(--bind "$PFAMDIR:$PFAMDIR:ro")
  log "launch $unit: build-stage --profile $PROF_COPY --root $ROOT ${bs_args[*]}"
  bash "$TREE/scripts/pod_run.sh" --profile "$PROF_COPY" --image "$IMAGE" --digest "$DIGEST" --mem "$MEM" \
    --threads "$Q" --bind /tmp/if2453 "${ro_binds[@]}" --env "PYTHONPATH=$PYPATH" \
    --log "$LOGS/$unit.log" "$unit" build-stage --profile "$PROF_COPY" --root "$ROOT" "${bs_args[@]}" >> "$LOG" 2>&1
  sleep 2
  watch_unit "$ROOT" "$rel" "$unit" "$WORK/run/$unit" "$WORK/estimate/before.json" "${UNIT_PREFIX}-rel-" "$SCRATCH"
  local wrc=$?
  [ "$wrc" = 3 ] && fail "STOP appeared in $ROOT during $unit; left in place"
  local rc
  rc=$(exit_of "$unit")
  [ "$rc" = 0 ] || fail "$unit exited '$rc' (see $LOGS/$unit.log)"
  grep -qE "build-stage: stage [^ ]+ and release $rel: all [0-9]+ steps done" "$LOGS/$unit.log" \
    || fail "$unit did not finish all steps (see $LOGS/$unit.log)"
  find "$ROOT/runs" -name '*.publishing.json' | grep -q . && fail "unfinished .publishing.json left in $ROOT/runs"

  # 5. After: verify new release and every earlier release (hashes unchanged)
  listing "$ROOT" "$WORK/after/listing.tsv"
  ok=1
  verify_rel "$ROOT" "$rel" "$WORK/after" || ok=0
  local new_hash h0 h1
  new_hash=$(hash_of "$WORK/after/verify-$rel.log")
  for r in $earlier_rels; do
    verify_rel "$ROOT" "$r" "$WORK/after" || ok=0
    h0=$(hash_of "$WORK/before/verify-$r.log")
    h1=$(hash_of "$WORK/after/verify-$r.log")
    [ -n "$h1" ] && [ "$h0" = "$h1" ] || { log "after: $r hash '$h1' != '$h0'"; ok=0; }
  done
  [ "$ok" = 1 ] || fail "post-build release verify failed or earlier hash changed"
  bb_ro cli releases --root "$ROOT" --profile "$PROF_COPY" > "$WORK/after/releases.txt" 2>&1 \
    || fail "bio-box releases failed (see $WORK/after/releases.txt)"

  # 6. Search pack in parts (W-35) + pack verify
  local punit="${UNIT_PREFIX}-pack" pack_model=()
  [ -f "$PFAM" ] && pack_model=(--model "pfam-38.2=$PFAM")
  bash "$TREE/scripts/pod_run.sh" --profile "$PROF_COPY" --image "$IMAGE" --digest "$DIGEST" --mem "$PACK_MEM" \
    --threads "$PACK_JOBS" --bind /tmp/if2453 --bind "$ROOT:$ROOT:ro" "${ro_binds[@]}" --env "PYTHONPATH=$PYPATH" \
    --log "$LOGS/$punit.log" "$punit" pack build --root "$ROOT" --release "$rel" "${pack_model[@]}" \
    --profile "$PROF_COPY" --parts "$PACK_PARTS" --jobs "$PACK_JOBS" >> "$LOG" 2>&1
  sleep 2
  while systemctl --user is-active --quiet "$punit"; do sleep 5; done
  rc=$(exit_of "$punit")
  [ "$rc" = 0 ] || fail "pack build exited '$rc' (see $LOGS/$punit.log)"
  bb_ro cli pack verify --root "$ROOT" --release "$rel" --profile "$PROF_COPY" > "$WORK/pack/verify.log" 2>&1 \
    || fail "pack verify failed (see $WORK/pack/verify.log)"
  local pack_id
  pack_id=$(sed -n 's/.*pack \([0-9a-f]\{32\}\).*/\1/p' "$LOGS/$punit.log" | tail -n 1)

  # 7. Worked examples release check
  rm -rf "$WORK/examples"
  TREE=$TREE PROFILE=$PROF_COPY IMAGE=$IMAGE DIGEST=$DIGEST ROOT=$ROOT OUT=$WORK/examples \
    UNIT="${UNIT_PREFIX}-ex" THREADS=$EX_THREADS PFAM=$PFAM \
    bash "$TREE/scripts/release_examples.sh" "$rel" > "$WORK/examples.log" 2>&1
  rc=$?
  log "examples: release_examples.sh $rel exit $rc: $(tail -n 1 "$WORK/examples.log" | cut -c1-400)"
  if [ "$rc" != 0 ] && [ "$ALLOW_EXAMPLES_FAIL" != 1 ]; then
    fail "release_examples.sh failed with exit $rc (see $WORK/examples.log)"
  fi

  # 8. Post-release check, root_added.tsv, and announcement draft
  local chk_args=(--root "$ROOT" --before "$WORK/before/listing.tsv" --after "$WORK/after/listing.tsv"
                  --release "$rel" --stages "$stages" --out "$WORK/check")
  [ -n "$latest" ] && chk_args+=(--earlier "$latest")
  [ -n "$new_hash" ] && chk_args+=(--manifest-hash "$new_hash")
  [ -n "$pack_id" ] && chk_args+=(--pack-id "$pack_id")
  [ -s "$WORK/examples/summary.tsv" ] && chk_args+=(--examples-summary "$WORK/examples/summary.tsv")
  python3 -I "$BIN/train_check.py" "${chk_args[@]}" > "$WORK/check/check.out"
  rc=$?
  cat "$WORK/check/check.out"
  [ "$rc" = 0 ] || fail "train_check.py failed (see $WORK/check/check.out)"
  cp "$WORK/check/root_added.tsv" "$WORK/root_added.tsv" || fail "cp root_added.tsv"
  cp "$WORK/check/announcement.md" "$WORK/announcement.md" || fail "cp announcement.md"
  log "release_train run complete: $rel (manifest $new_hash, pack $pack_id)"
  echo "$rel" > "$WORK/DONE"
}

do_rehearse() {
  local rh=$WORK/rh
  mkdir -p "$rh/tmp" "$rh/durable" "$rh/bulk" "$rh/inputs" "$rh/logs"
  cat > "$rh/profile.yaml" << EOF
max_threads_day: 40
max_threads_night: 56
mem_scope: 16G
tmp: $rh/tmp
durable: $rh/durable
bulk: $rh/bulk
nice_monolithic: 19
nice_fanout: 15
day_start_hour: 8
day_end_hour: 20
cpu_min_per_process: 300
cluster_cpu_min_max: 20000
small_file_floor_bytes: 1000
require_mem_cap: true
EOF
  # Prepare 3 distinct stages (S0=lambda, S1=M. genitalium, S2=T4) + incomplete S3 (stored only)
  if [ ! -s "$rh/inputs/objects-S2.json" ]; then
    IMG=$TESTIMAGE ROOT=$REAL_ROOT bb_ro py "$BIN/c24_prep.py" "$TREE" "$rh" >> "$LOG" 2>&1 || fail "rehearse: c24_prep.py"
    python3 -I - "$rh/inputs" <<'PY'
import json, sys
from pathlib import Path
d = Path(sys.argv[1])
g0 = json.loads((d / "genomes-S0.json").read_text())
o0 = json.loads((d / "objects-S0.json").read_text())
g1 = json.loads((d / "genomes-S1.json").read_text())
o1 = json.loads((d / "objects-S1.json").read_text())
# Split S0 (2 genomes) into S0 (first genome) and S1 (second genome); move original S1 (T4) to S2; S3 mirrors S0
(d / "genomes-S0a.json").write_text(json.dumps([g0[0]], indent=1))
(d / "objects-S0a.json").write_text(json.dumps([o0[0]], indent=1))
(d / "genomes-S1a.json").write_text(json.dumps([g0[1]], indent=1))
(d / "objects-S1a.json").write_text(json.dumps([o0[1]], indent=1))
(d / "genomes-S2.json").write_text(json.dumps(g1, indent=1))
(d / "objects-S2.json").write_text(json.dumps(o1, indent=1))
(d / "genomes-S3.json").write_text(json.dumps([g0[0]], indent=1))
(d / "objects-S3.json").write_text(json.dumps([o0[0]], indent=1))
PY
  fi

  local r_root=$rh/root
  local r_logs=$rh/logs
  mkdir -p "$r_root" "$r_logs"

  rh_launch() {
    local u=$1
    shift
    bash "$TREE/scripts/pod_run.sh" --profile "$rh/profile.yaml" --image "$IMAGE" --digest "$DIGEST" --mem 8G \
      --threads 4 --bind /tmp/if2453 --bind "$REAL_ROOT:$REAL_ROOT:ro" --bind "$PFAMDIR:$PFAMDIR:ro" \
      --env "PYTHONPATH=$PYPATH" --log "$r_logs/$u.log" "$u" "$@" >> "$LOG" 2>&1
    sleep 2
    while systemctl --user is-active --quiet "$u"; do sleep 2; done
    local rc
    rc=$(sed -n "s/^pod_run: $u exited \([0-9]*\)$/\1/p" "$r_logs/$u.log" | tail -n 1)
    [ "$rc" = 0 ] || fail "rehearse: unit $u exited '$rc' (see $r_logs/$u.log)"
  }

  # 1. Build S0 (stage only), then S1 + initial release r26.11.1 over S0,S1 (called depth, fast, <=4 threads)
  if [ ! -e "$r_root/runs/stage=S0/tapes.json" ]; then
    rh_launch c46-rh-s0 build-stage --profile "$rh/profile.yaml" --root "$r_root" --stage S0 --depth called \
      --genomes "$rh/inputs/genomes-S0a.json" --objects "$rh/inputs/objects-S0a.json" --threads 2
  fi
  if [ ! -e "$r_root/runs/stage=r26.11.1/catalog.json" ]; then
    rh_launch c46-rh-r1 build-stage --profile "$rh/profile.yaml" --root "$r_root" --stage S1 --depth called \
      --release r26.11.1 --stages S0,S1 --genomes "$rh/inputs/genomes-S1a.json" \
      --objects "$rh/inputs/objects-S1a.json" --threads 2 --scan-jobs 1 --scan-cpus 2 --cluster-parts 2
  fi
  # 2. Build S2 (complete, called depth) and S3 (stored only -> must be named and left out!)
  if [ ! -e "$r_root/runs/stage=S2/tapes.json" ]; then
    rh_launch c46-rh-s2 build-stage --profile "$rh/profile.yaml" --root "$r_root" --stage S2 --depth called \
      --genomes "$rh/inputs/genomes-S2.json" --objects "$rh/inputs/objects-S2.json" --threads 2
  fi
  if [ ! -e "$r_root/runs/stage=S3/ingest.json" ]; then
    rh_launch c46-rh-s3 build-stage --profile "$rh/profile.yaml" --root "$r_root" --stage S3 --depth stored \
      --genomes "$rh/inputs/genomes-S3.json" --objects "$rh/inputs/objects-S3.json" --threads 2
  fi

  # 3. Dry-run before cut (must name r26.11.2, stages S0,S1,S2, exclude S3, write nothing)
  ROOT=$r_root WORK=$rh/train PROFILE=$rh/profile.yaml TREE=$TREE IMAGE=$IMAGE DIGEST=$DIGEST \
    bash "$0" --dry-run --date 2026-11-19 --parts 2 --no-live-cpu > "$rh/dry-run-before.txt" || fail "rehearse: dry-run before"
  grep -q "Next release: r26.11.2" "$rh/dry-run-before.txt" || fail "rehearse: dry-run did not name r26.11.2"
  grep -q "Stages (3): S0,S1,S2" "$rh/dry-run-before.txt" || fail "rehearse: dry-run did not pick S0,S1,S2"
  grep -q "S3 (stored)" "$rh/dry-run-before.txt" || fail "rehearse: dry-run did not exclude S3"

  # 4. Full release_train.sh run cutting r26.11.2 (lightweight: <=4 threads, 2 cluster parts, 2 pack parts)
  ROOT=$r_root WORK=$rh/train PROFILE=$rh/profile.yaml TREE=$TREE IMAGE=$IMAGE DIGEST=$DIGEST \
    T=2 J=1 CPUS=2 Q=4 MEM=8G PACK_MEM=8G PARTS=2 PACK_PARTS=2 PACK_JOBS=2 EX_THREADS=2 \
    TMP_MIN_GB=10 START_MIN_GB=10 ALLOW_EXAMPLES_FAIL=1 UNIT_PREFIX=c46-rh \
    bash "$0" run --date 2026-11-19 --parts 2 --no-live-cpu || fail "rehearse: release_train run"

  # 5. Dry-run after cut: 0 new stages -> noop, and --same-stages shows B-49 reuse of r26.11.2!
  ROOT=$r_root WORK=$rh/train PROFILE=$rh/profile.yaml TREE=$TREE IMAGE=$IMAGE DIGEST=$DIGEST \
    bash "$0" --dry-run --date 2026-11-19 --parts 2 --no-live-cpu > "$rh/dry-run-after.txt" || fail "rehearse: dry-run after"
  grep -q "Next release: r26.11.3" "$rh/dry-run-after.txt" || fail "rehearse: dry-run after did not name r26.11.3"
  grep -q "noop" "$rh/dry-run-after.txt" || fail "rehearse: dry-run after did not report noop"

  ROOT=$r_root WORK=$rh/train PROFILE=$rh/profile.yaml TREE=$TREE IMAGE=$IMAGE DIGEST=$DIGEST \
    bash "$0" --dry-run --date 2026-11-19 --parts 2 --same-stages --no-live-cpu > "$rh/dry-run-reuse.txt" || fail "rehearse: dry-run reuse"
  grep -q "REUSE from r26.11.2" "$rh/dry-run-reuse.txt" || fail "rehearse: dry-run --same-stages did not detect B-49 reuse"

  touch "$rh/OK"
  echo "rehearsal: OK (r26.11.2 cut over S0,S1,S2; S3 excluded; r26.11.1 unchanged; pack verified; B-49 reuse verified)"
}

case "$MODE" in
  dry-run) do_dry_run ;;
  run) do_run ;;
  rehearse) do_rehearse ;;
  *) echo "usage: release_train.sh [--dry-run | run | rehearse] [options]" >&2; exit 2 ;;
esac
