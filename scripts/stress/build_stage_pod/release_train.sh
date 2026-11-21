#!/bin/bash
# Card C-46: the release train (DEC-52, DEC-53, DEC-54, B-49, W-35).
#
# Cuts the next rolling release from every complete stage (`called` or `annotated`) in the data root,
# preserving the latest defined release's stage order and appending newly complete stages in `inputs.json`
# timestamp order. Incomplete and `stored`-only stages are named and left out, never waited for.
# Never deletes or moves anything of an earlier release.
#
# Card A-15: `run` in a $WORK whose release it launched and left unfinished (its plan.json names a release with run
# records in the root and no manifest, and the plan names it again with the same stages) resumes it: it keeps the first
# launch's before/listing.tsv, before/verify-*.log and estimate/before.json, verifies the earlier releases again into
# $WORK/resume-<n>/ (each manifest hash must equal the first launch's), launches the same build-stage command (checked
# against the first launch's line in driver.log, with the same profile copy) as unit <prefix>-rel-<n>, a name not used
# before in this $WORK, and judges the rate guard against the first launch's estimate. Another release or other
# stages while such a release is unfinished are refused. `--dry-run` says `Action: resume ...` in that state.
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

stage_confs_present() {  # card A-12: a tree without stages/*.conf would leave the stage driver's units out of the lock
  compgen -G "$BIN/stages/*.conf" > /dev/null
}

build_units() {
  # the older names, and the build units <stage>-<n> of every stage of card A-12's stage driver (stages/<STAGE>.conf)
  local pats=('r50-*' 'c24-h1-*' 'c25-*' 'c26-pl1-*' 'c33-*' 'c34-*' 'c35-rel-*' 'c35-pack*' 'b47-*' 'c47-c5-*'
              'c48c-*' "${UNIT_PREFIX}-rel-*") f s
  for f in "$BIN"/stages/*.conf; do
    [ -e "$f" ] || continue
    s=$(basename "$f" .conf)
    [[ "$s" =~ ^[A-Z][A-Za-z0-9]{0,15}$ ]] && pats+=("${s,,}-[0-9]*")
  done
  systemctl --user list-units --type=service --state=active,activating,deactivating --no-legend --plain "${pats[@]}" \
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

# Card A-15: a rerun of `run` resumes the release that this $WORK launched and left unfinished. `train_plan.py` plans a
# release with no `release_def` and no manifest again under the same name and stages, and `build-stage` skips the
# steps already published; what a resume must not do is take the root's listing, the earlier releases' hashes and the
# estimate again on a root that holds half the release, or launch a unit under a name whose log already ends with an
# exit line (pod_run.sh appends to the log: card A-14's watch would read the old unit's exit line as the new one's).
RS_ACTION="" RS_REL="" RS_RECORDS="" RS_UNIT="" RS_HOW="" RS_REASON="" RS_PROBLEMS="" RS_VERIFY_N=0 RS_CMD=""
BS_ARGS=()

resume_state() {  # NEWPLAN_JSON: first launch, resume or refusal, from this $WORK, the root and the new plan (RS_*).
  # Reads only. The unit number is one more than any `<prefix>-rel-<n>` recorded under $WORK/run/ or $WORK/logs/.
  local out k v
  out=$(python3 -I - "$WORK" "$ROOT" "$UNIT_PREFIX" "$1" <<'PY'
import json
import re
import sys
from pathlib import Path

work, root, prefix = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
new = json.loads(sys.argv[4])


def say(**kv):
    for k, v in kv.items():
        print(f"{k}\t{' '.join(str(v).split())}")


def run_dir(rel):
    return root / "runs" / f"stage={rel}"


def records(rel):  # run records published
    d = run_dir(rel)
    return [p for p in d.glob("*.json") if not p.name.endswith(".publishing.json")] if d.is_dir() else []


def started(rel):  # a run record or a run-step intent: the root holds something of the release
    d = run_dir(rel)
    return d.is_dir() and any(d.glob("*.json"))


def finished(rel):
    return (root / "manifests" / f"release={rel}.json").exists()


def defined(rel):
    return (root / "release_def" / f"release={rel}").is_dir()


def plan_of(p):
    return f"release {p['release']} (stages {','.join(p['stages']) or 'none'})"


def earlier(p):
    return p.get("earlier_manifest_releases") or []


name = re.compile(re.escape(prefix) + r"-rel-([1-9][0-9]*)")
seen = [0]
for d, suffix in ((work / "run", ""), (work / "logs", ".log")):
    for p in (d.iterdir() if d.is_dir() else ()):
        if suffix and not p.name.endswith(suffix):
            continue
        m = name.fullmatch(p.name[:len(p.name) - len(suffix)])
        if m:
            seen.append(int(m[1]))
unit = max(seen) + 1

try:
    old = json.loads((work / "plan.json").read_text())
    old = old if isinstance(old.get("release"), str) and isinstance(old.get("stages"), list) else None
except (OSError, ValueError, AttributeError):
    old = None

if old is not None and started(old["release"]) and not finished(old["release"]):
    rel, k = old["release"], len(records(old["release"]))
    if new["release"] == rel and new["stages"] == old["stages"] and earlier(new) == earlier(old):
        how = "the same plan"
    elif (defined(rel) and new.get("latest_release") == rel and new["stages"] == old["stages"]
          and not new.get("new_stages") and earlier(new) == earlier(old)):
        how = "the same stages; its release_def was written since the first launch"
    else:
        say(action="refuse", unit=unit, reason=(
            f"this $WORK ({work}) launched {plan_of(old)}, and the root holds {k} run records of it and no manifest: "
            f"it is unfinished; the plan now names {plan_of(new)}. Resume {rel} from a tree and root that plan it "
            f"again, or set it aside by the PI's decision"))
        sys.exit(0)
    say(action="resume", release=rel, records=k, unit=unit, how=how)
else:
    foreign = sorted(d.name.removeprefix("stage=") for d in (root / "runs").glob("stage=*")
                     if d.is_dir() and not (d / "inputs.json").exists() and any(d.glob("*.json"))
                     and not finished(d.name.removeprefix("stage=")) and not defined(d.name.removeprefix("stage=")))
    if foreign:
        say(action="refuse", unit=unit, reason=(
            f"the root holds run records of {', '.join(foreign)} (no release_def, no manifest: unfinished), which this "
            f"$WORK ({work}) did not launch (its plan.json names {old['release'] if old else 'no release'}); the plan "
            f"names {plan_of(new)}. Resume it with the $WORK that launched it"))
        sys.exit(0)
    say(action="first", release=new["release"], records=len(records(new["release"])), unit=unit)
PY
) || return 1
  RS_ACTION="" RS_REL="" RS_RECORDS="" RS_UNIT="" RS_HOW="" RS_REASON=""
  while IFS=$'\t' read -r k v; do
    case "$k" in
      action) RS_ACTION=$v ;; release) RS_REL=$v ;; records) RS_RECORDS=$v ;; unit) RS_UNIT=$v ;;
      how) RS_HOW=$v ;; reason) RS_REASON=$v ;;
    esac
  done <<< "$out"
  [ -n "$RS_ACTION" ] && [ -n "$RS_UNIT" ]
}

make_bs_args() {  # PLAN_JSON: build-stage's arguments for the plan's release into BS_ARGS (the first launch's, on a resume)
  local line
  BS_ARGS=()
  while IFS= read -r line; do
    BS_ARGS+=("$line")
  done < <(python3 -I - "$1" <<'PY'
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
  BS_ARGS+=(--threads "$T" --scan-jobs "$J" --scan-cpus "$CPUS")
  [ -n "$PARTS" ] && BS_ARGS+=(--cluster-parts "$PARTS")
  return 0
}

first_launch_cmd() {  # REL: the build-stage command of the last launch of REL that $WORK/driver.log records
  [ -f "$LOG" ] || return 0
  grep -F -- " --release $1 " "$LOG" | sed -n "s/^.* launch ${UNIT_PREFIX}-rel-[0-9]*: \(build-stage .*\)$/\1/p" \
    | tail -n 1
}

resume_checks() {  # REL: what a resume keeps from the first launch, and that it launches the same command; the
  # problems into RS_PROBLEMS (empty when it may go on). Reads only; BS_ARGS made from $WORK/plan.json first.
  local rel=$1 r n=0 first now
  RS_PROBLEMS=""
  [ -s "$WORK/before/listing.tsv" ] || RS_PROBLEMS+="; the first launch's $WORK/before/listing.tsv is missing or empty"
  [ -s "$WORK/estimate/before.json" ] || RS_PROBLEMS+="; the first launch's $WORK/estimate/before.json is missing or empty"
  for r in $(python3 -I -c 'import json,sys; print(" ".join(json.load(open(sys.argv[1])).get("earlier_manifest_releases") or []))' "$WORK/plan.json"); do
    n=$((n + 1))
    [ -n "$(hash_of "$WORK/before/verify-$r.log" 2> /dev/null)" ] \
      || RS_PROBLEMS+="; the first launch's $WORK/before/verify-$r.log is missing or holds no manifest hash"
  done
  if [ ! -s "$PROF_COPY" ]; then
    RS_PROBLEMS+="; the first launch's profile copy $PROF_COPY is missing"
  elif ! sed "s#^tmp:.*#tmp: $WORK/tmp#" "$REALPROF" | cmp -s - "$PROF_COPY"; then
    RS_PROBLEMS+="; the profile copy made from $REALPROF would differ from the first launch's $PROF_COPY (its cut and scratch must not change)"
  fi
  first=$(first_launch_cmd "$rel")
  now="build-stage --profile $PROF_COPY --root $ROOT ${BS_ARGS[*]}"
  if [ -z "$first" ]; then
    RS_PROBLEMS+="; $LOG records no launch of $rel to compare the command with"
  elif [ "$first" != "$now" ]; then
    RS_PROBLEMS+="; the command would differ from the first launch's: now '$now', first '$first'"
  fi
  RS_PROBLEMS=${RS_PROBLEMS#; }
  RS_VERIFY_N=$n
  RS_CMD=$now
  [ -z "$RS_PROBLEMS" ]
}

action_line() {  # LINE: the plan's text (stdin) with its `Action:` line replaced by LINE
  local l
  while IFS= read -r l; do
    case "$l" in "Action: "*) printf '%s\n' "$1" ;; *) printf '%s\n' "$l" ;; esac
  done
}

do_dry_run() {
  plan_args
  local newplan txt rc
  newplan=$(bb_ro py "$BIN/train_plan.py" "${PLAN_ARGS[@]}" --json) || return $?
  txt=$(bb_ro py "$BIN/train_plan.py" "${PLAN_ARGS[@]}")
  rc=$?
  if ! resume_state "$newplan"; then
    printf '%s\n' "$txt"
    echo "release_train: could not tell a resume from a first launch (the resume helper failed)" >&2
    return 1
  fi
  case "$RS_ACTION" in
    first) printf '%s\n' "$txt"; return $rc ;;
    refuse) action_line "Action: refused: $RS_REASON" <<< "$txt"; return 1 ;;
  esac
  action_line "Action: resume $RS_REL ($RS_RECORDS run records published; unit $RS_UNIT)" <<< "$txt"
  make_bs_args "$WORK/plan.json"
  resume_checks "$RS_REL"
  local ok=$?
  echo "Resume of $RS_REL ($RS_HOW), launched by this \$WORK ($WORK): unit ${UNIT_PREFIX}-rel-$RS_UNIT, its own log $LOGS/${UNIT_PREFIX}-rel-$RS_UNIT.log and state $WORK/run/${UNIT_PREFIX}-rel-$RS_UNIT"
  echo "  kept from the first launch: before/listing.tsv ($(wc -l 2> /dev/null < "$WORK/before/listing.tsv" || echo 0) files), $RS_VERIFY_N verify logs, estimate/before.json (the rate guard's estimate)"
  echo "  the $RS_VERIFY_N earlier releases are verified again into $WORK/resume-$RS_UNIT; each manifest hash must equal the first launch's"
  echo "  command: $RS_CMD"
  if [ "$ok" = 0 ]; then
    echo "  the same as the first launch's in $LOG"
  else
    echo "  refused: $RS_PROBLEMS"
  fi
  [ "$rc" = 0 ] && [ "$ok" = 0 ]
}

# The watch follows card A-12's stage driver (stage_driver.sh unit_ended, watch_answer, stop_unit, bgcall, on_exit).
RUNNING_UNIT=""
RUNNING_LOGS=""
BGPID=""

unit_ended() {  # UNIT LOGS: 0 on the unit's own exit line or a real answer of systemd that it is inactive or failed;
  # 1 on a real answer that it runs; 2 on no answer (the user bus or systemctl failed): it counts as running
  local s
  [ -n "$(sed -n "s/^pod_run: $1 exited \([0-9]*\)$/\1/p" "$2/$1.log" 2> /dev/null | tail -n 1)" ] && return 0
  if ! s=$(systemctl --user show -p ActiveState --value "$1" 2> /dev/null) || [ -z "$s" ]; then
    log "no answer from systemd about $1; it counts as running"
    return 2
  fi
  case "$s" in inactive|failed) return 0 ;; *) return 1 ;; esac
}

# How long the train goes on without one readable answer about its unit before it gives up (the user bus may be gone):
# it then says the unit's state is unknown and that it may still be running, and exits non-zero. Never claims it
# stopped. Same bound as stage_driver.sh.
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
    fail "no readable answer from systemd about unit $2 for ${WATCH_LIMIT_S} s: its state is unknown and it may still be running (systemctl --user list-units '${UNIT_PREFIX}-rel-*'; a STOP file in the root ends it at a step boundary)"
  fi
  return 1
}

stop_unit() {  # UNIT WHY: stop it and see it gone before saying so; fail loudly at its bound, never claiming it
  local unit=$1 reason=$2 k=0 t0 rc
  t0=$(date +%s)
  log "stopping $unit: $reason"
  while :; do
    unit_ended "$unit" "$LOGS"
    rc=$?
    [ "$rc" = 0 ] && break
    watch_answer "$rc" "$unit"
    if [ $(( $(date +%s) - t0 )) -ge "$STOP_LIMIT_S" ]; then
      RUNNING_UNIT=""
      fail "$reason; then unit $unit was asked to stop for ${STOP_LIMIT_S} s and is still not seen gone: it may still be running"
    fi
    if [ $((k % 15)) = 0 ]; then
      systemctl --user stop "$unit" || log "systemctl --user stop $unit failed; asking again"
    fi
    k=$((k + 1))
    sleep 2
  done
  RUNNING_UNIT=""
  fail "$reason (unit $unit stopped; root left at last published step; rerun to resume)"
}

bgcall() {  # OUT CMD...: CMD in the background, its output into OUT, waited for here, so that a TERM to the train runs
  # its trap at once (bash runs a trap only after a foreground command returns)
  local out=$1
  shift
  "$@" > "$out" 2>> "$LOG" &
  BGPID=$!
  wait "$BGPID"
}

on_exit() {  # the train is leaving while its unit runs (a signal, a crash): stop the unit, as far as it can
  local rc=$? k=0
  [ -n "$BGPID" ] && kill "$BGPID" 2> /dev/null
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

watch_unit() {
  local root=$1 rel=$2 unit=$3 st=$4 est=$5 prefix=$6 scr=$7 k=0 n nrec=-1 line rc foreign=0
  mkdir -p "$st"
  while :; do
    unit_ended "$unit" "$LOGS"
    rc=$?
    [ "$rc" = 0 ] && break
    watch_answer "$rc" "$unit"
    if [ "$foreign" = 0 ] && [ -e "$root/STOP" ]; then
      foreign=1
      log "STOP appeared in $root: $(head -c 300 "$root/STOP"); waiting for $unit to end at step boundary"
    fi
    if [ $((k % 6)) = 0 ]; then
      bgcall "$st/guard.out" python3 -I "$BIN/r50_watch.py" --root "$root" --scratch "$scr" --since "$st/since" --state "$st" \
        --unit-prefix "$prefix" --tmp-min-gb "$TMP_MIN_GB" --scratch-max-gb "$SCRATCH_MAX_GB" --cpu-max "$CPU_MAX" \
        --stages "$rel"
      rc=$?
      line=$(cat "$st/guard.out" 2> /dev/null)
      echo "$line" >> "$st/guard.jsonl"
      [ "$rc" = 0 ] || stop_unit "$unit" "guard (r50_watch.py exit $rc; 3: a stop of the process guard): $(echo "$line" | cut -c1-1000)"
    fi
    n=$(find "$root/runs/stage=$rel" -maxdepth 1 -name '*.json' ! -name '*.publishing.json' 2> /dev/null | wc -l)
    if [ "$n" != "$nrec" ]; then
      nrec=$n
      bgcall "$st/rates.out" python3 -I "$BIN/train_est.py" rates --root "$root" --release "$rel" --est "$est"
      rc=$?
      line=$(cat "$st/rates.out" 2> /dev/null)
      echo "$line" >> "$st/rates.jsonl"
      [ "$rc" = 0 ] || stop_unit "$unit" "rate guard (train_est.py rates exit $rc; 3: a judged step >20% over its estimate, DEC-43, DEC-53; other: the guard itself failed): $(echo "$line" | cut -c1-1000)"
    fi
    k=$((k + 1))
    sleep 5
  done
  RUNNING_UNIT=""
  line=$(python3 -I "$BIN/train_est.py" rates --root "$root" --release "$rel" --est "$est" --out "$st/rates-end.json" 2>> "$LOG")
  rc=$?
  log "rates at end of $unit (exit $rc): $(echo "$line" | cut -c1-900)"
  [ "$rc" = 0 ] || fail "rate guard at end of $unit (train_est.py rates exit $rc; 3: a judged step >20% over its estimate, DEC-43, DEC-53): $(echo "$line" | cut -c1-1000)"
  [ "$foreign" = 1 ] && return 3
  return 0
}

do_run() {
  mkdir -p "$WORK" "$LOGS" "$WORK/before" "$WORK/estimate" "$WORK/run" "$WORK/after" "$WORK/pack" "$WORK/check"
  # $WORK/DONE is the success marker of a run (card A-11): written last, only when every step passed. A waiter tests
  # it, never a log line or an exit status alone; a run withdraws the marker of an earlier one before it starts.
  rm -f "$WORK/DONE"
  # the profile copy is written once the run is known to be a first launch; a resume keeps the first launch's (A-15)
  mkdir -p "$WORK/tmp" "$LOGS"
  log "release_train run: root=$ROOT work=$WORK profile=$REALPROF"

  [ -e "$ROOT/STOP" ] && fail "$ROOT/STOP exists: $(head -c 300 "$ROOT/STOP")"
  stage_confs_present || fail "this tree has no stage confs ($BIN/stages/*.conf): its lock could not see a stage driver's build units (card A-12); run the train from a tree that has them"
  local act free
  act=$(build_units)
  [ -z "$act" ] || fail "active build unit(s): $act"
  free=$(free_gb)
  [ "$free" -ge "$START_MIN_GB" ] || fail "/tmp has ${free} GB free, under ${START_MIN_GB} GB"

  # 1. Plan (rules 1 & 2). Card A-15: the new plan is written beside the first launch's until the run is known to be a
  # first launch (it then replaces it) or a resume of the release this $WORK launched (it goes to $WORK/resume-<n>/,
  # and everything that follows reads the first launch's plan.json).
  plan_args
  local plan_json=$WORK/plan.json next_json=$WORK/plan.next.json next_txt=$WORK/plan.next.txt resume=0 rdir=""
  bb_ro py "$BIN/train_plan.py" "${PLAN_ARGS[@]}" --json > "$next_json" || exit $?
  bb_ro py "$BIN/train_plan.py" "${PLAN_ARGS[@]}" > "$next_txt" || { cat "$next_txt"; fail "train_plan.py failed"; }
  resume_state "$(cat "$next_json")" || fail "could not tell a resume from a first launch (the resume helper failed)"
  case "$RS_ACTION" in
    first)
      mv -f "$next_json" "$plan_json" || fail "mv $next_json"
      mv -f "$next_txt" "$WORK/plan.txt" || fail "mv $next_txt"
      make_profile
      cat "$WORK/plan.txt"
      ;;
    resume)
      resume=1
      rdir=$WORK/resume-$RS_UNIT
      mkdir -p "$rdir"
      mv -f "$next_json" "$rdir/plan.json" || fail "mv $next_json"
      mv -f "$next_txt" "$rdir/plan.txt" || fail "mv $next_txt"
      action_line "Action: resume $RS_REL ($RS_RECORDS run records published; unit $RS_UNIT)" < "$rdir/plan.txt"
      log "resume $RS_REL ($RS_HOW): $RS_RECORDS run records published; unit ${UNIT_PREFIX}-rel-$RS_UNIT"
      make_bs_args "$plan_json"
      resume_checks "$RS_REL" || fail "resume of $RS_REL refused, nothing launched: $RS_PROBLEMS"
      ;;
    *) fail "refused, nothing changed: $RS_REASON" ;;
  esac

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

  # 2. Before: verify every earlier manifest release & snapshot listing. A resume keeps the first launch's listing and
  # verify logs (taken before anything of the release was written) and verifies again into its own directory: each
  # manifest hash must equal the first launch's.
  local r ok=1 h0 h1
  if [ "$resume" = 1 ]; then
    for r in $earlier_rels; do
      verify_rel "$ROOT" "$r" "$rdir" || ok=0
      h0=$(hash_of "$WORK/before/verify-$r.log")
      h1=$(hash_of "$rdir/verify-$r.log")
      [ -n "$h1" ] && [ "$h0" = "$h1" ] || { log "resume: $r hash '$h1' != the first launch's '$h0'"; ok=0; }
    done
    [ "$ok" = 1 ] || fail "an earlier release did not verify before the resume, or its manifest hash is not the first launch's (see $rdir/verify-*.log and $WORK/before/verify-*.log)"
    log "resume: keeps the first launch's before/listing.tsv ($(wc -l < "$WORK/before/listing.tsv") files), its verify logs and estimate/before.json"
  else
    for r in $earlier_rels; do
      verify_rel "$ROOT" "$r" "$WORK/before" || ok=0
    done
    [ "$ok" = 1 ] || fail "an earlier release did not verify before launch (see $WORK/before/verify-*.log)"
    listing "$ROOT" "$WORK/before/listing.tsv"
    log "before: $(wc -l < "$WORK/before/listing.tsv") files in $ROOT"
  fi

  # 3. Estimate before launch (a resume is judged against the first launch's)
  if [ "$resume" != 1 ]; then
    local est_args=(before --root "$ROOT" --release "$rel" --stages "$stages" --profile "$PROF_COPY" --tmp "$WORK/tmp/est")
    [ -n "$latest" ] && est_args+=(--earlier "$latest")
    [ -n "$PARTS" ] && est_args+=(--parts "$PARTS")
    bb_ro py "$BIN/train_est.py" "${est_args[@]}" > "$WORK/estimate/before.json" 2> "$WORK/estimate/before.err" \
      || fail "train_est.py before failed (see $WORK/estimate/before.err)"
  fi

  # 4. Build release steps via build-stage (stage steps for last_stage all skip; on a resume the release's published
  # steps skip too). The unit's name is new to this $WORK (A-15): its log holds no other unit's exit line.
  [ "$resume" = 1 ] || make_bs_args "$plan_json"
  local bs_args=("${BS_ARGS[@]}")

  local unit="${UNIT_PREFIX}-rel-$RS_UNIT"
  [ -e "$LOGS/$unit.log" ] && fail "$LOGS/$unit.log exists: another unit's exit line could be read as $unit's"
  mkdir -p "$WORK/run/$unit"
  touch "$WORK/run/$unit/since"
  local ro_binds=()
  [ -d "$REAL_ROOT" ] && [ "$ROOT" != "$REAL_ROOT" ] && ro_binds+=(--bind "$REAL_ROOT:$REAL_ROOT:ro")
  [ -d "$PFAMDIR" ] && ro_binds+=(--bind "$PFAMDIR:$PFAMDIR:ro")
  log "launch $unit: build-stage --profile $PROF_COPY --root $ROOT ${bs_args[*]}"
  RUNNING_UNIT=$unit  # the trap stops it if the train itself is stopped while it runs (card A-14)
  RUNNING_LOGS=$LOGS
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
  [ "$RS_UNIT" = 1 ] || punit="${UNIT_PREFIX}-pack-$RS_UNIT"
  [ -e "$LOGS/$punit.log" ] && fail "$LOGS/$punit.log exists: another unit's exit line could be read as $punit's"
  [ -f "$PFAM" ] && pack_model=(--model "pfam-38.2=$PFAM")
  bash "$TREE/scripts/pod_run.sh" --profile "$PROF_COPY" --image "$IMAGE" --digest "$DIGEST" --mem "$PACK_MEM" \
    --threads "$PACK_JOBS" --bind /tmp/if2453 --bind "$ROOT:$ROOT:ro" "${ro_binds[@]}" --env "PYTHONPATH=$PYPATH" \
    --log "$LOGS/$punit.log" "$punit" pack build --root "$ROOT" --release "$rel" "${pack_model[@]}" \
    --profile "$PROF_COPY" --parts "$PACK_PARTS" --jobs "$PACK_JOBS" >> "$LOG" 2>&1
  sleep 2
  while :; do  # the watch's rule (A-14): its own exit line or a real inactive answer ends the wait, a failed answer not
    unit_ended "$punit" "$LOGS"
    rc=$?
    [ "$rc" = 0 ] && break
    watch_answer "$rc" "$punit"
    sleep 5
  done
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
