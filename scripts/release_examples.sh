#!/bin/bash
# Run every worked example under examples/ on one release and gather their checks (card E-11, AGENTS 5.7). Run ON THE
# POD HOST (not in the image), before a milestone release is announced (RUNBOOK section 6):
#
#   TREE=... PROFILE=... IMAGE=... DIGEST=sha256:... ROOT=... OUT=/tmp/if2453/<you>/examples-rN.3 \
#     bash "$TREE/scripts/release_examples.sh" [--dry-run | --check-only] RELEASE
#
# The environment is each example's (examples/README.md): TREE (a checkout of the code, also the image's
# PYTHONPATH), PROFILE, IMAGE (the .sif, with matplotlib), DIGEST (sha256: of the image), ROOT (the data root, bound
# read-only), OUT (a new or empty directory outside ROOT), RELEASE (optional; when set it must be the argument).
# Optional, passed on to every example: THREADS (2), CAP, BIND, MEM, PFAM; UNIT (unit name prefix, default rex; an
# example's units are named UNIT-<example>-<step>-<pid>).
#
# Every directory examples/<name>/ with a run.sh and an expected.tsv is an example (examples/_shared/ is not). Each
# runs in turn, read-only on ROOT, into OUT/<name>/, each of its launches waiting for room under the account's thread
# cap (examples/_shared/launch.sh). Then the checking step: the shared check.py on each example's
# OUT/<name>/observed.tsv against its expected.tsv and recorded/RELEASE.tsv (an example that left no observed.tsv
# is checked against an empty one, so its invariants are missing), and one summary:
#
#   OUT/summary.tsv   example, check, kind, expected, observed, verdict (ok, FAIL, differs, new, missing, unexpected)
#   OUT/note.md       a measurement-note skeleton: the run's facts, per-example counts and costs, every check, and
#                     TODO where a person writes the reading
#   OUT/checks/<name>.txt, OUT/checks/<name>.tsv   check.py's lines and table per example
#   OUT/logs/<name>.out, <name>.err                each run.sh's own output
#
# --dry-run prints what would run and writes nothing. --check-only runs only the checking step, on an OUT that
# already holds OUT/<name>/observed.tsv (a run of this script, or files made by hand); there LOCAL_PYTHON, when set,
# runs check.py with that interpreter instead of in the image (the tests do; on the pod leave it unset).
#
# Exit: 0 when every invariant of every example holds; 1 when one fails or is missing (a failed invariant stops the
# release's announcement); 2 on a usage error, a malformed expected file, or a STOP file in ROOT.
set -euo pipefail

ME=release_examples.sh
usage() { echo "usage: $ME [--dry-run | --check-only] RELEASE (environment: TREE PROFILE IMAGE DIGEST ROOT OUT)" >&2; exit 2; }

mode=run
case "${1:-}" in
  --dry-run) mode=dry; shift ;;
  --check-only) mode=check; shift ;;
  -*) usage ;;
esac
[ $# -eq 1 ] || usage
if [ -n "${RELEASE:-}" ] && [ "$RELEASE" != "$1" ]; then
  echo "$ME: RELEASE is $RELEASE in the environment, not $1" >&2
  exit 2
fi
RELEASE=$1
: "${TREE:?set TREE}" "${OUT:?set OUT}"
if [ "$mode" != check ] || [ -z "${LOCAL_PYTHON:-}" ]; then
  : "${PROFILE:?set PROFILE}" "${IMAGE:?set IMAGE}" "${DIGEST:?set DIGEST}" "${ROOT:?set ROOT}"
fi
ROOT=${ROOT:-} PROFILE=${PROFILE:-} IMAGE=${IMAGE:-} DIGEST=${DIGEST:-}
THREADS=${THREADS:-2}
BIND=${BIND:-/tmp/if2453}
UNIT=${UNIT:-rex}
MEM=${MEM:-16G}
SHARED="$TREE/examples/_shared"
CHECK="$SHARED/check.py"
export TREE PROFILE IMAGE DIGEST ROOT RELEASE THREADS BIND MEM
[ -z "${CAP:-}" ] || export CAP
[ -z "${PFAM:-}" ] || export PFAM
[ -f "$SHARED/launch.sh" ] || { echo "$ME: TREE has no $SHARED/launch.sh" >&2; exit 2; }
. "$SHARED/launch.sh"

names=()
for d in "$TREE"/examples/*/; do
  name=$(basename "$d")
  case "$name" in _*) continue ;; esac
  [ -f "$d/run.sh" ] && [ -f "$d/expected.tsv" ] || continue
  names+=("$name")
done
[ "${#names[@]}" -gt 0 ] || die "no example under $TREE/examples"

if [ "$mode" = dry ]; then
  echo "$ME: release $RELEASE, ${#names[@]} examples, root $ROOT (read-only), out $OUT, $THREADS threads per launch"
  for name in "${names[@]}"; do
    echo "run    $name: OUT=$OUT/$name UNIT=$UNIT-$name bash $TREE/examples/$name/run.sh" \
      "> $OUT/logs/$name.out 2> $OUT/logs/$name.err"
  done
  for name in "${names[@]}"; do
    echo "check  $name: python -I $CHECK --expected $TREE/examples/$name/expected.tsv" \
      "--recorded $TREE/examples/$name/recorded/$RELEASE.tsv --table $OUT/checks/$name.tsv $OUT/$name/observed.tsv"
  done
  echo "write  $OUT/summary.tsv and $OUT/note.md"
  exit 0
fi

case "$OUT/" in "$ROOT"/*) [ -z "$ROOT" ] || die "OUT is under ROOT" ;; esac
[ ! -e "$OUT/summary.tsv" ] || die "$OUT already has a summary (summary.tsv); give a new OUT"
if [ "$mode" = run ]; then
  [ ! -e "$ROOT/STOP" ] || die "a STOP file exists in $ROOT; not starting (AGENTS 7.4)"
  if [ -e "$OUT" ] && [ -n "$(ls -A "$OUT")" ]; then
    die "OUT $OUT is not empty"
  fi
fi
mkdir -p "$OUT/logs" "$OUT/checks"
STARTED=$(date '+%F %T %Z')

declare -A ran
if [ "$mode" = run ]; then
  for name in "${names[@]}"; do
    [ ! -e "$ROOT/STOP" ] || die "a STOP file appeared in $ROOT; $name and the examples after it were not run"
    note "$name: starting (OUT=$OUT/$name)"
    set +e
    OUT="$OUT/$name" UNIT="$UNIT-$name" bash "$TREE/examples/$name/run.sh" > "$OUT/logs/$name.out" \
      2> "$OUT/logs/$name.err"
    ran[$name]=$?
    set -e
    note "$name: run.sh exited ${ran[$name]}"
  done
fi

py() {  # py LOG ARGS...: python -I ARGS, in the image (a capped scope) unless LOCAL_PYTHON names an interpreter
  local log=$1
  shift
  if [ -n "${LOCAL_PYTHON:-}" ]; then
    "$LOCAL_PYTHON" -I "$@"
  else
    tool "$log" python -I "$@"
  fi
}

# The checking step.
worst=0
declare -A rc
for name in "${names[@]}"; do
  ex="$TREE/examples/$name"
  obs="$OUT/$name/observed.tsv"
  if [ ! -f "$obs" ]; then
    obs="$OUT/checks/$name.none.tsv"
    printf 'check_id\tvalue\n' > "$obs"
  fi
  set +e
  py "$OUT/logs/$name.check.time" "$CHECK" --expected "$ex/expected.tsv" --recorded "$ex/recorded/$RELEASE.tsv" \
    --table "$OUT/checks/$name.tsv" "$obs" > "$OUT/checks/$name.txt"
  rc[$name]=$?
  set -e
  [ "${rc[$name]}" -le "$worst" ] || worst=${rc[$name]}
done

# One summary over the examples.
printf 'example\tcheck\tkind\texpected\tobserved\tverdict\n' > "$OUT/summary.tsv"
for name in "${names[@]}"; do
  [ -f "$OUT/checks/$name.tsv" ] || continue
  awk -v ex="$name" 'NR > 1 { print ex "\t" $0 }' "$OUT/checks/$name.tsv" >> "$OUT/summary.tsv"
done
bad=$(awk -F'\t' 'NR > 1 && $3 == "invariant" && $6 != "ok" { n++ } END { print n + 0 }' "$OUT/summary.tsv")
per=""
for name in "${names[@]}"; do
  n=$(awk -F'\t' -v ex="$name" 'NR > 1 && $1 == ex && $3 == "invariant" && $6 != "ok" { n++ } END { print n + 0 }' \
    "$OUT/summary.tsv")
  per="$per${per:+, }$name: $n"
done

cpu_s() {  # CPU-seconds of an example's launches: the units' CPUUsageNSec plus the scopes' user and system time
  local d=$1
  awk '
    /^pod_run: .* usage / { for (i = 1; i <= NF; i++) if ($i ~ /^CPUUsageNSec=/) { split($i, a, "="); s += a[2] / 1e9 } }
    /User time \(seconds\):|System time \(seconds\):/ { s += $NF }
    END { printf "%.1f", s }' $(find "$d" -maxdepth 2 \( -name '*.log' -o -name '*.time' \) 2> /dev/null | sort) /dev/null
}

{
  echo "# Worked examples on \`$RELEASE\` (release check, AGENTS 5.7)"
  echo
  echo "*Measured*, $(hostname -s), started $STARTED, by \`scripts/release_examples.sh $RELEASE\`"
  echo "($mode mode): tree \`$TREE\` (commit: TODO), image \`$IMAGE\` (\`$DIGEST\`), data root \`$ROOT\`"
  echo "(read-only), outputs \`$OUT\`, $THREADS threads per launch. TODO where a person writes the reading."
  echo
  echo "## Result"
  echo
  if [ "$worst" -eq 0 ]; then
    echo "Every invariant of the ${#names[@]} examples holds. TODO: one paragraph."
  else
    echo "**$bad invariants failed or are missing** ($per): the release is not announced (RUNBOOK section 6). TODO."
  fi
  echo
  echo "## Per example"
  echo
  echo "| Example | run.sh exit | check exit | ok | FAIL | differs | new | missing | unexpected | CPU-s (*derived*) |"
  echo "|---|---|---|---|---|---|---|---|---|---|"
  for name in "${names[@]}"; do
    counts=$(awk -F'\t' -v ex="$name" 'NR > 1 && $1 == ex { n[$6]++ }
      END { printf "%d | %d | %d | %d | %d | %d", n["ok"], n["FAIL"], n["differs"], n["new"], n["missing"], n["unexpected"] }' \
      "$OUT/summary.tsv")
    echo "| \`$name\` | ${ran[$name]:--} | ${rc[$name]} | $counts | $(cpu_s "$OUT/$name") |"
  done
  echo
  echo "CPU-s is the sum of the units' \`CPUUsageNSec\` (\`pod_run: ... usage\` lines) and of \`/usr/bin/time -v\` user"
  echo "plus system time of the scopes, from the example's logs."
  echo
  echo "## Differences from the recorded values"
  echo
  if awk -F'\t' 'NR > 1 && $6 == "differs" { f = 1 } END { exit !f }' "$OUT/summary.tsv"; then
    echo "| Example | Check | Recorded | Observed | Cause |"
    echo "|---|---|---|---|---|"
    awk -F'\t' 'NR > 1 && $6 == "differs" { print "| `" $1 "` | `" $2 "` | " $4 " | " $5 " | TODO |" }' "$OUT/summary.tsv"
  else
    echo "None."
  fi
  echo
  echo "## Every check"
  echo
  echo "| Example | Check | Kind | Expected | Observed | Verdict |"
  echo "|---|---|---|---|---|---|"
  awk -F'\t' 'NR > 1 { print "| `" $1 "` | `" $2 "` | " $3 " | " $4 " | " $5 " | " $6 " |" }' "$OUT/summary.tsv"
} > "$OUT/note.md"

if [ "$worst" -eq 0 ]; then
  echo "$ME: $RELEASE: ${#names[@]} examples, every invariant holds; summary $OUT/summary.tsv" >&2
else
  echo "$ME: $RELEASE: ${#names[@]} examples, $bad invariant(s) failed or missing ($per); summary $OUT/summary.tsv" >&2
fi
exit "$worst"
