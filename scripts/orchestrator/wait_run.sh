#!/bin/sh
# Wait until one GitHub Actions run is completed, then print its conclusion and each job's failed steps.
# `gh run watch` can return before a run has finished (seen 2026-11-12 and 2026-11-13), so this polls the status.
# Usage: wait_run.sh <repo> <run id> [interval_s]
# Exit 0 on success, 1 on any other conclusion, 4 when the run is not completed after 4 hours.
REPO=$1
ID=$2
INTERVAL=${3:-120}
w=0
st=""
while [ $((w * INTERVAL)) -lt 14400 ]; do
  st=$(gh run view "$ID" --repo "$REPO" --json status --jq .status 2>/dev/null)
  [ "$st" = "completed" ] && break
  w=$((w + 1))
  sleep "$INTERVAL"
done
if [ "$st" != "completed" ]; then
  echo "run $ID NOT COMPLETED after 4 hours (status '$st')"
  exit 4
fi
c=$(gh run view "$ID" --repo "$REPO" --json conclusion --jq .conclusion 2>/dev/null)
echo "run $ID: $c"
gh run view "$ID" --repo "$REPO" --json jobs \
  --jq '.jobs[] | .name + ": " + .conclusion + " / failed steps: " + ([.steps[] | select(.conclusion=="failure") | .name] | join(", "))'
[ "$c" = "success" ]
