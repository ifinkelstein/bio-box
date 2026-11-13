#!/bin/sh
# Run the `ci` workflow N times on one branch by re-running its latest push run; stop at the first failure.
# Usage: ci_repeat.sh <repo> <branch> <n>
REPO=$1
BRANCH=$2
N=$3
i=0
while [ $i -lt 30 ]; do
  id=$(gh run list --repo "$REPO" --workflow ci --branch "$BRANCH" --limit 1 --json databaseId --jq '.[0].databaseId' 2>/dev/null)
  [ -n "$id" ] && break
  i=$((i + 1))
  sleep 20
done
if [ -z "$id" ]; then echo "NO RUN FOUND for $BRANCH"; exit 2; fi
k=1
while [ $k -le "$N" ]; do
  # `gh run watch` can return before the run has finished (seen 2026-11-12, card C-16): poll the status itself
  # and read the conclusion only once it is completed; give up after 2 hours.
  w=0
  st=""
  while [ $w -lt 240 ]; do
    st=$(gh run view "$id" --repo "$REPO" --json status --jq .status 2>/dev/null)
    [ "$st" = "completed" ] && break
    w=$((w + 1))
    sleep 30
  done
  if [ "$st" != "completed" ]; then echo "attempt $k of $N: NOT COMPLETED after 2 hours (run $id, status '$st')"; exit 4; fi
  c=$(gh run view "$id" --repo "$REPO" --json conclusion --jq .conclusion 2>/dev/null)
  echo "attempt $k of $N: $c (run $id)"
  if [ "$c" != "success" ]; then
    gh run view "$id" --repo "$REPO" --log-failed 2>/dev/null | grep -E "FAILED|AssertionError|Error|passed|failed" | tail -n 20
    exit 1
  fi
  k=$((k + 1))
  [ $k -le "$N" ] || break
  gh run rerun "$id" --repo "$REPO" > /dev/null 2>&1 || { echo "RERUN REFUSED at attempt $k"; exit 3; }
  sleep 30
done
echo "ALL $N GREEN on $BRANCH"
