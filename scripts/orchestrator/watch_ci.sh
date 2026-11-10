#!/bin/sh
# Wait for the workflows of one commit to finish; print each conclusion and the fixture manifest hash.
# Usage: watch_ci.sh <repo dir> <sha>
cd "$1" || exit 2
SHA=$2
i=0
while [ $i -lt 90 ]; do
  out=$(gh run list --commit "$SHA" --json name,status,conclusion,databaseId 2>/dev/null) || out=""
  n=$(printf '%s' "$out" | python3 -c 'import json,sys; print(len(json.load(sys.stdin)))' 2>/dev/null || echo 0)
  pending=$(printf '%s' "$out" | python3 -c 'import json,sys; d=json.load(sys.stdin); print(sum(1 for r in d if r["status"]!="completed"))' 2>/dev/null || echo 1)
  if [ "$n" != "0" ] && [ "$pending" = "0" ]; then
    printf '%s' "$out" | python3 -c 'import json,sys
for r in json.load(sys.stdin): print(r["name"], r["conclusion"], r["databaseId"])'
    id=$(printf '%s' "$out" | python3 -c 'import json,sys
print(next((r["databaseId"] for r in json.load(sys.stdin) if r["name"]=="cross-host"), ""))')
    if [ -n "$id" ]; then
      gh run view "$id" --log 2>/dev/null | grep -o '"manifest_hash": *"[0-9a-f]*"' | sort -u
      gh run view "$id" --log-failed 2>/dev/null | tail -n 25
    fi
    exit 0
  fi
  i=$((i + 1))
  sleep 60
done
echo "CI TIMEOUT after 90 min for $SHA"
exit 3
