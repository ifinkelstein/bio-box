#!/usr/bin/env bash
# Render every diagrams/*.mmd to SVG with mermaid-cli, using the system Chrome.
# Usage: bash presentations/render.sh   (from anywhere)
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
d="$here/diagrams"
fail=0
for f in "$d"/d*.mmd; do
  if ! mmdc -q -i "$f" -o "${f%.mmd}.svg" \
       -c "$d/mermaid-config.json" -p "$d/puppeteer-config.json" -b transparent; then
    echo "FAIL $f" >&2
    fail=1
  fi
done
exit "$fail"
