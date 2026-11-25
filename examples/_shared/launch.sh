# The launcher of the worked examples (card E-11): sourced, never run, by each `examples/<name>/run.sh` and by
# `scripts/release_examples.sh`, ON THE POD HOST. Written for the CRISPR-Cas example's run.sh and moved here
# unchanged but for the prefix of its messages and the optional Pfam bind.
#
#   . "$TREE/examples/_shared/launch.sh"
#
# Defines die, note, cap, room, search_mem_gb, bb and tool. They read, when called: ROOT (the data root, bound read-only), OUT
# (the caller's scratch directory; `bb` writes OUT/logs/NAME.log and NAME.launch), TREE, PROFILE, IMAGE, DIGEST,
# THREADS (threads of one launch), MEM (pod_run.sh --mem), SEARCH_MEM_GB (optional override for `search_mem_gb`),
# BIND (the directory bound read-write into the image,
# holding TREE, OUT, the pack directory and the image), UNIT (unit name prefix), CAP (optional: the account's thread
# cap; default 40 Mon-Fri 08:00-20:00, else 80 when no other login user uses > 0.5 CPUs; AGENTS 7.2, DEC-74, C-28),
# PFAM_DIR (optional: a directory bound read-only into `tool`, the Pfam file's) and ME (the name in messages; default
# run.sh). Nothing here writes under ROOT.

die() { echo "${ME:-run.sh}: $*" >&2; exit 2; }
note() { echo "${ME:-run.sh}: $(date '+%F %T') $*" >&2; }  # stderr: stdout of `tool` may be an output file

others_busy() {  # 0 (true) if other login users (uid > 999, != self) use > 0.5 CPUs (> 50% CPU) or ps fails (C-28)
  local my_uid=${LAUNCH_UID:-$(id -u)} ps_out
  if [ -n "${LAUNCH_PS_OUT+x}" ]; then
    ps_out=$LAUNCH_PS_OUT
  else
    ps_out=$(ps -eo uid=,pcpu= 2>/dev/null) || return 0
  fi
  awk -v me="$my_uid" '
    NF >= 2 && $1 ~ /^[0-9]+$/ && $1 + 0 > 999 && $1 + 0 != me + 0 { s += $2 + 0 }
    END { exit (s > 50.0 ? 0 : 1) }
  ' <<< "$ps_out"
}

cap() {  # the account's thread cap of this hour (AGENTS 7.2, DEC-74, C-28)
  if [ -n "${CAP:-}" ]; then echo "$CAP"; return; fi
  local h=${LAUNCH_HOUR:-$(date +%-H)}
  local d=${LAUNCH_DOW:-$(date +%u)}
  if [ "$d" -le 5 ] && [ "$h" -ge 8 ] && [ "$h" -lt 20 ]; then
    echo 40
  elif others_busy; then
    echo 40
  else
    echo 80
  fi
}

room() {  # wait until this account's CPU (ps %CPU, summed) plus THREADS fits under the cap
  local total v c
  while true; do
    [ ! -e "$ROOT/STOP" ] || die "a STOP file appeared in $ROOT"
    total=0
    for v in $(ps -u "$USER" -o pcpu=); do
      total=$((total + ${v%.*}))
    done
    c=$(cap)
    if [ $(((total + 99) / 100 + THREADS)) -le "$c" ]; then
      return
    fi
    note "waiting for room: account at ${total}% CPU, cap $c threads"
    sleep 60
  done
}

search_mem_gb() {  # derive `search hmm --mem-gb` from MEM (2G headroom, min 1); SEARCH_MEM_GB overrides (A-25b)
  if [ -n "${SEARCH_MEM_GB:-}" ]; then
    echo "$SEARCH_MEM_GB"
    return
  fi
  local raw=${MEM:-16G} gb
  if [[ "$raw" =~ ^([0-9]+)[Gg]([Ii]?[Bb]?)?$ ]]; then
    gb=$((10#${BASH_REMATCH[1]}))
  elif [[ "$raw" =~ ^([0-9]+)[Mm]([Ii]?[Bb]?)?$ ]]; then
    gb=$((10#${BASH_REMATCH[1]} / 1024))
  elif [[ "$raw" =~ ^([0-9]+)[Tt]([Ii]?[Bb]?)?$ ]]; then
    gb=$((10#${BASH_REMATCH[1]} * 1024))
  else
    echo 14
    return
  fi
  if [ "$gb" -gt 2 ]; then
    echo $((gb - 2))
  else
    echo 1
  fi
}

bb() {  # bb NAME ARGS...: one bio-box command in its own unit through the runbook's launcher; waits for it
  local name=$1
  shift
  local unit="$UNIT-$name-$$"
  local log="$OUT/logs/$name.log"
  room
  bash "$TREE/scripts/pod_run.sh" --profile "$PROFILE" --image "$IMAGE" --digest "$DIGEST" --mem "$MEM" \
    --threads "$THREADS" --bind "$BIND" --bind "$ROOT:$ROOT:ro" --env PYTHONPATH="$TREE/src" \
    --env TMPDIR="$OUT/tmp" --log "$log" "$unit" "$@" > "$OUT/logs/$name.launch"
  until grep -q "^pod_run: $unit exited" "$log"; do
    sleep 5
  done
  grep -q "^pod_run: $unit exited 0$" "$log" || die "$name failed; see $log"
  note "$name: $(grep "^pod_run: $unit usage" "$log" || true)"
}

tool() {  # tool LOG ARGS...: one non-bio-box command of the image (HMMER, python) in a capped scope
  local log=$1
  shift
  room
  systemd-run --user --scope --quiet -p MemoryMax=8G -p MemorySwapMax=0 -p CPUQuota=$((THREADS * 100))% \
    nice -n 19 /usr/bin/time -v -a -o "$log" singularity exec --cleanenv --bind "$BIND" --bind "$ROOT:$ROOT:ro" \
    ${PFAM_DIR:+--bind "$PFAM_DIR:$PFAM_DIR:ro"} "$IMAGE" "$@"
}
