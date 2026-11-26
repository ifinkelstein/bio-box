#!/bin/bash
# Copy files from the PI's NAS `tower` to the pod, streaming through this laptop (card C-73): the mirror of
# relay_to_tower.sh (card C-49) with the roles swapped. The pod and tower do not see each other, so the laptop relays;
# it has little free disk, so nothing lands on it: each file goes `ssh tower cat` | `ssh pod tee`.
#
#   relay_to_pod.sh [options] SRC_DIR DST_DIR LIST
#
# SRC_DIR  absolute directory on tower (read only: `stat` and `cat`; nothing on tower is written or removed)
# DST_DIR  absolute directory on the pod (made if absent), strictly under /tmp/if2453/ or, for files of 64 MiB and
#          more only, under /stor/scratch/Wilke/if2453/bio-box-staging/ (AGENTS §7.3). Anything else: exit 2.
# LIST     a local file, one file per line: `<path relative to SRC_DIR><TAB><sha256>`, or the SHA256SUMS form
#          `<sha256>  <path>` (`sha256sum` output, `*` marker allowed), such as tower's SHA256SUMS. Blank lines and
#          lines starting `#` are skipped. Every path and checksum is checked before the first remote call; a path
#          may hold only letters, digits and . _ = + , - in components joined by /, no component . or .. or starting
#          with -; the names SHA256SUMS, STOP and *.partial are the relay's own. Anything else: exit 2.
#
# The plan reads the sizes on tower, then what the pod holds. For DST_DIR under /stor, a file to copy of less than
# 64 MiB (67,108,864 bytes) refuses the whole run (exit 2) before any call to the pod.
# Per file: a file already at DST_DIR/<path> is skipped and never overwritten. Otherwise tower's bytes are streamed
# into DST_DIR/<dir>/.<name>.partial while the pod hashes the same stream (`tee .partial | sha256sum`); the hash must
# equal the record and the .partial's size tower's, then `mv -n` puts it in place (refused if the name appeared
# meanwhile). A failed try removes its .partial; a .partial left by a killed run is removed and copied again.
# A STOP file in DST_DIR (AGENTS §7.4) ends the run before the next file: nothing more is streamed, exit 1.
# At the end, for DST_DIR under /tmp/if2453/ only, DST_DIR/SHA256SUMS is written or extended with the files verified
# in this run (copied, or read back with --reverify); a file in place that was not verified is not added. A line
# already there with another sha256 is a conflict: reported, and SHA256SUMS is left as it was. Under /stor no
# SHA256SUMS is written (a file under 64 MiB); the log holds each file's sha256.
#
# Options:
#   --streams N  files copied at once (default 1: two streams measured 8.4 MB/s together against 9.70 for one in the
#                other direction, card C-49, so a second stream is not expected to add)
#   --reverify   also read every file back on the pod and hash it again before it is put in place, and hash the
#                files found in place (doubles the pod's reads)
#   --dry-run    print what would be copied and the bytes; no write anywhere, no log
#   --tries N    tries per file (default 3); --pause S seconds between tries (default 30)
#   --log FILE   append the log here (default ${XDG_STATE_HOME:-$HOME/.local/state}/bio-box/relay_to_pod.log)
#
# Environment (each checked: letters, digits and . _ @ : - only):
#   BIOBOX_POD             user@host of the pod (default if2453@wilkcomp02.ccbb.utexas.edu)
#   BIOBOX_POD_JUMP        user@host of a jump host for the pod, off campus (as scripts/pod_pytest.sh, card A-09)
#   BIOBOX_TOWER           user@host of the NAS (default root@tower; the laptop resolves the name)
#   BIOBOX_TOWER_IDENTITY  the laptop's key for tower (default $HOME/.ssh/id_ed25519); no key is made or copied
#   BIOBOX_RELAY_SSH_INCLUDE  an ssh config included first when it exists (default $HOME/.ssh/config), for host
#                          names and keys the laptop already uses
#
# Exit: 0 every file in place and verified or skipped; 1 a file failed, was missing on tower, a STOP file ended the
# run, or a conflict in SHA256SUMS; 2 a refusal before any write (the list, the directories and the hosts before any
# remote call; a file under 64 MiB bound for /stor before any call to the pod). Run it again after a failure or a
# stop: what is in place is skipped.
set -u
export LC_ALL=C

usage() { sed -n '2,/^set -u/p' "$0" | sed '$d' | sed 's/^# \{0,1\}//'; }
die() { echo "relay_to_pod: $*" >&2; exit 2; }

STREAMS=1 DRY=0 REVERIFY=0 TRIES=3 PAUSE=30
LOG=${XDG_STATE_HOME:-${HOME:-/tmp}/.local/state}/bio-box/relay_to_pod.log
while [ $# -gt 0 ]; do
  case "$1" in
    --streams) [ $# -ge 2 ] || die "--streams needs a value"; STREAMS=$2; shift 2 ;;
    --dry-run) DRY=1; shift ;;
    --reverify) REVERIFY=1; shift ;;
    --tries) [ $# -ge 2 ] || die "--tries needs a value"; TRIES=$2; shift 2 ;;
    --pause) [ $# -ge 2 ] || die "--pause needs a value"; PAUSE=$2; shift 2 ;;
    --log) [ $# -ge 2 ] || die "--log needs a value"; LOG=$2; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    --) shift; break ;;
    -*) die "unknown option $1 (see --help)" ;;
    *) break ;;
  esac
done
[ $# -eq 3 ] || die "usage: relay_to_pod.sh [options] SRC_DIR DST_DIR LIST (see --help)"
SRC=${1%/} DST=${2%/} LIST=$3
case "$STREAMS" in ''|*[!0-9]*) die "--streams must be a whole number" ;; esac
[ "$STREAMS" -ge 1 ] && [ "$STREAMS" -le 8 ] || die "--streams must be 1 to 8"
case "$TRIES" in ''|*[!0-9]*) die "--tries must be a whole number" ;; esac
[ "$TRIES" -ge 1 ] || die "--tries must be at least 1"
case "$PAUSE" in ''|*[!0-9]*) die "--pause must be a whole number of seconds" ;; esac

# --- checks, all before the first remote call --------------------------------------------------------------------
COMP='[A-Za-z0-9._=+,-]+'
PATH_RE="^$COMP(/$COMP)*\$"
ABS_RE="^(/$COMP)+\$"
HOST_RE='^[A-Za-z0-9._@:-]+$'
SHA_RE='^[0-9a-f]{64}$'
POD_TMP=/tmp/if2453
POD_STOR=/stor/scratch/Wilke/if2453/bio-box-staging
MIB64=67108864

safe_components() {  # no component . or .. or starting with -
  case "/$1/" in */./*|*/../*|*/-*) return 1 ;; esac
  return 0
}
[[ $SRC =~ $ABS_RE ]] && safe_components "$SRC" || die "SRC_DIR must be an absolute path of safe characters: $SRC"
[[ $DST =~ $ABS_RE ]] && safe_components "$DST" || die "DST_DIR must be an absolute path of safe characters: $DST"
case "$DST" in
  "$POD_TMP"/?*) ON_STOR=0 ;;
  "$POD_STOR"/?*) ON_STOR=1 ;;
  *) die "DST_DIR must lie under $POD_TMP/ or $POD_STOR/ (AGENTS §7.3): $DST" ;;
esac
POD=${BIOBOX_POD:-if2453@wilkcomp02.ccbb.utexas.edu}
JUMP=${BIOBOX_POD_JUMP:-}
TOWER=${BIOBOX_TOWER:-root@tower}
IDENT=${BIOBOX_TOWER_IDENTITY:-${HOME:-}/.ssh/id_ed25519}
INCLUDE=${BIOBOX_RELAY_SSH_INCLUDE-${HOME:-}/.ssh/config}
[[ $POD =~ $HOST_RE ]] || die "BIOBOX_POD must be one user@host"
[ -z "$JUMP" ] || [[ $JUMP =~ $HOST_RE ]] || die "BIOBOX_POD_JUMP must be one user@host"
[[ $TOWER =~ $HOST_RE ]] || die "BIOBOX_TOWER must be one user@host"
[[ $IDENT =~ $ABS_RE ]] || die "BIOBOX_TOWER_IDENTITY must be an absolute path of safe characters"
[ -z "$INCLUDE" ] || [[ $INCLUDE =~ $ABS_RE ]] || die "BIOBOX_RELAY_SSH_INCLUDE must be an absolute path of safe characters"
[ -f "$LIST" ] || die "no list file: $LIST"

REL=() SHA=()
n=0
tab=$(printf '\t')
while IFS= read -r line || [ -n "$line" ]; do
  n=$((n + 1))
  case "$line" in ''|'#'*) continue ;; esac
  case "$line" in
    *"$tab"*"$tab"*) die "list line $n: more than two fields" ;;
    *"$tab"*) rel=${line%%"$tab"*}; sum=${line#*"$tab"} ;;
    *) sum=${line%% *}; rest=${line#* }
       case "$rest" in ' '*|'*'*) rel=${rest#?} ;; *) die "list line $n: neither <path><TAB><sha256> nor <sha256>  <path>" ;; esac ;;
  esac
  [[ $sum =~ $SHA_RE ]] || die "list line $n: the checksum is not 64 lower-case hex digits"
  [[ $rel =~ $PATH_RE ]] && safe_components "$rel" || die "list line $n: refused path (only letters, digits, . _ = + , - and /; no . or .. or leading -)"
  case "$rel" in SHA256SUMS|.SHA256SUMS*|STOP|*.partial) die "list line $n: the name $rel is the relay's own" ;; esac
  REL+=("$rel") SHA+=("$sum")
done < "$LIST"
N=${#REL[@]}
[ "$N" -gt 0 ] || die "the list names no file"
dup=$(printf '%s\n' "${REL[@]}" | sort | uniq -d | head -n 1)
[ -z "$dup" ] || die "the list names $dup twice"

# --- the ssh config: hosts and the identity from the environment, nothing else ------------------------------------
TMPD=$(mktemp -d "${TMPDIR:-/tmp}/relay.XXXXXX") || die "cannot make a temporary directory"
trap 'rm -rf "$TMPD"' EXIT
trap 'echo "relay_to_pod: stopped by a signal" >&2; [ "$DRY" = 1 ] || say "stopped by a signal"; kill $(jobs -p) 2>/dev/null; exit 130' INT TERM
CFG=$TMPD/ssh_config
{
  [ -n "$INCLUDE" ] && [ -f "$INCLUDE" ] && echo "Include $INCLUDE"
  echo "Host pod"
  echo "  HostName ${POD#*@}"
  case "$POD" in *@*) echo "  User ${POD%%@*}" ;; esac
  [ -z "$JUMP" ] || echo "  ProxyJump $JUMP"
  echo "  BatchMode yes"
  echo "  ServerAliveInterval 30"
  echo "Host tower"
  echo "  HostName ${TOWER#*@}"
  case "$TOWER" in *@*) echo "  User ${TOWER%%@*}" ;; esac
  echo "  IdentityFile $IDENT"
  echo "  IdentitiesOnly yes"
  echo "  BatchMode yes"
  echo "  ServerAliveInterval 30"
} > "$CFG"

pod() { ssh -n -F "$CFG" pod "$1"; }        # no stdin
tower() { ssh -n -F "$CFG" tower "$1"; }    # no stdin
say() { mkdir -p "$(dirname "$LOG")" && echo "$(date '+%F %T') $*" >> "$LOG"; }
mbs() { awk -v b="$1" -v s="$2" 'BEGIN { if (s > 0) printf "%.2f", b / 1e6 / s; else printf "n/a" }'; }
stopped() {  # true unless the pod answers that DST_DIR holds no STOP file
  [ "$(pod "if [ -e '$DST/STOP' ]; then echo stop; else echo go; fi")" != go ]
}

# --- the plan: sizes on tower, then what the pod holds, two calls ---------------------------------------------------
printf '%s\n' "${REL[@]}" > "$TMPD/paths"
ssh -F "$CFG" tower "cd '$SRC' || exit 3; while IFS= read -r p; do if [ -f \"\$p\" ] && [ -r \"\$p\" ]; then stat -c %s -- \"\$p\"; else echo missing; fi; done" < "$TMPD/paths" > "$TMPD/sizes"
rc=$?
[ "$rc" = 0 ] || { echo "relay_to_pod: reading sizes on tower failed (exit $rc; 3: no directory $SRC)" >&2; exit 1; }
SIZE=()
while IFS= read -r v; do SIZE+=("$v"); done < "$TMPD/sizes"
[ "${#SIZE[@]}" = "$N" ] || { echo "relay_to_pod: the plan's answers do not match the list" >&2; exit 1; }
for i in $(seq 0 $((N - 1))); do
  case "${SIZE[$i]}" in missing) ;; ''|*[!0-9]*) echo "relay_to_pod: tower gave no size for ${REL[$i]}" >&2; exit 1 ;; esac
done
# a file under 64 MiB bound for /stor refuses the run before any call to the pod
if [ "$ON_STOR" = 1 ]; then
  small=0
  for i in $(seq 0 $((N - 1))); do
    if [ "${SIZE[$i]}" != missing ] && [ "${SIZE[$i]}" -lt "$MIB64" ]; then
      echo "relay_to_pod: refused ${REL[$i]} (${SIZE[$i]} bytes): under 64 MiB, bound for /stor (AGENTS §7.3)" >&2
      small=$((small + 1))
    fi
  done
  [ "$small" = 0 ] || die "$small file(s) under 64 MiB bound for $DST; send them under $POD_TMP/ instead"
fi
ssh -F "$CFG" pod "if [ -e '$DST/STOP' ]; then echo stop; else echo go; fi; if cd '$DST' 2>/dev/null; then while IFS= read -r p; do if [ -e \"\$p\" ]; then echo yes; else echo no; fi; done; else while IFS= read -r p; do echo no; done; fi" < "$TMPD/paths" > "$TMPD/has"
rc=$?
[ "$rc" = 0 ] || { echo "relay_to_pod: reading $DST on the pod failed (exit $rc)" >&2; exit 1; }
HAS=()
STOP_AT_START=
{ IFS= read -r STOP_AT_START; while IFS= read -r v; do HAS+=("$v"); done; } < "$TMPD/has"
[ "${#HAS[@]}" = "$N" ] && { [ "$STOP_AT_START" = go ] || [ "$STOP_AT_START" = stop ]; } \
  || { echo "relay_to_pod: the plan's answers do not match the list" >&2; exit 1; }

todo=0 todo_bytes=0 inplace=0 missing=0
for i in $(seq 0 $((N - 1))); do
  if [ "${HAS[$i]}" = yes ]; then inplace=$((inplace + 1))
  elif [ "${SIZE[$i]}" = missing ]; then missing=$((missing + 1))
  else todo=$((todo + 1)); todo_bytes=$((todo_bytes + SIZE[i])); fi
done

if [ "$DRY" = 1 ]; then
  for i in $(seq 0 $((N - 1))); do
    if [ "${HAS[$i]}" = yes ]; then echo "in place ${REL[$i]}"
    elif [ "${SIZE[$i]}" = missing ]; then echo "missing on tower ${REL[$i]}"
    else echo "copy ${REL[$i]} ${SIZE[$i]}"; fi
  done
  took=$(awk -v b="$todo_bytes" 'BEGIN { s = b / 9.70e6; if (s < 120) printf "%.0f s", s; else printf "%.1f minutes", s / 60 }')
  echo "dry run: $todo files to copy, $todo_bytes bytes, about $took at 9.70 MB/s (estimate: one stream pod to tower, measured 2026-11-20, card C-49; this direction not measured); $inplace in place; $missing missing on tower"
  [ "$STOP_AT_START" = go ] || echo "dry run: a STOP file is in $DST on the pod; a run would copy nothing"
  [ "$missing" = 0 ] && [ "$STOP_AT_START" = go ] || exit 1
  exit 0
fi

say "start (pid $$): $N files, $todo to copy ($todo_bytes bytes), $inplace in place, $missing missing on tower; $STREAMS stream(s); tower:$SRC -> pod:$DST${JUMP:+ (pod through $JUMP)}"
if [ "$STOP_AT_START" != go ]; then
  say "STOP file in $DST on the pod: nothing copied (remove it and run again)"
  tail -n 2 "$LOG"
  exit 1
fi
for i in $(seq 0 $((N - 1))); do
  [ "${SIZE[$i]}" != missing ] || [ "${HAS[$i]}" = yes ] || say "${REL[$i]}: FAILED: not a readable file on tower"
done

# --- one file -----------------------------------------------------------------------------------------------------
one() {  # one <index> <result file>
  local i=$1 out=$2 rel=${REL[$1]} want=${SHA[$1]} size=${SIZE[$1]} dir base part dest try got got2 t0 t1 st psize
  dest="$DST/$rel"
  case "$rel" in */*) dir="$DST/${rel%/*}"; base=${rel##*/} ;; *) dir=$DST; base=$rel ;; esac
  part="$dir/.$base.partial"
  if [ "${HAS[$i]}" = yes ]; then
    if [ "$REVERIFY" = 1 ]; then
      got=$(pod "sha256sum -- '$dest'"); got=${got%% *}
      if [ "$got" = "$want" ]; then
        say "$rel: already in place; read back, sha256 $got equal to the record"
        printf 'ok\t%s\t%s\n' "$want" "$rel" >> "$out"
      else
        say "$rel: FAILED: already in place, and read back its sha256 ${got:-?} differs from the record $want (left as it is)"
        printf 'fail\t\t%s\n' "$rel" >> "$out"
      fi
    else
      say "$rel: already in place on the pod; skipped, not read"
      printf 'inplace\t\t%s\n' "$rel" >> "$out"
    fi
    return
  fi
  for try in $(seq 1 "$TRIES"); do
    if ! pod "mkdir -p -- '$dir' && rm -f -- '$part'"; then
      say "$rel: try $try: cannot prepare $dir on the pod"
    else
      say "$rel: try $try: copy started ($size bytes)"
      t0=$(date +%s)
      ssh -n -F "$CFG" tower "cat -- '$SRC/$rel'" | ssh -F "$CFG" pod "tee -- '$part' | sha256sum" > "$TMPD/hash.$i"
      st="${PIPESTATUS[0]} ${PIPESTATUS[1]}"
      t1=$(date +%s)
      got=$(cat "$TMPD/hash.$i"); got=${got%% *}
      if [ "$st" != "0 0" ]; then
        say "$rel: try $try: the stream failed (tower cat, pod tee: exit $st)"
      elif [ "$got" != "$want" ]; then
        say "$rel: try $try: sha256 of the stream on the pod $got differs from the record $want"
      else
        psize=$(pod "stat -c %s -- '$part'")
        if [ "$psize" != "$size" ]; then
          say "$rel: try $try: $part holds ${psize:-?} bytes, tower $size"
        else
          got2=$got
          if [ "$REVERIFY" = 1 ]; then got2=$(pod "sha256sum -- '$part'"); got2=${got2%% *}; fi
          if [ "$got2" != "$want" ]; then
            say "$rel: try $try: read back on the pod, sha256 $got2 differs from the record $want"
          elif pod "mv -n -- '$part' '$dest' && [ ! -e '$part' ] && [ -f '$dest' ]"; then
            say "$rel: OK: $size bytes in $((t1 - t0)) s ($(mbs "$size" $((t1 - t0))) MB/s), sha256 $want equal to the record$([ "$REVERIFY" = 1 ] && echo ', read back again')"
            printf 'ok\t%s\t%s\n' "$want" "$rel" >> "$out"
            return
          else
            pod "rm -f -- '$part'"
            say "$rel: FAILED: $dest appeared on the pod during the copy; not overwritten, the copy removed"
            printf 'fail\t\t%s\n' "$rel" >> "$out"
            return
          fi
        fi
      fi
      pod "rm -f -- '$part'" || say "$rel: try $try: could not remove $part"
    fi
    [ "$try" = "$TRIES" ] || sleep "$PAUSE"
  done
  say "$rel: FAILED after $TRIES tries"
  printf 'fail\t\t%s\n' "$rel" >> "$out"
}

# --- queues: largest first, each file to the stream with the fewest bytes -----------------------------------------
for i in $(seq 0 $((N - 1))); do
  [ "${HAS[$i]}" = yes ] && { printf '%s\t%s\n' 0 "$i"; continue; }
  [ "${SIZE[$i]}" = missing ] || printf '%s\t%s\n' "${SIZE[$i]}" "$i"
done | sort -t "$(printf '\t')" -k1,1nr -k2,2n | awk -F '\t' -v n="$STREAMS" -v d="$TMPD" '
  { k = 0; for (j = 1; j < n; j++) if (load[j] < load[k]) k = j; load[k] += $1; print $2 > (d "/queue." k) }'
for k in $(seq 0 $((STREAMS - 1))); do
  [ -f "$TMPD/queue.$k" ] || continue
  ( while IFS= read -r i; do
      if [ "${HAS[$i]}" != yes ] && stopped; then  # AGENTS §7.4: a STOP file ends the run before the next copy
        say "STOP file in $DST on the pod (or the check failed): this stream stops before ${REL[$i]}"
        printf 'stop\t\t%s\n' "${REL[$i]}" >> "$TMPD/result.$k"
        while IFS= read -r j; do printf 'stop\t\t%s\n' "${REL[$j]}" >> "$TMPD/result.$k"; done
        break
      fi
      one "$i" "$TMPD/result.$k"
    done < "$TMPD/queue.$k" ) &
done
wait

# --- the tally, and SHA256SUMS on the pod -------------------------------------------------------------------------
cat "$TMPD"/result.* > "$TMPD/results" 2>/dev/null || : > "$TMPD/results"
ok=$(awk -F '\t' '$1 == "ok"' "$TMPD/results" | wc -l | tr -d ' ')
failed=$(awk -F '\t' '$1 == "fail"' "$TMPD/results" | wc -l | tr -d ' ')
skipped=$(awk -F '\t' '$1 == "inplace"' "$TMPD/results" | wc -l | tr -d ' ')
halted=$(awk -F '\t' '$1 == "stop"' "$TMPD/results" | wc -l | tr -d ' ')
RC=0
[ "$failed" = 0 ] && [ "$missing" = 0 ] && [ "$halted" = 0 ] || RC=1

if [ "$ok" -gt 0 ] && [ "$ON_STOR" = 1 ]; then
  say "SHA256SUMS: not written: $DST is on /stor, where no file under 64 MiB is written (AGENTS §7.3); each OK line above holds the sha256"
elif [ "$ok" -gt 0 ]; then
  if pod "if [ -f '$DST/SHA256SUMS' ]; then cat -- '$DST/SHA256SUMS'; fi" > "$TMPD/sums.old"; then
    awk -F '\t' '$1 == "ok" { print $2 "  " $3 }' "$TMPD/results" > "$TMPD/sums.new"
    # one line per path; the old file's `*` marker is dropped; a path with two sha256 values is a conflict
    awk 'length($1) == 64 && length($0) > 66 { print substr($0, 67) "\t" $1 }' "$TMPD/sums.old" "$TMPD/sums.new" \
      | sort -u > "$TMPD/sums.pairs"
    conflicts=$(cut -f1 "$TMPD/sums.pairs" | uniq -d)
    if [ -n "$conflicts" ]; then
      say "SHA256SUMS: CONFLICT, left as it was: these paths have another sha256 there: $(echo $conflicts)"
      RC=1
    else
      awk -F '\t' '{ print $2 "  " $1 }' "$TMPD/sums.pairs" > "$TMPD/sums.merged"
      if cmp -s "$TMPD/sums.merged" "$TMPD/sums.old"; then
        say "SHA256SUMS: unchanged ($(wc -l < "$TMPD/sums.merged" | tr -d ' ') lines)"
      elif ssh -F "$CFG" pod "cat > '$DST/.SHA256SUMS.tmp' && mv -f -- '$DST/.SHA256SUMS.tmp' '$DST/SHA256SUMS'" < "$TMPD/sums.merged"; then
        say "SHA256SUMS: written, $(wc -l < "$TMPD/sums.merged" | tr -d ' ') lines ($ok verified in this run)"
      else
        say "SHA256SUMS: FAILED to write on the pod"
        RC=1
      fi
    fi
  else
    say "SHA256SUMS: FAILED to read on the pod"
    RC=1
  fi
fi
[ "$skipped" = 0 ] || say "note: $skipped file(s) in place were not read; they enter SHA256SUMS only when verified (--reverify)"
if [ "$RC" = 0 ]; then
  say "DONE: $ok verified, $skipped in place and skipped, 0 failed"
else
  say "ENDED WITH FAILURES OR A STOP: $ok verified, $skipped skipped, $failed failed, $missing missing on tower, $halted not copied after a STOP; read this log"
fi
tail -n 3 "$LOG"
exit $RC
