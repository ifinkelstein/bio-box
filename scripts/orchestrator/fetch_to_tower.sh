#!/bin/bash
# Fetch files by URL into the PI's NAS staging area, run on `tower` itself (card R-08). It is the script of the three
# fetch batches of 2026-11-25 (`_fetch_one.sh` and `_fetch_batch2.sh`, which existed only on tower), with one rule added:
# a URL that no entry of sources.yaml covers is refused before anything is fetched (AGENTS.md section 8).
#
#   fetch_to_tower.sh [options] JOBS_TSV [JOBS_TSV ...]
#
# JOBS_TSV  a job list, one file per line, five tab-separated fields:
#             <name> <TAB> <file> <TAB> <url> <TAB> <md5 or -> <TAB> <licence>
#           name     a directory under DEST, letters, digits . _ - in components joined by / (no . or .. component)
#           file     the name the file is kept under in DEST/<name>/: letters, digits, space and . _ , = + ( ) -;
#                    no /, no leading -, not ending .partial, not SOURCE_URLS.tsv or SHA256SUMS
#           url      https, or ftp where no https route exists (DEC-81 (5), card R-15); must match a URL form of an
#                    entry of sources.yaml (below), and an ftp:// URL matches only a form that is itself ftp://
#           md5      the provider's md5 (32 lower-case hex digits) or -; a file whose md5 differs is not kept
#           licence  as the job list states it, recorded unchanged (a literal {} stays {}); empty is recorded as -
#         Blank lines and lines starting # are skipped. The four-field lists of batch 1 (no licence) are accepted.
#         Every line of every list is checked, and every URL matched, before the first request: any refusal, exit 2.
#
# The URL forms are read from sources.yaml at run time: every top-level `url` or `url_<what>:` line of an entry whose value
# starts with https:// (or http://, or ftp://) is a form. {name} stands for one path segment (no /, ? or #); {path} and
# {x_path} for a path of one or more segments. Everything else is literal, the scheme included, so an https URL matches
# only an https:// form and an ftp URL only an ftp:// form; the URL must equal the form from its first to its last
# character, and the host part of a form may hold no variable. An http URL is refused whatever is listed. curl is held to
# the job's own protocol (--proto and --proto-redir =https, or =ftp for an ftp job); redirects are followed within it, but
# only the URL given is matched. The default User-Agent is `curl` (-A curl); no key, name, email address or project
# detail is sent (an ftp login is curl's own anonymous one).
#
# Per file (as run on 2026-11-25): a file already in DEST/<name>/ and in its SHA256SUMS is skipped, never fetched again;
# else curl writes DEST/<name>/.<file>.partial (a .partial is not a landed file); the partial must be non-empty, equal the
# md5 if one is given, and pass the archive test of its extension (gzip -t, unzip -tq, tar -tf, a FASTA's first byte >)
# where that tool exists (a file whose md5 equalled the provider's is kept when only that test fails, and is reported as
# OK plus a second stdout line WARN_<TEST>\t<name>\t<file>; a file with no md5 that fails it is removed, FAIL_<TEST>);
# then the sha256 and size are taken, `mv -n` puts it in place, and one line is appended to
# DEST/<name>/SHA256SUMS (`<sha256>  <file>`) and to DEST/<name>/SOURCE_URLS.tsv (filename, url, bytes, sha256, license,
# date_utc; the header is written when the file is made). A failed file is reported (FAIL_*, on stderr), its partial is
# removed, and the run goes on with the next. If DEST/STOP exists the run ends cleanly before the next file.
#
# Several connections (--jobs N, card R-22): the run, still the one holder of the lock, deals the jobs out in turn (the
# i-th job of the lists to worker i mod N) and starts N worker processes; each runs the steps above on its own share,
# one file at a time, and pauses --pause seconds after each file it lands. A destination is named once in the lists (a
# second naming is refused before any request, as with one connection), so no two workers ever fetch the same file.
# What workers share is written by one worker at a time, under a lock directory in the run's temporary directory: a
# directory's SHA256SUMS and SOURCE_URLS.tsv (made, read for the skip, appended), the move into place, and every line
# of output, so each line is whole; lines come in the order files end, not in list order. A failed file is one FAIL_
# line and its worker goes on. Each worker hands its counts to the run, which prints the one DONE line. A STOP file
# ends each worker before its next file. TERM or INT ends every worker (a transfer at once; a move into place and its
# two record lines are finished first), each removes its .partial, and then the lock is released. A worker that ends
# without finishing its share (killed) stops the run: the others are ended, exit 1, no DONE line. The lock names the
# workers and each open transfer (curl.<K>): a transfer left by a killed worker is ended by the run, and the lock is
# kept while it lives; a lock whose run has ended is not taken while a worker or a transfer of it lives; a worker whose
# run has ended stops before its next file, or while it waits for the write lock. Ctrl-C (INT to the process group) is
# acted on by the run alone. --jobs above 1 needs bash 4.3 or later (wait -n). With --jobs 1 nothing of this applies.
#
# Options:
#   --dest DIR     where files go (default /mnt/user/research-bulk/bio-box/sources_incoming)
#   --sources FILE the sources.yaml to match against (default: the one two directories above this script)
#   --check        read and match only: print `match <entry> <url>` per job and a total; no request, nothing written
#   --pause S      seconds between files, per worker (default 1)
#   --jobs N       files fetched at a time, each by its own worker over its own connection, 1 to 8 (default 1)
#
# Exit: 0 every file landed or skipped (or, with --check, every URL matched); 1 a file failed, or a worker ended without
# finishing its share; 2 a refusal before any request; 3 another run, or a worker of an ended run, holds
# DEST/.fetch_to_tower.lock.
set -u
export LC_ALL=C

usage() { sed -n '2,/^set -u/p' "$0" | sed '$d' | sed 's/^# \{0,1\}//'; }
die() { echo "fetch_to_tower: $*" >&2; exit 2; }

DEST=/mnt/user/research-bulk/bio-box/sources_incoming
HERE=${0%/*}
[ "$HERE" != "$0" ] || HERE=.
SOURCES=$HERE/../../sources.yaml
CHECK=0 PAUSE=1 JOBS=1
while [ $# -gt 0 ]; do
  case "$1" in
    --dest) [ $# -ge 2 ] || die "--dest needs a value"; DEST=$2; shift 2 ;;
    --sources) [ $# -ge 2 ] || die "--sources needs a value"; SOURCES=$2; shift 2 ;;
    --check) CHECK=1; shift ;;
    --pause) [ $# -ge 2 ] || die "--pause needs a value"; PAUSE=$2; shift 2 ;;
    --jobs) [ $# -ge 2 ] || die "--jobs needs a value"; JOBS=$2; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    --) shift; break ;;
    -*) die "unknown option $1 (see --help)" ;;
    *) break ;;
  esac
done
[ $# -ge 1 ] || die "usage: fetch_to_tower.sh [options] JOBS_TSV [JOBS_TSV ...] (see --help)"
case "$PAUSE" in ''|*[!0-9]*) die "--pause must be a whole number of seconds" ;; esac
case "$JOBS" in [1-8]) ;; *) die "--jobs must be a whole number from 1 to 8: $JOBS" ;; esac
if [ "$JOBS" -gt 1 ] && { [ "${BASH_VERSINFO[0]}" -lt 4 ] || { [ "${BASH_VERSINFO[0]}" = 4 ] && [ "${BASH_VERSINFO[1]}" -lt 3 ]; }; }; then
  die "--jobs above 1 needs bash 4.3 or later (wait -n); this is bash $BASH_VERSION: run with --jobs 1"
fi
DEST=${DEST%/}
case "$DEST" in /*) ;; *) die "--dest must be an absolute path" ;; esac
case "$DEST" in *[!A-Za-z0-9._/=+,-]*|*/./*|*/../*|*/-*) die "--dest holds a character or component that is refused: $DEST" ;; esac
[ -f "$SOURCES" ] || die "no sources file: $SOURCES"
for f in "$@"; do [ -f "$f" ] || die "no job list: $f"; done

TMPD=$(mktemp -d "${TMPDIR:-/tmp}/fetch_to_tower.XXXXXX") || die "cannot make a temporary directory"
LOCK=
PAR=0  # 1 while workers of --jobs above 1 may run: the run then ends them before it removes anything
KEEPLOCK=0  # 1 when a transfer of a killed worker still runs: the lock then stays, naming it
cleanup() { [ "$PAR" = 0 ] || stop_workers; rm -rf "$TMPD"; [ -z "$LOCK" ] || [ "$KEEPLOCK" = 1 ] || rm -rf "$LOCK"; }
trap cleanup EXIT
trap 'echo "fetch_to_tower: stopped by a signal" >&2; exit 130' INT TERM

# --- the URL forms of sources.yaml: one line `<entry><TAB><form>` -------------------------------------------------------
awk -v q="'" '
  /^  - name:/ { e = $0; sub(/^  - name:[ \t]*/, "", e); next }
  /^    url[A-Za-z0-9_]*:[ \t]+/ {
    v = $0; sub(/^    url[A-Za-z0-9_]*:[ \t]+/, "", v)
    sub("^[\"" q "]", "", v)
    split(v, a, /[ \t]/); u = a[1]; sub("[\"" q "]$", "", u)
    if (u ~ /^(https?|ftp):\/\//) print e "\t" u
  }' "$SOURCES" > "$TMPD/forms"
[ -s "$TMPD/forms" ] || die "$SOURCES holds no URL form"

# --- check every line of every list and match every URL, before any request -------------------------------------------
: > "$TMPD/jobs"
awk -F '\t' -v forms="$TMPD/forms" -v out="$TMPD/jobs" '
  function to_ere(f,   o, i, c, j, name, host) {
    host = f; sub(/^(https?|ftp):\/\//, "", host); sub(/\/.*/, "", host)
    if (host == "" || index(host, "{") || index(host, "}")) return ""
    o = ""; i = 1
    while (i <= length(f)) {
      c = substr(f, i, 1)
      if (c == "{") {
        j = index(substr(f, i), "}")
        if (j == 0) return ""
        name = substr(f, i + 1, j - 2)
        if (name !~ /^[A-Za-z_][A-Za-z0-9_]*$/) return ""
        o = o ((name == "path" || name ~ /_path$/) ? "[^?#]+" : "[^/?#]+")
        i += j
      } else {
        if (index("\\.?+*()[]^$|{}", c)) o = o "\\"
        o = o c; i++
      }
    }
    return "^" o "$"
  }
  function bad(msg) { printf "fetch_to_tower: %s line %d: %s\n", FILENAME, FNR, msg > "/dev/stderr"; nbad++ }
  BEGIN {
    while ((getline line < forms) > 0) {
      split(line, p, "\t"); r = to_ere(p[2])
      if (r != "") { nf++; ere[nf] = r; ent[nf] = p[1] }
    }
    if (nf == 0) { print "fetch_to_tower: no usable URL form in sources.yaml" > "/dev/stderr"; exit 2 }
  }
  /^[ \t]*$/ { next }
  /^#/ { next }
  {
    if (NF != 4 && NF != 5) { bad("expected 4 or 5 tab-separated fields, found " NF); next }
    name = $1; file = $2; url = $3; md5 = $4; lic = (NF == 5 ? $5 : "")
    ok = 1
    if (name !~ /^[A-Za-z0-9_][A-Za-z0-9_.-]*(\/[A-Za-z0-9_][A-Za-z0-9_.-]*)*$/ || name ~ /(^|\/)\.\.?(\/|$)/) { bad("refused name: " name); ok = 0 }
    if (file !~ /^[A-Za-z0-9_., ()=+-]+$/ || file ~ /^-/ || file == "." || file == ".." || file ~ /\.partial$/ || file == "SOURCE_URLS.tsv" || file == "SHA256SUMS") { bad("refused file name: " file); ok = 0 }
    if (url !~ /^(https?|ftp):\/\// || url !~ /^[A-Za-z0-9._~:\/?#@%&=+,;!*()-]+$/) { bad("refused URL: " url); ok = 0 }
    else if (url !~ /^(https|ftp):\/\//) { bad("not https: " url); ok = 0 }
    if (!(md5 == "-" || (length(md5) == 32 && md5 ~ /^[0-9a-f]+$/))) { bad("md5 is neither - nor 32 lower-case hex digits: " md5); ok = 0 }
    if (lic !~ /^[A-Za-z0-9._ {}()=+,\/-]*$/) { bad("refused licence text: " lic); ok = 0 }
    if (lic == "") lic = "-"
    key = name "/" file
    if (key in seen) { bad("the list names " key " twice (first at line " seen[key] ")"); ok = 0 } else seen[key] = FNR
    if (!ok) next
    hit = ""
    for (k = 1; k <= nf; k++) if (url ~ ere[k]) { hit = ent[k]; break }
    if (hit == "") { bad("URL not covered by any sources.yaml entry: " url); next }
    printf "%s\t%s\t%s\t%s\t%s\t%s\n", name, file, url, md5, lic, hit >> out
    njobs++
  }
  END { if (nbad > 0) { printf "fetch_to_tower: %d refusal(s); nothing was fetched\n", nbad > "/dev/stderr"; exit 2 } }
' "$@"
rc=$?
[ "$rc" = 0 ] || exit 2
N=$(wc -l < "$TMPD/jobs")
N=$((N + 0))
[ "$N" -gt 0 ] || die "the lists name no file"

if [ "$CHECK" = 1 ]; then
  awk -F '\t' '{ printf "match %s %s\n", $6, $3 }' "$TMPD/jobs"
  echo "check: $N URLs, all covered by sources.yaml; nothing fetched"
  exit 0
fi

# --- one run at a time on a destination ------------------------------------------------------------------------------
mkdir -p -- "$DEST" || die "cannot make $DEST"
LOCKDIR=$DEST/.fetch_to_tower.lock
if ! mkdir "$LOCKDIR" 2>/dev/null; then
  other=$(cat "$LOCKDIR/pid" 2>/dev/null)
  if [ -n "$other" ] && kill -0 "$other" 2>/dev/null; then
    echo "fetch_to_tower: another run (pid $other) holds $LOCKDIR" >&2
    exit 3
  fi
  for w in $(cat "$LOCKDIR/workers" 2>/dev/null); do
    if kill -0 "$w" 2>/dev/null; then
      echo "fetch_to_tower: a worker (pid $w) of an ended run still holds $LOCKDIR" >&2
      exit 3
    fi
  done
  for c in "$LOCKDIR"/curl.*; do
    [ -f "$c" ] || continue
    t=$(cat "$c" 2>/dev/null)
    if [ -n "$t" ] && kill -0 "$t" 2>/dev/null; then
      echo "fetch_to_tower: a transfer (pid $t) of an ended run still runs; $LOCKDIR is kept until it ends" >&2
      exit 3
    fi
  done
  rm -rf "$LOCKDIR"
  mkdir "$LOCKDIR" 2>/dev/null || { echo "fetch_to_tower: cannot take $LOCKDIR" >&2; exit 3; }
fi
LOCK=$LOCKDIR
echo $$ > "$LOCK/pid"

sha256_of() { if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | awk '{print $1}'; else shasum -a 256 "$1" | awk '{print $1}'; fi; }
md5_of() { if command -v md5sum >/dev/null 2>&1; then md5sum "$1" | awk '{print $1}'; else md5 -q "$1"; fi; }
in_sums() { [ -f "$1" ] && awk -v f="$2" '{ l = $0; sub(/^[0-9a-f]+  /, "", l); if (l == f) found = 1 } END { exit !found }' "$1"; }

OK=0 SKIP=0 FAIL=0
CPID= INBOOK=0 STOPPING=0 cur= FORKING=0 BEFORE= WK=
# hold and release take and give the run's write lock when --jobs is above 1; with one connection they do nothing. The
# lock is a symlink whose target is the holder's pid (ln -s makes it with its content at once, or fails), so a stop
# can tell its own lock from another's. While the run lives, a holder that dies is a worker that dies, and the run
# stops every worker; once the run has ended, a waiter ends when the holder is dead. A stop that reaches a holder is
# acted on at release.
hold() {
  [ "$JOBS" = 1 ] && return 0
  until ln -s "$BASHPID" "$TMPD/write.lock" 2>/dev/null; do
    if ! kill -0 $$ 2>/dev/null; then  # the run has ended: nobody stops a holder that died, so end if it has
      h=$(readlink "$TMPD/write.lock" 2>/dev/null)
      [ -z "$h" ] || kill -0 "$h" 2>/dev/null || worker_end
    fi
    sleep 0.1 2>/dev/null || sleep 1
  done
  INBOOK=1
}
release() { [ "$JOBS" = 1 ] && return 0; rm -f "$TMPD/write.lock"; INBOOK=0; [ "$STOPPING" = 0 ] || worker_end; }
fail() { hold; printf '%s\t%s\t%s\n' "$1" "$2" "$3" >&2; release; FAIL=$((FAIL + 1)); }
get() {  # curl: in the foreground with --jobs 1; else in the background, so that a stop ends the transfer at once
  if [ "$JOBS" = 1 ]; then curl "$@"; return; fi
  # the transfer writes its own pid into the lock before curl starts (exec keeps the pid), so that a transfer whose
  # worker is killed is still named there: the run ends it, and no run takes the lock while it lives
  BEFORE=${!:-} FORKING=1
  ( echo "$BASHPID" > "$LOCK/curl.$WK"; exec curl "$@" ) &
  CPID=$!
  FORKING=0
  wait "$CPID"
  local rc=$?
  rm -f -- "$LOCK/curl.$WK"
  CPID=
  return "$rc"
}

fetch_jobs() {  # fetch_jobs < JOBS: the steps of each file, one file at a time (the loop of --jobs 1, and of each worker)
while IFS=$'\t' read -r name file url md5 lic entry; do
  if [ -e "$DEST/STOP" ]; then hold; echo "STOP file found in $DEST: ended cleanly before $name/$file"; release; break; fi
  if [ "$JOBS" != 1 ] && ! kill -0 $$ 2>/dev/null; then
    hold; echo "fetch_to_tower: the run (pid $$) has ended: a worker stops before $name/$file" >&2; release; break
  fi
  dir=$DEST/$name
  mkdir -p -- "$dir"
  cd "$dir" || { fail FAIL_DIR "$name" "$file"; continue; }
  hold
  [ -f SOURCE_URLS.tsv ] || printf 'filename\turl\tbytes\tsha256\tlicense\tdate_utc\n' > SOURCE_URLS.tsv
  [ -f SHA256SUMS ] || : > SHA256SUMS
  if [ -s "$file" ] && in_sums SHA256SUMS "$file"; then
    echo "SKIP	$name	$file"
    release
    SKIP=$((SKIP + 1))
    continue
  fi
  release
  partial=.$file.partial
  cur=$dir/$partial
  rm -f -- "$partial"
  proto=https
  case "$url" in ftp://*) proto=ftp ;; esac
  if ! get -fsSL -A curl --proto "=$proto" --proto-redir "=$proto" --retry 3 --retry-delay 5 --connect-timeout 30 -o "$partial" "$url"; then
    rm -f -- "$partial"; fail FAIL_CURL "$name/$file" "$url"; continue
  fi
  if [ ! -s "$partial" ]; then rm -f -- "$partial"; fail FAIL_EMPTY "$name" "$file"; continue; fi
  if [ "$md5" != "-" ]; then
    got=$(md5_of "$partial")
    if [ "$got" != "$md5" ]; then rm -f -- "$partial"; fail FAIL_MD5 "$name/$file" "expected=$md5 got=$got"; continue; fi
  fi
  bad=
  case "$file" in
    *.zip) command -v unzip >/dev/null 2>&1 && ! unzip -tq "$partial" >/dev/null 2>&1 && bad=FAIL_ZIP ;;
    *.tar.gz|*.tgz|*.gz) command -v gzip >/dev/null 2>&1 && ! gzip -t "$partial" 2>/dev/null && bad=FAIL_GZIP ;;
    *.tar) command -v tar >/dev/null 2>&1 && ! tar -tf "$partial" >/dev/null 2>&1 && bad=FAIL_TAR ;;
    *.fasta|*.fna|*.fa) [ "$(head -c 1 "$partial")" = ">" ] || bad=FAIL_FASTA ;;
  esac
  warn=
  if [ -n "$bad" ] && [ "$md5" != "-" ]; then warn=WARN_${bad#FAIL_}; bad=; fi
  if [ -n "$bad" ]; then rm -f -- "$partial"; fail "$bad" "$name" "$file"; continue; fi
  bytes=$(wc -c < "$partial")
  bytes=$((bytes + 0))
  sha=$(sha256_of "$partial")
  dt=$(date -u +%Y-%m-%dT%H:%M:%SZ)
  hold
  mv -n -- "$partial" "$file"
  if [ ! -f "$file" ] || [ -f "$partial" ]; then rm -f -- "$partial"; release; fail FAIL_MV "$name" "$file"; continue; fi
  printf '%s  %s\n' "$sha" "$file" >> SHA256SUMS
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$file" "$url" "$bytes" "$sha" "$lic" "$dt" >> SOURCE_URLS.tsv
  printf 'OK\t%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$file" "$bytes" "$sha" "$lic" "$dt"
  [ -z "$warn" ] || printf '%s\t%s\t%s\n' "$warn" "$name" "$file"
  release
  OK=$((OK + 1))
  [ "$PAUSE" = 0 ] || sleep "$PAUSE"
done
}

[ "$JOBS" -le "$N" ] || JOBS=$N
if [ "$JOBS" = 1 ]; then
  fetch_jobs < "$TMPD/jobs"
else
  # --- --jobs above 1: the jobs dealt out in turn to JOBS workers, each the loop above on its own share ---------------
  worker_end() {  # a stop reached a worker: end its transfer, give back the write lock if it holds it, remove its .partial
    # a stop between the fork of a transfer and CPID=$!: $! is then already the new transfer
    if [ -z "$CPID" ] && [ "$FORKING" = 1 ] && [ "${!:-}" != "$BEFORE" ]; then CPID=$!; fi
    [ -z "$CPID" ] || { kill -TERM "$CPID" 2>/dev/null; wait "$CPID" 2>/dev/null; rm -f -- "$LOCK/curl.$WK"; }
    [ "$(readlink "$TMPD/write.lock" 2>/dev/null)" != "$BASHPID" ] || rm -f "$TMPD/write.lock"
    [ -z "$cur" ] || rm -f -- "$cur"
    exit 143
  }
  worker() {  # worker K: the jobs of share K; its counts to tally.K when its share is done
    trap - EXIT
    WK=$1
    trap '' INT  # a Ctrl-C reaches the whole process group: only the run acts on it, and it ends the workers by TERM
    trap 'STOPPING=1; [ "$INBOOK" = 1 ] || worker_end' TERM
    fetch_jobs < "$TMPD/jobs.$1"
    printf '%s %s %s\n' "$OK" "$SKIP" "$FAIL" > "$TMPD/tally.$1"
    exit 0
  }
  stop_workers() {  # end every worker that has not handed in its counts (and is not known gone), then wait for all
    trap '' INT TERM
    k=0
    while [ "$k" -lt "$JOBS" ]; do
      p=${WPID[$k]:-}
      [ -z "$p" ] || [ -f "$TMPD/tally.$k" ] || [ -f "$TMPD/gone.$k" ] || kill -TERM "$p" 2>/dev/null
      k=$((k + 1))
    done
    [ "$LAUNCHED" -ge "$JOBS" ] || [ -z "${!:-}" ] || kill -TERM "$!" 2>/dev/null  # a stop while starting them
    wait
    # A transfer still named in the lock belongs to a worker that was killed: it is no child of the run and cannot be
    # waited for. End it, give it 2 s, and keep the lock (naming it) while it lives, so no run starts beside it.
    t=0
    while :; do
      live=
      for c in "$LOCK"/curl.*; do
        [ -f "$c" ] || continue
        p=$(cat "$c" 2>/dev/null)
        if [ -n "$p" ] && kill -0 "$p" 2>/dev/null; then
          live="$live $p"
          [ "$t" != 0 ] || kill -TERM "$p" 2>/dev/null
        else
          rm -f -- "$c"
        fi
      done
      [ -n "$live" ] || break
      if [ "$t" -ge 20 ]; then
        KEEPLOCK=1
        echo "fetch_to_tower: transfer(s)$live of a killed worker still run; $LOCK is kept until they end" >&2
        break
      fi
      t=$((t + 1))
      sleep 0.1 2>/dev/null || sleep 1
    done
  }
  awk -v n="$JOBS" -v d="$TMPD" '{ print > (d "/jobs." ((NR - 1) % n)) }' "$TMPD/jobs"
  LAUNCHED=0 PAR=1
  while [ "$LAUNCHED" -lt "$JOBS" ]; do
    worker "$LAUNCHED" &
    WPID[$LAUNCHED]=$!
    echo "$!" >> "$LOCK/workers"  # at once: a run killed while it starts them still names every worker
    LAUNCHED=$((LAUNCHED + 1))
  done
  while :; do
    wait -n
    wrc=$?  # 127: no worker left to wait for
    left=0 lost=
    k=0
    while [ "$k" -lt "$JOBS" ]; do
      if [ ! -f "$TMPD/tally.$k" ]; then
        if ! kill -0 "${WPID[$k]}" 2>/dev/null; then lost="$lost $k"; : > "$TMPD/gone.$k"
        elif [ "$wrc" = 127 ]; then lost="$lost $k"
        else left=$((left + 1)); fi
      fi
      k=$((k + 1))
    done
    if [ -n "$lost" ]; then
      echo "fetch_to_tower: worker(s)$lost ended without finishing a share; the others are stopped; a rerun fetches what is missing" >&2
      exit 1
    fi
    [ "$left" -gt 0 ] || break
  done
  wait  # every worker has handed in its counts: let each one exit before they are read
  PAR=0
  k=0
  while [ "$k" -lt "$JOBS" ]; do
    read -r o s f < "$TMPD/tally.$k"
    OK=$((OK + o)) SKIP=$((SKIP + s)) FAIL=$((FAIL + f))
    k=$((k + 1))
  done
fi

echo "DONE	ok=$OK	skip=$SKIP	fail=$FAIL"
[ "$FAIL" = 0 ] || exit 1
exit 0
