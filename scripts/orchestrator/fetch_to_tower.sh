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
#           url      https only; must match a URL form of an entry of sources.yaml (below)
#           md5      the provider's md5 (32 lower-case hex digits) or -; a file whose md5 differs is not kept
#           licence  as the job list states it, recorded unchanged (a literal {} stays {}); empty is recorded as -
#         Blank lines and lines starting # are skipped. The four-field lists of batch 1 (no licence) are accepted.
#         Every line of every list is checked, and every URL matched, before the first request: any refusal, exit 2.
#
# The URL forms are read from sources.yaml at run time: every top-level `url` or `url_<what>:` line of an entry whose value
# starts with https:// (or http://) is a form. {name} stands for one path segment (no /, ? or #); {path} and {x_path} for
# a path of one or more segments. Everything else is literal, the URL must equal the form from its first to its last
# character, and the host part of a form may hold no variable. Redirects are followed (https only) but only the URL given
# is matched. The default User-Agent is `curl` (-A curl); no key, name, email address or project detail is sent.
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
# Options:
#   --dest DIR     where files go (default /mnt/user/research-bulk/bio-box/sources_incoming)
#   --sources FILE the sources.yaml to match against (default: the one two directories above this script)
#   --check        read and match only: print `match <entry> <url>` per job and a total; no request, nothing written
#   --pause S      seconds between files (default 1)
#
# Exit: 0 every file landed or skipped (or, with --check, every URL matched); 1 a file failed; 2 a refusal before any
# request; 3 another run holds DEST/.fetch_to_tower.lock.
set -u
export LC_ALL=C

usage() { sed -n '2,/^set -u/p' "$0" | sed '$d' | sed 's/^# \{0,1\}//'; }
die() { echo "fetch_to_tower: $*" >&2; exit 2; }

DEST=/mnt/user/research-bulk/bio-box/sources_incoming
HERE=${0%/*}
[ "$HERE" != "$0" ] || HERE=.
SOURCES=$HERE/../../sources.yaml
CHECK=0 PAUSE=1
while [ $# -gt 0 ]; do
  case "$1" in
    --dest) [ $# -ge 2 ] || die "--dest needs a value"; DEST=$2; shift 2 ;;
    --sources) [ $# -ge 2 ] || die "--sources needs a value"; SOURCES=$2; shift 2 ;;
    --check) CHECK=1; shift ;;
    --pause) [ $# -ge 2 ] || die "--pause needs a value"; PAUSE=$2; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    --) shift; break ;;
    -*) die "unknown option $1 (see --help)" ;;
    *) break ;;
  esac
done
[ $# -ge 1 ] || die "usage: fetch_to_tower.sh [options] JOBS_TSV [JOBS_TSV ...] (see --help)"
case "$PAUSE" in ''|*[!0-9]*) die "--pause must be a whole number of seconds" ;; esac
DEST=${DEST%/}
case "$DEST" in /*) ;; *) die "--dest must be an absolute path" ;; esac
case "$DEST" in *[!A-Za-z0-9._/=+,-]*|*/./*|*/../*|*/-*) die "--dest holds a character or component that is refused: $DEST" ;; esac
[ -f "$SOURCES" ] || die "no sources file: $SOURCES"
for f in "$@"; do [ -f "$f" ] || die "no job list: $f"; done

TMPD=$(mktemp -d "${TMPDIR:-/tmp}/fetch_to_tower.XXXXXX") || die "cannot make a temporary directory"
LOCK=
cleanup() { rm -rf "$TMPD"; [ -z "$LOCK" ] || rm -rf "$LOCK"; }
trap cleanup EXIT
trap 'echo "fetch_to_tower: stopped by a signal" >&2; exit 130' INT TERM

# --- the URL forms of sources.yaml: one line `<entry><TAB><form>` -------------------------------------------------------
awk -v q="'" '
  /^  - name:/ { e = $0; sub(/^  - name:[ \t]*/, "", e); next }
  /^    url[A-Za-z0-9_]*:[ \t]+/ {
    v = $0; sub(/^    url[A-Za-z0-9_]*:[ \t]+/, "", v)
    sub("^[\"" q "]", "", v)
    split(v, a, /[ \t]/); u = a[1]; sub("[\"" q "]$", "", u)
    if (u ~ /^https?:\/\//) print e "\t" u
  }' "$SOURCES" > "$TMPD/forms"
[ -s "$TMPD/forms" ] || die "$SOURCES holds no URL form"

# --- check every line of every list and match every URL, before any request -------------------------------------------
: > "$TMPD/jobs"
awk -F '\t' -v forms="$TMPD/forms" -v out="$TMPD/jobs" '
  function to_ere(f,   o, i, c, j, name, host) {
    host = f; sub(/^https?:\/\//, "", host); sub(/\/.*/, "", host)
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
    if (url !~ /^https?:\/\// || url !~ /^[A-Za-z0-9._~:\/?#@%&=+,;!*()-]+$/) { bad("refused URL: " url); ok = 0 }
    else if (url !~ /^https:\/\//) { bad("not https: " url); ok = 0 }
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
  rm -rf "$LOCKDIR"
  mkdir "$LOCKDIR" 2>/dev/null || { echo "fetch_to_tower: cannot take $LOCKDIR" >&2; exit 3; }
fi
LOCK=$LOCKDIR
echo $$ > "$LOCK/pid"

sha256_of() { if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | awk '{print $1}'; else shasum -a 256 "$1" | awk '{print $1}'; fi; }
md5_of() { if command -v md5sum >/dev/null 2>&1; then md5sum "$1" | awk '{print $1}'; else md5 -q "$1"; fi; }
in_sums() { [ -f "$1" ] && awk -v f="$2" '{ l = $0; sub(/^[0-9a-f]+  /, "", l); if (l == f) found = 1 } END { exit !found }' "$1"; }

OK=0 SKIP=0 FAIL=0
fail() { printf '%s\t%s\t%s\n' "$1" "$2" "$3" >&2; FAIL=$((FAIL + 1)); }

while IFS=$'\t' read -r name file url md5 lic entry; do
  if [ -e "$DEST/STOP" ]; then echo "STOP file found in $DEST: ended cleanly before $name/$file"; break; fi
  dir=$DEST/$name
  mkdir -p -- "$dir"
  cd "$dir" || { fail FAIL_DIR "$name" "$file"; continue; }
  [ -f SOURCE_URLS.tsv ] || printf 'filename\turl\tbytes\tsha256\tlicense\tdate_utc\n' > SOURCE_URLS.tsv
  [ -f SHA256SUMS ] || : > SHA256SUMS
  if [ -s "$file" ] && in_sums SHA256SUMS "$file"; then
    echo "SKIP	$name	$file"
    SKIP=$((SKIP + 1))
    continue
  fi
  partial=.$file.partial
  rm -f -- "$partial"
  if ! curl -fsSL -A curl --proto =https --proto-redir =https --retry 3 --retry-delay 5 --connect-timeout 30 -o "$partial" "$url"; then
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
  mv -n -- "$partial" "$file"
  if [ ! -f "$file" ] || [ -f "$partial" ]; then rm -f -- "$partial"; fail FAIL_MV "$name" "$file"; continue; fi
  printf '%s  %s\n' "$sha" "$file" >> SHA256SUMS
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$file" "$url" "$bytes" "$sha" "$lic" "$dt" >> SOURCE_URLS.tsv
  printf 'OK\t%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$file" "$bytes" "$sha" "$lic" "$dt"
  [ -z "$warn" ] || printf '%s\t%s\t%s\n' "$warn" "$name" "$file"
  OK=$((OK + 1))
  [ "$PAUSE" = 0 ] || sleep "$PAUSE"
done < "$TMPD/jobs"

echo "DONE	ok=$OK	skip=$SKIP	fail=$FAIL"
[ "$FAIL" = 0 ] || exit 1
exit 0
