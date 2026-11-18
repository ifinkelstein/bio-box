#!/bin/sh
# Card C-41: fetch the GlobDB release 226 metadata tables and checksum file (entry globdb_r226 in sources.yaml,
# written under DEC-50 before any fetch) into /tmp/if2453/s3/globdb/dl/, one file at a time over one connection.
#
# Run ON THE POD in a capped unit (at most 2 threads, nice 19):
#
#     systemd-run --user --unit=c41-fetch-tables --collect -p MemoryMax=4G -p MemorySwapMax=0 -p CPUQuota=100% \
#       nice -n 19 ionice -c 2 -n 7 sh /tmp/if2453/ci/C-41/scripts/stress/build_stage_pod/c41_fetch.sh tables
#
# Modes:
#   tables: README.txt, globdb_r226_md5sums.txt, globdb_r226_tax_plus_stats.tsv.gz, globdb_r226_taxonomy.tsv.gz,
#           globdb_r226_dictionaries.tar.gz, and environmental_metadata/high_level_env__{pct,value}.tsv (54.4 MB).
#           Every top-level file listed in globdb_r226_md5sums.txt is checked against that md5, and every file is
#           checked against the Content-Length byte count recorded in sources.yaml.
#
# No key, email or identifying user agent (`-A curl`). Stop after 5 refusals in a row (HTTP 403, 429, 503) or at once
# on 404, 410, 451. Checks STOP in /tmp/if2453/r10/root and /tmp/if2453/s3/globdb before each file.
set -u
W=/tmp/if2453/s3/globdb
D=$W/dl
mkdir -p "$D/environmental_metadata" "$W/logs"
exec >> "$W/logs/fetch.log" 2>&1 < /dev/null
rm -f "$D/DONE" "$D/FAILED"
REFUSED=0
BASE=https://fileshare.lisc.univie.ac.at/globdb/globdb_r226

now() { date -u +%FT%TZ; }
fail() { echo "FAILED: $1 $(now)"; echo "$1" > "$D/FAILED"; exit 1; }

provider_md5() {  # file -> md5 from globdb_r226_md5sums.txt, or -
  f=$1
  if [ -f "$D/globdb_r226_md5sums.txt" ]; then
    m=$(awk -v t="$f" '$2 == t || $2 == "./" t {print $1}' "$D/globdb_r226_md5sums.txt")
    if [ -n "$m" ]; then echo "$m"; return 0; fi
  fi
  echo "-"
}

fetch() {  # name relpath bytes
  name=$1; rel=$2; bytes=$3
  url="$BASE/$rel"
  dest="$D/$rel"
  for f in /tmp/if2453/r10/root/STOP "$W/STOP"; do
    [ -e "$f" ] && fail "STOP present: $f"
  done
  if [ -e "$dest" ]; then echo "$name: $dest present, not fetched again"; return 0; fi
  n=0
  t0=$(date +%s)
  while :; do
    n=$((n + 1))
    echo "$name try $n start $(now) $url"
    curl -A curl --fail --location --silent --show-error --http1.1 --connect-timeout 60 \
      --speed-limit 1024 --speed-time 300 -D "$dest.headers.$n" -C - -o "$dest.part" "$url"
    rc=$?
    echo "$name try $n rc=$rc bytes=$(stat -c %s "$dest.part" 2>/dev/null || echo 0) $(now)"
    if [ $rc = 0 ]; then REFUSED=0; break; fi
    st=$(grep -E '^HTTP/' "$dest.headers.$n" 2>/dev/null | tail -1 | awk '{print $2}')
    echo "$name try $n last HTTP status: ${st:-none}"
    case "$st" in
      404|410|451) fail "$name refused: HTTP $st" ;;
      403|429|503) REFUSED=$((REFUSED + 1)); [ $REFUSED -ge 5 ] && fail "5 refusals in a row (last $name HTTP $st)" ;;
    esac
    if [ $rc = 33 ] || [ $rc = 36 ]; then echo "$name: no resume, next try from byte 0"; rm -f "$dest.part"; fi
    [ $n -ge 3 ] && fail "$name: gave up after $n tries (last rc=$rc)"
    sleep $((60 * n))
  done
  t1=$(date +%s)
  got=$(stat -c %s "$dest.part")
  m=$(md5sum "$dest.part" | cut -d' ' -f1)
  want=$(provider_md5 "$rel")
  if [ "$got" != "$bytes" ]; then mv "$dest.part" "$dest.bad"; fail "$name: $got bytes, sources.yaml says $bytes"; fi
  if [ "$want" != "-" ] && [ "$m" != "$want" ]; then
    mv "$dest.part" "$dest.bad"
    fail "$name: md5 $m, globdb_r226_md5sums.txt says $want"
  fi
  mv "$dest.part" "$dest"
  sha=$(sha256sum "$dest" | cut -d' ' -f1)
  echo "$rel $got bytes md5 $m sha256 $sha provider_md5 $want wall_s $((t1 - t0)) $(now)" | tee -a "$D/SUMS"
  sleep 1
}

tables() {
  fetch md5sums globdb_r226_md5sums.txt 526
  fetch readme README.txt 10941
  fetch tax_plus_stats globdb_r226_tax_plus_stats.tsv.gz 11811303
  fetch taxonomy globdb_r226_taxonomy.tsv.gz 5832515
  fetch dictionaries globdb_r226_dictionaries.tar.gz 3176715
  fetch env_high_pct environmental_metadata/high_level_env__pct.tsv 20736725
  fetch env_high_val environmental_metadata/high_level_env__value.tsv 12894501
}

echo "c41_fetch start $(now) host $(hostname) args: $*"
case "${1:-}" in
  tables) tables ;;
  *) fail "usage: c41_fetch.sh tables" ;;
esac
now > "$D/DONE"
echo "c41_fetch done $(now)"
