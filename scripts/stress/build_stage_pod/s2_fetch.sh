#!/bin/sh
# Fetch the two S2 sources the PI approved on 2026-11-11 (sources.yaml: hpgc_huge_phages, plsdb) to
# /tmp/if2453/s2/dl, one file at a time over one connection. Run ON THE POD inside a capped unit:
#
#     systemd-run --user --unit=s2-fetch --collect -p MemoryMax=8G -p MemorySwapMax=0 -p CPUQuota=200% \
#       nice -n 19 ionice -c 2 -n 7 sh /tmp/if2453/s2/bin/s2_fetch.sh
#
# Only the three URLs of those two entries, exactly as written there. No key, email or identifying
# user agent (`-A curl`). Each try resumes the partial file (`-C -`); a server that does not resume
# makes the next try start again from byte 0. Up to 6 tries with backoff 30, 60, 120, 240, 480 s.
# HTTP 403, 404, 410 or 451 stops the run at once (a refusal). A STOP file in the data root or in
# /tmp/if2453/s2 stops it before each file. A finished file is renamed from <name>.part; its size,
# md5 and sha256 go to dl/SUMS. At the end dl/DONE, or dl/FAILED with the reason. Log: logs/fetch.log.
set -u
W=/tmp/if2453/s2
D=$W/dl
exec >> $W/logs/fetch.log 2>&1 < /dev/null
mkdir -p $D
rm -f $D/DONE $D/FAILED

now() { date -u +%FT%TZ; }
fail() { echo "FAILED: $1 $(now)"; echo "$1" > $D/FAILED; exit 1; }

fetch() {  # name url dest
  name=$1; url=$2; dest=$D/$3
  for f in /tmp/if2453/r10/root/STOP $W/STOP; do
    [ -e "$f" ] && fail "STOP present: $f"
  done
  if [ -e "$dest" ]; then echo "$name: $dest present, not fetched again"; return 0; fi
  n=0
  while :; do
    n=$((n + 1))
    echo "$name try $n start $(now) $url"
    curl -A curl --fail --location --silent --show-error --http1.1 --connect-timeout 60 \
      --speed-limit 1024 --speed-time 300 -D "$dest.headers.$n" -C - -o "$dest.part" "$url"
    rc=$?
    echo "$name try $n rc=$rc bytes=$(stat -c %s "$dest.part" 2>/dev/null || echo 0) $(now)"
    [ $rc = 0 ] && break
    st=$(grep -E '^HTTP/' "$dest.headers.$n" 2>/dev/null | tail -1 | awk '{print $2}')
    echo "$name try $n last HTTP status: ${st:-none}"
    case "$st" in 403|404|410|451) fail "$name refused: HTTP $st";; esac
    if [ $rc = 33 ] || [ $rc = 36 ]; then echo "$name: no resume, next try from byte 0"; rm -f "$dest.part"; fi
    [ $n -ge 6 ] && fail "$name: gave up after $n tries (last rc=$rc)"
    sleep $((30 * (1 << (n - 1))))
  done
  mv "$dest.part" "$dest"
  echo "$name $(stat -c %s "$dest") bytes md5 $(md5sum "$dest" | cut -d' ' -f1) sha256 $(sha256sum "$dest" | cut -d' ' -f1) $(now)" | tee -a $D/SUMS
}

echo "s2_fetch start $(now) host $(hostname)"
fetch hpgc https://ndownloader.figshare.com/files/59221811 HPGC.7230.genomes.final.fasta.gz
fetch plsdb_fasta https://ccb-microbe.cs.uni-saarland.de/plsdb2025/download_fasta plsdb_download_fasta
fetch plsdb_meta https://ccb-microbe.cs.uni-saarland.de/plsdb2025/download_meta.tar.gz download_meta.tar.gz
now > $D/DONE
echo "s2_fetch done $(now)"
