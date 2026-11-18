#!/bin/sh
# Card C-40: fetch the two plasmid sources written into sources.yaml under DEC-50 (gut_plasmids_yu2024, plasmidscope)
# to /tmp/if2453/s3/plasmids/dl/<source>/, one file at a time over one connection. Run ON THE POD in a capped unit:
#
#     systemd-run --user --unit=c40-fetch --collect -p MemoryMax=4G -p MemorySwapMax=0 -p CPUQuota=100% \
#       nice -n 19 ionice -c 2 -n 7 sh /tmp/if2453/s3/plasmids/bin/c40_fetch.sh
#
# Only the URL forms of those two entries, exactly as written there. No key, email or identifying user agent
# (`-A curl`). A file already present is not fetched again. Each file: at most 3 tries (resuming with -C -); HTTP 403,
# 429 or 503 counts as a refusal and 5 refusals in a row stop the run; 404, 410 or 451 stop it at once. A STOP file in
# the data root or in /tmp/if2453/s3/plasmids stops it before each file. A finished file is renamed from <name>.part
# and checked: its size against the entry's bytes, its md5 against the provider's where one is published; a mismatch
# stops the run (the file is kept as <name>.bad). Size, md5, sha256 and the UTC time go to dl/SUMS. At the end
# dl/DONE, or dl/FAILED with the reason. Nothing fetched is executed: the archives are unpacked by c40_unpack.sh into
# new, empty directories. Log: logs/fetch.log.
set -u
W=/tmp/if2453/s3/plasmids
D=$W/dl
mkdir -p "$D/gut" "$D/plasmidscope" "$W/logs"
exec >> "$W/logs/fetch.log" 2>&1 < /dev/null
rm -f "$D/DONE" "$D/FAILED"
REFUSED=0

now() { date -u +%FT%TZ; }
fail() { echo "FAILED: $1 $(now)"; echo "$1" > "$D/FAILED"; exit 1; }

fetch() {  # name url dest bytes md5|-
  name=$1; url=$2; dest=$D/$3; bytes=$4; md5=$5
  for f in /tmp/if2453/r10/root/STOP "$W/STOP"; do
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
  got=$(stat -c %s "$dest.part")
  m=$(md5sum "$dest.part" | cut -d' ' -f1)
  if [ "$got" != "$bytes" ]; then mv "$dest.part" "$dest.bad"; fail "$name: $got bytes, the entry says $bytes"; fi
  if [ "$md5" != "-" ] && [ "$m" != "$md5" ]; then mv "$dest.part" "$dest.bad"; fail "$name: md5 $m, the provider's $md5"; fi
  mv "$dest.part" "$dest"
  echo "$3 $got bytes md5 $m sha256 $(sha256sum "$dest" | cut -d' ' -f1) provider_md5 $md5 $(now)" | tee -a "$D/SUMS"
  sleep 2
}

Z=https://zenodo.org/api/records/5843600/files
P=https://plasmidapi.deepomics.org/api/database/files
echo "c40_fetch start $(now) host $(hostname)"
fetch gut_summary "$Z/predicted_plasmids_summary.txt.gz/content" gut/predicted_plasmids_summary.txt.gz 14952780 e12cd9fee4ad63015a50f467bbd793c8
fetch gut_clusters "$Z/predicted_plasmids_clusters.txt.gz/content" gut/predicted_plasmids_clusters.txt.gz 1873405 ce6b77fd22c83ee9f2a56cf6370e8542
fetch gut_systems "$Z/predicted_plasmids_systems.txt.gz/content" gut/predicted_plasmids_systems.txt.gz 293764 21bf763a32925bc6eeb0dfdead45a389
fetch gut_selfrec "$Z/predicted_plasmids_assembly_self_recruitment.txt.gz/content" gut/predicted_plasmids_assembly_self_recruitment.txt.gz 2364939 6e3dc23190e70c8c2030d9971f5fa379
fetch gut_fasta "$Z/predicted_plasmids.fa.gz/content" gut/predicted_plasmids.fa.gz 252900436 3252f40283f594e1f703199a7f7e7dd8
while read -r src fb mb; do
  fetch "ps_${src}_meta" "$P/$src/data/$src.plasmid_list.download.tsv/" "plasmidscope/$src.plasmid_list.download.tsv" "$mb" -
  [ "$fb" = "-" ] || fetch "ps_${src}_fasta" "$P/$src/$src.fasta.tar.gz/" "plasmidscope/$src.fasta.tar.gz" "$fb" -
done <<EOF
ALL - 105861535
TPA 215293 906
Kraken2 15183583 106106
DDBJ 151612202 585010
ENA 181390959 675330
COMPASS 388863677 1435170
mMGE 313613870 9424875
PLSDB 1619480622 5797812
RefSeq 2148947321 10012206
GenBank 2328467541 10578194
IMG-PR 5235080408 81722123
EOF
now > "$D/DONE"
echo "c40_fetch done $(now)"
