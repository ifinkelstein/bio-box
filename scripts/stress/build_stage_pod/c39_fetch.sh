#!/bin/sh
# Card C-39: fetch the sources.yaml entries uhgv, gsv_atlas and vire_tables to /tmp/if2453/s3/dl/<source>/, one file
# at a time over one connection. Run ON THE POD in a capped unit:
#
#     systemd-run --user --unit=c39-fetch --collect -p MemoryMax=4G -p MemorySwapMax=0 -p CPUQuota=100% \
#       nice -n 19 ionice -c 2 -n 7 sh /tmp/if2453/s3/bin/c39_fetch.sh uhgv vire
#
# The arguments name the sources to fetch, in order (uhgv, gsv, vire).
#
# Only the URLs of those three entries, exactly as written there. No key, email or identifying user agent (`-A curl`).
# Each try resumes the partial file (`-C -`) unless the server does not resume. At most 5 tries a file; HTTP 403,
# 404, 410 or 451 stops the run at once, and 5 refusals in a row (any HTTP status of 400 or more) stop it. A STOP
# file in the data root or in /tmp/if2453/s3 stops it before each file. A finished file is renamed from <name>.part;
# its size, md5 and sha256 go to dl/SUMS, and for a file with a provider md5 the two are compared (a mismatch
# stops the run). The GSV archive's size is read from a HEAD request first (the entry asks for it). At the end
# dl/DONE, or dl/FAILED with the reason. Log: logs/fetch.log. Nothing fetched is executed or unpacked here.
set -u
W=/tmp/if2453/s3
D=$W/dl
exec >> $W/logs/fetch.log 2>&1 < /dev/null
mkdir -p $D/uhgv $D/gsv $D/vire
rm -f $D/DONE $D/FAILED
REFUSALS=0

now() { date -u +%FT%TZ; }
fail() { echo "FAILED: $1 $(now)"; echo "$1" > $D/FAILED; exit 1; }

fetch() {  # name url dest [md5]
  name=$1; url=$2; dest=$D/$3; want=${4:-}
  for f in /tmp/if2453/r10/root/STOP $W/STOP; do
    [ -e "$f" ] && fail "STOP present: $f"
  done
  if [ -e "$dest" ]; then echo "$name: $dest present, not fetched again"; return 0; fi
  n=0
  while :; do
    n=$((n + 1))
    echo "$name try $n start $(now) $url"
    curl -A curl --fail --location --silent --show-error --http1.1 --connect-timeout 60 \
      --speed-limit 1024 --speed-time 600 --max-filesize 20000000000 -D "$dest.headers.$n" -C - -o "$dest.part" "$url"
    rc=$?
    echo "$name try $n rc=$rc bytes=$(stat -c %s "$dest.part" 2>/dev/null || echo 0) $(now)"
    [ $rc = 0 ] && break
    st=$(grep -E '^HTTP/' "$dest.headers.$n" 2>/dev/null | tail -1 | awk '{print $2}')
    echo "$name try $n last HTTP status: ${st:-none}"
    case "$st" in 403|404|410|451) fail "$name refused: HTTP $st";; esac
    case "$st" in 4*|5*) REFUSALS=$((REFUSALS + 1));; esac
    [ $REFUSALS -ge 5 ] && fail "$name: 5 refusals in a row (last HTTP $st)"
    if [ $rc = 33 ] || [ $rc = 36 ]; then echo "$name: no resume, next try from byte 0"; rm -f "$dest.part"; fi
    [ $n -ge 5 ] && fail "$name: gave up after $n tries (last rc=$rc)"
    sleep $((30 * (1 << (n - 1))))
  done
  REFUSALS=0
  mv "$dest.part" "$dest"
  md5=$(md5sum "$dest" | cut -d' ' -f1)
  echo "$name $3 $(stat -c %s "$dest") bytes md5 $md5 sha256 $(sha256sum "$dest" | cut -d' ' -f1) $(now)" | tee -a $D/SUMS
  if [ -n "$want" ]; then
    [ "$md5" = "$want" ] || fail "$name: md5 $md5 is not the provider's $want"
    echo "$name: md5 equals the provider's"
  fi
}

Z=https://zenodo.org/records/17402089/files
G=https://g-83fdd0.1beed.03c0.data.globus.org/cartuids/35846_b8b9149f2f59ff950549b5c89b8b35e2.tar

uhgv() {
  fetch uhgv_votus_fna "$Z/votus_full.fna.zst?download=1" uhgv/votus_full.fna.zst 9ad996fea08c6087a452194f0c2f9d33
  fetch uhgv_votus_meta "$Z/votus_metadata_extended.tsv.zst?download=1" uhgv/votus_metadata_extended.tsv.zst \
    6702355349e933fb605e62a8bbaa50a9
  fetch uhgv_biosample "$Z/source_biosample_metadata.tsv.zst?download=1" uhgv/source_biosample_metadata.tsv.zst \
    9934204a79cdac60ce9681a6233282b6
  fetch uhgv_host_range "$Z/host_range_breadth.tsv.zst?download=1" uhgv/host_range_breadth.tsv.zst \
    25960f854fa411956ef1353213133bfc
  fetch uhgv_genome_meta "$Z/uhgv_metadata.tsv.zst?download=1" uhgv/uhgv_metadata.tsv.zst f897bdf92b8b93cf3941f6d2c3d55efb
}

gsv() {  # run once on 2026-11-16: the HEAD was reset at the TLS hello and the README answered HTTP 400 three times
  if [ ! -e $D/gsv/gsv_atlas.tar ]; then
    echo "gsv HEAD $(now)"
    curl -A curl --fail --location --silent --show-error --http1.1 --connect-timeout 60 -I "$G" > $D/gsv/head.txt
    echo "gsv HEAD rc=$? $(grep -i -E '^(HTTP/|content-length|last-modified|etag)' $D/gsv/head.txt | tr -d '\r' | tr '\n' ' ')"
  fi
  fetch gsv_readme "https://data.pnnl.gov/sites/default/files/2023-12/README_Global%20Soil%20Virus%20%28GSV%29%20Atlas.txt" \
    gsv/README_GSV_Atlas.txt
  fetch gsv_archive "$G" gsv/gsv_atlas.tar
}

uhgv_faa() {  # the provider's proteins of the vOTU representatives, for the estimate of distinct proteins
  fetch uhgv_votus_faa "$Z/votus_full.faa.zst?download=1" uhgv/votus_full.faa.zst 1eff9626a8c3990a98d2b4ea58ff0596
}

vire() {
  fetch vire_meta "https://spire.embl.de/vire/api/download_metadata/?environment=all_samples&data_type=metadata&format=txt" \
    vire/vire_all_samples_metadata.txt
}

echo "c39_fetch start $(now) host $(hostname), sources: $*"
for s in "$@"; do
  case "$s" in
    uhgv) uhgv ;;
    uhgv_faa) uhgv_faa ;;
    gsv) gsv ;;
    vire) vire ;;
    *) fail "unknown source $s" ;;
  esac
done
now > $D/DONE
echo "c39_fetch done $(now)"
