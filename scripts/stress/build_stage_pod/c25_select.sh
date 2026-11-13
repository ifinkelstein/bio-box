#!/bin/sh
# Card C-25 (1): the S2C selection on the pod, one CPU (planning/measurements/2026-11-12_s2c_selection.md).
#
#   systemd-run --user --unit c25-select --collect -p MemoryMax=8G -p MemorySwapMax=0 -p CPUQuota=100% \
#     -p StandardOutput=append:/tmp/if2453/c25/logs/select.log -p StandardError=append:/tmp/if2453/c25/logs/select.log \
#     nice -n 19 sh /tmp/if2453/c25/bin/c25_select.sh
#
# Reads the durable GTDB R226 metadata archive on /stor (read only), unpacks it under /tmp/if2453/c25/gtdb, runs
# c25_select.py into /tmp/if2453/c25/out, and reruns phase 1's r50_select.py (whose reader c25_select.py shares)
# into /tmp/if2453/c25/r50_recheck to show that its all.tsv is unchanged. Writes nothing outside /tmp/if2453/c25.
set -eu
W=/tmp/if2453/c25
TAR=/stor/scratch/Wilke/if2453/all-the-data/provisional/s1_inputs/gtdb_r226_metadata.tar
NOT_SERVED=GCA_900637225.1,GCA_902153285.2,GCA_902502775.2,GCA_902502895.2,GCA_963667335.1,GCA_963921365.2,GCF_040719245.2
echo "start $(date -Is) host $(hostname -s) python $(python3 --version 2>&1)"
sha256sum "$TAR"
tar -xf "$TAR" -C "$W/gtdb"
sha256sum "$W/gtdb/bac120_metadata_r226.tsv.gz" "$W/gtdb/ar53_metadata_r226.tsv.gz"
rm -rf "$W/out.tmp"
/usr/bin/time -v python3 -I "$W/bin/c25_select.py" "$W/gtdb" "$W/phase1/all.tsv" "$W/out.tmp" \
  --cap 20 --phase1-not-served "$NOT_SERVED" > "$W/logs/selection.stdout.json"
rm -rf "$W/out"
mv "$W/out.tmp" "$W/out"
sha256sum "$W/out/s2c_selected.tsv" "$W/out/remainder.tsv" "$W/out/by_phylum.tsv" "$W/out/selection.json"
rm -rf "$W/r50_recheck"
python3 -I "$W/bin/r50_select.py" "$W/gtdb" "$W/ladder_sel" /tmp/if2453/ladder/genomes "$W/r50_recheck" > /dev/null
cmp "$W/r50_recheck/all.tsv" "$W/phase1/all.tsv"
echo "r50_select all.tsv unchanged"
echo "end $(date -Is)"
