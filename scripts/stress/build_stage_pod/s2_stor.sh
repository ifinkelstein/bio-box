#!/bin/sh
# Permanent copies of the S2 downloads and of stage H1's inputs (run ON THE POD):
#
#     systemd-run --user --unit=s2-stor --collect -p MemoryMax=8G -p MemorySwapMax=0 -p CPUQuota=100% \
#       nice -n 19 ionice -c 2 -n 7 sh /tmp/if2453/s2/bin/s2_stor.sh
#
# Into /stor/scratch/Wilke/if2453/all-the-data/provisional/s2_inputs/: new files only (an existing name is
# never overwritten), each at least 64 MB, each written as .<name>.tmp, its sha256 read back with
# O_DIRECT (from the server, not this host's page cache) and compared with the source's, then renamed
# with `mv -n`. One line per file in MANIFEST.sha256 (written last, the same way). Nothing under /stor
# is deleted: a copy that fails leaves its .tmp and stops the run. Log: logs/stor.log.
#
# s2_inputs.tar bundles /tmp/if2453/s2/{inputs, genomes/H1, work/{hpgc,plsdb,stage_uids.tsv}, logs, bin,
# dl/SUMS, dl/*.headers.*}: H1's inputs with the FASTA object they name, and the records of this work.
set -u
W=/tmp/if2453/s2
D=/stor/scratch/Wilke/if2453/all-the-data/provisional/s2_inputs
S=$W/stor_stage
exec >> $W/logs/stor.log 2>&1 < /dev/null
for f in /tmp/if2453/r10/root/STOP $W/STOP; do [ -e "$f" ] && { echo "STOP present: $f"; exit 3; }; done
echo "stor start $(date)"
mkdir -p $S
if [ ! -e $S/s2_inputs.tar ]; then
  cd $W || exit 1
  tar --sort=name --exclude=logs/stor.log --exclude=__pycache__ -cf $S/s2_inputs.tar.tmp inputs genomes/H1 work/hpgc work/plsdb work/stage_uids.tsv logs bin \
    dl/SUMS dl/*.headers.* || exit 1
  mv $S/s2_inputs.tar.tmp $S/s2_inputs.tar
  tar -tvf $S/s2_inputs.tar
fi
mkdir -p $D
[ -e $D/MANIFEST.sha256 ] && { echo "refused: $D/MANIFEST.sha256 exists"; exit 1; }
: > $S/MANIFEST.sha256

copy() {  # source, name on /stor
  src=$1; name=$2
  if [ -e "$D/$name" ] || [ -e "$D/.$name.tmp" ]; then echo "refused: $D/$name or its .tmp exists"; exit 1; fi
  size=$(stat -c %s "$src")
  [ "$size" -ge 67108864 ] || { echo "refused: $src is $size bytes, under 64 MB"; exit 1; }
  want=$(sha256sum "$src" | cut -d' ' -f1)
  echo "copy $src -> $D/$name ($size bytes, sha256 $want) $(date)"
  cp "$src" "$D/.$name.tmp" || exit 1
  got=$(dd if="$D/.$name.tmp" iflag=direct bs=8M status=none | sha256sum | cut -d' ' -f1)
  [ "$got" = "$want" ] || { echo "FAILED: $D/.$name.tmp sha256 $got, source $want"; exit 1; }
  mv -n "$D/.$name.tmp" "$D/$name"
  [ -e "$D/.$name.tmp" ] && { echo "FAILED: rename of $D/.$name.tmp refused"; exit 1; }
  echo "$got  $name" >> $S/MANIFEST.sha256
  echo "verified $name sha256 $got (read back with O_DIRECT) $(date)"
}

copy $W/dl/HPGC.7230.genomes.final.fasta.gz HPGC.7230.genomes.final.fasta.gz
copy $W/dl/plsdb_download_fasta plsdb_2024_05_31_v2_sequences.fasta
copy $W/dl/download_meta.tar.gz plsdb_2024_05_31_v2_meta.tar.gz
copy $S/s2_inputs.tar s2_inputs.tar
cp $S/MANIFEST.sha256 $D/.MANIFEST.sha256.tmp && mv -n $D/.MANIFEST.sha256.tmp $D/MANIFEST.sha256
cd $D && sha256sum -c MANIFEST.sha256
echo "stor done rc=$? $(date)"
