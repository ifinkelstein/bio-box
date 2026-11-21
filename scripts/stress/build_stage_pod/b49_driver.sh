#!/bin/bash
# Card B-49 (DEC-45 (2)) on the pod: a scratch release over the seven stages of rM.2 that takes rM.2's clustering.
# Host side, bash; one mode per call, each a short command or one detached unit:
#
#   bash b49_driver.sh prep     check free space (stop under 205 GB), make the view C/view (b49_view.py: file links
#                               over ROOT for S1, P1 to P6, rM.1 and rM.2) and the scratch profile C/profile.yaml
#                               (the tree's wilkcomp02 profile with tmp C/tmp)
#   bash b49_driver.sh build    one unit b49-build (pod_run.sh: the runtime image of faec1f4 by digest, the tree's
#                               src first, ROOT bound read-only, --threads T, --mem 32G): build-stage of stage S1 and
#                               release rB49.1 over the seven stages in C/view, --cluster-parts 24 as rM.2's plan
#                               (its cut also holds the profile's cluster_cpu_min_max), --threads T. The stage steps
#                               are passed over (their run records are linked); the release steps run: cluster (the
#                               reuse), tokens, postings, dim, manifest, catalog
#   bash b49_driver.sh verify   `release verify` of rB49.1 in C/view, in the image, in a capped scope (ROOT read-only)
#   bash b49_driver.sh check    b49_check.py: the copies' sha256 against rM.2's run record, the run records, the
#                               ledger and the reuse record, the view's real bytes; then `du` of C
#   bash b49_driver.sh clean    remove C/view, C/tmp (C/logs and C/*.json stay)
# Writes only under C. Nothing in ROOT: it is bound read-only into every command run in the image.
set -u

C=/tmp/if2453/b49
ROOT=/tmp/if2453/r10/root
TREE=${TREE:-/tmp/if2453/ci/B-49}
IMAGE=/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif
DIGEST=sha256:d67300aa1321b571fa8da79e74dfccdd87a8dca2e56518f324902efa6b9f098a
REL=rB49.1
STAGES=S1,P1,P2,P3,P4,P5,P6
T=${T:-4}
MODE=${1:-}
GCV=pyrodigal-3.7.1-gv-0.3.2-p2

img() {  # a short command in the image, capped, nice 19, ROOT read-only
  systemd-run --user --scope --quiet -p MemoryMax=16G -p MemorySwapMax=0 -p CPUQuota=200% nice -n 19 \
    singularity exec --cleanenv --env PYTHONPATH="$TREE/src" --env OMP_NUM_THREADS=2 --env PYTHONUNBUFFERED=1 \
    --bind /tmp/if2453 --bind "$ROOT:$ROOT:ro" --bind /sys/fs/cgroup:/sys/fs/cgroup:ro "$IMAGE" "$@"
}

case "$MODE" in
  prep)
    free=$(df -B1 --output=avail /tmp | tail -n 1)
    echo "free on /tmp: $free bytes"
    [ "$free" -ge 205000000000 ] || { echo "under 205 GB free: stop"; exit 1; }
    mkdir -p "$C/tmp" "$C/logs"
    python3 -I "$C/bin/b49_view.py" "$ROOT" "$C/view" S1,P1,P2,P3,P4,P5,P6 rM.1,rM.2 || exit 1
    sed "s#^tmp: .*#tmp: $C/tmp#" "$TREE/profiles/wilkcomp02.yaml" > "$C/profile.yaml"
    grep -n '^tmp:' "$C/profile.yaml"
    ;;
  build)
    [ -d "$C/view" ] || { echo "no view: run prep"; exit 1; }
    bash "$TREE/scripts/pod_run.sh" --profile "$C/profile.yaml" --image "$IMAGE" --digest "$DIGEST" --mem 32G \
      --threads "$((T + 2))" --bind /tmp/if2453 --bind "$ROOT:$ROOT:ro" --env "PYTHONPATH=$TREE/src" \
      --log "$C/logs/b49-build.log" b49-build \
      build-stage --profile "$C/profile.yaml" --root "$C/view" --stage S1 --release "$REL" --stages "$STAGES" \
      --genomes /dev/null --objects /dev/null --library pfam 38.2 /dev/null gathering --threads "$T" \
      --cluster-parts 24 --tmp "$C/tmp/scratch" --gene-call-version "$GCV"
    ;;
  verify)
    img python -m biobox.cli release verify --root "$C/view" --release "$REL"
    ;;
  check)
    img python -I "$C/bin/b49_check.py" "$C/view" "$REL" rM.2 || exit $?  # not `du`'s status (card A-11)
    du -sb "$C"
    ;;
  clean)
    rm -rf "$C/view" "$C/tmp"
    du -sb "$C"
    ;;
  *)
    echo "usage: b49_driver.sh prep|build|verify|check|clean"
    exit 2
    ;;
esac
