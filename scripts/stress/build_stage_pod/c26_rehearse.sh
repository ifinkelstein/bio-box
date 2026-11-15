#!/bin/bash
# Card C-26 part 1: the rehearsal of the PLSDB stage on a sample of its records, in a scratch root that links stages
# P1 and P2 of the data root (bound read-only in every launch); nothing is written in the data root. Copies of this
# script and c26_view.py, c26_sample.py, c26_rh_check.py, make_inputs.py run from /tmp/if2453/c26/bin. Each mode in
# its own small unit, for example
#
#   systemd-run --user --unit=c26-rh-driver-run --collect -p MemoryMax=4G -p MemorySwapMax=0 \
#     -p StandardOutput=append:/tmp/if2453/c26/driver.out -p StandardError=append:/tmp/if2453/c26/driver.out \
#     nice -n 19 bash /tmp/if2453/c26/bin/c26_rehearse.sh run
#
# How the stage-only launch ends (card C-26 as amended 2026-11-15). Once every other stage of the release is
# gene-called, build-stage's cluster gate (card B-36, `_cluster_gate`) no longer refuses, so a launch of the stage
# under the nine-stage release would go on into the release's clustering. No stage step's argv names the release
# (`stage.steps`), so the stage's files do not depend on which release the launch names. The launch is ended by a
# STOP file in the root, written by this driver as soon as the stage's last step (the filtered CRISPR producer,
# crispr-$NEAR) has published: its run record or its intent file (`<step>.publishing.json`, which run-step writes only
# after its own post-command STOP check, `runstep.run_step`) appears. Between that and the first release step's launch
# build-stage computes the cluster gate (protein counts from every stage's footers), dim's ownership check
# (`dim.owner_problem`, DuckDB over every stage's contig_index and genome_input), the parts decision, and in `one()`
# reads the CPU allowance (`runstep.allowance`, a sampling of the host's CPU) before run-step checks STOP once more
# before writing anything. So the STOP lands in that window; build-stage says `STOP exists; stopping before
# <release>.cluster...` and exits 75, and no file of the release is written. Polling every POLL seconds.
#
# Modes:
#   prep    The view (c26_view.py ROOT RH/root P1,P2), the sample (c26_sample.py from the stage's real inputs), its
#           inputs (make_inputs.py plasmid ... --root RH/root) and a scratch profile (8 threads, the real budget).
#   run     The launch, the real command's form for the sample stage PLR under release rC26.1 over P1,P2,PLR, with
#           --library pfam 38.2 <the real Pfam file> gathering --crispr; the STOP watcher. Expects exit 75 at the
#           release's cluster step, every PLR stage step recorded and nothing of rC26.1 in RH/root; then removes the
#           STOP it wrote (its own content only). Writes RH/run/OK.
#   dim     The dimension tables of a scratch release rC26.2 over P1,P2,PLR, through run-step (`release dim`), in the
#           view: one row per contig and per genome, the second genomes of the shared plasmids (DEC-27).
#   check   c26_rh_check.py: DEC-26 and DEC-27 on the rehearsal's files, read only.
set -u

C=/tmp/if2453/c26
RH=$C/rh
ROOT=/tmp/if2453/r10/root
TREE=${TREE:-/tmp/if2453/ci/C-26}
IMAGE=/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif
DIGEST=sha256:d67300aa1321b571fa8da79e74dfccdd87a8dca2e56518f324902efa6b9f098a
PFAMDIR=/stor/scratch/Wilke/if2453/all-the-data/l0/refdb/pfam/38.2
PFAM=$PFAMDIR/Pfam-A.hmm.gz
S2=/tmp/if2453/s2
IN=$S2/inputs/PL1
META=$S2/plsdb/meta
STAGE=PLR
REL=rC26.1
DIMREL=rC26.2
ALL=P1,P2,PLR
NEAR=minced-0.4.2-w8-nr2-r23x47-s26x50-cas1kb-caspfam1-pfam-38.2
T=${T:-8}
J=${J:-2}
CPUS=${CPUS:-4}
POLL=${POLL:-0.5}
BIN=$C/bin
LOGS=$C/logs
LOG=$C/driver.log
mkdir -p "$LOGS"

log() { echo "$(date '+%F %T.%3N %Z') $*" >> "$LOG"; }
fail() { log "FAILED: $*"; mkdir -p "$RH/$MODE"; echo "$*" > "$RH/$MODE/FAILED"; exit 1; }

py() {  # python in the runtime image, the tree's src on the path, the data root and the metadata bound read-only
  systemd-run --user --scope --quiet -p MemoryMax=16G -p MemorySwapMax=0 -p CPUQuota=200% nice -n 19 \
    singularity exec --cleanenv --env PYTHONPATH="$TREE/src" --env OMP_NUM_THREADS=2 --env PYTHONUNBUFFERED=1 \
    --bind /tmp/if2453 --bind "$ROOT:$ROOT:ro" --bind "$META:$META:ro" --bind /sys/fs/cgroup:/sys/fs/cgroup:ro \
    "$IMAGE" python "$@"
}

exit_of() { sed -n "s/^pod_run: $1 exited \([0-9]*\)$/\1/p" "$LOGS/$1.log" | tail -n 1; }

prep() {
  [ -e "$RH/inputs/genomes.json" ] && { log "prep: inputs present; nothing to do"; return; }
  python3 -I "$BIN/c26_view.py" "$ROOT" "$RH/root" P1,P2 >> "$LOG" 2>&1 || fail "prep: view"
  mkdir -p "$RH/inputs" "$RH/tmp" "$RH/durable" "$RH/bulk"
  python3 -I "$BIN/c26_sample.py" "$IN" "$META/nuccore.csv" "$IN/plsdb_2024_05_31_v2_sequences.fasta" \
    "$RH/inputs/sample.fasta" "$RH/inputs/sample.tsv" P1,P2 >> "$LOG" 2>&1 || fail "prep: sample"
  py "$BIN/make_inputs.py" plasmid "$RH/inputs/sample.fasta" "$META/nuccore.csv" "$RH/inputs" --root "$RH/root" \
    > "$RH/inputs/report.json" 2>> "$LOG" || fail "prep: make_inputs.py"
  log "prep: sample inputs: $(cat "$RH/inputs/report.json")"
  cat > "$RH/profile.yaml" << EOF
max_threads_day: $T
max_threads_night: $T
mem_scope: 32G
tmp: $RH/tmp
durable: $RH/durable
bulk: $RH/bulk
nice_monolithic: 19
nice_fanout: 15
day_start_hour: 8
day_end_hour: 20
cpu_min_per_process: 300
cluster_cpu_min_max: 20000
small_file_floor_bytes: 1000
require_mem_cap: true
EOF
  log "prep: done"
}

run() {
  local unit=c26-rh-1 rc seen="" t0 n
  [ -e "$RH/run/OK" ] && { log "run: OK present; nothing to do"; return; }
  [ -e "$RH/root/STOP" ] && fail "$RH/root/STOP exists before the launch"
  mkdir -p "$RH/run"
  find "$RH/root" ! -type d -printf '%P\t%s\t%y\n' | LC_ALL=C sort > "$RH/run/listing-before.tsv"  # links not followed
  log "launch $unit: build-stage --root $RH/root --stage $STAGE --release $REL --stages $ALL (Pfam 38.2, --crispr)"
  bash "$TREE/scripts/pod_run.sh" --profile "$RH/profile.yaml" --image "$IMAGE" --digest "$DIGEST" --mem 32G \
    --threads "$T" --bind /tmp/if2453 --bind "$ROOT:$ROOT:ro" --bind "$PFAMDIR:$PFAMDIR:ro" \
    --env "PYTHONPATH=$TREE/src" --log "$LOGS/$unit.log" "$unit" \
    build-stage --profile "$RH/profile.yaml" --root "$RH/root" --stage "$STAGE" --release "$REL" --stages "$ALL" \
    --genomes "$RH/inputs/genomes.json" --objects "$RH/inputs/objects.json" --library pfam 38.2 "$PFAM" gathering \
    --crispr --threads "$T" --scan-jobs "$J" --scan-cpus "$CPUS" >> "$LOG" 2>&1
  sleep 2
  while systemctl --user is-active --quiet "$unit"; do
    if [ -z "$seen" ] && { [ -e "$RH/root/runs/stage=$STAGE/crispr-$NEAR.json" ] \
         || [ -e "$RH/root/runs/stage=$STAGE/crispr-$NEAR.publishing.json" ]; }; then
      echo "c26 rehearsal watcher $(date '+%F %T.%N')" > "$RH/root/STOP"
      seen=1
      log "run: the last stage step's record or intent file appeared; wrote $RH/root/STOP"
    fi
    sleep "$POLL"
  done
  rc=$(exit_of "$unit")
  log "run: $unit exit '$rc'; $(grep -E 'build-stage: .*STOP|refused' "$LOGS/$unit.log" | tail -n 2 | tr '\n' ' ')"
  log "run: usage $(grep -E 'pod_run: .* usage' "$LOGS/$unit.log" | tail -n 1)"
  [ -n "$seen" ] || fail "the watcher never saw crispr-$NEAR's record (exit '$rc')"
  [ "$rc" = 75 ] || fail "$unit exit '$rc', not 75"
  grep -q "STOP exists; stopping before $REL\.cluster" "$LOGS/$unit.log" || fail "no 'stopping before $REL.cluster' line"
  for step in inputs ingest genecall "crispr-minced-0.4.2-w8-nr3-r23x47-s26x50" igr tapes domains-pfam-38.2 "crispr-$NEAR"; do
    [ -e "$RH/root/runs/stage=$STAGE/$step.json" ] || fail "no run record $STAGE/$step"
  done
  n=$(find "$RH/root" "$RH/tmp" -path '*release=rC26.1*' -o -path '*stage=rC26.1*' -o -path '*rC26.1.*' | wc -l)
  log "run: paths naming $REL in the view and its scratch: $n"
  [ "$n" = 0 ] || fail "$n paths of $REL"
  [ -e "$RH/root/release_def/release=$REL" ] && fail "$REL has a release_def"
  grep -q "^c26 rehearsal watcher" "$RH/root/STOP" && rm "$RH/root/STOP"
  log "run: removed the watcher's STOP"
  touch "$RH/run/OK"
  log "run: OK"
}

dim() {
  local unit=c26-rh-dim
  [ -e "$RH/root/runs/stage=$DIMREL/dim.json" ] && { log "dim: record present"; return; }
  mkdir -p "$RH/tmp/scratch/$DIMREL.dim"
  bash "$TREE/scripts/pod_run.sh" --profile "$RH/profile.yaml" --image "$IMAGE" --digest "$DIGEST" --mem 16G \
    --threads 2 --bind /tmp/if2453 --bind "$ROOT:$ROOT:ro" --env "PYTHONPATH=$TREE/src" --log "$LOGS/$unit.log" \
    "$unit" run-step --profile "$RH/profile.yaml" --root "$RH/root" --stage "$DIMREL" --step dim --threads 2 -- \
    python -m biobox.cli release dim --root "$RH/root" --release "$DIMREL" --stages "$ALL" \
    --tmp "$RH/tmp/scratch/$DIMREL.dim" >> "$LOG" 2>&1
  sleep 5
  while systemctl --user is-active --quiet "$unit"; do sleep 5; done
  log "dim: $unit exit '$(exit_of "$unit")': $(grep -E 'published|failed|refused' "$LOGS/$unit.log" | tail -n 1)"
  [ "$(exit_of "$unit")" = 0 ] || fail "dim: exit '$(exit_of "$unit")'"
}

check() {
  py "$BIN/c26_rh_check.py" --root "$RH/root" --data-root "$ROOT" --inputs "$RH/inputs" --stage "$STAGE" \
    --release "$DIMREL" --linked P1,P2 --log "$LOGS/c26-rh-1.log" --out "$RH/check" > "$RH/check.out" 2>&1
  local rc=$?
  log "check: exit $rc: $(tail -n 1 "$RH/check.out" | cut -c1-600)"
  return $rc
}

MODE=${1:-}
log "c26 rehearsal: mode '$MODE' (pid $$), tree $TREE, image $IMAGE"
case "$MODE" in
  prep) prep ;;
  run) run ;;
  dim) dim ;;
  check) check ;;
  *) echo "usage: c26_rehearse.sh prep|run|dim|check" >&2; exit 2 ;;
esac
log "c26 rehearsal: mode '$MODE' done"
