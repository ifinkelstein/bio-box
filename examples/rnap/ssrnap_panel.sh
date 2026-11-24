#!/bin/bash
# Companion two-namespace single-subunit & single-polypeptide viral RNAP super-panel check (card W-40, DEC-69).
# Leaves the locked PF00940 release check (`run.sh`, `expected.tsv`, `recorded/*.tsv`) untouched.
#
# Usage on wilkcomp02 (read-only against ROOT):
#   TREE=... PROFILE=... IMAGE=... ROOT=/tmp/if2453/r10/root RELEASE=r26.11.2 \
#   HMMS=/tmp/if2453/w40/hmms OUT=/tmp/if2453/w40/panel_out \
#     bash "$TREE/examples/rnap/ssrnap_panel.sh"
set -euo pipefail

: "${TREE:?set TREE}" "${PROFILE:?set PROFILE}" "${IMAGE:?set IMAGE}"
: "${ROOT:?set ROOT}" "${RELEASE:?set RELEASE}" "${HMMS:?set HMMS}" "${OUT:?set OUT}"
STAGES=${STAGES:-S1,H1,PL1,P1,C5,V1}
THREADS=${THREADS:-4}
MEM_GB=${MEM_GB:-28}
CUTOFFS=${CUTOFFS:-$TREE/examples/rnap/ssrnap_cutoffs.tsv}

[ ! -e "$ROOT/STOP" ] || { echo "STOP file present in $ROOT" >&2; exit 2; }
mkdir -p "$OUT/runs" "$OUT/work"

exec singularity exec --cleanenv \
  --env "PYTHONPATH=$TREE/src" \
  --env "OMP_NUM_THREADS=$THREADS" \
  --bind "$ROOT:$ROOT:ro" \
  --bind /tmp/if2453/bio-box:/tmp/if2453/bio-box:ro \
  --bind "$TREE:$TREE:ro" \
  --bind "$HMMS:$HMMS:ro" \
  --bind "$OUT:$OUT" \
  --bind /stor/scratch/Wilke/if2453/all-the-data:/stor/scratch/Wilke/if2453/all-the-data:ro \
  --bind /sys/fs/cgroup:/sys/fs/cgroup:ro \
  "$IMAGE" \
  python -m biobox.cli search panel \
    --root "$ROOT" \
    --release "$RELEASE" \
    --profile "$PROFILE" \
    --runs "$OUT/runs" \
    --work "$OUT/work" \
    --stages "$STAGES" \
    --library palm@pfam:PF00940,PF14700,PF29426,PF21769,PF21867,PF21894,PF21624 \
    --hmm "palm:$HMMS/vogdb_rnap_palm_clean.hmm" \
    --hmm "palm:$HMMS/phrogs_rnap_palm_clean.hmm" \
    --hmm "palm:$HMMS/seeds_palm.hmm" \
    --library dpbb_single@pfam:PF05183 \
    --hmm "dpbb_single:$HMMS/seeds_dpbb.hmm" \
    --decoy-library decoy_pol1@pfam:PF00476 \
    --decoy-library decoy_cell_rnap@pfam:PF04997,PF00623,PF04983,PF05000,PF04998,PF04992,PF04990,PF04563,PF04561,PF04565,PF04566,PF10385,PF04567,PF00562,PF04560 \
    --cutoff gathering \
    --cutoffs-tsv "$CUTOFFS" \
    --positive palm:PF00940.25=2eef240afefddcd2c8ac7660eb50e384 \
    --positive palm:PF00940.25=eb8377fda717b96aac266fda63a2b096 \
    --positive palm:PF00940.25=42c5546422cdf8424b5a0db0d0bdc6f8 \
    --positive palm:PF00940.25=bf8803f4abe8a136c0d132ad2f31882b \
    --positive palm:PF21769.5=1e396b2d4e3d7a5280a6835bcf18e10d \
    --positive palm:N4_RNAP2_gp15=a17b8fbf735e80abf6546c9a641ee2fb \
    --positive palm:N4_RNAP2_gp16=eda1349c1f642a9294c36fda6fd0793a \
    --positive dpbb_single:crass_vrnap_dpbb=835ad638bcf19577465c5c2792dcb759 \
    --positive dpbb_single:yono_dpbb=82bbc9c36fa5629026f5ea0a668ba1be \
    --negative shuffled_seed \
    --max-gap-genes 2 \
    --min-models 1 \
    --neighbors 5 \
    --cpus "$THREADS" \
    --threads "$THREADS" \
    --mem-gb "$MEM_GB" \
    --out "$OUT/panel"
