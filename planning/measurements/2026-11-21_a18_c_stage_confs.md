# Stage driver configurations for `C1` to `C4` and `C6` (card A-18, DEC-62)

Date: 2026-11-21, wilkcomp02. Every number carries a label: *measured* (host, date, command or file), *derived* (formula and inputs) or *estimate*. **No build**: `/tmp/if2453/r10/root` was bound read-only (`:ro`) into every command, and nothing was written to the data root or `/stor`.

## 1. Read-only overlap reconciliation against `/tmp/if2453/r10/root`

Stages stored in `/tmp/if2453/r10/root/l0/contig_index` (*measured*, wilkcomp02, 2026-11-21, `a18_overlap.py` inside `bio-box_faec1f4.sif`, `is_default_variant` rows): `S1` (27,258), `P1` (539), `P2` (5,018), `P3` (9,107), `P4` (4,005), `P5` (4,111), `P6` (4,227), `H1` (6,714), `PL1` (61,647), `V1` (168,358), `C5` (2,653).

Reconciliation of `/tmp/if2453/s2/inputs/C1`, `C2`, `C3`, `C4`, `C6` (`genomes.json`) against all 11 stored stages and against `planning/measurements/2026-11-15_s2c_phase2_inputs.md` §3 (*measured*, wilkcomp02, 2026-11-21):

| Stage | Distinct contigs | `P1` | `P2` | `P3` | `P4` | `P5` | `P6` | `PL1` | `S1`/`H1`/`V1`/`C5` | Total held | Contigs added |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `C1` | 2,918 | 0 | 3 | 6 | 1 | 0 | 3 | 1,032 | 0 | 1,045 | 1,873 |
| `C2` | 6,675 | 0 | 0 | 0 | 3 | 5 | 4 | 3,220 | 0 | 3,232 | 3,443 |
| `C3` | 5,579 | 1 | 4 | 4 | 6 | 4 | 6 | 2,236 | 0 | 2,261 | 3,318 |
| `C4` | 4,571 | 0 | 5 | 5 | 4 | 0 | 0 | 1,918 | 0 | 1,932 | 2,639 |
| `C6` | 4,845 | 0 | 5 | 0 | 1 | 0 | 7 | 2,288 | 0 | 2,301 | 2,544 |

- **Contigs held by multiple stored stages**: 0 across all five stages (*measured*).
- **Reconciliation with `2026-11-15_s2c_phase2_inputs.md` §3**:
  - Every `P1`–`P6` and `PL1` count matches `2026-11-15_s2c_phase2_inputs.md` §3 (`analyze/summary.json` `stored_elsewhere_before_it.by_owner` and `pl1_stored`) exactly.
  - **`V1`**: holds 0 contigs of `C1`, `C2`, `C3`, `C4`, or `C6` (*measured*).
  - **`C5`**: holds 0 contigs of `C1`, `C2`, `C3`, `C4`, or `C6` (*measured*). The 1 contig (820 bp, `contig_uid` `368474c3475d9c3a80470a3bde270aea`, `NZ_AP025650.1`) shared between `C4` and `C5`'s inputs in `2026-11-15_s2c_phase2_inputs.md` §3 (`analyze/C4_contigs_elsewhere.tsv` line 397) is owned by `PL1`, which was built before `C5`; `C5` therefore skipped it at ingest (DEC-26) and `PL1` remains its sole owner (already counted in `PL1`'s 1,918 contigs of `C4`).
  - **Pairwise overlap among `C1`, `C2`, `C3`, `C4`, `C6`**: 0 contigs (*measured*). No `AFTER=` ordering is required among the five stages.
- **Input checksums**: `genomes.json` and `objects.json` sha256 in `/tmp/if2453/s2/inputs/C1`, `C2`, `C3`, `C4`, `C6` all match `2026-11-15_s2c_phase2_inputs.md` §2 (`inputs/SUMS`) byte for byte (*measured*, wilkcomp02, 2026-11-21).

## 2. `stage_driver.sh <STAGE> estimate` results

Each stage was checked read-only against `/tmp/if2453/r10/root` on `wilkcomp02` (2026-11-21, `STAGE_DIR=/tmp/if2453/a18/stages/c<k> bash /tmp/if2453/ci/A-18/scripts/stress/build_stage_pod/stage_driver.sh C<k> estimate`, inside `systemd-run --user --scope` with `MemoryMax=4G`, `CPUQuota=200%`, `nice -n 19`, at most 2 stages at once so <= 4 CPUs total; `stage_driver.sh` runs `stage_est.py before` inside `bio-box_faec1f4.sif` with `/tmp/if2453/r10/root:ro`).

All five stages exited 0 and wrote `estimate/DONE` with their conf's sha256 (`parked_owners: []`, `refusals: []`, `contigs_held_by_two_stages: 0`, every `EXPECT=held_by` equal to `owned_by_stage`).

| Stage | FASTA / `gz_bytes` (*measured*) | Held bases (*measured*) | Bases stored (*measured*) | Ingest CPU-min (batches?) (*derived*) | Gene calling CPU-min (batches?; called / not) (*derived*) | CRISPR default CPU-min (*derived*) | CRISPR near-Cas CPU-min (*derived*) | `estimate/DONE` conf sha256 (*measured*) |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| `C1` | 9,487,676,919 | 240,103,498 | 9,114,195,377 | 40.544 (no) | 151.864 (no; 1,764 / 12) | 25.535 | 29.365 | `2e2f68c2…d1a1` |
| `C2` | 9,489,508,963 | 227,238,138 | 9,126,873,544 | 40.605 (no) | 152.455 (no; 3,094 / 2) | 25.570 | 29.406 | `98c34e36…7a53` |
| `C3` | 9,621,864,080 | 635,727,847 | 8,808,991,551 | 39.505 (no) | 153.807 (no; 2,447 / 9) | 24.680 | 28.382 | `81120a73…052a` |
| `C4` | 9,347,446,902 | 411,856,959 | 8,806,059,807 | 39.226 (no) | 150.223 (no; 1,968 / 2) | 24.672 | 28.372 | `774963ad…9078` |
| `C6` | 9,426,417,771 | 137,715,892 | 9,166,710,932 | 40.697 (no) | 151.739 (no; 2,117 / 0) | 25.682 | 29.534 | `8d9ff10f…dce6` |

Formulas (*derived* by `stage_est.py before` from `biobox.stage` and `biobox.l1.genecall`, `margin_cpu_min = 240.0`):
- **Ingest** (card B-61 / A-17, all objects `.fna.gz`, no `min_contig_len`): `0.07 + max(0, fasta_bytes - held_bases) * 4.11 / 1e9 + gz_bytes * 0.26 / 1e9`. All five are 39.2–40.7 CPU-min, one process (`ingest_in_batches: false`), well under 240 CPU-min and above DEC-43's 24 CPU-min rate-guard floor.
- **Gene calling** (`genecall.cut_plan`, single mode, `STEP_OVERHEAD = 1.23`): 150.2–153.8 CPU-min in 1 process each (`merge_estimate_cpu_min` 1.025–1.043), identical to `2026-11-15_s2c_phase2_inputs.md` §4 (`genes_expected`: 8.45 M–8.65 M; largest genome 0.134–0.238 CPU-min).
- **CRISPR** (`stage.crispr_cpu_min`, `circular_stored = 0`, `flank_bases = 0`): default producer 24.672–25.682 CPU-min and near-Cas producer 28.372–29.534 CPU-min; both exceed DEC-43's 24 CPU-min threshold on every stage and will be judged by the rate guard.
- **Pfam 38.2 full scan (DEC-62)**: 9,961–12,428 CPU-min per stage in 42–52 shards (*estimate*, `2026-11-15_s2c_phase2_inputs.md` §4; before card B-47b's copy of proteins already measured in earlier stages).
