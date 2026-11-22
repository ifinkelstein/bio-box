# Card B-66: `O(1)`-stage `SeqStore` point lookups and per-stored-contig CRISPR CPU-minute term

Date: 2026-11-22, `wilkcomp02`. Every number carries a label: *measured* (host, date, command/file), *derived* (formula and inputs), or *estimate*. **Nothing was built in the data root**: `/tmp/if2453/r10/root` and `/tmp/if2453/s2/inputs/SPX1` were bound read-only (`:ro`), and all scratch ran under `/tmp/if2453/b66/` (removed after verification).

## 1. Summary

1. **`O(1)`-stage `SeqStore` point lookups and single-locate `get`**: `KeyLookup.find_first(key, predicate=None)` stops at the first matching file/row group without decoding later stages, and `SeqStore.get` locates `raw_digest` once via `_where(rd)` instead of calling `length_of(rd)` + `raw(rd)`. On 20,000 `G1` contigs (`165,408,771` bp) with all 12 root stages opened (`b66_crispr_contig_arms.py seqstore`, *measured*, wilkcomp02, 2026-11-22), `SeqStore.get` dropped from **`1.4881` to `0.5627` ms/contig (`-62.2%`)**, cutting `contig_index` row groups read from `12` to `1` and `seq_locator` row groups read from `22` to `10` (the 10 row groups of `G1` alone, zero reads of the other 11 stages).
2. **Per-stored-contig CRISPR term (`CRISPR_CPU_S_PER_K_CONTIGS = 0.50` CPU-s / 1,000 stored contigs = `8.333` CPU-min / M contigs)**: fitted from post-fix end-to-end `crispr.run` and `crispr.run_near_cas` (`cpus=4`) on `G1` 1-in-20 (`61,361` contigs, `510.04` Mbp, 12 stages) and `SPX1` 1-in-40 (`49,153` contigs, `730.72` Mbp, 13 stages). The per-contig term is **not** multiplied by `CRISPR_NEAR_CAS_COST = 1.15` because per-contig sequence/FASTA I/O is identical in both producers (confirmed by `G1` near-Cas `128.034` vs `121.072` CPU-s default, `+5.8%`, giving `-1.0%` error without multiplying vs `-4.4%` with multiplying, and improving `V1` near-Cas from `-9.4%` to `-8.6%`).
3. **Fit across all stages**:
   - All 10 non-`C5` complete stages (`S1, P1..P6, H1, PL1, V1`, 20 records) sit within **`10.3%`** (`0.897×` to `1.081×`, well inside the `15%` target). `C5` remains documented under `DEC-53` (`0.836×` default, `0.758×` near-Cas).
   - `G1` (`1,227,214` contigs, `10.043` Gbp): new estimate **`38.363` CPU-min** default / **`42.583` CPU-min** near-Cas (vs pre-fix default record `66.117` CPU-min; post-fix sample projects `39.92` CPU-min default `+4.1%` and `42.17` CPU-min near-Cas `-1.0%`).
   - `SPX1` (`2,829,368` contigs, `35.885` Gbp): new estimate **`124.116` CPU-min** (`< 240` CPU-min margin); 1-in-40 sample runs at `2.3149` vs `2.4569` CPU-min (`0.942×`, `-5.8%`, within `15%`), projecting `~116.0` CPU-min full-stage.
   - `PL2` (`112,579` contigs, 11-stage view) and `PL3` (`747,644` contigs, 11-stage view): **`8.486` / `9.618` CPU-min** (`PL2`) and **`60.110` / `68.192` CPU-min** (`PL3`).

---

## 2. Before / after `SeqStore.get` on `G1` (12 stages opened)

Command: `python scripts/stress/b66_crispr_contig_arms.py seqstore --n 20000` (unit `b66-seqstore2`, `CPUQuota=200%`, `MemoryMax=16G`, `nice -n 19`, `/tmp/if2453/r10/root:ro`, 2026-11-22, *measured*):

Stages opened (12): `['G1', 'C5', 'H1', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'PL1', 'S1', 'V1']`.
Contigs queried: first `20,000` default-variant contigs of `G1` (`165,408,771` bp).

| Mode | Wall (s) | CPU (s) | ms / contig | `index.row_groups_read` | `locator.row_groups_read` | Projected `G1` (`1.227M` contigs) |
|---|---:|---:|---:|---:|---:|---:|
| Before (`get` + double `_where`) | 29.7694 (*measured*) | 29.7626 (*measured*) | 1.4881 (*derived*) | 12 (*measured*) | 22 (*measured*) | 30.44 CPU-min (*derived*) |
| After (`find_first` + single `_where`) | 11.2569 (*measured*) | 11.2546 (*measured*) | 0.5627 (*derived*) | 1 (*measured*) | 10 (*measured*) | 11.51 CPU-min (*derived*) |
| **Savings** | **−18.5125 s (−62.2%)** | **−18.5080 s (−62.2%)** | **−0.9254 ms** | **−11** | **−12** | **−18.93 CPU-min** (*derived*) |

`G1`'s own `seq_locator` has 10 row groups (`1,227,214 / 122,880`); 20,000 contigs hit all 10 of `G1`'s row groups and zero row groups of the other 11 stages.

---

## 3. Post-fix end-to-end CRISPR scratch arms (`G1` 1-in-20 and `SPX1` 1-in-40)

Commands (each `cpus=4`, `CPUQuota=400%`, `MemoryMax=16G`, `nice -n 19`, `/tmp/if2453/r10/root:ro`, 2026-11-22, *measured*):
- `b66-g1`: `python scripts/stress/b66_crispr_contig_arms.py g1 --stride 20 --cpus 4` (12 stages opened)
- `b66-spx1`: `python scripts/stress/b66_crispr_contig_arms.py spx1 --stride 40 --max-objects 8 --cpus 4` (13 stages opened; `min_contig_len=5000`; `8` of `304` `SPX1` objects)

| Arm | Stages | Contigs stored | Bases stored (bp) | Producer | Wall (s) | CPU (s) | CPU-min | Base term (CPU-s, *derived*) | Residual (CPU-s / k-contigs, *derived*) | New est (CPU-min, `c=0.50`) | Ratio (`meas/est`, *derived*) |
|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|
| `G1` 1-in-20 | 12 | 61,361 (*measured*) | 510,039,671 (*measured*) | default | 50.797 (*measured*) | 121.072 (*measured*) | 2.0179 (*measured*) | 85.738 | 0.5758 | 1.9403 (*derived*) | **1.040** (+4.0%) |
| `G1` 1-in-20 | 12 | 61,361 (*measured*) | 510,039,671 (*measured*) | near-Cas | 52.915 (*measured*) | 128.034 (*measured*) | 2.1339 (*measured*) | 98.598 | 0.4797 | 2.1546 (*derived*) | **0.990** (−1.0%) |
| `SPX1` 1-in-40 | 13 | 49,153 (*measured*) | 730,724,623 (*measured*) | default | 45.246 (*measured*) | 138.896 (*measured*) | 2.3149 (*measured*) | 122.835 | 0.3268 | 2.4569 (*derived*) | **0.942** (−5.8%) |

### Fit and `near_cas` decision

- The three post-fix residual per-contig rates are `0.5758`, `0.4797`, and `0.3268` CPU-s / 1,000 contigs (mean `0.461`, relative-error least-squares fit `~0.48`, *derived*). We choose **`CRISPR_CPU_S_PER_K_CONTIGS = 0.50`** CPU-s / 1,000 stored contigs (`8.333` CPU-min / M contigs), slightly conservative against `DEC-64` over-estimate stops while keeping every arm and all 10 non-`C5` built stages well within `15%`.
- **`near_cas` does not multiply the contig term**: on `G1` 1-in-20, `near_cas - default = +6.962 CPU-s` (`+5.8%`), matching `+15%` of the base term alone (`+12.861 CPU-s`, total error `-1.0%`) much better than multiplying both terms (`+18.161 CPU-s`, total error `-4.4%`). Mechanistically, `SeqStore.get`, BGZF decompression, `ContigOwners`, and FASTA batch formatting are identical in `crispr.run` and `crispr.run_near_cas`.

---

## 4. All 11 complete stages (`S1..V1, C5`) with `CRISPR_CPU_S_PER_K_CONTIGS = 0.50`

Run records read from `/tmp/if2453/r10/root/runs/stage=<S>/crispr-*.json` (`b66_crispr_contig_arms.py summary`, wilkcomp02, 2026-11-22, *measured*). Estimates *derived* from `stage.crispr_cpu_min(bases, flank, contigs=n_stored, near_cas=...)`.

| Stage | Contigs stored (*measured*) | Scanned Gbp (*measured*) | Default meas (CPU-min) | Default est (CPU-min) | Default ratio | Near-Cas meas (CPU-min) | Near-Cas est (CPU-min) | Near-Cas ratio |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `S1` | 27,258 | 1.664 | 5.085 | 4.889 | 1.040 | 5.487 | 5.589 | 0.982 |
| `P1` | 539 | 1.017 | 3.085 | 2.853 | 1.081 | 3.383 | 3.280 | 1.031 |
| `P2` | 5,018 | 9.596 | 26.489 | 26.927 | 0.984 | 31.794 | 30.960 | 1.027 |
| `P3` | 9,107 | 9.949 | 29.195 | 27.950 | 1.045 | 33.386 | 32.131 | 1.039 |
| `P4` | 4,005 | 9.976 | 27.585 | 27.984 | 0.986 | 32.606 | 32.177 | 1.013 |
| `P5` | 4,111 | 9.968 | 27.505 | 27.960 | 0.984 | 32.361 | 32.149 | 1.007 |
| `P6` | 4,227 | 9.963 | 27.938 | 27.947 | 1.000 | 32.796 | 32.134 | 1.021 |
| `H1` | 6,714 | 1.820 | 4.925 | 5.155 | 0.955 | 5.876 | 5.919 | 0.993 |
| `PL1` | 61,647 | 12.130 | 30.938 | 34.498 | 0.897 | 35.513 | 39.595 | 0.897 |
| `V1` | 168,358 | 6.997 | 20.560 | 21.007 | 0.979 | 21.893 | 23.948 | 0.914 |
| `C5` (`DEC-53`) | 2,653 | 9.109 | 21.361 | 25.543 | 0.836 | 22.257 | 29.371 | 0.758 |

All 10 non-`C5` complete stages (`S1..V1`, 20 records) lie between **`0.897×` and `1.081×`** (within `10.3%`, inside the `15%` target). `C5` stays documented under `DEC-53`.

---

## 5. Upcoming stages: `G1`, `SPX1`, `PL2`, `PL3`

All estimates *derived* from `stage.crispr_cpu_min(bases, flank, contigs=n_stored, near_cas=...)` with `CRISPR_CPU_S_PER_GBP = 168.1`, `CRISPR_NEAR_CAS_COST = 1.15`, `CRISPR_CPU_S_PER_K_CONTIGS = 0.50`.

| Stage | Contigs stored | Bases stored (bp) | Flank bases (bp) | Old B-53 est (def / nc, CPU-min) | New B-66 est (def / nc, CPU-min) | Notes |
|---|---:|---:|---:|---:|---:|---|
| `G1` | 1,227,214 (*measured*, root) | 10,042,502,859 (*measured*, root) | 0 | 28.136 / 32.356 | **38.363 / 42.583** | Pre-fix default record `66.117` (*measured*, root) included `~18.9` CPU-min of O(12)-stage `SeqStore.get` overhead; post-fix 1-in-20 sample projects **`39.92` default (`1.041×`) / `42.17` near-Cas (`0.990×`)** (*derived*) |
| `SPX1` | 2,829,368 (*measured*, staged) | 35,885,000,000 (*measured*, staged) | 0 | 100.538 / — | **124.116 / —** | `called` depth (default producer only); `< 240` CPU-min margin; 1-in-40 sample runs at `0.942×` (`-5.8%`, within `15%`) and projects **`~116.0` CPU-min** (*derived*) |
| `PL2` (11-stage view) | 112,579 (*measured*, C-52 `before.json`) | 1,976,283,378 (*measured*) | 717,808,262 (*measured*) | 7.548 / 8.680 | **8.486 / 9.618** | Build order (`112,562` contigs, `1.975` Gbp) differs by `< 0.01` CPU-min (*derived*) |
| `PL3` (11-stage view) | 747,644 (*measured*, C-52 `before.json`) | 16,409,762,458 (*measured*) | 2,821,603,614 (*measured*) | 53.880 / 61.962 | **60.110 / 68.192** | Build order (`745,871` contigs, `16.283` Gbp) ~`59.7 / 67.8` CPU-min (*derived*); well under `240` CPU-min |
