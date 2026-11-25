# `W-38a` Pre-Build Protocol: GPU pLM & `3Di` Remote-Homology Bake-Off Gate (`SD-21` Extension)

**Status:** Ready-to-execute protocol & harness specification (zero GPU compute executed; awaits PI approval of `W-38a` / `DEC-PLM-1..2`).  
**Target Hardware:** `wilkcomp03.ccbb.utexas.edu` (`4× NVIDIA RTX 6000 Ada Generation`, `48 GB` VRAM each = `192 GB` total VRAM; `*measured*`, `2026-11-19`, `nvidia-smi -L`).  
**Target Corpus Release:** Parameterized via `--release` (default `r26.11.2`; compatible with `r26.11.3` and any subsequent rolling release). **Note (`DEC-66`):** `r26.11.1` is scheduled for archive retirement under `S-4` (`N = 2` newest rolling releases kept hot on `/tmp`), so this protocol never hardcodes `r26.11.1` and reads release tables strictly through `biobox.catalog.ReleaseFiles.of(root, release)` (`AGENTS.md` §3.9).  
**Harness & Unit Tests:** [`scripts/stress/w38a_plm_bakeoff.py`](../../scripts/stress/w38a_plm_bakeoff.py) and [`tests/r1/test_w38a_plm_bakeoff.py`](../../tests/r1/test_w38a_plm_bakeoff.py).  
**Primary Design Context:** [`planning/design/pLM_remote_system_search.md`](pLM_remote_system_search.md) (§1–§5), [`planning/search/SEARCH_PLAN.md`](../search/SEARCH_PLAN.md) (`SD-21`), and [`planning/tasks/TASKS.md`](../tasks/TASKS.md) (`W-38`, `W-38a`).

**Number Labels (`AGENTS.md` §6):**
- `*published*` — cited paper, repository source line, or vendor specification;
- `*measured*` — observed directly on `wilkcomp02` / `wilkcomp03`;
- `*derived*` — explicit arithmetic formula from labeled inputs;
- `*estimate*` — stated modeling assumption.

---

## 1. Purpose & Falsifiability Contract

Before committing TACC GPU allocations (`Tier 2` `0.92 Tbp` / `Tier 3` `14.6 Tbp`) or registering new production tables in `src/biobox/registry/schemas.py`, `W-38a` executes a single bounded GPU bake-off (`<10` GPU-hours on `wilkcomp03`) that tests three concrete architectural decisions against quantitative Go/No-Go gates:

1. **Gate 1 (`Option A`, `SD-21`):** Does `3Di` (`ProstT5` via Foldseek v10 vs `ESM-2 3B 3Di` via Domainator) recover `>= +20` percentage points higher `Recall@1%FPR` than `3-round jackhmmer` on `<20%` AA ID & `hmmsearch E > 10` remote homologs and multi-domain fusions at `<= 1.5×` Foldseek-AFDB latency?
2. **Gate 2 (`Option B`, Corpus-Wide Dense ANN):** Does `ESM-C 300M` pooled ANN (`fp16` / `IVF-PQ 64B`, full-ORF vs `segments.faa`) recover `>= +10%` marginal unique true-positive systems at `FPR <= 0.1%`/genome **beyond** `jackhmmer + HHblits + 3Di-Foldseek` with `< 5%` fusion recall loss? If **no**, Option B is **rejected permanently**.
3. **Gate 3 (`Option C` Tape-Context Scorer vs `Option D` On-Demand Rounds):** Do supervised operon-context scorers (`DefensePredictor` patched for circular wraparound & 1-CDS `NaN` bug; `ESM_DF + GeneCLR_DF`) achieve `>= 2×` `bio-box nbr` guilt-by-association recall at `0.1%` FPR on novel-fold dark systems (`Foldseek E > 1`), with `< 15%` relative recall drop on `1.5–5 kb` truncated contigs and `< 3%` AP drop under `NR90` representative substitution?

---

## 2. Manifest-Backed Data Access & Invariants (`AGENTS.md` §3–§4)

Every data read in `w38a_plm_bakeoff.py` obeys `bio-box`'s hard data-layer rules:
- **Manifest-only reads (`AGENTS.md` §3.9):** All Parquet and BGZF shard paths are resolved through `biobox.catalog.ReleaseFiles.of(root, release)` and `biobox.l0seq.store.SeqStore(root, files)` — **zero filesystem globs** across stages or releases.
- **Stable keys (`AGENTS.md` §3.2–§3.4):** Every record carries its 16-byte `feature_uid` and `protein_uid` (hex-encoded in JSONL/TSV). Tape `ordinal` is never used as a key without `feature_uid` beside it; release-local `token_id`s never cross the benchmark boundary.
- **Bounded memory (`AGENTS.md` §4.1):** Release tables (`protein_seq`, `contig_tapes`, `l1_domains`, `protein_occurrence`) are streamed one row group at a time; background sampling uses order-invariant keyed bottom-`k` reservoir sampling (`blake3(seed || protein_uid)`) so peak heap is `O(k)` regardless of corpus size.
- **Determinism (`AGENTS.md` §4.3):** Every split, background sample, and truncation window is keyed by content hashes (`blake3`) + explicit `--seed` (`20261125` default), yielding byte-identical outputs independent of thread count or Parquet read order.

---

## 3. Strict Holdout Construction & Three Hard Stratifications

### 3.1 Positive Set & Dual-Gate Leakage Filter (Closing the `30%` Gap)

1. **Candidate positives:** Map curated anti-phage defense proteins from:
   - `DefenseFinder` (`>= 152` canonical systems),
   - Millman et al. 2022 (`Cell Host & Microbe`), and
   - Laub lab experimentally validated `DS-1..DS-42` (`DeWeirdt, Mahoney & Laub 2026`, *Science* `PMC13092281` Table S7)
   onto the live release (`--release r26.11.2` or `r26.11.3`) at `>= 90%` AA identity and `>= 80%` alignment coverage.
2. **Mandatory dual-gate leakage filter (`build-holdout`):** Unlike `DefensePredictor`'s CV split (`>30%` AA ID and `qcov > 0.8`, `05_04_02_split_data.ipynb`), a held-out test protein is admitted **only if both hold simultaneously**:
   - **Sequence/profile identity gate:** `max_train_aa_id < 0.20` (`< 20%` AA identity against the entire training/query split under `mmseqs search -s 7.5 --num-iterations 3`), **AND**
   - **HMM sensitivity gate:** `min_train_hmm_evalue > 10.0` under `hmmsearch` at fixed `Z = 1e8`, `domZ = 1e8` (`SD-04`) against all training-split HMMs.
   - **Fail-closed rule:** Missing or `NaN` homology fields are rejected; proteins with zero hits above reporting threshold record `max_train_aa_id = 0.0` and `min_train_hmm_evalue = inf`.
3. **Structural novelty sub-stratification:**
   - **`remote_homolog`:** Passes the dual `<20%` ID & `E > 10` gate and has `Foldseek E <= 1.0` against training structures (tests Gate 1 & Gate 2).
   - **`novel_fold`:** Passes the dual gate **and** has `Foldseek E > 1.0` against all training structures (tests Gate 3).

### 3.2 Stratification (i): Multi-Domain Fusions (`>= 600` aa & `> 1,022` aa)

- Held-out positives are partitioned into:
  - **`monodomain`**: `< 600` aa or single-domain architecture;
  - **`fusion`**: `>= 600` aa with `>= 2` distinct non-overlapping domain envelopes (or `>= 1` defense effector fused to `>= 1` non-defense/promiscuous partner, e.g., `NACHT+effector`, `TniQ+nuclease`, `helicase+PD-(D/E)XK`);
  - **`fusion_gt_1022`**: subset of `fusion` with `aa_len > 1,022` aa, specifically isolating `ESM2`'s `1,022`-aa hard truncation (`defense_predictor/core.py:101,125-128`).
- Exposes global-pooling dilution (`Option B`) and hard-truncation blindness (`DefensePredictor`), while verifying local alignment resilience (`3Di` Foldseek, `pLM-BLAST`).

### 3.3 Stratification (ii): Closed `P1..P6` Operons vs Empirical `1.5–5 kb` `SPIRE` Truncation

- Take every held-out operon residing on closed `P1..P6` (`S2C`) chromosomes (`is_circular = True` or complete linear chromosome) and generate a paired artificial-fragmentation set via `w38a_plm_bakeoff.py truncate-loci`:
  - Draw a deterministic fragment length `L ∈ [1500, 5000]` bp from `SPIRE`'s empirical contig-length distribution (`C-42`; or the built-in right-skewed empirical CDF when `--spire-lengths` is omitted), keyed by `blake3(seed || contig_uid || anchor_feature_uid)`.
  - **Circular-aware slicing (`AGENTS.md` §3.1):** If the parent contig is circular and the operon wraps the origin (either a slot has `end > contig_len` or the shortest arc spanning the operon crosses `0`), unwrap slots modulo `contig_len` before slicing, then project intersecting slots into linear fragment coordinates `[0, L)`.
  - Mark every fragment `is_circular = False`, recompute `dist_to_contig_end_bp = min(local_start, L - local_end)`, `orfs_to_contig_start`, `orfs_to_contig_end`, and `is_partial` (True when clipped at `0` or `L`), while retaining exact `parent_contig_uid`, `parent_feature_uid`, and `protein_uid`.
- Measuring `truncation_relative_recall_drop = (recall_closed - recall_truncated) / recall_closed` directly quantifies degradation when `±2` / `±10` neighbors are censored and genome-wide GC z-scoring is unavailable.

### 3.4 Stratification (iii): `1,000,000` Real Background Genes (FP Trap Cohorts)

- Draw `1,000,000` real non-defense genes from the live release (excluding all known `DefenseFinder` / Millman / Laub positives and any gene with `>= 30%` AA ID to them) using deterministic hash-keyed reservoir sampling across five stratified cohorts:
  1. **`ntpase`** (`150,000` genes): housekeeping P-loop NTPases, AAA+ ATPases, and DEAD/DEAH helicases (e.g., `PF00004`, `PF00005`, `PF00270`, `PF00271`, `PF05729`);
  2. **`methyltransferase`** (`150,000` genes): housekeeping Rossmann-fold / SAM-dependent methyltransferases (`PF01555`, `PF08241`, `PF08242`, `PF13649`, `PF13847`);
  3. **`transposase`** (`150,000` genes): solitary IS elements, DDE transposases, and integrases (`PF00665`, `PF01527`, `PF01609`, `PF02371`, `PF13358`);
  4. **`hydrolase`** (`100,000` genes): metabolic alpha/beta hydrolases, phosphatases, and non-defensive nucleases (`PF00561`, `PF00702`, `PF12697`);
  5. **`general`** (`450,000` genes): unstratified random background genes across the release.
- Reports both overall `FPR <= 1%` / `0.1%` and per-cohort trap FPRs so promiscuous catalytic superfamilies cannot hide behind balanced-dataset AUROC.

---

## 4. Competitors & Mandatory Bug Fixes for `DefensePredictor`

| Group | Method Key | Configuration | Notes |
|---|---|---|---|
| **CPU Baselines** | `jackhmmer_3iter` | 3 rounds, `Z = 1e8`, `domZ = 1e8`, `incE = 1e-3` | Gate 1 baseline |
| **CPU Baselines** | `hhblits_3iter` | 3 iterations vs Pfam + PDB70 (`-e 1e-3`) | Gate 2 baseline |
| **CPU Baselines** | `nbr_gba` | `bio-box nbr` guilt-by-association (`F-01`/`F-03`, `±10` slots) | Gate 3 baseline |
| **Option A (`SD-21`)** | `prostt5_foldseek` | Foldseek v10 built-in `--prostt5-model` + `Foldseek` (`--sort-by-structure-bits 0`) | Preferred if within `2` pp of `ESM-2 3B 3Di` |
| **Option A (`SD-21`)** | `esm2_3b_3di_foldseek` | Domainator `ESM-2 3B 3Di` (`esmologs`, `<= 2,500` aa) + `Foldseek` / `reseek` | Compare recall & throughput vs `ProstT5` |
| **Option B (ANN)** | `esmc_300m_full_fp16` | `ESM-C 300M` (`960`d, `bf16` -> `fp16`) full-ORF mean-pool cosine | Tests multi-domain dilution |
| **Option B (ANN)** | `esmc_300m_segments_ivfpq` | `ESM-C 300M` over `segments.faa` (`SD-06`) + `FAISS IVF-PQ 64B` | Tests whether segment pooling saves Option B |
| **Option C/D** | `defense_predictor_stock` | Unmodified `PeterDeWeirdt/defense_predictor` (`ESM2-150M`, `±2` CDS, 5-fold LGBM) | Reference baseline showing unpatched fragmentation drop |
| **Option C/D** | `defense_predictor_fixed` | Patched for the two §5 audit bugs (below) | Gate 3 candidate |
| **Option C/D** | `esm_df_geneclr_df` | `ESM_DF` (rsLoRA) + `GeneCLR_DF` (12L DS-Attn Transformer, `asinh(0.03*d)`), edge-censored | Gate 3 candidate |
| **Option D (Refinement)** | `plm_blast_top1k` | Per-residue `ProtT5` + `pLM-BLAST` local SW on top-1,000 candidates per query | Candidate-scale refinement only |

**Mandatory `defense_predictor_fixed` Patches (`pLM_remote_system_search.md` §5):**
1. **Fix 1-CDS `NaN -> +1.1 SD` bug (`core.py:158,166,451`):** Never z-score GC/dinucleotides using single-contig `std(ddof=1)` on contigs with `< 5` CDSs; instead center/scale against the contig's `contig_dim` genus/species median & IQR (or pass `0.0` / explicit missing mask when `< 10 kb` unbinned).
2. **Fix circular origin wraparound (`core.py:70,314-320`):** When `is_circular = True`, wrap `±2` neighbor indices modulo `n_slots` across `(0, 1) <-> (N-2, N-1)` and compute inter-gene distances using `biobox.l3.geometry.pair_geometry` (`end > contig_len`, `AGENTS.md` §3.1) rather than linear slicing + `end < start`.

---

## 5. Exact Quantitative Go/No-Go Gates (Gates 1–3)

Evaluated automatically by `w38a_plm_bakeoff.py check-gates`:

| Gate | Primary Metric(s) | Exact Pass Criterion | Action on Pass | Action on Fail |
|---|---|---|---|---|
| **Gate 1 (`Option A` `3Di`)** | (1) `Recall@1%FPR` on `remote_homolog` + `fusion` vs `jackhmmer_3iter`; (2) query latency ratio vs `foldseek_afdb` | `delta_recall_pp >= +20.0` **AND** `latency_ratio <= 1.50` | Proceed to `W-38b`/`W-38c` (`3di_nr50/`); select `prostt5_foldseek` if within `2.0` pp of `esm2_3b_3di_foldseek`, else `esm2_3b_3di_foldseek`. | Hold corpus-wide `3Di`; keep `HHblits` + external Foldseek only. |
| **Gate 2 (`Option B` Dense ANN)** | (1) Marginal unique TP systems at `FPR <= 0.1%` beyond `union(jackhmmer, hhblits, best_3di)`; (2) `fusion` vs `monodomain` recall loss | `unique_tp_gain >= +0.10` (`+10%`) **AND** `fusion_recall_loss < 0.05` (`5%`) | Approve `segments.faa` `IVF-PQ` pack at `Tier 1/2` only. | **Reject Option B permanently**; zero corpus-wide dense ANN index. |
| **Gate 3 (`Option C` vs `Option D`)** | (a) `Recall@0.1%FPR` on `novel_fold` (`Foldseek E > 1`) vs `nbr_gba`; (b) Stratification (ii) relative recall drop on `1.5–5 kb` contigs; (c) `NR90` rep substitution AP drop | (a) `ratio >= 2.0` **AND** (b) `truncation_drop < 0.15` **AND** (c) `nr90_ap_drop < 0.03` | Approve `W-38e` ephemeral `contig_tapes` scorer (`>= 5` genes, `0%` NR90 error via VRAM streaming) + `W-38d` on-demand rounds. | Keep **`W-38d` (`Option D` on-demand rounds) only**; no corpus-wide operon scorer. |

---

## 6. `wilkcomp03` (`4× RTX 6000 Ada`, `192 GB` VRAM) Runtime & VRAM Budget

### 6.1 Hardware & Throughput Constants
- **GPU hardware (`*measured*`, `wilkcomp03`, `2026-11-19`, `nvidia-smi -L`):** `4× NVIDIA RTX 6000 Ada Generation`, `48 GB` GDDR6 ECC VRAM each (`192 GB` total), FP16/BF16 tensor throughput `~1.7×` A100-40GB (`*estimate*`, `pLM_remote_system_search.md` §5.4).
- **Benchmark workload size (`*derived*`):** `1,000,000` background genes + `~5,000` held-out positives + `~20,000` `±2` flanking context genes + `~5,000` truncated-locus genes = **`1.03 × 10^6` total ORFs** (`~1.05M` with `segments.faa` `1.5×` factor on `esmc_300m_segments_ivfpq` -> `~1.55 × 10^6` segments). Mean length `324` aa (`*measured*`, `P2` mean `323.6` aa, `C-38` N16) -> **`3.34 × 10^8` aa** (`5.02 × 10^8` aa for `segments.faa`).
- **Derived single-RTX-6000-Ada throughputs (`1.7×` A100-40GB constants from `pLM_remote_system_search.md` §3.3):**
  - `ProstT5` (Foldseek v10, `bf16`): `77 × 1.7 = 131` prot/s (`42.4k` res/s, `*derived*`).
  - `ESM-2 3B 3Di` (`bf16`): `~25 × 1.7 = 42.5` prot/s (`13.8k` res/s, `*derived*` from `~5k–10k` res/s A100).
  - `ESM-C 300M` (`bf16` FlashAttention2): `247 × 1.7 = 420` prot/s (`136k` res/s, `*derived*`).
  - `ESM2-150M` (`DefensePredictor`, `bf16`): `247 × 1.7 = 420` prot/s (`*derived*`).
  - `ESM2-35M` (`GeneCLR_DF`, `bf16`): `350 × 1.7 = 595` prot/s (`*derived*`).

### 6.2 Per-Step Runtime & VRAM Table (`*derived*`)

| Step / Competitor | Sequences (`N`) | Peak VRAM / GPU (`*estimate*`) | 1× Ada Wall-Time (`*derived*`) | 4× Ada Wall-Time (`*derived*`, `0.90` scaling) | Disk Footprint (`*derived*`) |
|---|---:|---:|---:|---:|---:|
| `ProstT5` `3Di` (`prostt5_foldseek`) | `1.03M` | `14 GB` | `1.03e6 / 131 / 3600` = **`2.18 h`** | **`0.61 h`** | `0.12 GB` zstd (`120 B/ORF`) |
| `ESM-2 3B 3Di` (`esm2_3b_3di_foldseek`) | `1.03M` | `26 GB` | `1.03e6 / 42.5 / 3600` = **`6.73 h`** | **`1.87 h`** | `0.12 GB` zstd |
| `ESM-C 300M` full-ORF + `segments.faa` | `1.03M` + `1.55M` | `12 GB` | `2.58e6 / 420 / 3600` = **`1.71 h`** | **`0.47 h`** | `4.95 GB` `fp16` / `0.17 GB` `IVF-PQ` |
| `DefensePredictor` (`ESM2-150M` + 5-fold LGBM) | `1.03M` | `10 GB` | `1.03e6 / 420 / 3600` = **`0.68 h`** + `0.25 h` LGBM = **`0.93 h`** | **`0.26 h`** | `1.32 GB` `fp16` (ephemeral) |
| `ESM_DF + GeneCLR_DF` (`ESM2-35M` + 12L DS-Attn) | `1.03M` | `10 GB` | `1.03e6 / 595 / 3600` = **`0.48 h`** + `0.20 h` attn = **`0.68 h`** | **`0.19 h`** | `0.99 GB` `fp16` (ephemeral) |
| `pLM-BLAST` (`ProtT5` per-residue on top-1k/query) | `~50k` unique | `16 GB` | `5.0e4 / 131 / 3600` = **`0.11 h`** + `0.30 h` SW = **`0.41 h`** | **`0.15 h`** | `~33 GB` ephemeral `/tmp` (deleted after scoring) |
| **Total GPU Bake-Off (`W-38a`)** | — | **`<= 26 GB / 48 GB`** | **`12.64` single-GPU Ada-h** (`~5.91 h` if `ESM-2 3B 3Di` runs on `200k` stratified subset per `SD-21`) | **`~3.55 h` wall** (`~1.7 h` on `200k` `3B` subset) | **`< 40 GB` peak `/tmp`** |

**Operational Guardrails (`AGENTS.md` §7):**
- Confirm `wilkcomp03` GPU etiquette with pod owner before launch (`DEC-PLM-1`); check `nvidia-smi` for idle GPUs.
- Run detached via `systemd-run --user --unit=w38a-bakeoff --collect -p MemoryMax=180G` at `nice 19`, checking for `<root>/STOP` before each step.
- Work strictly under `/tmp/if2453/bio-box/w38a/`; never write to `/stor` except the final measurement note (`planning/measurements/<date>_plm_bakeoff_<release>.md`).

---

## 7. Separate CUDA Container Draft (`container/environment-cuda.yml` Outline)

Per `DEC-PLM-2`, `bio-box.sif` (`container/environment.yml`) remains **CPU-only** and untouched. Upon `W-38b` approval, a separate `bio-box-cuda.sif` is built from the following pinned outline (`container/environment-cuda.yml` — **draft only; not yet created**):

```yaml
# Draft specification for bio-box-cuda.sif (DEC-PLM-2; separate from CPU container/environment.yml).
# Do NOT merge into container/environment.yml. Pin exact build hashes at W-38b approval.
channels:
  - pytorch
  - nvidia
  - conda-forge
  - bioconda
dependencies:
  - python=3.12
  - cuda-version=12.4
  - pytorch=2.5.*
  - pytorch-cuda=12.4
  - foldseek>=10.941cd33      # built-in ProstT5 3Di translator (--prostt5-model) + GPU search
  - mmseqs2=18.8cc5c         # matches container/environment.yml
  - hmmer=3.4                # matches container/environment.yml
  - hhsuite=3.3.0            # HHblits 3-iter baseline
  - lightgbm=4.5.*           # DefensePredictor 5-fold ensemble
  - faiss-gpu=1.9.*          # Option B IVF-PQ 64B evaluation
  - pyarrow>=17
  - duckdb>=1.1
  - numpy>=2
  - blake3>=1
  # Pip/uv layer inside bio-box-cuda.sif:
  #   esm (EvolutionaryScale ESM-C 300M), transformers + sentencepiece (ProtT5 / pLM-BLAST),
  #   fair-esm (DefensePredictor ESM2-150M), peft (ESM_DF rsLoRA), flash-attn (Ada sm_89)
```

---

## 8. Reproducible Execution Commands (`scripts/stress/w38a_plm_bakeoff.py`)

Once PI approval is granted and `r26.11.2` / `r26.11.3` is live on disk:

```bash
# 1. Build deterministic strict holdout (<20% AA ID AND hmmsearch E > 10) + 1M trap-enriched background
python scripts/stress/w38a_plm_bakeoff.py build-holdout \
  --root /tmp/if2453/r10/root \
  --release r26.11.2 \
  --positives /tmp/if2453/bio-box/w38a/curated_positives.jsonl \
  --out /tmp/if2453/bio-box/w38a/holdout \
  --bg-size 1000000 \
  --seed 20261125

# 2. Generate paired 1.5-5 kb SPIRE-truncated loci (Stratification ii) preserving circular wraparound & feature_uid
python scripts/stress/w38a_plm_bakeoff.py truncate-loci \
  --loci /tmp/if2453/bio-box/w38a/holdout/closed_loci.jsonl \
  --out /tmp/if2453/bio-box/w38a/holdout/truncated_loci.jsonl \
  --seed 20261125 \
  --min-bp 1500 \
  --max-bp 5000

# 3. After running CPU baselines + GPU competitors, evaluate stratified PR-AUC, Recall@1%FPR, Recall@0.1%FPR, Recall@100
python scripts/stress/w38a_plm_bakeoff.py evaluate \
  --holdout /tmp/if2453/bio-box/w38a/holdout/holdout.jsonl \
  --predictions /tmp/if2453/bio-box/w38a/scores/manifest.json \
  --out /tmp/if2453/bio-box/w38a/eval_report.json

# 4. Evaluate falsifiable Go/No-Go Gates 1-3
python scripts/stress/w38a_plm_bakeoff.py check-gates \
  --eval-json /tmp/if2453/bio-box/w38a/eval_report.json \
  --out /tmp/if2453/bio-box/w38a/gates_report.json
```
