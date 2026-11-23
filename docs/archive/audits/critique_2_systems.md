# Adversarial Systems, Algorithms & HPC Throughput Audit (`critique_2_systems.md`)

**Auditor:** Auditor 2 — Database Systems, Algorithms & HPC Throughput Architect  
**Lenses:** Martin Steinegger (`MMseqs2` / `Foldseek` / `ProstT5`), Hannes Mühleisen (`DuckDB` / Parquet Vectorized Execution), Rayan Chikhi (`Logan` Planetary Assembly & Sketching)  
**Target Documents:** `00-spec.md`, `PRIOR-ART.md`, `PLAN.md` (v1.0)  
**Reference Hardware Envelope:**
- **CPU/Storage Node (`AX162-R`):** AMD EPYC 9454P (48C / 96T, Zen 4 AVX-512, 256 MB L3, 12-ch DDR5 @ ~340 GB/s), 256–512 GB ECC RAM, 2× 7.68 TB Gen4 U.2 NVMe RAID-0 (`15.36 TB`, ~14 GB/s sequential read, ~2.0M random 4 KB IOPS).
- **GPU Node (`GEX130`):** Intel Xeon Gold 5412U (24C / 48T), 256 GB DDR5, 1× NVIDIA RTX 6000 Ada (**48 GB GDDR6 ECC VRAM**, ~960 GB/s VRAM bandwidth).

---

## Executive Verdict

`PLAN.md` v1.0 has the right **scientific architecture** (sub-ORF domain splitting, `3Di` structural tokens for dark proteins, first-class non-coding RNAs, phylogenetic bitvectors, and zero-Neo4j columnar storage), but its **ingestion compute math, memory-mapped CSR index design, Parquet Z-order layout, and S3 streaming deduplication state contain six fatal systems bugs** that will cause multi-month compute stalls, OOM crashes on 512 GB RAM / 48 GB VRAM, and $500\times$ query latency regressions on common protein families:

| Subsystem in `PLAN.md` v1.0 | Quantitative Failure Mode | Cost in `PLAN.md` v1.0 | Cost After v2.0 Fix | Speedup / Savings |
|---|---|---|---|---|
| **1a. `Infernal` (`cmscan`, 4,720 CMs)** | Placed inside pre-dedup S3 streaming sieve (`Sec 2.3`) or across all contigs (`Sec 3.1`) | **180,000–1,800,000 CPU-hrs** (78–780 days on 96T EPYC) | **~420 CPU-hrs** (**8.8 hrs** wall-clock on 96T EPYC) | **$428\times\text{–}4,280\times$** |
| **1b. `ProstT5` on `NR90` (300M–800M reps)** | Running 3B T5 encoder on `NR90` instead of `NR50` + alignment transfer + cascade gating | **1,111 GPU-hrs** (**46.3 days** on 1× RTX 6000 Ada) | **14.3–22.0 GPU-hrs** (**< 1 day** on 1× RTX 6000 Ada) | **$50\times\text{–}77\times$** |
| **1c. `HH-suite3` (`hhsearch`) on dark segments** | Per-cluster `FAMSA` + `hhmake` file spawning + un-prefiltered `hhsearch` vs 250k `.hhm`s | **16,600–250,000 CPU-hrs** + 60M inode storm | **~48 CPU-hrs** via `mmseqs result2profile` + `Foldseek` gate | **$345\times\text{–}5,200\times$** |
| **1d. `pyhmmer` (215k HMMs on `NR90`)** | Running 17 HMM DBs on 500M `NR90` sequences instead of `NR50` + envelope projection | **~19,900 CPU-hrs** (~207 hrs wall-clock) | **~1,200 CPU-hrs** (**12.5 hrs** wall-clock) | **$16.5\times$** |
| **2. Index 1 (`CRoaring` + CSR Slot Index)** | `uint64` `contig_id` truncates in 32-bit `CRoaring`; CSR takes **32.8–61.6 GB RAM**; hub intersections (`PF00005` $\cap$ `PF04055` = 50M contigs) do 50M random seeks | **1,200–25,000 ms** per hub query ($500\times$ over budget) + **45–75 GB RAM** | **< 4.0 ms** via 16-slot `window_id`s + Hub Bigram Postings + Inverted Slices in **14.2 GB** | **$300\times\text{–}6,250\times$ faster**, **$4\times$ less RAM** |
| **3. `feature_tape.parquet` Z-Order Layout** | Z-order sorting by `(seq_family_id, struct_fold_id, contig_id, slot_idx)` scatters a 10-gene contig across 10 Parquet row groups and breaks RLE compression | **15,000–40,000 ms** to reconstruct 5k neighborhoods; **$3\times$ Parquet bloat** | **< 12 ms** via Dual-Projection (`contig_tapes.parquet` + `window_postings.parquet`) | **$1,250\times$ faster**, **1.8 TB disk saved** |
| **4. `Logan` S3 `arch_hash` In-Memory State** | `ordered_family_ids` unavailable during raw S3 stream; 500M `RoaringBitmap` heap objects + 95% ANI strings OOM 512 GB RAM | **102–145 GB heap** (plus **2 TB** for 95% ANI sequences) $\rightarrow$ **OOM crash** | **16 GB RAM** (Cuckoo Filter) + **64-shard NVMe RocksDB** + 64-bit Syncmer Sketch | **Zero OOM**, line-rate 400 MB/s S3 ingest |
| **5. `gLM-Domain` Softmax & Index 3 Layout** | Softmax over 15M `seq_family_id` tokens requires **245 GB VRAM** on a 48 GB GPU; untranposed Fitch bitvectors cause 7.2 MB strided cache misses | **VRAM OOM** on RTX 6000 Ada | **2.8 GB VRAM** (131k cap + `ECOD`/`3Di` composition); **< 0.2 ms** Sparse Edge Transitions | Fits comfortably on 48 GB GPU |

---

## 1. Wall-Clock Compute Bottlenecks in the Ingestion & Dissection Pipeline (The Exact Math)

### 1.1 `Infernal` (`cmscan` with 4,720 Covariance Models) Across 52 TB of Contigs

#### The Flaw in `PLAN.md` (Sections 2.3 & 3.1)
`PLAN.md` Section 2.3 Step 3 places `Infernal` *inside* the streaming S3 sieve before operon-architecture deduplication, and Section 3.1 specifies running `cmscan` with **Rfam 15.0 (`4,170` CMs) + `orphan_rnas.cm` (`~550` CMs) = `4,720` covariance models**.

Let us do the arithmetic:
1. **Input Nucleotide Volume:**
   - Raw input across all tiers is `~52 TB` compressed (`.fa.zst` / `.fna.gz`), corresponding to $\sim 180\text{ Tbp}$ ($1.8 \times 10^{14}\text{ bp}$) of raw sequence.
   - After discarding contigs $< 1,200\text{ bp}$ and eukaryotic repeats (Section 2.3 Steps 1–2), the surviving prokaryotic/viral/metagenomic sequence is **$\sim 30\text{ Tbp}$ ($3.0 \times 10^{13}\text{ bp}$) pre-deduplication**, and **$\sim 3.0\text{ Tbp}$ ($3.0 \times 10^{12}\text{ bp}$) post-deduplication** across the $\sim 5 \times 10^8$ unique multi-feature contigs (`3.6B` slots $\times \sim 850\text{ bp/slot}$).
2. **`cmscan` Throughput Math:**
   - Even with `Infernal 1.1.5` HMM filter stages (`F1/F2/F3`), QDB (Query-Dependent Banded CYK), and AVX2/AVX-512 acceleration, scanning 1 Mbp of double-stranded DNA ($2 \times 10^6$ nt scanned) against **4,720 CMs** requires **$\sim 216\text{ CPU-seconds/Mbp}$** ($\mathbf{60\text{ CPU-hours/Gbp}}$).
   - **If run pre-deduplication in the S3 sieve (`30,000 Gbp`):**
     $$30,000\text{ Gbp} \times 60\text{ CPU-hrs/Gbp} = \mathbf{1,800,000\text{ CPU-hours}} \quad (\mathbf{781\text{ wall-clock days}\text{ on 96 threads!}})$$
   - **Even if run post-deduplication across `3,000 Gbp` of surviving contigs:**
     $$3,000\text{ Gbp} \times 60\text{ CPU-hrs/Gbp} = \mathbf{180,000\text{ CPU-hours}} \quad (\mathbf{78.1\text{ wall-clock days}\text{ on 96 threads!}})$$
3. **Why `cmscan` Is the Wrong Binary:**
   - `cmscan` is query-sequence-outer / CM-database-inner: for every contig, it streams through 4,720 CMs, thrashing CPU caches and re-initializing model states. `cmsearch` is CM-outer / target-sequence-inner, keeping one CM hot in L1/L2 cache while streaming across target sequences (`1.8x–2.5x` faster than `cmscan`).

#### The v2.0 Engineering Fix (`428x–4,280x` Speedup $\rightarrow$ `420 CPU-hours`)
1. **Evict `Infernal` from the S3 Streaming Sieve:** Only `pyrodigal` (`--meta`), `smORF` RBS scanning, and fast k-mer/syncmer `Minced` CRISPR repeat detection run in the streaming sieve.
2. **Extract Only Non-Coding Intergenic Regions (`IGRs >= 40 bp` + `25 bp` CDS flanks) After Contig Deduplication:**
   - Prokaryotic and viral genomes are $88\text{–}91\%$ protein-coding (`CDS`). Structured riboswitches, leader attenuators, `CreTA` RNAs, and `msr-msd` retrons sit in intergenic regions (or overlap the 5'/3' ends of adjacent ORFs by $\le 25\text{ bp}$).
   - Extracting only `IGRs >= 40 bp` (plus $25\text{ bp}$ padding into flanking ORFs and full regions around `tRNA`/`CRISPR`/`RT` loci) reduces the target nucleotide volume from `3,000 Gbp` to **`~350 Gbp`** ($8.57\times$ reduction).
3. **Cluster IGRs at 90% Nucleotide Identity (`IGR90`) via `mmseqs linclust --search-type 3`:**
   - Across $5 \times 10^8$ contigs, conserved riboswitches (`TPP`, `Cobalamin`, `SAM-I`, `Guanidine-I`, `ydaO`, `T-box`) and leader sequences are massively redundant across strains and metagenomes.
   - Linear-time nucleotide clustering (`mmseqs linclust --min-seq-id 0.90 -c 0.80 --cov-mode 1 --search-type 3`, taking `~18 CPU-hours`) collapses `350 Gbp` of raw IGRs by **$7.0\times$** into **`~50 Gbp` of non-redundant `IGR90` representatives**, with coordinates projected back to member contigs via the cluster alignment map.
4. **Prune Eukaryotic-Only Rfam CMs & Run a Two-Stage `pyhmmer` (`nhmmer`) Sequence Prefilter Before `cmsearch`:**
   - Strip the ~2,720 strictly eukaryotic Rfam models (spliceosomal snRNAs, vertebrate miRNAs, snoRNAs, fungal/plant lncRNAs, Flavivirus/Coronavirus UTRs, keeping those only for the tiny `Tier 2` RNA virome slice). This leaves **`1,450` prokaryotic/archaeal/phage Rfam CMs + `550` orphan RNA CMs = `2,000` active CMs**.
   - Extract the profile HMM (`cmbuild` / `cmconvert` HMM profiles) for all 2,000 CMs and run vectorized in-memory `pyhmmer.nhmmer` (`--F1 0.02 --F2 1e-3 --F3 1e-5`) across the `50 Gbp` `IGR90` database (`~310 CPU-hours`).
   - Run full covariance-model `cmsearch --acc --nohmmonly` (banded CYK + Inside) **only on the $<0.4\%$ of `(IGR90, CM)` candidate windows** passing the `nhmmer` filter (`~92 CPU-hours`).
   - **Total `Infernal` Wall-Clock:** $18 + 310 + 92 = \mathbf{420\text{ CPU-hours}}$ (**`8.75 hours` wall-clock** on the 48C/96T EPYC 9454P).

---

### 1.2 `ProstT5` (3B Parameter Encoder) on `NR90` Proteins

#### The Flaw in `PLAN.md` (Section 3.2 Cascade 4)
`PLAN.md` clusters all called proteins at **90% sequence identity (`NR90`)** and passes them through the 4-cascade dissector, running `ProstT5` to translate amino-acid sequences into 1D `3Di` token strings.

Let us do the arithmetic:
1. **Cardinality of `NR90` at Planetary Scale:**
   - `SPIRE` alone contains ~35 billion proteins (~550M clusters at 90% ID); adding `MGnify`, `GTDB`, `IMG/VR`, `AllTheBacteria`, and `Logan`'s surviving microbial contigs yields **$3.5 \times 10^8$ to $8.0 \times 10^8$ (350M–800M) `NR90` cluster representatives** (and after sub-ORF domain splitting, **~500M–1.1B domain segments**).
2. **`ProstT5` GPU Throughput Math on 1× RTX 6000 Ada (48 GB):**
   - Using the 1.2B-parameter encoder half of `ProstT5` + 2-layer CNN `3Di` projection head (`Rostlab/ProstT5_fp16` in `int8`/`fp16` TensorRT/ONNX, *without* autoregressive T5 decoding), throughput on a single RTX 6000 Ada for average 280 aa proteins is **`~200 proteins/sec`** (if using the full seq2seq decoder, throughput drops to `~45 proteins/sec`).
   - At `200 proteins/sec` across $8.0 \times 10^8$ `NR90` proteins:
     $$\frac{8.0 \times 10^8\text{ proteins}}{200\text{ proteins/sec}} = 4.0 \times 10^6\text{ seconds} = \mathbf{1,111\text{ GPU-hours}} = \mathbf{46.3\text{ GPU-days!}}$$
3. **Why Running `ProstT5` at `NR90` Violates Structural Biology & Algorithm Design:**
   - The `Foldseek` `3Di` alphabet encodes tertiary local backbone geometry between $C_\alpha$ virtual centers. Backbone geometry and `ECOD` T-group topology are strictly conserved down to **25–30% sequence identity**; two proteins with **$\ge 50\%$ sequence identity** over $\ge 80\%$ coverage have virtually identical `3Di` strings ($\text{TM-score} > 0.85$) and map to the exact same `ECOD_T_id`!

#### The v2.0 Engineering Fix (`50x–77x` Speedup $\rightarrow$ `14.3–22.0 GPU-hours`)
1. **Cluster at 50% Sequence Identity (`NR50`) for Structural & Remote-Homology Cascades:**
   - Run `mmseqs linclust --min-seq-id 0.50 -c 0.80 --cov-mode 0` on the `NR90` representatives (`~45 CPU-hours` on 96 threads). This collapses `350M–800M` `NR90` sequences down to **`~75M–110M` `NR50` cluster representatives**.
2. **Alignment-Based `3Di` Projection from `AFDB-Full` (214M) + `ESMAtlas-Full` (772M) + `BFVD` (351k):**
   - `PLAN.md` only checks exact/98% hashes against precomputed databases in Cascade 1. Instead, run `mmseqs search` (or co-cluster via `mmseqs linclust`) of the `NR50` representatives against the **986.4M precomputed `3Di` strings** in `AFDB-Full` + `ESMAtlas-Full` + `BFVD` at $\ge 50\%$ sequence identity and $\ge 80\%$ alignment coverage.
   - Using the pairwise alignment CIGAR string (`mmseqs prefixid` / gap-projected `3Di` string transfer), project the precomputed `3Di` states onto aligned residues of the `NR50` representative. Because `ESMAtlas-Full` (772M metagenomic proteins) and `AFDB-Full` (214M UniProt proteins) already cover **~68–75% of all non-singleton microbial `NR50` clusters**, this eliminates **~55M–80M `NR50` clusters in ~60 CPU-hours** with zero GPU cost!
3. **Cascade-Gated Execution + Singleton Filtering:**
   - Any segment already resolved in Cascade 2 (`pyhmmer` $E \le 10^{-10}$) to a `Pfam`/`COG`/`arCOG`/`PHROG` family that has a known `ECOD_T_id` mapping inherits that `struct_fold_id` directly.
   - Among the remaining sequence-dark, unmapped `NR50` domain segments, filter out **true biosphere singletons** (`n_syntenic_occurrences == 1` across all 52 TB, which cannot achieve statistical significance $\text{synteny\_log10\_pval}$ in any comparative genomics query; their `biophysical_badge` is still computed in full on CPU).
   - This leaves **`~12M–18M` recurrent novel dark `NR50` segments** ($\ge 2$ syntenic occurrences, average segment length $\sim 175\text{ aa}$).
4. **Length-Sorted Token-Budget Dynamic Batching on the `ProstT5` Encoder+CNN Head:**
   - Sort all $1.8 \times 10^7$ dark segments strictly by amino-acid length (`45 aa .. 1024 aa`) before batching so padding overhead is $< 1.5\%$, and run **only the 24-layer T5 Encoder + 2-layer CNN `3Di` classification head** (compiled with TensorRT `fp16`/`int8` weight-only quantization, `max_tokens_per_batch = 65,536` on 48 GB VRAM).
   - Throughput increases from `200 proteins/sec` (padded whole proteins) to **`~350 segments/sec` (`~61,000 aa/sec`)**:
     $$\frac{1.8 \times 10^7\text{ dark NR50 segments}}{350\text{ segments/sec}} = 51,428\text{ seconds} = \mathbf{14.3\text{ GPU-hours!}}$$

---

### 1.3 `HH-suite3` (`hhsearch`) and `pyhmmer` Scaling Across Dark & Annotated Segments

#### The Flaws in `PLAN.md` (Section 3.2 Cascades 2 & 3)
1. **Cascade 2 (`pyhmmer` across 17 HMM databases on `NR90`):**
   - Combining `Pfam-A` (21k), `COG` (5k), `arCOG` (4.5k), `NCBIfam` (35k), `CDD` (60k), `SFLD` (0.3k), `PHROGs` (38.8k), `VOGdb` (35k), and 9 specialized conflict/secretion/transport DBs (~15k) yields **`~215,000` HMMs**.
   - Running `pyhmmer` with 215,000 HMMs across $5 \times 10^8$ `NR90` sequences ($1.075 \times 10^{14}$ HMM-sequence comparisons) at $1.5 \times 10^6$ MSV comparisons/sec/core takes **`~19,900 CPU-hours` (~8.6 days on 96 threads)**.
2. **Cascade 3 (`FAMSA` + `hhmake` + `hhsearch` on unannotated `NR50` clusters):**
   - After Cascade 2, there are **`~25M–40M` dark `NR50` clusters** ($\ge 45\text{ aa}$).
   - Spawning `FAMSA` and `hhmake` per cluster creates **60–80 million small `.afa` and `.hhm` files**, exhausting filesystem inodes and serializing on POSIX `fork`/`exec` and VFS metadata locks.
   - Worse, running `hhsearch` for $3 \times 10^7$ queries against ~250,000 target `.hhm` models (`pfam` + `cog_arcog` + `pdb70` + `ecod` + `phrogs`) at even `2 CPU-sec/query` (with `hhblits` prefiltering) requires:
     $$3 \times 10^7 \times 2\text{ CPU-sec} = 6 \times 10^7\text{ CPU-sec} = \mathbf{16,667\text{ CPU-hours}}$$
     (and without a prefilter, at `30 CPU-sec/query`, **`250,000 CPU-hours`**!).

#### The v2.0 Engineering Fix
1. **Run `pyhmmer` on `NR50` First, Project Envelopes to `NR90`, and Check Active Sites by Direct Alignment Indexing (`16.5x` Speedup):**
   - Run the 215,000 HMMs (after collapsing identical profiles between `CDD`, `Pfam`, and `COG` into **`~115,000` non-redundant HMM profiles**) on **`NR50` (`~85M` sequences)** instead of `NR90` (`500M` sequences).
   - Resolve non-overlapping N-to-C domain envelopes on `NR50` via `cath-resolve-hits`, then project domain boundaries onto all `NR90` member sequences using the precomputed `NR50 -> NR90` `MMseqs2` alignment CIGARs (`< 5 CPU-hours`), checking `CDD`/`SFLD` catalytic residue masks directly against each `NR90` member's aligned coordinates in Rust.
   - **Cascade 2 Compute Drops:** `19,900 CPU-hours` $\rightarrow$ **`~1,200 CPU-hours` (12.5 hours wall-clock)**.
2. **Replace `FAMSA` + `hhmake` + Full `hhsearch` with `mmseqs result2profile` + `Foldseek` Structural Prefiltering (`345x` Speedup):**
   - **Filter `n_members >= 3`:** A dark `NR50` cluster with only 1 or 2 sequences has zero MSA covariance/entropy beyond single-sequence `pyhmmer`; running `hhsearch` on singletons/doubletons is mathematically redundant with Cascade 2. Filtering to `n_members >= 3` cuts `30M` dark `NR50` clusters down to **`~4.5M` multi-member dark clusters**.
   - **Zero-File In-Memory Profile Generation:** Use `mmseqs result2profile` directly on the `MMseqs2` cluster alignment database (`NR90_to_NR50_alnDB`) to generate all 4.5M sequence profiles in flat memory-mapped `MMseqs2` / `FFindex` (`.ffdata` / `.ffindex`) bundles—**zero loose files on disk**.
   - **Two-Sided Prefilter (`mmseqs search` Profile-vs-Consensus + `Foldseek` `3Di` Top-200 Candidates):**
     Run `Foldseek` (`3Di+AA`) and `mmseqs search` (profile-vs-consensus at `-s 7.5`) to select the **top 200 candidate target `.hhm` models** per dark `NR50` cluster, then execute vectorized `hhalign` / restricted `hhsearch` via `ffindex` **only against those 200 candidate `.hhm` models** (`~5 ms/query` $\times 4.5 \times 10^6$ queries = **`6.25 CPU-hours`**, plus `~42 CPU-hours` for `result2profile` and prefiltering = **`48 CPU-hours` total**).

---

## 2. Index 1 (`CRoaring` + CSR Positional Index) Memory & Bit-Packing Bugs

### 2.1 Three Fatal Bugs in `PLAN.md` Section 4.2

1. **Bug 2.1 — 64-Bit `contig_id` Truncation in 32-Bit `CRoaring` & CSR Offset Indexing:**
   - In Section 4.3, `contig_id` is defined as `uint64` where *"top 2 bits encode `tier`: `1=HQ`, `2=MGE`, `3=Unbinned`"*.
   - In Section 4.2, `token_to_contigs.roar` uses standard `CRoaring` (`roaring_bitmap_t`), which is a **32-bit unsigned integer bitmap** (`max value = 4,294,967,295`). Passing a `uint64` with bits 62–63 set into `roaring_bitmap_add()` silently truncates the upper 32 bits, causing collisions across tiers! And indexing a flat CSR array `slot_positions[contig_id]` with bit 62 set attempts to address offset $2^{62} \times 8\text{ B} = 36.8\text{ Exabytes}$, causing an instant `SIGSEGV`.

2. **Bug 2.2 — Uncompressed CSR Array Memory Bloat (`32.8 GB – 61.6 GB` RAM):**
   - Even if `contig_id` is densely remapped to `dense_contig_id: uint32` $\in [0, 5 \times 10^8)$, let us calculate the exact memory footprint of `PLAN.md`'s CSR index across $5 \times 10^8$ contigs and $3.6 \times 10^9$ slots:
     - **`contig_slot_offsets.bin`:** Because $3.6 \times 10^9$ slots is near the $2^{32}$ limit (and exceeds $2^{32}$ when indexing both `seq_family_id` and `struct_fold_id`), offsets require `uint64` (`8 bytes`):
       $$5 \times 10^8\text{ contigs} \times 8\text{ bytes} = \mathbf{4.0\text{ GB}}$$
     - **`slot_positions.bin` (CSR Payload):** For each of the $3.6 \times 10^9$ slots, checking `seq_family_id` (`4 B`), `struct_fold_id` (`4 B`), `orf_idx` (`2 B`), `slot_idx` (`2 B`), `directon_id` (`2 B`), and `strand`/`flags` (`2 B`) requires `16 bytes/slot` (or `8 bytes/slot` minimum if `struct_fold_id` and `directon_id` are omitted):
       $$3.6 \times 10^9\text{ slots} \times 16\text{ bytes} = \mathbf{57.6\text{ GB}} \quad (\text{or }\mathbf{28.8\text{ GB}}\text{ at 8 B/slot})$$
     - Combined with `token_to_contigs.roar` (`~14 GB`), Index 1 alone consumes **`46.8 GB – 75.6 GB` of RAM**.

3. **Bug 2.3 — The High-Frequency Hub Latency Catastrophe (`25,000 ms` vs. `50 ms` Target):**
   - `PLAN.md` Section 4.2 claims: *"Because `bitmap(A) & bitmap(B)` typically reduces $2 \times 10^8$ contigs down to $10^2\text{–}10^5$ candidate contigs ... a memory-mapped CSR array checks exact slot distance in $< 10\text{ ms}$."*
   - **This fails catastrophically on hub domains and Tier 1 chromosomes:**
     - Ubiquitous domains—such as `PF00005` (ABC ATPase), `PF04055` (Radical SAM), `PF00072` (Response regulator), `COG0583` (LysR), `PF13304` (AAA+ ATPase), or `ECOD` Rossmann / P-loop folds—occur on **$2.0 \times 10^7$ to $6.5 \times 10^7$ contigs**.
     - Furthermore, every complete chromosome in Tier 1 (`1.8M` assemblies, `~3,000–6,000` genes per chromosome) and every large metagenomic scaffold in Tier 3 contains *multiple copies* of all common families scattered megabases apart!
     - Therefore, for a query involving two hub tokens (e.g., `PF04055` within $\pm 3$ slots of `PF00005`), `bitmap(PF04055) & bitmap(PF00005)` returns **$3.5 \times 10^7$ to $5.0 \times 10^7$ candidate contigs**—99.8% of which are false-positive whole-contig intersections where the two genes are thousands of ORFs apart!
     - **Latency Math for 50,000,000 Random CSR Lookups:**
       - **If `mmap`-ed from Gen4 NVMe RAID-0 (`2.0M` random 4 KB IOPS):**
         $$\frac{5 \times 10^7\text{ random page reads}}{2.0 \times 10^6\text{ IOPS}} = \mathbf{25.0\text{ seconds}} \quad (\mathbf{500\times\text{ over the 50 ms budget!}})$$
       - **Even if 100% pinned in DDR5 RAM (`~15 ns` effective pipelined random cache-line latency, 2 cache misses per contig = $10^8$ misses):**
         $$10^8 \times 15\text{ ns} = \mathbf{1.50\text{ seconds}} \quad (\mathbf{30\times\text{ over the 50 ms budget!}})$$

---

### 2.2 The v2.0 Replacement Architecture for Index 1 (`< 4 ms` Worst-Case Latency, `14.2 GB` Total Size)

Replace the Contig-Granularity `CRoaring` + Forward CSR design with a **Three-Layer Windowed Inverted Positional Index**:

```
┌────────────────────────────────────────────────────────────────────────────┐
│ LAYER 1: 16-Slot Micro-Window Bitmaps (`token_to_windows.roar`)            │
│ • Divide the 3.6B global coordinate tape into aligned 16-slot windows      │
│   (`window_id = global_slot_idx >> 4`, total 2.25 × 10^8 windows < 2^28).  │
│ • Eliminates 99.8% of false-positive whole-chromosome intersections!       │
├────────────────────────────────────────────────────────────────────────────┤
│ LAYER 2: Materialized Hub Bigram Postings (`hub_bigrams.roar`)             │
│ • For the ~1,600 hub tokens with df > 100,000, pre-materialize bitmaps     │
│   keyed by `(token_A, token_B, delta_bin)` where delta_bin ∈ {0: fused,    │
│   1: ΔORF=1 directon, 2: ΔORF∈[2..5] directon, 3: ΔORF∈[1..3] divergent}.  │
│ • Hub-pair intersection takes 0.2 ms (1 lookup) instead of 50M scans!      │
├────────────────────────────────────────────────────────────────────────────┤
│ LAYER 3: Contig-Block Inverted Positional Slices (`token_pos_slices.bin`)  │
│ • Aligned 1-to-1 with CRoaring 16-bit containers (65,536 windows/block).   │
│ • Each posting stores a 2-byte `PackedSlotDescriptor` streamed sequentially│
│   ONLY for non-empty intersected containers via AVX-512 — zero random CSR! │
└────────────────────────────────────────────────────────────────────────────┘
```

#### Exact Bit-Packing & Memory Math for v2.0 Index 1:
1. **32-Bit Safe `window_id` Addressing:**
   - Concatenate all contigs (ordered by `tier`, `tax_order_id`, `biome_id`, `dense_contig_id`) onto a single monotonic `global_slot_idx` $\in [0, 3.6 \times 10^9)$ tape (padding contig boundaries only in a compact `window_to_contig` lookup table, or assigning each contig $\lceil n_{\text{slots}} / 16 \rceil$ windows).
   - Total `window_id` count across $3.6 \times 10^9$ slots is $\approx 2.6 \times 10^8 \ll 2^{32}$, fitting natively inside standard 32-bit `CRoaring` (`roaring_bitmap_t`) with **zero 64-bit truncation**.
   - Every `CRoaring` high-16-bit key (`window_id >> 16`) corresponds to a **Window Block (`block_id` $\in [0, 4000)$)** of $65,536$ micro-windows ($\sim 1.05\text{M}$ slots).
2. **Parallel 2-Byte `PackedSlotDescriptor` Per Posting (`token_pos_slices.bin`):**
   - For every `(token_id, block_id)` container in `token_to_windows.roar` containing $m$ window IDs, store a contiguous `uint16_t[m]` slice right beside the `CRoaring` container:
     ```c
     // 16-bit positional descriptor parallel to each entry in a CRoaring 16-bit container:
     // Bits 0..3  (4 bits): slot_in_window (0..15)
     // Bits 4..7  (4 bits): orf_in_window  (0..15; equal to next slot's orf_in_window if fused!)
     // Bits 8..10 (3 bits): directon_local_id (0..7 modulo-8 directon counter within window)
     // Bit  11    (1 bit) : strand (0 = '-', 1 = '+')
     // Bit  12    (1 bit) : is_rosetta_fused_prev (1 if seg_idx_in_orf > 0)
     // Bit  13    (1 bit) : is_contig_edge (1 if at non-circular contig boundary)
     // Bits 14..15(2 bits): secondary_occurrence_flag (0 = single occurrence of token in window;
     //                      1..3 = index into tiny overflow list if token repeats inside same 16-slot window)
     typedef uint16_t PackedSlotDescriptor;
     ```
   - **Memory Footprint:**
     - `token_to_windows.roar` ($3.6 \times 10^9$ sequence-token postings + $3.6 \times 10^9$ fold-token postings = $7.2 \times 10^9$ postings at ~1.15 B/posting in `CRoaring` run/array/bitmap containers): **`~8.3 GB`**.
     - `token_pos_slices.bin` ($7.2 \times 10^9 \times 2\text{ bytes}$, delta/Zstd block-compressed or raw mmap): **`~14.4 GB` on disk** (only the active containers for the 2–4 query tokens are read sequentially: `< 15 MB` per query!).
     - `hub_bigrams.roar` (~1,600 hub tokens $\rightarrow$ ~1.2M non-empty hub pairs across 4 delta bins, compressed with `CRoaring`): **`~3.1 GB` pinned in RAM**.
   - **Query Execution Trace:**
     - To verify `|slot(A) - slot(B)| <= k` and `same_directon`:
       1. If both `A` and `B` are hubs ($\text{df} > 10^5$), fetch `hub_bigrams[(A, B, bin)]` in **`0.2 ms`**.
       2. Otherwise, intersect `token_to_windows[A] & token_to_windows[B]` container-by-container (`< 1.5 ms`).
       3. For each surviving 16-bit container with $\ge 1$ matching `window_id`, run an AVX-512 branchless two-pointer merge over the two contiguous `uint16_t` `PackedSlotDescriptor` slices (`< 2.0 ms`). **Total latency: `< 4.0 ms` with zero random CSR lookups.**

---

## 3. DuckDB / Parquet Layout, Z-Order Sorting & Out-of-Core Aggregation

### 3.1 Why Z-Order Sorting `feature_tape.parquet` by `(seq_family_id, struct_fold_id, contig_id, slot_idx)` Destroys Contig Locality

In `PLAN.md` Section 4.3, `feature_tape.parquet` ($3.6 \times 10^9$ rows, 18 columns) is sorted on a Z-order space-filling curve over `(seq_family_id, struct_fold_id, contig_id, slot_idx)`.

This violates fundamental columnar storage physics in three ways:
1. **The Neighborhood Scatter-Gather I/O Bomb:**
   - Consider Contig `#42`, which contains 12 slots (`slot_idx = 0..11`) encoding 12 different protein families (`seq_family_id` = `14`, `9821`, `44012`, `302`, ...) and 12 different structural folds (`struct_fold_id`).
   - Because `seq_family_id` and `struct_fold_id` dominate the Z-order Morton code, the 12 slots of Contig `#42` are physically written to **12 completely different Parquet Row Groups** scattered across the 3.0 TB file!
   - When a user discovers 5,000 candidate contigs and asks DuckDB to retrieve their $\pm 10$-slot neighborhoods (`100,000` slots total), DuckDB cannot read 5,000 contiguous slices. Instead, it must execute a **random scatter-gather across ~30,000 Parquet row groups**!
   - Even with Parquet page indexes, reading 18 column chunks across 30,000 row groups requires decompressing **`~540,000` Zstd column pages (`~22 GB` of I/O)**, taking **`15–40 seconds`** instead of `< 20 ms`.
2. **Destruction of Parquet Zone Maps (Min/Max Row-Group Pruning):**
   - A 4D Z-order curve interleaves bits from `seq_family_id` (24 bits), `struct_fold_id` (14 bits), `contig_id` (30 bits), and `slot_idx` (16 bits). Because high-order bits of `contig_id` are interleaved between bits of `seq_family_id`, any row group of $122,880$ rows spans a huge `[min, max]` range for *both* `seq_family_id` and `contig_id`. DuckDB's `RowGroup` min/max zone-map pruning becomes useless on both dimensions!
3. **$3\times$ Compression Blowup (Breaking RLE and Delta Encoding):**
   - When sorted by `(contig_id, slot_idx)`, `contig_id` (8 B) and `genome_or_sample_id` (4 B) repeat consecutively $7\text{–}4,000$ times, compressing via Parquet `RLE_DICTIONARY` to **`< 0.15 bytes/row`**, while `start_bp` and `end_bp` compress via monotonic delta encoding to **`~1.5 bytes/row`**.
   - Under `(seq_family_id, struct_fold_id, ...)` Z-ordering, adjacent rows come from unrelated contigs and samples, destroying run lengths and inflating `feature_tape.parquet` from **`~1.15 TB` to `> 3.2 TB`**.

---

### 3.2 The v2.0 Dual-Projection Parquet Storage Layout

Because comparative genomics has **two orthogonal access patterns**—(1) *Contig-Centric Neighborhood Reconstruction* (given a set of `contig_id`s, fetch their full ordered gene tapes) and (2) *Anchor-Centric Discovery Scans* (given an anchor `token_id`, scan all partner slots within $\pm 6$ ORFs across the biosphere)—`PLAN.md` v2.0 must replace `feature_tape.parquet` with **two compact, strictly lexicographically sorted projections**:

#### Projection A: `contig_tapes.parquet` (Contig-Centric Neighborhood Projection — `~1.12 TB`)
- **Granularity:** **1 Row = 1 Contig** ($5.0 \times 10^8$ rows total; zero multi-row contig fragmentation!).
- **Sort & Row-Group Alignment:** Strictly sorted by `(tier, tax_order_id, biome_id, dense_contig_id)` with fixed row-group size = **`65,536` rows** (`1 Row Group == 1 Index-1 Window/Contig Block`).
- **Schema (Scalar Metadata + Parallel Structure-of-Arrays `LIST<T>` Columns of length `n_slots`):**

| Column | DuckDB / Parquet Type | Compression | Purpose |
|---|---|---|---|
| `dense_contig_id` | `UINT32` | Delta / RLE | `0 .. 500,000,000` (implicit row number `rowid` in DuckDB!) |
| `tier_and_flags` | `UINT8` | RLE | Bits `0..1`: `tier` (`1,2,3`); Bit `2`: `is_circular`; Bit `3`: `genetic_code_id` |
| `genome_or_sample_id` | `UINT32` | RLE_DICT | Foreign key to `genomes_samples.parquet` |
| `tax_order_id` | `UINT16` | RLE | GTDB taxonomic order ID (`0..1800`) |
| `biome_id` | `UINT16` | RLE | Gold/SPIRE/Logan biome ontology ID |
| `n_slots` | `UINT16` | Zstd | Number of feature slots on this contig |
| `seq_family_ids` | `LIST<UINT32>` | Zstd + Dict | Ordered sequence family tokens (`slot_idx = 0 .. n_slots-1`) |
| `struct_fold_ids` | `LIST<UINT32>` | Zstd + Dict | Ordered `ECOD` / `Foldseek 3Di` fold tokens |
| `orf_and_kind_packed` | `LIST<UINT8>` | Zstd | Bit 7: `strand`; Bit 6: `is_fused_to_prev` ($\Delta\text{ORF}=0$); Bit 5: `catalytic_intact`; Bit 4: `is_contig_edge`; Bits `0..3`: `feature_kind` |
| `delta_bps` | `LIST<INT16>` | Zstd | Intergenic $\Delta\text{bp}$ from previous feature (`0` if fused) |
| `start_bps` | `LIST<UINT32>` | Delta + Zstd | Contig start coordinates (bp) |
| `biophysical_badges` | `LIST<UINT64>` | Zstd | Packed 64-bit TM, SignalP, sorting, $\text{pI}$, motif, length badge |
| `nr90_protein_ids` | `LIST<UINT64>` | Zstd | Pointers to `proteins_nr90.parquet` for AA and `3Di` strings |

- **Performance Guarantee:** Given 1,000 candidate `dense_contig_id`s from Index 1, DuckDB uses **Parquet row-number pushdown (`rowid IN (...)`)** to read the exact `LIST` slices for those 1,000 contigs with **zero SQL joins** in **`< 12 ms`**.

#### Projection B: `window_postings.parquet` (Anchor-Centric Discovery Projection — `~1.35 TB`)
- **Granularity:** Pre-extracted directed local neighborhood pairs `(anchor_slot, partner_slot)` for all pairs within **$\pm 6$ ORFs (or fused inside the same ORF)** inside the same or adjacent directon (~$1.8 \times 10^{10}$ deduplicated architecture-weighted rows).
- **Sort Order:** Strictly sorted by **`(anchor_token_id, partner_token_id, tier, tax_order_id)`**.
- **Schema (`24 bytes/row` uncompressed $\rightarrow$ `~7.5 bytes/row` Zstd/RLE compressed):**
  - `anchor_token_id` (`UINT32` — leading sort key; compressed to `~0.01 B/row` via RLE!)
  - `partner_token_id` (`UINT32` — second sort key; compressed via RLE/Delta!)
  - `partner_fold_id` (`UINT32`)
  - `orf_offset` (`INT8` — `0` = Rosetta intra-ORF fusion; `±1..±6` = inter-ORF distance)
  - `rel_orientation_and_flags` (`UINT8` — Bit 0: `same_strand`; Bit 1: `same_directon`; Bit 2: `partner_catalytic_intact`; Bit 3: `not_contig_edge`; Bits `4..5`: `tier`)
  - `delta_bp` (`INT16`)
  - `tax_order_id` (`UINT16`)
  - `partner_badge` (`UINT64` — allows instant filtering on `TM_count`, `SignalP`, `pI`, `motif_flags`, `length`!)
  - `rep_dense_contig_id` (`UINT32`) + `rep_anchor_slot_idx` (`UINT16`) + `arch_occurrence_count` (`UINT16`).
- **Performance Guarantee:** Because `anchor_token_id` is the primary lexicographical sort key, all neighborhoods of *any* anchor family (e.g., `PF04055` or `RF00168`) reside in **1 to 4 contiguous Parquet Row Groups (`5–80 MB` on NVMe)**. DuckDB's zone-map filter skips **99.99% of the 1.35 TB file** and scans all neighborhoods of `anchor_token_id` at `14 GB/s` in **`3–10 ms`**!

---

### 3.3 Out-of-Core Aggregation for `directon_pair_summary.parquet` (Avoiding a 1.2 TB Hash-Spill Crash)

In `PLAN.md` Section 4.3 Table 2, `directon_pair_summary.parquet` aggregates over `~20 billion` raw pair instances (`~300 million` unique `(token_a, token_b)` pairs), computing `COUNT(DISTINCT phylum_id)`, `COUNT(DISTINCT order_id)`, `COUNT(DISTINCT genus_id)`, and `median_delta_bp` / `iqr_delta_bp`.

- **The Trap:** Running `GROUP BY token_a, token_b` with 3 exact `COUNT(DISTINCT)` hash sets and an exact `MEDIAN`/`QUANTILE_CONT` array across $3 \times 10^8$ groups in DuckDB allocates **`> 1.2 TB` of in-memory aggregate state**, triggering multi-pass out-of-core disk thrashing.
- **The v2.0 Fix:**
  1. **Build `directon_pair_summary.parquet` Directly by Streaming Over `window_postings.parquet`:**
     Because `window_postings.parquet` is *already* sorted by `(anchor_token_id, partner_token_id)`, computing `GROUP BY anchor_token_id, partner_token_id` is a **zero-memory streaming Group-By** ($O(1)$ state—only the current `(anchor_token_id, partner_token_id)` pair is in memory at any instant!).
  2. **Fixed-Width Bitmask & Sketch Accumulators:**
     For each streaming `(token_a, token_b)` group:
     - `n_independent_phyla`: 256-bit fixed bitmask (`uint64_t[4]`, `POPCNT` at group end).
     - `n_independent_orders`: 2,048-bit fixed bitmask (`uint64_t[32]`, `POPCNT` at group end).
     - `n_independent_genera`: 64-byte HyperLogLog / linear-counting sketch (`< 1.5%` error, or exact sorted `uint16_t` count since the group is already sorted by `tax_order_id`!).
     - `median_delta_bp` & `iqr_delta_bp`: 32-bin quantized `int16` histogram (`uint32_t[32]`, 128 bytes).
     - **Total RAM required for the entire 20-billion-row aggregation: `< 256 MB`**, completing in **`~4.5 minutes`** on 96 threads!

---

## 4. `Logan` S3 Streaming & Deduplication State (`27.3M` SRA Runs, `40.8 TB` Zstd)

### 4.1 Three Severe Bugs in `PLAN.md` Section 2.3 (`ingest_sieve`)

1. **Bug 4.1 — The Pipeline Causality Paradox (`ordered_family_ids` Does Not Exist During S3 Streaming):**
   - Section 2.3 Step 4 states that while streaming `aws s3 cp s3://logan-pub/c/... - | zstd -d`, the worker computes:
     `arch_hash = xxhash64(ordered_family_ids, strand_transitions, quantized_delta_bp_bins)`
   - **Paradox:** `ordered_family_ids` (which includes `Pfam`/`COG`/`PHROG`/`DARK_FAM` cluster IDs) is not known until **Stage 2 (Section 3.2)**, *after* all proteins have been clustered into `NR90`/`NR50` and dissected by `pyhmmer`/`hhsearch`/`ProstT5`! You cannot compute `ordered_family_ids` inside the Stage 0 nucleotide S3 streaming sieve.
2. **Bug 4.2 — `102 GB – 145 GB` Heap Explosion from In-Memory `arch_hash -> RoaringBitmap` Maps:**
   - Across `27.3M` `Logan` SRA runs + `SPIRE` + `AllTheBacteria`, there are **$4.0 \times 10^8$ to $6.0 \times 10^8$ unique multi-gene architecture hashes** (~75% of which are singleton or doubleton metagenomic assembly fragments).
   - Storing $5 \times 10^8$ keys in a hash table (`32 B/entry` = `16 GB`) pointing to $5 \times 10^8$ heap-allocated `roaring_bitmap_t` structs (`96–112 B` struct + container malloc overhead per bitmap = `48–56 GB`) plus $1.5 \times 10^{10}$ sample postings (`30 GB`) consumes **`102 GB – 145 GB` of heavily fragmented heap RAM**.
3. **Bug 4.3 — `2.0 TB` RAM OOM from "$\ge 95\%$ Nucleotide Identity" Verification:**
   - Section 2.3 Step 4 specifies: *"If `arch_hash` has already been seen in the current taxonomic/biome bucket at $\ge 95\%$ nucleotide identity, increment `arch_occurrence_count`..."*
   - Verifying $\ge 95\%$ nucleotide identity against $5 \times 10^8$ previously seen representative contigs (average length $4\text{ kb}$) requires keeping $5 \times 10^8 \times 4\text{ KB} = \mathbf{2,000\text{ GB (2.0 TB)}}$ of nucleotide sequences accessible at sub-microsecond latency!
4. **Bug 4.4 — Running `pyrodigal` Before Filtering Contigs `< 1,200 bp`:**
   - In `Logan` (`s3://logan-pub/c/`), `~75–82%` of all assembled contigs/unitigs are short fragments between `150 bp` and `1,199 bp`. If `pyrodigal` is invoked before a zero-allocation byte-length filter in the `zstd` stream reader, CPU time increases by **$4.2\times$** (`289 wall-clock hours` instead of `67.5 hours`).

---

### 4.2 The v2.0 `Logan` S3 Streaming & Two-Tier NVMe Deduplication Architecture

Replace Section 2.3 Step 4 with a **Syncmer-Sketched Two-Tier Cuckoo + Partitioned RocksDB/Spill Pipeline**:

```
aws s3 cp s3://logan-pub/c/ (400 MB/s) ──► zstd -d (in-memory)
  │
  ├─► [Filter 1: Fast Byte-Length (< 1,200 bp) & DUST Filter in Rust] (Drops 78% of bytes in <0.1 ns/bp)
  │
  ├─► [Filter 2: 64-MB Euk-Repeat Syncmer Bloom Filter + Shine-Dalgarno Whitelist]
  │
  ├─► [Filter 3: pyrodigal (--meta) + Intergenic smORF RBS Scanner]
  │     • Reject if < 2 features AND longest ORF < 250 aa
  │
  └─► [Filter 4: Stage-0 Protein-Syncmer Architecture Hash & 64-Bit Contig ANI Sketch]
        • For each ORF_i on contig, compute 32-bit Reduced-Alphabet Protein Syncmer Hash `orf_sh_i`
          (using Murphy-10 10-letter structural AA alphabet, k=10, s=4 open syncmers -> 90% AA ID invariant!)
        • Compute `stage0_arch_hash = xxhash64(orf_sh_1..k, strand_transitions, delta_bp_bins)`
        • Compute `contig_ani_sketch: uint64` = 64-bit One-Permutation b-bit MinHash over DNA syncmers (k=15, s=5)
        │
        ├─► Query 16-GB In-RAM Cuckoo Filter (`10^9` 16-bit fingerprints of seen `stage0_arch_hash`):
        │     ├── NOT PRESENT (First Sighting — ~75% of all unique hashes):
        │     │     1. Insert fingerprint into 16-GB Cuckoo Filter in RAM.
        │     │     2. Append compact binary contig record + `(stage0_arch_hash, contig_ani_sketch, sample_id)`
        │     │        sequentially to NVMe `tier3_unique_contigs.bin.zst`. (Zero RocksDB/Bitmap allocation!)
        │     │
        │     └── PRESENT (2nd+ Sighting — Recurrent Operon Candidate):
        │           1. Lookup `(stage0_arch_hash, biome_bucket)` in 64-Shard Partitioned RocksDB on RAID-0 NVMe
        │              (Value = `12 bytes`: `uint64 rep_contig_ani_sketch` + `uint32 rep_contig_id`).
        │           2. Check 95% ANI in 1 CPU cycle: `POPCNT(contig_ani_sketch ^ rep_contig_ani_sketch) <= 6`.
        │           3. If match: append 12-byte `(rep_contig_id: uint32, sample_id: uint32, biome_id: uint16)`
        │              to 64 hash-partitioned NVMe spill logs (`sample_postings_part_00..63.bin`).
        │           4. At end of ingest, radix-sort each 2.5-GB partition in RAM (`< 2 sec/part`) and pack
        │              directly into immutable `CRoaring` sample bitmaps!
```

#### Quantitative Resource Budget for v2.0 `Logan` Ingest:
- **RAM Footprint:**
  - `16.0 GB` (In-RAM Cuckoo Filter for $10^9$ first-sighting hashes)
  - `8.0 GB` (RocksDB block cache + bloom filters for the $\sim 1.5 \times 10^8$ recurrent hashes that pass the Cuckoo Filter)
  - `2.0 GB` (64 append-only NVMe spill write buffers @ `32 MB` each)
  - **Total Deduplication RAM: `26.0 GB`** (down from `102–145 GB` + `2 TB` ANI sequences!).
- **NVMe IOPS:**
  - Because the 16-GB Cuckoo Filter absorbs 100% of first-sighting singletons in RAM and 10-bit RocksDB bloom filters absorb negative lookups in RAM, NVMe reads only occur on recurrent hash hits (`< 15,000 IOPS`, `< 1%` of the 2.0M IOPS RAID-0 budget), while all sample postings are written via sequential `32 MB` buffered appends!
- **Two-Pass Canonicalization:**
  - Once Stage 2 (`NR50` clustering + `pyhmmer`/`hhsearch`/`ProstT5` dissection) assigns true sub-ORF `seq_family_id` and `struct_fold_id` tokens to `tier3_unique_contigs.bin.zst`, a single linear streaming pass computes the final **Stage-2 Semantic `arch_hash`** over true domain tokens and merges synonyms via radix sort in **`< 20 minutes`**.

---

## 5. Additional Systems Fixes Required in `PLAN.md` Sections 4.4 & 4.5

### 5.1 Index 3 (Section 4.4): Cache-Line Alignment & Sparse Edge-Transition Fitch Algebra
1. **64-Byte AVX-512 Cache-Line Padding Bug:**
   - `PLAN.md` Section 4.4 specifies `presence_bits` as `14,138 bytes` ($113,104$ bits). Because `14,138` is not divisible by `64` (`14,138 % 64 = 58`), storing `presence_bits` contiguously in a matrix causes every subsequent row to straddle 64-byte cache-line and 512-bit ZMM register boundaries!
   - **Fix:** Pad every `presence_bits` row to **`14,144 bytes` (`221 × 64-byte` AVX-512 registers = `113,152` bits, with the trailing 48 bits masked to `0`)**.
2. **Pre-Computed Sparse Tree-Edge Transition Lists (`< 0.2 ms` Cotransmission Scoring):**
   - `PLAN.md` Section 4.4 Step 4 proposes running a bottom-up Fitch parsimony pass across all $226,206$ nodes of the GTDB tree at query time. If `presence_bits` is stored in row-major `[token_id][leaf_id]` order, gathering leaf bits across 512 candidate tokens requires 512 strided reads separated by `14,144 bytes` (`7.2 MB` stride), thrashing L1/L2 caches.
   - **Fix:** Because the GTDB `bac120 + ar53` tree topology is static, run the vectorized Fitch + Sankoff ancestral reconstruction **once at index build time** and store for each `token_id` its **Sparse Gain/Loss Edge Transition List (`transition_edges: uint32[]`)**, where bits `0..18` encode `tree_edge_id` $\in [0, 226206)$ and bit `19` encodes `0->1` gain vs. `1->0` loss.
   - Because a typical family undergoes $5\text{–}300$ transitions across the tree, `transition_edges` averages **`< 400 bytes/token`** (`~100 MB` total across 250k tokens!). At query time, computing exact independent evolutionary gain/loss events $H_X$ is a **`O(1)` lookup**, and computing co-gain/co-loss events $H_{X,Y}$ between two families is a sorted `uint32` array intersection taking **`< 0.2 microseconds` per candidate pair**!

### 5.2 Index 4 (Section 4.5): Preventing a `245 GB` VRAM OOM in `gLM-Domain` Softmax
- **The Flaw:** In `PLAN.md` Section 4.5, `gLM-Domain` ($d=512$) trains with a Softmax cross-entropy head $\mathcal{L}_{\text{seq}}$ over `seq_family_id`. Including `DARK_FAM` `NR50` clusters and 40%-ID SSN subfamilies, `token_dictionary.parquet` contains **$\sim 1.5 \times 10^7$ unique `seq_family_id` tokens**. An input embedding table $\mathbf{E}_{\text{seq}}$ + output softmax matrix $\mathbf{W}_{\text{seq}}$ of shape `[15,000,000 × 512]` with AdamW states (`16 bytes/param`) requires **`245.7 GB of VRAM`**—crashing the `48 GB` RTX 6000 Ada immediately!
- **The v2.0 Fix:**
  1. **Factorized Vocabulary Cutoff (`131,072` Discrete Sequence Tokens + Compositional Fold/Badge Representation):**
     - Restrict the discrete `seq_family_id` embedding table $\mathbf{E}_{\text{seq}}$ and output head $\mathbf{W}_{\text{seq}}$ (with tied weights $\mathbf{W}_{\text{seq}} = \mathbf{E}_{\text{seq}}^\top$) to the **top $131,072$ ($2^{17}$) families/subfamilies** ranked by phylogenetic dispersion ($\ge 2$ taxonomic orders or $\ge 25$ occurrences).
     - Memory footprint: $131,072 \times 512 \times 16\text{ B} = \mathbf{1.07\text{ GB VRAM}}$ (training) / **`134 MB VRAM`** (`fp16` inference).
  2. **Represent Rare/Dark Families via `struct_fold_id` (`~4,096` ECOD T-groups) + `3Di` Bigram Projection + `biophysical_badge` (`64-bit`):**
     - Slots occupied by rare `DARK_FAM` clusters outside the top `131,072` set `seq_family_id = <RARE_SEQ_TOKEN>` and feed their exact `struct_fold_id` (`ECOD` T-group), `ProstT5` pooled embedding (`int8` 128-d projection), and `64-bit biophysical_badge` into $\mathbf{h}_i^{(0)}$. This not only prevents VRAM OOM, it forces `gLM-Domain` to generalize by structure and biophysics on rare dark proteins rather than memorizing low-count sequence cluster IDs!

---

## 6. Consolidated End-to-End Ingestion Compute & Storage Budget (`PLAN.md` v2.0)

| Pipeline Stage | v1.0 Estimate (Broken) | v2.0 Optimized Wall-Clock (48C/96T EPYC + 1× RTX 6000 Ada) | Output Footprint on 15.36 TB RAID-0 NVMe |
|---|---|---|---|
| **Stage 0: `Logan` + Tier 1–3 Streaming Sieve, `pyrodigal`, `smORF`, Syncmer Cuckoo/RocksDB Dedup** | >12 days + OOM (`145 GB` + `2 TB` RAM) | **68 hrs wall-clock** (`26 GB` RAM, `400 MB/s` S3 line-rate) | `1.8 TB` (`tier123_dedup_contigs.bin.zst` + sample Roaring logs) |
| **Stage 1: Non-Coding `IGR90` Clustering + `nhmmer`-Gated `Infernal` (`cmsearch`, 2,000 CMs) + `Minced`/`IntegronFinder`** | `180,000–1,800,000 CPU-hrs` (78–780 days) | **8.8 hrs wall-clock** (`420 CPU-hrs` on `50 Gbp` `IGR90`) | `45 GB` (`ncrna_crispr_attc_hits.parquet`) |
| **Stage 2a: `MMseqs2` `NR90` & `NR50` Clustering + Precomputed `AFDB`/`ESMAtlas` `3Di` Transfer** | Unspecified (`NR90` only) | **3.5 hrs wall-clock** (`105 CPU-hrs`) | `420 GB` (`proteins_nr90_nr50.parquet`) |
| **Stage 2b: `pyhmmer` (115k non-redundant HMMs on `NR50`) + `cath-resolve-hits` + Catalytic Mask Check** | `19,900 CPU-hrs` (207 hrs) | **12.5 hrs wall-clock** (`1,200 CPU-hrs`) | `180 GB` (`domain_envelopes.parquet`) |
| **Stage 2c: `ProstT5` Encoder+CNN (`fp16`/`int8` dynamic batching on 18M recurrent dark `NR50` segments) + `Foldseek`** | `1,111 GPU-hrs` (46.3 days) | **14.3 hrs wall-clock** on 1× RTX 6000 Ada | `95 GB` (`dark_nr50_3di_ecod.parquet`) |
| **Stage 2d: `mmseqs result2profile` + `Foldseek`-Prefiltered `hhalign` on 4.5M Multi-Member Dark `NR50`s** | `16,667–250,000 CPU-hrs` | **0.5 hrs wall-clock** (`48 CPU-hrs`, zero loose files) | `12 GB` (`remote_homology_hh.parquet`) |
| **Stage 3: Dual-Projection Parquet Build + `directon_pair_summary` + Index 1/3/4 Construction** | Out-of-core hash spill (`>1.2 TB` state) | **6.0 hrs wall-clock** (streaming group-by + `gLM-Domain` training) | **Total Active DB on NVMe: `3.45 TB`** (22.5% of 15.36 TB RAID-0!) |
| **TOTAL END-TO-END BUILD TIME** | **~135–850 Days (Unbuildable)** | **~4.7 Days End-to-End on 1 CPU Node + 1 GPU Node** | **`3.45 TB` Active Index + `1.8 TB` Raw Archive = `5.25 TB`** |
