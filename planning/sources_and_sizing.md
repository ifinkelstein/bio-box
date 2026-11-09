# Source ladder and sizing for the staged build (wilkcomp02 first, Hetzner later)

Status: plan only. Nothing here has been run. Every rate is **unmeasured – measure in S0**: each one is given as a formula so the S0 microbenchmarks can replace it. Counts marked "est." are extrapolations, and the source of each is named. Container/CLI name: `bio-box`.

## 1. The envelope on wilkcomp02

From `WILKECOMP.org` (probed 2026-09-30) and `fragment-defense/old-plan.2026-10-18.bak` §9.3a/b.

| Resource | Usable value | Planning consequence |
|---|---|---|
| Cores | 56 physical Haswell E7-4850 v3, 2.2 GHz, AVX2 (no AVX-512). 40 threads by day, 56 at night and weekends | Realistic throughput ≈ 14 h × 40 + 10 h × 56 = 1,120 core-h/day peak. **Plan with 800 core-h/day** (70% duty). |
| Per-core speed | About 1/3 of a modern desktop core | Modern published rates are multiplied by 3 below (`H = 3`) |
| RAM | 1,007 GB physical; the norm is a 200 GB `systemd-run` cgroup (≤20% by day, ≤65% only overnight for a few hours) | **Hard design limit: 200 GB per job.** |
| Local `/tmp` | ~940 GB ext4, volatile, not shared | All MMseqs2/HMMER/Infernal working sets go here. **Hard limit: about 800 GB per job, leaving headroom** |
| NFS `/stor/scratch/Wilke` | 30 TB, 1.2 TB used, not backed up | Holds the sequence store and stage outputs. Read and written as few large files only |
| NFS `/stor/work/Wilke` | 75 TB shared lab quota (62% used), backed up | Manifests, ledgers and small durable outputs only (GB, not TB) |
| NFS I/O | One server, 10 GbE, shared by three hosts | No per-genome files on NFS. Stream tarballs and write sharded Parquet/bgzip of ≥1 GB |
| Watchdog | Renices any single PID after ~369 CPU-minutes and mails rctf-support | Monolithic multithreaded tools (MMseqs2, Foldseek, HMMER on large shards) **launch at `nice 19`**. Fan-out workers run at `nice 15` and are sized so no PID exceeds ~300 CPU-min |
| Scheduler | None; singularity present; docker group use unclear | `bio-box` ships as an Apptainer `.sif`. The orchestrator is Snakemake with a local executor and explicit `--cores` caps |
| Open files | Hard limit 4,096 (8,192 allowed via `ulimit -n`) | Sharded outputs, not millions of small files |

## 2. Unit rates (all unmeasured – replace in S0)

`H = 3` converts a modern-core rate to wilkcomp02. Each rate is given as core-hours (core-h) per unit of input.

| Tool / step | Formula (core-h) | Published anchor | Estimate on wilkcomp02 |
|---|---|---|---|
| pyrodigal, single mode (isolate/complete genomes) | `Gbp × 1e3 / (R_pyro × 3600)` | Not found numerically. The pyrodigal benchmarks page reports SIMD gains of ⅓–½ over Prodigal but gives no Mbp/s figure ([docs](https://pyrodigal.readthedocs.io/en/v3.2.1/benchmarks.html)) | `R_pyro` = 0.1–0.3 Mbp/s/thread → **0.9–2.8 core-h/Gbp**. Meta mode and pyrodigal-gv: assume ×2–4 until measured |
| pyhmmer hmmsearch, Pfam-A | `Gaa × r_pfam` | pyhmmsearch README: 4,977 proteins vs Pfam in 144 s on 1 thread (≈1.6 Maa assumed at 330 aa/protein) ([pyhmmsearch](https://github.com/jolespin/pyhmmsearch)); PyHMMER paper ([PMC10159651](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10159651/)) | `r_pfam` ≈ 0.024 core-h/Maa modern × H → **≈73 core-h/Gaa** |
| pyhmmer, NCBIfam | `Gaa × 1.3 r_pfam` | Scaled by model positions (assumption) | **≈95 core-h/Gaa** |
| pyhmmer, PHROGs (only on viral, plasmid and prophage proteins) | `Gaa_mge × 1.5 r_pfam` | Scaled by model count (38,880; assumption) | **≈110 core-h/Gaa** |
| Infernal, full Rfam on whole genomes | `Mbp × 2 strands × N_cm × s_cm / 3600` | Infernal 1.1 paper: 0.6 s/Mb/query; about 1 h per 2–5 Mb genome against Rfam 11 (2,208 CMs) on a quad-core desktop ([PMC3810854](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3810854/)) | Rfam 15.x has about 4,200 CMs → **≈1,400–5,600 core-h/Gbp**. Run this mode only in S0, as the recall reference |
| Infernal, production mode (IGRs ≈12% of bp, `--rfam` filters, prokaryotic/viral subset) | Same × f_IGR × f_subset × f_rfam | Assumes f_IGR 0.12, f_subset 0.5, f_rfam 0.3–1 | **≈30–340 core-h/Gbp of genome**. The widest and most consequential range; S0 must measure it |
| geNomad end-to-end | `Gbp × g` | 10,000 scaffolds (18.3 Mb) in 242 s on c5.12xlarge with **thread count not stated** ([PMC11324519](https://pmc.ncbi.nlm.nih.gov/articles/PMC11324519/table/Tab1)) | **11–530 core-h/Gbp** (1 to 48 threads × H). Also test `--disable-nn-classification` |
| CheckM2 | `n_genomes × c` | No published per-genome time found ([Biowulf doc](https://hpc.nih.gov/apps/checkm2.html)) | 0.05–0.3 core-h/genome, assumed. **Skip for GTDB genomes**: GTDB metadata already carries CheckM2 values (verify field in R226 metadata) |
| MMseqs2 `linclust` (90% or exact-ish) | `N_seq/1e6 × l` | Linclust paper: ~1.6B sequences in ~10 h on one server (Steinegger & Söding 2018; recalled, verify) | `l` ≈ 0.25 core-h per 1e6 seqs per pass |
| MMseqs2 cascaded `cluster` to 50% | `N_seq/1e6 × 10–30 l` | Assumption | 2.5–7.5 core-h per 1e6 seqs |
| MMseqs2 linclust memory | `N_seq × kmer_per_seq (21) × 16 B` | MMseqs2 wiki memory formula (recalled, verify) | **≈336 B/seq.** This sets the cgroup ceiling: 200 GB ≈ 6e8 sequences per linclust run |

**Conversion constants (assumptions, measure in S0):**
- Proteins: aa ≈ 0.30 × bp (88% coding ÷ 3). 330 aa per protein on average.
- IGRs ≥ 40 bp: ≈12% of bp in prokaryotes, ≈8% in phages.

**Base-scan unit ([SEARCH_PLAN SD-05](search/SEARCH_PLAN.md)):** library scans (the pyhmmer rows above) run on exact-unique proteins through S2C, and on NR90 representatives from S3 on, with hits projected to members through `cluster_member_aln`. Members whose projected bitscore is within 10 bits of the cutoff, or whose alignment covers < 90% of the domain envelope, are re-scored exactly. So `Gaa` in the pyhmmer formulas is exact-unique residues at S0–S2C and new NR90-rep residues plus the re-score fraction (unmeasured) from S3 on. Stored hits are per protein after projection.

## 3. Sequence store (new requirement): every contig, protein and IGR kept and retrievable

- **What is kept:** raw contig nucleotide sequences for **every** contig, including dedup members, which are never discarded. Predicted proteins and IGRs are kept too. All are reachable by any downstream project through a stable sequence UID.
- **Format:** primary is bgzip FASTA shards (`.fai` + `.gzi`, ~2^30 bases, ~0.33 GB each) + Parquet `protein_seq` per `data_layer_and_tests.md` §2A, decided by the S0 bake-off K-13 (DEC-17). zstd-seekable is not used. The size model below holds for either format.
- **Lookup table:** `seq_index.parquet` maps `seq_uid` (BLAKE3-128 of the uppercased sequence; exact-byte `raw_digest` is sha256) to `{source, accession, version, file, record_name, length, md5}`. Proteins and IGRs also carry `{parent_contig_uid, start, end, strand}`.
- **Size model:**

| Sequence type | Compressed size | Index overhead |
|---|---|---|
| Nucleotide (contigs, IGRs) | **0.309 B/bp** (*measured*, golden set, 2026-11-09; shards with indexes; the estimate was 0.25–0.30, using 0.27) | `.fai` ≈ 50 B/record, `.gzi` ≈ 0.25 MB/GB |
| Protein | **0.73 B/aa** (*measured*, golden set, 2026-11-09; `protein_seq`, 227 bytes per distinct protein; the estimate was 0.55) | Same as above |

- **Per-stage store size:**
  - `store_GB = 0.309·Gbp_nt + 0.73·Gaa + 0.309·0.12·Gbp_nt + 50 B·N_records`
  - That is ≈ **0.565 GB per Gbp of genome** (*derived* from the measured constants above: 0.309 nt + 0.219 aa (0.73 × 0.30) + 0.037 IGR (0.309 × 0.12); the estimate was 0.46: 0.27 nt + 0.165 aa + 0.032 IGR).
- **Location:** on `/stor/scratch` (not backed up). Durability comes from each source's pinned checksum and from a second copy of the store on the lab storage server or an object bucket once it exceeds ~1 TB (S6). The `seq_index` and manifests go to `/stor/work`, which is backed up.

## 4. The ladder

### Changes from the candidate order

1. **IMG/VR HQ vOTU reps (and IMG/PR) move to S4, ahead of all GTDB reps.**
   - Size: they are small (231,408 HQ vOTUs, est. ~12 Gbp).
   - Order logic: completing the viral and plasmid layer before the big bacterial expansion means geNomad/CheckV regime calls and vehicle dedup can be validated on MGEs first.
   - Licensing: JGI policy makes them **local-only**, never in a redistributable pack. The constraint is the license, not compute.
2. **"RefSeq complete" is folded into S6 as a flag, not a separate source.**
   - Overlap: GTDB genomes are drawn from RefSeq/GenBank, so most RefSeq complete genomes are already GTDB members.
   - What RefSeq adds: the `assembly_level = Complete Genome` flag from `assembly_summary_refseq.txt`, used for `hq_closed` masks.
   - Pangenome planes need the multi-strain set, which is S6.
3. **S2C (complete prokaryotic genomes) is inserted between S2 and S3** (2026-10-29, [SEARCH_PLAN SD-16, §6.1](search/SEARCH_PLAN.md)).
   - Selection: GTDB R226 genomes with `ncbi_assembly_level = Complete Genome` in `bac120_metadata_r226.tsv` / `ar53_metadata_r226.tsv`, ≤20 per species (CheckM2 quality, then ANI diversity). The count is an S−1 task.
   - Complete species reps enter at S2C and are not re-ingested at S3 or S5; complete non-rep members enter at S2C and are not re-ingested at S6. S5 and S6 totals are unchanged; their deltas shrink.
   - The RefSeq `assembly_level` flag (item 2) still arrives at S6 for the `hq_closed` masks; S2C selects on the GTDB metadata field, not on `assembly_summary_refseq.txt`.

### Stage table

All rates are unmeasured; "Core-h" is the range from §2 at the stage's scale.

| Stage | Sources (version) | n genomes | Gbp | Download GB | Core-h, est. | Disk: seq store / everything else | License flag | New capability |
|---|---|---|---|---|---|---|---|---|
| **S0** | Golden set (critique_8 table: MG1655, B. sub 168, T4, λ, *M. genitalium*, SR1, Lak, *M. acetivorans*, *H. volcanii*, *S. acidocaldarius*, H37Rv, *B. theta*, SF370, Obelisk) + CAMI II marine subset | 15 + CAMI subset (*measured*, golden set, 2026-11-09) | 0.0365 (36,475,789 bp, *measured*, golden set, 2026-11-09) (+ CAMI subset, est. 1–3) | <2 | 100–1,500 (dominated by full-Rfam reference runs and microbenchmarks) | 0.03 GB (+ ~1 GB CAMI) / <20 GB | Open (RefSeq/INSDC; CAMI CC-BY) | Calibrates every rate above. End-to-end truth: operons vs RegulonDB, frameshift/Sec/intron/code-4/code-25 handling, leaderless starts. DuckDB/Parquet schema and tests at toy scale |
| **S1** | INPHARED 14Apr2025 (or newer) complete phages | est. 30–40k (14,244 at Jan 2021; current count must be read from `*_data.tsv`) | est. 1.7–2.2 | ~0.6 | 250–2,500 | 0.9 GB / ~15 GB | Sequences INSDC (open). INPHARED scripts AGPL-3.0 | Phage gene calling (pyrodigal-gv, codes 11/15), PHROGs, introns/HEGs, frameshifts (gpG-T), viral regime class; first `slot_pairs` and `pair_stats` |
| **S2** | PLSDB 2024_05_31_v2 | 72,360 plasmids | est. 5–7 ("~7 GB of sequence information") | ~2–7 | 800–6,000 | 2.8 GB / ~40 GB | Records from NCBI Nucleotide (open). Database license not stated → cite, redistribute derived tables only after check | Plasmid backbone/cargo, MOB/Rep, conjugation, integrons/attC, AMR passengers. Benchmark-5 trap material appears; vehicle dedup on real plasmids |
| **S2C** | GTDB R226 metadata (`bac120/ar53_metadata_r226`, `ncbi_assembly_level = Complete Genome`, ≤20 per species) + NCBI datasets for the sequences | 25,458 (*measured*, GTDB R226 metadata, 2026-11-09; 56,387 contigs, 97.7 M proteins; the estimate was 15–30k) | 106.4 (*measured*, GTDB R226 metadata, 2026-11-09; the estimate was 60–120) | 29 (*derived*: 0.273 GB/Gbp, the S5 tarball ratio 127 GB / 465.5 Gbp, × 106.4 Gbp; the estimate was 17–38) | 4,600–18,300 *derived*: 106.4 Gbp × 43–172 core-h/Gbp (S5 rate: 20,000/465.5 and 80,000/465.5); SEARCH_PLAN §6.1's 7,500–20,000 lies inside this range. Library scans on exact-unique proteins (§2) | 60 GB (*derived*: 0.565 GB/Gbp × 106.4 Gbp) / ~80–240 GB (S5 ratio 1.3–2.0 GB/Gbp; computed from the superseded estimate; recompute); search pack ~50–105 GB (below; computed from the superseded estimate; recompute) | Genomes INSDC (open); GTDB CC BY-SA (verify) | Closed chromosomes with exact neighborhoods (no edge censoring); multi-strain dependence-aware counts; search tier ST2 (user tokens, on-demand anchor statistics, `nbr score`, local HH-suite) |
| **S3** | GTDB R226 species reps < 2 Mbp (Patescibacteria ≈6,435 species; DPANN; Mycoplasmatota; endosymbionts), selected from `bac120/ar53_metadata_r226` `genome_size` | 38,750 (*measured*, GTDB R226 metadata, 2026-11-09; 4.17 M contigs, 53.2 M proteins; the earlier estimate was about a quarter of this) | 52.4 (*measured*, GTDB R226 metadata, 2026-11-09; the earlier estimate was about a fifth of this) | extract from the S5 tarball, or ~4 via NCBI datasets (computed from the superseded estimate; recompute) | 1,000–10,000 (computed from the superseded estimate; recompute) | 5.2 GB / ~60 GB (computed from the superseded estimate; recompute) | GTDB: CC BY-SA 4.0 (verify) for taxonomy; genomes INSDC | Genetic code detection (code 25 in BD1-5/SR1, code 4 in Mollicutes), reduced genomes, leaderless archaea, small-genome NOGD sanity. First tree-ordered bitvectors (sparse) |
| **S4** | IMG/VR v4 `-hc` HQ vOTU reps (231,408 vOTUs with HQ member; 495,576 UViGs ≥90% complete) + optionally IMG/PR (699,973 plasmids) | 0.23M (+0.7M) | est. 12 (+ est. 10–30) | est. 5–15 | 4,500–17,000 (HMM + Rfam on IGRs; geNomad/CheckV not rerun, because IMG/VR ships them) | 5.5–19 GB / ~150 GB | **Gated: JGI Data Usage Policy.** Local-only layer, flagged `redistributable = false`. IMG/PR moved to the JGI Data Portal (old Genome Portal sunset mid-2025) | Uncultivated viral diversity, host predictions, prophage/provirus boundaries, viral spacer links |
| **S5** | GTDB R226 all species reps (`gtdb_genomes_reps_r226.tar.gz`, 127.30 GB); complete reps already loaded at S2C and S3 reps are not re-ingested | 143,614 total (remaining = 143,614 − S3 − complete reps at S2C; 24.95 M contigs, 437 M proteins, *measured*, GTDB R226 metadata, 2026-11-09) | 465.5 (*measured*, GTDB R226 metadata, 2026-11-09; the estimate was 400–450) | 127 | **20,000–80,000** (25–100 days at 800 core-h/day; computed from the superseded estimate; recompute). Long poles: Rfam on IGRs, geNomad (until measured), Pfam on NR90 | 207 GB / ~0.6–0.8 TB on scratch (computed from the superseded estimate; recompute); **local `/tmp` peak 500–900 GB** | GTDB CC BY-SA (verify); genomes INSDC | Tree-wide phylogenetic dispersion, species-level presence bitvectors (core plane only), clade-local NOGD, full token dictionary (`seq_family`, `subfamily_50`) |
| **S6** | GTDB R226 remaining non-rep members (588,861 minus the complete members loaded at S2C) via NCBI datasets, + RefSeq `Complete Genome` flag | 732,475 total (111.5 M contigs, 2,436 M proteins, *measured*, GTDB R226 metadata, 2026-11-09) | 2,588 cumulative (*measured*, GTDB R226 metadata, 2026-11-09; the estimate was 2,700–3,100) | est. 800–950 (computed from the superseded estimate; recompute) | est. 60,000–1,600,000 (geNomad dominates if it is not cheap; computed from the superseded estimate; recompute); Pfam on new NR90 reps only, via `clusterupdate` | **+1.34 TB (cumulative ≈1.56 TB)**; the delta shrinks by the S2C members' share, the cumulative does not (computed from the superseded estimate; recompute) / +2–4 TB | INSDC (open) | Pangenome planes (`core_chr`, `pan_any`, `mge_accessory`), `hq_closed_or_multi_strain` mask, strain-level vehicle dedup |
| **S7+** | Off box: AllTheBacteria (2.44M, MIT), SPIRE, GlobDB extras, Logan (AWS us-east-1 compute-to-data) | Millions | Tbp–Pbp | TB–100s TB | Hetzner/AWS | Hetzner object storage | Mixed; Logan has no stated license → NIH/NCBI policy | Tier 3: unbinned diversity, rare architectures |

### Disk per layer for S5 (the last stage that fits on the box)

| Layer | Size (est.) | Where |
|---|---|---|
| Raw download (genome tarball; GTDB protein tarballs not needed because we re-call) | 127 GB | scratch, read-only |
| Sequence store: nt + aa + IGR, bgzip + fai | ~207 GB | scratch |
| Per-genome annotations: pyrodigal GFF → Parquet; Pfam/NCBIfam hits ≈1.5 hits/protein × 60 B over 437 M proteins (*measured*, GTDB R226 metadata, 2026-11-09; scanned on NR90 reps, projected to every protein, §2) | ~50 GB (computed from the superseded 4.8e8 estimate; recompute) | scratch |
| Clusters: MMseqs2 DBs (~2× aa bytes) + NR90/NR50 membership Parquet (36.5 B per protein per level, *measured*, golden set, 2026-11-09; the estimate was 16 B/protein) | ~300 GB DB (transient; computed from the superseded estimate; recompute) + ~32 GB membership (*derived*: 36.5 B × 2 levels × 437 M proteins; the estimate was ~8 GB) | DB on `/tmp` during runs, membership on scratch |
| Parquet tapes (est. 6e8 slots (superseded estimate; the S5 protein count is 437 M) × 30.4 B per slot, *measured*, golden set, 2026-11-09; the estimate was ~40 B) + `slot_pairs` + `pair_stats` | ~25–60 GB (computed from the superseded estimate; recompute) | scratch |
| DuckDB catalog / views | <1 GB | scratch |

### S5 protein residues (reconciled 2026-10-29)

Three values were in use ([`search/06_infra_inventory.md`](search/06_infra_inventory.md) §3 item 26):

| Value | Derivation | Status |
|---|---|---|
| **120–135 Gaa** | 0.30 aa/bp × 400–450 Gbp (§2 constant, the same formula used at every stage; computed from the superseded Gbp estimate, S5 is 465.5 Gbp *measured*, GTDB R226 metadata, 2026-11-09; recompute) | **Working range** |
| 145 Gaa | Stated in §5 (`/tmp` line) without a derivation | Superseded |
| 158 Gaa | ~4.8×10⁸ proteins × 330 aa | Superseded: it multiplies the conservative protein count by the mean length, so it double-counts the margin |

- The protein count consistent with the working range is 3.6–4.1×10⁸ (*derived*: 120–135 Gaa ÷ 330 aa). The ~4.8×10⁸ figure is kept as the **upper bound for memory and `/tmp` sizing** (linclust 160 GB, §5), not as a residue source.
- Protein store at S5: 0.73 B/aa × 120–135 Gaa = 88–99 GB (*derived*; the 120–135 Gaa is computed from the superseded Gbp estimate; recompute; matches the 0.219 GB/Gbp term in §3; the estimate was 0.55 B/aa and 66–74 GB). `data_layer_and_tests.md` §2A.4 gives 45 GB and needs the same update.
- S0 measures aa/bp and mean protein length on the golden set; S2C (106.4 Gbp of complete genomes, *measured*, GTDB R226 metadata, 2026-11-09) is the first stage large enough to refit both before S5.

### Search pack sizing ([SEARCH_PLAN §3.2](search/SEARCH_PLAN.md))

Built once per release by `bio-box pack build`, as ≥256 MB shards on NFS with a manifest (SD-06); at S0–S2 each target is a single file, below the shard floor. Figures are cumulative for the release. All are *estimates* unless marked *derived*; NR90 and NR50 counts are **estimates (±2×)** from `06_infra_inventory.md` §4 and are not given in any other plan file. They are replaced by S0–S3 measurements.

Byte constants (*derived* from the SEARCH_PLAN §3.2 S5 rows): zstd protein FASTA ≈ 0.6 B/aa (20–35 GB / 33–56 Gaa); MMseqs2 `--prefilter-mode 1` index ≈ 1.1 B/aa (35–60 GB / 33–56 Gaa); DIAMOND DB ≈ 0.8–1.1 B/aa (100–150 GB / 120–135 Gaa; 1.1 used at S0–S2C); zstd nucleotide FASTA ≈ 0.27 B/bp (14 GB / 48–54 Gbp). Dark ORFs ≈ 0.07 Gaa per Gbp of genome (*estimate*: SEARCH_PLAN ~30 Gaa / 400–450 Gbp at S5, ~200 Gaa / 2.7–3.1 Tbp at S6; measure at S0).

| Target (tier it starts) | S0 (ST0) | S1–S2 (ST1) | S2C (ST2) | S5 (ST3) | S6 (ST4) |
|---|---|---|---|---|---|
| All / exact-unique proteins (`unique.dmnd`; the S0 target, recall check after) | 15 Maa, all proteins (4.5×10⁴), no reduction; <0.1 GB | ≤2.0–2.8 Gaa (all proteins; exact-unique not estimated); 2.2–3.1 GB | 8.6–16 Gaa: S0–S2 plus S2C 2–4×10⁷ exact-unique × 330 aa; 9.5–18 GB | ~120–135 Gaa (reps are nearly all unique) plus S2C members; ~100–150 GB | ~1×10⁹ proteins, ~330 Gaa; ~330 GB, Hetzner only |
| NR90 reps (`nr90.faa`) | not built (ST0 is unreduced) | 1.4–4.1×10⁶ → 0.4–1.2 Gaa; 0.25–0.7 GB | 2.7–9.1 Gaa; 1.6–5.5 GB | 2.4–3.4×10⁸ → 80–110 Gaa; 48–66 GB | 3.5–6×10⁸ → 115–200 Gaa; 70–120 GB |
| NR50 reps (`nr50.faa` + MMseqs2 index) | not built | 0.65–2.6×10⁶ → 0.2–0.75 Gaa; 0.1–0.45 GB + 0.2–0.8 GB index | 1.0–4.1 Gaa; 0.6–2.5 GB + 1.1–4.5 GB index | 1.0–1.7×10⁸ → 33–56 Gaa; 20–35 GB + 35–60 GB index | 1.2–2.5×10⁸ → 40–80 Gaa; 24–48 GB + 44–88 GB index |
| `segments.faa` (≈1.1× NR90 residues) | not built | 0.45–1.3 Gaa; 0.3–0.8 GB | 3.0–10 Gaa; 1.8–6.0 GB | 88–121 Gaa; 53–73 GB | not estimated (blank in SEARCH_PLAN §3.2; excluded from the total) |
| `igr.fna` + `igr90.fna` | ~6 Mbp; <0.01 GB | 0.74–1.0 Gbp; 0.2–0.3 GB | 7.9–15.4 Gbp; 2.1–4.2 GB | 48–54 Gbp; ~14 GB | 324–372 Gbp; 87–100 GB |
| `dark_orf.faa` | ~3.5 Maa; <0.01 GB (measure here) | 0.47–0.64 Gaa; 0.3–0.4 GB | 4.7–9.0 Gaa; 2.8–5.4 GB | ~30 Gaa; ~18 GB | ~200 Gaa; ~120 GB |
| `spacers.fna`, `models/` | small | small | small | small | small |
| `hh/` (Pfam-A HHM, PDB70; fixed reference) | not built | not built (from ST2) | 30–60 GB | 30–60 GB | 30–60 GB |
| `3di_nr50/` (ProstT5) | not built | not built | not built (from ST3) | ~50 GB | 60–120 GB |
| **Total** (*derived* sum of rows) | **<0.1 GB** | **~3.5–6.5 GB** | **~50–105 GB** (20–45 GB without `hh/`) | **~370–530 GB** | **~0.76–1.0 TB** (~0.43–0.65 TB without `unique.dmnd`) |

Derivations:
- **S1–S2 NR90/NR50:** S1 NR90 0.5–1.5×10⁶ × 215 aa (phage mean 200–230 aa) = 0.11–0.32 Gaa, plus S2 NR90 0.9–2.6×10⁶ × 330 aa = 0.30–0.86 Gaa. S1 NR50 0.2–1.0×10⁶ × 215 aa plus S2 NR50 (0.1–0.25 × 4.5–6.4×10⁶) × 330 aa. Cross-stage merging is ignored, so these are upper ends.
- **S2C proteins:** 106.4 Gbp × 0.30 = 31.9 Gaa; ÷ 330 aa = 9.7×10⁷ proteins (*derived*; the S2C protein count is 97.7 M, *measured*, GTDB R226 metadata, 2026-11-09). The S2C column of the table above was computed from the superseded 60–120 Gbp estimate; recompute. Exact-unique 2–4×10⁷ applies the S6 unique fraction (35%, §5) because S2C holds up to 20 strains per species; it may be higher, since many species have a single complete genome. SEARCH_PLAN §6 ST2 "~2–4×10⁷ proteins" is this exact-unique count.
- **S2C NR90/NR50:** the S6 ratios to exact-unique (NR90 3.5–6×10⁸ / 1×10⁹ = 0.35–0.6; NR50 1.2–2.5×10⁸ / 1×10⁹ = 0.12–0.25) applied to 2–4×10⁷, × 330 aa, plus the S0–S2 values.
- **S2C IGR and dark ORFs:** 12% × 106.4 Gbp = 12.8 Gbp; 0.07 × 106.4 = 7.4 Gaa; plus S0–S2 (*derived*).
- **S5 and S6 rows:** SEARCH_PLAN §3.2 values, on the conservative 4.8×10⁸ protein count (NR50 and NR90 there are 0.2–0.35× and 0.5–0.7× of it). On the reconciled 120–135 Gaa they would be ~15–25% lower.
- **S5 total:** the sum of SEARCH_PLAN's own S5 rows is ~370–530 GB, which SEARCH_PLAN now uses (it first quoted "about 200–350 GB"). Either way it is under 2% of the 30 TB scratch. Per-run `/tmp` working set (NR50 FASTA + index + IGRs): ~70–110 GB at S5, ~4–11 GB at S2C, inside the 150 GB search-lane quota (SD-15).
- **Base scans versus search targets:** library scans use exact-unique proteins through S2C and NR90 reps from S3 (§2, SD-05); the pack's `unique.dmnd` is a recall-check target from S5 on, not the scan unit.

## 5. Where the wilkcomp02 envelope stops being enough

- **S0–S4, including S2C, fit with margin.** Each is under 20k core-h (S2C's high end is ~18.3k, about 23 days; *derived*: 18,300 ÷ 800), under 200 GB of working set, and under 1 TB of disk.
- **S5 is the boundary stage.** It fits on disk (~1 TB of scratch), but it is borderline on three axes:
  1. **Memory.** MMseqs2 linclust over ~4.8e8 proteins (the upper bound, §4 reconciliation) ≈ 4.8e8 × 336 B ≈ **160 GB**, at the edge of the 200 GB cgroup.
     - Mitigation: exact-sequence dedup first. Then either shard linclust by phylum and merge with `clusterupdate` (MMseqs2 ≥ 18-8cc5c, per the pod notes), or get a written overnight exception to 400 GB from Claus Wilke / rctf-support.
  2. **Local `/tmp`.** MMseqs2 DB + index + tmp for 120–135 Gaa (§4 reconciliation; 145 Gaa superseded) is est. 500–900 GB, against ~940 GB free. The search lane's 150 GB `/tmp` quota (SEARCH_PLAN SD-15) comes out of the same disk.
     - Mitigation: `--split-memory-limit`, shard-wise processing, delete tmp between passes.
  3. **Wall-clock and the watchdog.** 25–100 days of 56-core nights.
     - Every MMseqs2 run must be at `nice 19`.
     - HMMER/Infernal work runs as ~10k short shard processes, each kept under ~300 CPU-min, to stay inside etiquette.
     - The range collapses once S0 measures `R_pyro`, `s_cm` and the geNomad rate.
- **S6 breaks the envelope.** Reasons:
  1. **Memory.** linclust over ~3e9 proteins needs ≈1 TB RAM, which is above the cgroup and close to physical RAM. Even after exact dedup (est. 35% unique → ~1e9 seqs ≈ 340 GB) it exceeds the 200 GB norm.
  2. **Local `/tmp`.** The working sets exceed 940 GB.
  3. **NFS.** ~1 TB of input downloads and ~1.3 TB of store writes go over the shared NFS server. That is acceptable only throttled (`ionice -c3`, one stream).
  4. **geNomad on 588k genomes.** At anything but the low end of its rate range, this is CPU-years.

  What can still run on box: pyrodigal plus exact-dedup and the sequence store; the per-species pangenome planes, which are embarrassingly parallel by species. What moves to Hetzner (AX162-class, 48 Zen 4 cores, AVX-512): the global clustering and geNomad.
- **S7+ is off-box by design.** The sizes are TB–PB, and Logan transfer only makes sense next to its bucket.

## 6. S0 measurement list (outputs replace §2)

All runs go through `hyperfine` or `/usr/bin/time -v` inside the `bio-box` `.sif`, on wilkcomp02, at nice 19, with logged thread counts:

1. **pyrodigal** single, meta and pyrodigal-gv on MG1655, T4, SR1 and a CAMI shard → `R_pyro` (Mbp/s/thread).
2. **pyhmmer** Pfam-A 37.x, NCBIfam and PHROGs on 1 Maa of MG1655 + T4 proteins at 1/8/40 threads → `r_pfam` and the scaling curve. Record Z (search space) for E-value comparability; ad hoc searches use the fixed `Z_prot = 1e8`, `Z_nt = 1e5` Mb (SEARCH_PLAN SD-04, confirmed at S−1).
3. **cmsearch/cmscan**:
   - Full Rfam on MG1655 and *B. subtilis* 168 whole genomes (the reference run).
   - The same with IGR-only, prokaryotic subset and `--rfam`.
   - Output: `s_cm` and recall of production mode vs reference (target ≥ 0.98 of reference hits with E < 1e-5).
4. **geNomad** on 100 Mbp of mixed S1/S2/S3 sequence at 1 and 8 threads, with and without NN → `g`.
5. **MMseqs2** linclust and cluster at 50% on 1e7 proteins → `l`, peak RSS and tmp bytes. Check the 336 B/seq memory formula.
6. **Sequence store:** bgzip vs zstd-seekable bytes/bp and random-access latency (1,000 random records from NFS) → checks the golden-set values 0.309 B/bp and 0.73 B/aa (*measured*, golden set, 2026-11-09; the estimates were 0.27 and 0.55).
7. **DuckDB:** Parquet scan and lookup on the S0 tapes, as the template for the data-layer tests at each stage.

## 7. Facts verified for this plan (2026-10-28)

- **GTDB R226:**
  - 732,475 genomes, 143,614 species clusters (136,646 bacterial, 6,968 archaeal) ([stats](https://gtdb.ecogenomic.org/stats/r226)).
  - `gtdb_genomes_reps_r226.tar.gz` = 127.30 G; proteins aa 87.49 G; proteins nt 124.19 G ([listing](https://data.gtdb.ecogenomic.org/releases/release226/226.0/genomic_files_reps/)).
  - Patescibacteriota ≈ 6,435 species clusters (read by position from a jumbled table; confirm).
- **INPHARED:** latest release file set dated 14Apr2025, hosted at `millardlab-inphared.s3.climb.ac.uk`; repo is AGPL-3.0 ([GitHub](https://github.com/RyanCook94/inphared)). The count was not available; 14,244 genomes at Jan 2021 ([preprint](https://www.biorxiv.org/content/10.1101/2021.05.01.442102.full.pdf)).
- **PLSDB 2025:**
  - 72,360 plasmids (v2024_05_31_v2) from NCBI Nucleotide, "about 7 GB of sequence information".
  - No database license stated; the article is CC BY 4.0 ([PMC11701622](https://pmc.ncbi.nlm.nih.gov/articles/PMC11701622/)).
- **IMG/VR v4:**
  - 15,663,652 UViGs; 8.7M vOTUs (8,606,551 in results); 231,408 vOTUs with a HQ representative; 495,576 UViGs ≥90% complete.
  - The `-hc` bundle holds 5,576,197 high-confidence UViGs. No license named ([PMC9825611](https://pmc.ncbi.nlm.nih.gov/articles/PMC9825611)).
  - JGI data policy applies ([JGI](https://jgi.doe.gov/data-policy-support/data-policy)).
- **IMG/PR:** 699,973 plasmids. The Genome Portal is being sunset mid-2025, so use the JGI Data Portal ([portal](https://genome.jgi.doe.gov/portal/IMG_PR/IMG_PR.home.html)).
- **Infernal 1.1 speed:** 0.6 s/Mb/query; about 1 h per bacterial genome against Rfam 11 on a quad-core ([PMC3810854](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3810854/)).
- **geNomad speed:** 242 s for 18.3 Mb, threads unstated ([PMC11324519](https://pmc.ncbi.nlm.nih.gov/articles/PMC11324519/table/Tab1)).

**Not verified, to check before S1:**
- Current INPHARED genome count.
- Total bp of PLSDB and IMG/VR-hc.
- GTDB taxonomy license.
- Whether R226 metadata carries CheckM2 fields.
- The Linclust throughput and memory formula.
- RefSeq `Complete Genome` bacterial/archaeal count (read from `assembly_summary_refseq.txt`).
- GTDB R226 complete-genome count, reps and members, and the ≤20-per-species selection (`ncbi_assembly_level` in `bac120/ar53_metadata_r226`). This pins S2C and is an S−1 task.
- AllTheBacteria precomputed annotation availability.
