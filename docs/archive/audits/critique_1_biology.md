# Audit Report 1: Comparative Genomics & Evolutionary Biology Critique of `all-the-data` (`PLAN.md` v1.0)

**Auditor:** Comparative Genomics & Evolutionary Biology Critic (Domain Lenses: *Eugene Koonin, L. Aravind, Valérie de Crécy-Lagard, Morgan Price, Ronald Breaker*)
**Target Documents Audited:** `00-spec.md`, `PRIOR-ART.md`, `PLAN.md` (v1.0 — Pre-Audit Draft)
**Overall Verdict:** The core vision of `all-the-data`—unifying sub-ORF domain/fold dissection, cis-regulatory RNAs, phylogenetic dispersion bitvectors, and experimental closure across Tiers 1–3—is a genuine paradigm shift over `STRING`, `AnnoTree`, `EFI-GNT`, and `Foldseek`. However, **`PLAN.md` v1.0 contains 10 severe biological and evolutionary modeling flaws** (`5 CRITICAL`, `4 MAJOR`, `1 MODERATE`) where simplifying computational assumptions collide with real prokaryotic, archaeal, and viral molecular evolution. Left unfixed, these flaws will systematically corrupt `feature_tape.parquet`, blind the NOGD bitvector engine to >80% of true pathway replacements, shatter riboswitch-regulated operons, and erase polymorphic toxin tips.

---

## Executive Summary Matrix of Biological Failure Modes

| ID | Audit Category | Specific Biological Failure Mode in `PLAN.md` v1.0 | Severity | Primary Affected Components |
|---|---|---|---|---|
| **BIO-01** | **Sub-ORF Domain Evolution** | Linear `cath-resolve-hits` N-to-C splitting shatters **discontinuous / nested domains** (`Cas9` `RuvC` split by `HNH`/`REC`; `AspRS` split by `GAD`; Radical SAMs) into fake tandem duplicates with `catalytic_intact = false` | **CRITICAL** | §3.2 Cascade 2, `feature_tape.parquet` |
| **BIO-02** | **Sub-ORF Domain Evolution** | Full-protein `NR90` pre-clustering (`90% ID / 90% cov`) before domain dissection **erases hyper-variable C-terminal toxin tips (`CTDs`)** in long polymorphic toxins (`Rhs`, `CdiA`, `VgrG`, `MafB`) and mis-pairs immunity proteins | **CRITICAL** | §3.2 `NR90` step, `directon_pair_summary` |
| **BIO-03** | **Sub-ORF Domain Evolution** | **Tandem repeat arrays** (`TPR`, `ANK`, `WD40`, `Pentapeptide`, `RHS` shell) expand into 15–45 consecutive `slot_idx` tokens per protein, breaking `±k` slot window searches (`Index 1`) and flooding `gLM-Domain` context windows | **MAJOR** | §3.2 Cascade 2, §4.2 Index 1, §4.5 Index 4 |
| **BIO-04** | **Sub-ORF Domain Evolution** | **Circular permutations (`CP`)** in methyltransferases, `PLD`/`PD-(D/E)xK` nucleases, and swap-domain effectors produce split partial HMM hits (`catalytic_intact = false`) or fake "dark domains" | **MAJOR** | §3.2 Cascade 2 & 4 |
| **BIO-05** | **Feature Tape & Gene Calling** | Unarbitrated overlap between **spurious `pyrodigal`/`smORF` calls and real structured RNAs (`Infernal`), CRISPR arrays (`Minced`), or `CreTA`/uORFs** injects fake antisense strand flips that sever `directon_id` between riboswitches and their target operons | **CRITICAL** | §3.1 Stage 1, §4.3 `feature_tape.parquet` |
| **BIO-06** | **Feature Tape & Gene Calling** | **Programmed ribosomal frameshifts (`prfB`, phage `gpG-gpGT`, IS `orfA-orfB`), selenoproteins (`UGA=Sec`), split inteins, and Group I/II self-splicing introns** shatter intact genes into pseudo-fragments separated by homing endonucleases | **MAJOR** | §3.1 Stage 1, §4.3 `feature_tape.parquet`, §4.4 NOGD |
| **BIO-07** | **Feature Tape & Gene Calling** | Hardcoded Shine-Dalgarno (`AGGAG`/`GGAGG`) requirement in the `smORF` scanner is **blind to leaderless transcripts** (40–70% of *Archaea* & *Actinomycetota*) and **SD-independent *Bacteroidota* / *Cyanobacteriota* initiation** | **MAJOR** | §3.1 Stage 1 (`smORF` scanner), §2.3 Gate |
| **BIO-08** | **NOGD & Phylo-Profiling** | Single-genome binary `presence_bits` across `113,104` GTDB species reps (65% MAGs at 90–93% completeness) causes **~8% assembly-dropout false holes per gene** and conflates **patchy plasmid/prophage MGEs (5% of strains)** with core chromosomes | **CRITICAL** | §4.4 Index 3 (`presence_bits`) |
| **BIO-09** | **NOGD & Phylo-Profiling** | Indexing `presence_bits` and `directon_pair_summary` on coarse `Pfam`/`COG` tokens instead of **isofunctional orthologous sub-clusters (`SSN` + active-site signature)** makes NOGD completely blind to multi-paralog families (`COG0720` `PTPS-I..VI`, `PF04055` Radical SAM, `SDR`) | **CRITICAL** | §4.1 `token_dictionary`, §4.3, §4.4 Index 3 |
| **BIO-10** | **Operon & Regulatory Topology** | Unidirectional `directon_id` (broken on any strand flip, `>300 bp` gap, or `TransTermHP` terminator) **severs `<-- -->` divergons (TetR/LysR regulators, R-M systems), `--> <--` antisense excludons, and riboswitch attenuator terminators** | **CRITICAL** | §4.3 `feature_tape.parquet`, `directon_pair_summary` |
| **BIO-11** | **Specificity & Deduplication** | `arch_hash` deduplication (`§2.3`) and `Infernal` family tokens collapse **riboswitch ligand-specificity switches (`ykkC` Guanidine-I vs. `ppGpp` vs. `PRPP`)** and **catalytically inactivated pseudo-enzyme scaffolds** into identical tokens | **MODERATE** | §2.3 Ingest Sieve, §3.1 `Infernal`, §5 Benchmark 1 |

---

## Detailed Audit Findings & Concrete v2.0 Engineering Fixes

---

### 1. Linear Sub-ORF Domain Splitting (`cath-resolve-hits`) vs. Real Protein Evolution

#### [BIO-01] Discontinuous / Nested ("Russian-Doll") Protein Domains Are Shattered or Dropped by Linear N-to-C Segmentation
- **Severity:** **CRITICAL**
- **Exact Biological Failure Mechanism & Real Examples:**
  In `00-spec.md` (§2 Invariant 2) and `PLAN.md` (§3.2 Cascade 2, line 75), `cath-resolve-hits` resolves HMM hits into strictly linear, non-overlapping N-to-C segments `(seg_1, seg_2, ..., seg_k)` placed sequentially onto `slot_idx`. In real protein evolution, newly recruited domains frequently **insert directly into a surface loop inside an existing catalytic fold**, creating discontinuous (nested) domain architectures:
  1. **CRISPR-Cas & Transposon RNA-Guided Nucleases (`Cas9`, `Cas12b`/`C2c1`, `IscB`, `TnpB`):** The `RuvC` (`RNase H`-fold) nuclease domain is split into three discontinuous sequence blocks (`RuvC-I`, `RuvC-II`, `RuvC-III`) interrupted by the α-helical recognition (`REC`) lobe and an inserted **`HNH` endonuclease domain (`PF01844`)** nested directly between `RuvC-II` and `RuvC-III`.
  2. **Class II Aminoacyl-tRNA Synthetases (`AspRS` `COG0173`, `LysRS`):** Prokaryotic AspRS contains a ~120 aa **GAD editing domain (`PF02938`)** nested inside the catalytic core (`PF00152`) between Motif 2 and Motif 3.
  3. **Metabolic TIM Barrels, Rossmann Folds, P-Loop NTPases, and Radical SAMs:** E.g., `RecBCD`/`AddAB` helicase-nucleases, `Lon` AAA+ protease, `RuBisCO`, and `B12`/SPASM-containing Radical SAM enzymes (`PF04055`) where auxiliary substrate-recognition or ferredoxin modules interrupt the core fold.
  - **How `PLAN.md` v1.0 breaks:**
    - If `cath-resolve-hits` runs in naive contiguous mode, it either discards the nested domain (`HNH` or `GAD`) because its bitscore is lower than the enclosing domain, or discards the enclosing domain.
    - If `cath-resolve-hits` splits the enclosing domain into two separate N-to-C segments around the insertion (`[RuvC_N] -> [HNH] -> [RuvC_C]`), `feature_tape.parquet` emits **two separate `slot_idx` rows with the same `seq_family_id` (`RuvC`)**. That creates a false tandem duplication signature (`RuvC + RuvC`), shifts all downstream `slot_idx` offsets, and causes the active-site checker (`PLAN.md` line 78) to evaluate `RuvC_N` (`D10` only) and `RuvC_C` (`E762/H983` only) separately—**falsely marking both halves `catalytic_intact = false`!**
- **Concrete Schema / Algorithmic Fix for `PLAN.md` v2.0:**
  1. **Discontinuous HMM Envelope Stitching in Stage 2 (`Cascade 2 & 3`):**
     - Run `pyhmmer` domain table extraction retaining multi-envelope alignments per target HMM. When two protein regions `[env_start_1..env_end_1]` and `[env_start_2..env_end_2]` separated by an intervening inserted domain `[ins_start..ins_end]` map monotonically in HMM profile coordinates (`hmm_to_1 < hmm_from_2 + 15 aa`), **stitch them into a single logical domain instance** before checking CDD/SFLD catalytic residues (`catalytic_intact` is evaluated on the stitched full alignment!).
  2. **Schema Update to `feature_tape.parquet` (§4.3 Table 1):**
     - Emit **one `slot_idx` row per logical domain** (preventing fake tandem split copies of the host domain).
     - Add `parent_seg_idx_in_orf: int8` (`-1` for top-level N-to-C domains; `>= 0` pointing to the enclosing host `seg_idx_in_orf` when a domain is nested inside a discontinuous host domain) and `is_discontinuous_host: bool`.
     - Add `nested_insertion_count: uint32` to `directon_pair_summary.parquet` so domain-in-domain nesting is explicitly queryable alongside N-to-C Rosetta-stone fusions!

---

#### [BIO-02] Whole-Protein `NR90` Pre-Clustering Erases Hyper-Variable C-Terminal Toxin Tips (`CTDs`) in Polymorphic Toxin Systems
- **Severity:** **CRITICAL**
- **Exact Biological Failure Mechanism & Real Examples:**
  In `PLAN.md` §3.2 (line 65), *"Every called protein is first mapped to the Non-Redundant Sequence Dictionary (`NR90`, `MMseqs2 linclust` at 90% ID / 90% coverage) so expensive profile and structure tools run only once per unique sequence cluster."*
  - In bacterial biological conflict (`Zhang, de Souza, Anantharaman, Iyer, Aravind 2012, Biol Direct`), **polymorphic toxins** (`Rhs`/`YD`-repeat toxins, `CdiA` contact-dependent inhibition proteins, `T6SS` evolved `VgrG`/`PAAR`, `T7SS` `LXG`, `MafB`, `FHA` filamentous hemagglutinins) consist of a huge conserved N-terminal delivery shaft (`1,200–3,500 aa`), a `pre-tox` cleavage motif (`PT-hint`, `VENN`, `RHS` core-box `PxxxxDPxGL`), and a **hyper-variable C-terminal toxin tip (`CTD`, `60–130 aa`)** that swaps rapidly by recombination, followed immediately on the chromosome by its **cognate immunity gene (`Imm`)** and often a downstream array of **orphan `[CTD_i - Imm_i]` replacement cassettes**.
  - **Consider the exact arithmetic:** Two `2,500 aa` `Rhs` or `CdiA` proteins in different strains share an identical `2,390 aa` N-terminal shaft (`99.8%` identical) but carry **completely unrelated `110 aa` C-terminal toxin tips** (e.g., an `EndoU` RNase tip vs. a `deaminase` tip). Across the full `2,500 aa` length, these two proteins have **$2390 / 2500 = 95.6\%$ sequence identity and $95.6\%$ alignment coverage**!
  - `MMseqs2 linclust --min-seq-id 0.9 -c 0.9` **collapses both proteins into a single `nr90_protein_id` representative**!
  - **Consequences:**
    1. All alternative C-terminal toxin tips sharing a conserved shaft are **permanently deleted before domain dissection**!
    2. In `feature_tape.parquet`, Genome 2's `Rhs` protein inherits Genome 1's `EndoU` toxin tip token, which is now placed immediately upstream of Genome 2's **`deaminase` immunity protein (`Imm_2`)**—**corrupting `directon_pair_summary` with false toxin-immunity pairings!**
    3. Furthermore, when a `2,500 aa` `Rhs` protein is passed to `ProstT5` in Cascade 4 (`PLAN.md` line 87), standard `ProstT5` token windows (`1,024` or `2,048` aa) truncate or dilute the C-terminal `100 aa` tip!
- **Concrete Schema / Algorithmic Fix for `PLAN.md` v2.0:**
  1. **Terminal-Window Guard on `NR90` Clustering (§3.2):**
     - When clustering proteins $> 350\text{ aa}$ into `NR90`, require **both** $\ge 90\%$ global sequence identity/coverage **AND $\ge 90\%$ sequence identity specifically across the terminal `120 aa` N-terminal and `120 aa` C-terminal windows** (implemented via a fast 64-bit rolling k-mer sketch / alignment check on `seq[-120:]` and `seq[:120]`). Two `Rhs`/`CdiA` proteins with different `CTD` tips are never collapsed into the same `nr90_protein_id`.
  2. **Dedicated Polymorphic Toxin `pre-tox` Slicer & Orphan Cassette Detector (§3.2 Cascade 2–4):**
     - Whenever a protein contains a delivery shaft (`RHS`, `VgrG`, `PAAR`, `LXG`, `WXG100`, `FHA`, `MafN`, `SpvB`) or a `pre-tox` motif (`PT-hint`, `VENN`, `DPxGL`), explicitly slice the downstream C-terminal segment (`35–180 aa` to C-terminus) as an independent segment (`seg_subkind = CTD_TOXIN_TIP`) and run `hhsearch` (Cascade 3) and `ProstT5`/`Foldseek` (Cascade 4) **on the isolated CTD sequence alone** (not the whole `2,500 aa` protein!).
     - Flag downstream short ORFs (`50–200 aa`) containing a `pre-tox` linker + `CTD` without an N-terminal delivery shaft as **`ORPHAN_CTD_CASSETTE`** so they are not mistaken for pseudogene assembly fragments.

---

#### [BIO-03] Tandem Repeat Arrays (`TPR`, `ANK`, `WD40`, `LRR`, `Pentapeptide`) Flood `slot_idx` and Break Operon Distance Windows
- **Severity:** **MAJOR**
- **Exact Biological Failure Mechanism & Real Examples:**
  In `Pfam` and `NCBIfam`, repeat-unit HMMs (`PF00515` `TPR_1`, `PF00023` `Ank`, `PF00400` `WD40`, `PF00805` `Pentapeptide`, `PF05593` `RHS_repeat`, `PF00560` `LRR_1`, `PF01839` `FG-GAP`) match single `25–42 aa` structural repeats rather than the entire solenoid or β-propeller domain.
  - A single `STAND` NTPase sensor (`Avs`, `bNACHT`, `MalT`-family regulator) or cyanobacterial signaling protein routinely contains **15 to 40 tandem `TPR`, `WD40`, or `ANK` repeat hits**.
  - Under `PLAN.md` §3.2 and §4.3, `cath-resolve-hits` resolves these into 15–40 separate non-overlapping segments `(seg_1, ..., seg_35)`, and `feature_tape.parquet` increments `slot_idx` for every segment!
  - **How `PLAN.md` v1.0 breaks:**
    1. **Index 1 (`contig_slot_offsets.bin`, §4.2 line 123):** Gapped operon queries check slot distance $|s_A - s_B| \le k$ (typically $k = 4\text{–}6$). A single `TPR`-containing protein sitting between Gene A and Gene B inflates $|s_A - s_B|$ from `2` to `36`, causing Index 1 to miss the operon!
    2. **Index 4 (`gLM-Domain`, §4.5 line 191):** Trained on windows of `4 to 24 slots long`. A single repeat protein fills the entire 24-slot transformer window with `[TPR, TPR, TPR, ..., TPR]`, blinding the language model to the surrounding operon!
    3. **Cascade 4 (`ProstT5` / `Foldseek`):** Individual `30 aa` repeat units are too short for meaningful `Foldseek` `ECOD` topology assignment, whereas the complete `12-repeat` solenoid or `7-blade` WD40 β-propeller maps cleanly to an `ECOD` `T-group`.
- **Concrete Schema / Algorithmic Fix for `PLAN.md` v2.0:**
  1. **Homotypic Repeat-Array Super-Segment Collapsing in Stage 2 (`§3.2`):**
     - Define a repeat-clan registry (`CL0020` TPR, `CL0465` ANK, `CL0186` WD40/β-propeller, `CL0022` LRR, `CL0193` Pentapeptide, `RHS`/`YD` shell, `SLH`, `Dockerin`/`Cohesin`).
     - Collapse $N \ge 2$ consecutive hits belonging to the same repeat clan (separated by $\le 35\text{ aa}$) into **a single `slot_idx` / `seg_idx_in_orf` super-segment** spanning from the start of Repeat 1 to the end of Repeat $N$.
     - Store `repeat_unit_count: uint8` in the segment's metadata (or replace 6 bits of `length_aa_bin` / add a `repeat_count` field in `feature_tape.parquet`), and pass the **complete collapsed solenoid/propeller region** to `Foldseek` (`Cascade 4`).
  2. **Index 1 Dual Coordinate Addressing (`§4.2`):**
     - Store each entry in `contig_slot_offsets.bin` as a packed `uint32` encoding `(orf_idx: uint16, seg_idx_in_orf: uint8, strand_and_kind: uint8)` so queries can filter on **both** inter-gene distance $|\text{orf}_A - \text{orf}_B| \le k_{\text{orf}}$ and intra-ORF fusion distance $|\text{seg}_A - \text{seg}_B| = 1$ without any interference from multi-domain complexity!

---

#### [BIO-04] Circular Permutations (`CP`) in Methyltransferases, Nucleases, and Effectors Produce Split False-Negative Hits
- **Severity:** **MAJOR**
- **Exact Biological Failure Mechanism & Real Examples:**
  Circular permutations—where the canonical N-terminal half `[A]` and C-terminal half `[B]` of a domain are rearranged into `[B] -> [A]` via tandem duplication and terminal trimming—are widespread across prokaryotic defense and modification enzymes:
  - **Class I Rossmann-fold DNA/RNA Methyltransferases (`PF01555`, `PF13649`, `BREX` `PglX`, Type I/II/III R-M MTases):** Frequently permute the relative order of the SAM-binding motif (`Motif I`) and catalytic target-base motif (`Motif IV`, `NPPY`/`DPPY`).
  - **Restriction-like `PD-(D/E)xK` Endonucleases, `PLD` Nucleases (`Cardiolipin synthase` / phage defense nucleases), `GIY-YIG`, and `Cas` Effectors.**
  - **How `PLAN.md` v1.0 breaks:** Linear HMM alignment (`pyhmmer`) aligns `[B]` at the protein's N-terminus (`hmm_coords = [130..240]`) and emits a second separate hit for `[A]` at the protein's C-terminus (`hmm_coords = [1..129]`). `cath-resolve-hits` either emits two partial domains of the same family (both failing `catalytic_intact`!) or drops the shorter half into Cascade 3 (`>= 45 aa` dark segment), hallucinating a fake "novel dark domain" fused to a truncated enzyme.
- **Concrete Schema / Algorithmic Fix for `PLAN.md` v2.0:**
  - Add a **Circular-Permutation (`CP`) Stitcher** immediately after `pyhmmer`/`hhsearch` envelope extraction:
    - When two adjacent hits in the same ORF map to the same `seq_family_id` or ` clan_or_cog_id` with inverted HMM coordinates—Hit 1 covering `[hmm_mid..hmm_L]` at protein `[p1..p2]` and Hit 2 covering `[1..hmm_mid]` at protein `[p3..p4]` ($p_3 > p_2$, HMM overlap $< 15\%$, combined HMM coverage $\ge 70\%$)—merge them into **one domain segment**, evaluate `catalytic_intact` on the combined alignment, and set `is_circularly_permuted: bool = true`.

---

### 2. Gene Calling, `smORF`, & Non-Coding Collisions on the Feature Tape

#### [BIO-05] Spurious ORF/`smORF` Collisions with Riboswitches (`Infernal`), CRISPR Arrays (`Minced`), and `CreTA` Corrupt `slot_idx` and Sever `directon_id`
- **Severity:** **CRITICAL**
- **Exact Biological Failure Mechanism & Real Examples:**
  In `PLAN.md` §3.1, four independent feature callers run on every contig: `pyrodigal` (CDS), the intergenic `smORF` scanner (`25–75 aa`), `Infernal cmscan` (`Rfam` + 550 orphan riboswitches/ncRNAs), and `Minced`/`IntegronFinder`/`tRNAscan-SE`. All called features are then placed sequentially onto `feature_tape.parquet` (`slot_idx = 0, 1, 2, ...`) where `directon_id` increments whenever `strand` flips (`PLAN.md` line 142).
  - **Three catastrophic collision modes occur because `PLAN.md` v1.0 lacks an overlap-arbitration hierarchy:**
    1. **Spurious Sense/Antisense ORFs Inside Riboswitches & Large Structured RNAs (`Breaker / Weinberg`):**
       High-GC structured RNAs—such as `Cobalamin` riboswitches (`RF00174`, ~210 nt), `T-box` leaders (~280 nt), `Lysine` / `Glycine` / `SAM-I/II` riboswitches, `OLE` RNA (`RF01071`, ~610 nt), `GOLLD` (~800 nt), `HEARO` (~350 nt), `RNase P`, and `tmRNA`—routinely contain G-rich base-paired stems that match the `smORF` scanner's Shine-Dalgarno seed (`AGGAG`/`GGAGG`, `PLAN.md` line 55) followed by an `ATG`/`GTG`/`TTG` and a `25–90 aa` open reading frame on either the `+1` or `-1` strand.
       - If a spurious `35 aa` ORF on `strand = -1` is called inside a `strand = +1` riboswitch sitting upstream of its regulated transporter (`strand = +1`), the coordinate-ordered tape becomes:
         `[Riboswitch (+1)] -> [Spurious_smORF (-1)] -> [Transporter (+1)]`.
       - Because `strand` flips (`+1 -> -1 -> +1`), **`directon_id` increments twice**, severing the riboswitch from the operon it regulates, corrupting `delta_bp_prev_orf`, and inserting a hallucinated `DARK_FAM` protein between the riboswitch and its target gene!
    2. **Spurious ORFs Across CRISPR Repeat-Spacer Arrays (`Minced`):**
       In ~33% of CRISPR arrays, the repeat+spacer periodicity (`e.g., 30 bp repeat + 36 bp spacer = 66 bp = 22 codons`) is a multiple of 3 nucleotides, creating an unbroken reading frame across the entire CRISPR array that `pyrodigal --meta` calls as a repetitive protein!
    3. **True Biologically Overlapping RNA + Peptide Elements (`CreTA`, Attenuator Leader `uORFs`, Retron `msr-msd`):**
       Conversely, deleting *all* ORFs that overlap non-coding RNAs would destroy **real dual-function loci**:
       - **Leader Attenuator Peptides (`trpL`, `hisL`, `leuL`, `ilvB`, `pheL`, `pyrL`, `ermC`, `speF`):** A true `14–55 aa` regulatory `uORF` is translated *inside* the 5' cis-regulatory attenuator RNA!
       - **`CreTA` Safeguard Modules (`Li / Koonin 2021 Science`):** The mini-toxin `CreT` (`15–50 aa` tRNA-sequestering ORF) sits immediately inside or adjacent to the `CreA`-flanked CRISPR-like RNA cassette.
       - **Retron `msr-msd` ncRNAs:** Often overlap the 5' end of the retron reverse transcriptase (`RT`) or encode a nested effector/toxin `smORF`.
- **Concrete Schema / Algorithmic Fix for `PLAN.md` v2.0:**
  1. **Add Stage 1.4: Coordinate Overlap Arbitration & Hierarchy Rules (`§3.1`):**
     - **Rule 1 (Hard Non-Coding Mask):** Any `pyrodigal` ORF or `smORF` overlapping $\ge 25\%$ of its length with a **CRISPR array (`Minced`), tRNA (`tRNAscan-SE`), rRNA, `RNase P`, `SRP` RNA, or `Integron` `attC` array** that lacks a high-confidence known protein HMM hit (`E < 1e-10`, excluding low-complexity repeat models) is **purged before `NR90` and tape construction**.
     - **Rule 2 (Riboswitch / Orphan ncRNA vs. ORF Arbitration):**
       - When a called ORF/`smORF` overlaps an `Infernal` covariance hit (`E < 1e-5` or bitscore > gathering cutoff):
         - **Antisense (`strand_orf != strand_rna`):** Purge the ORF unless it has a validated `Pfam`/`COG`/`PHROG`/`AFDB` hit (`E < 1e-10`) AND is $> 90\text{ aa}$ (or matches a known Type I TA toxin profile like `hok`/`tisB`/`symE`/`txpA`, in which case retain both as a **`SENSE_ANTISENSE_TA_PAIR`**!).
         - **Sense (`strand_orf == strand_rna`):** If the ORF is $\le 75\text{ aa}$ and nested within the `Infernal` 5'-leader region, classify it as `feature_kind = 8` (**`uORF_cis_leader`**) or `CreT` mini-ORF, **attach it as a child of the riboswitch/leader (`parent_slot_offset = -1`)**, and **compute `delta_bp_prev_orf` for the downstream structural gene from the 3' end of the `Infernal` RNA element**, so an embedded `uORF` never shifts or breaks the riboswitch-to-operon linkage!

---

#### [BIO-06] Programmed Ribosomal Frameshifts (`prfB`, `gpG-gpGT`, IS Elements), Selenoproteins (`Sec`), Split Inteins, and Group I/II Self-Splicing Introns Shatter Intact Proteins
- **Severity:** **MAJOR**
- **Exact Biological Failure Mechanism & Real Examples:**
  In `PLAN.md` §3.1, `pyrodigal` assumes every protein is encoded by a single uninterrupted open reading frame from start codon to stop codon (plus whole-contig genetic code switching for codes `4, 15, 25`). In real bacterial, archaeal, and especially **Tier 2 phage (`IMG/VR`)** genomes, intact functional proteins are frequently split at the DNA level by four mechanisms:
  1. **Self-Splicing Group I & Group II Introns & Inteins (Rampant in Phages, Thermophiles, Cyanobacteria, and Archaea):**
     - Core phage and cellular replication/metabolism genes—including ribonucleotide reductase (`nrdA`, `nrdB`, `nrdD`), thymidylate synthase (`td`/`thyA`), DNA polymerase (`polA`/`gp43`), large terminase (`terL`), `recA`, and cyanobacterial `dnaE` (split intein)—are routinely interrupted by $300\text{–}2,500\text{ bp}$ Group I (`RF00028`) or Group II (`RF00029`) introns containing an internal **homing endonuclease (`LAGLIDADG` `PF00961`, `GIY-YIG` `PF01541`, `HNH` `PF01844`)** or **intron maturase/RT (`PF00078`)**, or by in-frame/split **inteins (`HINT` `PF00665`)**.
     - `pyrodigal` calls `[nrdA_exon1] -> [HNH_endonuclease] -> [nrdA_exon2]`.
     - **Failure cascade:** `nrdA` is shattered into two truncated ORFs separated by $>1,000\text{ bp}$ (breaking `directon_id`!), both exons fail `catalytic_intact` (`false`), `directon_pair_summary` learns a spurious functional association between `nrdA` and `HNH` nucleases, and in Tier 1 genomes where exons are short, the genome is falsely scored as **lacking `nrdA` in NOGD `presence_bits`**!
  2. **Programmed `-1` / `+1` Ribosomal Frameshifts (`PRF`):**
     - **Release Factor 2 (`prfB`, `COG1186`):** Uses a conserved `+1` frameshift (`CUU UGA C`) across ~70% of bacteria. `pyrodigal` misses the N-terminal 25–30 aa fragment (`PF02978`), calls a truncated ORF, and miscalculates $\Delta\text{bp}$ to the upstream gene (`lysS`/`herC`).
     - **Tailed Phage Tail Assembly Chaperones (`gpG-gpGT`, `PF06199`):** Conserved across >70% of *Caudoviricetes* (`Tier 2`); uses a `-1` slippery heptamer to fuse `gpG` and `gpGT`.
     - **IS Transposases (`IS1` `insA-insB`, `IS3` `orfA-orfB`):** Use a `-1` frameshift (`AAAAAAC`) to fuse the HTH DNA-binding domain (`orfA`) with the `DDE` catalytic domain (`orfB`).
  3. **Selenocysteine (`UGA = Sec`) & Pyrrolysine (`UAG = Pyl`) Readthrough:**
     - Whole-contig genetic code switching (`PLAN.md` line 52) **cannot detect `Sec` or `Pyl`** because only 1–5 genes in the genome (`fdhF`, `fdoG`, `grdA/B`, `selD`, `mtmB`/`mtbB`/`mttB`) read through a single internal in-frame `TGA` (`Sec`) or `TAG` (`Pyl`) codon while the rest of the genome uses `TGA`/`TAG` as standard stops! `pyrodigal` cuts every selenoprotein in half right at the catalytic `Sec` residue!
- **Concrete Schema / Algorithmic Fix for `PLAN.md` v2.0:**
  - **Add Stage 1.5: Co-Linear Split-Gene, Frameshift, Intron, and `Sec`/`Pyl` Reconstructor:**
    1. **Selenocysteine / Pyrrolysine / `+1`/`-1` Frameshift Stitcher:**
       - Scan adjacent same-strand ORF pairs `(orf_i, orf_{i+1})` with $-60\text{ bp} \le \Delta\text{bp} \le +90\text{ bp}$ (including pairs separated by a single in-frame `TGA`/`TAG` codon or overlapping in `+1`/`-1` frames) against `pyhmmer`.
       - When `orf_i` aligns to `[1..m]` and `orf_{i+1}` aligns to `[m+1..L]` of the **same single-polypeptide HMM (`COG`/`NCBIfam`/`Pfam`/`PHROG`)** (`prfB`, `fdhF`, `selD`, `gpG-gpGT`, `IS3`, etc.), reconstruct the merged full-length virtual protein sequence (`U` at in-frame `TGA`, `O` at `TAG`, or frameshift-stitched AA sequence) for `NR90` and `catalytic_intact` scoring, and merge them into a single `orf_idx` with a flag bit in `biophysical_badge` (`bit 40: has_frameshift_or_sec_readthrough`, `bit 41: has_spliced_intron_or_intein`).
    2. **Group I/II Intron & Split-Intein Exon Bridge:**
       - When two same-strand ORFs within $\le 3,500\text{ bp}$ (separated only by a Group I/II intron RNA `RF00028`/`RF00029`, an `HNH`/`LAGLIDADG`/`GIY-YIG` homing endonuclease, an intron maturase `PF00078`, or a `HINT` intein `PF00665`) map monotonically to consecutive N- and C-terminal halves `[1..m]` and `[m+1..L]` of the same host HMM (`nrdA/B/D`, `td`, `polA`, `terL`, `recA`, `dnaE`), stitch the host exons for `NR90`, `catalytic_intact`, and `presence_bits`, keep both exons in the same `directon_id`, and flag the intervening HEG/intron slot with `is_intronic_nested = true` so it does not pollute host operon synteny!

---

#### [BIO-07] Canonical Shine-Dalgarno (`AGGAG`/`GGAGG`) Requirement in `smORF` Scanner Drops Leaderless (*Archaea*/*Actinomycetota*) and SD-Independent (*Bacteroidota*/*Cyanobacteriota*) `smORFs`
- **Severity:** **MAJOR**
- **Exact Biological Failure Mechanism & Real Examples:**
  In `PLAN.md` §3.1 (line 55), the intergenic `smORF` scanner (`25–75 aa`) requires `ATG`/`GTG`/`TTG` to be *"preceded within 4–14 bp by a canonical ribosome-binding site (`AGGAG` / `GGAGG`) or within 0–5 bp of an upstream stop codon (translational coupling `ATGA` / `TGATG`)."*
  - While this works for interior translationally coupled `smORFs`, **it fails for operon-initial (first-in-operon) and standalone `smORFs`** across four huge branches of microbial life:
    1. **Archaea (*Halobacteriota*, *Thermoproteota*, *Methanobacteriota*, *Asgardarchaeota*) & *Actinomycetota* (*Mycobacterium*, *Streptomyces*, *Corynebacterium*):** **40% to 70% of all transcripts are leaderless mRNAs** where transcription initiates directly at the `+1` `ATG`/`GTG` start codon (driven by an archaeal `TATA` box + `BRE` at `-26..-30 bp` or an actinobacterial `-10` `TANNNT` hexamer at `-12..-7 bp`) with **zero 5'-UTR and zero Shine-Dalgarno motif**!
    2. ***Bacteroidota* (*Bacteroides*, *Flavobacterium*, *Cytophaga*) & *Cyanobacteriota*:** Ribosomes in *Bacteroidota* have their anti-SD 16S rRNA tail sequestered by ribosomal protein `bS21`/`bS1` and **do not use `AGGAG`/`GGAGG` RBSs at all**—instead, translation initiation relies on an **AT-rich, unstructured `-20..-4 bp` upstream window** (`Accetto / Fredrick / Mankin 2020–2023`)!
    3. **Eukaryotic-Host / Metatranscriptomic Viral Elements (`Obelisks`, RNA viruses in Tier 2/3):** Lack bacterial `AGGAG` motifs; single-ORF or two-ORF `Obelisks` (~1,000 nt circular RNAs) will fail both the `smORF` scanner and the `>= 2-Feature Gate` (`§2.3` line 41) if `Oblin-2` (`~50 aa`) is missed!
- **Concrete Schema / Algorithmic Fix for `PLAN.md` v2.0:**
  1. **Use `pyrodigal`'s Per-Contig Learned Initiation Table + 4-Mode Initiation Scoring in the `smORF` Scanner (`§3.1`):**
     - `pyrodigal` already infers whether a contig uses Shine-Dalgarno, leaderless, or AT-rich non-SD initiation (`training_info.uses_sd`). Reuse `pyrodigal`'s learned translation initiation weight matrix down to `--min-gene 60` (`20 aa`) in intergenic regions, accepting an operon-initial `smORF` (`20–75 aa`) if it satisfies **any** of:
       - **(a) Canonical SD:** `AGGAG` / `GGAGG` / $\ge 4\text{-nt}$ 16S anti-SD match in `-15..-4 bp`,
       - **(b) Translational Coupling:** `-11 to +5 bp` from an upstream same-strand stop codon,
       - **(c) Leaderless Promoter Match:** `ATG`/`GTG` located `6–11 bp` downstream of a `-10` `TANNNT` box (*Actinomycetota*/*Deinococcota*) or `22–32 bp` downstream of an archaeal `TATA` box (`TTTAWA`/`TTTAAA`),
       - **(d) *Bacteroidota* Ribosomal `bS1` Window:** $\ge 78\%$ `A/T` content and predicted $\Delta G_{\text{fold}} > -1.5\text{ kcal/mol}$ in `-18..-4 bp` on contigs classified as *Bacteroidota* / *Cyanobacteriota* (`uses_sd == false`),
       - **(e) Profile / Biophysical / Comparative Rescue:** Any `18–75 aa` ORF matching a known HMM (`TADB`, `AcrDB`, `RiPP` precursor leader, `CreT`, `PHROGs`), containing a predicted `TM` helix (`TMbed`) or `CxxC..CxxC` metal-binding motif, or conserved across $\ge 3$ contigs with $d_N/d_S < 0.5$.
  2. **Whitelist `Obelisks` & Viroid-like Elements in the `§2.3` 2-Feature Gate:**
     - Whitelist any circular or RNA-virome contig $\ge 750\text{ bp}$ encoding an `Oblin`-family hit or `hammerhead`/`HDV`/`twister` self-cleaving ribozyme (`Infernal`), even if it contains only 1–2 short dark ORFs.

---

### 3. NOGD & Phylogenetic Profiling Failure Modes

#### [BIO-08] Single-Genome Species Representative Bitvectors (`113,104` GTDB Reps) Collapse Under 90%-Complete MAG Dropouts and Patchy MGEs
- **Severity:** **CRITICAL**
- **Exact Biological Failure Mechanism & Real Examples:**
  In `PLAN.md` §4.4 (`Index 3: Tree-Ordered Phylogenetic Bitvectors`), each `token_id` is represented by a single `113,104`-bit vector (`presence_bits`) over the `113,104` GTDB r220 species representatives, and NOGD pathway holes are computed via:
  $$\text{Hole\_Genomes} = \text{bits}(A) \ \& \ \text{bits}(B) \ \& \ \text{bits}(D) \ \& \ \sim\text{bits}(C)$$
  - **Two fatal biological confounds invalidate this equation on raw GTDB species reps:**
    1. **The "90% Complete MAG = 8% Missing Core Genes" Math Catastrophe:**
       Roughly 65% of the `113,104` GTDB r220 species representatives are **unclosed metagenome-assembled genomes (MAGs) or draft assemblies** with completeness between `90%` and `95%` (`Tier 1` admits all genomes $\ge 90\%$ complete, `PLAN.md` line 29) split across `50–250` contigs.
       - Consider a universal 4-step biosynthetic pathway $(A, B, C, D)$ present in `20,000` species where $A, B, D$ sit in one locus and $C$ sits elsewhere (or at the end of the operon).
       - In `13,000` MAG species reps with average completeness `92%`, gene $C$ falls into an unsequenced assembly gap or contig break in **$8\%$ of MAGs $\approx 1,040$ species representatives**!
       - Meanwhile, true evolutionary Non-Orthologous Gene Displacement (NOGD) of $C$ might occur in only `150` species!
       - Thus, **$>87\%$ ($1,040 / 1,190$) of the bits set to `1` in $\text{Hole\_Genomes}$ are false-positive assembly dropouts in incomplete MAGs, NOT evolutionary gene displacements!**
       - Worse, because assembly gaps are distributed randomly across the phylogenetic tree, Fitch parsimony (`PLAN.md` line 184) scores those `1,040` MAG dropouts as **1,040 independent evolutionary loss events ($H_C = 1,040$)**!
    2. **Plasmid / Prophage / ICE Accessory Genes (`5%` of Strains) vs. Core Chromosome Conflation:**
       - If a gene $X$ lives on a conjugative plasmid, ICE, or prophage present in only 5% of *E. coli* strains:
         - If the single GTDB *E. coli* species representative happens to lack that plasmid, `presence_bits[E_coli] = 0` (even though 5,000 sequenced *E. coli* genomes in Tier 1/2 carry it!).
         - Conversely, if the single species rep happens to carry a prophage encoding `ThyX` or a defense system, `presence_bits` marks the entire species `1`, even though it is a transient MGE absent from 95% of the species!
- **Concrete Schema / Algorithmic Fix for `PLAN.md` v2.0:**
  1. **Replace the Single Binary Bitvector with a 3-Plane Pangenome-Aggregated Bitvector Tuple per Token (`§4.4 Index 3`):**
     Aggregate across **all `596,859` GTDB QC genomes + RefSeq Complete + `AllTheBacteria` species clusters + host-assigned Tier 2 plasmids/phages** for each of the `113,104` species leaves (`3 × 14,138 bytes = 42.4 KB` per token):
     - **`core_chr_bits` (`14,138 B`):** `1` iff the token is present on the chromosome in $\ge 50\%$ of high-completeness genomes of species $s$ (or present in the species rep and not flagged as MGE).
     - **`pan_any_bits` (`14,138 B`):** `1` iff the token is present in **any** genome assembly or host-linked plasmid/prophage of species $s$ (prevents false `~bits(C)` calls when the single species rep has a MAG assembly gap or lost a plasmid, while a sister strain of the same species has `C`!).
     - **`mge_accessory_bits` (`14,138 B`):** `1` iff the token occurs in $< 30\%$ of strains of species $s$ or resides on a Tier 2 plasmid/phage/integron/ICE.
  2. **Global Assembly-Quality & Contig-Edge Masks for NOGD Queries (`§4.4`):**
     - Pre-compute two global bitvectors (`14,138 B` each):
       - `hq_closed_or_multi_strain_mask`: `1` iff species $s$ either has a **closed complete chromosome** (`RefSeq Complete`), OR has $\ge 3$ independent sequenced genomes/MAGs (so `pan_any_bits` has $< 0.08^3 = 0.05\%$ chance of random assembly dropout!), OR has completeness $\ge 97\%$ with $\le 15$ contigs.
       - Per-pathway `operon_not_at_contig_edge_bits`: `0` if any gene of $(A, B, D)$ has `is_contig_edge == true` in species $s$.
     - **Updated NOGD Pathway-Hole Equation (0% MAG-Dropout Corruption):**
       $$\text{Verified\_Hole} = \text{core\_chr}(A) \ \& \ \text{core\_chr}(B) \ \& \ \text{core\_chr}(D) \ \& \ \sim\text{pan\_any}(C) \ \& \ \text{hq\_closed\_or\_multi\_strain\_mask} \ \& \ \text{operon\_not\_at\_contig\_edge}$$
     - Additionally, require that any candidate NOGD clade without a closed chromosome shows **phylogenetic coherence** (the hole $\sim\text{pan\_any}(C)$ is shared by $\ge 2$ consecutive leaves in the same genus/family on the DFS-ordered bitvector, tested in 1 clock cycle via bitwise shift `hole & (hole >> 1)`!).

---

#### [BIO-09] Multi-Paralog Superfamilies (`COG0720` PTPS, `PF04055` Radical SAM, `SDR`, `ABC`) Blind NOGD Unless Sub-Partitioned into Isofunctional Orthogroups Across All Indices
- **Severity:** **CRITICAL**
- **Exact Biological Failure Mechanism & Real Examples:**
  In `PLAN.md` §4.1–4.4, `token_dictionary.parquet`, `feature_tape.parquet`, `directon_pair_summary.parquet`, and `Index 3` (`presence_bits`) are indexed by `seq_family_id` (`Pfam` / `COG` / `arCOG` / `NCBIfam` / `PHROG`). While **Benchmark 4 (`PLAN.md` line 222)** briefly mentions splitting families into 40%-identity `subfamily_token_id`s for `family_closure_card`, **`subfamily_token_id` is completely absent from `feature_tape.parquet` (§4.3), `directon_pair_summary.parquet` (§4.3), `Index 1` (`token_to_contigs.roar`), and `Index 3` (`presence_bits` NOGD)!**
  - **Why NOGD and Synteny fail without isofunctional orthogroup sub-tokens:**
    1. **The `COG0720` (`PTPS` T-fold) Trap (`de Crécy-Lagard`):** `COG0720` / `PF01242` encompasses at least six distinct isofunctional enzyme families: `PTPS-I` (`QueD`, queuosine synthesis, present in >85% of bacteria), `PTPS-II` (biopterin synthesis), `PTPS-III` (the non-orthologous replacement for `FolB` in folate biosynthesis!), and `PTPS-IV..VI`. Because `QueD` (`PTPS-I`) is in almost all bacteria, the coarse family bitvector `bits(COG0720)` is **`1` in >85% of all species**! When searching for the NOGD replacement of `FolB`, `COG0720` shows **zero mutual exclusivity with `FolB`** at the family level because `FolB+` genomes still have `COG0720` (`QueD`)! Only the sub-cluster **`COG0720_PTPS_III`** is mutually exclusive with `FolB` and syntenic with `folKP`!
    2. **Mega-Families (`PF04055` Radical SAM, `PF00106` SDR Dehydrogenases, `PF00583` GNAT, `PF00005` ABC ATPases, `COG0500` SAM MTases, `PF13419` HAD Hydrolases):** Almost every prokaryotic genome encodes `5–60` paralogs of these families. Their family-level `presence_bits` is a solid wall of `113,104` ones (`100%` of species), making family-level phylogenetic profiling and NOGD algebra mathematically useless for the very enzyme superfamilies that drive >60% of metabolic NOGDs!
    3. **Asymmetric / Multi-Solution NOGDs (`MetE`/`MetH`, `ThyA`/`ThyX`, `NrdAB`/`NrdDG`):** `PLAN.md` line 183 searches for candidates $X$ with global mutual exclusivity against $C$. In reality: (a) a single pathway hole $\sim C$ across the tree of life is often filled by **3–5 different non-homologous displacements $X_1, X_2, X_3$** in different phyla (so each $X_i$ covers only 15–40% of $\text{Hole\_Genomes}$), and (b) facultative organisms often retain **both** $C$ and $X$ as condition-specific isoenzymes (e.g., $B_{12}$-dependent `MetH` + $B_{12}$-independent `MetE`).
- **Concrete Schema / Algorithmic Fix for `PLAN.md` v2.0:**
  1. **Two-Tier Token Hierarchy (`family_token_id` + `orthogroup_sub_id`) Across All 4 Indices (`§4.1–4.5`):**
     - Partition every `Pfam`/`COG`/`arCOG`/`PHROG`/`ECOD` family with $> 250$ sequences or mean copy number $> 1.15\text{/genome}$ into **Isofunctional Orthologous Sub-Clusters (`orthogroup_sub_id: uint32`)** using a combined **Sequence + Active-Site Pocket + Operon-Context** splitter:
       - Step A: Cluster domain sequences at `40% AA identity / 80% coverage` (`MMseqs2 cluster`).
       - Step B: Sub-split any 40%-ID cluster whose members diverge at **CDD/SFLD/AlphaFill active-site/ligand-pocket positions (`pocket_signature_hash`)** (cleanly separating `PTPS-I`, `PTPS-II`, `PTPS-III`, etc. even when overall sequence identity is ~35–40%).
       - Step C: Sub-split any cluster whose members stably reside in $\ge 2$ non-overlapping conserved operon contexts (`flanking_family_jaccard < 0.1` across $\ge 5$ orders).
     - Add `orthogroup_sub_id: uint32` as a primary indexed column in `feature_tape.parquet`, `directon_pair_summary.parquet`, `Index 1` (`token_to_contigs.roar`), `Index 3` (`presence_bits`), and `family_closure_card.parquet`.
  2. **Clade-Local Asymmetric NOGD Scoring Function (`§4.4`):**
     - Score candidate replacements $X$ (`family_token_id` or `orthogroup_sub_id`) not by global Jaccard over all holes, but by **Synteny-Seeded Independent Clade Displacement**:
       - First, identify all candidates $X$ that physically sit inside or adjacent to the $(A, B, D)$ locus (`directon` or `divergon`) in **at least 2 independent taxonomic orders** lacking $C$ (`Verified_Hole = 1`).
       - Second, score $X$ by the Fitch parsimony count of **independent evolutionary clades where gain of $X$ coincides with loss of $C$ ($H_{X^+ \land C^-}$)**, allowing partial coverage of global holes (so $X_1, X_2, X_3$ in different phyla are each recovered at the top of the list!) and tolerating $[C^+, X^+]$ dual-isoenzyme retention when $X$ is syntenically embedded in the $(A, B, D)$ operon.

---

### 4. Operon Topology Beyond Simple Unidirectional Directons

#### [BIO-10] `directon_id` Severs Divergons (`<-- -->`), Convergent Antisense Excludons (`--> <--`), and Riboswitch-Attenuated Operons (`>300 bp` / `TransTermHP`)
- **Severity:** **CRITICAL**
- **Exact Biological Failure Mechanism & Real Examples:**
  In `PLAN.md` §4.3 Table 1 (line 142), `directon_id` is defined as:
  *"Increments whenever strand flips or `delta_bp_prev_orf > 300 bp` or a terminator occurs."*
  And in `directon_pair_summary.parquet` (line 159), relative orientation is summarized solely by `same_directon_frac: float32`.
  - **This definition destroys three fundamental classes of prokaryotic transcriptional architecture:**
    1. **Divergons (`<-- -->`, Bidirectional Shared Promoter-Operator Pairs):**
       Across bacteria and archaea, **25–35% of cis-acting transcription factors and two-component or defense regulators are transcribed divergently (`<-- -->`, separated by $30\text{–}250\text{ bp}$) from the operon they regulate** (`Korbel & Bork 2004`; `Price & Arkin`):
       - **Regulator-Target Divergons:** `TetR`/`AcrR` $\leftarrow \ \rightarrow$ `tetA`/`acrAB` efflux pumps; `LysR` (`ilvY <-> ilvC`, `oxyR <-> katG`); `MarR`, `IclR`, `AraC`, `MerR`, `ArsR` metalloregulators.
       - **Restriction-Modification (`R-M`) & Defense Systems:** Type II R-M systems routinely arrange the Methyltransferase (`<--`) and Restriction Endonuclease (`-->`) as a divergent `<-- -->` pair controlled by a central C-protein operator so methylation is established before restriction upon horizontal transfer!
       - **Metabolic Divergons:** `SmpB <-> RnfH`, `MalT <-> malPQ`, and bidirectional siderophore/antibiotic BGC gene clusters.
       - Under `PLAN.md` v1.0, `directon_id` breaks at the `<-- -->` junction, and `same_directon_frac = 0.0`, making a tightly coupled `<-- -->` divergon indistinguishable from two unrelated genes!
    2. **Convergons, Excludons, & Antisense-Overlapping TA Modules (`--> <--`):**
       - **Excludons (`Sesto, Wurtzel, Sorek, Cossart 2013 Nat Rev Microbiol`; `Toledo-Arana 2009 Nature`):** Adjacent genes/operons with mutually exclusive functions (e.g., flagellar biosynthesis vs. `MogR` motility repressor in *Listeria*; alternative carbon or virulence pathways) are arranged in **convergent (`--> <--`) or divergent antisense-overlapping orientation** (`delta_bp < 0`, where the 3'-UTR or 5'-UTR of one operon extends across the opposing ORF), providing mutual transcriptional exclusion via RNA polymerase collision and `RNase III` dsRNA cleavage.
       - **Type I Toxin-Antitoxin Systems (`hok/sok`, `tisB/istR`, `symE/symR`, `txpA/ratA`):** A short hydrophobic toxin `smORF` (`-->`) overlaps or converges with its antisense non-coding antitoxin RNA (`<--`).
    3. **Riboswitch Attenuator Terminators & Long Leader Regions (`> 300 bp` + `TransTermHP`):**
       - Over **70% of riboswitches and T-boxes in *Bacillota* (`Firmicutes`)** regulate downstream gene expression via **conditional Rho-independent transcription terminator hairpins (attenuators)** located inside the 5' leader region between the riboswitch aptamer and the start codon of the first ORF (or inside intergenic regions of complex operons like the `cobalamin` or `S-box`/`SAM` operons)!
       - Furthermore, large riboswitches (`Cobalamin` ~220 bp, `T-box` ~280 bp, tandem `Glycine` or `SAM+AdoCbl` riboswitches `320–500 bp`, `OLE` RNA ~600 bp) frequently create an **ORF-to-ORF distance `delta_bp_prev_orf > 300 bp`** even though the RNA element itself bridges the gap with only `10–25 bp` on either side!
       - **How `PLAN.md` v1.0 breaks:** Because `TransTermHP` calls the riboswitch's intrinsic attenuator hairpin as a terminator (`feature_kind = 7`) and/or `delta_bp_prev_orf > 300 bp`, **`directon_id` increments right after the riboswitch, severing almost every transcription-attenuated riboswitch in *Bacillota* from the operon it regulates!**
- **Concrete Schema / Algorithmic Fix for `PLAN.md` v2.0:**
  1. **Fix `directon_id` Continuity Rules for Cis-Regulatory RNAs & Attenuators (`§4.3 Table 1`):**
     - Rename `delta_bp_prev_orf` to `delta_bp_prev_slot: int16` (measuring gap from the **preceding feature slot**, whether CDS or `Infernal` ncRNA, while also storing `delta_bp_prev_cds: int16`). A `350 bp` riboswitch sitting `20 bp` upstream of a CDS now has `delta_bp_prev_slot = 20 bp` (not `370 bp`), keeping `directon_id` intact!
     - **Attenuator Terminator Exemption:** Distinguish **3'-operon-stop terminators** (`feature_kind = 7`, located `-10 to +120 bp` downstream of a CDS stop codon with no overlapping `Infernal` leader) from **5'-cis-attenuator terminators** (`feature_kind = 9`, overlapping or within `80 bp` downstream of an `Infernal` riboswitch/leader RNA or `uORF_cis_leader` on the same strand). **Never increment `directon_id` on a `5'-cis-attenuator` (`feature_kind = 9`) or internal leaky stoichiometry attenuator when the downstream same-strand CDS begins within $\le 150\text{ bp}$!**
  2. **Add `regulon_locus_id: uint32` and `junction_topology_prev: uint8` to `feature_tape.parquet` (`§4.3 Table 1`):**
     - Define `junction_topology_prev: uint8` for every transition between adjacent features:
       - `0 = INTRA_ORF_FUSION` (`seg_idx_in_orf > 0`, same polypeptide)
       - `1 = CO_DIRECTIONAL_OPERON` (`same strand`, $\Delta\text{bp} \le 300\text{ bp}$)
       - `2 = DIVERGON_PROMOTER` (`-1 -> +1` head-to-head flip with shared intergenic promoter $-20\text{ bp} \le \Delta\text{bp} \le 250\text{ bp}$)
       - `3 = CONVERGON_EXUDON` (`+1 -> -1` tail-to-tail flip with antisense overlap or tight convergence $-600\text{ bp} \le \Delta\text{bp} \le +60\text{ bp}$)
       - `4 = SENSE_ANTISENSE_NESTED` (antisense RNA/uORF nested within feature coordinates)
       - `5 = INTER_LOCUS_BREAK` ($\Delta\text{bp} > 300\text{ bp}$ or uncoupled strand break)
     - Define `regulon_locus_id: uint32`, which increments **only** on `5 = INTER_LOCUS_BREAK` (thus keeping an entire `<-- [Directon A]` + `--> [Directon B]` divergon and any convergent excludon inside the same `regulon_locus_id`!).
  3. **Expand `directon_pair_summary.parquet` (`§4.3 Table 2`) to Index All 4 Topologies:**
     - Replace the single `same_directon_frac` scalar with explicit topological fractions and distances:
       - `codir_frac: float32`, `divergon_frac: float32`, `convergon_excludon_frac: float32`, `antisense_nested_frac: float32`
       - `median_codir_bp: int16`, `median_divergon_bp: int16`, `median_convergon_bp: int16`.
     - Include `junction_topology_prev` in `gLM-Domain`'s positional embedding (`§4.5` line 194) and sample full `regulon_locus_id` windows (including `<-- -->` divergons!) for `gLM-Domain` training.

---

#### [BIO-11] Riboswitch Ligand-Specificity Subtypes (`ykkC` `Guanidine-I` / `ppGpp` / `PRPP`) and Pseudo-Enzyme Scaffolds Are Collapsed by Coarse Family Tokens & `arch_hash`
- **Severity:** **MODERATE**
- **Exact Biological Failure Mechanism & Real Examples:**
  1. **Riboswitch Aptamer Specificity Switches (`Breaker Lab`):**
     Look at **Benchmark 1 (`PLAN.md` §5 line 213)**: the benchmark explicitly requires recovering `ykkC (Guanidine-I)`, `ykkC-2a (ppGpp)`, and `ykkC-2b (PRPP)` as distinct systems!
     - However, in `Rfam`, `ykkC` / `PRPP` / `ppGpp` aptamers share the same global stem-loop covariance scaffold (`RF00442` / `CL00116`) and differ only at **3 to 5 specific nucleotide positions in the bulge/P3 ligand-binding pocket** (`Sherlock, Sudarsan, Breaker 2017–2018 PNAS / eLife`)! Similarly, `ydaO` (`c-di-AMP`), `GEMM` (`c-di-GMP` vs `c-AMP-GMP`), and purine riboswitches (`Guanine` vs `Adenine` vs `2'-dG` via single `C74` vs `U74` mutation!) switch ligand identity via point mutations at conserved binding-core columns while matching the same `Infernal` covariance model.
     - If Stage 1 (`§3.1`) records only the top `Infernal` model accession (`seq_family_id = RF00442`), `Guanidine-I`, `ppGpp`, and `PRPP` riboswitches receive the **exact same `token_id`**, making it impossible for `all-the-data` to pass its own Benchmark 1!
  2. **Catalytically Inactivated Pseudo-Enzymes vs. Active Enzymes (`Koonin / Aravind`):**
     - In `PLAN.md` §2.3 (line 43), `arch_hash = xxhash64(ordered_family_ids, strand_transitions, quantized_delta_bp_bins)` deduplicates contigs, and `directon_pair_summary` aggregates by `(token_a, token_b)` without distinguishing whether `catalytic_intact` is `true` or `false`. In many signaling and defense operons, a duplicated enzyme copy loses its catalytic residues (`catalytic_intact = false`) and functions as an allosteric regulatory scaffold alongside an active copy.
- **Concrete Schema / Algorithmic Fix for `PLAN.md` v2.0:**
  1. **Covariance Alignment Pocket-Column Classifier for `Infernal` Hits (`§3.1`):**
     - After `cmscan` aligns a sequence to an `Infernal` `.cm` model, inspect the aligned nucleotides at curated **ligand-specificity consensus columns** (`cm_pocket_rules.json`: e.g., `ykkC` pocket columns separating `Guanidine-I`, `ykkC-2a/ppGpp`, `ykkC-2b/PRPP`, `ykkC-2c/NDP`; purine riboswitch position `74` separating `G`/`A`/`2'-dG`; `c-di-GMP` vs `c-AMP-GMP` `GEMM-I` columns) to assign fine-grained **`orthogroup_sub_id`** tokens for every riboswitch subtype!
  2. **Include `orthogroup_sub_id` and `catalytic_intact` in `arch_hash` and `directon_pair_summary`:**
     - Compute `arch_hash = xxhash64(ordered_orthogroup_sub_ids, catalytic_intact_bits, junction_topology_prev, quantized_delta_bp_bins)` so operons with specificity switches or pseudo-enzyme neofunctionalizations are never prematurely deduplicated.

---

## Summary Checklist of Required `PLAN.md` v2.0 Revisions

1. **§2.3 Streaming Sieve:** Whitelist `Obelisks` / circular ribozyme-viroids (`>= 750 bp`) in the 2-Feature Gate; defer `arch_hash` deduplication until `orthogroup_sub_id`, `catalytic_intact`, and `CTD` slicing are computed, and never drop CRISPR spacer or polymorphic toxin cassette diversity.
2. **§3.1 Stage 1 (Gene, `smORF`, & ncRNA Calling):**
   - Replace hardcoded `AGGAG`/`GGAGG` `smORF` seed with 5-mode initiation (`pyrodigal` learned initiation table + SD + translational coupling + archaeal/actinobacterial leaderless `-10`/`TATA` + *Bacteroidota* AT-rich window + HMM/TM/dN-dS rescue).
   - Add **Stage 1.4 Overlap Arbitration** (purge spurious ORFs overlapping CRISPR arrays/tRNAs/rRNAs/riboswitches; retain and parent-link sense `uORF_cis_leader`, `CreTA`, and antisense Type I TA pairs).
   - Add **Stage 1.5 Split-Gene / Frameshift / Intron / `Sec`-`Pyl` Reconstructor** (`prfB`, `gpG-gpGT`, IS elements, selenoproteins, and Group I/II intron + HEG exon bridging).
   - Add **Riboswitch Binding-Pocket Subtyper** (`ykkC-Guanidine-I / ppGpp / PRPP`, purine `G / A / dG`, `c-di-GMP / c-AMP-GMP`).
3. **§3.2 Stage 2 (`NR90` & 4-Cascade Domain Dissector):**
   - Guard `NR90` clustering with a **120-aa N/C-terminal local identity check** so polymorphic toxins (`Rhs`, `CdiA`, `VgrG`, `LXG`, `MafB`) with different `CTD` tips are never collapsed.
   - Slice `pre-tox`-bounded C-terminal toxin tips (`CTDs`) into standalone segments for `hhsearch` and `ProstT5`/`Foldseek`.
   - Stitch **discontinuous / nested domains** (`parent_seg_idx_in_orf`) and **circular permutations (`is_circularly_permuted`)** *before* evaluating `catalytic_intact`.
   - Collapse **homotypic tandem repeat arrays** (`TPR`, `ANK`, `WD40`, `LRR`, `Pentapeptide`) into single super-segments with `repeat_unit_count`.
4. **§4.1–4.3 Schema & Indices (`feature_tape.parquet` & `directon_pair_summary.parquet`):**
   - Promote `orthogroup_sub_id: uint32` (`40% ID` + active-site pocket signature + syntenic context split) to a first-class column across all tables and bitvectors.
   - Add `regulon_locus_id: uint32` and `junction_topology_prev: uint8` (`INTRA_ORF`, `CO_DIRECTIONAL`, `DIVERGON`, `CONVERGON_EXUDON`, `ANTISENSE_NESTED`, `BREAK`) and exempt `5'-cis-attenuator` terminators from breaking `directon_id`.
5. **§4.4 Index 3 (Phylogenetic Bitvectors & NOGD):**
   - Replace single-genome species-rep bitvectors with the 3-plane pangenome tuple (`core_chr_bits`, `pan_any_bits`, `mge_accessory_bits`) + `hq_closed_or_multi_strain_mask` + `operon_not_at_contig_edge` veto + clade-local asymmetric NOGD scoring.
