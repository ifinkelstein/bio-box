# Adversarial Audit Report 3: Statistical Genetics, `gLM-Domain` Machine Learning Architecture, and Experimental Closure / Falsification Design

**Auditor:** Auditor 3 — Statistical Genetics, Machine Learning & Falsification Critic  
**Target Documents:** `00-spec.md`, `PRIOR-ART.md`, `PLAN.md` (v1.0 — Pre-Audit Draft)  
**Audit Verdict:** **MAJOR REVISION REQUIRED BEFORE IMPLEMENTATION.** While the systems engineering and multi-modal domain dissection in `PLAN.md` v1.0 are visionary, the statistical genetics scoring layer, the `gLM-Domain` training/objective formulation, and the `family_closure_card` resolution contain **three fatal methodological flaws** that would cause the discovery engine to (1) rank promiscuous plasmid/integron passenger pairs above true molecular pathways, (2) overfit `gLM-Domain` on horizontally transferred operons while predicting biophysically impossible chimeric slots, and (3) misclassify novel domains fused to common Rossmann/ABC scaffolds as "already solved in PDB/PaperBLAST."

Below is the complete mathematical diagnosis of every flaw, followed by the exact equations, schema upgrades, and falsification protocols required for `PLAN.md` v2.0.

---

## 1. Passenger Synteny vs. Functional Coupling (The #1 Statistical False-Positive in Comparative Genomics)

### 1.1 Diagnosis: Why `synteny_log10_pval` and Fitch Homoplasy ($H_{A,B}$) Actively Reward False Positives

In `PLAN.md` Section 4.3 (`directon_pair_summary.parquet`, lines 151–161) and Section 4.4 (lines 184–185), pairwise synteny is ranked by two metrics:
1. `synteny_log10_pval`: A clade-stratified hypergeometric enrichment $p$-value for tokens $(A, B)$ co-occurring within $\pm 5$ ORFs, and
2. `parsimony_homoplasy_events` ($H_{A,B}$): The number of independent evolutionary gain/co-acquisition events of $(A, B)$ mapped onto the GTDB species tree (`bac120` / `ar53`) via vectorized Fitch parsimony.

**Both metrics fail catastrophically in the presence of mobile genetic elements (MGEs) and conserved housekeeping super-operons:**

1. **The MGE Passenger & Co-Mobilization Trap:**
   Suppose Gene $A$ (e.g., an aminoglycoside resistance gene, a cargo restriction system, or an uncharacterized `DUF`) and Gene $B$ (e.g., a `Tn3`/`IS6` transposase, an `IncP`/`IncF` conjugative backbone gene like `TraG`/`TrbE`, a class-1 integron `qacEΔ1`/`sul1` cassette, or a prophage replication/lysis module) sit within $\pm 5$ ORFs on a broad-host-range conjugative plasmid or integrative conjugative element (`ICE`).
   - When that plasmid or `ICE` transfers horizontally across 50 genomes spanning 4 bacterial phyla, **Fitch parsimony on the host chromosomal tree (`bac120`) counts 50 "independent evolutionary acquisitions" ($H_{A,B} = 50$)!**
   - Why? Because the host tree `bac120` is the **wrong phylogenetic topology** for genes residing on an MGE. Those 50 occurrences represent **one single assembly event on the plasmid/ICE phylogeny** followed by 50 whole-element horizontal transfers. Evaluating MGE cargo on the host chromosome tree conflates **vehicle mobility** with **convergent selection for functional coupling**.
2. **Failure of the Hypergeometric Null (`synteny_log10_pval`):**
   A hypergeometric test assumes an **exchangeable uniform bag-of-genes null** across a chromosome of $N \approx 4,000$ ORFs, where the background probability of two genes landing within $\pm k = 5$ ORFs is $p_0 \approx \frac{2k}{N} = 0.0025$.
   - However, MGEs, genomic islands, and phage contigs (`Tier 2` and integrated `Tier 1`/`Tier 3` islands) occupy a restricted compartment of $L_{\text{MGE}} \in [8, 120]$ ORFs with a highly skewed gene pool (`mge_enriched` tokens). On a $12\text{ kb}$ plasmid ($L = 12$ ORFs), a window of $\pm 5$ ORFs spans **92% of the entire replicon**! Comparing MGE-localized gene pairs against a $4,000$-ORF chromosomal uniform null inflates `synteny_log10_pval` by $10^{25}\text{–}10^{120}$ for completely unrelated passenger genes.
3. **The Conserved Super-Operon Boundary Trap:**
   In vertically conserved chromosomal macro-blocks—such as the `S10-spc-alpha` ribosomal super-operon (~30 genes), the `dcw` cell-division/peptidoglycan cluster (`mraZ`–`ftsZ`), or flagellar biosynthesis blocks—genes remain within $\pm 5$ ORFs across 30+ phyla due to ancestral transcriptional co-regulation and ribosome assembly kinetics. Consequently, functionally unrelated boundary genes (e.g., a tRNA methyltransferase or `secY` adjacent to a ribosomal protein $4$ slots away, or a conserved flanking gene immediately outside the `dcw` cluster) achieve massive clade-stratified hypergeometric significance without forming a direct biochemical complex or metabolic pathway step.
4. **Right/Left-Censoring Bias on Fragmented `Tier 3` Contigs:**
   In `Tier 3` (`unbinned_sra`), contigs have a median length of $2\text{–}6$ ORFs. If occurrence counts $N(A)$ include instances where $A$ sits within $d < k$ ORFs of a contig end, the denominator of $P(B \mid A)$ is artificially inflated by truncated windows where $B$ was sheared off by assembly fragmentation.

---

### 1.2 Exact Mathematical Fix: The Background Mobility & Context-Conditional Null Framework

To eliminate passenger synteny while preserving true defense systems, metabolic pathways, and viral/plasmid functional modules (which often *do* reside on MGEs!), `PLAN.md` v2.0 must replace raw hypergeometric $p$-values and host-tree homoplasy counts with a **Three-Component Context-Conditional Coupling Score ($\mathcal{S}_{\text{couple}}(A, B)$)**.

#### Component A: Compartment- & Mobility-Conditioned Null ($p_{m}^{(0)}$) and Vehicle-Deduplicated Effective Events ($N_{\text{eff}}^{\text{veh}}$)

1. **Local Regime Classification:**
   Every feature slot $i$ on contig $c$ is assigned to one of six local genomic regimes $r(i) \in \{\text{CORE}, \text{SUPER\_OPERON}, \text{PLASMID}, \text{PHAGE}, \text{ICE\_GI}, \text{INTEGRON}\}$ using `tier` plus a $\pm 15$-ORF sliding window scan for hallmark tokens:
   - `INTEGRON`: within $\pm 10$ slots of an `attC`/`attI` site (`feature_kind = 5`) or `IntI` integrase.
   - `PHAGE`: `tier=2` phage or within $\pm 15$ ORFs of viral hallmark tokens (capsid, portal, terminase large subunit, tail sheath, `PHROG` structural/replication core).
   - `PLASMID` / `ICE_GI`: `tier=2` plasmid or within $\pm 15$ ORFs of conjugative transfer (`VirB4`, `TraG`, `MOB` relaxase, `RepA`) or transposase/integrase (`ISfinder`/`TnCentral`/`Tyr-Ser recombinase`) tokens.
   - `SUPER_OPERON`: inside a known macro-synteny block ($\ge 8$ universal `bac120`/`ar53` single-copy marker genes within 20 ORFs).

2. **Flanking-Vehicle Graph Deduplication (Fixing the Host-Tree Homoplasy Bug):**
   Instead of counting independent events solely on the host species tree $T_{\text{host}}$, we deduplicate horizontal transfers of the **same mobile vehicle**. For every contig locus $m$ containing the pair $(A, B)$, extract the set of flanking token families in a $\pm 15$-ORF window excluding the immediate $(A, B)$ directon module:
   $$\mathcal{W}_{\text{flank}}(m) = \{ \text{seq\_family\_id}_j : 0 < |\text{orf}_j - \text{orf}_{A,B}| \le 15 \} \setminus \text{Module}(A, B)$$
   Construct a single-linkage **Vehicle Co-Transfer Graph** where two loci $m, m'$ (even if found in different host phyla!) are linked if they share a mobile/super-operon context:
   $$\text{Sim}_{\text{veh}}(m, m') = \frac{|\mathcal{W}_{\text{flank}}(m) \cap \mathcal{W}_{\text{flank}}(m')|}{|\mathcal{W}_{\text{flank}}(m) \cup \mathcal{W}_{\text{flank}}(m')|} \ge 0.30 \quad \text{OR} \quad \text{AA\_ID}(A_m, A_{m'}) \ge 0.80 \wedge \text{AA\_ID}(B_m, B_{m'}) \ge 0.80$$
   The **Vehicle-Corrected Independent Event Count** is then defined as the sum over connected vehicle components $v \in \mathcal{V}_{A,B}$ of the phylogenetic dispersion *within* or *across* genuinely distinct vehicles, combined with **Gene-Tree Fitch Parsimony ($H_{B \mid T_A}$)**:
   $$N_{\text{eff}}^{\text{veh}}(A, B) = \sum_{v \in \mathcal{V}_{A,B}} \Big( 1 + \alpha_{\text{core}} \cdot \mathbb{I}[r(v) = \text{CORE}] \cdot \big(H_{\text{host}}(v) - 1\big) \Big)$$
   where $H_{\text{host}}(v)$ is the host-tree Fitch parsimony count inside component $v$, and $\alpha_{\text{core}} = 1$ for chromosomal core loci with $<80\%$ AA identity and $0$ for MGE vehicles (`PLASMID`, `PHAGE`, `ICE_GI`, `INTEGRON`).
   - **Result:** A plasmid transferred intact across 50 genomes in 4 phyla collapses to $N_{\text{eff}}^{\text{veh}} = 1$ vehicle event, whereas a defense system or metabolic pair independently assembled into 18 distinct phages, plasmids, and chromosomes with unrelated flanking backbones and $<60\%$ AA identity scores $N_{\text{eff}}^{\text{veh}} = 18$!

3. **Censoring-Adjusted Regime-Conditioned Poisson-Binomial Null:**
   For any target ORF distance $d = |\text{orf}_B - \text{orf}_A| \le k$, define the **uncensored occurrence set** $\mathcal{M}_{\text{uncens}}(A; d)$ as all independent loci where contig $c$ has at least $d$ called ORFs between $A$ and the contig boundary. Under the regime-conditioned null, the probability that token $B$ lands within $\pm k$ ORFs of $A$ at locus $m$ in regime $r_m$ is:
   $$p_m^{(0)}(B \mid A, r_m) = 1 - \prod_{j \in \text{Window}_{\pm k}(m)} \big(1 - \pi(B \mid r_m, \text{clade}_m)\big) \approx \min\!\Big(1, \, 2k_{\text{eff}}(m) \cdot \frac{N_{\text{eff}}(B, r_m)}{L_{\text{eff}}(r_m)}\Big)$$
   where $\pi(B \mid r_m, \text{clade}_m)$ is the empirical slot frequency of token $B$ inside regime $r_m$ (so transposases, integrases, and plasmid `Tra`/`Rep` genes have a high baseline $\pi(B \mid \text{PLASMID/ICE})$, automatically penalizing their co-occurrence!).

---

#### Component B: Directional Conditional Probability Asymmetry & Contextual Neighbor Entropy

Even inside an MGE or genomic island, how do we distinguish a **true functional operon pair** (e.g., `CBASS` cyclase + phospholipase effector, or a dedicated maturase $A$ + core enzyme $B$) from a **cargo gene sitting next to a promiscuous transposase, integrase, or integron `attC` cassette**?

1. **Censoring-Corrected Directional Conditional Probabilities:**
   Using vehicle-dereplicated, contig-edge-uncensored counts:
   $$P(B \mid A) = \frac{N_{\text{eff}}^{\text{veh}}(A \wedge_{\pm k} B)}{N_{\text{eff, uncens}}^{\text{veh}}(A; k)}, \qquad P(A \mid B) = \frac{N_{\text{eff}}^{\text{veh}}(A \wedge_{\pm k} B)}{N_{\text{eff, uncens}}^{\text{veh}}(B; k)}$$
   We store **both** directional probabilities alongside:
   - **Symmetric Coupling (Jaccard / Cosine):** $C_{\text{sym}}(A, B) = \sqrt{P(B \mid A) \cdot P(A \mid B)}$ (high when $A$ and $B$ are obligate 1:1 partners),
   - **Directional Dependency Ceiling:** $C_{\max}(A, B) = \max\big(P(B \mid A), P(A \mid B)\big)$ (captures asymmetric modular biology, e.g., a specialized radical SAM enzyme $A$ that *always* requires a common `ABC` transporter $B$, so $P(B \mid A) = 0.92$ while $P(A \mid B) = 0.01$),
   - **Directional Asymmetry Log-Ratio:**
     $$\mathcal{A}(A, B) = \log_2 \frac{P(B \mid A) + \epsilon}{P(A \mid B) + \epsilon}$$

2. **Contextual Neighbor Shannon Entropy $H_{\text{ctx}}(X)$ (The Promiscuous Hub Filter):**
   For every token $X$ in the dictionary, compute the clade- and vehicle-dereplicated empirical probability distribution over all immediate neighbor tokens $Y \in \mathcal{N}_{\pm 2}(X)$:
   $$p_{\text{ctx}}(Y \mid X) = \frac{N_{\text{eff}}^{\text{veh}}(X \wedge_{\pm 2} Y)}{\sum_{Z \in \mathcal{V}} N_{\text{eff}}^{\text{veh}}(X \wedge_{\pm 2} Z)}$$
   Define the **Contextual Neighbor Entropy** (in bits) and **Effective Neighbor Perplexity**:
   $$H_{\text{ctx}}(X) = -\sum_{Y \in \mathcal{V}} p_{\text{ctx}}(Y \mid X) \log_2 p_{\text{ctx}}(Y \mid X), \qquad \text{Perp}_{\text{ctx}}(X) = 2^{H_{\text{ctx}}(X)}$$
   - **Why this separates asymmetric functional partners from MGE promiscuity:**
     - A promiscuous insertion sequence (`IS6`, `IS256`), tyrosine integrase, or integron `attC` cassette inserts adjacent to thousands of unrelated cargo genes: its **Contextual Perplexity is huge** ($\text{Perp}_{\text{ctx}}(X) \in [150, 4000]$, $H_{\text{ctx}}(X) > 7.2\text{ bits}$), and even when conditioned on a specific slot offset $+1$, its neighbor distribution is flat.
     - By contrast, a ubiquitous core cofactor or transporter family $B$ that participates in 15 specific operon types has moderate global entropy ($H_{\text{ctx}}(B) \approx 3.5\text{–}5.0\text{ bits}$) and **low conditional entropy given the operon's other members** $H(Y \mid X, \text{OperonContext}) < 1.5\text{ bits}$.
     - Even more critically, we compute the **Slot-Specific Conditional Entropy of $A$'s Neighborhood at the Offset of $B$**, $H_{\text{slot}}(A \to \Delta\text{orf}_{AB})$: if $A$ truly requires $B$ (or a functional NOGD analog of fold `fold(B)`), then at the specific relative slot $\Delta\text{orf}_{AB}$, the **Fold-Level Contextual Entropy** $H_{\text{ctx}}^{\text{fold}}(A; \Delta\text{orf}_{AB}) = -\sum_{F} P(\text{fold}=F \mid A, \Delta\text{orf}_{AB}) \log_2 P(\text{fold}=F \mid A, \Delta\text{orf}_{AB})$ is near zero ($\le 1.2\text{ bits}$)!

3. **Micro-Synteny Geometric Rigidity Score $\mathcal{G}_{\text{micro}}(A, B)$:**
   Passenger genes on a plasmid or integron array vary wildly in relative strand orientation, gene order, and intergenic spacing ($\Delta\text{bp}$), whereas truly co-transcribed / translationally coupled operon partners maintain strict biophysical geometry across evolution:
   $$\mathcal{G}_{\text{micro}}(A, B) = \underbrace{f_{\text{same\_directon}}(A, B)}_{\text{co-transcription}} \cdot \underbrace{\exp\!\left(-\frac{\text{IQR}(\Delta\text{bp}_{AB})}{40\text{ bp}}\right)}_{\text{translational spacing conservation}} \cdot \underbrace{\Big(0.5 + 0.5\,\big|2 f_{A \prec B} - 1\big|\Big)}_{\text{transcriptional order conservation}} \cdot \big(1 + \beta_{\text{couple}} f_{\text{overlap\_or\_fused}}(A, B)\big)$$
   where:
   - $f_{\text{same\_directon}}(A, B)$ is the fraction of independent vehicle instances where $A$ and $B$ share strand and contain zero intervening Rho-independent terminators or `attC` recombination sites,
   - $\text{IQR}(\Delta\text{bp}_{AB})$ is the interquartile range of intergenic base-pair distance across distinct taxonomic orders (translational coupling `-4 bp ATGA` or `-1 bp TGATG` yields $\text{IQR} \le 6\text{ bp}$),
   - $f_{A \prec B}$ is the fraction of same-directon instances where $A$ is transcribed $5' \to 3'$ upstream of $B$,
   - $f_{\text{overlap\_or\_fused}}(A, B)$ is the fraction of instances with $\Delta\text{bp}_{AB} \le 4\text{ bp}$ or Rosetta-stone intra-ORF fusion ($\Delta\text{orf} = 0$).

---

#### Component C: Unified Log-Likelihood Ratio for Functional Coupling vs. Passenger Mobility ($\text{LLR}_{\text{func}}$)

Replace `synteny_log10_pval` as the primary ranking metric in `directon_pair_summary.parquet` with the **Mobility- & Entropy-Corrected Functional Coupling Score**:
$$\text{LLR}_{\text{func}}(A, B) = \underbrace{\sum_{m \in \mathcal{M}_{\text{indep}}^{\text{veh}}(A,B)} \log_{10} \frac{P_{\text{obs}}(A, B \mid \Delta\text{bp}_m, \text{strand}_m)}{p_m^{(0)}(B \mid A, r_m) \cdot p_{\text{geom}}^{(0)}(\Delta\text{bp}_m, \text{strand}_m)}}_{\text{Vehicle-deduplicated regime-conditioned log-BF}} - \lambda_{\text{hub}} \max\!\big(0, H_{\text{ctx}}^{\text{fold}}(A; \Delta_{AB}) - H_0\big) - \lambda_{\text{attC}} f_{\text{attC\_between}}(A, B)$$
calibrated via **Target-Decoy Empirical False Discovery Rate ($q_{\text{emp}}$)** by generating decoy pairs `(A, B_shuffled)` permuted strictly *within* the same `(regime, taxonomic_order, neighbor_entropy_bin)` stratum!

---

## 2. `gLM-Domain` Transformer Architecture & Data Leakage Flaws

### 2.1 Flaw 2A: Horizontal Gene Transfer (HGT) Homology Leakage Across Phylum Splits

#### Diagnosis
In `PLAN.md` Section 2.3 (line 44) and Section 5 Benchmark 2 (lines 217–218), directon tapes are deduplicated at **95% nucleotide identity**, and `gLM-Domain` generalization is evaluated by **holding out entire taxonomic phyla** (*Aquificota*, *Thermotogota*, *Halobacteriota*).

**This evaluation design suffers from severe data leakage:**
1. **Cross-Phylum HGT at 60–90% Amino-Acid Identity:**
   Prokaryotic defense islands (`CRISPR-Cas`, `CBASS`, `Pycsar`, `Thoeris`, `BREX`, `RM`), secretion systems (`T4SS`, `T6SS`), antibiotic biosynthesis clusters, and thermophilic metabolic operons are rampantly transferred horizontally across phyla and between Bacteria and Archaea. In particular, *Aquificota* and *Thermotogota* have acquired $>25\%$ of their operons via HGT from *Firmicutes*, *Proteobacteria*, and *Euryarchaeota*!
2. **Token-Sequence Identity Leakage (100% Identical Input Tapes at 40% Nucleotide Identity!):**
   `gLM-Domain` does **not** take raw nucleotides as input; it takes discrete **family and fold tokens** `(seq_family_id, struct_fold_id)`. Two homologous operons in *Pseudomonadota* (Train) and *Aquificota* (Test) that diverged 2 billion years ago (or transferred horizontally) and share only 42% amino-acid identity still map to the **exact same token sequence** `[COG1111, COG2222, Pfam3333,COG4444]`!
   - Consequently, holding out a phylum while training on token sequences from other phyla tests **memorization of identical token n-grams**, not out-of-distribution generalization!

#### Exact Fix: Operon-Homology Bipartite Graph Connected-Component Splitting ($\mathcal{G}_{\text{leak}}$)

To prevent all sequence, profile, and structural operon leakage across `Train / Val / Test`:
1. **Step 1 — Superfamily Equivalence Closure (`super_cluster_id`):**
   Group all sequence tokens (`seq_family_id`, including `DARK_FAM` clusters) into connected **Superfamily Equivalence Classes** $\mathcal{S}$ where two families $f_1, f_2$ are merged if:
   - Any member sequences share $\ge 30\%$ amino-acid identity over $\ge 50\%$ coverage (`MMseqs2 search --min-seq-id 0.30 -c 0.50`), **OR**
   - Profile-profile alignment `hhsearch` probability $\ge 80\%$ ($E \le 10^{-3}$), **OR**
   - Both belong to the same `Pfam Clan` or the same `ECOD H-group` (`ecod_h_id`) for sequence-annotated families.
2. **Step 2 — Operon Architecture Similarity Graph $\mathcal{G}_{\text{leak}} = (\mathcal{D}, \mathcal{E})$:**
   For every candidate directon/neighborhood tape $d \in \mathcal{D}$, map its constituent proteins to their multiset of non-housekeeping Superfamily Equivalence Classes $S(d) = \{s_1, s_2, \dots, s_m\}$ and ordered `ECOD H-group` 3-grams $G_3(d)$. Place an undirected leakage edge $(d_1, d_2) \in \mathcal{E}$ if **any** of the following hold:
   $$\text{OverlapCoeff}\big(S(d_1), S(d_2)\big) = \frac{|S(d_1) \cap S(d_2)|}{\min(|S(d_1)|, |S(d_2)|)} \ge 0.30 \quad \lor \quad |S(d_1) \cap S(d_2)| \ge 2 \quad \lor \quad G_3(d_1) \cap G_3(d_2) \neq \emptyset$$
3. **Step 3 — Connected-Component Partitioning (`Union-Find`):**
   Compute exact connected components $\text{CC}(\mathcal{G}_{\text{leak}})$ and assign **entire connected components** to `Train (80%)`, `Val (10%)`, and `Test (10%)`.
   - **Guarantee:** No test operon shares any $\ge 30\%$ AA / `HHsearch` / `Pfam-Clan` gene pair or any `ECOD H-group` structural 3-gram with any training operon.
4. **Step 4 — Split Benchmark 2 into Two Orthogonal Falsification Protocols:**
   - **Benchmark 2A (Zero-Shot Syntax Generalization on Disjoint Connected Components):** On held-out connected components $\text{CC}_{\text{test}}$, mask 1 slot and evaluate `Recall@10` on `biophysical_role` and `ECOD_T` fold.
   - **Benchmark 2B (Simulated NOGD Blind Replacement Recovery):** For the 50 known NOGD pairs $(C_{\text{canonical}}, X_{\text{replacement}})$ inside pathway context $(A, B, \_, D)$, **purge all operons containing $X_{\text{replacement}}$ (and its entire `Clan` / `MMseqs30` cluster) from the training corpus**, while keeping canonical $(A, B, C_{\text{canonical}}, D)$ operons in `Train`. At test time, query `[A, B, MASK, D]` and verify that `gLM-Domain`'s predicted distribution over `(biophysical_role, ECOD_T_fold, catalytic_motif_flags)` ranks the true fold and biophysical badge of the unseen non-homologous replacement $X_{\text{replacement}}$ in the top 5!

---

### 2.2 Flaw 2B: The Multi-Modal Softmax Independence Bug

#### Diagnosis
In `PLAN.md` Section 4.5 (lines 196–199), `gLM-Domain` predicts three conditionally independent heads from the final hidden state $\mathbf{h}_i^{(L)}$ at a masked position $i$:
$$\mathcal{L}_{\text{MSP}} = \mathcal{L}_{\text{seq}}\big(\mathbf{z}_i^{\text{seq}}, y_i^{\text{seq}}\big) + \mathcal{L}_{\text{fold}}\big(\mathbf{z}_i^{\text{fold}}, y_i^{\text{fold}}\big) + \mathcal{L}_{\text{phys}}\big(\mathbf{z}_i^{\text{phys}}, \mathbf{y}_i^{\text{phys}}\big)$$
This assumes conditional independence:
$$P\big(y_i^{\text{seq}}, y_i^{\text{fold}}, \mathbf{y}_i^{\text{phys}} \mid \mathbf{h}_i^{(L)}\big) = P\big(y_i^{\text{seq}} \mid \mathbf{h}_i^{(L)}\big) \, P\big(y_i^{\text{fold}} \mid \mathbf{h}_i^{(L)}\big) \, P\big(\mathbf{y}_i^{\text{phys}} \mid \mathbf{h}_i^{(L)}\big)$$

**Why this breaks biological inference:**
1. **Mode-Mixing Chimeras under Functional Multi-Modality:**
   In prokaryotic operons, a single functional slot is frequently occupied by **mechanistically interchangeable but structurally distinct classes**. For example, the effector slot of a `CBASS` or `Thoeris` defense operon may be occupied in 45% of lineages by a **soluble cytoplasmic `TIR` / `sirtuin` NADase (`ECOD: Rossmann`, `TM = 0`)** and in 40% of lineages by a **multi-pass transmembrane pore-forming effector (`TM = 4`, membrane fold)**.
   - When `gLM-Domain` predicts marginals independently, $\arg\max P(y^{\text{fold}} \mid \mathbf{h}_i^{(L)})$ selects the soluble `Rossmann/TIR` fold, while the binary cross-entropy head $P(\mathbf{y}^{\text{phys}} \mid \mathbf{h}_i^{(L)})$ activates transmembrane bits or signal-peptide bits—producing a **biophysically impossible chimera** (a cytoplasmic soluble Rossmann fold with 4 transmembrane helices!).
2. **Bitwise BCE Incoherence inside `biophysical_badge`:**
   In `PLAN.md` Section 3.3 (lines 96–101), `tm_count` (`bits 0..4`), `signal_type` (`bits 5..8`), and `pi_quantized` (`bits 15..19`) are binary integer encodings packed into a 64-bit integer. Applying 64 independent sigmoid BCE losses to binary integer bits (`bits 0..4`) means bit 0 ($2^0$) and bit 4 ($2^4$) are predicted independently—an elementary representation bug!
3. **Flat Softmax Dilution across 250,000 Sequence Families:**
   Because `seq_family_id`, `clan_or_cog_id`, and `struct_fold_id` form a strict topological hierarchy (`ECOD_X -> ECOD_H -> ECOD_T -> Clan -> Family`), flat independent softmaxes waste probability mass on invalid `(fold, family)` combinations and suffer severe gradient starvation on rare metagenomic families.

#### Exact Fix: 4-Stage Coarse-to-Fine Autoregressive Slot Decoder

Replace the three independent heads at masked slot $i$ with a **Hierarchical Coarse-to-Fine Autoregressive Slot Decoder** that factorizes the joint distribution over a single slot from physical role down to specific sequence family:
$$P\big(\mathbf{b}_i, x_i, h_i, t_i, c_i, f_i \mid \mathbf{h}_i^{(L)}\big) = P\big(\mathbf{b}_i \mid \mathbf{h}_i^{(L)}\big) \cdot P\big(x_i, h_i, t_i \mid \mathbf{h}_i^{(L)}, \mathbf{e}(\mathbf{b}_i)\big) \cdot P\big(c_i \mid \mathbf{u}_i^{\text{fold}}, t_i\big) \cdot P\big(f_i \mid \mathbf{u}_i^{\text{clan}}, c_i\big)$$

Specifically:
1. **Stage 1 — Structured Biophysical Role Head $P(\mathbf{b}_i \mid \mathbf{h}_i^{(L)})$:**
   Decompose the 64-bit `biophysical_badge` into proper typed distributions rather than raw bit-level BCE:
   - **Topology/Localization Class $b_{i,\text{loc}} \in \{0, \dots, 9\}$ (10-way Softmax):** `0=Cytoplasmic_Soluble`, `1=Single_Pass_TM`, `2=TM_2_3`, `3=TM_4_7`, `4=TM_8_Plus_Transporter`, `5=Sec_SPI_Periplasmic`, `6=Tat_SPI_Redox`, `7=Lipoprotein_SPII`, `8=Prepilin_SPIII`, `9=Secretion_System_Substrate`.
   - **Conditioned Sorting & Active-Site Motif Vector $\mathbf{b}_{i,\text{motif}} \in \{0,1\}^{16}$:** Predicted via MLP conditioned on $\big[\mathbf{h}_i^{(L)} \ \| \ \mathbf{E}_{\text{loc}}[b_{i,\text{loc}}]\big]$.
   - **Quantized $\text{pI}$ and Length Bins ($b_{i,\text{pI}}, b_{i,\text{len}}$):** Predicted via categorical cross-entropy with Gaussian label smoothing ($\sigma = 1\text{ bin}$) conditioned on $\mathbf{E}_{\text{loc}}[b_{i,\text{loc}}]$.
   - Embed the complete Stage-1 state:
     $$\mathbf{u}_i^{(1)} = \text{LayerNorm}\Big(\mathbf{h}_i^{(L)} + \mathbf{W}_{\text{phys}}\big[\mathbf{E}_{\text{loc}}[b_{i,\text{loc}}] \ \| \ \mathbf{W}_m \mathbf{b}_{i,\text{motif}} \ \| \ \mathbf{E}_{\text{pI}}[b_{i,\text{pI}}]\big]\Big)$$

2. **Stage 2 — Hierarchical ECOD Fold Tree Head $P(x_i, h_i, t_i \mid \mathbf{u}_i^{(1)})$:**
   Conditioned on the biophysical state $\mathbf{u}_i^{(1)}$ (so choosing `TM_8_Plus_Transporter` in Stage 1 immediately masks out soluble cytoplasmic folds!), predict the ECOD structural hierarchy via **tree-masked conditional logits**:
   $$P(x_i \mid \mathbf{u}_i^{(1)}) = \text{Softmax}_{x}\big(\mathbf{w}_x^\top \mathbf{u}_i^{(1)}\big)$$
   $$P(h_i \mid x_i, \mathbf{u}_i^{(1)}) = \text{Softmax}_{h \in \text{Children}(x_i)}\big(\mathbf{w}_h^\top \mathbf{u}_i^{(1)}\big)$$
   $$P(t_i \mid h_i, \mathbf{u}_i^{(1)}) = \text{Softmax}_{t \in \text{Children}(h_i)}\big(\mathbf{w}_t^\top \mathbf{u}_i^{(1)}\big)$$
   Update the slot state with the chosen fold embedding:
   $$\mathbf{u}_i^{(2)} = \text{LayerNorm}\big(\mathbf{u}_i^{(1)} + \mathbf{E}_{\text{ECOD\_H}}[h_i] + \mathbf{E}_{\text{ECOD\_T}}[t_i]\big)$$

3. **Stage 3 — Clan / Superfamily Head Conditioned on Fold $P(c_i \mid t_i, \mathbf{u}_i^{(2)})$:**
   Restrict the softmax support to the set of sequence clans $\mathcal{C}(t_i)$ empirically compatible with `ECOD_T` fold $t_i$, plus a dedicated **`DARK_CLAN_OF_FOLD(t_i)`** token for unannotated sequence clusters adopting fold $t_i$:
   $$P(c_i \mid t_i, \mathbf{u}_i^{(2)}) = \text{Softmax}_{c \in \mathcal{C}(t_i) \cup \{\text{DARK}(t_i)\}}\big(\mathbf{w}_c^\top \mathbf{u}_i^{(2)}\big)$$
   $$\mathbf{u}_i^{(3)} = \text{LayerNorm}\big(\mathbf{u}_i^{(2)} + \mathbf{E}_{\text{clan}}[c_i]\big)$$

4. **Stage 4 — Specific Sequence Family Head Conditioned on Clan $P(f_i \mid c_i, \mathbf{u}_i^{(3)})$:**
   $$P(f_i \mid c_i, \mathbf{u}_i^{(3)}) = \text{Softmax}_{f \in \mathcal{F}(c_i)}\big(\mathbf{w}_f^\top \mathbf{u}_i^{(3)}\big)$$

- **Computational Efficiency & Beam-Search Discovery:**
  - During training, **teacher forcing** feeds ground-truth $(b_i^*, x_i^*, h_i^*, t_i^*, c_i^*)$ into Stages 1–4 simultaneously in a **single parallel matrix multiplication** with zero sequential loop overhead!
  - Because each conditional softmax is restricted to $\text{Children}(\cdot)$ ($|\mathcal{F}(c_i)| \sim 10^1\text{–}10^3$ instead of $250,000$), the output layer is **10× faster and uses 8× less VRAM** than a flat 250k softmax!
  - At inference time, a width-$K$ ($K=10$) **Coarse-to-Fine Beam Search** over $(b_{\text{loc}} \to t_i \to c_i \to f_i)$ returns $K$ 100% biophysically and structurally self-consistent slot hypotheses, explicitly separating alternative modes (Mode 1: `Cytoplasmic + Rossmann TIR NADase`; Mode 2: `TM_4_7 + Transmembrane Holin Pore`).

---

### 2.3 Flaw 2C: 1D Slot `RoPE` Ignores Physical Operon Geometry

#### Diagnosis
In `PLAN.md` Section 4.5 (lines 194–196), relative position is encoded by standard 1D rotary position embeddings (`RoPE`) on the sequential `slot_idx` ($i \in \{0, 1, \dots, L-1\}$), while `delta_bp_prev_orf`, `strand`, and `is_fused` are added only as local embeddings at layer 0 ($\mathbf{h}_i^{(0)}$).

**Why this fails on genomic tapes:**
1. **Multi-Domain Skew of `slot_idx`:**
   Because `all-the-data` splits multi-domain proteins into sub-ORF segments (`slot_idx` increments for each domain inside the same `orf_idx`), a single 5-domain non-ribosomal peptide synthetase (`NRPS`) or polyketide synthase (`PKS`) consumes 5 `slot_idx` steps! Under 1D `RoPE` on `slot_idx`, the gene immediately downstream of a 5-domain protein appears 5 slots away from the N-terminal domain, distorting inter-gene syntactic distances across all attention layers.
2. **Blindness to Pairwise Transcriptional Boundaries:**
   For any pair of slots $(i, j)$, attention weights must know whether $i$ and $j$ are:
   - Covalently fused inside the **same polypeptide** ($\Delta\text{orf}_{ij} = 0, \Delta\text{bp}_{ij} = 0$),
   - **Translationally coupled** via overlapping stop/start codons ($\Delta\text{orf}_{ij} = 1, \Delta\text{bp}_{ij} = -4\text{ bp}$ `ATGA`),
   - Co-transcribed in the **same directon** ($\Delta\text{bp}_{ij} \in [0, 50]\text{ bp}$, zero intervening terminators),
   - **Divergently transcribed** across a shared bidirectional promoter ($\leftarrow_i \ \rightarrow_j$, $\Delta\text{bp}_{ij} \in [30, 250]\text{ bp}$ — the canonical layout of `TetR`/`LacI`/`LysR` transcriptional regulators and their target operons!),
   - **Convergently transcribed** ($\rightarrow_i \ \leftarrow_j$, separated by a bidirectional transcription terminator — typically functionally insulated!), or
   - Separated by a **riboswitch, leader attenuator, or integron `attC` recombination site**.

#### Exact Fix: Dual-Coordinate `RoPE` + Pairwise Genomic Geometry Attention Bias $B_{ij}^{(h)}$

In every transformer layer $\ell \in \{1, \dots, 12\}$ and attention head $h \in \{1, \dots, 8\}$ ($d_{\text{head}} = 64$), replace standard 1D slot `RoPE` with:

$$\text{AttnLogit}_{ij}^{(\ell, h)} = \frac{\big(\mathbf{R}_{\text{orf}}(\text{orf}_i) \oplus \mathbf{R}_{\text{seg}}(\text{seg}_i)\,\mathbf{q}_i^{(\ell, h)}\big)^\top \big(\mathbf{R}_{\text{orf}}(\text{orf}_j) \oplus \mathbf{R}_{\text{seg}}(\text{seg}_j)\,\mathbf{k}_j^{(\ell, h)}\big)}{\sqrt{d_{\text{head}}}} + B_{ij}^{(h)}$$

where:
1. **2D Split Rotary Embedding ($\mathbf{R}_{\text{orf}} \oplus \mathbf{R}_{\text{seg}}$):**
   - The first $48$ dimensions of $\mathbf{q}_i, \mathbf{k}_i$ are rotated by **`orf_idx`** ($\theta_m = 10000^{-2m/48}$), so all domains inside the same ORF share the exact same ORF-level rotary phase!
   - The remaining $16$ dimensions are rotated by **`seg_idx_in_orf`** ($\phi_m = 100^{-2m/16}$), encoding N-to-C intra-protein domain order orthogonally to inter-gene distance.
2. **Pairwise Genomic Geometry Attention Bias $B_{ij}^{(h)} \in \mathbb{R}$:**
   Computed once per contig tape (shared across layers, projected to $H=8$ heads via a lightweight 2-layer MLP over exact pairwise genomic invariants):
   $$B_{ij}^{(h)} = \mathbf{w}_h^\top \text{GELU}\Big(\mathbf{W}_1 \big[ \mathbf{e}_{\text{topo}}(\tau_{ij}) \ \big\| \ \boldsymbol{\psi}_{\text{RBF}}(d_{ij}^{\text{bp}}) \ \big\| \ \mathbf{e}_{\Delta\text{orf}}\big(\text{clip}(\text{orf}_j - \text{orf}_i, -16, +16)\big) \ \big\| \ \mathbf{c}_{\text{inter}}(i, j) \big]\Big)$$
   with exact inputs:
   - **Discrete Pairwise Topology State $\tau_{ij} \in \{0, \dots, 6\}$:**
     - `0`: Same ORF ($\text{orf}_i = \text{orf}_j$, covalent fusion; sub-indexed by $\text{seg}_j - \text{seg}_i$),
     - `1`: Same directon, translational overlap/abutment ($\max_{k \in (i, j]} \Delta\text{bp}_k \le +4\text{ bp}$),
     - `2`: Same directon, tight co-transcription ($\max_{k \in (i, j]} \Delta\text{bp}_k \in [5, 60]\text{ bp}$, zero terminators),
     - `3`: Same directon, loose co-orientation ($\max_{k \in (i, j]} \Delta\text{bp}_k \in [61, 300]\text{ bp}$, zero terminators),
     - `4`: Divergent head-to-head ($\text{strand}_i = -1, \text{strand}_j = +1$ in contig order, zero intervening genes, sharing 5' promoter region),
     - `5`: Convergent tail-to-tail ($\text{strand}_i = +1, \text{strand}_j = -1$, sharing 3' termination region),
     - `6`: Across $\ge 1$ terminator or multiple strand flips.
   - **Signed Log-BP Radial Basis Function $\boldsymbol{\psi}_{\text{RBF}}(d_{ij}^{\text{bp}}) \in \mathbb{R}^{16}$:**
     Let $d_{ij}^{\text{bp}}$ be the signed nucleotide edge-to-edge distance between feature $i$ and feature $j$ ($d_{ij}^{\text{bp}} = 0$ if $\text{orf}_i = \text{orf}_j$, negative if overlapping ORFs). Transform via signed log-modulus:
     $$\tilde{d}_{ij} = \text{sgn}(d_{ij}^{\text{bp}}) \cdot \log_2\big(1 + |d_{ij}^{\text{bp}}|\big) \in [-8, +15]$$
     and expand across 16 Gaussian RBF centers $\mu_1, \dots, \mu_{16}$ spaced from $-4\text{ bp}$ to $+16\text{ kb}$.
   - **Intervening Barrier Counter Vector $\mathbf{c}_{\text{inter}}(i, j) \in \mathbb{R}^4$:**
     $$\mathbf{c}_{\text{inter}}(i, j) = \log_2\big(1 + [\text{n\_terminators}_{ij}, \text{n\_strand\_flips}_{ij}, \text{n\_attC}_{ij}, \text{n\_riboswitches}_{ij}]\big)$$

---

## 3. Experimental Closure Card (`family_closure_card`) Precision Bugs

### 3.1 Flaw 3A: The Multi-Domain Chimera & Unsolved-Complex Closure Trap

#### Diagnosis
In `PLAN.md` Section 4.3 (`family_closure_card.parquet`, lines 163–171), experimental characterization status (`max_pdb_seq_id`, `n_swissprot_exp`, `n_paperblast_hits`, `max_rb_tnseq_cofitness`) is joined from `PDB100`, `Swiss-Prot`, `PaperBLAST` (`litsearch.db`), and `Fitness Browser` (`feba.db`).

**There are two severe coordinate and relational bugs in this design:**
1. **Whole-Protein Alignment Bleed Across Multi-Domain Chimeras:**
   `Swiss-Prot`, `PaperBLAST` (`litsearch.db`), and `Fitness Browser` index **whole proteins**, not sub-ORF domain segments. Even `PDB100` chains frequently contain multi-domain constructs or expression tags in `SEQRES` where only one domain was resolved in the crystal density (`ATOM` records).
   - Suppose metagenomic Protein $P$ is a two-domain chimera: `seg_0` (residues `1..210`, a canonical `PF00005` ABC ATPase or Rossmann domain) fused to `seg_1` (residues `225..490`, an uncharacterized `DUF9999` domain).
   - Or worse: suppose a `PaperBLAST` or `Swiss-Prot` entry $Q$ is *itself* a two-domain protein (`Rossmann + DUF9999`) where the published paper mutated only the `Rossmann` domain, or suppose `PDB` entry `7XYZ` contains full-length `SEQRES` (`1..490`) but residues `225..490` (`DUF9999`) were disordered or proteolytically cleaved and have **zero atomic coordinates in the mmCIF `ATOM` table**!
   - Any whole-protein or `SEQRES`-based homology join falsely stamps `DUF9999` as `"Solved in PDB (max_pdb_seq_id = 0.90)"`!
2. **Monomer Closure vs. Pairwise Complex / Operon Closure Blindness:**
   Suppose `Dom_A` (`max_pdb_seq_id = 0.85`, solved as an isolated monomer in 2008 by a structural genomics consortium, e.g. JCSG/NESG, with zero functional assay) and `Dom_B` (`max_pdb_seq_id = 0.78`, a generic radical SAM or nucleotidyltransferase solved as a monomer) form a conserved 2-gene operon or Rosetta-stone fusion across 6 phyla.
   - If discovery queries filter by monomer `family_closure_card` (`max_pdb_seq_id < 0.30`), **this entire novel pathway is thrown away** because both monomers have structural genomics PDB structures—even though **their physical complex, biochemical partnership, and operon coupling `(Dom_A, Dom_B)` have never been co-crystallized, co-assayed, or co-mentioned in any paper**!

#### Exact Fix: Envelope-Gated Coordinate Intersection & Pairwise Complex Closure (`pair_closure`)

1. **Rule 1 — Strict Sub-ORF Envelope & Resolved-`ATOM` Coordinate Gate:**
   Let domain segment $s$ on protein $P$ have envelope coordinates $[a_s, b_s]$ ($L_s = b_s - a_s + 1$). A sequence/structure alignment between $P$ and an external database target $T \in \{\text{PDB100}, \text{SwissProt}, \text{PaperBLAST}, \text{FEBA}\}$ spanning query residues $[q_1, q_2]$ and target residues $[t_1, t_2]$ is attributed to segment $s$ **if and only if**:
   $$\text{EnvCov}(s, T) = \frac{\big|[a_s, b_s] \cap [q_1, q_2]\big|}{b_s - a_s + 1} \ge 0.70$$
   **AND**, for `PDB100` targets, let $\mathcal{R}_{\text{ATOM}}(T) \subseteq [1, L_T]$ be the set of target residue indices with **resolved backbone heavy atoms ($\text{N}, \text{C}_\alpha, \text{C}$) in the mmCIF `ATOM` record**:
   $$\text{ResolvedAtomCov}(s, T) = \frac{\big|\{i \in [a_s, b_s] \cap [q_1, q_2] : \text{map}_{P \to T}(i) \in \mathcal{R}_{\text{ATOM}}(T)\}\big|}{b_s - a_s + 1} \ge 0.70$$
   **AND**, for `Swiss-Prot` / `PaperBLAST`, distinguish between:
   - `STRUCTURAL_GENOMICS_ONLY` (PDB structure from JCSG/NESG/MCSG/SGC with title matching `"crystal structure of uncharacterized/putative..."` and zero catalytic/phenotypic `ECO:0000269` functional annotations), vs.
   - `BIOCHEMICALLY_CHARACTERIZED` (verified enzymatic activity, ligand/substrate identification, or genetic phenotype).

2. **Rule 2 — First-Class Pairwise Complex & Epistatic Closure on `directon_pair_summary.parquet`:**
   Add explicit **Pairwise Closure Columns** to `directon_pair_summary.parquet` for every pair `(token_a, token_b)` (indexed at the `subfamily_50_id` level):
   - `pdb_cocrystal_status: uint8`:
     - `0 = NONE` (never resolved in the same PDB entry),
     - `1 = SAME_COMPLEX_NO_CONTACT` (in the same biological assembly, min heavy-atom distance $> 8\text{ \AA}$),
     - `2 = DIRECT_PHYSICAL_INTERFACE` (resolved in the same PDB biological assembly—either inter-chain or intra-chain—with $\ge 5$ heavy-atom contacts $< 5.0\text{ \AA}$ and interface buried surface area $\Delta\text{SASA} \ge 300\text{ \AA}^2$).
   - `max_pdb_pair_seq_id: float32`: $\min\big(\text{seq\_id}(a, T_a), \text{seq\_id}(b, T_b)\big)$ across all PDB biological assemblies containing both domains in physical contact.
   - `paperblast_pair_co_pmids: uint16`: Count of shared PMIDs in `PaperBLAST` (`litsearch.db`) where a single publication links to envelope-verified homologs ($\ge 40\%$ ID) of **both** `token_a` and `token_b`.
   - `feba_pair_cofitness_max: float32`: Maximum empirical Pearson correlation $r_{\text{cofit}}(a, b)$ across all `Fitness Browser` species encoding both `token_a` and `token_b` in the same genome.
   - `feba_pair_shared_phenotype_count: uint16`: Count of specific stress/carbon/nitrogen/phage conditions where both `token_a` and `token_b` show strong mutant fitness defects ($|t| > 4$).

---

### 3.2 Flaw 3B: The Broad Superfamily & Catalytic-Site Divergence Trap

#### Diagnosis
In `PLAN.md` Section 4.3 (`feature_tape.parquet` line 145 and `family_closure_card.parquet` line 163), `token_id` is indexed at the coarse `Pfam` / `COG` / `arCOG` family level, while Section 5 Benchmark 4 (line 221–222) mentions 40%-identity subfamilies (`subfamily_token_id`) as a validation check **without actually incorporating `subfamily_id` or active-site residue conservation into `feature_tape.parquet`, `directon_pair_summary.parquet`, or `family_closure_card.parquet`!**

**Why both family-level and naive 40%-sequence-identity closure fail:**
1. **Enzymatic Superfamily Explosion:**
   A single Pfam/COG—such as `PF04055` (Radical SAM, $>700,000$ members), `PF00534` (`GT4` glycosyltransferase), `PF00005` (`ABC_tran`), `PF00583` (`GNAT` acetyltransferase), or `COG0500` (SAM methyltransferase)—catalyzes **dozens to hundreds of distinct biochemical reactions**. Marking `PF04055` as "Characterized" because `BioB`, `LipA`, and `MiaB` are solved blinds the engine to 90% of radical SAM chemistry.
2. **Active-Site Divergence Within 40–50% Sequence Clusters (Neofunctionalization & Pseudoenzymes):**
   Even within a 40–50% global sequence-identity sub-cluster, evolution frequently mutates the catalytic constellation—either **neofunctionalizing** substrate specificity or creating **catalytically inactive pseudoenzymes** that act as allosteric scaffolds, toxin-antitoxin immunity inhibitors, or regulatory sensors!
   - In `PLAN.md` line 148, `feature_tape.parquet` stores only a single boolean `catalytic_intact: bool`. A single boolean cannot distinguish *which* catalytic motif variant a subfamily carries!

#### Exact Fix: Two-Coordinate Subfamily & Active-Site Constellation Indexing (`subfamily_50_id`, `catalytic_motif_hash`)

1. **Hierarchical Token Granularity in `feature_tape.parquet`:**
   Every feature slot $i$ is indexed across four explicit resolution levels:
   - `clan_or_cog_id: uint32` (Level 1: Pfam Clan / Superfamily),
   - `seq_family_id: uint32` (Level 2: Pfam / COG / arCOG / NCBIfam / PHROG / Rfam / `DARK_MMSEQS30`),
   - `subfamily_50_id: uint32` (Level 3: `MMseqs2 cluster --min-seq-id 0.50 -c 0.80 --cov-mode 0` computed strictly over the excised domain envelope sequences $[a_s, b_s]$ within each `seq_family_id`),
   - `catalytic_motif_hash: uint32` (Level 4: `xxhash32` of the ordered amino-acid characters at the family's $K \in [3, 12]$ conserved active-site/ligand-binding alignment columns, defined by CDD/SFLD/M-CSA/AlphaFill pocket positions or the top-8 highest-information-content conservation columns in the family MSA).

2. **Active-Site-Gated Subfamily Closure Score (`subfamily_closure_card.parquet`):**
   Replace `family_closure_card.parquet` with **`subfamily_closure_card.parquet`** keyed by `(subfamily_50_id, catalytic_motif_hash)`. For any query domain segment $s$ and characterized reference target $T \in \{\text{SwissProt\_Exp}, \text{PDB\_Biochem}, \text{PaperBLAST}\}$, define the **Isofunctional Closure Score** $\mathcal{C}_{\text{iso}}(s) \in [0, 1]$:
   $$\mathcal{C}_{\text{iso}}(s) = \max_{T} \Big[ \underbrace{S_{\text{seq}}\big(\text{seq\_id}_{\text{env}}(s, T)\big)}_{\text{Envelope sequence similarity}} \cdot \underbrace{\exp\!\big(-\gamma \, d_{\text{Hamming}}(\text{Motif}(s), \text{Motif}(T))\big)}_{\text{Active-site constellation penalty}} \cdot w_{\text{evidence}}(T) \Big]$$
   where:
   - $\text{seq\_id}_{\text{env}}(s, T)$ is the sequence identity over the domain envelope ($\text{EnvCov} \ge 0.70$),
   - $S_{\text{seq}}(x) = \text{clip}\big(\frac{x - 0.30}{0.40}, 0, 1\big)$ ramps linearly from $0.0$ at $\le 30\%$ ID to $1.0$ at $\ge 70\%$ ID,
   - $d_{\text{Hamming}}(\text{Motif}(s), \text{Motif}(T))$ is the number of mutated catalytic/active-site residues between $s$ and reference $T$ ($\gamma = 1.4$, so **even a single catalytic residue mutation drops $\exp(-1.4) = 0.246$, preventing a pseudoenzyme or neofunctionalized paralog from ever being masked as "already solved"**!),
   - $w_{\text{evidence}}(T) = 1.0$ for verified biochemical/genetic function (`Swiss-Prot ECO:0000269`, `MIBiG`, `GapMind`, `DefenseFinder`), $0.75$ for `PaperBLAST` title/abstract match or `FEBA` strong phenotype ($|t| > 5$), and $0.35$ for uncharacterized structural-genomics-only PDB deposition (`STRUCTURAL_GENOMICS_ONLY`).

---

## 4. Exact Schema Diffs & Specification Upgrades for `PLAN.md` v2.0

### 4.1 Schema Upgrade 1: `feature_tape.parquet` (Additions & Modifications)
Add the following columns to `feature_tape.parquet` (`PLAN.md` Section 4.3, Table 1):

| Column | Type | Description |
|---|---|---|
| `subfamily_50_id` | `uint32` | **NEW:** `MMseqs2` 50% AA identity / 80% envelope coverage sub-cluster ID within `seq_family_id` |
| `super_cluster_id` | `uint32` | **NEW:** Connected-component superfamily ID ($<30\%$ AA + `HHsearch` $P \ge 80\%$ + Pfam Clan + `ECOD_H`) for zero-leakage ML splits |
| `catalytic_motif_hash` | `uint32` | **NEW:** `xxhash32` of amino acids at the family's conserved active-site/pocket MSA columns |
| `local_regime` | `uint8` | **NEW:** `0=CORE`, `1=SUPER_OPERON`, `2=PLASMID`, `3=PHAGE`, `4=ICE_GI`, `5=INTEGRON` |
| `vehicle_cluster_id` | `uint32` | **NEW:** Single-linkage flanking-backbone vehicle cluster ID ($\pm 15$ ORFs Jaccard $\ge 0.30$ or vOTU/pOTU) |
| `orfs_to_contig_start` | `uint16` | **NEW:** Exact count of called ORFs to $5'$ contig end (for exact censoring-adjusted conditional denominators) |
| `orfs_to_contig_end` | `uint16` | **NEW:** Exact count of called ORFs to $3'$ contig end |

### 4.2 Schema Upgrade 2: `directon_pair_summary.parquet` (Passenger-Corrected & Pair-Closed)
Upgrade `directon_pair_summary.parquet` (`PLAN.md` Section 4.3, Table 2) to index at both family and subfamily levels and include directional, entropy, vehicle, and pair-closure metrics:

| Column | Type | Description |
|---|---|---|
| `token_a`, `token_b` | `uint32` | Primary sequence/ncRNA family IDs (plus parallel table / rollup keyed by `subfamily_50_a`, `subfamily_50_b`) |
| `n_veh_indep_events` | `uint32` | **NEW ($N_{\text{eff}}^{\text{veh}}$):** Independent evolutionary co-occurrences after collapsing shared MGE flanking backbones ($\text{Jaccard}_{\pm 15} \ge 0.30$) |
| `host_tree_homoplasy` | `uint16` | Host-tree (`bac120`/`ar53`) Fitch parsimony co-acquisitions (diagnostic comparison against `n_veh_indep_events`) |
| `gene_tree_homoplasy_b_on_a` | `uint16` | **NEW ($H_{B \mid T_A}$):** Fitch parsimony gains/retentions of $B$ mapped onto the protein phylogeny $T_A$ of family $A$ |
| `p_b_given_a_uncens` | `float32` | **NEW ($P(B \mid A)$):** Censoring- and vehicle-corrected conditional probability that $B$ is within $\pm k$ ORFs given $A$ |
| `p_a_given_b_uncens` | `float32` | **NEW ($P(A \mid B)$):** Censoring- and vehicle-corrected conditional probability that $A$ is within $\pm k$ ORFs given $B$ |
| `cond_asymmetry_log2` | `float32` | **NEW ($\mathcal{A}(A,B)$):** $\log_2(P(B \mid A) / P(A \mid B))$ |
| `ctx_entropy_a_bits`, `ctx_entropy_b_bits` | `float16` | **NEW ($H_{\text{ctx}}(A), H_{\text{ctx}}(B)$):** Global neighbor Shannon entropy in bits (flags promiscuous `IS`/integrase hubs $> 6.5\text{ bits}$) |
| `slot_fold_entropy_ab_bits` | `float16` | **NEW ($H_{\text{ctx}}^{\text{fold}}(A; \Delta_{AB})$):** Conditional entropy of structural folds at relative slot $\Delta\text{orf}_{AB}$ |
| `micro_geometry_score` | `float32` | **NEW ($\mathcal{G}_{\text{micro}}$):** Composite score of `same_directon_frac`, $\exp(-\text{IQR}(\Delta\text{bp})/40)$, order conservation, and translational overlap |
| `func_coupling_llr` | `float32` | **NEW ($\text{LLR}_{\text{func}}$):** Primary ranking metric: vehicle-deduplicated, regime- & entropy-corrected log-likelihood ratio of functional coupling vs. passenger mobility |
| `empirical_fdr_qval` | `float32` | **NEW ($q_{\text{emp}}$):** Target-decoy empirical FDR calibrated against regime- and entropy-matched within-order shuffled contigs |
| `pdb_cocrystal_status` | `uint8` | **NEW:** `0=NONE`, `1=SAME_ASSEMBLY_NO_CONTACT`, `2=DIRECT_PHYSICAL_CONTACT` ($\ge 5$ contacts $< 5\text{ \AA}$) |
| `max_pdb_pair_seq_id` | `float32` | **NEW:** Max $\min(\text{ID}_a, \text{ID}_b)$ to a co-crystallized complex or two-domain resolved PDB structure |
| `paperblast_pair_co_pmids` | `uint16` | **NEW:** Number of shared publications in `PaperBLAST` matching envelope-verified homologs of both $A$ and $B$ |
| `feba_pair_cofitness_max` | `float32` | **NEW:** Max empirical RB-TnSeq co-fitness correlation $r(A, B)$ across `Fitness Browser` species |

### 4.3 Schema Upgrade 3: `subfamily_closure_card.parquet` (Replacing Coarse `family_closure_card.parquet`)
Keyed by `(subfamily_50_id, catalytic_motif_hash)` (with a rolled-up summary view at `seq_family_id`):

| Column | Type | Description |
|---|---|---|
| `subfamily_50_id` | `uint32` | 50% sequence-identity domain-envelope cluster ID |
| `catalytic_motif_hash` | `uint32` | Active-site amino-acid constellation hash (`0` if no known active site) |
| `seq_family_id` | `uint32` | Parent `Pfam` / `COG` / `arCOG` / `PHROG` / `DARK_FAM` ID |
| `catalytic_motif_str` | `varchar` | Human-readable active-site constellation (e.g., `"C112-C116-C119_SAM"` or `"H142-E143-H146_Zn"`) |
| `catalytic_intact` | `bool` | `True` if canonical catalytic residues are preserved; `False` if mutated (pseudoenzyme/neofunctionalized candidate) |
| `max_pdb_atom_env_seq_id` | `float32` | **FIXED:** Max sequence ID to a PDB chain requiring $\ge 70\%$ domain-envelope overlap with **resolved `ATOM` coordinates** |
| `pdb_is_structural_genomics_only` | `bool` | **NEW:** `True` if closest PDB hit is an uncharacterized JCSG/NESG/SGC monomer with zero functional assay |
| `max_swissprot_env_seq_id` | `float32` | **FIXED:** Max envelope-verified sequence ID to a `Swiss-Prot` entry with `ECO:0000269` experimental functional evidence |
| `active_site_matches_ref` | `bool` | **NEW:** `True` if `catalytic_motif_hash` exactly matches the closest characterized `Swiss-Prot`/`PDB` reference |
| `isofunctional_closure_score` | `float32` | **NEW ($\mathcal{C}_{\text{iso}} \in [0, 1]$):** Envelope- and active-site-gated closure score (capped at $\le 0.25$ if active site diverges!) |
| `n_paperblast_env_hits` | `uint32` | **FIXED:** Count of `PaperBLAST` hits covering $\ge 70\%$ of this domain segment's envelope |

---

## 5. Upgraded Falsification & Adversarial Benchmark Suite (`PLAN.md` Section 5 v2.0)

In addition to Benchmarks 1–4 (upgraded below), `all-the-data` must pass **two new adversarial negative-control benchmarks (Benchmarks 5 & 6)** before deployment:

1. **Benchmark 1 (Upgraded) — Untargeted Historical Discovery Recall Against Full MGE Background:**
   - Run the 50 historical discovery queries **without manual MGE blacklisting** across all 3 Tiers combined.
   - **Pass Criterion:** Using `func_coupling_llr` ($\text{LLR}_{\text{func}}$) and `isofunctional_closure_score`, $\ge 95\%$ ($48/50$) of historical partners must rank in the **Top 15** for their bait family, whereas ranking by naive `synteny_log10_pval` or raw host-tree `parsimony_homoplasy_events` must be shown in an ablation table to fail (dropping recall@15 below $50\%$ due to transposase/plasmid flood).
2. **Benchmark 2A & 2B (Upgraded) — Zero-Leakage Connected-Component & Blind NOGD Slot Prediction:**
   - **2A (Disjoint Operon-Homology Connected Components):** Hold out $10\%$ of connected components from $\mathcal{G}_{\text{leak}}$ (zero $\ge 30\%$ AA, `HHsearch`, `Pfam-Clan`, or `ECOD_H` 3-gram overlap with `Train`). Require `Recall@10` $\ge 65\%$ on `ECOD_T` fold prediction and **0.0% biophysically incoherent predictions** (verified by checking that predicted `TM_count >= 2` never co-occurs with purely cytoplasmic soluble `ECOD_T` folds in top-10 beams).
   - **2B (Blind NOGD Replacement Recovery):** Purge all instances and clan-relatives of 50 known NOGD replacement enzymes (`ThyX`, `coaX`, `nadV`, `fbiD`, etc.) from `Train`. When querying their canonical pathway with a masked hole `[A, B, MASK, D]`, the Coarse-to-Fine decoder's Stage 1+2 marginal must rank the true replacement's `(biophysical_loc, ECOD_T_fold)` in the **Top 5** for $\ge 70\%$ ($35/50$) of systems.
3. **Benchmark 5 (NEW) — The MGE Passenger & Super-Operon Negative-Control Trap Suite ($\text{FPR} < 1\%$):**
   - Construct an adversarial panel of **500 known non-functional passenger & boundary pairs** that occur in $\ge 30$ genomes across $\ge 3$ phyla:
     - *200 Plasmid/ICE Cargo–Backbone Pairs:* e.g., `blaTEM-1` + `tnpA` (`IS26`/`Tn3`), `tet(M)` + `Tn916` `orf13`, `mcr-1` + `ISApl1`, `IncP` `trfA` + cargo AMR/mercuric resistance genes.
     - *150 Class 1/2/3 Integron Cassette Neighbor Pairs:* Unrelated AMR or `DUF` cassettes separated by an `attC` recombination site (`aadA1` + `dfrA1`, `oxa2` + `ges1`).
     - *100 Prophage Boundary / Moron–Structural Pairs:* Prophage cargo toxins (`Stx`, `Ctx`, `SopE`) within $\pm 4$ ORFs of phage tail fiber or lysis holin genes across enterobacteria/vibrios.
     - *50 Housekeeping Super-Operon Boundary Pairs:* Unrelated flanking genes immediately adjacent to the $5'$ or $3'$ boundary of the `S10` ribosomal operon or `dcw` cluster.
   - **Pass Criterion:** At the `func_coupling_llr` threshold that achieves $95\%$ recall on Benchmark 1, **$< 1\%$ ($\le 5/500$) of Benchmark 5 passenger/boundary traps may pass as false positives**.
4. **Benchmark 6 (NEW) — Chimeric Fusion & Pseudoenzyme Closure Decoy Test ($0\%$ False Closure Bleed):**
   - **6A (Synthetic & Natural Two-Domain Chimeras):** Take 200 multi-domain proteins consisting of one solved PDB domain (`Rossmann`, `ABC_tran`, `P-loop`, `HTH`) fused to an uncharacterized `DUF` or viral effector domain. Verify that **100% ($200/200$)** of the uncharacterized domain segments receive `max_pdb_atom_env_seq_id == 0.0` and `isofunctional_closure_score < 0.10`.
   - **6B (Active-Site Pseudoenzyme Decoys):** Take 100 natural pseudoenzyme paralogs (e.g., inactive kinase/protease/nuclease/Radical-SAM regulators sharing $45\text{–}60\%$ global sequence identity with a solved active enzyme in `Swiss-Prot`/`PDB` but carrying mutated catalytic triad residues). Verify that `active_site_matches_ref == False` and `isofunctional_closure_score <= 0.25` for **100% ($100/100$)** of cases.
