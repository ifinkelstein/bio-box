# Audit Report 5: Data Science, Statistical Inference, Calibration & Validation Design

**Auditor:** Applied statistician / ML practitioner (ranking models on biological data; calibration, FDR, leakage, power)
**Target Documents:** `00-spec.md`, `PLAN.md` (v2.0), `DEPLOYMENT-AND-MALLEABILITY.md`, `PRIOR-ART.md`, `audits/critique_1..3`, `audits/SYNTHESIS.md`
**Verdict:** The v2.0 plan is an engineering specification wearing the clothes of a statistical one. Every quantity that will decide which hypotheses a human looks at — `func_coupling_llr`, `empirical_fdr_qval`, `isofunctional_closure_score`, `Recall@10` — is **defined by assertion rather than by a model, has no estimation procedure, no calibration target, no power analysis, and no control that could fail**. Two of the headline numbers are arithmetically wrong (the 95% ANI sketch test, the "90% AA-ID invariant" Stage-0 hash), one headline guarantee is structurally impossible (zero-leakage single-linkage splitting at planetary scale), and the entire discovery benchmark is circular (the systems to be "rediscovered" are the sources of the HMM libraries used as input). All are fixable, and the fixes are cheaper than the architecture already specified.

---

## Step A — Deepening the audit task before executing it

**(1) What a top professional would do that was not asked for.**
Not "list statistical objections." They would (a) write down the *generative model* the plan implicitly assumes and show where the plan's statistic is not the likelihood ratio of that model; (b) **do the arithmetic** on every numeric claim that a statistician can check from the page alone (sketch resolution, detection floor, multiplicity, CI width on the pass criteria) and report where it breaks; (c) replace assertion-thresholds with an *estimation procedure* (fit the weights on a labeled set, lock, evaluate out-of-sample); (d) design the **controls that could falsify the system**, especially the negative/ablation controls that the plan omits — a context-shuffled control for the language model and a frequency-prior baseline, without which `Recall@10 ≥ 65%` is uninterpretable; (e) insist on a **pre-registration + temporal holdout** protocol, because a retrospective "recall of known discoveries" on a database built from those discoveries measures nothing; (f) specify what **sufficient statistics must be persisted** so the score can be recomputed under new parameters without a 4.7-day rebuild — a statistical requirement that is also a storage-layout requirement.

**(2) The rawest material they would start from.**
Not the three existing audits (already-digested conclusions), and not the aggregate tables the plan proposes. The raw material is: the **primary methods papers** for each borrowed statistic (Li & König b-bit minwise hashing; Ondov/Mash ANI–Jaccard identity; Shaw & Yu skani; Efron's empirical null/`locfdr`; Storey π0; Yekutieli hierarchical FDR; Strona's curveball null; Pagel's correlated-evolution test; Teufel's GraphPart); a **1%-sample pilot build** of real `window_postings` rows from a few thousand GTDB genomes plus one SPIRE biome, which is the only way to learn the component-size distribution, the score null distribution, and the dedup collapse ratio; the **date-stamped historical database snapshots** (Pfam 27.0, COG 2014, Rfam 11, PDB/SwissProt/PaperBLAST frozen at a cut date) that make a real temporal holdout possible; and **experimentally verified operon/complex ground truth** (RegulonDB/ODB verified operons, M-CSA catalytic-residue mutagenesis, PDB biological assemblies, RB-TnSeq cofitness) rather than HMM-library labels.

**(3) Specific methods and tools.**
Poisson-binomial score test / conditional logistic (Cochran–Mantel–Haenszel stratified by regime and clade) for the coupling statistic; cluster-bootstrap over vehicle clusters for variance; Efron empirical-null local fdr + Storey π0 + parametric tail with permutation validation in the reachable range; curveball/fastball degree-preserving swaps and block-circular permutation with block length > window; Pagel/BayesTraits-style correlated-evolution LRT or Mk/Dollo transition modeling instead of bit-shift adjacency; Beta-binomial shrinkage for pangenome plane calls; monotone-constrained gradient boosting (XGBoost `monotone_constraints`) or regularized logistic regression to *learn* `λ_hub, λ_attC, α_core, γ, β_couple` instead of asserting them; Sobol/Morris global sensitivity analysis over a `params.yaml`; temperature scaling + reliability diagrams + ECE + risk–coverage curves for `gLM-Domain`; GraphPart-style balanced homology partitioning instead of single-linkage; Zombi/ALF genome-evolution simulation with HGT and rearrangement to generate data with **known** coupling truth; skani for ANI verification; Coinfinder as an external baseline method to beat.

**(4) What, other than the auditor or the author, would show the result is right.**
Four things, none of which is an opinion: (i) **simulation with known ground truth** — Zombi/ALF genomes with planted co-evolving modules and planted promiscuous MGE passengers, where observed FDR can be compared to nominal `q_emp` on a reliability diagram; (ii) **temporal holdout** — a pre-2016 token dictionary recovering post-2016 published discoveries, which no amount of tuning can fake; (iii) **orthogonal experimental evidence not used in fitting** — AlphaFold3/Boltz interface pTM on top-ranked pairs, RB-TnSeq cofitness released after the freeze date, co-expression from independent metatranscriptomes, PDB assemblies deposited after the freeze; (iv) **beating a published baseline** (Coinfinder, STRING's own scores, gLM2) on a shared, pre-registered task. Plus a pre-registered prediction registry with a BLAKE3-hashed timestamp, so the hit rate is auditable by a third party later.

**Rewritten prompt.**
> Audit this plan as the statistician who will have to defend its numbers to a skeptical reviewer. Do not enumerate objections. Instead: (1) write down the generative model the plan implicitly assumes for pairwise gene coupling, and prove whether `LLR_func` is the likelihood ratio of that model; if it is not, derive the statistic that is, in a form computable in one streaming pass. (2) Check the arithmetic of every statistical claim checkable from the page: the 64-bit sketch's actual ANI operating point and its false-negative rate; the mutation-invariance of the Stage-0 architecture hash; the multiplicity correction implied by ~1e8–1e10 tested pairs and the resulting minimum detectable number of co-occurrences per mobility regime; the Wilson/exact-binomial confidence intervals on every "≥95%", "0%", and "100%" pass criterion. (3) Determine whether the zero-leakage single-linkage operon split can exist at this scale (percolation/giant-component argument) and specify the partitioning method that can. (4) Identify the circularity in the discovery benchmark and replace it with a dated temporal-holdout protocol naming the snapshot versions and per-system publication dates. (5) Specify the missing controls: frequency-prior baseline, context-shuffled ablation, calibration (ECE/reliability), power/`min_detectable_effect` per pair, and global sensitivity analysis over every hard-coded constant. (6) Specify what sufficient statistics must be persisted so any parameter change is a recompute, not a rebuild. (7) Name the external evidence — simulation with known truth, temporal holdout, orthogonal experiment, published baseline — that would show the result is right without the author or the auditor asserting it. Show all arithmetic; cite primary sources.

---

## Findings Summary

| ID | Topic | Severity | Primary location |
|---|---|---|---|
| **DS-01** | `LLR_func` is not a likelihood ratio; numerator undefined, units incoherent, scale confounded with sample size | **CRITICAL** | `PLAN.md` §4.3 Table 2C, `critique_3` §1.2C |
| **DS-02** | Seven free hyperparameters with no estimation procedure; ranking is hand-tuned | **CRITICAL** | `PLAN.md` §4.3, §4.3 Table 2D |
| **DS-03** | Stage-0 95% ANI test is arithmetically a ~99.3% ANI test with ~50% FNR; 64 bits is ~16× too few | **CRITICAL** | `PLAN.md` §2.3 Filter 4, `critique_2` §4.2 |
| **DS-04** | "90% AA-ID invariant" exact architecture hash is not mutation-invariant; collapse ratios and all storage estimates unfounded | **CRITICAL** | `PLAN.md` §2.3, `00-spec.md` Invariant 2 |
| **DS-05** | Benchmark 1 is circular: the defense/RiPP systems to be "rediscovered" are the provenance of the input HMM libraries | **CRITICAL** | `PLAN.md` §5 B1, B4 |
| **DS-06** | Zero-leakage single-linkage split collapses to one giant component; "0.0% leakage" unattainable and unverified | **CRITICAL** | `PLAN.md` §4.5(1), `critique_3` §2.1 |
| **DS-07** | Target-decoy FDR is a category error; permutation null cannot reach the 1e-10 tail required by multiplicity; no π0, no local fdr, no hierarchy correction across 4 token levels | **CRITICAL** | `PLAN.md` §4.3, §5 |
| **DS-08** | No power analysis; statistic is powerless exactly on small replicons (Tier 2 MGEs), where p0 → 1 | **MAJOR** | `PLAN.md` §4.3, `critique_3` §1.2A |
| **DS-09** | Non-independence beyond vehicles: pair-within-operon, BioProject/study pseudo-replication, sampling-effort bias in `arch_occurrence_count` | **MAJOR** | `PLAN.md` §2.3, §4.2 |
| **DS-10** | Pass criteria stated as 0%/100%/≥95% on n=50–500 without confidence intervals; traps not pre-registered or blinded | **MAJOR** | `PLAN.md` §5 B1–B6 |
| **DS-11** | No ML controls: no frequency-prior baseline, no context-shuffled ablation, no calibration, no frequency-stratified metrics | **MAJOR** | `PLAN.md` §4.5, §5 B2A |
| **DS-12** | Annotation uncertainty (HMM E-values, predicted 3Di, gene calls) treated as ground truth; no propagation, no measured ProstT5→ECOD error rate on novel folds | **MAJOR** | `PLAN.md` §3.2c, §4.2 |
| **DS-13** | Pangenome plane thresholds (≥50%, <30%) undefined at small n; bit-shift "phylogenetic coherence" is not a test | **MAJOR** | `PLAN.md` §4.4, `critique_1` BIO-08 |
| **DS-14** | Hard-coded constants (0.30, 1.4, 0.70, ±6, 16, 6.5 bits, 0.95) with no sensitivity analysis | **MAJOR** | throughout |
| **DS-15** | No simulation-with-known-truth tier; nothing can measure actual FDR | **MAJOR** | `PLAN.md` §5 |
| **DS-16** | Aggregates destroy the sufficient statistics needed to recompute scores; "malleability" covers new columns but not new statistics | **MAJOR** | `DEPLOYMENT-AND-MALLEABILITY.md` §3 |
| **DS-17** | No external baseline comparison; no prospective prediction registry | **MODERATE** | `PLAN.md` §5, `PRIOR-ART.md` |

---

## DS-01 — `LLR_func` is not a likelihood ratio

**Severity: CRITICAL.** `PLAN.md` §4.3 Table 2C (`func_coupling_llr`), formula in `critique_3` §1.2 Component C.

The statistic is written as

$$\text{LLR}_{\text{func}}(A,B) = \sum_{m} \log_{10} \frac{P_{\text{obs}}(A,B \mid \Delta\text{bp}_m, \text{strand}_m)}{p_m^{(0)}(B \mid A, r_m)\, p_{\text{geom}}^{(0)}} - \lambda_{\text{hub}}\max(0, H^{\text{fold}}_{\text{ctx}} - H_0) - \lambda_{\text{attC}} f_{\text{attC}}$$

Four defects, each independently fatal to interpretation:

1. **The numerator is never defined.** `P_obs` is called an "observed probability" of an event that was observed — i.e. 1 at every locus in the sum, in which case the statistic is $-\sum_m \log_{10} p_m^{(0)}$, a monotone function of $N_{\text{eff}}^{\text{veh}}$ and nothing else. Any alternative must have a **free parameter estimated from the data**; there is none. As written there is no alternative hypothesis, so there is no likelihood ratio.
2. **Scale is confounded with sample size.** Because the sum runs over loci, a weakly coupled pair seen in 10,000 vehicles outranks an obligate pair seen in 30. Effect size and evidence strength are collapsed into one number, which is precisely the failure mode the audit claims to fix for hypergeometric p-values. A ranked list sorted by this quantity is a list sorted mostly by abundance.
3. **Units are incoherent.** A $\log_{10}$ Bayes factor has additive penalties in **bits** (`H_ctx`) and in **unitless fractions** (`f_attC`) subtracted from it. $\lambda$ cannot repair this; it only hides it.
4. **Hits are not independent even after vehicle dedup** (see DS-09).

**Fix.** Specify the generative model, then use its score test — which is streaming-computable from exactly the accumulators `SYS-06` already maintains.

Let $y_m \in \{0,1\}$ indicate that $B$ falls within $\pm k$ ORFs of $A$ at independent locus $m$ (vehicle-deduplicated, uncensored). Null: $y_m \sim \text{Bernoulli}(p_m^{(0)})$ with $p_m^{(0)}$ from the regime- and clade-conditioned background. Alternative: a single coupling parameter $\beta$ on the logit scale,

$$\text{logit}\,\pi_m(\beta) = \text{logit}\,p_m^{(0)} + \beta .$$

Then report **three separate numbers** rather than one:

- **Effect size:** $\hat\beta$ = MLE of the conditional logistic model (interpretable as a log odds ratio of adjacency), with cluster-robust SE.
- **Evidence:** the Poisson-binomial / CMH score statistic
 $$Z_{AB} = \frac{\sum_m (y_m - p_m^{(0)})}{\sqrt{\sum_m p_m^{(0)}(1-p_m^{(0)})}}\;\xrightarrow{H_0}\;\mathcal{N}(0,1),$$
 whose two accumulators $\sum p_m^{(0)}$ and $\sum p_m^{(0)}(1-p_m^{(0)})$ are $O(1)$ state per group — the same streaming budget as now.
- **Geometry/promiscuity features kept separate:** $\mathcal{G}_{\text{micro}}$, $H^{\text{fold}}_{\text{ctx}}$, $f_{\text{attC}}$, regime fractions stored as columns and combined by a *fitted* model (DS-02), never by subtraction from a log-BF.

Store `beta_hat`, `beta_se`, `z_score`, `n_eff_veh`, `sum_p0`, `sum_p0q0`. Keep `func_coupling_llr` only as a derived view so existing text need not be rewritten wholesale.

**Validation.** On Zombi/ALF-simulated genomes with planted modules, $Z$ must be ~N(0,1) on unplanted pairs (QQ plot slope 1.00 ± 0.02 over the central 99%) and the rank correlation between $\hat\beta$ and planted coupling strength must exceed 0.8.

---

## DS-02 — Seven asserted hyperparameters, no estimation procedure

**Severity: CRITICAL.** `λ_hub`, `λ_attC`, `H_0`, `α_core`, `γ = 1.4`, `β_couple`, the Jaccard cut 0.30, the IQR scale 40 bp, and the ramp $S_{\text{seq}}$ endpoints 0.30/0.70 are all stated as constants. Nowhere does the plan say how they are chosen, on what data, or what would change them. A ranking function whose weights were chosen by its author and evaluated on a benchmark chosen by the same author has no out-of-sample guarantee at all.

`γ = 1.4` is the clearest case: it is reverse-engineered so that $e^{-1.4} = 0.2466$ lands just under the 0.25 pass line of the plan's own Benchmark 6B. That is fitting the parameter to the test.

**Fix.**
1. Put every constant in a single versioned `params.yaml` (hashed into `manifest.blake3.json`), with provenance for each: `asserted`, `fitted`, or `derived`.
2. Replace the hand-weighted combination with a **monotone-constrained model** (XGBoost `monotone_constraints`, or L2 logistic regression) over interpretable features $[\,\hat\beta, Z, N^{\text{veh}}_{\text{eff}}, \mathcal{G}_{\text{micro}}, H_{\text{ctx}}, H^{\text{fold}}_{\text{ctx}}, f_{\text{attC}}, \text{regime fractions}, \text{topology fractions}\,]$, trained on a **split** gold set: positives = experimentally verified co-functional pairs (RegulonDB/ODB verified operons, PDB biological assemblies with a real interface, MetaCyc/GapMind same-pathway adjacent steps, MIBiG same-BGC), negatives = the pre-registered trap panel. Fit on half, lock, evaluate on the other half, report AUPRC with a cluster-bootstrap CI over vehicle clusters.
3. Calibrate `γ` against labeled data instead of a target number: M-CSA plus Swiss-Prot `ECO:0000269` site-directed-mutagenesis annotations give thousands of mutations labeled activity-abolishing vs tolerated. Fit $\gamma$ by maximum likelihood on that set, and make the penalty **information-weighted and length-normalized**, since one mismatch out of a 3-residue triad is not one mismatch out of 12 columns:
 $$\exp\Big(-\gamma \textstyle\sum_j w_j \mathbb{1}[\text{mismatch}_j]\Big), \qquad w_j = \frac{\text{IC}_j}{\frac{1}{K}\sum_{j'} \text{IC}_{j'}}$$
 with $\text{IC}_j$ the information content of pocket column $j$ in the family MSA.

---

## DS-03 — The 95% ANI sketch test is wrong by ~4 percentage points of ANI and has ~50% FNR

**Severity: CRITICAL.** `PLAN.md` §2.3 Filter 4 step 2 and `critique_2` §4.2: `contig_ani_sketch: uint64` = "64-bit One-Permutation b-bit MinHash over DNA syncmers (k=15, s=5)"; "Check ≥ 95% ANI in 1 CPU cycle: `POPCNT(sketch_a ^ sketch_b) <= 6`."

**The arithmetic.** With one-permutation hashing into 64 bins at $b = 1$ bit per bin, the probability that two sketches agree in a given bin is the standard $b$-bit minwise collision probability
$$P_{\text{agree}} = J + \tfrac{1}{2}(1-J) = \tfrac{1+J}{2}, \qquad P_{\text{mismatch}} = \tfrac{1-J}{2}.$$
Converting 95% ANI to $k$-mer Jaccard at $k = 15$ via the Mash identity $J = \frac{r^k}{2 - r^k}$:
$$r^k = 0.95^{15} = e^{15\ln 0.95} = e^{-0.7694} = 0.4634, \qquad J = \frac{0.4634}{2 - 0.4634} = 0.3016.$$
So $P_{\text{mismatch}} = (1 - 0.3016)/2 = 0.3492$ and
$$\mathbb{E}[\text{POPCNT}] = 64 \times 0.3492 = \mathbf{22.3}, \qquad \text{sd} = \sqrt{64 \times 0.3492 \times 0.6508} = 3.81.$$
The threshold 6 sits $(6 - 22.3)/3.81 = -4.3$ sd below the mean of a genuine 95%-ANI pair. $P(\text{Binom}(64, 0.3492) \le 6) \approx 10^{-5}$. **A true 95% ANI pair passes this test essentially never.**

Inverting the threshold to find the ANI it actually encodes: $\mathbb{E}[\text{POPCNT}] = 6 \Rightarrow (1-J)/2 = 6/64 \Rightarrow J = 0.8125 \Rightarrow r^k = \frac{2J}{1+J} = 0.8966 \Rightarrow r = 0.8966^{1/15} = \mathbf{0.9927}$. The test is a **~99.3% ANI test**, not 95%. Worse, at that operating point sd $= \sqrt{64 \times 0.09375 \times 0.90625} = 2.33$, so a threshold placed at the mean gives a **~50% false-negative rate** even for pairs at the ANI it was designed for.

**Resolution requirement.** To separate 95% ANI ($J = 0.3016$, mismatch 0.3492) from 90% ANI ($0.9^{15} = 0.2059 \Rightarrow J = 0.1148$, mismatch 0.4426) with 3 sd of margin on each side:
$$|p_1 - p_2| = 0.0934 \ \ge\ 6\sqrt{\tfrac{p(1-p)}{n}} \ \Rightarrow\ n \ \ge\ \frac{0.2475}{(0.01557)^2} \approx \mathbf{1{,}021\ \text{bits}}.$$
So **≥1,024 bits (128 bytes)** per contig, not 64 bits — a 16× increase. This is affordable where it matters: only 2nd-and-later sightings are verified, so the RocksDB value grows from 12 bytes to ~132 bytes for the ~1.5×10⁸ recurrent hashes (~20 GB on NVMe, not RAM), and the Cuckoo-filter path is untouched.

**Fix.** Either (a) widen to a 1,024-bit sketch and recompute the threshold from the target ANI via the Mash identity rather than asserting an integer, deriving it as $\tau = n\frac{1-J^*}{2} - z_{\alpha}\sqrt{n p(1-p)}$ for an explicit one-sided error rate; or (b) **drop the fixed-width sketch for verification entirely and call `skani`** on candidate pairs — it was built for exactly this case (fragmented, incomplete assemblies, where plain sketching is known to underestimate ANI), is >20× faster than FastANI, and is reliable above ~82% ANI. Because only recurrent-hash candidates need verification (<1% of the stream by the plan's own numbers), the cost is negligible. Recommended: (b) for correctness, with (a) as the in-RAM prefilter.

Either way, `00-spec.md` Invariant 2 and `PLAN.md` §2.3 must stop claiming a 1-cycle POPCNT gives 95% ANI.

---

## DS-04 — The Stage-0 architecture hash is not mutation-invariant, so no collapse ratio in the plan is founded

**Severity: CRITICAL.** `PLAN.md` §2.3 Filter 4: `orf_sh_i` is a "32-bit Reduced-Alphabet Protein Syncmer Hash (Murphy-10, k=10, s=4 open syncmers → 90% AA ID invariant)", and `stage0_arch_hash = xxhash64(orf_sh_1..k, ctd_sh_1..k, strand_transitions, delta_bp_bins)`.

**Exact hashing cannot be identity-invariant.** `xxhash64` over a concatenation matches only on **byte-exact** equality of every input. For the composite hash to collide, every ORF's `orf_sh_i` must be bit-identical, every CTD hash identical, and all quantized spacings in the same bins. The claim "→ 90% AA ID invariant" would require each `orf_sh_i` to be a function only of features conserved at 90% identity. Open syncmers are not: density is $2/(k-s+1) = 2/7 = 0.286$, and one substitution perturbs every syncmer covering it, ≈ $10 \times 0.286 = 2.9$ syncmers per substitution.

For a 300-aa ORF at 90% AA identity (30 substitutions): expected syncmers perturbed $\approx \min(30 \times 2.9,\, 0.286 \times 291) \approx 83$ of ~83 — i.e. **all of them**. Even as a *set* comparison the Jaccard is low: fraction of unmutated 10-mers $= 0.9^{10} = 0.349 \Rightarrow J = 0.349/(2-0.349) = 0.211$. Murphy-10 helps only for substitutions internal to a group; granting a generous halving of the effective rate gives $0.95^{10} = 0.599 \Rightarrow J = 0.428$. Neither is 1.0, and only 1.0 produces a hash collision.

**Consequence.** Stage 0 collapses only near-identical architectures (≳99.9% AA), not 90%-identical ones. Therefore:
- the "~75% first sighting" Cuckoo-filter hit rate,
- the collapse of ~52 TB input to 2.5×10⁹ contigs / 3.6×10⁹ slots,
- `arch_occurrence_count` as a measure of recurrence,
- and every downstream footprint (1.80 TB archive, 3.45 TB index, the 68-hour Stage-0 wall clock, the three edition sizes)

are **unvalidated**, and all of them move in the direction of *more* data surviving than budgeted. The plan's RAM/disk/time table is not a prediction, it is an assumption.

**Fix.**
1. Replace exact-hash equality with **banded MinHash LSH** at the architecture level: sketch each contig as an ORF-level MinHash set, hash into $b$ bands × $r$ rows, and treat a shared band as a candidate; tune $(b, r)$ to the desired $(s_1, s_2)$ sensitivity curve, with verification by skani (DS-03). This gives a *tunable, stated* similarity threshold instead of an implicit one.
2. State the invariance level honestly in `00-spec.md` Invariant 2 — "near-exact architecture dedup (≳99% AA)" — if exact hashing is kept.
3. **Measure the collapse ratio on a 1% pilot before committing to any storage number.** Replace every absolute footprint in `PLAN.md` §1.3 and the edition sizes in `DEPLOYMENT-AND-MALLEABILITY.md` §4.1 with a formula plus a pilot-measured constant and a stated uncertainty (e.g. "3.45 TB ± 40%, pending pilot"), and make the build plan degrade gracefully (sharded Tier 3 by biome) if the ratio is worse than assumed.

---

## DS-05 — Benchmark 1 is circular

**Severity: CRITICAL.** `PLAN.md` §5 Benchmark 1 requires recovering CBASS, Thoeris, Pycsar, DarTG, DISARM, BREX, Zorya, Septu, Gabija, CreTA, toxSAS, pAgo, MftA-F, PqqA-E, ykkC subtypes, etc. `PLAN.md` §3.2 Stage 2b scans against "PADLOC, DefenseFinder, AcrDB, TADB, TXSScan, MIBiG, PHROGs" HMMs, and §3.1 Step 3b uses "500+ Weinberg/Breaker orphan RNA `.sto` models", and `catalytic_motif_hash` is defined on CDD/SFLD/M-CSA curated active-site columns.

Those HMM libraries **exist because those systems were discovered**. PADLOC and DefenseFinder are curated catalogs of CBASS/Thoeris/Gabija/Zorya; the Breaker `.sto` set contains the `ykkC` subtype alignments the benchmark asks the engine to resolve. "Recall@15 ≥ 95%" therefore measures whether the index can retrieve its own input labels. This is a lookup test presented as a discovery test, and it will pass trivially while saying nothing about the engine's ability to find anything new.

**Fix — dated temporal holdout.** Rebuild the token dictionary and all closure tables from **pre-discovery snapshots**, and require the system to rank partners of discoveries published *after* the freeze. Concretely, define two configurations:

- `dict@2013`: Pfam 27.0, COG 2003/2014, Rfam 11.0, TIGRFAM 15, PDB/Swiss-Prot/PaperBLAST/FEBA frozen 2013-12-31; **no** PADLOC, DefenseFinder, AcrDB, TADB, PHROGs, MIBiG, no Breaker orphan `.sto` set, no `cm_pocket_rules.json` entries for post-2013 subtypes.
- `dict@now`: the full v2.0 library.

Score a panel with per-system publication dates (CBASS 2019; Thoeris 2018/2021; Pycsar 2021; Gabija/Septu/Zorya/Hachiman 2018; DarTG as antiphage 2022; `ykkC`-ppGpp/PRPP 2017–2018; Obelisks 2024; MftA-F mycofactocin 2011+; `ThyX` 2002; `coaX`/`nadV`/`fbiD`/`DpdA-C` with their own dates). Pass criterion: **recall@15 on `dict@2013` for post-2013 systems**, reported with a Wilson CI, with the `dict@now` run shown only as the (uninformative) ceiling and the naive-hypergeometric run as the floor. Any system whose HMM is in the frozen set is excluded from the numerator, not quietly counted.

This is the one benchmark that cannot be gamed by tuning, and it is also the only one that tests the plan's actual thesis.

---

## DS-06 — Zero-leakage single-linkage splitting cannot exist at this scale

**Severity: CRITICAL.** `PLAN.md` §4.5(1) builds `super_cluster_id` by **single-linkage transitive closure** over (≥30% AA ID) ∪ (HHsearch P ≥ 80%) ∪ (shared Pfam clan) ∪ (shared ECOD_H), then links operons sharing ≥2 SuperClusters **or any ECOD_H 3-gram**, then assigns connected components 80/10/10 and claims "0.0% HGT operon leakage and 0.0% remote-homology token leakage."

Two independent reasons this fails:

1. **Single linkage at 30% AA identity chains superfamilies together.** P-loop NTPases, Rossmann folds, HTH, TIM barrels and the ECOD_H level are exactly the levels at which everything is transitively connected. Single linkage has no stopping condition: one 30%-identity bridge merges two clans permanently. The `super_cluster_id` partition will contain a handful of mega-clusters covering most slots.
2. **Percolation on the operon graph.** Nearly every prokaryotic operon contains at least one ubiquitous family (ABC_tran, HTH, response regulator, radical SAM, GNAT — the plan's own §4.1 lists 2,048 hub tokens with >2×10⁶ occurrences each). An edge rule that fires on *any* shared ECOD_H 3-gram gives these hubs degree in the 10⁷–10⁸ range. Mean degree ≫ 1 guarantees a giant component; with hub families present in a large majority of operons, the expected giant component contains **>90% of all operons**. An 80/10/10 split of components where one component is 90% of the data is not a split.

Note also an internal contradiction: `critique_3` §2.1 Step 2 restricts $S(d)$ to **non-housekeeping** SuperClasses, while `PLAN.md` §4.5(1) drops that qualifier and adds the ECOD_H 3-gram rule unconditionally. The stricter version is the one that percolates.

**Fix.**
1. **Measure before claiming.** Make the component-size histogram a *build artifact and a gate*: if the largest component exceeds 20% of operons, the split is rejected and the edge rule is retightened. Publish $|C_{\max}|/|\mathcal{D}|$ in the manifest.
2. **Drop single linkage.** Use average/complete-linkage or profile-level clustering for `super_cluster_id`, and exclude the hub/housekeeping families (and ECOD_H levels above a frequency cut) from the *edge definition* while keeping them as features.
3. **Use a partitioning algorithm with a cut objective, not component assignment** — GraphPart-style balanced homology partitioning (which keeps more data than homology reduction and reports achieved separation), or METIS/Louvain on the operon graph followed by community-to-fold assignment under a bounded cross-partition edge weight.
4. **Replace "0.0%" with a measured bound.** Report the achieved leakage ("≤0.5% of test operons share any SuperCluster with train; 0 share ≥2") with the measurement script. Zero-by-assertion is not auditable; a measured small number is.
5. **Add a leakage-proof second split**: temporal (train on `dict@2013` corpus) and taxonomic-plus-homology (hold out an entire phylum *and* its homology components), reporting both.

---

## DS-07 — The FDR design is misnamed, cannot reach the required tail, and ignores hierarchy

**Severity: CRITICAL.** `PLAN.md` §4.3 (`empirical_fdr_qval`) and `critique_3` §1.2C specify "Target-Decoy Empirical FDR ($q_{\text{emp}}$) by generating decoy pairs `(A, B_shuffled)` permuted strictly within the same (regime, taxonomic_order, neighbor_entropy_bin) stratum", with "within-order, within-regime, entropy-stratified circular block permutation decoys."

**(a) It is not target-decoy.** Target-decoy competition (proteomics) requires a decoy *database* of entities that compete head-to-head with targets for each observation and are known-false by construction. What is described is a **permutation null**. The correct machinery is permutation p-values → Storey π0 → BH, or Efron's empirical null with local fdr. Keeping the name invites applying TDC's error-rate guarantees, which do not hold here.

**(b) Multiplicity exceeds what permutation can resolve.** The plan tests ~3×10⁸ unique family-level pairs, plus a subfamily rollup (potentially 10⁹–10¹⁰), **times** four token levels (clan / family / subfamily / catalytic motif), **times** multiple window sizes $k$ and four topologies. A Bonferroni-equivalent threshold is $\alpha \approx 0.05/3{\times}10^8 = 1.7{\times}10^{-10}$, and 10–20× worse with the secondary dimensions. A permutation p-value is floored at $1/(\#\text{permutations}+1)$: resolving $10^{-10}$ by permutation requires $>10^{10}$ decoys **per pair**. This is not a tuning problem; it is impossible as specified.

**Fix:** a two-tier null. Use the parametric score-test tail (DS-01) for the far tail, **validated** against permutation in the reachable range ($10^{-3}$ to $10^{-6}$) via a QQ plot, with the deviation at $10^{-6}$ reported as the extrapolation error. Add Efron's **empirical null**: fit $\mathcal{N}(\mu_0, \sigma_0^2)$ to the central 90% of the $Z$ histogram, estimate $\pi_0$, and report **local fdr** per pair in addition to tail $q$ — local fdr is what a user ranking the top 100 actually needs.

**(c) The permutation scheme is anti-conservative for exactly the pairs tested.** Circular block permutation preserves within-block adjacency. If block length is comparable to the ±6-ORF window, the decoys retain the local co-occurrence the test is measuring, deflating the null and inflating significance. **Require block length ≥ 4× the window (≥ 25 ORFs) and require the pair to be split across blocks.**

**(d) Add a degree-preserving content null.** Block permutation does not preserve family prevalence and per-genome gene counts simultaneously. Use **curveball/fastball** edge swaps on the genome × family incidence matrix — the standard null in co-occurrence analysis — as a second, independent null, and report results under both. Divergence between nulls is itself a finding.

**(e) Correct across the token hierarchy.** The four levels are nested hypotheses, not 4× independent tests. Use **hierarchical FDR** (Yekutieli's procedure for hierarchically ordered families): test family level first, descend to subfamily only within rejected families, spending α down the tree. Pre-register one **primary** level (recommend `subfamily_50_id`) and label the rest secondary.

**(f) Phylogenetic non-independence in the null.** Co-occurrence across related genomes is not exchangeable. Use lineage-aware testing (as Coinfinder does by discounting clade-driven co-occurrence) and, for presence/absence coupling, a correlated-evolution LRT (Pagel-style dependent vs independent two-state model) rather than Fitch counts (see DS-13).

---

## DS-08 — No power analysis; the statistic is powerless on small replicons

**Severity: MAJOR.** The plan states detection *thresholds* but never a **detection floor**, so users cannot distinguish "no coupling" from "no power."

**Chromosomal regime.** With $N \approx 4{,}000$ ORFs and $k = 6$, $p_0 \approx 2k/N = 0.003$. For a pair with $N_{\text{eff}}^{\text{veh}} = 20$ independent loci, co-occurrence count under the null is Poisson($\lambda = 0.06$):
$$P(X \ge 5) = 6.1\times10^{-9}, \qquad P(X \ge 6) = 6.5\times10^{-11}.$$
At $\alpha = 1.7\times10^{-10}$, **≥6 vehicle-independent co-occurrences** are needed. For $N_{\text{eff}} = 100$ ($\lambda = 0.3$): $P(X\ge8) = 1.3\times10^{-9}$, $P(X\ge9) \approx 4\times10^{-11}$ → **≥9 needed**. So the floor is on the *count*, ~6–9, and is reachable — this should be documented, not discovered by users.

**MGE regime — the statistic dies.** `critique_3` §1.2A itself notes that on a 12-ORF plasmid a ±5-ORF window spans 92% of the replicon. Then $p_0 \to 1$, the null probability of co-occurrence is ~certainty, and $Z_{AB} \to 0$ regardless of biology. **No pair confined to small replicons can ever be significant under this null.** Since Tier 2 (17M phage/plasmid/MGE contigs) is a third of the scientific motivation, this is a hole in the middle of the plan: the regime-conditioning that fixes the false-positive problem creates a total false-negative problem.

**Fix.**
1. Store `min_detectable_count` and a boolean `powered` per pair, computed from $\sum p_m^{(0)}$ at the operative α. Queries must be able to filter on `powered = true`, and the UI must distinguish "tested, no evidence" from "untestable."
2. For small replicons, **change the null rather than accept zero power**: condition on replicon gene content and test *geometry* (adjacency, order, Δbp, strand) against a replicon-matched background pool — i.e. given that both genes are on this replicon, is their arrangement more rigid than in matched replicons carrying both? This moves the MGE evidence from co-presence (uninformative) to micro-synteny (informative), which the plan already computes as $\mathcal{G}_{\text{micro}}$.
3. State explicitly in §4.3 that for low-$N_{\text{eff}}$ pairs the frequentist route is abandoned in favor of structural/complex evidence (pair closure, interface prediction), and rank those separately rather than mixing incomparable scores in one list.

---

## DS-09 — Residual non-independence: within-operon pairs, study pseudo-replication, sampling effort

**Severity: MAJOR.** `STAT-01`'s vehicle dedup addresses one dependence (shared MGE backbone) and declares the problem solved. Three remain:

1. **Pair-within-module dependence.** A single 5-gene operon generates $\binom{5}{2} = 10$ pairs, all driven by one biological fact and one evolutionary event. Ranking pairs independently will fill the top of the list with 10 near-duplicate rows per real discovery, and any FDR computed over pairs treats 10 correlated tests as 10 independent ones. **Fix:** after pair scoring, run a module/community step (threshold the pair graph, extract dense modules) and report **modules** as the discovery unit with pairs as supporting evidence; compute FDR at the module level; use cluster-bootstrap over modules for all CIs.
2. **Study / BioProject pseudo-replication.** `arch_occurrence_count` counts SRA occurrences. SRA sampling is wildly non-uniform — one deeply sequenced time-series project can contribute thousands of near-identical assemblies. Vehicle and phylogenetic dedup do not touch this. **Fix:** cap contribution at **one count per BioProject per vehicle cluster** for inferential purposes (keep the raw count as a separate descriptive column), and store `n_distinct_bioprojects` and `n_distinct_biomes` alongside.
3. **Sampling-effort bias in any prevalence claim.** "Occurs in N samples" conflates abundance with how often that biome was sequenced. **Fix:** report effort-normalized prevalence (occurrences per sequenced Gbp per biome) and never use raw sample counts in a ranking statistic.

---

## DS-10 — Pass criteria without confidence intervals; traps neither pre-registered nor blinded

**Severity: MAJOR.** `PLAN.md` §5 states "≥95% (48/50)", "0.0% false hole calls", "<1% (≤5/500)", "100% (200/200)", "100% (100/100)".

**Arithmetic.** 48/50 has a 95% Wilson interval of **[0.868, 0.985]** — the benchmark cannot distinguish 95% performance from 87% performance. 0 failures in 200 trials bounds the failure rate only at **≤1.5%** (rule of three, 3/200), not 0%. 0 failures in 100 bounds it at **≤3%**. 5/500 has a Wilson interval of roughly **[0.004, 0.023]**. Writing "0%" and "100%" as achievements on these sample sizes is a statistical error that a reviewer will find immediately.

**Fix.**
1. Restate every criterion as a one-sided bound with the method named: "≤1.5% false closure rate (95% upper bound, 0/200 observed)". Never print 0% or 100%.
2. **Enlarge the positive set** from 50 to ≥300 by auto-generating co-functional pairs from MetaCyc/GapMind same-pathway adjacent steps, MIBiG same-BGC pairs, and PDB biological assemblies with real interfaces — then 95% recall has a CI of ±2.5% rather than ±9%.
3. **Pre-register and blind the traps.** The 500-trap panel must be curated by someone other than the person tuning the score, frozen before any scoring run, and its **BLAKE3 hash committed to the repo** with a timestamp. Allow an `UNKNOWN` label: some nominal negatives (integron cassettes, toxin-immunity pairs, prophage morons) are genuinely coupled, and mislabeled negatives put a ceiling on the measurable FPR. Report FPR on confident negatives only, with the `UNKNOWN` count disclosed.
4. State the **ablation table** as a required deliverable with CIs, not a claim: naive hypergeometric and raw host-tree homoplasy must be *run*, not asserted to fail.

---

## DS-11 — The ML evaluation lacks every control that would make it meaningful

**Severity: MAJOR.** `PLAN.md` §5 Benchmark 2A requires `Recall@10 ≥ 65%` on `ECOD_T` fold prediction with "0.0% biophysically incoherent predictions." Four controls are missing, and without them the number is uninterpretable:

1. **Frequency-prior baseline.** Predicting the 10 most common T-groups ignoring all context is a strong baseline when the fold distribution is heavy-tailed. Report it. If `gLM-Domain` does not beat it by a wide margin, the model is memorizing marginals.
2. **Context-shuffled ablation — the single most important missing control.** Re-run inference with the neighbor slots randomly permuted (or replaced by random slots) while keeping the masked slot's position. If Recall@10 does not drop substantially, the model has learned *nothing about synteny* and the whole Index 4 premise is unsupported. This costs one inference pass and is the only direct test of the plan's core scientific claim.
3. **Calibration.** No calibration target appears anywhere in the plan. Require reliability diagrams and **ECE < 0.05** after temperature scaling on validation, for each of the four decoder stages, plus a **risk–coverage curve** so users can choose a confidence cut. Add NLL/Brier alongside Recall@K — a top-K metric alone cannot detect overconfidence, and overconfident fold calls on dark proteins are exactly how this system would mislead a biologist.
4. **Frequency-stratified metrics.** The 131,072-token vocabulary cap means >99% of families collapse to `<RARE_SEQ_TOKEN>`; Stage-4 metrics will be dominated by the common head while the scientific interest is in the tail. Report all metrics by frequency decile, and note that Stage 4 *cannot* name a novel family by construction.

Also: the vocabulary cap is "ranked by phylogenetic dispersion" computed on the full corpus, i.e. including held-out components — a mild but real leakage. Compute selection statistics on the training split only.

---

## DS-12 — Annotation uncertainty is treated as ground truth

**Severity: MAJOR.** Every token in the tape is an *estimate*: `pyrodigal` gene calls, HMM hits at an E-value cut, `catalytic_intact` from alignment columns, and — critically — `struct_fold_id` from **ProstT5-predicted** 3Di for the 18M dark segments that matter most. The plan then computes exact set algebra and exact p-values on these estimates, propagating zero uncertainty. The `catalytic_intact` boolean in particular will be wrong at an unmeasured rate whenever the alignment is marginal, and it gates the closure score that decides what counts as "novel."

Worse, the ProstT5→ECOD path is **out-of-distribution exactly where it is relied upon**: predicted structures for sequence-dark, remote-homology-free proteins are the least reliable, and these are the segments the engine will surface as discoveries.

**Fix.**
1. Store per-slot confidence: HMM E-value/bitscore, `cath-resolve-hits` margin, ProstT5 head entropy, Foldseek E-value and TM-score, and a `fold_source` enum (`transferred_3Di`, `predicted_3Di`, `hmm_inherited`).
2. **Propagate by multiple imputation**: sample $K = 10$ annotation realizations from the top-$k$ per slot, recompute $Z_{AB}$ and the ranking under each, and store the **rank standard deviation** per pair. Report "rank 7 ± 4" instead of "rank 7."
3. **Measure the dark-fold error rate on a prospective set**: evaluate ProstT5+Foldseek ECOD_T assignment against PDB/AFDB entries deposited after the ProstT5 training cutoff and absent from its training data; publish the confusion rate by sequence-identity bin and by segment length. A single number ("ECOD_T top-1 accuracy on novel folds = X%") should appear in the spec, because every downstream claim about dark matter is conditional on it.

---

## DS-13 — Pangenome plane thresholds are undefined at small n; the coherence test is not a test

**Severity: MAJOR.** `PLAN.md` §4.4 sets `core_chr_bits` = present in ≥50% of QC assemblies, `mge_accessory_bits` = present in <30% of strains. For a species with 1, 2 or 3 genomes these fractions are 0, 0.5 or 0.33 with enormous variance, yet they are thresholded as if measured. Most GTDB species clusters are small, so most bits are set by a coin flip.

`critique_1` BIO-08's proposed coherence test — `hole & (hole >> 1)` on a DFS-ordered bitvector, "tested in 1 clock cycle" — is a 1-bit adjacency check that (a) assumes DFS adjacency equals clade membership, which holds only inside a clade and breaks at every subtree boundary, and (b) is not a statistical test: it has no null, no p-value, and no error rate.

**Fix.**
1. Replace hard fractions with **Beta-binomial shrinkage**: estimate species-level presence probability as a posterior mean with a prior fitted at the genus level, store `n_genomes_in_species` as the explicit denominator, and set plane bits from the posterior with a stated credible-interval rule. Species below a power threshold get a third state (`UNKNOWN`), not a 0.
2. Make the bitvectors **three-valued** (present / absent / unknown) — a two-state encoding of a quantity that is often unmeasured is the root cause of the false-hole problem that BIO-09 only partially patches.
3. Replace the shift trick with a real phylogenetic test: model gain/loss with a continuous-time Mk (or Dollo) process on the GTDB tree, and test coupling with a **correlated-evolution likelihood-ratio test** (dependent vs independent two-state models, Pagel-style), reporting a p-value. Precompute per-token transition lists as already planned, but use them as sufficient statistics for the LRT rather than as raw counts.
4. Benchmark 3's "0.0% false hole calls" under 8% simulated dropout is impossible in principle (masking reduces, never eliminates). Restate as: "false-hole rate ≤1% (95% upper bound) **and** nominal $q$ matches observed FDR within 2× on the reliability diagram."

---

## DS-14 — No sensitivity analysis for ~15 hard-coded constants

**Severity: MAJOR.** Jaccard 0.30; $\gamma = 1.4$; EnvCov 0.70 and ResolvedAtomCov 0.70; 50% ID / 80% coverage; ±6 ORFs; 16-slot windows; $H_{\text{ctx}} > 6.5$ bits; IQR scale 40 bp; 95% ANI; 120-aa CTD window; 35-aa nested-domain gap; ≥30% AA leakage cut; 131,072 vocabulary; 0.25/0.10 closure pass lines. Each is a scientific decision presented as a constant. A result that flips when 0.30 becomes 0.35 is not a result.

**Fix.**
1. All constants into `params.yaml`, hashed into the manifest, each tagged with provenance and a plausible range.
2. **One-factor-at-a-time plus global (Morris screening → Sobol indices)** sensitivity analysis on a 1% stratified sample, reporting for each constant the Kendall τ between the top-100 ranking at nominal and at ±50%.
3. Gate: **τ ≥ 0.9** for every constant, or the constant is declared fragile and published as a caveat with the ranking reported under both settings. Fragility is acceptable; undisclosed fragility is not.

---

## DS-15 — No simulation tier, so actual FDR can never be measured

**Severity: MAJOR.** Every benchmark in §5 uses real data with uncertain labels. Nothing in the plan can answer "when the engine reports q = 0.01, what fraction of those are actually false?" — the one question that determines whether the output is usable.

**Fix: add Benchmark 0 (Simulation with Known Truth), to be run first and re-run on every parameter change.** Use **Zombi** (which uniquely models genome evolution in extinct lineages, so transfers from unsampled donors occur — the realistic confounder here) and/or **ALF** (which models rearrangement, LGT, duplication, loss, fusion and fission jointly). Generate ~5,000 genomes with: planted co-evolving modules at known coupling strengths; planted promiscuous MGE passengers with high transfer rates; planted NOGD replacements; realistic rearrangement; then simulate assembly fragmentation and MAG incompleteness.

Deliverables: a **reliability diagram of observed vs nominal FDR**, a power curve (detection probability vs $N_{\text{eff}}$ and coupling strength per regime), and a demonstration that planted passengers are *not* recovered. This is the only artifact in the entire plan that can falsify the statistical layer, and it is cheap relative to a 4.7-day planetary build.

---

## DS-16 — Aggregation destroys the sufficient statistics, so no parameter can be changed without a rebuild

**Severity: MAJOR.** `DEPLOYMENT-AND-MALLEABILITY.md` §3 defines malleability as new rows (epochs), new columns (sidecars), new token IDs (remap), new models (scan). **Changing a statistical parameter is in none of these categories.** Yet `SYS-06` collapses 1.8×10¹⁰ postings into `directon_pair_summary` with scalars like `func_coupling_llr` and `empirical_fdr_qval` baked in. Re-fitting $\lambda_{\text{hub}}$, moving the Jaccard cut, or adopting the score test of DS-01 then requires re-streaming `window_postings.parquet` — which the `lite` and `core` editions **do not even ship**. So the editions most users will run have statistics that are permanently frozen at the vendor's parameter choices, which contradicts the document's own premise.

**Fix.** Persist **per-pair sufficient statistics** rather than only final scores, chosen so that every parameterized score is recomputable in closed form:
- $\sum_m p_m^{(0)}$, $\sum_m p_m^{(0)}(1-p_m^{(0)})$, co-occurrence counts $\sum_m y_m$ — **per regime, per tier, per taxonomic order** (6 × 3 × ~1,800 is too wide; store per regime × tier, plus a 2,048-bit order bitmask as already planned);
- the 32-bin Δbp histogram, topology fractions, and strand/order counts already specified;
- vehicle-cluster count and BioProject count (DS-09);
- the neighbor-distribution sketch needed to recompute $H_{\text{ctx}}$ at a different window.

With these, re-scoring under new parameters is a scan over a ~10 GB table, shippable in `lite`. Add a `bio-box rescore --params params.yaml` subcommand to §2.3's CLI list, and state in §3 that **statistical malleability** is the fifth mechanism. Also record in the manifest: snapshot IDs of every reference DB, RNG seeds, `params.yaml` hash, and code commit — none of which the current manifest carries, so no score is currently reproducible.

---

## DS-17 — No external baseline, no prospective registry

**Severity: MODERATE.** `PRIOR-ART.md` argues at length that nothing comparable exists, and then the plan never compares against anything. "Better than nothing" is unfalsifiable.

**Fix.**
1. **Beat a published baseline on a shared task.** On the subset of genomes each tool can handle: STRING's own combined score, **Coinfinder** (lineage-aware pangenome association/dissociation), a phylogenetic logistic regression, and **gLM2** (the OMG-corpus mixed-modality model, 3.1T bp / 3.3B CDS, which already does masked-slot prediction over genomic context and reports protein-task benchmarks) for Index 4. Report AUPRC on the pre-registered positive/negative sets with cluster-bootstrap CIs. If gLM-Domain does not beat gLM2 on masked-slot fold prediction, the custom model is not justified.
2. **Prospective prediction registry.** Before any publication, commit a hashed, timestamped file of the top 50 novel predicted pairs/modules, with the specific claim for each (physical complex / same pathway / same defense system). Score it 12–24 months later against: newly deposited PDB assemblies, newly published RB-TnSeq cofitness, new DefenseFinder/PADLOC releases, and AlphaFold3/Boltz interface pTM computed by a third party. This converts the system's output from an assertion into a measurable hit rate, and it is the only evidence that does not depend on the authors.

---

## Required document edits (minimum set)

1. `00-spec.md` Invariant 2 and `PLAN.md` §2.3: remove the "95% ANI in 1 CPU cycle" and "90% AA ID invariant" claims; specify 1,024-bit sketch + skani verification + banded-LSH architecture matching; mark all footprint and wall-clock figures as pilot-pending with stated uncertainty (DS-03, DS-04).
2. `PLAN.md` §4.3 Table 2C: replace `func_coupling_llr` as primary metric with `beta_hat`, `beta_se`, `z_score`, `local_fdr`, `q_hier`, `powered`, `min_detectable_count`, `rank_sd`, and the sufficient-statistic columns; keep `func_coupling_llr` as a derived view (DS-01, DS-07, DS-08, DS-12, DS-16).
3. `PLAN.md` §4.3 / new §4.6: add the fitted-weights protocol, `params.yaml`, and the sensitivity-analysis gate (DS-02, DS-14).
4. `PLAN.md` §4.5(1): replace single-linkage components with GraphPart-style balanced partitioning; add the component-size gate; replace "0.0% leakage" with a measured bound; add temporal split (DS-06).
5. `PLAN.md` §4.4: three-valued plane bits with Beta-binomial shrinkage; replace the bit-shift coherence check with a correlated-evolution LRT (DS-13).
6. `PLAN.md` §5: add **Benchmark 0** (Zombi/ALF simulation, FDR reliability diagram, power curves); rewrite **Benchmark 1** as a dated temporal holdout on `dict@2013`; add to **Benchmark 2A** the frequency-prior baseline, context-shuffled ablation, ECE/reliability, and frequency-stratified metrics; restate all criteria with exact-binomial bounds; add trap pre-registration with a committed hash; add the external-baseline comparison and the prospective prediction registry (DS-05, DS-10, DS-11, DS-15, DS-17).
7. `DEPLOYMENT-AND-MALLEABILITY.md` §3: add "Malleability Mechanism 5 — Statistical Rescoring" with `bio-box rescore`, the persisted sufficient statistics, and full provenance (DB snapshot IDs, seeds, params hash, commit) in the manifest (DS-16).

---

## Sources

- [Zombi: a phylogenetic simulator of trees, genomes and sequences that accounts for dead lineages (Davín et al., *Bioinformatics* 2020)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7031779/)
- [ALF — a simulation framework for genome evolution (Dalquen, Anisimova, Gonnet, Dessimoz, *MBE* 29(4):1115–1123, 2012)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3341827)
- [GraphPart: homology partitioning for biological sequence analysis (Teufel et al., *NAR Genomics and Bioinformatics* 5(4), 2023)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10578201/)
- [Fast-Part: fast and accurate data partitioning for biological sequence analysis (bioRxiv 2024)](https://www.biorxiv.org/content/10.1101/2024.11.13.623463.full.pdf)
- [Coinfinder: detecting significant associations and dissociations in pangenomes (Whelan, Rusilowicz, McInerney, *Microbial Genomics* 6, 2020)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7200068/)
- [Fast and robust metagenomic sequence comparison through sparse chaining with skani (Shaw & Yu, *Nature Methods* 20:1661–1665, 2023)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10630134/)
- [The OMG dataset and gLM2: an Open MetaGenomic corpus for mixed-modality genomic language modeling (Cornman et al., ICLR 2025)](https://proceedings.iclr.cc/paper_files/paper/2025/hash/dbb40375daa3c39f1c098b12608aed12-Abstract-Conference.html)
