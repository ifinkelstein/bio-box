# `synthworld`: synthetic-world generator spec

Status: **planning only** (2026-10-29). Nothing is installed or implemented. This file settles `TEST_ARCHITECTURE.md` §11 item 1 ("Zombi fit at micro-corpus size, ALF fallback") and specifies the generator behind §6.1 (micro-corpus), §6.2 (synthetic worlds, null worlds), MR-18, and `PLAN.md` §5 rows B-SIM and B-NULL.

Conventions used throughout come from the store (`PLAN.md` §4.0–4.1, `planning/data_layer_and_tests.md` §2):
- Coordinates are 0-based half-open on the stored contig. Strand is `int8` (+1/−1).
- `contig_uid` is BLAKE3-128 of the uppercase nucleotide sequence, stored as `FIXED_LEN_BYTE_ARRAY(16)`. `protein_uid` is BLAKE3-128 of the amino-acid sequence. `igr_uid` is BLAKE3-128 of the sliced IGR.

Items marked **verify** are not confirmed from a primary source. They are checked in the S−1 fit spike (§1.4).

---

## 1. Verdict

### 1.1 Decision

**Use ZOMBI2 for the evolutionary backbone, and write the rest as custom code.** ZOMBI2 is the 2026 rewrite of Zombi by the same author. It is pinned to an exact version and used at its *ordered* resolution, for the species tree, gene families, gene trees, gene order, strand, HGT, rearrangements and replicons. The custom `synthworld` code realizes nucleotides from that output: protein and codon sequences, IGRs, RNA motifs, CRISPR arrays, phages, fusions, assembly artifacts and truth tables.

**Neither Zombi 1.0.1 nor ALF is used.** The fallback changes from ALF to two options:
- an adapter that reads Zombi 1 `G`-mode gene orders into the same interface;
- if that is also unusable, a small in-house ordered-genome simulator that writes ZOMBI2's file schema.

Either works because `synthworld` consumes only a narrow set of backbone files (§2.2).

### 1.2 Evidence

| criterion | Zombi 1.0.1 (AADavin/Zombi) | ZOMBI2 0.54.0 (AADavin/zombi2) | ALF 0.93 (DessimozLab/ALF) | custom only |
|---|---|---|---|---|
| Python 3 | README asks for "python 3.6" with ETE3, numpy, networkx and pyvolve ([README](https://github.com/AADavin/Zombi)). The S mode broke once on Python 3.8 because pyvolve imported the removed `Bio.Alphabet` ([issue #7](https://github.com/AADavin/Zombi/issues/7)). Pyvolve was archived read-only on 2025-08-07 ([pyvolve](https://github.com/sjspielman/pyvolve)). | `Requires-Python >=3.10`; dependencies are only `numpy>=1.24,<3` and `tqdm` ([PyPI JSON](https://pypi.org/pypi/zombi2/json)). | Not Python. Written in Darwin (`.drw` files in the [repo](https://github.com/DessimozLab/ALF); language "Other (Darwin)" per [GSR](https://surveillance.cancer.gov/genetic-simulation-resources/packages/alf/)). Darwin is closed source ([Wikipedia](https://en.wikipedia.org/wiki/Darwin_(programming_language))). Its binaries were "distributed in object form free of charge" from an archived ETH page last updated 2003 ([biorecipes](https://www.biorecipes.com/darwin/index.html)). | yes |
| gene order with strand | `G` mode writes `<genome>_GENOME.tsv` with POSITION, GENE_FAMILY, ORIENTATION, GENE_ID ([GenomeSimulator.py](https://github.com/AADavin/Zombi/blob/master/GenomeSimulator.py)). | `gene_order.tsv` (lineage, chromosome, topology, position, strand, family, copy) for every node, ancestors included ([output files](https://aadavin.github.io/zombi2/docs/output-files/)). | Rearrangements are modeled: inversion, translocation, inverted translocation ([paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC3341827)). The gene-order output format is not documented in the sources read (**verify**). | must be written |
| genome events | Duplication, transfer, loss, inversion, transposition, origination, with geometric extents and replacement transfer ([GenomeParameters.tsv](https://raw.githubusercontent.com/AADavin/Zombi/master/Parameters/GenomeParameters.tsv)). | The same events, plus translocation, fission, fusion, chromosome origination and loss, and circular or linear topology ([ordered guide](https://aadavin.github.io/zombi2/docs/guide/genomes-ordered/)). **Groups** move as one block, with `together`, `adjacent` and `fates` dials ([CHANGELOG 0.53.0](https://github.com/AADavin/zombi2/blob/main/CHANGELOG.md)). Hub families and `Recipients` decoys are available, and transfers can be distance-biased or clade-restricted ([API](https://aadavin.github.io/zombi2/docs/api/genomes/)). | Duplication, loss, LGT, gene fusion and fission, inversion, translocation ([paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC3341827)). | must be written |
| intergenic regions | `Gf` mode tracks intergenic *lengths* (flat-Dirichlet initial, `_LENGTHS.tsv`) ([Gf wiki](https://github.com/AADavin/Zombi/wiki/Full-genome-(Gf))). In our reading of `SequenceSimulator.py`, `run_s` is a stub and `generate_intergenic_sequences` is never called ([source](https://github.com/AADavin/Zombi/blob/master/SequenceSimulator.py); **verify**). | The nucleotide resolution tracks spacer blocks with breakpoints only between genes, and writes `blocks.tsv`, `genes.tsv`, per-node GFF, BED and FASTA ([nucleotide guide](https://aadavin.github.io/zombi2/docs/guide/genomes-nucleotide/), [output files](https://aadavin.github.io/zombi2/docs/output-files/)). Spacer evolves under `intergene_model` (default JC69) at 3× speed ([sequences guide](https://aadavin.github.io/zombi2/docs/guide/sequences/)). | Sequences are coding sequences plus indels. IGRs are not described in the sources read (**verify**). | must be written |
| sequence models | Pyvolve with nucleotide, amino-acid (WAG, LG, …) or codon (MG, GY) models, and per-family multipliers in `Su` ([SequenceParameters.tsv](https://github.com/AADavin/Zombi/blob/master/Parameters/SequenceParameters.tsv), [source](https://github.com/AADavin/Zombi/blob/master/SequenceSimulator.py)). | JC69/K80/HKY85/GTR and Poisson/Dayhoff/JTT/WAG/LG, Γ+I, per-family multipliers, relaxed clocks, and a `divergence` target ([sequences guide](https://aadavin.github.io/zombi2/docs/guide/sequences/)). **There is no codon model**; the roadmap lists it as "considering" ([roadmap](https://aadavin.github.io/zombi2/docs/roadmap/)). | Codon models (M0–M8, empirical), 7 aa and 4 nt models, indels, Γ, GC amelioration ([paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC3341827)). | custom (§2.3) |
| determinism | Seeds `random` and `numpy.random` only when `SEED≠0`. The default `SEED 0` means unseeded ([GenomeSimulator.py](https://github.com/AADavin/Zombi/blob/master/GenomeSimulator.py), [GenomeParameters.tsv](https://raw.githubusercontent.com/AADavin/Zombi/master/Parameters/GenomeParameters.tsv)). | The same version, arguments and seed give the same dataset, with identical output across operating systems, CPU architectures and Python 3.10–3.13, checked by digest tests on six CI configurations. Seeds are **not** stable across versions. The parallel engine differs from the serial one. A `run.zombi2` provenance report records the version, the arguments and the numpy version ([reproducibility](https://aadavin.github.io/zombi2/docs/reproducibility/)). | Not documented (**verify**). | ours |
| speed | About 1 min for 2,000 taxa with 500-gene root genomes on a laptop ([paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC7031779/)). | Ordered genomes: 0.68 s at 1,000 tips and 6.9 s at 3,000 tips (Apple Silicon). The author reports about 210× over ZOMBI1 on a family task ([performance](https://aadavin.github.io/zombi2/docs/performance/), [repo](https://github.com/AADavin/zombi2)). The genome size in that benchmark is not stated (**verify** at 2,000–4,000 genes). | 4.5 h for 20 *E. coli* genomes (4,352 genes) under an empirical codon model on one 2.26 GHz core ([paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC3341827)). That is over budget for a 15-min micro-corpus. | fast if numpy |
| maintenance | Last commit 2025-04-03, a reverted optimization ([commits API](https://api.github.com/repos/AADavin/Zombi/commits?per_page=5)). One open issue (#6, 2020: `MIN_LINEAGES` not implemented) ([issues](https://github.com/AADavin/Zombi/issues)). The author now develops ZOMBI2. | Very active: 1,815 commits; releases from 0.2.0 (2026-07-07) to 0.54.0 (2026-10-25); CI and a `tests/` directory. Development Status 4 (Beta); single maintainer; 4 stars ([PyPI JSON](https://pypi.org/pypi/zombi2/json), [repo](https://github.com/AADavin/zombi2), [CHANGELOG](https://github.com/AADavin/zombi2/blob/main/CHANGELOG.md)). | Repo created and last pushed 2020-04-27 ([API](https://api.github.com/repos/DessimozLab/ALF)). Open issue #2 (2021, `SE_CodonMutate` error) ([issues](https://github.com/DessimozLab/ALF/issues)). The last release was about 14 years ago ([GSR](https://surveillance.cancer.gov/genetic-simulation-resources/packages/alf/)). | ours |
| license | GPL-3.0-or-later ([README](https://github.com/AADavin/Zombi)). | MIT ([docs](https://aadavin.github.io/zombi2/), [PyPI JSON](https://pypi.org/pypi/zombi2/json)). | No license in the repo (`license: null`, [API](https://api.github.com/repos/DessimozLab/ALF)). The Darwin runtime is closed. | ours |

### 1.3 Deciding reasons

1. **Python 3 and dependencies.** ZOMBI2 needs only numpy and runs on Python 3.10–3.13. Zombi 1's sequence path depends on pyvolve, which is now archived. ALF needs a closed Darwin binary whose official distribution page is an unmaintained archive. That binary cannot be pinned or rebuilt for the `bio-box-test` image.
2. **Gene-order output fits the store.** ZOMBI2 writes signed gene order at every node, with chromosome topology. It also has *groups*, families that move as one block with tunable co-transfer, adjacency and fates. Groups are exactly the planted-module mechanism B-SIM needs, and they are absent from Zombi 1 and ALF.
3. **Determinism is documented and tested** across platforms, with a provenance report. Zombi 1's default is unseeded. ALF's seeding is undocumented.
4. **Speed.** At the ordered resolution the micro-corpus backbone takes seconds. ALF's codon simulation would by itself exceed the 15-min budget.
5. **License.** MIT fits a test image without copyleft questions. ALF has no license, so redistributing it in an image is not clearly permitted.
6. **What ZOMBI2 lacks, we write ourselves:**
   - a codon model and HMM-profile-aware protein evolution, so planted families are recognizable by their source Pfam HMM at a set divergence;
   - IGR structure with promoters, RBS and RNA motifs;
   - CRISPR and protospacer coupling, fusions and splits, and assembly and caller artifacts.

   None of these exist in ALF or Zombi 1 either, in the form the truth tables need.

**Main risk:** ZOMBI2 is beta, has a single maintainer, and has released about 50 versions in 4 months. Mitigations:
- pin `zombi2==0.54.0` and vendor the wheel into the test image;
- consume only the five backbone files in §2.2, behind an adapter with its own schema tests;
- record ZOMBI2's version and `run.zombi2` in the world manifest (§6).

An upgrade changes worlds; it is a `synthworld` MAJOR bump (§6.3).

### 1.4 S−1 fit spike (replaces "check Zombi fit")

Run once, on a laptop, before sign-off. Results go into this file.

1. Install `zombi2==0.54.0` in a clean venv. Run `simulate_genomes_ordered` on a 46-tip tree with `initial_families ≈ 2,000` and the micro preset rates (§5.1). Record wall time and peak memory. Pass: under 60 s.
2. Confirm that `gene_order.tsv` `copy` maps one-to-one to gene-tree leaves `n<species>_g<copy>` ([output files](https://aadavin.github.io/zombi2/docs/output-files/); **verify**).
3. Confirm group semantics:
   - members stay contiguous and co-oriented in every extant genome after inversions;
   - `together<1` transfers subsets;
   - `fates` per-family losses work;
   - `replacement` is refused with groups, as documented ([API](https://aadavin.github.io/zombi2/docs/api/genomes/)).
4. Confirm that the same seed gives byte-identical outputs on macOS and Linux CI. Confirm that `simulate_genomes_ordered` has no `parallel` argument, so it always runs the serial engine (the [API](https://aadavin.github.io/zombi2/docs/api/genomes/) excerpt read lists none; **verify**).
5. Confirm that genome size stays within a ±30% band of target across the tree when origination and loss rates are balanced.
6. If any of 2–4 fail and cannot be fixed by an upstream issue within the S−1 window, switch to the fallback in §1.1.

### 1.5 Consequential edits to other files (not made here)

- `TEST_ARCHITECTURE.md` §0.3, §6.1, §6.2, §8 and §11.1: replace "Zombi with Pyvolve; ALF fallback" with "ZOMBI2 (pinned) + `synthworld` realization". Replace "sequences by `hmmemit` at set divergence" with the profile-constrained evolver in §2.3 (`hmmemit` stays as a cross-check).
- `DEPLOYMENT-AND-MALLEABILITY.md` §2.1 image table: the `bio-box-test` layer carries `zombi2==0.54.0`, `blake3` and `synthworld`. It drops Zombi 1, Pyvolve and ALF.
- `PLAN.md` §5 B-SIM: "Zombi and/or ALF" becomes "`synthworld` (ZOMBI2 backbone)".
- `TEST_ARCHITECTURE.md` §8 `test_run`: add the columns in §6.4.

### 1.6 Other tools considered

| tool | what it does | why not the backbone |
|---|---|---|
| SimPhy | Gene trees inside locus trees inside a species tree, with ILS, DTL and gene conversion. Sequences come via INDELible ([paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4748750/), [repo](https://github.com/adamallo/SimPhy), GPL-2.0). | No gene order or genome structure. |
| SaGePhy | Species, gene and domain trees, with sub-gene (domain) duplication, transfer and loss, and additive or replacing transfer. Java, BSD ([site](https://compbio.engr.uconn.edu/software/SaGePhy)). | No gene order. Its domain-level events could inform fusion planting later. |
| EvolSimulator | Gene-family evolution with LGT under a birth–death model; amino-acid only; no user root input ([GSR](https://surveillance.cancer.gov/genetic-simulation-resources/packages/evolsimulator), [NU-IN note](https://pmc.ncbi.nlm.nih.gov/articles/PMC3161368/)). | Old; protein-only. |
| sgEvolver (Mauve) | Nucleotide genomes: HKY substitutions via seq-gen, indels, gain and loss from a donor pool, inversions ([Mauve dev guide](https://asap.genetics.wisc.edu/software/mauve/mauve-developer-guide/evaluating-alignment-quality-and-stress-testing-the-aligner.php), [arXiv 1207.6964](https://arxiv.org/pdf/1207.6964)). | No gene families or gene trees. Duplications are not tracked. |
| Sim-it | Structural-variant and long-read simulator (**verify**; not fetched). | It is a read simulator, not a genome-evolution simulator. |
| INDELible, pyvolve, AliSim alone | Sequence evolution along a given tree. INDELible and AliSim are **verify**; pyvolve is archived ([pyvolve](https://github.com/sjspielman/pyvolve)). | No genome model. Our evolver needs per-site HMM profiles and codon-usage coupling, which none of them offer as-is. |
| custom on dendropy/ete | Everything ours. | DTL with dead lineages, segmental events and groups are what ZOMBI2 already does well and fast. Rewriting them is the fallback, not the plan. |

---

## 2. Architecture

### 2.1 Stages and tools

`synthworld build --preset <name> --seed <u64> --out worlds/<world_id>/` runs these stages in order. Each stage reads only the outputs of earlier stages and its own child RNG stream (§6.1).

| # | stage | tool | output |
|---|---|---|---|
| S1 | **Species tree.** Birth–death with extinct lineages kept for transfers ([ZOMBI2 docs](https://aadavin.github.io/zombi2/)). Optional clonal bursts: short-branch clusters grafted onto chosen tips, to give uneven sampling. Synthetic taxonomy from depth cuts (phylum … species). | ZOMBI2 `species` + custom grafting | `species_complete.nwk`, extant tree, `truth_genome.tax_lineage` |
| S2 | **Genome evolution with gene order.** Families with per-family rates (core vs accessory). Segmental duplication, transfer and loss with geometric extents. Inversion, transposition, translocation. Chromosome plus 0–2 plasmid replicons. Groups for planted modules, vehicles (chassis + cargo) and Cas operons. Hub families. `Recipients` decoys. Distance-decay transfer ([ordered guide](https://aadavin.github.io/zombi2/docs/guide/genomes-ordered/), [API](https://aadavin.github.io/zombi2/docs/api/genomes/)). | ZOMBI2 `simulate_genomes_ordered` (serial engine) | `gene_order.tsv`, gene trees, `genome_events.tsv`, `groups.tsv`, `chromosome_events.tsv` |
| S3 | **Family design.** Each family gets a domain architecture of 1–3 Pfam HMMs from the world's HMM panel, or a *dark* profile from a Pfam HMM held out of the scan set. Plus N/C tails and linkers. Lengths are fitted so the protein-length distribution matches the golden genomes (§7). Promiscuous-domain traps (HTH, P-loop) are prepended to many unrelated families. | custom; pyhmmer to parse HMMs (`match_emissions`, `transition_probabilities`; [pyhmmer docs](https://pyhmmer.readthedocs.io/en/v0.8.1/api/plan7.html)) | family table |
| S4 | **Family sequences.** Profile-constrained codon evolution along each family's *complete* gene tree, including dead lineages and transfer branches. Lineage codon usage, GC and genetic code. Details in §2.3. | custom numpy | per-copy CDS |
| S5 | **IGRs and regulatory elements.** Each gene copy carries an upstream cassette: optional σ70-like promoter PWM, a 5′ UTR, optional RNA motif, and an RBS 5–10 nt before the start. Operon-last genes may also carry a downstream terminator stem-loop. Cassettes evolve along the gene tree with the gene. The IGR is assembled per junction type (§2.4). RNA motifs are sampled with `cmemit` from Rfam CMs, once per gene-tree clade (§2.5). | custom + Infernal `cmemit` ([man](https://manpages.debian.org/trixie/infernal/cmemit.1.en.html)) | replicon sequences, `truth_igr`, `truth_rna` |
| S6 | **Phages and CRISPR.** A separate small ZOMBI2 ordered run with linear topology gives a related phage set, realized by S3–S5. Some phages are integrated as prophages in host clades. Arrays are seeded at species-tree nodes next to a Cas group. Spacers are acquired at the leader end along branches and deleted internally. Spacers are slices of phage genomes next to a PAM, with mismatches. | ZOMBI2 + custom | phage contigs, arrays, `truth_spacer_target` |
| S7 | **Fusions, splits, NOGD.** Fuse adjacent co-directional pairs, or split multi-domain ORFs, in every extant genome of a chosen clade. Non-orthologous displacement relabels one module member to a new family in a clade. | custom | `truth_fusion`, `truth_module` (NOGD rows) |
| S8 | **Gene-call error injection.** Frameshift indels, upstream and downstream alternative starts, planted smORFs (10–50 aa), internal-stop pseudogenes, and clades using genetic code 4 or 25 (assigned in S4). | custom | `truth_callerr` |
| S9 | **Assembly and fragmentation.** Per genome, an assembly class: complete, draft (target N50) or MAG (target N50 + completeness 70–95%). Breakpoints prefer repeats (vehicle chassis copies, synthetic rRNA operons, recent duplicates ≥ 98% identical, ≥ 1 kb), with random breaks added to reach the N50. Contigs under 500 bp are dropped. Each contig is reverse-complemented with p = 0.5. Complete circular replicons are rotated to a random origin, and at least k origin-spanning genes are planted (MR-04). | custom | contigs, `truth_contig`, `truth_edges` |
| S10 | **Contamination and duplicates.** Foreign contigs from a distant clade or an outgroup world with a different GC are added to chosen assemblies. Exact duplicates are written as separate "sources" under new accessions and BioProjects. Near-duplicates at 99.9% ANI (a mutation fraction drawn per genome, biased to synonymous and IGR sites). A mock FCS-GX sidecar flags the contaminant contigs (MR-19). | custom | `truth_contam`, `truth_dup`, `sidecars/fcs_gx_mock.parquet` |
| S11 | **Emit.** Per-source FASTA and an ingest manifest (accession, BioProject, `source_class`, `replicon_class`, `molecule_type`, synthetic `tax_lineage` with an optional mislabel rate, biome). Truth Parquet tables. A convenience truth GFF. `manifest.json` with digests (§6.2). | custom | world directory |

### 2.2 Backbone interface (the only ZOMBI2 coupling)

`synthworld.backbone` reads five ZOMBI2 outputs into typed tables: species tree, `gene_order.tsv`, gene trees (`gene_trees.tsv` with `--compact`), `genome_events.tsv` and `groups.tsv` ([output files](https://aadavin.github.io/zombi2/docs/output-files/)).

An R1 schema test pins the columns. Any adapter (Zombi 1 `_GENOME.tsv`, or the in-house fallback) must produce the same tables. The raw ZOMBI2 directory, including `run.zombi2`, is kept under `worlds/<id>/zombi2/` for audit.

### 2.3 Family sequence model (replaces "hmmemit at set divergence")

**Why `hmmemit` alone is not enough.** It samples each sequence independently from the profile ([man](https://www.mankier.com/1/hmmemit)). Its options include `-N`, `-a`, `-c`, `-C`, `-p`, `-L` and `--seed`. There is no divergence, temperature or identity option ([Debian man, HMMER 3.4](https://manpages.debian.org/trixie/hmmer/hmmemit.1.en.html)). Independent samples form a star phylogeny at the profile's saturation identity, so they carry neither tree structure nor controllable divergence. Divergence therefore comes from evolving along the gene tree.

**Profile-F81 codon evolver** (custom numpy):
- **Sites.** A family's protein is a concatenation of segments: N-tail, domain₁, linker, …, domainₖ, C-tail.
  - A domain site *s* has the equilibrium distribution πₛ, its source HMM's match-state emission row (pyhmmer `match_emissions`; [docs](https://pyhmmer.readthedocs.io/en/v0.8.1/api/plan7.html)).
  - Insert-state residues sampled at the root, and all tail and linker sites, use the lineage's background amino-acid composition.
- **Root.** The root is drawn from πₛ per match site, with insert lengths from the HMM transitions. This is a pure-numpy equivalent of `hmmemit` from the core model.
- **Substitution.** On a branch of length *t*, site *s* has a substitution event with probability 1 − exp(−μ_f · r_s · c_b · t). Here μ_f is the family rate, r_s a Γ rate category, and c_b a lineage relaxed-clock multiplier. On an event, the new amino acid is drawn from πₛ (F81 with site-specific frequencies). This is exact, needs no matrix exponentials, and vectorizes over sites and branches.
- **Divergence control in closed form.** For two copies separated by path length T, expected identity is mean over *s* of [e^{−μ r_s T} + (1 − e^{−μ r_s T})·hₛ], with hₛ = Σₐ πₛₐ². `synthworld` solves for μ_f to hit a target mean pairwise identity, or a target bitscore band (§5.3). At T → ∞ the members become independent `hmmemit`-like draws.
- **Codons.** A changed amino acid gets a codon drawn from the lineage's codon usage, conditioned on the amino acid and the clade's genetic code. Synonymous events occur at rate κ_syn·μ_f and redraw from the same conditional. GC and codon usage follow a Brownian process on the species tree, so transferred genes carry donor usage and ameliorate (as ALF models GC amelioration; [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC3341827)).
- **Indels.** v1 has indels only in tails and linkers (1–10 aa), so domain envelopes stay full-length and hit tests stay interpretable. Domain-internal indels are a v2 knob.
- **Each gene-tree node's species** is read from internal labels `<event>_n<species>` ([output files](https://aadavin.github.io/zombi2/docs/output-files/)). That species gives the node's codon usage and genetic code.

**Cross-check (R1).** For 20 HMMs, the score distribution of `hmmsearch` on (a) our t → ∞ samples and (b) `hmmemit -N 500 --seed <n>` samples must agree (two-sample KS p > 0.01 at fixed seeds). This keeps the numpy sampler honest against the binary.

**Circularity guard.** Planted families scanned with the same HMMs that generated them give optimistic recall. Two knobs:
- `root_source = real_seed_member`: start from a real Pfam seed-alignment sequence instead of an HMM sample;
- `dark_fraction`: families generated from Pfam HMMs excluded from the scan panel.

### 2.4 IGR assembly per junction

For adjacent genes i, i+1 on a replicon, the junction type comes from the strands: co-directional, divergent (←→) or convergent (→←). For co-directional junctions it also comes from the operon class:
- **tight** for every within-group pair;
- **tight** with probability p_tight for background pairs, decided by a hash of (family_i, family_{i+1}, orientation). The same adjacency in different lineages therefore gets the same class, which mimics inheritance without a planted coupling.

| junction | sequence | length source |
|---|---|---|
| co-directional, tight | RBS part of U_{i+1}; overlaps −1 (`TAATG`) and −4 (`ATGA`) realized by templates | golden intra-operon distribution |
| co-directional, loose | U_{i+1} with its promoter | golden inter-directon co-directional distribution |
| divergent | revcomp(U_i) + U_{i+1} | golden divergent distribution |
| convergent | D_i + revcomp(D_{i+1}) | golden convergent distribution |

Each copy's cassette has a *preferred length*. It is inherited along the gene tree with small multiplicative noise, and the realized IGR is trimmed or padded from the cassette to the target. Orthologous IGRs in close relatives are therefore near-identical, which matters for ANI, dedup and motif conservation.

Cassette nucleotides evolve at 3× the gene's synonymous rate. Elements (promoter boxes, RBS, motif paired positions) evolve at a lower rate. This is the same idea as ZOMBI2's `intergene_speed` default of 3 ([sequences guide](https://aadavin.github.io/zombi2/docs/guide/sequences/)).

### 2.5 RNA motifs

For each (Rfam CM, attached family) link:
1. Cut the family's gene tree at height h into clades.
2. Run `cmemit -N 1 -a --dna --seed <derived> [--exp x] model.cm` once per clade, keeping the Stockholm consensus structure ([cmemit man, Infernal 1.1.5](https://manpages.debian.org/trixie/infernal/cmemit.1.en.html)).
3. Within a clade, evolve the motif along the gene tree. Substitutions hit unpaired positions; paired positions change only by compensatory swaps among WC/GU pairs.

Difficulty knobs:
- `--exp` < 1 gives sequences harder for the CM to distinguish ([man](https://manpages.debian.org/trixie/infernal/cmemit.1.en.html));
- `-l` gives local-mode emission;
- h sets clade depth.

`cmemit --seed 0` means an arbitrary seed, so `synthworld` always passes a derived non-zero seed ([man](https://manpages.debian.org/trixie/infernal/cmemit.1.en.html)). The same rule applies to `hmmemit` ([man](https://www.mankier.com/1/hmmemit)).

### 2.6 Runtime budget (micro preset)

| part | estimate |
|---|---|
| ZOMBI2 ordered run | seconds (**verify**, §1.4) |
| S4 evolver: about 3,000 families × ~300 sites × ~200 tree nodes, vectorized | under 2 min single core (estimate, measured at S−1) |
| `cmemit` | about 50 calls |
| all of `synthworld` | target ≤ 3 min, leaving ≥ 12 min of the 15-min micro-corpus budget for the bio-box build (`TEST_ARCHITECTURE.md` §6.1) |

---

## 3. Truth tables

### 3.1 General rules

- **Format.** All tables are Parquet with schemas registered in the schema registry (`tests/registry/`, LY-01). Type names: `uid16` = `FIXED_LEN_BYTE_ARRAY(16)`. Coordinates are `uint32`, 0-based half-open on the **stored contig**. Strand is `int8` (+1/−1).
- **Common columns.** Every table also carries `world_id string` and `release_tag string` (the world's tag, not a store `release_id`).
- **Join keys.** Content keys are `contig_uid`, `protein_uid` and `igr_uid`, computed exactly as the store computes them. Exact duplicates share `contig_uid`. The per-occurrence key is therefore (`source_label`, `accession`), which joins to the store's `contig_occurrence`.
- **UID cross-check.** `synthworld` computes UIDs with the `blake3` package. An R1 test asserts equality with the store's UID function on every contig and protein, so a store UID bug cannot hide behind a shared implementation.
- **Slot join rule.** The store's `slot_ordinal` depends on `gene_call_version`, so truth never stores it. Truth genes join to called slots by (`contig_uid`, `strand`, `stop_bp`): the stop-anchored match standard in prokaryotic gene-call evaluation. Start agreement is then scored separately. `stop_bp` is the coordinate just past the stop codon on the + strand, and the stop codon's first base on the − strand. Both are expressed in contig coordinates.
- **Protein UIDs.** `protein_uid` hashes the true translation without the terminal `*`. It must match the store's convention (**verify** against `data_layer_and_tests.md` §2A).

### 3.2 Supporting tables (needed for every join)

**`truth_genome`** (one row per genome per source copy)

| column | type | notes |
|---|---|---|
| genome_id | string | world-unique, e.g. `g017` |
| species_node | string | ZOMBI2 tip label |
| source_label | string | `src_a`, `src_b`, …; duplicates get a new label |
| accession | string | synthetic assembly accession |
| bioproject | string | synthetic; bursts share one |
| tax_lineage | string | synthetic, depth-cut; possibly mislabeled |
| tax_lineage_true | string | before mislabeling |
| genetic_code | uint8 | 11, 4 or 25 |
| assembly_class | string | `complete`, `draft`, `mag` |
| completeness_true | float32 | 1.0 unless MAG |
| gc_true | float32 | |
| is_phage | bool | |

**`truth_contig`** (one row per emitted contig occurrence)

| column | type | notes |
|---|---|---|
| contig_uid | uid16 | |
| accession | string | contig accession |
| source_label | string | |
| genome_id | string | |
| length | uint32 | |
| replicon_id | string | `chr`, `p1`, `phage_07`, … |
| replicon_class | string | chromosome, plasmid, phage, prophage-bearing |
| is_circular | bool | true only for unfragmented circular replicons |
| is_contaminant | bool | |
| is_dropped | bool | MAG dropout; the sequence is kept in `dropped/` and not emitted |

**`truth_gene`** (one row per true gene on an emitted contig)

| column | type | notes |
|---|---|---|
| gene_id | string | world-unique |
| contig_uid | uid16 | |
| start | uint32 | |
| end | uint32 | |
| strand | int8 | |
| stop_bp | uint32 | slot join key (§3.1) |
| partial_5 | bool | truncated at a contig edge |
| partial_3 | bool | truncated at a contig edge |
| origin_spanning | bool | |
| feature_kind | string | SO-mapped: CDS, smORF, pseudogene, rRNA |
| family_id | string | |
| copy_label | string | ZOMBI2 leaf `n<sp>_g<copy>` |
| protein_uid | uid16 | null for non-CDS |
| domain_arch | list<struct<hmm_acc string, aa_start uint16, aa_end uint16>> | aa coordinates 0-based half-open |
| module_id | string | null if none |
| vehicle_id | string | null if none |
| operon_id | string | truth transcription unit |
| genetic_code | uint8 | |
| hgt_on_lineage | uint16 | transfer events on the root-to-leaf path |

**`truth_igr`** (one row per junction)

| column | type | notes |
|---|---|---|
| contig_uid | uid16 | |
| start | uint32 | negative-length overlaps are encoded with `start = end` plus `overlap_bp` |
| end | uint32 | |
| overlap_bp | uint8 | |
| igr_uid | uid16 | null when overlap |
| left_gene_id | string | |
| right_gene_id | string | |
| junction | string | `codir_tight`, `codir_loose`, `divergent`, `convergent` |
| elements | list<struct<kind string, start uint32, end uint32, strand int8>> | promoter_m35, promoter_m10, rbs, terminator, rna_motif |

### 3.3 The ten required tables

**`truth_family`** (one row per family; members via `truth_gene.family_id`)

| column | type | notes |
|---|---|---|
| family_id | string | |
| zombi_family | uint32 | ZOMBI2 family id |
| kind | string | `pfam`, `dark`, `phage`, `cas`, `vehicle_chassis`, `cargo`, `hub`, `decoy` |
| domain_hmms | list<string> | Pfam accessions with version |
| in_scan_panel | bool | false for dark families |
| root_source | string | `hmm_sample`, `real_seed_member` |
| rate_multiplier | float32 | μ_f |
| target_identity | float32 | |
| realized_mean_identity | float32 | over extant pairs |
| realized_min_identity | float32 | |
| ref_bitscore_median | float32 | source-HMM score of members, computed at build time by pyhmmer (difficulty label, not truth of detection) |
| n_extant_copies | uint32 | |
| n_genomes_present | uint32 | |
| n_transfers | uint32 | from `genome_events.tsv` |
| gene_tree_extant | string | Newick with `gene_id` leaves |

**`truth_module`** (one row per module member)

| column | type | notes |
|---|---|---|
| module_id | string | |
| kind | string | `planted`, `vehicle_chassis`, `cas_operon`, `profile_only` (`adjacent=False`) |
| family_id | string | |
| member_role | string | `core`, `nogd_replacement`, `nogd_replaced` |
| nogd_clade_node | string | null unless NOGD |
| order_index | uint8 | position in the group |
| together | float32 | ZOMBI2 co-transfer probability |
| fates | string | ZOMBI2 group fates |
| adjacent | bool | |
| origin_node | string | species-tree node of group origination |
| strength_class | string | `weak`, `medium`, `strong` (§5.3) |
| realized_loci | uint32 | extant loci carrying ≥ 2 members adjacent |
| realized_independent_transfers | uint32 | group transfer events with an extant descendant |
| realized_vehicle_independent_loci | uint32 | loci in distinct flank contexts (same ±15-ORF rule as `vehicle_cluster_id`) |
| powered_expected | bool | realized loci ≥ the `min_detectable_count` rule (`PLAN.md` §4.3), so recall is scored only on modules the data can show |

A view `truth_module_pairs(module_id, family_a, family_b, kind)` enumerates within-module pairs for pair-level scoring. Module-level scoring follows `PLAN.md` §4.3 (modules are the discovery unit).

**`truth_fusion`** (one row per extant occurrence)

| column | type | notes |
|---|---|---|
| fusion_id | string | event id |
| kind | string | `fusion`, `split` |
| family_a | string | |
| family_b | string | for a split, the two halves' parent family appears in both columns, with `split_domain_index` |
| split_domain_index | uint8 | |
| clade_node | string | |
| genome_id | string | |
| contig_uid | uid16 | |
| start | uint32 | |
| end | uint32 | |
| strand | int8 | of the fused ORF, or of the split pair |
| junction_aa | uint16 | protein coordinate where A ends (fusion) |
| linker_aa | uint8 | |
| gene_ids | list<string> | 1 for fusion, 2 for split |
| expected_arch | string | `fused`, `split_adjacent` (MR-11) |
| in_module | string | module_id or null |

**`truth_rna`** (one row per planted motif instance)

| column | type | notes |
|---|---|---|
| rna_id | string | |
| rfam_acc | string | with Rfam release |
| cm_name | string | |
| clade_emission_id | string | shared within a gene-tree clade |
| exp | float32 | `cmemit --exp` |
| local_mode | bool | |
| attached_gene_id | string | |
| attached_family_id | string | |
| contig_uid | uid16 | |
| start | uint32 | |
| end | uint32 | |
| strand | int8 | |
| dist_to_start | int32 | bp from motif end to start codon |
| n_subst_from_emission | uint16 | |
| n_compensatory | uint16 | |
| ref_bitscore | float32 | cmsearch score of the motif alone at build (difficulty label) |

**`truth_spacer_target`** (one row per spacer × target)

| column | type | notes |
|---|---|---|
| array_id | string | |
| genome_id | string | |
| contig_uid | uid16 | |
| array_start | uint32 | |
| array_end | uint32 | |
| array_strand | int8 | leader side = 5′ |
| repeat_consensus | string | |
| cas_module_id | string | |
| spacer_index | uint16 | 0 = leader-proximal |
| spacer_start | uint32 | |
| spacer_end | uint32 | |
| spacer_seq | string | |
| acquisition_node | string | species-tree branch |
| target_contig_uid | uid16 | |
| target_kind | string | `phage_contig`, `prophage`, `plasmid`, `none` (decoy spacer) |
| proto_start | uint32 | |
| proto_end | uint32 | |
| proto_strand | int8 | |
| n_mismatches | uint8 | |
| pam_seq | string | |
| pam_present | bool | PAM rule per Cas type (**verify** PAM conventions before planting) |

**`truth_dup`** (one row per member contig)

| column | type | notes |
|---|---|---|
| dup_group_id | string | |
| kind | string | `exact`, `near` |
| ani_target | float32 | |
| ani_realized | float32 | |
| rep_genome_id | string | |
| member_genome_id | string | |
| member_source_label | string | |
| member_accession | string | |
| rep_contig_uid | uid16 | |
| member_contig_uid | uid16 | equal to `rep_contig_uid` for exact duplicates (MR-06) |
| n_mutations | uint32 | |
| expected_collapse | bool | for counting (MR-06, MR-07) |

**`truth_contam`** (one row per contaminant contig)

| column | type | notes |
|---|---|---|
| contig_uid | uid16 | |
| host_genome_id | string | |
| donor_genome_id | string | |
| donor_world | string | `main` or `outgroup` |
| tree_distance | float32 | host to donor, null for outgroup |
| length | uint32 | |
| gc | float32 | |
| kind | string | `cross_clade`, `outgroup`, `vector_like` |
| in_mock_mask | bool | flagged in `fcs_gx_mock` (MR-19) |

**`truth_edges`** (one row per emitted contig)

| column | type | notes |
|---|---|---|
| contig_uid | uid16 | |
| genome_id | string | |
| replicon_id | string | |
| replicon_start | uint32 | in the *pre-rotation, forward* replicon frame |
| replicon_end | uint32 | in the same frame |
| orientation | int8 | −1 if the contig is the reverse complement |
| rotation_offset | uint32 | for circular complete replicons |
| left_edge | string | `replicon_end`, `break_repeat`, `break_random`, `none` |
| right_edge | string | same values |
| repeat_family_left | string | null if not a repeat break |
| repeat_family_right | string | null if not a repeat break |
| n_partial_genes | uint8 | |

A helper view maps any replicon coordinate to (`contig_uid`, contig coordinate, strand), so MR-05 and S-09 edge tests can compute expected censoring.

**`truth_callerr`** (one row per injected error)

| column | type | notes |
|---|---|---|
| callerr_id | string | |
| kind | string | `frameshift_ins`, `frameshift_del`, `alt_start_upstream`, `alt_start_downstream`, `smorf`, `internal_stop`, `alt_code` |
| genome_id | string | |
| contig_uid | uid16 | |
| true_start | uint32 | the intended gene model |
| true_end | uint32 | |
| strand | int8 | |
| error_pos | uint32 | |
| true_gene_id | string | |
| true_protein_uid | uid16 | null for pseudogenes |
| expected_calls | list<struct<start uint32, end uint32, strand int8>> | what a correct standard-code caller should emit |
| expected_arch | string | e.g. `not_split` for a frameshift: O-ARCH must not call a split (F-08) |
| genetic_code | uint8 | |

**`truth_traps`** (one row per trap × family)

| column | type | notes |
|---|---|---|
| trap_id | string | |
| kind | string | `promiscuous_domain`, `mge_cargo`, `hub_family`, `decoy_profile`, `clonal_burst`, `tandem_paralog` |
| family_id | string | |
| vehicle_id | string | for cargo |
| partner_family_ids | list<string> | families it co-occurs with by construction |
| expected_label | string | `must_not_call`, `UNKNOWN` |
| rationale | string | e.g. "cargo always shares chassis flanks; vehicle dedup must collapse" |

Trap pairs are disjoint from `truth_module_pairs` by construction, except vehicle chassis, which are `truth_module.kind = vehicle_chassis` and excluded from both recall and null evaluation.

---

## 4. Null worlds

### 4.1 Constructions

**N0, generative null (primary; B-NULL).** N0 uses the same preset, seed and stage code as the signal world, with **zero planted groups**:
- `truth_module` has no `planted`, `profile_only` or `cas_operon` rows;
- no `Recipients` decoys and no NOGD.

N0 keeps unchanged:
- the species tree, including dead lineages and clonal bursts;
- all DTL, rearrangement, extent and replicon rates, and hub families;
- vehicles (chassis + cargo as a ZOMBI2 group);
- phages, prophages and CRISPR arrays;
- promiscuous domains, background tight junctions (hash rule, §2.4);
- fragmentation, call errors, contamination and duplicates.

Two variants:
- **N0-full** (default): as above. Null evaluation covers every pair not in `truth_module`; chassis pairs are excluded.
- **N0-strict**: also has no vehicles and no Cas group, so every pair is null. It gives a clean uniformity test.

Because ZOMBI2 streams are per level ([reproducibility](https://aadavin.github.io/zombi2/docs/reproducibility/)), the species tree is identical between a signal world and its N0 at the same seed. Gene-level draws differ once groups change (**verify**).

**N1, profile-preserving relabel (MR-18).** Take a *signal* world. For each planted-module family A, choose a non-module family A′ whose phyletic profile is within Hamming distance d (default 10% of genomes) and whose transfer count is within ±25%. Swap the labels of A and A′ in every genome, without moving anything. Presence profiles, phylogeny and HGT pattern are kept per label; the adjacency coupling between module labels is broken. O-STATS p-values for the relabeled pairs must be uniform.

**N2, multiple seeds.** Calibration uses ≥ 20 N0 seeds at micro scale (KS power), and the R5 counts in §5.2.

### 4.2 Why N0 is not too easy (F-07)

The null keeps every non-coupling reason for two genes to co-occur or sit adjacent:

1. **Vertical inheritance of order.** Rearrangement rates are tuned so the adjacency-conservation curve over tree distance matches real data (§7).
2. **Segmental HGT.** Geometric extents with a mean of 2–4 genes move neighbors together, including from dead lineages ([Zombi paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC7031779/)). Transfers decay with distance and are clade-restricted ([API](https://aadavin.github.io/zombi2/docs/api/genomes/)).
3. **Uneven sampling.** Clonal bursts, shared BioProjects, exact and 99.9% duplicates.
4. **Hub families, vehicles with hitchhiking cargo, prophages, plasmids.**
5. **Promiscuous domains, tandem paralogs, dark families.**
6. **Measurement noise.** Fragmentation and edge censoring, MAG dropout, call errors, contamination, taxonomy mislabels.

### 4.3 Difficulty audit (a generator gate, run on every N0 build at R2)

A null that a naive method passes is too easy.

- **Positive control for confounding.** On N0-full, a naive statistic must be visibly miscalibrated. The naive statistic is a hypergeometric co-occurrence or adjacency test with no clade conditioning and no vehicle or BioProject dedup. Miscalibrated means genomic inflation λ ≥ 1.5, or ≥ 5× excess of pairs at p < 10⁻³. If the naive test is calibrated, the world lacks structure and the build fails.
- **Ablation.** Switching off each confounder group (2, 3 and 4 in §4.2) one at a time must lower the naive inflation. This shows each one contributes.
- **The real statistic.** O-STATS, O-ANCHORSTATS and O-NBRSCORE must then be calibrated on the same N0. The B-NULL criteria are a KS test, QQ against permutation, and FDR within binomial bands (`PLAN.md` §5 B-NULL).

---

## 5. Parameters and presets

All parameters live in `presets/<name>.yaml`. Each one is tagged `asserted`, `fitted` (from golden genomes, §7) or `derived`, as `params.yaml` is in `PLAN.md` §4.3.1. The resolved YAML is hashed into the manifest.

### 5.1 `micro` (R2, every merge; deterministic; synthworld ≤ 3 min, full build < 15 min on a laptop)

**Genomes and families**

| parameter | value |
|---|---|
| synthetic genomes | 46 extant tips, including one clonal burst of 4. Plus 2 exact and 2 near (99.9%) duplicates, for ~50 genome records. The three real genomes (λ, T4, *M. genitalium*) are added by the corpus builder, not by `synthworld`. |
| genome size | chromosome mean ≈ 1,500 genes, range 1,000–4,000 (tip-level origination/loss balance); 0–2 plasmids of 20–200 genes |
| model panel | 200 Pfam HMMs for planted families, plus 100 held-out HMMs for dark families (`dark_fraction` 0.3); 20 Rfam CMs |
| genetic codes | one clade with code 4, one tip with code 25 |

**Planted features**

| feature | count |
|---|---|
| modules | 25 (sizes 2–6; strength grid 3 × 3, §5.3); plus 3 `profile_only` and 2 NOGD |
| vehicles | 3, each with 3–8 cargo families |
| hubs, decoys | 5 hubs, 5 decoys |
| promiscuous domains | 3 |
| fusions, splits | 15 fusions, 5 splits |
| RNA motif links | 20 (10 CMs × 2 families) |
| phages | 12 (40–80 genes each); 4 integrated as prophages |
| CRISPR arrays | 10 arrays, 5–40 spacers each |
| call errors | ~200 |
| origin-spanning genes | ≥ 5 |

**Assembly and artifacts**

| parameter | value |
|---|---|
| assembly mix | 40% complete; 40% draft (N50 50 or 200 kb); 20% MAG (N50 20 kb, completeness 70–95%) |
| contamination | 3 genomes |

**Seeds.** The default seed is fixed and committed. CI also builds one rotating seed per nightly run (R3), so seed-specific luck is caught.

### 5.2 `stage` (R5 calibration at S2C, S5, S6; B-SIM)

| parameter | value |
|---|---|
| genomes | 5,000 (B-SIM), or matched to the stage's genome count when calibrating there |
| genome size | 3,000 genes mean |
| model panel | full Pfam and Rfam panels; families ~30k |
| planted modules | 500 over the strength grid |
| traps | scaled ×20 |
| seeds | 5 signal + 10 N0 |
| runtime | ZOMBI2 for 5,000 tips is not benchmarked; 3,000 tips took 6.9 s ([performance](https://aadavin.github.io/zombi2/docs/performance/)) |
| S4 | parallel per family, deterministic (§6.1); estimated 1–3 core-hours (**verify** at S−1) |
| host | stage host |

### 5.3 Difficulty knobs (all in YAML, all logged in truth)

**Sequence**
- `target_identity` per family, as a distribution: easy ≥ 60%, medium 35–60%, hard 25–35%.
- `root_source`; `dark_fraction`.
- κ_syn; GC drift σ.

**Order and transfer**
- Rearrangement rates.
- Transfer rate and extent; distance decay; replacement probability.
- Hub transfer multiplier.

**Module strength**
- Grid: `together` ∈ {0.3, 0.7, 1.0} × `fates` ∈ {per-family loss, whole-group} × origin depth ∈ {young, old}.
- Strength classes are binned from these.

**Sampling**
- Burst size and number; mislabel rate.

**Assembly and calling**
- N50; MAG completeness; contamination rate.
- Call-error rate per kind.

**Elements**
- Motif `--exp` and clade height h.
- Spacer mismatches 0–5; PAM presence.

---

## 6. Determinism and versioning

### 6.1 Seeds

- One master `seed` (u64) per world.
- Each stage and entity gets `numpy.random.SeedSequence(seed, spawn_key=blake3("<stage>/<entity_id>")[:16] as ints)`. Results therefore do not depend on worker count or iteration order. This is the same property ZOMBI2 gets from per-family streams ([gene-family guide](https://aadavin.github.io/zombi2/docs/guide/genomes/)).
- External tools get derived non-zero 31-bit seeds:
  - ZOMBI2 `seed=`;
  - `cmemit --seed`;
  - `hmmemit --seed` in cross-checks.

  Seed 0 is never passed, because it means arbitrary in HMMER and Infernal ([hmmemit](https://www.mankier.com/1/hmmemit), [cmemit](https://manpages.debian.org/trixie/infernal/cmemit.1.en.html)).
- ZOMBI2 runs in one fixed engine. Ordered runs are serial (§1.4 item 4); the phage run uses the same engine. Parallel and serial ZOMBI2 runs differ at the same seed ([reproducibility](https://aadavin.github.io/zombi2/docs/reproducibility/)).

### 6.2 Outputs and hashing

- FASTA is written uncompressed with 80-column wrapping. Any gzip uses `mtime=0`.
- Parquet uses fixed row-group sizes, sorted rows and no wall-clock metadata.
- `manifest.json` records:
  - `synthworld_version` and git commit;
  - the preset name and resolved-params digest; the master seed;
  - tool versions: ZOMBI2, numpy, pyhmmer, Infernal, HMMER, Pfam release, Rfam release;
  - the embedded `run.zombi2` report;
  - sha256 for every file;
  - `world_digest` = BLAKE3-128 over the sorted (path, sha256) list.
- **Determinism tests (R1/R2).**
  - Build twice: identical `world_digest`.
  - Build with 1 and 8 workers: identical.
  - The macOS-vs-Linux digest is checked nightly.

### 6.3 Versioning

`synthworld` uses semver with a world-specific meaning:
- **MAJOR**: any change that maps a seed to a different world, including any RNG-stream change or a ZOMBI2 upgrade. ZOMBI2 itself warns that seeds change across versions ([reproducibility](https://aadavin.github.io/zombi2/docs/reproducibility/)).
- **MINOR**: new features that are off by default and leave existing presets byte-identical.
- **PATCH**: no output change.

An R1 test pins `world_digest` for (`micro`, default seed) per MAJOR.MINOR. If the digest changes without a version bump, that is a failure.

### 6.4 Entry into `test_run`

Add three columns to `test_run` (`TEST_ARCHITECTURE.md` §8; proposed, not edited there):
- `world_id string`;
- `world_digest FIXED_LEN_BYTE_ARRAY(16)`;
- `generator_version string` (e.g. `synthworld 1.2.0+zombi2-0.54.0`).

The `params_hash` includes `world_digest`. A gate query can then require that all blocking tests passed on the same world digest. A regression diff separates generator changes from code changes.

---

## 7. Validating the generator itself

### 7.1 Realism against the golden genomes

Comparisons use the prokaryotic chromosomes of the S0 golden set (`STAGES.md` S0): MG1655, *B. subtilis* 168, *M. genitalium*, *M. acetivorans*, *H. volcanii*, *S. acidocaldarius*, H37Rv, *B. thetaiotaomicron*, *S. pyogenes* SF370. Each statistic is computed per genome.

**Pass rule ("golden envelope").** A preset passes a statistic when both hold:
- the median of its synthetic genomes lies within [min − 10%, max + 10%] of the golden values;
- the Wasserstein-1 distance from each synthetic distribution to its nearest golden genome is ≤ the largest pairwise golden distance.

Bands are *asserted* now and fixed at S0 after the first measurement.

| statistic | definition |
|---|---|
| protein / gene length | distribution of CDS lengths (aa) |
| IGR length | distribution per junction type (co-directional, divergent, convergent), with the −1 and −4 overlap frequencies |
| coding density | fraction of replicon in CDS |
| directon size | run lengths of co-directional genes |
| strand switches | per 100 genes |
| GC, GC3 | per genome and per gene (spread) |
| codon usage | RSCU distance to the nearest golden genome; ENC |
| aa composition | Jensen–Shannon distance |
| RBS, promoter | fraction of genes with an SD-like motif upstream, as detected by pyrodigal's RBS motif output |
| pangenome shape | gene-frequency U-shape and core fraction, compared with published ranges (**verify** source at S−1) |
| adjacency decay | fraction of adjacent pairs conserved vs tree distance. Golden pairs are too few, so S−1 compares with a published conserved-adjacency curve (**verify** source) and S2C refits against real pairs. |

### 7.2 Tool parity (tools behave on synthetic data as on golden data)

| tool | statistic |
|---|---|
| pyrodigal | stop-match recall and start accuracy against `truth_gene` |
| hmmsearch | fraction of CDS with a hit to the scan panel. Planted Pfam families should be in the golden range of 60–85%; dark families should mostly have no hit. |
| MinCED | recovers `truth_spacer_target` arrays |
| cmsearch | recall of `truth_rna` by difficulty bin |

A synthetic world on which these tools do much better than on golden data is too clean.

### 7.3 Internal consistency (R1, every build)

- Every `truth_gene` slice translates to its `protein_uid`. This reuses the O-SLICE and O-TRANSLATE oracles.
- Recomputed `contig_uid` and `igr_uid` match.
- Module members with `adjacent=True` are contiguous in every extant carrier.
- Each spacer equals its protospacer with exactly `n_mismatches` differences.
- `ani_realized` lies within ±0.02% of target.
- Every `truth_edges` row maps back to the replicon.
- No trap pair appears in `truth_module_pairs`.

### 7.4 Optional realism adversary

A 4-mer and codon-pair logistic classifier trained to separate synthetic from golden 5 kb windows. Its AUC is reported, not gated. A falling AUC across versions shows realism improving.

---

## 8. Open questions and risks

### 8.1 Open questions (**verify** in the S−1 spike, §1.4)

1. **ZOMBI2 group behavior.**
   - Do group members stay co-oriented after inversions?
   - Does `copy` map to gene-tree leaves?
   - Is `parallel` really absent at the ordered resolution?
   - Are groups unavailable at the nucleotide resolution (only documented for ordered)?
2. **Genome-size stability.** Does genome size stay in band at 1,000–4,000 genes per genome with balanced rates? What is the runtime at 5,000 tips × 3,000 genes?
3. **ZOMBI2 nucleotide resolution (v2 option).** It tracks spacer homology natively ([nucleotide guide](https://aadavin.github.io/zombi2/docs/guide/genomes-nucleotide/)). It could replace §2.4's cassette model if:
   - its blocks can be painted with our codon sequences;
   - groups are supported there.

   Its genes evolve under nucleotide models only, with no codon model ([sequences guide](https://aadavin.github.io/zombi2/docs/guide/sequences/), [roadmap](https://aadavin.github.io/zombi2/docs/roadmap/)). Its letters therefore cannot be used directly for CDS.
4. **Store conventions to confirm.**
   - The `protein_uid` convention (terminal `*` or not).
   - `stop_bp` as the slot-join key.
   - The ingest manifest columns `synthworld` must emit.
5. **External calibration sources.** Which published curves serve for pangenome shape and adjacency decay at S−1.
6. **CRISPR details.** PAM rules and repeat sources for planted arrays (golden-genome repeats vs synthetic palindromes), and whether MinCED and CRISPRCasFinder treat them as real.

### 8.2 Risks

1. **ZOMBI2 churn.** It is beta, has a single maintainer, and releases several times a week ([CHANGELOG](https://github.com/AADavin/zombi2/blob/main/CHANGELOG.md)). Mitigations: an exact pin, a vendored wheel, the adapter in §2.2, and the fallback in §1.1. Its upstream test suite and cross-platform digest tests ([reproducibility](https://aadavin.github.io/zombi2/docs/reproducibility/)) lower but do not remove the risk.
2. **Circularity.** Families generated from Pfam HMMs and scanned with the same HMMs inflate recall. Mitigations: `dark_fraction`, `root_source = real_seed_member`, and reporting recall by difficulty bin (§2.3).
3. **Model simplifications.**
   - Profile-F81 ignores epistasis, covarions and structure.
   - Domain-internal indels are absent in v1.
   - Replication-strand gene-orientation bias is not modeled (ZOMBI2 has no replication origin model; **verify**).
   - Terminators are a crude stem-loop model.

   Any of these could make some statistic easier or harder than on real data. §7 measures the gap; it does not close it.
4. **Null realism.** Calibration on N0 holds only as far as N0's confounders match real ones. §4.3 checks that confounding is present; it cannot prove it is complete. B5 traps and published replays (S-01) remain necessary complements.
5. **Defining coupling in a simulation.** Planted modules are groups (physical, selected cohesion). Cargo is physically linked too. The only difference the store can see is context diversity: modules land at many loci, while cargo always shares chassis flanks. This matches the vehicle-dedup design (`PLAN.md` §4.3). It should be stated in the B-SIM report so nobody reads the simulation as testing more than that.
6. **Stage-scale cost and storage.** A 5,000-genome world is about 15 Gb of FASTA plus truth. Worlds are cached by `world_digest` on the stage host and regenerated rather than archived.
7. **Licenses.**
   - ZOMBI2 is MIT. The Zombi 1 fallback is GPL-3.0-or-later ([README](https://github.com/AADavin/Zombi)); it is acceptable only as a separate tool in the test image.
   - ALF has no license and is not used.
   - Pfam and Rfam reuse terms for emitted sequences need confirming (**verify**; believed CC0).
