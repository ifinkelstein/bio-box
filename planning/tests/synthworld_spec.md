# `synthworld`: synthetic-world generator spec

Status: **planning only** (2026-10-29). Nothing is installed or implemented. **Update 2026-11-12 (card F-02 session 1):** the backbone (§2.1 S1–S2, §2.2, §6.1, the backbone half of §5.1) is implemented in `src/biobox/synthworld/` on `zombi2==0.54.0` under DEC-15; §2.1, §2.2, §2.6, §5.1 and §6.1 carry the spike's corrections (`synthworld_spike.md` §7). **Update 2026-11-12 (card F-02 session 2):** the rest of S1–S2 is implemented: the synthetic taxonomy, plasmids, vehicles, hubs, decoys, distance-decay transfer, lineage-varying rates and the typed event and gene-tree tables (§2.1, §2.2, §5.1, §5.3 updated). **Update 2026-11-12 (card F-02 session 3):** family design and protein-coding sequences (§2.1 S3–S4, §2.3) are implemented; §2.3, §2.6, §3.2, §5.1 and §6.1 carry "as built" notes. The rest of this file is still a plan. This file settles `TEST_ARCHITECTURE.md` §11 item 1 ("Zombi fit at micro-corpus size, ALF fallback") and specifies the generator behind §6.1 (micro-corpus), §6.2 (synthetic worlds, null worlds), MR-18, and `PLAN.md` §5 rows B-SIM and B-NULL.

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

**S1–S2 as built (F-02 session 1; DEC-15; spike §7 items 1, 8 and 10).**
- *Clonal burst.* ZOMBI2 grows `n_extant − (burst − 1)` tips; `synthworld` grafts a random binary subtree of `burst` tips onto one extant tip (chosen from stream `S1/clonal_burst`), rooted `burst_depth` × crown height before the present.
- *Root scramble (DEC-15 (c)).* ZOMBI2 lays the root out with every gene on + and all declared families in one block at the end of each chromosome (spike §2.3), which plants adjacency between modules that a null does not have. `synthworld` replaces the species tree's stem with a stem of `root_scramble.stem_length` on which only inversions and transpositions run (`changing_at` skyline; every other rate, the declared families' own rates included, is 0 there). Groups move whole and are never cut. Every build checks the result against ZOMBI2's `initial_genome.tsv`: + fraction in [0.4, 0.6], ≤ 5% of non-unit initial adjacencies kept, ≤ 50% of declared genes' outside neighbours declared (derived thresholds; `micro` default seed measured 0.488, 0.29% and 3.0%).
- *Groups (DEC-15 (a), (b)).* Adjacency is guaranteed only for `fates="group"`: with `fates="family"` units break, mostly by tandem copies landing inside them (spike §2.2; on `micro` at the default seed 264 of 3,842 multi-member units, 6.9%, measured on wilkcomp02 2026-11-12 by `tests/r2/test_synthworld_backbone.py`). An adjacency-scored group must have `fates="group"`; a group with `fates="family"` is allowed only unscored. Truth is computed from the realized order in every case. `replacement` is 0 in every world with groups (ZOMBI2 refuses any value above 0 with groups).
- *`together`* is the probability that a transfer starting on a member carries the whole group. At extent > 1 the realized co-transfer share differs from the dial (spike §2.3), so truth counts `event` ids from `genome_events.tsv` and never derives transfers from the dial.

**S1–S2 as built (F-02 session 2; spike §7 items 2 and 7).** Each feature has an R2 test that recovers it from ZOMBI2's raw files (`tests/r2/test_synthworld_backbone.py`); the numbers below are from the `micro` default seed (measured, wilkcomp02, 2026-11-12, `scripts/pod_pytest.sh F-02-2 -s tests/r2/test_synthworld_backbone.py`).
- *Synthetic taxonomy.* A genome's taxon at a rank is the species-tree branch alive at that rank's cut, `crown + cut × (present − crown)` with the crown at the first speciation (the scramble stem is not counted). `micro` cuts phylum … species at 0.08, 0.2, 0.35, 0.5, 0.7, 0.9 of the crown age, so the clonal burst (rooted at 0.98) is one species. Taxa are named `p__Phylum01` … `s__Species12`, numbered per rank by the first genome that holds them, under `d__Bacteria`; `taxa.tsv` keeps the crossing branch, cut time, genome count and parent taxon as truth. This is `tax_lineage_true` of §3.2; mislabelling belongs to S11.
- *Plasmids are excised groups, not ZOMBI2 chromosomes.* ZOMBI2 deals the root genes round-robin over its chromosomes, and lands every origination and every arriving transfer on a uniformly chosen chromosome (`zombi2/genomes/ordered.py`, `_do_transfer`); translocation moves blocks between them at a per-copy rate, and `chromosome_loss` can take the main chromosome. Its replicons therefore converge in size, and it has no per-replicon rates. So each plasmid type is an adjacent group (`fates="group"`, `together=1`) of its preset size, laid out whole in the root genome; its members carry the per-plasmid transfer and loss rates divided by the size (rates are per copy, and a group is hit once per member), no duplication rate of their own, and a recipient rule that gives a lineage already carrying the plasmid weight 0 (entry exclusion; a ZOMBI2 joint run). ZOMBI2 runs one circular chromosome; `synthworld` cuts every extant unit of a plasmid group out into a circular replicon of its own (replicon 1, 2, …) and joins the chromosome around it, so the realized order and its truth have the plasmid on its own ring. Consequences: a plasmid's gene content is fixed (nothing is gained or lost inside it, ZOMBI2 inserts nothing into a unit); all its genes share one strand (a group is laid out on +); a background run that starts on a neighbour can still copy or carry it. `micro`: two types of 120 and 40 genes, transfer 0.3 and loss 0.7 per plasmid per unit time; sizes are exact (tolerance 0); genomes carry 0, 1 or 2 plasmids (17, 20 and 9 of 46).
- *Vehicles* are adjacent groups (`fates="group"`, `together=1`) of a chassis family, the cargo families and a second chassis family, with per-vehicle transfer 0.6 and loss 0.2. Every transfer event that moves a vehicle gene moves the whole vehicle, and cargo never occurs outside a unit flanked by its chassis (§8.2 risk 5). Vehicles, like plasmids, are unscored groups; whether their pairs are excluded from null evaluation is §3 and §4's to settle.
- *Hubs* are single declared families transferring at 20× the genome's per-copy rate (2.0) with loss 0.1, so they fill `max_family_size` (10 copies per genome). Their neighbour diversity, the distinct families next to any copy, is 6.0–6.6× the median root family's present in ≥ 40 genomes (233–259 against 39; the R2 test requires 3×). Per copy the gap is small (about 1.1× in the session-2 scans), because a transfer of extent 3 that starts on a hub carries its right-hand neighbours along. A hub's copy line can die out early, with probability about loss/transfer (0.05), so another seed can lose a hub.
- *Decoys* are single declared families whose transfers are weighted by the presence of a planted module's first family in the recipient (`Recipients().weighted_by("genomes:mod<m>.0", {present: 1, absent: 0})`, a ZOMBI2 joint run), so their phyletic profile follows the module without adjacency or co-transfer. Every decoy transfer lands in a lineage that carries the partner; carriers were a mean 0.64 of the other contemporaneous lineages at those moments (0.92 with the widespread partners first tried, modules 0, 5, 10, 15, 20, so the partners are the patchiest modules of the default seed).
- *Distance-decay transfer*: ZOMBI2 `Distance(decay=4)`, recipient weight exp(−4 d / mean root-to-tip time). A family with its own recipient rule (plasmids, decoys) uses that rule instead of the run's for its own genes, and a segment's rules multiply. Summed over transfer events, the recipient is closer to the donor than the other contemporaneous lineages 3.9× as often as it is farther (2,911 against 756 over 6,227 events); under `transfer_to: uniform` the ratio is 1.03 (1,460 against 1,423).
- *Lineage-varying rates* (spike §7 item 7): two ZOMBI2 axes `Random('lineages', Drift(LogNormal(0, 0.8)))`, one on loss and one shared by duplication, transfer and origination. With the spike's rates a genome's size relaxes over about 20 time units against a crown of about 4, so lineage rate differences cannot move sizes far: none of the configurations scanned at the spike's rates, or with σ ≤ 0.6, reached both a chromosome of ≤ 1,100 genes and one of ≥ 3,600. `micro` therefore doubles duplication, transfer, loss and origination (0.04, 0.1, 0.24, 0.1 × root families). Chromosomes then hold 969–3,798 genes, mean 1,873 (+25% on the §5.1 mean of 1,500; the session-1 test's band is ±30%). Root families present in ≥ 40 of 46 genomes number 243.
- *Typed tables* (§2.2): `events.tsv` is `genome_events.tsv` without its speciation rows, with the acting lineage and a transfer's recipient read from the copy labels; `gene_tree_nodes.tsv` is every complete gene tree, one row per node in preorder.

### 2.2 Backbone interface (the only ZOMBI2 coupling)

`synthworld.backbone` reads five ZOMBI2 outputs into typed tables: species tree, `gene_order.tsv`, gene trees (`gene_trees.tsv` with `--compact`), `genome_events.tsv` and `groups.tsv` ([output files](https://aadavin.github.io/zombi2/docs/output-files/)).

An R1 schema test pins the columns. Any adapter (Zombi 1 `_GENOME.tsv`, or the in-house fallback) must produce the same tables. The raw ZOMBI2 directory is kept under `worlds/<id>/zombi2/` for audit.

**As built (F-02 session 1; spike §7 items 4–6).** `biobox.synthworld.backbone` is the only module that imports ZOMBI2. `build(preset, seed, out)` writes `out/zombi2/` (`OrderedGenomesResult.write(compact=True)`) and `out/backbone/`: `genomes.tsv`, `genes.tsv` (the realized order of every extant genome with strand, family, copy label, unit), `families.tsv`, `groups.tsv`, `truth_units.tsv`, `truth_adjacency.tsv` and `backbone.json`; the module docstring fixes the columns. Pinned ZOMBI2 0.54.0 headers (`backbone.PINNED_HEADERS`, checked on every build):
- `gene_order.tsv`: `lineage chromosome topology position strand family copy`; `copy` is `g<id>`, `family` an integer id;
- `genome_events.tsv`: 11 columns, `time kind family parents children chromosome start length dest_chromosome dest_position event`; transfer kinds are `transfer_additive` and `transfer_replacing`;
- `rearrangement_events.tsv`: `time kind lineage chromosome start length dest_chromosome dest_position flipped cuts`;
- `groups.tsv`: `group unit family gene`, `gene` a bare integer; `adjacent=False` groups write no rows; declared names map to ids through `genome_summary.json` `families.names`;
- `gene_trees.tsv`: `family complete extant`; the leaf label is `<lineage>_<copy>`;
- `initial_genome.tsv`: `chromosome topology position strand family copy`.

**As built (F-02 session 2).** The format is `synthworld-backbone/2`. Added files: `replicons.tsv` (chromosome and excised plasmids), `taxonomy.tsv` and `taxa.tsv` (§2.1 S1), `events.tsv` (`genome_events.tsv` typed, speciation rows left out: the gene trees hold them) and `gene_tree_nodes.tsv` (the complete gene trees of `gene_trees.tsv`, one row per node in preorder, leaves marked `extant` or `lost`). `families.tsv` gains a last column, `partner` (a decoy's partner family); the existing columns of every session-1 file are unchanged, and `genes.tsv` now numbers plasmids as replicons 1, 2, … . The module docstring fixes every column.

Extancy comes from the species tree's fates (`gene_order.tsv` rows of extant lineages) or from the `extant` gene trees, never from an `n` prefix in a `complete` gene tree: those also label lost copies. Copy ids are re-minted at every speciation; a unit id is kept across speciations, so a unit is keyed by (genome, unit). Groups can be declared only from the Python API, which writes no `run.zombi2`; `synthworld` records its own provenance in `backbone.json` (ZOMBI2 and numpy versions, the resolved preset and its digest, the master seed and derived tool seeds). Where a CLI `run.zombi2` or `.log` exists it holds timestamps and absolute paths and is not byte-stable, so it is not used.

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

**As built (F-02 session 3; spike §7 item 6).** `src/biobox/synthworld/{profiles,codons,sequences}.py`; `sequences.build(world)` reads the backbone (never ZOMBI2's files: the species tree comes through `backbone.species_branches`) and writes `world/sequences/` (format `synthworld-sequences/1`: `cds.fna`, `proteins.faa`, `roots.faa` with `.fai`, and `families.tsv`, `domains.tsv`, `genes.tsv`, `lineages.tsv`, `panel.tsv`, `sequences.json`); the module docstring fixes the model and every column. Where it settles or departs from the text above:
- *Panel.* `profiles.source: pfam` reads the panel from a Pfam HMM file with pyhmmer (HMMs selected in file order by a per-name draw, so only the head of the file is read; dark and promiscuous members by further draws). The `micro` default is `synthetic`, because CI has no Pfam file: per-site rows Dirichlet(A × background) with match lengths and A fitted to Pfam 38.2 (measured, `micro.yaml` gives the numbers), and Pfam's median transitions. The 3 promiscuous traps are panel members, not named HTH or P-loop models.
- *Root.* `sample_domain` walks the core model with inserts *and deletes* from the transitions, as `hmmemit` does; "per match site" above would leave deletes out and score above `hmmemit` (the §2.3 cross-check is the test). Insert residues come from the background. Tails and linkers are geometric (means 20, 25, 20 in `micro`); the initiator M is a site of its own that never substitutes.
- *Substitution and codons* as written. c_b is an uncorrelated log-normal per species branch (σ 0.3, mean 1); Γ is Yang's 4-category mean discrete gamma (α 1); κ_syn 3 applies to every codon, the initiation and stop codons included (redrawn from fixed ATG/GTG/TTG and TAA/TGA/TAG weights). Codon usage per lineage: an amino acid's codons weighted GC3 or 1 − GC3 by their third base, with logit GC3 Brownian on the species tree (σ 0.25 per √time from 0.5 at the stem), constant along a branch. Genetic code 11 only (codes 4 and 25 stay with S8). The process runs on the gene tree pruned to the ancestors of extant copies (the extant copies' distribution is unchanged; dead branches are not drawn).
- *Indels* as written (tails and linkers, 1–10 aa), rate 0.05 per site per expected substitution; within a branch substitutions come before indels.
- *Divergence control.* mu_f solves the closed form over (up to 256) pairs of extant copies, or copy-to-root for a single copy, over the root's tail, linker and domain sites; a target below saturation is raised to saturation + 0.01. The "target bitscore band" alternative is not built: see the next item.
- *Detectability does not fall with divergence under this model.* The equilibrium of every site is its emission row and the root is drawn from it, so every copy is on its own a draw from the profile, whatever its divergence; divergence lowers identity between copies, not a copy's score against the source HMM. No choice of mu_f can therefore place a family in a bitscore band, and a Pfam family's copies stay detectable at any depth: tests/r2/test_synthworld_sequences_pfam.py measures it. Making detection fall with divergence needs a different model (for example a root sharper than the profile, an equilibrium relaxed toward the background, or `real_seed_member` roots); that is a decision for the PI. Measured (wilkcomp02, 2026-11-12, `micro` default seed with `profiles.source: pfam` on Pfam 38.2, `BIOBOX_POD_PYDEPS=… pytest -s tests/r2/test_synthworld_sequences_pfam.py` in `bio-box-test_74b3e05.sif` with Pfam bound read-only): over 1,928 copies of 40 Pfam-seeded families, `hmmsearch --cut_ga` found 596/698 copies within 0.5 substitutions per site of the root, 174/239, 54/87, 157/178, 131/179, 184/208 in the bins to 3.0, and 261/339 beyond (overall 81%, inside §7.2's golden 60–85%); 29 of 40 roots pass GA; the far half of a family's copies scored lower than the near half in 19 of 40 families; a copy agrees with its root on detection in 0.883 of the nearest third and 0.807 of the farthest, and differs from the root's score by a median 5.4 against 14.3 bits. Cross-check: for 20 scan-panel HMMs the numpy root sampler and `hmmemit -N 500` give KS p from 0.012 to 0.995 (all above 0.01).
- *Spike §7 item 6.* Extant copies are the backbone's `extant` gene-tree leaves (from `genes.tsv`), never an `n` prefix; copy ids are re-minted at each speciation, so `copy_label` names the extant leaf only, and a copy's path is read from the gene tree.
- *Measured on `micro`* (wilkcomp02, 2026-11-12, default seed, `scripts/pod_pytest.sh F-02-3 -s tests/r2/test_synthworld_sequences.py`): 2,698 families and 89,102 genes; realized pairwise domain identity against the closed form in quintiles of expectation 0.286/0.286, 0.456/0.458, 0.618/0.619, 0.775/0.774, 0.915/0.916 over 6,355 pairs, Spearman(mu T, identity) −0.987; unrelated families' identity above chance by 0.0001; non-domain composition at Jensen–Shannon 0.00002 nats from the background; protein length median 277 aa (5% 88, 95% 677); genome GC3 against lineage GC3 r 0.988.

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
| ZOMBI2 ordered run | seconds (**verify**, §1.4). Measured: the whole `micro` backbone (species tree, burst graft, scramble stem, ZOMBI2 run and write, adapter tables, digests) took 17.5 s and 19.6 s wall in two runs (wilkcomp02, 2026-11-12, one process inside `bio-box-test_74b3e05.sif` at nice 19 with the host near its 40-thread guest cap; `scripts/pod_pytest.sh F-02 -s tests/r2/test_synthworld_backbone.py`, `test_micro_backbone_runtime_within_budget`). After session 2 (plasmids, vehicles, hubs, decoys in a ZOMBI2 joint run, doubled turnover rates, typed event and gene-tree tables): 36.2 s and 36.5 s wall in two runs, same host and command (`F-02-2`), same conditions |
| S4 evolver: about 3,000 families × ~300 sites × ~200 tree nodes, vectorized | under 2 min single core (estimate, measured at S−1). Measured (F-02 session 3, S3 and S4 together, synthetic panel): 2,698 families, 206,583 pruned gene-tree branches, 89,102 genes in 75.3 s and 78.8 s wall in two runs (wilkcomp02, 2026-11-12, one process at nice 19 inside `bio-box-test_74b3e05.sif`, the account near its guest cap; a direct `sequences.build` call and `test_micro_sequences_runtime_within_budget`); peak RSS of the stage in a child process 430 MB. With the backbone, synthworld so far takes about 115 s against the 3 min target |
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
| copy_label | string | ZOMBI2 extant leaf `n<sp>_g<copy>`; extancy from the extant trees or `gene_order.tsv`, never the prefix (lost copies carry it too), and copy ids change at each speciation (spike §7 item 6) |
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
| genome size | chromosome mean ≈ 1,500 genes; 0–2 plasmids of 20–200 genes. The range 1,000–4,000 first written here is not what one root genome with uniform rates gives: the spike measured about ±9% across tips (spike §7 item 7), and the F-02 `micro` backbone at the default seed has mean 1,509, range 1,393–1,708 (measured, wilkcomp02, 2026-11-12, `test_micro_backbone_runtime_within_budget`). Session 2 reaches the range with lineage-varying rates and doubled turnover rates (§2.1 S1–S2 as built, session 2): chromosomes of 969–3,798 genes, mean 1,873, at the default seed (measured, wilkcomp02, 2026-11-12); the R2 test requires ≤ 1,100 and ≥ 3,600 (within 10% of the range's ends) and a mean within 1,500 ± 30%. Plasmids are excised ZOMBI2 groups of exact preset size (§2.1). |
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

**Backbone rates (F-02 session 1).** §5.1 first gave no DTL or rearrangement rates (spike §7 item 7). `src/biobox/synthworld/presets/micro.yaml` now holds them, each tagged: the spike's rates (duplication 0.02, transfer 0.05, loss 0.12, inversion 0.05, transposition 0.02 per copy; origination 0.05 × root families per lineage; geometric extents 2, 2, 3, 5, 3; `inversion_probability` 0.5; uniform recipients; replacement 0), 1,000 root families on one circular chromosome, the root scramble of §2.1, and the modules: 25 adjacent groups of sizes 2–6 at `fates="group"` on a 3 × 3 grid of `together` {0.3, 0.7, 1.0} × member transfer rate {0.05, 0.15, 0.30} (member loss 0.03), plus 3 `profile_only` groups (`adjacent=False`) of 3. Because DEC-15 (a) fixes `fates`, the member transfer rate replaces §5.3's "per-family loss" axis in the grid; the "origin depth" axis is not available for adjacent groups (spike §7 item 3). Vehicles, hubs, decoys and NOGD are not in the backbone yet.

**Backbone parameters (F-02 session 2).** Added to `micro.yaml`, each tagged: the taxonomy cuts; `transfer_to: distance` with decay 4; two lineage axes at σ 0.8 (loss; duplication, transfer and origination); duplication, transfer, loss and origination doubled from the spike's (0.04, 0.1, 0.24, 0.1 × root families), because at the spike's rates no lineage spread tried reached the size range (§2.1); plasmids of 120 and 40 genes (transfer 0.3, loss 0.7 per plasmid, entry exclusion); 3 vehicles with 3, 5 and 8 cargo families (transfer 0.6, loss 0.2 per vehicle); 5 hubs (transfer 20× the genome rate, loss 0.1); 5 decoys on modules 1, 9, 11, 19 and 22 (transfer and loss 0.4, weight 0 without the partner). NOGD stays for S7. The strength grid is unchanged, waiting for the PI.

**Sequence parameters (F-02 session 3).** `micro.yaml` gains a `sequences` section, each value tagged: the panel (source `synthetic` by default, 200 scan + 100 dark + 3 promiscuous, the Pfam selection rule, and the synthetic panel fitted to Pfam 38.2 measurements), the architecture (1, 2 or 3 domains at 0.6, 0.3, 0.1; `dark_fraction` 0.3; promiscuous fraction 0.1; tail and linker means), the evolution (§5.3 difficulty classes at weights 0.5, 0.35, 0.15, easy up to 0.95; Γ α 1 in 4 categories; clock σ 0.3; κ_syn 3; indel rate 0.05; rate bounds) and the codons (code 11, GC3 0.5 at the stem with σ 0.25, start and stop codon weights near *E. coli*'s).

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
- Transfer rate and extent; distance decay (`genomes.transfer_decay`). There is no replacement-probability knob in practice: ZOMBI2 refuses `replacement` above 0 in any run with groups (spike §7 item 2), and since session 2 plasmids and vehicles are groups, so only a world without modules, plasmids, vehicles and Cas groups could replace. N0-strict (§4.1) drops vehicles and Cas groups but keeps plasmids, so it too has `replacement` 0 unless it also drops plasmids, which would make it differ from N0-full in more than the planted groups.
- Hub transfer multiplier (`hubs.transfer_multiplier`); lineage rate spread (`lineage_rates` σ).

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
- ZOMBI2 runs in one fixed engine. `simulate_genomes_ordered` has no `parallel` argument and the CLI refuses `--parallel` at the ordered resolution, so ordered runs are always serial (spike §4 and §7 item 9); the phage run uses the same engine. Parallel and serial ZOMBI2 runs differ at the same seed only at the family resolution ([reproducibility](https://aadavin.github.io/zombi2/docs/reproducibility/)). macOS-vs-Linux digest equality is still unchecked.
- **As built** (`src/biobox/synthworld/seeds.py`): "as ints" means the 16 bytes read as four big-endian uint32 words; a tool seed is the stream's first uint32 taken mod (2³¹ − 1), plus 1, so it lies in [1, 2³¹ − 1]. Component paths in use: `S1/species_tree` (ZOMBI2 species seed), `S1/clonal_burst` (burst tip and split times), `S2/genomes` (ZOMBI2 genomes seed). Renaming a path changes every world (a MAJOR bump, §6.3). Session 2 adds no path: the lineage multipliers are drawn by ZOMBI2 from its `S2/genomes` seed, and the taxonomy and plasmid excision draw nothing. The `micro` backbone is byte-identical across two builds and across a subprocess with another `PYTHONHASHSEED` and thread count (`tests/r2/test_synthworld_backbone.py`). Session 3 adds `S3/panel`, `S3/panel/<NAME>`, `S3/family/<family_id>`, `S4/lineages` and `S4/family/<family_id>` (listed in `seeds.py`); written floats are rounded to 6 significant digits and every dictionary of the preset is read in key order, so the order a preset was loaded in (YAML or the JSON copy in `backbone.json`) cannot change a draw. The sequences are byte-identical in a subprocess with another `PYTHONHASHSEED` and thread count (`tests/r2/test_synthworld_sequences.py`).

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
