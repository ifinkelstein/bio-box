# 0-08 Non-coding RNA: which Rfam models, where, with which flags, and which overlapping hits to keep

Proposal for the PI's decision card 0-08 (TASKS.md), 2026-11-11. Design only: nothing under `src/` or
`tests/` changes here. Edit the rules and tick the choices in section 6; the cards in section 6 follow
from the ticks.

Sources. The measurement this rests on is 0-09c, `planning/measurements/2026-11-11_infernal_intergenic.md`
(Infernal 1.1.5, Rfam 15.1, 8 golden genomes, 18,106,231 bp). Its per-model CPU lines and truth hits are on
the pod under `/tmp/if2453/m09c/out` and were read, not re-run. Files beside this note:

- `0-08_rfam15.1_models.tsv`: the header of every model in Rfam 15.1 `Rfam.cm` (ACC, NAME, CLEN, W, GA,
  TC, NC, NSEQ, consensus base pairs, filter-HMM MAXL, DESC). *Measured* on wilkcomp02, 2026-11-11, by
  `0-08_cmheaders.py` in `bio-box_80d4c55.sif` inside a `systemd-run --user --scope -p MemoryMax=4G -p
  MemorySwapMax=0 -p CPUQuota=100%` scope at nice 19, reading `/tmp/if2453/m09c/lib/Rfam.cm` (the
  gunzipped refdb.yaml file). 4,227 models, all with GA, 347 without base pairs (as 0-09c found); sha256
  `40655bea60e5435aa06b00c2a0786b1bf5015eb9ec4f606ab995aae917673c95`.
- `0-08_model_costs.tsv`: per model, its truth hits in the 8 genomes and its own `# CPU time` line (user +
  system) in the 0-09c runs `wgA_default` + `wgB_default` (whole genome), `igr_m200_notrunc` and
  `igr_m0_notrunc`. *Measured* values, joined by `0-08_model_rules.py join`.
- `0-08_model_rules.py rules`: applies the rules below and prints every count and cost in this note.

Every cost of a model set below is *derived* by one formula:
**CPU-hours per 10 Gbp = (sum of the set's per-model CPU-s in the named run) / 18,106,231 × 1e10 / 3,600.**
The intergenic runs cover the same genomes, so the margin's saving is inside the number. The per-model
lines sum to 0.4% less than the processes' CPU (12,403.2 against 12,450.1 CPU-s), so the full library reads
1,903 here and 1,910 in 0-09c. 0-09c checked that per-model lines predict a subset run within 0.7%; for the
intergenic runs that prediction is untested (no subset was run on intergenic records).

## 1. The cost problem

1. Today the producer (card E-01) scans whole genomes with all 4,227 models: **1,910 CPU-hours per 10 Gbp**
   (*derived* in 0-09c §1 from 1,454 bases per CPU-second, *measured*), 4.3 times Pfam's 443 (*derived*,
   card 0-09).
2. Cost adds up model by model: each model's own CPU line in the full run predicts the 258-model run within
   0.7% and the 20-model run within 1.4% (*measured*, 0-09c §3b). So a model set can be priced before it runs.
3. **3,969 of the 4,227 models hit nothing in the 8 genomes and carry 84% of the CPU** (*derived*: 1 −
   2,029.4 / 12,403.2 CPU-s, whole-genome per-model lines).
4. Region is the second lever: intergenic regions with a 200-base margin and `--notrunc` cost 796 at 98.0%
   hit recall (0-09c §2). The levers multiply only in part, because the long rRNA models that any set keeps
   cost the most per intergenic record.

## 2. The model set, as rules

A model takes the first rule that matches it, in the order (a), (b), (c), (d). Rules (b) and (c) read only
the model's NAME and DESC (a heuristic; see "Where descriptions mislead" below). The rule code is in
`0-08_model_rules.py` (`INCLUDE`, `EXCLUDE`); the PI can edit a rule and run it again.

- **(a) A hit in prokaryotes or phages in the measurement:** every model with a hit at or above GA in the 8
  genomes: 258 models (Appendix A: name, accession, description, hits, CPU share). **Exception (a−):** drop
  the four eukaryotic rRNA models, SSU_rRNA_eukarya RF01960, SSU_rRNA_microsporidia RF02542,
  LSU_rRNA_eukarya RF02543 and 5_8S_rRNA RF00002. Their 87 hits all lie on loci that a bacterial or archaeal
  rRNA model also hits (*derived* from the truth hits: every one is in a same-strand overlap cluster with
  SSU_rRNA_bacteria/archaea or LSU_rRNA_bacteria/archaea). They add no locus, and they cost 467.9 of the
  set's CPU-s: 72 CPU-hours per 10 Gbp whole genome and 49 intergenic.
- **(b) Included whatever 8 genomes show** (Appendix B lists every accession): (b1) rRNA and (b2) tRNA and
  tRNA-Sec (all already in (a)); (b3) RNase P, bacterial and archaeal; (b4) tmRNA and its variants; (b5)
  bacterial and archaeal SRP; (b6) 6S; (b7) CRISPR direct repeats and tracrRNA; (b8) group I and II introns;
  (b9) riboswitches and T-box; (b10) sRNAs whose description names a model prokaryote (*E. coli*,
  *Salmonella*, *Shigella*, *Bacillus*, *Mycobacterium*, *Staphylococcus*, *Pseudomonas*, *Vibrio*,
  *Listeria*, *Streptococcus*, *Caulobacter*, *Clostridium*, *Neisseria* and a few more); (b11) the archaeal
  box C/D sRNAs, which Rfam describes as "Small nucleolar RNA sR*n*".
- **(c) Excluded, eukaryote-only by description:** (c1) snoRNA and scaRNA; (c2) microRNA; (c3) spliceosomal
  snRNA and U7; (c4) lncRNA (named human and mouse lncRNAs); (c5) vault; (c6) telomerase. **Optional:**
  (c7) other eukaryotic housekeeping RNAs (nuclear RNase P, MRP, eukaryotic SRP, 7SK, Y RNA); (c8) elements
  of eukaryotic viruses; (c9) mitochondrial and plastid RNAs.
- **(d) Long tail:** every other model, 1,040. Most are bacterial sRNAs and motifs, by description. They
  cannot be excluded on what we know, and scanning them all costs a quarter of the library. Proposal: scan
  the whole library once on a sample of the corpus (measurement M1, section 7) and add every (d) model that
  hits a prokaryote or phage. The same scan checks (c): a (c) model that hits there shows a rule is wrong.

| rule | models | CPU-s, whole genome | share | CPU-h/10 Gbp, whole genome | intergenic 200, `--notrunc` | intergenic 0, `--notrunc` |
|---|---|---|---|---|---|---|
| full library | 4,227 | 12,403.2 | 100.0% | 1,903 | 792 | 308 |
| (a) a hit in the 8 genomes | 258 | 2,029.4 | 16.4% | 311 | 176 | 112 |
| (b3) RNase P, bacterial and archaeal | 1 | 3.5 | 0.0% | 1 | 0 | 0 |
| (b4) tmRNA | 4 | 91.5 | 0.7% | 14 | 11 | 9 |
| (b5) SRP, bacterial and archaeal | 1 | 1.4 | 0.0% | 0 | 0 | 0 |
| (b6) 6S | 2 | 4.5 | 0.0% | 1 | 0 | 0 |
| (b7) CRISPR direct repeats, tracrRNA | 61 | 60.3 | 0.5% | 9 | 4 | 2 |
| (b8) group I and II introns | 6 | 16.6 | 0.1% | 3 | 1 | 0 |
| (b9) riboswitches, T-box | 33 | 58.1 | 0.5% | 9 | 3 | 1 |
| (b10) sRNAs named for a model prokaryote | 188 | 1,010.4 | 8.1% | 155 | 51 | 14 |
| (b11) archaeal box C/D sRNAs | 38 | 43.0 | 0.3% | 7 | 3 | 1 |
| (c1) snoRNA, scaRNA | 733 | 2,048.7 | 16.5% | 314 | 129 | 40 |
| (c2) microRNA | 1,597 | 3,038.2 | 24.5% | 466 | 194 | 64 |
| (c3) snRNA | 16 | 48.7 | 0.4% | 7 | 3 | 1 |
| (c4) lncRNA | 68 | 78.6 | 0.6% | 12 | 5 | 2 |
| (c5) vault | 1 | 1.9 | 0.0% | 0 | 0 | 0 |
| (c6) telomerase | 4 | 46.7 | 0.4% | 7 | 2 | 0 |
| (c7) other eukaryotic housekeeping | 12 | 76.0 | 0.6% | 12 | 4 | 1 |
| (c8) eukaryotic viruses | 129 | 329.0 | 2.7% | 50 | 19 | 6 |
| (c9) organellar | 35 | 99.5 | 0.8% | 15 | 7 | 3 |
| (d) long tail | 1,040 | 3,317.3 | 26.7% | 509 | 178 | 53 |

The set at each step (same formula; the last row adds section 4's whole-genome scan of RNAs inside coding
sequence):

| step | models | CPU-h/10 Gbp, whole genome | intergenic 200, `--notrunc` |
|---|---|---|---|
| 0. today: full library | 4,227 | 1,903 (0-09c: 1,910) | 792 (0-09c: 796) |
| 1. (a) | 258 | 311 | 176 |
| 2. (a) + (b) | 592 | 509 | 250 |
| **2'. (a−) + (b): proposed** | **588** | **437** | **201** |
| 2''. as 2' without (b10) | 400 | 282 | 150 |
| 3. full library without (c1)-(c6) | 1,808 | 1,095 | 459 |
| 3'. full library without (c1)-(c9) | 1,632 | 1,018 | 428 |
| **2' intergenic, with its 23 in-CDS models whole genome (section 4)** | **588** | | **247** |

After M1, step 2' gains the (d) models that hit. They add between 0 and 178 CPU-hours per 10 Gbp
intergenic (*derived*: the (d) row), so the proposed set ends between 247 and about 425.

**Where descriptions mislead (found while writing this note).** As first written, rule (c1) also matched
the 48 archaeal box C/D sRNAs (sR1 to sR60, snoPyro_CD), 10 of which hit *S. acidocaldarius* here. A
"Y RNA" pattern in (c7) matched *E. coli* GlmY ("Glm Y RNA"), the bacterial Y RNA-like YrlA and YrlB, and
the *Deinococcus* Y RNA. Both are fixed by hand in the rules above (b11; c7 matches only the model named
Y_RNA). The error also runs the other way: (b10) misses enterobacterial sRNAs whose description names no
organism. isrG ("isrG Hfq binding RNA", 164.6 CPU-s, the second most expensive model in the library, no
hit) stays in (d). So the description rules are not exact, and that is the point of M1. Two Rfam files
would make them exact: `family.txt` gives each family's type (the `TP` field: Gene; rRNA, tRNA, sRNA,
snoRNA, miRNA, riboswitch, ...), and `Rfam.clanin` gives clan membership (needed for section 5). Neither is
fetched; section 6 lists them for `refdb.yaml`.

## 3. tRNA and rRNA: Rfam's models or dedicated tools

Neither tRNAscan-SE nor barrnap is in the image: `container/environment.yml` pins python, uv, mmseqs2,
hmmer 3.4, infernal 1.1.5, samtools and prodigal (an oracle only). Adding either is a PI change to that file
and a new image.

- **tRNA.** RF00005 and RF01852 together cost 9.0 CPU-s of 12,403 (0.07%; *measured* lines), so a dedicated
  tool cannot save time. It could buy sensitivity and annotation. The E-02 panel loses 2 of its 5 planted
  tRNA draws (RF00005 at 18758-18837 and 29073-29145) to Infernal's filters at the fixed Z (DEC-04); both
  are found with `--max` (*measured*, `tests/r3/test_g04_spikeins.py`; 0-09c §4 shows `--rfam` changes
  nothing). On real genomes the loss may be smaller. Truth hits per genome (*measured*, 0-09c): *E. coli*
  86 tRNA plus 1 tRNA-Sec, *B. subtilis* 86, *M. tuberculosis* 44, *S. pyogenes* 60, *S. acidocaldarius* 47,
  T4 8, LAK-A1 32, lambda 0. The counts usually reported for MG1655 and strain 168 are 86 each (*estimate*:
  recalled, not fetched, unverified). tRNAscan-SE would add the anticodon, isotype and pseudogene calls that
  Rfam's single tRNA model does not give. Its cost here is not measured.
- **rRNA.** The 9 rRNA models cost 857.1 CPU-s of 12,403 whole genome (6.9%) and 630.8 of 5,162.7 at
  intergenic margin 200 with `--notrunc` (12.2%) (*measured* lines). LSU_rRNA_eukarya alone costs 315.0
  CPU-s, the most of any model. In the proposed set 2' the five bacterial and archaeal rRNA models are 13.7%
  of the whole-genome cost and 23.8% of the intergenic cost (*derived*). barrnap, an HMM search with nhmmer
  (*estimate*, recalled), would likely be cheaper (not measured), but it would give HMM rather than CM
  boundaries and add a second ncRNA producer with its own score scale. (0-09c's class heuristic counted
  tracrRNA and four xrRNA models as rRNA: 214 hits and 866 CPU-s there; the rRNA models proper have 213 hits
  and 857.1 CPU-s.)

**Recommendation.** Keep Rfam's models for both now. Drop the four eukaryotic rRNA models (a−), which saves
72 CPU-hours per 10 Gbp whole genome (49 intergenic) and loses no locus in these 8 genomes. Before any
switch, measure M2 (how many real tRNAs the filters drop, with `--max`) and, only if the PI first adds the
tools to `environment.yml`, M3 (tRNAscan-SE and barrnap against Rfam on the 8 genomes and the panel).

## 4. Where to scan, and with which flags

From 0-09c §2 (*measured* recall against the whole-genome truth, 1,120 hits in 258 families; CPU *derived*):

| region | CPU-h/10 Gbp (full library) | hit recall | family recall | what is lost |
|---|---|---|---|---|
| whole genome | 1,910 | 1 | 1 | (3 of 100 spike-ins, to the filters) |
| intergenic, margin 0 | 913 | 0.906 | 0.860 | 27 of 364 tRNAs, 7 riboswitches/T-box, 62 sRNAs; 36 families |
| intergenic, margin 50 | 1,241 | 0.935 | 0.922 | 17 tRNAs, 5 riboswitches, 44 sRNAs; 20 families |
| intergenic, margin 200 | 1,231 | 0.983 | 0.965 | 1 tRNA, 18 sRNAs; 9 families |
| intergenic, margin 0, `--notrunc` | 311 | 0.868 | 0.837 | 27 tRNAs, 39 rRNA model hits, 69 sRNAs; 42 families |
| intergenic, margin 200, `--notrunc` | 796 | 0.980 | 0.957 | 1 tRNA, 22 sRNAs; 11 families |

Margin 0 loses tRNAs and riboswitch leaders to short spurious ORFs called across them (14 of its 27 lost
tRNAs are in the archaeon), so a gene caller's error becomes an ncRNA loss. Margin 200 recovers every
riboswitch, T-box, leader and rRNA hit. `--notrunc` stops the re-search of each short record's ends for
truncated hits, which costs more than the main search at margin 0 (0-09c §2). **Proposal: intergenic
regions with a 200-base margin, `--notrunc`.**

**RNA inside coding sequence.** Margin 200 with `--notrunc` loses 11 families entirely: DnaX (the *dnaX*
frameshift element), the *M. tuberculosis* cis-antisense RNAs AS1726, ASdes, ASpks and Ms_AS-5, and
5_ureB_sRNA, MicL, AsdA, AsrC, Erse and Ysr276. All lie inside or across coding sequence. They are a real
class: antisense regulators, recoding signals, and sRNAs that overlap mRNAs. No margin finds them.
Proposal: scan the whole genome for a small **in-CDS class** within the set: those 11 families plus every
model whose name or description says antisense, asRNA, frameshift, recoding or SECIS. In set 2' that is 23
models: RF00382 DnaX, RF01779 AS1726, RF01781 ASdes, RF01782 ASpks, RF01843 neisseria_FSE, RF01988 SECIS_2,
RF01989 SECIS_3, RF02194 HPnc0260, RF02376 antisense_SR1, RF02465 Ms_AS-5, RF02514 5_ureB_sRNA, RF02592 Spy392987, RF02654 MicL, RF02662 BLi_r0872, RF02674 AsdA,
RF02736 Rev11, RF02737 Rev13, RF02738 Rev24, RF02741 Rev44, RF02742 Rev72, RF02746 AsrC, RF02757 Erse,
RF02850 Ysr276. Cost: 55.6 CPU-hours per 10 Gbp whole genome instead of 9.6 intergenic, so set 2' costs
**247.4** (*derived*: 201.4 − 9.6 + 55.6; formula above). With only the 11 lost families whole genome it
costs 222.8. Either way, 17 of margin 200's 23 lost hits come back. The 6 still lost are DicF in
*B. subtilis* (inside a gene), two IS009 hits, one naRNA4 hit, one *E. coli* tRNA and one RtT hit, each cut
at a region boundary. That gives hit recall 1,027 of the 1,033 truth hits of set 2' models (99.4%), every
family found (*derived* from `igr_m200_notrunc.compare.json`). The 87 hits of the four dropped eukaryotic
rRNA models are left out by design; their loci are covered by the bacterial and archaeal models. The class is
defined by description and by what 8 genomes showed; M1 should list which (d) and (b) models hit inside
genes on the sample.

Model set and region change the hits, so both must be part of the producer's name (for example
`infernal-1.1.5-rfam-15.1-prok1-igr200nt`). Today `library_version` would name a subset CM only by its
sha256 prefix.

## 5. Overlapping hits (card 0-08's first question)

**Today** (`src/biobox/l1/ncrna.py`): every hit at or above GA is kept. Two models with identical
coordinates share one `feature` row; any other overlap is a separate feature and tape slot. On the 8
genomes, 1,120 hits become 1,095 features on 971 same-strand loci (*derived* from 0-09c's truth hits).
Of those loci, 56 carry more than one model (205 hits):

- 25 SSU loci, each with 4 SSU models;
- 25 LSU loci, each with 3 LSU models. The 5.8S model also hits 12 of them, at the 5' end of the 23S gene.
  In 6 *S. pyogenes* operons the 5S gene joins the cluster, because the LSU_rRNA_eukarya hit runs 48 bases
  past the 23S gene into the 5S gene;
- 4 tRNA plus RtT loci (*E. coli*);
- 1 Bacteria_large_SRP plus Bacteria_small_SRP locus (*B. subtilis*);
- 1 GlmZ plus GlmY locus (*E. coli*).

On opposite strands, 12 pairs of hits overlap: naRNA4 on both strands of 6 palindromic repeats,
Bacillaceae-1 antisense to the 16S gene 4 times, HPnc0260 across a tRNA, and SdsR across SraC.

**Rfam's rule** is clan competition: among overlapping hits of models in one clan, keep the best. In
Infernal this is `cmscan --clanin <Rfam.clanin> --oclan` (*estimate*, recalled; `cmsearch` has no clan
option, so the producer would apply it after the scan). It cannot be applied without `Rfam.clanin`, which
this project has not fetched.

Options:

1. **Keep all, let readers choose (today).** Nothing is lost. But the tape carries 124 extra slots on
   these 8 genomes (1,095 features for 971 loci), an rRNA operon looks like 9 genes to a neighbour query,
   and every reader must repeat the choice.
2. **Best-scoring hit per overlapping region, regardless of clan.** No new file, but it is wrong for nested
   and adjacent real features. On these genomes, on the same strand it drops the 6 *S. pyogenes* 5S rRNAs
   (chained through the eukaryotic LSU model's overhang), the 4 RtT hits and GlmY. Across strands it also
   drops one partner of each antisense pair (SdsR or SraC, HPnc0260 or the tRNA, Bacillaceae-1 or the 16S).
3. **Fetch `Rfam.clanin` and apply clan competition** (card 0-08 option 1b). Hits on the same contig and
   strand that overlap by at least one base, from models of one clan: the best one (by bit score; ties by
   accession, so it is deterministic) is kept and becomes the feature. The others stay in `l1_ncrna`,
   flagged, off the tape. Models in no clan never compete. Expected on these genomes: one SSU and one LSU
   feature per operon. 5S, tRNA, RtT and the antisense pairs are kept, unless Rfam puts them in one clan
   with their neighbour (clan memberships are unverified here).

**Recommendation: 3**, together with (a−), which removes the 5S chaining at its source. The rule needs two
`l1_ncrna` columns appended before `rS0.1` is signed (AGENTS §3.8): the model's clan (nullable) and whether
the hit won. It also needs one registry statement: only winning hits have a `feature` row, so a losing row's
`feature_uid` names its coordinates (AGENTS §3.3) but no feature. The test that holds it uses constructed
truth (AGENTS §5.6), the E-02 panel extended with an operon. Plant, at recorded coordinates, `cmemit` draws
of SSU_rRNA_bacteria, a tRNA, LSU_rRNA_bacteria and 5S_rRNA in operon order with realistic spacers (5S
starting 50 to 120 bases after the 23S), plus an antisense pair (SdsR on one strand across SraC on the
other). Each draw must be detectable alone (the panel's rule). The test asserts:

- exactly one kept hit of each planted rRNA gene's clan covers at least half of it;
- the 5S, the tRNA and both antisense partners are kept;
- no two kept hits of one clan overlap;
- every losing hit is in `l1_ncrna` with its flag;
- tape slots equal planted genes.

A unit test with a hand-written clan file checks the selection against a brute-force pairwise loop written
separately (the oracle).

## 6. The decision

Tick one per line:

1. **Model rule:** [ ] (a) only, 258 models · [ ] (a) + (b) · [x] **(a−) + (b), 588 models, plus the (d)
   models that hit in M1** · [ ] as before without (b10) · [ ] full library without (c). Optional exclusions:
   [ ] (c7) [ ] (c8) [ ] (c9) (they only matter if (d) is scanned whole).
2. **tRNA and rRNA:** [x] **Rfam's models, without the 4 eukaryotic rRNA models** · [ ] add tRNAscan-SE to
   `environment.yml` (then M3) · [ ] add barrnap to `environment.yml` (then M3).
3. **Region and flags:** [x] **intergenic, margin 200, `--notrunc`, and the in-CDS class (23 models) on
   whole genomes** · [ ] the same with only the 11 lost families whole genome · [ ] margin 0 · [ ] whole
   genome for everything.
4. **Overlaps:** [ ] keep all · [ ] best per region · [x] **clan competition from `Rfam.clanin`, losers
   kept in `l1_ncrna` and off the tape**.
5. **Files for `refdb.yaml` (same release, so no conflict with DEC-12; not fetched; URL form follows the
   `rfam` entry, existence unverified):** [ ] `Rfam.clanin`,
   `https://ftp.ebi.ac.uk/pub/databases/Rfam/15.1/Rfam.clanin` (needed for 4) · [ ] `family.txt.gz`,
   `https://ftp.ebi.ac.uk/pub/databases/Rfam/15.1/database_files/family.txt.gz` (optional: it would replace
   the description rules with Rfam's type field).

Proposed configuration: **247 CPU-hours per 10 Gbp** (*derived*: 247.4), against 1,910 today, 7.7 times less,
and below Pfam's 443. It rises to at most about 425 once M1's long-tail families are added.

**Cards these choices spawn** (each **M**, one per session; tests first, AGENTS §5):

- **E-07, amended: ncRNA model set, regions, shards and `build-stage`.** Depends: E-01, C-12, B-20, this
  card. Files: `src/biobox/l1/ncrna.py`, `src/biobox/l1/rfam_sets/prok-1.tsv`+ (accession and whether
  whole genome; generated from this note's rules and frozen), `src/biobox/stage.py`, `src/biobox/cli.py`
  (`stage ncrna plan|shard|merge`, `build-stage --ncrna-cm`), `tests/r2/test_ncrna_regions.py`+,
  `tests/r2/test_ncrna_shards.py`+, `RUNBOOK.md`. Tests first:
  - regions from a fixture's gene calls against a separately written per-base loop (margin, merging,
    linear ends, origin wrap);
  - hits mapped back to contig coordinates, including across the origin;
  - the subset's hits equal the full fixture library's hits restricted to the subset (oracle: the
    `cmsearch` binary on the full fixture);
  - in-CDS models are scanned whole genome and nowhere twice;
  - sharded output is byte-identical to unsharded.

  Producer name carries set and region. Shards are cut by genome bases: about 202 Mbp per process at 247.4
  CPU-hours per 10 Gbp under 300 CPU-minutes (*derived*: 18,000 / 3,600 / 247.4 × 1e10). Memory: 0-09c
  measured at most 4.8 GiB per `cmsearch --cpu 2` whole genome and 3.8 GiB intergenic with `--notrunc`;
  the card's 12.3 GiB was at `--cpu 8`. Accept: tests pass; on the 8 golden genomes (Pod: yes) the step
  reproduces 0-09c's margin-200 `--notrunc` hits restricted to the set, plus the whole-genome hits of the
  in-CDS class, at a CPU cost within 20% of 247.4 per 10 Gbp; `build-stage` with and without
  `--ncrna-cm` gives two tape versions with byte-identical gene-call tables.
- **E-08, new: clan competition.** Depends: E-01, and the PI adding `Rfam.clanin` to `refdb.yaml` (H).
  Files: `src/biobox/l1/ncrna.py`, `src/biobox/registry/schemas.py` (`l1_ncrna`: append nullable `clan` and
  `clan_winner`), the schema locks by `scripts/update_schema_lock.py`, `tests/r2/test_ncrna_clans.py`+,
  `tests/r3/test_g04b_operon.py`+. Tests first: the unit test and the operon panel of section 5. Accept:
  both pass; on the 8 golden genomes, one kept SSU and one kept LSU feature per operon (25 and 25), with
  every other truth hit kept.
- **0-08a, measurement: the long tail and the exclusions on a corpus sample** (M1 below). Pod: yes. Files:
  `scripts/stress/annotation_cost/rfam_sample.py`+, `planning/measurements/<date>_rfam_sample.md`+.
- **0-08b, measurement: tRNA filter loss, and set 2' as one scan** (M2, M4, M5). Pod: yes. Files: the same
  scripts directory, `planning/measurements/<date>_ncrna_set.md`+.
- **0-08c, only if 2 adds a tool:** the PI's change to `container/environment.yml` (H), then M3 (M, Pod:
  yes).

## 7. What this note could not establish, and the measurement that would

Each fits in about a day on the pod or less, at the thread cap.

- **M1. How many more families hit beyond 8 genomes, and which (d) models to keep.** Scan the full library
  on whole genomes of a 100-genome sample of the 2,497 prokaryotes of the 10 Gbp set, stratified by phylum,
  with archaea and a few hundred phages added. Report:
  - families with hits as genomes are added (8 to 100);
  - the (d) models that hit;
  - any (c) model that hits (a wrong exclusion);
  - hits lying inside genes, by family.

  Cost: about 81 CPU-hours (*derived*: 1,903 × 0.0425; the sample's 425 Mbp is an *estimate* at 4.25 Mbp
  per genome, STATUS.md). That is 2 to 3 hours of wall time at 40 threads in about 20 processes of under 300
  CPU-minutes each (*derived*).
- **M2. Real tRNAs lost to the filters.** RF00005 and RF01852 alone with `--max` on the 8 genomes; count
  loci at or above GA found only with filters off. Cost not measured; *estimate* under 10 CPU-hours (two
  short models).
- **M3. Dedicated tools** (needs the PI's environment change first): tRNAscan-SE and barrnap against the
  Rfam hits on the 8 genomes and the E-02 panel: recall, boundaries, CPU.
- **M4. Set 2' as one scan.** Run the 588 models at margin 200 with `--notrunc`, and the 23 in-CDS models
  whole genome, on the 8 genomes. This checks the derived 247.4 (per-model lines were checked as predictors
  only for whole-genome subsets). Cost about 0.45 CPU-hours (*derived*: 247.4 × 18,106,231 / 1e10).
- **M5. Short contigs.** The same scan on about 50 Mbp of phage and metagenome-like contigs from the
  10 Gbp set, because intergenic cost follows the number of records (0-09c §2). Cost about 1.2 CPU-hours if
  the golden rate holds (*estimate*).
- Not measurable without fetching: exact family types and clans (`family.txt`, `Rfam.clanin`), and whether
  Rfam's clans join any of the nested pairs in section 5.

## Appendix A. The 258 models with a hit in the 8 genomes (rule a)

Hits: truth hits at or above GA in the 8 genomes (*measured*, 0-09c `wg_default`). CPU-s: the model's own
`# CPU time` line summed over `wgA_default` and `wgB_default` (*measured*); share: of all models' 12,403.2
CPU-s. The four eukaryotic rRNA models (RF00002, RF01960, RF02542, RF02543) are dropped by (a-).

| accession | name | description | hits | CPU-s | share |
|---|---|---|---|---|---|
| RF00001 | 5S_rRNA | 5S ribosomal RNA | 26 | 1.64 | 0.013% |
| RF00002 | 5_8S_rRNA | 5.8S ribosomal RNA | 12 | 5.97 | 0.048% |
| RF00005 | tRNA | tRNA | 363 | 4.51 | 0.036% |
| RF00010 | RNaseP_bact_a | Bacterial RNase P class A | 2 | 12.75 | 0.103% |
| RF00011 | RNaseP_bact_b | Bacterial RNase P class B | 2 | 7.14 | 0.058% |
| RF00013 | 6S | 6S / SsrS RNA | 4 | 2.31 | 0.019% |
| RF00014 | DsrA | DsrA RNA | 1 | 1.73 | 0.014% |
| RF00018 | CsrB | CsrB/RsmB RNA family | 1 | 6.66 | 0.054% |
| RF00021 | Spot_42 | Spot 42 RNA | 1 | 2.22 | 0.018% |
| RF00022 | GcvB | GcvB RNA | 1 | 4.85 | 0.039% |
| RF00023 | tmRNA | transfer-messenger RNA | 5 | 45.57 | 0.367% |
| RF00028 | Intron_gpI | Group I catalytic intron | 2 | 24.31 | 0.196% |
| RF00029 | Intron_gpII | Group II catalytic intron | 2 | 0.91 | 0.007% |
| RF00033 | MicF | MicF RNA | 1 | 5.59 | 0.045% |
| RF00034 | RprA | RprA RNA | 1 | 2.43 | 0.020% |
| RF00035 | OxyS | OxyS RNA | 1 | 1.75 | 0.014% |
| RF00039 | DicF | DicF RNA | 2 | 1.67 | 0.013% |
| RF00040 | rne5 | RNase E 5' UTR element | 1 | 8.35 | 0.067% |
| RF00050 | FMN | FMN riboswitch aptamer (RFN element) | 4 | 2.17 | 0.017% |
| RF00057 | RyhB | RyhB RNA | 1 | 1.58 | 0.013% |
| RF00059 | TPP | TPP riboswitch aptamer (THI element) | 11 | 1.64 | 0.013% |
| RF00077 | SraB | SraB RNA | 1 | 4.20 | 0.034% |
| RF00078 | MicA | MicA sRNA | 1 | 2.14 | 0.017% |
| RF00079 | OmrA-B | OmrA-B family | 2 | 1.86 | 0.015% |
| RF00080 | Manganese | Manganese Mn2+ riboswitch aptamer (yybP-ykoY) | 1 | 2.67 | 0.022% |
| RF00081 | ArcZ | ArcZ RNA | 1 | 2.62 | 0.021% |
| RF00082 | SraG | SraG RNA | 1 | 10.87 | 0.088% |
| RF00083 | GlmZ_SraJ | GlmZ RNA activator of glmS mRNA | 1 | 2.55 | 0.021% |
| RF00084 | CsrC | CsrC RNA family | 1 | 5.17 | 0.042% |
| RF00101 | SraC_RyeA | SraC/RyeA RNA | 1 | 2.53 | 0.020% |
| RF00110 | RybB | RybB RNA | 1 | 1.59 | 0.013% |
| RF00111 | SdsR_RyeB | SdsR_RyeB RNA | 1 | 1.67 | 0.013% |
| RF00112 | CyaR_RyeE | CyaR/Rye RNA | 1 | 1.90 | 0.015% |
| RF00113 | SIB_RNA | Short Intergenic Abundant RNA | 5 | 4.86 | 0.039% |
| RF00114 | S15 | Ribosomal S15 leader | 1 | 10.75 | 0.087% |
| RF00115 | McaS | McaS/IsrA RNA | 1 | 1.55 | 0.012% |
| RF00116 | C0465 | C0465 RNA | 1 | 2.98 | 0.024% |
| RF00117 | C0719 | C0719 RNA | 1 | 5.15 | 0.042% |
| RF00118 | rydB | rydB RNA | 1 | 1.05 | 0.008% |
| RF00119 | C0299 | C0299 RNA | 1 | 1.82 | 0.015% |
| RF00121 | MicC | MicC RNA | 1 | 9.89 | 0.080% |
| RF00122 | GadY | GadY | 1 | 1.79 | 0.014% |
| RF00124 | IS102 | IS102 RNA | 1 | 6.08 | 0.049% |
| RF00125 | IS128 | IS128 RNA | 1 | 18.62 | 0.150% |
| RF00126 | ryfA | ryfA RNA | 1 | 11.55 | 0.093% |
| RF00127 | t44 | t44 RNA | 1 | 1.14 | 0.009% |
| RF00128 | GlmY_tke1 | Glm Y RNA activator of glmS mRNA | 2 | 1.81 | 0.015% |
| RF00140 | Alpha_RBS | Alpha operon ribosome binding site | 1 | 2.54 | 0.020% |
| RF00162 | SAM | SAM riboswitch aptamer (S box leader) | 11 | 2.15 | 0.017% |
| RF00167 | Purine | Purine riboswitch aptamer | 6 | 1.72 | 0.014% |
| RF00168 | Lysine | Lysine riboswitch aptamer | 3 | 3.24 | 0.026% |
| RF00169 | Bacteria_small_SRP | Bacterial small signal recognition particle RNA | 4 | 1.60 | 0.013% |
| RF00174 | Cobalamin | Cobalamin riboswitch aptamer | 4 | 2.47 | 0.020% |
| RF00177 | SSU_rRNA_bacteria | Bacterial small subunit ribosomal RNA | 25 | 43.20 | 0.348% |
| RF00230 | T-box | T-box leader | 18 | 6.10 | 0.049% |
| RF00234 | glmS | glmS glucosamine-6-phosphate activated ribozyme aptamer | 1 | 2.19 | 0.018% |
| RF00368 | sroB | sroB RNA | 1 | 3.51 | 0.028% |
| RF00369 | sroC | sroC RNA | 1 | 1.83 | 0.015% |
| RF00370 | sroD | sroD RNA | 1 | 2.23 | 0.018% |
| RF00371 | sroE | sroE RNA | 1 | 1.28 | 0.010% |
| RF00372 | sroH | sroH RNA | 1 | 3.42 | 0.028% |
| RF00373 | RNaseP_arch | Archaeal RNase P | 1 | 6.25 | 0.050% |
| RF00379 | c-di-AMP | c-di-AMP riboswitch aptamer (ydaO/yuaA leader) | 3 | 1.54 | 0.012% |
| RF00380 | Magnesium | M-box riboswitch aptamer (ykoK leader) | 3 | 2.65 | 0.021% |
| RF00382 | DnaX | DnaX ribosomal frameshifting element | 1 | 2.03 | 0.016% |
| RF00391 | RtT | RtT RNA | 13 | 2.59 | 0.021% |
| RF00442 | Guanidine-I | Guanidine-I riboswitch aptamer | 2 | 1.46 | 0.012% |
| RF00504 | Glycine | Glycine riboswitch aptamer | 5 | 1.77 | 0.014% |
| RF00505 | RydC | RydC RNA | 1 | 1.75 | 0.014% |
| RF00506 | Thr_leader | Threonine operon leader | 1 | 9.90 | 0.080% |
| RF00512 | Leu_leader | Leucine operon leader | 1 | 4.97 | 0.040% |
| RF00513 | Trp_leader | Tryptophan operon leader | 1 | 3.89 | 0.031% |
| RF00514 | His_leader | Histidine operon leader | 1 | 2.71 | 0.022% |
| RF00515 | PyrR | PyrR binding site | 5 | 5.71 | 0.046% |
| RF00516 | ylbH | ylbH leader | 1 | 3.72 | 0.030% |
| RF00522 | PreQ1 | PreQ1 (pre queuosine) riboswitch aptamer | 1 | 1.11 | 0.009% |
| RF00534 | SgrS | SgrS RNA | 1 | 4.94 | 0.040% |
| RF00552 | rncO | rncO | 1 | 4.74 | 0.038% |
| RF00555 | L13_leader | Ribosomal protein L13 leader | 2 | 1.02 | 0.008% |
| RF00556 | L19_leader | Ribosomal protein L19 leader | 1 | 1.17 | 0.009% |
| RF00557 | L10_leader | Ribosomal protein L10 leader | 2 | 3.24 | 0.026% |
| RF00558 | L20_leader | Ribosomal protein L20 leader | 2 | 4.73 | 0.038% |
| RF00559 | L21_leader | Ribosomal protein L21 leader | 2 | 1.05 | 0.008% |
| RF00630 | P26 | Pseudomonas sRNA P26 | 1 | 1.74 | 0.014% |
| RF00634 | SAM-IV | SAM (S-adenosyl methionine) riboswitch aptamer | 1 | 4.32 | 0.035% |
| RF01054 | preQ1-II | PreQ1-II (pre queuosine) riboswitch aptamer | 1 | 2.65 | 0.021% |
| RF01055 | Molybdenum | Molybdenum/Tungsten riboswitch aptamer (Moco RNA motif) | 1 | 1.78 | 0.014% |
| RF01056 | Magnesium-II | Magnesium riboswitch aptamer (Magnesium responsive RNA element) | 1 | 2.19 | 0.018% |
| RF01065 | 23S-methyl | 23S methyl RNA motif | 1 | 2.06 | 0.017% |
| RF01066 | 6C | 6C RNA | 1 | 1.09 | 0.009% |
| RF01068 | Guanidine-II | Guanidine-II riboswitch aptamer (mini-ykkC) | 1 | 1.10 | 0.009% |
| RF01125 | sR4 | Small nucleolar RNA sR4 | 1 | 0.60 | 0.005% |
| RF01136 | sR28 | Small nucleolar RNA sR28 | 1 | 0.68 | 0.005% |
| RF01137 | sR21 | Small nucleolar RNA sR21 | 1 | 0.67 | 0.005% |
| RF01140 | sR20 | Small nucleolar RNA sR20 | 1 | 2.55 | 0.021% |
| RF01141 | sR18 | Small nucleolar RNA sR18 | 1 | 0.68 | 0.005% |
| RF01144 | sR17 | Small nucleolar RNA sR17 | 1 | 0.68 | 0.005% |
| RF01145 | sR14 | Small nucleolar RNA sR14 | 1 | 0.67 | 0.005% |
| RF01304 | sR5 | Small nucleolar RNA sR5 | 2 | 3.96 | 0.032% |
| RF01312 | sR9 | Small nucleolar RNA sR9 | 1 | 1.07 | 0.009% |
| RF01317 | CRISPR-DR4 | CRISPR RNA direct repeat element | 1 | 0.75 | 0.006% |
| RF01335 | CRISPR-DR22 | CRISPR RNA direct repeat element | 6 | 1.14 | 0.009% |
| RF01338 | CRISPR-DR25 | CRISPR RNA direct repeat element | 16 | 1.06 | 0.009% |
| RF01400 | istR | istR Hfq binding RNA | 1 | 3.39 | 0.027% |
| RF01401 | rseX | rseX Hfq binding RNA | 1 | 3.81 | 0.031% |
| RF01407 | STnc560 | STnc560 Hfq binding RNA | 1 | 9.50 | 0.077% |
| RF01408 | sraL | sraL Hfq binding RNA | 1 | 2.43 | 0.020% |
| RF01410 | BsrC | BsrC | 2 | 1.54 | 0.012% |
| RF01411 | BsrF | BsrF | 1 | 2.15 | 0.017% |
| RF01412 | BsrG | BsrG | 3 | 13.52 | 0.109% |
| RF01470 | rli38 | Listeria sRNA rli38 | 1 | 9.40 | 0.076% |
| RF01492 | rli28 | Listeria sRNA rli28 | 1 | 10.60 | 0.085% |
| RF01517 | iscRS | iscR stability element | 1 | 1.54 | 0.012% |
| RF01688 | Actino-pnp | Actino-pnp RNA | 1 | 1.49 | 0.012% |
| RF01690 | Bacillaceae-1 | Bacillaceae-1 RNA | 9 | 1.52 | 0.012% |
| RF01691 | Bacillus-plasmid | Bacillus-plasmid RNA | 1 | 1.04 | 0.008% |
| RF01707 | JUMPstart | JUMPstart RNA | 3 | 0.91 | 0.007% |
| RF01708 | L17DE | L17 ribosomal protein downstream element | 1 | 1.29 | 0.010% |
| RF01709 | Lacto-rpoB | Lacto-rpoB RNA | 1 | 1.25 | 0.010% |
| RF01732 | asd | MarS sRNA | 1 | 2.11 | 0.017% |
| RF01735 | epsC | epsC RNA | 1 | 1.40 | 0.011% |
| RF01746 | mraW | mraW RNA | 1 | 1.67 | 0.013% |
| RF01748 | nuoG | nuoG RNA | 1 | 1.04 | 0.008% |
| RF01749 | pan | pan motif | 1 | 3.48 | 0.028% |
| RF01764 | azaaromatic | azaaromatic riboswitch aptamer (yjdF RNA) | 1 | 1.29 | 0.010% |
| RF01766 | cspA | cspA thermoregulator | 7 | 23.00 | 0.185% |
| RF01769 | greA | Enterobacteria greA leader | 1 | 3.54 | 0.029% |
| RF01770 | rimP | Gammaprotebacteria rimP leader | 1 | 9.90 | 0.080% |
| RF01771 | rnk_leader | Enterobacteria rnk leader | 1 | 2.95 | 0.024% |
| RF01776 | RatA | RNA anti-toxin A | 1 | 2.19 | 0.018% |
| RF01779 | AS1726 | AS1726 sRNA | 1 | 11.33 | 0.091% |
| RF01780 | AS1890 | AS1890 sRNA | 1 | 3.56 | 0.029% |
| RF01781 | ASdes | ASdes TB sRNA | 2 | 3.93 | 0.032% |
| RF01782 | ASpks | ASpks TB sRNA | 5 | 4.42 | 0.036% |
| RF01783 | b55 | Mycobacterium B11 | 1 | 1.24 | 0.010% |
| RF01791 | F6 | F6 TB sRNA | 1 | 1.62 | 0.013% |
| RF01794 | sok | sok antitoxin (CssrC) | 5 | 3.17 | 0.026% |
| RF01796 | frnS | Fumarate/nitrate reductase regulator sRNA | 1 | 4.59 | 0.037% |
| RF01798 | g2 | G2 | 1 | 1.03 | 0.008% |
| RF01804 | Lambda_thermo | Lambda phage CIII thermoregulator element | 1 | 1.62 | 0.013% |
| RF01809 | symR | SymR antitoxin | 1 | 2.20 | 0.018% |
| RF01813 | rdlD | rdlD antitoxin | 4 | 1.25 | 0.010% |
| RF01820 | RsaE | RNA Staph. aureus E (RoxS) | 1 | 3.59 | 0.029% |
| RF01829 | sR6 | sR6 snoRNA | 1 | 1.41 | 0.011% |
| RF01830 | StyR-44 | Salmonella enterica Typhi npcRNA 44 | 7 | 2.59 | 0.021% |
| RF01832 | ROSE_2 | Repression of heat shock gene expression (ROSE) element | 1 | 1.69 | 0.014% |
| RF01852 | tRNA-Sec | Selenocysteine transfer RNA | 1 | 4.48 | 0.036% |
| RF01854 | Bacteria_large_SRP | Bacterial large signal recognition particle RNA | 1 | 3.65 | 0.029% |
| RF01857 | Archaea_SRP | Archaeal signal recognition particle RNA | 1 | 4.68 | 0.038% |
| RF01859 | Phe_leader | Phenylalanine leader peptide | 2 | 16.63 | 0.134% |
| RF01959 | SSU_rRNA_archaea | Archaeal small subunit ribosomal RNA | 25 | 48.88 | 0.394% |
| RF01960 | SSU_rRNA_eukarya | Eukaryotic small subunit ribosomal RNA | 25 | 102.00 | 0.822% |
| RF01988 | SECIS_2 | Selenocysteine insertion sequence 2 | 1 | 2.12 | 0.017% |
| RF01989 | SECIS_3 | Selenocysteine insertion sequence 3 | 2 | 1.34 | 0.011% |
| RF02003 | group-II-D1D4-4 | Group II catalytic intron D1-D4-4 | 2 | 1.64 | 0.013% |
| RF02029 | sraA | sraA | 1 | 4.58 | 0.037% |
| RF02030 | tp2 | tp2 | 1 | 2.53 | 0.020% |
| RF02031 | tpke11 | tpke11 | 1 | 1.26 | 0.010% |
| RF02051 | STnc450 | Enterobacterial sRNA STnc450 | 1 | 1.21 | 0.010% |
| RF02052 | STnc630 | Enterobacterial sRNA STnc630 | 1 | 6.93 | 0.056% |
| RF02053 | STnc430 | Enterobacterial sRNA STnc430 | 1 | 12.29 | 0.099% |
| RF02057 | STnc40 | Salmonella enterica sRNA STnc40 | 1 | 2.66 | 0.021% |
| RF02060 | STnc410 | Enterobacterial sRNA STnc410 | 1 | 11.68 | 0.094% |
| RF02064 | STnc370 | Enterobacterial sRNA STnc370 | 1 | 6.10 | 0.049% |
| RF02068 | STnc480 | Enterobacterial sRNA STnc480 | 1 | 2.13 | 0.017% |
| RF02069 | STnc70 | Enterobacterial sRNA STnc70 | 1 | 1.53 | 0.012% |
| RF02074 | STnc240 | Enterobacterial sRNA STnc240 | 1 | 3.65 | 0.029% |
| RF02079 | STnc180 | Enterobacterial sRNA STnc180 | 1 | 13.46 | 0.109% |
| RF02081 | STnc550 | Enterobacterial sRNA STnc550 | 1 | 17.33 | 0.140% |
| RF02082 | STnc540 | Enterobacterial sRNA STnc540 | 1 | 5.34 | 0.043% |
| RF02083 | OrzO-P | OrzO-P RNA antitoxin family | 1 | 3.58 | 0.029% |
| RF02084 | STnc130 | Enterobacteria sRNA STnc130 | 1 | 8.28 | 0.067% |
| RF02099 | rivX | rivX sRNA | 1 | 5.39 | 0.043% |
| RF02111 | IS009 | IS009 | 5 | 14.06 | 0.113% |
| RF02194 | HPnc0260 | Bacterial antisense RNA HPnc0260 | 1 | 20.88 | 0.168% |
| RF02273 | FsrA | FsrA | 1 | 1.77 | 0.014% |
| RF02341 | ncrMT1302 | Mycobacterium sRNA ncrMT1302 | 1 | 1.72 | 0.014% |
| RF02348 | tracrRNA | Trans-activating crRNA | 1 | 2.27 | 0.018% |
| RF02370 | Trp_leader_2 | Bacillus tryptophan operon leader | 1 | 3.54 | 0.029% |
| RF02371 | PyrG_leader | PyrG leader | 1 | 4.53 | 0.037% |
| RF02372 | PyrC_leader | PyrC leader | 1 | 1.40 | 0.011% |
| RF02373 | PyrD_leader | PyrD leader | 1 | 1.44 | 0.012% |
| RF02376 | antisense_SR1 | SR1 sRNA | 1 | 8.39 | 0.068% |
| RF02377 | SurA | SurA sRNA | 1 | 8.65 | 0.070% |
| RF02378 | SurC | SurC sRNA | 1 | 4.96 | 0.040% |
| RF02384 | FasX | FasX small RNA | 1 | 11.98 | 0.097% |
| RF02401 | ClpQY_promoter | ClpQY promoter | 1 | 0.95 | 0.008% |
| RF02419 | Spd-sr37 | Streptococcus sRNA Spd-sr37 | 1 | 1.79 | 0.014% |
| RF02449 | ncr1015 | Bacillus sRNA ncr1015 | 1 | 4.98 | 0.040% |
| RF02450 | ncr1175 | Bacillus sRNA ncr1175 | 1 | 3.87 | 0.031% |
| RF02451 | ncr1241 | Bacillus sRNA ncr1241 | 1 | 1.13 | 0.009% |
| RF02452 | ncr1575 | Bacillus sRNA ncr1575 | 1 | 5.13 | 0.041% |
| RF02453 | ncr952 | Bacillus sRNA ncr952 | 1 | 6.00 | 0.048% |
| RF02454 | ncr982 | Bacillus sRNA ncr982 | 1 | 1.78 | 0.014% |
| RF02465 | Ms_AS-5 | Mycobacterium sRNA Ms_AS-5 | 1 | 1.40 | 0.011% |
| RF02495 | ohsC_RNA | Oppression of Hydrophobic ORF by sRNA | 1 | 2.48 | 0.020% |
| RF02514 | 5_ureB_sRNA | 5' ureB small RNA | 2 | 41.42 | 0.334% |
| RF02524 | sagA | Streptococcus sRNA sagA | 1 | 2.06 | 0.017% |
| RF02525 | SSRC30 | Streptococcus sRNA SSRC30 | 1 | 3.22 | 0.026% |
| RF02526 | SSRC34_1 | Streptococcus sRNA SSRC34 | 1 | 7.95 | 0.064% |
| RF02527 | SSRC38 | Streptococcus sRNA SSRC38 | 1 | 10.55 | 0.085% |
| RF02528 | SSRC41 | Streptococcus sRNA SSRC41 | 1 | 4.88 | 0.039% |
| RF02529 | SSRC8_1 | Streptococcus sRNA SSRC8 | 1 | 5.28 | 0.043% |
| RF02540 | LSU_rRNA_archaea | Archaeal large subunit ribosomal RNA | 25 | 153.05 | 1.234% |
| RF02541 | LSU_rRNA_bacteria | Bacterial large subunit ribosomal RNA | 25 | 142.49 | 1.149% |
| RF02542 | SSU_rRNA_microsporidia | Microsporidia small subunit ribosomal RNA | 25 | 44.89 | 0.362% |
| RF02543 | LSU_rRNA_eukarya | Eukaryotic large subunit ribosomal RNA | 25 | 314.99 | 2.540% |
| RF02564 | naRNA4 | Nucleoid-associated noncoding RNA 4 (CssrE) | 96 | 16.33 | 0.132% |
| RF02566 | Ms1 | Mycobacterium smegmatis small RNA 1 | 1 | 7.40 | 0.060% |
| RF02571 | Mcr7 | mcr7 sRNA | 1 | 14.69 | 0.118% |
| RF02589 | Spy779816 | S. pyogenes small RNA 779816 | 3 | 1.63 | 0.013% |
| RF02590 | Spy1186876 | S. pyogenes small RNA 1186876 | 1 | 13.29 | 0.107% |
| RF02591 | Spy1786666 | S. pyogenes small RNA 1786666 | 1 | 1.50 | 0.012% |
| RF02592 | Spy392987 | S. pyogenes antisense RNA 392987 | 1 | 1.09 | 0.009% |
| RF02616 | SSRC8_2 | Streptococcus sRNA 8 | 5 | 5.56 | 0.045% |
| RF02617 | SSRC10 | Streptococcus sRNA 10 | 1 | 2.18 | 0.018% |
| RF02618 | SSRC34_2 | Streptococcus sRNA 34 | 1 | 10.58 | 0.085% |
| RF02640 | MOSES4 | S. pyogenes small RNA MOSES4 | 1 | 2.64 | 0.021% |
| RF02641 | Spy490483c | S. pyogenes small RNA Spy490483c | 1 | 2.13 | 0.017% |
| RF02642 | Spy491311c | S. pyogenes small RNA Spy491311c | 1 | 1.66 | 0.013% |
| RF02643 | Spy491738 | S. pyogenes small RNA Spy491738 | 1 | 1.35 | 0.011% |
| RF02644 | Spy490380c | S. pyogenes small RNA Spy490380c | 1 | 6.29 | 0.051% |
| RF02654 | MicL | MicL sRNA | 1 | 11.01 | 0.089% |
| RF02659 | ncRv12659 | ncRv12659 sRNA | 1 | 1.91 | 0.015% |
| RF02674 | AsdA | antisense RNA of dnaA mRNA | 1 | 2.80 | 0.023% |
| RF02694 | RalA | RalR antitoxin | 1 | 5.16 | 0.042% |
| RF02700 | HtrA_thermometer | HtrA 5' UTR thermometer | 1 | 1.16 | 0.009% |
| RF02744 | Rev39_5p_UTR | Soft rot Enterobacteriaceae Rev 39 5'UTR | 1 | 6.06 | 0.049% |
| RF02746 | AsrC | antisense RNA of rseC mRNA | 1 | 41.51 | 0.335% |
| RF02750 | ES003 | ES003 sRNA | 1 | 1.62 | 0.013% |
| RF02751 | ES036 | ES036 (CssrF) sRNA | 27 | 1.44 | 0.012% |
| RF02752 | ES056 | ES056 sRNA | 1 | 2.03 | 0.016% |
| RF02753 | ES173 | ES173 sRNA | 1 | 1.16 | 0.009% |
| RF02754 | ES205 | ES205 sRNA | 1 | 1.71 | 0.014% |
| RF02755 | ES222 | ES222 sRNA | 1 | 1.61 | 0.013% |
| RF02756 | ES239 | ES239 sRNA | 1 | 3.74 | 0.030% |
| RF02757 | Erse | Erse small RNA | 1 | 27.70 | 0.223% |
| RF02768 | Ysr155_RyfD | Yersinia sRNA 155(RyfD) | 1 | 1.87 | 0.015% |
| RF02773 | TrxA_thermometer | TrxA 5' UTR thermometer | 1 | 1.53 | 0.012% |
| RF02782 | CpoB_ybgF_thermometer | CpoB/ybgF ICR thermometer | 1 | 2.26 | 0.018% |
| RF02797 | PssrA | Plasmid-encoded Shigella sRNA A | 1 | 1.66 | 0.013% |
| RF02798 | CssrA | Chromosome-encoded Shigella sRNA A | 1 | 1.34 | 0.011% |
| RF02822 | Srn266 | Streptococcus RNA 266 | 2 | 5.23 | 0.042% |
| RF02850 | Ysr276 | Yersinia sRNA 276 | 1 | 13.80 | 0.111% |
| RF02855 | Ysr251 | Yersinia sRNA 251 | 1 | 5.05 | 0.041% |
| RF02886 | npcTB_6715 | Mycobacterium sRNA 6715 | 1 | 24.92 | 0.201% |
| RF02891 | AgrB | ArsR-gov region gene B | 2 | 2.31 | 0.019% |
| RF02892 | SR6_antitoxin | Bacillus SR6 antitoxin | 1 | 2.01 | 0.016% |
| RF02893 | RpsF_leader | RpsF leader | 1 | 1.41 | 0.011% |
| RF02913 | pemK | pemK RNA | 1 | 0.94 | 0.008% |
| RF02929 | algC | algC RNA | 1 | 1.13 | 0.009% |
| RF03039 | Peptidase-S11 | Peptidase-S11 RNA | 1 | 2.10 | 0.017% |
| RF03075 | DUF3800-VIII | DUF3800-VIII RNA | 1 | 0.84 | 0.007% |
| RF03131 | L25-Gammaproteobacteria | L25-Gammaproteobacteria ribosomal protein leader | 1 | 0.86 | 0.007% |
| RF03136 | S10-Clostridia | S10-Clostridia ribosomal protein leader | 2 | 1.61 | 0.013% |
| RF03156 | L31-Firmicutes | L31-Firmicutes ribosomal protein leader | 1 | 0.74 | 0.006% |
| RF03157 | L31-Gammaproteobacteria | L31-Gammaproteobacteria ribosomal protein leader | 1 | 1.00 | 0.008% |
| RF04312 | trpB | trpB ncRNA motif | 1 | 2.97 | 0.024% |

## Appendix B. Models added by rule (b) (no hit in the 8 genomes)

(b1) and (b2) add nothing: every rRNA and tRNA model is already in (a).

- (b3) RNase P, bacterial and archaeal, 1: RF02357 RNaseP-T
- (b4) tmRNA, 4: RF01849 alpha_tmRNA, RF01850 beta_tmRNA, RF01851 cyano_tmRNA, RF04321 tmRNA_permuted
- (b5) SRP, bacterial and archaeal, 1: RF04183 Bacteroidales_small_SRP
- (b6) 6S, 2: RF01685 6S-Flavo, RF01693 Bacteroidales-1
- (b7) CRISPR direct repeats, tracrRNA, 61: RF01315 CRISPR-DR2, RF01316 CRISPR-DR3, RF01318 CRISPR-DR5, RF01319 CRISPR-DR6, RF01320 CRISPR-DR7, RF01321 CRISPR-DR8, RF01322 CRISPR-DR9, RF01323 CRISPR-DR10, RF01324 CRISPR-DR11, RF01325 CRISPR-DR12, RF01326 CRISPR-DR13, RF01327 CRISPR-DR14, RF01328 CRISPR-DR17, RF01329 CRISPR-DR15, RF01330 CRISPR-DR16, RF01331 CRISPR-DR18, RF01332 CRISPR-DR19, RF01333 CRISPR-DR20, RF01334 CRISPR-DR21, RF01336 CRISPR-DR23, RF01337 CRISPR-DR24, RF01339 CRISPR-DR27, RF01340 CRISPR-DR29, RF01341 CRISPR-DR30, RF01342 CRISPR-DR32, RF01343 CRISPR-DR33, RF01344 CRISPR-DR34, RF01345 CRISPR-DR35, RF01346 CRISPR-DR36, RF01347 CRISPR-DR37, RF01348 CRISPR-DR38, RF01349 CRISPR-DR40, RF01350 CRISPR-DR41, RF01351 CRISPR-DR42, RF01352 CRISPR-DR43, RF01353 CRISPR-DR44, RF01354 CRISPR-DR45, RF01355 CRISPR-DR26, RF01356 CRISPR-DR28, RF01357 CRISPR-DR31, RF01358 CRISPR-DR39, RF01359 CRISPR-DR46, RF01360 CRISPR-DR47, RF01361 CRISPR-DR48, RF01362 CRISPR-DR49, RF01363 CRISPR-DR50, RF01364 CRISPR-DR51, RF01365 CRISPR-DR52, RF01366 CRISPR-DR53, RF01367 CRISPR-DR54, RF01368 CRISPR-DR55, RF01369 CRISPR-DR56, RF01370 CRISPR-DR57, RF01371 CRISPR-DR58, RF01373 CRISPR-DR60, RF01374 CRISPR-DR61, RF01375 CRISPR-DR62, RF01376 CRISPR-DR63, RF01377 CRISPR-DR64, RF01378 CRISPR-DR65, RF01379 CRISPR-DR66
- (b8) group I and II introns, 6: RF01998 group-II-D1D4-1, RF01999 group-II-D1D4-2, RF02001 group-II-D1D4-3, RF02004 group-II-D1D4-5, RF02005 group-II-D1D4-6, RF02012 group-II-D1D4-7
- (b9) riboswitches, T-box, 33: RF00521 SAM_alpha, RF01051 c-di-GMP-I, RF01057 SAH, RF01482 AdoCbl, RF01510 2dG-I, RF01689 AdoCbl-II, RF01704 Glutamine-II, RF01725 SAM-I-IV, RF01726 SAM-II, RF01727 SAM-SAH, RF01734 Fluoride, RF01739 Glutamine, RF01750 ZMP-ZTP, RF01763 Guanidine-III, RF01767 SAM-III, RF01786 c-di-GMP-II, RF01826 SAM_V, RF01831 THF, RF02680 PreQ1-III, RF02683 NiCo, RF02885 SAM_VI, RF02912 AAC_AAD_leader, RF02977 THF-II, RF03013 nadA, RF03054 Xanthine, RF03071 Sodium, RF03164 NiCo-AnGGAG, RF03165 2dG-II, RF03166 glycine-GGAnGA, RF03167 c-di-GMP-I-GGC, RF03168 c-di-GMP-I-UAU, RF03169 c-di-GMP-II-GAG, RF03170 c-di-GMP-II-GCG
- (b10) sRNAs named for a model prokaryote, 188: RF00615 LhrA, RF00616 LhrC, RF00623 P1, RF00624 P9, RF00625 P11, RF00627 P15, RF00629 P24, RF01456 VrrA, RF01457 rli22, RF01458 rli23, RF01459 rliE, RF01460 rliH, RF01461 rli24, RF01462 rli26, RF01463 rli27, RF01464 rliA, RF01465 rli31, RF01466 rli34, RF01467 rli36, RF01468 rli32, RF01469 rli33, RF01471 rliB, RF01472 rli40, RF01473 rli41, RF01475 rli45, RF01476 rliF, RF01477 rli43, RF01478 rli47, RF01479 rli48, RF01480 rli52, RF01481 rli53, RF01483 rli56, RF01484 rli59, RF01485 rli61, RF01486 rli62, RF01487 rliI, RF01488 rli49, RF01489 sbrA, RF01490 rli51, RF01491 rli54, RF01493 rli37, RF01494 rliD, RF01519 CC0196, RF01520 CC0734, RF01521 CC1840, RF01527 CrfA, RF01528 CC3513, RF01530 CC3664, RF01665 P13, RF01668 P10, RF01669 P14, RF01670 P17, RF01671 P18, RF01674 P27, RF01675 CrcZ, RF01676 P31, RF01677 P35, RF01678 P37, RF01679 P36, RF01680 P5, RF01681 P4, RF01682 P8, RF01683 P6, RF01699 Clostridiales-1, RF01772 rnk_pseudo, RF01773 rpsL_pseudo, RF01808 MicX, RF01843 neisseria_FSE, RF01867 CC2171, RF02048 STnc30, RF02049 STnc460, RF02050 STnc470, RF02054 STnc420, RF02055 STnc380, RF02056 STnc390, RF02059 STnc50, RF02062 STnc361, RF02063 STnc350, RF02065 STnc340, RF02066 STnc320, RF02067 STnc310, RF02070 STnc300, RF02071 STnc280, RF02072 STnc590, RF02073 STnc260, RF02075 STnc230, RF02077 STnc220, RF02078 STnc210, RF02080 STnc170, RF02088 STnc510, RF02270 nse_sRNA, RF02385 sau-13, RF02386 sau-19, RF02387 sau-27, RF02388 sau-30, RF02389 sau-31, RF02390 sau-41, RF02392 sau-53, RF02393 sau-59, RF02394 sau-63, RF02395 sau-66, RF02396 sau-5949, RF02397 sau-5971, RF02398 sau-6072, RF02404 P33, RF02405 P34, RF02414 rli60, RF02415 rliG, RF02418 Spd-sr07, RF02425 SpF01_sRNA, RF02426 SpF03_sRNA, RF02427 SpF10_sRNA, RF02428 SpF11_sRNA, RF02429 SpF14_sRNA, RF02430 SpF19_sRNA, RF02431 SpF22_sRNA, RF02432 SpF25_sRNA, RF02433 SpF36_sRNA, RF02434 SpF39_sRNA, RF02435 SpF41_sRNA, RF02436 SpF43_sRNA, RF02437 SpF44_sRNA, RF02438 SpF51_sRNA, RF02439 SpF56_sRNA, RF02440 SpF59_sRNA, RF02441 SpF61_sRNA, RF02442 SpF66_sRNA, RF02443 SpR08_sRNA, RF02444 SpR10_sRNA, RF02445 SpR14_sRNA, RF02446 SpR18_sRNA, RF02447 SpR19_sRNA, RF02448 SpR20_sRNA, RF02463 Ms_AS-1, RF02467 Ms_IGR-2, RF02468 Ms_IGR-4, RF02470 Ms_IGR-8, RF02537 TarA, RF02538 TarB, RF02567 Va-907, RF02568 UptR, RF02577 tsr24, RF02578 tsr25, RF02579 tsr26, RF02580 tsr31, RF02581 tsr32, RF02582 tsr33, RF02583 Teg23, RF02627 Ssr1, RF02652 StyR-3, RF02653 StyR-143, RF02662 BLi_r0872, RF02676 Esr41, RF02736 Rev11, RF02737 Rev13, RF02738 Rev24, RF02739 Rev41, RF02740 Fwd6_3p_UTR, RF02741 Rev44, RF02742 Rev72, RF02745 Rev42_5p_UTR, RF02799 CssrB, RF02813 PA5194_thermometer, RF02818 V_AS5, RF02819 V_AS7, RF02820 V_AS9, RF02821 V_IGR5, RF02834 VqmR, RF02861 IR_84, RF02877 NmsRb, RF02887 mgtC_leader, RF02888 BtsR1, RF02889 Sr006, RF02894 S35, RF02895 S414, RF02896 S774, RF02897 S808, RF02908 rli30, RF02909 rli117, RF02937 Clostridiales-2, RF02942 Clostridiales-3, RF03091 Clostridium-PBP, RF03112 Staphylococcus-1, RF03140 S4-Clostridia, RF03531 n00280, RF03533 SQ1002_RNA, RF03534 sRNA71, RF03535 TEG147
- (b11) archaeal box C/D sRNAs (described as snoRNA), 38: RF00095 snoPyro_CD, RF01119 sR32, RF01120 sR33, RF01121 sR38, RF01122 sR39, RF01123 sR35, RF01124 sR36, RF01126 sR41, RF01127 sR42, RF01128 sR43, RF01129 sR44, RF01130 sR46, RF01131 sR47, RF01132 sR48, RF01133 sR3, RF01134 sR30, RF01135 sR24, RF01138 sR23, RF01139 sR2, RF01142 sR19, RF01143 sR16, RF01146 sR15, RF01147 sR12, RF01149 sR10, RF01150 sR11, RF01152 sR1, RF01273 sR34, RF01274 sR45, RF01275 sR22, RF01276 sR53, RF01297 sR40, RF01303 sR49, RF01305 sR51, RF01306 sR52, RF01307 sR55, RF01308 sR58, RF01309 sR60, RF01310 sR7
