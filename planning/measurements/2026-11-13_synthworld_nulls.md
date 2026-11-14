# synthworld null worlds (§4.1) on `micro`

Card F-02, session 9 (branch `task/F-02-9`). Code: `src/biobox/synthworld/nulls.py` (its docstring states what each
null world keeps and destroys and the decisions below), the null hooks of `src/biobox/synthworld/emit.py`, preset
section `nulls` of `presets/micro.yaml`. Tests: `tests/r1/test_synthworld_nulls_units.py` (null presets leaf by leaf,
the N1 matcher against a brute-force eligibility loop on constructed families, N2 seeds, world ids),
`tests/r2/test_synthworld_nulls.py` (one fixture builds the signal world and seven null worlds in child processes).
Small tables beside this note in `2026-11-13_synthworld_nulls/`.

Host and commands for every *measured* number: wilkcomp02, 2026-11-13, default test image `bio-box-test_64f1bf0.sif`,
`micro` default seed 20261112 (and 20261113 for the second seed), the account at about 38 threads of its 56-thread
night cap before the run. **Run R2**: `sh scripts/pod_pytest.sh F-02-9 -s tests/r1/test_synthworld_nulls_units.py
tests/r2/test_synthworld_nulls.py` (20:47, 3 failures in the tests' own oracles, fixed; 20:55, 25 passed). **Dev**: the same worlds built once into a scratch directory by a
development test (not committed) at 20:35, four at once, for the first look that set the tests' design.

## What was built

| world | how | shares with the signal world |
|---|---|---|
| N0-full | `nulls.build("N0-full", micro, seed)`: preset `modules.count` 0, `profile_only.count` 0, `decoys.partner_modules` [], `fusions.nogd` 0; S1-S11 by the signal world's stage code | species tree, taxonomy, lineage tables, phage set and every per-genome draw (exact, below); no stage directory (ZOMBI2's gene-level draws change with the groups) |
| N0-strict | as N0-full, and `vehicles.cargo` [], `crispr.arrays` 0 | as N0-full |
| N1 | `nulls.build_n1(signal, out)`: S1-S10 of the signal world hard-linked, the relabel drawn, S11 written | S1-S10 byte for byte (every file), the emitted FASTA, `ingest.tsv`, `genomes.json`, `truth.gff`, and every truth table but `truth_module`, `truth_gene.module_id`, `truth_fusion.in_module` |
| N2 | `nulls.n2_seeds(seed, 20)`: 20 u64 seeds from stream `N2/seeds`; each an N0 world | not built in the tests (20 x about 3 min) |

## Determinism and the signal world

| quantity | value | label |
|---|---|---|
| signal world, every stage | backbone `025fdd0d…`, sequences `34ebce13…`, replicons `55bda971…`, mobile `0bd7f0db…`, fusions `db09abf3…`, assembly `2f425158…`, sources `770007bd…`, emit `a2f0e448…`, world `f0262027e4ccdcbc72cd7d90f23fbb5c`: the session-8 digests, unchanged | measured, run R2 |
| N0-full, default seed | world `c2621dba6b989d3c7d94932b64fbf4b6`, byte-identical in two child processes (PYTHONHASHSEED 919 and 1 thread; 0 and 4 threads) | measured, run R2 and dev |
| N0-strict, default seed | world `856547d47e316eb18bc6351943042e10`, byte-identical in the two processes | measured, run R2 and dev |
| N1, default seed | world `4d6648c73d39b91970cb693a8a178978`, byte-identical in the two processes | measured, run R2 and dev |
| second seed 20261113 | N0-full `3c9123841cf4a755b238b216e8562028`, N0-strict `787161057efcfc7271bec00b3e8cb45f`; N1 with null seed 20261113 `b705a3c7381f5cdd6d001dbd8e66bcf0`; all seven world digests distinct | measured, run R2 |

## N0: kept and destroyed

Kept exactly (both variants, default seed; every comparison file by file or row by row against the signal world):
the species tree (`zombi2/species_complete.nwk`), taxonomy and taxa, the S4 and S5 lineage tables, burst tips and
ZOMBI2 tool seeds, the whole phage set (`mobile/phages/`, `phages.fna`, `phage_genes.tsv`), the 16 prophage instances
(phage, clade node, genome, strand), each genome's genetic code, assembly class, N50 and completeness targets, the 200
planned call errors (id, kind, genome), sources (label, kind, accession, BioProject), contamination host-donor pairs,
duplicate roles, and every label of `truth_genome`; the background tight rule gives the 199 (N0-full) and 179 (N0-strict)
background co-directional family pairs met in both worlds the same class in both. The recorded presets differ from the
signal's in exactly the leaves §4.1 names (plus the `null` marker).

Kept in distribution (test in the R2 docstring of `test_n0_keeps_the_background_process_in_distribution`):

| quantity | signal | N0-full | N0-strict | test | label |
|---|---|---|---|---|---|
| root + originated genes per genome, mean | 1,318.0 | 1,273.8 | 1,256.1 | geometric mean of per-genome ratio in [0.9, 1.1]: 0.958, 0.943 | measured |
| same, Pearson r with the signal over 46 genomes | - | 0.976 | 0.985 | >= 0.9 | measured |
| root families' genomes-present, KS D against the signal | - | 0.0487 | 0.0336 | < 0.0728 (alpha 0.01) | measured |
| + strand share of genes | 0.529 | 0.511 | 0.537 | in [0.4, 0.6] | measured |
| tight share of distinct background co-directional family pairs | - | 0.5497 (8,478 pairs) | 0.5499 (8,291) | within 4 binomial s.e. of p_tight 0.55 | measured |
| genome GC, |difference| from the same genome in the signal | - | mean 0.0013, max 0.0045 | mean 0.0016, max 0.0034 | max <= 0.01 | measured |

**Not kept: total genome size.** All genes per genome average 1,937 in the signal world, 1,453 in N0-full and 1,388 in
N0-strict (measured, dev). The signal world's planted families hold 380.7 genes per genome on average, profile-only
families 29.6 and decoys 6.7 (measured, dev): 21.5% of its genes (*derived*: 417.0 / 1,937). Module members carry
loss 0.03 against the genome's 0.24 and transfer up to 0.30 against 0.10, so they grow to many copies; N0 removes them
by definition. At the second seed (crown age 8.3 against about 4 at the default seed) the signal world averages 701
genes per genome and N0-full about 210 (measured: dev build of the signal at 20261113, 105.6 s; run R2's N0 worlds),
so planted families are about 70% of that signal world's genes (*derived*) and both plasmid types have died out. For
the PI: the `micro` preset is fitted at the default seed only; N2's 20 seeds and the nightly rotating seed (§5.1) will
meet worlds of a few hundred genes per genome dominated by planted families unless the module rates or the genome
turnover are refitted (a change to the signal world, outside this card's session).

Destroyed, measured by brute force on the emitted contigs (pairs within one vehicle, plasmid or Cas locus left out in
both worlds; tau the 99.9th percentile of the signal world's background pairs; one-sided exact conditional binomial
test of N0's rate against the signal background's):

| measure | tau | signal planted | signal background | N0-full | N0-strict | label |
|---|---|---|---|---|---|---|
| genomes in which a pair is adjacent | 31 | 61 / 111 | 18 / 21,344 | 15 / 16,764 (p 0.50) | 26 / 16,205 (p 0.024) | measured, run R2 |
| distinct flank contexts | 20 | 37 / 111 | 18 / 21,344 | 15 / 16,764 (p 0.50) | 21 / 16,205 (p 0.12) | measured, run R2 |

Truth: no planted, profile-only or NOGD rows, no decoy family or trap, no `group` junction. N0-strict: no vehicle
or Cas family, no array, no spacer row, `truth_module` empty; plasmids and the 4 prophages kept.

§4.3 positive control (naive hypergeometric co-occurrence over the S8 genomes, no conditioning): share of pairs at
p < 1e-3 is 0.0465 in N0-full (2,388,205 pairs of 2,186 families) and 0.0486 in N0-strict, 46.5x and 48.6x the 0.001 of
a calibrated test (criterion >= 5x); lambda from the median p is 0.544 and 0.579, below 1 because the discrete test is
conservative (measured, run R2). The signal world gives 42.2x (measured, dev). The ablation of §4.3 is not done.

## N1: kept and destroyed

| quantity | default null seed | null seed 20261113 | label |
|---|---|---|---|
| module members relabeled | 27 of 100 (12 of 25 modules) | 28 of 100 (13 modules) | measured, run R2 |
| evaluated label pairs (at least one relabeled member) | 47 | 52 | measured |
| Hamming distance of a replacement's profile | <= 4 genomes (floor 0.1 x 46) | same | measured |
| mean genomes adjacent, evaluated label pairs | 0.277 | 0.038 | measured |
| same, matched random pairs (2,000 draws) | 0.133 (one-sided p 0.170) | 0.150 (p 0.715) | measured |
| signal module pairs (positive control), 175 pairs | 16.33 against 0.181 (p 0.0005, the floor) | same | measured |

The 73 unmatched members are module families for which no free root or originated family has both a profile within
4 genomes and a transfer count within 25%; their transfer counts run from 4 to 234, 51 of them 43 or more (measured,
run R2), while the 27 replacements drawn transfer 0 to 52 times (measured, the relabel table; mod22, in 3 genomes, is replaced
by four originated families present in 1 genome each).
§4.1's ±25% therefore relabels about a quarter of the members, mostly of the weak modules. For the PI: keep the
spec's bounds (N1 then calibrates on 47 to 52 pairs per world), or relax the transfer tolerance. The relabel table of
the default seed is `2026-11-13_synthworld_nulls/n1_relabel_default_seed.tsv`; N0-full's `truth_null` is
`n0full_truth_null_default_seed.tsv`.

## Time and memory (run R2 twice, 20:47 and 20:55; seven children at once, then three)

| world | wall, run 1 / run 2 | peak RSS (VmHWM), run 1 / run 2 | label |
|---|---|---|---|
| signal S1-S11 | 231.5 s / 239.3 s | 382 / 383 MB | measured |
| N0-full (two processes each run) | 183.9, 185.8 s / 183.3, 184.5 s | 303, 305 / 304, 304 MB | measured |
| N0-strict (two processes each run) | 169.9, 160.6 s / 169.1, 168.6 s | 285, 288 / 290, 286 MB | measured |
| N0-full, N0-strict at 20261113 | 41.6, 28.7 s / 46.3, 30.1 s | 169, 105 / 173, 104 MB | measured |
| N1 (link S1-S10, relabel, S11), three builds | 30.3, 32.2, 29.9 s / 27.6, 26.5, 27.0 s | 342 to 343 / 164 to 168 MB | measured (run 2 reads backbone.json instead of the whole backbone) |
| R2 fixture, all ten | 265.4 s / 268.3 s | - | measured |

The brief's estimate for a null world was about 230 s (a full world). N0 takes 160 to 186 s, 20 to 30% less (*derived*),
because its genomes lack the 21.5% of genes that planted families hold; N1 shares S1-S10 and costs one S11.

## Decisions taken in this session (for the PI to confirm)

- Decoys and NOGD are removed from both N0 variants (§4.1); NOGD rows are also left out of N1's truth.
- CRISPR-Cas loci: N0-full keeps them with their `cas_operon` rows and excludes within-locus Cas pairs from null
  evaluation (read of §4.1, whose N0 bullet says no `cas_operon` rows while its variants put the Cas group in N0-full
  only); N0-strict drops whole loci at S6 (`crispr.arrays` 0), the backbone untouched by that change.
- Vehicle pairs: the chassis pair is excluded; pairs with cargo are null traps (must_not_call).
- Plasmid pairs are null in both N0 variants.
- N1 is applied to the truth at S11 (the genomes do not move); the module definitions take the replacements, which is
  the spec's label swap up to renaming; N1 evaluates only label pairs with a relabeled member.
- Values: N1 bounds 0.1 of genomes and ±25% transfers (spec defaults); replacement kinds root and originated; N2 has
  20 seeds from stream `N2/seeds`; N0 at a seed uses that seed for every stage.

## What remains of §4

- §4.3: the positive control is a test (above); the ablation (switching off confounder groups 2, 3 and 4 one at a
  time) needs three more worlds per run and is not built; nor is a build-time gate that refuses an N0 world.
- N2: the seeds are defined; the 20 worlds are not built here (about 20 x 3 min), and at other seeds the preset gives
  small, module-dominated worlds (above).
- The rest of the card: item (8), §6.2 to §6.4 and §7.
