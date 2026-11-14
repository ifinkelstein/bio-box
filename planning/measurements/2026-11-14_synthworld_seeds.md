# synthworld `micro` over 21 seeds: why it holds its size and planted share only at the default seed

Card F-02, session 10 (branch `task/F-02-10`). **Outcome: stopped and reported (AGENTS §9; the brief's stop rule).** The
measurement shows that the main cause is not the one the brief described, and no change tried that is small puts the 21
seeds inside the spec's ranges while keeping the strength grid. The generator is unchanged; no digest moves. This note
holds the evidence and the options for the PI.

Code: `scripts/stress/f0210_seed_scan.py` (backbone-only scan over seeds, every number counted from the emitted
tables; `--crown` and `--rule` are measurement devices that patch ZOMBI2's calls in the worker process, never the
generator), `scripts/stress/f0210_world_summary.py` (a full world S1-S11 with stage times, and a summary counted from
its files). Tables beside this note in `2026-11-14_synthworld_seeds/`: `<label>.tsv` (one row per seed),
`<label>_grid.tsv` (one row per seed and planted module), `<label>.json` (the overrides), `default_world_summary.json`
(the current default world, S1-S11). A draft of the test that would have caught this is
`2026-11-14_synthworld_seeds/draft_test_synthworld_preset_seeds.py` (not under `tests/`: it fails on the current preset
by design, and no fix lands in this session).

Host and commands for every *measured* number: wilkcomp02, 2026-11-14 00:00-01:00 CST (night; the account at about
36 to 42 of its 56 threads with this session's 2 to 6), image `bio-box-test_64f1bf0.sif`, units `f0210-*` at nice 19,
one thread per process. **Scan**: `python scripts/stress/f0210_seed_scan.py scan /tmp/if2453/f0210/<dir> <label>
--jobs <2..5> --seeds <...> [options]`, options per label in the table of §2. **World**: `python
scripts/stress/f0210_world_summary.py build /tmp/if2453/f0210/world-before 20261112`, then `summary`.

Seeds: the default 20261112 and the 20 N2 seeds `nulls.n2_seeds(20261112, 20)`. Session 9's second seed 20261113 is
not an N2 seed; it is in the 6-seed cause experiments. "Crown age" is the time from the first speciation of the complete
tree to the present (the scramble stem of 1.0 not counted); session 9's "about 4 at the default seed" is the tree height
with the stem (3.56).

## 1. The 21 seeds on the current preset (step 1)

Spec ranges used (each stated where it comes from): §5.1 chromosome mean 1,500 +- 30% (1,050-1,950;
`test_micro_backbone_parameters`); the size range's ends, smallest chromosome <= 1,100 and largest >= 3,600
(`test_lineage_varying_rates_reach_the_genome_size_range`, at the default seed); §2.1 S2 / §5.1 0-2 plasmids per genome
with both types realized (`test_plasmids_are_whole_replicons_of_their_preset_size`); every planted module with all its
members in at least 9 genomes (PLAN.md §4.3 `min_detectable_count`, 6-9 co-occurrences, the upper end; the spec has no
stricter number); DEC-15 (c) thresholds of `backbone.check_scrambled` and every scored unit intact (DEC-15 (a)). The
spec states no planted-share range; the draft test derives one: 25 modules of sizes 2-6 (100 families), 3 profile-only
groups of 3 and 5 decoys are 114 families, 7.6% of a 1,500-gene chromosome at one copy each; the bound is twice that,
15.2% (*derived*, for the PI to confirm).

All *measured* (scan `base`, `--jobs 5`, seeds default,n2; table `base.tsv`). Shares are of all genes; "planted" is the
planted-module families only (the planted + profile-only + decoy share is `share_planted_all` in the TSV).

| seed | crown age | tree height | genomes | genes per genome mean (min-max) | chromosome mean (min-max) | planted / profile-only / decoy / vehicle / hub / plasmid / background share | plasmid genomes pls0, pls1 (0/1/2 per genome) | modules all members in >= 9 genomes | copies per present member | DEC-15 (c) plus, kept, clustered; scored units intact | ranges broken |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 20261112 | 2.5632 | 3.5632 | 46 | 1937.0 (969-3918) | 1872.7 (969-3798) | 19.7 / 1.5 / 0.3 / 4.6 / 2.5 / 3.3 / 68.0 % | 18, 20 (17/20/9) | 24 | 3.762 | 0.482, 0.0019, 0.021; 3309/3309 | modules<9: 1 |
| 1001685937451957930 | 5.2636 | 6.2636 | 46 | 404.4 (250-1685) | 404.4 (250-1685) | 59.5 / 0.8 / 0.0 / 3.2 / 12.0 / 0.0 / 24.5 % | 0, 0 (46/0/0) | 8 | 3.127 | 0.429, 0.0010, 0.000; 1770/1770 | mean, max<3600, plasmid lost, modules<9: 17 |
| 6329801112132740259 | 4.1528 | 5.1528 | 46 | 991.8 (343-1857) | 987.5 (343-1857) | 35.2 / 1.9 / 0.4 / 7.9 / 4.8 / 0.4 / 49.2 % | 0, 5 (41/5/0) | 23 | 3.779 | 0.469, 0.0019, 0.032; 2832/2832 | mean, max<3600, plasmid lost, modules<9: 2 |
| 6036017085589997374 | 3.9808 | 4.9808 | 46 | 14151.2 (357-38757) | 14089.5 (357-38757) | 3.8 / 0.4 / 0.1 / 0.7 / 0.3 / 0.4 / 94.3 % | 22, 5 (22/21/3) | 20 | 6.121 | 0.573, 0.0057, 0.032; 4161/4161 | mean, modules<9: 5 |
| 4383110225690451409 | 6.156 | 7.156 | 46 | 2152.1 (319-4016) | 2112.9 (319-3896) | 22.5 / 2.3 / 0.2 / 5.0 / 1.8 / 1.8 / 66.5 % | 15, 0 (31/15/0) | 20 | 4.901 | 0.459, 0.0000, 0.032; 2990/2990 | mean, plasmid lost, modules<9: 5 |
| 10975598825187167680 | 4.76 | 5.76 | 46 | 3429.8 (538-8135) | 3421.9 (538-8135) | 15.5 / 1.0 / 0.2 / 1.8 / 1.4 / 0.2 / 79.8 % | 3, 0 (44/1/1) | 24 | 5.741 | 0.436, 0.0029, 0.010; 3799/3799 | mean, plasmid lost, modules<9: 1 |
| 5467273169533009601 | 3.7091 | 4.7091 | 46 | 1319.2 (662-2233) | 1240.1 (542-2088) | 30.9 / 3.5 / 0.2 / 4.2 / 3.7 / 6.0 / 51.5 % | 28, 7 (14/29/3) | 23 | 4.16 | 0.454, 0.0029, 0.021; 3274/3274 | max<3600, modules<9: 2 |
| 4734920099712499969 | 4.8146 | 5.8146 | 46 | 1738.8 (349-2025) | 1738.8 (349-2025) | 22.9 / 3.2 / 0.1 / 8.3 / 2.8 / 0.0 / 62.6 % | 0, 0 (46/0/0) | 21 | 4.257 | 0.554, 0.0019, 0.032; 3153/3153 | max<3600, plasmid lost, modules<9: 4 |
| 1570278518166000558 | 5.3348 | 6.3348 | 46 | 3259.5 (2261-3681) | 3257.8 (2221-3681) | 18.5 / 1.8 / 0.0 / 3.5 / 1.5 / 0.1 / 74.6 % | 0, 2 (44/2/0) | 23 | 5.845 | 0.534, 0.0010, 0.010; 4300/4300 | mean, min>1100, plasmid lost, modules<9: 2 |
| 3830344756309592292 | 5.5021 | 6.5021 | 46 | 3006.2 (477-6580) | 2982.7 (477-6540) | 18.2 / 1.9 / 0.3 / 4.2 / 1.6 / 0.8 / 72.9 % | 0, 27 (19/27/0) | 25 | 5.885 | 0.554, 0.0000, 0.000; 4178/4178 | mean, plasmid lost |
| 9254619518473077140 | 4.6847 | 5.6847 | 46 | 3759.6 (464-7204) | 3759.6 (464-7204) | 11.2 / 0.8 / 0.0 / 3.1 / 0.8 / 0.0 / 84.1 % | 0, 0 (46/0/0) | 16 | 3.844 | 0.417, 0.0000, 0.010; 3401/3401 | mean, plasmid lost, modules<9: 9 |
| 4701969716356831982 | 5.9082 | 6.9082 | 46 | 1000.0 (533-1674) | 971.3 (533-1674) | 35.5 / 1.2 / 0.4 / 3.6 / 4.9 / 2.9 / 51.4 % | 11, 0 (35/11/0) | 15 | 3.965 | 0.450, 0.0029, 0.021; 2453/2453 | mean, max<3600, plasmid lost, modules<9: 10 |
| 14843108908847119356 | 4.3471 | 5.3471 | 46 | 15445.2 (702-39465) | 15445.2 (702-39465) | 4.6 / 0.4 / 0.1 / 1.0 / 0.3 / 0.0 / 93.5 % | 0, 0 (46/0/0) | 24 | 7.025 | 0.452, 0.0000, 0.043; 5395/5395 | mean, plasmid lost, modules<9: 1 |
| 1758546030705738202 | 6.5563 | 7.5563 | 46 | 1727.5 (402-2454) | 1693.6 (379-2454) | 31.3 / 2.9 / 0.0 / 1.7 / 2.8 / 2.0 / 59.3 % | 13, 0 (33/13/0) | 20 | 5.396 | 0.538, 0.0029, 0.032; 3849/3849 | max<3600, plasmid lost, modules<9: 5 |
| 1995372076782357328 | 6.2258 | 7.2258 | 46 | 2896.3 (296-5021) | 2878.1 (296-4981) | 17.3 / 1.5 / 0.0 / 1.3 / 1.7 / 0.6 / 77.5 % | 0, 21 (25/21/0) | 22 | 5.867 | 0.554, 0.0000, 0.021; 3586/3586 | mean, plasmid lost, modules<9: 3 |
| 13514522956726874792 | 7.6234 | 8.6234 | 46 | 446.9 (325-563) | 446.9 (325-563) | 63.6 / 7.3 / 0.0 / 9.9 / 8.7 / 0.0 / 10.5 % | 0, 0 (46/0/0) | 10 | 2.887 | 0.480, 0.0010, 0.021; 1716/1716 | mean, max<3600, plasmid lost, modules<9: 15 |
| 12621855678475765395 | 6.7251 | 7.7251 | 46 | 636.7 (254-1570) | 636.7 (254-1570) | 39.5 / 3.7 / 0.0 / 12.4 / 7.8 / 0.0 / 36.6 % | 0, 0 (46/0/0) | 10 | 2.944 | 0.429, 0.0000, 0.021; 1477/1477 | mean, max<3600, plasmid lost, modules<9: 15 |
| 7689157899245528466 | 3.7641 | 4.7641 | 46 | 1284.8 (338-1936) | 1239.6 (298-1776) | 28.4 / 3.3 / 0.5 / 5.1 / 3.8 / 3.5 / 55.4 % | 7, 31 (15/24/7) | 24 | 4.043 | 0.537, 0.0038, 0.032; 3122/3122 | max<3600, modules<9: 1 |
| 13748185708697958051 | 6.72 | 7.72 | 46 | 11192.2 (462-21685) | 11192.2 (462-21685) | 6.1 / 0.6 / 0.0 / 1.1 / 0.4 / 0.0 / 91.7 % | 0, 0 (46/0/0) | 24 | 6.84 | 0.473, 0.0038, 0.077; 4426/4426 | mean, plasmid lost, modules<9: 1 |
| 15390672246448918260 | 4.9484 | 5.9484 | 46 | 1322.6 (508-2971) | 1290.5 (488-2942) | 29.2 / 3.0 / 0.1 / 10.0 / 2.9 / 2.4 / 52.2 % | 12, 1 (33/13/0) | 15 | 3.901 | 0.450, 0.0019, 0.021; 3190/3190 | max<3600, modules<9: 10 |
| 11714203776604047358 | 5.6616 | 6.6616 | 46 | 675.4 (436-1168) | 669.3 (436-1128) | 55.1 / 3.9 / 0.5 / 1.1 / 7.3 / 0.9 / 31.2 % | 0, 7 (39/7/0) | 21 | 4.012 | 0.447, 0.0029, 0.032; 2833/2833 | mean, max<3600, plasmid lost, modules<9: 4 |

Summary (*derived* from the table): chromosome mean inside 1,050-1,950 at 6 of 21 seeds (404 to 15,445); both ends of
the size range reached at 9; both in band and reaching the ends only at the default seed. Both plasmid types alive at 5
of 21 seeds. All 25 modules powered (>= 9 genomes) at 1 of 21 (3830344756309592292); the default seed misses one (mod22,
3 genomes, as session 9 found). Planted + profile-only + decoy share 4.3% to 70.8% (median 26.3%); planted-module genes
per genome 241 to 709 (median 407), at 2.9 to 7.0 copies per genome of a present member. DEC-15 holds at every seed:
+ fraction 0.417 to 0.573, initial adjacency kept <= 0.0057, declared clustered <= 0.077, every scored unit intact.
The default seed has the shortest crown age of the 21 (2.56; the N2 seeds 3.71 to 7.62, median of all 21 5.26).

Realized strength grid over the 21 seeds (`base_grid.tsv`; medians over the 2-3 modules of a cell and the seeds;
*measured*): member transfer rate 0.05 / 0.15 / 0.30 gives 11.5-20.7 / 35-131 / 97-198 transfer events per member
family, and realized together (share of a module's transfer events that carry every member) is 0.12-0.15, 0.42-0.50
and 1.00 for together 0.3, 0.7 and 1.0. Seeds whose cell has no module with all members in >= 9 genomes: 3, 3 and 1 of
21 in the three cells at member transfer 0.05; 2, 1, 0 at 0.15; 0 at 0.30.

## 2. What causes it (step 2): six seeds, one factor at a time

Seeds: 20261112, 20261113, 1001685937451957930, 10975598825187167680, 4383110225690451409, 6329801112132740259 (the
default, session 9's second seed, the most collapsed and an exploding N2 seed, the deepest and one more). All
*measured* (tables `<label>.tsv`).

| label | change from the preset | chromosome mean per seed (in the order above) | planted+profile-only+decoy share | modules >= 9 genomes | seeds losing a plasmid type | smallest/largest chromosome reaching 1,100/3,600 |
|---|---|---|---|---|---|---|
| base | none (seed 20261113 from session 9: 701 genes per genome, measured there) | 1873, -, 404, 3422, 2113, 988 | 16.7-60.3% | 8-24 | 4 of 5 | 3 of 5 |
| c_crown | tree rescaled to crown age 2.5632 (`--crown 2.5632`) | 1873, 654, 352, 2190, 1715, 1055 | 20.5-45.9% | 12-25 | 2 of 6 | 2 of 6 |
| noaxes | no lineage-rate axes (`--axes none`), native trees | 1717, 1548, 1699, 1752, 1644, 1825 | 25.4-36.1% | 16-25 | 5 of 6 | 0 of 6 |
| c_crown_noaxes | both | 1717, 1527, 1502, 1713, 1635, 1681 | 24.0-26.8% | 18-25 | 2 of 6 | 0 of 6 |
| c_modbg | planted and profile-only members at the run's transfer 0.1 and loss 0.24 (`--rule modbg`) | 1386, 331, 144, 3040, 1591, 612 | 1.8-3.6% | 1-8 | 3 of 6 | - |
| c_modbal | member loss = member transfer + run duplication (`--rule modbal`) | 1477, 243, 168, 3131, 1557, 616 | 5.0-13.9% | 1-18 | 3 of 6 | - |

What each explains:

- **Genome size is set by the lineage-rate axes, not by tree depth.** With the axes off the chromosome mean is 1,548 to
  1,825 at all six seeds although their crowns run from 2.56 to 8.28, and background genes per genome are 891 to 1,119
  (`noaxes.tsv`). Rescaled to the default's crown age with the axes on, the means still run from 352 to 2,190. Over the
  21 seeds the Spearman correlation of crown age with chromosome mean is -0.24 (*derived*). The mechanism (*derived*
  from the preset): a lineage's genes change at about origination x c_b + (duplication 0.04 x 2 + transfer 0.1 x 3) x
  c_b x N - loss 0.24 x 2 x c_l x N per unit time, with c_l the loss axis and c_b the axis shared by duplication,
  transfer and origination, each a mean-1 log-normal drift with sigma 0.8 per speciation. The gains are 0.38/0.48 of the
  loss at c = 1, so a lineage whose c_b/c_l exceeds 1.26 grows without a bound short of `max_family_size` (up to
  39,465 genes measured) and one below it shrinks toward 100 c_b/(0.48 c_l - 0.38 c_b) genes. The first splits are
  shared by many genomes, so whole worlds go one way: means of 11,192 to 15,445 at three N2 seeds, 404 to 669 at four.
  The axes are also what gives the default seed its spread: with them off, no seed reaches either end of the range
  (`noaxes`: smallest chromosome 1,232 to 1,418, largest 1,685 to 2,248).
- **The planted share comes from the module rate imbalance, and depth adds little to it.** Planted-module genes per
  genome do not follow crown age (Spearman -0.02 over 21 seeds, *derived*): members gain copies at duplication 0.04 x c_b
  plus transfer 0.05-0.30 per copy against loss 0.03, so they fill `max_family_size` (copies per present member 2.9 to
  7.0; up to 9.6 in the strongest cell). The share then varies mostly with the background: 4.3% where the background
  exploded, 70.8% where it collapsed. Without axes it is 24-36% at every seed. Balancing the members' rates removes the
  excess (share 1.8-13.9%) but leaves 1 to 18 of 25 modules with all members in 9 genomes: a module starts as one unit
  in the root, and a critical or subcritical copy number loses it in most lineages (losses also come from background
  runs that end in a unit, about 0.24 per unit per unit time, *derived*: loss rate x (mean extent - 1)).
- **Plasmids depend on depth and on early chance.** With transfer 0.3 below loss 0.7 per plasmid and entry exclusion, a
  type's carrier share falls about as exp(-0.4 x crown age) (*derived*): 0.36 at the default's 2.56, 0.05 at 7.6. At a
  fixed crown age a type still dies out at 2 of 6 seeds (c_crown_noaxes), by losses in the few lineages near the crown.
- **Module loss at some seeds is a property of the tree.** At 20261113 one daughter of the first speciation has no
  extant descendant (43 of 43 ZOMBI2 tips descend from the other; at the default seed 3 and 40, measured with
  `zombi2.species.simulate_species_tree` at the S1 seeds), and 7 to 24 modules fall below 9 genomes under every
  variant tried. (The crown age the taxonomy and these tables use is measured from that first speciation, so at such
  seeds the extant crown is younger than the "crown age".)

## 3. What was tried toward a fix, and why it stops here (step 3)

| label | change | chromosome mean | share | modules >= 9 | plasmid type lost | ends reached | strength axis |
|---|---|---|---|---|---|---|---|
| o_crown_noaxes_excl | crown 2.5632, no axes, members with entry exclusion (a lineage holding the group's first family receives none of its families by transfer) and no duplication of their own (`--rule excl+dup0`) | 1089-1324 | 5.6-8.0% | 16-25 | 2 of 6 | 0 of 6 | transfers per member 4.0-11.0 in every cell (base: 11.5 to 198) |
| o_crown_orig16_excl | as above, and one axis on origination alone, sigma 1.6 (`--axes origination:1.6`) | 919-1375 | 4.9-9.9% | 17-25 | 4 of 6 | 0 of 6 | as above |

(Earlier probes with the generator patched in a scratch tree, not kept: crown 2.5632 with origination and the other
three rates on separate axes at sigma 0.8, means 1,914 to 2,556; origination alone at sigma 1.2, smallest chromosome
above 1,100 and largest below 3,600 at the default seed. Their numbers are not used here.)

Every element of the brief's candidate list was measured, and none suffices:

1. *A fixed crown age* (rescaled tree; it leaves the taxonomy cut fractions, the burst depth and every per-split lineage
   draw unchanged, and only stretches time) fixes the plasmid decay and part of the module loss, but not genome size
   (352 to 2,190).
2. *Per-family rates per unit of expected substitutions* is the same change as 1 for the backbone (one time scale).
3. *Module loss balancing the gain* gives a stationary mean copy number but loses most modules (1 to 18 of 25 powered).
4. *Capping the planted share by design* (entry exclusion, no own duplication) keeps every grid cell powered at all six
   seeds (`o_crown_noaxes_excl_grid.tsv`: no seed without a powered module in any cell) and the share at 5.6-8.0%, but a
   transfer into a genome already holding the module is refused, so the member transfer axis collapses: 4.0 to 11.0
   transfer events per member family in every cell, against 11.5 / 35-131 / 97-198 for 0.05 / 0.15 / 0.30 now. The
   PI's grid (2026-11-12) would no longer mean what it says.

And genome size needs the axes redesigned: the only mechanism measured that reaches both ends of the size range at the
default seed (pinned by `test_lineage_varying_rates_reach_the_genome_size_range`, which this session may not loosen) is
the divergent one that also sends other seeds to 404 or 15,445 genes. A linear one (origination on its own axis) keeps
the mean but did not reach the ends at sigma up to 1.6 within the 2.56 time units of the default crown.

A fix therefore touches four coupled things at once: the time scale, the lineage axes (and with them the size-range
requirement), the module dynamics (and with them the grid's transfer axis) and the plasmid rates, each needing its own
21-seed scan and refit. That is the large redesign the brief rules out for this session.

## 4. Options for the PI

1. **Decide what "inside the range at any seed" means for genome size.** (a) Per seed: chromosome mean 1,500 +- 30% only,
   with the 1,000-4,000 spread asked of the default seed or dropped. Then the axes can go or become linear (origination
   alone), and `noaxes` / `c_crown_noaxes` show the mean holding at 1,502-1,825 at six seeds. (b) Per seed, both the mean
   and the spread: needs a new size mechanism (for example origination on its own axis with a longer fixed crown age so
   sizes relax within the tree, 1 - exp(-0.1 x T) of the way), refitted over the 21 seeds.
2. **Decide what the member transfer axis must mean.** (a) Keep raw rates 0.05/0.15/0.30 and accept saturation, with a
   genome cap for modules instead of exclusion (ZOMBI2 0.54.0 has only the global `max_family_size`; a per-family cap
   would have to be written in synthworld). (b) Keep the cap by entry exclusion and redefine the axis as realized
   transfer events per member (then refit rates so the three cells differ). (c) Balance members' loss partly (between
   0.03 and transfer + duplication) and accept fewer powered modules, scored with `powered_expected` (§3.3).
3. **Fix the time scale** (`species.crown_age`, a tree rescaled to it): needed under every option for plasmids and module
   loss; the value is a choice (the default seed's 2.5632 keeps the rates' fits; the 21-seed median 5.26 lets sizes
   relax further but needs plasmid loss refitted, since their prevalence goes as exp(-(loss - transfer) x T)).
4. **Plasmids**: with entry exclusion a type survives only if early losses miss; lower turnover (for example transfer
   0.15, loss 0.55 at the same net decline) or a fixed crown age with a 21-seed scan of plasmid loss.
5. **Seeds whose first speciation has one extinct side** (20261113, 6329801112132740259 of the six): either measure the
   crown from the extant MRCA (it moves the taxonomy cuts) or accept fewer powered modules there.

The signal world and every digest stay as sessions 1-9 left them. The draft test
(`draft_test_synthworld_preset_seeds.py`) states the per-seed checks a fix must pass and fails on the current preset:
3 of 3 failed in 648 s (measured, wilkcomp02, 2026-11-14, `sh scripts/pod_pytest.sh F-02-10
tests/r2/test_synthworld_preset_seeds.py` with the draft in `tests/r2/`). Its own counts agree with the scan's (chromosome
means 1,872.7, 404.4, 987.5, 14,089.5 and 2,112.9) and give session 9's second seed 20261113 a chromosome mean of 700.9,
a planted + profile-only + decoy share of 57.3%, no plasmid, and a grid cell (together 0.3, transfer 0.15) whose three
modules are in no genome whole; at 1001685937451957930 three cells have no module whole in 9 genomes.

## 5. The current default world, for the before-and-after table of the next session

`default_world_summary.json` (S1-S11 at 20261112, *measured*, the build command above): stage digests equal to the
pinned ones (backbone `025fdd0d…`, sequences `34ebce13…`, replicons `55bda971…`, mobile `0bd7f0db…`, fusions
`db09abf3…`, assembly `2f425158…`, sources `770007bd…`, emit `a2f0e448…`, world `f0262027e4ccdcbc72cd7d90f23fbb5c`);
stage times 37.1, 86.7, 48.0, 11.5, 5.3, 18.9, 6.6 and 24.5 s, 238.6 s in all; genes per genome 1,937.0 (969-3,918);
planted share 21.5%; plasmid genomes pls0 18, pls1 20; with every TSV's row count, family counts by kind, phages,
prophage instances, arrays, fusion events, call errors, assembly classes and sources by kind.
