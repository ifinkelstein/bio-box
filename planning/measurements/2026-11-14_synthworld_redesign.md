# synthworld `micro` redesigned to hold at every seed (DEC-36, DEC-37)

Card F-02, session 11 (branch `task/F-02-11`). The redesign of `planning/measurements/2026-11-14_synthworld_seeds.md`
§4 as the PI decided it on 2026-11-14 (DEC-36: entry exclusion and a realized transfer axis; DEC-37: genome size as
a mean, one origination axis, a fixed crown age from the living MRCA, refitted plasmid rates). Tables beside this note
in `2026-11-14_synthworld_redesign/`: `<label>.tsv` (one row per seed), `<label>_grid.tsv` (one row per seed and
planted module), `<label>.json` (the overrides), `crowns.tsv` (the species trees alone), `after*.tsv` (the preset as
committed), `default_world_summary_after.json` (the new default world, S1-S11).

Host and commands for every *measured* number: wilkcomp02, 2026-11-14, 11:00-19:00 CST (Saturday; the account at
48 to 51 of its 56 threads with the release build `r50-rM2-14` and another agent, this session's 4), image
`bio-box-test_faec1f4.sif`, units `f0211-*` at nice 19, one thread per process. **Scan**: `python
scripts/stress/f0210_seed_scan.py scan /tmp/if2453/f0211/<dir> <label> --jobs <2|4> --seeds default,20261113,n2
[--n0] [--axes origination:S] [--set a.b=JSON ...]` (overrides per label in `<label>.json`). **Crowns**: `python
scripts/stress/f0210_seed_scan.py crowns /tmp/if2453/f0211/scan/crowns.tsv`. Seeds: the default 20261112, session 9's
second seed 20261113 and the 20 N2 seeds `nulls.n2_seeds(20261112, 20)` (22 rows; "the 21 seeds" are the default and
the N2 seeds).

**Outcome: the redesign is in the generator and the preset; three existing assertions cannot hold under it and are
reported, not loosened (§4).** DEC-36 and DEC-37 (1) and (2) hold at all 22 seeds; DEC-37 (3) holds as a stated share
(both plasmid types at 15 of 22 seeds). The stage digests of the default world changed once and are re-pinned (§5).

## 1. What the generator does now

- **Fixed crown age from the living MRCA** (`species.crown_age` 5.0; `backbone._run_zombi2`, `extant_crown`): the drawn
  tree is rescaled (`zombi2.tree.rescale`, every branch times one factor) so that the most recent common ancestor of the
  living lineages lies 5.0 before the present. At 8 of the 22 seeds one side of the first speciation has no living
  genome (`crowns.tsv`: `mrca_is_root` False), and the first speciation then lies 5.0 to 9.49 before the present.
- **Rates start at the crown** (`species.rate_onset: crown`): duplication, transfer, loss and origination (and every
  family's own rates) are 0 from the first speciation to the MRCA, so every world evolves its gene content for exactly
  5.0 time units; inversions and transpositions run from the first speciation as before (they move no gene content).
- **What the rescaling does** (*derived* from the code): the taxonomy cut fractions (0.08 ... 0.9) are unchanged and are
  now fractions of the time from the MRCA to the present, so the cut times scale with the factor and move to the MRCA
  at the 8 seeds above; the clonal burst is rooted 0.02 x 5.0 = 0.1 before the present at every seed; the lineage
  multipliers are drift steps drawn at each speciation (`Drift(LogNormal(0, sigma))`), so their number and distribution
  do not change, only how long each acts. At the default seed the factor is 5.0 / 2.5632 = 1.95: its world is twice as
  old as before.
- **One lineage axis, on origination alone** (`lineage_rates.axes` [{origination, sigma 0.2}]).
- **Entry exclusion** (DEC-36; `modules.entry_exclusion`, `profile_only.entry_exclusion`): ZOMBI2 0.54.0 has no
  per-group knob, so the adapter gives every member family the recipient rule `Recipients().weighted_by("genomes:<group>.0",
  {present: 0, absent: 1})`; ZOMBI2 normalises the weights over the contemporaneous lineages and makes a transfer with no
  eligible recipient a no-op (`zombi2/genomes/_transfer.py`). This rule replaces the run's distance decay for the
  members' own transfers (one weighting per family in 0.54.0); a segment that also carries background genes multiplies
  the rules. Members start no duplication of their own (`member_duplication` 0).
- **The truth records the realized strength** (`emit.py`, `truth_module`): `realized_transfers_per_member` (transfer
  rows of the backbone's `events.tsv` on the group's families over its size) and `realized_together` (share of the
  distinct transfer events touching a member that carry every member), on every member row of a planted, profile-only
  or vehicle group.

Tested from the emitted tables, never from the generator's summary: `test_dec36_no_module_family_is_transferred_into_a_lineage_that_holds_the_module`
replays ZOMBI2's `genome_events.tsv` (all lineages, speciation rows included) and counts, at the default seed, 2,117
module-family transfer rows, 0 into a lineage holding the module's first family (728 into a lineage holding some other
member, which partial transfers at together < 1 allow); the same world without the rule has 10,643 rows and 10,216
into a holder (*measured*, `pod_pytest.sh F-02-11 -s tests/r2/test_synthworld_backbone.py`). The per-seed test replays
the rule at six seeds and the scan at all 22 (`module_transfers_into_first_holder` 0 of 644 to 2,123 at every seed).
`test_species_tree_has_the_preset_crown_age_from_the_living_mrca` reads the crown from the Newick with a parser of its
own; `test_truth_module_records_the_realized_strength` recounts the two truth columns from `events.tsv`.

## 2. The crown age: 5.0

DEC-37 (2) asks for about the 21-seed median. Over the default and the 20 N2 seeds the crown measured from the first
speciation (session 10's definition) has median 5.26 (*measured*, session 10, `base.tsv`), and the crown measured from
the living MRCA, as DEC-37 (3) now defines it, has median 4.52 (min 2.56, max 6.72; *measured*, `crowns.tsv`, the
`crowns` command above). 5.0 lies between the two medians and is a round value; every fit below is at 5.0.

## 3. The fits (22 seeds each; one row per scan, every candidate tried)

Columns (*measured*, from each `<label>.tsv` and `_grid.tsv` by `fitsumm.py`-style counting): chromosome mean range over
the seeds and seeds inside 1,050-1,950; seeds with both plasmid types; seeds with a genome of more than 2 plasmids;
pooled medians of realized transfer events per member family by member-transfer level; seeds whose three level medians
strictly increase; N0-full against the signal, Pearson r of background genome sizes at the default seed / 20261113 /
seeds with r >= 0.9 (`--n0`); seeds with all nine grid cells powered (a module with every member in >= 9 genomes);
planted + profile-only + decoy share. Batch C (`c_*`) is the earlier worker's (origination 0.15, sigma 0.3 unless named,
member transfer 0.02/0.15/0.6 unless named, plasmids 0.3/0.7 unless named); E (`e_*`: member transfer 0.01/0.06/0.4,
sigma 0.7 unless named, plasmids 0.3/0.7); F (`f_*`: sigma 0.25, member transfer 0.01/0.06/0.4, plasmids 1.0/0.5
unless named); G (`g_*`, `after`: plasmids 0.8/0.4, member transfer 0.01/0.06/0.4); H (`h_*`: the committed preset and
the spike's background rates).

| scan | chromosome mean | in band | both plasmid types | a genome > 2 plasmids | realized transfers per member, pooled medians | levels increasing | N0 r default / 20261113 / >= 0.9 | all cells powered | planted share |
|---|---|---|---|---|---|---|---|---|---|
| c_fit1_o015_s0.2 | 1145-1733 | 22 | 3 | 2 | 0.02: 5.3, 0.15: 19.0, 0.6: 19.2 | 13 | 0.71 / 0.06 / 7 | 16 | 3.7-8.8% |
| c_fit1_o015_s0.35 | 986-2029 | 20 | 4 | 0 | 0.02: 5.7, 0.15: 17.0, 0.6: 19.4 | 16 | 0.53 / 0.12 / 12 | 18 | 3.5-11.2% |
| e_s0.5 | 954-2557 | 16 | 5 | 0 | 0.01: 5.0, 0.06: 10.0, 0.4: 20.8 | 21 | 0.72 / 0.69 / 12 | 12 | 3.0-8.6% |
| e_s0.7 | 775-3384 | 13 | 4 | 0 | 0.01: 5.0, 0.06: 9.9, 0.4: 19.6 | 20 | 0.85 / 0.50 / 11 | 20 | 2.1-9.5% |
| e_s1.0 | 684-4858 | 11 | 2 | 0 | 0.01: 4.7, 0.06: 9.8, 0.4: 20.6 | 20 | 0.66 / 0.21 / 9 | 17 | 1.4-11.9% |
| f_s0.2 | 1136-1808 | 22 | 10 | 6 | 0.01: 5.0, 0.06: 11.0, 0.4: 18.4 | 20 | 0.68 / 0.66 / 7 | 13 | 4.0-6.8% |
| f_s0.25 | 1116-1834 | 22 | 14 | 7 | 0.01: 5.0, 0.06: 11.0, 0.4: 19.0 | 20 | 0.73 / 0.77 / 9 | 15 | 3.4-7.3% |
| f_s0.3 | 986-2000 | 20 | 13 | 10 | 0.01: 5.5, 0.06: 9.7, 0.4: 20.2 | 20 | 0.85 / 0.68 / 9 | 17 | 3.4-7.5% |
| g_s0.25_o0.15 (sigma 0.25, origination 0.15) | 1024-1844 | 21 | 11 | 9 | 0.01: 5.0, 0.06: 9.8, 0.4: 19.7 | 20 | 0.60 / -0.21 / 8 | 14 | 3.2-7.3% |
| g_s0.25_o0.16 | 1134-1963 | 21 | 12 | 6 | 0.01: 5.0, 0.06: 10.8, 0.4: 19.0 | 20 | 0.78 / 0.83 / 9 | 15 | 3.4-7.6% |
| **g_s0.2_o0.155 (committed)** | 1149-1795 | 22 | 15 | 11 | 0.01: 5.0, 0.06: 10.9, 0.4: 19.0 | 22 | 0.79 / 0.74 / 8 | 18 | 3.2-7.6% |
| c_fit3_mt_0.01,0.1,0.5 | 1054-1984 | 21 | 3 | 0 | 0.01: 5.5, 0.1: 16.0, 0.5: 20.2 | 18 | - | 18 | 3.4-8.2% |
| c_fit3_mt_0.02,0.15,0.8 | 1113-1878 | 22 | 3 | 1 | 0.02: 5.3, 0.15: 17.7, 0.8: 20.2 | 14 | - | 18 | 4.7-8.8% |
| c_fit3_mt_0.02,0.2,1.5 | 1073-2051 | 21 | 2 | 1 | 0.02: 6.0, 0.2: 20.0, 1.5: 20.6 | 12 | - | 21 | 4.2-9.7% |
| e_mt_0.005,0.05,0.4 | 837-3564 | 12 | 3 | 0 | 0.005: 4.4, 0.05: 9.5, 0.4: 18.5 | 20 | - | 20 | 2.8-10.2% |
| e_mt_0.01,0.04,0.25 | 780-3471 | 15 | 6 | 0 | 0.01: 5.0, 0.04: 8.0, 0.25: 18.0 | 17 | - | 14 | 2.1-9.7% |
| e_mt_0.003,0.03,0.3 | 842-3437 | 14 | 5 | 0 | 0.003: 4.1, 0.03: 8.0, 0.3: 21.0 | 21 | - | 18 | 2.2-8.7% |
| c_fit2_pls_0.5-0.2 | 952-1965 | 20 | 18 | 16 | 0.02: 5.7, 0.15: 17.6, 0.6: 21.0 | 15 | - | 21 | 4.2-8.9% |
| c_fit2_pls_0.6-0.3 | 1028-1934 | 21 | 13 | 10 | 0.02: 5.1, 0.15: 19.0, 0.6: 20.8 | 16 | - | 17 | 3.5-9.7% |
| c_fit2_pls_0.8-0.4 | 1101-1957 | 21 | 16 | 6 | 0.02: 5.8, 0.15: 18.2, 0.6: 20.9 | 16 | - | 19 | 4.3-8.2% |
| e_pls_0.2-0.25 | 863-3333 | 13 | 14 | 5 | 0.01: 4.5, 0.06: 10.0, 0.4: 20.0 | 20 | - | 17 | 2.1-9.0% |
| f_pls_0.8-0.3 | 1103-1919 | 22 | 10 | 5 | 0.01: 5.9, 0.06: 9.3, 0.4: 20.8 | 19 | - | 16 | 3.1-7.4% |
| f_pls_1.0-0.5 | 1116-1834 | 22 | 14 | 7 | 0.01: 5.0, 0.06: 11.0, 0.4: 19.0 | 20 | - | 15 | 3.4-7.3% |
| f_pls_1.2-0.6 | 1114-1881 | 22 | 10 | 6 | 0.01: 5.2, 0.06: 10.0, 0.4: 20.4 | 21 | - | 19 | 3.9-7.5% |
| f_pls_1.0-0.3 | 1106-1865 | 22 | 14 | 11 | 0.01: 4.6, 0.06: 10.4, 0.4: 20.0 | 22 | - | 15 | 3.5-7.2% |
| f_pls_1.5-0.8 | 1067-1880 | 22 | 12 | 9 | 0.01: 4.6, 0.06: 10.0, 0.4: 20.6 | 21 | - | 18 | 3.4-7.0% |
| h_spike (committed, background 0.02/0.05/0.12) | 1471-2196 | 16 | 17 | 10 | 0.01: 3.0, 0.06: 9.0, 0.4: 14.0 | 20 | - | 17 | 2.8-5.7% |

`f_s0.25` and `f_pls_1.0-0.5` are the same settings run twice in two units and gave identical tables (`h_ref` and
`g_s0.2_o0.155` likewise): the scan is deterministic across processes.

What each fit decided:

- **Origination axis** (DEC-37 (1)): the mean leaves the band at 1 to 2 seeds from sigma 0.25 (depending on the draw)
  and at 6, 9 and 11 seeds at 0.5, 0.7 and 1.0; at sigma 0.2 with origination 0.155 x root families every seed holds
  (1,149-1,795; mean of the 22 means 1,454, *derived*). Size still moves with the axis: at the default seed a
  genome's originated genes rank with the origination multiplier its line carried since the MRCA at Spearman 0.80
  (permutation p 0.0005; *measured*, the migrated test).
- **Member transfer** (DEC-36): under exclusion a member is transferred only into a lineage without the module, so the
  realized count has a floor near 4 (background transfer runs that carry a member) and saturates near 20 (the holes that
  losses make, about 0.27 per unit per unit time with background loss runs, *derived*: 0.03 + 0.24); the old rates
  0.05/0.15/0.30 sit at about 5/18/20 (c_fit3_mt_0.02,0.15,0.8 with 0.8 for 0.30). 0.01/0.06/0.4 gives pooled medians
  5.0/10.9/19.0, a factor 2 between levels, and strictly increasing level medians at 22 of 22 seeds.
- **Plasmids** (DEC-37 (3)): see §4.

## 4. What cannot be met, and why (stop and report)

1. **Both plasmid types at every seed** (DEC-37 (3)): reached at 15 of 22 seeds with the committed 0.8/0.4, and at 10 to
   18 of 22 at every pair tried from 0.2/0.25 to 1.5/0.8. A plasmid unit is cured at its own rate plus about 0.24 per
   unit per unit time by background loss runs that end in it (ZOMBI2: "a run that ends inside the group takes the whole
   group"; *derived*: 0.24 x sum over k >= 2 of P(extent >= k) = 0.24 at geometric mean 2), and a lineage that carries
   it can also die; with two lineages at the crown a type is lost early with probability about ((d + L + 0.24) / (b +
   T))^2 per type (*derived*, a branching process of carrier lineages: 0.3 death, 1.0 birth), about 0.2 at 0.8/0.4.
   Making T large enough to push that below 1/42 saturates every genome, and background duplication runs that end in
   the unit then give some genomes a second copy of a type. So survival at every seed and §5.1's "0-2 plasmids per
   genome" with some genomes carrying none (asserted at the default seed by
   `test_plasmids_are_whole_replicons_of_their_preset_size`) exclude each other in this model: at the committed
   rates the default seed has 44 and 46 carriers, no genome without a plasmid (0/2/44), so that assertion fails.
2. **N0's background sizes correlate with the signal's at r >= 0.9** (`test_n0_keeps_the_background_process_in_distribution`):
   the assertion rests on "size is set by the per-lineage rates, which ZOMBI2 draws with the tree" (its docstring), true
   of the two runaway axes DEC-37 (1) removed. With one origination axis the per-lineage spread is small against the
   stochastic spread of a genome's size, so a new draw (N0) of the same lineage rates correlates weakly: r at the
   default seed 0.53 to 0.85 over sigma 0.2 to 1.0, never 0.9, and at sigma 1.0 the default seed's N0 size ratio is
   1.20 (outside [0.9, 1.1]) and 11 of 22 seed means leave the band. No sigma meets both DEC-37 (1) and r >= 0.9 at the
   default seed. The committed preset gives r 0.79 at the default seed and 0.74 at 20261113.
3. **The core genome at the default seed** (`test_hub_families_have_more_diverse_neighbours` needs at least 100 root
   families present in >= 40 of the 46 genomes to compute its comparator): 9 now, 243 in session 2. The default seed's
   tree is stretched 1.95x to the fixed crown, and session 2 had doubled every turnover rate for a crown of 2.56 (the
   spread DEC-37 (1) no longer asks for); a root gene's expected survival on a lineage falls as exp(-0.48 x T). With the
   spike's rates (duplication 0.02, transfer 0.05, loss 0.12) the core is 36 at the default seed (36 to 237 over 22
   seeds, 11 seeds >= 100) and the means rise to 1,471-2,196 (`h_spike`): halving the turnover is the direction, but it
   needs its own 21-seed refit of origination, plasmids and member rates, and a decision on the floor of 100.

Every other range holds at every seed: chromosome mean (22 of 22), crown age 5.0 from the MRCA (22 of 22), entry
exclusion (0 breaches at 22 of 22), planted + profile-only + decoy share 3.2-7.6% (<= 15.2%), copies per present member
0.81 to 1.15 (session 10: 2.9 to 7.0), DEC-15 (c) and every scored unit intact (22 of 22).

## 5. The committed preset over the 22 seeds (scan `g_s0.2_o0.155`, identical to the preset as committed)

| seed | crown age (MRCA) / from first speciation | genomes | genes per genome mean (min-max) | chromosome mean (min-max) | planted / profile-only / decoy / vehicle / hub / plasmid / background share % | plasmid genomes pls0, pls1 (0/1/2 per genome) | modules all members in >= 9 genomes | copies per present member | DEC-15 (c) plus, kept, clustered; scored units intact | realized transfers per member, level medians | N0 background size ratio, r | module transfers into a holder |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 20261112 | 5 / 5 | 46 | 1415.8 (1071-1761) | 1261.0 (911-1601) | 5.4 / 0.7 / 0.6 / 8.6 / 3.5 / 10.9 / 70.3 | 44, 46 (0/2/44) | 23 | 1.145 | 0.460, 0.0000, 0.021; 704/704 | 9.0, 17.0, 32.67 | 1.0263, 0.7905 | 0 of 2117 |
| 20261113 | 5 / 8.435 | 46 | 1188.0 (945-1388) | 1148.9 (945-1348) | 6.5 / 0.9 / 0.3 / 11.9 / 4.1 / 3.3 / 73.1 | 0, 45 (3/41/2) | 22 | 1.044 | 0.459, 0.0019, 0.021; 580/580 | 3.8, 15.0, 22.33 | 1.104, 0.7416 | 0 of 1503 |
| 1001685937451957930 | 5 / 5 | 46 | 1572.1 (1186-1848) | 1524.2 (1026-1805) | 2.6 / 0.4 / 0.2 / 5.4 / 3.1 / 3.0 / 85.2 | 3, 46 (0/43/3) | 14 | 0.812 | 0.524, 0.0048, 0.021; 440/440 | 1.0, 2.0, 16.5 | 1.0395, 0.8489 | 0 of 722 |
| 6329801112132740259 | 5 / 5.136 | 46 | 1887.0 (1134-2212) | 1740.9 (974-2149) | 3.2 / 0.5 / 0.3 / 5.2 / 2.6 / 7.7 / 80.5 | 40, 48 (0/6/38) | 19 | 1.034 | 0.477, 0.0010, 0.032; 638/638 | 8.0, 15.0, 21.0 | 1.0444, 0.9564 | 0 of 1657 |
| 6036017085589997374 | 5 / 5 | 46 | 1607.5 (1305-1902) | 1444.0 (1145-1742) | 3.3 / 0.5 / 0.1 / 8.5 / 2.4 / 10.2 / 74.9 | 47, 47 (0/2/40) | 17 | 0.859 | 0.412, 0.0019, 0.021; 485/485 | 6.0, 11.0, 19.0 | 1.0476, 0.8838 | 0 of 1295 |
| 4383110225690451409 | 5 / 5 | 46 | 1652.2 (1295-2019) | 1485.3 (1135-1859) | 4.0 / 0.6 / 0.2 / 8.7 / 3.0 / 10.1 / 73.5 | 50, 42 (0/7/32) | 18 | 0.994 | 0.432, 0.0010, 0.065; 633/633 | 5.0, 10.67, 21.0 | 1.0421, 0.9018 | 0 of 1234 |
| 10975598825187167680 | 5 / 5 | 46 | 1338.2 (1092-1872) | 1301.6 (1085-1832) | 4.2 / 0.4 / 0.2 / 8.3 / 3.7 / 2.7 / 80.5 | 0, 42 (6/38/2) | 17 | 1.055 | 0.522, 0.0029, 0.021; 460/460 | 3.0, 4.0, 19.25 | 1.0381, 0.8405 | 0 of 945 |
| 5467273169533009601 | 5 / 5.03 | 46 | 1376.0 (890-1689) | 1230.8 (850-1529) | 5.1 / 0.5 / 0.3 / 8.2 / 3.6 / 10.5 / 71.9 | 44, 35 (1/11/34) | 25 | 1.091 | 0.426, 0.0019, 0.043; 661/661 | 7.0, 8.5, 13.33 | 0.9043, 0.5604 | 0 of 1089 |
| 4734920099712499969 | 5 / 5 | 46 | 1535.2 (1102-1783) | 1372.6 (852-1623) | 4.4 / 0.4 / 0.3 / 7.6 / 3.2 / 10.6 / 73.7 | 48, 43 (0/3/41) | 20 | 0.953 | 0.439, 0.0038, 0.032; 649/649 | 5.0, 15.0, 17.0 | 1.0421, 0.9408 | 0 of 1218 |
| 1570278518166000558 | 5 / 6.18 | 46 | 1947.1 (1213-2350) | 1795.0 (1053-2190) | 3.1 / 0.5 / 0.5 / 7.3 / 2.5 / 7.8 / 78.2 | 44, 43 (0/6/39) | 18 | 0.963 | 0.454, 0.0010, 0.054; 547/547 | 2.4, 11.5, 14.5 | 0.9943, 0.9074 | 0 of 1234 |
| 3830344756309592292 | 5 / 5 | 46 | 1508.0 (971-1826) | 1441.9 (931-1701) | 4.9 / 0.4 / 0.2 / 2.3 / 3.2 / 4.4 / 84.7 | 10, 46 (0/36/10) | 20 | 1.028 | 0.435, 0.0010, 0.043; 640/640 | 6.33, 17.67, 18.33 | 0.9779, 0.9589 | 0 of 1665 |
| 9254619518473077140 | 5 / 7.657 | 46 | 1525.0 (1219-1773) | 1365.0 (1060-1613) | 5.0 / 0.5 / 0.5 / 7.4 / 3.2 / 10.5 / 72.9 | 45, 49 (0/1/42) | 22 | 1.048 | 0.525, 0.0010, 0.043; 697/697 | 12.0, 17.8, 27.0 | 1.0255, 0.8569 | 0 of 1822 |
| 4701969716356831982 | 5 / 5.087 | 46 | 1400.7 (1125-1641) | 1400.7 (1125-1641) | 4.8 / 0.5 / 0.2 / 4.1 / 3.5 / 0.0 / 86.9 | 0, 0 (46/0/0) | 19 | 0.969 | 0.428, 0.0000, 0.021; 636/636 | 5.0, 6.0, 18.0 | 1.0156, 0.8266 | 0 of 962 |
| 14843108908847119356 | 5 / 5 | 46 | 1759.5 (1332-2236) | 1759.5 (1332-2236) | 4.0 / 0.5 / 0.2 / 5.2 / 2.8 / 0.0 / 87.2 | 0, 0 (46/0/0) | 19 | 1.015 | 0.550, 0.0029, 0.054; 644/644 | 10.0, 18.0, 25.75 | 1.015, 0.8805 | 0 of 1787 |
| 1758546030705738202 | 5 / 5 | 46 | 1932.7 (930-2261) | 1782.2 (770-2101) | 2.9 / 0.4 / 0.0 / 6.2 / 2.5 / 7.8 / 80.2 | 44, 41 (0/8/37) | 16 | 0.884 | 0.563, 0.0000, 0.032; 509/509 | 1.67, 6.0, 8.0 | 1.0098, 0.9571 | 0 of 644 |
| 1995372076782357328 | 5 / 5 | 46 | 1755.2 (1216-2611) | 1589.9 (1056-2451) | 3.8 / 0.7 / 0.6 / 4.7 / 2.8 / 9.4 / 78.1 | 48, 46 (0/0/45) | 20 | 1.118 | 0.483, 0.0048, 0.043; 633/633 | 3.17, 7.6, 25.0 | 1.0052, 0.9679 | 0 of 1260 |
| 13514522956726874792 | 5 / 9.488 | 46 | 1434.6 (1183-1750) | 1299.8 (1067-1481) | 3.9 / 0.9 / 0.6 / 4.1 / 3.4 / 9.4 / 77.8 | 36, 47 (0/12/31) | 18 | 1.034 | 0.444, 0.0019, 0.054; 540/540 | 5.0, 7.5, 20.67 | 1.018, 0.6872 | 0 of 1082 |
| 12621855678475765395 | 5 / 7.443 | 46 | 1563.4 (1296-1782) | 1519.9 (1256-1742) | 3.2 / 0.8 / 0.1 / 6.7 / 3.1 / 2.8 / 83.3 | 0, 50 (1/43/0) | 15 | 1.038 | 0.468, 0.0019, 0.032; 519/519 | 4.0, 7.75, 16.67 | 1.0375, 0.575 | 0 of 979 |
| 7689157899245528466 | 5 / 5 | 46 | 1500.5 (1088-1942) | 1352.7 (928-1782) | 4.9 / 0.7 / 0.2 / 8.2 / 3.2 / 9.8 / 72.9 | 44, 38 (1/8/37) | 23 | 1.065 | 0.546, 0.0019, 0.010; 719/719 | 12.0, 16.0, 24.0 | 1.0521, 0.8928 | 0 of 2123 |
| 13748185708697958051 | 5 / 5 | 46 | 1314.8 (961-1701) | 1314.8 (961-1701) | 3.5 / 0.6 / 0.0 / 3.6 / 3.7 / 0.0 / 88.6 | 0, 0 (46/0/0) | 16 | 0.971 | 0.550, 0.0010, 0.021; 422/422 | 2.5, 14.0, 19.0 | 1.0317, 0.9501 | 0 of 1127 |
| 15390672246448918260 | 5 / 6.455 | 46 | 1621.9 (1336-1916) | 1580.2 (1296-1876) | 3.1 / 0.4 / 0.3 / 8.6 / 3.0 / 2.6 / 82.0 | 0, 48 (0/44/2) | 18 | 0.926 | 0.416, 0.0000, 0.021; 533/533 | 4.33, 10.67, 25.0 | 0.9532, 0.7738 | 0 of 1091 |
| 11714203776604047358 | 5 / 5.052 | 46 | 1427.2 (1113-1763) | 1272.4 (953-1603) | 4.5 / 0.7 / 0.2 / 11.1 / 3.4 / 10.8 / 69.2 | 44, 46 (1/2/41) | 19 | 1.068 | 0.473, 0.0029, 0.021; 571/571 | 10.0, 12.0, 17.67 | 1.0963, 0.7318 | 0 of 1477 |

Realized strength grid (medians over the modules of a cell and the 22 seeds, range in parentheses; realized together =
share of a module's transfer events carrying every member; powered = seeds at which the cell has a module with every
member in >= 9 genomes; *measured*):

| member transfer \ together | 0.3 | 0.7 | 1 | level, pooled |
|---|---|---|---|---|
| 0.01 | 4.1 (0-19); together 0.17; powered at 21/22 | 5.3 (0-23); together 0.59; powered at 22/22 | 7.0 (0-17); together 1.00; powered at 22/22 | 5.0 |
| 0.06 | 8.8 (0-21); together 0.18; powered at 22/22 | 12.1 (0-31.67); together 0.55; powered at 20/22 | 12.5 (0-33); together 1.00; powered at 22/22 | 10.9 |
| 0.4 | 16.4 (0-43.83); together 0.25; powered at 22/22 | 18.8 (0-40); together 0.62; powered at 21/22 | 23.5 (3-36); together 1.00; powered at 22/22 | 19.0 |

## 6. The default world before and after (S1-S11 at 20261112; `f0210_world_summary.py`)

Before: `2026-11-14_synthworld_seeds/default_world_summary.json` (session 10); after:
`2026-11-14_synthworld_redesign/default_world_summary_after.json` (*measured*, wilkcomp02, 2026-11-14, `python
scripts/stress/f0210_world_summary.py build /tmp/if2453/f0211/world-after 20261112`, then `summary`).

| quantity | before | after |
|---|---|---|
| arrays | 10 | 10 |
| chromosome_genes.max | 3798 | 1601 |
| chromosome_genes.mean | 1872.7 | 1261.0 |
| chromosome_genes.min | 969 | 911 |
| families_with_extant_copies_by_kind.cargo | 16 | 16 |
| families_with_extant_copies_by_kind.decoy | 5 | 5 |
| families_with_extant_copies_by_kind.hub | 5 | 5 |
| families_with_extant_copies_by_kind.originated | 1399 | 5212 |
| families_with_extant_copies_by_kind.planted | 100 | 96 |
| families_with_extant_copies_by_kind.plasmid | 160 | 160 |
| families_with_extant_copies_by_kind.profile_only | 9 | 9 |
| families_with_extant_copies_by_kind.root | 998 | 954 |
| families_with_extant_copies_by_kind.vehicle_chassis | 6 | 6 |
| genes_per_genome.max | 3918 | 1761 |
| genes_per_genome.mean | 1937.0 | 1415.8 |
| genes_per_genome.min | 969 | 1071 |
| genes_per_genome_by_kind.background | 1318.0 | 995.0 |
| genes_per_genome_by_kind.decoy | 6.7 | 8.3 |
| genes_per_genome_by_kind.hub | 49.3 | 48.9 |
| genes_per_genome_by_kind.planted | 380.7 | 76.0 |
| genes_per_genome_by_kind.plasmid | 64.3 | 154.8 |
| genes_per_genome_by_kind.profile_only | 29.6 | 10.5 |
| genes_per_genome_by_kind.vehicle | 88.3 | 122.2 |
| genomes | 46 | 46 |
| kinds.assembly/callerr.tsv:kind.alt_code | 6551 | 4971 |
| kinds.assembly/callerr.tsv:kind.alt_start_downstream | 30 | 30 |
| kinds.assembly/callerr.tsv:kind.alt_start_upstream | 30 | 30 |
| kinds.assembly/callerr.tsv:kind.frameshift_del | 30 | 30 |
| kinds.assembly/callerr.tsv:kind.frameshift_ins | 30 | 30 |
| kinds.assembly/callerr.tsv:kind.internal_stop | 40 | 40 |
| kinds.assembly/callerr.tsv:kind.smorf | 40 | 40 |
| kinds.assembly/edits.tsv:kind.alt_start_downstream | 30 | 30 |
| kinds.assembly/edits.tsv:kind.alt_start_upstream | 88 | 86 |
| kinds.assembly/edits.tsv:kind.frameshift_del | 30 | 30 |
| kinds.assembly/edits.tsv:kind.frameshift_ins | 30 | 30 |
| kinds.assembly/edits.tsv:kind.internal_stop | 40 | 40 |
| kinds.assembly/edits.tsv:kind.recode_codon | 26039 | 19588 |
| kinds.assembly/edits.tsv:kind.recode_overlap | 334 | 220 |
| kinds.assembly/edits.tsv:kind.recode_stop | 2189 | 1620 |
| kinds.assembly/edits.tsv:kind.smorf | 40 | 40 |
| kinds.assembly/genomes.tsv:assembly_class.complete | 18 | 18 |
| kinds.assembly/genomes.tsv:assembly_class.draft | 18 | 18 |
| kinds.assembly/genomes.tsv:assembly_class.mag | 10 | 10 |
| kinds.assembly/replicons.tsv:kind.chromosome | 46 | 46 |
| kinds.assembly/replicons.tsv:kind.plasmid | 38 | 90 |
| kinds.fusions/edits.tsv:kind.fusion | 52 | 52 |
| kinds.fusions/edits.tsv:kind.nogd | 15 | 7 |
| kinds.fusions/edits.tsv:kind.split | 14 | 14 |
| kinds.fusions/events.tsv:kind.fusion | 15 | 15 |
| kinds.fusions/events.tsv:kind.nogd | 2 | 2 |
| kinds.fusions/events.tsv:kind.split | 5 | 5 |
| kinds.fusions/fusions.tsv:kind.fusion | 52 | 52 |
| kinds.fusions/fusions.tsv:kind.split | 14 | 14 |
| kinds.fusions/replicons.tsv:kind.chromosome | 46 | 46 |
| kinds.fusions/replicons.tsv:kind.plasmid | 38 | 90 |
| kinds.mobile/insertions.tsv:kind.crispr | 30 | 30 |
| kinds.mobile/insertions.tsv:kind.prophage | 16 | 16 |
| kinds.mobile/replicons.tsv:kind.chromosome | 46 | 46 |
| kinds.mobile/replicons.tsv:kind.plasmid | 38 | 90 |
| kinds.sources/contam.tsv:kind.cross_clade | 10 | 10 |
| kinds.sources/dup.tsv:kind.exact | 108 | 71 |
| kinds.sources/dup.tsv:kind.near | 87 | 124 |
| kinds.sources/sources.tsv:kind.exact_dup | 2 | 2 |
| kinds.sources/sources.tsv:kind.near_dup | 2 | 2 |
| kinds.sources/sources.tsv:kind.original | 46 | 46 |
| phages.genes | [12 items] | [12 items] |
| phages.n | 12 | 12 |
| planted_share | 0.2153 | 0.067 |
| plasmid_genomes_by_type.pls0 | 18 | 44 |
| plasmid_genomes_by_type.pls1 | 20 | 46 |
| prophage_instances | 16 | 16 |
| seconds.S1-S11 | 238.61 | 208.53 |
| seconds.S1-S2 backbone | 37.1 | 48.59 |
| seconds.S10 sources | 6.56 | 5.27 |
| seconds.S11 emit | 24.52 | 18.36 |
| seconds.S3-S4 sequences | 86.73 | 73.8 |
| seconds.S5 replicons | 47.97 | 37.17 |
| seconds.S6 mobile | 11.47 | 9.02 |
| seconds.S7 fusions | 5.34 | 3.21 |
| seconds.S8-S9 assembly | 18.91 | 13.1 |
| stage_digests.assembly | 2f425158dc38f2302c1d02e9ebed30da | a9157aa3b65d4d67758721921f383011 |
| stage_digests.backbone | 025fdd0d27fdc65e678782216ecd21fb | 9de0eaf2a8ce45448527454c0b8cba81 |
| stage_digests.emit | a2f0e448484dab63d41281a5ef6a00c0 | 8b04777a3532b447baa4151876d711a3 |
| stage_digests.fusions | db09abf3fe79ae34a73e3ab0664a9b13 | f8424b142a12035c30660f0a5fe60de9 |
| stage_digests.mobile | 0bd7f0db5691de58ea5e1a0df1e6e570 | c2ebc655933e39f906a802a9711c96b6 |
| stage_digests.replicons | 55bda97155d75148e5cdca7b90a1d7be | 955c6dfb47f812debf149f3d045b2678 |
| stage_digests.sequences | 34ebce1324413a261acd874d7ee4fe3e | ce0b49129909e5b134dde252e64ab3c9 |
| stage_digests.sources | 770007bde8615f41d9dd5f2eb912e753 | df3b6f60dbb01a857acec61ab2be3e25 |
| tsv_rows.assembly/callerr.tsv | 6751 | 5171 |
| tsv_rows.assembly/contigs.tsv | 2431 | 1907 |
| tsv_rows.assembly/edits.tsv | 28820 | 21684 |
| tsv_rows.assembly/gene_contigs.tsv | 92488 | 68037 |
| tsv_rows.assembly/genes.tsv | 90282 | 66305 |
| tsv_rows.assembly/genomes.tsv | 46 | 46 |
| tsv_rows.assembly/repeats.tsv | 1344 | 1840 |
| tsv_rows.assembly/replicons.tsv | 84 | 136 |
| tsv_rows.backbone/events.tsv | 50898 | 116199 |
| tsv_rows.backbone/families.tsv | 2698 | 6467 |
| tsv_rows.backbone/gene_tree_nodes.tsv | 262548 | 270731 |
| tsv_rows.backbone/genes.tsv | 89102 | 65125 |
| tsv_rows.backbone/genomes.tsv | 46 | 46 |
| tsv_rows.backbone/groups.tsv | 33 | 33 |
| tsv_rows.backbone/replicons.tsv | 84 | 136 |
| tsv_rows.backbone/taxa.tsv | 82 | 82 |
| tsv_rows.backbone/taxonomy.tsv | 46 | 46 |
| tsv_rows.backbone/truth_adjacency.tsv | 89102 | 65125 |
| tsv_rows.backbone/truth_units.tsv | 3912 | 1522 |
| tsv_rows.emit/ingest.tsv | 2327 | 1896 |
| tsv_rows.emit/truth/truth_callerr.tsv | 8509 | 6459 |
| tsv_rows.emit/truth/truth_contam.tsv | 10 | 10 |
| tsv_rows.emit/truth/truth_contig.tsv | 2648 | 2124 |
| tsv_rows.emit/truth/truth_dup.tsv | 195 | 195 |
| tsv_rows.emit/truth/truth_edges.tsv | 2217 | 1825 |
| tsv_rows.emit/truth/truth_family.tsv | 2810 | 6575 |
| tsv_rows.emit/truth/truth_fusion.tsv | 69 | 68 |
| tsv_rows.emit/truth/truth_gene.tsv | 93883 | 69641 |
| tsv_rows.emit/truth/truth_genome.tsv | 62 | 62 |
| tsv_rows.emit/truth/truth_igr.tsv | 92090 | 67883 |
| tsv_rows.emit/truth/truth_module.tsv | 175 | 175 |
| tsv_rows.emit/truth/truth_rna.tsv | 1110 | 895 |
| tsv_rows.emit/truth/truth_spacer_target.tsv | 3887 | 3630 |
| tsv_rows.emit/truth/truth_traps.tsv | 213 | 473 |
| tsv_rows.fusions/edits.tsv | 81 | 73 |
| tsv_rows.fusions/events.tsv | 22 | 22 |
| tsv_rows.fusions/fusions.tsv | 66 | 66 |
| tsv_rows.fusions/genes.tsv | 90242 | 66265 |
| tsv_rows.fusions/nogd.tsv | 15 | 7 |
| tsv_rows.fusions/nogd_families.tsv | 2 | 2 |
| tsv_rows.fusions/replicons.tsv | 84 | 136 |
| tsv_rows.mobile/arrays.tsv | 30 | 30 |
| tsv_rows.mobile/cas.tsv | 40 | 40 |
| tsv_rows.mobile/genes.tsv | 90280 | 66303 |
| tsv_rows.mobile/insertions.tsv | 46 | 46 |
| tsv_rows.mobile/phage_genes.tsv | 730 | 730 |
| tsv_rows.mobile/phages.tsv | 12 | 12 |
| tsv_rows.mobile/prophages.tsv | 16 | 16 |
| tsv_rows.mobile/replicons.tsv | 84 | 136 |
| tsv_rows.mobile/spacer_targets.tsv | 3462 | 3363 |
| tsv_rows.mobile/spacers.tsv | 316 | 316 |
| tsv_rows.replicons/elements.tsv | 171872 | 133881 |
| tsv_rows.replicons/emissions.tsv | 252 | 163 |
| tsv_rows.replicons/genes.tsv | 89102 | 65125 |
| tsv_rows.replicons/junctions.tsv | 89102 | 65125 |
| tsv_rows.replicons/lineages.tsv | 117 | 117 |
| tsv_rows.replicons/links.tsv | 20 | 20 |
| tsv_rows.replicons/models.tsv | 10 | 10 |
| tsv_rows.replicons/replicons.tsv | 84 | 136 |
| tsv_rows.replicons/rna.tsv | 1078 | 852 |
| tsv_rows.sequences/domains.tsv | 4329 | 10214 |
| tsv_rows.sequences/families.tsv | 2698 | 6463 |
| tsv_rows.sequences/genes.tsv | 89102 | 65125 |
| tsv_rows.sequences/lineages.tsv | 117 | 117 |
| tsv_rows.sequences/panel.tsv | 300 | 300 |
| tsv_rows.sources/contam.tsv | 10 | 10 |
| tsv_rows.sources/dup.tsv | 195 | 195 |
| tsv_rows.sources/mutations.tsv | 3135 | 2838 |
| tsv_rows.sources/sidecars/fcs_gx_mock.tsv | 10 | 10 |
| tsv_rows.sources/source_contigs.tsv | 2315 | 1884 |
| tsv_rows.sources/sources.tsv | 50 | 50 |
| world_digest | f0262027e4ccdcbc72cd7d90f23fbb5c | 9caab69fa23f48573bd1aacdc6c8e81e |

## 7. Re-pinned digests and migrated assertions

Every stage digest of the default world (seed 20261112) changed once (*measured*, the build above). Pinned in
`tests/r2/test_synthworld_emit.py` (DIGESTS, EMIT_DIGEST, WORLD_DIGEST), `test_synthworld_nulls.py` (SIGNAL,
SIGNAL_WORLD), `test_synthworld_assembly.py` (DIGESTS, ASSEMBLY_DIGEST) and `test_synthworld_fusions.py` (S6_DIGEST):

| stage | before | after |
|---|---|---|
| backbone | 025fdd0d27fdc65e678782216ecd21fb | 9de0eaf2a8ce45448527454c0b8cba81 |
| sequences | 34ebce1324413a261acd874d7ee4fe3e | ce0b49129909e5b134dde252e64ab3c9 |
| replicons | 55bda97155d75148e5cdca7b90a1d7be | 955c6dfb47f812debf149f3d045b2678 |
| mobile | 0bd7f0db5691de58ea5e1a0df1e6e570 | c2ebc655933e39f906a802a9711c96b6 |
| fusions | db09abf3fe79ae34a73e3ab0664a9b13 | f8424b142a12035c30660f0a5fe60de9 |
| assembly | 2f425158dc38f2302c1d02e9ebed30da | a9157aa3b65d4d67758721921f383011 |
| sources | 770007bde8615f41d9dd5f2eb912e753 | df3b6f60dbb01a857acec61ab2be3e25 |
| emit | a2f0e448484dab63d41281a5ef6a00c0 | 8b04777a3532b447baa4151876d711a3 |
| world | f0262027e4ccdcbc72cd7d90f23fbb5c | 9caab69fa23f48573bd1aacdc6c8e81e |

Assertions migrated or changed (each with its reason; no other assertion was changed):

- `test_lineage_varying_rates_reach_the_genome_size_range` became `test_one_origination_axis_varies_genome_size_without_a_runaway`
  (DEC-37 (1)). Before: at the default seed the smallest chromosome <= 1,100 and the largest >= 3,600 genes, and > 10
  distinct ZOMBI2 multipliers for each rate of the two axes. After: the axes are exactly [[origination]]; > 10 distinct
  origination multipliers and duplication, transfer and loss multipliers all 1; the chromosome mean within 1,500 +-
  30%; a genome's originated genes rank with the origination multiplier its line carried since the MRCA (Spearman > 0,
  one-sided permutation p < 0.01, 2,000 permutations, fixed seed); the largest chromosome below 4,000.
- `test_taxonomy_is_cut_from_the_species_tree_at_the_preset_depths`: the crown is the living lineages' MRCA, read from
  the Newick (before: the root's end, the first speciation) (DEC-37 (3)).
- `test_dec15_a_fates_family_splits_a_unit_and_the_truth_says_so`: both bounds on within-unit adjacencies add one wrap
  pair per excised plasmid replicon (a plasmid is a unit and a circular replicon of its own; truth_adjacency includes
  the wrap pair). Forced by the plasmid count (38 replicons before, about 90 now): without it the upper bound failed,
  11,478 > 11,426.
- `test_n1_shares_s1_to_s10_and_the_emitted_contigs_with_its_signal_world` (nulls): N1's S1-S10 digests are compared
  with its own signal world's instead of the pinned SIGNAL, so the module runs at a second seed; at the default seed
  the first test still pins SIGNAL, so the assertion there is unchanged in effect. `test_signal_world_keeps_every_stage_digest`
  skips at a seed other than the default (its digests are the default seed's). `BIOBOX_SYNTHWORLD_NULL_SEED` selects
  the seed.
- New tests: `test_dec36_no_module_family_is_transferred_into_a_lineage_that_holds_the_module`,
  `test_species_tree_has_the_preset_crown_age_from_the_living_mrca` (backbone), `test_truth_module_records_the_realized_strength`
  (emit), `test_extant_crown_is_the_living_lineages_mrca_not_the_first_speciation`,
  `test_micro_preset_carries_the_redesign_of_dec36_and_dec37` and three `check_preset` refusals (r1), and
  `tests/r2/test_synthworld_preset_seeds.py` (six seeds: default, 20261113 and the first four N2 seeds). Its ranges,
  each with its source: chromosome mean 1,050-1,950 (DEC-37 (1)); planted + profile-only + decoy share <= 2 x 114 /
  1,500 (session 10's derivation from §5.1's counts); both plasmid types at >= half the seeds (DEC-37 (3)'s stated share;
  15 of 22 reached); crown age 5.0 from the MRCA within 1e-6 (DEC-37 (2), (3)); at most one grid cell per seed without a
  module whole in 9 genomes (PLAN.md §4.3 min_detectable_count; 18 of 22 seeds have all nine, the rest one fewer);
  pooled level medians each >= 1.5 x the level below and, per seed, the top level's median above the bottom's (DEC-36);
  entry exclusion replayed (0 breaches) and DEC-15 (every scored unit intact).

## 8. Null worlds

NULLS

## 9. The test that failed first

`tests/r2/test_synthworld_preset_seeds.py` as session 10 drafted it, run on the generator before this session (the base
commit fe9d12e in a scratch copy): 3 failed in 589.90 s (*measured*, wilkcomp02, 2026-11-14, `env
BIOBOX_TREE=<scratch copy> BIOBOX_POD_CI=/tmp/if2453/f0211/ci sh scripts/pod_pytest.sh F-02-11
tests/r2/test_synthworld_preset_seeds.py`; output in `failing_first_old_generator.txt`): chromosome means 1,872.7,
700.9, 404.4, 987.5, 14,089.5 and 2,112.9 with planted shares up to 60.3%; plasmid types lost at four of six seeds;
grid cells with no module whole in 9 genomes at two seeds.
