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

**Outcome: the redesign is in the generator and the preset; seven existing assertions fail under it and are
reported, not loosened (§4, §8).** DEC-36 and DEC-37 (1) and (2) hold at all 22 seeds; DEC-37 (3) holds as a stated share
(both plasmid types at 14 of 22 seeds). The stage digests of the default world changed once and are re-pinned (§5).

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
replays ZOMBI2's `genome_events.tsv` (all lineages, speciation rows included) and counts, at the default seed, 2,065
module-family transfer rows, 0 into a lineage holding the module's first family (710 into a lineage holding some other
member, which partial transfers at together < 1 allow); the same world without the rule has 10,664 rows and 10,241
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
| g_s0.2_o0.155 (burst_depth 0.02) | 1149-1795 | 22 | 15 | 11 | 0.01: 5.0, 0.06: 10.9, 0.4: 19.0 | 22 | 0.79 / 0.74 / 8 | 18 | 3.2-7.6% |
| **final (committed: g_s0.2_o0.155 with burst_depth 0.01)** | 1152-1799 | 22 | 14 | 10 | 0.01: 5.0, 0.06: 10.6, 0.4: 19.1 | 22 | 0.81 / 0.77 / 10 | 18 | 3.2-7.7% |
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
  genome's originated genes rank with the origination multiplier its line carried since the MRCA at Spearman 0.805
  (permutation p 0.0005; *measured*, the migrated test).
- **Member transfer** (DEC-36): under exclusion a member is transferred only into a lineage without the module, so the
  realized count has a floor near 4 (background transfer runs that carry a member) and saturates near 20 (the holes that
  losses make, about 0.27 per unit per unit time with background loss runs, *derived*: 0.03 + 0.24); the old rates
  0.05/0.15/0.30 sit at about 5/18/20 (c_fit3_mt_0.02,0.15,0.8 with 0.8 for 0.30). 0.01/0.06/0.4 gives pooled medians
  5.0/10.9/19.0, a factor 2 between levels, and strictly increasing level medians at 22 of 22 seeds.
- **Plasmids** (DEC-37 (3)): see §4.

## 4. What cannot be met, and why (stop and report)

1. **Both plasmid types at every seed** (DEC-37 (3)): reached at 14 of 22 seeds with the committed 0.8/0.4, and at 10 to
   18 of 22 at every pair tried from 0.2/0.25 to 1.5/0.8. A plasmid unit is cured at its own rate plus about 0.24 per
   unit per unit time by background loss runs that end in it (ZOMBI2: "a run that ends inside the group takes the whole
   group"; *derived*: 0.24 x sum over k >= 2 of P(extent >= k) = 0.24 at geometric mean 2), and a lineage that carries
   it can also die; with two lineages at the crown a type is lost early with probability about ((d + L + 0.24) / (b +
   T))^2 per type (*derived*, a branching process of carrier lineages: 0.3 death, 1.0 birth), about 0.2 at 0.8/0.4.
   Making T large enough to push that below 1/42 saturates every genome, and background duplication runs that end in
   the unit then give some genomes a second copy of a type. So survival at every seed and §5.1's "0-2 plasmids per
   genome" with some genomes carrying none (asserted at the default seed by
   `test_plasmids_are_whole_replicons_of_their_preset_size`) exclude each other in this model: at the committed
   rates the default seed has 43 and 43 carriers, no genome without a plasmid (0/6/40), so that assertion fails.
2. **N0's background sizes correlate with the signal's at r >= 0.9** (`test_n0_keeps_the_background_process_in_distribution`):
   the assertion rests on "size is set by the per-lineage rates, which ZOMBI2 draws with the tree" (its docstring), true
   of the two runaway axes DEC-37 (1) removed. With one origination axis the per-lineage spread is small against the
   stochastic spread of a genome's size, so a new draw (N0) of the same lineage rates correlates weakly: r at the
   default seed 0.53 to 0.85 over sigma 0.2 to 1.0, never 0.9, and at sigma 1.0 the default seed's N0 size ratio is
   1.20 (outside [0.9, 1.1]) and 11 of 22 seed means leave the band. No sigma meets both DEC-37 (1) and r >= 0.9 at the
   default seed. The committed preset gives r 0.81 at the default seed and 0.77 at 20261113, where the size ratio is 1.101 (just
   outside [0.9, 1.1]).
3. **The core genome at the default seed** (`test_hub_families_have_more_diverse_neighbours` needs at least 100 root
   families present in >= 40 of the 46 genomes to compute its comparator): 7 now (7 to 105 over the 22 seeds, `final`),
   243 in session 2. The default seed's
   tree is stretched 1.95x to the fixed crown, and session 2 had doubled every turnover rate for a crown of 2.56 (the
   spread DEC-37 (1) no longer asks for); a root gene's expected survival on a lineage falls as exp(-0.48 x T). With the
   spike's rates (duplication 0.02, transfer 0.05, loss 0.12) the core is 36 at the default seed (36 to 237 over 22
   seeds, 11 seeds >= 100) and the means rise to 1,471-2,196 (`h_spike`): halving the turnover is the direction, but it
   needs its own 21-seed refit of origination, plasmids and member rates, and a decision on the floor of 100.

4. **Sequence-level tests at the default seed, whose world is now twice as old** (DEC-37 (2); the sequence model is
   not this session's): `test_orthologous_igrs_are_near_identical_in_close_relatives` measured burst identity 0.864
   with burst_depth 0.02 (the burst rooted 0.1 before the present) and 0.895 with 0.01 (0.05, as before; the test asks
   0.9; *measured*, the pod_pytest runs below); `test_every_planted_motif_is_recovered_at_its_truth_coordinates_by_an_independent_scan`
   fails on one instance at 27% divergence from its emission whose best window lies elsewhere (motif clades are rooted
   at 0.5 of the crown age, now 2.5 time units instead of 1.28 at the default seed; halving `rna_motifs` clade height
   would keep them as young, a world-changing step not taken here).
5. **N1's matched-random comparator** (`test_n1_label_pairs_are_no_more_adjacent_than_matched_random_family_pairs`):
   it draws, for each planted family, root or originated families present in the same number of genomes +- 2; with
   the small core of item 3 a planted family in nearly every genome has none (`ValueError: high <= 0`, *measured* at the
   default seed). Same cause as item 3.

Proposed next steps for the PI: (a) the background turnover at the fixed crown (the spike's rates halve session 2's
and bring the core back toward the old world, `h_spike`), refitted with origination, plasmids and member rates over
the 21 seeds; (b) whether §5.1's "0-2 plasmids per genome" or DEC-37 (3)'s survival at every seed gives way (§4 item 1);
(c) replacing N0's r >= 0.9 with a test that does not assume size is set by lineage rates (item 2); (d) motif clade
height and the burst identity bound (item 4).

Every other range holds at every seed: chromosome mean (22 of 22), crown age 5.0 from the MRCA (22 of 22), entry
exclusion (0 breaches at 22 of 22), planted + profile-only + decoy share 3.2-7.6% (<= 15.2%), copies per present member
0.81 to 1.15 (session 10: 2.9 to 7.0), DEC-15 (c) and every scored unit intact (22 of 22).

## 5. The committed preset over the 22 seeds (scan `final`)

`final.tsv`: the preset as committed (sigma 0.2, origination 0.155, member transfer 0.01/0.06/0.4, plasmids 0.8/0.4,
burst_depth 0.01; *measured*, the scan command above with no overrides, `--n0`). Session 10's `base.tsv` columns
plus the crown from the MRCA, the realized levels, N0 and the entry-exclusion replay:

| seed | crown age (MRCA) / from first speciation | genomes | genes per genome mean (min-max) | chromosome mean (min-max) | planted / profile-only / decoy / vehicle / hub / plasmid / background share % | plasmid genomes pls0, pls1 (0/1/2 per genome) | modules all members in >= 9 genomes | copies per present member | DEC-15 (c) plus, kept, clustered; scored units intact | realized transfers per member, level medians | N0 background size ratio, r | module transfers into a holder |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 20261112 | 5 / 5 | 46 | 1412.9 (1106-1667) | 1263.3 (946-1547) | 5.3 / 0.7 / 0.6 / 8.7 / 3.5 / 10.6 / 70.6 | 43, 43 (0/6/40) | 22 | 1.148 | 0.460, 0.0000, 0.021; 698/698 | 9.0, 18.75, 30.67 | 1.0219, 0.8096 | 0 of 2065 |
| 20261113 | 5 / 8.435 | 46 | 1186.0 (955-1384) | 1152.1 (955-1344) | 6.6 / 0.9 / 0.3 / 11.8 / 4.2 / 2.9 / 73.5 | 0, 39 (7/39/0) | 22 | 1.049 | 0.459, 0.0019, 0.021; 591/591 | 3.8, 14.0, 21.67 | 1.101, 0.7721 | 0 of 1517 |
| 1001685937451957930 | 5 / 5 | 46 | 1550.6 (1063-1841) | 1510.6 (1023-1801) | 2.7 / 0.4 / 0.2 / 5.4 / 3.1 / 2.6 / 85.7 | 0, 46 (0/46/0) | 14 | 0.808 | 0.524, 0.0048, 0.021; 445/445 | 1.0, 2.0, 16.0 | 1.0476, 0.8729 | 0 of 716 |
| 6329801112132740259 | 5 / 5.136 | 46 | 1887.0 (1134-2212) | 1740.9 (974-2149) | 3.2 / 0.5 / 0.3 / 5.2 / 2.6 / 7.7 / 80.5 | 40, 48 (0/6/38) | 19 | 1.034 | 0.477, 0.0010, 0.032; 638/638 | 8.0, 15.0, 21.0 | 1.0444, 0.9564 | 0 of 1657 |
| 6036017085589997374 | 5 / 5 | 46 | 1591.6 (1322-1903) | 1431.6 (1157-1718) | 3.3 / 0.5 / 0.1 / 8.6 / 2.5 / 10.1 / 74.9 | 47, 43 (0/6/36) | 17 | 0.859 | 0.412, 0.0019, 0.021; 489/489 | 5.33, 11.75, 20.0 | 1.0529, 0.8688 | 0 of 1301 |
| 4383110225690451409 | 5 / 5 | 46 | 1652.2 (1295-2019) | 1485.3 (1135-1859) | 4.0 / 0.6 / 0.2 / 8.7 / 3.0 / 10.1 / 73.5 | 50, 42 (0/7/32) | 18 | 0.994 | 0.432, 0.0010, 0.065; 633/633 | 5.0, 10.67, 21.0 | 1.0421, 0.9018 | 0 of 1234 |
| 10975598825187167680 | 5 / 5 | 46 | 1338.2 (1092-1872) | 1301.6 (1085-1832) | 4.2 / 0.4 / 0.2 / 8.3 / 3.7 / 2.7 / 80.5 | 0, 42 (6/38/2) | 17 | 1.055 | 0.522, 0.0029, 0.021; 460/460 | 3.0, 4.0, 19.25 | 1.0381, 0.8405 | 0 of 945 |
| 5467273169533009601 | 5 / 5.03 | 46 | 1376.0 (890-1689) | 1230.8 (850-1529) | 5.1 / 0.5 / 0.3 / 8.2 / 3.6 / 10.5 / 71.9 | 44, 35 (1/11/34) | 25 | 1.091 | 0.426, 0.0019, 0.043; 661/661 | 7.0, 8.5, 13.33 | 0.9043, 0.5604 | 0 of 1089 |
| 4734920099712499969 | 5 / 5 | 46 | 1535.2 (1102-1783) | 1372.6 (852-1623) | 4.4 / 0.4 / 0.3 / 7.6 / 3.2 / 10.6 / 73.7 | 48, 43 (0/3/41) | 20 | 0.953 | 0.439, 0.0038, 0.032; 649/649 | 5.0, 15.0, 17.0 | 1.0421, 0.9408 | 0 of 1218 |
| 1570278518166000558 | 5 / 6.18 | 46 | 1951.3 (1239-2309) | 1799.1 (1079-2149) | 3.1 / 0.5 / 0.5 / 7.3 / 2.5 / 7.8 / 78.2 | 43, 46 (0/3/43) | 19 | 0.97 | 0.454, 0.0010, 0.054; 543/543 | 3.4, 11.33, 15.0 | 0.9921, 0.9315 | 0 of 1268 |
| 3830344756309592292 | 5 / 5 | 46 | 1512.6 (974-1848) | 1441.3 (934-1721) | 4.9 / 0.4 / 0.2 / 2.4 / 3.2 / 4.7 / 84.2 | 12, 46 (0/34/12) | 21 | 1.029 | 0.435, 0.0010, 0.043; 650/650 | 6.33, 17.8, 21.33 | 0.9869, 0.959 | 0 of 1689 |
| 9254619518473077140 | 5 / 7.657 | 46 | 1526.4 (1240-1745) | 1369.8 (1080-1585) | 4.9 / 0.5 / 0.5 / 7.4 / 3.2 / 10.2 / 73.1 | 45, 45 (0/3/42) | 22 | 1.047 | 0.525, 0.0010, 0.043; 696/696 | 11.0, 17.8, 30.33 | 1.0206, 0.8541 | 0 of 1846 |
| 4701969716356831982 | 5 / 5.087 | 46 | 1398.2 (1099-1585) | 1398.2 (1099-1585) | 4.9 / 0.5 / 0.2 / 4.1 / 3.5 / 0.0 / 86.7 | 0, 0 (46/0/0) | 19 | 0.968 | 0.428, 0.0000, 0.021; 641/641 | 4.33, 6.67, 18.33 | 1.0185, 0.7793 | 0 of 1001 |
| 14843108908847119356 | 5 / 5 | 46 | 1756.7 (1341-2240) | 1756.7 (1341-2240) | 4.0 / 0.5 / 0.2 / 5.2 / 2.8 / 0.0 / 87.3 | 0, 0 (46/0/0) | 19 | 1.006 | 0.550, 0.0029, 0.054; 639/639 | 8.0, 19.0, 26.5 | 1.0124, 0.9118 | 0 of 1793 |
| 1758546030705738202 | 5 / 5 | 46 | 1935.9 (975-2286) | 1786.4 (815-2126) | 2.8 / 0.4 / 0.1 / 6.2 / 2.5 / 7.7 / 80.2 | 43, 43 (0/9/34) | 17 | 0.892 | 0.563, 0.0000, 0.032; 510/510 | 1.5, 5.75, 9.67 | 1.0065, 0.9667 | 0 of 611 |
| 1995372076782357328 | 5 / 5 | 46 | 1755.2 (1216-2611) | 1589.9 (1056-2451) | 3.8 / 0.7 / 0.6 / 4.7 / 2.8 / 9.4 / 78.1 | 48, 46 (0/0/45) | 20 | 1.118 | 0.483, 0.0048, 0.043; 633/633 | 3.17, 7.6, 25.0 | 1.0052, 0.9679 | 0 of 1260 |
| 13514522956726874792 | 5 / 9.488 | 46 | 1436.1 (1191-1711) | 1300.5 (1088-1438) | 3.8 / 0.9 / 0.6 / 4.0 / 3.4 / 9.4 / 77.8 | 37, 45 (1/10/33) | 17 | 1.037 | 0.444, 0.0019, 0.054; 532/532 | 5.4, 6.75, 20.33 | 1.0187, 0.6126 | 0 of 1102 |
| 12621855678475765395 | 5 / 7.443 | 46 | 1567.1 (1324-1789) | 1523.6 (1284-1749) | 3.2 / 0.8 / 0.1 / 6.9 / 3.1 / 2.8 / 83.2 | 0, 50 (0/44/0) | 15 | 1.041 | 0.468, 0.0019, 0.032; 525/525 | 4.0, 7.0, 17.0 | 1.0321, 0.6107 | 0 of 962 |
| 7689157899245528466 | 5 / 5 | 46 | 1496.0 (1015-1964) | 1349.1 (895-1804) | 4.9 / 0.7 / 0.2 / 8.0 / 3.3 / 9.8 / 73.2 | 44, 37 (0/11/35) | 22 | 1.064 | 0.546, 0.0019, 0.010; 724/724 | 11.8, 16.0, 24.0 | 1.0531, 0.9158 | 0 of 2113 |
| 13748185708697958051 | 5 / 5 | 46 | 1314.8 (961-1701) | 1314.8 (961-1701) | 3.5 / 0.6 / 0.0 / 3.6 / 3.7 / 0.0 / 88.6 | 0, 0 (46/0/0) | 16 | 0.971 | 0.550, 0.0010, 0.021; 422/422 | 2.5, 14.0, 19.0 | 1.0317, 0.9501 | 0 of 1127 |
| 15390672246448918260 | 5 / 6.455 | 46 | 1627.5 (1354-1858) | 1587.5 (1298-1818) | 3.1 / 0.4 / 0.3 / 8.7 / 3.0 / 2.5 / 82.0 | 0, 46 (1/44/1) | 18 | 0.932 | 0.416, 0.0000, 0.021; 539/539 | 5.5, 10.67, 22.5 | 0.9486, 0.8139 | 0 of 1100 |
| 11714203776604047358 | 5 / 5.052 | 46 | 1427.0 (1099-1728) | 1273.1 (939-1568) | 4.5 / 0.7 / 0.2 / 11.2 / 3.4 / 10.8 / 69.2 | 45, 42 (0/7/37) | 19 | 1.072 | 0.473, 0.0029, 0.021; 563/563 | 10.0, 11.0, 18.67 | 1.0943, 0.7271 | 0 of 1473 |

Realized strength grid (medians over the modules of a cell and the 22 seeds, range in parentheses; realized together =
share of a module's transfer events carrying every member; powered = seeds at which the cell has a module with every
member in >= 9 genomes; *measured*):

| member transfer \ together | 0.3 | 0.7 | 1 | level, pooled |
|---|---|---|---|---|
| 0.01 | 4.1 (0-19); together 0.19; powered at 21/22 | 5.2 (0-23); together 0.61; powered at 22/22 | 7.0 (0-17); together 1.00; powered at 22/22 | 5.0 |
| 0.06 | 8.9 (0-22.33); together 0.19; powered at 22/22 | 10.9 (0-31.5); together 0.54; powered at 20/22 | 12.5 (0-33); together 1.00; powered at 22/22 | 10.6 |
| 0.4 | 15.8 (0-45); together 0.24; powered at 22/22 | 19.1 (0-40); together 0.62; powered at 21/22 | 23.0 (3-35); together 1.00; powered at 22/22 | 19.1 |

## 6. The default world before and after (S1-S11 at 20261112; `f0210_world_summary.py`)

Before: `2026-11-14_synthworld_seeds/default_world_summary.json` (session 10); after:
`2026-11-14_synthworld_redesign/default_world_summary_after.json` (*measured*, wilkcomp02, 2026-11-14, `python
scripts/stress/f0210_world_summary.py build /tmp/if2453/f0211/world-final 20261112`, then `summary`).

| quantity | before | after |
|---|---|---|
| arrays | 10 | 10 |
| chromosome_genes.max | 3798 | 1547 |
| chromosome_genes.mean | 1872.7 | 1263.3 |
| chromosome_genes.min | 969 | 946 |
| families_with_extant_copies_by_kind.cargo | 16 | 16 |
| families_with_extant_copies_by_kind.decoy | 5 | 5 |
| families_with_extant_copies_by_kind.hub | 5 | 5 |
| families_with_extant_copies_by_kind.originated | 1399 | 5193 |
| families_with_extant_copies_by_kind.planted | 100 | 96 |
| families_with_extant_copies_by_kind.plasmid | 160 | 160 |
| families_with_extant_copies_by_kind.profile_only | 9 | 9 |
| families_with_extant_copies_by_kind.root | 998 | 954 |
| families_with_extant_copies_by_kind.vehicle_chassis | 6 | 6 |
| genes_per_genome.max | 3918 | 1667 |
| genes_per_genome.mean | 1937.0 | 1412.9 |
| genes_per_genome.min | 969 | 1106 |
| genes_per_genome_by_kind.background | 1318.0 | 997.5 |
| genes_per_genome_by_kind.decoy | 6.7 | 8.4 |
| genes_per_genome_by_kind.hub | 49.3 | 48.9 |
| genes_per_genome_by_kind.planted | 380.7 | 75.2 |
| genes_per_genome_by_kind.plasmid | 64.3 | 149.6 |
| genes_per_genome_by_kind.profile_only | 29.6 | 10.5 |
| genes_per_genome_by_kind.vehicle | 88.3 | 122.9 |
| genomes | 46 | 46 |
| kinds.assembly/callerr.tsv:kind.alt_code | 6551 | 5028 |
| kinds.assembly/callerr.tsv:kind.alt_start_downstream | 30 | 30 |
| kinds.assembly/callerr.tsv:kind.alt_start_upstream | 30 | 30 |
| kinds.assembly/callerr.tsv:kind.frameshift_del | 30 | 30 |
| kinds.assembly/callerr.tsv:kind.frameshift_ins | 30 | 30 |
| kinds.assembly/callerr.tsv:kind.internal_stop | 40 | 40 |
| kinds.assembly/callerr.tsv:kind.smorf | 40 | 40 |
| kinds.assembly/edits.tsv:kind.alt_start_downstream | 30 | 30 |
| kinds.assembly/edits.tsv:kind.alt_start_upstream | 88 | 75 |
| kinds.assembly/edits.tsv:kind.frameshift_del | 30 | 30 |
| kinds.assembly/edits.tsv:kind.frameshift_ins | 30 | 30 |
| kinds.assembly/edits.tsv:kind.internal_stop | 40 | 40 |
| kinds.assembly/edits.tsv:kind.recode_codon | 26039 | 19877 |
| kinds.assembly/edits.tsv:kind.recode_overlap | 334 | 232 |
| kinds.assembly/edits.tsv:kind.recode_stop | 2189 | 1588 |
| kinds.assembly/edits.tsv:kind.smorf | 40 | 40 |
| kinds.assembly/genomes.tsv:assembly_class.complete | 18 | 18 |
| kinds.assembly/genomes.tsv:assembly_class.draft | 18 | 18 |
| kinds.assembly/genomes.tsv:assembly_class.mag | 10 | 10 |
| kinds.assembly/replicons.tsv:kind.chromosome | 46 | 46 |
| kinds.assembly/replicons.tsv:kind.plasmid | 38 | 86 |
| kinds.fusions/edits.tsv:kind.fusion | 52 | 52 |
| kinds.fusions/edits.tsv:kind.nogd | 15 | 7 |
| kinds.fusions/edits.tsv:kind.split | 14 | 14 |
| kinds.fusions/events.tsv:kind.fusion | 15 | 15 |
| kinds.fusions/events.tsv:kind.nogd | 2 | 2 |
| kinds.fusions/events.tsv:kind.split | 5 | 5 |
| kinds.fusions/fusions.tsv:kind.fusion | 52 | 52 |
| kinds.fusions/fusions.tsv:kind.split | 14 | 14 |
| kinds.fusions/replicons.tsv:kind.chromosome | 46 | 46 |
| kinds.fusions/replicons.tsv:kind.plasmid | 38 | 86 |
| kinds.mobile/insertions.tsv:kind.crispr | 30 | 30 |
| kinds.mobile/insertions.tsv:kind.prophage | 16 | 16 |
| kinds.mobile/replicons.tsv:kind.chromosome | 46 | 46 |
| kinds.mobile/replicons.tsv:kind.plasmid | 38 | 86 |
| kinds.sources/contam.tsv:kind.cross_clade | 10 | 10 |
| kinds.sources/dup.tsv:kind.exact | 108 | 76 |
| kinds.sources/dup.tsv:kind.near | 87 | 118 |
| kinds.sources/sources.tsv:kind.exact_dup | 2 | 2 |
| kinds.sources/sources.tsv:kind.near_dup | 2 | 2 |
| kinds.sources/sources.tsv:kind.original | 46 | 46 |
| phages.genes | [12 items] | [12 items] |
| phages.n | 12 | 12 |
| planted_share | 0.2153 | 0.0666 |
| plasmid_genomes_by_type.pls0 | 18 | 43 |
| plasmid_genomes_by_type.pls1 | 20 | 43 |
| prophage_instances | 16 | 16 |
| seconds.S1-S11 | 238.61 | 203.87 |
| seconds.S1-S2 backbone | 37.1 | 53.94 |
| seconds.S10 sources | 6.56 | 5.47 |
| seconds.S11 emit | 24.52 | 17.84 |
| seconds.S3-S4 sequences | 86.73 | 69.06 |
| seconds.S5 replicons | 47.97 | 34.59 |
| seconds.S6 mobile | 11.47 | 8.44 |
| seconds.S7 fusions | 5.34 | 2.96 |
| seconds.S8-S9 assembly | 18.91 | 11.57 |
| stage_digests.assembly | 2f425158dc38f2302c1d02e9ebed30da | e8976e16089cc6b0558475d790d70f61 |
| stage_digests.backbone | 025fdd0d27fdc65e678782216ecd21fb | f51833fa9a0a605cd90b19201f871c00 |
| stage_digests.emit | a2f0e448484dab63d41281a5ef6a00c0 | d46c15d998991020a53f35ec03de11ab |
| stage_digests.fusions | db09abf3fe79ae34a73e3ab0664a9b13 | d338a5a84ff766709b4a58d99fccf033 |
| stage_digests.mobile | 0bd7f0db5691de58ea5e1a0df1e6e570 | be536abec0056be18811b6886b89e3c1 |
| stage_digests.replicons | 55bda97155d75148e5cdca7b90a1d7be | d74442d291b8cdd57459585bfa91023a |
| stage_digests.sequences | 34ebce1324413a261acd874d7ee4fe3e | e4667ebdc498d851c2153cb6eaec1569 |
| stage_digests.sources | 770007bde8615f41d9dd5f2eb912e753 | 80a3e2e9f87f7ab8d02594f4becf74ec |
| tsv_rows.assembly/callerr.tsv | 6751 | 5228 |
| tsv_rows.assembly/contigs.tsv | 2431 | 1909 |
| tsv_rows.assembly/edits.tsv | 28820 | 21942 |
| tsv_rows.assembly/gene_contigs.tsv | 92488 | 67898 |
| tsv_rows.assembly/genes.tsv | 90282 | 66174 |
| tsv_rows.assembly/genomes.tsv | 46 | 46 |
| tsv_rows.assembly/repeats.tsv | 1344 | 1846 |
| tsv_rows.assembly/replicons.tsv | 84 | 132 |
| tsv_rows.backbone/events.tsv | 50898 | 115800 |
| tsv_rows.backbone/families.tsv | 2698 | 6448 |
| tsv_rows.backbone/gene_tree_nodes.tsv | 262548 | 270213 |
| tsv_rows.backbone/genes.tsv | 89102 | 64994 |
| tsv_rows.backbone/genomes.tsv | 46 | 46 |
| tsv_rows.backbone/groups.tsv | 33 | 33 |
| tsv_rows.backbone/replicons.tsv | 84 | 132 |
| tsv_rows.backbone/taxa.tsv | 82 | 82 |
| tsv_rows.backbone/taxonomy.tsv | 46 | 46 |
| tsv_rows.backbone/truth_adjacency.tsv | 89102 | 64994 |
| tsv_rows.backbone/truth_units.tsv | 3912 | 1522 |
| tsv_rows.emit/ingest.tsv | 2327 | 1872 |
| tsv_rows.emit/truth/truth_callerr.tsv | 8509 | 6540 |
| tsv_rows.emit/truth/truth_contam.tsv | 10 | 10 |
| tsv_rows.emit/truth/truth_contig.tsv | 2648 | 2125 |
| tsv_rows.emit/truth/truth_dup.tsv | 195 | 194 |
| tsv_rows.emit/truth/truth_edges.tsv | 2217 | 1795 |
| tsv_rows.emit/truth/truth_family.tsv | 2810 | 6556 |
| tsv_rows.emit/truth/truth_fusion.tsv | 69 | 69 |
| tsv_rows.emit/truth/truth_gene.tsv | 93883 | 69468 |
| tsv_rows.emit/truth/truth_genome.tsv | 62 | 62 |
| tsv_rows.emit/truth/truth_igr.tsv | 92090 | 67766 |
| tsv_rows.emit/truth/truth_module.tsv | 175 | 175 |
| tsv_rows.emit/truth/truth_rna.tsv | 1110 | 873 |
| tsv_rows.emit/truth/truth_spacer_target.tsv | 3887 | 3668 |
| tsv_rows.emit/truth/truth_traps.tsv | 213 | 474 |
| tsv_rows.fusions/edits.tsv | 81 | 73 |
| tsv_rows.fusions/events.tsv | 22 | 22 |
| tsv_rows.fusions/fusions.tsv | 66 | 66 |
| tsv_rows.fusions/genes.tsv | 90242 | 66134 |
| tsv_rows.fusions/nogd.tsv | 15 | 7 |
| tsv_rows.fusions/nogd_families.tsv | 2 | 2 |
| tsv_rows.fusions/replicons.tsv | 84 | 132 |
| tsv_rows.mobile/arrays.tsv | 30 | 30 |
| tsv_rows.mobile/cas.tsv | 40 | 40 |
| tsv_rows.mobile/genes.tsv | 90280 | 66172 |
| tsv_rows.mobile/insertions.tsv | 46 | 46 |
| tsv_rows.mobile/phage_genes.tsv | 730 | 730 |
| tsv_rows.mobile/phages.tsv | 12 | 12 |
| tsv_rows.mobile/prophages.tsv | 16 | 16 |
| tsv_rows.mobile/replicons.tsv | 84 | 132 |
| tsv_rows.mobile/spacer_targets.tsv | 3462 | 3363 |
| tsv_rows.mobile/spacers.tsv | 316 | 316 |
| tsv_rows.replicons/elements.tsv | 171872 | 133356 |
| tsv_rows.replicons/emissions.tsv | 252 | 164 |
| tsv_rows.replicons/genes.tsv | 89102 | 64994 |
| tsv_rows.replicons/junctions.tsv | 89102 | 64994 |
| tsv_rows.replicons/lineages.tsv | 117 | 117 |
| tsv_rows.replicons/links.tsv | 20 | 20 |
| tsv_rows.replicons/models.tsv | 10 | 10 |
| tsv_rows.replicons/replicons.tsv | 84 | 132 |
| tsv_rows.replicons/rna.tsv | 1078 | 834 |
| tsv_rows.sequences/domains.tsv | 4329 | 10191 |
| tsv_rows.sequences/families.tsv | 2698 | 6444 |
| tsv_rows.sequences/genes.tsv | 89102 | 64994 |
| tsv_rows.sequences/lineages.tsv | 117 | 117 |
| tsv_rows.sequences/panel.tsv | 300 | 300 |
| tsv_rows.sources/contam.tsv | 10 | 10 |
| tsv_rows.sources/dup.tsv | 195 | 194 |
| tsv_rows.sources/mutations.tsv | 3135 | 2842 |
| tsv_rows.sources/sidecars/fcs_gx_mock.tsv | 10 | 10 |
| tsv_rows.sources/source_contigs.tsv | 2315 | 1860 |
| tsv_rows.sources/sources.tsv | 50 | 50 |
| world_digest | f0262027e4ccdcbc72cd7d90f23fbb5c | 31085c4d1be9fa94df9238f785a375d6 |

## 7. Re-pinned digests and migrated assertions

Every stage digest of the default world (seed 20261112) changed once (*measured*, the build above). Pinned in
`tests/r2/test_synthworld_emit.py` (DIGESTS, EMIT_DIGEST, WORLD_DIGEST), `test_synthworld_nulls.py` (SIGNAL,
SIGNAL_WORLD), `test_synthworld_assembly.py` (DIGESTS, ASSEMBLY_DIGEST) and `test_synthworld_fusions.py` (S6_DIGEST):

| stage | before | after |
|---|---|---|
| backbone | 025fdd0d27fdc65e678782216ecd21fb | f51833fa9a0a605cd90b19201f871c00 |
| sequences | 34ebce1324413a261acd874d7ee4fe3e | e4667ebdc498d851c2153cb6eaec1569 |
| replicons | 55bda97155d75148e5cdca7b90a1d7be | d74442d291b8cdd57459585bfa91023a |
| mobile | 0bd7f0db5691de58ea5e1a0df1e6e570 | be536abec0056be18811b6886b89e3c1 |
| fusions | db09abf3fe79ae34a73e3ab0664a9b13 | d338a5a84ff766709b4a58d99fccf033 |
| assembly | 2f425158dc38f2302c1d02e9ebed30da | e8976e16089cc6b0558475d790d70f61 |
| sources | 770007bde8615f41d9dd5f2eb912e753 | 80a3e2e9f87f7ab8d02594f4becf74ec |
| emit | a2f0e448484dab63d41281a5ef6a00c0 | d46c15d998991020a53f35ec03de11ab |
| world | f0262027e4ccdcbc72cd7d90f23fbb5c | 31085c4d1be9fa94df9238f785a375d6 |

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
  14 of 22 reached); crown age 5.0 from the MRCA within 1e-6 (DEC-37 (2), (3)); at most one grid cell per seed without a
  module whole in 9 genomes (PLAN.md §4.3 min_detectable_count; 18 of 22 seeds have all nine, the rest one fewer);
  pooled level medians each >= 1.5 x the level below and, per seed, the top level's median above the bottom's (DEC-36);
  entry exclusion replayed (0 breaches) and DEC-15 (every scored unit intact).

## 8. Null worlds

`tests/r2/test_synthworld_nulls.py` at the default seed and at 20261113 (`BIOBOX_SYNTHWORLD_NULL_SEED`), each in a
detached unit at 2 CPUs in the image from the pod_pytest mirror (*measured*, wilkcomp02, 2026-11-14, `python -m pytest
-q -s tests/r2/test_synthworld_nulls.py`): **default seed 2 failed, 12 passed in 1,388 s; 20261113 1 failed, 12
passed, 1 skipped (the pinned-digest test) in 1,371 s.** No assertion was loosened; no threshold there depended on a
count that changed (the planted-pair counts in docstrings are printed, not asserted).

| test | default seed | 20261113 |
|---|---|---|
| N0 background size: ratio in [0.9, 1.1], r >= 0.9 | fails: ratio 1.022, r 0.810 | fails: ratio 1.101, r 0.772 |
| N1 label pairs against matched random pairs | fails: no root or originated family within +-2 genomes of a planted family's prevalence (§4 item 5) | passes |
| every other test (byte identity in two processes, N0 keeps the tree and per-genome draws, no planted signal, N0-strict, another seed, truth_null, naive statistic miscalibrated, N1 shares S1-S10, N1 relabel bounds, times) | pass | pass |

N1 relabels **11 of 100** module members at the default seed (4 of 25 modules, 29 evaluated label pairs) and 17 of 100
at 20261113 (7 modules, 41 pairs); session 9 relabeled 27 or 28. Fewer non-module families now match a module
family's profile within 4 genomes and its transfers within 25% (entry exclusion keeps the members' transfer counts at 3
to 48, *measured*, while the root families' profiles spread with the doubled time). The naive co-occurrence statistic
stays miscalibrated on both N0 worlds (21.0x and 20.5x the calibrated share at p < 1e-3, default seed). Fixes in this
file: N1's relabel check reads a module family present in no genome as an empty profile (4 of the 100 are absent at
the default seed; the test had read `pres[a]` and stopped with KeyError 'F001035').

## 9. The test that failed first

`tests/r2/test_synthworld_preset_seeds.py` as session 10 drafted it, run on the generator before this session (the base
commit fe9d12e in a scratch copy): 3 failed in 589.90 s (*measured*, wilkcomp02, 2026-11-14, `env
BIOBOX_TREE=<scratch copy> BIOBOX_POD_CI=/tmp/if2453/f0211/ci sh scripts/pod_pytest.sh F-02-11
tests/r2/test_synthworld_preset_seeds.py`; output in `failing_first_old_generator.txt`): chromosome means 1,872.7,
700.9, 404.4, 987.5, 14,089.5 and 2,112.9 with planted shares up to 60.3%; plasmid types lost at four of six seeds;
grid cells with no module whole in 9 genomes at two seeds.

## 10. Test runs on the committed tree

Every synthworld test outside the null module and `tests/r0` together: **4 failed, 500 passed, 6 skipped in 2,141 s**
(*measured*, wilkcomp02, 2026-11-14, `env BIOBOX_TREE=<worktree> sh scripts/pod_pytest.sh F-02-11 -s tests/r0
tests/r1/test_synthworld_*.py tests/r2/test_synthworld_{backbone,preset_seeds,sequences,sequences_pfam,intergenic,
rna_rfam,mobile,fusions,assembly,emit}.py`). The four: `test_plasmids_are_whole_replicons_of_their_preset_size` (§4
item 1), `test_hub_families_have_more_diverse_neighbours` (item 3), `test_orthologous_igrs_are_near_identical_in_close_relatives`
and `test_every_planted_motif_is_recovered_at_its_truth_coordinates_by_an_independent_scan` (item 4). The per-seed test
passes (both plasmid types at 4 of its 6 seeds; one unpowered cell at 4383110225690451409; pooled levels 5.6, 12.9,
20.5), as do the new DEC-36 and crown tests and every re-pinned digest. The null module: §8.

## 11. Follow-up for the PI: the crown at the default seed's old age (measurement only)

**Hypothesis tested** (the orchestrator's, not adopted): four of the five causes in §4 follow from the default world
being 1.95x older than the one every rate since session 2 was fitted in. **Variant**: `species.crown_age` 2.5632 measured
from the living lineages' MRCA, the definition of §1. At the default seed the MRCA is the first speciation, so this is
its old age under both definitions, the rescale factor there is 1 and its tree is the one sessions 1 to 10 used. Also
`burst_depth` back to 0.02 (0.051 time units, as before), the background turnover unchanged (duplication 0.04,
transfer 0.1, loss 0.24; this session never changed it), motif clade height unchanged as a fraction (0.5 of the crown,
so 1.28 time units at the default seed, as before). Refitted only origination, plasmids and member transfer. Set by
`--set` overrides in the scans and, for the tests, in a scratch copy of the branch (`/tmp/f0211-local/variant`,
mirrored to `/tmp/if2453/f0211/ci/F-02-11` on the pod); `presets/micro.yaml`, the generator and the tests on the branch
are unchanged and nothing is re-pinned.

All *measured*, wilkcomp02, 2026-11-14 19:00-21:00 CST, image `bio-box-test_faec1f4.sif`, units `f0211-*`, at most
8 threads: scans `python scripts/stress/f0210_seed_scan.py scan /tmp/if2453/f0211/scan<V|W|X|VF> <label> --jobs 7 --seeds
default,20261113,n2 --set species.crown_age=2.5632 --set species.burst_depth=0.02 [...]` (overrides per label in
`<label>.json`).

**Fits** (22 seeds each; columns as named; "core" = root families present in >= 40 of the 46 genomes):

| scan | chromosome mean (in band) | default mean | planted share % | copies per present member | exclusion breaches | both plasmid types | a genome > 2 plasmids | genomes without a plasmid at default | all cells powered | levels pooled (increasing) | core >= 40 genomes: default, range (>= 100) | N0 ratio, r at default |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| v_o0.25_s0.3 | 1205-2120 (20/22) | 1205 | 3.6-7.4 | 0.96-1.13 | 0 | 18/22 | 9 | 1 | 22/22 | 4.0, 7.6, 11.1 (22/22) | 29, 29-200 (9) | - |
| v_o0.3_s0.3 | 1271-2313 (17/22) | 1293 | 3.2-7.2 | 0.89-1.12 | 0 | 17/22 | 9 | 0 | 22/22 | 4.0, 7.7, 11.0 (20/22) | 29, 29-206 (10) | - |
| v_o0.35_s0.3 | 1325-2564 (15/22) | 1367 | 3.6-6.2 | 0.92-1.14 | 0 | 17/22 | 11 | 0 | 21/22 | 4.0, 7.8, 11.8 (20/22) | 30, 30-216 (12) | - |
| v_o0.3_s0.5 | 1125-2827 (17/22) | 1222 | 3.0-7.6 | 0.96-1.20 | 0 | 18/22 | 8 | 0 | 21/22 | 3.5, 7.3, 10.8 (20/22) | 26, 26-238 (7) | - |
| v_mt_0.01,0.06,0.4 | 1234-2399 (17/22) | 1297 | 3.3-6.8 | 0.89-1.14 | 0 | 21/22 | 12 | 0 | 22/22 | 3.0, 6.1, 11.0 (20/22) | 42, 42-212 (10) | - |
| v_mt_0.03,0.15,1.2 | 1242-2280 (17/22) | 1353 | 3.3-7.4 | 0.94-1.15 | 0 | 16/22 | 10 | 0 | 21/22 | 4.3, 9.0, 11.1 (17/22) | 37, 37-201 (8) | - |
| v_pls_0.5-0.2 | 1224-2289 (16/22) | 1319 | 3.9-7.1 | 0.87-1.16 | 0 | 20/22 | 14 | 0 | 22/22 | 4.0, 8.0, 11.0 (19/22) | 36, 36-222 (9) | - |
| v_pls_1.2-0.6 | 1248-2272 (17/22) | 1313 | 3.2-7.9 | 0.85-1.16 | 0 | 15/22 | 9 | 1 | 22/22 | 3.7, 7.9, 11.0 (17/22) | 38, 38-258 (9) | - |
| v_pls_0.3-0.3 | 1268-2348 (17/22) | 1324 | 3.7-7.2 | 0.92-1.11 | 0 | 18/22 | 6 | 1 | 22/22 | 4.0, 7.5, 11.0 (19/22) | 31, 31-203 (9) | - |
| w_o0.2_s0.3 | 1228-1879 (22/22) | 1228 | 4.5-7.1 | 0.93-1.12 | 0 | 20/22 | 11 | 0 | 21/22 | 3.0, 6.1, 11.0 (20/22) | 38, 38-199 (10) | - |
| w_o0.15_s0.3 | 1132-1645 (22/22) | 1132 | 3.8-7.4 | 0.92-1.10 | 0 | 18/22 | 10 | 0 | 20/22 | 3.0, 5.5, 10.9 (19/22) | 40, 38-219 (10) | - |
| w_o0.2_s0.2 | 1222-1734 (22/22) | 1319 | 4.7-7.0 | 0.93-1.15 | 0 | 19/22 | 14 | 0 | 22/22 | 3.0, 6.0, 10.1 (20/22) | 41, 38-213 (9) | - |
| w_o0.15_s0.2 | 1187-1521 (22/22) | 1200 | 4.9-8.1 | 0.77-1.13 | 0 | 17/22 | 7 | 0 | 19/22 | 3.0, 6.0, 11.0 (19/22) | 47, 47-239 (9) | - |

Chosen: origination 0.2 x root families at sigma 0.3 (the only candidate with the mean in band at 22 of 22 and the most
plasmid survival, 20 of 22), member transfer 0.01/0.06/0.4 (the committed levels; pooled 3.0/6.1/11.0, each >= 1.8x
the one below; 0.02/0.1/0.8 gave 4.0/7.6/11.1 with the top step 1.46x), plasmids 0.8/0.4 (the committed pair; survival
15 to 21 of 22 over the pairs tried at this crown). The fitted variant run again with `--n0` (`variant.tsv`) gives the
same table as `w_o0.2_s0.3`.

**The variant over the 22 seeds** (`variant.tsv`): chromosome mean 1,228-1,879 (22 of 22 in band; default 1,228);
planted + profile-only + decoy share 4.5-7.1%; copies per present member 0.93-1.12; exclusion breaches 0; both plasmid
types at 20 of 22; a genome with more than 2 plasmids at 11; genomes without a plasmid at the default seed 0; all nine
cells powered at 21 of 22; realized levels pooled 3.0, 6.1, 11.0, strictly increasing at 20 of 22 seeds; core at the
default seed 38 (38 to 199 over the seeds, >= 100 at 10); N0 at the default seed ratio 0.998, r 0.846 (r 0.17 at
20261113, >= 0.9 at 10 of 22).

**Why the core stays small** (default seed only, batch X, `x_*`):

| scan | chromosome mean (in band) | default mean | planted share % | copies per present member | exclusion breaches | both plasmid types | a genome > 2 plasmids | genomes without a plasmid at default | all cells powered | levels pooled (increasing) | core >= 40 genomes: default, range (>= 100) | N0 ratio, r at default |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| x_oldaxes | 1392-1392 (1/1) | 1392 | 7.2-7.2 | 1.06-1.06 | 0 | 1/1 | 0 | 2 | 1/1 | 2.0, 7.8, 9.0 (1/1) | 204, 204-204 (1) | - |
| x_noaxes | 1646-1646 (1/1) | 1646 | 5.1-5.1 | 1.07-1.07 | 0 | 1/1 | 1 | 0 | 1/1 | 5.0, 7.5, 13.3 (1/1) | 35, 35-35 (0) | - |
| x_oldaxes_noexcl | 1873-1873 (1/1) | 1873 | 21.5-21.5 | 3.76-3.76 | 7708 | 1/1 | 0 | 17 | 1/1 | 11.0, 39.0, 127.7 (1/1) | 243, 243-243 (1) | - |

`x_oldaxes_noexcl` is the old preset in every backbone value (two axes, origination 0.1, no exclusion, old member and
plasmid rates) run through this session's generator at crown 2.5632: it reproduces session 10's default world (mean
1,873, planted share 21.5%, core 243), so the variant's machinery is sound. With DEC-36 and the old axes the core is
204; with the one origination axis it is 26 to 47 at crown 2.5632 over every draw above, and 35 with no axis. **The
core is set by the two lineage axes that DEC-37 (1) removed, not by the crown**: a mean-corrected log-normal multiplier
at sigma 0.8 drifted over the splits leaves most lineages well below 1 on loss (the median falls as exp(-0.32 k) after k
splits, *derived*), and those lineages keep the root families, so many are in >= 40 genomes.

**The six assertions at the default seed** (variant: the six tests themselves on the scratch copy, `pod_pytest.sh` with
`BIOBOX_TREE=/tmp/f0211-local/variant BIOBOX_POD_CI=/tmp/if2453/f0211/ci`, "4 failed, 34 passed in 595.27s" for the
backbone and intergenic modules, `variant_backbone_intergenic_tests.txt`; the null module in a detached unit, "2 failed,
12 passed in 796.98s" with the pinned-digest test failing as expected and ignored, `variant_nulls_default_seed.txt`):

| assertion (default seed unless named) | bound | crown 5.0 (committed) | crown 2.5632 (variant) | crown 5.0, halved turnover (`h_spike`, not refitted) |
|---|---|---|---|---|
| plasmids: some genome without one, at most 2 | >= 1 genome with none | fails: 0 (0/6/40) | fails: 0 | fails: 0 (scan) |
| hub test: root families in >= 40 genomes | >= 100 | fails: 7 | fails: 38 | fails: 36 |
| burst orthologous IGR identity | >= 0.9 | fails: 0.895 | fails: 0.882 | not measured |
| planted motif found at its coordinates | best window or divergence > 0.3 | fails (one instance at 27%) | fails (one instance at 8%, another window 1 mismatch away) | not measured |
| N0 background size ratio, r | [0.9, 1.1], >= 0.9 | fails: 1.022, 0.810 | fails: 0.998, 0.846 | not measured |
| N1 matched random comparator | candidates exist | fails (none within +-2 genomes) | **passes** | not measured |
| both plasmid types at every seed (22) | 22 | 14 | 20 | 17 |
| chromosome mean within 30% at every seed (22) | 22 | 22 | 22 | 16 (origination not refitted) |

What each would still need decided: **crown 5.0**: (a)-(d) of §4. **Crown 2.5632**: an amendment of DEC-37 (2) (a crown
of 2.56 is the shortest of the 21, not their median), and still (b) plasmids, (c) N0's r, the core (the axes, below)
and the two sequence-level tests. **Halved turnover at 5.0**: a refit of origination and the rest, and still the core
(36 at the default seed), (b) and (c).

**Verdict: the hypothesis is wrong for the core genome and for the two sequence-level tests, and right only for N1's
comparator and, in part, plasmid survival** (20 of 22 at the shorter crown against 14 of 22, *measured*; the early
loss of a type is less likely when the first splits come sooner in absolute time). The core comes from the lineage
axes (above). The burst identity is 0.882 at the old burst age and the old crown, below 0.895 at 5.0, so it is not the
burst's age; my best account, not tested, is the same axes: the old shared multiplier on duplication and transfer was
below 1 in most lineages, so most junctions were single-copy, whereas now paralogous copies of a junction (same left
family, right family and type, a different copy) enter the comparison (the low pairs, 0.67 to 0.81, look like such
pairs). The motif failure is one instance with a near-identical window elsewhere in the same region, which a paralogous
copy would also explain. N0's r stays below 0.9 at either crown (0.846), as §4 item 2 says.

## 12. The background rates refitted for one lineage axis (question (a))

The orchestrator withdrew the crown hypothesis after §11: the core came from the old lineage axes, and DEC-37 (1)
removed those axes but left session 2's base rates (duplication 0.04, transfer 0.1, loss 0.24), chosen for a world in
which most lineages ran far below them. Those base rates are free preset values; this section refits them at the
decided crown (5.0) with the single origination axis. All *measured*, wilkcomp02, 2026-11-14 20:30 to 2026-11-15
03:00 CST, at most 8 threads, `f0210_seed_scan.py scan` with the overrides in each `<label>.json`.

**1. What the old world's lineages saw** (`f0210_seed_scan.py effrates`, new: each extant genome's ZOMBI2 multiplier
averaged over time along its line of descent from the rate onset to the present; the old preset run through this
tree's generator and kept, `y_oldworld`, which reproduces session 10's table, e.g. 1,873 at the default seed and 11,192
at 13748185708697958051; per seed in `y_oldworld_effrates.tsv`):

| multiplier (time-averaged per genome) | q10 | q25 | median | q75 | mean |
|---|---|---|---|---|---|
| duplication = transfer = origination, 22 seeds (1,012 genomes) | 0.060 | 0.157 | 0.361 | 0.893 | 0.867 |
| loss, 22 seeds | 0.149 | 0.465 | 0.849 | 1.654 | 1.573 |
| duplication = transfer = origination, default seed | 0.026 | 0.035 | 0.060 | 0.061 | 0.217 |
| loss, default seed | 0.287 | 0.433 | 0.451 | 0.679 | 0.577 |

So a typical old lineage saw duplication 0.04 x 0.36 = 0.014, transfer 0.036 and loss 0.24 x 0.85 = 0.20 over the
seeds, and at the default seed duplication 0.0024, transfer 0.006 and loss 0.108 over its crown of 2.56 (*derived*).
(This corrects §11's reading that most lineages drifted to a loss multiplier well below 1: it is the shared
duplication and transfer multiplier that was low; the loss median over the seeds is 0.85.) The old world had a core of
at least 100 root families at only 8 of the 22 seeds (0 to 743; 243 at the default seed, `y_oldworld`).

**2. The fits** (22 seeds, crown 5.0, one origination axis at sigma 0.2; columns as in §11 plus N1's comparator
(planted families present with no root or originated family within +-2 genomes of their prevalence, read from the
backbone) and the burst's paralogous junction keys (per burst genome, chromosome adjacency keys of left family, right
family and strands occurring more than once)):

| scan | chromosome mean (in band) | default mean | planted share % | copies per present member | exclusion breaches | both plasmid types | a genome > 2 plasmids | genomes without a plasmid at default | all cells powered | levels pooled (increasing) | core >= 40 genomes: default, range (>= 100) | N0 ratio, r at default | N1 unmatched planted at default | burst multicopy junction keys at default |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| y_f0.2_o0.065 | 1354-1783 (22/22) | 1436 | 4.2-6.1 | 0.92-1.10 | 0 | 19/22 | 5 | 0 | 21/22 | 2.1, 6.7, 10.1 (20/22) | 203, 203-464 (22) | - | 0 | 159 of 5157 |
| y_f0.2_o0.08 | 1452-1812 (22/22) | 1558 | 4.2-6.4 | 0.90-1.10 | 0 | 21/22 | 6 | 0 | 22/22 | 2.0, 7.0, 10.0 (20/22) | 247, 216-488 (22) | - | 0 | 287 of 5614 |
| y_f0.3_o0.077 | 1340-1856 (22/22) | 1532 | 4.0-5.9 | 0.92-1.14 | 0 | 14/22 | 5 | 0 | 22/22 | 3.0, 8.0, 11.8 (19/22) | 129, 117-365 (22) | - | 0 | 384 of 5699 |
| y_f0.3_o0.09 | 1408-1882 (22/22) | 1472 | 3.1-7.1 | 0.83-1.08 | 0 | 16/22 | 4 | 1 | 21/22 | 2.7, 7.8, 11.0 (20/22) | 111, 111-352 (22) | - | 0 | 311 of 5575 |
| y_f0.4_o0.088 | 1289-1818 (22/22) | 1422 | 3.9-6.5 | 0.79-1.09 | 0 | 19/22 | 7 | 0 | 20/22 | 3.0, 8.7, 12.0 (18/22) | 69, 69-263 (18) | - | 0 | 284 of 5125 |
| y_f0.4_o0.1 | 1396-1849 (22/22) | 1408 | 3.1-5.9 | 0.71-1.10 | 0 | 13/22 | 8 | 0 | 20/22 | 3.0, 8.3, 12.0 (22/22) | 74, 74-264 (17) | - | 0 | 201 of 5371 |
| y_committed_default | 1263-1263 (1/1) | 1263 | 6.7-6.7 | 1.15-1.15 | 0 | 1/1 | 0 | 0 | 1/1 | 9.0, 18.8, 30.7 (1/1) | 7, 7-7 (0) | - | 13 | 301 of 4392 |
| y_oldworld | 404-15445 (6/22) | 1873 | 4.2-70.7 | 2.90-7.03 | 19131 | 5/22 | 0 | 17 | 16/22 | 15.2, 75.0, 131.0 (20/22) | 243, 0-743 (8) | - | 0 | 618 of 6911 |
| z_mt_0.01,0.03,0.4 | 1367-1760 (22/22) | 1448 | 4.4-6.5 | 0.88-1.30 | 0 | 15/22 | 4 | 0 | 21/22 | 2.3, 4.0, 10.3 (21/22) | 253, 228-513 (22) | - | 0 | 203 of 5301 |
| z_mt_0.005,0.03,0.4 | 1396-1775 (22/22) | 1482 | 4.1-6.2 | 0.87-1.30 | 0 | 12/22 | 3 | 1 | 21/22 | 2.0, 4.5, 10.0 (21/22) | 238, 236-486 (22) | - | 0 | 186 of 5359 |
| z_mt_0.01,0.035,0.8 | 1361-1737 (22/22) | 1526 | 4.3-6.1 | 0.86-1.19 | 0 | 14/22 | 1 | 0 | 22/22 | 2.3, 4.6, 10.2 (21/22) | 252, 235-490 (22) | - | 0 | 243 of 5568 |
| z_n0 | 1354-1783 (22/22) | 1436 | 4.2-6.1 | 0.92-1.10 | 0 | 19/22 | 5 | 0 | 21/22 | 2.1, 6.7, 10.1 (20/22) | 203, 203-464 (22) | 1.0417, 0.3748 | 0 | 159 of 5157 |
| b_t0.05_o0.01 | 1559-1946 (22/22) | 1778 | 3.8-5.5 | 0.85-1.09 | 0 | 15/22 | 2 | 0 | 22/22 | 3.8, 6.0, 9.7 (19/22) | 420, 381-590 (22) | - | 0 | 371 of 6559 |
| b_t0.05_o0.02 | 1678-2060 (21/22) | 1860 | 3.5-4.8 | 0.89-1.09 | 0 | 15/22 | 5 | 0 | 22/22 | 3.5, 5.6, 10.0 (19/22) | 415, 396-611 (22) | - | 0 | 383 of 6940 |
| b_t0.045_o0.02 | 1514-1900 (22/22) | 1646 | 3.6-5.8 | 0.79-1.09 | 0 | 20/22 | 2 | 0 | 22/22 | 3.0, 5.4, 10.0 (19/22) | 342, 342-560 (22) | - | 0 | 405 of 6312 |
| b_t0.045_o0.03 | 1585-1988 (21/22) | 1815 | 3.6-5.4 | 0.82-1.07 | 0 | 17/22 | 2 | 1 | 21/22 | 3.5, 5.0, 10.3 (17/22) | 365, 350-538 (22) | - | 0 | 376 of 6687 |
| b_t0.05_l0.06_o0.01 | 1415-1847 (22/22) | 1631 | 4.0-5.9 | 0.80-1.11 | 0 | 18/22 | 2 | 0 | 21/22 | 3.1, 5.0, 11.3 (17/22) | 252, 252-477 (22) | - | 0 | 305 of 6139 |
| b_t0.05_l0.07_o0.02 | 1398-1799 (22/22) | 1555 | 3.7-8.1 | 0.86-1.35 | 0 | 17/22 | 4 | 0 | 22/22 | 3.6, 6.0, 11.0 (18/22) | 222, 201-442 (22) | - | 0 | 282 of 5649 |
| b_t0.05_l0.08_o0.03 | 1242-1593 (22/22) | 1480 | 3.7-6.3 | 0.80-1.09 | 0 | 19/22 | 3 | 0 | 21/22 | 4.0, 5.6, 12.0 (19/22) | 174, 143-368 (22) | - | 0 | 224 of 5303 |
| c_t0.05_l0.07_o0.02_mt0.06 | 1380-1715 (22/22) | 1513 | 3.3-7.7 | 0.89-1.34 | 0 | 19/22 | 5 | 0 | 22/22 | 3.5, 8.0, 13.0 (20/22) | 193, 180-465 (22) | 1.0446, 0.5348 | 0 | 263 of 5724 |
| final4 | 1367-1760 (22/22) | 1448 | 4.4-6.5 | 0.88-1.30 | 0 | 15/22 | 4 | 0 | 21/22 | 2.3, 4.0, 10.3 (21/22) | 253, 228-513 (22) | 1.0233, 0.2826 | 0 | 203 of 5301 |
| z_pls_0.15-0.3 | 1392-1729 (22/22) | 1405 | 4.1-7.2 | 0.91-1.26 | 0 | 15/22 | 1 | 10 | 22/22 | 2.1, 5.0, 10.0 (19/22) | 211, 206-496 (22) | - | 0 | 208 of 5547 |
| z_pls_0.1-0.3 | 1389-1776 (22/22) | 1481 | 4.5-6.7 | 0.94-1.39 | 0 | 17/22 | 2 | 4 | 21/22 | 2.4, 4.0, 10.1 (18/22) | 248, 226-494 (22) | - | 0 | 193 of 5675 |
| z_pls_3.0-0.1 | 1319-1690 (22/22) | 1496 | 4.2-6.0 | 0.95-1.13 | 0 | 22/22 | 7 | 0 | 21/22 | 2.0, 4.2, 11.0 (18/22) | 240, 233-486 (22) | - | 0 | 272 of 5372 |

Labels: `y_f<F>_o<O>`: duplication, transfer and loss at F x session 2's, hub multiplier 20 / F, origination O;
`y_committed_default`: the preset before this section, default seed only; `z_mt_*`: F 0.2, O 0.065 with member transfer
levels as named; `b_t<T>_o<O>`, `b_t<T>_l<L>_o<O>`: duplication 0.008, transfer T, loss L (0.048 unless named),
member transfer 0.01/0.03/0.4; `c_*`: the same with 0.01/0.06/0.4; `final4`: the preset as committed; `z_pls_*`: it with
plasmid rates as named.

**Full worlds at the default seed** (the sequence-level tests read S3-S5, which the scan does not build; `pod_pytest.sh`
on scratch copies or the branch, outputs in this directory):

| background (duplication, transfer, loss; origination; member transfer) | burst IGR identity (>= 0.9) | equal-length cross-phylum IGR pairs (>= 200) | planted motif found | other failures | file |
|---|---|---|---|---|---|
| 0.04, 0.1, 0.24; 0.155; 0.01/0.06/0.4 (before) | 0.895 | 263 | no | hub core, plasmids | §10 |
| 0.008, 0.02, 0.048; 0.065; 0.01/0.06/0.4 | 0.934 | 118 | yes | plasmids, levels 1.40x over six seeds | `tests_f0.2.txt` |
| 0.012, 0.03, 0.072; 0.077; 0.01/0.06/0.4 | 0.934 | 169 | yes | (intergenic only) | `intergenic_f0.3.txt` |
| 0.008, 0.035, 0.048; 0.04; 0.01/0.06/0.4 | 0.943 | 188 | yes | (intergenic only) | `intergenic_t0.035.txt` |
| 0.008, 0.05, 0.048; 0.01; 0.01/0.06/0.4 | 0.957 | 262 | no (an exact copy of the emission 257 bases away) | (intergenic only) | `intergenic_t0.05.txt` |
| 0.008, 0.05, 0.07; 0.02; 0.01/0.06/0.4 | 0.947 | 236 | no (an exact copy 511 bases away) | junction lengths (codir_tight q90 33 against 38, tolerance 4.56), plasmids, levels 1.40x | `tests_t0.05_l0.07.txt` |
| **0.008, 0.02, 0.048; 0.065; 0.01/0.03/0.4 (committed)** | 0.931 | 150 | yes | see §13 | `tests_final_refit.txt` |

The paralogous-junction account of §11 is refuted by counting: the burst's chromosomes carry adjacency keys occurring
more than once in 301 of 4,392 keys at the preset before this section (6.9%), 618 of 6,911 in the old world (8.9%), 159
of 5,157 at factor 0.2 (3.1%); the old world had more repeated junctions than the failing one, and burst identity rose
to 0.93 at every lower turnover. The cross-phylum pair count instead follows the background transfer (118, 169, 188,
262 at 0.02, 0.03, 0.035, 0.05): a vertically inherited junction's IGR changes length over the 4.6 time units since the
phyla split, so equal-length pairs across phyla come mostly from junctions carried between them by transfer. The
planted-motif test failed at both settings with transfer 0.05 because another, exact copy of the same emission sat in
the same intergenic region (two samples; a mechanism is not established).

**Chosen**: duplication 0.008, transfer 0.02, loss 0.048 (session 2's x 0.2), origination 0.065 x root families, hub
multiplier 100 (hub transfer 2.0 kept), member transfer 0.01/0.03/0.4 (refitted because the 0.06 and 0.4 levels fell
within 1.4x over the per-seed test's six seeds at this background). Transfer 0.05 would meet the cross-phylum count but
failed the motif test twice and the junction-length test once, and put the levels within 1.4x again; it is not taken.

**The committed preset over the 22 seeds** (`final4.tsv`, with `--n0`):

| seed | crown age (MRCA) / from first speciation | genomes | genes per genome mean (min-max) | chromosome mean (min-max) | planted / profile-only / decoy / vehicle / hub / plasmid / background share % | plasmid genomes pls0, pls1 (0/1/2 per genome) | modules all members in >= 9 genomes | copies per present member | DEC-15 (c) plus, kept, clustered; scored units intact | realized transfers per member, level medians | N0 background size ratio, r | module transfers into a holder |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 20261112 | 5 / 5 | 46 | 1569.2 (1363-1741) | 1447.5 (1243-1606) | 4.6 / 0.6 / 0.1 / 6.7 / 3.1 / 7.8 / 77.1 | 45, 5 (0/42/4) | 25 | 1.034 | 0.460, 0.0000, 0.021; 763/763 | 2.0, 4.0, 13.0 | 1.0233, 0.2826 | 0 of 759 |
| 20261113 | 5 / 8.435 | 46 | 1544.2 (1293-1684) | 1422.5 (1265-1524) | 4.8 / 0.5 / 0.2 / 8.1 / 3.2 / 7.9 / 75.3 | 32, 44 (1/16/27) | 21 | 1.019 | 0.459, 0.0019, 0.021; 716/716 | 2.6, 4.0, 10.25 | 1.0616, 0.5 | 0 of 690 |
| 1001685937451957930 | 5 / 5 | 46 | 1557.0 (1372-1880) | 1517.0 (1332-1840) | 4.7 / 0.5 / 0.3 / 6.7 / 3.2 / 2.6 / 82.0 | 0, 46 (0/46/0) | 20 | 0.988 | 0.524, 0.0048, 0.021; 788/788 | 2.67, 3.25, 5.0 | 1.0907, 0.7628 | 0 of 407 |
| 6329801112132740259 | 5 / 5.136 | 46 | 1719.0 (1413-1933) | 1644.2 (1394-1814) | 3.9 / 0.5 / 0.4 / 2.9 / 2.9 / 4.3 / 85.0 | 15, 41 (5/26/15) | 20 | 1.0 | 0.477, 0.0010, 0.032; 682/682 | 4.0, 8.0, 11.0 | 1.0011, 0.7523 | 0 of 809 |
| 6036017085589997374 | 5 / 5 | 46 | 1663.1 (1491-1833) | 1502.2 (1336-1673) | 3.9 / 0.5 / 0.3 / 3.5 / 3.0 / 9.7 / 79.1 | 46, 47 (0/3/41) | 20 | 1.035 | 0.412, 0.0019, 0.021; 716/716 | 3.0, 4.67, 7.0 | 1.024, 0.7624 | 0 of 596 |
| 4383110225690451409 | 5 / 5 | 46 | 1615.0 (1425-1846) | 1475.9 (1265-1686) | 4.2 / 0.5 / 0.3 / 8.1 / 3.1 / 8.6 / 75.3 | 41, 37 (1/12/33) | 22 | 1.003 | 0.432, 0.0010, 0.065; 743/743 | 2.5, 3.0, 10.0 | 1.0897, 0.9424 | 0 of 627 |
| 10975598825187167680 | 5 / 5 | 46 | 1722.4 (1526-1884) | 1557.2 (1366-1724) | 4.0 / 0.4 / 0.2 / 6.5 / 2.9 / 9.6 / 76.4 | 48, 46 (0/0/44) | 23 | 1.036 | 0.522, 0.0029, 0.021; 722/722 | 2.0, 3.0, 6.33 | 0.9931, 0.1894 | 0 of 459 |
| 5467273169533009601 | 5 / 5.03 | 46 | 1524.5 (1307-1643) | 1367.2 (1147-1483) | 4.3 / 0.6 / 0.2 / 6.7 / 3.2 / 10.3 / 74.7 | 45, 46 (0/3/41) | 22 | 1.004 | 0.426, 0.0019, 0.043; 724/724 | 1.0, 3.4, 12.0 | 1.0826, 0.2821 | 0 of 591 |
| 4734920099712499969 | 5 / 5 | 46 | 1741.6 (1447-1905) | 1601.6 (1327-1745) | 4.5 / 0.4 / 0.2 / 10.1 / 2.8 / 8.0 / 73.9 | 46, 23 (0/23/23) | 21 | 1.053 | 0.439, 0.0038, 0.032; 778/778 | 2.33, 4.8, 12.0 | 1.0111, 0.7247 | 0 of 713 |
| 1570278518166000558 | 5 / 6.18 | 46 | 1891.2 (1495-2049) | 1759.9 (1335-1889) | 3.9 / 0.5 / 0.1 / 9.1 / 2.6 / 6.9 / 77.0 | 35, 46 (0/11/35) | 22 | 1.029 | 0.454, 0.0010, 0.054; 840/840 | 2.5, 4.0, 11.67 | 1.0206, 0.7384 | 0 of 652 |
| 3830344756309592292 | 5 / 5 | 46 | 1669.7 (1403-1897) | 1567.9 (1283-1777) | 4.0 / 0.5 / 0.3 / 6.6 / 3.0 / 6.1 / 79.6 | 39, 0 (7/39/0) | 23 | 0.979 | 0.435, 0.0010, 0.043; 765/765 | 3.0, 3.33, 11.5 | 1.0575, 0.8838 | 0 of 603 |
| 9254619518473077140 | 5 / 7.657 | 46 | 1647.0 (1465-1762) | 1489.6 (1305-1618) | 4.2 / 0.6 / 0.2 / 5.7 / 3.0 / 9.6 / 76.7 | 45, 46 (0/1/45) | 24 | 1.001 | 0.525, 0.0010, 0.043; 777/777 | 3.0, 5.8, 10.5 | 1.0333, 0.7984 | 0 of 724 |
| 4701969716356831982 | 5 / 5.087 | 46 | 1771.7 (1658-1880) | 1635.2 (1498-1741) | 4.1 / 0.8 / 0.3 / 6.6 / 2.8 / 7.7 / 77.7 | 37, 46 (0/9/37) | 20 | 0.882 | 0.428, 0.0000, 0.021; 683/683 | 1.0, 2.5, 10.67 | 1.0089, 0.4756 | 0 of 615 |
| 14843108908847119356 | 5 / 5 | 46 | 1665.1 (1378-1948) | 1626.0 (1349-1948) | 3.9 / 0.4 / 0.2 / 5.5 / 3.0 / 2.4 / 84.6 | 0, 45 (2/43/1) | 21 | 1.026 | 0.550, 0.0029, 0.054; 720/720 | 3.0, 7.0, 14.33 | 1.0399, 0.8765 | 0 of 845 |
| 1758546030705738202 | 5 / 5 | 46 | 1861.0 (1302-2025) | 1741.0 (1182-1905) | 3.7 / 0.5 / 0.3 / 6.5 / 2.7 / 6.5 / 80.0 | 46, 0 (0/46/0) | 21 | 0.964 | 0.563, 0.0000, 0.032; 739/739 | 0.0, 3.0, 7.0 | 0.984, 0.9036 | 0 of 336 |
| 1995372076782357328 | 5 / 5 | 46 | 1854.0 (1610-2386) | 1718.3 (1457-2226) | 5.7 / 0.5 / 0.3 / 8.8 / 2.7 / 7.3 / 74.7 | 38, 42 (1/10/35) | 20 | 1.296 | 0.483, 0.0048, 0.043; 684/684 | 3.67, 1.0, 10.0 | 0.9923, 0.9867 | 0 of 1111 |
| 13514522956726874792 | 5 / 9.488 | 46 | 1530.7 (1411-1647) | 1519.4 (1411-1628) | 5.3 / 0.5 / 0.6 / 8.2 / 3.2 / 0.7 / 81.5 | 0, 13 (33/13/0) | 21 | 1.062 | 0.444, 0.0019, 0.054; 866/866 | 1.0, 4.0, 9.17 | 1.0239, 0.2815 | 0 of 562 |
| 12621855678475765395 | 5 / 7.443 | 46 | 1801.8 (1626-1911) | 1663.6 (1511-1751) | 4.4 / 0.7 / 0.2 / 8.9 / 2.7 / 7.7 / 75.5 | 38, 45 (0/9/37) | 25 | 1.05 | 0.468, 0.0019, 0.032; 799/799 | 2.0, 4.0, 8.0 | 1.0131, 0.5687 | 0 of 571 |
| 7689157899245528466 | 5 / 5 | 46 | 1424.7 (1220-1621) | 1385.6 (1180-1581) | 4.4 / 0.8 / 0.3 / 4.0 / 3.5 / 2.8 / 84.2 | 0, 45 (1/45/0) | 23 | 1.051 | 0.546, 0.0019, 0.010; 658/658 | 3.0, 5.0, 14.5 | 1.0693, 0.534 | 0 of 834 |
| 13748185708697958051 | 5 / 5 | 46 | 1489.7 (1339-1749) | 1429.7 (1219-1670) | 4.4 / 0.6 / 0.2 / 2.6 / 3.3 / 4.0 / 84.9 | 23, 0 (23/23/0) | 23 | 0.99 | 0.550, 0.0010, 0.021; 760/760 | 2.0, 4.0, 10.0 | 1.0798, 0.6227 | 0 of 632 |
| 15390672246448918260 | 5 / 6.455 | 46 | 1843.1 (1627-1963) | 1715.2 (1507-1803) | 4.6 / 0.4 / 0.5 / 7.6 / 2.7 / 6.9 / 77.2 | 35, 42 (3/9/34) | 24 | 1.091 | 0.416, 0.0000, 0.021; 743/743 | 4.0, 4.33, 14.0 | 0.9732, 0.3059 | 0 of 833 |
| 11714203776604047358 | 5 / 5.052 | 46 | 1634.3 (1357-1744) | 1478.7 (1317-1584) | 4.6 / 0.7 / 0.3 / 5.8 / 3.0 / 9.5 / 76.0 | 45, 44 (0/3/43) | 25 | 1.043 | 0.473, 0.0029, 0.021; 767/767 | 3.0, 6.17, 12.33 | 0.991, 0.5056 | 0 of 935 |

Realized strength grid (`final4_grid.tsv`):

| member transfer \ together | 0.3 | 0.7 | 1 | level, pooled |
|---|---|---|---|---|
| 0.01 | 1.9 (0-7.17); together 0.20; powered at 22/22 | 2.3 (0-10); together 0.54; powered at 21/22 | 3.0 (0-10); together 1.00; powered at 22/22 | 2.3 |
| 0.03 | 3.5 (0-8.5); together 0.24; powered at 22/22 | 4.3 (0-10.75); together 0.54; powered at 22/22 | 5.0 (0-17); together 1.00; powered at 22/22 | 4.0 |
| 0.4 | 9.8 (1.5-111.17); together 0.30; powered at 22/22 | 10.3 (3.25-18.67); together 0.70; powered at 22/22 | 14.0 (2-24); together 1.00; powered at 22/22 | 10.3 |

## 13. Where the card stands after the refit

Committed tree, default seed unless named (*measured*): every synthworld test outside the null module and `tests/r0`:
"4 failed, 500 passed, 6 skipped in 1881.21s" (`tests_final_refit.txt`); the null module "2 failed, 12 passed in
1353.09s" at the default seed and "1 failed, 12 passed, 1 skipped in 1292.63s" at 20261113
(`nulls-final-*.txt`).

| assertion | bound | value now | cause | owner |
|---|---|---|---|---|
| `test_plasmids_are_whole_replicons_of_their_preset_size` | a genome without a plasmid | 0 genomes without (pls0 in 45, pls1 in 5) | §4 item 1 | PI (b) |
| `test_n0_keeps_the_background_process_in_distribution` | ratio [0.9, 1.1], r >= 0.9 | 1.023, 0.283 (20261113: 1.062, 0.500) | §4 item 2 | PI (c) |
| `test_orthologous_igrs_are_near_identical_in_close_relatives` | >= 200 cross-phylum pairs | 150 (identity 0.931 passes) | needs background transfer >= about 0.045 (above) | this card: conflicts with the motif result at 0.05 |
| `test_dec15_a_fates_family_splits_a_unit_and_the_truth_says_so` | a unit broken under fates "family" | none broken | at duplication 0.008 and member duplication 0, no tandem copy lands in a unit in that world | this card |
| `test_decoy_transfers_land_only_where_the_partner_is` | >= 3 transfers per decoy | decoy 1297: 2 (others 6-54) | its partner module is rare at this draw | this card |
| `test_n1_label_pairs_are_no_more_adjacent_than_matched_random_family_pairs` | relabeled pairs not more adjacent than random, p | observed 0.585 against 0.113, p 0.0055 (passes at 20261113) | 41 of 100 members relabeled (60 at factor 0.2 with the 0.06 level, where it passed); account not tested | this card |

Resolved by the refit: the hub test's core (253 root families in >= 40 genomes at the default seed, 228-513 over the
22 seeds), N1's comparator exists (0 unmatched), the burst identity (0.931), the planted motif, the per-seed test (both
plasmid types at 5 of 6, levels 2.55/4.08/9.83).

**For the PI, (b)** (`z_pls_*`, `final4`; 22 seeds): keeping §5.1's range (at most 2 per genome, some genomes
without): plasmids 0.15/0.3 give both types at 15 of 22 seeds with some genome without at 22 of 22 and more than 2 at 1;
0.1/0.3 give 17 of 22, 22 of 22 and 2. Survival at every seed: 3.0/0.1 gives both types at 22 of 22 with every genome a
carrier at every seed (no genome without at any seed) and more than 2 plasmids in some genome at 7. The committed 0.8/0.4:
15 of 22, some genome without at 10, more than 2 at 4. **(c)** (`final4`, `--n0`): N0's background size ratio is in
[0.9, 1.1] at 22 of 22 seeds (0.973-1.091; default 1.023); r >= 0.9 at 3 of 22 (0.19-0.99; default 0.283); a two-sample
Kolmogorov-Smirnov test of the 46 background sizes (asymptotic p) gives D 0.17-0.72, p >= 0.01 at 11 of 22 seeds
(default D 0.283, p 0.041; 20261113 D 0.652, p < 1e-4): N0's sizes are shifted from the signal's at about half the
seeds, so a KS replacement would fail there too; the ratio alone holds everywhere.

## 14. Last pass: the four failures of this card (measurement; the preset is unchanged)

All *measured*, wilkcomp02, 2026-11-15 00:00-03:30 CST, at most 8 threads. The preset on the branch keeps the values of
§12 (duplication 0.008, transfer 0.02, loss 0.048, origination 0.065, member transfer 0.01/0.03/0.4); no candidate below
passed every assertion other than (b) and (c), so none is committed.

**Backbone scans** (22 seeds; transfer 0.04 / 0.045 with duplication 0.012 / 0.015, loss derived so that
2 loss - 2 dup - 3 transfer = -0.01, origination 0.035, hub multiplier 2.0 / transfer; `d_committed` is the branch
preset; new column: `decoy_transfers_min`, the fewest transfer rows of any decoy):

| scan | chromosome mean (in band) | default mean | planted share % | copies per present member | exclusion breaches | both plasmid types | a genome > 2 plasmids | genomes without a plasmid at default | all cells powered | levels pooled (increasing) | core >= 40 genomes: default, range (>= 100) | N0 ratio, r at default | N1 unmatched planted at default | burst multicopy junction keys at default |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| d_committed | 1367-1760 (22/22) | 1448 | 4.4-6.5 | 0.88-1.30 | 0 | 15/22 | 4 | 0 | 21/22 | 2.3, 4.0, 10.3 (21/22) | 253, 228-513 (22) | - | 0 | 203 of 5301 |
| d_t0.04_d0.012_l0.067 | 1401-1764 (22/22) | 1569 | 4.2-5.6 | 0.92-1.19 | 0 | 20/22 | 5 | 2 | 20/22 | 3.0, 5.2, 10.7 (20/22) | 210, 201-437 (22) | - | 0 | 291 of 5500 |
| d_t0.04_d0.015_l0.0685 | 1398-1724 (22/22) | 1523 | 4.2-5.7 | 0.87-1.14 | 0 | 13/22 | 5 | 0 | 22/22 | 3.5, 5.2, 10.3 (19/22) | 195, 185-425 (22) | - | 0 | 283 of 5486 |
| d_t0.045_d0.012_l0.076 | 1250-1711 (22/22) | 1388 | 4.3-6.6 | 0.92-1.14 | 0 | 15/22 | 7 | 3 | 21/22 | 3.0, 5.5, 11.7 (18/22) | 168, 143-423 (22) | - | 0 | 204 of 5121 |
| d_t0.045_d0.015_l0.0775 | 1320-1786 (22/22) | 1460 | 3.8-6.9 | 0.87-1.32 | 0 | 15/22 | 5 | 0 | 20/22 | 3.2, 5.4, 11.2 (19/22) | 156, 156-390 (22) | - | 0 | 263 of 5280 |

Decoys with at least 3 transfers each (the backbone test's bound), seeds out of 22: committed 7, candidates 8, 11, 10, 8
(default seed minimum 2, 11, 37, 8, 9). **The decoy assertion fails at about two thirds of the seeds at every setting:
it is a property of the preset (decoy transfer 0.4 against loss 0.4 is critical, so a decoy's copy line often dies or
stays small, and the partner modules are the patchiest ones by design), not a draw at the default seed.**

**Full worlds at the default seed** (scratch copies of the branch with the candidate values; `tests_candidate_*.txt`,
`cnulls-*.txt`):

| assertion | committed (§13) | A: transfer 0.045, dup 0.015, loss 0.0775 | B: transfer 0.045, dup 0.012, loss 0.076 |
|---|---|---|---|
| cross-phylum IGR pairs >= 200 | 150 | **200** | 195 |
| DEC-15 (a): a unit broken under fates "family" | none | broken units present (passes) | passes |
| decoys >= 3 transfers | 2 | 9-49 | 8-90 |
| N1 label pairs no more adjacent than matched random | p 0.0055 | passes | passes |
| hub core >= 100 | 253 | 156 | 168 |
| N1 comparator, burst identity | yes, 0.931 | yes, 0.952 | yes, 0.945 |
| planted motif found | yes | **no**: twin (below) | yes |
| junction lengths, per-seed test | pass | pass (levels 3.9/5.9/13.8) | pass (3.5/6.0/12.7) |
| plasmids (b) | fails | fails | **passes** (3 genomes without) |
| N0 r (c) | 0.283 | 0.587 | 0.102 |
| totals | 4 + 2 nulls | backbone/intergenic/per-seed "2 failed, 42 passed"; nulls "2 failed, 12 passed" (digest pin, N0) | "1 failed, 43 passed"; nulls "2 failed, 12 passed" (digest pin, N0) |

**The motif failure is a real second planted instance, not a coincidence.** In the transfer-0.05 world the winning window
at 1,626,402 is `n49_g114798.rna`, an instance of the same emission (L05.8, family F000199, genome g002) on the gene copy
`n49_g114798`; the failing instance is on `n49_g114799`. Both copies come from one tandem duplication in lineage n49
(`events.tsv` event 6336: `duplication 199 n49 ... n49_g98824 -> n49_g114798;n49_g114799`), so each carries the clade's
emission, and the twin, with fewer substitutions, wins the scan of the shared intergenic region. Candidate A's failure is
the same (`n89_g193995` and `n89_g193996`, emission L10.3, 8 against 7 substitutions). **The test's premise (one planted
instance per region, so the instance at its coordinates is the best window) is broken by tandem duplication of a
motif-carrying gene, a thing the world does.** Same-emission twins within 1 kb (a scratch count over `replicons/rna.tsv`, not in the tree):
session 10's world 4 (1 closer to the emission), committed 0, transfer 0.035 28 (10 closer), A 16 (4 closer), B 6 (1
closer); the test failed in 4 of the 9 full worlds built this session. It will fail at transfer 0.04-0.045 with some
probability at any seed whenever duplication is not zero.

**N1** relabels 11 of 100 members at the preset before §12, 60 at factor 0.2 with the 0.06 level, 41 at the committed
preset, 39 (A) and 35 (B); session 9 had 27 or 28. Lower turnover gives more background families within 4 genomes of a
member's profile and 25% of its transfers. The adjacency excess at the committed preset (a scratch script, not in the tree, on the kept null
worlds): of the 41 replacements, 1 is adjacent to a member of its own module in any genome; of the 82 evaluated pairs,
the 32 (moved, unmoved) pairs are adjacent in 2 genome instances, the 50 (moved, moved) pairs in 47 (6 pairs adjacent
somewhere): the excess comes from pairs of replacement families that are adjacent to each other, two neighbouring
background families whose shared history gives both a module-like profile. The test passes at A and B and at
20261113.

**What cannot hold together under one origination axis** (measured trade-off): the cross-phylum IGR count needs the
background transfer at about 0.045 or more (118, 169, 188, 195-200, 236-262 at 0.02, 0.03, 0.035, 0.045, 0.05); the
other assertions of this card hold there (A: all but the motif and plasmids; B: all but the IGR count by 5). The motif
test is not a rate trade-off: it fails whenever a tandem twin of a planted instance is closer to the emission, in 4 of
9 worlds this session, at duplication from 0.008 to 0.04. The decoy test holds at a third of the seeds at every
setting. So the decisions left are about the tests' premises, not the rates: (i) whether the motif test should allow a
same-emission twin, (ii) whether the decoy test's 3 transfers is a per-seed requirement for decoys at critical rates,
(iii) (b) plasmids and (iv) (c) N0's r, with the numbers of §13.


## 15. Session 12: the preset of DEC-44, re-pinned

Card F-02, session 12 (branch `task/F-02-12`: session 11's `task/F-02-11` at `5de6adc` with `plan/task-cards` merged,
which brought DEC-43 to DEC-46; the merge was clean). DEC-44 is the specification. Host and commands for every
*measured* number in this section: wilkcomp02, 2026-11-15 13:25-18:00 CST (Sunday; the account at 18 to 26 of its 56
threads with other agents, this session's at most 8), image `bio-box-test_faec1f4.sif`, units `f0212-*` at nice 19,
one thread per process. Scans: `python scripts/stress/f0210_seed_scan.py scan /tmp/if2453/f0212/<dir> <label> --jobs
<5|6|7> --seeds default,20261113,n2 [--n0] [--set a.b=JSON ...]` (overrides per label in `session12/<label>.json`);
full worlds at the default seed: `session12/scripts/igr.py` (builds S1-S5 with the overrides and calls the two
sequence-level test functions of `tests/r2/test_synthworld_intergenic.py`, the migrated motif test included); the
scripts are in `session12/scripts/`, the scan tables and logs in `session12/`. Row summaries: `session12/scripts/summ.py`
(session 11's `vsumm.py` columns plus the DEC-44 ones).

**Outcome: the preset of DEC-44 is in and re-pinned, and two existing assertions outside the migrations fail under it
(stop and report, §15.6).** The preset is session 11's candidate B with plasmids 0.1/0.3 and decoy transfer 8.0;
background transfer did not have to move: with DEC-44's plasmid and decoy rates the default world at 0.045 has 206
equal-length cross-phylum IGR pairs. Clauses (1), (3), (4) and (5) are met, the three migrations are in (§15.2), and
every stage digest changed once and is re-pinned (§15.3). Clause (2) stops: the decoy rate that holds "3 transfers per
decoy" at every seed (8.0) makes the null module's root-presence KS fail at both of its seeds (it holds at 22 of 22
seeds with decoy 0.4), and 16.0, which holds both there, leaves the default world 7 cross-phylum IGR pairs short of
clause (5). Separately, the default world under the final preset has two contigs of one genome with identical
sequence, so two truth genes share one `feature_uid` (`test_every_truth_gene_reads_its_s8_gene_on_the_emitted_contig`).
Neither assertion was changed.

### 15.1 The fits (22 seeds each, one row per scan)

Labels `e_t<T>_l<L>_dec<D>`: candidate B (duplication 0.012, origination 0.035) with background transfer T, loss L =
(-0.01 + 2 x 0.012 + 3 T) / 2 (session 11's rule, *derived*), hub multiplier 2.0 / T, plasmids 0.1/0.3 and decoy
transfer D. `final12` is the preset as committed, no overrides, with `--n0`.

| scan | chromosome mean (in band) | default mean | planted share % | copies per present member | exclusion breaches | both plasmid types (default) | some genome without a plasmid | a genome > 2 plasmids | all cells powered | levels pooled (increasing) | core >= 40 genomes: default, range (>= 100) | decoys >= 3 transfers: seeds, fewest, default | N0 ratio range (in [0.9, 1.1]) | N1 unmatched planted at default |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| e_t0.045_l0.076_dec0.4 | 1317-1854 (22/22) | 1460 | 4.4-6.4 | 0.87-1.11 | 0 | 14/22 (yes) | 21/22 | 2 | 21/22 | 4.0, 5.8, 11.5 (18/22) | 136, 136-393 (22) | 9/22, 0, 3 | - | 0 |
| e_t0.045_l0.076_dec0.8 | 1328-1754 (22/22) | 1536 | 3.9-7.6 | 0.87-1.15 | 0 | 16/22 (yes) | 22/22 | 1 | 22/22 | 3.3, 6.0, 12.0 (19/22) | 177, 156-392 (22) | 15/22, 0, 106 | - | 0 |
| e_t0.045_l0.076_dec1.2 | 1354-1822 (22/22) | 1492 | 4.7-8.2 | 0.84-1.09 | 0 | 16/22 (yes) | 22/22 | 2 | 20/22 | 3.4, 5.0, 11.0 (18/22) | 156, 156-391 (22) | 15/22, 0, 125 | - | 0 |
| e_t0.045_l0.076_dec4.0 | 1448-1808 (22/22) | 1654 | 5.6-8.4 | 0.91-1.11 | 0 | 18/22 (yes) | 22/22 | 3 | 22/22 | 3.7, 6.3, 12.0 (19/22) | 195, 182-413 (22) | 21/22, 1, 263 | - | 0 |
| **e_t0.045_l0.076_dec8.0** | 1403-1815 (22/22) | 1570 | 5.4-8.3 | 0.75-1.11 | 0 | 14/22 (yes) | 21/22 | 0 | 22/22 | 4.0, 6.0, 12.4 (19/22) | 181, 175-406 (22) | **22/22, 14, 220** | - | 0 |
| e_t0.045_l0.076_dec16.0 | 1420-1739 (22/22) | 1554 | 5.9-8.7 | 0.93-1.22 | 0 | 11/22 (yes) | 22/22 | 0 | 21/22 | 4.5, 7.0, 11.0 (19/22) | 173, 173-374 (22) | 22/22, 22, 314 | - | 0 |
| e_t0.05_l0.082_dec0.8 | 1381-1740 (22/22) | 1509 | 4.9-7.0 | 0.84-1.13 | 0 | 16/22 (yes) | 22/22 | 3 | 21/22 | 3.6, 5.8, 12.5 (18/22) | 158, 143-412 (22) | 12/22, 0, 149 | - | 0 |
| e_t0.05_l0.082_dec1.2 | 1301-1788 (22/22) | 1559 | 5.1-8.0 | 0.87-1.10 | 0 | 13/22 (yes) | 21/22 | 0 | 21/22 | 3.5, 6.2, 13.3 (19/22) | 143, 143-376 (22) | 12/22, 0, 0 | - | 0 |
| e_t0.05_l0.082_dec8.0 | 1450-1770 (22/22) | 1600 | 6.3-8.2 | 0.89-1.12 | 0 | 14/22 (yes) | 22/22 | 2 | 21/22 | 4.5, 6.5, 12.4 (16/22) | 192, 167-365 (22) | 21/22, 0, 269 | - | 0 |
| e_t0.055_l0.0895_dec8.0 | 1376-1813 (22/22) | 1613 | 4.8-8.4 | 0.84-1.24 | 0 | 13/22 (**no**) | 22/22 | 0 | 20/22 | 4.3, 6.9, 13.0 (18/22) | 177, 129-418 (22) | 21/22, 1, 227 | - | 0 |
| **final12 (committed)** | 1403-1815 (22/22) | 1570 | 5.4-8.3 | 0.75-1.11 | 0 | 14/22 (yes) | 21/22 | 0 | 22/22 | 4.0, 6.0, 12.4 (19/22) | 181, 175-406 (22) | 22/22, 14, 220 | 0.893-1.086 (21/22) | 0 |

`final12` reproduces `e_t0.045_l0.076_dec8.0` row for row (the overrides and the committed YAML give the same worlds).

**Decoys (DEC-44 (2)).** Why a decoy has fewer than 3 transfers (`decoys.py`, a replay of ZOMBI2's
`genome_events.tsv` on worlds kept at the seven failing seeds of `dec1.2`, scan `z_dec*` in `session12/`, output in
`decoys_kept_worlds.txt`): either both crown copies of the decoy die (loss rows, or the lineage dies) before either
transfers (4 of 7: 0 transfers, 2 or 3 losses, while the partner module is held by many lineages), or the partner
module is held by at most one lineage from t 0.56 to 1.25 on the tree's clock, so there is no recipient (3 of 7; the
partners are modules 1, 9, 11, 19 and 22, all at the two lowest member transfer levels, chosen as the patchiest in
session 2). Both are races against an early clock, so only a much higher rate wins them at every seed: 0.4, 0.8, 1.2,
4.0, 8.0, 12.0 and 16.0 give 9, 15, 15, 21, 22, 20 and 22 of 22 seeds (12.0 in §15.6 (A)). 8.0 is the smallest
tried that holds at every seed (the fewest at any seed 14) with every backbone property of this table and the other
DEC-44 clauses; the test's other two assertions (every decoy transfer lands where the partner is; mean partner share
among the other lineages <= 0.9) hold in the test itself at the default seed (§15.5); the null module's root-presence
KS does not (§15.6 (A)). Consequence (*measured*, the default world below): a decoy is now a multi-copy family in its
partner's carriers, 40.8 decoy genes per genome (2.2 before; 2.5% of the genes), still inside the planted-share bound
(planted + profile-only + decoy 5.4-8.3% over the 22 seeds against 15.2%).

**Background transfer (DEC-44 (5)).** Full default worlds with plasmids 0.1/0.3 and decoy 8.0 (`igr_t*.log`):

| background transfer, loss, hub multiplier | equal-length cross-phylum IGR pairs (>= 200) | burst IGR identity (>= 0.9) | across phyla | planted motif test (migrated) | best windows another planted instance of the emission |
|---|---|---|---|---|---|
| **0.045, 0.076, 44.444444** | **206** | 0.946 | 0.641 | passes (1,150 of 1,159 unique best) | 4 |
| 0.05, 0.082, 40 | 177 (fails) | 0.943 | 0.682 | passes (1,141 of 1,153) | 3 |
| 0.055, 0.0895, 36.363636 | 221 | 0.944 | 0.690 | passes (1,079 of 1,095) | 8 |

The count is not monotone in T at one seed (session 11: 195 at 0.045 with its plasmid and decoy rates, 236-262 at
0.05): it is one draw of the world. 0.045 is the smallest transfer that holds at the default seed, and at 0.05 and
0.055 a decoy has fewer than 3 transfers at one seed and, at 0.055, the default seed loses a plasmid type (table
above). So background transfer stays at 0.045, and loss (0.076) and origination (0.035) stay: the chromosome mean is
inside 1,500 +- 30% at every seed (1,403-1,815) and the core is >= 100 root families in >= 40 genomes at every seed
(175-406). The member transfer levels stay (0.01/0.03/0.4): pooled medians 4.0, 6.0, 12.4 over 22 seeds, 4.1, 6.2,
16.0 over the per-seed test's six seeds (`six.py`: steps 1.52x and 2.58x against its 1.5x), so they are still distinct.

**Plasmids (DEC-44 (3)).** At 0.1/0.3 under the final preset both types live at the default seed (pls0 in 35 genomes,
pls1 in 2; 9 genomes without a plasmid, none with more than 2), and at 14 of 22 seeds (session 11 measured 17 of 22 at
its preset); some genome is without a plasmid at 21 of 22 seeds, and no genome carries more than 2 at any seed.

**N0 (DEC-44 (4)).** The background size ratio is in [0.9, 1.1] at 21 of 22 seeds (0.893-1.086): at
10975598825187167680 it is 0.893. The null module runs at the default seed (0.946) and 20261113 (1.026), where it
holds; the one seed outside is reported here, not asserted (the module does not run there).

### 15.2 The three migrations (DEC-44 (1), (3), (4)); no other assertion changed

- `tests/r2/test_synthworld_intergenic.py::test_every_planted_motif_is_recovered_at_its_truth_coordinates_by_an_independent_scan`
  (DEC-44 (1)). Before: for every planted instance, `(start, strand) in winners or divergence > 0.3`, where `winners`
  are the windows of the instance's intergenic region (both strands) at the lowest Hamming distance from the clade's
  emission. After: `(start, strand) in winners or twin or divergence > 0.3`, where `twin` holds only when the instance
  is not among the winners and *every* winner is another planted instance of the same emission on the same replicon,
  at that instance's own start and strand (read from `rna.tsv` by a loop in the test, `planted_at`). Unchanged: the
  per-instance checks before the scan, the scan itself, the divergence escape, and the count of instances whose own
  window is the unique best (still >= 0.95 of them). A best window that is no planted instance of the emission fails as
  before. The test now prints how many instances have only other planted instances of their emission as best windows
  (4 at the default seed under the final preset; 3 and 8 in the transfer-0.05 and 0.055 worlds of §15.1).
- `tests/r2/test_synthworld_preset_seeds.py`: `test_both_plasmid_types_live_at_the_share_the_refit_reaches` became
  `test_both_plasmid_types_live_at_the_default_seed_and_are_reported_per_seed` (DEC-44 (3)). Before: both types at no
  fewer than half of the six seeds (`PLASMID_BOTH_SHARE` 0.5, removed). After: both types at the default seed; the
  carriers of each type and the most plasmids in one genome are printed per seed. §5.1's per-genome range is still
  asserted at the default seed by `test_plasmids_are_whole_replicons_of_their_preset_size` (unchanged).
- `tests/r2/test_synthworld_nulls.py::test_n0_keeps_the_background_process_in_distribution` (DEC-44 (4)). Before:
  `0.9 <= ratio <= 1.1 and r >= 0.9`. After: `0.9 <= ratio <= 1.1`; r is still computed, printed and recorded as a
  property. Every other assertion of the null module is unchanged (KS on root presence, + strand share, tight share,
  GC, and every other test).

Re-pinned, not migrated: the digests (§15.3). The spec: `synthworld_spec.md` §2.1's plasmid paragraph still says
"transfer 0.3 and loss 0.7 per plasmid ... (17, 20 and 9 of 46)" and §5.1's backbone paragraph the session-2 rates;
DEC-44 and the preset supersede both (not edited here, as the brief says).

### 15.3 Re-pinned digests and the default world before and after

Every stage digest of the default world (seed 20261112) changed once (*measured*, wilkcomp02, 2026-11-15, `python
scripts/stress/f0210_world_summary.py build /tmp/if2453/f0212/world-final 20261112` on the committed tree, then
`summary`; `session12/default_world_summary_after12.json`). Pinned in `tests/r2/test_synthworld_emit.py` (DIGESTS,
EMIT_DIGEST, WORLD_DIGEST), `test_synthworld_nulls.py` (SIGNAL, SIGNAL_WORLD), `test_synthworld_assembly.py` (DIGESTS,
ASSEMBLY_DIGEST) and `test_synthworld_fusions.py` (S6_DIGEST); the pinned tests pass (§15.5):

| stage | pre-redesign (session 10) | session 11 (`89cfb46`) | session 12 |
|---|---|---|---|
| backbone | 025fdd0d27fdc65e678782216ecd21fb | 545240f325b552247e0665459ccb1c0e | f65094aa0fb1c9bed885bc6bc2f4e825 |
| sequences | 34ebce1324413a261acd874d7ee4fe3e | 294beafcfe7efa64b582276549e353b9 | e9bf8833a898541f1bc9606c97d55463 |
| replicons | 55bda97155d75148e5cdca7b90a1d7be | 98b0a2fada22354766bdf4d12df6ade6 | 1930863ca47dcc67a5fde9deb1d12f84 |
| mobile | 0bd7f0db5691de58ea5e1a0df1e6e570 | bbe7289e32176b4e9874c706be6ccf03 | 1596365eedc377773ba85ef5aa7ec593 |
| fusions | db09abf3fe79ae34a73e3ab0664a9b13 | 818031fea3ed35abab093554605474ad | d2f80ac770d8cd78346467aed2efc9c7 |
| assembly | 2f425158dc38f2302c1d02e9ebed30da | 74ea9afac18fdab3fe8470728679f544 | 03b3e1b97dad03a8b2b7d78d755eebaf |
| sources | 770007bde8615f41d9dd5f2eb912e753 | a333a1e97bb6c0ba75b0c870f3ad28d8 | 5366fca6f2796b3956265d248b6529d9 |
| emit | a2f0e448484dab63d41281a5ef6a00c0 | 0f1ae8e62576ee05dc1c2fa9fe78f378 | 267f3a0200ddc3f7639e3c2fd6ffc134 |
| world | f0262027e4ccdcbc72cd7d90f23fbb5c | 765b20a0d7d0051d62c5036143d6ced5 | 684e539e11bd1ecbc5d3df6fe6405e7f |

The default world (S1-S11 at 20261112; `f0210_world_summary.py summary`; pre-redesign:
`2026-11-14_synthworld_seeds/default_world_summary.json`, session 11: `default_world_summary_refit.json`, session 12:
`session12/default_world_summary_after12.json`; all *measured*):

| quantity | pre-redesign (session 10) | session 11 (89cfb46) | session 12 |
|---|---|---|---|
| arrays | 10 | 10 | 10 |
| chromosome_genes.max | 3798 | 1606 | 1801 |
| chromosome_genes.mean | 1872.7 | 1447.5 | 1570.5 |
| chromosome_genes.min | 969 | 1243 | 1321 |
| families_with_extant_copies_by_kind.cargo | 16 | 16 | 16 |
| families_with_extant_copies_by_kind.decoy | 5 | 5 | 5 |
| families_with_extant_copies_by_kind.hub | 5 | 5 | 5 |
| families_with_extant_copies_by_kind.originated | 1399 | 2908 | 1415 |
| families_with_extant_copies_by_kind.planted | 100 | 100 | 100 |
| families_with_extant_copies_by_kind.plasmid | 160 | 160 | 160 |
| families_with_extant_copies_by_kind.profile_only | 9 | 9 | 9 |
| families_with_extant_copies_by_kind.root | 998 | 999 | 999 |
| families_with_extant_copies_by_kind.vehicle_chassis | 6 | 6 | 6 |
| genes_per_genome.max | 3918 | 1741 | 1921 |
| genes_per_genome.mean | 1937.0 | 1569.2 | 1663.6 |
| genes_per_genome.min | 969 | 1363 | 1432 |
| genes_per_genome_by_kind.background | 1318.0 | 1209.9 | 1254.8 |
| genes_per_genome_by_kind.decoy | 6.7 | 2.2 | 40.8 |
| genes_per_genome_by_kind.hub | 49.3 | 49.4 | 49.3 |
| genes_per_genome_by_kind.planted | 380.7 | 71.5 | 77.6 |
| genes_per_genome_by_kind.plasmid | 64.3 | 121.7 | 93.0 |
| genes_per_genome_by_kind.profile_only | 29.6 | 9.1 | 8.2 |
| genes_per_genome_by_kind.vehicle | 88.3 | 105.3 | 139.9 |
| genomes | 46 | 46 | 46 |
| kinds.assembly/callerr.tsv:kind.alt_code | 6551 | 5342 | 6300 |
| kinds.assembly/callerr.tsv:kind.alt_start_downstream | 30 | 30 | 30 |
| kinds.assembly/callerr.tsv:kind.alt_start_upstream | 30 | 30 | 30 |
| kinds.assembly/callerr.tsv:kind.frameshift_del | 30 | 30 | 30 |
| kinds.assembly/callerr.tsv:kind.frameshift_ins | 30 | 30 | 30 |
| kinds.assembly/callerr.tsv:kind.internal_stop | 40 | 40 | 40 |
| kinds.assembly/callerr.tsv:kind.smorf | 40 | 40 | 40 |
| kinds.assembly/edits.tsv:kind.alt_start_downstream | 30 | 30 | 30 |
| kinds.assembly/edits.tsv:kind.alt_start_upstream | 88 | 87 | 78 |
| kinds.assembly/edits.tsv:kind.frameshift_del | 30 | 30 | 30 |
| kinds.assembly/edits.tsv:kind.frameshift_ins | 30 | 30 | 30 |
| kinds.assembly/edits.tsv:kind.internal_stop | 40 | 40 | 40 |
| kinds.assembly/edits.tsv:kind.recode_codon | 26039 | 20986 | 24347 |
| kinds.assembly/edits.tsv:kind.recode_overlap | 334 | 264 | 294 |
| kinds.assembly/edits.tsv:kind.recode_stop | 2189 | 1816 | 2053 |
| kinds.assembly/edits.tsv:kind.smorf | 40 | 40 | 40 |
| kinds.assembly/genomes.tsv:assembly_class.complete | 18 | 18 | 18 |
| kinds.assembly/genomes.tsv:assembly_class.draft | 18 | 18 | 18 |
| kinds.assembly/genomes.tsv:assembly_class.mag | 10 | 10 | 10 |
| kinds.assembly/replicons.tsv:kind.chromosome | 46 | 46 | 46 |
| kinds.assembly/replicons.tsv:kind.plasmid | 38 | 50 | 37 |
| kinds.fusions/edits.tsv:kind.fusion | 52 | 52 | 52 |
| kinds.fusions/edits.tsv:kind.nogd | 15 | 7 | 7 |
| kinds.fusions/edits.tsv:kind.split | 14 | 14 | 14 |
| kinds.fusions/events.tsv:kind.fusion | 15 | 15 | 15 |
| kinds.fusions/events.tsv:kind.nogd | 2 | 2 | 2 |
| kinds.fusions/events.tsv:kind.split | 5 | 5 | 5 |
| kinds.fusions/fusions.tsv:kind.fusion | 52 | 52 | 52 |
| kinds.fusions/fusions.tsv:kind.split | 14 | 14 | 14 |
| kinds.fusions/replicons.tsv:kind.chromosome | 46 | 46 | 46 |
| kinds.fusions/replicons.tsv:kind.plasmid | 38 | 50 | 37 |
| kinds.mobile/insertions.tsv:kind.crispr | 30 | 30 | 30 |
| kinds.mobile/insertions.tsv:kind.prophage | 16 | 16 | 16 |
| kinds.mobile/replicons.tsv:kind.chromosome | 46 | 46 | 46 |
| kinds.mobile/replicons.tsv:kind.plasmid | 38 | 50 | 37 |
| kinds.sources/contam.tsv:kind.cross_clade | 10 | 10 | 10 |
| kinds.sources/dup.tsv:kind.exact | 108 | 105 | 110 |
| kinds.sources/dup.tsv:kind.near | 87 | 120 | 138 |
| kinds.sources/sources.tsv:kind.exact_dup | 2 | 2 | 2 |
| kinds.sources/sources.tsv:kind.near_dup | 2 | 2 | 2 |
| kinds.sources/sources.tsv:kind.original | 46 | 46 | 46 |
| phages.genes | [12 items] | [12 items] | [12 items] |
| phages.n | 12 | 12 | 12 |
| planted_share | 0.2153 | 0.0528 | 0.0761 |
| plasmid_genomes_by_type.pls0 | 18 | 45 | 35 |
| plasmid_genomes_by_type.pls1 | 20 | 5 | 2 |
| prophage_instances | 16 | 16 | 16 |
| seconds.S1-S11 | 238.61 | 210.1 | 215.9 |
| seconds.S1-S2 backbone | 37.1 | 31.44 | 51.97 |
| seconds.S10 sources | 6.56 | 5.89 | 7.56 |
| seconds.S11 emit | 24.52 | 21.97 | 21.87 |
| seconds.S3-S4 sequences | 86.73 | 73.79 | 66.72 |
| seconds.S5 replicons | 47.97 | 43.67 | 39.15 |
| seconds.S6 mobile | 11.47 | 13.96 | 9.48 |
| seconds.S7 fusions | 5.34 | 4.51 | 4.02 |
| seconds.S8-S9 assembly | 18.91 | 14.86 | 15.14 |
| stage_digests.assembly | 2f425158dc38f2302c1d02e9ebed30da | 74ea9afac18fdab3fe8470728679f544 | 03b3e1b97dad03a8b2b7d78d755eebaf |
| stage_digests.backbone | 025fdd0d27fdc65e678782216ecd21fb | 545240f325b552247e0665459ccb1c0e | f65094aa0fb1c9bed885bc6bc2f4e825 |
| stage_digests.emit | a2f0e448484dab63d41281a5ef6a00c0 | 0f1ae8e62576ee05dc1c2fa9fe78f378 | 267f3a0200ddc3f7639e3c2fd6ffc134 |
| stage_digests.fusions | db09abf3fe79ae34a73e3ab0664a9b13 | 818031fea3ed35abab093554605474ad | d2f80ac770d8cd78346467aed2efc9c7 |
| stage_digests.mobile | 0bd7f0db5691de58ea5e1a0df1e6e570 | bbe7289e32176b4e9874c706be6ccf03 | 1596365eedc377773ba85ef5aa7ec593 |
| stage_digests.replicons | 55bda97155d75148e5cdca7b90a1d7be | 98b0a2fada22354766bdf4d12df6ade6 | 1930863ca47dcc67a5fde9deb1d12f84 |
| stage_digests.sequences | 34ebce1324413a261acd874d7ee4fe3e | 294beafcfe7efa64b582276549e353b9 | e9bf8833a898541f1bc9606c97d55463 |
| stage_digests.sources | 770007bde8615f41d9dd5f2eb912e753 | a333a1e97bb6c0ba75b0c870f3ad28d8 | 5366fca6f2796b3956265d248b6529d9 |
| tsv_rows.assembly/callerr.tsv | 6751 | 5542 | 6500 |
| tsv_rows.assembly/contigs.tsv | 2431 | 2102 | 2239 |
| tsv_rows.assembly/edits.tsv | 28820 | 23323 | 26942 |
| tsv_rows.assembly/gene_contigs.tsv | 92488 | 75260 | 79762 |
| tsv_rows.assembly/genes.tsv | 90282 | 73364 | 77704 |
| tsv_rows.assembly/genomes.tsv | 46 | 46 | 46 |
| tsv_rows.assembly/repeats.tsv | 1344 | 1512 | 2005 |
| tsv_rows.assembly/replicons.tsv | 84 | 96 | 83 |
| tsv_rows.backbone/events.tsv | 50898 | 39792 | 57310 |
| tsv_rows.backbone/families.tsv | 2698 | 4208 | 2715 |
| tsv_rows.backbone/gene_tree_nodes.tsv | 262548 | 213418 | 242674 |
| tsv_rows.backbone/genes.tsv | 89102 | 72184 | 76524 |
| tsv_rows.backbone/genomes.tsv | 46 | 46 | 46 |
| tsv_rows.backbone/groups.tsv | 33 | 33 | 33 |
| tsv_rows.backbone/replicons.tsv | 84 | 96 | 83 |
| tsv_rows.backbone/taxa.tsv | 82 | 82 | 82 |
| tsv_rows.backbone/taxonomy.tsv | 46 | 46 | 46 |
| tsv_rows.backbone/truth_adjacency.tsv | 89102 | 72184 | 76524 |
| tsv_rows.backbone/truth_units.tsv | 3912 | 1497 | 1637 |
| tsv_rows.emit/ingest.tsv | 2327 | 2088 | 2243 |
| tsv_rows.emit/truth/truth_callerr.tsv | 8509 | 6933 | 8162 |
| tsv_rows.emit/truth/truth_contam.tsv | 10 | 10 | 10 |
| tsv_rows.emit/truth/truth_contig.tsv | 2648 | 2349 | 2509 |
| tsv_rows.emit/truth/truth_dup.tsv | 195 | 225 | 248 |
| tsv_rows.emit/truth/truth_edges.tsv | 2217 | 1982 | 2131 |
| tsv_rows.emit/truth/truth_family.tsv | 2810 | 4320 | 2827 |
| tsv_rows.emit/truth/truth_fusion.tsv | 69 | 69 | 69 |
| tsv_rows.emit/truth/truth_gene.tsv | 93883 | 76797 | 81649 |
| tsv_rows.emit/truth/truth_genome.tsv | 62 | 62 | 62 |
| tsv_rows.emit/truth/truth_igr.tsv | 92090 | 75140 | 79887 |
| tsv_rows.emit/truth/truth_module.tsv | 175 | 175 | 175 |
| tsv_rows.emit/truth/truth_rna.tsv | 1110 | 1022 | 1215 |
| tsv_rows.emit/truth/truth_spacer_target.tsv | 3887 | 3695 | 3781 |
| tsv_rows.emit/truth/truth_traps.tsv | 213 | 295 | 207 |
| tsv_rows.fusions/edits.tsv | 81 | 73 | 73 |
| tsv_rows.fusions/events.tsv | 22 | 22 | 22 |
| tsv_rows.fusions/fusions.tsv | 66 | 66 | 66 |
| tsv_rows.fusions/genes.tsv | 90242 | 73324 | 77664 |
| tsv_rows.fusions/nogd.tsv | 15 | 7 | 7 |
| tsv_rows.fusions/nogd_families.tsv | 2 | 2 | 2 |
| tsv_rows.fusions/replicons.tsv | 84 | 96 | 83 |
| tsv_rows.mobile/arrays.tsv | 30 | 30 | 30 |
| tsv_rows.mobile/cas.tsv | 40 | 40 | 40 |
| tsv_rows.mobile/genes.tsv | 90280 | 73362 | 77702 |
| tsv_rows.mobile/insertions.tsv | 46 | 46 | 46 |
| tsv_rows.mobile/phage_genes.tsv | 730 | 730 | 730 |
| tsv_rows.mobile/phages.tsv | 12 | 12 | 12 |
| tsv_rows.mobile/prophages.tsv | 16 | 16 | 16 |
| tsv_rows.mobile/replicons.tsv | 84 | 96 | 83 |
| tsv_rows.mobile/spacer_targets.tsv | 3462 | 3363 | 3363 |
| tsv_rows.mobile/spacers.tsv | 316 | 316 | 316 |
| tsv_rows.replicons/elements.tsv | 171872 | 153623 | 161405 |
| tsv_rows.replicons/emissions.tsv | 252 | 224 | 248 |
| tsv_rows.replicons/genes.tsv | 89102 | 72184 | 76524 |
| tsv_rows.replicons/junctions.tsv | 89102 | 72184 | 76524 |
| tsv_rows.replicons/lineages.tsv | 117 | 117 | 117 |
| tsv_rows.replicons/links.tsv | 20 | 20 | 20 |
| tsv_rows.replicons/models.tsv | 10 | 10 | 10 |
| tsv_rows.replicons/replicons.tsv | 84 | 96 | 83 |
| tsv_rows.replicons/rna.tsv | 1078 | 984 | 1159 |
| tsv_rows.sequences/domains.tsv | 4329 | 6655 | 4344 |
| tsv_rows.sequences/families.tsv | 2698 | 4208 | 2715 |
| tsv_rows.sequences/genes.tsv | 89102 | 72184 | 76524 |
| tsv_rows.sequences/lineages.tsv | 117 | 117 | 117 |
| tsv_rows.sequences/panel.tsv | 300 | 300 | 300 |
| tsv_rows.sources/contam.tsv | 10 | 10 | 10 |
| tsv_rows.sources/dup.tsv | 195 | 225 | 248 |
| tsv_rows.sources/mutations.tsv | 3135 | 3037 | 3488 |
| tsv_rows.sources/sidecars/fcs_gx_mock.tsv | 10 | 10 | 10 |
| tsv_rows.sources/source_contigs.tsv | 2315 | 2076 | 2231 |
| tsv_rows.sources/sources.tsv | 50 | 50 | 50 |
| world_digest | f0262027e4ccdcbc72cd7d90f23fbb5c | 765b20a0d7d0051d62c5036143d6ced5 | 684e539e11bd1ecbc5d3df6fe6405e7f |

### 15.4 The committed preset over the 22 seeds (scan `final12`)

`final12.tsv`, `final12_grid.tsv` (*measured*, the scan command above with no overrides, `--n0`; `tablegen.py`).
Session 11's `final.tsv` columns plus the fewest transfers of any decoy and both plasmid types; N0's size ratio is in
the N0 column. Ranges over the 22 seeds: chromosome mean 1,403-1,815 (22 of 22 inside 1,050-1,950); planted +
profile-only + decoy share 5.4-8.3%; copies per present member 0.75-1.11; module transfers into a holder 0 at every
seed; fewest decoy transfers 14-274 (>= 3 at 22 of 22); both plasmid types at 14 of 22 (the default among them), a
genome without a plasmid at 21 of 22, more than 2 in a genome at 0; N0 size ratio 0.893-1.086 (in [0.9, 1.1] at 21 of
22; 0.893 at 10975598825187167680); N0 r -0.08 to 0.97; all nine cells powered at 22 of 22; level medians strictly
increasing at 19 of 22; root families in >= 40 genomes 175-406.

| seed | crown age (MRCA) / from first speciation | genomes | genes per genome mean (min-max) | chromosome mean (min-max) | planted / profile-only / decoy / vehicle / hub / plasmid / background share % | plasmid genomes pls0, pls1 (0/1/2 per genome) | modules all members in >= 9 genomes | copies per present member | DEC-15 (c) plus, kept, clustered; scored units intact | realized transfers per member, level medians | N0 background size ratio, r | module transfers into a holder | fewest decoy transfers | both plasmid types |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 20261112 | 5 / 5 | 46 | 1663.6 (1432-1921) | 1570.5 (1321-1801) | 4.7 / 0.5 / 2.5 / 8.4 / 3.0 / 5.6 / 75.4 | 35, 2 (9/37/0) | 23 | 1.075 | 0.460, 0.0000, 0.021; 760/760 | 7.0, 7.5, 19.0 | 0.9457, 0.4382 | 0 of 1198 | 220 | yes |
| 20261113 | 5 / 8.435 | 46 | 1531.3 (1304-1768) | 1436.6 (1208-1648) | 4.3 / 0.6 / 2.5 / 7.3 / 3.2 / 6.2 / 75.9 | 36, 1 (9/37/0) | 21 | 1.047 | 0.459, 0.0019, 0.021; 644/644 | 4.0, 4.4, 16.67 | 1.0261, 0.5746 | 0 of 840 | 229 | yes |
| 1001685937451957930 | 5 / 5 | 46 | 1576.6 (1338-1811) | 1507.9 (1338-1810) | 4.2 / 0.5 / 2.5 / 5.0 / 3.1 / 4.4 / 80.4 | 24, 7 (19/23/4) | 20 | 0.976 | 0.524, 0.0048, 0.021; 720/720 | 2.0, 3.5, 11.0 | 1.0012, 0.4935 | 0 of 573 | 33 | yes |
| 6329801112132740259 | 5 / 5.136 | 46 | 1864.6 (1390-2154) | 1808.1 (1350-2114) | 3.8 / 0.5 / 2.2 / 7.4 / 2.6 / 3.0 / 80.4 | 11, 32 (8/33/5) | 21 | 1.088 | 0.477, 0.0010, 0.032; 627/627 | 7.0, 9.0, 20.0 | 0.942, 0.8717 | 0 of 1187 | 252 | yes |
| 6036017085589997374 | 5 / 5 | 46 | 1679.5 (1339-1875) | 1663.9 (1339-1835) | 5.2 / 0.3 / 2.8 / 7.9 / 2.9 / 0.9 / 79.9 | 0, 18 (28/18/0) | 25 | 1.066 | 0.412, 0.0019, 0.021; 884/884 | 5.5, 8.8, 16.0 | 0.9053, 0.7889 | 0 of 1095 | 274 | no |
| 4383110225690451409 | 5 / 5 | 46 | 1809.3 (1615-2023) | 1678.9 (1496-1863) | 4.2 / 0.6 / 2.2 / 8.5 / 2.7 / 7.2 / 74.5 | 39, 33 (0/20/26) | 25 | 1.084 | 0.432, 0.0010, 0.065; 786/786 | 3.33, 5.33, 19.0 | 1.0082, 0.2193 | 0 of 924 | 209 | yes |
| 10975598825187167680 | 5 / 5 | 46 | 1727.2 (1527-2046) | 1645.4 (1483-1926) | 4.8 / 0.5 / 2.3 / 7.4 / 2.9 / 4.7 / 77.5 | 28, 10 (12/30/4) | 24 | 1.106 | 0.522, 0.0029, 0.021; 764/764 | 2.33, 6.0, 10.0 | 0.8927, 0.6439 | 0 of 762 | 169 | yes |
| 5467273169533009601 | 5 / 5.03 | 46 | 1589.1 (1377-1779) | 1515.2 (1298-1659) | 4.2 / 0.6 / 2.6 / 7.7 / 3.1 / 4.7 / 77.2 | 24, 13 (13/29/4) | 21 | 0.952 | 0.426, 0.0019, 0.043; 645/645 | 4.0, 5.0, 11.5 | 0.9718, 0.0351 | 0 of 807 | 138 | yes |
| 4734920099712499969 | 5 / 5 | 46 | 1548.3 (1229-1729) | 1448.3 (1156-1579) | 5.4 / 0.7 / 1.8 / 5.4 / 3.2 / 6.5 / 77.1 | 33, 16 (9/25/12) | 23 | 1.112 | 0.439, 0.0038, 0.032; 814/814 | 4.0, 6.33, 13.33 | 1.086, 0.8333 | 0 of 923 | 14 | yes |
| 1570278518166000558 | 5 / 6.18 | 46 | 1823.4 (1327-2004) | 1814.7 (1327-2004) | 4.1 / 0.3 / 1.8 / 8.0 / 2.7 / 0.5 / 82.6 | 0, 10 (36/10/0) | 23 | 1.024 | 0.454, 0.0010, 0.054; 786/786 | 5.67, 7.0, 16.0 | 0.9818, 0.33 | 0 of 947 | 87 | no |
| 3830344756309592292 | 5 / 5 | 46 | 1685.3 (1418-1950) | 1637.5 (1299-1843) | 3.5 / 0.5 / 2.1 / 4.1 / 2.9 / 2.8 / 84.0 | 15, 10 (25/17/4) | 21 | 1.085 | 0.435, 0.0010, 0.043; 603/603 | 4.0, 6.33, 10.0 | 0.9427, 0.6206 | 0 of 696 | 150 | yes |
| 9254619518473077140 | 5 / 7.657 | 46 | 1684.3 (1389-1886) | 1630.4 (1356-1824) | 3.9 / 0.6 / 2.4 / 9.5 / 2.9 / 3.2 / 77.5 | 17, 11 (20/24/2) | 21 | 1.089 | 0.525, 0.0010, 0.043; 645/645 | 5.0, 6.4, 16.0 | 0.971, 0.8377 | 0 of 1026 | 217 | yes |
| 4701969716356831982 | 5 / 5.087 | 46 | 1426.8 (1289-1725) | 1403.3 (1289-1565) | 4.5 / 0.6 / 2.7 / 2.4 / 3.5 / 1.7 / 84.8 | 2, 21 (25/19/2) | 16 | 0.922 | 0.428, 0.0000, 0.021; 613/613 | 3.17, 2.0, 8.67 | 0.9695, 0.6317 | 0 of 507 | 21 | yes |
| 14843108908847119356 | 5 / 5 | 46 | 1740.6 (1485-2020) | 1685.8 (1485-1900) | 4.2 / 0.5 / 2.1 / 6.3 / 2.9 / 3.1 / 80.9 | 21, 0 (25/21/0) | 23 | 1.031 | 0.550, 0.0029, 0.054; 792/792 | 6.0, 9.0, 12.5 | 1.0172, 0.6274 | 0 of 992 | 103 | no |
| 1758546030705738202 | 5 / 5 | 46 | 1737.1 (1041-1956) | 1737.1 (1041-1956) | 3.1 / 0.7 / 1.6 / 7.1 / 2.8 / 0.0 / 84.7 | 0, 0 (46/0/0) | 16 | 0.753 | 0.563, 0.0000, 0.032; 543/543 | 1.5, 1.0, 5.0 | 0.9575, 0.9733 | 0 of 379 | 43 | no |
| 1995372076782357328 | 5 / 5 | 46 | 1691.8 (1534-1874) | 1646.6 (1494-1790) | 4.1 / 0.5 / 1.9 / 9.8 / 2.9 / 2.7 / 78.2 | 5, 37 (9/32/5) | 20 | 0.966 | 0.483, 0.0048, 0.043; 736/736 | 2.83, 8.83, 12.75 | 0.9617, 0.8447 | 0 of 759 | 35 | yes |
| 13514522956726874792 | 5 / 9.488 | 46 | 1592.5 (1362-1758) | 1576.8 (1362-1758) | 4.5 / 0.6 / 2.4 / 11.4 / 2.5 / 1.0 / 77.6 | 0, 18 (29/16/1) | 24 | 0.984 | 0.444, 0.0019, 0.054; 754/754 | 4.17, 5.0, 12.33 | 0.973, 0.1002 | 0 of 829 | 189 | no |
| 12621855678475765395 | 5 / 7.443 | 46 | 1640.1 (1386-1767) | 1601.0 (1386-1734) | 3.8 / 0.6 / 2.0 / 4.8 / 3.0 / 2.4 / 83.4 | 12, 9 (26/19/1) | 19 | 1.041 | 0.468, 0.0019, 0.032; 619/619 | 1.17, 4.0, 7.67 | 0.9447, -0.0787 | 0 of 507 | 84 | yes |
| 7689157899245528466 | 5 / 5 | 46 | 1530.2 (1221-1749) | 1504.9 (1204-1709) | 5.1 / 0.6 / 2.6 / 4.2 / 3.2 / 1.7 / 82.7 | 1, 26 (20/25/1) | 23 | 1.077 | 0.546, 0.0019, 0.010; 758/758 | 6.0, 10.0, 18.0 | 0.9371, 0.5681 | 0 of 1439 | 224 | yes |
| 13748185708697958051 | 5 / 5 | 46 | 1604.2 (1336-1860) | 1578.1 (1265-1860) | 4.5 / 0.6 / 1.9 / 4.6 / 3.1 / 1.6 / 83.6 | 10, 0 (36/10/0) | 21 | 1.028 | 0.550, 0.0010, 0.021; 642/642 | 6.0, 5.2, 15.67 | 0.9154, 0.6402 | 0 of 1015 | 18 | no |
| 15390672246448918260 | 5 / 6.455 | 46 | 1570.8 (1433-1676) | 1553.4 (1433-1676) | 4.2 / 0.7 / 2.0 / 6.4 / 3.2 / 1.1 / 82.5 | 0, 20 (27/18/1) | 22 | 0.933 | 0.416, 0.0000, 0.021; 647/647 | 2.67, 6.0, 11.0 | 1.0021, 0.4566 | 0 of 927 | 24 | no |
| 11714203776604047358 | 5 / 5.052 | 46 | 1488.5 (1348-1729) | 1485.9 (1348-1683) | 4.8 / 0.6 / 1.8 / 9.4 / 3.3 / 0.2 / 79.9 | 0, 3 (44/1/1) | 24 | 1.023 | 0.473, 0.0029, 0.021; 729/729 | 3.5, 9.67, 16.33 | 0.996, 0.2587 | 0 of 912 | 88 | no |

| member transfer \ together | 0.3 | 0.7 | 1 | level, pooled |
|---|---|---|---|---|
| 0.01 | 4.0 (0-13.6); together 0.10; powered at 22/22 | 4.2 (0-15.5); together 0.39; powered at 22/22 | 4.0 (0-13); together 1.00; powered at 22/22 | 4.0 |
| 0.03 | 5.1 (0-25.67); together 0.13; powered at 22/22 | 6.1 (0-16.75); together 0.38; powered at 22/22 | 7.0 (0-24); together 1.00; powered at 22/22 | 6.0 |
| 0.4 | 11.2 (0-47.5); together 0.24; powered at 22/22 | 12.4 (0-20.25); together 0.69; powered at 22/22 | 17.0 (5-31); together 1.00; powered at 22/22 | 12.4 |

### 15.5 Test runs on the committed tree

All *measured*, wilkcomp02, 2026-11-15, image `bio-box-test_faec1f4.sif`:

- **Test that failed first**: `tests/r2/test_synthworld_backbone.py::test_decoy_transfers_land_only_where_the_partner_is`
  on the merged branch before any change (decoy transfer 0.4): `AssertionError: Counter({(1298, True): 54, (1296, True):
  24, (1299, True): 12, (1300, True): 6, (1297, True): 2})`, "1 failed, 23 deselected in 38.19s" (`first.log`; a
  detached unit on a mirror of the tree, `python -m pytest -q tests/r2/test_synthworld_backbone.py -k
  decoy_transfers_land_only_where_the_partner_is`). After the change: decoys 220, 365, 410, 432 and 293 transfers, all
  where the partner is, mean partner share 0.786 (bound 0.9).
- **Every synthworld test outside the null module and `tests/r0`**: "1 failed, 503 passed, 6 skipped in 2286.64s"
  (`env BIOBOX_TREE=<worktree> sh scripts/pod_pytest.sh F-02-12 -s tests/r0 tests/r1/test_synthworld_*.py
  tests/r2/test_synthworld_{backbone,preset_seeds,sequences,sequences_pfam,intergenic,rna_rfam,mobile,fusions,assembly,emit}.py`;
  `tests_final12.txt`). The one failure is §15.6 (B). Passing among them: the migrated motif test (1,150 of 1,159
  instances the unique best window, 4 accepted as another planted instance of the emission), the cross-phylum IGR test
  (206 pairs, burst identity 0.946), the hub test (181 root families in >= 40 genomes, hubs 304-342 neighbour families
  against a root median of 60), DEC-15 (a) with fates "family", the plasmid replicon test (9 genomes without a plasmid,
  none with more than 1), the decoy test, every per-seed test (both plasmid types at 5 of 6 seeds, the default among
  them; levels pooled 4.08, 6.21, 16.0), and every re-pinned digest.
- **The null module** at the default seed: "1 failed, 13 passed in 1296.44s"; at 20261113: "1 failed, 12 passed, 1
  skipped in 1267.37s" (the pinned-digest test skips off the default seed) (`session12/scripts/nulls.sh <seed> 200`, a
  detached unit at 2 CPUs: `python -m pytest -q -s tests/r2/test_synthworld_nulls.py` with
  `BIOBOX_SYNTHWORLD_NULL_SEED`; `nulls-20261112.log`, `nulls-20261113.log`). The one failure at each seed is §15.6
  (A). N1's label pairs pass at both seeds: default 37 of 100 members relabeled, 83 pairs, p 0.992 (`n1`) and 1.000
  (`n1_s2`); 20261113 45 relabeled, 90 pairs, p 0.760 and 0.174. The naive co-occurrence statistic stays
  miscalibrated on both N0 worlds (27.0x and 25.5x the calibrated share at p < 1e-3, default seed). N0's size ratio:
  0.946 at the default seed (from the scan; the test stops at the KS line before printing), 1.026 at 20261113, r 0.438
  and 0.575 (printed only).
- **The whole suite** (`bash /tmp/if2453/ci/F-02-12/scripts/orchestrator/detached_suite.sh F-02-12 2 --deselect
  tests/r2/test_synthworld_nulls.py`, the null module run above): shard 0 "1201 passed, 12 skipped, 18 warnings in
  3979.98s", shard 1 "1 failed, 1101 passed, 12 skipped, 14 deselected, 48 warnings in 6791.63s" (the failure is
  §15.6 (B)). GitHub Actions `ci` on `task/F-02-12` at the note's first commit (run 43771558605, Linux): "2 failed, 2110
  passed, 226 skipped", the two of §15.6 with the same values as on the pod (D 0.0917; `feature_uid`
  a2873401f69e506e1a0f1d0923da0992), so the default world is the same on both hosts.

### 15.6 What fails, and why (stop and report; no assertion changed)

**(A) N0's root-presence KS** (`test_n0_keeps_the_background_process_in_distribution`, the line after the migrated
one: two-sample KS of the genomes-present counts of the root families, D below 1.628 sqrt((n + m) / (n m)), the
alpha = 0.01 critical value). Under the committed preset it fails at both seeds the null module runs: default seed
N0-full D 0.0917 against 0.0728 (N0-strict 0.034 passes), 20261113 N0-strict D 0.0785 against 0.0736 (N0-full 0.052
passes) (*measured*, the null module, `nulls-2026111{2,3}.log`). Session 11 never reached this line at its presets:
the r >= 0.9 assertion before it failed first, so a failure here was masked until DEC-44 (4) removed r. Over the 22
seeds (scan `--n0`, which now also builds N0-strict and computes this KS from the same `families.tsv` columns; it
reproduces the test's D at both seeds; `kssumm.py`; all *measured*, wilkcomp02, 2026-11-15 15:40-17:15 CST):

| scan | N0-full KS below bound | N0-strict KS below bound | default: full D, strict D (bound) | 20261113: full D, strict D (bound) | null - signal mean root presence, full; strict | N0-full size ratio (in [0.9, 1.1]) | decoys >= 3 transfers: seeds, fewest |
|---|---|---|---|---|---|---|---|
| n_dec0.4 (committed, decoy 0.4) | 22/22 | 22/22 | 0.0556, 0.0266 (0.0729) | 0.0561, 0.064 (0.0733) | -1.99 to 1.37; -1.75 to 1.34 | 0.946-1.140 (21/22) | 9/22, 0 |
| n_s11 (session 11's `89cfb46` preset) | 19/22 | 19/22 | 0.0363, 0.0422 (0.0728) | 0.0343, 0.0195 (0.0731) | -1.02 to 1.92; -1.39 to 2.54 | 0.973-1.091 (22/22) | 7/22 (session 11) |
| **n_final12 (committed, decoy 8.0)** | 19/22 | 18/22 | **0.0917**, 0.034 (0.0728) | 0.052, **0.0785** (0.0736) | -2.21 to 1.82; -2.19 to 2.03 | 0.893-1.086 (21/22) | 22/22, 14 |
| n_dec12.0 | 21/22 | 21/22 | 0.0607, 0.027 (0.0728) | 0.0499, 0.0297 (0.0733) | -2.30 to 1.33; -1.57 to 1.60 | 0.899-1.079 (21/22) | 20/22, 0 |
| n_dec16.0 | 22/22 | 19/22 | 0.0685, 0.0331 (0.0728) | 0.0384, 0.0614 (0.0733) | -1.53 to 1.04; -2.19 to 1.37 | 0.862-1.060 (21/22) | 22/22, 22 |

Every rate changes the whole draw, so neither column is monotone in the decoy rate: 12.0 loses the decoy test at two
seeds where 8.0 held it. The mean presence of a root family is not shifted one way (null minus signal -2.3 to +1.9
genomes at every rate), so this reads as the KS being sensitive to one world's draw (the 1,000 root families of one
world share one tree and are not independent samples, so the alpha = 0.01 bound is exceeded at 0 to 4 of 22 seeds
rather than at about 1%) more than as a bias from the decoys. 16.0 holds the decoy test at every seed and the KS at both
seeds the module runs, so it was tried as the clause allows; its full default world (`igr_dec16.log`, background
transfer 0.045) has **193** equal-length cross-phylum IGR pairs, below DEC-44 (5)'s 200, so 16.0 would reopen the
transfer refit, whose every step redraws every one of these per-seed properties again. The session stops here rather
than search further draws: the decoy rate that holds the decoy test at every seed with the other DEC-44 clauses (8.0)
breaks an existing null assertion at the null module's seeds, which is DEC-44 (2)'s stop. The preset stays at 8.0
(the smallest rate that holds the clause), the test is not touched. For the PI: (i) keep decoy 8.0 and decide the KS
(a bound over the 22 seeds, or a presence test that conditions on the tree), or (ii) accept a decoy rate that does not
hold at every seed (4.0: 21 of 22, fewest 1), or (iii) choose decoy partners among modules at the higher member transfer
levels (the preset tags `partner_modules` as fitted on session 2's default seed; a decoy fails when its partner module
is soon held by one lineage, §15.1), which should make a modest rate hold everywhere.

**(B) Two truth genes with one `feature_uid`** (`tests/r2/test_synthworld_emit.py::test_every_truth_gene_reads_its_s8_gene_on_the_emitted_contig`,
assertion `r["feature_uid"] not in keys`). In the default world under the final preset, genome g036 carries an exact,
young tandem triplication of a two-gene segment (F000364, F000733: `assembly/repeats.tsv`); S9 breaks the chromosome at
repeats, so contigs `g036_c0058` (830,022-831,449) and `g036_c0060` (832,931-834,358) are 1,427 bp each with
identical sequence (`assembly/contigs.tsv`: both `contig_uid` 90a034d056ccbea31575afc643504c8f). `contig_uid` is
BLAKE3 of the sequence, so both contigs are one contig to the truth (and to `make_inputs`' "distinct contigs"), and the
3' halves of the copies `n105_g238748` and `n105_g238750` both sit at [0, 708) on + on it: one `feature_uid`, two
truth rows (*measured*, `world-final/emit/truth/truth_gene.tsv` lines 64,303-64,304). DEC-27 and the emit docstring
cover an exact duplicate across genomes ("an exact duplicate shares the contig_uid, so its rows are the original's")
but not identical contigs inside one genome. It is a property of this world's draw, not of a rate: duplication 0.012
makes young exact tandem copies, and S9's repeat breaks cut them into identical contigs. What the truth should do
here (one row with both gene ids, or one row and a note of the copies, or S9 joining identical contigs) is a design
decision for `emit.py`/DEC-27, not this card's; the test is not changed.

### 15.7 What remains of card F-02

- The two decisions of §15.6: (A) the null module's root-presence KS against decoy 8.0 (keep the rate and decide the
  test, accept a rate that misses some seeds, or move the decoy partners), and (B) how the truth keys genes on two
  identical contigs of one genome (`emit.py`/DEC-27). Until they are decided the branch is not mergeable on its own
  tests: §15.5 has exactly these two failures (one in each run), everything else passes.
- Session 13 (DEC-39): detectability against divergence in the sequence model.
- §4: the §4.3 ablation worlds and build-time gate, the 20 N2 worlds; §6.2 to §6.4 (macOS against Linux digests, a
  CLI, `test_run` columns) and §7 (realism against the golden genomes, tool parity), as the card lists.
- The spec's plasmid paragraph (§2.1) and §5.1's backbone paragraph still give session 2's rates; the preset and
  DEC-44 supersede them (not edited here).
