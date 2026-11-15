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
