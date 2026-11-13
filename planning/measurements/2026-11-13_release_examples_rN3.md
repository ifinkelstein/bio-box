# Worked examples on `rN.3` (release check, AGENTS 5.7)

*Measured*, wilkcomp02, started 2026-11-13 14:06:32 CST, by `scripts/release_examples.sh rN.3`
(run mode): tree `/tmp/if2453/e11/tree` (commit `8ea227f` of `task/E-11`, mirrored with `rsync`), the runtime image
`/tmp/if2453/images/bio-box_64f1bf0.sif` (`sha256:cb7353f85316981be5268c2212100322b5fbf1086c32b07e2f9b47802b6afba3`,
RUNBOOK section 0; matplotlib 3.11.2), data root `/tmp/if2453/r10/root` (read-only), outputs `/tmp/if2453/e11/acc`, 2
threads per launch, started as the unit `e11-acc-driver` (`MemoryMax=2G`, nice 19), each step its own unit
`e11-acc-<example>-<step>-<pid>`. The card's acceptance run (card E-11; AGENTS §5.7): the first run of the rule. This
file is the script's own `note.md`, with the reading written in where it said TODO.

## Result

Every invariant of the 2 examples holds (exit 0, last line `release_examples.sh: rN.3: 2 examples, every invariant
holds`), and every recorded value equals its `recorded/rN.3.tsv`: 0 `differs`, 0 `new`, 0 `missing`.

- **CRISPR-Cas** (`examples/crispr_cas_neighbourhoods/`): 19 invariants and 31 recorded values ok. The run's
  `observed.tsv` is byte-identical to the example's committed `results/observed.tsv` (`cmp`), and so are
  `tniq_cas_loci.tsv`, the 10 kb summary (`results/cas7f_10kb.tsv`), `uids_tniq_cas.txt` and the two figures that are
  committed from `run.sh`'s steps (`figures/tniq_cas_pfam.svg`, `figures/cas7f_all_nr50.svg`; `sha256sum`). That run
  (`/tmp/if2453/ex2/rN.3-b`) used the test image `bio-box-test_64f1bf0.sif` and the launcher inside `run.sh`; this
  one the runtime image and the shared launcher (`examples/_shared/launch.sh`). Its committed `results/check.txt`
  differs from this run's `check.txt` only in the recorded rows, `new` there (no `recorded/rN.3.tsv` existed yet) and
  `ok` here.
- **RNA polymerase** (`examples/rnap/`, new): 18 invariants and 24 recorded values ok. Its recorded values were taken
  from a first run the same day (`/tmp/if2453/e11/rnap1`, 13:47 to 14:05, `recorded_new.tsv` copied); this run's
  `observed.tsv` is byte-identical to that run's, and so are `hits.tsv`, `hits.faa`, `neighbors.tsv`, the 2,581 locus
  tables and the 432 `nbr flanks` groups with their figures (`diff -r -q`: only `RUN` and the HMM path in
  `recipe.yaml` differ, as the run id and `OUT` differ).

**Against the measurement notes** (the RNA polymerase numbers, all on the same genomes): 1,957 proteins (1,931 in
phages, 26 in prokaryotes) on 2,581 genes, 2,526 contigs, 463 censored loci, 24,473 neighbours, T7 gene 1 at 367.7
bits (`2026-11-11_trial_search_rN1_second.md` section 3; `2026-11-11_build_r50.md` §A3); 2,569 hit genes with a Pfam
architecture and `Phage_T7_1_1` beside 1,032 loci (§A3); 432 flank groups, 1,558 leaves, 598 context types, none
holding phage and prokaryote loci (`2026-11-12_flanks_rN3.md`). Every one is reproduced; none differs. One reading to
note: T7 gene 1's own flank group is `c1bf639f...` (162 loci, 10 types), not the note's "T7-like group" `4a0062a1...`
(306 loci, 11 types), which is another NR90 cluster of T7-like polymerases and is present here too.

**Cost against the notes.** CRISPR-Cas: 1,371.2 CPU-seconds (22.85 CPU-minutes) against 1,366 (22.8) in its README
(+0.4%). RNA polymerase: 1,106.4 CPU-seconds (18.4 CPU-minutes); its search unit 11.25 CPU-minutes in the first run
against 10.78 in §A3 (+4%), `nbr flanks --figure` 4.66 against 4.57 (`2026-11-12_flanks_rN3.md`, W-27; +2%)
(`examples/rnap/README.md` section 4). The run's wall time with 2 search CPUs, 394.1 s, is 40% above §A3's 280.5 s
with 4: not a like-for-like number. Wall time of the whole check: 56 min 55 s (*derived*: 14:06:32 to 15:03:27 and the
checking step), of which at least 26 minutes waiting for room under the day cap (*derived*: 14 and 12 `waiting for
room` lines, 60 s each or more): the account held 38 to 41 CPUs for most of the hour. With room it would take about
half an hour (*estimate*: 57 - 26 min).

**What the examples found missing in bio-box** is in their READMEs (CRISPR section 5; nothing new from the RNA
polymerase example). Not run here: the 50 Gbp release `rM.2`, which has no search pack yet; card E-11 asks for it
when it exists and verifies.

## Per example

| Example | run.sh exit | check exit | ok | FAIL | differs | new | missing | unexpected | CPU-s (*derived*) |
|---|---|---|---|---|---|---|---|---|---|
| `crispr_cas_neighbourhoods` | 0 | 0 | 50 | 0 | 0 | 0 | 0 | 0 | 1371.2 |
| `rnap` | 0 | 0 | 42 | 0 | 0 | 0 | 0 | 0 | 1106.4 |

CPU-s is the sum of the units' `CPUUsageNSec` (`pod_run: ... usage` lines) and of `/usr/bin/time -v` user
plus system time of the scopes, from the example's logs. Wall: CRISPR-Cas 30 min 56 s (14:06:32 to 14:37:28), RNA
polymerase 25 min 59 s (14:37:28 to 15:03:27), from the driver's log.

## Differences from the recorded values

None.

## Every check

| Example | Check | Kind | Expected | Observed | Verdict |
|---|---|---|---|---|---|
| `crispr_cas_neighbourhoods` | `tniq_control_bits` | invariant | >= 27 | 61.3 | ok |
| `crispr_cas_neighbourhoods` | `cas7f_control_bits` | invariant | >= 25 | 138.5 | ok |
| `crispr_cas_neighbourhoods` | `tniq_control_order` | invariant | == IMPDH tnsA+tnsB tniB . . . . . . . ^ cas5f cas7f cas6f . . . . . . . | IMPDH tnsA+tnsB tniB . . . . . . . ^ cas5f cas7f cas6f . . . . . . . | ok |
| `crispr_cas_neighbourhoods` | `pdam_complete_if_window_has_cas1_and_cas3` | invariant | == true | true | ok |
| `crispr_cas_neighbourhoods` | `icp1_cas7f_window_has_cas1_and_cas3` | invariant | == true | true | ok |
| `crispr_cas_neighbourhoods` | `tniq_cas7f_loci_with_cas1_cas2_cas3` | invariant | == 0 | 0 | ok |
| `crispr_cas_neighbourhoods` | `cas7f_tniq_loci_with_cas1_cas2_cas3` | invariant | == 0 | 0 | ok |
| `crispr_cas_neighbourhoods` | `bp_tniq_cas7f_loci_with_cas1_cas2_cas3` | invariant | == 0 | 0 | ok |
| `crispr_cas_neighbourhoods` | `bp_cas7f_tniq_loci_with_cas1_cas2_cas3` | invariant | == 0 | 0 | ok |
| `crispr_cas_neighbourhoods` | `tniq_cas7f_share_cas5f_at_plus1` | invariant | >= 0.5 | 0.8 | ok |
| `crispr_cas_neighbourhoods` | `tniq_pfam_proteins_not_hit` | invariant | == 0 | 0 | ok |
| `crispr_cas_neighbourhoods` | `cas7f_pfam_proteins_not_hit` | invariant | == 0 | 0 | ok |
| `crispr_cas_neighbourhoods` | `max_bits_diff_search_vs_pfam` | invariant | <= 0.1 | 0.0 | ok |
| `crispr_cas_neighbourhoods` | `oracle_tniq_only_search` | invariant | == 0 | 0 | ok |
| `crispr_cas_neighbourhoods` | `oracle_tniq_only_hmmsearch` | invariant | == 0 | 0 | ok |
| `crispr_cas_neighbourhoods` | `oracle_tniq_max_bits_diff` | invariant | <= 0.1 | 0.00 | ok |
| `crispr_cas_neighbourhoods` | `oracle_cas7f_only_search` | invariant | == 0 | 0 | ok |
| `crispr_cas_neighbourhoods` | `oracle_cas7f_only_hmmsearch` | invariant | == 0 | 0 | ok |
| `crispr_cas_neighbourhoods` | `oracle_cas7f_max_bits_diff` | invariant | <= 0.1 | 0.00 | ok |
| `crispr_cas_neighbourhoods` | `tniq_hit_proteins` | recorded | 587 | 587 | ok |
| `crispr_cas_neighbourhoods` | `cas7f_hit_proteins` | recorded | 114 | 114 | ok |
| `crispr_cas_neighbourhoods` | `tniq_hits_not_in_pfam` | recorded | 13 | 13 | ok |
| `crispr_cas_neighbourhoods` | `cas7f_hits_not_in_pfam` | recorded | 0 | 0 | ok |
| `crispr_cas_neighbourhoods` | `tniq_loci_virus` | recorded | 0 | 0 | ok |
| `crispr_cas_neighbourhoods` | `tniq_loci_not_virus` | recorded | 635 | 635 | ok |
| `crispr_cas_neighbourhoods` | `cas7f_loci_virus` | recorded | 46 | 46 | ok |
| `crispr_cas_neighbourhoods` | `cas7f_loci_not_virus` | recorded | 111 | 111 | ok |
| `crispr_cas_neighbourhoods` | `tniq_cas7f_loci` | recorded | 20 | 20 | ok |
| `crispr_cas_neighbourhoods` | `tniq_cas7f_loci_virus` | recorded | 0 | 0 | ok |
| `crispr_cas_neighbourhoods` | `tniq_cas7f_contigs` | recorded | 18 | 18 | ok |
| `crispr_cas_neighbourhoods` | `tniq_cas7f_genomes` | recorded | 16 | 16 | ok |
| `crispr_cas_neighbourhoods` | `tniq_cas7f_loci_with_tnsA_tnsB_tnsC` | recorded | 6 | 6 | ok |
| `crispr_cas_neighbourhoods` | `tniq_cas7f_loci_with_cas8f_csy1` | recorded | 0 | 0 | ok |
| `crispr_cas_neighbourhoods` | `tniq_minimal_non_if_loci` | recorded | 3 | 3 | ok |
| `crispr_cas_neighbourhoods` | `tniq_minimal_ib_loci` | recorded | 0 | 0 | ok |
| `crispr_cas_neighbourhoods` | `tniq_minimal_ic_loci` | recorded | 3 | 3 | ok |
| `crispr_cas_neighbourhoods` | `tniq_other_cas_loci` | recorded | 4 | 4 | ok |
| `crispr_cas_neighbourhoods` | `cas7f_not_virus_loci_without_cas1_cas3` | recorded | 22 | 22 | ok |
| `crispr_cas_neighbourhoods` | `cas7f_not_virus_loci_without_cas1_cas3_with_tniq` | recorded | 19 | 19 | ok |
| `crispr_cas_neighbourhoods` | `cas7f_virus_loci_with_cas1_and_cas3` | recorded | 46 | 46 | ok |
| `crispr_cas_neighbourhoods` | `tniq_cas7f_loci_near_IMPDH` | recorded | 1 | 1 | ok |
| `crispr_cas_neighbourhoods` | `tniq_cas7f_loci_near_4HBT` | recorded | 0 | 0 | ok |
| `crispr_cas_neighbourhoods` | `tniq_cas7f_loci_near_SIS` | recorded | 0 | 0 | ok |
| `crispr_cas_neighbourhoods` | `tniq_cas7f_loci_near_Mg_chelatase` | recorded | 0 | 0 | ok |
| `crispr_cas_neighbourhoods` | `bp_tniq_loci` | recorded | 622 | 622 | ok |
| `crispr_cas_neighbourhoods` | `bp_cas7f_loci` | recorded | 157 | 157 | ok |
| `crispr_cas_neighbourhoods` | `bp_tniq_cas7f_loci` | recorded | 21 | 21 | ok |
| `crispr_cas_neighbourhoods` | `bp_cas7f_tniq_loci` | recorded | 19 | 19 | ok |
| `crispr_cas_neighbourhoods` | `bp_tniq_cas7f_loci_with_tnsA_tnsB_tnsC` | recorded | 12 | 12 | ok |
| `crispr_cas_neighbourhoods` | `bp_tniq_cas7f_loci_with_tnsB` | recorded | 18 | 18 | ok |
| `rnap` | `t7_control_bits` | invariant | >= 30 | 367.7 | ok |
| `rnap` | `t7_control_architecture` | invariant | == RPOL_N+RNA_pol | RPOL_N+RNA_pol | ok |
| `rnap` | `t7_control_order` | invariant | == ocr Phage_T7_GP04 Gp0_5 Phage_T7_0_6 Phage_T7_kinase ^ Phage_T7_1_1 T7-like_gp12 DNA_ligase_A_M+DNA_ligase_C DUF5464 BC10 | ocr Phage_T7_GP04 Gp0_5 Phage_T7_0_6 Phage_T7_kinase ^ Phage_T7_1_1 T7-like_gp12 DNA_ligase_A_M+DNA_ligase_C DUF5464 BC10 | ok |
| `rnap` | `t7_window_is_the_contig_order` | invariant | == true | true | ok |
| `rnap` | `t7_positive_recovered` | invariant | == true | true | ok |
| `rnap` | `coord_t7_kinase_offset` | invariant | == -1 | -1 | ok |
| `rnap` | `coord_t7_ligase_offset` | invariant | == 3 | 3 | ok |
| `rnap` | `rnap_negative_control_hits` | invariant | == 0 | 0 | ok |
| `rnap` | `rnap_loci_virus` | invariant | >= 1 | 2555 | ok |
| `rnap` | `rnap_loci_not_virus` | invariant | >= 1 | 26 | ok |
| `rnap` | `rnap_loci_missing_windows` | invariant | == 0 | 0 | ok |
| `rnap` | `rnap_hits_below_cutoff` | invariant | == 0 | 0 | ok |
| `rnap` | `rnap_pfam_proteins_at_cutoff_not_hit` | invariant | == 0 | 0 | ok |
| `rnap` | `max_bits_diff_search_vs_pfam` | invariant | <= 0.1 | 0.0 | ok |
| `rnap` | `oracle_rnap_only_search` | invariant | == 0 | 0 | ok |
| `rnap` | `oracle_rnap_only_hmmsearch` | invariant | == 0 | 0 | ok |
| `rnap` | `oracle_rnap_max_bits_diff` | invariant | <= 0.1 | 0.00 | ok |
| `rnap` | `rnap_flank_loci_missing` | invariant | == 0 | 0 | ok |
| `rnap` | `rnap_hmm_sha256` | recorded | 32e6f0f0b1ef15d23d8ead23095f0b2a5a3323551b899bba99218d75f830e87e | 32e6f0f0b1ef15d23d8ead23095f0b2a5a3323551b899bba99218d75f830e87e | ok |
| `rnap` | `rnap_hit_proteins` | recorded | 1957 | 1957 | ok |
| `rnap` | `rnap_hit_genes` | recorded | 2581 | 2581 | ok |
| `rnap` | `rnap_hit_contigs` | recorded | 2526 | 2526 | ok |
| `rnap` | `rnap_hit_genomes` | recorded | 2526 | 2526 | ok |
| `rnap` | `rnap_proteins_virus` | recorded | 1931 | 1931 | ok |
| `rnap` | `rnap_proteins_not_virus` | recorded | 26 | 26 | ok |
| `rnap` | `rnap_genes_with_architecture` | recorded | 2569 | 2569 | ok |
| `rnap` | `rnap_hits_not_in_pfam` | recorded | 11 | 11 | ok |
| `rnap` | `rnap_pfam_proteins` | recorded | 1957 | 1957 | ok |
| `rnap` | `rnap_pfam_proteins_below_cutoff` | recorded | 11 | 11 | ok |
| `rnap` | `rnap_summary_loci` | recorded | 2581 | 2581 | ok |
| `rnap` | `rnap_summary_contigs` | recorded | 2526 | 2526 | ok |
| `rnap` | `rnap_summary_censored_loci` | recorded | 463 | 463 | ok |
| `rnap` | `rnap_summary_neighbours` | recorded | 24473 | 24473 | ok |
| `rnap` | `rnap_top_pfam_neighbour` | recorded | Phage_T7_1_1 | Phage_T7_1_1 | ok |
| `rnap` | `rnap_top_pfam_neighbour_loci` | recorded | 1032 | 1032 | ok |
| `rnap` | `rnap_flank_groups` | recorded | 432 | 432 | ok |
| `rnap` | `rnap_flank_leaves` | recorded | 1558 | 1558 | ok |
| `rnap` | `rnap_flank_context_types` | recorded | 598 | 598 | ok |
| `rnap` | `rnap_flank_types_mixed_kind` | recorded | 0 | 0 | ok |
| `rnap` | `t7_flank_group_loci` | recorded | 162 | 162 | ok |
| `rnap` | `t7_flank_group_types` | recorded | 10 | 10 | ok |
| `rnap` | `t7_flank_type_kinds` | recorded | virus | virus | ok |
