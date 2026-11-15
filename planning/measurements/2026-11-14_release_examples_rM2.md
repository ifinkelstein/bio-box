# Worked examples on `rM.2` (release check, AGENTS 5.7)

*Measured*, wilkcomp02, started 2026-11-14 22:38:09 CST, by `scripts/release_examples.sh rM.2`
(run mode).

- Tree: `/tmp/if2453/w36/tree`, commit `a123da1` of `task/W-36`, mirrored with `rsync`. It holds card W-36's reader
  fix; without the fix this check stopped at its first step (`2026-11-14_partial_labels_rM2.md`).
- Image: the runtime image `/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif`
  (`sha256:d67300aa1321b571fa8da79e74dfccdd87a8dca2e56518f324902efa6b9f098a`, RUNBOOK section 0).
- Data root: `/tmp/if2453/r10/root`, read-only. Search pack: `26e10f05401faeefb78769e83f7f6cb4` (card W-35).
- Outputs: `/tmp/if2453/w36/ex-rM2`.
- Started as the unit `w36-ex-rM2` (`MemoryMax=2G`, nice 19), 2 threads per launch. Each step ran as its own unit,
  `w36-ex-rM2-<example>-<step>-<pid>`.

This file is the script's own `note.md`, with the reading written in where it said TODO.

## Result

Every invariant of the 2 examples holds. The exit was 0, and the last line was `release_examples.sh: rM.2: 2 examples,
every invariant holds; summary /tmp/if2453/w36/ex-rM2/summary.tsv`.

- Invariants: 19 for the CRISPR-Cas example and 18 for the RNA polymerase example, 37 in all. These are the 37 that
  were missing in the stopped check `ex-rM2`.
- Recorded values: all 55 are `new` (31 and 24), since no `recorded/rM.2.tsv` existed. They are copied from each
  example's `recorded_new.tsv` into `examples/<example>/recorded/rM.2.tsv` in the same commit.
- Arrays: `rM.2` holds none, so the CRISPR example's `[requires l1_crispr]` checks are `n/a`, as on `rN.3`. They wait
  for the release of card C-33.

**Against `rN.3`** (`recorded/rN.3.tsv`). `rM.2` adds the prokaryote stages P3 to P6 to `rN.3`'s S1, P1 and P2, and
the phage stage S1 is the same. So every virus count is unchanged and only the counts outside viruses grow:

- RNA polymerase:
  - Unchanged: 1,931 proteins in phages, T7 gene 1 at 367.7 bits, T7's own flank group (162 loci) and its kind
    `virus`, `Phage_T7_1_1` as the top neighbour.
  - Grown: proteins in prokaryotes 26 to 108, hit proteins 1,957 to 2,036, genes 2,581 to 2,663, flank groups 432
    to 505.
  - One context type now holds both phage and prokaryote loci (`rnap_flank_types_mixed_kind` 1, against 0).
- CRISPR-Cas:
  - Unchanged: 46 phage Cas7f loci, all with Cas1 and Cas3, and no TniQ locus in a phage.
  - Grown: TniQ loci 635 to 2,882 and Cas7f loci 111 to 502, both outside phages. TniQ loci beside Cas7f grew 20 to
    60, on 58 contigs of 55 genomes; 22 of them have TnsA, TnsB and TnsC, against 6.
  - The invariants that no TniQ-Cas7f locus carries Cas1, Cas2 and Cas3, and that at least half share Cas5f at +1
    (0.867), hold on the larger release.

**Cost.** CRISPR-Cas 4,870.4 CPU-seconds and RNA polymerase 3,082.0 (81.2 and 51.4 CPU-minutes, *derived* by the
script from the units' usage lines). That is 3.6 and 2.8 times `rN.3`'s (1,355.1 and 1,114.6 the same night, run
`w36-ex-rN3-new`), against 4.3 times the proteins (*derived*: 45,815,843 / 10,620,684).

- The largest step was the RNA polymerase search: 2,456.2 CPU-seconds (40.9 CPU-minutes), then the Cas7f search at
  33.8 and the TniQ search at 29.4 (*measured*, `CPUUsageNSec`).
- No step came near the 240 checked before the launch. The bound was each example's total on `rN.3` times 4.31, the
  protein ratio: at most 98.5 CPU-minutes per step (*derived*: 22.85 × 4.31).
- The RNA polymerase search estimate was 48.5 CPU-minutes (*derived*: 11.25 × 4.31); it measured 16% under that.
- The TniQ search unit's `MemoryPeak` was 17,179,869,184 bytes, exactly its 16 GiB `MEM` cap. It ended with exit 0:
  a cgroup's peak counts the page cache of the files it reads (here the 25.8 GB pack), and the search's own resident
  memory was not separated from it here. The other searches peaked at 15.3 and 15.7 GB. A card should say whether
  the search's own memory nears the cap on a 50 Gbp release.
- Wall time was 88 minutes (22:38:09 to 00:06:15), beside two runs on `rN.3` that ended at 23:10.

## Per example

| Example | run.sh exit | check exit | ok | FAIL | differs | new | missing | unexpected | CPU-s (*derived*) |
|---|---|---|---|---|---|---|---|---|---|
| `crispr_cas_neighbourhoods` | 0 | 0 | 19 | 0 | 0 | 31 | 0 | 0 | 4870.4 |
| `rnap` | 0 | 0 | 18 | 0 | 0 | 24 | 0 | 0 | 3082.0 |

CPU-s is the sum of the units' `CPUUsageNSec` (`pod_run: ... usage` lines) and of `/usr/bin/time -v` user
plus system time of the scopes, from the example's logs.

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
| `crispr_cas_neighbourhoods` | `tniq_cas7f_share_cas5f_at_plus1` | invariant | >= 0.5 | 0.867 | ok |
| `crispr_cas_neighbourhoods` | `tniq_pfam_proteins_not_hit` | invariant | == 0 | 0 | ok |
| `crispr_cas_neighbourhoods` | `cas7f_pfam_proteins_not_hit` | invariant | == 0 | 0 | ok |
| `crispr_cas_neighbourhoods` | `max_bits_diff_search_vs_pfam` | invariant | <= 0.1 | 0.0 | ok |
| `crispr_cas_neighbourhoods` | `oracle_tniq_only_search` | invariant | == 0 | 0 | ok |
| `crispr_cas_neighbourhoods` | `oracle_tniq_only_hmmsearch` | invariant | == 0 | 0 | ok |
| `crispr_cas_neighbourhoods` | `oracle_tniq_max_bits_diff` | invariant | <= 0.1 | 0.00 | ok |
| `crispr_cas_neighbourhoods` | `oracle_cas7f_only_search` | invariant | == 0 | 0 | ok |
| `crispr_cas_neighbourhoods` | `oracle_cas7f_only_hmmsearch` | invariant | == 0 | 0 | ok |
| `crispr_cas_neighbourhoods` | `oracle_cas7f_max_bits_diff` | invariant | <= 0.1 | 0.00 | ok |
| `crispr_cas_neighbourhoods` | `tniq_hit_proteins` | recorded | - | 2655 | new |
| `crispr_cas_neighbourhoods` | `cas7f_hit_proteins` | recorded | - | 504 | new |
| `crispr_cas_neighbourhoods` | `tniq_hits_not_in_pfam` | recorded | - | 60 | new |
| `crispr_cas_neighbourhoods` | `cas7f_hits_not_in_pfam` | recorded | - | 0 | new |
| `crispr_cas_neighbourhoods` | `tniq_loci_virus` | recorded | - | 0 | new |
| `crispr_cas_neighbourhoods` | `tniq_loci_not_virus` | recorded | - | 2882 | new |
| `crispr_cas_neighbourhoods` | `cas7f_loci_virus` | recorded | - | 46 | new |
| `crispr_cas_neighbourhoods` | `cas7f_loci_not_virus` | recorded | - | 502 | new |
| `crispr_cas_neighbourhoods` | `tniq_cas7f_loci` | recorded | - | 60 | new |
| `crispr_cas_neighbourhoods` | `tniq_cas7f_loci_virus` | recorded | - | 0 | new |
| `crispr_cas_neighbourhoods` | `tniq_cas7f_contigs` | recorded | - | 58 | new |
| `crispr_cas_neighbourhoods` | `tniq_cas7f_genomes` | recorded | - | 55 | new |
| `crispr_cas_neighbourhoods` | `tniq_cas7f_loci_with_tnsA_tnsB_tnsC` | recorded | - | 22 | new |
| `crispr_cas_neighbourhoods` | `tniq_cas7f_loci_with_cas8f_csy1` | recorded | - | 0 | new |
| `crispr_cas_neighbourhoods` | `tniq_minimal_non_if_loci` | recorded | - | 5 | new |
| `crispr_cas_neighbourhoods` | `tniq_minimal_ib_loci` | recorded | - | 0 | new |
| `crispr_cas_neighbourhoods` | `tniq_minimal_ic_loci` | recorded | - | 4 | new |
| `crispr_cas_neighbourhoods` | `tniq_other_cas_loci` | recorded | - | 31 | new |
| `crispr_cas_neighbourhoods` | `cas7f_not_virus_loci_without_cas1_cas3` | recorded | - | 80 | new |
| `crispr_cas_neighbourhoods` | `cas7f_not_virus_loci_without_cas1_cas3_with_tniq` | recorded | - | 59 | new |
| `crispr_cas_neighbourhoods` | `cas7f_virus_loci_with_cas1_and_cas3` | recorded | - | 46 | new |
| `crispr_cas_neighbourhoods` | `tniq_cas7f_loci_near_IMPDH` | recorded | - | 3 | new |
| `crispr_cas_neighbourhoods` | `tniq_cas7f_loci_near_4HBT` | recorded | - | 0 | new |
| `crispr_cas_neighbourhoods` | `tniq_cas7f_loci_near_SIS` | recorded | - | 0 | new |
| `crispr_cas_neighbourhoods` | `tniq_cas7f_loci_near_Mg_chelatase` | recorded | - | 1 | new |
| `crispr_cas_neighbourhoods` | `bp_tniq_loci` | recorded | - | 2817 | new |
| `crispr_cas_neighbourhoods` | `bp_cas7f_loci` | recorded | - | 548 | new |
| `crispr_cas_neighbourhoods` | `bp_tniq_cas7f_loci` | recorded | - | 61 | new |
| `crispr_cas_neighbourhoods` | `bp_cas7f_tniq_loci` | recorded | - | 59 | new |
| `crispr_cas_neighbourhoods` | `bp_tniq_cas7f_loci_with_tnsA_tnsB_tnsC` | recorded | - | 42 | new |
| `crispr_cas_neighbourhoods` | `bp_tniq_cas7f_loci_with_tnsB` | recorded | - | 51 | new |
| `rnap` | `t7_control_bits` | invariant | >= 30 | 367.7 | ok |
| `rnap` | `t7_control_architecture` | invariant | == RPOL_N+RNA_pol | RPOL_N+RNA_pol | ok |
| `rnap` | `t7_control_order` | invariant | == ocr Phage_T7_GP04 Gp0_5 Phage_T7_0_6 Phage_T7_kinase ^ Phage_T7_1_1 T7-like_gp12 DNA_ligase_A_M+DNA_ligase_C DUF5464 BC10 | ocr Phage_T7_GP04 Gp0_5 Phage_T7_0_6 Phage_T7_kinase ^ Phage_T7_1_1 T7-like_gp12 DNA_ligase_A_M+DNA_ligase_C DUF5464 BC10 | ok |
| `rnap` | `t7_window_is_the_contig_order` | invariant | == true | true | ok |
| `rnap` | `t7_positive_recovered` | invariant | == true | true | ok |
| `rnap` | `coord_t7_kinase_offset` | invariant | == -1 | -1 | ok |
| `rnap` | `coord_t7_ligase_offset` | invariant | == 3 | 3 | ok |
| `rnap` | `rnap_negative_control_hits` | invariant | == 0 | 0 | ok |
| `rnap` | `rnap_loci_virus` | invariant | >= 1 | 2555 | ok |
| `rnap` | `rnap_loci_not_virus` | invariant | >= 1 | 108 | ok |
| `rnap` | `rnap_loci_missing_windows` | invariant | == 0 | 0 | ok |
| `rnap` | `rnap_hits_below_cutoff` | invariant | == 0 | 0 | ok |
| `rnap` | `rnap_pfam_proteins_at_cutoff_not_hit` | invariant | == 0 | 0 | ok |
| `rnap` | `max_bits_diff_search_vs_pfam` | invariant | <= 0.1 | 0.0 | ok |
| `rnap` | `oracle_rnap_only_search` | invariant | == 0 | 0 | ok |
| `rnap` | `oracle_rnap_only_hmmsearch` | invariant | == 0 | 0 | ok |
| `rnap` | `oracle_rnap_max_bits_diff` | invariant | <= 0.1 | 0.00 | ok |
| `rnap` | `rnap_flank_loci_missing` | invariant | == 0 | 0 | ok |
| `rnap` | `rnap_hmm_sha256` | recorded | - | 32e6f0f0b1ef15d23d8ead23095f0b2a5a3323551b899bba99218d75f830e87e | new |
| `rnap` | `rnap_hit_proteins` | recorded | - | 2036 | new |
| `rnap` | `rnap_hit_genes` | recorded | - | 2663 | new |
| `rnap` | `rnap_hit_contigs` | recorded | - | 2602 | new |
| `rnap` | `rnap_hit_genomes` | recorded | - | 2601 | new |
| `rnap` | `rnap_proteins_virus` | recorded | - | 1931 | new |
| `rnap` | `rnap_proteins_not_virus` | recorded | - | 108 | new |
| `rnap` | `rnap_genes_with_architecture` | recorded | - | 2651 | new |
| `rnap` | `rnap_hits_not_in_pfam` | recorded | - | 11 | new |
| `rnap` | `rnap_pfam_proteins` | recorded | - | 2038 | new |
| `rnap` | `rnap_pfam_proteins_below_cutoff` | recorded | - | 13 | new |
| `rnap` | `rnap_summary_loci` | recorded | - | 2663 | new |
| `rnap` | `rnap_summary_contigs` | recorded | - | 2602 | new |
| `rnap` | `rnap_summary_censored_loci` | recorded | - | 466 | new |
| `rnap` | `rnap_summary_neighbours` | recorded | - | 25284 | new |
| `rnap` | `rnap_top_pfam_neighbour` | recorded | - | Phage_T7_1_1 | new |
| `rnap` | `rnap_top_pfam_neighbour_loci` | recorded | - | 1035 | new |
| `rnap` | `rnap_flank_groups` | recorded | - | 505 | new |
| `rnap` | `rnap_flank_leaves` | recorded | - | 1639 | new |
| `rnap` | `rnap_flank_context_types` | recorded | - | 682 | new |
| `rnap` | `rnap_flank_types_mixed_kind` | recorded | - | 1 | new |
| `rnap` | `t7_flank_group_loci` | recorded | - | 162 | new |
| `rnap` | `t7_flank_group_types` | recorded | - | 11 | new |
| `rnap` | `t7_flank_type_kinds` | recorded | - | virus | new |
