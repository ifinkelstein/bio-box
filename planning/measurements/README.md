# Measurement and Build Records (`planning/measurements/`)

Topical index of all **108** empirical measurement notes, stage build reports, and release records (*measured* count of `.md` files in `planning/measurements/`, 2026-11-23). Supporting TSV/JSON/SVG artifacts sit in matching date-prefixed subdirectories alongside their `.md` note.

See [`../../DOCS.md`](../../DOCS.md) for the repository-wide documentation map and [`../DECISIONS.md`](../DECISIONS.md) for standing decisions.

## Releases & Worked-Example Release Gates (8)

End-to-end builds, verifications, search packs, worked-example invariant checks, per-release source-diversity figures ([`release_diversity.json`](release_diversity.json) rendered by `scripts/plot_release_diversity.py`), and durable backups of published releases.

| Date / File | Summary Heading |
|---|---|
| [`2026-11-24_release_r26.11.3.md`](2026-11-24_release_r26.11.3.md) | Release `r26.11.3` over the thirteen stages (wilkcomp02, 2026-11-23 to 2026-11-25) + [`source_diversity_r26.11.3.png`](2026-11-24_release_r26.11.3/source_diversity_r26.11.3.png) |
| [`2026-11-21_release_r26.11.2.md`](2026-11-21_release_r26.11.2.md) | Release `r26.11.2` over the eleven stages (wilkcomp02, 2026-11-20 to 2026-11-22) |
| [`2026-11-19_release_train.md`](2026-11-19_release_train.md) | Card C-46: Release train rehearsal and real-root dry-run (wilkcomp02, 2026-11-19) |
| [`2026-11-16_release_r26.11.1.md`](2026-11-16_release_r26.11.1.md) | Release `r26.11.1` over the nine stages (wilkcomp02, 2026-11-16 to 2026-11-19) |
| [`2026-11-14_release_examples_rM2.md`](2026-11-14_release_examples_rM2.md) | Worked examples on `rM.2` (release check, AGENTS 5.7) |
| [`2026-11-13_release_examples_rN3.md`](2026-11-13_release_examples_rN3.md) | Worked examples on `rN.3` (release check, AGENTS 5.7) |
| [`2026-11-11_build_r50.md`](2026-11-11_build_r50.md) | The current-code release over the 10 Gbp stages, and the 50 Gbp build (wilkcomp02, 2026-11-11) |
| [`2026-11-10_import_ladder/rL3.md`](2026-11-10_import_ladder/rL3.md) | Release `rL.3` finished and synced (wilkcomp02, 2026-11-10) |

## Built Stages (`H1`, `PL1`, `C5`, `V1`, `G1`, `SPX1`) (8)

Stage build records, runtime/rate measurements, `stage_check.py` audits, and `/stor/work` + `tower` backups.

| Date / File | Summary Heading |
|---|---|
| [`2026-11-23_stage_SPX1.md`](2026-11-23_stage_SPX1.md) | Stage SPX1 (SPIRE v1 unbinned contigs, tranche 1: 304 aquatic and terrestrial samples) built into the data root |
| [`2026-11-23_stage_G1.md`](2026-11-23_stage_G1.md) | Stage G1 (GlobDB r226 non-GTDB archaeal species representatives) built into the data root |
| [`2026-11-21_a18_c_stage_confs.md`](2026-11-21_a18_c_stage_confs.md) | Stage driver configurations for `C1` to `C4` and `C6` (card A-18, DEC-62) |
| [`2026-11-20_stage_V1.md`](2026-11-20_stage_V1.md) | Stage V1 (UHGV v1.0 vOTU representatives) built into the data root (card C-48c part 1) |
| [`2026-11-20_stage_C5.md`](2026-11-20_stage_C5.md) | Stage C5 (GTDB R226 complete genomes, Enterobacterales to *Vibrio*) built into the data root (card C-47 part 1) |
| [`2026-11-15_stage_PL1.md`](2026-11-15_stage_PL1.md) | Stage PL1 (PLSDB plasmids) built into the data root (card C-26 part 2) |
| [`2026-11-15_stage_H1.md`](2026-11-15_stage_H1.md) | Stage H1, the huge phages, built into the data root (card C-24 as amended on 2026-11-15) |
| [`2026-11-10_build_stage_pod.md`](2026-11-10_build_stage_pod.md) | First real build with `build-stage` on the pod (card C-10, wilkcomp02, 2026-11-10) |

## Source Surveys, Crosswalks, Staging & Storage Relays (25)

External sequence/metadata inventories, JGI/SPIRE/GlobDB/PlasmidScope/virus staging, `DEC-65` terminal-repeat trims, deferred-accession registers, and `tower` NAS mirroring.

| Date / File | Summary Heading |
|---|---|
| [`2026-11-22_v2_gv1_dec65.md`](2026-11-22_v2_gv1_dec65.md) | Stage `V2` re-staged under DEC-65 and `V2.conf` / `GV1.conf` updated for build order after `G1` and `SPX1` (card C-53) |
| [`2026-11-21_pl2_pl3_dec63.md`](2026-11-21_pl2_pl3_dec63.md) | Stages PL2 and PL3 under DEC-63 and DEC-65: circular flags, terminal repeats, and final staging (card C-52) |
| [`2026-11-21_c51_deferred_list.md`](2026-11-21_c51_deferred_list.md) | The complete genomes that S2C phase 2's selection left out, by accession (card C-51, DEC-59) |
| [`2026-11-20_virus_v2_gv1_inputs.md`](2026-11-20_virus_v2_gv1_inputs.md) | IMG/VR and Inovirus as stage V2; giant-virus MAGs (GV1) measured but not staged (card C-48c, part 2) |
| [`2026-11-20_virus_gv1_inputs.md`](2026-11-20_virus_gv1_inputs.md) | GVMAGs V2 as stage GV1: inputs, overlap and cost (card C-48c, part 3) |
| [`2026-11-20_tower_relay.md`](2026-11-20_tower_relay.md) | The NAS `tower` and the relay through the laptop: capacities and rates (card C-49) |
| [`2026-11-20_spire_tranche_1.md`](2026-11-20_spire_tranche_1.md) | SPIRE tranche 1 (stage `SPX1`): 304 aquatic and terrestrial samples chosen, fetched and staged |
| [`2026-11-19_jgi_mags_iuig_dedup.md`](2026-11-19_jgi_mags_iuig_dedup.md) | Deduplication and Staging of JGI `GEMs`, `IMG-UMAG` (Prokaryote & Eukaryote MAGs), and `IUIG` after `GTDB` R226 and `GlobDB` r226 (Card C-48d) |
| [`2026-11-19_gold_img_crosswalk.md`](2026-11-19_gold_img_crosswalk.md) | JGI GOLD and IMG/M Environmental Metadata Crosswalk & DuckDB Sidecar (Card C-48b) |
| [`2026-11-18_spire_pilot.md`](2026-11-18_spire_pilot.md) | SPIRE 1,000-sample unbinned-contig pilot: habitat composition, depth-1 costs, NR50 novelty and on-the-hit annotation |
| [`2026-11-18_jgi_harvest_inventory.md`](2026-11-18_jgi_harvest_inventory.md) | 2026-11-18 JGI, NERSC and Zenodo Harvest Inventory (Card C-48a) |
| [`2026-11-18_globdb_inputs.md`](2026-11-18_globdb_inputs.md) | GlobDB release 226 species absent from GTDB: metadata, stage proposal, `G1` inputs, and full-set cost (card C-41) |
| [`2026-11-16_virus_votu_inputs.md`](2026-11-16_virus_votu_inputs.md) | Viruses from metagenomes as vOTU representatives: inputs of stage V1 and what they share (card C-39) |
| [`2026-11-16_run_crosswalk.md`](2026-11-16_run_crosswalk.md) | One assembly per sequencing run: SPIRE, MGnify, INSDC (NCBI WGS) and what IMG/M would add |
| [`2026-11-16_proposed_sources.md`](2026-11-16_proposed_sources.md) | Proposed source entries for cards C-39, C-40 and C-43 |
| [`2026-11-16_pl2_inputs.md`](2026-11-16_pl2_inputs.md) | Stage PL2 (plasmids from metagenomes and PlasmidScope): inputs, overlap and build estimates (card C-40) |
| [`2026-11-15_s2c_phase2_inputs.md`](2026-11-15_s2c_phase2_inputs.md) | S2C phase 2: the inputs of stages C1 to C6, and what they share (card C-36) |
| [`2026-11-15_plsdb_inputs.md`](2026-11-15_plsdb_inputs.md) | Stage PL1 (PLSDB plasmids): its inputs, what it adds, and how its build will end (card C-26, part 1) |
| [`2026-11-13_s2c_skani_selection.md`](2026-11-13_s2c_skani_selection.md) | S2C phase 2: skani run and the ANI-diversity selection (card C-27 step 4; 2026-11-13) |
| [`2026-11-13_s2c_skani_fetch.md`](2026-11-13_s2c_skani_fetch.md) | S2C phase 2: skani pinned, ANI selection code, the genome fetch (card C-27 steps 1 to 3; 2026-11-12/13) |
| [`2026-11-13_s2c_phase2_fetch.md`](2026-11-13_s2c_phase2_fetch.md) | S2C phase 2: fetch of the rest of the remainder (card C-25 step 2; 2026-11-13) |
| [`2026-11-12_s2c_selection.md`](2026-11-12_s2c_selection.md) | S2C selection: the ≤20-per-species complete genomes (card C-25, step 1; wilkcomp02, 2026-11-12) |
| [`2026-11-12_s2c_ani.md`](2026-11-12_s2c_ani.md) | S2C ANI diversity: does GTDB R226 publish the values? (card C-25, step 1b; 2026-11-12) |
| [`2026-11-11_s2_inputs.md`](2026-11-11_s2_inputs.md) | S2 inputs: HPGC huge phages and PLSDB plasmids (wilkcomp02, 2026-11-11) |
| [`2026-11-11_r50_inputs.md`](2026-11-11_r50_inputs.md) | Inputs for the 50 Gbp build (wilkcomp02, 2026-11-11) |

## Clustering, Sequence Store, Postings & Scaling Refits (25)

Empirical stress tests, memory/row-group bounds, split MMseqs2 clustering, representative Pfam evaluation, and rate-model calibrations.

| Date / File | Summary Heading |
|---|---|
| [`2026-11-22_b66_seqstore_and_crispr_contigs.md`](2026-11-22_b66_seqstore_and_crispr_contigs.md) | Card B-66: `O(1)`-stage `SeqStore` point lookups and per-stored-contig CRISPR CPU-minute term |
| [`2026-11-22_b65_protein_seqs_bounds.md`](2026-11-22_b65_protein_seqs_bounds.md) | Card B-65: Bounded Row-Group-Pruning `export.protein_seqs` and `hits_faa` at Multi-Stage Scale |
| [`2026-11-21_ingest_terms.md`](2026-11-21_ingest_terms.md) | The ingest estimate prices gzip input and bases read but not stored (card B-61) |
| [`2026-11-21_c5_ingest_cost.md`](2026-11-21_c5_ingest_cost.md) | What stage C5's ingest spent its CPU-minutes on (card B-58, amended 2026-11-21) |
| [`2026-11-21_b64_cluster_factors.md`](2026-11-21_b64_cluster_factors.md) | Card B-64: Release Clustering Calibration Factors Across `rM.2`, `r26.11.1`, and `r26.11.2` |
| [`2026-11-21_b62_postings_bounds.md`](2026-11-21_b62_postings_bounds.md) | Card B-62: Bounded memory and open files in `biobox.l3.postings` |
| [`2026-11-20_peak_helper_evidence.md`](2026-11-20_peak_helper_evidence.md) | Who kills the peak helper: card A-10's diagnosis and the evidence card A-10a keeps (2026-11-20) |
| [`2026-11-20_c5_pfam_by_representatives.md`](2026-11-20_c5_pfam_by_representatives.md) | Pfam by NR90 representatives on stage C5: DEC-35 (6) stops on missed pairs (card C-47 part 2) |
| [`2026-11-20_audit_of_1117_1120.md`](2026-11-20_audit_of_1117_1120.md) | Audit of the 2026-11-17 to 20 campaign: what was claimed, and what was checked again |
| [`2026-11-19_cluster_survivor_estimates.md`](2026-11-19_cluster_survivor_estimates.md) | Refitting intermediate cascade (`input_step1` / `input_step2`) and `kmermatcher` estimates (wilkcomp02, 2026-11-19) |
| [`2026-11-16_ingest_estimate.md`](2026-11-16_ingest_estimate.md) | The ingest estimate counts what a stage stores (card B-52, DEC-47 (2)) |
| [`2026-11-15_label_files_old_stages.md`](2026-11-15_label_files_old_stages.md) | Label files for the stages ingested before card B-27 (card C-34) |
| [`2026-11-15_cluster_reuse.md`](2026-11-15_cluster_reuse.md) | A release takes an earlier release's clustering (card B-49, DEC-45 (2)) |
| [`2026-11-12_split_clustering_rN1b.md`](2026-11-12_split_clustering_rN1b.md) | Split clustering reproduces `rN.1b` (card B-42 step 1) |
| [`2026-11-11_representatives_only.md`](2026-11-11_representatives_only.md) | Scanning cluster representatives only: saving and inference error (card 0-09b, wilkcomp02, 2026-11-11) |
| [`2026-11-11_mmseqs_profile_vs_hmmsearch.md`](2026-11-11_mmseqs_profile_vs_hmmsearch.md) | MMseqs2 profile search against the producer's `hmmsearch` at Pfam GA (card 0-09a, wilkcomp02, 2026-11-10/11) |
| [`2026-11-10_protein_stars.md`](2026-11-10_protein_stars.md) | Stop codons in `protein_seq` (card W-11) |
| [`2026-11-10_postings_memory.md`](2026-11-10_postings_memory.md) | Postings builder memory (card B-29, wilkcomp02, 2026-11-10) |
| [`2026-11-10_postings_bytes.md`](2026-11-10_postings_bytes.md) | `token_postings` bytes per column (card B-26) |
| [`2026-11-10_phase_b_exit.md`](2026-11-10_phase_b_exit.md) | Phase-B exit measurement, 2026-11-10 (card B-18) |
| [`2026-11-10_label_and_membership_lookups.md`](2026-11-10_label_and_membership_lookups.md) | Contig labels and cluster membership by their lookup keys (card B-27), 2026-11-10 |
| [`2026-11-10_cluster_memory.md`](2026-11-10_cluster_memory.md) | Clustering memory before and after card B-24 (wilkcomp02, 2026-11-10) |
| [`2026-11-10_cluster_id_churn.md`](2026-11-10_cluster_id_churn.md) | Cluster id churn between releases when each release clusters from scratch (wilkcomp02, 2026-11-10) |
| [`2026-11-09_design_stress_test.md`](2026-11-09_design_stress_test.md) | Design stress test, 2026-11-09 |
| [`2026-10-30_k13_bakeoff.md`](2026-10-30_k13_bakeoff.md) | K-13 bake-off, 2026-10-30: Parquet `seq_chunks` vs bgzip + faidx, scored by LY-06 |

## Ad Hoc Search, Search Packs & Locus Flank Comparison (12)

Custom-HMM search trials, memory profiles, sharded pack builds, positive controls, and `nbr flanks` synteny grouping.

| Date / File | Summary Heading |
|---|---|
| [`2026-11-14_pack_build_in_parts.md`](2026-11-14_pack_build_in_parts.md) | The search pack built in parts: `rN.3` reproduced, `rM.2` built (card W-35, 2026-11-14) |
| [`2026-11-14_flanks_figure_defaults.md`](2026-11-14_flanks_figure_defaults.md) | `nbr flanks` figure defaults of DEC-41 (card W-32) |
| [`2026-11-13_flanks_context_types.md`](2026-11-13_flanks_context_types.md) | `nbr flanks` on mobile elements: context types and figure colours (card W-31) |
| [`2026-11-12_positive_lookup_rN3.md`](2026-11-12_positive_lookup_rN3.md) | Resolving a positive control by accession on `rN.3` (card W-15), 2026-11-12 |
| [`2026-11-12_flanks_rN3.md`](2026-11-12_flanks_rN3.md) | `nbr flanks` on rN.3: the RNA polymerase loci grouped by their flanks (card W-25) |
| [`2026-11-11_trial_search_rN1_second.md`](2026-11-11_trial_search_rN1_second.md) | Second trial: a custom-HMM search with neighbourhoods on `rN.1`, end to end (2026-11-11) |
| [`2026-11-11_search_memory_rN1.md`](2026-11-11_search_memory_rN1.md) | A search on `rN.1` after card W-19: memory, time and the first `search status` (2026-11-11) |
| [`2026-11-11_pack_build_rN1.md`](2026-11-11_pack_build_rN1.md) | The search pack of `rN.1`, built by the streaming builder (card W-16, 2026-11-10 to 11) |
| [`2026-11-10_trial_search_rN1.md`](2026-11-10_trial_search_rN1.md) | Trial: a custom-HMM search with neighbourhoods on `rN.1` (2026-11-10) |
| [`2026-11-10_token_lookup.md`](2026-11-10_token_lookup.md) | `token_of` by stable key, and catalog connections (card B-34), 2026-11-10 |
| [`2026-11-10_search_pack_size.md`](2026-11-10_search_pack_size.md) | Search pack size, and why the pack is kept off shared storage (card W-03, 2026-11-10) |
| [`2026-11-10_neighbour_pruning.md`](2026-11-10_neighbour_pruning.md) | Neighbour and lookup macros that prune (card B-31), 2026-11-10 |

## Gene Calling, CRISPR Arrays & Non-Coding RNA (`Infernal`) (19)

Prodigal/pyrodigal batching and meta-mode scaling, default 3-repeat and filtered 2-repeat near-Cas CRISPR producers, and curated `prok-2` Rfam intergenic scanning.

| Date / File | Summary Heading |
|---|---|
| [`2026-11-21_crispr_estimate_refit.md`](2026-11-21_crispr_estimate_refit.md) | Can a term bring both CRISPR estimates within 15% of the eleven built stages? (card B-60) |
| [`2026-11-21_b63_meta_genecall.md`](2026-11-21_b63_meta_genecall.md) | Meta-mode gene-calling cost and batch sizing before `PL3` (card B-63) |
| [`2026-11-17_crispr_genecall_estimates.md`](2026-11-17_crispr_genecall_estimates.md) | Refitted estimates for CRISPR on circular plasmids and gene calling on circular viruses (card B-53) |
| [`2026-11-15_genecall_batches.md`](2026-11-15_genecall_batches.md) | Gene calling in batch steps under the CPU budget (card B-48) |
| [`2026-11-15_crispr_near_cas_phages.md`](2026-11-15_crispr_near_cas_phages.md) | The filtered CRISPR producer on phages: its chance rate on S1 (card E-16) |
| [`2026-11-15_crispr_arrays_on_built_stages.md`](2026-11-15_crispr_arrays_on_built_stages.md) | CRISPR arrays of both producers on the seven built stages of the data root (card C-33 part 1, DEC-40) |
| [`2026-11-14_crispr_two_repeat_producer.md`](2026-11-14_crispr_two_repeat_producer.md) | The producer of two-repeat CRISPR arrays beside Cas genes (card E-15, DEC-40) |
| [`2026-11-14_crispr_min_repeats.md`](2026-11-14_crispr_min_repeats.md) | CRISPR arrays at a minimum of two repeats against three (card E-13) |
| [`2026-11-14_arrays_in_neighbourhoods.md`](2026-11-14_arrays_in_neighbourhoods.md) | CRISPR arrays beside a locus (card E-14) |
| [`2026-11-13_crispr_arrays.md`](2026-11-13_crispr_arrays.md) | CRISPR arrays by MinCED on stages P1, P2 and S1 (card E-12) |
| [`2026-11-12_rfam_sample.md`](2026-11-12_rfam_sample.md) | Which long-tail Rfam models hit on a sample of the corpus, 2026-11-12 (card 0-08a) |
| [`2026-11-12_ncrna_set.md`](2026-11-12_ncrna_set.md) | The ncRNA set as one scan, tRNAs lost to filters, and short contigs, 2026-11-12 (card 0-08b) |
| [`2026-11-12_ncrna_prok2.md`](2026-11-12_ncrna_prok2.md) | The ncRNA producer with `prok-2` on the 8 genomes of 0-09c, 2026-11-12 (card E-09) |
| [`2026-11-12_context_free_genecall.md`](2026-11-12_context_free_genecall.md) | Context-free gene calling (policy p3) on stage P1, 2026-11-12 (card C-21f) |
| [`2026-11-11_ncrna_producer.md`](2026-11-11_ncrna_producer.md) | The ncRNA producer of card E-07 on the 8 genomes of 0-09c, 2026-11-11 |
| [`2026-11-11_ncrna_clans.md`](2026-11-11_ncrna_clans.md) | Clan competition for ncRNA hits (card E-08) on the 8 genomes of 0-09c, 2026-11-11 |
| [`2026-11-11_infernal_intergenic.md`](2026-11-11_infernal_intergenic.md) | Infernal on intergenic regions against whole-genome scanning, 2026-11-11 (card 0-09c) |
| [`2026-11-10_ncrna_rate.md`](2026-11-10_ncrna_rate.md) | ncRNA scan rate with full Rfam 15.1, 2026-11-10 (card E-01) |
| [`2026-11-10_genecall_short_contigs.md`](2026-11-10_genecall_short_contigs.md) | Gene calling of short contigs, 2026-11-10 (card B-32) |

## Synthetic Worlds, Null Models, Images & Test Harness (10)

Constructed-truth synthetic genome worlds (`F-02`), S0 golden-set gates, container images, and CI diagnostics.

| Date / File | Summary Heading |
|---|---|
| [`2026-11-14_synthworld_seeds.md`](2026-11-14_synthworld_seeds.md) | synthworld `micro` over 21 seeds: why it holds its size and planted share only at the default seed |
| [`2026-11-13_synthworld_nulls.md`](2026-11-13_synthworld_nulls.md) | synthworld null worlds (§4.1) on `micro` |
| [`2026-11-13_synthworld_mobile.md`](2026-11-13_synthworld_mobile.md) | synthworld S6 (phages, prophages, CRISPR arrays) on `micro`: time, memory and realized values |
| [`2026-11-13_synthworld_fusions.md`](2026-11-13_synthworld_fusions.md) | synthworld S7 (gene fusions, splits, NOGD) on `micro`: time, memory and realized values |
| [`2026-11-13_synthworld_emit.md`](2026-11-13_synthworld_emit.md) | synthworld S11 (emit) on `micro`, and the emitted world through the normal builders |
| [`2026-11-13_synthworld_assembly.md`](2026-11-13_synthworld_assembly.md) | synthworld S8-S10 (call errors, codes 4 and 25, assembly, contamination, duplicates) on `micro`: time, memory and realized values |
| [`2026-11-13_image_faec1f4.md`](2026-11-13_image_faec1f4.md) | The images of `faec1f4` on the pod (card C-32) |
| [`2026-11-10_ci_flakes.md`](2026-11-10_ci_flakes.md) | Two search tests that failed intermittently on GitHub Actions, 2026-11-10 (card A-07) |
| [`2026-10-30_wilkcomp02.md`](2026-10-30_wilkcomp02.md) | wilkcomp02, 2026-10-30: first run from the image |
| [`2026-10-30_s0_gate.md`](2026-10-30_s0_gate.md) | S0 gate, first report (2026-10-30) |

## Other Notes (2)

| Date / File | Summary Heading |
|---|---|
| [`2026-11-15_b47b_copy_plan_H1.md`](2026-11-15_b47b_copy_plan_H1.md) | Copy measured Pfam rows instead of scanning again: the plan of stage H1 (card B-47b) |
| [`2026-11-14_partial_labels_rM2.md`](2026-11-14_partial_labels_rM2.md) | A release whose stages are only partly labelled: lookups by accession on `rM.2` (card W-36) |
