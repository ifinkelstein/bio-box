# The goal, the order, and the register of what is deferred

Started 2026-11-20 at the PI's word (DEC-59). This file is the one place that says what the corpus is meant to
hold in the end, why sources are taken in the order they are, and, for every selection made on the way, what was
left out and how to go back for it. **Nothing in this file is a decision to leave a source out for good.** Keep it
current: a card that selects, filters or dereplicates adds its row here in the same commit as its note.

Numbers carry the labels of AGENTS §6. Those marked *(card)* are copied from the card or note named beside them and
were not measured again for this file.

## 1. The goal

In the PI's words (2026-11-20): "the ultimate goal is a COMPLETE catalog across all useful seuqeced data", and "i do
want to go and pick up less diverse sources in full later".

So:

- **The end state is completeness**: every useful sequenced genome, MAG, plasmid, virus and assembled metagenome
  contig, with its neighbourhood, not a diverse sample of them.
- **Diversity first is an order of work, not a scope.** Sources are ranked by new protein families (NR50 clusters)
  per Gbp because storage and CPU are short now (STATUS open item 19). A source ranked low is taken later, in full.
- **Every selection is reversible on paper.** A step that takes representatives, dereplicates, filters by length or
  stops a habitat early records what it did not take, with identifiers sufficient to fetch and ingest the remainder
  later (accessions, sample ids, the archive member, the filter's parameters), in its measurement note and in §4
  below.
- **What "complete" covers (DEC-60, the PI, 2026-11-20).** For unbinned metagenome assemblies: every contig of
  5,000 bp or more; contigs under 5 kb are outside the corpus and are not deferred data. RNA viruses are in scope
  (by a design card, C-50). Protein-only catalogues are deferred ("can wait"). Raw reads are not taken (Logan's
  assemblies stand in for them, §4.3).
- **Where it lives** is the main constraint (STATUS open item 19). The pod's local disk holds the working root; the
  PI's NAS `tower` is offered as a place for tranches as the corpus grows (DEC-60 (3), card C-49).

## 2. The order of work, and why

Measured yield, new NR50 clusters per Gbp (*derived* from the totals in each note; the SPIRE figures come from
`linclust` on the pilot and the others from the release clustering, so compare them loosely):

| Source type | New NR50 clusters per Gbp | From |
|---|---|---|
| Huge phages (`H1`) | about 489,000 | 890,552 new clusters over 1,819,831,338 bp (`r26.11.1` against `rM.2`; card C-35) |
| SPIRE unbinned contigs, aquatic | 465,323 | card C-42, `novelty/novelty_by_habitat1.tsv` (*measured*, contigs of 1.5 kb and more) |
| SPIRE unbinned, anthropogenic | 425,040 | the same |
| SPIRE unbinned, terrestrial | 407,439 | the same |
| SPIRE unbinned, host-associated | 269,797 | the same |
| Plasmids from isolates (`PL1`) | about 34,000 | 200,698 new clusters over 5,924,552,875 bp (card C-35) |

Of the pilot's 41,282,757 new clusters, 9,466,693 (23%) have more than one member (*measured*, card C-42): part of
the unbinned novelty is fragments. Yield fell between the first and second half of the pilot's samples in every
large habitat (for example terrestrial 465,288 to 349,345 per Gbp), so each habitat has a point where the next
tranche stops paying; at that point the habitat is **deferred, not finished** (§4).

The order that follows (the orchestrator's, 2026-11-20, for the PI to change):

1. Everything already on disk (§3), with rolling releases.
2. Viral vOTU representatives in full (IMG/VR v4.0, on disk) at depth `called`.
3. Unbinned metagenome contigs of 5 kb and more (DEC-57), by habitat, aquatic and terrestrial first: SPIRE, then
   IMG/M's unrestricted metagenomes that SPIRE lacks, then MGnify's and INSDC's assemblies outside both.
4. Species representatives of cellular genomes: GTDB R226's and GlobDB's.
5. Then the less diverse remainder of every source, in full (§4).

## 3. What is in, and what is on disk and queued

| Source | State on 2026-11-20 | Stage |
|---|---|---|
| INPHARED phages (14Apr2025) | built | `S1` |
| GTDB R226, best complete genome per species (11,966 fetched) | built | `P1` to `P6` |
| Huge Phage Genome Collection | built | `H1` |
| PLSDB 2024_05_31_v2 | built | `PL1` |
| UHGV v1.0 vOTU representatives (168,536 genomes; 168,358 contigs stored, 5.18 Gbp, 7,364,696 genes) | built 2026-11-20; in no release yet | `V1` |
| GTDB R226 complete genomes, second selection (13,477 genomes, 56.18 Gbp) | `C5` built 2026-11-20 (2,029 genomes, 9.11 Gbp stored, fully scanned; in release `r26.11.2`); `C1` to `C4` and `C6` (46.79 Gbp) held since 2026-11-21 for the PI's decision on how their Pfam is made (STATUS open item 21) | `C1` to `C6` |
| Plasmids from metagenomes (Yu 2024, PlasmidScope, IMG/PR; 859,036 new contigs) | inputs staged | `PL2` |
| IMG/VR v4.1 high-confidence vOTU representatives (224,520) and Inovirus (5,964 records): 230,098 genomes, 261,626 contigs, 9.34 Gbp, of which 248,864 contigs and 8.36 Gbp are new (*measured*, card C-48c part 2) | inputs staged | `V2` |
| GVMAGs V2: 18,499 genomes, 281,804 contigs, 4.32 Gbp, of which 264,134 contigs and 3.89 Gbp are new after `V2` (*measured*, card C-48c part 3) | inputs staged | `GV1` |
| GlobDB r226 non-GTDB archaea (7,139 genomes) | inputs staged | `G1` |
| SPIRE v1 unbinned contigs, tranche 1: 304 aquatic and terrestrial samples, 2,829,368 contigs of 5 kb and more, 35.885 Gbp (*measured*, card C-48e part 1); to be built at depth `called` after `G1`, which owns the 12,052 contigs they share | inputs staged | `SPX1` |
| IUIG novel prokaryotes (about 4,235 representatives, *derived*, card C-48d) | archive on disk; not staged | `J1` |
| IMG-UMAG eukaryote bins (691 bins, 14.44 Gbp) | archive on disk; waits on a `"eukaryote"` genome kind | `E1`, `E2` |

## 4. The register of what is deferred

Each row is something a selection left out, or a source not yet fetched. "How to pick it up" names what a later
card needs. Rows are added, and a row is struck through with its stage's name when its remainder is in the corpus.

### 4.1 Left out by a selection inside a source we hold

| Source | Taken | Left for later, in full | Where the remainder is recorded | How to pick it up |
|---|---|---|---|---|
| GTDB R226 genomes (732,475 rows *(card C-48d)*) | 11,966 best complete genomes per species (`P1` to `P6`); 13,477 further complete genomes (`C1` to `C6`) | every other GTDB genome: the species representatives that are not complete genomes (plan stages S3 and S5; 143,607 representatives, 465.47 Gbp, are also inside GlobDB's archive *(card C-41)*) and all non-representative genomes (plan stage S6) | `planning/measurements/2026-11-11_r50_inputs.md`, `2026-11-15_s2c_phase2_inputs.md` (the selections and their `remainder.tsv`, `rest.tsv`); `sources.yaml` | fetch by assembly accession from NCBI datasets (cards C-25, C-27 give the form); needs item 19 answered: S6 is 111.5 M contigs by the plan's count |
| GTDB R226 complete genomes of stage `C5`'s 306 species (Enterobacterales to *Vibrio*) | `C5`: 2,032 accessions (2,029 genomes), up to 20 per species by ANI diversity (DEC-34), built 2026-11-20 | the complete genomes of the 52 species over the cap that the selection did not pick: **10,318** (*measured*, card C-51, 2026-11-21, equal to card C-47's *derived* 11,358 complete genomes in those species, less 988 in `C5` and 52 in phase 1; `planning/measurements/2026-11-20_stage_C5.md` §11); 7 of them were not served by NCBI and are marked `not_served = y` | **per accession**: `planning/measurements/2026-11-21_c51_deferred_list/left_out.tsv` (all six `C` stages together, 23,866 rows, 1,374,419 bytes, sha256 `d0f10fd5e7a0fdc6aa2b3513aa4caea8230ee459dc53d19f152ef4f2de4e55a1`; `C5`'s species are the rows whose species `species_check.tsv` gives stage `C5`; 10,318 rows); note `planning/measurements/2026-11-21_c51_deferred_list.md`; per species: `2026-11-13_s2c_skani_selection/ani_species.tsv` (agrees in all 243 species); the picks: `2026-11-15_s2c_phase2_inputs/stages/stages.tsv`. Copy on `tower` `deferred_lists/C/`: relay prepared by card C-51 (dry run only), the transfer is the orchestrator's | candidates not in the picks, by assembly accession from NCBI datasets as card C-27 fetched them; a further `C` stage with `contig_uid` dedup against `C5` |
| GTDB R226 complete genomes of stage `C1`'s 579 species (Archaea to *Anoxybacillus_A*), `C1` not built | `C1` inputs: 1,783 genomes, 19 species over the cap, up to 20 per species by ANI diversity (DEC-34) | the complete genomes of those 19 species that the selection did not pick: 1,018 (*measured*, card C-51) | the same list and note as the `C5` row; stage per species in `species_check.tsv` | as the `C5` row |
| GTDB R226 complete genomes of stage `C2`'s 526 species (*Geobacillus* to *Hungateiclostridium*), `C2` not built | `C2` inputs: 3,102 genomes, 67 species over the cap | 5,979 (*measured*, card C-51), of which 2 not served | the same | as the `C5` row |
| GTDB R226 complete genomes of stage `C3`'s 775 species (*Aristaeella* to *Bradyrhizobium diazoefficiens*), `C3` not built | `C3` inputs: 2,470 genomes, 34 species over the cap | 1,712 (*measured*, card C-51) | the same | as the `C5` row |
| GTDB R226 complete genomes of stage `C4`'s 458 species (*Bradyrhizobium elkanii* to *Edwardsiella piscicida*), `C4` not built | `C4` inputs: 1,972 genomes, 33 species over the cap | 2,116 (*measured*, card C-51) | the same | as the `C5` row |
| GTDB R226 complete genomes of stage `C6`'s 541 species (*Vibrio zhugei* to *Akkermansia*), `C6` not built | `C6` inputs: 2,118 genomes, 38 species over the cap | 2,723 (*measured*, card C-51), of which 17 not served | the same | as the `C5` row |
| Genomes NCBI did not serve for phase 2 | none | the 34 accessions of DEC-34 (3) (26 candidates, 8 others), not retried | `planning/measurements/2026-11-13_s2c_phase2_fetch/not_served.txt` (8 accessions) and the C-27 note (the rest). The 26 candidates are inside the per-accession list of the rows above, marked `not_served = y` (`C2` 2, `C5` 7, `C6` 17); the 8 failed picks are not in that list (species at or under the cap): by stage of their species `C2` 1, `C3` 2, `C5` 1, `C6` 1 and 3 in species with no stage (*measured*, card C-51, `2026-11-21_c51_deferred_list.md` §4) | fetch again by accession when NCBI serves them |
| UHGV v1.0 | one representative per vOTU (168,536; stage `V1`, built 2026-11-20) | the 705,458 non-representative genomes (*derived*: 873,994 − 168,536, card C-48c part 1) | `planning/measurements/2026-11-16_virus_votu_inputs.md` §3 (`uhgv_metadata.tsv` has every genome's vOTU) | fetch the provider's full genome FASTA; stage as a further `V` stage, `contig_uid` dedup against `V1` |
| IMG/VR v4.1 high-confidence (5,576,197 UViGs, 2,917,521 vOTUs, *measured*, equal to the published counts; card C-48c part 2) | `V2`: 224,520 vOTU representatives, one per vOTU that has a high-quality member (9.29 Gbp), by highest completeness, then lowest contamination, then longest (the provider names no representative) | 538,684 other members of those vOTUs (15.3 Gbp; 262,171 of them high-quality: strain-level copies) and 2,693,001 vOTUs without a high-quality member (4.81 M UViGs, 57.7 Gbp) | `/tmp/if2453/s2/inputs/V2/V2_left_for_later_hc.tsv.gz` (5,351,677 rows, sha256 `ace2ec98…5ad5b62`; **on `/tmp` only until mirrored**, made again by `c48c_v2_gv1_left.py`); `planning/measurements/2026-11-20_virus_v2_gv1_inputs.md` §7 | no fetch needed: the same member of `jgi_img_vr_v4.1_hc_20261118.zip`, by id; `--depth called`; the vOTU representatives of this remainder take the root past 7.86 M contigs, so card B-56's design (open item 16) comes first |
| IMG/VR v4.0 full, the putative part (8,606,551 distinct vOTU ids in its table, *measured*; `sources.yaml` says 9,169,185, *published*) | none | 10,101,426 putative UViGs (77.8 Gbp; 5,689,030 vOTUs that are putative only); at 10 kb and more, 2.15 M vOTUs and 45.0 Gbp (*derived*, the same note §6) | `/tmp/if2453/s2/inputs/V2/V2_left_for_later_v4.0_putative.tsv.gz` (sha256 `d607ba78…0c80751`; on `/tmp` only until mirrored) | no fetch needed: `jgi_img_vr_v4_full_20261118.zip`, by id; `--depth called` with a least contig length (geNomad did not call these high-confidence) |
| Inovirus (Roux 2019) | `V2`: all 5,964 records | none (364 repeat an IMG/VR representative or each other and are listed as second genomes, DEC-27) | `V2/repeated_records.tsv` in the note's directory | nothing to pick up |
| GVMAGs V2 (18,726 giant-virus genomes; no quality table in the archive) | `GV1`, staged whole after `V2` (the orchestrator's choice of 2026-11-20, DEC-59: no tier exists to select by): 18,499 genomes (227 of the 18,726 repeat an earlier genome's contig set and are listed once), 281,804 contigs, 4.32 Gbp, 3.89 Gbp of it new (*measured*, card C-48c part 3) | none: the provider's 8,508-genome `fna_files` subset holds no contig `GV1` lacks (8,066 genomes identical, 442 with 746 contigs removed). 17,664 contigs are owned by `V2`'s IMG/VR giant-virus representatives and 1,116 genomes have every contig held earlier; they stay genomes of `GV1` | `planning/measurements/2026-11-20_virus_gv1_inputs.md` §4, §6; its `repeated_records.tsv`, `skipped_members.tsv` | nothing to pick up |
| Plasmid sets of `PL2` | 859,036 new contigs | 15,135 circular gut contigs that are not representatives (PI question, card C-40); records skipped by the stage's rules | `planning/measurements/2026-11-16_pl2_inputs.md`; `/tmp/if2453/s2/inputs/PL2/skipped_records.tsv` (on `/tmp` only: copy into the note's directory before `PL2` is built) | a further plasmid stage from the same FASTA files |
| IUIG (143,796 genome rows *(card C-48d)*) | about 4,235 novel prokaryote representatives proposed (`J1`) | 120,654 prokaryote genomes that match a GTDB or GlobDB species (strain copies); 15,544 virus, 1,185 plasmid and 721 eukaryote rows | `planning/measurements/2026-11-19_jgi_mags_iuig_dedup.md` | no fetch needed: `jgi_iuig_genomes_20261118.zip` (524.8 GB) is on `/stor/work`; stream by member |
| IMG-UMAG prokaryote bins (18,912 metagenomes) | none yet; novel representatives proposed (`J2` to `J15`, 120 to 140 Gbp, *estimate*) | bins that match a held species | the same note | no fetch needed for the bins; their dereplication wants GlobDB's and GTDB's representatives on disk, or a decision to ingest without it |
| GEMs (52,515 MAGs, 45,599 OTUs *(card C-48d)*) | none as genomes (every OTU is already judged by GlobDB); the 24,313 VirSorter prophage coordinates as a sidecar | all 52,515 MAGs as genomes | the same note | no fetch needed: `jgi_gems_52515_fna_20261118.tar` is on `/stor/work` |
| SPIRE pilot (990 assemblies, 341.5 Gbp raw) | measurements only; nothing built into the root | the assemblies themselves (unlinked after hashing); 12 of its samples were fetched again and are in `SPX1` | `planning/measurements/2026-11-18_spire_pilot.md`, `select/pilot_selection.tsv` | fetch again by sample id (card C-48e) |
| SPIRE v1 unbinned contigs, tranche 1 (card C-48e part 1) | 304 aquatic and terrestrial samples (19 strata, 189 studies), staged as `SPX1`: 2,829,368 contigs of 5 kb and more, 35.885 Gbp (*measured*); not built yet | the other 17,968 aquatic and terrestrial samples with a published size (8.50 Tbp assembled; about 83.3 M contigs and 958.6 Gbp at 5 kb and more, *estimate*) and 361 without one; every sample of the other habitats (host-associated 63,909, anthropogenic 9,662, unlabelled 3,846, air 17) | `planning/measurements/2026-11-20_spire_tranche_1/out/remainder_aquatic_terrestrial.tsv.gz` (in the order a next tranche takes them), `out/other_habitats.tsv`, `select/by_stratum.tsv` | `c48e_select.py run --exclude <earlier selection.tsv>`, then `c48e_fetch.py` by sample id and `c48e_inputs.py`; a second tranche waits on open items 16 and 19 |

### 4.2 Left out by a filter

| Filter | Decision | What it leaves out | How to pick it up |
|---|---|---|---|
| Least contig length 5,000 bp on unbinned metagenome stages | DEC-57, DEC-60 (1) | contigs under 5 kb: 83.9% of assembled bases and 99.05% of contigs in the pilot (*derived*: 1 − 55.050 / 341.501 Gbp; 1 − 4.29 M / 453.49 M; **corrected 2026-11-20**: this line and DEC-57 first said 4.81 M contigs of 5 kb and more, a figure from the earlier orchestrator's summary; the pilot's note and its per-sample table say 4.29 M, 4,292,779). In tranche 1 as staged: 99.13% of records and 85.23% of bases (*measured*, card C-48e part 1) | **Not deferred: outside the corpus** (DEC-60 (1): "in full means contigs 5 kb and longer"). Listed here so that the share is on record; each tranche's note still says how many contigs and bases the filter dropped |
| A habitat stopped when its yield falls | card C-48e's stopping rule | the habitat's remaining samples | the tranche note lists the samples not taken; card C-43's crosswalk has every sample |
| Depth `called` in place of `annotated` | DEC-54 | Pfam domains outside the windows that rounds annotate | `bio-box annotate round`, or a later stage step: no data is lost |
| Depth `stored` for eukaryote bins | card C-48d | genes: no intron-aware caller is in the image | a later card adds a caller; the sequences are held |

### 4.3 Sources not yet fetched

| Source | Size | Why not yet | Needs |
|---|---|---|---|
| SPIRE v1 unbinned contigs | 96,067 samples; about 23.35 Tbp assembled, 3.76 Tbp at 5 kb and more (*derived*, card C-42) | first tranche next (card C-48e part 1) | item 19 beyond two or three tranches |
| IMG/M unrestricted metagenomes not in SPIRE | 36,580 metagenomes, 17.307 Tbp *(card C-48b)* | after SPIRE's first tranches | the PI's JGI session for each fetch; item 19 |
| MGnify assemblies not in SPIRE | 37,628 assemblies *(card C-48e)* | after the two above | item 19 |
| INSDC WGS metagenome assemblies outside all three | 2,306 to 4,948 *(card C-48e)* | after them | |
| GlobDB r226 genome archive | 251,371,669,660 bytes; 162,653 non-GTDB species (410.33 Gbp, 48,999,006 contigs) and 143,607 GTDB representatives (465.47 Gbp) *(card C-41)* | it was streamed once and only `G1` kept; a second fetch would stay on `/stor` | item 19 |
| Strain-level sets (AllTheBacteria's 2.44 M isolates; the MAG sets behind GlobDB's species, such as mOTUs' 2.83 M MAGs) *(card C-45)* | not sized here | least new diversity per Gbp | DEC-59 reads them as "in, later"; item 19 |
| Logan (assemblies of all of SRA) | not sized here | taken on the hit only (card C-44) | a decision that "complete" includes it in bulk |
| Global Soil Virus Atlas | not fetched | the provider's host refused the pod five times *(STATUS)* | another route to the files |
| JGI sets on tape (`SMC_2025` 750 MB table and 6.11 TB of GFFs; `IMG_NR` 58.14 GB) *(card C-48e)* | | JGI's tape restore was down on 2026-11-18 | a check that it is back |
| RNA viruses (LucaProt's 161,979 species, Serratus' RdRps, RVMT) *(card C-45)* | not sized here | **in scope (DEC-60 (2))**; everything built assumes DNA genomes | card C-50: a gene caller and token model for RNA genomes, the sources' sizes and terms |
| Protein-only catalogues (MGnify's proteins, the global gene catalogues) *(card C-45)* | not sized here | deferred: "can wait" (DEC-60 (2)); proteins without neighbours serve homology, not context | a search-only pack beside the corpus, when the PI calls for it |
