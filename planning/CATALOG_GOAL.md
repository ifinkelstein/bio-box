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
- What "useful" excludes is the PI's to say (card C-45: RNA viruses, protein-only catalogues, raw reads). Until the
  PI says, those are listed in §4 as deferred, not as out.

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
| UHGV v1.0 vOTU representatives (168,536) | building | `V1` |
| GTDB R226 complete genomes, second selection (13,477 genomes, 56.18 Gbp) | `C5` resuming; `C1` to `C4`, `C6` queued | `C1` to `C6` |
| Plasmids from metagenomes (Yu 2024, PlasmidScope, IMG/PR; 859,036 new contigs) | inputs staged | `PL2` |
| IMG/VR v4.1 high-confidence HQ vOTU representatives, Inovirus, GVMAGs V2 | inputs being staged | `V2`, `GV1` |
| GlobDB r226 non-GTDB archaea (7,139 genomes) | inputs staged | `G1` |
| IUIG novel prokaryotes (about 4,235 representatives, *derived*, card C-48d) | archive on disk; not staged | `J1` |
| IMG-UMAG eukaryote bins (691 bins, 14.44 Gbp) | archive on disk; waits on a `"eukaryote"` genome kind | `E1`, `E2` |

## 4. The register of what is deferred

Each row is something a selection left out, or a source not yet fetched. "How to pick it up" names what a later
card needs. Rows are added, and a row is struck through with its stage's name when its remainder is in the corpus.

### 4.1 Left out by a selection inside a source we hold

| Source | Taken | Left for later, in full | Where the remainder is recorded | How to pick it up |
|---|---|---|---|---|
| GTDB R226 genomes (732,475 rows *(card C-48d)*) | 11,966 best complete genomes per species (`P1` to `P6`); 13,477 further complete genomes (`C1` to `C6`) | every other GTDB genome: the species representatives that are not complete genomes (plan stages S3 and S5; 143,607 representatives, 465.47 Gbp, are also inside GlobDB's archive *(card C-41)*) and all non-representative genomes (plan stage S6) | `planning/measurements/2026-11-11_r50_inputs.md`, `2026-11-15_s2c_phase2_inputs.md` (the selections and their `remainder.tsv`, `rest.tsv`); `sources.yaml` | fetch by assembly accession from NCBI datasets (cards C-25, C-27 give the form); needs item 19 answered: S6 is 111.5 M contigs by the plan's count |
| UHGV v1.0 | one representative per vOTU (168,536) | the non-representative genomes of each vOTU | `planning/measurements/2026-11-16_virus_votu_inputs.md` §3 (`uhgv_metadata.tsv` has every genome's vOTU) | fetch the provider's full genome FASTA; stage as a further `V` stage, `contig_uid` dedup against `V1` |
| IMG/VR v4.1 and v4.0 (5,576,197 genomes, 2,917,521 vOTUs high-confidence; 9.17 M vOTUs in v4.0, *published*, *(card C-48c)*) | HQ high-confidence vOTU representatives (`V2`, being staged) | all other vOTU representatives (next in §2's order), then non-representative genomes | card C-48c part 2's note (`planning/measurements/2026-11-20_virus_v2_gv1_inputs.md`, when written); both archives are on `/stor/work` | no fetch needed: stage from `jgi_img_vr_v4_full_20261118.zip` at depth `called`; crosses 7.86 M contigs in the root, so card B-56's design (open item 16) comes first |
| Plasmid sets of `PL2` | 859,036 new contigs | 15,135 circular gut contigs that are not representatives (PI question, card C-40); records skipped by the stage's rules | `planning/measurements/2026-11-16_pl2_inputs.md`; `/tmp/if2453/s2/inputs/PL2/skipped_records.tsv` (on `/tmp` only: copy into the note's directory before `PL2` is built) | a further plasmid stage from the same FASTA files |
| IUIG (143,796 genome rows *(card C-48d)*) | about 4,235 novel prokaryote representatives proposed (`J1`) | 120,654 prokaryote genomes that match a GTDB or GlobDB species (strain copies); 15,544 virus, 1,185 plasmid and 721 eukaryote rows | `planning/measurements/2026-11-19_jgi_mags_iuig_dedup.md` | no fetch needed: `jgi_iuig_genomes_20261118.zip` (524.8 GB) is on `/stor/work`; stream by member |
| IMG-UMAG prokaryote bins (18,912 metagenomes) | none yet; novel representatives proposed (`J2` to `J15`, 120 to 140 Gbp, *estimate*) | bins that match a held species | the same note | no fetch needed for the bins; their dereplication wants GlobDB's and GTDB's representatives on disk, or a decision to ingest without it |
| GEMs (52,515 MAGs, 45,599 OTUs *(card C-48d)*) | none as genomes (every OTU is already judged by GlobDB); the 24,313 VirSorter prophage coordinates as a sidecar | all 52,515 MAGs as genomes | the same note | no fetch needed: `jgi_gems_52515_fna_20261118.tar` is on `/stor/work` |
| SPIRE pilot (990 assemblies, 341.5 Gbp raw) | measurements only; nothing built into the root | the assemblies themselves (unlinked after hashing) | `planning/measurements/2026-11-18_spire_pilot.md`, `select/pilot_selection.tsv` | fetch again by sample id (card C-48e) |

### 4.2 Left out by a filter

| Filter | Decision | What it leaves out | How to pick it up |
|---|---|---|---|
| Least contig length 5,000 bp on unbinned metagenome stages | DEC-57 | contigs under 5 kb: 83.9% of assembled bases and 98.9% of contigs in the pilot (*derived*: 1 − 55.050 / 341.501 Gbp; 1 − 4.81 M / 453.49 M) | each tranche's note lists its samples and the filter; the raw assemblies are removed after staging, so the remainder is fetched again by sample id and ingested with a lower `--min-contig-len` as a further stage. **Open, the PI's**: whether "in full" includes contigs under 5 kb, under 1.5 kb, or neither (DEC-54 (11) says a dropped contig is outside the corpus; DEC-59 makes that "outside for now") |
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
| RNA viruses (LucaProt, Serratus, RVMT); protein-only catalogues *(card C-45)* | not sized here | everything built assumes DNA genomes with neighbourhoods | the PI's word on whether "useful" includes them |
