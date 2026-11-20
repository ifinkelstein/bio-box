# IMG/VR and Inovirus as stage V2; giant-virus MAGs (GV1) measured but not staged (card C-48c, part 2)

Date: 2026-11-20, wilkcomp02, 09:54 to 12:00 CST (the pod's clock). Every number carries a label: *measured* (host,
date, command or file), *derived* (formula and inputs), *estimate*, or *published* (the provider's figure, with where it
is stated). Host and date of every *measured* number below are wilkcomp02, 2026-11-20, unless said. **No build**:
nothing was written in the data root `/tmp/if2453/r10/root` (bound read-only into every command that read it) and
nothing under `/stor` (the archives bound read-only). No external fetch (source freeze).

Small results are beside this file in [`2026-11-20_virus_v2_gv1_inputs/`](2026-11-20_virus_v2_gv1_inputs/): `look/`
(what the archives hold), `V2/` (`inputs_SUMS`, `report.json`, `check.json`, `objects.json`, `shared_records.tsv`,
`repeated_records.tsv`, `ani.json`, `est.json`, the summaries of the lists left for later and their
`left_for_later_SUMS`), `GV1/`
(`GV1_left_for_later.tsv.gz`, its summary, `gv1.json`), `imgvr_v4.0_full/imgvr-full.json`, `overlap.json`, `run.log`.
Scripts: `scripts/stress/build_stage_pod/c48c_virus_inputs.py` (the selection and the staging; tests
`tests/r1/test_c48c_virus_inputs.py`), `c48c_v2_gv1_look.py`, `c48c_v2_gv1_overlap.py`, `c48c_v2_gv1_ani.py`,
`c48c_v2_gv1_est.py`, `c48c_v2_gv1_left.py`, each run by `c48c_v2_gv1_run.sh <mode>` in one detached unit
`c48c2-<mode>` (pod copy of the tree `/tmp/if2453/c48c2/tree`, log `/tmp/if2453/c48c2/logs/run.log`, outputs
`/tmp/if2453/c48c2/work/<mode>/`).

## Summary

- **Stage V2** is staged at `/tmp/if2453/s2/inputs/V2/` on the pod: **230,098 genomes, 261,626 contigs,
  9,338,542,974 bp, kind virus (gv-meta), redistributable false** (*measured*, `V2/report.json`): one representative
  for each of IMG/VR v4.1's **224,520** high-confidence vOTUs that have a high-quality member (published 231,408: 3.0%
  fewer), chosen by a stated rule since IMG/VR names none, plus the 5,964 Inovirus genomes of Roux et al. 2019. Every
  archive member read has the md5 of `sources.yaml`. Dry check clean; a second run gives every file byte for byte.
- **Shared by `contig_uid`**: 12,762 contigs (974 Mbp) are held by an earlier stage: V1's inputs 10,933 (UHGV took
  IMG/VR v4 genomes), H1 1,677, PL1 68, S1 64, P2 to P5 20 (*measured*). **V2 adds 248,864 contigs, 8,364,370,889 bp.**
  At the vOTU threshold (skani, a sample of 11,434 representatives) 17.9% meet it with a virus of S1, H1 or V1, 11.0%
  near-identically (*measured*); the method's calibration says these are upper bounds (§4).
- **GV1 is not staged: its selection rule is not determined by the data.** GVMAGs V2 on disk is 18,726 genomes
  (4.58 Gbp) and an 8,508-genome subset, with **no table of quality, completeness or tiers** and no README (*measured*),
  so "the MAGs the provider calls high or medium quality" cannot be applied. Gene-call mode, if staged: gv-meta (kind
  virus), the registry's giant-virus models. The PI's choice is needed (§2.3); everything else about GV1 (overlap,
  cost, the list of what is left) is measured on both candidate sets.
- **Cost** (the code's functions and a measured gene-calling sample, *derived*): V2 at `--depth annotated` about
  21,030 CPU-min, of which the Pfam scan is 20,220 in 85 shards (2.54 G residues of about 11.4 M distinct proteins);
  `--depth called` 781, so `called` saves 96%. GV1 (all 18,726) annotated 9,640, Pfam 9,239 in 39 shards; called 389.
- **The rest of IMG/VR v4.0** (8.38 M vOTUs V2 does not take, 80.9 Gbp as one longest member each, *derived*) at
  `--depth called` with `--min-contig-len 10000` would add 2.15 M genomes and 45.0 Gbp for about 4,190 CPU-min and
  about 32 GB (*derived*, §6).
- **What is left for later is written down** (DEC-59, §7): 5,351,677 rows of the HC table and 10,101,426 putative rows
  of v4.0 as sorted lists beside V2's inputs, and GV1's 18,726 genomes in this note's directory.

## 1. What the archives hold (step 1, before any rule)

All read from the archives as streams of their members (`c48c_v2_gv1_look.py`, units `c48c2-look-vr`, `-look-ino`,
`-look-gv`, 10:00 to 10:06), *measured* unless said.

**IMG/VR v4.1 high-confidence** (`jgi_img_vr_v4.1_hc_20261118.zip`, 5 members: the JGI file manifest, the README,
`IMGVR_all_Host_information-high_confidence.tsv` 289,219,376 bytes, `IMGVR_all_nucleotides-high_confidence.fna.gz`
24,792,200,438 bytes stored, `IMGVR_all_Sequence_information-high_confidence.tsv` 1,888,622,625 bytes; sizes equal to
`sources.yaml`; the two tables' md5 equal to it, `look/imgvr_summary.json`; the FASTA's md5 equal to it,
`V2/report.json`). The README names the table's columns: UViG, IMG taxon and scaffold, coordinates (`whole` or a
provirus's range), GOLD ecosystem, vOTU, length, topology, geNomad score, confidence, CheckV completeness and
contamination, MIUViG quality, gene counts, virus taxonomy, host taxonomy and method, origin. The host table has 3
columns (UViG, method, host taxonomy) and one row per UViG.

| | count | bp |
|---|---|---|
| UViGs (rows) | 5,576,197 (published 5,576,197, `sources.yaml`) | 82,353,017,673 |
| vOTUs | 2,917,521 (published 2,917,521) | |
| Confidence `High-confidence` | 5,576,197 (every row) | |
| MIUViG `High-quality` | 486,691 | 19,563,453,979 |
| MIUViG `Genome fragment` | 4,994,325 | 62,156,546,259 |
| MIUViG `Unsure` / `Reference` (RefSeq) | 81,210 / 13,971 | 194,021,225 / 438,996,210 |
| **vOTUs with a high-quality member** | **224,520** (published 231,408, the card: 3.0% fewer) | |
| topology Linear / Provirus / DTR / ITR / GVMAG / Concatemer | 5,137,565 / 297,178 / 131,222 / 8,171 / 2,055 / 6 | |
| origin geNomad / RVMT / RefSeq / Inoviruses / GVMAGs (Schulz 2020) | 5,202,181 / 348,762 / 13,971 / 9,228 / 2,055 | |
| ecosystem Environmental / Host-associated / Engineered / none | 4,297,497 / 879,432 / 269,549 / 129,719 | |
| host prediction none / CRISPR / isolate / k-mer | 4,811,771 / 308,482 / 235,856 / 220,088 | |

**IMG/VR names no vOTU representative**: neither the README nor the table has a representative column, so the rule
below chooses one. Of the 224,520 vOTUs with a high-quality member, 152,504 have one high-quality member and 101,767
have one member at all. A GVMAG UViG is a bin of several FASTA records (5,621,398 records for 5,576,197 UViGs).

**Inovirus** (`Cusrus/data/Gb_files_inoviruses.zip` inside `jgi_inovirus_umag_iuig_meta_20261118.zip`, md5 equal to
`sources.yaml`): 5,964 GenBank files, one record each, `Family/Subfamily/Sp_<id>.gb` (species representatives), 49,235,851
bp, 5,958 distinct sequences, LOCUS topology circular 845 and linear 5,119, lengths 726 to 29,995; families
Protoinoviridae 3,614, Vespertilinoviridae 1,042, Tandem 398, Amplinoviridae 357, Densinoviridae 295, Paulinoviridae
226, Photinoviridae 32 (`look/inovirus_summary.json`). The card quotes no count. IMG/VR's table already holds 9,228
UViGs of origin "Inoviruses".

**GVMAGs V2** (`gvmags_v2_20261118.tar`: `all_18kgenomes_fna.tar.gz`, `fna_files.tar.gz`, `long_ortho_8000.tsv`;
`look/gvmags_summary.json`):

| set | genomes | records | bp | median contigs / bp per genome | under 100 kb |
|---|---|---|---|---|---|
| `confirmed/` in all_18kgenomes_fna | 18,726 (published "18k", `sources.yaml`) | 288,536 | 4,576,307,473 | 11 / 200,420 | 3,698 |
| fna_files | 8,508, every one named in the 18k set | 165,849 | 2,444,841,138 | 14 / 251,515 | 1,132 |

- 8,066 of the 8,508 are byte-identical to their 18k copy; the other 442 have fewer records in fna_files (166.6 against
  161.3 Mbp), as if a contig was removed. `long_ortho_8000.tsv` (1,534,752 rows) names 8,501 genomes, all in fna_files.
- 227 of the 18,726 genomes repeat another's contig set exactly (for example an NCBI assembly under two names), and
  6,626 contigs are in more than one genome (`overlap.json`).
- Headers are `<genome>|<contig>`; genomes include IMG bins (`IMGM…_BIN…`), TARA bins, SRA-assembled bins and isolate
  references (`GCA_…`, named viruses).
- **There is no quality, completeness or tier column anywhere in the archive**, nor a README; the one non-FASTA member
  is an editor backup of a chunking shell script (`confirmed/chunk_maker.sh~`, read as text, never run).

## 2. Selection rules (step 2; tested first)

The tests (`tests/r1/test_c48c_virus_inputs.py`, 9 tests) build small synthetic archives in `tmp_path` and check the
rule against a brute-force sort written in the test, the contig uids against `blake3` directly, md5s against `hashlib`.
They failed first on the pod (8 failed, 1 passed: the GenBank parser had been written for step 1 before its test) and
pass now.

### 2.1 V2: one representative per high-quality, high-confidence vOTU

`c48c_virus_inputs.imgvr_representatives`, a pure function of the table (its result does not depend on the rows'
order, tested over 20 shuffles):

1. Every row must be `High-confidence` (the file is IMG/VR's high-confidence set; another value is refused).
2. A vOTU is taken when at least one of its rows has MIUViG quality `High-quality` (estimated ≥ 90% complete).
   `Reference`, `Genome fragment` and `Unsure` rows are never representatives.
3. Among the vOTU's high-quality rows: the highest CheckV completeness, **capped at 100** (CheckV values above 100 say
   the sequence is longer than its reference genomes, not more complete; 306,588 of the 486,691 are at or above 100),
   then the lowest contamination, then the longest, then the least UViG id. A high-quality row whose completeness or
   contamination is not a number is refused (none was).

Plus every Inovirus GenBank record, files in name order.

### 2.2 How the genomes are listed

Kind virus, gv-meta mode (`genecall.mode_of`, as V1); genetic code 11 recorded (neither source gives one; gv-meta takes
each gene's code from its models); assembly quality `vOTU` (the schema's value for a vOTU representative; the MIUViG
tier, completeness, contamination, topology, ecosystem, host and origin are in `records.tsv`). An IMG/VR genome is every
FASTA record of its UViG (a GVMAG bin: 1,415 genomes of several records), whose bases must equal the table's length
(checked; none differed). Circular only for a one-record genome whose topology is IMG/VR's `Direct terminal repeat` or
GenBank's `circular` and whose sequence no longer carries the repeat (card C-39's test): **IMG/VR keeps the repeat in
59,502 of its 59,608 DTR representatives, which stay linear**; 106 IMG/VR and 478 Inovirus genomes are circular
(*measured* and *derived*: 584 − 106).

### 2.3 GV1: not determined by the data (stopped)

The brief's rule, "the MAGs at the quality tiers the provider calls high or medium", has nothing to read in the
archive (§1). What the data does offer, for the PI to choose from:

| option | genomes | bp | contigs added after V2 | what it rests on |
|---|---|---|---|---|
| (a) all of `confirmed/` | 18,726 | 4,576,307,473 | 264,134; 3,890,643,054 bp | the provider's own folder name: every genome a confirmed giant virus |
| (b) `fna_files` | 8,508 | 2,444,841,138 | 153,020; 2,091,751,441 bp | the subset with orthogroups (long_ortho_8000); 442 genomes one contig shorter than in (a) |
| (c) a completeness we compute | | | | a marker-gene check (e.g. NCLDV core genes) that no card or registry entry provides |

Gene calling, whichever is chosen: the registry offers kinds prokaryote, virus and plasmid and modes single, meta and
gv-meta; kind virus is called in gv-meta mode, pyrodigal-gv's `ViralGeneFinder(meta=True)`, which carries prodigal-gv's
models for giant viruses (and alternative codes). That is the fit, and no registry change is needed. Assembly quality
would be `mag_hq` or `mag_mq` per genome only with a tier; without one, `contig` or a single value the PI names. The
1,174 genomes all of whose contigs V2's IMG/VR GVMAG representatives already hold would be second genomes of those
contigs (DEC-27).

## 3. Stage V2's inputs (step 3)

`c48c_virus_inputs.py v2 HC_ZIP META_ZIP /tmp/if2453/s2/inputs/V2 --root /tmp/if2453/r10/root --root-stages
S1,P1,P2,P3,P4,P5,P6,H1,PL1,C5 --owner-inputs V1=/tmp/if2453/s2/inputs/V1` (unit `c48c2-inputs`, 10:14:50 to 10:44:42,
one thread; peak RSS 354 MB at 10:24, *measured* from `/proc`; the unit's 9.4 GB was page cache). Two passes over the
HC zip's members as streams (the table, then the 24.8 GB FASTA gzip through a hashing reader), the Inovirus zip held
(40 MB); written as `V2.tmp` and renamed. Owners: the root's `contig_index` default variants of the listed stages
(explicit part files; never V1's files in the root, which another agent is building), then V1's staged `records.tsv`.

| File (`/tmp/if2453/s2/inputs/V2/`) | Content | bytes | sha256 |
|---|---|---|---|
| `genomes.json` | 230,098 genomes, 261,626 contigs, kind virus, code 11, 584 circular | 32,448,192 | `bc46f5f954564948ed68aabf3418bd6fa2b11a10d9a5d1eb874be40c09385105` |
| `objects.json` | two objects, `jgi_img_vr` and `jgi_inovirus`, JGI terms, redistributable false | 957 | `583de85a2b1b86665c61fcbad9873df7fdeee6208a034775753f7e4b9f19eee7` |
| `records.tsv` | every record: source, genome, vOTU, contig uid, length, MIUViG, completeness, contamination, topology, code, circular, ecosystem, host method and taxonomy, origin | 67,920,005 | `3baabeda905114da9297607da2be8ec54c8861eae1759da7aeb68893f84de0bb` |
| `repeated_records.tsv` | 386 genomes whose contig set an earlier genome has | 42,730 | `c4b8078d3c3c9729e9d03fb0d0614a08f20364b5e1adcd9a11473716e01de44d` |
| `shared_records.tsv` | 12,765 records whose contig an earlier stage holds, with the owner | 917,801 | `4b4309b87e13ea61019873b896eb442f58697395d6aa0e36f313069e5e9b2c1e` |
| `imgvr_v4.1_hc_hq_votu_reps.fna` | 256,049 records of the 224,520 representatives, headers and line width as in the source | 9,464,839,653 | `69bf5d0984a4c430b25299607d3bacab1b0eea373a4204524515d650032d7e0d` |
| `inovirus_roux2019.fna` | 5,964 records from GenBank, header `<LOCUS> <family>/<subfamily>`, 60 bases a line | 50,500,520 | `209f5a171bab5b2d608ea236e64c959206729a3cb22e2a9c39fd487dc5d22539` |

*Measured* (`V2/inputs_SUMS`, `V2/report.json`). Objects: `retrieved_at` 2026-11-18T23:15:37Z (IMG/VR) and
23:16:08Z (Inovirus), the archive times of `sources.yaml` (the JGI download times were not recorded more finely);
licence "JGI Data Utilization Policy (dataset marked Unrestricted …); held local, not redistributed (DEC-50, DEC-51)".

- **Dry check** (`s2_check_inputs.py V2`, unit `c48c2-check`): exit 0, 230,098 genomes, 261,626 `genome_input` rows,
  230,098 distinct genome uids, 584 circular rows, kinds `[virus]`, codes `[11]`, quality `[vOTU]`, gene-call mode
  gv-meta, no problem (*measured*, `V2/check.json`). No new kind, so `s2_check_inputs.py` is unchanged.
- **Determinism**: the same command into scratch (unit `c48c2-rerun`, 10:45:22 to 11:10:00) gave all seven files byte
  for byte (`cmp`; `objects.json` after its own paths) (*measured*, `run.log`).
- **Representatives** (*measured*): topology Linear 123,336, DTR 59,608, Provirus 38,301, ITR 1,853, GVMAG 1,422; origin
  geNomad 195,158, RVMT 25,586, Inoviruses 2,354, GVMAGs 1,422; ecosystem Environmental 128,240, Host-associated
  69,147, Engineered 13,593, none 13,540; host predicted for 66,634 (CRISPR 35,613, isolate 26,498, k-mer 4,523).
- **Repeats** (DEC-27, *measured*): 386 genomes are listed once under an earlier genome: 22 IMG/VR representatives
  whose bytes equal another vOTU's representative, 358 Inovirus records whose bytes equal an IMG/VR representative, 6
  Inovirus records repeating another Inovirus record. One GVMAG bin lists one contig twice (counted, listed once).

## 4. What V2 and GV1 share by `contig_uid` (step 4)

`c48c_v2_gv1_overlap.py` (unit `c48c2-overlap`, 10:45:22 to 10:46:20), DuckDB on each stage's `contig_index` (owners,
default variants) and `genome_input` (members); GV1's candidate contigs streamed from the tar (281,804 distinct
contigs, 4,316,061,173 bp). All *measured* (`overlap.json`); each stage counted on its own, then "held earlier" is the
union in build order (stages, V1, then V2 for GV1):

| | S1 | P2–P5 | H1 | PL1 | C5 | V1 inputs | V2 | held earlier | **adds** |
|---|---|---|---|---|---|---|---|---|---|
| V2 (IMG/VR) contigs | 21 | 20 | 1,677 | 68 | 0 (2 listed, not owned) | 10,996 | | 12,719; 973.9 Mbp | 243,307; 8,317,869,614 bp |
| V2 (Inovirus) contigs | 43 | 0 | 0 | 0 | 0 | 0 | 358 in IMG/VR reps | 43; 0.3 Mbp | 5,915; 48,864,483 bp |
| **V2** | 64 | 20 | 1,677 | 68 | 0 | 10,996 | | **12,762; 974,172,085 bp** | **248,864; 8,364,370,889 bp** |
| GV1 (a) 18,726 | 0 | 0 | 6 | 0 | 0 | 0 | 17,664; 423,175,502 bp | 17,670; 425.4 Mbp | 264,134; 3,890,643,054 bp |
| GV1 (b) 8,508 | 0 | 0 | 2 | 0 | 0 | 0 | | 12,941; 325.5 Mbp | 153,020; 2,091,751,441 bp |

- P1 and P6 hold none. V2's genomes are held whole or not at all: 12,765 genomes all of whose contigs an earlier stage
  holds (they are second genomes of those contigs, DEC-26/27: no calls of V2's for them), none in part.
- **V1 and V2**: 10,996 contigs (512,003,810 bp) are byte-identical: UHGV took 13,730 IMG/VR v4 genomes (card C-39).
- **V2 and GV1**: 17,567 of the 17,664 shared contigs are contigs of V2's 1,422 IMG/VR GVMAG representatives (which hold
  32,950 contigs): GVMAGs V2 contains much of Schulz et al. 2020's set. Under GV1 (a), 1,174 genomes would be held
  whole and 682 in part.
- Proposed order: V2 before GV1, so that IMG/VR's GVMAG representatives own their contigs and GVMAGs V2's bins of them
  are second genomes; or the PI may prefer GV1 to own every giant-virus contig, in which case V2's 1,422 GVMAG
  representatives would be the second genomes. Either is DEC-27-correct; the counts above are for V2 first.

### At the vOTU threshold (skani)

`c48c_v2_gv1_ani.py` (unit `c48c2-ani`, 3 threads, 10:48:33 to 11:17:24), skani 0.4.0 of the image with card C-39's
settings (`dist -c 30 -m 200 --min-af 0 --qi --ri`). Samples, deterministic: every 20th V2 genome of one record that
repeats no other (11,434; 436.6 Mbp), every 20th stored contig of S1 and H1 (read from their sequence shards in the
root: 1,362 and 335) and of V1's FASTA (8,426). Run A: V2 sample against all of S1 (27,258 stored contigs), H1 (6,714)
and V1 (168,536); run B: the four samples against all of V2 (262,013 records). Rules as card C-39 (`shorter`: ANI ≥ 95,
the shorter ≥ 85% aligned; `both`; `near_identical`: ANI ≥ 99, both ≥ 90%). All *measured* (`V2/ani.json`):

| queries | n | shorter | both | near-identical |
|---|---|---|---|---|
| V2 sample → S1 | 11,434 | 145 | 69 | 29 |
| V2 sample → H1 | 11,434 | 113 | 97 | 96 |
| V2 sample → V1 | 11,434 | 1,891 | 1,453 | 1,149 |
| **V2 sample → S1, H1 or V1** | 11,434 | **2,045 (17.9%)** | **1,585 (13.9%)** | **1,263 (11.0%)** |
| S1 sample → V2 | 1,362 | 311 (22.8%) | 251 (18.4%) | 56 (4.1%) |
| H1 sample → V2 | 335 | 167 (49.9%) | 147 (43.9%) | 122 (36.4%) |
| V1 sample → V2 | 8,426 | 2,355 (27.9%) | 1,423 (16.9%) | 1,121 (13.3%) |
| calibration: V2 sample → other V2 | 11,434 | 1,626 (14.2%) | 998 (8.7%) | 244 (2.1%) |

- What it shows: V2's overlap with the root is mostly with V1 (UHGV's IMG/VR genomes and their relatives); about
  11% of V2's representatives are near-identical to a virus V1, S1 or H1 holds, against 4.8% byte-identical (*derived*:
  10,996 / 230,098), since proviruses cut at other coordinates and rotated or trimmed genomes are not the same bytes.
  Extrapolated (*derived*: fraction × 230,098 genomes): about 25,400 (near-identical) to 41,200 (shorter rule) V2
  representatives are vOTUs of a virus the root or V1 holds. Half of H1's huge phages have an IMG/VR high-quality vOTU.
- What it does not show: the calibration (V2 representatives against other V2 records) meets the 95/85 rule in 14.2%
  even though each vOTU contributes one representative; part of that is Inovirus records against the IMG/VR
  representative of the same virus (both in V2 by design), the rest is skani's estimate near the line, as card C-39
  found for UHGV. So the `shorter` and `both` counts are upper bounds and `near_identical` the firmer figure. GV1 was not
  sampled (it is not staged; giant-virus genomes of several contigs need genome-level ANI, `--qi` off).
- CPU: run A 1,771.9 s, run B 2,563.2 s of skani (*measured*, getrusage of the children).

## 5. Dry check (step 5)

See §3: `s2_check_inputs.py V2` exit 0, no problem. GV1 has no inputs to check.

## 6. What building each stage would cost (step 6)

`c48c_v2_gv1_est.py` (units `c48c2-est-v2`, `c48c2-est-gv1`, one thread each): the code's own functions from the tree.
The distinct proteins are not known before calling, so a sample (every 40th V2 genome called; every 50th GV1 genome,
name order) was called here with the step's own path (`genecall._call`, one thread) and its genes, distinct proteins
and residues per Mbp scaled to the stage's stored bases: *derived*, and an over-estimate of distinct proteins (a sample
repeats fewer proteins than the whole).

| Step | V2 (230,098 genomes; 8.364 Gbp stored) | GV1 (a) 18,726 (3.891 Gbp stored) | GV1 (b) 8,508 (2.092 Gbp) | How |
|---|---|---|---|---|
| ingest | 35.2 (39.2 counting every byte) | 18.4 | 9.8 | *derived*: `stage.ingest_cpu_min(FASTA bytes, held bases)`; V2's objects 9,515,340,173 bytes, 974,172,085 bases held |
| gene calling, gv-meta | 722.2 in 4 batches | 359.9 in 2 | 191.1 in 1 | *derived*: `genecall.one_process_cpu_min` over the bases read (V2: 8,354,292,426 linear, 26,054,307 circular with flanks; every genome with a contig of its own called whole) |
| genes | 11.48 M (code rate 12.05 M) | 4.07 M | 2.19 M | *derived*: sample genes per Mbp × stored Mbp |
| distinct proteins; residues | 11.40 M; 2,540.6 M | 4.01 M; 1,160.6 M | 2.15 M; 624.0 M | *derived*: sample per Mbp × stored Mbp (upper bound) |
| Pfam 38.2 scan | **20,223.5 in 85 shards** | **9,238.6 in 39** | 4,967.0 in 21 | *derived*: residues × 7.96 (`SCAN_CPU_MIN_PER_M`); shards of 30,150,753 residues |
| CRISPR, default and near-Cas | 23.5 and 27.0 | 10.9 and 12.5 | 5.9 and 6.7 | *derived*: `stage.crispr_cpu_min(stored, flanks)`; V2's circular flanks 18,951,436 bp |
| **in all, `--depth annotated`** | **21,031** (351 CPU-hours) | **9,640** | 5,181 | *derived*: the sum |
| **`--depth called`** | **781** | **389** | 207 | *derived*: ingest, calls, default CRISPR |
| what `called` saves | 20,250 (96%) | 9,251 (96%) | 4,974 | *derived*: Pfam and the near-Cas CRISPR producer |
| files in the root | about 6.0 to 6.1 GB | about 2.8 GB | about 1.5 GB | *derived*: stored Gbp × 0.7175 (S1) to 0.7316 (H1) GB per Gbp (card C-39, *measured* there) |

**The samples** (one thread, `genecall._call`, *measured*): V2 every 40th called genome, 5,442 genomes, 210,261,281 bp:
288,695 genes (1,373.0 per Mbp), 286,610 distinct proteins (1,363.1 per Mbp; 40,891 of them already in the root's
`protein_seq`, which since today includes V1's), 303,746 residues per Mbp, calls 973.8 CPU-s = 4.631 s per Mbp. GV1
every 50th genome of the 18,726, 374 genomes, 90,017,161 bp: 94,106 genes (1,045.4 per Mbp), 92,670 distinct
(1,029.5 per Mbp; 482 in the root), 298,313 residues per Mbp, 2.605 s per Mbp (`V2/est.json`, `GV1/gv1.json`).

**Against the code's rates.** The code's gv-meta linear rate is 4.2019 s per Mbp (the largest measured, rounded up,
`genecall.CALL_S_PER_MBP`). V2's sample ran at 4.631, **10.2% above it** (*derived*), while four of our processes and the
V1 build shared the host; the step's estimate multiplies the calls by `STEP_OVERHEAD` 1.23, so the batch sizes still
hold (5.17 per Mbp with the overhead, *derived*). GV1's sample ran at 2.605, **38% under** the code's rate: giant-virus
genomes are long and gene-dense in long genes (1,045 genes per Mbp against V2's 1,373), so the caller spends less per
base. The code errs on the safe side for GV1; AGENTS §6's 20% is passed, reported here and in §9, not acted on (a code
card's call).

No single process passes 300 CPU-minutes: gene calling runs in batches of at most 240 (code), each Pfam shard at most
240, ingest 39.2 at most.

**The rest of IMG/VR v4.0** (*derived* from its table alone, `imgvr_v4.0_full/imgvr-full.json`; the table copied out
of `jgi_img_vr_v4_full_20261118.zip`, md5 equal to `sources.yaml`; its 48 GB FASTA was not read): 15,677,623 UViGs
(published 15,677,623), 8,606,551 vOTUs (`sources.yaml` says 9,169,185: 6.1% fewer), 160.2 Gbp; 10,101,426 putative rows
(77.8 Gbp). The vOTUs V2 does not take number 8,382,031; one longest member each is 80.9 Gbp. Staging them at `--depth
called` (ingest, gv-meta calling, CRISPR; no Pfam) with `--min-contig-len` L would add, every genome linear:

| L | vOTUs left | bp | CPU-min (ingest + calls + CRISPR) | gene-call batches | root GB at 0.71 GB/Gbp |
|---|---|---|---|---|---|
| 5 kb | 5,913,093 | 70.68 G | 6,582 | 26 | 50.2 |
| 10 kb | 2,149,481 | 44.95 G | 4,186 | 17 | 31.9 |
| 20 kb | 694,776 | 25.15 G | 2,342 | 10 | 17.9 |

The longest member over-counts a representative's bases, and 0.71 GB per Gbp is DEC-49's figure for the called depth
(*derived* there from H1); the CPU figures use the code's gv-meta linear rate, which the samples here run under (§6).

## 7. What this selection leaves for later (DEC-59)

Nothing was dropped silently: every sequence of the three archives that V2 does not stage is in one of these lists,
with the provider's id, vOTU, quality tier, length, reason, archive and member. The two large lists live beside V2's
inputs on the pod (not committed); GV1's is in this note's directory. Each list is sorted (by reason, then id), written
through a temporary name with gzip at mtime 0; a second run of the HC list gave the same bytes (unit
`c48c2-left-rerun`, *measured*). Columns: `sequence_id votu quality length records reason archive member`.

| list | rows | bytes | sha256 | where |
|---|---|---|---|---|
| IMG/VR v4.1 HC rows V2 does not take | 5,351,677 | 58,947,247 | `ace2ec9845b0cea6e6f4c9f81c9baaf7b7b348f90440e0dda89d69de86ad5b62` | `/tmp/if2453/s2/inputs/V2/V2_left_for_later_hc.tsv.gz` |
| IMG/VR v4.0 putative rows | 10,101,426 | 106,448,343 | `d607ba78f5ed09342f6095e4b471b9f7703374f4757a4b070f4f3861f0c80751` | `/tmp/if2453/s2/inputs/V2/V2_left_for_later_v4.0_putative.tsv.gz` |
| GVMAGs V2 genomes (GV1 undecided) | 18,726 | 237,910 | `adb42e8d3dc5e05ea6c0f33d2d25f4e968cb59db30ac9ec379126237c12a7f4c` | `GV1/GV1_left_for_later.tsv.gz` here |

By reason and tier (*measured*, `V2/left_for_later_*.json`, `GV1/left_for_later.json`):

| reason | tier | rows | bp |
|---|---|---|---|
| HC: `not_representative` (another member of a vOTU V2 takes) | High-quality | 262,171 | 10,269,399,888 |
| | Genome fragment | 268,574 | 4,994,208,132 |
| | Unsure / Reference | 7,053 / 886 | 12,573,259 / 59,017,545 |
| HC: `votu_without_hq_member` (2,693,001 vOTUs) | Genome fragment | 4,725,751 | 57,162,338,127 |
| | Unsure / Reference | 74,157 / 13,085 | 181,447,966 / 379,978,665 |
| v4.0: `putative` (6,251,664 vOTUs with a putative row; 5,689,030 putative only) | Genome fragment | 8,280,573 | 66,015,149,311 |
| | High-quality / Unsure | 51,203 / 1,769,650 | 901,445,511 / 10,931,848,757 |
| GVMAGs V2: `gv1_selection_undecided` | none in the archive | 18,726 (288,536 records) | 4,576,307,473 |

Check (*derived*): HC rows left 5,351,677 + taken 224,520 = 5,576,197, and bases 73,058,963,582 + 9,294,054,091 =
82,353,017,673, the table's. `records` is 1 per IMG/VR row; a GVMAG UViG's FASTA has several records.

Skipped by the staging's own rules (in V2's inputs, as V1's): `repeated_records.tsv` (386 genomes listed under an
earlier genome of the same bytes), `shared_records.tsv` (12,765 records whose contig an earlier stage owns; they stay
genomes of V2, DEC-26/27). No record was malformed or refused; nothing else was skipped.

**Picking each remainder up later** (*derived* from the lists, *estimate* where said):

- **HC non-representative members** (538,684 UViGs, 15.3 Gbp): the same member of the same zip
  (`IMGVR_all_nucleotides-high_confidence.fna.gz`), selected by the list's ids (`reason = not_representative`); one
  streaming pass as V2's (30 min, one thread, *measured* for V2). The 262,171 high-quality members (10.3 Gbp) are
  strain-level copies of V2's vOTUs; at `--depth called` about 960 CPU-min for them (*derived*: 10.27 Gbp at the code's
  gv-meta linear rate × 1.23, 885, plus ingest 43 and CRISPR 29, every genome linear). Suggested: a stage of the high-quality members at `called`, the fragments with
  `--min-contig-len`.
- **HC vOTUs without a high-quality member** (2,693,001 vOTUs, 4.81 M rows, 57.7 Gbp): same member, by
  `reason = votu_without_hq_member`; one representative each by V2's rule extended to all tiers, then `--depth called`
  with `--min-contig-len` (§6's table counts these together with the putative-only vOTUs).
- **v4.0 putative** (10.1 M rows, 77.8 Gbp): `jgi_img_vr_v4_full_20261118.zip`, member
  `Cus_VR/IMG_VR_2022-12-19_7/IMGVR_all_nucleotides.fna.gz` (48.3 GB), by the list's ids; `--depth called` with
  `--min-contig-len` only (geNomad did not call them high-confidence).
- **GVMAGs V2**: `gvmags_v2_20261118.tar`, member `gvmags_v2/all_18kgenomes_fna.tar.gz` (or `fna_files.tar.gz`), each
  genome one `.fna`; 18,726 genomes, 3.89 Gbp new after V2; `annotated` costs about 9,640 CPU-min, `called` about 390
  (§6). Needs the PI's selection rule first (§2.3).

## 8. Proposed `sources.yaml` lines (not applied)

Under `jgi_img_vr`:

```yaml
    staged:
      - stage: V2
        object: /tmp/if2453/s2/inputs/V2/imgvr_v4.1_hc_hq_votu_reps.fna
        bytes: 9464839653
        sha256: 69bf5d0984a4c430b25299607d3bacab1b0eea373a4204524515d650032d7e0d
        derived_from: jgi_img_vr_v4.1_hc_20261118.zip, members IMGVR_all_Sequence_information-high_confidence.tsv
          (md5 3c516db128082fa29dc2c2f60520da1b) and IMGVR_all_nucleotides-high_confidence.fna.gz (md5 83301c9c6dfefea3305a53ee2a41bac3)
        rule: c48c_virus_inputs.py v2 (one representative per vOTU with a High-quality member; card C-48c)
        left_for_later: /tmp/if2453/s2/inputs/V2/V2_left_for_later_hc.tsv.gz (sha256 ace2ec98…5ad5b62),
          /tmp/if2453/s2/inputs/V2/V2_left_for_later_v4.0_putative.tsv.gz (sha256 d607ba78…0c80751)
```

Under `jgi_inovirus`:

```yaml
    staged:
      - stage: V2
        object: /tmp/if2453/s2/inputs/V2/inovirus_roux2019.fna
        bytes: 50500520
        sha256: 209f5a171bab5b2d608ea236e64c959206729a3cb22e2a9c39fd487dc5d22539
        derived_from: jgi_inovirus_umag_iuig_meta_20261118.zip, member Cusrus/data/Gb_files_inoviruses.zip
          (md5 a075958a353ec9a36c01264513f329e1), every GenBank record as FASTA
```

Under `gvmags_v2`: `staged: none (card C-48c part 2: no quality tiers in the archive; selection waits on the PI)`.

## 9. Open questions

1. **GV1's selection** (§2.3): (a) all 18,726, (b) the 8,508 of `fna_files`, or (c) a completeness computed by us; and
   the assembly quality to record. And whether V2 or GV1 owns the 17,664 shared giant-virus contigs (§4).
2. **Circular DTR genomes**: IMG/VR keeps the terminal repeat in 59,502 of 59,608 DTR representatives, so they are
   listed linear (the C-39 rule); a gene across the repeat junction is then split. Trim the repeat (a change of bytes,
   so of `contig_uid`) or accept linear.
3. **Quality value**: `vOTU` for every V2 genome; the MIUViG tier is only in `records.tsv`.
4. **The code's gv-meta rate** (§6): V2's sample ran 10.2% above the code's calls-only rate (inside the step overhead),
   GV1's 38% under it. A code card may add a giant-virus rate; nothing here depends on it.
5. `sources.yaml` says 9,169,185 vOTUs for v4.0; its table has 8,606,551 distinct vOTU ids (6.1% fewer).

## 10. Pod use and scratch

| Unit | CST | What | Threads |
|---|---|---|---|
| `c48c2-look-vr` | 10:00 (failed in its summary), 10:04 to 10:06 | the HC tables copied out, md5, summary | 2 |
| `c48c2-look-ino`, `-look-gv` | 10:00 to 10:06 (ino stopped once for a slow nested-zip read and rerun) | Inovirus, GVMAGs | 1 each |
| `c48c2-inputs` | 10:14:50 to 10:44:42 | V2's inputs | 1 |
| `c48c2-est-full` | 10:19:54 to 10:21:39 | v4.0's table, estimate | 2 |
| `c48c2-left-put`, `-left-gv1` | 10:29 (failed: an API name), 10:32:57 to 10:35:55 | lists left for later | 1 each |
| `c48c2-check`, `-rerun`, `-overlap`, `-left-v2` | 10:45:21 to 11:10:00 | dry check, determinism, overlap, HC list | 1 each |
| `c48c2-ani` | 10:48:33 to 11:17:24 | skani | 3 |
| `c48c2-est-gv1`, `-left-rerun`, `-est-v2` | 11:10:38 to 11:39:12 | estimates, list determinism | 1 each |
| `pod_pytest.sh C-48c-2 --shards 2` | from 11:40 | the whole suite | 2 shards |

Every unit: `systemd-run --user`, `nice 19`, `ionice` 7, a memory cap (8 to 32 GB), runtime image
`bio-box_faec1f4.sif`, `PYTHONPATH` the tree's `src` alone, the data root, the archives and V1's inputs bound read-only;
at most 4 threads of ours at once. No `STOP` in the data root at any step. The other agent's units (`c48c-*`) were not
touched.

Disk (*measured*, `du -sb`, 11:41): V2's inputs with the two lists 9,782,065,448 bytes; scratch `/tmp/if2453/c48c2`
8,423,810,377 (the copied HC and v4.0 tables 7.2 GB, skani's samples and tables, the tree): 18.2 GB in all, under the
brief's 100 GB. Left on the pod: V2's inputs and the two lists beside them, the
scratch `/tmp/if2453/c48c2` (the copied tables and the skani samples). Nothing on `/stor`.

## 11. Proposed register rows (for `planning/CATALOG_GOAL.md` §4.1)

| Source | Taken | Left for later, in full | Where the remainder is recorded | How to pick it up |
|---|---|---|---|---|
| IMG/VR v4.1 high-confidence (`jgi_img_vr`) | V2: 224,520 vOTU representatives (one per vOTU with a high-quality member), 9.29 Gbp | 538,684 other members of those vOTUs (15.3 Gbp; 262,171 high-quality) and 2,693,001 vOTUs without a high-quality member (4.81 M UViGs, 57.7 Gbp) | `/tmp/if2453/s2/inputs/V2/V2_left_for_later_hc.tsv.gz` (sha256 `ace2ec98…`); this note §7 | same zip member by id; `--depth called`, `--min-contig-len` for fragments |
| IMG/VR v4.0 full, putative (`jgi_img_vr`) | none | 10,101,426 putative UViGs (77.8 Gbp; 5,689,030 putative-only vOTUs) | `/tmp/if2453/s2/inputs/V2/V2_left_for_later_v4.0_putative.tsv.gz` (sha256 `d607ba78…`); §6, §7 | `jgi_img_vr_v4_full_20261118.zip` FASTA member by id; `--depth called` with `--min-contig-len` |
| Inovirus (`jgi_inovirus`) | V2: all 5,964 records (5,600 genomes listed; 364 repeat an IMG/VR representative or each other, DEC-27) | none | `V2/repeated_records.tsv` | — |
| GVMAGs V2 (`gvmags_v2`) | none (no quality tiers in the archive) | 18,726 genomes, 4.58 Gbp (3.89 Gbp new after V2) | `2026-11-20_virus_v2_gv1_inputs/GV1/GV1_left_for_later.tsv.gz` (sha256 `adb42e8d…`); §2.3 | PI picks the set; `all_18kgenomes_fna.tar.gz` per genome; kind virus, gv-meta |
