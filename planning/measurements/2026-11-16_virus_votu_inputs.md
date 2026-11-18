# Viruses from metagenomes as vOTU representatives: inputs of stage V1 and what they share (card C-39)

Dates: 2026-11-16, wilkcomp02, 13:44 to 15:33 CST, and 2026-11-17 08:15 (the pod's clock). Every number carries a
label: *measured* (host, date, command or file), *derived* (formula and inputs), *estimate*, or *published* (the
provider's figure, with where it is stated). **No build**: nothing was written in the data root
`/tmp/if2453/r10/root` (bound read-only into every command that read it) and nothing under `/stor`.

Small results are beside this file in [`2026-11-16_virus_votu_inputs/`](2026-11-16_virus_votu_inputs/):
`download_SUMS` (every fetched file: bytes, md5, sha256, time), `fetch.log`, `vire.json` (`c39_vire.py`), and `V1/`
(`report.json` of `c39_inputs.py`, `check.json` of `s2_check_inputs.py`, `inputs_SUMS`, `objects.json`,
`shared_records.tsv`, `overlap.json` of `c39_overlap.py`, `ani.json` of `c39_ani.py`, `est.json` of `c39_est.py`,
`run.log`). Scripts: `scripts/stress/build_stage_pod/c39_*`.

## Summary

- **Fetched** (pod, 2026-11-16, `c39_fetch.sh`): UHGV v1.0's vOTU representatives and tables from Zenodo, every md5
  equal to the record's; VIRE's all-samples genome table from EMBL. **The Global Soil Virus Atlas was not fetched**:
  PNNL's hosts refused the pod (TLS reset at the client hello on the Globus URL; HTTP 400 to curl on the README), five
  refusals in a row, so the fetch stopped as the entry's limits say. IMG/VR was not tried (JGI, DEC-50).
- **Stage V1** (UHGV): `/tmp/if2453/s2/inputs/V1/` on the pod, **168,536 genomes, one contig each,
  5,222,288,776 bp, kind virus** (gv-meta), 22,517 circular, 24,404 `complete` and 144,132 `vOTU` (*measured*,
  `V1/report.json`). Dry check clean; a second run gives every file byte for byte (*measured*).
- **Shared by `contig_uid`** with the stored stages: **178 contigs** (H1 177, S1 1; 43,795,699 bp); none with P1 to P6
  or PL1 (*measured*). The stage adds **168,358 contigs, 5,178,493,077 bp**.
- **At the vOTU threshold** (skani, a sample of 8,426 V1 representatives): 163 (1.9%) meet ANI ≥ 95% over ≥ 85% of
  the shorter genome with a phage of S1 or H1, 67 (0.8%) when both genomes must be ≥ 85% covered, 28 (0.3%)
  near-identical (*measured*). Of S1's phages 18.9% (shorter-genome rule) or 13.5% (both covered) fall in a UHGV vOTU;
  of H1's, 19.7% or 15.8% (*measured*, samples of 1,703 and 361). The method's own calibration (V1 against V1) says
  the 95/85 rules are loose with skani near the threshold (below).
- **Cost of building V1** (the code's own functions, *derived*): ingest 23.4 CPU-min; gene calling 687.5 CPU-min in
  3 batches (the measured sample gives 550, 20.0% under the code: see §5); CRISPR 15.7 and 18.7; Pfam about 10,400
  CPU-min in 44 shards (1,307 M residues of 5.60 M distinct proteins, *measured* on the provider's proteins); 3.7
  to 4.2 GB in the root.
- **Proposal: two stages, `V1` (UHGV) and `V2` (the soil atlas), not one** (§2). V2 has no inputs yet.
- **VIRE** (table only): 1,784,510 genomes in 697,136 species clusters, of which 255,274 human gut, 33,492 soil and
  408,370 from other habitats; 567,044 clusters have no member VIRE aligned to IMG/VR (*measured*, `vire.json`).

## 1. What was fetched, and what was not

Entries (`sources.yaml`): `uhgv` and `gsv_atlas` as approved by the PI, with two corrections made here, and two
additions written **before** their fetch under DEC-50 (commits `45d6ad0` and `8185dac` on `task/C-39`, before the
fetches at 19:49 and 20:31 UTC): `uhgv_metadata.tsv.zst` and `votus_full.faa.zst` in `uhgv`, and the new entry
`vire_tables`. One connection, `-A curl`, no name, address or project detail in any request; nothing fetched was
executed; zstd files were decompressed by Python's `compression.zstd` in the image, never unpacked as archives.

| Source | File | Bytes | Check | Fetched (UTC) |
|---|---|---|---|---|
| uhgv | `votus_full.fna.zst` | 964,680,445 | md5 = Zenodo's | 19:49:18 to 19:53:06 |
| uhgv | `votus_metadata_extended.tsv.zst` | 14,048,300 | md5 = Zenodo's | 19:53 |
| uhgv | `source_biosample_metadata.tsv.zst` | 152,063 | md5 = Zenodo's | 19:53 |
| uhgv | `host_range_breadth.tsv.zst` | 490,562 | md5 = Zenodo's | 19:53 |
| uhgv | `uhgv_metadata.tsv.zst` (added) | 38,978,745 | md5 = Zenodo's | 19:54 |
| uhgv | `votus_full.faa.zst` (added) | 1,022,105,075 | md5 = Zenodo's | 20:31 to 20:32 |
| vire_tables | `vire_all_samples_metadata.txt` | 903,989,625 | none published; ours recorded | 20:00:18 to 20:01:32 |
| gsv_atlas | archive and README | none | refused | 19:54 to 19:57 |

All *measured* (`download_SUMS`, `fetch.log`); sha256 values are in `sources.yaml`. UHGV's files need a zstd decoder
with a long window (`window_log_max` 31): the first inputs run failed on it ("Frame requires too much memory").

**Corrections to the approved entries.** `uhgv`: the URL form `https://zenodo.org/records/17402089/files/<name>?download=1`
was read from the record page and worked (the entry had it as unconfirmed). `gsv_atlas`: the README URL on
`datahub.pnnl.gov` redirects (302) to the site's front page; the same path on `data.pnnl.gov` serves it to a page
reader, and was written into the entry.

**The soil atlas refused us.** The HEAD of the archive URL (`g-83fdd0.1beed.03c0.data.globus.org`, 192.101.102.34)
timed out in the TLS handshake (curl 28), and a second try was reset by the peer at the client hello (curl 35); the
README on `data.pnnl.gov` answered HTTP 400 (247 bytes of HTML) three times to curl from the pod (*measured*,
`fetch.log`). Five failures in a row: the fetch stopped there, and nothing was retried with another client, user
agent or TLS setting, since that would be working around a refusal. The page reader (not the pod) could read both the
DataHub page and the README; what this note says of the atlas's contents is from the README (*published*). What the
PI can decide: fetch the archive from another network (a browser on campus, or a Globus transfer from the PI's
account) to `/tmp/if2453/s3/dl/gsv/`, or let an agent fetch it from the laptop. With the archive, V2's inputs are one
more mode of `c39_inputs.py` (the README names a vOTU column; the representative must be chosen by a rule if the
table does not name it).

## 2. One stage or two

Proposed: **`V1` for UHGV and `V2` for the soil atlas** (kind `virus` both). The terms do not force it (CC BY 4.0 and
CC0, both redistributable). The reasons:

- The schema carries no habitat or source per genome (`genome_dim`: kind, code, quality, tier, species-rep); the
  source rides only on `contig_occurrence` through the source object. The stage is the one label a release reader can
  filter genomes by without a join to occurrences, so "human gut viruses" and "soil viruses" stay separable.
- V1 is ready and V2 waits on the PI (§1). DEC-52's rolling releases take each complete stage when it is done.
- Each source's quirks stay in its own inputs (UHGV's long-window zstd and DTR trimming; the atlas's per-contig
  table without named representatives).
- Size does not force it: gene calling is batched now (card B-48), and V1 alone takes 3 batches.

IMG/VR, when the PI fetches it, is a third stage with `redistributable = false` (S4 of `STAGES.md`).

## 3. Stage V1's inputs

Script: `c39_inputs.py uhgv FNA_ZST VOTU_META_ZST OUT_DIR RETRIEVED_AT --root ROOT`, run by `c39_run.sh inputs` (unit
`c39-inputs`, 14:33:19 to 14:34:22 CST). It decompresses the FASTA byte for byte into the inputs directory as
`uhgv_v1.0_votus_full.fna` (the one source object, 5,228,187,536 bytes), joins each record to its row of
`votus_metadata_extended.tsv` by the header's first token (refusing a record without a row, a row without a record,
a length that differs from `genome_length`, or a value outside the maps below), and reads every stage's
`contig_index` for the owners (`make_inputs.stage_owners`, DEC-26).

| File (`/tmp/if2453/s2/inputs/V1/`) | Content | sha256 |
|---|---|---|
| `genomes.json` | 168,536 genomes, one contig each, kind virus, codes 11, 15 and 4, 22,517 circular | `752ac47ff81cca4432d88574a0183427fe76cb2d1b7148e90a0e4ce85c398240` |
| `objects.json` | the FASTA, source `uhgv`, release `zenodo:10.5281/zenodo.17402089 v1.0`, CC BY 4.0, redistributable, retrieved 2026-11-16T19:53:06Z (the end of its fetch) | `dc7bed72ef15ccb93f41da05a8d600f8b05cad9312993fe2723c2e5761864d6c` |
| `records.tsv` | per record: ordinal, genome, vOTU, contig_uid, length, CheckV quality and method, code, circular | `31ef5a640048d4bae72e947c380fb00ac0268f057114250c613199155e87a7cf` |
| `repeated_records.tsv` | records whose bytes repeat an earlier record's: none | `4142a62fae1678e3272cf972ff05eaa52c5ddaa4ea49e196d68529365e109cf6` |
| `shared_records.tsv` | the 178 records whose contig a stage owns, with the stage | `96216d41981855ecc6d1aad0094628d6421844180a5c4e96c757374bd9c8634d` |
| `uhgv_v1.0_votus_full.fna` | the FASTA | `9c595d262d89100f25ebe40964ee31912b26d9b64c6b5f453daf73b931071356` |

*Measured* (`V1/inputs_SUMS`). **Dry check** (`s2_check_inputs.py V1`, unit `c39-check`): exit 0, 168,536 genomes,
168,536 `genome_input` rows and distinct genome uids, 22,517 circular rows, kinds `[virus]`, codes `[4, 11, 15]`,
quality `[complete, vOTU]`, gene-call mode gv-meta, no problem (*measured*, `V1/check.json`). **Determinism**: the
same command again into scratch (unit `c39-rerun`) gave every file byte for byte (`cmp`; `objects.json` after its own
path) (*measured*, `V1/run.log`).

The FASTA holds exactly the 168,536 representatives the table names (`votu_representative` Yes in
`uhgv_metadata.tsv`), and the vOTU table has one row per record (*measured*, `V1/overlap.json`). UHGV's record says
873,995 genomes; its table has 873,994 rows (*measured*; one fewer, not reconciled; the representatives are complete).

**Fields and the rules applied** (all from `votus_metadata_extended.tsv`):

- **Genetic code**: `genetic_code` (prodigal-gv's pick): 11 166,548; 15 1,804; 4 184 (*measured*). Recorded in
  `genome_input`; in gv-meta mode the code of each gene comes from the caller's models (`genecall._finder`), so this
  value is only recorded.
- **Assembly quality**: `checkv_quality`: Complete 24,404 → `complete`; High-quality 34,257, Medium-quality 43,106,
  Low-quality 66,491 and Not-determined 278 → `vOTU`, the schema's value for a vOTU representative (*measured*). The
  CheckV tier itself is not kept in the stage (the enum has no tiers for viruses); it is in `records.tsv`.
- **Topology**: `checkv_completeness_method`: DTR 23,995, ITR 1,280, AAI-based 133,312, HMM-based 9,508, Provirus 163,
  NULL 278 (*measured*). A DTR (direct terminal repeat) genome is circular or terminally redundant. The gene caller
  wraps a circular contig end to start (`genecall._genes`), which is right only when the second copy of the repeat
  is gone. The script tests each DTR record (its first 20 bases found again in its last 2,000): **22,517 no longer
  carry the repeat** (UHGV trimmed it) and are listed circular (988,370,473 bp); 1,478 still carry it and are listed
  linear; every other record is linear (*measured*). The test is a heuristic of this card, not UHGV's statement;
  S1 and H1 list every phage linear.
- Not carried (the inputs format cannot): vOTU id, host lineage, lifestyle, taxonomy, source study, CheckV
  completeness. They are in UHGV's tables (fetched) for a later dimension card.

**What the representatives stand for** (*measured*, `V1/overlap.json` from `uhgv_metadata.tsv`): 91,867 of the
168,536 vOTUs have one genome; the median vOTU has 1, the largest 7,588. By CheckV tier and geNomad confidence:

| CheckV | Confident: contigs, bp | Uncertain: contigs, bp |
|---|---|---|
| Complete | 24,157; 1,080,437,045 | 247; 13,653,244 |
| High-quality | 33,357; 1,467,750,391 | 900; 39,038,515 |
| Medium-quality | 40,840; 1,328,578,136 | 2,266; 68,983,595 |
| Low-quality | 51,522; 972,738,909 | 14,969; 247,471,768 |
| Not-determined | 51; 363,090 | 227; 3,274,083 |

Low-quality and not-determined representatives are 66,769 contigs and 1,223,847,850 bp (*derived*: the sum of their
rows), 40% of the contigs and 23% of the bases. Card C-38's filter at ingest (item 4) would act here; this card keeps
every representative, as the card says.

**Source studies** of the representatives (*measured*): MGV 59,384, GPD 38,990, mMGE 14,231, IMGVR4 13,730, Hadza
13,610, CHVD 10,653, GVD 5,607, HuVirDB 4,910, COPSAC 4,478, DEVoC 1,591, GFD 759, Benler 593. **13,730
representatives (387,561,224 bp) come from IMG/VR v4** inside a CC BY 4.0 record: an open question for the PI (§8).

## 4. What V1 shares

### By `contig_uid`

The store's identity (`ids.contig_uid`, computed by `c39_inputs.py`; no strand or rotation canonicalisation),
against every stage's `contig_index` (owners, default variants) and `genome_input` (members), `c39_overlap.py`
(unit `c39-overlap2`), all *measured*:

| Stage | V1 contigs it owns | bp | V1 contigs its genomes list | of them not owned by it |
|---|---|---|---|---|
| S1 | 1 | 61,892 | 1 | 0 |
| P1 to P6 | 0 | 0 | 0 | 0 |
| H1 | 177 | 43,733,807 | 177 | 0 |
| PL1 | 0 | 0 | 0 | 0 |
| **all** | **178** | **43,795,699** | | |

- The one S1 contig is an IMG/VR v4 genome (High-quality); the 177 H1 contigs are huge phages UHGV took from MGV 85,
  Hadza 70, GPD 7, IMGVR4 8, GVD 4, CHVD 2, COPSAC 1 (10 Complete, 131 High, 36 Medium) (*measured*).
- Under DEC-26 and DEC-27 these 178 are second genomes of contigs the earlier stage keeps: V1 writes no calls for them.
- No two V1 records share bytes (no repeats) (*measured*).
- **Between the sources by `contig_uid`**: not possible here. The soil atlas was not fetched, and VIRE's table holds
  no sequence or hash.

### At the vOTU threshold (skani)

`c39_ani.py` (unit `c39-ani`, 5 threads, 14:35:13 to 14:38:54; tables summarised again by `c39-ani2`, `c39-ani3`),
skani 0.4.0 of the image, `dist -c 30 -m 200 --min-af 0 --qi --ri` (skani's small-genome settings, without
`--faster-small`, which drops short genomes from the screen). Samples, deterministic: every 20th V1 contig (8,426,
260,992,837 bp), every 20th record of S1's input FASTA (1,703) and of H1's (361). Run A: the V1 sample against all of
S1 and H1; run B: the three samples against all of V1. Three rules on skani's ANI and aligned fractions: **shorter**
(the vOTU threshold: ANI ≥ 95, the shorter genome ≥ 85% aligned), **both** (ANI ≥ 95, both ≥ 85%), **near-identical**
(ANI ≥ 99, both ≥ 90%). All *measured* (`V1/ani.json`):

| Queries | n | shorter | both | near-identical |
|---|---|---|---|---|
| V1 sample → S1 | 8,426 | 107 | 39 | 2 |
| V1 sample → H1 | 8,426 | 59 | 28 | 26 |
| V1 sample → S1 or H1 | 8,426 | **163 (1.93%)** | **67 (0.80%)** | **28 (0.33%)** |
| S1 sample → V1 | 1,703 | 322 (18.9%) | 230 (13.5%) | 25 (1.5%) |
| H1 sample → V1 | 361 | 71 (19.7%) | 57 (15.8%) | 33 (9.1%) |
| calibration: V1 sample → other V1 | 8,426 | 2,922 (34.7%) | 728 (8.6%) | 2 (0.02%) |

- **Calibration.** UHGV's representatives are, by construction, different vOTUs, yet 34.7% of them meet the 95/85
  rule of the shorter genome with another representative under skani, and 8.6% with both genomes covered; the pairs
  sit just over the line (ANI 95.1 to 97.1, fractions 86 to 98, in the 13 pairs read). UHGV clustered with MCL on an
  ANI graph (its GitHub page), which does not promise that representatives stay apart at 95/85, and skani's ANI and
  fractions are estimates from sketches and chained anchors, not alignments. So the 95/85 counts above are **upper bounds** on "is a vOTU of", and the
  near-identical column, which the calibration passes (2 of 8,426), is the lower one.
- Extrapolated to all of V1 (*derived*: fraction × 168,536): between about 560 (near-identical) and 3,260 (shorter
  rule) representatives are vOTUs of a phage the root already holds, 1,340 by the both-covered rule; the 178 exact
  shares are among them. By CheckV tier, the near-identical V1 matches are Complete 16, High 10, Medium 2, Low 0.
- The other way, 13.5% to 18.9% of S1's isolate phages and 15.8% to 19.7% of H1's huge phages have a UHGV vOTU
  (the 95/85 rules), 1.5% and 9.1% near-identically: H1's huge phages were in large part assembled from the same gut
  catalogues.
- Not done: the soil atlas against V1 and the root (not fetched); a second, alignment-based oracle (MMseqs2 or
  BLAST ANI over HSPs, as MIUViG's own tool) to replace the skani calibration.
- CPU: run A 106.3 s, run B 370.0 s of skani (*measured*, getrusage of the children).

## 5. What building V1 would cost

`c39_est.py` (unit `c39-est`, one core, 14:35:24 to 15:33:32): the code's own functions from the tree
(`biobox.stage`, `biobox.l1.genecall`), a sample called with the step's own call path, and the provider's proteins.

| Step | Estimate | How |
|---|---|---|
| ingest (one object of 5,228,187,536 bytes) | **23.44 CPU-min** | *derived*: `stage.ingest_cpu_min` |
| gene calling, gv-meta, every genome whole | **687.5 CPU-min, 3 batches** at 240 | *derived*: `genecall.one_process_cpu_min` over 4,233,918,303 linear bases and 2,494,857,527 bases read for the 22,517 circular genomes (their 988,370,473 bp with the 50 kb wrap flanks); the build runs it in batches (card B-48) |
| the same at the rates measured here | 549.9 CPU-min | *derived*: (3.7222 × 4,233.918 + 11.1929 × 988.370) / 60 × 1.23 (`STEP_OVERHEAD`); rates *measured* below |
| genes | about 7.39 M | *derived*: 1,426.36 genes per Mbp (*measured*) × 5,178.49 Mbp; the code's 1,441.1 gives 7.46 M; the provider's 7,426,124 proteins over all 5,222.29 Mbp (*published*, counted here) |
| distinct proteins | **5.60 M**, 0.17 M of them already in the root | *measured*: 5,597,468 distinct sequences of the provider's 7,426,124 (`ids.protein_uid`, trailing `*` removed), 174,820 of them in the root's `protein_seq` (3.1%); our calls are pyrodigal-gv, the provider's prodigal-gv |
| residues of the distinct proteins | 1,307 M | *measured*: 1,307,249,946 |
| Pfam 38.2 scan | **about 10,400 CPU-min in 44 shards** | *derived*: 1,307.25 M × 7.96 (`SCAN_CPU_MIN_PER_M`) = 10,406; shards of 30,150,753 residues (300 × 0.8 × 10⁶ / 7.96): 43.4, so 44. Without the proteins the root has: about 10,080 and 43 (*estimate*, scaled by 5,422,648 / 5,597,468) |
| CRISPR, default and filtered producers | **15.71 and 18.69 CPU-min** | *derived*: `stage.crispr_cpu_min(5,178,493,077)`, × `CRISPR_NEAR_CAS_COST` |
| IGR, tapes | about 0.9 and 1.0 CPU-min | *derived*: H1's 0.267 and 0.318 (*measured*, card C-24) × 7.39 M / 2.28 M genes |
| files in the root | **about 3.7 to 4.2 GB** | *derived*: 5.178 Gbp × 0.7175 (S1) to 0.7316 (H1) GB per Gbp of stored bases (*measured* here, `du`-equivalent of every `stage=` file over `contig_index` bases), raised by up to a tenth because V1 has more distinct proteins per Mbp (1,072 against H1's 884) |
| in all | **about 11,150 CPU-min** (186 CPU-hours) | *derived*: the sum, with the code's gene-calling figure; about 5 hours at 36 threads (*estimate*) |

**Storage per Gbp** (*measured*, `V1/est.json`, bytes of every file under a `stage=` path in the root over the stage's
stored bases): S1 0.7175 and H1 0.7316 GB per Gbp (the virus stages), P1 to P6 0.775 to 0.780, PL1 0.598.

**The gene-calling sample** (every 30th contig the stage stores, FASTA order: 5,611 genomes, 172,030,086 bp, called
one at a time with `genecall._call` as `genomes.json` gives them, one thread; *measured*): 245,377 genes, 1,426.36
per Mbp; calls alone 898.5 CPU-s, **3.7222 s per Mbp on the 4,833 linear genomes and 11.1929 s per Mbp of contig on
the 778 circular ones** (34.56 Mbp); 238,910 distinct proteins, 7,761 of them in the root.

**Against the code's rates.** Linear: 3.72 against the code's 4.2019 (11% under). Circular: the code has no measured
gv-meta circular rate and takes the linear rate times meta mode's ratio (6.3075 s per Mbp read, *estimate* in
`genecall.py`), which is 15.92 s per Mbp of contig for V1's circular genomes (*derived*: 6.3075 × 2,494.86 / 988.37);
measured 11.19, 30% under. The whole step at the measured rates is 549.9 against 687.5, **20.0% under the code's
estimate**, at the edge of AGENTS §6's 20%. The code errs on the safe side (it sizes batches and refusals), so the
build is not at risk; the measured circular gv-meta rate (11.19 s per contig Mbp on 778 genomes of mean 44.4 kb) is
the first such measurement and could replace the estimate in `CALL_S_PER_MBP` (a code card, not this one).

No single process passes 300 CPU-minutes: gene calling runs in 3 batches of at most 240 (code), each Pfam shard at most
240, ingest 23.4.

## 6. What VIRE would add (its table alone)

`c39_vire.py` over VIRE's all-samples genome table (unit `c39-vire`, 38 s), all *measured* (`vire.json`) unless said:

- 1,784,510 genomes (the paper's figure, *published*), 68,027,010,177 bp; CheckV Complete 384,035, High 417,105,
  Medium 983,370 (no low-quality genomes, unlike UHGV's representatives); geNomad topology DTR 393,553, ITR 28,247;
  codes 11 1,752,608, 15 25,999, 4 5,903; every row `source` metagenome.
- 697,136 species clusters (`species_cluster`; the paper's 706,281 includes reference genomes, *published*), 491,505
  of one genome; 520,463 genus clusters.
- **Habitat** (`microntology`; human gut = the labels `...:human host` and `...:digestive tract:intestine`, soil = a
  label naming soil): human gut 950,399 genomes (the paper's figure exactly) in 255,274 clusters; soil 48,679 in
  33,492; other 785,432 in 408,370 (aquatic 308,797 genomes, other animal guts and the mouth, built environment,
  wastewater and more).
- **Against IMG/VR** (VIRE's own `alignment_to_imgvr`; True 743,710, False 990,069, None 50,731 counted as not
  aligned): clusters with no member aligned: human gut 208,244, soil 20,514, other 339,286; 567,044 in all.
- **A stage of VIRE's species representatives** would be about 26.3 Gbp (*estimate*: the longest member of each cluster
  summed, 9.73 Gbp gut, 1.28 soil, 15.33 other, since the table names no representative), five times V1; the
  clusters no member of which aligns to IMG/VR are about 20.6 Gbp of it.
- **What it adds over V1 and V2**: not measurable from tables, since nothing joins VIRE's genomes to UHGV's or the
  atlas's. In outline: its 408,370 clusters outside the human gut and soil are habitats neither source covers (and that
  IMG/VR covers only in part); its 255,274 gut clusters are medium quality or better, where 40% of V1's representatives
  are low quality; its 33,492 soil clusters are close to the
  atlas's 38,508 vOTUs (*published*) in number, from other samples in part. A sequence comparison needs VIRE's genome
  archives (`all_samples_genomes.tar`, 20 GB on its downloads page, *published*), which the entry does not cover.

## 7. What the virus stages lack without IMG/VR

IMG/VR v4.1 is the planned stage S4 (231,408 high-quality vOTU representatives, *published* in `STAGES.md`). Without
it: viruses of every IMG biome other than the human gut and (once V2 is fetched) soil: marine, freshwater,
engineered and host-associated non-human habitats; IMG/VR's proviruses with their boundaries (S4's capability list);
and IMG/VR's host predictions. From the measured overlaps: 13,730 of V1's representatives already come from IMG/VR v4
(§3), and VIRE marks 743,710 of its genomes as aligned to IMG/VR, so IMG/VR overlaps what VIRE holds by about 42% of
its genomes (*derived*: 743,710 / 1,784,510).

**What is needed from the PI for IMG/VR** (DEC-50 (2)): from JGI's portal, with the PI's login and acceptance of
JGI's data usage policy, release `IMG_VR_2022-12-19_7.1`'s high-confidence nucleotide FASTA and its sequence
information table (the `.tsv` that gives vOTU, quality, topology, host and ecosystem), copied to the pod (for example
`/tmp/if2453/s3/dl/imgvr/`) with JGI's checksum if it gives one; then an entry `imgvr` in `sources.yaml`
(`redistributable: false`), and a stage of its own.

## 8. Open questions

1. **The soil atlas**: PNNL refuses the pod (§1). Fetch the archive from another network, or allow a fetch from the
   laptop and a copy to the pod?
2. **IMG/VR genomes inside UHGV**: 13,730 V1 representatives (387.6 Mbp) are IMG/VR v4 sequences published by UHGV
   under CC BY 4.0. Whether JGI's terms follow them, and so whether V1 is `redistributable` as a whole, is the PI's
   call (DEC-51 lets us use them either way; the inputs mark the object redistributable, as the record states).
3. **Low-quality representatives**: 40% of V1's contigs (23% of bases) are CheckV low-quality or not determined.
   Keep them (as now) or drop them under card C-38's filter at ingest.
4. **Circular DTR genomes**: listed circular when the terminal repeat was trimmed (22,517). Accept the rule, or list
   every V1 genome linear as S1 and H1 are.
5. **Quality value**: `vOTU` for every non-complete representative (the schema's value); CheckV's tier is not kept.
6. **The code's circular gv-meta rate** is about 30% above what was measured (§5): a code card may replace it.

## 9. Pod use and scratch

| Unit | CST | What | Threads |
|---|---|---|---|
| `c39-fetch` | 13:49 to about 13:58 (stopped after the refusals) | UHGV, then the soil atlas's refusals | 1 |
| `c39-fetch-vire` | 14:00 to 14:01 | VIRE's table | 1 |
| `c39-inputs` | 14:28 (failed: zstd window), 14:28:57 to 14:30:15, 14:33:19 to 14:34:22 | the inputs (the second run before the topology rule; removed) | 1 |
| `c39-vire` | 14:29 to 14:30, 14:33 to 14:34 | VIRE's table (the second with the human-gut rule fixed) | 1 |
| `c39-fetch-faa` | 14:31 to 14:32 | UHGV's proteins | 1 |
| `c39-check` | 14:34:32 to 14:34:57 | the dry check, sha256 of the inputs | 1 |
| `c39-overlap`, `c39-overlap2` | 14:35 to 14:39 | overlap by `contig_uid` | 2 |
| `c39-ani`, `-ani2`, `-ani3` | 14:35 to 14:43 | skani; the summaries again | 5 |
| `c39-est` | 14:35:24 to 15:33:32 | estimates, sample, proteins (at most one core for 58 minutes) | 1 |
| `c39-rerun` | 14:43:14 to 14:44:40 | determinism | 1 |

Every unit: `systemd-run --user`, `nice 19`, `ionice` 7, a memory cap (4 to 32 GB), runtime image
`bio-box_faec1f4.sif`, `PYTHONPATH` the tree's `src` alone, the data root and the downloads bound read-only. At most
8 threads of ours at once; no `c35-*` unit was active at any check; no `STOP` in the data root at any step. Tests:
`pod_pytest.sh C-39 tests/r1/test_make_inputs.py` (10 passed; `make_inputs.py` unchanged).

Disk (*measured*, `du -sb`, 2026-11-17 08:15): downloads `/tmp/if2453/s3/dl` 2,944,454,023 bytes; work
`/tmp/if2453/s3/work` 654,412,833; inputs `/tmp/if2453/s2/inputs/V1` 5,268,969,342 (the FASTA 5,228,187,536). About
8.9 GB in all; `/tmp` had 568 GB free at 08:15 on 2026-11-17 (`df -h`). Left on the pod: the downloads, the inputs,
the work directory. Nothing on `/stor`: the archive of the downloads and inputs onto `/stor` (as `provisional/s2_inputs`
for H1) is for the build card.
