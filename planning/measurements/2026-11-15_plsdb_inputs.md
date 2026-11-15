# Stage PL1 (PLSDB plasmids): its inputs, what it adds, and how its build will end (card C-26, part 1)

Date: 2026-11-15, wilkcomp02, 09:51 to 10:35 CST (the pod's clock). Every number carries a label: *measured* (host,
date, command or file), *derived* (formula and inputs) or *estimate*. **No build of the stage**: nothing was written
in the data root `/tmp/if2453/r10/root` (bound read-only into every command that read it) and nothing under `/stor`.
Stage H1 (card C-24) is half-built and stopped (STATUS open item 7); its units, scratch and files were not touched.

Small results are committed beside this file in [`2026-11-15_plsdb_inputs/`](2026-11-15_plsdb_inputs/): `PL1/`
(the inputs' `objects.json`, `repeated_records.tsv`, `shared_records.tsv`, the script's `report.json`, the dry
`check.json`, the inputs' `SUMS`, `c26_analyze.py`'s `summary.json`, `c26_rate.py`'s `rate.json`, the pod's
`driver.log`), `rehearsal/` (the sample's list and inputs records, the three units' logs, the check's output).

## Summary

- **Inputs**: `/tmp/if2453/s2/inputs/PL1/` on the pod, made by `make_inputs.py plasmid` (new mode, this card):
  **72,556 records, 68,693 genomes, 68,693 contigs listed, 7,297,669,596 bases in records** (*measured*,
  `c26/inputs/report.json`). A record whose bytes repeat an earlier record's is the same genome (`genome_uid` hashes
  the contig set): 3,863 records, listed once, every one still ingested as an occurrence. Kind `plasmid`, genetic code
  11, meta-mode gene calling.
- **Topology** from `NUCCORE_Topology`: circular 66,511 records, linear 6,045, **empty 0** (*measured*); by genome
  62,797 circular, 5,896 linear (*derived*: 68,693 − 62,797).
- **Quality** from `NUCCORE_Completeness`: `complete` 68,703 records → `complete`; empty 3,853 → `contig`; no other
  value (*measured*); by genome 64,936 complete, 3,757 contig.
- **Shared with an earlier stage** (the stage whose `contig_index` holds the default variant owns the contig,
  DEC-26): 7,139 records, 7,046 contigs, 1,177,838,866 bases, all in P1 to P6 and S1; **none in H1** (*measured*).
  These are second genomes of contigs the earlier stage keeps (DEC-27).
- **What the stage adds**: **61,647 contigs, 5,924,552,875 bases** (*measured*); expected about 6.4 M genes in meta
  mode and 1.4 to 1.8 M distinct proteins (*derived*, below).
- **Cost** (*derived* from the code's own rates unless said): ingest 33.3 CPU-minutes; **gene calling about 860
  CPU-minutes in one step** (*derived* from a rate measured here); Pfam 3,261 to 4,149 CPU-minutes in 14 to 18 shards;
  the two CRISPR producers 18.0 and 21.4; root files about 4.6 GB at most, scratch about 3.8 GB. **The gene-call step
  would pass the 300 CPU-minute rule about threefold, most of it in one process** (section 4): the build cannot run as
  the card is written until that is decided. Every other step fits.
- **How the stage-only launch ends**: once every other stage of `rM.3` is gene-called the cluster gate no longer
  stops `build-stage`; a STOP file written by the driver as soon as the stage's last step has published ends the
  launch at the release's first step (exit 75), with no file of the release written. **Rehearsed**: it held.
- **Rehearsal** (315 PLSDB records as stage `PLR` in a scratch view of P1 and P2, the real command's form, image
  `bio-box_faec1f4.sif` by digest): DEC-26 held (the 103 contigs P1 and P2 own got no row of PLR's), DEC-27 held
  (their 103 second genomes are in `genome_input` and `genome_dim`; `contig_dim` one row per contig, `n_genomes` 2),
  counts by locus do not fan out (9,136 genes by locus; a join to membership gives 18,272).

## 1. Setup

- Branch `task/C-26` from `origin/plan/task-cards` at `f5621cf`. Tree synced to `/tmp/if2453/ci/C-26` by
  `scripts/pod_pytest.sh C-26 tests/r1/test_make_inputs.py`; scripts run from copies in `/tmp/if2453/c26/bin`.
- Python ran in the runtime image `/stor/home/if2453/bio-box/images/bio-box_faec1f4.sif` (the rehearsal's launches
  through `pod_run.sh --digest sha256:d67300aa…f098a`, which checked it), `PYTHONPATH` the tree's `src` alone; the data
  root, PLSDB's metadata and the download bound read-only.
- Units `c26-inputs`, `c26-check`, `c26-analyze`, `c26-analyze2`, `c26-rate` (one or two CPUs, 8 to 12 GB caps,
  nice 19, ionice 7), the rehearsal's `c26-rh-driver-*`, `c26-rh-1`, `c26-rh-2`, `c26-rh-dim` (at most 7 threads,
  32 GB). Never more than 8 threads of ours at once; nothing else of ours computed (`top` at 09:51 and 10:09: load
  0.22 and 1.45, only our own processes). `/tmp` had 225 GB free at 09:59.
- Nothing fetched from outside: the PLSDB files are the ones of the S2 inputs note (`/tmp/if2453/s2/dl/`,
  `/tmp/if2453/s2/plsdb/meta/`), archived under `provisional/s2_inputs/` on `/stor` (not touched in this session).

## 2. The inputs

Script: `scripts/stress/build_stage_pod/make_inputs.py plasmid FASTA META_CSV OUT_DIR --root ROOT`, run by
`c26_inputs.sh inputs` (unit `c26-inputs`, 09:59:48 to 10:01:33 CST). It first checked the FASTA's sha256 against
`sources.yaml` (`e654c106…a587a`, equal), then hard-linked it into the inputs directory as
`plsdb_2024_05_31_v2_sequences.fasta`, the name of its copy on `/stor` (same inode as `/tmp/if2453/s2/dl/
plsdb_download_fasta`, no second copy), so the source object's name in `l0_manifest` will be the archive's.

| File (`/tmp/if2453/s2/inputs/PL1/`) | Content | sha256 |
|---|---|---|
| `genomes.json` | 68,693 genomes, one contig each, kind plasmid, code 11, 62,797 circular, quality complete or contig | `69ec5621edcef675398d1c3cbd07989ccdd80b96a3dd94110ab387f4494229e9` |
| `objects.json` | one object: the FASTA, 7,425,668,755 bytes, `source_name` plsdb, `source_release` 2024_05_31_v2, licence "INSDC records; PLSDB states no licence", redistributable, retrieved 2026-11-11T23:14:30Z (the end of its fetch) | `c70d2ec6470394707f976c5f5aef86f05c8c04f8e0e6d86934ff84707b51150a` |
| `repeated_records.tsv` | the 3,863 records whose genome is an earlier record's: ordinal, accession, contig, the first record, whether the two agree on topology and quality | `3d04f91367f2d582dff9837c5ba98c2a457844e293e03512d8fe5880d20396bb` |
| `shared_records.tsv` | the 7,139 records whose contig an earlier stage owns, with that stage | `4567837327e4709e3189e99bbbd97ad2a6ed034ca0b9d28bf6217d37589ef297` |
| `plsdb_2024_05_31_v2_sequences.fasta` | the PLSDB FASTA (hard link) | `e654c106d8a15a2d5d55a8b680087377cb84e654e63d41bca181fc69a23e587a` |

The H1 inputs hold `dropped_records.tsv` (records the old rule did not list); PL1 drops none, so its two record files
are the repeats (one genome, several records) and the shared ones (a genome whose contig another stage keeps).

Dry check (`s2_check_inputs.py PL1`, unit `c26-check`, the inputs step's parser and `_rows` checks, writes nothing):
exit 0, 68,693 genomes, 68,693 `genome_input` rows and distinct genome uids, 62,797 circular rows, kinds `[plasmid]`,
codes `[11]`, quality `[complete, contig]`, the object a file of 7,425,668,755 bytes, no problem (*measured*,
`PL1/check.json`).

### Topology: `NUCCORE_Topology`

The only topology field in the metadata (`nuccore.csv`, 21 columns; the FASTA headers carry none). Values by record
(*measured*, `c26_analyze.py`, `metadata_values_by_record`): **circular 66,511, linear 6,045, empty 0**. The script maps
circular to `circular` and linear (or empty, none here) to linear; any other value would be refused. By genome (the
first record of each sequence): 62,797 circular. Of the 3,863 repeats, **112 disagree with their first record on
topology** (*measured*); the genome takes the first record's (FASTA order), which is a rule of the script, not a fact
of the molecule. For the 7,046 contigs an earlier stage owns, `contig_dim` takes the topology from the owning stage
(`dim._views`: rows of the owning stage), which lists every contig linear (the prok inputs carry no topology); 6,249 of
them PLSDB calls circular (*measured*). Circular contigs are called with the wrap flanks, so a gene across the origin
is one gene (`genecall._genes`).

### Quality: `NUCCORE_Completeness`

By record: `complete` 68,703, empty 3,853 (*measured*). Mapped to the schema's `assembly_quality` enum:
**complete → `complete`, empty → `contig`** (a record NCBI does not call complete is an assembled sequence of unknown
completeness; the enum's `plasmid` value would repeat the kind and say nothing of quality). Any other value would be
refused. By genome: 64,936 complete, 3,757 contig; on the contigs the stage stores, 58,259 and 3,388 (*measured*). 99
repeats disagree with their first record on quality (*measured*).

### Genetic code and gene-calling mode

- **Code 11** for every genome (the card's value; PLSDB's metadata has no code per record; `proteins.csv` gives PGAP's
  `transl_table` per protein, not read here). It matters little: for kind `plasmid` the code calls in **meta mode at
  any length** (`genecall._finder`: `if genome.kind == "plasmid" or ... < SINGLE_MODE_MIN_BP: GeneFinder(meta=True)`;
  measured on a 4 kb and a 400 kb plasmid genome: `meta`, `c26_analyze.py`), and meta mode takes the code per gene
  from its models, so `genetic_code` is only recorded.
- Contigs another stage owns were called there in its mode (single mode, trained on the whole P genome); PL1 calls them
  again whole in meta mode, for nothing: the calls are dropped before writing (DEC-26), as the rehearsal showed.
- PLSDB marks 7,658 records `VIRAL_putative` True (`VIRAL_category` likely_viral_1 7,413, likely_viral_2 235,
  likely_viral_1-2 10) (*measured*). They stay kind `plasmid` (the card's kind); were they phage–plasmids to be called
  as viruses, the inputs format could say so (`kind` virus, gv-meta), at the cost of the stage's one kind. Not done.

### What the stage's `genome_dim` could carry, and what the inputs format holds

`genome_dim` will hold, per PL1 genome: `n_contigs` 1, `genome_len`, `kind` plasmid, `genetic_code` 11,
`assembly_quality`; `tier`, `tax_order` null, `is_species_rep` false (no step fills them). The accession and NCBI title
of every record ride in `contig_occurrence` and `contig_label` (from the FASTA header). PLSDB's metadata has more that
a plasmid genome's row would want, and **the inputs format (`GenomeInput`: contigs, code, kind, circular, quality,
tier, species-rep) can hold none of it**: host taxonomy (`TAXONOMY_UID` → `taxonomy.csv`; Bacteria 71,949, Archaea
607 records, S2 note), the host assembly (`ASSEMBLY_UID`, 57,392 records linked: it would tie a plasmid to the P-stage
genome it came from), `NUCCORE_Source` (refseq 59,471, insd 13,085), `STATUS` (196 suppressed records, kept: the card
drops none), `VIRAL_putative`, the BioSample ecology, MOB-typer mobility and replicon types, PlasmidFinder replicons
and AMR genes. `genome_dim.tax_order` exists but no step writes it. No schema change made; reported.

## 3. Records against the earlier stages

Ownership read from every stage's `contig_index` on disk (default variants; `make_inputs.py --root`), membership
from every stage's `genome_input` (`c26_analyze.py`, unit `c26-analyze`, 10:02:04 to 10:03:44). The identity is
`ids.contig_uid` (BLAKE3-128 of the uppercase bytes), the store's; no strand or rotation canonicalisation.

| Stage | records whose contig it owns | contigs | bases | contigs its genomes list but it does not own |
|---|---|---|---|---|
| S1 | 2 | 2 | 202,881 | 0 |
| P1 | 120 | 118 | 18,650,102 | 0 |
| P2 | 1,294 | 1,267 | 227,508,917 | 0 |
| P3 | 1,189 | 1,174 | 160,535,406 | 1 |
| P4 | 1,513 | 1,503 | 261,879,972 | 1 |
| P5 | 1,713 | 1,692 | 303,733,594 | 1 |
| P6 | 1,308 | 1,290 | 205,327,994 | 1 |
| **H1** (half-built; its `contig_index` published) | **0** | **0** | 0 | 0 |
| **all** | **7,139** | **7,046** | **1,177,838,866** | 4 |

- All *measured*; the sums are *derived* and equal the script's own totals (7,139 records; genome bases 7,102,391,741 −
  new bases 5,924,552,875 = 1,177,838,866).
- The four contigs listed by a P stage that another stage owns are the four of card C-20 (P1/P3 and P2/P4 plasmids,
  and the two S1 phage–plasmids also in P5 and P6). That is why the S2 note's per-stage counts by membership (P3
  1,191, P4 1,516, P5 1,714, P6 1,309) are one to two higher than these by ownership; the total, 7,139, is the same.
- No PLSDB record is an H1 contig (H1's 6,714 stored contigs are phages over 200 kbp). H1's state does not change any
  number here, since its `contig_index` is published and will not change.

**Records that repeat each other within PLSDB** (*measured*, `c26_analyze.py` `repeats`): 2,003 sequences are
carried by more than one record, 5,866 records in all (3,863 beyond the first of each); groups of 2 records 1,382,
of 3 308, of 4 127, …, the largest 57 records. **What DEC-27 makes of them with the code as it is: one contig, one
genome** (the PL1 genome of a single-contig record is identified by its contig set, so identical records are the same
genome, `ids.genome_uid`; `genome_input._rows` refuses a genome listed twice), and every record an occurrence of the
contig (`contig_occurrence`, `contig_label.n_records`), so no accession is lost. The "several genomes" of DEC-27 are
across stages here: a PL1 genome {plasmid} beside the P genome {chromosome, plasmid, …} that holds the same contig.

## 4. What the stage adds, and what it will cost

| Quantity | Value | Label and source |
|---|---|---|
| contigs the stage stores | 61,647 | *measured*, `report.json` `new_contigs` |
| their bases | 5,924,552,875 | *measured*, `new_bases` |
| genomes, of which second genomes of contigs another stage owns | 68,693; 7,046 | *measured* |
| bases gene calling reads (every genome called whole) | 7,102,391,741 | *measured*, `genome_bases` |
| genes on the stored contigs, meta mode | **about 6.37 M** | *derived*: 5,924.55 Mbp × 1,074.7 genes per Mbp, *measured* by `c26_rate.py` (the code's `genecall._call`, one thread, on every 30th stored record: 2,054 genomes, 186,854,199 bases, 200,816 genes) |
| the same at the density of the plasmids P1 to P6 already called | about 5.65 M | *derived*: × 954.0 genes per Mbp, the 7,044 PLSDB contigs P1 to P6 own (1,123,460 genes on 1,177.6 Mbp, *measured* from their `feature_protein`; single mode, genome-trained) |
| the same from PGAP | about 5.61 M | *derived*: 5,005,221 PGAP CDS on the 59,767 annotated new records (5,286.3 Mbp; `proteins.csv`, *measured*; 946.8 per Mbp) × 5,924.55 / 5,286.31 |
| distinct proteins of the stage | 1.38 to 1.76 M | *measured* lower bound: 1,383,442 distinct PGAP protein sequences of the 59,767 annotated new records (`proteins.fasta`, hashed with `ids.protein_uid`); *derived*: × 1.1207 for the unannotated 10.8% of bases (1.55 M), × 1.1351 more for meta mode's genes over PGAP's (1,074.7 / 946.8 per Mbp; 1.76 M) |
| their residues (what the Pfam scan reads) | 409.7 to 521.2 M | the same three ways (409,706,637 *measured* on PGAP's; 459.2 M; 521.2 M) |

Meta mode calls 12.7% more genes per Mbp than the P stages' genome-trained calls of the same kind of molecule
(1,074.7 against 954.0) and 13.5% more than PGAP (946.8), both inside AGENTS §6's 20%. The P stages' own plasmids give
0.95 distinct proteins per gene (1,067,518 of 1,123,460, one plasmid set per
species), PGAP's PLSDB 0.28 (1,383,442 of 5,005,221): PLSDB holds many near-identical plasmids of the same species
(*Escherichia*, *Klebsiella*), so its proteins repeat. Gene counts in meta mode will differ from both (the rehearsal's
185 short plasmids, all under 60 kbp: 3,953 genes on 3.34 Mbp, 1,183 per Mbp, *measured*; small plasmids are denser).
How many of the stage's proteins are new against `protein_seq` of the earlier stages is measured after gene calling
(card C-26's build session).

**Cost of each step** (the code's rates; `stage.py` constants, measured on wilkcomp02 2026-11-10 unless said):

| Step | Estimate | How |
|---|---|---|
| ingest (one process; one object of 7,425,668,755 bytes, under the 240 that would batch it) | **33.3 CPU-min** | *derived*: `ingest_cpu_min` = 0.07 + 7,425,668,755 × 42.9 / 9,596,180,328 |
| gene calling (one step, workers inside) | **about 860 CPU-min** | *derived*: 7,102,391,741 bases (every genome is called whole, the 7,046 whose contig another stage owns too) × 7.26 s per Mbp / 60 = 859.4. The rate is *measured* (`c26_rate.py`, unit `c26-rate`): 1,356.6 CPU-s for 186.85 Mbp, calls only (the step's writing and sorting come on top); circular genomes 7.70 s/Mbp (they are called with the wrap flanks, `WRAP_FLANK` 50 kb or the whole contig on each side, so a small circular plasmid is called on three times its length), linear 2.94. For comparison: meta finding on P1's short contigs 3.285 s/Mbp (C-21f note), H1's whole gv-meta step 4.47 s/Mbp (146.091 CPU-min for 1,961.5 Mbp, linear), P2's single-mode step 0.90 |
| default CRISPR producer | **18.0 CPU-min** | *derived*: `crispr_cpu_min(5,924,552,875)` = bases × 182.0 s/Gbp / 60 |
| filtered CRISPR producer | **21.4 CPU-min** | *derived*: × 1.19 (`CRISPR_NEAR_CAS_COST`) |
| IGR, tapes | about 0.7 and 0.8 CPU-min | *derived*: H1's 0.267 and 0.318 (*measured*, C-24 note) × 5.6 M / 2.28 M genes |
| Pfam 38.2 scan | **3,261 to 4,149 CPU-min in 14 to 18 shards** | *derived*: 409.7 to 521.2 M residues × 7.96 CPU-min per M (`SCAN_CPU_MIN_PER_M`); shards of at most 30,150,753 residues (`shard_residues`: 300 × 0.8 × 10⁶ / 7.96); each full shard estimated at 240 and measured on H1 at 0.887 to 0.920 of that, so about 213 to 221; the remainder shard costs about 1.5 CPU-min of model loading whatever its size (C-24 note, shard 13) |
| files added to the root | at most about 4.6 GB | *derived*: 5,924.55 Mbp × 0.780 bytes per stored base, stage P3's files (7,762,670,119 B for 9,949,077,522 bases, *measured* `c26_analyze.py`); PL1's proteins repeat four times more, so its protein and domain tables will be smaller |
| scratch | about 3.8 GB at most | *derived*: 5.39 GB, the largest stage-step scratch of P3 to P6 at about 10 Gbp (C-24 note), × 7.10 / 10 |

At the build's 36 threads (9 shards of 4 for the scan) the wall time is about 3 hours (*estimate*: ingest about 35
minutes on one core, gene calling about 25, the scan 3,261 to 4,149 / 36 ≈ 90 to 115 minutes plus its tail, the rest
minutes). In all about 4,200 to 5,100 CPU-minutes (*derived*: the sum of the rows).

**Whether any single process would pass 300 CPU-minutes: yes, gene calling.** The code makes no estimate of it and has
no way to split it: `build-stage` sizes the scans, ingest and clustering against the budget and refuses a CRISPR step
above the margin, but gene calling is one `run-step` step whatever its size. At about 860 CPU-minutes it is 2.9 times
the 300 of AGENTS §7.2 (the run record would say `over_cpu_budget: true`), and it is not spread over processes: a
genome of at least 10,000 bases (`THREAD_MIN_BASES`) is called in a thread of the step's own process, shorter ones in
worker processes (`genecall` module docstring, `_Workers`), and nearly all of PLSDB's bases are in genomes over 10 kbp
(the median record is 53,666 bp, S2 note), so almost all of the 860 CPU-minutes land in one process (*estimate*). No
earlier stage came near this: P2's gene calling was 143.3 CPU-min for 9.6 Gbp in single mode, H1's 146.1 in gv-meta
(*measured*, run records). PL1 is the first stage called wholly in meta mode, at eight times P2's rate per base, and
mostly circular. Ingest (33.3), each Pfam shard (at most 240 estimated), each CRISPR step (18.0, 21.4) and IGRs and
tapes fit. **This is for the orchestrator and the PI before the build** (AGENTS §9; not this session's to settle).
The ways out I can see, none tried: (a) the PLSDB records as four stages of about 1.8 Gbp of genomes each (about 215
CPU-minutes of gene calling each, *derived*: 859.4 / 4), which changes card C-35's release from nine stages to twelve;
(b) a code card that runs gene calling as batch steps under the budget, as card B-30 did for ingest; (c) an exemption
for the step, as DEC-29 gives one clustering call. A smaller saving, whichever is chosen: 143 of the 860 CPU-minutes
(*derived*: 1,177.84 Mbp × 7.26 / 60) call the 7,046 genomes whose one contig another stage owns, in meta mode, where
nothing is trained, and the calls are then dropped (DEC-26); skipping a genome all of whose contigs another stage owns
would be a code change.

## 5. How the stage-only launch ends before the release's cluster step

From the code (`src/biobox/stage.py`):

- `build-stage` has no stage-only mode. Its step list (`steps`) ends with the release's `cluster`, `tokens`,
  `postings`, `dim`, `manifest`, `catalog`, and **no stage step's argv names the release or its stages** (`st(...)`
  carries `--root`, `--gene-call-version`, `--stage` and the step's own options), so the stage's files are the same
  whichever release the launch names.
- Before the cluster step, `_cluster_gate` refuses (exit 2) only while a stage of the release has no gene-call run
  record (`release_proteins`). Stage H1 was gene-called on 2026-11-15 (C-24 note); PL1 is gene-called in its own
  launch. So **a launch of PL1 under `rM.3` over the nine stages passes the gate**: it estimates the clustering of all
  nine stages' proteins (about 47.4 M plus PL1's, *derived*: 45,815,843 of `rM.2` + H1's 1,608,023, *measured*),
  routes it into parts (`_cluster_split`, above 240 CPU-minutes), finds no refusal in `cluster_refusal` (no cut is
  published) or `dim.owner_problem`, and would start `rM.3.cluster.plan`. That is the release's work (card C-35), and
  it must not start here.
- The amendment names two ways: a STOP at a step boundary written by the driver when the stage's last step is
  recorded, or the steps run one by one through `run-step`. The second means driving the Pfam scan's plan, its 14 to
  16 shards and their merge by hand, which build-stage does with its allowance, rescans and lost-scratch recovery;
  not chosen.
- **The STOP, and why it lands in time.** The stage's last step is the filtered CRISPR producer,
  `crispr-minced-0.4.2-w8-nr2-r23x47-s26x50-cas1kb-caspfam1-pfam-38.2`. run-step checks STOP once after the command
  (a STOP present then discards the step's work) and only then writes the intent file
  `runs/stage=PL1/<step>.publishing.json` and the record. A watcher that writes `ROOT/STOP` when either appears is
  therefore never early. After the record, build-stage checks STOP at the top of its loop, computes the cluster gate
  (footers of every stage's `protein_seq`), `dim.owner_problem` (DuckDB over every stage's `contig_index` and
  `genome_input`), the parts decision, checks STOP again in `one()`, reads the CPU allowance (`runstep.allowance`, a
  5-second sampling, `CPU_WINDOW_S`), and calls `run_step`, which checks STOP before it writes anything. The window
  from the record to that last check is at least those 5 seconds.
- **On the rehearsal** (*measured*, `rehearsal/c26-rh-1.log`, driver log): the filtered producer published at
  10:13:07; build-stage logged the parts decision at 10:13:07 and passed its own check; the watcher (polling every
  0.5 s) wrote STOP at 10:13:08.023; at 10:13:13 `run-step: …/STOP exists; not starting rC26.1.cluster.plan`, then
  `build-stage: rC26.1.cluster.plan exited 75; stopping here`; `pod_run: c26-rh-1 exited 75`. No path naming `rC26.1`
  in the view or its scratch, no `release_def` of it. The same launch again (`c26-rh-2`, every stage step recorded):
  the watcher wrote STOP 2.5 s after the launch, `build-stage: …/STOP exists; not starting stage PLR`, exit 75, the
  view's listing (path, size, mtime) unchanged.
- **So the build session runs**: the real command (below) with a watcher that (1) writes `ROOT/STOP` with a marker line
  as soon as `runs/stage=PL1/crispr-<filtered producer>.json` or its `.publishing.json` exists, polling at most every
  0.5 s; (2) stops the unit (`systemctl --user stop`) if any path `runs/stage=rM.3` or `release=rM.3` appears, as a
  second guard; (3) after the unit ends, requires exit 75, a line `STOP exists; (stopping before|not starting)
  rM.3.cluster`, every PL1 stage step recorded and no `rM.3` path, then removes only the STOP holding its marker.
  `c26_rehearse.sh` (`launch`, `verify`) is that driver for the rehearsal. The STOP is in the shared data root for
  seconds to minutes; nothing else may build there meanwhile (one build at a time, and C-35 waits for this card).

  ```
  pod_run.sh --profile <profile, tmp under /tmp/if2453/c26> --image .../bio-box_faec1f4.sif \
    --digest sha256:d67300aa1321b571fa8da79e74dfccdd87a8dca2e56518f324902efa6b9f098a --threads 40 \
    --bind /tmp/if2453 --bind $PFAMDIR:$PFAMDIR:ro --env PYTHONPATH=$TREE/src --log ... c26-pl1-1 \
    build-stage --profile <same> --root /tmp/if2453/r10/root --stage PL1 --release rM.3 \
    --stages S1,P1,P2,P3,P4,P5,P6,H1,PL1 --genomes /tmp/if2453/s2/inputs/PL1/genomes.json \
    --objects /tmp/if2453/s2/inputs/PL1/objects.json --library pfam 38.2 $PFAMDIR/Pfam-A.hmm.gz gathering \
    --crispr --threads 36 --scan-jobs 9 --scan-cpus 4
  ```

- An alternative without a race, for the orchestrator to weigh (not in the amendment, not rehearsed): name in the PL1
  launch a stage list with one stage that has no gene calls (`--stages <the nine>,PLX`), so the gate refuses exit 2 as
  it did for H1; no `release_def` is written before the cluster step (C-24 note), so nothing on disk would name the
  ten stages, but the command would name `rM.3` with a list that is not its definition.

## 6. The rehearsal

`c26_rehearse.sh` modes `prep`, `run`, `dim`, `check` (units `c26-rh-driver-*`), 10:06:46 to 10:18:13.

- **The view** `/tmp/if2453/c26/rh/root` (`c26_view.py ROOT … P1,P2`): real directories with 14 symlinks to stage
  directories of the data root (P1 and P2's `contig_index`, `seq_locator`, `seq_shards`, `genome_input`,
  `feature_protein` and `protein_seq` of `pyrodigal-3.7.1-gv-0.3.2-p2`, `runs/stage=P*`) and `store_params.json`;
  every launch bound the data root read-only.
- **The sample** (`c26_sample.py`, FASTA order, records copied byte for byte): 315 records, 288 genomes — 100 whose
  contig P1 or P2 owns, 3 repeat groups P1 or P2 own (7 records), 15 repeat groups no stage owns, 150 new records
  spread over the file (every 50th under 60 kbp), 10 linear and 10 with an empty completeness; none that S1, P3 to P6
  or H1 owns. Its inputs (`make_inputs.py plasmid … --root <view>`, `rehearsal/report.json`): 288 genomes, 220
  circular, 27 repeats (1 disagreeing on topology), 107 shared records (P1 42, P2 65) on 103 contigs (41, 62), 185
  new contigs of 3,340,530 bases, 6 suppressed records.
- **The launch** `c26-rh-1`: `build-stage --stage PLR --release rC26.1 --stages P1,P2,PLR --library pfam 38.2 <the real
  Pfam-A.hmm.gz> gathering --crispr --threads 7 --scan-jobs 1 --scan-cpus 4`, profile with the real budget (300 CPU-min,
  cluster cap 20,000), 8 threads. Every stage step published (*measured*, run records):

  | Step | wall s | CPU-min | peak RSS MiB | result |
  |---|---|---|---|---|
  | inputs | 2.8 | 0.06 | 324 | 288 genomes, 288 contigs |
  | ingest | 3.7 | 0.07 | 324 | 315 records, 185 new sequences |
  | genecall | 20.0 | 2.01 | 399 | 3,953 genes, 3,426 distinct proteins, mode meta, 64 wrapping genes, 103 skipped contigs |
  | crispr default | 4.4 | 0.10 | 324 | 0 arrays, 103 skipped |
  | igr, tapes | 2.7, 2.2 | 0.06, 0.05 | 325 | 3,209 IGRs; 185 tapes, 3,953 slots |
  | domains-pfam-38.2 (one process) | 74.2 | 7.85 | 365 | |
  | crispr filtered | 4.7 | 0.12 | 325 | 0 kept (42 candidates far from Cas), 103 skipped |

  Then the STOP (section 5), exit 75; the unit's totals `MemoryPeak=1831735296 CPUUsageNSec=631365819000` (10.5
  CPU-minutes). Gene calling's skips: exactly the 103 shared contigs with their owners (41 P1, 62 P2), in the result
  and in the ledger (`owned_by:P1`, `owned_by:P2`), among them `ee816c57…`, the P1/P3 contig of card C-20.
- **The dimension tables** of a scratch release `rC26.2` over P1, P2, PLR (`run-step … release dim`, unit
  `c26-rh-dim`, 3.4 s, 0.06 CPU-min): 5,742 contigs (539 + 5,018 + 185), 2,785 genomes, 103 shared contigs, 153
  circular.
- **The checks** (`c26_rh_check.py`, every oracle the sample's own inputs or the linked stages' files; `problems: []`,
  *measured*, `rehearsal/check.json`):
  - DEC-26: PLR's `contig_index` holds the 185 new contigs and none of the 103; no row on a shared contig in any of
    PLR's feature, `feature_protein`, `protein_occurrence`, IGR, tape or CRISPR files (only `contig_occurrence`, 107
    rows, `contig_label`, 103, and `genome_input`, 103, name them, as they should); P1 and P2's own gene calls cover
    them (4,281 genes on 41 contigs, 4,855 on 62).
  - DEC-27: PLR's `genome_input` lists all 288 genomes, the 103 second genomes among them; `contig_dim` 5,742 rows,
    5,742 distinct, equal to the membership's contigs, `n_genomes` and `genome_uid` equal to the membership's count and
    least genome for every contig; the 103 shared contigs have `n_genomes` 2; `genome_dim` 2,785 rows, one per genome
    of the membership, the 103 second genomes kind plasmid with one contig.
  - Counts on the 103 shared contigs: **by locus 9,136 distinct features, equal to the owning stages' 9,136 gene rows;
    a join of those genes to the membership gives 18,272**, the fan-out a count by locus must avoid (card C-21c's
    readers count distinct); **by genome 156 = the 53 P genomes that hold them + 103 PL genomes**.
  - The catalog-level readers (`genome_member`, `q1_genomes`, `bio-box releases`) need a clustered release and are
    not run here (C-21c's tests cover them on fixtures; card C-35 checks them on `rM.3`).

## 7. Pod use and scratch

| Unit | CST | What | Threads |
|---|---|---|---|
| `c26-inputs` | 09:59:48 to 10:01:33 | sha256 of the FASTA, the hard link, `make_inputs.py plasmid --root` | 1 |
| `c26-check` | 10:01:51 to 10:01:56 | `s2_check_inputs.py PL1`, the inputs' sha256 | 1 |
| `c26-analyze` | 10:02:04 to 10:03:44 | `c26_analyze.py` (owners, members, densities, PGAP, root bytes) | 2 |
| `c26-rate` | 10:06:45 to 10:30:14 | `c26_rate.py`, one thread | 1 |
| `c26-rh-driver-prep` | 10:06:46 to 10:08:20 | the view, the sample, its inputs | 1 |
| `c26-rh-driver-run`, `c26-rh-1` | 10:09:4x to 10:13:13 | the rehearsal's launch and its STOP watcher | 7 |
| `c26-rh-driver-verify1`, `-run2`, `c26-rh-2` | 10:16:07 to 10:16:31 | the first launch's checks; the same launch again | 7 |
| `c26-rh-driver-dim`, `c26-rh-dim` | 10:17:3x to 10:17:38 | `release dim` of `rC26.2` in the view | 2 |
| `c26-rh-driver-check` | 10:18:10 to 10:18:13 | `c26_rh_check.py` | 2 |
| `c26-analyze2` | 10:19:20 to 10:20:54 | `c26_analyze.py` again, with the metadata's value counts | 2 |

At most 8 threads of ours at once (the rehearsal's 7 beside the one-thread rate). Tests: `pod_pytest.sh C-26` (32 GB
scope, the test image of `faec1f4`).

Scratch: `/tmp/if2453/c26` 17 MB at the end (`du -sh`, 10:31), its peak about the same (the view's own files are
small, the stage directories are links); the inputs directory 7.0 GB by `du`, of which 7,425,668,755 bytes are the
FASTA's hard link (no new blocks) and 12.6 MB the four files the script wrote. `/tmp` free: 225 GB at 09:59 and at
10:31 (`df -h`). Nothing in
the data root, nothing under `/stor`, no H1 unit or file touched. Left on the pod: `/tmp/if2453/s2/inputs/PL1/` (the
inputs) and `/tmp/if2453/c26/` (the rehearsal's view, scratch, logs).
